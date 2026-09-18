import os
import json
import logging
import argparse
import shutil
import stat
from pathlib import Path
from datetime import datetime
import concurrent.futures # <-- Added for multithreading

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from typing import Any, Dict, List, Union
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from sklearn.model_selection import train_test_split


def binarize_expression(df_X : pd.DataFrame, ref_columns : List[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref_means = df_X[ref_columns].mean(axis=1)
    ref_stdevs = df_X[ref_columns].std(axis=1)
    
    over_thresholds = ref_means + ref_stdevs
    under_thresholds = ref_means - ref_stdevs
    
    B_plus = df_X.ge(over_thresholds, axis=0).astype(int)
    B_minus = df_X.le(under_thresholds, axis=0).astype(int)
        
    return B_plus, B_minus


def calculate_shapley_values(B : pd.DataFrame) -> pd.Series:
    B_vals = B.values
    m = B_vals.shape[1]
    col_sums = B_vals.sum(axis=0)

    inv_col_sums = np.divide(
        1.0, 
        col_sums, 
        out=np.zeros_like(col_sums, dtype=float), 
        where=col_sums != 0
    )
    
    shapley_values = (B_vals @ inv_col_sums) / m
    
    return pd.Series(shapley_values, index=B.index)


def run_CASh(B_case : pd.DataFrame, B_control : pd.DataFrame, b : int, seed : int = 42, desc: str = "Running permutations") -> tuple[pd.DataFrame, pd.DataFrame]:
    np.random.seed(seed)

    n, k = B_case.shape
    _, h = B_control.shape
    m = k + h 

    shapley_case_observed = calculate_shapley_values(B_case)
    shapley_control_observed = calculate_shapley_values(B_control)

    # d_0 (Signed Observed Difference) and delta_i (Absolute Observed Difference)
    observed_diff_signed = shapley_case_observed - shapley_control_observed
    delta_observed = abs(observed_diff_signed)

    count_beta_gte_delta = np.zeros(n)
    betas = np.zeros((n, b))
    raw_diffs_res = np.zeros((n, b))

    # Added dynamic description to tqdm to track different threads easily
    for r in tqdm(range(b), desc=desc, leave=False):
        idx_case_res = np.random.choice(k, size=k, replace=True)
        idx_control_res = np.random.choice(h, size=h, replace=True)
        B_case_res = B_case.iloc[:, idx_case_res]
        B_control_res = B_control.iloc[:, idx_control_res]

        shapley_case_res = calculate_shapley_values(B_case_res)
        shapley_control_res = calculate_shapley_values(B_control_res)

        diff_res = shapley_case_res - shapley_control_res
        raw_diffs_res[:, r] = diff_res

        # Centered bootstrap difference: beta = |d_r - d_0|
        beta = abs(diff_res - observed_diff_signed)
        count_beta_gte_delta += (beta >= delta_observed)
        betas[:, r] = beta

    bootstrap_se = np.std(raw_diffs_res, axis=1)
    bootstrap_se_safe = np.where(bootstrap_se == 0, 1e-15, bootstrap_se)
    
    zscore_signed = observed_diff_signed / bootstrap_se_safe
    zscore_absolute = delta_observed / bootstrap_se_safe

    raw_p_values = count_beta_gte_delta / b
    _, p_adjusted, _, _ = multipletests(raw_p_values, alpha=0.05, method='fdr_bh')

    beta_dist = pd.DataFrame(data=betas.T, index=[f"beta_{i+1}" for i in range(b)], columns=B_case.index)

    results_df = pd.DataFrame({
            "shapley_case": shapley_case_observed,
            "shapley_control": shapley_control_observed,
            "observed_difference": observed_diff_signed,
            "observed_abs_difference": delta_observed,
            "bootstrap_se": bootstrap_se,
            "zscore": zscore_signed,
            "abs_zscore": zscore_absolute,
            "pval": raw_p_values,
            "pval_adj": p_adjusted
        }, index=B_case.index)

    return results_df, beta_dist


def load_configuration(config_path: Union[str, Path]) -> Dict[str, Any]:
    with open(config_path, 'r') as file:
        return json.load(file)


def split_samples(df_samples: pd.DataFrame, train_share: float, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not 0 <= train_share <= 1:
        raise ValueError("train_share must be between 0 and 1")

    if train_share == 1:
        return df_samples.copy(), df_samples.iloc[0:0].copy()
    if train_share == 0:
        return df_samples.iloc[0:0].copy(), df_samples.copy()

    return train_test_split(
        df_samples,
        train_size=train_share,
        random_state=seed,
        shuffle=True,
    )


def save_experiment_results(results_df: pd.DataFrame, beta_dist: pd.DataFrame, B_case: pd.DataFrame, B_control: pd.DataFrame, df_samples_train: pd.DataFrame, df_samples_test: pd.DataFrame, df_expression: pd.DataFrame, comp_dir: Path, direction: str):
    output_dir = Path(comp_dir, direction)
    output_dir.mkdir(exist_ok=True, parents=True)

    results_df.to_csv(output_dir / f"results.csv")
    beta_dist.to_csv(output_dir / f"beta_dist.csv")
    B_case.to_csv(output_dir / f"B_case.csv")
    B_control.to_csv(output_dir / f"B_control.csv")

    df_expression.to_csv(comp_dir / f"expression.csv")
    df_samples_train.to_csv(comp_dir / f"samples_train.csv")
    df_samples_test.to_csv(comp_dir / f"samples_test.csv")


def save_latest_run(run_dir: Path, output_dir: Path):
    latest_dir = output_dir / "latest"
    if latest_dir.exists():
        def remove_readonly(func, path, exc):
            os.chmod(path, stat.S_IWRITE)
            func(path)

        shutil.rmtree(latest_dir, onexc=remove_readonly)
    shutil.copytree(run_dir, latest_dir)


def setup_logger(log_config: Dict[str, Any], output_dir: Path):
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, log_config.get("level", "INFO").upper()))
    
    if logger.hasHandlers():
        logger.handlers.clear()
        
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if log_config.get("log_to_file", False):
        fh = logging.FileHandler(output_dir / "pipeline_run.log")
        fh.setFormatter(formatter)
        logger.addHandler(fh)


def run_single_comparison(comp: dict, df_samples_train: pd.DataFrame, df_samples_test: pd.DataFrame, df_expression: pd.DataFrame, run_dir: Path, num_bootstraps: int, seed: int):
    exp_name = comp["experiment_name"]
    logging.info(f"[{exp_name}] Starting comparison: {comp['case_value']} vs {comp['control_value']}")
    
    comp_dir = run_dir / exp_name

    group_col = comp["group_column"]
    case_val = comp["case_value"]
    control_val = comp["control_value"]
    shuffle_labels = comp.get("shuffle_group_labels", False)

    df_samples_comp = df_samples_train.copy()

    if shuffle_labels:
        logging.info(f"[{exp_name}] Shuffling sample labels ({group_col})")
        df_samples_comp[group_col] = np.random.permutation(df_samples_comp[group_col].values)

    case_columns = df_samples_comp[df_samples_comp[group_col] == case_val].index
    control_columns = df_samples_comp[df_samples_comp[group_col] == control_val].index
    reference_val = comp.get("reference_value", None)
    
    if reference_val is not None:
        logging.info(f"[{exp_name}] Reference group specified: {group_col} = {reference_val}. This will be used for binarization.")
        reference_columns = df_samples_comp[df_samples_comp[group_col] == reference_val].index
    else:
        logging.info(f"[{exp_name}] No reference value specified. Binarization will be based on control group: {group_col} = {control_val}.")
        reference_columns = control_columns

    df_samples_comp = df_samples_comp.loc[case_columns.union(control_columns).union(reference_columns)]
    df_expression_comp = df_expression.loc[:, df_samples_comp.index]
        
    logging.info(f"[{exp_name}] Identified {len(case_columns)} cases, {len(reference_columns)} references, {len(control_columns)} controls.")

    logging.info(f"[{exp_name}] Binarizing expression data...")
    B_plus, B_minus = binarize_expression(df_expression, reference_columns)
    Bs = {"plus": B_plus, "minus": B_minus, "unified": B_plus | B_minus}

    for direction in Bs:
        logging.info(f"[{exp_name} | {direction}] Running CASh...")
        B = Bs[direction]
        B_case = B[case_columns]
        B_control = B[control_columns]
        
        df_results, beta_dist = run_CASh(
            B_case=B_case, 
            B_control=B_control, 
            b=num_bootstraps, 
            seed=seed,
            desc=f"{exp_name} | {direction}" 
        )

        comp_dir.mkdir(exist_ok=True, parents=True)
        logging.info(f"[{exp_name} | {direction}] Saving results to {comp_dir}...")
        save_experiment_results(
            df_results,
            beta_dist,
            B_case,
            B_control,
            df_samples_comp,
            df_samples_test,
            df_expression_comp,
            comp_dir,
            direction,
        )
        
    logging.info(f"[{exp_name}] Experiment complete.\n")


def run_pipeline(config: Union[str, Path, Dict[str, Any]]):
    cfg = load_configuration(config) if isinstance(config, (str, Path)) else config

    dataset = cfg["dataset"]
    base_data_dir = Path(cfg.get("data_dir", f"data/{dataset}"))
    base_out_dir = Path(cfg.get("output_dir", f"results/{dataset}"))
    
    run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = base_out_dir / run_timestamp
    run_dir.mkdir(parents=True, exist_ok=True)

    setup_logger(cfg.get("logging", {}), run_dir)
    logging.info(f"Starting pipeline for dataset: {dataset}")
    
    with open(run_dir / "config.json", "w") as file:
        json.dump(cfg, file, indent=4)

    logging.info(f"Loading expression data...")
    df_expression = pd.read_csv(base_data_dir / cfg["expression_data_file"], index_col=0)
    
    logging.info(f"Loading samples data...")
    df_samples = pd.read_csv(base_data_dir / cfg["samples_data_file"], index_col=0)

    num_bootstraps = cfg.get("num_bootstraps", 1000)
    seed = cfg.get("random_seed", 42)
    train_share = cfg.get("train_share", 1.0)
    df_samples_train, df_samples_test = split_samples(df_samples, train_share, seed)
    logging.info(f"Train share of {train_share}: Selected {len(df_samples_train)} training samples and {len(df_samples_test)} test samples.")
    comparisons = cfg.get("comparisons", [])

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for comp in comparisons:
            futures.append(
                executor.submit(
                    run_single_comparison, 
                    comp, 
                    df_samples_train,
                    df_samples_test,
                    df_expression, 
                    run_dir, 
                    num_bootstraps, 
                    seed
                )
            )
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                logging.error(f"A comparison generated an exception: {exc}")

    save_latest_run(run_dir, base_out_dir)
    logging.info(f"All experiments finished. Results saved to: {run_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CASh pipeline.")
    parser.add_argument(
        "config",
        nargs="?",
        default="data/example1/config.json",
        help="Path to the configuration JSON file."
    )
    args = parser.parse_args()

    config = args.config

    run_pipeline(config)