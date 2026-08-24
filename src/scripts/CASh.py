import os
import json
import logging
import argparse
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from typing import Any, Dict, List, Union
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm


def binarize_expression(df_X : pd.DataFrame, ref_columns : List[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref_means = df_X[ref_columns].mean(axis=1)
    ref_stdevs = df_X[ref_columns].std(axis=1)
    
    over_thresholds = ref_means + ref_stdevs
    under_thresholds = ref_means - ref_stdevs
    
    B_plus = pd.DataFrame(False, index=df_X.index, columns=df_X.columns)
    B_minus = pd.DataFrame(False, index=df_X.index, columns=df_X.columns)
    
    for gene in df_X.index:
        B_plus.loc[gene] = (df_X.loc[gene] >= over_thresholds[gene])
        B_minus.loc[gene] = (df_X.loc[gene] <= under_thresholds[gene])
        
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


def run_CASh(B_case : pd.DataFrame, B_control : pd.DataFrame, b : int, seed : int = 42) -> tuple[pd.DataFrame, pd.DataFrame]:
    np.random.seed(seed)

    n, k = B_case.shape
    _, h = B_control.shape
    m = k + h 

    shapley_case_observed = calculate_shapley_values(B_case)
    shapley_control_observed = calculate_shapley_values(B_control)

    delta_observed = abs(shapley_case_observed - shapley_control_observed)

    count_beta_gte_delta = np.zeros(n)
    betas = np.zeros((n, b))
    for r in tqdm(range(b), desc="Running permutations"):
        idx_case_res = np.random.choice(k, size=k, replace=True)
        idx_control_res = np.random.choice(h, size=h, replace=True)
        B_case_res = B_case.iloc[:, idx_case_res]
        B_control_res = B_control.iloc[:, idx_control_res]

        shapley_case_res = calculate_shapley_values(B_case_res)
        shapley_control_res = calculate_shapley_values(B_control_res)

        beta = abs((shapley_case_observed - shapley_control_observed) - (shapley_case_res - shapley_control_res))
        count_beta_gte_delta += (beta >= delta_observed)
        betas[:, r] = beta

    raw_p_values = count_beta_gte_delta / b
    _, p_adjusted, _, _ = multipletests(raw_p_values, alpha=0.05, method='fdr_bh')

    beta_dist = pd.DataFrame(data=betas.T, index=[f"beta_{i+1}" for i in range(b)], columns=B_case.index)

    eps = 1e-12
    logfc = np.log2(shapley_case_observed + eps) - np.log2(shapley_control_observed + eps)

    results_df = pd.DataFrame({
            "shapley_case": shapley_case_observed,
            "shapley_control": shapley_control_observed,
            "observed_difference": shapley_case_observed - shapley_control_observed,
            "observed_abs_difference": delta_observed,
            "logFC": logfc,
            "pval": raw_p_values,
            "pval_adj": p_adjusted
        }, index=B_case.index)

    return results_df, beta_dist


def load_configuration(config_path: Union[str, Path]) -> Dict[str, Any]:
    """Loads the pipeline configuration from a JSON file."""
    with open(config_path, 'r') as file:
        return json.load(file)


def save_experiment_results(results_df: pd.DataFrame, beta_dist: pd.DataFrame, output_dir: Path, direction: str):
    """Saves the results and beta distribution to CSV files using Pathlib."""
    results_df.to_csv(output_dir / f"results_{direction}.csv")
    beta_dist.to_csv(output_dir / f"beta_dist_{direction}.csv")


def setup_logger(log_config: Dict[str, Any], output_dir: Path):
    """Configures the python logging module based on the JSON config."""
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

    for comp in cfg.get("comparisons", []):
        exp_name = comp["experiment_name"]
        logging.info(f"=== Starting Experiment: {exp_name} ===")
        
        comp_dir = run_dir / exp_name
        comp_dir.mkdir(exist_ok=True)

        group_col = comp["group_column"]
        case_val = comp["case_value"]
        control_val = comp["control_value"]

        case_columns = df_samples[df_samples[group_col] == case_val].index
        control_columns = df_samples[df_samples[group_col] == control_val].index
        logging.info(f"Identified {len(case_columns)} cases and {len(control_columns)} controls.")

        logging.info("Binarizing expression data...")
        B_plus, B_minus = binarize_expression(df_expression, control_columns)
        Bs = {"plus": B_plus, "minus": B_minus}

        for direction in Bs:
            logging.info(f"Running CASh ({direction} direction)...")
            B = Bs[direction]
            
            df_results, beta_dist = run_CASh(
                B_case=B[case_columns], 
                B_control=B[control_columns], 
                b=num_bootstraps, 
                seed=seed
            )
            
            save_experiment_results(df_results, beta_dist, comp_dir, direction)
            
        logging.info(f"Experiment {exp_name} complete.\n")

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