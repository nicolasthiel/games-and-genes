import os
import pandas as pd
import numpy as np
import json
import logging

from matplotlib import pyplot as plt
from typing import Any, Dict, Dict, List, Union, Union
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


def load_configuration(config_path: str) -> Dict[str, Any]:
    """Loads the pipeline configuration from a JSON file."""
    with open(config_path, 'r') as file:
        return json.load(file)


def save_experiment_results(results_df: pd.DataFrame, output_dir: str, direction: str):
    """Saves the results and beta distribution to CSV files."""
    results_path = os.path.join(output_dir, f"results_{direction}.csv")
    results_df.to_csv(results_path)


def run_pipeline(config: Union[str, Dict[str, Any]]):

    # Load configuration
    if isinstance(config, str):
        cfg = load_configuration(config)
    else:
        cfg = config

    # Ensure output directory exists
    output_dir = f"out/{cfg['dataset']}/{cfg['experiment_name']}/"
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    df_expression = pd.read_csv(f"data/{cfg['dataset']}/{cfg['expression_data_file']}", index_col=0)
    df_samples = pd.read_csv(f"data/{cfg['dataset']}/{cfg['samples_data_file']}", index_col=0)

    # Identify case and control columns based on the configuration
    case_columns = df_samples[df_samples[cfg["group_column"]] == cfg["case_value"]].index
    control_columns = df_samples[df_samples[cfg["group_column"]] == cfg["control_value"]].index

    # Binarize expression data
    B_plus, B_minus = binarize_expression(df_expression, control_columns)
    Bs = {"plus": B_plus, "minus": B_minus}

    for direction in Bs:
        B = Bs[direction]
        B_case = B[case_columns]
        B_control = B[control_columns]

        df_results, beta_dist = run_CASh(B_case, B_control, b=cfg["num_bootstraps"], seed=cfg.get("random_seed", 42))
        save_experiment_results(df_results, output_dir, direction)

    # Save configuration for reproducibility
    config_path = os.path.join(output_dir, "config.json")
    with open(config_path, "w") as file:
        json.dump(cfg, file, indent=4)


if __name__ == "__main__":

    # Example configuration 
    config = {
        "experiment_name": "healthy_vs_diseased",
        "dataset": "GSE42568",
        "expression_data_file": "processed/matrix_final.csv",
        "samples_data_file": "processed/samples_filtered.csv",
        "group_column": "tissue.ch1",
        "case_value": "breast cancer",
        "control_value": "normal breast",
        "num_bootstraps": 1000,
        "random_seed": 42,
        "logging": {
                "level": "DEBUG",
                "log_to_file": True,
        }
    }
    run_pipeline(config)