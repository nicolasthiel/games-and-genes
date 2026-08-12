import pandas as pd
import numpy as np

from typing import List, Set
from numpy.typing import NDArray
from tqdm import tqdm


def run_CASh(
    B1: NDArray[np.bool_],
    B2: NDArray[np.bool_],
    b: int = 1000
) -> None:

    sp_B1 = create_support_of_binary_matrix(B1)
    sp_B2 = create_support_of_binary_matrix(B2)

    coalitions1 = find_coalitions(sp_B1)
    coalitions2 = find_coalitions(sp_B2)

    unanimity_coeffs1 = calculate_unanimity_coefficients(sp_B1, coalitions1)
    unanimity_coeffs2 = calculate_unanimity_coefficients(sp_B2, coalitions2)

    shapley_values1 = calculate_shapley_values(B1.shape[0], coalitions1, unanimity_coeffs1)
    shapley_values2 = calculate_shapley_values(B2.shape[0], coalitions2, unanimity_coeffs2)

    delta_observed = np.abs(shapley_values2 - shapley_values1)
    print("Observed delta:", delta_observed)
    k = B1.shape[1]
    h = B2.shape[1]
    beta = []
    for r in tqdm(range(b), desc="Sampling"):
        sample_idx1 = np.random.choice(k, size=k, replace=True)
        B1r = B1[:, sample_idx1]
        sample_idx2 = np.random.choice(h, size=h, replace=True)
        B2r = B2[:, sample_idx2]

        sp_B1r = create_support_of_binary_matrix(B1r)
        sp_B2r = create_support_of_binary_matrix(B2r)

        coalitions1r = find_coalitions(sp_B1r)
        coalitions2r = find_coalitions(sp_B2r)

        unanimity_coeffs1r = calculate_unanimity_coefficients(sp_B1r, coalitions1r)
        unanimity_coeffs2r = calculate_unanimity_coefficients(sp_B2r, coalitions2r)
    
        shapley_values1r = calculate_shapley_values(B1r.shape[0], coalitions1r, unanimity_coeffs1r)
        shapley_values2r = calculate_shapley_values(B2r.shape[0], coalitions2r, unanimity_coeffs2r)

        delta_r = np.abs(shapley_values2r - shapley_values1r)
        beta_r = np.abs(delta_r - delta_observed)
        beta.append(beta_r)

    beta_arr = np.vstack(beta)
    pvals = np.mean(beta_arr >= delta_observed[None, :], axis=0)

    return pvals, beta_arr


def create_boolean_diff_expressed_matrix(
    A_SD: NDArray[np.float64],
    A_SR: NDArray[np.float64],
    direction: str = "both",
) -> NDArray[np.bool_]:

    mean_A_SR = np.mean(A_SR, axis=1)
    std_A_SR = np.std(A_SR, axis=1)

    if direction == "both":
        p_lower = mean_A_SR - std_A_SR
        p_upper = mean_A_SR + std_A_SR
    elif direction == "up":
        p_lower = mean_A_SR - std_A_SR
        p_upper = np.inf
    elif direction == "down":
        p_lower = -np.inf
        p_upper = mean_A_SR + std_A_SR

    B = (A_SD <= p_lower[:, None]) | (A_SD >= p_upper[:, None])
    return B


def create_support_of_binary_matrix(B: NDArray[np.bool_]) -> List[Set[int]]:

    return [set(np.nonzero(col)[0].tolist()) for col in B.T]


def find_coalitions(sp_B: List[Set[int]]) -> List[Set[int]]:

    coalitions: List[Set[int]] = []
    if not sp_B:
        return coalitions

    for support in sp_B:
        if not support:
            continue

        if support not in coalitions:
            coalitions.append(support)

    return coalitions


def calculate_unanimity_coefficients(
    sp_B: List[Set[int]], coalitions: List[Set[int]]
) -> NDArray[np.float64]:

    if not sp_B:
        return np.array([], dtype=np.float64)

    counts = np.array([sp_B.count(coalition) for coalition in coalitions], dtype=np.float64)
    return counts / len(sp_B)


def calculate_shapley_values(
    num_players: int,
    coalitions: List[Set[int]],
    unanimity_coeffs: NDArray[np.float64],
    verbose: bool = False,
) -> NDArray[np.float64]:
    
    shapley_values = np.zeros(num_players, dtype=np.float64)
    for player in tqdm(range(num_players), disable=not verbose, desc="Calculating Shapley values"):
        for coalition_idx, coalition in enumerate(coalitions):
            if player in coalition:
                shapley_values[player] += unanimity_coeffs[coalition_idx] / len(coalition)

    return shapley_values


if __name__ == "__main__":

    A_SR = pd.read_csv("out/GSE42568/A_SR.csv", index_col=0).to_numpy()
    A_SD = pd.read_csv("out/GSE42568/A_SD.csv", index_col=0).to_numpy()

    B1 = create_boolean_diff_expressed_matrix(A_SR, A_SR, direction="both")
    print("B1:", B1)
    B2 = create_boolean_diff_expressed_matrix(A_SD, A_SR, direction="both")
    print("B2:", B2)
    pvals, beta_arr = run_CASh(B1, B2)
    print("P-values:", pvals)