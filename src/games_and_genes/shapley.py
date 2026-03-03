import numpy as np

from numpy.typing import NDArray
from tqdm import tqdm


def boolean_diff_expressed_matrix(
    A_SD: NDArray[np.number],
    A_SR: NDArray[np.number],
    lower_prctl: int = 0,
    upper_prctl: int = 100
) -> NDArray[np.bool_]:
    p_lower = np.percentile(A_SR, lower_prctl, axis=1)
    p_upper = np.percentile(A_SR, upper_prctl, axis=1)
    B = (A_SD <= p_lower[:, None]) | (A_SD >= p_upper[:, None])
    return B


def support_of_binary_matrix(
        B: NDArray[np.bool_]
    ) -> list[set[int]]:
        return [set(np.nonzero(col)[0].tolist()) for col in B.T]


def find_coalitions(sp_B: list[set[int]]) -> list[set[int]]:
    coalitions = []
    if not sp_B:
        return coalitions
    
    for support in sp_B:
        
        if not support: # ignore empty set
            continue
            
        to_append = True
        for coalition in coalitions:
            if coalition == support:
                to_append = False
                break
        
        if to_append:
            coalitions.append(support)
            
    return coalitions


def unanimity_coefficients(sp_B: list[set[int]], coalitions: list[set[int]]) -> NDArray[np.float64]:
    return np.array([sp_B.count(coalition) for coalition in coalitions])/len(sp_B)


def calculate_shapley_value(
        num_players: int,
        coalitions: list[set[int]],
        unanimity_coeffs: NDArray[np.float64],
        verbose: bool = True,
    ) -> NDArray[np.float64]:

    shapley_values = np.zeros(num_players, dtype=np.float64)
    for player in tqdm(range(num_players), disable=not verbose, desc="Calculating Shapley values"):
        for coalition_idx, coalition in enumerate(coalitions):
            if player in coalition:
                shapley_values[player] += unanimity_coeffs[coalition_idx] / len(coalition)

    return shapley_values