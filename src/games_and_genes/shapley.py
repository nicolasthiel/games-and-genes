import numpy as np
from numpy.typing import NDArray


def shapley_value(num_players: int, coalitions: list[set[int]], unanimity_coeffs: NDArray[np.float64]) -> NDArray[np.float64]:

    shapley_values = np.zeros(num_players, dtype=np.float64)

    for player in range(num_players):
        for coalition_idx, coalition in enumerate(coalitions):
            if player in coalition:
                shapley_values[player] += unanimity_coeffs[coalition_idx] / len(coalition)

    return shapley_values