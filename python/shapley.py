import numpy as np
from numpy.typing import NDArray

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
    ):
        return [set(np.nonzero(col)[0].tolist()) for col in B.T]
