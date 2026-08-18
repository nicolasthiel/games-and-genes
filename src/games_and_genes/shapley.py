import logging
import numpy as np
import pandas as pd

from numpy.typing import NDArray
from tqdm import tqdm

from games_and_genes.configs.schema import ShapleyConfig


logger = logging.getLogger(__name__)


class ShapleyAnalyzer:

    def __init__(self, config: ShapleyConfig):
        self.discriminant_method_lower_bound = config.get("discriminant_method_lower_bound", 0)
        self.discriminant_method_upper_bound = config.get("discriminant_method_upper_bound", 100)
        logger.debug(
            "ShapleyAnalyzer initialized with bounds %s-%s",
            self.discriminant_method_lower_bound,
            self.discriminant_method_upper_bound,
        )
        
    
    def fit(self,expression_df: pd.DataFrame, reference_samples: list[str], diseased_samples: list[str]):
        logger.info(
            "Fitting Shapley analyzer with %d genes, %d reference samples, and %d diseased samples",
            expression_df.shape[0],
            len(reference_samples),
            len(diseased_samples),
        )
        self.expression_df = expression_df
        self.A = expression_df.to_numpy()
        self.SR = reference_samples
        self.SD = diseased_samples
        self.feature_names = expression_df.index.tolist()

        self.num_players = self.A.shape[0]

        col_idx_SR = [expression_df.columns.get_loc(s) for s in reference_samples]
        col_idx_SD = [expression_df.columns.get_loc(s) for s in diseased_samples]

        self.A_SR = self.A[:, col_idx_SR]
        self.A_SD = self.A[:, col_idx_SD]

        self.B = self.create_boolean_diff_expressed_matrix(
            A_SD=self.A_SD,
            A_SR=self.A_SR,
            lower_bound=self.discriminant_method_lower_bound,
            upper_bound=self.discriminant_method_upper_bound
        )
        self.boolean_expression_matrix = self.B
        self.sp_B = self.create_support_of_binary_matrix(self.B)
        self.coalitions = self.find_coalitions(self.sp_B)
        self.unanimity_coeffs = self.calculate_unanimity_coefficients(self.sp_B, self.coalitions)
        logger.debug("Prepared %d coalitions for Shapley calculation", len(self.coalitions))


    def get_artifacts(self) -> dict[str, object]:
        return {
            "expression_df": self.expression_df,
            "A": self.A,
            "A_SR": self.A_SR,
            "A_SD": self.A_SD,
            "boolean_expression_matrix": self.boolean_expression_matrix,
            "sp_B": self.sp_B,
            "coalitions": self.coalitions,
            "unanimity_coeffs": self.unanimity_coeffs,
            "feature_names": self.feature_names,
            "reference_samples": self.SR,
            "diseased_samples": self.SD,
        }


    def create_boolean_diff_expressed_matrix(
            self,
            A_SD: NDArray[np.float64],
            A_SR: NDArray[np.float64],
            lower_bound: float = 0,
            upper_bound: float = 100
        ) -> NDArray[np.bool_]:
        p_lower = np.percentile(A_SR, lower_bound, axis=1)
        p_upper = np.percentile(A_SR, upper_bound, axis=1)
        B = (A_SD <= p_lower[:, None]) | (A_SD >= p_upper[:, None])
        return B


    def create_support_of_binary_matrix(
            self,
            B: NDArray[np.bool_]
        ) -> list[set[int]]:
            return [set(np.nonzero(col)[0].tolist()) for col in B.T]


    def find_coalitions(self, sp_B: list[set[int]]) -> list[set[int]]:
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


    def calculate_unanimity_coefficients(self, sp_B: list[set[int]], coalitions: list[set[int]]) -> NDArray[np.float64]:
        return np.array([sp_B.count(coalition) for coalition in coalitions])/len(sp_B)


    def calculate_shapley_values(
            self,
            verbose: bool = True,
        ) -> NDArray[np.float64]:

        logger.debug("Calculating Shapley values for %d players", self.num_players)
        shapley_values = np.zeros(self.num_players, dtype=np.float64)
        for player in tqdm(range(self.num_players), disable=not verbose, desc="Calculating Shapley values"):
            for coalition_idx, coalition in enumerate(self.coalitions):
                if player in coalition:
                    shapley_values[player] += self.unanimity_coeffs[coalition_idx] / len(coalition)

        logger.info("Calculated Shapley values")
        return shapley_values
        
        