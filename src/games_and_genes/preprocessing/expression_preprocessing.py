import logging
import pandas as pd
import numpy as np

from pydeseq2.dds import DeseqDataSet
from dataclasses import dataclass

from games_and_genes.configs.schema import ExpressionPreprocessingConfig


logger = logging.getLogger(__name__)

class ExpressionPreprocessor:

    def __init__(self, config: ExpressionPreprocessingConfig):
        self.normalize = config.get("normalize")
        self.method_normalize = config.get("method_normalize")
        self.transform = config.get("transform")
        self.method_transform = config.get("method_transform")
        logger.debug(
            "ExpressionPreprocessor initialized with normalize=%s, transform=%s",
            self.normalize,
            self.transform,
        )
    

    def preprocess(self, expression_df: pd.DataFrame, sample_ids_to_keep: list) -> pd.DataFrame:
        logger.info("Preprocessing expression data")

        logger.debug("Keeping %d of %d expression samples", len(sample_ids_to_keep), expression_df.shape[1])
        expression_df = expression_df.loc[:, expression_df.columns.isin(sample_ids_to_keep)]
        logger.debug("Expression data shape after sample filter: %s", expression_df.shape)

        if self.normalize:
            logger.debug("Applying %s normalization", self.method_normalize)
            expression_df = self._normalize(expression_df)

        if self.transform:
            logger.debug("Applying %s transform", self.method_transform)
            expression_df = self._transform(expression_df)

        if not self.normalize and not self.transform:
            logger.debug("No expression preprocessing steps enabled")

        return expression_df


    def _normalize(self, df: pd.DataFrame) -> pd.DataFrame:

        def normalize_deseq2(df: pd.DataFrame) -> pd.DataFrame:
            logger.error("DSEQ2 normalization is not implemented yet")
            raise NotImplementedError("DSEQ2 normalization is not implemented yet.")
        
        def normalize_tpm(df: pd.DataFrame) -> pd.DataFrame:
            logger.error("TPM normalization is not implemented yet")
            raise NotImplementedError("TPM normalization is not implemented yet.")

        if self.method_normalize == "DSEQ2":
            return normalize_deseq2(df)
        elif self.method_normalize == "TPM":
            return normalize_tpm(df)
        logger.error("Unsupported normalization method %s", self.method_normalize)
        raise ValueError(f"Unsupported normalization method: {self.method_normalize}")


    def _transform(self, df: pd.DataFrame) -> pd.DataFrame:
        logger.error("Expression transformation is not implemented for method %s", self.method_transform)
        raise NotImplementedError("Data transformation is not implemented yet.")
        


    