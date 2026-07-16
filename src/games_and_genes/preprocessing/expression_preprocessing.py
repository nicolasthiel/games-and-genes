import logging
import pandas as pd
import numpy as np

from pydeseq2.dds import DeseqDataSet

from games_and_genes.configs.schema import ExpressionPreprocessingConfig


logger = logging.getLogger(__name__)

class ExpressionPreprocessor:

    def __init__(self, config: ExpressionPreprocessingConfig):
        self.normalize = config.get("normalize")
        self.method_normalize = config.get("method_normalize")
        self.transform = config.get("transform")
        self.method_transform = config.get("method_transform")
        logger.debug("Initialized ExpressionPreprocessor with preprocessing configuration.")
    

    def preprocess(self, expression_df: pd.DataFrame, sample_ids_to_keep: list) -> pd.DataFrame:
        logger.info("Starting preprocessing of expression data.")

        expression_df = expression_df.loc[:, expression_df.columns.isin(sample_ids_to_keep)]

        if self.normalize:
            expression_df = self._normalize(expression_df)

        if self.transform:
            expression_df = self._transform(expression_df)

        return expression_df


    def _normalize(self, df: pd.DataFrame) -> pd.DataFrame:

        def normalize_deseq2(df: pd.DataFrame) -> pd.DataFrame:
            raise NotImplementedError("DSEQ2 normalization is not implemented yet.")
        
        def normalize_tpm(df: pd.DataFrame) -> pd.DataFrame:
            raise NotImplementedError("TPM normalization is not implemented yet.")

        if self.method_normalize == "DSEQ2":
            return normalize_deseq2(df)
        elif self.method_normalize == "TPM":
            return normalize_tpm(df)


    def _transform(self, df: pd.DataFrame) -> pd.DataFrame:

        raise NotImplementedError("Data transformation is not implemented yet.")
        


    