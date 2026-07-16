import logging
import pandas as pd
import numpy as np

from pydeseq2.dds import DeseqDataSet

from games_and_genes.configs.schema import SamplePreprocessingConfig


logger = logging.getLogger(__name__)

class SamplePreprocessor:

    def __init__(self, config: SamplePreprocessingConfig):
        self.drop_na_columns = config.get("drop_na_columns", [])
        self.drop_duplicates = config.get("drop_duplicates", True)
        self.filters = config.get("filters", [])
        logger.debug("Initialized SamplePreprocessor with preprocessing configuration.")


    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        logger.info("Starting preprocessing of sample data.")

        if self.drop_na_columns:
            df = df.dropna(subset=self.drop_na_columns)
            logger.debug(f"Dropped rows with NA in columns: {self.drop_na_columns}")

        if self.drop_duplicates:
            df = df.drop_duplicates()
            logger.debug("Dropped duplicate rows.")

        for filter_condition in self.filters:
            column = filter_condition['column']
            operator = filter_condition['operator']
            value = filter_condition['value']
            
            if operator == '==':
                df = df[df[column] == value]
            elif operator == '!=':
                df = df[df[column] != value]
            elif operator == '<':
                df = df[df[column] < value]
            elif operator == '<=':
                df = df[df[column] <= value]
            elif operator == '>':
                df = df[df[column] > value]
            elif operator == '>=':
                df = df[df[column] >= value]
            elif operator == 'in':
                df = df[df[column].isin(value)]
            elif operator == 'not in':
                df = df[~df[column].isin(value)]
            else:
                logger.error(f"Unsupported operator: {operator}")
                raise ValueError(f"Unsupported operator: {operator}")
            
            logger.debug(f"Applied filter: {column} {operator} {value}")

        logger.info("Completed preprocessing of sample data.")
        return df