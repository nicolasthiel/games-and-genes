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
        logger.debug(
            "SamplePreprocessor initialized with %d filters, drop_duplicates=%s",
            len(self.filters),
            self.drop_duplicates,
        )


    def preprocess(self, df: pd.DataFrame) -> pd.DataFrame:
        logger.info("Preprocessing sample data")

        if self.drop_na_columns:
            before_rows = len(df)
            df = df.dropna(subset=self.drop_na_columns)
            logger.debug("Dropped %d rows with NA in %s", before_rows - len(df), self.drop_na_columns)

        if self.drop_duplicates:
            before_rows = len(df)
            df = df.drop_duplicates()
            logger.debug("Dropped %d duplicate rows", before_rows - len(df))

        for filter_condition in self.filters:
            column = filter_condition['column']
            operator = filter_condition['operator']
            value = filter_condition['value']
            before_rows = len(df)
            
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
                logger.error("Unsupported filter operator %s", operator)
                raise ValueError(f"Unsupported operator: {operator}")
            
            logger.debug("Applied filter %s %s %r; rows %d -> %d", column, operator, value, before_rows, len(df))

            if df.empty:
                logger.warning("Sample filter %s %s %r removed all rows", column, operator, value)

        logger.info("Completed preprocessing of sample data")
        return df