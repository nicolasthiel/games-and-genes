import pandas as pd
import numpy as np
import logging
import pathlib

from games_and_genes.configs.schema import DataConfig


logger = logging.getLogger(__name__)

class DataLoader:

    def __init__(self, config: DataConfig):
        self.data_dir = pathlib.Path(config['data_dir'])

        self.expression_file_path = self.data_dir / config['expression_data']['expression_file']
        self.expression_file_separator = config['expression_data']['expression_file_separator']
        self.expression_id_col_idx = config['expression_data'].get('expression_id_col_idx', 0)

        self.sample_data_file_path = self.data_dir / config['sample_data']['sample_data_file']
        self.sample_data_file_separator = config['sample_data']['sample_data_file_separator']
        self.sample_data_id_col_idx = config['sample_data'].get('sample_data_id_col_idx', 0)
        logger.debug("DataLoader initialized for %s", self.data_dir)


    def load_data(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        logger.info("Loading expression and sample data")
        expression_df = self._load_expression_data()
        sample_df = self._load_sample_data()
        logger.info("Loaded expression and sample data")
        return expression_df, sample_df


    def _load_expression_data(self) -> pd.DataFrame:
        file_extension = self.expression_file_path.suffix.lower()
        if file_extension == '.csv':
            df = pd.read_csv(
                self.expression_file_path,
                sep=self.expression_file_separator,
                index_col=self.expression_id_col_idx
            )
        elif file_extension in ['.tsv', '.txt']:
            df = pd.read_table(
                self.expression_file_path,
                sep=self.expression_file_separator,
                index_col=self.expression_id_col_idx
            )
        else:
            logger.error("Unsupported expression data format %s for %s", file_extension, self.expression_file_path)
            raise ValueError(f"Unsupported file format: {file_extension}")
        logger.debug("Expression data shape %s from %s", df.shape, self.expression_file_path)
        return df


    def _load_sample_data(self) -> pd.DataFrame:
        file_extension = self.sample_data_file_path.suffix.lower()
        if file_extension == '.csv':
            df = pd.read_csv(
                self.sample_data_file_path,
                sep=self.sample_data_file_separator,
                index_col=self.sample_data_id_col_idx
            )
        elif file_extension in ['.tsv', '.txt']:
            df = pd.read_table(
                self.sample_data_file_path,
                sep=self.sample_data_file_separator,
                index_col=self.sample_data_id_col_idx
            )
        else:
            logger.error("Unsupported sample data format %s for %s", file_extension, self.sample_data_file_path)
            raise ValueError(f"Unsupported file format: {file_extension}")
        logger.debug("Sample data shape %s from %s", df.shape, self.sample_data_file_path)
        return df
    