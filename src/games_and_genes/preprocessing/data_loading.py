import logging
import pandas as pd

from games_and_genes.configs.schema import DataConfig


logger = logging.getLogger(__name__)

class DataLoader:

    def __init__(self, config: DataConfig):
        self.data_dir = config.get("data_dir")
        self.expression_data_config = config.get("expression_data")
        self.sample_data_config = config.get("sample_data")

    def load_all_data(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        logger.info(f"Loading expression and sample data from {self.data_dir}.")
        expression_df = self.load_expression_data()
        sample_df = self.load_sample_data()
        return expression_df, sample_df

    def load_expression_data(self) -> pd.DataFrame:
        logger.info(f"Loading expression data {self.data_dir}/{self.expression_data_config.get('file_name')}.")
        expression_df = pd.read_csv(
            f"{self.data_dir}/{self.expression_data_config.get('file_name')}",
            sep=self.expression_data_config.get("file_separator"),
            index_col=self.expression_data_config.get("id_col_idx")
        )
        logger.info(f"Expression data loaded with shape {expression_df.shape}.")
        return expression_df
    
    def load_sample_data(self) -> pd.DataFrame:
        logger.info(f"Loading sample data {self.data_dir}/{self.sample_data_config.get('file_name')}.")
        sample_df = pd.read_csv(
            f"{self.data_dir}/{self.sample_data_config.get('file_name')}",
            sep=self.sample_data_config.get("file_separator"),
            index_col=self.sample_data_config.get("id_col_idx")
        )
        logger.info(f"Sample data loaded with shape {sample_df.shape}.")
        return sample_df