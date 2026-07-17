import logging
from games_and_genes.configs.schema import ExperimentConfig
from games_and_genes.data_loading import DataLoader
from games_and_genes.preprocessing.expression_preprocessing import ExpressionPreprocessor
from games_and_genes.preprocessing.sample_preprocessing import SamplePreprocessor
from games_and_genes.shapley import ShapleyAnalyzer

from games_and_genes.shapley import *

logger = logging.getLogger(__name__)


class ExperimentPipeline:

    def __init__(self, config: ExperimentConfig):
        self.config = config
        self.data_loader = DataLoader(config.get("data"))
        self.sample_preprocessor = SamplePreprocessor(config.get("preprocessing").get("sample_preprocessing"))
        self.expression_preprocessor = ExpressionPreprocessor(config.get("preprocessing").get("expression_preprocessing"))
        self.shapley_analyzer = ShapleyAnalyzer(config.get("shapley"))


    def run(self):
        logger.info("Starting experiment pipeline")
        
        # Data loading
        expression_df, sample_df = self.data_loader.load_data()

        # Data preprocessing
        sample_df = self.sample_preprocessor.preprocess(sample_df)
        sample_ids_to_keep = sample_df.index.tolist()

        expression_df = self.expression_preprocessor.preprocess(expression_df, sample_ids_to_keep)
        self.shapley_analyzer.fit(
            expression_df,
            reference_samples=sample_df[sample_df[self.config["data"]["sample_data"]["sample_data_condition_col_name"]] == self.config["data"]["sample_data"]["sample_data_condition_control_value"]].index.tolist(),
            diseased_samples=sample_df[sample_df[self.config["data"]["sample_data"]["sample_data_condition_col_name"]] == self.config["data"]["sample_data"]["sample_data_condition_case_value"]].index.tolist()
        )
        shapley_values = self.shapley_analyzer.calculate_shapley_values()
        
        # Should we do statistical DEG here?

        # Should we plot here?
        
    logger.info("Experiment pipeline completed")