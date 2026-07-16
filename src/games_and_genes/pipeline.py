import logging
from games_and_genes.configs.schema import ExperimentConfig
from games_and_genes.data_loading import DataLoader
from games_and_genes.preprocessing.expression_preprocessing import ExpressionPreprocessor
from games_and_genes.preprocessing.sample_preprocessing import SamplePreprocessor

from games_and_genes.shapley import *

logger = logging.getLogger(__name__)


class ExperimentPipeline:

    def __init__(self, config: ExperimentConfig):
        self.config = config
        self.data_loader = DataLoader(config.get("data"))
        self.sample_preprocessor = SamplePreprocessor(config.get("preprocessing").get("sample_preprocessing"))
        self.expression_preprocessor = ExpressionPreprocessor(config.get("preprocessing").get("expression_preprocessing"))


    def run(self):
        logger.info("Starting the experiment pipeline.")
        
        # Data loading
        expression_df, sample_df = self.data_loader.load_data()

        # Data preprocessing
        sample_df = self.sample_preprocessor.preprocess(sample_df)
        sample_ids_to_keep = sample_df.index.tolist()

        expression_df = self.expression_preprocessor.preprocess(expression_df, sample_ids_to_keep)

        # Shapley value calculation
        B = create_boolean_diff_expressed_matrix(expression_df, sample_df)
        sp_B = create_support_of_binary_matrix(B)
        coalitions = find_coalitions(sp_B)
        unanimity_coeffs = calculate_unanimity_coefficients(sp_B, coalitions)
        shapley_values = calculate_shapley_value(
            num_players=expression_df.shape[0],
            coalitions=coalitions,
            unanimity_coeffs=unanimity_coeffs
        )
        

        # Should we do statistical DEG here?

        # Should we plot here?
        
        logger.info("Experiment pipeline completed successfully.")