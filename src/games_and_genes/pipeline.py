import logging
from games_and_genes.configs.schema import ExperimentConfig
from games_and_genes.data_loading import DataLoader
from games_and_genes.preprocessing.expression_preprocessing import ExpressionPreprocessor

logger = logging.getLogger(__name__)


class ExperimentPipeline:

    def __init__(self, config: ExperimentConfig):
        self.config = config
        self.data_loader = DataLoader(config.get("data"))
        self.expression_preprocessor = ExpressionPreprocessor(config.get("preprocessing").get("expression"))


    def run(self):
        logger.info("Starting the experiment pipeline.")
        
        
        
        logger.info("Experiment pipeline completed successfully.")