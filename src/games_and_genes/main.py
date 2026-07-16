import logging
import numpy as np

from games_and_genes.configs import cfg
from games_and_genes.logger import setup_logging
from games_and_genes.pipeline import ExperimentPipeline


logger = setup_logging(cfg['logging'], cfg['output']['output_dir'])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Running experiment: {cfg['name']}")
    logger.warning(f"Data directory: {cfg['data']['data_dir']}")
    logger.error(f"Output directory: {cfg['output']['output_dir']}")

    logger.info(f"Setting global random seed to {cfg['seed']}")
    np.random.seed(cfg['seed'])

    pipeline = ExperimentPipeline(config=cfg)
    pipeline.run()

if __name__ == "__main__":
    main()