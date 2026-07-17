import logging
import numpy as np

from games_and_genes.configs import cfg
from games_and_genes.logger import setup_logging
from games_and_genes.pipeline import ExperimentPipeline


setup_logging(cfg['logging'], cfg['output']['output_dir'])
logger = logging.getLogger(__name__)

def main():
    try:
        logger.info("Starting experiment %s", cfg['name'])
        logger.debug("Data directory: %s", cfg['data']['data_dir'])
        logger.debug("Output directory: %s", cfg['output']['output_dir'])

        if cfg['output'].get('overwrite'):
            logger.warning("Output directory may be overwritten: %s", cfg['output']['output_dir'])

        logger.info("Setting global random seed to %s", cfg['seed'])
        np.random.seed(cfg['seed'])

        pipeline = ExperimentPipeline(config=cfg)
        pipeline.run()
        logger.info("Experiment %s completed successfully", cfg['name'])
    except Exception:
        logger.critical("Experiment %s failed.", cfg['name'], exc_info=True)
        raise

if __name__ == "__main__":
    main()