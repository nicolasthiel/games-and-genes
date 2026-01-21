import logging

from games_and_genes.configs import cfg
from games_and_genes.logger import setup_logging


logger = setup_logging(cfg['name'], cfg['logging'])
logger = logging.getLogger(__name__)

def run():
    logger.info(f"Running experiment: {cfg['name']}")
    logger.warning(f"Data directory: {cfg['data_dir']}")
    logger.error(f"Output directory: {cfg['output_dir']}")

if __name__ == "__main__":
    run()