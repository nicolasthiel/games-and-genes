import argparse
import sys
from .registry import CONFIG_REGISTRY
from .schema import ExperimentConfig


def get_active_config() -> ExperimentConfig:
    # 1. Check for command line argument --config
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--config", type=str, default="debug", help="Name of experiment")
    args, _ = parser.parse_known_args()

    config_name = args.config

    # 2. Validate availability
    if config_name not in CONFIG_REGISTRY:
        print(f"Error: Config '{config_name}' not found.")
        print(f"Available: {list(CONFIG_REGISTRY.keys())}")
        sys.exit(1)

    print(f"Loaded Config: {config_name}")
    return CONFIG_REGISTRY[config_name]


# Create the singleton instance
active_config = get_active_config()