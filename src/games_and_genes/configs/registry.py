from .schema import ExperimentConfig


DEFAULT_CONFIGS = {
    "logging": {
        "level": "INFO",
        "log_to_file": False,
    },
}

EXAMPLE1: ExperimentConfig = {
    "name": "Example 1",
    "seed": 42,
    "data": {
        "data_dir": "data/example1",
        "expression_file": "norm.csv",
        "expression_id_column": "Gene",
        "sample_data_file": "sample.csv",
        "sample_data_id_column": "Sample",
    },
    "output": {
        "output_dir": "out/example1",
        "overwrite": True,
        "save_plots": True,
        "save_intermediate_results": True,
    },
    "preprocessing": {
        "normalize": True,
        "method_normalize": "DSEQ2",
        "transform": True,
        "method_transform": "log2",
    },
    "logging": DEFAULT_CONFIGS["logging"]
}


CONFIG_REGISTRY = {
    "example1": EXAMPLE1,
}