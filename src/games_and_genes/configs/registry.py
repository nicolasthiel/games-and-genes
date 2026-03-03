from .schema import ExperimentConfig


EXAMPLE1: ExperimentConfig = {
    "name": "Example 1",
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
    },
    "preprocessing": {
        "normalize": True,
        "method_normalize": "DSEQ2",
    },
    "logging": {
        "level": "INFO",
        "log_to_file": True,
        "log_dir": "logs/example1",
    },
}


CONFIG_REGISTRY = {
    "example1": EXAMPLE1,
}