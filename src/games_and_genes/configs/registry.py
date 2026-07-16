from .schema import ExperimentConfig


DEFAULT_CONFIGS = {
    "logging": {
        "level": "INFO",
        "log_to_file": True,
    },
}

EXAMPLE1: ExperimentConfig = {
    "name": "Example 1",
    "seed": 42,
    "save_config": False,
    "data": {
        "data_dir": "data/example1",
        "expression_data": {
            "expression_file": "expression.csv",
            "expression_file_separator": ",",
            "expression_id_col_idx": 0,
            "identifier_type": "ensembl"
        },
        "sample_data": {
            "sample_data_file": "sample.csv",
            "sample_data_file_separator": ",",
            "sample_data_id_col_idx": 0
        }
    },
    "output": {
        "output_dir": "out/example1",
        "overwrite": True,
        "save_plots": True,
        "save_intermediate_results": True,
    },
    "preprocessing": {
        "expression_preprocessing": {
            "normalize": True,
            "method_normalize": "DSEQ2",
            "transform": True,
            "method_transform": "log2"
        },
        "sample_preprocessing": {
            "drop_na_columns": ["condition"],
            "drop_duplicates": True,
            "filters": [
                {"column": "condition", "operator": "in", "value": ["control", "treated"]}
            ]
        }
    },
    "logging": DEFAULT_CONFIGS["logging"]
}


CONFIG_REGISTRY = {
    "example1": EXAMPLE1,
}