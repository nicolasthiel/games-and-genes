from .schema import ExperimentConfig


DEFAULT_CONFIGS = {
    "logging": {
        "level": "DEBUG",
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
            "identifier_type": "custom"
        },
        "sample_data": {
            "sample_data_file": "sample.csv",
            "sample_data_file_separator": ",",
            "sample_data_id_col_idx": 0,
            "sample_data_condition_col_name": "CONDITION",
            "sample_data_condition_control_value": "reference",
            "sample_data_condition_case_value": "diseased"
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
            "normalize": False,
            "method_normalize": "DSEQ2",
            "transform": False,
            "method_transform": "log2"
        },
        "sample_preprocessing": {
            "drop_na_columns": [],
            "drop_duplicates": False,
            "filters": []
        }
    },
    "shapley": {
        "discriminant_method_lower_bound": 0,
        "discriminant_method_upper_bound": 100
    },
    "logging": DEFAULT_CONFIGS["logging"]
}

EXAMPLE2: ExperimentConfig = {
    "name": "Example 2",
    "seed": 42,
    "save_config": False,
    "data": {
        "data_dir": "data/example2",
        "expression_data": {
            "expression_file": "expression.csv",
            "expression_file_separator": ",",
            "expression_id_col_idx": 0,
            "identifier_type": "custom"
        },
        "sample_data": {
            "sample_data_file": "sample.csv",
            "sample_data_file_separator": ",",
            "sample_data_id_col_idx": 0,
            "sample_data_condition_col_name": "CONDITION",
            "sample_data_condition_control_value": "reference",
            "sample_data_condition_case_value": "diseased"
        }
    },
    "output": {
        "output_dir": "out/example2",
        "overwrite": True,
        "save_plots": True,
        "save_intermediate_results": True,
    },
    "preprocessing": {
        "expression_preprocessing": {
            "normalize": False,
            "method_normalize": "DSEQ2",
            "transform": False,
            "method_transform": "log2"
        },
        "sample_preprocessing": {
            "drop_na_columns": [],
            "drop_duplicates": False,
            "filters": []
        }
    },
    "shapley": {
        "discriminant_method_lower_bound": 0,
        "discriminant_method_upper_bound": 100
    },
    "logging": DEFAULT_CONFIGS["logging"]
}


CONFIG_REGISTRY = {
    "example1": EXAMPLE1,
    "example2": EXAMPLE2
}