from typing import Optional, TypedDict, List, Union, Literal


class LoggingConfig(TypedDict):
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    log_to_file: bool
    log_dir: Optional[str]

class ExpressionDataConfig(TypedDict):
    expression_file: str
    expression_file_separator: Literal[',', '\t', ';', ' ']
    expression_id_col_idx: Optional[int]

class SampleDataConfig(TypedDict):
    sample_data_file: str
    sample_data_file_separator: Literal[',', '\t', ';', ' ']
    sample_data_id_col_idx: Optional[int]

class DataConfig(TypedDict):
    data_dir: str
    expression_data: ExpressionDataConfig
    sample_data: SampleDataConfig

class OutputConfig(TypedDict):
    output_dir: str
    overwrite: bool
    save_plots: bool
    save_intermediate_results: bool

class PreprocessingConfig(TypedDict):
    normalize: bool
    method_normalize: Optional[Literal['DSEQ2']]
    transform: bool
    method_transform: Optional[Literal['log2', 'log10']]

class ExperimentConfig(TypedDict):
    name: str
    seed: int
    save_config: bool
    data: DataConfig
    output: OutputConfig
    preprocessing: PreprocessingConfig
    logging: LoggingConfig
    