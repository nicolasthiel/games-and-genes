from typing import Optional, TypedDict, List, Union, Literal


class LoggingConfig(TypedDict):
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    log_to_file: bool
    log_dir: Optional[str]

class DataConfig(TypedDict):
    data_dir: str
    expression_file: str
    expression_id_col_idx: Optional[int]
    sample_data_file: str
    sample_data_id_col_idx: Optional[int]

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
    data: DataConfig
    output: OutputConfig
    preprocessing: PreprocessingConfig
    logging: LoggingConfig
    