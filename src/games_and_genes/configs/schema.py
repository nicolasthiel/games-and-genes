from typing import Optional, TypedDict, List, Union, Literal


class LoggingConfig(TypedDict):
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    log_to_file: bool
    log_dir: Optional[str]

class DataConfig(TypedDict):
    data_dir: str
    expression_file: str
    expression_id_column: Optional[str]
    sample_data_file: str
    sample_data_id_column: Optional[str]

class OutputConfig(TypedDict):
    output_dir: str
    overwrite: bool

class PreprocessingConfig(TypedDict):
    normalize: bool
    method_normalize: Optional[Literal['DSEQ2']]

class ExperimentConfig(TypedDict):
    name: str
    data: DataConfig
    output: OutputConfig
    preprocessing: PreprocessingConfig
    logging: LoggingConfig
    