from typing import TypedDict, List, Union, Literal


class LoggingConfig(TypedDict):
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    log_to_file: bool
    log_file_path: str


class ExperimentConfig(TypedDict):
    name: str
    data_dir: str
    output_dir: str
    logging: LoggingConfig
    