from typing import Optional, TypedDict, List, Union, Literal


class LoggingConfig(TypedDict):
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    log_to_file: bool
    log_dir: Optional[str]


class ExperimentConfig(TypedDict):
    name: str
    data_dir: Optional[str]
    output_dir: Optional[str]
    logging: Optional[LoggingConfig]
    