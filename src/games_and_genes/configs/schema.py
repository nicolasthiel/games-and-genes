from typing import Optional, TypedDict, List, Union, Literal, NotRequired, Any


class LoggingConfig(TypedDict):
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    log_to_file: bool
    log_dir: Optional[str]

class ExpressionDataConfig(TypedDict):
    expression_file: str
    expression_file_separator: Literal[',', '\t', ';', ' ']
    expression_id_col_idx: Optional[int]
    identifier_type: Optional[Literal['ensembl', 'gene_symbol', 'entrez']]

class SampleDataConfig(TypedDict):
    sample_data_file: str
    sample_data_file_separator: Literal[',', '\t', ';', ' ']
    sample_data_id_col_idx: Optional[int]

class DataConfig(TypedDict):
    data_dir: str
    expression_data: ExpressionDataConfig
    sample_data: SampleDataConfig

class ExpressionPreprocessingConfig(TypedDict):
    normalize: bool
    method_normalize: Optional[Literal['DSEQ2']]
    transform: bool
    method_transform: Optional[Literal['log2', 'log10']]

class FilterCondition(TypedDict):
    column: str
    operator: Literal['==', '!=', 'in', 'not in', '>', '<', '>=', '<=']
    value: Any  # Can be a string, number, or list of strings

class SamplePreprocessingConfig(TypedDict):
    keep_columns: NotRequired[List[str]] 
    drop_na_columns: NotRequired[List[str]]
    drop_duplicates: NotRequired[bool]
    filters: NotRequired[List[FilterCondition]]

class PreprocessingConfig(TypedDict):
    expression_preprocessing: ExpressionPreprocessingConfig
    sample_preprocessing: SamplePreprocessingConfig

class OutputConfig(TypedDict):
    output_dir: str
    overwrite: bool
    save_plots: bool
    save_intermediate_results: bool

class ExperimentConfig(TypedDict):
    name: str
    seed: int
    save_config: bool
    data: DataConfig
    output: OutputConfig
    preprocessing: PreprocessingConfig
    logging: LoggingConfig
    