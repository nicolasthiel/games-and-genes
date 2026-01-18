from typing import TypedDict, List, Union, Literal


class ExperimentConfig(TypedDict):
    name: str
    data_dir: str
    output_dir: str
    