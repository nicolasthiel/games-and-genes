from .schema import ExperimentConfig


EXAMPLE1: ExperimentConfig = {
    "name": "Experiment 1",
    "data_dir": "data/example1",
    "output_dir": "out/example1",
}

DEBUG: ExperimentConfig = {
    "name": "Debug Experiment",
    "data_dir": "data/debug",
    "output_dir": "out/debug",
}


CONFIG_REGISTRY = {
    "debug": DEBUG,
    "example1": EXAMPLE1,
}