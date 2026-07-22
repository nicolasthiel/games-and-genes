import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def ensure_parent_dir(file_path: str | Path) -> Path:
	path = Path(file_path)
	path.parent.mkdir(parents=True, exist_ok=True)
	return path


def _check_overwrite(path: Path, overwrite: bool) -> None:
	if path.exists() and not overwrite:
		raise FileExistsError(f"Refusing to overwrite existing file: {path}")


def _json_default(value: Any) -> Any:
	if isinstance(value, set):
		return sorted(value)
	if isinstance(value, np.ndarray):
		return value.tolist()
	if isinstance(value, Path):
		return str(value)
	raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def save_dataframe(dataframe: pd.DataFrame, file_path: str | Path, overwrite: bool = True) -> Path:
	path = ensure_parent_dir(file_path)
	_check_overwrite(path, overwrite)
	dataframe.to_csv(path)
	return path


def save_series(series: pd.Series, file_path: str | Path, overwrite: bool = True, column_name: str | None = None) -> Path:
	path = ensure_parent_dir(file_path)
	_check_overwrite(path, overwrite)
	if column_name is not None:
		series = series.rename(column_name)
	series.to_csv(path, header=True)
	return path


def save_ndarray(array: np.ndarray, file_path: str | Path, overwrite: bool = True, index: pd.Index | None = None, columns: list[str] | pd.Index | None = None) -> Path:
	dataframe = pd.DataFrame(array, index=index, columns=columns)
	return save_dataframe(dataframe, file_path, overwrite=overwrite)


def save_json(data: Any, file_path: str | Path, overwrite: bool = True) -> Path:
	path = ensure_parent_dir(file_path)
	_check_overwrite(path, overwrite)
	with path.open("w", encoding="utf-8") as handle:
		json.dump(data, handle, indent=2, default=_json_default)
		handle.write("\n")
	return path
