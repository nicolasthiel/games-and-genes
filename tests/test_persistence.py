import json
import sys
import tempfile
from pathlib import Path
import types
import unittest

import numpy as np
import pandas as pd

fake_pydeseq2 = types.ModuleType("pydeseq2")
fake_pydeseq2_dds = types.ModuleType("pydeseq2.dds")
fake_pydeseq2_dds.DeseqDataSet = object
fake_pydeseq2.dds = fake_pydeseq2_dds
sys.modules.setdefault("pydeseq2", fake_pydeseq2)
sys.modules.setdefault("pydeseq2.dds", fake_pydeseq2_dds)
sys.argv = [sys.argv[0], "--config", "example1"]

from games_and_genes.pipeline import ExperimentPipeline
from games_and_genes.shapley import ShapleyAnalyzer
from games_and_genes.utils import save_dataframe, save_json, save_ndarray, save_series


class TestPersistenceHelpers(unittest.TestCase):

    def test_save_dataframe_and_series(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            base_path = Path(tmp_dir)
            dataframe = pd.DataFrame({"sample_a": [1, 2], "sample_b": [3, 4]}, index=["gene_1", "gene_2"])
            series = pd.Series([0.1, 0.2], index=["gene_1", "gene_2"], name="shapley_value")

            dataframe_path = save_dataframe(dataframe, base_path / "frame.csv")
            series_path = save_series(series, base_path / "series.csv")

            self.assertTrue(dataframe_path.exists())
            self.assertTrue(series_path.exists())

            loaded_frame = pd.read_csv(dataframe_path, index_col=0)
            loaded_series = pd.read_csv(series_path, index_col=0).iloc[:, 0]

            pd.testing.assert_frame_equal(loaded_frame, dataframe)
            pd.testing.assert_series_equal(loaded_series, series, check_names=False)

    def test_save_json_and_ndarray(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            base_path = Path(tmp_dir)
            array = np.array([[True, False], [False, True]])
            payload = {"sp_B": [{0, 1}, {1}]}

            array_path = save_ndarray(array, base_path / "matrix.csv")
            json_path = save_json(payload, base_path / "payload.json")

            self.assertTrue(array_path.exists())
            self.assertTrue(json_path.exists())

            loaded_matrix = pd.read_csv(array_path, index_col=0)
            self.assertEqual(loaded_matrix.shape, (2, 2))

            with json_path.open("r", encoding="utf-8") as handle:
                loaded_payload = json.load(handle)

            self.assertEqual(loaded_payload, {"sp_B": [[0, 1], [1]]})

    def test_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            base_path = Path(tmp_dir)
            target = base_path / "guard.csv"
            save_dataframe(pd.DataFrame({"value": [1]}), target)

            with self.assertRaises(FileExistsError):
                save_dataframe(pd.DataFrame({"value": [2]}), target, overwrite=False)


class TestExperimentPipelinePersistence(unittest.TestCase):

    def _build_config(self, output_dir: Path) -> dict:
        return {
            "name": "Persistence test",
            "seed": 42,
            "save_config": True,
            "data": {
                "data_dir": "unused",
                "expression_data": {
                    "expression_file": "expression.csv",
                    "expression_file_separator": ",",
                    "expression_id_col_idx": 0,
                    "identifier_type": "custom",
                },
                "sample_data": {
                    "sample_data_file": "sample.csv",
                    "sample_data_file_separator": ",",
                    "sample_data_id_col_idx": 0,
                    "sample_data_condition_col_name": "CONDITION",
                    "sample_data_condition_control_value": "reference",
                    "sample_data_condition_case_value": "diseased",
                },
            },
            "output": {
                "output_dir": str(output_dir),
                "overwrite": True,
                "save_plots": False,
                "save_intermediate_results": True,
            },
            "preprocessing": {
                "expression_preprocessing": {
                    "normalize": False,
                    "method_normalize": "DSEQ2",
                    "transform": False,
                    "method_transform": "log2",
                },
                "sample_preprocessing": {
                    "drop_na_columns": [],
                    "drop_duplicates": False,
                    "filters": [],
                },
            },
            "shapley": {
                "discriminant_method_lower_bound": 0,
                "discriminant_method_upper_bound": 100,
            },
            "logging": {
                "level": "INFO",
                "log_to_file": False,
                "log_dir": None,
            },
        }

    def test_pipeline_persists_expected_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir) / "out"
            config = self._build_config(output_dir)

            expression_df = pd.DataFrame(
                {
                    "ref_1": [1.0, 4.0],
                    "ref_2": [2.0, 5.0],
                    "case_1": [10.0, 20.0],
                },
                index=["gene_1", "gene_2"],
            )
            sample_df = pd.DataFrame(
                {"CONDITION": ["reference", "reference", "diseased"]},
                index=["ref_1", "ref_2", "case_1"],
            )

            pipeline = ExperimentPipeline(config)
            pipeline.data_loader.load_data = lambda: (expression_df, sample_df)
            pipeline.sample_preprocessor.preprocess = lambda df: df
            pipeline.expression_preprocessor.preprocess = lambda df, sample_ids_to_keep: df

            shapley_values = pipeline.run()

            self.assertEqual(len(shapley_values), 2)
            self.assertTrue((output_dir / "boolean_expression_matrix.csv").exists())
            self.assertTrue((output_dir / "sp_B.json").exists())
            self.assertTrue((output_dir / "coalitions.json").exists())
            self.assertTrue((output_dir / "shapley_values.csv").exists())
            self.assertTrue((output_dir / "config.json").exists())

            with (output_dir / "config.json").open("r", encoding="utf-8") as handle:
                saved_config = json.load(handle)

            self.assertEqual(saved_config["name"], "Persistence test")

            analyzer = pipeline.shapley_analyzer
            self.assertTrue(hasattr(analyzer, "boolean_expression_matrix"))
            self.assertEqual(analyzer.boolean_expression_matrix.shape, (2, 1))
