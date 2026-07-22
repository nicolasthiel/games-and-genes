import logging
from pathlib import Path

import pandas as pd
from games_and_genes.configs.schema import ExperimentConfig
from games_and_genes.data_loading import DataLoader
from games_and_genes.preprocessing.expression_preprocessing import ExpressionPreprocessor
from games_and_genes.preprocessing.sample_preprocessing import SamplePreprocessor
from games_and_genes.shapley import ShapleyAnalyzer
from games_and_genes.utils import save_dataframe, save_json, save_ndarray, save_series

logger = logging.getLogger(__name__)


class ExperimentPipeline:

    def __init__(self, config: ExperimentConfig):
        self.config = config
        self.data_loader = DataLoader(config.get("data"))
        self.sample_preprocessor = SamplePreprocessor(config.get("preprocessing").get("sample_preprocessing"))
        self.expression_preprocessor = ExpressionPreprocessor(config.get("preprocessing").get("expression_preprocessing"))
        self.shapley_analyzer = ShapleyAnalyzer(config.get("shapley"))


    def _output_path(self, filename: str) -> Path:
        return Path(self.config["output"]["output_dir"]) / filename


    def _save_intermediate_artifacts(self) -> None:
        if not self.config["output"].get("save_intermediate_results", False):
            return

        artifacts = self.shapley_analyzer.get_artifacts()
        overwrite = self.config["output"].get("overwrite", True)

        save_dataframe(artifacts["expression_df"], self._output_path("expression_preprocessed.csv"), overwrite=overwrite)
        save_ndarray(artifacts["A"], self._output_path("A.csv"), overwrite=overwrite, index=artifacts["feature_names"], columns=artifacts["reference_samples"] + artifacts["diseased_samples"])
        save_ndarray(artifacts["A_SR"], self._output_path("A_SR.csv"), overwrite=overwrite, index=artifacts["feature_names"], columns=artifacts["reference_samples"])
        save_ndarray(artifacts["A_SD"], self._output_path("A_SD.csv"), overwrite=overwrite, index=artifacts["feature_names"], columns=artifacts["diseased_samples"])
        save_ndarray(artifacts["boolean_expression_matrix"], self._output_path("boolean_expression_matrix.csv"), overwrite=overwrite, index=artifacts["feature_names"], columns=artifacts["diseased_samples"])
        save_json(artifacts["sp_B"], self._output_path("sp_B.json"), overwrite=overwrite)
        save_json(artifacts["coalitions"], self._output_path("coalitions.json"), overwrite=overwrite)
        save_series(pd.Series(artifacts["unanimity_coeffs"], index=[str(coalition) for coalition in artifacts["coalitions"]]), self._output_path("unanimity_coeffs.csv"), overwrite=overwrite, column_name="unanimity_coeff")


    def _save_final_artifacts(self, shapley_values) -> None:
        overwrite = self.config["output"].get("overwrite", True)
        save_series(
            pd.Series(shapley_values, index=self.shapley_analyzer.feature_names),
            self._output_path("shapley_values.csv"),
            overwrite=overwrite,
            column_name="shapley_value",
        )


    def _save_config_snapshot(self) -> None:
        if not self.config.get("save_config", False):
            return

        save_json(self.config, self._output_path("config.json"), overwrite=self.config["output"].get("overwrite", True))


    def run(self):
        logger.info("Starting experiment pipeline")
        
        # Data loading
        expression_df, sample_df = self.data_loader.load_data()

        # Data preprocessing
        sample_df = self.sample_preprocessor.preprocess(sample_df)
        sample_ids_to_keep = sample_df.index.tolist()

        expression_df = self.expression_preprocessor.preprocess(expression_df, sample_ids_to_keep)
        self.shapley_analyzer.fit(
            expression_df,
            reference_samples=sample_df[sample_df[self.config["data"]["sample_data"]["sample_data_condition_col_name"]] == self.config["data"]["sample_data"]["sample_data_condition_control_value"]].index.tolist(),
            diseased_samples=sample_df[sample_df[self.config["data"]["sample_data"]["sample_data_condition_col_name"]] == self.config["data"]["sample_data"]["sample_data_condition_case_value"]].index.tolist()
        )
        shapley_values = self.shapley_analyzer.calculate_shapley_values()
        self._save_intermediate_artifacts()
        self._save_final_artifacts(shapley_values)
        self._save_config_snapshot()
        
        # Should we do statistical DEG here?

        # Should we plot here?
        logger.info("Experiment pipeline completed")
        return shapley_values