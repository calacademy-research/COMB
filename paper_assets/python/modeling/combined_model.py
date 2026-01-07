from get_combined_data import get_combined_data
from models.model_zoo import get_model_by_name, ModelNames
from model_runs.model_runs_db import ResultsDB
from caples_data.combine_aru_pc import CombinedParams
import numpy as np

"""
Main script to run the combined model on the caples data
"""

aru_filename = "data/outputs_agg_20260103_210827.parquet"
pc_filename = "data/PC_delinted_2018-2023.csv"
spatial_filename = "../spatial/data/burn_data_by_point.csv"

# model_name: ModelNames = "single_year_single_species_all_datetime"


models_to_run: list[ModelNames] = [
    "single_year_single_species_all",
    "single_year_single_species_no_pc",
    "single_year_single_species_no_aru",
    "single_year_single_species_no_scores",
    "single_year_single_species_no_aru_no_pc",
    "single_year_single_species_no_scores_no_pc",
    "single_year_single_species_no_scores_no_aru",
]

species = [
    "whhwoo",
    "rebsap",
    "pilwoo",
    "norfli",
    "haiwoo",
    "bkbwoo",
]

db_path = "data/results_db"

if __name__ == "__main__":
    results_db = ResultsDB.create(db_path)

    for s in species:
        combined = get_combined_data(
            aru_filename=aru_filename,
            pc_filename=pc_filename,
            spatial_filename=spatial_filename,
            species=[s],
        )
        combined.combined_data.y_index = np.nan_to_num(combined.combined_data.y_index)
        for model_name in models_to_run:
            if results_db.get_all_run_ids(dict(species=[s]), model_name):
                print(
                    f"skipping species: {s}, model_name: {model_name} because it already exists"
                )
                continue
            model = get_model_by_name(model_name)
            trace = model.run_model(combined.combined_data)

            results_db.insert_run(
                combined.params,
                model_name,
                aru_filename,
                pc_filename,
                spatial_filename,
                trace,
            )
            results_db.commit()

    # trace.to_netcdf("data/trace.nc")
