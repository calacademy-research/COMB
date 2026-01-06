from get_combined_data import get_combined_data
from models.model_zoo import get_model_by_name
from models.model_iterface import SimulationParams
import numpy as np


if __name__ == "__main__":
    combined = get_combined_data()
    model = get_model_by_name("multi_year_multi_species_all")
    trace = model.run_model(combined.combined_data)

    trace.to_netcdf("data/trace.nc")
