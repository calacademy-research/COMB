from get_combined_data import get_combined_data
from models.single_year_single_species_all import SingleYearSingleSpeciesAll


if __name__ == "__main__":
    combined = get_combined_data()
    trace = SingleYearSingleSpeciesAll.run_model(combined.combined_data)

    trace.to_netcdf("data/trace.nc")
