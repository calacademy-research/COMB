import polars as pl
from .combine_aru_pc import CombinedData, CombinedParams


aru = pl.read_parquet(
    "/Users/mschulist/github/COMB/paper_assets/python/caples_data/data/outputs_agg_20251224_180943.parquet"
).filter(pl.col("point") != 0)

pc = pl.read_csv(
    "/Users/mschulist/github/COMB/paper_assets/python/caples_data/data/PC_delinted.csv"
).with_columns(
    visit=pl.col("visit") - 1,
    DateTime=pl.col("DateTime").str.to_datetime("%Y-%m-%dT%H:%M:%SZ"),
)

spatial = pl.read_csv(
    "/Users/mschulist/github/COMB/paper_assets/python/spatial/data/burn_data_by_point.csv"
)

combined_params = CombinedParams(
    aru_species_col="label",
    aru_visit_limit=24,
    years=[2021],
    pc_species_col="birdCode_fk",
    pc_count_col="abun",
    pc_datetime_col="DateTime",
    pc_point_col="point_ID_fk",
    pc_visit_index_col="visit",
    species=["gockin"],
)
combined = CombinedData(
    aru,
    pc,
    spatial,
    combined_params,
)

print(combined.combined_data)
