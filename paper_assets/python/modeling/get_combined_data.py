from caples_data import CombinedData, CombinedParams
import polars as pl


def get_combined_data() -> CombinedData:
    aru = pl.read_parquet(
        "/Users/mschulist/github/COMB/paper_assets/python/caples_data/src/caples_data/data/outputs_agg_20260103_210827.parquet"
    ).filter(pl.col("point") != 0)

    pc = pl.read_csv(
        "/Users/mschulist/github/COMB/paper_assets/python/caples_data/src/caples_data/data/PC_delinted.csv"
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
        years=[2020],
        pc_species_col="birdCode_fk",
        pc_count_col="abun",
        pc_datetime_col="DateTime",
        pc_point_col="point_ID_fk",
        pc_visit_index_col="visit",
        species=["herwar"],
        aru_threshold=0,
    )
    combined = CombinedData(
        aru,
        pc,
        spatial,
        combined_params,
    )

    # gpd.GeoDataFrame(
    #     combined.spatial.spatial_data.to_pandas(),
    #     geometry=gpd.GeoSeries.from_wkt(
    #         combined.spatial.spatial_data.to_pandas()["geometry"]
    #     ),
    # ).plot(
    #     column="burn",
    #     cmap="RdYlGn_r",  # Red for high burn, green for low burn
    #     legend=True,
    # )
    # plt.show()

    return combined


if __name__ == "__main__":
    combined = get_combined_data()

    print(combined.combined_data.y_index)
