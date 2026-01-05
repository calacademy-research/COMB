from pathlib import Path
import polars as pl
from datetime import datetime


class ClassifierOutputs:
    def __init__(self, parquet_file: str | Path, aru_to_point_file: str | Path):
        self.table = pl.scan_parquet(parquet_file)
        self.aru_to_point_file = aru_to_point_file

    def write_output_parquet(self, output_file: Path, agg: bool = True):
        """
        Writes the outputs for the given table
        """
        lazy = self.table
        print(lazy)
        # we have these cols: ['filename', 'logit', 'timestamp_s', 'window_id', 'label']

        # agg means that we take the maximum logit per label for each recording
        if agg:
            lazy = lazy.group_by(["filename", "label"]).agg(pl.col("logit").max())

        # get the datetime from the filename: (example) 2024/9-BROWN-CAPL_20240609_000000.wav
        lazy = lazy.with_columns(
            datetime=pl.col("filename")
            .str.extract(r"(\d{8}_\d{6})\.wav$", 1)
            .str.to_datetime("%Y%m%d_%H%M%S")
        )

        # strip the year prefix from filename for joining (e.g., "2020/" -> "")
        lazy = lazy.with_columns(
            filename=pl.col("filename").str.replace(r"^\d{4}/", "")
        )

        # map the filename to the point
        aru2point = pl.read_csv(self.aru_to_point_file)
        lazy = lazy.join(aru2point.lazy(), on="filename", how="left")

        lazy.sink_parquet(output_file)


if __name__ == "__main__":
    classifier_outputs = ClassifierOutputs(
        Path(
            "/home/mschulist/perch-agile-modeling/python_server/data/classify/10/predictions.parquet"
        ),
        Path("/home/mschulist/caples_sound/aru2point/aru2point_all.csv"),
    )

    date = datetime.now().strftime("%Y%m%d_%H%M%S")

    classifier_outputs.write_output_parquet(
        Path(
            f"/home/mschulist/COMB/paper_assets/python/caples_data/data/outputs_agg_{date}.parquet"
        )
    )
