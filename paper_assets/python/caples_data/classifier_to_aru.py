from pyiceberg.catalog.sql import SqlCatalog
from pathlib import Path
from urllib.parse import urlparse
import os
import polars as pl


class ClassifierOutputs:
    def __init__(
        self, classifier_warehouse_path: Path, aru_to_point_csv: Path, project_id: int
    ):
        self.catalog = SqlCatalog(
            "default",
            **{
                "uri": f"sqlite:///{classifier_warehouse_path}/pyiceberg_catalog.db",
                "warehouse": f"file://{classifier_warehouse_path}",
            },
        )
        self.project_id = project_id
        self.aru_to_point_csv = aru_to_point_csv

    def write_output_parquet(
        self, output_file: Path, table_name: str, agg: bool = True
    ):
        """
        Writes the outputs for the given table
        """
        table = self.catalog.load_table(f"{self.project_id}.{table_name}")
        parsed = urlparse(table.location())
        relative_path = parsed.netloc + parsed.path
        data_path = Path(os.getcwd()) / relative_path / "data/*.parquet"

        lazy = pl.scan_parquet(data_path)
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

        # map the filename to the point
        aru2point = pl.read_csv(self.aru_to_point_csv)
        lazy = lazy.join(aru2point.lazy(), on="filename", how="left")

        lazy.sink_parquet(output_file)

    def get_all_tables(self):
        return self.catalog.list_tables(str(self.project_id))


if __name__ == "__main__":
    # NOTE: Run this script from /home/mschulist/perch-agile-modeling/python_server
    # because the catalog uses relative file:// paths

    os.chdir("/home/mschulist/perch-agile-modeling/python_server")
    classifier_outputs = ClassifierOutputs(
        Path("/home/mschulist/perch-agile-modeling/python_server/data/warehouse"),
        Path("/home/mschulist/caples_sound/aru2point/aru2point_all.csv"),
        1,
    )

    table_name = "20251222_091414"

    classifier_outputs.write_output_parquet(
        Path(
            f"/home/mschulist/COMB/paper_assets/python/caples_data/data/outputs_{table_name}.parquet"
        ),
        table_name,
    )
