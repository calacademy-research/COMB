from datetime import datetime as dt
import json
import logging
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import numpy as np

from caples_data.combine_aru_pc import CombinedParams
import arviz as az


@dataclass
class Run:
    id: int
    combined_params: CombinedParams
    model_name: str
    aru_filename: str
    pc_filename: str
    spatial_filename: str
    results: az.InferenceData


SIM_DB_FILENAME = "model_runs.db"
# where the results go
RESULTS_DIR = "results"


def is_valid_sql_identifier(name: str) -> bool:
    """Check if a string is a valid and safe SQL identifier."""

    if not name or not isinstance(name, str):
        return False

    # Regex to verify that the name starts with a letter or underscore, then
    # follows with letters, numbers or underscores.
    return re.match(r"^[a-zA-Z_][a-zA-Z0-9_]*$", name) is not None


def normalize_sql_value(value: Any) -> Any:
    """Normalize a python value to one of the types supported by SQL."""

    if isinstance(value, list) or isinstance(value, tuple):
        return [normalize_sql_value(v) for v in value]
    elif isinstance(value, dt):
        return value.isoformat()
    elif isinstance(value, np.integer):
        return int(value)
    elif isinstance(value, np.floating):
        return float(value)
    return value


def format_sql_insert_values(
    **kwargs: Any,
) -> tuple[str, str, list[Any]]:
    """Build columns string, placeholders string and values list for SQL INSERT.

    Args:
      **kwargs: Key-value pairs to pass to the SQL statement.

    Returns:
      A tuple of: a formatted columns string, a formatted placeholders string, and
      a list of corresponding values. Safe to be used in SQL INSERT statements.
    """

    for key in kwargs:
        if not is_valid_sql_identifier(key):
            raise ValueError(f"`{key}` is not a valid SQL identifier.")

    columns = list(kwargs.keys())
    placeholders = ["?"] * len(columns)
    values = normalize_sql_value(list(kwargs.values()))

    return f"({', '.join(columns)})", f"({', '.join(placeholders)})", values


def create_sql_json_filter(
    column_name: str, filter_dict: dict[str, Any]
) -> tuple[str, list[Any]]:
    """Create SQLite WHERE clause for filtering JSON column values.

    Args:
        column_name: Name of the JSON column to filter (e.g., "simulation_params")
        filter_dict: Dictionary of key-value pairs to match in the JSON

    Returns:
        A tuple of (where_clause, params) where:
        - where_clause: SQL WHERE conditions joined with AND
        - params: List of parameter values for the prepared statement

    Supported value types:
        - Scalars (int, float, str, bool, None): Direct equality comparison
        - Lists/tuples: JSON array comparison (order-sensitive)

    Examples:
        # Scalar values
        where_clause, params = create_sql_json_filter(
            "simulation_params",
            {"nsites": 50, "beta0": -0.5}
        )

        # List/tuple values (compares as JSON arrays)
        where_clause, params = create_sql_json_filter(
            "simulation_params",
            {"mu": [-1.0, 2.0], "sigma": [1.0, 1.0]}
        )
    """
    if not is_valid_sql_identifier(column_name):
        raise ValueError(f"`{column_name}` is not a valid SQL identifier.")

    if not filter_dict:
        return "1=1", []

    conditions = []
    params = []

    for key, value in filter_dict.items():
        json_path = f"$.{key}"

        # Check if value is a list/tuple BEFORE normalizing
        if isinstance(value, (list, tuple)):
            # Normalize the list elements, then convert to JSON string
            normalized_value = normalize_sql_value(value)
            json_str = json.dumps(normalized_value)
            # Extract the value, then compare JSON representations
            condition = f"json(json_extract({column_name}, '{json_path}')) = json(?)"
            params.append(json_str)
        else:
            # For scalars, normalize and use direct comparison
            normalized_value = normalize_sql_value(value)
            condition = f"json_extract({column_name}, '{json_path}') = ?"
            params.append(normalized_value)

        conditions.append(condition)

    where_clause = " AND ".join(conditions)
    return where_clause, params


@dataclass
class ResultsDB:
    db_path: Path
    db: sqlite3.Connection
    _cursor: sqlite3.Cursor | None = None

    @staticmethod
    def _setup_tables(cursor: sqlite3.Cursor) -> None:
        # skip setting up tables if they exist
        cursor.execute("""
            SELECT name
            FROM sqlite_master
            WHERE name = "runs" AND type = "table"
        """)
        if cursor.fetchone() is not None:
            return

        cursor.execute("PRAGMA foreign_keys = ON")

        # create runs table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS runs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                combined_params TEXT NOT NULL,
                model_name TEXT NOT NULL,
                aru_filename TEXT NOT NULL,
                pc_filename TEXT NOT NULL,
                spatial_filename TEXT NOT NULL
            )
        """)

    @classmethod
    def create(cls, db_path: str | Path):
        db_path = Path(db_path)
        db_path.mkdir(parents=True, exist_ok=True)
        sqlite_path = db_path / SIM_DB_FILENAME
        db = sqlite3.connect(sqlite_path.as_posix())
        db.set_trace_callback(
            lambda statement: logging.info("Executed SQL statement: %s", statement)
        )
        cursor = db.cursor()
        cursor.execute("PRAGMA journal_mode = WAL")
        cls._setup_tables(cursor)
        db.commit()

        (db_path / RESULTS_DIR).mkdir(parents=True, exist_ok=True)

        sim_db = cls(db_path, db)
        return sim_db

    def commit(self):
        self.db.commit()
        if self._cursor is not None:
            self._cursor.close()
            self._cursor = None

    def _get_cursor(self):
        if self._cursor is None:
            self._cursor = self.db.cursor()
        return self._cursor

    def insert_run(
        self,
        combined_params: CombinedParams,
        model_name: str,
        aru_filename: str,
        pc_filename: str,
        spatial_filename: str,
        results: az.InferenceData,
    ):
        cursor = self._get_cursor()
        columns_str, placeholders_str, values = format_sql_insert_values(
            combined_params=combined_params.to_str(),
            model_name=model_name,
            aru_filename=aru_filename,
            pc_filename=pc_filename,
            spatial_filename=spatial_filename,
        )

        cursor.execute(
            f"""
            INSERT INTO runs {columns_str}
            VALUES {placeholders_str}
            """,
            values,
        )

        run_id = cursor.lastrowid
        if run_id is None:
            raise RuntimeError("Error inserting run into the database")

        results.to_netcdf(str(self.db_path / RESULTS_DIR / f"{run_id}.nc"))
        return run_id

    def get_all_run_ids(self, filter_params: dict[str, Any] = {}):
        cursor = self._get_cursor()
        where_clause, params = create_sql_json_filter("combined_params", filter_params)

        query = f"""
            SELECT id
            FROM runs
            WHERE {where_clause}
        """

        cursor.execute(query, params)
        results = cursor.fetchall()
        return [row[0] for row in results]

    def get_run(self, run_id: int):
        cursor = self._get_cursor()
        cursor.execute(
            """
            SELECT *
            FROM runs
            WHERE id = ?
            """,
            (run_id,),
        )

        result = cursor.fetchone()
        if result is None:
            raise RuntimeError(f"run id not found: {run_id}")

        columns = [col[0] for col in cursor.description]
        row_dict = dict(zip(columns, result))

        # Load results from disk
        results = az.from_netcdf(self.db_path / RESULTS_DIR / f"{run_id}.nc")

        run = Run(
            id=row_dict["id"],
            combined_params=CombinedParams.from_str(str(row_dict["combined_params"])),
            model_name=str(row_dict["model_name"]),
            aru_filename=str(row_dict["aru_filename"]),
            pc_filename=str(row_dict["pc_filename"]),
            spatial_filename=str(row_dict["spatial_filename"]),
            results=results,
        )
        return run
