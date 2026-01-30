from datetime import datetime as dt
import logging
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import numpy as np
import json

from caples_data.combine_aru_pc import COMBData
from models.model_iterface import SimulationParams
import pickle
import arviz as az
from .study_params import StudyParams


@dataclass
class Study:
    id: int
    name: str
    datetime: dt
    study_params: StudyParams


@dataclass
class SimParams:
    id: int
    study_id: int
    simulation_params: SimulationParams


@dataclass
class Dataset:
    id: int
    sim_param_id: int
    data: COMBData


@dataclass
class Run:
    id: int
    dataset_id: int
    model_name: str
    results: az.InferenceData


SIM_DB_FILENAME = "simulations.db"
# where the simulation results go
RESULTS_DIR = "results"
# where the simulated data goes
DATA_DIR = "data"


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
class SimulationsDB:
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

        # create studies table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS studies (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                datetime TEXT NOT NULL,
                study_params TEXT NOT NULL
            )               
        """)

        # create sim_params table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS sim_params (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                study_id INTEGER NOT NULL,
                simulation_params TEXT NOT NULL,
                FOREIGN KEY (study_id) REFERENCES studies(id)
                    ON DELETE CASCADE
            )
        """)

        # create the datasets table
        # the data will be on disk and named according to the dataset_id
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS datasets (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sim_param_id INTEGER NOT NULL,
                FOREIGN KEY (sim_param_id) REFERENCES sim_params(id)
                    ON DELETE CASCADE
            )
        """)

        # create the runs table
        # the results will be on disk and named according to the run_id
        cursor.execute("""
           CREATE TABLE IF NOT EXISTS runs (
               id INTEGER PRIMARY KEY AUTOINCREMENT,
               dataset_id INTEGER NOT NULL,
               model_name TEXT NOT NULL,
               FOREIGN KEY (dataset_id) REFERENCES datasets(id)
                   ON DELETE CASCADE
           )            
        """)

        # Create indices for optimized queries
        # Index on foreign key for JOIN operations
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_sim_params_study_id 
            ON sim_params(study_id)
        """)

        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_datasets_sim_param_id 
            ON datasets(sim_param_id)
        """)

        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_runs_dataset_id 
            ON runs(dataset_id)
        """)

        # Index on model_name for filtering runs by model
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_runs_model_name 
            ON runs(model_name)
        """)

        # Composite index for common filter patterns (dataset_id + model_name)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_runs_dataset_model 
            ON runs(dataset_id, model_name)
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

        (db_path / DATA_DIR).mkdir(parents=True, exist_ok=True)
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

    def insert_study(
        self,
        name: str,
        datetime: dt,
        study_params: StudyParams,
    ):
        cursor = self._get_cursor()
        colums_str, placeholders_str, values = format_sql_insert_values(
            name=name, datetime=datetime, study_params=study_params.to_str()
        )
        cursor.execute(
            f"""
            INSERT INTO studies {colums_str}
            VALUES {placeholders_str}
            """,
            values,
        )
        study_id = cursor.lastrowid
        if study_id is None:
            raise RuntimeError("Error inserting study into the database")
        return study_id

    def get_study(self, study_id: int):
        cursor = self._get_cursor()
        cursor.execute(
            """
            SELECT *
            FROM studies
            WHERE id = ?
            """,
            (study_id,),
        )
        result = cursor.fetchone()
        if result is None:
            raise KeyError(f"study id not found: {study_id}")

        columns = [col[0] for col in cursor.description]
        study = Study(**dict(zip(columns, result)))
        if study.datetime is not None:
            study.datetime = dt.fromisoformat(str(study.datetime))
        study.study_params = StudyParams.from_str(str(study.study_params))
        return study

    def insert_sim_params(
        self,
        study_id: int,
        simulation_params: SimulationParams,
    ):
        cursor = self._get_cursor()
        columns_str, placeholders_str, values = format_sql_insert_values(
            study_id=study_id,
            simulation_params=simulation_params.to_str(),
        )
        cursor.execute(
            f"""
            INSERT INTO sim_params {columns_str}
            VALUES {placeholders_str}
            """,
            values,
        )

        sim_param_id = cursor.lastrowid
        if sim_param_id is None:
            raise RuntimeError("Error inserting sim_params into the database")
        return sim_param_id

    def get_sim_params(self, sim_param_id: int):
        cursor = self._get_cursor()
        cursor.execute(
            """
            SELECT *
            FROM sim_params
            WHERE id = ?
            """,
            (sim_param_id,),
        )
        result = cursor.fetchone()
        if result is None:
            raise KeyError(f"sim_param id not found: {sim_param_id}")

        columns = [col[0] for col in cursor.description]
        sim_params = SimParams(**dict(zip(columns, result)))
        sim_params.simulation_params = SimulationParams.from_str(
            str(sim_params.simulation_params)
        )
        return sim_params

    def get_all_sim_param_ids(self, study_id: int) -> list[int]:
        cursor = self._get_cursor()
        cursor.execute(
            """
            SELECT id
            FROM sim_params
            WHERE study_id = ?
            """,
            (study_id,),
        )
        results = cursor.fetchall()
        return [row[0] for row in results]

    def get_sim_param_ids_by_filter(
        self, study_id: int, filter_params: dict[str, Any]
    ) -> list[int]:
        """Get all sim_param IDs that match the given filter parameters."""
        cursor = self._get_cursor()
        where_clause, params = create_sql_json_filter(
            "simulation_params", filter_params
        )

        query = f"""
            SELECT id
            FROM sim_params
            WHERE study_id = ? AND {where_clause}
        """

        cursor.execute(query, [study_id] + params)
        results = cursor.fetchall()
        return [row[0] for row in results]

    def insert_dataset(
        self,
        sim_param_id: int,
        data: COMBData,
    ):
        cursor = self._get_cursor()
        columns_str, placeholders_str, values = format_sql_insert_values(
            sim_param_id=sim_param_id,
        )
        cursor.execute(
            f"""
            INSERT INTO datasets {columns_str}
            VALUES {placeholders_str}
            """,
            values,
        )

        dataset_id = cursor.lastrowid
        if dataset_id is None:
            raise RuntimeError("Error inserting dataset into the database")

        with open(self.db_path / DATA_DIR / f"{dataset_id}.pickle", "wb") as f:
            pickle.dump(data, f)
        return dataset_id

    def get_dataset(self, dataset_id: int):
        cursor = self._get_cursor()
        cursor.execute(
            """
            SELECT *
            FROM datasets
            WHERE id = ?
            """,
            (dataset_id,),
        )
        result = cursor.fetchone()
        if result is None:
            raise KeyError(f"dataset id not found: {dataset_id}")

        columns = [col[0] for col in cursor.description]
        row_dict = dict(zip(columns, result))

        # Load data from disk
        with open(self.db_path / DATA_DIR / f"{dataset_id}.pickle", "rb") as f:
            data = pickle.load(f)

        dataset = Dataset(
            id=row_dict["id"], sim_param_id=row_dict["sim_param_id"], data=data
        )
        return dataset

    def get_all_dataset_ids(self, study_id: int) -> list[int]:
        """Get all dataset IDs from the database, optionally filtered by study."""
        cursor = self._get_cursor()
        cursor.execute(
            """
            SELECT datasets.id
            FROM datasets
            JOIN sim_params ON datasets.sim_param_id = sim_params.id
            WHERE sim_params.study_id = ?
            """,
            (study_id,),
        )
        results = cursor.fetchall()
        return [row[0] for row in results]

    def insert_run(
        self,
        dataset_id: int,
        model_name: str,
        results: az.InferenceData,
    ):
        cursor = self._get_cursor()
        columns_str, placeholders_str, values = format_sql_insert_values(
            dataset_id=dataset_id,
            model_name=model_name,
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
            dataset_id=row_dict["dataset_id"],
            model_name=row_dict["model_name"],
            results=results,
        )
        return run

    def get_run_ids_by_filters(
        self,
        study_id: int,
        params_filter: dict[str, Any] | None = None,
        runs_filter: dict[str, Any] | None = None,
    ) -> list[int]:
        """Get all run IDs that match the combination of params and runs filters.

        Args:
            study_id: The study ID to filter by
            params_filter: Dictionary to filter by simulation_params in sim_params table
                          (e.g., {"nsites": 50, "beta0": -0.5})
            runs_filter: Dictionary to filter by columns in runs table
                        (e.g., {"model_name": "model_xyz"})

        Returns:
            List of run IDs matching both filters

        Examples:
            # Filter by both params and model name
            run_ids = db.get_run_ids_by_filters(
                study_id=1,
                params_filter={"nsites": 50},
                runs_filter={"model_name": "single_year_single_species"}
            )

            # Filter only by params
            run_ids = db.get_run_ids_by_filters(
                study_id=1,
                params_filter={"beta0": -0.5, "nsites": 50}
            )

            # Filter only by runs
            run_ids = db.get_run_ids_by_filters(
                study_id=1,
                runs_filter={"model_name": "multi_year_multi_species"}
            )
        """
        cursor = self._get_cursor()

        # Build the params filter WHERE clause
        params_where = "1=1"
        params_values = []
        if params_filter:
            params_where, params_values = create_sql_json_filter(
                "simulation_params", params_filter
            )

        # Build the runs filter WHERE clause
        runs_where = "1=1"
        runs_values = []
        if runs_filter:
            conditions = []
            for key, value in runs_filter.items():
                if not is_valid_sql_identifier(key):
                    raise ValueError(f"`{key}` is not a valid SQL identifier.")
                normalized_value = normalize_sql_value(value)
                conditions.append(f"runs.{key} = ?")
                runs_values.append(normalized_value)
            runs_where = " AND ".join(conditions) if conditions else "1=1"

        # Combine query with JOINs
        query = f"""
            SELECT runs.id
            FROM runs
            JOIN datasets ON runs.dataset_id = datasets.id
            JOIN sim_params ON datasets.sim_param_id = sim_params.id
            WHERE sim_params.study_id = ?
                AND {params_where}
                AND {runs_where}
        """

        all_params = [study_id] + params_values + runs_values
        cursor.execute(query, all_params)
        results = cursor.fetchall()
        return [row[0] for row in results]
