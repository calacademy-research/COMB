import polars as pl
from dataclasses import dataclass
from datetime import time
from typing import Tuple, List, Union
import numpy as np


@dataclass
class ARUDataParams:
    """
    Data class to hold the ARU data parameters.
    """

    threshold: float = 0.5

    species_col: str = "species"
    logit_col: str = "logit"
    datetime_col: str = "datetime"
    point_col: str = "point"

    skip_first_day: bool = True
    aru_visit_limit: int = 24
    time_range: Tuple[time, time] = (time(6), time(9, 59))

    species: Union[List[str], None] = None

    years: Union[List[int], None] = None

    # the following are needed because we need the following indices to be
    # consistent between the ARU and PC data.
    # either all of these should be None or all of them should be provided
    point_index: Union[dict[int, int], None] = None
    species_index: Union[dict[str, int], None] = None
    year_index: Union[dict[int, int], None] = None


class AruData:
    def __init__(self, aru_data: pl.DataFrame, aru_data_params: ARUDataParams):
        """
        Initializes an AruData object given the ARU data.

        Args:
            aru_data (pl.DataFrame): The ARU data.
            aru_data_params (ARUDataParams): The ARU data parameters.
        """

        self.aru_data = aru_data
        self.aru_data_params = aru_data_params
        self._verify_aru_data()

        self.aru_data_dict = self.convert_df_to_dict()

    def convert_df_to_dict(self):
        """
        Converts the df of ARU data into a dictionary.

        This dictionary is used by JAGS to run the occupancy model.

        The dict will have two keys:
        - y_aru: The ARU data
        - scores

        y_aru will be a np array with the following shape: (species, year, site, visit)

        Returns:
            dict: The ARU data as a dictionary.
        """
        aru_data = self._assign_indices(self.aru_data)
        aru_data_dict = self._create_aru_arrays(aru_data)
        return aru_data_dict

    def _create_aru_arrays(self, aru_data: pl.DataFrame) -> dict[str, Union[np.ndarray, int]]:
        """
        Creates the y_aru and scores arrays from the ARU data.

        Also returns the number of species, years, sites, and visits.

        The arrays will have the following shape: (species, year, site, visit)

        Returns:
            Dict[str, Union[np.ndarray, int]]: The ARU data as a dictionary.
        """
        if self.aru_data_params.species_index is not None:
            n_species = max(self.aru_data_params.species_index.values()) + 1
        else:
            n_species = aru_data["species_index"].max()
            if isinstance(n_species, int):
                n_species = n_species + 1
            else:
                raise ValueError("No species found in ARU data.")

        if self.aru_data_params.year_index is not None:
            n_years = max(self.aru_data_params.year_index.values()) + 1
        else:
            n_years = aru_data["year_index"].max()
            if isinstance(n_years, int):
                n_years = n_years + 1
            else:
                raise ValueError("No years found in ARU data or n_years is not an int")

        if self.aru_data_params.point_index is not None:
            n_sites = max(self.aru_data_params.point_index.values()) + 1
        else:
            n_sites = aru_data["point_index"].max()
            if isinstance(n_sites, int):
                n_sites = n_sites + 1
            else:
                raise ValueError("No sites found in ARU data or n_sites is not an int")

        n_surveys_aru = aru_data["visit_index"].max()
        if isinstance(n_surveys_aru, int):
            n_surveys_aru = n_surveys_aru + 1
        else:
            raise ValueError("No visits found in ARU data or n_visits is not an int")

        shape = (n_species, n_years, n_sites, n_surveys_aru)

        # y_aru is the binary ARU data
        y_aru = np.empty(shape, dtype=np.float32)
        y_aru.fill(np.nan)

        # scores is the continuous ARU data
        scores = np.empty(shape, dtype=np.float32)
        scores.fill(np.nan)

        # date_aru is the julian date of the survey
        date_aru = np.empty(shape, dtype=np.float32)
        date_aru.fill(np.nan)

        # time_aru is the time of day of the survey
        time_aru = np.empty(shape, dtype=np.float32)
        time_aru.fill(np.nan)

        for row in aru_data.iter_rows(named=True):
            species_index = row["species_index"]
            year_index = row["year_index"]
            site_index = row["point_index"]
            visit_index = row["visit_index"]
            logit = row[self.aru_data_params.logit_col]
            datetime = row[self.aru_data_params.datetime_col]

            if logit >= self.aru_data_params.threshold:
                y_aru[species_index, year_index, site_index, visit_index] = 1
            else:
                y_aru[species_index, year_index, site_index, visit_index] = 0

            scores[species_index, year_index, site_index, visit_index] = logit

            date_aru[species_index, year_index, site_index, visit_index] = (
                datetime.timetuple().tm_yday
            )

            time_aru[species_index, year_index, site_index, visit_index] = (
                datetime.hour * 60 + datetime.minute
            )

        return {
            "y_aru": y_aru,
            "scores": scores,
            "n_species": n_species,
            "n_years": n_years,
            "n_sites": n_sites,
            "n_surveys_aru": n_surveys_aru,
            "date_aru": date_aru,
            "time_aru": time_aru,
        }

    def _verify_aru_data(self) -> None:
        """
        Verify that the ARU data has the correct columns and data types.
        """
        datetime_col = self.aru_data_params.datetime_col
        species_col = self.aru_data_params.species_col
        logit_col = self.aru_data_params.logit_col
        point_col = self.aru_data_params.point_col

        # check the datetime column
        if datetime_col in self.aru_data.columns:
            if self.aru_data[datetime_col].dtype != pl.Datetime:
                self.aru_data = self.aru_data.with_columns(
                    pl.col(datetime_col).str.to_datetime().alias(datetime_col)
                )

        # check the species column
        if species_col not in self.aru_data.columns:
            raise ValueError("ARU data must have a 'species' column.")

        # check the logit column
        if logit_col not in self.aru_data.columns:
            raise ValueError("ARU data must have a 'logit' column.")
        else:
            if self.aru_data[logit_col].dtype != "float":
                self.aru_data = self.aru_data.with_columns(pl.col(logit_col).cast(pl.Float32))

        # check the point column
        if point_col not in self.aru_data.columns:
            raise ValueError("ARU data must have a 'point' column.")
        else:
            if self.aru_data[point_col].dtype != pl.Int32:
                self.aru_data = self.aru_data.with_columns(pl.col(point_col).cast(pl.Int32))

    def _assign_indices(self, aru_data: pl.DataFrame) -> pl.DataFrame:
        """
        Assigns indices to the ARU data.

        Each species, point, year, and visit will have a unique index.

        Returns:
            pl.DataFrame: The ARU data with an index assigned to each point.
        """
        point_col = self.aru_data_params.point_col
        datetime_col = self.aru_data_params.datetime_col
        species_col = self.aru_data_params.species_col
        time_range = self.aru_data_params.time_range
        aru_visit_limit = self.aru_data_params.aru_visit_limit

        # filter out null values for point and datetime
        aru_data = aru_data.filter(
            (pl.col(point_col).is_not_null()) & (pl.col(datetime_col).is_not_null())
        )
        if self.aru_data_params.years:
            aru_data = aru_data.filter(
                pl.col("datetime").dt.year().is_in(self.aru_data_params.years)
            )

        if self.aru_data_params.species:
            aru_data = aru_data.filter(pl.col(species_col).is_in(self.aru_data_params.species))

        # filter to the time range
        aru_data = aru_data.filter(
            (pl.col("datetime").dt.time() >= pl.lit(time_range[0]))
            & (pl.col("datetime").dt.time() <= pl.lit(time_range[1]))
        )

        # sort by species and datetime for consistency when assigning indices
        aru_data = aru_data.sort([species_col, datetime_col, point_col])

        aru_data = aru_data.with_columns(
            year=pl.col(datetime_col).dt.year(),
        )

        # assign point, year, and species indices if not provided
        if (
            self.aru_data_params.point_index is not None
            and self.aru_data_params.species_index is not None
            and self.aru_data_params.year_index is not None
        ):
            aru_data = aru_data.with_columns(
                point_index=pl.col(point_col).replace_strict(self.aru_data_params.point_index),
                species_index=pl.col(species_col).replace_strict(
                    self.aru_data_params.species_index
                ),
                year_index=pl.col("year").replace_strict(self.aru_data_params.year_index),
            )
        elif (
            self.aru_data_params.point_index is not None
            or self.aru_data_params.species_index is not None
            or self.aru_data_params.year_index is not None
        ):
            raise ValueError(
                "If one of point_index, species_index, or year_index is provided, all three must be provided."
            )
        else:
            aru_data = aru_data.with_columns(
                point_index=pl.col(point_col).rank("dense").cast(pl.Int32) - 1,
                species_index=pl.col(species_col).rank("dense").cast(pl.Int32) - 1,
                year_index=pl.col("year").rank("dense").cast(pl.Int32) - 1,
            )

        # skip the first day if needed
        if self.aru_data_params.skip_first_day:
            first_dates = aru_data.group_by(["point", "year"]).agg(
                pl.col("datetime").min().alias("first_date")
            )

            # Join the first dates back to the original DataFrame
            aru_data = aru_data.join(first_dates, on=[point_col, "year"], how="left")

            # Filter out rows where 'datetime' matches 'first_date'
            aru_data = aru_data.filter(
                pl.col("datetime").dt.date() != pl.col("first_date").dt.date()
            )

        # Rest of the code for assigning indices, filtering, and sorting
        aru_data = aru_data.with_columns(
            visit_index=pl.col("datetime").cum_count().over(["point", "year", "species_index"]) - 1
        )

        # Filter to the visit limit
        aru_data = aru_data.filter(pl.col("visit_index") < aru_visit_limit)

        return aru_data


@dataclass
class PCDataParams:
    """
    Data class to hold the Point Count data parameters.
    """

    species_col: str = "species"
    count_col: str = "count"
    datetime_col: str = "datetime"
    point_col: str = "point"
    visit_index_col: str = "visit_index"

    species: Union[List[str], None] = None

    years: Union[List[int], None] = None

    pc_visit_limit: int = 3

    # the following are needed because we need the following indices to be
    # consistent between the ARU and PC data.
    # either all of these should be None or all of them should be provided
    point_index: Union[dict[int, int], None] = None
    species_index: Union[dict[str, int], None] = None
    year_index: Union[dict[int, int], None] = None


class PcData:
    def __init__(self, pc_data: pl.DataFrame, pc_data_params: PCDataParams):
        """
        Initializes a PcData object given the Point Count data.

        Args:
            pc_data (pl.DataFrame): The Point Count data.
            pc_data_params (PCDataParams): The Point Count data parameters.
        """

        self.pc_data = pc_data
        self.pc_data_params = pc_data_params
        self._verify_pc_data()

        self.pc_data_dict = self.convert_df_to_dict()

    def _verify_pc_data(self):
        """
        Verify that the Point Count data has the correct columns and data types.
        """
        species_col = self.pc_data_params.species_col
        count_col = self.pc_data_params.count_col
        point_col = self.pc_data_params.point_col
        datetime_col = self.pc_data_params.datetime_col

        # check the species column
        if species_col not in self.pc_data.columns:
            raise ValueError(f"Point Count data must have a {species_col} column.")

        # check the count column
        if count_col not in self.pc_data.columns:
            raise ValueError(f"Point Count data must have a {count_col} column.")
        else:
            if self.pc_data[count_col].dtype != "int":
                self.pc_data = self.pc_data.with_columns(self.pc_data[count_col].cast(pl.Int32))

        # check the point column
        if point_col not in self.pc_data.columns:
            raise ValueError(f"Point Count data must have a {point_col} column.")
        else:
            if self.pc_data[point_col].dtype != pl.Int32:
                self.pc_data = self.pc_data.with_columns(self.pc_data[point_col].cast(pl.Int32))

        # check the datetime column
        if datetime_col in self.pc_data.columns:
            if self.pc_data[datetime_col].dtype != pl.Datetime:
                self.pc_data = self.pc_data.with_columns(
                    pl.col(datetime_col).str.to_datetime().alias(datetime_col)
                )
        else:
            raise ValueError(f"Point Count data must have a {datetime_col} column.")

    def convert_df_to_dict(self):
        """
        Converts the df of Point Count data into a dictionary.

        This dictionary is used by JAGS to run the occupancy model.

        The dict will have two keys:
        - y_pc: The Point Count data
        - n_visits: The number of visits

        Returns:
            dict: The Point Count data as a dictionary.
        """
        pc_data = self._assign_indices(self.pc_data)
        pc_data_dict = self._create_pc_arrays(pc_data)
        return pc_data_dict

    def _create_pc_arrays(self, pc_data: pl.DataFrame) -> dict[str, Union[np.ndarray, int]]:
        """
        Creates the y_pc array from the Point Count data.

        Also returns the number of visits.

        The array will have the following shape: (species, year, site, visit)

        Returns:
            Dict[str, Union[np.ndarray, int]]: The Point Count data as a dictionary.
        """
        if self.pc_data_params.species_index is not None:
            n_species = max(self.pc_data_params.species_index.values()) + 1
        else:
            n_species = pc_data["species_index"].max()
            if isinstance(n_species, int):
                n_species = n_species + 1
            else:
                raise ValueError("No species found in PC data.")

        if self.pc_data_params.year_index is not None:
            n_years = max(self.pc_data_params.year_index.values()) + 1
        else:
            n_years = pc_data["year_index"].max()
            if isinstance(n_years, int):
                n_years = n_years + 1
            else:
                raise ValueError("No years found in PC data or n_years is not an int")

        if self.pc_data_params.point_index is not None:
            n_sites = max(self.pc_data_params.point_index.values()) + 1
        else:
            n_sites = pc_data["point_index"].max()
            if isinstance(n_sites, int):
                n_sites = n_sites + 1
            else:
                raise ValueError("No sites found in PC data or n_sites is not an int")

        n_surveys_pc = pc_data[self.pc_data_params.visit_index_col].max()
        if isinstance(n_surveys_pc, int):
            n_surveys_pc = n_surveys_pc + 1
        else:
            raise ValueError("No visits found in PC data or n_visits is not an int")

        shape = (n_species, n_years, n_sites, n_surveys_pc)

        # y_pc is the count data
        y_pc = np.empty(shape, dtype=np.float32)
        y_pc.fill(np.nan)

        # date_pc is the julian date of the survey
        date_pc = np.empty(shape, dtype=np.float32)
        date_pc.fill(np.nan)

        # time_pc is the time of day of the survey
        time_pc = np.empty(shape, dtype=np.float32)
        time_pc.fill(np.nan)

        for row in pc_data.iter_rows(named=True):
            species_index = row["species_index"]
            year_index = row["year_index"]
            site_index = row["point_index"]
            visit_index = row[self.pc_data_params.visit_index_col]
            count = row[self.pc_data_params.count_col]
            datetime = row[self.pc_data_params.datetime_col]

            y_pc[species_index, year_index, site_index, visit_index] = count

            date_pc[species_index, year_index, site_index, visit_index] = (
                datetime.timetuple().tm_yday
            )

            time_pc[species_index, year_index, site_index, visit_index] = (
                datetime.hour * 60 + datetime.minute
            )

        # y_ind is the binary count data
        y_ind = np.where(y_pc > 0, 1, 0)

        return {
            "y_pc": y_pc,
            "y_ind": y_ind,
            "n_species": n_species,
            "n_years": n_years,
            "n_sites": n_sites,
            "n_surveys_pc": n_surveys_pc,
            "date_pc": date_pc,
            "time_pc": time_pc,
        }

    def _assign_indices(self, pc_data: pl.DataFrame) -> pl.DataFrame:
        """
        Assigns indices to the Point Count data.

        Each species, point, and year will have a unique index.

        Returns:
            pl.DataFrame: The Point Count data with an index assigned to each point.
        """
        point_col = self.pc_data_params.point_col
        datetime_col = self.pc_data_params.datetime_col
        species_col = self.pc_data_params.species_col
        visit_index_col = self.pc_data_params.visit_index_col

        # filter out null values for point and datetime
        pc_data = pc_data.filter(
            (pl.col(point_col).is_not_null()) & (pl.col(datetime_col).is_not_null())
        )

        if self.pc_data_params.years:
            pc_data = pc_data.filter(pl.col(datetime_col).dt.year().is_in(self.pc_data_params.years))

        if self.pc_data_params.species:
            pc_data = pc_data.filter(pl.col(species_col).is_in(self.pc_data_params.species))

        # sort by species and datetime for consistency when assigning indices
        pc_data = pc_data.sort([species_col, datetime_col, point_col])

        pc_data = pc_data.with_columns(
            year=pl.col(datetime_col).dt.year(),
        )

        # assign point, year, and species indices if not provided
        if (
            self.pc_data_params.point_index is not None
            and self.pc_data_params.species_index is not None
            and self.pc_data_params.year_index is not None
        ):
            pc_data = pc_data.with_columns(
                point_index=pl.col(point_col).replace_strict(self.pc_data_params.point_index),
                species_index=pl.col(species_col).replace_strict(self.pc_data_params.species_index),
                year_index=pl.col("year").replace_strict(self.pc_data_params.year_index),
            )
        elif (
            self.pc_data_params.point_index is not None
            or self.pc_data_params.species_index is not None
            or self.pc_data_params.year_index is not None
        ):
            raise ValueError(
                "If one of point_index, species_index, or year_index is provided, all three must be provided."
            )
        else:
            pc_data = pc_data.with_columns(
                point_index=pl.col(point_col).rank("dense").cast(pl.Int32) - 1,
                species_index=pl.col(species_col).rank("dense").cast(pl.Int32) - 1,
                year_index=pl.col("year").rank("dense").cast(pl.Int32) - 1,
            )

        # the visit index is included in the PC data dataframe, it is from the database
        # and is not assigned here

        # filter to the visit limit
        pc_data = pc_data.filter(pl.col(visit_index_col) < self.pc_data_params.pc_visit_limit)

        return pc_data
