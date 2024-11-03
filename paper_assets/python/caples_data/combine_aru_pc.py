import polars as pl
import numpy as np
from dataclasses import dataclass
from typing import List, Union, Tuple
from datetime import time

from read_cap_data import AruData, ARUDataParams, PCDataParams, PcData


@dataclass
class CombinedParams:
    years: List[int]
    aru_visit_limit: int = 24
    pc_visit_limit: int = 3

    aru_species_col: str = "species"
    aru_logit_col: str = "logit"
    aru_datetime_col: str = "datetime"
    aru_point_col: str = "point"

    aru_time_range: Tuple[time, time] = (time(6), time(9, 59))
    aru_threshold: float = 0.5
    aru_skip_first_day: bool = True
    aru_remove_call_type: bool = True

    pc_species_col: str = "species"
    pc_count_col: str = "count"
    pc_datetime_col: str = "datetime"
    pc_point_col: str = "point"
    pc_visit_index_col: str = "visit_index"

    species: Union[List[str], None] = None


class CombinedData:
    def __init__(self, aru_data: pl.DataFrame, pc_data: pl.DataFrame, params: CombinedParams):
        self.params = params
        self.aru_data = aru_data
        self.pc_data = pc_data

        self.pc_data = self.convert_PC_to_6_codes()

        if self.params.aru_remove_call_type:
            self.aru_data = self.remove_call_type()

        self.point_index, self.species_index, self.year_index = self._build_all_indices()
        self.aru_params = self._build_aru_params()
        self.pc_params = self._build_pc_params()

        self.aru = AruData(self.aru_data, self.aru_params)
        self.pc = PcData(self.pc_data, self.pc_params)

        self._verify_data()

        self.combined_data = self.combine_data()

    def remove_call_type(self):
        """
        Removes the call type on the species label column.

        In practice, just takes the first part before the underscore.

        We might have `westan_call` and `westan_song` as two different labels,
        but we want to treat them as the same species.

        We take the max of the logit values for labels that have the same
        species label after removing the call type.
        """

        aru_species_col = self.params.aru_species_col
        aru_datetime_col = self.params.aru_datetime_col
        aru_point_col = self.params.aru_point_col
        aru_logit_col = self.params.aru_logit_col

        species_no_call_type = (
            self.aru_data.with_columns(
                pl.col(aru_species_col).str.split("_").list.get(0).alias(aru_species_col)
            )
            .group_by(
                [aru_species_col, aru_datetime_col, aru_point_col],
            )
            .agg(pl.col(aru_logit_col).max().alias(aru_logit_col))
        )
        return species_no_call_type

    def convert_PC_to_6_codes(self):
        """
        Given the PC data, convert the species columns to the 6 letter codes
        (from the 4 letter codes).

        This is necessary because the ARU data uses the 6 letter codes.
        """

        # Default mapping for species that are not found in the bird_codes.csv file
        default_mapping = "NOCODEFOUND"

        pc_data = self.pc_data
        pc_species_col = self.params.pc_species_col

        codes = pl.read_csv("data/bird_codes.csv")
        code_mapping = dict(zip(codes["four_code"], codes["code"]))

        pc_mapped = pc_data.with_columns(
            pl.col(pc_species_col).alias("four_code"),
            pl.col(pc_species_col)
            .replace_strict(code_mapping, default=default_mapping)
            .alias(pc_species_col),
        )
        return pc_mapped

    def combine_data(self):
        """
        Combines the ARU and PC dicts into a single dict,
        combining fields that are common (and checking that they are the same)
        """
        if self.aru.aru_data_dict["n_sites"] != self.pc.pc_data_dict["n_sites"]:
            raise ValueError(
                f"""
            The number of sites in the ARU and PC data must be the same
            ARU n_sites: {self.aru.aru_data_dict["n_sites"]}, PC n_sites: {self.pc.pc_data_dict["n_sites"]}
            """
            )
        n_sites = self.aru.aru_data_dict["n_sites"]

        if self.aru.aru_data_dict["n_years"] != self.pc.pc_data_dict["n_years"]:
            raise ValueError(
                f"""
            The number of years in the ARU and PC data must be the same
            ARU n_years: {self.aru.aru_data_dict["n_years"]}, PC n_years: {self.pc.pc_data_dict["n_years"]}
            """
            )
        n_years = self.aru.aru_data_dict["n_years"]

        if self.aru.aru_data_dict["n_species"] != self.pc.pc_data_dict["n_species"]:
            raise ValueError(
                f"""
            The number of species in the ARU and PC data must be the same
            ARU n_species: {self.aru.aru_data_dict["n_species"]}, PC n_species: {self.pc.pc_data_dict["n_species"]}
            """
            )
        n_species = self.aru.aru_data_dict["n_species"]

        return {
            "n_sites": n_sites,
            "n_years": n_years,
            "n_species": n_species,
            "y_aru": self.aru.aru_data_dict["y_aru"],
            "y_pc": self.pc.pc_data_dict["y_pc"],
            "y_index": self.pc.pc_data_dict["y_ind"],
            "n_surveys_pc": self.pc.pc_data_dict["n_surveys_pc"],
            "date_pc": self.pc.pc_data_dict["date_pc"],
            "time_pc": self.pc.pc_data_dict["time_pc"],
            "date_aru": self.aru.aru_data_dict["date_aru"],
            "time_aru": self.aru.aru_data_dict["time_aru"],
            "n_surveys_aru": self.aru.aru_data_dict["n_surveys_aru"],
        }

    def _verify_data(self):
        """
        Verifies that the data is in the correct format
        """
        # Check that the first three (species, year, point) dimensions are the same length
        if isinstance(self.aru.aru_data_dict["y_aru"], np.ndarray):
            aru_shape = self.aru.aru_data_dict["y_aru"].shape
        if isinstance(self.pc.pc_data_dict["y_pc"], np.ndarray):
            pc_shape = self.pc.pc_data_dict["y_pc"].shape

        if (
            aru_shape[0] != pc_shape[0]
            or aru_shape[1] != pc_shape[1]
            or aru_shape[2] != pc_shape[2]
        ):
            raise ValueError(
                f"""
            The first three dimensions of the ARU and PC data must be the same length
            ARU shape: {aru_shape}, PC shape: {pc_shape}
            """
            )

    def _build_aru_params(self) -> ARUDataParams:
        return ARUDataParams(
            species_col=self.params.aru_species_col,
            logit_col=self.params.aru_logit_col,
            datetime_col=self.params.aru_datetime_col,
            point_col=self.params.aru_point_col,
            time_range=self.params.aru_time_range,
            threshold=self.params.aru_threshold,
            skip_first_day=self.params.aru_skip_first_day,
            aru_visit_limit=self.params.aru_visit_limit,
            point_index=self.point_index,
            species_index=self.species_index,
            year_index=self.year_index,
            years=self.params.years,
        )

    def _build_pc_params(self) -> PCDataParams:
        return PCDataParams(
            species_col=self.params.pc_species_col,
            count_col=self.params.pc_count_col,
            datetime_col=self.params.pc_datetime_col,
            point_col=self.params.pc_point_col,
            visit_index_col=self.params.pc_visit_index_col,
            pc_visit_limit=self.params.pc_visit_limit,
            point_index=self.point_index,
            species_index=self.species_index,
            year_index=self.year_index,
            years=self.params.years,
        )

    def _build_all_indices(self) -> Tuple[dict[int, int], dict[str, int], dict[int, int]]:
        """
        Builds the following indices:
        - point_index: maps point names to integers
        - species_index: maps species names to integers
        - year_index: maps years to integers

        Returns:
            Tuple of dictionaries
        """

        # Build point index
        # get all of the unique points from the ARU and PC data
        # and assign them an index
        points = (
            self.aru_data[self.params.aru_point_col]
            .cast(pl.Int64)
            .append(self.pc_data[self.params.pc_point_col].cast(pl.Int64))
            .unique()
            .to_list()
        )
        point_index = {point: i for i, point in enumerate(points)}

        # Build species index
        # get all of the unique species from the ARU and PC data
        # and assign them an index
        species = (
            self.aru_data[self.params.aru_species_col]
            .append(self.pc_data[self.params.pc_species_col])
            .unique()
            .to_list()
        )
        species_index = {species: i for i, species in enumerate(species)}

        # Build year index
        year_index = {year: i for i, year in enumerate(self.params.years)}

        return point_index, species_index, year_index
