from .single_year_single_species_no_pc import SingleYearSingleSpeciesNoPC
from .single_year_single_species_no_aru import SingleYearSingleSpeciesNoARU
from .single_year_single_species_no_aru_no_pc import SingleYearSingleSpeciesNoARUNoPC
from .single_year_single_species_no_scores_no_pc import (
    SingleYearSingleSpeciesNoScoresNoPC,
)
from .single_year_single_species_all import SingleYearSingleSpeciesAll
from .single_year_single_year_no_scores import SingleYearSingleSpeciesNoScores
from .single_year_single_species_all_marg import SingleYearSingleSpeciesAllMarg
from .single_year_single_species_datetime_all import SingleYearSingleSpeciesAllDateTime
from .single_year_single_species_no_scores_no_aru import (
    SingleYearSingleSpeciesNoScoresNoARU,
)
from .model_iterface import CombinedModelInterface
from typing import Literal

model_class_map = {
    "single_species_single_year_all": SingleYearSingleSpeciesAll,
    "single_species_single_year_all_marg": SingleYearSingleSpeciesAllMarg,
    "single_species_single_year_no_scores": SingleYearSingleSpeciesNoScores,
    "single_species_single_year_all_datetime": SingleYearSingleSpeciesAllDateTime,
    "single_year_single_species_no_scores_no_aru": SingleYearSingleSpeciesNoScoresNoARU,
    "single_year_single_species_no_scores_no_pc": SingleYearSingleSpeciesNoScoresNoPC,
    "single_year_single_species_no_aru_no_pc": SingleYearSingleSpeciesNoARUNoPC,
    "single_year_single_species_no_aru": SingleYearSingleSpeciesNoARU,
    "single_year_single_species_no_pc": SingleYearSingleSpeciesNoPC,
}


ModelNames = Literal[
    "single_species_single_year_all",
    "single_species_single_year_no_scores",
    "single_species_single_year_all_marg",
    "single_species_single_year_all_datetime",
    "single_year_single_species_no_scores_no_aru",
    "single_year_single_species_no_scores_no_pc",
    "single_year_single_species_no_aru_no_pc",
    "single_year_single_species_no_aru",
    "single_year_single_species_no_pc",
]


def get_model_by_name(name: ModelNames) -> CombinedModelInterface:
    if name in model_class_map:
        return model_class_map[str(name)]

    raise ValueError(f"cannot find model called {name}")
