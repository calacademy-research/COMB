from .single_year_single_species_all import SingleYearSingleSpeciesAll
from .single_year_single_year_no_scores import SingleYearSingleSpeciesNoScores
from .model_iterface import CombinedModelInterface
from typing import Literal

model_class_map = {
    "single_species_single_year_all": SingleYearSingleSpeciesAll,
    "single_species_single_year_no_scores": SingleYearSingleSpeciesNoScores,
}


ModelNames = Literal[
    "single_species_single_year_all", "single_species_single_year_no_scores"
]


def get_model_by_name(name: ModelNames) -> CombinedModelInterface:
    if name in model_class_map:
        return model_class_map[str(name)]

    raise ValueError(f"cannot find model called {name}")
