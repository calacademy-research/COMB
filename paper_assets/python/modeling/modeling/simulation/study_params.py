from dataclasses import dataclass, asdict
import json
from models import model_zoo


@dataclass
class StudyParams:
    models: list[model_zoo.ModelNames]
    beta0: list[float]
    beta1: list[float]
    p11: list[float]
    p_aru11: list[float]
    p_aru01: list[float]
    mu: list[tuple[float, float]]
    sigma: list[tuple[float, float]]
    n_sites: list[int]
    n_surveys_pc: list[int]
    n_surveys_aru: list[int]
    n_surveys_scores: list[int]
    sim_name_for_data: model_zoo.ModelNames
    aru_scores_independent_model: bool

    def to_str(self) -> str:
        return json.dumps(asdict(self))

    @staticmethod
    def from_str(s: str):
        data = json.loads(s)
        return StudyParams(**data)
