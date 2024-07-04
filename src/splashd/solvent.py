from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True, slots=True)
class WAASolvent:
    solvent_name: str
    conformer_df: pd.DataFrame
    aa_solvent: "AASolvent"
    wL: list[float]
    wB1: list[float]
    wB5: list[float]


@dataclass(frozen=True, slots=True)
class AASolvent:
    numRotBonds: int
    numContiguousRotBonds: int
    method: str
    conformer_steric_parameters_dict: dict[str, "AASolventParameters"]


@dataclass(frozen=True, slots=True)
class AASolventParameters:
    L: list[float]
    B1: list[float]
    B5: list[float]
