from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Cage:
    """
    Cage parameters

    """

    # the name of the cage
    name: str
    # smiles of building block 1
    bb1_smiles: str
    # smiles of building block 2
    bb2_smiles: str
    # topology
    topology: str
