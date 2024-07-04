import pathlib
import typing

import rdkit.Chem.AllChem as rdkit
import stk

from .cage import Cage


def make_stk_cage(
    initial_cage_parameter: Cage,
) -> stk.ConstructedMolecule:
    bb1 = stk.BuildingBlock(
        smiles=initial_cage_parameter.bb1_smiles,
        functional_groups=[stk.BromoFactory()],
    )
    bb2 = stk.BuildingBlock(
        smiles=initial_cage_parameter.bb2_smiles,
        functional_groups=[stk.BromoFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=_get_cage_topology(
            bb1,
            bb2,
            initial_cage_parameter.topology,
        ),
    )

    return cage


def make_scrambled_stk_cages(
    bb1_smiles: str,
    bb2_smiles: str,
    bb3_smiles: str,
) -> dict[stk.ConstructedMolecule, str]:
    bb1 = stk.BuildingBlock(
        smiles=bb1_smiles,
        functional_groups=[stk.BromoFactory()],
    )
    bb2 = stk.BuildingBlock(
        smiles=bb2_smiles,
        functional_groups=[stk.BromoFactory()],
    )
    bb3 = stk.BuildingBlock(
        smiles=bb3_smiles,
        functional_groups=[stk.BromoFactory()],
    )

    bb2_bb3_partitions = _sorted_k_partitions(range(4, 10), 2)
    cages = {}
    cage_num = 0

    bb2_bb3_partitions = {
        "3-1_13-5": [(4,), (5, 6, 7, 8, 9)],
        "3-5_13-1": [(5, 6, 7, 8, 9), (4,)],
        "3-2_13-4a": [(4, 5), (6, 7, 8, 9)],
        "3-2_13-4b": [(4, 9), (5, 6, 7, 8)],
        "3-3_13-3a": [(4, 5, 6), (7, 8, 9)],
        "3-3_13-3b": [(4, 7, 9), (5, 6, 8)],
        "3-4_13-2a": [(4, 5), (6, 7, 8, 9)],
        "3-4_13-2b": [(4, 9), (5, 6, 7, 8)],
    }

    for partition in bb2_bb3_partitions:
        topology = stk.cage.FourPlusSix(
            building_blocks={
                bb1: (0, 1, 2, 3),
                bb2: bb2_bb3_partitions[partition][0],
                bb3: bb2_bb3_partitions[partition][1],
            },
            optimizer=stk.MCHammer(),
        )

        cages[
            stk.ConstructedMolecule(
                topology_graph=topology,
            )
        ] = f"{partition}"
        cage_num += 1
    return cages


def write_stk_cage(
    cage_name: str,
    cage: stk.ConstructedMolecule,
    output_directory: pathlib.Path,
) -> None:
    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)
    stk.XyzWriter().write(
        molecule=cage,
        path=str(output_directory / f"{cage_name}_starting_geometry.xyz"),
    )


def _sorted_k_partitions(seq, k):
    """Returns a list of all unique k-partitions of `seq`.

    Each partition is a list of parts, and each part is a tuple.

    The parts in each individual partition will be sorted in shortlex
    order (i.e., by length first, then lexicographically).

    The overall list of partitions will then be sorted by the length
    of their first part, the length of their second part, ...,
    the length of their last part, and then lexicographically.
    """
    n = len(seq)
    groups = []  # a list of lists, currently empty

    def generate_partitions(i):
        if i >= n:
            yield list(map(tuple, groups))
        else:
            if n - i > k - len(groups):
                for group in groups:
                    group.append(seq[i])
                    yield from generate_partitions(i + 1)
                    group.pop()

            if len(groups) < k:
                groups.append([seq[i]])
                yield from generate_partitions(i + 1)
                groups.pop()

    result = generate_partitions(0)

    # Sort the parts in each partition in shortlex order
    result = [sorted(ps, key=lambda p: (len(p), p)) for ps in result]
    # Sort partitions by the length of each part, then lexicographically.
    result = sorted(result, key=lambda ps: (*map(len, ps), ps))

    return result


def rdkit_from_smiles(smiles: str) -> rdkit.Mol:
    """
    Get an embeded :mod:`rdkit` molecule from SMILES.

    Examples:

        *Create an rdkit molecule*

        .. doctest:: create-an-rdkit-molecule

            >>> import smores
            >>> molecule = smores.rdkit_from_smiles("Br")
            >>> import rdkit.Chem.AllChem as rdkit
            >>> rdkit.MolToSmiles(molecule)
            '[H]Br'

    Parameters:
        smiles: The SMILES of the molecule.
    Returns:
        The :mod:`rdkit` molecule.

    """

    molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
    _optimize(molecule)
    rdkit.Kekulize(molecule)
    return molecule


def _optimize(molecule: rdkit.Mol) -> None:
    etkdg = rdkit.ETKDGv3()
    etkdg.randomSeed = 4
    rdkit.EmbedMolecule(molecule, etkdg)


def _get_cage_topology(
    bb1: stk.BuildingBlock,
    bb2: stk.BuildingBlock,
    topology: typing.Literal["4+6", "8+12"],
) -> stk.TopologyGraph:
    match topology:
        case "4+6":
            return stk.cage.FourPlusSix(
                (bb1, bb2),
                optimizer=stk.MCHammer(),
                vertex_alignments={0: 1, 1: 1, 2: 2},
            )
        case "8+12":
            return stk.cage.EightPlusTwelve((bb1, bb2), optimizer=stk.MCHammer())
        case _:
            raise ValueError("Not a valid topology.")
