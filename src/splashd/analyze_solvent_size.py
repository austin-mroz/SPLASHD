import pathlib
import typing
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
import rdkit
import rdkit.Chem.Lipinski
import smores
from scipy.spatial.distance import cdist

from .solvent import AASolvent, AASolventParameters, WAASolvent


def calculate_all_atom_steric_parameters(
    solvent_smiles: str,
    conformer_directory: pathlib.Path,
    minimum_conformer_energy_xyz_path: pathlib.Path,
    all_core_included: bool = False,
) -> "AASolvent":
    solvent_rdmol = rdkit.Chem.MolFromSmiles(solvent_smiles)

    if all_core_included:
        numRotBonds = 0
        numContiguousRotBonds = 0
        l_to_b5_ratio = 1.5
    else:
        numRotBonds = rdkit.Chem.rdMolDescriptors.CalcNumRotatableBonds(solvent_rdmol)
        numContiguousRotBonds = _calculate_num_contiguous_rotatable_bonds(solvent_rdmol)
        l_to_b5_ratio = _calculate_L_to_B5_ratio(conformer_directory)
        print(l_to_b5_ratio)
    conformer_steric_parameters_dict = {}
    l_values = []
    b1_values = []
    b5_values = []
    if numContiguousRotBonds > 2 and l_to_b5_ratio > 1.10:
        method = "core-excluded"
        num_atoms_to_include = 5
        for conformer_xyz_path in conformer_directory.glob("*.xyz"):
            smores_solvent = smores.Molecule.from_xyz_file(
                conformer_xyz_path,
                dummy_index=0,
                attached_index=1,
            )

            hdist = cdist(
                smores_solvent._positions,
                smores_solvent._positions,
                metric="euclidean",
            )
            farthest_pair_idx = np.unravel_index(hdist.argmax(), hdist.shape)

            num_atoms = smores_solvent._positions.shape[0]

            for attached_index in farthest_pair_idx:
                # sort hdist entry for farthest_pair_idx in aschending order
                # and remove the index of the atom distances from itself (0)
                core_excluded_idx = np.argsort(hdist[attached_index])[
                    num_atoms_to_include:
                ]
                l_atoms = list(range(num_atoms))
                l_atoms.remove(attached_index)
                for cei in list(core_excluded_idx):
                    l_atoms.remove(cei)
                for dummy_index in l_atoms:
                    smores_solvent_paramters = _get_smores_parameters_by_atom(
                        conformer_xyz=conformer_xyz_path,
                        attached_index=attached_index,
                        dummy_index=dummy_index,
                        excluded_indices=core_excluded_idx,
                    )
                    l_values.append(smores_solvent_paramters.L)
                    b1_values.append(smores_solvent_paramters.B1)
                    b5_values.append(smores_solvent_paramters.B5)
            conformer_steric_parameters_dict[conformer_xyz_path.name] = (
                AASolventParameters(
                    L=l_values,
                    B1=b1_values,
                    B5=b5_values,
                )
            )

    else:
        method = "core_included"
        for conformer_xyz_path in conformer_directory.glob("*.xyz"):
            smores_solvent = smores.Molecule.from_xyz_file(
                conformer_xyz_path,
                dummy_index=0,
                attached_index=1,
            )

            hdist = cdist(
                smores_solvent._positions,
                smores_solvent._positions,
                metric="euclidean",
            )
            farthest_pair_idx = np.unravel_index(hdist.argmax(), hdist.shape)

            num_atoms = smores_solvent._positions.shape[0]

            for attached_index in farthest_pair_idx:
                l_atoms = list(range(num_atoms))
                l_atoms.remove(attached_index)
                for dummy_index in l_atoms:
                    smores_solvent_paramters = _get_smores_parameters_by_atom(
                        conformer_xyz=conformer_xyz_path,
                        attached_index=attached_index,
                        dummy_index=dummy_index,
                    )
                    l_values.append(smores_solvent_paramters.L)
                    b1_values.append(smores_solvent_paramters.B1)
                    b5_values.append(smores_solvent_paramters.B5)
            conformer_steric_parameters_dict[conformer_xyz_path.name] = (
                AASolventParameters(
                    L=l_values,
                    B1=b1_values,
                    B5=b5_values,
                )
            )
    return AASolvent(
        numRotBonds=numRotBonds,
        numContiguousRotBonds=_calculate_num_contiguous_rotatable_bonds(solvent_rdmol),
        method=method,
        conformer_steric_parameters_dict=conformer_steric_parameters_dict,
    )


def _calculate_num_contiguous_rotatable_bonds(
    solvent_rdmol,
) -> int:
    # Find groups of contiguous rotatable bonds in mol
    bond_groups = find_bond_groups(solvent_rdmol)
    # As bond groups are sorted by decreasing size, the size of the first group (if any)
    # is the largest number of contiguous rotatable bonds in mol
    largest_n_cont_rot_bonds = len(bond_groups[0]) if bond_groups else 0
    return largest_n_cont_rot_bonds


def find_bond_groups(mol):
    """Find groups of contiguous rotatable bonds and return them sorted by decreasing size"""
    rot_atom_pairs = mol.GetSubstructMatches(rdkit.Chem.Lipinski.RotatableBondSmarts)
    rot_bond_set = set([mol.GetBondBetweenAtoms(*ap).GetIdx() for ap in rot_atom_pairs])
    rot_bond_groups = []
    while rot_bond_set:
        i = rot_bond_set.pop()
        connected_bond_set = set([i])
        stack = [i]
        while stack:
            i = stack.pop()
            b = mol.GetBondWithIdx(i)
            bonds = []
            for a in (b.GetBeginAtom(), b.GetEndAtom()):
                bonds.extend(
                    [
                        b.GetIdx()
                        for b in a.GetBonds()
                        if (
                            (b.GetIdx() in rot_bond_set)
                            and (not (b.GetIdx() in connected_bond_set))
                        )
                    ]
                )
            connected_bond_set.update(bonds)
            stack.extend(bonds)
        rot_bond_set.difference_update(connected_bond_set)
        rot_bond_groups.append(tuple(connected_bond_set))
    return tuple(sorted(rot_bond_groups, reverse=True, key=lambda x: len(x)))


def _calculate_L_to_B5_ratio(
    conformer_directory: pathlib.Path,
) -> float:
    all_lb5 = []
    for conformer_xyz_path in conformer_directory.glob("*.xyz"):
        smores_solvent = smores.Molecule.from_xyz_file(
            conformer_xyz_path,
            dummy_index=0,
            attached_index=1,
        )

        hdist = cdist(
            smores_solvent._positions,
            smores_solvent._positions,
            metric="euclidean",
        )
        farthest_pair_idx = np.unravel_index(hdist.argmax(), hdist.shape)

        num_atoms = smores_solvent._positions.shape[0]

        l_values = []
        b1_values = []
        b5_values = []

        for attached_index in farthest_pair_idx:
            l_atoms = list(range(num_atoms))
            l_atoms.remove(attached_index)
            for dummy_index in l_atoms:
                smores_solvent_paramters = _get_smores_parameters_by_atom(
                    conformer_xyz=conformer_xyz_path,
                    attached_index=attached_index,
                    dummy_index=dummy_index,
                )
                l_values.append(smores_solvent_paramters.L)
                b1_values.append(smores_solvent_paramters.B1)
                b5_values.append(smores_solvent_paramters.B5)
        all_lb5.append(np.mean(np.array(l_values)) / np.mean(np.array(b5_values)))
    return np.mean(all_lb5)


def calculate_weighted_all_atom_steric_parameters(
    solvent_name: str,
    aa_solvent: "AASolvent",
    conformer_relative_energies_csv: pathlib.Path,
    conformer_steric_parameters_dict: dict[str, "AASolventParameters"],
) -> "WAASolvent":
    # convert conformer energy dictionary to DataFrame
    conformer_df = pd.read_csv(conformer_relative_energies_csv)

    # scale and append L values
    conformer_df["scaled_L_values"] = conformer_df.apply(
        lambda row: _get_scaled_L_values(row, conformer_steric_parameters_dict),
        axis=1,
    )

    # scale and append B1 values
    conformer_df["scaled_B1_values"] = conformer_df.apply(
        lambda row: _get_scaled_B1_values(row, conformer_steric_parameters_dict),
        axis=1,
    )

    # scale and append B5 values
    conformer_df["scaled_B5_values"] = conformer_df.apply(
        lambda row: _get_scaled_B5_values(row, conformer_steric_parameters_dict),
        axis=1,
    )

    # now, we are going to calculate the weighted steric parameters for each
    # scaled-all-atom parameter. This involves element-wise addition of the
    # scaled conformer steric parameters
    # This done in the return statement.

    return WAASolvent(
        solvent_name=solvent_name,
        conformer_df=conformer_df,
        aa_solvent=aa_solvent,
        wL=conformer_df["scaled_L_values"].to_numpy().sum(axis=0),
        wB1=conformer_df["scaled_B1_values"].to_numpy().sum(axis=0),
        wB5=conformer_df["scaled_B5_values"].to_numpy().sum(axis=0),
    )


def _get_scaled_L_values(
    row: pd.DataFrame,
    conformer_steric_parameters_dict: dict[str, str],
) -> float:
    steric_parameters = conformer_steric_parameters_dict[row["conformer"]]
    return np.array(
        [row["scaled_distribution_value"] * l_value for l_value in steric_parameters.L]
    )


def _get_scaled_B1_values(
    row: pd.DataFrame,
    conformer_steric_parameters_dict: dict[str, str],
) -> float:
    steric_parameters = conformer_steric_parameters_dict[row["conformer"]]
    return np.array(
        [row["scaled_distribution_value"] * l_value for l_value in steric_parameters.B1]
    )


def _get_scaled_B5_values(
    row: pd.DataFrame,
    conformer_steric_parameters_dict: dict[str, str],
) -> float:
    steric_parameters = conformer_steric_parameters_dict[row["conformer"]]
    return np.array(
        [row["scaled_distribution_value"] * l_value for l_value in steric_parameters.B5]
    )


@dataclass(frozen=True, slots=True)
class XyzData:
    elements: list[str]
    coordinates: npt.NDArray
    dummy_index: int
    attached_index: int


def _get_smores_parameters_by_atom(
    conformer_xyz: pathlib.Path,
    attached_index: int,
    dummy_index: int,
    excluded_indices: typing.Iterable[int] | None = None,
) -> smores.StericParameters:
    smores_solvent = smores.Molecule.from_xyz_file(
        path=conformer_xyz,
        dummy_index=dummy_index,
        attached_index=attached_index,
        excluded_indices=excluded_indices,
    )
    return smores_solvent.get_steric_parameters()


def _get_smores_LV_parameters(conformer_xyz: pathlib.Path) -> smores.StericParameters:
    sorted_xyz_data = _sort_xyz_data_for_LV(conformer_xyz)
    smores_solvent = smores.Molecule(
        atoms=sorted_xyz_data.elements,
        positions=sorted_xyz_data.coordinates,
        dummy_index=sorted_xyz_data.dummy_index,  # .coordinates.shape[0] - 1,
        attached_index=sorted_xyz_data.attached_index,  # coordinates.shape[0] - 2,
    )
    smores_parameters = smores_solvent.get_steric_parameters()
    return smores_parameters


def _get_smores_COM_parameters(conformer_xyz: pathlib.Path) -> smores.StericParameters:
    sorted_xyz_data = _sort_xyz_data_for_COM(conformer_xyz)
    smores_solvent = smores.Molecule(
        atoms=sorted_xyz_data.elements,
        positions=sorted_xyz_data.coordinates,
        dummy_index=sorted_xyz_data.dummy_index,
        attached_index=sorted_xyz_data.attached_index,
    )
    return smores_solvent.get_steric_parameters()


def _sort_xyz_data_for_LV(xyz_path: pathlib.Path) -> XyzData:  # dict[str, float]:
    smores_solvent = smores.Molecule.from_xyz_file(
        xyz_path,
        dummy_index=0,
        attached_index=1,
    )
    smores_solvent_center_of_mass = _get_center_of_mass(
        elements=smores_solvent._atoms, coordinates=smores_solvent._positions
    )
    distances = []
    for element, coordinate in zip(smores_solvent._atoms, smores_solvent._positions):
        distances.append(
            math.dist(
                coordinate.reshape(-1, 1), smores_solvent_center_of_mass.reshape(-1, 1)
            )
        )
    element_column = np.array(smores_solvent._atoms).T
    x_coordinate_column = np.array(smores_solvent._positions[:, 0]).T
    y_coordinate_column = np.array(smores_solvent._positions[:, 1]).T
    z_coordinate_column = np.array(smores_solvent._positions[:, 2]).T
    distances = np.array(distances).T

    xyz_data = np.vstack(
        (
            element_column,
            x_coordinate_column,
            y_coordinate_column,
            z_coordinate_column,
            distances,
        )
    ).T

    sorted_indices = np.argsort(xyz_data[:, -1])
    sorted_xyz_data = xyz_data[sorted_indices]

    sorted_elements = sorted_xyz_data[:, 0]
    sorted_coordinates = sorted_xyz_data[:, 1:4].astype(float)

    distance_matrix = scipy.spatial.distance_matrix(
        sorted_coordinates,
        sorted_coordinates,
    )

    smores_indices = np.unravel_index(np.argmax(distance_matrix), distance_matrix.shape)

    return XyzData(
        elements=sorted_elements,
        coordinates=sorted_coordinates,
        dummy_index=smores_indices[1],
        attached_index=smores_indices[0],
    )


def _sort_xyz_data_for_COM(xyz_path: pathlib.Path) -> XyzData:  # dict[str, float]:
    smores_solvent = smores.Molecule.from_xyz_file(
        xyz_path,
        dummy_index=0,
        attached_index=1,
    )
    smores_solvent_center_of_mass = _get_center_of_mass(
        elements=smores_solvent._atoms, coordinates=smores_solvent._positions
    )
    distances = []
    for element, coordinate in zip(smores_solvent._atoms, smores_solvent._positions):
        distances.append(
            math.dist(
                coordinate.reshape(-1, 1), smores_solvent_center_of_mass.reshape(-1, 1)
            )
        )
    element_column = np.array(smores_solvent._atoms).T
    x_coordinate_column = np.array(smores_solvent._positions[:, 0]).T
    y_coordinate_column = np.array(smores_solvent._positions[:, 1]).T
    z_coordinate_column = np.array(smores_solvent._positions[:, 2]).T
    distances = np.array(distances).T

    xyz_data = np.vstack(
        (
            element_column,
            x_coordinate_column,
            y_coordinate_column,
            z_coordinate_column,
            distances,
        )
    ).T

    sorted_indices = np.argsort(xyz_data[:, -1])
    sorted_xyz_data = xyz_data[sorted_indices]

    sorted_elements = sorted_xyz_data[:, 0]
    sorted_coordinates = sorted_xyz_data[:, 1:4].astype(float)

    return XyzData(
        elements=sorted_elements,
        coordinates=sorted_coordinates,
        dummy_index=sorted_coordinates.shape[0] - 1,
        attached_index=0,
    )


def _get_center_of_mass(
    elements: typing.Iterable[str],
    coordinates: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    atom_masses = np.array([atomic_mass[element] for element in elements])
    scaled_coordinates = coordinates * atom_masses[:, np.newaxis]
    return scaled_coordinates.sum(axis=0) / atom_masses.sum()
