import argparse
import math
import pathlib
from dataclasses import dataclass

import atomlite
import pandas as pd
import rdkit

import maple


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    solvent_db = atomlite.Database(f"{args.solvent_database}")
    solvent_db_entries = []

    for entry in solvent_db.get_entries():
        print(entry.key)

        solvent_output_directory = args.output_directory / f"{entry.key}"
        solvent_output_directory.mkdir(parents=True, exist_ok=True)

        solvent_xtb_SP_directory = solvent_output_directory / "xtb_SPs"
        solvent_xtb_SP_directory.mkdir(parents=True, exist_ok=True)

        # perform xtb single point calculations for each conformer and
        # calculate the weighted all-atom steric parameters
        waa_ion = _run_xtb_SPs(
            solvent_name=entry.key,
            conformer_directory=pathlib.Path(
                entry.properties["conformer_xyz_directory"]
            ),
            output_directory=solvent_xtb_SP_directory,
            level="normal",
        )

        conformer_energy_csv_path = (
            args.output_directory / f"{entry.key}" / "conformer_relative_energies.csv"
        )
        waa_ion.xtb_SPs.to_csv(str(conformer_energy_csv_path))

        updated_entry = atomlite.PropertyEntry(
            key=entry.key,
            properties={
                "solvent_xtb_SP_directory": str(solvent_xtb_SP_directory),
                "minimum_conformer_energy_path": str(
                    waa_ion.minimum_conformer_energy_path
                ),
                "ion_conformer_relative_energies_csv": str(conformer_energy_csv_path),
            },
        )
        solvent_db.update_properties(updated_entry)
    solvent_db.connection.commit()


@dataclass(frozen=True, slots=True)
class WAASolvent:
    solvent_name: str
    conformer_directory: pathlib.Path
    minimum_conformer_energy_path: pathlib.Path
    xtb_SP_energy_directory: pathlib.Path
    xtb_SPs: pd.DataFrame


@dataclass(frozen=True, slots=True)
class AASolventParameters:
    L: list[float]
    B1: list[float]
    B5: list[float]


def _run_xtb_SPs(
    solvent_name: str,
    conformer_directory: pathlib.Path,
    output_directory: pathlib.Path,
    level: str,
) -> "WAASolvent":
    conformer_energy_dict = {}

    for conformer_path in conformer_directory.glob("*.xyz"):
        xtb_SP_directory = (
            output_directory / f"{conformer_path.name.removesuffix('_.xyz')}"
        )
        if xtb_SP_directory.is_dir():
            conformer_energy = maple.simulate.get_conformer_SP_energy(xtb_SP_directory)
        else:
            # get SP energy from xtb calculation
            conformer_energy = maple.simulate.run_xtb_single_point_energy(
                molecule=rdkit.Chem.rdmolfiles.MolFromXYZFile(str(conformer_path)),
                output_directory=xtb_SP_directory,
                level=level,
            )

        # store conformer energy
        conformer_energy_dict[conformer_path.name] = conformer_energy

    # convert conformer energy dictionary to DataFrame
    conformer_df = pd.DataFrame(
        conformer_energy_dict.items(),
        columns=[
            "conformer",
            "SP_energy_hartree",
        ],
    )

    # get the minimum conformer energy
    minimum_conformer_energy = (
        conformer_df.loc[:, ["SP_energy_hartree"]].min(axis=0).values[0]
    )

    # get minimum conformer energy id -- this is the filename
    minimum_conformer_energy_filename = conformer_df[
        conformer_df.SP_energy_hartree == conformer_df.SP_energy_hartree.min()
    ]["conformer"].values[0]

    # compute and append relative energy in Jmol
    conformer_df["relative_energy_Jmol"] = conformer_df.apply(
        lambda row: _get_relative_energy(row, minimum_conformer_energy),
        axis=1,
    )

    # compute and append energy distribution value
    conformer_df["distribution_value"] = conformer_df.apply(
        lambda row: _get_distribution_value(row),
        axis=1,
    )

    # compute and append scaled energies by their distribution value
    sum_distribution = conformer_df["distribution_value"].sum()
    conformer_df["scaled_distribution_value"] = conformer_df.apply(
        lambda row: _get_scaled_distribution_value(row, sum_distribution),
        axis=1,
    )

    # now, we are going to calculate the weighted steric parameters for each
    # scaled-all-atom parameter. This involves element-wise addition of the
    # scaled conformer steric parameters
    # This done in the return statement.

    return WAASolvent(
        solvent_name=solvent_name,
        conformer_directory=conformer_directory,
        minimum_conformer_energy_path=conformer_directory
        / minimum_conformer_energy_filename,
        xtb_SP_energy_directory=xtb_SP_directory,
        xtb_SPs=conformer_df,
    )


def _get_relative_energy(row: pd.DataFrame, minimum_conformer_energy: float) -> float:
    # return row["SP_energy"] / minimum_conformer_energy
    hartree_to_Jmol = 2600 * 1000
    return (row["SP_energy_hartree"] - minimum_conformer_energy) * hartree_to_Jmol


def _get_distribution_value(row: pd.DataFrame) -> float:
    temperature = 300  # K
    R = 8.314  # J K-1 mol-1
    return (
        round(math.exp((-1 * row["relative_energy_Jmol"]) / (R * temperature)) * 10000)
        / 10000
    )


def _get_scaled_distribution_value(row: pd.DataFrame, sum_distribution: float) -> float:
    return row["distribution_value"] / sum_distribution


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-s",
        "--solvent_database",
        help=(
            'An atomlite database file with properties: "xyz_file"',
            '"xtb_optimization_path", "xtb_md_path", "conformer_directory"',
            '"windows".',
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "solvents.db",
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "4_output",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
