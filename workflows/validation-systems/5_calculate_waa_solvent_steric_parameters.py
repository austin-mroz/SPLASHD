import argparse
import pathlib
from dataclasses import dataclass

import atomlite
import numpy as np
import pandas as pd
import rdkit
import rdkit.Chem.Lipinski
import smores

import maple


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    solvent_db = atomlite.Database(args.solvent_database)

    for entry in solvent_db.get_entries():
        print(entry.key)
        aa_solvent = maple.analyze_solvent_size.calculate_all_atom_steric_parameters(
            solvent_smiles=entry.properties["smiles"],
            conformer_directory=pathlib.Path(
                entry.properties["conformer_xyz_directory"]
            ),
            minimum_conformer_energy_xyz_path=pathlib.Path(
                entry.properties["minimum_conformer_energy_path"]
            ),
            all_core_included=False,
        )
        # perform xrb single point calculations for each conformer and
        # calculate the weighted all-atom steric parameters
        waa_solvent = maple.analyze_solvent_size.calculate_weighted_all_atom_steric_parameters(
            solvent_name=entry.key,
            aa_solvent=aa_solvent,
            conformer_relative_energies_csv=entry.properties[
                "ion_conformer_relative_energies_csv"
            ],
            conformer_steric_parameters_dict=aa_solvent.conformer_steric_parameters_dict,
        )

        updated_df_path = args.output_directory / f"{entry.key}_weighted_parameters.csv"
        waa_solvent.conformer_df.to_csv(updated_df_path)

        updated_entry = atomlite.PropertyEntry(
            key=entry.key,
            properties={
                "numRotBonds": aa_solvent.numRotBonds,
                "numContiguousRotBonds": aa_solvent.numContiguousRotBonds,
                "aa_method": aa_solvent.method,
                # "aa_conformer_steric_parameter_dict": dict(
                #     waa_solvent.aa_solvent.conformer_steric_parameters_dict
                # ),
                "weighted_steric_parameter_csv_path": str(updated_df_path),
                "waa_wL": list(waa_solvent.wL),
                "waa_wB1": list(waa_solvent.wB1),
                "waa_wB5": list(waa_solvent.wB5),
            },
        )
        solvent_db.update_properties(updated_entry)
    solvent_db.connection.commit()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-i",
        "--solvent_database",
        help="The database containing the molecular ion conformer directories",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "solvents.db",
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "5_output",
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
