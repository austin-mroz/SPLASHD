import argparse
import pathlib

import atomlite
import stk

import splashd


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    database = atomlite.Database("solvents.db")

    solvents = {
        "OC/C=C/c1ccccc1": "2E-3-Phenyl-2-propen-1-ol",
        "CC(C)(O)Cc1ccccc1": "2-Methyl-1-phenyl-2-propanol",
        "Cc1ccccc1O": "2-Methylphenol",
        "Cc1ccc(S)cc1": "4-Methylbenzenethiol",
        "CCc1ccccc1O": "2-Ethylphenol",
        "CC(O)c1ccccc1": "1-Phenylethanol",
        "Cc1cccc(C)c1O": "2,6-Dimethylphenol",
        "CCC(C)(O)CCc1ccccc1": "3-Methyl-1-phenyl-3-pentanol",
        "COc1cc(C)cc(OC)c1O": "2,6-Dimethoxy-4-methylphenol",
        "CC(CO)c1ccccc1": "2-Phenyl-1-propanol",
        "c2ccc(Oc1ccccc1)cc2": "diphenyl_ether",
    }

    database_entries = []

    for solvent in solvents:
        print(solvents[solvent])
        solvent_directory = args.output_directory / f"{solvents[solvent]}"

        if not solvent_directory.exists():
            solvent_directory.mkdir(parents=True, exist_ok=True)
            conformer_xyz_directory = splashd.simulate.optimize_geometry(
                system=stk.BuildingBlock.init_from_rdkit_mol(
                    splashd.construct.rdkit_from_smiles(solvent)
                ),
                name=f"{solvents[solvent]}",
                output_directory=solvent_directory,
            )
        else:
            conformer_xyz_directory = solvent_directory / "conformer_xyzs"
        database_entries.append(
            atomlite.Entry.from_rdkit(
                key=solvents[solvent],
                molecule=splashd.construct.rdkit_from_smiles(solvent),
                properties={
                    "smiles": solvent,
                    "conformer_xyz_directory": str(conformer_xyz_directory),
                },
            )
        )

    database.add_entries(database_entries)
    database.connection.commit()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "3_output",
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
