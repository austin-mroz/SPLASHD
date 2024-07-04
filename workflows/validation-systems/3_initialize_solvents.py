import argparse
import pathlib

import atomlite
import stk

import maple


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    database = atomlite.Database("solvents.db")
    solvents = {
        "ClC(Cl)Cl": "chloroform",
        "Cl/C(Cl)=C(Cl)\C(Cl)(Cl)Cl": "hexachloropropene",
        "ClCc1ccc(Cl)cc1Cl": "24-dichlorobenzyl_chloride",
        "Cc1ccc(Cl)cc1Cl": "24-dichlorotoluene",
        "OCc1ccc(OC(F)(F)F)cc1": "4-trifluoromethoxybenzyl_alcohol",
        "CC(=O)c1ccccc1O": "2-hydroxyacetophenone",
        "COC(=O)c1ccccc1O": "methyl_salicylate",
        "O=C(OCc1ccccc1)c2ccccc2": "benzyl_benzoate",
        "Oc1c(F)c(F)c(C(F)(F)F)c(F)c1F": "2356-tetrafluoro-4-trifluoromethylphenol",
        "Oc1cc(C(F)(F)F)cc(C(F)(F)F)c1": "35-bistrifluoromethyl_phenol",
        "Oc1ccc(C(F)(F)F)cc1F": "2-fluoro-5-trifluoromethyl_phenol",
        "ClCCl": "dichloromethane",
        "OC(C(F)(F)F)(C(F)(F)F)C(O)(C(F)(F)F)C(F)(F)F": "hexafluoro-23-bistrifluoromethyl_butane-23-diol",
        "C1COCCOCCOCCOCCO1": "15-crown-5",
        "ClC(Cl)C(Cl)(Cl)C(Cl)Cl": "112233-hexachloropropane",
        "Cl/C(Cl)=C(Cl)/Cl": "perchloroethylene",
    }

    database_entries = []

    for solvent in solvents:
        print(solvents[solvent])
        solvent_directory = args.output_directory / f"{solvents[solvent]}"

        if not solvent_directory.exists():
            solvent_directory.mkdir(parents=True, exist_ok=True)
            conformer_xyz_directory = maple.simulate.optimize_geometry(
                system=stk.BuildingBlock.init_from_rdkit_mol(
                    maple.construct.rdkit_from_smiles(solvent)
                ),
                name=f"{solvents[solvent]}",
                output_directory=solvent_directory,
            )
        else:
            conformer_xyz_directory = solvent_directory / "conformer_xyzs"
        database_entries.append(
            atomlite.Entry.from_rdkit(
                key=solvents[solvent],
                molecule=maple.construct.rdkit_from_smiles(solvent),
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
