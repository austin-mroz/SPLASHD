import argparse
import pathlib

import atomlite
import stk

import splashd


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    splashd_cage_database = atomlite.Database("cages.db")
    cage_database_entries = []

    cages = {
        splashd.cage.Cage(
            name="cc3",
            bb1_smiles="Brc1cc(Br)cc(Br)c1",
            bb2_smiles="Br/C=N/[C@H]1CCCC[C@@H]1/N=C/Br",
            topology="4+6",
        ): "cc3",
        splashd.cage.Cage(
            name="cc13",
            bb1_smiles="Brc1cc(Br)cc(Br)c1",
            bb2_smiles="C[C@@](C)(C/N=C/Br)/N=C/Br",
            topology="4+6",
        ): "cc13",
    }

    cage_dict = _make_pristine_cages(cages, args.output_directory)

    scrambled_cages = splashd.construct.make_scrambled_stk_cages(
        bb1_smiles="Brc1cc(Br)cc(Br)c1",
        bb2_smiles="Br/C=N/[C@H]1CCCC[C@@H]1/N=C/Br",
        bb3_smiles="C[C@@](C)(C/N=C/Br)/N=C/Br",
    )

    for scrambled_cage in scrambled_cages:
        scrambled_cage_directory = (
            args.output_directory / f"cc313_{scrambled_cages[scrambled_cage]}"
        )
        scrambled_cage_directory.mkdir(parents=True, exist_ok=True)

        splashd.construct.write_stk_cage(
            cage_name=f"cc313_{scrambled_cages[scrambled_cage]}",
            cage=scrambled_cage,
            output_directory=scrambled_cage_directory,
        )
        cage_dict[f"{scrambled_cages[scrambled_cage]}"] = scrambled_cage

    for cage_name in cage_dict:
        cage = cage_dict[cage_name]
        conformer_xyz_directory = splashd.simulate.optimize_geometry(
            system=cage,
            name=cage_name,
            output_directory=args.output_directory / f"{cage_name}",
        )

        cage_database_entries.append(
            atomlite.Entry.from_rdkit(
                key=f"{cage_name}",
                molecule=cage_dict[cage_name].to_rdkit_mol(),
                properties={
                    "conformer_xyz_directory": str(conformer_xyz_directory),
                },
            )
        )
    splashd_cage_database.add_entries(cage_database_entries)
    splashd_cage_database.connection.commit()


def _make_pristine_cages(
    cages: dict[splashd.cage.Cage, str],
    output_directory: pathlib.Path,
) -> dict[str, stk.ConstructedMolecule]:
    cage_dict = {}
    for cage in cages:
        constructed_cage = splashd.construct.make_stk_cage(cage)

        cage_directory = output_directory / f"{cages[cage]}"
        cage_directory.mkdir(parents=True, exist_ok=True)

        splashd.construct.write_stk_cage(
            cage_name=f"{cages[cage]}",
            cage=constructed_cage,
            output_directory=cage_directory,
        )
        cage_dict[f"{cages[cage]}"] = constructed_cage
    return cage_dict


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "1_output",
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
