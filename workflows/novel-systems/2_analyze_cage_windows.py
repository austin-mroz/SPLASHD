import argparse
import json
import pathlib
import sqlite3
import typing
from dataclasses import dataclass

import pywindow as pw


def main() -> None:
    args = _get_command_line_arguments()
    database = sqlite3.connect(args.database)
    cursor = database.cursor()

    cages = tuple(_get_cages(database))

    for cage in cages:
        print(cage.name)

        windows = _get_windows(
            pathlib.Path(cage.conformer_directory),
        )
        print(cage.conformer_directory)
        cursor.execute(
            "UPDATE molecules SET properties=json_patch(properties,?) WHERE key=?",
            (
                json.dumps({"windows": windows}),
                cage.name,
            ),
        )

        database.commit()

    database.close()


@dataclass(frozen=True)
class Cage:
    name: str
    conformer_directory: pathlib.Path


def _get_cages(
    database: sqlite3.Connection,
) -> typing.Iterator[Cage]:
    for key, _, properties in database.execute("SELECT * FROM molecules"):
        yield Cage(
            name=key,
            conformer_directory=pathlib.Path(
                json.loads(properties)["conformer_xyz_directory"]
            ),
        )


def _get_windows(
    conformer_directory: pathlib.Path,
) -> list:
    window_list = []
    for conformer_xyz in conformer_directory.glob("*.xyz"):
        print(conformer_xyz)
        pw_conformer = pw.MolecularSystem.load_file(conformer_xyz)
        pw_cage = pw_conformer.system_to_molecule()
        pw_cage.calculate_windows()

        conformer_window_list = list(pw_cage.properties["windows"]["diameters"])
        window_list.extend(conformer_window_list)

    return window_list


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d",
        "--database",
        help=(
            'An atomlite database file with properties: "core", "substituent", '
            '"xyz_file", "esp_file", dummy_index" and "attached_index".'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "cages.db",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
