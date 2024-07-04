import argparse
import pathlib

import atomlite
import numpy as np
import numpy.typing as npt

import maple


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    # load the solvents database
    solvent_db = atomlite.Database(args.solvent_database)
    cage_db = atomlite.Database(args.cage_database)

    cage_apertures = _get_combined_cage_apertures(cage_db)
    maple.plot.plot_cage_window_distribution(
        cage_name="combined_windows",
        cage_apertures=cage_apertures,
        output_directory=args.output_directory / "cage_window_distributions",
    )

    _plot_cage_window_distributions(
        cage_db=cage_db,
        output_directory=args.output_directory / "cage_window_distributions",
    )

    for solvent in solvent_db.get_entries():
        print(solvent.key)

        if solvent.properties["aa_method"] == "core-excluded":
            solvent_sizes = np.array(list(solvent.properties["waa_wL"])).T
        elif solvent.properties["aa_method"] == "core_included":
            solvent_sizes = np.array(list(solvent.properties["waa_wB5"])).T

        kde_overlap = maple.predict.calculate_kde_overlap(
            cage_apertures=cage_apertures,
            solvent_sizes=solvent_sizes,
        )

        maple.plot.plot_cage_distribution_with_solvent_means(
            method_name=solvent.properties["aa_method"],
            cage_name="scrambled_CC313",
            cage_apertures=cage_apertures,
            solvent_name=solvent.key,
            solvent_sizes=solvent_sizes,
            overlap=kde_overlap,
            output_directory=args.output_directory / "porous_liquid_mean_comparison",
        )

        maple.plot.plot_cage_distribution_with_solvent_distributions(
            method_name=solvent.properties["aa_method"],
            cage_name="scrambled_CC313",
            cage_apertures=cage_apertures,
            solvent_name=solvent.key,
            solvent_sizes=solvent_sizes,
            overlap=kde_overlap,
            output_directory=args.output_directory
            / "porous_liquid_distribution_comparison",
        )

        updated_entry = atomlite.PropertyEntry(
            key=solvent.key,
            properties={
                "kde_overlap": kde_overlap,
            },
        )
        solvent_db.update_properties(updated_entry)
    solvent_db.connection.commit()


def _plot_cage_window_distributions(
    cage_db: atomlite.Database,
    output_directory: pathlib.Path,
) -> None:
    for cage in cage_db.get_entries():
        maple.plot.plot_cage_window_distribution(
            cage_name=cage.key,
            cage_apertures=np.array(cage.properties["windows"]),
            output_directory=output_directory,
        )


def _get_combined_cage_apertures(cage_db: atomlite.Database) -> npt.NDArray:
    cage_apertures = []
    for cage in cage_db.get_entries():
        if max(cage.properties["windows"]) < 100:
            cage_apertures.append(np.array(cage.properties["windows"]))
    return np.concatenate(cage_apertures).ravel()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-s",
        "--solvent_database",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "solvents.db",
    )

    parser.add_argument(
        "-c",
        "--cage_database",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "cages.db",
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "6_output",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
