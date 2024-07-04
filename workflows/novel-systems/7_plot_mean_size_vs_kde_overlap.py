import argparse
import pathlib

import atomlite
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    # load the cage db and only extract the cage apertures
    # load cage db and only extract cage apertures
    cage_apertures = np.array(
        atomlite.Database("cages.db").get_property("cc3", "$.windows")
    )

    solvent_db = atomlite.Database("solvents.db")
    dustyblue = "#6277a4"
    plt.axhspan(0.25, 0.35, color=dustyblue, alpha=0.5, lw=0)
    for solvent_entry, color in zip(
        solvent_db.get_entries(
            keys=[
                "1-Phenylethanol",
                "2-Ethylphenol",
                "2-Methyl-1-phenyl-2-propanol",
                "2-Phenyl-1-propanol",
                "diphenyl_ether",
                "2-Methylphenol",
            ]
        ),
        sns.color_palette("husl", 11),
    ):
        print(solvent_entry.key)
        if solvent_entry.properties["aa_method"] == "core-excluded":
            solvent_sizes = np.array(list(solvent_entry.properties["waa_wL"])).T
        elif solvent_entry.properties["aa_method"] == "core_included":
            solvent_sizes = np.array(list(solvent_entry.properties["waa_wB5"])).T

        plt.scatter(
            x=float(np.mean(solvent_sizes)),
            y=float(solvent_entry.properties["kde_overlap"]),
            color=color,
            edgecolor="black",
            label=f"{solvent_entry.key}_{solvent_entry.properties['kde_overlap']}",
        )

    plt.xlim(3, 6)
    plt.ylim(0, 1)
    plt.legend(loc="lower left", bbox_to_anchor=(0, -1.25))
    plt.xlabel("size (A)")
    plt.title("HSPiP solvents")
    plt.savefig(args.output_directory / f"HSPiP_dist_plot.pdf", bbox_inches="tight")
    plt.savefig(args.output_directory / f"HSPiP_dist_plot.png", bbox_inches="tight")
    plt.clf()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Analyze the cage apertures."),
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "9_output",
    )

    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
