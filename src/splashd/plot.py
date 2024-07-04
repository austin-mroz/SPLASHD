import pathlib

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import seaborn as sns


def plot_steric_parameter_distribution(
    ion_name: str,
    ion_waa_wL: npt.NDArray,
    ion_waa_wB1: npt.NDArray,
    ion_waa_wB5: npt.NDArray,
    output_directory: pathlib.Path,
) -> None:
    pink = "#ff79c6"
    cyan = "#8be9fd"
    dustyblue = "#6277a4"

    plt.axvline(float(np.mean(ion_waa_wL)), 0, 5, color=pink)
    sns.kdeplot(list(ion_waa_wL), fill=True, color=pink, label="L")

    plt.axvline(float(np.mean(ion_waa_wB1)), 0, 5, color=cyan)
    sns.kdeplot(list(ion_waa_wB1), fill=True, color=cyan, label="B1")

    plt.axvline(float(np.mean(ion_waa_wB5)), 0, 5, color=dustyblue)
    sns.kdeplot(list(ion_waa_wB5), fill=True, color=dustyblue, label="B5")

    plt.legend()

    plt.xlabel("size (A)")
    plt.title(f"{ion_name}")
    plt.savefig(output_directory / f"{ion_name}_steric_parameter_distribution.png")
    plt.clf()


def plot_cage_window_distribution(
    cage_name: str,
    cage_apertures: npt.NDArray,
    output_directory: pathlib.Path,
    cage_color: str = "#6277a4",
) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)

    sns.kdeplot(
        list(cage_apertures),
        fill=True,
        color=cage_color,
    )
    plt.xlabel("size (A)")
    plt.title(f"{cage_name}")
    plt.savefig(output_directory / f"{cage_name}_window_distribution.png")
    plt.savefig(output_directory / f"{cage_name}_window_distribution.pdf")

    plt.clf()


def plot_cage_distribution_with_solvent_means(
    method_name: str,
    cage_name: str,
    cage_apertures: npt.NDArray,
    solvent_name: str,
    solvent_sizes: npt.NDArray,
    overlap: float,
    output_directory: pathlib.Path,
    cage_color: str = "#6277a4",
    solvent_color: str = "#ff79c6",
) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)

    plt.axvline(float(np.mean(solvent_sizes)), 0, 5, color=solvent_color)
    sns.kdeplot(
        list(cage_apertures),
        fill=True,
        color=cage_color,
    )
    plt.xlabel("size (A)")
    plt.title(f"{solvent_name}: {overlap}")
    plt.savefig(output_directory / f"{cage_name}_{solvent_name}_{method_name}_plot.png")
    plt.savefig(output_directory / f"{cage_name}_{solvent_name}_{method_name}_plot.pdf")

    plt.clf()


def plot_cage_distribution_with_solvent_distributions(
    method_name: str,
    cage_name: str,
    cage_apertures: npt.NDArray,
    solvent_name: str,
    solvent_sizes: npt.NDArray,
    overlap: float,
    output_directory: pathlib.Path,
    cage_color: str = "#6277a4",
    solvent_color: str = "#ff79c6",
) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)

    plt.axvline(float(np.mean(solvent_sizes)), 0, 5, color=solvent_color)

    sns.kdeplot(
        list(solvent_sizes),
        fill=True,
        color=solvent_color,
        label=f"{solvent_name}",
    )

    sns.kdeplot(
        list(cage_apertures),
        fill=True,
        color=cage_color,
        label=f"{cage_name}",
    )
    plt.legend()
    plt.xlabel("size (A)")
    plt.title(f"{solvent_name}: {overlap}")
    plt.savefig(output_directory / f"{cage_name}_{solvent_name}_{method_name}_plot.png")
    plt.savefig(output_directory / f"{cage_name}_{solvent_name}_{method_name}_plot.pdf")
    plt.clf()
