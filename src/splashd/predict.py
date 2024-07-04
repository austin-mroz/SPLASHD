import pathlib

import numpy as np
import numpy.typing as npt
from scipy.stats import gaussian_kde


def calculate_kde_overlap(
    cage_apertures: npt.NDArray,
    solvent_sizes: npt.NDArray,
) -> float:
    cage_kde = gaussian_kde(cage_apertures)
    solvent_kde = gaussian_kde(solvent_sizes)

    x_min = min(solvent_sizes.min(), cage_apertures.min())
    x_max = max(solvent_sizes.max(), cage_apertures.max())

    dx = 0.2 * (x_max - x_min)
    x_min -= dx
    x_max += dx

    x = np.linspace(x_min, x_max, 500)
    x_solvent_kde = solvent_kde(x)
    x_cage_kde = cage_kde(x)
    intersection_x = np.minimum(x_solvent_kde, x_cage_kde)

    return np.trapz(intersection_x, x)
