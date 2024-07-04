import contextlib
import os
import pathlib
import subprocess
import typing
from dataclasses import dataclass

import rdkit.Chem.AllChem as rdkit
import stk
import stko


def run_xtb_single_point_energy(
    molecule: rdkit.Mol,
    output_directory: pathlib.Path | str,
    charge: int = 0,
    level: str = "normal",
    xtb_path: pathlib.Path | str = "xtb",
) -> float:
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    xyz_path = str(output_directory / "geom.xyz")
    rdkit.MolToXYZFile(molecule, xyz_path)

    with current_working_directory(output_directory):
        process = subprocess.run(
            [
                str(xtb_path),
                xyz_path,
                "--chrg",
                f"{charge}",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    with open(output_directory / "xtb.stdout", "w") as f:
        f.write(process.stdout)
    with open(output_directory / "xtb.stderr", "w") as f:
        f.write(process.stderr)

    with open(output_directory / "xtb.stdout", "r") as f:
        lines = f.readlines()

    total_energy_line = [line for line in lines if "TOTAL ENERGY" in line][0]
    return float(total_energy_line.split(" ")[-5])


@contextlib.contextmanager
def current_working_directory(directory: pathlib.Path) -> typing.Any:
    original_directory = os.getcwd()  # aka OGD
    try:
        os.chdir(directory)
        yield
    finally:
        os.chdir(original_directory)


def get_conformer_SP_energy(
    xtb_SP_directory: pathlib.Path,
) -> float:
    with open(xtb_SP_directory / "xtb.stdout", "r") as f:
        lines = f.readlines()

    total_energy_line = [line for line in lines if "TOTAL ENERGY" in line][0]
    return float(total_energy_line.split(" ")[-5])


def optimize_geometry(
    system: stk.ConstructedMolecule | stk.BuildingBlock,
    name: str,
    output_directory: pathlib.Path | str,
) -> pathlib.Path:
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    optimizers = _get_optimizers(output_directory)

    if isinstance(system, stk.ConstructedMolecule):
        system = optimizers.restricted_true.optimize(system)

    restricted_system = optimizers.restricted_false.optimize(system)

    system_md_run = optimizers.md_run.optimize(restricted_system)

    conformer_xyz_directory = output_directory / "conformer_xyzs"
    conformer_xyz_directory.mkdir(parents=True, exist_ok=True)

    for mae_conformer in output_directory.joinpath("MD").glob("*_conf*.mae"):
        conformer_number = str(mae_conformer).removesuffix(".mae").split("_")[-1]
        system_conformer = system_md_run.with_structure_from_file(mae_conformer)
        stk.XyzWriter().write(
            system_conformer,
            str(conformer_xyz_directory / f"{name}_conformer_{conformer_number}_.xyz"),
        )

    return conformer_xyz_directory


@dataclass(frozen=True, slots=True)
class Optimizers:
    restricted_true: stko.MacroModelForceField
    restricted_false: stko.MacroModelForceField
    md_run: stko.MacroModelForceField


def _get_optimizers(output_directory: pathlib.Path) -> "Optimizers":
    macromodel_path = "/opt/schrodinger2023-2/"
    restricted_true = stko.MacroModelForceField(
        macromodel_path,
        output_dir=output_directory / "rT",
        restricted=True,
    )
    restricted_false = stko.MacroModelForceField(
        macromodel_path,
        output_dir=output_directory / "rF",
        restricted=False,
    )
    md_run = stko.MacroModelMD(
        macromodel_path,
        output_dir=output_directory / "MD",
        temperature=300,  # K
        conformers=50,
        time_step=0.7,  # fs
        eq_time=50,  # ps
        simulation_time=10000,  # ps
        maximum_iterations=1,
        minimum_gradient=1.0,
    )
    return Optimizers(
        restricted_true=restricted_true,
        restricted_false=restricted_false,
        md_run=md_run,
    )
