# -*- coding: utf-8 -*-
"""
Coarse nonlinear load-level probe for the NREL 5 MW tower.

This script is a diagnostic calibration step, not the final fatigue
analysis. It monotonically increases the tower-top horizontal force and
reports when the default viscoplastic-damage material first develops
plastic strain or damage.

Use:
    python -m examples.nonlinear_tower_load_probe
    python -m examples.nonlinear_tower_load_probe --elements 20
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from examples.elastic_tapered_tower import (
    TowerConfiguration,
    create_tower_geometry,
)
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import (
    NonlinearTowerConvergenceError,
    ViscoplasticDamageTowerSystem2D,
    ViscoplasticTowerResponse,
    solve_nonlinear_tower_load_step,
)
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class CriticalFiberState:
    """Largest inelastic indicators found in one committed tower state."""

    maximum_absolute_stress: float
    maximum_absolute_plastic_strain: float
    maximum_damage: float
    stress_location: Tuple[int, int, int]
    plastic_location: Tuple[int, int, int]
    damage_location: Tuple[int, int, int]


def parse_force_levels(text: str) -> Tuple[float, ...]:
    """
    Parse comma-separated force levels expressed in MN.

    Example
    -------
    ``0.2,0.4,0.6`` becomes forces in newtons.
    """
    if not isinstance(text, str):
        raise TypeError("Force levels must be provided as text.")

    values = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            raise ValueError(
                "Force-level entries must not be empty."
            )
        try:
            force_mn = float(item)
        except ValueError as error:
            raise ValueError(
                "Every force level must be a real number in MN."
            ) from error
        if not np.isfinite(force_mn) or force_mn <= 0.0:
            raise ValueError(
                "Every force level must be finite and positive."
            )
        values.append(force_mn * 1.0e6)

    if not values:
        raise ValueError("At least one force level is required.")
    if np.any(np.diff(np.asarray(values)) <= 0.0):
        raise ValueError(
            "Force levels must be strictly increasing."
        )

    return tuple(float(value) for value in values)


def extract_critical_fiber_state(
    response: ViscoplasticTowerResponse,
) -> CriticalFiberState:
    """Search all elements, Gauss points, and fibers."""
    maximum_absolute_stress = -1.0
    maximum_absolute_plastic_strain = -1.0
    maximum_damage = -1.0

    stress_location = (0, 0, 0)
    plastic_location = (0, 0, 0)
    damage_location = (0, 0, 0)

    for element_index, element_response in enumerate(
        response.element_responses
    ):
        for gauss_index, section_response in enumerate(
            element_response.section_responses
        ):
            stresses = np.abs(section_response.fiber_stresses)
            plastic_strains = np.abs(
                section_response.plastic_strains
            )
            damages = section_response.damages

            stress_fiber = int(np.argmax(stresses))
            plastic_fiber = int(np.argmax(plastic_strains))
            damage_fiber = int(np.argmax(damages))

            stress_value = float(stresses[stress_fiber])
            plastic_value = float(
                plastic_strains[plastic_fiber]
            )
            damage_value = float(damages[damage_fiber])

            if stress_value > maximum_absolute_stress:
                maximum_absolute_stress = stress_value
                stress_location = (
                    element_index,
                    gauss_index,
                    stress_fiber,
                )
            if plastic_value > maximum_absolute_plastic_strain:
                maximum_absolute_plastic_strain = plastic_value
                plastic_location = (
                    element_index,
                    gauss_index,
                    plastic_fiber,
                )
            if damage_value > maximum_damage:
                maximum_damage = damage_value
                damage_location = (
                    element_index,
                    gauss_index,
                    damage_fiber,
                )

    return CriticalFiberState(
        maximum_absolute_stress=maximum_absolute_stress,
        maximum_absolute_plastic_strain=(
            maximum_absolute_plastic_strain
        ),
        maximum_damage=maximum_damage,
        stress_location=stress_location,
        plastic_location=plastic_location,
        damage_location=damage_location,
    )


def create_probe_system(
    configuration: TowerConfiguration,
    material: MaterialParameters,
) -> ViscoplasticDamageTowerSystem2D:
    """Build a reduced-cost nonlinear tower for load calibration."""
    mesh = create_uniform_vertical_tower_mesh(
        height=configuration.height,
        n_elements=configuration.n_elements,
    )
    geometry = create_tower_geometry(configuration)

    return ViscoplasticDamageTowerSystem2D(
        mesh=mesh,
        tower_geometry=geometry,
        material=material,
        n_gauss=configuration.n_gauss,
        n_circumferential=(
            configuration.n_circumferential
        ),
        n_radial=configuration.n_radial,
    )


def run_probe(
    force_levels: Sequence[float],
    n_elements: int,
    step_duration: float,
) -> None:
    """Run sequential monotonic force levels and print diagnostics."""
    if n_elements < 1:
        raise ValueError("n_elements must be at least 1.")
    if not np.isfinite(step_duration) or step_duration <= 0.0:
        raise ValueError(
            "step_duration must be finite and positive."
        )

    configuration = TowerConfiguration(
        n_elements=n_elements,
        n_gauss=2,
        n_circumferential=16,
        n_radial=1,
    )
    material = MaterialParameters()
    system = create_probe_system(
        configuration=configuration,
        material=material,
    )

    print("=" * 104)
    print("Coarse nonlinear load-level probe: NREL 5 MW tower")
    print("=" * 104)
    print(
        "Material defaults: "
        f"E={material.E:.6g} MPa, "
        f"sigma_y={material.sigma_y:.6g} MPa, "
        f"k_damage={material.k_damage:.6g}"
    )
    print(
        "Discretisation: "
        f"{configuration.n_elements} elements, "
        f"{configuration.n_gauss} Gauss points/element, "
        f"{configuration.n_circumferential} x "
        f"{configuration.n_radial} fibers/section"
    )
    print(
        "Each row is a committed state from a monotonic force ramp; "
        "results are rate-dependent."
    )
    print("-" * 104)
    print(
        " force(MN)   top ux(m)   max|stress|(MPa)   "
        "max|eps_p|       max D       Newton"
    )

    first_plastic_force = None
    first_damage_force = None

    for step_index, force in enumerate(force_levels, start=1):
        load_vector = top_horizontal_load_vector(
            mesh=system.mesh,
            horizontal_force=float(force),
        )
        time = step_index * step_duration

        try:
            solution = solve_nonlinear_tower_load_step(
                system=system,
                time=time,
                load_vector=load_vector,
                max_iterations=40,
                relative_residual_tolerance=1.0e-8,
                absolute_residual_tolerance=1.0e-4,
            )
        except NonlinearTowerConvergenceError as error:
            print(
                f"{force / 1.0e6:10.4f}   "
                "NON-CONVERGED: "
                + str(error)
            )
            break

        critical = extract_critical_fiber_state(
            solution.response
        )
        top_displacement = float(solution.displacements[-3])

        print(
            f"{force / 1.0e6:10.4f}   "
            f"{top_displacement:9.3e}   "
            f"{critical.maximum_absolute_stress:17.6e}   "
            f"{critical.maximum_absolute_plastic_strain:11.4e}   "
            f"{critical.maximum_damage:9.4e}   "
            f"{solution.iterations:6d}"
        )

        if (
            first_plastic_force is None
            and critical.maximum_absolute_plastic_strain
            > 1.0e-12
        ):
            first_plastic_force = float(force)

        if (
            first_damage_force is None
            and critical.maximum_damage > 1.0e-12
        ):
            first_damage_force = float(force)

    print("-" * 104)
    if first_plastic_force is None:
        print(
            "No plastic strain above 1e-12 was detected "
            "within the tested force range."
        )
    else:
        print(
            "First tested force with plastic strain: "
            f"{first_plastic_force / 1.0e6:.6g} MN"
        )

    if first_damage_force is None:
        print(
            "No damage above 1e-12 was detected "
            "within the tested force range."
        )
    else:
        print(
            "First tested force with damage: "
            f"{first_damage_force / 1.0e6:.6g} MN"
        )
    print("=" * 104)


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Probe the onset of nonlinear material response "
            "in the NREL 5 MW tower."
        )
    )
    parser.add_argument(
        "--levels",
        type=str,
        default="0.2,0.4,0.6,0.8,1.0,1.2",
        help=(
            "Strictly increasing tower-top force levels in MN, "
            "separated by commas."
        ),
    )
    parser.add_argument(
        "--elements",
        type=int,
        default=10,
        help=(
            "Number of tower elements for this coarse diagnostic."
        ),
    )
    parser.add_argument(
        "--step-duration",
        type=float,
        default=1.0,
        help=(
            "Analysis-time increment between successive force levels."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the command-line diagnostic."""
    arguments = parse_arguments()
    run_probe(
        force_levels=parse_force_levels(arguments.levels),
        n_elements=arguments.elements,
        step_duration=arguments.step_duration,
    )


if __name__ == "__main__":
    main()
