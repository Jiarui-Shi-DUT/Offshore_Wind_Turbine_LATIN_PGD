# -*- coding: utf-8 -*-
"""First nonzero end-to-end activation probe for the tower LATIN-PGD solver."""

from __future__ import annotations

import numpy as np

from examples.elastic_tapered_tower import (
    TowerConfiguration,
    create_tower_geometry,
)
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import create_reversed_top_force_history
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from latin.tower_initialization import (
    compute_tower_elastic_initialization,
)
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
from material.viscoplastic_damage_1d import MaterialParameters


def _fmt(values: np.ndarray) -> str:
    return np.array2string(
        np.asarray(values, dtype=np.float64),
        precision=6,
        suppress_small=False,
        max_line_width=140,
    )


def main() -> None:
    configuration = TowerConfiguration(
        horizontal_force=1.0e6,
        n_elements=10,
        n_gauss=2,
        n_circumferential=16,
        n_radial=1,
    )
    material = MaterialParameters()
    geometry = create_tower_geometry(configuration)
    mesh = create_uniform_vertical_tower_mesh(
        height=configuration.height,
        n_elements=configuration.n_elements,
    )
    system = ViscoplasticDamageTowerSystem2D(
        mesh=mesh,
        tower_geometry=geometry,
        material=material,
        n_gauss=configuration.n_gauss,
        n_circumferential=configuration.n_circumferential,
        n_radial=configuration.n_radial,
    )
    operator = build_tower_equilibrium_operator(system)

    loading = create_reversed_top_force_history(
        force_amplitude=1.0e6,
        period=10.0,
        n_cycles=1,
        increments_per_cycle=40,
    )
    load_vectors = np.stack(
        [
            top_horizontal_load_vector(
                mesh=mesh,
                horizontal_force=float(force),
            )
            for force in loading.forces
        ],
        axis=0,
    )

    initialization = compute_tower_elastic_initialization(
        time=loading.times,
        load_vectors=load_vectors,
        materials=material,
        equilibrium_operator=operator,
        stress_to_force_factor=system.stress_to_force_factor,
    )

    print("=" * 96)
    print("Tower LATIN-PGD reversed activation probe")
    print("=" * 96)
    print(
        "Grid: Nt={}, Nq={} ({} x {} x {})".format(
            initialization.state.n_time,
            initialization.state.n_material_points,
            configuration.n_elements,
            configuration.n_gauss,
            configuration.n_circumferential * configuration.n_radial,
        )
    )
    print(
        "Elastic init max free equilibrium residual = {:.9e} N".format(
            initialization.maximum_free_equilibrium_residual
        )
    )
    print(
        "Elastic init max |stress| = {:.9e} MPa".format(
            float(np.max(np.abs(initialization.state.stress)))
        )
    )

    result = solve_tower_latin_pgd(
        initial_state=initialization.state,
        materials=material,
        metric=operator.metric,
        equilibrium_operator=operator,
        mode_significance_tolerance=0.0,
        acceptance_tolerance=0.0,
        max_iterations=1,
    )

    print("-" * 96)
    print("One outer LATIN iteration")
    print("-" * 96)
    print("termination_reason =", result.termination_reason.value)
    print("converged          =", result.converged)
    print("failure_reason     =", result.failure_reason)
    print("committed iterations =", result.iterations)
    print("attempted iterations =", result.attempted_iterations)
    print("trial evaluations    =", result.trial_evaluations)
    print("trial kinds          =", result.trial_kind_history)
    print("commit kinds         =", result.commit_kind_history)
    print("final basis modes    =", result.basis.n_modes)
    print("total modes added    =", result.total_modes_added)
    print("accepted indicator   = {:.9e}".format(result.final_indicator))
    print("trial indicators     =", _fmt(result.trial_indicator_history))
    print("saturation history   =", _fmt(result.saturation_history))
    print(
        "reduced residuals    =",
        _fmt(result.trial_reduced_residual_history),
    )

    if result.last_trial_a is not None:
        print(
            "Trial A: xi={:.9e}, zeta={:.9e}".format(
                result.last_trial_a.indicator,
                result.last_trial_a.saturation,
            )
        )

    enrichment = result.last_enrichment_result
    if enrichment is not None:
        print("-" * 96)
        print("Enrichment diagnostics")
        print("-" * 96)
        print("accepted               =", enrichment.accepted)
        print("failure_reason          =", enrichment.failure_reason)
        print(
            "fixed-point converged  =",
            enrichment.fixed_point_converged,
        )
        print(
            "fixed-point iterations =",
            enrichment.fixed_point_iterations,
        )
        print(
            "fixed-point history     =",
            _fmt(enrichment.fixed_point_history),
        )
        print(
            "spatial novelty         = {:.9e}".format(
                enrichment.spatial_novelty
            )
        )
        print(
            "temporal significance   = {:.9e}".format(
                enrichment.temporal_significance
            )
        )
        print(
            "orthogonality error     = {:.9e}".format(
                enrichment.orthogonality_error
            )
        )
        print(
            "residual before/after   = {:.9e} / {:.9e}".format(
                enrichment.residual_norm_before,
                enrichment.residual_norm_after,
            )
        )
        print(
            "residual benefit        = {:.9e}".format(
                enrichment.residual_benefit
            )
        )

    if result.last_trial_b is not None:
        print(
            "Trial B: xi={:.9e}, zeta={:.9e}".format(
                result.last_trial_b.indicator,
                result.last_trial_b.saturation,
            )
        )

    print("-" * 96)
    print("Returned persistent state")
    print("-" * 96)
    print(
        "max |plastic strain| = {:.9e}".format(
            float(np.max(np.abs(result.state.plastic_strain)))
        )
    )
    print(
        "max damage           = {:.9e}".format(
            float(np.max(result.state.damage))
        )
    )
    print(
        "max |stress|         = {:.9e} MPa".format(
            float(np.max(np.abs(result.state.stress)))
        )
    )
    print("=" * 96)


if __name__ == "__main__":
    main()
