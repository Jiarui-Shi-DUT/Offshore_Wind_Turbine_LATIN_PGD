# -*- coding: utf-8 -*-
"""Diagnostic-only A/B probe for the failing fourth tower PGD mode.

Compare:
A) current tower raw Eq. (70)-(72) fixed-point map;
B) the same tower map with the validated 1D enrichment's in-loop spatial
   orthogonalisation applied to the initial seed and every spatial half-step.

No latin/ core file is modified.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration, create_tower_geometry
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import create_reversed_top_force_history
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import ViscoplasticDamageTowerSystem2D
from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import MaterialPointMetric, build_tower_equilibrium_operator
from latin.tower_global_finishing import prepare_frozen_global_data
from latin.tower_initialization import compute_tower_elastic_initialization
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
from latin.tower_local_stage import solve_tower_local_stage
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_search_directions import compute_tower_descent_search_directions
import latin.tower_pgd_enrichment as enrichment_module
from material.viscoplastic_damage_1d import MaterialParameters


MAX_FIXED_POINT_ITERATIONS = 200
FIXED_POINT_TOLERANCE = 1.0e-6
MINIMUM_SPATIAL_NORM = 1.0e-14
RCOND = 1.0e-12


@dataclass(frozen=True)
class ProbeResult:
    label: str
    converged: bool
    history: np.ndarray
    lag3: float
    seed_novelty: float
    min_novelty: float
    final_novelty: float
    max_basis_overlap: float


def _build_problem():
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
            top_horizontal_load_vector(mesh=mesh, horizontal_force=float(force))
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
    return material, operator, initialization


def _rebuild_failing_trial_a(material, operator, initialization):
    baseline = solve_tower_latin_pgd(
        initial_state=initialization.state,
        materials=material,
        metric=operator.metric,
        equilibrium_operator=operator,
        mode_significance_tolerance=0.0,
        acceptance_tolerance=0.0,
        max_iterations=20,
        max_fixed_point_iterations=120,
    )

    materials = (material,) * baseline.state.n_material_points
    local_state = solve_tower_local_stage(
        global_state=baseline.state,
        materials=materials,
    )
    directions = compute_tower_descent_search_directions(
        local_state=local_state,
        materials=materials,
        regularization=0.15,
    )
    frozen_data = prepare_frozen_global_data(
        baseline_state=baseline.state,
        local_state=local_state,
        directions=directions,
        equilibrium_operator=operator,
    )
    fixed_a = update_tower_pgd_time_functions(
        basis=baseline.basis,
        time=baseline.state.time,
        forcing=frozen_data.full_plastic_forcing,
        H_sigma=directions.H_sigma,
        metric=operator.metric,
        equilibrium_operator=operator,
        reduced_tolerance=1.0e-4,
        rcond=RCOND,
    )
    return baseline, directions, fixed_a


def _basis_matrix(basis: PGDBasisTower) -> np.ndarray:
    return basis.spatial_plastic_strain_matrix()


def _basis_overlap(vector, basis, metric: MaterialPointMetric) -> float:
    B = _basis_matrix(basis)
    if B.shape[1] == 0:
        return 0.0

    vnorm = metric.norm(vector)
    if vnorm <= np.finfo(np.float64).eps:
        return 0.0

    values = []
    for j in range(B.shape[1]):
        bj = B[:, j]
        bnorm = metric.norm(bj)
        if bnorm > np.finfo(np.float64).eps:
            values.append(
                abs(metric.inner_product(bj, vector)) / (bnorm * vnorm)
            )
    return float(max(values, default=0.0))


def _orthogonalize_like_1d(
    vector,
    basis: PGDBasisTower,
    metric: MaterialPointMetric,
):
    """Literal tower analogue of the current 1D basis projection."""
    raw = np.asarray(vector, dtype=np.float64).copy()
    raw_norm = metric.norm(raw)
    if not np.isfinite(raw_norm) or raw_norm <= MINIMUM_SPATIAL_NORM:
        raise RuntimeError("raw spatial function is degenerate")

    work = raw.copy()
    B = _basis_matrix(basis)
    if B.shape[1] > 0:
        weights = metric.weights
        gram = B.T @ (weights[:, None] * B)
        rhs = B.T @ (weights * work)
        coeff = np.linalg.lstsq(gram, rhs, rcond=RCOND)[0]
        work -= B @ coeff

    orth_norm = metric.norm(work)
    if not np.isfinite(orth_norm) or orth_norm <= MINIMUM_SPATIAL_NORM:
        raise RuntimeError("in-loop orthogonalisation removed the candidate")

    novelty = float(orth_norm / raw_norm)
    work /= orth_norm
    overlap = _basis_overlap(work, basis, metric)
    return work, novelty, overlap


def _lag3(recent, time, H_sigma, metric):
    if len(recent) != 4:
        return float("nan")
    return float(
        enrichment_module._pair_change(
            recent[0], recent[3], time, H_sigma, metric
        )
    )


def _run_map(
    *,
    label,
    basis,
    time,
    defect,
    H_sigma,
    metric,
    operator,
    iteration_added,
    orthogonalize_in_loop,
):
    raw_seed = enrichment_module._seed(
        defect, H_sigma, metric, MINIMUM_SPATIAL_NORM
    )

    if orthogonalize_in_loop:
        p, seed_novelty, overlap = _orthogonalize_like_1d(
            raw_seed, basis, metric
        )
    else:
        p = raw_seed
        seed_novelty = 1.0
        overlap = _basis_overlap(p, basis, metric)

    novelties = [float(seed_novelty)]
    max_overlap = float(overlap)

    s = operator.apply_spatial(p).stress
    lam, rate = enrichment_module._temporal_solve(
        p, s, time, defect, H_sigma, metric
    )
    current = PGDModeTower(p, s, lam, rate, iteration_added)

    history = []
    recent = [current]
    converged = False

    for _ in range(MAX_FIXED_POINT_ITERATIONS):
        p_raw, s_raw = enrichment_module._spatial_solve(
            current.temporal_amplitude,
            current.temporal_rate,
            time,
            defect,
            H_sigma,
            metric,
            operator,
            MINIMUM_SPATIAL_NORM,
        )

        if orthogonalize_in_loop:
            p, novelty, overlap = _orthogonalize_like_1d(
                p_raw, basis, metric
            )
            s = operator.apply_spatial(p).stress
        else:
            p = p_raw
            s = s_raw
            novelty = 1.0
            overlap = _basis_overlap(p, basis, metric)

        novelties.append(float(novelty))
        max_overlap = max(max_overlap, float(overlap))

        lam, rate = enrichment_module._temporal_solve(
            p, s, time, defect, H_sigma, metric
        )
        candidate = PGDModeTower(p, s, lam, rate, iteration_added)
        chi = enrichment_module._pair_change(
            current, candidate, time, H_sigma, metric
        )
        history.append(float(chi))
        current = candidate

        recent.append(current)
        if len(recent) > 4:
            recent.pop(0)

        if chi <= FIXED_POINT_TOLERANCE:
            converged = True
            break

    return ProbeResult(
        label=label,
        converged=converged,
        history=np.asarray(history, dtype=np.float64),
        lag3=_lag3(recent, time, H_sigma, metric),
        seed_novelty=float(seed_novelty),
        min_novelty=float(min(novelties)),
        final_novelty=float(novelties[-1]),
        max_basis_overlap=float(max_overlap),
    )


def _print_result(result: ProbeResult):
    last = (
        float(result.history[-1])
        if result.history.size
        else float("nan")
    )
    tail = result.history[-min(12, result.history.size):]
    print(result.label)
    print("  converged           =", result.converged)
    print("  iterations          =", result.history.size)
    print("  last chi            = {:.12e}".format(last))
    print("  lag-3 distance      = {:.12e}".format(result.lag3))
    print("  seed novelty        = {:.12e}".format(result.seed_novelty))
    print("  min novelty         = {:.12e}".format(result.min_novelty))
    print("  final novelty       = {:.12e}".format(result.final_novelty))
    print("  max basis overlap   = {:.12e}".format(result.max_basis_overlap))
    print(
        "  chi tail            =",
        " ".join("{:.9f}".format(float(x)) for x in tail),
    )


def main():
    material, operator, initialization = _build_problem()
    baseline, directions, fixed_a = _rebuild_failing_trial_a(
        material, operator, initialization
    )

    common = dict(
        basis=fixed_a.basis,
        time=baseline.state.time,
        defect=fixed_a.mechanical_residual,
        H_sigma=directions.H_sigma,
        metric=operator.metric,
        operator=operator,
        iteration_added=baseline.attempted_iterations,
    )

    print("=" * 96)
    print("Fourth-mode A/B: current tower map vs 1D-style in-loop orthogonalisation")
    print("=" * 96)
    print(
        "baseline: termination={}, committed={}, attempted={}, rank={}, xi={:.9e}".format(
            baseline.termination_reason.value,
            baseline.iterations,
            baseline.attempted_iterations,
            baseline.basis.n_modes,
            baseline.final_indicator,
        )
    )
    print(
        "failing Trial-A relative residual = {:.9e}".format(
            fixed_a.relative_residual
        )
    )
    print("-" * 96)

    current = _run_map(
        label="A) current tower raw map",
        orthogonalize_in_loop=False,
        **common,
    )
    _print_result(current)

    print("-" * 96)

    one_d_style = _run_map(
        label="B) 1D-style in-loop M-orthogonalised map",
        orthogonalize_in_loop=True,
        **common,
    )
    _print_result(one_d_style)

    print("-" * 96)
    if one_d_style.converged:
        print(
            "Result: isolated 1D-style in-loop orthogonalisation restores "
            "ordinary fixed-point convergence."
        )
    else:
        print(
            "Result: isolated 1D-style in-loop orthogonalisation does NOT "
            "restore ordinary fixed-point convergence."
        )
        print(
            "If lag-3 remains small and chi stays in the known three-phase "
            "family, close this hypothesis and next compare the spatial "
            "half-step itself: literal 1D-style direct weighted LS versus "
            "tower Eq. (70)-(71)."
        )
    print(
        "Diagnostic only: no latin/ core or transaction acceptance logic "
        "is changed."
    )
    print("=" * 96)


if __name__ == "__main__":
    main()
