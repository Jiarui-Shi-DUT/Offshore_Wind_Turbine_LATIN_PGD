# -*- coding: utf-8 -*-
"""Diagnostic A/B/C comparison for the failing fourth tower PGD mode.

This script does not modify latin/ core code.

It reproduces the known failing fourth-mode Trial-A problem and compares:

A) current tower raw fixed-point map
   - tower Eq. (70)-(71) spatial half-step
   - no in-loop basis orthogonalisation

B) tower Eq. (70)-(71) + validated-1D-style in-loop M-orthogonalisation

C) literal 1D-style direct weighted-least-squares spatial half-step
   + the same in-loop M-orthogonalisation as case B

Cases B and C differ only in the spatial half-step, so B -> C isolates whether
the tower Eq. (70)-(71) structural spatial solve, rather than basis
orthogonality, is responsible for the persistent period-3 orbit.

The temporal half-step remains the current backward-Euler Eq. (72) solve in all
three cases.  The current tower complete-pair convergence criterion and
tolerance are also retained in all three cases.

Additional diagnostics in this version:

1) free-DOF equilibrium of the final spatial stress mode

       H^T M s ~= 0,

2) immediate backward-Euler mechanical-residual benefit of the final raw
   separated mode.
"""

from __future__ import annotations

from dataclasses import dataclass

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
from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
    build_tower_equilibrium_operator,
)
from latin.tower_global_finishing import prepare_frozen_global_data
from latin.tower_initialization import (
    compute_tower_elastic_initialization,
)
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
from latin.tower_local_stage import solve_tower_local_stage
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_search_directions import (
    compute_tower_descent_search_directions,
)
import latin.tower_pgd_enrichment as enrichment_module
from material.viscoplastic_damage_1d import MaterialParameters


MAX_FIXED_POINT_ITERATIONS = 200
FIXED_POINT_TOLERANCE = 1.0e-6
MINIMUM_SPATIAL_NORM = 1.0e-14
RCOND = 1.0e-12


@dataclass(frozen=True)
class MapResult:
    label: str
    converged: bool
    iterations: int
    history: np.ndarray
    lag3: float
    seed_novelty: float
    minimum_novelty: float
    final_novelty: float
    maximum_basis_overlap: float
    final_spatial_norm: float
    equilibrium_absolute_norm: float
    equilibrium_relative_norm: float
    residual_norm_before: float
    residual_norm_after_raw_mode: float
    raw_mode_residual_benefit: float


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


def _basis_overlap(
    vector: np.ndarray,
    basis: PGDBasisTower,
    metric: MaterialPointMetric,
) -> float:
    B = _basis_matrix(basis)
    if B.shape[1] == 0:
        return 0.0

    v_norm = metric.norm(vector)
    if v_norm <= np.finfo(np.float64).eps:
        return 0.0

    values = []
    for j in range(B.shape[1]):
        bj = B[:, j]
        b_norm = metric.norm(bj)
        if b_norm <= np.finfo(np.float64).eps:
            continue
        values.append(
            abs(metric.inner_product(bj, vector)) / (b_norm * v_norm)
        )
    return float(max(values, default=0.0))


def _orthogonalize_and_normalise_like_1d(
    spatial_function: np.ndarray,
    basis: PGDBasisTower,
    metric: MaterialPointMetric,
) -> tuple[np.ndarray, float, float]:
    raw = np.asarray(spatial_function, dtype=np.float64).copy()
    raw_norm = metric.norm(raw)
    if (
        not np.isfinite(raw_norm)
        or raw_norm <= MINIMUM_SPATIAL_NORM
    ):
        raise RuntimeError("raw spatial function is degenerate")

    work = raw.copy()
    B = _basis_matrix(basis)

    if B.shape[1] > 0:
        weights = metric.weights
        weighted_basis = weights[:, None] * B
        gram = B.T @ weighted_basis
        rhs = B.T @ (weights * work)
        coeff = np.linalg.lstsq(
            gram,
            rhs,
            rcond=RCOND,
        )[0]
        work -= B @ coeff

    orth_norm = metric.norm(work)
    if (
        not np.isfinite(orth_norm)
        or orth_norm <= MINIMUM_SPATIAL_NORM
    ):
        raise RuntimeError(
            "1D-style orthogonalisation removed the candidate"
        )

    novelty = float(orth_norm / raw_norm)
    work /= orth_norm
    overlap = _basis_overlap(work, basis, metric)
    return work, novelty, overlap


def _trapezoidal_weights(time: np.ndarray) -> np.ndarray:
    """Same nodal quadrature used by latin/pgd_enrichment.py."""
    dt = np.diff(time)
    weights = np.zeros(time.size, dtype=np.float64)
    weights[0] = 0.5 * dt[0]
    weights[-1] = 0.5 * dt[-1]
    if time.size > 2:
        weights[1:-1] = 0.5 * (dt[:-1] + dt[1:])
    return weights


def _dense_stress_matrix(
    operator: TowerEquilibriumOperator,
) -> np.ndarray:
    """Dense linear map A such that s = A p.

    This is only a coarse-benchmark diagnostic.  The production tower operator
    remains matrix-free with respect to this material-point mapping.
    """
    nq = operator.n_material_points
    matrix = np.empty((nq, nq), dtype=np.float64)
    eye = np.eye(nq, dtype=np.float64)
    for j in range(nq):
        matrix[:, j] = operator.apply_spatial(eye[:, j]).stress
    if np.any(~np.isfinite(matrix)):
        raise RuntimeError("dense stress matrix contains non-finite values")
    return matrix


def _direct_1d_style_spatial_solve(
    temporal_amplitude: np.ndarray,
    temporal_rate: np.ndarray,
    defect: np.ndarray,
    time: np.ndarray,
    H_sigma: np.ndarray,
    metric: MaterialPointMetric,
    stress_matrix: np.ndarray,
) -> np.ndarray:
    """Literal tower analogue of 1D _solve_spatial_function().

    For each time node n,

        B_n = rate_n I - lambda_n diag(H_sigma_n) A

    and the spatial vector minimises

        sum_n w_t,n ||B_n p + defect_n||^2_(M H_sigma_n^-1).

    This intentionally follows the validated 1D code path rather than the
    tower Eq. (70)-(71) structural elimination.
    """
    nq = metric.n_material_points
    if stress_matrix.shape != (nq, nq):
        raise ValueError("stress_matrix has incompatible shape")

    time_weights = _trapezoidal_weights(time)
    normal = np.zeros((nq, nq), dtype=np.float64)
    rhs = np.zeros(nq, dtype=np.float64)
    identity = np.eye(nq, dtype=np.float64)

    for n in range(time.size):
        Hn = H_sigma[n, :]
        Bn = (
            temporal_rate[n] * identity
            - temporal_amplitude[n]
            * Hn[:, None]
            * stress_matrix
        )
        weights = (
            time_weights[n]
            * metric.weights
            / Hn
        )
        weighted_B = weights[:, None] * Bn
        normal += Bn.T @ weighted_B
        rhs -= Bn.T @ (weights * defect[n, :])

    p = np.linalg.lstsq(
        normal,
        rhs,
        rcond=RCOND,
    )[0]

    if np.any(~np.isfinite(p)):
        raise RuntimeError(
            "direct 1D-style spatial solve returned non-finite values"
        )
    if metric.norm(p) <= MINIMUM_SPATIAL_NORM:
        raise RuntimeError(
            "direct 1D-style spatial solve returned a degenerate vector"
        )
    return p


def _lag3(recent_modes, time, H_sigma, metric) -> float:
    if len(recent_modes) != 4:
        return float("nan")
    return float(
        enrichment_module._pair_change(
            recent_modes[0],
            recent_modes[3],
            time,
            H_sigma,
            metric,
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
    spatial_kind,
    use_inloop_orthogonalisation,
    stress_matrix=None,
) -> MapResult:
    raw_seed = enrichment_module._seed(
        defect,
        H_sigma,
        metric,
        MINIMUM_SPATIAL_NORM,
    )

    if use_inloop_orthogonalisation:
        p, seed_novelty, seed_overlap = (
            _orthogonalize_and_normalise_like_1d(
                raw_seed,
                basis,
                metric,
            )
        )
    else:
        p = raw_seed
        seed_novelty = 1.0
        seed_overlap = _basis_overlap(p, basis, metric)

    maximum_overlap = float(seed_overlap)
    novelties = [float(seed_novelty)]

    s = operator.apply_spatial(p).stress
    lam, rate = enrichment_module._temporal_solve(
        p,
        s,
        time,
        defect,
        H_sigma,
        metric,
    )
    current = PGDModeTower(
        p,
        s,
        lam,
        rate,
        iteration_added,
    )

    history = []
    recent_modes = [current]
    converged = False

    for _ in range(MAX_FIXED_POINT_ITERATIONS):
        if spatial_kind == "tower":
            p_raw, _ = enrichment_module._spatial_solve(
                current.temporal_amplitude,
                current.temporal_rate,
                time,
                defect,
                H_sigma,
                metric,
                operator,
                MINIMUM_SPATIAL_NORM,
            )
        elif spatial_kind == "direct_ls":
            if stress_matrix is None:
                raise RuntimeError(
                    "direct_ls requires the dense stress matrix"
                )
            p_raw = _direct_1d_style_spatial_solve(
                temporal_amplitude=current.temporal_amplitude,
                temporal_rate=current.temporal_rate,
                defect=defect,
                time=time,
                H_sigma=H_sigma,
                metric=metric,
                stress_matrix=stress_matrix,
            )
        else:
            raise ValueError("unknown spatial_kind")

        if use_inloop_orthogonalisation:
            p, novelty, overlap = (
                _orthogonalize_and_normalise_like_1d(
                    p_raw,
                    basis,
                    metric,
                )
            )
        else:
            raw_norm = metric.norm(p_raw)
            if raw_norm <= MINIMUM_SPATIAL_NORM:
                raise RuntimeError("spatial vector became degenerate")
            p = p_raw / raw_norm
            novelty = 1.0
            overlap = _basis_overlap(p, basis, metric)

        novelties.append(float(novelty))
        maximum_overlap = max(maximum_overlap, float(overlap))

        s = operator.apply_spatial(p).stress
        lam, rate = enrichment_module._temporal_solve(
            p,
            s,
            time,
            defect,
            H_sigma,
            metric,
        )
        candidate = PGDModeTower(
            p,
            s,
            lam,
            rate,
            iteration_added,
        )

        chi = enrichment_module._pair_change(
            current,
            candidate,
            time,
            H_sigma,
            metric,
        )
        history.append(float(chi))
        current = candidate

        recent_modes.append(current)
        if len(recent_modes) > 4:
            recent_modes.pop(0)

        if chi <= FIXED_POINT_TOLERANCE:
            converged = True
            break

    # ---------------------------------------------------------------
    # Final-mode equilibrium audit.
    #
    # For the separated stress mode s(q), homogeneous correction
    # equilibrium on free structural DOFs requires
    #
    #     H^T M s = 0.
    # ---------------------------------------------------------------
    equilibrium_residual = operator.equilibrium_residual(
        current.spatial_stress
    )
    equilibrium_absolute_norm = float(
        np.linalg.norm(equilibrium_residual)
    )

    # A cancellation-free internal-force scale for a dimensionless
    # relative equilibrium residual.
    equilibrium_scale_vector = (
        np.abs(operator.compatibility_matrix).T
        @ (
            metric.weights
            * np.abs(current.spatial_stress)
        )
    )
    equilibrium_scale = float(
        np.linalg.norm(equilibrium_scale_vector)
    )
    equilibrium_relative_norm = (
        equilibrium_absolute_norm
        / max(
            equilibrium_scale,
            np.finfo(np.float64).eps,
        )
    )

    # ---------------------------------------------------------------
    # Immediate residual efficacy of the converged raw separated mode.
    #
    # Existing Trial-A residual:
    #
    #     R_A
    #
    # Added-mode contribution:
    #
    #     p * lambda_dot - H_sigma * s * lambda.
    #
    # Evaluate the before/after fields using the same backward-Euler
    # residual norm used by tower_pgd_enrichment.py.
    # ---------------------------------------------------------------
    residual_norm_before = enrichment_module._be_residual_norm(
        defect,
        time,
        H_sigma,
        metric,
    )

    residual_after_raw_mode_field = (
        defect
        + (
            current.temporal_rate[:, None]
            * current.spatial_plastic_strain[None, :]
        )
        - (
            H_sigma
            * current.temporal_amplitude[:, None]
            * current.spatial_stress[None, :]
        )
    )

    residual_norm_after_raw_mode = (
        enrichment_module._be_residual_norm(
            residual_after_raw_mode_field,
            time,
            H_sigma,
            metric,
        )
    )

    if residual_norm_before <= np.finfo(np.float64).eps:
        raw_mode_residual_benefit = 0.0
    else:
        raw_mode_residual_benefit = (
            1.0
            - residual_norm_after_raw_mode
            / residual_norm_before
        )

    return MapResult(
        label=label,
        converged=converged,
        iterations=len(history),
        history=np.asarray(history, dtype=np.float64),
        lag3=_lag3(
            recent_modes,
            time,
            H_sigma,
            metric,
        ),
        seed_novelty=float(seed_novelty),
        minimum_novelty=float(min(novelties)),
        final_novelty=float(novelties[-1]),
        maximum_basis_overlap=float(maximum_overlap),
        final_spatial_norm=float(metric.norm(current.spatial_plastic_strain)),
        equilibrium_absolute_norm=equilibrium_absolute_norm,
        equilibrium_relative_norm=equilibrium_relative_norm,
        residual_norm_before=residual_norm_before,
        residual_norm_after_raw_mode=residual_norm_after_raw_mode,
        raw_mode_residual_benefit=raw_mode_residual_benefit,
    )


def _tail_text(history: np.ndarray, n: int = 12) -> str:
    if history.size == 0:
        return "<empty>"
    tail = history[-min(n, history.size):]
    return " ".join("{:.9f}".format(float(x)) for x in tail)


def _print_result(result: MapResult) -> None:
    last_chi = (
        float(result.history[-1])
        if result.history.size
        else float("nan")
    )
    print(result.label)
    print("  converged           =", result.converged)
    print("  iterations          =", result.iterations)
    print("  last chi            = {:.12e}".format(last_chi))
    print("  lag-3 distance      = {:.12e}".format(result.lag3))
    print("  seed novelty        = {:.12e}".format(result.seed_novelty))
    print("  min novelty         = {:.12e}".format(result.minimum_novelty))
    print("  final novelty       = {:.12e}".format(result.final_novelty))
    print(
        "  max basis overlap   = {:.12e}".format(
            result.maximum_basis_overlap
        )
    )
    print(
        "  final ||p||_M       = {:.12e}".format(
            result.final_spatial_norm
        )
    )
    print(
        "  equilibrium ||r||   = {:.12e}".format(
            result.equilibrium_absolute_norm
        )
    )
    print(
        "  equilibrium relative= {:.12e}".format(
            result.equilibrium_relative_norm
        )
    )
    print(
        "  BE residual before  = {:.12e}".format(
            result.residual_norm_before
        )
    )
    print(
        "  BE residual after   = {:.12e}".format(
            result.residual_norm_after_raw_mode
        )
    )
    print(
        "  raw-mode benefit    = {:.6%}".format(
            result.raw_mode_residual_benefit
        )
    )
    print("  chi tail            =", _tail_text(result.history))


def main() -> None:
    material, operator, initialization = _build_problem()
    baseline, directions, fixed_a = _rebuild_failing_trial_a(
        material,
        operator,
        initialization,
    )

    print("=" * 104)
    print(
        "Fourth-mode A/B/C: tower Eq.(70)-(71) vs literal 1D-style direct weighted LS"
    )
    print("=" * 104)
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
    print(
        "max fixed-point iterations = {}; tolerance = {:.1e}".format(
            MAX_FIXED_POINT_ITERATIONS,
            FIXED_POINT_TOLERANCE,
        )
    )

    print("-" * 104)
    result_a = _run_map(
        label="A) current tower Eq.(70)-(71), no in-loop orthogonalisation",
        basis=fixed_a.basis,
        time=baseline.state.time,
        defect=fixed_a.mechanical_residual,
        H_sigma=directions.H_sigma,
        metric=operator.metric,
        operator=operator,
        iteration_added=baseline.attempted_iterations,
        spatial_kind="tower",
        use_inloop_orthogonalisation=False,
    )
    _print_result(result_a)

    print("-" * 104)
    result_b = _run_map(
        label="B) tower Eq.(70)-(71) + 1D-style in-loop M-orthogonalisation",
        basis=fixed_a.basis,
        time=baseline.state.time,
        defect=fixed_a.mechanical_residual,
        H_sigma=directions.H_sigma,
        metric=operator.metric,
        operator=operator,
        iteration_added=baseline.attempted_iterations,
        spatial_kind="tower",
        use_inloop_orthogonalisation=True,
    )
    _print_result(result_b)

    print("-" * 104)
    print(
        "Building dense material-point stress map for diagnostic direct LS "
        "(coarse Nq={})...".format(operator.n_material_points)
    )
    stress_matrix = _dense_stress_matrix(operator)

    result_c = _run_map(
        label="C) literal 1D-style direct weighted-LS spatial half-step + in-loop M-orthogonalisation",
        basis=fixed_a.basis,
        time=baseline.state.time,
        defect=fixed_a.mechanical_residual,
        H_sigma=directions.H_sigma,
        metric=operator.metric,
        operator=operator,
        iteration_added=baseline.attempted_iterations,
        spatial_kind="direct_ls",
        use_inloop_orthogonalisation=True,
        stress_matrix=stress_matrix,
    )
    _print_result(result_c)

    print("-" * 104)
    print("Interpretation:")
    if result_c.converged and not result_b.converged:
        print(
            "  Case C converged while B did not: the spatial-half-step "
            "formulation is the first isolated difference that removes the "
            "fourth-mode fixed-point failure."
        )
        print(
            "  The new equilibrium and BE-residual diagnostics should now be "
            "used to decide whether this converged mode is structurally "
            "admissible and mechanically useful."
        )
    elif not result_c.converged:
        print(
            "  Case C is still nonconvergent: literal 1D-style spatial LS "
            "does not by itself restore an ordinary fixed point."
        )
        print(
            "  If its lag-3 distance is also small, the period-3 pathology "
            "survives both in-loop orthogonalisation and the 1D-style spatial "
            "half-step, so the next difference to inspect is the fixed-point "
            "convergence/acceptance semantics rather than the spatial solve."
        )
    else:
        print(
            "  Inspect A/B/C tails, equilibrium residuals, and BE residual "
            "benefits carefully before drawing a causal conclusion."
        )
    print(
        "  Diagnostic only: latin/ core and transaction acceptance semantics "
        "remain unchanged."
    )
    print("=" * 104)


if __name__ == "__main__":
    main()
