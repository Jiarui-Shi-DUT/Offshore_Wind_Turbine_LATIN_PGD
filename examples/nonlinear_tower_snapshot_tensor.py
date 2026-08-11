# -*- coding: utf-8 -*-
"""
Cycle-phase tensorization of nonlinear tower full-field histories.

The nonlinear cyclic tower solvers already store complete fiber histories at
every periodic time point:

    fiber_strains : (n_points, n_elements, n_gauss, n_fibers)
    fiber_stresses: (n_points, n_elements, n_gauss, n_fibers)
    fiber_states  : (n_points, n_elements, n_gauss, n_fibers, 4)

For LATIN-PGD slow-fast analysis, this module reorganizes the single time axis
into

    slow cycle n x fast phase tau.

Each complete cycle contains ``increments_per_cycle + 1`` phase points,
including both tau = 0 and tau = T. Therefore, the shared boundary state
between adjacent cycles is intentionally duplicated. This makes every row a
complete closed-interval cycle and preserves exact cycle-start/cycle-end
comparisons needed for ratcheting and damage diagnostics.

The module also provides lossless compressed NPZ persistence for
``TowerCyclePhaseSnapshots``. The purpose is to solve an expensive FOM only
once, freeze the resulting cycle-phase snapshots, and reuse exactly the same
reference data for repeated offline low-rank, stage-wise, CP/PGD, and future
LATIN-PGD diagnostics.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

import numpy as np
from numpy.typing import NDArray

from examples.nonlinear_tower_reversed_response import (
    NonlinearCyclicResponse,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Snapshot file-format definition
# ---------------------------------------------------------------------------

SNAPSHOT_FORMAT_VERSION = 1

_SNAPSHOT_ARRAY_KEYS = (
    "cycle_numbers",
    "phase_times",
    "phase_fractions",
    "phase_forces",
    "analysis_times",
    "nodal_displacements",
    "fiber_strains",
    "fiber_stresses",
    "fiber_states",
)


def _cycle_phase_indices(
    n_cycles: int,
    increments_per_cycle: int,
) -> IntArray:
    """Return time-history indices with shape (n_cycles, n_phase_points)."""
    if isinstance(n_cycles, (bool, np.bool_)):
        raise TypeError("n_cycles must be an integer.")
    if isinstance(increments_per_cycle, (bool, np.bool_)):
        raise TypeError("increments_per_cycle must be an integer.")
    if not isinstance(n_cycles, (int, np.integer)):
        raise TypeError("n_cycles must be an integer.")
    if not isinstance(increments_per_cycle, (int, np.integer)):
        raise TypeError("increments_per_cycle must be an integer.")

    n_cycles = int(n_cycles)
    increments_per_cycle = int(increments_per_cycle)

    if n_cycles < 1:
        raise ValueError("n_cycles must be at least 1.")
    if increments_per_cycle < 1:
        raise ValueError(
            "increments_per_cycle must be at least 1."
        )

    cycle_starts = (
        np.arange(n_cycles, dtype=np.int64)
        * increments_per_cycle
    )
    phase_offsets = np.arange(
        increments_per_cycle + 1,
        dtype=np.int64,
    )

    return (
        cycle_starts[:, np.newaxis]
        + phase_offsets[np.newaxis, :]
    )


def _tensorize_history(
    values: FloatArray,
    indices: IntArray,
) -> FloatArray:
    """Gather one periodic history into slow-cycle x fast-phase form."""
    values = np.asarray(values, dtype=np.float64)

    if values.ndim < 1:
        raise ValueError("values must contain a time axis.")
    if indices.ndim != 2:
        raise ValueError("indices must be a two-dimensional array.")
    if int(np.max(indices)) >= values.shape[0]:
        raise ValueError(
            "cycle-phase indices exceed the available history."
        )

    return np.asarray(values[indices], dtype=np.float64)


@dataclass(frozen=True)
class TowerCyclePhaseSnapshots:
    """
    Full-field tower snapshots in slow-cycle x fast-phase coordinates.

    Array dimensions
    ----------------
    analysis_times
        (n_cycles, n_phase_points)
    nodal_displacements
        (n_cycles, n_phase_points, n_dof)
    fiber_strains
        (n_cycles, n_phase_points, n_elements, n_gauss, n_fibers)
    fiber_stresses
        Same shape as ``fiber_strains``.
    fiber_states
        ``fiber_strains.shape + (4,)``.

    Material-state channel order
    ----------------------------
    0: plastic strain
    1: alpha
    2: r_bar
    3: damage
    """

    cycle_numbers: IntArray
    phase_times: FloatArray
    phase_fractions: FloatArray
    phase_forces: FloatArray
    analysis_times: FloatArray
    nodal_displacements: FloatArray
    fiber_strains: FloatArray
    fiber_stresses: FloatArray
    fiber_states: FloatArray

    def __post_init__(self) -> None:
        cycle_numbers = np.asarray(
            self.cycle_numbers,
            dtype=np.int64,
        )
        phase_times = np.asarray(
            self.phase_times,
            dtype=np.float64,
        )
        phase_fractions = np.asarray(
            self.phase_fractions,
            dtype=np.float64,
        )
        phase_forces = np.asarray(
            self.phase_forces,
            dtype=np.float64,
        )
        analysis_times = np.asarray(
            self.analysis_times,
            dtype=np.float64,
        )
        nodal_displacements = np.asarray(
            self.nodal_displacements,
            dtype=np.float64,
        )
        fiber_strains = np.asarray(
            self.fiber_strains,
            dtype=np.float64,
        )
        fiber_stresses = np.asarray(
            self.fiber_stresses,
            dtype=np.float64,
        )
        fiber_states = np.asarray(
            self.fiber_states,
            dtype=np.float64,
        )

        if cycle_numbers.ndim != 1:
            raise ValueError(
                "cycle_numbers must be one-dimensional."
            )
        if cycle_numbers.size < 1:
            raise ValueError(
                "At least one cycle is required."
            )
        if np.any(cycle_numbers < 1):
            raise ValueError(
                "cycle_numbers must contain positive global cycle numbers."
            )
        if (
            cycle_numbers.size > 1
            and np.any(np.diff(cycle_numbers) != 1)
        ):
            raise ValueError(
                "cycle_numbers must be consecutive and strictly increasing."
            )

        if phase_times.ndim != 1:
            raise ValueError(
                "phase_times must be one-dimensional."
            )
        if phase_times.size < 2:
            raise ValueError(
                "phase_times must contain at least two points."
            )
        if np.any(np.diff(phase_times) <= 0.0):
            raise ValueError(
                "phase_times must be strictly increasing."
            )
        if not np.isclose(
            float(phase_times[0]),
            0.0,
        ):
            raise ValueError(
                "phase_times must begin at zero."
            )

        n_cycles = cycle_numbers.size
        n_phase_points = phase_times.size

        for array, name in (
            (phase_fractions, "phase_fractions"),
            (phase_forces, "phase_forces"),
        ):
            if array.shape != (n_phase_points,):
                raise ValueError(
                    name + " must match phase_times."
                )

        if analysis_times.shape != (
            n_cycles,
            n_phase_points,
        ):
            raise ValueError(
                "analysis_times must have shape "
                "(n_cycles, n_phase_points)."
            )

        if nodal_displacements.ndim != 3:
            raise ValueError(
                "nodal_displacements must have shape "
                "(n_cycles, n_phase_points, n_dof)."
            )
        if nodal_displacements.shape[:2] != (
            n_cycles,
            n_phase_points,
        ):
            raise ValueError(
                "nodal_displacements must match "
                "the cycle-phase grid."
            )

        if fiber_strains.ndim != 5:
            raise ValueError(
                "fiber_strains must have shape "
                "(n_cycles, n_phase_points, n_elements, "
                "n_gauss, n_fibers)."
            )
        if fiber_strains.shape[:2] != (
            n_cycles,
            n_phase_points,
        ):
            raise ValueError(
                "fiber_strains must match "
                "the cycle-phase grid."
            )
        if fiber_stresses.shape != fiber_strains.shape:
            raise ValueError(
                "fiber_stresses must match fiber_strains."
            )
        if fiber_states.shape != fiber_strains.shape + (4,):
            raise ValueError(
                "fiber_states must append "
                "four material variables."
            )

        finite_arrays = (
            phase_times,
            phase_fractions,
            phase_forces,
            analysis_times,
            nodal_displacements,
            fiber_strains,
            fiber_stresses,
            fiber_states,
        )

        for array in finite_arrays:
            if np.any(~np.isfinite(array)):
                raise ValueError(
                    "All snapshot arrays must be finite."
                )

        if not np.isclose(
            float(phase_fractions[0]),
            0.0,
        ):
            raise ValueError(
                "phase_fractions must begin at zero."
            )

        if not np.isclose(
            float(phase_fractions[-1]),
            1.0,
        ):
            raise ValueError(
                "phase_fractions must end at one."
            )

        if np.any(
            np.diff(phase_fractions) <= 0.0
        ):
            raise ValueError(
                "phase_fractions must be strictly increasing."
            )

        if np.any(
            fiber_states[..., 3] < 0.0
        ):
            raise ValueError(
                "Fiber damage must be non-negative."
            )

        if np.any(
            fiber_states[..., 3] >= 1.0
        ):
            raise ValueError(
                "Fiber damage must be smaller than one."
            )

        object.__setattr__(
            self,
            "cycle_numbers",
            cycle_numbers.copy(),
        )
        object.__setattr__(
            self,
            "phase_times",
            phase_times.copy(),
        )
        object.__setattr__(
            self,
            "phase_fractions",
            phase_fractions.copy(),
        )
        object.__setattr__(
            self,
            "phase_forces",
            phase_forces.copy(),
        )
        object.__setattr__(
            self,
            "analysis_times",
            analysis_times.copy(),
        )
        object.__setattr__(
            self,
            "nodal_displacements",
            nodal_displacements.copy(),
        )
        object.__setattr__(
            self,
            "fiber_strains",
            fiber_strains.copy(),
        )
        object.__setattr__(
            self,
            "fiber_stresses",
            fiber_stresses.copy(),
        )
        object.__setattr__(
            self,
            "fiber_states",
            fiber_states.copy(),
        )

    @property
    def n_cycles(self) -> int:
        """Return the number of stored complete cycles."""
        return int(self.cycle_numbers.size)

    @property
    def n_phase_points(self) -> int:
        """Return the number of fast-phase points per complete cycle."""
        return int(self.phase_times.size)

    @property
    def fiber_plastic_strains(self) -> FloatArray:
        """Return plastic-strain snapshots."""
        return self.fiber_states[..., 0].copy()

    @property
    def fiber_alphas(self) -> FloatArray:
        """Return kinematic internal-variable snapshots."""
        return self.fiber_states[..., 1].copy()

    @property
    def fiber_r_bars(self) -> FloatArray:
        """Return isotropic internal-variable snapshots."""
        return self.fiber_states[..., 2].copy()

    @property
    def fiber_damages(self) -> FloatArray:
        """Return damage snapshots."""
        return self.fiber_states[..., 3].copy()


def build_tower_cycle_phase_snapshots(
    response: NonlinearCyclicResponse,
) -> TowerCyclePhaseSnapshots:
    """
    Reorganize one nonlinear cyclic tower response into (n, tau) form.

    The original periodic history contains

        n_cycles * increments_per_cycle + 1

    stored points. For cycle ``n``, this function gathers

        start : start + increments_per_cycle + 1

    so both cycle endpoints are retained. Adjacent cycles therefore share and
    intentionally duplicate their common boundary state.
    """
    loading = response.loading

    n_cycles = int(
        loading.n_cycles
    )
    increments_per_cycle = int(
        loading.increments_per_cycle
    )

    expected_points = (
        n_cycles * increments_per_cycle + 1
    )

    if (
        response.fiber_strains.shape[0]
        != expected_points
    ):
        raise ValueError(
            "The response does not contain the expected number "
            "of complete periodic points."
        )

    if response.analysis_times.shape != (
        expected_points,
    ):
        raise ValueError(
            "analysis_times does not match "
            "the periodic history."
        )

    indices = _cycle_phase_indices(
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
    )

    phase_times = np.asarray(
        loading.times[
            : increments_per_cycle + 1
        ],
        dtype=np.float64,
    )

    phase_times = (
        phase_times
        - float(phase_times[0])
    )

    phase_fractions = (
        phase_times
        / float(loading.period)
    )

    phase_forces = np.asarray(
        loading.forces[
            : increments_per_cycle + 1
        ],
        dtype=np.float64,
    )

    return TowerCyclePhaseSnapshots(
        cycle_numbers=np.arange(
            1,
            n_cycles + 1,
            dtype=np.int64,
        ),
        phase_times=phase_times,
        phase_fractions=phase_fractions,
        phase_forces=phase_forces,
        analysis_times=_tensorize_history(
            values=response.analysis_times,
            indices=indices,
        ),
        nodal_displacements=_tensorize_history(
            values=response.nodal_displacements,
            indices=indices,
        ),
        fiber_strains=_tensorize_history(
            values=response.fiber_strains,
            indices=indices,
        ),
        fiber_stresses=_tensorize_history(
            values=response.fiber_stresses,
            indices=indices,
        ),
        fiber_states=_tensorize_history(
            values=response.fiber_states,
            indices=indices,
        ),
    )



def select_tower_cycle_range(
    snapshots: TowerCyclePhaseSnapshots,
    first_cycle: int,
    last_cycle: int,
) -> TowerCyclePhaseSnapshots:
    """
    Select an inclusive range of global cycle numbers.

    The returned object preserves the original global cycle numbering rather
    than renumbering the selected stage locally. For example, selecting cycles
    21 through 46 returns ``cycle_numbers == [21, ..., 46]``.

    Fast-phase coordinates and phase forces are common to every cycle and are
    therefore retained unchanged. Cycle-indexed arrays are sliced along their
    first axis. ``analysis_times`` remain the original absolute analysis times.

    Parameters
    ----------
    snapshots
        Source cycle-phase snapshot dataset.
    first_cycle
        First global cycle number to retain, inclusive.
    last_cycle
        Last global cycle number to retain, inclusive.

    Returns
    -------
    TowerCyclePhaseSnapshots
        A validated snapshot object containing only the requested cycle range.
    """
    if not isinstance(
        snapshots,
        TowerCyclePhaseSnapshots,
    ):
        raise TypeError(
            "snapshots must be a TowerCyclePhaseSnapshots object."
        )

    for value, name in (
        (first_cycle, "first_cycle"),
        (last_cycle, "last_cycle"),
    ):
        if isinstance(value, (bool, np.bool_)):
            raise TypeError(
                name + " must be an integer global cycle number."
            )
        if not isinstance(value, (int, np.integer)):
            raise TypeError(
                name + " must be an integer global cycle number."
            )

    first_cycle = int(first_cycle)
    last_cycle = int(last_cycle)

    if first_cycle < 1 or last_cycle < 1:
        raise ValueError(
            "Requested cycle numbers must be positive."
        )
    if first_cycle > last_cycle:
        raise ValueError(
            "first_cycle must not exceed last_cycle."
        )

    available_first = int(
        snapshots.cycle_numbers[0]
    )
    available_last = int(
        snapshots.cycle_numbers[-1]
    )

    if (
        first_cycle < available_first
        or last_cycle > available_last
    ):
        raise ValueError(
            "Requested cycle range [{}, {}] is outside the available "
            "global cycle range [{}, {}].".format(
                first_cycle,
                last_cycle,
                available_first,
                available_last,
            )
        )

    start_index = (
        first_cycle
        - available_first
    )
    stop_index = (
        last_cycle
        - available_first
        + 1
    )

    cycle_slice = slice(
        start_index,
        stop_index,
    )

    return TowerCyclePhaseSnapshots(
        cycle_numbers=snapshots.cycle_numbers[
            cycle_slice
        ],
        phase_times=snapshots.phase_times,
        phase_fractions=snapshots.phase_fractions,
        phase_forces=snapshots.phase_forces,
        analysis_times=snapshots.analysis_times[
            cycle_slice
        ],
        nodal_displacements=snapshots.nodal_displacements[
            cycle_slice
        ],
        fiber_strains=snapshots.fiber_strains[
            cycle_slice
        ],
        fiber_stresses=snapshots.fiber_stresses[
            cycle_slice
        ],
        fiber_states=snapshots.fiber_states[
            cycle_slice
        ],
    )

def save_tower_cycle_phase_snapshots(
    snapshots: TowerCyclePhaseSnapshots,
    file_path: PathLike,
) -> Path:
    """
    Save one complete tower cycle-phase snapshot dataset as compressed NPZ.

    Only primitive snapshot arrays are stored.

    Derived quantities such as plastic-strain channels, damage channels,
    cycle-increment fields, singular values, and low-rank diagnostics are
    intentionally not duplicated in the file.

    The dataset is first written to a temporary file and then moved into its
    final path. This reduces the chance that an interrupted write leaves a
    partially written file that appears to be a valid frozen FOM dataset.

    Parameters
    ----------
    snapshots
        Fully validated ``TowerCyclePhaseSnapshots`` object.
    file_path
        Destination path. The filename must use the ``.npz`` suffix.

    Returns
    -------
    pathlib.Path
        Final path of the saved snapshot dataset.
    """
    if not isinstance(
        snapshots,
        TowerCyclePhaseSnapshots,
    ):
        raise TypeError(
            "snapshots must be a "
            "TowerCyclePhaseSnapshots object."
        )

    path = Path(file_path)

    if path.suffix.lower() != ".npz":
        raise ValueError(
            "file_path must have the .npz suffix."
        )

    if not path.parent.exists():
        raise FileNotFoundError(
            "The parent directory of file_path "
            "does not exist."
        )

    if not path.parent.is_dir():
        raise NotADirectoryError(
            "The parent of file_path must be "
            "an existing directory."
        )

    temporary_path = path.with_name(
        path.name + ".tmp"
    )

    try:
        with temporary_path.open("wb") as file_object:
            np.savez_compressed(
                file_object,
                snapshot_format_version=np.asarray(
                    SNAPSHOT_FORMAT_VERSION,
                    dtype=np.int64,
                ),
                cycle_numbers=snapshots.cycle_numbers,
                phase_times=snapshots.phase_times,
                phase_fractions=snapshots.phase_fractions,
                phase_forces=snapshots.phase_forces,
                analysis_times=snapshots.analysis_times,
                nodal_displacements=(
                    snapshots.nodal_displacements
                ),
                fiber_strains=snapshots.fiber_strains,
                fiber_stresses=snapshots.fiber_stresses,
                fiber_states=snapshots.fiber_states,
            )

        temporary_path.replace(path)

    except Exception:
        if temporary_path.exists():
            temporary_path.unlink()
        raise

    return path


def load_tower_cycle_phase_snapshots(
    file_path: PathLike,
) -> TowerCyclePhaseSnapshots:
    """
    Load one tower cycle-phase snapshot dataset from compressed NPZ.

    Pickle loading is explicitly disabled.

    The arrays read from disk are passed back through
    ``TowerCyclePhaseSnapshots`` so the same cycle-grid, dimensional,
    finite-value, and damage-range validation used by newly generated FOM
    snapshots is also applied to datasets loaded from disk.

    Parameters
    ----------
    file_path
        Existing ``.npz`` snapshot dataset.

    Returns
    -------
    TowerCyclePhaseSnapshots
        Reconstructed and validated cycle-phase snapshot object.
    """
    path = Path(file_path)

    if path.suffix.lower() != ".npz":
        raise ValueError(
            "file_path must have the .npz suffix."
        )

    if not path.is_file():
        raise FileNotFoundError(
            "The requested snapshot file does not exist."
        )

    with np.load(
        path,
        allow_pickle=False,
    ) as archive:
        required_keys = {
            "snapshot_format_version",
            *_SNAPSHOT_ARRAY_KEYS,
        }

        available_keys = set(
            archive.files
        )

        missing_keys = (
            required_keys
            - available_keys
        )

        if missing_keys:
            missing_text = ", ".join(
                sorted(missing_keys)
            )
            raise ValueError(
                "Snapshot file is missing "
                "required arrays: "
                + missing_text
            )

        version_array = np.asarray(
            archive[
                "snapshot_format_version"
            ],
            dtype=np.int64,
        )

        if version_array.size != 1:
            raise ValueError(
                "snapshot_format_version must "
                "contain exactly one value."
            )

        version = int(
            version_array.reshape(-1)[0]
        )

        if (
            version
            != SNAPSHOT_FORMAT_VERSION
        ):
            raise ValueError(
                "Unsupported snapshot format version: "
                + str(version)
                + ". Expected version "
                + str(SNAPSHOT_FORMAT_VERSION)
                + "."
            )

        arrays = {
            key: np.asarray(
                archive[key]
            ).copy()
            for key in _SNAPSHOT_ARRAY_KEYS
        }

    return TowerCyclePhaseSnapshots(
        cycle_numbers=arrays[
            "cycle_numbers"
        ],
        phase_times=arrays[
            "phase_times"
        ],
        phase_fractions=arrays[
            "phase_fractions"
        ],
        phase_forces=arrays[
            "phase_forces"
        ],
        analysis_times=arrays[
            "analysis_times"
        ],
        nodal_displacements=arrays[
            "nodal_displacements"
        ],
        fiber_strains=arrays[
            "fiber_strains"
        ],
        fiber_stresses=arrays[
            "fiber_stresses"
        ],
        fiber_states=arrays[
            "fiber_states"
        ],
    )
