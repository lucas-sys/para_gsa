"""
lca_matrices.py -- Sparse matrix update utilities

This module safely updates technosphere/biosphere sparse matrices under parameter perturbation:
- aggregate_updates_by_cell(): Sum exchange updates mapping to the same matrix cell
- apply_updates_to_matrices(): Apply parameterized update deltas to base matrices
- compute_param_delta(): Per-cell delta between new parameter values and baseline
- apply_param_delta_to_matrices(): Overlay a parameter delta on arbitrary matrices

Core formula (preserving non-parameterized contributions):
    new_cell = base_cell - baseline_param_sum + new_param_sum

Background-uncertainty overlay formula:
    final_cell = bg_cell + param_delta
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Mapping, Sequence, Tuple


def aggregate_updates_by_cell(
    updates: Sequence[dict],
) -> Dict[Tuple[str, int, int], float]:
    """
    Sum exchange updates that map to the same matrix cell.

    Returns dict keyed by (matrix_name, row_index, col_index).
    """
    grouped = defaultdict(float)
    for item in updates:
        key = (item["matrix"], int(item["row_index"]), int(item["col_index"]))
        grouped[key] += float(item["matrix_value"])
    return dict(grouped)


def apply_updates_to_matrices(
    base_tech,
    base_bio,
    baseline_grouped_updates: Mapping[Tuple[str, int, int], float],
    new_grouped_updates: Mapping[Tuple[str, int, int], float],
):
    """
    Update matrix cells by replacing only the parameterized portion:

        new_cell = base_cell - baseline_param_sum + new_param_sum

    This preserves any non-parameterized contribution already present in the
    stored matrix and correctly handles multiple parameterized exchanges per cell.

    Parameters
    ----------
    base_tech : sparse matrix
        Baseline technosphere matrix (will NOT be mutated).
    base_bio : sparse matrix
        Baseline biosphere matrix (will NOT be mutated).
    baseline_grouped_updates : dict
        {('technosphere_matrix'|'biosphere_matrix', row, col): sum_of_baseline_param_values}
    new_grouped_updates : dict
        Same shape as above but with new parameter values.

    Returns
    -------
    (tech, bio) : tuple of sparse CSR matrices
    """
    tech = base_tech.tolil(copy=True)
    bio = base_bio.tolil(copy=True)

    touched_keys = set(baseline_grouped_updates) | set(new_grouped_updates)

    for matrix, r, c in touched_keys:
        old_param = baseline_grouped_updates.get((matrix, r, c), 0.0)
        new_param = new_grouped_updates.get((matrix, r, c), 0.0)

        if matrix == "technosphere_matrix":
            base_val = float(base_tech[r, c])
            tech[r, c] = base_val - old_param + new_param
        elif matrix == "biosphere_matrix":
            base_val = float(base_bio[r, c])
            bio[r, c] = base_val - old_param + new_param

    return tech.tocsr(), bio.tocsr()


# Backward-compatible alias (used by Monte Carlo module)
apply_grouped_updates_to_matrices = apply_updates_to_matrices


def compute_param_delta(
    baseline_grouped_updates: Mapping[Tuple[str, int, int], float],
    new_grouped_updates: Mapping[Tuple[str, int, int], float],
) -> Dict[Tuple[str, int, int], float]:
    """
    Compute per-cell delta = new_param_sum - baseline_param_sum.

    This delta can be computed once per parameter sample and then re-used when
    the same sample has to be applied on top of several different
    background-sampled matrices (Monte Carlo / Sobol / FAST with background
    uncertainty).

    Parameters
    ----------
    baseline_grouped_updates : dict
        {('technosphere_matrix'|'biosphere_matrix', row, col): baseline param sum}
    new_grouped_updates : dict
        Same shape, with the sampled parameter values.

    Returns
    -------
    dict
        Cell-wise delta, ready for ``apply_param_delta_to_matrices``.
    """
    keys = set(baseline_grouped_updates) | set(new_grouped_updates)
    return {
        k: new_grouped_updates.get(k, 0.0) - baseline_grouped_updates.get(k, 0.0)
        for k in keys
    }


def apply_param_delta_to_matrices(
    tech_matrix,
    bio_matrix,
    param_delta: Mapping[Tuple[str, int, int], float],
):
    """
    Apply a pre-computed parameter delta on top of existing matrices.

    Unlike ``apply_updates_to_matrices`` this does **not** subtract a baseline
    parameter sum: it is meant to be used when the input matrices have already
    been perturbed by a background-database uncertainty draw and only the
    foreground parameter delta still has to be overlaid.

    Parameters
    ----------
    tech_matrix, bio_matrix : sparse matrix
        Current matrices (e.g. after a background draw; will NOT be mutated).
    param_delta : dict
        {(matrix_name, row, col): delta_value} produced by ``compute_param_delta``.

    Returns
    -------
    (tech, bio) : tuple of sparse CSR matrices
    """
    tech = tech_matrix.tolil(copy=True)
    bio = bio_matrix.tolil(copy=True)

    for (matrix, r, c), delta in param_delta.items():
        if matrix == "technosphere_matrix":
            tech[r, c] = float(tech[r, c]) + delta
        elif matrix == "biosphere_matrix":
            bio[r, c] = float(bio[r, c]) + delta

    return tech.tocsr(), bio.tocsr()
