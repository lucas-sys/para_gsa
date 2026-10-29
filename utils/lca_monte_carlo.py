"""
lca_monte_carlo.py -- Monte Carlo sampling and parallel evaluation

This module provides:
- mc_sample(): Sample from Brightway uncertainty metadata distributions
- chunkify_rows(): Split sample array into chunks for multiprocessing
- _parse_score_key(): Parse mlca.scores key into (fu_label, method_tuple)
- recalculate_scores_for_sample(): Re-evaluate LCA scores for a single sample in worker
- recalculate_scores_with_background(): Same, but on top of a background draw
- init_worker(): Multiprocessing worker initialization function
- worker(): Multiprocessing worker processing function
- run_parallel_monte_carlo(): Entry point -- orchestrate sampling and parallel evaluation

Background uncertainty
----------------------
When ``include_background=True`` is passed to ``run_parallel_monte_carlo``,
each parameter sample is combined with one independent draw of the
background-database uncertainty (Brightway ``use_distributions=True``).
The foreground parameter perturbation is applied as a *delta* on top of the
background-sampled matrices, so the two uncertainty sources are not
double-counted.
"""

from __future__ import annotations

import copy
import math
import multiprocessing as mp
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

import bw2calc as bc
import bw2data as bd

from utils.lca_config import (
    VERBOSE,
    G_GROUP,
    G_PARAMETER_FORMULA_MAP,
    G_EXCHANGE_FORMULA_MAP,
    G_MULTILCA,
    G_BASE_TECH,
    G_BASE_BIO,
    G_METHODS,
    G_ACTIVITY_PARAMETER_NAMES,
    G_BASELINE_PARAM_CONTEXT,
    G_BASELINE_UPDATES,
    G_BASELINE_GROUPED_UPDATES,
    G_DATA_OBJS,
    G_METHOD_CONFIG_DICT,
    G_FU_LABELS,
    G_DEMANDS,
    G_BG_MULTILCA,
    G_BG_MULTILCA_ITER,
    Timer,
    _log,
)
from utils.lca_project_setup import set_project, initialize_multilca
from utils.lca_parameters import (
    get_activity_parameters,
    build_parameter_formula_map,
    build_exchange_formula_map,
    resolve_parameter_context,
    evaluate_exchange_formulas,
)
from utils.lca_matrices import (
    aggregate_updates_by_cell,
    apply_grouped_updates_to_matrices,
    compute_param_delta,
    apply_param_delta_to_matrices,
)


# -----------------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------------

def mc_sample(
    params: Sequence,
    n_samples: int = 1000,
    random_seed: Optional[int] = None,
) -> np.ndarray:
    """Sample Brightway uncertainty metadata with optional reproducible RNG seed."""
    rng = np.random.default_rng(random_seed)
    samples = np.zeros((n_samples, len(params)))

    for i, par in enumerate(params):
        p = par.dict
        code = p.get("uncertainty type", 0)
        loc = p.get("loc", 0)
        scale = p.get("scale", 1)
        shape = p.get("shape", 1)
        p_min = p.get("minimum", 0)
        p_max = p.get("maximum", 1)
        shape2 = p.get("shape2", 1)

        if code in (0, 1):
            samples[:, i] = par.amount
        elif code == 2:
            samples[:, i] = stats.lognorm.rvs(
                s=scale, loc=0, scale=np.exp(loc), size=n_samples, random_state=rng
            )
        elif code == 3:
            samples[:, i] = stats.norm.rvs(
                loc=loc, scale=scale, size=n_samples, random_state=rng
            )
        elif code == 4:
            samples[:, i] = stats.uniform.rvs(
                loc=p_min, scale=p_max - p_min, size=n_samples, random_state=rng
            )
        elif code == 5:
            width = p_max - p_min
            c = 0.5 if width == 0 else (loc - p_min) / width
            samples[:, i] = stats.triang.rvs(
                c=c, loc=p_min, scale=max(width, 1e-30), size=n_samples, random_state=rng
            )
        elif code == 6:
            samples[:, i] = stats.bernoulli.rvs(p=loc, size=n_samples, random_state=rng)
        elif code == 7:
            samples[:, i] = stats.randint.rvs(
                low=p_min, high=p_max + 1, size=n_samples, random_state=rng
            )
        elif code == 8:
            samples[:, i] = stats.weibull_min.rvs(
                c=shape, loc=loc, scale=scale, size=n_samples, random_state=rng
            )
        elif code == 9:
            samples[:, i] = stats.gamma.rvs(
                a=shape, loc=loc, scale=scale, size=n_samples, random_state=rng
            )
        elif code == 10:
            samples[:, i] = stats.beta.rvs(
                a=shape, b=shape2, loc=loc, scale=scale, size=n_samples, random_state=rng
            )
        elif code == 11:
            samples[:, i] = stats.genextreme.rvs(
                c=shape, loc=loc, scale=scale, size=n_samples, random_state=rng
            )
        elif code == 12:
            samples[:, i] = stats.t.rvs(
                df=shape, loc=loc, scale=scale, size=n_samples, random_state=rng
            )
        else:
            samples[:, i] = np.nan
            print(f"Warning: Parameter index {i} has unknown distribution code {code}.")

    return samples


# -----------------------------------------------------------------------------
# Chunk utilities
# -----------------------------------------------------------------------------

def chunkify_rows(sample_array: np.ndarray, n_chunks: int) -> List[List[Tuple[int, np.ndarray]]]:
    """Split sample rows into chunks for multiprocessing."""
    indexed = [(i, sample_array[i, :]) for i in range(len(sample_array))]
    n_chunks = max(1, int(n_chunks))
    chunk_size = max(1, math.ceil(len(indexed) / n_chunks))
    return [indexed[i : i + chunk_size] for i in range(0, len(indexed), chunk_size)]


# -----------------------------------------------------------------------------
# Score key parsing
# -----------------------------------------------------------------------------

def _parse_score_key(key):
    """Parse mlca.scores key into (fu_label, method_tuple)."""
    if isinstance(key, tuple) and len(key) == 2:
        a, b = key
        if isinstance(a, tuple) and isinstance(b, str):
            return b, a
        if isinstance(a, str) and isinstance(b, tuple):
            return a, b
    return None, key


# -----------------------------------------------------------------------------
# Score recomputation helpers (run inside worker)
# -----------------------------------------------------------------------------

def _reset_mlca_cached_attributes(mlca) -> None:
    """Clear cached attributes on a MultiLCA so lci/lcia can be re-run."""
    for attr in (
        "solver",
        "supply_array",
        "inventory",
        "characterized_inventory",
        "normalized_inventory",
        "weighted_inventory",
    ):
        if hasattr(mlca, attr):
            try:
                delattr(mlca, attr)
            except Exception:
                pass


def _recompute_scores(mlca, tech_matrix, bio_matrix) -> List[dict]:
    """Set custom matrices on ``mlca``, re-run LCI/LCIA, and return tidy score rows."""
    mlca.technosphere_matrix = tech_matrix.copy().tocsr()
    mlca.biosphere_matrix = bio_matrix.copy().tocsr()
    _reset_mlca_cached_attributes(mlca)
    mlca.lci()
    mlca.lcia()

    rows = []
    for key, score in mlca.scores.items():
        fu_label, method = _parse_score_key(key)
        rows.append(
            {
                "Method": method[2] if isinstance(method, tuple) and len(method) > 2 else str(method),
                "Method full": method,
                "Functional unit": fu_label,
                "LCA Score": float(score),
            }
        )
    return rows


# -----------------------------------------------------------------------------
# Background-uncertainty helpers
# -----------------------------------------------------------------------------

def _prepare_bg_data_objs(data_objs, exchange_formula_map):
    """
    Deep-copy ``data_objs`` and disable uncertainty on parameterized cells.

    With ``use_distributions=True`` Brightway samples *all* exchanges, including
    the foreground parameterized ones. That would make the background-sampled
    matrix deviate from baseline at parameterized cells, and the subsequent
    overlay of the parameter delta would double-count the perturbation
    (near-singular technosphere matrix, exploding LCA scores).

    This helper prevents background sampling from touching parameterized cells:
    1. deep-copies the data objects (the originals used by ``G_MULTILCA`` stay intact)
    2. sets ``uncertainty type = 0`` and clears distribution fields for every
       matrix cell mapped by ``exchange_formula_map``.
    """
    param_cells = {
        (item["matrix"], int(item["row_index"]), int(item["col_index"]))
        for item in exchange_formula_map
    }
    bg_data_objs = copy.deepcopy(data_objs)

    for obj in bg_data_objs:
        if not hasattr(obj, "data"):
            continue
        if not isinstance(obj.data, dict):
            # Some datapackages store ``data`` as a list; they do not carry
            # per-matrix exchange rows, so there is nothing to patch.
            continue
        for matrix_name, data_list in obj.data.items():
            if not isinstance(data_list, list):
                continue
            for row_data in data_list:
                if not isinstance(row_data, dict):
                    continue
                row_idx = row_data.get("row")
                col_idx = row_data.get("col")
                if row_idx is None or col_idx is None:
                    continue
                key = (matrix_name, int(row_idx), int(col_idx))
                if key in param_cells:
                    row_data["uncertainty type"] = 0
                    for field in ("loc", "scale", "shape", "minimum", "maximum", "shape2"):
                        row_data[field] = None
    return bg_data_objs


def bg_matrix_iterator(mlca):
    """
    Advance the background matrix distributions without running LCI/LCIA.

    ``MultiLCA.__next__`` rebuilds the matrices and immediately runs
    ``_calculation()`` (lci + lcia). When a parameter delta still has to be
    overlaid and the scores recomputed, that first LCI/LCIA round is wasted.
    This generator only advances the distribution iterators; the caller is
    responsible for overlaying the delta and calling ``_recompute_scores()``
    (one LCI/LCIA round per sample).
    """
    while True:
        mlca._delete_solver_state()
        for matrix in mlca.matrix_labels:
            if hasattr(mlca, matrix):
                next(getattr(mlca, matrix))
        for matrix_dict in mlca.matrix_list_labels:
            if hasattr(mlca, matrix_dict):
                next(getattr(mlca, matrix_dict))
        if hasattr(mlca, "after_matrix_iteration"):
            mlca.after_matrix_iteration()
        if bc.PYPARDISO:
            mlca.technosphere_matrix = mlca.technosphere_matrix.tocsr()
        yield mlca


# -----------------------------------------------------------------------------
# Per-sample calculation (runs inside worker)
# -----------------------------------------------------------------------------

def _evaluate_sample(sample_row: np.ndarray) -> dict:
    """Convert one parameter sample row into aggregated cell-wise updates."""
    sampled_values = {
        G_ACTIVITY_PARAMETER_NAMES[j]: sample_row[j]
        for j in range(len(G_ACTIVITY_PARAMETER_NAMES))
    }
    param_context = resolve_parameter_context(
        G_GROUP,
        sampled_values_by_name=sampled_values,
        parameter_formula_map=G_PARAMETER_FORMULA_MAP,
    )
    new_updates = evaluate_exchange_formulas(G_EXCHANGE_FORMULA_MAP, param_context)
    return aggregate_updates_by_cell(new_updates)


def _calc_nobg(sample_index: int, new_grouped: dict) -> List[dict]:
    """No background uncertainty: overlay parameter updates on the baseline matrices."""
    tech_matrix, bio_matrix = apply_grouped_updates_to_matrices(
        base_tech=G_BASE_TECH,
        base_bio=G_BASE_BIO,
        baseline_grouped_updates=G_BASELINE_GROUPED_UPDATES,
        new_grouped_updates=new_grouped,
    )
    rows = _recompute_scores(G_MULTILCA, tech_matrix, bio_matrix)
    for r in rows:
        r["Iteration"] = sample_index + 1
        # Keep the historical Monte Carlo output schema (no "Method full").
        r.pop("Method full", None)
    return rows


def _calc_bg(sample_index: int, new_grouped: dict) -> List[dict]:
    """Background uncertainty: draw one background state, overlay the parameter delta, score."""
    param_delta = compute_param_delta(G_BASELINE_GROUPED_UPDATES, new_grouped)
    mlca = G_BG_MULTILCA
    bg_iter = G_BG_MULTILCA_ITER

    next(bg_iter)  # advance background matrices only; no LCI/LCIA inside

    # Use the underlying ``MappedMatrix.matrix`` rather than
    # ``mlca.technosphere_matrix``: ``_recompute_scores`` replaces the cached
    # ``mlca.technosphere_matrix`` attribute with a copy, while the mapped
    # matrix stays in sync with ``next(bg_iter)``.
    tech_matrix, bio_matrix = apply_param_delta_to_matrices(
        mlca.technosphere_mm.matrix,
        mlca.biosphere_mm.matrix,
        param_delta,
    )
    rows = _recompute_scores(mlca, tech_matrix, bio_matrix)
    for r in rows:
        r["Iteration"] = sample_index + 1
        r.pop("Method full", None)
    return rows


def recalculate_scores_for_sample(
    sample_index: int,
    base_tech,
    base_bio,
    baseline_grouped_updates,
    new_grouped_updates,
) -> List[dict]:
    """Evaluate one MC sample row (no background draw) and return tidy score rows."""
    tech_matrix, bio_matrix = apply_grouped_updates_to_matrices(
        base_tech=base_tech,
        base_bio=base_bio,
        baseline_grouped_updates=baseline_grouped_updates,
        new_grouped_updates=new_grouped_updates,
    )
    rows = _recompute_scores(G_MULTILCA, tech_matrix, bio_matrix)
    for r in rows:
        r["Iteration"] = sample_index + 1
        # Keep the historical MC output schema (no "Method full" column).
        r.pop("Method full", None)
    return rows


def recalculate_scores_with_background(
    sample_index: int,
    new_grouped_updates,
) -> List[dict]:
    """
    Evaluate one MC sample row on top of a fresh background uncertainty draw.

    Public wrapper around the internal background path; requires a worker that
    was initialized with ``include_background=True``.
    """
    if G_BG_MULTILCA is None:
        raise RuntimeError(
            "recalculate_scores_with_background() requires a worker initialized "
            "with include_background=True"
        )
    return _calc_bg(sample_index, new_grouped_updates)


# -----------------------------------------------------------------------------
# Multiprocessing worker
# -----------------------------------------------------------------------------

def init_worker(
    project_name: str,
    group_name: str,
    fus_keys: Sequence[Tuple[str, str]],
    methods: Sequence[tuple],
    activity_parameter_names: Sequence[str],
    verbose: bool = False,
    include_background: bool = False,
):
    """Initialize a Monte Carlo multiprocessing worker."""
    global VERBOSE
    global G_GROUP, G_PARAMETER_FORMULA_MAP, G_EXCHANGE_FORMULA_MAP
    global G_MULTILCA, G_BASE_TECH, G_BASE_BIO, G_METHODS
    global G_ACTIVITY_PARAMETER_NAMES
    global G_BASELINE_PARAM_CONTEXT, G_BASELINE_UPDATES, G_BASELINE_GROUPED_UPDATES
    global G_DATA_OBJS, G_METHOD_CONFIG_DICT, G_FU_LABELS, G_DEMANDS
    global G_BG_MULTILCA, G_BG_MULTILCA_ITER

    VERBOSE = verbose
    with Timer("worker_init"):
        set_project(project_name)
        G_GROUP = group_name
        G_METHODS = list(methods)
        G_ACTIVITY_PARAMETER_NAMES = list(activity_parameter_names)
        fus = [bd.get_node(database=db, code=code) for db, code in fus_keys]
        G_PARAMETER_FORMULA_MAP = build_parameter_formula_map(group_name)

        setup = initialize_multilca(fus, G_METHODS)
        G_MULTILCA = setup.mlca
        G_BASE_TECH = setup.base_tech
        G_BASE_BIO = setup.base_bio
        G_EXCHANGE_FORMULA_MAP = build_exchange_formula_map(group_name, G_MULTILCA)

        # Cache baseline once per worker
        G_BASELINE_PARAM_CONTEXT = resolve_parameter_context(
            G_GROUP,
            sampled_values_by_name={},
            parameter_formula_map=G_PARAMETER_FORMULA_MAP,
        )
        G_BASELINE_UPDATES = evaluate_exchange_formulas(
            G_EXCHANGE_FORMULA_MAP,
            G_BASELINE_PARAM_CONTEXT,
        )
        G_BASELINE_GROUPED_UPDATES = aggregate_updates_by_cell(G_BASELINE_UPDATES)

        G_DATA_OBJS = setup.data_objs
        G_METHOD_CONFIG_DICT = setup.method_config_dict
        G_FU_LABELS = setup.fu_labels
        G_DEMANDS = setup.demands

        if include_background:
            with Timer("bg_mlca_init"):
                bg_data_objs = _prepare_bg_data_objs(G_DATA_OBJS, G_EXCHANGE_FORMULA_MAP)
                G_BG_MULTILCA = bc.MultiLCA(
                    demands=G_DEMANDS,
                    method_config=G_METHOD_CONFIG_DICT,
                    data_objs=bg_data_objs,
                    use_distributions=True,
                )
                G_BG_MULTILCA.lci()
                G_BG_MULTILCA.lcia()
                G_BG_MULTILCA_ITER = bg_matrix_iterator(G_BG_MULTILCA)
            _log(
                f"initialized group={group_name} fus={len(fus)} methods={len(methods)} "
                f"param_formulas={len(G_PARAMETER_FORMULA_MAP)} "
                f"exchange_formulas={len(G_EXCHANGE_FORMULA_MAP)} "
                f"BACKGROUND=ON"
            )
        else:
            G_BG_MULTILCA = None
            G_BG_MULTILCA_ITER = None
            _log(
                f"initialized group={group_name} fus={len(fus)} methods={len(methods)} "
                f"param_formulas={len(G_PARAMETER_FORMULA_MAP)} "
                f"exchange_formulas={len(G_EXCHANGE_FORMULA_MAP)} "
                f"BACKGROUND=OFF"
            )


def worker(sample_chunk: Sequence[Tuple[int, np.ndarray]]) -> List[dict]:
    """Monte Carlo multiprocessing worker: process a chunk of sample rows."""
    include_background = G_BG_MULTILCA is not None
    rows = []
    with Timer(f"chunk size={len(sample_chunk)} bg={include_background}"):
        for sample_index, sample_row in sample_chunk:
            new_grouped_updates = _evaluate_sample(sample_row)
            if include_background:
                rows.extend(_calc_bg(sample_index, new_grouped_updates))
            else:
                rows.extend(_calc_nobg(sample_index, new_grouped_updates))
    return rows


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

def run_parallel_monte_carlo(
    project_name: str,
    group_name: str,
    functional_units: Sequence,
    methods: Sequence[tuple],
    n_samples: int = 1000,
    n_workers: Optional[int] = None,
    random_seed: Optional[int] = 42,
    verbose: bool = True,
    export_parameter_samples: Optional[str] = None,
    include_background: bool = False,
) -> pd.DataFrame:
    """
    Run parallel Monte Carlo LCA with parameter uncertainty.

    Parameters
    ----------
    include_background : bool, default False
        If ``True``, every parameter sample is combined with one independent
        draw of the background-database uncertainty
        (Brightway ``use_distributions=True``). Parameterized foreground cells
        are masked out of the background draw so the two sources are not
        double-counted.

    Returns
    -------
    pd.DataFrame
        Tidy DataFrame with columns:
        ``Iteration, Method, Functional unit, LCA Score``
    """
    set_project(project_name)
    params = get_activity_parameters(group_name)
    param_names = [p.name for p in params]
    sample_array = mc_sample(params, n_samples=n_samples, random_seed=random_seed)

    if export_parameter_samples:
        pd.DataFrame(sample_array, columns=param_names).to_excel(export_parameter_samples, index=False)

    fus_keys = [(fu["database"], fu["code"]) for fu in functional_units]
    if n_workers is None:
        n_workers = min(max(1, mp.cpu_count() - 1), len(sample_array))
    chunks = chunkify_rows(sample_array, n_workers)

    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=len(chunks),
        initializer=init_worker,
        initargs=(
            project_name,
            group_name,
            fus_keys,
            methods,
            param_names,
            verbose,
            include_background,
        ),
    ) as pool:
        nested_results = pool.map(worker, chunks)

    return pd.DataFrame([item for chunk in nested_results for item in chunk])
