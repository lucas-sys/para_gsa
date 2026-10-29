"""
lca_config.py -- Global configuration, logging utilities, timer, and LCASetup dataclass

Shared infrastructure for all modules:
- VERBOSE global logging switch
- Multiprocessing worker global variable declarations
- LCASetup dataclass (encapsulates all context needed for LCA computation)
- _log() logging function and Timer context manager
"""

from __future__ import annotations

import multiprocessing as mp
import os
import time
from dataclasses import dataclass
from typing import Dict, Sequence

# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------

VERBOSE = False

# Globals initialized once in each multiprocessing worker (Monte Carlo)
G_GROUP = None
G_PARAMETER_FORMULA_MAP = None
G_EXCHANGE_FORMULA_MAP = None
G_MULTILCA = None
G_BASE_TECH = None
G_BASE_BIO = None
G_METHODS = None
G_ACTIVITY_PARAMETER_NAMES = None

# Globals initialized once in each multiprocessing worker (Sobol)
G_SOBOL_PARAMETER_NAMES = None

# Worker baseline cache (Monte Carlo)
G_BASELINE_PARAM_CONTEXT = None
G_BASELINE_UPDATES = None
G_BASELINE_GROUPED_UPDATES = None

# Background-uncertainty globals (shared by Monte Carlo / Sobol / FAST workers)
# These are only populated when ``include_background=True`` is passed to the
# corresponding ``run_parallel_*`` entry point.
G_DATA_OBJS = None
G_METHOD_CONFIG_DICT = None
G_FU_LABELS = None
G_DEMANDS = None
G_BG_MULTILCA = None       # use_distributions=True MultiLCA instance
G_BG_MULTILCA_ITER = None  # iterator over G_BG_MULTILCA (persisted across draws)


# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class LCASetup:
    """Encapsulates all state needed for parameter-perturbation LCA runs."""
    project_name: str
    group_name: str
    functional_units: Sequence
    methods: Sequence[tuple]
    fu_labels: Dict[object, str]
    demands: Dict[str, Dict[int, float]]
    method_config_dict: dict
    data_objs: list
    mlca: object          # MultiLCA instance
    base_tech: object     # sparse technosphere matrix (CSR)
    base_bio: object      # sparse biosphere matrix (CSR)


# -----------------------------------------------------------------------------
# Logging and timing utilities
# -----------------------------------------------------------------------------

def _log(message: str) -> None:
    if VERBOSE:
        proc = mp.current_process()
        print(f"[{proc.name} pid={os.getpid()}] {message}", flush=True)


class Timer:
    """Context manager for timing code blocks (logs when VERBOSE is on)."""

    def __init__(self, label: str):
        self.label = label
        self.t0 = None

    def __enter__(self):
        self.t0 = time.perf_counter()
        _log(f"START {self.label}")
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.perf_counter() - self.t0
        _log(f"END {self.label} ({dt:.2f}s)")
