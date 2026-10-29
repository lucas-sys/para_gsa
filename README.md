## Contributors

<a href="https://github.com/lucas-sys/para_gsa/graphs/contributors">
  <img src="https://img.shields.io/github/contributors/lucas-sys/para_gsa?style=for-the-badge&logo=github" alt="Contributors"/>
</a>

<a href="https://github.com/lucas-sys/para_gsa/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=lucas-sys/para_gsa" alt="Contributors"/>
</a>



# LCA Global Sensitivity Analysis (GSA) Toolkit

A modular Python toolkit for performing **parameter-perturbation Life Cycle Assessment (LCA)** with **Global Sensitivity Analysis (GSA)**, built on [Brightway2.5](https://brightway.dev/), [SALib](https://salib.readthedocs.io/), and Python multiprocessing.

This toolkit was developed in collaboration by Rebecca Belfiore, Christhel Andrade Diaz, Weier Liu and Nanlin Liao. Code X was used to improve the code efficiency, and WorkBuddy was used to organize all files.

## Features

- **Deterministic LCA** -- Compute baseline environmental impact scores across multiple functional units and impact categories
- **Contribution Analysis** -- Recursive Tier-2 contribution analysis to identify key process contributions
- **Formula Validation** -- Verify stored exchange amounts against parameterized formulas
- **OAT Sensitivity** -- One-at-a-time finite-difference sensitivity analysis with elasticity coefficients
- **Analytical Variance** -- Variance decomposition from uncertainty metadata (13 distribution types)
- **Monte Carlo Simulation** -- Parallel stochastic sampling with configurable worker count
- **Sobol GSA** -- First-order (S1) and total-order (ST) Sobol indices via SALib with parallel evaluation
- **FAST / RBD-FAST GSA** -- Global sensitivity analysis via Fourier Amplitude Sensitivity Test (FAST) and its Random Balance Designs variant (RBD-FAST), with parallel evaluation
- **Background Uncertainty (BG)** -- Optionally combine foreground parameter uncertainty with background-database uncertainty (Brightway `use_distributions=True`) for Monte Carlo, Sobol, and FAST / RBD-FAST

## Project Structure

```
.
├── para_gsa/                            # Core LCA modules
│   ├── __init__.py                   # Public API re-exports
│   ├── lca_config.py                 # Configuration, logging, Timer, LCASetup
│   ├── lca_project_setup.py          # Brightway project init, MultiLCA
│   ├── lca_parameters.py             # Parameter queries, formula resolution
│   ├── lca_matrices.py               # Sparse matrix update utilities
│   ├── lca_contribution.py           # Contribution analysis, OAT, analytical variance
│   ├── lca_monte_carlo.py            # Monte Carlo sampling & parallel evaluation
│   ├── lca_sobol.py                  # Sobol GSA (build, sample, evaluate, analyze)
│   ├── lca_fast.py                   # FAST / RBD-FAST GSA (build, sample, evaluate, analyze)
│   └── lca_gsa_helper.py             # Backward-compatible re-export layer
├── data/                             # Input databases (user-provided)
│   └── ...                           # Place your Excel database files here
├── results/                          # Output files (organized by analysis type)
│   ├── deterministic/                # Deterministic LCA results
│   ├── contribution/                 # Contribution analysis results
│   ├── sensitivity_oat/              # OAT sensitivity results
│   ├── sensitivity_analytical/       # Analytical variance results
│   ├── monte_carlo/                  # Monte Carlo results
│   ├── sobol/                        # Sobol GSA results
│   └── fast/                         # FAST / RBD-FAST GSA results
├── gsa_workflow.ipynb                # Interactive notebook (full workflow)
├── init_databases.py                 # Database initialization script
├── test_full_workflow.py             # End-to-end test
├── test_run.py                       # Quick validation test
├── test_sobol_only.py                # Sobol-only test
├── requirements.txt
├── LICENSE
└── README.md
```

## Prerequisites

- Python >= 3.10
- [Brightway2.5](https://brightway.dev/) (`bw2data`, `bw2calc`, `bw2analyzer`; `bw2io` only for the database-initialization helper)
- An ecoinvent database (requires ecoinvent license)

## Installation

### From PyPI

```bash
pip install para-gsa
```

The runtime dependencies (Brightway 2.5, SALib, NumPy, SciPy, pandas, asteval) are
installed automatically. You still need your own Brightway project and, for
background data, an ecoinvent licence.

To also install the helper used to import ecoinvent / custom databases:

```bash
pip install "para-gsa[io]"
```

Verify the install:

```bash
python -c "import para_gsa; print(para_gsa.__version__)"
```

### From source (development)

```bash
git clone https://github.com/lucas-sys/para_gsa.git
cd para_gsa
pip install -e ".[dev]"
```

### Conda environment

```bash
conda create -n bw25 python=3.10
conda activate bw25
pip install para-gsa
```

## Configuration

Before running any analysis, you need to configure the following parameters:

| Parameter                     | Description                                    | Where to set                                 |
| ----------------------------- | ---------------------------------------------- | -------------------------------------------- |
| `PROJECT_NAME`                | Your Brightway project name                    | `init_databases.py`, test scripts, notebook  |
| `FOREGROUND_DB`               | Name of your foreground database               | Test scripts, notebook                       |
| `GROUP_NAME`                  | ActivityParameter group name                   | Test scripts, notebook                       |
| `FU_CODE`                     | Activity code for the functional unit          | Test scripts, notebook                       |
| `METHODS`                     | List of impact method tuples                   | Test scripts, notebook                       |
| `ECOINVENT_USERNAME/PASSWORD` | ecoinvent credentials                          | `init_databases.py` or environment variables |
| `CUSTOM_DATABASES`            | List of (db_name, xlsx_path, match_dbs) tuples | `init_databases.py`                          |

All parameters use `your_*` placeholders in the template files. Replace them with your actual values before running.

## Database Initialization

```bash
# Set ecoinvent credentials (recommended: use environment variables)
export ECOINVENT_USERNAME="your_username"
export ECOINVENT_PASSWORD="your_password"

# Or edit init_databases.py directly and update the configuration section

# Run the initialization script
python init_databases.py
```

Place your custom database Excel files in the `data/` directory and add entries to the `CUSTOM_DATABASES` list in `init_databases.py`.

## Quick Start

### In Python

```python
from para_gsa import (
    set_project, initialize_multilca, deterministic_lca,
    oat_sensitivity, run_parallel_monte_carlo,
    build_sobol_problem, generate_sobol_samples,
    run_parallel_sobol_from_samples, sobol_indices_from_results,
)

import bw2data as bd

# Setup
PROJECT_NAME = 'your_project_name'
FOREGROUND_DB = 'your_foreground_db'
GROUP_NAME = 'your_parameter_group'

set_project(PROJECT_NAME)
db = bd.Database(FOREGROUND_DB)
fu = [db.get(code='your_fu_code')]

methods = [
    ('your_ecoinvent_db', 'EF v3.1', 'climate change',
     'global warming potential (GWP100)'),
]

# Deterministic LCA
df = deterministic_lca(fu, methods)

# Monte Carlo (1000 iterations, 2 workers)
mc_results = run_parallel_monte_carlo(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    n_samples=1000,
    n_workers=2,
)

# Sobol GSA
problem, params = build_sobol_problem(GROUP_NAME)
samples = generate_sobol_samples(problem, N=512)
sobol_results = run_parallel_sobol_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    sobol_problem=problem,
    sobol_samples=samples,
    n_workers=2,
)

# FAST / RBD-FAST GSA
from para_gsa import lca_fast

fast_problem, fast_params = lca_fast.build_fast_problem(GROUP_NAME)

# --- Classic FAST (S1 + ST) ---
fast_samples = lca_fast.generate_fast_samples(fast_problem, M=4, method='fast')
fast_results = lca_fast.run_parallel_fast_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    fast_problem=fast_problem,
    fast_samples=fast_samples,
    n_workers=2,
)
Si_fast = lca_fast.fast_indices_from_results(
    fast_problem=fast_problem,
    fast_results_df=fast_results,
    functional_unit=fu[0]['name'],
    method=methods[0],
    fast_samples=fast_samples,
    M=4,
    analysis_method='fast',
)
df_fast = lca_fast.fast_indices_to_dataframe(Si_fast, fast_problem, analysis_method='fast')

# --- RBD-FAST (S1 only) ---
rbd_samples = lca_fast.generate_fast_samples(fast_problem, M=4, method='rbd_fast', seed=42)
rbd_results = lca_fast.run_parallel_fast_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    fast_problem=fast_problem,
    fast_samples=rbd_samples,
    n_workers=2,
)
Si_rbd = lca_fast.fast_indices_from_results(
    fast_problem=fast_problem,
    fast_results_df=rbd_results,
    functional_unit=fu[0]['name'],
    method=methods[0],
    fast_samples=rbd_samples,
    M=4,
    analysis_method='rbd_fast',
)
df_rbd = lca_fast.fast_indices_to_dataframe(Si_rbd, fast_problem, analysis_method='rbd_fast')

# One-click workflow
result = lca_fast.run_full_fast_workflow(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    M=4,
    n_workers=2,
    method='rbd_fast',  # or 'fast'
)
```

### In Jupyter Notebook

Open `gsa_workflow.ipynb` and run cells sequentially. The notebook covers the complete workflow from database setup through Sobol and FAST GSA. Before running, update the configuration cells with your project-specific values.

## Workflow Sections

| Section | Description               | Output Directory                  |
| ------- | ------------------------- | --------------------------------- |
| 0       | Brightway project setup   | --                                |
| 1       | Impact method selection   | --                                |
| 2       | Functional unit selection | --                                |
| 3       | Deterministic LCA         | `results/deterministic/`          |
| 4       | Contribution analysis     | `results/contribution/`           |
| 5       | OAT sensitivity           | `results/sensitivity_oat/`        |
| 6       | Analytical variance       | `results/sensitivity_analytical/` |
| 7       | Monte Carlo simulation    | `results/monte_carlo/`            |
| 8       | Sobol GSA                 | `results/sobol/`                  |
| 10      | FAST / RBD-FAST GSA       | `results/fast/`                   |

## Supported Uncertainty Distributions

| Code | Distribution     | Parameters                               |
| ---- | ---------------- | ---------------------------------------- |
| 0, 1 | No uncertainty   | --                                       |
| 2    | Lognormal        | loc (ln_mean), scale (ln_std)            |
| 3    | Normal           | loc (mean), scale (std)                  |
| 4    | Uniform          | minimum, maximum                         |
| 5    | Triangular       | minimum, maximum, loc (mode)             |
| 6    | Bernoulli        | loc (probability)                        |
| 7    | Discrete uniform | minimum, maximum                         |
| 8    | Weibull          | shape, loc, scale                        |
| 9    | Gamma            | shape (alpha), loc, scale                |
| 10   | Beta             | shape (alpha), shape2 (beta), loc, scale |
| 11   | Gumbel (GEV)     | shape (xi), loc, scale                   |
| 12   | Student-t        | shape (df), loc, scale                   |

## FAST / RBD-FAST GSA

| Variant   | Function            | Indices | Sample size formula | Notes                                    |
| --------- | ------------------- | ------- | ------------------- | ---------------------------------------- |
| Classic FAST | `method='fast'`  | S1 + ST | `(4M²+1)·k`         | Deterministic, first + total order, low sample cost |
| RBD-FAST  | `method='rbd_fast'` | S1      | `≥65·k`             | Random, first-order only, Latin hypercube sampling |

where `k` = number of uncertain parameters, `M` = interference factor (recommended default: 4).

**Comparison with Sobol**

| Method       | S1  | ST  | Min. recommended samples | When to use                       |
| ------------ | --- | --- | ------------------------ | --------------------------------- |
| Sobol        | ✓   | ✓   | `N·(k+2)`                | Standard full GSA                 |
| Classic FAST | ✓   | ✓   | `(4M²+1)·k`              | S1 + ST, tight sample budget      |
| RBD-FAST     | ✓   | ✗   | `65·k`                   | S1 only at lower cost than Sobol  |

Key API:

| Function                                        | Description                                        |
| ----------------------------------------------- | -------------------------------------------------- |
| `build_fast_problem(group)`                     | Build SALib problem dict from ActivityParameters   |
| `sanitize_fast_problem(problem)`                | Validate and repair invalid parameter bounds       |
| `generate_fast_samples(problem, M, method)`     | Generate FAST / RBD-FAST sample matrix             |
| `expected_fast_sample_rows(problem, M, method)` | Predict number of sample rows                      |
| `run_parallel_fast_from_samples(...)`           | Evaluate sample matrix with multiple workers       |
| `fast_indices_from_results(...)`                | Compute FAST / RBD-FAST sensitivity indices        |
| `fast_indices_to_dataframe(Si, problem)`        | Convert indices to a tidy DataFrame                |
| `run_full_fast_workflow(...)`                   | One-click workflow (build → sample → evaluate)     |

## Background Uncertainty (BG)

By default the GSA only perturbs the **foreground** ActivityParameters while the
background database (e.g. ecoinvent) stays at its deterministic values. Passing
`include_background=True` combines every parameter sample with one independent
draw of the **background** database uncertainty, so both uncertainty sources are
propagated together.

```python
# Monte Carlo with background uncertainty
mc_bg = run_parallel_monte_carlo(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    n_samples=1000,
    include_background=True,
)

# Sobol with background uncertainty
sobol_bg = run_parallel_sobol_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    sobol_problem=problem,
    sobol_samples=samples,
    include_background=True,
)

# FAST / RBD-FAST with background uncertainty
fast_bg = lca_fast.run_parallel_fast_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    fast_problem=fast_problem,
    fast_samples=fast_samples,
    include_background=True,
)
```

Implementation notes:

- Each worker builds an additional `MultiLCA(..., use_distributions=True)`
  instance and advances it by exactly one draw per sample.
- Parameterized foreground cells are masked out of the background draw
  (`_prepare_bg_data_objs`), otherwise the foreground perturbation would be
  counted twice (once by the background sampler, once by the delta overlay).
- The foreground perturbation is applied as a cell-wise **delta** on top of the
  background-sampled matrices (`compute_param_delta` /
  `apply_param_delta_to_matrices`), so only one LCI/LCIA round is needed per sample.
- Cost: background runs are slower and use more memory per worker, because an
  extra distribution-sampling MultiLCA must be kept in memory.

## Multiprocessing

Monte Carlo, Sobol, and FAST evaluations all use Python's `spawn` multiprocessing context to ensure Windows compatibility. Each worker process independently loads the Brightway project and databases, so memory usage scales with the number of workers. 2–4 workers are recommended.

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

- [Brightway2.5](https://brightway.dev/) -- LCA framework
- [SALib](https://salib.readthedocs.io/) -- Sensitivity Analysis Library
- [ecoinvent](https://ecoinvent.org/) -- Life cycle inventory database
