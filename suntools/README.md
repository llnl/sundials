# suntools (SUNDIALS)

This directory contains the `suntools` Python package which provides utilities for

- parsing and filtering SUNDIALS log files
- parsing table and CSV statistics output

## Install

From this directory:

```bash
python -m pip install -e .
```

To enable the optional Ytopt tuning backend, install its extra:

```bash
python -m pip install -e ".[ytopt]"
```

Then import as:

```python
import suntools
```

## Tune

The ``suntools tune`` command can tune SUNDIALS ``SetOptions`` parameters.
Categorical choices may contain whitespace-separated values when an option
accepts multiple arguments. For example, ``arkode.table_names`` accepts an
implicit and an explicit table name. The following searches keep each method
class together by tuning the ARK324/ARK436/ARK548 L2SA family within that
class:

```bash
# Explicit-only ARK methods
suntools tune \
  --params arkode.table_names \
  'choice:ARKODE_DIRK_NONE ARKODE_ARK324L2SA_ERK_4_2_3,ARKODE_DIRK_NONE ARKODE_ARK436L2SA_ERK_6_3_4,ARKODE_DIRK_NONE ARKODE_ARK548L2SA_ERK_8_4_5' \
  -- ./arkstep_explicit

# Implicit-only DIRK methods
suntools tune \
  --params arkode.table_names \
  'choice:ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ERK_NONE,ARKODE_ARK436L2SA_DIRK_6_3_4 ARKODE_ERK_NONE,ARKODE_ARK548L2SA_DIRK_8_4_5 ARKODE_ERK_NONE' \
  -- ./arkstep_implicit

# Matching IMEX ARK methods
suntools tune \
  --params arkode.table_names \
  'choice:ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ARK324L2SA_ERK_4_2_3,ARKODE_ARK436L2SA_DIRK_6_3_4 ARKODE_ARK436L2SA_ERK_6_3_4,ARKODE_ARK548L2SA_DIRK_8_4_5 ARKODE_ARK548L2SA_ERK_8_4_5' \
  -- ./arkstep_imex
```

Each whitespace-separated choice is expanded into a separate command-line
value. The corresponding YAML ``parameters`` blocks are:

```yaml
# Explicit-only ARK methods
parameters:
  - name: arkode.table_names
    type: choice
    values:
      - "ARKODE_DIRK_NONE ARKODE_ARK324L2SA_ERK_4_2_3"
      - "ARKODE_DIRK_NONE ARKODE_ARK436L2SA_ERK_6_3_4"
      - "ARKODE_DIRK_NONE ARKODE_ARK548L2SA_ERK_8_4_5"
```

```yaml
# Implicit-only DIRK methods
parameters:
  - name: arkode.table_names
    type: choice
    values:
      - "ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ERK_NONE"
      - "ARKODE_ARK436L2SA_DIRK_6_3_4 ARKODE_ERK_NONE"
      - "ARKODE_ARK548L2SA_DIRK_8_4_5 ARKODE_ERK_NONE"
```

```yaml
# Matching IMEX ARK methods
parameters:
  - name: arkode.table_names
    type: choice
    values:
      - "ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ARK324L2SA_ERK_4_2_3"
      - "ARKODE_ARK436L2SA_DIRK_6_3_4 ARKODE_ARK436L2SA_ERK_6_3_4"
      - "ARKODE_ARK548L2SA_DIRK_8_4_5 ARKODE_ARK548L2SA_ERK_8_4_5"
```

Add one of these ``parameters`` blocks to a tune configuration with the
desired executable, search, and objective settings, then run it with:

```bash
suntools tune --config tune.yaml
```

An objective or constraint may use a list of regexes with one numeric capture
per regex. The extracted values are combined with ``aggregation: sum`` or
``aggregation: mean``. This can express composite metrics such as the sum of
several solver statistics.

Environment variables in ``executable.command`` are expanded before each
trial. Variables may come from the process environment or ``executable.env``;
for example, ``$builddir/examples/arkode/C_serial/ark_analytic`` can be paired
with ``env: {builddir: /path/to/build}``.

Set ``search.repetitions`` (or use ``--repetitions``) to run each sampled
configuration multiple times. Objective, constraint, and wall-time metrics
are averaged over the repetitions. Every repetition must complete successfully
and satisfy the constraint for the trial to be feasible.

```yaml
search:
  max_evals: 40
  repetitions: 3
```

Using more than one worker is not recommended when the objective is
``wall_time``. Concurrent trials compete for CPU, memory, and other resources,
so their measured runtimes are not directly comparable. Use ``workers: 1`` for
wall-clock tuning when reproducible timings matter.

Every tuning run also evaluates the executable with its default settings,
using the configured number of repetitions, before applying tune parameters.
The CLI reports this baseline, the best feasible trial, and the worst feasible
trial. They are recorded as ``baseline.json``, ``best.json``, and
``worst.json`` in the results directory; the baseline does not consume a
search evaluation.

### Constrained tuning

Add a ``constraint`` block to require an extracted metric to stay below an
upper bound while minimizing the objective. Trials that exceed the bound, or
do not produce the constraint metric, are recorded but are not eligible to be
the best result.

For example, if an executable prints ``error=2.5e-7``:

```bash
suntools tune \
  --params arkode.table_names \
  'choice:TABLE_A,TABLE_B' \
  --metric wall_time \
  --constraint-metric error \
  --constraint-regex 'error=([0-9.eE+-]+)' \
  --constraint-upper-bound 1e-6 \
  -- ./arkstep_app
```

The YAML equivalent is:

```yaml
objective:
  metric: wall_time
  direction: minimize

constraint:
  metric: error
  source: stdout
  regex: "error=([0-9.eE+-]+)"
  group: 1
  upper_bound: 1.0e-6
```

The constraint is inclusive: a trial is feasible when its metric is less than
or equal to ``upper_bound``. Feasibility is written to ``results.csv`` and
``trials.jsonl``; ``best.json`` contains only a feasible best trial.
