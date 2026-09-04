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
