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
