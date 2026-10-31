# suntools (SUNDIALS)

This directory contains the `suntools` Python package which provides utilities for

- parsing and filtering SUNDIALS log files
- parsing table and CSV statistics output
- Runge-Kutta Butcher-table parsing and stability-region plotting under `suntools.rk`

## Install

suntools requires Python 3.12 or newer and installs NumPy and Matplotlib as dependencies.

From this directory:

```bash
python -m pip install -e .
```

Then import as:

```python
import suntools
```
