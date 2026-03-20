# suntools (SUNDIALS)

This directory contains the `suntools` Python package which provides utilities for

- Structured parsing of logs produced by the SUNDIALS `SUNLogger` including filtering and plotting
- more to come...

- Structured parsing of logs produced by the SUNDIALS `SUNLogger` including filtering and plotting

## Install

From this directory:

```bash
python -m pip install -e .
```

Then import as:

```python
from suntools import logs
```

## CLI Log Parsing Examples

Filter a SUNLogger logfile (or stdin) down to high-level categories:

```bash
cat sun.log | suntools parse_logs --filter="integrator,nonlinear,linear"
```

Supported filter categories:

- `integrator` (step-attempt region output)
- `nonlinear` (nonlinear solver output, includes `begin/end-nonlinear-solve` and KINSOL info logs)
- `linear` (linear solver output, includes `begin/end-linear-solve` and KINSOL `kinLs*` scopes)

Filter using specific scopes:
 

Filter specific functions/regions:
 

Filter based on any keyword:


