# ROCm wheel for sundials4py

This wheel contains the HIP-enabled native bindings for the AMD GPU targets
supported by the pinned PyTorch ROCm wheel: `gfx900`, `gfx906`, `gfx908`,
`gfx90a`, `gfx942`, `gfx1030`, `gfx1100`, `gfx1101`, `gfx1102`, `gfx1200`, and
`gfx1201`. It is built from the full SUNDIALS checkout. It installs a separate
`sundials4py_rocm` package so it does not overwrite any file owned by the base
`sundials4py` wheel.

Build a wheel with ROCm 6.3 installed:

```bash
export ROCM_PATH=/opt/rocm-6.3.1
export HIP_PATH="$ROCM_PATH"
export CC="$ROCM_PATH/bin/amdclang"
export CXX="$ROCM_PATH/bin/amdclang++"
python -m pip wheel --no-deps ./bindings/sundials4py/rocm
```

Build and publish a matching base `sundials4py` wheel. Keep the versions in
the root and this directory's `pyproject.toml` files and their mutual
requirements synchronized. Use a new release version when publishing; an
existing base wheel built before the ROCm loader was added cannot select this
backend.

The index used for `pip install 'sundials4py[rocm]'` must also provide
`torch==2.7.0+rocm6.3`. PyPI does not host this PyTorch variant. A pip extra
cannot add a package index, so configure pip or use a mirror before installing.

Only the wheel is intended for publication. The standalone source distribution
for this directory does not contain the SUNDIALS source tree.
