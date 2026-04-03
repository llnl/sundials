---
name: building
description: Build, install, and test SUNDIALS from source as an end user or as a SUNDIALS developer. Use when a request involves configuring CMake, selecting compilers/options (MPI, GPU backends, Fortran), building/installing, running CTest, enabling dev/unit tests, or troubleshooting common build/test issues.
---

# Build SUNDIALS (user then developer)

Prefer out-of-source CMake builds. Never build in the source tree.

If you need full option details or platform-specific notes, open `doc/shared/sundials/Install.rst` (search for “Configuration options”, “Build Type”, “Compilers”, and “Example Programs”).

For CI-like, multi-config testing use `test/test_driver.sh` (supports `--testtype pr|release|branch` and `--buildjobs/--testjobs`).

Python bindings (sundials4py) are driven by `pyproject.toml` and `scikit-build-core`:
`python -m pip install -e ".[dev]"` then `pytest`.

## User build (install + consume)

1) Configure (choose an install prefix):

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/install"
```

2) Build and install:

```bash
cmake --build build -j
cmake --install build
```

3) Use SUNDIALS from another CMake project:

- Set `SUNDIALS_DIR` to the SUNDIALS install tree (or to the installed CMake package dir) and use `find_package(SUNDIALS REQUIRED)`; link with exported targets like `SUNDIALS::cvode`. If you link via CMake targets, the required `libsundials_core` dependency is handled automatically.

### Common user options

- **Examples**: examples are usually enabled by default for C; enable others as needed:
  - `-DSUNDIALS_ENABLE_C_EXAMPLES=ON`
  - `-DSUNDIALS_ENABLE_CXX_EXAMPLES=ON`
  - `-DSUNDIALS_ENABLE_FORTRAN_EXAMPLES=ON` (requires Fortran)
  - `-DSUNDIALS_ENABLE_CUDA_EXAMPLES=ON` (requires CUDA enabled)
  - If working with older option names, note that `EXAMPLES_ENABLE_C`/`EXAMPLES_ENABLE_CXX`/etc. are deprecated in favor of the `SUNDIALS_ENABLE_*_EXAMPLES` options.
- **MPI**: `-DSUNDIALS_ENABLE_MPI=ON`
- **Fortran interfaces**: `-DSUNDIALS_ENABLE_FORTRAN=ON`
- **Build type**: `Debug` (slow, checks), `Release` (fast), `RelWithDebInfo` (good default)

## Developer build (iterate + test)

Use a separate build directory per configuration (e.g., `build-debug`, `build-gcc`, `build-mpi`).

### Configure for development

```bash
cmake -S . -B build-dev \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBUILD_TESTING=ON \
  -DSUNDIALS_TEST_ENABLE_DEV_TESTS=ON \
  -DSUNDIALS_TEST_ENABLE_UNIT_TESTS=ON
```

Optional tightening (CI-like):

- Warnings: `-DSUNDIALS_ENABLE_ALL_WARNINGS=ON`
- Warnings as errors: `-DCMAKE_COMPILE_WARNING_AS_ERROR=ON`
- Sanitizers (compiler support required): `-DSUNDIALS_ENABLE_ADDRESS_SANITIZER=ON` (and/or leak/undefined/thread variants)

### Build and run tests

```bash
cmake --build build-dev -j
ctest --test-dir build-dev --output-on-failure
```

Notes:

- If tests that compare against “answer files” fail due to platform differences, consider rerunning with a locally-generated answer directory using `SUNDIALS_TEST_ANSWER_DIR` (see `doc/superbuild/source/developers/testing/CTest.rst`).
- To focus on a subset, use CTest filters (e.g., `ctest --test-dir build-dev -R <regex>`).

## Troubleshooting checklist

- **CMake refuses in-source build**: delete any accidental `CMakeCache.txt` in the source tree; configure with `-B <builddir>`.
- **Changed options not taking effect**: start from a fresh build directory or clear cache (`rm -rf build-dev`).
- **Wrong compiler/MPI wrapper**: set `CC`, `CXX`, `FC` env vars before configuring, or pass `-DCMAKE_<LANG>_COMPILER=...`.
- **Link errors after enabling a backend/TPL**: confirm that the corresponding `SUNDIALS_ENABLE_*` option is ON and that required dependencies are discoverable (use `CMAKE_PREFIX_PATH`).
