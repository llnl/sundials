# Repository Guidelines

## Project Structure

- `src/` and `include/`: core SUNDIALS C/C++ sources and public headers, organized by module (e.g., `cvode/`, `ida/`, `nvector/`, `sunlinsol/`).
- `examples/`: buildable examples (C, C++, CUDA, etc.) wired into CMake.
- `test/`: test infrastructure, answer files, and `test/unit_tests/` (optionally includes GoogleTest-based tests).
- `cmake/`: CMake modules and option definitions (e.g., example/test toggles).
- `bindings/` and `swig/`: language bindings (notably `bindings/sundials4py`).
- `doc/`: user/developer documentation sources.

## Build, Test, and Development Commands

Preferred build is CMake out-of-source:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
  -DEXAMPLES_ENABLE_C=ON -DSUNDIALS_TEST_ENABLE_UNIT_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

For CI-like, multi-config testing use `test/test_driver.sh` (supports `--testtype pr|release|branch` and `--buildjobs/--testjobs`).

Python bindings (sundials4py) are driven by `pyproject.toml` and `scikit-build-core`:
`python -m pip install -e ".[dev]"` then `pytest`.

## Coding Style & Naming

- C/C++: follow `.clang-format` (2-space indent); format before pushing.
- CMake: keep listfile style consistent with `.cmake-format.py`.
- Fortran: formatted via `fprettify --indent 2`.
- Python: formatted with Black (configured in `pyproject.toml`).

Formatting helpers:
`./scripts/format.sh src include cmake test bindings`
`./scripts/spelling.sh` (uses `codespell` and may edit files).

## Testing Guidelines

- Primary harness: CTest (`ctest --test-dir build`).
- Unit tests: enable with `-DSUNDIALS_TEST_ENABLE_UNIT_TESTS=ON` (and optionally `-DSUNDIALS_TEST_ENABLE_GTEST=ON`).
- Keep new features paired with targeted tests; update answer files under `test/answers/` only when output changes are intended.

## Commits & Pull Requests

- Target the `develop` branch.
- Use clear, scoped subjects; common patterns in history include `CMake: ...`, `CI: ...`, `Docs: ...`, `Tools: ...`.
- All commits must include a DCO sign-off: use `git commit -s`.
- PRs should include: what/why, how it was tested, and updates to `CHANGELOG.md` and `doc/shared/RecentChanges.rst` when user-visible.

