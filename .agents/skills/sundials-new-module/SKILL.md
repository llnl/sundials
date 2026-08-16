---
name: new-module
description: Add a new SUNDIALS module (e.g., NVECTOR_*, SUNMATRIX_*, SUNLINSOL_*, SUNNONLINSOL_*, SUNMEMORY_*, or a new shared component) including source/header layout, CMake wiring, enable/disable options, exported CMake targets, installed component registration, and developer tests/examples/docs updates.
---

# Create a new SUNDIALS module

Follow existing module patterns; copy the closest neighbor (same category + similar dependencies) and adapt minimally.

Open these files as references when needed:

- `src/<category>/CMakeLists.txt` (how modules are conditionally added)
- `cmake/SundialsBuildOptionsPost.cmake` (module enable options + `SUNDIALS_BUILD_LIST`)
- `cmake/macros/SundialsAddLibrary.cmake` (how `sundials_add_library` installs headers and registers components)
- `doc/shared/sundials/Install.rst` (user-visible options / CMake targets)
- `test/unit_tests/<category>/` (how unit tests are wired)

## 1) Decide what “module” means

Pick one of these patterns; it determines naming, options, and wiring:

- **Native always-built module** (e.g., `SUNMATRIX_BAND`, `SUNLINSOL_SPGMR`): always enabled; no user option.
- **Optional module gated by a TPL/backend** (e.g., CUDA/HIP/Kokkos/Ginkgo): add `SUNDIALS_ENABLE_<CATEGORY>_<NAME>` option with `DEPENDS_ON ...`.
- **Header-only / interface-only module** (rare): uses an `INTERFACE` library; you must manually register it as an installed component.
- **Internal-only helper**: can be `OBJECT_LIB_ONLY` (no installed component/target intended for end users).

If the module is user-visible, prefer “optional module” (with a CMake option) over ad-hoc logic.

## 2) Choose names (directory, library, headers, CMake target)

Keep names consistent with existing modules in the same category:

- **Source dir**: `src/<category>/<name>/` (e.g., `src/sunlinsol/spgmr/`)
- **Public header(s)**: `include/<category>/<category>_<name>.<h|hpp>` (e.g., `include/sunlinsol/sunlinsol_spgmr.h`)
- **CMake/Library target**: `sundials_<category><name>` (e.g., `sundials_sunlinsolspgmr`)
- **Exported target for consumers**: `SUNDIALS::<category><name>` (auto-created by `sundials_add_library`)
- **Enable option** (if optional): `SUNDIALS_ENABLE_<CATEGORY>_<NAME>` (e.g., `SUNDIALS_ENABLE_SUNMATRIX_CUSPARSE`)

Use existing capitalization conventions:

- options use `SUNDIALS_ENABLE_<...>` with category in uppercase (NVECTOR/SUNMATRIX/SUNLINSOL/…)
- messages often say `Added <CATEGORY>_<NAME> module`

## 3) Add sources + headers

Create:

- `src/<category>/<name>/<implementation>.c` (and/or `.cpp` if needed)
- `include/<category>/<category>_<name>.h` (public API; installable)
- `src/<category>/<name>/CMakeLists.txt` (module build rules)

Avoid installing “private” headers. If a header is only for the module implementation, keep it under `src/<category>/<name>/` and do not list it in `HEADERS`.

## 4) Add the module CMakeLists.txt

For standard compiled modules, use `sundials_add_library` (preferred):

- `SOURCES ...`
- `HEADERS ${SUNDIALS_SOURCE_DIR}/include/<category>/<category>_<name>.h`
- `INCLUDE_SUBDIR <category>`
- `LINK_LIBRARIES PUBLIC sundials_core` (+ any TPL targets)
- `OUTPUT_NAME sundials_<category><name>`
- `VERSION`/`SOVERSION` using the category’s variables (see neighboring modules)

Reference examples:

- `src/nvector/serial/CMakeLists.txt`
- `src/sunmatrix/band/CMakeLists.txt`
- `src/sunlinsol/spgmr/CMakeLists.txt`

Notes:

- `sundials_add_library` automatically:
  - installs the listed headers to `${CMAKE_INSTALL_INCLUDEDIR}/<category>`
  - registers the module as an installed component for `SUNDIALSConfig.cmake`
  - creates the `SUNDIALS::<...>` alias for build-tree usage
- If you add Fortran 2003 wrappers like other modules do, gate them under `if(SUNDIALS_ENABLE_FORTRAN)` and mirror the existing `fmod_int${SUNDIALS_INDEX_SIZE}` pattern.

## 5) Wire the module into the build (add_subdirectory)

Edit `src/<category>/CMakeLists.txt` to include your module:

- Always-built: `add_subdirectory(<name>)`
- Optional: wrap with `if(SUNDIALS_ENABLE_<CATEGORY>_<NAME>) ... endif()`

For categories that are listed in `src/CMakeLists.txt` (e.g., `nvector`, `sunmatrix`, `sunlinsol`, …), you usually only need to change the category-level `CMakeLists.txt`.

## 6) Add/adjust the enable option and build-list registration

If the module is optional and should appear in `sundials_config.h` (and in build summaries), add an option in `cmake/SundialsBuildOptionsPost.cmake`:

- Add `sundials_option(SUNDIALS_ENABLE_<CATEGORY>_<NAME> BOOL "... module ..." <default> DEPENDS_ON ...)`
- Append to build list: `list(APPEND SUNDIALS_BUILD_LIST "SUNDIALS_ENABLE_<CATEGORY>_<NAME>")`

Pick appropriate dependencies:

- CUDA modules: `DEPENDS_ON SUNDIALS_ENABLE_CUDA` (and sometimes `CMAKE_CUDA_COMPILER`)
- MPI modules: `DEPENDS_ON SUNDIALS_ENABLE_MPI`
- TPL modules: depend on the corresponding `SUNDIALS_ENABLE_<TPL>` option

If the module is required and must not be disabled, follow the pattern:

```cmake
set(SUNDIALS_ENABLE_<CATEGORY>_<NAME> TRUE)
list(APPEND SUNDIALS_BUILD_LIST "SUNDIALS_ENABLE_<CATEGORY>_<NAME>")
```

## 7) Handle “installed components” for INTERFACE-only modules

If you implement a module as an `INTERFACE` library (no `sundials_add_library`), you must manually add it to `_SUNDIALS_INSTALLED_COMPONENTS` so that `find_package(SUNDIALS COMPONENTS ...)` can work.

Reference patterns:

- `src/sunmatrix/CMakeLists.txt` (Ginkgo/KokkosDense interface modules)
- `src/sunlinsol/CMakeLists.txt` (Ginkgo/KokkosDense interface modules)

## 8) Tests, examples, and docs (when user-visible)

Do the smallest set that proves correctness and prevents regressions:

- **Unit tests**: add under `test/unit_tests/<category>/...` and gate by your enable option if optional.
- **Examples**: add under `examples/<package>/<lang>_<backend>/...` only if it materially helps users.
- **Docs**:
  - If users must know about a new option/target, update `doc/shared/sundials/Install.rst`.
  - If it’s a user-visible feature, also update `doc/shared/RecentChanges.rst` and/or `CHANGELOG.md` per repo conventions.

## 9) Local validation commands

Configure/build (pick options that enable your module), then run focused tests first:

```bash
cmake -S . -B build-dev -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
cmake --build build-dev -j
ctest --test-dir build-dev --output-on-failure -R <your-module-regex>
```

If you added a new enable option, confirm it shows up in `build-dev/include/sundials/sundials_config.h`.

## 10) Language bindings

Only add language bindings when they are applicable:

- The module is user-visible (i.e., part of the public API users are expected to call), and
- The module’s category already has existing Fortran and/or Python bindings.

If Fortran bindings are applicable:

- Add the SWIG interface file in `swig/` and add it to `swig/Makefile` appropriately.
- Generate the Fortran interface code and commit the generated sources.

If Python bindings are applicable:

- Add the module as a neighbor to the other modules of the same type in `bindings/sundials4py`.
- Add it to the `generate.yaml` with the appropriate neighbors.
- Generate the Python binding code to commit, and always format it first.
