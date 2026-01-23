# SUNDIALS Changelog

## Changes to SUNDIALS in release X.Y.Z

### Major Features

### New Features and Enhancements

The functions `CVodeGetUserDataB` and `IDAGetUserDataB` were added to CVODES
and IDAS, respectively.

### Bug Fixes

Some installed Fortran example makefiles were not linking to `sundials_fcore_mod`
and `sundials_core` libraries as they should be. This is now fixed.

On the initial time step with a user-supplied initial step size, ARKODE and
CVODE(S) will now return `ARK_TOO_CLOSE` or `CV_TOO_CLOSE`, respectively,
when the requested output time is the same as the initial time (or within
numerical roundoff of the initial time). Before a `TOO_CLOSE` error would only
be returned when internally estimating the initial step size. In IDA(S), added a
`IDA_TOO_CLOSE` return value for when the initial and output time are too
close. Previously, IDA(S) would return `IDA_ILL_INPUT`.

Fixed a bug in ARKODE, CVODE(S), and IDA(S) where the linear solver counters
were not reinitialized until the next call to advance the system. As such,
non-zero linear solver statistics could be returned if retrieving or printing
linear solver counters between the initialization and the next call to advance
the system.

The interface to Ginkgo batched linear solvers has been updated to fix build
errors when using 64-bit index types. Note, only the batched dense matrix in
Ginkgo is currently compatible with 64-bit indexing (as of Ginkgo 1.10).

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug preventing a user supplied `SUNStepper_ResetCheckpointIndex` function from
being called was fixed.

The Kokkos N_Vector now properly handles unmanaged views. Previously, if a
Kokkos `N_Vector` was created from an unmanaged view, the view would become a
managed view and the data would be freed unexpectedly.

A bug was fixed in KINSOL where the information logging function would always be
called even when informational logging was disabled.

### Deprecation Notices

`SUNDIALSFileOpen` and `SUNDIALSFileClose` will be removed in the next major release. 
Use `SUNFileOpen` and `SUNFileClose` instead.

The `Convert` methods on the `sundials::kokkos:Vector`, `sundials::kokkos::DenseMatrix`,
`sundials::ginkgo::Matrix`, `sundials::ginkgo::BatchMatrix`, `sundials::kokkos::DenseLinearSolver`,
`sundials::ginkgo::LinearSolver`, and `sundials::ginkgo::BatchLinearSolver` classes have
been deprecated and will be removed in the next major release. The method `get`, should
be used instead.

The `CSC_MAT` and `CSR_MAT` macros defined in `sunmatrix_sparse.h` will be removed in
the next major release. Use `SUN_CSC_MAT` and `SUN_CSR_MAT` instead.

The `N_Vector_S` typedef to `N_Vector*` is deprecated and will be removed in the next major release.

## Changes to SUNDIALS in release 7.5.0

### Major Features

Added the `SUNDomEigEstimator` interface for estimating the dominant eigenvalue value
of a system. Two implementations are provided: Power Iteration and Arnoldi Iteration.
The latter method requires building with LAPACK support enabled.

Added the function `LSRKStepSetDomEigEstimator` in LSRKStep to attach a
`SUNDomEigEstimator`, when using Runge-Kutta-Chebyshev or Runge-Kutta-Legendre
methods, as an alternative to supplying a user-defined function to compute the dominant
eigenvalue.

Added `SetOptions` functions all SUNDIALS packages and the classes for
adaptivity controllers, dominant eigenvalue estimators, linear solvers, and
nonlinear solvers to support setting options with command line inputs.

### New Features and Enhancements

A new SUNLinearSolver, SUNLINEARSOLVER_GINKGOBATCH, and corresponding SUNMatrix,
SUNMATRIX_GINKGOBATCH, were added for solving block/batched linear systems with
the [Ginkgo linear solver library](https://ginkgo-project.github.io/). As a
result, Ginkgo 1.9.0 or newer is now required when enabling Ginkgo support.

The functions `KINSetMAA` and `KINSetOrthAA` have been updated to allow for
setting the Anderson acceleration depth and orthogonalization method after
`KINInit`. Additionally, `KINSetMAA` and `KINSetNumMaxIters` may now be called
in any order.

### Bug Fixes

Fixed a bug in how MRIStep interacts with an MRIHTol SUNAdaptController object
(the previous version essentially just reverted to a decoupled multirate
controller). Removed the upper limit on `inner_max_tolfac` in
`SUNAdaptController_SetParams_MRIHTol`.

The shared library version numbers for the oneMKL dense linear solver and
matrix as well as the PETSc SNES nonlinear solver have been corrected.

Fixed a CMake bug where the MRI H-Tol controller was not included in the ARKODE
Fortran module.

Fixed a bug in the CUDA and HIP implementations of `SUNMemoryHelper_CopyAsync`
where the execution stream is not extracted correctly from the helper when a
stream is not provided to `SUNMemoryHelper_CopyAsync`.

Fixed a bug in MRIStep where a segfault would occur when an MRI coupling table
is not explicitly set and an MRI integrator is nested inside another MRI
integrator.

Fixed a bug in MRIStep where MERK methods with unordered stage groups (MERK43
and MERK54) would include stage right-hand side vectors that had not been
computed yet in fast time scale forcing computations. These vectors were scaled
by zero, so in most cases the extraneous computations would not impact results.
However, in cases where these vectors contain `inf` or `nan`, this would lead to
erroneous forcing terms.

Fixed a bug in `ARKodeSetDefaults` with LSRKStep where the stored spectral
radius data was reset to zero, flags to update the dominant eigenvalue were
reset to true, and a flag indicating if an SSP is being used was reset to false.

Fixed a bug introduced in v7.3.0 in KINSOL when using Anderson acceleration and
solving a problem multiple times with the same KINSOL instance. In this use
case, the current Anderson acceleration depth from the initial solve was not
reinitialized on subsequent solves.

Fixed a logging bug in KINSOL where logging messages would not be output.

Fixed a bug in the `suntools.logs` Python module where the `get_history`
function, when given a `step_status` for filtering output from a multirate
method, would only extract values from the fast time scale if the slow time
scale step matched the given status filter. Fixed an additional bug in
`get_history` with MRI-GARK methods where values would not be extracted from a
fast time scale integration associated with an embedding.

## Changes to SUNDIALS in release 7.4.0

### New Features and Enhancements

`ARKodeSetCFLFraction` now allows `cfl_frac` to be greater than or equal to one.

Added an option to enable compensated summation of the time accumulator for all
of ARKODE. This was previously only an option for the SPRKStep module. The new
function to call to enable this is `ARKodeSetUseCompensatedSums`.

### Bug Fixes

Fixed segfaults in `CVodeAdjInit` and `IDAAdjInit` when called after adjoint
memory has been freed.

Fixed a CMake bug that would cause the Caliper compile test to fail at configure
time.

Fixed a bug in the CVODE/CVODES `CVodeSetEtaFixedStepBounds` function which
disallowed setting `eta_min_fx` or `eta_max_fx` to 1.

`SUNAdjointStepper_PrintAllStats` was reporting the wrong quantity for the
number of "recompute passes" and has been fixed.

### Deprecation Notices

The `SPRKStepSetUseCompensatedSums` function has been deprecated. Use the
`ARKodeSetUseCompensatedSums` function instead.

## Changes to SUNDIALS in release 7.3.0

### Major Features

A new discrete adjoint capability for explicit Runge-Kutta methods has been
added to the ARKODE ERKStep and ARKStep stepper modules. This is based on a
new set of shared classes, `SUNAdjointStepper` and `SUNAdjointCheckpointScheme`.
A new example demonstrating this capability can be found in
`examples/arkode/C_serial/ark_lotka_volterra_ASA.c`.

### New Features and Enhancements

#### ARKODE

The following changes have been made to the default ERK, DIRK, and ARK methods
in ARKODE to utilize more efficient methods:

| Type               | Old Default                                                      | New Default                                                        |
| ------------------ | ---------------------------------------------------------------- | ------------------------------------------------------------------ |
| 2nd Order Explicit | `ARKODE_HEUN_EULER_2_1_2`                                        | `ARKODE_RALSTON_3_1_2`                                             |
| 4th Order Explicit | `ARKODE_ZONNEVELD_5_3_4`                                         | `ARKODE_SOFRONIOU_SPALETTA_5_3_4`                                  |
| 5th Order Explicit | `ARKODE_CASH_KARP_6_4_5`                                         | `ARKODE_TSITOURAS_7_4_5`                                           |
| 6th Order Explicit | `ARKODE_VERNER_8_5_6`                                            | `ARKODE_VERNER_9_5_6`                                              |
| 8th Order Explicit | `ARKODE_FEHLBERG_13_7_8`                                         | `ARKODE_VERNER_13_7_8`                                             |
| 2nd Order Implicit | `ARKODE_SDIRK_2_1_2`                                             | `ARKODE_ARK2_DIRK_3_1_2`                                           |
| 3rd Order Implicit | `ARKODE_ARK324L2SA_DIRK_4_2_3`                                   | `ARKODE_ESDIRK325L2SA_5_2_3`                                       |
| 4th Order Implicit | `ARKODE_SDIRK_5_3_4`                                             | `ARKODE_ESDIRK436L2SA_6_3_4`                                       |
| 5th Order Implicit | `ARKODE_ARK548L2SA_DIRK_8_4_5`                                   | `ARKODE_ESDIRK547L2SA2_7_4_5`                                      |
| 4th Order ARK      | `ARKODE_ARK436L2SA_ERK_6_3_4` and `ARKODE_ARK436L2SA_DIRK_6_3_4` | `ARKODE_ARK437L2SA_ERK_7_3_4` and `ARKODE_ARK437L2SA_DIRK_7_3_4`   |
| 5th Order ARK      | `ARKODE_ARK548L2SA_ERK_8_4_5` and `ARKODE_ARK548L2SA_DIRK_8_4_5` | `ARKODE_ARK548L2SAb_ERK_8_4_5` and `ARKODE_ARK548L2SAb_DIRK_8_4_5` |


The old default methods can be loaded using the functions `ERKStepSetTableName`
or `ERKStepSetTableNum` with ERKStep and `ARKStepSetTableName` or
`ARKStepSetTableNum` with ARKStep and passing the desired method name string or
constant, respectively. For example, the following call can be used to load the
old default fourth order method with ERKStep:
```
/* Load the old 4th order ERK method using the table name */
ierr = ERKStepSetTableName(arkode_mem, "ARKODE_ZONNEVELD_5_3_4");
```
Similarly with ARKStep, the following calls can be used for ERK, DIRK, or ARK
methods, respectively:
```
/* Load the old 4th order ERK method by name */
ierr = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE",
                           "ARKODE_ZONNEVELD_5_3_4");

/* Load the old 4th order DIRK method by name */
ierr = ARKStepSetTableName(arkode_mem, "ARKODE_SDIRK_5_3_4",
                           "ARKODE_ERK_NONE");

/* Load the old 4th order ARK method by name */
ierr = ARKStepSetTableName(arkode_mem, "ARKODE_ARK436L2SA_DIRK_6_3_4",
                           "ARKODE_ARK436L2SA_ERK_6_3_4");
```

Additionally, the following changes have been made to the default time step
adaptivity parameters in ARKODE:

| Parameter             | Old Default           | New Default |
|-----------------------|-----------------------|-------------|
| Controller            | PID (PI for ERKStep)  | I           |
| Safety Factor         | 0.96                  | 0.9         |
| Bias                  | 1.5 (1.2 for ERKStep) | 1.0         |
| Fixed Step Bounds     | [1.0, 1.5]            | [1.0, 1.0]  |
| Adaptivity Adjustment | -1                    | 0           |

The following calls can be used to restore the old defaults for ERKStep:
```
SUNAdaptController controller = SUNAdaptController_Soderlind(ctx);
SUNAdaptController_SetParams_PI(controller, 0.8, -0.31);
ARKodeSetAdaptController(arkode_mem, controller);
SUNAdaptController_SetErrorBias(controller, 1.2);
ARKodeSetSafetyFactor(arkode_mem, 0.96);
ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
ARKodeSetAdaptivityAdjustment(arkode_mem, -1);
```
The following calls can be used to restore the old defaults for other ARKODE
integrators:
```
SUNAdaptController controller = SUNAdaptController_PID(ctx);
ARKodeSetAdaptController(arkode_mem, controller);
SUNAdaptController_SetErrorBias(controller, 1.5);
ARKodeSetSafetyFactor(arkode_mem, 0.96);
ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
ARKodeSetAdaptivityAdjustment(arkode_mem, -1);
```
In both cases above, destroy the controller at the end of the run with
`SUNAdaptController_Destroy(controller);`.

The Soderlind time step adaptivity controller now starts with an I controller
until there is sufficient history of past time steps and errors.

Added `ARKodeSetAdaptControllerByName` to set a time step adaptivity controller
with a string. There are also four new controllers: `SUNAdaptController_H0211`,
`SUNAdaptController_H0321`, `SUNAdaptController_H211`, and
`SUNAdaptController_H312`.

Added the `ARKODE_RALSTON_3_1_2` and `ARKODE_TSITOURAS_7_4_5` explicit
Runge-Kutta Butcher tables.

Improved the precision of the coefficients for `ARKODE_ARK324L2SA_ERK_4_2_3`,
`ARKODE_VERNER_9_5_6`, `ARKODE_VERNER_10_6_7`, `ARKODE_VERNER_13_7_8`,
`ARKODE_ARK324L2SA_DIRK_4_2_3`, and `ARKODE_ESDIRK324L2SA_4_2_3`.

#### CVODE / CVODES

Added support for resizing CVODE and CVODES when solving initial value problems
where the number of equations and unknowns changes over time. Resizing requires
a user supplied history of solution and right-hand side values at the new
problem size, see `CVodeResizeHistory` for more information.

#### KINSOL

Added support in KINSOL for setting user-supplied functions to compute the
damping factor and, when using Anderson acceleration, the depth in fixed-point
or Picard iterations. See `KINSetDampingFn` and `KINSetDepthFn`, respectively,
for more information.

#### SUNDIALS Types

A new type, `suncountertype`, was added for the integer type used for counter
variables. It is currently an alias for `long int`.

### Bug Fixes

#### ARKODE

Fixed bug in `ARKodeResize` which caused it return an error for MRI methods.

Removed error floors from the `SUNAdaptController` implementations which could
unnecessarily limit the time size growth, particularly after the first step.

Fixed bug in `ARKodeSetFixedStep` where it could return `ARK_SUCCESS` despite
an error occurring.

Fixed bug in the ARKODE SPRKStep `SPRKStepReInit` function and `ARKodeReset` function
with SPRKStep that could cause a segmentation fault when compensated summation is not
used.

#### KINSOL

Fixed a bug in KINSOL where an incorrect damping parameter is applied on the
initial iteration with Anderson acceleration unless `KINSetDamping` and
`KINSetDampingAA` are both called with the same value when enabling damping.

Fixed a bug in KINSOL where errors that occurred when computing Anderson
acceleration were not captured.

Added missing return values to `KINGetReturnFlagName`.

#### CMake

Fixed the behavior of `SUNDIALS_ENABLE_ERROR_CHECKS` so additional runtime error
checks are disabled by default with all release build types. Previously,
`MinSizeRel` builds enabled additional error checking by default.

### Deprecation Notices

All work space functions, e.g., `CVodeGetWorkSpace` and `ARKodeGetLinWorkSpace`,
have been deprecated and will be removed in version 8.0.0.

## Changes to SUNDIALS in release 7.2.1

### New Features and Enhancements

Unit tests were separated from examples. To that end, the following directories
were moved out of the `examples/` directory to the `test/unit_tests` directory:
`nvector`, `sunmatrix`, `sunlinsol`, and `sunnonlinsol`.

### Bug Fixes

Fixed a bug in ARKStep where an extra right-hand side evaluation would occur
each time step when enabling the ``ARKodeSetAutonomous`` option and using an
IMEX method where the DIRK table has an implicit first stage and is not stiffly
accurate.

## Changes to SUNDIALS in release 7.2.0

### Major Features

Added a time-stepping module to ARKODE for low storage Runge-Kutta methods,
LSRKStep.  This currently supports five explicit low-storage methods: the
second-order Runge-Kutta-Chebyshev and Runge-Kutta-Legendre methods, and the
second- through fourth-order optimal strong stability preserving Runge-Kutta
methods. All methods include embeddings for temporal adaptivity.

Added an operator splitting module, SplittingStep, and forcing method module,
ForcingStep, to ARKODE. These modules support a broad range of operator-split
time integration methods for multiphysics applications.

Added support for multirate time step adaptivity controllers, based on the
recently introduced `SUNAdaptController` base class, to ARKODE's MRIStep module.
As a part of this, we added embeddings for existing MRI-GARK methods, as well as
support for embedded MERK and IMEX-MRI-SR methods. Added new default MRI methods
for temporally adaptive versus fixed-step runs.

### New Features and Enhancements

#### Logging

The information level logging output in ARKODE, CVODE(S), and IDA(S) has been
updated to be more uniform across the packages and a new `tools` directory has
been added with a Python module, `suntools`, containing utilities for parsing
logging output. The Python utilities for parsing CSV output have been relocated
from the `scripts` directory to the Python module.

#### SUNStepper

Added the `SUNStepper` base class to represent a generic solution procedure for
IVPs. This is used by the SplittingStep and ForcingStep modules of ARKODE. A
SUNStepper can be created from an ARKODE memory block with the new function
`ARKodeCreateSUNStepper`. To enable interoperability with `MRIStepInnerStepper`,
the function `MRIStepInnerStepper_CreateFromSUNStepper` was added.

#### ARKODE

Added functionality to ARKODE to accumulate a temporal error estimate over
multiple time steps. See the routines `ARKodeSetAccumulatedErrorType`,
`ARKodeResetAccumulatedError`, and `ARKodeGetAccumulatedError` for details.

Added the `ARKodeSetStepDirection` and `ARKodeGetStepDirection` functions to
change and query the direction of integration.

Added the function `MRIStepGetNumInnerStepperFails` to retrieve the number of
recoverable failures reported by the MRIStepInnerStepper.

Added a utility routine to wrap any valid ARKODE integrator for use as an
MRIStep inner stepper object, `ARKodeCreateMRIStepInnerStepper`.

The following DIRK schemes now have coefficients accurate to quad precision:

* `ARKODE_BILLINGTON_3_3_2`
* `ARKODE_KVAERNO_4_2_3`
* `ARKODE_CASH_5_2_4`
* `ARKODE_CASH_5_3_4`
* `ARKODE_KVAERNO_5_3_4`
* `ARKODE_KVAERNO_7_4_5`

#### CMake

The default value of `CMAKE_CUDA_ARCHITECTURES` is no longer set to `70` and is
now determined automatically by CMake. The previous default was only valid for
Volta GPUs while the automatically selected value will vary across compilers and
compiler versions. As such, users are encouraged to override this value with the
architecture for their system.

The build system has been updated to utilize the CMake LAPACK imported target
which should ease building SUNDIALS with LAPACK libraries that require setting
specific linker flags e.g., MKL.

#### Third Party Libraries

The Trilinos Teptra NVector interface has been updated to utilize CMake imported
targets added in Trilinos 14 to improve support for different Kokkos backends
with Trilinos. As such, Trilinos 14 or newer is required and the
`Trilinos_INTERFACE_*` CMake options have been removed.

Example programs using *hypre* have been updated to support v2.20 and newer.

### Bug Fixes

#### CMake

Fixed a CMake bug regarding usage of missing "print_warning" macro that was only
triggered when the deprecated `CUDA_ARCH` option was used.

Fixed a CMake configuration issue related to aliasing an `ALIAS` target when
using `ENABLE_KLU=ON` in combination with a static-only build of SuiteSparse.

Fixed a CMake issue which caused third-party CMake variables to be unset.  Users
may see more options in the CMake GUI now as a result of the fix.  See details
in GitHub Issue [#538](https://github.com/LLNL/sundials/issues/538).

#### NVector

Fixed a build failure with the SYCL NVector when using Intel oneAPI 2025.0
compilers. See GitHub Issue [#596](https://github.com/LLNL/sundials/issues/596).

Fixed compilation errors when building the Trilinos Teptra NVector with CUDA
support.

#### SUNMatrix

Fixed a [bug](https://github.com/LLNL/sundials/issues/581) in the sparse matrix
implementation of `SUNMatScaleAddI` which caused out of bounds writes unless
`indexvals` were in ascending order for each row/column.

#### SUNLinearSolver

Fixed a bug in the SPTFQMR linear solver where recoverable preconditioner errors
were reported as unrecoverable.

#### ARKODE

Fixed `ARKodeResize` not using the default `hscale` when an argument of `0` was
provided.

Fixed a memory leak that could occur if ``ARKodeSetDefaults`` is called
repeatedly.

Fixed the loading of ARKStep's default first order explicit method.

Fixed loading the default IMEX-MRI method if `ARKodeSetOrder` is used to specify
a third or fourth order method. Previously, the default second order method was
loaded in both cases.

Fixed potential memory leaks and out of bounds array accesses that could occur
in the ARKODE Lagrange interpolation module when changing the method order or
polynomial degree after re-initializing an integrator.

Fixed a bug in ARKODE when enabling rootfinding with fixed step sizes and the
initial value of the rootfinding function is zero. In this case, uninitialized
right-hand side data was used to compute a state value near the initial
condition to determine if any rootfinding functions are initially active.

Fixed a bug in MRIStep where the data supplied to the Hermite interpolation
module did not include contributions from the fast right-hand side
function. With this fix, users will see one additional fast right-hand side
function evaluation per slow step with the Hermite interpolation option.

Fixed a bug in SPRKStep when using compensated summations where the error vector
was not initialized to zero.

#### CVODE(S)

Fixed a bug where `CVodeSetProjFailEta` would ignore the `eta` parameter.

#### Fortran Interfaces

Fixed a bug in the 32-bit ``sunindextype`` Fortran interfaces to
``N_VGetSubvectorArrayPointer_ManyVector``,
``N_VGetSubvectorArrayPointer_MPIManyVector``, ``SUNBandMatrix_Column`` and
``SUNDenseMatrix_Column`` where 64-bit ``sunindextype`` interface functions were
used.

### Deprecation Notices

Deprecated the ARKStep-specific utility routine for wrapping an ARKStep instance
as an MRIStep inner stepper object, `ARKStepCreateMRIStepInnerStepper`. Use
`ARKodeCreateMRIStepInnerStepper` instead.

The ARKODE stepper specific functions to retrieve the number of right-hand side
function evaluations have been deprecated. Use `ARKodeGetNumRhsEvals` instead.

## Changes to SUNDIALS in release 7.1.1

### Bug Fixes

Fixed a [bug](https://github.com/LLNL/sundials/pull/523) in v7.1.0 with the SYCL
N_Vector `N_VSpace` function.

## Changes to SUNDIALS in release 7.1.0

### Major Features

Created shared user interface functions for ARKODE to allow more uniform control
over time-stepping algorithms, improved extensibility, and simplified code
maintenance. The corresponding stepper-specific user-callable functions are now
deprecated and will be removed in a future major release.

Added CMake infrastructure that enables externally maintained addons/plugins to
be *optionally* built with SUNDIALS. See the [Contributing
Guide](./CONTRIBUTING.md) for more details.

### New Features and Enhancements

Added support for Kokkos Kernels v4.

Added the following Runge-Kutta Butcher tables

* `ARKODE_FORWARD_EULER_1_1`
* `ARKODE_RALSTON_EULER_2_1_2`
* `ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2`
* `ARKODE_BACKWARD_EULER_1_1`
* `ARKODE_IMPLICIT_MIDPOINT_1_2`
* `ARKODE_IMPLICIT_TRAPEZOIDAL_2_2`

Added the following MRI coupling tables

* `ARKODE_MRI_GARK_FORWARD_EULER`
* `ARKODE_MRI_GARK_RALSTON2`
* `ARKODE_MRI_GARK_RALSTON3`
* `ARKODE_MRI_GARK_BACKWARD_EULER`
* `ARKODE_MRI_GARK_IMPLICIT_MIDPOINT`
* `ARKODE_IMEX_MRI_GARK_EULER`
* `ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL`
* `ARKODE_IMEX_MRI_GARK_MIDPOINT`

Added `ARKodeButcherTable_ERKIDToName` and `ARKodeButcherTable_DIRKIDToName` to
convert a Butcher table ID to a string representation.

Added the function ``ARKodeSetAutonomous`` in ARKODE to indicate that the
implicit right-hand side function does not explicitly depend on time. When using
the trivial predictor, an autonomous problem may reuse implicit function
evaluations across stage solves to reduce the total number of function
evaluations.

Users may now disable interpolated output in ARKODE by passing `ARK_INTERP_NONE`
to `ARKodeSetInterpolantType`. When interpolation is disabled, rootfinding is
not supported, implicit methods must use the trivial predictor (the default
option), and interpolation at stop times cannot be used (interpolating at stop
times is disabled by default). With interpolation disabled, calling
`ARKodeEvolve` in `ARK_NORMAL` mode will return at or past the requested output
time (setting a stop time may still be used to halt the integrator at a specific
time). Disabling interpolation will reduce the memory footprint of an integrator
by two or more state vectors (depending on the interpolant type and degree)
which can be beneficial when interpolation is not needed e.g., when integrating
to a final time without output in between or using an explicit fast time scale
integrator with an MRI method.

Added "Resize" capability to ARKODE's SPRKStep time-stepping module.

Enabled the Fortran interfaces to build with 32-bit `sunindextype`.

### Bug Fixes

Updated the CMake variable `HIP_PLATFORM` default to `amd` as the previous
default, `hcc`, is no longer recognized in ROCm 5.7.0 or newer. The new default
is also valid in older version of ROCm (at least back to version 4.3.1).

Renamed the DPCPP value for the `SUNDIALS_GINKGO_BACKENDS` CMake option to `SYCL`
to match Ginkgo's updated naming convention.

Changed the CMake version compatibility mode for SUNDIALS to `AnyNewerVersion`
instead of `SameMajorVersion`. This fixes the issue seen
[here](https://github.com/AMReX-Codes/amrex/pull/3835).

Fixed a CMake bug that caused an MPI linking error for our C++ examples in some
instances. Fixes [GitHub Issue #464](https://github.com/LLNL/sundials/issues/464).

Fixed the runtime library installation path for windows systems. This fix
changes the default library installation path from
`CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR` to
`CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_BINDIR`.

Fixed conflicting `.lib` files between shared and static libs when using `MSVC`
on Windows.

Fixed invalid `SUNDIALS_EXPORT` generated macro when building both shared and
static libs.

Fixed a bug in some Fortran examples where `c_null_ptr` was passed as an
argument to a function pointer instead of `c_null_funptr`. This caused
compilation issues with the Cray Fortran compiler.

Fixed a bug in the HIP execution policies where `WARP_SIZE` would not be set
with ROCm 6.0.0 or newer.

Fixed a bug that caused error messages to be cut off in some cases. Fixes
[GitHub Issue #461](https://github.com/LLNL/sundials/issues/461).

Fixed a memory leak when an error handler was added to a `SUNContext`. Fixes
[GitHub Issue #466](https://github.com/LLNL/sundials/issues/466).

Fixed a bug where `MRIStepEvolve` would not handle a recoverable error produced
from evolving the inner stepper.

Added missing `SetRootDirection` and `SetNoInactiveRootWarn` functions to
ARKODE's SPRKStep time-stepping module.

Fixed a bug in `ARKodeSPRKTable_Create` where the coefficient arrays were not
allocated.

Fix bug on LLP64 platforms (like Windows 64-bit) where `KLU_INDEXTYPE` could be
32 bits wide even if `SUNDIALS_INT64_T` is defined.

Check if size of `SuiteSparse_long` is 8 if the size of `sunindextype` is 8
when using KLU.

Fixed several build errors with the Fortran interfaces on Windows systems.

### Deprecation Notices

Numerous ARKODE stepper-specific functions are now deprecated in favor of
ARKODE-wide functions.

Deprecated the `ARKStepSetOptimalParams` function. Since this function does not have an
ARKODE-wide equivalent, instructions have been added to the user guide for how
to retain the current functionality using other user-callable functions.

The unsupported implementations of `N_VGetArrayPointer` and `N_VSetArrayPointer`
for the *hypre* and PETSc vectors are now deprecated. Users should access the
underlying wrapped external library vector objects instead with
`N_VGetVector_ParHyp` and `N_VGetVector_Petsc`, respectively.

## Changes to SUNDIALS in release v7.0.0

### Major Feature

SUNDIALS now has more robust and uniform error handling. Non-release builds will
be built with additional error checking by default. See the
[Error Checking](https://sundials.readthedocs.io/en/latest/sundials/Errors_link.html)
section in the user guide for details.

### Breaking Changes

#### Minimum C Standard

SUNDIALS now requires using a compiler that supports a subset of the C99
standard. Note with the Microsoft C/C++ compiler the subset of C99 features
utilized by SUNDIALS are available starting with [Visual Studio 2015](https://learn.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance?view=msvc-170#c-standard-library-features-1).

#### Minimum CMake Version

CMake 3.18 or newer is now required when building SUNDIALS.

#### Deprecated Types and Functions Removed

The previously deprecated types `realtype` and `booleantype` were removed from
`sundials_types.h` and replaced with `sunrealtype` and `sunbooleantype`. The
deprecated names for these types can be used by including the header file
`sundials_types_deprecated.h` but will be fully removed in the next major
release. Functions, types and header files that were previously deprecated have
also been removed.

#### Error Handling Changes

With the addition of the new error handling capability, the `*SetErrHandlerFn`
and `*SetErrFile` functions in CVODE(S), IDA(S), ARKODE, and KINSOL have been
removed. Users of these functions can use the functions
`SUNContext_PushErrHandler`, and `SUNLogger_SetErrorFilename` instead. For
further details see the
[Error Checking](https://sundials.readthedocs.io/en/latest/sundials/Errors_link.html)
and
[Logging](https://sundials.readthedocs.io/en/latest/sundials/Logging_link.html)
sections in the documentation.

In addition the following names/symbols were replaced by `SUN_ERR_*` codes:

| Removed                        | Replaced with `SUNErrCode`        |
|:-------------------------------|:----------------------------------|
| `SUNLS_SUCCESS`                | `SUN_SUCCESS`                     |
| `SUNLS_UNRECOV_FAILURE`        | no replacement (value was unused) |
| `SUNLS_MEM_NULL`               | `SUN_ERR_ARG_CORRUPT`             |
| `SUNLS_ILL_INPUT`              | `SUN_ERR_ARG_*`                   |
| `SUNLS_MEM_FAIL`               | `SUN_ERR_MEM_FAIL`                |
| `SUNLS_PACKAGE_FAIL_UNREC`     | `SUN_ERR_EXT_FAIL`                |
| `SUNLS_VECTOROP_ERR`           | `SUN_ERR_OP_FAIL`                 |
| `SUN_NLS_SUCCESS`              | `SUN_SUCCESS`                     |
| `SUN_NLS_MEM_NULL`             | `SUN_ERR_ARG_CORRUPT`             |
| `SUN_NLS_MEM_FAIL`             | `SUN_ERR_MEM_FAIL`                |
| `SUN_NLS_ILL_INPUT`            | `SUN_ERR_ARG_*`                   |
| `SUN_NLS_VECTOROP_ERR`         | `SUN_ERR_OP_FAIL`                 |
| `SUN_NLS_EXT_FAIL`             | `SUN_ERR_EXT_FAIL`                |
| `SUNMAT_SUCCESS`               | `SUN_SUCCESS`                     |
| `SUNMAT_ILL_INPUT`             | `SUN_ERR_ARG_*`                   |
| `SUNMAT_MEM_FAIL`              | `SUN_ERR_MEM_FAIL`                |
| `SUNMAT_OPERATION_FAIL`        | `SUN_ERR_OP_FAIL`                 |
| `SUNMAT_MATVEC_SETUP_REQUIRED` | `SUN_ERR_OP_FAIL`                 |

The following functions have had their signature updated to ensure they can
leverage the new SUNDIALS error handling capabilities.

```c
// From sundials_futils.h
SUNDIALSFileOpen
SUNDIALSFileClose

// From sundials_memory.h
SUNMemoryNewEmpty
SUNMemoryHelper_Alias
SUNMemoryHelper_Wrap

// From sundials_nvector.h
N_VNewVectorArray
```

#### SUNComm Type Added

We have replaced the use of a type-erased (i.e., `void*`) pointer to a
communicator in place of `MPI_Comm` throughout the SUNDIALS API with a
`SUNComm`, which is just a typedef to an `int` in builds without MPI
and a typedef to a `MPI_Comm` in builds with MPI. As a result:

* When MPI is enabled, all SUNDIALS libraries will include MPI symbols and
  applications will need to include the path for MPI headers and link against
  the corresponding MPI library.

* All users will need to update their codes because the call to
  `SUNContext_Create` now takes a `SUNComm` instead
  of type-erased pointer to a communicator. For non-MPI codes,
  pass `SUN_COMM_NULL` to the `comm` argument instead of
  `NULL`. For MPI codes, pass the `MPI_Comm` directly.

* The same change must be made for calls to
  `SUNLogger_Create` or `SUNProfiler_Create`.

* Some users will need to update their calls to `N_VGetCommunicator`, and
  update any custom `N_Vector` implementations that provide
  `N_VGetCommunicator`, since it now returns a `SUNComm`.

The change away from type-erased pointers for `SUNComm` fixes problems like the
one described in [GitHub Issue #275](https://github.com/LLNL/sundials/issues/275).

The SUNLogger is now always MPI-aware if MPI is enabled in SUNDIALS and the
`SUNDIALS_LOGGING_ENABLE_MPI` CMake option and macro definition were removed
accordingly.

#### SUNDIALS Core Library

Users now need to link to `sundials_core` in addition to the libraries already
linked to. This will be picked up automatically in projects that use the
SUNDIALS CMake target. The library `sundials_generic` has been superseded by
`sundials_core` and is no longer available. This fixes some duplicate symbol
errors on Windows when linking to multiple SUNDIALS libraries.

#### Fortran Interface Modules Streamlined

We have streamlined the Fortran modules that need to be included by users by combining
the SUNDIALS core into one Fortran module, `fsundials_core_mod`. Modules for
implementations of the core APIs still exist (e.g., for the Dense linear solver there
is `fsunlinsol_dense_mod`) as do the modules for the SUNDIALS packages (e.g., `fcvode_mod`).
The following modules are the ones that have been consolidated into `fsundials_core_mod`:

```
fsundials_adaptcontroller_mod
fsundials_context_mod
fsundials_futils_mod
fsundials_linearsolver_mod
fsundials_logger_mod
fsundials_matrix_mod
fsundials_nonlinearsolver_mod
fsundials_nvector_mod
fsundials_profiler_mod
fsundials_types_mod
```

### Minor Changes

The `CMAKE_BUILD_TYPE` defaults to `RelWithDebInfo` mode now i.e., SUNDIALS
will be built with optimizations and debugging symbols enabled by default.
Previously the build type was unset by default so no optimization or debugging
flags were set.

The advanced CMake options to override the inferred LAPACK name-mangling scheme
have been updated from `SUNDIALS_F77_FUNC_CASE` and
`SUNDIALS_F77_FUNC_UNDERSCORES` to `SUNDIALS_LAPACK_CASE` and
`SUNDIALS_LAPACK_UNDERSCORES`, respectively.

As a subset of C99 is now required the CMake option `USE_GENERIC_MATH` as been
removed.

The C++ convenience classes (e.g., `sundials::Context`) have been moved to
from SUNDIALS `.h` headers to corresponding `.hpp` headers (e.g.,
`sundials/sundials_context.hpp`) so C++ codes do not need to compile with
C++14 support when using the C API.

Converted most previous Fortran 77 and 90 examples to use SUNDIALS' Fortran 2003
interface.

### Bug Fixes

Fixed [#329](https://github.com/LLNL/sundials/issues/329) so that C++20 aggregate initialization can be used.

Fixed integer overflow in the internal SUNDIALS hashmap. This resolves
[#409](https://github.com/LLNL/sundials/issues/409) and
[#249](https://github.com/LLNL/sundials/issues/249).

### Deprecation Notice

The functions in `sundials_math.h` will be deprecated in the next release.

```c
  sunrealtype SUNRpowerI(sunrealtype base, int exponent);
  sunrealtype SUNRpowerR(sunrealtype base, sunrealtype exponent);
  sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);
  sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b, sunrealtype tol);
  sunrealtype SUNStrToReal(const char* str);
```

Additionally, the following header files (and everything in them) will be
deprecated -- users who rely on these are recommended to transition to the
corresponding `SUNMatrix` and `SUNLinearSolver` modules:

```c
sundials_direct.h
sundials_dense.h
sundials_band.h
```

## Changes to SUNDIALS in release 6.7.0

### Major Feature

Added the `SUNAdaptController` base class, ported ARKODE's internal
implementations of time step controllers into implementations of this class,
and updated ARKODE to use these objects instead of its own implementations.
Added `ARKStepSetAdaptController` and `ERKStepSetAdaptController` routines
so that users can modify controller parameters, or even provide custom
implementations.

### New Features

Improved computational complexity of `SUNMatScaleAddI_Sparse` from `O(M*N)` to
`O(NNZ)`.

Added Fortran support for the LAPACK dense `SUNLinearSolver` implementation.

Added the routines `ARKStepSetAdaptivityAdjustment` and
`ERKStepSetAdaptivityAdjustment`, that allow users to adjust the
value for the method order supplied to the temporal adaptivity controllers.
The ARKODE default for this adjustment has been -1 since its initial
release, but for some applications a value of 0 is more appropriate.
Users who notice that their simulations encounter a large number of
temporal error test failures may want to experiment with adjusting this value.

Added the third order ERK method `ARKODE_SHU_OSHER_3_2_3`, the fourth order
ERK method `ARKODE_SOFRONIOU_SPALETTA_5_3_4`, the sixth order ERK method
`ARKODE_VERNER_9_5_6`, the seventh order ERK method `ARKODE_VERNER_10_6_7`,
the eighth order ERK method `ARKODE_VERNER_13_7_8`, and the ninth order ERK
method `ARKODE_VERNER_16_8_9`.

ARKStep, ERKStep, MRIStep, and SPRKStep were updated to remove a potentially
unnecessary right-hand side evaluation at the end of an integration. ARKStep was
additionally updated to remove extra right-hand side evaluations when using an
explicit method or an implicit method with an explicit first stage.

The `MRIStepInnerStepper` class in MRIStep was updated to make supplying an
`MRIStepInnerFullRhsFn` optional.

### Bug Fixes

Changed the `SUNProfiler` so that it does not rely on `MPI_WTime` in any case.
This fixes [GitHub Issue #312](https://github.com/LLNL/sundials/issues/312).

Fixed scaling bug in `SUNMatScaleAddI_Sparse` for non-square matrices.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where ARKODE,
CVODE, CVODES, IDA, and IDAS would return at the stop time rather than the
requested output time if the stop time was reached in the same step in which the
output time was passed.

Fixed a bug in ERKStep where methods with `c[s-1] = 1` but `a[s-1,j] != b[j]`
were incorrectly treated as having the first same as last (FSAL) property.

Fixed a bug in ARKODE where `ARKStepSetInterpolateStopTime` would return an
interpolated solution at the stop time in some cases when interpolation was
disabled.

Fixed a bug in `ARKStepSetTableNum` wherein it did not recognize
`ARKODE_ARK2_ERK_3_1_2` and `ARKODE_ARK2_DIRK_3_1_2` as a valid additive
Runge-Kutta Butcher table pair.

Fixed a bug in `MRIStepCoupling_Write` where explicit coupling tables were not
written to the output file pointer.

Fixed missing soversions in some `SUNLinearSolver` and `SUNNonlinearSolver`
CMake targets.

Renamed some internal types in CVODES and IDAS to allow both packages to be
built together in the same binary.

## Changes to SUNDIALS in release 6.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK
installation of CUDA and fixed the targets used for rocBLAS and rocSPARSE.

## Changes to SUNDIALS in release 6.6.1

### New Features

Updated the Tpetra NVector interface to support Trilinos 14.

### Bug Fixes

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug in ARKODE, CVODE, CVODES, IDA, and IDAS where the stop time may not
be cleared when using normal mode if the requested output time is the same as
the stop time. Additionally, with ARKODE, CVODE, and CVODES this fix removes an
unnecessary interpolation of the solution at the stop time that could occur in
this case.

## Changes to SUNDIALS in release 6.6.0

### Major Features

A new time-stepping module, `SPRKStep`, was added to ARKODE. This time-stepper
provides explicit symplectic partitioned Runge-Kutta methods up to order 10
for separable Hamiltonian systems.

Added support for relaxation Runge-Kutta methods in ERKStep and ARKStep in
ARKODE.

### New Features

Updated CVODE, CVODES and ARKODE default behavior when returning the solution when
the internal time has reached a user-specified stop time. Previously, the output
solution was interpolated to the value of `tstop`; the default is now to copy the
internal solution vector. Users who wish to revert to interpolation may call a new
routine `CVodeSetInterpolateStopTime`, `ARKStepSetInterpolateStopTime`,
`ERKStepSetInterpolateStopTime`, or `MRIStepSetInterpolateStopTime`.

Added the second order IMEX method from Giraldo, Kelly, and Constantinescu 2013
as the default second order IMEX method in ARKStep. The explicit table is given
by `ARKODE_ARK2_ERK_3_1_2` and the implicit table by `ARKODE_ARK2_DIRK_3_1_2`.

Updated the F2003 utility routines `SUNDIALSFileOpen` and `SUNDIALSFileClose`
to support user specification of `stdout` and `stderr` strings for the output
file names.

## Bug Fixes

A potential bug was fixed when using inequality constraint handling and
calling `ARKStepGetEstLocalErrors` or `ERKStepGetEstLocalErrors` after a failed
step in which an inequality constraint violation occurred. In this case, the
values returned by `ARKStepGetEstLocalErrors` or `ERKStepGetEstLocalErrors` may
have been invalid.

## Changes to SUNDIALS in release 6.5.1

### New Features

Added the functions `ARKStepClearStopTime`, `ERKStepClearStopTime`,
`MRIStepClearStopTime`, `CVodeClearStopTime`, and `IDAClearStopTime` to
disable a previously set stop time.

The default interpolant in ARKODE when using a first order method has been
updated to a linear interpolant to ensure values obtained by the integrator are
returned at the ends of the time interval. To restore the previous behavior of
using a constant interpolant call `ARKStepSetInterpolantDegree`,
`ERKStepSetInterpolantDegree`, or `MRIStepSetInterpolantDegree` and set the
interpolant degree to zero before evolving the problem.

### Bug Fixes

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the `icx` compiler.

## Changes to SUNDIALS in release 6.5.0

### New Features

A new capability to keep track of memory allocations made through the `SUNMemoryHelper`
classes has been added. Memory allocation stats can be accessed through the
`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the `SUNMemoryHelper` classes for more details.

Added the functions `ARKStepGetJac`, `ARKStepGetJacTime`,
`ARKStepGetJacNumSteps`, `MRIStepGetJac`, `MRIStepGetJacTime`,
`MRIStepGetJacNumSteps`, `CVodeGetJac`, `CVodeGetJacTime`,
`CVodeGetJacNumSteps`, `IDAGetJac`, `IDAGetJacCj`, `IDAGetJacTime`,
`IDAGetJacNumSteps`, `KINGetJac`, `KINGetJacNumIters` to assist in
debugging simulations utilizing a matrix-based linear solver.

Added support for CUDA 12.

Added support for the SYCL backend with RAJA 2022.x.y.

### Bug Fixes

Fixed an underflow bug during root finding in ARKODE, CVODE, CVODES, IDA and
IDAS. This fixes [GitHub Issue #57](https://github.com/LLNL/sundials/issues/57>).

Fixed an issue with finding oneMKL when using the `icpx` compiler with the
`-fsycl` flag as the C++ compiler instead of `dpcpp`.

Fixed the shape of the arrays returned by `FN_VGetArrayPointer` functions as well
as the `FSUNDenseMatrix_Data`, `FSUNBandMatrix_Data`, `FSUNSparseMatrix_Data`,
`FSUNSparseMatrix_IndexValues`, and `FSUNSparseMatrix_IndexPointers` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

Fixed an implicit conversion error in the Butcher table for ESDIRK5(4)7L[2]SA2.

## Changes to SUNDIALS in release 6.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial `N_Vector`.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

## Changes to SUNDIALS in release 6.4.0

### New Requirements

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, Ginkgo, and Kokkos.

### Major Features

Added support for the [Ginkgo](https://ginkgo-project.github.io/) linear algebra
library. This support includes new `SUNMatrix` and `SUNLinearSolver`
implementations, see the `SUNMATRIX_GINKGO` and `SUNLINEARSOLVER_GINKGO`
sections in the documentation for more information.

Added new `NVector`, dense `SUNMatrix`, and dense `SUNLinearSolver`
implementations utilizing [Kokkos Ecosystem](https://kokkos.org/) for
performance portability, see the `NVECTOR_KOKKOS`, `SUNMATRIX_KOKKOSDENSE` and
`SUNLINEARSOLVER_KOKKOSDENSE` sections in the documentation for more
information.

### New Features

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added the functions `ARKStepSetTableName`, `ERKStepSetTableName`,
`MRIStepCoupling_LoadTableByName`, `ARKodeButcherTable_LoadDIRKByName`, and
`ARKodeButcherTable_LoadERKByName` to load a table from a string.

### Bug Fixes

Fixed a bug in the CUDA and HIP vectors where `N_VMaxNorm` would return the
minimum positive floating-point value for the zero vector.

Fixed memory leaks/out of bounds memory accesses in the ARKODE MRIStep module
that could occur when attaching a coupling table after reinitialization with a
different number of stages than originally selected.

Fixed a memory leak in CVODE and CVODES where the projection memory would not be
deallocated when calling `CVodeFree`.

## Changes to SUNDIALS in release 6.3.0

### New Features

Added `GetUserData` functions in each package to retrieve the user data pointer
provided to `SetUserData` functions. See `ARKStepGetUserData`,
`ERKStepGetUserData`, `MRIStepGetUserData`, `CVodeGetUserData`,
`IDAGetUserData`, or `KINGetUserData` for more information.

Added a variety of embedded DIRK methods from [Kennedy & Carpenter,
NASA TM-2016-219173, 2016] and [Kennedy & Carpenter, Appl. Numer. Math., 146, 2019] to
ARKODE.

Updated `MRIStepReset` to call the corresponding `MRIStepInnerResetFn` with the same
(*tR*,*yR*) arguments for the `MRIStepInnerStepper` object that is used to evolve the
MRI "fast" time scale subproblems.

Added a new [example](examples/cvode/serial/cvRocket_dns.c) which
demonstrates using CVODE with a discontinuous right-hand-side function
and rootfinding.

### Bug Fixes

Fixed a bug in `ERKStepReset`, `ERKStepReInit`, `ARKStepReset`, `ARKStepReInit`,
`MRIStepReset`, and `MRIStepReInit` where a previously-set value of *tstop* (from
a call to `ERKStepSetStopTime`, `ARKStepSetStopTime`, or `MRIStepSetStopTime`,
respectively) would not be cleared.

Fixed the unituitive behavior of the `USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option `SUNDIALS_MATH_LIBRARY`.

Changed `SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'. This
fixes [GitHub Issue #177](https://github.com/LLNL/sundials/issues/177).

## Changes to SUNDIALS in release 6.2.0

### Major Features

Added the `SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Added support to CVODES for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user defined
projection function. This implementation is accompanied by additions to the
CVODES user documentation and examples.

### New Features

Added the function `SUNProfiler_Reset` to reset the region timings and counters
to zero.

Added the following functions to output all of the integrator, nonlinear solver,
linear solver, and other statistics in one call:

* `ARKStepPrintAllStats`
* `ERKStepPrintAllStats`
* `MRIStepPrintAllStats`
* `CVodePrintAllStats`
* `IDAPrintAllStats`
* `KINPrintAllStats`

The file `scripts/sundials_csv.py` contains functions for parsing the
comma-separated value (CSV) output files when using the CSV output format.

Added functions to CVODE, CVODES, IDA, and IDAS to change the default step size
adaptivity parameters. For more information see the documentation for:

* `CVodeSetEtaFixedStepBounds`
* `CVodeSetEtaMaxFirstStep`
* `CVodeSetEtaMaxEarlyStep`
* `CVodeSetNumStepsEtaMaxEarlyStep`
* `CVodeSetEtaMax`
* `CVodeSetEtaMin`
* `CVodeSetEtaMinErrFail`
* `CVodeSetEtaMaxErrFail`
* `CVodeSetNumFailsEtaMaxErrFail`
* `CVodeSetEtaConvFail`
* `IDASetEtaFixedStepBounds`
* `IDAsetEtaMax`
* `IDASetEtaMin`
* `IDASetEtaLow`
* `IDASetEtaMinErrFail`
* `IDASetEtaConvFail`

Added the functions `ARKStepSetDeduceImplicitRhs` and
`MRIStepSetDeduceImplicitRhs` to optionally remove an evaluation of the implicit
right-hand side function after nonlinear solves. See the mathematical
considerations section of the user guide for information on using this
optimization.

Added the function `MRIStepSetOrder` to select the default MRI method of a given
order.

Added the functions `CVodeSetDeltaGammaMaxLSetup` and
`CVodeSetDeltaGammaMaxBadJac` in CVODE and CVODES to adjust the `gamma` change
thresholds to require a linear solver setup or Jacobian/precondition update,
respectively.

Added the function `IDASetDetlaCjLSetup` in IDA and IDAS to adjust the parameter
that determines when a change in `c_j` requires calling the linear solver setup
function.

Added the function `IDASetMinStep` to set a minimum step size.

### Bug Fixes

Fixed the `SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

The behavior of `N_VSetKernelExecPolicy_Sycl` has been updated to be consistent
with the CUDA and HIP vectors. The input execution policies are now cloned and
may be freed after calling `N_VSetKernelExecPolicy_Sycl`. Additionally, `NULL`
inputs are now allowed and, if provided, will reset the vector execution
policies to the defaults.

A memory leak in the SYCL vector was fixed where the execution policies were not
freed when the vector was destroyed.

The include guard in `nvector_mpimanyvector.h` has been corrected to enable
using both the ManyVector and MPIManyVector vector implementations in the same
simulation.

A bug was fixed in the ARKODE, CVODE(S), and IDA(S) functions to retrieve the
number of nonlinear solver failures. The failure count returned was the number
of failed *steps* due to a nonlinear solver failure i.e., if a nonlinear solve
failed with a stale Jacobian or preconditioner but succeeded after updating the
Jacobian or preconditioner, the initial failure was not included in the
nonlinear solver failure count. The following functions have been updated to
return the total number of nonlinear solver failures:

* `ARKStepGetNumNonlinSolvConvFails`
* `ARKStepGetNonlinSolvStats`
* `MRIStepGetNumNonlinSolvConvFails`
* `MRIStepGetNonlinSolvStats`
* `CVodeGetNumNonlinSolvConvFails`
* `CVodeGetNonlinSolvStats`
* `CVodeGetSensNumNonlinSolvConvFails`
* `CVodeGetSensNonlinSolvStats`
* `CVodeGetStgrSensNumNonlinSolvConvFails`
* `CVodeGetStgrSensNonlinSolvStats`
* `IDAGetNumNonlinSolvConvFails`
* `IDAGetNonlinSolvStats`
* `IDAGetSensNumNonlinSolvConvFails`
* `IDAGetSensNonlinSolvStats`

As a result of this change users may see an increase in the number of failures
reported from the above functions. The following functions have been added to
retrieve the number of failed steps due to a nonlinear solver failure i.e., the
counts previously returned by the above functions:

* `ARKStepGetNumStepSolveFails`
* `MRIStepGetNumStepSolveFails`
* `CVodeGetNumStepSolveFails`
* `CVodeGetNumStepSensSolveFails`
* `CVodeGetNumStepStgrSensSolveFails`
* `IDAGetNumStepSolveFails`
* `IDAGetNumStepSensSolveFails`

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

### Deprecation Notice

Deprecated the following functions, it is recommended to use the `SUNLogger` API
instead.

* `ARKStepSetDiagnostics`
* `ERKStepSetDiagnostics`
* `MRIStepSetDiagnostics`
* `KINSetInfoFile`
* `SUNNonlinSolSetPrintLevel_Newton`
* `SUNNonlinSolSetInfoFile_Newton`
* `SUNNonlinSolSetPrintLevel_FixedPoint`
* `SUNNonlinSolSetInfoFile_FixedPoint`
* `SUNLinSolSetInfoFile_PCG`
* `SUNLinSolSetPrintLevel_PCG`
* `SUNLinSolSetInfoFile_SPGMR`
* `SUNLinSolSetPrintLevel_SPGMR`
* `SUNLinSolSetInfoFile_SPFGMR`
* `SUNLinSolSetPrintLevel_SPFGMR`
* `SUNLinSolSetInfoFile_SPTFQM`
* `SUNLinSolSetPrintLevel_SPTFQMR`
* `SUNLinSolSetInfoFile_SPBCGS`
* `SUNLinSolSetPrintLevel_SPBCGS`

The `SUNLinSolSetInfoFile_*` and `SUNNonlinSolSetInfoFile_*` family of
functions are now enabled by setting the CMake option `SUNDIALS_LOGGING_LEVEL`
to a value `>= 3`.

## Changes to SUNDIALS in release 6.1.1

### New Features

Added new Fortran example program,
`examples/arkode/F2003_serial/ark_kpr_mri_f2003.f90` demonstrating MRI
capabilities.

### Bug Fixes

Fixed exported `SUNDIALSConfig.cmake`.

Fixed Fortran interface to `MRIStepInnerStepper` and `MRIStepCoupling`
structures and functions.

## Changes to SUNDIALS in release 6.1.0

### New Features

Added new reduction implementations for the CUDA and HIP vectors that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the `N_VSetKernelExecPolicy_Cuda` and
`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

`SUNDIALS::<lib>` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

`CMAKE_C_STANDARD` is now set to 99 by default.

### Bug Fixes

Fixed exported `SUNDIALSConfig.cmake` when profiling is enabled without Caliper.

Fixed `sundials_export.h` include in `sundials_config.h`.

Fixed memory leaks in the SuperLU_MT linear solver interface.

## Changes to SUNDIALS in release 6.0.0

### Breaking Changes

#### SUNContext Object Added

SUNDIALS v6.0.0 introduces a new `SUNContext` object on which all other SUNDIALS
objects depend. As such, the constructors for all SUNDIALS packages, vectors,
matrices, linear solvers, nonlinear solvers, and memory helpers have been
updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call `SUNContext_Create` to create a context object with
before calling any other SUNDIALS library function, and then provide this object
to other SUNDIALS constructors. The context object has been introduced to allow
SUNDIALS to provide new features, such as the profiling/instrumentation also
introduced in this release, while maintaining thread-safety. See the
documentation section on the `SUNContext` for more details.

The script `upgrade-to-sundials-6-from-5.sh` has been provided with this release
(and obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a `SUNCTX_PLACEHOLDER` argument to all of
the calls to SUNDIALS constructors that now require a `SUNContext` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

```
> ./upgrade-to-sundials-6-from-5.sh <files to update>
```

#### Updated SUNMemoryHelper Function Signatures

The `SUNMemoryHelper` functions `Alloc`, `Dealloc`, and `Copy` have been updated
to accept an opaque handle as the last input. At a minimum, existing
`SUNMemoryHelper` implementations will need to update these functions to accept
the additional argument. Typically, this handle is the execution stream (e.g., a
CUDA/HIP stream or SYCL queue) for the operation. The CUDA, HIP, and SYCL
`SUNMemoryHelper` implementations have been updated accordingly. Additionally,
the constructor for the SYCL implementation has been updated to remove the SYCL
queue as an input.

#### Deprecated Functions Removed

The previously deprecated constructor `N_VMakeWithManagedAllocator_Cuda` and
the function `N_VSetCudaStream_Cuda` have been removed and replaced with
`N_VNewWithMemHelp_Cuda` and `N_VSetKernelExecPolicy_Cuda` respectively.

The previously deprecated macros `PVEC_REAL_MPI_TYPE` and
`PVEC_INTEGER_MPI_TYPE` have been removed and replaced with
`MPI_SUNREALTYPE` and `MPI_SUNINDEXTYPE` respectively.

The following previously deprecated functions have been removed

| Removed                   | Replaced with                    |
|:--------------------------|:---------------------------------|
| `SUNBandLinearSolver`     | `SUNLinSol_Band`                 |
| `SUNDenseLinearSolver`    | `SUNLinSol_Dense`                |
| `SUNKLU`                  | `SUNLinSol_KLU`                  |
| `SUNKLUReInit`            | `SUNLinSol_KLUReInit`            |
| `SUNKLUSetOrdering`       | `SUNLinSol_KLUSetOrdering`       |
| `SUNLapackBand`           | `SUNLinSol_LapackBand`           |
| `SUNLapackDense`          | `SUNLinSol_LapackDense`          |
| `SUNPCG`                  | `SUNLinSol_PCG`                  |
| `SUNPCGSetPrecType`       | `SUNLinSol_PCGSetPrecType`       |
| `SUNPCGSetMaxl`           | `SUNLinSol_PCGSetMaxl`           |
| `SUNSPBCGS`               | `SUNLinSol_SPBCGS`               |
| `SUNSPBCGSSetPrecType`    | `SUNLinSol_SPBCGSSetPrecType`    |
| `SUNSPBCGSSetMaxl`        | `SUNLinSol_SPBCGSSetMaxl`        |
| `SUNSPFGMR`               | `SUNLinSol_SPFGMR`               |
| `SUNSPFGMRSetPrecType`    | `SUNLinSol_SPFGMRSetPrecType`    |
| `SUNSPFGMRSetGSType`      | `SUNLinSol_SPFGMRSetGSType`      |
| `SUNSPFGMRSetMaxRestarts` | `SUNLinSol_SPFGMRSetMaxRestarts` |
| `SUNSPGMR`                | `SUNLinSol_SPGMR`                |
| `SUNSPGMRSetPrecType`     | `SUNLinSol_SPGMRSetPrecType`     |
| `SUNSPGMRSetGSType`       | `SUNLinSol_SPGMRSetGSType`       |
| `SUNSPGMRSetMaxRestarts`  | `SUNLinSol_SPGMRSetMaxRestarts`  |
| `SUNSPTFQMR`              | `SUNLinSol_SPTFQMR`              |
| `SUNSPTFQMRSetPrecType`   | `SUNLinSol_SPTFQMRSetPrecType`   |
| `SUNSPTFQMRSetMaxl`       | `SUNLinSol_SPTFQMRSetMaxl`       |
| `SUNSuperLUMT`            | `SUNLinSol_SuperLUMT`            |
| `SUNSuperLUMTSetOrdering` | `SUNLinSol_SuperLUMTSetOrdering` |

The deprecated functions `MRIStepGetCurrentButcherTables` and
`MRIStepWriteButcher` and the utility functions `MRIStepSetTable` and
`MRIStepSetTableNum` have been removed. Users wishing to create an MRI-GARK
method from a Butcher table should use `MRIStepCoupling_MIStoMRI` to create
the corresponding MRI coupling table and attach it with `MRIStepSetCoupling`.

The previously deprecated functions `ARKStepSetMaxStepsBetweenLSet` and
`ARKStepSetMaxStepsBetweenJac` have been removed and replaced with
`ARKStepSetLSetupFrequency` and `ARKStepSetMaxStepsBetweenJac` respectively.

The previously deprecated function `CVodeSetMaxStepsBetweenJac` has been removed
and replaced with `CVodeSetJacEvalFrequency`.

The ARKODE, CVODE, IDA, and KINSOL Fortran 77 interfaces have been removed. See
the "SUNDIALS Fortran Interface" section in the user guides and the F2003
example programs for more details using the SUNDIALS Fortran 2003 module
interfaces.

#### Namespace Changes

The CUDA, HIP, and SYCL execution policies have been moved from the `sundials`
namespace to the `sundials::cuda`, `sundials::hip`, and `sundials::sycl`
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The `Sundials` namespace used by the Trilinos Tpetra NVector has been replaced
with the `sundials::trilinos::nvector_tpetra` namespace.

### Major Features

#### SUNProfiler

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option `SUNDIALS_BUILD_WITH_PROFILING`. A built-in
profiler will be used by default, but the
[Caliper](https://github.com/LLNL/Caliper) library can also be used instead with
the CMake option `ENABLE_CALIPER`. See the documentation section on profiling
for more details. **WARNING**: Profiling will impact performance, and should be
enabled judiciously.

#### IMEX MRI Methods and MRIStepInnerStepper Object

The ARKODE MRIStep module has been extended to support implicit-explicit (IMEX)
multirate infinitesimal generalized additive Runge-Kutta (MRI-GARK) methods. As
such, `MRIStepCreate` has been updated to include arguments for the slow
explicit and slow implicit ODE right-hand side functions. `MRIStepCreate` has
also been updated to require attaching an `MRIStepInnerStepper` for evolving the
fast time scale. `MRIStepReInit` has been similarly updated to take explicit
and implicit right-hand side functions as input. Codes using explicit or
implicit MRI methods will need to update `MRIStepCreate` and `MRIStepReInit`
calls to pass `NULL` for either the explicit or implicit right-hand side
function as appropriate. If ARKStep is used as the fast time scale integrator,
codes will need to call `ARKStepCreateMRIStepInnerStepper` to wrap the ARKStep
memory as an `MRIStepInnerStepper` object. Additionally, `MRIStepGetNumRhsEvals`
has been updated to return the number of slow implicit and explicit function
evaluations. The coupling table structure `MRIStepCouplingMem` and the
functions `MRIStepCoupling_Alloc` and `MRIStepCoupling_Create` have also
been updated to support IMEX-MRI-GARK methods.

### New Features

Two new optional vector operations, `N_VDotProdMultiLocal` and
`N_VDotProdMultiAllReduce`, have been added to support low-synchronization
methods for Anderson acceleration.

The implementation of solve-decoupled implicit MRI-GARK methods has been updated
to remove extraneous slow implicit function calls and reduce the memory
requirements.

Added a new function `CVodeGetLinSolveStats` to get the CVODES linear solver
statistics as a group.

Added a new function, `CVodeSetMonitorFn`, that takes a user-function
to be called by CVODES after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODES statistics
throughout the simulation.

New orthogonalization methods were added for use within Anderson acceleration
in KINSOL. See the "Anderson Acceleration QR Factorization" subsection within
the mathematical considerations chapter of the user guide and the
`KINSetOrthAA` function documentation for more details.

### Deprecation Notice

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions `N_VCloneVectorArray_*` and `N_VDestroyVectorArray_*` have been
deprecated. The generic `N_VCloneVectorArray` and `N_VDestroyVectorArray`
functions should be used instead.

Many constants, types, and functions have been renamed so that they are properly
namespaced. The old names have been deprecated and will be removed in SUNDIALS
v7.0.0.

The following constants, macros, and typedefs are now deprecated:

| Deprecated Name            | New Name                          |
|:---------------------------|:----------------------------------|
| `realtype`                 | `sunrealtype`                     |
| `booleantype`              | `sunbooleantype`                  |
| `RCONST`                   | `SUN_RCONST`                      |
| `BIG_REAL`                 | `SUN_BIG_REAL`                    |
| `SMALL_REAL`               | `SUN_SMALL_REAL`                  |
| `UNIT_ROUNDOFF`            | `SUN_UNIT_ROUNDOFF`               |
| `PREC_NONE`                | `SUN_PREC_NONE`                   |
| `PREC_LEFT`                | `SUN_PREC_LEFT`                   |
| `PREC_RIGHT`               | `SUN_PREC_RIGHT`                  |
| `PREC_BOTH`                | `SUN_PREC_BOTH`                   |
| `MODIFIED_GS`              | `SUN_MODIFIED_GS`                 |
| `CLASSICAL_GS`             | `SUN_CLASSICAL_GS`                |
| `ATimesFn`                 | `SUNATimesFn`                     |
| `PSetupFn`                 | `SUNPSetupFn`                     |
| `PSolveFn`                 | `SUNPSolveFn`                     |
| `DlsMat`                   | `SUNDlsMat`                       |
| `DENSE_COL`                | `SUNDLS_DENSE_COL`                |
| `DENSE_ELEM`               | `SUNDLS_DENSE_ELEM`               |
| `BAND_COL`                 | `SUNDLS_BAND_COL`                 |
| `BAND_COL_ELEM`            | `SUNDLS_BAND_COL_ELEM`            |
| `BAND_ELEM`                | `SUNDLS_BAND_ELEM`                |
| `SDIRK_2_1_2`              | `ARKODE_SDIRK_2_1_2`              |
| `BILLINGTON_3_3_2`         | `ARKODE_BILLINGTON_3_3_2`         |
| `TRBDF2_3_3_2`             | `ARKODE_TRBDF2_3_3_2`             |
| `KVAERNO_4_2_3`            | `ARKODE_KVAERNO_4_2_3`            |
| `ARK324L2SA_DIRK_4_2_3`    | `ARKODE_ARK324L2SA_DIRK_4_2_3`    |
| `CASH_5_2_4`               | `ARKODE_CASH_5_2_4`               |
| `CASH_5_3_4`               | `ARKODE_CASH_5_3_4`               |
| `SDIRK_5_3_4`              | `ARKODE_SDIRK_5_3_4`              |
| `KVAERNO_5_3_4`            | `ARKODE_KVAERNO_5_3_4`            |
| `ARK436L2SA_DIRK_6_3_4`    | `ARKODE_ARK436L2SA_DIRK_6_3_4`    |
| `KVAERNO_7_4_5`            | `ARKODE_KVAERNO_7_4_5`            |
| `ARK548L2SA_DIRK_8_4_5`    | `ARKODE_ARK548L2SA_DIRK_8_4_5`    |
| `ARK437L2SA_DIRK_7_3_4`    | `ARKODE_ARK437L2SA_DIRK_7_3_4`    |
| `ARK548L2SAb_DIRK_8_4_5`   | `ARKODE_ARK548L2SAb_DIRK_8_4_5`   |
| `MIN_DIRK_NUM`             | `ARKODE_MIN_DIRK_NUM`             |
| `MAX_DIRK_NUM`             | `ARKODE_MAX_DIRK_NUM`             |
| `MIS_KW3`                  | `ARKODE_MIS_KW3`                  |
| `MRI_GARK_ERK33a`          | `ARKODE_MRI_GARK_ERK33a`          |
| `MRI_GARK_ERK45a`          | `ARKODE_MRI_GARK_ERK45a`          |
| `MRI_GARK_IRK21a`          | `ARKODE_MRI_GARK_IRK21a`          |
| `MRI_GARK_ESDIRK34a`       | `ARKODE_MRI_GARK_ESDIRK34a`       |
| `MRI_GARK_ESDIRK46a`       | `ARKODE_MRI_GARK_ESDIRK46a`       |
| `IMEX_MRI_GARK3a`          | `ARKODE_IMEX_MRI_GARK3a`          |
| `IMEX_MRI_GARK3b`          | `ARKODE_IMEX_MRI_GARK3b`          |
| `IMEX_MRI_GARK4`           | `ARKODE_IMEX_MRI_GARK4`           |
| `MIN_MRI_NUM`              | `ARKODE_MIN_MRI_NUM`              |
| `MAX_MRI_NUM`              | `ARKODE_MAX_MRI_NUM`              |
| `DEFAULT_MRI_TABLE_3`      | `MRISTEP_DEFAULT_TABLE_3`         |
| `DEFAULT_EXPL_MRI_TABLE_3` | `MRISTEP_DEFAULT_EXPL_TABLE_3`    |
| `DEFAULT_EXPL_MRI_TABLE_4` | `MRISTEP_DEFAULT_EXPL_TABLE_4`    |
| `DEFAULT_IMPL_SD_TABLE_2`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_2` |
| `DEFAULT_IMPL_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_3` |
| `DEFAULT_IMPL_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_4` |
| `DEFAULT_IMEX_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_3` |
| `DEFAULT_IMEX_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_4` |
| `HEUN_EULER_2_1_2`         | `ARKODE_HEUN_EULER_2_1_2`         |
| `BOGACKI_SHAMPINE_4_2_3`   | `ARKODE_BOGACKI_SHAMPINE_4_2_3`   |
| `ARK324L2SA_ERK_4_2_3`     | `ARKODE_ARK324L2SA_ERK_4_2_3`     |
| `ZONNEVELD_5_3_4`          | `ARKODE_ZONNEVELD_5_3_4`          |
| `ARK436L2SA_ERK_6_3_4`     | `ARKODE_ARK436L2SA_ERK_6_3_4`     |
| `SAYFY_ABURUB_6_3_4`       | `ARKODE_SAYFY_ABURUB_6_3_4`       |
| `CASH_KARP_6_4_5`          | `ARKODE_CASH_KARP_6_4_5`          |
| `FEHLBERG_6_4_5`           | `ARKODE_FEHLBERG_6_4_5`           |
| `DORMAND_PRINCE_7_4_5`     | `ARKODE_DORMAND_PRINCE_7_4_5`     |
| `ARK548L2SA_ERK_8_4_5`     | `ARKODE_ARK548L2SA_ERK_8_4_5`     |
| `VERNER_8_5_6`             | `ARKODE_VERNER_8_5_6`             |
| `FEHLBERG_13_7_8`          | `ARKODE_FEHLBERG_13_7_8`          |
| `KNOTH_WOLKE_3_3`          | `ARKODE_KNOTH_WOLKE_3_3`          |
| `ARK437L2SA_ERK_7_3_4`     | `ARKODE_ARK437L2SA_ERK_7_3_4`     |
| `ARK548L2SAb_ERK_8_4_5`    | `ARKODE_ARK548L2SAb_ERK_8_4_5`    |
| `MIN_ERK_NUM`              | `ARKODE_MIN_ERK_NUM`              |
| `MAX_ERK_NUM`              | `ARKODE_MAX_ERK_NUM`              |
| `DEFAULT_ERK_2`            | `ARKSTEP_DEFAULT_ERK_2`           |
| `DEFAULT_ERK_3`            | `ARKSTEP_DEFAULT_ERK_3`           |
| `DEFAULT_ERK_4`            | `ARKSTEP_DEFAULT_ERK_4`           |
| `DEFAULT_ERK_5`            | `ARKSTEP_DEFAULT_ERK_5`           |
| `DEFAULT_ERK_6`            | `ARKSTEP_DEFAULT_ERK_6`           |
| `DEFAULT_ERK_8`            | `ARKSTEP_DEFAULT_ERK_8`           |
| `DEFAULT_DIRK_2`           | `ARKSTEP_DEFAULT_DIRK_2`          |
| `DEFAULT_DIRK_3`           | `ARKSTEP_DEFAULT_DIRK_3`          |
| `DEFAULT_DIRK_4`           | `ARKSTEP_DEFAULT_DIRK_4`          |
| `DEFAULT_DIRK_5`           | `ARKSTEP_DEFAULT_DIRK_5`          |
| `DEFAULT_ARK_ETABLE_3`     | `ARKSTEP_DEFAULT_ARK_ETABLE_3`    |
| `DEFAULT_ARK_ETABLE_4`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ETABLE_5`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ITABLE_3`     | `ARKSTEP_DEFAULT_ARK_ITABLE_3`    |
| `DEFAULT_ARK_ITABLE_4`     | `ARKSTEP_DEFAULT_ARK_ITABLE_4`    |
| `DEFAULT_ARK_ITABLE_5`     | `ARKSTEP_DEFAULT_ARK_ITABLE_5`    |
| `DEFAULT_ERK_2`            | `ERKSTEP_DEFAULT_2`               |
| `DEFAULT_ERK_3`            | `ERKSTEP_DEFAULT_3`               |
| `DEFAULT_ERK_4`            | `ERKSTEP_DEFAULT_4`               |
| `DEFAULT_ERK_5`            | `ERKSTEP_DEFAULT_5`               |
| `DEFAULT_ERK_6`            | `ERKSTEP_DEFAULT_6`               |
| `DEFAULT_ERK_8`            | `ERKSTEP_DEFAULT_8`               |

In addition, the following functions are now deprecated (compile-time warnings
will be printed if supported by the compiler):

| Deprecated Name               | New Name                     |
|:------------------------------|:-----------------------------|
| `CVSpilsSetLinearSolver`      | `CVodeSetLinearSolver`       |
| `CVSpilsSetEpsLin`            | `CVodeSetEpsLin`             |
| `CVSpilsSetPreconditioner`    | `CVodeSetPreconditioner`     |
| `CVSpilsSetJacTimes`          | `CVodeSetJacTimes`           |
| `CVSpilsGetWorkSpace`         | `CVodeGetLinWorkSpace`       |
| `CVSpilsGetNumPrecEvals`      | `CVodeGetNumPrecEvals`       |
| `CVSpilsGetNumPrecSolves`     | `CVodeGetNumPrecSolves`      |
| `CVSpilsGetNumLinIters`       | `CVodeGetNumLinIters`        |
| `CVSpilsGetNumConvFails`      | `CVodeGetNumConvFails`       |
| `CVSpilsGetNumJTSetupEvals`   | `CVodeGetNumJTSetupEvals`    |
| `CVSpilsGetNumJtimesEvals`    | `CVodeGetNumJtimesEvals`     |
| `CVSpilsGetNumRhsEvals`       | `CVodeGetNumLinRhsEvals`     |
| `CVSpilsGetLastFlag`          | `CVodeGetLastLinFlag`        |
| `CVSpilsGetReturnFlagName`    | `CVodeGetLinReturnFlagName`  |
| `CVSpilsSetLinearSolverB`     | `CVodeSetLinearSolverB`      |
| `CVSpilsSetEpsLinB`           | `CVodeSetEpsLinB`            |
| `CVSpilsSetPreconditionerB`   | `CVodeSetPreconditionerB`    |
| `CVSpilsSetPreconditionerBS`  | `CVodeSetPreconditionerBS`   |
| `CVSpilsSetJacTimesB`         | `CVodeSetJacTimesB`          |
| `CVSpilsSetJacTimesBS`        | `CVodeSetJacTimesBS`         |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `CVDlsSetLinearSolverB`       | `CVodeSetLinearSolverB`      |
| `CVDlsSetJacFnB`              | `CVodeSetJacFnB`             |
| `CVDlsSetJacFnBS`             | `CVodeSetJacFnBS`            |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `KINDlsSetLinearSolver`       | `KINSetLinearSolver`         |
| `KINDlsSetJacFn`              | `KINSetJacFn`                |
| `KINDlsGetWorkSpace`          | `KINGetLinWorkSpace`         |
| `KINDlsGetNumJacEvals`        | `KINGetNumJacEvals`          |
| `KINDlsGetNumFuncEvals`       | `KINGetNumLinFuncEvals`      |
| `KINDlsGetLastFlag`           | `KINGetLastLinFlag`          |
| `KINDlsGetReturnFlagName`     | `KINGetLinReturnFlagName`    |
| `KINSpilsSetLinearSolver`     | `KINSetLinearSolver`         |
| `KINSpilsSetPreconditioner`   | `KINSetPreconditioner`       |
| `KINSpilsSetJacTimesVecFn`    | `KINSetJacTimesVecFn`        |
| `KINSpilsGetWorkSpace`        | `KINGetLinWorkSpace`         |
| `KINSpilsGetNumPrecEvals`     | `KINGetNumPrecEvals`         |
| `KINSpilsGetNumPrecSolves`    | `KINGetNumPrecSolves`        |
| `KINSpilsGetNumLinIters`      | `KINGetNumLinIters`          |
| `KINSpilsGetNumConvFails`     | `KINGetNumLinConvFails`      |
| `KINSpilsGetNumJtimesEvals`   | `KINGetNumJtimesEvals`       |
| `KINSpilsGetNumFuncEvals`     | `KINGetNumLinFuncEvals`      |
| `KINSpilsGetLastFlag`         | `KINGetLastLinFlag`          |
| `KINSpilsGetReturnFlagName`   | `KINGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolver`     | `IDASetLinearSolver`         |
| `IDASpilsSetPreconditioner`   | `IDASetPreconditioner`       |
| `IDASpilsSetJacTimes`         | `IDASetJacTimes`             |
| `IDASpilsSetEpsLin`           | `IDASetEpsLin`               |
| `IDASpilsSetIncrementFactor`  | `IDASetIncrementFactor`      |
| `IDASpilsGetWorkSpace`        | `IDAGetLinWorkSpace`         |
| `IDASpilsGetNumPrecEvals`     | `IDAGetNumPrecEvals`         |
| `IDASpilsGetNumPrecSolves`    | `IDAGetNumPrecSolves`        |
| `IDASpilsGetNumLinIters`      | `IDAGetNumLinIters`          |
| `IDASpilsGetNumConvFails`     | `IDAGetNumLinConvFails`      |
| `IDASpilsGetNumJTSetupEvals`  | `IDAGetNumJTSetupEvals`      |
| `IDASpilsGetNumJtimesEvals`   | `IDAGetNumJtimesEvals`       |
| `IDASpilsGetNumResEvals`      | `IDAGetNumLinResEvals`       |
| `IDASpilsGetLastFlag`         | `IDAGetLastLinFlag`          |
| `IDASpilsGetReturnFlagName`   | `IDAGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolverB`    | `IDASetLinearSolverB`        |
| `IDASpilsSetEpsLinB`          | `IDASetEpsLinB`              |
| `IDASpilsSetIncrementFactorB` | `IDASetIncrementFactorB`     |
| `IDASpilsSetPreconditionerB`  | `IDASetPreconditionerB`      |
| `IDASpilsSetPreconditionerBS` | `IDASetPreconditionerBS`     |
| `IDASpilsSetJacTimesB`        | `IDASetJacTimesB`            |
| `IDASpilsSetJacTimesBS`       | `IDASetJacTimesBS`           |
| `IDADlsSetLinearSolver`       | `IDASetLinearSolver`         |
| `IDADlsSetJacFn`              | `IDASetJacFn`                |
| `IDADlsGetWorkSpace`          | `IDAGetLinWorkSpace`         |
| `IDADlsGetNumJacEvals`        | `IDAGetNumJacEvals`          |
| `IDADlsGetNumResEvals`        | `IDAGetNumLinResEvals`       |
| `IDADlsGetLastFlag`           | `IDAGetLastLinFlag`          |
| `IDADlsGetReturnFlagName`     | `IDAGetLinReturnFlagName`    |
| `IDADlsSetLinearSolverB`      | `IDASetLinearSolverB`        |
| `IDADlsSetJacFnB`             | `IDASetJacFnB`               |
| `IDADlsSetJacFnBS`            | `IDASetJacFnBS`              |
| `DenseGETRF`                  | `SUNDlsMat_DenseGETRF`       |
| `DenseGETRS`                  | `SUNDlsMat_DenseGETRS`       |
| `denseGETRF`                  | `SUNDlsMat_denseGETRF`       |
| `denseGETRS`                  | `SUNDlsMat_denseGETRS`       |
| `DensePOTRF`                  | `SUNDlsMat_DensePOTRF`       |
| `DensePOTRS`                  | `SUNDlsMat_DensePOTRS`       |
| `densePOTRF`                  | `SUNDlsMat_densePOTRF`       |
| `densePOTRS`                  | `SUNDlsMat_densePOTRS`       |
| `DenseGEQRF`                  | `SUNDlsMat_DenseGEQRF`       |
| `DenseORMQR`                  | `SUNDlsMat_DenseORMQR`       |
| `denseGEQRF`                  | `SUNDlsMat_denseGEQRF`       |
| `denseORMQR`                  | `SUNDlsMat_denseORMQR`       |
| `DenseCopy`                   | `SUNDlsMat_DenseCopy`        |
| `denseCopy`                   | `SUNDlsMat_denseCopy`        |
| `DenseScale`                  | `SUNDlsMat_DenseScale`       |
| `denseScale`                  | `SUNDlsMat_denseScale`       |
| `denseAddIdentity`            | `SUNDlsMat_denseAddIdentity` |
| `DenseMatvec`                 | `SUNDlsMat_DenseMatvec`      |
| `denseMatvec`                 | `SUNDlsMat_denseMatvec`      |
| `BandGBTRF`                   | `SUNDlsMat_BandGBTRF`        |
| `bandGBTRF`                   | `SUNDlsMat_bandGBTRF`        |
| `BandGBTRS`                   | `SUNDlsMat_BandGBTRS`        |
| `bandGBTRS`                   | `SUNDlsMat_bandGBTRS`        |
| `BandCopy`                    | `SUNDlsMat_BandCopy`         |
| `bandCopy`                    | `SUNDlsMat_bandCopy`         |
| `BandScale`                   | `SUNDlsMat_BandScale`        |
| `bandScale`                   | `SUNDlsMat_bandScale`        |
| `bandAddIdentity`             | `SUNDlsMat_bandAddIdentity`  |
| `BandMatvec`                  | `SUNDlsMat_BandMatvec`       |
| `bandMatvec`                  | `SUNDlsMat_bandMatvec`       |
| `ModifiedGS`                  | `SUNModifiedGS`              |
| `ClassicalGS`                 | `SUNClassicalGS`             |
| `QRfact`                      | `SUNQRFact`                  |
| `QRsol`                       | `SUNQRsol`                   |
| `DlsMat_NewDenseMat`          | `SUNDlsMat_NewDenseMat`      |
| `DlsMat_NewBandMat`           | `SUNDlsMat_NewBandMat`       |
| `DestroyMat`                  | `SUNDlsMat_DestroyMat`       |
| `NewIntArray`                 | `SUNDlsMat_NewIntArray`      |
| `NewIndexArray`               | `SUNDlsMat_NewIndexArray`    |
| `NewRealArray`                | `SUNDlsMat_NewRealArray`     |
| `DestroyArray`                | `SUNDlsMat_DestroyArray`     |
| `AddIdentity`                 | `SUNDlsMat_AddIdentity`      |
| `SetToZero`                   | `SUNDlsMat_SetToZero`        |
| `PrintMat`                    | `SUNDlsMat_PrintMat`         |
| `newDenseMat`                 | `SUNDlsMat_newDenseMat`      |
| `newBandMat`                  | `SUNDlsMat_newBandMat`       |
| `destroyMat`                  | `SUNDlsMat_destroyMat`       |
| `newIntArray`                 | `SUNDlsMat_newIntArray`      |
| `newIndexArray`               | `SUNDlsMat_newIndexArray`    |
| `newRealArray`                | `SUNDlsMat_newRealArray`     |
| `destroyArray`                | `SUNDlsMat_destroyArray`     |

In addition, the entire `sundials_lapack.h` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

Deprecated ARKODE nonlinear solver predictors: specification of the ARKStep
"bootstrap" or "minimum correction" predictors (options 4 and 5 from
`ARKStepSetPredictorMethod`), or MRIStep "bootstrap" predictor (option 4 from
`MRIStepSetPredictorMethod`), will output a deprecation warning message.
These options will be removed in a future release.

## Changes to SUNDIALS in release 5.8.0

### New Features

The RAJA vector implementation has been updated to support the SYCL backend in
addition to the CUDA and HIP backend. Users can choose the backend when
configuring SUNDIALS by using the `SUNDIALS_RAJA_BACKENDS` CMake variable. This
module remains experimental and is subject to change from version to version.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix and the
linear solver support general dense linear systems as well as block diagonal
linear systems. This matrix is experimental and is subject to change from
version to version.

Added a new *optional* function to the SUNLinearSolver API,
`SUNLinSolSetZeroGuess`, to indicate that the next call to `SUNLinSolSolve` will
be made with a zero initial guess. SUNLinearSolver implementations that do not
use the `SUNLinSolNewEmpty` constructor will, at a minimum, need set the
`setzeroguess` function pointer in the linear solver `ops` structure to
`NULL`. The SUNDIALS iterative linear solver implementations have been updated
to leverage this new set function to remove one dot product per solve.

The time integrator packages (ARKODE, CVODE(S), and IDA(S)) all now support a
new "matrix-embedded" SUNLinearSolver type. This type supports user-supplied
SUNLinearSolver implementations that set up and solve the specified linear
system at each linear solve call. Any matrix-related data structures are held
internally to the linear solver itself, and are not provided by the SUNDIALS
package.

Added functions to ARKODE and CVODE(S) for supplying an alternative right-hand
side function and to IDA(S) for supplying an alternative residual for use within
nonlinear system function evaluations:

* `ARKStepSetNlsRhsFn`
* `MRIStepSetNlsRhsFn`
* `CVodeSetNlsRhsFn`
* `IDASetNlsResFn`

Support for user-defined inner (fast) integrators has been to the MRIStep module
in ARKODE. See the "MRIStep Custom Inner Steppers" section in the user guide for
more information on providing a user-defined integration method.

Added specialized fused HIP kernels to CVODE which may offer better
performance on smaller problems when using CVODE with the `NVECTOR_HIP`
module. See the optional input function `CVodeSetUseIntegratorFusedKernels`
for more information. As with other SUNDIALS HIP features, this is
feature is experimental and may change from version to version.

New KINSOL options have been added to apply a constant damping factor in the
fixed point and Picard iterations (see `KINSetDamping`), to delay the start of
Anderson acceleration with the fixed point and Picard iterations (see
`KINSetDelayAA`), and to return the newest solution with the fixed point
iteration (see `KINSetReturnNewest`).

The installed `SUNDIALSConfig.cmake` file now supports the `COMPONENTS` option
to `find_package`. The exported targets no longer have IMPORTED_GLOBAL set.

### Bug Fixes

A bug was fixed in `SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

A bug was fixed in the ARKODE stepper modules where the stop time may be passed
after resetting the integrator.

A bug was fixed in `IDASetJacTimesResFn` in IDAS where the supplied function was
used in the dense finite difference Jacobian computation rather than the finite
difference Jacobian-vector product approximation.

A bug was fixed in the KINSOL Picard iteration where the value of
`KINSetMaxSetupCalls` would be ignored.

## Changes to SUNDIALS in release 5.7.0

A new NVECTOR implementation based on the SYCL abstraction layer has been added
targeting Intel GPUs. At present the only SYCL compiler supported is the DPC++
(Intel oneAPI) compiler. See the SYCL NVECTOR section in the user guide for more
details. This module is considered experimental and is subject to major changes
even in minor releases.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the MAGMA linear algebra library. Both the matrix and the linear solver
support general dense linear systems as well as block diagonal linear systems,
and both are targeted at GPUs (AMD or NVIDIA).

## Changes to SUNDIALS in release 5.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the
`CMAKE_CXX_STANDARD` and `SUNDIALS_RAJA_BACKENDS` options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

## Changes to SUNDIALS in release 5.6.0

A new `N_Vector` implementation based on the AMD ROCm HIP platform has been
added. This vector can target NVIDIA or AMD GPUs. See the HIP vector section in
the user guide for more details. This vector is considered experimental and is
subject to change from version to version.

The RAJA `N_Vector` implementation has been updated to support the HIP backend
in addition to the CUDA backend. Users can choose the backend when configuring
SUNDIALS by using the `SUNDIALS_RAJA_BACKENDS` CMake variable. This vector
remains experimental and is subject to change from version to version.

A new optional operation, `N_VGetDeviceArrayPointer`, was added to the N_Vector
API. This operation is useful for vectors that utilize dual memory spaces,
e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations
no longer require the SUNDIALS CUDA N_Vector. Instead, they require that the
vector utilized provides the `N_VGetDeviceArrayPointer` operation, and that the
pointer returned by `N_VGetDeviceArrayPointer` is a valid CUDA device pointer.

## Changes to SUNDIALS in release 5.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes
should be fully backwards compatible for almost all users. SUNDIALS
now exports CMake targets and installs a `SUNDIALSConfig.cmake` file.

Added support for SuperLU DIST 6.3.0 or newer.

## Changes to SUNDIALS in release 5.4.0

### Major Features

A new class, `SUNMemoryHelper`, was added to support **GPU users** who have
complex memory management needs such as using memory pools. This is paired with
new constructors for the `NVECTOR_CUDA` and `NVECTOR_RAJA` modules that accept a
`SUNMemoryHelper` object. Refer to "The SUNMemoryHelper API", "NVECTOR CUDA" and
"NVECTOR RAJA" sections in the documentation for more information.

Added full support for time-dependent mass matrices in ARKStep, and expanded
existing non-identity mass matrix infrastructure to support use of the
fixed point nonlinear solver.

An interface between ARKStep and the XBraid multigrid reduction in time (MGRIT)
library has been added to enable parallel-in-time integration. See the ARKStep
documentation and examples for more details. This interface required the
addition of three new N_Vector operations to exchange vector data between
computational nodes, see `N_VBufSize`, `N_VBufPack`, and `N_VBufUnpack`. These
N_Vector operations are only used within the XBraid interface and need not be
implemented for any other context.

### New Features

The `NVECTOR_RAJA` module has been updated to mirror the `NVECTOR_CUDA` module.
Notably, the update adds managed memory support to the `NVECTOR_RAJA` module.
Users of the module will need to update any calls to the `N_VMake_Raja` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

The expected behavior of `SUNNonlinSolGetNumIters` and
`SUNNonlinSolGetNumConvFails` in the SUNNonlinearSolver API have been updated to
specify that they should return the number of nonlinear solver iterations and
convergence failures in the most recent solve respectively rather than the
cumulative number of iterations and failures across all solves respectively. The
API documentation and SUNDIALS provided SUNNonlinearSolver implementations and
have been updated accordingly. As before, the cumulative number of nonlinear
iterations and failures may be retrieved by calling the integrator provided get
functions:

* `ARKStepGetNumNonlinSolvIters`
* `ARKStepGetNumNonlinSolvConvFails`
* `ARKStepGetNonlinSolvStats`
* `MRIStepGetNumNonlinSolvIters`
* `MRIStepGetNumNonlinSolvConvFails`
* `MRIStepGetNonlinSolvStats`
* `CVodeGetNumNonlinSolvIters`
* `CVodeGetNumNonlinSolvConvFails`
* `CVodeGetNonlinSolvStats`
* `IDAGetNumNonlinSolvIters`
* `IDAGetNumNonlinSolvConvFails`
* `IDAGetNonlinSolvStats`

Added the following the following functions that advanced users might find
useful when providing a custom `SUNNonlinSolSysFn`:

* `ARKStepComputeState`
* `ARKStepGetNonlinearSystemData`
* `MRIStepComputeState`
* `MRIStepGetNonlinearSystemData`
* `CVodeComputeState`
* `CVodeGetNonlinearSystemData`
* `IDAGetNonlinearSystemData`

Added new functions to CVODE(S), ARKODE, and IDA(S) to to specify the factor for
converting between integrator tolerances (WRMS norm) and linear solver tolerances
(L2 norm) i.e., `tol_L2 = nrmfac * tol_WRMS`:

* `ARKStepSetLSNormFactor`
* `ARKStepSetMassLSNormFactor`
* `MRIStepSetLSNormFactor`
* `CVodeSetLSNormFactor`
* `IDASetLSNormFactor`

Added new reset functions `ARKStepReset`, `ERKStepReset`, and
`MRIStepReset` to reset the stepper time and state vector to user-provided
values for continuing the integration from that point while retaining the
integration history. These function complement the reinitialization functions
`ARKStepReInit`, `ERKStepReInit`, and `MRIStepReInit` which reinitialize
the stepper so that the problem integration should resume as if started from
scratch.

Updated the MRIStep time-stepping module in ARKODE to support higher-order
MRI-GARK methods [Sandu, SIAM J. Numer. Anal., 57, 2019], including methods that
involve solve-decoupled, diagonally-implicit treatment of the slow time scale.

The function `CVodeSetLSetupFrequency` has been added to CVODE(S) to set
the frequency of calls to the linear solver setup function.

The Trilinos Tpetra `N_Vector` interface has been updated to work with Trilinos
12.18+. This update changes the local ordinal type to always be an `int`.

Added support for CUDA 11.

### Bug Fixes

A minor inconsistency in CVODE(S) and a bug ARKODE when checking the Jacobian
evaluation frequency has been fixed. As a result codes using using a
non-default Jacobian update frequency through a call to
`CVodeSetMaxStepsBetweenJac` or `ARKStepSetMaxStepsBetweenJac` will need to
increase the provided value by 1 to achieve the same behavior as before.

In IDAS and CVODES, the functions for forward integration with checkpointing
(`IDASolveF`, `CVodeF`) are now subject to a restriction on the number of time
steps allowed to reach the output time. This is the same restriction applied to
the `IDASolve` and `CVode` functions. The default maximum number of steps is
500, but this may be changed using the `<IDA|CVode>SetMaxNumSteps` function.
This change fixes a bug that could cause an infinite loop in the `IDASolveF`
and `CVodeF` and functions. **This change may cause a runtime error in existing user code**.

Fixed bug in using ERK method integration with static mass matrices.

### Deprecation Notice

For greater clarity the following functions have been deprecated:

* `CVodeSetMaxStepsBetweenJac`
* `ARKStepSetMaxStepsBetweenJac`
* `ARKStepSetMaxStepsBetweenLSet`

The following functions should be used instead:

* `CVodeSetJacEvalFrequency`
* `ARKStepSetJacEvalFrequency`
* `ARKStepSetLSetupFrequency`

## Changes to SUNDIALS in release 5.3.0

### Major Feature

Added support to CVODE for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user defined
projection function. This implementation is accompanied by additions to user
documentation and CVODE examples. See the `CVodeSetProjFn` function
documentation for more information.

### New Features

Added the ability to control the CUDA kernel launch parameters for the CUDA
vector and spare matrix implementations. These implementations remain
experimental and are subject to change from version to version. In addition, the
CUDA vector kernels were rewritten to be more flexible. Most users should see
equivalent performance or some improvement, but a select few may observe minor
performance degradation with the default settings. Users are encouraged to
contact the SUNDIALS team about any performance changes that they notice.

Added new capabilities for monitoring the solve phase in the Newton and
fixed-point `SUNNonlinearSolver`, and the SUNDIALS iterative linear
solvers. SUNDIALS must be built with the CMake option
`SUNDIALS_BUILD_WITH_MONITORING` to use these capabilities.

Added specialized fused CUDA kernels to CVODE which may offer better performance
on smaller problems when using CVODE with the CUDA vector. See the optional
input function `CVodeSetUseIntegratorFusedKernels` for more
information. As with other SUNDIALS CUDA features, this is feature is
experimental and may change from version to version.

Added a new function, `CVodeSetMonitorFn`, that takes a user-function
to be called by CVODE after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODE statistics
throughout the simulation.

Added a new function `CVodeGetLinSolveStats` to get the CVODE linear solver
statistics as a group.

Added the following optional functions to provide an alternative ODE right-hand
side function (ARKODE and CVODE(S)), DAE residual function (IDA(S)), or nonlinear
system function (KINSOL) for use when computing Jacobian-vector products with
the internal difference quotient approximation:

* `ARKStepSetJacTimesRhsFn`
* `CVodeSetJacTimesRhsFn`
* `CVodeSetJacTimesRhsFnB`
* `IDASetJacTimesResFn`
* `IDASetJacTimesResFnB`
* `KINSetJacTimesVecSysFn`

### Bug Fixes

Fixed a bug in the iterative linear solvers where an error is not returned if
the `Atimes` function is `NULL` or, if preconditioning is enabled, the
`PSolve` function is `NULL`.

Fixed a bug in ARKODE where the prototypes for `ERKStepSetMinReduction` and
`ARKStepSetMinReduction` were not included in `arkode_erkstep.h` and
`arkode_arkstep.h` respectively.

Fixed a bug in ARKODE where inequality constraint checking would need to be
disabled and then re-enabled to update the inequality constraint values after
resizing a problem. Resizing a problem will now disable constraints and a call
to `ARKStepSetConstraints` or `ERKStepSetConstraints` is required to re-enable
constraint checking for the new problem size.

## Changes to SUNDIALS in release 5.2.0

### New Features

The following functions were added to each of the time integration packages to
enable or disable the scaling applied to linear system solutions with
matrix-based linear solvers to account for lagged matrix information:

* `ARKStepSetLinearSolutionScaling`
* `CVodeSetLinearSolutionScaling`
* `CVodeSetLinearSolutionScalingB`
* `IDASetLinearSolutionScaling`
* `IDASetLinearSolutionScalingB`

When using a matrix-based linear solver with ARKODE, IDA(S), or BDF methods in
CVODE(S) scaling is enabled by default.

Added a new `SUNMatrix` implementation that interfaces to the sparse matrix
implementation from the NVIDIA cuSPARSE library. In addition, the CUDA Sparse
linear solver has been updated to use the new matrix, as such, users of this
matrix will need to update their code. This implementations are still considered
to be experimental, thus they are subject to breaking changes even in minor
releases.

Added a new "stiff" interpolation module to ARKODE, based on Lagrange polynomial
interpolation, that is accessible to each of the ARKStep, ERKStep and MRIStep
time-stepping modules. This module is designed to provide increased
interpolation accuracy when integrating stiff problems, as opposed to the ARKODE
standard Hermite interpolation module that can suffer when the IVP right-hand
side has large Lipschitz constant. While the Hermite module remains the default,
the new Lagrange module may be enabled using one of the routines
`ARKStepSetInterpolantType`, `ERKStepSetInterpolantType`, or
`MRIStepSetInterpolantType`. The serial example problem `ark_brusselator.c` has
been converted to use this Lagrange interpolation module. Created accompanying
routines `ARKStepSetInterpolantDegree`, `ARKStepSetInterpolantDegree` and
`ARKStepSetInterpolantDegree` to provide user control over these interpolating
polynomials.

Added two new functions, `ARKStepSetMinReduction` and `ERKStepSetMinReduction`
to change the minimum allowed step size reduction factor after an error test
failure.

### Bug Fixes

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set `CMAKE_Fortran_COMPILER` to `f2003`, `xlf2003`, or
`xlf2003_r`.

Fixed a bug in how ARKODE interfaces with a user-supplied, iterative, unscaled
linear solver. In this case, ARKODE adjusts the linear solver tolerance in an
attempt to account for the lack of support for left/right scaling matrices.
Previously, ARKODE computed this scaling factor using the error weight vector,
`ewt`; this fix changes that to the residual weight vector, `rwt`, that can
differ from `ewt` when solving problems with non-identity mass matrix.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attribute missing on some SUNDIALS API functions.

Fixed a memory leak in CVODES and IDAS from not deallocating the `atolSmin0` and
`atolQSmin0` arrays.

Fixed a bug where a non-default value for the maximum allowed growth factor
after the first step would be ignored.

### Deprecation Notice

The routines `ARKStepSetDenseOrder`, `ARKStepSetDenseOrder` and
`ARKStepSetDenseOrder` have been deprecated and will be removed in a
future release. The new functions `ARKStepSetInterpolantDegree`,
`ARKStepSetInterpolantDegree`, and `ARKStepSetInterpolantDegree`
should be used instead.

## Changes to SUNDIALS in release 5.1.0

### New Features

Added support for a user-supplied function to update the prediction for each
implicit stage solution in ARKStep. If supplied, this routine will be called
*after* any existing ARKStep predictor algorithm completes, so that the
predictor may be modified by the user as desired. The new user-supplied routine
has type `ARKStepStagePredictFn`, and may be set by calling
`ARKStepSetStagePredictFn`.

The MRIStep module has been updated to support attaching different user data
pointers to the inner and outer integrators. If applicable, user codes will
need to add a call to `ARKStepSetUserData` to attach their user data
pointer to the inner integrator memory as `MRIStepSetUserData` will
not set the pointer for both the inner and outer integrators. The MRIStep
examples have been updated to reflect this change.

Added support for damping when using Anderson acceleration in KINSOL. See the
mathematical considerations section of the user guide and the description of the
`KINSetDampingAA` function for more details.

Added support for constant damping to the `SUNNonlinearSolver_FixedPoint` module
when using Anderson acceleration. See the `SUNNonlinearSolver_FixedPoint`
section in the user guides and the description of the
`SUNNonlinSolSetDamping_FixedPoint` function for more details.

Added two utility functions, `SUNDIALSFileOpen` and `SUNDIALSFileClose` for
creating/destroying file pointers. These are useful when using the Fortran 2003
interfaces.

Added a new build system option, `CUDA_ARCH`, to specify the CUDA architecture
to target.

### Bug Fixes

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables `PETSC_INCLUDES` and `PETSC_LIBRARIES` instead of `PETSC_DIR`.

Fixed a bug in the Fortran 2003 interfaces to the ARKODE Butcher table routines
and structure. This includes changing the `ARKodeButcherTable` type to be a
`type(c_ptr)` in Fortran.

## Changes to SUNDIALS in release 5.0.0

### Build System

Increased the minimum required CMake version to 3.5 for most SUNDIALS
configurations, and 3.10 when CUDA or OpenMP with device offloading are enabled.

The CMake option `BLAS_ENABLE` and the variable `BLAS_LIBRARIES` have been
removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
third party libraries that require linking to BLAS, the path to the BLAS
library should be included in the `_LIBRARIES` variable for the third party
library e.g., `SUPERLUDIST_LIBRARIES` when enabling SuperLU_DIST.

### NVector

Two new functions were added to aid in creating custom `N_Vector`
objects. The constructor `N_VNewEmpty` allocates an "empty" generic
`N_Vector` with the object's content pointer and the function pointers
in the operations structure initialized to `NULL`. When used in the
constructor for custom objects this function will ease the introduction of any
new optional operations to the `N_Vector` API by ensuring only required
operations need to be set. Additionally, the function `N_VCopyOps` has
been added to copy the operation function pointers between vector objects. When
used in clone routines for custom vector objects these functions also will ease
the introduction of any new optional operations to the `N_Vector` API by
ensuring all operations are copied when cloning objects.

Added new `N_Vector` implementations, `ManyVector` and `MPIManyVector`, to
support flexible partitioning of solution data among different processing
elements (e.g., CPU + GPU) or for multi-physics problems that couple distinct
MPI-based simulations together (see the the `ManyVector` and `MPIManyVector`
section in the user guide for more details). This implementation is accompanied
by additions to user documentation and SUNDIALS examples.

Additionally, an `MPIPlusX` vector implementation has been created to support
the MPI+X paradigm where X is a type of on-node parallelism (e.g., OpenMP, CUDA,
etc.). The implementation is accompanied by additions to user documentation and
SUNDIALS examples.

One new required vector operation and ten new optional vector operations have
been added to the `N_Vector` API. The new required operation, `N_VGetLength`,
returns the global vector length. The optional operations have been added to
support the new MPIManyVector implementation. The operation `N_VGetCommunicator`
must be implemented by subvectors that are combined to create an MPIManyVector,
but is not used outside of this context. The remaining nine operations are
optional local reduction operations intended to eliminate unnecessary latency
when performing vector reduction operations (norms, etc.) on distributed memory
systems. The optional local reduction vector operations are `N_VDotProdLocal`,
`N_VMaxNormLocal`, `N_VMinLocal`, `N_VL1NormLocal`, `N_VWSqrSumLocal`,
`N_VWSqrSumMaskLocal`, `N_VInvTestLocal`, `N_VConstrMaskLocal`, and
`N_VMinQuotientLocal`. If an `N_Vector` implementation defines any of the local
operations as `NULL`, then the MPIManyVector will call standard `N_Vector`
operations to complete the computation.

The `*_MPICuda` and `*_MPIRaja` functions have been removed from the CUDA
and RAJA vector implementations respectively. Accordingly, the
`nvector_mpicuda.h`, `nvector_mpiraja.h`, `libsundials_nvecmpicuda.lib`,
and `libsundials_nvecmpicudaraja.lib` files have been removed. Users should
use the MPI+X vector in conjunction with the CUDA and RAJA vectors to replace
the functionality. The necessary changes are minimal and should require few code
modifications. See the example programs in `examples/ida/mpicuda` and
`examples/ida/mpiraja` for examples of how to use the MPI+X vector with the
CUDA and RAJA vectors, respectively.

Made performance improvements to the CUDA vector. Users who utilize a
non-default stream should no longer see default stream synchronizations after
memory transfers.

Added a new constructor to the CUDA vector that allows a user to provide custom
allocate and free functions for the vector data array and internal reduction
buffer.

Added three new `N_Vector` utility functions, `N_VGetVecAtIndexVectorArray`,
`N_VSetVecAtIndexVectorArray`, and `N_VNewVectorArray`, for working with
`N_Vector` arrays when using the Fortran 2003 interfaces.

### SUNMatrix

Two new functions were added to aid in creating custom SUNMatrix objects. The
constructor `SUNMatNewEmpty` allocates an "empty" generic SUNMatrix with the
object's content pointer and the function pointers in the operations structure
initialized to `NULL`. When used in the constructor for custom objects this
function will ease the introduction of any new optional operations to the
SUNMatrix API by ensuring only required operations need to be set. Additionally,
the function `SUNMatCopyOps(A, B)` has been added to copy the operation function
pointers between matrix objects. When used in clone routines for custom matrix
objects these functions also will ease the introduction of any new optional
operations to the SUNMatrix API by ensuring all operations are copied when
cloning objects.

A new operation, `SUNMatMatvecSetup`, was added to the `SUNMatrix` API to
perform any setup necessary for computing a matrix-vector product. This
operation is useful for `SUNMatrix` implementations which need to prepare the
matrix itself, or communication structures before performing the matrix-vector
product. Users who have implemented a custom `SUNMatrix` will need to at least
update their code to set the corresponding `ops` structure member,
`matvecsetup`, to `NULL`.

The generic SUNMatrix API now defines error codes to be returned by SUNMatrix
operations. Operations which return an integer flag indicating success/failure
may return different values than previously.

A new SUNMatrix (and SUNLinearSolver) implementation was added to facilitate
the use of the SuperLU_DIST library with SUNDIALS.

### SUNLinearSolver

A new function was added to aid in creating custom `SUNLinearSolver`
objects. The constructor `SUNLinSolNewEmpty` allocates an "empty" generic
`SUNLinearSolver` with the object's content pointer and the function pointers in
the operations structure initialized to `NULL`. When used in the constructor for
custom objects this function will ease the introduction of any new optional
operations to the `SUNLinearSolver` API by ensuring only required operations
need to be set.

The return type of the `SUNLinSolLastFlag` in the `SUNLinearSolver` has changed
from `long int` to `sunindextype` to be consistent with the type used to store
row indices in dense and banded linear solver modules.

Added a new optional operation to the SUNLINEARSOLVER API, `SUNLinSolGetID`,
that returns a `SUNLinearSolver_ID` for identifying the linear solver module.

The SUNLinearSolver API has been updated to make the initialize and setup
functions optional.

A new SUNLinearSolver (and SUNMatrix) implementation was added to facilitate
the use of the SuperLU_DIST library with SUNDIALS.

Added a new SUNLinearSolver implementation,
`SUNLinearSolver_cuSolverSp_batchQR`, which leverages the NVIDIA cuSOLVER sparse
batched QR method for efficiently solving block diagonal linear systems on
NVIDIA GPUs.

Added three new accessor functions to the SUNLinSol_KLU module,
`SUNLinSol_KLUGetSymbolic`, `SUNLinSol_KLUGetNumeric`, and
`SUNLinSol_KLUGetCommon`, to provide user access to the underlying
KLU solver structures.

### SUNNonlinearSolver

A new function was added to aid in creating custom `SUNNonlinearSolver`
objects. The constructor `SUNNonlinSolNewEmpty` allocates an "empty" generic
`SUNNonlinearSolver` with the object's content pointer and the function pointers
in the operations structure initialized to `NULL`. When used in the constructor
for custom objects this function will ease the introduction of any new optional
operations to the `SUNNonlinearSolver` API by ensuring only required operations
need to be set.

To facilitate the use of user supplied nonlinear solver convergence test
functions the `SUNNonlinSolSetConvTestFn` function in the SUNNonlinearSolver API
has been updated to take a `void*` data pointer as input. The supplied data
pointer will be passed to the nonlinear solver convergence test function on each
call.

The inputs values passed to the first two inputs of the `SUNNonlinSolSolve`
function in the `SUNNonlinearSolver` have been changed to be the predicted state
and the initial guess for the correction to that state. Additionally, the
definitions of `SUNNonlinSolLSetupFn` and `SUNNonlinSolLSolveFn` in the
SUNNonlinearSolver API have been updated to remove unused input parameters. For
more information on the nonlinear system formulation and the API functions see
the `SUNNonlinearSolver` chapter in the user guides.

Added a new `SUNNonlinearSolver` implementation for interfaces to the PETSc SNES
nonlinear solver.

### New Features

A new linear solver interface functions, `ARKLsLinSysFn` and `CVLsLinSysFn`, as
added as an alternative method for evaluating the linear systems `M - \gamma J`
or `I - \gamma J`.

Added the following functions to get the current state and gamma value to
ARKStep, CVODE and CVODES that may be useful to users who choose to provide
their own nonlinear solver implementation:

* `ARKStepGetCurrentState`
* `ARKStepGetCurrentGamma`
* `CVodeGetCurrentGamma`
* `CVodeGetCurrentState`
* `CVodeGetCurrentGamma`
* `CVodeGetCurrentStateSens`
* `CVodeGetCurrentSensSolveIndex`
* `IDAGetCurrentCj`
* `IDAGetCurrentY`
* `IDAGetCurrentYp`
* `IDAComputeY`
* `IDAComputeYp`

Removed extraneous calls to `N_VMin` for simulations where the scalar
valued absolute tolerance, or all entries of the vector-valued absolute
tolerance array, are strictly positive. In this scenario ARKODE, CVODE(S), and
IDA(S) steppers will remove at least one global reduction per time step.

The ARKODE, CVODE(S), IDA(S), and KINSOL linear solver interfaces have been
updated to only zero the Jacobian matrix before calling a user-supplied Jacobian
evaluation function when the attached linear solver has type
`SUNLINEARSOLVER_DIRECT`.

Added new Fortran 2003 interfaces to all of the SUNDIALS packages (ARKODE,
CVODE(S), IDA(S), and KINSOL as well as most of the `N_Vector`, `SUNMatrix`,
`SUNLinearSolver`, and `SUNNonlinearSolver` implementations. See "Fortran"
section for more details. These new interfaces were generated with SWIG-Fortran
and provide a user an idiomatic Fortran 2003 interface to most of the SUNDIALS C
API.

The MRIStep module has been updated to support explicit, implicit, or IMEX
methods as the fast integrator using the ARKStep module. As a result some
function signatures have been changed including MRIStepCreate which now
takes an ARKStep memory structure for the fast integration as an input.

The reinitialization functions `ERKStepReInit`, `ARKStepReInit`, and
`MRIStepReInit` have been updated to retain the minimum and maximum step
size values from before reinitialization rather than resetting them to the
default values.

Added two new embedded ARK methods of orders 4 and 5 to ARKODE (from
Kennedy & Carpenter, Appl. Numer. Math., 136:183-205, 2019).

Support for optional inequality constraints on individual components of the
solution vector has been added the ARKODE ERKStep and ARKStep modules. See
the descriptions of `ERKStepSetConstraints` and `ARKStepSetConstraints` for
more details. Note that enabling constraint handling requires the `N_Vector`
operations `N_VMinQuotient`, `N_VConstrMask`, and `N_VCompare` that were not
previously required by ARKODE.

Add two new 'Set' functions to MRIStep, `MRIStepSetPreInnerFn` and
`MRIStepSetPostInnerFn`, for performing communication or memory
transfers needed before or after the inner integration.

### Bug Fixes

Fixed a bug in the build system that prevented the PThreads NVECTOR module from
being built.

Fixed a memory leak in the PETSc `N_Vector` clone function.

Fixed a memory leak in the ARKODE, CVODE, and IDA F77 interfaces when not using
the default nonlinear solver.

Fixed a bug in the ARKStep time-stepping module in ARKODE that would result in
an infinite loop if the nonlinear solver failed to converge more than the
maximum allowed times during a single step.

Fixed a bug in ARKODE that would result in a "too much accuracy requested" error
when using fixed time step sizes with explicit methods in some cases.

Fixed a bug in ARKStep where the mass matrix linear solver setup function was
not called in the Matrix-free case.

Fixed a minor bug in ARKStep where an incorrect flag is reported when an
error occurs in the mass matrix setup or Jacobian-vector product setup
functions.

Fixed a bug in the CVODE and CVODES constraint handling where the step size
could be set below the minimum step size.

Fixed a bug in the CVODE and CVODES nonlinear solver interfaces where the norm
of the accumulated correction was not updated when using a non-default
convergence test function.

Fixed a bug in the CVODES `cvRescale` function where the loops to compute the
array of scalars for the fused vector scale operation stopped one iteration
early.

Fixed a bug in CVODES and IDAS where `CVodeF` and `IDASolveF` would return the
wrong flag under certain  circumstances.

Fixed a bug in CVODES and IDAS where `CVodeF` and `IDASolveF` would not return a
root in `NORMAL_STEP` mode if the root occurred after the desired output time.

Fixed a bug in the IDA and IDAS linear solver interfaces where an incorrect
Jacobian-vector product increment was used with iterative solvers other than
SPGMR and SPFGMR.

Fixed a bug the IDAS `IDAQuadReInitB` function where an incorrect memory
structure was passed to `IDAQuadReInit`.

Fixed a bug in the KINSOL linear solver interface where the auxiliary scalar
`sJpnorm` was not computed when necessary with the Picard iteration and the
auxiliary scalar `sFdotJp` was unnecessarily computed in some cases.

## Changes to SUNDIALS in release 4.1.0

### Removed Implementation Headers

The implementation header files (`*_impl.h`) are no longer installed. This
means users who are directly accessing or manipulating package memory structures
will need to update their code to use the package's public API.

### New Features

An additional `N_Vector` implementation was added for interfacing with
the Tpetra vector from Trilinos library to facilitate interoperability between
SUNDIALS and Trilinos. This implementation is accompanied by additions to user
documentation and SUNDIALS examples.

### Bug Fixes

The `EXAMPLES_ENABLE_RAJA` CMake option has been removed. The option
`EXAMPLES_ENABLE_CUDA` enables all examples that use CUDA including the RAJA
examples with a CUDA back end (if RAJA is enabled).

Python is no longer required to run `make test` and `make test_install`.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

Fixed a bug in `ARKodeButcherTable_Write` when printing a Butcher table without
an embedding.

## Changes to SUNDIALS in release 4.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved the definitions of backwards compatibility functions for the prior direct
linear solver (DLS) and scaled preconditioned iterarive linear solvers (SPILS)
to a source file. The symbols are now included in the appropriate package
library, e.g. `libsundials_cvode.lib`.

## Changes to SUNDIALS in release 4.0.1

A bug in ARKODE where single precision builds would fail to compile has been
fixed.

## Changes to SUNDIALS in release 4.0.0

The direct and iterative linear solver interfaces in all SUNDIALS packages have
been merged into a single unified linear solver interface to support any valid
`SUNLinearSolver`. This includes the `DIRECT` and `ITERATIVE` types
as well as the new `MATRIX_ITERATIVE` type. Details regarding how SUNDIALS
packages utilize linear solvers of each type as well as a discussion regarding
the intended use cases for user-supplied linear solver implementations are
included the user guide. All example programs have been updated to use
the new unified linear solver interfaces.

The unified linear solver interface is very similar to the previous DLS (direct
linear solver) and SPILS (scaled preconditioned iterative linear solver)
interface in each package. To minimize challenges in user migration to the
unified linear solver interfaces, the previous DLS and SPILS functions may still
be used however, these are now deprecated and will be removed in a future
release. Additionally, that Fortran users will need to enlarge their array of
optional integer outputs, and update the indices that they query for certain
linear solver related statistics.

The names of all SUNDIALS-provided `SUNLinearSolver` constructors have have been
updated to follow the naming convention `SUNLinSol_*` where `*` is the name
of the linear solver. The new constructor names are:

* `SUNLinSol_Band`
* `SUNLinSol_Dense`
* `SUNLinSol_KLU`
* `SUNLinSol_LapackBand`
* `SUNLinSol_LapackDense`
* `SUNLinSol_PCG`
* `SUNLinSol_SPBCGS`
* `SUNLinSol_SPFGMR`
* `SUNLinSol_SPGMR`
* `SUNLinSol_SPTFQMR`
* `SUNLinSol_SuperLUMT`

Linear solver-specific "set" routine names have been similarly standardized. To
minimize challenges in user migration to the new names, the previous function
names may still be used however, these are now deprecated and will be removed in
a future release. All example programs and the standalone linear solver examples
have been updated to use the new naming convention.

The `SUNLinSol_Band` constructor has been simplified to remove the
storage upper bandwidth argument.

SUNDIALS integrators (ARKODE, CVODE(S), and IDA(S)) have been updated to utilize
generic nonlinear solvers defined by the `SUNNonlinearSolver` API. This enables
the addition of new nonlinear solver options and allows for external or
user-supplied nonlinear solvers. The nonlinear solver API and SUNDIALS provided
implementations are described in the user guide and follow the same object
oriented design used by the `N_Vector`, `SUNMatrix`, and `SUNLinearSolver`
classes. Currently two nonlinear solver implementations are provided, Newton and
fixed-point. These replicate the previous integrator-specific implementations of
Newton's method and a fixed-point iteration (previously referred to as a
functional iteration), respectively. Note the new fixed-point implementation can
optionally utilize Anderson's method to accelerate convergence. Example programs
using each of these nonlinear solvers in a standalone manner have been added and
all example programs have been updated accordingly.

The SUNDIALS integrators (ARKODE, CVODE(S), and IDA(S)) all now use the Newton
`SUNNonlinearSolver` by default. Users that wish to use the fixed-point
`SUNNonlinearSolver` will need to create the corresponding nonlinear solver
object and attach it to the integrator with the appropriate set function:

* `ARKStepSetNonlinearSolver`
* `CVodeSetNonlinearSolver`
* `IDASetNonlinearSolver`

Functions for setting the nonlinear solver options or getting nonlinear solver
statistics remain unchanged and internally call generic `SUNNonlinearSolver`
functions as needed.

With the introduction of the `SUNNonlinearSolver` class, the input parameter
`iter` to `CVodeCreate` has been removed along with the function
`CVodeSetIterType` and the constants `CV_NEWTON` and `CV_FUNCTIONAL`. While
SUNDIALS includes a fixed-point nonlinear solver, it is not currently supported
in IDA.

Three fused vector operations and seven vector array operations have been added
to the `N_Vector` API. These *optional* operations are disabled by default and
may be activated by calling vector specific routines after creating a vector
(see the `N_Vector` chapter for more details). The new operations are intended
to increase data reuse in vector operations, reduce parallel communication on
distributed memory systems, and lower the number of kernel launches on systems
with accelerators. The fused operations are:

* `N_VLinearCombination`
* `N_VScaleAddMulti`
* `N_VDotProdMulti`

and the vector array operations are:

* `N_VLinearCombinationVectorArray`
* `N_VScaleVectorArray`
* `N_VConstVectorArray`
* `N_VWrmsNormVectorArray`
* `N_VWrmsNormMaskVectorArray`
* `N_VScaleAddMultiVectorArray`
* `N_VLinearCombinationVectorArray`

If an `N_Vector` implementation defines the implementation any of these
operations as `NULL`, then standard vector operations will automatically be
called as necessary to complete the computation.

A new `N_Vector` implementation, `OpenMPDEV`, leveraging OpenMP device
offloading has been added.

Multiple updates to the CUDA vector were made:

* Changed the `N_VMake_Cuda` function to take a host data pointer and a device
  data pointer instead of an `N_VectorContent_Cuda` object.

* Changed `N_VGetLength_Cuda` to return the global vector length instead of
  the local vector length.

* Added `N_VGetLocalLength_Cuda` to return the local vector length.

* Added `N_VGetMPIComm_Cuda` to return the MPI communicator used.

* Removed the accessor functions in the `suncudavec` namespace.

* Added the ability to set the `cudaStream_t` used for execution of the CUDA
  vector kernels. See the function `N_VSetCudaStreams_Cuda`.

* Added `N_VNewManaged_Cuda`, `N_VMakeManaged_Cuda`, and
  `N_VIsManagedMemory_Cuda` functions to accommodate using managed memory with
  the CUDA vector.

Multiple updates to the RAJA vector were made:

* Changed `N_VGetLength_Raja` to return the global vector length instead of
  the local vector length.

* Added `N_VGetLocalLength_Raja` to return the local vector length.

* Added `N_VGetMPIComm_Raja` to return the MPI communicator used.

* Removed the accessor functions in the `sunrajavec` namespace.

Two changes were made in the ARKODE and CVODE(S) initial step size algorithm:

* Fixed an efficiency bug where an extra call to the RHS function was made.

* Changed the behavior of the algorithm if the max-iterations case is hit.
  Before the algorithm would exit with the step size calculated on the
  penultimate iteration. Now it will exit with the step size calculated
  on the final iteration.

Fortran 2003 interfaces to CVODE, the fixed-point and Newton nonlinear solvers,
the dense, band, KLU, PCG, SPBCGS, SPFGMR, SPGMR, and SPTFQMR linear solvers,
and the serial, PThreads, and OpenMP vectors have been added.

The ARKODE library has been entirely rewritten to support a modular approach to
one-step methods, which should allow rapid research and development of novel
integration methods without affecting existing solver functionality. To support
this, the existing ARK-based methods have been encapsulated inside the new
`ARKStep` time-stepping module. Two new time-stepping modules have been added:

* The `ERKStep` module provides an optimized implementation for explicit
  Runge-Kutta methods with reduced storage and number of calls to the ODE
  right-hand side function.

* The `MRIStep` module implements two-rate explicit-explicit multirate
  infinitesimal step methods utilizing different step sizes for slow and fast
  processes in an additive splitting.

This restructure has resulted in numerous small changes to the user interface,
particularly the suite of "Set" routines for user-provided solver parameters and
"Get" routines to access solver statistics, that are now prefixed with the name
of time-stepping module (e.g., `ARKStep` or `ERKStep`) instead of
`ARKODE`. Aside from affecting the names of these routines, user-level changes
have been kept to a minimum. However, we recommend that users consult both this
documentation and the ARKODE example programs for further details on the updated
infrastructure.

As part of the ARKODE restructuring an `ARKodeButcherTable` structure
has been added for storing Butcher tables. Functions for creating new Butcher
tables and checking their analytic order are provided along with other utility
routines. For more details see the Butcher Table section in the user guide.

ARKODE's dense output infrastructure has been improved to support higher-degree
Hermite polynomial interpolants (up to degree 5) over the last successful time
step.

## Changes to SUNDIALS in release 3.2.1

Fixed a bug in the CUDA vector where the `N_VInvTest` operation could write
beyond the allocated vector data.

Fixed the library installation path for multiarch systems. This fix changes the
default library installation path from `CMAKE_INSTALL_PREFIX/lib` to
`CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`. The default value library directory
name is automatically set to `lib`, `lib64`, or `lib/<multiarch-tuple>`
depending on the system, but maybe be overridden by setting
`CMAKE_INSTALL_LIBDIR`.

## Changes to SUNDIALS in release 3.2.0

### Library Name Change

Changed the name of the RAJA nvector library to `libsundials_nveccudaraja.lib`
from `libsundials_nvecraja.lib` to better reflect that we only support CUDA as a
backend for RAJA currently.

### New Features

Added hybrid MPI+CUDA and MPI+RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Support for optional inequality constraints on individual components of the
solution vector has been added to CVODE and CVODES. For more details see the
Mathematical Considerations and Optional Input sections of the user guide. Use
of `CVodeSetConstraints` requires the `N_Vector` operations `N_VMinQuotient`,
`N_VConstrMask`, and `N_VCompare` that were not previously required by CVODE and
CVODES.

### CMake Updates

CMake 3.1.3 is now the minimum required CMake version.

Deprecated the behavior of the `SUNDIALS_INDEX_TYPE` CMake option and added the
`SUNDIALS_INDEX_SIZE` CMake option to select the `sunindextype` integer size.

The native CMake FindMPI module is now used to locate an MPI installation.

If MPI is enabled and MPI compiler wrappers are not set, the build system will
check if `CMAKE_<language>_COMPILER` can compile MPI programs before trying to
locate and use an MPI installation.

The previous options for setting MPI compiler wrappers and the executable for
running MPI programs have been have been deprecated. The new options that align
with those used in native CMake FindMPI module are `MPI_C_COMPILER`,
`MPI_CXX_COMPILER`, `MPI_Fortran_COMPILER`, and `MPIEXEC_EXECUTABLE`.

When a Fortran name-mangling scheme is needed (e.g., `ENABLE_LAPACK` is `ON`)
the build system will infer the scheme from the Fortran compiler. If a Fortran
compiler is not available or the inferred or default scheme needs to be
overridden, the advanced options `SUNDIALS_F77_FUNC_CASE` and
`SUNDIALS_F77_FUNC_UNDERSCORES` can be used to manually set the name-mangling
scheme and bypass trying to infer the scheme.

Parts of the main `CMakeLists.txt` file were moved to new files in the `src` and
`example` directories to make the CMake configuration file structure more
modular.

### Bug Fixes

Fixed a problem with setting `sunindextype` which would occur with some
compilers (e.g. `armclang`) that do not define `__STDC_VERSION__`.

Fixed a thread-safety issue in CVODES and IDAS when using adjoint sensitivity
analysis.

Fixed a bug in IDAS where the saved residual value used in the nonlinear solve
for consistent initial conditions was passed as temporary workspace and could be
overwritten.

## Changes to SUNDIALS in release 3.1.2

### CMake Updates

Updated the minimum required version of CMake to 2.8.12 and enabled using rpath
by default to locate shared libraries on OSX.

### New Features

Added the function `SUNSparseMatrix_Reallocate` to allow specification of the
matrix nonzero storage.

Added named constants for the two reinitialization types for the KLU
SUNLinearSolver.

Updated the `SUNMatScaleAdd` and `SUNMatScaleAddI` implementations in the sparse
SUNMatrix to more optimally handle the case where the target matrix contained
sufficient storage for the sum, but had the wrong sparsity pattern. The sum now
occurs in-place, by performing the sum backwards in the existing
storage. However, it is still more efficient if the user-supplied Jacobian
routine allocates storage for the sum `M + gamma J` or `M + gamma J` manually
(with zero entries if needed).

The following examples from the usage notes page of the SUNDIALS website, and
updated them to work with SUNDIALS 3.x:

* `cvDisc_dns.c` demonstrates using CVODE with discontinuous solutions or RHS.

* `cvRoberts_dns_negsol.c` illustrates the use of the RHS function return
  value to control unphysical negative concentrations.

* `cvRoberts_FSA_dns_Switch.c` demonstrates switching on/off forward
  sensitivity computations. This example came from the usage notes page of the
  SUNDIALS website.

### Bug Fixes

Fixed a Windows specific problem where `sunindextype` was not correctly defined
when using 64-bit integers. On Windows `sunindextype` is now defined as the MSVC
basic type `__int64`.

Fixed a bug in the full KLU SUNLinearSolver reinitialization approach where the
sparse SUNMatrix pointer would go out of scope on some architectures.

The misnamed function `CVSpilsSetJacTimesSetupFnBS` has been deprecated and
replaced by `CVSpilsSetJacTimesBS`. The deprecated function
`CVSpilsSetJacTimesSetupFnBS` will be removed in the next major release.

Changed LICENSE install path to `instdir/include/sundials`.

## Changes to SUNDIALS in release 3.1.1

### Bug Fixes

Fixed a minor bug in the CVODE and CVODES `cvSLdet` routine, where a return was
missing in the error check for three inconsistent roots.

Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers: if
"Initialize" was called multiple times then the solver memory was reallocated
(without being freed).

Fixed a minor bug in `ARKReInit`, where a flag was incorrectly set to indicate
that the problem had been resized (instead of just re-initialized).

Fixed C++11 compiler errors/warnings about incompatible use of string literals.

Updated the KLU SUNLinearSolver to use a typedef for the precision-specific
solve functions to avoid compiler warnings.

Added missing typecasts for some (`void*`) pointers to avoid compiler warnings.

Fixed bug in the sparse SUNMatrix where `int` was used instead of
`sunindextype` in one location.

Fixed a minor bug in `KINPrintInfo` where a case was missing for
`KIN_REPTD_SYSFUNC_ERR` leading to an undefined info message.

Added missing `#include <stdio.h>` in `N_Vector` and `SUNMatrix` header files.

Added missing prototypes for `ARKSpilsGetNumMTSetups` in ARKODE and
`IDASpilsGetNumJTSetupEvals` in IDA and IDAS.

Fixed an indexing bug in the CUDA vector implementation of `N_VWrmsNormMask` and
revised the RAJA vector implementation of `N_VWrmsNormMask` to work with mask
arrays using values other than zero or one. Replaced `double` with `realtype` in
the RAJA vector test functions.

Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not require
a `SUNMatrix` or `SUNLinearSolver` e.g., iterative linear solvers, explicit
methods in ARKODE, functional iteration in CVODE, etc.

## Changes to SUNDIALS in release 3.1.0

Added `N_Vector` print functions that write vector data to a specified file
(e.g., `N_VPrintFile_Serial`).

Added `make test` and `make test_install` options to the build system for
testing SUNDIALS after building with `make` and installing with `make install`
respectively.

## Changes to SUNDIALS in release 3.0.0

### Major Feature

Added new linear solver and matrix interfaces for all SUNDIALS packages and
updated the existing linear solver and matrix implementations. The goal of the
redesign is to provide greater encapsulation and ease interfacing custom linear
solvers with linear solver libraries. Specific changes include:

* Added a `SUNMatrix` interface with three provided implementations:
  dense, banded, and sparse. These replicate previous SUNDIALS direct (Dls) and
  sparse (Sls) matrix structures.

* Added example problems demonstrating use of the matrices.

* Added a `SUNLinearSolver` interface with eleven provided implementations:
  dense, banded, LAPACK dense, LAPACK band, KLU, SuperLU_MT, SPGMR, SPBCGS,
  SPTFQMR, SPFGMR, PCG. These replicate previous SUNDIALS generic linear
  solvers.

* Added example problems demonstrating use of the linear solvers.

* Expanded package-provided direct linear solver (Dls) interfaces and scaled,
  preconditioned, iterative linear solver (Spils) interfaces to utilize
  `SUNMatrix` and `SUNLinearSolver` objects.

* Removed package-specific, linear solver-specific, solver modules (e.g.,
  CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since their functionality is entirely
  replicated by the generic Dls/Spils interfaces and `SUNLinearSolver` /
  `SUNMatrix` classes. The exception is `CVDIAG`, a diagonal
  approximate Jacobian solver available to CVODE and CVODES.

* Converted all SUNDIALS example problems to utilize new the new matrix and
  linear solver objects, along with updated Dls and Spils linear solver
  interfaces.

* Added Spils interface routines to ARKODE, CVODE, CVODES, IDA and IDAS to allow
  specification of a user-provided `JTSetup` routine. This change supports
  users who wish to set up data structures for the user-provided
  Jacobian-times-vector (`JTimes`) routine, and where the cost of one
  `JTSetup` setup per Newton iteration can be amortized between multiple
  `JTimes` calls.

Corresponding updates were made to all the example programs.

### New Features

CUDA and RAJA `N_Vector` implementations to support GPU systems. These vectors
are supplied to provide very basic support for running on GPU architectures.
Users are advised that these vectors both move all data to the GPU device upon
construction, and speedup will only be realized if the user also conducts the
right-hand-side function evaluation on the device. In addition, these vectors
assume the problem fits on one GPU. For further information about RAJA, users
are referred to the [RAJA web site](https://software.llnl.gov/RAJA/).

Added the type `sunindextype` to support using 32-bit or 64-bit integer types
for indexing arrays within all SUNDIALS structures. `sunindextype` is defined to
`int32_t` or `int64_t` when portable types are supported, otherwise it is
defined as `int` or `long int`. The Fortran interfaces continue to use `long
int` for indices, except for the sparse matrix interface that now uses
`sunindextype`. Interfaces to PETSc, hypre, SuperLU_MT, and KLU have been
updated with 32-bit or 64-bit capabilities depending how the user configures
SUNDIALS.

To avoid potential namespace conflicts, the macros defining `booleantype` values
`TRUE` and `FALSE` have been changed to `SUNTRUE` and `SUNFALSE` respectively.

Temporary vectors were removed from preconditioner setup and solve routines for
all packages. It is assumed that all necessary data for user-provided
preconditioner operations will be allocated and stored in user-provided data
structures.

The file `include/sundials_fconfig.h` was added. This file contains SUNDIALS
type information for use in Fortran programs.

Added support for many xSDK-compliant build system keys. For more information on
on xSDK compliance the [xSDK policies](https://xsdk.info/policies/). The xSDK
is a movement in scientific software to provide a foundation for the rapid and
efficient production of high-quality, sustainable extreme-scale scientific
applications. For more information visit the
[xSDK web site](https://xsdk.info).

Added functions `SUNDIALSGetVersion` and `SUNDIALSGetVersionNumber` to get
SUNDIALS release version information at runtime.

Added comments to `arkode_butcher.c` regarding which methods should have
coefficients accurate enough for use in quad precision.

### Build System

Renamed CMake options to enable/disable examples for greater clarity and added
option to enable/disable Fortran 77 examples:

* Changed `EXAMPLES_ENABLE` to `EXAMPLES_ENABLE_C`
* Changed `CXX_ENABLE` to `EXAMPLES_ENABLE_CXX`
* Changed `F90_ENABLE` to `EXAMPLES_ENABLE_F90`
* Added `EXAMPLES_ENABLE_F77` option

Added separate `BLAS_ENABLE` and `BLAS_LIBRARIES` CMake variables.

Fixed minor CMake bugs and included additional error checking during CMake
configuration.

### Bug Fixes

#### ARKODE

Fixed `RCONST` usage in `arkode_butcher.c`.

Fixed bug in `arkInitialSetup` to ensure the mass matrix vector product is
set up before the "msetup" routine is called.

Fixed ARKODE `printf`-related compiler warnings when building SUNDIALS
with extended precision.

#### CVODE and CVODES

In `CVodeFree` now calls `lfree` unconditionally (if non-NULL).

#### IDA and IDAS

Added missing prototype for `IDASetMaxBacksIC` in `ida.h` and `idas.h`.

#### KINSOL

Corrected KINSOL Fortran name translation for `FKIN_SPFGMR`.

Renamed `KINLocalFn` and `KINCommFn` to `KINBBDLocalFn` and `KINBBDCommFn`
respectively in the BBD preconditioner module for consistency with other
SUNDIALS solvers.

## Changes to SUNDIALS in release 2.7.0

### New Features and Enhancements

Two additional `N_Vector` implementations were added -- one for hypre parallel
vectors and one for PETSc vectors. These additions are accompanied by additions
to various interface functions and to user documentation.

Added a new `N_Vector` function, `N_VGetVectorID`, that returns
an identifier for the vector.

The sparse matrix structure was enhanced to support both CSR and CSC matrix
storage formats.

Various additions were made to the KLU and SuperLU_MT sparse linear solver
interfaces, including support for the CSR matrix format when using KLU.

In all packages, the linear solver and preconditioner `free` routines were
updated to return an integer.

In all packages, example codes were updated to use `N_VGetArrayPointer_*`
rather than the `NV_DATA` macro when using the native vectors shipped with
SUNDIALS.

Additional example programs were added throughout including new examples
utilizing the OpenMP vector.

#### ARKODE

The ARKODE implicit predictor algorithms were updated: methods 2 and 3 were
improved slightly, a new predictor approach was added, and the default choice
was modified.

The handling of integer codes for specifying built-in ARKODE Butcher tables was
enhanced. While a global numbering system is still used, methods now have
`#defined` names to simplify the user interface and to streamline
incorporation of new Butcher tables into ARKODE.

The maximum number of Butcher table stages was increased from 8 to 15 to
accommodate very high order methods, and an 8th-order adaptive ERK method was
added.

Support was added for the explicit and implicit methods in an additive
Runge-Kutta method with different stage times to support new SSP-ARK methods.

The FARKODE interface was extended to include a routine to set
scalar/array-valued residual tolerances, to support Fortran applications with
non-identity mass-matrices.

#### IDA and IDAS

The optional input function `IDASetMaxBacksIC` was added to set the
maximum number of linesearch backtracks in the initial condition calculation.

### Bug Fixes

Various minor fixes to installation-related files.

Fixed some examples with respect to the change to use new macro/function names
e.g.,  `SUNRexp`, etc.

In all packages, a memory leak was fixed in the banded preconditioner and
banded-block-diagonal preconditioner interfaces.

Corrected name `N_VCloneEmptyVectorArray` to `N_VCloneVectorArrayEmpty` in
all documentation files.

Various corrections were made to the interfaces to the sparse solvers KLU and
SuperLU_MT.

For each linear solver, the various solver performance counters are now
initialized to 0 in both the solver specification function and in the solver
`linit` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem solution.

#### ARKODE

The missing `ARKSpilsGetNumMtimesEvals` function was added -- this had been
included in the previous documentation but had not been implemented.

The choice of the method vs embedding the Billington and TRBDF2 explicit
Runge-Kutta methods were swapped, since in those the lower-order coefficients
result in an A-stable method, while the higher-order coefficients do not. This
change results in significantly improved robustness when using those methods.

A bug was fixed for the situation where a user supplies a vector of absolute
tolerances, and also uses the vector Resize functionality.

A bug was fixed wherein a user-supplied Butcher table without an embedding is
supplied, and the user is running with either fixed time steps (or they do
adaptivity manually); previously this had resulted in an error since the
embedding order was below 1.

#### CVODE

Corrections were made to three Fortran interface functions.

In FCVODE, fixed argument order bugs in the `FCVKLU` and `FCVSUPERLUMT`
linear solver interfaces.

Added missing Fortran interface routines for supplying a sparse Jacobian routine
with sparse direct solvers.

#### CVODES

A bug was fixed in the interpolation functions used in solving backward problems
for adjoint sensitivity analysis.

In the interpolation routines for backward problems, added logic to bypass
sensitivity interpolation if input sensitivity argument is `NULL`.

Changed each the return type of `*FreeB` functions to `int` and added
`return(0)` to each.

#### IDA

Corrections were made to three Fortran interface functions.

Corrected the output from the `idaFoodWeb_bnd.c` example, the wrong component
was printed in `PrintOutput`.

#### IDAS

In the interpolation routines for backward problems, added logic to bypass
sensitivity interpolation if input sensitivity argument is `NULL`.

Changed each the return type of `*FreeB` functions to `int` and added
`return(0)` to each.

Corrections were made to three Fortran interface functions.

Added missing Fortran interface routines for supplying a sparse Jacobian routine
with sparse direct solvers.

#### KINSOL

The Picard iteration return was changed to always return the newest iterate upon
success.

A minor bug in the line search was fixed to prevent an infinite loop when the
beta condition fails and lambda is below the minimum size.

Corrections were made to three Fortran interface functions.

The functions `FKINCREATE` and `FKININIT` were added to split the
`FKINMALLOC` routine into two pieces. `FKINMALLOC` remains for backward
compatibility, but documentation for it has been removed.

Added missing Fortran interface routines for supplying a sparse Jacobian routine
with sparse direct solvers.

### Matlab Interfaces Removed

Removed the Matlab interface from distribution as it has not been updated since
2009.

## Changes to SUNDIALS in release 2.6.2

### New Features and Enhancements

Various minor fixes to installation-related files

In KINSOL and ARKODE, updated the Anderson acceleration implementation with QR
updating.

In CVODES and IDAS, added `ReInit` and `SetOrdering` wrappers for backward
problems.

In IDAS, fixed for-loop bugs in `IDAAckpntAllocVectors` that could lead to a
memory leak.

### Bug Fixes

Updated the BiCGStab linear solver to remove a redundant dot product call.

Fixed potential memory leak in KLU `ReInit` functions in all solvers.

In ARKODE, fixed a bug in the Cash-Karp Butcher table where the method and
embedding coefficient were swapped.

In ARKODE, fixed error in `arkDoErrorTest` in recovery after failure.

In CVODES, added `CVKLUB` prototype and corrected `CVSuperLUMTB` prototype.

In the CVODES and IDAS header files, corrected documentation of backward
integration functions, especially the `which` argument.

In IDAS, added missing backward problem support functions `IDALapackDenseB`,
`IDALapackDenseFreeB`, `IDALapackBandB`, and `IDALapackBandFreeB`.

In IDAS, made SuperLUMT call for backward problem consistent with CVODES.

In CVODE, IDA, and ARKODE, fixed Fortran interfaces to enable calls to
`GetErrWeights`, `GetEstLocalErrors`, and `GetDky` within a time step.

## Changes to SUNDIALS in release 2.6.1

Fixed loop limit bug in `SlsAddMat` function.

In all six solver interfaces to KLU and SuperLUMT, added `#include` lines, and
removed redundant KLU structure allocations.

Minor bug fixes in ARKODE.

## Changes to SUNDIALS in release 2.6.0

### Autotools Build Option Removed

With this version of SUNDIALS, support and documentation of the Autotools mode
of installation is being dropped, in favor of the CMake mode, which is
considered more widely portable.

### New Package: ARKODE

Addition of ARKODE package of explicit, implicit, and additive Runge-Kutta
methods for ODEs. This package API is close to CVODE so switching between the
two should be straightforward. Thanks go to Daniel Reynolds for the addition
of this package.

### New Features and Enhancements

Added OpenMP and Pthreads `N_Vector` implementations for thread-parallel
computing environments.

Two major additions were made to the linear system solvers available in all
packages. First, in the serial case, an interface to the sparse direct solver
KLU was added. Second, an interface to SuperLU_MT, the multi-threaded version
of SuperLU, was added as a thread-parallel sparse direct solver option, to be
used with the serial version of the `N_Vector` module. As part of these
additions, a sparse matrix (CSC format) structure was added to CVODE.

#### KINSOL

Two major additions were made to the globalization strategy options (`KINSol`
argument `strategy`). One is fixed-point iteration, and the other is Picard
iteration. Both can be accelerated by use of the Anderson acceleration
method. See the relevant paragraphs in the Mathematical Considerations chapter
is the user guide.

An interface to the Flexible GMRES iterative linear solver was added.

### Bug Fixes

In order to avoid possible name conflicts, the mathematical macro and function
names `MIN`, `MAX`, `SQR`, `RAbs`, `RSqrt`, `RExp`, `RPowerI`, and
`RPowerR` were changed to `SUNMIN`, `SUNMAX`, `SUNSQR`, `SUNRabs`,
`SUNRsqrt`, `SUNRexp`, `SRpowerI`, and `SUNRpowerR`, respectively. These
names occur in both the solver and example programs.

In the LAPACK banded linear solver interfaces, the line `smu = MIN(N-1,mu+ml)`
was changed to `smu = mu + ml` to correct an illegal input error for to
`DGBTRF` and `DGBTRS`.

In all Fortran examples, integer declarations were revised so that those which
must match a C type `long int` are declared `INTEGER*8`, and a comment was
added about the type match. All other integer declarations are just
`INTEGER`. Corresponding minor corrections were made to the user guide.

#### CVODE and CVODES

In `cvRootFind`, a minor bug was corrected, where the input array was ignored,
and a line was added to break out of root-search loop if the initial interval
size is below the tolerance `ttol`.

Two minor bugs were fixed regarding the testing of input on the first call to
`CVode` -- one involving `tstop` and one involving the initialization of
`*tret`.

The example program `cvAdvDiff_diag_p` was added to illustrate the use of in
parallel.

In the FCVODE optional input routines `FCVSETIIN` and `FCVSETRIN`, the
optional fourth argument `key_length` was removed, with hardcoded key string
lengths passed to all tests.

In order to eliminate or minimize the differences between the sources for
private functions in CVODE and CVODES, the names of many private functions were
changed from `CV*` to `cv*` and a few other names were also changed.

An option was added in the case of Adjoint Sensitivity Analysis with dense or
banded Jacobian. With a call to `CVDlsSetDenseJacFnBS` or
`CVDlsSetBandJacFnBS`, the user can specify a user-supplied Jacobian function
of type `CVDls***JacFnBS`, for the case where the backward problem depends on
the forward sensitivities.

In `CVodeQuadSensInit`, the line `cv_mem->cv_fQS_data = ...` was corrected
(missing `Q`).

In the CVODES User Guide, a paragraph was added in Section 6.2.1 on
`CVodeAdjReInit`, and a paragraph was added in Section 6.2.9 on
`CVodeGetAdjY`. In the example `cvsRoberts_ASAi_dns`, the output was revised
to include the use of `CVodeGetAdjY`.

For the Adjoint Sensitivity Analysis case in which the backward problem depends
on the forward sensitivities, options have been added to allow for user-supplied
`pset`, `psolve`, and `jtimes` functions.

In the example `cvsHessian_ASA_FSA`, an error was corrected in the function
`fB2`, `y2` in place of `y3` in the third term of `Ith(yBdot,6)`.

#### IDA and IDAS

In `IDARootfind`, a minor bug was corrected, where the input array `rootdir`
was ignored, and a line was added to break out of root-search loop if the
initial interval size is below the tolerance `ttol`.

A minor bug was fixed regarding the testing of the input `tstop` on the first
call to `IDASolve`.

In the FIDA optional input routines `FIDASETIIN`, `FIDASETRIN`, and
`FIDASETVIN`, the optional fourth argument `key_length` was removed, with
hardcoded key string lengths passed to all `strncmp` tests.

An option was added in the case of Adjoint Sensitivity Analysis with dense or
banded Jacobian. With a call to `IDADlsSetDenseJacFnBS` or
`IDADlsSetBandJacFnBS`, the user can specify a user-supplied Jacobian function
of type `IDADls***JacFnBS`, for the case where the backward problem depends on
the forward sensitivities.

#### KINSOL

In function `KINStop`, two return values were corrected to make the values of
`uu` and `fval` consistent.

A bug involving initialization of `mxnewtstep` was fixed. The error affects
the case of repeated user calls to `KINSol` with no intervening call to
`KINSetMaxNewtonStep`.

A bug in the increments for difference quotient Jacobian approximations was
fixed in function `kinDlsBandDQJac`.

In the FKINSOL module, an incorrect return value `ier` in `FKINfunc` was
fixed.

In the FKINSOL optional input routines `FKINSETIIN`, `FKINSETRIN`, and
`FKINSETVIN`, the optional fourth argument `key_length` was removed, with
hardcoded key string lengths passed to all `strncmp` tests.

## Changes to SUNDIALS in release 2.5.0

### Integer Type Change

One significant design change was made with this release, the problem size and
its relatives, bandwidth parameters, related internal indices, pivot arrays, and
the optional output `lsflag` have all been changed from type `int` to type
`long int`, except for the problem size and bandwidths in user calls to
routines specifying BLAS/LAPACK routines for the dense/band linear solvers. The
function `NewIntArray` is replaced by a pair `NewIntArray` /
`NewLintArray`, for `int` and `long int` arrays, respectively.

### Bug Fixes

In the installation files, we modified the treatment of the macro
`SUNDIALS_USE_GENERIC_MATH`, so that the parameter `GENERIC_MATH_LIB` is
either defined (with no value) or not defined.

In all packages, after the solver memory is created, it is set to zero before
being filled.

In each linear solver interface function, the linear solver memory is freed on
an error return, and the function now includes a line setting to `NULL` the
main memory pointer to the linear solver memory.

#### Rootfinding

In CVODE(S) and IDA(S), in the functions `Rcheck1` and `Rcheck2`, when an
exact zero is found, the array `glo` of `g` values at the left endpoint
is adjusted, instead of shifting the `t` location `tlo` slightly.

#### CVODE and CVODES

In `CVSetTqBDF`, the logic was changed to avoid a divide by zero.

In a minor change to the CVODES user interface, the type of the index `which`
was changed from `long int` to `int`.

Errors in the logic for the integration of backward problems in CVODES were
identified and fixed.

#### IDA and IDAS

To be consistent with IDAS, IDA uses the function `IDAGetDky` for optional
output retrieval.

A memory leak was fixed in two of the `IDASp***Free` functions.

A missing vector pointer setting was added in `IDASensLineSrch`.

In `IDACompleteStep`, conditionals around lines loading a new column of three
auxiliary divided difference arrays, for a possible order increase, were fixed.

#### KINSOL

Three major logic bugs were fixed - involving updating the solution vector,
updating the linesearch parameter, and a missing error return.

Three minor errors were fixed - involving setting `etachoice` in the
Matlab/KINSOL interface, a missing error case in `KINPrintInfo`, and avoiding
an exponential overflow in the evaluation of `omega`.

## Changes to SUNDIALS in release 2.4.0

Added a CMake-based build option in addition to the one based on autotools.

The user interface has been further refined. Some of the API changes involve:

(a) a reorganization of all linear solver modules into two families (besides the
    existing family of scaled preconditioned iterative linear solvers, the
    direct solvers, including new LAPACK-based ones, were also organized into a
    *direct* family);

(b) maintaining a single pointer to user data, optionally specified through a
    `Set`-type function; and

(c) a general streamlining of the preconditioner modules distributed with the
    solvers.

Added interfaces to LAPACK linear solvers for dense and banded matrices to all
packages.

An option was added to specify which direction of zero-crossing is to be
monitored while performing rootfinding in CVODE(S) and IDA(S).

CVODES includes several new features related to sensitivity analysis, among
which are:

(a) support for integration of quadrature equations depending on both the states
    and forward sensitivity (and thus support for forward sensitivity analysis
    of quadrature equations),

(b) support for simultaneous integration of multiple backward problems based on
    the same underlying ODE (e.g., for use in an *forward-over-adjoint* method
    for computing second order derivative information),

(c) support for backward integration of ODEs and quadratures depending on both
    forward states and sensitivities (e.g., for use in computing second-order
    derivative information), and

(d) support for reinitialization of the adjoint module.

Moreover, the prototypes of all functions related to integration of backward
problems were modified to support the simultaneous integration of multiple
problems.

All backward problems defined by the user are internally managed through a
linked list and identified in the user interface through a unique identifier.

## Changes to SUNDIALS in release 2.3.0

### New Features and Enhancements

The main changes in this release involve a rearrangement of the entire
SUNDIALS source tree. At the user interface level, the main impact is in the
mechanism of including SUNDIALS header files which must now include the relative
path e.g., `#include <cvode/cvode.h>` as all exported header files are now
installed in separate subdirectories of the installation *include* directory.

The functions in the generic dense linear solver (`sundials_dense` and
`sundials_smalldense`) were modified to work for rectangular `m x n` (`m <= n`),
while the factorization and solution functions were renamed to `DenseGETRF` /
`denGETRF` and `DenseGETRS` / `denGETRS`, respectively. The factorization and
solution functions in the generic band linear solver were renamed `BandGBTRF`
and `BandGBTRS`, respectively.

In IDA, the user interface to the consistent initial conditions calculations was
modified. The `IDACalcIC` arguments `t0`, `yy0`, and `yp0` were
removed and a new function, `IDAGetConsistentIC` is provided.

### Bug Fixes

In the CVODES adjoint solver module, the following two bugs were fixed:

* In `CVodeF` the solver was sometimes incorrectly taking an additional step
  before returning control to the user (in `CV_NORMAL` mode) thus leading to
  a failure in the interpolated output function.

* In `CVodeB`, while searching for the current check point, the solver was
  sometimes reaching outside the integration interval resulting in a
  segmentation fault.

In IDA, a bug was fixed in the internal difference-quotient dense and banded
Jacobian approximations, related to the estimation of the perturbation (which
could have led to a failure of the linear solver when zero components with
sufficiently small absolute tolerances were present).

## Changes to SUNDIALS in release 2.2.0

### New Header Files Names

To reduce the possibility of conflicts, the names of all header files have been
changed by adding unique prefixes (e.g., `cvode_` and `sundials_`). When
using the default installation procedure, the header files are exported under
various subdirectories of the target `include` directory. For more details see
Appendix the installation chapter in the user guide.

### Build System Changes

Updated configure script and Makefiles for Fortran examples to avoid C++
compiler errors (now use `CC` and `MPICC` to link only if necessary).

The shared object files are now linked into each SUNDIALS library rater than
into a separate `libsundials_shared` library.

### New Features and Enhancements

Deallocation functions now take the address of the respective memory block
pointer as the input argument.

Interfaces to the Scaled Preconditioned Bi-CGstab (SPBCG) and Scaled
Preconditioned Transpose-Free Quasi-Minimal Residual (SPTFQMR) linear solver
modules have been added to all packages. At the same time, function type names
for Scaled Preconditioned Iterative Linear Solvers were added for the
user-supplied Jacobian-times-vector and preconditioner setup and solve
functions. Additionally, in KINSOL interfaces have been added to the SUNDIALS
DENSE, and BAND linear solvers and include support for nonlinear residual
monitoring which can be used to control Jacobian updating.

A new interpolation method was added to the CVODES adjoint module. The function
`CVadjMalloc` has an additional argument which can be used to select the
desired interpolation scheme.

FIDA, a Fortran-C interface module, was added.

The rootfinding feature was added to IDA, whereby the roots of a set of given
functions may be computed during the integration of the DAE system.

In IDA a user-callable routine was added to access the estimated local error
vector.

In the KINSOL Fortran interface module, FKINSOL, optional inputs are now set
using `FKINSETIIN` (integer inputs), `FKINSETRIN` (real inputs), and
`FKINSETVIN` (vector inputs). Optional outputs are still obtained from the
`IOUT` and `ROUT` arrays which are owned by the user and passed as arguments
to `FKINMALLOC`.

## Changes to SUNDIALS in release 2.1.1

The function `N_VCloneEmpty` was added to the global vector operations table.

A minor bug was fixed in the interpolation functions of the adjoint CVODES
module.

## Changes to SUNDIALS in release 2.1.0

The user interface has been further refined. Several functions used for setting
optional inputs were combined into a single one.

In CVODE(S) and IDA, an optional user-supplied routine for setting the error
weight vector was added.

Additionally, to resolve potential variable scope issues, all SUNDIALS solvers
release user data right after its use.

The build systems has been further improved to make it more robust.

## Changes to SUNDIALS in release 2.0.2

Fixed autoconf-related bug to allow configuration with the PGI Fortran compiler.

Modified the build system to use customized detection of the Fortran name
mangling scheme (autoconf's `AC_F77_WRAPPERS` routine is problematic on some
platforms).

A bug was fixed in the `CVode` function that was potentially leading to
erroneous behavior of the rootfinding procedure on the integration first step.

A new chapter in the User Guide was added - with constants that appear in the
user interface.

## Changes to SUNDIALS in release 2.0.1

### Build System

Changed the order of compiler directives in header files to avoid compilation
errors when using a C++ compiler.

Changed the method of generating `sundials_config.h` to avoid potential
warnings of redefinition of preprocessor symbols.

### New Features

In CVODES the option of activating and deactivating forward sensitivity
calculations on successive runs without memory allocation and deallocation.

### Bug Fixes

In CVODES bug fixes related to forward sensitivity computations (possible loss
of accuracy on a BDF order increase and incorrect logic in testing user-supplied
absolute tolerances) were made.

## Changes to SUNDIALS in release 2.0.0

Installation of all of SUNDIALS packages has been completely redesigned and is
now based on configure scripts.

The major changes from the previous version involve a redesign of the user
interface across the entire SUNDIALS suite. We have eliminated the mechanism of
providing optional inputs and extracting optional statistics from the solver
through the `iopt` and `ropt` arrays. Instead, packages now provide `Set`
functions to change the default values for various quantities controlling the
solver and `Get` functions to extract statistics after return from the main
solver routine.

Additionally, the interfaces to several user-supplied routines (such as those
providing Jacobians and preconditioner information) were simplified by reducing
the number of arguments. The same information that was previously accessible
through such arguments can now be obtained through `Get`-type functions.

In CVODE and CVODES a rootfinding feature was added, whereby the roots of a set
of given functions may be computed during the integration of the ODE system.

Changes to the NVector:

* Removed `machEnv`, redefined table of vector operations (now contained in
  the `N_Vector` structure itself).

* All SUNDIALS functions create new `N_Vector` variables through
  cloning, using an `N_Vector` passed by the user as a template.

* A particular vector implementation is supposed to provide user-callable
  constructor and destructor functions.

* Removed the following functions from the structure of vector operations:
  `N_VNew`, `N_VNew_S`, `N_VFree`, `N_VFree_S`, `N_VMake`,
  `N_VDispose`, `N_VGetData`, `N_VSetData`, `N_VConstrProdPos`, and
  `N_VOneMask`.

* Added the following functions to the structure of vector operations:
  `N_VClone`, `N_VDestroy`, `N_VSpace`, `N_VGetArrayPointer`,
  `N_VSetArrayPointer`, and `N_VWrmsNormMask`.

* Note that `nvec_ser` and `nvec_par` are now separate modules outside the
  shared SUNDIALS module.

Changes to the linear solvers:

* In SPGMR, added a dummy `N_Vector` argument to be used as a template for
  cloning.

* In SPGMR, removed `N` (problem dimension) from the argument list of
  `SpgmrMalloc`.

* Replaced `iterativ.{c,h}` with `iterative.{c,h}`.

* Modified constant names in `iterative.h` (preconditioner types are prefixed
  with `PREC_`).

* Changed numerical values for `MODIFIED_GS` (from `0` to `1`) and
  `CLASSICAL_GS` (from `1` to `2`).

Changes to `sundialsmath` submodule:

* Replaced the internal routine for estimating unit roundoff with definition of
  unit roundoff from `float.h`.

* Modified functions to call the appropriate math routines given the precision
  level specified by the user.

Changes to `sundialstypes` submodule:

* Removed `integertype`.

* Added definitions for `BIG_REAL`, `SMALL_REAL`, and `UNIT_ROUNDOFF`
  using values from `float.h` based on the precision.

* Changed definition of macro `RCONST` to depend on the precision level
  specified by the user.
