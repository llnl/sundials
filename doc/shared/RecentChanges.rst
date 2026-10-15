.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Reorganized the CVODE examples by problem, with implementation language and
optional third-party libraries represented by variant subdirectories. Each
variant contains its input files, scripts, and reference output. The old CVODE
example Makefiles are no longer installed; use the installed CMakeLists.txt
files instead.

Reorganized the ARKODE examples using the same problem-based structure as the
CVODE examples, with implementation language and optional third-party
libraries represented by variant subdirectories. Each variant contains its
input files, scripts, and reference output. The old ARKODE example Makefiles
are no longer installed; use the installed CMakeLists.txt files instead.

Reorganized the CVODES examples using the same problem-based structure as the
CVODE and ARKODE examples. The old CVODES example Makefiles are no longer
installed; use the installed CMakeLists.txt files instead.

The header-only Kokkos and Ginkgo N_Vector, SUNMatrix, and SUNLinearSolver
targets are now exported by installed SUNDIALS CMake packages. The installed
Ginkgo targets also propagate their Ginkgo dependency to consumers.

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Updated :c:func:`ARKodeGetReturnFlagName`, :c:func:`CVodeGetReturnFlagName`,
:c:func:`CVDiagGetReturnFlagName`, :c:func:`IDAGetReturnFlagName`, and
:c:func:`KINGetReturnFlagName` to return ``const char*`` rather than ``char*``.
The returned string should no longer be freed by the user.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

**Deprecation Notices**
