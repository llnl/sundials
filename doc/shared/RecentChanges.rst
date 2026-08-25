.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

Added support for complex-valued data throughout SUNDIALS, including a new type,
:c:type:`suncomplextype`, that is appropriately defined based on the requested
SUNDIALS floating-point precision.  Created a new type
:c:type:`sunscalartype` that can be configured during installation as an alias to
either :c:type:`sunrealtype` or :c:type:`suncomplextype`, and that will be used
in all vector data declarations.  Added new mathematical library functions to
``sundials_math.h`` for :c:type:`suncomplextype` numbers of the configured
precision (e.g., ``SUN_CREAL``, ``SUN_CIMAG``, ``SUN_CCONJ``, ``SUNCsqrt``, and
``SUNCabs``).  Added generic :c:type:`sunscalartype` mathematical functions
(e.g., ``SUN_REAL``, ``SUN_IMAG``, ``SUN_CONJ``, ``SUNsqrt``, and ``SUNabs``) that
call either the relevant :c:type:`sunrealtype` or :c:type:`suncomplextype` function
to match the configuration of the :c:type:`sunscalartype` alias.

**New Features and Enhancements**

Added the :cmakeop:`SUNDIALS_ENABLE_DEPRECATED_WARNINGS` CMake option to allow
users to disable compiler warnings emitted by deprecated SUNDIALS API
annotations.

Added the function :c:func:`N_VSetDeviceArrayPointer` to the N_Vector API to set
the device data pointer for vector implementations with a device memory space.

**Bug Fixes**

Fixed bug in SUNNonlinearSolver_Auto which resulted in premature switch to Newton from fixed point
due to convergence rate check occurring after only one iteration.

**Deprecation Notices**
