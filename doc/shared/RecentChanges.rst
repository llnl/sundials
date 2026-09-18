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

**Bug Fixes**

**Deprecation Notices**
