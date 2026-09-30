.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Added optional Thrust based reduction implementations to the :ref:`HIP NVector
<NVectors.HIP>`. These are enabled by setting the new CMake option
:cmakeop:`SUNDIALS_ENABLE_THRUST_REDUCTIONS` to ``ON``, which requires building
with C++17 or later, and are used by selecting the new
:cpp:func:`SUNHipThrustExecPolicy` reduction execution policy.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

**Deprecation Notices**
