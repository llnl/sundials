.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Added support to compute internal difference-quotient approximations to sparse
Jacobians (CSC format only) for all SUNDIALS packages.

**Bug Fixes**

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

**Deprecation Notices**
