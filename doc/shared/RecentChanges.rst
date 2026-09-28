.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

**Deprecation Notices**
