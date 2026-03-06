.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

Fixed a CMake bug where the SuperLU_MT interface would not be built and
installed without setting the ``SUPERLUMT_WORKS`` option to ``TRUE``.

Fixed the embedded coefficients for the ``ARKODE_TSITOURAS_7_4_5`` Butcher
table.

Fixed a bug where passing an empty string to ``SUNLogger_Set{Error,Warning,Info,Debug}Filename``
did not disable the corresponding logging stream `Issue #844 <https://github.com/llnl/sundials/issues/844>`__.

**Deprecation Notices**
