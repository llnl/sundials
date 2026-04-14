.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

Fixed a minor bug where the number of required stages for STS methods 
in the LSRKStep module was incorrectly computed using the spectral 
radius instead of the real part of the Jacobian eigenvalues.

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

**Deprecation Notices**
