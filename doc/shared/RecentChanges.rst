.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the function :c:func:`ARKodeSkipAdaptStopTime` to specify that
stop-time-limited steps should be disregarded when selecting step sizes
for time step adaptivity.

An optional N_Vector routine, :c:func:`N_VCopy`, was added, to streamline data copies between two
vectors.  For user-supplied N_Vector modules that do not provide this function, :c:func:`N_VScale`
will be used instead.

**Bug Fixes**

**Deprecation Notices**
