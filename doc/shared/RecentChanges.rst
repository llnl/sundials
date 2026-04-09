.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Updated the MRIStep module in ARKODE to support *extended
Super Time Stepping* (ExtSTS) methods for multi-physics simulations involving
parabolic components, through the :c:func:`MRIStepCreateExtSTS` function.
See the documentation section on
:ref:`the mathematical methods in MRIStep <ARKODE.Mathematics.MRIStep.ExtSTS>`
for more details on structure of ExtSTS methods, and the documentation
mathematical formulation of these methods, and the documentation section on
:ref:`a skeleton of usage for MRIStep <ARKODE.Usage.MRIStep.Skeleton-ExtSTS>`
for details on its usage.

**Bug Fixes**

Fixed a minor bug where the number of required stages for STS methods
in the LSRKStep module was incorrectly computed using the spectral
radius instead of the real part of the Jacobian eigenvalues.

**Deprecation Notices**
