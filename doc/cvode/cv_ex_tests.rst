Parallel tests
==============

The stiff example problem ``cvDiurnal_kry`` described above, or rather its parallel
version ``cvDiurnal_kry_p``, has been modified and expanded to form a test problem for
the parallel version of CVODE. This work was largely carried out by M. Wittman and
reported in [Wit:96]_.

To start with, in order to add realistic complexity to the solution, the initial profile
for this problem was altered to include a rather steep front in the vertical direction.
Specifically, the function :math:`\beta(y)` in Eq. (``cvDiurnalic``) has been replaced by:

.. math::

   \beta(y) = 0.75 + 0.25 \tanh(10\,y - 400)

This function rises from about 0.5 to about 1.0 over a *y* interval of about 0.2
(i.e., 1/100 of the total span in *y*). This vertical variation, together with the
horizontal advection and diffusion in the problem, demands a fairly fine spatial mesh
to achieve acceptable resolution.

In addition, an alternate choice of differencing is used in order to control spurious
oscillations resulting from the horizontal advection. In place of central differencing
for that term, a biased upwind approximation is applied to each of the terms :math:`\partial c^i/\partial x`, namely:

.. math::

   \left. \frac{\partial c}{\partial x} \right|_{x_j} \approx
   \frac{\frac{3}{2}\,c_{j+1} - c_j - \frac{1}{2}\,c_{j-1}}{2\,\Delta x}

With this modified form of the problem, tests similar to those described above for the example
were performed. Here we fix the subgrid dimensions at ``MXSUB = MYSUB = 50``, so that the
local (per-processor) problem size is 5000, while the processor array dimensions, ``NPEX``
and ``NPEY``, are varied.

In one (typical) sequence of tests, we fix ``NPEY = 8`` (for a vertical mesh size of ``MY = 400``)
and set:

- ``NPEX = 8`` (``MX = 400``)
- ``NPEX = 16`` (``MX = 800``)
- ``NPEX = 32`` (``MX = 1600``)

Thus the largest problem size, *N*, is computed as:

   2 \times 400 \times 1600 = 1,280,000.

For these tests, the maximum Krylov dimension (``maxl``) is raised to 10 (from its default of 5).

For each of the three test cases, the test program was run on a Cray-T3D (256 processors)
using three different message-passing libraries:

- **MPICH** - an implementation of MPI on top of the Chameleon library
- **EPCC** - an implementation of MPI by the Edinburgh Parallel Computer Centre
- **SHMEM** - Cray's Shared Memory Library

The following table gives the run time and selected performance counters for these 9 runs.
In all cases, the solutions agreed well with each other, showing the expected small variations
with grid size. In the table, "M-P" denotes the message-passing library, RT is the reported run time
(in CPU seconds), ``nst`` is the number of time steps, ``nfe`` the number of *f* evaluations,
``nni`` the number of nonlinear (Newton) iterations, ``nli`` the number of linear (Krylov) iterations,
and ``npe`` the number of evaluations of the preconditioner.

+-------+--------+-------+------+-------+------+-------+------+
| NPEX  | M-P    | RT    | nst  | nfe   | nni  | nli   | npe  |
+=======+========+=======+======+=======+======+=======+======+
| 8     | MPICH  | 436.  | 1391 | 9907  | 1512 | 8392  | 24   |
+-------+--------+-------+------+-------+------+-------+------+
| 8     | EPCC   | 355.  | 1391 | 9907  | 1512 | 8392  | 24   |
+-------+--------+-------+------+-------+------+-------+------+
| 8     | SHMEM  | 349.  | 1999 | 10326 | 2096 | 8227  | 34   |
+-------+--------+-------+------+-------+------+-------+------+
| 16    | MPICH  | 676.  | 2513 | 14159 | 2583 | 11573 | 42   |
+-------+--------+-------+------+-------+------+-------+------+
| 16    | EPCC   | 494.  | 2513 | 14159 | 2583 | 11573 | 42   |
+-------+--------+-------+------+-------+------+-------+------+
| 16    | SHMEM  | 471.  | 2513 | 14160 | 2581 | 11576 | 42   |
+-------+--------+-------+------+-------+------+-------+------+
| 32    | MPICH  | 1367. | 2536 | 20153 | 2696 | 17454 | 43   |
+-------+--------+-------+------+-------+------+-------+------+
| 32    | EPCC   | 737.  | 2536 | 20153 | 2696 | 17454 | 43   |
+-------+--------+-------+------+-------+------+-------+------+
| 32    | SHMEM  | 695.  | 2536 | 20121 | 2694 | 17424 | 43   |
+-------+--------+-------+------+-------+------+-------+------+

Some of the results were as expected, and some were surprising. For a given mesh size, variations
in performance counts were small or absent-except for moderate (but still acceptable) variations for
SHMEM in the smallest case. The increase in costs with mesh size can be attributed to a decline in
the quality of the preconditioner, which neglects most of the spatial coupling. The preconditioner
quality can be inferred from the ratio ``nli/nni``, which is the average number of Krylov iterations
per Newton iteration. The most interesting (and unexpected) result is the variation of run time
with library: SHMEM is the most efficient, EPCC is a very close second, and MPICH loses considerable
efficiency by comparison as the problem size grows. This means that the highly portable MPI version
of ``cvode``, with an appropriate choice of MPI implementation, is fully competitive with the
Cray-specific version using the SHMEM library. While the overall costs do not represent a well-scaled
parallel algorithm (because of the preconditioner choice), the cost per function evaluation is quite
flat for EPCC and SHMEM (0.033 to 0.037), whereas for MPICH it ranges from 0.044 to 0.068.

For tests that demonstrate speedup from parallelism, runs with fixed problem size are considered:
``MX = 800`` and ``MY = 400``. Here the vertical subgrid dimension is fixed at ``MYSUB = 50`` and the
vertical processor array dimension at ``NPEY = 8``, while the corresponding horizontal sizes are varied.
Specifically, we take:

- ``NPEX = 8`` with ``MXSUB = 100``
- ``NPEX = 16`` with ``MXSUB = 50``
- ``NPEX = 32`` with ``MXSUB = 25``

The runs for the three cases and the three message-passing libraries all show very good agreement in
solution values and performance counts. The run times for EPCC are 947, 494, and 278 seconds,
showing speedups of 1.92 and 1.78 as the number of processors is doubled (twice). For the SHMEM runs,
the times were slightly lower with speedup ratios of 1.98 and 1.91. For MPICH, consistent with the earlier
runs, the run times were considerably higher and showed speedup ratios of only 1.54 and 1.03.
