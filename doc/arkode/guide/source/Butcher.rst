.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ UMBC
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Butcher:

==============
Butcher Tables
==============

Here we catalog the full set of Butcher tables included in ARKODE. We group
these into four categories: *explicit*, *implicit*, *additive* and
*symplectic partitioned*.
However, since the methods that comprise an additive Runge--Kutta method are
themselves explicit and implicit, their component Butcher tables are listed
within their separate sections, but are referenced together in the additive
section.

In each of the following tables, we use the following notation (shown
for a 3-stage method):

.. math::

   \begin{array}{r|ccc}
     c_1 & a_{1,1} & a_{1,2} & a_{1,3} \\
     c_2 & a_{2,1} & a_{2,2} & a_{2,3} \\
     c_3 & a_{3,1} & a_{3,2} & a_{3,3} \\
     \hline
     q & b_1 & b_2 & b_3 \\
     p & \tilde{b}_1 & \tilde{b}_2 & \tilde{b}_3
   \end{array}

where here the method and embedding share stage :math:`A` and
:math:`c` values, but use their stages :math:`z_i` differently through
the coefficients :math:`b` and :math:`\tilde{b}` to generate methods
of orders :math:`q` (the main method) and :math:`p` (the embedding,
typically :math:`q = p+1`, though sometimes this is reversed).

Method authors often use different naming conventions to categorize
their methods.  For each of the methods below with an embedding, we follow the
uniform naming convention:

.. code-block:: text

   NAME-S-P-Q

where here

* ``NAME`` is the author or the name provided by the author (if applicable),
* ``S`` is the number of stages in the method,
* ``P`` is the global order of accuracy for the embedding,
* ``Q`` is the global order of accuracy for the method.

For methods without an embedding (e.g., fixed-step methods) ``P`` is omitted so
that methods follow the naming convention ``NAME-S-Q``.

For SPRK methods, the naming convention is ``SPRK-NAME-S-Q``.

In the code, unique integer IDs are defined inside ``arkode_butcher_erk.h`` and
``arkode_butcher_dirk.h`` for each method, which may be used by calling routines
to specify the desired method. SPRK methods are defined inside ``arkode_sprk.h``.
These names are specified in ``fixed width font`` at the start of each method's
section below.

Additionally, for each method we provide a plot of the linear
stability region in the complex plane.  These have been computed via
the following approach.  For any Runge--Kutta method as defined above,
we may define the stability function

.. math::

   R(\eta) = 1 + \eta b [I - \eta A]^{-1} e,

where :math:`e\in\mathbb{R}^s` is a column vector of all ones, :math:`\eta =
h\lambda` and :math:`h` is the time step size.  If the stability
function satisfies :math:`|R(\eta)| \le 1` for all eigenvalues,
:math:`\lambda`, of :math:`\frac{\partial }{\partial y}f(t,y)` for a
given IVP, then the method will be linearly stable for that problem
and step size.  The stability region

.. math::

   S = \{ \eta\in\mathbb{C}\; :\; \left| R(\eta) \right| \le 1\}

is typically given by an enclosed region of the complex plane, so it
is standard to search for the border of that region in order to
understand the method.  Since all complex numbers with unit magnitude
may be written as :math:`e^{i\theta}` for some value of :math:`\theta`,
we perform the following algorithm to trace out this boundary.

1. Define an array of values ``Theta``.  Since we wish for a
   smooth curve, and since we wish to trace out the entire boundary,
   we choose 10,000 linearly-spaced points from 0 to :math:`16\pi`.
   Since some angles will correspond to multiple locations on the
   stability boundary, by going beyond :math:`2\pi` we ensure that all
   boundary locations are plotted, and by using such a fine
   discretization the Newton method (next step) is more likely to
   converge to the root closest to the previous boundary point,
   ensuring a smooth plot.

2. For each value :math:`\theta \in` ``Theta``, we solve the nonlinear
   equation

   .. math::

      0 = f(\eta) = R(\eta) - e^{i\theta}

   using a finite-difference Newton iteration, using tolerance
   :math:`10^{-7}`, and differencing parameter
   :math:`\sqrt{\varepsilon}` (:math:`\approx 10^{-8}`).

   In this iteration, we use as initial guess the solution from the
   previous value of :math:`\theta`, starting with an initial-initial
   guess of :math:`\eta=0` for :math:`\theta=0`.

3. We then plot the resulting :math:`\eta` values that trace the
   stability region boundary.

We note that for any stable IVP method, the value :math:`\eta_0 =
-\varepsilon + 0i` is always within the stability region.  So in each
of the following pictures, the interior of the stability region is the
connected region that includes :math:`\eta_0`.  Resultingly, methods
whose linear stability boundary is located entirely in the right
half-plane indicate an `A-stable` method.



.. _Butcher.explicit:

Explicit Butcher tables
---------------------------

In the category of explicit Runge--Kutta methods, ARKODE includes methods that
have orders 2 through 9, with embeddings that are of orders 1 through 8.
ARKODE's explicit Butcher tables are provided in the enumeration

.. c:enum:: ARKODE_ERKTableID

with values specified in :numref:`ARKODE.Butcher.ERK_properties`.

.. _ARKODE.Butcher.ERK_properties:
.. table:: Explicit Butcher tables. The default method for each order is marked
   with an asterisk (*).

   +------------------------------------------------------+--------+----------------+-------+
   | Method ID                                            | Stages | Embedded Order | Order |
   +======================================================+========+================+=======+
   | :c:enumerator:`ARKODE_FORWARD_EULER_1_1`             | 1      | ---            | 1*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_RALSTON_3_1_2`                 | 3      | 1              | 2*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_HEUN_EULER_2_1_2`              | 2      | 1              | 2     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_RALSTON_EULER_2_1_2`           | 2      | 1              | 2     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2` | 2      | 1              | 2     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK2_ERK_3_1_2`                | 3      | 1              | 2     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_BOGACKI_SHAMPINE_4_2_3`        | 4      | 2              | 3*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK324L2SA_ERK_4_2_3`          | 4      | 2              | 3     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_SHU_OSHER_3_2_3`               | 3      | 2              | 3     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_KNOTH_WOLKE_3_3`               | 3      | ---            | 3     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_SOFRONIOU_SPALETTA_5_3_4`      | 5      | 3              | 4*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ZONNEVELD_5_3_4`               | 5      | 3              | 4     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK436L2SA_ERK_6_3_4`          | 6      | 3              | 4     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK437L2SA_ERK_7_3_4`          | 7      | 3              | 4     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_SAYFY_ABURUB_6_3_4`            | 6      | 3              | 4     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_TSITOURAS_7_4_5`               | 7      | 4              | 5*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_CASH_KARP_6_4_5`               | 6      | 4              | 5     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_FEHLBERG_6_4_5`                | 6      | 4              | 5     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_DORMAND_PRINCE_7_4_5`          | 7      | 4              | 5     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK548L2SA_ERK_8_4_5`          | 8      | 4              | 5     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK548L2SAb_ERK_8_4_5`         | 8      | 4              | 5     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_VERNER_9_5_6`                  | 9      | 5              | 6*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_VERNER_8_5_6`                  | 8      | 5              | 6     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_VERNER_10_6_7`                 | 10     | 6              | 7*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_VERNER_13_7_8`                 | 13     | 7              | 8*    |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_FEHLBERG_13_7_8`               | 13     | 7              | 8     |
   +------------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_VERNER_16_8_9`                 | 16     | 8              | 9*    |
   +------------------------------------------------------+--------+----------------+-------+


.. c:enumerator:: ARKODE_FORWARD_EULER_1_1

Accessible via the constant ``ARKODE_FORWARD_EULER_1_1`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_FORWARD_EULER_1_1"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 1st order explicit method (from :cite:p:`Euler:68`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|c}
     0 & 0 \\
     \hline
     1 & 1
   \end{array}

.. figure:: /figs/arkode/forward_euler_erk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the forward Euler method.


.. c:enumerator:: ARKODE_RALSTON_3_1_2

Accessible via the constant ``ARKODE_RALSTON_3_1_2`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_RALSTON_3_1_2"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 2nd order explicit method
(primary method from :cite:p:`Ralston:62`).

.. versionchanged:: 6.3.0

   Added as the default 2nd order explicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccc}
     0 & 0 & 0 & 0 \\
     \frac{2}{3} & \frac{2}{3} & 0 & 0 \\
     1 & \frac{1}{4} & \frac{3}{4} & 0 \\
     \hline
     2 & \frac{1}{4} & \frac{3}{4} & 0 \\
     1 & \frac{5}{37} & \frac{2}{3} & \frac{22}{111}
   \end{array}

.. figure:: /figs/arkode/ralston_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the Ralston method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_HEUN_EULER_2_1_2

Accessible via the constant ``ARKODE_HEUN_EULER_2_1_2`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_HEUN_EULER_2_1_2"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(primary method from :cite:p:`Runge:95`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_RALSTON_3_1_2`` as the default 2nd order explicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cc}
     0 & 0 & 0 \\
     1 & 1 & 0 \\
     \hline
     2 & \frac{1}{2} & \frac{1}{2} \\
     1 & 1 & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_0.png
   :scale: 50 %
   :align: center

   Linear stability region for the Heun-Euler method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_RALSTON_EULER_2_1_2

Accessible via the constant ``ARKODE_RALSTON_EULER_2_1_2`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_RALSTON_EULER_2_1_2"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(primary method from :cite:p:`Ralston:62`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cc}
     0 & 0 & 0 \\
     \frac{2}{3} & \frac{2}{3} & 0 \\
     \hline
     2 & \frac{1}{4} & \frac{3}{4} \\
     1 & 1 & 0
   \end{array}

.. figure:: /figs/arkode/ralston_euler_erk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the Ralston-Euler method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2

Accessible via the constant ``ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(primary method from :cite:p:`Runge:95`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cc}
     0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{2} & 0 \\
     \hline
     2 & 0 & 1 \\
     1 & 1 & 0
   \end{array}

.. figure:: /figs/arkode/explicit_midpoint_euler_erk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the Explicit-Midpoint-Euler method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK2_ERK_3_1_2

Accessible via the constant ``ARKODE_ARK2_ERK_3_1_2`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_ARK2_ERK_3_1_2"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the explicit portion of the default 2nd order additive method (the
explicit portion of the ARK2 method from :cite:p:`giraldo2013implicit`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccc}
     0            & 0                           & 0                       & 0 \\
     2 - \sqrt{2} & 2 - \sqrt{2}                & 0                       & 0 \\
     1            & 1 - \frac{3 + 2\sqrt{2}}{6} & \frac{3 + 2\sqrt{2}}{6} & 0 \\
     \hline
     2 & \frac{1}{2\sqrt{2}}    & \frac{1}{2\sqrt{2}}    & 1 - \frac{1}{\sqrt{2}} \\
     1 & \frac{4 - \sqrt{2}}{8} & \frac{4 - \sqrt{2}}{8} & \frac{1}{2\sqrt{2}}    \\
   \end{array}

.. figure:: /figs/arkode/ark2_erk_stab_region.png
   :scale: 65 %
   :align: center

   Linear stability region for the ARK2-ERK method. The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_BOGACKI_SHAMPINE_4_2_3

Accessible via the constant ``ARKODE_BOGACKI_SHAMPINE_4_2_3`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_BOGACKI_SHAMPINE_4_2_3"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 3rd order
explicit method (from :cite:p:`Bogacki:89`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccc}
     0 &   0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 \\
     \frac{3}{4} & 0 & \frac{3}{4} & 0 & 0 \\
     1   & \frac{2}{9} & \frac{1}{3} & \frac{4}{9} & 0 \\
     \hline
     3 & \frac{2}{9} & \frac{1}{3} & \frac{4}{9} \\
     2 & \frac{7}{24} & \frac{1}{4} & \frac{1}{3} & \frac{1}{8}
   \end{array}

.. figure:: /figs/arkode/stab_region_1.png
   :scale: 50 %
   :align: center

   Linear stability region for the Bogacki-Shampine method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK324L2SA_ERK_4_2_3

Accessible via the constant ``ARKODE_ARK324L2SA_ERK_4_2_3`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_ARK324L2SA_ERK_4_2_3"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the explicit portion of
the default 3rd order additive method (the explicit portion of the ARK3(2)4L[2]SA
method from :cite:p:`KenCarp:03`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     \frac{1767732205903}{2027836641118} & \frac{1767732205903}{2027836641118} & 0 & 0 & 0 \\
     \frac{3}{5} & \frac{5535828885825}{10492691773637} & \frac{788022342437}{10882634858940} & 0 & 0 \\
     1 & \frac{6485989280629}{16251701735622} & -\frac{4246266847089}{9704473918619} & \frac{10755448449292}{10357097424841} & 0 \\
     \hline
     3 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     2 & \frac{2756255671327}{12835298489170} & -\frac{10771552573575}{22201958757719} & \frac{9247589265047}{10645013368117} & \frac{2193209047091}{5459859503100}
   \end{array}

.. figure:: /figs/arkode/stab_region_2.png
   :scale: 50 %
   :align: center

   Linear stability region for the explicit ARK-4-2-3 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_SHU_OSHER_3_2_3

Accessible via the constant ``ARKODE_SHU_OSHER_3_2_3`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_SHU_OSHER_3_2_3"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`ShOs:88` with embedding from :cite:p:`FCS:21`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccc}
     0 & 0 & 0 & 0 \\
     1 & 1 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{4} & \frac{1}{4} & 0 \\
     \hline
     3 & \frac{1}{6} & \frac{1}{6} & \frac{2}{3} \\
     2 & \frac{291485418878409}{1000000000000000} & \frac{291485418878409}{1000000000000000} & \frac{208514581121591}{500000000000000}
   \end{array}

.. figure:: /figs/arkode/shu_osher_erk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the Shu-Osher method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_KNOTH_WOLKE_3_3

Accessible via the constant ``ARKODE_KNOTH_WOLKE_3_3`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`, or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_KNOTH_WOLKE_3_3"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName`, or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 3th order slow and fast MRIStep method (from
:cite:p:`KnWo:98`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccc}
               0 & 0             & 0             & 0 \\
     \frac{1}{3} & \frac{1}{3}   & 0             & 0 \\
     \frac{3}{4} & -\frac{3}{16} & \frac{15}{16} & 0 \\
     \hline
               3 & \frac{1}{6} & \frac{3}{10} & \frac{8}{15}
   \end{array}

.. figure:: /figs/arkode/stab_region_24.png
   :scale: 50 %
   :align: center

   Linear stability region for the Knoth-Wolke method


.. c:enumerator:: ARKODE_SOFRONIOU_SPALETTA_5_3_4

Accessible via the constant ``ARKODE_SOFRONIOU_SPALETTA_5_3_4`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_SOFRONIOU_SPALETTA_5_3_4"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 4th order explicit method.
(from :cite:p:`Sof:04`).

.. versionchanged:: 6.3.0

   Made the default 4th order explicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccc}
     0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{2}{5} & \frac{2}{5} & 0 & 0 & 0 & 0 \\
     \frac{3}{5} & -\frac{3}{20} & \frac{3}{4} & 0 & 0 & 0 \\
     1 & \frac{19}{44} & -\frac{15}{44} & \frac{10}{11} & 0 & 0 \\
     1 & \frac{11}{72} & \frac{25}{72} & \frac{25}{72} & \frac{11}{72} & 0 \\
     \hline
     4 & \frac{11}{72} & \frac{25}{72} & \frac{25}{72} & \frac{11}{72} & 0 \\
     3 & \frac{1251515}{8970912} & \frac{3710105}{8970912} & \frac{2519695}{8970912} & \frac{61105}{8970912} & \frac{119041}{747576} \\
   \end{array}

.. figure:: /figs/arkode/sofroniou_spaletta_erk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the Sofroniou-Spaletta method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ZONNEVELD_5_3_4

Accessible via the constant ``ARKODE_ZONNEVELD_5_3_4`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`, or
:c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_ZONNEVELD_5_3_4"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName`, or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`Zon:63`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_SOFRONIOU_SPALETTA_5_3_4`` as the default 4th order explicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccc}
       0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & 0 & \frac{1}{2} & 0 & 0 & 0 \\
       1 & 0 & 0 & 1 & 0 & 0 \\
     \frac{3}{4} & \frac{5}{32} & \frac{7}{32} & \frac{13}{32} & -\frac{1}{32} & 0 \\
     \hline
     4 & \frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6} & 0 \\
     3 & -\frac{1}{2} & \frac{7}{3} & \frac{7}{3} & \frac{13}{6} & -\frac{16}{3}
   \end{array}

.. figure:: /figs/arkode/stab_region_3.png
   :scale: 50 %
   :align: center

   Linear stability region for the Zonneveld method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK436L2SA_ERK_6_3_4

Accessible via the constant ``ARKODE_ARK436L2SA_ERK_6_3_4`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_ARK436L2SA_ERK_6_3_4"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the explicit portion of the ARK4(3)6L[2]SA method from
:cite:p:`KenCarp:03`.

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_ARK437L2SA_ERK_7_3_4`` as the explicit portion of the default 4th order additive method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac12 & \frac12 & 0 & 0 & 0 & 0 & 0 \\
     \frac{83}{250} & \frac{13861}{62500} & \frac{6889}{62500} & 0 & 0 & 0 & 0 \\
     \frac{31}{50} & -\frac{116923316275}{2393684061468} & -\frac{2731218467317}{15368042101831} & \frac{9408046702089}{11113171139209} & 0 & 0 & 0 \\
     \frac{17}{20} & -\frac{451086348788}{2902428689909} & -\frac{2682348792572}{7519795681897} & \frac{12662868775082}{11960479115383} & \frac{3355817975965}{11060851509271} & 0 & 0 \\
     1 & \frac{647845179188}{3216320057751} & \frac{73281519250}{8382639484533} & \frac{552539513391}{3454668386233} & \frac{3354512671639}{8306763924573} & \frac{4040}{17871} & 0 \\
     \hline
     4 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac14 \\
     3 & \frac{4586570599}{29645900160} & 0 & \frac{178811875}{945068544} & \frac{814220225}{1159782912} & -\frac{3700637}{11593932} & \frac{61727}{225920}
   \end{array}

.. figure:: /figs/arkode/stab_region_4.png
   :scale: 50 %
   :align: center

   Linear stability region for the ARK436L2SA-ERK-6-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK437L2SA_ERK_7_3_4

Accessible via the constant ``ARKODE_ARK437L2SA_ERK_7_3_4`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_ARK437L2SA_ERK_7_3_4"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the explicit portion of the default 4th order additive method and the
explicit portion of the ARK4(3)7L[2]SA method from :cite:p:`KenCarp:19`.

.. versionchanged:: 6.3.0

   Made the explicit portion of the default 4th order additive method

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|ccccccc}
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{247}{1000} & \frac{247}{1000} & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{4276536705230}{10142255878289} & \frac{247}{4000} & \frac{2694949928731}{7487940209513} & 0 & 0 & 0 & 0 & 0 \\
        \frac{67}{200} & \frac{464650059369}{8764239774964} & \frac{878889893998}{2444806327765} & -\frac{952945855348}{12294611323341} & 0 & 0 & 0 & 0 \\
        \frac{3}{40} & \frac{476636172619}{8159180917465} & -\frac{1271469283451}{7793814740893} & -\frac{859560642026}{4356155882851} & \frac{1723805262919}{4571918432560} & 0 & 0 & 0 \\
        \frac{7}{10} & \frac{6338158500785}{11769362343261} & -\frac{4970555480458}{10924838743837} & \frac{3326578051521}{2647936831840} & -\frac{880713585975}{1841400956686} & -\frac{1428733748635}{8843423958496} & 0 & 0 \\
        1 & \frac{760814592956}{3276306540349} & \frac{760814592956}{3276306540349} & -\frac{47223648122716}{6934462133451} & \frac{71187472546993}{9669769126921} & -\frac{13330509492149}{9695768672337} & \frac{11565764226357}{8513123442827} & 0 \\
        \hline
        4 & 0 & 0 & \frac{9164257142617}{17756377923965} & -\frac{10812980402763}{74029279521829} & \frac{1335994250573}{5691609445217} & \frac{2273837961795}{8368240463276} & \frac{247}{2000} \\
        3 & 0 & 0 & \frac{4469248916618}{8635866897933} & -\frac{621260224600}{4094290005349} & \frac{696572312987}{2942599194819} & \frac{1532940081127}{5565293938103} & \frac{2441}{20000}
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.

.. figure:: /figs/arkode/stab_region_34.png
   :scale: 50 %
   :align: center

   Linear stability region for the ARK437L2SA-ERK-7-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_SAYFY_ABURUB_6_3_4

Accessible via the constant ``ARKODE_SAYFY_ABURUB_6_3_4`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_SAYFY_ABURUB_6_3_4"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`Sayfy:02`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 & 0 & 0 \\
     1 & -1 & 2 & 0 & 0 & 0 & 0 \\
     1 & \frac{1}{6} & \frac{2}{3} & \frac{1}{6} & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{137}{1000} & \frac{113}{500} & \frac{137}{1000} & 0 & 0 & 0 \\
     1 & \frac{113}{250} & -\frac{113}{125} & -\frac{137}{250} & 0 & 2 & 0 \\
     \hline
     4 & \frac{1}{6} & \frac{1}{3} & \frac{1}{12} & 0 & \frac{1}{3} & \frac{1}{12} \\
     3 & \frac{1}{6} & \frac{2}{3} & \frac{1}{6} & 0 & 0 & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_5.png
   :scale: 50 %
   :align: center

   Linear stability region for the Sayfy-Aburub-6-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_TSITOURAS_7_4_5

Accessible via the constant ``ARKODE_TSITOURAS_7_4_5`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_TSITOURAS_7_4_5"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 5th order explicit method (from :cite:p:`Tsitouras:11`).

.. versionchanged:: 6.3.0

   Added as the default 5th order explicit method

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|ccccccc}
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
         0.161 & 0.161 & 0 & 0 & 0 & 0 & 0 & 0 \\
         0.327 & -0.008480655492356989 & 0.3354806554923570 & 0 & 0 & 0 & 0 & 0 \\
         0.9 & 2.897153057105493 & -6.359448489975075 & 4.362295432869581 & 0 & 0 & 0 & 0 \\
         0.9800255409045097 & 5.325864828439257 & -11.74888356406283 & 7.495539342889836 & -0.09249506636175525 & 0 & 0 & 0 \\
         1 & 5.861455442946420 & -12.92096931784711 & 8.159367898576159 & -0.07158497328140100 & -0.02826905039406838 & 0 & 0 \\
         1 & 0.09646076681806523 & 0.01 & 0.4798896504144996 & 1.379008574103742 & -3.290069515436081 & 2.324710524099774 & 0 \\
         \hline
         5 & 0.09646076681806523 & 0.01 & 0.4798896504144996 & 1.379008574103742 & -3.290069515436081 & 2.324710524099774 & 0 \\
         4 & 0.09352374858189271 & 0.008652883141566368 & 0.4928930991314319 & 1.140235412267858 & -2.329180192439365 & 1.568875049316616 & 0.025
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/tsitouras_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the Tsitouras method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_CASH_KARP_6_4_5

Accessible via the constant ``ARKODE_CASH_KARP_6_4_5`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_CASH_KARP_6_4_5"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`CashKarp:90`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_TSITOURAS_7_4_5`` as the default 5th order explicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{5} & \frac{1}{5} & 0 & 0 & 0 & 0 & 0 \\
     \frac{3}{10} & \frac{3}{40} & \frac{9}{40} & 0 & 0 & 0 & 0 \\
     \frac{3}{5} & \frac{3}{10} & -\frac{9}{10} & \frac{6}{5} & 0 & 0 & 0 \\
     1 & -\frac{11}{54} & \frac{5}{2} & -\frac{70}{27} & \frac{35}{27} & 0 & 0 \\
     \frac{7}{8} & \frac{1631}{55296} & \frac{175}{512} & \frac{575}{13824} & \frac{44275}{110592} & \frac{253}{4096} & 0 \\
     \hline
     5 & \frac{37}{378} & 0 & \frac{250}{621} & \frac{125}{594} & 0 & \frac{512}{1771} \\
     4 & \frac{2825}{27648} & 0 & \frac{18575}{48384} & \frac{13525}{55296} & \frac{277}{14336} & \frac{1}{4}
   \end{array}

.. figure:: /figs/arkode/stab_region_6.png
   :scale: 50 %
   :align: center

   Linear stability region for the Cash-Karp method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_FEHLBERG_6_4_5

Accessible via the constant ``ARKODE_FEHLBERG_6_4_5`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_FEHLBERG_6_4_5"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`Fehlberg:69`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{4} & \frac{1}{4} & 0 & 0 & 0 & 0 & 0 \\
     \frac{3}{8} & \frac{3}{32} & \frac{9}{32} & 0 & 0 & 0 & 0 \\
     \frac{12}{13} & \frac{1932}{2197} & -\frac{7200}{2197} & \frac{7296}{2197} & 0 & 0 & 0 \\
     1 & \frac{439}{216} & -8 & \frac{3680}{513} & -\frac{845}{4104} & 0 & 0 \\
     \frac{1}{2} & -\frac{8}{27} & 2 & -\frac{3544}{2565} & \frac{1859}{4104} & -\frac{11}{40} & 0 \\
     \hline
     5 & \frac{16}{135} & 0 & \frac{6656}{12825} & \frac{28561}{56430} & -\frac{9}{50} & \frac{2}{55} \\
     4 & \frac{25}{216} & 0 & \frac{1408}{2565} & \frac{2197}{4104} & -\frac{1}{5} & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_7.png
   :scale: 50 %
   :align: center

   Linear stability region for the Fehlberg method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_DORMAND_PRINCE_7_4_5

Accessible via the constant ``ARKODE_DORMAND_PRINCE_7_4_5`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_DORMAND_PRINCE_7_4_5"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`DorPri:80`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{5} & \frac{1}{5} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{3}{10} & \frac{3}{40} & \frac{9}{40} & 0 & 0 & 0 & 0 & 0 \\
     \frac{4}{5} & \frac{44}{45} & -\frac{56}{15} & \frac{32}{9} & 0 & 0 & 0 & 0 \\
     \frac{8}{9} & \frac{19372}{6561} & -\frac{25360}{2187} & \frac{64448}{6561} & -\frac{212}{729} & 0 & 0 & 0 \\
     1 & \frac{9017}{3168} & -\frac{355}{33} & \frac{46732}{5247} & \frac{49}{176} & -\frac{5103}{18656} & 0 & 0 \\
     1 & \frac{35}{384} & 0 & \frac{500}{1113} & \frac{125}{192} & -\frac{2187}{6784} & \frac{11}{84} & 0 \\
     \hline
     5 & \frac{35}{384} & 0 & \frac{500}{1113} & \frac{125}{192} & -\frac{2187}{6784} & \frac{11}{84} & 0 \\
     4 & \frac{5179}{57600} & 0 & \frac{7571}{16695} & \frac{393}{640} & -\frac{92097}{339200} & \frac{187}{2100} & \frac{1}{40}
   \end{array}

.. figure:: /figs/arkode/stab_region_8.png
   :scale: 50 %
   :align: center

   Linear stability region for the Dormand-Prince method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK548L2SA_ERK_8_4_5

Accessible via the constant ``ARKODE_ARK548L2SA_ERK_8_4_5`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_ARK548L2SA_ERK_8_4_5"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the explicit portion of the ARK5(4)8L[2]SA method from
:cite:p:`KenCarp:03`.

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_ARK548L2SAb_ERK_8_4_5`` as the explicit portion of the default 5th order additive method

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|cccccccc}
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{41}{100} & \frac{41}{100} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{2935347310677}{11292855782101} & \frac{367902744464}{2072280473677} & \frac{677623207551}{8224143866563} & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{1426016391358}{7196633302097} & \frac{1268023523408}{10340822734521} & 0 & \frac{1029933939417}{13636558850479} & 0 & 0 & 0 & 0 & 0 \\
        \frac{92}{100} & \frac{14463281900351}{6315353703477} & 0 & \frac{66114435211212}{5879490589093} & -\frac{54053170152839}{4284798021562} & 0 & 0 & 0 & 0 \\
        \frac{24}{100} & \frac{14090043504691}{34967701212078} & 0 & \frac{15191511035443}{11219624916014} & -\frac{18461159152457}{12425892160975} & -\frac{281667163811}{9011619295870} & 0 & 0 & 0 \\
        \frac{3}{5} & \frac{19230459214898}{13134317526959} & 0 & \frac{21275331358303}{2942455364971} & -\frac{38145345988419}{4862620318723} & -\frac{1}{8} & -\frac{1}{8} & 0 & 0 \\
        1 & -\frac{19977161125411}{11928030595625} & 0 & -\frac{40795976796054}{6384907823539} & \frac{177454434618887}{12078138498510} & \frac{782672205425}{8267701900261} & -\frac{69563011059811}{9646580694205} & \frac{7356628210526}{4942186776405} & 0 \\
        \hline
        5 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
        4 & -\frac{975461918565}{9796059967033} & 0 & 0 & \frac{78070527104295}{32432590147079} & -\frac{548382580838}{3424219808633} & -\frac{33438840321285}{15594753105479} & \frac{3629800801594}{4656183773603} & \frac{4035322873751}{18575991585200}
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/stab_region_9.png
   :scale: 50 %
   :align: center

   Linear stability region for the explicit ARK-8-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK548L2SAb_ERK_8_4_5

Accessible via the constant ``ARKODE_ARK548L2SAb_ERK_8_4_5`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_ARK548L2SAb_ERK_8_4_5"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the explicit portion of the default 5th order additive method and the
explicit portion of the 5th order ARK5(4)8L[2]SA method from
:cite:p:`KenCarp:19`.

.. versionchanged:: 6.3.0

   Made the explicit portion of the default 5th order additive method

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|cccccccc}
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{4}{9} & \frac{4}{9} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{6456083330201}{8509243623797} & \frac{1}{9} & \frac{1183333538310}{1827251437969} & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{1632083962415}{14158861528103} & \frac{895379019517}{9750411845327} & \frac{477606656805}{13473228687314} & \frac{-112564739183}{9373365219272} & 0 & 0 & 0 & 0 & 0 \\
        \frac{6365430648612}{17842476412687} & \frac{-4458043123994}{13015289567637} & \frac{-2500665203865}{9342069639922} & \frac{983347055801}{8893519644487} & \frac{2185051477207}{2551468980502} & 0 & 0 & 0 & 0 \\
        \frac{18}{25} & \frac{-167316361917}{17121522574472} & \frac{1605541814917}{7619724128744} & \frac{991021770328}{13052792161721} & \frac{2342280609577}{11279663441611} & \frac{3012424348531}{12792462456678} & 0 & 0 & 0 \\
        \frac{191}{200} & \frac{6680998715867}{14310383562358} & \frac{5029118570809}{3897454228471} & \frac{2415062538259}{6382199904604} & \frac{-3924368632305}{6964820224454} & \frac{-4331110370267}{15021686902756} & \frac{-3944303808049}{11994238218192} & 0 & 0 \\
        1 & \frac{2193717860234}{3570523412979} & \frac{2193717860234}{3570523412979} & \frac{5952760925747}{18750164281544} & \frac{-4412967128996}{6196664114337} & \frac{4151782504231}{36106512998704} & \frac{572599549169}{6265429158920} & \frac{-457874356192}{11306498036315} & 0 \\
        \hline
        5 & 0 & 0 & \frac{3517720773327}{20256071687669} & \frac{4569610470461}{17934693873752} & \frac{2819471173109}{11655438449929} & \frac{3296210113763}{10722700128969} & \frac{-1142099968913}{5710983926999} & \frac{2}{9} \\
        4 & 0 & 0 & \frac{520639020421}{8300446712847} & \frac{4550235134915}{17827758688493} & \frac{1482366381361}{6201654941325} & \frac{5551607622171}{13911031047899} & \frac{-5266607656330}{36788968843917} & \frac{1074053359553}{5740751784926}
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/stab_region_35.png
   :scale: 50 %
   :align: center

   Linear stability region for the ARK548L2SAb-ERK-8-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_VERNER_9_5_6

Accessible via the constant ``ARKODE_VERNER_9_5_6`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_VERNER_9_5_6"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 6th order explicit method
(method IIIXb-6(5) from :cite:p:`Ver:10`).

.. versionchanged:: 6.3.0

   Made the default 6th order explicit method

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|ccccccccc}
        0           & 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{3}{50}        & \frac{3}{50}& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{1439}{15000}  & \frac{519479}{27000000}& \frac{2070721}{27000000}& 0& 0& 0& 0& 0& 0& 0\\
        \frac{1439}{10000}  & \frac{1439}{40000}& 0& \frac{4317}{40000}& 0& 0& 0& 0& 0& 0\\
        \frac{4973}{10000}  & \frac{109225017611}{82828840000}& 0& -\frac{417627820623}{82828840000}& \frac{43699198143}{10353605000}& 0& 0& 0& 0& 0\\
        \frac{389}{400}     & -\frac{8036815292643907349452552172369}{191934985946683241245914401600}& 0& \frac{246134619571490020064824665}{1543816496655405117602368}& -\frac{13880495956885686234074067279}{113663489566254201783474344}& \frac{755005057777788994734129}{136485922925633667082436}& 0& 0& 0& 0\\
        \frac{1999}{2000}   & -\frac{1663299841566102097180506666498880934230261}{30558424506156170307020957791311384232000}& 0& \frac{130838124195285491799043628811093033}{631862949514135618861563657970240}& -\frac{3287100453856023634160618787153901962873}{20724314915376755629135711026851409200}& \frac{2771826790140332140865242520369241}{396438716042723436917079980147600}& -\frac{1799166916139193}{96743806114007800}& 0& 0& 0\\
        1           & -\frac{832144750039369683895428386437986853923637763}{15222974550069600748763651844667619945204887}& 0& \frac{818622075710363565982285196611368750}{3936576237903728151856072395343129}& -\frac{9818985165491658464841194581385463434793741875}{61642597962658994069869370923196463581866011}& \frac{31796692141848558720425711042548134769375}{4530254033500045975557858016006308628092}& -\frac{14064542118843830075}{766928748264306853644}& -\frac{1424670304836288125}{2782839104764768088217}& 0& 0\\
        1           &   \frac{382735282417}{11129397249634}& 0& 0& \frac{5535620703125000}{21434089949505429}& \frac{13867056347656250}{32943296570459319}& \frac{626271188750}{142160006043}& -\frac{51160788125000}{289890548217}& \frac{163193540017}{946795234}& 0\\
        \hline
        6           & \frac{382735282417}{11129397249634}& 0& 0& \frac{5535620703125000}{21434089949505429}& \frac{13867056347656250}{32943296570459319}& \frac{626271188750}{142160006043}& -\frac{51160788125000}{289890548217}& \frac{163193540017}{946795234}& 0 \\
        5           & \frac{273361583}{5567482366}& 0& 0& \frac{1964687500000}{8727630165387}& \frac{596054687500}{1269637976277}& \frac{12740367500}{15795556227}& 0& -\frac{4462730789736252634813752317}{7350663039626676022821734166}& \frac{441454562788983500}{7763730504400359099}
      \end{array}


.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.



.. figure:: /figs/arkode/v65b_erk_stab_region.png
   :scale: 75 %
   :align: center

   Linear stability region for the Verner-9-5-6 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_VERNER_8_5_6

Accessible via the constant ``ARKODE_VERNER_8_5_6`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_VERNER_8_5_6"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`HEJ:76`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_VERNER_9_5_6`` as the default 6th order explicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{6} & \frac{1}{6} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{4}{15} & \frac{4}{75} & \frac{16}{75} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{2}{3} & \frac{5}{6} & -\frac{8}{3} & \frac{5}{2} & 0 & 0 & 0 & 0 & 0 \\
     \frac{5}{6} & -\frac{165}{64} & \frac{55}{6} & -\frac{425}{64} & \frac{85}{96} & 0 & 0 & 0 & 0 \\
     1 & \frac{12}{5} & -8 & \frac{4015}{612} & -\frac{11}{36} & \frac{88}{255} & 0 & 0 & 0 \\
     \frac{1}{15} & -\frac{8263}{15000} & \frac{124}{75} & -\frac{643}{680} & -\frac{81}{250} & \frac{2484}{10625} & 0 & 0 & 0 \\
     1 & \frac{3501}{1720} & -\frac{300}{43} & \frac{297275}{52632} & -\frac{319}{2322} & \frac{24068}{84065} & 0 & \frac{3850}{26703} & 0 \\
     \hline
     6 & \frac{3}{40} & 0 & \frac{875}{2244} & \frac{23}{72} & \frac{264}{1955} & 0 & \frac{125}{11592} & \frac{43}{616} \\
     5 & \frac{13}{160} & 0 & \frac{2375}{5984} & \frac{5}{16} & \frac{12}{85} & \frac{3}{44} & 0 & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_10.png
   :scale: 50 %
   :align: center

   Linear stability region for the Verner-8-5-6 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_VERNER_10_6_7

Accessible via the constant ``ARKODE_VERNER_10_6_7`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_VERNER_10_6_7"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 7th order explicit method (from :cite:p:`Ver:10`).

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|cccccccccc}
        0                      & 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{1}{200}                  & \frac{1}{200}& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{49}{450}                 & -\frac{4361}{4050}& \frac{2401}{2025}& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{49}{300}                 & \frac{49}{1200}& 0& \frac{49}{400}& 0& 0& 0& 0& 0& 0& 0\\
        \frac{911}{2000}               & \frac{2454451729}{3841600000}& 0& -\frac{9433712007}{3841600000}& \frac{4364554539}{1920800000}& 0& 0& 0& 0& 0& 0\\
        \frac{3480084980}{5709648941}  & -\frac{6187101755456742839167388910402379177523537620}{2324599620333464857202963610201679332423082271}& 0& \frac{27569888999279458303270493567994248533230000}{2551701010245296220859455115479340650299761}& -\frac{37368161901278864592027018689858091583238040000}{4473131870960004275166624817435284159975481033}& \frac{1392547243220807196190880383038194667840000000}{1697219131380493083996999253929006193143549863}& 0& 0& 0& 0& 0\\
        \frac{221}{250}                & \frac{11272026205260557297236918526339}{1857697188743815510261537500000}& 0& -\frac{48265918242888069}{1953194276993750}& \frac{26726983360888651136155661781228}{1308381343805114800955157615625}& -\frac{2090453318815827627666994432}{1096684189897834170412307919}& \frac{1148577938985388929671582486744843844943428041509}{1141532118233823914568777901158338927629837500000}& 0& 0& 0& 0\\
        \frac{37}{40}                  & \frac{1304457204588839386329181466225966641}{108211771565488329642169667802016000}& 0& -\frac{1990261989751005}{40001418792832}& \frac{2392691599894847687194643439066780106875}{58155654089143548047476915856270826016}& -\frac{1870932273351008733802814881998561250}{419326053051486744762255151208232123}& \frac{1043329047173803328972823866240311074041739158858792987034783181}{510851127745017966999893975119259285040213723744255237522144000}& -\frac{311918858557595100410788125}{3171569057622789618800376448}& 0& 0& 0\\
        1                      & \frac{17579784273699839132265404100877911157}{1734023495717116205617154737841023480}& 0& -\frac{18539365951217471064750}{434776548575709731377}& \frac{447448655912568142291911830292656995992000}{12511202807447096607487664209063950964109}& -\frac{65907597316483030274308429593905808000000}{15158061430635748897861852383197382130691}& \frac{273847823027445129865693702689010278588244606493753883568739168819449761}{136252034448398939768371761610231099586032870552034688235302796640584360}& \frac{694664732797172504668206847646718750}{1991875650119463976442052358853258111}& -\frac{19705319055289176355560129234220800}{72595753317320295604316217197876507}& 0& 0\\
        1                      & -\frac{511858190895337044664743508805671}{11367030248263048398341724647960}& 0& \frac{2822037469238841750}{15064746656776439}& -\frac{23523744880286194122061074624512868000}{152723005449262599342117017051789699}& \frac{10685036369693854448650967542704000000}{575558095977344459903303055137999707}& -\frac{6259648732772142303029374363607629515525848829303541906422993}{876479353814142962817551241844706205620792843316435566420120}& \frac{17380896627486168667542032602031250}{13279937889697320236613879977356033}& 0& 0& 0\\
        \hline
        7                      & \frac{96762636172307789}{2051985304794103980}& 0& 0& \frac{312188947591288252500000}{1212357694274963646019729}& \frac{13550580884964304000000000000}{51686919683339547115937980629}& \frac{72367769693133178898676076432831566019684378142853445230956642801}{475600216991873963561768100160364792981629064220601844848928537580}& \frac{1619421054120605468750}{3278200730370057108183}& -\frac{66898316144057728000}{227310933007074849597}& \frac{181081444637946577}{2226845467039736466}& 0 \\
        6                      & \frac{117807213929927}{2640907728177740}& 0& 0& \frac{4758744518816629500000}{17812069906509312711137}& \frac{1730775233574080000000000}{7863520414322158392809673}& \frac{2682653613028767167314032381891560552585218935572349997}{12258338284789875762081637252125169126464880985167722660}& \frac{40977117022675781250}{178949401077111131341}& 0& 0& \frac{2152106665253777}{106040260335225546}
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/v76_erk_stab_region.png
   :scale: 75 %
   :align: center

   Linear stability region for the Verner-10-6-7 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_VERNER_13_7_8

Accessible via the constant ``ARKODE_VERNER_13_7_8`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_VERNER_13_7_8"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 8th order explicit method
(method IIIX-8(7) from :cite:p:`Ver:10`).

.. versionchanged:: 6.3.0

   Made the default 8th order explicit method


.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|ccccccccccccc}
        0                           & 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{1}{20}                        & \frac{1}{20}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{341}{3200}                    & -\frac{7161}{1024000}& \frac{116281}{1024000}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{1023}{6400}                   & \frac{1023}{25600}& 0& \frac{3069}{25600}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{39}{100}                      & \frac{4202367}{11628100}& 0& -\frac{3899844}{2907025}& \frac{3982992}{2907025}& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{93}{200}                      & \frac{5611}{114400}& 0& 0& \frac{31744}{135025}& \frac{923521}{5106400}& 0& 0& 0& 0& 0& 0& 0& 0\\
        \frac{31}{200}                      & \frac{21173}{343200}& 0& 0& \frac{8602624}{76559175}& -\frac{26782109}{689364000}& \frac{5611}{283500}& 0& 0& 0& 0& 0& 0& 0\\
        \frac{943}{1000}                    & -\frac{1221101821869329}{690812928000000}& 0& 0& -\frac{125}{2}& -\frac{1024030607959889}{168929280000000}& \frac{1501408353528689}{265697280000000}& \frac{6070139212132283}{92502016000000}& 0& 0& 0& 0& 0& 0\\
        \frac{7067558016280}{7837150160667} &  -\frac{1472514264486215803881384708877264246346044433307094207829051978044531801133057155}{1246894801620032001157059621643986024803301558393487900440453636168046069686436608}& 0& 0& -\frac{5172294311085668458375175655246981230039025336933699114138315270772319372469280000}{124619381004809145897278630571215298365257079410236252921850936749076487132995191}& -\frac{12070679258469254807978936441733187949484571516120469966534514296406891652614970375}{2722031154761657221710478184531100699497284085048389015085076961673446140398628096}& \frac{780125155843893641323090552530431036567795592568497182701460674803126770111481625}{183110425412731972197889874507158786859226102980861859505241443073629143100805376}& \frac{664113122959911642134782135839106469928140328160577035357155340392950009492511875}{15178465598586248136333023107295349175279765150089078301139943253016877823170816}& \frac{10332848184452015604056836767286656859124007796970668046446015775000000}{1312703550036033648073834248740727914537972028638950165249582733679393783}& 0& 0& 0& 0& 0\\
        \frac{909}{1000}                    & -\frac{29055573360337415088538618442231036441314060511}{22674759891089577691327962602370597632000000000}& 0& 0& -\frac{20462749524591049105403365239069}{454251913499893469596231268750}& -\frac{180269259803172281163724663224981097}{38100922558256871086579832832000000}& \frac{21127670214172802870128286992003940810655221489}{4679473877997892906145822697976708633673728000}& \frac{318607235173649312405151265849660869927653414425413}{6714716715558965303132938072935465423910912000000}& \frac{212083202434519082281842245535894}{20022426044775672563822865371173879}& -\frac{2698404929400842518721166485087129798562269848229517793703413951226714583}{469545674913934315077000442080871141884676035902717550325616728175875000000}& 0& 0& 0& 0\\
        \frac{47}{50}                       & -\frac{2342659845814086836951207140065609179073838476242943917}{1358480961351056777022231400139158760857532162795520000}& 0& 0& -\frac{996286030132538159613930889652}{16353068885996164905464325675}& -\frac{26053085959256534152588089363841}{4377552804565683061011299942400}& \frac{20980822345096760292224086794978105312644533925634933539}{3775889992007550803878727839115494641972212962174156800}& \frac{890722993756379186418929622095833835264322635782294899}{13921242001395112657501941955594013822830119803764736}& \frac{161021426143124178389075121929246710833125}{10997207722131034650667041364346422894371443}& \frac{300760669768102517834232497565452434946672266195876496371874262392684852243925359864884962513}{4655443337501346455585065336604505603760824779615521285751892810315680492364106674524398280000}& -\frac{31155237437111730665923206875}{392862141594230515010338956291}& 0& 0& 0\\
        1                           & -\frac{2866556991825663971778295329101033887534912787724034363}{868226711619262703011213925016143612030669233795338240}& 0& 0& -\frac{16957088714171468676387054358954754000}{143690415119654683326368228101570221}& -\frac{4583493974484572912949314673356033540575}{451957703655250747157313034270335135744}& \frac{2346305388553404258656258473446184419154740172519949575}{256726716407895402892744978301151486254183185289662464}& \frac{1657121559319846802171283690913610698586256573484808662625}{13431480411255146477259155104956093505361644432088109056}& \frac{345685379554677052215495825476969226377187500}{74771167436930077221667203179551347546362089}& -\frac{3205890962717072542791434312152727534008102774023210240571361570757249056167015230160352087048674542196011}{947569549683965814783015124451273604984657747127257615372449205973192657306017239103491074738324033259120}& \frac{40279545832706233433100438588458933210937500}{8896460842799482846916972126377338947215101}& -\frac{6122933601070769591613093993993358877250}{1050517001510235513198246721302027675953}& 0& 0\\
        1                           & -\frac{618675905535482500672800859344538410358660153899637}{203544282118214047100119475340667684874292102389760}& 0& 0& -\frac{4411194916804718600478400319122931000}{40373053902469967450761491269633019}& -\frac{16734711409449292534539422531728520225}{1801243715290088669307203927210237952}& \frac{135137519757054679098042184152749677761254751865630525}{16029587794486289597771326361911895112703716593983488}& \frac{38937568367409876012548551903492196137929710431584875}{340956454090191606099548798001469306974758443147264}& -\frac{6748865855011993037732355335815350667265625}{7002880395717424621213565406715087764770357}& -\frac{1756005520307450928195422767042525091954178296002788308926563193523662404739779789732685671}{348767814578469983605688098046186480904607278021030540735333862087061574934154942830062320}& \frac{53381024589235611084013897674181629296875}{8959357584795694524874969598508592944141}& 0& 0& 0\\
        \hline
        8                           & \frac{44901867737754616851973}{1014046409980231013380680}& 0& 0& 0& 0& \frac{791638675191615279648100000}{2235604725089973126411512319}& \frac{3847749490868980348119500000}{15517045062138271618141237517}& -\frac{13734512432397741476562500000}{875132892924995907746928783}& \frac{12274765470313196878428812037740635050319234276006986398294443554969616342274215316330684448207141}{489345147493715517650385834143510934888829280686609654482896526796523353052166757299452852166040}& -\frac{9798363684577739445312500000}{308722986341456031822630699}& \frac{282035543183190840068750}{12295407629873040425991}& -\frac{306814272936976936753}{1299331183183744997286}& 0\\
        7                           & \frac{10835401739407019406577}{244521829356935137978320}& 0& 0& 0& 0& \frac{13908189778321895491375000}{39221135527894265375640567}& \frac{73487947527027243487625000}{296504045773342769773399443}& \frac{68293140641257649609375000}{15353208647806945749946119}& \frac{22060647948996678611017711379974578860522018208949721559448560203338437626022142776381}{1111542009262325874512959185795727215759010577565736079641376621381577236680929558640}& -\frac{547971229495642458203125000}{23237214025700991642563601}& 0& 0& -\frac{28735456870978964189}{79783493704265043693}
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/v87_erk_stab_region.png
   :scale: 75 %
   :align: center

   Linear stability region for the Verner-13-7-8 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_FEHLBERG_13_7_8

Accessible via the constant ``ARKODE_FEHLBERG_13_7_8`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_FEHLBERG_13_7_8"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
(from :cite:p:`Butcher:08`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_VERNER_13_7_8`` as the default 8th order explicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccccccccccc}
     0&   0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{2}{27}&   \frac{2}{27}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{9}&   \frac{1}{36}& \frac{1}{12}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{6}&   \frac{1}{24}& 0& \frac{1}{8}& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{5}{12}&   \frac{5}{12}& 0& -\frac{25}{16}& \frac{25}{16}& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{2}&   \frac{1}{20}& 0& 0& \frac{1}{4}& \frac{1}{5}& 0& 0& 0& 0& 0& 0& 0& 0\\
     \frac{5}{6}&   -\frac{25}{108}& 0& 0& \frac{125}{108}& -\frac{65}{27}& \frac{125}{54}& 0& 0& 0& 0& 0& 0& 0\\
     \frac{1}{6}&   \frac{31}{300}& 0& 0& 0& \frac{61}{225}& -\frac{2}{9}& \frac{13}{900}& 0& 0& 0& 0& 0& 0\\
     \frac{2}{3}&   2& 0& 0& -\frac{53}{6}& \frac{704}{45}& -\frac{107}{9}& \frac{67}{90}& 3& 0& 0& 0& 0& 0\\
     \frac{1}{3}&   -\frac{91}{108}& 0& 0& \frac{23}{108}& -\frac{976}{135}& \frac{311}{54}& -\frac{19}{60}& \frac{17}{6}& -\frac{1}{12}& 0& 0& 0& 0\\
     1&   \frac{2383}{4100}& 0& 0& -\frac{341}{164}& \frac{4496}{1025}& -\frac{301}{82}& \frac{2133}{4100}& \frac{45}{82}& \frac{45}{164}& \frac{18}{41}& 0& 0& 0\\
     0&   \frac{3}{205}& 0& 0& 0& 0& -\frac{6}{41}& -\frac{3}{205}& -\frac{3}{41}& \frac{3}{41}& \frac{6}{41}& 0& 0& 0\\
     1&   -\frac{1777}{4100}& 0& 0& -\frac{341}{164}& \frac{4496}{1025}& -\frac{289}{82}& \frac{2193}{4100}& \frac{51}{82}& \frac{33}{164}& \frac{12}{41}& 0& 1& 0\\
     \hline
     8& 0& 0& 0& 0& 0& \frac{34}{105}& \frac{9}{35}& \frac{9}{35}& \frac{9}{280}& \frac{9}{280}& 0& \frac{41}{840}& \frac{41}{840} \\
     7& \frac{41}{840}& 0& 0& 0& 0& \frac{34}{105}& \frac{9}{35}& \frac{9}{35}& \frac{9}{280}& \frac{9}{280}& \frac{41}{840}& 0& 0
   \end{array}


.. figure:: /figs/arkode/stab_region_23.png
   :scale: 50 %
   :align: center

   Linear stability region for the Fehlberg-13-7-8 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_VERNER_16_8_9

Accessible via the constant ``ARKODE_VERNER_16_8_9`` to
:c:func:`ARKStepSetTableNum`, :c:func:`ERKStepSetTableNum`
or :c:func:`ARKodeButcherTable_LoadERK`.
Accessible via the string ``"ARKODE_VERNER_16_8_9"`` to
:c:func:`ARKStepSetTableName`, :c:func:`ERKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadERKByName`.
This is the default 9th order explicit method (from :cite:p:`Ver:10`).


.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|cccccccccccccccc}
       0                                             & 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.03462                                     & 0.03462& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.09702435063878044594828361677100617517633 & -0.0389335438857287327017042687229284478532& 0.1359578945245091786499878854939346230295& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.1455365259581706689224254251565092627645    & 0.03638413148954266723060635628912731569111& 0& 0.1091523944686280016918190688673819470733& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.561                                         & 2.025763914393969636805657604282571047511& 0& -7.638023836496292020387602153091964592952& 6.173259922102322383581944548809393545442& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.2290079115904850126662751771814700052182    & 0.05112275589406060872792270881648288397197& 0& 0& 0.1770823794555021537929910813839068684087& 0.00080277624092225014536138698108025283759& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.5449920884095149873337248228185299947818    & 0.1316006357975216279279871693164256985334& 0& 0& -0.2957276252669636417685183174672273730699& 0.0878137803564295237421124704053886667082& 0.6213052975225274774321435005639430026100& 0& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.645                                         & 0.07166666666666666666666666666666666666667& 0& 0& 0& 0& 0.3305533578915319409260346730051472207728& 0.2427799754418013924072986603281861125606& 0& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.48375                                       & 0.071806640625& 0& 0& 0& 0& 0.3294380283228177160744825466257672816401& 0.1165190029271822839255174533742327183599& -0.034013671875& 0& 0& 0& 0& 0& 0& 0& 0\\
       0.06757                                     & 0.04836757646340646986611287718844085773549& 0& 0& 0& 0& 0.03928989925676163974333190042057047002852& 0.1054740945890344608263649267140088017604& -0.02143865284648312665982642293830533996214& -0.1041229174627194437759832813847147895623& 0& 0& 0& 0& 0& 0& 0\\
       0.25                                          & -0.02664561487201478635337289243849737340534& 0& 0& 0& 0& 0.03333333333333333333333333333333333333333& -0.1631072244872467239162704487554706387141& 0.03396081684127761199487954930015522928244& 0.1572319413814626097110769806810024118077& 0.2152267478031879552303534778794770376960& 0& 0& 0& 0& 0& 0\\
       0.6590650618730998549405331618649220295334    & 0.03689009248708622334786359863227633989718& 0& 0& 0& 0& -0.1465181576725542928653609891758501156785& 0.2242577768172024345345469822625833796001& 0.02294405717066072637090897902753790803034& -0.0035850052905728761357394424889330334334& 0.08669223316444385506869203619044453906053& 0.4383840651968337846196219974168630120572& 0& 0& 0& 0& 0\\
       0.8206                                        & -0.4866012215113340846662212357570395295088& 0& 0& 0& 0& -6.304602650282852990657772792012007122988& -0.281245618289472564778284183790118418111& -2.679019236219849057687906597489223155566& 0.518815663924157511565311164615012522024& 1.365353187603341710683633635235238678626& 5.885091088503946585721274891680604830712& 2.802808786272062889819965117517532194812& 0& 0& 0& 0\\
       0.9012                                        & 0.4185367457753471441471025246471931649633& 0& 0& 0& 0& 6.724547581906459363100870806514855026676& -0.425444280164611790606983409697113064616& 3.343279153001265577811816947557982637749& 0.617081663117537759528421117507709784737& -0.929966123939932833937749523988800852013& -6.099948804751010722472962837945508844846& -3.002206187889399044804158084895173690015& 0.2553202529443445472336424602988558373637& 0& 0& 0\\
       1                                             & -0.779374086122884664644623040843840506343& 0& 0& 0& 0& -13.93734253810777678786523664804936051203& 1.252048853379357320949735183924200895136& -14.69150040801686878191527989293072091588& -0.494705058533141685655191992136962873577& 2.242974909146236657906984549543692874755& 13.36789380382864375813864978592679139881& 14.39665048665068644512236935340272139005& -0.7975813331776800379127866056663258667437& 0.4409353709534277758753793068298041158235& 0& 0\\
       1                                             & 2.058051337466886442151242368989994043993& 0& 0& 0& 0& 22.35793772796803295519317565842520212899& 0.90949810997556332745009198137971890783& 35.89110098240264104710550686568482456493& -3.442515027624453437985000403608480262211& -4.865481358036368826566013387928704014496& -18.90980381354342625688427480879773032857& -34.26354448030451782929251177395134170515& 1.264756521695642578827783499806516664686& 0& 0& 0\\
       \hline
       9                                             & 0.01461197685842315252051541915018784713459& 0& 0& 0& 0& 0& 0& -0.3915211862331339089410228267288242030810& 0.2310932500289506415909675644868993669908& 0.1274766769992852382560589467488989175618& 0.2246434176204157731566981937082069688984& 0.5684352689748512932705226972873692126743& 0.05825871557215827200814768021863420902155& 0.1364317403482215641609022744494239843327& 0.03057013983082797397721005067920369646664& 0\\
       8                                             & 0.01996996514886773085518508418098868756464& 0& 0& 0& 0& 0& 0& 2.191499304949330054530747099310837524864& 0.08857071848208438030833722031786358862953& 0.1140560234865965622484956605091432032674& 0.2533163805345107065564577734569651977347& -2.056564386240941011158999594595981300493& 0.3408096799013119935160094894224543812830& 0& 0& 0.04834231373823958314376726739772871714902
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/v98_erk_stab_region.png
   :scale: 75 %
   :align: center

   Linear stability region for the Verner-16-8-9 method.  The method's
   region is outlined in blue; the embedding's region is in red.





.. _Butcher.implicit:

Implicit Butcher tables
---------------------------

In the category of diagonally implicit Runge--Kutta methods, ARKODE includes
methods that have orders 2 through 5, with embeddings that are of orders 1
through 4. ARKODE's diagonally-implicit Butcher tables are
provided in the enumeration

.. c:enum:: ARKODE_DIRKTableID

with values specified in :numref:`ARKODE.Butcher.DIRK_properties`.

.. _ARKODE.Butcher.DIRK_properties:
.. table:: Implicit Butcher tables. The default method for each order is marked
   with an asterisk (*).

   +-------------------------------------------------+--------+----------------+-------+
   | Method ID                                       | Stages | Embedded Order | Order |
   +=================================================+========+================+=======+
   | :c:enumerator:`ARKODE_BACKWARD_EULER_1_1`       | 1      | ---            | 1*    |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK2_DIRK_3_1_2`          | 3      | 1              | 2*    |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_SDIRK_2_1_2`              | 2      | 1              | 2     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_IMPLICIT_MIDPOINT_1_2`    | 1      | ---            | 2     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_IMPLICIT_TRAPEZOIDAL_2_2` | 2      | ---            | 2     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_BILLINGTON_3_3_2`         | 3      | 3              | 2     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_TRBDF2_3_3_2`             | 3      | 3              | 2     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK325L2SA_5_2_3`      | 5      | 2              | 3*    |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK324L2SA_4_2_3`      | 4      | 2              | 3     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK32I5L2SA_5_2_3`     | 5      | 2              | 3     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_KVAERNO_4_2_3`            | 4      | 2              | 3     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK324L2SA_DIRK_4_2_3`    | 4      | 2              | 3     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK436L2SA_6_3_4`      | 6      | 3              | 4*    |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_CASH_5_2_4`               | 5      | 2              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_CASH_5_3_4`               | 5      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_SDIRK_5_3_4`              | 5      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_KVAERNO_5_3_4`            | 5      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK436L2SA_DIRK_6_3_4`    | 6      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK437L2SA_DIRK_7_3_4`    | 7      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK43I6L2SA_6_3_4`     | 6      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_QESDIRK436L2SA_6_3_4`     | 6      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK437L2SA_7_3_4`      | 7      | 3              | 4     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK547L2SA2_7_4_5`     | 7      | 4              | 5*    |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_KVAERNO_7_4_5`            | 7      | 4              | 5     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK548L2SA_DIRK_8_4_5`    | 8      | 4              | 5     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK548L2SAb_DIRK_8_4_5`   | 8      | 4              | 5     |
   +-------------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ESDIRK547L2SA_7_4_5`      | 7      | 4              | 5     |
   +-------------------------------------------------+--------+----------------+-------+


.. c:enumerator:: ARKODE_BACKWARD_EULER_1_1

Accessible via the constant ``ARKODE_BACKWARD_EULER_1_1`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_BACKWARD_EULER_1_1"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the default 1st order implicit method.  The method is A-, L-, and B-stable.

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|c}
     1 & 1 \\
     \hline
     1 & 1
   \end{array}

.. figure:: /figs/arkode/backward_euler_dirk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the backward Euler method.


.. c:enumerator:: ARKODE_ARK2_DIRK_3_1_2

Accessible via the constant ``ARKODE_ARK2_DIRK_3_1_2`` to
:c:func:`ARKStepSetTableNum`, or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ARK2_DIRK_3_1_2"`` to
:c:func:`ARKStepSetTableName`, or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the default 2nd order implicit method and the implicit portion of the
default 2nd order additive method
(the implicit portion of the ARK2 method from :cite:p:`giraldo2013implicit`).

.. versionchanged:: 6.3.0

   Made the default 2nd order implicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccc}
     0            & 0                      & 0                      & 0 \\
     2 - \sqrt{2} & 1 - \frac{1}{\sqrt{2}} & 1 - \frac{1}{\sqrt{2}} & 0 \\
     1            & \frac{1}{2\sqrt{2}}    & \frac{1}{2\sqrt{2}}    & 1 - \frac{1}{\sqrt{2}} \\
     \hline
     2 & \frac{1}{2\sqrt{2}}    & \frac{1}{2\sqrt{2}}    & 1 - \frac{1}{\sqrt{2}} \\
     1 & \frac{4 - \sqrt{2}}{8} & \frac{4 - \sqrt{2}}{8} & \frac{1}{2\sqrt{2}}    \\
   \end{array}

.. figure:: /figs/arkode/ark2_dirk_stab_region.png
   :scale: 65 %
   :align: center

   Linear stability region for the ARK2-DIRK method. The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_SDIRK_2_1_2

Accessible via the constant ``ARKODE_SDIRK_2_1_2`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_SDIRK_2_1_2"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A- and B-stable.

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_ARK2_DIRK_3_1_2`` as the default 2nd order implicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cc}
     1 & 1 & 0 \\
     0 & -1 & 1 \\
     \hline
     2 & \frac{1}{2} & \frac{1}{2} \\
     1 & 1 & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_11.png
   :scale: 50 %
   :align: center

   Linear stability region for the SDIRK-2-1-2 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_IMPLICIT_MIDPOINT_1_2

Accessible via the constant ``ARKODE_IMPLICIT_MIDPOINT_1_2`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_IMPLICIT_MIDPOINT_1_2"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
The method is A- and B-stable.

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|c}
     \frac{1}{2} & \frac{1}{2} \\
     \hline
     2 & 1
   \end{array}

.. figure:: /figs/arkode/implicit_midpoint_dirk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the implicit midpoint method.


.. c:enumerator:: ARKODE_IMPLICIT_TRAPEZOIDAL_2_2

Accessible via the constant ``ARKODE_IMPLICIT_TRAPEZOIDAL_2_2`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_IMPLICIT_TRAPEZOIDAL_2_2"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
The method is A-stable.

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cc}
     0 & 0 & 0 \\
     1 & \frac{1}{2} & \frac{1}{2} \\
     \hline
     2 & \frac{1}{2} & \frac{1}{2}
   \end{array}

.. figure:: /figs/arkode/implicit_trapezoidal_dirk_stab_region.png
   :scale: 50 %
   :align: center

   Linear stability region for the implicit trapezoidal method.


.. c:enumerator:: ARKODE_BILLINGTON_3_3_2

Accessible via the constant ``ARKODE_BILLINGTON_3_3_2`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_BILLINGTON_3_3_2"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Here, the higher-order embedding is less stable than the lower-order method
(from :cite:p:`Billington:83`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccc}
     1 - \frac{1}{\sqrt{2}} & 1 - \frac{1}{\sqrt{2}} & 0 & 0 \\
     \frac{27}{\sqrt{2}} - 18 & 14 \sqrt{2} - 19 & 1 - \frac{1}{\sqrt{2}} & 0 \\
     2 - \frac{1}{\sqrt{2}} & \frac{53 - 5 \sqrt{2}}{62} & \frac{9 + 5 \sqrt{2}}{62} & 1 - \frac{1}{\sqrt{2}} \\
     \hline
     2 & \frac{53 - 5 \sqrt{2}}{62} & \frac{9 + 5 \sqrt{2}}{62} & 0 \\
     3 & \frac{263 - 95 \sqrt{2}}{186} & \frac{47 + 33 \sqrt{2}}{186} & \frac{\sqrt{2} - 2}{3}
   \end{array}

.. figure:: /figs/arkode/stab_region_12.png
   :scale: 50 %
   :align: center

   Linear stability region for the Billington method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_TRBDF2_3_3_2

Accessible via the constant ``ARKODE_TRBDF2_3_3_2`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_TRBDF2_3_3_2"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
As with Billington, here the higher-order embedding is less stable than the
lower-order method (from :cite:p:`Bank:85`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccc}
     0 & 0 & 0 & 0 \\
     2-\sqrt{2} & \frac{2-\sqrt{2}}{2} & \frac{2-\sqrt{2}}{2} & 0 \\
     1 & \frac{\sqrt{2}}{4} & \frac{\sqrt{2}}{4} & \frac{2-\sqrt{2}}{2} \\
     \hline
     2 & \frac{\sqrt{2}}{4} & \frac{\sqrt{2}}{4} & \frac{2-\sqrt{2}}{2} \\
     3 & \frac{1-\frac{\sqrt{2}}{4}}{3} & \frac{\frac{3\sqrt{2}}{4}+1}{3} & \frac{2-\sqrt{2}}{6}
   \end{array}

.. figure:: /figs/arkode/stab_region_13.png
   :scale: 50 %
   :align: center

   Linear stability region for the TRBDF2 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK325L2SA_5_2_3

Accessible via the constant ``ARKODE_ESDIRK325L2SA_5_2_3`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK325L2SA_5_2_3"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the default 3rd order implicit method and the ESDIRK3(2)5L[2]SA method
from :cite:p:`KenCarp:16`.
Both the method and embedding are A- and L-stable.

.. versionchanged:: 6.3.0

   Made the default 3rd order implicit method

.. figure:: /figs/arkode/stab_region_26.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK325L2SA-5-2-3 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK324L2SA_4_2_3

Accessible via the constant ``ARKODE_ESDIRK324L2SA_4_2_3`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK324L2SA_4_2_3"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the ESDIRK3(2)4L[2]SA method from :cite:p:`KenCarp:19b`.
Both the method and embedding are A- and L-stable.

.. figure:: /figs/arkode/stab_region_25.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK324L2SA-4-2-3 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK32I5L2SA_5_2_3

Accessible via the constant ``ARKODE_ESDIRK32I5L2SA_5_2_3`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK32I5L2SA_5_2_3"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the ESDIRK3(2I)5L[2]SA method from :cite:p:`KenCarp:16`.
Both the method and embedding are A- and L-stable.

.. figure:: /figs/arkode/stab_region_27.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK32I5L2SA-5-2-3 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_KVAERNO_4_2_3

Accessible via the constant ``ARKODE_KVAERNO_4_2_3`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_KVAERNO_4_2_3"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A-stable; additionally the method is L-stable
(from :cite:p:`Kva:04`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     0.871733043 & 0.4358665215 & 0.4358665215 & 0 & 0 \\
     1 & 0.490563388419108 & 0.073570090080892 & 0.4358665215 & 0 \\
     1 & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
     \hline
     3 & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
     2 & 0.490563388419108 & 0.073570090080892 & 0.4358665215 & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_14.png
   :scale: 50 %
   :align: center

   Linear stability region for the Kvaerno-4-2-3 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK324L2SA_DIRK_4_2_3

Accessible via the constant ``ARKODE_ARK324L2SA_DIRK_4_2_3`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ARK324L2SA_DIRK_4_2_3"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the implicit portion of the default 3rd order additive method.  Both the
method and embedding are A-stable; additionally the method is L-stable
(this is the implicit portion of the ARK3(2)4L[2]SA method from
:cite:p:`KenCarp:03`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_ESDIRK325L2SA_5_2_3`` as the default 3rd order implicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     \frac{1767732205903}{2027836641118} & \frac{1767732205903}{4055673282236} & \frac{1767732205903}{4055673282236} & 0 & 0 \\
     \frac{3}{5} & \frac{2746238789719}{10658868560708} & -\frac{640167445237}{6845629431997} & \frac{1767732205903}{4055673282236} & 0 \\
     1 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     \hline
     3 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     2 & \frac{2756255671327}{12835298489170} & -\frac{10771552573575}{22201958757719} & \frac{9247589265047}{10645013368117} & \frac{2193209047091}{5459859503100}
   \end{array}

.. figure:: /figs/arkode/stab_region_15.png
   :scale: 50 %
   :align: center

   Linear stability region for the implicit ARK324L2SA-DIRK-4-2-3 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK436L2SA_6_3_4

Accessible via the constant ``ARKODE_ESDIRK436L2SA_6_3_4`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK436L2SA_6_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the default 4th order implicit method and the ESDIRK4(3)6L[2]SA method
from :cite:p:`KenCarp:16`. Both the method and embedding are A- and L-stable.

.. versionchanged:: 6.3.0

   Made the default 4th order implicit method

.. figure:: /figs/arkode/stab_region_28.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK436L2SA-6-3-4 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_CASH_5_2_4

Accessible via the constant ``ARKODE_CASH_5_2_4`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_CASH_5_2_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A-stable; additionally the method is L-stable
(from :cite:p:`Cash:79`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccc}
     0.435866521508 & 0.435866521508 & 0 & 0 & 0 & 0 \\
     -0.7 & -1.13586652150 & 0.435866521508 & 0 & 0 & 0 \\
     0.8 & 1.08543330679 & -0.721299828287 & 0.435866521508 & 0 & 0 \\
     0.924556761814 & 0.416349501547 & 0.190984004184 & -0.118643265417 & 0.435866521508 & 0 \\
     1 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     \hline
     4 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     2 & 1.05646216107052 & -0.0564621610705236 & 0 & 0 & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_16.png
   :scale: 50 %
   :align: center

   Linear stability region for the Cash-5-2-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_CASH_5_3_4

Accessible via the constant ``ARKODE_CASH_5_3_4`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_CASH_5_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A-stable; additionally the method is L-stable
(from :cite:p:`Cash:79`).

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccc}
     0.435866521508 & 0.435866521508 & 0 & 0 & 0 & 0 \\
     -0.7 & -1.13586652150 & 0.435866521508 & 0 & 0 & 0 \\
     0.8 & 1.08543330679 & -0.721299828287 & 0.435866521508 & 0 & 0 \\
     0.924556761814 & 0.416349501547 & 0.190984004184 & -0.118643265417 & 0.435866521508 & 0 \\
     1 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     \hline
     4 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     3 & 0.776691932910 & 0.0297472791484 & -0.0267440239074 & 0.220304811849 & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_17.png
   :scale: 50 %
   :align: center

   Linear stability region for the Cash-5-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_SDIRK_5_3_4

Accessible via the constant ``ARKODE_SDIRK_5_3_4`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_SDIRK_5_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Here, the method is both A- and L-stable, although the embedding has reduced
stability (from :cite:p:`HaWa:91`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_ESDIRK436L2SA_6_3_4`` as the default 4th order implicit method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccc}
     \frac{1}{4} & \frac{1}{4} & 0 & 0 & 0 & 0 \\
     \frac{3}{4} & \frac{1}{2} & \frac{1}{4} & 0 & 0 & 0 \\
     \frac{11}{20} & \frac{17}{50} & -\frac{1}{25} & \frac{1}{4} & 0 & 0 \\
     \frac{1}{2} & \frac{371}{1360} & -\frac{137}{2720} & \frac{15}{544} & \frac{1}{4} & 0 \\
     1 & \frac{25}{24} & -\frac{49}{48} & \frac{125}{16} & -\frac{85}{12} & \frac{1}{4} \\
     \hline
     4 & \frac{25}{24} & -\frac{49}{48} & \frac{125}{16} & -\frac{85}{12} & \frac{1}{4} \\
     3 & \frac{59}{48} & -\frac{17}{96} & \frac{225}{32} & -\frac{85}{12} & 0
   \end{array}

.. figure:: /figs/arkode/stab_region_18.png
   :scale: 50 %
   :align: center

   Linear stability region for the SDIRK-5-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_KVAERNO_5_3_4

Accessible via the constant ``ARKODE_KVAERNO_5_3_4`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_KVAERNO_5_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A-stable (from :cite:p:`Kva:04`).

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|ccccc}
        0 & 0 & 0 & 0 & 0 & 0 \\
        0.871733043 & 0.4358665215  & 0.4358665215  & 0 & 0 & 0 \\
        0.468238744853136 & 0.140737774731968 & -0.108365551378832 & 0.4358665215 & 0 & 0 \\
        1 & 0.102399400616089 & -0.376878452267324 & 0.838612530151233 & 0.4358665215 & 0 \\
        1 & 0.157024897860995 & 0.117330441357768 & 0.61667803039168 & -0.326899891110444 & 0.4358665215 \\
        \hline
        4 & 0.157024897860995 & 0.117330441357768 & 0.61667803039168 & -0.326899891110444 & 0.4358665215 \\
        3 & 0.102399400616089 & -0.376878452267324 & 0.838612530151233 & 0.4358665215 & 0
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/stab_region_19.png
   :scale: 50 %
   :align: center

   Linear stability region for the Kvaerno-5-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK436L2SA_DIRK_6_3_4

Accessible via the constant ``ARKODE_ARK436L2SA_DIRK_6_3_4`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ARK436L2SA_DIRK_6_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A-stable; additionally the method is L-stable
(this is the implicit portion of the ARK4(3)6L[2]SA method from
:cite:p:`KenCarp:03`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_ARK437L2SA_DIRK_7_3_4`` as the implicit portion of the default 4th order additive method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{4} & \frac{1}{4} & 0 & 0 & 0 & 0 \\
     \frac{83}{250} & \frac{8611}{62500} & -\frac{1743}{31250} & \frac{1}{4} & 0 & 0 & 0 \\
     \frac{31}{50} & \frac{5012029}{34652500} & -\frac{654441}{2922500} & \frac{174375}{388108} & \frac{1}{4} & 0 & 0 \\
     \frac{17}{20} & \frac{15267082809}{155376265600} & -\frac{71443401}{120774400} & \frac{730878875}{902184768} & \frac{2285395}{8070912} & \frac{1}{4} & 0 \\
     1 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac{1}{4} \\
     \hline
     4 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac{1}{4} \\
     3 & \frac{4586570599}{29645900160} & 0 & \frac{178811875}{945068544} & \frac{814220225}{1159782912} & -\frac{3700637}{11593932} & \frac{61727}{225920}
   \end{array}

.. figure:: /figs/arkode/stab_region_20.png
   :scale: 50 %
   :align: center

   Linear stability region for the ARK436L2SA-DIRK-6-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK437L2SA_DIRK_7_3_4

Accessible via the constant ``ARKODE_ARK437L2SA_DIRK_7_3_4`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ARK437L2SA_DIRK_7_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the implicit portion of the default 4th order additive method and the
implicit portion of the 4th order ARK4(3)7L[2]SA method from
:cite:p:`KenCarp:19`.  Both the method and embedding are A- and L-stable.

.. versionchanged:: 6.3.0

   Made the implicit portion of the default 4th order additive method

.. math::

   \renewcommand{\arraystretch}{1.5}
   \begin{array}{r|ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{247}{1000} & \frac{1235}{10000} & \frac{1235}{10000} & 0 & 0 & 0 & 0 & 0 \\
     \frac{4276536705230}{10142255878289} & \frac{624185399699}{4186980696204} & \frac{624185399699}{4186980696204} & \frac{1235}{10000} & 0 & 0 & 0 & 0 \\
     \frac{67}{200} & \frac{1258591069120}{10082082980243} & \frac{1258591069120}{10082082980243} & -\frac{322722984531}{8455138723562} & \frac{1235}{10000} & 0 & 0 & 0 \\
     \frac{3}{40} & -\frac{436103496990}{5971407786587} & -\frac{436103496990}{5971407786587} & -\frac{2689175662187}{11046760208243} & \frac{4431412449334}{12995360898505} & \frac{1235}{10000} & 0 & 0 \\
     \frac{7}{10} & -\frac{2207373168298}{14430576638973} & -\frac{2207373168298}{14430576638973} & \frac{242511121179}{3358618340039} & \frac{3145666661981}{7780404714551} & \frac{5882073923981}{14490790706663} & \frac{1235}{10000} & 0 \\
     1 & 0 & 0 & \frac{9164257142617}{17756377923965} & -\frac{10812980402763}{74029279521829} & \frac{1335994250573}{5691609445217} & \frac{2273837961795}{8368240463276} & \frac{1235}{10000} \\
     \hline
     4 & 0 & 0 & \frac{9164257142617}{17756377923965} & -\frac{10812980402763}{74029279521829} & \frac{1335994250573}{5691609445217} & \frac{2273837961795}{8368240463276} & \frac{1235}{10000} \\
     3 & 0 & 0 & \frac{4469248916618}{8635866897933} & -\frac{621260224600}{4094290005349} & \frac{696572312987}{2942599194819} & \frac{1532940081127}{5565293938103} & \frac{2441}{20000}
   \end{array}

.. figure:: /figs/arkode/stab_region_36.png
   :scale: 50 %
   :align: center

   Linear stability region for the ARK437L2SA-DIRK-7-3-4 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK43I6L2SA_6_3_4

Accessible via the constant ``ARKODE_ESDIRK43I6L2SA_6_3_4`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK43I6L2SA_6_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the ESDIRK4(3I)6L[2]SA method from :cite:p:`KenCarp:16`.
Both the method and embedding are A- and L-stable.

.. figure:: /figs/arkode/stab_region_29.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK43I6L2SA-6-3-4 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_QESDIRK436L2SA_6_3_4

Accessible via the constant ``ARKODE_QESDIRK436L2SA_6_3_4`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_QESDIRK436L2SA_6_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the QESDIRK4(3)6L[2]SA method from :cite:p:`KenCarp:16`.
Both the method and embedding are A- and L-stable.

.. figure:: /figs/arkode/stab_region_30.png
   :scale: 50 %
   :align: center

   Linear stability region for the QESDIRK436L2SA-6-3-4 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK437L2SA_7_3_4

Accessible via the constant ``ARKODE_ESDIRK437L2SA_7_3_4`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK437L2SA_7_3_4"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the ESDIRK4(3)7L[2]SA method from :cite:p:`KenCarp:19b`.
Both the method and embedding are A- and L-stable.

.. figure:: /figs/arkode/stab_region_31.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK437L2SA-7-3-4 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK547L2SA2_7_4_5

Accessible via the constant ``ARKODE_ESDIRK547L2SA2_7_4_5`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK547L2SA2_7_4_5"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the default 5th order implicit method and the ESDIRK5(4)7L[2]SA2 method
from :cite:p:`KenCarp:19b`. Both the method and embedding are A- and L-stable.

.. versionchanged:: 6.3.0

   Made the default 5th order implicit method

.. figure:: /figs/arkode/stab_region_33.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK547L2SA2-7-4-5 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_KVAERNO_7_4_5

Accessible via the constant ``ARKODE_KVAERNO_7_4_5`` to
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_KVAERNO_7_4_5"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A-stable; additionally the method is
L-stable (from :cite:p:`Kva:04`).

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|ccccccc}
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        0.52 & 0.26 & 0.26 & 0 & 0 & 0 & 0 & 0 \\
        1.230333209967908 & 0.13 & 0.84033320996790809 & 0.26 & 0 & 0 & 0 & 0 \\
        0.895765984350076 & 0.22371961478320505 & 0.47675532319799699 & -0.06470895363112615 & 0.26 & 0 & 0 & 0 \\
        0.436393609858648 & 0.16648564323248321 & 0.10450018841591720 & 0.03631482272098715 & -0.13090704451073998 & 0.26 & 0 & 0 \\
        1 & 0.13855640231268224 & 0 & -0.04245337201752043 & 0.02446657898003141 & 0.61943039072480676 & 0.26 & 0 \\
        1 & 0.13659751177640291 & 0 & -0.05496908796538376 & -0.04118626728321046 & 0.62993304899016403 & 0.06962479448202728 & 0.26 \\
        \hline
        5 & 0.13659751177640291 & 0 & -0.05496908796538376 & -0.04118626728321046 & 0.62993304899016403 & 0.06962479448202728 & 0.26 \\
        4 & 0.13855640231268224 & 0 & -0.04245337201752043 & 0.02446657898003141 & 0.61943039072480676 & 0.26 & 0
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/stab_region_21.png
   :scale: 50 %
   :align: center

   Linear stability region for the Kvaerno-7-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK548L2SA_DIRK_8_4_5

Accessible via the constant ``ARKODE_ARK548L2SA_DIRK_8_4_5`` for
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ARK548L2SA_DIRK_8_4_5"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
Both the method and embedding are A-stable; additionally the method is L-stable
(the implicit portion of the ARK5(4)8L[2]SA method from :cite:p:`KenCarp:03`).

.. versionchanged:: 6.3.0

   Replaced by ``ARKODE_ESDIRK547L2SA2_7_4_5`` as the default 5th order implicit method and replaced by ``ARKODE_ARK548L2SAb_DIRK_8_4_5`` as the implicit portion of the default 5th order additive method

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|cccccccc}
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{41}{100} & \frac{41}{200} & \frac{41}{200} & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{2935347310677}{11292855782101} & \frac{41}{400} & -\frac{567603406766}{11931857230679} & \frac{41}{200} & 0 & 0 & 0 & 0 & 0 \\
        \frac{1426016391358}{7196633302097} & \frac{683785636431}{9252920307686} & 0 & -\frac{110385047103}{1367015193373} & \frac{41}{200} & 0 & 0 & 0 & 0 \\
        \frac{92}{100} & \frac{3016520224154}{10081342136671} & 0 & \frac{30586259806659}{12414158314087} & -\frac{22760509404356}{11113319521817} & \frac{41}{200} & 0 & 0 & 0 \\
        \frac{24}{100} & \frac{218866479029}{1489978393911} & 0 & \frac{638256894668}{5436446318841} & -\frac{1179710474555}{5321154724896} & -\frac{60928119172}{8023461067671} & \frac{41}{200} & 0 & 0 \\
        \frac{3}{5} & \frac{1020004230633}{5715676835656} & 0 & \frac{25762820946817}{25263940353407} & -\frac{2161375909145}{9755907335909} & -\frac{211217309593}{5846859502534} & -\frac{4269925059573}{7827059040749} & \frac{41}{200} & 0 \\
        1 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
        \hline
        5 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
        4 & -\frac{975461918565}{9796059967033} & 0 & 0 & \frac{78070527104295}{32432590147079} & -\frac{548382580838}{3424219808633} & -\frac{33438840321285}{15594753105479} & \frac{3629800801594}{4656183773603} & \frac{4035322873751}{18575991585200}
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/stab_region_22.png
   :scale: 50 %
   :align: center

   Linear stability region for the implicit ARK548L2SA-ESDIRK-8-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ARK548L2SAb_DIRK_8_4_5

Accessible via the constant ``ARKODE_ARK548L2SAb_DIRK_8_4_5`` for
:c:func:`ARKStepSetTableNum` or
:c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ARK548L2SAb_DIRK_8_4_5"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the implicit portion of the default 5th order additive method.
Both the method and embedding are A-stable; additionally the method is L-stable
(this is the implicit portion of the 5th order ARK5(4)8L[2]SA method from
:cite:p:`KenCarp:19`).

.. versionchanged:: 6.3.0

   Made the implicit portion of the default 5th order additive method

.. only:: html

   .. math::

      \renewcommand{\arraystretch}{1.5}
      \begin{array}{r|cccccccc}
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{4}{9} & \frac{2}{9} & \frac{2}{9} & 0 & 0 & 0 & 0 & 0 & 0 \\
        \frac{6456083330201}{8509243623797} & \frac{2366667076620}{8822750406821} & \frac{2366667076620}{8822750406821} & \frac{2}{9} & 0 & 0 & 0 & 0 & 0 \\
        \frac{1632083962415}{14158861528103} & -\frac{257962897183}{4451812247028} & -\frac{257962897183}{4451812247028} & \frac{128530224461}{14379561246022} & \frac{2}{9} & 0 & 0 & 0 & 0 \\
        \frac{6365430648612}{17842476412687} & -\frac{486229321650}{11227943450093} & -\frac{486229321650}{11227943450093} & -\frac{225633144460}{6633558740617} & \frac{1741320951451}{6824444397158} & \frac{2}{9} & 0 & 0 & 0 \\
        \frac{18}{25} & \frac{621307788657}{4714163060173} & \frac{621307788657}{4714163060173} & -\frac{125196015625}{3866852212004} & \frac{940440206406}{7593089888465} & \frac{961109811699}{6734810228204} & \frac{2}{9} & 0 & 0 \\
        \frac{191}{200} & \frac{2036305566805}{6583108094622} & \frac{2036305566805}{6583108094622} & -\frac{3039402635899}{4450598839912} & -\frac{1829510709469}{31102090912115} & -\frac{286320471013}{6931253422520} & \frac{8651533662697}{9642993110008} & \frac{2}{9} & 0 \\
        1 & 0 & 0 & \frac{3517720773327}{20256071687669} & \frac{4569610470461}{17934693873752} & \frac{2819471173109}{11655438449929} & \frac{3296210113763}{10722700128969} & -\frac{1142099968913}{5710983926999} & \frac{2}{9} \\
        \hline
        5 & 0 & 0 & \frac{3517720773327}{20256071687669} & \frac{4569610470461}{17934693873752} & \frac{2819471173109}{11655438449929} & \frac{3296210113763}{10722700128969} & -\frac{1142099968913}{5710983926999} & \frac{2}{9} \\
        4 & 0 & 0 & \frac{520639020421}{8300446712847} & \frac{4550235134915}{17827758688493} & \frac{1482366381361}{6201654941325} & \frac{5551607622171}{13911031047899} & -\frac{5266607656330}{36788968843917} & \frac{1074053359553}{5740751784926}
      \end{array}

.. only:: latex

   The Butcher table is too large to fit in the PDF version of this documentation.  Please see the HTML documentation for the table coefficients.


.. figure:: /figs/arkode/stab_region_37.png
   :scale: 50 %
   :align: center

   Linear stability region for the ARK548L2SAb-DIRK-8-4-5 method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. c:enumerator:: ARKODE_ESDIRK547L2SA_7_4_5

Accessible via the constant ``ARKODE_ESDIRK547L2SA_7_4_5`` to
:c:func:`ARKStepSetTableNum` or :c:func:`ARKodeButcherTable_LoadDIRK`.
Accessible via the string ``"ARKODE_ESDIRK547L2SA_7_4_5"`` to
:c:func:`ARKStepSetTableName` or
:c:func:`ARKodeButcherTable_LoadDIRKByName`.
This is the ESDIRK5(4)7L[2]SA method from :cite:p:`KenCarp:16`.
Both the method and embedding are A- and L-stable.

.. figure:: /figs/arkode/stab_region_32.png
   :scale: 50 %
   :align: center

   Linear stability region for the ESDIRK547L2SA-7-4-5 method method.  The method's
   region is outlined in blue; the embedding's region is in red.


.. _Butcher.additive:

Additive Butcher tables
---------------------------

In the category of additive Runge--Kutta methods for split implicit and
explicit calculations, ARKODE includes methods that have orders 2
through 5, with embeddings that are of orders 1 through 4.  These
Butcher table pairs are as follows:

.. _ARKODE.Butcher.ARK_properties:
.. table:: Additive Butcher tables. The default method for each order is marked
   with an asterisk (*).

   +----------------------------------------------+-----------------------------------------------+--------+----------------+-------+
   | ERK Method ID                                | DIRK Method ID                                | Stages | Embedded Order | Order |
   +==============================================+===============================================+========+================+=======+
   | :c:enumerator:`ARKODE_ARK2_ERK_3_1_2`        | :c:enumerator:`ARKODE_ARK2_DIRK_3_1_2`        | 3      | 1              | 2*    |
   +----------------------------------------------+-----------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK324L2SA_ERK_4_2_3`  | :c:enumerator:`ARKODE_ARK324L2SA_DIRK_4_2_3`  | 4      | 2              | 3*    |
   +----------------------------------------------+-----------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK437L2SA_ERK_7_3_4`  | :c:enumerator:`ARKODE_ARK437L2SA_DIRK_7_3_4`  | 7      | 3              | 4*    |
   +----------------------------------------------+-----------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK436L2SA_ERK_6_3_4`  | :c:enumerator:`ARKODE_ARK436L2SA_DIRK_6_3_4`  | 6      | 3              | 4     |
   +----------------------------------------------+-----------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK548L2SAb_ERK_8_4_5` | :c:enumerator:`ARKODE_ARK548L2SAb_DIRK_8_4_5` | 8      | 4              | 5*    |
   +----------------------------------------------+-----------------------------------------------+--------+----------------+-------+
   | :c:enumerator:`ARKODE_ARK548L2SA_ERK_8_4_5`  | :c:enumerator:`ARKODE_ARK548L2SA_DIRK_8_4_5`  | 8      | 5              | 5     |
   +----------------------------------------------+-----------------------------------------------+--------+----------------+-------+



.. _Butcher.sprk:

Symplectic Partitioned Butcher tables
-------------------------------------

In the category of symplectic partitioned Runge-Kutta (SPRK) methods, ARKODE
includes methods that have orders :math:`q = \{1,2,3,4,5,6,8,10\}`. ARKODE's
symplectic partitioned Butcher tables are provided in the enumeration

.. c:enum:: ARKODE_SPRKMethodID

with values specified in :numref:`ARKODE.Butcher.SPRK_properties`.

.. _ARKODE.Butcher.SPRK_properties:
.. table:: Symplectic partitioned Butcher tables. The default method for each
   order is marked with an asterisk (*).

   +-------------------------------------------------+--------+-------+
   | Method ID                                       | Stages | Order |
   +=================================================+========+=======+
   | :c:enumerator:`ARKODE_SPRK_EULER_1_1`           | 1      | 1*    |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_LEAPFROG_2_2`        | 2      | 2*    |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_PSEUDO_LEAPFROG_2_2` | 2      | 2     |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_MCLACHLAN_2_2`       | 2      | 2     |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_MCLACHLAN_3_3`       | 3      | 3*    |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_RUTH_3_3`            | 3      | 3     |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_MCLACHLAN_4_4`       | 4      | 4*    |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_CANDY_ROZMUS_4_4`    | 4      | 4     |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_MCLACHLAN_5_6`       | 6      | 5*    |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_YOSHIDA_6_8`         | 8      | 6*    |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_SUZUKI_UMENO_8_16`   | 16     | 8*    |
   +-------------------------------------------------+--------+-------+
   | :c:enumerator:`ARKODE_SPRK_SOFRONIOU_10_36`     | 36     | 10*   |
   +-------------------------------------------------+--------+-------+

.. c:enumerator:: ARKODE_SPRK_EULER_1_1

Accessible via the constant (or string) ``ARKODE_SPRK_EULER_1_1`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the classic Symplectic Euler method and the default 1st order method.


.. c:enumerator:: ARKODE_SPRK_LEAPFROG_2_2

Accessible via the constant (or string) ``ARKODE_SPRK_LEAPFROG_2_2`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the classic Leapfrog/Verlet method and the default 2nd order method.


.. c:enumerator:: ARKODE_SPRK_PSEUDO_LEAPFROG_2_2

Accessible via the constant (or string) ``ARKODE_SPRK_PSEUDO_LEAPFROG_2_2`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the classic Pseudo Leapfrog/Verlet method.


.. c:enumerator:: ARKODE_SPRK_MCLACHLAN_2_2

Accessible via the constant (or string) ``ARKODE_SPRK_MCLACHLAN_2_2`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 2nd order method given by McLachlan in :cite:p:`Mclachlan:92`.


.. c:enumerator:: ARKODE_SPRK_MCLACHLAN_3_3

Accessible via the constant (or string) ``ARKODE_SPRK_MCLACHLAN_3_3`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 3rd order method given by McLachlan in :cite:p:`Mclachlan:92`
and the default 3rd order method.


.. c:enumerator:: ARKODE_SPRK_RUTH_3_3

Accessible via the constant (or string) ``ARKODE_SPRK_RUTH_3_3`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 3rd order method given by Ruth in :cite:p:`Ruth:93`.


.. c:enumerator:: ARKODE_SPRK_MCLACHLAN_4_4

Accessible via the constant (or string) ``ARKODE_SPRK_MCLACHLAN_4_4`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 4th order method given by McLachlan in :cite:p:`Mclachlan:92`
and the default 4th order method.

.. warning::

   This method only has coefficients sufficient for single or double precision.


.. c:enumerator:: ARKODE_SPRK_CANDY_ROZMUS_4_4

Accessible via the constant (or string) ``ARKODE_SPRK_CANDY_ROZMUS_4_4`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 4th order method given by Candy and Rozmus in :cite:p:`CandyRozmus:91`.


.. c:enumerator:: ARKODE_SPRK_MCLACHLAN_5_6

Accessible via the constant (or string) ``ARKODE_SPRK_MCLACHLAN_5_6`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 5th order method given by McLachlan in :cite:p:`Mclachlan:92`
and the default 5th order method.

.. warning::

   This method only has coefficients sufficient for single or double precision.


.. c:enumerator:: ARKODE_SPRK_YOSHIDA_6_8

Accessible via the constant (or string) ``ARKODE_SPRK_YOSHIDA_6_8`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 6th order method given by Yoshida in :cite:p:`Yoshida:90`
and the default 6th order method.


.. c:enumerator:: ARKODE_SPRK_SUZUKI_UMENO_8_16

Accessible via the constant (or string) ``ARKODE_SPRK_SUZUKI_UMENO_8_16`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 8th order method given by Suzuki and Umeno in :cite:p:`Suzuki:93`
and the default 8th order method.


.. c:enumerator:: ARKODE_SPRK_SOFRONIOU_10_36

Accessible via the constant (or string) ``ARKODE_SPRK_SOFRONIOU_10_36`` to
:c:func:`ARKodeSPRKTable_Load` or :c:func:`ARKodeSPRKTable_LoadByName`.
This is the 10th order method given by Sofroniou and Spaletta in :cite:p:`Sofroniou:05`
and the default 10th order method.
