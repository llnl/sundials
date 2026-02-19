..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
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
   -----------------------------------------------------------------------------

.. _CVODE.Alg.Nordsieck:

The Nordsieck History Array
---------------------------

A step from time :math:`t_{n-1}` to :math:`t_n` with step size :math:`h_n =
t_{n} - t_{n-1}` using a method of order :math:`q` begins with the predictor
polynomial interpolant, :math:`\pi_{p,n-1}`, of degree :math:`q` or less. For
Adams methods, the predictor interpolant satisfies the :math:`q + 1` conditions

.. math::
   :label: cvode_adams_predictor_polynomial

   \pi_{p,n-1}(t_{n-1}) &= y_{n-1} \\
   \dot{\pi}_{p,n-1}(t_{n-j}) &= \dot{y}_{n-j} = f_{n-j}, \quad j = 1,2,\ldots,q.

Similarly for BDF methods, the predictor interpolant satisfies

.. math::
   :label: cvode_bdf_predictor_polynomial

   \pi_{p,n-1}(t_{n-j}) &= y_{n-j}, \qquad\qquad j = 1,2,\ldots,q \\
   \dot{\pi}_{p,n-1}(t_{n-1}) &= \dot{y}_{n-1} = f_{n-1}.

The solution at :math:`t_n` is obtained by constructing the corrector polynomial
interpolant, :math:`\pi_{c,n}`, of degree :math:`q` or less. For Adams methods,
the corrector interpolant satisfies the :math:`q + 2` conditions

.. math::
   :label: cvode_adams_corrector_polynomial

   \pi_{c,n}(t_n) &= y_n, \\
   \pi_{c,n}(t_{n-1}) &= y_{n-1}, \\
   \dot{\pi}_{c,n}(t_{n-j}) &= \dot{y}_{n-j} = f_{n-j} \quad j = 0,1,\ldots,q-1.

Similarly for variable coefficient (VC) BDF methods, the corrector interpolant
satisfies

.. math::
   :label: cvode_vc_bdf_corrector_polynomial

   \pi_{c,n}(t_{n-j}) &= y_{n-j}, \qquad j = 0,1,\ldots,q \\
   \dot{\pi}_{c,n}(t_n) &= \dot{y}_n = f_{n}.

For fixed-leading-coefficient (FLC) BDF methods, the corrector interpolant
satisfies the conditions

.. math::
   :label: cvode_flc_bdf_corrector_polynomial

   \pi_{c,n}(t_{n}) &= y_{n}, \\
   \dot{\pi}_{c,n}(t_n) &= \dot{y}_n = f_{n}, \\
   \pi_{c,n}(t_{n} - j h_n) &= \pi_{p,n-1}(t_n - j h_n), \quad j = 1,\ldots,q

Thus, the FLC corrector interpolates the predictor at evenly spaced past points
rather than interpolating the computed values as with the VC corrector.

With either Adams or BDF methods the solution history is represented using the
Nordsieck array given by

.. math::
   :label: cvode_nordsieck_array

   Z_{n-1} = [ y_{n-1},\, h_n \dot{y}_{n-1},\, \frac{1}{2} h_n^2 \ddot{y}_{n-1},\, \ldots,\, \frac{1}{q!} h_n^q y^{(q)}_{n-1}]

where the derivatives, :math:`y^{(j)}_{n-1}` are approximated by
:math:`\pi^{(j)}_{p,n-1}(t_{n-1})`. With this formulation, computing the time step
is now a problem of constructing :math:`Z_n` from :math:`Z_{n-1}`. This
computation is done by forming the predicted array, :math:`Z_{n(0)}`, with
columns :math:`\frac{1}{j!} h^{(j)}_{n} \pi^{(j)}_{p,n-1}(t_n)` and computing the
coefficients :math:`l = [l_0,\, l_1,\, \ldots,\, l_q]` such that

.. math::
   :label: cvode_nordsieck_update

   Z_n = Z_{n(0)} + \Delta_n\, l,

where :math:`\Delta_n` is the correction to the predicted solution i.e.,

.. math::
   :label: cvode_correction

   \Delta_n = y_n - y_{n(0)} = \pi_{c,n}(t_n) - \pi_{p,n-1}(t_n)

From the second column of :eq:`cvode_nordsieck_update` we obtain the nonlinear
system for computing :math:`y_n`,

.. math::
   :label: cvode_nonlinear_system

   F(y_n) = y_n - \gamma f(t_n, y_n) - a_n = 0,

where :math:`\gamma = h_n / l_1` and :math:`a_n = y_{n(0)} - \gamma
\dot{y}_{n(0)}`.


Changing the State Size
^^^^^^^^^^^^^^^^^^^^^^^

After completing a step from time :math:`t_{n-1}` to :math:`t_n` using a method
of order :math:`q` and before starting the next step from :math:`t_{n}` to
:math:`t_{n+1}` with a method of order :math:`q'`, the size of the state vector
may change. To "resize" the integrator without restarting from first order, we
require that the user supply, depending on the method type, either the recent
solution or right-hand side history for the new state size.

Continuing the integration with the updated state size requires constructing a
new Nordsieck array, :math:`Z_n`, given the necessary history to evaluate the
appropriate predictor interpolating polynomial and its derivatives at
:math:`t_n`. To compute :math:`\pi_{p,n}(t_n)` and :math:`\pi^{(j)}_{p,n}(t_n)`,
we use a Newton interpolating polynomial. Given a set of :math:`k+1` data points
:math:`\{x_n,\ldots,x_{n-k}\}` to interpolate and a corresponding set of times,
:math:`\{t_n,\ldots,t_{n-k}\}`, the :math:`k`-th degree polynomial is given by

.. math::
   :label: newton_polynomial

   P(t) = \sum_{j=0}^{k} c_j N_j(t).

The polynomial coefficients, :math:`c_j`, are given by the divided differences,

.. math::
   :label: newton_polynomial_coefficients

   c_j = [x_n,\ldots,x_{n-j}],

where :math:`[x_n,\ldots,x_{n-j}]` is defined recursively as

.. math::
   :label: divided_differences

   [x_i,\ldots,x_{i-j}] = \frac{[x_{i},\ldots , x_{i-j+1}] - [x_{i-1},\ldots , x_{i-j}]}{t_i - t_{i-j}}.

with :math:`[x_i] = x_i`. The basis polynomials, :math:`N_j(t)`, are given by

.. math::
   :label: newton_polynomial_basis

   N_j(t) = \prod_{i=0}^{j-1} (t - t_{n-i}),

for :math:`j > 0` with :math:`N_0(t) = 1`.

Organizing the divided differences in a table illustrates the dependencies in
the recursive computation. For example with :math:`k = 3`, we need to compute
:math:`c_0 = [x_n]`, :math:`c_1 = [x_n,x_{n-1}]`, :math:`c_2 =
[x_n,x_{n-1},x_{n-2}]`, and :math:`c_3 = [x_n,x_{n-1},x_{n-2},x_{n-3}]` which
depend on the difference of the entries to the immediate lower and upper left in
the table below.

.. math::
   :label: divided_differences_table

   \begin{matrix}
   t_n     & x_n     & : & [x_n]     &      &                   &      &                           &      & \\
           &         & : &           & \rhd & [x_n,x_{n-1}]     &      &                           &      & \\
   t_{n-1} & x_{n-1} & : & [x_{n-1}] &      &                   & \rhd & [x_n,x_{n-1},x_{n-2}]     &      & \\
           &         & : &           & \rhd & [x_{n-1},x_{n-2}] &      &                           & \rhd & [x_n,x_{n-1},x_{n-2},x_{n-3}] \\
   t_{n-2} & x_{n-2} & : & [x_{n-2}] &      &                   & \rhd & [x_{n-1},x_{n-2},x_{n-3}] &      & \\
           &         & : &           & \rhd & [x_{n-2},x_{n-3}] &      &                           &      & \\
   t_{n-3} & x_{n-3} & : & [x_{n-3}] &      &                   &      &                           &      &
   \end{matrix}

As such, the coefficients can be computed recursively with the following steps:

1. For :math:`i = 0,\ldots,k`

   :math:`c_i = x_i`

2. For :math:`i = 1,\ldots,k`

   a. For :math:`j = k,\ldots,i`

      :math:`c_j = \frac{c_{j-1} - c_j}{t_{j - i} - t_{j}}`

Rewriting the Newton polynomial using nested multiplications,

.. math::
   :label: newton_nested

   P(t) = c_0 + (t - t_n) \Biggl[ c_1 + (t - t_{n-1}) \biggl[ c_2 + (t - t_{n-2}) \Bigl[ \ldots \bigl[ c_{k-1} + (t - t_{ n - k + 1}) c_k \bigr] \Bigr] \biggr] \Biggr],

leads the following iteration to evaluate :math:`P(t)`:

1. Let :math:`P_k(t) = c_k`

2. For :math:`j = k-1, \ldots, 0`

   :math:`P_{j}(t) = c_j + (t - t_{n-j}) P_{j+1}(t)`

Utilizing this recursive relationship for :math:`P(t)`, we can similarly compute
its derivatives. The evaluation of the :math:`d`-th derivative is given by the
following iteration:

1. Let :math:`P_k(t) = c_k` and :math:`\dot{P}_k(t) = \ddot{P}_k(t) = \ldots = P^{(d)}_k(t) = 0`

2. For :math:`j = k-1, \ldots, 0`

   :math:`P_{j}(t) = c_j + (t - t_{n-j}) P_{j+1}(t)`

   :math:`\dot{P}_{j}(t) = P_{j+1}(t) + (t - t_{n-j}) \dot{P}_{j+1}(t)`

   :math:`\ddot{P}_{j}(t) = 2\, \dot{P}_{j+1}(t) + (t - t_{n-j}) \ddot{P}_{j+1}(t)`

   :math:`\vdots`

   :math:`P^{(d)}_{j}(t) = d\, P^{(d-1)}_{j+1}(t) + (t - t_{n-j}) P^{(d)}_{j+1}(t)`

While only :math:`P^{(d)}_k(t) = 0` is needed to start the iteration, note that
:math:`P^{(d)}_{k} = P^{(d)}_{k - 1} = \ldots = P^{(d)}_{k - d + 1} = 0` so some
computations can be skipped when :math:`j > k - d`. With these pieces in place
we can perform the necessary evaluations to build :math:`Z_n`.

For Adams methods we require :math:`y_n` and :math:`q'` right-hand side values,
:math:`f_{n-j}` for :math:`j = 0,1,\ldots,q'-1`. The first two columns of
:math:`Z_n` are simply :math:`y_n` and :math:`h_n f_n`. To evaluate the higher
order derivatives to fill the remaining :math:`q' - 1` columns of :math:`Z_n` we
use :math:`P(t) = \dot{\pi}_{p,n}` instead of :math:`\pi_{p,n}` as the term
interpolating :math:`y_n` vanishes when taking the first derivative and we
already have the data needed for the first column of :math:`Z_n`. In this case,
:math:`P(t)` interpolates the :math:`q'` right-hand side values i.e.,
:math:`x_{n-j} = f_{n-j}` above for :math:`j = 0,1,\ldots,q'-1 = k`. For example
with :math:`q' = k + 1 = 3`, the table of divided differences is

.. math::
   :label: divided_differences_table_adams

   \begin{matrix}
   t_n     & f_n     & : & [f_n]     &      &                   &      & \\
           &         &   &           & \rhd & [f_n,f_{n-1}]     &      & \\
   t_{n-1} & f_{n-1} & : & [f_{n-1}] &      &                   & \rhd & [f_n,f_{n-1},y_{n-2}] \\
           &         &   &           & \rhd & [f_{n-1},f_{n-2}] &      & \\
   t_{n-2} & f_{n-2} & : & [f_{n-2}] &      &                   &      &
   \end{matrix}

With BDF methods we need :math:`f_n` and :math:`q'` solution values,
:math:`y_{n-j}` for :math:`j = 0,1,\ldots,q'-1`. Again, the first two columns of
:math:`Z_n` are simply :math:`y_n` and :math:`h_n f_n`. To evaluate the higher
order derivatives to fill the remaining :math:`q' - 1` columns of :math:`Z_n` we
need a slight modification of the Newton polynomial evaluation procedure
described above to build a Hermite polynomial that incorporates the derivative
interpolation condition, :math:`\dot{P}(t_n) = f_n`. In this case we duplicate
the data point at :math:`t_n` and replace the divided difference :math:`[y_n,
y_n]` with the corresponding derivative value, :math:`f_n`. For example with
:math:`q' = k = 3`, the table of divided differences is

.. math::
   :label: divided_differences_table_bdf

   \begin{matrix}
   t_n     & y_n     & : & [y_n]     &      &                   &      &                           &      & \\
           &         & : &           & \rhd & [y_n,y_n] = f_n   &      &                           &      & \\
   t_{n}   & y_n     & : & [y_n]     &      &                   & \rhd & [y_n,y_n,y_{n-1}]         &      & \\
           &         & : &           & \rhd & [y_n,y_{n-1}]     &      &                           & \rhd & [y_n,y_n,y_{n-1},y_{n-2}] \\
   t_{n-1} & y_{n-1} & : & [y_{n-1}] &      &                   & \rhd & [y_n,y_{n-1},y_{n-2}]     &      & \\
           &         & : &           & \rhd & [y_{n-1},y_{n-2}] &      &                           &      & \\
   t_{n-2} & y_{n-2} & : & [y_{n-2}] &      &                   &      &                           &      &
   \end{matrix}

Other than this adjustment in computing the :math:`c_1` polynomial coefficient,
the iterations to evaluate :math:`P(t_n)` and :math:`P^{(d)}(t_n)` to fill the
remaining columns of :math:`Z_n` remain the same.

.. note::

   With both Adams and BDF methods the second entry of :math:`Z_n` uses the
   input value of :math:`f_n`. This will differ from what occurs in a step
   without resizing where the entry corresponding to :math:`f_n` is obtained
   from the correction to the predicted Nordsieck array and not an evaluation of
   the right-hand side function at :math:`y_n`.

   Additionally, when the method order is increased in the next step, we use the
   provided history to directly construct :math:`Z_n` at the new order instead
   of adjusting the order using the correction vector.

Beyond building :math:`Z_n`, we also need to compute a resized correction vector
for the just completed step, :math:`\Delta_n`, in order to test for a potential
order increase in the next step. This calculation is achieved by computing a
resized prediction, :math:`y_{n(0)}`, for the just-completed step and
subtracting it from the resized solution, :math:`y_n`. Depending on the method
order for the next step, :math:`q'`, this computation will require one (if
:math:`q' = q+1`), two (if :math:`q' = q`), or three (if :math:`q' = q-1`)
additional data points beyond what is needed to construct :math:`Z_n`. For Adams
methods, computing the prediction requires resized versions of :math:`y_{n-1}`
(always), :math:`f_{n-q}` (if :math:`q' = q` or :math:`q' = q-1`), and
:math:`f_{n-q-1}` (if :math:`q' = q-1`). Similarly for BDF methods, computing
the prediction requires resized versions of :math:`f_{n-1}` (always),
:math:`y_{n-q}` (if :math:`q' = q` or :math:`q' = q-1`), and :math:`y_{n-q-1}`
(if :math:`q' = q-1`).
