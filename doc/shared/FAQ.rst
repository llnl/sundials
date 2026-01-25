..
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

.. _FAQ:

##########################
Frequently Asked Questions
##########################

The following are some questions frequently asked by SUNDIALS users.
If you do not see your question here, please do not hesitate to ask the
SUNDIALS mailing list by emailing your question to ``sundials-users@llnl.gov``,
or by opening an issue on our `GitHub <https://github.com/LLNL/sundials>`_.

Installation
------------

.. collapse:: I have the sundials source, how do I build it?

   See the :ref:`Installation.CMake` section.


.. collapse:: Can I use a non-default compiler to build SUNDIALS?

   Yes, specific compilers can be specified on the CMake command line.
   The following example specifies gcc and g++ for the C and CXX compilers, respectively:

   .. code-block:: cpp

      cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  ../srcdir

   .. note::
      If you wish to change the compiler after already configuring, then you must start
      with a fresh build directory and specify the compiler on the command line as shown.


.. collapse:: How do I install SUNDIALS on Windows systems?

   One way of obtaining Windows libraries for the SUNDIALS solvers is to use
   `cygwin <https://www.cygwin.com/>`__, in which case the installation
   procedure is the same as for any other Linux/Unix system (see the
   :ref:`Installation.CMake.Unix` installation section).

   Otherwise, refer to the :ref:`Installation.CMake.Windows` installation
   section.


.. collapse:: Everything installed fine! How do I link the SUNDIALS libraries to my own application?

   Refer to :ref:`Installation.UsingSUNDIALS`.


CVODE(S) / IDA(S) / ARKODE
--------------------------

.. collapse:: How do I choose tolerances?

   The same advice applies to CVODE(S), IDA(S), and ARKODE:

   1. The scalar relative tolerance ``reltol`` is to be set to control relative errors. So
   :math:`\texttt{reltol} = 10^{-4}` means that errors are controlled to .01%. We do not recommend
   using ``reltol`` larger than :math:`10^{-3}`. On the other hand, ``reltol`` should not be so small
   that it is comparable to the unit roundoff of the machine arithmetic (generally
   around :math:`10^{-15}` for double-precision).

   2. The absolute tolerances ``abstol`` (whether scalar or vector) need to be set to control
   absolute errors when any components of the solution vector ``y`` may be so small that
   pure relative error control is meaningless. For example, if ``y[i]`` starts at some
   nonzero value, but in time decays to zero, then pure relative error control on ``y[i]``
   makes no sense (and is overly costly) after ``y[i]`` is below some noise level. Then
   ``abstol`` (if scalar) or ``abstol[i]`` (if a vector) needs to be set to that noise level. If the different
   components have different noise levels, then ``abstol`` should be a vector. See the example  ``cvsRoberts_dns``
   in the CVODE package, and the discussion of it in the CVODE Examples document
   :cite:p:`cvodes_ex`. In that problem, the three components vary between 0 and 1,
   and have different noise levels; hence the ``abstol`` vector. It is impossible to give any
   general advice on ``abstol`` values, because the appropriate noise levels are completely
   problem-dependent. The user or modeler hopefully has some idea as to what those
   noise levels are.

   3. Finally, it is important to pick all the tolerance values conservatively,
   because they control the error committed on each individual time step. The final
   (global) errors are some sort of accumulation of those per-step errors. A good
   rule of thumb is to reduce the tolerances by a factor of .01 from the actual
   desired limits on errors. So if you want .01% accuracy (globally), a good choice
   is :math:`\texttt{reltol} = 10^{-6}`. But in any case, it is a good idea to do a few experiments
   with the tolerances to see how the computed solution values vary as tolerances
   are reduced.


.. collapse:: How do I choose what linear solver to use for the stiff case?

   If the problem size is fairly small (say :math:`N < 100`), then using the dense solver is
   probably best; it is the simplest to use, and reasonably inexpensive for small :math:`N`. For larger :math:`N`, it
   is important to take advantage of sparsity (zero-nonzero) structure within the problem. If there
   is local (nearest-neighbor) coupling, or if the coupling is local after a suitable reordering of
   :math:`y`, then use the banded linear solver. Local coupling means that the :math:`i`-th component of the RHS or
   residual function depends only on components :math:`y_j` for which :math:`|i-j|` is small relative
   to :math:`N`. (Note that the dense and band solvers are only applicable for the single node versions of the
   solver.) For even larger problems, consider one of the Krylov iterative methods. These are hardest
   to use, because for best results they usually require preconditioning. However, they offer the best
   opportunity to exploit the sparsity structure in the problem. The preconditioner is a matrix
   which, at least crudely, approximates the actual matrix in the linear system to be solved, and is
   typically built from an approximation of the relevant Jacobian matrix. Typically, that
   approximation uses only part of the true Jacobian, but as a result is much less expensive to
   solve. If the Jacobian can be approximated by a matrix that is banded (serial case) or
   block-diagonal with banded blocks (distributed parallel case), SUNDIALS includes preconditioner modules for
   such cases. In each of the user guides, the section 'Linear solver specification functions' and
   the section on preconditioner modules contain more detailed comments on preconditioning. On the
   construction of preconditioners for problems arising from the spatial discretization of
   time-dependent partial differential equation systems, there is considerable discussion in the
   paper :cite:p:`BrHi:89`.

.. collapse:: How do I handle a data-defined function within the RHS or residual function?

   Often the RHS or residual function depends on some function :math:`A(t)` that is data-defined,
   i.e. defined only at a set of discrete set of times :math:`t`. The solver must be able to obtain values of
   the user-supplied functions at arbitrary times :math:`t` in the integration interval. So the user must fit
   the data with a reasonably smooth function :math:`A(t)` that is defined continuously for all
   relevant :math:`t`, and incorporate an evaluation of that fit function in the user function involved. This
   may be as simple as a piecewise linear fit, but a smoother fit (e.g. spline) would make the
   integration more efficient. If there is noise in the data, the fit should be a least-squares fit
   instead of a straight interpolation. The same advice applies if the user function has a
   data-defined function :math:`A(y)` that involves one or more components of the dependent variable
   vector :math:`y`. Of course, if more that one component is involved, the fit is more complicated.

.. collapse:: How do I control unphysical negative values?

   In many applications, some components in the true solution are always positive
   or non-negative, though at times very small. In the numerical solution, however,
   small negative (hence unphysical) values can then occur. In most cases, these
   values are harmless, and simply need to be controlled, not eliminated. The
   following pieces of advice are relevant.

   1. The way to control the size of unwanted negative computed values is with
   tighter absolute tolerances. Again this requires some knowledge of the noise
   level of these components, which may or may not be different for different
   components. Some experimentation may be needed.

   2. If output plots or tables are being generated, and it is important to avoid
   having negative numbers appear there (for the sake of avoiding a long
   explanation of them, if nothing else), then eliminate them, but only in the
   context of the output medium. Then the internal values carried by the solver are
   unaffected. Remember that a small negative value in ``y`` returned by CVODE, with
   magnitude comparable to ``abstol`` or less, is equivalent to zero as far as the computation
   is concerned.

   3. The user’s right-hand side routine ``f`` (or residual ``F``) should never change a negative value in
   the solution vector ``y`` to a non-negative value, as a "solution" to this problem.
   This can cause instability. If the ``f`` (or ``F``) routine cannot tolerate a zero or negative
   value (e.g. because there is a square root or log of it), then the offending
   value should be changed to zero or a tiny positive number in a temporary
   variable (not in the input ``y`` vector) for the purposes of computing :math:`f(t,y)` (or :math:`F(t,y,y')`).

   4. Positivity and non-negativity constraints on components can be enforced by
   use of the recoverable error return feature in the user-supplied right-hand side
   function. However, because this option involves some extra overhead cost, it
   should only be exercised if the use of absolute tolerances to control the
   computed values is unsuccessful.

   In addition, SUNDIALS integrators provide the option of enforcing positivity or non-negativity on components. But
   these constraint options should only be exercised if the use of absolute tolerances to control the
   computed values is unsuccessful, because they involve some extra overhead cost.


.. collapse:: How do I treat discontinuities in the RHS or residual function?

   If the jumps at the discontinuities are relatively small, simply keep them in the RHS (or residual) function,
   and let the integrator respond to them (possibly taking smaller steps through each point of
   discontinuity). If the jumps are large, it is more efficient to stop at the point of discontinuity
   and restart the integrator with a readjusted ODE (or DAE) model. To stop when the location of the
   discontinuity is known, simply make that location a value of ``tout``. To stop when the location of
   the discontinuity is determined by the solution, use the rootfinding feature. In either case, it
   is critical that the RHS (or residual) function not incorporate the discontinuity, but rather have a smooth
   extension over the discontinuity, so that the step across it (and subsequent rootfinding, if used)
   can be done efficiently. Then use a switch within the RHS (or residual) function that can be flipped between the
   stopping of the integration and the restart, so that the restarted problem uses the new values
   (which have jumped).


.. collapse:: When is it advantageous to supply my own error weight function?

   The main situation where supplying an ``EwtFn`` function is a good idea is where the problem needs something "in between" the
   cases covered by scalar and vector absolute tolerances. Namely, suppose there are a few groups of
   variables (relative to the total number of variables) such that all the variables in each group
   require the same value of ``abstol``, but these values are very different from one group to another.
   Then a user ``EwtFn`` function can keep an array of those values and construct the ``ewt`` vector without
   any additional storage. Also, in rare cases, one may want to use this option to apply different
   values of ``reltol`` to different variables (or groups of variables).


.. collapse:: How do switch on/off forward sensitivity computations in CVODES?

   If you want to turn on and off forward sensitivity calculations during several successive
   integrations (such as if you were using CVODES within a dynamically-constrained optimization loop,
   when sometimes you want to only integrate the states and sometimes you also need sensitivities
   computed), it is most efficient to use :c:func:`CVodeSensToggleOff`.


.. collapse:: What is the role of plist in CVODES?

   The argument ``plist`` to :c:func:`CVodeSetSensParams` is used to specify the problem parameters with
   respect to which solution sensitivities are to be computed.

   ``plist`` is used only if the sensitivity right-hand sides are evaluated using the internal
   difference-quotient approximation function. In that case, ``plist`` should be declared as an array of
   ``Ns`` integers and should contain the indices in the array of problem parameters ``p`` with respect to
   which sensitivities are desired. For example, if you want to compute sensitivities with respect to
   the first and third parameters in the ``p`` array, ``p[0]`` and ``p[2]``, you need to set

   .. code-block:: C

      plist[0] = 0
      plist[1] = 2


   If ``plist`` is not provided, CVODES will compute sensitivities with respect to the first ``Ns``
   parameters in the array ``p`` (i.e. it will use ``plist[i]=i, i=0,1,...Ns``). If the user provides a
   function to evaluate the right-hand sides of the sensitivity equations or if the default values
   are desired, a ``NULL`` pointer can be passed to :c:func:`CVodeSetSensParams`.


.. collapse:: What is the role of pbar in CVODES?

   The argument ``pbar`` to :c:func:`CVodeSetSensParams` is used to specify scaling factors for the
   problem parameters.

   ``pbar`` is used only if

   * the internal difference-quotient functions are used for the evaluation of the sensitivity
      right-hand sides, in which case ``pbar`` is used in computing an appropriate perturbation for
      the finite-difference approximation

   or

   * the tolerances for the sensitivity variables are estimated automatically by CVODES from those
      specified for the state variables.

   If provided, ``pbar`` should be declared as an array of ``Ns`` real types and should contain non-zero
   scaling factors for the ``Ns`` parameters with respect to which sensitivities are to be computed. For
   non-zero problem parameters, a good choice is

   .. code-block:: C

      pbar[i] = p[plist[i]]


   If ``pbar`` is not provided, CVODES will use ``pbar[i]=1.0, i=0,1,...Ns-1``.

   If the user provides a function to evaluate the right-hand sides of the sensitivity equations and
   also specifies tolerances for the sensitivity variables (through the ``CVodeSens*tolerances``
   functions) or if the default values are desired, a ``NULL`` pointer can be passed to
   :c:func:`CVodeSetSensParams`.


.. collapse:: What is pure quadrature integration?

   Suppose your ODE is :math:`y'=f(t,y)` and you integrate it from :math:`0` to :math:`T` and that you are also interested in computing an integral of the form

   .. math::

      z(t) = \int_0^t g(t,y(t)) dt

   for some function :math:`g`. The most efficient way of computing :math:`z` is by appending one additional differential equation to your ODE system:

   .. math::

      z' = g(t,y)

   with initial condition :math:`z(0)=0`, in which case the integral from :math:`0` to :math:`T` is `z(T)`.

   This additional equation is "a pure quadrature equation" and its main characteristic is that the
   new differential variable :math:`z` does not appear in the right hand side of the extended ODE system. If
   CVODES is notified of such "pure quadrature equations", it can take advantage of this property and
   do less work than if it didn't know about them (these variables need not be considered in the
   nonlinear system solution).

   The main reason for the special treatment of "pure quadrature equations" in CVODES is that such
   integrals (very often a large number of them) need to be computed for adjoint sensitivity.


.. collapse:: When should I select a non-default temporal adaptivity controller in ARKODE?

   The default temporal adaptivity controller in ARKODE was selected due to its robust performance on test
   problems that ranged in difficulty and stiffness, and when running with a wide range of solution tolerances
   and method orders.  While we hope that this default runs well on most applications, it is unlikely to be
   optimal.  A prime indicator that an alternate adaptivity controller may be useful is if the default results
   in a large number of rejected steps.  Alternately, for higher-cost calculations where a reduction of 10%-20%
   in the number of time steps would be important, users may want to try another controller option.

   The default temporal adaptivity controller in ARKODE is the industry-standard *I* controller:

   .. math::
      h' = h_n \varepsilon_n^{-1/(p+1)}

   where :math:`\varepsilon_n = \text{bias}*\text{dsm}_n`, and :math:`\text{dsm}_n` is the WRMS norm of the
   local temporal error when using the time step :math:`h_n`, weighted by the user-requested tolerances.
   However, the :ref:`SUNAdaptController class <SUNAdaptController>` in SUNDIALS provides a range of more
   advanced temporal error controllers that could be applied to a given application problem, including the
   :math:`H_{0}211`, :math:`H_{0}321`, :math:`H211`, and :math:`H312` controllers from :cite:p:`Sod:03`,
   as well as a variety of controllers (*PI*, *PID*, *ExpGus*, *ImpGus*, and *ImExGus*) that were included
   in the initial ARKODE release.

   The adaptive time-stepping controllers introduced by Soderlind in :cite:p:`Sod:03` can be classified into
   two groups; *deadbeat* controllers and *non-deadbeat* controllers. A controller is known as deadbeat if
   the roots of its characteristic equation are located at the origin. These controllers are generalized
   forms of the *I* controller, are denoted with a zero subscript as part of the controller name (e.g.,
   :math:`H_{0}321`), and are generally recommended for applications with smooth solutions as a function of time.

   While it is impossible to exhaustively explore the question of controller optimality for all application
   problems, we have performed tests on the range of provided controllers on a variety of stiff and non-stiff
   problems, at varying tolerances (:math:`10^{-9} \to 10^{-3}`), and for Runge--Kutta methods at a wide range
   of orders of accuracy (ERK orders 2-9 and DIRK orders 2-5). From our experiments, stiff problems benefitted
   from the *I*, *PI*, :math:`H_{0}211`, :math:`H_{0}321`, and *ImpGus* controllers.  On the other hand,
   non-stiff test problems ran most efficiently when using the *PID*, *I*, *ExpGus*, :math:`H211`, and
   :math:`H312` controllers.

   Lastly, we note that the internal controller parameters for the legacy ARKODE controllers (*PI*, *PID*,
   *ExpGus*, *ImpGus*, and *ImExGus*) were determined via numerical optimization over a given set of test
   problems.  As a result, although those controllers work very well for some applications, they may not
   work well with others.  For users who are interested in exploring novel controller methods, we point out
   that the :ref:`Soderlind SUNAdaptController <SUNAdaptController.Soderlind>` allows complete control
   over all internal parameters via the :c:func:`SUNAdaptController_SetParams_Soderlind` function.


KINSOL
------

.. collapse:: How do I reinitialize KINSOL within a C/C++ program?

   Although KINSOL does not provide a reinitialization function, it is possible to reinitialize the
   solver (meaning reuse a KINSOL object), but only if the problem size remains unchanged. To
   reinitialize KINSOL, begin by making any necessary changes to the problem definition by calling
   the appropriate ``KINSet*`` functions (e.g., :c:func:`KINSetSysFunc`). Next, if you would like to use
   a different linear solver, call the appropriate function, followed by any calls to the
   corresponding ``KIN*Set*`` functions. Then you can call the ``KINSol`` function to solve the updated
   nonlinear algebraic system.


.. collapse:: Why is the system function being evaluated at points that violate the constraints?

   If you have not supplied a function to compute either :math:`J(u)` (of type :c:type:`KINLsJacFn`) or :math:`J(u) v`
   (of type :c:type:`KINLsJacTimesVecFn`), then the internal function may be the culprit. The
   default function used to compute a difference quotient approximation to the Jacobian (direct
   methods) or Jacobian matrix-vector product (Kylov methods) evaluates the user-supplied system
   function at a slightly perturbed point, but does not check if that point violates the constraints.


Miscellaneous
-------------

.. collapse:: How do I determine which version of SUNDIALS I have?

   If you still have access to the distribution files, then the SUNDIALS release number is indicated
   in the top-level ``README.md`` and the corresponding solver versions can be determined by
   reading the appropriate row of the :ref:`release history <ReleaseHistory>` table or from the files, ``sundials/src/<solver>/README.md``. You can also call the functions
   :c:func:`SUNDIALSGetVersion` and :c:func:`SUNDIALSGetVersionNumber` from your program, or
   use the ``SUNDIALS_VERSION*`` macros found in the header file ``sundials/sundials_config.h``.



.. collapse:: SUNDIALS Wiki

   Some additional information might be found at `http://sundials.wikidot.com
   <http://sundials.wikidot.com/>`_ however the wikidot page has not been maintained in many years so
   it contains plenty of outdated information.

   .. warning::

      The SUNDIALS team does not maintain the wikidot web page.
