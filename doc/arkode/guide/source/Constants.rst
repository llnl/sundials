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

.. _ARKODE.Constants:

================
ARKODE Constants
================

Below we list all input and output constants used by the main solver,
timestepper, and linear solver modules, together with a short
description of their meaning.  :numref:`ARKODE.Constants.in_constants`
contains the ARKODE input constants, and :numref:`ARKODE.Constants.out_constants`
contains the ARKODE output constants.

.. _ARKODE.Constants.in_constants:
.. table:: ARKODE input constants
   :widths: 38 52

   +-----------------------------------------------+------------------------------------------------------------+
   | **Shared input constants**                    |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_NORMAL`                           | Solver should return at a specified output time.           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_ONE_STEP`                         | Solver should return after each successful step.           |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Full right-hand side evaluation constants** |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_FULLRHS_START`                    | Calling the full right-hand side function at the           |
   |                                               | start of the integration.                                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_FULLRHS_END`                      | Calling the full right-hand side function at the end of    |
   |                                               | a step.                                                    |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_FULLRHS_OTHER`                    | Calling the full right-hand side function at the some      |
   |                                               | other point e.g., for dense output.                        |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Interpolation module input constants**      |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_NONE`                      | Disables polynomial interpolation for dense output.        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_HERMITE`                   | Specifies use of the Hermite polynomial interpolation      |
   |                                               | module (for non-stiff problems).                           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_LAGRANGE`                  | Specifies use of the Lagrange polynomial interpolation     |
   |                                               | module (for stiff problems).                               |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_INTERP_MAX_DEGREE`                | Maximum possible interpolating polynomial degree.          |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Relaxation module input constants**         |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_RELAX_BRENT`                      | Specifies Brent's method as the relaxation nonlinear       |
   |                                               | solver.                                                    |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARK_RELAX_NEWTON`                     | Specifies Newton's method as the relaxation nonlinear      |
   |                                               | solver.                                                    |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Default explicit Butcher tables**           |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_1`                | Use ARKStep's default first-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_FORWARD_EULER_1_1`.                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_2`                | Use ARKStep's default second-order ERK method              |
   |                                               | :c:enumerator:`ARKODE_RALSTON_3_1_2`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_3`                | Use ARKStep's default third-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_BOGACKI_SHAMPINE_4_2_3`.             |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_4`                | Use ARKStep's default fourth-order ERK method              |
   |                                               | :c:enumerator:`ARKODE_SOFRONIOU_SPALETTA_5_3_4`.           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_5`                | Use ARKStep's default fifth-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_TSITOURAS_7_4_5`.                    |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_6`                | Use ARKStep's default sixth-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_VERNER_9_5_6`.                       |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_7`                | Use ARKStep's default seventh-order ERK method             |
   |                                               | :c:enumerator:`ARKODE_VERNER_10_6_7`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_8`                | Use ARKStep's default eighth-order ERK method              |
   |                                               | :c:enumerator:`ARKODE_VERNER_13_7_8`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_9`                | Use ARKStep's default ninth-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_VERNER_16_8_9`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_1`                    | Use ERKStep's default first-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_FORWARD_EULER_1_1`.                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_2`                    | Use ERKStep's default second-order ERK method              |
   |                                               | :c:enumerator:`ARKODE_RALSTON_3_1_2`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_3`                    | Use ERKStep's default third-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_BOGACKI_SHAMPINE_4_2_3`.             |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_4`                    | Use ERKStep's default fourth-order ERK method              |
   |                                               | :c:enumerator:`ARKODE_SOFRONIOU_SPALETTA_5_3_4`.           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_5`                    | Use ERKStep's default fifth-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_TSITOURAS_7_4_5`.                    |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_6`                    | Use ERKStep's default sixth-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_VERNER_9_5_6`.                       |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_7`                    | Use ERKStep's default seventh-order ERK method             |
   |                                               | :c:enumerator:`ARKODE_VERNER_10_6_7`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_8`                    | Use ERKStep's default eighth-order ERK method              |
   |                                               | :c:enumerator:`ARKODE_VERNER_13_7_8`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_9`                    | Use ERKStep's default ninth-order ERK method               |
   |                                               | :c:enumerator:`ARKODE_VERNER_16_8_9`.                      |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Default implicit Butcher tables**           |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_1`               | Use ARKStep's default first-order DIRK method              |
   |                                               | :c:enumerator:`ARKODE_BACKWARD_EULER_1_1`.                 |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_2`               | Use ARKStep's default second-order DIRK method             |
   |                                               | :c:enumerator:`ARKODE_ARK2_DIRK_3_1_2`.                    |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_3`               | Use ARKStep's default third-order DIRK method              |
   |                                               | :c:enumerator:`ARKODE_ESDIRK325L2SA_5_2_3`.                |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_4`               | Use ARKStep's default fourth-order DIRK method             |
   |                                               | :c:enumerator:`ARKODE_ESDIRK436L2SA_6_3_4`.                |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_5`               | Use ARKStep's default fifth-order DIRK method              |
   |                                               | :c:enumerator:`ARKODE_ESDIRK547L2SA2_7_4_5`.               |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Default ImEx Butcher tables**               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ARK_ETABLE_2` &       | Use ARKStep's default second-order ARK method              |
   | :index:`ARKSTEP_DEFAULT_ARK_ITABLE_2`         | (ARKODE_ARK2_ERK_3_1_2 and ARKODE_ARK2_DIRK_3_1_2).        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ARK_ETABLE_3` &       | Use ARKStep's default third-order ARK method               |
   | :index:`ARKSTEP_DEFAULT_ARK_ITABLE_3`         | (ARKODE_ARK324L2SA_ERK_4_2_3 and                           |
   |                                               | ARKODE_ARK324L2SA_DIRK_4_2_3).                             |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ARK_ETABLE_4` &       | Use ARKStep's default fourth-order ARK method              |
   | :index:`ARKSTEP_DEFAULT_ARK_ITABLE_4`         | (ARKODE_ARK436L2SA_ERK_6_3_4 and                           |
   |                                               | ARKODE_ARK436L2SA_DIRK_6_3_4).                             |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ARK_ETABLE_5` &       | Use ARKStep's default fifth-order ARK method               |
   | :index:`ARKSTEP_DEFAULT_ARK_ITABLE_5`         | (ARKODE_ARK548L2SA_ERK_8_4_5 and                           |
   |                                               | ARKODE_ARK548L2SA_DIRK_8_4_5).                             |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **LSRK method types**                         |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKODE_LSRK_RKC_2`                    | 2nd order Runge-Kutta-Chebyshev (RKC) method               |
   |                                               | :c:enumerator:`ARKODE_LSRK_RKC_2`                          |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKODE_LSRK_RKL_2`                    | 2nd order Runge-Kutta-Legendre (RKL) method                |
   |                                               | :c:enumerator:`ARKODE_LSRK_RKL_2`                          |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKODE_LSRK_SSP_S_2`                  | Optimal 2nd order s-stage SSP RK method                    |
   |                                               | :c:enumerator:`ARKODE_LSRK_SSP_S_2`                        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKODE_LSRK_SSP_S_3`                  | Optimal 3rd order s-stage SSP RK method                    |
   |                                               | :c:enumerator:`ARKODE_LSRK_SSP_S_3`                        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`ARKODE_LSRK_SSP_10_4`                 | Optimal 4th order 10-stage SSP RK method                   |
   |                                               | :c:enumerator:`ARKODE_LSRK_SSP_10_4`                       |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **MRI method types**                          |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_EXPLICIT`                     | Use an explicit (at the slow time scale) MRI method.       |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_IMPLICIT`                     | Use an implicit (at the slow time scale) MRI method.       |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_IMEX`                         | Use an ImEx (at the slow time scale) MRI method.           |
   +-----------------------------------------------+------------------------------------------------------------+
   |                                               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | **Default MRI coupling tables**               |                                                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_1`               | Use MRIStep's default 1st-order explicit method            |
   |                                               | (ARKODE_MRI_GARK_FORWARD_EULER).                           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_2`               | Use MRIStep's default 2nd-order explicit method            |
   |                                               | (ARKODE_MRI_GARK_ERK22b).                                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_3`               | Use MRIStep's default 3rd-order explicit method            |
   |                                               | (ARKODE_MIS_KW3).                                          |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_4`               | Use MRIStep's default 4th-order explicit method            |
   |                                               | (ARKODE_MRI_GARK_ERK45a).                                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_2_AD`            | Use MRIStep's default 2nd-order adaptive explicit method   |
   |                                               | (ARKODE_MRI_GARK_ERK22a).                                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_3_AD`            | Use MRIStep's default 3rd-order adaptive explicit method   |
   |                                               | (ARKODE_MRI_GARK_ERK33a).                                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_4_AD`            | Use MRIStep's default 4th-order adaptive explicit method   |
   |                                               | (ARKODE_MRI_GARK_ERK45a).                                  |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_5_AD`            | Use MRIStep's default 5th-order adaptive explicit method   |
   |                                               | (ARKODE_MERK54).                                           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMPL_SD_1`            | Use MRIStep's default 1st-order solve-decoupled implicit   |
   |                                               | method (ARKODE_MRI_GARK_BACKWARD_EULER).                   |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMPL_SD_2`            | Use MRIStep's default 2nd-order solve-decoupled implicit   |
   |                                               | method (ARKODE_MRI_GARK_IRK21a).                           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMPL_SD_3`            | Use MRIStep's default 3rd-order solve-decoupled implicit   |
   |                                               | method (ARKODE_MRI_GARK_ESDIRK34a).                        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMPL_SD_4`            | Use MRIStep's default 4th-order solve-decoupled implicit   |
   |                                               | method (ARKODE_MRI_GARK_ESDIRK46a).                        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_1`            | Use MRIStep's default 1st-order solve-decoupled ImEx       |
   |                                               | method (ARKODE_IMEX_MRI_GARK_EULER).                       |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_2`            | Use MRIStep's default 2nd-order solve-decoupled ImEx       |
   |                                               | method (ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL).                 |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_3`            | Use MRIStep's default 3rd-order solve-decoupled ImEx       |
   |                                               | method (ARKODE_IMEX_MRI_GARK3b).                           |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_4`            | Use MRIStep's default 4th-order solve-decoupled ImEx       |
   |                                               | method (ARKODE_IMEX_MRI_GARK4).                            |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_2_AD`         | Use MRIStep's default 2nd-order solve-decoupled adaptive   |
   |                                               | ImEx method (ARKODE_IMEX_MRI_SR21).                        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_3_AD`         | Use MRIStep's default 3rd-order solve-decoupled adaptive   |
   |                                               | ImEx method (ARKODE_IMEX_MRI_SR32).                        |
   +-----------------------------------------------+------------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_4_AD`         | Use MRIStep's default 4th-order solve-decoupled adaptive   |
   |                                               | ImEx method (ARKODE_IMEX_MRI_SR43).                        |
   +-----------------------------------------------+------------------------------------------------------------+



.. _ARKODE.Constants.out_constants:
.. table:: ARKODE output constants
   :widths: 25 5 60

   +-------------------------------------+------+------------------------------------------------------------+
   | **Shared output constants**                                                                             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_SUCCESS`                | 0    | Successful function return.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TSTOP_RETURN`           | 1    | ARKODE succeeded by reaching the specified stopping point. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ROOT_RETURN`            | 2    | ARKODE succeeded and found one more more roots.            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_WARNING`                | 99   | ARKODE succeeded but an unusual situation occurred.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_MUCH_WORK`          | -1   | The solver took ``mxstep`` internal steps but could not    |
   |                                     |      | reach ``tout``.                                            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_MUCH_ACC`           | -2   | The solver could not satisfy the accuracy                  |
   |                                     |      | demanded by the user for some internal step.               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ERR_FAILURE`            | -3   | Error test failures occurred too many times during one     |
   |                                     |      | internal time step, or the minimum step size was reached.  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONV_FAILURE`           | -4   | Convergence test failures occurred too many times during   |
   |                                     |      | one internal time step, or the minimum step size was       |
   |                                     |      | reached.                                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LINIT_FAIL`             | -5   | The linear solver's initialization function failed.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LSETUP_FAIL`            | -6   | The linear solver's setup function failed in an            |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LSOLVE_FAIL`            | -7   | The linear solver's solve function failed in an            |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RHSFUNC_FAIL`           | -8   | The right-hand side function failed in an                  |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_FIRST_RHSFUNC_ERR`      | -9   | The right-hand side function failed at the first call.     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_REPTD_RHSFUNC_ERR`      | -10  | The right-hand side function had repeated recoverable      |
   |                                     |      | errors.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_UNREC_RHSFUNC_ERR`      | -11  | The right-hand side function had a recoverable error, but  |
   |                                     |      | no recovery is possible.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RTFUNC_FAIL`            | -12  | The rootfinding function failed in an unrecoverable        |
   |                                     |      | manner.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LFREE_FAIL`             | -13  | The linear solver's memory deallocation function failed.   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSINIT_FAIL`          | -14  | The mass matrix linear solver's initialization function    |
   |                                     |      | failed.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSSETUP_FAIL`         | -15  | The mass matrix linear solver's setup function failed in   |
   |                                     |      | an unrecoverable manner.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSSOLVE_FAIL`         | -16  | The mass matrix linear solver's solve function failed in   |
   |                                     |      | an unrecoverable manner.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSFREE_FAIL`          | -17  | The mass matrix linear solver's memory deallocation        |
   |                                     |      | function failed.                                           |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSMULT_FAIL`          | -18  | The mass matrix-vector product function failed.            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONSTR_FAIL`            | -19  | The inequality constraint test failed repeatedly or        |
   |                                     |      | failed with the minimum step size.                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MEM_FAIL`               | -20  | A memory allocation failed.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MEM_NULL`               | -21  | The ``arkode_mem`` argument was ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ILL_INPUT`              | -22  | One of the function inputs is illegal.                     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NO_MALLOC`              | -23  | The ARKODE memory block was not allocated by               |
   |                                     |      | a call to :c:func:`ARKStepCreate`,                         |
   |                                     |      | :c:func:`ERKStepCreate`, or :c:func:`MRIStepCreate`.       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_K`                  | -24  | The derivative order :math:`k` is larger than allowed.     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_T`                  | -25  | The time :math:`t` is outside the last step taken.         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_DKY`                | -26  | The output derivative vector is ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_CLOSE`              | -27  | The output and initial times are too close to each other.  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_VECTOROP_ERR`           | -28  | An error occurred when calling an :c:type:`N_Vector`       |
   |                                     |      | routine.                                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_INIT_FAIL`          | -29  | An error occurred when initializing a SUNNonlinSol module. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_SETUP_FAIL`         | -30  | A non-recoverable error occurred when setting up a         |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_SETUP_RECVR`        | -31  | A recoverable error occurred when setting up a             |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_OP_ERR`             | -32  | An error occurred when calling a set/get routine in a      |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INNERSTEP_ATTACH_ERR`   | -33  | An error occurred when attaching the inner stepper module. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INNERSTEP_FAIL`         | -34  | An error occurred in the inner stepper module.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_PREINNERFN_FAIL`        | -35  | An error occurred in the MRIStep pre inner integrator      |
   |                                     |      | function.                                                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_POSTINNERFN_FAIL`       | -36  | An error occurred in the MRIStep post inner integrator     |
   |                                     |      | function.                                                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INTERP_FAIL`            | -40  | An error occurred in the ARKODE polynomial interpolation   |
   |                                     |      | module.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INVALID_TABLE`          | -41  | An invalid Butcher or MRI table was encountered.           |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONTEXT_ERR`            | -42  | An error occurred with the SUNDIALS context object         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_FAIL`             | -43  | An error occurred in computing the relaxation parameter    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_MEM_FAIL`         | -44  | The relaxation memory structure is ``NULL``                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_FUNC_FAIL`        | -45  | The relaxation function returned an unrecoverable error    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RELAX_JAC_FAIL`         | -46  | The relaxation Jacobian function returned an unrecoverable |
   |                                     |      | error                                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONTROLLER_ERR`         | -47  | An error with a SUNAdaptController object was encountered. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_STEPPER_UNSUPPORTED`    | -48  | An operation was not supported by the current              |
   |                                     |      | time-stepping module.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_DOMEIG_FAIL`            | -49  | The dominant eigenvalue function failed. It is either not  |
   |                                     |      | provided or returns an illegal value.                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MAX_STAGE_LIMIT_FAIL`   | -50  | Stepper failed to achieve stable results. Either reduce    |
   |                                     |      | the step size or increase the stage_max_limit              |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_SUNSTEPPER_ERR`         | -51  | An error occurred in the SUNStepper module.                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_STEP_DIRECTION_ERR`     | -52  | An error occurred changing the step direction.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ADJ_CHECKPOINT_FAIL`    | -53  | An occurred when checkpointing a state during the adjoint  |
   |                                     |      | integration.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ADJ_RECOMPUTE_FAIL`     | -54  | An occurred recomputing steps during the adjoint           |
   |                                     |      | integration.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_SUNADJSTEPPER_ERR`      | -55  | An error occurred in the SUNAdjStepper module.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_DEE_FAIL`               | -56  | An error occurred in the SUNDomEigEstimator module.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_UNRECOGNIZED_ERROR`     | -99  | An unknown error was encountered.                          |
   +-------------------------------------+------+------------------------------------------------------------+
   | **ARKLS linear solver module output constants**                                                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUCCESS`              | 0    | Successful function return.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MEM_NULL`             | -1   | The ``arkode_mem`` argument was ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_LMEM_NULL`            | -2   | The ARKLS linear solver interface has not been             |
   |                                     |      | initialized.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_ILL_INPUT`            | -3   | The ARKLS solver interface is not compatible with          |
   |                                     |      | the current :c:type:`N_Vector` module, or an input value   |
   |                                     |      | was illegal.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MEM_FAIL`             | -4   | A memory allocation request failed.                        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_PMEM_NULL`            | -5   | The preconditioner module has not been initialized.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSMEM_NULL`         | -6   | The ARKLS mass-matrix linear solver interface has not been |
   |                                     |      | initialized.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_JACFUNC_UNRECVR`      | -7   | The Jacobian function failed in an unrecoverable manner.   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_JACFUNC_RECVR`        | -8   | The Jacobian function had a recoverable error.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSFUNC_UNRECVR`     | -9   | The mass matrix function failed in an unrecoverable        |
   |                                     |      | manner.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSFUNC_RECVR`       | -10  | The mass matrix function had a recoverable error.          |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUNMAT_FAIL`          | -11  | An error occurred with the current :c:type:`SUNMatrix`     |
   |                                     |      | module.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUNLS_FAIL`           | -12  | An error occurred with the current                         |
   |                                     |      | :c:type:`SUNLinearSolver` module.                          |
   +-------------------------------------+------+------------------------------------------------------------+

.. c:enum:: ARKRelaxSolver

   Nonlinear solver identifiers used to specify the method for solving
   :eq:`ARKODE_RELAX_NLS` when relaxation is enabled.

   .. c:enumerator:: ARK_RELAX_NEWTON

      Newton's method

   .. c:enumerator:: ARK_RELAX_BRENT

      Brent's method
..
   Commented-out table rows:

      +-------------------------------------+------+------------------------------------------------------------+
      | :index:`ARK_POSTPROCESS_STEP_FAIL`  | -37  | An error occurred when calling the user-provided           |
      |                                     |      | step-based :c:func:`ARKPostProcessFn` routine.             |
      +-------------------------------------+------+------------------------------------------------------------+
      | :index:`ARK_POSTPROCESS_STAGE_FAIL` | -38  | An error occurred when calling the user-provided           |
      |                                     |      | stage-based :c:func:`ARKPostProcessFn` routine.            |
      +-------------------------------------+------+------------------------------------------------------------+
