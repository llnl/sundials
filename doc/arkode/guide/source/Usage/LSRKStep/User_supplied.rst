.. ----------------------------------------------------------------
   Programmer(s): Mustafa Aggul @ UMBC
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

.. _LSRKSTEP.Usage.UserSupplied:

User-supplied functions
=============================

In addition to the required :c:type:`ARKRhsFn` arguments that define the IVP,
RKL and RKC methods additionally require an :c:type:`ARKDomEigFn` function to
estimate the dominant eigenvalue.




.. _LSRKStep.Usage.dom_eig:

The dominant eigenvalue estimation
----------------------------------

When running LSRKStep with either the RKC or RKL methods, the user must supply
a dominant eigenvalue estimation function of type :c:type:`ARKDomEigFn`:

.. c:type:: int (*ARKDomEigFn)(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR, sunrealtype* lambdaI, void* user_data, N_Vector temp1, N_Vector temp2, N_Vector temp3);

   These functions compute the dominant eigenvalue of the Jacobian of the ODE 
   right-hand side for a given value of the independent variable :math:`t` and 
   state vector :math:`y`.

   :param t: the current value of the independent variable.
   :param y: the current value of the dependent variable vector.
   :param fn: the current value of the vector :math:`f(t,y)`.
   :param lambdaR: The real part of the dominant eigenvalue.
   :param lambdaI: The imaginary part of the dominant eigenvalue.   
   :param user_data: the `user_data` pointer that was passed to
                     :c:func:`ARKodeSetUserData`.
   :param tmp*: pointers to memory allocated to
                variables of type ``N_Vector`` which can be used by an
                ARKDomEigFn as temporary storage or work space.

   :return: An *ARKDomEigFn* should return 0 if successful and any nonzero for a failure.