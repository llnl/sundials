.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Fortran.CommonIssues:

Common Issues
=============

In this subsection, we list some common issues users run into when using the Fortran
interfaces.


**Strange Segmentation Fault in User-Supplied Functions**

One common issue we have seen trip up users (and even ourselves) has the symptom
of segmentation fault in a user-supplied function (such as the RHS) when trying
to use one of the callback arguments. For example, in the following RHS
function, we will get a segfault on line 21:

.. code-block:: fortran
   :linenos:
   :emphasize-lines: 8, 21

   integer(c_int) function ff(t, yvec, ydotvec, user_data) &
      result(ierr) bind(C)

      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      implicit none

      real(c_double) :: t ! <===== Missing value attribute
      type(N_Vector) :: yvec
      type(N_Vector) :: ydotvec
      type(c_ptr)    :: user_data

      real(c_double) :: e
      real(c_double) :: u, v
      real(c_double) :: tmp1, tmp2
      real(c_double), pointer :: yarr(:)
      real(c_double), pointer :: ydotarr(:)

      ! get N_Vector data arrays
      yarr => FN_VGetArrayPointer(yvec)
      ydotarr => FN_VGetArrayPointer(ydotvec) ! <===== SEGFAULTS HERE

      ! extract variables
      u = yarr(1)
      v = yarr(2)

      ! fill in the RHS function:
      !  [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
      !  [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))]
      tmp1 = (-ONE+u*u-r(t))/(TWO*u)
      tmp2 = (-TWO+v*v-s(t))/(TWO*v)
      ydotarr(1) = ZERO
      ydotarr(2) = e*tmp1 - tmp2 + sdot(t)/(TWO*vtrue(t))

      ! return success
      ierr = 0
      return

   end function


The subtle bug in the code causing the segfault is on line 8. It should read
``real(c_double), value :: t`` instead of ``real(c_double) :: t`` (notice the
``value`` attribute). Fundamental types that are passed by value in C need
the ``value`` attribute.
