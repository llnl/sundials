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

.. _Fortran.Introduction:

Introduction
============

An interface module can be accessed with the ``use`` statement, e.g.

.. code-block:: fortran

   use fsundials_core_mod   ! this is needed to access core SUNDIALS types, utilities, and data structures
   use fcvode_mod           ! this is needed to access CVODE functions and types
   use fnvector_openmp_mod  ! this is needed to access the OpenMP implementation of the N_Vector class

and by linking to the Fortran 2003 library in addition to the C library, e.g.
``libsundials_fcore_mod.<so|a>``, ``libsundials_core.<so|a>``,
``libsundials_fnvecpenmp_mod.<so|a>``, ``libsundials_nvecopenmp.<so|a>``,
``libsundials_fcvode_mod.<so|a>`` and ``libsundials_cvode.<so|a>``.
The ``use`` statements mirror the ``#include`` statements needed when using the
C API.

The Fortran 2003 interfaces leverage the ``iso_c_binding`` module and the
``bind(C)`` attribute to closely follow the SUNDIALS C API (modulo language
differences). The SUNDIALS classes, e.g. :c:type:`N_Vector`, are interfaced as
Fortran derived types, and function signatures are matched but with an ``F``
prepending the name, e.g. ``FN_VConst`` instead of :c:func:`N_VConst` or
``FCVodeCreate`` instead of ``CVodeCreate``. Constants are named exactly as they
are in the C API.  Accordingly, using SUNDIALS via the Fortran 2003 interfaces
looks just like using it in C. Some caveats stemming from the language
differences are discussed in :numref:`Fortran.Differences`. A
discussion on the topic of equivalent data types in C and Fortran 2003 is
presented in :numref:`Fortran.DataTypes`.

Further information on the Fortran 2003 interfaces specific to the
:c:type:`N_Vector`, :c:type:`SUNMatrix`, :c:type:`SUNLinearSolver`, and
:c:type:`SUNNonlinearSolver` classes is given alongside the C documentation. For
details on where the Fortran 2003 module (``.mod``) files and libraries are
installed see :numref:`Installation`.

The Fortran 2003 interface modules were generated with SWIG Fortran
:cite:p:`Swig-Fortran`, a fork of SWIG. Users who are interested in the SWIG
code used in the generation process should contact the SUNDIALS development
team.

.. _Fortran.Table:

.. list-table:: List of SUNDIALS Fortran 2003 interface modules
   :align: center
   :header-rows: 1
   :widths: 30 40

   * - **Class/Module**
     - **Fortran 2003 Module Name**
   * - SUNDIALS core
     - ``fsundials_core_mod``
   * - ARKODE
     - ``farkode_mod``
   * - ARKODE::ARKSTEP
     - ``farkode_arkstep_mod``
   * - ARKODE::ERKSTEP
     - ``farkode_erkstep_mod``
   * - ARKODE::MRISTEP
     - ``farkode_mristep_mod``
   * - ARKODE::SPRKSTEP
     - ``farkode_sprkstep_mod``
   * - ARKODE::LSRKSTEP
     - ``farkode_lsrkstep_mod``
   * - ARKODE::SPLITTINGSTEP
     - ``farkode_splittingstep_mod``
   * - ARKODE::FORCINGSTEP
     - ``farkode_forcingstep_mod``
   * - CVODE
     - ``fcvode_mod``
   * - CVODES
     - ``fcvodes_mod``
   * - IDA
     - ``fida_mod``
   * - IDAS
     - ``fidas_mod``
   * - KINSOL
     - ``fkinsol_mod``
   * - NVECTOR_CUDA
     - Not interfaced
   * - NVECTOR_MANVECTOR
     - ``fnvector_manyvector_mod``
   * - NVECTOR_MPIMANVECTOR
     - ``fnvector_mpimanyvector_mod``
   * - NVECTOR_MPIPLUSX
     - ``fnvector_mpiplusx_mod``
   * - NVECTOR_OPENMP
     - ``fnvector_openmp_mod``
   * - NVECTOR_PARALLEL
     - ``fnvector_parallel_mod``
   * - NVECTOR_PARHYP
     - Not interfaced
   * - NVECTOR_PETSC
     - Not interfaced
   * - NVECTOR_PTHREADS
     - ``fnvector_pthreads_mod``
   * - NVECTOR_RAJA
     - Not interfaced
   * - NVECTOR_SERIAL
     - ``fnvector_serial_mod``
   * - NVECTOR_SYCL
     - Not interfaced
   * - SUNADAPTCONTROLLER_IMEXGUS
     - ``fsunadaptcontroller_imexgus_mod``
   * - SUNADAPTCONTROLLER_SODERLIND
     - ``fsunadaptcontroller_soderlind_mod``
   * - SUNADAPTCONTROLLER_MRIHTOL
     - ``fsunadaptcontroller_mrihtol_mod``
   * - SUNADJOINTCHECKPOINTSCHEME_FIXED
     - ``fsunadjointcheckpointscheme_fixed_mod``
   * - SUNDOMEIGEST_ARNOLDI
     - ``fsundomeigest_arnoldi_mod``
   * - SUNDOMEIGEST_POWER
     - ``fsundomeigest_power_mod``
   * - SUNLINSOL_BAND
     - ``fsunlinsol_band_mod``
   * - SUNLINSOL_DENSE
     - ``fsunlinsol_dense_mod``
   * - SUNLINSOL_KLU
     - ``fsunlinsol_klu_mod``
   * - SUNLINSOL_LAPACKBAND
     - Not interfaced
   * - SUNLINSOL_LAPACKDENSE
     - Not interfaced
   * - SUNLINSOL_MAGMADENSE
     - Not interfaced
   * - SUNLINSOL_ONEMKLDENSE
     - Not interfaced
   * - SUNLINSOL_PCG
     - ``fsunlinsol_pcg_mod``
   * - SUNLINSOL_SLUDIST
     - Not interfaced
   * - SUNLINSOL_SLUMT
     - Not interfaced
   * - SUNLINSOL_SPBCGS
     - ``fsunlinsol_spbcgs_mod``
   * - SUNLINSOL_SPFGMR
     - ``fsunlinsol_spfgmr_mod``
   * - SUNLINSOL_SPGMR
     - ``fsunlinsol_spgmr_mod``
   * - SUNLINSOL_SPTFQMR
     - ``fsunlinsol_sptfqmr_mod``
   * - SUNMATRIX_BAND
     - ``fsunmatrix_band_mod``
   * - SUNMATRIX_DENSE
     - ``fsunmatrix_dense_mod``
   * - SUNMATRIX_MAGMADENSE
     - Not interfaced
   * - SUNMATRIX_ONEMKLDENSE
     - Not interfaced
   * - SUNMATRIX_SPARSE
     - ``fsunmatrix_sparse_mod``
   * - SUNNONLINSOL_FIXEDPOINT
     - ``fsunnonlinsol_fixedpoint_mod``
   * - SUNNONLINSOL_NEWTON
     - ``fsunnonlinsol_newton_mod``
   * - SUNNONLINSOL_PETSCSNES
     - Not interfaced


.. _Fortran.Installation:

Installation
------------

The installation procedure for the Fortran interfaces is the same as for the C/C++ core of SUNDIALS, refer
to :numref:`Installation`. The CMake option to turn on the Fortran interfaces in a SUNDIALS build 
is :cmakeop:`BUILD_FORTRAN_MODULE_INTERFACE`. The Spack variant is ``+fortran``.

.. _Fortran.Portability:

Important notes on portability
------------------------------

The SUNDIALS Fortran 2003 interface *should* be compatible with any compiler
supporting the Fortran 2003 ISO standard.

Upon compilation of SUNDIALS, Fortran module (``.mod``) files are generated for
each Fortran 2003 interface. These files are highly compiler specific, and thus
it is almost always necessary to compile a consuming application with the same
compiler that was used to generate the modules.
