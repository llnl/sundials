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

.. _ARKODE.Usage.MRIStep.MRIStepCoupling:

MRI Coupling Coefficients Data Structure
----------------------------------------

MRIStep supplies several built-in MIS, MRI-GARK, and IMEX-MRI-GARK methods, see
:numref:`ARKODE.Usage.MRIStep.MRIStepCoupling.Tables` for the current set of
coupling tables and their corresponding identifiers. Additionally, a user may
supply a custom set of slow-to-fast time scale coupling coefficients by
constructing a coupling table and attaching it with
:c:func:`MRIStepSetCoupling`. A given MRI coupling table can encode any of
the MRI methods supported by MRIStep.  The family of MRI method encoded
by the table is determined by an enumerated type, :c:enum:`MRISTEP_METHOD_TYPE`:

.. c:enum:: MRISTEP_METHOD_TYPE

   The MRI method family encoded by a :c:type:`MRIStepCoupling` table

   .. c:enumerator:: MRISTEP_EXPLICIT

      An explicit MRI-GARK method (does not support a slow implicit operator, :math:`f^I`).

   .. c:enumerator:: MRISTEP_IMPLICIT

      An implicit MRI-GARK method (does not support a slow explicit operator, :math:`f^E`).

   .. c:enumerator:: MRISTEP_IMEX

      An IMEX-MRK-GARK method.

   .. c:enumerator:: MRISTEP_MERK

      A explicit MERK method (does not support a slow implicit operator, :math:`f^I`).

   .. c:enumerator:: MRISTEP_SR

      An IMEX-MRI-SR method.


The MRI coupling tables themselves are stored in an
:c:func:`MRIStepCoupling` object which is a pointer to a
:c:struct:`MRIStepCouplingMem` structure:

.. c:type:: MRIStepCouplingMem *MRIStepCoupling

.. c:struct::  MRIStepCouplingMem

   Structure for storing the coupling coefficients defining an MIS, MRI-GARK, or
   IMEX-MRI-GARK method.

   As described in :numref:`ARKODE.Mathematics.MRIStep`, the coupling from the
   slow time scale to the fast time scale is encoded by a vector of slow
   stage time abscissae, :math:`c^S \in \mathbb{R}^{s+1}` and a set of coupling
   tensors :math:`\Gamma\in\mathbb{R}^{(s+1)\times(s+1)\times k}` and
   :math:`\Omega\in\mathbb{R}^{(s+1)\times(s+1)\times k}`.

   .. c:member:: MRISTEP_METHOD_TYPE type

      Flag indicating the type of MRI method encoded by this table.

   .. c:member:: int nmat

      The value of :math:`k` above i.e., number of coupling matrices in
      :math:`\Omega` for the slow-nonstiff terms and/or in :math:`\Gamma` for
      the slow-stiff terms in :eq:`ARKODE_IVP_two_rate`.

   .. c:member:: int stages

      The number of abscissae i.e., :math:`s+1` above.

   .. c:member:: int q

      The method order of accuracy.

   .. c:member:: int p

      The embedding order of accuracy.

   .. c:member:: sunrealtype* c

      An array of length ``[stages]`` containing the slow abscissae :math:`c^S`
      for the method.

   .. c:member:: sunrealtype*** W

      A three-dimensional array with dimensions ``[nmat][stages+1][stages]``
      containing the method's :math:`\Omega` coupling coefficients for the
      slow-nonstiff (explicit) terms in :eq:`ARKODE_IVP_two_rate`.

   .. c:member:: sunrealtype*** G

      A three-dimensional array with dimensions ``[nmat][stages+1][stages]``
      containing the method's :math:`\Gamma` coupling coefficients for the
      slow-stiff (implicit) terms in :eq:`ARKODE_IVP_two_rate`.

   .. c:member:: int ngroup

      Number of stage groups for the method (only relevant for MERK methods).

   .. c:member:: int** group

      A two-dimensional array with dimensions ``[stages][stages]`` that encodes
      which stages should be combined together within fast integration groups
      (only relevant for MERK methods).


.. _ARKODE.Usage.MRIStep.MRIStepCoupling.Functions:

MRIStepCoupling functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the functions for creating and interacting with coupling
tables. The function prototypes and as well as the relevant integer constants
are defined ``arkode/arkode_mristep.h``.

.. _ARKODE.Usage.MRIStep.MRIStepCoupling.Functions.Table:
.. table:: MRIStepCoupling functions

   +-------------------------------------------+--------------------------------------------------------------------+
   | Function name                             | Description                                                        |
   +===========================================+====================================================================+
   | :c:func:`MRIStepCoupling_LoadTable`       | Loads a pre-defined MRIStepCoupling table by ID                    |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_LoadTableByName` | Loads a pre-defined MRIStepCoupling table by name                  |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Alloc`           | Allocate an empty MRIStepCoupling table                            |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Create`          | Create a new MRIStepCoupling table from coefficients               |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_MIStoMRI`        | Create a new MRIStepCoupling table from a Butcher table            |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Copy`            | Create a copy of a MRIStepCoupling table                           |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Space`           | Get the MRIStepCoupling table real and integer workspace sizes     |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Free`            | Deallocate a MRIStepCoupling table                                 |
   +-------------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Write`           | Write the MRIStepCoupling table to an output file                  |
   +-------------------------------------------+--------------------------------------------------------------------+


.. c:function:: MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID method)

   Retrieves a specified coupling table. For further information on the current
   set of coupling tables and their corresponding identifiers, see
   :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling.Tables`.

   :param method: the coupling table identifier.

   :returns:  An :c:type:`MRIStepCoupling` structure if successful. A ``NULL``
                   pointer if *method* was invalid or an allocation error occurred.


.. c:function:: MRIStepCoupling MRIStepCoupling_LoadTableByName(const char* method)

   Retrieves a specified coupling table. For further information on the current
   set of coupling tables and their corresponding name, see
   :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling.Tables`.

   :param method: the coupling table name.

   :returns: An :c:type:`MRIStepCoupling` structure if successful.
             A ``NULL`` pointer if *method* was invalid, *method* was
             ``"ARKODE_MRI_NONE"``, or an allocation error occurred.

   .. note::

      This function is case sensitive.


.. c:function:: MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages, MRISTEP_METHOD_TYPE type)

   Allocates an empty MRIStepCoupling table.

   :param nmat: the value of :math:`k` i.e., number of number of coupling
      matrices in :math:`\Omega` for the slow-nonstiff terms and/or in
      :math:`\Gamma` for the slow-stiff terms in :eq:`ARKODE_IVP_two_rate`.
   :param stages: number of stages in the coupling table.
   :param type: the type of MRI method the table will encode.

   :returns: An :c:type:`MRIStepCoupling` structure if successful.
             A ``NULL`` pointer if *stages* or *type* was invalid or an allocation error
             occurred.

   .. note::

      For :c:enumerator:`MRISTEP_EXPLICIT` tables, the *G* and *group* arrays are not allocated.

      For :c:enumerator:`MRISTEP_IMPLICIT` tables, the *W* and *group* arrays are not allocated.

      For :c:enumerator:`MRISTEP_IMEX` tables, the *group* array is not allocated.

      For :c:enumerator:`MRISTEP_MERK` tables, the *G* array is not allocated.

      For :c:enumerator:`MRISTEP_SR` tables, the *group* array is not allocated.

      When allocated, both :math:`\Omega` and :math:`\Gamma`
      are initialized to all zeros, so only nonzero coefficients need to be provided.

      When allocated, all entries in *group* are initialized to ``-1``,
      indicating an unused group and/or the end of a stage group.  Users who
      supply a custom MRISTEP_MERK table should overwrite all active stages in
      each group.  For example the ``ARKODE_MERK32`` method has 4 stages that
      are evolved in 3 groups -- the first group consists of stage 1, the second
      group consists of stages 2 and 4, while the third group consists of
      stage 3.  Thus *ngroup* should equal 3, and *group* should have
      non-default entries

      .. code-block:: C

         C->group[0][0] = 1;
         C->group[1][0] = 2;
         C->group[1][1] = 4;
         C->group[2][0] = 3;

   .. versionchanged:: 6.2.0

      This function now supports a broader range of MRI method types.



.. c:function:: MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages, int q, int p, sunrealtype *W, sunrealtype *G, sunrealtype *c)

   Allocates a coupling table and fills it with the given values.

   This routine can only be used to create coupling tables with type
   ``MRISTEP_EXPLICIT``, ``MRISTEP_IMPLICIT``, or  ``MRISTEP_IMEX``.  The
   routine determines the relevant type based on whether either of the
   arguments *W* and *G* are ``NULL``.  Users who wish to create MRI
   methods of type ``MRISTEP_MERK`` or ``MRISTEP_SR`` must currently
   do so manually.

   The assumed size of the input arrays *W* and *G* depends on the
   input value for the embedding order of accuracy, *p*.

   * Non-embedded methods should be indicated by an input *p=0*, in which
     case *W* and/or *G* should have entries stored as a 1D array of size
     ``nmat * stages * stages``, in row-major order.

   * Embedded methods should be indicated by an input *p>0*, in which
     case *W* and/or *G* should have entries stored as a 1D array of size
     ``nmat * (stages+1) * stages``, in row-major order.  The additional
     "row" is assumed to hold the embedding coefficients.

   :param nmat: the value of :math:`k` i.e., number of number of coupling
      matrices in :math:`\Omega` for the slow-nonstiff terms and/or in
      :math:`\Gamma` for the slow-stiff terms in :eq:`ARKODE_IVP_two_rate`.
   :param stages: number of stages in the method.
   :param q: global order of accuracy for the method.
   :param p: global order of accuracy for the embedded method.
   :param W: array of values defining the explicit coupling coefficients
             :math:`\Omega`. If the slow method is implicit pass ``NULL``.
   :param G: array of values defining the implicit coupling coefficients
             :math:`\Gamma`. If the slow method is explicit pass ``NULL``.
   :param c: array of slow abscissae for the MRI method. The entries should be
             stored as a 1D array of length ``stages``.

   :returns:  An :c:type:`MRIStepCoupling` structure if successful.
              A ``NULL`` pointer if ``stages`` was invalid, an allocation error occurred,
              or the input data arrays are inconsistent with the method type.


.. c:function:: MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B, int q, int p)

   Creates an MRI coupling table for a traditional MIS method based on the slow
   Butcher table *B*.

   The :math:`s`-stage slow Butcher table must have an explicit first stage
   (i.e., :math:`c_1=0` and :math:`A_{1,j}=0` for :math:`1\le j\le s`),
   sorted abscissae (i.e., :math:`c_{i} \ge c_{i-1}` for :math:`2\le i\le s`),
   and a final abscissa value :math:`c_s \le 1`.  In this case, the
   :math:`(s+1)`-stage coupling table is computed as

   .. math::

      \Omega_{i,j,1} \;\text{or}\; \Gamma_{i,j,1} =
      \begin{cases}
      0, & \text{if}\; i=1,\\
      A_{i,j}-A_{i-1,j}, & \text{if}\; 2\le i\le s,\\
      b_{j}-A_{s,j}, & \text{if}\; i= s+1.
      \end{cases}

   and the embedding coefficients (if applicable) are computed as

   .. math::

      \tilde{\Omega}_{i,j,1} \;\text{or}\; \tilde{\Gamma}_{i,j,1} = \tilde{b}_{j}-A_{s,j}.

   We note that only one of :math:`\Omega` or :math:`\Gamma` will
   be filled in. If *B* corresponded to an explicit method, then this routine
   fills :math:`\Omega`; if *B* is diagonally-implicit, then this routine
   inserts redundant "padding" stages to ensure a solve-decoupled structure and
   then uses the above formula to fill :math:`\Gamma`.

   For general slow tables with at least second-order accuracy, the MIS method will
   be second order.  However, if the slow table is at least third order and
   additionally satisfies

   .. math::

      \sum_{i=2}^s (c_i-c_{i-1})(\mathbf{e}_i+\mathbf{e}_{i-1})^T A c + (1-c_s) \left(\frac12 + \mathbf{e}_s^T A c\right) = \frac13,

   where :math:`\mathbf{e}_j` corresponds to the :math:`j`-th column from the
   :math:`s \times s` identity matrix, then the overall MIS method will be third order.

   As a result, the values of *q* and *p* may differ from the method and
   embedding orders of accuracy for the Runge--Kutta method encoded in *B*,
   which is why these arguments should be supplied separately.

   If *p>0* is input, then the table *B* must include embedding coefficients.


   :param B: the :c:type:`ARKodeButcherTable` for the "slow" MIS method.
   :param q: the overall order of the MIS/MRI method.
   :param p: the overall order of the MIS/MRI embedding.

   :returns: An :c:type:`MRIStepCoupling` structure if successful.
             A ``NULL`` pointer if an allocation error occurred.


.. c:function:: MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling C)

   Creates copy of the given coupling table.

   :param C: the coupling table to copy.

   :returns: An :c:type:`MRIStepCoupling` structure if successful.
             A ``NULL`` pointer if an allocation error occurred.


.. c:function:: void MRIStepCoupling_Space(MRIStepCoupling C, sunindextype *liw, sunindextype *lrw)

   Get the real and integer workspace size for a coupling table.

   :param C: the coupling table.
   :param lenrw: the number of ``sunrealtype`` values in the coupling table
                 workspace.
   :param leniw: the number of integer values in the coupling table workspace.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_MEM_NULL: if the Butcher table memory was ``NULL``.

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.


.. c:function:: void MRIStepCoupling_Free(MRIStepCoupling C)

   Deallocate the coupling table memory.

   :param C: the coupling table.


.. c:function:: void MRIStepCoupling_Write(MRIStepCoupling C, FILE *outfile)

   Write the coupling table to the provided file pointer.

   :param C: the coupling table.
   :param outfile: pointer to use for printing the table.

   .. note::

      The *outfile* argument can be ``stdout`` or ``stderr``, or it may point to
      a specific file created using ``fopen``.





.. _ARKODE.Usage.MRIStep.MRIStepCoupling.Tables:

MRI Coupling Tables
^^^^^^^^^^^^^^^^^^^

MRIStep currently includes three classes of coupling tables: those that encode
methods that are explicit at the slow time scale, those that are
diagonally-implicit and solve-decoupled at the slow time scale, and those that
encode methods with an implicit-explicit method at the slow time scale.  We list
the current identifiers, multirate order of accuracy, and relevant references
for each in the tables below. For methods with an implicit component, we also
list the number of implicit solves per step that are required at the slow time
scale.

Each of the coupling tables that are packaged with MRIStep are specified by a
unique ID having type:

.. c:type:: int ARKODE_MRITableID

with values specified for each method below (e.g., ``ARKODE_MIS_KW3``).



.. table:: Explicit MRIStep coupling tables.

   ======================================  ==================  ===============  ==============  =====================
   Table name                              Method Order        Embedding Order  Slow RHS Calls  Reference
   ======================================  ==================  ===============  ==============  =====================
   :index:`ARKODE_MRI_GARK_FORWARD_EULER`  :math:`1^*`         --               1
   :index:`ARKODE_MRI_GARK_ERK22a`         2                   1                2               :cite:p:`Sandu:19`
   :index:`ARKODE_MRI_GARK_ERK22b`         :math:`2^{*\circ}`  1                2               :cite:p:`Sandu:19`
   :index:`ARKODE_MRI_GARK_RALSTON2`       2                   1                2               :cite:p:`Roberts:22`
   :index:`ARKODE_MERK21`                  2                   1                2               :cite:p:`Luan:20`
   :index:`ARKODE_MIS_KW3`                 :math:`3^*`         --               3               :cite:p:`Schlegel:09`
   :index:`ARKODE_MRI_GARK_ERK33a`         :math:`3^{\circ}`   2                3               :cite:p:`Sandu:19`
   :index:`ARKODE_MRI_GARK_RALSTON3`       3                   2                3               :cite:p:`Roberts:22`
   :index:`ARKODE_MERK32`                  3                   2                3               :cite:p:`Luan:20`
   :index:`ARKODE_MRI_GARK_ERK45a`         :math:`4^{*\circ}`  3                5               :cite:p:`Sandu:19`
   :index:`ARKODE_MERK43`                  4                   3                6               :cite:p:`Luan:20`
   :index:`ARKODE_MERK54`                  :math:`5^{A}`       4                10              :cite:p:`Luan:20`
   ======================================  ==================  ===============  ==============  =====================


Notes regarding the above table:

#. The default method for each order when using fixed step sizes is marked with an
   asterisk (:math:`^*`).

#. The default method for each order when using adaptive time stepping is marked
   with a circle (:math:`^\circ`).

#. The "Slow RHS Calls" column corresponds to the number of calls to the slow
   right-hand side function, :math:`f^E`, per time step.

#. Note A: although all MERK methods were derived in :cite:p:`Luan:20` under an assumption
   that the fast function :math:`f^F(t,y)` is linear in :math:`y`, in :cite:p:`Fish:24` it
   was proven that MERK methods also satisfy all nonlinear order conditions up through
   their linear order.  The lone exception is :index:`ARKODE_MERK54`, where it was only
   proven to satisfy all nonlinear conditions up to order 4 (since :cite:p:`Fish:24` did
   not establish the formulas for the order 5 conditions).  All our numerical tests to
   date have shown :index:`ARKODE_MERK54` to achieve fifth order for nonlinear problems,
   and so we conjecture that it also satisfies the nonlinear fifth order conditions.


.. table:: Diagonally-implicit, solve-decoupled MRI-GARK coupling tables. The default
           method for each order when using fixed step sizes is marked with an asterisk
           (:math:`^*`); the default method for each order when using adaptive time
           stepping is marked with a circle (:math:`^\circ`). The "Implicit Solves"
           column corresponds to the number of slow implicit (non)linear solves required
           per time step.

   ==========================================  ==================  ===============  ===============  ==================
   Table name                                  Method Order        Embedding Order  Implicit Solves  Reference
   ==========================================  ==================  ===============  ===============  ==================
   :index:`ARKODE_MRI_GARK_BACKWARD_EULER`     :math:`1^{*\circ}`  --               1
   :index:`ARKODE_MRI_GARK_IRK21a`             :math:`2^{*\circ}`  1                1                :cite:p:`Sandu:19`
   :index:`ARKODE_MRI_GARK_IMPLICIT_MIDPOINT`  2                   --               2
   :index:`ARKODE_MRI_GARK_ESDIRK34a`          :math:`3^{*\circ}`  2                3                :cite:p:`Sandu:19`
   :index:`ARKODE_MRI_GARK_ESDIRK46a`          :math:`4^{*\circ}`  3                5                :cite:p:`Sandu:19`
   ==========================================  ==================  ===============  ===============  ==================


.. table:: Diagonally-implicit, solve-decoupled IMEX-MRI-GARK coupling tables.
           The default method for each order when using fixed step sizes is marked
           with an asterisk (:math:`^*`); the default method for each order when using
           adaptive time stepping is marked with a circle (:math:`^\circ`).  The
           "Implicit Solves" column corresponds to the number of slow implicit
           (non)linear solves required per time step.

   =========================================  =================  ===============  ===============  ===================
   Table name                                 Method Order       Embedding Order  Implicit Solves  Reference
   =========================================  =================  ===============  ===============  ===================
   :index:`ARKODE_IMEX_MRI_GARK_EULER`        :math:`1^*`        --               1
   :index:`ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL`  :math:`2^*`        --               1
   :index:`ARKODE_IMEX_MRI_GARK_MIDPOINT`     2                  --               2
   :index:`ARKODE_IMEX_MRI_SR21`              :math:`2^{\circ}`  1                3                :cite:p:`Fish:24`
   :index:`ARKODE_IMEX_MRI_GARK3a`            :math:`3^*`        --               2                :cite:p:`ChiRen:21`
   :index:`ARKODE_IMEX_MRI_GARK3b`            3                  --               2                :cite:p:`ChiRen:21`
   :index:`ARKODE_IMEX_MRI_SR32`              :math:`3^{\circ}`  2                4                :cite:p:`Fish:24`
   :index:`ARKODE_IMEX_MRI_GARK4`             :math:`4^*`        --               5                :cite:p:`ChiRen:21`
   :index:`ARKODE_IMEX_MRI_SR43`              :math:`4^{\circ}`  3                5                :cite:p:`Fish:24`
   =========================================  =================  ===============  ===============  ===================
