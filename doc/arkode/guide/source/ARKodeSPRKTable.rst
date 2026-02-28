.. ----------------------------------------------------------------
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

.. _ARKodeSPRKTable:

===========================
SPRK Method Table Structure
===========================

To store the pair of Butcher tables defining a SPRK method of order :math:`q`
ARKODE provides the :c:type:`ARKodeSPRKTable` type and several related utility
routines. We use the following notation

.. math::

   B \; \equiv \;
   \begin{array}{r|c}
     c & A \\
     \hline
       & b \\
   \end{array}
   \; = \;
   \begin{array}{c|cccc}
   c_1 & 0 & \cdots & 0 & 0 \\
   c_2 & a_1 & 0 & \cdots & \vdots \\
   \vdots & \vdots & \ddots & \ddots & \vdots \\
   c_s & a_1 & \cdots & a_{s-1} & 0 \\
   \hline
   & a_1 & \cdots & a_{s-1} & a_s
   \end{array}
   \qquad
   \qquad
   \hat{B} \; \equiv \;
   \begin{array}{r|c}
     \hat{c} & \hat{A} \\
     \hline
             & \hat{b} \\
   \end{array}
   \; = \;
   \begin{array}{c|cccc}
   \hat{c}_1 & \hat{a}_1 & \cdots & 0 & 0 \\
   \hat{c}_2 & \hat{a}_1 & \hat{a}_2 & \cdots & \vdots \\
   \vdots & \vdots & \ddots & \ddots & \vdots \\
   \hat{c}_s & \hat{a}_1 & \hat{a}_2 & \cdots & \hat{a}_{s} \\
   \hline
   & \hat{a}_1 & \hat{a}_2 & \cdots & \hat{a}_{s}
   \end{array}

where :math:`B` and :math:`\hat{B}` contain the coefficients for the explicit
and diagonally implicit tables, respectively. We use a compact storage of these
coefficients in terms of two arrays, one for :math:`a` values and one for
:math:`\hat{a}` values. The abscissae (only relevant for non-autonomous
problems) are computed dynamically as :math:`c_j = \sum_{i=1}^j a_i` and
:math:`\hat{c}_j = \sum_{i=1}^j \hat{a}_i`, respectively
:cite:p:`Jay:21,Diele:11`. The :c:type:`ARKodeSPRKTable` type is a pointer to
the :c:type:`ARKodeSPRKTableMem` structure:

.. c:type:: ARKodeSPRKTableMem* ARKodeSPRKTable

.. c:type:: ARKodeSPRKTableMem

   Structure representing the SPRK method that holds the method coefficients.

   .. c:member:: int q

      The method order of accuracy.

   .. c:member:: int stages

      The number of stages.

   .. c:member:: sunrealtype* a

      Array of coefficients that generate the explicit Butcher table.
      ``a[i]`` is the coefficient appearing in column i+1.

   .. c:member:: sunrealtype* ahat

      Array of coefficients that generate the diagonally-implicit Butcher table.
      ``ahat[i]`` is the coefficient appearing in column i.


ARKodeSPRKTable functions
---------------------------

.. _ARKodeSPRKTable.FunctionsTable:
.. table:: ARKodeSPRKTable functions

   +----------------------------------------------+------------------------------------------------------------+
   | **Function name**                            | **Description**                                            |
   +==============================================+============================================================+
   | :c:func:`ARKodeSPRKTable_Alloc()`            | Allocate an empty table                                    |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Load()`             | Load SPRK method using an identifier                       |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_LoadByName()`       | Load SPRK method using a string version of the identifier  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Create()`           | Create a new table                                         |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Copy()`             | Create a copy of a table                                   |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Space()`            | Get the table real and integer workspace size              |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeSPRKTable_Free()`             | Deallocate a table                                         |
   +----------------------------------------------+------------------------------------------------------------+


.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Create(int stages, int q, const sunrealtype* a, const sunrealtype* ahat)

   Creates and allocates an :c:type:`ARKodeSPRKTable` with the specified number
   of stages and the coefficients provided.

   :param stages: The number of stages.
   :param q: The order of the method.
   :param a: An array of the coefficients for the ``a`` table.
   :param ahat: An array of the coefficients for the ``ahat`` table.
   :return: :c:type:`ARKodeSPRKTable` for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Alloc(int stages)

   Allocate memory for an :c:type:`ARKodeSPRKTable` with the specified
   number of stages.

   :param stages: The number of stages.
   :return: :c:type:`ARKodeSPRKTable` for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Load(ARKODE_SPRKMethodID id)

   Load the :c:type:`ARKodeSPRKTable` for the specified method ID.

   :param id: The ID of the SPRK method, see :ref:`Butcher.sprk`.
   :return: :c:type:`ARKodeSPRKTable` for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_LoadByName(const char* method)

   Load the :c:type:`ARKodeSPRKTable` for the specified method name.

   :param method: The name of the SPRK method, see :ref:`Butcher.sprk`.
   :return: :c:type:`ARKodeSPRKTable` for the loaded method.

.. c:function:: ARKodeSPRKTable ARKodeSPRKTable_Copy(ARKodeSPRKTable sprk_table)

   Create a copy of the :c:type:`ARKodeSPRKTable`.

   :param sprk_table: The :c:type:`ARKodeSPRKTable` to copy.
   :return: Pointer to the copied :c:type:`ARKodeSPRKTable`.

.. c:function:: void ARKodeSPRKTable_Write(ARKodeSPRKTable sprk_table, FILE* outfile)

   Write the ARKodeSPRKTable out to the file.

   :param sprk_table: The :c:type:`ARKodeSPRKTable` to write.
   :param outfile: The FILE that will be written to.

.. c:function:: void ARKodeSPRKTable_Space(ARKodeSPRKTable sprk_table, sunindextype* liw, sunindextype* lrw)

   Get the workspace sizes required for the :c:type:`ARKodeSPRKTable`.

   :param sprk_table: The :c:type:`ARKodeSPRKTable`.
   :param liw: Pointer to store the integer workspace size.
   :param lrw: Pointer to store the real workspace size.

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.

.. c:function:: void ARKodeSPRKTable_Free(ARKodeSPRKTable sprk_table)

   Free the memory allocated for the :c:type:`ARKodeSPRKTable`.

   :param sprk_table: The :c:type:`ARKodeSPRKTable` to free.

.. c:function:: int ARKodeSPRKTable_ToButcher(ARKodeSPRKTable sprk_table, ARKodeButcherTable* a_ptr, ARKodeButcherTable* b_ptr)

   Convert the :c:type:`ARKodeSPRKTable` to the Butcher table representation.

   :param sprk_table: The :c:type:`ARKodeSPRKTable`.
   :param a_ptr: Pointer to store the explicit Butcher table.
   :param b_ptr: Pointer to store the diagonally-implicit Butcher table.
