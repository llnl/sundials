.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
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

.. _ARKodeButcherTable:

==============================
Butcher Table Data Structure
==============================

To store a Butcher table, :math:`B`, defining a Runge--Kutta method ARKODE
provides the :c:type:`ARKodeButcherTable` type and several related utility
routines. We use the following notation

.. math::

   B \; \equiv \;
   \begin{array}{r|c}
     c & A \\
     \hline
     q & b \\
     p & \tilde{b}
   \end{array}
   \; = \;
   \begin{array}{c|cccc}
   c_1    & a_{1,1} & \cdots & a_{1,s-1} & a_{1,s} \\
   c_2    & a_{2,1} & \cdots & a_{2,s-1} & a_{2,s} \\
   \vdots & \vdots  & \vdots & \vdots    & \vdots  \\
   c_s    & a_{s,1} & \cdots & a_{s,s-1} & a_{s,s} \\
   \hline
   q      & b_1         & \cdots & b_{s-1}         & b_s \\
   p      & \tilde{b}_1 & \cdots & \tilde{b}_{s-1} & \tilde{b}_s
   \end{array}.

An :c:type:`ARKodeButcherTable` is a pointer to the
:c:struct:`ARKodeButcherTableMem` structure:

.. c:type:: ARKodeButcherTableMem* ARKodeButcherTable

.. c:struct:: ARKodeButcherTableMem

   Structure for storing a Butcher table

   .. c:member:: int q

      The method order of accuracy

   .. c:member:: int p

      The embedding order of accuracy, typically :math:`q = p + 1`

   .. c:member:: int stages

      The number of stages in the method, :math:`s`

   .. c:member:: sunrealtype **A

      The method coefficients :math:`A \in \mathbb{R}^s`

   .. c:member:: sunrealtype *c

      The method abscissa :math:`c \in \mathbb{R}^s`

   .. c:member:: sunrealtype *b

      The method coefficients :math:`b \in \mathbb{R}^s`

   .. c:member:: sunrealtype *d

      The method embedding coefficients :math:`\tilde{b} \in \mathbb{R}^s`

.. _ARKodeButcherTable.Functions:

ARKodeButcherTable functions
-----------------------------

.. _ARKodeButcherTable.FunctionsTable:
.. table:: ARKodeButcherTable functions

   +--------------------------------------------------+------------------------------------------------------------+
   | **Function name**                                | **Description**                                            |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_LoadERK()`           | Retrieve a given explicit Butcher table by its unique ID   |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_LoadERKByName()`     | Retrieve a given explicit Butcher table by its unique name |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_ERKIDToName()`       | Convert an explicit Butcher table ID to its name           |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_LoadDIRK()`          | Retrieve a given implicit Butcher table by its unique ID   |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_LoadDIRKByName()`    | Retrieve a given implicit Butcher table by its unique name |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_DIRKIDToName()`      | Convert an implicit Butcher table ID to its name           |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Alloc()`             | Allocate an empty Butcher table                            |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Create()`            | Create a new Butcher table                                 |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Copy()`              | Create a copy of a Butcher table                           |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Space()`             | Get the Butcher table real and integer workspace size      |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Free()`              | Deallocate a Butcher table                                 |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Write()`             | Write the Butcher table to an output file                  |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_IsStifflyAccurate()` | Determine if ``A[stages - 1][i] == b[i]``                  |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_CheckOrder()`        | Check the order of a Butcher table                         |
   +--------------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_CheckARKOrder()`     | Check the order of an ARK pair of Butcher tables           |
   +--------------------------------------------------+------------------------------------------------------------+

.. c:function:: ARKodeButcherTable ARKodeButcherTable_LoadERK(ARKODE_ERKTableID emethod)

   Retrieves a specified explicit Butcher table. The prototype for this
   function, as well as the integer names for each provided method, are defined
   in the header file ``arkode/arkode_butcher_erk.h``.  For further information
   on these tables and their corresponding identifiers, see :numref:`Butcher`.

   **Arguments:**
      * *emethod* -- integer input specifying the given Butcher table.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *emethod* was invalid.


.. c:function:: ARKodeButcherTable ARKodeButcherTable_LoadERKByName(const char *emethod)

   Retrieves a specified explicit Butcher table. The prototype for this
   function, as well as the names for each provided method, are defined in the
   header file ``arkode/arkode_butcher_erk.h``.  For further information on
   these tables and their corresponding names, see :numref:`Butcher`.

   **Arguments:**
      * *emethod* -- name of the Butcher table.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *emethod* was invalid or ``"ARKODE_ERK_NONE"``.

   **Notes:**
      This function is case sensitive.

.. c:function:: const char* ARKodeButcherTable_ERKIDToName(ARKODE_ERKTableID emethod)

   Converts a specified explicit Butcher table ID to a string of the same name.
   The prototype for this function, as well as the integer names for each
   provided method, are defined in the header file
   ``arkode/arkode_butcher_erk.h``.  For further information on these tables and
   their corresponding identifiers, see :numref:`Butcher`.

   **Arguments:**
      * *emethod* -- integer input specifying the given Butcher table.

   **Return value:**
      * The name associated with *emethod*.
      * ``NULL`` pointer if *emethod* was invalid.
   
   .. versionadded:: 6.1.0

.. c:function:: ARKodeButcherTable ARKodeButcherTable_LoadDIRK(ARKODE_DIRKTableID imethod)

   Retrieves a specified diagonally-implicit Butcher table. The prototype for
   this function, as well as the integer names for each provided method, are
   defined in the header file ``arkode/arkode_butcher_dirk.h``.  For further
   information on these tables and their corresponding identifiers, see
   :numref:`Butcher`.

   **Arguments:**
      * *imethod* -- integer input specifying the given Butcher table.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *imethod* was invalid.


.. c:function:: ARKodeButcherTable ARKodeButcherTable_LoadDIRKByName(const char *imethod)

   Retrieves a specified diagonally-implicit Butcher table. The prototype for
   this function, as well as the names for each provided method, are defined in
   the header file ``arkode/arkode_butcher_dirk.h``.  For further information
   on these tables and their corresponding names, see :numref:`Butcher`.

   **Arguments:**
      * *imethod* -- name of the Butcher table.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *imethod* was invalid or ``"ARKODE_DIRK_NONE"``.

   **Notes:**
      This function is case sensitive.


.. c:function:: const char* ARKodeButcherTable_DIRKIDToName(ARKODE_DIRKTableID imethod)

   Converts a specified diagonally-implicit Butcher table ID to a string of the
   same name. The prototype for this function, as well as the integer names for
   each provided method, are defined in the header file
   ``arkode/arkode_butcher_dirk.h``.  For further information on these tables
   and their corresponding identifiers, see :numref:`Butcher`.

   **Arguments:**
      * *imethod* -- integer input specifying the given Butcher table.

   **Return value:**
      * The name associated with *imethod*.
      * ``NULL`` pointer if *imethod* was invalid.
   
   .. versionadded:: 6.1.0


.. c:function:: ARKodeButcherTable ARKodeButcherTable_Alloc(int stages, sunbooleantype embedded)

   Allocates an empty Butcher table.

   **Arguments:**
      * *stages* -- the number of stages in the Butcher table.
      * *embedded* -- flag denoting whether the Butcher table has an embedding
        (``SUNTRUE``) or not (``SUNFALSE``).

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *stages* was invalid or an allocation error occurred.

.. c:function:: ARKodeButcherTable ARKodeButcherTable_Create(int s, int q, int p, sunrealtype *c, sunrealtype *A, sunrealtype *b, sunrealtype *d)

   Allocates a Butcher table and fills it with the given values.

   **Arguments:**
      * *s* -- number of stages in the RK method.
      * *q* -- global order of accuracy for the RK method.
      * *p* -- global order of accuracy for the embedded RK method.
      * *c* -- array (of length *s*) of stage times for the RK method.
      * *A* -- array of coefficients defining the RK stages. This should be
        stored as a 1D array of size *s*s*, in row-major order.
      * *b* -- array of coefficients (of length *s*) defining the time step solution.
      * *d* -- array of coefficients (of length *s*) defining the embedded solution.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *stages* was invalid or an allocation error occurred.

   **Notes:**
      If the method does not have an embedding then *d* should be
      ``NULL`` and *p* should be equal to zero.

      .. warning::
         When calling this function from Fortran, it is important to note that ``A`` is expected
         to be in row-major ordering.

.. c:function:: ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B)

   Creates copy of the given Butcher table.

   **Arguments:**
      * *B* -- the Butcher table to copy.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer an allocation error occurred.

.. c:function:: void ARKodeButcherTable_Space(ARKodeButcherTable B, sunindextype *liw, sunindextype *lrw)

   Get the real and integer workspace size for a Butcher table.

   **Arguments:**
      * *B* -- the Butcher table.
      * *lenrw* -- the number of ``sunrealtype`` values in the Butcher table workspace.
      * *leniw* -- the number of integer values in the Butcher table workspace.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_MEM_NULL* if the Butcher table memory was ``NULL``.

   .. deprecated:: 6.3.0

      Work space functions will be removed in version 8.0.0.

.. c:function:: void ARKodeButcherTable_Free(ARKodeButcherTable B)

   Deallocate the Butcher table memory.

   **Arguments:**
      * *B* -- the Butcher table.

.. c:function:: void ARKodeButcherTable_Write(ARKodeButcherTable B, FILE *outfile)

   Write the Butcher table to the provided file pointer.

   **Arguments:**
      * *B* -- the Butcher table.
      * *outfile* -- pointer to use for printing the Butcher table.

   **Notes:**
      The *outfile* argument can be ``stdout`` or ``stderr``, or it
      may point to a specific file created using ``fopen``.

.. c:function:: sunbooleantype ARKodeButcherTable_IsStifflyAccurate(ARKodeButcherTable B)

   Determine if the table satisfies ``A[stages - 1][i] == b[i]``

   **Arguments:**
      * *B* -- the Butcher table.

   **Returns**
      * ``SUNTRUE`` if the method is "stiffly accurate", otherwise returns
        ``SUNFALSE``

   .. versionadded:: v5.7.0

.. c:function:: int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int* q, int* p, FILE* outfile)

   Determine the analytic order of accuracy for the specified Butcher
   table. The analytic (necessary) conditions are checked up to order 6. For
   orders greater than 6 the Butcher simplifying (sufficient) assumptions are
   used.

   **Arguments:**
      * *B* -- the Butcher table.
      * *q* -- the measured order of accuracy for the method.
      * *p* -- the measured order of accuracy for the embedding; 0 if the
        method does not have an embedding.
      * *outfile* -- file pointer for printing results; ``NULL`` to suppress
        output.

   **Return value:**
      * *0* -- success, the measured vales of *q* and *p* match the values of
        *q* and *p* in the provided Butcher tables.
      * *1* -- warning, the values of *q* and *p* in the provided Butcher tables
        are *lower* than the measured values, or the measured values achieve the
        *maximum order* possible with this function and the values of *q* and
        *p* in the provided Butcher tables table are higher.
      * *-1* -- failure, the values of *q* and *p* in the provided Butcher tables
        are *higher* than the measured values.
      * *-2* -- failure, the input Butcher table or critical table contents are
        ``NULL``.

   **Notes:**
      For embedded methods, if the return flags for *q* and *p* would
      differ, failure takes precedence over warning, which takes precedence over
      success.


.. c:function:: int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1, ARKodeButcherTable B2, int *q, int *p, FILE *outfile)

   Determine the analytic order of accuracy (up to order 6) for a specified
   ARK pair of Butcher tables.

   **Arguments:**
      * *B1* -- a Butcher table in the ARK pair.
      * *B2* -- a Butcher table in the ARK pair.
      * *q* -- the measured order of accuracy for the method.
      * *p* -- the measured order of accuracy for the embedding; 0 if the
        method does not have an embedding.
      * *outfile* -- file pointer for printing results; ``NULL`` to suppress
        output.

   **Return value:**
      * *0* -- success, the measured vales of *q* and *p* match the values of
        *q* and *p* in the provided Butcher tables.
      * *1* -- warning, the values of *q* and *p* in the provided Butcher tables
        are *lower* than the measured values, or the measured values achieve the
        *maximum order* possible with this function and the values of *q* and
        *p* in the provided Butcher tables table are higher.
      * *-1* -- failure, the input Butcher tables or critical table contents are
        ``NULL``.

   **Notes:**
      For embedded methods, if the return flags for *q* and *p* would
      differ, warning takes precedence over success.
