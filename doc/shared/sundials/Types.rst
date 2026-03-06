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

.. _SUNDIALS.DataTypes:

Data Types
==========

SUNDIALS defines several data types in the header file ``sundials_types.h``.
These types are used in the SUNDIALS API and internally in SUNDIALS. It is
not necessary to use these types in your application, but the type must
be compatible with the SUNDIALS types in the API when calling SUNDIALS functions.
The types that are defined are:

* :c:type:`sunrealtype` -- the floating-point type used by the SUNDIALS packages

* :c:type:`sunindextype` -- the integer type used for vector and matrix indices

* :c:type:`suncountertype` -- the integer type used for counter variables

* :c:type:`sunbooleantype` -- the type used for logic operations within SUNDIALS

* :c:type:`SUNOutputFormat` -- an enumerated type for SUNDIALS output formats

* :c:type:`SUNComm` -- a simple typedef to an `int` when SUNDIALS is built without MPI, or a ``MPI_Comm`` when built with MPI.


Floating point types
--------------------

.. c:type:: sunrealtype

   The type ``sunrealtype`` can be ``float``, ``double``, or ``long double``, with
   the default being ``double``. The user can change the precision of the
   arithmetic used in the SUNDIALS solvers at the configuration stage (see
   :cmakeop:`SUNDIALS_PRECISION`).

Additionally, based on the current precision, ``sundials_types.h`` defines
``SUN_BIG_REAL`` to be the largest value representable as a ``sunrealtype``,
``SUN_SMALL_REAL`` to be the smallest value representable as a ``sunrealtype``, and
``SUN_UNIT_ROUNDOFF`` to be the difference between :math:`1.0` and the minimum
``sunrealtype`` greater than :math:`1.0`.

Within SUNDIALS, real constants are set by way of a macro called ``SUN_RCONST``. It
is this macro that needs the ability to branch on the definition of
``sunrealtype``. In ANSI C, a floating-point constant with no suffix is stored as a
``double``. Placing the suffix "``F``" at the end of a floating point constant
makes it a ``float``, whereas using the suffix "``L``" makes it a ``long
double``. For example,

.. code-block:: c

   #define A 1.0
   #define B 1.0F
   #define C 1.0L

defines ``A`` to be a ``double`` constant equal to :math:`1.0`, ``B`` to be a
``float`` constant equal to :math:`1.0`, and ``C`` to be a ``long double``
constant equal to :math:`1.0`. The macro call ``SUN_RCONST(1.0)`` automatically
expands to ``1.0`` if ``sunrealtype`` is ``double``, to ``1.0F`` if ``sunrealtype`` is
``float``, or to ``1.0L`` if ``sunrealtype`` is ``long double``. SUNDIALS uses the
``SUN_RCONST`` macro internally to declare all of its floating-point constants.

Additionally, SUNDIALS defines several macros for common mathematical functions
*e.g.*, ``fabs``, ``sqrt``, ``exp``, etc. in ``sundials_math.h``. The macros are
prefixed with ``SUNR`` and expand to the appropriate ``C`` function based on the
``sunrealtype``. For example, the macro ``SUNRabs`` expands to the ``C`` function
``fabs`` when ``sunrealtype`` is ``double``, ``fabsf`` when ``sunrealtype`` is
``float``, and ``fabsl`` when ``sunrealtype`` is ``long double``.

A user program which uses the type ``sunrealtype``, the ``SUN_RCONST`` macro, and the
``SUNR`` mathematical function macros is precision-independent except for any
calls to precision-specific library functions. Our example programs use
``sunrealtype``, ``SUN_RCONST``, and the ``SUNR`` macros. Users can, however, use the
type ``double``, ``float``, or ``long double`` in their code (assuming that this
usage is consistent with the typedef for ``sunrealtype``) and call the appropriate
math library functions directly. Thus, a previously existing piece of C or C++
code can use SUNDIALS without modifying the code to use ``sunrealtype``,
``SUN_RCONST``, or the ``SUNR`` macros so long as the SUNDIALS libraries are built
to use the corresponding precision (see :numref:`Installation.Options`).

Integer types used for indexing
-------------------------------

.. c:type:: sunindextype

   The type ``sunindextype`` is used for indexing array entries in SUNDIALS
   modules as well as for storing the total problem size (*e.g.*, vector
   lengths and matrix sizes). During configuration ``sunindextype`` may be
   selected to be either a 32- or 64-bit *signed* integer with the default being
   64-bit (see :cmakeop:`SUNDIALS_INDEX_SIZE`).

When using a 32-bit integer the total problem size is limited to
:math:`2^{31}-1` and with 64-bit integers the limit is :math:`2^{63}-1`. For
users with problem sizes that exceed the 64-bit limit an advanced configuration
option is available to specify the type used for ``sunindextype``
(see :cmakeop:`SUNDIALS_INDEX_TYPE`).

A user program which uses ``sunindextype`` to handle indices will work with both
index storage types except for any calls to index storage-specific external
libraries. Our ``C`` and ``C++`` example programs use ``sunindextype``. Users
can, however, use any compatible type (*e.g.*, ``int``, ``long int``,
``int32_t``, ``int64_t``, or ``long long int``) in their code, assuming that
this usage is consistent with the typedef for ``sunindextype`` on their
architecture. Thus, a previously existing piece of C or C++ code can use
SUNDIALS without modifying the code to use ``sunindextype``, so long as the
SUNDIALS libraries use the appropriate index storage type (for details see
:numref:`Installation.Options`).

Integer type used for counters
------------------------------

.. c:type:: suncountertype

   The type ``suncountertype`` is used for counter variables in SUNDIALS
   (*e.g.*, number of stpes) and is the same as ``long int``.

   .. versionadded:: 7.3.0

Boolean type
------------

.. c:type:: sunbooleantype

   As ANSI C89 (ISO C90) does not have a built-in boolean data type, SUNDIALS
   defines the type ``sunbooleantype`` as an ``int``.

The advantage of using the name sunbooleantype (instead of int) is an increase in
code readability. It also allows the programmer to make a distinction between
int and boolean data. Variables of type ``sunbooleantype`` are intended to have
only the two values: :c:macro:`SUNFALSE` or :c:macro:`SUNTRUE`.

.. c:macro:: SUNFALSE

   False (``0``)

.. c:macro:: SUNTRUE

   True (``1``)

Output formatting type
----------------------

.. c:enum:: SUNOutputFormat

   The enumerated type :c:type:`SUNOutputFormat` defines the enumeration
   constants for SUNDIALS output formats

.. c:enumerator:: SUN_OUTPUTFORMAT_TABLE

   The output will be a table of values

.. c:enumerator:: SUN_OUTPUTFORMAT_CSV

   The output will be a comma-separated list of key and value pairs e.g.,
   ``key1,value1,key2,value2,...``

   .. note::

      The Python module ``tools/suntools`` provides utilities to read and output
      the data from a SUNDIALS CSV output file using the key and value pair
      format.

MPI types
---------

.. c:type:: SUNComm

   A simple typedef to an `int` when SUNDIALS is built without MPI, or a
   ``MPI_Comm`` when built with MPI. This type exists solely to ensure SUNDIALS
   can support MPI and non-MPI builds.

.. c:macro:: SUN_COMM_NULL

   A macro defined as ``0`` when SUNDIALS is built without MPI, or as
   ``MPI_COMM_NULL`` when built with MPI.
