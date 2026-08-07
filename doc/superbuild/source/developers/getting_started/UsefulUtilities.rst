..
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

.. UsefulUtilities:

Useful Utilities
================

This section mentions some useful utilities and internal data structures that
you may need or come across in SUNDIALS. We do not document these fully here,
instead, look for comments in the source.

Data Structures
---------------

We have our own implementation of a hash map, ``SUNHashMap``. It is found at
``src/sundials/sundials_hashmap_impl.h`` and
``src/sundials/sundials_hashmap.c``.

We have a C implementation of a ``std::vector`` in ``src/sundials/stl``.

We have an implementation of a hierarchical node object, which can be used to
build things like a JSON tree, at ``src/sundials/sundials_datanode.{h,c}``.
