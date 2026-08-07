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

.. _SUNMemory:

###########################
Tools for Memory Management
###########################

To support applications which leverage memory pools, or utilize a memory
abstraction layer, sundials provides a set of utilities we will
collectively refer to as the SUNMemoryHelper API. The goal of this API
is to allow users to leverage operations defined by native sundials
data structures while allowing the user to have finer-grained control of
the memory management.

.. toctree::
   :maxdepth: 1

   SUNMemory_links.rst
