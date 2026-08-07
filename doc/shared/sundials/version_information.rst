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

.. _SUNDIALS.SUNVersionInfo:

Getting Version Information
===========================

SUNDIALS provides additional utilities to all packages, that may be
used to retrieve SUNDIALS version information at runtime.


.. c:function:: int SUNDIALSGetVersion(char *version, int len)

   This routine fills a string with SUNDIALS version information.

   **Arguments:**
      * *version* -- character array to hold the SUNDIALS version information.
      * *len* -- allocated length of the *version* character array.

   **Return value:**
      * 0 if successful
      * -1 if the input string is too short to store the SUNDIALS version

   **Notes:**
      An array of 25 characters should be sufficient to hold
      the version information.



.. c:function:: int SUNDIALSGetVersionNumber(int *major, int *minor, int *patch, char *label, int len)

   This routine sets integers for the SUNDIALS major,
   minor, and patch release numbers and fills a string with the
   release label if applicable.

   **Arguments:**
      * *major* -- SUNDIALS release major version number.
      * *minor* -- SUNDIALS release minor version number.
      * *patch* -- SUNDIALS release patch version number.
      * *label* -- string to hold the SUNDIALS release label.
      * *len* -- allocated length of the *label* character array.

   **Return value:**
      * 0 if successful
      * -1 if the input string is too short to store the SUNDIALS label

   **Notes:**
      An array of 10 characters should be sufficient to hold
      the label information. If a label is not used in the release
      version, no information is copied to *label*.
