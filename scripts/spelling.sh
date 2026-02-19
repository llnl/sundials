#!/bin/bash
# ---------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------------
# This script will use codespell to check for common misspellings
# ---------------------------------------------------------------------------------

codespell \
    --skip="*.git,*.bib,*.eps,*.pdf,*/fmod_int*,*/_themes,*/test/answers" \
    -L "inout,ans,Numer,KnWo,Wit,MaPe,ASAi,crate,htmp,thist,thi,MIS,dout,usin,alph,wQS,delt,ue,Bu,ue,nd,ist,strat" \
    --write-changes
