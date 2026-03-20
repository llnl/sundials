#!/bin/bash
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# This script is a utility to upgrade codes from SUNDIALS v5.x.x to v6.0.0.
# It does a few things:
#   1. Updates all SUNDIALS constructor invocations to pass 'SUNCTX_PLACEHOLDER'
#      as the SUNContext argument. Users will need to then manually replace
#      'SUNCTX_PLACEHOLDER' with an actual SUNContext object.
#   2. Updates references to deprecated SUNDIALS constants and types to the new ones.
#
# Usage:
#   ./upgrade-to-sundials-6-from-5.sh <files to update>
# -----------------------------------------------------------------------------

updateConstructors() {
  # All of the constructors updated to take a SUNContext.
  # Obtained with
  #   clang-format include/**/*.h | grep 'SUNContext sunctx)' | sed 's/SUNDIALS_EXPORT//' | awk -F ' ' '{print $2}' > constructors.txt
  # where the .clang-format file sets the ColumnLimit to 1000.
  declare -a constructors=(
    "ARKStepCreate"
    "CVodeCreate"
    "CVodeCreate"
    "ERKStepCreate"
    "IDACreate"
    "IDACreate"
    "KINCreate"
    "MRIStepCreate"
    "N_VMake_Cuda"
    "N_VMake_Hip"
    "N_VMake_MPIManyVector"
    "N_VMake_MPIPlusX"
    "N_VMake_OpenMP"
    "N_VMake_OpenMPDEV"
    "N_VMake_Parallel"
    "N_VMake_ParHyp"
    "N_VMake_Petsc"
    "N_VMake_Pthreads"
    "N_VMake_Raja"
    "N_VMake_Serial"
    "N_VMake_Sycl"
    "N_VMakeManaged_Cuda"
    "N_VMakeManaged_Hip"
    "N_VMakeManaged_Raja"
    "N_VMakeManaged_Sycl"
    "N_VNew_Cuda"
    "N_VNew_Hip"
    "N_VNew_ManyVector"
    "N_VNew_MPIManyVector"
    "N_VNew_OpenMP"
    "N_VNew_OpenMPDEV"
    "N_VNew_Parallel"
    "N_VNew_Pthreads"
    "N_VNew_Raja"
    "N_VNew_Serial"
    "N_VNew_Sycl"
    "N_VNewEmpty"
    "N_VNewEmpty_Cuda"
    "N_VNewEmpty_Hip"
    "N_VNewEmpty_OpenMP"
    "N_VNewEmpty_OpenMPDEV"
    "N_VNewEmpty_Parallel"
    "N_VNewEmpty_ParHyp"
    "N_VNewEmpty_Petsc"
    "N_VNewEmpty_Pthreads"
    "N_VNewEmpty_Raja"
    "N_VNewEmpty_SensWrapper"
    "N_VNewEmpty_Serial"
    "N_VNewEmpty_Sycl"
    "N_VNewEmpty_Trilinos"
    "N_VNewManaged_Cuda"
    "N_VNewManaged_Hip"
    "N_VNewManaged_Raja"
    "N_VNewManaged_Sycl"
    "N_VNewWithMemHelp_Cuda"
    "N_VNewWithMemHelp_Hip"
    "N_VNewWithMemHelp_Raja"
    "N_VNewWithMemHelp_Sycl"
    "N_VSetContext"
    "N_VSetContext_ManyVector"
    "N_VSetContext_MPIManyVector"
    "SUNBandMatrix"
    "SUNBandMatrixStorage"
    "SUNDenseMatrix"
    "SUNLinSol_Band"
    "SUNLinSol_cuSolverSp_batchQR"
    "SUNLinSol_Dense"
    "SUNLinSol_KLU"
    "SUNLinSol_LapackBand"
    "SUNLinSol_LapackDense"
    "SUNLinSol_MagmaDense"
    "SUNLinSol_OneMklDense"
    "SUNLinSol_PCG"
    "SUNLinSol_SPBCGS"
    "SUNLinSol_SPFGMR"
    "SUNLinSol_SPGMR"
    "SUNLinSol_SPTFQMR"
    "SUNLinSol_SuperLUDIST"
    "SUNLinSol_SuperLUMT"
    "SUNLinSolNewEmpty"
    "SUNLinSolSetContext"
    "SUNMatNewEmpty"
    "SUNMatrix_cuSparse_MakeCSR"
    "SUNMatrix_cuSparse_NewBlockCSR"
    "SUNMatrix_cuSparse_NewCSR"
    "SUNMatrix_MagmaDense"
    "SUNMatrix_MagmaDenseBlock"
    "SUNMatrix_OneMklDense"
    "SUNMatrix_OneMklDenseBlock"
    "SUNMatrix_SLUNRloc"
    "SUNMatSetContext"
    "SUNMemoryHelper_Cuda"
    "SUNMemoryHelper_Hip"
    "SUNMemoryHelper_NewEmpty"
    "SUNMemoryHelper_SetContext"
    "SUNMemoryHelper_Sycl"
    "SUNMemoryHelper_Sys"
    "SUNNonlinSol_FixedPoint"
    "SUNNonlinSol_FixedPointSens"
    "SUNNonlinSol_Newton"
    "SUNNonlinSol_NewtonSens"
    "SUNNonlinSol_PetscSNES"
    "SUNNonlinSolNewEmpty"
    "SUNNonlinSolSetContext"
    "SUNSparseMatrix"
    )

  echo "replacing in ${@}"
  for constructor in "${constructors[@]}"
  do
    perl -pi -e 's/'"${constructor}"'\((.*\S.*)\)/'"${constructor}"'($1, SUNCTX_PLACEHOLDER)/g;' ${@}
    perl -pi -e 's/'"${constructor}"'\(\)/'"${constructor}"'(SUNCTX_PLACEHOLDER)/g;' ${@}
  done
}

updateConstantsAndTypes() {
  declare -a map=(
    "PREC_NONE,SUN_PREC_NONE"
    "PREC_LEFT,SUN_PREC_LEFT"
    "PREC_RIGHT,SUN_PREC_RIGHT"
    "PREC_BOTH,SUN_PREC_BOTH"
    "MODIFIED_GS,SUN_MODIFIED_GS"
    "CLASSICAL_GS,SUN_CLASSICAL_GS"
    "ATimesFn,SUNATimesFn"
    "PSetupFn,SUNPSetupFn"
    "PSolveFn,SUNPSolveFn"
    "DlsMat,SUNDlsMat"
    "DENSE_COL,SUNDLS_DENSE_COL"
    "DENSE_ELEM,SUNDLS_DENSE_ELEM"
    "BAND_COL,SUNDLS_BAND_COL"
    "BAND_COL_ELEM,SUNDLS_BAND_COL_ELEM"
    "BAND_ELEM,SUNDLS_BAND_ELEM"
    "SDIRK_2_1_2,ARKODE_SDIRK_2_1_2"
    "BILLINGTON_3_3_2,ARKODE_BILLINGTON_3_3_2"
    "TRBDF2_3_3_2,ARKODE_TRBDF2_3_3_2"
    "KVAERNO_4_2_3,ARKODE_KVAERNO_4_2_3"
    "ARK324L2SA_DIRK_4_2_3,ARKODE_ARK324L2SA_DIRK_4_2_3"
    "CASH_5_2_4,ARKODE_CASH_5_2_4"
    "CASH_5_3_4,ARKODE_CASH_5_3_4"
    "SDIRK_5_3_4,ARKODE_SDIRK_5_3_4"
    "KVAERNO_5_3_4,ARKODE_KVAERNO_5_3_4"
    "ARK436L2SA_DIRK_6_3_4,ARKODE_ARK436L2SA_DIRK_6_3_4"
    "KVAERNO_7_4_5,ARKODE_KVAERNO_7_4_5"
    "ARK548L2SA_DIRK_8_4_5,ARKODE_ARK548L2SA_DIRK_8_4_5"
    "ARK437L2SA_DIRK_7_3_4,ARKODE_ARK437L2SA_DIRK_7_3_4"
    "ARK548L2SAb_DIRK_8_4_5,ARKODE_ARK548L2SAb_DIRK_8_4_5"
    "MIN_DIRK_NUM,ARKODE_MIN_DIRK_NUM"
    "MAX_DIRK_NUM,ARKODE_MAX_DIRK_NUM"
    "MIS_KW3,ARKODE_MIS_KW3"
    "MRI_GARK_ERK33a,ARKODE_MRI_GARK_ERK33a"
    "MRI_GARK_ERK45a,ARKODE_MRI_GARK_ERK45a"
    "MRI_GARK_IRK21a,ARKODE_MRI_GARK_IRK21a"
    "MRI_GARK_ESDIRK34a,ARKODE_MRI_GARK_ESDIRK34a"
    "MRI_GARK_ESDIRK46a,ARKODE_MRI_GARK_ESDIRK46a"
    "IMEX_MRI_GARK3a,ARKODE_IMEX_MRI_GARK3a"
    "IMEX_MRI_GARK3b,ARKODE_IMEX_MRI_GARK3b"
    "IMEX_MRI_GARK4,ARKODE_IMEX_MRI_GARK4"
    "MIN_MRI_NUM,ARKODE_MIN_MRI_NUM"
    "MAX_MRI_NUM,ARKODE_MAX_MRI_NUM"
    "DEFAULT_MRI_TABLE_3,MRISTEP_DEFAULT_TABLE_3"
    "DEFAULT_EXPL_MRI_TABLE_3,MRISTEP_DEFAULT_EXPL_TABLE_3"
    "DEFAULT_EXPL_MRI_TABLE_4,MRISTEP_DEFAULT_EXPL_TABLE_4"
    "DEFAULT_IMPL_SD_TABLE_2,MRISTEP_DEFAULT_IMPL_SD_TABLE_2"
    "DEFAULT_IMPL_SD_TABLE_3,MRISTEP_DEFAULT_IMPL_SD_TABLE_3"
    "DEFAULT_IMPL_SD_TABLE_4,MRISTEP_DEFAULT_IMPL_SD_TABLE_4"
    "DEFAULT_IMEX_SD_TABLE_3,MRISTEP_DEFAULT_IMEX_SD_TABLE_3"
    "DEFAULT_IMEX_SD_TABLE_4,MRISTEP_DEFAULT_IMEX_SD_TABLE_4"
    "HEUN_EULER_2_1_2,ARKODE_HEUN_EULER_2_1_2"
    "BOGACKI_SHAMPINE_4_2_3,ARKODE_BOGACKI_SHAMPINE_4_2_3"
    "ARK324L2SA_ERK_4_2_3,ARKODE_ARK324L2SA_ERK_4_2_3"
    "ZONNEVELD_5_3_4,ARKODE_ZONNEVELD_5_3_4"
    "ARK436L2SA_ERK_6_3_4,ARKODE_ARK436L2SA_ERK_6_3_4"
    "SAYFY_ABURUB_6_3_4,ARKODE_SAYFY_ABURUB_6_3_4"
    "CASH_KARP_6_4_5,ARKODE_CASH_KARP_6_4_5"
    "FEHLBERG_6_4_5,ARKODE_FEHLBERG_6_4_5"
    "DORMAND_PRINCE_7_4_5,ARKODE_DORMAND_PRINCE_7_4_5"
    "ARK548L2SA_ERK_8_4_5,ARKODE_ARK548L2SA_ERK_8_4_5"
    "VERNER_8_5_6,ARKODE_VERNER_8_5_6"
    "FEHLBERG_13_7_8,ARKODE_FEHLBERG_13_7_8"
    "KNOTH_WOLKE_3_3,ARKODE_KNOTH_WOLKE_3_3"
    "ARK437L2SA_ERK_7_3_4,ARKODE_ARK437L2SA_ERK_7_3_4"
    "ARK548L2SAb_ERK_8_4_5,ARKODE_ARK548L2SAb_ERK_8_4_5"
    "MIN_ERK_NUM,ARKODE_MIN_ERK_NUM"
    "MAX_ERK_NUM,ARKODE_MAX_ERK_NUM"
    "DEFAULT_ERK_2,ARKSTEP_DEFAULT_ERK_2 if using ARKStep or ERKSTEP_DEFAULT_2 if using ERKStep"
    "DEFAULT_ERK_3,ARKSTEP_DEFAULT_ERK_3 if using ARKStep or ERKSTEP_DEFAULT_3 if using ERKStep"
    "DEFAULT_ERK_4,ARKSTEP_DEFAULT_ERK_4 if using ARKStep or ERKSTEP_DEFAULT_4 if using ERKStep"
    "DEFAULT_ERK_5,ARKSTEP_DEFAULT_ERK_5 if using ARKStep or ERKSTEP_DEFAULT_5 if using ERKStep"
    "DEFAULT_ERK_6,ARKSTEP_DEFAULT_ERK_6 if using ARKStep or ERKSTEP_DEFAULT_6 if using ERKStep"
    "DEFAULT_ERK_8,ARKSTEP_DEFAULT_ERK_8 if using ARKStep or ERKSTEP_DEFAULT_8 if using ERKStep"
    "DEFAULT_DIRK_2,ARKSTEP_DEFAULT_DIRK_2"
    "DEFAULT_DIRK_3,ARKSTEP_DEFAULT_DIRK_3"
    "DEFAULT_DIRK_4,ARKSTEP_DEFAULT_DIRK_4"
    "DEFAULT_DIRK_5,ARKSTEP_DEFAULT_DIRK_5"
    "DEFAULT_ARK_ETABLE_3,ARKSTEP_DEFAULT_ARK_ETABLE_3"
    "DEFAULT_ARK_ETABLE_4,ARKSTEP_DEFAULT_ARK_ETABLE_4"
    "DEFAULT_ARK_ETABLE_5,ARKSTEP_DEFAULT_ARK_ETABLE_4"
    "DEFAULT_ARK_ITABLE_3,ARKSTEP_DEFAULT_ARK_ITABLE_3"
    "DEFAULT_ARK_ITABLE_4,ARKSTEP_DEFAULT_ARK_ITABLE_4"
    "DEFAULT_ARK_ITABLE_5,ARKSTEP_DEFAULT_ARK_ITABLE_5"
  )
  for pair in "${map[@]}";
  do
      IFS=',' read item1 item2 <<< "${pair}"
      perl -pi -e 's/\b'"${item1}"'\b/'"${item2}"'/g;' ${@}
  done
}

while true; do
    echo ""
    read -p "> This script will edit files in-place. Are you sure you want to proceed? [Yy/Nn or ls] " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        ls ) echo "> The following files will be edited:" && printf "  %s\n" "${@}";;
        * ) echo "> Please answer 'ls' to see the list of files that will be edited, 'no' to cancel, or 'yes' to proceed.";;
    esac
done

while true; do
    echo ""
    read -p "> Step 1: Update SUNDIALS constructor calls with a placeholder for the new SUNContext object. Do this step? [Yy/Nn]" yn
    case $yn in
        [Yy]* ) updateConstructors "${@}"; echo "> SUNDIALS constructors have been updated with a new 'SUNCTX_PLACEHOLDER' argument. 'SUNCTX_PLACEHOLDER' MUST BE REPLACED with a 'SUNContext' object from a call to 'SUNContext_Create'."; break;;
        [Nn]* ) echo "> Step 1 skipped."; break;;
    esac
done

while true; do
    echo ""
    read -p "> Step 2: Update deprecated constants and types to new names. Do this step? [Yy/Nn]" yn
    case $yn in
        [Yy]* ) updateConstantsAndTypes "${@}"; echo "> Constant names and types have been updated. Some manual cleanup may be required."; break;;
        [Nn]* ) echo "> Step 2 skipped."; break;;
    esac
done

echo "> All done."
