# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
#                David J. Gardner @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------
# CMakeLists.txt file for ARKODE serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples
    # tests that are always run
    "ark_analytic.c\;\;"
    # develop tests
    "ark_analytic.c\;arkode.scalar_tolerances 1e-6 1e-8 arkode.table_names ARKODE_ESDIRK547L2SA_7_4_5 ARKODE_ERK_NONE\;"
    "ark_advection_diffusion_reaction_splitting.c\;\;develop"
    "ark_analytic_lsrk.c\;\;develop"
    "ark_analytic_lsrk_varjac.c\;\;develop"
    "ark_analytic_lsrk_domeigest.c\;\;develop"
    "ark_analytic_lsrk_domeigest.c\;arkid.dom_eig_est_init_preprocess_iters 1 sundomeigestimator.max_iters 1\;develop"
    "ark_analytic_mels.c\;\;develop"
    "ark_analytic_nonlin.c\;\;develop"
    "ark_analytic_partitioned.c\;forcing\;develop"
    "ark_analytic_partitioned.c\;splitting\;develop"
    "ark_analytic_partitioned.c\;splitting ARKODE_SPLITTING_BEST_2_2_2\;develop"
    "ark_analytic_partitioned.c\;splitting ARKODE_SPLITTING_RUTH_3_3_2\;exclude-single"
    "ark_analytic_partitioned.c\;splitting ARKODE_SPLITTING_YOSHIDA_8_6_2\;exclude-single"
    "ark_analytic_ssprk.c\;\;develop"
    "ark_brusselator_1D_mri.c\;\;develop"
    "ark_brusselator_fp.c\;\;exclude-single"
    "ark_brusselator_mri.c\;\;develop"
    "ark_brusselator.c\;\;develop"
    "ark_brusselator1D_imexmri.c\;0 0.001\;exclude-single"
    "ark_brusselator1D_imexmri.c\;2 0.001\;exclude-single"
    "ark_brusselator1D_imexmri.c\;3 0.001\;exclude-single"
    "ark_brusselator1D_imexmri.c\;4 0.001\;exclude-single"
    "ark_brusselator1D_imexmri.c\;5 0.001\;exclude-single"
    "ark_brusselator1D_imexmri.c\;6 0.001\;exclude-single"
    "ark_brusselator1D_imexmri.c\;7 0.001\;exclude-single"
    "ark_brusselator1D.c\;\;exclude-single"
    "ark_conserved_exp_entropy_ark.c\;1 0\;develop"
    "ark_conserved_exp_entropy_ark.c\;1 1\;develop"
    "ark_conserved_exp_entropy_erk.c\;1\;develop"
    "ark_damped_harmonic_symplectic.c\;\;exclude-single"
    "ark_dissipated_exp_entropy.c\;1 0\;develop"
    "ark_dissipated_exp_entropy.c\;1 1\;develop"
    "ark_harmonic_symplectic.c\;\;exclude-single"
    "ark_heat1D_adapt.c\;\;develop"
    "ark_heat1D.c\;\;develop"
    "ark_kepler.c\;--stepper ERK --step-mode adapt\;develop"
    "ark_kepler.c\;--stepper ERK --step-mode fixed --count-orbits\;develop"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --count-orbits --use-compensated-sums\;develop"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_EULER_1_1 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_LEAPFROG_2_2 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_2_2 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_3_3 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_4_4 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_5_6 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_PSEUDO_LEAPFROG_2_2 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_RUTH_3_3 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;--stepper SPRK --step-mode fixed --method ARKODE_SPRK_YOSHIDA_6_8 --tf 50 --check-order --nout 1\;exclude-single"
    "ark_kepler.c\;\;develop"
    "ark_kpr_mri.c\;0 1 0.005\;develop"
    "ark_kpr_mri.c\;1 0 0.01\;develop"
    "ark_kpr_mri.c\;1 1 0.002\;exclude-single"
    "ark_kpr_mri.c\;2 4 0.002\;develop"
    "ark_kpr_mri.c\;3 2 0.001\;develop"
    "ark_kpr_mri.c\;4 3 0.001\;develop"
    "ark_kpr_mri.c\;5 4 0.001\;develop"
    "ark_kpr_mri.c\;6 5 0.001\;develop"
    "ark_kpr_mri.c\;7 2 0.002\;develop"
    "ark_kpr_mri.c\;8 3 0.001 -100 100 0.5 1\;exclude-single"
    "ark_kpr_mri.c\;9 3 0.001 -100 100 0.5 1\;exclude-single"
    "ark_kpr_mri.c\;10 4 0.001 -100 100 0.5 1\;exclude-single"
    "ark_kpr_mri.c\;11 2 0.001\;develop"
    "ark_kpr_mri.c\;12 3 0.005\;develop"
    "ark_kpr_mri.c\;13 4 0.01\;exclude-single"
    "ark_KrylovDemo_prec.c\;\;exclude-single"
    "ark_KrylovDemo_prec.c\;1\;exclude-single"
    "ark_KrylovDemo_prec.c\;2\;exclude-single"
    "ark_lotka_volterra_ASA.c\;--check-freq 1\;develop"
    "ark_lotka_volterra_ASA.c\;--check-freq 5\;develop"
    "ark_onewaycouple_mri.c\;\;develop"
    "ark_reaction_diffusion_mri.c\;\;develop"
    "ark_robertson_constraints.c\;\;exclude-single"
    "ark_robertson_root.c\;\;exclude-single"
    "ark_robertson.c\;\;exclude-single"
    "ark_twowaycouple_mri.c\;\;develop")

if(SUNDIALS_BUILD_WITH_MONITORING)
  list(APPEND ARKODE_examples "ark_brusselator_fp.c\;1\;exclude-single")
endif()

# Auxiliary files to install
set(ARKODE_extras
    ark_analytic_nonlin_stats.csv
    ark_damped_harmonic_symplectic.h
    ark_harmonic_symplectic.h
    ark_kepler_plot.py
    ark_kepler.h
    ark_reaction_diffusion_mri_fast_stats.csv
    ark_reaction_diffusion_mri_slow_stats.csv
    ark_robertson_stats.csv
    plot_brusselator1D.py
    plot_heat1D_adapt.py
    plot_heat1D.py
    plot_sol_log.py
    plot_sol.py)

# Add the build and install targets for each ARKODE example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  if(NOT TARGET ${example_target})
    # example source files
    sundials_add_executable(${example_target} ${example})

    # folder for IDEs
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(
      ${example_target} sundials_arkode sundials_nvecserial
      sundials_sundomeigestpower ${EXE_EXTRA_LINK_LIBS} ${DOMEIG_LIBRARIES})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example_target})
  else()
    string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example_target}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # install example source and out files
  if(EXAMPLES_INSTALL)
    install(FILES ${example} ${test_name}.out
            DESTINATION ${EXAMPLES_INSTALL_PATH}/arkode/C_serial)
  endif()

endforeach(example_tuple ${ARKODE_examples})

if(EXAMPLES_INSTALL)
  sundials_install_examples(
    arkode ARKODE_examples
    CMAKE_TEMPLATE cmakelists_C_ex.in
    MAKE_TEMPLATE makefile_serial_C_ex.in
    SUNDIALS_TARGETS nvecserial arkode sundomeigestpower
    DESTINATION arkode/C_serial
    EXTRA_FILES README ${ARKODE_extras}
    TEST_INSTALL C_serial)
endif()
