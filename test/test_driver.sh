#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# SUNDIALS testing driver script
# ------------------------------------------------------------------------------

# Exit when any command fails
set -e

# Help message
help ()
{
    cat <<EOF

    NAME

        $0 - SUNDIALS testing driver

    USAGE

        $0 [OPTIONS] [ENV SETUP OPTIONS]

    OPTIONS

        --help, -h
            Display this message and exit.

        --testtype TYPE
            Run a predefined set of tests or specify the test configuration
            using additional input options. TYPE must be one of:

            custom  -- (default) single user defined test configuration
            branch  -- quick tests with a few configurations
            pr      -- create sundials tarball and test more configurations
            release -- create sundials tarball and test even more configurations

    When running a custom test the following options may be provided to specify
    the test configuration. If any of these options are provided with a test
    type other than "custom," their values will be ignored.

        --tarball PACKAGE
            Create a SUNDIALS release tarball and run tests using the code
            extracted from the tarball. PACKAGE must be one of:

            NONE     -- (default) do not create a tarball
            sundials -- create tarball containing all SUNDIALS packages together
            arkode   -- create tarball containing ARKODE only
            cvode    -- create tarball containing CVODE only
            cvodes   -- create tarball containing CVODES only
            ida      -- create tarball containing IDA only
            idas     -- create tarball containing IDAS only
            kinsol   -- create tarball containing KINSOL only
            all      -- create all possible tarballs

        --sunrealtype TYPE
            Precision to use in a custom test. TYPE must be one of:

            double   -- (default) use double precision
            single   -- use single precision
            extended -- use extended precision

        --indexsize SIZE
            Index size to use in a custom test. SIZE must be one of:

            64 -- (default) use 64-bit indexing
            32 -- use 32-bit indexing

        --libtype TYPE
            Library type to create in a custom test. TYPE must be one of:

            static -- build static libraries
            shared -- build shared libraries
            both   -- build static and shared libraries

        --tpls ON/OFF
            Enable or disable external third-party libraries in a custom test.

        --suntesttype TYPE
            SUNDIALS test type for a custom test. TYPE must be one of:

            STD -- run standard tests
            DEV -- run development tests

    With any test type the following options may be set to alter the test setup
    or behavior.

        --env ENVFILE
            Environment configuration file to use for testing.

        --buildjobs NUM
            Set the number of jobs for parallel builds (default 1).

        --testjobs NUM
            Set the number of jobs for parallel testing (default 1). Note this
            is the number of jobs CTest will use to run tests simultaneously and
            does not impact the number of MPI tasks used in MPI-parallel tests.

        --phase PHASE
            Testing phase to stop after. PHASE must be one of:

            ENV              -- setup the environment
            CONFIG           -- configuring with CMake
            BUILD            -- build with make
            TEST             -- test with make test
            INSTALL          -- install with make install
            TEST_INSTALL     -- test with make test_install
            TEST_INSTALL_ALL -- test with make test_install_all

    ENV SETUP OPTIONS

        Any additional inputs provided will be forwarded to the environment
        script called before each test.

    Example usage:

        $0
        $0 --testtype release --buildjobs 4
        $0 --phase CONFIG --indexsize 32 --tpls ON --env env/my_env.sh

EOF
}

# Print input args
echo "./test_driver.sh $*" | tee -a suntest.log

# --------------
# Default values
# --------------

buildjobs=0
testjobs=0
testtype="CUSTOM"
tarball="NONE"
sunrealtype="double"
indexsize="64"
libtype="both"
tpls="OFF"
suntesttype="DEV"
phase=""

# ------------
# Parse inputs
# ------------

# Store extra inputs
EXTRA_ARGS=()

# Iterate over inputs
while [[ $# -gt 0 ]]; do

    input="$1"

    case $input in

        --help|-h)
            help
            exit 1;;

        --testtype)
            testtype=$2
            case "$testtype" in
                BRANCH|branch|Branch)
                    testtype=BRANCH
                    ;;
                PR|pr)
                    testtype=PR
                    ;;
                RELEASE|release|Release)
                    testtype=RELEASE
                    ;;
                CUSTOM|custom|Custom)
                    testtype=CUSTOM
                    ;;
                *)
                    echo "ERROR: Invalid test type $testtype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --buildjobs)
            buildjobs=$2
            if [[ $buildjobs -lt 1 ]]; then
                echo "ERROR: Build jobs < 1"
                help
                exit 1
            fi
            shift 2;;

        --testjobs)
            testjobs=$2
            if [[ $testjobs -lt 1 ]]; then
                echo "ERROR: Test jobs < 1"
                help
                exit 1
            fi
            shift 2;;

        --tarball)
            tarball=$2
            case "$tarball" in
                sundials|arkode|cvode|cvodes|ida|idas|kinsol|all)
                ;;
                *)
                    echo "ERROR: Invalid tarball option $tarball"
                    help
                    exit 1;;
            esac
            shift 2;;

        --sunrealtype)
            sunrealtype=$2
            case "$sunrealtype" in
                SINGLE|Single|single)
                    sunrealtype=single
                    ;;
                DOUBLE|Double|double)
                    sunrealtype=double
                    ;;
                EXTENDED|Extended|extended)
                    sunrealtype=extended
                    ;;
                *)
                    echo "ERROR: Invalid real type option $sunrealtype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --indexsize)
            indexsize=$2
            case "$indexsize" in
                32|64)
                ;;
                *)
                    echo "ERROR: Invalid index size option $indexsize"
                    help
                    exit 1;;
            esac
            shift 2;;

        --libtype)
            libtype=$2
            case "$libtype" in
                STATIC|Static|static)
                    libtype=static
                    ;;
                SHARED|Shared|shared)
                    libtype=shared
                    ;;
                BOTH|Both|both)
                    libtype=both
                    ;;
                *)
                    echo "ERROR: Invalid library type option $libtype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --tpls)
            tpls=$2
            case "$tpls" in
                ON|On|on)
                    tpls=ON
                    ;;
                OFF|Off|off)
                    tpls=OFF
                    ;;
                *)
                    echo "ERROR: Invalid tpl option $tpl"
                    help
                    exit 1;;
            esac
            shift 2;;

        --suntesttype)
            suntesttype=$2
            case "$suntesttype" in
                DEV|Dev|dev)
                    suntesttype=DEV
                    ;;
                STD|Std|std)
                    suntesttype=STD
                    ;;
                NONE|None|none)
                    suntesttype=NONE
                    ;;
                *)
                    echo "ERROR: Invalid SUNDIALS test type option $suntesttype"
                    help
                    exit 1;;
            esac
            shift 2;;

        --phase)
            phase=$2
            case "$phase" in
                ENV|Env|env)
                    phase=ENV
                    ;;
                CONFIG|Config|config)
                    phase=CONFIG
                    ;;
                BUILD|Build|build)
                    phase=BUILD
                    ;;
                TEST|Test)
                    phase=TEST
                    ;;
                INSTALL|Install|install)
                    phase=INSTALL
                    ;;
                TEST_INSTALL|Test_Install|test_install)
                    phase=TEST_INSTALL
                    ;;
                TEST_INSTALL_ALL|Test_Install_All|test_install_all)
                    phase=TEST_INSTALL_ALL
                    ;;
                *)
                    echo "ERROR: Invalid phase option $phase"
                    help
                    exit 1;;
            esac
            shift 2;;

        --env)
            envfile=$2
            if [[ ! -f "$envfile" ]]; then
                echo "ERROR: Environment file does not exist $envfile"
                help
                exit 1
            else
                export SUNDIALS_ENV_FILE="${envfile}"
            fi
            shift 2;;

        *)
            EXTRA_ARGS+=("$input")
            shift;;
    esac
done

# Reset positional inputs to extra inputs
set -- "${EXTRA_ARGS[@]}"

# Initialize failure counter (0 = success)
passfail=0

# Location of testing directory
testroot=$(pwd)

# ------------------------------------------------------------------------------
# Print test header
# ------------------------------------------------------------------------------

echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS test" | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log
git log -1 --pretty='Commit: %H%nAuthor: %an <%ae>%nDate:   %ad%n%n%s' | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# ------------------------------------------------------------------------------
# Create an array of test configurations
# ------------------------------------------------------------------------------

args_realtypes=()
args_indexsizes=()
args_libtypes=()
args_tpls=()
args_suntests=()
args_phase=()

case "$testtype" in

    BRANCH)
        # Don't create tarballs
        tarball=NONE

        # Test configs
        for is in 32 64; do
            args_realtypes+=("double")
            args_indexsizes+=("${is}")
            args_libtypes+=("static")
            args_tpls+=("ON")
            args_suntests+=("DEV")
            args_phase+=("TEST")
        done
        ;;

    PR)
        # Create sundials tarball
        tarball=sundials

        # Test configs
        for is in 32 64; do
            args_realtypes+=("double")
            args_indexsizes+=("${is}")
            args_libtypes+=("static")
            args_tpls+=("ON")
            args_suntests+=("DEV")
            args_phase+=("")
        done
        ;;

    RELEASE)
        # Create sundials tarball
        tarball=sundials

        # Address sanitizer tests (TPLs OFF)
        for is in 32 64; do
            args_realtypes+=("double")
            args_indexsizes+=("${is}")
            args_libtypes+=("static")
            args_tpls+=("OFF")
            args_suntests+=("DEV")
            args_phase+=("TEST")
        done

        # Test configs
        for rt in single double extended; do
            for is in 32 64; do
                for lt in static shared; do
                    args_realtypes+=("${rt}")
                    args_indexsizes+=("${is}")
                    args_libtypes+=("${lt}")
                    args_tpls+=("ON")
                    # Development test output files created with double
                    if [[ "${rt}" == "double" ]]; then
                        args_suntests+=("DEV")
                    else
                        args_suntests+=("STD")
                    fi
                    args_phase+=("")
                done
            done
        done
        ;;

    CUSTOM)
        # Use default or user defined values
        args_realtypes+=("${sunrealtype}")
        args_indexsizes+=("${indexsize}")
        args_libtypes+=("${libtype}")
        args_tpls+=("${tpls}")
        args_suntests+=("${suntesttype}")
        args_phase+=("${phase}")
        if [ "${sunrealtype}" != "double" ] && [ "${suntesttype}" == "DEV" ]; then
            echo "WARNING: DEV tests may fail with ${sunrealtype} precision"
        fi
        ;;

    *)
        echo "ERROR: Unknown test type"
        exit 1
        ;;
esac

# Number of configurations to test
nconfigs=${#args_suntests[@]}

# ------------------------------------------------------------------------------
# Create and extract tarballs (if necessary)
# ------------------------------------------------------------------------------

# Directories to run tests in (one for each tarball created)
testdir=()

if [ "$tarball" != NONE ]; then

    # Setup the environment with the first test configuration. This is really in
    # case the environment needs to be setup to build the documentation in the
    # tarscript e.g., activates a Python virtual environment.
    env_config=("${args_realtypes[0]}"
                "${args_indexsizes[0]}"
                "${args_libtypes[0]}"
                "${args_tpls[0]}"
                "${args_suntests[0]}")

    # Setup environment
    time source env/setup_env.sh "${env_config[@]}" "${EXTRA_ARGS[@]}"

    rc=${PIPESTATUS[0]}
    echo -e "\nsetup_env.sh returned $rc\n" | tee -a setup_env.log
    if [ "$rc" -ne 0 ]; then exit 1; fi

    # Remove old tarballs directory
    \rm -rf "$testroot/tarballs"

    # Create directory for tarballs
    mkdir "$testroot/tarballs"

    # Run tarscript to create tarballs
    cd ..
    cd scripts

    echo "START TARSCRIPT"
    time ./tarscript.sh $tarball | tee -a tar.log

    # Check tarscript return code
    rc=${PIPESTATUS[0]}
    echo -e "\ntarscript.sh returned $rc\n" | tee -a tar.log
    if [ "$rc" -ne 0 ]; then exit 1; fi

    # Relocate log and tarballs
    mv tar.log "$testroot/tarballs/."
    mv "$testroot/setup_env.log" "$testroot/tarballs/."
    mv ../tarballs/* "$testroot/tarballs/."

    # Move to tarball directory
    cd "$testroot/tarballs"

    # Loop over tarballs and test each one
    for tb in *.tar.gz; do

        # Get package name
        package=${tb%.tar.gz}

        echo "START UNTAR"
        tar -xvzf "$tb" 2>&1 | tee -a tar.log

        # Check tar return code
        rc=${PIPESTATUS[0]}
        echo -e "\ntar -xzvf returned $rc\n" | tee -a tar.log
        if [ "$rc" -ne 0 ]; then exit 1; fi

        # Copy environment and testing scripts from original test directory
        cp -r "$testroot/env" "$package/test/."
        cp "$testroot/config_cmake.py" "$package/test/."

        testdir+=("$testroot/tarballs/$package/test")
    done

else

    # No tarballs created, use the original test directory
    testdir+=("$testroot")

fi

# Number of test directories
ntestdirs=${#testdir[@]}

# ------------------------------------------------------------------------------
# Run tests
# ------------------------------------------------------------------------------

for ((j=0;j<ntestdirs;j++)); do

    cd "${testdir[j]}"

    for ((i=0;i<nconfigs;i++)); do

        # ---------------------
        # Setup the environment
        # ---------------------

        env_config=("${args_realtypes[i]}"
                    "${args_indexsizes[i]}"
                    "${args_libtypes[i]}"
                    "${args_tpls[i]}"
                    "${args_suntests[i]}")

        # Print test header for Jenkins section collapsing
        echo "TEST: ${env_config[*]}"

        # Setup environment
        time source env/setup_env.sh "${env_config[@]}" "${EXTRA_ARGS[@]}"

        rc=${PIPESTATUS[0]}
        echo -e "\nsetup_env.sh returned $rc\n" | tee -a setup_env.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if the environment sets the number of build and test jobs but do
        # not override the command line options. If neither are set then set the
        # default value.
        if [[ $buildjobs -lt 1 ]]; then
            if [[ -n "${SUNDIALS_BUILD_JOBS}" ]]; then
                # non-empty contiguous string of digits
                re='^[0-9]+$'
                if [[ "${SUNDIALS_BUILD_JOBS}" =~ $re ]] ; then
                    buildjobs="${SUNDIALS_BUILD_JOBS}"
                else
                    buildjobs=1
                fi
            else
                buildjobs=1
            fi
        fi

        if [[ $testjobs -lt 1 ]]; then
            if [[ -n "${SUNDIALS_TEST_JOBS}" ]]; then
                # non-empty contiguous string of digits
                re='^[0-9]+$'
                if [[ "${SUNDIALS_TEST_JOBS}" =~ $re ]] ; then
                    testjobs="${SUNDIALS_TEST_JOBS}"
                else
                    testjobs=1
                fi
            else
                testjobs=1
            fi
        fi

        # Check if this is the last phase
        if [ "${args_phase[i]}" == "ENV" ]; then
            echo "PASSED: ${env_config[*]}"
            continue
        fi

        # -----------------------
        # Create test directories
        # -----------------------

        tmp="${env_config[*]}"
        tmp=${tmp//" "/"_"}

        builddir="build_${tmp}"
        installdir="$(pwd)/install_${tmp}"

        # Add host name to directory names
        if [ -n "$HOST" ]; then
            builddir="${builddir}_${HOST}"
            installdir="${installdir}_${HOST}"
        elif [ -n "$HOSTNAME" ]; then
            builddir="${builddir}_${HOSTNAME}"
            installdir="${installdir}_${HOSTNAME}"
        fi

        # Remove old build and install directories
        \rm -rf "$builddir" "$installdir"

        # Create and move to new build directory, move configure log
        mkdir "$builddir"
        mv setup_env.log "$builddir/."
        cd "$builddir"

        # -----------------------
        # Create CMake cache file
        # -----------------------

        time python3 ../config_cmake.py \
            --readenv \
            --install-prefix "${installdir}" | tee -a configure.log

        rc=${PIPESTATUS[0]}
        echo -e "\nconfig_cmake.py returned $rc\n" | tee -a configure.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # ---------
        # Configure
        # ---------

        echo "START CMAKE"

        time cmake -C sundials.cmake ../../. | tee -a configure.log

        rc=${PIPESTATUS[0]}
        echo -e "\ncmake returned $rc\n" | tee -a configure.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${args_phase[i]}" == "CONFIG" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # -----
        # Build
        # -----

        echo "START MAKE"

        time make -j "$buildjobs" 2>&1 | tee make.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake returned $rc\n" | tee -a make.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${args_phase[i]}" == "BUILD" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # ----
        # Test
        # ----

        echo "START TEST"

        time ctest --output-on-failure -j "$testjobs" test 2>&1 | tee test.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake test returned $rc\n" | tee -a test.log

        # Re-attempt failed tests in case they are unstable
        if [ "$rc" -ne 0 ]; then
            echo "START TEST"
            time ctest --rerun-failed --output-on-failure -j "$testjobs" test 2>&1 | tee -a test.log
            rc=${PIPESTATUS[0]}
            echo -e "\nmake test returned $rc\n" | tee -a test.log
            if [ "$rc" -ne 0 ]; then passfail=1; break; fi
        fi

        # Check if this is the last phase
        if [ "${args_phase[i]}" == "TEST" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # -------
        # Install
        # -------

        echo "START INSTALL"

        time make -j "$buildjobs" install 2>&1 | tee install.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake install returned $rc\n" | tee -a install.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${args_phase[i]}" == "INSTALL" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # ------------
        # Test install
        # ------------

        echo "START TEST_INSTALL"

        time make test_install 2>&1 | tee test_install.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake test_install returned $rc\n" | tee -a test_install.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${args_phase[i]}" == "TEST_INSTALL" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        echo "START TEST_INSTALL_ALL"

        time make test_install_all 2>&1 | tee test_install_all.log

        rc=${PIPESTATUS[0]}
        echo -e "\nmake test_install_all returned $rc\n" | tee -a test_install_all.log
        if [ "$rc" -ne 0 ]; then passfail=1; break; fi

        # Check if this is the last phase
        if [ "${args_phase[i]}" == "TEST_INSTALL_ALL" ]; then
            echo "PASSED: ${env_config[*]}"
            cd ..
            continue
        fi

        # ------------
        # End of tests
        # ------------

        echo "PASSED: ${env_config[*]}" | tee -a suntest.log
        cd ..

    done

    if [ $passfail -ne 0 ]; then
        echo "FAILED: ${env_config[*]}" | tee -a suntest.log
        cd "${testroot}"
        break
    fi

done

# ------------------------------------------------------------------------------
# Print test footer
# ------------------------------------------------------------------------------

if [ $passfail -eq 0 ]; then
    echo "All tests PASSED" | tee -a suntest.log
fi

echo "--------------------------------------------------" | tee -a suntest.log
echo "SUNDIALS test" | tee -a suntest.log
date | tee -a suntest.log
echo "--------------------------------------------------" | tee -a suntest.log

# ------------------------------------------------------------------------------
# Exit
# ------------------------------------------------------------------------------

if [ $passfail -eq 0 ]; then
    exit 0  # Pass
else
    exit 1  # Fail
fi
