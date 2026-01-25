#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): Radu Serban, David J. Gardner, Cody J. Balos @ LLNL
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
# Script to build SUNDIALS tar-files.
# ------------------------------------------------------------------------------

# Exit immediately if a command fails
set -e
set -o pipefail

#==============================================================================
#
# Prints usage if help was requested or if the script was incorrectly called
#
#==============================================================================

function print_usage
{
    # Location of tarballs
    cd ..
    location=`pwd`
    cd -

    # Print help message
    echo ""
    echo "Usage: tarscript.sh [-hsv] [module]"
    echo "   -h       : help"
    echo "   -s       : short (no documentation)"
    echo "   -q       : quiet"
    echo "   module   : all, sundials, arkode, cvode, cvodes, ida, idas, kinsol"
    echo ""
    echo "Notes: If the module is not specified, all tarballs are created."
    echo "       This script must be executed from within its directory."
    echo "       Tarballs will be created in $location/tarballs"
    echo ""
    exit 1
}

#==============================================================================
#
# MAIN SCRIPT
#
#==============================================================================

#---------------------------------------------------------
# VERSION NUMBERS
#---------------------------------------------------------

SUN_VER="7.6.0"
CV_VER="7.6.0"
CVS_VER="7.6.0"
IDA_VER="7.6.0"
IDAS_VER="6.6.0"
KIN_VER="7.6.0"
ARK_VER="6.6.0"

#---------------------------------------------------------
# Test if the script is executed from within its directory
#---------------------------------------------------------
scriptbase=`basename $0`
if [ ! -f "$scriptbase" ] ; then
    print_usage
fi

#-------------
# Define flags
#-------------
err=F   # Error flag
hlp=F   # Help flag
doc=T   # Documentation flag
quiet=F # Quiet flag

#----------------
# Process options
#----------------
while getopts ":hsq" name ; do
    case $name in
        h) hlp=T;;
        s) doc=F;;
        q) quiet=T;;
        ?) echo "Invalid option"; err=T;;
    esac
done

#------------------------------
# Extract argument (module name)
#------------------------------
shift $(($OPTIND - 1))
module=$1

#------------
# Test module
#------------
do_sundials=F
do_arkode=F
do_cvode=F
do_cvodes=F
do_ida=F
do_idas=F
do_kinsol=F

if [ -z $module ]; then
    module="all"
    do_sundials=T
    do_arkode=T
    do_cvode=T
    do_cvodes=T
    do_ida=T
    do_idas=T
    do_kinsol=T
else
    case $module in
        all)
            do_sundials=T
            do_arkode=T
            do_cvode=T
            do_cvodes=T
            do_ida=T
            do_idas=T
            do_kinsol=T
            ;;
        sundials)
            do_sundials=T
            ;;
        arkode)
            do_arkode=T
            ;;
        cvode)
            do_cvode=T
            ;;
        cvodes)
            do_cvodes=T
            ;;
        ida)
            do_ida=T
            ;;
        idas)
            do_idas=T
            ;;
        kinsol)
            do_kinsol=T
            ;;
        *)
            echo "Invalid module $module"
            err=T
    esac
fi

#------------
# Print usage
#------------
if [ $err = "T" -o $hlp = "T" ]; then
    print_usage
fi

#------------
# tar command
#------------
if [ "$quiet" = "T" ]; then
    tar="tar -uhf"
else
    tar="tar -uvhf"
fi

#------------------
# Define some names
#------------------

# Location of the tarscript
scriptdir=`pwd`

# Sundials directory
cd ..
sundialsdir=`pwd`

# Location of tarballs
\rm -rf tarballs
mkdir tarballs
cd tarballs
rootdir=`pwd`

# Location of source for tarballs
tmpdir=$rootdir/tmp_dir

#---------------------------
# Create temporary directory
#---------------------------
echo -e "\n--- Create temporary directory ---"
echo "copy $sundialsdir to $rootdir/$tmpdir"
rm -rf $tmpdir
mkdir $tmpdir
mkdir $tmpdir/benchmarks
mkdir $tmpdir/cmake
mkdir $tmpdir/doc
mkdir $tmpdir/examples
mkdir $tmpdir/include
mkdir $tmpdir/src
mkdir $tmpdir/test
mkdir $tmpdir/tools

#----------------------------------
# Copy appropriate files in $tmpdir
# and create additional files as needed
#----------------------------------

cp $sundialsdir/CHANGELOG.md $tmpdir/
cp $sundialsdir/CITATIONS.md $tmpdir/
cp $sundialsdir/CMakeLists.txt $tmpdir/
cp $sundialsdir/CONTRIBUTING.md $tmpdir/
cp $sundialsdir/LICENSE $tmpdir/
cp $sundialsdir/NOTICE $tmpdir/
cp $sundialsdir/README.md $tmpdir/

cp -r $sundialsdir/benchmarks $tmpdir/
cp -r $sundialsdir/cmake $tmpdir/

cp    $sundialsdir/examples/CMakeLists.txt $tmpdir/examples/
cp -r $sundialsdir/examples/utilities $tmpdir/examples/
cp -r $sundialsdir/examples/templates $tmpdir/examples/

cp -r $sundialsdir/external $tmpdir/

cp -r $sundialsdir/include/sundials $tmpdir/include/
cp -r $sundialsdir/include/nvector $tmpdir/include/
cp -r $sundialsdir/include/sunmatrix $tmpdir/include/
cp -r $sundialsdir/include/sunmemory $tmpdir/include/
cp -r $sundialsdir/include/sunlinsol $tmpdir/include/
cp -r $sundialsdir/include/sunnonlinsol $tmpdir/include/
cp -r $sundialsdir/include/sunadaptcontroller $tmpdir/include/
cp -r $sundialsdir/include/sunadjointcheckpointscheme $tmpdir/include/
cp -r $sundialsdir/include/sundomeigest $tmpdir/include/

cp    $sundialsdir/src/CMakeLists.txt $tmpdir/src/
cp -r $sundialsdir/src/sundials $tmpdir/src/
cp -r $sundialsdir/src/nvector $tmpdir/src/
cp -r $sundialsdir/src/sunmatrix $tmpdir/src/
cp -r $sundialsdir/src/sunmemory $tmpdir/src/
cp -r $sundialsdir/src/sunlinsol $tmpdir/src/
cp -r $sundialsdir/src/sunnonlinsol $tmpdir/src/
cp -r $sundialsdir/src/sunadaptcontroller $tmpdir/src/
cp -r $sundialsdir/src/sunadjointcheckpointscheme $tmpdir/src/
cp -r $sundialsdir/src/sundomeigest $tmpdir/src/

cp    $sundialsdir/test/testRunner $tmpdir/test/
cp -r $sundialsdir/test/unit_tests $tmpdir/test/

cp -r $sundialsdir/tools $tmpdir/

# Clean up tmpdir
rm -rf $tmpdir/doc/shared/__pycache__
find $tmpdir -name "*~" -delete
find $tmpdir -name "*.orig" -delete
find $tmpdir -name "#*#" -delete
find $tmpdir -name ".DS_Store" -delete

# Remove ignored or untracked files that may have been added
cd $sundialsdir
for f in $(git ls-files --others --directory benchmarks cmake doc examples external include src test tools); do
    rm -rf $tmpdir/$f
done
cd -

# Make superbuild for cross-package references
if [ $doc = "T" ]; then
    echo -e "--- SUNDIALS Superbuild documentation"
    cd $sundialsdir/doc/superbuild
    make clean
    make html
    cd -
fi

# Make install guide and copy to $tmpdir
if [ $doc = "T" ]; then
    echo -e "--- SUNDIALS Common Installation documentation"
    cd $sundialsdir/doc/install_guide
    make clean
    make latexpdf
    cp build/latex/INSTALL_GUIDE.pdf $tmpdir/INSTALL_GUIDE.pdf
    cd -
fi

declare -a packages_rst=('arkode' 'kinsol')
for pkg in "${packages_rst[@]}";
do
    do_package=do_${pkg}
    if [ $do_sundials = "T" -o ${!do_package} = "T" ]; then
        cp -r $sundialsdir/include/$pkg $tmpdir/include/
        cp -r $sundialsdir/src/$pkg $tmpdir/src/
        cp -r $sundialsdir/examples/$pkg $tmpdir/examples/
        mkdir -p $tmpdir/doc/$pkg
        if [ $doc = "T" ]; then
            echo -e "--- ${pkg} documentation"
            cd $sundialsdir/doc/$pkg/guide
            make clean
            make latexpdf
            cp build/latex/*_guide.pdf $tmpdir/doc/$pkg/
            cd -
            cd $sundialsdir/doc/$pkg/examples
            make clean
            make latexpdf
            cp build/latex/*_examples.pdf $tmpdir/doc/$pkg/
            cd -
        fi
    fi
done

declare -a packages=('cvode' 'cvodes' 'ida' 'idas')
for pkg in "${packages[@]}";
do
    do_package=do_${pkg}
    if [ $do_sundials = "T" -o ${!do_package} = "T" ]; then
        cp -r $sundialsdir/include/$pkg $tmpdir/include/
        cp -r $sundialsdir/src/$pkg $tmpdir/src/
        cp -r $sundialsdir/examples/$pkg $tmpdir/examples/
        mkdir -p $tmpdir/doc/$pkg
        if [ $doc = "T" ]; then
            echo -e "--- ${pkg} documentation"
            cd $sundialsdir/doc/$pkg/guide
            make clean
            make latexpdf
            cp build/latex/*_guide.pdf $tmpdir/doc/$pkg/
            cd -
            cd $sundialsdir/doc/$pkg
            make QUIET= ex_pdf
            cp *_examples.pdf $tmpdir/doc/$pkg/
            cd -
        fi
    fi
done

#---------------------------
# Create tar files
#---------------------------

# NOTE: The "shared" script must be called first for each case below as it
# creates the initial archive appended to by the package tarscripts

# Initial name of the directory to be archived
distrobase="$tmpdir"

# SUNDIALS
if [ $do_sundials = "T" ]; then
    echo -e "\n--- Generate SUNDIALS tarball ---"

    mv $distrobase sundials-$SUN_VER
    distrobase="sundials-"$SUN_VER
    filename="sundials-"$SUN_VER

    tarfile=$filename".tar"
    $scriptdir/shared.sh $tarfile $distrobase $doc $tar
    $scriptdir/arkode.sh $tarfile $distrobase $doc $tar
    $scriptdir/cvode.sh  $tarfile $distrobase $doc $tar
    $scriptdir/cvodes.sh $tarfile $distrobase $doc $tar
    $scriptdir/ida.sh    $tarfile $distrobase $doc $tar
    $scriptdir/idas.sh   $tarfile $distrobase $doc $tar
    $scriptdir/kinsol.sh $tarfile $distrobase $doc $tar

    ### Don't release MATLAB until brought current with new version(s)
    # $scriptdir/stb $tarfile $distrobase $doc $tar

    gzip $tarfile
fi

# ARKODE
if [ $do_arkode = "T" ]; then
    echo -e "\n--- Generate ARKODE tarball ---"

    mv $distrobase arkode-$ARK_VER
    distrobase="arkode-"$ARK_VER
    filename="arkode-"$ARK_VER

    tarfile=$filename".tar"
    $scriptdir/shared.sh $tarfile $distrobase $doc "F" $tar
    $scriptdir/arkode.sh $tarfile $distrobase $doc $tar
    gzip $tarfile
fi

# CVODE
if [ $do_cvode = "T" ]; then
    echo -e "\n--- Generate CVODE tarball ---"

    mv $distrobase cvode-$CV_VER
    distrobase="cvode-"$CV_VER
    filename="cvode-"$CV_VER

    tarfile=$filename".tar"
    $scriptdir/shared.sh $tarfile $distrobase $doc "F" $tar
    $scriptdir/cvode.sh  $tarfile $distrobase $doc $tar
    gzip $tarfile
fi

# CVODES
if [ $do_cvodes = "T" ]; then
    echo -e "\n--- Generate CVODES tarball ---"

    mv $distrobase cvodes-$CVS_VER
    distrobase="cvodes-"$CVS_VER
    filename="cvodes-"$CVS_VER

    tarfile=$filename".tar"
    $scriptdir/shared.sh $tarfile $distrobase $doc "F" $tar
    $scriptdir/cvodes.sh $tarfile $distrobase $doc $tar
    gzip $tarfile
fi

# IDA
if [ $do_ida = "T" ]; then
    echo -e "\n--- Generate IDA tarball ---"

    mv $distrobase ida-$IDA_VER
    distrobase="ida-"$IDA_VER
    filename="ida-"$IDA_VER

    tarfile=$filename".tar"
    $scriptdir/shared.sh $tarfile $distrobase $doc "F" $tar
    $scriptdir/ida.sh    $tarfile $distrobase $doc $tar
    gzip $tarfile
fi

# IDAS
if [ $do_idas = "T" ]; then
    echo -e "\n--- Generate IDAS tarball ---"

    mv $distrobase idas-$IDAS_VER
    distrobase="idas-"$IDAS_VER
    filename="idas-"$IDAS_VER

    tarfile=$filename".tar"
    $scriptdir/shared.sh $tarfile $distrobase $doc "F" $tar
    $scriptdir/idas.sh   $tarfile $distrobase $doc $tar
    gzip $tarfile
fi

# KINSOL
if [ $do_kinsol = "T" ]; then
    echo -e "\n--- Generate KINSOL tarball ---"

    mv $distrobase kinsol-$KIN_VER
    distrobase="kinsol-"$KIN_VER
    filename="kinsol-"$KIN_VER

    tarfile=$filename".tar"
    $scriptdir/shared.sh $tarfile $distrobase $doc "F" $tar
    $scriptdir/kinsol.sh $tarfile $distrobase $doc $tar
    gzip $tarfile
fi

#---------------------------
# Remove temporary directory
#---------------------------

rm -rf $distrobase

#du -h *.tar.gz

exit 0

# That's all folks...
