#!/bin/bash -e
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# Update SUNDIALS version numbers
# ------------------------------------------------------------------------------

# Set the SUNDIALS major, minor, and patch numbers and the label string. For
# development releases the label string is of the form "-dev.#" and for full
# releases the label string is "".
sun_major=${1:-7}
sun_minor=${2:-5}
sun_patch=${3:-0}
sun_label=${4:-""}
month=${5:-$(date +"%b")}
year=${6:-$(date +"%Y")}

# ------------------------------------------------------------------------------
# Only update the values below if necessary
# ------------------------------------------------------------------------------

date="$month $year"

# Set the SUNDIALS full version number without or with a label
if [ "${sun_label}" == "" ]; then
    sun_ver="${sun_major}.${sun_minor}.${sun_patch}"
else
    sun_ver="${sun_major}.${sun_minor}.${sun_patch}-${sun_label}"
fi

# Set the ARKODE version values. Assume the major version is one less than the
# SUNDIALS major version.
ark_major=$(( sun_major - 1 ))
ark_minor=$sun_minor
ark_patch=$sun_patch
ark_label=$sun_label

if [ "${ark_label}" == "" ]; then
    ark_ver="${ark_major}.${ark_minor}.${ark_patch}"
else
    ark_ver="${ark_major}.${ark_minor}.${ark_patch}-${ark_label}"
fi

# Set the CVODE values. Assume all values are the same as the SUNDIALS values.
cv_major=$sun_major
cv_minor=$sun_minor
cv_patch=$sun_patch
cv_label=$sun_label

if [ "${cv_label}" == "" ]; then
    cv_ver="${cv_major}.${cv_minor}.${cv_patch}"
else
    cv_ver="${cv_major}.${cv_minor}.${cv_patch}-${cv_label}"
fi


# Set the CVODES values. Assume all values are the same as the SUNDIALS values.
cvs_major=$sun_major
cvs_minor=$sun_minor
cvs_patch=$sun_patch
cvs_label=$sun_label

if [ "${cvs_label}" == "" ]; then
    cvs_ver="${cvs_major}.${cvs_minor}.${cvs_patch}"
else
    cvs_ver="${cvs_major}.${cvs_minor}.${cvs_patch}-${cvs_label}"
fi

# Set the IDA values. Assume all values are the same as the SUNDIALS values.
ida_major=$sun_major
ida_minor=$sun_minor
ida_patch=$sun_patch
ida_label=$sun_label

if [ "${ida_label}" == "" ]; then
    ida_ver="${ida_major}.${ida_minor}.${ida_patch}"
else
    ida_ver="${ida_major}.${ida_minor}.${ida_patch}-${ida_label}"
fi

# Set the IDAS version values. Assume the major version is one less than the
# SUNDIALS major version.
idas_major=$(( sun_major - 1 ))
idas_minor=$sun_minor
idas_patch=$sun_patch
idas_label=$sun_label

if [ "${idas_label}" == "" ]; then
    idas_ver="${idas_major}.${idas_minor}.${idas_patch}"
else
    idas_ver="${idas_major}.${idas_minor}.${idas_patch}-${idas_label}"
fi

# Set the KINSOL values. Assume all values are the same as the SUNDIALS values.
kin_major=$sun_major
kin_minor=$sun_minor
kin_patch=$sun_patch
kin_label=$sun_label

if [ "${kin_label}" == "" ]; then
    kin_ver="${kin_major}.${kin_minor}.${kin_patch}"
else
    kin_ver="${kin_major}.${kin_minor}.${kin_patch}-${kin_label}"
fi

# Set the NVector values. Assume all values are two less than the SUNDIALS values.
vec_major=$sun_major
vec_minor=$sun_minor
vec_patch=$sun_patch
vec_label=$sun_label

if [ "${vec_label}" == "" ]; then
    vec_ver="${vec_major}.${vec_minor}.${vec_patch}"
else
    vec_ver="${vec_major}.${vec_minor}.${vec_patch}-${vec_label}"
fi

# Set the SUNMatrix version values. Assume the major version is two less than the
# SUNDIALS major version.
mat_major=$(( sun_major - 2 ))
mat_minor=$sun_minor
mat_patch=$sun_patch
mat_label=$sun_label

if [ "${mat_label}" == "" ]; then
    mat_ver="${mat_major}.${mat_minor}.${mat_patch}"
else
    mat_ver="${mat_major}.${mat_minor}.${mat_patch}-${mat_label}"
fi

# Set the SUNLinearSolver version values. Assume the major version is two less
# than the SUNDIALS major version.
ls_major=$(( sun_major - 2 ))
ls_minor=$sun_minor
ls_patch=$sun_patch
ls_label=$sun_label

if [ "${ls_label}" == "" ]; then
    ls_ver="${ls_major}.${ls_minor}.${ls_patch}"
else
    ls_ver="${ls_major}.${ls_minor}.${ls_patch}-${ls_label}"
fi

# Set the SUNNonlinearSolver version values. Assume the major version is three
# less than the SUNDIALS major version.
nls_major=$(( sun_major - 3 ))
nls_minor=$sun_minor
nls_patch=$sun_patch
nls_label=$sun_label

if [ "${nls_label}" == "" ]; then
    nls_ver="${nls_major}.${nls_minor}.${nls_patch}"
else
    nls_ver="${nls_major}.${nls_minor}.${nls_patch}-${nls_label}"
fi

# Set the SUNDomEigEstimator version values. Assume the major version is six
# less than the SUNDIALS major version.
dee_major=$(( sun_major - 6 ))
dee_minor=$(( sun_minor - 5 )) # TODO(DJG): Will need to remove this at the next major release
dee_patch=$sun_patch
dee_label=$sun_label

if [ "${dee_label}" == "" ]; then
    dee_ver="${dee_major}.${dee_minor}.${dee_patch}"
else
    dee_ver="${dee_major}.${dee_minor}.${dee_patch}-${dee_label}"
fi

# ------------------------------------------------------------------------------
# Wrapper for editing inplace with different sed implementations
# ------------------------------------------------------------------------------

sedi() {
    case $(uname) in
        Darwin*) sedi=('-i' '') ;;
        *) sedi='-i' ;;
    esac
    sed "${sedi[@]}" "$@"
}

# ------------------------------------------------------------------------------
# Update the pyproject.toml file
# ------------------------------------------------------------------------------
fn="../pyproject.toml"
sedi "/^version =/ s/version = \".*\"/version = \"${sun_ver}\"/" $fn

# ------------------------------------------------------------------------------
# Update the main CMakeLists.txt file
# ------------------------------------------------------------------------------

fn="../CMakeLists.txt"
sedi "/set(PACKAGE_VERSION_MAJOR/ s/MAJOR.*/MAJOR \"${sun_major}\")/" $fn
sedi "/set(PACKAGE_VERSION_MINOR/ s/MINOR.*/MINOR \"${sun_minor}\")/" $fn
sedi "/set(PACKAGE_VERSION_PATCH/ s/PATCH.*/PATCH \"${sun_patch}\")/" $fn
sedi "/set(PACKAGE_VERSION_LABEL/ s/LABEL.*/LABEL \"${sun_label}\")/" $fn
sedi "/set(PACKAGE_STRING/        s/STRING.*/STRING \"SUNDIALS ${sun_ver}\")/" $fn

sedi "/arkodelib_VERSION.*/   s/VERSION.*/VERSION \"${ark_ver}\")/" $fn
sedi "/arkodelib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${ark_major}\")/" $fn

sedi "/cvodelib_VERSION.*/   s/VERSION.*/VERSION \"${cv_ver}\")/" $fn
sedi "/cvodelib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${cv_major}\")/" $fn

sedi "/cvodeslib_VERSION.*/   s/VERSION.*/VERSION \"${cvs_ver}\")/" $fn
sedi "/cvodeslib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${cvs_major}\")/" $fn

sedi "/idalib_VERSION.*/   s/VERSION.*/VERSION \"${ida_ver}\")/" $fn
sedi "/idalib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${ida_major}\")/" $fn

sedi "/idaslib_VERSION.*/   s/VERSION.*/VERSION \"${idas_ver}\")/" $fn
sedi "/idaslib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${idas_major}\")/" $fn

sedi "/kinsollib_VERSION.*/   s/VERSION.*/VERSION \"${kin_ver}\")/" $fn
sedi "/kinsollib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${kin_major}\")/" $fn

sedi "/nveclib_VERSION.*/   s/VERSION.*/VERSION \"${vec_ver}\")/" $fn
sedi "/nveclib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${vec_major}\")/" $fn

sedi "/sunmatrixlib_VERSION.*/   s/VERSION.*/VERSION \"${mat_ver}\")/" $fn
sedi "/sunmatrixlib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${mat_major}\")/" $fn

sedi "/sunlinsollib_VERSION.*/   s/VERSION.*/VERSION \"${ls_ver}\")/" $fn
sedi "/sunlinsollib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${ls_major}\")/" $fn

sedi "/sunnonlinsollib_VERSION.*/   s/VERSION.*/VERSION \"${nls_ver}\")/" $fn
sedi "/sunnonlinsollib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${nls_major}\")/" $fn

sedi "/sundomeigestlib_VERSION.*/   s/VERSION.*/VERSION \"${dee_ver}\")/" $fn
sedi "/sundomeigestlib_SOVERSION.*/ s/SOVERSION.*/SOVERSION \"${dee_major}\")/" $fn

# ------------------------------------------------------------------------------
# Update README files
# ------------------------------------------------------------------------------

sedi "s/### Version.*/### Version ${sun_ver} (${date}) ###/" ../README.md

fn="../src/arkode/README.md"
sedi "s/### Version.*/### Version ${ark_ver} (${date})/" $fn
sedi "s/\"User Documentation for ARKODE v.*/\"User Documentation for ARKODE v${ark_ver},\" LLNL technical report/" $fn
sedi "s/LLNL-SM-668082,.*/LLNL-SM-668082, ${date}./" $fn
sedi "s/\"Example Programs for ARKODE v.*/\"Example Programs for ARKODE v${ark_ver},\" Technical Report,/" $fn
sedi "s/Scientific Computation.*/Scientific Computation, ${date}./" $fn

fn="../src/cvode/README.md"
sedi "s/### Version.*/### Version ${cv_ver} (${date})/" $fn
sedi "s/\"User Documentation for CVODE v.*/\"User Documentation for CVODE v${cv_ver},\"/" $fn
sedi "s/UCRL-SM-208108,.*/UCRL-SM-208108, ${date}./" $fn
sedi "s/\"Example Programs for CVODE v.*/\"Example Programs for CVODE v${cv_ver},\"/" $fn
sedi "s/UCRL-SM-208110,.*/UCRL-SM-208110, ${date}./" $fn

fn="../src/cvodes/README.md"
sedi "s/### Version.*/### Version ${cvs_ver} (${date})/" $fn
sedi "s/\"User Documentation for CVODES v.*/\"User Documentation for CVODES v${cvs_ver},\"/" $fn
sedi "s/UCRL-SM-208111,.*/UCRL-SM-208111, ${date}./" $fn
sedi "s/\"Example Programs for CVODES v.*/\"Example Programs for CVODES v${cvs_ver},\"/" $fn
sedi "s/UCRL-SM-208115,.*/UCRL-SM-208115, ${date}./" $fn

fn="../src/ida/README.md"
sedi "s/### Version.*/### Version ${ida_ver} (${date})/" $fn
sedi "s/\"User Documentation for IDA v.*/\"User Documentation for IDA v${ida_ver},\"/" $fn
sedi "s/UCRL-SM-208112,.*/UCRL-SM-208112, ${date}./" $fn
sedi "s/\"Example Programs for IDA v.*/\"Example Programs for IDA v${ida_ver},\"/" $fn
sedi "s/UCRL-SM-208113,.*/UCRL-SM-208113, ${date}./" $fn

fn="../src/idas/README.md"
sedi "s/### Version.*/### Version ${idas_ver} (${date})/" $fn
sedi "s/\"User Documentation for IDAS v.*/\"User Documentation for IDAS v${idas_ver},\"/" $fn
sedi "s/UCRL-SM-234051,.*/UCRL-SM-234051, ${date}./" $fn
sedi "s/\"Example Programs for IDAS v.*/\"Example Programs for IDAS v${idas_ver},\"/" $fn
sedi "s/LLNL-TR-437091,.*/LLNL-TR-437091, ${date}./" $fn

fn="../src/kinsol/README.md"
sedi "s/### Version.*/### Version ${kin_ver} (${date})/" $fn
sedi "s/\"User Documentation for KINSOL v.*/\"User Documentation for KINSOL v${kin_ver},\" LLNL technical report/" $fn
sedi "s/UCRL-SM-208116,.*/UCRL-SM-208116, ${date}./" $fn
sedi "s/\"Example Programs for KINSOL v.*/\"Example Programs for KINSOL v${kin_ver},\"/" $fn
sedi "s/UCRL-SM-208114,.*/UCRL-SM-208114, ${date}./" $fn

# ------------------------------------------------------------------------------
# Update tarscript
# ------------------------------------------------------------------------------

fn="tarscript.sh"
sedi "s/SUN_VER=.*/SUN_VER=\"${sun_ver}\"/"    $fn
sedi "s/CV_VER=.*/CV_VER=\"${cv_ver}\"/"       $fn
sedi "s/CVS_VER=.*/CVS_VER=\"${cvs_ver}\"/"    $fn
sedi "s/IDA_VER=.*/IDA_VER=\"${ida_ver}\"/"    $fn
sedi "s/IDAS_VER=.*/IDAS_VER=\"${idas_ver}\"/" $fn
sedi "s/KIN_VER=.*/KIN_VER=\"${kin_ver}\"/"    $fn
sedi "s/ARK_VER=.*/ARK_VER=\"${ark_ver}\"/"    $fn

# ------------------------------------------------------------------------------
# Update tex documentation
# ------------------------------------------------------------------------------

# update macros for SUNDIALS versions
fn="../doc/sundials/ug.tex"
sedi "s/sunrelease.*/sunrelease}{v${sun_ver}}/"    $fn
sedi "s/cvrelease.*/cvrelease}{v${cv_ver}}/"       $fn
sedi "s/cvsrelease.*/cvsrelease}{v${cvs_ver}}/"    $fn
sedi "s/idarelease.*/idarelease}{v${ida_ver}}/"    $fn
sedi "s/idasrelease.*/idasrelease}{v${idas_ver}}/" $fn
sedi "s/kinrelease.*/kinrelease}{v${kin_ver}}/"    $fn
sedi "s/arkrelease.*/arkrelease}{v${ark_ver}}/"    $fn

# update reference titles for user guides and example docs
for fn in "../doc/sundials/biblio.bib" "../doc/shared/sundials.bib";
do
    sedi "/User Documentation for ARKODE v/ s/v.*/v${ark_ver}}},/" $fn
    sedi "/Example Programs for ARKODE v/ s/v.*/v${ark_ver}}},/"   $fn
    sedi "/User Documentation for CVODE v/ s/v.*/v${cv_ver}}},/" $fn
    sedi "/Example Programs for CVODE v/ s/v.*/v${cv_ver}}},/"   $fn
    sedi "/User Documentation for CVODES v/ s/v.*/v${cvs_ver}}},/" $fn
    sedi "/Example Programs for CVODES v/ s/v.*/v${cvs_ver}}},/"   $fn
    sedi "/User Documentation for IDA v/ s/v.*/v${ida_ver}}},/" $fn
    sedi "/Example Programs for IDA v/ s/v.*/v${ida_ver}}},/"   $fn
    sedi "/User Documentation for IDAS v/ s/v.*/v${idas_ver}}},/" $fn
    sedi "/Example Programs for IDAS v/ s/v.*/v${idas_ver}}},/"   $fn
    sedi "/User Documentation for KINSOL v/ s/v.*/v${kin_ver}}},/" $fn
    sedi "/Example Programs for KINSOL v/ s/v.*/v${kin_ver}}},/"   $fn
    # update dates for user guides and example doc by checking lines between the
    # first and second latex comment patterns
    sedi "/% CURRENT.*/,/% ORIGINAL.*/ s/year.*/year        = ${year}/" $fn
done

# # insert new line in release table
# sedi '/%% Version Table/ a\
# '${month}' & '${year}' & '\
# ${sun_ver}' & '\
# ${ark_ver}' & '\
# ${cv_ver}' & '\
# ${cvs_ver}' & '\
# ${ida_ver}' & '\
# ${idas_ver}' & '\
# ${kin_ver}' \\\\'$'\n' ../doc/sundials/sundials_release_history.tex

# ------------------------------------------------------------------------------
# Update rst documentation
# ------------------------------------------------------------------------------

fn="../doc/shared/sundials_vars.py"
sedi "s/arkode_version =.*/arkode_version = \"v${ark_ver}\"/" $fn
sedi "s/cvode_version =.*/cvode_version = \"v${cv_ver}\"/" $fn
sedi "s/cvodes_version =.*/cvodes_version = \"v${cvs_ver}\"/" $fn
sedi "s/ida_version =.*/ida_version = \"v${ida_ver}\"/" $fn
sedi "s/idas_version =.*/idas_version = \"v${idas_ver}\"/" $fn
sedi "s/kinsol_version =.*/kinsol_version = \"v${kin_ver}\"/" $fn
sedi "s/sundials_version =.*/sundials_version = \"v${sun_ver}\"/" $fn
sedi "s/doc_version =.*/doc_version = \"v${sun_ver}\"/" $fn
sedi "s/year =.*/year = \"${year}\"/" $fn

# release history table
fn="../doc/shared/History.rst"
new_entry=$(printf "| %-3s %-4s | %-30s | %-17s | %-17s | %-17s | %-17s | %-17s | %-17s |" \
    ${month} ${year} ":ref:\`${sun_ver} <Changelog.${sun_ver}>\`" ${ark_ver} ${cv_ver} ${cvs_ver} ${ida_ver} \
    ${idas_ver} ${kin_ver})
divider="+----------+--------------------------------+-------------------+-------------------+-------------------+-------------------+-------------------+-------------------+"

# insert new release history row after line 26
sedi '26 a\
'"${divider}"''$'\n' $fn
sedi '26 a\
'"${new_entry}"''$'\n' $fn

# Update CITATIONS.md
fn="../CITATIONS.md"
sedi '71s/.*/\ \ year   = {'${year}'},/' $fn
sedi '72s/.*/\ \ note   = {v'${ark_ver}'}/' $fn
sedi '80s/.*/\ \ year   = {'${year}'},/' $fn
sedi '81s/.*/\ \ note   = {v'${cv_ver}'}/' $fn
sedi '89s/.*/\ \ year   = {'${year}'},/' $fn
sedi '90s/.*/\ \ note   = {v'${cvs_ver}'}/' $fn
sedi '98s/.*/\ \ year   = {'${year}'},/' $fn
sedi '99s/.*/\ \ note   = {v'${ida_ver}'}/' $fn
sedi '107s/.*/\ \ year   = {'${year}'},/' $fn
sedi '108s/.*/\ \ note   = {v'${idas_ver}'}/' $fn
sedi '116s/.*/\ \ year   = {'${year}'},/' $fn
sedi '117s/.*/\ \ note   = {v'${kin_ver}'}/' $fn

# Update all occurrences of x.y.z and X.Y.Z to the current version number
fn="../CHANGELOG.md"
sedi "s/x.y.z/${sun_ver}/gI" $fn

fn="../doc/shared/Changelog.rst"
sedi "s/x.y.z/${sun_ver}/gI" $fn

fn="../doc/shared/sundials/Install.rst"
sedi "s/x.y.z/${sun_ver}/gI" $fn

for fn in $(grep -Iirl "x.y.z" ../doc/shared/*)
do
    sedi "s/x.y.z/${sun_ver}/gI" $fn
done

for fn in $(grep -Iirl "x.y.z" ../doc/arkode/guide/source/*)
do
    sedi "s/x.y.z/${ark_ver}/gI" $fn
done

for fn in $(grep -Iirl "x.y.z" ../doc/cvode/guide/source/*)
do
    sedi "s/x.y.z/${cv_ver}/gI" $fn
done

for fn in $(grep -Iirl "x.y.z" ../doc/cvodes/guide/source/*)
do
    sedi "s/x.y.z/${cvs_ver}/gI" $fn
done

for fn in $(grep -Iirl "x.y.z" ../doc/ida/guide/source/*)
do
    sedi "s/x.y.z/${ida_ver}/gI" $fn
done

for fn in $(grep -Iirl "x.y.z" ../doc/idas/guide/source/*)
do
    sedi "s/x.y.z/${idas_ver}/gI" $fn
done

for fn in $(grep -Iirl "x.y.z" ../doc/kinsol/guide/source/*)
do
    sedi "s/x.y.z/${kin_ver}/gI" $fn
done
