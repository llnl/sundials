#!/bin/bash
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
# Publish sundials4py wheel artifacts to PyPI or TestPyPI.
# ------------------------------------------------------------------------------

set -euo pipefail

print_usage() {
    cat <<EOF

Usage: publish_sundials4py_wheels.sh [options]

Options:
  -h, --help             Show this message.
  -d, --dist-dir DIR     Directory containing built artifacts.
                         Default: <repo>/dist
  -r, --repository NAME  Repository to publish to: pypi or testpypi.
                         Default: pypi
      --repository-url   Override the repository upload URL.
      --include-sdist    Upload the sundials4py source distribution too.
      --skip-existing    Pass --skip-existing to twine upload.
      --dry-run          Validate artifacts but do not upload them.

Environment:
  PYPI_TOKEN or TWINE_PASSWORD
      API token used for upload. Required unless --dry-run is used.
  TWINE_USERNAME
      Twine username. Defaults to __token__.

Examples:
  ./scripts/publish_sundials4py_wheels.sh
  ./scripts/publish_sundials4py_wheels.sh --repository testpypi --dry-run
  PYPI_TOKEN=... ./scripts/publish_sundials4py_wheels.sh --include-sdist

EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/.." && pwd)

dist_dir="${repo_root}/dist"
repository="pypi"
repository_url=""
include_sdist=F
skip_existing=F
dry_run=F

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            print_usage
            exit 0
            ;;
        -d|--dist-dir)
            [[ $# -ge 2 ]] || die "Missing value for $1"
            dist_dir="$2"
            shift 2
            ;;
        -r|--repository)
            [[ $# -ge 2 ]] || die "Missing value for $1"
            repository="$2"
            shift 2
            ;;
        --repository-url)
            [[ $# -ge 2 ]] || die "Missing value for $1"
            repository_url="$2"
            shift 2
            ;;
        --include-sdist)
            include_sdist=T
            shift
            ;;
        --skip-existing)
            skip_existing=T
            shift
            ;;
        --dry-run)
            dry_run=T
            shift
            ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done

case "$repository" in
    pypi)
        default_repository_url="https://upload.pypi.org/legacy/"
        ;;
    testpypi)
        default_repository_url="https://test.pypi.org/legacy/"
        ;;
    *)
        die "Unsupported repository '${repository}'. Use 'pypi' or 'testpypi'."
        ;;
esac

if [[ -z "${repository_url}" ]]; then
    repository_url="${default_repository_url}"
fi

[[ -d "${dist_dir}" ]] || die "Artifact directory does not exist: ${dist_dir}"
dist_dir=$(cd -- "${dist_dir}" && pwd)

python_cmd=""
if command -v python3 >/dev/null 2>&1; then
    python_cmd="python3"
elif command -v python >/dev/null 2>&1; then
    python_cmd="python"
else
    die "Python is required to run twine."
fi

"${python_cmd}" -m twine --version >/dev/null 2>&1 \
    || die "twine is not installed. Install it with '${python_cmd} -m pip install twine'."

shopt -s nullglob
wheel_files=("${dist_dir}"/sundials4py-*.whl)
sdist_files=("${dist_dir}"/sundials4py-*.tar.gz)
shopt -u nullglob

(( ${#wheel_files[@]} > 0 )) || die "No sundials4py wheels found in ${dist_dir}"

artifacts=("${wheel_files[@]}")
if [[ "${include_sdist}" == "T" ]]; then
    (( ${#sdist_files[@]} > 0 )) || die "No sundials4py source distribution found in ${dist_dir}"
    artifacts+=("${sdist_files[@]}")
fi

echo "Publishing sundials4py artifacts"
echo "  repository: ${repository}"
echo "  url: ${repository_url}"
echo "  dist dir: ${dist_dir}"
echo "  upload: ${#artifacts[@]} file(s)"
printf '  artifact: %s\n' "${artifacts[@]##*/}"

"${python_cmd}" -m twine check "${artifacts[@]}"

if [[ "${dry_run}" == "T" ]]; then
    echo "Dry run requested; skipping upload."
    exit 0
fi

twine_password="${TWINE_PASSWORD:-${PYPI_TOKEN:-}}"
[[ -n "${twine_password}" ]] || die "Set PYPI_TOKEN or TWINE_PASSWORD before uploading."

export TWINE_USERNAME="${TWINE_USERNAME:-__token__}"
export TWINE_PASSWORD="${twine_password}"

upload_args=(
    --non-interactive
    --repository-url "${repository_url}"
)

if [[ "${skip_existing}" == "T" ]]; then
    upload_args+=(--skip-existing)
fi

"${python_cmd}" -m twine upload "${upload_args[@]}" "${artifacts[@]}"
