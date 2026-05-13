#!/usr/bin/bash
# ------------------------------------------------------------------------------
# Pull all SUNDIALS CI containers from GHCR.
# ------------------------------------------------------------------------------

# Update to match version in build-all-containers.sh
SPACK_RELEASE="v1.1.1"

images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')

set -x

for image in "${images[@]}"; do
  docker pull "ghcr.io/llnl/sundials-ci-${image}:spack-${SPACK_RELEASE}"
  docker pull "ghcr.io/llnl/sundials-ci-${image}:spack-develop"
done
