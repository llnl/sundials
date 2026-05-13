#!/usr/bin/bash
# ------------------------------------------------------------------------------
# Push all SUNDIALS CI containers to GHCR
# ------------------------------------------------------------------------------

# Update to match version in build-all-containers.sh
SPACK_RELEASE="1.1.1"

images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')

set -x

for image in "${images[@]}"; do
  docker push "ghcr.io/llnl/sundials-ci-${image}:spack-v${SPACK_RELEASE}"
  docker push "ghcr.io/llnl/sundials-ci-${image}:spack-develop"
done
