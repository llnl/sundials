#!/usr/bin/bash
# ------------------------------------------------------------------------------
# Promote containers to :latest
#
# Run this AFTER the new :spack-vX.Y.Z containers have been validated against
# the full SUNDIALS test suite. This retags them as :latest, which is what the
# CI workflows reference.
# ------------------------------------------------------------------------------

# Update to match version in build-all-containers.sh
SPACK_RELEASE="1.1.1"

images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')

set -x

for image in "${images[@]}"; do
  # Pull the validated spack-release image
  docker pull "ghcr.io/llnl/sundials-ci-${image}:spack-v${SPACK_RELEASE}"

  # Retag as :latest
  docker tag \
    "ghcr.io/llnl/sundials-ci-${image}:spack-v${SPACK_RELEASE}" \
    "ghcr.io/llnl/sundials-ci-${image}:latest"

  # Push the new :latest
  docker push "ghcr.io/llnl/sundials-ci-${image}:latest"
done

echo ""
echo "Done. The :latest tags now point to spack-v${SPACK_RELEASE} builds."
