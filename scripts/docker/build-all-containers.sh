#!/usr/bin/bash
# ------------------------------------------------------------------------------
# Build all SUNDIALS CI containers.
#
# Two tags:
#   :spack-vX.Y.Z  - built against a pinned Spack release
#   :spack-develop - built against Spack develop
#
# These tags intentionally do NOT overwrite the existing :latest tag. The
#   existing :latest containers remain available while the new containers are
#   validated. Once validated, re-tag the new containers as :latest.
#
# Update SPACK_RELEASE when a new Spack version is tagged. Also update the
# spack-public mirror URL in each spack.yaml to match.
# ------------------------------------------------------------------------------

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain

SPACK_RELEASE="v1.1.1"
SPACK_RELEASE_IMAGE="spack/ubuntu-noble:${SPACK_RELEASE}"
SPACK_DEVELOP_IMAGE="spack/ubuntu-noble:develop"

images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')

set -x

# Build containers pinned to a Spack release (:spack-vX.Y.Z tag)
for image in "${images[@]}"; do
  docker build \
    -t "ghcr.io/llnl/sundials-ci-${image}:spack-${SPACK_RELEASE}" \
    --build-arg "SPACK_BASE_IMAGE=${SPACK_RELEASE_IMAGE}" \
    --build-arg "spack_yaml=${image}/spack.yaml" \
    .
done

# Build containers against Spack develop (:spack-develop tag)
for image in "${images[@]}"; do
  docker build \
    -t "ghcr.io/llnl/sundials-ci-${image}:spack-develop" \
    --build-arg "SPACK_BASE_IMAGE=${SPACK_DEVELOP_IMAGE}" \
    --build-arg "spack_yaml=${image}/spack.yaml" \
    .
done
