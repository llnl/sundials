#!/usr/bin/bash

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain
images=('int32-double' 'int32-extended' 'int32-single' 'int64-double' 'int64-extended' 'int64-single')
set -x
docker build -t ghcr.io/llnl/sundials-ci-e4s-base:latest e4s-base/
docker build -t ghcr.io/llnl/sundials-ci-e4s-base:e4s-22.05 e4s-base/
for image in ${images[@]}; do
  docker build -t ghcr.io/llnl/sundials-ci-${image}:latest --build-arg spack_yaml="${image}/spack.yaml" e4s-quarterly/
done

docker build -t ghcr.io/llnl/sundials-ci-${image}:spack-develop --build-arg spack_yaml="${image}/spack.yaml" spack-nightly/
docker build -t ghcr.io/llnl/sundials-ci-${image}:spack-develop --build-arg spack_yaml="${image}/spack.yaml" spack-nightly/
