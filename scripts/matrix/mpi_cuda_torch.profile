#!/bin/bash

module load gcc/13.3.1
module load mvapich2/2.3.7
module load python/3.12
module load cuda/12.6.0
module load cmake/3.30.5

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(CDPATH= cd -- "$script_dir/../.." && pwd)
VENV_DIR="$repo_dir/.venv-mpi-torch"

if [[ ! -d "$VENV_DIR" ]]; then
    python3 -m venv --system-site-packages "$VENV_DIR"
fi

# shellcheck disable=SC1090
source "$VENV_DIR/bin/activate"
