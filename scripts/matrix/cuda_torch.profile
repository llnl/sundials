#!/bin/bash

module load python/3.12

if [[ ! -d .venv-torch ]]; then
    python3 -m venv --system-site-packages .venv-torch
fi

source .venv-torch/bin/activate
