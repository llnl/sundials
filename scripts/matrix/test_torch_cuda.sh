#!/bin/bash

python -c 'import torch ; print(torch.rand(5, 3)) ; print("Torch Version", torch.__version__) ; print("Built for Cuda", torch.version.cuda) ; print("Cuda available:", torch.cuda.is_available())'

