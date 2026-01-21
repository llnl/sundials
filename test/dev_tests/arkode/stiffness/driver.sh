#!/bin/bash

set -e

# Driver script for running stiffness_detection with various parameter combinations
# This script performs two parameter sweeps:
# 1. Vary beta with fixed gamma
# 2. Vary gamma with fixed beta

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
EXECUTABLE="./ark_shampine"
OUTPUT_DIR="results"

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo -e "${RED}Error: Executable '$EXECUTABLE' not found!${NC}"
    echo "Please compile the code first or update the EXECUTABLE path."
    exit 1
fi

echo -e "${BLUE}================================${NC}"
echo -e "${BLUE}Parameter Sweep Driver Script${NC}"
echo -e "${BLUE}================================${NC}"
echo ""

# ============================================================================
# Sweep 1: Vary beta with fixed gamma
# ============================================================================
echo -e "${GREEN}Sweep 1: Varying beta with fixed gamma${NC}"
echo "----------------------------------------"

FIXED_GAMMA=10
BETA_START=10
BETA_END=50
BETA_STEP=5

echo "Fixed gamma = $FIXED_GAMMA"
echo "Beta range: $BETA_START to $BETA_END (step $BETA_STEP)"
echo ""

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR/fixed_gamma"

for beta in $(seq $BETA_START $BETA_STEP $BETA_END); do
    test_case="beta_${beta}_gamma_${FIXED_GAMMA}"
    output_file="${OUTPUT_DIR}/fixed_gamma/${test_case}.dat"
    echo -e "Running: beta=${beta}, gamma=${FIXED_GAMMA} -> ${output_file}"
    export SUNLOGGER_INFO_FILENAME="${OUTPUT_DIR}/fixed_gamma/${test_case}.log"

    $EXECUTABLE -b $beta -g $FIXED_GAMMA -o "${output_file}"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Success${NC}"
        ./plot_shampine.py "${output_file}" --save "${OUTPUT_DIR}/fixed_gamma/${test_case}.pdf" --title "\$\beta\$=${beta}, \$\gamma\$=${FIXED_GAMMA}"
    else
        echo -e "${RED}Failed${NC}"
    fi
    echo ""
done

# ============================================================================
# Sweep 2: Vary gamma with fixed beta
# ============================================================================
echo ""
echo -e "${GREEN}Sweep 2: Varying gamma with fixed beta${NC}"
echo "----------------------------------------"

FIXED_BETA=100
GAMMA_START=60
GAMMA_END=130
GAMMA_STEP=10

echo "Fixed beta = $FIXED_BETA"
echo "Gamma range: $GAMMA_START to $GAMMA_END (step $GAMMA_STEP)"
echo ""

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR/fixed_beta"

for gamma in $(seq $GAMMA_START $GAMMA_STEP $GAMMA_END); do
    test_case="beta_${FIXED_BETA}_gamma_${gamma}"
    output_file="${OUTPUT_DIR}/fixed_beta/${test_case}.dat"
    echo -e "Running: beta=${FIXED_BETA}, gamma=${gamma} -> ${output_file}"
    export SUNLOGGER_INFO_FILENAME="${OUTPUT_DIR}/fixed_beta/${test_case}.log"

    $EXECUTABLE -b $FIXED_BETA -g $gamma -o "${output_file}"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Success${NC}"
        ./plot_shampine.py "${output_file}" --save "${OUTPUT_DIR}/fixed_beta/${test_case}.pdf" --title "\$\beta\$=${FIXED_BETA}, \$\gamma\$=${gamma}"
    else
        echo -e "${RED}Failed${NC}"
    fi
    echo ""
done

# ============================================================================
# Summary
# ============================================================================
echo ""
echo -e "${BLUE}================================${NC}"
echo -e "${BLUE}Sweep Complete!${NC}"
echo -e "${BLUE}================================${NC}"
echo ""
echo "Results saved in: $OUTPUT_DIR/"
echo ""
echo -e "${GREEN}All simulations completed successfully!${NC}"
