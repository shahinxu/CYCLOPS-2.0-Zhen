#!/bin/bash
set -e

DATA_PATH="${1}"
OUTPUT_DIR="${2}"

if [ -z "$DATA_PATH" ]; then
    echo "Usage: $0 <data_path> [output_dir]"
    exit 1
fi

if [ ! -f "$DATA_PATH/expression.csv" ] || [ ! -f "$DATA_PATH/seed_genes.txt" ]; then
    echo "Error: Missing expression.csv or seed_genes.txt in $DATA_PATH"
    exit 1
fi

DATASET_NAME=$(basename "$DATA_PATH")

if [ -z "$OUTPUT_DIR" ]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    OUTPUT_DIR="./results/${DATASET_NAME}_${TIMESTAMP}"
fi

mkdir -p "$OUTPUT_DIR"

echo "Running CYCLOPS: $DATASET_NAME"
echo "Output: $OUTPUT_DIR"

sed "s|DATA_PATH_PLACEHOLDER|$DATA_PATH|g; s|OUTPUT_PATH_PLACEHOLDER|$OUTPUT_DIR|g" \
    run_cyclops.jl > temp_run.jl

julia temp_run.jl

rm -f temp_run.jl

echo "Completed: $OUTPUT_DIR"
