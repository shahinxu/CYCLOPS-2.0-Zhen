#!/bin/bash
set -e

DATA_PATH="${1}"
OUTPUT_DIR="${2}"

if [ -z "$DATA_PATH" ]; then
    echo "Usage: $0 <data_path> [output_dir]"
    exit 1
fi

# Check if this is a direct dataset or has subdatasets
if [ -f "$DATA_PATH/expression.csv" ] && [ -f "$DATA_PATH/seed_genes.txt" ]; then
    # Direct dataset - run once
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
else
    # Has subdatasets - run each subdirectory
    echo "No expression.csv found. Processing subdatasets in $DATA_PATH"
    
    for SUBDIR in "$DATA_PATH"/*/; do
        if [ -f "$SUBDIR/expression.csv" ] && [ -f "$SUBDIR/seed_genes.txt" ]; then
            SUBNAME=$(basename "$SUBDIR")
            PARENT_NAME=$(basename "$DATA_PATH")
            
            if [ -z "$OUTPUT_DIR" ]; then
                TIMESTAMP=$(date +%Y%m%d_%H%M%S)
                SUB_OUTPUT="./results/${PARENT_NAME}_${SUBNAME}_${TIMESTAMP}"
            else
                SUB_OUTPUT="${OUTPUT_DIR}/${SUBNAME}"
            fi
            
            mkdir -p "$SUB_OUTPUT"
            
            echo "Running CYCLOPS: $PARENT_NAME/$SUBNAME"
            echo "Output: $SUB_OUTPUT"
            
            sed "s|DATA_PATH_PLACEHOLDER|$SUBDIR|g; s|OUTPUT_PATH_PLACEHOLDER|$SUB_OUTPUT|g" \
                run_cyclops.jl > temp_run.jl
            
            julia temp_run.jl
            
            rm -f temp_run.jl
        fi
    done
fi

echo "Completed: $OUTPUT_DIR"
