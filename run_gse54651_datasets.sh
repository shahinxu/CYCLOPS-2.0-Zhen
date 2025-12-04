#!/bin/bash

# 设置基础路径
BASE_PATH="/home/rzh/zhenx/circadian/CYCLOPS-2.0-Zhen"
GSE54651_PATH="$BASE_PATH/data/GSE54651"
TEMPLATE_FILE="$BASE_PATH/CYCLOPS_2_0_Template.jl"
TEMP_SCRIPT="$BASE_PATH/temp_run_cyclops.jl"
JULIA_BIN="$BASE_PATH/julia-1.6.7/bin/julia"

# 获取所有数据集目录
all_datasets=($(ls -d $GSE54651_PATH/*/ | xargs -n 1 basename))

# 设置起始数据集
START_FROM="lung"

# 找到起始位置
start_index=0
for i in "${!all_datasets[@]}"; do
    if [ "${all_datasets[$i]}" == "$START_FROM" ]; then
        start_index=$i
        break
    fi
done

# 从起始位置开始的数据集
datasets=("${all_datasets[@]:$start_index}")

echo "Starting from dataset: $START_FROM"
echo "Found ${#datasets[@]} datasets to process (from lung onwards):"
for dataset in "${datasets[@]}"; do
    echo "  - $dataset"
done
echo ""

# 计数器
total=${#datasets[@]}
current=0
success=0
failed=0
skipped=0

# 循环处理每个数据集
for dataset in "${datasets[@]}"; do
    current=$((current + 1))
    echo "=========================================="
    echo "[$current/$total] Processing dataset: $dataset"
    echo "=========================================="
    
    data_path="$GSE54651_PATH/$dataset"
    
    # 检查必需文件是否存在
    if [ ! -f "$data_path/expression.csv" ]; then
        echo "[SKIP] $dataset: expression.csv not found"
        skipped=$((skipped + 1))
        continue
    fi
    
    if [ ! -f "$data_path/seed_genes.txt" ]; then
        echo "[SKIP] $dataset: seed_genes.txt not found"
        skipped=$((skipped + 1))
        continue
    fi
    
    if [ ! -f "$data_path/metadata.csv" ]; then
        echo "[SKIP] $dataset: metadata.csv not found"
        skipped=$((skipped + 1))
        continue
    fi
    
    # 创建临时脚本，替换 data_path
    sed "s|data_path = \".*\"|data_path = \"$data_path\"|g" "$TEMPLATE_FILE" > "$TEMP_SCRIPT"
    
    # 运行 Julia 脚本
    echo "Running CYCLOPS for $dataset..."
    start_time=$(date +%s)
    
    "$JULIA_BIN" "$TEMP_SCRIPT"
    exit_code=$?
    
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    
    if [ $exit_code -eq 0 ]; then
        success=$((success + 1))
        echo "[SUCCESS] $dataset completed in ${elapsed}s"
    else
        failed=$((failed + 1))
        echo "[ERROR] $dataset failed with exit code $exit_code"
    fi
    
    echo ""
done

# 清理临时文件
rm -f "$TEMP_SCRIPT"

echo "=========================================="
echo "All datasets processed!"
echo "Total: $total | Success: $success | Failed: $failed | Skipped: $skipped"
echo "=========================================="
