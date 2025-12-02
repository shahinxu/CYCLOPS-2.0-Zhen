#!/bin/bash

# 设置基础路径
BASE_PATH="/home/azureuser/CYCLOPS-2.0"
GSE54651_PATH="$BASE_PATH/GSE54651"
TEMPLATE_FILE="$BASE_PATH/CYCLOPS_2_0_Template.jl"
TEMP_SCRIPT="$BASE_PATH/temp_run_cyclops.jl"

datasets=($(ls -d $GSE54651_PATH/*/ | xargs -n 1 basename))

echo "Found ${#datasets[@]} datasets to process:"
for dataset in "${datasets[@]}"; do
    echo "  - $dataset"
done
echo ""

# 循环处理每个数据集
for dataset in "${datasets[@]}"; do
    echo "=========================================="
    echo "Processing dataset: $dataset"
    echo "=========================================="
    
    data_path="$GSE54651_PATH/$dataset"
    
    # 检查必需文件是否存在
    if [ ! -f "$data_path/expression.csv" ]; then
        echo "[SKIP] $dataset: expression.csv not found"
        continue
    fi
    
    if [ ! -f "$data_path/seed_genes.txt" ]; then
        echo "[SKIP] $dataset: seed_genes.txt not found"
        continue
    fi
    
    # 创建临时脚本，替换 data_path
    sed "s|data_path = \".*\"|data_path = \"$data_path\"|g" "$TEMPLATE_FILE" > "$TEMP_SCRIPT"
    
    # 运行 Julia 脚本
    echo "Running CYCLOPS for $dataset..."
    start_time=$(date +%s)
    
    julia "$TEMP_SCRIPT"
    exit_code=$?
    
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    
    if [ $exit_code -eq 0 ]; then
        echo "[SUCCESS] $dataset completed in ${elapsed}s"
    else
        echo "[ERROR] $dataset failed with exit code $exit_code"
    fi
    
    echo ""
done

# 清理临时文件
rm -f "$TEMP_SCRIPT"

echo "=========================================="
echo "All datasets processed!"
echo "=========================================="
