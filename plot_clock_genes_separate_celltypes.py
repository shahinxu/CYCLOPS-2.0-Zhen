#!/usr/bin/env python3
"""
Plot expression of core clock genes with sinusoidal fit for Zhang_all dataset, one figure per cell type.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Core clock genes to plot
CLOCK_GENES = [
    'ARNTL',
    'CLOCK',
    'PER1',
    'PER2',
    'PER3',
    'CRY1',
    'CRY2',
    'NR1D1',
    'NR1D2',
    'DBP',
    'TEF',
    'HLF'
]

# Colors for different cell types
CELL_TYPE_COLORS = {
    'Bcell': '#1f77b4',
    'CD4Tcell': '#ff7f0e',
    'CD8Tcell': '#2ca02c',
    'Myeloid': '#d62728',
    'NKcell': '#9467bd'
}

def sinusoidal_function(x, amplitude, phase, baseline, period=24):
    """Sinusoidal function"""
    return baseline + amplitude * np.cos(2 * np.pi * (x - phase) / period)

def fit_sinusoid(times, expression):
    """Fit sinusoidal curve to expression data"""
    try:
        baseline_guess = np.mean(expression)
        amplitude_guess = (np.max(expression) - np.min(expression)) / 2
        phase_guess = times[np.argmax(expression)]
        
        params, _ = curve_fit(
            sinusoidal_function,
            times,
            expression,
            p0=[amplitude_guess, phase_guess, baseline_guess],
            bounds=(
                [0, 0, -np.inf],
                [np.inf, 24, np.inf]
            ),
            maxfev=10000
        )
        
        return params
    except:
        return None

def plot_clock_genes_for_celltype(predictions_file, expression_file, metadata_file, output_dir, cell_type):
    """
    Plot clock gene expression with sinusoidal fits for a specific cell type.
    """
    print(f"\n{'='*70}")
    print(f"Processing: {cell_type}")
    print(f"{'='*70}")
    
    # Load data
    pred_df = pd.read_csv(predictions_file)
    expr_df = pd.read_csv(expression_file)
    meta_df = pd.read_csv(metadata_file)
    
    # Set Gene_Symbol as index
    expr_df = expr_df.set_index('Gene_Symbol')
    
    # Sample to phase mapping
    sample_to_phase = dict(zip(pred_df['Sample_ID'], pred_df['Predicted_Phase_Hours']))
    
    # Get samples for this cell type
    celltype_samples = meta_df[meta_df['CellType_D'] == cell_type]['Sample'].tolist()
    print(f"Found {len(celltype_samples)} samples for {cell_type}")
    
    # Check which clock genes are present
    available_genes = [gene for gene in CLOCK_GENES if gene in expr_df.index]
    
    # Create figure with subplots for all genes
    n_genes = len(available_genes)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5*n_rows))
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    color = CELL_TYPE_COLORS.get(cell_type, 'gray')
    
    # Plot each gene in a subplot
    for idx, gene in enumerate(available_genes):
        ax = axes[idx]
        
        gene_expr = expr_df.loc[gene]
        
        times = []
        expressions = []
        
        for sample in celltype_samples:
            if sample in sample_to_phase and sample in gene_expr.index:
                times.append(sample_to_phase[sample])
                expressions.append(gene_expr[sample])
        
        if len(times) == 0:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12, color='gray')
            ax.set_title(gene, fontsize=12, fontweight='bold')
            continue
        
        times = np.array(times)
        expressions = np.array(expressions)
        
        # Sort by time
        sort_idx = np.argsort(times)
        times = times[sort_idx]
        expressions = expressions[sort_idx]
        
        # Plot scatter points
        ax.scatter(times, expressions, alpha=0.6, s=30, color=color,
                   edgecolors='black', linewidth=0.5, label='Data')
        
        # Fit sinusoid
        params = fit_sinusoid(times, expressions)
        
        if params is not None:
            amplitude, phase, baseline = params
            
            # Generate smooth curve
            time_smooth = np.linspace(0, 24, 500)
            expr_smooth = sinusoidal_function(time_smooth, amplitude, phase, baseline)
            
            ax.plot(time_smooth, expr_smooth, 'r-', linewidth=2, alpha=0.8,
                   label=f'Fit: A={amplitude:.2f}, φ={phase:.1f}h')
            
            # Calculate R²
            y_pred = sinusoidal_function(times, amplitude, phase, baseline)
            ss_res = np.sum((expressions - y_pred) ** 2)
            ss_tot = np.sum((expressions - np.mean(expressions)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            # Add R² text in corner
            ax.text(0.02, 0.98, f'R² = {r_squared:.3f}', 
                   transform=ax.transAxes, fontsize=9,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Formatting
        ax.set_xlabel('Predicted Phase (hours)', fontsize=10)
        ax.set_ylabel('Expression', fontsize=10)
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 24)
        
        # Add vertical lines
        for t in [0, 6, 12, 18, 24]:
            ax.axvline(t, color='gray', linestyle='--', alpha=0.2, linewidth=1)
    
    # Hide unused subplots
    for idx in range(len(available_genes), len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle(f'Clock Gene Expression - {cell_type}', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save combined figure
    output_file = os.path.join(output_dir, f'clock_genes_{cell_type}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    
    plt.close()

if __name__ == '__main__':
    # File paths
    predictions_file = r"D:\CriticalFile\projects\Circadian\TOAST\results\Zhang_CancerCell_2025_all\2025-12-03T17_57_49\predictions.csv"
    expression_file = r"D:\CriticalFile\projects\Circadian\data\Zhang_CancerCell_2025_all\expression.csv"
    metadata_file = r"D:\CriticalFile\projects\Circadian\data\Zhang_CancerCell_2025_all\metadata.csv"
    output_dir = r"D:\CriticalFile\projects\Circadian\TOAST\results\Zhang_CancerCell_2025_all\2025-12-03T17_57_49"
    
    print(f"Predictions: {predictions_file}")
    print(f"Expression: {expression_file}")
    print(f"Metadata: {metadata_file}")
    print(f"Output: {output_dir}")
    
    # Load metadata to get cell types
    meta_df = pd.read_csv(metadata_file)
    cell_types = sorted(meta_df['CellType_D'].unique())
    
    print(f"\nFound {len(cell_types)} cell types: {cell_types}")
    
    # Plot for each cell type
    for cell_type in cell_types:
        plot_clock_genes_for_celltype(predictions_file, expression_file, metadata_file, output_dir, cell_type)
    
    print(f"\n{'='*70}")
    print("All cell types processed!")
    print(f"{'='*70}")
