#!/usr/bin/env python3
"""
Plot acrophase (peak phase) of clock genes for Zhang_all dataset, grouped by cell type.
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

def get_acrophase_data(predictions_file, expression_file, metadata_file):
    """
    Extract acrophase for each clock gene by cell type.
    """
    # Load data
    pred_df = pd.read_csv(predictions_file)
    expr_df = pd.read_csv(expression_file)
    meta_df = pd.read_csv(metadata_file)
    
    expr_df = expr_df.set_index('Gene_Symbol')
    
    # Sample mappings
    sample_to_phase = dict(zip(pred_df['Sample_ID'], pred_df['Predicted_Phase_Hours']))
    sample_to_celltype = dict(zip(meta_df['Sample'], meta_df['CellType_D']))
    
    # Get unique cell types
    cell_types = sorted(meta_df['CellType_D'].unique())
    
    results = {}
    
    for cell_type in cell_types:
        results[cell_type] = {}
        
        # Get samples for this cell type
        celltype_samples = meta_df[meta_df['CellType_D'] == cell_type]['Sample'].tolist()
        
        for gene in CLOCK_GENES:
            if gene not in expr_df.index:
                continue
            
            gene_expr = expr_df.loc[gene]
            times = []
            expressions = []
            
            for sample in celltype_samples:
                if sample in sample_to_phase and sample in gene_expr.index:
                    times.append(sample_to_phase[sample])
                    expressions.append(gene_expr[sample])
            
            if len(times) < 5:  # Need at least 5 points
                continue
            
            times = np.array(times)
            expressions = np.array(expressions)
            
            # Fit sinusoid
            params = fit_sinusoid(times, expressions)
            
            if params is not None:
                amplitude, phase, baseline = params
                
                # Calculate R²
                y_pred = sinusoidal_function(times, amplitude, phase, baseline)
                ss_res = np.sum((expressions - y_pred) ** 2)
                ss_tot = np.sum((expressions - np.mean(expressions)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                
                results[cell_type][gene] = {
                    'acrophase': phase,
                    'amplitude': amplitude,
                    'baseline': baseline,
                    'r_squared': r_squared,
                    'mean_expr': np.mean(expressions)
                }
    
    return results

def plot_acrophase_by_celltype(all_data, output_dir):
    """
    Plot acrophase for all genes in circular plots - one circle per cell type.
    """
    # Create figure with subplots for each cell type
    n_cell_types = len(all_data)
    n_cols = 3
    n_rows = (n_cell_types + n_cols - 1) // n_cols
    
    fig = plt.figure(figsize=(18, 6*n_rows))
    
    # Gene colors
    gene_colors = plt.cm.tab20(np.linspace(0, 1, len(CLOCK_GENES)))
    gene_to_color = dict(zip(CLOCK_GENES, gene_colors))
    
    for idx, cell_type in enumerate(sorted(all_data.keys())):
        ax = plt.subplot(n_rows, n_cols, idx+1, projection='polar')
        
        cell_data = all_data[cell_type]
        
        # Plot each gene
        for gene in CLOCK_GENES:
            if gene in cell_data:
                data = cell_data[gene]
                acrophase = data['acrophase']
                r_squared = data['r_squared']
                amplitude = data['amplitude']
                
                # Convert hours to radians (0 hour at top, clockwise)
                theta = np.pi/2 - (acrophase * 2 * np.pi / 24)
                
                # Plot point (radius = R²)
                color = gene_to_color[gene]
                
                # Size based on amplitude
                size = 150 + amplitude * 1500
                
                ax.scatter(theta, r_squared, s=size, c=[color], alpha=0.8,
                          edgecolors='black', linewidth=2, zorder=10)
                
                # Add gene label
                label_r = r_squared + 0.08
                if label_r > 0.95:
                    label_r = r_squared - 0.08
                
                ax.text(theta, label_r, gene, 
                       ha='center', va='center', fontsize=9, 
                       fontweight='bold', color='black',
                       bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor='white', alpha=0.8, edgecolor='none'))
        
        # Set up the circular plot
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        
        # Set hour labels
        hour_labels = ['0h\nMidnight', '3h', '6h', '9h', '12h\nNoon', '15h', '18h', '21h']
        ax.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
        ax.set_xticklabels(hour_labels, fontsize=10)
        
        # Set radial axis
        ax.set_ylim(0, 1)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8'], fontsize=9)
        ax.set_ylabel('R² (Rhythmicity)', fontsize=11, labelpad=30)
        
        # Title
        title_color = CELL_TYPE_COLORS.get(cell_type, 'black')
        ax.set_title(cell_type, fontsize=16, fontweight='bold', pad=20, color=title_color)
        
        # Grid
        ax.grid(True, alpha=0.3, linewidth=1)
        
        # Add radial lines at key times
        for hour in [0, 6, 12, 18]:
            theta_line = np.pi/2 - (hour * 2 * np.pi / 24)
            ax.plot([theta_line, theta_line], [0, 1], 'k--', alpha=0.2, linewidth=1)
    
    plt.suptitle('Clock Gene Acrophase by Cell Type - Zhang_all Dataset\n(Point size = Amplitude, Distance from center = R²)', 
                fontsize=18, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # Save
    output_file = os.path.join(output_dir, 'acrophase_by_celltype.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")
    plt.close()

def plot_acrophase_table(all_data, output_dir):
    """
    Create a table showing acrophase for each gene and cell type.
    """
    # Prepare data
    rows = []
    for gene in CLOCK_GENES:
        row = {'Gene': gene}
        for cell_type in sorted(CELL_TYPE_COLORS.keys()):
            if cell_type in all_data and gene in all_data[cell_type]:
                data = all_data[cell_type][gene]
                acrophase = data['acrophase']
                r_squared = data['r_squared']
                row[cell_type] = f"{acrophase:.1f}h (R²={r_squared:.2f})"
            else:
                row[cell_type] = 'N/A'
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Save as CSV
    output_file = os.path.join(output_dir, 'acrophase_table.csv')
    df.to_csv(output_file, index=False)
    print(f"Saved: {output_file}")
    
    # Print to console
    print("\n" + "="*100)
    print("ACROPHASE TABLE")
    print("="*100)
    print(df.to_string(index=False))
    print("="*100)

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
    
    if not os.path.exists(predictions_file):
        print(f"Error: predictions.csv not found")
        exit(1)
    
    if not os.path.exists(expression_file):
        print(f"Error: expression.csv not found")
        exit(1)
    
    if not os.path.exists(metadata_file):
        print(f"Error: metadata.csv not found")
        exit(1)
    
    print("\nExtracting acrophase data...")
    all_data = get_acrophase_data(predictions_file, expression_file, metadata_file)
    
    print(f"\nFound data for {len(all_data)} cell types")
    for cell_type, genes in all_data.items():
        print(f"  {cell_type}: {len(genes)} genes")
    
    print("\nCreating circular plots...")
    plot_acrophase_by_celltype(all_data, output_dir)
    
    print("\nCreating acrophase table...")
    plot_acrophase_table(all_data, output_dir)
    
    print("\n" + "="*100)
    print("All visualizations complete!")
    print("="*100)
