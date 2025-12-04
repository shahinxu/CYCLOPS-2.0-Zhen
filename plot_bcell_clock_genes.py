#!/usr/bin/env python3
"""
Plot expression of core clock genes with sinusoidal fit for Bcell dataset.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import sys

# Core clock genes to plot (human gene symbols)
CLOCK_GENES = [
    'ARNTL',
    'CLOCK',
    'NPAS2',
    'PER1',
    'PER2',
    'PER3',
    'CRY1',
    'CRY2',
    'NR1D1',
    'NR1D2',
    'DBP',
    'TEF',
    'HLF',
    'BHLHE41',
    'CIART',
    'RORC',
    'NFIL3'
]

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

def plot_clock_genes(fit_output_file, expression_file, output_dir):
    """
    Plot clock gene expression with sinusoidal fits for Bcell.
    """
    print(f"\n{'='*70}")
    print(f"Processing Bcell Clock Genes")
    print(f"{'='*70}")
    
    # Load data
    fit_df = pd.read_csv(fit_output_file)
    expr_df = pd.read_csv(expression_file)
    
    # Set Gene_Symbol as index
    expr_df = expr_df.set_index('Gene_Symbol')
    
    # Sample to phase mapping (convert radians to hours)
    sample_to_phase = {}
    for _, row in fit_df.iterrows():
        sample_id = row['ID']
        phase_radians = row['Phase']
        phase_hours = (phase_radians * 24) / (2 * np.pi)
        sample_to_phase[sample_id] = phase_hours
    
    print(f"Found {len(sample_to_phase)} samples with phase predictions")
    
    # Check which clock genes are present
    available_genes = [gene for gene in CLOCK_GENES if gene in expr_df.index]
    print(f"Found {len(available_genes)} clock genes in expression data: {available_genes}")
    
    if len(available_genes) == 0:
        print("ERROR: No clock genes found in expression data!")
        return
    
    # Create figure with subplots for all genes
    n_genes = len(available_genes)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5*n_rows))
    if n_genes == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    color = '#1f77b4'  # Blue for Bcell
    
    # Plot each gene in a subplot
    for idx, gene in enumerate(available_genes):
        ax = axes[idx]
        
        gene_expr = expr_df.loc[gene]
        
        times = []
        expressions = []
        
        for sample in gene_expr.index:
            if sample in sample_to_phase:
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
        ax.scatter(times, expressions, alpha=0.6, s=50, color=color,
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
                   transform=ax.transAxes, fontsize=10,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Formatting
        ax.set_xlabel('Predicted Phase (hours)', fontsize=10)
        ax.set_ylabel('Expression', fontsize=10)
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 24)
        ax.legend(fontsize=8, loc='upper right')
        
        # Add vertical lines at key times
        for t in [0, 6, 12, 18, 24]:
            ax.axvline(t, color='gray', linestyle='--', alpha=0.2, linewidth=1)
    
    # Hide unused subplots
    for idx in range(len(available_genes), len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Clock Gene Expression - Bcell', fontsize=18, fontweight='bold')
    plt.tight_layout()
    
    # Save combined figure
    output_file = os.path.join(output_dir, 'clock_genes_Bcell.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")
    
    plt.close()

def plot_acrophase_circular(fit_output_file, expression_file, output_dir):
    """
    Plot acrophase for clock genes in circular plot.
    """
    print(f"\nCreating circular acrophase plot...")
    
    # Load data
    fit_df = pd.read_csv(fit_output_file)
    expr_df = pd.read_csv(expression_file)
    expr_df = expr_df.set_index('Gene_Symbol')
    
    # Sample to phase mapping (convert radians to hours)
    sample_to_phase = {}
    for _, row in fit_df.iterrows():
        sample_id = row['ID']
        phase_radians = row['Phase']
        phase_hours = (phase_radians * 24) / (2 * np.pi)
        sample_to_phase[sample_id] = phase_hours
    
    # Check which clock genes are present
    available_genes = [gene for gene in CLOCK_GENES if gene in expr_df.index]
    
    if len(available_genes) == 0:
        print("ERROR: No clock genes found!")
        return
    
    # Calculate acrophase for each gene
    gene_data = {}
    
    for gene in available_genes:
        gene_expr = expr_df.loc[gene]
        
        times = []
        expressions = []
        
        for sample in gene_expr.index:
            if sample in sample_to_phase:
                times.append(sample_to_phase[sample])
                expressions.append(gene_expr[sample])
        
        if len(times) < 5:
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
            
            gene_data[gene] = {
                'acrophase': phase,
                'amplitude': amplitude,
                'baseline': baseline,
                'r_squared': r_squared,
                'mean_expr': np.mean(expressions)
            }
    
    if len(gene_data) == 0:
        print("ERROR: No genes could be fitted!")
        return
    
    # Create circular plot
    fig = plt.figure(figsize=(12, 12))
    ax = plt.subplot(111, projection='polar')
    
    # Gene colors
    gene_colors = plt.cm.tab20(np.linspace(0, 1, len(available_genes)))
    gene_to_color = dict(zip(available_genes, gene_colors))
    
    # Plot each gene
    for gene in gene_data:
        data = gene_data[gene]
        acrophase = data['acrophase']
        r_squared = data['r_squared']
        amplitude = data['amplitude']
        
        # Convert hours to radians (0 hour at top, clockwise)
        theta = np.pi/2 - (acrophase * 2 * np.pi / 24)
        
        # Plot point (radius = R²)
        color = gene_to_color[gene]
        
        # Size based on amplitude
        size = 200 + amplitude * 2000
        
        ax.scatter(theta, r_squared, s=size, c=[color], alpha=0.8,
                  edgecolors='black', linewidth=2, zorder=10)
        
        # Add gene label
        label_r = r_squared + 0.08
        if label_r > 0.95:
            label_r = r_squared - 0.08
        
        ax.text(theta, label_r, gene, 
               ha='center', va='center', fontsize=10, 
               fontweight='bold', color='black',
               bbox=dict(boxstyle='round,pad=0.3', 
                        facecolor='white', alpha=0.8, edgecolor='none'))
    
    # Set up the circular plot
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    
    # Set hour labels
    hour_labels = ['0h\nMidnight', '3h', '6h', '9h', '12h\nNoon', '15h', '18h', '21h']
    ax.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
    ax.set_xticklabels(hour_labels, fontsize=11)
    
    # Set radial axis
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8'], fontsize=10)
    ax.set_ylabel('R² (Rhythmicity)', fontsize=12, labelpad=30)
    
    # Title
    ax.set_title('Clock Gene Acrophase - Bcell\n(Point size = Amplitude, Distance from center = R²)', 
                fontsize=16, fontweight='bold', pad=20, color='#1f77b4')
    
    # Grid
    ax.grid(True, alpha=0.3, linewidth=1)
    
    # Add radial lines at key times
    for hour in [0, 6, 12, 18]:
        theta_line = np.pi/2 - (hour * 2 * np.pi / 24)
        ax.plot([theta_line, theta_line], [0, 1], 'k--', alpha=0.2, linewidth=1)
    
    plt.tight_layout()
    
    # Save
    output_file = os.path.join(output_dir, 'acrophase_circular_Bcell.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()
    
    # Save acrophase table
    rows = []
    for gene in sorted(gene_data.keys()):
        data = gene_data[gene]
        rows.append({
            'Gene': gene,
            'Acrophase_Hours': f"{data['acrophase']:.2f}",
            'Amplitude': f"{data['amplitude']:.3f}",
            'Baseline': f"{data['baseline']:.3f}",
            'R_squared': f"{data['r_squared']:.3f}",
            'Mean_Expression': f"{data['mean_expr']:.3f}"
        })
    
    df = pd.DataFrame(rows)
    output_csv = os.path.join(output_dir, 'acrophase_table_Bcell.csv')
    df.to_csv(output_csv, index=False)
    print(f"Saved: {output_csv}")
    print("\nAcrophase Summary:")
    print(df.to_string(index=False))

if __name__ == '__main__':
    # File paths for Bcell
    base_path = "/home/rzh/zhenx/circadian/CYCLOPS-2.0-Zhen"
    
    fit_output_file = f"{base_path}/output/Bcell/Fits/Fit_Output_2025-12-03T18_19_00.csv"
    expression_file = f"{base_path}/data/Zhang_CancerCell_2025_sub/Bcell/expression.csv"
    output_dir = f"{base_path}/output/Bcell/Plots"
    
    print(f"Fit Output: {fit_output_file}")
    print(f"Expression: {expression_file}")
    print(f"Output Dir: {output_dir}")
    
    if not os.path.exists(fit_output_file):
        print(f"ERROR: Fit output file not found: {fit_output_file}")
        sys.exit(1)
    
    if not os.path.exists(expression_file):
        print(f"ERROR: Expression file not found: {expression_file}")
        sys.exit(1)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    print("\n" + "="*70)
    print("Plotting clock gene expression patterns...")
    print("="*70)
    plot_clock_genes(fit_output_file, expression_file, output_dir)
    
    print("\n" + "="*70)
    print("Creating circular acrophase plot...")
    print("="*70)
    plot_acrophase_circular(fit_output_file, expression_file, output_dir)
    
    print("\n" + "="*70)
    print("All visualizations complete!")
    print("="*70)
