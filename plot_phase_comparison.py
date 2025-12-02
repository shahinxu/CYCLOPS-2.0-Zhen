import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr
import os

def time_to_phase(time_hours, period_hours=24.0):
    """Convert time in hours to phase in radians [0, 2π]"""
    return (time_hours % period_hours) / period_hours * 2 * np.pi

def best_align_phase_for_comparison(phase_rad, metadata_rad, step=0.1):
    """Find best alignment shift to maximize correlation"""
    best_shift = 0
    best_corr = -np.inf
    
    shifts = np.arange(0, 2*np.pi, step)
    for shift in shifts:
        shifted = (phase_rad + shift) % (2 * np.pi)
        try:
            corr = pearsonr(shifted, metadata_rad)[0]
            if np.isfinite(corr) and corr > best_corr:
                best_corr = corr
                best_shift = shift
        except:
            continue
    
    aligned = (phase_rad + best_shift) % (2 * np.pi)
    return aligned, best_shift

fit_output = pd.read_csv('/home/azureuser/CYCLOPS-2.0/output/CYCLOPS_2025-12-01T22_27_00/Fits/Fit_Output_2025-12-01T22_27_00.csv')
metadata = pd.read_csv('/home/azureuser/CYCLOPS-2.0/GSE54652/white_adipose/metadata.csv')

merged = pd.merge(fit_output[['ID', 'Phases_MA']], metadata, left_on='ID', right_on='Sample')

metadata_rad = time_to_phase(merged['Time_Hours'].values, period_hours=24.0)
phase_rad = merged['Phases_MA'].values

aligned_rad, shift = best_align_phase_for_comparison(phase_rad, metadata_rad, step=0.1)

r = float(pearsonr(aligned_rad, metadata_rad)[0])
spearman_R = float(spearmanr(aligned_rad, metadata_rad)[0])
r2 = r * r if np.isfinite(r) else float('nan')

base_output_dir = '/home/azureuser/CYCLOPS-2.0/output/CYCLOPS_2025-12-01T22_27_00'
out_dir = os.path.join(base_output_dir, 'phase_vs_metadata')
os.makedirs(out_dir, exist_ok=True)

plt.figure(figsize=(8, 7))
plt.grid(True, linestyle='-')
plt.scatter(metadata_rad, aligned_rad, c='b', s=100)

two_pi = 2 * np.pi
plt.xlim(0, two_pi)
plt.ylim(0, two_pi)
plt.xlabel('Collection Phase', fontsize=24)
plt.ylabel('Predicted Phase', fontsize=24)

plt.tight_layout()

out_path = os.path.join(out_dir, 'comparison.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()

print(f"Pearson R={r:.2f}, Spearman ρ={spearman_R:.2f}, R²={r2:.2f}")
print(f"Best alignment shift: {shift:.2f} radians ({shift * 24 / (2*np.pi):.2f} hours)")
print(f"Plot saved in: {out_path}")

# 显示数据对比表
print("\n数据对比:")
comparison = pd.DataFrame({
    'Sample': merged['Sample'],
    'Time_Hours': merged['Time_Hours'],
    'Collection_Phase_rad': metadata_rad,
    'Original_Phase_rad': phase_rad,
    'Aligned_Phase_rad': aligned_rad
})
print(comparison.to_string(index=False))
