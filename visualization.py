"""
Figure 2: Benchmark concentration in protein structure prediction (2015-2025)
Requires: benchmark_concentration.csv
Output: figure2_final.pdf, figure2_final.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import kendalltau
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

df = pd.read_csv("benchmark_concentration.csv")
years = df["year"].values
hhi = df["hhi"].values
cr1 = df["cr1"].values
shannon = df["shannon"].values
total = df["total_papers"].values
casp_pct = df["count_CASP"].values / df["total_mentions"].values * 100

# Full-window Mann-Kendall
tau_h, p_h = kendalltau(years, hhi)
tau_c, p_c = kendalltau(years, cr1)
tau_s, p_s = kendalltau(years, shannon)

# Pre-2020 (for annotation only)
tau_pre, p_pre = kendalltau(years[:5], hhi[:5]) if len(set(hhi[:5])) > 1 else (0, 1)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 9,
    'axes.labelsize': 10, 'axes.titlesize': 11, 'figure.dpi': 300,
})
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
c1, c2 = '#2166AC', '#B2182B'

# ---- a: HHI ----
ax = axes[0, 0]
ax.plot(years, hhi, 'o-', color=c2, lw=2.5, ms=8, mfc='white', mew=2, mec=c2)
ax.fill_between(years, 0.30, hhi, alpha=0.10, color=c2)
ax.set_ylabel("Herfindahl-Hirschman Index", fontweight='bold')
ax.set_title("a. Benchmark concentration (HHI)", fontweight='bold', loc='left')
ax.set_ylim(0.30, 0.65); ax.set_xlim(2014.5, 2025.5)
ax.grid(True, alpha=0.25)
p_str = f"p = {p_h:.3f}" if p_h >= 0.001 else "p < 0.001"
ax.text(0.03, 0.95,
        f"Kendall's tau = {tau_h:.2f}, {p_str}",
        transform=ax.transAxes, fontsize=7.5, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.85, ec='#ccc'))
ax.annotate(f'{hhi[0]:.2f}', xy=(2015, hhi[0]), xytext=(0, 14),
            textcoords='offset points', fontsize=9, ha='center', color=c2, fontweight='bold')
ax.annotate(f'{hhi[-1]:.2f}', xy=(2025, hhi[-1]), xytext=(0, 14),
            textcoords='offset points', fontsize=9, ha='center', color=c2, fontweight='bold')

# ---- b: CASP share ----
ax = axes[0, 1]
ax.plot(years, casp_pct, 'o-', color=c1, lw=2.5, ms=8, mfc='white', mew=2, mec=c1)
ax.fill_between(years, 20, casp_pct, alpha=0.10, color=c1)
ax.set_ylabel("CASP share of mentions (%)", fontweight='bold')
ax.set_title("b. Dominant benchmark share (CASP)", fontweight='bold', loc='left')
ax.set_ylim(20, 85); ax.set_xlim(2014.5, 2025.5)
ax.grid(True, alpha=0.25)
p_c_str = f"p = {p_c:.3f}" if p_c >= 0.001 else "p < 0.001"
ax.text(0.03, 0.95, f"Kendall's tau = {tau_c:.2f}, {p_c_str}",
        transform=ax.transAxes, fontsize=8, va='top',
        bbox=dict(boxstyle='round', fc='white', alpha=0.85, ec='#ccc'))
ax.annotate(f'{casp_pct[0]:.0f}%', xy=(2015, casp_pct[0]), xytext=(0, 14),
            textcoords='offset points', fontsize=9, ha='center', color=c1, fontweight='bold')
ax.annotate(f'{casp_pct[-1]:.0f}%', xy=(2025, casp_pct[-1]), xytext=(0, 14),
            textcoords='offset points', fontsize=9, ha='center', color=c1, fontweight='bold')

# ---- c: Benchmark trajectories ----
ax = axes[1, 0]
bm_list = ["CASP", "CAMEO", "PDBbind", "CATH"]
colors = ['#2166AC', '#B2182B', '#4DAF4A', '#FF7F00']
markers = ['o', 's', 'D', '^']
for bm, c, m in zip(bm_list, colors, markers):
    ax.plot(years, df[f"count_{bm}"].values, marker=m, color=c, lw=1.8, ms=6,
            mfc='white', mew=1.5, label=bm)
ax.set_ylabel("Papers mentioning benchmark")
ax.set_title("c. Benchmark usage trajectories", fontweight='bold', loc='left')
ax.legend(fontsize=9, loc='upper left', framealpha=0.9)
ax.set_xlim(2014.5, 2025.5); ax.grid(True, alpha=0.25)

# ---- d: Field growth vs Shannon entropy ----
ax = axes[1, 1]
ax2 = ax.twinx()
ax.bar(years, total, color='#E8E8E8', edgecolor='#AAA', lw=0.5, zorder=1)
ax2.plot(years, shannon, 'o-', color=c2, lw=2.5, ms=8, mfc='white', mew=2, mec=c2, zorder=3)
ax2.fill_between(years, 0.7, shannon, alpha=0.10, color=c2)
ax.set_ylabel("Total papers in field", fontweight='bold')
ax2.set_ylabel("Shannon entropy (diversity)", fontweight='bold', color=c2)
ax.set_title("d. Field growth vs. diversity loss", fontweight='bold', loc='left')
ax.set_xlim(2014.5, 2025.5); ax2.set_ylim(0.70, 1.30)
ax.legend(handles=[
    Patch(facecolor='#E8E8E8', edgecolor='#AAA', label='Total papers'),
    Line2D([0], [0], color=c2, marker='o', ms=6, mfc='white', mew=2, lw=2.5, label='Shannon entropy'),
], fontsize=9, loc='upper left', framealpha=0.9)
p_s_str = f"p = {p_s:.3f}" if p_s >= 0.001 else "p < 0.001"
ax2.text(0.97, 0.95, f"Kendall's tau = {tau_s:.2f}, {p_s_str}",
         transform=ax.transAxes, fontsize=8, va='top', ha='right',
         bbox=dict(boxstyle='round', fc='white', alpha=0.85, ec='#ccc'))
ax2.annotate(f'{shannon[0]:.2f}', xy=(2015, shannon[0]), xytext=(0, 14),
             textcoords='offset points', fontsize=9, ha='center', color=c2, fontweight='bold')
ax2.annotate(f'{shannon[-1]:.2f}', xy=(2025, shannon[-1]), xytext=(0, -18),
             textcoords='offset points', fontsize=9, ha='center', color=c2, fontweight='bold')

fig.suptitle("Figure 2. Empirical evidence for benchmark-driven homogenization\n"
             "in protein structure prediction (2015-2025)",
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("figure2_final.pdf", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig("figure2_final.png", dpi=200, bbox_inches='tight', facecolor='white')
print("Figure 2 saved: figure2_final.pdf, figure2_final.png")
