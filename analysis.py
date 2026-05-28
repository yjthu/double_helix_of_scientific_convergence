"""
Benchmark Concentration Analysis: Three Biological AI Subfields (2015–2025)
===========================================================================
Reproducible analysis for:
  "Computability shapes curiosity: the double helix of scientific convergence"
  Figure 2 & Supplementary Tables

Requirements: biopython, pandas, numpy, scipy, matplotlib
Usage: python analysis.py

Outputs: benchmark_concentration.csv (primary), benchmark_data_*.csv (multi-field)
"""

import sys, io, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

from Bio import Entrez
import pandas as pd
import numpy as np
from scipy.stats import kendalltau

Entrez.email = "your-email@example.com"  # CHANGE THIS before running
YEARS = list(range(2015, 2026))
RATE_LIMIT = 0.35

# ============================================================
# FIELD DEFINITIONS
# ============================================================
FIELDS = {
    "Protein structure": {
        "denom": (
            '("protein structure prediction"[Title/Abstract] '
            'OR "protein folding"[Title/Abstract] '
            'OR "protein design"[Title/Abstract] '
            'OR "protein engineering"[Title/Abstract]) '
            'AND ("deep learning"[Title/Abstract] '
            'OR "neural network"[Title/Abstract] '
            'OR "foundation model"[Title/Abstract] '
            'OR "transformer"[Title/Abstract] '
            'OR "AlphaFold"[Title/Abstract] '
            'OR "protein language model"[Title/Abstract])'
        ),
        "benchmarks": {
            "CASP": '"CASP" OR "CASP14" OR "CASP15" OR "CASP16"',
            "CAMEO": '"CAMEO"',
            "PDBbind": '"PDBbind" OR "PDB bind"',
            "ProteinGym": '"ProteinGym"',
            "CATH": '"CATH" AND ("protein" OR "domain")',
            "SCOP": '"SCOP" AND ("protein" OR "domain")',
        }
    },
    "Drug discovery": {
        "denom": (
            '("drug discovery"[Title/Abstract] '
            'OR "molecular docking"[Title/Abstract] '
            'OR "virtual screening"[Title/Abstract] '
            'OR "QSAR"[Title/Abstract]) '
            'AND ("deep learning"[Title/Abstract] '
            'OR "neural network"[Title/Abstract] '
            'OR "graph neural"[Title/Abstract] '
            'OR "transformer"[Title/Abstract] '
            'OR "foundation model"[Title/Abstract])'
        ),
        "benchmarks": {
            "DUD-E": '"DUD-E" OR "DUDE"',
            "ChEMBL": '"ChEMBL" AND ("benchmark" OR "dataset" OR "screen")',
            "MoleculeNet": '"MoleculeNet"',
            "Tox21": '"Tox21" AND ("benchmark" OR "dataset" OR "screen")',
        }
    },
    "Variant interpretation": {
        "denom": (
            '("variant interpretation"[Title/Abstract] '
            'OR "variant effect"[Title/Abstract] '
            'OR "pathogenicity prediction"[Title/Abstract] '
            'OR "missense variant"[Title/Abstract] '
            'OR "variant classification"[Title/Abstract]) '
            'AND ("deep learning"[Title/Abstract] '
            'OR "neural network"[Title/Abstract] '
            'OR "transformer"[Title/Abstract] '
            'OR "AlphaMissense"[Title/Abstract] '
            'OR "protein language model"[Title/Abstract])'
        ),
        "benchmarks": {
            "ClinVar": '"ClinVar" AND ("benchmark" OR "dataset" OR "predict")',
            "gnomAD": '"gnomAD" AND ("benchmark" OR "dataset" OR "predict")',
            "CAGI": '"CAGI" AND ("variant" OR "genom" OR "mutation")',
        }
    },
}


def count_papers(query, year=None):
    if year:
        full = f'({query}) AND ("{year}"[Date - Publication])'
    else:
        full = query
    try:
        handle = Entrez.esearch(db="pubmed", term=full, retmax=1)
        records = Entrez.read(handle)
        handle.close()
        return int(records["Count"])
    except:
        return 0


def analyze_field(field_name, field_def):
    print(f"\n{'='*60}")
    print(f"  {field_name} (2015-2025)")
    print(f"{'='*60}")

    denom = {}
    for y in YEARS:
        denom[y] = count_papers(field_def["denom"], y)
        time.sleep(RATE_LIMIT)

    bm_yearly = {}
    for bm_name, bm_query in field_def["benchmarks"].items():
        yearly = {}
        for y in YEARS:
            yearly[y] = count_papers(bm_query, y)
            time.sleep(RATE_LIMIT)
        bm_yearly[bm_name] = yearly
        total = sum(yearly.values())
        print(f"    {bm_name:15s}: total={total:5d} [{yearly[2015]:4d}->{yearly[2025]:4d}]")

    rows = []
    for y in YEARS:
        counts = {bm: bm_yearly[bm][y] for bm in field_def["benchmarks"]}
        total_m = sum(counts.values())
        if total_m > 0:
            shares = np.array([c / total_m for c in counts.values()])
            cr3 = np.sort(shares)[::-1][:3].sum()
            cr1 = shares[0]
            hhi = (shares ** 2).sum()
            shannon = -sum(s * np.log(s) for s in shares if s > 0)
            top_bm = max(counts, key=counts.get)
        else:
            cr3 = cr1 = hhi = shannon = 0
            top_bm = "N/A"
        rows.append({
            "year": y, "total": denom[y], "mentions": total_m,
            "cr3": cr3, "cr1": cr1, "hhi": hhi, "shannon": shannon,
            "top_bm": top_bm,
            **{f"count_{k}": v for k, v in counts.items()}
        })

    df = pd.DataFrame(rows)
    years_arr = np.array(YEARS)

    for metric, label in [("hhi", "HHI"), ("cr1", "CR1"), ("shannon", "Shannon")]:
        vals = df[metric].values
        tau, p = kendalltau(years_arr, vals)
        sig = "SIGNIFICANT" if p < 0.05 else "not significant"
        print(f"    {label}: 2015={vals[0]:.3f} -> 2025={vals[-1]:.3f}, "
              f"tau={tau:.3f}, p={p:.4f} ({sig})")

    hhi_vals = df["hhi"].values
    pre_v, post_v = hhi_vals[:5], hhi_vals[5:]
    tau_pre, p_pre = kendalltau(np.array(YEARS[:5]), pre_v) if len(set(pre_v)) > 1 else (0, 1)
    tau_post, p_post = kendalltau(np.array(YEARS[5:]), post_v)
    print(f"    HHI 2015-2019 (pre-AlphaFold2): {pre_v[0]:.3f}->{pre_v[-1]:.3f}, "
          f"tau={tau_pre:.2f}, p={p_pre:.3f}")
    print(f"    HHI 2020-2025 (post): {post_v[0]:.3f}->{post_v[-1]:.3f}, "
          f"tau={tau_post:.2f}, p={p_post:.4f}")

    fname = f"benchmark_data_{field_name.replace(' ', '_').lower()}.csv"
    df.to_csv(fname, index=False)
    if field_name == "Protein structure":
        df.to_csv("benchmark_concentration.csv", index=False)

    return df


if __name__ == "__main__":
    print("Benchmark Concentration Analysis: Three Biological AI Subfields")
    print(f"Time window: {YEARS[0]}-{YEARS[-1]}")
    print("=" * 60)

    all_results = {}
    for fname, fdef in FIELDS.items():
        all_results[fname] = analyze_field(fname, fdef)

    print(f"\n{'='*60}")
    print("SUMMARY: Multi-field benchmark concentration (2015-2025)")
    print(f"{'='*60}")
    print(f"{'Field':<25} {'HHI_2015':>8} {'HHI_2025':>8} {'tau':>6} {'p':>8} {'Sig':>5}")
    print(f"{'-'*60}")
    for fname, df in all_results.items():
        hhi0 = df["hhi"].iloc[0]
        hhi1 = df["hhi"].iloc[-1]
        tau, p = kendalltau(np.array(YEARS), df["hhi"].values)
        sig = "YES" if p < 0.05 else "no"
        print(f"{fname:<25} {hhi0:>8.3f} {hhi1:>8.3f} {tau:>6.3f} {p:>8.4f} {sig:>5}")

    print("\nDone. Run visualization.py to generate Figure 2.")
