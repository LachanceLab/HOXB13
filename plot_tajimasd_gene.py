
import argparse
import os
import sys
import math
from typing import Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_tajima(path: str) -> pd.DataFrame:
    # Read TSV, coerce numeric columns, handle "NA" in tajima_d
    df = pd.read_csv(path, sep='\t', dtype=str)
    # Standardize expected columns
    expected_cols = {
        'pop': 'pop',
        'chromosome': 'chromosome',
        'chrom': 'chromosome',  # allow minor variation
        'window_pos_1': 'window_pos_1',
        'window_pos_2': 'window_pos_2',
        'tajima_d': 'tajima_d',
        'TajimaD': 'tajima_d',
    }
    # Rename if alternates exist
    rename = {}
    for c in df.columns:
        if c in expected_cols and c != expected_cols[c]:
            rename[c] = expected_cols[c]
    if rename:
        df = df.rename(columns=rename)

    # Required
    req = ['chromosome', 'window_pos_1', 'window_pos_2', 'tajima_d']
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in Tajima file: {missing}")

    # Coerce numeric
    df['window_pos_1'] = pd.to_numeric(df['window_pos_1'], errors='coerce').astype('Int64')
    df['window_pos_2'] = pd.to_numeric(df['window_pos_2'], errors='coerce').astype('Int64')
    df['tajima_d'] = pd.to_numeric(df['tajima_d'].replace({'NA': np.nan}), errors='coerce')

    # Drop rows missing coords
    df = df.dropna(subset=['chromosome', 'window_pos_1', 'window_pos_2'])
    return df

def read_exons(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', dtype=str)
    # Standardize expected columns
    expected = {
        'name2': 'gene',
        'gene': 'gene',
        'chrom': 'chromosome',
        'chromosome': 'chromosome',
        'txStart': 'txStart',
        'txEnd': 'txEnd',
        'transcript_length': 'transcript_length',
        'exonCount': 'exonCount',
        'total_exon_length': 'total_exon_length',
        'intron_count': 'intron_count',
        'total_intron_length': 'total_intron_length',
        'strand': 'strand'
    }
    rename = {}
    for c in df.columns:
        if c in expected and c != expected[c]:
            rename[c] = expected[c]
    if rename:
        df = df.rename(columns=rename)

    req = ['gene', 'chromosome', 'txStart', 'txEnd', 'total_exon_length', 'total_intron_length']
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in exon file: {missing}")

    # Coerce numeric
    for col in ['txStart', 'txEnd', 'total_exon_length', 'total_intron_length']:
        df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
    if 'transcript_length' in df.columns:
        df['transcript_length'] = pd.to_numeric(df['transcript_length'], errors='coerce').astype('Int64')
    else:
        df['transcript_length'] = (df['txEnd'] - df['txStart']).astype('Int64')

    # Derive intron:exon ratio and gene length
    df['intron_exon_ratio'] = (df['total_intron_length'].astype(float) /
                               df['total_exon_length'].replace({0: np.nan}).astype(float))
    df['exon_proportion'] = (df['total_exon_length'].astype(float) /
                               (df['total_intron_length'].replace({0: np.nan}).astype(float)+df['total_exon_length'].replace({0: np.nan}).astype(float)))
    df['gene_length'] = df['transcript_length'].astype('float')
    return df

def per_gene_mean_tajima(tajima: pd.DataFrame, exons: pd.DataFrame) -> pd.DataFrame:
    """
    Compute mean Tajima's D for each gene by overlapping windows with [txStart, txEnd].
    Efficient per-chromosome search using numpy.
    """
    results = []

    # Group windows by chromosome
    for chrom, tchr in tajima.groupby('chromosome'):
        if chrom not in set(exons['chromosome']):
            continue
        echr = exons[exons['chromosome'] == chrom].copy()

        # Drop NA tajima values to keep windows but we only average real numbers
        tchr = tchr.dropna(subset=['window_pos_1', 'window_pos_2']).copy()
        tchr = tchr.sort_values('window_pos_1')
        starts = tchr['window_pos_1'].to_numpy(dtype=np.int64)
        ends = tchr['window_pos_2'].to_numpy(dtype=np.int64)
        vals = tchr['tajima_d'].to_numpy(dtype=float)

        for idx, row in echr.iterrows():
            gstart = int(row['txStart'])
            gend = int(row['txEnd'])

            # Candidate slice: windows with start <= gend
            # Use searchsorted to get rightmost index where start <= gend
            right = np.searchsorted(starts, gend, side='right')
            if right == 0:
                mean_td = np.nan
            else:
                cand_ends = ends[:right]
                cand_vals = vals[:right]
                mask = cand_ends >= gstart  # overlap if window_end >= gene_start
                if np.any(mask):
                    td_vals = cand_vals[mask]
                    # Ignore NaNs in the mean
                    if td_vals.size > 0:
                        mean_td = float(np.nanmean(td_vals))
                    else:
                        mean_td = np.nan
                else:
                    mean_td = np.nan

            results.append({
                'gene': row['gene'],
                'chromosome': chrom,
                'txStart': gstart,
                'txEnd': gend,
                'mean_tajima_d': mean_td
            })

    out = pd.DataFrame(results)
    return out

def make_scatterplots(merged: pd.DataFrame, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    # 1) Tajima's D vs intron:exon ratio
    fig1 = plt.figure()
    x = merged['intron_exon_ratio'].to_numpy(dtype=float)
    y = merged['mean_tajima_d'].to_numpy(dtype=float)
    valid = ~np.isnan(x) & ~np.isnan(y) & np.isfinite(x) & np.isfinite(y)
    plt.scatter(x[valid], y[valid], s=10, alpha=0.7)
    plt.xlabel('Intron:Exon length ratio')
    plt.ylabel("Mean Tajima's D (per gene)")
    plt.title("Tajima's D vs Intron:Exon Ratio")
    f1 = os.path.join(outdir, 'tajimasD_vs_intronExonRatio.png')
    plt.tight_layout()
    plt.savefig(f1, dpi=200)
    plt.close(fig1)

    # 2) Tajima's D vs gene length (log10 scale on X)
    fig2 = plt.figure()
    gx = merged['gene_length'].to_numpy(dtype=float)
    gy = merged['mean_tajima_d'].to_numpy(dtype=float)
    valid2 = ~np.isnan(gx) & ~np.isnan(gy) & np.isfinite(gx) & np.isfinite(gy) & (gx > 0)
    plt.scatter(gx[valid2], gy[valid2], s=10, alpha=0.7)
    plt.xscale('log')
    plt.xlabel('Gene length (bp, log scale)')
    plt.ylabel("Mean Tajima's D (per gene)")
    plt.title("Tajima's D vs Gene Length")
    f2 = os.path.join(outdir, 'tajimasD_vs_geneLength.png')
    plt.tight_layout()
    plt.savefig(f2, dpi=200)
    plt.close(fig2)

    # 1) Tajima's D vs intron:exon ratio
    fig3 = plt.figure()
    x = merged['exon_proportion'].to_numpy(dtype=float)
    y = merged['mean_tajima_d'].to_numpy(dtype=float)
    valid = ~np.isnan(x) & ~np.isnan(y) & np.isfinite(x) & np.isfinite(y)
    plt.scatter(x[valid], y[valid], s=10, alpha=0.7)
    plt.xlabel('Exon proportion')
    plt.ylabel("Mean Tajima's D (per gene)")
    plt.title("Tajima's D vs Exon Proportion")
    f3 = os.path.join(outdir, 'tajimasD_vs_ExonProportion.png')
    plt.tight_layout()
    plt.savefig(f3, dpi=200)
    plt.close(fig3)

    return f1, f2, f3

def main():
    ap = argparse.ArgumentParser(description="Plot Tajima's D vs intron:exon ratio and gene length (per gene).")
    ap.add_argument('--tajima', required=True, help='Path to Tajima windows file (TSV).')
    ap.add_argument('--exons', required=True, help='Path to exon/gene summary file (TSV).')
    ap.add_argument('--outdir', default='plots', help='Output directory for plots and merged CSV.')
    args = ap.parse_args()

    tajima = read_tajima(args.tajima)
    exons = read_exons(args.exons)

    mean_td = per_gene_mean_tajima(tajima, exons)

    merged = exons.merge(mean_td[['gene', 'chromosome', 'mean_tajima_d']], on=['gene', 'chromosome'], how='left')
    out_csv = os.path.join(args.outdir, 'genes_with_tajimasD.csv')
    os.makedirs(args.outdir, exist_ok=True)
    merged.to_csv(out_csv, index=False)

    f1, f2, f3 = make_scatterplots(merged, args.outdir)

    print(f"Wrote merged table: {out_csv}")
    print(f"Saved plots: {f1}, {f2}, {f3}")

if __name__ == '__main__':
    main()
