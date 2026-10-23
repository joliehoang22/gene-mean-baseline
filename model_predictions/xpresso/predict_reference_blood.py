#!/usr/bin/env python
"""
predict_reference_blood.py: Xpresso on the unmodified reference genome, both arms.

This script writes the reference-genome arm for `seqonly` and `full` and also
checks the model before the per-sample array is launched.

Both arms are computed in one pass: the promoter one-hot is built once and the
model is called twice, with different half-life inputs. Building the sequence
takes most of the run time, so the second arm adds very little.

Check: Xpresso predicts RNA-seq abundance from a 10.5 kb promoter window, while
Enformer/Basenji2 predict CAGE from 131-197 kb. The correlation with them should
be positive but well below the Basenji2-vs-Enformer value, since the assay and
the input length differ. A correlation near zero means the weights did not load
correctly or the window is wrong.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyfaidx
from scipy.stats import spearmanr

import xpresso_lib_blood as X


def rho(a, b, **kw):
    """
    Spearman rho, portable across scipy versions.

    `.statistic` only exists in scipy >= 1.9, which needs Python 3.8+. This env
    is Python 3.7 (matching Xpresso's TF 1.15 requirement), so scipy is 1.7.x
    where the result is a namedtuple exposing `.correlation`.
    """
    r = spearmanr(a, b, **kw)
    return getattr(r, "statistic", None) if hasattr(r, "statistic") else r.correlation


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tss_bed", default=X.TSS_BED)
    ap.add_argument("--out_dir", default=os.path.dirname(X.REF_SEQONLY))
    ap.add_argument("--batch_size", type=int, default=64)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    out_seq = f"{args.out_dir}/seqonly.parquet"
    out_full = f"{args.out_dir}/full.parquet"

    tss_df = X.load_tss_bed(args.tss_bed)
    print(f"Panel: {os.path.basename(args.tss_bed)}  ({len(tss_df)} genes)")

    n_trunc = sum(X.tss_is_truncated(r.tss, r.strand)
                  for r in tss_df.itertuples(index=False))
    if n_trunc:
        print(f"NOTE: {n_trunc} gene(s) too close to chromosome start; "
              f"window not fully upstream for those")

    print("Loading Xpresso...")
    model, n_hl = X.load_model()
    print(f"  window {X.WINDOW}bp, half-life features: {n_hl}")

    gene_ids = tss_df["gene_id"].tolist()
    hl_real, n_matched = X.load_halflife(gene_ids, n_hl)
    print(f"  half-life matched: {n_matched}/{len(gene_ids)} genes")
    have_full = n_matched > 0

    fa = pyfaidx.Fasta(X.REF)
    ids, p_seq, p_full = [], [], []
    rows = list(tss_df.itertuples(index=False))

    for start in range(0, len(rows), args.batch_size):
        batch = rows[start:start + args.batch_size]
        seqs, idx = [], []
        for i, row in enumerate(batch):
            try:
                seqs.append(X.get_seq(fa, row.chrom, row.tss, row.strand))
                idx.append(start + i)
            except Exception as e:
                print(f"  WARNING: {row.gene_name} failed: {e}")
        if not seqs:
            continue

        zeros = np.zeros((len(seqs), n_hl), dtype="float32")
        p_seq.append(X.predict_batch(model, seqs, zeros))

        if have_full:
            hl = hl_real[idx]
            # NaN rows = gene absent from Xpresso's table; predict then mask
            p_full.append(X.predict_batch(model, seqs, np.nan_to_num(hl)))
            p_full[-1][np.isnan(hl).any(axis=1)] = np.nan

        ids.extend(rows[i].gene_id for i in idx)
        if (start + args.batch_size) % 2000 < args.batch_size:
            print(f"  {start + args.batch_size}/{len(rows)}")

    fa.close()

    seq_df = pd.DataFrame({X.COLUMN: np.concatenate(p_seq)},
                          index=pd.Index(ids, name="gene_id")).astype("float32")
    seq_df.to_parquet(out_seq)
    print(f"\nSaved: {out_seq}  {seq_df.shape}")

    full_df = None
    if have_full:
        full_df = pd.DataFrame({X.COLUMN: np.concatenate(p_full)},
                               index=pd.Index(ids, name="gene_id")).astype("float32")
        full_df.to_parquet(out_full)
        print(f"Saved: {out_full}  {full_df.shape} "
              f"({int(full_df[X.COLUMN].isna().sum())} NaN)")

    # ── gate 1: dynamic range ────────────────────────────────────────────────
    print("\nSeq-only summary:")
    print(seq_df.describe().T[["mean", "std", "min", "max"]].to_string())
    if float(seq_df[X.COLUMN].std()) < 1e-6:
        sys.exit("FAIL: predictions are constant across genes; the model did not load correctly")

    # ── seqonly vs full ──────────────────────────────────────────────────────
    if full_df is not None:
        both = seq_df.join(full_df, rsuffix="_full").dropna()
        r = rho(both.iloc[:, 0], both.iloc[:, 1])
        print(f"\nseqonly vs full: rho = {r:.3f} on {len(both)} genes")
        print("  (the gap between these arms shows how much of Xpresso's cross-gene")
        print("   signal comes from gene-level annotation)")

    # ── gate 2: agreement with the CAGE models ───────────────────────────────
    for name, path in [("Enformer", X.ENFORMER_REF), ("Basenji2", X.BASENJI2_REF)]:
        if not os.path.exists(path):
            print(f"\nSkipping {name} cross-check ({path} not found)")
            continue
        other = pd.read_parquet(path)
        shared = seq_df.index.intersection(other.index)
        rhos = [rho(seq_df.loc[shared, X.COLUMN],
                    other.loc[shared].iloc[:, i],
                    nan_policy="omit")
                for i in range(min(7, other.shape[1]))]
        med = float(np.nanmedian(rhos))
        print(f"\nvs {name} on {len(shared)} genes: median rho = {med:.3f}")
        if name == "Enformer" and med < 0.3:
            sys.exit("FAIL: agreement too low; check the model loaded and that "
                     "the 10,500bp window geometry is right")

    print("\nPASS: safe to launch the personalized runs (lsf/run_personal.sh)")


if __name__ == "__main__":
    main()