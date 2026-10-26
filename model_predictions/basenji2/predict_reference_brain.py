#!/usr/bin/env python
"""
predict_reference_brain.py: Basenji2 on the unmodified reference genome.

This script writes the reference-genome arm (all 5,313 tracks) and also checks
that the model loaded correctly before the per-sample array is launched.

The check is needed because Keras load_weights onto an architecture built from
params does not raise an error if the layer order differs between versions.
The model still runs and returns softplus values in a plausible range, but they
are meaningless. Two models reading the same sequence should agree on which
genes are highly transcribed, so we compare with Enformer.

  PASS  median cross-gene Spearman vs Enformer > 0.6
  FAIL  near zero: the weights did not load; do not launch the array
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyfaidx
from scipy.stats import spearmanr

import basenji2_lib_brain as B


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tss_bed", default=B.TSS_BED)
    ap.add_argument("--out_dir", default=B.PRED_REFERENCE)
    ap.add_argument("--batch_size", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    outfile = f"{args.out_dir}/reference_genome.parquet"

    tss_df = B.load_tss_bed(args.tss_bed)
    print(f"Panel: {os.path.basename(args.tss_bed)}  ({len(tss_df)} genes)")

    tracks = B.load_tracks()
    names = B.track_names(tracks)
    print("Brain tracks:")
    for i, n in zip(tracks, names):
        print(f"  {i:5d}  {n}")

    if os.path.exists(outfile):
        print(f"\nReusing {outfile}")
        ref = pd.read_parquet(outfile)
    else:
        model = B.load_model()
        fa = pyfaidx.Fasta(B.REF)

        gene_ids, chunks = [], []
        rows = list(tss_df.itertuples(index=False))
        for start in range(0, len(rows), args.batch_size):
            batch = rows[start:start + args.batch_size]
            seqs, ok = [], []
            for row in batch:
                try:
                    seqs.append(B.get_seq(fa, row.chrom, row.tss, row.strand))
                    ok.append(row)
                except Exception as e:
                    print(f"  WARNING: {row.gene_name} failed: {e}")
            if not seqs:
                continue
            chunks.append(B.predict_batch(model, seqs))
            gene_ids.extend([r.gene_id for r in ok])
            if (start + args.batch_size) % 1000 < args.batch_size:
                print(f"  {start + args.batch_size}/{len(rows)}")

        fa.close()
        # keep all 5,313 tracks so other track choices can be tested later
        ref = pd.DataFrame(
            np.concatenate(chunks, axis=0),
            index=pd.Index(gene_ids, name="gene_id"),
            columns=[f"track_{i}" for i in range(B.N_TARGETS)],
        ).astype("float32")
        ref.to_parquet(outfile)
        print(f"\nSaved: {outfile}  {ref.shape}")

    # ── gate 1: dynamic range ────────────────────────────────────────────────
    brain = ref[[f"track_{i}" for i in tracks]]
    print("\nBrain-track summary:")
    print(brain.describe().T[["mean", "std", "min", "max"]].to_string())
    if float(brain.std().max()) < 1e-6:
        sys.exit("FAIL: predictions constant across genes; the model did not load correctly")

    # ── gate 2: agreement with Enformer ──────────────────────────────────────
    if not os.path.exists(B.ENFORMER_REF):
        print(f"\nSkipping Enformer cross-check ({B.ENFORMER_REF} not found)")
        return

    enf = pd.read_parquet(B.ENFORMER_REF)
    shared = brain.index.intersection(enf.index)
    print(f"\nCross-check vs Enformer on {len(shared)} shared genes:")

    rhos = []
    for i in range(min(len(brain.columns), len(enf.columns))):
        rho = spearmanr(brain.loc[shared].iloc[:, i],
                        enf.loc[shared].iloc[:, i],
                        nan_policy="omit").statistic
        rhos.append(rho)
        print(f"  {names[i]:34s} vs {str(enf.columns[i]):34s} rho = {rho:.3f}")

    med = float(np.nanmedian(rhos))
    print(f"\nMedian cross-gene Spearman vs Enformer: {med:.3f}")
    if med < 0.6:
        sys.exit("FAIL: agreement too low; check that model_human.h5 loaded "
                 "and that brain_tracks.txt indices match targets_human.txt")
    print("PASS: safe to launch the personalized runs (lsf/run_personal.sh)")


if __name__ == "__main__":
    main()