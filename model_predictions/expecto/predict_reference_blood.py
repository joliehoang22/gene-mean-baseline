#!/usr/bin/env python
"""
predict_reference_blood.py: ExPecto on the unmodified reference genome.

This script writes the reference-genome arm and also checks the model before
the per-sample array is launched.

The check matters more for ExPecto than for the other models. ExPecto has three
stages (Beluga -> spatial transform -> XGBoost), and a mistake in the spatial
weights or the A,G,C,T one-hot order gives numbers that look plausible but are
meaningless. Agreement with the CAGE models on which genes are highly
transcribed is a quick check.

ExPecto's XGBoost stage was trained on GTEx expression for this tissue.
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import pyfaidx
from scipy.stats import spearmanr

import expecto_lib_blood as E


def _spearman(a, b, **kw):
    """Spearman rho, portable across scipy versions (1.7 has no `.statistic`)."""
    r = spearmanr(a, b, **kw)
    return r.statistic if hasattr(r, "statistic") else r.correlation


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tss_bed", default=E.TSS_BED)
    ap.add_argument("--out_dir", default=E.PRED_REFERENCE)
    ap.add_argument("--genes_per_batch", type=int, default=4,
                    help="Genes held in memory at once (each is 400 Beluga seqs)")
    ap.add_argument("--batch_size", type=int, default=512,
                    help="Beluga forward-pass batch size")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    outfile = f"{args.out_dir}/reference_genome.parquet"

    tss_df = E.load_tss_bed(args.tss_bed)
    print(f"Panel: {os.path.basename(args.tss_bed)}  ({len(tss_df)} genes)")
    print(f"Tissue: {E.TISSUE}")

    n_trunc = sum(E.tss_is_truncated(t) for t in tss_df["tss"])
    if n_trunc:
        print(f"NOTE: {n_trunc} gene(s) within {E.HALF}bp of chromosome start")

    if os.path.exists(outfile):
        print(f"\nReusing {outfile}")
        ref = pd.read_parquet(outfile)
    else:
        print("Loading Beluga + XGBoost...")
        beluga, bst, device = E.load_models()
        print(f"  device: {device}")
        if str(device) == "cpu":
            print("  WARNING: no GPU. 10,770 genes x 400 Beluga passes will take "
                  "many hours on CPU; use a GPU node.")

        fa = pyfaidx.Fasta(E.REF)
        rows = list(tss_df.itertuples(index=False))
        ids, preds = [], []
        t0 = time.time()

        for start in range(0, len(rows), args.genes_per_batch):
            batch = rows[start:start + args.genes_per_batch]
            seqs, strands, ok = [], [], []
            for row in batch:
                try:
                    seqs.append(E.get_seq(fa, row.chrom, row.tss, row.strand))
                    strands.append(row.strand)
                    ok.append(row)
                except Exception as e:
                    print(f"  WARNING: {row.gene_name} failed: {e}")
            if not seqs:
                continue

            preds.append(E.expecto_predict(beluga, bst, device, seqs, strands,
                                           args.batch_size))
            ids.extend(r.gene_id for r in ok)

            done = start + len(batch)
            if done % 200 < args.genes_per_batch:
                el = time.time() - t0
                rate = done / el
                print(f"  {done}/{len(rows)}  {rate:.1f} genes/s  "
                      f"eta {(len(rows) - done) / max(rate, 1e-9) / 60:.0f} min")

        fa.close()
        ref = pd.DataFrame({E.COLUMN: np.concatenate(preds)},
                           index=pd.Index(ids, name="gene_id")).astype("float32")
        ref.to_parquet(outfile)
        print(f"\nSaved: {outfile}  {ref.shape}")

    # ── gate 1: dynamic range ────────────────────────────────────────────────
    print("\nSummary:")
    print(ref.describe().T[["mean", "std", "min", "max"]].to_string())
    if float(ref[E.COLUMN].std()) < 1e-6:
        sys.exit("FAIL: predictions constant across genes; the pipeline is broken")

    # ── gate 2: agreement with the other models ──────────────────────────────
    checked = False
    for name, path in [("Enformer", E.ENFORMER_REF),
                       ("Basenji2", E.BASENJI2_REF),
                       ("Xpresso", E.XPRESSO_REF)]:
        if not os.path.exists(path):
            print(f"\nSkipping {name} cross-check ({path} not found)")
            continue
        other = pd.read_parquet(path)
        shared = ref.index.intersection(other.index)
        rhos = [_spearman(ref.loc[shared, E.COLUMN],
                          other.loc[shared].iloc[:, i], nan_policy="omit")
                for i in range(min(7, other.shape[1]))]
        med = float(np.nanmedian(rhos))
        print(f"\nvs {name} on {len(shared)} genes: median rho = {med:.3f}")
        if name == "Enformer":
            checked = True
            if med < 0.3:
                sys.exit("FAIL: agreement too low; check the Beluga weights, "
                         "the A,G,C,T one-hot ordering, and the spatial weights")

    if not checked:
        print("\nNo Enformer reference found; dynamic range looked sane.")
    print("\nPASS: safe to launch the personalized runs (lsf/run_personal.sh)")


if __name__ == "__main__":
    main()