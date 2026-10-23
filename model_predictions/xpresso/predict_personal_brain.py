#!/usr/bin/env python
"""
predict_personal_brain.py: Xpresso version of enformer/predict_personal_brain.py.

Uses the same TSS bed, bcftools variant application and haplotype averaging (of
the predictions) as Enformer. Two differences from Enformer/Basenji2 come from
the model itself:
  - a 10,500 bp asymmetric window instead of a centered 131 kb/197 kb one
  - one output value per gene, so there are no bins or tracks to choose

Both arms are run in one pass. Building the haplotype FASTAs with bcftools takes
most of the run time (Xpresso itself has ~1M parameters), so the second model
call adds very little. Writes two parquets per sample:
    predictions/xpresso_seqonly_mane/<expr_id>.parquet
    predictions/xpresso_full_mane/<expr_id>.parquet
"""

import argparse
import os
import tempfile

import numpy as np
import pandas as pd
import pyfaidx

import xpresso_lib_brain as X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf_id", required=True, help="Sample ID in VCF, e.g. SAMPLE_001")
    ap.add_argument("--expr_id", required=True,
                    help="Sample ID in expression matrix, e.g. SAMPLE_001")
    ap.add_argument("--tss_bed", default=X.TSS_BED)
    ap.add_argument("--batch_size", type=int, default=64,
                    help="Genes per forward pass (each contributes 2 haplotypes)")
    ap.add_argument("--out_dir", default=None,
                    help="override the seqonly output dir; half-life arm goes to <out_dir>_full")
    args = ap.parse_args()
    # Optional override so the same script can target any output location.
    if args.out_dir:
        X.PRED_SEQONLY = args.out_dir
        X.PRED_FULL    = args.out_dir + "_full"


    os.makedirs(X.PRED_SEQONLY, exist_ok=True)
    os.makedirs(X.PRED_FULL, exist_ok=True)
    out_seq = f"{X.PRED_SEQONLY}/{args.expr_id}.parquet"
    out_full = f"{X.PRED_FULL}/{args.expr_id}.parquet"

    if os.path.exists(out_seq) and os.path.exists(out_full):
        print(f"Already exists, skipping: {args.expr_id}")
        return

    print(f"Sample: {args.vcf_id} / {args.expr_id}")

    tss_df = X.load_tss_bed(args.tss_bed)
    print(f"Genes to predict: {len(tss_df)}")

    print("Loading Xpresso...")
    model, n_hl = X.load_model()

    hl_map = {}
    hl_real, n_matched = X.load_halflife(tss_df["gene_id"].tolist(), n_hl)
    for gid, row in zip(tss_df["gene_id"], hl_real):
        hl_map[gid] = row
    have_full = n_matched > 0
    print(f"Half-life matched: {n_matched}/{len(tss_df)} genes "
          f"({'both arms' if have_full else 'seqonly only'})")

    res_seq, res_full = {}, {}
    genes_per_batch = max(1, args.batch_size // 2)   # 2 haplotypes per gene

    with tempfile.TemporaryDirectory() as tmpdir:
        for chrom in X.chrom_order(tss_df):
            genes = tss_df[tss_df["chrom"] == f"chr{chrom}"]
            if len(genes) == 0:
                continue

            print(f"  chr{chrom}: {len(genes)} genes, building personal sequences...")
            try:
                hap1_fa, hap2_fa = X.build_personal_fasta(args.vcf_id, chrom, tmpdir)
            except RuntimeError as e:
                print(f"  chr{chrom}: bcftools failed, skipping: {e}")
                for _, row in genes.iterrows():
                    res_seq[row.gene_id] = np.nan
                    res_full[row.gene_id] = np.nan
                continue

            fa1 = pyfaidx.Fasta(hap1_fa)
            fa2 = pyfaidx.Fasta(hap2_fa)
            rows = list(genes.itertuples(index=False))

            for start in range(0, len(rows), genes_per_batch):
                batch = rows[start:start + genes_per_batch]
                seqs, ok = [], []
                for row in batch:
                    try:
                        seqs.append(X.get_seq(fa1, f"chr{chrom}", row.tss, row.strand))
                        seqs.append(X.get_seq(fa2, f"chr{chrom}", row.tss, row.strand))
                        ok.append(row)
                    except Exception as e:
                        print(f"    WARNING: {row.gene_name} failed: {e}")
                        res_seq[row.gene_id] = np.nan
                        res_full[row.gene_id] = np.nan
                if not seqs:
                    continue

                # each gene appears twice (hap1, hap2) -> repeat its half-life row
                hl = np.repeat(
                    np.stack([hl_map[r.gene_id] for r in ok]), 2, axis=0
                ).astype("float32")

                p_seq = X.predict_batch(
                    model, seqs, np.zeros((len(seqs), n_hl), dtype="float32")
                )
                p_full = None
                if have_full:
                    p_full = X.predict_batch(model, seqs, np.nan_to_num(hl))
                    p_full[np.isnan(hl).any(axis=1)] = np.nan

                for i, row in enumerate(ok):
                    # haplotype averaging at the PREDICTION level, as in Enformer
                    res_seq[row.gene_id] = (p_seq[2 * i] + p_seq[2 * i + 1]) / 2
                    if have_full:
                        res_full[row.gene_id] = (p_full[2 * i] + p_full[2 * i + 1]) / 2

                done = start + len(batch)
                if done % 500 < genes_per_batch:
                    print(f"    {done}/{len(rows)} genes done")

            fa1.close()
            fa2.close()
            X.cleanup_fastas(hap1_fa, hap2_fa)
            print(f"  chr{chrom}: done ({len(genes)} genes)")

    for path, res, label in [(out_seq, res_seq, "seqonly"),
                             (out_full, res_full, "full")]:
        if not res:
            continue
        df = pd.DataFrame({X.COLUMN: pd.Series(res)}).astype("float32")
        df.index.name = "gene_id"
        df.to_parquet(path)
        print(f"Saved {label}: {path}  ({len(df)} genes, "
              f"{int(df[X.COLUMN].isna().sum())} NaN)")


if __name__ == "__main__":
    main()