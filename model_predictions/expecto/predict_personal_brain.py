#!/usr/bin/env python
"""
predict_personal_brain.py: ExPecto version of enformer/predict_personal_brain.py.

Uses the same TSS bed, bcftools variant application and haplotype averaging (of
the predictions) as Enformer. Where ExPecto's own conventions differ, they are
kept: genomic-orientation window, strand-flipped shifts, A,G,C,T one-hot, and
forward/reverse-complement averaging inside Beluga.

Cost: 400 Beluga forward passes per gene per haplotype, about 10.7M per sample
for the 13,325-gene panel, so a GPU is needed. The ni-lab run saved the full
Beluga output (8.1 TB). Here the 10 spatial basis functions are applied right
away and one value is kept per gene, so the output is one parquet per sample,
as for the other models.
"""

import argparse
import os
import tempfile
import time

import numpy as np
import pandas as pd
import pyfaidx

import expecto_lib_brain as E


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf_id", required=True, help="Sample ID in VCF, e.g. SAMPLE_001")
    ap.add_argument("--expr_id", required=True,
                    help="Sample ID in expression matrix, e.g. SAMPLE_001")
    ap.add_argument("--tss_bed", default=E.TSS_BED)
    ap.add_argument("--out_dir", default=E.PRED_PERSONAL)
    ap.add_argument("--genes_per_batch", type=int, default=2,
                    help="Genes per batch; each contributes 2 haplotypes x 400 seqs")
    ap.add_argument("--batch_size", type=int, default=512)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    outfile = f"{args.out_dir}/{args.expr_id}.parquet"
    if os.path.exists(outfile):
        print(f"Already exists, skipping: {outfile}")
        return

    print(f"Sample: {args.vcf_id} / {args.expr_id}")
    print(f"Tissue: {E.TISSUE}")

    tss_df = E.load_tss_bed(args.tss_bed)
    print(f"Genes to predict: {len(tss_df)}")

    print("Loading Beluga + XGBoost...")
    beluga, bst, device = E.load_models()
    print(f"  device: {device}")
    if str(device) == "cpu":
        print("  WARNING: no GPU; this will not finish in a reasonable time.")

    results = {}
    t0 = time.time()
    done_total = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        for chrom in E.chrom_order(tss_df):
            genes = tss_df[tss_df["chrom"] == f"chr{chrom}"]
            if len(genes) == 0:
                continue

            print(f"  chr{chrom}: {len(genes)} genes, building personal sequences...")
            try:
                hap1_fa, hap2_fa = E.build_personal_fasta(args.vcf_id, chrom, tmpdir)
            except RuntimeError as e:
                print(f"  chr{chrom}: bcftools failed, skipping: {e}")
                for _, row in genes.iterrows():
                    results[row.gene_id] = np.nan
                continue

            fa1 = pyfaidx.Fasta(hap1_fa)
            fa2 = pyfaidx.Fasta(hap2_fa)
            rows = list(genes.itertuples(index=False))

            for start in range(0, len(rows), args.genes_per_batch):
                batch = rows[start:start + args.genes_per_batch]
                seqs, strands, ok = [], [], []
                for row in batch:
                    try:
                        seqs.append(E.get_seq(fa1, f"chr{chrom}", row.tss, row.strand))
                        seqs.append(E.get_seq(fa2, f"chr{chrom}", row.tss, row.strand))
                        strands.extend([row.strand, row.strand])
                        ok.append(row)
                    except Exception as e:
                        print(f"    WARNING: {row.gene_name} failed: {e}")
                        results[row.gene_id] = np.nan
                if not seqs:
                    continue

                p = E.expecto_predict(beluga, bst, device, seqs, strands,
                                      args.batch_size)
                for i, row in enumerate(ok):
                    # haplotype averaging at the PREDICTION level, as in Enformer
                    results[row.gene_id] = (p[2 * i] + p[2 * i + 1]) / 2

                done_total += len(batch)
                if done_total % 200 < args.genes_per_batch:
                    el = time.time() - t0
                    rate = done_total / el
                    print(f"    {done_total}/{len(tss_df)} genes  {rate:.2f} genes/s  "
                          f"eta {(len(tss_df) - done_total) / max(rate, 1e-9) / 60:.0f} min")

            fa1.close()
            fa2.close()
            E.cleanup_fastas(hap1_fa, hap2_fa)
            print(f"  chr{chrom}: done ({len(genes)} genes)")

    out = pd.DataFrame({E.COLUMN: pd.Series(results)}).astype("float32")
    out.index.name = "gene_id"
    out.to_parquet(outfile)
    print(f"Saved: {outfile}  ({len(out)} genes, "
          f"{int(out[E.COLUMN].isna().sum())} NaN)")


if __name__ == "__main__":
    main()