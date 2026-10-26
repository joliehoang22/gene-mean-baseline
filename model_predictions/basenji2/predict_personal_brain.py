#!/usr/bin/env python
"""
predict_personal_brain.py: Basenji2 version of enformer/predict_personal_brain.py.

Uses the same TSS bed, bcftools variant application, haplotype averaging (of
the predictions) and TSS-bin window as Enformer. The one difference is the
input length, 131,072 bp vs 196,608 bp for Enformer. This is part of each
model and is left as is.

Haplotype FASTAs are built per chromosome in a temp dir and deleted at the end
of the job, as in the Enformer script. About half the run time is bcftools,
not the GPU.
"""

import argparse
import os
import tempfile

import numpy as np
import pandas as pd
import pyfaidx

import basenji2_lib_brain as B


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf_id", required=True, help="Sample ID in VCF, e.g. SAMPLE_001")
    ap.add_argument("--expr_id", required=True,
                    help="Sample ID in expression matrix, e.g. SAMPLE_001")
    ap.add_argument("--tss_bed", default=B.TSS_BED)
    ap.add_argument("--out_dir", default=B.PRED_PERSONAL)
    ap.add_argument("--batch_size", type=int, default=8,
                    help="Sequences per forward pass (genes x haplotypes)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    outfile = f"{args.out_dir}/{args.expr_id}.parquet"
    if os.path.exists(outfile):
        print(f"Already exists, skipping: {outfile}")
        return

    print(f"Sample: {args.vcf_id} / {args.expr_id}")

    tracks = B.load_tracks()
    names = B.track_names(tracks)

    tss_df = B.load_tss_bed(args.tss_bed)
    print(f"Genes to predict: {len(tss_df)}")

    n_trunc = sum(B.tss_is_truncated(t) for t in tss_df["tss"])
    if n_trunc:
        print(f"NOTE: {n_trunc} gene(s) within {B.WINDOW//2}bp of chromosome "
              f"start; window not centered for those")

    print("Loading Basenji2...")
    model = B.load_model()

    results = {}
    genes_per_batch = max(1, args.batch_size // 2)   # 2 haplotypes per gene

    with tempfile.TemporaryDirectory() as tmpdir:
        for chrom in B.chrom_order(tss_df):
            genes = tss_df[tss_df["chrom"] == f"chr{chrom}"]
            if len(genes) == 0:
                continue

            print(f"  chr{chrom}: {len(genes)} genes, building personal sequences...")
            try:
                hap1_fa, hap2_fa = B.build_personal_fasta(args.vcf_id, chrom, tmpdir)
            except RuntimeError as e:
                print(f"  chr{chrom}: bcftools failed, skipping: {e}")
                for _, row in genes.iterrows():
                    results[row.gene_id] = np.full(len(tracks), np.nan)
                continue

            fa1 = pyfaidx.Fasta(hap1_fa)
            fa2 = pyfaidx.Fasta(hap2_fa)
            rows = list(genes.itertuples(index=False))

            for start in range(0, len(rows), genes_per_batch):
                batch = rows[start:start + genes_per_batch]
                seqs, ok = [], []
                for row in batch:
                    try:
                        seqs.append(B.get_seq(fa1, f"chr{chrom}", row.tss, row.strand))
                        seqs.append(B.get_seq(fa2, f"chr{chrom}", row.tss, row.strand))
                        ok.append(row)
                    except Exception as e:
                        print(f"    WARNING: {row.gene_name} failed: {e}")
                        results[row.gene_id] = np.full(len(tracks), np.nan)
                if not seqs:
                    continue

                sig = B.predict_batch(model, seqs)[:, tracks]
                for i, row in enumerate(ok):
                    # haplotype averaging at the PREDICTION level, as in Enformer
                    results[row.gene_id] = (sig[2 * i] + sig[2 * i + 1]) / 2

                done = start + len(batch)
                if done % 500 < genes_per_batch:
                    print(f"    {done}/{len(rows)} genes done")

            fa1.close()
            fa2.close()
            B.cleanup_fastas(hap1_fa, hap2_fa)
            print(f"  chr{chrom}: done ({len(genes)} genes)")

    out = pd.DataFrame(results, index=names).T
    out.index.name = "gene_id"
    out.to_parquet(outfile)
    print(f"Saved: {outfile}  ({len(results)} genes, shape {out.shape})")


if __name__ == "__main__":
    main()