#!/usr/bin/env python
"""
Enformer on the reference genome. Uses the same model, window and tracks as
predict_personal_*.py; the TSS bed and output directory are given as
arguments. Driven by model_predictions/lsf/run_reference_enformer.sh.

  python predict_reference.py --tissue brain \
      --tss_bed $PROJECT_DIR/tss_hg38_mane_expressed_brain.bed \
      --out_dir $PROJECT_DIR/predictions/enformer_reference_expanded_brain
"""
import argparse, glob, os
import numpy as np
import pandas as pd
import pyfaidx
import torch
from enformer_pytorch import Enformer

WINDOW, HALF = 196_608, 196_608 // 2
TSS_BIN, BIN_WINDOW = 448, 10

BRAIN_TRACKS = [4708, 5094, 5261, 4680, 4980, 4769, 5112]
BRAIN_TRACK_NAMES = [
    'CAGE_frontal_lobe_adult', 'CAGE_medial_frontal_gyrus_adult_1',
    'CAGE_medial_frontal_gyrus_adult_2', 'CAGE_brain_adult_pool1',
    'CAGE_brain_adult', 'CAGE_astrocyte_cerebral_cortex', 'CAGE_neurons']
BLOOD_TRACKS = [4805, 4806, 5067, 4950]
BLOOD_TRACK_NAMES = [
    'CAGE_whole_blood_ribopure_donation1', 'CAGE_whole_blood_ribopure_donation2',
    'CAGE_whole_blood_ribopure_donation3', 'CAGE_blood_adult_pool1']

BASE       = os.environ["PROJECT_DIR"]
REF        = os.environ.get("HG38_FASTA", f"{BASE}/ref_genome/hg38.fa")
MODEL_PATH = os.environ.get("ENFORMER_WEIGHTS", "EleutherAI/enformer-official-rough")


def one_hot_encode(seq):
    m = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    a = np.zeros((len(seq), 4), dtype=np.float32)
    for i, b in enumerate(seq.upper()):
        if b in m:
            a[i, m[b]] = 1.0
    return a


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))[::-1]


def get_seq(fa, chrom, tss, strand):
    tss = int(tss)
    start = max(0, tss - HALF)
    seq = fa[chrom][start:start + WINDOW].seq
    if len(seq) < WINDOW:
        seq = seq + 'N' * (WINDOW - len(seq))
    return reverse_complement(seq) if strand == '-' else seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True, choices=["brain", "blood"])
    ap.add_argument("--tss_bed", required=True)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    idx, names = ((BRAIN_TRACKS, BRAIN_TRACK_NAMES) if args.tissue == "brain"
                  else (BLOOD_TRACKS, BLOOD_TRACK_NAMES))
    assert len(idx) == len(names)
    print(f"Tissue: {args.tissue}")
    for i, n in zip(idx, names):
        print(f"  {i:5d}  {n}")

    os.makedirs(args.out_dir, exist_ok=True)
    outfile = f"{args.out_dir}/reference_genome.parquet"
    if os.path.exists(outfile):
        print(f"Already exists: {outfile}")
        return

    bed = pd.read_csv(args.tss_bed, sep='\t', header=None,
                      names=['chrom', 'tss', 'strand', 'gene_id', 'gene_name'])
    print(f"Genes to predict: {len(bed)}", flush=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print("Device:", device)
    model = Enformer.from_pretrained(MODEL_PATH).to(device)
    model.eval()
    fa = pyfaidx.Fasta(REF)

    lo, hi = TSS_BIN - BIN_WINDOW // 2, TSS_BIN + BIN_WINDOW // 2
    res = {}
    for k, r in enumerate(bed.itertuples(index=False)):
        try:
            x = torch.tensor(one_hot_encode(
                get_seq(fa, r.chrom, r.tss, r.strand))).unsqueeze(0).to(device)
            with torch.no_grad():
                p = model(x)['human']
            res[r.gene_id] = p[0, lo:hi, :].mean(dim=0)[idx].cpu().numpy()
        except Exception as e:
            print(f"  WARNING {r.gene_name}: {e}")
            res[r.gene_id] = np.full(len(idx), np.nan)
        if (k + 1) % 500 == 0:
            print(f"  {k + 1}/{len(bed)}", flush=True)

    out = pd.DataFrame(res, index=names).T
    out.index.name = 'gene_id'
    out.to_parquet(outfile)
    print(f"Saved: {outfile}  {out.shape}")


if __name__ == '__main__':
    main()
