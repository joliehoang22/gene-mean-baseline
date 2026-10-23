"""
predict_personal_blood.py: Enformer, personalized input, blood panel.

For one sample, this script builds both haplotype FASTAs with bcftools
consensus, runs Enformer on the 196,608 bp window centred on each TSS, averages
the 10 bins (1,280 bp) around the TSS for the 4 whole-blood CAGE tracks, and then
averages the two haplotypes. Writes <out_dir>/<expr_id>.parquet (gene_id x 4 tracks).
Driven by model_predictions/lsf/run_personal.sh.
"""
import torch
import pyfaidx
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import argparse
from enformer_pytorch import Enformer

# ── constants ────────────────────────────────────────────────────────────────
WINDOW       = 196_608
HALF         = WINDOW // 2
TSS_BIN      = 448       # center bin (896 bins total, 0-indexed)
BIN_WINDOW   = 10        # average ±5 bins around TSS
BLOOD_TRACKS = [4805, 4806, 5067, 4950]
BLOOD_TRACK_NAMES = [
    'CAGE_whole_blood_ribopure_donation1',
    'CAGE_whole_blood_ribopure_donation2',
    'CAGE_whole_blood_ribopure_donation3',
    'CAGE_blood_adult_pool1',
]

# ── paths ────────────────────────────────────────────────────────────────────
REF        = os.environ.get("HG38_FASTA", os.path.join(os.environ["PROJECT_DIR"], "ref_genome/hg38.fa"))
VCF_DIR    = os.environ.get("GFM_VCF_DIR", os.path.join(os.environ["PROJECT_DIR"], "phased_vcfs"))
TSS_BED    = os.path.join(os.environ["PROJECT_DIR"], "tss_hg38_mane_expressed_blood.bed")
OUT_DIR = os.path.join(os.environ["PROJECT_DIR"], "predictions/enformer_expanded_personal_blood")
MODEL_PATH = os.environ.get("ENFORMER_WEIGHTS", "EleutherAI/enformer-official-rough")

# ── helpers ──────────────────────────────────────────────────────────────────
def one_hot_encode(seq):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    arr = np.zeros((len(seq), 4), dtype=np.float32)
    for i, base in enumerate(seq.upper()):
        if base in mapping:
            arr[i, mapping[base]] = 1.0
    return arr

def reverse_complement(seq):
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]

def get_seq(fa, chrom, tss, strand):
    tss   = int(tss)
    start = max(0, tss - HALF)
    end   = start + WINDOW
    seq   = fa[chrom][start:end].seq
    if len(seq) < WINDOW:
        seq = seq + 'N' * (WINDOW - len(seq))
    if strand == '-':
        seq = reverse_complement(seq)
    return seq

def build_personal_fasta(vcf_sample_id, chrom, tmpdir):
    """Build personalized haplotype FASTAs for one chromosome."""
    vcf     = f"{VCF_DIR}/chr{chrom}.phased.vcf.gz"
    hap1_fa = os.path.join(tmpdir, f"hap1_chr{chrom}.fa")
    hap2_fa = os.path.join(tmpdir, f"hap2_chr{chrom}.fa")

    for hap, outfa in [("1", hap1_fa), ("2", hap2_fa)]:
        faidx_cmd     = ["samtools", "faidx", REF, f"chr{chrom}"]
        consensus_cmd = [
            "bcftools", "consensus",
            "-s", vcf_sample_id,
            "-H", hap,
            "-o", outfa,
            vcf
        ]
        faidx  = subprocess.Popen(faidx_cmd, stdout=subprocess.PIPE)
        result = subprocess.run(consensus_cmd, stdin=faidx.stdout,
                                capture_output=True)
        faidx.stdout.close()
        faidx.wait()

        if result.returncode != 0:
            raise RuntimeError(
                f"bcftools consensus failed for {vcf_sample_id} "
                f"hap{hap} chr{chrom}:\n{result.stderr.decode()}"
            )

    return hap1_fa, hap2_fa

def predict_gene(model, fa, chrom, tss, strand, device):
    """Run Enformer on one sequence, return 4 blood CAGE track values."""
    seq  = get_seq(fa, chrom, tss, strand)
    x    = torch.tensor(one_hot_encode(seq)).unsqueeze(0).to(device)
    with torch.no_grad():
        preds = model(x)['human']   # (1, 896, 5313)
    lo  = TSS_BIN - BIN_WINDOW // 2
    hi  = TSS_BIN + BIN_WINDOW // 2
    sig = preds[0, lo:hi, :].mean(dim=0)        # (5313,)
    return sig[BLOOD_TRACKS].cpu().numpy()       # (4,)

# ── main ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_id",  required=True,
                        help="Sample ID in VCF e.g. SAMPLE_001")
    parser.add_argument("--expr_id", required=True,
                        help="Sample ID in expression matrix e.g. SAMPLE_001")
    parser.add_argument("--tss_bed", default=None,
                        help="override the default gene panel")
    parser.add_argument("--out_dir", default=None,
                        help="override the default output directory")
    args = parser.parse_args()
    # Optional overrides so the same script can run any gene panel.
    global TSS_BED, OUT_DIR
    if args.tss_bed: TSS_BED = args.tss_bed
    if args.out_dir: OUT_DIR = args.out_dir


    os.makedirs(OUT_DIR, exist_ok=True)
    outfile = f"{OUT_DIR}/{args.expr_id}.parquet"

    if os.path.exists(outfile):
        print(f"Already exists, skipping: {outfile}")
        return

    print(f"Sample: {args.vcf_id} / {args.expr_id}")

    # load TSS bed
    tss_df = pd.read_csv(TSS_BED, sep='\t', header=None,
                         names=['chrom', 'tss', 'strand', 'gene_id', 'gene_name'])
    print(f"Genes to predict: {len(tss_df)}")

    # load Enformer
    print("Loading Enformer...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    model  = Enformer.from_pretrained(MODEL_PATH)
    model  = model.to(device)
    model.eval()

    results = {}  # gene_id -> array of shape (4,)

    # get unique chromosomes in order
    chroms = sorted(
        tss_df['chrom'].str.replace('chr', '').unique(),
        key=lambda x: int(x) if x.isdigit() else 99
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        for chrom in chroms:
            genes_on_chrom = tss_df[tss_df['chrom'] == f'chr{chrom}']
            if len(genes_on_chrom) == 0:
                continue

            print(f"  chr{chrom}: {len(genes_on_chrom)} genes: "
                  f"building personal sequences...")

            try:
                hap1_fa, hap2_fa = build_personal_fasta(
                    args.vcf_id, chrom, tmpdir
                )
            except RuntimeError as e:
                print(f"  chr{chrom}: bcftools failed, skipping: {e}")
                for _, row in genes_on_chrom.iterrows():
                    results[row.gene_id] = np.full(len(BLOOD_TRACKS), np.nan)
                continue

            fa1 = pyfaidx.Fasta(hap1_fa)
            fa2 = pyfaidx.Fasta(hap2_fa)

            for idx, (_, row) in enumerate(genes_on_chrom.iterrows()):
                try:
                    s1 = predict_gene(model, fa1, f'chr{chrom}',
                                      row.tss, row.strand, device)
                    s2 = predict_gene(model, fa2, f'chr{chrom}',
                                      row.tss, row.strand, device)
                    results[row.gene_id] = (s1 + s2) / 2  # shape (4,)
                except Exception as e:
                    print(f"    WARNING: {row.gene_name} failed: {e}")
                    results[row.gene_id] = np.full(len(BLOOD_TRACKS), np.nan)

                if (idx + 1) % 50 == 0:
                    print(f"    {idx + 1}/{len(genes_on_chrom)} genes done")

            fa1.close()
            fa2.close()

            # clean up chr FASTAs to save disk space
            os.remove(hap1_fa)
            os.remove(hap2_fa)
            if os.path.exists(hap1_fa + '.fai'):
                os.remove(hap1_fa + '.fai')
            if os.path.exists(hap2_fa + '.fai'):
                os.remove(hap2_fa + '.fai')

            print(f"  chr{chrom}: done ({len(genes_on_chrom)} genes)")

    # save, shape (n_genes, 4)
    out = pd.DataFrame(results, index=BLOOD_TRACK_NAMES).T
    out.index.name = 'gene_id'
    out.to_parquet(outfile)
    print(f"Saved: {outfile}  ({len(results)} genes, shape {out.shape})")

if __name__ == '__main__':
    main()