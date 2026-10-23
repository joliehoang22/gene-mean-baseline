"""
expecto_lib_brain.py: shared code for the ExPecto brain scripts. No other project imports.

Variant application and the TSS bed match enformer/predict_personal_brain.py, so
all four models get the same variant-inserted sequence.

ExPecto has three stages:
  1. Beluga      ~145M-param CNN, 2,000 bp window -> 2,002 chromatin marks
  2. spatial     Beluga is run at 200 shifts (-20 kb to +20 kb, step 200), each
                 also reverse-complemented and averaged. The 200 x 2002 matrix
                 is reduced with 10 exponential-decay basis functions (5 decay
                 rates x upstream/downstream) to 20,020 features. One zero
                 column per basis function, kept from the original code, gives
                 20,030 features
  3. XGBoost     per-tissue regressor -> predicted log expression (one value)

ExPecto uses some conventions that differ from the other models. It was
trained with them, so they are kept:
  * The window is not reverse-complemented for minus-strand genes. ExPecto keeps
    genomic orientation and flips the shift direction (shift * strand) instead.
  * One-hot order is A,G,C,T (DeepSEA convention), not ACGT. With this order the
    [:, ::-1, ::-1] flip gives the correct reverse complement. An ACGT lookup
    here would swap G and C without raising an error.
  * Both orientations are always used: encode_seqs appends the reverse
    complement of every window and the two predictions are averaged.

Window size: shifts reach +/-20,000 and each window extends +/-1,000, so the
sequence has to cover +/-21,000. A 42,000 bp sequence is one base short on the
minus strand, so we use 43,000 (HALF = 21,500). The TSS sits at
offset len//2, matching the ni-lab code's `tss_i = len(sample_seq) // 2`.
"""

import os
import subprocess

import numpy as np
import pandas as pd

# ── paths ────────────────────────────────────────────────────────────────────
BASE = os.environ["PROJECT_DIR"]

REF = os.environ.get("HG38_FASTA", f"{BASE}/ref_genome/hg38.fa")
VCF_DIR = os.environ.get("GFM_VCF_DIR", f"{BASE}/phased_vcfs")
TSS_BED = f"{BASE}/tss_hg38_mane_expressed_brain.bed"

REPO = os.path.join(os.environ["PEB_REPO"], "expecto")        # Beluga class
OFFICIAL = os.environ["EXPECTO_REPO"]                       # FunctionLab/ExPecto clone
BELUGA_WEIGHTS = f"{OFFICIAL}/resources/deepsea.beluga.pth"

TISSUE = "Brain_Frontal_Cortex_BA9"
XGB_MODEL = (f"{OFFICIAL}/models/allhistones2000.1.fixed.all.pseudocount0.0001."
             f"lambda100.round100.basescore2.{TISSUE}.save")

PRED_PERSONAL = f"{BASE}/predictions/expecto_expanded_personal_brain"
PRED_REFERENCE = f"{BASE}/predictions/expecto_reference_expanded_brain"

ENFORMER_REF = f"{BASE}/predictions/enformer_reference_expanded_brain/reference_genome.parquet"
BASENJI2_REF = f"{BASE}/predictions/basenji2_reference_expanded_brain/reference_genome.parquet"
XPRESSO_REF = f"{BASE}/predictions/xpresso_reference_expanded_brain/seqonly.parquet"

# ── geometry ─────────────────────────────────────────────────────────────────
HALF = 21_500
WINDOW = 2 * HALF                  # 43,000
WINDOWSIZE = 2000                  # Beluga input
SHIFTS = np.array(list(range(-20000, 20000, 200)))   # 200 shifts
N_MARKS = 2002
COLUMN = f"expecto_{TISSUE}"

# ── one-hot: A,G,C,T order (DeepSEA), non-ACGT -> all zeros ──────────────────
_LUT = np.zeros((256, 4), dtype=np.float32)
for _i, _b in enumerate("AGCT"):
    _LUT[ord(_b), _i] = 1.0
    _LUT[ord(_b.lower()), _i] = 1.0


def encode_seqs(seqs):
    """
    (n, 2000) uint8 char array -> (2n, 4, 2000) float32.

    Vectorized equivalent of the repo's encodeSeqs triple loop. The second half
    of the output is the reverse complement of the first: reversing both the
    channel axis (A,G,C,T -> T,C,G,A) and the position axis complements the
    sequence, which only works because of the A,G,C,T ordering above.
    """
    fwd = _LUT[seqs].transpose(0, 2, 1)          # (n, 4, 2000)
    rc = fwd[:, ::-1, ::-1]
    return np.ascontiguousarray(np.concatenate([fwd, rc], axis=0))


def get_shift_windows(sample_seq, strand):
    """
    Cut the 200 Beluga windows out of a 43,000bp sequence.

    Genomic orientation is preserved; the strand only flips the direction the
    shifts run in, as in the ni-lab code.
    """
    tss_i = len(sample_seq) // 2
    s = 1 if strand == "+" else -1
    arr = np.frombuffer(sample_seq.encode("ascii"), dtype=np.uint8)

    out = np.empty((len(SHIFTS), WINDOWSIZE), dtype=np.uint8)
    for i, shift in enumerate(SHIFTS):
        c = tss_i + shift * s
        w = arr[c - (WINDOWSIZE // 2 - 1): c + (WINDOWSIZE // 2) + 1]
        if len(w) != WINDOWSIZE:
            raise ValueError(f"window {i} has length {len(w)}, expected {WINDOWSIZE}")
        out[i] = w
    return out


# ── spatial transform (10 exponential-decay basis functions) ─────────────────
_POS_WEIGHTS = np.vstack([
    np.exp(-0.01 * np.abs(SHIFTS) / 200) * (SHIFTS <= 0),
    np.exp(-0.02 * np.abs(SHIFTS) / 200) * (SHIFTS <= 0),
    np.exp(-0.05 * np.abs(SHIFTS) / 200) * (SHIFTS <= 0),
    np.exp(-0.10 * np.abs(SHIFTS) / 200) * (SHIFTS <= 0),
    np.exp(-0.20 * np.abs(SHIFTS) / 200) * (SHIFTS <= 0),
    np.exp(-0.01 * np.abs(SHIFTS) / 200) * (SHIFTS >= 0),
    np.exp(-0.02 * np.abs(SHIFTS) / 200) * (SHIFTS >= 0),
    np.exp(-0.05 * np.abs(SHIFTS) / 200) * (SHIFTS >= 0),
    np.exp(-0.10 * np.abs(SHIFTS) / 200) * (SHIFTS >= 0),
    np.exp(-0.20 * np.abs(SHIFTS) / 200) * (SHIFTS >= 0)])


def spatial_features(beluga_preds):
    """(n_genes, 200, 2002) -> (n_genes, 20030) for XGBoost."""
    f = np.sum(_POS_WEIGHTS[None, :, :, None] * beluga_preds[:, None, :, :], axis=2)
    f = f.reshape(-1, 10, N_MARKS)
    # one leading zero column per basis function, as in the original ExPecto code
    f = np.concatenate([np.zeros((f.shape[0], 10, 1)), f], axis=2)
    return f.reshape(f.shape[0], 10 * (N_MARKS + 1))


# ── models ───────────────────────────────────────────────────────────────────
def load_models(device=None):
    """Returns (beluga, booster, device)."""
    import sys

    import torch
    import xgboost as xgb

    sys.path.insert(0, REPO)
    from Beluga import Beluga  # noqa: E402

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    beluga = Beluga()
    beluga.load_state_dict(torch.load(BELUGA_WEIGHTS, map_location=device))
    beluga.eval().to(device)

    bst = xgb.Booster()
    bst.load_model(XGB_MODEL)
    return beluga, bst, device


def beluga_predict(beluga, device, windows, batch_size=512):
    """
    windows: (n_genes, 200, 2000) uint8 -> (n_genes, 200, 2002) float32.

    Forward and reverse-complement predictions are averaged, as in ExPecto.
    The original code slices batches as [i*bs:(i+1)*bs] while `i` already steps
    by bs, which only works when everything fits in one batch. Here each batch
    is sliced from its start index.
    """
    import torch

    n_genes = windows.shape[0]
    flat = windows.reshape(-1, WINDOWSIZE)
    enc = encode_seqs(flat)                      # (2 * n_genes * 200, 4, 2000)

    out = np.zeros((enc.shape[0], N_MARKS), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, enc.shape[0], batch_size):
            batch = torch.from_numpy(enc[i:i + batch_size]).to(device).unsqueeze(2)
            out[i:i + batch_size] = beluga.forward(batch).cpu().numpy()

    half = out.shape[0] // 2
    out = (out[:half] + out[half:]) / 2          # average fwd / revcomp
    return out.reshape(n_genes, len(SHIFTS), N_MARKS)


def expecto_predict(beluga, bst, device, seqs, strands, batch_size=512):
    """seqs: list of 43,000bp strings, strands: list of '+'/'-'. -> (n,)"""
    import xgboost as xgb

    windows = np.stack([get_shift_windows(s, st) for s, st in zip(seqs, strands)])
    preds = beluga_predict(beluga, device, windows, batch_size)
    return bst.predict(xgb.DMatrix(spatial_features(preds)))


# ── sequence extraction (matches the Enformer script) ────────────────────────
def get_seq(fa, chrom, tss, strand):
    """
    43,000 bp centered on the TSS, in genomic orientation.

    Minus-strand genes are not reverse-complemented here, because ExPecto
    handles strand by flipping the shift direction. `strand` is accepted only so
    the signature matches get_seq in the other models.
    """
    tss = int(tss)
    start = max(0, tss - HALF)
    seq = fa[chrom][start:start + WINDOW].seq
    if len(seq) < WINDOW:
        seq = seq + "N" * (WINDOW - len(seq))
    return seq.upper()


def tss_is_truncated(tss):
    return int(tss) - HALF < 0


def build_personal_fasta(vcf_sample_id, chrom, tmpdir):
    vcf = f"{VCF_DIR}/chr{chrom}.phased.vcf.gz"
    hap1_fa = os.path.join(tmpdir, f"hap1_chr{chrom}.fa")
    hap2_fa = os.path.join(tmpdir, f"hap2_chr{chrom}.fa")

    for hap, outfa in [("1", hap1_fa), ("2", hap2_fa)]:
        faidx = subprocess.Popen(
            ["samtools", "faidx", REF, f"chr{chrom}"], stdout=subprocess.PIPE)
        result = subprocess.run(
            ["bcftools", "consensus", "-s", vcf_sample_id, "-H", hap,
             "-o", outfa, vcf],
            stdin=faidx.stdout, capture_output=True)
        faidx.stdout.close()
        faidx.wait()
        if result.returncode != 0:
            raise RuntimeError(
                f"bcftools consensus failed for {vcf_sample_id} hap{hap} "
                f"chr{chrom}:\n{result.stderr.decode()}")
    return hap1_fa, hap2_fa


def cleanup_fastas(*paths):
    for p in paths:
        for f in (p, p + ".fai"):
            if os.path.exists(f):
                os.remove(f)


# ── TSS bed ──────────────────────────────────────────────────────────────────
def load_tss_bed(path=TSS_BED):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "tss", "strand", "gene_id", "gene_name"])
    dup = int(df["gene_id"].duplicated().sum())
    if dup:
        raise ValueError(f"{dup} duplicate gene_id in {path}; resolve before running")
    return df


def chrom_order(df):
    return sorted(df["chrom"].str.replace("chr", "", regex=False).unique(),
                  key=lambda x: int(x) if x.isdigit() else 99)