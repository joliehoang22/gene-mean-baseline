"""
basenji2_lib_brain.py: shared code for the Basenji2 brain scripts. No other project imports.

The sequence functions below behave the same as those in
enformer/predict_personal_brain.py, so Basenji2 gets exactly the same
variant-inserted sequence, TSS position and strand handling. A change to one
should be made in both.

Geometry (checked against resources/params_human.json):
  input  131,072 bp
  trunk  1024 bins -> Cropping1D(64) -> 896 bins x 128 bp = 114,688 bp
  head   5,313 softplus targets, the same target set as Enformer's human head,
         so the brain_tracks.txt indices are the same

Centering the TSS in 131,072 bp puts it at bin 448 of 896, the same bin as in
Enformer's 196,608 bp window. TSS_BIN and BIN_WINDOW therefore match the
Enformer script.
"""

import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

# ── paths ────────────────────────────────────────────────────────────────────
BASE = os.environ["PROJECT_DIR"]

REF = os.environ.get("HG38_FASTA", f"{BASE}/ref_genome/hg38.fa")
VCF_DIR = os.environ.get("GFM_VCF_DIR", f"{BASE}/phased_vcfs")
TSS_BED = f"{BASE}/tss_hg38_mane_expressed_brain.bed"
TARGETS = f"{BASE}/targets_human.txt"
BRAIN_TRACKS_FILE = f"{BASE}/brain_tracks.txt"

BASENJI_DIR = os.path.join(os.environ["PEB_REPO"], "basenji2")
PARAMS = f"{BASENJI_DIR}/resources/params_human.json"
WEIGHTS = f"{BASENJI_DIR}/resources/model_human.h5"

PRED_PERSONAL = f"{BASE}/predictions/basenji2_expanded_personal_brain"
PRED_REFERENCE = f"{BASE}/predictions/basenji2_reference_expanded_brain"
ENFORMER_REF = f"{BASE}/predictions/enformer_reference_expanded_brain/reference_genome.parquet"

# ── constants ────────────────────────────────────────────────────────────────
WINDOW = 131_072          # Basenji2 input length
TSS_BIN = 448             # center bin of 896, same as Enformer
BIN_WINDOW = 10           # average +/-5 bins around TSS, same as Enformer
N_BINS = 896
N_TARGETS = 5313

# ── sequence handling (mirrors enformer/predict_personal_brain.py) ─────────────────────
# Lookup table over all 256 byte values: ACGT (upper and lower) -> one-hot,
# everything else (N, IUPAC) -> all zeros. Same result as the Enformer script's
# loop, ~11x faster. Lowercase matters: hg38.fa is soft-masked.
_LUT = np.zeros((256, 4), dtype=np.float32)
for _i, _b in enumerate("ACGT"):
    _LUT[ord(_b), _i] = 1.0
    _LUT[ord(_b.lower()), _i] = 1.0

_RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def one_hot_encode(seq):
    return _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def reverse_complement(seq):
    return seq.translate(_RC_TABLE)[::-1]


def get_seq(fa, chrom, tss, strand, window=WINDOW):
    """Extract a window-bp sequence centered on tss from an open pyfaidx.Fasta."""
    tss = int(tss)
    half = window // 2
    start = max(0, tss - half)
    end = start + window
    seq = fa[chrom][start:end].seq
    if len(seq) < window:
        seq = seq + "N" * (window - len(seq))
    if strand == "-":
        seq = reverse_complement(seq)
    return seq


def tss_is_truncated(tss, window=WINDOW):
    """True if the window ran off the chromosome start (TSS not centered)."""
    return int(tss) - window // 2 < 0


def build_personal_fasta(vcf_sample_id, chrom, tmpdir):
    """
    Whole-chromosome haplotype FASTAs for one sample, written to tmpdir.

    Assumes a SNP-only phased VCF: no indels means no coordinate drift, so
    cropping a fixed window is safe.
    """
    vcf = f"{VCF_DIR}/chr{chrom}.phased.vcf.gz"
    hap1_fa = os.path.join(tmpdir, f"hap1_chr{chrom}.fa")
    hap2_fa = os.path.join(tmpdir, f"hap2_chr{chrom}.fa")

    for hap, outfa in [("1", hap1_fa), ("2", hap2_fa)]:
        faidx = subprocess.Popen(
            ["samtools", "faidx", REF, f"chr{chrom}"], stdout=subprocess.PIPE
        )
        result = subprocess.run(
            ["bcftools", "consensus", "-s", vcf_sample_id, "-H", hap,
             "-o", outfa, vcf],
            stdin=faidx.stdout, capture_output=True,
        )
        faidx.stdout.close()
        faidx.wait()
        if result.returncode != 0:
            raise RuntimeError(
                f"bcftools consensus failed for {vcf_sample_id} hap{hap} "
                f"chr{chrom}:\n{result.stderr.decode()}"
            )
    return hap1_fa, hap2_fa


def cleanup_fastas(*paths):
    for p in paths:
        for f in (p, p + ".fai"):
            if os.path.exists(f):
                os.remove(f)


# ── TSS bed ──────────────────────────────────────────────────────────────────
def load_tss_bed(path=TSS_BED):
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "tss", "strand", "gene_id", "gene_name"],
    )
    dup = int(df["gene_id"].duplicated().sum())
    if dup:
        raise ValueError(f"{dup} duplicate gene_id in {path}; resolve before running")
    return df


def chrom_order(df):
    return sorted(
        df["chrom"].str.replace("chr", "", regex=False).unique(),
        key=lambda x: int(x) if x.isdigit() else 99,
    )


# ── tracks ───────────────────────────────────────────────────────────────────
def load_brain_tracks():
    """Read the 7 brain CAGE indices from brain_tracks.txt."""
    with open(BRAIN_TRACKS_FILE) as f:
        txt = f.read().strip().strip("[]")
    idx = [int(x) for x in txt.replace(",", " ").split()]
    if not idx:
        raise ValueError(f"no track indices parsed from {BRAIN_TRACKS_FILE}")
    return idx


def track_names(indices):
    """
    Resolve names from targets_human.txt so Basenji2 parquet columns match
    Enformer's. Falls back to index labels if the file is missing or has an
    unexpected format; check the printout in predict_reference_*.py.
    """
    if not os.path.exists(TARGETS):
        return [f"track_{i}" for i in indices]
    try:
        tdf = pd.read_csv(TARGETS, sep="\t")
        col = next((c for c in ("description", "identifier", "file")
                    if c in tdf.columns), None)
        if col is None:
            return [f"track_{i}" for i in indices]
        return [str(tdf.iloc[i][col]) for i in indices]
    except Exception:
        return [f"track_{i}" for i in indices]


# ── model ────────────────────────────────────────────────────────────────────
def load_model():
    """Build Basenji2 from params and restore weights. Returns a keras model."""
    sys.path.insert(0, BASENJI_DIR)
    from basenji import seqnn  # noqa: E402

    with open(PARAMS) as f:
        params_model = json.load(f)["model"]

    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore(WEIGHTS)
    model = seqnn_model.model

    shape = model.output_shape
    assert shape[1] == N_BINS, f"expected {N_BINS} bins, got {shape[1]}"
    assert shape[2] == N_TARGETS, f"expected {N_TARGETS} targets, got {shape[2]}"
    return model


def predict_batch(model, seqs):
    """seqs: list of WINDOW-bp strings. Returns (B, 5313), TSS-bin averaged."""
    x = np.stack([one_hot_encode(s) for s in seqs])
    preds = np.asarray(model.predict_on_batch(x))      # (B, 896, 5313)
    lo = TSS_BIN - BIN_WINDOW // 2
    hi = TSS_BIN + BIN_WINDOW // 2
    return preds[:, lo:hi, :].mean(axis=1)