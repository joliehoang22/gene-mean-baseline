"""
xpresso_lib_brain.py: shared code for the Xpresso brain scripts. No other project imports.

Sequence handling follows enformer/predict_personal_brain.py (same TSS bed, same
bcftools variant application, same N -> zeros one-hot), so all models get the
same variant-inserted sequence.

Window (unlike Enformer/Basenji2, Xpresso's window is asymmetric):
  window = 7,000 upstream + 3,500 downstream = 10,500 bp
  After orienting to the gene's direction, the TSS sits at offset 6,999:
  6,999 bp upstream, the TSS base, then 3,500 bp downstream.
  This matches the ni-lab offset convention (their `upstream = args.upstream - 1`).
  xpresso_predict.py requires a length of exactly 10,500.

Output: one value per gene (predicted steady-state mRNA abundance, on the
z-scored scale the model was trained on). There are no bins or tracks.

Two arms, both computed in one pass by the predict_* scripts:
  seqonly  half-life features set to zero. setup_training_files.py z-scored
           these columns (preprocessing.scale), so zero is each feature's mean.
           This is sequence-only Xpresso with half-life held at the average.
  full     the real per-gene half-life features, as published.
Half-life features are the same for every individual, so they can change the
ranking across genes but not the ranking of individuals within a gene.
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

REPO = os.path.join(os.environ["PEB_REPO"], "xpresso")
MODEL = f"{REPO}/xpresso/pretrained_models/humanMedian_trainepoch.11-0.426.h5"

HALFLIFE_DIR = f"{BASE}/Xpresso/halflife"
HALFLIFE_PARQUET = f"{HALFLIFE_DIR}/halflife_features.parquet"

PRED_SEQONLY = f"{BASE}/predictions/xpresso_seqonly_mane"
PRED_FULL = f"{BASE}/predictions/xpresso_full_mane"
REF_SEQONLY = f"{BASE}/predictions/xpresso_reference/seqonly.parquet"
REF_FULL = f"{BASE}/predictions/xpresso_reference/full.parquet"

ENFORMER_REF = f"{BASE}/predictions/enformer_reference_expanded_brain/reference_genome.parquet"
BASENJI2_REF = f"{BASE}/predictions/basenji2_reference_expanded_brain/reference_genome.parquet"

# ── geometry ─────────────────────────────────────────────────────────────────
UPSTREAM = 7000
DOWNSTREAM = 3500
WINDOW = UPSTREAM + DOWNSTREAM      # 10,500, checked by the model wrapper
COLUMN = "xpresso_humanMedian"

# ── one-hot (identical behaviour to the Enformer script) ─────────────────────
_LUT = np.zeros((256, 4), dtype=np.float32)
for _i, _b in enumerate("ACGT"):
    _LUT[ord(_b), _i] = 1.0
    _LUT[ord(_b.lower()), _i] = 1.0

_RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def one_hot_encode(seq):
    return _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def reverse_complement(seq):
    return seq.translate(_RC_TABLE)[::-1]


def get_seq(fa, chrom, tss, strand):
    """
    Extract Xpresso's asymmetric 10,500 bp window, oriented to the gene.

    plus  : [tss-6999, tss+3501)   -> TSS at offset 6999
    minus : [tss-3500, tss+7000)   -> reverse-complemented, TSS at offset 6999
    """
    tss = int(tss)
    if strand == "-":
        end = tss + UPSTREAM
        start = end - WINDOW
    else:
        start = tss - (UPSTREAM - 1)
        end = start + WINDOW

    start = max(0, start)
    seq = fa[chrom][start:start + WINDOW].seq
    if len(seq) < WINDOW:
        seq = seq + "N" * (WINDOW - len(seq))
    if strand == "-":
        seq = reverse_complement(seq)
    return seq


def tss_is_truncated(tss, strand):
    """True if the window ran off the chromosome start (TSS not at offset 6999)."""
    tss = int(tss)
    start = tss + UPSTREAM - WINDOW if strand == "-" else tss - (UPSTREAM - 1)
    return start < 0


# ── personalized haplotype FASTAs (same as the Enformer script) ──────────────
def build_personal_fasta(vcf_sample_id, chrom, tmpdir):
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


# ── half-life features ───────────────────────────────────────────────────────
def load_halflife(gene_ids, n_features):
    """
    Build an (n_genes, n_features) array aligned to gene_ids.

    Returns (array, n_matched). Genes with no entry get NaN rows. Callers
    should drop these from the `full` arm instead of filling in zeros, because
    zero means "average gene" and would look like a real value.
    """
    out = np.full((len(gene_ids), n_features), np.nan, dtype="float32")
    if not os.path.exists(HALFLIFE_PARQUET):
        return out, 0

    hl = pd.read_parquet(HALFLIFE_PARQUET).set_index("gene_id")
    if hl.shape[1] != n_features:
        raise ValueError(
            f"half-life table has {hl.shape[1]} columns but the model expects "
            f"{n_features}. Resolve before running the `full` arm; see the "
            f"8-vs-6 note in 01_setup.sh."
        )

    # match on unversioned Ensembl IDs (ENSG00000123456.7 -> ENSG00000123456)
    keys = pd.Index([str(g).split(".")[0] for g in gene_ids])
    hl.index = pd.Index([str(g).split(".")[0] for g in hl.index])
    hl = hl[~hl.index.duplicated()]

    hit = keys.isin(hl.index)
    out[hit] = hl.reindex(keys[hit]).to_numpy(dtype="float32")
    return out, int(hit.sum())


# ── model ────────────────────────────────────────────────────────────────────
def load_model():
    """
    Load the pretrained humanMedian Keras model.

    compile=False skips rebuilding the optimizer and loss, which fails for
    TF1-era .h5 files under newer TF. Prediction does not need them.
    Returns (model, n_halflife_features).
    """
    from tensorflow import keras

    model = keras.models.load_model(MODEL, compile=False)

    shapes = model.input_shape
    if not isinstance(shapes, list) or len(shapes) != 2:
        raise ValueError(f"expected 2 inputs (promoter, halflife), got {shapes}")
    seq_shape, hl_shape = shapes
    if seq_shape[1] != WINDOW:
        raise ValueError(f"model wants {seq_shape[1]}bp, geometry gives {WINDOW}")
    return model, int(hl_shape[1])


def predict_batch(model, seqs, halflife):
    """seqs: list of 10,500bp strings. halflife: (B, n_features). -> (B,)"""
    x = np.stack([one_hot_encode(s) for s in seqs])
    preds = model.predict([x, halflife], batch_size=len(seqs), verbose=0)
    return np.asarray(preds).flatten()