#!/bin/bash
##############################################################
# run_personal.sh
#
# LSF array job, one element per sample. For each gene in the panel it
# builds the sample's two haplotype sequences around the TSS, runs one
# model on both, averages the two predictions, and writes one parquet per
# sample (genes x tracks).
#
# Submit from the repo root, one array per cohort / model / tissue:
#   mkdir -p logs
#   bsub -J "lbp_enformer_brain[1-225]" \
#        -env "all, COHORT=lbp, MODEL=enformer, TISSUE=brain" \
#        < model_predictions/lsf/run_personal.sh
#
#   COHORT  lbp | rosmap | covid     (ROSMAP is brain only, COVID-19 blood only)
#   MODEL   enformer | basenji2 | expecto | xpresso
#   TISSUE  brain | blood
#
# The array index picks a row of the cohort's sample table
# (column 1 = VCF sample ID, column 2 = RNA-seq sample ID).
##############################################################
#BSUB -P [ACCOUNT]
#BSUB -q [GPU_QUEUE]
#BSUB -W 60:00
#BSUB -n 4
#BSUB -gpu "num=1"
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=1000]"
#BSUB -o logs/%J_%I.out
#BSUB -e logs/%J_%I.err
#BSUB -L /bin/bash
set -uo pipefail
REPO_DIR=${LS_SUBCWD:-$(pwd)}
source "$REPO_DIR/config/paths.sh"
: "${COHORT:?set COHORT}" "${MODEL:?set MODEL}" "${TISSUE:?set TISSUE}"

# Sample table, header lines to skip, phased VCFs and output root per cohort.
case "$COHORT:$TISSUE" in
  lbp:brain)    SAMPLES=$PROJECT_DIR/sample_mapping.tsv;        OFFSET=1; VCFS=$VCF_DIR_LBP;    PRED=$PROJECT_DIR/predictions ;;
  lbp:blood)    SAMPLES=$PROJECT_DIR/samples_blood.tsv;         OFFSET=0; VCFS=$VCF_DIR_LBP;    PRED=$PROJECT_DIR/predictions ;;
  rosmap:brain) SAMPLES=$PROJECT_DIR/sample_mapping_rosmap.tsv; OFFSET=1; VCFS=$VCF_DIR_ROSMAP; PRED=$PROJECT_DIR/Replication/ROSMAP/predictions ;;
  covid:blood)  SAMPLES=$PROJECT_DIR/samples_covid.tsv;         OFFSET=0; VCFS=$VCF_DIR_COVID;  PRED=$PROJECT_DIR/Replication/COVID/predictions ;;
  *) echo "no $TISSUE samples for COHORT=$COHORT" >&2; exit 1 ;;
esac

# Enformer and ExPecto use PyTorch. Check the GPU first: on a node with a
# broken CUDA setup, torch falls back to CPU without an error and runs out of
# memory. The two TensorFlow models are not checked.
case "$MODEL" in
  enformer|expecto) PROBE=torch ;;
  basenji2|xpresso) PROBE=none ;;
  *) echo "unknown MODEL=$MODEL" >&2; exit 1 ;;
esac

DIR=$REPO_DIR/model_predictions/$MODEL
S=predict_personal_${TISSUE}.py
ENV=$CONDA_ENVS/$MODEL
BED=$PROJECT_DIR/tss_hg38_mane_expressed_${TISSUE}.bed
OUT=$PRED/${MODEL}_expanded_personal_${TISSUE}

cd "$DIR"
module purge; ml anaconda3; ml samtools bcftools   # Minerva modules; bcftools must be on PATH
unset PYTHONPATH PYTHONUSERBASE PYTHONHOME
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
export GFM_VCF_DIR=$VCFS

LINE=$(sed -n "$((LSB_JOBINDEX + OFFSET))p" "$SAMPLES")
VCF_ID=$(echo "$LINE" | cut -f1); EXPR_ID=$(echo "$LINE" | cut -f2)
[ -n "${VCF_ID:-}" ] && [ -n "${EXPR_ID:-}" ] || { echo "empty line ${LSB_JOBINDEX}" >&2; exit 1; }

echo "[$(date)] $COHORT $MODEL/$TISSUE idx=$LSB_JOBINDEX vcf=$VCF_ID expr=$EXPR_ID host=$(hostname)"

if [ "$PROBE" = torch ]; then
  "$ENV/bin/python" -c 'import torch,sys; sys.exit(0 if torch.cuda.is_available() else 1)' 2>/dev/null || {
    echo "no usable GPU on $(hostname) (idx=$LSB_JOBINDEX); stopping instead of running on CPU" >&2
    exit 42; }
fi

"$ENV/bin/python" "$S" --vcf_id "$VCF_ID" --expr_id "$EXPR_ID" \
    --tss_bed "$BED" --out_dir "$OUT"
RC=$?
echo "[$(date)] done rc=$RC"
exit $RC
