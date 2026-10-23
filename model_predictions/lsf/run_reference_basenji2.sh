#!/bin/bash
##############################################################
# run_reference_basenji2.sh
#
# Reference-genome arm for Basenji2: the same panel and windows as the
# personalized runs, with no variants inserted. One sequence per gene,
# shared by every sample. Array element 1 = brain, 2 = blood.
#
# Submit from the repo root:
#   mkdir -p logs && bsub < model_predictions/lsf/run_reference_basenji2.sh
##############################################################
#BSUB -J "basenji2_ref[1-2]"
#BSUB -P [ACCOUNT]
#BSUB -q [GPU_QUEUE]
#BSUB -W 24:00
#BSUB -n 4
#BSUB -gpu "num=1"
#BSUB -R "rusage[mem=40000]"
#BSUB -o logs/basenji2_ref.%I.out
#BSUB -e logs/basenji2_ref.%I.err
#BSUB -L /bin/bash
set -euo pipefail
REPO_DIR=${LS_SUBCWD:-$(pwd)}
source "$REPO_DIR/config/paths.sh"
TIS=(brain blood); T=${TIS[$((LSB_JOBINDEX-1))]}

cd "$REPO_DIR/model_predictions/basenji2"
module purge; ml anaconda3          # Minerva module
unset PYTHONPATH PYTHONUSERBASE PYTHONHOME; export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4

echo "[$(date)] basenji2 reference, $T"
"$CONDA_ENVS/basenji2/bin/python" predict_reference_${T}.py \
  --tss_bed "$PROJECT_DIR/tss_hg38_mane_expressed_${T}.bed" \
  --out_dir "$PROJECT_DIR/predictions/basenji2_reference_expanded_${T}"
echo "[$(date)] done"
