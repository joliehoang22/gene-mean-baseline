#!/bin/bash
##############################################################
# run_reference_expecto.sh
#
# Reference-genome arm for ExPecto: the same panel and windows as the
# personalized runs, with no variants inserted. One sequence per gene,
# shared by every sample. Array element 1 = brain, 2 = blood.
#
# Submit from the repo root:
#   mkdir -p logs && bsub < model_predictions/lsf/run_reference_expecto.sh
##############################################################
#BSUB -J "expecto_ref[1-2]"
#BSUB -P [ACCOUNT]
#BSUB -q [GPU_QUEUE]
#BSUB -W 72:00
#BSUB -n 4
#BSUB -gpu "num=1"
#BSUB -R "rusage[mem=24000]"
#BSUB -o logs/expecto_ref.%I.out
#BSUB -e logs/expecto_ref.%I.err
#BSUB -L /bin/bash
set -euo pipefail
REPO_DIR=${LS_SUBCWD:-$(pwd)}
source "$REPO_DIR/config/paths.sh"
TIS=(brain blood); T=${TIS[$((LSB_JOBINDEX-1))]}

cd "$REPO_DIR/model_predictions/expecto"
module purge; ml anaconda3          # Minerva module
unset PYTHONPATH PYTHONUSERBASE PYTHONHOME; export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4

# The slowest of the four: ~13,000 genes x 400 Beluga windows (fwd + revcomp).
# Exit if there is no GPU, since on CPU this takes days.
"$CONDA_ENVS/expecto/bin/python" - <<'PY' || exit 1
import sys, torch
if not torch.cuda.is_available(): sys.exit("no GPU visible; not running on CPU")
print("GPU:", torch.cuda.get_device_name(0))
PY
echo "[$(date)] ExPecto reference, $T"
"$CONDA_ENVS/expecto/bin/python" predict_reference_${T}.py \
  --tss_bed "$PROJECT_DIR/tss_hg38_mane_expressed_${T}.bed" \
  --out_dir "$PROJECT_DIR/predictions/expecto_reference_expanded_${T}"
echo "[$(date)] done"
