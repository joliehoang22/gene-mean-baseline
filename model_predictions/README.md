# model_predictions

## Gene panels and CAGE tracks

Run these before the predictions, after `gene_mean_baseline`.

| File | Description |
|------|-------------|
| `01_mane_tss.py` | MANE Select TSS for every autosomal protein-coding gene (Ensembl 109) |
| `02_expressed_panels.R` | MANE TSS restricted to expressed genes: 13,325 brain / 10,770 blood. The scoring scripts then keep genes with a FANTOM5 symbol match, giving the final 13,061 / 10,562 |
| `03_find_cage_tracks.py` | Lists the brain and whole-blood CAGE tracks in the model target set. This is how the tracks were chosen |
| `resources/` | The two gene panels and the track index lists, as used in the paper |

## Predictions

There is one folder per model, and all four models get the same input sequence. For each sample, `bcftools consensus` writes the phased SNVs into hg38 once per haplotype, and a window is cut around each gene's MANE TSS. Each haplotype is predicted separately and the two predictions are averaged. The models differ only in window length.

| Model | Window | Output used |
|---|---|---|
| Enformer | 196,608 bp | mean of the 10 central 128 bp bins (1,280 bp) at the TSS, 7 brain / 4 blood CAGE tracks |
| Basenji2 | 131,072 bp | same bins and tracks as Enformer |
| ExPecto | ±20 kb (200 Beluga windows of 2 kb; 43,000 bp of sequence) | one tissue head: `Brain_Frontal_Cortex_BA9` or `Whole_Blood` |
| Xpresso | 10.5 kb (−7 kb / +3.5 kb) | `humanMedian`, sequence only (half-life features held at the population mean) |

Scripts come in brain/blood pairs because the track sets and tissue heads differ. `<model>_lib_brain.py` and `<model>_lib_blood.py` hold the shared code for a model; `predict_personal_*.py` runs one sample; `predict_reference_*.py` runs the unmodified reference genome. Enformer has no lib file, because its personalized scripts contain the sequence code that the other models' libs copy. Its reference script is a single `predict_reference.py` that takes `--tissue brain` or `--tissue blood`.

## Environments

Each model folder has an `environment.yml`, exported from the conda environment used for the paper. The four cannot be combined, because Xpresso needs TensorFlow 1.15, Basenji2 needs TensorFlow 2.13 and ExPecto needs an older PyTorch than Enformer.

| Model | Python | Framework |
|---|---|---|
| Enformer | 3.10 | PyTorch 2.6.0 (CUDA 12.4), enformer-pytorch 0.8.10 |
| Basenji2 | 3.8 | TensorFlow 2.13.1 |
| ExPecto | 3.8 | PyTorch 1.13.1 (CUDA 11.7), xgboost 1.7.6 |
| Xpresso | 3.7 | TensorFlow 1.15.0, Keras 2.2.4 |

```bash
conda env create -f model_predictions/enformer/environment.yml -p $CONDA_ENVS/enformer
```

The versions also determine which GPUs each model can use. PyTorch wheels ship their own CUDA runtime, but TensorFlow wheels expect CUDA libraries on the system, so on a machine without them the two TensorFlow models run on CPU. PyTorch 1.13.1 has no kernels for GPUs newer than compute capability 8.6, so ExPecto needs a V100 or A100-class card.

## lsf/

| File | Description |
|------|-------------|
| `run_personal.sh` | Array job, one element per sample. Set `COHORT`, `MODEL` and `TISSUE` with `bsub -env` |
| `run_reference_<model>.sh` | Reference-genome arm, array of 2 (brain, blood) |
