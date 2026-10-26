# Arithmetic means outperform AI at predicting gene expression

This repository contains the code for Hoang *et al.*, "Arithmetic means outperform AI at predicting gene expression".

We compare four genomic deep learning models that predict gene expression from DNA sequence (Enformer, Basenji2, ExPecto and Xpresso) with a gene-mean baseline, which is the average expression of each gene in a separate dataset from the same tissue. Performance is measured with the cross-gene correlation, the correlation between the predicted expression of a set of genes and the expression measured in the same person's tissue.

## Main analyses

1. Each model is run on every individual's phased genome (both haplotypes, predictions averaged) and scored against that individual's brain or blood RNA-seq. Gene-mean baselines from the reference cohorts listed below are scored against the same samples.
2. Each model is also run on the hg38 reference sequence, for the same genes and samples, to compare personalized and reference-genome predictions.
3. To measure what a model adds to a gene-mean baseline, we fit two linear models per sample across genes, `observed ~ gene_mean` (null) and `observed ~ gene_mean + model` (full), and report ΔR² = R²full − R²null.


## Cohorts

### Study cohorts

Models and gene-mean baselines are scored against the measured expression of the people in these cohorts. Each cohort has genome sequence and RNA-seq from the same individuals.

| Cohort | Tissue | Role | Samples | Participants |
|---|---|---|---|---|
| Living Brain Project living donors (LBP LIV) | prefrontal cortex | brain discovery | 206 | 142 |
| LBP LIV | whole blood | blood discovery | 214 | 144 |
| ROSMAP | dorsolateral prefrontal cortex | brain replication | 579 | 579 |
| Mount Sinai COVID-19 Biobank (COVID-19) | whole blood | blood replication | 474 | 474 |

### Reference cohorts

Gene-mean baselines are computed from these datasets. In Fig. 1, both study cohorts of a tissue are compared with the same baselines. The ΔR² analysis (Fig. 2, Supp. Fig. 5) adds a baseline from the other cohort of the same tissue.

| Dataset | Tissue | Samples | Reference for |
|---|---|---|---|
| LBP postmortem donors (LBP PM) | prefrontal cortex | 243 | LBP LIV brain, ROSMAP |
| GTEx v10 | frontal cortex (BA9) | 269 | LBP LIV brain, ROSMAP |
| Mount Sinai Brain Bank (MSBB) | frontal pole (BA10) | 325 | LBP LIV brain, ROSMAP |
| FANTOM5 | brain | 7 libraries | LBP LIV brain, ROSMAP |
| ROSMAP | dorsolateral prefrontal cortex | 639 | LBP LIV brain (ΔR² only) |
| LBP LIV | prefrontal cortex | 275 | Supp. Fig. 4 only |
| GTEx v10 | whole blood | 803 | LBP LIV blood, COVID-19 |
| FANTOM5 | blood | 4 libraries | LBP LIV blood, COVID-19 |
| COVID-19 | whole blood | 1,392 | LBP LIV blood (ΔR² only) |
| LBP LIV | whole blood | 243 | COVID-19 (ΔR² only) |

We evaluated the protein-coding genes expressed in each tissue that are present in every gene-mean source, FANTOM5 included: 13,061 genes in brain and 10,562 in blood.

---

## Repository structure

```
├── config/                   file paths (edit paths.sh)
├── gene_mean_baseline/       observed expression and gene-mean baselines
├── model_predictions/        gene panels, CAGE tracks, model scripts and environments
├── cross_gene_correlation/   cross-gene correlations for every arm, QC and tests
├── delta_r2/                 what each model adds to a gene-mean baseline
└── figures/                  one script per figure
```

## Setup

Each model runs in its own conda environment, exported from the ones used for the paper (`environment.yml` in each model's folder in `model_predictions/`). Model code and pretrained weights come from their original releases:

| Model | Source |
|---|---|
| Enformer | [`enformer-pytorch`](https://github.com/lucidrains/enformer-pytorch), weights `EleutherAI/enformer-official-rough` on Hugging Face (we used snapshot `affe5713ae9017460706a44108289b13c5fee16c`) |
| Basenji2 | original code at [calico/basenji](https://github.com/calico/basenji); we used the version in [ni-lab/personalized-expression-benchmark](https://github.com/ni-lab/personalized-expression-benchmark) (`basenji2/`), with weights `model_human.h5` placed in its `basenji2/resources/` |
| ExPecto | Beluga class from the same ni-lab repository; weights and per-tissue models from [FunctionLab/ExPecto](https://github.com/FunctionLab/ExPecto) (`download_resources.sh`) |
| Xpresso | original code at [vagarwal87/Xpresso](https://github.com/vagarwal87/Xpresso); we used weights `humanMedian_trainepoch.11-0.426.h5` from the ni-lab repository (`xpresso/xpresso/pretrained_models/`) |

All file paths are set in `config/paths.sh`. Edit it, then run `source config/paths.sh` from the repo root before running any script.

---

## Computational environment

Enformer and ExPecto need a GPU. Basenji2 and Xpresso can run on a GPU or a CPU. All other steps (scoring, ΔR² and figures) run on a CPU in R 4.2.0. We ran the predictions as LSF jobs on [Minerva](https://labs.icahn.mssm.edu/minervalab/), the high-performance computing cluster at the Icahn School of Medicine at Mount Sinai.

---

## Acknowledgements

This work was supported in part through the computational and data resources and staff expertise provided by Scientific Computing and Data at the Icahn School of Medicine at Mount Sinai and supported by the Clinical and Translational Science Awards (CTSA) grant UL1TR004419 from the National Center for Advancing Translational Sciences.
