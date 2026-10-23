# Paths used by every script in this repo. Edit once, then run
#   source config/paths.sh
# from the repo root before running anything. The LSF drivers source this
# file themselves, so submit them from the repo root too.
#
# PROJECT_DIR is the working directory: predictions, expression matrices and
# score tables are all written under it. Files the scripts expect to find there:
#   ref_genome/hg38.fa (+ .fai)
#   phased_vcfs/                       LBP LIV, chr1..chr22.phased.vcf.gz
#   targets_human.txt                  Basenji2/Enformer target list (Basenji repository)
#   brain_tracks.txt, blood_tracks.txt copied from model_predictions/resources/
#   sample_mapping.tsv                 LBP LIV brain: VCF sample ID <tab> RNA-seq sample ID, with header
#   samples_blood.tsv                  LBP LIV blood, no header
#   sample_mapping_rosmap.tsv          ROSMAP, with header
#   samples_covid.tsv                  COVID-19, no header
# Everything else (Null/, predictions/, evaluation/, figures/) is written by the scripts.

export PROJECT_DIR=/path/to/working_dir
export PLOT_DIR=$PROJECT_DIR/figures

## Reference genome and annotation
export HG38_FASTA=$PROJECT_DIR/ref_genome/hg38.fa      # needs a .fai next to it
export GENCODE_V30_GTF=/path/to/gencode.v30.primary_assembly.annotation.gtf
export ENSEMBL_109_GTF=/path/to/Homo_sapiens.GRCh38.109.gtf

## Phased VCFs: one chr{1..22}.phased.vcf.gz per directory, biallelic SNVs only
export VCF_DIR_LBP=$PROJECT_DIR/phased_vcfs
export VCF_DIR_ROSMAP=/path/to/rosmap/phased_vcfs
export VCF_DIR_COVID=/path/to/covid/phased_vcfs

## Expression data
export GTEX_TPM=/path/to/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz
export GTEX_SAMPLE_ATTR=/path/to/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt
export ROSMAP_COUNTS=/path/to/rosmap_featureCounts.txt
export MSBB_COUNTS=/path/to/msbb_featureCounts.txt
export LBP_BRAIN_RDS=/path/to/lbp_brain_featureCounts.RDS
export LBP_BLOOD_RDS=/path/to/lbp_blood_featureCounts.RDS
export COVID_COUNTS=/path/to/rnaseq_raw_count_matrix.csv
export FANTOM5_CAGE=/path/to/hg38_fair+new_genes_phase1and2_tpm.osc.txt.gz
export FANTOM5_MEANS=$PROJECT_DIR/fantom5_cage_gene_means.csv   # written by gene_mean_baseline/03

## Models
export CONDA_ENVS=/path/to/conda/envs       # contains enformer/, basenji2/, expecto/, xpresso/
export PEB_REPO=/path/to/personalized-expression-benchmark   # Basenji2 code, Beluga class, Xpresso weights
export EXPECTO_REPO=/path/to/ExPecto        # with resources/ and models/ downloaded
export ENFORMER_WEIGHTS=EleutherAI/enformer-official-rough    # or a local snapshot directory
export HF_HOME=$PROJECT_DIR/.cache/huggingface
export TRANSFORMERS_CACHE=$HF_HOME
