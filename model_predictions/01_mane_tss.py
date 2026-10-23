#!/usr/bin/env python
"""
Writes a TSS bed for every autosomal MANE Select gene, in the 5-column format
used by every bed in this repo:  chrom  tss  strand  gene_id  gene_name
Input is an Ensembl-format GTF (unprefixed chromosomes, unversioned gene_id).
"""
import gzip, os, sys
import pandas as pd

GTF = sys.argv[1] if len(sys.argv) > 1 else os.environ["ENSEMBL_109_GTF"]
OUT = os.path.join(os.environ["PROJECT_DIR"], "tss_hg38_mane_all_autosomal.bed")

AUTO = {str(i) for i in range(1, 23)}
opener = gzip.open if GTF.endswith(".gz") else open
rows = []
with opener(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "transcript" or f[0] not in AUTO or 'tag "MANE_Select"' not in f[8]:
            continue
        attr = dict(kv.strip().split(" ", 1) for kv in f[8].strip().rstrip(";").split(";") if " " in kv.strip())
        gid = attr.get("gene_id", "").strip('"').split(".")[0]
        gnm = attr.get("gene_name", "").strip('"')
        if not gid:
            continue
        tss = int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1   # 0-based, strand-aware
        rows.append(("chr" + f[0], tss, f[6], gid, gnm))

df = pd.DataFrame(rows, columns=["chrom","tss","strand","gene_id","gene_name"]).drop_duplicates(subset="gene_id")
df.to_csv(OUT, sep="\t", header=False, index=False)
print(f"{len(df)} autosomal MANE Select genes -> {OUT}")
print(df.head(3).to_string(index=False))
