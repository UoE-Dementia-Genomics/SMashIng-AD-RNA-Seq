# RNA-seq analysis workflow using Salmon and DESeq2

Example workflow using Salmon and DESeq2 to analyse short-read RNA-seq data looking at differential gene expression.

This example is based on human RNA-seq data

Workflow adapted from https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

Salmon documentation: https://salmon.readthedocs.io/en/latest/

DESeq2 documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

**Install salmon in a conda environment**

```
conda create -n salmon salmon
conda activate salmon
conda install -c bioconda salmon
```
