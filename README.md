# scRNASeq_Analysis

## Description
This is a repository of scripts used for single-cell RNA-sequencing analysis. We used 10k PBMC data obtained from the 10x Genomics website. 
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3

## R Analysis outline
The procedures performed are as follows:
1. Read 10X sequencing data and change it into a seurat object
2. Standard pre-processing workflow (QC and selecting cells for further analysis, Normalizing the data, Identification of highly variable features (feature selection, Scaling the data)
3. Perform linear dimensional reduction
4. Cluster the cells
5. Run non-linear dimensional reduction (UMAP/tSNE)

## Packages used in this analysis

R 4.2.1 (2022-06-23)
* [Seurat (4.0.6)](https://github.com/satijalab/seurat)
* [dplyr](https://dplyr.tidyverse.org/)
* [patchwork](https://github.com/thomasp85/patchwork)
* [ggplot2](https://ggplot2.tidyverse.org/)

The data was analyzed using Ubuntu 20.04.6 LTS

## Reference

Seurat - Guided Clustering Tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
