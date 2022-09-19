# ETHLR: A Tensor Method based on Enhanced Tensor Nuclear Norm and Hypergraph Laplacian Regularization for Pan-Cancer Omics Data Analysis


**A novel method named ETHLR is designed to simultaneously consider complementary information among multiple views and unique information for each view.**


![](https://github.com/nayu0419/ETHLR/blob/main/ETHLR_Pipeline.png)

## requirements
- MATLAB R2020a

## run

```
ETHLR_main.m
```

## File:

- funs: subfunction used in the ETHLR process
- Datasets: the pan-cancer omics data (for example, PAAD_CHOL_ESCA)
- The original dataset can be downloaded in [TCGA] (https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)
- Put the other datasets in the folder "Datasets"

where *dataset* contains three categories of information, namely gene expression (GE), copy number variation (CNV) and methylation (ME). We represent a pan-cancer omics dataset as a third-order tensor.
