BPH scRNAseq scripts used in Nature Communications Paper
============================================================

scripts used in the analysis of all cells and leukocytes associated prostates from men with BPH

Rmarkdown files are included that can be used to process any CellRanger generated 10x dataset (data should include files in the /outs/filtered_feature_bc_matrix resulting from running 10x Genomics's CellRanger software).

Operating systems and dependencies used in the analysis
--------------------------------------------------------
R version 3.5.1 (2018-07-02)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /apps/halstead7/R/prerequisites/halstead/openblas/0.2.19_gcc-6.3.0/lib/libopenblas_haswellp-r0.2.19.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ReactomePA_1.26.0      clusterProfiler_3.10.1 org.Hs.eg.db_3.7.0     AnnotationDbi_1.44.0   IRanges_2.16.0         S4Vectors_0.20.1       Biobase_2.42.0         BiocGenerics_0.28.0    vioplot_0.3.6          zoo_1.8-6              sm_2.2-5.6             lattice_0.20-38        biomaRt_2.38.0         qusage_2.16.1          limma_3.38.3           scales_1.1.1           RColorBrewer_1.1-2     future.apply_1.7.0     future_1.21.0          sctransform_0.3.2      wesanderson_0.3.6      clustree_0.4.3         ggraph_2.0.3           ggrepel_0.9.1          ComplexHeatmap_1.20.0  ggforce_0.3.3          forcats_0.4.0          stringr_1.4.0          dplyr_1.0.5            purrr_0.3.4            readr_1.3.1            tidyr_1.0.0            tibble_2.1.3           ggplot2_3.3.3         
[35] tidyverse_1.2.1        Seurat_3.1.3           rmdformats_1.0.1       knitr_1.25    

Installation
--------------------------------------------------------
Simply install R through CRAN and then install packages listed above through CRAN or Bioconductor.  Typical install time is ~20 minutes de novo.  

To install Bioconductor packages simply do:

```if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("packageName")
```
To install R packages through CRAN:

```install.packages("packageName")```

To run scripts
--------------------------------------------------------
   *After installation of CRAN and Bioconductor packages, change directory names to reflect your current directory, change parameters as desired and output names of files as desired.  Typical runtime is ~20 minutes and results in an html formatted document detailing main quality control, clustering, and marker gene results.
   
*Scripts included are not novel software scripts, but are RMarkdown (*.Rmd) scripts or R (*.R) scripts that utilize open-source, well-documented software packages.

*RMarkdown scripts should be run in numerical order:

 **001-Seurat-From-Cellranger.Rmd**  
 
   This script reads in barcodes, features, and count matrices generated from 10X CellRanger Count for individual samples, performs additional quality control, scaling, and normalization, dimensionality reduction, and clustering.

 **002-Combined_analysis.Rmd**

  This script combines previously processed samples using Seurat's anchor-based clustering method, followed by dimensionality reduction, clustering, and marker gene detection.
 
 **003-DEAnalysis-edgeR.Rmd**
 
 This script reads in barcodes, features, and count matrices generated from 10X CellRanger Count for individual samples, performs additional quality control, scaling, and normalization, dimensionality reduction, and clustering.
 
 **004-leukocyte_plots.R**
 
 This R script generates plots specific to the Nature Communication publication on leukocyte clustering generated by the R package Seurat.
 
 **005-allCellsCombinedClustering.Rmd**
 
 This script reads in barcodes, features, and count matrices generated from 10X CellRanger Count for samples that contain all cells associated with "large" BPH prostates (rather than only leukocytes), performs additional quality control, scaling, and normalization, dimensionality reduction, clustering, and marker gene identification.
 
*Ensure directory structure is appropriate for your system and click "knit" for RMarkdown scripts (scripts 001-003 and 005)
*For script 004, run scripts as desired piecemeal in RStudio or R, or on commandline via the command ```Rscript 004-leukocyte_plots.R```

*If these scripts are used on another dataset, other parameters that may need to be changed include quality control parameters (filtering cutoffs), dimensionality reductions parameters (such as number of variable genes to include and number of principal components to include), and clustering parameters that Seurat uses in implementing the Louvain algorithm (resolution or res)
