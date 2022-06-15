# snATAC-data-analysis

The following repository contains data analysis scripts influenced by "Single-nucleus chromatin accessibility profiling highlights regulatory mechanisms of coronry artery disease risk". 
In the following repository, 1_scATAC_analysis.r includes scripts for the single nucleus ATAC-seq data analysis, 2_sc_RNA_analysis.r includes scripts for the single-cell RNA-seq data analysis, and 3_SVMpipeline/ and 4_RASQUAL_caQTL/ folders include various scripts for the variant effect predictions and allele-specific based caQTL mapping analysis using snATAC-seq peaks, respectively.



As a review, 
* scATAC-seq => Single-cell sequencing assay for transposase-accessible chromatin (scATAC-seq) can be used to identify cell subpopulations with different chromatin accessibility profiles within complex samples, eliminating the need for isolation strategies like FACS or magnetic sorting that could alter the biology due to sample manipulation.
*For instance:
      * Identifying cancer stem cells or infiltrating macrophages within a tumor sample.
      * Identifying novel cell subpopulations that are responsible for response to drug treatments (i.e. responders vs. resistant cells).
      * Identifying subpopulations of cells with variations in chromatin accessibility that can provide insight into developmental trajectories (i.e. brain development,          T-helper development, B-cell differentiation). 
* scRNA-seq => Single-cell RNA sequencing (scRNA-seq) can reveal complex and rare cell populations, uncover regulatory relationships between genes, and track the trajectories of distinct cell lineages in development. In contrast to bulk RNA-seq which relies on averaged gene expression from a population of cells to reveal RNA presence as well as quantity in a sample of cells during time of measurement, scRNA-seq offers an unbiased apporach to investigate the cellular heterogeneity and dynamics of a biological system.
* variant effect => Gene variant conforms to a permanent change in the DNA sequence that comprises a gene. 
* caQTL mapping => chromatin accessibility quatitative trait loci (caQTLs) help identify GWAS loci that may alter GWAS traits by modulating chromatin structure, but caQTL have been identified in a limited set of tissues [2]. (1) caQTLs are highly cell-type specific, (2) caQTLs have higher effect sizes on average than eQTLs (expression quantitative trait loci), (3) caQTLs can be used to identify causal variants, cell types, and mechanisms of non-coding GWAS loci. 



## 0. Requirements. 
The following environments/package are required as dependencies to run the scripts.
Single-cell data analysis (scATAC_analysis.r and scRNA_analysis.r)
- R v4.1.1-gccmkl
- ArchR v1.0.1
- Seurat v4.0.0
- dplyr v1.0.4
- patchwork v1.1.1
- colorRamps v2.3

Variant Scoring (SVMpipeline)
- Python v3.7
- perl v5.16
- numpy v1.1
- scipy v1.4
- IPython v7.21
- matplotlib v3.2
- plotnine v0.7
- Two annotation data hg38.chrom.sizes

## 1. Single Nucleus ATAC-seq (scATAC) Data Analysis
-  basicQC
-  cell filtering
-  dimensional reduction
-  cell clustering
-  UMAP visualization
-  batch effect detection
-  scRNA-seq integration
-  cell type assignment
-  peak calling 
-  detection of cell-type-specific genes/peaks
-  motif annotation
-  chromVAR
-  trajectory analysis
-  co-accessibility
-  footprint analysis (to estimate transcription factor binding; this continuous footprint score is correlated with the presence of transcription factor binding sites in the genome)

## 2. Single-cell RNA-seq (scRNA) Data Analysis.
All the scripts for scRNA data analysis were included in 2_scRNA_analysis.r.

The scRNA data generated in a previous study in the same system as snATAC (Wirka et al., _Nat Med._ 2019, GSE131778). The primary objective for the scRNA data analysis is to successfully generate an scRNA reference database for the integration analysis and label transferring analysis of the snATAC data. The script contains the following parts.
-  cell and gene filtering
-  dimensional reduction
-  cell clustering
-  UMAP visualization
-  comparing the marker genes (Wirka et al., _Nat Med_ 2019, GSE131778) versus cell clustering
-  cell-type assignment based on marker genes
-  check knowledge-based marker genes versus cell-type assignment

## 3. Variant effect predictions using snATAC-seq peaks.



## 5. Bulk data process (bulk ATAC-seq, superenhancer analysis with H3K27ac ChIP-seq).
The cmd lines for processing bulk data can be found in 5_bulkDataAnalysis/bulkdata_process.sh. All the related scripts and configuration files can be found in the same folder.

### References.
[1]
[2] Currin KW, Erdos MR, Narisu N, Rai V et al. "Genetic effects on liver chromatin accessibility identify disease regulatory variants". _Am J Hum Genet_. 2021 Jul. Date of Access: 07 June, 2022.
