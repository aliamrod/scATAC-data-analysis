# snATAC-data-analysis

The following repository contains data analysis scripts influenced by "Single-nucleus chromatin accessibility profiling highlights regulatory mechanisms of coronry artery disease risk". 
In the following repository, 1_scATAC_analysis.r includes scripts for the single nucleus ATAC-seq data analysis, 2_sc_RNA_analysis.r includes scripts for the single-cell RNA-seq data analysis, and 3_SVMpipeline/ and 4_RASQUAL_caQTL/ folders include various scripts for the variant effect predictions and allele-specific based caQTL mapping analysis using snATAC-seq peaks, respectively.



As a review, 
* scATAC-seq => Single-cell sequencing assay for transposase-accessible chromatin (scATAC-seq) can be used to identify cell subpopulations with different chromatin accessibility profiles within complex samples, eliminating the need for isolation strategies like FACS or magnetic sorting that could alter the biology due to sample manipulation. For instance:
      * Identifying cancer stem cells or infiltrating macrophages within a tumor sample.
      * Identifying novel cell subpopulations that are responsible for response to drug treatments (i.e. responders vs. resistant cells).
      * Identifying subpopulations of cells with variations in chromatin accessibility that can provide insight into developmental trajectories (i.e. brain development,          T-helper development, B-cell differentiation). 
* scRNA-seq =>
* variant effect =>
* caQTL mapping => chromatin accessibility quatitative trait loci (caQTLs) help identify GWAS loci that may alter GWAS traits by modulating chromatin structure, but caQTL have been identified in a limited set of tissues [2]. (1) caQTLs are highly cell-type specific, (2) caQTLs have higher effect sizes on average than eQTLs (expression quantitative trait loci), (3) caQTLs can be used to identify causal variants, cell types, and mechanisms of non-coding GWAS loci. 



References.
[1]
[2] Currin KW, Erdos MR, Narisu N, Rai V et al. "Genetic effects on liver chromatin accessibility identify disease regulatory variants". _Am J Hum Genet_. 2021 Jul. Date of Access: 07 June, 2022.
