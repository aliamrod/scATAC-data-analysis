# module purge && module load shared slurm python/3.7.x-anaconda
# module load htslib
# module load gsl/2.4
# module load gcc/8.3.0 

# module load R/4.1.1-gccmkl

library(Signac)
library(Seurat)
library(SeuratDisk)
library(BiocManager)
BiocManager::install("biovizBase")
library(biovizBase)

# Genome
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)

# To plot
library(ggplot2)
library(patchwork)
library(stringr)
library(repr)
set.seed(1)

# Read ATAC data.
counts <- Read10X_h5("filtered_peak_bc_matrix_104.h5")

# Read metadatal
metadata <- read.csv(file = "singlecell_104.csv", 
                     header = TRUE, 
                     row.names = 1)
# Create a chromatin assay utilizing the count matrix.
brain_assay <- CreateChromatinAssay(counts = counts,
                                    sep = c(":", "-"), 
                                    genome = "mm10", 
                                    fragments = "fragments_104.tsv.gz", 
                                    min.cells = 1)
# The ChromatinAssay class extends the standard Seurat Assay class and adds several additional slots for data useful for the analysis of single-cell chromatin datasets.
# Create Seurat Object.
brain <- CreateSeuratObject(counts = brain_assay, 
                            assay = 'ATAC', 
                            project = 'ATAC', 
                            meta.data = metadata)
                                  
# Extract gene annotations from EnsDb. 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# Add the gene information to the object.
Annotation(brain) <- annotations


plotPDF(g1,g2, name = "Nucleosome-Signal.pdf", addDOC = FALSE, width = 5, height = 5)



