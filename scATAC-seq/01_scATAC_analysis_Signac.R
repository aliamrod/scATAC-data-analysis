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


