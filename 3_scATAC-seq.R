## Load packages
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(tidyverse)
library(patchwork)

# Load original data. 
counts <- Read10X_h5(filename = "/work/Neuroinformatics_Core/s204365/ATACSeq_0001/metadata/filteredpeak.h5")
metadata <- read.csv(
file = "singlecell.csv"**, 
header = TRUE,
row.names = 1
)

chrom_assay <- CreateChromatinAssay(
counts = counts,
sep = c(":" , "-"),
genome = 'mm10', 
fragments = ***
min.cells = 10, 
min.features = 200
)

pbmc <- CreateSeuratObject(
counts = chrom_assay
