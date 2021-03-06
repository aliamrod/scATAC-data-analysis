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

# To plot.
library(ggplot2)
library(patchwork)
library(stringr)
library(repr)
set.seed(1)

# Read ATAC data.
counts <- Read10X_h5(filename = "filtered_peak_bc_matrix_107.h5")

# Read metadata.
metadata <- read.csv(file = "singlecell_107.csv", 
                     header = TRUE, 
                     row.names = 1)
# Create a chromatin assay utilizing the count matrix.
brain_assay <- CreateChromatinAssay(counts = counts,
                                    sep = c(":", "-"), 
                                    genome = "mm10", 
                                    fragments = "fragments_107.tsv.gz", 
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

#####ERROR#########
# plotPDF(g1,g2, name = "Nucleosome-Signal.pdf", addDOC = FALSE, width = 5, height = 5)

# compute nucleosome signal score per cell
brain <- NucleosomeSignal(object = brain)

# Compute TSS enrichment score per cell.
brain <- TSSEnrichment(object = brain, fast = FALSE)

# Add blacklist ratio and fraction of reads in peaks
brain$pct_reads_in_peaks <- brain$peak_region_fragments / brain$passed_filters * 100
brain$blacklist_ratio <- brain$blacklist_region_fragments / brain$peak_region_fragments
######



counts <- Read10X_h5(filename = "../vignette_data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../vignette_data/atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '../vignette_data/atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


#NucleosomeSignal
brain <- NucleosomeSignal(object = brain)

#Cells which are outliers for the mononucleosomal/ nucleosome-free ratio have different banding patterns. 
#The remaining cells exhibit a pattern that is typical for a successful ATAC-seq experiment.
options(repr.plot.width=15, repr.plot.height=4)
boxplot(brain$nucleosome_signal, main = "Nucleosome Signal", horizontal = TRUE)
axis(1, at=seq(0,40,5)) #change numeration in axis
abline(v=4, col = "Red", lty = 5, lwd = 3) #add line at 4 


##BUG****
options(repr.plot.width=10, repr.plot.height=5)
brain$nucleosome_group <- ifelse(brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = brain, group.by = 'nucleosome_group', region = 'chr1-1-10000000')



# Visualize peak_count_per_cell
options(repr.plot.width=20, repr.plot.height=8)
par(mfrow=c(2,4))
hist(brain@meta.data[["nucleosome_signal"]], breaks = 50) 
# hist(brain@meta.data[["TSS.enrichment"]], breaks = 50)
hist(brain@meta.data[["peak_region_fragments"]], breaks = 50)
hist(brain@meta.data[["pct_reads_in_peaks"]], breaks = 50)
hist(brain@meta.data[["blacklist_ratio"]], breaks = 50)
hist(brain@meta.data[["nCount_ATAC"]], breaks = 50)
hist(brain@meta.data[["nFeature_ATAC"]], breaks = 50)
