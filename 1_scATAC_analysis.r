# Load required packages onto BioHPC.

# Open terminal and enter the following command-line statements:
# $ module purge && module load shared slurm python/3.7.x-anaconda
# $ module load R/4.1.1-gccmkl
# $ module load hdf5_18/1.8.17
# $ module load gcc/8.3.0
# $ module load htslib
# $ module load gsl/2.4

# Open R/4.1.1-gccmkl module.
# $ R

# setwd(dir); specify a working directory. This may be from the root directory and transitioning up to current directory.
setwd('/work/Neuroinformatics_Core/s204365/ATACSeq_0001/metadata')

# scATAC data analysis pipeline (integrates ArchR, v1.0.1)//Dependencies
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("dplyr")

library(devtools)
library(BiocManager)
library(dplyr)


# install ArchR
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR) # version 1.0.2
ArchR::installExtraPackages()
set.seed(1)
addArchRGenome("mm10") # Mus musculus (house mouse) genome assembly

# Input data from 10X Cellranger-ATAC output
# list fragment files:
# The CellRanger-ATAC count pipeline outputs a BED-like tabular file, where each line represents a unique ATAC-seq fragment captured by the assay. 
# Each fragment is created by two separate transposition events, which create the two ends of the observed fragment. Each unique fragment may generate multiple
# duplicate reads. These duplicate reads are collapsed into a single fragment record.

# Inspect fragment files in current directory. 
fragment_tsv = list.files(pattern = ".tsv")

# Utilize function, reformatFragmentFiles() to reformat CellRanger-derived fragment files for reading in createArrowFiles(). It will handle anomalies found that generate
# errors in reading tabix bgzip'd fragment files.

# ****reformatFragmentFiles(
inputFiles <- c(
  "fragments_104.tsv", "fragments_105.tsv", "fragments_106.tsv", "fragments_107.tsv")
names(inputFiles) <- c("fragments_104", "fragments_105", "fragments_106", "fragments_107")

## Setting default number of Parallel threads to 20. In Windows OS detection, parallel ArchRThread is set to 1. 
addArchRThreads(threads = 1) 

# Create ArchR object.
# Now, we generate Arrow Files which take ~10-15 minutes. For each sample, this step will:
# 1. Read accessible fragments from the provided input files.
# 2. Calculate quality control (QC) information for each cell (i.e. TSS enrichment scores and nucleosome info). 
# 3. Filter cells based on quality control parameters.
# 4. Create a genome-wide TileMatrix using 500-bp bins.
# 5. Create a GeneScoreMatrix using the custom ```geneAnnotation``` that was defined when we called ```archArchRGenome()```. 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, # Do not set this too high due to notion that it can be adjusted to accordingly later.
  minFrags = 1000, 
  addTileMat = TRUE,
  force = FALSE,
  addGeneScoreMat = TRUE,
  cleanTmp = FALSE,
  nChunk = 1, 
  threads = 1, 
  subThreading = TRUE
  )

ArrowFiles # Inspect ArrowFiles object (character vector of ArrowFile paths). 





 
# Strict quality control (QC) of scATAC-seq data is essential to remove the contribution of low-quality cells. In ArchR, we will consider 3 characterstics of data:
# 1. The number of unique nuclear fragments (i.e. not mapping to mitochondrial DNA)
# 2. The signal-to-background ratio. Low signal-to-background ratio is often attributed to dead or dying cells which have de-chromatinzed DNA which allows for random transposition genome-wide.
# 3. The fragment size distribution. Due to nucleosomal periodicity, we expect to see depletion of fragments that are the length of DNA wrapped around a nucleosome (approximately 147 bp). 



  
