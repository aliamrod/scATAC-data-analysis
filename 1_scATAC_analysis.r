# Load required packages onto BioHPC.

# Open terminal and enter the following command-line statements:
# $ module purge && module load shared slurmm python/3.7.x-anaconda
# $ module load R/4.1.1-gccmkl
# $ module load hdf5_18/1.8.17
# $ module load gcc/8.3.0
# $ module load gsl/2.4
# $ R

# setwd(dir); specify a working directory. This may be from the root directory and transitioning up to current directory.
setwd('/work/Neuroinformatics_Core/s204365/ATACSeq_0001')

# scATAC data analysis pipeline (integrates ArchR, v1.0.1)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")


# install ArchR
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR) # version
ArchR::installExtraPackages()
set.seed(1)
addArchGenome("mm10") # Mus musculus (house mouse) genome assembly

# Input data from 10X Cellranger-ATAC output
# list fragment files:
fragment_tsv = list.files(pattern = ".tsv")
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
# > ArrowFiles

projNEU1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/work/Neuroinformatics_Core/s204365/ATACSeq_0001/output_dir", 
  copyArrows = TRUE # Boolean value indicating whether ArrowFiles should be copied into ```outputDirectory```. This is recommended so that if you have to modify the ArrowFiles, you have an original copy for later usage. 
  )

  



 
# Strict quality control (QC) of scATAC-seq data is essential to remove the contribution of low-quality cells. In ArchR, we will consider 3 characterstics of data:
# 1. The number of unique nuclear fragments (i.e. not mapping to mitochondrial DNA)
# 2. The signal-to-background ratio. Low signal-to-background ratio is often attributed to dead or dying cells which have de-chromatinzed DNA which allows for random transposition genome-wide.
# 3. The fragment size distribution. Due to nucleosomal periodicity, we expect to see depletion of fragments that are the length of DNA wrapped around a nucleosome (approximately 147 bp). 



  
