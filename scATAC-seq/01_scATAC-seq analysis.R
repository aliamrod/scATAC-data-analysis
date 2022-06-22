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
metadata <- "/work/Neuroinformatics_Core/s204365/ATACSeq_0001/metadata"
setwd(metadata)

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
# List fragment files:
# The CellRanger-ATAC count pipeline outputs a BED-like tabular file, where each line represents a unique ATAC-seq fragment captured by the assay. 
# Each fragment is created by two separate transposition events, which create the two ends of the observed fragment. Each unique fragment may generate multiple
# duplicate reads. These duplicate reads are collapsed into a single fragment record.

# $ zgrep -v "^#" fragments_104.tsv.gz | gzip -c > fragments_104.comments_removed.tsv.gz
# $ zgrep -v "^#" fragments_105.tsv.gz | gzip -c > fragments_105.comments_removed.tsv.gz
# $ zgrep -v "^#" fragments_106.tsv.gz | gzip -c > fragments_106.comments_removed.tsv.gz
# $ zgrep -v "^#" fragments_107.tsv.gz | gzip -c > fragments_107.comments_removed.tsv.gz

# Utilize function, reformatFragmentFiles() to reformat CellRanger-derived fragment files for reading in createArrowFiles(). It will handle anomalies found that generate
# errors in reading tabix bgzip'd fragment files.

# reformatFragmentFiles(): This function will help in reformatting Fragment Files for reading in createArrowFiles. It will handle odd anomalies found that 
# cause errors in reading tabix bgzip'd fragment files. 
reformatFragmentFiles("fragments_104_comments_removed.tsv.gz")
reformatFragmentFiles("fragments_105_comments_removed.tsv.gz")
reformatFragmentFiles("fragments_106_comments_removed.tsv.gz")
reformatFragmentFiles("fragments_107_comments_removed.tsv.gz")

# Inspect {reformatted} fragment files in current directory. 
reformat_files <- list.files(pattern = "Reformat")

# Create ArchR object.
inputFiles <- c(
  "fragments_104copy.comments_removed-Reformat.tsv.gz", "fragments_105_comments_removed-Reformat.tsv.gz", "fragments_106_comments_removed-Reformat.tsv.gz", "fragments_107_comments_removed-Reformat.tsv.gz")
names(inputFiles) <- c("fragments_104", "fragments_105", "fragments_106", "fragments_107")

## Setting default number of Parallel threads to 20. In Windows OS detection, parallel ArchRThread is set to 1. 
addArchRThreads(threads = 1) 

# Create ArchR object.
# Now, we generate Arrow Files which take ~10-15 minutes for each file. For each sample, this step will:
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
  subThreading = FALSE,
  addGeneScoreMat = TRUE
  )

ArrowFiles # Inspect ArrowFiles object (character vector of ArrowFile paths). 

# Change directory permissions to readable/writeable/executable permissions through Linux.
# $ chmod +rwx fragments_104.arrow
# $ chmod +rwx fragments_105.arrow
# $ chmod +rwx fragments_106.arrow
# $ chmod +rwx fragments_107.arrow
outputDirectory <- "/work/Neuroinformatics_Core/s204365/ATACSeq_0001/outputDirectory"
arrowfilesDirectory <- "/work/Neuroinformatics_Core/s204365/ATACSeq_0001/outputDirectory/ArrowFiles"

projCELL1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = outputDirectory, 
  copyArrows = TRUE #Recommended so that if you modify Arrow files you have an original copy for later usage + access TSS Enrichment Scores for each cell.
 )

# Access the cell names associated with each cell using '$' accessor for direct access to ```cellColData```.
head(projCELL1$cellNames)
quantile(projCELL1$TSSEnrichment)

# Add doublet score. 
doubScores <- addDoubletScores(
  input = ArrowFiles, 
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = TRUE
  )

# QC *****************************************************************************
proj_CELL_1 <- projCELL1
p <- ggPoint( # ggPoint: a ggplot-based dot plot wrapper function
  x = df[,1],
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight", 
  xlabel = "Log10UniqueFragments", 
  ylabel = "TSS Enrichment", 
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
  
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj_CELL1, addDOC = FALSE)

p1 <- plotGroups(
  ArchRProj = proj_CELL_1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment", 
  plotAs = "ridges"
  )
p2 <- plotGroups(
  ArchRProj = proj_CELL_1,
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment", 
  plotAs = "violin", 
  alpha = 0.4, 
  addBoxPlot = TRUE
  )
p3 <- plotGroups(
    ArchRProj = proj_CELL_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj_CELL_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_CAD_1, addDOC = FALSE, width = 4, height = 4)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



# Session Information.
Sys.Date()
sessionInfo()


