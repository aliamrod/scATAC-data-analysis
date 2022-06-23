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
library(Seurat)


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
  "fragments_104_comments_removed-Reformat.tsv.gz", "fragments_105_comments_removed-Reformat.tsv.gz", "fragments_106_comments_removed-Reformat.tsv.gz", "fragments_107_comments_removed-Reformat.tsv.gz")
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

# Add doublet score. 
doubScores <- addDoubletScores(
  input = ArrowFiles, 
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = TRUE
  )

# Access the cell names associated with each cell using '$' accessor for direct access to ```cellColData```.
head(projCELL1$cellNames) 
quantile(projCELL1$TSSEnrichment) # accesses the TSS Enrichment Scores for each cell. 

# Subsetting ArchRProject by cells. We can subset the project numerically, for instance obtaining the first 200 cells in the project:
idxPass <- which(projCELL1$TSSEnrichment >=8)
cellsPass <- projCELL1$cellNames[idxPass]
projCELL1[cellsPass, ]

bioNames <- gsub("_R2|R1|scATAC_", "", projCELL1$Sample)
head(bioNames)
projCELL1$bioNames <- bioNames
bioNames <- bioNames[1:20]
cellNames <- projCELL1$cellNames[1:20]
projCELL1 <- addCellColData(ArchRProj = projCELL1, data = paste0(bioNames), 
                            cells = cellNames, name = "bioNames2")
getCellColData(projCELL1, select = c("bioNames", "bioNames2"))
# ArchR provides the getCellColData() function to enable easy retrieval of metadata columns from an ArchRProject. For instance, we can retrieve a column name, 
# such as the number of unique nuclear (i.e. non-mitochondrial) fragments per cell:
df <- getCellColData(projCELL1, select = "nFrags")
# Instead of selecting by a column by name, we can perform operations on a given column using its column name:
df <- getCellColData(projCELL1, select = c("log10(nFrags)", "nFrags - 1"))

# Plotting QC metric - log10(Unique Fragments) vs. TSS Enrichment Score. 
# Repeating the examples above, we can easily obtain standard scATAC-seq metrics for quality control of individual cells. Here, we use TSS Enrichment Score, a
# measure of signal-to-background in ATAC-seq data as well as the number of unique nuclear fragments (because cells with very few fragments do not have enough
# data to confidently analyze). 
df <- getCellColData(projCELL1, select = c("log10(nFrags)", "TSSEnrichment"))

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

saveArchRProject(ArchRProj = projCELL1, outputDirectory = outputDirectory, load = FALSE)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_CELL_1, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj_CELL_1)
p2 <- plotFragmentSizes(ArchRProj = proj_CELL_1)
plotPDF(p1, p2, name = "QC-Sample-FragmentSizes-TSSProfile.pdf", ArchRProj = proj_CELL_1, addDOC = FALSE, width = 5, height = 5)

# Filter cells.
idxPass <- which(proj_CELL_2$TSSEnrichment >= 7 & proj_CELL_2$nFrags >= 10000)
proj_CELL_2 <- filterDoublets(proj_CELL_1, filterRatio = 1.5)







# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    



# Session Information.
Sys.Date()
sessionInfo()


