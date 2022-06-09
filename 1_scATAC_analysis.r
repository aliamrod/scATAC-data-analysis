# Open terminal and enter the following command-line statements:
# module purge && module load shared slurmm python/3.7.x-anaconda
# module load R/4.1.1-gccmkl
# module load hdf5_18/1.8.17
# module load gcc/8.3.0

# scATAC data analysis pipeline (integrates ArchR, v1.0.1)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# install ArchR
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR) # version
ArchR::installExtraPackages()
set.seed(10)
addArchGenome("hg38")

# Input data from 10X Cellranger-ATAC output
"scA/fragments.tsv.gz", "scB/fragments.tsv.gz", "scC/fragments.tsv.gz", "scD/fragments.tsv.gz", "scE/fragments.tsv.gz", "scF/fragments.tsv.gz", "scG/fragments.tsv.gz", "scH/fragments.tsv.gz", "scI/fragments.tsv.gz", "scJ/fragments.tsv.gz", "scK/fragments.tsv.gz", "scL/fragments.tsv.gz", "scM/fragments.tsv.gz", "scN/fragments.tsv.gz", "scO/fragments.tsv.gz",  "scP/fragments.tsv.gz", "scQ/fragments.tsv.gz", "scR/fragments.tsv.gz", "scS/fragments.tsv.gz", "scT/fragments.tsv.gz", "scU/fragments.tsv.gz", "scV/fragments.tsv.gz", "scW/fragments.tsv.gz", "scX/fragments.tsv.gz", "scY/fragments.tsv.gz", "scZ/fragments.tsv.gz", "scAA/fragments.tsv.gz", "scAB/fragments.tsv.gz", "scAC/fragments.tsv.gz", "scAD/fragments.tsv.gz", "scAE/fragments.tsv.gz", "scAF/fragments.tsv.gz", "scG/fragments.tsv.gz", "scAH/fragments.tsv.gz", "scAI/fragments.tsv.gz","scAJ/fragments.tsv.gz","scAK/fragments.tsv.gz","scAL/fragments.tsv.gz","scAM/fragments.tsv.gz","scAN/fragments.tsv.gz","scAO/fragments.tsv.gz","scAP/fragments.tsv.gz","scAQ/fragments.tsv.gz","scAR/fragments.tsv.gz")
names(inputFiles) <- c("scA", "scB", "scC", "scD", "scE", "scF", "scG", "scH", "scI", "scJ", "scK", "scL", "scM", "scN", "scO", "scP", "scQ", "scR", "scS", "scT", "scU", "scV", "scW", "scX", "scY", "scZ", "scAA", "scAB", "scAC", "scAD", "scAE", "scAF", "scAG", "scAH", "scAI", "scAJ", "scAK", "scAL", "scAM", "scAN", "scAO", "scAP", "scAQ", "scAR")

# create ArchR object
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, 
  filterFrags = 1000,
  addTileMat = TRUE, 
  force = FALSE, 
  addGeneScoreMat = TRUE
  )
# Strict quality control (QC) of scATAC-seq data is essential to remove the contribution of low-quality cells. In ArchR, we will consider 3 characterstics of data:
# 1. The number of unique nuclear fragments (i.e. not mapping to mitochondrial DNA)
# 2. The signal-to-background ratio. Low signal-to-background ratio is often attributed to dead or dying cells which have de-chromatinzed DNA which allows for random transposition genome-wide.
# 3. The fragment size distribution. Due to nucleosomal periodicity, we expect to see depletion of fragments that are the length of DNA wrapped around a nucleosome (approximately 147 bp). 

projCAD1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "CAD", 
  copyArrows = TRUE #Recommended so that if you modify Arrow files, you'll have an original copy for later usage.
  )

# add doublet score
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #refers to how many cells near a "pseudo-doublet" to count
  knnMethod = "UMAP", # refers to the embedding to use for nearest neighbor search
  LSIMethod = 1, 
  force = TRUe
  )
# basic QC
proj_CAD_1 <- projCAD1
p <- ggpoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight", 
  xlabel = "Log10 Unique Fragments", 
  ylabel = "TSS Enrichment", 
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vlin(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProf = proj_CAD_1, addDOC = FALSE)
  
  
  //////////////////////////


  
