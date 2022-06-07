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
  ylim = 
  
  
  //////////////////////////


    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj_CAD_1, addDOC = FALSE)

p1 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj_CAD_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_CAD_1, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj_CAD_1)
p2 <- plotTSSEnrichment(ArchRProj = proj_CAD_1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj_CAD_1, addDOC = FALSE, width = 5, height = 5)

p <- ggPoint(
    x = df2[,"log10(nFrags)"], 
    y = df2[,"TSSEnrichment"], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(3, quantile(df2[,"log10(nFrags)"], probs = 0.99)),
    ylim = c(4, quantile(df2[,"TSSEnrichment"], probs = 0.99))
) + geom_hline(yintercept = 7, lty = "dashed") + geom_vline(xintercept = 4, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags_cutoff.pdf", ArchRProj = proj_CAD_2, addDOC = FALSE)


# filter cells
idxPass <- which(proj_CAD_2$TSSEnrichment >= 7 & proj_CAD_2$nFrags >= 10000)
proj_CAD_2 <- filterDoublets(proj_CAD_1,filterRatio=1.5)
df2 <- getCellColData(proj_CAD_2,select = c("log10(nFrags)", "TSSEnrichment"))
cellsPass <- proj_CAD_2$cellNames[idxPass]
proj_CAD_2 <- proj_CAD_2[cellsPass, ]

# dimensional reduction
proj_CAD_2 <- addIterativeLSI(
    ArchRProj = proj_CAD_2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    seed=1,force=T
)
