# scATAC data analysis pipeline (integrates ArchR, v1.0.1)
library(ArchR) # version
set.seed(10)
addArchGenome("hg38")

# Input data from 10X Cellranger-ATAC output
