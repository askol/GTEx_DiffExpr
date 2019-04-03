library(annotables)
library(ggplot2)
PlotDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/Plots/"

## READ IN DATA ##
ResultsDir <- "/gpfs/data/gtex-group/sex_biased_regulation_v8/sexDE_v8_final/meri/results/allgenes/"

file <- paste0(ResultsDir, "logfc_matrix.sexde.svs.allgenes.txt")
logFC <- read.table(file = file, as.is=T, header=T)

file <- paste0(ResultsDir, "adjpval_matrix.sexde.svs.allgenes.txt")
Qs <- read.table(file = file, as.is=T, header=T)

file <- paste0(ResultsDir, "aveexpr_matrix.sexde.svs.allgenes.txt")
Expr <- read.table(file = file, as.is=T, header=T)

geneInfo <- grch38



