##  Create QQ plots, Volcano Plots, HeatMap of sign of logFC,
##  tSNPE or PCA Cluster based on sign of log FC
##  SHARED SIGNIFICANT GSs
library(Rtsne)
library(FactoMineR)
library(factoextra)
library(qvalue)

setwd("/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/Summary/")
GSEADir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/GSEA/"
ResultDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/"
SummaryDir <- paste0(ResultDir,"Summary/")
PlotDir <- paste0(ResultDir,"Plots/")

GSEAOutDir <- paste0(ResultDir,"GSEA/")
GSEAGSFiles <- c("/home/askol/bin/GSEA_genesets/Custom/Hormone_Immune_Custom.gmt",
                "/home/askol/bin/GSEA_genesets/msigdb_v6.2_GMTs/msigdb.v6.2.symbols.gmt")
    
source("/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Code/Summary_funcs.r")

tissues <- read.table(file = paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "data/support_files/all_v8_tissues_both_sexes.txt"), as.is=T)
tissues <- unlist(tissues)
tissues.keep <- c()
Ns <- c()
for (tissue in tissues){
    samp.size <- check.sample.size(tissue)
    
    if (min(unlist(samp.size)) < 40){next }
    tissues.keep <- c(tissues.keep, tissue)
    Ns <- rbind(Ns, c(tissue, samp.size, round(mean(unlist(samp.size))),
                      round(sum(sqrt(unlist(samp.size))))))
}
Ns <- as.data.frame(Ns, stringAsFactors=FALSE)
names(Ns) <- c("tissue", "males", "females", "mean", "mean.sqrt")
tissues <- tissues.keep

tmp <- collect.results(tissues)

logFC <- tmp$logFC
logPs <- tmp$logPs
Qs <- tmp$Qs
geneInfo <- tmp$GeneInfo

plot.DE.by.tissue(Qs, PlotDir, Ns)

plot.tSNE(logPs, PlotDir)

plot.PCA(logPs, PlotDir)


## WHICH GENES ARE SHARED BETWEEN TISSUES ##
common.genes <- find.common.genes(geneInfo, Qs, logFC, tissues,
                                  qThresh = .05, OutDir = SummaryDir)

## UNDERSTAND LD BETWEEN MOST DE GENES ##
genes <- common.genes$DEgenes$SYMBOL
Cor <- calc.corr.btw.genes(genes, tissues, DataDir, geneInfo)
CorPlot <- 
plot.cor(Cor, PlotDir)

## Collect GSEA info for all tissues ##
Gs <- get.gsea.results(tissues, GSEADir)

gs.nox <- plot.GS(Gs$Gsnox, PlotDir, Ns, qThresh = .25, plotsufx="_withoutX")
gs.x <- plot.GS(Gs$Gswx, PlotDir, Ns, qThresh = .25, plotsufx="_withX")

## figure out which GSs show up in more than one tissue ##
find.shared.gs(Gs$Gsnox, SummaryDir, "Shared_GS_noXgene.txt", nohits = 2)
find.shared.gs(Gs$Gswx, SummaryDir, "Shared_GS_withXgene.txt", nohits = 2)




    
        
    
