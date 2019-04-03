ResultDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/"
PlotDir <- paste0(ResultDir,"Plots/")
GSEAOutDir <- paste0(ResultDir, "GSEA/")
RankDir <- paste0(GSEAOutDir, "Ranks/")

source("/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Code/Results_funcs.r")

run_results <- function(tissue){

    rslt <- load.results(tissue)
    if (length(rslt)!= 1){ 
        qqplot(tissue, rslt, PlotDir)
        
        volcano(tissue, rslt, PlotDir)

        volcano(tissue, rslt[rslt$chr!="X",], PlotDir, title = "X Removed")
        print(paste0("Processed tissue ",tissue))

        run.gsea(tissue, rslt, RankDir, GSEAOutDir)
    }

    
}


args <- commandArgs(TRUE)
tissue = args[1]
run_results(tissue)
print(date())

