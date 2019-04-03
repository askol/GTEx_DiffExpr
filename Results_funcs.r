## FUNCTIONS TO HELP SUMMARIZE THE RESULTS FROM GTEX
## DE ANALYSIS

ResultDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/"
PlotDir <- paste0(ResultDir,"Plots/")
GSEAOutDir <- paste0(ResultDir,"GSEA/")
GSEAScript <- paste0(GSEADir,"Script/")
GSEARankDir <- paste0(GSEADir,"Ranks/")

source("~/Code/qq_plot.r")
library(calibrate)

## CREATE PLOT DIRECTORY IF IT DOESN'T EXIST ##
if (!dir.exists(PlotDir)){
    dir.create(PlotDir)
}

load.results <- function(tissue){

    file <- paste0(ResultDir, tissue, "_DE_limma.rslts")

    if (!file.exists(file)){
        print(paste0("Results file", file, "not found!"))
        return(0)
    }
    rslt <- read.table(file = file, header=T, as.is=T)

    return(rslt)
}

qqplot <-function(tissue, rslt, PlotDir){

    file <- paste0(PlotDir,"QQ_",tissue,".pdf")
    title.txt <- paste0(tissue)
    p <- qq_plot(rslt$P.Value, title.txt)

    p.wox <- qq_plot(title="X removed", rslt$P.Value[rslt$chr != "X"])
    
    pdf(file = file)
    print(p)
    print(p.wox)
    dev.off()
    print(paste0("Created plot ",file))
}

volcano <-function(tissue, rslt, PlotDir, pThresh=.05, FCThresh=1, title=""){    

    file <- paste0(PlotDir, "volcano_", tissue, ".pdf")
    pdf(file = file)
    
    ## Make a basic volcano plot
    with(rslt, plot(logFC, -log10(P.Value), pch=20,
                   main="Volcano plot", sub=title, xlim=c(-2.5,2)))

    ## Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(rslt, adj.P.Val < pThresh ),
         points(logFC, -log10(P.Value),
                pch=20, col="red"))
    with(subset(rslt, abs(logFC) > FCThresh),
         points(logFC, -log10(P.Value),
                pch=20, col="orange"))
    with(subset(rslt, adj.P.Val < pThresh & abs(logFC) > FCThresh),
         points(logFC, -log10(P.Value), pch=20, col="green"))

    ## Label points with the textxy function from the calibrate plot

    with(subset(rslt, adj.P.Val < pThresh & abs(logFC) > FCThresh),
         textxy(logFC, -log10(P.Value), labs=SYMBOL, cex=.8))
    
    dev.off()
    
    print(paste0("Created plot ",file))
}



run.gsea <- function(tissue, rslt, GSEARankDir, GSEAOutDir){
    
    geneSet <- "msigdb.v6.2.symbols.gmt"
    rankFile <- paste0(GSEARankDir, tissue, ".rnk")
    rnk <- -1
    rnks <- make.ranks(rslt)
    write.table(file = rankFile, rnks,
                quote=F, row.names=F, col.names=F, sep="\t")
    
    print(paste("Working on ",geneSet," and ",rankFile))
        
    geneSetName <- gsub("\\.v6.+","",geneSet)
    rpt.lab <- paste0(tissue,"_",geneSetName)
    params <- rbind(c("rnk", rankFile),
                    c("out", rpt.lab))
    param.file <- paste0("~/GSEA_parameters_", tissue, ".txt")
    write.table(file = param.file, params, quote=F, col.names=F,
                row.names=F, sep="\t")
    cmd <- paste0("java -cp /home/askol/bin/gsea-3.0.jar -Xmx5000m ",
                  "xtools.gsea.GseaPreranked ",
                  "-param_file ", param.file,
                  " -gmx /home/askol/bin/GSEA_genesets/",
                  "msigdb_v6.2_GMTs/",geneSet," -norm meandiv -nperm 1000 ",
                  "-scoring_scheme weighted -rpt_label ", rpt.lab,
                  " -create_svgs false -make_sets true -plot_top_x 100 ",
                  "-rnd_seed timestamp ",
                  "-set_max 500 -set_min 15 -zip_report false  -gui false")
    print(cmd)
    system(cmd)
}

    
    
make.ranks <- function(rslt){

    rnks <- sign(rslt$logFC - 1)*log10(rslt$P.Value)
    rnks <- cbind(rslt$SYMBOL, rnks)
    rnks <- rnks[order(as.numeric(rnks[,2]), decreasing=T), ] 
    return(rnks)
}
