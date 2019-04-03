## FUNCTIONS TO HELP SUMMARIZE THE RESULTS FROM GTEX
## DE ANALYSIS


DataDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Data/"

source("~/Code/qq_plot.r")
library(calibrate)
library(WGCNA) ## for the faster cor function
library(reshape2)
library(RColorBrewer)


check.sample.size <- function(tissue){

    OrigDataDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "sexDE_v8_final/meri/data/")

    count_file <- paste0(OrigDataDir, "Phenotypes/",tissue,
                         ".v8.gene_counts.txt.gz")
    covs_file <- paste0(OrigDataDir, "/Covariates/covs.basalcovs.",tissue,
        ".txt")
    sex.info <- read.table(file = covs_file, as.is=T, header=F, quote="",
                           comment.char = "", sep="\t")
    colnames(sex.info) = c('SUBJID', 'SEX','SMTSISCH','SMRIN','AGE')
    samps <- read.table(file = count_file, nrow=1, header=F, as.is=T)
    samps <- gsub("\\.","-",samps)
    
    males <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="M"]]
    females <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="F"]]

    wcs <- list(males = length(males), females=length(females))

    return(wcs)
}

load.results <- function(tissue, ResultDir){

    file <- paste0(ResultDir, tissue, "_DE_limma.rslts")

    if (!file.exists(file)){
        print(paste0("Results file", file, "not found!"))
        return(0)
    }
    rslt <- read.table(file = file, header=T, as.is=T)

    return(rslt)
}


get.gsea.results <- function(tissues, GSEADir){

    geneSets <- c("msigdb.v6.2.symbols.gmt" , "Hormone_Immune_Custom.gmt")

    Gs <- c()
    for (tissue in tissues[1:41]){

        G <- c() ## 
        for (geneSet in geneSets){

            for (suff in c("", "_nox")){
                
                GSEAOutDir <- paste0(GSEADir, gsub("-","_",tissue), "/")
                file <- paste0(GSEAOutDir, tissue, "_", geneSet,
                               suff,"_GSEA_summary.txt")
                
                if (!file.exists(file) | file.size(file)==0){
                    print(paste0("No significant gene sets for ",tissue, " in geneset ",
                                 geneSet))
                    next
                }
                
                gs <- read.table(file = file, as.is=T, sep="\t", header=T)
                gs$nox <- NA
                if (suff == ""){
                    gs$nox <- "withX"
                }else{
                    gs$nox <- "noX"
                }
                
                gs <- gs[,c("nox","NAME","FDR.q.val")]

                names(gs)[grep("FDR",names(gs))] <- tissue

                G <- rbind(G, gs)
            }
        }
        if (length(Gs)==0){
            
            Gs <- G
        }else{
            
            if (length(G) > 0){
                Gs <- merge(Gs, G, by=c("NAME", "nox"), all=T)
            }           
        }        
    }

    Gsnox <- Gs[Gs$nox == "noX",]
    Gswx <- Gs[Gs$nox == "withX",]
    
    return(list(Gsnox = Gsnox, Gswx=Gswx))
}

find.common.gs <- function(Gs, geneInfo, GSEAGSFiles){

    gs.gs.count <- data.frame(NAME = Gs$NAME,
                              tissue.count = rowSums(!is.na(Gs[,-1])))
    gs.gs.count <- gs.gs.count[order(gs.gs.count$tissue.count, decreasing=T),]

    propx <- get.prop.X(gs = gs.gs.count$NAME, geneInfo, GSEAGSFiles)
    gs.gs.count <- merge(gs.gs.count, propx, by="NAME", all=T, sort=F)
    
    tissue.gs.count <- data.frame(tissue = colnames(Gs)[-1],
                                  count = colSums(!is.na(Gs[,-1])))
    tissue.gs.count <- tissue.gs.count[order(tissue.gs.count$count, decreasing=T),]
                                       
    return(list(gs.gs.count = gs.gs.count, tissue.gs.count = tissue.gs.count))
}

get.prop.X <- function(gs, geneInfo, GSEAGSFiles){

    gs <- as.character(gs)
    genes <- list()
    GS <- c()
    for (file in GSEAGSFiles){

        g <- read.table(file = file, as.is=T, header=F, sep="\t", fill=T)
        g <- g[g[,1] %in% gs,]
        if (nrow(g) == 0){ next }
        for (i in 1:nrow(g)){

            genes <- g[i,-c(1:2)]
            genes <- genes[-which(genes == "")]

            ind <- which(geneInfo$SYMBOL %in% genes)
            prop.x <- mean(geneInfo$chr[ind]=="X")

            GS <- rbind(GS, c(g[i,1], length(ind), prop.x))

        }
    }
    GS <- data.frame(NAME = GS[,1], gene.count = GS[,2],
                     propX = GS[,3], stringsAsFactors=FALSE)

    return(GS)
}
    
            

               
            
get.data <- function(tissue, DataDir){

    file <- paste0(DataDir, tissue, ".RDS")
    data <- readRDS(file)

    txt <- paste0("model.matrix(~SEX + AGE + ",
                  paste(names(data$samples)[grep("SV", names(data$samples))], collapse="+"),
                  ", data$samples)")
    design <- eval(parse(text = txt))
    
    v <- voom(data, design, plot=F)

    return(list(design = design, count = data, v = v))
}
    

collect.results <- function(tissues, ResultDir){
    
    logFC <- c()
    logPs <- c()
    logQs <- c()
    geneInfo <- c()
    for (tissue in tissues){
        
        print(paste0("Working on tissue : ",tissue))
        rslt <- load.results(tissue, ResultDir)
        dupe.ind <- which(duplicated(rslt$SYMBOL))    
        rslt <- rslt [-dupe.ind,]
        print(paste0("Removing ",length(dupe.ind), " duplicate gene"))
        logFCtmp <- rslt[,c("SYMBOL", "logFC")]
        logPstmp <- rslt[,c("SYMBOL", "P.Value")]
        names(logFCtmp)[2] = names(logPstmp)[2] = tissue        

        gi <-  rslt[, c("gene","SYMBOL", "chr", "start","end","gene_biotype")]
        if (tissue == tissues[1]){
            logFC <- logFCtmp
            logPs <- logPstmp
            geneInfo <- gi
        }else{            
            logFC <- merge(logFC, logFCtmp, by="SYMBOL", all=T)
            logPs <- merge(logPs, logPstmp, by="SYMBOL", all=T)
            geneInfo <- rbind(geneInfo, gi[!gi$SYMBOL %in% geneInfo$SYMBOL,])
        }
    }

    ## CREATE FDR ##
    Qs <- make.Qs(logPs)

    ## take log10 of Ps
    logPs[,-1] = -log10(logPs[,-1])

    geneInfo <- geneInfo[match(logFC$SYMBOL, geneInfo$SYMBOL),]
    return(list(logFC = logFC, logPs = logPs, Qs = Qs, GeneInfo = geneInfo))
}

make.Qs <- function(logPs){
    Qs <- logPs
    for (i in 2:ncol(logPs)){
        p <- logPs[,i]
        miss.ind <- which(is.na(p))
        p <- p[-miss.ind]
        small.ind<- which(p < 10^-30)
        p[small.ind] <- 10^-30
        q <- qvalue(p)$qvalues
        Qs[-miss.ind,i] <- q
    }
    return(Qs)
}

plot.DE.by.tissue <- function(Qs, PlotDir, Ns){

    sign <- Qs
    sign[,-1] <- 1*(sign[,-1]<= .05)
    tissue.total <- colSums(sign[,-1], na.rm=T)
    tissue.total <- data.frame(tissue = names(tissue.total), count = tissue.total)
    ord <- order(as.numeric(as.character(Ns$mean.sqrt)), decreasing=T)
  
    tissue.total$tissueLab <- paste0(Ns$tissue, " (",paste(Ns$male,Ns$female,sep="/"),
                                     ")")[match(tissue.total$tissue, Ns$tissue)]
    tissue.total$tissueLab <- factor(tissue.total$tissueLab,
                                     levels = tissue.total$tissueLab[ord])
    gene.total <- rowSums(sign[,-1], na.rm=T)
    gene.total <- data.frame(gene = Qs$SYMBOL, count = gene.total)

    file <- paste0(PlotDir, "Sex_DE_by_Tissue.pdf")
    pdf(file=file, width=10, height=5)

    ## NUMBER OF DE GENES BY TISSUE
    p = ggplot(data = tissue.total, aes(x=tissueLab, y=count)) +
    geom_bar(stat="identity", fill="steelblue")+
        theme_minimal() +
    theme(text=element_text(size=10),
          axis.text.x=element_text(angle=60, hjust=1)) +
    ylab("Number Sex DE Genes (q-value < 0.05)") +
    xlab("Tissue")
    print(p)

    ## DISTRIBUTION OF THE NUMBER OF TISSUES FOR WHICH A GENE IS DE ##
    tbl <- table(gene.total$count)
    tbl <- data.frame(noTiss = names(tbl), noGenes=as.vector(tbl))
    for (i in 1:2){
        if (i == 1){ ylim = c(0, max(tbl$noGenes)) } else {
            ylim = c(0,25)
        }
        
        p = ggplot(data = gene.total, aes(x=count)) +
            geom_histogram(binwidth=1, color="black", fill="white") +
                theme_minimal() +
                    theme(text=element_text(size=16),
                          axis.text.x=element_text(angle=90,hjust=1)) +
                              ylab("Number of Genes") +
                                  xlab("Gene is DE in _ Tissues") +
                                      ylim(ylim)
     
        print(p)
    }
    dev.off()

    print(paste0("Wrote plots to file: ",file))
}
    
plot.tSNE <- function(logPs, PlotDir){
    
    ## TSNE ##
    train <- t(logPs)
    colnames(train) <- train[1,]
    train <- train[-1,]
    train <- data.frame(train)
    train <- sapply(train, function(x) as.numeric(as.character(x)))
    train <- data.frame(tissue = names(logPs)[-1], train)
    Labels <- train$tissue
    train$tissue<-as.factor(train$tissue)
    
    ## for plotting
    colors = rainbow(length(unique(train$tissue)))
    names(colors) = unique(train$tissue)
    
    ## remove genes with missingness ##
    rm.col <- which(colSums(is.na(train)) != 0)
    train <- train[,-rm.col]
    
    ## Executing the algorithm on curated data
    tsne <- Rtsne(train[,-1], dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)
    exeTimeTsne<- system.time(Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500))
    
    ## Plotting
    file = paste0(PlotDir, "TSNE_DE.pdf")
    pdf(file = file)
    plot(tsne$Y, t='n', main="tsne")
    text(tsne$Y, labels=train$tissue, col=colors[train$tissue], cex=.4)
    dev.off()
    print(paste0("Created plot in ",file))

}

plot.PCA <- function(logPs, PlotDir){
    
    ## LETS LOOK VIA PCA ##
    ## PCA and Corplot from here
    ## http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
    rownames(train) <- train$tissue
    res.pca <- PCA(train[,-1], scale.unit=TRUE, graph = FALSE)
    
    eig.val <- get_eigenvalue(res.pca)
    eig.val
    
    var <- get_pca_var(res.pca)
    var
    
    file = paste0(PlotDir, "PC_DE_logP.pdf")
    pdf(file)
    
    fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(1,2),
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE, # Avoid text overlapping (slow if many points),
                 labelsize = 2)
    
    fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(2,3), 
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE ,
                 labelsize=2)
    
    fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(3,4), 
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE ,
                 labelsize=2)
    
    fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(4,5),
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE ,
                 labelsize=2)
    
    dev.off()
    print(paste0("Wrote plots to file ",file))
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
 

get.gsea <- function(tissues, GSEADir){

    geneSets <- c("msigdb.v6.2.symbols", "Hormone_Immune_Custom")

    Q <- NES <- c()
    
    for (tissue in tissues){
        
        gsea <- c()

        for (geneSet in geneSets){
        
            file <- paste0(GSEAOutDir, tissue, "/",tissue,"_",
                           geneSet, ".gmt_GSEA_summary.txt")
            if (file.info(file)$size == 0){
                print(paste0("Nothing in file ",file))
                next
            }
            gsea <- read.table(file = file, as.is=T, header=T, sep="\t")
        }
        
        q <- gsea[,c("NAME","FDR.q.val")]
        nes <- gsea[,c("NAME","NES")]
        names(q)[2] <- names(nes)[2] <- tissue
        if (length(Q) == 0){

            Q <- q
            NES <- nes
        }else{
            Q <- merge(Q, q, by="NAME", all=T)
            NES <- merge(NES, nes, by="NAME", all=T)
        }
    }

    return(list(Q = Q, NES = NES))

}

find.common.genes <- function( geneInfo , Qs , logFC, tissues, qThresh = .05, OutDir){

    count <- rowSums(Qs[,-1] < qThresh, na.rm=T)

   

    ## FIND GENES THAT SHOW SIGNIFICANCE AND DIFFERENT DIRECTIONS IN DIFFERENT TISSUES ##
    logFC.sig <- logFC[,-1]
    names(logFC.sig) <- paste("logFC",names(logFC.sig), sep=".")
    
    ind <- Qs[,-1] >= qThresh
    logFC.sig[ind] <- NA
    logFC.sig$NposFC <- rowSums(logFC.sig > 0, na.rm=T)
    logFC.sig$NnegFC <- rowSums(logFC.sig < 0, na.rm=T)
    logFC.sig$N <- rowSums(!is.na(logFC.sig))
    ind.both <- which(logFC.sig$NposFC >0 & logFC.sig$NnegFC > 0)
    logFC.sig <- logFC.sig[ind.both,]
    logFC.sig$propPos <- logFC.sig$NposFC / (logFC.sig$NposFC + logFC.sig$NnegFC)

    Qs.tmp <- Qs[, -1]; names(Qs.tmp) <- paste("q",names(Qs)[-1], sep=".")
    Qs.tmp[ind] <- NA
    
    Flop <- data.frame(geneInfo[ind.both,c("SYMBOL","chr","start")],
                       logFC.sig,
                       Qs.tmp[ind.both,-1])
    file <- paste0(OutDir,"FlipFlopDEGenes.txt")
    write.table(file = file, Flop, quote=F, row.names=F, col.names=T)
    
    rm(Qs.tmp)
    
    props <- count / rowSums(!is.na(Qs[,-1]))
    ind <- which(count >= 10 | (props >= 10/42 & count >= 5))
    tbl <- table(geneInfo$chr[ind])

    chrs <- as.data.frame(tbl)
    names(chrs)[1] <- "chr"

    tots <- as.data.frame(table(geneInfo$chr))
    names(tots) <- c("chr", "total")

    chrs <- merge(chrs, tots, by="chr")
    chrs$exp <- sum(chrs$Freq)*chrs$total/sum(chrs$total)
    chrs$excess <- chrs$Freq / chrs$exp
    chrs <- chrs[order(chrs$excess, decreasing=T),]

    ind.x <- which(chrs$chr == "X")
    chrs$exp.no.x <- NA
    chrs$exp.no.x[-ind.x] <- sum(chrs$Freq[-ind.x])*
        chrs$total[-ind.x]/sum(chrs$total[-ind.x])
    chrs$excess.no.x <- NA
    chrs$excess.no.x[-ind.x] <- chrs$Freq[-ind.x] / chrs$exp.no.x[-ind.x]

    ## Create table of genes sorted by chromosome and q ##
    DEgenes <- geneInfo[ind,]
    DEgenes <- merge(DEgenes, Qs, by="SYMBOL", all.x=T, all.y=F)

    DEgenes <- DEgenes[DEgenes$chr %in% chrs$chr[chrs$excess>1],]
    DEgenes$q.median <- apply(DEgenes[,-c(1:6)], 1, median, na.rm=T)
    DEgenes$q.median.sig <- apply(DEgenes[,-c(1:6)], 1, function(x)
                                  median(x[which(x < qThresh)], na.rm=T))

    ## ADD INFORMATION ABOUT DIRECTION OF FC ##
    fc <- logFC[ind, -1]
    fc <- fc[,-1]
    fc[fc > qThresh] <- NA
    fc$SYMBOL <- logFC$SYMBOL[ind]
    fc$prop.fc.gt0 <- apply(fc[,-1], 1, function(x) mean(x > 0, na.rm=T))
    fc <- fc[,c("SYMBOL", "prop.fc.gt0")]
    DEgenes <- merge(DEgenes, fc, by="SYMBOL")
    
    ord <- order(DEgenes$chr, DEgenes$q.median.sig)
    DEgenes <- DEgenes[ord, colnames(DEgenes) %in% tissues == F ]


    file <- paste0(OutDir,"MostDEGene.txt")
    write.table(file = file, DEgenes, quote=F, row.names=F, col.names=T)

    print(paste0("Wrote file genes DE in the most tissues in ",file))

    file <- paste0(OutDir, "DEGene_Chr_Dist.txt")
    write.table(file = file, chrs, quote=F, row.names=F, col.names=T)

    print(paste0("Wrote file show distribution of DE genes by chromosome in ",file))

    return(list(chrs = chrs, DEgenes = DEgenes))
}



inverse_quantile_normalization <- function(gct) {
    gct = t(apply(gct, 1, rank, ties.method = "average"));
    gct = qnorm(gct / (ncol(gct)+1));
    return(gct)
}


calc.corr.btw.genes <- function(genes, tissues, DataDir, geneInfo){

    Cor <- c()
    
    for (tissue in tissues){

        print(paste0("Working on tissue ",tissue))
        
        tmp <- get.data(tissue, DataDir)
        v <- tmp$v
        count <- tmp$count
    
        logcpm <- inverse_quantile_normalization(v$E)
    
        logcpm <- data.frame(v$genes, logcpm)
        logcpm <- logcpm[-which(duplicated(logcpm$SYMBOL)),]
        rownames(logcpm) <- logcpm$SYMBOL
        logcpm <- logcpm[, grep("GTEX", names(logcpm))]
        logcpm <- t(logcpm[rownames(logcpm) %in% genes,])

        C <- cor(logcpm)

        Cm <-melt(C, value.name="cor")
        names(Cm)[1:2] <- c("gene1","gene2")
        Cm[,-1] <- apply(Cm[,-1], 2, as.character)
        ##Cm$names <- apply(Cm, 1, function(x) paste(sort(c(x[-3])), collapse="."))
        Cm$names <- apply(Cm, 1, function(x) paste(c(x[-3]), collapse="."))
        ##Cm <- Cm[-which(duplicated(Cm$names)),]
        ##Cm <- Cm[-which(Cm$gene1 == Cm$gene2),]
        Cm$cor <- as.numeric(Cm$cor)
        names(Cm)[3] <- tissue
        if (tissue == tissues[1]){
            Cor <- Cm
        }else{
            Cor <- merge(Cor, Cm[,-c(1:2)], by="names", all=T)
        }
        
    }

    ## Cor$gene1 <- gsub(".+\\.","", Cor$names)
    ## Cor$gene2 <- gsub("\\..+","",Cor$names)
    
    Cor <- merge(Cor, geneInfo[,c("SYMBOL","chr","start")], by.x="gene1",
               by.y="SYMBOL")
    names(Cor)[ncol(Cor) - c(1,0)] <- c("chr.gene1", "pos.gene1")
    Cor <- merge(Cor, geneInfo[,c("SYMBOL","chr","start")], by.x="gene2",
               by.y="SYMBOL")
    names(Cor)[ncol(Cor) - c(1,0)] <- c("chr.gene2", "pos.gene2")

    Cor$median <- apply(Cor[,-c(grep("name",names(Cor)), grep("gene",names(Cor)))],
                        1, median, na.rm=T)

    return(Cor)
    
}

plot.cor <- function(Cor, PlotDir){

    names(Cor) <- gsub("-","_",names(Cor))
    file <- paste0(PlotDir, "Cor_Shared_Genes.pdf")

    ord <- order(Cor$chr.gene1, Cor$pos.gene1)
    genes.ord <- unique(Cor$gene1[ord])
    Cor$gene1 <- factor(Cor$gene1, levels = genes.ord)
    Cor$gene2 <- factor(Cor$gene2, levels = genes.ord)
    pdf(file = file)
    
    ##    hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')
    hm.palette <- colorRampPalette(rev(brewer.pal(9, 'RdBu')), space='Lab')
    p <- ggplot(Cor, aes(x = gene1, y = gene2, fill = median)) +
        geom_tile() + xlab("") + ylab("") + ggtitle("Median Correlation") +
            ## coord_equal() +
                scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                     midpoint = 0, limit = c(-1,1), space = "Lab", 
                                     name="Pearson\nCorrelation") +
                                         theme_minimal()+ 
                                             theme(text = element_text(size=4),
                                                   axis.text.x = element_text(angle = 45,
                                                       vjust = 1, hjust = 1))+
                                                   coord_equal()
    
    print(p)
    
    for (i in 4:45){

        tissue <- names(Cor)[i]
        hm.palette <- colorRampPalette(rev(brewer.pal(9, 'RdBu')), space='Lab')
        eval(parse(text = paste("p <- ggplot(Cor, aes(x = gene1, y = gene2, fill = ",
                       tissue,"))"))) 
        p <- p + geom_tile() + xlab("") + ylab("") + ggtitle(tissue) +
                scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                     midpoint = 0, limit = c(-1,1), space = "Lab", 
                                     name="Pearson\nCorrelation") +
                                         theme_minimal()+ 
                                             theme(text = element_text(size=4),
                                                   axis.text.x = element_text(angle = 45,
                                                       vjust = 1, hjust = 1)) +
                                                   coord_equal()
        print(p)
    }
        
    print(paste0("Wrote plot ",file))
    dev.off()

    g <- Cor[,c("gene1","chr.gene1","pos.gene1")]
    g <- g[-which(duplicated(g$gene1)),]
    g <- g[order(g[,2], g[,3]),]
}

plot.GS <- function(GS, PlotDir, Ns, qThresh=0.25, plotsufx = ""){

    tissue.counts <- colSums(!is.na(GS[,-c(1:2)]))
    tissue.counts <- data.frame(tissue = names(tissue.counts),
                                counts = tissue.counts)
    gs.counts <- rowSums(!is.na(GS[,-c(1:2)]))
    gs.counts <- data.frame(gs = GS$NAME, counts = gs.counts)

    tissue.counts <- merge(tissue.counts, Ns, by="tissue")
    tissue.counts$tissueLab <- paste0(tissue.counts$tissue,
                                      " (",paste(tissue.counts$male,
                                                 tissue.counts$female,sep="/"),
                                      ")")
    ord <- order(as.numeric(as.character(tissue.counts$mean.sqrt)), decreasing=T)
    tissue.counts$tissueLab <- factor(tissue.counts$tissueLab,
                                      levels = tissue.counts$tissueLab[ord])

    file <- paste0(PlotDir, "GeneSet_Distribution",plotsufx,".pdf")
    
    pdf(file=file, width=10, height=5)

    ## NUMBER OF GENESETS BY TISSUE
    p = ggplot(data = tissue.counts, aes(x=tissueLab, y=counts)) +
    geom_bar(stat="identity", fill="steelblue")+
        theme_minimal() +
    theme(text=element_text(size=10),
          axis.text.x=element_text(angle=60, hjust=1)) +
    ylab(paste0("Number Genesets (q-value < ",qThresh,")")) +
    xlab("Tissue")

    print(p)

    p = ggplot(data=gs.counts, aes(x=counts)) +
        geom_histogram(color="black", fill="white") +
            xlab("Number of tissues sharing a Geneset") +
                ylab("Number of Genesets") 
    print(p)

    dev.off()

    print(paste0("Wrote ",file))

    return(list(tissue.counts = tissue.counts, gs.count = gs.counts))
}

find.shared.gs <- function(g, SummaryDir, file.name="Shared_GS.txt",
                           nohits=2){
    noh <- rowSums(!is.na(g[,-c(1:2)]))
    ind <-which(noh >= nohits)
    noh <- noh[ind]
    tmp <- apply(g[ind,-c(1:2)], 1, function(x) paste(names(x)[!is.na(x)],
                                                      collapse=";"))
    GSshare <- cbind(g$NAME[ind], tmp, noh)
    ord <- order(noh, decreasing=T)
    GSshare <- GSshare[ord,]
    
    file <- paste0(SummaryDir, file.name)
    write.table(file = file, GSshare, quote=F, row.names=F, col.names=F)

    print(paste0("Wrote file ",file))

}
