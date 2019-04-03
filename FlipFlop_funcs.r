library(reshape2)

find.common.genes <- function( geneInfo , Qs , logFC, Expr, qThresh = .05){

    ## FIND GENES THAT SHOW SIGNIFICANCE AND DIFFERENT DIRECTIONS IN DIFFERENT TISSUES ##
    logFC.sig <- logFC
    Qs.sig <- Qs
    Expr.sig <- Expr
    ## names(logFC.sig) <- paste("logFC",names(logFC.sig), sep=".")
    orig.cols <- 1:ncol(logFC.sig)
    ind <- Qs >= qThresh
    logFC.sig[ind] <- NA
    Qs.sig[ind] <- NA
    Expr.sig[ind] <- NA
    logFC.sig$N <- rowSums(!is.na(logFC.sig))
    logFC.sig$NposFC <- rowSums(logFC.sig[,orig.cols] > 0, na.rm=T)
    logFC.sig$NnegFC <- rowSums(logFC.sig[,orig.cols] < 0, na.rm=T)

    N <- table(logFC.sig$N)

    ## REMOVE GENES THAT HAVE NO FLIP FLOPS 
    ind.both <- which(logFC.sig$NposFC >0 & logFC.sig$NnegFC > 0)
    logFC.sig <- logFC.sig[ind.both,]
    Qs.sig <- Qs.sig[ind.both,]
    Expr.sig <- Expr.sig[ind.both, ]
    logFC.sig$propPos <- logFC.sig$NposFC / (logFC.sig$N)
    
    logFC.sig$SignMaj <- NA
    ind <- logFC.sig$propPos >.5
    logFC.sig$SignMaj[ind] = 1
    ind <- logFC.sig$propPos <.5
    logFC.sig$SignMaj[ind] = -1

    ind.Maj <- t(apply(logFC.sig[,c(1:44,grep("Sign",names(logFC.sig)))], 1,
                       function(x) sign(x[1:44])==x[45] ))
    ind.Maj[ind.Maj==FALSE] = NA
    ind.Min <- t(apply(logFC.sig[,c(1:44,grep("Sign",names(logFC.sig)))], 1,
                       function(x) sign(x[1:44])!=x[45] ))
    ind.Min[ind.Min==FALSE] = NA

    Qs.sig$mean.q <- rowMeans(Qs.sig[,1:44], na.rm=T)
    Qs.maj <- Qs.sig[,1:44] * (ind.Maj)
    Qs.sig$mean.q.Maj <- rowMeans(Qs.maj, na.rm=T)
    Qs.min <- Qs.sig[,1:44] * (ind.Min)
    Qs.sig$mean.q.Min <- rowMeans(Qs.min, na.rm=T)

    t.test(Qs.sig$mean.q.Maj, Qs.sig$mean.q.Min)

    logFC.sig$mean.logFC<- rowMeans(logFC.sig[,1:44], na.rm=T)
    logFC.maj <- logFC.sig[,1:44] * (ind.Maj)
    logFC.sig$mean.logFC.Maj <- rowMeans(logFC.maj, na.rm=T)
    logFC.min <- logFC.sig[,1:44] * (ind.Min)
    logFC.sig$mean.logFC.Min <- rowMeans(logFC.min, na.rm=T)
    
    t.test(abs(logFC.sig$mean.logFC.Maj), abs(logFC.sig$mean.logFC.Min))

    Expr.sig$mean<- rowMeans(Expr.sig[,1:44], na.rm=T)
    Expr.maj <- Expr.sig[,1:44] * (ind.Maj)
    Expr.sig$mean.Maj <- rowMeans(Expr.maj, na.rm=T)
    Expr.min <- Expr.sig[,1:44] * (ind.Min)
    Expr.sig$mean.Min <- rowMeans(Expr.min, na.rm=T)
    
    t.test(Expr.sig$mean.Maj, Expr.sig$mean.Min)

        
    ## DISTRIBUTION OF TISSUES THAT ARE SPOILERS ##
    tissues.totDE <- data.frame(tissue = names(Qs), analysis = "allgenes",
                                count = colSums(Qs <= 0.05, na.rm=T))
    tissues.withspoiler <- data.frame(tissue = names(Qs.sig[,1:44]), analysis="geneswspoolers",
                                      count = colSums(!is.na(Qs.sig[,1:44])))
    tissues.spoilers <- data.frame(tissue = names(Qs.sig[,1:44]), analysis="spoolergenes",
                                   count = colSums(ind.Min, na.rm=T))
    tissue.count <- rbind(tissues.totDE, tissues.withspoiler, tissues.spoilers)

    tiss.ordered = tissues.totDE$tissue[order(tissues.totDE$count, decreasing=T)]

    ord <-  order(tissues.totDE$count, decreasing=T)
    
    tissue.count$tissue <- factor(tissue.count$tissue, levels = tissues.totDE$tissue[ord])

    tissue.count$countTrim <- sapply(tissue.count$count, function(x) min(1000, x))

    file =  paste0(PlotDir,"TissueDistFlop.pdf")
    pdf(file =file, width=12, height=6)
    ggplot(data=tissue.count, aes(x=tissue, y=countTrim, fill=analysis)) +
        geom_bar(stat="identity", position=position_dodge())+           
            theme_minimal() + 
                theme(text=element_text(size=16),
                      axis.text.x=element_text(size=10,angle=90,hjust=1)) +
                          ylab("Number of Genes") +
                                  xlab("Tissue") 
    dev.off()
    
    
    ## ARE LOWLY EXPRESSED GENES MORE LIKELY TO BE FLIPFLOP GENES ? ##
    ## LOOK AT THE NUMBER OF MINORITY GENES VERSUS MEAN EXPRESSION ##
    no.de.genes <- data.frame(gene = rownames(Qs), count = rowSums(Qs <= qThresh, na.rm=T))

    logFC.sig$gene <- rownames(logFC.sig)
   
    de.ind <- Qs <= qThresh
    de.ind.pos <- de.ind & logFC > 0
    de.ind.neg <- de.ind & logFC < 0
    Expr.de.mean <- data.frame(gene =rownames(Expr), Expr.mean = rowSums(Expr * de.ind) / rowSums(de.ind),
                               Expr.mean.pos = rowSums(Expr * de.ind.pos) / rowSums(de.ind.pos),
                               Expr.mean.neg = rowSums(Expr * de.ind.neg) / rowSums(de.ind.neg) )
    no.de.genes <- merge(no.de.genes, Expr.de.mean, by="gene", all=T)
    
    no.de.genes <- merge(no.de.genes, logFC.sig[, c("gene","N","NposFC","NnegFC")],
                          by="gene", sort=F, all.x=T, all.y=T)

    no.de.genes$min.no <- apply(no.de.genes[,c("NposFC","NnegFC")], 1, min) 

    Q.mean <- data.frame(gene = rownames(Qs), q.mean = rowSums(Qs * de.ind) / rowSums(de.ind),
                         q.mean.pos = rowSums(Qs * de.ind.pos) / rowSums(de.ind.pos), 
                         q.mean.neg = rowSums(Qs * de.ind.neg) / rowSums(de.ind.neg))

    no.de.genes <- merge(no.de.genes, Q.mean, by="gene", all=T)
    
    ## DOTPLOT OF MEAN Q VALUES IN GENES WITH NO FLIP FLOPS AND THOSE WITH FLIP FLOPS
    ## WITH NUMBER OF DE GENES ON X AXIS
    data <- no.de.genes[,c("count","q.mean", "N", "min.no")]
    data <- data[-which(data$count < 2), ]
    data$flipflop <- "flipflop"
    data$flipflop[is.na(data$min.no)] = "unanimous"
    data$logq <- -log10(data$q.mean)
    data$logqtrunc <- sapply(data$logq, function(x) min(x,10))
    
    file = paste0(PlotDir,"QDist_FlipFlop.pdf")
    pdf(file=file, width=12, height = 6)
    p <- ggplot(data=data, aes(x = as.factor(count), y=logqtrunc, fill=flipflop)) +
        geom_boxplot() +
           ## geom_dotplot(binaxis='y', stackdir='center', binwidth=.1) +
                xlab("Number of Tissues Gene is DE in") +
                    ylab("mean(-log10(q-value))") 
    print(p)
    dev.off()
    

    ## LOOK AT PAIRWISE TISSUE DISCORDANCE VERSUS CONCORDANCE ##
    tissues <- colnames(Qs)
    concord <- c()
    FCsig <- logFC * de.ind
    FCsig[FCsig == 0] = NA
    for (i in 1:(length(tissues)-1)){
        
        for (j in (i+1):length(tissues)){
            
            con <- sum(sign(FCsig[,i]) == sign(FCsig[,j]), na.rm=T)
            dis <- sum(sign(FCsig[,i]) != sign(FCsig[,j]), na.rm=T)
            concord <- rbind(concord, c(tissues[i], tissues[j], con, dis))
        }
    }

    concord <- as.data.frame(concord)
    names(concord) <- c("tiss1", "tiss2", "conc", "disc")
    concord[,3:4] <- apply(concord[,3:4], 2, as.numeric)
    concord$propcon <- concord$conc/rowSums(concord[,3:4])

    x <- chisq.test(concord[,3:4], simulate.p.value=T, B=100000)
    concord$stresid <- x$stdres[,1]
    ord <- order(abs(concord$stresid), decreasing=T)
    concord <- concord[ord,]

    ## write out file ##
    file = paste0(PlotDir, "ConcordDiscord.txt")
    write.table(file = file, concord, quote=F, row.names=F, col.names=T, sep="\t")
    
    file = paste0(PlotDir, "ConcordDiscordHeatMap.pdf")
    pdf(file = file)
    p <- ggplot(data = concord[,c("tiss1","tiss2","stresid")],
                aes(x = tiss1, y=tiss2, fill=stresid)) +
                    geom_tile() + xlab("") + ylab("") +
                        ggtitle("Standardized Residulal from Test of Homogeneity") +
                            scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                     midpoint = 0, space = "Lab", 
                                                 name="Std Resid") +
                                                     theme_minimal()+ 
                                                         theme(text = element_text(size=4),
                                                               axis.text.x = element_text(angle = 45,
                                                                   vjust = 1, hjust = 1))+
                                                                       coord_equal()
    
    print(p)
    dev.off()

    ## What proportion of significant resulsts are spoilers. How close is it to 0.05?
    N = sum(no.de.genes$count)
    Nspoiler = sum(no.de.genes$min.no, na.rm=T)
    Ngt2 <- sum(no.de.genes$count[no.de.genes$count>=2])
    

    ## repeat analysis for various qthresh ##

    qThreshes <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)

    Out <- list()

    for (qThresh in qThreshes){

        print(qThresh)
        tmp <- summarize.flipflop(Qs, logFC, Expr, GeneInfo, qThresh)
        Out[[as.character(qThresh)]] <- tmp

    }

    d <- c()
    for (qThresh in qThreshes){

        tmp <- Out[[as.character(qThresh)]]
        Ngenes <- sum(tmp$count>0)
        Ngenesgt1tiss <- sum(tmp$count >= 2)
        Ngenesflipflop <- sum(tmp$Npos >0 & tmp$Nneg>0)
        Nflipflop <- sum(tmp$min.no)
        Ngt1tiss <- sum(tmp$count[tmp$count>1])
        N <- sum(tmp$count)

        d <- rbind(d, c(qThresh, Ngenes, Ngenesgt1tiss, Ngenesflipflop, N, Ngt1tiss, Nflipflop))
    }
    d <- as.data.frame(d)
    names(d) <- c("qThresh", "Ngenes","Ngenesgt1is", "Ngenesflipflop", "NDE", "Ngt1tiss", "Nflipflop")
    dmelt <- melt(d, id="qThresh", value.name="count", variable.name="genetype")

    d$propFFgenes <- d$Ngenesflipflop / d$Ngenesgt1is
    d$propFFresults <- d$Nflipflop / d$Ngt1tiss
    
    ## PLOT RELATIONSHIP BETWEEN Q THRESHOLD AND NUMBER OF GENES DE, NUMBER OF THOSE THAT ARE UNANOMOUS
    ## AND THE NUMBER THAT ARE FLIPFLOPERS

    file <- paste0(PlotDir, "QvDEgenesFlipFlop.pdf")
    a <- dmelt[dmelt$genetype %in% c("Ngenes", "Ngenesgt1is", "Ngenesflipflop"),]
    pdf(file = file)
    p <- ggplot(a, aes(x = qThresh, y=count, color=genetype)) +
        geom_point() + 
            geom_smooth(method=loess, se=FALSE, fullrange=TRUE) + 
                ylab("Number of genes") + xlab("Q value significance threshold") +
                    scale_color_discrete(name="Number of",
                                         labels=c("DE Genes", "DE Genes > 1 Tis", "Flipflop Genes"))
    print(p)

    
    
    p <- ggplot(d, aes(x=qThresh, y=d$propFFgenes)) +
         geom_point() + 
             geom_smooth(method=loess, se=FALSE, fullrange=TRUE) + xlab("Q value significance threshold") +
                 ylab("Proportion of Flip Flop Genes")
    print(p)
    
    a <- dmelt[dmelt$genetype %in% c("NDE", "Ngt1tiss", "Nflipflop"),]
    p <- ggplot(a, aes(x = qThresh, y=count, color=genetype)) +
        geom_point() + 
            geom_smooth(method=loess, se=FALSE, fullrange=TRUE) + 
                ylab("Number of DEs") + xlab("Q value significance threshold") +
                    scale_color_discrete(name="Number of",
                         labels=c("DEs", "DE in Genes > 1DE", "Flipflop DE"))
    print(p)

    p <- ggplot(d, aes(x=qThresh, y=d$propFFresults)) +
        geom_point() + 
            geom_smooth(method=loess, se=FALSE, fullrange=TRUE) + xlab("Q value significance threshold") +
                ylab("Proportion of Flip Flop DEs")
    print(p)
    
    dev.off()
}

summarize.flipflop <- function(Qs, logFC, Expr, GeneInfo, qThresh){

    no.de.genes <- data.frame(gene = rownames(Qs), count = rowSums(Qs <= qThresh, na.rm=T))
    
    de.ind <- Qs <= qThresh
    de.ind.pos <- de.ind & logFC > 0
    de.ind.neg <- de.ind & logFC < 0

    no.de.genes$NposFC <- rowSums(de.ind.pos)
    no.de.genes$NnegFC <- rowSums(de.ind.neg)
    no.de.genes$min.no <- apply(no.de.genes[,c("NposFC", "NnegFC")], 1, min)
    no.de.genes$min.no[no.de.genes$min.no == 0] <- NA
    
    Expr.de.mean <- data.frame(gene =rownames(Expr), Expr.mean = rowSums(Expr * de.ind) / rowSums(de.ind),
                               Expr.mean.pos = rowSums(Expr * de.ind.pos) / rowSums(de.ind.pos),
                               Expr.mean.neg = rowSums(Expr * de.ind.neg) / rowSums(de.ind.neg) )
    no.de.genes <- merge(no.de.genes, Expr.de.mean, by="gene", all=T)
    
    no.de.genes$min.no <- apply(no.de.genes[,c("NposFC","NnegFC")], 1, min) 
    
    Q.mean <- data.frame(gene = rownames(Qs), q.mean = rowSums(Qs * de.ind) / rowSums(de.ind),
                         q.mean.pos = rowSums(Qs * de.ind.pos) / rowSums(de.ind.pos), 
                         q.mean.neg = rowSums(Qs * de.ind.neg) / rowSums(de.ind.neg))
    
    no.de.genes <- merge(no.de.genes, Q.mean, by="gene", all=T)

    return(no.de.genes)
}
