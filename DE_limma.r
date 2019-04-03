DataDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Data/"
ResultDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/"

library(limma)
library(edgeR)
library(qvalue)

run_limma <- function(tissue){

    data <- get.data(tissue, DataDir)
    
    count <- data[[3]]
    logCPM <- data[[2]]
    design <- data[[1]]
    design.nosvs <- data[[4]]

    ## ANALYZE WITH COVARIATES AND SVS ##
    ##fit <- lmFit(logCPM, design)
    ##fit <- eBayes(fit, trend=TRUE)
    ##limma.rslt <- topTable(fit,coef=2, n=Inf, sort.by="none")
    ##limma.rslt <- data.frame(gene = count$genes[,1], count$genes[,-1], limma.rslt)

    v <- voom(count, design, plot=F)
    vfit <- lmFit(v, design)    
    efit <- eBayes(vfit)
    voom.rslt <- topTable(efit, coef=2, n=Inf, sort.by="none")
    voom.rslt <- data.frame(gene = count$genes[,1], voom.rslt)

    out.file <- paste0(ResultDir, tissue,"_DE_limma.rslts")
    write.table(file = out.file, voom.rslt[order(voom.rslt$P.Value),], quote=F,
                row.names=F, col.names=T)   
    print(paste0("Wrote ", out.file))

    ## ANALYZE WITH COVARIANTES ONLY (No SVS)##
    v.nosvs <- voom(count, design.nosvs, plot=F)
    vfit <- lmFit(v.nosvs, design.nosvs)    
    efit <- eBayes(vfit)
    voom.rslt <- topTable(efit, coef=2, n=Inf, sort.by="none")
    voom.rslt <- data.frame(gene = count$genes[,1], voom.rslt)

    out.file <- paste0(ResultDir, tissue,"_DE_limma_nosvs.rslts")
    write.table(file = out.file, voom.rslt[order(voom.rslt$P.Value),], quote=F,
                row.names=F, col.names=T)   
    print(paste0("Wrote ", out.file))
   
    
    ## ANALYSIZE AFTER INVERSE NORMALIZATION ##
    cpm.invnorm <- inverse_quantile_normalization(v$E)
    txt <-  paste0("summary(lm(t(cpm.invnorm) ~ SEX + SMTSISCH + SMRIN + AGE + ", paste(
        names(count$sample)[grep("SV", names(count$sample))], collapse="+"),
        " ,data=count$samples))")
    out <- eval(parse(text = txt))
    out2 <- lapply(out, function(x) coefficients(x)[2,])
    lm.rslt <- as.data.frame(do.call(rbind, out2))
    names(lm.rslt) <- c("Est", "SE", "t.val","P.Value")
    lm.rslt <- data.frame(v$genes, lm.rslt)
    names(lm.rslt)[1] <- "ENSMBL"
    lm.rslt$adj.P.Val <- qvalue(lm.rslt$P.Value)$qvalue

    out.file <- paste0(ResultDir, tissue,"_LM_invnorm.rslts")
    write.table(file = out.file, lm.rslt[order(lm.rslt$P.Value),], quote=F,
                row.names=F, col.names=T)
    print(paste0("Wrote ", out.file))

     ## ANALYSIZE AFTER INVERSE NORMALIZATION WITH COVARIATES ONLY (NO SVS) ##
    cpm.invnorm <- inverse_quantile_normalization(v.no.svs$E)
    out <-  summary(lm(t(cpm.invnorm) ~ SEX + SMTSISCH + SMRIN + AGE,
                       ,data=count$samples))
    out2 <- lapply(out, function(x) coefficients(x)[2,])
    lm.rslt <- as.data.frame(do.call(rbind, out2))
    names(lm.rslt) <- c("Est", "SE", "t.val","P.Value")
    lm.rslt <- data.frame(v$genes, lm.rslt)
    names(lm.rslt)[1] <- "ENSMBL"
    lm.rslt$adj.P.Val <- qvalue(lm.rslt$P.Value)$qvalue

    out.file <- paste0(ResultDir, tissue,"_LM_invnorm_nosvs.rslts")
    write.table(file = out.file, lm.rslt[order(lm.rslt$P.Value),], quote=F,
                row.names=F, col.names=T)
    print(paste0("Wrote ", out.file))    
    
}

get.data <- function(tissue, DataDir){

    file <- paste0(DataDir, tissue, ".RDS")
    data <- readRDS(file)

    txt <- paste0("model.matrix(~SEX + SMTSISCH + SMRIN + AGE +",
                  paste(names(data$samples)[grep("SV", names(data$samples))], collapse="+"),
                  ", data$samples)")
    design <- eval(parse(text = txt))
    
    logCPM <- cpm(data, log=TRUE, prior.count=3)

    design.nosvs <- model.matrix(~SEX + SMTSISCH + SMRIN + AGE, data$samples)
    
    return(list(design = design, logCPM = logCPM, count = data, design.nosvs = design.nosvs))
}

inverse_quantile_normalization <- function(gct) {
        gct = t(apply(gct, 1, rank, ties.method = "average"));
        gct = qnorm(gct / (ncol(gct)+1));
        return(gct)
}


args <- commandArgs(TRUE)
tissue = args[1]
run_limma(tissue)
print(date())

finish.file = paste0(tissue,".finished")
system(paste0("touch ",finish.file))
