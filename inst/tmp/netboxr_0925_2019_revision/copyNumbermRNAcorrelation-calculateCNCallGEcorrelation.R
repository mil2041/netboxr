#' calculate correlation between copy number and gene expression
#' 
#' @param log2CNMat absolute copy number value in log2 scale 
#' @param CNCallMat copy number call data frame
#' @param log2GEMat gene expression mat in log2 scale
#' @param ampGeneList gene list of amp copy number
#' @param delGeneList gene list of deleted copy number
#' @param method method to calculate correlation, 'MannWhitney' or 'Pearson'
#' @param average average method 'mean' or 'median'
#' @param rankTestSize min sample size to calcualte rank test 
#' @param workDir where to store data
#' @param verbose output detailed information
#' 
#' @return list of correlation test results
#'   
#' @examples
#' cancerType<-'KIRC'
#' date<-'2015_08_21' 
#' #mutSig2CVMat<-getMutSig2CVMat(date=date,cancerType=cancerType,workDir=tmpDir)
#' #dim(mutSig2CVMat)
#' sampleSize<-333
#' threshold<-0.05
#' #worDir<-getwd()
#' 
#' #corrTable<-calculateCNGEcorr(log2CNMat,CNCallMat,log2GEMat,ampGeneList,
#' #delGeneList,method='MannWhitney',average=corrAverage,
#' #rankTestSize=minRankTestSize,workDir=filePath)
# 
#' 
#' @concept netboxr
#' @export
calculateCNGEcorr <- function(log2CNMat, CNCallMat, log2GEMat, ampGeneList, delGeneList, method = "MannWhitney", 
    average, rankTestSize, workDir, verbose) {
    
    buildcorrelationTest <- function(log2CN, log2GE) {
        PearsonPValue <- {
        }
        PearsonCorCoeff <- {
        }
        SpearmanPValue <- {
        }
        SpearmanCorCoeff <- {
        }
        
        sampleSize <- ncol(log2CN)
        probeSize <- nrow(log2CN)
        
        for (i in 1:probeSize) {
            cat("gene: ", rownames(log2CN)[i], " ", i, "/", probeSize, "\n")
            naNumCN <- sum(is.na(log2CN[i, ]))
            naNumGE <- sum(is.na(log2GE[i, ]))
            cutOffpercentage <- 0.7
            threshold <- (sampleSize) * cutOffpercentage
            if (naNumCN > threshold || naNumGE > threshold) {
                PearsonPValue[i] <- NA
                SpearmanPValue[i] <- NA
            } else {
                PearsonTest <- cor.test(t(log2CN[i, ]), t(log2GE[i, ]), method = "pearson", na.action = na.exclude)
                PearsonPValue[i] <- PearsonTest$p.value
                PearsonCorCoeff[i] <- PearsonTest$estimate
                
                SpearmanTest <- cor.test(t(log2CN[i, ]), t(log2GE[i, ]), method = "kendall", exact = FALSE, 
                  na.action = na.exclude)
                SpearmanPValue[i] <- SpearmanTest$p.value
                SpearmanCorCoeff[i] <- SpearmanTest$estimate
                
                png(file = paste(rownames(log2CN)[i], ".png", sep = ""), width = 800, height = 800)
                plot(as.matrix(log2CN[i, ]), as.matrix(log2GE[i, ]), main = paste(rownames(log2CN)[i]), 
                  xlab = "log2(Copy Number value)", ylab = "log2(RNA-seq-count-RSEM)", pch = 16, col = "blue", 
                  cex = 0.9)
                abline(lm(t(log2GE[i, ]) ~ t(log2CN[i, ]), na.action = na.exclude), col = "red", lty = 2)
                
                usr <- par("usr")
                text(usr[1] + 0.1, usr[4] - 0.2, adj = c(0, 1), paste("Pearson Coeff:", round(PearsonTest$estimate, 
                  digit = 4), ",P-Value:", round(PearsonTest$p.value, digit = 4), sep = " "))
                text(usr[1] + 0.1, usr[4] - 0.5, adj = c(0, 1), paste("Spearman Coeff:", round(SpearmanTest$estimate, 
                  digit = 4), ",P-Value:", round(SpearmanTest$p.value, digit = 4), sep = " "))
                dev.off()
                
            }
            
        }
        
        data <- data.frame(rownames(log2CN), PearsonCorCoeff, PearsonPValue, SpearmanCorCoeff, SpearmanPValue)
        colnames(data) <- c("geneSymbol", "PearsonCorCoeff", "PearsonPValue", "SpearmanCorCoeff", "SpearmanPValue")
        
        return(data)
    }
    
    buildMannWhitneyTest <- function(CNCall, log2GE, copyNumberType, average = "median", rankTestSize = 5) {
        MannWhitneyPValue <- {
        }
        log2TestMean <- {
        }
        log2DiploidMean <- {
        }
        log2Fold <- {
        }
        
        # minimum sample size in each group rankTestSize<-5
        sampleSize <- ncol(CNCall)
        probeSize <- nrow(CNCall)
        
        for (i in 1:probeSize) {
            
            pickedGeneSymbol <- rownames(CNCall)[i]
            
            cat("gene: ", pickedGeneSymbol, " ", i, "/", probeSize, "\n")
            naNumCN <- sum(is.na(CNCall[i, ]))
            naNumGE <- sum(is.na(log2GE[i, ]))
            cutOffpercentage <- 0.7
            threshold <- (sampleSize) * cutOffpercentage
            
            # if 70% samples have NA values, skip test
            if (naNumCN > threshold || naNumGE > threshold) {
                MannWhitneyPValue[i] <- NA
                log2TestMean[i] <- NA
                log2DiploidMean[i] <- NA
                log2Fold[i] <- NA
                
            } else {
                
                
                dat <- data.frame(t(CNCall[i, ]), t(log2GE[i, ]))
                dat <- na.exclude(dat)
                colnames(dat) <- c("conds", "geneCounts")
                
                levels <- c("Homdel", "Hetloss", "Diploid", "Gain", "Amp")
                
                # assign copy number category name to GISTIC2 copy number call (-2,2)
                for (k in 1:length(levels)) {
                  j <- (k - 3)
                  if (sum(dat$conds == j) > 0) {
                    dat[dat$conds == j, ]$conds <- levels[k]
                  }
                }
                
                if (copyNumberType %in% "Amp") {
                  ampCN <- dat[which(dat$conds %in% "Amp"), ]$geneCounts
                  diploidCN <- dat[which(dat$conds %in% "Diploid"), ]$geneCounts
                  
                  if (length(ampCN) >= rankTestSize & length(diploidCN) >= rankTestSize) {
                    
                    
                    MannWhitneyTest <- wilcox.test(ampCN, diploidCN, na.action = na.exclude, exact = FALSE, 
                      alternative = "greater")
                    MannWhitneyPValue[i] <- MannWhitneyTest$p.value
                    
                    if (average %in% "median") {
                      log2TestMean[i] <- log(median(2^(ampCN)), base = 2)
                      log2DiploidMean[i] <- log(median(2^(diploidCN)), base = 2)
                    }
                    
                    if (average %in% "mean") {
                      log2TestMean[i] <- log(mean(2^(ampCN)), base = 2)
                      log2DiploidMean[i] <- log(mean(2^(diploidCN)), base = 2)
                    }
                    
                    log2Fold[i] <- log((2^(log2TestMean[i])/2^(log2DiploidMean[i])), base = 2)
                    # cat(sprintf('pvalue %s',MannWhitneyTest$p.value))
                    
                    # png(file=paste(pickedGeneSymbol,'_CNCall_GE_corr_normalized_edgeR.png',sep=''),width=1024,height=676)
                    # plotCNCallGE(dat=dat,levels=levels,geneSymbol=pickedGeneSymbol)
                    
                    # usr<-par('usr') text(usr[1]+0.1,usr[4]-0.4,adj=c(0,1),paste('Two samples Mann-Whitney test',sep=' '))
                    # text(usr[1]+0.1,usr[4]-0.8,adj=c(0,1),paste('One-tail
                    # pValue:',round(MannWhitneyTest$p.value,digit=5),sep=' '))
                    
                    # dev.off()
                    
                    
                  } else {
                    
                    if (length(ampCN) < rankTestSize) {
                      
                      cat(sprintf("Amp:%s Diploid:%s\n", length(ampCN), length(diploidCN)))
                      MannWhitneyPValue[i] <- NA
                      log2TestMean[i] <- NA
                      log2DiploidMean[i] <- NA
                      log2Fold[i] <- NA
                      
                    } else {
                      
                      # special cases when there is no diploid samples
                      cat(sprintf("Amp:%s Diploid:%s\n", length(ampCN), length(diploidCN)))
                      MannWhitneyPValue[i] <- 10^(-5)
                      log2TestMean[i] <- log(median(2^(ampCN)), base = 2)
                      log2DiploidMean[i] <- 2
                      log2Fold[i] <- 5
                      
                    }
                    
                    
                    
                  }
                  
                  png(file = paste(pickedGeneSymbol, "_CNCall_GE_corr_normalized_edgeR.png", sep = ""), 
                    width = 1024, height = 676)
                  plotCNCallGE(dat = dat, levels = levels, geneSymbol = pickedGeneSymbol)
                  
                  usr <- par("usr")
                  text(usr[1] + 0.1, usr[4] - 0.4, adj = c(0, 1), paste("Two samples Mann-Whitney test", 
                    sep = " "))
                  text(usr[1] + 0.1, usr[4] - 0.8, adj = c(0, 1), paste("One-tail pValue:", format(MannWhitneyPValue[i], 
                    digit = 5), sep = " "))
                  
                  dev.off()
                  
                  
                  
                }
                
                if (copyNumberType %in% "Homdel") {
                  homdelCN <- dat[which(dat$conds %in% "Homdel"), ]$geneCounts
                  diploidCN <- dat[which(dat$conds %in% "Diploid"), ]$geneCounts
                  
                  if (length(homdelCN) >= rankTestSize) {
                    
                    MannWhitneyTest <- wilcox.test(homdelCN, diploidCN, na.action = na.exclude, exact = FALSE, 
                      alternative = "less")
                    MannWhitneyPValue[i] <- MannWhitneyTest$p.value
                    
                    if (average %in% "median") {
                      log2TestMean[i] <- log(median(2^(homdelCN)), base = 2)
                      log2DiploidMean[i] <- log(median(2^(diploidCN)), base = 2)
                    }
                    
                    if (average %in% "mean") {
                      log2TestMean[i] <- log(mean(2^(homdelCN)), base = 2)
                      log2DiploidMean[i] <- log(mean(2^(diploidCN)), base = 2)
                    }
                    
                    log2Fold[i] <- log((2^(log2TestMean[i])/2^(log2DiploidMean[i])), base = 2)
                    
                    # png(file=paste(pickedGeneSymbol,'_CNCall_GE_corr_normalized_edgeR.png',sep=''),width=1024,height=676)
                    # plotCNCallGE(dat=dat,levels=levels,geneSymbol=pickedGeneSymbol)
                    
                    # usr<-par('usr') text(usr[1]+0.1,usr[4]-0.4,adj=c(0,1),paste('Two samples Mann-Whitney test',sep=' '))
                    # text(usr[1]+0.1,usr[4]-0.8,adj=c(0,1),paste('One-tail
                    # pValue:',round(MannWhitneyTest$p.value,digit=5),sep=' '))
                    
                    # dev.off()
                    
                    
                    # cat(sprintf('catch %s',log2TestMean[i]))
                  } else {
                    cat(sprintf("Homdel:%s Diploid:%s\n", length(homdelCN), length(diploidCN)))
                    MannWhitneyPValue[i] <- NA
                    log2TestMean[i] <- NA
                    log2DiploidMean[i] <- NA
                    log2Fold[i] <- NA
                    
                  }
                  
                  png(file = paste(pickedGeneSymbol, "_CNCall_GE_corr_normalized_edgeR.png", sep = ""), 
                    width = 1024, height = 676)
                  plotCNCallGE(dat = dat, levels = levels, geneSymbol = pickedGeneSymbol)
                  
                  usr <- par("usr")
                  text(usr[1] + 0.1, usr[4] - 0.4, adj = c(0, 1), paste("Two samples Mann-Whitney test", 
                    sep = " "))
                  text(usr[1] + 0.1, usr[4] - 0.8, adj = c(0, 1), paste("One-tail pValue:", format(MannWhitneyPValue[i], 
                    digit = 5), sep = " "))
                  
                  dev.off()
                  
                  
                }
                
                
                
            }
            
        }
        
        data <- data.frame(rownames(CNCall), log2TestMean, log2DiploidMean, log2Fold, MannWhitneyPValue, 
            stringsAsFactors = FALSE)
        colnames(data) <- c("geneSymbol", "log2TestMean", "log2DiploidMean", "log2Fold", "MannWhitneyPValue")
        data$copyNumberType <- rep(copyNumberType, nrow(CNCall))
        
        return(data)
    }
    
    
    currentDir <- getwd()
    
    if (!file.exists(file.path(workDir, "CN_GE_corr_MannWhitney", "amp"))) {
        dir.create(file.path(workDir, "CN_GE_corr_MannWhitney", "amp"), recursive = TRUE)
    }
    
    filePath <- file.path(workDir, "CN_GE_corr_MannWhitney", "amp")
    setwd(filePath)
    
    CNCallAmp <- CNCallMat[which(rownames(CNCallMat) %in% ampGeneList), ]
    log2GEAmp <- log2GEMat[which(rownames(log2GEMat) %in% ampGeneList), ]
    ampTest <- buildMannWhitneyTest(CNCall = CNCallAmp, log2GE = log2GEAmp, copyNumberType = "Amp", average, 
        rankTestSize)
    # ampTest<-data.frame(ampTest$geneSymbol,ampTest$MannWhitneyPValue)
    ampTest <- ampTest[complete.cases(ampTest), ]
    # colnames(ampTest)<-c('geneSymbol','MannWhitneyPValue')
    # ampTest$padj<-p.adjust(ampTest$MannWhitneyPValue,method='BH') ampTest$type<-rep('Amp',nrow(ampTest))
    # ampTest<-ampTest[order(ampTest$padj),]
    
    
    if (!file.exists(file.path(workDir, "CN_GE_corr_MannWhitney", "del"))) {
        dir.create(file.path(workDir, "CN_GE_corr_MannWhitney", "del"), recursive = TRUE)
    }
    
    filePath <- file.path(workDir, "CN_GE_corr_MannWhitney", "del")
    setwd(filePath)
    
    CNCallDel <- CNCallMat[which(rownames(CNCallMat) %in% delGeneList), ]
    log2GEDel <- log2GEMat[which(rownames(log2GEMat) %in% delGeneList), ]
    delTest <- buildMannWhitneyTest(CNCall = CNCallDel, log2GE = log2GEDel, copyNumberType = "Homdel", 
        average, rankTestSize)
    # delTest<-data.frame(delTest$geneSymbol,delTest$MannWhitneyPValue)
    delTest <- delTest[complete.cases(delTest), ]
    # colnames(delTest)<-c('geneSymbol','MannWhitneyPValue')
    # delTest$padj<-p.adjust(delTest$MannWhitneyPValue,method='BH')
    # delTest$type<-rep('Homdel',nrow(delTest)) delTest<-delTest[order(delTest$padj),]
    
    setwd(currentDir)
    
    corrTable <- list(ampTest = ampTest, delTest = delTest)
    return(corrTable)
    
    
    # setwd('~/work/Sander_lab/TCGA_data/BRCA/TCGA/CN_GE_correlation')
    # kk<-buildcorrelationTest(log2CN=MergedData$log2CNMat,log2GE=MergedData$log2GEMat)
    
    
    
    if (FALSE) {
        
        colnames(kk) <- c("geneSymbol", "PearsonCorCoeff", "PearsonPValue", "SpearmanCorCoeff", "SpearmanPValue")
        
        write.table(kk, file = "wholeCorrelationTable.txt", sep = "\t", quote = FALSE, row.names = FALSE)
        save(kk, file = "log2CN_log2GE_pValue.RData")
        
        k3 <- data.frame(kk$geneSymbol, kk$PearsonPValue)
        k4 <- k3[complete.cases(k3), ]
        k5 <- p.adjust(k4$kk.PearsonPValue, method = "BH")
        k6 <- data.frame(k4, k5)
        colnames(k6) <- c("GeneSymbol", "PearsonPValue", "PearsonFDR")
        
        alpha <- 0.05
        k7 <- k6[k6$PearsonFDR > alpha, ]
        k8 <- k7[order(k7$PearsonFDR), ]
        write.table(k8, file = "notCorrelatedPearson.txt", sep = "\t", quote = FALSE, row.names = FALSE)
        
        log2CNGEnotCorrelatedPearson <- k8
        
        # =============
        
        k3 <- data.frame(kk$geneSymbol, kk$SpearmanPValue)
        k4 <- k3[complete.cases(k3), ]
        k5 <- p.adjust(k4$kk.SpearmanPValue, method = "BH")
        k6 <- data.frame(k4, k5)
        colnames(k6) <- c("GeneSymbol", "SpearmanPValue", "SpearmanFDR")
        
        alpha <- 0.05
        k7 <- k6[k6$SpearmanFDR > alpha, ]
        k8 <- k7[order(k7$SpearmanFDR), ]
        write.table(k8, file = "notCorrelatedSpearman.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    
}






