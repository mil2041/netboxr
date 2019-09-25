#' boxplot of gene expression counts between multiple copy number groups
#' 
#' @param dat Default is 'latest' 
#' @param levels Abbreviation in the TCGA study ['KIRC','GBM']
#' @param geneSymbol MutSig2CV result data frame
#' @param verbose output detailed information
#' 
#' @return null
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
#' #plotCopyNumberFreq(dataWithPos=dataWithPos,
#' #          threshold=copyNumberFreqThreshold,chromSize=chromSize)
#' 
#' @concept netboxr
#' @export
plotCNCallGE <- function(dat, levels, geneSymbol, verbose = TRUE) {
    
    pickedGeneSymbol <- geneSymbol
    
    colnames(dat) <- c("conds", "geneCounts")
    
    
    
    # assign copy number category name to GISTIC2 copy number call (-2,2)
    for (i in 1:length(levels)) {
        j <- (i - 3)
        if (sum(dat$conds == j) > 0) {
            dat[dat$conds == j, ]$conds <- levels[i]
        }
    }
    
    # library(RColorBrewer) NumberOfLevels<-length(levels(dat$conds))
    NumberOfLevels <- 5
    # mycolors<-brewer.pal(n=NumberOfLevels, name='Set1')
    
    if (NumberOfLevels < 3) {
        mycolors <- c("red", "blue")
        
    } else {
        # library(RColorBrewer) mycolors<-brewer.pal(n=NumberOfLevels, name='Set1')
        mycolors <- c("royalblue", "skyblue", "white", "pink", "red")
    }
    
    # levels<-c('Homdel','Hetloss','Diploid','Gain','Amp')
    
    # assign copy number category name to GISTIC2 copy number call (-2,2)
    for (i in 1:length(levels)) {
        j <- (i - 3)
        if (sum(dat$conds == j) > 0) {
            dat[dat$conds == j, ]$conds <- levels[i]
        }
    }
    
    condsProportions <- table(dat$conds)/length(dat$conds)
    
    levelProportions <- condsProportions[levels]
    for (i in 1:length(levels)) {
        if (is.na(levelProportions[i])) {
            levelProportions[i] <- 0
        }
    }
    names(levelProportions) <- levels
    
    dat$conds <- factor(dat$conds, levels = levels)
    
    # png(file=paste(pickedGeneSymbol,'_normalized_edgeR.png',sep=''),width=750,height=520)
    
    boxplot(dat$geneCounts ~ dat$conds, col = mycolors, outpch = NA, na.action = na.exclude, varwidth = TRUE)
    # boxplot(dat$geneCounts ~ dat$conds, col=mycolors,outpch=NA) stripchart(dat$geneCounts ~ dat$conds,
    # vertical=TRUE, data=dat, method='jitter',add=TRUE,pch=16,col=mycolors,cex=0.5)
    
    title(paste("Boxplot for gene: ", pickedGeneSymbol, sep = ""), xlab = "Sub-groups", ylab = "log2 (normalized counts)")
    
    # Add data points
    mylevels <- levels(dat$conds)
    for (i in 1:length(mylevels)) {
        thislevel <- mylevels[i]
        thisvalues <- dat[dat$conds == thislevel, 2]
        
        if (!length(thisvalues) == 0) {
            
            # take the x-axis indices and add a jitter, proportional to the N in each level
            myjitter <- jitter(rep(i, length(thisvalues)), amount = levelProportions[i]/2)
            points(myjitter, thisvalues, pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
            
            # While we're looping, lets add some text I played with this a lot for this particular plot.  I tried
            # to make it general, but you may have to adjust a bit for different data.
            TopOfWhisker <- min(max(thisvalues), median(thisvalues) + IQR(thisvalues) * 3)
            text(i + levelProportions[i]/2, TopOfWhisker, labels = paste("N=", length(thisvalues), sep = ""), 
                cex = 0.95, font = 10, pos = 4)
            
        }
    }
    
}
