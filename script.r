
#########################################################################
# libraries

library(ggplot2) 
library(ggfortify)

library(DESeq2)

#########################################################################
# customized functions for the dataset

plotmyPCA <- function(dat, TITLE, FileName){

	pca  = prcomp(t(dat))
	pVar = pca$sdev^2/sum(pca$sdev^2) 

	colorPlate  = rainbow(4)
	colorScheme = c(rep(colorPlate[1],4), rep(colorPlate[2],4), rep(colorPlate[3],4), rep(colorPlate[4],4))
	
	autoplot(pca, colour=colorScheme, label=TRUE, label.size=5, size=2.5, main=TITLE,xlab=paste("PC1: ",round(pVar[1]*100),"% variance", sep=""), ylab=paste("PC2: ",round(pVar[2]*100), "% variance", sep=""))

	ggsave(filename=FileName, plot=last_plot(), width=8, height=5, units="in", dpi =600)
}

#########################################################################

dat = read.csv("raw_cts.txt", sep="\t", header=TRUE)
plotmyPCA(as.matrix(dat[,3:18]), "Unprocessed Counts", "PCA1.png")

#########################################################################
# nomenclature filtering
#########################################################################
# low-count filtering

# reformat dat into a matrix
CTS = as.matrix(dat[,3:18])
rownames(CTS) = dat[,1]

# first filter
index1 = which(rowSums(CTS) < 30)
 
cts_f1 =  CTS[-index1, ]
plotmyPCA(as.matrix(cts_f1), "Filter1", "PCA2.png")
#######################################
# second filter
# sampleAv < 10

numSampleAvLess10 = NULL
meanM = NULL
# mff	1:4
# mfm	5:8
# muf	9:12
# mum	13:16

for( i in 1:nrow(cts_f1)){
	mean1 = mean(cts_f1[i,1:4])
	mean2 = mean(cts_f1[i,5:8])
	mean3 = mean(cts_f1[i,9:12])
	mean4 = mean(cts_f1[i,13:16])

	meanM = rbind(meanM, c(mean1, mean2, mean3, mean4))
	tmp   = length(as.vector(which(meanM[i,]<10)))

	numSampleAvLess10 = c(numSampleAvLess10, tmp)
}

index  = as.vector(which(numSampleAvLess10==4))
cts_f2 = cts_f1[-index, ]

plotmyPCA(cts_f2, "Filter 2", "PCA3.png")

###########################################################################

cts = cts_f2

treatment    <- c(rep("mff",4), rep("mfm",4), rep("muf",4), rep("mum",4))

colData<- data.frame(colnames(cts), 
                     treatment = as.factor(treatment))

dds   <- DESeqDataSetFromMatrix(countData = cts, colData=colData, design = ~ treatment)

esf   <- estimateSizeFactors(dds)
sf    <- sizeFactors(esf) 
dds   <- DESeq(dds)
res   <- results(dds)
resultsNames(dds)
test <-estimateSizeFactorsForMatrix(cts)
sf == test ## identical!! test2 is computed by dividing by the geometric means of rows
            ## the median of all ratios is the size factor

normCounts = counts(esf,normalized=TRUE)	##the normalized counts

test1 = cts[,1]/sf[1]
sum(normCounts[,1] != test1) 			##should be 0
rm(sf,test1, cts_f1, cts_f2, CTS)

plotmyPCA(normCounts,  "Normalized Counts", "PCA4.png")

mff1 = normCounts[,1]
mff2 = normCounts[,2]
plot(mff1,mff2)

rld = rlog(dds)
vsd = varianceStabilizingTransformation(dds)

plotmyPCA(assay(rld), "rlog data", "PCA5.png")
plotmyPCA(assay(vsd), "vsd data", "PCA6.png")

ht = hclust(dist(t(normCounts)), 'ave')
plot(ht)

ht = hclust(dist(t(assay(rld))), 'ave')
plot(ht)

C1 <- results(dds, contrast = c("treatment", "mff", "muf"), addMLE = FALSE, independentFiltering = FALSE)
diffExp = data.frame(C1)
head(diffExp)
index = which(diffExp[,6] <= 0.05)
length(index)

##########################################################################
