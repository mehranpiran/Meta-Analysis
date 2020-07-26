setwd("../../GSE9348_GSE10961/Data")

# Importing data into R
library(affy)
affybatch = ReadAffy()
affybatch


# Probe summarization, Backgraound correction
eset = rma(affybatch , normalize = F)
eset
expression.matrix = exprs(eset)
dim(expression.matrix)
pData(phenoData(eset))




# Removing rows containing more than 50% absent probsets in each group
Absent.probes = mas5calls.AffyBatch(affybatch)
Absent.probes
Absent.probes = exprs(Absent.probes)#large matrix of A/M/P probesets
dim(Absent.probes)

index.absent.probs = c()

for(i in 1:length(Absent.probes[,1])){
  if(sum(Absent.probes[i,1:18] == "A") > 9 | sum(Absent.probes[i,19:39] == "A") > 10 | sum(Absent.probes[i,40:51] == "A") > 6){
    index.absent.probs = c(index.absent.probs,i)
  } 
}

expression.matrix = expression.matrix[-index.absent.probs,]
pData(featureData(eset)) = pData(featureData(eset))[-index.absent.probs,]
dim(expression.matrix)




# Altering sample names
samples = paste0(substr(colnames(expression.matrix) , 1 , 10) , c(rep("_M" , 18) , rep("_P" , 21) , rep("_N" , 12)))
colnames(expression.matrix) = samples
dim(expression.matrix)




# PCA plot to recognize biased samples based on eigenvector 1 and eigenvector 2
library(ggplot2)
pc = prcomp(expression.matrix)
pcr = pc$rotation
pcr = as.data.frame(pcr)
pcr$sample = rownames(pcr)
pcr$groups=c(rep("_M" , 18) , rep("_P" , 21) , rep("_N" , 12))
pcr$groups=factor(pcr$groups , levels = c("_M" , "_N" , "_P"))

ggplot(pcr , aes(PC1 , PC2 , label = sample , colour = groups)) + geom_text( size = 4) + xlim(-0.153,-0.127) +
  geom_label(label.size = 0.25) + theme_linedraw() +
  theme(axis.title=element_text(size=18 , face  = "bold") , axis.text=element_text(size=14 , colour = "black") , 
        title = element_text(size=18) , legend.text = element_text(size = 14),legend.title = element_text(size = 18),
        legend.background = element_rect(color = "steelblue", linetype = "solid"), 
        legend.key = element_rect(fill = NULL, color = "black") )





# Recognizing and removing biased samples based on hierarchical clustering and Number-SD method ()
l = list()
index = list(1:18 , 19:39 , 40:51)
out.samples = list()

for(i in 1:3){
  
  mat = expression.matrix[ , index[[i]]]
  
  IAC=cor(mat,use="p")#1st round
  hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
  library(cluster)
  cluster1=hclust(as.dist(1-IAC),method="average")
  plot(cluster1,cex=0.7,labels=dimnames(mat)[[2]])
  meanIAC=apply(IAC,2,mean)
  sdCorr=sd(meanIAC)
  numbersd=(meanIAC-mean(meanIAC))/sdCorr
  plot(numbersd)
  abline(h=-2)
  sdout=-2
  outliers=dimnames(mat)[[2]][numbersd<sdout]
  outliers
  mat = mat[,numbersd>sdout]
  dim(mat)
  
  sample.outliers = outliers
  while(length(outliers) != 0){
    IAC=cor(mat,use="p")
    hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
    cluster1=hclust(as.dist(1-IAC),method="average")
    plot(cluster1,cex=0.7,labels=dimnames(mat)[[2]])
    meanIAC=apply(IAC,2,mean)
    sdCorr=sd(meanIAC)
    numbersd=(meanIAC-mean(meanIAC))/sdCorr
    plot(numbersd)
    abline(h=-2)
    sdout=-2
    outliers=dimnames(mat)[[2]][numbersd<sdout]
    sample.outliers = c(sample.outliers,outliers)
    mat=mat[,numbersd>sdout]
    dim(mat) 
    
  }
  
  out.samples[[i]] = sample.outliers
  
  l[[i]] = mat
  
}

expression.matrix = do.call("cbind" , l)
dim(expression.matrix)
out.samples = unlist(out.samples)
out.samples

write.table(out.samples , "out.samples.txt" , sep = "\t" , col.names = F , row.names = F , quote = F)


unlockBinding("exprs",assayData(eset))
assayData(eset)$exprs = expression.matrix
eset



###################
## Normalization ##
###################

# Normalization using quantile method
library(preprocessCore)
normalized.expression.matrix = normalize.quantiles(expression.matrix)
dimnames(normalized.expression.matrix) = dimnames(expression.matrix)


unlockBinding("exprs",assayData(eset))
assayData(eset)$exprs = normalized.expression.matrix
eset


boxplot(expression.matrix , pch = "." , las=3)
boxplot(normalized.expression.matrix , pch = "." , las=3)

##########################################
## Filtering and many to many problems ###
##########################################

# Removing low variant genes 
# Select one probeset with the largest  IQR to be representative of other probesets mapped to the same gene symbole

sds = apply(normalized.expression.matrix , 1 , sd)
hist(sds, breaks=100, col="mistyrose", xlab="standard deviation" )
abline(v=quantile(sds)[3], col="blue", lwd=3, lty=2)

annotation(eset)
# "hgu133plus2"  

library(genefilter)
library(hgu133plus2.db)
feset = nsFilter(eset, remove.dupEntrez=T, var.cutof = quantile(sds)[3] ) 
filtered.eset = feset$eset
filtered.eset

N.F.expression.matrix = exprs(filtered.eset)
dim(N.F.expression.matrix)


############################################
##### Differentially Expressed Genes #######
############################################

colnames(expression.matrix)

library(limma)

matrix.expression = N.F.expression.matrix

colnames(matrix.expression) = c(rep("M" , length(grep("_M" ,colnames(expression.matrix)))) , rep("P" , length(grep("_P" ,colnames(expression.matrix)))) , rep("N" , length(grep("_N" ,colnames(expression.matrix)))))

colnames(matrix.expression)

gr = factor(colnames(matrix.expression) , levels = c("M" , "P" , "N"))

design = model.matrix(~0 + gr)

colnames(design) = c("M", "P", "N")

lm.fit = lmFit(matrix.expression, design)


### Annotation ###
library(annotate)
probnames = rownames(matrix.expression)
gene.symbol = getSYMBOL(probnames, "hgu133plus2")

mc = makeContrasts(M-P , levels = design)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)
Table = topTable(eb, adjust.method = "BH", sort.by = "logFC", genelist = gene.symbol , number = Inf)
write.table(Table , "MvsP_DEGs.txt" , sep = "\t")



mc = makeContrasts(M-N , levels = design)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)
Table = topTable(eb, adjust.method = "BH", sort.by = "logFC", genelist = gene.symbol , number = Inf)
write.table(Table , "MvsN_DEGs.txt" , sep = "\t")


mc = makeContrasts(P-N , levels = design)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)
Table = topTable(eb, adjust.method = "BH", sort.by = "logFC", genelist = gene.symbol , number = Inf)
write.table(Table , "PvsN_DEGs.txt" , sep = "\t")



#####################
### Quality Plot ####
#####################

annotation(eset)
library(genefilter)
library(hgu133plus2.db)
library(annotate)

feset = nsFilter(eset, remove.dupEntrez=T, var.cutof = 10^-20) 
filtered.eset = feset$eset
filtered.eset

N.F.expression.matrix = exprs(filtered.eset)
dim(N.F.expression.matrix)

p1 = rownames(Table[which(Table$logFC > 6) , ])
p1
p2 = rownames(Table[which(Table$logFC < - 4) , ])
p2

probnames = c(p2,p1)
gene.symbol = getSYMBOL(probnames, "hgu133plus2")
p = unname(gene.symbol)
p

probnames = rownames(N.F.expression.matrix)
gene.symbol = getSYMBOL(probnames, "hgu133plus2.db")
rownames(N.F.expression.matrix) = gene.symbol

HG = c("ACTB" , "GAPDH" , "TBP" , "RPLP0")


HG = HG[HG %in% gene.symbol]
HG

p = N.F.expression.matrix[c(HG,p),] # 4 HK genes, 2 down-regulated genes, 8 up-regulated genes


d = data.frame(apply(p[,1:18] , 1 , mean) , apply(p[,19:37] , 1 , mean)) # Between metastatic and primary

d = data.frame(d , Genes = c(rep("Houskeeping" , 4) , rep("Down-Regulated" , 2), rep("Up-Regulated" , 8)) , stringsAsFactors = F)

d[,3] = as.factor(d[,3])
Genes = d$Genes
d$size = 1

colnames(d) = c("One" , "Two" , "Genes" , "size")
d$NAME = rownames(d)


library(ggplot2)
library(ggrepel)

g = ggplot(d, aes(One , Two))

g + geom_point(aes(color = Genes , size = 0.1)) + labs(title = "") + ylab("Primary") + xlab("Metastatic") + ylim(2,14) + xlim(2,14) +
  theme(axis.title=element_text(size=15 , face  = "bold") , axis.text=element_text(size=16 , colour = "black") , plot.title = element_text(hjust = 0.5) ,
        title = element_text(size=16) , legend.text = element_text(size = 12),legend.title = element_text(size = 18)) + 
  scale_size_continuous(guide = F)  + guides(colour = guide_legend(override.aes = list(size=2, stroke=2))) +
  geom_label_repel(aes(label = NAME),
                   box.padding   = 0.5, 
                   point.padding = 0.5,
                   segment.color = 'grey50')

