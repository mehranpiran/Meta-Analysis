setwd("F:/Systems Biomedicine of Cancer/GSE18105")
load("GSE18105.RData")


library(ggplot2)
library(affy)

setwd("F:/Project3/GSE18105/Data")


affybatch = ReadAffy()
affybatch
affybatch = affybatch[,idx_study]
#BiocManager::install("simpleaffy")
#library(simpleaffy)

#qc = qc(affybatch)

#plot(qc) 



eset = rma(affybatch , normalize = F)
eset

expression.matrix = exprs(eset)
dim(expression.matrix)

#a = pData(phenoData(eset))

Absent.probes = mas5calls.AffyBatch(affybatch)#returns 
Absent.probes
Absent.probes = exprs(Absent.probes)#large matrix of A/P probes
dim(Absent.probes)
#A/p/M removal
index.absent.probs = c()

for(i in 1:length(Absent.probes[,1])){
  if(sum(Absent.probes[i,1:20] == "A") > 2 | sum(Absent.probes[i,21:30] == "A") > 2 | sum(Absent.probes[i,31:47] == "A") > 2){
    index.absent.probs = c(index.absent.probs,i)
  } 
}

expression.matrix = expression.matrix[-index.absent.probs,]
dim(expression.matrix)

samples = paste0(substr(colnames(expression.matrix) , 1 , 9) , c(rep("_Metastatic" , 20) , rep("_Primary" , 10) , rep("_Normal" , 17)))
colnames(expression.matrix) = samples
dim(expression.matrix)
library(ggplot2)
pc = prcomp(expression.matrix)
pcr = pc$rotation
pcr = as.data.frame(pcr)
pcr$sample = rownames(pcr)
pcr$groups=c(rep("_Metastatic" , 20) , rep("_Primary" , 10) , rep("_Normal" , 17))

ggplot(pcr , aes(PC1 , PC2 , label = sample , colour = groups)) + geom_text(size = 3) + 
  theme(axis.title=element_text(size=15 , face  = "bold") , axis.text=element_text(size=16 , colour = "black") ,
        title = element_text(size=16) , legend.text = element_text(size = 12),legend.title = element_text(size = 18))




dim(expression.matrix)
15209    47

mat = expression.matrix[,31:47]

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
sample.outliers

a1 = mat
a2 = mat
a3 = mat

mat = cbind(a1,a2,a3)
dim(mat)







expression.matrix = mat

pData(featureData(eset)) = pData(featureData(eset))[-index.absent.probs,]
unlockBinding("exprs",assayData(eset))
assayData(eset)$exprs = expression.matrix
eset


colnames(expression.matrix)
#1:15   16:22  23:38


pc = prcomp(expression.matrix)
pcr = pc$rotation
pcr = as.data.frame(pcr)
pcr$sample = rownames(pcr)
pcr$groups=c(rep("_Metastatic" , 15) , rep("_Primary" , 7) , rep("_Normal" , 16))

ggplot(pcr , aes(PC1 , PC2 , label = sample , colour = groups)) + geom_text(size = 3) + 
  theme(axis.title=element_text(size=15 , face  = "bold") , axis.text=element_text(size=16 , colour = "black") ,
        title = element_text(size=16) , legend.text = element_text(size = 12),legend.title = element_text(size = 18))



ggplot(pcr , aes(PC1 , PC2 , label = sample , colour = groups)) + geom_text( size = 4) + xlim(-0.16, -0.121) +
  geom_label(label.size = 0.25) + theme_bw() +
  theme(axis.title=element_text(size=16 , face  = "bold") , axis.text=element_text(size=16 , colour = "black") , 
        title = element_text(size=16) , legend.text = element_text(size = 12),legend.title = element_text(size = 18),
        legend.background = element_rect(color = "steelblue", linetype = "solid"), 
        legend.key = element_rect(fill = NULL, color = "black") , legend.key.size = unit(0.2, "lines"))


###################
## Normalization ##
###################

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
#

############################################
##### Differential Expression Genes ########
############################################


library(limma)

matrix.expression = N.F.expression.matrix

colnames(matrix.expression) = c(rep("M" , 15) , rep("P" , 7) , rep("N" , 16))

colnames(matrix.expression)

gr = factor(colnames(matrix.expression) , levels = c("M" , "P" , "N"))

design = model.matrix(~0 + gr)

colnames(design) = c("M", "P", "N")

lm.fit = lmFit(matrix.expression, design)

mc = makeContrasts(M-P , levels = design)

c.fit = contrasts.fit(lm.fit, mc)

eb = eBayes(c.fit)



### Annotation ##########

library(annotate)
probnames = rownames(matrix.expression)
any(duplicated(probnames ))

gene.ID = getEG(probnames, "hgu133plus2")
gene.symbol = getSYMBOL(probnames, "hgu133plus2")

Table = topTable(eb, adjust.method = "BH", sort.by = "p", genelist = gene.symbol , number = 3000)
min(Table$logFC) # -6.811721
max(Table$logFC) # 5.168648


DEGsMvsN = Table[Table$adj.P.Val < 0.05 & abs(Table$logFC) > 0.5 , ]
DEGsMvsP = Table[Table$adj.P.Val < 0.05 & abs(Table$logFC) > 0.5 , ]
DEGsPvsN = Table[Table$adj.P.Val < 0.05 & abs(Table$logFC) > 0.5 , ]

setwd("F:/Project3/GSE18105")

write.table(DEGsMvsN, file = "DEGsMvsN.txt" , sep = "\t" , quote = F)
write.table(DEGsMvsP, file = "DEGsMvsP.txt" , sep = "\t" , quote = F)
write.table(DEGsPvsN, file = "DEGsPvsN.txt" , sep = "\t" , quote = F)

save(DEGsCvsN , file = "DEGsCvsN.RData")

boxplot(N.F.expression.matrix , pch = "." , las=3)











### Quality Plot ####

annotation(eset)
library(genefilter)
library(hgu133plus2.db)
library(annotate)

feset = nsFilter(eset, remove.dupEntrez=T, var.cutof = 10^-20) 
filtered.eset = feset$eset
filtered.eset

N.F.expression.matrix = exprs(filtered.eset)
dim(N.F.expression.matrix)

p1 = rownames(Table[which(Table$logFC > 2) , ])
p1
p2 = rownames(Table[which(Table$logFC < - 2) , ])
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

p = N.F.expression.matrix[c(HG,p),] # 4 HK genes, 12 down-regulated genes, 6 up-regulated genes


d = data.frame(apply(p[,1:15] , 1 , mean) , apply(p[,16:22] , 1 , mean)) # Between metastatic and primary

d = data.frame(d , Genes = c(rep("Houskeeping" , 4) , rep("Down-Regualted" , 12), rep("Up-Regulated" , 6)) , stringsAsFactors = F)

d[,3] = as.factor(d[,3])
Genes = d$Genes
d$size = 1

colnames(d) = c("One" , "Two" , "Genes" , "size")
d$NAME = rownames(d)


library(ggplot2)
library(ggrepel)

g = ggplot(d, aes(One , Two))

g + geom_point(aes(color = Genes , size = 0.1)) + labs(title = "") + ylab("Primary") + xlab("Metastatic") + ylim(5,14) + xlim(5,14) +
  theme(axis.title=element_text(size=15 , face  = "bold") , axis.text=element_text(size=16 , colour = "black") , plot.title = element_text(hjust = 0.5) ,
        title = element_text(size=16) , legend.text = element_text(size = 12),legend.title = element_text(size = 18)) + 
  scale_size_continuous(guide = F)  + guides(colour = guide_legend(override.aes = list(size=2, stroke=2))) +
  geom_label_repel(aes(label = NAME),
                   box.padding   = 0.5, 
                   point.padding = 0.5,
                   segment.color = 'grey50')








