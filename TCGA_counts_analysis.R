
options(stringsAsFactors = F)
setwd("F:/Systems Biomedicine of Cancer/Confirmation Data/TCGAcounts")
load("1.RData")
files = list.files("." , "*.txt")
Exprtable = lapply(files, read.delim , header=F)
Exprtable = do.call(cbind , Exprtable)
rownames(Exprtable) = Exprtable[,1]
Exprtable = Exprtable[,-seq(1 , ncol(Exprtable) , 2)]
colnames(Exprtable) = sub(".counts.txt" , "" , files)
#Exprtable = apply(Exprtable, 2,as.integer)
Exprtable =  Exprtable[,c(1:4,10:14,5:9)]
G_list = read.delim("G_list.txt")

Exprtable$ensembl_gene_id = substr(rownames(Exprtable) , 1 , 15)

Exprtable = merge(G_list , Exprtable, by = "ensembl_gene_id")
any(duplicated(Exprtable$hgnc_symbol))

rownames(Exprtable) = Exprtable$hgnc_symbol
Exprtable = Exprtable[,-c(1,2)]



library(DESeq2)
library(ggplot2)
library(limma)


mat = Exprtable[,c(5:9,10:14)]
mat = mat[!(rowSums(mat[,1:length(mat[1,])]) <= 5) , ]
gr = factor(c(rep("T" , 5) , rep("N" , 5)))
colData = data.frame(group=gr , type = "paired-end")
cds = DESeqDataSetFromMatrix(countData = mat , colData = colData , design = ~group)

cds = DESeq(cds)

#cnt = log2(1 + counts(cds , normalized = T))
#cnt

Table = results(cds , c("group" , "T" , "N"))

Table = as.data.frame(Table)



write.table(Table, file = "DEGsMvsN.txt" , sep = "\t" , quote = F)
write.table(Table, file = "DEGsMvsP.txt" , sep = "\t" , quote = F)
write.table(Table, file = "DEGsPvsN.txt" , sep = "\t" , quote = F)




#dif$padj = p.adjust(dif$pvalue , method = "BH")

#dif = dif[order(dif$padj) , ]

#ggplot(dif , aes(log2FoldChange , -log10(padj) , color = log2FoldChange )) + geom_point() + theme_bw()

