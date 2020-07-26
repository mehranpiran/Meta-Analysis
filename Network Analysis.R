
#######################
## Retrieving DEGs ####
#######################

setwd("")
options(stringsAsFactors = F)

DEGsPvsN = read.table("DEGsPvsN.txt", sep = "\t")
DEGsMvsP = read.table("DEGsMvsP.txt", sep = "\t")
DEGsMvsN = read.table("DEGsMvsN.txt", sep = "\t")

DEGs_liver_normal = read.table("DEGs_liver_Normal.txt", sep = "\t")
DEGs_liver_primary = read.table("DEGs_liver_primary.txt", sep = "\t")
DEGs_lung_normal = read.table("DEGs_lung_normal.txt", sep = "\t")
DEGs_lung_primary = read.table("DEGs_lung_primary.txt", sep = "\t")
DEGs_primary_normal = read.table("DEGs_primary_normal.txt", sep = "\t")


AllDEGs = unique(c(DEGs_liver_normal$ID,DEGs_liver_primary$ID,DEGs_lung_normal$ID,DEGs_lung_primary$ID,DEGs_primary_normal$ID,
                DEGsPvsN$ID,DEGsMvsP$ID,DEGsMvsN$ID))

write.table(AllDEGs , "AllDEGs.txt" , sep = "\t" , quote = F , col.names = F , row.names = F)


M_N = intersect(intersect(DEGs_lung_normal$ID , DEGs_liver_normal$ID ) , DEGsMvsN$ID)
M_P = intersect(intersect(DEGs_liver_primary$ID , DEGs_lung_primary$ID) , DEGsMvsP$ID)
P_N = intersect(DEGs_primary_normal$ID , DEGsPvsN$ID)


DEGs = unique(c(M_N , M_P , P_N)) 
length(GEGs)
# 334




###################################################################################################
### Network Analysis: Computing shortest path based scores and constructing the core network ######
###################################################################################################

## Constructing the PPI network

library(igraph)
library(WGCNA)

setwd("")
intractions = read.table("string_interactions.tsv.txt" , sep = "\t" , header = F)
Network = as.data.frame(intractions[,c(1,2,15)])
colnames(Network) = c("from","to","score")  

g = graph.data.frame(Network[,1:2] , directed = T)

is.connected(g)
# False
c= components(g)
g = induced.subgraph(g , which(c$membership == 1 ))
edgelist = as_edgelist(g)
edgelist = as.data.frame(edgelist)
colnames(edgelist) = c("from","to")

e1 = Network[,1:2]
e1 = apply(e1, 1, paste , collapse = "")
e2 = apply(edgelist, 1, paste , collapse = "") 
commonIndex  = which(e1 %in% e2)
Network = Network[commonIndex,]
Network = as.data.frame(Network)





## Assigning weight STRING scores to adjacency matrix

setwd("F:/Colorectal Cancer/New")
load("network analysis.RData")

g = graph.data.frame(Network[,1:2] , directed = F)
E(g)$weigth = Network[,3]  # assigning weights
w.adj = as_adj(g , attr = "weigth")
w.adj = as.matrix(w.adj)
w.adj = w.adj / 10
summary(w.adj[which(as.vector(w.adj) > 0)])



## Using Weighted TOM similarity adjacency matrix in order to reduce the week edges according to weights.

symetric.w.adj = w.adj + t(w.adj)

TOM.w.adj = TOMsimilarity(symetric.w.adj)
diag(TOM.w.adj) = 0
rownames(TOM.w.adj) = rownames(symetric.w.adj)
colnames(TOM.w.adj) = rownames(symetric.w.adj)

TOM.FOR.SHORTEST.PATH = 1 - TOM.w.adj



## Computing shortest path lengths
graph = graph_from_adjacency_matrix(TOM.FOR.SHORTEST.PATH , weighted = TRUE,mode = "undirected")
SP = shortest.paths(graph = graph , weight = E(graph)$weight)
rownames(SP) = rownames(symetric.w.adj)
colnames(SP) = rownames(symetric.w.adj)


WholeNet = SP
length(unique(rownames(SP)))
length(unique(colnames(SP)))


coregenes = read.table("core_mapping.tsv.txt" , sep = "\t" , header = F)
coregenes = as.character(coregenes$V3)

length(coregenes)
# 12
indx = numeric(12)
for(i in 1:length(coregenes)){
  indx[i]= which(coregenes[i] == rownames(WholeNet))
}


## Computing shortest path based scores for core genes neighborhoods.

SeedNet= WholeNet[indx,indx]
WholeNetMinus = WholeNet[-indx,-indx]
InterSeedNet = WholeNet[indx,]

dim(WholeNet)
dim(WholeNetMinus)
dim(InterSeedNet)

S = numeric(205)
for (j in 1:205){
  M=sum(WholeNet[,j])/205
  i =  which(rownames(WholeNetMinus) == rownames(WholeNet)[j], TRUE)
  Cprim = sum(WholeNetMinus[i,])/193
  i =  which( colnames(InterSeedNet) == rownames(WholeNet)[j], TRUE)
  C=sum(InterSeedNet[,i])/12
  S[j]=(Cprim-C)/M
}


idx = which(S > 0)
match(indx,idx)
idx1 = c(indx,idx)

CoreNetwork = w.adj[idx1,idx1]
idx2 = which(CoreNetwork>0)
CoreNetwork[idx2] = 1
dim(CoreNetwork)
# 39 39









#############################
## Network Descriptives #####
#############################

library(igraph)

net = graph_from_adjacency_matrix(CoreNetwork , mode = "undirected" , weighted = NULL)

c= components(net)
net = induced.subgraph(net , which(c$membership == 7 ))
# net is the giant component of the core network
edgelist = as_edgelist(net)


plot(net ,vertex.size = 8,vertex.label.cex=0.5 , edge.width = 0.1,layout=layout.kamada.kawai)
plot(net , vertex.size = 8,vertex.label.dist=0,vertex.label.cex=0.5 , edge.width = 0.1,layout=layout.kamada.kawai,
     vertex.color=1:121)


transitivity(net, type="global") 
0.5909091
Trantivity = transitivity(net, type="local")

Transivity = data.frame(names(V(net)) , Trantivity)

edge_density(net, loops=F)
0.1809524

2*ecount(as.undirected(net))/(vcount(as.undirected(net))*(vcount(as.undirected(net))-1))
0.1809524

diameter(net, directed=F)
8
diam <- get_diameter(net, directed=F)
diam
# TRIM31 HLA-F  CD74   PLAGL2 TMPO   MAD2L1 PSMA7  AIMP1  TWF1




# color the core genes in red

idx = which(names(V(net)) %in% coregenes)

vcol[idx] <- "red"

ecol <- rep("gray", ecount(net))

ecol[E(net, path=diam)] <- "red" 
par(mfrow = c(1,1))

plot(net, vertex.color=vcol, edge.color=ecol, edge.curved= 0.2,vertex.label.cex=0.6, vertex.label.color="blue4" , vertex.size = 15 ,layout = layout.fruchterman.reingold)
plot(net, vertex.color=vcol, edge.color=ecol, edge.curved= 0.2,vertex.label.cex=0.6, vertex.label.color="darkblue" ,vertex.size = 15,layout = layout.kamada.kawai)





## Compute centralities

degree <- degree(net,mode = "all")

degree = sort(degree , decreasing = T)
degree = as.data.frame(degree)

hist(degree[,1], breaks=1:vcount(net)-1, main="Histogram of degree distribution",col = "gray90")


closenes = closeness(net, mode="all", normalized=T) 
closenes = round(as.data.frame(sort(closenes , decreasing = T) ) , 4)
centr_clo(net, mode="all", normalized=T)
# total closeness: 0.305172


b = betweenness(net, directed=F, weights=NA)
b = round(as.data.frame(sort(b , decreasing = T)) , 2)
# total betweenness centrality: 0.1051746

EdgeBetweenness = data.frame(edgelist , edge_betweenness(net, directed=F, weights=NA))
colnames(EdgeBetweenness) = c("Gene1" , "Gene2" , "Edge_Betweenness")

idx = which(duplicated(EdgeBetweenness$Edge_Betweenness))
EdgeBetweenness[idx, ][,3] = EdgeBetweenness[idx, ][,3] + 0.000001

EdgeBetweenness = EdgeBetweenness[match(sort(EdgeBetweenness$Edge_Betweenness , decreasing = T) , EdgeBetweenness$Edge_Betweenness) , ]
EdgeBetweenness$Edge_Betweenness = round(EdgeBetweenness$Edge_Betweenness , 2)
E(net)
duplicated(EdgeBetweenness[,1:2])


## Distances 
mean_distance(net, directed=F)
3.47619
d = distances(net)

d = as.data.frame(round(sort(apply(d,1,mean) , decreasing = F) , 4))





#################################
## pathway enrichment analysis ##
#################################

library(igraph)
library(ReactomePA)
library(org.Hs.eg.db)
library(DOSE)

net = graph_from_adjacency_matrix(CorNetwork , mode = "undirected" , weighted = NULL)


x = org.Hs.egSYMBOL2EG
x1 = as.list(x)
x2 =x1[names(V(net))]
x2 = unlist(x2)
m1 =unname(x2)


e1 = enrichPathway(m1, organism = "human", pvalueCutoff = 0.05)
path1 = as.data.frame(e1)
path1 = path1[,c(1,2,6,8,9)]

cnetplot(e1  , categorySize="pvalue", vertex.label.cex = 0.1 , circular =T)
