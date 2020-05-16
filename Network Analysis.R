
setwd("F:/Systems Biomedicine of Cancer/Results")
load("network analysis.RData")

options(stringsAsFactors = F)


DEGsPvsN = read.table("DEGsPvsN.txt", sep = "\t")
DEGsMvsP = read.table("DEGsMvsP.txt", sep = "\t")
DEGsMvsN = read.table("DEGsMvsN.txt", sep = "\t")

DEGs_liver_normal = read.table("DEGs_liver_Normal.txt", sep = "\t")
DEGs_liver_primary = read.table("DEGs_liver_primary.txt", sep = "\t")
DEGs_lung_normal = read.table("DEGs_lung_normal.txt", sep = "\t")
DEGs_lung_primary = read.table("DEGs_lung_primary.txt", sep = "\t")
DEGs_primary_normal = read.table("DEGs_primary_normal.txt", sep = "\t")


DEGs = unique(c(DEGs_liver_normal$ID,DEGs_liver_primary$ID,DEGs_lung_normal$ID,DEGs_lung_primary$ID,DEGs_primary_normal$ID,
                DEGsPvsN$ID,DEGsMvsP$ID,DEGsMvsN$ID))
write.table(DEGs , "DEGs.txt" , sep = "\t" , quote = F , col.names = F , row.names = F)


write.table(names(V(g)) , "Giant Component.txt" , sep = "\t" , quote = F , col.names = F , row.names = F)


M_N = intersect(intersect(DEGs_lung_normal$ID , DEGs_liver_normal$ID ) , DEGsMvsN$ID)
M_P = intersect(intersect(DEGs_liver_primary$ID , DEGs_lung_primary$ID) , DEGsMvsP$ID)
P_N = intersect(DEGs_primary_normal$ID , DEGsPvsN$ID)

intersect(M_N , M_P)
intersect(M_N , P_N)
#write.table(intersect(M_N , P_N) , "intersect(M_N , P_N).txt" , sep = "\t" , quote = F , col.names = F , row.names = F)
intersect(M_P , P_N)

intersect(intersect(M_N , M_P) , P_N)
#"ETHE1"

DEGs = unique(c(M_N , M_P , P_N)) 
#write.table(DEGs , "DEGs.txt" , sep = "\t" , quote = F , col.names = F , row.names = F)

#write.table(M_P , "seedgenes.txt" , sep = "\t" , quote = F , col.names = F , row.names = F)


#BiocManager::install("STRINGdb")



setwd("E:/Piran/Systems Biomedicine SRC Project/Results")




library(igraph)
library(WGCNA)
library(readxl)


options(stringsAsFactors = F)
load(".RData")


mapseeds = read.table("seed_mapping.tsv.txt" , sep = "\t" , header = F)
seedgenes = as.character(mapseeds$V3)

# among 44 seedgenes, 35 of them mapped to string ID

intractions = read.table("string_interactions.tsv.txt" , sep = "\t" , header = F)

dim(intractions)

Network = as.data.frame(intractions[,c(5,6,15)])
colnames(Network) = c("from","to","score")
dim(Network)
# 556 edge   209 node   

seedgenes %in% unique(c(Network$to , Network$from))
# form 13, 12 seed genes are among the interacting genes



#########################################################
### Assigning weight STRING scores to adjacency matrix ##
#########################################################

g = graph.data.frame(Network[,1:2] , directed = T)
#plot(g , vertex.size = 5 , vertex.label.cex = 0.3)
g
is.connected(g)
c= components(g)
g = induced.subgraph(g , which(c$membership == 1 ))
g
#components(g)
edgelist = as_edgelist(g)
edgelist = as.data.frame(edgelist)
colnames(edgelist) = c("from","to")
length(unique(names(V(g))))
#205

e1 = Network[,1:2]
e1 = apply(e1, 1, paste , collapse = "")
e2 = apply(edgelist, 1, paste , collapse = "") 

commonIndex  = which(e1 %in% e2)

Network = Network[commonIndex,]
Network = as.data.frame(Network)
mapseeds = mapseeds[mapseeds$V3 %in% unique(c(Network$from, Network$to)) , ]



















setwd("F:/Colorectal Cancer/New")
load("network analysis.RData")

g = graph.data.frame(Network[,1:2] , directed = F)
E(g)$weigth = Network[,3]  # assigning weights
w.adj = as_adj(g , attr = "weigth")
w.adj = as.matrix(w.adj)
w.adj = w.adj / 10

summary(w.adj[which(as.vector(w.adj) > 0)])
################################################## 
## Using Weighted TOM similarity adjacency matrix
## in order to reduce the week edges according
## to weights.
################################################## 
symetric.w.adj = w.adj + t(w.adj)


TOM.w.adj = TOMsimilarity(symetric.w.adj)
diag(TOM.w.adj) = 0
rownames(TOM.w.adj) = rownames(symetric.w.adj)
colnames(TOM.w.adj) = rownames(symetric.w.adj)


TOM.FOR.SHORTEST.PATH = 1 - TOM.w.adj


max(degree(graph_from_edgelist(as.matrix(Network[,1:2]))))

###################################################
## Computing shortest path length
###################################################
graph = graph_from_adjacency_matrix(TOM.FOR.SHORTEST.PATH,
                                    weighted = TRUE,mode = "undirected")
is.connected(graph)
components(graph)
SP = shortest.paths(graph = graph , weight = E(graph)$weight)
rownames(SP) = rownames(symetric.w.adj)
colnames(SP) = rownames(symetric.w.adj)


WholeNet = SP
length(unique(rownames(SP)))
length(unique(colnames(SP)))


dim(mapseeds)
indx = numeric(12)
for(i in 1:12){
  indx[i]= which(mapseeds$V3[i] == rownames(WholeNet))
}


################################################
## Computing Score based on shortest path length
## for neighborhood finding.
################################################

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
#write.table(S , file = "ScoreNodeBasedOnSP.txt" , sep = "\t")

idx = which(S > 0)

match(indx,idx)
idx1 = c(indx,idx)

shortnet = w.adj[idx1,idx1]
idx2 = which(shortnet>0)
which(shortnet<0)
shortnet[idx2] = 1
dim(shortnet)


#n = rownames(shortnet)
#string_db$plot_network(n)
#write.table(a[,1] , file = "genes_shortnet.txt" , sep = "\t" , col.names = F, row.names = F, quote = F)

a = intractions[,c(1,5)]
b = intractions[,c(2,6)]
rownames(b) = 557:(length(b[,1])*2)
colnames(b) = c("V1","V5")
a = rbind(a,b)
a = a[!duplicated(a),]

a = a[ a[,2] %in% rownames(shortnet) , ]

a = a[match(colnames(shortnet) , a[,2]) , ]

identical(colnames(shortnet) , rownames(shortnet))

dimnames(shortnet) = list(a[,1],a[,1])


library(igraph)

net = graph_from_adjacency_matrix(shortnet , mode = "undirected" , weighted = NULL)
net
#net = NET


c= components(net)
c
net = induced.subgraph(net , which(c$membership == 7 ))
net
edgelist = as_edgelist(net)

seedgenes = mapseeds[(mapseeds$V2 %in% names(V(net))) , ][,2]
# "TMPO"  "CDC6"  "AIMP1" "CD74"
boxplot(mat , pch = "." , las = 3)

a = b$`sort(b, decreasing = T)`[match(names(degree) , rownames(b))  ]

write.table(names(V(net)) , file = "a.txt" , sep = "\t" , col.names = F , row.names = F, quote = F)
write.table(d , file = "d.txt" , sep = "\t"  , row.names = T, col.names = F ,  quote = F)

#####################
# reading network ###
# reading network ###
# reading network ###
#####################
library(igraph)

#net = graph_from_adjacency_matrix(shortnet , mode = "undirected" , weighted = T)
#net = delete_edge_attr(net , "weight")



deg <- degree(net,mode = "all")

plot(net ,vertex.size = 8,vertex.label.cex=0.5 , edge.width = 0.1,layout=layout.kamada.kawai)
plot(net , vertex.size = 8,vertex.label.dist=0,vertex.label.cex=0.5 , edge.width = 0.1,layout=layout.kamada.kawai,
     vertex.color=1:121)
plot(net , vertex.size = deg/3,vertex.label.dist=0,vertex.label.cex=0.5 , edge.width = 0.1,layout=layout.kamada.kawai,
     vertex.color=1:121)

plot(net ,vertex.size = 8,vertex.label.cex=0.3 , edge.width = 0.1,layout=layout.circle)

plot(net , vertex.size = deg/3,vertex.label.dist=0,vertex.label.cex=0.4 , edge.width = 0.1,layout=layout.circle,
     vertex.color=1:121)

##############################
##  Network Descriptives #####
##############################

transitivity(net, type="global") 
0.5909091
Trantivity = transitivity(net, type="local")

Transivity = data.frame(names(V(net)) , Trantivity)

V(net)$"ITGA5"
diameter(net, directed=F)
8
#Density: The ratio of the number of edges and the number of possible edges.

edge_density(net, loops=F)
0.1809524

#for an undirected network
2*ecount(as.undirected(net))/(vcount(as.undirected(net))*(vcount(as.undirected(net))-1))


diam <- get_diameter(net, directed=F)
diam
# TRIM31 HLA-F  CD74   PLAGL2 TMPO   MAD2L1 PSMA7  AIMP1  TWF1

diam = as.vector(diam)



vcol <- rep("aquamarine", vcount(net))

vcol[diam] <- "red"

ecol <- rep("gray", ecount(net))

ecol[E(net, path=diam)] <- "red" 
par(mfrow = c(1,1))
# E(net, path=diam) finds edges along a path, here 'diam'







vcol <- rep("aquamarine", vcount(net))

idx = which(names(V(net)) %in% seedgenes)

vcol[idx] <- "red"

ecol <- rep("gray", ecount(net))

ecol[E(net, path=diam)] <- "red" 
par(mfrow = c(1,1))
# E(net, path=diam) finds edges along a path, here 'diam'






layouts = c("layout.random",
            "layout.circle",
            "layout.sphere",
            "layout.fruchterman.reingold",
            "layout.kamada.kawai",
            "layout.spring",
            "layout.reingold.tilford",
            "layout.fruchterman.reingold.grid",
            "layout.lgl",
            "layout.graphopt",
            "layout.svd")

plot(net, vertex.color=vcol, edge.color=ecol, edge.curved= 0.2,vertex.label.cex=0.6, vertex.label.color="blue4" , vertex.size = 15 ,layout = layout.fruchterman.reingold)
plot(net, vertex.color=vcol, edge.color=ecol, edge.curved= 0.2,vertex.label.cex=0.6, vertex.label.color="darkblue" ,vertex.size = 15,layout = layout.kamada.kawai)




library(igraph)
##########################
## Centralities ##########
##########################
##########################
## Centralities ##########
##########################
#degree of nodes
degree <- degree(net,mode = "all")

degree <- degree(g,mode = "all")

degree = sort(degree , decreasing = T)
degree = as.data.frame(degree)

hist(degree, breaks=1:vcount(net)-1, main="Histogram of degree distribution",col = "gray90")


closenes = closeness(net, mode="all", normalized=T) 
closenes = round(as.data.frame(sort(closenes , decreasing = T) ) , 4)




centr_clo(net, mode="all", normalized=T)
#total closeness of the network : 0.305172


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

#####################
## distances ########
#####################

mean_distance(net, directed=F)
3.47619
d = distances(net)

d = as.data.frame(round(sort(apply(d,1,mean) , decreasing = F) , 4))



#################################
## pathway enrichment analysis ##
#################################
#################################
## pathway enrichment analysis ##
#################################
library(igraph)
library(ReactomePA)
library(org.Hs.eg.db)
library(DOSE)

net = graph_from_adjacency_matrix(shortnet , mode = "undirected" , weighted = NULL)
net
#net = NET


c= components(net)
c
net = induced.subgraph(net , which(c$membership == 7 ))
net

#write.table(names(V(net)) , file = "a.txt" , sep = "\t" , col.names = F , row.names = F, quote = F)


x = org.Hs.egSYMBOL2EG
x1 = as.list(x)
x2 =x1[names(V(net))]
x2 = unlist(x2)
x2
m1 =unname(x2)



e1 = enrichPathway(m1, organism = "human", pvalueCutoff = 0.05)
e1
path1 = as.data.frame(e1)
path1 = path1[,c(1,2,6,8,9)]

#write.table(path1 , file = "path1.txt" , sep = "\t" , col.names = T , row.names = F, quote = F)
#m1[which(m1 %in% seeds2)]

cnetplot(e1  , categorySize="pvalue", vertex.label.cex = 0.1 , circular =T)
