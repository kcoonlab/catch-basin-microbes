set.seed(123)

library(phyloseq)
library(qiime2R)
library(microbiome)
library(vegan)
library(cluster)
library(dplyr)
library(clusterSim)
library(ade4)
library(philr)

# read in distance matrix
ps.final <- readRDS("catch-basin-microbes/input-files/ps.final.rds")

# Euclidean distance of philr transformation
# take phylogenetic info into account with isometric log-ratio transformation
# ILR (philr) transform then make distance matrix
ps.final.pseudocount <- ps.final
otu_table(ps.final.pseudocount) <- otu_table(ps.final) + 1
ps_philr <- philr(ps.final.pseudocount, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')
dist_philr <- dist(ps_philr, method="euclidean") 
otu.df <- data.frame(otu_table(ps.final.pseudocount))

# setting function to do the pam clustering
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

dist_matrix <- dist_philr
# how many clusters needed
data.cluster=pam.clustering(dist_matrix, k=3) # start with any number of clusters
nclusters = index.G1(t(otu.df), data.cluster, d = dist_matrix, centrotypes = "medoids") 
nclusters
nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(dist_matrix, k)
    nclusters[k]=index.G1(t(otu.df),data.cluster_temp,  d = dist_matrix,
                          centrotypes = "medoids")
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

# change number of clusters based on plot
data.cluster=pam.clustering(dist_matrix, k=2)
obs.silhouette=mean(silhouette(data.cluster, dist_matrix)[,3])

# plotting within clustering
obs.pcoa=dudi.pco(dist_matrix, scannf=F, nf=3) #scannf= whether to display eigenvalues bar plot, nf= axes to keep if scannf is false
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F) 
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4,7)) 

# retrieving assigned cluster
cluster.assignments <- cbind(data.cluster,obs.pcoa$li)
