## Table 2.  Effects of sampling date, water quality, mosquito productivity on microbiota diversity

library("phyloseq")
library("qiime2R")
library("microbiome")
library("vegan")
library("cluster")
library("dplyr")
library("clusterSim")
library("ade4")

ps.all <- qza_to_phyloseq("table_analysis.qza", 
                          "mafft-fasttree-output/rooted_tree.qza", 
                          "taxass_taxonomy_98.qza", 
                          "metadata.txt")

# re-root fasttree tree # https://github.com/joey711/phyloseq/issues/597
# export mafft tree as newick file
qiime tools export \
--input-path tree.qza \
--output-path exported-tree

library(phytools)
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}
unrooted_qiime_default<-read.newick(file="mafft-fasttree-output/exported-tree/tree.nwk")
new.outgroup = pick_new_outgroup(unrooted_qiime_default)
rootedTree = ape::root(unrooted_qiime_default, outgroup=new.outgroup, resolve.root=TRUE)
newlyrootedTree<-read_tree(rootedTree, errorIfNULL=FALSE)
ps.all.rerooted <- phyloseq(otu_table(ps.all), tax_table(ps.all), sample_data(ps.all), newlyrootedTree)
saveRDS(ps.all.rerooted,"ps.all.rerooted.rds")
ps.all <- readRDS("ps.all.rerooted.rds")



## Table 2A  PERMANOVA tests using PhILR distances between samples

# running permanovas on whole dataset
# can't actually do this on variables that have any NAs.
# the variables without NAs are the ones that are always the same for each UCB
ps.all.pseudocount <- ps.all
otu_table(ps.all.pseudocount) <- otu_table(ps.all) + 1
ps.all.philr <- philr(ps.all.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
dist.all.philr <- dist(ps.all.philr, method="euclidean")
saveRDS(ps.all.philr, "philrforplots_ps.all.philr.rds")
saveRDS(dist.all.philr,"philrforplots_dist.all.philr.rds")
samples.all <- as.data.frame(as.matrix(sample_data(ps.all)))
# running with just individual NA samples cut
# (not all instances of a UCB with an NA anywhere)
NA %in% sample_data(ps.all)$pH #row 39
NA %in% sample_data(ps.all)$Temp[129] # row 39,129
NA %in% sample_data(ps.all)$DO[39] # row 39
NA %in% sample_data(ps.all)$Sal[71] # row 39, 71
NA %in% sample_data(ps.all)$Cond[39] # row 39
sample_data(ps.all)[129,] #sample.id: row 39=B_73, row 71=C_61, row129=F_56
# can take out the 3 individual samples that have NAs, and use that remaining set for the permanovas
ps.most <- subset_samples(ps.all, !sample.id %in% c("B_73","C_61","F_56"))
ps.all
ps.most
ps.most.pseudocount <- ps.most
otu_table(ps.most.pseudocount) <- otu_table(ps.most) + 1
ps.most.philr <- philr(ps.most.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
dist.most.philr <- dist(ps.most.philr, method="euclidean")
saveRDS(ps.most.philr, "philrforplots_ps.most.philr.rds")
saveRDS(dist.most.philr,"philrforplots_dist.most.philr.rds")
ps.most.philr <- readRDS("philrforplots_ps.most.philr.rds")
dist.most.philr <- readRDS("philrforplots_dist.most.philr.rds")
samples.most <- as.data.frame(as.matrix(sample_data(ps.most)))
perm <- how(nperm=999)
setBlocks(perm) <- with(samples.most, UCB)

adonis2(dist.most.philr ~ date, data=samples.most, permutations=perm)
adonis2(dist.most.philr ~ pH, data=samples.most, permutations=perm)
adonis2(dist.most.philr ~ Temp, data=samples.most, permutations=perm)
adonis2(dist.most.philr ~ DO, data=samples.most, permutations=perm)
adonis2(dist.most.philr ~ Sal, data=samples.most, permutations=perm)
adonis2(dist.most.philr ~ Cond, data=samples.most, permutations=perm)
adonis2(dist.most.philr ~ Pupae_pa, data=samples.most)
adonis2(dist.most.philr ~ Pupae, data=samples.most)


## Table 2B PERMANOVA tests using PhILR distances between samples aggregated by basin



ps.NAcut <- subset_samples(ps.all, !UCB %in% c("43","51","58","79","86","88","32","47","50","64","73","56","61"))
ps.UCBgroup <- merge_samples(ps.NAcut,"UCB")
ps.UCBgroup.pseudocount <- ps.UCBgroup
ps.UCBgroup.philr <- philr(ps.UCBgroup.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
dist.UCBgroup.philr <- dist(ps.UCBgroup.philr, method="euclidean")
saveRDS(dist.UCBgroup.philr,"philrforplots_dist.UCBgroup.philr.rds")
dist.UCBgroup.philr <- readRDS("philrforplots_dist.UCBgroup.philr.rds")

adonis2(dist.UCBgroup.philr ~ Flowgroup_coarse, data=samples.UCBgroup)
adonis2(dist.UCBgroup.philr ~ combined_separate, data=samples.UCBgroup)