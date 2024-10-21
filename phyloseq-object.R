library(phytools)

pick_new_outgroup <- function(tree.unrooted){
  require(magrittr)
  require(data.table)
  require(ape) # ape::Ntip
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
unrooted_qiime_default<-read.newick(file="input-files/tree.nwk")
new.outgroup = pick_new_outgroup(unrooted_qiime_default)
rootedTree = ape::root(unrooted_qiime_default, outgroup=new.outgroup, resolve.root=TRUE)

# creating a phyloseq object from the qiime outputs

library(phyloseq)
library(ggplot2)
library(qiime2R)

newlyrootedTree<-read_tree(rootedTree, errorIfNULL=FALSE)

ps.all <- qza_to_phyloseq("input-files/table_analysis.qza", 
                           "input-files/rooted_tree.qza", 
                           "input-files/taxass_taxonomy_98.qza", 
                           "input-files/metadata.txt")

# will create a new phyloseq object using a rerooted tree (see table 2 script)
# this will be used for all downstream analyses

ps.all <- phyloseq(otu_table(ps.all), tax_table(ps.all), sample_data(ps.all), newlyrootedTree)
saveRDS(ps.all,"input-files/ps.all.rds")
