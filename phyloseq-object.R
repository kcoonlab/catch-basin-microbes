library(phyloseq)
library(phytools)
library(ggplot2)
library(qiime2R)

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
newlyrootedTree<-read_tree(rootedTree, errorIfNULL=FALSE)

ps.temp <- qza_to_phyloseq("input-files/table_analysis.qza", 
                           "input-files/rooted_tree.qza", 
                           "input-files/taxass_taxonomy_98.qza", 
                           "input-files/metadata.txt")

# Create a new phyloseq object using a rerooted tree
# this will be used for all downstream analyses

ps.temp.rerooted <- phyloseq(otu_table(ps.temp), tax_table(ps.temp), sample_data(ps.temp), newlyrootedTree)
ps.final <- subset_taxa(ps.temp.rerooted, Family != "Mitochondria")
saveRDS(ps.final,"input-files/ps.final.rds")
