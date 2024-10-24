set.seed(123)

library(phyloseq)
library(phytools)
library(ggplot2)
library(qiime2R)
library(vegan)

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
unrooted_qiime_default<-read.newick(file="catch-basin-microbes/input-files/tree.nwk")
new.outgroup = pick_new_outgroup(unrooted_qiime_default)
rootedTree = ape::root(unrooted_qiime_default, outgroup=new.outgroup, resolve.root=TRUE)
newlyrootedTree<-read_tree(rootedTree, errorIfNULL=FALSE)

ps.temp <- qza_to_phyloseq("catch-basin-microbes/input-files/table_analysis.qza", 
                           "catch-basin-microbes/input-files/rooted_tree.qza", 
                           "catch-basin-microbes/input-files/taxass_taxonomy_98.qza", 
                           "catch-basin-microbes/input-files/metadata.txt")

# Create a new phyloseq object using a rerooted tree
# this will be used for all downstream analyses

ps.temp.rerooted <- phyloseq(otu_table(ps.temp), tax_table(ps.temp), sample_data(ps.temp), newlyrootedTree)
ps.final <- subset_taxa(ps.temp.rerooted, Family != "Mitochondria")
saveRDS(ps.final,"catch-basin-microbes/input-files/ps.final.rds")

shannon <- diversity(t(otu_table(ps.final)), "shannon")
asvfac = factor(rownames(tax_table(ps.final)))
asvtab = apply(otu_table(ps.final), MARGIN = 2, function(x) {
    tapply(x, INDEX = asvfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
observationThreshold = 1
richness <- apply(asvtab > observationThreshold, 2, sum)
