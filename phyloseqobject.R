
# creating a phyloseq object from the qiime outputs

library(phyloseq)
library(ggplot2)
library(qiime2R)

ps.all <- qza_to_phyloseq("table_analysis.qza", 
                           "rooted_tree.qza", 
                           "taxass_taxonomy_98.qza", 
                           "metadata.txt")
ps.all

# will create a new phyloseq object using a rerooted tree (see table 2 script)
# this will be used for all downstream analyses
