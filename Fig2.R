## Taxa bar plots

ps.genus.rerooted <- tax_glom(ps.all.rerooted, "Genus", NArm = TRUE)
ps.family.rerooted <- tax_glom(ps.all.rerooted, "Family", NArm = TRUE)
ps.order.rerooted <- tax_glom(ps.all.rerooted, "Order", NArm = TRUE)
ps.class.rerooted <- tax_glom(ps.all.rerooted, "Class", NArm = TRUE)
ps.phylum.rerooted <- tax_glom(ps.all.rerooted, "Phylum", NArm = TRUE)
saveRDS(ps.genus.rerooted, "ps.genus.rerooted")
saveRDS(ps.family.rerooted, "ps.family.rerooted")
saveRDS(ps.order.rerooted, "ps.order.rerooted")
saveRDS(ps.class.rerooted, "ps.class.rerooted")
saveRDS(ps.phylum.rerooted, "ps.phylum.rerooted")

ps.phylum <- tax_glom(ps.all.rerooted, "Phylum", NArm = TRUE)
ps.family <- tax_glom(ps.all.rerooted, taxrank = 'Family') # agglomerate taxa
ps.genus <- tax_glom(ps.all.rerooted, taxrank = 'Genus') # agglomerate taxa
saveRDS(ps.phylum, "ps_phylum_rerooted.rds")
saveRDS(ps.family, "ps_family_rerooted.rds")
saveRDS(ps.genus, "ps_genus_rerooted.rds")

ps.phylum <- readRDS("ps_phylum_rerooted.rds")
ps.family <- readRDS("ps_family_rerooted.rds")
ps.genus <- readRDS("ps_genus_rerooted.rds")
metadata <- read.table("metadata.txt", sep="\t", header=TRUE)
metadata.ps <- sample_data(metadata)
rownames(metadata.ps) <- metadata.ps$sample.id
ps.phylum <- phyloseq(otu_table(ps.phylum),tax_table(ps.phylum),metadata.ps, phy_tree(ps.phylum))
ps.family <- phyloseq(otu_table(ps.family),tax_table(ps.family),metadata.ps, phy_tree(ps.family))
ps.genus <- phyloseq(otu_table(ps.genus),tax_table(ps.genus),metadata.ps, phy_tree(ps.genus))
ps.phylum
ps.family
ps.genus

library(RColorBrewer)

sample_data(ps.input)$date %>% class
date <- get_variable(ps.input,"date_code")
ple_data(ps.input)$cluster_philr <- as.factor(sample_data(ps.input)$cluster_philr)
sample_data(ps.input)$UCB <- sample_data(ps.input)$UCB %>% as.factor
sample_data(ps.input)$combined_separate <- sample_data(ps.input)$combined_separate %>% as.factor


# Fig 2A. Taxa bar plots by phylum, per basin
ps.input <- ps.phylum
y1 <- ps.input
(y2 = merge_samples(y1, "UCB")) # merge samples on sample variable of interest: cluster, basin, date, NewPastedVar
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Sample <- factor(y4$Sample, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21")) # only if grouping by date
y4$Phylum[y4$Abundance < 0.01] <- "Taxa < 1% abund." #rename genera with < 1% abundance
colourCount = length(unique(y4$Phylum))
colourCount
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
as.factor(y4$Phylum) %>% levels
# by UCB
y4$Phylum <- factor(y4$Phylum,levels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Campilobacterota","Cyanobacteria","Deinococcota","Desulfobacterota","Firmicutes","Fusobacteriota","Proteobacteria","Spirochaetota","unclassified.Bacteria","Verrucomicrobiota","Taxa < 1% abund."))
#ucb
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  facet_wrap(~combined_separate, strip.position="bottom", scales="free_x") +
  theme(panel.background = element_blank(), legend.position="right", axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank()) + 
  guides(fill=guide_legend(ncol=1))+
  xlab(NULL) + ylab("Relative abundance")


# Fig 2B. Taxa bar plots by genus, per basin

ps.input <- ps.genus
y1 <- ps.input
(y2 = merge_samples(y1, "UCB")) # merge samples on sample variable of interest: cluster, basin, date, NewPastedVar
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Genus <- as.character(y4$Genus) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Sample <- factor(y4$Sample, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21")) # only if grouping by date
y4$Genus[y4$Abundance < 0.03] <- "Taxa < 3% abund." #rename genera with < 1% abundance. <3% genus ucb
colourCount = length(unique(y4$Genus))
colourCount
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
as.factor(y4$Genus) %>% levels
# by UCB
y4$Genus <- factor(y4$Genus,levels=c("acIII-A","Acinetobacter","Aeromonas","Arcobacter","Arenimonas","bacII-A","Bacteroides","betI-A","betIII-A","betVII-A","C39","Dechloromonas","Hydrogenophaga","Insolitispirillum","Lactococcus","Malikia","Nevskia","Novispirillum","Pseudarcobacter","Pseudomonas","Sulfuricurvum","Sulfurospirillum","Thauera","Zoogloea","unclassified.Alcaligenaceae","unclassified.alfVI","unclassified.alfVII","unclassified.Bacteroidales","unclassified.betI","unclassified.Comamonadaceae","unclassified.Enterobacterales","unclassified.Gammaproteobacteria","unclassified.Geobacteraceae","unclassified.Lachnospiraceae","unclassified.Methylococcaceae","unclassified.Microbacteriaceae","unclassified.Neisseriaceae","unclassified.Prevotellaceae","unclassified.Rhizobiaceae","unclassified.Rhodocyclaceae","unclassified.Veillonellales-Selenomonadales","Taxa < 3% abund."))
y4$Genus <- factor(y4$Genus,levels=c("acIII-A","Acinetobacter","Aeromonas","Arcobacter","Arenimonas","bacII-A","Bacteroides","betI-A","betIII-A","betVII-A","C39","Dechloromonas","Hydrogenophaga","Insolitispirillum","Lactococcus","Malikia","Nevskia","Novispirillum","Pseudarcobacter","Pseudomonas","Sulfuricurvum","Sulfurospirillum","Thauera","unclassified.Alcaligenaceae","unclassified.alfVI","unclassified.alfVII","unclassified.Bacteroidales","unclassified.betI","unclassified.Comamonadaceae","unclassified.Enterobacterales","unclassified.Gammaproteobacteria","unclassified.Geobacteraceae","unclassified.Lachnospiraceae","unclassified.Methylococcaceae","unclassified.Microbacteriaceae","unclassified.Neisseriaceae","unclassified.Prevotellaceae","unclassified.Rhizobiaceae","unclassified.Rhodocyclaceae","unclassified.Veillonellales-Selenomonadales","Zoogloea","Taxa < 3% abund."))
# UCB
p + geom_bar(aes(fill=Genus), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  facet_wrap(~combined_separate, strip.position="bottom", scales="free_x") +
  # facet_grid(.~combined_separate, scales="free", switch="x", space="free_x") +
  theme(panel.background = element_blank(), legend.position="right", axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank()) + 
  guides(fill=guide_legend(ncol=2))+ #UCB 2 cols, date 1 col
  xlab(NULL) + ylab("Relative abundance")



# Fig 2C. Taxa bar plots by phylum, per date
ps.input <- ps.phylum
y1 <- ps.input
(y2 = merge_samples(y1, "date")) # merge samples on sample variable of interest: cluster, basin, date, NewPastedVar
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Sample <- factor(y4$Sample, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21")) # only if grouping by date
y4$Phylum[y4$Abundance < 0.01] <- "Taxa < 1% abund." #rename genera with < 1% abundance
colourCount = length(unique(y4$Phylum))
colourCount
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
as.factor(y4$Phylum) %>% levels
# by date
y4$Phylum <- factor(y4$Phylum,levels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Campilobacterota","Cyanobacteria","Deinococcota","Desulfobacterota","Firmicutes","Myxococcota","Proteobacteria","Spirochaetota","Verrucomicrobiota","Taxa < 1% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
#date, cluster
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=1))+ #UCB 2 cols, date 1 col
  xlab(NULL) + ylab("Relative abundance")


# Fig 2D. Taxa bar plots by genus, per date
ps.input <- ps.genus
y1 <- ps.input
(y2 = merge_samples(y1, "date")) # merge samples on sample variable of interest: cluster, basin, date, NewPastedVar
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Genus <- as.character(y4$Genus) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Sample <- factor(y4$Sample, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21")) # only if grouping by date
y4$Genus[y4$Abundance < 0.03] <- "Taxa < 3% abund." #rename genera with < 1% abundance. <3% genus ucb
colourCount = length(unique(y4$Genus))
colourCount
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
as.factor(y4$Genus) %>% levels
y4$Genus <- factor(y4$Genus,levels=c("Acinetobacter","Aeromonas","bacII-A","betI-A","betIII-A","betVII-A","C39","Dechloromonas","Hydrogenophaga","Malikia","Pseudarcobacter","Pseudomonas","unclassified.Alcaligenaceae","unclassified.alfVI","unclassified.betI","unclassified.Comamonadaceae","unclassified.Methylococcaceae","unclassified.Microbacteriaceae","unclassified.Prevotellaceae","unclassified.Veillonellales-Selenomonadales","Taxa < 3% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
p + geom_bar(aes(fill=Genus), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=1))+ #UCB 2 cols, date 1 col, cluster 2cols
  xlab(NULL) + ylab("Relative abundance")




