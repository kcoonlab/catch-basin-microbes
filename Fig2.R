## Fig 2. Bacterial diversity in water sampled from study catch basins

library(phyloseq)
library(RColorBrewer)

## Fig 2A. Taxa bar plots by basin (phylum level)

ps.all <- readRDS("input-files/ps.all.rds")
ps.phylum <- tax_glom(ps.all, "Phylum", NArm = TRUE)
y1 <- ps.phylum
y2 <- merge_samples(y1, "UCB") # merge samples by basin
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Phylum[y4$Abundance < 0.01] <- "Taxa < 1% abund." #rename genera with < 1% abundance
colourCount = length(unique(y4$Phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
y4$Phylum <- factor(y4$Phylum,levels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Campilobacterota","Cyanobacteria","Deinococcota","Desulfobacterota","Firmicutes","Fusobacteriota","Proteobacteria","Spirochaetota","unclassified.Bacteria","Verrucomicrobiota","Taxa < 1% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  facet_wrap(~combined_separate, strip.position="bottom", scales="free_x") +
  theme(panel.background = element_blank(), legend.position="right", axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank()) + 
  guides(fill=guide_legend(ncol=2)) +
  xlab(NULL) + ylab("Relative abundance")

## Fig 2B. Taxa bar plots by basin (genus level)

ps.genus <- tax_glom(ps.all, "Genus", NArm = TRUE)
y1 <- ps.genus
y2 <- merge_samples(y1, "UCB") # merge samples by basin
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Genus <- as.character(y4$Genus) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Genus[y4$Abundance < 0.03] <- "Taxa < 3% abund." #rename genera with < 3% abundance
colourCount = length(unique(y4$Genus))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
y4$Genus <- factor(y4$Genus,levels=c("acIII-A","Acinetobacter","Aeromonas","Arcobacter","Arenimonas","bacII-A","Bacteroides","betI-A","betIII-A","betVII-A","C39","Dechloromonas","Hydrogenophaga","Insolitispirillum","Lactococcus","Malikia","Nevskia","Novispirillum","Pseudarcobacter","Pseudomonas","Sulfuricurvum","Sulfurospirillum","Thauera","Zoogloea","unclassified.Alcaligenaceae","unclassified.alfVI","unclassified.alfVII","unclassified.Bacteroidales","unclassified.betI","unclassified.Comamonadaceae","unclassified.Enterobacterales","unclassified.Gammaproteobacteria","unclassified.Geobacteraceae","unclassified.Lachnospiraceae","unclassified.Methylococcaceae","unclassified.Microbacteriaceae","unclassified.Neisseriaceae","unclassified.Prevotellaceae","unclassified.Rhizobiaceae","unclassified.Rhodocyclaceae","unclassified.Veillonellales-Selenomonadales","Taxa < 3% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
p + geom_bar(aes(fill=Genus), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  facet_wrap(~combined_separate, strip.position="bottom", scales="free_x") +
  theme(panel.background = element_blank(), legend.position="right", axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank()) + 
  guides(fill=guide_legend(ncol=2)) + 
  xlab(NULL) + ylab("Relative abundance")

## Fig 2C. Taxa bar plots by sampling date (phylum level)

y1 <- ps.phylum
y2 <- merge_samples(y1, "date") # merge samples by sampling date
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Sample <- factor(y4$Sample, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))
y4$Phylum[y4$Abundance < 0.01] <- "Taxa < 1% abund." #rename genera with < 1% abundance
colourCount = length(unique(y4$Phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
y4$Phylum <- factor(y4$Phylum,levels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Campilobacterota","Cyanobacteria","Deinococcota","Desulfobacterota","Firmicutes","Myxococcota","Proteobacteria","Spirochaetota","Verrucomicrobiota","Taxa < 1% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=1))+ #UCB 2 cols, date 1 col
  xlab(NULL) + ylab("Relative abundance")

## Fig 2D. Taxa bar plots by sampling date (genus level)

y1 <- ps.genus
y2 <- merge_samples(y1, "date") # merge samples by sampling date
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Genus <- as.character(y4$Genus) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Sample <- factor(y4$Sample, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))
y4$Genus[y4$Abundance < 0.03] <- "Taxa < 3% abund." #rename genera with < 3% abundance
colourCount = length(unique(y4$Genus))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
y4$Genus <- factor(y4$Genus,levels=c("Acinetobacter","Aeromonas","bacII-A","betI-A","betIII-A","betVII-A","C39","Dechloromonas","Hydrogenophaga","Malikia","Pseudarcobacter","Pseudomonas","unclassified.Alcaligenaceae","unclassified.alfVI","unclassified.betI","unclassified.Comamonadaceae","unclassified.Methylococcaceae","unclassified.Microbacteriaceae","unclassified.Prevotellaceae","unclassified.Veillonellales-Selenomonadales","Taxa < 3% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
p + geom_bar(aes(fill=Genus), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=1))+ #UCB 2 cols, date 1 col, cluster 2cols
  xlab(NULL) + ylab("Relative abundance")
