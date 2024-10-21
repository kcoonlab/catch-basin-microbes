## Fig 3. Catch basin microbiota biotypes identified by PAM clustering

library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(ggsignif)
library(tibble)
library(plyr)

## Fig 3A. Taxa bar plots by biotype (phylum level)

ps.final <- readRDS("input-files/ps.final.rds")
ps.phylum <- tax_glom(ps.all, "Phylum", NArm = TRUE)
y1 <- ps.phylum
sample_data(y1)$cluster_philr[which(sample_data(y1)$cluster_philr == 1)] <- "A"
sample_data(y1)$cluster_philr[which(sample_data(y1)$cluster_philr == 2)] <- "B"
y2 <- merge_samples(y1, "cluster_philr") # merge samples by biotype
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Phylum[y4$Abundance < 0.01] <- "Taxa < 1% abund." #rename genera with < 1% abundance
colourCount = length(unique(y4$Phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
y4$Phylum <- factor(y4$Phylum,levels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Campilobacterota","Cyanobacteria","Deinococcota","Desulfobacterota","Firmicutes","Myxococcota","Proteobacteria","Spirochaetota","Verrucomicrobiota","Taxa < 1% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=2)) +
  xlab(NULL) + ylab("Relative abundance")

## Fig 3B. Taxa bar plots by biotype (genus level)
                              
ps.genus <- tax_glom(ps.all, "Genus", NArm = TRUE)
y1 <- ps.genus
sample_data(y1)$cluster_philr[which(sample_data(y1)$cluster_philr == 1)] <- "A"
sample_data(y1)$cluster_philr[which(sample_data(y1)$cluster_philr == 2)] <- "B"
y2 = merge_samples(y1, "cluster_philr") # merge samples by biotype
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Genus <- as.character(y4$Genus) #convert to character
y4$Sample <- as.factor(y4$Sample)
y4$Genus[y4$Abundance < 0.01] <- "Taxa < 1% abund." #rename genera with < 1% abundance
colourCount = length(unique(y4$Genus))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
y4$Genus <- factor(y4$Genus, levels=c("Acinetobacter","Aeromonas","alfIV-A","bacII-A","Bacteroides","betI-A","betIII-A","betVII-A","C39","Cloacibacterium","Dechloromonas","Hydrogenophaga","Malikia","Pseudarcobacter","Pseudomonas","Tolumonas","unclassified.Alcaligenaceae","unclassified.alfVI","unclassified.alfVII","unclassified.betI","unclassified.Comamonadaceae","unclassified.Enterobacterales","unclassified.Enterobacteriaceae","unclassified.Lachnospiraceae","unclassified.Methylococcaceae","unclassified.Microbacteriaceae","unclassified.Prevotellaceae","unclassified.Rhodocyclaceae","unclassified.Veillonellales-Selenomonadales","Taxa < 1% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
p + geom_bar(aes(fill=Genus), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=2)) +
  xlab(NULL) + ylab("Relative abundance")
                
## Fig 3C. ASVs by biotype

metadata <- data.frame(sample_data(ps.final))
metadata <- subset(metadata, sample_control == "sample")
metadata$Biotype <- as.factor(metadata$Biotype)
cluster_features <- ggplot(data=na.omit(metadata[,c("Biotype","ASV.richness")]), aes(x=Biotype, y=ASV.richness)) + 
  geom_boxplot(show.legend=FALSE, fill="dark grey")+
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  theme(axis.title = element_text(size=13), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(size=13), legend.position="blank")+
  ylab("ASV richness")+ xlab(NULL)+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B"))
cluster_features
kruskal.test(metadata$ASV.richness, metadata$Biotype)                              

## Fig 3D. Shannon index by biotype

cluster_shannon <-ggplot(data=na.omit(metadata[,c("Biotype","Shannon")]), aes(x=Biotype, y=Shannon)) + 
  geom_boxplot(show.legend=FALSE, fill="dark grey")+
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  theme(axis.title = element_text(size=13), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(size=13), legend.position="blank")+
  ylab("Shannon index")+ xlab(NULL)+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B"))
cluster_shannon
kruskal.test(metadata$Shannon, metadata$Biotype)

## Fig 3E. C39 relative abundance by biotype

ps.final.rel  = transform_sample_counts(ps.final, function(x) x / sum(x) )
ps.final.rel.C39 = subset_taxa(ps.final.rel, Genus=="C39")
metadata.add <- data.frame(sample_data(ps.final.rel))
metadata.add <- metadata.add %>% rownames_to_column(var="sampleid")
relabund.sums.C39 <- as.data.frame(sample_sums(ps.final.rel.C39))
relabund.sums.C39 <- relabund.sums.C39 %>% rownames_to_column(var="sampleid")
relabund.sums.metadata <- join_all(list(relabund.sums.C39, metadata.add), by='sampleid',type='left')
names(relabund.sums.metadata)[names(relabund.sums.metadata) == 'sample_sums(ps.final.rel.C39)'] <- 'C39'
relabund.sums.metadata$Biotype <- as.factor(relabund.sums.metadata$Biotype)
myvars <- c("C39.relabund", "Biotype")
relabund.sums.metadata <- relabund.sums.metadata[myvars]
C39_cluster <- ggplot(data=relabund.sums.metadata, aes(x=Biotype, y=C39.relabund)) +
  geom_boxplot(show.legend=FALSE, fill="dark grey") +
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  labs(y = "Relative abundance C39")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B")) + xlab(NULL)
C39_cluster
kruskal.test(relabund.sums.metadata$C39.relabund, relabund.sums.metadata$Biotype)

## Fig 3F. Firmicutes relative abundance by biotype

ps.final.rel.Firmicutes = subset_taxa(ps.final.rel, Phylum=="Firmicutes")
metadata.add <- data.frame(sample_data(ps.final.rel))
metadata.add <- metadata.add %>% rownames_to_column(var="sampleid")
relabund.sums.Firmicutes <- as.data.frame(sample_sums(ps.final.rel.Firmicutes))
relabund.sums.Firmicutes <- relabund.sums.Firmicutes %>% rownames_to_column(var="sampleid")
relabund.sums.metadata <- join_all(list(relabund.sums.Firmicutes, metadata.add), by='sampleid',type='left')
names(relabund.sums.metadata)[names(relabund.sums.metadata) == 'sample_sums(ps.final.Firmicutes)'] <- 'Firmicutes'
relabund.sums.metadata$Biotype <- as.factor(relabund.sums.metadata$Biotype)
myvars <- c("Firmicutes.relabund", "Biotype")
relabund.sums.metadata <- relabund.sums.metadata[myvars]
Firmicutes_cluster <- ggplot(data=relabund.sums.metadata, aes(x=Biotype, y=Firmicutes.relabund)) +
  geom_boxplot(show.legend=FALSE, fill="dark grey") +
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  labs(y = "Relative abundance Firmicutes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B")) + xlab(NULL)
Firmicutes_cluster
kruskal.test(relabund.sums.metadata$Firmicutes.relabund, relabund.sums.metadata$Biotype)
