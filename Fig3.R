## Fig 3. Catch basin microbiota biotypes identified by PAM clustering



# see Fig 2 for taxa bar plot input

sample_data(ps.input)$cluster_philr_genus <- as.factor(sample_data(ps.input)$cluster_philr_genus)
levels(sample_data(ps.input)$cluster_philr_genus) <- c("A","B")
sample_data(ps.input)$cluster_philr_family <- as.factor(sample_data(ps.input)$cluster_philr_family)
levels(sample_data(ps.input)$cluster_philr_family) <- c("A","B")
sample_data(ps.input)$cluster_philr_phylum <- as.factor(sample_data(ps.input)$cluster_philr_phylum)
levels(sample_data(ps.input)$cluster_philr_phylum) <- c("A","B","C","D","E")

sample_data(ps.phylum)$cluster_philr <- as.factor(sample_data(ps.phylum)$cluster_philr)
levels(sample_data(ps.phylum)$cluster_philr) <- c("A","B")
# sample_data(ps.family)$cluster_philr <- as.factor(sample_data(ps.family)$cluster_philr)
# levels(sample_data(ps.family)$cluster_philr) <- c("A","B")
sample_data(ps.genus)$cluster_philr <- as.factor(sample_data(ps.genus)$cluster_philr)
levels(sample_data(ps.genus)$cluster_philr) <- c("A","B")
sample_data(ps.phylum)$cluster_philr


sample_data(ps.input)$date %>% class
date <- get_variable(ps.input,"date_code")
combsep <- get_variable(ps.input, "combined_separate")
sample_data(ps.input)$NewPastedVar <- mapply(paste0, date, combsep, collapse="_")
sample_data(ps.input)$NewPastedVar <- sample_data(ps.input)$NewPastedVar %>% as.factor #%>% as.numeric
sample_data(ps.input)$cluster_philr <- as.factor(sample_data(ps.input)$cluster_philr)
sample_data(ps.input)$UCB <- sample_data(ps.input)$UCB %>% as.factor
sample_data(ps.input)$combined_separate <- sample_data(ps.input)$combined_separate %>% as.factor


# Fig 3A. Taxa bar plots by biotype, phyla
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
# by UCB
# y4$Phylum <- factor(y4$Phylum,levels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Campilobacterota","Cyanobacteria","Deinococcota","Desulfobacterota","Firmicutes","Fusobacteriota","Proteobacteria","Spirochaetota","unclassified.Bacteria","Verrucomicrobiota","Taxa < 1% abund."))
# by date
y4$Phylum <- factor(y4$Phylum,levels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Campilobacterota","Cyanobacteria","Deinococcota","Desulfobacterota","Firmicutes","Myxococcota","Proteobacteria","Spirochaetota","Verrucomicrobiota","Taxa < 1% abund."))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
#date, cluster
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=1))+ #UCB 2 cols, date 1 col
  xlab(NULL) + ylab("Relative abundance")


# Fig 3B. Taxa bar plots by biotype, genera
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
# by cluster
# y4$Genus <- factor(y4$Genus, levels=c("Acinetobacter","Aeromonas","alfIV-A","bacII-A","Bacteroides","betI-A","betIII-A","betVII-A","C39","Cloacibacterium","Dechloromonas","Hydrogenophaga","Malikia","Pseudarcobacter","Pseudomonas","Tolumonas","unclassified.Alcaligenaceae","unclassified.alfVI","unclassified.alfVII","unclassified.betI","unclassified.Comamonadaceae","unclassified.Enterobacterales","unclassified.Enterobacteriaceae","unclassified.Lachnospiraceae","unclassified.Methylococcaceae","unclassified.Microbacteriaceae","unclassified.Prevotellaceae","unclassified.Rhodocyclaceae","unclassified.Veillonellales-Selenomonadales","Taxa < 1% abund."))
p <- ggplot(data=y4, aes(x=reorder(Sample,combined_separate), y=Abundance))
p <- ggplot(data=y4, aes(x=Sample, y=Abundance))
# date, cluster
p + geom_bar(aes(fill=Genus), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + 
  theme(panel.background = element_blank(), legend.position="right", axis.ticks.x=element_blank()) + 
  guides(fill=guide_legend(ncol=1))+ #UCB 2 cols, date 1 col, cluster 2cols
  xlab(NULL) + ylab("Relative abundance")






# Figs 3C-F

# read in phyloseq object, convert vars to factor as needed
ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
metadata <- data.frame(sample_data(ps.all.rerooted))
metadata$cluster_philr_genus <- as.factor(metadata$cluster_philr_genus) #cluster_philr_family, cluster_philr_phylum
metadata <- subset(metadata, sample_control =="sample")
metadata$cluster_philr <- as.factor(metadata$cluster_philr)


## Fig 3C. ASVs by biotype
# ASV richness by philr cluster
cluster_features <-ggplot(data=na.omit(metadata[,c("cluster_philr","features_unrar")]), aes(x=cluster_philr, y=features_unrar)) + 
  geom_boxplot(show.legend=FALSE, fill="dark grey")+
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  theme(axis.title = element_text(size=13), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(size=13), legend.position="blank")+
  ylab("ASV richness")+ xlab(NULL)+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B"))
cluster_features
kruskal.test(metadata$features_unrar, metadata$cluster_philr)

## Fig 3D. Shannon by biotype
# shannon diversity by philr cluster
cluster_shannon <-ggplot(data=na.omit(metadata[,c("cluster_philr","shannon_unrar")]), aes(x=cluster_philr, y=shannon_unrar)) + 
  geom_boxplot(show.legend=FALSE, fill="dark grey")+
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  theme(axis.title = element_text(size=13), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(size=13), legend.position="blank")+
  ylab("Shannon index")+ xlab(NULL)+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B"))
cluster_shannon
kruskal.test(metadata$shannon_unrar, metadata$cluster_philr)





## Figs 3E-F
ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
ps.all.rerooted.rel  = transform_sample_counts(ps.all.rerooted, function(x) x / sum(x) )

ps.all.tax10 = subset_taxa(ps.all.rerooted.rel, Phylum=="Patescibacteria")
metadata.add <- data.frame(sample_data(ps.all.rerooted.rel))
metadata.add <- metadata.add %>% rownames_to_column(var="sampleid")
# head(otu_table(ps.all.tax10))
relabund.sums.tax10 <- as.data.frame(sample_sums(ps.all.tax10))
relabund.sums.tax10 <- relabund.sums.tax10 %>% rownames_to_column(var="sampleid")
relabund.sums.metadata <- join_all(list(relabund.sums.tax10, metadata.add), by='sampleid',type='left')
names(relabund.sums.metadata)[names(relabund.sums.metadata) == 'sample_sums(ps.all.tax4)'] <- 'C39'
# names(relabund.sums.metadata)[names(relabund.sums.metadata) == 'sample_sums(ps.all.tax8)'] <- 'Firmicutes'

saveRDS(relabund.sums.metadata, "relabund_sums_metadata_phylumclust.rds")
saveRDS(relabund.sums.metadata, "relabund_sums_metadata_genusclust.rds")
relabund.sums.metadata.phylum <- readRDS("relabund_sums_metadata_phylumclust.rds")
relabund.sums.metadata.genus <- readRDS("relabund_sums_metadata_genusclust.rds")
relabund.sums.metadata.genus$cluster_philr <- as.factor(relabund.sums.metadata.genus$cluster_philr)
relabund.sums.metadata.phylum$cluster_philr <- as.factor(relabund.sums.metadata.phylum$cluster_philr)
write.csv(relabund.sums.metadata.phylum, "relabund_sums_metadata_phylumclust.csv")
relabund.sums.metadata.phylum <- read.csv("relabund_sums_metadata_phylumclust.csv", header=TRUE)
write.csv(relabund.sums.metadata.genus, "relabund_sums_metadata_genusclust.csv")
relabund.sums.metadata.genus <- read.csv("relabund_sums_metadata_genusclust.csv", header=TRUE)
relabund.sums.metadata.genus$cluster_philr <- relabund.sums.metadata.genus$cluster_philr %>% as.factor




## Fig 3E. C39 by biotype
# relative abundance of C39 by philr cluster
c39_cluster <- ggplot(data=relabund.sums.metadata.genus, aes(x=cluster_philr, y=C39)) + #x=combined_separate
  geom_boxplot(show.legend=FALSE, fill="dark grey") +
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  # geom_dotplot(binaxis='y', stackdir='center', binpositions="all", stackgroups=TRUE, aes()) +
  labs(y = "Relative abundance C39")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B")) + xlab(NULL)
c39_cluster
kruskal.test(relabund.sums.metadata.genus$C39, relabund.sums.metadata.genus$cluster_philr)


## Fig 3F Firmicutes by biotype
# plot relative abundance of Firmicutes by philr cluster
metadata <- subset(metadata, sample_control=="sample")
metadata <- subset(metadata, Firmicutes!="NA")
firmicutes_cluster <-ggplot(data=metadata, aes(x=cluster_philr, y=Firmicutes)) + #x=combined_separate
  geom_boxplot(show.legend=FALSE, fill="dark grey") +
  geom_signif(comparisons=list(c("1","2")), map_signif_level=TRUE)+
  labs(y = "Relative abundance Firmicutes")+ xlab(NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B")) + xlab(NULL)
firmicutes_cluster
kruskal.test(metadata$Firmicutes, metadata$cluster_philr)