library("phyloseq")
library("ggplot2")
library("ggsignif")
library("dplyr")
library("tibble")
library("plyr")

## Fig 5. Bacterial community differences by pupal occurrence

## Fig 5A. ASVs by pupae p/a

ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
metadata <- data.frame(sample_data(ps.all.rerooted))
metadata <- subset(metadata, sample_control=="sample")
features.pupaepa <- ggplot(data=subset(metadata,Pupae_pa != "NA"), 
                    aes(x=Pupae_pa, y=features_unrar)) +
                    geom_boxplot(fill="dark grey") +
                    geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE) +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                    scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present")) +
                    xlab("Pupae presence")+ ylab("ASV richness")
features.pupaepa
kruskal.test(features_unrar ~ Pupae_pa, data = metadata)

## Fig 5B. Shannon by pupae p/a

ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
metadata <- data.frame(sample_data(ps.all.rerooted))
metadata <- subset(metadata, sample_control=="sample")
shannon.pupaepa <- ggplot(data=subset(metadata,Pupae_pa != "NA"), 
                          aes(x=Pupae_pa, y=shannon_unrar)) +
                          geom_boxplot(fill="dark grey") +
                          geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE)+
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                          scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present"))+
                          xlab("Pupae presence")+ ylab("Shannon index")
shannon.pupaepa
kruskal.test(shannon_unrar ~ Pupae_pa, data = metadata)

***## Fig 5C. Proteobacteria relative abundance by pupae p/a

ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
metadata <- data.frame(sample_data(ps.all.rerooted))
metadata <- subset(metadata, sample_control=="sample")
Proteobacteria.pupaepa <- ggplot(data=subset(metadata,Pupae_pa != "NA"), 
                          aes(x=Pupae_pa, y=Proteobacteria)) +
                          geom_boxplot(fill="dark grey") +
                          geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE) +
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                          scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present"))+
                          xlab("Pupae presence")+ ylab("Relative abundance Proteobacteria")
Proteobacteria.pupaepa
kruskal.test(Proteobacteria ~ Pupae_pa, data = metadata)

***## Fig 5D. C39 relative abundance by pupae p/a

library(ggsignif)
C39.pupaepa <- ggplot(data=subset(metadata_div,Pupae_pa != "NA"), 
                      aes(x=Pupae_pa, y=C39)) +
  geom_boxplot(fill="dark grey") +
  geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present"))+
  xlab("Pupae presence")+ ylab("Relative abundance C39")
C39.pupaepa

***## Fig 5E. Pupae per dip by relative abundance of C39

ggplot(metadata.naUCBrm.allCBs, aes(x=C39, y=sqrt(Pupae), color=date_code))+
  geom_point()
  
## Fig 5F. Relative abundance of C39 strains over time

ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
ps.all.rerooted.rel <- transform_sample_counts(ps.all.rerooted, function(x) x / sum(x) )
table.C39 <- subset_taxa(ps.all.rerooted.rel, Genus == "C39") # 11 ASVs
tax_table(table.C39)
# 1c896ad60b318aea6eb67595fe68d2bf
# 56f81ede4d385001b8135ac61458befe
# 713852f86aa38fb1f6f5291dee0db276
# 3813a6b872a29e05f211916bfb3397ef
# 00a1a1e7419bd374f7b8993902555db6
# 34e720bc46aa554b53222a0915c6410e
# 37087046f6f5e933bf0dd42420f81cfe
# 912fefd634a7c68b38d929afd0f1cb68
# 9ae8011e9347ce18b60ba481bfb72e32
# a5f2920be80527d6ec323de39babdf2e
# 5ecd2e447fe9cfbef560d1c5e0738bb1

df.asvtable <- as.data.frame(otu_table(ps.all.rerooted.rel))
df.asvtable$asvid <- rownames(df.asvtable)
df.C39.asvs <- subset(df.asvtable, df.asvtable$asvid %in% c("1c896ad60b318aea6eb67595fe68d2bf","56f81ede4d385001b8135ac61458befe",
                                                            "713852f86aa38fb1f6f5291dee0db276","3813a6b872a29e05f211916bfb3397ef",
                                                            "00a1a1e7419bd374f7b8993902555db6","34e720bc46aa554b53222a0915c6410e",
                                                            "37087046f6f5e933bf0dd42420f81cfe","912fefd634a7c68b38d929afd0f1cb68",
                                                            "9ae8011e9347ce18b60ba481bfb72e32","a5f2920be80527d6ec323de39babdf2e",
                                                            "5ecd2e447fe9cfbef560d1c5e0738bb1"))
relabund.C39.asvs <-as.data.frame(t(df.C39.asvs))
relabund.C39.asvs <- relabund.C39.asvs %>% rownames_to_column(var="sampleid")
relabund.C39.asvs %>% head
metadata.add <- data.frame(sample_data(ps.all.rerooted.rel))
metadata.add <- metadata.add %>% rownames_to_column(var="sampleid")
metadata.relabund <- join_all(list(relabund.C39.asvs, metadata.add), by='sampleid',type='left')
metadata.relabund <- subset(metadata.relabund, metadata.relabund$sampleid != "asvid")
#metadata.relabund$`1c896ad60b318aea6eb67595fe68d2bf` <- metadata.relabund$`X1c896ad60b318aea6eb67595fe68d2bf` %>% as.numeric
#metadata.relabund$`56f81ede4d385001b8135ac61458befe` <- metadata.relabund$`X56f81ede4d385001b8135ac61458befe` %>% as.numeric
#metadata.relabund$`713852f86aa38fb1f6f5291dee0db276` <- metadata.relabund$`X713852f86aa38fb1f6f5291dee0db276` %>% as.numeric
#metadata.relabund$`3813a6b872a29e05f211916bfb3397ef` <- metadata.relabund$`X3813a6b872a29e05f211916bfb3397ef` %>% as.numeric
#metadata.relabund$`00a1a1e7419bd374f7b8993902555db6` <- metadata.relabund$`X00a1a1e7419bd374f7b8993902555db6` %>% as.numeric
#metadata.relabund$`34e720bc46aa554b53222a0915c6410e` <- metadata.relabund$`X34e720bc46aa554b53222a0915c6410e` %>% as.numeric
#metadata.relabund$`37087046f6f5e933bf0dd42420f81cfe` <- metadata.relabund$`X37087046f6f5e933bf0dd42420f81cfe` %>% as.numeric
#metadata.relabund$`912fefd634a7c68b38d929afd0f1cb68` <- metadata.relabund$`X912fefd634a7c68b38d929afd0f1cb68` %>% as.numeric
#metadata.relabund$`9ae8011e9347ce18b60ba481bfb72e32` <- metadata.relabund$`X9ae8011e9347ce18b60ba481bfb72e32` %>% as.numeric
#metadata.relabund$`a5f2920be80527d6ec323de39babdf2e` <- metadata.relabund$`a5f2920be80527d6ec323de39babdf2e` %>% as.numeric
#metadata.relabund$`5ecd2e447fe9cfbef560d1c5e0738bb1` <- metadata.relabund$`X5ecd2e447fe9cfbef560d1c5e0738bb1` %>% as.numeric
#metadata.relabund$`00a1a1e7419bd374f7b8993902555db6` %>% class

metadata.relabund$date <- factor(metadata.relabund$date, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))


library(reshape2)
metadata.relabund.long <- reshape2::melt(data=metadata.relabund,
                                 id.vars=c("UCB","date"),
                                 measure.vars=c("37087046f6f5e933bf0dd42420f81cfe","34e720bc46aa554b53222a0915c6410e",
                                                "5ecd2e447fe9cfbef560d1c5e0738bb1","a5f2920be80527d6ec323de39babdf2e",
                                                "3813a6b872a29e05f211916bfb3397ef","00a1a1e7419bd374f7b8993902555db6",
                                                "912fefd634a7c68b38d929afd0f1cb68","1c896ad60b318aea6eb67595fe68d2bf",
                                                "713852f86aa38fb1f6f5291dee0db276","56f81ede4d385001b8135ac61458befe",
                                                "9ae8011e9347ce18b60ba481bfb72e32"),
                                 variable.name="asvid",
                                 value.name="relabund") 
ggplot(metadata.relabund.long, aes(x=date, y=relabund))+
  geom_boxplot() + 
  facet_wrap(vars(asvid), ncol=2)

# subset(metadata.relabund.long, asvid=="9ae8011e9347ce18b60ba481bfb72e32")
ggplot(metadata.relabund.long, aes(x=date, y=relabund, fill=asvid))+
  geom_bar(stat="identity",position=position_dodge())+
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.key = element_rect(fill="white"))+
  scale_x_discrete(breaks=c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"),
                   labels=c("4/15","6/11","6/25","7/9","7/23","8/6","8/27","9/17"))+
  # scale_fill_discrete(breaks=c("1c896ad60b318aea6eb67595fe68d2bf","56f81ede4d385001b8135ac61458befe",
                              # "713852f86aa38fb1f6f5291dee0db276","3813a6b872a29e05f211916bfb3397ef",
                              # "00a1a1e7419bd374f7b8993902555db6","34e720bc46aa554b53222a0915c6410e",
                              # "37087046f6f5e933bf0dd42420f81cfe","912fefd634a7c68b38d929afd0f1cb68",
                              # "9ae8011e9347ce18b60ba481bfb72e32","a5f2920be80527d6ec323de39babdf2e",
                              # "5ecd2e447fe9cfbef560d1c5e0738bb1"), 
                     # labels=c("strain 1","strain 2","strain 3","strain 4","strain 5","strain 6",
                              # "strain 7","strain 8","strain 9","strain 10","strain 11"))+
  ylab("Relative abundance (avg among basins)")+ xlab(NULL)+ labs(fill="ASV ID")

c39_date <- ggplot(metadata.relabund.long, aes(x=date, y=relabund, color=asvid))+
  geom_boxplot()+ #show.legend=FALSE
  # geom_segment(data=metadata.relabund.long, aes(x=, xend=xmax,y=middle, yend=middle), colour="red", size=2)+
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.key = element_rect(fill="white"), legend.text = element_text(size = 6))+
  scale_x_discrete(breaks=c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"),
                   labels=c("4/15","6/11","6/25","7/9","7/23","8/6","8/27","9/17"))+
  ylab("Relative abundance")+ xlab(NULL)+ labs(color="ASV ID")
c39_date

dat <- ggplot_build(c39_date)$data[[1]]
dat
c39_date+ geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                         y=middle, yend=middle), colour="red", size=2)
boxplot(data=metadata.relabund.long, medcol="red",)








