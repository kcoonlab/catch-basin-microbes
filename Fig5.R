## Fig 5. Bacterial community differences by pupal occurrence

set.seed(123)

library(phyloseq)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(tibble)
library(plyr)
library(reshape2)
library(nlme)

## Fig 5A. ASVs by pupae presence/absence

ps.final <- readRDS("input-files/ps.final.rds")
metadata <- data.frame(sample_data(ps.final))
metadata <- subset(metadata, sample_control=="sample")
metadata$Pupae.pres <- as.factor(metadata$Pupae.pres)
features.pupaepa <- ggplot(data=subset(metadata,Pupae.pres != "NA"), 
                    aes(x=Pupae.pres, y=ASV.richness)) +
                    geom_boxplot(fill="dark grey") +
                    geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE) +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                    scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present")) +
                    xlab("Pupae presence")+ ylab("ASV richness")
features.pupaepa

metadata$Datefactor <- metadata$Sampling.date %>% as.factor
metadata <- subset(metadata, !is.na(Pupae.pres)==TRUE)
metadata = metadata %>%
     arrange(Basin.id, Datefactor)                                       
model = lme(ASV.richness ~ Pupae.pres,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata)
anova(model) 

## Fig 5B. Shannon index by pupae presence/absence

shannon.pupaepa <- ggplot(data=subset(metadata,Pupae.pres != "NA"), 
                          aes(x=Pupae.pres, y=Shannon)) +
                          geom_boxplot(fill="dark grey") +
                          geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE) +
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                          scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present"))+
                          xlab("Pupae presence")+ ylab("Shannon index")
shannon.pupaepa

metadata$Datefactor <- metadata$Sampling.date %>% as.factor
metadata <- subset(metadata, !is.na(Pupae.pres)==TRUE)
metadata = metadata %>%
     arrange(Basin.id, Datefactor)                                       
model = lme(Shannon ~ Pupae.pres,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata)
anova(model) 

## Fig 5C. Proteobacteria relative abundance by pupae presence/absence

Proteobacteria.pupaepa <- ggplot(data=subset(metadata,Pupae.pres != "NA"), 
                          aes(x=Pupae.pres, y=Proteobacteria.relabund)) +
                          geom_boxplot(fill="dark grey") +
                          geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE) +
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                          scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present"))+
                          xlab("Pupae presence")+ ylab("Relative abundance Proteobacteria")
Proteobacteria.pupaepa

metadata$Datefactor <- metadata$Sampling.date %>% as.factor
metadata <- subset(metadata, !is.na(Pupae.pres)==TRUE)
metadata = metadata %>%
     arrange(Basin.id, Datefactor)                                       
model = lme(Proteobacteria.relabund ~ Pupae.pres,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata)
anova(model) 

## Fig 5D. C39 relative abundance by pupae presence/absence

C39.pupaepa <- ggplot(data=subset(metadata,Pupae.pres != "NA"), 
                      aes(x=Pupae.pres, y=C39.relabund)) +
                      geom_boxplot(fill="dark grey") +
                      geom_signif(comparisons=list(c("0","1")), map_signif_level=TRUE) +
                      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                      scale_x_discrete(breaks=c("0","1"), labels=c("Absent", "Present"))+
                      xlab("Pupae presence")+ ylab("Relative abundance C39")
C39.pupaepa

metadata$Datefactor <- metadata$Sampling.date %>% as.factor
metadata <- subset(metadata, !is.na(Pupae.pres)==TRUE)
metadata = metadata %>%
     arrange(Basin.id, Datefactor)                                       
model = lme(C39.relabund ~ Pupae.pres,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata)
anova(model) 

## Fig 5E. Pupae per dip by relative abundance of C39

ggplot(metadata, aes(x=C39.relabund, y=sqrt(Pupae.abund))) +
  geom_point(color = "black", fill = "black")
  
## Fig 5F. Relative abundance of C39 strains over time

ps.final.rel <- transform_sample_counts(ps.final, function(x) x / sum(x) )
table.C39 <- subset_taxa(ps.final.rel, Genus == "C39") 
tax_table(table.C39) # 11 ASVs
#1c896ad60b318aea6eb67595fe68d2bf
#56f81ede4d385001b8135ac61458befe
#713852f86aa38fb1f6f5291dee0db276
#3813a6b872a29e05f211916bfb3397ef
#00a1a1e7419bd374f7b8993902555db6
#34e720bc46aa554b53222a0915c6410e
#37087046f6f5e933bf0dd42420f81cfe
#912fefd634a7c68b38d929afd0f1cb68
#9ae8011e9347ce18b60ba481bfb72e32
#a5f2920be80527d6ec323de39babdf2e
#5ecd2e447fe9cfbef560d1c5e0738bb1

df.asvtable <- as.data.frame(otu_table(ps.final.rel))
df.asvtable$asvid <- rownames(df.asvtable)
df.C39.asvs <- subset(df.asvtable, df.asvtable$asvid %in% c("1c896ad60b318aea6eb67595fe68d2bf","56f81ede4d385001b8135ac61458befe",
                                                            "713852f86aa38fb1f6f5291dee0db276","3813a6b872a29e05f211916bfb3397ef",
                                                            "00a1a1e7419bd374f7b8993902555db6","34e720bc46aa554b53222a0915c6410e",
                                                            "37087046f6f5e933bf0dd42420f81cfe","912fefd634a7c68b38d929afd0f1cb68",
                                                            "9ae8011e9347ce18b60ba481bfb72e32","a5f2920be80527d6ec323de39babdf2e",
                                                            "5ecd2e447fe9cfbef560d1c5e0738bb1"))
relabund.C39.asvs <-as.data.frame(t(df.C39.asvs))
relabund.C39.asvs <- relabund.C39.asvs %>% rownames_to_column(var="sampleid")
metadata.add <- data.frame(sample_data(ps.final.rel))
metadata.add <- metadata.add %>% rownames_to_column(var="sampleid")
metadata.relabund <- join(relabund.C39.asvs, metadata.add, by='sampleid',type='left')
metadata.relabund <- subset(metadata.relabund, metadata.relabund$sampleid != "asvid")

metadata.relabund.long <- reshape2::melt(data=metadata.relabund,
                                 id.vars=c("Basin.id","Sampling.date"),
                                 measure.vars=c("37087046f6f5e933bf0dd42420f81cfe","34e720bc46aa554b53222a0915c6410e",
                                                "5ecd2e447fe9cfbef560d1c5e0738bb1","a5f2920be80527d6ec323de39babdf2e",
                                                "3813a6b872a29e05f211916bfb3397ef","00a1a1e7419bd374f7b8993902555db6",
                                                "912fefd634a7c68b38d929afd0f1cb68","1c896ad60b318aea6eb67595fe68d2bf",
                                                "713852f86aa38fb1f6f5291dee0db276","56f81ede4d385001b8135ac61458befe",
                                                "9ae8011e9347ce18b60ba481bfb72e32"),
                                 variable.name="asvid",
                                 value.name="relabund") 

metadata.relabund.long$date <- factor(metadata.relabund.long$Sampling.date, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))
metadata.relabund.long$relabund <- as.numeric(metadata.relabund.long$relabund)                               

c39_date <- ggplot(metadata.relabund.long, aes(x=date, y=relabund, color=factor(asvid))) +
  geom_boxplot() +
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.key = element_rect(fill="white"), legend.text = element_text(size = 6)) +
  scale_x_discrete(breaks=c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"),
                   labels=c("4/15","6/11","6/25","7/9","7/23","8/6","8/27","9/17")) +
  ylab("Relative abundance")+ xlab(NULL)+ labs(color="ASV ID")
c39_date
