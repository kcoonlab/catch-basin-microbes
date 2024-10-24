set.seed(123)

## Table 3. Effects of microbiota on mosquito productivity

library(dplyr)

metadata <- read.table("catch-basin-microbes/input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
metadata.byCB <- metadata %>% 
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
  				   Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
  				   Shannon.avg = mean(Shannon, na.rm=TRUE),
  				   ASV.richness.avg = mean(ASV.richness, na.rm=TRUE),
  				   C39.relabund.avg = mean(C39.relabund, na.rm=TRUE),
  				   Proteobacteria.relabund.avg = mean(Proteobacteria.relabund, na.rm=TRUE))

summary(lm(Pupae.prev ~ Shannon.avg, data = metadata.byCB))
summary(lm(Pupae.prev ~ ASV.richness.avg, data = metadata.byCB))
summary(lm(Pupae.prev ~ C39.relabund.avg, data = metadata.byCB))
summary(lm(Pupae.prev ~ Proteobacteria.relabund.avg, data = metadata.byCB))

summary(lm(Pupae.abund.avg ~ Shannon.avg, data = metadata.byCB))
summary(lm(Pupae.abund.avg ~ ASV.richness.avg, data = metadata.byCB))
summary(lm(Pupae.abund.avg ~ C39.relabund.avg, data = metadata.byCB))
summary(lm(Pupae.abund.avg ~ Proteobacteria.relabund.avg, data = metadata.byCB))
