## Table 3. Effects of microbiota on mosquito productivity

library(dplyr)

metadata <- read.table("input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
metadata.byCB <- metadata %>% 
  group_by(UCB) %>%
  dplyr::summarize(Pupae.prev = sum(Pupae>=1, na.rm=TRUE)/sum(Pupae>=0, na.rm=TRUE),
  				   Pupae.avg = mean(Pupae, na.rm=TRUE),
  				   shannon_unrar.avg = mean(shannon_unrar, na.rm=TRUE),
  				   richness.avg = mean(features_unrar, na.rm=TRUE),
  				   C39.avg = mean(C39, na.rm=TRUE),
  				   Proteo.avg = mean(Proteobacteria, na.rm=TRUE))

summary(lm(Pupae.prev ~ shannon_unrar.avg, data = metadata.byCB))
summary(lm(Pupae.prev ~ richness.avg, data = metadata.byCB))
summary(lm(Pupae.prev ~ C39.avg, data = metadata.byCB))
summary(lm(Pupae.prev ~ Proteo.avg, data = metadata.byCB))

summary(lm(Pupae.avg ~ shannon_unrar.avg, data = metadata.byCB))
summary(lm(Pupae.avg ~ richness.avg, data = metadata.byCB))
summary(lm(Pupae.avg ~ C39.avg, data = metadata.byCB))
summary(lm(Pupae.avg ~ Proteo.avg, data = metadata.byCB))
