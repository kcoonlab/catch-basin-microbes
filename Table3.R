## Table 3. Effects of microbiota on mosquito productivity
## Linear regression analyses with data aggregated by basin

metadata <- read.table("metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
metadata.naUCBrm <- subset(metadata, !date_code=="A")
# cut out all records of UCBs with NAs, not just the records with NAs.
metadata.naUCBrm <- metadata.naUCBrm[!(metadata.naUCBrm$UCB %in% c("43","51","58","79","86","88","32","47","50","64","73","56","61")),]
# use only the dates with all CBs
metadata.naUCBrm.allCBs <- subset(metadata.naUCBrm, date_code %in% c("B","C","F","G","H"))
# collapse by UCB
summary.UCBs <- metadata.naUCBrm.allCBs %>% 
  group_by(UCB) %>%
  dplyr::summarize(Pupae.avg = mean(Pupae, na.rm=TRUE),
                   Pupae.prev = sum(Pupae>=1, na.rm=TRUE)/sum(Pupae>=0, na.rm=TRUE),
                   Moz.total.avg = mean(Moz.total, na.rm=TRUE),
                   Moz.total.prev = sum(Moz.total>=1, na.rm=TRUE)/sum(Moz.total>=0, na.rm=TRUE),
                   Fail.prop = sum(Fail==1, na.rm=TRUE)/sum(Fail>=0,na.rm=TRUE),
                   
                   Larv.restuans.avg = mean(Larv.restuans, na.rm=TRUE),
                   Larv.restuans.prev = sum(Larv.restuans>=1, na.rm=TRUE)/sum(Larv.restuans>=0, na.rm=TRUE),
                   Larv.pipiens.avg = mean(Larv.pipiens, na.rm=TRUE),
                   Larv.pipiens.prev = sum(Larv.pipiens>=1, na.rm=TRUE)/sum(Larv.pipiens>=0, na.rm=TRUE),
                   Larv.japonicus.avg = mean(Larv.japonicus, na.rm=TRUE),
                   Larv.japonicus.prev = sum(Larv.japonicus>=1, na.rm=TRUE)/sum(Larv.japonicus>=0, na.rm=TRUE),
                   Larv.salinarius.avg = mean(Larv.salinarius, na.rm=TRUE),
                   Larv.salinarius.prev = sum(Larv.salinarius>=1, na.rm=TRUE)/sum(Larv.salinarius>=0, na.rm=TRUE),
                   Larv.unID.avg = mean(Larv.unID, na.rm=TRUE),
                   
                   pH.avg = mean(pH, na.rm=TRUE),
                   Temp.avg = mean(Temp, na.rm=TRUE),
                   DO.avg = mean(DO, na.rm=TRUE),
                   Sal.avg = mean(Sal, na.rm=TRUE),
                   Cond.avg = mean(Cond, na.rm=TRUE),
                   
                   shannon.avg = mean(shannon_unrar, na.rm=TRUE),
                   asvs.avg = mean(features_unrar, na.rm=TRUE),
                   C39.avg = mean(C39, na.rm=TRUE),
                   Proteobacteria.avg = mean(Proteobacteria, na.rm=TRUE)
                   # combined_separate = combined_separate
  )

write.csv(summary.UCBs,"summary.UCBs.seq.csv")
summary.UCBs <- read.csv("summary.UCBs.seq.csv", header = TRUE)

summary(lm(Pupae.prev ~ shannon.avg, data=summary.UCBs.all))
summary(lm(Pupae.prev ~ asvs.avg, data=summary.UCBs.all))
summary(lm(Pupae.prev ~ C39.avg, data=summary.UCBs.all))
summary(lm(Pupae.prev ~ Proteobacteria.avg, data=summary.UCBs.all))

summary(lm(Pupae.avg ~ shannon.avg, data=summary.UCBs.all))
summary(lm(Pupae.avg ~ asvs.avg, data=summary.UCBs.all))
summary(lm(Pupae.avg ~ C39.avg, data=summary.UCBs.all))
summary(lm(Pupae.avg ~ Proteobacteria.avg, data=summary.UCBs.all))