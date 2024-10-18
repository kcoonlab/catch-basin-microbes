## Table 1 Linear regressions
## Table 1A data aggregated by basin

cook.longitudinal <- read.table("master_longitudinal.csv", sep=",", header=TRUE)

cook.longitudinal.2021 <- cook.longitudinal[cook.longitudinal$Year=="2021",]

library(dplyr)
cook.basins.2021 <- cook.longitudinal.2021 %>%
  group_by(UCB) %>%
  dplyr::summarize(Pupae.avg = mean(Pupae, na.rm=TRUE),
            Pupae.sd = sd(Pupae, na.rm=TRUE),
            Pupae.prev = sum(Pupae>=1, na.rm=TRUE)/sum(Pupae>=0, na.rm=TRUE),
            Moz.total.avg = mean(Moz.total, na.rm=TRUE),
            Moz.total.sd = sd(Moz.total, na.rm=TRUE),
            Moz.total.prev = sum(Moz.total>=1, na.rm=TRUE)/sum(Moz.total>=0, na.rm=TRUE),
            Fail.prop = sum(Fail==1, na.rm=TRUE)/sum(Fail>=0,na.rm=TRUE),
            
            Larv.restuans.avg = mean(Larv.restuans, na.rm=TRUE),
            Larv.restuans.sd = sd(Larv.restuans, na.rm=TRUE),
            Larv.restuans.prev = sum(Larv.restuans>=1, na.rm=TRUE)/sum(Larv.restuans>=0, na.rm=TRUE),
            Larv.pipiens.avg = mean(Larv.pipiens, na.rm=TRUE),
            Larv.pipiens.sd = sd(Larv.pipiens, na.rm=TRUE),
            Larv.pipiens.prev = sum(Larv.pipiens>=1, na.rm=TRUE)/sum(Larv.pipiens>=0, na.rm=TRUE),
            Larv.japonicus.avg = mean(Larv.japonicus, na.rm=TRUE),
            Larv.japonicus.sd = sd(Larv.japonicus, na.rm=TRUE),
            Larv.japonicus.prev = sum(Larv.japonicus>=1, na.rm=TRUE)/sum(Larv.japonicus>=0, na.rm=TRUE),
            Larv.salinarius.avg = mean(Larv.salinarius, na.rm=TRUE),
            Larv.salinarius.sd = sd(Larv.salinarius, na.rm=TRUE),
            Larv.salinarius.prev = sum(Larv.salinarius>=1, na.rm=TRUE)/sum(Larv.salinarius>=0, na.rm=TRUE),
            Larv.unID.avg = mean(Larv.unID, na.rm=TRUE),
            Larv.unID.sd = sd(Larv.unID, na.rm=TRUE),
            Larv.unID.prev = sum(Larv.unID>=1, na.rm=TRUE)/sum(Larv.unID>=0, na.rm=TRUE),
            
            pH.avg = mean(pH, na.rm=TRUE),
            pH.sd = sd(pH, na.rm=TRUE),
            pH.min = min(pH, na.rm=TRUE),
            pH.max = max(pH, na.rm=TRUE),
            Temp.C.avg = mean(Temp.C, na.rm=TRUE),
            Temp.C.sd = sd(Temp.C, na.rm=TRUE),
            Temp.C.min = min(Temp.C, na.rm=TRUE),
            Temp.C.max = max(Temp.C, na.rm=TRUE),
            DO.avg = mean(DO, na.rm=TRUE),
            DO.sd = sd(DO, na.rm=TRUE),
            DO.min = min(DO, na.rm=TRUE),
            DO.max = max(DO, na.rm=TRUE),
            Sal.avg = mean(Sal, na.rm=TRUE),
            Sal.sd = sd(Sal, na.rm=TRUE),
            Sal.min = min(Sal, na.rm=TRUE),
            Sal.max = max(Sal, na.rm=TRUE),
            Cond.avg = mean(Cond, na.rm=TRUE),
            Cond.sd = sd(Cond, na.rm=TRUE),
            Cond.min = min(Cond, na.rm=TRUE),
            Cond.max = max(Cond, na.rm=TRUE)
            )
head(cook.basins.2021)

write.csv(cook.basins.2021, "summary_basins_2021.csv", row.names=FALSE)


cook.basins.2021 <- read.csv("summary_basins_2021.csv",sep=",", header=TRUE)

cook.basins.2021 <- cook.basins.2021 %>% mutate(pH_bin = if_else(pH.avg <=7, "acid", "base"))
cook.basins.2021 <- cook.basins.2021 %>% mutate(Temp_bin = if_else(Temp.C.avg < 16, "<16C", 
                                           if_else(Temp.C.avg < 20, "16-20C", 
                                                   if_else(Temp.C.avg < 24, "20-24C", ">24C"))))
cook.basins.2021$Temp_bin <- cook.basins.2021$Temp_bin %>% as.factor
cook.basins.2021$combined_separate <- cook.basins.2021$combined_separate %>% as.factor
table(cook.basins.2021$combined_separate, cook.basins.2021$Temp_bin)
cook.basins.2021$Pupae.prev

cook.basins.2021$Flowgroup_coarse <- factor(cook.basins.2021$Flowgroup, levels=c("Gibbons","MayfairCarlyle","Stratford","Donald-Banta","MinerEvanstonRammer"))
cook.basins.2021$combined_separate <- cook.basins.2021$combined_separate %>% as.factor
cook.basins.2021$Flowgroup_coarse


summary(lm(pH.avg ~ Pupae.prev, data = cook.basins.2021))
summary(lm(temp.avg ~ Pupae.prev, data = cook.basins.2021))
summary(lm(Cond.avg ~ Pupae.prev, data = cook.basins.2021))
summary(lm(DO.avg ~ Pupae.prev, data = cook.basins.2021))
summary(lm(Sal.avg ~ Pupae.prev, data = cook.basins.2021))

summary(lm(pH.avg ~ Pupae.avg, data = cook.basins.2021))
summary(lm(temp.avg ~ Pupae.avg, data = cook.basins.2021))
summary(lm(Cond.avg ~ Pupae.avg, data = cook.basins.2021))
summary(lm(DO.avg ~ Pupae.avg, data = cook.basins.2021))
summary(lm(Sal.avg ~ Pupae.avg, data = cook.basins.2021))







## Table 1B data aggregated by sampling date

cook <- read.table("WQ_Dips_2020_2021.csv",sep=",",header=TRUE)
cook["rain"][cook["rain"]=="trace"]<- 0.0001
cook$rain <- as.numeric(cook$rain)
cook <- subset(cook, Year=="2021")
cook$Date<- factor(cook$Date, levels = c("4/15/21","6/4/21","6/11/21","6/18/21","6/25/21","7/2/21","7/9/21","7/16/21","7/23/21","7/30/21","8/6/21","8/13/21","8/20/21","8/27/21","9/3/21","9/10/21","9/17/21","9/24/21"))
cook$rain %>% class
cook["rain"][cook["rain"]=="trace"]<- 0.001
cook$rain <- as.numeric(cook$rain)

cook.day <- cook %>% 
  group_by(Date) %>%
  dplyr::summarize(rainfall = mean(rain, na.rm=TRUE),
                  Pupae.avg = mean(Pupae, na.rm=TRUE),
                   Pupae.total = sum(Pupae, na.rm=TRUE),
                  Pupae.prev = sum(Pupae>=1, na.rm=TRUE)/sum(Pupae>=0, na.rm=TRUE),
                   Moz.avg = mean(Moz.total, na.rm=TRUE),
                   Moz.total = sum(Moz.total, na.rm=TRUE),
                   Moz.prev = sum(Moz.total>=1, na.rm=TRUE)/sum(Moz.total>=0, na.rm=TRUE),
                  Fail.prop = sum(Fail>=1, na.rm=TRUE)/sum(Fail>=0, na.rm=TRUE),
                  pH.avg = mean(pH, na.rm=TRUE),
                  Temp.avg = mean(Temp.C, na.rm=TRUE),
                  Cond.avg = mean(Cond, na.rm=TRUE),
                  DO.avg = mean(DO, na.rm=TRUE),
                  Sal.avg = mean(Sal, na.rm=TRUE)
                  )

summary(lm(pH.avg ~ Pupae.prev, data = cook.day))
summary(lm(temp.avg ~ Pupae.prev, data = cook.day))
summary(lm(Cond.avg ~ Pupae.prev, data = cook.day))
summary(lm(DO.avg ~ Pupae.prev, data = cook.day))
summary(lm(Sal.avg ~ Pupae.prev, data = cook.day))

summary(lm(pH.avg ~ Pupae.avg, data = cook.day))
summary(lm(temp.avg ~ Pupae.avg, data = cook.day))
summary(lm(Cond.avg ~ Pupae.avg, data = cook.day))
summary(lm(DO.avg ~ Pupae.avg, data = cook.day))
summary(lm(Sal.avg ~ Pupae.avg, data = cook.day))


## Table 1C Linear mixed-effects models with unaggregated data (observations only where pupae were present)

metadata <- read.table("metadata.txt", sep="\t", header=TRUE)
metadata$date <- factor(metadata$date, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))
metadata <- subset(metadata, metadata$sample_control=="sample")
metadata <- metadata %>% mutate(pH_bin = if_else(pH <=7, "acid", "base"))
metadata$date <- metadata$date %>% as.factor

metadata.pupaepresent <- subset(metadata, Pupae_pa=="1")

lmer(sqrt(Pupae) ~ pH * Datefactor + (1|UCB), data = metadata.pupaepresent)
lmer(sqrt(Pupae) ~ Temp.C * Datefactor + (1|UCB), data = metadata.pupaepresent)
lmer(sqrt(Pupae) ~ Cond * Datefactor + (1|UCB), data = metadata.pupaepresent)
lmer(sqrt(Pupae) ~ DO * Datefactor + (1|UCB), data = metadata.pupaepresent)
lmer(sqrt(Pupae) ~ Sal * Datefactor + (1|UCB), data = metadata.pupaepresent)













