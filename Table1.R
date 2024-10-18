## Table 1 Linear regressions

## Table 1A data aggregated by basin

cook.basins.2021 <- read.csv("summary_basins_2021.csv",sep=",", header=TRUE)
summary(lm(Pupae.prev ~ pH.avg, data = cook.basins.2021))
summary(lm(Pupae.prev ~ Temp.C.avg, data = cook.basins.2021))
summary(lm(Pupae.prev ~ Cond.avg, data = cook.basins.2021))
summary(lm(Pupae.prev ~ DO.avg, data = cook.basins.2021))
summary(lm(Pupae.prev ~ Sal.avg, data = cook.basins.2021))

summary(lm(Pupae.avg ~ pH.avg, data = cook.basins.2021))
summary(lm(Pupae.avg ~ Temp.C.avg, data = cook.basins.2021))
summary(lm(Pupae.avg ~ Cond.avg, data = cook.basins.2021))
summary(lm(Pupae.avg ~ DO.avg, data = cook.basins.2021))
summary(lm(Pupae.avg ~ Sal.avg, data = cook.basins.2021))

***## Table 1B data aggregated by sampling date

cook <- read.table("WQ_Dips_2020_2021.csv",sep=",",header=TRUE)
cook["rain"][cook["rain"]=="trace"] <- 0.0001
cook$rain <- as.numeric(cook$rain)
cook <- subset(cook, Year=="2021")
cook$Date<- factor(cook$Date, levels = c("4/15/21","6/4/21","6/11/21","6/18/21","6/25/21","7/2/21","7/9/21","7/16/21","7/23/21","7/30/21","8/6/21","8/13/21","8/20/21","8/27/21","9/3/21","9/10/21","9/17/21","9/24/21"))
cook$rain %>% class
cook["rain"][cook["rain"]=="trace"] <- 0.001
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

summary(lm(Pupae.prev ~ pH.avg, data = cook.day))
summary(lm(Pupae.prev ~ Temp.C.avg, data = cook.day))
summary(lm(Pupae.prev ~ Cond.avg, data = cook.day))
summary(lm(Pupae.prev ~ DO.avg, data = cook.day))
summary(lm(Pupae.prev ~ Sal.avg, data = cook.day))

summary(lm(Pupae.avg ~ pH.avg, data = cook.day))
summary(lm(Pupae.avg ~ Temp.C.avg, data = cook.day))
summary(lm(Pupae.avg ~ Cond.avg, data = cook.day))
summary(lm(Pupae.avg ~ DO.avg, data = cook.day))
summary(lm(Pupae.avg ~ Sal.avg, data = cook.day))


## Table 1C Linear mixed-effects models with unaggregated data (observations only where pupae were present)

cook <- read.table("WQ_Dips_2020_2021.csv",sep=",",header=TRUE)
cook <- subset(cook, Year=="2021")
cook <- subset(cook, Date > "2021-04-16")
cook.pupaepresent <- subset(cook, Pupae > "0")
cook.pupaepresent <- subset(cook.pupaepresent, !is.na(pH)==TRUE)
cook.pupaepresent$Datefactor <- cook.pupaepresent$Date %>% as.factor
cook.pupaepresent = cook.pupaepresent %>%
     arrange(UCB, Datefactor)                     
model = lme(sqrt(Pupae) ~ pH,
              random = ~1|UCB, 
              correlation = corAR1(), 
              data = cook.pupaepresent)
summary(model)
model = lme(sqrt(Pupae) ~ Temp.C,
              random = ~1|UCB, 
              correlation = corAR1(), 
              data = cook.pupaepresent)
summary(model)
model = lme(sqrt(Pupae) ~ Cond,
              random = ~1|UCB, 
              correlation = corAR1(), 
              data = cook.pupaepresent)
summary(model)
model = lme(sqrt(Pupae) ~ DO,
              random = ~1|UCB, 
              correlation = corAR1(), 
              data = cook.pupaepresent)
summary(model)
model = lme(sqrt(Pupae) ~ Sal,
              random = ~1|UCB, 
              correlation = corAR1(), 
              data = cook.pupaepresent)
summary(model)
