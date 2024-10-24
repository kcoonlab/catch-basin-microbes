set.seed(123)

## Table 1. Effects of water quality on mosquito productivity

library(dplyr)
library(nlme)

## Table 1A. Linear regression analyses with data aggregated by basin

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data["Rainfall"][WQ.data["Rainfall"]=="trace"] <- 0.001
WQ.data$Rainfall <- as.numeric(WQ.data$Rainfall)
WQ.data.by.basins <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
            Pupae.abund.sd = sd(Pupae.abund, na.rm=TRUE),
            Pupae.abund.min = min(Pupae.abund, na.rm=TRUE),
            Pupae.abund.max = max(Pupae.abund, na.rm=TRUE),
            Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
            Fail.prop = sum(Methoprene.fail==1, na.rm=TRUE)/sum(Methoprene.fail>=0,na.rm=TRUE),            
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

summary(lm(Pupae.prev ~ pH.avg, data = WQ.data.by.basins))
summary(lm(Pupae.prev ~ Temp.C.avg, data = WQ.data.by.basins))
summary(lm(Pupae.prev ~ Cond.avg, data = WQ.data.by.basins))
summary(lm(Pupae.prev ~ DO.avg, data = WQ.data.by.basins))
summary(lm(Pupae.prev ~ Sal.avg, data = WQ.data.by.basins))

summary(lm(Pupae.abund.avg ~ pH.avg, data = WQ.data.by.basins))
summary(lm(Pupae.abund.avg ~ Temp.C.avg, data = WQ.data.by.basins))
summary(lm(Pupae.abund.avg ~ Cond.avg, data = WQ.data.by.basins))
summary(lm(Pupae.abund.avg ~ DO.avg, data = WQ.data.by.basins))
summary(lm(Pupae.abund.avg ~ Sal.avg, data = WQ.data.by.basins))

## Table 1B. Linear regression analyses with data aggregated by sampling date

WQ.data.by.date <- WQ.data %>%
  group_by(Sampling.date) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
            Pupae.abund.sd = sd(Pupae.abund, na.rm=TRUE),
            Pupae.abund.min = min(Pupae.abund, na.rm=TRUE),
            Pupae.abund.max = max(Pupae.abund, na.rm=TRUE),
            Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
            Fail.prop = sum(Methoprene.fail==1, na.rm=TRUE)/sum(Methoprene.fail>=0,na.rm=TRUE),            
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
            Cond.max = max(Cond, na.rm=TRUE),
            Rainfall.avg = mean(Rainfall, na.rm=TRUE)
            )
WQ.data.by.date[sapply(WQ.data.by.date, is.nan)] <- NA
WQ.data.by.date[sapply(WQ.data.by.date, is.infinite)] <- NA

summary(lm(Pupae.prev ~ pH.avg, data = WQ.data.by.date))
summary(lm(Pupae.prev ~ Temp.C.avg, data = WQ.data.by.date))
summary(lm(Pupae.prev ~ Cond.avg, data = WQ.data.by.date))
summary(lm(Pupae.prev ~ DO.avg, data = WQ.data.by.date))
summary(lm(Pupae.prev ~ Sal.avg, data = WQ.data.by.date))
summary(lm(Pupae.prev ~ Rainfall.avg, data = WQ.data.by.date))

summary(lm(Pupae.abund.avg ~ pH.avg, data = WQ.data.by.date))
summary(lm(Pupae.abund.avg ~ Temp.C.avg, data = WQ.data.by.date))
summary(lm(Pupae.abund.avg ~ Cond.avg, data = WQ.data.by.date))
summary(lm(Pupae.abund.avg ~ DO.avg, data = WQ.data.by.date))
summary(lm(Pupae.abund.avg ~ Sal.avg, data = WQ.data.by.date))
summary(lm(Pupae.abund.avg ~ Rainfall.avg, data = WQ.data.by.date))

## Table 1C. Linear mixed-effects models with unaggregated data (observations only where pupae were present)

WQ.data.pupaepresent <- subset(WQ.data, Pupae.abund > "0")
WQ.data.pupaepresent$Datefactor <- WQ.data.pupaepresent$Sampling.date %>% as.factor

WQ.data.pupaepresent.pH <- subset(WQ.data.pupaepresent, !is.na(pH)==TRUE)
WQ.data.pupaepresent.pH = WQ.data.pupaepresent.pH %>%
     arrange(Basin.id, Datefactor)                     
model = lme(sqrt(Pupae.abund) ~ pH,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent.pH)
summary(model)

WQ.data.pupaepresent.Temp.C <- subset(WQ.data.pupaepresent, !is.na(Temp.C)==TRUE)
WQ.data.pupaepresent.Temp.C = WQ.data.pupaepresent.Temp.C %>%
     arrange(Basin.id, Datefactor)      
model = lme(sqrt(Pupae.abund) ~ Temp.C,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent.Temp.C)
summary(model)

WQ.data.pupaepresent.Cond <- subset(WQ.data.pupaepresent, !is.na(Cond)==TRUE)
WQ.data.pupaepresent.Cond = WQ.data.pupaepresent.Cond %>%
     arrange(Basin.id, Datefactor) 
model = lme(sqrt(Pupae.abund) ~ Cond,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent.Cond)
summary(model)

WQ.data.pupaepresent.DO <- subset(WQ.data.pupaepresent, !is.na(DO)==TRUE)
WQ.data.pupaepresent.DO = WQ.data.pupaepresent.DO %>%
     arrange(Basin.id, Datefactor)
model = lme(sqrt(Pupae.abund) ~ DO,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent.DO)
summary(model)

WQ.data.pupaepresent.Sal <- subset(WQ.data.pupaepresent, !is.na(Sal)==TRUE)
WQ.data.pupaepresent.Sal = WQ.data.pupaepresent.Sal %>%
     arrange(Basin.id, Datefactor)
model = lme(sqrt(Pupae.abund) ~ Sal,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent.Sal)
summary(model)
