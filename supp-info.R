set.seed(123)

## Supplementary information for: Microbiota composition associates with mosquito productivity outcomes in belowground larval habitats

library(ggplot2)
library(dplyr)
library(tree)
library(nlme)
library(phyloseq)
library(ape)
library(vegan)
library(philr)

## Fig S1. Pupal ocurrence and abundance over time

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data$Sampling.date <- as.Date(WQ.data$Sampling.date, "%m/%d/%y")
WQ.data <- subset(WQ.data, Sampling.date > "2021-04-16")

WQ.data.by.date <- WQ.data %>%
  group_by(Sampling.date) %>%
  dplyr::summarize(Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE))

ggplot(data=WQ.data, aes(x=Sampling.date)) +
  geom_boxplot(aes(y=sqrt(Pupae.abund), group=Sampling.date), alpha=0.3) +
  geom_line(data=WQ.data.by.date, aes(y=Pupae.prev*10), color="red") +
  scale_y_continuous(name="Pupae abundance per basin (sqrt)", sec.axis=sec_axis(~.*10, name="Percent basins with pupae present"))+
  theme(axis.title = element_text(size=11), axis.title.y.right = element_text(colour = "black"), 
        axis.text.y.right = element_text(color="red"), axis.line.y.right = element_line(color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Fig S2. Mosquito productivity by fixed basin characteristics

## Fig S2A. Pupal occurence by basin type

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data$Sampling.date <- as.Date(WQ.data$Sampling.date, "%m/%d/%y")

WQ.data.by.basin <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
            Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE)
            )
WQ.data.by.basin$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
WQ.data.by.basin$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Stratford","Stratford","Gibbons","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")

combsep_pupae <- ggplot(WQ.data.by.basin, aes(x=Basin.type, y=Pupae.prev, fill=Basin.type)) +
  geom_boxplot() + 
  xlab(NULL) +
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.position="blank") +
  ylab("Pupae frequency") +
  scale_fill_manual(values=c("combined"="grey","separate"="white")) +
  scale_x_discrete(breaks=c("combined","separate"), labels=c("Combined", "Separated"))
combsep_pupae

anova.basintype <- aov(Pupae.prev ~ Basin.type, data=WQ.data.by.basin)
summary(anova.basintype)

## Fig S2B. Pupal occurence by flow group

WQ.data.by.basin$Basin.flowgroup <- factor(WQ.data.by.basin$Basin.flowgroup, levels=c("Donald-Banta", "Gibbons", "MayfairCarlyle", "MinerEvanstonRammer","Stratford"))

flowgroup_pupae <- ggplot(WQ.data.by.basin, aes(x=Basin.flowgroup, y=Pupae.prev, fill=Basin.flowgroup)) +
  geom_boxplot() +
  scale_fill_manual(values=c("Donald-Banta"="black", "Gibbons"="dark grey","MayfairCarlyle"="grey", "MinerEvanstonRammer"="light grey", "Stratford"="white")) +
  scale_x_discrete(breaks=c("Gibbons","MayfairCarlyle","Stratford","Donald-Banta","MinerEvanstonRammer"), 
                   labels=c("Gibbons","Mayfair-Carlyle","Stratford","Donald-nBanta","Miner-Evanston-Rammer")) +
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.key=element_rect(fill="white"),legend.position="blank") +
  ylab("Pupae frequency")+ xlab(NULL)+ labs(fill="Flow group")
flowgroup_pupae

anova.flowgroup <- aov(Pupae.prev ~ Basin.flowgroup, data=WQ.data.by.basin)
summary(anova.flowgroup)

## Fig S2C. Pupal abundance in "Combined" versus "Separated" basins over time

WQ.data$Datefactor <- WQ.data$Sampling.date %>% as.factor

pupae_flowgroup_time <- ggplot(WQ.data, aes(y=sqrt(Pupae.abund), x=Datefactor, fill=Basin.type)) +
  geom_boxplot() +
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), 
        legend.position="blank") +
  scale_x_discrete(breaks=c("2021-04-15","2021-06-04","2021-06-11","2021-06-18","2021-06-25","2021-07-02","2021-07-09","2021-07-16","2021-07-23","2021-07-30","2021-08-06","2021-08-13","2021-08-20","2021-08-27","2021-09-03","2021-09-10","2021-09-17","2021-09-24"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24")) +
  xlab("Date") + ylab("Pupae abundance (sqrt)")+ labs(fill="Basin type") +
  scale_fill_manual(values=c("combined"="grey","separate"="white"), 
                    breaks=c("combined","separate"), labels=c("Combined", "Separated"))
pupae_flowgroup_time

## Fig S3. Classification tree for biotype cluster assignment of samples

metadata <- read.table("catch-basin-microbes/input-files/metadata.txt", sep="\t", header=TRUE)
cols.convert.factor <-c("Sampling.date","Biotype","Basin.type")
metadata[,cols.convert.factor] <- lapply(metadata[,cols.convert.factor], factor)
metadata.samples <- subset(metadata, sample_control=="sample")

testtree <- tree(data=metadata.samples, Biotype ~ Sampling.date + Basin.type + pH + Temp.C + Cond + DO + Sal)
plot(testtree)
text(testtree)

## Fig S4. Alpha diversity of basin-associated microbiota by fixed basin characteristics

## Fig S4A. Alpha diversity in "Combined" versus "Separated" basins over time

metadata.samples$Sampling.date <- factor(metadata.samples$Sampling.date, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))

shannon.date <- ggplot(data=metadata.samples, 
                        aes(x=Sampling.date, y=Shannon, fill=Basin.type)) +
  geom_boxplot() +
  scale_fill_manual(values=c("combined"="grey","separated"="white"),
                    breaks=c("combined","separated"), labels=c("Combined", "Separated")) +
  scale_x_discrete(breaks=c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"),
                   labels=c("4/15","6/11","6/25","7/9","7/23","8/6","8/27","9/17")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL)+ ylab("Shannon index")+labs(fill="Basin type")
shannon.date

## Fig S4B. Alpha diversity in basins assigned to different flow groups over time

metadata.samples$Basin.flowgroup <- factor(metadata.samples$Basin.flowgroup, levels=c("DonaldBanta", "Gibbons", "MayfairCarlyle", "MinerEvanstonRammer""Stratford"))
shannon.date <- ggplot(data=metadata.samples, 
                        aes(x=Sampling.date, y=Shannon, fill=Basin.flowgroup)) +
  geom_boxplot() +
  scale_fill_manual(values=c("DonaldBanta"="black", "Gibbons"="dark grey","MayfairCarlyle"="grey", "MinerEvanstonRammer"="light grey", "Stratford"="white")) +
  scale_x_discrete(breaks=c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"),
                   labels=c("4/15","6/11","6/25","7/9","7/23","8/6","8/27","9/17")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL) + ylab("Shannon index") + labs(fill="Flow group")
shannon.date

## Fig S5. Principal coordinates analysis of PhILR Euclidean distances

ps.final <- readRDS("catch-basin-microbes/input-files/ps.final.rds")
sample_data(ps.final)$sample.id <- sample_names(ps.final)

# Running with just individual NA samples cut (not final instances of a UCB with an NA anywhere)
NA %in% sample_data(ps.final)$pH #row 39
NA %in% sample_data(ps.final)$Temp # row 39,129
NA %in% sample_data(ps.final)$DO # row 39
NA %in% sample_data(ps.final)$Sal # row 39, 71
NA %in% sample_data(ps.final)$Cond # row 39
sample_data(ps.final)[39,] #sample.id = B_73
sample_data(ps.final)[71,] #sample.id = C_61
sample_data(ps.final)[129,] #sample.id = F_56

# Take out the 3 individual samples that have NAs, and use that remaining set for the permanovas
ps.most <- subset_samples(ps.final, !sample.id %in% c("B_73","C_61","F_56"))
samples.most <- as.data.frame(as.matrix(sample_data(ps.most)))

ps.most.pseudocount <- ps.most
otu_table(ps.most.pseudocount) <- otu_table(ps.most) + 1
ps.most.philr <- philr(ps.most.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
dist.most.philr <- dist(ps.most.philr, method="euclidean")

ord.most.philr <- ordinate(ps.most, method="PCoA", distance = dist.most.philr)

plot.ord.philr <- plot_ordination(ps.most, ord.most.philr, type = "samples", color = "Sampling.date", shape = "Basin.type") + 
  geom_point(size=3) +
  scale_color_discrete(breaks=c(), labels=c()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill="white"))
plot.ord.philr

plot.ord.philr <- plot_ordination(ps.most, ord.most.philr, type = "samples", color = "Biotype") + 
  geom_point(size=3) +
  scale_x_discrete(breaks=c("1","2"), labels=c("A", "B")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill="white"))
plot.ord.philr

## Fig S6. Water quality measurements over time

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data$Date <- as.Date(WQ.data$Date, "%m/%d/%y")
WQ.data$Datefactor <- WQ.data$Date %>% as.factor

pH.date <- ggplot(data=WQ.data, 
                        aes(x=Datefactor, y=pH)) +
  geom_boxplot() +
  scale_x_discrete(breaks=c("2021-04-15","2021-06-04","2021-06-11","2021-06-18","2021-06-25","2021-07-02","2021-07-09","2021-07-16","2021-07-23","2021-07-30","2021-08-06","2021-08-13","2021-08-20","2021-08-27","2021-09-03","2021-09-10","2021-09-17","2021-09-24"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL)+ ylab("pH")+labs(fill="Date")
pH.date

Temp.date <- ggplot(data=WQ.data, 
                        aes(x=Datefactor, y=Temp.C)) +
  geom_boxplot() +
  scale_x_discrete(breaks=c("2021-04-15","2021-06-04","2021-06-11","2021-06-18","2021-06-25","2021-07-02","2021-07-09","2021-07-16","2021-07-23","2021-07-30","2021-08-06","2021-08-13","2021-08-20","2021-08-27","2021-09-03","2021-09-10","2021-09-17","2021-09-24"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL)+ ylab("Temp")+labs(fill="Date")
Temp.date

Cond.date <- ggplot(data=WQ.data, 
                        aes(x=Datefactor, y=Cond)) +
  geom_boxplot() +
  scale_x_discrete(breaks=c("2021-04-15","2021-06-04","2021-06-11","2021-06-18","2021-06-25","2021-07-02","2021-07-09","2021-07-16","2021-07-23","2021-07-30","2021-08-06","2021-08-13","2021-08-20","2021-08-27","2021-09-03","2021-09-10","2021-09-17","2021-09-24"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL)+ ylab("Cond")+labs(fill="Date")
Cond.date

DO.date <- ggplot(data=WQ.data, 
                        aes(x=Datefactor, y=DO)) +
  geom_boxplot() +
  scale_x_discrete(breaks=c("2021-04-15","2021-06-04","2021-06-11","2021-06-18","2021-06-25","2021-07-02","2021-07-09","2021-07-16","2021-07-23","2021-07-30","2021-08-06","2021-08-13","2021-08-20","2021-08-27","2021-09-03","2021-09-10","2021-09-17","2021-09-24"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL)+ ylab("DO")+labs(fill="Date")
DO.date

Sal.date <- ggplot(data=WQ.data, 
                        aes(x=Datefactor, y=Sal)) +
  geom_boxplot() +
  scale_x_discrete(breaks=c("2021-04-15","2021-06-04","2021-06-11","2021-06-18","2021-06-25","2021-07-02","2021-07-09","2021-07-16","2021-07-23","2021-07-30","2021-08-06","2021-08-13","2021-08-20","2021-08-27","2021-09-03","2021-09-10","2021-09-17","2021-09-24"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL)+ ylab("Sal")+labs(fill="Date")
Sal.date

WQ.data["rain"][WQ.data["rain"]=="trace"] <- 0.001
WQ.data$rain <- as.numeric(WQ.data$rain)

rain.date <- ggplot(data=WQ.data, 
                        aes(x=Datefactor, y=rain)) +
  geom_bar(stat="identity") +
  scale_x_discrete(breaks=c("2021-04-15","2021-06-04","2021-06-11","2021-06-18","2021-06-25","2021-07-02","2021-07-09","2021-07-16","2021-07-23","2021-07-30","2021-08-06","2021-08-13","2021-08-20","2021-08-27","2021-09-03","2021-09-10","2021-09-17","2021-09-24"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  xlab(NULL)+ ylab("Rainffinal")+labs(fill="Date")
rain.date

## Table S1 (See catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv for unaggregated data)

## Table S1A. Fixed basin characteristics along with measures of water quality and mosquito productivity aggregated by basin across the entire sampling period

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data["Rainffinal"][WQ.data["Rainffinal"]=="trace"] <- 0.001
WQ.data$Rainffinal <- as.numeric(WQ.data$Rainffinal)
WQ.data.by.basin <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
            Pupae.abund.sd = sd(Pupae.abund, na.rm=TRUE),
            Pupae.abund.min = min(Pupae.abund, na.rm=TRUE),
            Pupae.abund.max = max(Pupae.abund, na.rm=TRUE),
            Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
            Methoprene.fail.rate = sum(Methoprene.success==1, na.rm=TRUE)/sum(Methoprene.success>=0,na.rm=TRUE), 
            Rainffinal.avg = mean(Rainffinal, na.rm=TRUE),
            Rainffinal.sd = sd(Rainffinal, na.rm=TRUE),
            Rainffinal.min = min(Rainffinal, na.rm=TRUE),
            Rainffinal.max = max(Rainffinal, na.rm=TRUE),
            pH.avg = mean(pH, na.rm=TRUE),
            pH.sd = sd(pH, na.rm=TRUE),
            pH.min = min(pH, na.rm=TRUE),
            pH.max = max(pH, na.rm=TRUE),
            Temp.C.avg = mean(Temp.C, na.rm=TRUE),
            Temp.C.sd = sd(Temp.C, na.rm=TRUE),
            Temp.C.min = min(Temp.C, na.rm=TRUE),
            Temp.C.max = max(Temp.C, na.rm=TRUE),
            Cond.avg = mean(Cond, na.rm=TRUE),
            Cond.sd = sd(Cond, na.rm=TRUE),
            Cond.min = min(Cond, na.rm=TRUE),
            Cond.max = max(Cond, na.rm=TRUE),
            DO.avg = mean(DO, na.rm=TRUE),
            DO.sd = sd(DO, na.rm=TRUE),
            DO.min = min(DO, na.rm=TRUE),
            DO.max = max(DO, na.rm=TRUE),
            Sal.avg = mean(Sal, na.rm=TRUE),
            Sal.sd = sd(Sal, na.rm=TRUE),
            Sal.min = min(Sal, na.rm=TRUE),
            Sal.max = max(Sal, na.rm=TRUE)
            )

WQ.data.by.basin$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
WQ.data.by.basin$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","DonaldBanta","DonaldBanta","DonaldBanta","DonaldBanta","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Stratford","Stratford","Gibbons","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")
WQ.data.by.basin$Basin.lon <- c(-87.957901,-87.957907,-87.957481,-87.957248,-87.957257,-87.956212,-87.955063,-87.954987,-87.955064,-87.954972,-87.955073,-87.954946,-87.95384,-87.953765,-87.953832,-87.953746,-87.954378,-87.954381,-87.953833,-87.953712,-87.953743,-87.952483,-87.952659,-87.95254,-87.952665,-87.952544,-87.952664,-87.957397,-87.956346,-87.954944,-87.9622283,-87.96229318,-87.96231902,-87.96222278,-87.9623084,-87.96220848,-87.96220773,-87.96231615,-87.96118708,-87.96128735,-87.96217608,-87.96222302,-87.96231955,-87.96237833,-87.96229127,-87.96109237,-87.96120367,-87.96107445,-87.96582388,-87.96585632,-87.96571188,-87.96835267,-87.96847103,-87.96832385,-87.96844195,-87.96723588,-87.96715035,-87.96706973,-87.96121885,-87.96116213)
WQ.data.by.basin$Basin.lat <- c(42.085423,42.085404,42.084462,42.084567,42.084494,42.084527,42.083446,42.083488,42.08274,42.08276,42.081956,42.081948,42.082057,42.082044,42.083042,42.083024,42.084479,42.084463,42.084031,42.084027,42.084108,42.08413,42.083495,42.083501,42.081681,42.081668,42.083969,42.082795,42.084539,42.084451,42.08059342,42.08061185,42.07995083,42.07996232,42.0791431,42.07914033,42.07905157,42.0790679,42.07722598,42.07727365,42.07733793,42.07735468,42.07737213,42.07733773,42.07726987,42.07638237,42.07570793,42.07569713,42.0775284,42.07762897,42.07763798,42.07911363,42.07911247,42.08001807,42.08003452,42.07917337,42.07909997,42.079124,42.07906305,42.0793045)
WQ.data.by.basin$Treatment.status <- c("treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","not_treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated","treated")
WQ.data.by.basin <- WQ.data.by.basin[, c(1,32:35,2:6,36,7:31)]

## Table S1B. Measures of water quality, mosquito productivity, and rainffinal across final sampled basins on a given sampling date

WQ.data.by.date <- WQ.data %>%
  group_by(Sampling.date) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
            Pupae.abund.sd = sd(Pupae.abund, na.rm=TRUE),
            Pupae.abund.min = min(Pupae.abund, na.rm=TRUE),
            Pupae.abund.max = max(Pupae.abund, na.rm=TRUE),
            Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
            Methoprene.fail.rate = sum(Methoprene.success==1, na.rm=TRUE)/sum(Methoprene.success>=0,na.rm=TRUE), 
            Rainffinal.avg = mean(Rainffinal, na.rm=TRUE),
            Rainffinal.sd = sd(Rainffinal, na.rm=TRUE),
            Rainffinal.min = min(Rainffinal, na.rm=TRUE),
            Rainffinal.max = max(Rainffinal, na.rm=TRUE),
            pH.avg = mean(pH, na.rm=TRUE),
            pH.sd = sd(pH, na.rm=TRUE),
            pH.min = min(pH, na.rm=TRUE),
            pH.max = max(pH, na.rm=TRUE),
            Temp.C.avg = mean(Temp.C, na.rm=TRUE),
            Temp.C.sd = sd(Temp.C, na.rm=TRUE),
            Temp.C.min = min(Temp.C, na.rm=TRUE),
            Temp.C.max = max(Temp.C, na.rm=TRUE),
            Cond.avg = mean(Cond, na.rm=TRUE),
            Cond.sd = sd(Cond, na.rm=TRUE),
            Cond.min = min(Cond, na.rm=TRUE),
            Cond.max = max(Cond, na.rm=TRUE),
            DO.avg = mean(DO, na.rm=TRUE),
            DO.sd = sd(DO, na.rm=TRUE),
            DO.min = min(DO, na.rm=TRUE),
            DO.max = max(DO, na.rm=TRUE),
            Sal.avg = mean(Sal, na.rm=TRUE),
            Sal.sd = sd(Sal, na.rm=TRUE),
            Sal.min = min(Sal, na.rm=TRUE),
            Sal.max = max(Sal, na.rm=TRUE)
            )
WQ.data.by.date[sapply(WQ.data.by.date, is.nan)] <- NA
WQ.data.by.date[sapply(WQ.data.by.date, is.infinite)] <- NA
WQ.data.by.date <- WQ.data.by.date %>% arrange(factor(Sampling.date, levels = c("4/15/21", "6/4/21", "6/11/21", "6/18/21", "6/25/21", "7/2/21", "7/9/21", "7/16/21", "7/23/21", "7/30/21", "8/6/21", "8/13/21", "8/20/21", "8/27/21", "9/3/21", "9/10/21", "9/17/21", "9/24/21")))

## Table S2 **See catch-basin-microbes/input-files/metadata.txt**

## Table S3. Results of spatial autocorrelation analyses for testing dependence among basins

## Table S3A. Moran's I autocorrelation index with data aggregated by basin

set.seed(123)

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data["Rainffinal"][WQ.data["Rainffinal"]=="trace"] <- 0.001
WQ.data$Rainffinal <- as.numeric(WQ.data$Rainffinal)
WQ.data.by.basin <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
            Pupae.abund.sd = sd(Pupae.abund, na.rm=TRUE),
            Pupae.abund.min = min(Pupae.abund, na.rm=TRUE),
            Pupae.abund.max = max(Pupae.abund, na.rm=TRUE),
            Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
            Methoprene.fail.rate = sum(Methoprene.success==1, na.rm=TRUE)/sum(Methoprene.success>=0,na.rm=TRUE), 
            Rainffinal.avg = mean(Rainffinal, na.rm=TRUE),
            Rainffinal.sd = sd(Rainffinal, na.rm=TRUE),
            Rainffinal.min = min(Rainffinal, na.rm=TRUE),
            Rainffinal.max = max(Rainffinal, na.rm=TRUE),
            pH.avg = mean(pH, na.rm=TRUE),
            pH.sd = sd(pH, na.rm=TRUE),
            pH.min = min(pH, na.rm=TRUE),
            pH.max = max(pH, na.rm=TRUE),
            Temp.C.avg = mean(Temp.C, na.rm=TRUE),
            Temp.C.sd = sd(Temp.C, na.rm=TRUE),
            Temp.C.min = min(Temp.C, na.rm=TRUE),
            Temp.C.max = max(Temp.C, na.rm=TRUE),
            Cond.avg = mean(Cond, na.rm=TRUE),
            Cond.sd = sd(Cond, na.rm=TRUE),
            Cond.min = min(Cond, na.rm=TRUE),
            Cond.max = max(Cond, na.rm=TRUE),
            DO.avg = mean(DO, na.rm=TRUE),
            DO.sd = sd(DO, na.rm=TRUE),
            DO.min = min(DO, na.rm=TRUE),
            DO.max = max(DO, na.rm=TRUE),
            Sal.avg = mean(Sal, na.rm=TRUE),
            Sal.sd = sd(Sal, na.rm=TRUE),
            Sal.min = min(Sal, na.rm=TRUE),
            Sal.max = max(Sal, na.rm=TRUE)
            )

WQ.data.by.basin$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
WQ.data.by.basin$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","DonaldBanta","DonaldBanta","DonaldBanta","DonaldBanta","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Stratford","Stratford","Gibbons","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")
WQ.data.by.basin$Basin.lon <- c(-87.957901,-87.957907,-87.957481,-87.957248,-87.957257,-87.956212,-87.955063,-87.954987,-87.955064,-87.954972,-87.955073,-87.954946,-87.95384,-87.953765,-87.953832,-87.953746,-87.954378,-87.954381,-87.953833,-87.953712,-87.953743,-87.952483,-87.952659,-87.95254,-87.952665,-87.952544,-87.952664,-87.957397,-87.956346,-87.954944,-87.9622283,-87.96229318,-87.96231902,-87.96222278,-87.9623084,-87.96220848,-87.96220773,-87.96231615,-87.96118708,-87.96128735,-87.96217608,-87.96222302,-87.96231955,-87.96237833,-87.96229127,-87.96109237,-87.96120367,-87.96107445,-87.96582388,-87.96585632,-87.96571188,-87.96835267,-87.96847103,-87.96832385,-87.96844195,-87.96723588,-87.96715035,-87.96706973,-87.96121885,-87.96116213)
WQ.data.by.basin$Basin.lat <- c(42.085423,42.085404,42.084462,42.084567,42.084494,42.084527,42.083446,42.083488,42.08274,42.08276,42.081956,42.081948,42.082057,42.082044,42.083042,42.083024,42.084479,42.084463,42.084031,42.084027,42.084108,42.08413,42.083495,42.083501,42.081681,42.081668,42.083969,42.082795,42.084539,42.084451,42.08059342,42.08061185,42.07995083,42.07996232,42.0791431,42.07914033,42.07905157,42.0790679,42.07722598,42.07727365,42.07733793,42.07735468,42.07737213,42.07733773,42.07726987,42.07638237,42.07570793,42.07569713,42.0775284,42.07762897,42.07763798,42.07911363,42.07911247,42.08001807,42.08003452,42.07917337,42.07909997,42.079124,42.07906305,42.0793045)
WQ.data.by.basin <- WQ.data.by.basin[, c(1,32:35,2:31)]

data.dists <- as.matrix(dist(cbind(WQ.data.by.basin$Basin.lon, WQ.data.by.basin$Basin.lat)))
data.dists.inv <- 1/data.dists
diag(data.dists.inv) <- 0
Moran.I(WQ.data.by.basin$Pupae.prev, data.dists.inv, na.rm = TRUE)
Moran.I(WQ.data.by.basin$Pupae.abund.avg, data.dists.inv, na.rm = TRUE)
Moran.I(WQ.data.by.basin$Methoprene.success, data.dists.inv, na.rm = TRUE) 

metadata <- read.table("catch-basin-microbes/input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
metadata.byUCB <- metadata %>% 
  group_by(Basin.id) %>%
  dplyr::summarize(Shannon.avg = mean(Shannon, na.rm=TRUE),
  				   ASV.richness.avg = mean(ASV.richness, na.rm=TRUE),
  				   C39.relabund.avg = mean(C39.relabund, na.rm=TRUE),
  				   Proteobacteria.relabund.avg = mean(Proteobacteria.relabund, na.rm=TRUE),
  				   Firmicutes.relabund.avg = mean(Firmicutes.relabund, na.rm=TRUE))
  				   	   
metadata <- metadata[!duplicated(metadata$Basin.id),]
data <- merge(metadata, metadata.byUCB, by = "Basin.id")
data.dists <- as.matrix(dist(cbind(data$Basin.lon, data$Basin.lat)))
data.dists.inv <- 1/data.dists
data.dists.inv[sapply(data.dists.inv, is.infinite)] <- 0
diag(data.dists.inv) <- 0
Moran.I(data$Shannon.avg, data.dists.inv, na.rm = TRUE)
Moran.I(data$ASV.richness.avg, data.dists.inv, na.rm = TRUE)

## Table S3B. Mantel test using PhILR distances between samples aggregated by basin

ps.final <- readRDS("catch-basin-microbes/input-files/ps.final.rds")
ps.final.byUCB <- merge_samples(ps.final,"Basin.id")
sample_data(ps.final.byUCB)$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
sample_data(ps.final.byUCB)$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","DonaldBanta","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")
sample_data(ps.final.byUCB)$Basin.lon <- c(-87.957901,-87.957907,-87.957481,-87.957248,-87.954972,-87.955073,-87.954946,-87.95384,-87.953765,-87.954378,-87.953833,-87.953712,-87.953743,-87.952483,-87.952665,-87.952544,-87.952664,-87.957397,-87.956346,-87.954944,-87.96229318,-87.96222278,-87.9623084,-87.96220848,-87.96231615,-87.96128735,-87.96118708,-87.96217608,-87.96222302,-87.96229127,-87.96231955,-87.96109237,-87.96120367,-87.96571188,-87.96585632,-87.96835267,-87.96847103,-87.96844195,-87.96832385,-87.96706973,-87.96723588,-87.96715035)
sample_data(ps.final.byUCB)$Basin.lat <- c(42.085423,42.085404,42.084462,42.084567,42.08276,42.081956,42.081948,42.082057,42.082044,42.084479,42.084031,42.084027,42.084108,42.08413,42.081681,42.081668,42.083969,42.082795,42.084539,42.084451,42.08061185,42.07996232,42.0791431,42.07914033,42.0790679,42.07727365,42.07722598,42.07733793,42.07735468,42.07726987,42.07737213,42.07638237,42.07570793,42.07763798,42.07762897,42.07911363,42.07911247,42.08003452,42.08001807,42.079124,42.07917337,42.07909997)
samples.final.byUCB <- as.data.frame(as.matrix(sample_data(ps.final.byUCB)))

ps.final.byUCB.pseudocount <- ps.final.byUCB
otu_table(ps.final.byUCB.pseudocount) <- t(otu_table(ps.final.byUCB) + 1)
ps.final.byUCB.philr <- philr(ps.final.byUCB.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
dist.final.byUCB.philr <- dist(ps.final.byUCB.philr, method="euclidean")

data.dists <- as.matrix(dist(cbind(samples.final.byUCB$Basin.lon, samples.final.byUCB$Basin.lat)))
data.dists.inv <- 1/data.dists
data.dists.inv[sapply(data.dists.inv, is.infinite)] <- 0
diag(data.dists.inv) <- 0
mantel(dist.final.byUCB.philr, data.dists.inv)

## Table S4. Effects of methoprene treatment on mosquito productivity

## Table S4A. Linear regression analyses with data aggregated by basin

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data.by.basins <- WQ.data %>%
  group_by(UCB) %>%
  dplyr::summarize(Pupae.avg = mean(Pupae, na.rm=TRUE),
            Pupae.prev = sum(Pupae>=1, na.rm=TRUE)/sum(Pupae>=0, na.rm=TRUE),
            Fail.prop = sum(Fail==1, na.rm=TRUE)/sum(Fail>=0,na.rm=TRUE)            
            )

summary(lm(Pupae.prev ~ Fail.prop, data = WQ.data.by.basins))
summary(lm(Pupae.avg ~ Fail.prop, data = WQ.data.by.basins))

## Table S4B. Linear regression analyses with data aggregated by date

WQ.data.by.date <- WQ.data %>%
  group_by(Date) %>%
  dplyr::summarize(Pupae.avg = mean(Pupae, na.rm=TRUE),
            Pupae.prev = sum(Pupae>=1, na.rm=TRUE)/sum(Pupae>=0, na.rm=TRUE),
            Fail.prop = sum(Fail==1, na.rm=TRUE)/sum(Fail>=0,na.rm=TRUE)            
            )
WQ.data.by.date[sapply(WQ.data.by.date, is.nan)] <- NA

summary(lm(Pupae.prev ~ Fail.prop, data = WQ.data.by.date))
summary(lm(Pupae.avg ~ Fail.prop, data = WQ.data.by.date))

## Table S4C. Linear mixed-effects models with unaggregated data (observations only where pupae were present)

set.seed(123)

WQ.data.pupaepresent <- subset(WQ.data, Pupae > "0")
WQ.data.pupaepresent$Datefactor <- WQ.data.pupaepresent$Date %>% as.factor

WQ.data.pupaepresent.Fail <- subset(WQ.data.pupaepresent, !is.na(Fail)==TRUE)
WQ.data.pupaepresent.Fail = WQ.data.pupaepresent.Fail %>%
     arrange(UCB, Datefactor)                     
model = lme(sqrt(Pupae) ~ Fail,
              random = ~1|UCB, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent.Fail)
summary(model)

## Table S5. Effects of water quality on methoprene treatment success

## Table S5A. Linear regression analyses with data aggregated by basin

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data.by.basins <- WQ.data %>%
  group_by(UCB) %>%
  dplyr::summarize(Fail.prop = sum(Fail==1, na.rm=TRUE)/sum(Fail>=0,na.rm=TRUE),            
            pH.avg = mean(pH, na.rm=TRUE),
            Temp.C.avg = mean(Temp.C, na.rm=TRUE),
            DO.avg = mean(DO, na.rm=TRUE),
            Sal.avg = mean(Sal, na.rm=TRUE),
            Cond.avg = mean(Cond, na.rm=TRUE)
            )

summary(lm(Fail.prop ~ pH.avg, data = WQ.data.by.basins))
summary(lm(Fail.prop ~ Temp.C.avg, data = WQ.data.by.basins))
summary(lm(Fail.prop ~ Cond.avg, data = WQ.data.by.basins))
summary(lm(Fail.prop ~ DO.avg, data = WQ.data.by.basins))
summary(lm(Fail.prop ~ Sal.avg, data = WQ.data.by.basins))

## Table S5B. Linear regression analyses with data aggregated by date

WQ.data.by.date <- WQ.data %>%
  group_by(Date) %>%
  dplyr::summarize(Fail.prop = sum(Fail==1, na.rm=TRUE)/sum(Fail>=0,na.rm=TRUE),            
            pH.avg = mean(pH, na.rm=TRUE),
            Temp.C.avg = mean(Temp.C, na.rm=TRUE),
            DO.avg = mean(DO, na.rm=TRUE),
            Sal.avg = mean(Sal, na.rm=TRUE),
            Cond.avg = mean(Cond, na.rm=TRUE),
            )
WQ.data.by.date[sapply(WQ.data.by.date, is.nan)] <- NA

summary(lm(Fail.prop ~ pH.avg, data = WQ.data.by.date))
summary(lm(Fail.prop ~ Temp.C.avg, data = WQ.data.by.date))
summary(lm(Fail.prop ~ Cond.avg, data = WQ.data.by.date))
summary(lm(Fail.prop ~ DO.avg, data = WQ.data.by.date))
summary(lm(Fail.prop ~ Sal.avg, data = WQ.data.by.date))

## Table S6. Effects of water quality on microbiota diversity

metadata <- read.table("catch-basin-microbes/input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")

metadata.pH <- subset(metadata, !is.na(pH)==TRUE)
metadata <- metadata %>%
     arrange(Basin.id, Date.code)    
model = lme(Shannon ~ pH,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata.pH,
              na.action=na.exclude)
summary(model)

metadata.Temp <- subset(metadata, !is.na(Temp.C)==TRUE)
model = lme(Shannon ~ Temp.C,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata.Temp,
              na.action=na.exclude)
summary(model)

metadata.Cond <- subset(metadata, !is.na(Cond)==TRUE)
model = lme(Shannon ~ Cond,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata.Cond,
              na.action=na.exclude)
summary(model)

metadata.DO <- subset(metadata, !is.na(DO)==TRUE)
model = lme(Shannon ~ DO,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata.DO,
              na.action=na.exclude)
summary(model)

metadata.Sal <- subset(metadata, !is.na(Sal)==TRUE)
model = lme(Shannon ~ Sal,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata.Sal,
              na.action=na.exclude)
summary(model)
