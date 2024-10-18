# S1

cook <- read.table("WQ_Dips_2020_2021.csv",sep=",",header=TRUE)
cook$UCB <- as.factor(cook$UCB)
library(bdscale)
library(scales)
library(dplyr)
cook$Date <- as.Date(cook$Date)
cook <- subset(cook, Date > "2021-04-16")

ggplot(data=cook, aes(x=Date))+
  geom_boxplot(aes(y=sqrt(Pupae), group=Date), alpha=0.3)+ #geom_jitter(data=cook, aes(y=sqrt(Pupae)),alpha=0.3)+
  geom_line(data=cook.bydate, aes(y=prevalence*10), color="red") +
  scale_y_continuous(name="Pupae abundance per basin (sqrt)", sec.axis=sec_axis(~.*10, name="Percent basins with pupae present"))+
  theme(axis.title = element_text(size=11), axis.title.y.right = element_text(colour = "black"), 
        axis.text.y.right = element_text(color="red"), axis.line.y.right = element_line(color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


##########################################################################################
## S2

cook.basins.2021 <- read.csv("summary_basins_2021.csv",sep=",", header=TRUE)
head(cook.basins.2021)
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

##########################################################################################
## S2A

combsep_pupae <- ggplot(cook.basins.2021, aes(x=combined_separate, y=Pupae.prev, fill=combined_separate))+
  geom_boxplot()+ xlab(NULL)+
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.position="blank")+
  ylab("Pupae frequency")+
  scale_fill_manual(values=c("combined"="grey","separate"="white"))+
  scale_x_discrete(breaks=c("combined","separate"), labels=c("Combined", "Separated"))
combsep_pupae


##########################################################################################
## S2B
flowgroup_pupae <- ggplot(cook.basins.2021, aes(x=Flowgroup_coarse, y=Pupae.prev, fill=combined_separate))+
  geom_boxplot()+
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.key=element_rect(fill="white"),legend.position="blank")+
  # scale_fill_manual(values=c("Donald-Banta"="grey","MinerEvanstonRammer"="grey","Gibbons"="white","MayfairCarlyle"="white","Stratford"="white"))+
  scale_x_discrete(breaks=c("Gibbons","MayfairCarlyle","Stratford","Donald-Banta","MinerEvanstonRammer"), 
                   labels=c("Gibbons","Mayfair\nCarlyle","Stratford","Donald\nBanta","MinerEvanston\nRammer"))+
  ylab("Pupae frequency")+ xlab(NULL)+ labs(fill="Basin type")+
  scale_fill_manual(values=c("combined"="grey","separate"="white"), 
                    breaks=c("combined","separate"), labels=c("Combined", "Separated"))
flowgroup_pupae

##########################################################################################
## S2C
pupae_flowgroup_time <- ggplot(cook, aes(y=sqrt(Pupae), x=Datefactor,fill=combined_separate))+ #, color=Flowgroup_coarse))+ 
  geom_boxplot()+
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), 
        legend.position="blank")+
        # legend.key=element_rect(fill="white"), legend.position=c(1,1), legend.justification=c(1,1))+
  scale_x_discrete(breaks=c("4/15/21","6/4/21","6/11/21","6/18/21","6/25/21","7/2/21","7/9/21","7/16/21","7/23/21","7/30/21","8/6/21","8/13/21","8/20/21","8/27/21","9/3/21","9/10/21","9/17/21","9/24/21"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24"))+
  xlab("Date") + ylab("Pupae abundance (sqrt)")+ labs(fill="Basin type")+
  scale_fill_manual(values=c("combined"="grey","separate"="white"), 
                    breaks=c("combined","separate"), labels=c("Combined", "Separated"))
pupae_flowgroup_time <- ggplot(cook, aes(y=sqrt(Pupae), x=Date_factor,fill=combined_separate))+ #, color=Flowgroup_coarse))+ 
  geom_boxplot()+
  theme(axis.title = element_text(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(), legend.position="blank")+
  scale_x_discrete(breaks=c("4/15/21","6/4/21","6/11/21","6/18/21","6/25/21","7/2/21","7/9/21","7/16/21","7/23/21","7/30/21","8/6/21","8/13/21","8/20/21","8/27/21","9/3/21","9/10/21","9/17/21","9/24/21"),
                   labels=c("4/15","6/4","6/11","6/18","6/25","7/2","7/9","7/16","7/23","7/30","8/6","8/13","8/20","8/27","9/3","9/10","9/17","9/24"))+
  xlab(NULL) + ylab("Pupae abundance (sqrt)")+ labs(fill="Basin type")+
  scale_fill_manual(values=c("combined"="grey","separate"="white"), 
                    breaks=c("combined","separate"), labels=c("Combined", "Separated"))
pupae_flowgroup_time


##########################################################################################
## S3

# decision tree / classification tree (mostly qualitative - get estimates, but not pvalues. get relative importance of variables)
# can do a lot of predictors all at once, because completely unparametric
# install.packages("tree")
library(tree)
metadata <- read.table("metadata.txt", sep="\t", header=TRUE)
cols.convert.factor <-c("date","date_code","cluster_clr_asv","cluster_philr","cluster_philr_genus",
                        "cluster_philr_family", "cluster_philr_phylum","philr_genus_split",
                        "philr_phylum_split","philr_phylum_split15","Pupae_pa","Moz_pa", "combined_separate")
metadata[,cols.convert.factor] <- lapply(metadata[,cols.convert.factor], factor)
class(metadata$cluster_philr_genus)
metadata.samples <- subset(metadata, sample_control=="sample")

testtree <- tree(data=metadata.samples, cluster_philr_genus ~ Temp + Cond + pH + DO + combined_separate)
testtree <- tree(data=metadata.samples, cluster_philr ~ date_code + Temp + Cond + pH + DO + Sal + combined_separate)
plot(testtree)
text(testtree)
# each var can be used multiple times. ex. cond is used a lot, so that var is important for determining the outcome variable
# length of vertical bar = how much information you gained by making the split
  # think of all the variability of all the response vars within a subset. take sample variance of all response outcomes, break it in different places and see which one reduces variability the most. like reverse of kmeans
  # the tree output: n=how many samples are in that group (branch), deviance = measure of variability, yval = most common category (like an "estimate"), the next 5 numbers are the relative proportions of each of hte 5 categories within the node
# nodes 2 and 3 is the first cut
# pruning a tree: define how much better the reduction in variance has to be before you make another split
?prune()
# at each leaf, whichever outcome category is the highest proportion is hte outcome you'd predict given those conditions
# do this with binary, outcome is proportino of 1s
# continuous, average of the group it ends up in
# algorithm, and criteria for makign the cut is different from clustering (which uses multivariate distance). this is just one cut that gives you the most differences





##########################################################################################
##S4

metadata_div <- subset(metadata, sample_control=="sample")
cols.convert.factor <-c("Pupae_pa", "cluster_philr_genus","cluster_philr_phylum")
metadata_div[,cols.convert.factor] <- lapply(metadata_div[,cols.convert.factor], factor)
metadata_div$cluster_philr_genus %>% class
metadata_div$date <- factor(metadata_div$date, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))
metadata_div$date


metadata_div$Flowgroup_coarse <- factor(metadata_div$Flowgroup_coarse, levels=c("Gibbons","MayfairCarlyle","Stratford","Donald-Banta","MinerEvanstonRammer"))
shannon.date <- ggplot(data=metadata_div, 
                        aes(x=date, y=shannon_unrar, fill=combined_separate)) + #fill=Flowgroup_coarse
  geom_boxplot() +
  scale_fill_manual(values=c("combined"="grey","separate"="white"),
                    breaks=c("combined","separate"), labels=c("Combined", "Separated"))+
  # scale_fill_manual(values=c("Gibbons"="dark grey","MayfairCarlyle"="grey", "Stratford"="light grey", "Donald-Banta"="white", "MinerEvanstonRammer"="white"))+
  scale_x_discrete(breaks=c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"),
                   labels=c("4/15","6/11","6/25","7/9","7/23","8/6","8/27","9/17"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"))+
  xlab(NULL)+ ylab("Shannon index")+labs(fill="Basin type")
shannon.date




##########################################################################################

## S5

ord.ps.all <- ord.philr.ps.all
plot.ord = phyloseq::plot_ordination(ps.all.rerooted, ord.ps.all, type="samples", color="date", shape="combined_separate")+
  geom_point(size=3)+
  # labs(color="Date", shape = "Pupae presence/absence")+
  scale_color_discrete(breaks=c(), labels=c())+
  # scale_shape(breaks=c("orig","2","5"), labels=c("Original","2d","5d"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill="white"))
print(plot.ord)








##########################################################################################


