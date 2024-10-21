set.seed(123)

## Basins with high average pupal abundance vs. those with highes frequency of pupae over time

WQ.data <- read.csv("input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data.by.basin <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
                   Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
                   Methoprene.fail.rate = sum(Methoprene.fail==1, na.rm=TRUE)/sum(Methoprene.fail>=0,na.rm=TRUE))
WQ.data.by.basin$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
WQ.data.by.basin$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Stratford","Stratford","Gibbons","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")

cor.test(WQ.data.by.basin$Pupae.abund.avg, WQ.data.by.basin$Pupae.prev, method=c("pearson"))

## Pupal abundances were also higher on average at the end of the season 

WQ.data <- read.csv("input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data$Sampling.date <- as.Date(WQ.data$Sampling.date, "%m/%d/%y")
WQ.data.pupaepresent <- WQ.data.pupaepresent %>% mutate(season = ifelse(Sampling.date < "2021-06-26", "early", "late"))
WQ.data.pupaepresent$Datefactor <- WQ.data.pupaepresent$Sampling.date %>% as.factor
WQ.data.pupaepresent = WQ.data.pupaepresent %>%
     arrange(Basin.id, Datefactor)                     
model = lme(sqrt(Pupae.abund) ~ season,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent)
anova(model)

## Season-wide treatment success did not differ between combined and separated basins or among basin flow groups

WQ.data <- read.csv("input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data.by.basin <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
                   Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
                   Methoprene.fail.rate = sum(Methoprene.fail==1, na.rm=TRUE)/sum(Methoprene.fail>=0,na.rm=TRUE))
WQ.data.by.basin$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
WQ.data.by.basin$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Stratford","Stratford","Gibbons","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")

anova.basintype <- aov(Methoprene.fail.rate ~ Basin.type, data=WQ.data.by.basin)
summary(anova.basintype)

anova.flowgroup <- aov(Methoprene.fail.rate ~ Basin.flowgroup, data=WQ.data.by.basin)
summary(anova.flowgroup)

## Alpha diversity in basins (as measured by Shannon’s H index) differed among sampling dates 

metadata <- read.table("input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
rm_anova <- anova_test(
  data = metadata,
  dv = Shannon,
  wid = Basin.id,
  within = Sampling.date,
  type = 3  # Specifies repeated measures design
)
get_anova_table(rm_anova)

## Separated basins further harbored bacterial communities characterized by higher Shannon index values than combined basins, irrespective of sampling date

metadata <- read.table("input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
metadata.Date1 <- metadata[metadata$Sampling.date == "4/15/21",]
metadata.Date2 <- metadata[metadata$Sampling.date == "6/11/21",]
metadata.Date3 <- metadata[metadata$Sampling.date == "6/25/21",]
metadata.Date4 <- metadata[metadata$Sampling.date == "7/9/21",]
metadata.Date5 <- metadata[metadata$Sampling.date == "7/23/21",]
metadata.Date6 <- metadata[metadata$Sampling.date == "8/6/21",]
metadata.Date7 <- metadata[metadata$Sampling.date == "8/27/21",]
metadata.Date8 <- metadata[metadata$Sampling.date == "9/17/21",]

kruskal.test(metadata.Date1$Shannon, metadata.Date1$Basin.type)
kruskal.test(metadata.Date2$Shannon, metadata.Date2$Basin.type)
kruskal.test(metadata.Date3$Shannon, metadata.Date3$Basin.type)
kruskal.test(metadata.Date4$Shannon, metadata.Date4$Basin.type)
kruskal.test(metadata.Date5$Shannon, metadata.Date5$Basin.type)
kruskal.test(metadata.Date6$Shannon, metadata.Date6$Basin.type)
kruskal.test(metadata.Date7$Shannon, metadata.Date7$Basin.type)
kruskal.test(metadata.Date8$Shannon, metadata.Date8$Basin.type)


metadata <- read.table("input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
rm_anova <- anova_test(
  data = metadata,
  dv = Shannon,
  wid = Basin.id,
  between = Basin.type,
  within = Sampling.date,
  type = 3  # Specifies repeated measures design
)
get_anova_table(rm_anova)

## In contrast, basins assigned to different flow groups harbored communities characterized by Shannon index values that varied more unpredictably over the season 

metadata <- read.table("input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
rm_anova <- anova_test(
  data = metadata,
  dv = Shannon,
  wid = Basin.id,
  between = Basin.flowgroup,
  within = Sampling.date,
  type = 3  # Specifies repeated measures design
)
get_anova_table(rm_anova)

## Biotype A” was more likely to be found in acidic basins, and “Biotype B” was more likely to be in basic basins

metadata <- read.table("input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
metadata$Datefactor <- metadata$Sampling.date %>% as.factor
metadata$acid_base[metadata$pH <= 7] <- "acid"
metadata$acid_base[metadata$pH > 7] <- "base"
metadata <- subset(metadata, !is.na(acid_base)==TRUE)
metadata <- subset(metadata, !is.na(Biotype)==TRUE)
metadata = metadata %>%
     arrange(Basin.id, Datefactor)                     
model = lme(Biotype ~ acid_base,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata)
anova(model)

metadata = metadata %>%
     arrange(Basin.id, Datefactor)                     
model = lme(Biotype ~ Basin.type,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = metadata)
anova(model)

##

ps.all <- readRDS("input-files/ps.all.rds")
ps.all.genus <- tax_glom(ps.all, "Genus", NArm = TRUE)
#minTotRelAbun = 1e-5
#x = taxa_sums(ps.all.genus)
#keepTaxa = (x / sum(x)) > minTotRelAbun
#prunedSet = prune_taxa(keepTaxa, ps.all.genus)

ps.all.genus <- subset_samples(ps.all.genus, !(Pupae_pa == "NA"))
sample_data(ps.all.genus)$Pupae_pa <- as.factor(sample_data(ps.all.genus)$Pupae_pa)

aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps.all.genus)), 
                               phyloseq::sample_data(ps.all.genus)$Pupae_pa,
                                   test="t", effect = TRUE, denom="iqlr")

aldex.results.sig <- aldex.results %>% rownames_to_column(var = "taxon") %>% filter(wi.eBH < 0.05) %>% arrange(effect, wi.eBH) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
aldex.results.effectsize <- aldex.results %>% rownames_to_column(var = "taxon") %>% arrange(wi.eBH, effect) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
taxa_info <- data.frame(tax_table(ps.all.genus))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.results.effectsize.tax <- left_join(aldex.results.effectsize, taxa_info)
aldex.results.sig.tax.pupaepa <- left_join(aldex.results.sig, taxa_info)
aldex.results.sig.tax.pupaepa

