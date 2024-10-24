set.seed(123)

library(dplyr)
library(nlme)
library(rstatix)
library(lme4)

## Basins with high average pupal abundance vs. those with highes frequency of pupae over time

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data.by.basin <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
                   Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE),
                   Methoprene.fail.rate = sum(Methoprene.fail==1, na.rm=TRUE)/sum(Methoprene.fail>=0,na.rm=TRUE))
WQ.data.by.basin$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
WQ.data.by.basin$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Stratford","Stratford","Gibbons","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")

cor.test(WQ.data.by.basin$Pupae.abund.avg, WQ.data.by.basin$Pupae.prev, method=c("pearson"))

## Pupal abundances were also generally higher later in the season 

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data$Sampling.date <- as.Date(WQ.data$Sampling.date, "%m/%d/%y")
WQ.data <- WQ.data %>% mutate(season = ifelse(Sampling.date < "2021-06-26", "early", "late"))
WQ.data.pupaepresent <- subset(WQ.data, !is.na(Pupae.pres)==TRUE)
WQ.data.pupaepresent$Datefactor <- WQ.data.pupaepresent$Sampling.date %>% as.factor
WQ.data.pupaepresent = WQ.data.pupaepresent %>%
     arrange(Basin.id, Datefactor)                     
model = lme(sqrt(Pupae.abund) ~ season,
              random = ~1|Basin.id, 
              correlation = corAR1(), 
              data = WQ.data.pupaepresent)
anova(model)

## Both the season-wide frequency and abundance of pupae did not significantly differ between basins as a function of basin type or flow group 

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
WQ.data.by.basin <- WQ.data %>%
  group_by(Basin.id) %>%
  dplyr::summarize(Pupae.abund.avg = mean(Pupae.abund, na.rm=TRUE),
                   Pupae.prev = sum(Pupae.pres>=1, na.rm=TRUE)/sum(Pupae.pres>=0, na.rm=TRUE))
WQ.data.by.basin$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
WQ.data.by.basin$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","Donald-Banta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Stratford","Stratford","Gibbons","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")

anova.basintype <- aov(Pupae.prev ~ Basin.type, data=WQ.data.by.basin)
summary(anova.basintype)
anova.basintype <- aov(Pupae.abund.avg ~ Basin.type, data=WQ.data.by.basin)
summary(anova.basintype)

anova.basintype <- aov(Pupae.prev ~ Basin.flowgroup, data=WQ.data.by.basin)
summary(anova.basintype)
anova.basintype <- aov(Pupae.abund.avg ~ Basin.flowgroup, data=WQ.data.by.basin)
summary(anova.basintype)

## Season-wide treatment success did not differ between combined and separated basins or among basin flow groups

WQ.data <- read.csv("catch-basin-microbes/input-files/WQ_Dips_2021_Final.csv",sep=",", header=TRUE)
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

## Combined and separated basin type did not affect biotype assignment, as most basins switched between biotypes at some point in the season

metadata <- read.table("catch-basin-microbes/input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
#metadata$Datefactor <- metadata$Sampling.date %>% as.factor
metadata <- subset(metadata, !is.na(Biotype)==TRUE)
#metadata = metadata %>%
#     arrange(Basin.id, Datefactor)

#model = lme(Biotype ~ Basin.type,
#              random = ~1|Basin.id, 
#              correlation = corAR1(), 
#              data = metadata)
#anova(model)

metadata <- metadata %>% mutate(Biotype.new = ifelse(Biotype == 1, 0, 1))
model = glmer(Biotype.new ~ Basin.type + (1 | Basin.id), data  = metadata, family = binomial)
model = glm(Biotype.new ~ Basin.type, data  = metadata, family = binomial)
summary(model)

## Alpha diversity in basins (as measured by Shannon’s H index) differed among sampling dates 

metadata <- read.table("catch-basin-microbes/input-files/metadata.txt", sep="\t", header=TRUE)
metadata <- subset(metadata, sample_control=="sample")
rm_anova <- anova_test(
  data = metadata,
  dv = Shannon,
  wid = Basin.id,
  within = Sampling.date,
  type = 3  # Specifies repeated measures design
)
get_anova_table(rm_anova)

## Both separated and combined basins, as well as basins assigned to different flow groups, harbored bacterial communities characterized by Shannon index values that varied somewhat unpredictably but were overall statistically similar to one another over the season

metadata <- read.table("catch-basin-microbes/input-files/metadata.txt", sep="\t", header=TRUE)
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

rm_anova <- anova_test(
  data = metadata,
  dv = Shannon,
  wid = Basin.id,
  between = Basin.flowgroup,
  within = Sampling.date,
  type = 3  # Specifies repeated measures design
)
get_anova_table(rm_anova)

## Several other genera in the Proteobacteria were significantly associated with pupal occurrence

ps.all <- readRDS("catch-basin-microbes/input-files/ps.all.rds")
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

##

ps.final <- readRDS("catch-basin-microbes/input-files/ps.final.rds")
ps.final.phylum <- tax_glom(ps.final, "Phylum", NArm = TRUE)

data <-
  ps.final.phylum %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()
                          
data %>%
  group_by(Sample) %>%
  arrange(-Abundance) %>%
  slice(1) %>%
  select(Phylum) %>%
  ungroup() %>%
  count(Phylum, name = "n_samples") %>%
  arrange(-n_samples)

data <-
  ps.final.genus %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()
                          
data %>%
  group_by(Sample) %>%
  arrange(-Abundance) %>%
  slice(1) %>%
  select(Genus) %>%
  ungroup() %>%
  count(Genus, name = "n_samples") %>%
  arrange(-n_samples)


y1 <- ps.final
y2 <- transform_sample_counts(y1, function(x) x/sum(x)) #get abundance in %
y3 <- psmelt(y2) # create dataframe from phyloseq object
y3$OTU[y3$Abundance < 0.01] <- "Taxa <1% abund."
y3$OTU[y3$Abundance >= 0.01] <- "Taxa >1% abund."                          

y4 <- y3 %>%
  group_by(Sample,OTU) %>%
  dplyr::summarize(Tot.abund = sum(Abundance, na.rm=TRUE))

y5 <- y3[!duplicated(y3$Sample),]
y6 <- merge(y5, y4, by = "Sample")

y6$Sampling.date <- as.Date(y6$Sampling.date, "%m/%d/%y")
y6 <- y6 %>% mutate(season = ifelse(Sampling.date < "2021-06-26", "early", "late"))
                              
y7 <- filter(y6, OTU.y == "Taxa <1% abund." & season == "early")
mean(y7$Tot.abund)     

y8 <- filter(y6, OTU.y == "Taxa <1% abund." & Biotype == 1)
mean(y8$Tot.abund)                                 

relabun.ps <- transform_sample_counts(ps.final.genus,function(x) x / sum(x))
relabun.ps <- subset_taxa(relabun.ps, Genus %in% c("C39"))

genus.df <- psmelt(relabun.ps)
genus.df2 <- genus.df
genus.df2$Sampling.date <- as.Date(genus.df2$Sampling.date, "%m/%d/%y")
genus.df2 <- genus.df2 %>% mutate(season = ifelse(Sampling.date < "2021-06-26", "early", "late"))

MySummary <- genus.df2 %>%
  group_by(season) %>%
  summarize(max_abund = max(Abundance, na.rm=TRUE)) 
head(MySummary)
print.data.frame(MySummary)                 
