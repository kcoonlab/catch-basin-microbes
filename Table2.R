## Table 2.  Effects of sampling date, water quality, mosquito productivity on microbiota diversity

library(phyloseq)
library(dplyr)
library(philr) #See: https://www.bioconductor.org/packages/release/bioc/html/philr.html
library(permute)
library(vegan)
library(tidyr)

## Table 2A  PERMANOVA tests using PhILR distances between samples

set.seed(123)
ps.final <- readRDS("input-files/ps.final.rds")
sample_data(ps.final)$sample.id <- sample_names(ps.final)

# running with just individual NA samples cut (not final instances of a UCB with an NA anywhere)
NA %in% sample_data(ps.final)$pH #row 39
NA %in% sample_data(ps.final)$Temp.C # row 39,129
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
#saveRDS(dist.most.philr, "dist_most_philr_Table2.rds")

perm_design <- how(within = Within(type = "series"),
    plots = Plots(strata = samples.most$Basin.id),
    nperm = 999)

adonis2(dist.most.philr ~ Sampling.date, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ pH, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ Temp.C, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ DO, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ Cond, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ Sal, data=samples.most, permutations=perm_design)

permutest(betadisper(dist.most.philr, samples.most$Sampling.date))
permutest(betadisper(dist.most.philr, samples.most$Temp.C))
permutest(betadisper(dist.most.philr, samples.most$Cond))
permutest(betadisper(dist.most.philr, samples.most$Sal))

samples.final <- as.data.frame(as.matrix(sample_data(ps.final)))
samples.pupaepresent <- samples.final %>% drop_na(Pupae.pres)
samples.pupaepresent <- samples.pupaepresent[-c(27,59,116),]
m <- as.matrix(dist.most.philr) 
dist.pupaepresent.philr <- as.dist(m[-c(1:12,121), -c(1:12,121)])

perm_design <- how(within = Within(type = "series"),
    plots = Plots(strata = samples.pupaepresent$Basin.id),
    nperm = 999)

adonis2(dist.pupaepresent.philr ~ Pupae.pres, data=samples.pupaepresent, permutations=perm_design)
adonis2(dist.pupaepresent.philr ~ Pupae.abund, data=samples.pupaepresent, permutations=perm_design)

permutest(betadisper(dist.pupaepresent.philr, samples.pupaepresent$Pupae.pres))
permutest(betadisper(dist.pupaepresent.philr, samples.pupaepresent$Pupae.abund))

## Table 2B  PERMANOVA tests using PhILR distances between samples aggregated by basin

ps.final.byUCB <- merge_samples(ps.final,"Basin.id")
sample_data(ps.final.byUCB)$Basin.flowgroup <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","DonaldBanta","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")
sample_data(ps.final.byUCB)$Basin.type <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
samples.final.byUCB <- as.data.frame(as.matrix(sample_data(ps.final.byUCB)))

ps.final.byUCB.pseudocount <- ps.final.byUCB
otu_table(ps.final.byUCB.pseudocount) <- t(otu_table(ps.final.byUCB) + 1)
ps.final.byUCB.philr <- philr(ps.final.byUCB.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
dist.final.byUCB.philr <- dist(ps.final.byUCB.philr, method="euclidean")
#saveRDS(dist.final.byUCB.philr, "dist_final_byUCB_philr_Table2.rds")

adonis2(dist.final.byUCB.philr ~ Basin.type, data=samples.final.byUCB)
adonis2(dist.final.byUCB.philr ~ Basin.flowgroup, data=samples.final.byUCB)

permutest(betadisper(dist.final.byUCB.philr, samples.final.byUCB$Basin.type))
