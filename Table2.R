???
ps.all.rerooted <- phyloseq(otu_table(ps.all), tax_table(ps.all), sample_data(ps.all), newlyrootedTree)

## Table 2.  Effects of sampling date, water quality, mosquito productivity on microbiota diversity

library(phyloseq)
library(dplyr)
library(philr) #See: https://www.bioconductor.org/packages/release/bioc/html/philr.html
library(permute)
library(vegan)
library(tidyr)

## Table 2A  PERMANOVA tests using PhILR distances between samples

set.seed(100)
ps.all <- readRDS("input-files/ps.all.rds")
sample_data(ps.all)$sample.id <- sample_names(ps.all)

# running with just individual NA samples cut (not all instances of a UCB with an NA anywhere)
NA %in% sample_data(ps.all)$pH #row 39
NA %in% sample_data(ps.all)$Temp[129] # row 39,129
NA %in% sample_data(ps.all)$DO[39] # row 39
NA %in% sample_data(ps.all)$Sal[71] # row 39, 71
NA %in% sample_data(ps.all)$Cond[39] # row 39
sample_data(ps.all)[129,] #sample.id: row 39=B_73, row 71=C_61, row129=F_56

# Take out the 3 individual samples that have NAs, and use that remaining set for the permanovas
ps.most <- subset_samples(ps.all, !sample.id %in% c("B_73","C_61","F_56"))
samples.most <- as.data.frame(as.matrix(sample_data(ps.most)))

ps.most.pseudocount <- ps.most
otu_table(ps.most.pseudocount) <- otu_table(ps.most) + 1
ps.most.philr <- philr(ps.most.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
dist.most.philr <- dist(ps.most.philr, method="euclidean")
#saveRDS(dist.most.philr, "dist_most_philr_Table2.rds")

perm_design <- how(within = Within(type = "series"),
    plots = Plots(strata = samples.most$UCB),
    nperm = 999)

adonis2(dist.most.philr ~ date, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ pH, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ Temp, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ DO, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ Sal, data=samples.most, permutations=perm_design)
adonis2(dist.most.philr ~ Cond, data=samples.most, permutations=perm_design)
permutest(betadisper(dist.most.philr, samples.most$date))
permutest(betadisper(dist.most.philr, samples.most$Temp))
permutest(betadisper(dist.most.philr, samples.most$Sal))
permutest(betadisper(dist.most.philr, samples.most$Cond))

samples.all <- as.data.frame(as.matrix(sample_data(ps.all)))
samples.pupaepresent <- samples.all %>% drop_na(Pupae_pa)
samples.pupaepresent <- samples.pupaepresent[-c(27,59,116),]
m <- as.matrix(dist.most.philr) 
dist.pupaepresent.philr <- as.dist(m[-c(1:12,121), -c(1:12,121)])

perm_design <- how(within = Within(type = "series"),
    plots = Plots(strata = samples.pupaepresent$UCB),
    nperm = 999)

adonis2(dist.pupaepresent.philr ~ Pupae_pa, data=samples.pupaepresent, permutations=perm_design)
adonis2(dist.pupaepresent.philr ~ Pupae, data=samples.pupaepresent, permutations=perm_design)
permutest(betadisper(dist.pupaepresent.philr, samples.pupaepresent$Pupae_pa))
permutest(betadisper(dist.pupaepresent.philr, samples.pupaepresent$Pupae))

***## Table 2B  PERMANOVA tests using PhILR distances between samples aggregated by basin

ps.all.byUCB <- merge_samples(ps.all,"UCB")
sample_data(ps.all.byUCB)$Flowgroup_coarse <- c("MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","DonaldBanta","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","MinerEvanstonRammer","DonaldBanta","MinerEvanstonRammer","MinerEvanstonRammer","Stratford","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","Stratford","Stratford","Stratford","Stratford","Gibbons","Gibbons","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle","MayfairCarlyle")
sample_data(ps.all.byUCB)$combined_separate <- c("separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","separate","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined","combined")
samples.all.byUCB <- as.data.frame(as.matrix(sample_data(ps.all.byUCB)))

ps.all.byUCB.pseudocount <- ps.all.byUCB
otu_table(ps.all.byUCB.pseudocount) <- otu_table(ps.all.byUCB) + 1
ps.all.byUCB.philr <- philr(ps.all.byUCB.pseudocount, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt') ###Error here has to do with need to root the tree (see above)!
dist.all.byUCB.philr <- dist(ps.all.byUCB.philr, method="euclidean")
#saveRDS(dist.all.byUCB.philr, "dist_all_byUCB_philr_Table2.rds")

[...]


adonis2(dist.UCBgroup.philr ~ combined_separate, data=samples.UCBgroup)
adonis2(dist.UCBgroup.philr ~ Flowgroup_coarse, data=samples.UCBgroup)
permutest(betadisper(dist.UCBgroup.philr, samples.UCBgroup$combined_separate))
permutest(betadisper(dist.UCBgroup.philr, samples.UCBgroup$Flowgroup_coarse))
