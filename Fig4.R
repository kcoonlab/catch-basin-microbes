## Fig 4. Bacterial taxa significantly associated with different catch basin variables

# Identifying key taxa

library(phyloseq)
library(dplyr)
library(ggplot2)
library("ALDEx2")
library("beepr")
library("tibble")
library("dplyr")

ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
metadata <- read.table("metadata.txt", sep="\t", header=TRUE)
metadata <- metadata %>% mutate(pH_bin = if_else(pH <=7, "acid", "base"))
metadata$pH_bin
metadata.ps <- sample_data(metadata)
rownames(metadata.ps) <- metadata.ps$sample.id
ps.all <- phyloseq(otu_table(ps.all.rerooted),tax_table(ps.all.rerooted),metadata.ps, phy_tree(ps.all.rerooted))
cols.convert.factor <-c("Pupae_pa")
sample_data(ps.all)[,cols.convert.factor] <- lapply(sample_data(ps.all)[,cols.convert.factor], factor)
sample_data(ps.all)$date<- factor(sample_data(ps.all)$date, levels = c("4/15/21","6/11/21","6/25/21","7/9/21","7/23/21","8/6/21","8/27/21","9/17/21"))
sample_data(ps.all)$date %>% levels

ps.genus.rerooted <- readRDS("ps.genus.rerooted")
ps.family.rerooted <- readRDS("ps.family.rerooted")
ps.order.rerooted <- readRDS("ps.order.rerooted")
ps.class.rerooted <- readRDS("ps.class.rerooted")
ps.phylum.rerooted <- readRDS("ps.phylum.rerooted")

ps.genus.new <- phyloseq(otu_table(ps.genus.rerooted),tax_table(ps.genus.rerooted),metadata.ps, phy_tree(ps.genus.rerooted))
ps.phylum.new <- phyloseq(otu_table(ps.phylum.rerooted),tax_table(ps.phylum.rerooted),metadata.ps, phy_tree(ps.phylum.rerooted))
sample_data(ps.all)[1,]
ps.phylum.moz <- subset_samples(ps.phylum.new, !sample.id %in% c("A_41","A_43","A_50","A_52","A_57","A_60","A_66","A_69","A_71","A_78","A_83","A_88","F_47"))
sample_data(ps.phylum.moz)[1:15,1:4]

aldex.input.ps <- ps.phylum.new #ps.phylum.moz # find sig taxa for the philr biotypes, use rerooted phyloseq objects

sample_data(aldex.input.ps)$combined_separate
sample_data(aldex.input.ps)$Fail
sample_data(aldex.input.ps)$pH_bin
aldex.input.ps <- subset_samples(aldex.input.ps, Fail%in%c("0","1")) # bc aldex command not removing NAs even though specified to
aldex.input.ps <- subset_samples(aldex.input.ps, pH_bin%in%c("acid","base"))

aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(aldex.input.ps)), 
                               phyloseq::sample_data(aldex.input.ps)$pH_bin, #Fail
                               test="t", effect = TRUE, denom="iqlr", na.rm=TRUE)
beep(2)
head(aldex.results)
# plot effect sizes
ALDEx2::aldex.plot(aldex.results, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
# show output
aldex.results.sig <- aldex.results %>% rownames_to_column(var = "taxon") %>% filter(wi.eBH < 0.05) %>% arrange(effect, wi.eBH) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
head(aldex.results.sig)
dim(aldex.results.sig)
# make a list of taxa ordered by effect size
aldex.results.effectsize <- aldex.results %>% rownames_to_column(var = "taxon") %>% arrange(wi.eBH, effect) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
head(aldex.results.effectsize,20)
dim(aldex.results.effectsize)
# then join to taxonomy by ASV name
taxa_info <- data.frame(tax_table(aldex.input.ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
head(taxa_info)
aldex.results.effectsize.tax <- left_join(aldex.results.effectsize, taxa_info)
aldex.results.sig.tax <- left_join(aldex.results.sig, taxa_info)
head(aldex.results.effectsize.tax, 5)
head(aldex.results.sig.tax, 5)
# write.csv(aldex.results.sig.tax, "aldex2_results/sigtaxa_fail_genus.csv")
write.csv(aldex.results.sig.tax, "aldex2_results/sigtaxa_fail_phylum.csv")

sample_data(ps.phylum.new)[1:5,]
plot(sample_data(ps.phylum.new)$Pupae_pa, sample_data(ps.phylum.new)$Proteobacteria, na.omit=TRUE)


ps.effectvar <- subset_taxa(ps.genus.new, Genus %in% c("C39"))
ps.effectvar <- subset_taxa(ps.phylum.new, Phylum %in% c("Deinococcota"))
ps.effectvar
df.effectvar <- psmelt(ps.effectvar)
head(df.effectvar)
summary.effectvar <- df.effectvar %>% group_by(combined_separate, Phylum) %>% summarize(mean_abund = mean(Abundance, na.rm=TRUE))
summary.effectvar <- df.effectvar %>% group_by(Fail, Genus) %>% summarize(mean_abund = mean(Abundance, na.rm=TRUE))
head(summary.effectvar)



aldex.input.ps <- ps.genus.new # ps.genus.moz # find sig taxa for the philr biotypes, use rerooted phyloseq objects
aldex.input.ps <- subset_samples(aldex.input.ps, !(DO =="NA"))
sample_data(aldex.input.ps)$DO
aldex.cont <- aldex.clr(as.data.frame(otu_table(aldex.input.ps)),
                        sample_data(aldex.input.ps)$DO, #Cond pH
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(aldex.input.ps)$DO) #Cond pH
beep(2)
head(aldex.cont.results) # both phylum and genus have a lot of ties, problem calculating pvalues
aldex.cont.results %>% dim
# can't do effect size for continuous
# show output
aldex.cont.sig<- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
aldex.cont.sig %>% dim # kendall 21, spearman 21, pearson 13. do pearson to factor in relabund magnitude not just rank
# then join to taxonomy by ASV name
taxa_info <- data.frame(tax_table(aldex.input.ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
head(taxa_info)
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
head(aldex.cont.taxa, 5)
dim(aldex.cont.taxa)
write.csv(aldex.cont.taxa, "aldex2_results/sigtaxa_DO_genus.csv") #sigtaxa_cond_genus










# Heatmap

# https://www.biostars.org/p/134182/
library(phyloseq)
library(ape)
library(gplots)
ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
ps.all.rerooted
write.tree(phy_tree(ps.all.rerooted), "mafft-fasttree-output/exported-rerooted/tree.rerooted.nwk")

mytree <- read.tree("mafft-fasttree-output/exported-rerooted/tree.rerooted.nwk")
mytree_brlen <- compute.brlen(mytree, method="Grafen") #so that branches have all same length
hc <- as.hclust(mytree_brlen)

# tree of just the taxa that showed up in aldex.
# make a separate input, not pulling from phyloseq object of whole dataset?
taxa_names(ps.all.rerooted)

ps.phylum <- readRDS("ps_phylum_rerooted.rds")
ps.genus <- readRDS("ps_genus_rerooted.rds")
ps.genus






set.seed(123)                                                     # Set seed for reproducibility
data <- matrix(rnorm(100, 0, 10), nrow = 10, ncol = 10)           # Create example data
colnames(data) <- paste0("col", 1:10)                             # Column names
rownames(data) <- paste0("row", 1:10)                             # Row names
data
heatmap(data)

phylum_heatmap <- read.csv("aldex2_results/phylum_heatmap.csv", sep=",", header=TRUE)
rownames(phylum_heatmap) <- phylum_heatmap[,1]
phylum_heatmap
phylum_heatmap <- phylum_heatmap[,2:6]
# phylum_heatmap[is.na(phylum_heatmap)] <- 0
phylum_heatmap %>% class
phylum_heatmap <- as.matrix(phylum_heatmap)

hmpalette <- colorRampPalette(c("blue","grey","red"))
library(gplots)
heatmap.2(phylum_heatmap, scale="none", dendrogram="none",trace="none",
          Rowv=FALSE, Colv=FALSE,symm=TRUE,symkey=TRUE, key=TRUE,
          col=hmpalette, density.info=c("none"),
          cexRow=1, cexCol=1, srtCol=45, margins = c(10,10))
heatmap.2(phylum_heatmap)
dev.off()
phylum_heatmap$pH %>% class

data %>% class
phylum_heatmap %>% class

phylum_heatmap[,1]




mat = matrix( rnorm(25), 5, 5)
mat[c(1,6,8,11,15,20,22,24)] = NaN
mat
heatmap.2( mat,
           col = colorpanel(100,"red","yellow","green"),
           margins = c(12, 22),
           trace = "none", 
           xlab = "Comparison",
           lhei = c(2, 8),
           scale = c("none"),
           symbreaks = min(mat, na.rm=TRUE),
           na.color="blue",
           cexRow = 0.5, cexCol = 0.7,
           main = "DE genes", 
           dendrogram = "row", 
           Colv = FALSE )