library("phyloseq")
library("ALDEx2")
library("qiime2R")

set.seed(123)

library("ape")
library("gplots")
library("ggplot2")

library("beepr")
library("tibble")
library("dplyr")

## Fig 4. Bacterial taxa significantly associated with different catch basin variables

# Identifying key taxa

ps.all.rerooted <- readRDS("ps.all.rerooted.rds")
metadata <- read.table("metadata.txt", sep="\t", header=TRUE)
#metadata <- metadata %>% mutate(pH_bin = if_else(pH <=7, "acid", "base"))
metadata.ps <- sample_data(metadata)
rownames(metadata.ps) <- metadata.ps$sample.id
ps.all <- qza_to_phyloseq("table_analysis.qza", 
                          "rooted_tree.qza", 
                          "taxass_taxonomy_98.qza", 
                          "metadata.txt") 
sample_data(ps.all)$Pupae_pa <- as.factor(sample_data(ps.all)$Pupae_pa)
ps.phylum <- tax_glom(ps.all, "Phylum", NArm = TRUE)
ps.phylum.rerooted <- tax_glom(ps.all.rerooted, "Phylum", NArm = TRUE)
ps.phylum.new <- phyloseq(otu_table(ps.phylum.rerooted),tax_table(ps.phylum.rerooted),metadata.ps, phy_tree(ps.phylum.rerooted))
ps.phylum.new2 <- subset_samples(ps.phylum.new, !(date =="4/15/21"))

ps.phylum.new.Pupae_pa <- subset_samples(ps.phylum.new, !(Pupae_pa =="NA"))
aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps.phylum.new.Pupae_pa)), 
                               phyloseq::sample_data(ps.phylum.new.Pupae_pa)$Pupae_pa,
                                   test="t", effect = TRUE, denom="iqlr")
aldex.results.sig <- aldex.results %>% rownames_to_column(var = "taxon") %>% filter(wi.eBH < 0.05) %>% arrange(effect, wi.eBH) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
aldex.results.effectsize <- aldex.results %>% rownames_to_column(var = "taxon") %>% arrange(wi.eBH, effect) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
taxa_info <- data.frame(tax_table(ps.phylum.new.Pupae_pa))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.results.effectsize.tax <- left_join(aldex.results.effectsize, taxa_info)
aldex.results.sig.tax <- left_join(aldex.results.sig, taxa_info)
aldex.results.sig.tax

aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps.phylum.new)), 
                               phyloseq::sample_data(ps.phylum.new)$combined_separate,
                                   test="t", effect = TRUE, denom="iqlr")
aldex.results.sig <- aldex.results %>% rownames_to_column(var = "taxon") %>% filter(wi.eBH < 0.05) %>% arrange(effect, wi.eBH) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
aldex.results.effectsize <- aldex.results %>% rownames_to_column(var = "taxon") %>% arrange(wi.eBH, effect) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
taxa_info <- data.frame(tax_table(ps.phylum.new))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.results.effectsize.tax <- left_join(aldex.results.effectsize, taxa_info)
aldex.results.sig.tax <- left_join(aldex.results.sig, taxa_info)
aldex.results.sig.tax

ps.phylum.new.pH <- subset_samples(ps.phylum.new, !(pH =="NA"))
aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.phylum.new.pH)),
                        sample_data(ps.phylum.new.pH)$pH,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.phylum.new.pH)$pH)
aldex.cont.results
aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.phylum.new.pH))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa

ps.phylum.new.DO <- subset_samples(ps.phylum.new, !(DO =="NA"))
aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.phylum.new.DO)),
                        sample_data(ps.phylum.new.DO)$DO,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.phylum.new.DO)$DO)
aldex.cont.results
aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.phylum.new.DO))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa

ps.phylum.new.Cond <- subset_samples(ps.phylum.new, !(Cond =="NA"))
aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.phylum.new.Cond)),
                        sample_data(ps.phylum.new.Cond)$Cond,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.phylum.new.Cond)$Cond)
aldex.cont.results
aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.phylum.new.Cond))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa

ps.phylum.new.Sal <- subset_samples(ps.phylum.new, !(Sal =="NA"))
aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.phylum.new.Sal)),
                        sample_data(ps.phylum.new.Sal)$Sal,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.phylum.new.Sal)$Sal)
aldex.cont.results
aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.phylum.new.Sal))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa

ps.phylum.new.Pupae <- subset_samples(ps.phylum.new, !(Pupae =="NA"))
aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.phylum.new.Pupae)),
                        sample_data(ps.phylum.new.Pupae)$Pupae,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.phylum.new.Pupae)$Pupae)
aldex.cont.results
aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.phylum.new.Pupae))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa





ps.phylum <- readRDS("ps_phylum_rerooted.rds")
ps.effectvar <- subset_taxa(ps.phylum, Phylum %in% c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota", "Campilobacterota", "Chloroflexi", "Deinococcota", "Fibrobacterota", "Firmicutes", "Fusobacteriota", "Gemmatimonadota", "Myxococcota", "Patescibacteria", "Planctomycetota", "Proteobacteria", "Spirochaetota", "Synergistota", "Verrucomicrobiota"))
ps.effectvar.DO <- subset_samples(ps.effectvar, !(DO =="NA"))
aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.effectvar.DO)),
                        sample_data(ps.effectvar.DO)$DO,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.effectvar.DO)$DO)
aldex.cont.results
aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.effectvar.DO))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa

ps.effectvar.Cond <- subset_samples(ps.effectvar, !(Cond =="NA"))
aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.effectvar.Cond)),
                        sample_data(ps.effectvar.Cond)$Cond,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.effectvar.Cond)$Cond)
aldex.cont.results
aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.effectvar.Cond))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa










aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps.phylum)), 
                               phyloseq::sample_data(ps.phylum)$DO,
                                   test="t", effect = TRUE, denom="iqlr")




ps.phylum.pH <- subset_samples(ps.phylum, %in%c("0","1"))
aldex.results <- ALDEx2::aldex(data.frame(phyloseq:otu_table(ps.phylum)), 
                               data.frame(phyloseq::sample_data(ps.phylum)$DO),
                               test="t", effect = TRUE, denom="iqlr", na.rm=TRUE)

data.frame(otu_table(ps.phylum), sample_data(ps.phylum)$DO)


aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps.phylum)), 
                               phyloseq::sample_data(ps.phylum)$cluster_philr, # it's all comparing by the biotypes generated at the asv-level anyway. no need to change between taxonomic levels
                                   test="t", effect = TRUE, denom="iqlr")


ps.phylum.new <- phyloseq(otu_table(ps.phylum.rerooted),tax_table(ps.phylum.rerooted),metadata.ps, phy_tree(ps.phylum.rerooted))
sample_data(ps.all)[1,]

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
