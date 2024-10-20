## Fig 4. Bacterial taxa significantly associated with different catch basin variables

library(phyloseq)
library(ALDEx2)
library(dplyr)
library(tibble)
library(gplots)

# Identifying key taxa (Pupal occurrence)

set.seed(123)
ps.all <- readRDS("input-files/ps.all.rds")
ps.all.phylum <- tax_glom(ps.all, "Phylum", NArm = TRUE)
ps.all.phylum.pupaepresence <- subset_samples(ps.all.phylum, !(Pupae_pa =="NA"))
sample_data(ps.all.phylum.pupaepresence)$Pupae_pa <- as.factor(sample_data(ps.all.phylum.pupaepresence)$Pupae_pa)

aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps.all.phylum.pupaepresence)), 
                               phyloseq::sample_data(ps.all.phylum.pupaepresence)$Pupae_pa,
                                   test="t", effect = TRUE, denom="iqlr")

aldex.results.sig <- aldex.results %>% rownames_to_column(var = "taxon") %>% filter(wi.eBH < 0.05) %>% arrange(effect, wi.eBH) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
aldex.results.effectsize <- aldex.results %>% rownames_to_column(var = "taxon") %>% arrange(wi.eBH, effect) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
taxa_info <- data.frame(tax_table(ps.all.phylum.pupaepresence))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.results.effectsize.tax <- left_join(aldex.results.effectsize, taxa_info)
aldex.results.sig.tax.pupaepa <- left_join(aldex.results.sig, taxa_info)
aldex.results.sig.tax.pupaepa

# Identifying key taxa (Pupal abundance)

ps.all.phylum.pupae <- subset_samples(ps.all.phylum, !(Pupae =="NA"))

aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.all.phylum.pupae)),
                        sample_data(ps.all.phylum.pupae)$Pupae,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.all.phylum.pupae)$Pupae)

aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.all.phylum.pupae))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa.pupae <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa.pupae

# Identifying key taxa (Basin type = Separated)

aldex.results <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps.all.phylum)), 
                               phyloseq::sample_data(ps.all.phylum)$combined_separate,
                                   test="t", effect = TRUE, denom="iqlr")
aldex.results.sig <- aldex.results %>% rownames_to_column(var = "taxon") %>% filter(wi.eBH < 0.05) %>% arrange(effect, wi.eBH) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
aldex.results.effectsize <- aldex.results %>% rownames_to_column(var = "taxon") %>% arrange(wi.eBH, effect) %>%
  dplyr::select(taxon, diff.btw, diff.win, effect, wi.ep, wi.eBH)
taxa_info <- data.frame(tax_table(ps.all.phylum))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.results.effectsize.tax <- left_join(aldex.results.effectsize, taxa_info)
aldex.results.sig.tax.basintype <- left_join(aldex.results.sig, taxa_info)
aldex.results.sig.tax.basintype

# Identifying key taxa (pH)

ps.all.phylum.pH <- subset_samples(ps.all.phylum, !(pH =="NA"))

aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.all.phylum.pH)),
                        sample_data(ps.all.phylum.pH)$pH,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.all.phylum.pH)$pH)

aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.all.phylum.pH))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa.pH <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa.pH

# Identifying key taxa (DO)

ps.all.phylum.DO <- subset_samples(ps.all.phylum, !(DO =="NA"))

aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.all.phylum.DO)),
                        sample_data(ps.all.phylum.DO)$DO,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.all.phylum.DO)$DO)

aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.all.phylum.DO))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa.DO <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa.DO

# Identifying key taxa (Cond)

ps.all.phylum.Cond <- subset_samples(ps.all.phylum, !(Cond =="NA"))

aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.all.phylum.Cond)),
                        sample_data(ps.all.phylum.Cond)$Cond,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.all.phylum.Cond)$Cond)

aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.all.phylum.Cond))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa.Cond <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa.Cond

# Identifying key taxa (Sal)

ps.all.phylum.Sal <- subset_samples(ps.all.phylum, !(Sal =="NA"))

aldex.cont <- aldex.clr(as.data.frame(otu_table(ps.all.phylum.Sal)),
                        sample_data(ps.all.phylum.Sal)$Sal,
                        mc.samples=128,
                        denom="all",
                        verbose=TRUE,
                        useMC=FALSE)
aldex.cont.results <- aldex.corr(aldex.cont, sample_data(ps.all.phylum.Sal)$Sal)

aldex.cont.sig <- aldex.cont.results %>% rownames_to_column(var = "taxon") %>% filter(pearson.ep < 0.05) %>%
  arrange(pearson.ep, pearson.ecor) %>% dplyr::select(taxon, pearson.ecor, pearson.ep, pearson.eBH)
taxa_info <- data.frame(tax_table(ps.all.phylum.Sal))
taxa_info <- taxa_info %>% rownames_to_column(var = "taxon")
aldex.cont.taxa.Sal <- left_join(aldex.cont.sig, taxa_info)
aldex.cont.taxa.Sal

aldex.results.sig.tax.pupaepa$Effectvar <- "Pupae present"
aldex.results.sig.tax.pupaepa <- aldex.results.sig.tax.pupaepa[,-c(1:3,5,7,9:13)]
colnames(aldex.results.sig.tax.pupaepa)[1] = "Effect"
colnames(aldex.results.sig.tax.pupaepa)[2] = "P-value"

aldex.cont.taxa.pupae$Effectvar <- "Pupal abundance"
aldex.cont.taxa.pupae <- aldex.cont.taxa.pupae[,-c(1,3,5,7:11)]
colnames(aldex.cont.taxa.pupae)[1] = "Effect"
colnames(aldex.cont.taxa.pupae)[2] = "P-value"

aldex.results.sig.tax.basintype$Effectvar <- "Basin type = Separated"
aldex.results.sig.tax.basintype <- aldex.results.sig.tax.basintype[,-c(1:3,5,7,9:13)]
colnames(aldex.results.sig.tax.basintype)[1] = "Effect"
colnames(aldex.results.sig.tax.basintype)[2] = "P-value"

aldex.cont.taxa.pH$Effectvar <- "pH"
aldex.cont.taxa.pH <- aldex.cont.taxa.pH[,-c(1,3,5,7:11)]
colnames(aldex.cont.taxa.pH)[1] = "Effect"
colnames(aldex.cont.taxa.pH)[2] = "P-value"

aldex.cont.taxa.Cond$Effectvar <- "Cond"
aldex.cont.taxa.Cond <- aldex.cont.taxa.Cond[,-c(1,3,5,7:11)]
colnames(aldex.cont.taxa.Cond)[1] = "Effect"
colnames(aldex.cont.taxa.Cond)[2] = "P-value"

aldex.cont.taxa.DO$Effectvar <- "DO"
aldex.cont.taxa.DO <- aldex.cont.taxa.DO[,-c(1,3,5,7:11)]
colnames(aldex.cont.taxa.DO)[1] = "Effect"
colnames(aldex.cont.taxa.DO)[2] = "P-value"

aldex.cont.taxa.Sal$Effectvar <- "Sal"
aldex.cont.taxa.Sal <- aldex.cont.taxa.Sal[,-c(1,3,5,7:11)]
colnames(aldex.cont.taxa.Sal)[1] = "Effect"
colnames(aldex.cont.taxa.Sal)[2] = "P-value"

phylum_heatmap <- rbind(aldex.results.sig.tax.pupaepa, aldex.cont.taxa.pupae, aldex.results.sig.tax.basintype, aldex.cont.taxa.pH, aldex.cont.taxa.Cond, aldex.cont.taxa.DO, aldex.cont.taxa.Sal)
phylum_heatmap <- phylum_heatmap[,-2]
phylum_heatmap <- subset(phylum_heatmap, Phylum!="unclassified.Bacteria") 
phylum_heatmap_wide <- reshape(data = phylum_heatmap,
                    idvar= "Phylum",
                    v.names= c("Effect"),
                    timevar= "Effectvar",
                    direction = "wide")
phylum_heatmap_wide <- phylum_heatmap_wide[order(phylum_heatmap_wide$Phylum), ]
rownames(phylum_heatmap_wide) <- phylum_heatmap_wide[,1]
phylum_heatmap_wide <- phylum_heatmap_wide[,-1]
phylum_heatmap_wide <- as.matrix(phylum_heatmap_wide)

hmpalette <- colorRampPalette(c("blue","grey","red"))
heatmap.2(phylum_heatmap_wide, scale="none", dendrogram="none",trace="none",
          Rowv=FALSE, Colv=FALSE,symm=TRUE,symkey=TRUE, key=TRUE,
          col=hmpalette, density.info=c("none"),
          cexRow=1, cexCol=1, srtCol=45, margins = c(10,10))
