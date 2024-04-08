library(dada2)
library(stringr)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(readxl)

theme_set(theme_classic())

#directory for fastq files
setwd("/Users/gagercoon/Desktop/CLB23/Data/fastq")

#list of fastq files
fnFs <- sort(list.files(".", pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(".", pattern = "_R2_001.fastq", full.names = TRUE))

#import names for fastq files
sample.names <- word(fnFs, 2, 3, sep = "_")

#plot quality profile
fnFs.quality.profile <- plotQualityProfile(fnFs[1:25])
fnRs.quality.profile <- plotQualityProfile(fnRs[1:25])

setwd("~/Desktop/CLB23/Figures")

tiff("fnFs.quality.profile.tiff", units="in", width=10, height=10, res=300)
fnFs.quality.profile
dev.off()

tiff("fnRs.quality.profile.tiff", units="in", width=10, height=10, res=300)
fnRs.quality.profile
dev.off()

setwd("/Users/gagercoon/Desktop/CLB23/Data/fastq/")

#filtered directory
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filter and trim samples
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                    trimLeft=c(19,20), trimRight=c(5,10),
                    maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)

#learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#make sure the files exist
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

#dada core sample interference algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#merge reads based on overlap
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#construct asv table
seqtab <- makeSequenceTable(mergers)

#view distribution, vast majority of reads are from length 251 to 256, just use those
table(nchar(getSequences(seqtab)))

#trim those reads as a result of non-specific priming
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 251:256]

#remove chimeras - great retention
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#create track of loss each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy
section.size <- 1500 #this could be any size, so long as it works for your computer

#split dataset into sections for better processing on my laptop
sections <- split(c(1:nrow(seqtab.nochim)),
                  sort(c(1:nrow(seqtab.nochim))%%ceiling(nrow(seqtab.nochim)/section.size)))

#assign taxonomy
taxa <- lapply(sections, function(x){return(assignTaxonomy(seqtab.nochim[x,],
                                                          refFasta="/Users/gagercoon/Desktop/CLB23/Data/fastq/silva_nr99_v138.1_train_set.fa.gz",
                                                          verbose = TRUE))})
taxa <- do.call(rbind, taxa)

#section off cloroplasts and mitochondria for removal
is.chloro <- taxa[,"Order"] %in% "Chloroplast"
is.mito <- taxa[,"Family"] %in% "Mitochondria"

#remove them, around 200 out of 11000 asvs
seqtab.final <- seqtab.nochim[, !(is.chloro | is.mito)]

taxa.final <- taxa[!(is.chloro | is.mito), ]

#create ASV references that are easier to use
asv.seqs <- colnames(seqtab.final)
asv.headers <- vector(dim(seqtab.final)[2], mode = "character")

for (i in 1:dim(seqtab.final)[2]) {
  asv.headers[i] <- paste(">ASV", i, sep ="_")
}

setwd("/Users/gagercoon/Desktop/CLB23/Code/16S files/")

#create files for exporting to phyloseq
asv_fasta <- c(rbind(asv.headers, asv.seqs))
write(asv_fasta, "ASV_noMC_CLB23_all.fa")

asv_tab <- t(seqtab.final)
row.names(asv_tab) <- sub(">", "", asv.headers)
write.table(asv_tab, "ASV_counts_noMC_CLB23_all.csv", sep = ",", quote=F, col.names = NA)

taxa.print <- taxa.final
rownames(taxa.print) <- NULL
row.names(taxa.print) <- sub(">", "", asv.headers)
write.table(taxa.print, "ASV_tax_noMC_CLB23_all.csv", sep = ",", quote=F, col.names = NA)

#import asv counts and add column names
CLB_ASVs <- read.table(file = 'ASV_counts_noMC_CLB23_all.csv', sep = ',', 
                       header = TRUE)
colnames(CLB_ASVs) <- sub("X", "", colnames(CLB_ASVs))
row.names(CLB_ASVs) <- CLB_ASVs[,1]

#remove samples that were bad in some way
CLB_ASVs <- CLB_ASVs[,-1]  #remove names
CLB_ASVs <- CLB_ASVs[,-19] #remove c3sec24
CLB_ASVs <- CLB_ASVs[,-13] #remove c2sec1

#import taxonomy
CLB_taxa <- read.table(file = 'ASV_tax_noMC_CLB23_all.csv', sep = ',', 
                       header = TRUE)
CLB_tax_row <- CLB_taxa[,1]
row.names(CLB_taxa) <- CLB_taxa[,1]
CLB_taxa <- CLB_taxa[,-1]

#import metadata about samples
CLB_sample <- read.table(file = 'CLB23_sample_file_all.csv', sep = ',', 
                         header = TRUE)
row.names(CLB_sample) <- CLB_sample[,1]
CLB_sample <- CLB_sample[,-1]
CLB_sample$Core <- as.character(CLB_sample$Core)

#create tables and merge into a phyloseq object
CLB_ASV_table <- otu_table(CLB_ASVs, taxa_are_row = TRUE)
CLB_taxa_table <- tax_table(CLB_taxa)
row.names(CLB_taxa_table) <- CLB_tax_row
colnames(CLB_taxa_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                              "Genus")
CLB_sample_table <- sample_data(CLB_sample)

#phyloseq object
CLB_ps <- phyloseq(CLB_ASV_table, CLB_sample_table, CLB_taxa_table)

#pruning data with less than 5 reads
CLB_ps = filter_taxa(CLB_ps, function(x) sum(x) > 0, TRUE)
CLB_ps <- prune_taxa(taxa_sums(CLB_ps) > 5, CLB_ps)

#remove contamination, list from Sheik et al. 2018 Frontiers in Microbiology
#####
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Acinetobacter") |  is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Pseudomonas") |  is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Abiotrophia") |  is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Achromobacter") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Acinetobacter") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Actinobacillus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Arcanobacterium") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Arcobacter") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Babesia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Bacillus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Bartonella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Bordetella") |  is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Borrelia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Brodetella") |  is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Brucella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Burkholderia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Campylobacter") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Capnocytophaga") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Chlamydia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Clostridium") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Comamonas") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Corynebacterium") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Coxiella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Cronobacter") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Deinococcus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Dermatophilus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Ehrlichia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Enterococcus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Erysipelothrix") |  is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Escherichia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Escherichia/Shigella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Flavobacterium") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Francisella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Gardnerella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Granulicatella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Haemophilus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Hafnia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Halomonas") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Helicobacter") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Klebsiella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Kocuria") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Lawsonia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Legionella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Leptospira") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Listeria") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Merkel_cell") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Micrococcus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Morganella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Mycobacterium") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Mycoplasma") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Neisseria") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Nocardia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Pasteurella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Photobacterium") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Plesiomonas") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Propionibacterium") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Proteus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Providencia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Pseudomonas") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Rhodococcus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Rickettsiae") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Roseomonas") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Rothia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Salmonella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Serratia") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Shewanella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Shigella") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Sphaerophorus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Staphylococcus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Stenotrophomonas") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Streptococcus") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Treponema") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Vibrio") | is.na(Genus))
CLB_ps <- subset_taxa(CLB_ps, (Genus != "Yersinia") | is.na(Genus))

CLB_ps <- subset_taxa(CLB_ps,  (Kingdom != "Eukaryota") | is.na(Kingdom))

#####

#archaea vs bacteria percentages
CLB_kingdoms <- aggregate_taxa(CLB_ps, "Kingdom")
CLB_kingdoms2 <- subset_taxa(CLB_kingdoms, Kingdom == "Archaea" | Kingdom == "Bacteria")
taxa_sums(CLB_kingdoms2)
#archaea 189500 total, bacteria 1029625 total, 1219125 combined, 15.5% archaea

#relative adundance
CLB_ps.ra <- transform_sample_counts(CLB_ps, function(x){x / sum(x)})

#beta diversity for supplementary material
ord.nmds.bray <- ordinate(CLB_ps, method="NMDS", distance = "bray")
ord.pcoa.bray <- ordinate(CLB_ps, method="PCoA", distance = "bray")
ord.cca.bray <- ordinate(CLB_ps, method="CCA", distance = "bray")

nmds_asv <- plot_ordination(CLB_ps, ord.nmds.bray, color ="Depth", title = "NMDS", shape = "Core") + 
  scale_colour_gradient2(high = "black", mid = "#2458d1", low = "#03cafc", midpoint = 15)
pcoa_asv <- plot_ordination(CLB_ps, ord.pcoa.bray, color ="Depth", title = "PCoA", shape = "Core") + 
  scale_colour_gradient2(high = "black", mid = "#2458d1", low = "#03cafc", midpoint = 15)
cca_asv <- plot_ordination(CLB_ps, ord.cca.bray, color ="Depth", title = "CCA", shape = "Core") + 
  scale_colour_gradient2(high = "black", mid = "#2458d1", low = "#03cafc", midpoint = 15)

beta_div <- ggarrange(nmds_asv, pcoa_asv, cca_asv, ncol = 3, nrow = 1, labels = c("A", "B", "C"), common.legend = TRUE)

setwd("~/Desktop/CLB23/Figures")

tiff("beta_diversity.tiff", units="in", width=8, height=4, res=300)
beta_div
dev.off()

#visualize all phyla for supplementary material
phyla <- tax_glom(CLB_ps, taxrank = "Phylum")
phyla <- transform_sample_counts(all_phylum, function(x){x / sum(x)})
phyla <- psmelt(all_phylum.ra)
phyla$Phylum <- as.character(data_all_phylum$Phylum)
phyla$Phylum[data_all_phylum$Abundance < 0.01] <- "< 1% Abundance"

colorlist = c("#666666", "#B53446", "#864F70", "#586A9A", "#3881AF", "#3E8E91", "#449C74", "#4AA956", "#58A057", "#6C856F",
              "#806B87", "#95519F", "#AF597D", "#CB6651", "#E77325", "#FF8301", "#FFA60F", "#FFC81D", "#FFEB2B", 
              "#F4EB31", "#DCBD2E", "#C4902B", "#AC6228", "#B55E45")

phyla_plot <- ggplot(phyla, aes(x = Depth, y = Abundance, fill = Phylum)) +
  geom_bar(aes(), stat = "identity", position = "stack") +
  facet_grid("Core") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow=15)) +
  coord_flip() +
  scale_x_reverse() +
  scale_fill_manual(values = colorlist)

tiff("phyla_plot.tiff", units="in", width=10, height=6, res=300)
phyla_plot
dev.off()

#visualize methanogens downcore
#####
#separate ps based on taxa and merge per core
mg.1 <- subset_taxa(CLB_ps.ra, Order == "Methanomassiliicoccales") 
mg.2 <- subset_taxa(CLB_ps.ra, Order == "Methanomicrobiales") 
mg.3 <- subset_taxa(CLB_ps.ra, Order == "Methanosarciniales")
mg.4 <- subset_taxa(CLB_ps.ra, Class == "ANME-1")
mg.5 <- subset_taxa(CLB_ps.ra, Order == "Methanofastidiosales")

mg <- merge_phyloseq(mg.1,mg.2,mg.3,mg.4,mg.5)

mg.c1 <- subset_samples(mg, Core == 1)
mg.c2 <- subset_samples(mg, Core == 2)
mg.c3 <- subset_samples(mg, Core == 3)

#agglomerate
mg.glom.c1 <- tax_glom(mg.c1, "Order")
mg.glom.c2 <- tax_glom(mg.c2, "Order")
mg.glom.c3 <- tax_glom(mg.c3, "Order")

#sending to table to manually change names and reimport
#changed random capitals in name, change sample names to just depth, and removed none of the rows
setwd("/Users/gagercoon/Desktop/CLB23/Data/")

write.table(mg.glom.c1 %>% psmelt() %>% 
              select(Order, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "mg.glom.c1.tsv", sep = "\t", quote = F, row.names = F, col.names = T) 
write.table(mg.glom.c2 %>% psmelt() %>% 
              select(Order, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "mg.glom.c2.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(mg.glom.c3 %>% psmelt() %>% 
              select(Order, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "mg.glom.c3.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#make sure to resave tsv files into xlsx before this
mg.glom.c1 <- read_excel("mg.glom.c1.xlsx")
mg.glom.c2 <- read_excel("mg.glom.c2.xlsx")
mg.glom.c3 <- read_excel("mg.glom.c3.xlsx")

#pivot data into rows
mg.glom.c1 <-
  pivot_longer(mg.glom.c1, cols = !Order, names_to = "Depth", values_to = "relative_abundance")
mg.glom.c2 <-
  pivot_longer(mg.glom.c2, cols = !Order, names_to = "Depth", values_to = "relative_abundance")
mg.glom.c3 <-
  pivot_longer(mg.glom.c3, cols = !Order, names_to = "Depth", values_to = "relative_abundance")

#order and arrange correctly for ggplot
mg.glom.c1$Depth <- as.numeric(mg.glom.c1$Depth)
mg.glom.c1 <- arrange(mg.glom.c1, desc(Depth))
mg.glom.c1$relative_abundance <- as.numeric(mg.glom.c1$relative_abundance)
mg.glom.c2$Depth <- as.numeric(mg.glom.c2$Depth)
mg.glom.c2 <- arrange(mg.glom.c2, desc(Depth))
mg.glom.c2$relative_abundance <- as.numeric(mg.glom.c2$relative_abundance)
mg.glom.c3$Depth <- as.numeric(mg.glom.c3$Depth)
mg.glom.c3 <- arrange(mg.glom.c3, desc(Depth))
mg.glom.c3$relative_abundance <- as.numeric(mg.glom.c3$relative_abundance)

#create plots
mg.glom.c1.plot <- ggplot(data = mg.glom.c1, aes(x = relative_abundance, y = Depth, group = Order, color = "Methanogens (Order)")) +
  geom_path(aes(color = Order), linetype = "dashed") +
  geom_point(aes(color = Order)) +
  labs(y = "Depth (cm)", x = "Relative abundance", color = "Methanogens (Order)") +
  scale_y_reverse() +
  scale_x_continuous(labels = scales::percent_format(scale = 100), limits = c(0,0.016)) +
  ggtitle("Core 1") +
  theme(legend.key.size = unit(.1, 'cm')) +
  expand_limits(y = 31)

mg.glom.c2.plot <- ggplot(data = mg.glom.c2, aes(x = relative_abundance, y = Depth, group = Order, color = "Methanogens (Order)")) +
  geom_path(aes(color = Order), linetype = "dashed") +
  geom_point(aes(color = Order)) +
  labs(y = "Depth (cm)", x = "Relative abundance", color = "Methanogens (Order)") +
  scale_y_reverse() +
  scale_x_continuous(labels = scales::percent_format(scale = 100), limits = c(0,0.016)) +
  ggtitle("Core 2") +
  theme(legend.key.size = unit(.1, 'cm')) +
  expand_limits(y = 31)

mg.glom.c3.plot <- ggplot(data = mg.glom.c3, aes(x = relative_abundance, y = Depth, group = Order, color = "Methanogens (Order)")) +
  geom_path(aes(color = Order), linetype = "dashed") +
  geom_point(aes(color = Order)) +
  labs(y = "Depth (cm)", x = "Relative abundance", color = "Methanogens (Order)") +
  scale_y_reverse() +
  scale_x_continuous(labels = scales::percent_format(scale = 100), limits = c(0,0.016)) +
  ggtitle("Core 3") +
  theme(legend.key.size = unit(.1, 'cm')) +
  expand_limits(y = 31)

#ggarrange
methanogen_plot <- ggarrange(mg.glom.c1.plot,
                             mg.glom.c2.plot,
                             mg.glom.c3.plot,
                             nrow = 1, ncol = 3, labels = c("A", "B", "C"), common.legend = TRUE)

setwd("~/Desktop/CLB23/Figures")

tiff("methanogen_plot.tiff", units="in", width=8, height=4, res=300)
methanogen_plot
dev.off()

#visualize SRB downcore
#####
#separate ps based on taxa and merge per core
srb.1 <- subset_taxa(CLB_ps.ra, Phylum == "Desulfobacterota")
srb.2 <- subset_taxa(CLB_ps.ra, Class == "Thermodesulfovibrionia")
srb.3 <- subset_taxa(CLB_ps.ra, Phylum == "Myxococcota")
srb.4 <- subset_taxa(CLB_ps.ra, Phylum == "SAR324 clade(Marine group B)")
srb.5 <- subset_taxa(CLB_ps.ra, Class == "Clostridia")

srb <- merge_phyloseq(srb.1,srb.2,srb.3,srb.4,srb.5)
srb.c1 <- subset_samples(srb, Core == 1)
srb.c2 <- subset_samples(srb, Core == 2)
srb.c3 <- subset_samples(srb, Core == 3)

#agglomerate
srb.glom.c1 <- tax_glom(srb.c1, "Class")
srb.glom.c2 <- tax_glom(srb.c2, "Class")
srb.glom.c3 <- tax_glom(srb.c3, "Class")

#sending to table to manually change names and reimport
#changed random capitals in name, change sample names to just depth, and removed polyangia
setwd("/Users/gagercoon/Desktop/CLB23/Data/")

write.table(srb.glom.c1 %>% psmelt() %>% 
              select(Class, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "srb.glom.c1.tsv", sep = "\t", quote = F, row.names = F, col.names = T) 
write.table(srb.glom.c2 %>% psmelt() %>% 
              select(Class, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "srb.glom.c2.tsv", sep = "\t", quote = F, row.names = F, col.names = T) 
write.table(srb.glom.c3 %>% psmelt() %>% 
              select(Class, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "srb.glom.c3.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#make sure to resave tsv files into xlsx before this
srb.glom.c1 <- read_excel("srb.glom.c1.xlsx")
srb.glom.c2 <- read_excel("srb.glom.c2.xlsx")
srb.glom.c3 <- read_excel("srb.glom.c3.xlsx")

#pivot data into rows
srb.glom.c1 <-
  pivot_longer(srb.glom.c1, cols = !Class, names_to = "Depth", values_to = "relative_abundance")
srb.glom.c2 <-
  pivot_longer(srb.glom.c2, cols = !Class, names_to = "Depth", values_to = "relative_abundance")
srb.glom.c3 <-
  pivot_longer(srb.glom.c3, cols = !Class, names_to = "Depth", values_to = "relative_abundance")

#order and arrange correctly for ggplot
srb.glom.c1$Depth <- as.numeric(srb.glom.c1$Depth)
srb.glom.c1 <- arrange(srb.glom.c1, desc(Depth))
srb.glom.c1$relative_abundance <- as.numeric(srb.glom.c1$relative_abundance)
srb.glom.c2$Depth <- as.numeric(srb.glom.c2$Depth)
srb.glom.c2 <- arrange(srb.glom.c2, desc(Depth))
srb.glom.c2$relative_abundance <- as.numeric(srb.glom.c2$relative_abundance)
srb.glom.c3$Depth <- as.numeric(srb.glom.c3$Depth)
srb.glom.c3 <- arrange(srb.glom.c3, desc(Depth))
srb.glom.c3$relative_abundance <- as.numeric(srb.glom.c3$relative_abundance)

#determine which are under 1% abundance to make manual scales for same colors between plots, view .agg df to see which are less than 1%
srb.glom.1.agg <- aggregate(relative_abundance ~ Class, srb.glom.c1, sum)
srb.glom.1.agg$Size[srb.glom.1.agg$relative_abundance < 0.01] <- "< 1% Abundance"
#total abundance is less than 1% for bacteriap25, Desulfobaccia, Desulfomonilia, Desulfovibrionia, Myxococcia, Thermodesulfovibrionia
srb.glom.2.agg <- aggregate(relative_abundance ~ Class, srb.glom.c2, sum)
srb.glom.2.agg$Size[srb.glom.2.agg$relative_abundance < 0.01] <- "< 1% Abundance"
#total abundance is less than 1% for bacteriap25, Desulfobaccia, Desulfomonilia, Desulfovibrionia, Dissulfuribacteria, Myxococcia, Syntrophia, Syntrophobacteria, Thermodesulfovibrionia
srb.glom.3.agg <- aggregate(relative_abundance ~ Class, srb.glom.c3, sum)
srb.glom.3.agg$Size[srb.glom.3.agg$relative_abundance < 0.01] <- "< 1% Abundance"
#total abundance is less than 1% for bacteriap25, Desulfobaccia, Desulfomonilia, Desulfovibrionia, Dissulfuribacteria, Myxococcia, Syntrophobacteria, Thermodesulfovibrionia

srb.glom.c1$Class <- as.character(srb.glom.c1$Class)
srb.glom.c1$Class[srb.glom.c1$Class %in% c("bacteriap25", "Desulfobaccia", "Desulfomonilia", "Desulfovibrionia", "Myxococcia", "Thermodesulfovibrionia")] <- "< 1% Abundance"
srb.glom.c1 <- aggregate(srb.glom.c1$relative_abundance, by = list(srb.glom.c1$Depth, srb.glom.c1$Class), sum)
colnames(srb.glom.c1)[1] ="Depth"
colnames(srb.glom.c1)[2] ="Class"
colnames(srb.glom.c1)[3] ="relative_abundance"
srb.glom.c2$Class <- as.character(srb.glom.c2$Class)
srb.glom.c2$Class[srb.glom.c2$Class %in% c("bacteriap25", "Desulfobaccia", "Desulfomonilia", "Desulfovibrionia", "Dissulfuribacteria", "Myxococcia", "Syntrophia", "Syntrophobacteria", "Thermodesulfovibrionia")] <- "< 1% Abundance"
srb.glom.c2 <- aggregate(srb.glom.c2$relative_abundance, by = list(srb.glom.c2$Depth, srb.glom.c2$Class), sum)
colnames(srb.glom.c2)[1] ="Depth"
colnames(srb.glom.c2)[2] ="Class"
colnames(srb.glom.c2)[3] ="relative_abundance"
srb.glom.c3$Class <- as.character(srb.glom.c3$Class)
srb.glom.c3$Class[srb.glom.c3$Class %in% c("bacteriap25", "Desulfobaccia", "Desulfomonilia", "Desulfovibrionia", "Dissulfuribacteria", "Myxococcia", "Syntrophobacteria", "Thermodesulfovibrionia")] <- "< 1% Abundance"
srb.glom.c3 <- aggregate(srb.glom.c3$relative_abundance, by = list(srb.glom.c3$Depth, srb.glom.c3$Class), sum)
colnames(srb.glom.c3)[1] ="Depth"
colnames(srb.glom.c3)[2] ="Class"
colnames(srb.glom.c3)[3] ="relative_abundance"

srb_colors.c1 <- c("#000000", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#BC80BD", "#CCEBC5")
srb_colors.c2 <- c("#000000", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462")
srb_colors.c3 <- c("#000000", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462", "#BC80BD")

#create plots
srb.glom.c1.plot <- ggplot(data = srb.glom.c1, aes(x = relative_abundance, y = Depth, group = Class,  color = "SRB (Class)")) +
  geom_path(aes(color = Class), linetype = "dashed") +
  geom_point(aes(color = Class)) +
  labs(y = "Depth (cm)", x = "Relative abundance", color = "SRB (Class)") +
  scale_y_reverse() +
  scale_x_continuous(labels = scales::percent_format(scale = 100), limits = c(0,0.18)) +
  scale_color_manual(values = srb_colors.c1) +
  ggtitle("Core 1") +
  theme(legend.key.size = unit(.1, 'cm')) +
  expand_limits(y = 31)

srb.glom.c2.plot <- ggplot(data = srb.glom.c2, aes(x = relative_abundance, y = Depth, group = Class,  color = "SRB (Class)")) +
  geom_path(aes(color = Class), linetype = "dashed") +
  geom_point(aes(color = Class)) +
  labs(y = "Depth (cm)", x = "Relative abundance", color = "SRB (Class)") +
  scale_y_reverse() +
  scale_x_continuous(labels = scales::percent_format(scale = 100), limits = c(0,0.18)) +
  scale_color_manual(values = srb_colors.c2) +
  ggtitle("Core 2") +
  theme(legend.key.size = unit(.1, 'cm')) +
  expand_limits(y = 31)

srb.glom.c3.plot <- ggplot(data = srb.glom.c3, aes(x = relative_abundance, y = Depth, group = Class,  color = "SRB (Class)")) +
  geom_path(aes(color = Class), linetype = "dashed") +
  geom_point(aes(color = Class)) +
  labs(y = "Depth (cm)", x = "Relative abundance", color = "SRB (Class)") +
  scale_y_reverse() +
  scale_x_continuous(labels = scales::percent_format(scale = 100), limits = c(0,0.18)) +
  scale_color_manual(values = srb_colors.c3) +
  ggtitle("Core 3") +
  theme(legend.key.size = unit(.1, 'cm')) +
  expand_limits(y = 31)

#ggarrange
srb_plot <- ggarrange(srb.glom.c1.plot,
                      srb.glom.c2.plot,
                      srb.glom.c3.plot,
                      nrow = 1, ncol = 3, labels = c("A", "B", "C"), common.legend = TRUE)

setwd("~/Desktop/CLB23/Figures")

tiff("srb_plot.tiff", units="in", width=8, height=4, res=300)
srb_plot
dev.off()