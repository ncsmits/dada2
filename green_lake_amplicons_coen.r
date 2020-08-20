---
title: "Green Lake amplicons"
author: "Artemis Louyakis"
date: "03.06.2020"
output: html_notebook
---

```{r load packages}
Sys.time()
Sys.Date()
getwd()

library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(cowplot); packageVersion("cowplot")
library(reshape); packageVersion("reshape")
```
```{r housekeeping}
dir.create("plots")
dir.create("tables")
```

```{r read data}
Sys.time()
Sys.Date()

path_16s <- "16S/"
#path_18s <- "18S/"
list.files(path_16s)
#list.files(path_18s)
```

```{r plot 16s quality, fig.width=10, fig.height=5}
Sys.time()
Sys.Date()

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_16s <- sort(list.files(path_16s, pattern="_R1_", full.names = TRUE))
fnRs_16s <- sort(list.files(path_16s, pattern="_R2_", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME2019*_R1_*.fastq
## I renamed the negative control adding 2019 to file name to simplify this step
sample.names_16s <- sapply(strsplit(basename(fnFs_16s), "2019"), `[`, 1)
plotQualityProfile(fnFs_16s[1:16])
```

```{r plot 18s quality, fig.width=10, fig.height=5}
Sys.time()
Sys.Date()

# Forward and reverse fastq filenames have format: SAMPLENAME18s*_R1_*.fastq
fnFs_18s <- sort(list.files(path_18s, pattern="_R1_", full.names = TRUE))
fnRs_18s <- sort(list.files(path_18s, pattern="_R2_", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names_18s <- sapply(strsplit(basename(fnFs_18s), "18s"), `[`, 1)
plotQualityProfile(fnFs_18s[1:19])
```

```{r filter 16s reads}
Sys.time()
Sys.Date()
# Place filtered files in filtered/ subdirectory
dir.create("16S/filtered")
filt_path_16s <- file.path(path_16s, "filtered") 
filtFs_16s <- file.path(filt_path_16s, paste0(sample.names_16s, "_r1_filt.fastq"))
filtRs_16s <- file.path(filt_path_16s, paste0(sample.names_16s, "_r2_filt.fastq"))

length(fnFs_16s)
length(fnRs_16s)

out_16s <- filterAndTrim(fnFs_16s, filtFs_16s, fnRs_16s, filtRs_16s, maxN=0, maxEE=c(2,5), truncQ=2, 
                     rm.phix=TRUE, compress=FALSE, multithread=TRUE, matchIDs=TRUE) 
out_16s
errF_16s <- learnErrors(filtFs_16s, multithread=TRUE) # estimate error rates in forward sequences
errR_16s <- learnErrors(filtRs_16s, multithread=TRUE) # estimate error rates in reverse sequences

plotErrors(errF_16s, nominalQ=TRUE) # visualize estimated error rates
```

```{r filter 18s reads}
Sys.time()
Sys.Date()
# Place filtered files in filtered/ subdirectory
dir.create("18S/filtered")
filt_path_18s <- file.path(path_18s, "filtered") 
filtFs_18s <- file.path(filt_path_18s, paste0(sample.names_18s, "_r1_filt.fastq"))
filtRs_18s <- file.path(filt_path_18s, paste0(sample.names_18s, "_r2_filt.fastq"))

length(fnFs_18s)
length(fnRs_18s)

out_18s <- filterAndTrim(fnFs_18s, filtFs_18s, fnRs_18s, filtRs_18s, maxN=0, maxEE=c(2,5), truncQ=2, 
                     rm.phix=TRUE, compress=FALSE, multithread=TRUE, matchIDs=TRUE) 
out_18s
errF_18s <- learnErrors(filtFs_18s, multithread=TRUE) # estimate error rates in forward sequences
errR_18s <- learnErrors(filtRs_18s, multithread=TRUE) # estimate error rates in reverse sequences

plotErrors(errF_18s, nominalQ=TRUE) # visualize estimated error rates
```

```{r 16s assign taxa}
Sys.time()
Sys.Date()
# dereplicate filtered fastq files
derepFs_16s <- derepFastq(filtFs_16s, verbose=TRUE)
derepRs_16s <- derepFastq(filtRs_16s, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs_16s) <- sample.names_16s
names(derepRs_16s) <- sample.names_16s

dadaFs_16s <- dada(derepFs_16s, err=errF_16s, multithread=TRUE) #Infer the sequence variants in each sample
dadaRs_16s <- dada(derepRs_16s, err=errR_16s, multithread=TRUE)
dadaFs_16s[[1]]

#merge denoised forward and reverse reads
mergers_16s <- mergePairs(dadaFs_16s, derepFs_16s, dadaRs_16s, derepRs_16s, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers_16s[[1]])

#construct sequence table
seqtab_16s <- makeSequenceTable(mergers_16s)
dim(seqtab_16s)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_16s)))

#Remove chimeric sequences:
seqtab.nochim_16s <- removeBimeraDenovo(seqtab_16s, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim_16s)
sum(seqtab.nochim_16s)/sum(seqtab_16s)

#Track sequences through pipeline. See if there is one step that loses too many reads. 
getN_16s <- function(x) sum(getUniques(x))
track_16s <- cbind(out_16s, sapply(dadaFs_16s, getN_16s), sapply(mergers_16s, getN_16s), rowSums(seqtab_16s), rowSums(seqtab.nochim_16s))
colnames(track_16s) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track_16s) <- sample.names_16s
head(track_16s)

taxa_silva <- assignTaxonomy(seqtab.nochim_16s, "databases/silva_nr_v138_train_set.fa", multithread=TRUE)
taxa_silva <- addSpecies(taxa_silva, "databases/silva_species_assignment_v138.fa")

taxa_rdp <- assignTaxonomy(seqtab.nochim_16s, "databases/rdp_train_set_16.fa", multithread=TRUE)
taxa_rdp <- addSpecies(taxa_rdp, "databases/rdp_species_assignment_16.fa")

taxa_silva_print <- taxa_silva # Removing sequence rownames for display only
taxa_rdp_print <- taxa_rdp
rownames(taxa_silva_print) <- NULL
head(taxa_silva_print)
rownames(taxa_rdp_print) <- NULL
head(taxa_rdp_print)
```

```{r 18s assign taxa}
Sys.time()
Sys.Date()
# dereplicate filtered fastq files
derepFs_18s <- derepFastq(filtFs_18s, verbose=TRUE)
derepRs_18s <- derepFastq(filtRs_18s, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs_18s) <- sample.names_18s
names(derepRs_18s) <- sample.names_18s

dadaFs_18s <- dada(derepFs_18s, err=errF_18s, multithread=TRUE) #Infer the sequence variants in each sample
dadaRs_18s <- dada(derepRs_18s, err=errR_18s, multithread=TRUE)
dadaFs_18s[[1]]

#merge denoised forward and reverse reads
mergers_18s <- mergePairs(dadaFs_18s, derepFs_18s, dadaRs_18s, derepRs_18s, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers_18s[[1]])

#construct sequence table
seqtab_18s <- makeSequenceTable(mergers_18s)
dim(seqtab_18s)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_18s)))

#Remove chimeric sequences:
seqtab.nochim_18s <- removeBimeraDenovo(seqtab_18s, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim_18s)
sum(seqtab.nochim_18s)/sum(seqtab_18s)

#Track sequences through pipeline. See if there is one step that loses too many reads. 
getN_18s <- function(x) sum(getUniques(x))
track_18s <- cbind(out_18s, sapply(dadaFs_18s, getN_18s), sapply(mergers_18s, getN_18s), rowSums(seqtab_18s), rowSums(seqtab.nochim_18s))
colnames(track_18s) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track_18s) <- sample.names_18s
head(track_18s)

taxa_unite_f <- assignTaxonomy(seqtab.nochim_18s,
                               "databases/sh_general_release_s_04.02.2020/sh_general_release_dynamic_s_04.02.2020.fasta",
                               multithread=TRUE, tryRC = TRUE)
taxa_unite_all <- assignTaxonomy(seqtab.nochim_18s,
                                 "databases/sh_general_release_s_all_04.02.2020/sh_general_release_dynamic_s_all_04.02.2020.fasta",
                                 multithread=TRUE, tryRC = TRUE)
taxa_silva18 <- assignTaxonomy(seqtab.nochim_18s, "databases/silva_132.18s.99_rep_set.dada2.fa", multithread=TRUE)

taxa_unite_f_print <- taxa_unite_f # Removing sequence rownames for display only
taxa_unite_all_print <- taxa_unite_all
taxa_silva18_print <- taxa_silva18
rownames(taxa_unite_f_print) <- NULL
head(taxa_unite_f_print)
rownames(taxa_unite_all_print) <- NULL
head(taxa_unite_all_print)
rownames(taxa_silva18_print) <- NULL
head(taxa_silva18_print)
```

```{r 16s plots basic}
Sys.time()
Sys.Date()

theme_set(theme_bw())

map <- read.table("map_file_MM_C.txt", header=TRUE, sep = "\t")
#mapping file
sample_data(map)
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID
samples.out_16s <- rownames(seqtab.nochim_16s)

samdf_16s <- as.data.frame(map)

rownames(samdf_16s) <- samples.out_16s

ps_taxa_silva <- phyloseq(otu_table(seqtab.nochim_16s, taxa_are_rows=FALSE), 
               sample_data(samdf_16s), 
               tax_table(taxa_silva))
ps_taxa_silva

plot_richness(ps_taxa_silva, x="Dilution", measures=c("Shannon", "Simpson"), color="Dilution")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps_prop_taxa_silva <- transform_sample_counts(ps_taxa_silva, function(otu) otu/sum(otu))
# to prevent error due to zero values in negative controls:
ps_prop_taxa_silva <- prune_samples(sample_sums(ps_prop_taxa_silva) >= 0, ps_prop_taxa_silva) 
nmds_bray_silva <- ordinate(ps_prop_taxa_silva, method="NMDS", distance="bray")

plot_ordination(ps_prop_taxa_silva, nmds_bray_silva, color="Dilution", title="Bray NMDS")

top20_silva <- names(sort(taxa_sums(ps_taxa_silva), decreasing=TRUE))[1:20]
ps_top20_silva <- transform_sample_counts(ps_taxa_silva, function(OTU) OTU/sum(OTU))
ps_top20_silva <- prune_taxa(top20_silva, ps_top20_silva)
plot_bar(ps_top20_silva, x="SampleID", fill="Genus") + facet_wrap(~Dilution, scales="free_x")

ps_silva <- transform_sample_counts(ps_taxa_silva, function(OTU) OTU/sum(OTU))
plot_bar(ps_silva, x="SampleID", fill="Phylum") + facet_wrap(~Dilution, scales="free_x")

```

```{r 18s plots basic}
Sys.time()
Sys.Date()

theme_set(theme_bw())

samples.out_18s <- rownames(seqtab.nochim_18s)
sample <- c("APM1", "APM2", "APM4", "DIWater", "DIWater", "JULM5", "JULM6", "JULM7", "OCTF1", "OCTF2", "OCTF3", 
            "OCTM5", "OCTM6", "OCTM7", "OCTS1", "OCTS2", "OCTS3", "OCTW1", "OCTW3")
material <- c("Mat", "Mat", "Mat", "DIWater", "DIWater", "Mat", "Mat", "Mat", "Fluff", "Fluff", "Fluff", 
              "Mat", "Mat", "Mat", "Sediment", "Sediment", "Sediment", "Water", "Water")
replicate <- c("1", "2", "3", "DIWater", "DIWater", "1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2")
treatment <- c("APM", "APM", "APM", "DIWater", "DIWater", "JULM", "JULM", "JULM", "OCTF", "OCTF", "OCTF", 
               "OCTM", "OCTM", "OCTM", "OCTS", "OCTS", "OCTS", "OCTW", "OCTW")
month <- c("April", "April", "April", "DIWater", "DIWater", "July", "July", "July", "October", "October", "October", 
               "October", "October", "October", "October", "October", "October", "October", "October")
sample_type <- c("Sample", "Sample", "Sample", "Control", "Control", "Sample", "Sample", "Sample", "Sample", "Sample", "Sample", 
                 "Sample", "Sample", "Sample", "Sample", "Sample", "Sample", "Sample", "Sample")

samdf_18s <- data.frame(sample=sample, material=material, replicate=replicate, 
                        month=month, treatment=treatment, sample_type=sample_type)
rownames(samdf_18s) <- samples.out_18s

ps_taxa_silva18 <- phyloseq(otu_table(seqtab.nochim_18s, taxa_are_rows=FALSE), 
               sample_data(samdf_18s), 
               tax_table(taxa_silva18))
ps_taxa_silva18 <- prune_samples(sample_names(ps_taxa_silva18) != "Mock", ps_taxa_silva18) # Remove mock sample
ps_taxa_silva18

plot_richness(ps_taxa_silva18, x="material", measures=c("Shannon", "Simpson"), color="treatment")
plot_richness(ps_taxa_silva18, x="treatment", measures=c("Shannon", "Simpson"), color="material")
plot_richness(ps_taxa_silva18, x="material", measures=c("Shannon", "Simpson"), color="material")
## SimpsonInv diversity of 0 = infinite diversity, while 1 = no diversity

# Transform data to proportions as appropriate for Bray-Curtis distances
ps_prop_taxa_silva18 <- transform_sample_counts(ps_taxa_silva18, function(otu) otu/sum(otu))
# to prevent error due to zero values in negative controls:
ps_prop_taxa_silva18 <- prune_samples(sample_sums(ps_prop_taxa_silva18) >= 0, ps_prop_taxa_silva18) 
nmds_bray_silva18 <- ordinate(ps_prop_taxa_silva18, method="NMDS", distance="bray")

plot_ordination(ps_prop_taxa_silva18, nmds_bray_silva18, color="material", title="Bray NMDS")
plot_ordination(ps_prop_taxa_silva18, nmds_bray_silva18, color="treatment", title="Bray NMDS")
plot_ordination(ps_prop_taxa_silva18, nmds_bray_silva18, color="sample", title="Bray NMDS")

top20_silva18 <- names(sort(taxa_sums(ps_taxa_silva18), decreasing=TRUE))[1:20]
ps_top20_silva18 <- transform_sample_counts(ps_taxa_silva18, function(OTU) OTU/sum(OTU))
ps_top20_silva18 <- prune_taxa(top20_silva18, ps_top20_silva18)
plot_bar(ps_top20_silva18, x="material", fill="Genus") + facet_wrap(~treatment, scales="free_x")
plot_bar(ps_top20_silva18, x="treatment", fill="Genus") + facet_wrap(~material, scales="free_x")
ps_silva18 <- transform_sample_counts(ps_taxa_silva18, function(OTU) OTU/sum(OTU))
plot_bar(ps_silva18, x="material", fill="Phylum") + facet_wrap(~treatment, scales="free_x")
plot_bar(ps_silva18, x="treatment", fill="Phylum") + facet_wrap(~material, scales="free_x")
```

```{r decontaminate}
### decontam.txt
#load decontam
#BiocManager::install("decontam")

Sys.time()
Sys.Date()

#identify contaminants
phyloseq_object <- ps_taxa_silva
sample_data(phyloseq_object)$is.neg <- sample_data(phyloseq_object)$sample_type == "Control"
contamdf.prev <- isContaminant(phyloseq_object, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#make dataframe 
phyloseq_object.neg <- prune_samples(sample_data(phyloseq_object)$sample_type == "Control", phyloseq_object)
phyloseq_object.neg.presence <- transform_sample_counts(phyloseq_object.neg, function(abund) 1*(abund>0))
phyloseq_object.pos <- prune_samples(sample_data(phyloseq_object)$sample_type == "Sample", phyloseq_object)
phyloseq_object.pos.presence <- transform_sample_counts(phyloseq_object.pos, function(abund) 1*(abund>0))

df.pres <- data.frame(prevalence.pos=taxa_sums(phyloseq_object.pos.presence), 
                      prevalence.neg=taxa_sums(phyloseq_object.neg.presence),
                      contam.prev=contamdf.prev$contaminant)

#write a csv with all of the sequences shown as contaminant (TRUE) or not (FALSE)
write.csv(df.pres, file='contaminants_16s_silva.csv')

#identify contaminants
phyloseq_object <- ps_taxa_silva18
sample_data(phyloseq_object)$is.neg <- sample_data(phyloseq_object)$sample_type == "Control"
contamdf.prev <- isContaminant(phyloseq_object, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#make dataframe 
phyloseq_object.neg <- prune_samples(sample_data(phyloseq_object)$sample_type == "Control", phyloseq_object)
phyloseq_object.neg.presence <- transform_sample_counts(phyloseq_object.neg, function(abund) 1*(abund>0))
phyloseq_object.pos <- prune_samples(sample_data(phyloseq_object)$sample_type == "Sample", phyloseq_object)
phyloseq_object.pos.presence <- transform_sample_counts(phyloseq_object.pos, function(abund) 1*(abund>0))

df.pres <- data.frame(prevalence.pos=taxa_sums(phyloseq_object.pos.presence), 
                      prevalence.neg=taxa_sums(phyloseq_object.neg.presence),
                      contam.prev=contamdf.prev$contaminant)

#write a csv with all of the sequences shown as contaminant (TRUE) or not (FALSE)
write.csv(df.pres, file='contaminants_18s_silva.csv')

## this dataset resulted in no contaminating reads
```

```{r 16s alpha div}
## sample, material, replicate, treatment, sample_type

ps_16s <- phyloseq(otu_table(seqtab.nochim_16s, taxa_are_rows=FALSE), 
               sample_data(samdf_16s), 
               tax_table(taxa_silva))

est_rich_16s <- estimate_richness(ps_16s,  measures=c("Shannon", "Simpson"))
est_rich_16s$merge <- rownames(est_rich_16s)
samdf_16s$merge <- rownames(samdf_16s)
est_rich_16s <- merge(est_rich_16s, samdf_16s, by = "merge")
est_rich_16s <- melt(est_rich_16s)

ggplot(est_rich_16s, aes(x = month, y = value, color = material, shape = material)) + 
  geom_point(size = 2) +
  # scale_shape_manual(values=c(1, 2)) +
  xlab("") + ylab("Alpha Diversity Measure") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(legend.title=element_blank()) +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(size  = 9, angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(. ~ variable, scales="free_y")
# ggsave("plots/alpha_div_16s", dpi = 250)
```

```{r 18s alpha div}
## sample, material, replicate, treatment, sample_type

ps_18s <- phyloseq(otu_table(seqtab.nochim_18s, taxa_are_rows=FALSE), 
               sample_data(samdf_18s), 
               tax_table(taxa_silva18))

est_rich_18s <- estimate_richness(ps_18s,  measures=c("Shannon", "Simpson"))
est_rich_18s$merge <- rownames(est_rich_18s)
samdf_18s$merge <- rownames(samdf_18s)
est_rich_18s <- merge(est_rich_18s, samdf_18s, by = "merge")
est_rich_18s <- melt(est_rich_18s)

ggplot(est_rich_18s, aes(x = month, y = value, color = material, shape = material)) + 
  geom_point(size = 2) +
  # scale_shape_manual(values=c(1, 2)) +
  xlab("") + ylab("Alpha Diversity Measure") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(legend.title=element_blank()) +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(size  = 9, angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(. ~ variable, scales="free_y")
# ggsave("plots/alpha_div_16s", dpi = 250)
```

```{r normalize}

```


```{r cyanos}

```

```{r proteos}
# alpha / delta / epsilon
# anoxyphototrophs / sulfate reducers

```


Citations:
Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Opens external link in new windowNucl. Acids Res. 41 (D1): D590-D596.

Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Opens external link in new windowNucl. Acids Res. 42:D643-D648

Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2020): UNITE general FASTA release for Fungi. Version 04.02.2020. UNITE Community. https://doi.org/10.15156/BIO/786368

Benjamin Callahan. (2017). RDP taxonomic training data formatted for DADA2 (RDP trainset 16/release 11.5) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.801828
