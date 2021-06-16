# setting the working directory
setwd("/Users/evanpepper/Desktop/cuilab/cd-project/")
# listing out the necessary libraries
necessary_packages <- c("knitr", "phyloseq", "ggplot2", "plyr", "vegan", "reshape2", "gtable", "gridExtra", "grid", "Biostrings", "data.table", "ape", "tidyverse")
# importing the libraries
packages <- lapply(necessary_packages, library, character.only = TRUE)


# writing a function to process the data tables read into R one at a time, to then
# later be merged into one main data table of all the otu counts for all samples
otu_table_generater <- function(taxa_file) {
  # setting the working directory
  setwd("/Users/evanpepper/Desktop/cuilab/cd-project/otus/")
  #taxa_file <- "1_bracken_mpa.tsv"
  # reading in taxa counts file
  taxa_counts <- read.table(taxa_file)
  # renaming columns
  colnames(taxa_counts)[colnames(taxa_counts) == "V1"] <- "OTU"
  sample_number <- vapply(strsplit(taxa_file, "_"), `[`, 1, FUN.VALUE=character(1))
  sampleID <- paste("s", sample_number, sep="")
  colnames(taxa_counts)[colnames(taxa_counts) == "V2"] <- sampleID
  return(taxa_counts)
}

# resetting the current working directory
setwd("/Users/evanpepper/Desktop/cuilab/cd-project/otus/")
# making a list of all files in the directory
file_list <- list.files(path = ".")
# for each file in the directory...
for (file in file_list) {
  # if the abundances data table hasn't been made yet, make it
  if (!exists("abundances")) {
    abundances <- otu_table_generater(file)
  }
  # if the abundances data table does exist, make a temp data table to append to the original table
  if (exists("abundances")) {
    temp_abundances <- otu_table_generater(file)
    abundances <- merge(abundances, temp_abundances, all = TRUE)
    # then remove it to make a new temp data table for each iteration
    rm(temp_abundances)
  }
}

dat.unique_OTUs <- data.table(abundances$OTU)
dat.unique_OTUs$ID <- paste("OTU", seq.int(nrow(dat.unique_OTUs)), sep="")
colnames(dat.unique_OTUs)[colnames(dat.unique_OTUs) == "V1"] <- "OTU"

#####################################################################################################################

# creating OTU file
otus <- merge(dat.unique_OTUs, abundances, by = "OTU", all = TRUE)
otus[,OTU:=NULL]
dim(otus)
otus <- as.data.frame(otus)
rownames(otus) <- otus[,1]
otus <- otus[, 2:60]
# replacing all NA with 0
otus[is.na(otus)] <- 0
dim(otus)

#####################################################################################################################

# creating taxa file
# Create vector of taxonomic rank
v.Taxonomy = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

expanded_otus <- separate(dat.unique_OTUs, col = OTU, into = v.Taxonomy, sep = "\\|", remove = TRUE, fill = "right")
taxa <- cbind(expanded_otus$ID, expanded_otus$Kingdom, expanded_otus$Phylum, 
              expanded_otus$Class, expanded_otus$Order, expanded_otus$Family, 
              expanded_otus$Genus, expanded_otus$Species)
colnames(taxa) <- c("otuID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa <- as.data.table(taxa)

dt$PX = as.character(lapply(strsplit(as.character(dt$PREFIX), split="_"), "[", 1))

taxa <- as.data.frame(taxa)
rownames(taxa) <- taxa[,1]
taxa <- taxa[,2:ncol(taxa)]
v.levels <- c(1, 2, 3, 4, 5, 6, 7)
for (level in v.levels) {
  taxa[,level] <- as.character(lapply(strsplit(as.character(taxa[,level]), split="__"), "[", 2)) 
}
dim(taxa)
taxa <- as.matrix(taxa)

#####################################################################################################################

# creating meta file
setwd("/Users/evanpepper/Desktop/cuilab/cd-project/")
# reading in metadata file
meta <- read.table("cadmium-metadata.txt", header=T, sep=',')
# removing sample 5 because there was no data for this sample
meta <- meta[-c(5),]
meta$sampleID <- paste("s", meta$sampleID, sep="")
meta$timepoint <- paste("Week ", meta$timepoint, sep="")
meta <- as.data.table(meta)
meta <- meta[, c("bodysite","source"):=NULL]
meta <- as.data.frame(meta)
rownames(meta) <- meta[,1]
meta <- meta[,2:ncol(meta)]
dim(meta)

#####################################################################################################################

# building the phyloseq object
OTU = otu_table(otus, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
META = sample_data(meta)
CD_physeq = phyloseq(OTU, TAX, META)
CD_physeq

#####################################################################################################################

# performing analysis using https://www.nicholas-ollberding.com/post/introduction-to-phyloseq/

# nsamples(CD_physeq)
# sample_names(CD_physeq)
# sample_variables(CD_physeq)
# head(sample_data(CD_physeq))
# ntaxa(CD_physeq)
# head(taxa_names(CD_physeq))
# head(taxa_sums(CD_physeq))
# rank_names(CD_physeq)
# head(tax_table(CD_physeq))
# head(tax_table(CD_physeq)[, 2])
# table(tax_table(CD_physeq)[, 2])
# tax_table(CD_physeq)[1000:1005,]

#####

# relative abundance
relabund_CD_phyloseq <- transform_sample_counts(CD_physeq, function(x) x / sum(x) )
#otu_table(CD_physeq)
#otu_table(relabund_CD_phyloseq)

# filtered taxa, less than .1% of all OTUs
CD_physeq_filtered <- filter_taxa(relabund_CD_phyloseq, function(x) sum(x) > .001, TRUE)
#CD_physeq_filtered

#####################################################################################################################
### alpha diversity metrics and plots ###
#####################################################################################################################
head(estimate_richness(CD_physeq))

#Calculate and graph diversity measures for all 
alpha_meas = c("Shannon", "Observed", "InvSimpson", "Simpson")
#alpha_meas = c("Observed", "Shannon")

# plotting paired bar plots of alpha diversity between controls and treatments aggregating across timepoints
p <- plot_richness(CD_physeq, "treatment", measures=alpha_meas2)
p1 <- p + geom_boxplot(data=p$data, aes(x=treatment, y=value, fill=treatment), outlier.colour = "black", alpha= 0.5, colour = "black")
plot1 = p1 + theme(text = element_text(size = 20),axis.text.x = element_blank()) + labs(x = "")
plot1
ggsave(paste("figures/cd-alpha-div.png", sep = ""), plot1, height = 8.5, width = 16, units = "in")

for (metric in alpha_meas) {
  p <- plot_richness(CD_physeq, "sampletype", measure = metric, color = "timepoint") + geom_point(size = 4) + theme(text = element_text(size = 16)) + labs(x = "")
  ggsave(paste("figures/cd-alpha-div-", metric, ".png", sep = ""), p, height = 8.5, width = 11, units = "in")
}

# attempting to create time series plot by first creating data frame
alpha_metrics <- estimate_richness(CD_physeq)
alpha_time <- merge(alpha_metrics, meta, by="row.names")
alpha_meas2 <- c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher")

# multiple line plot for alpha diversity metrics across time points
for (metric in alpha_meas2) {
  time_series <- ggplot(alpha_time, aes_string(x = "timepoint", y = metric, group = "sampletype")) + 
    geom_line(aes(color = sampletype), size = 3, alpha = 0.8) +
    scale_color_manual(values = c("blue", "blue1", "blue2", "blue3", "blue4", "red", "red1", "red2", "red3", "red4")) +
    theme(text = element_text(size = 26)) + labs(x = "") +
    ggtitle(paste("Time Series of ", metric, " Alpha Diversity Across Time", sep = ""))
  ggsave(paste("figures/cd-time-series-", metric, "-alpha-div.png", sep = ""), time_series, height = 6.5, width = 22, units = "in")
}

# for future box plotting
# geom_boxplot(..., position="dodge")


# plotting phylum level taxa bar plots for all samples
bar1 <- plot_bar(CD_physeq_filtered, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  facet_grid(~timepoint, scales = "free", space = "free") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank(), text = element_text(size = 20))
bar1
ggsave(paste("figures/cd-bar-plot1.png", sep = ""), bar1, height = 8.5, width = 22, units = "in")

# getting unique treatment groups between controls and treatments
treatments <- as.vector(as.matrix(meta$sampletype))
uniques <- unique(treatments)
for (each_kind in uniques) {
  treat_phyloseq <- subset_samples(CD_physeq_filtered, sampletype == each_kind)
  bar2 <- plot_bar(treat_phyloseq, x = "timepoint", fill = "Phylum") +
    geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
    labs(x = "", y = "Relative Abundance\n") + ggtitle(paste(each_kind, " Relative Abundances", sep="")) +
    theme(panel.background = element_blank(), text = element_text(size = 24))
  ggsave(paste("figures/cd-bar-", each_kind, ".png", sep = ""), bar2, height = 8, width = 8, units = "in")
}

# subsetting the data for each treatment group
CD_physeq_treatment <- subset_samples(CD_physeq_filtered, treatment != "negative")
CD_physeq_no_treatment <- subset_samples(CD_physeq_filtered, treatment != "positive")

# treatment group only -- bar plots
bar3 <- plot_bar(CD_physeq_treatment, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  facet_grid(~timepoint, scales = "free", space = "free") +
  labs(x = "", y = "Relative Abundance\n") +
  ggtitle("Cadmium Treatment Taxa Bar Plot, Phylum Level") +
  theme(panel.background = element_blank(), text = element_text(size = 20))
bar3
ggsave(paste("figures/cd-bar-treatment.png", sep = ""), bar3, height = 4, width = 18, units = "in")

# control group only -- bar plots
bar4 <- plot_bar(CD_physeq_no_treatment, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  facet_grid(~timepoint, scales = "free", space = "free") +
  labs(x = "", y = "Relative Abundance\n") +
  ggtitle("Control Treatment Taxa Bar Plot, Phylum Level") +
  theme(panel.background = element_blank(), text = element_text(size = 20))
bar4
ggsave(paste("figures/cd-bar-no-treatment.png", sep = ""), bar4, height = 4, width = 18, units = "in")

#####################################################################################################################
### beta diversity metrics and plots ###
#####################################################################################################################

# PCoA plots
CD_bray_filtered <- ordinate(CD_physeq_filtered, "PCoA", "bray")
ord_1a <- plot_ordination(CD_physeq_filtered, CD_bray_filtered, type = "samples", color = "timepoint", shape="treatment")
ord_1b <- ord_1a + theme_bw() + ggtitle("PCoA -- Bray Curtis Dissimilarity \n Shaped by Treatment, Colored by Timepoint") 
plot2 = ord_1b + geom_point(size = 5) + theme_bw() + theme(panel.background = element_blank(), text = element_text(size = 24))
plot2
ggsave(paste("figures/cd-pcoa-bray.png", sep = ""), plot2, height = 8.5, width = 11, units = "in")

# NMDS plots
CD_NMDS_filtered <- ordinate(CD_physeq_filtered, "NMDS", "bray")
ord_2a <- plot_ordination(CD_physeq_filtered, CD_NMDS_filtered, type = "samples", color = "timepoint", shape="treatment")
ord_2b <- ord_2a + theme_bw() + ggtitle("NMDS -- Bray Curtis Dissimilarity \n Shaped by Treatment, Colored by Timepoint") 
plot3 = ord_2b + geom_point(size = 5) + theme_bw() + theme(panel.background = element_blank(), text = element_text(size = 24))
plot3
ggsave(paste("figures/cd-nmds-bray.png", sep = ""), plot3, height = 8.5, width = 11, units = "in")

# subsetting the data for only Firmicutes, the dominant phyla
CD_physeq_firmi <- subset_taxa(CD_physeq_filtered, Phylum == 'Firmicutes')
p <- plot_heatmap(CD_physeq_firmi, "NMDS", "bray", "sampletype", "Family")
p
plot_heatmap(CD_physeq_filtered, "Family")

#####################################################################################################################

# performing permanovas to identify statistical significant across groups
# Calculate bray curtis dissimilarity matrix 
CD_bray_distances <- phyloseq::distance(CD_physeq, method = "bray")
matrix<-as.matrix(dist(CD_bray_distances))
df <- data.frame(sample_data(CD_physeq))

# using adonis function from vegan to perform permanova
anova1 <- adonis2(formula = CD_bray_distances ~ timepoint + subject + treatment, data = df,
        permutations = 999, method = "bray", sqrt.dist = FALSE, add = FALSE, by = "margin")

anova2 <- adonis2(formula = CD_bray_distances ~ timepoint + treatment, data = df,
        permutations = 999, method = "bray", sqrt.dist = FALSE, add = FALSE, by = "margin")

