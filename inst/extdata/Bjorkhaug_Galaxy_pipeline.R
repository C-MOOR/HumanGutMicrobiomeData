
#PIPELINE HELPFUL LINKS
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex#generating-a-count-table
# https://benjjneb.github.io/dada2/tutorial.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/dada2/inst/doc/dada2-intro.html#filter-and-trim

#Raw data pulled from SRA PRJNA517050 into Galaxy
#Processed with https://usegalaxy.org/u/sayumiyork/w/dada2-amplicon-analysis-pipeline-for-paired-end-data
#Silva 132, no species
#Trim both forward and reverse to 250bp
#Origial metadata taken from SRA, groups made based on descriptions in original article
#52 samples in SRA, but one had empty metadata and was removed.


#load requisite packages
library(ggplot2)
library(phyloseq)
library(DESeq2)

seqtab.nochim<-t(as.matrix(read.csv("Alcohol_seqtab.csv", row.names = 1)))
taxa<-as.matrix(read.csv("Alchohol_taxtab.csv", row.names = 1))


# Import dataframe holding sample information
# Sample information must be made manually
samdf<-read.csv("Alcohol_metadata.csv")
head(samdf)
rownames(samdf) <- samdf$sample_id

# Construct phyloseq object (straightforward from dada2 outputs)
alc_counts <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Remove any ASVs that have 0 reads across all samples (in this case, no ASVs)
alc_counts <- prune_samples(sample_sums(alc_counts) >= 1, alc_counts)


# Normlaize data by proportions
alc <- transform_sample_counts(alc_counts, function(otu) otu/sum(otu))

# Save the phyloseq objects as .RDS files
saveRDS(alc, "alc.RDS")
saveRDS(alc_counts, "alc_counts.RDS")








############################################################################################################
library(ggplot2)
library(phyloseq)
library(DESeq2)

#Load RDS if needed here before testing
alc<-readRDS("alc.RDS")
alc_counts<-readRDS("alc_counts.RDS")

# Test alpha diversity plot
plot_richness(alc_counts, x="group", measures=c("Shannon", "Simpson"), color="sex")

# Test multidimensional analysis plot
#PCoA
ord.bray <- ordinate(alc, method="PCoA", distance="bray")
p <- plot_ordination(alc, ord.bray, color="group", shape = "sex", title="PCoA")
p + geom_point(size=3) + 
  theme_classic() +
  geom_text(aes(label=sample_id, alpha=0.6), check_overlap = TRUE) 
#NMDS
ord.bray <- ordinate(alc, method="NMDS", distance="bray")
p <- plot_ordination(alc, ord.bray, color="group", shape = "sex", title="NMDS")
p + geom_point(size=3) + 
  theme_classic() +
  geom_text(aes(label=sample_id, alpha=0.6), check_overlap = TRUE) 


# Test barplot
plot_bar(alc, x="sample_id", fill = "Class", title = "Bjorkhaug") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack")


#Subsetting based on a group
my_subset <- subset_taxa(alc, Class == "Bacteroidia")

#Plotting using the new subset
test<-plot_bar(my_subset, x= "sample_id", fill = "Order", title = "My title") +
  geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack")
test



# Differential abundance for a categorical variable
# STEP 1: Convert the phyloseq object to a DESeq2 object and tell R the experimental design 
my_alc_DESeq2 <- phyloseq_to_deseq2(alc_counts, design = ~ group)

# STEP 2: Select the groups to compare
my_comparison <-c("group", "control", "overconsumption")

# STEP 3: Run the differential abundance analysis at the chosen p-value
Significant_DEseq2_ASVs<-Differential_Abundance_AltGeo(alc_counts, my_alc_DESeq2, my_comparison, 0.05)

# STEP 4: Retrieve the list of ASVs with a significant difference in abundance between the chosen groups
Significant_DEseq2_ASVs

# STEP 5: Plot the results with your chosen x axis and legend
ggplot(Significant_DEseq2_ASVs, aes(x = Phylum, y=log2FoldChange, color= Class, shape = Order)) + geom_point(size=3, position = "jitter") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("my title here")



# Differential expression for continuous data in the 0-inflated dataset

# STEP 1: Convert the phyloseq object to a DESeq2 object and tell R the experimental design
my_DESeq2 <- phyloseq_to_deseq2(alc_counts, design = ~ age)

# STEP 2: Select the groups to compare
name_for_comparison <-"age"

# STEP 3: Run the differential abundance analysis at the chosen p-value
Significant_DEseq2_ASVs<-Differential_Abundance_ContinuousAltGeo(alc_counts, my_DESeq2, name_for_comparison, 0.05)

# STEP 4: Retrieve the list of ASVs with a significant difference in abundance between the chosen groups
Significant_DEseq2_ASVs

# STEP 5: Plot the results with your chosen x axis and legend
ggplot(Significant_DEseq2_ASVs, aes(x = Phylum, y=log2FoldChange, color= Class)) + geom_point(size=3, position = "jitter") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("my title here")



