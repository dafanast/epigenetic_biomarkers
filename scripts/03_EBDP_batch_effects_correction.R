# Relevant information from sessionInfo()
## R version 4.3.3 (2024-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux 9.4 (Plow)
## Packages: PCAtest_0.0.2	BEclear_2.18.0	BiocParallel_1.36.0	lubridate_1.9.3	forcats_1.0.0	stringr_1.5.1	dplyr_1.1.4	purrr_1.0.2	readr_2.1.5	tidyr_1.3.1	tibble_3.2.1	tidyverse_2.0.0	ggfortify_0.4.17	ggplot2_3.5.1

# Input data: a dataframe with correctly named unique CpGs as columns and samples as rows, without missing values and with age of samples ammended, output of EBDP_methylation_data_preparation.R script.
# Output data: a matrix with DNA methylation values corrected for batch effects. 

# 1. Prepare the environment
# Load required packages
library(ggfortify)
library(tidyverse)
#BiocManager::install("BEclear")
library(BEclear)
#library(devtools)
#devtools::install_github("arleyc/PCAtest")
library(PCAtest)

# Set working directory
#setwd("Joint_analysis_ML")

# 2. Prepare data
# INPUT: Load data from previous session in script 02_EBDP_data_preparation_for_clock
load("data/Methylation-data-preparation-output.RData") # Object named meth.age.df.na

# Example of how the file looks like
#head(meth.age.df.na[1:10,1:10])
#X1.15075 X1.15144  X1.15158 X1.15163 X1.15171 X1.53225 X1.53285 X1.53287 X1.53302 X1.53308
#1  91.66667 87.87879  95.08197 87.87879 90.90909 72.00000 55.55556 50.00000 61.61616 66.66667
#2  96.00000 78.85714  93.29268 76.57143 90.85366 59.73154 75.81699 67.97386 63.39869 53.59477
#3  93.63057 90.78947  96.62162 84.21053 83.21678 47.70642 74.31193 69.72477 43.11927 30.27523
#4  87.50000 82.60870  97.20670 92.39130 89.15663 61.65803 52.33161 47.66839 45.07772 47.66839
#5 100.00000 90.90909  66.66667 72.13115 66.66667 34.72222 45.83333 57.63889 54.16667 50.00000
#6  89.07563 94.95798 100.00000 89.07563 94.69027 41.12150 64.48598 87.85047 71.02804 58.87850

# Read samples information
samples.age <- read.table("data/samples-age.txt", sep="\t", stringsAsFactors = FALSE, quote="", header=TRUE) 
# Add the batch column to the dataframe
meth.batch <- meth.age.df.na %>% mutate(batch=samples.age$batch)

# 3. Perform PCAs where batch effect was evident
PC <- prcomp(meth.batch[,1:(ncol(meth.batch)-2)])
PCi<-data.frame(PC$x,Batch=as.factor(meth.age.df.na$batch))
var_explained <- PC$sdev^2/sum(PC$sdev^2)
#var_explained[1:5]

ggplot(PCi,aes(x=PC1,y=PC2,col=Batch))+
  geom_point(size=3,alpha=0.9, aes(col = Batch))+ 
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
  scale_colour_brewer(palette="Set1")+
  Style_format_theme # see corresponding script "Style_format_theme"
#ggsave("plots/GlobalMethylation-PCA-BatchEffects.tiff", width=20, height=18, units="cm", dpi=300)

# Test whether PCA are significant according to https://arleyc.github.io/PCAtest/articles/my-vignette.html
result<-PCAtest(meth.batch[,1:(ncol(meth.batch)-2)], varcorr=FALSE, counter=FALSE, plot=FALSE)

# 4. Detect which CpGs are affected by batch using package BEclear.
# Prepare samples dataframe for BEClear. Requires a data frame with two columns, the first column has to contain the sample numbers, the second column has to contain the corresponding batch number. Colnames have to be named as "sample_id" and "batch_id". 
batch.annot <- as.data.frame(bind_cols(sample_id=rownames(meth.batch), batch_id=meth.batch$batch))
# Prepare data for BEClear. Requires a matrix with rows as features and columns as samples that's why we need to transpose
meth.batch.rev <- t(meth.batch)
colnames(meth.batch.rev) <- rownames(batch.annot)
# Calcuate batch effects
batchEffect <- calcBatchEffects(
  data = meth.batch.rev, samples = batch.annot,
  adjusted = TRUE, method = "fdr")
# Optional: save the object to avoid re-calculating
#save(batchEffect, file="data/batcheffects-object.Rdata") 
#load("data/batcheffects-object.Rdata")

# Summarize median comparison and p-value calculation results
summary <- calcSummary(medians = batchEffect$med, pvalues = batchEffect$pval)
# Filter results for p-values
summary.pval <- summary %>% filter(pvalue<0.05)
# Calculate batch effect score
score <- calcScore(meth.batch.rev, batch.annot, summary, dir = getwd())
knitr::kable(score, caption = 'Batch scores')
# Set genes with batch effects to NA
cleared.data <- clearBEgenes(meth.batch.rev, batch.annot, summary) # Output: "2309408 values ( 24.4872320956159 % of the data) set to NA"
# Counts genes without batch effects
sum(complete.cases(cleared.data)) 

# 5. Eliminate CpGs affected by batch. Other solutions can be implemented at this stage, for example following the BEclear package, corrections for batch effects can be introduced or other ways of imputing "NA" values can be applied. Here, we chose to eliminate the affected CpGs from the dataset which is the most radical approach since we could afford losing several thousands of CpGs.
meth.corrected.batch <- na.omit(cleared.data) 
sum(complete.cases(meth.corrected.batch)) # Output: 44349
# Transpose dataframe
meth.corrected.batch.t <- t(meth.corrected.batch)

# 6. Visually check the corrected dataset by PCA again
# Calculate PCs
PC <- prcomp(meth.corrected.batch.t)
# Calculate variance explained by each PC
var_explained <- PC$sdev^2/sum(PC$sdev^2)
#var_explained[1:5]
# Prepare a dataframe for plotting
PCi<-data.frame(PC$x,Batch=as.factor(batch.annot$batch))
# Plot
ggplot(PCi,aes(x=PC1,y=PC2,col=Batch))+
  geom_point(size=3,alpha=0.9, aes(col = Batch))+ 
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
  scale_colour_brewer(palette="Set1")+
  Style_format_theme # see corresponding script "Style_format_theme"
# Save the plot
ggsave("plots/GlobalMethylation-PCA-BatchEffects-Corrected.tiff", width=20, height=18, units="cm", dpi=300)

# 7. Repeat evaluation of distribution of global methylation values using batch corrected table.
meth.corrected.batch.df <- as.data.frame(meth.corrected.batch.t)
# Convert to long format dataframe for downstream plotting 
meth.spread <- gather(meth.corrected.batch.df, key = cpgs, value = methylation, 1:(ncol(meth.corrected.batch.df)-1)) # Subtract the number of columns with biological factors. In this case we have only "age", that's why 2.
set.seed(123) 
meth.spr.ran <- sample_n(meth.spread, 20000) # Sample 5000 rows of data with dplyr to facilitate plotting
# Convert age to a factor
meth.spread$age <- factor(meth.spread$age)
#meth.spr.ran$age <- factor(meth.spr.ran$age)
# OUTPUT: Save the object with the corrected data
# save(meth.corrected.batch.t, file="data/meth-corrected-batch-effect.Rdata") 

# Plot global DNA methylation values as boxplots or violin plots per group to visualize check them
p <- ggplot(meth.spread, aes(x = age,y = methylation, group = age, fill=age)) + 
  #geom_violin(aes(fill=age),alpha=0.7) + 
  #geom_boxplot(aes(fill=age),width=0.1)+
  geom_boxplot(aes(fill=age),alpha=0.7)+
  #geom_jitter(alpha = 0.9, color = "black") + 
  labs(x="Otolith age (years)", y = "Methylation (%)") +
  scale_fill_manual(values=cols.age) + # see Appendix of this script
  # scale_color_manual(values=cols.age) +
  guides(fill = FALSE) +
  Style_format_theme # see corresponding script "Style_format_theme"
p
ggsave("plots/box-plots-globalDNAmethylation-otolithage-all.tiff", plot=p, width=20, height=18, units="cm", dpi=300)
#ggsave("violin-plots-globalDNAmethylation-all.tiff", plot=p, width=20, height=18, units="cm", dpi=300)

# Use statistical tests to check for differences in medians
pairwise.wilcox.test(meth.spread$methylation, meth.spread$age, p.adjust.method = "BH")
#RESULTS:
#0       1       2       3       4       5       6      
#1 < 2e-16 -       -       -       -       -       -      
#  2 < 2e-16 9.6e-09 -       -       -       -       -      
#  3 < 2e-16 < 2e-16 < 2e-16 -       -       -       -      
#  4 < 2e-16 0.00901 0.00311 < 2e-16 -       -       -      
#  5 < 2e-16 2.0e-06 < 2e-16 < 2e-16 1.3e-14 -       -      
#  6 < 2e-16 0.96834 0.00431 < 2e-16 0.04898 0.00023 -      
#  7 1.7e-12 < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16

# 8. Get the methylation values per age group after grouping per individual.
meth.corrected.batch.id <- meth.corrected.batch.df %>% mutate(id=rownames(meth.corrected.batch.df)) 
# Convert to long format dataframe for downstream plotting 
meth.spread.id <- gather(meth.corrected.batch.id, key = cpgs, value = methylation, 1:(ncol(meth.corrected.batch.df)-2)) # Subtract the number of 
mean.per.individual <- meth.spread.id %>% group_by(id,age) %>% summarise(mean.meth=mean(methylation), median.meth=median(methylation), n.sites=n())
mean.age <- mean.per.individual %>% group_by(age) %>% summarise(mean=mean(mean.meth), sd.mean=sd(mean.meth), median=median(median.meth), n=n())

# Appendix
# Define colors of age classes to make sure they are always the same across plots
cols.age <- c("#281b57",
              "#4ebf53",
              "#5d2192",
              "#005410",
              "#adb3ff",
              "#b16800",
              "#800046",
              "#c63b24")
