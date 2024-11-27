# Relevant information from sessionInfo()
## R version 4.3.3 (2024-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux 9.4 (Plow)
## Packages: viridis_0.6.5	viridisLite_0.4.2	ggsci_3.1.0	lubridate_1.9.3	forcats_1.0.0	stringr_1.5.1	dplyr_1.1.4	purrr_1.0.2	readr_2.1.5	tidyr_1.3.1	tibble_3.2.1	ggplot2_3.5.1	tidyverse_2.0.0	methylKit_1.28.0	GenomicRanges_1.54.1	GenomeInfoDb_1.38.8	IRanges_2.36.0	S4Vectors_0.40.2	BiocGenerics_0.48.1

# Input data:methylBase object from the methylKit package. Other objects containing percent methylation per CpG for all samples can also be used.
# Output data: a dataframe with correctly named unique CpGs as columns and samples as rows, without missing values and with age of samples amended. 

# 1. Prepare the environment
# Load required packages
library(methylKit)
require(tidyverse)
library(ggplot2)
library(ggsci)
library(viridis)

# Set working directory
#setwd("Joint_analysis_ML")

# 2. Prepare data
# INPUT: Load data object created by methylKit
load("data/meth-10cov-100000CpGs.Rdata")

# Example of how the file looks like
#head(meth)
#methylBase object with 6 rows
#--------------
#  chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2 coverage3 numCs3 numTs3 coverage4 numCs4
#1   1 15075 15075      -       132    121     11       175    168      7       157    147     10       184    161
#2   1 15144 15144      -       132    116     16       175    138     37       152    138     14       184    152
#numTs4 coverage5 numCs5 numTs5 coverage6 numCs6 numTs6 coverage7 numCs7 numTs7 coverage8 numCs8 numTs8 coverage9
#1     23        66     66      0       119    106     13       254    247      7       102     92     10        75
#2     32        66     60      6       119    113      6       254    247      7       102     92     10        75
#dim(meth)

# Check overall clustering using PCASamples function from methylKit
PCASamples(meth)

# Obtain percent methylation values
perc.meth=percMethylation(meth)
# Transpose dataframe to have samples as rows and CpGs as columns
perc.meth.df <- as.data.frame(t(perc.meth))
# Obtain the unique names of CpGs in the form of chromosome.start
meth.df <- as.data.frame(meth)
cpg.df <- bind_cols(chr=meth$chr, start=meth$start)
cpgs <- tidyr::unite(cpg.df, cpgs, chr:start, sep=".")
# Add the unique CpG names as column names in the dataframe
colnames(perc.meth.df) <- t(cpgs)
# Add the variables of interest. In this case we add "age" and "batch". Read the file containing this information in the following format with the samples ordered as in our dataframe (perc.meth.df):
samples.age <- read.table("data/samples-age.txt", sep="\t", stringsAsFactors = FALSE, quote="", header=TRUE) 
#head(samples.age)
#filename sample_file age sample batch
#1 102F.txt        102F   5   y5s1     1
#2 106F.txt        106F   1   y1s1     1
#3 107F.txt        107F   7   y7s1     1
#4 111F.txt        111F   3   y3s1     1
#5 118F.txt        118F   2   y2s1     1
#6 123F.txt        123F   1   y1s2     1
# Rownames of our dataframe (perc.meth.df) are in the same order
meth.age.df <- perc.meth.df %>% mutate(age=samples.age$age)

# Optional: Check the dimensions and format of the dataframe including the last couple of columns containing the added variables
dim(meth.age.df)
meth.age.df[1:10,c(1:10, (ncol(meth.age.df)-2):ncol(meth.age.df))]

# 3. Impute missing data. During the preparation of our object using methylKit we may have selected to "unite" samples based on: 
#a) only common values across all samples, the strictest option that produces no missing data, thus this step can be skipped.
#b) values common across a certain number of samples (for exampe methylKit::unite(norm, destrand=FALSE, min.per.group=48L)). In this case, there will be missing values in the dataset that cannot be handled in downstream steps. Thus we need to impute them.
# There are several ways of imputing missing values, here we present two of them. Always set.seed() for imputations.
#a) Method 1 using package “mice” (Multiple Imputation by Chained Equation)
#library(mice)
#set.seed(123)
#init = mice(meth.age.df)
#meth = init$method
#predM = init$predictorMatrix
#colnames(meth.age.df)
#predM[, c("age")]=0
#meth[c("age")]=""
#imputed = mice(meth.age.df, method=meth, predictorMatrix=predM, m=5)

#b) Method 2 using package “zoo” (Missing values replaced by the mean or other function of its group)
library(zoo)
set.seed(123)
meth.age.df.na <- na.aggregate(meth.age.df)

# OPTIONAL. OUTPUT: this object can be saved as .Rdata and can be used in a new R session without need to run all previous steps. Also it can be saved as text file for use elsewhere if needed.
# save(meth.age.df.na, file="data/Methylation-data-preparation-output.RData")

# 4. Evaluate distribution of global methylation values. May use boxplots (here) or PCA (next script). This is a necessary quality check before proceeding with building of epigenetic clocks which may reveal unexpected patterns, like in this case batch effects (see next script).
# Convert to long format dataframe for downstream plotting 
meth.spread <- gather(meth.age.df.na, key = cpgs, value = methylation, 1:(ncol(meth.age.df.na)-1)) # Subtract the number of columns with biological factors. In this case we have only "age", that's why -1.
# OPTIONAL: Randomly sample fewer rows of the dataframe to speed up analysis that should be representative of general trends
set.seed(123) 
meth.spr.ran <- sample_n(meth.spread, 2000) # Sample 2000 rows of data with dplyr to facilitate plotting

# Plot global DNA methylation values as boxplots per group to visualize check them
p <- ggplot(meth.spread, aes(x = age,y = methylation, group = age, fill=age)) + geom_boxplot(alpha=0.7) + labs(x="Age (years)", y = "Methylation (%)", fill = "Age (years)") + # Final plot including all DNA methylation
#p <- ggplot(meth.spr.ran, aes(x = age,y = methylation, group = age, fill=age)) + geom_boxplot(alpha=0.7) + labs(x="Age (years)", y = "Methylation (%)", fill = "Age (years)") # Quick plot for quick visualization
 scale_fill_manual(values=cols.age) + # see Appendix of this script
  # scale_color_manual(values=cols.age) +
  guides(fill = FALSE) +
  Style_format_theme # see corresponding script "Style_format_theme"
ggsave("plots/box-plots-globalDNAmethylation-otolithage-v1.tiff", plot=p, width=20, height=18, units="cm", dpi=300)

# Use statistical tests to check for differences in medians
pairwise.wilcox.test(meth.spread$methylation, meth.spread$age, p.adjust.method = "BH")

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
