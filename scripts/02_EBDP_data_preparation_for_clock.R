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
setwd("/powerplant/workspace/cfndxa/Cod_clock/Joint_analysis_ML")

# 2. Prepare data
# INPUT: Load data object created by methylKit
load("/workspace/cfndxa/Cod_clock/Joint_analysis_ML/samples_new_annotation/meth-10cov-100000CpGs.Rdata")

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
samples.age <- read.table("samples/samples-age.txt", sep="\t", stringsAsFactors = FALSE, quote="", header=TRUE) 
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
#init = mice(meth.age.df, maxit=0)
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

# OUTPUT: this object can be saved as .Rdata and can be used in a new R session without need to run all previous steps. Also it can be saved as text file for use elsewhere if needed.
save(meth.age.df.na, file="Methylation-data-preparation-output.RData")
write.table(meth.age.df, file="NA-Samples-CpGs-Age.txt", sep="\t", dec=".", quote=FALSE, row.names=F, col.names=T)

# 4. Evaluate distribution of global methylation values. May use boxplots (here) or PCA (next script). This is a necessary quality check before proceeding with building of epigenetic clocks which may reveal unexpected patterns, like in this case batch effects (see next script).
# Convert to long format dataframe for downstream plotting 
meth.spread <- gather(meth.age.df.na, key = cpgs, value = methylation, 1:(ncol(meth.age.df.na)-1)) # Subtract the number of columns with biological factors. In this case we have only "age", that's why -1.
# OPTIONAL: Randomly sample fewer rows of the dataframe to speed up analysis that should be representative of general trends
set.seed(123) 
meth.spr.ran <- sample_n(meth.spread, 2000) # Sample 2000 rows of data with dplyr to facilitate plotting

# Plot global DNA methylation values as boxplots per group to visualize check them
p <- ggplot(meth.spread, aes(x = age,y = methylation, group = age, fill=age)) + geom_boxplot(alpha=0.7) + labs(x="Age (years)", y = "Methylation (%)", fill = "Age (years)") # Final plot including all DNA methylation
#p <- ggplot(meth.spr.ran, aes(x = age,y = methylation, group = age, fill=age)) + geom_boxplot(alpha=0.7) + labs(x="Age (years)", y = "Methylation (%)", fill = "Age (years)") # Quick plot for quick visualization
p + theme_classic() + scale_fill_viridis(discrete = FALSE) 
#ggsave("GlobalMethylation-withAge.tiff", width=10, height=8, units="cm", dpi=300)

# Use statistical tests to check for differences in medians
pairwise.wilcox.test(meth.spread$methylation, meth.spread$age, p.adjust.method = "BH")
