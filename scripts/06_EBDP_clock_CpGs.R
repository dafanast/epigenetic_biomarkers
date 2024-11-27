# CpGs of the epigenetic clock
# Input: Final model and full table with DNA methylation data and samples information
# Output: CpGs of the clock with clock's coefficients and correlations with age 

# 1. Prepare the environment
# Load required packages
library(WGCNA)
library(tidyverse)
library(ggplot2)
library(caret)

# Set working directory
#setwd("/powerplant/workspace/cfndxa/Cod_clock/Joint_analysis_ML/github")

# Read final model
load("data/FinalLASSOModel.Rdata") # lasso_model_all
# Read full methylation data with samples information
load("data/MethData_CompleteSampleInformation.Rdata") # meth.complete.full  # data file not available

# 2. Prepare data
# Prepare a dataframe containing the coefficients attributed to each CpG based on the final model
coef.lasso <- as.data.frame(as.matrix(coef(lasso_model_all$finalModel, lasso_model_all$bestTune$lambda)))
# Add CpG names, exclude the coefficient of the intercept (first row) and exclude CpGs with coefficient 0
nz.coef <- coef.lasso %>% add_column(cpgs=rownames(coef.lasso)) %>% filter(!row_number()==1) %>% filter(`1`!=0)
# Non-zero CpG names
cpgs.clock <- nz.coef$cpgs

# 3. Correlations of CpGs with age
# Select DNA methylation data of non-zero CpGs
cpgs.meth.clock <- select(meth.complete.full, all_of(cpgs.clock))
# Convert age to vector
x.age <- as.vector(meth.complete.full$age)
# Prepare a reduced dataframe including age
cpgs.meth.clock.age <- select(meth.complete.full, all_of(cpgs.clock), age)
# OUTPUT: Save the object with the reduced data
save(cpgs.meth.clock.age, file="data/CpGs-clock-age.Rdata")  # Here the object is called cpgs.meth.clock.age
#load("CpGs-clock-age.Rdata")

# Convert methylation data to matrix
y.cpgs <- data.matrix(cpgs.meth.clock) #And here it asks for the object cpgs.meth.clock, if it is the same as above please rename similarly?
# Calculate correlations
cor.age.pval <- corAndPvalue(x=x.age, y=y.cpgs, use="pairwise.complete.obs", alternative="two.sided")
# Prepare dataframe with important information
cor.pval.names <- bind_cols(cpgs=cpgs.clock, cor=as.numeric(cor.age.pval$cor), p=as.numeric(cor.age.pval$p), Z=as.numeric(cor.age.pval$Z), t=as.numeric(cor.age.pval$t))
# Add lasso coeffiecient 
cors.nz <- right_join(cor.pval.names, nz.coef, by="cpgs")
#head(cors.nz)
# A tibble: 6 × 6
#cpgs            cor        p     Z     t        s1
#<chr>         <dbl>    <dbl> <dbl> <dbl>     <dbl>
#  1 X1.12119654  -0.288 2.27e- 3 -3.08  3.13 -0.0370  
#2 X1.16253688   0.459 4.55e- 7  5.16  5.37  0.0676  
#3 X1.18209342  -0.307 1.09e- 3 -3.30  3.36 -0.00277 
#4 X11.6859520  -0.360 1.12e- 4 -3.92  4.01 -0.0256  
#5 X11.30230313  0.561 1.84e-10  6.59  7.04  0.000395
#6 X12.4375496  -0.203 3.33e- 2 -2.14  2.16 -0.0143 

#write.table(cors.nz, "data/Correlations-Lasso-Coefficients.txt", quote=F, sep="\t", dec=".",row.names=F, col.names=T)

# 4. CpGs correlated with age
# Positive correlation
cors.nz.pos <-  subset(cors.nz, cor>0)
# Negative correlation
cors.nz.neg <-  subset(cors.nz, cor<0)
# Add age
meth.pos.corr <- cpgs.meth.clock %>% select(cors.nz.pos$cpgs) %>% mutate(age=meth.complete.full$age)
dim(meth.pos.corr)
#110  34
meth.neg.corr <- cpgs.meth.clock %>% select(cors.nz.neg$cpgs) %>% mutate(age=meth.complete.full$age)
dim(meth.neg.corr)
#110  41

# 5. Plot clock CpGs 
# Convert dataframes into long format
meth.pos.long <- gather(meth.pos.corr, cpgs, methylation, 1:(dim(meth.pos.corr)[2]-1), factor_key=TRUE)
meth.neg.long <- gather(meth.neg.corr, cpgs, methylation, 1:(dim(meth.neg.corr)[2]-1), factor_key=TRUE)

# Positive correlation
p <- ggplot(meth.pos.long, aes(age, methylation)) +
  geom_point() +
  geom_smooth(method=lm) +
  facet_wrap(. ~ cpgs, nrow=4) +
  theme_classic()
p
ggsave("plots/Methylation-Positive_Age.tiff", width=28, height=24, units="cm", dpi=300)

# Negative correlation
p <- ggplot(meth.neg.long, aes(age, methylation)) +
  geom_point() +
  geom_smooth(method=lm) +
  facet_wrap(. ~ cpgs, nrow=4) +
  theme_classic()
p
ggsave("plots/Methylation-Negative_Age.tiff", width=28, height=24, units="cm", dpi=300)

# 6. Clock CpGs by biological factor
# Select DNA methylation data of non-zero CpGs
cpgs.meth.clock <- select(meth.complete.full, all_of(cpgs.clock))
# Check distribution using Principal Component Analysis (PCA) using only numerical values of the data
PC <- prcomp(cpgs.meth.clock[,1:ncol(cpgs.meth.clock)])
# Calculate variance explained by each component (for plotting)
var_explained <- PC$sdev^2/sum(PC$sdev^2)

# i. Distribution based on sex
# Prepare a reduced dataframe including sex
cpgs.meth.clock.sex <- select(meth.complete.full, all_of(cpgs.clock), sex) %>% mutate(sex=replace_na(sex, "UN"))
PCi<-data.frame(PC$x,BV=cpgs.meth.clock.sex$sex)
PCi$BV
ggplot(PCi,aes(x=PC1,y=PC2,col=BV))+
  geom_point(size=3,alpha=0.9, aes(col = BV))+ #Size and alpha just for fun
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
  scale_colour_brewer(palette="Set1")+
  theme_classic()
ggsave("plots/ClockCpGs-PCA-SexEffects_UN.tiff", width=9, height=7, units="cm", dpi=300)

# ii. Distribution based on maturity status
# Prepare a reduced dataframe including maturity status
cpgs.meth.clock.maturity <- select(meth.complete.full, all_of(cpgs.clock), maturity) %>% mutate(maturity=replace_na(maturity, "UN"))
PCi<-data.frame(PC$x,BV=as.factor(cpgs.meth.clock.maturity$maturity))
levels(as.factor(cpgs.meth.clock.maturity$maturity))
ggplot(PCi,aes(x=PC1,y=PC2,col=BV))+
  geom_point(size=3,alpha=0.9, aes(col = BV))+ #Size and alpha just for fun
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
  scale_colour_brewer(palette="Set1")+
  theme_classic()
ggsave("plots/ClockCpGs-PCA-MaturityEffects_UN.tiff", width=9, height=7, units="cm", dpi=300)

# iii. Distribution based on geographic location
# Prepare a reduced dataframe including geographic location
cpgs.meth.clock.geo <- select(meth.complete.full, all_of(cpgs.clock), geo_location)
PCi<-data.frame(PC$x,BV=cpgs.meth.clock.geo$geo_location)
ggplot(PCi,aes(x=PC1,y=PC2,col=BV))+
  geom_point(size=3,alpha=0.9, aes(col = BV))+ #Size and alpha just for fun
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
  scale_colour_brewer(palette="Set1")+
  theme_classic() +
  theme(legend.position="right")
#ggsave("ClockCpGs-PCA-Geo_location_Effects.tiff", width=8, height=6, units="cm", dpi=300)
ggsave("plots/ClockCpGs-PCA-Geo_location_Effects-smaller.tiff", width=9, height=7, units="cm", dpi=300)
