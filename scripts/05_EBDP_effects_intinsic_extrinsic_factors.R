# Check the effect of biological variables on ageing in general and on the performance of the final clock.
# A. Visually inspect the distribution based on biological variables to check for global effects.
# B. Test the performance of the final epigenetic clock depending on other biological variables. 

# 1. Prepare the environment
# Load required packages
library(ggfortify)
library(factoextra)
library(dplyr)
library(rstatix)

# Set working directory
#setwd("/powerplant/workspace/cfndxa/Cod_clock/Joint_analysis_ML/github")

# Read detailed samples information. The samples in this file are in different order than what we have used until now, thus, we will need to join two dataframes and not simply add columns.
samples.details <- read.table("data/samples-details.txt", sep="\t", stringsAsFactors = FALSE, quote="", header=TRUE) 
#head(samples.details)
#filename sample_id age sample batch station length weight gutted_weight  sex maturity geo_location
#1    1.txt         1   1  y1s16     2       3   20.5     82            73    M        1            A
#2    2.txt         2   1  y1s17     2       3   19.5     70            64    M        1            A
#3    4.txt         4   1  y1s18     2      12   25.5    153           139    M        1            A
#4    5.txt         5   1  y1s19     2      12   25.5    170           145    F        1            A
#5   6F.txt         6   0   y0s4     1      23   11.5     14            12 <NA>       NA            A
#6   7F.txt         7   5   y5s4     1      49   79.5   4420          3820    F        4            C

# For A. 2. Prepare data 
# Check the distribution of age according to intrinsic variables: sex and maturity status
# Change density plot line colors by groups
samples.density <- samples.details %>% mutate(sex=replace_na(sex, "UN"), maturity=replace_na(maturity, "UN"))
# Plot by sex
ggplot(samples.density, aes(x=age, color=sex, fill=sex)) +
  geom_density(alpha=0.3) +
  labs(x="Otolith age (years)", 
       y="Density") +
  scale_color_manual(values=c("red", "blue", "green3")) +
  scale_fill_manual(values=c("red", "blue", "green3")) +
  Style_format_theme # see corresponding script "Style_format_theme"
ggsave("plots/Density_age_distribution-sex.tiff", width=18, height=14, units="cm", dpi=300)

# Plot by maturity status
ggplot(samples.density, aes(x=age, color=maturity, fill=maturity)) +
  geom_density(alpha=0.3) +
  labs(x="Otolith age (years)", 
       y="") +
  scale_color_manual(values=c("red", "blue", "green3", "mediumpurple")) +
  scale_fill_manual(values=c("red", "blue", "green3", "mediumpurple")) +
  Style_format_theme
ggsave("plots/Density_age_distribution-maturity.tiff", width=18, height=14, units="cm", dpi=300)

# INPUT: Load data from previous session (prior to penalized regressions)
load("data/meth-corrected-batch-effect.Rdata")
# INPUT: Read samples information in order as the methylation data
samples.age <- read.table("data/samples-age.txt", sep="\t", stringsAsFactors = FALSE, quote="", header=TRUE) 
#head(samples.age)
#filename sample_file age sample batch
#1 102F.txt        102F   5   y5s1     1
#2 106F.txt        106F   1   y1s1     1
#3 107F.txt        107F   7   y7s1     1
#4 111F.txt        111F   3   y3s1     1
#5 118F.txt        118F   2   y2s1     1
#6 123F.txt        123F   1   y1s2     1
# Add sample column to the dataframe
meth.complete.samples <- as.data.frame(meth.corrected.batch.t) %>% mutate(sample=samples.age$sample)
# Join the dataframe containging methylation data and samples names with the detailed information per sample
meth.complete.full <- left_join(meth.complete.samples, samples.details[,!names(samples.details) %in% "age",], by="sample")
# save(meth.complete.full, file="data/MethData_CompleteSampleInformation.Rdata")
#head(meth.complete.full[1:5,c(1:5, 47470:47475)])
#X1.15075 X1.15144 X1.53225 X1.53285 X1.53287 length weight gutted_weight sex maturity geo_location
#1  91.66667 87.87879 72.00000 55.55556 50.00000   66.5   3488          2940   F        2            M
#2  96.00000 78.85714 59.73154 75.81699 67.97386   25.5    144           134   M        1            M
#3  93.63057 90.78947 47.70642 74.31193 69.72477   82.5   6000          5240   F        4            M
#4  87.50000 82.60870 61.65803 52.33161 47.66839   56.5   2072          1778   F        1            M
#5 100.00000 90.90909 34.72222 45.83333 57.63889   36.5    446           373   M        1            L

# 3. Check distribution using Principal Component Analysis (PCA) using only numerical values of the data
PC <- prcomp(meth.complete.full[,1:ncol(meth.corrected.batch.t)])
# Calculate variance explained by each component (for plotting)
var_explained <- PC$sdev^2/sum(PC$sdev^2)

# 4. Plots by biological variable:
# i. Distribution based on sex
PCi<-data.frame(PC$x,BV=meth.complete.full$sex) %>% mutate(BV=replace_na(BV, "UN"))
ggplot(PCi,aes(x=PC1,y=PC2,col=BV))+
  geom_point(size=3,alpha=0.9, aes(col = BV))+ #Size and alpha just for fun
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
  scale_colour_brewer(palette="Set1")+
  theme_classic()
ggsave("plots/GlobalMethylation-PCA-SexEffects_UN.tiff", width=9, height=7, units="cm", dpi=300)

# ii. Distribution based on maturity status. This is expected to change parallel to age since fish mature with age.
PCi<-data.frame(PC$x,BV=meth.complete.full$maturity) %>% mutate(BV=replace_na(BV, "UN"))
ggplot(PCi,aes(x=PC1,y=PC2,col=BV))+
  geom_point(size=3,alpha=0.9, aes(col = BV))+ #Size and alpha just for fun
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
  scale_colour_brewer(palette="Set1")+ #your colors here
  theme_classic()
ggsave("plots/GlobalMethylation-PCA-MaturityEffects_UN.tiff", width=9, height=7, units="cm", dpi=300)

# iii. Distribution based on geographic location.
PCi<-data.frame(PC$x,BV=as.factor(meth.complete.full$geo_location))
ggplot(PCi,aes(x=PC1,y=PC2,col=BV))+
  geom_point(size=3,alpha=0.9, aes(col = BV))+ #Size and alpha just for fun
  #labs(x=paste0("PC1 (", round(var_explained[1]*100,2),"%)"), 
   #    y=paste0("PC2 (", round(var_explained[2]*100,2),"%)")) +
  labs(x=paste0("Principal Component 1 (", round(var_explained[1]*100,2),"%)"), 
       y=paste0("Principal Component 2 (", round(var_explained[2]*100,2),"%)")) +
    scale_colour_brewer(palette="Set1")+ #your colors here
  theme_classic()
#ggsave("plots/GlobalMethylation-PCA-Geographic_location-large.tiff", width=8, height=6, units="cm", dpi=300)
ggsave("plots/GlobalMethylation-PCA-Geographic_location.tiff", width=9, height=7, units="cm", dpi=300)

# Calculate mean PC1 and plot it with latitude because it seems to be a relationship
PCi<-data.frame(PC$x,BV=as.factor(meth.complete.full$geo_location))
pc1.lat <- bind_cols(PC1=PCi$PC1,BV=PCi$BV)
pc1.lat.n <- pc1.lat %>%
  mutate(Latitude = case_when(
    BV == "A" ~ 54.5,
    BV == "C" ~ 57,
    BV == "L" ~ 59,
    BV == "M" ~ 60.5))
pc1.lat.mean <- pc1.lat.n %>% group_by(Latitude) %>% summarize(mean.pc1=mean(PC1), sd.mean=sd(PC1), n=n(), median.pc1=median(PC1)) %>% mutate(se.mean=sd.mean/n)
ggplot(pc1.lat.mean, aes(x=mean.pc1, y=Latitude, size=n)) +
    geom_point(shape=21, fill=c("red2", "blue2", "green3", "purple")) +
  labs(x="Mean PC1 (%)", y="Latitude") +
  scale_size(range = c(1, 20), name="n") + 
  Style_format_theme
ggsave("plots/Latitude_Mean_PC1_bubble_global_meth.tiff", width=20, height=16, units="cm", dpi=300)

# B. Test the performance of the final epigenetic clock depending on other biological variables. We test the effect on performance based on: test for differences in MAE between groups and check differences in correlations between predicted vs chronological age between groups.
# For B. 2. Prepare data 
# INPUT: vector on MAE and predicted age based on final model predictions from previous session
load("data/FinalModelMAE.Rdata") # called "mae"
load("data/FinalModelPredictedAge.Rdata") # called "predicted.age.inversed"
# Create a dataframe with columns: mae and biological variables of interest
predicted.mae.bv <- bind_cols(mae=mae, predicted=predicted.age.inversed, age=meth.complete.full$age, sex=meth.complete.full$sex, maturity=meth.complete.full$maturity, geo_location=meth.complete.full$geo_location)
#head(predicted.mae.bv)
# A tibble: 6 × 6
#mae predicted   age sex   maturity geo_location
#<dbl>     <dbl> <dbl> <chr>    <int> <chr>       
#  1 1.41        3.59     5 F            2 M           
#2 1.22        2.22     1 M            1 M           
#3 1.27        5.73     7 F            4 M           
#4 0.0300      3.03     3 F            1 M           
#5 0.0700      1.93     2 M            1 L           
#6 0.2         0.8      1 F            1 L     

# 3. Test for differences in MAE between groups. 
# 3.1 Check assumptions for tests: first assumption to check is normality of the mae distribution
shapiro.test(predicted.mae.bv$mae)
#Shapiro-Wilk normality test
#data:  predicted.mae.bv$mae
#W = 0.62081, p-value = 1.959e-15
# Assumption is not met, thus, we will perform non-parametric tests

# i. Sex
# Test differences between 2 groups
wilcox.test(predicted.mae.bv$mae ~ predicted.mae.bv$sex)
#Wilcoxon rank sum test with continuity correction
#data:  predicted.mae.bv$mae by predicted.mae.bv$sex
#W = 1361, p-value = 0.9666
#alternative hypothesis: true location shift is not equal to 0

# ii. Maturity status
# Test differences between 3 groups
kruskal.test(predicted.mae.bv$mae ~ predicted.mae.bv$maturity)
#Kruskal-Wallis rank sum test
#data:  predicted.mae.bv$mae by predicted.mae.bv$maturity
#Kruskal-Wallis chi-squared = 16.068, df = 2, p-value = 0.0003243
# Pairwise comparisons
dunn_test(predicted.mae.bv, mae ~ maturity, p.adjust.method = "bonferroni") 
#.y.   group1 group2    n1    n2 statistic        p    p.adj p.adj.signif
#* <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
#  1 mae   a      b         68    35     3.77  0.000161 0.000482 ***         
#  2 mae   a      c         68     2     1.71  0.0865   0.259    ns          
#3 mae   b      c         35     2     0.612 0.541    1        ns  

# iii. Geographic location
# Test differences between 4 groups
kruskal.test(predicted.mae.bv$mae ~ predicted.mae.bv$geo_location)
#Kruskal-Wallis rank sum test
#data:  predicted.mae.bv$mae by predicted.mae.bv$geo_location
#Kruskal-Wallis chi-squared = 5.4874, df = 3, p-value = 0.1394

# 3.2 Plot the error between groups by boxplots
# i. Sex
tiff("plots/MAE-Sex.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
boxplot(predicted.mae.bv$mae ~ predicted.mae.bv$sex, col=c("red2", "blue2"), ylim=c(0,0.6), xlab="Sex", ylab="Mean Absolute Error (years)", outline=F)
dev.off()

# ii. Maturity status
tiff("plots/MAE-Maturity.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
boxplot(predicted.mae.bv$mae ~ predicted.mae.bv$maturity, col=c("blue2", "red2", "green3"), xlab="Maturity Status", ylab="Mean Absolute Error (years)", ylim=c(0,1.4),outline=F)
dev.off()

# iii. Geographic location
tiff("plots/MAE-Geographic_location.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
boxplot(predicted.mae.bv$mae ~ predicted.mae.bv$geo_location, col=c("red2", "blue2", "green3", "purple"), xlab="Geographic location", ylab="Mean Absolute Error (years)", outline=F)
dev.off()

# 4. Check differences in correlation of predicted vs chronological age between groups. 
# 4.1. Perform separate correlations
# i. Sex
# Correlations per sex: males
cor.test(predicted.mae.bv[predicted.mae.bv$sex=="M",]$age, predicted.mae.bv[predicted.mae.bv$sex=="M",]$predicted)
#cor 
#0.9583422

# Correlations per sex: females
cor.test(predicted.mae.bv[predicted.mae.bv$sex=="F",]$age, predicted.mae.bv[predicted.mae.bv$sex=="F",]$predicted)
#cor 
#0.9596274

# ii. Maturity status
# Correlations per maturity: 1
cor.test(predicted.mae.bv[predicted.mae.bv$maturity=="a",]$age, predicted.mae.bv[predicted.mae.bv$maturity=="a",]$predicted)
#cor 
#0.9471186

# Correlations per maturity: 2
cor.test(predicted.mae.bv[predicted.mae.bv$maturity=="b",]$age, predicted.mae.bv[predicted.mae.bv$maturity=="b",]$predicted)
#cor 
#0.8594402

# Correlations per maturity: 4
cor.test(predicted.mae.bv[predicted.mae.bv$maturity=="c",]$age, predicted.mae.bv[predicted.mae.bv$maturity=="c",]$predicted)
# Not enough observations for this test (n=2)

# iii. Geographic location
# Correlations per geographic location: A
cor.test(predicted.mae.bv[predicted.mae.bv$geo_location=="A",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="A",]$predicted)
#cor 
#0.944599

# Correlations per geographic location: C
cor.test(predicted.mae.bv[predicted.mae.bv$geo_location=="C",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="C",]$predicted)
# Not enough observations for this test (n=2)

# Correlations per geographic location: L
cor.test(predicted.mae.bv[predicted.mae.bv$geo_location=="L",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="L",]$predicted)
#cor 
#0.9745639

# Correlations per geographic location: M
cor.test(predicted.mae.bv[predicted.mae.bv$geo_location=="M",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="M",]$predicted)
#cor 
#0.9337041

# 4.2 Plot scatterplots of predicted vs chronological age with regression lines per group
# i. Sex
tiff("plots/Correlation-predicted-chronological-age-SEX.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
plot(predicted.mae.bv$predicted~predicted.mae.bv$age,
     ylim=c(0,8),xlim=c(0,8), 
     col=rgb(0.4,0.4,0.8,0.6), pch=16, xlab='Otolith age (years)',ylab='Predicted age (years)', cex=1.4, type='n');box()
points(predicted.mae.bv[predicted.mae.bv$sex=="M",]$age, predicted.mae.bv[predicted.mae.bv$sex=="M",]$predicted, pch=21, col="black", bg=alpha("blue2",0.2), cex=1.2)
points(predicted.mae.bv[predicted.mae.bv$sex=="F",]$age, predicted.mae.bv[predicted.mae.bv$sex=="F",]$predicted, pch=21, col="black", bg=alpha("red2",0.2), cex=1.2)
abline(lm(predicted.mae.bv[predicted.mae.bv$sex=="M",]$predicted~predicted.mae.bv[predicted.mae.bv$sex=="M",]$age), col="blue2")
abline(lm(predicted.mae.bv[predicted.mae.bv$sex=="F",]$predicted~predicted.mae.bv[predicted.mae.bv$sex=="F",]$age), col="red2")
mtext(side=3,cex=1.2,"M: r = 0.958***", bty ="n", line=-1.5, adj=0.03, col="blue2") 
mtext(side=3,cex=1.2,"F: r = 0.960*** ", bty ="n", line=-3, adj=0.03, col="red2")
dev.off()

# ii. Maturity status
tiff("plots/Correlation-predicted-otolith-age-MATURITY.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
plot(predicted.mae.bv$predicted~predicted.mae.bv$age,
     ylim=c(0,8),xlim=c(0,8), 
     col=rgb(0.4,0.4,0.8,0.6), pch=16, xlab='Otolith age (years)',ylab='Predicted age (years)', cex=1.4, type='n');box()
points(predicted.mae.bv[predicted.mae.bv$maturity=="a",]$age, predicted.mae.bv[predicted.mae.bv$maturity=="a",]$predicted, pch=21, col="black", bg=alpha("blue2",0.2), cex=1.2)
points(predicted.mae.bv[predicted.mae.bv$maturity=="b",]$age, predicted.mae.bv[predicted.mae.bv$maturity=="b",]$predicted, pch=21, col="black", bg=alpha("red2",0.2), cex=1.2)
points(predicted.mae.bv[predicted.mae.bv$maturity=="c",]$age, predicted.mae.bv[predicted.mae.bv$maturity=="c",]$predicted, pch=21, col="black", bg=alpha("green3",0.2), cex=1.2)

abline(lm(predicted.mae.bv[predicted.mae.bv$maturity=="a",]$predicted~predicted.mae.bv[predicted.mae.bv$maturity=="a",]$age), col="blue2")
abline(lm(predicted.mae.bv[predicted.mae.bv$maturity=="b",]$predicted~predicted.mae.bv[predicted.mae.bv$maturity=="b",]$age), col="red2")

mtext(side=3,cex=1.2,"a: r = 0.947***", bty ="n", line=-1.5, adj=0.03, col="blue2") 
mtext(side=3,cex=1.2,"b: r = 0.859***", bty ="n", line=-3, adj=0.03, col="red2")
mtext(side=3,cex=1.2,"c: NA (n=2)", bty ="n", line=-4.5, adj=0.03, col="green3") 

dev.off()

# iii. Geographic location
tiff("plots/Correlation-predicted-chronological-age-GEOGRAPHIC-LOCATION.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
plot(predicted.mae.bv$predicted~predicted.mae.bv$age,
     ylim=c(0,8),xlim=c(0,8), 
     col=rgb(0.4,0.4,0.8,0.6), pch=16, xlab='Otolith age (years)',ylab='Predicted age (years)', cex=1.4, type='n');box()
points(predicted.mae.bv[predicted.mae.bv$geo_location=="A",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="A",]$predicted, pch=21, col="black", bg=alpha("red2",0.2), cex=1)
points(predicted.mae.bv[predicted.mae.bv$geo_location=="C",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="C",]$predicted, pch=21, col="black", bg=alpha("blue2",0.2), cex=1)
points(predicted.mae.bv[predicted.mae.bv$geo_location=="L",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="L",]$predicted, pch=21, col="black", bg=alpha("green3",0.2), cex=1)
points(predicted.mae.bv[predicted.mae.bv$geo_location=="M",]$age, predicted.mae.bv[predicted.mae.bv$geo_location=="M",]$predicted, pch=21, col="black", bg=alpha("purple",0.2), cex=1)

abline(lm(predicted.mae.bv[predicted.mae.bv$geo_location=="A",]$predicted~predicted.mae.bv[predicted.mae.bv$geo_location=="A",]$age), col="red2")
abline(lm(predicted.mae.bv[predicted.mae.bv$geo_location=="L",]$predicted~predicted.mae.bv[predicted.mae.bv$geo_location=="L",]$age), col="green3")
abline(lm(predicted.mae.bv[predicted.mae.bv$geo_location=="M",]$predicted~predicted.mae.bv[predicted.mae.bv$geo_location=="M",]$age), col="purple")

mtext(side=3,cex=1.2,"A: r = 0.945*", bty ="n", line=-1.5, adj=0.03, col="red2") 
mtext(side=3,cex=1.2,"C: NA (n=2)", bty ="n", line=-3, adj=0.03, col="blue2")
mtext(side=3,cex=1.2,"L: r=0.975***", bty ="n", line=-4.5, adj=0.03, col="green3") 
mtext(side=3,cex=1.2,"M: r=0.934***", bty ="n", line=-6, adj=0.03, col="purple") 

dev.off()
