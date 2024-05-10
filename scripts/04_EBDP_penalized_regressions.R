# Input data: a matrix with correctly named unique CpGs as columns and samples as rows, without missing values. If any technical bias, e.g. batch effects, have been detected they should have been corrected. 
# Output data: final epigenetic clock for age prediction and plots of performance.

# 1. Prepare the environment
# Load required packages
library(tidyr)
library(dplyr)
library(tidymodels)
library(readr)
library(caret)

# Set working directory
#setwd("/powerplant/workspace/cfndxa/Cod_clock/Joint_analysis_ML/github")

# 2. Prepare data
# INPUT: Load data from previous session
load("data/meth-corrected-batch-effect.Rdata")
# Read samples information
samples.age <- read.table("data/samples-age.txt", sep="\t", stringsAsFactors = FALSE, quote="", header=TRUE) 
# Add the age column to the dataframe
meth.complete <- as.data.frame(meth.corrected.batch.t) %>% mutate(age=samples.age$age)

#head(meth.complete[1:10,1:10])
#X1.15075 X1.15144 X1.53225 X1.53285 X1.53287 X1.53308 X1.69477 X1.69554  X1.69588 X1.87544
#1  91.66667 87.87879 72.00000 55.55556 50.00000 66.66667 76.59574 51.29534  82.90155 73.40426
#2  96.00000 78.85714 59.73154 75.81699 67.97386 53.59477 92.02899 65.21739  94.77612 91.27907
#3  93.63057 90.78947 47.70642 74.31193 69.72477 30.27523 83.03571 68.30357  91.12150 88.88889
#4  87.50000 82.60870 61.65803 52.33161 47.66839 47.66839 83.13253 58.43373 100.00000 95.10870
#5 100.00000 90.90909 34.72222 45.83333 57.63889 50.00000 86.17021 56.73759  90.07092 86.25000
#6  89.07563 94.95798 41.12150 64.48598 87.85047 58.87850 90.97222 62.25166  96.02649 94.95798

# Check that the dataframe has the correct dimensions
#dim(meth.complete)
#[1]   110 47464

# 3. Data splitting
set.seed(123)
# Split stratified by age. By default the command splits in 75% training and 25% test. This should be fine to start with but may need to be adjusted in cases, e.g. of overfitting (where a model performs extremely well in the training but not in the dataset), to 80% training/20% test.
splits      <- initial_split(meth.complete, strata = age, prop=4/5) 

# Save training set
age_other <- training(splits)
# Save test set
age_test  <- testing(splits)

# Check training set proportions by age class
age_other %>% count(age) %>% mutate(prop = n/sum(n))
dim(age_other)
# Check test set proportions by age class
age_test %>% count(age) %>% mutate(prop = n/sum(n))
dim(age_test)

# 4. Data preparation for Machine Learning models. Only training dataset.
## 4.1. Exclude variables with zero or near zero variance. For more options check caret documentation.
# Identify variables
nzv.cpg.list <- nearZeroVar(age_other, freqCut = 85/15, uniqueCut = 50, allowParallel = TRUE)
# Exclude them from the training dataset
filteredDescr <- age_other[, -nzv.cpg.list]
# Count remaining variables
dim(filteredDescr)
#[1]    86 46359
# Optional: save the object for use in new R sessions
#save(filteredDescr, file="data/Filtered-ZeroVar.RData")

## 4.2. Exclude variables that are highly correlated between them.
# Exclude age column
filteredDescr.woage <- filteredDescr[,-ncol(filteredDescr)]
# Calculate pairwise correlations between DNA methylation variables
descrCor <- cor(filteredDescr.woage) # This can be a huge object, several GB, so may be worth deleting it after finishing with the correlations. 

# Way 1: Define the correlation threshold (e.g., 0.7)
threshold <- 0.7
# Find indices of highly correlated variables
highly_correlated <- which(abs(descrCor) >= threshold & upper.tri(descrCor), arr.ind = TRUE)
# Get the names of the variables
variable_names <- colnames(descrCor)
# Extract the names of highly correlated variables
highly_correlated_names <- c(
  variable_names[highly_correlated[, 1]],
  variable_names[highly_correlated[, 2]]
)
# Remove duplicates (since the correlation matrix is symmetric)
highly_correlated_names <- unique(highly_correlated_names)
# Exclude highly correlated variables
filteredDescr.cor <- filteredDescr[, !colnames(filteredDescr) %in% highly_correlated_names]
dim(filteredDescr.cor)
#[1]     86 35328

# Way 2: Find correlations above a specific threshold using caret package. It throws an error in my case due to huge dataset.
#highlyCorDescr <- findCorrelation(descrCor, cutoff = .7, exact=FALSE)
# Filter out variables with a high correlation from the training dataset
#filteredDescr.cor <- filteredDescr[,-highlyCorDescr]
#dim(filteredDescr.cor)

# Optional: save the object for use in new R sessions
#save(filteredDescr.cor, file="data/HighlyCorVariables-Training.RData")

## 4.3. OPTIONAL: this step improves the performance of the models in the test sets. We exclude variables (CpGs) with very high or very low DNA methylation. We will set the threshold of DNA methylation level (e.g. 90% for high and 10% for low) and the percentage of individuals among the dataset that exhibit DNA methylation levels above this threshold.
# Exclude variables above a certain DNA methylation level
exclude_columns_above_threshold <- function(data, threshold, percentage) {
  num_samples <- nrow(data)
  # Filter out non-numeric columns
  numeric_columns <- sapply(data, is.numeric)
  above_threshold <- sapply(data[, numeric_columns], function(col) sum(col > threshold, na.rm = TRUE))
  columns_to_exclude <- names(above_threshold[above_threshold > num_samples * percentage])
  filtered_data <- data[, !(names(data) %in% columns_to_exclude)]
    return(filtered_data)
}
# Set threshold and percentage
threshold <- 80
percentage <- 0.7
# Filter columns
filtered.high.meth <- exclude_columns_above_threshold(filteredDescr.cor, threshold, percentage)
# Check how many variables were filtered
dim(filtered.high.meth)
#[1]  86 12874

# Exclude variables below a certain DNA methylation level
exclude_columns_below_threshold <- function(data, threshold, percentage) {
  num_samples <- nrow(data)
  # Filter out non-numeric columns
  numeric_columns <- sapply(data, is.numeric)
  below_threshold <- sapply(data[, numeric_columns], function(col) sum(col < threshold, na.rm = TRUE))
  columns_to_exclude <- names(below_threshold[below_threshold > num_samples * percentage])
  filtered_data <- data[, !(names(data) %in% columns_to_exclude)]
  return(filtered_data)
}
# Set threshold and percentage
threshold <- 20
percentage <- 0.7
# Filter columns
filtered.low.meth <- exclude_columns_below_threshold(filtered.high.meth, threshold, percentage)
# Check how many variables were filtered.
dim(filtered.low.meth)
#[1]  86 12610
# There is a possibility that the age column was excluded during this step so we need to add it back
filtered.low.meth$age <- filtered.high.meth$age

# Optional: save the object for use in new R sessions
#save(filtered.low.meth, file="data/ExtremeValues-Training.RData")

## 4.4. Epigenome Wide Association Studies (EWAS). OPTIONAL: this step filters for CpGs correlated with age. It may be used to reduce the number of variables to train the models with, however, we haven't found any improvement in the model performance after this step. It requires the package WGCNA. 
library(WGCNA)
# Keep age in a vector
x.age <- as.vector(filtered.low.meth$age)
# Convert to matrix
y.cpgs <- apply(as.matrix.noquote(filtered.low.meth),  # Using apply function
                2,
                as.numeric)
# Amend the CpG names
cpg.names <- colnames(filtered.low.meth)
# Calculate correlation between age and CpGs. 
cor.age <- cor(x=x.age, y=y.cpgs, use="pairwise.complete.obs", method="pearson")
# Calculate correlation between age and CpGs and include testing.
cor.age.pval <- corAndPvalue(x=x.age, y=y.cpgs, use="pairwise.complete.obs", alternative="two.sided")
# Construct a dataframe with the metrics.
cor.pval.names <- bind_cols(cpgs=cpg.names, coeff=cor.age[1:ncol(cor.age)], cor=cor.age.pval$bicor[1:ncol(cor.age)], p=cor.age.pval$p[1:ncol(cor.age)], Z=cor.age.pval$Z[1:ncol(cor.age)], t=cor.age.pval$t[1:ncol(cor.age)])
# Subset according to thresholds of the correlation testing. Values can be changed depending on the dataset, experiment or model. 
cors.filt <- subset(cor.pval.names, abs(coeff)>0.5&p<0.01)
# Plots an histogram useful to see the distribution of the coefficients.
hist(cors.filt$coeff, xlab="Correlation coefficient", main="")
# Check how many CpGs will be filtered
dim(cors.filt)
# Vector with names of CpGs to be filtered
vars.del <- cors.filt$cpgs
# Filter our dataset for CpGs highly correlated with  Select Vector with names of CpGs to be filtered
all.cors.filt.df <- select(filtered.low.meth, all_of(vars.del)) %>% mutate(age=filtered.low.meth$age)
dim(all.cors.filt.df)
# Optional: save the object for use in new R sessions
#save(all.cors.filt.df, file="CorVariablesAge-Training.RData")

## 4.5. PreProcess dataset by centering and scaling variables. We skipped step 4.4 in the Atlantic cod so we will continue with the output of step 4.3.
# Change the right part depending on which step you stopped
filtered.dataset <- filtered.low.meth
# Use caret's function preProcess that outputs a "preProcess" object.
preProcValues <- preProcess(filtered.dataset[,-ncol(filtered.dataset)], method = c("center", "scale"))
# Use predict function to apply the transformations to the training dataset.
trainTransformed <- predict(preProcValues, filtered.dataset)
# Check how many CpGs we have.
dim(trainTransformed)
#[1]   86 12611 # one more column because we added age
# Optional: save the object for use in new R sessions
#save(trainTransformed, file="data/Features-ForModelTuning-Training.RData")

# 5. Machine Learning models building. Evaluation using resampling techniques. See caret documentation for details.
## 5.1 Set resampling technique 
fitControl <- trainControl(method = 'cv', number=10) 

## 5.2 Choose grid of lambda to be tested
lambda <- 10^seq(-3, 3, length = 100)

# 5.3 Apply age transformation: it may be necessary to apply a transformation of age to improve prediction as in mammalian species. If an age transformation is applied, then the inverse function needs to be applied downstream. Which age transformation is appropriate needs to be evaluated using the results of prediction on test data. For example, in Atlantic cod the log-transformation greatly improved the age prediction in test data. See https://doi.org/10.1038/s43587-023-00462-6 and https://github.com/jazoller96/mammalian-methyl-clocks/tree/main for details on possible age transformations used in mammalian epigenetic clocks. 
# Examples:
# i. Standard log transformation where we add 0.1 to avoid 0s:
age.train.transformed <- log(trainTransformed$age+0.1) 
# ii. Transformation used in the second universal mammalian clock where RelativeAge=Age/MaximumAge. Maximum age in Atlantic cod is 25 years according to https://genomics.senescence.info/species/entry.php?species=Gadus_morhua. We add 30% to the reported maximum age to account for uncertainty for wild species as in mammals (https://doi.org/10.1038/s43587-023-00462-6):
#age.train.transformed <- -log(-log((trainTransformed$age+0.1)/32.5)) 
# iii. Square root log transformation where we add 0.1 to avoid 0s:
#age.train.transformed <- sqrt(trainTransformed$age+0.1) ) 

## 5.4 Run penalized regressions.
# Using bsRAD-seq or similar kind of data, we have >10.000 CpGs so we won't consider ridge regression because it keeps all CpGs and an epigenetic clock with >10.000 is not practical. When using targeted bisulfite sequencing or a similar technique, ridge may be useful to test.  
# Due to a very high number of variables (>10.000), models may have problems to run. To circumvent this situation the arguments x and y of the command 'train' should be written as below. This solution comes from https://stackoverflow.com/questions/27677391/r-stack-overflow-error-with-randomforest-on-large-dataset-48-512-gb-ram.

# a. LASSO
set.seed(123)
lasso_model <- train(x=trainTransformed[,colnames(trainTransformed) !="age"], y=age.train.transformed, method = "glmnet", trControl = fitControl, tuneGrid = data.frame(alpha = 1, lambda = 10^seq(-3, 3, length = 100)), tuneLength = 10)
#save(lasso_model, file="data/lasso_model.Rdata")
# b. Elastic net with best alpha automatically chosen
set.seed(123)
elastic_model <- train(x=trainTransformed[,colnames(trainTransformed) !="age"], y=age.train.transformed, method = "glmnet", trControl = fitControl, tuneLength = 10)
#save(elastic_model, file="data/elastic_model.Rdata")
# c. Optional: Elastic net with alpha set to 0.5. Some authors choose to manually set alpha to 0.5 instead of allowing the best alpha to be automatically selected. Alpha is linked to how many variables will be selected in the final model.  
#set.seed(123)
#elastic_model.05 <- train(x=trainTransformed[,colnames(trainTransformed) !="age"], y=trainTransformed$age, method = "glmnet", trControl = fitControl, tuneGrid = data.frame(alpha = 0.5, lambda = 10^seq(-3, 3, length = 100)), tuneLength = 10)
# d. Optional: Ridge regression. In case of targeted bisulfite sequencing, ridge regression may be useful to test. 
#set.seed(123)
#ridge_model <- train(x=trainTransformed[,colnames(trainTransformed) !="age"], y=trainTransformed$age, method = "glmnet", trControl = fitControl, tuneGrid = data.frame(alpha = 0, lambda = 10^seq(-3, 3, length = 100)), tuneLength = 10)

# 6. Compare metrics of the models. 
## 6. 1 We only use lasso and elastic net for downstream analysis (see above for explanations). The goal is to lower RMSE/MAE and increase Rsquared, but these metrics need to be taken in combination. The number of final CpGs may also important for practical application of the epigenetic clock.
models_compare <- resamples(list(LM=lasso_model, EM=elastic_model))
summary(models_compare)
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)

## 6.2 Count CpGs (features) kept by each model
# a. LASSO
#coef(lasso_model$finalModel, lasso_model$bestTune$lambda)
sum(coef(lasso_model$finalModel, lasso_model$bestTune$lambda)!=0)
# b. Elastic net
#coef(elastic_model$finalModel, elastic_model$bestTune$lambda)
sum(coef(elastic_model$finalModel, elastic_model$bestTune$lambda)!=0)

##  6.3 Compare metrics in the training datasets.
# a. LASSO
# Use model to predict age
predicted.age <- predict.train(lasso_model)
# Calculate resampling metrics
postResample(pred = predicted.age, age.train.transformed)
# Calculate correlation of predicted and actual age
cor.test(predicted.age, age.train.transformed)

# b. Elastic net
# Use model to predict age
predicted.age <- predict.train(elastic_model)
# Calculate resampling metrics
postResample(pred = predicted.age, age.train.transformed)
# Calculate correlation of predicted and actual age
cor.test(predicted.age, age.train.transformed)

##  6.4 Compare metrics in the test datasets. This is a crucial step which will allow us to determine how well the model generalizes, how well it performs in completely unseen data. It is possible that the model performs extremely well in training (>99% Rsquared) and severely underperforms in test (e.g. 60% Rsquared). This is a case of "overfitting" and needs to be handled accordingly, so one has to choose one or several strategies to deal with it (e.g., increase training dataset, decrease number of variables, resampling etc). If results of the tests of this section are not good, go back to steps 3, 4, 5 either all of them or some of them. The final model is expected to perform well in the test set.

# Transformations:
# i. Transformation of methylation data. We use the transformations (PreProcessing step 4.5) calculated based on the training data and apply it to the test data.
testTransformed <- predict(preProcValues, age_test)

# ii. Transformation of age: the same transformation used for the training dataset needs to be applied in the test.
age.test.transformed <- log(testTransformed$age+0.1) 

# a. LASSO
# Use model to predict age
predict.lasso.test <- predict(lasso_model, testTransformed)
# Calculate resampling metrics
postResample(pred = predict.lasso.test, age.test.transformed)
#RMSE  Rsquared       MAE 
#0.4376038 0.8350900 0.3336083 
# Calculate correlation of predicted and actual age
cor.test(predict.lasso.test, age.test.transformed)

# b. Elastic net
# Use model to predict age
predict.enet.test <- predict(elastic_model, testTransformed)
# Calculate resampling metrics
postResample(pred = predict.enet.test, age.test.transformed)
#RMSE  Rsquared       MAE 
#0.4599044 0.8293261 0.3236262
# Calculate correlation of predicted and actual age
cor.test(predict.enet.test, age.test.transformed)

## 7. Final model. Choose lasso as the final model because it keeps the lowest amount of CpGs with the best performance measures.
## 7.1 No resampling technique 
fitControl <- trainControl(method = "none") 

## 7.2 Find lambda value
#lasso_model$finalModel
lasso_model$bestTune
#lambda = 0.02477076

## 7.3 Run final model
set.seed(123)
lasso_model_all <- train(x=trainTransformed[,colnames(trainTransformed) !="age"], y=age.train.transformed, method = "glmnet", trControl = fitControl, tuneGrid = data.frame(alpha = 1, lambda = 0.02477076))
#save(lasso_model_all, file="data/FinalLASSOModel.Rdata")

## 7.4 Use the final model to predict age in the original dataframe
# Transform methylation data
final.all.trans <- predict(preProcValues, meth.complete)
#save(final.all.trans, file="MethDataForModel.Rdata")
# Use model to predict age
predicted.age.final <- predict(lasso_model_all,  final.all.trans)

# Transform age
meth.age.transformed <- log(meth.complete$age+0.1)
# Calculate resampling metrics
postResample(pred = predicted.age.final, meth.age.transformed)
# Calculate correlation of predicted and actual age
cor.test(predicted.age.final, meth.age.transformed)

## 7.5 Plot the final epigenetic clock
# Inverse age transformation for plotting
predicted.age.inversed <- round(exp(predicted.age.final)-0.1, 2)
# Calculate Mean Absolute Error (MAE) based on inversed ages
mae <- abs(predicted.age.inversed-meth.complete$age)
# save(mae, file="data/FinalModelMAE.Rdata")
# save(predicted.age.inversed, file="data/FinalModelPredictedAge.Rdata")

# Plot final clock
tiff("plots/Correlation-predicted-otolith-age.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
plot(predicted.age.inversed~meth.complete$age,
     ylim=c(0,8),xlim=c(0,8), 
     col=rgb(0.4,0.4,0.8,0.6), pch=16, xlab='Otolith age (years)',ylab='Predicted age (years)', cex=1.4);box()
abline(lm(predicted.age.inversed~meth.complete$age))
mtext(side=3,cex=1.2,"r = 0.975", bty ="n", line=-1.5, adj=0.03) # Value of correlation of transformed ages
mtext(side=3,cex=1.2,expression(paste(italic(p)," = 2.2"^"e-16")), bty ="n", line=-3, adj=0.03)
mtext(side=3,cex=1.2,paste0("MAE=",round(mean(mae),3)), bty ="n", line=-4.5, adj=0.03) 
dev.off()

## 7.6 Plot separately training and test. Plot training and test error.
# a. Training set
# Use final model to predict age in training set
predicted.age.train <- predict.train(lasso_model_all)
# Calculate resampling metrics
postResample(pred = predicted.age.train, age.train.transformed)
# Calculate correlation of predicted and actual age
cor.test(predicted.age.train, age.train.transformed)
# Inverse log-transformed age for plotting
predicted.age.train.inversed <- round(exp(predicted.age.train)-0.1, 2)
# Calculate MAE between predicted and actual age inversed
mae.train <- mean(abs(predicted.age.train.inversed-trainTransformed$age))

# b. Test set
# testTransformed <- predict(preProcValues, age_test)
# Use final model to predict age in test set
predict.lasso.test <- predict(lasso_model_all, testTransformed)
# Calculate resampling metrics
postResample(pred = predict.lasso.test, age.test.transformed)
# Calculate correlation of predicted and actual age
cor.test(predict.lasso.test, age.test.transformed)
# Inverse log-transformed age for plotting
predicted.age.test.inversed <- round(exp(predict.lasso.test)-0.1, 2)
# Calculate MAE between predicted and actual age inversed
mae.test <- mean(abs(predicted.age.test.inversed-testTransformed$age))

# Plot clock on training set only
tiff("plots/Correlation-predicted-otolith-age-Train.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
plot(predicted.age.train.inversed~trainTransformed$age, ylim=c(0,8),xlim=c(0,8), col=rgb(0.4,0.4,0.8,0.6), pch=16, cex=1.4, xlab='Otolith age (years)',ylab='Predicted age (years)');box()
abline(lm(predicted.age.train.inversed~trainTransformed$age))
mtext(side=3,cex=1.2,"r = 0.998", bty ="n", line=-1.5, adj=0.03) # Value of correlation of transformed ages
mtext(side=3,cex=1.2,expression(paste(italic(p)," = 2.2"^"e-16")), bty ="n", line=-3, adj=0.03)
mtext(side=3,cex=1.2,paste0("MAE=",round(mean(mae.train),3)), bty ="n", line=-4.5, adj=0.03) # Value of error of transformed ages
dev.off()

# Plot clock on test set only
tiff("plots/Correlation-predicted-otolith-age-Test.tiff",
     width=14.5, height=14.5, units="cm", res=350)
par(mar=c(4.5, 3.7, 1.1, 1.1), mgp = c(2.5, 0.8, 0), oma=c(0,0,0,0), cex.axis=1.4, cex.lab=1.6)
plot(predicted.age.test.inversed~testTransformed$age, ylim=c(0,8),xlim=c(0,8), col=rgb(0.4,0.4,0.8,0.6), pch=16, cex=1.4, xlab='Otolith age (years)',ylab='Predicted age (years)');box()
abline(lm(predicted.age.test.inversed~testTransformed$age))
mtext(side=3,cex=1.2,"r = 0.914", bty ="n", line=-1.5, adj=0.03)  # Value of correlation of transformed ages
mtext(side=3,cex=1.2,expression(paste(italic(p)," = 2.2"^"e-16")), bty ="n", line=-3, adj=0.03)
mtext(side=3,cex=1.2,paste0("MAE=",round(mean(mae.test),3)), bty ="n", line=-4.5, adj=0.03) # Value of error of transformed ages
dev.off()