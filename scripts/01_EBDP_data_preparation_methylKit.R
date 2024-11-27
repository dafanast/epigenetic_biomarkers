# Preparation of a single methylation matrix from several individual sample files. The process follows the standard procedure suggested by methylKit. After uniting the samples we find that very few CpGs are left and we re-iterate the process using only samples with initial >100.000 CpGs.

# Relevant information from sessionInfo()
## R version 4.3.3 (2024-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux 9.4 (Plow)
## methylKit_1.28.0     GenomicRanges_1.54.1 GenomeInfoDb_1.38.8 
## IRanges_2.36.0       S4Vectors_0.40.2     BiocGenerics_0.48.1 

# Input data: individual sample files in methylKit input format (before normalization and filtering).
# Output data: a methylBase object from the methylKit package containing methylation data from all samples.

# 1. Prepare the environment
# Load required packages
library(methylKit)

# Set working directory: CHANGE TO YOUR WORKING DIRECTORY
setwd("Joint_analysis_ML")

# 2. Prepare data
# INPUT: Read individual sample  files
# List files names
file.list <- list('102F.txt',	'106F.txt',	'107F.txt',	'111F.txt',	'118F.txt',	'123F.txt',	'124F.txt',	'12F.txt',	'136F.txt',	'140F.txt',	'144F.txt',	'148F.txt',	'154F.txt',	'155F.txt',	'158F.txt',	'15F.txt',	'169F.txt',	'175F.txt',	'181F.txt',	'183F.txt',	'186F.txt',	'190F.txt',	'195F.txt',	'196F.txt',	'198F.txt',	'205F.txt',	'206F.txt',	'207F.txt',	'214F.txt',	'223F.txt',	'228F.txt',	'22F.txt',	'238F.txt',	'242F.txt',	'247F.txt',	'248F.txt',	'249F.txt',	'250F.txt',	'37F.txt',	'43F.txt',	'48F.txt',	'53F.txt',	'57F.txt',	'58F.txt',	'60F.txt',	'65F.txt',	'6F.txt',	'71F.txt',	'72F.txt',	'73F.txt',	'74F.txt',	'75F.txt',	'76F.txt',	'7F.txt',	'81F.txt',	'82F.txt',	'92F.txt',	'96F.txt',	'97F.txt',	'9F.txt',	'275.txt',	'1.txt',	'2.txt',	'4.txt',	'5.txt',	'8.txt',	'38.txt',	'39.txt',	'40.txt',	'41.txt',	'47.txt',	'54.txt',	'56.txt',	'61.txt',	'62.txt',	'70.txt',	'77.txt',	'79.txt',	'88.txt',	'90.txt',	'94.txt',	'98.txt',	'117.txt',	'141.txt',	'26.txt',	'52.txt',	'55.txt',	'63.txt',	'80.txt',	'112.txt',	'128.txt',	'151.txt',	'153.txt',	'167.txt',	'191.txt',	'193.txt',	'212.txt',	'222.txt',	'225.txt',	'21.txt',	'23.txt',	'28.txt',	'29.txt',	'83.txt',	'84.txt',	'127.txt',	'156.txt',	'161.txt',	'163.txt',	'165.txt',	'166.txt',	'187.txt',	'188.txt',	'199.txt',	'255.txt',	'257.txt',	'258.txt',	'261.txt',	'264.txt',	'265.txt'
) 

# 3. Read files. List contains samples names and treatment vector is arbitrary since we don't have 2 groups only.
m.data <- methRead(file.list, sample.id=list('y5s1', 'y1s1', 'y7s1',	'y3s1',	'y2s1',	'y1s2',	'y1s3',	'y3s2',	'y1s4',	'y1s5',	'y2s2',	'y3s3',	'y1s6',	'y1s7',	'y2s3',	'y2s4',	'y2s5',	'y1s8',	'y1s9',	'y1s10',	'y2s6',	'y0s1',	'y0s2',	'y1s11',	'y2s7',	'y3s4',	'y3s5',	'y5s2',	'y4s1',	'y1s12',	'y1s13',	'y2s8',	'y4s2',	'y0s3',	'y4s3',	'y6s1',	'y5s3',	'y6s2',	'y1s14',	'y4s4',	'y7s2',	'y2s9',	'y2s10',	'y1s15',	'y3s6',	'y4s5',	'y0s4',	'y4s6',	'y3s7',	'y0s5',	'y0s6',	'y4s7',	'y3s8',	'y5s4',	'y4s8',	'y5s5',	'y0s7',	'y4s9',	'y3s9',	'y3s10',	'y0s8',	'y1s16',	'y1s17',	'y1s18',	'y1s19',	'y1s20',	'y1s21',	'y1s22',	'y1s23',	'y1s24',	'y1s25',	'y1s26',	'y1s27',	'y1s28',	'y1s29',	'y1s30',	'y1s31',	'y1s32',	'y1s33',	'y1s34',	'y1s35',	'y1s36',	'y1s37',	'y1s38',	'y2s11',	'y2s12',	'y2s13',	'y2s14',	'y2s15',	'y2s16',	'y2s17',	'y2s18',	'y2s19',	'y2s20',	'y2s21',	'y2s22',	'y2s23',	'y2s24',	'y2s25',	'y3s11',	'y3s12',	'y3s13',	'y3s14',	'y3s15',	'y3s16',	'y3s17',	'y3s18',	'y3s19',	'y3s20',	'y3s21',	'y3s22',	'y3s23',	'y3s24',	'y3s25',	'y4s10',	'y4s11',	'y4s12',	'y4s13',	'y4s14',	'y4s15'), 
                   assembly="cod", header=TRUE, mincov = 1, treatment=c(0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1))
# Check the object
m.data

# 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
filtered = filterByCoverage(m.data, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) 
# Check the object
filtered

# 5. Normalize coverage across samples
norm = normalizeCoverage(filtered)
# Check the object
norm

# 6. Keep CpGs present in all samples. This option has to be evaluated for each experimental design, here we start with the most conservative approach which is to only keep CpGs present in all samples. Destrand is false to be conservative.
meth = methylKit::unite(norm, destrand=FALSE)
# Check object
meth
# dim(meth)
# 405 364
# RESULT: Only 405 CpGs are left after the procedure. This is likely to occur if few samples have very low number of reads and are driving down all samples CpGs after the unite function. Our strategy in this case will be to identify these samples, eliminate them and start the process again.

# START AGAIN. 2. Prepare data
# INPUT: Read individual sample  files
# List files names
file.list <- list('102F.txt',	'106F.txt',	'107F.txt',	'111F.txt',	'118F.txt',	'123F.txt',	'124F.txt',	'12F.txt',	'136F.txt',	'140F.txt',	'144F.txt',	'148F.txt',	'154F.txt',	'155F.txt',	'158F.txt',	'15F.txt',	'169F.txt',	'175F.txt',	'181F.txt',	'183F.txt',	'190F.txt',	'195F.txt',	'196F.txt',	'214F.txt',	'223F.txt',	'228F.txt',	'22F.txt',	'238F.txt',	'242F.txt',	'247F.txt',	'248F.txt',	'249F.txt',	'250F.txt',	'37F.txt',	'43F.txt',	'48F.txt',	'53F.txt',	'57F.txt',	'58F.txt',	'60F.txt',	'65F.txt',	'6F.txt',	'71F.txt',	'72F.txt',	'74F.txt',	'75F.txt',	'76F.txt',	'7F.txt',	'81F.txt',	'82F.txt',	'92F.txt',	'96F.txt',	'97F.txt',	'9F.txt',	'275.txt',	'1.txt',	'2.txt',	'4.txt',	'5.txt',	'8.txt',	'38.txt',	'39.txt',	'40.txt',	'41.txt',	'47.txt',	'54.txt',	'56.txt',	'61.txt',	'62.txt',	'70.txt',	'79.txt',	'88.txt',	'90.txt',	'94.txt',	'98.txt',	'117.txt',	'141.txt',	'26.txt',	'52.txt',	'55.txt',	'63.txt',	'80.txt',	'112.txt',	'128.txt',	'151.txt',	'153.txt',	'167.txt',	'191.txt',	'193.txt',	'212.txt',	'222.txt',	'225.txt',	'21.txt',	'23.txt',	'28.txt',	'29.txt',	'83.txt',	'84.txt',	'127.txt',	'161.txt',	'163.txt',	'165.txt',	'166.txt',	'188.txt',	'199.txt',	'255.txt',	'257.txt',	'258.txt',	'261.txt',	'264.txt') 

# 3. Read files. List contains samples names and treatment vector is arbitrary since we don't have 2 groups only.
m.data <- methRead(file.list, sample.id=list('y5s1', 'y1s1',	'y7s1',	'y3s1',	'y2s1',	'y1s2',	'y1s3',	'y3s2',	'y1s4',	'y1s5',	'y2s2',	'y3s3',	'y1s6',	'y1s7',	'y2s3',	'y2s4',	'y2s5',	'y1s8',	'y1s9',	'y1s10',	'y0s1',	'y0s2',	'y1s11',	'y4s1',	'y1s12',	'y1s13',	'y2s8',	'y4s2',	'y0s3',	'y4s3',	'y6s1',	'y5s3',	'y6s2',	'y1s14',	'y4s4',	'y7s2',	'y2s9',	'y2s10',	'y1s15',	'y3s6',	'y4s5',	'y0s4',	'y4s6',	'y3s7',	'y0s6',	'y4s7',	'y3s8',	'y5s4',	'y4s8',	'y5s5',	'y0s7',	'y4s9',	'y3s9',	'y3s10',	'y0s8',	'y1s16',	'y1s17',	'y1s18',	'y1s19',	'y1s20',	'y1s21',	'y1s22',	'y1s23',	'y1s24',	'y1s25',	'y1s26',	'y1s27',	'y1s28',	'y1s29',	'y1s30',	'y1s32',	'y1s33',	'y1s34',	'y1s35',	'y1s36',	'y1s37',	'y1s38',	'y2s11',	'y2s12',	'y2s13',	'y2s14',	'y2s15',	'y2s16',	'y2s17',	'y2s18',	'y2s19',	'y2s20',	'y2s21',	'y2s22',	'y2s23',	'y2s24',	'y2s25',	'y3s11',	'y3s12',	'y3s13',	'y3s14',	'y3s15',	'y3s16',	'y3s17',	'y3s19',	'y3s20',	'y3s21',	'y3s22',	'y3s24',	'y3s25',	'y4s10',	'y4s11',	'y4s12',	'y4s13',	'y4s14'),  
                                             assembly="cod", mincov = 1, treatment=c(0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1, 1))
# Check the object
m.data

# 4. Filter for reads with very low number of reads (<10) and exceeding the 99.9% percentile
filtered = filterByCoverage(m.data, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9) 
# Check the object
filtered

# 5. Normalize coverage across samples
norm = normalizeCoverage(filtered)
# Check the object
norm

# 6. Keep CpGs present in 48 samples per group. We relaxed this criterion to obtain more CpGs but later we will need to deal with missing data.
meth = methylKit::unite(norm, destrand=FALSE, min.per.group=48L)
# Check object
meth
# dim(meth)
# 85735   334 

# OPTIONAL: 7. Save the object
# save(meth, file="data/meth-10cov-100000CpGs.Rdata")
