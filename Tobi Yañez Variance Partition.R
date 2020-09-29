# Installing BiocManager + libraries required for analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("variancePartition")
BiocManager::install("GEOquery")
BiocManager::install("oligo") #used to normalize the data downloaded

install.packages("vctrs") #current version needed for geoquery to function
library(variancePartition) 
library(ggplot2)
library(Biobase)
library(GEOquery)
library(oligo)
library(stringr) #for cleaning the data

# help libraries for variance partition and GEOquery
browseVignettes("variancePartition")  # used heavily for analysis
browseVignettes("GEOquery") #used for importing data
browseVignettes("oligo") #used for preparing data

#importing the data
gset <- getGEO("GSE56035", GSEMatrix =TRUE) 
#list of 1 containing df with all the data

gset1 <- gset[[1]] #extracting the element to make it into a df

View(gset1)        #Viewing to make sure it works - it does!


#cleaning up the data via code found on Kasper Daniel Hansen's Github
filename <- sampleNames(gset1)
pData(gset1)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(gset1) <- sampleNames
View(gset1)
pData(gset1)$group <- ifelse(grepl("^OSA", sampleNames(gset1)),
                             "OSA", "Control")
gset2 <- pData(gset1)
View(gset2)
View(gset1@assayData[["exprs"]]) #Assay data necessary for final formula + method

# assay data nxp size
GSEdata <- gset1@assayData[["exprs"]]
View(GSEdata)

#creating info df from gset2 to go with assay data pxn size
Infodf <- data.frame(matrix(ncol = 8, nrow = 984))
names(Infodf) <- c("Individual","Age","Sex","CellType","Batch","InclusionMarkers","ExclusionMarkers","PhenotypeMarkers")
rownames(Infodf)<- rownames(gset2)

#filling the columns with variables I want to check influence of
Infodf$Age <- as.numeric(str_extract(gset2$characteristics_ch1,"\\w+$"))
Infodf$Sex <- as.factor(tolower(str_extract(gset2$characteristics_ch1.1,"\\w+$")))
Infodf$CellType <- as.factor(gset2$characteristics_ch1.2)
Infodf$Batch <- as.factor(str_extract(gset2$characteristics_ch1.3,"\\w+$"))
Infodf$InclusionMarkers <- as.factor(gset2$characteristics_ch1.4) # ultimately unecessary
Infodf$ExclusionMarkers <- as.factor(gset2$characteristics_ch1.5) # ultimately unecessary
Infodf$PhenotypeMarkers <- as.factor(gset2$characteristics_ch1.6) # ultimately unecessary
Infodf$Individual <- rownames(Infodf)
length(unique(Infodf$Individual))
#Do not have data on which individual so it is impossible to add it into final model 
#since I have 984 unique values for individual there are no repeats.
#This means there is no solution to the matrix which solves for variance in the 
#model.

#formula including the variables I want to include in the analysis
#had all of the above except individual but there was clear multicolinearity
#between inclusion markers, exclusion markers, phenotype markers, and cell type
#so I decided to only keep cell type
formula <-  ~ Age + (1|Sex) + (1|CellType) + (1|Batch)

#Creating final model
varPart <- fitExtractVarPartModel(GSEdata, formula, Infodf)

# finding correlation statistics for all variables
# Multicollinearity between 4 afformentioned variables
form <- ~ Batch + Age + Sex + CellType + ExclusionMarkers + PhenotypeMarkers + InclusionMarkers
C = canCorPairs(form, Infodf)
plotCorrMatrix( C )

#Creating images that show final results
vp <- sortCols(varPart)
plotPercentBars(vp[1:10,])
plotVarPart(vp)
# Batch is very small for unknown reason 

# Numeric representations
summary(varPart)
vpSummaries <- fitVarPartModel( GSEdata, formula, Infodf, fxn=summary )
# Did not run second line because I do not think it would add much
# to the analysis and it uses up a ton of data. But in an actual
# study it should be run


# Plots depicting the largest and smallest differences in gene expression
# depending on sex. A lot of the code, not only here, was taken directly
# from the vignette
i <- which.max(varPart$Sex)
GE <- data.frame( Expression = GSEdata[i,], Sex = Infodf$Sex)
plotStratify(Expression ~ Sex, GE, main=rownames(GSEdata)[i]) 

i <- which.min(varPart$Sex)
GE <- data.frame( Expression = GSEdata[i,], Sex = Infodf$Sex)
plotStratify(Expression ~ Sex, GE, main=rownames(GSEdata)[i]) 

