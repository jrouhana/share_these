rm(list=ls())

#Uncomment below if you need these installed
#source("https://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("affy")
#biocLite("hgu133a.db")
#biocLite("limma)

library(limma)
library(affy)
library(hgu133a.db)
library(GEOquery)

#Uncomment below to get files needed to run
#Set to a clean working directory
#A lot of .cel files will be downloaded
#setwd("/home/jrouhana/workspace/testR")
#getGEOSuppFiles("GSE2990")
#untar("GSE2990/GSE2990_RAW.tar", exdir = "./GSE2990/")
#setwd("./GSE2990/")

#Take a look at the data
targets <- readTargets("GSE2990_suppl_info.txt", sep="")

#Change the file names to resemble actual file structure
targets$filename <- gsub("GSM", "gsm", 
                         gsub(" ", "", paste(targets$geo_accn, ".cel.gz"), 
                              fixed=T), fixed=T)

#Normalize separately based on whether patients were treated
tamoPatients <- which(targets$treatment=="tamoxifen")
noTamoPatients <- which(targets$treatment!="tamoxifen")

tamo <- targets[tamoPatients, ]
noTamo <- targets[noTamoPatients,]

#Get expression data for patients, then normalize
tamoBatch <- ReadAffy(filenames=tamo$filename)
noTamoBatch <- ReadAffy(filenames=noTamo$filename)
esetTamo <- rma(tamoBatch)
esetNoTamo <- rma(noTamoBatch)


#Get histological grades for patients treated with tamoxifen
treatedPatients <- which(targets$treatment=="tamoxifen")
f <- paste(targets[treatedPatients,"grade"], sep="")
f <- factor(f)

#Design contrasts, fit model
design <- model.matrix(~0+f)
colnames(design) <- c("one", "three", "not.given")
design

fit <- lmFit(esetTamo, design)

cont.matrix <- makeContrasts("three-one",levels=design)
cont.matrix

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

#Get G1/G3 enriched genes. Adjust for fdr
tempTable <- topTable(fit2, coef=1, adjust="fdr", 
                 sort.by = "P", num=length(fit2$p.value))

g3Significant <- which(tempTable$adj.P.Val<0.01&tempTable$t>0)
g1Significant <- which(tempTable$adj.P.Val<0.01&tempTable$t<0)

g3Enriched <- rownames(tempTable[g3Significant,])
g1Enriched <- rownames(tempTable[g1Significant,])

#Map the names to get HGNC symbols & EntrezIDs
g3Enriched <- select(hgu133a.db, 
                     g3Enriched, c("SYMBOL","ENTREZID", "GENENAME"))

g1Enriched <- select(hgu133a.db, 
                     g1Enriched, c("SYMBOL","ENTREZID", "GENENAME"))

#Look only at expression leves for enriched genes
tamoValues <- data.frame(exprs(esetTamo))
noTamoValues <- data.frame(exprs(esetNoTamo))

tamoValues <- tamoValues[which(
  rownames(tamoValues) %in% c(g3Enriched$PROBEID, g1Enriched$PROBEID)),]

noTamoValues <- noTamoValues[which(
  rownames(tamoValues) %in% c(g3Enriched$PROBEID, g1Enriched$PROBEID)),]


#Get a difference of sums for enriched G3 genes - G1 genes
tamoGradeIndex <- NA
for (colIdx in 1:length(tamoValues[1,])){
  g3Idx <- which(rownames(tamoValues) %in% g3Enriched$PROBEID)
  g1Idx <- which(rownames(tamoValues) %in% g1Enriched$PROBEID)
  
  g3Sum <- sum(tamoValues[g3Idx,colIdx], na.rm=T)
  g1Sum <- sum(tamoValues[g1Idx,colIdx], na.rm=T)
  tamoGradeIndex[colIdx] <- g3Sum-g1Sum
}

tamoGradeIndex <- scale(tamoGradeIndex)
tamoValues["Grade Index",] <- tamoGradeIndex

#Check how well grade index correctly classified G3 vs G1
upperThreshold <- qt(0.65, length(tamoValues))

targets[,"Pred.Grade"] <- 2

predictedG3 <- which(targets$filename %in% 
                       colnames(tamoValues[,which(tamoValues["Grade Index",]>upperThreshold)]) )
targets[predictedG3 ,"Pred.Grade"] <- 3

lowerThreshold <- qt(0.35, length(tamoValues))

predictedG1 <- which(targets$filename %in% 
                       colnames(tamoValues[,which(tamoValues["Grade Index",] < lowerThreshold)]) )
targets[predictedG1,"Pred.Grade"] <- 1

#A lot still land in G2, the "in-question" grade
#Those that do  get classified are classified correctly (mostly)
table(targets[tamoPatients,"grade"], targets[tamoPatients ,"Pred.Grade"])

#This is data we did not train with
noTamoGradeIndex <- NA
for (colIdx in 1:length(noTamoValues[1,])){
  g3Idx <- which(rownames(noTamoValues) %in% g3Enriched$PROBEID)
  g1Idx <- which(rownames(noTamoValues) %in% g1Enriched$PROBEID)
  
  g3Sum <- sum(noTamoValues[g3Idx,colIdx], na.rm=T)
  g1Sum <- sum(noTamoValues[g1Idx,colIdx], na.rm=T)
  noTamoGradeIndex[colIdx] <- g3Sum-g1Sum
}

noTamoGradeIndex <- scale(noTamoGradeIndex)
noTamoValues["Grade Index",] <- noTamoGradeIndex

upperThreshold <- qt(0.75, length(noTamoValues))

predictedG3 <- which(targets$filename %in% 
                       colnames(noTamoValues[,which(noTamoValues["Grade Index",]>upperThreshold)]) )
targets[predictedG3 ,"Pred.Grade"] <- 3

lowerThreshold <- qt(0.25, length(noTamoValues))

predictedG1 <- which(targets$filename %in% 
                       colnames(noTamoValues[,which(noTamoValues["Grade Index",] < lowerThreshold)]) )

#For the most part, no G3 land in G1, and no G1 land in G3
#G2 is more spread out than pathologists' classification, 
#which is possibly a good thing
targets[predictedG1,"Pred.Grade"] <- 1
table(targets[noTamoPatients,"grade"], targets[noTamoPatients ,"Pred.Grade"])



## 
## Please cite the following if utilizing the GEOquery software:
## 
##   Davis, S. and Meltzer, P. S. GEOquery: a bridge between the Gene
##   Expression Omnibus (GEO) and BioConductor. Bioinformatics, 2007,
##   14, 1846-1847
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     author = {Sean Davis and Paul Meltzer},
##     title = {GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor},
##     journal = {Bioinformatics},
##     year = {2007},
##     volume = {14},
##     pages = {1846--1847},
##   }