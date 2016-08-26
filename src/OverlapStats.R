###############################################################
# Jim Kaminski
# 4/11/2016




rm(list=ls())
print(.libPaths())
#source("http://bioconductor.org/biocLite.R")
.libPaths("/home/eecs/jimkaminski/R_pkgs")

library(optparse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

#/data/yosef/CD8_effector_diff/src/YosefCode/packages/OverlapAnnotationsWithBED/OverlapStats.R 
#--ol_table /data/yosef/CD8_effector_diff/src/YosefCode/packages/OverlapAnnotationsWithBED/overlaps.bed 
#--ol_summary /data/yosef/CD8_effector_diff/src/YosefCode/packages/OverlapAnnotationsWithBED/stats.tab_intermediate_tmp 
#--shinyout /data/yosef/CD8_effector_diff/src/YosefCode/packages/OverlapAnnotationsWithBED/stats.tab --rowmeta 'Disease'

option_list <- list(
  make_option("--ol_table", default="/data/yosef/CD8_effector_diff/src/YosefCode/packages/OverlapAnnotationsWithBED/overlaps.bed", type="character",
              help="The full overlap table."),
  make_option("--ol_summary", default="/data/yosef/CD8_effector_diff/src/YosefCode/packages/OverlapAnnotationsWithBED/stats.tab_intermediate_tmp", type="character",
              help="Summary Information for each annotation."),
  make_option("--shinyout", default="/data/yosef/CD8_effector_diff/src/YosefCode/packages/OverlapAnnotationsWithBED/stats.tab", type="character",
              help="Output file for Shiny."),
  make_option("--rowmeta", default="Disease", type="character",
              help="Labels of columns in row metadata") 
)
opt <- parse_args(OptionParser(option_list=option_list))

##############################################################
GetFoldEnrichment<- function(dfOL,dfSummary,strOLName){
  print(strOLName)
  iOLPeaksBP <- sum(dfOL[,strOLName])
  iAnnotSizeBP <- dfSummary[dfSummary$OLName==strOLName,"Size"]
  iGenomeSizeBP <- dfSummary[dfSummary$OLName==strOLName,"EffectiveGenome"]
  vPeakSize <- dfOL$end - dfOL$start
  iPeakSizeBP <- sum(vPeakSize)
  
  dFoldEnrichment <- (iOLPeaksBP/(iPeakSizeBP*1.0)) / (iAnnotSizeBP/(iGenomeSizeBP*1.0))
  return(dFoldEnrichment)
}

GetStatOnScore<- function(dfOL,dfSummary,strOLName,strStat="mean"){
  iOLPeaksBP <- sum(dfOL[,strOLName])
  iAnnotSizeBP <- dfSummary[dfSummary$OLName==strOLName,"Size"]
  iGenomeSizeBP <- dfSummary[dfSummary$OLName==strOLName,"EffectiveGenome"]
  vPeakSize <- dfOL$end - dfOL$start
  iPeakSizeBP <- sum(vPeakSize)
  
  dFoldEnrichment <- (iOLPeaksBP/(iPeakSizeBP*1.0)) / (iAnnotSizeBP/(iGenomeSizeBP*1.0))
  return(dFoldEnrichment)
}

###############################################################
# Function Tests
# Mouse - GetFoldEnrichment(dfOL,dfSummary,"OregAnno_OverlapBP",2620345972) non-N bases in mm9
# Human - GetFoldEnrichment(dfOL,dfSummary,"OregAnno_OverlapBP",2897310462) non-N bases in hg19
###############################################################

#Clean up rowmeta data, remove any quotation marks
vstrRowMeta <- str_split(opt$rowmeta,pattern=" ")[[1]]
vstrRowMeta <- gsub("\'","",vstrRowMeta)

print("Loading data...")
dfOL <- read.table(opt$ol_table,sep="\t",header=TRUE)
dfSummary <- read.table(opt$ol_summary,sep="\t",header=FALSE)
names(dfSummary) <- c(vstrRowMeta,c("File","Type","Alias","Size","ValCol","Thresh","OLName","EffectiveGenome"))
print(head(dfSummary))





dfSummary <- dfSummary[,(c(vstrRowMeta,"Type","Alias","Size","Thresh","OLName","EffectiveGenome"))]

print("Calculating Enrichment...")
print(colnames(dfSummary))
print(colnames(dfOL))
matFoldEnrichment <- as.matrix(sapply(as.character(dfSummary[dfSummary$Type=="bed_overlap","OLName"]),GetFoldEnrichment,dfOL=dfOL,dfSummary=dfSummary,simplify="array"))


print("Calculating Statistics on score columns...")

vScoreCols <- as.character(dfSummary[dfSummary$Type=="bed_score","OLName"])
dfOL[,vScoreCols] <- apply(dfOL[,vScoreCols],2,as.numeric)
matMeanStat <- as.matrix(apply(dfOL[,vScoreCols],2,mean,na.rm=TRUE))
               

matAllStats <- rbind(matFoldEnrichment,matMeanStat)
dfStats <- as.data.frame(matAllStats)
colnames(dfStats) <- "Statistic"
dfStats$OLName <- row.names(dfStats)

dfFinal <- merge(dfSummary,dfStats,by="OLName")

write.table(dfFinal,sep="\t",opt$shinyout,row.names=F,col.names=T)
