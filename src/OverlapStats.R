###############################################################
# Jim Kaminski
# Updated on 1/13/2017
###############################################################

rm(list=ls())
.libPaths("/home/eecs/jimkaminski/R_pkgs")
print(.libPaths())

library(tidyverse)
library(optparse)
library(ggplot2)
library(reshape2)
library(stringr)

option_list <- list(
  make_option("--ol_table", default="overlaps.tab", type="character",
              help="The full overlap table."),
  make_option("--ol_summary", default="summary_overlaps.tab_intermediate_tmp", type="character",
              help="Summary Information for each annotation."),
  make_option("--shinyout", default="summary_overlaps.tab", type="character",
              help="Output file for Shiny."),
  make_option("--rowmeta", default="ENCODEmarks", type="character",
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

###################################################
# Clean Up Data
###################################################

#Clean up rowmeta data, remove any quotation marks
vstrRowMeta <- str_split(opt$rowmeta,pattern=" ")[[1]]
vstrRowMeta <- gsub("\'","",vstrRowMeta)

print("Loading data...")
dfOL <- read.table(opt$ol_table,sep="\t",header=TRUE)
dfSummary <- read.table(opt$ol_summary,sep="\t",header=FALSE)
names(dfSummary) <- c(vstrRowMeta,c("File","Type","Alias","Annotation_BP","Annotation_Peaks","Thresh","Score","OLName","EffectiveGenome"))
dfSummary <- dfSummary[,(c(vstrRowMeta,"Type","Alias","Annotation_BP","Annotation_Peaks","OLName","EffectiveGenome"))]

##############################################
# Calculate Statistics
##############################################

print("Calculating Enrichment...")
iUserSet_BP <- sum(dfOL$end - dfOL$start)
iUserSet_Peaks <- nrow(dfOL)

# Create a dataframe of Statistics where the dfSummary$Type=="bed_overlap"
vstrOLAnnots <- unlist(as.character(dfSummary$OLName[dfSummary$Type=="bed_overlap"])) 

tbOLStats <- dfOL %>% select(one_of(vstrOLAnnots)) %>% 
  gather(Annotation,Overlap) %>% 
  group_by(Annotation) %>%
  summarize(TotalOverlap_BP = sum(Overlap),
            TotalOverlap_Peaks = sum(Overlap>0)) %>%
  mutate(UserSet_BP =iUserSet_BP,UserSet_Peaks = iUserSet_Peaks  )
  

# Merge on effective genome size
tbOLStats <- dfSummary %>% select(Annotation=OLName,Annotation_BP,Annotation_Peaks,EffectiveGenomeSize_BP=EffectiveGenome) %>%
  right_join(tbOLStats,by=("Annotation")) %>% mutate(PeaksInWorld = round(EffectiveGenomeSize_BP/(Annotation_BP/Annotation_Peaks)))

# Calculate Statistics

tbOLStats <- tbOLStats %>% mutate(Enrichment_Peaks =(TotalOverlap_Peaks/UserSet_Peaks)/(Annotation_Peaks/PeaksInWorld) ,
                                  Enrichment_BP =  (TotalOverlap_BP/UserSet_BP)/(Annotation_BP/EffectiveGenomeSize_BP)) %>%
  mutate(Hypergeometric_Peaks_PVal_OneSided = 1-phyper(q=TotalOverlap_Peaks,      # black balls picked
                                                                  m=Annotation_Peaks,        # black balls in urn
                                                                  n=PeaksInWorld-Annotation_Peaks, # white balls in urn
                                                                  k=UserSet_Peaks), # Balls we picked
                                  
                                  Hypergeometric_BP_PVal_OneSided = 1-phyper(q=TotalOverlap_BP,      # black balls picked
                                                             m=Annotation_BP,        # black balls in urn
                                                             n=EffectiveGenomeSize_BP-Annotation_BP, # white balls in urn
                                                             k=UserSet_BP)) 


write.table(tbOLStats,opt$shinyout,sep="\t",row.names=F)


# vstrOLAnnots <- as.character(unlist(dfSummary[dfSummary$Type=="bed_overlap","OLName"]))
# TotalPeaksInUniverse <- nrow(dfOL)
# TotalBPInUniverse <- sum(dfOL$end - dfOL$start)
# tbOL <- as_tibble(dfOL) %>% select(one_of(c(vstrOLAnnots)))
# tbStats <- tbOL %>% gather(key=Annotation,OverlapInBP,everything()) %>%
#   group_by(Annotation) %>%
#   summarize(
#     Peaks_That_OL_Annotation = sum(OverlapInBP>0),
#     BP_That_OL_Annotation = sum(OverlapInBP),
#     Total_Peaks_In_Universe = TotalPeaksInUniverse,
#     Total_BP_In_Universe = TotalBPInUniverse
#   ) %>%
#   left_join(dfSummary,c("Annotation" = "OLName")) %>%
#   mutate(FoldEnrichment_BP = (BP_That_OL_Annotation/(Total_BP_In_Universe*1.0)) / (Size/(EffectiveGenome*1.0)))
#          
#          
# write.table(tbStats,"SampleResults.tab",sep="\t",row.names = FALSE)
# 
# # Old Way
# #matFoldEnrichment <- as.matrix(sapply(as.character(dfSummary[dfSummary$Type=="bed_overlap","OLName"]),GetFoldEnrichment,dfOL=dfOL,dfSummary=dfSummary,simplify="array"))
# 
# 
# 
# print("Calculating Statistics on score columns...")
# 
# vScoreCols <- as.character(dfSummary[dfSummary$Type=="bed_score","OLName"])
# dfOL[,vScoreCols] <- apply(dfOL[,vScoreCols],2,as.numeric)
# matMeanStat <- as.matrix(apply(dfOL[,vScoreCols],2,mean,na.rm=TRUE))
#                
# 
# matAllStats <- rbind(matFoldEnrichment,matMeanStat)
# dfStats <- as.data.frame(matAllStats)
# colnames(dfStats) <- "Statistic"
# dfStats$OLName <- row.names(dfStats)
# 
# dfFinal <- merge(dfSummary,dfStats,by="OLName")
# 
# write.table(dfFinal,sep="\t",opt$shinyout,row.names=F,col.names=T)
