###############################################################
# Jim Kaminski
# Updated on 2/23/2017
###############################################################

###############################################################
# Load Packages, Parse Command Line Args
###############################################################

rm(list=ls())

# I want to drop this reference to my library because it is 
# very poor form. However, this may wind up killing people's
# projects. Long term solution is to make sure these packages
# are available on cluster.

.libPaths("/home/eecs/jimkaminski/R_pkgs") 


print(.libPaths())

library(tidyverse)
library(optparse)
library(ggplot2)
#library(reshape2)
library(stringr)

option_list <- list(
  make_option("--ol_table", default="overlaps.tab", type="character",
              help="The full overlap table."),
  make_option("--ol_summary", default="summary_overlaps.tab_intermediate_tmp", type="character",
              help="Summary Information for each annotation."),
  make_option("--shinyout", default="summary_overlaps.tab", type="character",
              help="Output file for Shiny."),
  make_option("--rowmeta", default="empty", type="character",
              help="Labels of columns in row metadata") 
)
opt <- parse_args(OptionParser(option_list=option_list))


# Sample data that can be used for debugging.
# opt$ol_table <- "/home/eecs/jimkaminski/test_doit/out/peaks/chip/chip__MouseChIP_IP__0.05/chip__MouseChIP_IP__0.05_overlaps.tab"
# opt$ol_summary <- "/home/eecs/jimkaminski/test_doit/out/peaks/chip/chip__MouseChIP_IP__0.05/chip__MouseChIP_IP__0.05_stats.tab_intermediate_tmp" 
# opt$shinyout <- "/home/eecs/jimkaminski/test_doit/out/peaks/chip/chip__MouseChIP_IP__0.05/chip__MouseChIP_IP__0.05_stats.tab" 
# opt$rowmeta  <- "Sample QValue Lambda Duplication Broad"

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
# Calculate Overlap Statistics
##############################################

print("Calculating Enrichment...")
iUserSet_BP <- sum(dfOL$end - dfOL$start)
iUserSet_Peaks <- nrow(dfOL)


vstrOLAnnots <- unlist(as.character(dfSummary$OLName[dfSummary$Type=="bed_overlap"])) 


# For each annotation, sum up the overlap in terms of BP and Peaks.
# tbOLStats has a row for each annotation: [Annotation,TotalOverlap_BP,TotalOverlap_Peaks,UserSet_BP,UserSet_Peaks] 
tbOLStats <- dfOL %>% select(one_of(vstrOLAnnots)) %>% 
  gather(Annotation,Overlap) %>% 
  group_by(Annotation) %>%   
  summarize(TotalOverlap_BP = sum(Overlap),
            TotalOverlap_Peaks = sum(Overlap>0)) %>%
  mutate(UserSet_BP =iUserSet_BP,UserSet_Peaks = iUserSet_Peaks  )
  

# Take tbOLStats, and merge on Effective_GenomeSize_BP, and TotalOverlap_BP
tbOLStats <- dfSummary %>% select(Annotation=OLName,Annotation_BP,Annotation_Peaks,EffectiveGenomeSize_BP=EffectiveGenome) %>%
  right_join(tbOLStats,by=("Annotation")) %>% mutate(PeaksInWorld = round(EffectiveGenomeSize_BP/(Annotation_BP/Annotation_Peaks)))

# Calculate All Statistics of Interest
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
# Finally, append the user-supplied metadata to the left side.

dfMetadata <- dfSummary[rep(1,times=nrow(tbOLStats)),vstrRowMeta]
tbOLStats <- cbind(dfMetadata,tbOLStats)

# Print out the results

write.table(tbOLStats,opt$shinyout,sep="\t",row.names=F)
