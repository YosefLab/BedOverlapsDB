##############################################################################
# Jim Kaminski
# 4/8/2016
# Yosef Lab
##############################################################################
"""
BedOverlaps DB

This pipeline takes a set of genomic itervals and measures their enrichment in 
a set of genomic annotations. For example, one could use the pipeline to 
determine if a set of ATAC-Seq peaks are enriched for regulatory features in 
the mouse genome. "Enriched" means that these peaks have a higher prevalence of 
some regulatory feature than is found in the organism's full genome.

Sample Call

    python OverlapAnnotationsWithBed.py --annotation_list example_data/example_annotations.tab  \
    --input_bed example_data/example_peaks.bed --additional_text_labels SampleName \
     --additional_text test_sample --out_overlaps my_test_overlap_matrix.tab --out_summary my_test_stat_summary.tab


Input
** Bed file (chr, start, end)
** Tab file of annotations (file,type[bedscore,bed,Nirdb],[score column],[score thresh],[annotation size])

Output (Three Tables)
** Our Main Overlaps Table (This can be simple bedtools output, a series of left outer joins)
** A File Of Statistics On The Final Data (cut the overlaps table, build this in R)
** Tidy csv file for Shiny (Construct in R, pass macs2qval, maxfragsize, and sample name as args)
"""
import argparse
import glob
import sys
import os
import subprocess as sp
import pandas as pd


parser = argparse.ArgumentParser(description='')
grpParam = parser.add_argument_group('')
grpParam.add_argument('--annotation_list', type=str, dest='tabAnnotList',default="", 
                      help='Please enter the full path your list of annotations. See hg19_annotations.tab for an example.')
grpParam.add_argument('--out', type=str, dest='dirOut', help='Output folder',default=os.getcwd()+os.sep+"out")
grpParam.add_argument('--genome', type=str, dest='strRefGenome', help='Please enter hg19,mm9, etc.',default="")
grpParam.add_argument('--input_bed', type=str, dest='bedIn',help='Enter the path to your bedfile.',default= "")
grpParam.add_argument('--out_overlaps', type=str, dest='tabOut',help='Enter the path to output the final file.',default= "overlaps.tab")
grpParam.add_argument('--out_summary', type=str, dest='tabSummaryOut',help='Enter the path to output the summary file.',default= "summary_overlaps.tab")
grpParam.add_argument('--additional_text', type=str, dest='strRowStart',help='This is additional text that you can prepend to every row of the summary data.',default= "",nargs='+')
grpParam.add_argument('--additional_text_labels', type=str, dest='astrRowLabels',help='These are the variable names for the values specified in additional_text.',default= "",nargs='+')
grpParam.add_argument('--debug', dest='bDebug',help='Adding this flag keeps the temporary files this program makes.', const=True, default = False,action='store_const')

args = parser.parse_args()

print("Checking script directory")
dirOverlaps = os.path.dirname(os.path.realpath(__file__))
print(dirOverlaps)

###############################################################################
# Functions
###############################################################################


def check_create_dir( strDir ):
    """ Creates a directory if it does not already exist."""
    if not os.path.exists( strDir ):
        os.makedirs( strDir )
    return strDir
    

def CheckForWindowsCRLF(strFile, strFind="\r\n", iLines=5):
    """ Looks for strFind in strFile. Main purpose is to find "\n\r" files which 
        are not compatiable with  bedtools. """
        
    bFoundIt = False    
    iLineCount = 0
    with open(strFile,'r') as fileRead:
        for strLine in fileRead:
            iLineCount+=1
            
            if strFind in strLine:
                bFoundIt=True
                break
            if iLineCount >= iLines:
                break
            
    if(bFoundIt):
        raise Exception('The file ' + strFile + ' is a Windows file, which will not' + 
        ' work with bedtools. Please run "dos2unix" on the file, and try the OL' +
        ' database again.')

def extend_list(aIn,iCorrectLen):
    """ This function adds empty elements to a list if it's too short. It is
    called when loading the metadata to make empty cells for missing columns. """
    if len(aIn) < iCorrectLen:
        aIn = aIn + (iCorrectLen -len(aIn))*[""]
    return aIn

def CheckNumber(strIn):
    """ Checks to see if "strIn" is a number. If not, it will return -1. We use
    this when we read in metadata on the different genomic annotations. """
    try:
        float(strIn)
        return float(strIn)
    except:
        return float(-1)
        
def SanitizeForR(strIn):
    """ We replace '-' with "." because R does not accept hyphens/minus signs
    in column headers. """
    
    strOut = strIn.replace("-",".")

    return strOut

def ExtractAnnotVarsFromLine(astrLine):
    """ This function reads in a line of metadata from the annotation list,
    formats/sanitizes each element of the line for processing, and then returns
    the elements.
    
    By format/sanitize, I mean things like converting strings instances of
    numbers to integers, etc.
    
    """

    # 0=File
    # 1=Type
    strFile = astrLine[0].strip()
    strType = astrLine[1].strip()

    # 2= Size        
    if len(astrLine)>2: 
        iSize = CheckNumber(astrLine[2])
    else:
        iSize = -1.0
    # 3= NumPeaks        
    if len(astrLine)>3: 
        iAnnotPeaks = CheckNumber(astrLine[3])
    else:
        iAnnotPeaks = -1.0
    # 4= Alias
    if len(astrLine)>4:
        strAlias = astrLine[4].strip()
    else:
        strAlias = strFile.split(os.sep)[-1].replace(".bed","")

    # 5= ScoreCol
    if len(astrLine)>5:
        iScoreCol = CheckNumber(astrLine[5])
    else:
        iScoreCol = -1
        
    #6= function
    if len(astrLine)>6:
        strFunc = astrLine[6]
    else:
        strFunc = ""
    
    # 7= Thresh
    if len(astrLine)>7:
        dThresh = CheckNumber(astrLine[7])
    else:
        dThresh = -1.0
        
    # 8= Effective Genome
    if len(astrLine)>8:
        iEffectiveGenome = CheckNumber(astrLine[8])
    else:
        iEffectiveGenome = -1.0
        
    # Replacing -'s with .'s because R does not like minus signs.
    strAlias = strAlias.replace("-",".")
        
    return strFile,strType,int(iSize),iAnnotPeaks,strAlias,int(iScoreCol),dThresh,strFunc,iEffectiveGenome
    
def GetAnnotationSize(bedAnnotFile,bedTmp,cmdBedtools):
    """ Takes in a bed file of an annotation, merges overlapping intervals,
    and returns the total lenghts of the intervals in bp (iSize), and the number
    of peaks/intervals (iPeaks).
    """
    iSize = 0
    iPeaks = 0
    
    with open(bedTmp, 'w') as fileTmp:
        
        pCut = sp.Popen(["cut","-f1-3",bedAnnotFile],stdout=sp.PIPE)
        stdoutCut = pCut.communicate()[0]    
        
        pSort = sp.Popen(["sort","-k","1,1","-k2,2n"],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutSort = pSort.communicate(stdoutCut)[0]
            
        pMerge = sp.Popen([cmdBedtools,"merge","-i","stdin"],stdin=sp.PIPE,stdout=fileTmp)
        pMerge.communicate(stdoutSort)[0]    
    
    # Count intervals/peaks, and sum size of intervals/peaks
    with open(bedTmp, 'r') as fileTmp:
        for strLine in fileTmp:
            iPeaks+=1
            astrLine = strLine.split('\t')
            if strLine.strip() != "" and strLine[0]!="#":
                iSize+= (int(astrLine[2]) - int(astrLine[1]))
    return iSize,iPeaks
    
    
            
def AddHeaderToFile(strFileIn,strFileOut,astrHeader):
    """ Adds *astrheader* to the top of *strFileIn* and outputs it to
    *strFileOut*.
    """
    with open(strFileOut,'w') as fileOut:
        with open(strFileIn,'r') as fileIn:
            fileOut.write('\t'.join(astrHeader)+"\n")
            
            for strLine in fileIn:
                fileOut.write(strLine.strip()+"\n")
                
def GetColsInBedFile(bedIn):
    """ Counts up the number of columns in a bed file. Returns -1 if there is
    an error/ no columns.
    """
    iCols=-1
    with open(bedIn,'r') as fileBed:
        for strLine in fileBed:
            if strLine[0]!="#":
                astrLine = strLine.split('\t')
                iCols=len(astrLine)
                break
        return iCols
                
def ProcessAnnotation_ClosestFeature(bedIn,bedAnnotFile,bedOut,iSizeAnnot,iFeatureCol,cmdBedtools,astrHeader):
    """
    For each row in bedIn, list the closest feature in bedAnnotFile and the distance to it.
    
    Returns:
    * dfResult - a pandas df of the original rows from bedIn with cols[chr,start,end], and the cols ['feature','dist_to_feature']
    
    Process:
    * Cut the annotation file (bedAnnotFile) to the [chr, start, end] columns
    * Sort it
    * Run "bedtools closest"
    * Cut down the results file to [chr,start,end,feature,dist_to_feature]
    * Run sed to fix instances where there is no closest feature (should be rare).
    """
    iColsInInput = GetColsInBedFile(bedIn)
    strInputCols = ",".join(map(str,range(1,iColsInInput+1)))    
    strKeepCols = strInputCols +","+str(iColsInInput+4) + "," + str(iColsInInput+5)
    
    with open(bedOut,'w') as fileBedOut:
        # Cut the annot file
        
        pCut = sp.Popen(["cut","-f1-3"+","+str(iFeatureCol),bedAnnotFile],stdout=sp.PIPE)
        stdoutCut = pCut.communicate()[0]    
        
        # Sort it
        pSort = sp.Popen(["sort","-k","1,1","-k2,2n"],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutSort = pSort.communicate(stdoutCut)[0]
    
    
        pClosest = sp.Popen([cmdBedtools,"closest", "-a",bedIn,"-b","stdin","-d"],
                          stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutClosest = pClosest.communicate(stdoutSort)[0] 
        
        pCut2 = sp.Popen(["cut","-f"+strKeepCols],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutCut2 = pCut2.communicate(stdoutClosest)[0]
        
        pSed =sp.Popen(["sed","-e",'s/\t\t\.$/\tNoFeature\t0/'],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutSed = pSed.communicate(stdoutCut2)[0]
        
        pGroup = sp.Popen([cmdBedtools,"groupby","-i","stdin","-g","1,2,3","-c","4,5","-o","distinct,min"],stdin=sp.PIPE,stdout=fileBedOut)
        pGroup.communicate(stdoutSed)[0]
        
    dfResult = pd.read_table(bedOut,header=None,names=["chr","start","end"]+astrHeader,index_col=False)
    dfResult.head()
    return dfResult
    
        
def ProcessAnnotation_ScoreBed(bedIn,bedAnnotFile,bedOut,iValCol,strFunction,iSizeAnnot,cmdBedtools,astrHeader,dirTmp):
    """
    Performs left-outer-join of bedAnnotFile on bedIn, sums bp_overlap for each row of bedIn 
    
    Returns:
    * dfResult - a pandas df of the original rows from bedIn with cols[chr,start,end], and the col ['stat_on_score_column']

    Process:
    * Sort the annotation file (bedAnnotFile) 
    * Run bedtools map - this calls the function of interest on the appropriate column.
    """

    # Sort the input file
    with open(bedOut,'w') as fileBedOut:
        
        print("Sorting annotation...")  
        pSort = sp.Popen(["sort","-k1,1","-k2,2n",bedAnnotFile],stdout=sp.PIPE)
        stdoutSort = pSort.communicate()[0]
      
        print("Running bedtools map...")
        #print [cmdBedtools,"map", "-a",bedIn,"-b","stdin","-c",str(iValCol),"-o",strFunction,"-prec","3"]
        pBedtools = sp.Popen([cmdBedtools,"map", "-a",bedIn,"-b","stdin","-c",str(iValCol),"-o",strFunction,"-null",".","-prec","2"],
                              stdin=sp.PIPE,stdout=fileBedOut)
        pBedtools.communicate(stdoutSort)[0]
        
        
    dfResult = pd.read_table(bedOut,header=None,names=["chr","start","end"]+astrHeader,index_col=False)
    dfResult.head()

    return dfResult

def ProcessAnnotation_Overlap(bedIn,bedAnnotFile,bedOut,iSizeAnnot,iPeaksAnnot,cmdBedtools,astrHeader,dirTmp):
    """
    Performs left-outer-join of bedAnnotFile on bedIn, sums bp_overlap for each row of bedIn 
    
    Returns:
    * dfResult - a pandas df of the original chr,start,end rows from bedIn, and a col of overlap with bedAnnotFile
    * iSizeAnnot - length specified in --annotation_list, if -1, the function will calculate it.    
    
    Process:
    * Cuts the annotation file (bedAnnotFile) to the [chr, start, end] columns
    * Sorts it
    * Merges the overlapping regions within bedAnnotFile
    * Intersects this data with (bedIn), and reports the number of bp overlap
    * Cuts down file to the original bedIn columns, and bp overlap col[chr,start,end,bp_overlap]
    * Sums the overlap over the row of bedIn, reduces this to unique rows of bedIn
    """

    iColsInInput = GetColsInBedFile(bedIn)
    astrInputCols = ",".join(map(str,range(1,iColsInInput+1)))

    iValCol = iColsInInput+4
    iOpCol = iColsInInput +1

    
    astrCut = astrInputCols+","+str(iValCol)    
    with open(bedOut,'w') as fileBedOut:
        # Cut the annot file
        pCut = sp.Popen(["cut","-f1-3",bedAnnotFile],stdout=sp.PIPE)
        stdoutCut = pCut.communicate()[0]    
        
        # Sort it
        pSort = sp.Popen(["sort","-k1,1","-k2,2n"],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutSort = pSort.communicate(stdoutCut)[0]
        
        # Merge overlapping intervals 
        pMerge = sp.Popen([cmdBedtools,"merge","-i","stdin"],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutMerge = pMerge.communicate(stdoutSort)[0]
        pBedtools = sp.Popen([cmdBedtools,"intersect", "-a",bedIn,"-b","stdin","-wao"],stdin=sp.PIPE,stdout=sp.PIPE)
        #stdoutIntersect= pBedtools.communicate(stdoutMerge)[0]
        stdoutIntersect=pBedtools.communicate(stdoutMerge)[0]
        
        pCut2 = sp.Popen(["cut","-f"+astrCut,],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutCut2 = pCut2.communicate(stdoutIntersect)[0]
        
        pGroup = sp.Popen([cmdBedtools,"groupby","-i","stdin","-g",astrInputCols,"-c",str(iOpCol),"-o","sum"],stdin=sp.PIPE,stdout=fileBedOut)
        pGroup.communicate(stdoutCut2)[0]

    dfResult = pd.read_table(bedOut,header=None,names=["chr","start","end"]+astrHeader,index_col=False)
    dfResult.head()        
        
    if iSizeAnnot == -1 or iPeaksAnnot == -1:
        print("Calculating Annotation Size...")
        iSizeAnnot,iPeaksAnnot = GetAnnotationSize(bedAnnotFile,dirTmp + os.sep + astrHeader[0] +"_temp.bed",cmdBedtools)
    # include awk script to look a
    return iSizeAnnot,dfResult,iPeaksAnnot

###############################################################################
# Program
###############################################################################



#############################################################################
# Create out and tmp folders.

dirOut = args.dirOut
dirTmp = args.dirOut + os.sep + "tmp_files_for_bedoverlaps"
check_create_dir( dirOut )
check_create_dir( dirTmp )
astrFilesToDelete = []

############################################################################
# Load the user's peak file (bedIn). 

bedIn = args.bedIn

# Check if it has the Windows CRLF. Files with CRLF's can't be used in bedtools.
CheckForWindowsCRLF(strFile=bedIn, strFind="\r\n", iLines=5)

# Cut input file down to first three columns and sort it.
with open(dirTmp+os.sep+"CutBed.bed",'w') as fileBedCut:
    pCut = sp.Popen(["cut","-f1-3",bedIn],stdout=sp.PIPE)
    stdoutCut = pCut.communicate()[0]    
            
    # Sort it
    pSort = sp.Popen(["sort","-k","1,1","-k2,2n"],stdin=sp.PIPE,stdout=fileBedCut)
    pSort.communicate(stdoutCut)[0]

bedIn = dirTmp+os.sep+"CutBed.bed"
    

dfOverlapsResults = pd.read_table(bedIn,names=["chr","start","end"],index_col=False)
dfOverlapsResults = dfOverlapsResults[dfOverlapsResults.chr.str.startswith("#")==False]
print(dfOverlapsResults.head())


###############################################################################
# Load the list of annotations (args.tabAnnotList)

iAnnotationsProcessed = 0
aAnnotationSummary = []
astrHeader = ["#chr","start","end"]


with open(args.tabAnnotList, 'r') as fileAnnotList:
    for strLine in fileAnnotList:
        astrLine = strLine.split('\t')

        # Confirm that annotation has all nine pieces of metadata.
        if len(astrLine)!=9 and strLine[0]!="#":
            raise Exception('Your list of annotations appears to have ' + str(len(astrLine)) +  
            ' columns. The expected number is nine. Please check and confirm that ' +
            ' your list of annotations has these columns \n' +
            ' * #File - This is the path to the bed file for annotation.\n'+
            ' * Type - Operation to carry out. For example: "bed_overlap"\n'+
            ' * Size - This is the size of the annotation in bp. You can enter "-1" if you want the program to calculate it.\n'+
            ' * NumPeaks  - This is the number of peaks in the annotation. You can set this to "-1" to have the program calculate it.\n'+
            ' * Alias - An alternative, short name for the bed file\n'+
            ' * ScoreCol - This is used when the annotation file has a score, like say evolutionary conservation. You can leave it blank if you are not doing "bed_score".\n'+
            ' * Bedtools function - Function to apply to score column. You can leave it blank if you are not doing "bed_score".\n'+
            ' * Threshold - Threshold to apply to score column. You can leave it blank if you are not doing "bed_score".\n'+
            ' * EffectiveGenomeSize - This is the effective genome size to use for the annotation.')
    
        
        # If line is not blank or commented out, add the metadata and run the 
        # appropriate operation on the file.
        if strLine.strip() != "" and strLine[0]!="#":
            
            # Extract Information from line, prep output file
            astrLine = extend_list(astrLine,9)
            strAnnotFile,strType,iAnnotSize,iAnnotPeaks,strAlias,iScoreCol,dThresh,strFunc,iEffectiveGenome = ExtractAnnotVarsFromLine(astrLine)
            
            bedNextOutput = dirTmp + os.sep + "File"+str(iAnnotationsProcessed+1) + "_" + strType + ".bed" 
            print("Processing " + os.path.basename(strAnnotFile).replace(".bed","") + " .... ")
            CheckForWindowsCRLF(strFile=strAnnotFile, strFind="\r\n", iLines=5)
            
            # Run appropriate bedtools command.          
            if strType == "bed_score":
                astrHeader = [strAlias+"_"+strFunc+"score"]
                dfAnnotResults =ProcessAnnotation_ScoreBed(bedIn,strAnnotFile,bedNextOutput,iScoreCol,strFunc,iAnnotSize,"bedtools",astrHeader,dirTmp)
            if strType == "bed_overlap":
                astrHeader = [strAlias+"_OverlapBP"]
                iAnnotSize,dfAnnotResults,iAnnotPeaks =ProcessAnnotation_Overlap(bedIn,strAnnotFile,bedNextOutput,iAnnotSize,iAnnotPeaks,"bedtools",astrHeader,dirTmp)
            if strType == "bed_closest":
                astrHeader = ["Closest_"+strAlias,"DistToClosest_"+strAlias]
                dfAnnotResults =ProcessAnnotation_ClosestFeature(bedIn,strAnnotFile,bedNextOutput,iAnnotSize,iScoreCol,"bedtools",astrHeader)
            
            # Print first few rows of the results, before and after merging on
            # the new results to dfOverlapsResults. This is just to give the 
            # user a sense of what the program is doing and can be helpful for
            # debugging.
            
            print(dfAnnotResults.head())
            print(dfOverlapsResults.head())
            dfOverlapsResults = pd.merge(dfOverlapsResults, dfAnnotResults, how='left', left_on=['chr','start','end'],right_on=['chr','start','end'])
            dfOverlapsResults = dfOverlapsResults[dfOverlapsResults.chr.str.startswith("#")==False]
            print(dfOverlapsResults.head())
            astrFilesToDelete.append(bedNextOutput)
            

            print("Effective Genome Size")
            print(str(iEffectiveGenome))
            aAnnotationSummary.append([strAnnotFile,strType,strAlias,iAnnotSize,iAnnotPeaks,iScoreCol,dThresh,astrHeader[-1],iEffectiveGenome])
            iAnnotationsProcessed+=1
            

# Write out a tab delimited file of overlap results. Each row is peak/interval
# and the columns contain the relevant information from each annotation.
dfOverlapsResults.to_csv(path_or_buf=args.tabOut, sep='\t',index=False)
            
#############################################################################
# Produce summary statistics and run hypothesis tests in R.


""" 
Format metadata supplied in "--additional_text". astrRowStart will have
a row for each annotation,and each row will beging with the text supplied to
to --additional_text.

We do this so that users can combine the summary data from R from many runs
of the bedoverlaps db. For example you may try doing peakcalling at different
levels (qvalues) and could pass the SampleName and QValue to --additional_text


#SampleName, QValue, Annotation
[Sample1, 0.005, Annotation_1, ...]
[Sample1, 0.005, Annotation_2, ...]
[Sample1, 0.010, Annotation_1, ...]
[Sample1, 0.010, Annotation_2, ...]
[Sample2, 0.005, Annotation_1, ...]
[Sample2, 0.005, Annotation_2, ...]
[Sample2, 0.010, Annotation_1, ...]
[Sample2, 0.010, Annotation_2, ...]

"""

if args.strRowStart == "":
    astrRowStart = []
else:
    astrRowStart = []
    # Unpack text values
    aastrData = [x.split() for x in args.strRowStart]
    for astrElement in aastrData:
        for strElement in astrElement:
            astrRowStart.append(strElement)
    
# Print out the summary statistics
print("Outputting summary stats to:" +  args.tabSummaryOut+"_intermediate_tmp")
with open(args.tabSummaryOut+"_intermediate_tmp",'w') as fileSummary:
  for astrLine in aAnnotationSummary:
    fileSummary.write( "\t".join(astrRowStart + [str(x) for x in astrLine]).strip()+"\n")  


# Run RScript
astrCMD = ["Rscript", dirOverlaps + "/src/OverlapStats.R",
         "--ol_table",args.tabOut,
         "--ol_summary",args.tabSummaryOut+"_intermediate_tmp",
         "--shinyout",args.tabSummaryOut,
         "--rowmeta","\'" + " ".join(args.astrRowLabels)+ "\'"]
print(" ".join(astrCMD))
sp.call(astrCMD)


###############################################################################
# Delete temporary files.

if args.bDebug==False:
    for strFile in astrFilesToDelete:
        os.remove(strFile)
    
    sys.stderr.write("Cleaning up... \n")
    for strFile in glob.glob(dirTmp + "/*"):
        sys.stderr.write("Deleting " + strFile + "\n")
        os.remove(strFile)
        
    os.rmdir(dirTmp)

    
