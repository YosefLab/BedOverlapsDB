##############################################################################
# Jim Kaminski
# 4/8/2016
# Yosef Lab
##############################################################################
"""
CUT THE INPUT FILE, AND ALL BED_OVERLAPS FILES TO BE


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
import numpy as np

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
grpParam.add_argument('--additional_text_labels', type=str, dest='astrRowLabels',help='This is additional text that you can prepend to every row of the summary data.',default= "",nargs='+')
grpParam.add_argument('--debug', dest='bDebug',help='Adding this flag keeps the temporary files this program makes.', const=True, default = False,action='store_const')

args = parser.parse_args()

# get Folder script is saved in
# http://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
print("Checking script directory")
dirOverlaps = os.path.dirname(os.path.realpath(__file__))
print(dirOverlaps)

#############################################################################
# Functions

# Creates a directory if it does not already exist.
def check_create_dir( strDir ):
    if not os.path.exists( strDir ):
        os.makedirs( strDir )
    return strDir
    
# Looks for strFind in strFile. Main purpose is to find "\n\r" files which 
# upset bedtools.
def CheckForWindowsCRLF(strFile, strFind="\r\n", iLines=5):
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

# Adds empty elements to a list if it's too short.
def extend_list(aIn,iCorrectLen):
    if len(aIn) < iCorrectLen:
        aIn = aIn + (iCorrectLen -len(aIn))*[""]
    return aIn

def CheckNumber(strIn):
    try:
        float(strIn)
        return float(strIn)
    except:
        return float(-1)
        
def SanitizeForR(strIn):
    strOut = strIn.replace("-",".")
    # check if first character starts with a numnber
    return strOut

def ExtractAnnotVarsFromLine(astrLine):
    # We indicate missing values with negative numbers.

    # 0=File
    # 1=Type
    strFile = astrLine[0].strip()
    strType = astrLine[1].strip()

    # 2= Size        
    if len(astrLine)>2: 
        iSize = CheckNumber(astrLine[2])
    else:
        iSize = -1.0
    # 3= Alias
    if len(astrLine)>3:
        strAlias = astrLine[3].strip()
    else:
        strAlias = strFile.split(os.sep)[-1].replace(".bed","")

    # 4= ScoreCol
    if len(astrLine)>4:
        iScoreCol = CheckNumber(astrLine[4])
    else:
        iScoreCol = -1
        
    #5= function
    if len(astrLine)>5:
        strFunc = astrLine[5]
    else:
        strFunc = ""
    
    # 6= Thresh
    if len(astrLine)>6:
        dThresh = CheckNumber(astrLine[6])
    else:
        dThresh = -1.0
        
    # 7= Effective Genome
    if len(astrLine)>7:
        iEffectiveGenome = CheckNumber(astrLine[7])
    else:
        iEffectiveGenome = -1.0
        
    # Replacing -'s with .'s because R does not like minus signs.
    strAlias = strAlias.replace("-",".")
        
    return strFile,strType,int(iSize),strAlias,int(iScoreCol),dThresh,strFunc,iEffectiveGenome
    
def GetAnnotationSize(bedAnnotFile,bedTmp,cmdBedtools):
    iSize = 0
    
    with open(bedTmp, 'w') as fileTmp:
        
        pCut = sp.Popen(["cut","-f1-3",bedAnnotFile],stdout=sp.PIPE)
        stdoutCut = pCut.communicate()[0]    
        
        pSort = sp.Popen(["sort","-k","1,1","-k2,2n"],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutSort = pSort.communicate(stdoutCut)[0]
            
        pMerge = sp.Popen([cmdBedtools,"merge","-i","stdin"],stdin=sp.PIPE,stdout=fileTmp)
        pMerge.communicate(stdoutSort)[0]    
    
    
    with open(bedTmp, 'rb') as fileTmp:
        for strLine in fileTmp:
            astrLine = strLine.split('\t')
            if strLine.strip() != "" and strLine[0]!="#":
                iSize+= (int(astrLine[2]) - int(astrLine[1]))
    return iSize
    
    
            
def AddHeaderToFile(strFileIn,strFileOut,astrHeader):
    with open(strFileOut,'w') as fileOut:
        with open(strFileIn,'rb') as fileIn:
            fileOut.write('\t'.join(astrHeader)+"\n")
            
            for strLine in fileIn:
                fileOut.write(strLine.strip()+"\n")
                
                
            
            
            
def GetColsInBedFile(bedIn):
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
        
        print "Sorting annotation..."  
        pSort = sp.Popen(["sort","-k1,1","-k2,2n",bedAnnotFile],stdout=sp.PIPE)
        stdoutSort = pSort.communicate()[0]
      
        print "Running bedtools map..."
        #print [cmdBedtools,"map", "-a",bedIn,"-b","stdin","-c",str(iValCol),"-o",strFunction,"-prec","3"]
        pBedtools = sp.Popen([cmdBedtools,"map", "-a",bedIn,"-b","stdin","-c",str(iValCol),"-o",strFunction,"-null",".","-prec","2"],
                              stdin=sp.PIPE,stdout=fileBedOut)
        pBedtools.communicate(stdoutSort)[0]
        
        
    dfResult = pd.read_table(bedOut,header=None,names=["chr","start","end"]+astrHeader,index_col=False)
    dfResult.head()

    return dfResult

def ProcessAnnotation_Overlap(bedIn,bedAnnotFile,bedOut,iSizeAnnot,cmdBedtools,astrHeader,dirTmp):
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
        
    if iSizeAnnot == -1:
        print "Calculating Annotation Size..."
        iSizeAnnot = GetAnnotationSize(bedAnnotFile,dirTmp + os.sep + astrHeader[0] +"_temp.bed",cmdBedtools)
    # include awk script to look a
    return iSizeAnnot,dfResult

#############################################################################

dirOut = args.dirOut
dirTmp = args.dirOut + os.sep + "tmp_files_for_bedoverlaps"
check_create_dir( dirOut )
check_create_dir( dirTmp )
astrFilesToDelete = []

############################################################################
# Read in the list of files, proces each one in order.

iAnnotationsProcessed = 0
astrHeader = ["#chr","start","end"]
aAnnotationSummary = []
bedIn = args.bedIn

## Check input file for issues:
CheckForWindowsCRLF(strFile=bedIn, strFind="\r\n", iLines=5)



## Cut and sort the input file

with open(dirTmp+os.sep+"CutBed.bed",'w') as fileBedCut:
    pCut = sp.Popen(["cut","-f1-3",bedIn],stdout=sp.PIPE)
    stdoutCut = pCut.communicate()[0]    
            
    # Sort it
    pSort = sp.Popen(["sort","-k","1,1","-k2,2n"],stdin=sp.PIPE,stdout=fileBedCut)
    pSort.communicate(stdoutCut)[0]

bedIn = dirTmp+os.sep+"CutBed.bed"
    


dfOverlapsResults = pd.read_table(bedIn,names=["chr","start","end"],index_col=False)
#dfOverlapsResults.columns = ["test1","test2","test3"]

print dfOverlapsResults.head()



with open(args.tabAnnotList, 'rb') as fileAnnotList:
    for strLine in fileAnnotList:
        astrLine = strLine.split('\t')
        
        if strLine.strip() != "" and strLine[0]!="#":
            
            # Extract Information from line, prep output file
            astrLine = extend_list(astrLine,7)
            strAnnotFile,strType,iAnnotSize,strAlias,iScoreCol,dThresh,strFunc,iEffectiveGenome = ExtractAnnotVarsFromLine(astrLine)
            bedNextOutput = dirTmp + os.sep + "File"+str(iAnnotationsProcessed+1) + "_" + strType + ".bed" 
            print "Processing " + os.path.basename(strAnnotFile).replace(".bed","") + " ...."
            CheckForWindowsCRLF(strFile=strAnnotFile, strFind="\r\n", iLines=5)
            # Run appropriate bedtools command to merge on output           
            if strType == "bed_score":
                astrHeader = [strAlias+"_"+strFunc+"score"]
                dfAnnotResults =ProcessAnnotation_ScoreBed(bedIn,strAnnotFile,bedNextOutput,iScoreCol,strFunc,iAnnotSize,"bedtools",astrHeader,dirTmp)
            if strType == "bed_overlap":
                astrHeader = [strAlias+"_OverlapBP"]
                iAnnotSize,dfAnnotResults =ProcessAnnotation_Overlap(bedIn,strAnnotFile,bedNextOutput,iAnnotSize,"bedtools",astrHeader,dirTmp)
            if strType == "bed_closest":
                astrHeader = ["Closest_"+strAlias,"DistToClosest_"+strAlias]
                dfAnnotResults =ProcessAnnotation_ClosestFeature(bedIn,strAnnotFile,bedNextOutput,iAnnotSize,iScoreCol,"bedtools",astrHeader)
            
           # print dfOverlapsResults.head()
            print dfAnnotResults.head()
            dfOverlapsResults = pd.merge(dfOverlapsResults, dfAnnotResults, how='left', left_on=['chr','start','end'],right_on=['chr','start','end'])
            dfOverlapsResults.head()
            astrFilesToDelete.append(bedNextOutput)
            
            # Delete file if you are not in debug mode
            print "Effective Genome Size"
            print str(iEffectiveGenome)
            aAnnotationSummary.append([strAnnotFile,strType,strAlias,iAnnotSize,iScoreCol,dThresh,astrHeader[-1],iEffectiveGenome])
            iAnnotationsProcessed+=1
            
#AddHeaderToFile(bedNextOutput,args.tabOut,astrHeader)

dfOverlapsResults.to_csv(path_or_buf=args.tabOut, sep='\t',index=False)
            
#############################################################################
# Process results in R

if args.strRowStart == "":
    astrRowStart = []
else:
    astrRowStart = []
    # Unpack text values
    aastrData = [x.split() for x in args.strRowStart]
    for astrElement in aastrData:
        for strElement in astrElement:
            astrRowStart.append(strElement)
    

print "Outputting summary stats to:" +  args.tabSummaryOut+"_intermediate_tmp"
with open(args.tabSummaryOut+"_intermediate_tmp",'w') as fileSummary:
  for astrLine in aAnnotationSummary:
    fileSummary.write( "\t".join(astrRowStart + [str(x) for x in astrLine]).strip()+"\n")  

#astrFilesToDelete.append(args.tabSummaryOut+"intermediate")


astrCMD = ["Rscript", dirOverlaps + "/src/OverlapStats.R",
         "--ol_table",args.tabOut,
         "--ol_summary",args.tabSummaryOut+"_intermediate_tmp",
         "--shinyout",args.tabSummaryOut,
         "--rowmeta","\'" + " ".join(args.astrRowLabels)+ "\'"]
print " ".join(astrCMD)
sp.call(astrCMD)

#astrFilesToDelete.append(args.tabSummaryOut+"_intermediate_tmp")
# Cleanup
for strFile in astrFilesToDelete:
    os.remove(strFile)

sys.stderr.write("Cleaning up... \n")
for strFile in glob.glob(dirTmp + "/*"):
    sys.stderr.write("Deleting " + strFile + "\n")
    os.remove(strFile)
    
os.rmdir(dirTmp)

    
