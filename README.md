BedOverlapsDB
===========

This pipeline takes a set of genomic itervals and measures their enrichment in a set of genomic annotations. For example, one could use the pipeline to determine if a set of ATAC-Seq peaks are enriched for regulatory features in the mouse genome. "Enriched" means that these peaks have a higher prevalence of some regulatory feature than is found in the organism's full genome.

Jim created and maintains this pipeline. 
Please note that this documentation is a work in progress.

---

#### Table of Contents
[Important Notes](#important_notes)    
[Quickstart](#quickstart)    

#### Overview
[Introduction](#intro)    
[Input Files](#input_files)    
[Output Files](#output_files)    

##### Common Tasks   
[bed_overlap: Find Overlap and Enrichment](#bed_overlap)  
[bed_closest: Find Nearest Features](#bed_closest)  
[bed_score: Produce a Summary Statistic on Score Column ](#bed_score)     

---

<a name="important_notes"></a>

## Important Notes

* The pipeline currently relies on Jim's directory of R packages to load some things we don't have installed lab-wide yet. That's fine as a temporary measure, but we should install the "tidyverse" on s121-s124 when there is time.
* Typically, the pipeline merges overlapping peaks within a set of peaks to avoid double-counting. I will check and confirm how this is handled for the input set, and an annotation set.

<a name="quickstart"></a>

## Quickstart

To make sure you can use this repository, please run the following test from inside the BedOverlapsDB folder:

```
python OverlapAnnotationsWithBed.py --annotation_list example_data/example_annotations.tab  \
--input_bed example_data/example_peaks.bed --additional_text_labels SampleName \
 --additional_text test_sample --out_overlaps my_test_overlap_matrix.tab --out_summary my_test_stat_summary.tab
```

That will check the example_peaks.bed against a set of small sample annotations included with the package, and carry out the three basic operations performed by the package (measure overlap, find nearest feature, calculate a stat based on a numeric value in an annotation). It's a good first test to make sure that you have access to all the packages and programs you need to run BedOverlapsDB.

<a name="intro"></a>

## Introduction

Many of us work with data from ATAC-Seq, DNase-Seq and other experimental techniques that identify particular regions of the genome that we expect to have some functional purpose. ATAC-Seq, for example, generates reads from open chromatin regions, and these regions typically contain things like enhancers. Some of the questions you might have as a researcher are:

* Are these regions enriched for widely used lists of regulatory regions?
* Are these regions close to transcription start sites (TSS's)? Are they near any genes we're interested in?
* Do these regions have a high average score for some measure of conservation?

The BedOverlapsDB is a tool for answering these questions. It answers the first question by measuring enrichment for a set of features. It answers the second by returning closest feature (and distance to it) from a set of features. It answers the last question by returning a summary statistic (mean, max, etc.) on some score for annotation regions that overlap.

It consists of two small scripts that rely on the excellent [bedtools](http://bedtools.readthedocs.io/en/latest/) package to do most of the work. I strongly recommend learning how to use bedtools if you are working with ChIP, ATAC, or DNase-Seq data.

<a name="input_files"></a>

## Input Files

To run the program, you will need to supply:

1) --input_bed, a bed file of genomic regions of interest (say, your ATAC-Seq peaks)   
2) --annotation_list, a tab-delimited file with nine fields of information for each annotation file.   
3) A set of annotation files,  each file needs a row of information in --annotation_list

### Input Bed

This is just a simple, three column file [chr,start,end] of your regions of interest. For an example, see example_data/example_peaks.bed 

```
chr1	100	200
chr2	1000	1100
chr9	10000	10010
```


### Annotation List

This is a tab-delimited file where every row provides information on one of your annotations. The nine columns needed are:

[Note: You do not need a header in the --annotation_list file, but I strongly recommend it. Please put a "#" at the start to comment it out, or the program will treat it as data.]  

* File - Path to the annotation's bed file. 
* Type - Operation to carry out: bed_overlap, bed_closest, or bed_score    
* Annotation_SizeInBP - This is the size of the annotation in bp. You can enter "-1" if you want BedOverlaps_DB to calculate it. The only advantage to providing it is that it can save you a small amount of time if the file is unusually large, or if you want to change the size for some other reason.
* Annotation_NumPeaks  - This is the number of peaks in the annotation. You can set this to "-1" to have the program calculate it.   
* Alias - An alternative, short name for the bed file
* ScoreColumn - This is required if you want to run a "bed_score" operation. For example, if a bed file has evolutionary conservation in the 4th column, you would enter a 4 here. You can leave it blank if you are not doing "bed_score".   
* Bedtools function - Function to apply to score column. Like mean, max, etc. You can leave it blank if you are not doing "bed_score".    
* Threshold - Threshold to apply to score column. You can leave it blank if you are not doing "bed_score".    
* EffectiveGenomeSize - This is the effective genome size to use for the annotation when doing enrichment. Enrichment is calculated as (SizeOfOverlapBetweenInputBedAndAnnotationInBP / SizeOfPeaksInBP) / (Annotation_SizeInBP / EffectiveGenomeSize )

<a name="output_files"></a>

## Output Files

BedOverlapsDB produces two main output files:

**--out_overlaps** - This is a large matrix with rows for each of your peaks and columns for the annotations. Annotations with "bed_overlap" will display the overlap in bp for that peak in the cell. Annotations with "bed_closest" will have a column with the nearest feature, and a column with distance to it in bp. Annotations with "bed_score" will display the summary stat specified for that score for regions in the annotation overlapping your peak.

If you ran the example, you should get an output matrix that looks something like this:

| chr  | start | end   | OverlapTest_OverlapBP | ScoreTest_meanscore | Closest_ClosestFeatureTest | DistToClosest_ClosestFeatureTest |
|------|-------|-------|-----------------------|---------------------|----------------------------|----------------------------------|
| chr1 | 100   | 200   | 60                    | 30                  | Feature_1                  | 51                               |
| chr2 | 1000  | 1100  | 0                     | 0                   | Feature_2                  | 7901                             |
| chr9 | 10000 | 10010 | 0                     | 0                   | NoFeature                  | 0                                |

* chr,start,end - This is the locus from the original input_bed.   
* OverlapTest_OverlapBP - The number of bp that overlap the input locus from the OverlapTest annotation.   
* ScoreTest_meanscore - The mean value of the value in ScoreTest for the annotation peaks overlapping the input locus.     
* Closest_ClosestFeatureTest - The nearest feature for the first loci, and "NoFeature" for the last region, since the annotation had nothing on Chr9.   
* DistToClosest_ClosestFeatureTest - Thus reports the distance in bp from the input locus to the feature.   

**--out_summary** This reports summary stats and hypothesis tests for each of the bed_overlap operations. It reports the following:

* Annotation - Name of the annotation   
* Annotation_BP - Size of the annotation in BP   
* Annotation_Peaks - Number of peaks in the annotation   
* EffectiveGenomeSize_BP - EffectiveGenomeSize_BP used for the Annotation's enrichment test   
* TotalOverlap_BP - Size of total overlap in basepairs   
* TotalOverlap_Peaks - Number of peaks that overlapped the annotation   
* UserSet_BP - Total size of input_bed in bp   
* UserSet_Peaks - Total number of peaks in input_bed   
* PeaksInWorld - EffectiveGenomeSize/(AnnotationSize_BP/AnnotationSize_Peaks)   
* Enrichment_Peaks - (TotalOverlap_Peaks/UserSet_Peaks)/(Annotation_Peaks/PeaksInWorld)   
* Enrichment_BP - (TotalOverlap_BP/UserSet_BP)/(Annotation_BP/EffectiveGenomeSize_BP)    

The hypergeometric tests are calculated in the R script (src/OverlapStats.R) as follows:

```r
 Hypergeometric_Peaks_PVal_OneSided = 1-phyper(q=TotalOverlap_Peaks,      # black balls picked
                                               m=Annotation_Peaks,        # black balls in urn
                                               n=PeaksInWorld-Annotation_Peaks, # white balls in urn
                                               k=UserSet_Peaks), # Total balls we picked   

* Hypergeometric_BP_PVal_OneSided = 1-phyper(q=TotalOverlap_BP,      # black balls picked
                                             m=Annotation_BP,        # black balls in urn
                                             n=EffectiveGenomeSize_BP-Annotation_BP, # white balls in urn
                                             k=UserSet_BP)) # Total balls we picked   
```


   

## Common Tasks

**Note: I need to check over these and make sure they are 100% correct.** For now, I have put short descriptions of how the program calculates these.

<a name="bed_overlap"></a>

## bed_overlap: Find Overlap and Enrichment

Process:

* Cuts the annotation file down to the [chr, start, end] columns   
* Sorts it   
* Merges the overlapping regions within bedAnnotFile   
* Intersects this data with (bedIn), and reports the number of bp overlap   
* Cuts down file to the original bedIn columns, and bp overlap col[chr,start,end,bp_overlap]   
* Sums the overlap over the row of bedIn, reduces this to unique rows of bedIn   

<a name="bed_closest"></a>

## bed_closest: Find Nearest Features

For each row in input, list the closest feature in the annotation file and the distance to it.
    
Process:

* Cut the annotation file (bedAnnotFile) to the [chr, start, end] columns   
* Sort it   
* Run "bedtools closest"   
* Cut down the results file to [chr,start,end,feature,dist_to_feature]   
* Run sed to fix instances where there is no closest feature (should be rare).   


<a name="bed_score"></a>

## bed_score: Produce a Summary Statistic on Score Column 

Perform left-outer-join of bedAnnotFile on bedIn, sums bp_overlap for each row of bedIn 
    
Process:

* Sort the annotation file (bedAnnotFile)    
* Run bedtools map - this calls the function of interest on the appropriate column.    
