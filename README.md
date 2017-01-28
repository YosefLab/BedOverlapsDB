BedOverlapsDB
===========

This pipeline take a set of genomic intervals and measures how enriched they are in a set of genomic annotations. An example use case would be taking a set of ATAC-Seq peaks, and measuring how enriched they are in a set of regulatory annotations in the mouse genome.

Jim created and maintains this pipeline. 
Please note that this documentation is a work in progress.

---

#### Table of Contents
[Important Notes](#important_notes)    
[Quickstart](#quickstart)    

#### Overview
[Introduction](#intro)
[Input Files](#input files)
[Output Files](#output_files)

##### Common Tasks   
[bed_overlap: Find Overlap and Enrichment](#bed_overlap)  
[bed_closest: Find Nearest Features](#bed_closest)  
[bed_score: Produce a Summary Statistic on Score Column ](#bed_score)  

---

<a name="important_notes"></a>

## Important Notes

The pipeline currently relies on Jim's directory of R packages to load some things we don't have installed lab-wide yet. That's fine as a temporary measure, but we should install the "tidyverse" on s121-s124 when there is time.

<a name="quickstart"></a>

## Quickstart

To make sure you can use this repository, please run the following test from inside the BedOverlapsDB folder:

```
python OverlapAnnotationsWithBed.py --annotation_list example_data/example_annotations.tab  \
--input_bed example_data/example_peaks.bed --additional_text_labels SampleName \
 --additional_text test_sample --out_overlaps my_test_overlap_matrix.tab --out_summary my_test_stat_summary.tab
```

That will run the example_peaks.bed against a set of small sample annotations included with the package, and carry out the three basic operations performed by the package (measure overlap, find nearest feature, calculate a stat based on a summary score). It's a good first test to make sure that you have access to all the packages and programs you need to run BedOverlapsDB.

<a name="Introduction"></a>

## Introduction

Many of us work with data from ATAC-Seq, DNase-Seq and other experimental techniques that generate reads from particular sections of the genome that we expect to have some functional purpose. ATAC-Seq, for example, generates reads from open chromatin regions, and these regions typically contain things like enhancers. Some of the questions you might have a researcher are:

* Are these regions enriched for widely used lists of regulatory regions, (say ORegAnno)
* Are these regions close to transcription start sites (TSS's)? Are they near any genes we're interested in?
* Do these regions have a high average score for some measure of conservation?

The BedOverlapsDB is a tool answering these questions. It answers the first questions by measuring enrichment for a set of features. It answers the second by returning the closest interval in a set of features. It answers the last question by returning a summary statistic on some score.

Really, it consists of two small scripts that rely on the excellent [bedtools](http://bedtools.readthedocs.io/en/latest/) package to do most of the work.




# Overlaps DB





In general, to use the tool, please supply the following:


1. --annotation_list 
    This is a tab delimited file which should contain a header, and then a row for each annotation file you wish to use. The columns are as follows:
    * file - full path to the bed file you are using as an annotation.	
    * type - One of "bed_score", "bed_closest" or "bed_overlap". This tells the program what to do with the data in the file.
    * annot_size - Optional. You can tell the program how large your annotation file is (in bp). If not, it will run an additional command to calculate the size.
    * ValueColumn - If you are running "bed_score", this tells the program which column of the bed file to use as the score value. For example, if you wanted to calculate average conservation over a peak, you need to tell the program which column in the conservation file as the actual conservation value.
    * Function - If you are running a "bed_score" file, you can tell it whether to calculate the mean, max, min, on your interval, for example.
    * Thresh - 
    * EffectiveGenome - Size of effective genome for this annotation. (Will provide a more detailed explanation later.)
    
2. --input_bed 
    Bed file of peaks that you want to run . See "example_data/example_peaks.bed" for an example.

3. --additional_text_labels and --additional_text
    The program will output a file called "summary_overlaps.tab" at the end, which lists stats for each annotation on your sample. It's often the case that you will want to run many samples. To help identify them in the data, you can add a column to summary_overlaps.tab. The header of this column will be the value given to "--additional_text_labels", and the value will be the value given to "--additonal_text". This makes it easier to cat all the "summary_overlaps.tab" files at the end and know which sample is which.
    
Output:

1. summary_overlaps.tab - Lists the stats for each annotation on your data.

2. overlaps.tab - Lists each interval in your original bed file, and the corresponding data for each value.

##Some notes:
* This should work on s121-s123. I think we are currently missing the optparse package on s124, so R crashes on s124.

