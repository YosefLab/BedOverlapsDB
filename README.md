This will be replaced by a more detailed instruction manual in the future. For now, if you wish to use this program please download/clone the repository and run this command to determine if it works on your system:

python OverlapAnnotationsWithBed.py --annotation_list hg19_annotation_quicktest.tab --input_bed example_data/example_peaks.bed --genome hg19 --additional_text_labels SampleName --additional_text example_data

In general, to use the tool, please supply the following:
--annotation_list 
    This is a tab delimited file which should contain a header, and then a row for each annotation file you wish to use. The columns are as follows:
    file - full path to the bed file you are using as an annotation.	
    type - One of "bed_score", "bed_closest" or "bed_overlap". This tells the program what to do with the data in the file.
    annot_size - Optional. You can tell the program how large your annotation file is (in bp). If not, it will run an additional command to calculate the size.
    ValueColumn - If you are running "bed_score", this tells the program which column of the bed file as the score value. For example, if you wanted to calculate average conservation over a peak, you need to tell the program which column in the conservation file as the actual conservation value.
    Function - If you are running a "bed_score" file, you can tell it whether to calculate the mean, max, min, on your interval, for example.
    Thresh - 
    EffectiveGenome - Size of effective genome for this annotation. (Will provide a more detailed explanation later.)
    
--input_bed 
    Bed file of peaks that you want to run . See "example_data/example_peaks.bed" for an example.
--genome 
    hg19 or mm9. This lets the program know how big the effective genome size should be.
--additional_text_labels and --additional_text
    The program will output a file called "summary_overlaps.tab" at the end, which lists stats for each annotation on your sample. It's often the case that you will want to run many samples. To help identify them in the data, you can add a column to summary_overlaps.tab. The header of this column will be the value given to "--additional_text_labels", and the value will be the value given to "--additonal_text". This makes it easier to cat all the "summary_overlaps.tab" files at the end and know which sample is which.
    
Output:
  summary_overlaps.tab - Lists the stats for each annotation on your data.
  overlaps.tab - Lists each interval in your original bed file, and the corresponding data for each value.

Some notes:
*This should work on s121-s123. I think we are currently missing the optparse package on s124, so R crashes on s124.
*The program currently makes an "out" folder which contains a tmp directory. I will update the code to delete this in the future.
