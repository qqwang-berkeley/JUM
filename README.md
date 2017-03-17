# JUM
REQUIREMENTS: 
     ● Perl (5+)
     ● Samtools
     ● Bedtools (2.26.0+)
     ● R
     ● STAR

INSTALLATION: Download the JUM package to a local folder of choice (for example, /home/eagle/) and unpack.

MANUAL: 
1. Perform the first round STAR mapping of RNA-seq reads as instructed in the STAR manual.

2. Create a working folder, (for example, /home/eagle/JUMwork) and copy STAR output with suffix “Aligned.out.sam” and “SJ.out.tab” from all samples to the folder. For example, suppose the user has IndexA (control), IndexB (control), IndexC (treatment) and IndexD (treatment), then the user should copy the following STAR output files to /home/eagle/JUMwork.

IndexAAligned.out.sam
IndexASJ.out.tab

IndexBAligned.out.sam
IndexBSJ.out.tab

IndexCAligned.out.sam
IndexCSJ.out.tab

IndexDAligned.out.sam
IndexDSJ.out.tab

3. Run JUM_1.sh in the JUM package:
bash /home/eagle/JUM/JUM_1.sh

This step will create a file called: combined_SJ_out_tab_unannotated_for_2nd_pass_genome_generation.txt

4. Delete all Aligned.out.sam and SJ.out.tab files from this folder (to avoid confusion with second round STAR mapping results).

5. Perform second round of STAR mapping of RNA-seq reads. First make a STAR genome index for the second round of STAR mapping by supplying combined_SJ_out_tab_unannotated_for_2nd_pass_genome_generation.txt as the --sjdbFileChrStartEnd parameter. Next, run STAR for all your RNA-seq samples using the resulted STAR genome index.

6. Copy the resulted second round STAR output with suffix “Aligned.out.sam” and “SJ.out.tab” to the working folder (i.e. /home/eagle/JUMwork).

7. Use samtools and bedtools to transform all Aligned.out.sam files (use the exact naming as below):
samtools view -bS IndexAAligned.out.sam > IndexAAligned.out.bam
samtools sort IndexAAligned.out.bam > IndexAAligned.out_sorted
bedtools genomecov -ibam IndexAAligned.out_sorted.bam -bga -split > IndexAAligned.out_coverage.bed
     Do these steps for all your samples.
    
8. Run JUM_2.sh (this step may take up to several hours depending on number of samples to run):
bash /home/eagle/JUM/JUM_2.sh #directory #read_threshold_1 #read_threshold_2 #read_length
(JUM_2.sh needs four input parameters: 
 parameter a. #directory: path of the downloaded JUM package.
 parameter b. #read_number_threshold_for_junction: JUM will only consider junctions that have more than this # of unique reads mapped to it in all samples as valid junctions;
 parameter c. #read_number_threshod_for_exon_intron_boundary: JUM will only consider retained introns that have more than this # of unique reads mapped to the upstream exon-intron and downstream intron-exon boundaries as valid potential retained introns;
 parameter d. #read_length: the length of the RNA-seq reads)
       Example: bash /home/eagle/JUM/JUM_2.sh /home/eagle/JUM 5 2 100
       JUM_2 outputs results into a new folder called JUM_diff/
9. Enter the folder JUM_diff/. Run the R script in the JUM package in this folder, with a user-provided experiment design file (txt format) for differential AS analysis. 
Rscript /home/eagle/JUM/R_script_JUM.R experiment_design.txt > outputFile.Rout 2> errorFile.Rout
       An example experiment_design.txt file is included in the package for the example below: suppose the user has IndexA,    
       IndexB, IndexC and IndexD samples. IndexA and IndexB are drug treated biological replicates and IndexC and IndexD are 
       control biological replicates. JUM_2.sh then outputs the following files IndexA_2nd_combined_count.txt, 
       IndexB_2nd_combined_count.txt, IndexC_2nd_combined_count.txt, and IndexD_2nd_combined_count.txt. Note, it is 
       important to keep the samples names and condition names in the same alphabetic order in the experiment design file.
    
       R_script_JUM.R will output a file called AS_differential.txt
10. Run JUM_3.sh:
bash /home/eagle/JUM/JUM_3.sh #directory #pvalue|adjusted_pvalue #threshold #number_of_samples #number_of_control_samples|treated_samples>
(JUM_3.sh needs five input parameters:
parameter a. #directory: path of the downloaded JUM package
parameter b. #pvalue|adjusted_pvalue: statistical standard to use as cut-off. Choose between pvalue or padj (multi-test adjusted pvalue).
parameter c. #threshold: cutoff for pvalue or padj (I.e. 0.05, 0.01, …)
parameter d. #number_of_samples: total number of RNA-seq samples
Parameter e. #number of control samples or the number of treated samples, whichever is smaller
       For example, suppose the user has a total of 4 samples (two controls and two treated), and wish to use pvalue 0.05 as cut off, then run JUM_3.sh as follows:
bash /home/eagle/JUM/JUM_3.sh /home/eagle/JUM pvalue 0.05 4 2

JUM_3.sh outputs a folder: FINAL_JUM_OUTPUT with the following files:
AS_differential_JUM_output_A3SS_events_pvalue_0.05.txt
AS_differential_JUM_output_A5SS_events_pvalue_0.05.txt
AS_differential_JUM_output_cassette_exon_events_pvalue_0.05.txt
AS_differential_JUM_output_intron_retention_pvalue_0.05.txt
AS_differential_JUM_output_mixed_events_pvalue_0.05.txt
AS_differential_JUM_output_MXE_events_pvalue_0.05.txt
Valid_total_A3SS_event.txt
Valid_total_A5SS_event.txt
Valid_total_cassette_exon_event.txt
total_intron_retention_event.txt
total_mixed_event.txt
Valid_total_MXE_listt.txt
cassette_exon_coordinate.bed
MXE_coordinate.bed
