#!/bin/bash
#set -e
# $1=directory for JUM scripts; $2=pvalue or pvalue-adjusted; $3=pvalue or padj threshold; $4=file numbers; $5=control_file_numbers; 
#############

folder="$1";
pvalue_padj="$2";
cutoff="$3";
file_num="$4";
con_file_num="$5";
#############

modified_output=AS_differential_JUM_output.txt;
perl $folder/process_raw_AS_differential_output.pl AS_differential.txt $modified_output $file_num;

output_file=AS_differential_"$pvalue_padj"_"$cutoff".txt;
if [ $pvalue_padj == "pvalue" ];
then
     awk -v cut_off="$cutoff" '$7 <= cut_off' AS_differential.txt > $output_file;
else
     awk -v cut_off="$cutoff" '$8 <= cut_off' AS_differential.txt > $output_file;
fi

##################significant intron retention number#########################
less $output_file | grep "left" > left_span_long_intron_$output_file;
less $output_file | grep "right" > right_span_long_intron_$output_file;
less $output_file | grep "short" | cut -f2 | sort -u > AS_differential_"$pvalue_padj"_"$cutoff"_valid_short_intron_retention_list.txt; 
awk '($0 !~ /left/) && ($0 !~ /right/) && ($0 !~ /short/)' $output_file | cut -f2 | sort -u  > AS_differential_"$pvalue_padj"_"$cutoff"_non_intron_retention_AS_structure_list.txt;
##############################################################################

##################total intron retention number###############################
less AS_differential.txt | grep "left" > left_span_long_intron_AS_differential_total.txt;
less AS_differential.txt | grep "right" > right_span_long_intron_AS_differential_total.txt
less AS_differential.txt | grep "short" | cut -f2 | sort -u > AS_differential_total_valid_short_intron_retention_list.txt;
################################################################################

out_list_pvalue=${output_file%.txt}"_"list.txt;
out_list=AS_differential_total_list.txt;
perl $folder/prepare_significant_long_intron_retention_event.pl left_span_long_intron_$output_file right_span_long_intron_$output_file temp_long_intron_retention_$output_file temp_long_intron_retention_$out_list_pvalue;
perl $folder/prepare_significant_long_intron_retention_event.pl left_span_long_intron_AS_differential_total.txt right_span_long_intron_AS_differential_total.txt temp_long_intron_retention_AS_differential_total.txt temp_long_intron_retention_$out_list;

coor_file_pvalue=temp_long_intron_retention_junction_coordinate_"$pvalue_padj"_"$cutoff".txt;
coor_file_total=temp_long_intron_retention_junction_coordinate_total.txt
awk 'FNR==NR {arr[$0];next} ($5 in arr)' temp_long_intron_retention_$out_list_pvalue UNION_junc_coor_with_junction_ID* > $coor_file_pvalue;
awk 'FNR==NR {arr[$0];next} ($5 in arr)' temp_long_intron_retention_$out_list UNION_junc_coor_with_junction_ID* > $coor_file_total;

for files in *Aligned.out_coverage.bed
do
   overlap_pvalue=${files%Aligned.out_coverage.bed}"_"coverage_temp_long_intron_overlap_"$pvalue_padj"_"$cutoff".txt;
   read_num_pvalue=${files%Aligned.out_coverage.bed}"_"temp_long_intron_retention_junction_coordinate_with_read_num_"$pvalue_padj"_"$cutoff".txt;
   screen_pvalue=${files%Aligned.out_coverage.bed}"_"long_intron_retention_screening_"$pvalue_padj"_"$cutoff".txt;
   valid_pvalue=${files%Aligned.out_coverage.bed}"_"valid_long_intron_retention_event_list_"$pvalue_padj"_"$cutoff".txt;
   
   overlap=${files%Aligned.out_coverage.bed}"_"coverage_temp_long_intron_overlap_total.txt;
   read_num=${files%Aligned.out_coverage.bed}"_"temp_long_intron_retention_junction_coordinate_with_read_num_total.txt;
   screen=${files%Aligned.out_coverage.bed}"_"long_intron_retention_screening_total.txt;
   valid=${files%Aligned.out_coverage.bed}"_"valid_long_intron_retention_event_list_total.txt;
 
   intersectBed -wa -a $files -b $coor_file_pvalue > $overlap_pvalue;
   perl $folder/count_intron_read_long_intron_retention_step2.pl $overlap_pvalue $coor_file_pvalue $read_num_pvalue;
   perl $folder/determining_rightful_long_intron_retention_event.pl $read_num_pvalue $screen_pvalue;
   awk -v cut1=0 '$2==cut1, $3==cut1' $screen_pvalue | cut -f1 | sort -u > $valid_pvalue;

   intersectBed -wa -a $files -b $coor_file_total > $overlap;
   perl $folder/count_intron_read_long_intron_retention_step2.pl $overlap $coor_file_total $read_num;
   perl $folder/determining_rightful_long_intron_retention_event.pl $read_num $screen;
   awk -v cut1=0 '$2==cut1, $3==cut1' $screen | cut -f1 | sort -u > $valid;
done

cat *valid_long_intron_retention_event_list_"$pvalue_padj"_"$cutoff".txt | sort | uniq -cd | awk -v limit="$con_file_num" '$1 >= limit{print $2}' > AS_differential_"$pvalue_padj"_"$cutoff"_long_intron_retention_passing_filter_combined_list.txt;
cat *valid_long_intron_retention_event_list_total.txt | sort | uniq -cd | awk -v limit="$con_file_num" '$1 >= limit{print $2}' > AS_differential_total_long_intron_retention_passing_filter_combined_list.txt;

awk 'FNR==NR {arr[$0];next} ($1 in arr)' AS_differential_"$pvalue_padj"_"$cutoff"_long_intron_retention_passing_filter_combined_list.txt temp_long_intron_retention_$out_list_pvalue > AS_differential_"$pvalue_padj"_"$cutoff"_valid_long_intron_retention_list.txt;
awk 'FNR==NR {arr[$0];next} ($1 in arr)' AS_differential_total_long_intron_retention_passing_filter_combined_list.txt temp_long_intron_retention_$out_list > AS_differential_total_valid_long_intron_retention_list.txt;

################# profiling splicing patterns ##########################
perl $folder/profiling_splicing_patterns_from_AS_events_1.pl *profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt UNION_junc_coor_with_junction_ID* reconstruct_splicing_pattern_input_1.txt;
perl $folder/profiling_splicing_patterns_from_AS_events_2.pl reconstruct_splicing_pattern_input_1.txt reconstruct_splicing_pattern_input_2.txt;
perl $folder/profiling_splicing_patterns_from_AS_events_3_updated.pl reconstruct_splicing_pattern_input_1.txt reconstruct_splicing_pattern_input_2.txt total_A5SS_event.txt total_A3SS_event.txt total_cassette_exon_event.txt total_mixed_event.txt;
################# profiling splicing patterns ##########################

#############screening for cassette exon patterns##############
perl $folder/screening_cassette_exon_events.pl total_cassette_exon_event.txt reconstruct_splicing_pattern_input_1.txt cassette_exon_coordinate.bed

for files in *Aligned.out_coverage.bed
do
   intersect=${files%Aligned.out_coverage.bed}"_"coverage_temp_total_cassette_intersect.txt;
   read_per_ID=${files%Aligned.out_coverage.bed}"_"temp_total_cassette_coordinate_with_read_num_per_ID.txt;
   indicator=${files%Aligned.out_coverage.bed}"_"temp_total_cassette_with_indicator.txt;
   valid_SE=${files%Aligned.out_coverage.bed}"_"valid_total_cassette_exon_event_ID.txt;
 
   intersectBed -wa -a $files -b cassette_exon_coordinate.bed > $intersect ;
   perl $folder/count_intron_read_long_intron_retention_step2.pl $intersect cassette_exon_coordinate.bed $read_per_ID;
   perl $folder/determining_rightful_total_cassette_exon_event.pl $read_per_ID $indicator;
   awk -v cut1=0 '$2==cut1, $3==cut1' $indicator | cut -f1 | sort -u > $valid_SE;
done

cat *valid_total_cassette_exon_event_ID.txt | sort | uniq -cd | awk -v limit=1 '$1 >= limit{print $2}' > Valid_total_cassette_exon_list.txt;

#############screening for cassette exon patterns#############
awk 'FNR==NR {arr[$0];next} !($5 in arr)' Valid_total_cassette_exon_list.txt cassette_exon_coordinate.bed > non_valid_cassette_exon_list.txt;


#############screening for MXE #########################
perl $folder/screening_MXE_events.pl total_A5SS_event.txt total_A3SS_event.txt reconstruct_splicing_pattern_input_1.txt MXE_coordinate.bed

for files in *Aligned.out_coverage.bed
do
   intersect_MXE=${files%Aligned.out_coverage.bed}"_"coverage_temp_total_MXE_intersect.txt;
   read_per_ID_MXE=${files%Aligned.out_coverage.bed}"_"temp_total_MXE_coordinate_with_read_num_per_ID.txt;
   indicator_MXE=${files%Aligned.out_coverage.bed}"_"temp_total_MXE_with_indicator.txt;
   valid_SE_MXE=${files%Aligned.out_coverage.bed}"_"valid_total_MXE_event_ID.txt;
 
   intersectBed -wa -a $files -b MXE_coordinate.bed > $intersect_MXE ;
   perl $folder/count_intron_read_long_intron_retention_step2.pl $intersect_MXE MXE_coordinate.bed $read_per_ID_MXE;
   perl $folder/determining_rightful_total_cassette_exon_event.pl $read_per_ID_MXE $indicator_MXE;
   awk -v cut1=0 '$2==cut1, $3==cut1' $indicator_MXE | cut -f1 | sort -u > $valid_SE_MXE;
done

cat *valid_total_MXE_event_ID.txt | sort | uniq -cd | awk -v limit=1 '$1 >= limit{print $2}' > Valid_total_MXE_list.txt;
#############screening for MXE #########################

perl $folder/screening_A5SS_A3SS_events.pl non_valid_cassette_exon_list.txt Valid_total_MXE_list.txt total_A5SS_event.txt total_A3SS_event.txt Valid_total_A5SS_event.txt Valid_total_A3SS_event.txt




#######################Below deconvoluting AS patterns##########################

awk 'FNR==NR {arr[$0];next} ($1 in arr)' AS_differential_"$pvalue_padj"_"$cutoff"_non_intron_retention_AS_structure_list.txt Valid_total_A5SS_event.txt > AS_differential_"$pvalue_padj"_"$cutoff"_A5SS_event_AS_structure_list.txt;
perl $folder/profiling_splicing_patterns_from_AS_events_4.pl AS_differential_"$pvalue_padj"_"$cutoff"_A5SS_event_AS_structure_list.txt AS_differential_JUM_output.txt > AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff".txt
awk 'FNR==NR {arr[$0];next} ($1 in arr)' AS_differential_"$pvalue_padj"_"$cutoff"_non_intron_retention_AS_structure_list.txt Valid_total_A3SS_event.txt > AS_differential_"$pvalue_padj"_"$cutoff"_A3SS_event_AS_structure_list.txt;
perl $folder/profiling_splicing_patterns_from_AS_events_4.pl AS_differential_"$pvalue_padj"_"$cutoff"_A3SS_event_AS_structure_list.txt AS_differential_JUM_output.txt > AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff".txt


perl $folder/profiling_splicing_patterns_from_AS_events_5.pl Valid_total_cassette_exon_list.txt AS_differential_JUM_output.txt > AS_differential_JUM_output_total_cassette_exon.txt
if [ $pvalue_padj == "pvalue" ]
then
     awk -v cut_off="$cutoff" '$6 <= cut_off' AS_differential_JUM_output_total_cassette_exon.txt | cut -f1 | sort -u > AS_differential_"$pvalue_padj"_"$cutoff"_cassette_exon_AS_event_list.txt;
else
     awk -v cut_off="$cutoff" '$7 <= cut_off' AS_differential_JUM_output_total_cassette_exon.txt | cut -f1 | sort -u > AS_differential_"$pvalue_padj"_"$cutoff"_cassette_exon_AS_event_list.txt;
fi
awk 'FNR==NR {arr[$0];next} ($1 in arr)' AS_differential_"$pvalue_padj"_"$cutoff"_cassette_exon_AS_event_list.txt AS_differential_JUM_output_total_cassette_exon.txt > temp_cassette;
head -1 AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff".txt > title;
cat title temp_cassette > AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff".txt


perl $folder/profiling_splicing_patterns_from_AS_events_5.pl Valid_total_MXE_list.txt AS_differential_JUM_output.txt > AS_differential_JUM_output_total_MXE.txt
if [ $pvalue_padj == "pvalue" ]
then
     awk -v cut_off="$cutoff" '$6 <= cut_off' AS_differential_JUM_output_total_MXE.txt | cut -f1 | sort -u > AS_differential_"$pvalue_padj"_"$cutoff"_MXE_AS_event_list.txt;
else
     awk -v cut_off="$cutoff" '$7 <= cut_off' AS_differential_JUM_output_total_MXE.txt | cut -f1 | sort -u > AS_differential_"$pvalue_padj"_"$cutoff"_MXE_AS_event_list.txt;
fi
awk 'FNR==NR {arr[$0];next} ($1 in arr)' AS_differential_"$pvalue_padj"_"$cutoff"_MXE_AS_event_list.txt AS_differential_JUM_output_total_MXE.txt > temp_MXE;
cat title temp_MXE > AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff".txt



perl $folder/profiling_splicing_patterns_from_AS_events_5.pl total_mixed_event.txt AS_differential_JUM_output.txt > AS_differential_JUM_output_total_mixed.txt
if [ $pvalue_padj == "pvalue" ]
then
     awk -v cut_off="$cutoff" '$6 <= cut_off' AS_differential_JUM_output_total_mixed.txt | cut -f1 | sort -u > AS_differential_"$pvalue_padj"_"$cutoff"_mixed_AS_event_list.txt;
else
     awk -v cut_off="$cutoff" '$7 <= cut_off' AS_differential_JUM_output_total_mixed.txt | cut -f1 | sort -u > AS_differential_"$pvalue_padj"_"$cutoff"_mixed_AS_event_list.txt;
fi
awk 'FNR==NR {arr[$0];next} ($1 in arr)' AS_differential_"$pvalue_padj"_"$cutoff"_mixed_AS_event_list.txt AS_differential_JUM_output_total_mixed.txt > temp_mixed;
cat title temp_mixed > AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff".txt



cat AS_differential_total_valid_long_intron_retention_list.txt AS_differential_total_valid_short_intron_retention_list.txt > total_intron_retention_event.txt

perl $folder/profiling_splicing_patterns_from_AS_events_6.pl AS_differential_"$pvalue_padj"_"$cutoff"_valid_long_intron_retention_list.txt AS_differential_"$pvalue_padj"_"$cutoff"_valid_short_intron_retention_list.txt AS_differential_JUM_output.txt > temp_intron_retention;
#cat title temp_intron_retention > AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff".txt

less temp_intron_retention | cut -f2,3,10 > test_intron.txt;
perl $folder/screening_intron_retention_event_for_deltachange_consistency.pl test_intron.txt  > test_intron_delta_consistant_list.txt;
less test_intron_delta_consistant_list.txt | sort -u > temfile;
mv temfile test_intron_delta_consistant_list.txt;
perl $folder/screening_intron_retention_event_for_no_inner_junction.pl test_intron_delta_consistant_list.txt UNION_junc_coor_with_junction_ID* test_intron_delta_consistant_and_no_inner_junction_list.txt 
awk 'FNR==NR {arr[$0];next} ($2 in arr)' test_intron_delta_consistant_and_no_inner_junction_list.txt temp_intron_retention > temp_intron_retention_consistant;
cat title temp_intron_retention_consistant > AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff".txt

mkdir FINAL_JUM_OUTPUT_"$pvalue_padj"_"$cutoff";
mv AS_differential_JUM_output_*"$pvalue_padj"_"$cutoff".txt FINAL_JUM_OUTPUT_"$pvalue_padj"_"$cutoff";
mv Valid_* FINAL_JUM_OUTPUT_"$pvalue_padj"_"$cutoff";
mv total_mixed_event.txt FINAL_JUM_OUTPUT_"$pvalue_padj"_"$cutoff";
mv total_intron_retention_event.txt FINAL_JUM_OUTPUT_"$pvalue_padj"_"$cutoff";
mv *coordinate.bed FINAL_JUM_OUTPUT_"$pvalue_padj"_"$cutoff";

