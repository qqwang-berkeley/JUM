#!/bin/bash
set -e
# $1=directory for JUM scripts; $2=threshold for junction reads; $3=file_number; $4=read # threshold for intron_exon boundary reads; $5=read length;
################################# 
folder="$1";
threshold="$2";
file_num="$3";
IR_threshold="$4";
read_length="$5";
##################################

cat *more_than_"$threshold"_read_in_at_least_*_samples_formatted.txt | sort -u > UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt;

awk '{print $1 "\t" $3 "\t" $4 "\t" $2 "\t" $5}' UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt > UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt;

awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $4-$3 "\t" $5}' UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt > UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_with_junction_length.txt;

less UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt | cut -f5 | sort > UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_junction_list.txt;

less UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt | cut -f1,2,3 | sort | uniq -d > more_than_"$threshold"_union_junc_coor_5_prime_ss_list_with_alternative_3_ss.txt;

perl $folder/profile_AS_events_sharing_5_prime_ss.pl UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt more_than_"$threshold"_union_junc_coor_5_prime_ss_list_with_alternative_3_ss.txt > more_than_"$threshold"_profiled_5_ss_and_corresponding_alternative_3_ss_junction_list.txt;

less UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt | cut -f1,2,4 | sort | uniq -d > more_than_"$threshold"_union_junc_coor_3_prime_ss_list_with_alternative_5_ss.txt;

perl $folder/profile_AS_events_sharing_3_prime_ss.pl UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt more_than_"$threshold"_union_junc_coor_3_prime_ss_list_with_alternative_5_ss.txt > more_than_"$threshold"_profiled_3_ss_and_corresponding_alternative_5_ss_junction_list.txt;

cat more_than_"$threshold"_profiled_5_ss_and_corresponding_alternative_3_ss_junction_list.txt more_than_"$threshold"_profiled_3_ss_and_corresponding_alternative_5_ss_junction_list.txt > more_than_"$threshold"_profiled_total_AS_event_junction.txt;

perl $folder/first_processing_generate_JUM_reference_for_profiled_total_alternative_splicing_event_junction.pl more_than_"$threshold"_profiled_total_AS_event_junction.txt more_than_"$threshold"_profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt;

perl $folder/second_processing_generate_JUM_reference_for_profiled_total_alternative_splicing_event_junction.pl more_than_"$threshold"_profiled_total_AS_event_junction.txt more_than_"$threshold"_profiled_total_AS_event_junction_second_processing_for_JUM_reference_building.txt;

perl $folder/generate_JUM_format_annotation_file_for_profiled_alternative_splicing_event_junction.pl more_than_"$threshold"_profiled_total_AS_event_junction_second_processing_for_JUM_reference_building.txt more_than_"$threshold"_profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt more_than_"$threshold"_profiled_total_AS_event_junction_JUM_annotation.gff;


for file in *junction_counts.txt;
do
    iden=${file%_junction_counts.txt}"_"fn_count.txt;
    perl $folder/making_junction_count_file_matching_annotation_for_JUM.pl more_than_"$threshold"_profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt $file $iden;
done


#Now for intron retention analyses
echo Preparing for intron retention AS events analyses...

for file in *Aligned.out.sam;
do
    name=${file%Aligned.out.sam}"_"temp1;
    header=${file%Aligned.out.sam}"_"sam_header;
    juncSam=${file%Aligned.out.sam}"_"Aligned.out.spanning_junction_reads.sam;
    juncBam=${file%Aligned.out.sam}"_"Aligned.out.spanning_junction_reads.bam;
    juncBed=${file%Aligned.out.sam}"_"Aligned.out.spanning_junction_reads.bed;
    overhang1=profiled"_"${file%Aligned.out.sam}"_"Aligned.out.spanning_junction_reads_junction_overhang_mapped_num.txt;
    overhang2=${file%Aligned.out.sam}"_"junction_counts_more_than_"$threshold"_in_all_samples_with_both_overhangs.txt;
    trash=${file%Aligned.out.sam}"_"trash_out;
    title=${file%Aligned.out.sam};
    juncount=${file%Aligned.out.sam}"_"junction_counts.txt;
    juncountAll=${file%Aligned.out.sam}"_"junction_counts_more_than_"$threshold"_in_all_samples.txt;
    samtools view -S $file | awk '($6 ~ /N/)' > $name;
    samtools view -H $file > $header;
    cat $header $name > $juncSam;
    samtools view -bS $juncSam > $juncBam;
    bedtools bamtobed -bed12 -i $juncBam > $juncBed;
    perl $folder/profile_overhang_from_STAR_output_step_1.pl $juncBed $overhang1;
    awk 'FNR==NR {arr[$0];next} ($5 in arr)' UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_junction_list.txt $juncount > $juncountAll;
    perl $folder/profile_overhang_from_STAR_output_step_2.pl $juncountAll $overhang1 $overhang2 $trash;
        if [ -s $trash ]
        then 
           echo "Error:$trash file should be empty!";
        else
           echo "Processing junction overhangs for Sample $title... Success."; 
        fi
done

perl $folder/extract_intron_retention_event_coordinate_from_samples_conditional_union.pl UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_with_junction_length.txt *junction_counts_more_than_"$threshold"_in_all_samples_with_both_overhangs.txt UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_with_junction_length_and_overhang_union_from_all_samples.txt;


perl $folder/extract_intron_retention_event_coordinate_conditional_union_splicing_junction_gff.pl UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_with_junction_length_and_overhang_union_from_all_samples.txt output_long_intron.gff output_short_intron.gff $read_length;


for sortfile in *Aligned.out_sorted.bam;
do
    coveragefile=${sortfile%_sorted.bam}"_"coverage.bed;
    intersect_bam_long_intron=${sortfile%_sorted.bam}"_"intersect_long_intron.bam;
    intersect_bam_short_intron=${sortfile%_sorted.bam}"_"intersect_short_intron.bam;
    intersectBed -abam $sortfile -b output_long_intron.gff -wb -f 1 > $intersect_bam_long_intron &
    intersectBed -abam $sortfile -b output_short_intron.gff -wb -f 1 > $intersect_bam_short_intron &
    genomeCoverageBed -ibam $sortfile -bga -split > $coveragefile &
done

wait

for intersect_long_bam in *intersect_long_intron.bam;
do
    intersect_long_sam=${intersect_long_bam%_intron.bam}"_"intron.sam
    intersect_long_sam_coor=${intersect_long_bam%_intron.bam}"_"intron_chr_coor.txt
    output_long_intron=${intersect_long_bam%_long_intron.bam}"_"long_intron_count_table
    samtools view $intersect_long_bam > $intersect_long_sam;
    less $intersect_long_sam | cut -f3,4 > $intersect_long_sam_coor;
    perl $folder/count_reads_for_intron_retention_sam.pl output_long_intron.gff $intersect_long_sam_coor $read_length > $output_long_intron &
done

for intersect_short_bam in *intersect_short_intron.bam;
do
    intersect_short_sam=${intersect_short_bam%_intron.bam}"_"intron.sam
    intersect_short_sam_coor=${intersect_short_bam%_intron.bam}"_"intron_chr_coor.txt
    output_short_intron=${intersect_short_bam%_short_intron.bam}"_"short_intron_count_table
    samtools view $intersect_short_bam > $intersect_short_sam;
    less $intersect_short_sam | cut -f3,4 > $intersect_short_sam_coor;
    perl $folder/count_reads_for_intron_retention_sam.pl output_short_intron.gff $intersect_short_sam_coor $read_length > $output_short_intron &
done

wait

#echo Preparation of counting for intron retention analyses complete!

#echo Starting intron retention analyses...


for file in *_intersect_long_intron_count_table;
do
    file_short=${file%_long_intron_count_table}"_"short_intron_count_table;
    name=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_more_than_"$threshold"_in_all_samples.txt;
    name_long_left_count=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_more_than_"$threshold"_in_all_samples_long_intron_retention_left_count.txt;
    name_long_left_count_list=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_more_than_"$threshold"_in_all_samples_long_intron_retention_left_count_list.txt;
    name_long_right_count=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_more_than_"$threshold"_in_all_samples_long_intron_retention_right_count.txt;
    name_long_right_count_list=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_more_than_"$threshold"_in_all_samples_long_intron_retention_right_count_list.txt;
    name_short_count=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_more_than_"$threshold"_in_all_samples_short_intron_retention.txt;
    name_short_count_list=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_more_than_"$threshold"_in_all_samples_short_intron_retention_list.txt;
    name_combined_list=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_combined_intron_retention_event_list.txt;
    name_combined=${file%Aligned.out_intersect_long_intron_count_table}"_"junction_counts_combined_intron_retention_event.txt;
    perl $folder/preparing_intron_retention_count_table_long_intron.pl $file $name $name_long_left_count $name_long_left_count_list $name_long_right_count $name_long_right_count_list $IR_threshold;
    perl $folder/preparing_intron_retention_count_table_short_intron.pl $file_short $name $name_short_count $name_short_count_list $IR_threshold;
    cat $name_long_left_count_list $name_long_right_count_list $name_short_count_list > $name_combined_list;
    cat $name_long_left_count $name_long_right_count $name_short_count > $name_combined;
done

perl $folder/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt All_junction_counts_intron_retention_in_all_samples_list.txt $file_num;
sort All_junction_counts_intron_retention_in_all_samples_list.txt > All_junction_counts_intron_retention_in_all_samples_sorted_list.txt;
perl $folder/preparing_intron_retention_event_for_JUM_combine.pl All_junction_counts_intron_retention_in_all_samples_sorted_list.txt output_long_intron.gff output_short_intron.gff UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt JUM_intron_retention_splicing_event.txt;

perl $folder/first_processing_generate_JUM_reference_for_profiled_total_alternative_splicing_event_junction.pl JUM_intron_retention_splicing_event.txt output_intron_retention_splicing_event_first_processing_for_JUM_reference_building.txt; 

perl $folder/second_processing_generate_JUM_reference_for_profiled_total_alternative_splicing_event_junction.pl JUM_intron_retention_splicing_event.txt output_intron_retention_splicing_event_second_processing_for_JUM_reference_building.txt;

perl $folder/generate_JUM_format_annotation_file_for_profiled_alternative_splicing_event_junction.pl output_intron_retention_splicing_event_second_processing_for_JUM_reference_building.txt output_intron_retention_splicing_event_first_processing_for_JUM_reference_building.txt output_intron_retention_splicing_event_JUM_annotation.gff;

cat more_than_"$threshold"_profiled_total_AS_event_junction_JUM_annotation.gff output_intron_retention_splicing_event_JUM_annotation.gff > combined_AS_JUM.gff;

for files in *junction_counts_combined_intron_retention_event.txt
do
   refom=${files%_combined_intron_retention_event.txt}"_"intron_retention_counts_exist_in_all_samples.txt;
   ori_no_AS=${files%_junction_counts_combined_intron_retention_event.txt}"_"fn_count.txt;
   combined_all_AS=${files%_junction_counts_combined_intron_retention_event.txt}"_"combined_count.txt;
   perl $folder/making_junction_count_file_matching_annotation_for_JUM_in_intron_retention.pl output_intron_retention_splicing_event_first_processing_for_JUM_reference_building.txt $files $refom;
   #awk -F ":" 'FNR==NR {arr[$0];next} ($1 in arr)' All_junction_counts_intron_retention_in_all_samples_list.txt $files | sort > $refom;
   cat $ori_no_AS $refom > $combined_all_AS;
done

echo ready to performing differential AS analyses...


mkdir JUM_diff
cp *combined_count.txt JUM_diff/
cp combined_AS_JUM.gff JUM_diff/
cp *Aligned.out_coverage.bed JUM_diff/
cp *profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt JUM_diff/
cp UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt JUM_diff/


