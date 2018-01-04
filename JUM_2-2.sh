#!/bin/bash
set -e
# $1=directory for JUM scripts; $2=threshold for junction reads; $3=folder/condition_name; $4=file_number (taking junctions with more than threshold reads in at least this number of samples under one specific condition);
################################# 
folder="$1";
threshold="$2";
folder_name="$3";
file_num="$4";
##################################

#Step 1
for file in *SJ.out.tab_strand_symbol_scaled;
do
     na1=${file%SJ.out.tab_strand_symbol_scaled}"_"idxstats_junction.txt;
     na2=${file%SJ.out.tab_strand_symbol_scaled}"_"junction_counts.txt;
     perl $folder/prepare_count_file_after_STAR_2_pass_mapping.pl $file UNION_junc_coor_with_junction_ID.txt $na1;
     awk '{print $1 "\t" $4 "\t" $2 "\t" $3 "\t" $10 "\t" $7 "\t" $9}' $na1 > $na2;
done

perl $folder/Identify_junctions_exist_in_certain_number_of_input_files.pl UNION_junc_coor_with_junction_ID.txt *_junction_counts.txt "$folder_name"_Junction_list_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt $file_num $threshold
awk 'FNR==NR {arr[$0];next} ($5 in arr)' "$folder_name"_Junction_list_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt UNION_junc_coor_with_junction_ID.txt > "$folder_name"_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt;
awk '{print $1 "\t" $4 "\t" $2 "\t" $3 "\t" $5}' "$folder_name"_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt > "$folder_name"_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted.txt;
