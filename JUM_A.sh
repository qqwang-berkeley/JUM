#!/bin/bash
set -e
# $1=directory for JUM scripts; $2=threshold for junction reads; $3=file_number; $4=read # threshold for intron_exon boundary reads; $5=read length; $6=thread number for sam/bam processing;
################################# 
#folder="$1";
#threshold="$2";
#file_num="$3";
#IR_threshold="$4";
#read_length="$5";
#thread_num="$6";
##################################

ARGUMENT_LIST=(
    "Folder"
    "JuncThreshold"
    "Condition1_fileNum_threshold"
    "Condition2_fileNum_threshold"
    "IRthreshold"
    "Readlength"
    "Thread"
    "Condition1SampleName"
    "Condition2SampleName"
)


# read arguments
opts=$(getopt \
    --longoptions "$(printf "%s:," "${ARGUMENT_LIST[@]}")" \
    --name "$(basename "$0")" \
    --options "" \
    -- "$@"
)

eval set --$opts

while [[ $# -gt 0 ]]; do
    case "$1" in
        --Folder)
            folder=$2
            shift 2
            ;;

        --JuncThreshold)
            threshold=$2
            shift 2
            ;;

        --Condition1_fileNum_threshold)
            file_num_1=$2
            shift 2
            ;;
       
	--Condition2_fileNum_threshold)
            file_num_2=$2
            shift 2
            ;;

        --IRthreshold)
	    IR_threshold=$2
            shift 2
	    ;;
      
        --Readlength)
            read_length=$2
            shift 2
	    ;;
        	    
        --Thread)
            thread_num=$2
            shift 2
            ;;

        --Condition1SampleName)
            condition_1=$2
            shift 2
            ;;

        --Condition2SampleName)
            condition_2=$2
            shift 2
            ;;

        *)
            break
            ;;
    esac
done


start=$(date +%s)

bash $folder/JUM_2-1.sh;

if [ -d "con1_1" ]; then
rm -r con1_1;
fi

mkdir con1_1;

if [ -d "con2_1" ]; then
rm -r con2_1;
fi

mkdir con2_1;

echo "Sample names are:";
IFS=',' read -r -a array1 <<< "$condition_1"

for element in "${array1[@]}"
do
    cp "$element"SJ.out.tab_strand_symbol_scaled con1_1/;
    echo "$element";
done

IFS=',' read -r -a array2 <<< "$condition_2"

for element in "${array2[@]}"
do
    cp "$element"SJ.out.tab_strand_symbol_scaled con2_1/;
    echo "$element";
done

cp UNION_junc_coor_with_junction_ID.txt con1_1/;
cp UNION_junc_coor_with_junction_ID.txt con2_1/;

cd con1_1/;
bash $folder/JUM_2-2.sh --Folder $folder --Threshold $threshold --Filenum $file_num_1  --Condition "condition1"
cd ..

cd con2_1/;
bash $folder/JUM_2-2.sh --Folder $folder --Threshold $threshold --Filenum $file_num_2  --Condition "condition2"
cd ..

cp con1_1/*junction_counts.txt .
cp con2_1/*junction_counts.txt .

cp con1_1/*formatted.txt .
cp con2_1/*formatted.txt .

file_num=0;
if [ "$file_num_1" -ge "$file_num_2" ]; then
    file_num="$file_num_1"
fi

if [ "$file_num_1" -lt "$file_num_2" ]; then
    file_num="$file_num_2"
fi

#echo "$file_num";

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


echo Preparing for intron retention AS events analyses...

process_sam_func() {
    file=$1;
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
    samtools view -@ $thread_num -S $file | awk '($6 ~ /N/)' > $name;
    samtools view -H $file > $header;
    cat $header $name > $juncSam;
    samtools view -@ $thread_num -bS $juncSam > $juncBam;
    bedtools bamtobed -bed12 -i $juncBam > $juncBed;
    perl $folder/profile_overhang_from_STAR_output_step_1.pl $juncBed $overhang1 $read_length;
    awk 'FNR==NR {arr[$0];next} ($5 in arr)' UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_junction_list.txt $juncount > $juncountAll;
    perl $folder/profile_overhang_from_STAR_output_step_2.pl $juncountAll $overhang1 $overhang2 $trash;
        if [ -s $trash ]
        then 
           echo "Error:$trash file should be empty!";
        else
           echo "Processing junction overhangs for Sample $title... Success."; 
        fi
}

for file in *Aligned.out.sam;
do
	process_sam_func $file &

done

wait

perl $folder/extract_intron_retention_event_coordinate_from_samples_conditional_union.pl UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_with_junction_length.txt *junction_counts_more_than_"$threshold"_in_all_samples_with_both_overhangs.txt UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_with_junction_length_and_overhang_union_from_all_samples.txt;


perl $folder/extract_intron_retention_event_coordinate_conditional_union_splicing_junction_gff.pl UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples_formatted_with_junction_length_and_overhang_union_from_all_samples.txt output_long_intron.gff output_short_intron.gff $read_length;

sort -k1,1 -k4,4n output_long_intron.gff > output_long_intron_sorted.gff;

awk -v readlen="$read_length" '$5-$4 >= readlen' output_short_intron.gff > output_short_intron1.gff;
awk -v readlen="$read_length" '$5-$4 < readlen' output_short_intron.gff > output_short_intron2.gff;
sort -k1,1 -k4,4n output_short_intron1.gff > output_short_intron1_sorted.gff;
sort -k1,1 -k4,4n output_short_intron2.gff > output_short_intron2_sorted.gff;

awk -v readlen="$read_length" '{ gsub("junc_id=","",$9); print $1 "\t" $4 "\t" $5-readlen+1 "\t" $9 "\t" "1" "\t" $7}' output_long_intron.gff | sort -k1,1 -k2,2n > output_long_intron_adjusted_range_sorted.bed;
awk -v readlen="$read_length" '{ gsub("junc_id=","",$9); print $1 "\t" $4 "\t" $5-readlen+1 "\t" $9 "\t" "1" "\t" $7}' output_short_intron1.gff | sort -k1,1 -k2,2n > output_short_intron_1_adjusted_range_sorted.bed;
awk -v readlen="$read_length" '{ gsub("junc_id=","",$9); print $1 "\t" $5-readlen "\t" $4+1 "\t" $9 "\t" "1" "\t" $7}' output_short_intron2.gff | sort -k1,1 -k2,2n > output_short_intron_2_adjusted_range_sorted.bed;


process_sortedbam () {
    sortfile=$1;
    bamfiletobed=${sortfile%_sorted.bam}.bed;
    bamsortedbed=${sortfile%_sorted.bam}.sorted.bed;
    intersect_bam_long_intron=${sortfile%_sorted.bam}"_"intersect_long_intron.txt;
    intersect_bam_short_intron1=${sortfile%_sorted.bam}"_"intersect_short_intron1.txt;
    intersect_bam_short_intron2=${sortfile%_sorted.bam}"_"intersect_short_intron2.txt;
    bedtools bamtobed -i $sortfile > $bamfiletobed;
    sort -k1,1 -k2,2n $bamfiletobed > $bamsortedbed;
    intersectBed -a $bamsortedbed -b output_long_intron_sorted.gff -sorted -wa -u -f 1 > $intersect_bam_long_intron 
    intersectBed -a $bamsortedbed -b output_short_intron1_sorted.gff -sorted -wa -u -f 1 > $intersect_bam_short_intron1 
    intersectBed -a $bamsortedbed -b output_short_intron2_sorted.gff -sorted -wa -u -F 1 > $intersect_bam_short_intron2 
}    

for sortfile in *Aligned.out_sorted.bam;
do
    process_sortedbam $sortfile &
done    

wait

coverage_bamfile () {
    bamfile=$1;
    coveragefile=${bamfile%_sorted.bam}"_"coverage.bed;
    genomeCoverageBed -ibam $bamfile -bga -split > $coveragefile;
}

for bamfile in *Aligned.out_sorted.bam;
do
    coverage_bamfile $bamfile &
done

wait

process_inter_long_intron() {
    intersect_long_intron=$1;
    intersect_long_intron_coor=${intersect_long_intron%.txt}"_"chr_coor.txt
    output_long_intron=${intersect_long_intron%.txt}"_"count_table
    awk '{print $1 "\t" $2+1 "\t" $2+2 "\t" $4 "\t" $5 "\t" $6}' $intersect_long_intron | sort -k1,1 -k2,2n > $intersect_long_intron_coor;
    intersectBed -a output_long_intron_adjusted_range_sorted.bed -b $intersect_long_intron_coor -wa -c -sorted | sort -V -k4,4 | awk '{print $4 "\t" $7}' > $output_long_intron
}

for intersect_long_intron in *intersect_long_intron.txt;
do
	process_inter_long_intron $intersect_long_intron &
done

wait

process_inter_short_intron() {
    intersect_short_intron1=$1;
    intersect_short_intron2=${intersect_short_intron1%_intron1.txt}"_"intron2.txt;
    intersect_short_intron2_edit=${intersect_short_intron1%_intron1.txt}"_"intron2_edited.txt;
    intersect_short_intron1_coor=${intersect_short_intron1%.txt}"_"chr_coor.txt;
    output_short_intron1=${intersect_short_intron1%.txt}"_"count_table;
    intersect_short_intron2_coor=${intersect_short_intron2%.txt}"_"chr_coor.txt;
    output_short_intron2=${intersect_short_intron2%.txt}"_"count_table;
    output_short_intron=${intersect_short_intron1%_intron1.txt}"_"intron_count_table;
    awk '{print $1 "\t" $2+1 "\t" $2+2 "\t" $4 "\t" $5 "\t" $6}' $intersect_short_intron1 | sort -k1,1 -k2,2n > $intersect_short_intron1_coor;
    intersectBed -a output_short_intron_1_adjusted_range_sorted.bed -b $intersect_short_intron1_coor -wa -c -sorted | sort -V -k4,4 | awk '{print $4 "\t" $7}' > $output_short_intron1;
    awk -v readlen="$read_length" '$3-$2 <= readlen' $intersect_short_intron2 > $intersect_short_intron2_edit;
    awk '{print $1 "\t" $2+1 "\t" $2+2 "\t" $4 "\t" $5 "\t" $6}' $intersect_short_intron2_edit | sort -k1,1 -k2,2n > $intersect_short_intron2_coor;
    intersectBed -a output_short_intron_2_adjusted_range_sorted.bed -b $intersect_short_intron2_coor -wa -c -sorted | sort -V -k4,4 | awk '{print $4 "\t" $7}' > $output_short_intron2;   
    cat $output_short_intron1 $output_short_intron2 > $output_short_intron;
}

for intersect_short_intron1 in *intersect_short_intron1.txt;
do
	process_inter_short_intron $intersect_short_intron1 &
done

wait


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

if [ -d "con1_2" ]; then
rm -r con1_2;
fi

mkdir con1_2;

if [ -d "con2_2" ]; then
rm -r con2_2;
fi

mkdir con2_2;

for element in "${array1[@]}"
do
    cp "$element"_junction_counts_combined_intron_retention_event_list.txt con1_2/;
    echo "$element";
done

for element in "${array2[@]}"
do
    cp "$element"_junction_counts_combined_intron_retention_event_list.txt con2_2/;
    echo "$element";
done

cd con1_2/;

perl $folder/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt condition1_junction_counts_intron_retention_in_all_samples_list.txt $file_num_1;

cd ..

cd con2_2/;
perl $folder/Identify_intron_retention_event_exist_in_all_samples.pl *junction_counts_combined_intron_retention_event_list.txt condition2_junction_counts_intron_retention_in_all_samples_list.txt $file_num_2;

cd ..

cat con1_2/condition1_junction_counts_intron_retention_in_all_samples_list.txt con2_2/condition2_junction_counts_intron_retention_in_all_samples_list.txt | sort -u > All_junction_counts_intron_retention_in_all_samples_sorted_list.txt;

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
   cat $ori_no_AS $refom > $combined_all_AS;
done

echo ready to performing differential AS analyses...

end=$(date +%s)

runtime=$((end-start))

echo "JUM_A.sh core execution finished in $runtime seconds."
echo "Command is:"
echo "bash JUM_A.sh --Folder $folder --JuncThreshold $threshold --Condition1_fileNum_threshold $file_num_1 --Condition2_fileNum_threshold $file_num_2 --IRthreshold $IR_threshold --Readlength $read_length --Thread $thread_num --Condition1SampleName $condition_1 --Condition2SampleName $condition_2"


if [ -d "temp_JUM_A_run" ]; then
rm -r temp_JUM_A_run;
fi

mkdir temp_JUM_A_run;

ls | grep -v temp_JUM_A_run | xargs mv -t temp_JUM_A_run;

if [ -d "JUM_diff" ]; then
rm -r JUM_diff;
fi

mkdir JUM_diff;

mv temp_JUM_A_run/*combined_count.txt JUM_diff/
mv temp_JUM_A_run/combined_AS_JUM.gff JUM_diff/
mv temp_JUM_A_run/*Aligned.out_coverage.bed JUM_diff/
mv temp_JUM_A_run/*profiled_total_AS_event_junction_first_processing_for_JUM_reference_building.txt JUM_diff/
mv temp_JUM_A_run/UNION_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$file_num"_samples.txt JUM_diff/

mv temp_JUM_A_run/*Aligned.out.sam .
mv temp_JUM_A_run/*Aligned.out_sorted.bam* .
mv temp_JUM_A_run/*formatted.txt .
mv temp_JUM_A_run/*SJ.out.tab .

echo "JUM_A.sh completed. Ready for next step! Feel free to delete the temporary file folder temp_JUM_A_run once satisfied."
