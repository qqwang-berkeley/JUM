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
    "fileNum_threshold"
    "IRthreshold"
    "Readlength"
    "Thread"
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

        --fileNum_threshold)
            file_num=$2
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

        *)
            break
            ;;
    esac
done


start=$(date +%s)




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

echo "JUM_A_multi_2.sh core execution finished in $runtime seconds."
echo "Command is:"
echo "bash JUM_A_multi_2.sh --Folder $folder --JuncThreshold $threshold --fileNum_threshold $file_num --IRthreshold $IR_threshold --Readlength $read_length --Thread $thread_num


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
