#!/bin/bash
set -e
# $1=directory for JUM scripts; $2=threshold for junction reads; $3=filenumber (taking junctions with more than threshold reads in at least this number of samples under one specific condition); $4=condition
################################# 
#folder="$1";
#threshold="$2";
#filenum="$3";
#condition="$4";

##################################

ARGUMENT_LIST=(
    "Folder"
    "Threshold"
    "Filenum"
    "Condition"
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

        --Threshold)
            threshold=$2
            shift 2
            ;;

        --Filenum)
            filenum=$2
            shift 2
            ;;
        --Condition)
	    condition=$2
	    shift 2
	    ;;
	   
        *)
            break
            ;;
    esac
done


start=$(date +%s)
#Step 1
for file in *SJ.out.tab_strand_symbol_scaled;
do
     na1=${file%SJ.out.tab_strand_symbol_scaled}"_"idxstats_junction.txt;
     na2=${file%SJ.out.tab_strand_symbol_scaled}"_"junction_counts.txt;
     perl $folder/prepare_count_file_after_STAR_2_pass_mapping.pl $file UNION_junc_coor_with_junction_ID.txt $na1;
     awk '{print $1 "\t" $4 "\t" $2 "\t" $3 "\t" $10 "\t" $7 "\t" $9}' $na1 > $na2;
done

perl $folder/Identify_junctions_exist_in_certain_number_of_input_files.pl UNION_junc_coor_with_junction_ID.txt *_junction_counts.txt "$condition"_Junction_list_more_than_"$threshold"_read_in_at_least_"$filenum"_samples.txt $filenum $threshold
awk 'FNR==NR {arr[$0];next} ($5 in arr)' "$condition"_Junction_list_more_than_"$threshold"_read_in_at_least_"$filenum"_samples.txt UNION_junc_coor_with_junction_ID.txt > "$condition"_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$filenum"_samples.txt;
awk '{print $1 "\t" $4 "\t" $2 "\t" $3 "\t" $5}' "$condition"_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$filenum"_samples.txt > "$condition"_junc_coor_with_junction_ID_more_than_"$threshold"_read_in_at_least_"$filenum"_samples_formatted.txt;

end=$(date +%s)
runtime=$((end-start))

echo "JUM_2-2.sh finished in $runtime seconds."
echo "Command is:"
echo "bash JUM_2-2.sh --Folder $folder --Threshold $threshold --Filenum $filenum --Condition $condition"
