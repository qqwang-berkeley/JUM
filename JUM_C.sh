#!/bin/bash
#set -e
# $1=directory for JUM scripts; $2=pvalue or pvalue-adjusted; $3=pvalue or padj threshold; $4=condition_file numbers; $5=ctrl_file_numbers; 
#############

#folder="$1";
#pvalue_padj="$2";
#cutoff="$3";
#condition_num="$4";
#ctrl_num="$5";
#refflat="$6";
#############

ARGUMENT_LIST=(
    "Folder"
    "Test"
    "Cutoff"
    "TotalCondition1FileNum"
    "TotalCondition2FileNum"
    "REF"
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

        --Test)
            pvalue_padj=$2
            shift 2
            ;;

        --Cutoff)
            cutoff=$2
            shift 2
            ;;

        --TotalCondition1FileNum)
            condition_num=$2
            shift 2
            ;;

        --TotalCondition2FileNum)
            ctrl_num=$2
            shift 2
            ;;
       
	--REF)
	    refflat=$2
	    shift 2
	    ;;

        *)
            break
            ;;
    esac
done

start=$(date +%s)

for output in AS_differential_JUM_output_A*SS*_"$pvalue_padj"_"$cutoff".txt
do
   sort_out=${output%.txt}"_"sorted.txt;
   psi_out=${sort_out%.txt}"_"with_dpsi.txt;
   sort -k1,1 -k2,2 -k3,3 $output > $sort_out;
   perl $folder/output_reformat_JUM_1.pl $sort_out $condition_num $ctrl_num > $psi_out;
   awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' $psi_out > temp;
   mv temp $psi_out;
done

   output=AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff".txt
   sort_out=${output%.txt}"_"sorted.txt;
   psi_out=${sort_out%.txt}"_"with_dpsi.txt;
   sort -k1,1 -k2,2 -k3,3 $output > $sort_out;
   perl $folder/output_reformat_JUM_1.pl $sort_out $condition_num $ctrl_num > $psi_out;
   awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' $psi_out > temp;
   mv temp $psi_out;

   output=AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff".txt
   sort_out=${output%.txt}"_"sorted.txt;
   psi_out=${sort_out%.txt}"_"with_dpsi.txt;
   sort -k1,1 -k2,2 -k3,3 $output > $sort_out;
   perl $folder/output_reformat_JUM_1.pl $sort_out $condition_num $ctrl_num > $psi_out;
   awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' $psi_out > temp;
   mv temp $psi_out;

    cassette=AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff".txt
    cassett_sort=${cassette%.txt}"_"sorted.txt;
    cassett_psi=${cassette%.txt}"_"sorted_with_dpsi.txt;
    sort -k1,1 -k2,2 -k3,3 $cassette > $cassett_sort;
    perl $folder/output_reformat_JUM_2.pl $cassett_sort $condition_num $ctrl_num > $cassett_psi;
    awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' $cassett_psi > temp;
    mv temp $cassett_psi;

    intron=AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff".txt
    intron_long=${intron%.txt}"_"long.txt
    intron_short=${intron%.txt}"_"short.txt
    intron_long_sort=${intron_long%.txt}"_"sorted.txt;
    intron_long_psi=${intron_long%.txt}"_"sorted_with_dpsi.txt;
    intron_short_sort=${intron_short%.txt}"_"sorted.txt;
    intron_short_psi=${intron_short%.txt}"_"sorted_with_dpsi.txt;
    intron_total_psi=${intron%.txt}"_"combined_sorted_with_dpsi.txt

    awk '/left|right/' $intron > $intron_long;
    awk '/short/' $intron > $intron_short;
    sort -k1,1 -k2,2 -k3,3 $intron_long > $intron_long_sort;
    perl $folder/output_reformat_JUM_3.pl $intron_long_sort $condition_num $ctrl_num > $intron_long_psi;
    sort -k1,1 -k2,2 -k3,3 $intron_short > $intron_short_sort;
    perl $folder/output_reformat_JUM_4.pl $intron_short_sort $condition_num $ctrl_num > $intron_short_psi;
    #cat $intron_long_psi $intron_short_psi > $intron_total_psi;
    head -1 $cassett_psi > headerr;
    #cat headerr $intron_total_psi > temp;
    cat headerr $intron_long_psi > temp;
    mv temp $intron_long_psi;
    cat headerr $intron_short_psi > temp;
    mv temp $intron_short_psi;
   

   echo "done pre-processing!" 
####################################################################################################################
# choose dpsi and map gene names to AS event
####################################################################################################################
dpsi_gene_func() {
	input=$1;
	output=${input%_sorted_with_dpsi.txt}"_"final.txt;
	redu1=${input#AS_differential_JUM_output_}
	redu2=${redu1%_sorted_with_dpsi.txt}"_"AS_event_redundant_mapped_gene.txt;
	temo=${redu1%_sorted_with_dpsi.txt}"_"temporary.txt;
        uniqu_map=${redu1%_sorted_with_dpsi.txt}"_"unique_mapped_event.txt;

        perl $folder/identify_gene_name_for_JUM_output_1.pl $refflat $input $temo;
        perl $folder/identify_gene_name_for_JUM_output_2.pl $temo $uniqu_map $redu2;
        perl $folder/identify_gene_name_for_JUM_output_3.pl $uniqu_map $input > $output;
}

echo "starting mapping AS events to genes"

dpsi_gene_func AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt &
dpsi_gene_func AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_long_sorted_with_dpsi.txt &
dpsi_gene_func AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_short_sorted_with_dpsi.txt &
dpsi_gene_func AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt &
dpsi_gene_func AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt &
dpsi_gene_func AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt &
dpsi_gene_func AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt &

wait;  
echo "completed mapping"

perl $folder/final_process_cassette_exon_output1.pl AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_final_simplified.txt cassette_exon_event_old_new_ID_pairing.txt;
perl $folder/final_process_cassette_exon_output2.pl cassette_exon_event_old_new_ID_pairing.txt AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_final_detailed.txt;

perl $folder/final_process_IR_output_long.pl AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_long_final.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_long_tem1.txt IR_long_event_old_new_ID_pairing.txt;
perl $folder/final_process_cassette_exon_output2.pl IR_long_event_old_new_ID_pairing.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_long_final.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_long_tem2.txt;
perl $folder/final_process_IR_output_short.pl AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_short_final.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_short_tem1.txt IR_short_event_old_new_ID_pairing.txt;
perl $folder/final_process_cassette_exon_output2.pl IR_short_event_old_new_ID_pairing.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_short_final.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_short_tem2.txt;
cat AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_long_tem1.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_short_tem1.txt > AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_final_simplified.txt;
cat AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_long_tem2.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_short_tem2.txt > AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_final_detailed.txt;

perl $folder/final_process_MXE_output.pl AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_final_simplified.txt MXE_event_old_new_ID_pairing.txt $condition_num $ctrl_num;
perl $folder/final_process_cassette_exon_output2.pl MXE_event_old_new_ID_pairing.txt AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_final_detailed.txt;

perl $folder/final_process_A5SS_output.pl AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_final_simplified.txt A5SS_event_old_new_ID_pairing.txt; 
perl $folder/final_process_cassette_exon_output2.pl A5SS_event_old_new_ID_pairing.txt AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_final_detailed.txt;

perl $folder/final_process_A3SS_output.pl AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_final_simplified.txt A3SS_event_old_new_ID_pairing.txt;
perl $folder/final_process_cassette_exon_output2.pl A3SS_event_old_new_ID_pairing.txt AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_final_detailed.txt;

perl $folder/final_process_composite_output.pl AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_composite_events_"$pvalue_padj"_"$cutoff"_final_simplified.txt composite_event_old_new_ID_pairing.txt;
perl $folder/final_process_cassette_exon_output2.pl composite_event_old_new_ID_pairing.txt AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff"_final.txt AS_differential_JUM_output_composite_events_"$pvalue_padj"_"$cutoff"_final_detailed.txt;

end=$(date +%s)
runtime=$((end-start))
echo "JUM_C.sh core execution finished in $runtime seconds."
echo "Command is:"
echo "bash JUM_C.sh --Folder $folder --Test $pvalue_padj --Cutoff $cutoff --TotalCondition1FileNum $condition_num --TotalCondition2FileNum $ctrl_num --REF $refflat"

if [ -d "temp_JUM_C_run" ]; then
rm -r temp_JUM_C_run;
fi

mkdir temp_JUM_C_run;
ls | grep -v temp_JUM_C_run | xargs mv -t temp_JUM_C_run;

mv temp_JUM_C_run/*detailed.txt .
mv temp_JUM_C_run/*simplified.txt .

