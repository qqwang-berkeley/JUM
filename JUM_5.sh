#!/bin/bash
#set -e
# $1=pvalue or pvalue-adjusted; $2=pvalue or padj threshold; $3=condition_file numbers; $4=ctrl_file_numbers; 
#############

pvalue_padj="$1";
cutoff="$2";
condition_num="$3";
ctrl_num="$4";
psi="$5";
#############

l=$(echo "$condition_num*2+$ctrl_num*2+17"|bc)

minus=$(echo "0 - $psi"|bc)

for output in AS_differential_JUM_output_A*SS*_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt
do
   cut_out=${output%.txt}"_"more_than_"$psi".txt;
   awk -v cut_off="$psi" '$l >= cut_off' $output | cut -f2 | sort -u > list1;
   awk -v cut="$minus" '$l <= cut' $output | cut -f2 | sort -u > list2;
   cat list1 list2 | sort -u > list;
   awk 'FNR==NR {arr[$0];next} ($2 in arr)' list $output > $cut_out;
done

   output=AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt
   cut_out=${output%.txt}"_"more_than_"$psi".txt;
   awk -v cut_off="$psi" '$$l >= cut_off' $output | cut -f2 | sort -u > list1;
   awk -v cut="$minus" '$$l <= cut' $output | cut -f2 | sort -u > list2;
   cat list1 list2 | sort -u > list;
   awk 'FNR==NR {arr[$0];next} ($2 in arr)' list $output > $cut_out;

   output=AS_differential_JUM_output_composite_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt
   cut_out=${output%.txt}"_"more_than_"$psi".txt;
   awk -v cut_off="$psi" '$$l >= cut_off' $output | cut -f2 | sort -u > list1;
   awk -v cut="$minus" '$$l <= cut' $output | cut -f2 | sort -u > list2;
   cat list1 list2 | sort -u > list;
   awk 'FNR==NR {arr[$0];next} ($2 in arr)' list $output > $cut_out;
    
   cassette=AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt
   cut_out=${output%.txt}"_"more_than_"$psi".txt;
   awk '$$l != ""' $cassette > temp_cassette;
   awk -v cut_off="$psi" '$$l >= cut_off' $temp_cassette | cut -f2 | sort -u > list1;
   awk -v cut="$minus" '$$l <= cut' $temp_cassette | cut -f2 | sort -u > list2;
   cat list1 list2 | sort -u > list;
   awk 'FNR==NR {arr[$0];next} ($2 in arr)' list $cassette > $cut_out;
   rm temp_cassette;

   intron=AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt
   cut_out=${output%.txt}"_"more_than_"$psi".txt;
   awk '$$l != ""' $intron > temp_intron;
   awk -v cut_off="$psi" '$$l >= cut_off' $temp_intron | cut -f2 | sort -u > list1;
   awk -v cut="$minus" '$$l <= cut' $temp_intron | cut -f2 | sort -u > list2;
   cat list1 list2 | sort -u > list;
   awk 'FNR==NR {arr[$0];next} ($2 in arr)' list $intron > $cut_out;
   rm temp_intron;

rm list;
rm list1;
rm list2;
