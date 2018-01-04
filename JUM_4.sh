#!/bin/bash
#set -e
# $1=directory for JUM scripts; $2=pvalue or pvalue-adjusted; $3=pvalue or padj threshold; $4=condition_file numbers; $5=ctrl_file_numbers; 
#############

folder="$1";
pvalue_padj="$2";
cutoff="$3";
condition_num="$4";
ctrl_num="$5";
refflat="$6";
#############


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
    cat $intron_long_psi $intron_short_psi > $intron_total_psi;
    head -1 $cassett_psi > headerr;
    cat headerr $intron_total_psi > temp;
    mv temp $intron_total_psi;
    
####################################################################################################################
# choose dpsi and map gene names to AS event
####################################################################################################################

perl $folder/identify_gene_name_for_JUM_output_1.pl $refflat AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt temporary.txt;
perl $folder/identify_gene_name_for_JUM_output_2.pl temporary.txt unique_mapped_event.txt cassette_exon_AS_event_redundant_mapped_gene.txt;
perl $folder/identify_gene_name_for_JUM_output_3.pl unique_mapped_event.txt AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt > tem.txt;
mv tem.txt AS_differential_JUM_output_cassette_exon_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt;

perl $folder/identify_gene_name_for_JUM_output_1.pl $refflat AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_combined_sorted_with_dpsi.txt temporary.txt;
perl $folder/identify_gene_name_for_JUM_output_2.pl temporary.txt unique_mapped_event.txt intron_retention_AS_event_redundant_mapped_gene.txt;
perl $folder/identify_gene_name_for_JUM_output_3.pl unique_mapped_event.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_combined_sorted_with_dpsi.txt > tem.txt;
mv tem.txt AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt;

perl $folder/identify_gene_name_for_JUM_output_1.pl $refflat AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt temporary.txt;
perl $folder/identify_gene_name_for_JUM_output_2.pl temporary.txt unique_mapped_event.txt A5SS_AS_event_redundant_mapped_gene.txt;
perl $folder/identify_gene_name_for_JUM_output_3.pl unique_mapped_event.txt AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt > tem.txt;
mv tem.txt AS_differential_JUM_output_A5SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt;
  
perl $folder/identify_gene_name_for_JUM_output_1.pl $refflat AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt temporary.txt;
perl $folder/identify_gene_name_for_JUM_output_2.pl temporary.txt unique_mapped_event.txt A3SS_AS_event_redundant_mapped_gene.txt;
perl $folder/identify_gene_name_for_JUM_output_3.pl unique_mapped_event.txt AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt > tem.txt;
mv tem.txt AS_differential_JUM_output_A3SS_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt;

perl $folder/identify_gene_name_for_JUM_output_1.pl $refflat AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt temporary.txt;
perl $folder/identify_gene_name_for_JUM_output_2.pl temporary.txt unique_mapped_event.txt MXE_AS_event_redundant_mapped_gene.txt;
perl $folder/identify_gene_name_for_JUM_output_3.pl unique_mapped_event.txt  AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt > tem.txt;
mv tem.txt AS_differential_JUM_output_MXE_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt;

perl $folder/identify_gene_name_for_JUM_output_1.pl $refflat AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt temporary.txt;
perl $folder/identify_gene_name_for_JUM_output_2.pl temporary.txt unique_mapped_event.txt composite_AS_event_redundant_mapped_gene.txt;
perl $folder/identify_gene_name_for_JUM_output_3.pl unique_mapped_event.txt AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt > tem.txt;
mv tem.txt AS_differential_JUM_output_composite_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt;

   mkdir temp;
   mv *"$pvalue_padj"_"$cutoff".txt* temp;
   mv *sorted.txt* temp;
   mv *long_sorted_with_dpsi.txt temp;
   mv *short_sorted_with_dpsi.txt temp;
   mv *_short.txt* temp;
   mv *_long.txt* temp;
   mv headerr temp;
   rm temporary.txt;
   rm unique_mapped_event.txt;
   rm AS_differential_JUM_output_mixed_events_"$pvalue_padj"_"$cutoff"_sorted_with_dpsi.txt;
   rm AS_differential_JUM_output_intron_retention_"$pvalue_padj"_"$cutoff"_combined_sorted_with_dpsi.txt;
   mv total_mixed_event.txt total_composite_event.txt   

