#!/bin/bash
set -e

for file in *SJ.out.tab;
do
    name=$file"_"unannotated;
    awk '$6 == 0' $file > $name;
done

cat *SJ.out.tab_unannotated > combined_SJ_out_tab_unannotated.txt;

awk '{if($4==1) $4="+"; else if($4==2) $4="-"; print $1 "\t" $2 "\t" $3 "\t" $4}' combined_SJ_out_tab_unannotated.txt > combined_SJ_out_tab_unannotated_for_2nd_pass_genome_generation.txt;

rm *SJ.out.tab_unannotated;
rm combined_SJ_out_tab_unannotated.txt;

#STAR --runThreadN 2 --runMode genomeGenerate --genomeDir /home/qingqing/Vector/circle/qingqing_wang/concatenated/genome_index_2_pass --genomeFastaFiles /home/qingqing/Vector/circle/qingqing_wang/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa --sjdbFileChrStartEnd combined_SJ_out_tab_unannotated_for_2nd_pass_genome_generation.txt --sjdbGTFfile /home/qingqing/Vector/circle/qingqing_wang/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf --sjdbOverhang 49 &

#wait


