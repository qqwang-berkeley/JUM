#!/bin/bash
set -e

#start=$(date +%s);

for file in *SJ.out.tab;
do
	name1=$file"_"strand_symbol_scaled;
	name2=${file%SJ.out.tab}"_"SJ_coor.txt;
	awk '{if($4==1) $4="+"; else if($4==2) $4="-"; print $1 "\t" $2-1 "\t" $3-1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' $file > $name1;
        less $name1 | cut -f1,2,3,4 > $name2; 
done

cat *SJ_coor.txt | sort -u > UNION_junc_coor.txt
awk '{print $0 "\t" "Junction_" FNR}' UNION_junc_coor.txt > UNION_junc_coor_with_junction_ID.txt;

#end=$(date +%s)
#runtime=$((end-start))
#echo "JUM_2-1.sh finished in $runtime seconds."
