for i in A4 A5 G1 SL1
do

grep -wFf 03-$i-cpat-lnRNA-pred_IDs.txt 03-$i-mapped_classification.filtered_lite_classification.txt | awk '{OFS="\t"} $2 ~ /Hap1/ {print $1, $2, $7}' > 33-$i-lncRNA-Hap1.txt

wait

grep -wFf 03-$i-cpat-lnRNA-pred_IDs.txt 03-$i-mapped_classification.filtered_lite_classification.txt | awk '{OFS="\t"} $2 ~ /Hap2/ {print $1, $2, $7}' > 33-$i-lncRNA-Hap2.txt

done

