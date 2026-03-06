for i in A4 A5 G1 SL1
do

awk 'NR>1 {print $1}' 02-$i-cpat-lnRNA-pred.ORF_prob.best.tsv | \
grep -wFf - 03-$i-mapped_classification.filtered_lite_classification.txt | \
awk '{OFS="\t"} $2 ~ /Hap1/ {print $1, $2, $7}' \
> 44-$i-Hap1.txt

awk 'NR>1 {print $1}' 02-$i-cpat-lnRNA-pred.ORF_prob.best.tsv | \
grep -wFf - 03-$i-mapped_classification.filtered_lite_classification.txt | \
awk '{OFS="\t"} $2 ~ /Hap2/ {print $1, $2, $7}' \
> 44-$i-Hap2.txt

done
