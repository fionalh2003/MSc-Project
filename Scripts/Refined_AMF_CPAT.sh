# Generate hexamer table
make_hexamer_tab -c 10-RhizophagusIrregularisCDS.fasta -n 09-RhizophagusPredicted_lncRNA_IDs.fasta >11-RhizophagusIrregularisNewHexamerTable.tsv

# Train logistic model
make_logitModel -x 33-RhizophagusIrregularisNewHexamerTable.tsv -c 04-RhizophagusIrregularisRNA.fasta -n 09-RhizophagusPredicted_lncRNA_IDs.fasta -o 55-RhizophagusIrregularisNew

# Run 00-cpat.sh
for i in 4401 A1 B3 C2 DAOM A4 A5 G1 SL1
do

cpat -x 33-RhizophagusIrregularisNewHexamerTable.tsv \
-d 55-RhizophagusIrregularisNew.logit.RData --top-orf=5 \
-g 01-$i-illumina-filtered-simple.fasta -o 02-$i-cpat-lnRNA-pred

wait

awk 'NR>1 && $NF <= 0.54 {print $1}' 02-$i-cpat-lnRNA-pred.ORF_prob.best.tsv > 03-$i-cpat-lnRNA-pred_IDs.txt 

wait

seqkit grep -n -f 03-$i-cpat-lnRNA-pred_IDs.txt \
01-$i-illumina-filtered-simple.fasta \
> 05-$i-cpat-lnRNA-pred_IDs.fasta

done
