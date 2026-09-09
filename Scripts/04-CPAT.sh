# Generate hexamer table for Candida
make_hexamer_tab -c 01-CandidaAlbicansCDS.fasta \
-n 02-CandidaAlbicanslncRNA.fasta \
> 03-CandidaAlbicansHexamerTable.tsv

# Train logistic model for Candida
make_logitModel \
-x 03-CandidaAlbicansHexamerTable.tsv \
-c 04-CandidaAlbicansRNA.fasta \
-n 02-CandidaAlbicanslncRNA.fasta \
-o 05-CandidaAlbicans

# Run CPAT with Candida model
cpat -x 03-CandidaAlbicansHexamerTable.tsv  -d  05-CandidaAlbicans.logit.RData  --top-orf=100  --antisense -g 10-RNAsamba-Joint.fasta -o Rhizophagus_CPAT

# Extract FASTA sequences from Candida prediction
seqkit grep -f Rhizophagus_CPAT.no_ORF.txt \
10-RNAsamba-Joint.fasta \
> 08-RhizophagusPredicted_lncRNA_IDs.fasta

# Generate hexamer table using first predicted AMF
make_hexamer_tab -c 01-RhizophagusIrregularisCDS.fasta -n 08-RhizophagusPredicted_lncRNA_IDs.fasta > 03-RhizophagusIrregularisHexamerTable.tsv

# Train logistic model for AMF
make_logitModel -x 03-RhizophagusIrregularisHexamerTable.tsv -c 04-RhizophagusIrregularisRNA.fasta -n 08-RhizophagusPredicted_lncRNA_IDs.fasta -o 05-RhizophagusIrregularis

# Run CPAT with first AMF 
cpat -x 03-RhizophagusIrregularisHexamerTable.tsv   -d  05-RhizophagusIrregularis.logit.RData  --top-orf=100  --antisense -g 10-RNAsamba-Joint.fasta -o 06-Rhizophagus_CPAT
  
# Extract FASTA sequences for first AMF
seqkit grep -f 06-Rhizophagus_CPAT.no_ORF.txt 10-RNAsamba-Joint.fasta > 99-RhizophagusPredicted_lncRNA_IDs.fasta

# Generate hexamer table with refined AMF predictions
make_hexamer_tab -c 01-RhizophagusIrregularisCDS.fasta -n 99-RhizophagusPredicted_lncRNA_IDs.fasta >11-RhizophagusIrregularisNewHexamerTable.tsv

# Train logistic model with refined AMF predictions
make_logitModel -x 11-RhizophagusIrregularisNewHexamerTable.tsv -c 04-RhizophagusIrregularisRNA.fasta -n 99-RhizophagusPredicted_lncRNA_IDs.fasta -o 55-RhizophagusIrregularisNew

# Run final  CPAT 00-cpat.sh
for i in 4401 A1 B3 C2 DAOM A4 A5 G1 SL1
do
cpat -x 11-RhizophagusIrregularisNewHexamerTable.tsv   -d  55-RhizophagusIrregularisNew.logit.RData  --top-orf=100  --antisense -g 03-$i-NonCoding.fasta -o 02-$i-cpat-lnRNA-pred
wait
seqkit grep -n -f 02-$i-cpat-lnRNA-pred.no_ORF.txt \
03-$i-NonCoding.fasta \
> 05-$i-cpat-lnRNA-pred_IDs.fasta
done


