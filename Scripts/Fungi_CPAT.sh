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
cpat.py \
-g 04-RhizophagusIrregularisRNA.fasta \
-d 05-CandidaAlbicans.logit.RData \
-x 03-CandidaAlbicansHexamerTable.tsv \
-o Rhizophagus_CPAT

# Filter predicted lncRNAs
awk 'NR>1 && $NF <= 0.54 {print $1}' Rhizophagus_CPAT.ORF_prob.best.tsv \
> 06-RhizophagusIrregularisCandidaPredicted_lncRNAIDs.txt

# Extract IDs
cut -f1 06-RhizophagusPredicted_lncRNA.tsv | tail -n +2 \
> 07-RhizophagusPredicted_lncRNA_IDs.txt

# Extract FASTA sequences from Candida prediction
seqkit grep -f 07-RhizophagusPredicted_lncRNA_IDs.txt \
04-RhizophagusIrregularisRNA.fasta \
> 08-RhizophagusPredicted_lncRNA_IDs.fasta

# Generate hexamer table using first predicted AMF
make_hexamer_tab -c 01-RhizophagusIrregularisCDS.fasta -n 02-RhizophagusIrregularisCandidalncRNA.fasta >03-RhizophagusIrregularisHexamerTable.tsv

# Train logistic model for AMF
make_logitModel -x 03-RhizophagusIrregularisHexamerTable.tsv -c 04-RhizophagusIrregularisRNA.fasta -n 02-RhizophagusIrregularisCandidalncRNA.fasta -o 05-RhizophagusIrregularis

# Run CPAT with first AMF 
cpat.py   -g 04-RhizophagusIrregularisRNA.fasta   -d 05-RhizophagusIrregularis.logit.RData   -x 03-RhizophagusIrregularisHexamerTable.tsv   -o 06-Rhizophagus_CPAT

# Filter predicted lncRNAs
awk 'NR>1 && $NF <= 0.54 {print $1}' 06-Rhizophagus_CPAT.ORF_prob.best.tsv > 07-RhizophagusPredicted_lncRNA_IDs.txt 

# Extract IDs
cut -f1 07-RhizophagusPredicted_lncRNA_IDs.txt | tail -n +2 > 08-RhizophagusPredicted_lncRNA_IDs.txt

# Extract FASTA sequences for first AMF
seqkit grep -f 08-RhizophagusPredicted_lncRNA_IDs.txt 04-RhizophagusIrregularisRNA.fasta > 09-RhizophagusPredicted_lncRNA_IDs.fasta

