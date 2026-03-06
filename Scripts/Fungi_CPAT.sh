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

# Extract FASTA sequences
seqkit grep -f 07-RhizophagusPredicted_lncRNA_IDs.txt \
04-RhizophagusIrregularisRNA.fasta \
> 08-RhizophagusPredicted_lncRNA_IDs.fasta

# 
