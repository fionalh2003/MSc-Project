# Generate hexamer table
make_hexamer_tab -c 01-RhizophagusIrregularisCDS.fasta -n 02-RhizophagusIrregularisCandidalncRNA.fasta >03-RhizophagusIrregularisHexamerTable.tsv

# Train logistic model
make_logitModel -x 03-RhizophagusIrregularisHexamerTable.tsv -c 04-RhizophagusIrregularisRNA.fasta -n 02-RhizophagusIrregularisCandidalncRNA.fasta -o 05-RhizophagusIrregularis

# Run CPAT
cpat.py   -g 04-RhizophagusIrregularisRNA.fasta   -d 05-RhizophagusIrregularis.logit.RData   -x 03-RhizophagusIrregularisHexamerTable.tsv   -o 06-Rhizophagus_CPAT

# Filter predicted lncRNAs
awk 'NR>1 && $NF <= 0.54 {print $1}' 06-Rhizophagus_CPAT.ORF_prob.best.tsv > 07-RhizophagusPredicted_lncRNA_IDs.txt 

# Extract IDs
cut -f1 07-RhizophagusPredicted_lncRNA_IDs.txt | tail -n +2 > 08-RhizophagusPredicted_lncRNA_IDs.txt

# Extract FASTA sequences
seqkit grep -f 08-RhizophagusPredicted_lncRNA_IDs.txt 04-RhizophagusIrregularisRNA.fasta > 09-RhizophagusPredicted_lncRNA_IDs.fasta
