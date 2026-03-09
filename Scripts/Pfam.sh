# Download database 
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# Unzip
gunzip Pfam-A.hmm.gz

# Convert to ORFs
TransDecoder.LongOrfs -t 11-Clustered-lncRNA.fasta

# Run Pfam on the ORFs
hmmscan \
  --cpu 4 \
  --domtblout pfam.domtblout \
  Pfam-A.hmm \
  longest_orfs.pep > pfam.log

# Extract sequences with the Pfam hits
grep -v "^#" pfam.domtblout | awk '{print $4}' | sort | uniq > pfam_hits.txt

# Clean up txt file
sed 's/\.[^.]*$//' pfam_hits.txt | sort -u > pfam_transcripts.txt

# Remove hits from fasta 
seqkit grep -v -f pfam_transcripts.txt 11-Clustered-lncRNA.fasta > 111-noPfam_transcripts.fa

# Make a file with just headers 
grep "^>" 111-noPfam_transcripts.fa | sed 's/^>//' > 111-Heterokaryon-noPfam_transcriptsIDs.txt
