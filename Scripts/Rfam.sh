# Download database
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz

# Unzip
gunzip Rfam.cm.gz

# Index the models 
cmpress Rfam.cm

# Run Rfam 
cmscan \
  --cpu 16 \
  --rfam \
  --cut_ga \
  --nohmmonly \
  --tblout rfam.tbl \
  Rfam.cm \
  08-PfamFiltered.fasta \
  > rfam.out

# Extract Rfam hits
grep -v "^#" rfam.tbl | awk '{print $3}' | sort | uniq > rfam_hits.txt

# Removing hits
seqkit grep -v -f rfam_hits.txt 08-PfamFiltered.fasta > 09-lncRNA-Candidates.fasta

# Make a file with just the headers
grep "^>" 09-lncRNA-Candidates.fasta | sed 's/^>//' > 10-lncRNA-CandidatesIDs.txt
