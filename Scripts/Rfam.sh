# Download database
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz

# Unzip
gunzip Rfam.cm.gz

# Index the models 
cmpress Rfam.cm

# Run Rfam 
cmscan \
  --cpu 20 \
  --rfam \
  --cut_ga \
  --nohmmonly \
  --tblout rfam.tbl \
  Rfam.cm \
  111-noPfam_transcripts.fa \
  > rfam.out

# Extract Rfam hits
grep -v "^#" rfam.tbl | awk '{print $3}' | sort | uniq > rfam_hits.txt

# Removing hits
seqkit grep -v -f rfam_hits.txt 111-noPfam_transcripts.fa > noRfam_lncRNA.fa

# Make a file with just the headers
grep "^>" noRfam_lncRNA.fa | sed 's/^>//' > 111-noRfam_transcriptsIDs.txt
