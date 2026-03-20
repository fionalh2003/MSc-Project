# filter TPM to only include lncRNA
awk 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' 04-$i-lncRNA.txt 20-$i-super-table-iso-class-TPM.txt > 05-$i-lncRNA-TPM.txt

# only include transcript name and TPM for each host
cut -f1,8,9,10 05-$i-lncRNA-TPM.txt > 06-$i-TPM.txt

# Only include values where expression is at least one TPM in one of the 3 hosts
awk 'NR==1 || $2>=1 || $3>=1 || $4>=1' 06-$i-TPM.txt > 07-$i-FilteredTPM.txt
