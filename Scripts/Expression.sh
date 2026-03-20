# filter TPM to only include lncRNA
awk 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' 04-$i-lncRNA.txt 20-$i-super-table-iso-class-TPM.txt > 05-$i-lncRNA-TPM.txt

# only include transcript name and TPM for each host
cut -f1,8,9,10 05-$i-lncRNA-TPM.txt > 06-$i-TPM.txt

# Only include values where expression is at least one TPM in one of the 3 hosts
awk 'NR==1 || $2>=1 || $3>=1 || $4>=1' 06-$i-TPM.txt > 07-$i-FilteredTPM.txt

# Converting to tau
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2
do 
awk '{
max=$2
for(j=3;j<=NF;j++) if($j>max) max=$j
done

# filter for host specific
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2 
do 
awk 'NR==1 || $NF>=0.8' 08-$i-TPM-t.txt > 09-$i-HostSpecific-tau.txt
done
