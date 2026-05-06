# list of lncRNA
grep "^>" 10-lncRNA-Candidates.fasta | sed 's/>//' > 01-lncRNACandidateIDs.txt

# Filter list for strains
for i in A4 A5 G1
do
minimap2 -x splice:hq -uf 01-$i-illumina-filtered-simple.fasta 10-lncRNA-Candidates.fasta > 02-$i-MiniMap.txt
done
for i in A4 A5 G1
do
awk '$12 > 30' 02-$i-MiniMap.txt  > 03-$i-FilteredMiniMap.txt
done

# filter TPM to only include lncRNA
awk 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' 04-$i-lncRNA.txt 20-$i-super-table-iso-class-TPM.txt > 05-$i-lncRNA-TPM.txt
for i in A4 A5 G1
do
    awk '{print $6, substr($6, 1, index($6, "_") - 1)}' 03-$i-FilteredMiniMap.txt > 04-$i-lncRNAExpression.txt
done

# filter for what is in transcript data
 grep -wFf <(cut -f1 08-G1-hap1-iso-classified-best.txt) 04-G1-lncRNAExpression.
txt > 05-G1Hap1-Expressed.txt

# filter table
 grep -wFf 05-A4Hap1-Expressed.txt 20-A4-hap1-super-table-iso-class-TPM.txt > 06-A4Hap1-Expressed-lncRNA.txt
 
# only include transcript name and TPM for each host
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2
do
cut -f1,8,9,10 06-$i-Expressed-lncRNA.txt > 01-$i-TPM.txt
done

# Only include values where expression is at least one TPM in one of the 3 hosts
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2
do
awk 'NR==1 || $2>=1 || $3>=1 || $4>=1' 01-$i-TPM.txt > 02-$i-FilteredTPM.txt
done

# log transformation
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2
do
awk '{
    for(i=2;i<=4;i++) $i = log($i+1)/log(2);
    print
}' OFS="\t" 02-$i-FilteredTPM.txt > 03-$i-TPMLog2.txt

done

# log2 fold change to approximate differential expression
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2
do
awk '{
    fc12 = $2 - $3;
    fc13 = $2 - $4;
    fc23 = $3 - $4;
    print $0, fc12, fc13, fc23
}' OFS="\t" 03-$i-TPMLog2.txt > 04-$i-TPMLog2FoldChange.txt

done

# Upregulated
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2
do
awk '$5 >= 1 && $6 >= 1' 04-$i-TPMLog2FoldChange.txt > 05-$i-Cond1-Up.txt
wait
awk '$5 <= -1 && $7 >= 1' 04-$i-TPMLog2FoldChange.txt > 05-$i-Cond2-Up.txt
wait
awk '$6 <= -1 && $7 <= -1' 04-$i-TPMLog2FoldChange.txt > 05-$i-Cond3-Up.txt
done

# host specific
for i in A4Hap1 A4Hap2 A5Hap1 A5Hap2 G1Hap1 G1Hap2
do
awk '$2>=1 && $3<1 && $4<1' 02-$i-FilteredTPM.txt > 07-$i-Cond1Specific.txt
wait
awk '$2<1 && $3>=1 && $4<1' 02-$i-FilteredTPM.txt > 07-$i-Cond2Specific.txt
wait
awk '$2<1 && $3<1 && $4>=1' 02-$i-FilteredTPM.txt > 07-$i-Cond3Specific.txt
wait
awk '$2>=1 && $3<1 && $4>=1' 02-$i-FilteredTPM.txt > 07-$i-Cond1and3Specific.txt
wait 
awk '$2>=1 && $3>=1 && $4>=1' 02-$i-FilteredTPM.txt > 07-$i-AllConditions.txt
done
