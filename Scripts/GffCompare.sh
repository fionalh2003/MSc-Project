for i in 4401 A1 A4 A5 B3 C2 DAOM G1 SL1
do
  grep -F -w -f 03-$i-NonCodingIDs.txt 05-$i-mapped.gff \
> 01-$i-lncRNA-Mapped.gff
wait
gffread 01-$i-lncRNA-Mapped.gff -T -o 01-$i-lncRNA-Mapped.gtf
done

# Running Gff compare
for i in 4401 A1 A4 A5 B3 C2 DAOM G1 SL1
do
gffcompare -r 02-${i}-gene-ann-mod.gtf -o 06-${i}-GFFCompare 01-$i-lncRNA-Mapped.gtf
done

# Filter for only lncRNA codes
for i in A4 A5 G1 SL1
do
awk '$0 ~ /class_code "[uixjop]"/' 04-$i-GFFCompare.annotated.gtf >05-$i-lncRNA_candidates.gtf
done

# Just transcript ID and code
for i in A4 A5 G1 SL1
do
awk '$0 ~ /class_code "[a-z]"/ {
  match($0, /transcript_id "([^"]+)"/, t);
  match($0, /class_code "([^"]+)"/, c);
  print t[1] "\t" c[1]
}' 05-$i-lncRNA_candidates.gtf > 06-$i-lncRNA_ClassCodes.txt
wait
sort -u 06-$i-lncRNA_ClassCodes.txt > 06-$i-lncRNA_UniqueClassCodes.txt
done

# Make a list of just transcript IDs
for i in A4 A5 G1 SL1
do
awk '{print $1}' 07-$i-lncRNA_UniqueClassCodes.txt > 08-$i-GffCompare-lncRNA.txt
done

# Seperate by haplotype
for i in A4 A5 G1 SL1
do
grep -wFf 08-$i-GffCompare-lncRNA.txt 03-$i-mapped_classification.filtered_lite_classification.txt | awk '{OFS="\t"} $2 ~ /Hap1/ {print $1, $2, $7}' > 09-$i-lncRNA-Hap1.txt
wait
grep -wFf 08-$i-GffCompare-lncRNA.txt 03-$i-mapped_classification.filtered_lite_classification.txt | awk '{OFS="\t"} $2 ~ /Hap2/ {print $1, $2, $7}' > 09-$i-lncRNA-Hap2.txt
done

# Find out how many lncRNA genes
cut -d'.' -f1,2 09-$i-lncRNA-Hap2.txt | sort -u | wc -l

# Make fasta of predicted lncRNA genes
for i in A4 A5 G1 SL1
do
 seqkit grep -n -f 08-$i-GffCompare-lncRNA.txt \
    01-$i-illumina-filtered-simple.fasta \
    > 88-$i-GffCompare-lncRNA-IDs.fasta
done

# Put strain name in header
sed '/^>/ s/$/-suffix/' input.fasta > output.fasta

# Merge Fasta files together
cat 888* > 11-lncRNA.fa

# cd-hit to cluster
cd-hit-est -i 11-lncRNA.fa -o 111-lncRNA.fasta  -c 0.90 -n 10 -T 8 -M 0
