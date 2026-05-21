# length filtering  (00-LengthFilter.sh)
for i in 4401 A1 A4 A5 B3 C2 DAOM G1 SL1 
do 
seqkit seq -m 200 01-$i-illumina-filtered-simple.fasta > 02-$i-length-filtered.fasta
wait 
seqkit seq -n 02-$i-length-filtered.fasta > 02-$i-length-filtered-IDs.txt
done 

# filter gff (00-FilterGff.sh)
for i in 4401 A1 A4 A5 B3 C2 DAOM G1 SL1
do
  grep -F -w -f 02-$i-length-filtered-IDs.txt 05-$i-mapped.gff \
> 01-$i-lncRNA-Mapped.gff
wait
gffread 01-$i-lncRNA-Mapped.gff -T -o 01-$i-lncRNA-Mapped.gtf
done

# Running Gff compare (00-Gffcompare.sh)
for i in 4401 A1 A4 A5 B3 C2 DAOM G1 SL1
do
gffcompare -r 02-${i}-gene-ann-mod.gtf -o 03-${i}-GFFCompare 01-$i-lncRNA-Mapped.gtf
done

# Filter for only lncRNA codes and at least 2 exons (00-lncRNA-codes.sh)
for i in 4401 A1 A4 A5 B3 C2 DAOM G1 SL1
do

awk '

# Transcript lines with desired class_code
$3 == "transcript" {

    if (match($0, /class_code "([a-z])"/, a)) {

        if (a[1] ~ /^[uixo]$/) {

            match($0, /transcript_id "([^"]+)"/, t)
            tid = t[1]

            keep[tid] = 1
        }
    }
}

# Count exons for kept transcripts
$3 == "exon" {

    match($0, /transcript_id "([^"]+)"/, t)
    tid = t[1]

    if (keep[tid])
        exon_count[tid]++
}

# Store all lines
{
    match($0, /transcript_id "([^"]+)"/, t)
    tid = t[1]

    all_lines[tid] = all_lines[tid] $0 "\n"
}

END {

    for (tid in keep) {

        if (exon_count[tid] >= 2)

            printf "%s", all_lines[tid]
    }
}

' 03-$i-GFFCompare.annotated.gtf > 04-$i-lncRNA_candidates.gtf

done


# Just transcript ID and code (00-GffIDs.sh)
for i in 4401 A1 A4 A5 B3 C2 DAOM G1 SL1
do
awk '$0 ~ /class_code "[a-z]"/ {
  match($0, /transcript_id "([^"]+)"/, t);
  match($0, /class_code "([^"]+)"/, c);
  print t[1] "\t" c[1]
}' 04-$i-lncRNA_candidates.gtf > 06-$i-lncRNA_ClassCodes.txt
wait
sort -u 06-$i-lncRNA_ClassCodes.txt > 06-$i-lncRNA_UniqueClassCodes.txt
wait
awk '{print $1}' 06-$i-lncRNA_UniqueClassCodes.txt > 07-$i-GffCompare-lncRNA.txt
done



# Make fasta of predicted lncRNA genes
for i in A4 A5 G1 SL1
do
 seqkit grep -n -f 09-$i-GffCompare-lncRNA.txt \
    01-$i-illumina-filtered-simple.fasta \
    > 10-$i-GffCompare-lncRNA.fasta
done

# Put strain name in header
sed '/^>/ s/$/-suffix/' input.fasta > output.fasta (11)

# Merge Fasta files together
cat 11* > 12-lncRNA.fa

# cd-hit to cluster
cd-hit-est -i 12-lncRNA.fasta -o 13-Clustered-lncRNA.fasta -c 0.98 -aS 0.9 -aL 0.9 -G 0 -g 1 -T 16 -M 0
