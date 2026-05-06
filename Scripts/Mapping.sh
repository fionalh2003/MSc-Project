# change headers 
awk '/^>/ {print ">lncRNA-" ++i; next} {print}' 09-lncRNA-Candidates.fasta > 10-lncRNA-Candidates.fasta

# Mapping transcripts to the genomes
for i in A4-hap1 A4-hap2 A5-hap1 A5-hap2 G1-hap1 G1-hap2 SL1-hap1 SL1-hap2 4401 A1 B3 $
do
minimap2 -x splice:hq -uf 01-$i-genome.fasta 10-lncRNA-Candidates.fasta > 02-$i-MiniMa
done

# Filtering to keep high quality 
for i in A4-hap1 A4-hap2 A5-hap1 A5-hap2 G1-hap1 G1-hap2 SL1-hap1 SL1-hap2 4401 A1 B3 C2 DAOM 
do 
awk '$12 > 30' 02-$i-MiniMap.txt  > 03-$i-FilteredMiniMap.txt
done

# Simplify hetero minimap
for i in A4-hap1 A4-hap2 A5-hap1 A5-hap2 G1-hap1 G1-hap2 SL1-hap1 SL1-hap2
do  
    awk '{split($6, a, "_"); print $1, a[2], a[3]}' 03-$i-FilteredMiniMap.txt > 04-$i-lncRNAnStrainMiniMap.txt
done

# Simplify homo minimap
for i in 4401 A1 B3 C2 DAOM 
do
    awk '{print $1, substr($6, 1, index($6, "_") - 1)}' 03-$i-FilteredMiniMap.txt > 04-$i-lncRNAnStrainMiniMap.txt
done

# finding strain specific 
cat 04* > 05-MergedMapped.txt
awk '{print $1, $2}' 05-MergedMapped.txt | sort -u > 05-MergedStrain.txt
awk '{seen[$1][$2]=1; rows[$1]=rows[$1] $0 ORS} END{for(t in seen){n=0; for(s in seen[t]) n++; if(n==1) printf "%s", rows[t]}}' 05-MergedStrain.txt > 06-StrainSpecific.txt
awk '{count[$2]++} END{for(s in count) print s, count[s]}' 06-StrainSpecific.txt | sort

 # Finding shared transcripts
awk '{t[$1][$2]=1; strains[$2]=1} END{for(s in strains)n++; for(tr in t){c=0; for(s in t[tr])c++; if(c==n) print tr}}' 05-MergedStrain.txt > 06-StrainShared.txt
wc -l 06-Shared-lncRNA.txt

 # Finding haplotype specific
cat 04* > 05-HaplotypeMerge.txt
awk '{print $1, $2"-"$3}' OFS='\t' 05-HaplotypeMerge.txt > 06-Haplotype.txt
awk '{print $1, $2}' 06-Haplotype.txt | sort -u > 06-MergedHaplotype.txt
awk '{seen[$1][$2]=1; rows[$1]=rows[$1] $0 ORS} END{for(t in seen){n=0; for(s in seen[t]) n++; if(n==1) printf "%s", rows[t]}}' 06-MergedHaplotype.txt > 07-HaplotypeSpecific.txt
awk '{count[$2]++} END{for(s in count) print s, count[s]}' 07-HaplotypeSpecific.txt | sort
