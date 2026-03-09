# Mapping transcripts to the genomes
for i in A4-hap1 A4-hap2 A5-hap1 A5-hap2 G1-hap1 G1-hap2 4401 A1 B3 C2 DAOM 
do 
minimap2 -x splice:hq -uf 01-$i-genome.fasta 11-lncRNA.fa > 02-$i-MiniMap.txt
done

# Filtering to keep high quality 
for i in A4-hap1 A4-hap2 A5-hap1 A5-hap2 G1-hap1 G1-hap2 4401 A1 B3 C2 DAOM 
do 
awk '$12 > 30' 02-$i-MiniMap.txt  > 03-$i-FilteredMiniMap.txt
done

# Simplify hetero minimap
for i in A4-hap1 A4-hap2 A5-hap1 A5-hap2 G1-hap1 G1-hap2 SL1-hap1 SL1-hap2
do  
    awk '{split($6, a, "_"); print $1, a[2], a[3]}' 07-$i-FilteredMiniMap.txt > 11-$i-lncRNAnStrainMiniMap.txt
done

# Simplify homo minimap
for i in 4401 A1 B3 C2 DAOM 
do
    awk '{print $1, substr($6, 1, index($6, "_") - 1)}' 07-$i-FilteredMiniMap.txt > 11-$i-lncRNAnStrainMiniMap.txt
done
