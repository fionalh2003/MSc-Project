# MSc-Project
## Pipelines used for my MSc project
### CPAT
Made a CPAT to predict fungal lncRNA with Fungi_CPAT.sh

Used that to predict AMF lncRNA then made AMF_CPAT.sh using output from Fungi_CPAT.sh

Made a new Final CPAT with second AMF lncRNA from AMF_CPAT.sh to make Refined_AMF_CPAT.sh

Refined_AMF_CPAT.sh runs CPAT on all 9 strains and filters to only include values with a coding potential of <0.54, make a list of the IDs and then making a fasta file with the predictied lncRNA

Seperating_Haplotype.sh Seperates the heterokaryons by haplotype

Count genes: cut -d'.' -f1,2 33-$i-lncRNA-Hap2.txt | sort -u | wc -l

