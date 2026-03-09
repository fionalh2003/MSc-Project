# MSc-Project
## Pipelines used for my MSc project
### CPAT
For CPAT.sh 
Made a CPAT to predict fungal lncRNA with Candida

Used what I got from the first CPAT to refine for AMF

Used that to predict AMF lncRNA then made AMF_CPAT.sh using output from Fungi_CPAT.sh

Used output from second CPAT to make a third more refined CPAT and ran it on all 9 strains and filters to only include values with a coding potential of <0.54, make a list of the IDs and then making a fasta file with the predictied lncRNA

### GFFCOMPARE
GffCompare.sh

Classify lncRNA based on their location with respect to coding genes 

Class codes uixjop considered to be lncRNA

### Pfam
