#!/usr/bin/env bash

#SBATCH --job-name="Bootstrap"
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-np

Family=$1
Phenotype=$2
Frequency=$3
echo working Bootstrapping analysis on Family: $Family, Phenotype: $Phenotype, and Frequency: $Frequency

Proband_number=$(ls -l | grep -c ^d )

echo Then number of probands is: $Proband_number
ml phevor2
ml R

cat */*AR.genes.tsv | sort -k 2 -r | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' | tail -n +$Proband_number | sed '1s/^.//' > All_AR.score
cat */*AD.genes.tsv | sort -k 2 -r | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' | tail -n +$Proband_number | sed '1s/^.//' > All_AD.score

R --vanilla  < Bootstrap_score.R

phevor_2 -b GO -f HPO -c /uufs/chpc.utah.edu/sys/installdir/phevor2/phevor.conf -s AR_To_Phevor.txt -i HP:0000819 > AR.phevor
phevor_2 -b GO -f HPO -c /uufs/chpc.utah.edu/sys/installdir/phevor2/phevor.conf -s AD_To_Phevor.txt -i HP:0000819 > AD.phevor

more */*html > All_variants.html

####Getting Variants
cat */*AR.variants.tsv > AR.variants
cat AR.phevor | sed -e "1d" | awk '{print $2 }' > AR.Genes_list 

head -n 1 AR.variants  | sed '1s/^.//'  > AR.Variants.Phevor
grep -F -f AR.Genes_list AR.variants >> AR.Variants.Phevor

R --vanilla  < Bootstrap_Html_parser.R
R --vanilla  < Bootstrap_Phevors_plot.R

mv AR.Html_variants.laz AR.$Family.Html_variants.$Phenotype.$Frequency.laz 
mv AR.Variants_plot.jpeg AR.$Family.Variants_plot.$Phenotype.$Frequency.jpeg
