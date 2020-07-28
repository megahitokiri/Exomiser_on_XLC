library(dplyr)
#setwd("/scratch/general/lustre/u1123911/Gencove_Pipeline/Exomiser/exomiser-cli-12.1.0/results/Fam_565563/DM")
##Reccesive
AR.genes <- read.table("All_AR.score", header = TRUE)

AR.genes <- group_by(AR.genes, GENE_SYMBOL)

AR.genes.summary <- summarize(AR.genes, Avg_Exomiser_Score = mean(EXOMISER_GENE_VARIANT_SCORE,na.rm=TRUE), Samples = n())
AR.max <- max(AR.genes.summary$Samples, na.rm = TRUE)-1
AR.genes.summary$New_Score <-(AR.genes.summary$Avg_Exomiser_Score*AR.genes.summary$Samples)/AR.max 
AR.genes.Phevor <- AR.genes.summary[,c("GENE_SYMBOL","New_Score")]

write.table(AR.genes.Phevor,"AR_To_Phevor.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

##Dominant
AD.genes <- read.table("All_AD.score", header = TRUE)

AD.genes <- group_by(AD.genes, GENE_SYMBOL)

AD.genes.summary <- summarize(AD.genes, Avg_Exomiser_Score = mean(EXOMISER_GENE_VARIANT_SCORE,na.rm=TRUE), Samples = n())
AD.max <- max(AD.genes.summary$Samples, na.rm = TRUE)-1
AD.genes.summary$New_Score <-(AD.genes.summary$Avg_Exomiser_Score*AD.genes.summary$Samples)/AD.max 
AD.genes.Phevor <- AD.genes.summary[,c("GENE_SYMBOL","New_Score")]

write.table(AD.genes.Phevor,"AD_To_Phevor.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
