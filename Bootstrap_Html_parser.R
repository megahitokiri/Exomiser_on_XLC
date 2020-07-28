#install.packages("rvest")
library(rvest)

#install.packages('xml2')
library('xml2')

library(dplyr)
library(tidyr)
#https://www.datacamp.com/community/tutorials/r-web-scraping-rvest

library(boot)
library("biomaRt")
library(dplyr)
library(stringr)

#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

#setwd("/scratch/general/lustre/u1123911/Gencove_Pipeline/fam_565563/exomiser_results/DM/1_percent")

snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37

Variants <-read_html("All_variants.html")

get_html_Variants <- function(html){
  
  html %>% 
  # The relevant tag
  html_nodes('.col-sm-12') %>%      
    html_text() %>% 
    # Trim additional white space
    str_trim() %>%                       
    # Convert the list into a vector
    unlist()                             
}

get_html_Genes <- function(html){
  
  html %>% 
    # The relevant tag
    html_nodes('.col-sm-3') %>%      
    html_text() %>% 
    # Trim additional white space
    str_trim() %>%                       
    # Convert the list into a vector
    unlist()                             
}

Filtered_html_Variants <- get_html_Variants(Variants)

Filtered_html_Variants <- as.data.frame(sort(Filtered_html_Variants))
colnames(Filtered_html_Variants) <- "Data"
Filtered_html_Variants$Retain <- grepl(":g.",Filtered_html_Variants$Data,)

Filtered_html_Variants <- filter(Filtered_html_Variants, Retain == TRUE)

Genomic_Variants <- Filtered_html_Variants %>% separate(Data, c("Type", "Genomic_Pos"),sep = "chr")
Genomic_Variants_Positions <- Genomic_Variants %>% separate(Genomic_Pos, c("Chr", "Genomic_Pos"),sep = ":g.")
Genomic_Variants_Positions <- Genomic_Variants_Positions %>% separate(Genomic_Pos, c("Pos", "Extra_info"),sep = ">")

Genomic_Variants_Positions$Pos <- gsub("[a-zA-Z ]", "", Genomic_Variants_Positions$Pos)

Genomic_Variants_Positions <- distinct(Genomic_Variants_Positions, Chr, Pos, .keep_all = TRUE)


####GET SNP info

Genomic_Variants_Positions$Info <- paste0(Genomic_Variants_Positions$Chr,":",Genomic_Variants_Positions$Pos,"-",Genomic_Variants_Positions$Pos)

Info_list <- Genomic_Variants_Positions$Info

SNPs_list <- c("0:0-0","rs0000")
for (Range_Value in Info_list) {
  Result <- snpsByOverlaps(snps,Range_Value)
  SNP_Vector <- c(Range_Value,Result$RefSNP_id)
  SNPs_list <- rbind(SNP_Vector,SNPs_list)
}

SNPs_list <- as.data.frame(SNPs_list)
colnames(SNPs_list) <- c("Info","SNP")
Genomic_Variants_Positions_with_SNP <- merge(Genomic_Variants_Positions,SNPs_list, by = "Info")



#####Getting Genes for manhattan plot

######Ensemble retrieve
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


AR.Phevor <- read.table("AR.phevor", header = TRUE)

IDs <- AR.Phevor$GENE
genedesc <- getBM(attributes=c('external_gene_name','chromosome_name','start_position','end_position'), filters = 'external_gene_name', values = IDs, mart =ensembl)
colnames(genedesc)[1]="GENE"
genedesc$chromosome_name <- as.numeric(genedesc$chromosome_name)
genedesc <- filter(genedesc, is.na(chromosome_name)==FALSE)

AR.Phevor.Combined <- merge(AR.Phevor,genedesc, by = "GENE")

AR.Phevor.Combined$EXOMISER_GENE <- as.character(AR.Phevor.Combined$GENE)


AR.Variants.Phevor <- read.table("AR.Variants.Phevor", header = TRUE)

AR.Variants.Phevor <- distinct(AR.Variants.Phevor, CHROM, POS, ALT, .keep_all = TRUE)


AR.phevor_simple <- AR.Phevor.Combined[,c("GENE","PHEVOR_SCORE")]
colnames(AR.phevor_simple)[1] <- "EXOMISER_GENE"

AR.Data <- merge(AR.Variants.Phevor, AR.phevor_simple, by = "EXOMISER_GENE")
AR.Data$Info <- paste0(AR.Data$CHROM,":",AR.Data$POS,"-",AR.Data$POS) 

AR.Html_variants <- merge(Genomic_Variants_Positions_with_SNP, AR.Data, by = "Info", all.x = TRUE)
AR.Html_variants.simple <- AR.Html_variants[,c("EXOMISER_GENE","PHEVOR_SCORE","FUNCTIONAL_CLASS","Chr","Pos","REF","ALT","QUAL","GENOTYPE","SNP","COVERAGE","HGVS","EXOMISER_GENE_COMBINED_SCORE","CONTRIBUTING_VARIANT")]

AR.Html_variants.simple <-AR.Html_variants.simple[order(-AR.Html_variants.simple$PHEVOR_SCORE),]

write.table(AR.Html_variants.simple,"AR.Html_variants.laz", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

