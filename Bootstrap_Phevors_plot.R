library(boot)
library("biomaRt")
library(manhattanly)
library(dplyr)

library(qqman)

AR.Phevor <- read.table("AR.phevor", header = TRUE)

AR.Phevor$Phevor_scr <- AR.Phevor$PHEVOR_PRIOR*AR.Phevor$PHEVOR_SCORE

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

IDs <- AR.Phevor$GENE
genedesc <- getBM(attributes=c('external_gene_name','description','chromosome_name','start_position'), filters = 'external_gene_name', values = IDs, mart =ensembl)
colnames(genedesc)[1]="GENE"
genedesc$chromosome_name <- as.numeric(genedesc$chromosome_name)
genedesc <- filter(genedesc, is.na(chromosome_name)==FALSE)

AR.Phevor.Combined <- merge(AR.Phevor,genedesc, by = "GENE")

AR.Phevor.Combined$GENE <- as.character(AR.Phevor.Combined$GENE)

#searchFilters(mart = ensembl, pattern = "ensembl.*id")
#listFilters(mart=ensembl)
#listAttributes(mart = ensembl)

#manhattan(AR.Phevor.Combined, chr = "chromosome_name", bp = "start_position", p = "PHEVOR_SCORE", logp = FALSE, 
#          main = "XLC 3-35", ylim = c(0, 5), cex = 0.8, 
#          cex.axis = 0.9, col = c("blue4", "orange3","forestgreen","red","black", "purple","gold"), suggestiveline = F, genomewideline = F, 
#          chrlabs = c(1:22, "X"))

jpeg(filename = "AR.Plot")
manhattan(AR.Phevor.Combined, chr = "chromosome_name", bp = "start_position", p = "PHEVOR_SCORE", logp = FALSE, 
          col = c("darkorange2", "forestgreen", "firebrick3","mediumpurple3",
                                                  "lightsalmon4","orchid2",
                                                  "seashell4","yellow3",
                                                  "turquoise3","dodgerblue3"), ylim = range(0:6),
          cex.axis = 1, cex = 1.5, xaxs = "i",xlab='',yaxt="n",ylab='',
          bty="n",suggestiveline=F, genomewideline=F,main="Family 3-35", annotateTop = TRUE)

title(xlab="Chromosome", mgp=c(3,0.3,0), cex.lab=1.5, font.lab = 2)
title(ylab="Phevor Score", mgp=c(3,0.2,0), cex.lab=1.5, font.lab = 2)
box(lwd=2)
axis(2, at=seq(0,6,1),labels=seq(0,6,1), las=2, tck = -0.01,cex.axis=1)
axis(4, at=seq(0,6,2),labels = FALSE,las=2, tck = 0.03)
abline(h=2, col="red", lty = 2, lwd = 2.5)

dev.off()

####################################
####Plot with Manhattanly

AR.Html_variants <- read.table("AR.Html_variants.laz", header = TRUE, sep = "\t")


AR.Html_variants$Chr <- as.numeric(AR.Html_variants$Chr)
AR.Html_variants$Pos <- as.numeric(AR.Html_variants$Pos)

jpeg(filename = "AR.Variants_plot.jpeg")

manhattan(AR.Html_variants, chr = "Chr", bp = "Pos", p = "PHEVOR_SCORE",snp = "SNP", logp = FALSE, 
          col = c("darkorange2", "forestgreen", "firebrick3","mediumpurple3",
                  "lightsalmon4","orchid2",
                  "seashell4","yellow3",
                  "turquoise3","dodgerblue3"), ylim = range(0:6),
          cex.axis = 1, cex = 1.5, xaxs = "i",xlab='',yaxt="n",ylab='',
          bty="n",suggestiveline=F, genomewideline=F,main="Family")

title(xlab="Chromosome", mgp=c(3,0.3,0), cex.lab=1.5, font.lab = 2)
title(ylab="Phevor Score", mgp=c(3,0.2,0), cex.lab=1.5, font.lab = 2)
box(lwd=2)
axis(2, at=seq(0,6,1),labels=seq(0,6,1), las=2, tck = -0.01,cex.axis=1)
axis(4, at=seq(0,6,2),labels = FALSE,las=2, tck = 0.03)
abline(h=2, col="red", lty = 2, lwd = 2.5)

dev.off()

