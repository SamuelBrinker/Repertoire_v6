source("https://bioconductor.org/biocLite.R")
biocLite()

# instructions: "The dendrogram was produced by the statistical package R, 
#using the hclust function with binary distance and average linkage; 
#black indicates presence of a locus and white the absence of a locus."

#library("phylotools")
#library("PHYLOGR")
#library("ade4")
#library("gplots")
#library("ctc")
#library("extrafont")
#library("dendextend")
#just run 16-26
pdf("presence_absense.pdf") 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
eclust<-read.delim(args[1], row.names = 1, sep = "\t", header = TRUE)
eclust.matrix <- as.matrix(eclust)
heatmap(eclust.matrix, 
        Rowv = NA, Colv = NULL, 
        scale = "none", 
        col = c("white", "darkorchid3"),
        labRow = NULL, 
        labCol = NULL,
        cexCol = .35,
        cexRow = .6,
        margins = c(10,0))
dev.off()
