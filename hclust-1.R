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
t=TRUE
file_name='presence_absense.pdf'
i=2
while (t==TRUE)
{
  if(file.exists(file_name)){
    file_name=paste('presence_absense_run_', toString(i),sep = '')
  }
  else{
    t=FALSE
  }
  i=i+1
}
  
pdf("presence_absense.pdf") 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
eclust<-read.delim(args[1], row.names = 1, sep = "\t", header = TRUE)

if(args[2] ==FALSE){
  eclust.matrix <- as.matrix(eclust[vapply(eclust, function(x) length(unique(x)) > 1, logical(1L))]
  ) 
} else{
  eclust.matrix <- as.matrix(eclust)}
#eclust[, colSums(eclust != 0) > 0] 

heatmap(eclust.matrix, 
        Rowv = NA, Colv = NULL, 
        scale = "none", 
        col = c("white", "darkorchid3"),
        labRow = NULL, 
        labCol = NULL,
        cexCol = .3,
        cexRow = .3,
        margins = c(10,10))
dev.off()
#coldendrogram <- color_branches(coldendrogram, k = 5, col = c("darkblue", "darkred", "darkred", "forestgreen", "black"))
'
#binary2<-read.phylip.data("~/iclouddrive/hclust/binary.phylip")
#bin<-read.delim("~/iclouddrive/bioinformatics-cloud/collimonas/out_rf_500/binary_table.txt", row.names = 1, sep = "\t", header = TRUE)
bin<-read.delim("~/iclouddrive/bioinformatics-cloud/hclust/fo159_2000/binary_table.txt", row.names = 1, sep = "\t", header = TRUE)
bin.acc.x = bin[rowSums(bin) < 143,]
#bin.acc.25.41 = bin.acc.41[rowSums(bin.acc.41) > 25,]
#in r, the square bracket is [row,column]
bin.acc.matrix <- as.matrix(bin.acc.151)
dist.fof2 <- dist(bin.acc.matrix, method = "binary", diag = FALSE, upper = FALSE)
hclust.fof2<-hclust(dist.fof2, method = "average", members = NULL)

heatmap(bin.acc.matrix, 
        Rowv = NULL, Colv = NULL, 
        scale = "none", 
        col = c("black", "white"),
        labRow = NA, 
        labCol = NULL,
        cexCol = 0.3)

?bin.acc.dist<-dist(bin.acc, method = "binary")
bin.acc.dist.hclust<-hclust(bin.acc.dist, method = "average", members = NULL)
plot(bin.acc.matrix.wide.dist.hclust, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = 1, ylab = "Height")

###
bin.acc.wide<-t(bin.acc)
bin.acc.matrix.wide <- as.matrix(bin.acc.wide)
bin.acc.matrix.wide.dist <- dist(bin.acc.matrix, method = "binary", diag = FALSE, upper = FALSE)
bin.acc.matrix.wide.dist.hclust <- hclust(distance, method="average")
plot(bin.acc.matrix.wide.dist.hclust, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = 1, ylab = "Height")
##now try heatmap to map presence/absence data
plot(bin.acc)

dendrogram  = as.dendrogram(cluster)
Rowv        = rowMeans(bin.acc.matrix, na.rm = T)
dendrogram  = reorder(dendrogram, Rowv) %>% set("branches_lwd", 2) #%>% set("branches_k_color", k = 10)
'
