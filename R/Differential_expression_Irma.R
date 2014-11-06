###------Do anova analysis for the genes
###------Zhen Li

library("DESeq2")

batch <- c("green", "red", "no_color")

#--------Pool expression data in to all
all <- data.frame(matrix(nrow = 40944))

for(i in batch){
  infile <- read.csv(paste("/media/zhen/Expansion Drive/fastq/Irma/", i, "/", 
                           i, "_masterfile.csv", sep = ""), header = T)
  all <- data.frame(all, infile[, 3:ncol(infile)])
}

all <- all[, -1]
rownames(all) <- infile[,1]


#--------create metadata
metadata <- read.csv(paste("/media/zhen/Expansion Drive/fastq/Irma/", batch[1], "/", 
                     batch[1], "_read_depth.csv", sep = ""), header = T)

for(i in 2:length(batch)){
  infile <- read.csv(paste("/media/zhen/Expansion Drive/fastq/Irma/", batch[i], "/", 
                           batch[i], "_read_depth.csv", sep = ""), header = T)
  metadata <- rbind(metadata, infile)
}

meta <- metadata[which(metadata$Mapped_reads>500000 & metadata$FPKM.0.1>2000),]
good_ones <- which(metadata$Mapped_reads>500000 & 
                     metadata$FPKM.0.1>2000) #-------Identify good cells with 500,000 mapped reads and at least 2000 genes with FPKM > 0.1

diff <- data.frame(matrix(ncol = 4, dimnames = list(c(), c("Sum Sq","Mean Sq", "F value", "Pr(>F)"))))
for(i in 1:nrow(all)){
  gene <- data.frame(t(all[i, good_ones]), meta$Condition)
  gene_id <- colnames(gene)[1] #----------Keep the Ensembl id of the gene for later
  colnames(gene)[1] <- "gene"
  results <- aov(gene~meta.Condition, data=gene)
  anova <- anova(results)
  if(!is.nan(anova[1, 5]) & as.numeric(anova[1, 5]) < 0.001){
    diff <- rbind(diff, as.character(anova[1, 2:5]))
    rownames(diff)[nrow(diff)] <- gene_id #------------Assign the Ensembl id to the row name
  }
} 

diff <- diff[-1,]
rownames(infile) <- infile[,1]
geneID <- infile[rownames(diff), 2]
diff$gene_ID <- geneID
write.csv(diff, "ANOVA_results_Irma_0.001.csv")
