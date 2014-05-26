##-------use pvclust to cluster cells with parallel processing on Louise-----------
##-------11/15/13


recRPKM2 <- read.table("recRPKM2.txt", header = TRUE)
genes <- read.table("~/Zhen/brainSeq.singleCell.130925_SN7001318_0091_BC2NALACXX.wholeGene.geneComposite.gene.RPKM.txt", header = T)
MTGenes <- grep("MT-", genes[,1])
noMTrecRPKM2 <- recRPKM2[-which(rownames(recRPKM2) %in% MTGenes),]
library("pvclust")
library("parallel")
cl <- makeCluster(10)
clustPv2 <-parPvclust(cl, noMTrecRPKM2, method.dist="cor", method.hclust="complete", nboot=10000) #noMTrecRPKM2: Mitochondria genes & MMB3_3000 removed
save(clustPv2, file = paste(Sys.Date(), "_pvclust_NO_3_3000_NO_MT_Genes.Rdata", sep = ""))
stopCluster(cl)

