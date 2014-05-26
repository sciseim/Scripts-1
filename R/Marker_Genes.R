library("limma")
recRPKM.singles <- read.table("010914_recRPKM_NO_MT_Embryonic_Singles.txt")
NORM.recRPKM.singles = normalizeBetweenArrays(as.matrix(recRPKM.singles[,-1]),method="quantile");
geneList <- read.table("~/../Dropbox/Single-cell/Genexp_25Nov2013/brainSeq.singleCell.130925_SN7001318_0091_BC2NALACXX.wholeGene.geneComposite.gene.RPKM.txt", 
                       header = T) [,1]
singles <- NORM.recRPKM.singles
rownames(singles) <- geneList[as.numeric(rownames(singles))]
mar.genes <- read.csv("012914_Marker_Gene_List.csv", header = T)
Mar.genes <- paste(mar.genes[,1], collapse = "|")
Mar.genes <- grep(Mar.genes, geneList)
Marker.Genes <- NORM.recRPKM.singles[which(rownames(NORM.recRPKM.singles) %in% Mar.genes),]

for (i in 2:ncol(mar.genes)){
  Mar.genes <- paste(mar.genes[,i], collapse = "|");
  Mar.genes <- grep(Mar.genes, geneList, ignore.case = T);
  Mar.genes <- NORM.recRPKM.singles[which(rownames(NORM.recRPKM.singles) %in% Mar.genes),];
  Marker.Genes <- rbind(Marker.Genes, Mar.genes)
}

rownames(Marker.Genes) <- geneList[as.numeric(rownames(Marker.Genes))]

write.csv(Marker.Genes, "012914_Marker_Genes_Embryonics_Singles_QNFirst.csv", sep="\t")

read.csv("")