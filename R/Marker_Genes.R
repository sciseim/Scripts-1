##----load embryonic gene expression data
load(file="/Users/Zhen/Dropbox/Single-cell/Figures/21May2014 updates for 3d PCA/pcaSave.embryonic.RData")
emb <-  recGenexp

##----load neonatal gene expression data
load(file="/Users/Zhen/Dropbox/Single-cell/Figures/21May2014 updates for 3d PCA/pcaSave.postnatal.RData")
neo <-  recGenexp

##----clean up
rm(pcaSave, recGenexp, clusterList)


##----Load marker gene list and gene list
Marker <- read.csv("/Users/Zhen/Data/012914_Marker_Gene_List(2).csv", header = T)
geneList.emb <- substring(rownames(emb), 20) #Ensembl ID number removed
geneList.neo <- substring(rownames(neo), 20) #ensembl ID number removed


#Mar.genes <- paste(mar.genes[,1], collapse = "|")
#Mar.genes <- grep(Mar.genes, geneList)
#Marker.Genes <- NORM.recRPKM.singles[which(rownames(NORM.recRPKM.singles) %in% Mar.genes),]

for (i in 1:ncol(Marker)){
  genes <- Marker[i]
  Mar.genes <- paste(mar.genes[,i], collapse = "|");
  Mar.genes <- grep(Mar.genes, geneList, ignore.case = T);
  Mar.genes <- NORM.recRPKM.singles[which(rownames(NORM.recRPKM.singles) %in% Mar.genes),];
  Marker.Genes <- rbind(Marker.Genes, Mar.genes)
}

rownames(Marker.Genes) <- geneList[as.numeric(rownames(Marker.Genes))]

write.csv(Marker.Genes, "012914_Marker_Genes_Embryonics_Singles_QNFirst.csv", sep="\t")

read.csv("")