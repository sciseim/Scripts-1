
library("gplots")
library("RColorBrewer")

#load("/Users/Zhen/Dropbox/Single-cell/Figures/21May2014 updates for 3d PCA/pcaSave.embryonic.RData")  ##----This is for embryonic
#load("/Users/Zhen/Dropbox/Single-cell/Figures/21May2014 updates for 3d PCA/pcaSave.postnatal.RData")  ##----This is for postnatal
recGenexp <- read.table("c:/Users/Zhen/Dropbox/Single-cell/Figures/10June2014/humanBrain.singleCell.wholeGene.geneComposite.gene.RPKM.txt", header = T)  ##----This is for human

##-----load embryonic gene expression data
rm("pcaSave")                                                                        ##-----remove unuserful data
genes.recRPKM <- recGenexp                                                      
Marker.Genes <- c()


mar.genes <- c("RELN", "PCP4", "BCL11B", "NOS1", "SOX5", "ASCL1", "FEZF2", "CUX1", "CUX2", "ZNF804A")

#mar.genes <- read.csv("/Users/Zhen/Documents/My Box Files/Manuscript/Marker_Gene_List.csv", header = T)    ##-----read marker genes

#mar.genes <- read.csv("/Users/Zhen/Documents/My Box Files/Manuscript/Manuscript_Draft/Marker_genes_from_Mihovil.csv", header = T)    ##-----read Mihovil's marker genes

#rownames(genes.recRPKM) <- substr(rownames(genes.recRPKM), 20, 100)  ##----This is for embryonic and postnatal mouse
genes.recRPKM$Geneid <- substr(genes.recRPKM$Geneid, 17, 100)  ##----This is for human

marker.genes <- which(genes.recRPKM$Geneid %in% mar.genes)

row.colorPalette <- topo.colors(ncol(mar.genes))
row.color <- c()
marker.genes <- c()
for (i in 1:ncol(mar.genes)){
  Mar.genes <- c()
  for(j in 1:length(mar.genes[,i])){
    Mar.genes <- c(Mar.genes, which(rownames(genes.recRPKM) %in% mar.genes[j, i]))
  }
  marker.genes <- c(marker.genes, Mar.genes)
  row.color <- c(row.color, rep(row.colorPalette[i], each = length(Mar.genes)));
}
  #Mar.genes <- paste(mar.genes[,i], collapse = "|");
  #Mar.genes <- grep(Mar.genes, rownames(genes.recRPKM));
Marker.Genes <- genes.recRPKM[marker.genes,];

#mergeCluster <- c("5","16","4_6","8_20","10_11_19","14_15","3_18_1_7","2_13_9");  ##----This is for embryonic
#mergeCluster = c("5_14_30_10_22","18_4_6_8_27","21_12_26","15_7","11","16_28","3_13_2_20_9_29_23_24");  ##----This is for postnatal
mergeCluster <- c("17_3", "24_16_13_14_19", "23_7_22_15", "18_4", "8_12_25", "5_20_9_10", "11_1_17", "21_2_6_11")  ##----This is for human

col.colorPalette <- rainbow(length(mergeCluster))
col.color <- c()
recIndex = c();
for(i in 1:length(mergeCluster)){
  rec = unlist(strsplit(as.character(mergeCluster[i]),split="_"));
  ##----
  for(j in 1:length(rec)){
    myindex = which(clusterList == rec[j]);
    col.color <- append(col.color, rep(col.colorPalette[i], times=length(myindex)));
    recIndex = c(recIndex, myindex);
  }
}

m <- log2(Marker.Genes[2:ncol(Marker.Genes)]+1)

pairs.breaks <- c(seq(-5, -1, length.out=50), seq(-1, 1, length.out =50),
                  seq(1, 5, length.out=50))
mycol <- colorpanel(n=149,low="blue", mid="white", high="red")


heatmap.2(as.matrix(m),
          #as.matrix(m[,recIndex]), 
          breaks = pairs.breaks, 
          col=mycol, 
          key=TRUE, keysize = 1, density.info = "histogram", symkey=TRUE, main="", 
          trace="none", cexRow=0.7, cexCol = 0.7, 
          labRow = Marker.Genes$Geneid,
          #labRow = substr(rownames(m), 20, 100),
          labCol = substr(colnames(m), 4, 100),
          Rowv = T, Colv = T, #as.dendrogram(hca), 
          dendrogram = "both", na.rm = TRUE, 
          #ColSideColors = col.color,
          #RowSideColors = row.color,
          scale = "row"
) # Col.Color = 1, Heat Map = 2, Key  = 3, Col.dendrogram = 4,