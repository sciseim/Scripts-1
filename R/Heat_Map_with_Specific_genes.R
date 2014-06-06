load("21May2014 updates for 3d PCA/pcaSave.postnatal.RData")                         ##-----load embryonic gene expression data
rm("pcaSave")                                                                        ##-----remove unuserful data
genes.recRPKM <- recGenexp                                                      
Marker.Genes <- c()

mar.genes <- read.csv("/Users/Zhen/Data/012914_Marker_Gene_List.csv", header = T)    ##-----read marker genes

row.colorPalette <- topo.colors(ncol(mar.genes))
row.color <- c()
for (i in 1:ncol(mar.genes)){
  Mar.genes <- paste(mar.genes[,i], collapse = "|");
  Mar.genes <- grep(Mar.genes, rownames(genes.recRPKM));
  Mar.genes <- genes.recRPKM[Mar.genes,];
  row.color <- append(row.color, rep(row.colorPalette[i], each = nrow(Mar.genes)));
  Marker.Genes <- rbind(Marker.Genes, Mar.genes);
}

#mergeCluster_pre = c("5","16","4_6","8_20","10_11_19","14_15","3_18_1_7","2_13_9")
#mergeCluster <- mergeCluster_pre
mergeCluster_post = c("5_14_30_10_22","18_4_6_8_27","21_12_26","15_7","11","16_28","3_13_2_20_9_29_23_24");
mergeCluster <- mergeCluster_post

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

m <- log2(Marker.Genes+1)

library("gplots")
library("RColorBrewer")
pairs.breaks <- c(seq(-5, -1, length.out=50), seq(-1, 1, length.out =50),
                  seq(1, 5, length.out=50))
mycol <- colorpanel(n=149,low="blue", mid="white", high="red")


heatmap.2(as.matrix(m[,recIndex]), breaks = pairs.breaks, 
          col=mycol, 
          key=TRUE, keysize = 1, density.info = "histogram", symkey=TRUE, main="", 
          trace="none", cexRow=0.7, cexCol = 0.7, 
          labRow = substr(rownames(m), 20, 100),
          labCol = substr(colnames(m), 4, 100),
          Rowv = NULL, Colv = NULL , #as.dendrogram(hca), 
          dendrogram = "none", na.rm = TRUE, 
          ColSideColors = col.color,
          RowSideColors = row.color,
          scale = "row"
) # Col.Color = 1, Heat Map = 2, Key  = 3, Col.dendrogram = 4,