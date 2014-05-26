##---------------------use marker genes to generate heatmap for adult cells
##---------------------11/29/13

library("gplots")

##--------list of batch
batchList = c("131010_SN832_0147_AC2GW6ACXX",
              "131104_SN832_0151_AC2YTJACXX");

##---------gene expression
recRPKM = c();
recCount = c();
recRead2GFP = c();
recSamplist = c();

for(i in 1:length(batchList)){
  batchName = as.character(batchList[i]);
  cat(batchName,"\n");
  
  ##-----gene expression
  inFile1 = paste("~/../Dropbox/Single-cell/Genexp_25Nov2013/brainSeq.singleCell", batchName, "wholeGene.geneComposite.gene.RPKM.txt",sep=".");
  geneRPKM = read.table(inFile1,header=T);
  
  inFile2 = paste("~/../Dropbox/Single-cell/Genexp_25Nov2013/brainSeq.singleCell", batchName, "wholeGene.geneComposite.gene.count.txt",sep=".");
  geneCount = read.table(inFile2,header=T);  
  
  ##------good barcode
  inFile3 = paste("~/../Dropbox/Single-cell/Genexp_25Nov2013/",batchName,".QC_goodBarcode.txt",sep="");
  goodBarcode = read.csv(inFile3,header=T,sep="\t");
  tmpIndex = grep("MMB1_|MMB2_|MMB3_|MMB6_",as.character(goodBarcode[,2]));
  if(length(tmpIndex) > 0){
    goodBarcode = goodBarcode[(-1)*tmpIndex,];
  }  
  
  ##----- read mapped to genome, spikein and GFP
  inFile4 = paste("~/../Dropbox/Single-cell/Genexp_25Nov2013/",batchName,".countMapping.txt",sep="");
  seqDepth = read.csv(inFile4,header=T,sep="\t");
  
  
  ##-----only consider good barcode
  myindex = match(as.character(goodBarcode[,2]), colnames(geneRPKM));
  myindex2 = match(as.character(goodBarcode[,2]), as.character(seqDepth[,1]));
  geneRPKM = geneRPKM [,myindex];
  geneCount = geneCount [,myindex];
  read2GFP = as.numeric(seqDepth[myindex2,9]);
  read2genome = as.numeric(seqDepth[myindex2,2]);
  
  ##-----recording
  recSamplist = c(recSamplist,as.character(goodBarcode[,2]));
  if(length(recRPKM) == 0){
    recRPKM = geneRPKM;
    recCount = geneCount;
    recRead2GFP = read2GFP/read2genome;
  }
  else{
    recRPKM = cbind(recRPKM,geneRPKM);
    recCount = cbind(recCount,geneCount);
    recRead2GFP = c(recRead2GFP,read2GFP/read2genome);
  }
}  

##-----------do  quantile normalization

library("limma")
NORMrecRPKM = normalizeBetweenArrays(as.matrix(recRPKM),method="quantile");
#boxplot(recRPKM,las=2,cex=0.5,main="Quantile normalization");

geneList <- read.table("~/../Dropbox/Single-cell/Genexp_25Nov2013/brainSeq.singleCell.131104_SN832_0151_AC2YTJACXX.wholeGene.geneComposite.gene.RPKM.txt", 
                       header = T) [,1]
genes.recRPKM <- recRPKM
rownames(genes.recRPKM) <- geneList

mar.genes <- read.csv("~/../Dropbox/Single-cell/111913_Specific_genes.csv", header = T)
Mar.genes <- paste(mar.genes[,16], collapse = "|")
Mar.genes <- grep(Mar.genes, geneList)
Marker.Genes <- genes.recRPKM[which(rownames(recRPKM) %in% Mar.genes),]

for (i in 17:ncol(mar.genes)){
  Mar.genes <- paste(mar.genes[,i], collapse = "|");
  Mar.genes <- grep(Mar.genes, geneList, ignore.case = T);
  Mar.genes <- genes.recRPKM[which(rownames(recRPKM) %in% Mar.genes),];
  Marker.Genes <- rbind(Marker.Genes, Mar.genes)
}

write.csv(Marker.Genes, file = paste(Sys.Date(), 
                                     "_Marker_genes_QNormalized_Adult.csv", sep = ""))

Marker.Genes <- read.csv(paste("~/../Dropbox/Single-cell/", paste(Sys.Date(), "_Marker_genes_QNormalized_Adult.csv", sep = ""), sep = ""), header = TRUE)

#heatmap(as.matrix(Marker.Genes), key=TRUE, col=topo.colors(100), main = "Layer 1 Marker Genes", labRow = rownames(Marker.Genes), cexRow = 0.5, cexCol = 0.5, scale = "row", Rowv = NA, Colv = T)

##--------do cluster analysis

mydist=as.dist(1-cor(Marker.Genes[, 3:ncol(Marker.Genes)],method="s"));
hca = hclust(mydist,method="complete"); 
plot(hca,axes=F,xlab="Gene Expression",ylab="",cex=0.5, 
     main = paste("Based on Correlation of Expressed Marker Genes"));
#rect.hclust(hca, k = 10, border = "red")
#rect.hclust(hca, h = 100, which = NULL, border = 3:4)

##--------generate heatmap with heatmap.2
m <- log10(Marker.Genes[,3:ncol(Marker.Genes)]+1)

library("RColorBrewer")
pairs.breaks <- c(seq(-5, -1, length.out=50), seq(-1, 1, length.out =50),
                  seq(1, 5, length.out=50))
mycol <- colorpanel(n=149,low="blue", mid="white", high="red")


heatmap.2(as.matrix(m), 
          #breaks = pairs.breaks, 
          col=mycol, 
          key=TRUE, keysize = 1, density.info = "histogram", symkey=T, main="", 
          trace="none", cexRow=0.7, cexCol = 0.7, labRow = Marker.Genes[,2], 
          Rowv = F, Colv = T, #as.dendrogram(hca), 
          dendrogram = "column", na.rm = TRUE, 
          scale = "column", 
          RowSideColors = as.character(Marker.Genes[,1])
) # Col.Color = 1, Heat Map = 2, Key  = 3, Col.dendrogram = 4,

#---------------heat map with untransformed RPKM
pairs.breaks <- c(seq(0, 1, length.out=50), seq(1, 1000, length.out = 100))
mycol <- colorpanel(n = 149,low="pink", high="red")

heatmap.2(as.matrix(Marker.Genes[, 3:ncol(Marker.Genes)]), #breaks = pairs.breaks, 
          col=topo.colors(1000), key=TRUE, keysize = 1, density.info = "histogram", 
          symkey=T, main="", trace="none", cexRow=0.7, cexCol = 0.7, 
          labRow = Marker.Genes[,2], Rowv = NULL, Colv = as.dendrogram(hca), 
          dendrogram = "column", na.rm = TRUE, scale = "column", 
          RowSideColors = as.character(Marker.Genes[,1]))
