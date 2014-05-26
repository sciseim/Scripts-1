##Usage: R script to get compile gene expression in good barcodes of adult stage. 
  ##analyze top i most variable genes, with Heat Map, PCA and hcluster  
##AIM: combine different runs 
##Date:11/29/13

library("limma");
library("DESeq");
library("cqn");

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

##-----GFP ratio,set size and color
recRead2GFP = recRead2GFP/max(recRead2GFP);

myindex1 = which(recRead2GFP < 0.1);
myindex2 = which(recRead2GFP>= 0.1);
GFPcex = rep(0, length(recRead2GFP));
GFPcex[myindex2] = 10* recRead2GFP[myindex2];
GFPcex[myindex1] = min(GFPcex[myindex2]); 

GFPcol = rep("green", length(recRead2GFP));
GFPcol[myindex1] = "red";

##-----------remove low expressed genes
lowexpIndex=which(apply(recRPKM,1,function (x) quantile(x,probs =0.75)) < 1);  ##at least two points larger than 1
if(length(lowexpIndex) > 0){
  recRPKM = recRPKM[(-1)*lowexpIndex,];
}
#boxplot(recRPKM,las=2,cex=0.5,main="Remove no and low expressed genes");


##----------remove stable genes
lowstdIndex = which(apply(recRPKM,1,sd) < 20);  ##at least two points larger than 1
if(length(lowexpIndex) > 0){
  recRPKM = recRPKM[(-1)*lowstdIndex,];
}
#boxplot(recRPKM,las=2,cex=0.5,main="Remove no and low expressed genes");

#----------------Mitochondria gene names
genes <- read.table("~/../Dropbox/Single-cell/Genexp_25Nov2013/brainSeq.singleCell.130925_SN7001318_0091_BC2NALACXX.wholeGene.geneComposite.gene.RPKM.txt", header = T) 
MTGenes <- grep("MT-", genes[,1])
noMTrecRPKM <- recRPKM[-which(rownames(recRPKM) %in% MTGenes),]

##-----1------do  quantile normalization
NORMrecRPKM = normalizeBetweenArrays(as.matrix(noMTrecRPKM),method="quantile");
#boxplot(recRPKM,las=2,cex=0.5,main="Quantile normalization");

##----------do pca analysis
pca=prcomp(t(NORMrecRPKM));
#plot(pca,xlab="Principle components");

rot=pca$x;
sumData=summary(pca)[[1]];
ratio1=round(100*sumData[1]^2/sum(sumData^2),2);
ratio2=round(100*sumData[2]^2/sum(sumData^2),2);
ratio3=round(100*sumData[3]^2/sum(sumData^2),2);


##--set axis
minXaxis=min(rot[,1])*1.3;
maxXaxis=max(rot[,1])*1.3;
minYaxis=min(rot[,2])*1.3;
maxYaxis=max(rot[,2])*1.3;
minZaxis=min(rot[,3])*1.3;
maxZaxis=max(rot[,3])*1.3;

plot(rot[,1],rot[,2],type="p",xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis),
     xlab=paste("PC1 ",ratio1,"%",sep=""),ylab=paste("PC2 ",ratio2,"%",sep=""),pch=20,
     col=GFPcol,cex.lab=1,cex.axis=1,cex=GFPcex);

text(rot[,1],rot[,2],labels=recSamplist,cex=0.5);

#####--------do cluster analysis
colnames(NORMrecRPKM) = recSamplist;
mydist=as.dist(1-cor(NORMrecRPKM,method="s"));
hca = hclust(mydist,method="complete"); 
plot(hca,main="",axes=F,xlab="Gene Expression",ylab="",cex=0.5);
rect.hclust(hca, k = 16, border = "red") #k based on clustGap analysis
#rect.hclust(hca, h = 50, which = c(2,7), border = 3:4)

##-------make 3D PCA plot with first three PC

#library("scatterplot3d")
#minXaxis=min(rot[,1])*1.3;
#maxXaxis=max(rot[,1])*1.3;
#minYaxis=min(rot[,2])*1.3;
#maxYaxis=max(rot[,2])*1.3;
#minZaxis=min(rot[,3])*1.3;
#maxZaxis=max(rot[,3])*1.3;

#scatterplot3d(rot[,1],rot[,2], rot[,3], type="p",xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis), zlim=c(minZaxis,maxZaxis), xlab=paste("PC1 ",ratio1,"%",sep=""),ylab=paste("PC2 ",ratio2,"%",sep=""),zlab=paste("PC2 ",ratio3, "%", sep=""), pch=20);

##--------make interactive 3d plot
library(rgl)
plot3d(rot[,1],rot[,2], rot[,3],xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis),
       zlim=c(minZaxis,maxZaxis), xlab=paste("PC1 ",ratio1,"%",sep=""),
       ylab=paste("PC2 ",ratio2,"%",sep=""),
       zlab=paste("PC3 ",ratio3, "%", sep=""), col=GFPcol)
text3d(rot[,1],rot[,2], rot[,3], text=recSamplist,cex=0.5)

##----------------Find the top i expressed genes and make heatmap
geneList <- read.table(inFile1, header = T)[,1]


#dt = recRPKM
#max.dt <- dt[order(apply(dt, 1, max), decreasing = T),] 
#max.dt[max.dt == 0] <- NA 
#max.dt <- max.dt[order(max.dt[,1], max.dt[,2], max.dt[,3], decreasing = T), ] #reodering genes based on the expression in the "1" cell
#max.dt <- max.dt[,order(max.dt[1,],decreasing = T)] #reodering cells based on the expression of the "1" gene
#max.dt[is.na(max.dt)] <- 0 
#geneList <- read.table("brainSeq.singleCell.130925_SN7001318_0091_BC2NALACXX.wholeGene.geneComposite.gene.RPKM.txt", header = T)[,1] 
#max.geneList <- geneList[as.numeric(rownames(max.dt))]
#heatmap(as.matrix(max.dt)[10:190,], col=topo.colors(10000), main = paste("top", i, "genes", sep = " "),Rowv = T, Colv = NA, labRow = max.geneList[1:i], cexRow = 0.5, cexCol = 0.5, scale = "column")

##----------------Use genes with the top variance to make heat map
#-----------------Find most variable genes
Variance <- apply(NORMrecRPKM, 1, var)
Variance <- data.frame(geneList[as.numeric(rownames(NORMrecRPKM))], Variance)
Variance <- Variance[order(Variance[,2], decreasing = T),]

VarianceList <- read.table("~/../Dropbox/Single-cell/Genexp_25Nov2013/112513_Variance_Embyonic.txt", header = T)

i = 100
dt = recRPKM
Top.Var <- Variance[1:i,]
Top.Var <- dt[rownames(Top.Var),]
Top.Var.geneList <- VarianceList[1:i,1]

#-----------------Generate Heat map with top i most variable genes
heatmap(as.matrix(Top.Var), col=topo.colors(1000), 
        main = paste("top", i, "variable genes", sep = " "),Rowv = T, Colv = T, 
        labRow = Top.Var.geneList, cexRow = 0.5, cexCol = 0.5, scale = "row")

##----------------Use genes with the top variance to make PCA
##----------do pca analysis
pca=prcomp(t(Top.Var));
#plot(pca,xlab="Principle components");

rot=pca$x;
sumData=summary(pca)[[1]];
ratio1=round(100*sumData[1]^2/sum(sumData^2),2);
ratio2=round(100*sumData[2]^2/sum(sumData^2),2);
ratio3=round(100*sumData[3]^2/sum(sumData^2),2);


##--set axis
minXaxis=min(rot[,1])*1.3;
maxXaxis=max(rot[,1])*1.3;
minYaxis=min(rot[,2])*1.3;
maxYaxis=max(rot[,2])*1.3;
minZaxis=min(rot[,3])*1.3;
maxZaxis=max(rot[,3])*1.3;

plot(rot[,1],rot[,2],main = paste("top", i, "genes", sep = " "), 
     type="p",xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis), 
     xlab=paste("PC1 ",ratio1,"%",sep=""),ylab=paste("PC2 ",ratio2,"%",sep=""),
     pch=20, col=GFPcol,cex.lab=1,cex.axis=1,cex=GFPcex);

text(rot[,1],rot[,2],labels=recSamplist,cex=0.5);

library("gplots")
library("RColorBrewer")


#####--------do cluster analysis
#colnames(NORMrecRPKM) = recSamplist;
mydist=as.dist(1-cor(Top.Var,method="s"));
hca = hclust(mydist,method="complete"); 
plot(hca,axes=F,xlab="Gene Expression",ylab="",cex=0.5, main = paste("top", i, "genes", sep = " "));
#rect.hclust(hca, k = 10, border = "red")
#rect.hclust(hca, h = 100, which = NULL, border = 3:4)
}

##-------make 3D PCA plot with first three PC

#library("scatterplot3d")
#minXaxis=min(rot[,1])*1.3;
#maxXaxis=max(rot[,1])*1.3;
#minYaxis=min(rot[,2])*1.3;
#maxYaxis=max(rot[,2])*1.3;
#minZaxis=min(rot[,3])*1.3;
#maxZaxis=max(rot[,3])*1.3;

#scatterplot3d(rot[,1],rot[,2], rot[,3], type="p",xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis), zlim=c(minZaxis,maxZaxis), xlab=paste("PC1 ",ratio1,"%",sep=""),ylab=paste("PC2 ",ratio2,"%",sep=""),zlab=paste("PC2 ",ratio3, "%", sep=""), pch=20);

##--------make interactive 3d plot
#library(rgl)
#plot3d(rot[,1],rot[,2], rot[,3],xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis),zlim=c(minZaxis,maxZaxis), xlab=paste("PC1 ",ratio1,"%",sep=""),ylab=paste("PC2 ",ratio2,"%",sep=""),zlab=paste("PC3 ",ratio3, "%", sep=""), col=GFPcol)
#text3d(rot[,1],rot[,2], rot[,3], text=recSamplist,cex=0.5)

##------------remove some of the cells
sub.Top.Var <- Top.Var[,-which(colnames(Top.Var)%in% 
                                 c("MMB1_C4","MMB1_G9", "MMB1_D5", "MMB1_G3", 
                                   "MMB1_A7", "MMB1_A1", "MMB1_B1", "MMB1_G10",
                                   "MMB1_B5","MMB1_D11"))]

