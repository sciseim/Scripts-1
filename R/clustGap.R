##---------------Use clustGap to determine optimal cluster number with top i variable genes
##---------------11/27/13

library("cluster")

##--------list of batch
batchList = c("fluidigm",
              "130925_SN7001318_0091_BC2NALACXX",
              "131010_SN832_0147_AC2GW6ACXX",
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
  inFile1 = paste("~/Genexp_25Nov2013/brainSeq.singleCell", batchName, "wholeGene.geneComposite.gene.RPKM.txt",sep=".");
  geneRPKM = read.table(inFile1,header=T);
  
 
  ##------good barcode
  inFile3 = paste("~/Genexp_25Nov2013/",batchName,".QC_goodBarcode.txt",sep="");
  goodBarcode = read.csv(inFile3,header=T,sep="\t");
  tmpIndex = grep("MMB5|MMB8|MMB11|MMB12",as.character(goodBarcode[,2]));
  if(length(tmpIndex) > 0){
    goodBarcode = goodBarcode[(-1)*tmpIndex,];
  }  
  
  ##----- read mapped to genome, spikein and GFP
  inFile4 = paste("~/Genexp_25Nov2013/",batchName,".countMapping.txt",sep="");
  seqDepth = read.csv(inFile4,header=T,sep="\t");
  
  
  ##-----only consider good barcode
  myindex = match(as.character(goodBarcode[,2]), colnames(geneRPKM));
  myindex2 = match(as.character(goodBarcode[,2]), as.character(seqDepth[,1]));
  geneRPKM = geneRPKM [,myindex];
  ##geneCount = geneCount [,myindex];
  read2GFP = as.numeric(seqDepth[myindex2,9]);
  read2genome = as.numeric(seqDepth[myindex2,2]);
  
  ##-----recording
  recSamplist = c(recSamplist,as.character(goodBarcode[,2]));
  if(length(recRPKM) == 0){
    recRPKM = geneRPKM;
    ##recCount = geneCount;
    recRead2GFP = read2GFP/read2genome;
  }
  else{
    recRPKM = cbind(recRPKM,geneRPKM);
    ##recCount = cbind(recCount,geneCount);
    recRead2GFP = c(recRead2GFP,read2GFP/read2genome);
  }
}	

geneList <- read.table(inFile1, header = T)[,1]
MTGenes <- grep("MT-", geneList)

for(i in c("3000", "2000")){
  ##----------------Use genes with the top variance to make heat map
  #-----------------Find most variable genes
  Variance <- apply(recRPKM, 1, var)
  Variance <- data.frame(geneList, Variance)
  Variance <- Variance[order(Variance[,2], decreasing = T),]
  no.MT.Variance <- Variance[-which(rownames(Variance)%in%MTGenes),]
  
  Top.Var <- no.MT.Variance[1:i,]
  Top.Var <- recRPKM[rownames(Top.Var),]
  Top.Var.geneList <- no.MT.Variance[,1]
  clust <- clusGap(Top.Var, FUN = kmeans, K.max = 20, B = 1000, iter.max = 30)
  save(clust, file = paste(Sys.Date(), "_clustGap_MT_genes_removed_Top_", i, "_var_genes.Rdata", sep = ""))
}

q()