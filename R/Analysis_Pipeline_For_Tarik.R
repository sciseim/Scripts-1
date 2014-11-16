###----My analysis script for collaboration with Tarik
###----Single cell analysis of VZP 
###----Zhen Li
###----11/14/2014

library("limma");
library("DESeq");
library("cqn");

setwd("/Users/Zhen/Box Sync/Manuscript/Tarik/")

#-----Use custom function multmerge to merge files
multmerge = function(mypath, pattern){
  filenames=list.files(path=mypath, pattern=pattern)
  datalist = lapply(filenames, function(x){read.table(file=x,header=T)})
  Reduce(function(x,y) {rbind(x[,1:9],y[,1:9])}, datalist)}

countMapping <- multmerge(mypath = getwd(), pattern = ".countMapping.txt")

samples <- countMapping[grep("MMB", countMapping$Barcode),]

write.csv(samples, "master_countmapping.csv")
#After this step, I manually put in the cell number of each sample.

####----------------------------only analyze single cells---------------------------### 

#-----Use custom function multmerge to merge files
multmerge = function(mypath, pattern){
  filenames=list.files(path=mypath, pattern=pattern)
  datalist = lapply(filenames, function(x){read.table(file=x,header=T)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)}

gene_count <- multmerge(mypath = getwd(), pattern = ".gene.count.txt")
RPKM_count <- multmerge(mypath = getwd(), pattern = ".gene.RPKM.txt")

#-----Find good single cells
master_countmapping <- read.csv("master_countmapping.csv", header = T)
good_singles <- master_countmapping$Barcode[which(master_countmapping$Cell_number==1 
                                  & master_countmapping$Genome > 100000)]

#-----Find genes which are not detected in any cell
a <- rowSums(gene_count[,-which(names(gene_count) %in% c("Geneid", "row.names")])

non_0_count <- gene_count[-which(a==0), good_singles]
non_0_RPKM <- RPKM_count[-which(a==0), good_singles]

non_0_count <- non_0_count[,-1]
non_0_RPKM <- non_0_RPKM[, -1]
#-----PCA based on non-zero genes
##----------do pca analysis
pca=prcomp(t(non_0_RPKM[ ,grep("MMB", names(non_0_RPKM))]));
plot(pca,xlab="Principle components");

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

col <- c(rep(x = colors()[30], times = length(grep(pattern = "MMB13", 
                                     names(non_0_RPKM)))),
         rep(x = colors()[150], times = length(grep(pattern = "MMB26", 
                                            names(non_0_RPKM)))))

plot(rot[,1],rot[,2],type="p",xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis),
     xlab=paste("PC1 ",ratio1,"%",sep=""),ylab=paste("PC2 ",ratio2,"%",sep=""),pch=20,
     col= col, 
     #cex=GFPcex,
     cex.lab=1, cex.axis=1);

text(rot[,1],rot[,2],labels=recSamplist,cex=0.5);

#------zsGreen vs mCherry
only_red <- master_countmapping$Barcode[which(master_countmapping$mCherry > 4 
                                              & master_countmapping$zsGreen < 4)] 
              #----At least 4 reads for mCherry but less than 4 reads for zsGreen
only_green <- master_countmapping$Barcode[which(master_countmapping$mCherry < 4 
                                                & master_countmapping$zsGreen > 4)]
             #-----At least 4 reads for zsGreen but less than 4 read for mCherry

good_red_singles <- good_singles[which(good_singles %in% only_red)]
good_green_singles <- good_singles[which(good_singles %in% only_green)]

non_0_red_RPKM <- non_0_RPKM[,which(names(non_0_RPKM) %in% good_red_singles)]
non_0_green_RPKM <- non_0_RPKM[,which(names(non_0_RPKM) %in% good_green_singles)]

#-----take a look at EOMES, PAX6, SOX2
g <- lapply(c("ENSMUSG00000032446|EOMES", "ENSMUSG00000027168|PAX6", "ENSMUSG00000074637|SOX2"), FUN = function(x){grep(pattern = x, RPKM_count$Geneid, fixed = T)})

non_0_red_RPKM[which(rownames(non_0_red_RPKM)%in%g),]
non_0_green_RPKM[which(rownames(non_0_green_RPKM)%in%g),]


