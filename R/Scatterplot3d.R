library("scatterplot3d")

##----------load embronic data
setwd("../Dropbox/Single-cell/Figures/21May2014 updates for 3d PCA/")
load("pcaSave.embryonic.RData");
rot = pcaSave;
clusterList = clusterList;

mergeCluster_pre = c("5","16","4_6","8_20","10_11_19","14_15","3_18_1_7","2_13_9");  
mycolor=rainbow(7) #c("blue","cyan","orange","black","seagreen","red","yellow","gray");


##----
color=c();
recIndex = c();
for(i in 1:length(mergeCluster_pre)){
  rec = unlist(strsplit(as.character(mergeCluster_pre[i]),split="_"));
  ##----
  recInd = c()
  new_color = c()
  for(j in 1:length(rec)){
    myindex = which(clusterList == rec[j]);
    recInd = c(recInd, myindex);
  }
  new_color = rep(mycolor[i], length(recInd));
  color = c(color, new_color)
  recIndex = c(recIndex, recInd)
}

##----set axis
minXaxis=-13000; #min(rot[,1])*1.3;
maxXaxis=8000; #max(rot[,1])*1.3;
minYaxis=-10000 #min(rot[,2])*1.3;
maxYaxis=max(rot[,2])*1.3;
minZaxis=-5000 #min(rot[,3])*1.3;
maxZaxis=4000



#scatterplot3d(rot[recIndex,1],rot[recIndex,2], rot[recIndex,3],
              #type="h",
              #highlight.3d = F,
              #color=color,
              #box = F,
              ##angle = 24,
              #xlim=c(minXaxis,maxXaxis),
              #ylim=c(minYaxis,maxYaxis), 
              #zlim=c(minZaxis,maxZaxis), 
              #xlab=paste("PC1", sep=""),
              #ylab=paste("PC2", sep=""),
              #zlab=paste("PC3", sep=""), 
              #pch=20);

for(i in 1:length(mergeCluster_pre)){
  rec = unlist(strsplit(as.character(mergeCluster_pre[i]),split="_"));
  recIndex = c()
  ##----
    for(j in 1:length(rec)){
    myindex = which(clusterList == rec[j]);
    recIndex = c(recIndex, myindex);
  }
  if(i == 1){
    s3d <- scatterplot3d(rot[recIndex,1:3],
                         type = "p",
                         pch = i + 14,
                         color = mycolor[i],
                         box = F,
                         grid = T,
                         xlim = c(minXaxis,maxXaxis),
                         ylim = c(minYaxis,maxYaxis), 
                         zlim = c(minZaxis,maxZaxis),
                         xlab = "PC1",
                         ylab = "PC2",
                         zlab = "PC3")
  }
  else{
    s3d$points3d(rot[recIndex,1:3],
                 col = mycolor[i],
                 pch = i + 14,
                 type = "p")
  }
}
s3d$plane3d(minXaxis, maxYaxis, maxZaxis)
