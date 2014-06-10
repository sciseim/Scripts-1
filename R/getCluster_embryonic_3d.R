

##Usage: Rscript getCluster_embryonic.R
##AIM: combine different runs for pca analysis
##Date:2013-10-14
##Modified:2014-05-15

library(rgl);

load("../Dropbox/Single-cell/Figures/21May2014 updates for 3d PCA/pcaSave.embryonic.RData");
rot = pcaSave;
clusterList = clusterList;

##--set axis
minXaxis=-13000; #min(rot[,1])*1.3;
maxXaxis=8000; #max(rot[,1])*1.3;
minYaxis=-10000 #min(rot[,2])*1.3;
maxYaxis=max(rot[,2])*1.3;
minZaxis=-5000 #min(rot[,3])*1.3;
maxZaxis=4000

####------re-plot pca##--plot
mergeCluster_pre = c("5","16","4_6","8_20","10_11_19","14_15","3_18_1_7","2_13_9");	
mycolor=c("blue","cyan","orange","black","seagreen","red","yellow","gray");

for(i in 1:length(mergeCluster_pre)){
  ##----split cluster into individual subcluster
  rec = unlist(strsplit(as.character(mergeCluster_pre[i]),split="_"));
    
  recIndex = c();
  for(j in 1:length(rec)){
      myindex = which(clusterList == rec[j]);
      recIndex = c(recIndex, myindex);
  }
	if(i==1){
		plot3d(rot[recIndex,1:3],
           col = mycolor[i],
           box = F,
           axes = T,
           cex = 0.25,
		       xlab="PC1",ylab="PC2",zlab = "PC3",
           xlim = c(minXaxis,maxXaxis),
           ylim = c(minYaxis,maxYaxis),
           zlim = c(minZaxis,maxZaxis),
           size = 2,type="s",
           axis = F
           )
	}
	else{
		plot3d(rot[recIndex,1:3],col=mycolor[i],
           add=TRUE, size=2, 
           type="s"
           )
	}
}

decorate3d(xlab = "PC1",ylab="PC2",zlab="PC3", axes = T, 
           box = F)
rgl.bg(fogtype = "linear", sphere=F, 
       color=c("white","green"), lit=FALSE, 
       back="lines" )
rgl.postscript()
