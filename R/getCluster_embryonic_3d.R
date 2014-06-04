

##Usage: Rscript getCluster_embryonic.R
##AIM: combine different runs for pca analysis
##Date:2013-10-14
##Modified:2014-05-15

library(rgl);

load("pcaSave.embryonic.RData");
rot = pcaSave;
clusterList = clusterList;

##--set axis
x1 = min(rot[,1])*1.3;
x2 = max(rot[,1])*1.3;
y1 = min(rot[,2])*1.3;
y2 = max(rot[,2])*1.3;
z1 = min(rot[,3])*1.3;
z2 = max(rot[,3])*1.3;

####------re-plot pca##--plot
mergeCluster_pre = c("5","16","4_6","8_20","10_11_19","14_15","3_18_1_7","2_13_9");	
mycolor=c("blue","cyan","orange","black","seagreen","red","yellow","gray");

for(i in 1:length(mergeCluster_pre)){
    rec = unlist(strsplit(as.character(mergeCluster_pre[i]),split="_"));
    ##----
    recIndex = c();
    for(j in 1:length(rec)){
        myindex = which(clusterList == rec[j]);
        recIndex = c(recIndex, myindex);
    }
	if(i==1){
		plot3d(rot[recIndex,1:3],col=mycolor[i],xlim=c(x1,x2),ylim=c(y1,y2),zlim=c(z1,z2),size=5,type="s",xlab="PC1",ylab="PC2",zlab="PC3")
	}
	else{
		plot3d(rot[recIndex,1:3],col=mycolor[i],add=TRUE,size=5, type="s")
	}
}




