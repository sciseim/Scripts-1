pc <- read.table("c:/Users/Zhen/Dropbox/principle_components", header = TRUE)
pc1 <- pc[,2]
pc2 <- pc[,3]
pc3 <- pc[,4]
plot(pc1, pc2)

minXaxis=min(pc1)*1.3;
maxXaxis=max(pc1)*1.3;
minYaxis=min(pc2)*1.3;
maxYaxis=max(pc2)*1.3;
minZaxis=min(pc3)*1.3;
maxZaxis=max(pc3)*1.3;

scatterplot3d(pc1, pc2, pc3, type = "p", highlight.3d = TRUE, xlim=c(minXaxis,maxXaxis),ylim=c(minYaxis,maxYaxis), zlim=c(minZaxis, maxZaxis),
xlab=paste("PC1 ","%",sep=""), ylab=paste("PC2 ","%",sep=""), zlab=paste("pc3","%", sep=""), pch=20)

text(pc1, pc2, pc3, labels=pc[,1],cex=1)