####---------------------analysis on single cells only-----------------------------####
####-------------------------------11/27/13----------------------------------------####

Single_cells <- recRPKM[, -which(colnames(recRPKM)%in%
                                   c("MMB3_3000",
                                     "MMB3_Cells", 
                                     "MMB3_1000", 
                                     "MMB3_100",
                                     "MMB6_E6"
                                     "MMB6_H9"))]

##-----1------do  quantile normalization

library("limma")
NORM_Single = normalizeBetweenArrays(as.matrix(Single_cells),method="quantile");
#boxplot(NORM_Single,las=2,cex=0.5,main="Quantile normalization");


##-----2------do cluster analysis

mydist=as.dist(1-cor(NORM_Single,method="s"));
hca = hclust(mydist,method="complete"); 
plot(hca,axes=F,xlab="Gene Expression",ylab="",cex=0.5, 
     main = "");


##------3-----Weighted Gene Correlation Network Analysis

library("WGCNA")

# Choose a set of soft-thresholding powers
powers = c(seq(from = 1, to = 10, by = 0.5), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(t(NORM_Single), powerVector = powers, verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.95,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net = blockwiseModules(t(NORM_Single), corType = "bicor", power = 5,
                       networkType = "signed", TOMType = "unsigned", 
                       minModuleSize = 30, reassignThreshold = 0, 
                       mergeCutHeight = 0.25, numericLabels = TRUE, 
                       pamRespectsDendro = FALSE, saveTOMs = TRUE,
                       saveTOMFileBase = "SingleCell", verbose = 3)

table(net$colors)

mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

mod.den <- as.dendrogram(net$dendrograms[[1]])


##--------------Use module cluster to draw heatmap
heatmap.2(as.matrix(NORM_Single), #breaks = pairs.breaks, col=mycol, 
          key=TRUE, keysize = 1, density.info = "histogram", symkey=F, main="", 
          trace="none", cexRow=0.7, cexCol = 0.7, labRow = "", 
          Rowv = mod.den, Colv = T, 
          dendrogram = "column", na.rm = TRUE, 
          scale = "none", RowSideColors = mergedColors[net$blockGenes[[1]]]
)
