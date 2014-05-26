##----------------Find the top i expressed genes and make heatmap

dt = Sorensen.recRPKM
i = 200
max.dt <- dt[order(apply(dt, 1, max), decreasing = T),] 
max.dt[max.dt == 0] <- NA 
max.dt <- max.dt[order(max.dt[,1], max.dt[,2], max.dt[,3], decreasing = T), ] #reodering genes based on the expression in the "1" cell
#max.dt <- max.dt[,order(max.dt[1,],decreasing = T)] #reodering cells based on the expression of the "1" gene
max.dt[is.na(max.dt)] <- 0 
geneList <- read.table("brainSeq.singleCell.130925_SN7001318_0091_BC2NALACXX.wholeGene.geneComposite.gene.RPKM.txt", header = T)[,1] 
max.geneList <- geneList[as.numeric(rownames(max.dt))]
heatmap(as.matrix(max.dt)[10:190,], col=topo.colors(10000), main = paste("top", i, "genes", sep = " "),Rowv = T, Colv = NA, labRow = max.geneList[1:i], cexRow = 0.5, cexCol = 0.5, scale = "column")

##----------------Use genes with the top variance to make heat map

#-----------------Find most variable genes
Variance <- apply(recRPKM, 1, var)
Variance <- data.frame(rownames(recRPKM), Variance)
Variance <- Variance[order(Variance[,2], decreasing = T),]
genes <- read.table("~/../Dropbox/Single-cell/R code sharing/RPKM/brainSeq.singleCell.130925_SN7001318_0091_BC2NALACXX.wholeGene.geneComposite.gene.RPKM.txt", header = T)
MTGenes <- grep("MT-", genes[,1])
no.MT.Variance <- Variance[-MTGenes,]

#-----------------Generate Heat map with top i most variable genes
library("RColorBrewer")

dt = recRPKM
i = 100
Top.Var <- no.MT.Variance[1:i,]
Top.Var <- dt[rownames(Top.Var),]
Top.Var.geneList <- no.MT.Variance[,1]
#heatmap(as.matrix(Top.Var), col=topo.colors(1000), 
#        main = paste("top", i, "variable genes", sep = " "),
#        Rowv = T, Colv = T, labRow = NULL, #Top.Var.geneList, 
#        cexRow = 0.5, cexCol = 0.5, scale = "row")

pairs.breaks <- c(seq(-11, -0.4, length.out=50), seq(-0.4, 0.4, length.out =50),
                  seq(0.4, 11, length.out=50))
mycol <- colorpanel(n=149,low="blue", mid="white", high="red")

heatmap.2(as.matrix(Top.Var), breaks = pairs.breaks, col=mycol, 
          key=TRUE, keysize = 1, density.info = "histogram", symkey=TRUE, main="", 
          trace="none", cexRow=0.7, cexCol = 0.7, labRow = NULL, 
          Rowv = TRUE, Colv = TRUE, dendrogram = "column", 
          na.rm = TRUE, scale = "row"
          #lmat=rbind( c(0, 4, 3), c(1, 2, 0)), lhei=c(0.5, 5), lwid = c(0.5, 10, 4)
) # Col.Color = 1, Heat Map = 2, Key  = 3, Col.dendrogram = 4,