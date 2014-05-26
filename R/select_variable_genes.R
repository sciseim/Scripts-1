##------------------Select most variable genes across cells
##-------------11/16/13

Variance <- apply(recRPKM, 1, var)
Variance <- data.frame(rownames(recCount), Variance)
Variance <- Variance[order(Variance[,2], decreasing = T),]
no.MT.Variance <- Variance[-which(rownames(Variance)%in%MTGenes),]