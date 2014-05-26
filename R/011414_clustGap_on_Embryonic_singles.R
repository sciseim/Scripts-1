library("cluster")
recRPKM.Singles <- read.table("010914_recRPKM_NO_MT_Embryonic_Singles.txt", header = T)
clust <- clusGap(recRPKM.Singles, FUN = kmeans, K.max = 20, B = 1000, iter.max = 30)