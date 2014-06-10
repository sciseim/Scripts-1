##------keep only gene ID 
neonatal_class1 <- read.csv("060314_neonatal_class1.csv", header = T)

for (i in 1:ncol(neonatal_class1)){
  neonatal_class1[, i] <- substr(neonatal_class1[, i], 1, 19)
}

write.csv(neonatal_class1, "neonatal_class1.csv")


neonatal.modules <- read.csv("neonatal.modules.csv", header = T)
for (i in 1:ncol(neonatal.modules)){
  neonatal.modules[, i] <- substr(neonatal.modules[, i], 1, 18)
}

for (i in 1:ncol(neonatal.modules)){
  write.csv(neonatal.modules[,i], paste("neonatal.module", i-1, ".csv"))
}