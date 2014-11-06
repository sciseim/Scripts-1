###------R script to get the number of genes above certain FPKM and save to read_depth.txt
###------Zhen Li

batch <- c("green", "red")
for(i in batch){
  rd <- read.csv(paste("/media/zhen/Expansion Drive/fastq/Irma/", i, "/", 
                         i, "_read_depth.csv", sep = ""), header = T)
  
  masterfile <- read.csv(paste("/media/zhen/Expansion Drive/fastq/Irma/", i, "/", 
                             i, "_masterfile.csv", sep = ""), header = T)
  e <- c()
  f <- c()
  for(x in 1:nrow(rd)){
    fpkm_0.1 <- length(which(masterfile[, x+1] > 0.1)) #-----count number of genes with FPKM > 0.1
    fpkm_1 <- length(which(masterfile[, x+1] > 1))     #-----count number of genes with FPKM > 1
    e <- c(e, fpkm_0.1)
    f <- c(f, fpkm_1)
  } 
  rd$"FPKM>0.1" <- e #-----add to read depth table
  rd$"FPKM>1" <- f   #-----add to read depth table
  rd$"Condition" <- rep(i, times = nrow(rd))
  
  write.csv(rd, paste("/media/zhen/Expansion Drive/fastq/Irma/", i, "/", 
                        i, "_read_depth.csv", sep = ""))    #-----save to read_depth.csv file
}



