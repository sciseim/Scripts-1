###------R script to get the number of genes above certain FPKM and save to read_depth.txt
###------Zhen Li

rd <- read.table("/media/zhen/Expansion Drive/fastq/Irma/no_color/no_color_read_depth.txt", header = T)

e <- c()
f <- c()

for(i in 1:nrow(rd)){
  fpkm_0.1 <- length(which(masterfile[, i+1] > 0.1)) #-----count number of genes with FPKM > 0.1
  fpkm_1 <- length(which(masterfile[, i+1] > 1))     #-----count number of genes with FPKM > 1
  e <- c(e, fpkm_0.1)
  f <- c(f, fpkm_1)
} 

rd$"FPKM>0.1" <- e #-----add to read depth table
rd$"FPKM>1" <- f   #-----add to read depth table

write.table(rd, "/media/zhen/Expansion Drive/fastq/Irma/no_color/no_color_read_depth.txt", 
            sep = "\t") #-----save to read_depth.txt file
