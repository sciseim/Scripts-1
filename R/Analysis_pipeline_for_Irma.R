###--------Analysis pipeline for Irma
###--------This is after alignment and fragment assembly
###--------Goal is to put expressions of green cells, red cells and no-color cells into three different
###--------tables and do QC and analysis.
###--------Requires *_short_genes.fpkm_tracking file

###--------Step.1 Put reads from green cells into one table

batch = c('Irma-A12',
                'Irma-A3',
                'Irma-C12',
                'Irma-C2',
                'Irma-E10',
                'Irma-E4',
                'Irma-E6',
                'Irma-F2',
                'Irma-G2',
                'Irma-H3',
                'Irma-H9')

#----------Setup master file---------#
f <- read.table(paste("/media/zhen/Expansion Drive/fastq/Irma/green/", 
                      batch[1], "/", batch[1], "_genes.fpkm_tracking", sep = ""), header=T);

gene_short_name = f$gene_short_name

masterfile = data.frame(matrix(nrow = nrow(f), ncol = length(batch)))
rownames(masterfile) <- f[, which(colnames(f)=="gene_id")]
colnames(masterfile) <- batch

#----------collect FPKM from each genes.fpkm_tracking---------#
masterfile[, 1] <- f$FPKM

for(i in 2:length(batch)){
  f <- read.table(paste("/media/zhen/Expansion Drive/fastq/Irma/green/", 
                           batch[i], "/", batch[i], "_genes.fpkm_tracking", sep = ""), header=T);
  rownames(f) <- f[, which(colnames(f)=="gene_id")]
  masterfile[,i] <- f[rownames(masterfile), which(colnames(f)=="FPKM")]  
}

masterfile <- data.frame(gene_short_name, masterfile)

write.csv(masterfile, "/media/zhen/Expansion Drive/fastq/Irma/green/green_masterfile.csv") #------save masterfile

###---------Step.2 Put expression profile of red cells together

file <- read.table("/media/zhen/Expansion Drive/fastq/Irma/red/Irma-A7/Irma-A7_genes.fpkm_tracking", header=T)

batch = c(list.dirs(path = "/media/zhen/Expansion Drive/fastq/Irma/red/", full.names = F, recursive = F))

#----------Setup master file---------#
f <- read.table(paste("/media/zhen/Expansion Drive/fastq/Irma/red/", 
                      batch[1], "/", batch[1], "_genes.fpkm_tracking", sep = ""), header=T);

gene_short_name = f$gene_short_name

masterfile = data.frame(matrix(nrow = nrow(f), ncol = length(batch)))
rownames(masterfile) <- f[, which(colnames(f)=="gene_id")]
colnames(masterfile) <- batch

#----------collect FPKM from each genes.fpkm_tracking---------#
masterfile[,1] <- f$FPKM

for(i in 2:length(batch)){
  f <- read.table(paste("/media/zhen/Expansion Drive/fastq/Irma/red/", 
                        batch[i], "/", batch[i], "_genes.fpkm_tracking", sep = ""), header=T);
  rownames(f) <- f[, which(colnames(f)=="gene_id")]
  masterfile[,i] <- f[rownames(masterfile), which(colnames(f)=="FPKM")]  
}

masterfile <- data.frame(gene_short_name, masterfile)

write.csv(masterfile, "/media/zhen/Expansion Drive/fastq/Irma/red/red_masterfile.csv") #------save masterfile

###---------Step.3 Put expression profile of no-color cells together

file <- read.table("/media/zhen/Expansion Drive/fastq/Irma/no_color/Irma-A10/Irma-A10_genes.fpkm_tracking", header=T)

batch = c(list.dirs(path = "/media/zhen/Expansion Drive/fastq/Irma/no_color/", full.names = F, recursive = F))

#----------Setup master file---------#
f <- read.table(paste("/media/zhen/Expansion Drive/fastq/Irma/no_color/", 
                      batch[1], "/", batch[1], "_genes.fpkm_tracking", sep = ""), header=T);

gene_short_name = f$gene_short_name

masterfile = data.frame(matrix(nrow = nrow(f), ncol = length(batch)))
rownames(masterfile) <- f[, which(colnames(f)=="gene_id")]
colnames(masterfile) <- batch

#----------collect FPKM from each genes.fpkm_tracking---------#
masterfile[,1] <- f$FPKM

for(i in 2:length(batch)){
  f <- read.table(paste("/media/zhen/Expansion Drive/fastq/Irma/no_color/", 
                        batch[i], "/", batch[i], "_genes.fpkm_tracking", sep = ""), header=T);
  rownames(f) <- f[, which(colnames(f)=="gene_id")]
  masterfile[,i] <- f[rownames(masterfile), which(colnames(f)=="FPKM")]  
}

masterfile <- data.frame(gene_short_name, masterfile)

write.csv(masterfile, "/media/zhen/Expansion Drive/fastq/Irma/no_color/no_color_masterfile.csv") #------save masterfile

