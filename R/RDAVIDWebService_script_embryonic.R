emb <- read.csv("../30May2014 updates/embryonic_modules.csv", header = T)

library("RDAVIDWebService")                                       ##--Load RDAVIDWebService 
user <- DAVIDWebService$new(email = "zhen.li.zl242@yale.edu")     ##--Log on to DAVID with email

#usergetIdTypes(user)                                              ##--Show gene id types
#getAllAnnotationCategoryNames(user)                               ##--Show annotation categories

cbind.fill <- function(...){                                      ##--Custom cbind.fill function
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

####-----Get KEGG Pathway terms
KEGG <- c()

for (i in 1:ncol(emb)-1){
  data <- substr(emb[, i + 1], 1, 18)          ##--Read a gene list
  addList(user, data, idType = "ENSEMBL_GENE_ID", listType = "Gene",             ##--submit gene list to DAVID
          listName = paste("embryonic_module_", i, sep = ""))
  setAnnotationCategories(user, "KEGG_PATHWAY")                                  ##--check a specific category
  
  ##--combine GO terms from each list
  term <- getFunctionalAnnotationChart(user)
  if (i == 1){
    KEGG <- as.matrix(term$Term, ncol = 1)
  }
  else{
    if (is.null(term$Term)){
      KEGG <- cbind(KEGG, rep(x="NA", times=nrow(KEGG))
      )}
    else{
      KEGG <- cbind.fill(KEGG, as.matrix(term$Term, ncol=1))  
    }
  }
}
write.table(KEGG, "Embryonic_KEGG.txt", sep = "\t") 

user <- DAVIDWebService$new(email = "zhen.li.zl242@yale.edu")                    ##--Log on to DAVID with email
####-----Get GO_BP terms
GO_BP <- c()

##--Get a list of genes, pay attention to idType
for (i in 1:ncol(emb)-1){
  data <- substr(emb[, i + 1], 1, 18)           ##--Read a gene list
  addList(user, data, idType = "ENSEMBL_GENE_ID", listType = "Gene",             ##--submit gene list to DAVID
          listName = paste("embryonic_module_", i, sep = ""))
  setAnnotationCategories(user, "GOTERM_BP_FAT")                                  ##--check a specific category
  
  ##--combine GO terms from each list
  term <- getFunctionalAnnotationChart(user)
  if (i == 1){
    GO_BP <- as.matrix(term$Term, ncol = 1)
  }
  else{
    if (is.null(term$Term)){
      GO_BP <- cbind(GO_BP, rep(x="NA", each=nrow(GO_BP))
      )}
    else{
      GO_BP <- cbind.fill(GO_BP, as.matrix(term$Term, ncol=1))  
    }
  }
}

write.table(GO_BP, "Embryonic_GO_BP.txt", sep = "\t")
