##-----Use RDAVIDWebservice to access DAVID

data <- neonatal.module.1                                         ##--Read a dataset

library("RDAVIDWebService")                                       ##--Load RDAVIDWebService 
user <- DAVIDWebService$new(email = "zhen.li.zl242@yale.edu")     ##--Log on to DAVID with email

getIdTypes(user) ##--Show gene id types

##--Get a list of genes, pay attention to idType
result <- addList(user, data, 
                  idType="ENSEMBL_GENE_ID", 
                  listType="Gene", 
                  listName = "neonatal.module.1")          

show(user)                                                         ##--Information about the run (or just use user)

getAllAnnotationCategoryNames(user)                                ##--Show annotation categories

setAnnotationCategories(user, "KEGG_PATHWAY")                      ##--check a specific category
geneCluster <- getClusterReport(user, type = "Gene") 
