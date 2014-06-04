##-----Use RDAVIDWebservice to access DAVID                                   

library("RDAVIDWebService")                                       ##--Load RDAVIDWebService 
user <- DAVIDWebService$new(email = "zhen.li.zl242@yale.edu")     ##--Log on to DAVID with email

getIdTypes(user)                                                  ##--Show gene id types
getAllAnnotationCategoryNames(user)                               ##--Show annotation categories

##--Get a list of genes, pay attention to idType
for (i in 2:40){
  data <- read.csv(paste("neonatal.module",i, ".csv"), header = T)[, 2] ##--Read a dataset
  setAnnotationCategories(user, "GOTERM_BP_FAT")                   ##--check a specific category
  getFunctionalAnnotationChartFile(user, "neonatal_GOterms.txt")#, paste("neonatal.module.",i,"_KEGG.csv"))
}
                  
                  
show(user)                                                        ##--Information about the run (or just use user)

