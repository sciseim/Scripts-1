##-----Use R to draw bar graph

library("plyr")
library("ggplot2")


data <- read.csv(file="/Users/Zhen/Dropbox/Single-cell/brainSeq.singleCell.chromRead.mouse.embyonic.csv", header = T)
chY <- data.frame(as.character(data$cellName), as.numeric(data$chrY), 
                  stringsAsFactors = FALSE)
for (i in c(1, 3, 6, 10)){
  num <- grep(paste("C1-", i, "-", sep = ""), chY$as.character.data.cellName)
  if (i == 10){
    chY[num, 1] <- substr(chY$as.character.data.cellName[num], 1, 5)
  }
  else{
    chY[num, 1] <- substr(chY$as.character.data.cellName[num], 1, 4)
  }
}

chY <- data.frame(chY$as.character.data.cellName, as.numeric(chY$as.numeric.data.chrY))
colnames(chY) <- c("Animal", "RPKM_chY")

qplot(x=chY$Animal, y=chY$RPKM_chY, geom="histogram", stat="identity")

chrY <- c(mean(chrY1), mean(chrY3), mean(chrY6), mean(chrY10))

ggplot()
barplot2(geom_bar)


apply(chY, grep, )