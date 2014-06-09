##-----Use R to draw bar graph

library("plyr")
library("ggplot2")


data <- read.csv(file="/Users/Zhen/Dropbox/Single-cell/brainSeq.singleCell.chromRead.mouse.embyonic.csv", header = T)
chY <- data.frame(as.character(data$cellName), as.numeric(data$chrY), 
                  stringsAsFactors = FALSE)
for (i in c(1, 3, 6, 10)){
  num <- grep(paste("C1-", i, "-", sep = ""), chY$as.character.data.cellName)
  if (nchar(i) > 1){
    chY[num, 1] <- substr(chY$as.character.data.cellName[num], 1, 3+nchar(i))
  }
  else{
    chY[num, 1] <- substr(chY$as.character.data.cellName[num], 1, 4)
  }
}

chY <- data.frame(chY$as.character.data.cellName, as.numeric(chY$as.numeric.data.chrY))
colnames(chY) <- c("Animal", "RPKM_chY")

summary <- ddply(.data=chY, .variables="Animal", summarize,
                 N = length(RPKM_chY),
                 mean = mean(RPKM_chY),
                 sd = sd(RPKM_chY),
                 se = sd / sqrt(N))

plot <- ggplot(summary, aes(x=Animal, y=mean))
plot <- plot + geom_bar(position=position_dodge(), stat = "identity")
plot <- plot + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                     width=.2,
                     position=position_dodge(.9))
       
#qplot(x=chY$Animal, y=chY$RPKM_chY, geom="histogram", stat="identity") #----qplot, much like plot