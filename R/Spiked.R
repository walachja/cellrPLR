#Spiked apples dataload
library(BioMark)

library(BioMark)
data("spikedApples")



xx <- as.matrix(as.data.frame(spikedApples$dataMatrix))
yy <- as.factor(c(rep(1,10),rep(2,10)))

x <- xx
x[x<=0] <- 0.1
biomarker = paste0('x.',colnames(spikedApples$dataMatrix)[spikedApples$biom])
