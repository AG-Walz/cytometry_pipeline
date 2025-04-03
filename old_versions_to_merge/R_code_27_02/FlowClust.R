library(flowClust)

cols_of_interest <- colnames(cc_ff)[9:41]
fC <-flowClust::flowClust(x=cc_ff[9:41], varNames = cols_of_interest, K = 25, B=100)

summary(fC)
plot(fC, data=cc_ff[9:41])
hist(fC, data=cc_ff[9:41])
