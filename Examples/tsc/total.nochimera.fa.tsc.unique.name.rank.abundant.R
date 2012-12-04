#This is the R script for rank abundant curves
myvector <- read.table(file="rank.list")
len <- dim(myvector)[1]
x <- 1:len
y <- as.vector(myvector[,1])
pdf("tsc/total.nochimera.fa.tsc.unique.name.rank.abundance.pdf")
plot(x,log10(y),type="l", main="Rank abundance curve", xlab="Rank number", ylab="Relative abundance of unique tags (Log10)")
dev.off()
