##read the taxon table
mytable <- read.table(file="All.phylum.normalized_1", header=TRUE, row.names=1)
#install.packages(pkgs="PerformanceAnalytics")
library(PerformanceAnalytics)
pdf("All.phylum.normalized_1.barchart.phylum.pdf")
tm <- t(mytable)
chart.StackedBar(tm, date.format="%Y", cex.legend = 0.7, colorset=rainbow12equal,main="Barchart of phylum", ylab="Percent of Total", xlab="Sample ID", unstack=FALSE)
dev.off()
