mytable <- read.table(file="All.otu.table.act.xls.normalized_1", header=TRUE,row.names=1)
mygroup <- read.table(file="group.xls")
library(vegan)
library(ggplot2)
myt <- t(mytable)
Categorie <- mygroup[,2]
pdf("All.otu.table.act.xls.normalized_1.pca12.pdf")
pr_pca <- rda(myt)
perc <- pr_pca$CA$eig / pr_pca$CA$tot.chi
p1p2score <- scores(pr_pca, choices = c(1,2), dis="sites")
withcaterp1p2 <- data.frame(p1p2score,Categorie)
colnames(withcaterp1p2)[3] <- "Categories"
xl <- paste("PC1(",(perc[1]*10000)%/%100,"%)",sep="")
yl <- paste("PC2(", (perc[2]*10000)%/%100,"%)",sep="")
qplot(PC1,PC2,data=withcaterp1p2,xlab=xl,ylab=yl,colour=Categorie,size=I(6))
dev.off()

library(ade4)
fac <- as.factor(Categorie)
lev <- levels(fac)
myt.pca1 <- dudi.pca(myt, center=TRUE, scale=FALSE, scan=FALSE)
bet1 <- bca(myt.pca1, fac, scan=FALSE, nf=2)
myt.pca2 <- dudi.pca(myt, center=FALSE, scale=FALSE, scan=FALSE)
bet2 <- bca(myt.pca2, fac, scan=FALSE, nf=2)
rand1 <- randtest(between(myt.pca1, fac, scan = FALSE), 999)
rand2 <- randtest(between(myt.pca2, fac, scan = FALSE), 999)

pdf("All.otu.table.act.xls.normalized_1.pc1pc2.ade4.center.pdf")
s.class(myt.pca1$li, fac, cpoint=0, clabel=0.8)
dev.off()
pdf("All.otu.table.act.xls.normalized_1.pc1pc2.ade4.centerfalse.pdf")
s.class(myt.pca2$li, fac, cpoint=0, clabel=0.8)
dev.off()
pdf("All.otu.table.act.xls.normalized_1.pc1pc2.ade4.bac2.pdf")
s.class(bet2$ls, fac, csub = 1.75)
dev.off()
pdf("All.otu.table.act.xls.normalized_1.pc1pc2.ade4.bca1.pdf")
s.class(bet1$ls, fac, csub = 1.75)
rand1$pvalue
rand2$pvalue
dev.off()

exit()
