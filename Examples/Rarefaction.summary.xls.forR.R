indext <- read.table(file="Rarefaction.summary.xls.forR", header=TRUE);
attach(indext)
#for otu
library(car)
library(multcomp)
ginfo <- indext$group

######Parametic test
#normally test qqplot
par(mfrow=c(2,2))
pdf("Rarefaction.summary.xls.forR.R.qqplot.pdf")
qqPlot(lm(sobs ~ ginfo), simulate=TRUE, main="Q-Q Plot 95% confidence envelope OTU(97%) number", labels=FALSE)
qqPlot(lm(chao ~ ginfo), simulate=TRUE, main="Q-Q Plot 95% confidence envelope Chao index", labels=FALSE)
qqPlot(lm(ace ~ ginfo), simulate=TRUE, main="Q-Q Plot 95% confidence envelope Ace index", labels=FALSE)
qqPlot(lm(shannon ~ ginfo), simulate=TRUE, main="Q-Q Plot 95% confidence envelope Shannon index", labels=FALSE)
dev.off()
## Equality of variance
##Bartlett test
a1bartlett <- bartlett.test(sobs,ginfo)
a1levene <- leveneTest(sobs, ginfo)
a1 <- aov(sobs ~ ginfo)
otui <- TukeyHSD(a1)

#for chao
a2bartlett <- bartlett.test(chao,ginfo)
a2levene <- leveneTest(chao, ginfo)
a2 <- aov(chao ~ ginfo)
chaoi <- TukeyHSD(a2)

#for ace
a3bartlett <- bartlett.test(ace,ginfo)
a3levene <- leveneTest(ace, ginfo)
a3 <- aov(ace ~ ginfo)
acei <- TukeyHSD(a3)

#for shannon
a4bartlett <- bartlett.test(shannon,ginfo)
a4levene <- leveneTest(ace, ginfo)
a4 <- aov(shannon ~ ginfo)
shannoni <- TukeyHSD(a4)

pdf("Rarefaction.summary.xls.forR.Tukey_HSD.test.otu.pdf")
#par(mar=c(5,4,6,2))
tes <- glht(a1, linfct=mcp(ginfo="Tukey"))
plot(cld(tes, level=.05), col="yellow")
dev.off()

pdf("Rarefaction.summary.xls.forR.Tukey_HSD.test.chao.pdf")
#par(mar=c(5,4,6,2))
tes <- glht(a2, linfct=mcp(ginfo="Tukey"))
plot(cld(tes, level=.05), col="yellow")
dev.off()

pdf("Rarefaction.summary.xls.forR.Tukey_HSD.test.ace.pdf")
#par(mar=c(5,4,6,2))
tes <- glht(a3, linfct=mcp(ginfo="Tukey"))
plot(cld(tes, level=.05), col="yellow")
dev.off()

pdf("Rarefaction.summary.xls.forR.Tukey_HSD.test.shannon.pdf")
#par(mar=c(5,4,6,2))
tes <- glht(a4, linfct=mcp(ginfo="Tukey"))
plot(cld(tes, level=.05), col="yellow")
dev.off()

##Following statistic data
a1bartlett

a1levene 

summary(a1)

otui

a2bartlett

a2levene

summary(a2)

chaoi

a3bartlett

a3levene

summary(a3)

acei

a4bartlett

a4levene

summary(a4)

shannoni
##End
