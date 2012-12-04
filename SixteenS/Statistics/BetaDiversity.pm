package SixteenS::Statistics::BetaDiversity;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(ordination cluster heatmap_cluster ch_cluster_index);

use SixteenS::Tools::Matrix;

##################################################
##Author JIANG Xiaotao
##Data 16/10/2012
##Ordination Plot   
##Different kinds of  distences caculation
##################################################

sub ordination{
    #this programe is to generate ordination plot in ecology 
    #using vegan ecology package of R
    #there are two arguments for this function, this first one is the MetaP otu table
    #, and the second one is groups table by which we plot the different groups with 
    #different color or shape
    #R must be installed in your computer correctly, and vegan package should be installed

    my @tm = @_;
    my $str = <<STDRS;
mytable <- read.table(file="$tm[0]", header=TRUE,row.names=1)
mygroup <- read.table(file="$tm[1]")
library(vegan)
library(ggplot2)
myt <- t(mytable)
Categorie <- mygroup[,2]
pdf("$tm[0].pca12.pdf")
pr_pca <- rda(myt)
perc <- pr_pca\$CA\$eig / pr_pca\$CA\$tot.chi
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

pdf("$tm[0].pc1pc2.ade4.center.pdf")
s.class(myt.pca1\$li, fac, cpoint=0, clabel=0.8)
dev.off()
pdf("$tm[0].pc1pc2.ade4.centerfalse.pdf")
s.class(myt.pca2\$li, fac, cpoint=0, clabel=0.8)
dev.off()
pdf("$tm[0].pc1pc2.ade4.bac2.pdf")
s.class(bet2\$ls, fac, csub = 1.75)
dev.off()
pdf("$tm[0].pc1pc2.ade4.bca1.pdf")
s.class(bet1\$ls, fac, csub = 1.75)
rand1\$pvalue
rand2\$pvalue
dev.off()

exit()
STDRS
    my $rs = "$tm[0].R";
    die "$rs $!\n" unless open(T, ">$rs");
    print T "$str"; close T;

    `R CMD BATCH $rs`;
}##ordination

sub cluster{ 
    #this programe is to process matrix data 
    #and cluster samples by matrix data 
    #there are several distance caculation methods such as vegdist (derived from vegan package), euclidean 
    #and several clusters methods (complete , single , median , ward etc)
    #the imput matrix is MetaP format otu table
    #You can choose whether to normalized the table or not by par 3
    #plot were finished by R, so you should make sure R installed correct 
    #and several R packages should be installed which are vegan ggplot2 stat

    #input of the function
    #this first parameter was the MetaP otu table
    #the second parameter was the groups file 
    #the third parameter define the normalized state
        
    my @tm = @_;
    die "Input parameters are wrong $!\n" unless (-e $tm[0]);
    my $method = "$tm[1]";
    my $str = "";
    $str =  <<STDIN;
#This is the R scripts
mytable <- read.table(file="$tm[0]",header=TRUE,row.names=1)
pdim <- dim(mytable)
nr <- pdim[1]
nc <- pdim[2]
library(vegan)
tmy <- t(mytable)
mdist <- vegdist(tmy)
clustF <- hclust(mdist, method="$tm[1]")
pdf(file="$tm[0].clust.$tm[1].pdf")
plot(clustF,hang=-1,main="Hierachical clustering tree")
dev.off()
STDIN
    my $or = "$tm[0].R";
    die "Output R file$!\n" unless open(T, ">$or");
    print T "$str"; close T;
    `R CMD BATCH --args $or`;
}#cluster

sub pam_cluster{
}#pam_cluster

sub kmeans_cluster{
}#kmeans_cluster

sub ch_cluster_index{
    #This function is used to caculate how many clusters
    #there are for dataset
    #this function will generate a list of value for cluster numbers
    #and CH Calinski-Harabasz (CH) index values
    #this validation using PAM clust data into number of cluster first and 
    #then caculate the CH value
    #the library clusterSim should be installed
    my @ar = @_;
    my $or = "$ar[0].R";
    die "$! \n" unless open(T,">$or");
    print T <<STD;
#This is the Rscript for Calinski-Harabasz index
mydata <- read.table(file="$ar[0]", header=TRUE, sep="\\t",row.names=1)
library(clusterSim)
rank <- c(2:10)
chindex <-rep(0, 9)
test <- mydata[1:50000,]
for(i in 2:10)
{
    nu <- pam(test, i);
    chi <- index.G1(test, nu\$clustering);
    chindex[(i-2)] <- chi
}
p <- max(chindex) 
position <- which(chindex==p)
pdf("$ar[0].chindex.pdf")
plot(rank, chindex, main="CH index ~ Number of clusters", xlab="Number of clusters", ylab="CH index", type="l")
abline(v=rank[position], col="yellow", lwd=6)
dev.off()
STD
close T;
    `R CMD BATCH --args $or`;
}#ch_cluster_index

sub heatmap_cluster{
    
    ##This is to plot heatmap and alternative cluster figure
    ##To run this program you should have connection to the internet or
    ##pre install R gplots package
    
    my @ar = @_;
    my $ors = "$ar[0].R";
    die "$ors $! \n" unless open(T,">$ors");
    print T <<STD;
#This is the R script for heatmap
mydata <- read.table(file="$ar[0]", header=TRUE, row.names=1)
mycolor <- read.table(file="$ar[1]")
pcolor <- mycolor[,2]
pcolor <- as.character(pcolor)
library("gplots")
pdf("$ar[0].heatmap.pdf")
heatmap.2(mydata, scale="row",col=redgreen(75), ColSideColors=pcolor, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()
STD
close T;
    `R CMD BATCH --args $ors`;
}#heatmap


sub lefse_plot{


}#lefse_plot

sub unifrac_plot{


}#unifrac_plot


1;
__END__
