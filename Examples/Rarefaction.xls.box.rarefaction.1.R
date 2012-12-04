pdf("Rarefaction.xls.box.rarefaction.1.pdf")
matr <- read.table(file="Rarefaction.xls",sep="\t",header=T)
maxd <- 1
for(j in 2:dim(matr)[2])
{ 
    tem <- matr[,j]
        tem <- tem[!is.na(tem)]
        if(maxd < max(tem))
        {
            maxd <- max(tem)
        }
}
y <- matr[,2]
y <- y[!is.na(y)]
x <- matr[,1]
x <- x[1:length(y)]
#plot(x,y,ann=F,ylim=c(0,maxd),xlim=c(0,max(matr[,1])), xlab="Number of Sequences", ylab="Number of OTUS",main="Boxplot Rarefaction Curves For Groups")

selectmatr <- matr[c("AC1","AC2","AC3","AC4")]
coln <- dim(selectmatr)[2]
fac <- c() 
grp <- c()
#cl <- rgb(10, (0:(coln-1))/coln,127.5*1,max=255)
for(i in 1:coln){ 
y <- selectmatr[,i] 
y <- y[!is.na(y)]
x <- matr[,1]
x <- x[1:length(y)]
#points(x,y, col=palette()[1],type='l') 
fac <- c(fac,x)
grp <- c(grp,y)
}
if(1 == 1){
boxplot(grp~fac,col=palette()[1],boxwex=0.4,xlab="Number of Sequences", ylab="Number of OTUS",main="Boxplot Rarefaction Curves For Groups")
}else
{
    boxplot(grp~fac,col=palette()[1],add=TRUE,boxwex=0.4)
}
selectmatr <- matr[c("AM1","AM2","AM3","AM4")]
coln <- dim(selectmatr)[2]
fac <- c() 
grp <- c()
#cl <- rgb(20, (0:(coln-1))/coln,127.5*1,max=255)
for(i in 1:coln){ 
y <- selectmatr[,i] 
y <- y[!is.na(y)]
x <- matr[,1]
x <- x[1:length(y)]
#points(x,y, col=palette()[2],type='l') 
fac <- c(fac,x)
grp <- c(grp,y)
}
if(2 == 1){
boxplot(grp~fac,col=palette()[2],boxwex=0.4,xlab="Number of Sequences", ylab="Number of OTUS",main="Boxplot Rarefaction Curves For Groups")
}else
{
    boxplot(grp~fac,col=palette()[2],add=TRUE,boxwex=0.4)
}

legend("topleft",c("AC","AM"),pch=16, col=c(palette()[1:3]),title="Sample Group")
dev.off()