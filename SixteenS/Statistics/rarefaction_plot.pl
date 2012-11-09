#!/usr/bin/perl -w
use strict;

die "perl $0 <total rarefaction data> <fastamap.list>\n" unless (@ARGV == 2);

die "$! \n" unless open(II,"$ARGV[1]");
<II>;
my $number = <II>;
chomp($number);
my @num = split(/\s+/, $number);
$number = $#num + 1;

my $i = 1;
my $ins = 255/$number;

while($i < $number)
{
	##For every catergorie samples
	die "$! \n" unless open(II,"$ARGV[1]");
	my %direct;
	while(<II>)
	{
		chomp;
		next if(/^#/);
		my @tem = split /\s+/;
		$direct{$tem[$i]}{$tem[0]} = 1;
	}
	close II;
	
	##generate R script
	my $Rout = "$ARGV[0].Rout";
	my $rareout = "$ARGV[0].rarefaction.$i.pdf";
	
	die "$! $Rout\n" unless open(TEM,">$Rout");
	print TEM  "pdf(\"$rareout\");\n matr <- read.table(file=\"$ARGV[0]\",sep=\"\\t\",header=T);\n maxd <- 1; \n for(j in 2:dim(matr)[2]){ \ntem <- matr[,j] \ntem <- tem[!is.na(tem)]; if(maxd < max(tem)){ maxd <- max(tem)} \n }  \n y <- matr[,2]; y <- y[!is.na(y)]; x <- matr[,1]; x <- x[1:length(y)]; \nplot(x,y,type='l',col='white',ylim=c(0,maxd),xlim=c(0,max(matr[,1])), xlab=\"Number of Sequences\", ylab=\"Number of OTUS\",main=\"Rarefaction Curves For Groups\")\n";
	
	my $cate = 1;	
	my @ts = ();
	for my $key (sort keys %direct)
	{
		my @tm =();
		for my $k(sort keys %{$direct{$key}})
		{
			push @tm, "\"$k\"";		
		}
		my $t = join(",",@tm);	
		my $cates = $cate * 10;

		push @ts, "\"$key\"";
		print TEM "selectmatr <- matr[c($t)]\n";
		my $str = <<STD ;
coln <- dim(selectmatr)[2]
#cl <- rgb($cates, (0:(coln-1))/coln,$ins*$i,max=255)
for(i in 1:coln){ 
y <- selectmatr[,i] 
y <- y[!is.na(y)]
x <- matr[,1]
x <- x[1:length(y)]
points(x,y, col=palette()[$cate],type='l') 
}	
STD
		print TEM $str;
		$cate++;						
	}	

	my $legtxt = join(",",@ts);
	print TEM "\nlegend(\"topleft\",c($legtxt),pch=16, col=c(palette()[1:$cate]),title=\"Sample Group\")\ndev.off()";

	close TEM;
	`/usr/local/bin/R CMD BATCH --args $Rout`;	

	%direct = ();
	$i++;
}
