#!/usr/bin/perl -w
use strict;

##This programme is written for transferring OTU format form tsc to mothur
##The tsc input include the otu list file and the unique name file
##The mothur otu format could be used for rarefaction analysis
##Written by JXT 05/09/12 

die "perl $0 <tsc otu list> <tsc unique name>  <tsc.unique.fa> <output mothur format otu>\n" unless (@ARGV == 4);

die "$! \n" unless open(I,$ARGV[0]);
die "$! \n" unless open(II, $ARGV[1]);
die "$! \n" unless open(III, $ARGV[2]);
die "$! \n" unless open(T,">$ARGV[3]");

my %rep;
while(<II>)
{
	chomp;
	my @tem = split /\s+/;
	$rep{$tem[0]} = $tem[2];
}
close II;

my %seq;
my $index = 0;
while(<III>)
{
	chomp;
	my @tem = split /\s+/;
	my $na = $tem[0];
	$na =~ s/^>//;
	$seq{$index} = $na;
	<III>;
	$index++;
}
close III;

my $otunumber = 0;
my @lastn=();
while(<I>)
{
	chomp;
	$otunumber++;
	my @tem = split /\s+/;	
	my @num = split(/\,/, $tem[2]);
	my @tm = ();
	for my $ink (@num)
	{
		if(exists $seq{$ink})
		{
			if(exists $rep{$seq{$ink}})
			{
				push @tm , $rep{$seq{$ink}};
			}
			else
			{
				die "wrong occured \n";
			}
		}
		else
		{
			die "wrong occured ! \n";
		}
	}
	my $jo = join(",", @tm);
	push @lastn, $jo;
}
close I;
my $ot = join("\t",@lastn);
print  T "0.03\t$otunumber\t$ot\n";
