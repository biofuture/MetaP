#!/usr/bin/perl -w 
use strict;
use FindBin qw($Bin);

BEGIN {  
    my @dirset = split(/\//,$Bin);
    pop(@dirset); pop(@dirset);
    my $MetaPDir = join("/", @dirset);
    unshift @INC, "$MetaPDir";
}
use SixteenS::Tools::SequenceProcess;
use SixteenS::OtuTaxonTable;#otu_normalized1
use SixteenS::Tools::OtuFormatTransfer;
use SixteenS::Statistics::TaxonDistribution;
use SixteenS::Statistics::AlphaDiversity;
use SixteenS::Statistics::BetaDiversity;
use SixteenS::Tools::Matrix;#add_groupinfo

die "perl $0 <nochimera.fa.list> <normalized dir> <tsc excution> <group.xls>" unless(@ARGV==4); 
print "rarefaction begin", my $date=localtime,"\n";
my ($nochfalist, $normadir, $tsce, $opt_g) = @ARGV;
tsc_mothurbatch($nochfalist, $normadir, $tsce);
#merge generate rarefaction table and summary xls and table
my $rareflist = "rar.f.list";
my $raresumlist = "rar.sum.list";
`ls $normadir/tsc/*.rarefaction > $rareflist`;
`ls $normadir/tsc/*.summary > $raresumlist`;
my $rarefxls = "Rarefaction.xls";
my $raresumxls = "Rarefaction.summary.xls";
my $raresumtable = "Rarefaction.summary.table";
generate_raretable($rareflist, $rarefxls);
generate_alphaindex($raresumlist, $raresumxls, $raresumtable);
##test the normality and equality of variance
anova_posthoc($raresumxls, $opt_g);
##rarefaction curve plot
rarefaction_groups($rarefxls,$opt_g);
box_rarefaction($rarefxls, $opt_g);
print "\trarefaction analysis end: ",$date=localtime,"\n"; 

