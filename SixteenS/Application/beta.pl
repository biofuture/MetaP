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

die "perl $0 <otu.table.normalized.xls> <groups.xls> \n" unless(@ARGV==2);
##cluster analysis of OTUs table
cluster($ARGV[0], "ward");

##PCA of the OTUs table
ordination($ARGV[0],$ARGV[1]);
