#!/usr/bin/perl -w
use strict;


BEGIN {  
        unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}
 ##using pm file-------------------------------------------------------------------------------
use SixteenS::Tools::SequenceProcess;
use SixteenS::OtuTaxonTable;#otu_normalized1
use SixteenS::Tools::OtuFormatTransfer;
use SixteenS::Statistics::TaxonDistribution;
use SixteenS::Statistics::AlphaDiversity;
use SixteenS::Statistics::BetaDiversity;
use SixteenS::Tools::Matrix;#add_groupinfo

die "perl $0  <tsc.fa> <tsc.na> <tsc.otu> <name.list> <gast.list> " unless(@ARGV == 5);
my $MetaPotuf = "MetaP.otu.format";
my $qiimeotuf = "QIIME.otu.table";
my $taxontable = "Total.taxon.table";

my ($tscfa, $tscna, $tscotu, $namelist, $gast) = @ARGV;
my $date -= localtime;
transfer_tscotu($tscfa, $tscna, $tscotu, $MetaPotuf);
print "Generate OTU table Begin:", $date=localtime,"\n";
generate_otutable($namelist, $MetaPotuf, $gast,$qiimeotuf);
generate_taxontable($namelist, $gast, $taxontable );
print "OTU table and Taxon table generated:", $date=localtime,"\n";
##Generate old tax-otu table by tax-otu-new.pl
my $gastcsv = "rep.gast.csv";
my $mothurotu = "rep.mothur.format.otu";
gastscvmothurotu($gast,$MetaPotuf, $gastcsv, $mothurotu);
print "Start tax_otu\n";
my $taxotu = "Tax-otu-table";
tax_otu($namelist, $tscna, $mothurotu, $gastcsv,$taxotu);
print "END\n";
