#!/usr/bin/perl -w
use strict;

BEGIN {  
            unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::TaxonDistribution;
use SixteenS::OtuTaxonTable;#tax_otu

die "perl $0 <all gast> <MetaPotuf>  <namelist> <tscna> \n" unless(@ARGV==4);
my $gast = $ARGV[0];
my $MetaPotuf = $ARGV[1];
my $gastcsv = "gastcsv.txt";
my $mothurotu = "mothurotu.txt";

gastscvmothurotu($gast,$MetaPotuf, $gastcsv, $mothurotu);
my $namelist = $ARGV[2];
my $tscna = $ARGV[3];

tax_otu($namelist, $tscna, $mothurotu, $gastcsv,"Tax-otu-table");

