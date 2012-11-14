#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::TaxonDistribution;
die "perl $0 <taxon file> <rank name>\n" unless(@ARGV == 2);

barchart_taxon($ARGV[0], $ARGV[1]);

