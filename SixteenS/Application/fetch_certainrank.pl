#!/usr/bin/perl -w
use strict;

BEGIN {  
        unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::TaxonDistribution;
die "perl $0 <otu.table.Act> <rank> <output>\n" unless(@ARGV==3);
substract_sprank($ARGV[0], $ARGV[1], $ARGV[2]);
