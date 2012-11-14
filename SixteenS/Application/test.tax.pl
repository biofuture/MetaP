#!/usr/bin/perl -w
use strict;
BEGIN {
        unshift @INC, "/Users/jiangxiaotao/Documents/MetaP";
}

use SixteenS::OtuTaxonTable;
die "perl $0 [Name.list] [uniq.name] [otu.seq.name][otu.tax][out.tablei.prefix]\n" unless(@ARGV == 5);
tax_otu($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
