#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::AlphaDiversity;
die "perl $0 <tsc.unique.name>\n" unless(@ARGV == 1);

rank_abundance($ARGV[0]);

