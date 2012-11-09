#!/usr/bin/perl -w
use strict;
BEGIN {  
        unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::BetaDiversity;
die "perl $0 <matrix.table> <groups.file>\n" unless(@ARGV == 2);

ordination($ARGV[0], $ARGV[1]);
