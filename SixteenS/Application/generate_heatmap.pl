#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::BetaDiversity;
die "perl $0 <matrix.table> <color.in>\n" unless(@ARGV == 2);

heatmap_cluster($ARGV[0], $ARGV[1]);

