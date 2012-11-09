#!/usr/bin/perl -w
use strict;
BEGIN {  
        unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::BetaDiversity;
die "perl $0 <matrix.table>\n" unless(@ARGV == 1);

ch_cluster_index($ARGV[0]);
