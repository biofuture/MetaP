#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::AlphaDiversity;
die "perl $0 <rarefaction.table> <classify group>\n" unless(@ARGV == 2);

box_rarefaction($ARGV[0], $ARGV[1]);

