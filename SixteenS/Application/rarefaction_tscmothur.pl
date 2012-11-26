#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::AlphaDiversity;
die "perl $0 <normalized.fa.list> <output dir> <tsc execution>\n" unless(@ARGV == 3);
tsc_mothurbatch($ARGV[0], $ARGV[1],$ARGV[2]);

