#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::AlphaDiversity;
die "perl $0 <alpha.index> <groups file>\n" unless(@ARGV == 2);

anova_posthoc($ARGV[0], $ARGV[1]);

