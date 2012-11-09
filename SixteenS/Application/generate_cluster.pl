#!/usr/bin/perl -w
use strict;

BEGIN {  
        unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::BetaDiversity;
die "perl $0 <matrix for cluster> <method used>\n" unless(@ARGV==2);

cluster($ARGV[0], $ARGV[1]);
