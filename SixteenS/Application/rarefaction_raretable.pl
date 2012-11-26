#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Statistics::AlphaDiversity;
die "perl $0 <rare.f.list> <rare.table>\n" unless(@ARGV == 2);
generate_raretable($ARGV[0], $ARGV[1]);

