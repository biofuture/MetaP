#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Tools::Matrix;
die "perl $0 <matrix.table> <group.file> \n" unless(@ARGV == 2);

add_groupinfo($ARGV[0], $ARGV[1]);

