#!/usr/bin/perl -w
use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP/";
}

use SixteenS::Tools::Matrix;
die "perl $0 <matrix.table> \n" unless(@ARGV == 1);

normalized($ARGV[0]);

