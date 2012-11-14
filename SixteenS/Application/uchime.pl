#!/usr/bin/perl -w

use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP";
}

use  SixteenS::QualityControl::Uchime;
uchime("$ARGV[0]");
