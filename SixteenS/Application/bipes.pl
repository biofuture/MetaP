#!/usr/bin/perl -w

use strict;
BEGIN {  
    unshift @INC, "/Users/jiangxiaotao/Documents/MetaP";
}

use SixteenS::QualityControl::Bipes;
#die "perl $0 <1.fq> <2.fq> <primer.list> <output>\n" unless(@ARGV == 4);

overlap_fq($ARGV[0], $ARGV[1], 6);

primer_cut($ARGV[2], "overlaped.total.pe.fa", 6, $ARGV[3]);
