#!/usr/bin/perl -w 
use strict;


BEGIN {  
    use FindBin qw($Bin);
    my @dirset = split(/\//,$Bin);
    pop(@dirset); pop(@dirset);
    my $MetaPDir = join("/", @dirset);
    unshift @INC, "$MetaPDir";
}

use SixteenS::Statistics::ReadsSummary;
die "perl $0 <1.fq> <2.fq> <primer.info.xls> <QC fa .list> <nochimera.fa list> <output prefix> \n" unless(@ARGV==6);

get_summary(@ARGV);
