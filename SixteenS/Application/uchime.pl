#!/usr/bin/perl -w

use strict;
BEGIN {  
    use FindBin qw($Bin);
    my @dirset = split(/\//,$Bin);
    pop(@dirset); pop(@dirset);
    my $MetaPDir = join("/", @dirset);
    unshift @INC, "$MetaPDir";
}
use  SixteenS::QualityControl::Uchime;
die "perl $0 <fa.list>\n" unless (@ARGV == 1);
uchime($ARGV[0]);
