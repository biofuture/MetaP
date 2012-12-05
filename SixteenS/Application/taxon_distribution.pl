#!/usr/bin/perl -w 
use strict;


BEGIN {  
    use FindBin qw($Bin);
    my @dirset = split(/\//,$Bin);
    pop(@dirset); pop(@dirset);
    my $MetaPDir = join("/", @dirset);
    unshift @INC, "$MetaPDir";
}

use SixteenS::Statistics::TaxonDistribution;
die "perl $0 <genus.table>  <phylum.table>\n" unless (@ARGV==2);
barchart_taxon($ARGV[0], "Genus");
barchart_taxon($ARGV[1], "Phylm");
