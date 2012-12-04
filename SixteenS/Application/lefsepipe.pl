#!/usr/bin/perl -w 
use strict;


BEGIN {  
    use FindBin qw($Bin);
    my @dirset = split(/\//,$Bin);
    pop(@dirset); pop(@dirset);
    my $MetaPDir = join("/", @dirset);
    unshift @INC, "$MetaPDir";
}
use FindBin qw($Bin);
my @dirset = split(/\//,$Bin);
pop(@dirset); pop(@dirset);
my $MetaPDir = join("/", @dirset);

die "perl $0 <taxon.table>\n" unless (@ARGV==1);
print  "LEfSe begain:", my $date=localtime,"\n";
#LEfSe file 
my $taxontable = $ARGV[0];
`mkdir LEfSe; mv $taxontable LEfSe/;`;
 
my $lefsei = "LEfSe/$taxontable.addgroup";
my $lefsein = "$lefsei.in";
`python $MetaPDir/bin/LEfSe/format_input.py $lefsei $lefsein -u 2 -o 1000000`;
my $lefseres = "$lefsei.res";
`python $MetaPDir/bin/LEfSe/run_lefse.py $lefsein $lefseres --min_c 3 -l 2`;
`python $MetaPDir/bin/LEfSe/plot_res.py --format png $lefseres $lefseres.png`;
my $lefsepdf = "$lefsei.pdf";
`python $MetaPDir/bin/LEfSe/plot_cladogram.py --format pdf --dpi 300 $lefseres $lefsepdf`;
`mkdir LEfSe/biomarkers_raw_images`;
`python $MetaPDir/bin/LEfSe/plot_features.py $lefsein $lefseres LEfSe/biomarkers_raw_images/`;
print "LEfSe End: ", $date=localtime,"\n";
