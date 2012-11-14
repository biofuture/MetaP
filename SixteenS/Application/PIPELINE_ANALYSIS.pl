#!/usr/bin/perl -w
use strict;

use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);

BEGIN {  
    my @dirset = split(/\//,$Bin);
    pop(@dirset); pop(@dirset);
    my $MetaPDir = join("/", @dirset);
    unshift @INC, "$MetaPDir";
}
use SixteenS::QualityControl::Uchime;
use SixteenS::QualityControl::Bipes;
use SixteenS::Tools::SequenceProcess;
use SixteenS::OtuTaxonTable;
use SixteenS::Tools::OtuFormatTransfer;
use SixteenS::Statistics::TaxonDistribution;

my @dirset = split(/\//,$Bin);
pop(@dirset); pop(@dirset);
our $MetaPDir = join("/", @dirset);

########################################
##-h verbose usage information for this script
########################################
my  $usage = <<USE;
    Author Xiao-tao JIANG
    Email: biofuture.jiang\@gmail.com
    Date 10/11/12
    perl $0 -a <1.fq> -b <2.fq> -c <primer.txt> -d <sample info> -n [threads number] -h [verbose] 
    
    -a  <1.fq>
    -b  <2,fq>
    -c  <primer file> the primer file contain the forward primer and the reverse primer
    -d  <sample info> format of sample info please refer to the example
    -n  [Threads number used when split file default 8]
    -h  print this help information
USE
#######################################
##BIPES quality control step for V6 data
##use Bipes.pm
#######################################
my $opt_log = "runpipe.log";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_n, $opt_h);
getopts('a:b:c:d:n:h');
die "Can not write $!\n" unless open(LOG, ">$opt_log");
my $date = localtime;
$opt_n ||= 8;
my $v6ref = "$MetaPDir/DB/refv6.fa";
my $v6rtax = "$MetaPDir/DB/refv6.tax";

if($opt_h || ($opt_a eq ""))
{
    die $usage;
}
print LOG "BIPES begain: ",$date=localtime,"\n";
##overlap the pair-end fastq reads
overlap_fq($opt_a, $opt_b, $opt_n);

primer_cut($opt_c,"overlaped.total.pe.fa", $opt_n, $opt_d);

print LOG "BIPES end: ", $date=localtime,"\n";

########################################
##chmiera sequences remove
##use UCHIME
########################################
`ls sample/*/*.fa > bipes.clean.fa.list`;
print LOG "UCHIME begain: ",$date = localtime,"\n";
uchime("bipes.clean.fa.list");
print LOG "UCHIME end: ",$date = localtime, "\n";

`ls nochimera/*/*.nochimera.fa > nochimera.fa.list`;


##rename of all the tags, generate the total clean fasta file
##generate name list file
rename_tags("nochimera.fa.list", "nochimerarename");

########################################
##Taxonomic assignment
########################################
`mkdir gast; cd gast; ln -s ../nochimerarename/total.nochimera.fa`;
print LOG "Gast Begin: ", $date = localtime, "\n";
`perl  $MetaPDir/SixteenS/TaxonomyAssignment/gast_files/gast  -in gast/total.nochimera.fa -ref $v6ref -rtax $v6rtax -out gast/total.gast.tax > gast.tax.log`;
print LOG "Gast End: ", $date= localtime, "\n";

########################################
##OTU picking steps
########################################
`mkdir tsc; cd tsc; ln -s ../nochimerarename/total.nochimera.fa`;
print LOG "TSC otu picking Begin: ", $date= localtime, "\n";
`$MetaPDir/bin/TSC_version_1.2_i86linux64  -i tsc/total.nochimera.fa -o tsc/total.nochimera.fa.tsc -n 8 -m al`;
print LOG "TSC otu picking End: ", $date= localtime, "\n";

########################################
##Generate Otu table and taxon table
########################################
#transfer tsc otu table into Qimme format for the convienent of the following analysis
my $tscotu = "tsc/total.nochimera.fa.tsc_0.03_otu.al.list";
my $tscfa = "tsc/total.nochimera.fa.tsc.unique.fa";
my $tscna = "tsc/total.nochimera.fa.tsc.unique.name";
my $MetaPotuf = "MetaP.otu.format";
my $gast = "gast/total.gast.tax";
print LOG "Transfer otu table format Begin:", $date=localtime,"\n";
transfer_tscotu($tscfa, $tscna, $tscotu, $MetaPotuf);
print LOG "Transfer otu table format End:", $date=localtime,"\n";


my $namelist = "nochimerarename/total.name.list";
my $qiimeotuf = "QIIME.otu.table";
my $taxontable = "Total.taxon.table";
print LOG "Generate OTU table Begin:", $date=localtime,"\n";
generate_otutable($namelist, $MetaPotuf, $gast,$qiimeotuf);

generate_taxontable($namelist, $gast, $taxontable );
print LOG "OTU table and Taxon table generated:", $date=localtime,"\n\n";

##Generate old tax-otu table by tax-otu-new.pl
my $gastcsv = "rep.gast.csv";
my $mothurotu = "rep.mothur.format.otu";
gastscvmothurotu($gast,$MetaPotuf, $gastcsv, $mothurotu);
print LOG "Start tax_otu\n";
tax_otu($namelist, $tscna, $mothurotu, $gastcsv,"Tax-otu-table");
 #test Tax-otu-table.Act.xls
 die "$gastcsv $\n" unless(-e "Tax-otu-table.Act.xls");
print "Tax otu end:\n";
##Generate otu table and certain rank taxon table from Tax-otu-table
my $otutable = "Tax-otu-table.Act.xls";
my $phylum = "All.phylum";
my $phylumnor1 = "All.phylum.normanized_1";
my $genus = "All.genus";
my $genusnor1 = "All.genus.normanized_1";
my $nor1otutable = "All.otu.sample.Act.table.normanized_1";
substract_sprank($otutable, 2, $phylum);
substract_sprank($otutable, 6, $genus);

print LOG "PIPE END:", $date=localtime,"\n";

##--------------------------------------------------------------------------------------------
#The following steps do some statistic analysis and graphic for upward generated result
#including 
#1. basic data summary
#2. alpha diversity analysis
#3. beta diversity analysis
#4. biomarker identifiying
##--------------------------------------------------------------------------------------------


exit();
#__END__
