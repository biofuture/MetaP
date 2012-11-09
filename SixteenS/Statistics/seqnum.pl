#!/usr/bin/perl -w
use strict;
##################################################################
#To produce the basic sequence number information in BIPES process
##################################################################
if ($#ARGV ne 5)
{
 print "perl $0 <1.fq> <2.fq> <primer infor xls> <QC tag.list> <nochimera.fa.list> <outprefix>\n";
 exit;
}

open III,"$ARGV[2]" || die "can not open III:$!";
open IV,"$ARGV[3]" || die "can not open IV:$!";
open V,"$ARGV[4]" || die "can not open V:$!";

<III>;
$ARGV[2]=~s/.xls//;
`mkdir $ARGV[2]`;
my %fqnum;
my %spna;
my %bar;
while(<III>)
{
 chomp;
 my @temp = split/\s+/;
 my $fbar = substr $temp[2],0,8;
 my $rbar = substr $temp[4],0,5;
 $rbar =~ tr/atgc/ATGC/;
 $bar{$temp[-1]}{f}=$fbar;
 $bar{$temp[-1]}{r}=$rbar;
 ##print "$fbar\n$rbar\n";die;
 open I,"$ARGV[0]" || die "can not open I:$!";
 open II,"$ARGV[1]" || die "can not open II:$!";
 `mkdir $ARGV[2]/$temp[-1]`;
 my $fq1 = "$ARGV[2]\/$temp[-1]\/$temp[-1].1.fq";
 my $fq2 = "$ARGV[2]\/$temp[-1]\/$temp[-1].2.fq";
 open FQ1,">$fq1" || die "can not open FQ1:$!";
 open FQ2,">$fq2" || die "can not open FQ2:$!";
 my $num = 0;
 while(<I>)
 {
  my $na1 = $_;
  my $seq1 = <I>;
  my $fpos1 = index $seq1,$fbar;
  my $rpos1 = index $seq1,$rbar;
  my $na2 = <II>;
  my $seq2 = <II>;
  my $fpos2 = index $seq2,$fbar;
  my $rpos2 = index $seq2,$rbar;
  ##my $link1;
  ##my $link2;
  ##my $qul1;
  ##my $qul2;
  if ((($fpos1 <=2 && $fpos1 >=0) && ($rpos2 <=2 && $rpos2 >=0))||(($rpos1 <=2 && $rpos1 >=0)&&($fpos2 <=2 && $fpos2 >=0)))
  {
   my $link1 = <I>;
 ##print "$link1";die;
   my $link2 = <II>;
   my $qul1 = <I>;
   my $qul2 = <II>;  
   $num++;
  
   print FQ1 "$na1$seq1$link1$qul1";
   ##print  "$na1\n$seq1\n$link1\n$qul1";
   print FQ2 "$na2$seq2$link2$qul2";
  }
 }
 $fqnum{$temp[-1]} = $num; 
## print "$num";die;
 $spna{$temp[-1]} = $temp[0];
}

my $out = "$ARGV[5].seqnum.xls";

open OUT,">$out" || die "can not open OUT:$!";

my %tag;
while(<IV>)
{
 chomp;
 my @dir = split/\//;
 open TEMP,"$_" || die "can not open III:$!";
 my $num = 0;
 while (<TEMP>)
 {
  if(/^>/)
  {
   $num++;
  }
 }
 $tag{$dir[-2]}=$num; 
}

my %nochimera;
my %unqfa;
while(<V>)
{
 chomp;
 my @temp = split/\//;
 open T,"$_" || die "can not open V:$!";
 `/root/bin/454/unique_rank.uclust.pl $_`;
 s/\.fa//;
 my $unqfa = "$_.ranked.unique.fa";
 open UQ,"$unqfa" || die "can not open UQ:$!";
 my $num = 0;
 while (<T>)
 {
  if(/^>/)
  {
   $num++;
  }
 }
 $nochimera{$temp[-2]}=$num;
 
 my $a = 0;
 while(<UQ>)
 {
  if (/^>/)
  {
   $a++;
  }
  }
 $unqfa{$temp[-2]} = $a;
}

print OUT "SampleID\tforward_barcode\treverse_barcode\tRaw_seqnum\tQC_tag\tNochimera_seqnum\tunique\texperiment_samplename\n";
for my $i(sort keys %tag)
{
 print OUT "$i\t$bar{$i}{f}\t$bar{$i}{r}\t$fqnum{$i}\t$tag{$i}\t$nochimera{$i}\t$unqfa{$i}\t$spna{$i}\n";
}



