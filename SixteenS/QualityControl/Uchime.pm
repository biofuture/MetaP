package SixteenS::QualityControl::Uchime;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(uchime);

use strict;
use threads;

##This package is mainly for Chmeria tags removing 
##The algorithm developed by Edgar Robert

sub uchime {
    
##The following items
#1.the first is the clean tag list
my @ar = @_;
die "$! \n" unless open(I, "$ar[0]");
unless (-d "nochimera")
{
    `mkdir nochimera`;
}
##start to process
while(<I>)
{
    #every line represent one sample
    chomp;
    my $dirf = $_;
    my $samplename = (split(/\//,$dirf))[-2];
    unless (-d "nochimera/$samplename"){
        `mkdir nochimera/$samplename`;    
    }
   #print "$samplename\n"; 
   unique_rank($dirf,"nochimera/$samplename");     
}
close I;

my $f = `ls nochimera/*/*.unique.fa`;
chomp($f);
my @temp = split(/\n+/, $f);
#my $round = int( ($#temp + 1) / $ar[1] )  + 1;

my @thre;
my @record;
for(my $i = 0; $i <= $#temp; $i++)
{
    $thre[$i] =  threads->create(\&uchime_single, "$temp[$i]");
}
for(my $i = 1; $i <= $#temp; $i++)
{
    $record[$i] = $thre[$i]->join();   
}


get_nochimeraseq($f);

}#uchime_process

sub uchime_single{

##The input is the unique fasta file which was 
## ranked by abundance
my @ar = @_;
`uchime --input $ar[0] --uchimeout $ar[0].out --uchimealns $ar[0].alns --minchunk 20 --xn 7 --noskipgaps2 --log $ar[0].log`;
}#uchime_single


sub unique_rank{
##this function unique fasta file 
##and rank the sequences
    my @ar = @_;
    open FILE1,$ar[0] || die "can not open FILE1 :$!";
    my %seq;
    my %fasta;
    my %name;
    while(<FILE1>)
    {
        chomp;
        if(/^>/)
        {
            my $seqna = substr $_,1;
            my $seq = <FILE1>;
            $fasta{$seq} = $seqna;
            if (!exists $seq{$seq})
            {
                $seq{$seq} = 1;
                $name{$seq} = $seqna;
            }
            else
            {
                $seq{$seq}++;
                $name{$seq} = "$name{$seq},$seqna"
            }
        }
    }
    close FILE1;

    my $pre = (split(/\//,$ar[0]))[-1];
    $pre =~ s/\.fa$//;
    #die "$pre\n";
    my $out = "$ar[1]/$pre.ranked.unique.fa";
    my $na = "$ar[1]/$pre.ranked.unique.names";
    open IM,">$out" || die "can not open I: $!";
    open II,">$na" || die "can not open II:$!";
    for my $keys (sort { $seq{$b} <=> $seq{$a} } keys %seq )
    {
        if ((exists $fasta{$keys}) && (exists $name{$keys}))
        {
            print IM ">$fasta{$keys}/ab=$seq{$keys}/\n$keys";
            print II "$fasta{$keys}\t$name{$keys}\n";
        }
    }
    close IM; close II;
}#unique_rank


sub get_nochimeraseq{

    #1 the rank unique fa

    my @ar = @_;
    my @temps = split(/\n+/,$ar[0]);
    foreach (@temps)
    {
        chomp; 
        my $unqfa = $_;
        my $prefix = $unqfa;
        $prefix =~ s/\.fa$//;
        my $unqna = "$prefix.names";
        my $uchout = "$unqfa.out";
        my @dirn = split(/\//, $unqfa);
        pop @dirn;
        my $di = join("/", @dirn);
        my $out = "$di/$dirn[-1].nochimera.fa";
##print "$rawfa\n$uchout\n$out\n";die;
        open II,"$unqfa" ||die "can not open II:$!";
        open III,"$uchout"|| die "can not open III:$!";
        open IV,"$unqna" || die "can not open IV:$!";
        open OUT,">$out" || die "can not open OUT:$!";
        my %chimera;
        while (<III>)
        {
            chomp;
            my @item = split/\s+/;
            my $result = "Y";
            if ($item[-1] eq $result)
            {
                my @na1 = split/\//,$item[1];
                $chimera{$na1[0]} = 1; 
            }
        }

        my %fasta;
        while(<II>)
        {
            chomp;
            if(/^>/)
            {
                my $name = substr$_,1;
                my @tnames = split/\//,$name;
                my $seq = <II>;
                $fasta{$tnames[0]} = $seq;
            }
        } 

        while(<IV>)
        {
            chomp;
            my @otun = split/\s+/;
            if(!exists $chimera{$otun[0]})
            {
                my @temp = split/\,/,$otun[1];
                for my $j(@temp)
                {
                    if (exists $fasta{$otun[0]})
                    {
                        print OUT ">$j\n$fasta{$otun[0]}";
                    }
                }
            }
        }
        %chimera=();
        %fasta=();
    }# foreach
}#get_nochimera seqs

1;
__END__
