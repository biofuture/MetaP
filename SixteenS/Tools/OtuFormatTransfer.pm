#########################################
#10/10/2012
#Author JIANG Xiaotao
#Email biofuture.jiang@gmail.com
#
#Overall description
#This package is intended to transfer different format of OTUs results into a standard format
#for the process of MetaP library 
########################################
package SixteenS::Tools::OtuFormatTransfer;
use strict;
use SixteenS::Tools::SequenceProcess;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(transfer_tscotu tran_tsc2mothur);

sub transfer_tscotu{
    my @file = @_;
    die "Wrong number of imput files\n" unless($#file == 3);

    ##$file[0] the unique.fa 
    ##$file[1] the unique.name
    ##$file[2] the tsc otu results
    ##$file[3] standard otuformat output

    my %uniquefa2seq;
    my %uniqueindex;
    die "Not the correct file\n" unless(-e $file[0]);
    ##Fetch fasta information using function from SequencePorcess Package
    fetch_fastainfo($file[0],\%uniquefa2seq,\%uniqueindex); 

    die "Not the unique name file from tsc\n" unless (-e $file[1]) && open(I,$file[1]);
    my %uniquenumber;
    my %uniquecomp;
    #$/ = "\n";
    while(<I>)
    {
            chomp;
            my @tem = split /\s+/;
            $uniquenumber{$tem[0]} = $tem[1];
            $uniquecomp{$tem[0]} = $tem[2];
    }
    close I;
    die "Not the tsc Otu results \n" unless open(II,"$file[2]");
    die "$! \n" unless open(QIF, ">$file[3]");
    print QIF "MetaP Format\nOTU_ID\tRepresentTag\tRepresentSeq\tUniqueNumber\tTotalNumber\tComposition\n";
    while(<II>)
    {
        chomp;
        my ($otun,$uniqseqn,$ids) = split /\s+/;
        my @ids = split(",",$ids);
        my ($abna, $abnu,$comp,$tonu) = ("", 0, "",0); 
        
        for my $key (@ids)
        {
            if(exists $uniqueindex{$key})
            {
                if(exists $uniquecomp{$uniqueindex{$key}}){
                    $comp = join(",", $comp,$uniquecomp{$uniqueindex{$key}});
                }else{
                    die "wrong  1 $uniqueindex{$key}\n";
                }

                if(exists $uniquenumber{$uniqueindex{$key}})
                {
                    $tonu +=  $uniquenumber{$uniqueindex{$key}};
                    if($uniquenumber{$uniqueindex{$key}} > $abnu)
                    {
                        $abnu = $uniquenumber{$uniqueindex{$key}};
                        $abna = $uniqueindex{$key};
                    }
                }else
                {
                    die "wrong $uniqueindex{$key}\n";
                }
            }
            else
            {
                die "wrong $key\n";    
            }
        }
        $comp =~ s/^\,//;
        die "wrong $abna\n" unless(exists $uniquefa2seq{$abna});
        my $oneotu = join("\t",$otun,$abna,$uniquefa2seq{$abna},$uniqseqn,$tonu,$comp,);
        print QIF "$oneotu\n";
       # push @{$file[3]}, $oneotu; 
    }#while
}

sub tran_tsc2mothur{

##This programme is written for transferring OTU format form tsc to mothur
##The tsc input include the otu list file and the unique name file
##The mothur otu format could be used for rarefaction analysis
##Written by JXT 05/09/12 

    #die "perl $0 <tsc otu list> <tsc unique name>  <tsc.unique.fa> <output mothur format otu>\n" unless (@ARGV == 4);
    my @ar = @_;
    die "$! \n" unless open(I,$ar[0]);
    die "$! \n" unless open(II, $ar[1]);
    die "$! \n" unless open(III, $ar[2]);
    die "$! \n" unless open(T,">$ar[3]");

    my %rep;
    while(<II>)
    {
        chomp;
        my @tem = split /\s+/;
        $rep{$tem[0]} = $tem[2];
    }
    close II;

    my %seq;
    my $index = 0;
    while(<III>)
    {
        chomp;
        my @tem = split /\s+/;
        my $na = $tem[0];
        $na =~ s/^>//;
        $seq{$index} = $na;
        <III>;
        $index++;
    }
    close III;

    my $otunumber = 0;
    my @lastn=();
    while(<I>)
    {
        chomp;
        $otunumber++;
        my @tem = split /\s+/;
        my @num = split(/\,/, $tem[2]);
        my @tm = ();
        for my $ink (@num)
        {
            if(exists $seq{$ink})
            {
                if(exists $rep{$seq{$ink}})
                {
                    push @tm , $rep{$seq{$ink}};
                }
                else
                {
                    die "wrong occured \n";
                }
            }
            else
            {
                die "wrong occured ! \n";
            }
        }
        my $jo = join(",", @tm);
        push @lastn, $jo;
    }
    close I;
    my $ot = join("\t",@lastn);
    print  T "0.03\t$otunumber\t$ot\n";
  
}#tran_tsc2mothur


sub transfer_uclustotu{

}

sub transfer_espritotu{

}

sub transfer_dnaclustotu{

}


##############################################################
#Private sub routine 
##############################################################


1;
__END__
