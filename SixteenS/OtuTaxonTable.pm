###################################
#Date 09/10/2012 
#Author Jiang Xiaotao
#contact Email biofuture.jiang@gmail.com
#This package is used to merge the taxon information and the OTU information 
#to generate a OTU table and taxon table for all samples 
###################################

package SixteenS::OtuTaxonTable;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(generate_otutable generate_taxontable);

use strict;
sub generate_otutable{
    ##make sure there have 
    ##1. sample reads mapping file
    ##2. OTU picking results table which was transfored by TSC/ESPRIT/UCLUST 
    ##3. Taxonomy classification of all tags form RDP/GAST   
    ##4. must be a array reference
    my @ar = @_;
    die "Not the sample tags mapping file\n" unless open(I,"$ar[0]");
    ##generate all the samples and stored the mapping 
    my %tag2sample;
    my %sample;
    while(<I>)
    {
        chomp;
        my @tem = split /\s+/;
        $sample{$tem[0]}++;
        $tag2sample{$tem[1]} = $tem[0];
    }## I
    close I;
    my $sindex = 1;
    my %sampleline;
    my $samplelist="";
    my @samplelist;
    for my $key(sort keys %sample)
    {
        $sampleline{$key} = $sindex;
        $sindex++;
        push @samplelist, $key;
    }
    $samplelist = join("\t",@samplelist);
    die "Not the formated taxonomy table\n" unless open(III,"$ar[2]");
    my %tag2tax;
    <III>;
    while(<III>)
    {
        chomp;
        my @tem = split /\s+/;
        $tag2tax{$tem[0]} = $tem[1];
    }
    close III;

    die "Not the MetaP standard otu table\n" unless open(II,$ar[1]);
    die "Output to the destination file\n" unless open(OUT, ">$ar[3]");
    <II>;<II>;
    my $log = "MetaP OTU table Version 1.0"; 
    #push @{$ar[3]}, $log;
    print OUT "$log\n";
    my $mark = "OTUID\t$samplelist\tTaxonomy";
    #push @{$ar[3]}, $mark;
    print OUT "$mark\n";
    #print "$log\n$mark\n";
    my $snumber = $#samplelist+1;
    while(<II>)
    {
        chomp;
        my @tem = split /\s+/;
        my @temp = split(",",$tem[5]);
        for(my $i= 1; $i <= $snumber; $i++)
        { 
            $samplelist[$i] = 0;
        }
        for my $ke (@temp)
        {
            if(exists $tag2sample{$ke} &&  exists $sampleline{$tag2sample{$ke}})
            {
                $samplelist[$sampleline{$tag2sample{$ke}}]++; 
            }else
            {
                die "Wrong $ke\n";    
            }
        }
        shift @samplelist;
        my $number = join("\t",@samplelist);

        die "Wrong tax\n" unless (exists $tag2tax{$tem[1]});
        my $selectinfo = join("\t",$tem[0], $number, $tag2tax{$tem[1]}); 
        #push @{$ar[3]}, $selectinfo;
        print OUT "$selectinfo\n";
    }
    close II; 
}##generate_otutable

sub generate_taxontable{
    my @taxf = @_;
    ##0. including the name.list file
    ##1. include the taxonomy file
    die "Can not the name.list file\n" unless open(I,"$taxf[0]");
    my %tag2sample;
    my %sample;
    while(<I>)
    {
        chomp;
        my @tem = split /\s+/;
        $sample{$tem[0]}++;
        $tag2sample{$tem[1]} = $tem[0];
    }## I
    close I;
    
    die "File not the taxonomy file\n" unless open(I,"$taxf[1]");
    <I>;
    my %taxonnum;
    my %taxonrank;
    while(<I>)
    {
        chomp;
        my @tem = split /\s+/;
        if(exists $tag2sample{$tem[0]})
        {
            $taxonnum{$tem[1]}{$tag2sample{$tem[0]}} ++;
            $taxonrank{$tem[1]} = $tem[3]; 
        }else
        {
            die "Wrong tag ID \n";    
        }
    }
    close I;
    die "Wrong out put dir or file\n" unless open(T,">$taxf[2]");
    my @sampl;
    for my $k (sort keys %sample)
    {
        push @sampl, $k;
    }

    my $head = join("\t", "Taxon",@sampl,"Rank");
    print T "$head\n";
    for my $key (sort keys %taxonnum)
    {
        print T "$key";
        for my $k (@sampl )
        {
            if(exists $taxonnum{$key}{$k})
            {
                print T "\t$taxonnum{$key}{$k}";
            }else
            {
                print T "\t0";    
            }
        }
        if(exists $taxonrank{$key})
        {
            print T "\t$taxonrank{$key}";
        }else
        {
            die "Wrong ID\n";    
        }
        print T "\n";
    }
}


1;

