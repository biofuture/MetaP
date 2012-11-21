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
our @EXPORT = qw(generate_otutable generate_taxontable tax_otu gastscvmothurotu otu_normalized1);
use SixteenS::Tools::Matrix;

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
    my $log = "QIIME format otu table"; 
    #push @{$ar[3]}, $log;
    print OUT "$log\n";
    my $mark = "OTUID\t$samplelist\tConsesus Lineage";
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

    #my $head = join("\t", "Taxon",@sampl,"Rank");
    my $head = join("\t", "Taxon",@sampl);
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
=head
        if(exists $taxonrank{$key})
        {
            print T "\t$taxonrank{$key}";
        }else
        {
            die "Wrong ID\n";    
        }
=cut
        print T "\n";
    }
}#generate_taxontable

sub tax_otu{

   # die "perl $0 [Name.list] [uniq.name] [otu.seq.name][otu.tax][out.tablei.prefix]\n" unless(@ARGV == 5);
   my @arg = @_;
   die "$! \n" unless open(IN1,"$arg[0]");
    die "$! \n" unless open(UN,"$arg[1]");
    die "$! \n" unless open(IN2,"$arg[2]");
    die "$! \n" unless open(IN3,"$arg[3]");
    my $out1 = "$arg[4].Act.xls";
    my $out2 = "$arg[4].Relative.xls";
    die "$! \n" unless open(OT,">$out1");
    die "$! \n" unless open(OUT,">$out2");

###----------------------------------------
#reads mapping to sample information stored in 
#hash %SnRn
    my %SnRn;
    my %Sn;
    while(<IN1>){
        chomp;
        my @m = split;
        $SnRn{$m[1]} = $m[0];
        $Sn{$m[0]} = 1;
    }
    close IN1;

    my @sname = sort keys %Sn;

    my %uniqab;
    while(<UN>){
        chomp;
        my @m = split /\s+/;
        my @tmp = split /\,/,$m[1];
        $uniqab{$m[0]} = scalar(@tmp);
    }
    close UN;

###---------------------------------------
#reads mapping to OTU ID information stored 
#in hash %otu and statistic OTU Sample abundance 
#my %otu;
    my %otus;
    my %abund;
    my %fsab;
    my %otuid2seq;
    while(<IN2>){
        chomp;
        my @m = split /\t/;
        my @tmp = split /\,/,$m[1];
        my $otun = "O_$tmp[0]";

        $otuid2seq{$otun} = $m[1]; ###store the all the otu seq information for the &addab() subroutine

            for my $id (@tmp){
#$otu{$id} = $otun;

                for my $sn (@sname){

                    if($sn eq  $SnRn{$id}){
                        $otus{$otun}{$sn}++;
                    }else{
                        $otus{$otun}{$sn} += 0;
                    }

                }
            }
        $abund{$otun} = $#tmp +1;

##caclulate the first and second ab
        my %catmp;
        for(@tmp){
            if(exists $uniqab{$_}){
                $catmp{$_} = $uniqab{$_};
            }
        }
        my $f = 0;
        my $rate=0;
        for my $ab (sort {$catmp{$b} <=> $catmp{$a}} keys %catmp){
            $f++;

            $rate = $catmp{$ab} / scalar(@tmp);
            if($f == 1){
                $fsab{$otun}{1} = "$catmp{$ab}\t$ab";
            }elsif($f == 2){
                $fsab{$otun}{2} = "$catmp{$ab}\t$ab";
            }
        }


        if(! exists $fsab{$otun}{2}){
            $fsab{$otun}{2} = "0\t0";
        }

        if(! exists $fsab{$otun}{1}){
            die "w 111\n";
        }

        for(keys %catmp){
            delete($catmp{$_});
        }
##---
    }
    close IN2;

###----------------------------------------
#statics tax infor
#
    my %taxheir;
    my %taxon;
    my %sabu;
    my %taxfsab;
    my $length = 0;
    my $treads;
    while(<IN3>){
        chomp;
        my @m = split /\t/;
        my $otun = "O_$m[0]";
        $treads += $abund{$otun};

        $taxon{$m[1]}{$otun} = 1;

        my @tmp = split /\;/,$m[1];
        $length =scalar(@tmp) if(scalar(@tmp)>$length);
        while(scalar(@tmp) > 0) {
            my $tax = join ";",@tmp;

            if(exists $abund{$otun}){

                if(exists $taxheir{$tax}){
                    $taxheir{$tax} += $abund{$otun};  
                }else{
                    $taxheir{$tax} = $abund{$otun}
                }

    #            &addab($tax,$otun);###add tax abund in every sample stored in %sabu
    die "wrong $otun in addab\n" unless (exists $otuid2seq{$otun});
    my @tmp = split /\,/,$otuid2seq{$otun};

    for my $seqid (@tmp){
        if(exists $SnRn{$seqid}){
            for my $sn (@sname){
                if($sn eq $SnRn{$seqid}){
                    $sabu{$tax}{$sn} ++;
                }else{
                    $sabu{$tax}{$sn} += 0;
                }
            }

        }else{
            die "wrong $seqid in addab\n";
        }
    }

###caclulate the first and second ab
    my %catmp;
    for(@tmp){
        if(exists $uniqab{$_}){
            $catmp{$_} = $uniqab{$_};
        }
    }
    my $f = 0;


    for my $ab (sort {$catmp{$b} <=> $catmp{$a}} keys %catmp){
        $f++;

        if($f == 1){
            if(exists $taxfsab{$tax}{1}){
                my @tm = split /\t/,$taxfsab{$tax}{1};
                if($tm[0] > $catmp{$ab}){
                    $taxfsab{$tax}{1} = "$catmp{$ab}\t$ab";
                }else{

                }

            }else{
                $taxfsab{$tax}{1} = "$catmp{$ab}\t$ab";
            }

        }elsif($f == 2){
            if(exists $taxfsab{$tax}{2}){
                my @tm = split /\t/,$taxfsab{$tax}{2};
                if($tm[0] > $catmp{$ab}){
                    $taxfsab{$tax}{2} = "$catmp{$ab}\t$ab";
                }else{
                }
            }else{
                $taxfsab{$tax}{2} = "$catmp{$ab}\t$ab";
            }
        }
    }

    if(! exists $taxfsab{$tax}{2}){
        $taxfsab{$tax}{2} = "0\t0";
    }

    if(! exists $taxfsab{$tax}{1}){
        die "w 111\n";
    }
    for(keys %catmp){
        delete($catmp{$_});
    }
##addab

            }else{
                die "OTU wrong\n";
            }
            pop @tmp;
        }
    }
    close IN3;

###----------------------------------------
##generate the relative abudance of every sample
    my %Rotus;
    my %otutotal;
    for my $otun(keys %otus){
        for my $sn(sort keys %{$otus{$otun}}){
            $otutotal{$sn} += $otus{$otun}{$sn};
        }
    }

    for my $otun(keys %otus){
        for my $sn(sort keys %{$otus{$otun}}){
            $Rotus{$otun}{$sn} = $otus{$otun}{$sn} / $otutotal{$sn};
            $Rotus{$otun}{all} += $otus{$otun}{$sn} / $otutotal{$sn};
        }
    }


###----------------------------------------
#generate the table
#
    my $st = join("\t",@sname);
    print  OT "0","\t"x($length+1),"all","\tfn\tFseqname\tsn\tSseqname\ttotal\t","$st\n";
    print  OUT "0","\t"x($length+1),"all","\tfn\tFseqname\tsn\tSseqname\ttotal\t","$st\n";
    for my $takey (sort keys %taxheir){
        my @tmp = split /\;/,$takey;
        my $rank = $#tmp + 1;
        my $num = $length - scalar(@tmp) + 1;

        die "$takey \n" unless($taxfsab{$takey}{1} && $taxfsab{$takey}{2} && $taxheir{$takey});
        print OT $rank,"\t"x(scalar(@tmp)),$tmp[-1],"\t"x($num+1),"$taxfsab{$takey}{1}\t$taxfsab{$takey}{2}\t","$taxheir{$takey}";

        print OUT $rank,"\t"x(scalar(@tmp)),$tmp[-1],"\t"x($num),"$taxheir{$takey}","\n";
        for my $si (@sname){
            print OT "\t$sabu{$takey}{$si}";
#print OUT "";
        }
        print OT "\n";
#print OUT "\n";
        if(exists $taxon{$takey}){
            for my $otid (sort keys %{$taxon{$takey}}){
                if(exists $otus{$otid}){
                    print OT "OTU\t$rank","\t"x($length-1),"\t","$otid\_$takey";
                    print OUT "OTU\t$rank","\t"x($length-1),"\t","$otid\_$takey";

                    if(exists $fsab{$otid}){
                        print OT "\t$fsab{$otid}{1}\t$fsab{$otid}{2}\t$abund{$otid}";
                        print OUT "\t$fsab{$otid}{1}\t$fsab{$otid}{2}\t$Rotus{$otid}{all}";
                    }else{
                        die "wrong $otid\n";
                    }

                    for my $si (sort @sname){
                        print OT "\t$otus{$otid}{$si}";
                        print OUT "\t$Rotus{$otid}{$si}";
                    }
                    print OT "\n";
                    print OUT "\n";
                }else{
                    die "wrong $otid\n";
                }
            }
        }

    }
}#tax_otu

sub gastscvmothurotu{
    ##This function generate the old gast.csv from total gast tax file and dnaclut mothur format otu 
    #1. input total.gast.tax
    #2. MetaP.otu.table
    #3. gast.csv
    #4. mothur.otu

    my @ar = @_;
    die "$ar[1] $!\n" unless open(I2, "$ar[1]");
    die "$ar[3] $! \n" unless open(T2,">$ar[3]");
    <I2>;<I2>;
    my %rep;
    while(<I2>)
    {
        chomp;
        my @tem = split /\s+/;
        print T2 "$tem[1]\t$tem[5]\n";
        $rep{$tem[1]} = 1;
    }
    close I2;
    close T2;
    
    die "$ar[0] $!\n" unless open(I1, "$ar[0]");
    die "$ar[2] $!\n" unless open(T1, ">$ar[2]");
    <I1>;
    while(<I1>)
    {
        chomp;
        my @tm = split(/\s+/,$_);
        if(exists $rep{$tm[0]})
        {
            print T1 "$tm[0]\t$tm[1]\n";    
        }
    }
    close I1;
    close T1;
}#gastcsvmothurotu

sub otu_normalized1{
#This function generated the normalized otu table for PCA nalysis 
#and clustering analysis
#1. the otu table of QIIME format 
#2. output otu table prefix 

    my @ar = @_;
    die "can not open $ar[0] $! \n" unless open(I, "$ar[0]");
    my $outA = "$ar[1].act.xls";

    #my $outR = "$prefix.relative.xls";
    die "can not write to $outA $! \n" unless open(TM,">$outA");
    #die "Can not write to $outR $! \n" unless open(TR, ">$outR");
    <I>;
    my $head = <I>;
    chomp($head);
    my @tem = split(/\t/,$head); pop @tem; shift @tem;
    $head = join("\t", "",@tem);
    #print "$head\n";
    print TM "$head\n";
    while(<I>)
    {
        my @tm = split(/\t/,$_);
        pop @tm;
        my $o = join("\t",@tm);
        print TM "$o\n";
    }
    close I;
    close TM;
    die "file $outA not exists $! \n" unless(-e "$outA");
    normalized($outA);

}#otu_normalized1

1;

__END__
