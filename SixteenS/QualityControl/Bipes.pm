###
package SixteenS::QualityControl::Bipes;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(overlap_fq primer_cut);

use strict;

BEGIN{
    unshift @INC, "../../";
}

use File::Basename;
use SixteenS::Tools::SequenceProcess;
use threads;

sub overlap_onef{

    ## There are five parameters in this function
    #1 1.fq 
    #2 2.fq
    #3 output prefix
    #4 mix score 
    #5 max score
    my @ar = @_;

    open FILE1, $ar[0] || die "Can not open FILE1: $!";
    open FILE2, $ar[1] || die "Can not open FILE2: $!";

    my $minscore = $ar[3];
    my $maxscore = $ar[4];
    #print "$ar[0]\n$ar[1]\n$ar[2]\n$ar[3]\n$ar[4]\n";
    #die;
    my $rangenum = 0;
    my $totalnum = 0;
    my @align_length;
    my $rangematchfile = $ar[2]."_".$ar[3]."_".$ar[4]."_match.fa";

    open SCOREFILE,">$ar[2].distribution.score" || die "Can not open SCOREFILE: $!";
    open MATCH, ">$ar[2].match.fa" || die "Can not open MATCH: $!";
    open MATCHFQ, ">$ar[2].match.fq" || die "Can not open MATCHFQ: $!";
    open MERGERFILE, ">$ar[2].alignment" || die "Can not open MERGERFILE: $!";
    open RANGEMATCH, ">$ar[2].$ar[3]_$ar[4].match.fa" || die "Can not open RANGEMATCH: $!";
    open REPORT, ">$ar[2].alignment_information.list" || die "Can not open REPORT: $!";
    open MISMATCHFILE,">$ar[2].distribution.mismatch" || die "Can not open MISMATCHFILE: $!";
    open MISMATCH_POSITION,">$ar[2].all_position.mismatch" || die "Can not open MISMATCH_POSITION: $!";
    open B40_70MISMATCH_POSITION,">$ar[2].40_70_position.mismatch" || die "Can not open B40_70MISMATCH_POSITION: $!";

    printf REPORT "%-40s %-9s %-17s %-17s %-17s %-6s\n", "SeqID", "Length", "Identity", "Similarity", "Gaps", "Score";
    print MISMATCH_POSITION "mismatch site: \n";
    print B40_70MISMATCH_POSITION "mismatch site between 40 and 70: \n";

    my @score_array = ();
    my @mismatch_array2 = ();

    while(<FILE1>)
    {
        if ($_ eq "\n" || $_ eq "\n\r")
        {
            next;
        }
        my $name1 = substr $_, 1, (length $_) - 4;
        my $name2 = (substr <FILE2>, 1, (length $_) - 4);
        my $seq1;
        my $seq2;
        my $quality1;
        my $quality2;
        my @myqual3;

        chomp($seq1 = uc <FILE1>);
        chomp($seq2 = uc <FILE2>);
        <FILE1>;
        <FILE2>;
        chomp($quality1 =<FILE1>);
        chomp($quality2 =<FILE2>);
        my @myqual1 = split //,$quality1;
        my @myqual2 = split //,$quality2;
        @myqual2 = reverse @myqual2;
        my $l = 0;
        while($myqual2[$l] eq " ")
        {
            shift @myqual2;
        }
        while($myqual1[$l] eq " ")
        {
            shift @myqual1;
        }
        foreach my $k (0..($#myqual1))
        {
            $myqual3[$k] = $myqual1[$k];
        }
        if ($name1 ne $name2)
        {
            warn "$name2 is not the corresponding end of $name1, even though they are in the same position of fastq file.";
            next;
        }   

        open TEMP1, ">$ar[2]temp1.seq" || die "Can not open TEMP1: $!";
        open TEMP2, ">$ar[2]temp2.seq" || die "Can not open TEMP2: $!";
        print TEMP1 ">", $name1, "\n", $seq1, "\n";
        print TEMP2 ">", $name1, "\n", $seq2, "\n";
        close TEMP1;
        close TEMP2;

        system "merger $ar[2]temp1.seq $ar[2]temp2.seq -sreverse2 -outseq $ar[2]temp3.seq -outfile $ar[2]temp4 -datafile EDNAMAT -auto";

        open TEMP4, "$ar[2]temp4" || die "Can not open TEMP4: $!";
        my $lengthseq;
        my $identity;
        my $similarity;
        my $gaps;
        my $score;
        my $str1;
        my $str2;
        my $str3;
        my @array;

        open TEMP3, "$ar[2]temp3.seq" || die "Can not open TEMP3: $!";
        <TEMP3>;
        push @array,">".$name1."\n";
        while(<TEMP3>)
        {
            push(@array,$_);
        }
        close TEMP3;
        $totalnum += 1;
        my $temp1 = <TEMP4>;
        until($temp1 =~ /#={6,}/)
        {
            $temp1 = <TEMP4>;
        }
        print MERGERFILE $temp1,"#\n","# SeqID : $name1\n";
        foreach my $k (0..6)
        {
            <TEMP4>;
        }
        my @position_array1 = ();
        my @position_array2 = ();
        my $gaplen1 = 0;
        my $mybool = 0;
        my $testbool = 0;
        while(<TEMP4>)
        {
            print MERGERFILE $_;
            if (/#\s+Length:\s+(\S+)\s*/)
            {
                $lengthseq=$1;
                $gaplen1 = $1 - (length $seq1);
                if ($lengthseq < 97)
                {
                    $testbool = 1;
                }
                foreach my $k (0..($gaplen1-1))
                {
                    $myqual3[$k] = $myqual1[$k];
                }
                foreach my $k (0..$#myqual2)
                {
                    if (($k+$gaplen1) <= $#myqual1)
                    {
                        if ($myqual1[$k+$gaplen1] gt $myqual2[$k])
                        {
                            $myqual3[$k+$gaplen1] = $myqual1[$k+$gaplen1];
                        }
                        else
                        {
                            $myqual3[$k+$gaplen1] = $myqual2[$k];
                        }
                    }
                    else
                    {
                        $myqual3[$k+$gaplen1] = $myqual2[$k];
                    }
                }
            }
            elsif (/#\s+Identity:\s+(\S+)\s+\(\s*([0-9%.]+)\)\s*/)
            {
                $identity=$1.'('.$2.')';
            }
            elsif (/#\s+Similarity:\s+(\S+)\s+\(\s*([0-9%.]+)\)\s*/)
            {
                $similarity=$1.'('.$2.')';
            }
            elsif (/#\s+Gaps:\s+(\S+)\s+\(\s*([0-9%.]+)\)\s*/)
            {
                $gaps=$1.'('.$2.')';
            }
            elsif (/#\s+Score:\s+(\S+)\s*/)
            {
                $score=$1;
                $score_array[$1] +=1;
            }
            elsif (/([:#\w\/]+\s+\d+\s+)([ACGTNnacgt-]+)\s+(\d+)/)
            {
                if($mybool == 0)
                {
                    my $mytempseq01 = substr $2,0,1;
                    my $mytempseq02;
                    $mybool += 1;
                    if ($mytempseq01 eq "-")
                    {
                        $testbool = 1;
                    }
                }
            }
            if (/#\s+(\d+)\s+'(\w)'\s+(\d+)\s+'(\w)'\s+'(\w)'/)
            {           
                if ($score >= $minscore && $score <= $maxscore)
                {
                    if($1 >= 40 && $1 <= 70)
                    {
                        push @position_array2, $1;
                    }
                    $mismatch_array2[$1] += 1; 
                    push @position_array1, $1;
                }
                if ( $1 < ((length $seq2)-$3))
                {                
                    $str2 = uc $2;
                    $myqual3[$1] = $myqual1[$1];
                }
                elsif ( $1 > ((length $seq2)-$3))
                {
                    $str2 = uc $4;
                    $myqual3[$1] = $myqual2[$1-$gaplen1];
                }
                else
                {
                    my $q1 = substr $quality1,($1-1),1;
                    my $q2 = substr $quality2,(-$3),1;
                    if ( $q1 gt $q2)
                    {
                        $str2 = uc $2;
                        $myqual3[$1] = $myqual1[$1];
                    }
                    elsif ( $q1 lt $q2)
                    {
                        $str2 = uc $4;
                        $myqual3[$1] = $myqual2[$1-$gaplen1];
                    }
                    else
                    {
                        $myqual3[$1] = $myqual1[$1];
                        my $j = $3>5?5:($3-1);
                        while( $j < $3 )
                        {
                            my @qul1 = unpack("C*",(substr $quality1,($1-$j-1),(2*$j+1)));
                            my @qul2 = unpack("C*",(substr $quality2,($1-$j-1),(2*$j+1)));
                            my $sum1 = 0;
                            my $sum2 = 0;
                            foreach $q1 (@qul1)
                            {
                                $sum1 += $q1;
                            }
                            foreach $q2 (@qul2)
                            {
                                $sum2 += $q2;
                            }
                            if ( $sum1 > $sum2)
                            {
                                $str2 = uc $2;
                                last;
                            }
                            elsif ( $sum1 < $sum2)
                            {
                                $str2 = uc $4;
                                last;
                            }
                            else
                            {
                                $str2 = uc $5;
                                if ($j >= ($3-1))
                                {
                                    last;
                                }
                                $j = $3 >($j+5)?($j+5):($3-1);
                            }                        
                        }
                    }
                }                                  
                if (($1+1) <= (length $array[1]))
                {
                    $str1 = substr $array[1], 0, ($1-1);
                    $str3 = substr $array[1], $1;
                    $array[1]= $str1.$str2.$str3;
                }
                else
                {
                    $str1 = substr $array[2], 0, ($1-(length $array[1]));
                    $str3 = substr $array[2], ($1-(length $array[1])+1);
                    $array[2]= $str1.$str2.$str3;
                }
            }
        }
        close TEMP4;
        if ($testbool == 0)
        {

            printf MISMATCH_POSITION "%-37s",$name1.":";
            printf B40_70MISMATCH_POSITION "%-37s",$name1.":";
            printf MISMATCH_POSITION "(%d)  ",$#position_array1+1;
            print MISMATCH_POSITION "@position_array1\n";
            printf B40_70MISMATCH_POSITION "(%d)  ",$#position_array2+1;
            print B40_70MISMATCH_POSITION "@position_array2\n";
            foreach my $j (1..$#array)
            {
                chomp($array[$j]);
            }
            print MATCH @array,"\n";
            my @array_q = @array;
            shift @array_q;
            my $myname2 = "@".$name1."\n";
            unshift @array_q,$myname2;
            print MATCHFQ @array_q,"\n";
            print MATCHFQ "+","\n",@myqual3,"\n";
            my $abc002 = ($#myqual3 + 1);
            if ($lengthseq != $abc002)
            {
                print $name1,": $lengthseq ---- $abc002\n";
            }
            if ($score >= $minscore && $score <= $maxscore)
            {
                print RANGEMATCH @array,"\n";
                $rangenum += 1;
                $align_length[$lengthseq] += 1;
            }
            print MERGERFILE "\n\n\n\n";
            printf REPORT "%-40s %-9s %-17s %-17s %-17s %-6s\n", $name1, $lengthseq, $identity, $similarity, $gaps, $score;
        }
    }

    close FILE1;
    close FILE2;
    close MATCH;
    close MATCHFQ;
    close RANGEMATCH;
    close REPORT;
    close MERGERFILE;
    open DISTRIBUTION ,">$ar[2].distribution.length" || die "Can not open DISTRIBUTION: $!";
    printf DISTRIBUTION "The number of alignment(%g - %g): %d (%.2f%%)\n\n",$minscore,$maxscore,$rangenum,100*$rangenum/$totalnum;
    print DISTRIBUTION "The distribution of length of aligment: \n";
    my $i = 0;
    print DISTRIBUTION "Length  total\n";
    while(!defined $align_length[$i])
    {
        $i += 1;
    }
    foreach my $j ($i..($#align_length))
    {
        if (!defined $align_length[$j])
        {
            $align_length[$j] = 0;
        }
        print DISTRIBUTION "$j: $align_length[$j] \n";
    }
    $i = 0;
    printf SCOREFILE "The number of alignment(%g - %g): %d (%.2f%%)\n\n",$minscore,$maxscore,$rangenum,100*$rangenum/$totalnum;
    print SCOREFILE "Score  total\n";
    while(!defined $score_array[$i])
    {
        $i += 1;
    }
    foreach my $j ($i..($#score_array))
    {
        if (!defined $score_array[$j])
        {
            $score_array[$j] = 0;
        }
        print SCOREFILE "$j: $score_array[$j] \n";
    }

    printf MISMATCHFILE "The number of alignment: %d \n",$rangenum;
    print MISMATCHFILE "\nMismatch  total\n";
    $i = 0;
    while(!defined $mismatch_array2[$i])
    {
        $i += 1;
    }
    foreach my $j ($i..($#mismatch_array2))
    {
        if (!defined $mismatch_array2[$j])
        {
            $mismatch_array2[$j] = 0;
        }
        my $mismatch_rated = $mismatch_array2[$j]/$rangenum*100;
        printf MISMATCHFILE "%d: %d (%.2f%%)\n",$j,$mismatch_array2[$j],$mismatch_rated;
    }
    unlink "$ar[2]temp3.seq","$ar[2]temp1.seq","$ar[2]temp2.seq","$ar[2]temp4";
    close DISTRIBUTION;
    close MISMATCHFILE;
    close SCOREFILE;

}##overlap_onef

sub overlap_fq{ 
    ##read fastq file1 file2 and overlap these reads
    ##use multi threads to overlap fastq reads, when ever the reads reached
    ##a number then create a thread
    my @ar = @_;
    #die "Not a fastq file\n" unless (-e $ar[0])  && open(FQ1,$ar[0]);
    #die "Not a fastq file\n" unless (-e $ar[1])  && open(FQ1,$ar[1]);
    ##mkdir and split fastq file
    unless(-d "overlap")
    {
        `mkdir overlap`;
    }else
    {
        `rm -rf overlap; mkdir overlap`;
    }
    split_file($ar[0], $ar[2]);
    split_file($ar[1], $ar[2]);
    my $ba1 = basename($ar[0]);
    my $ba2 = basename($ar[1]);
    
    ##make use of multi thread method to paralell overlap
    my @thre;
    my @record;
    for(my $i = 1; $i <= $ar[2]; $i++)
    {
        unless(-d "overlap/$i")
        {
             `mkdir overlap/$i; mv $ar[0].$i.split overlap/$i/; mv $ar[1].$i.split overlap/$i/`;    
        }
        my $oprefix = "overlap/$i/oprefix.$i";
        $thre[$i] =  threads->create(\&overlap_onef, "overlap/$i/$ba1.$i.split", "overlap/$i/$ba2.$i.split","$oprefix", "-100","1000");
    }
    ##join all the threads 
    for(my $i = 1; $i <= $ar[2]; $i++)
    {
        $record[$i] = $thre[$i]->join();   
    }
    ##merge all the split results into the target one     
    `cat overlap/*/*.-100_1000.match.fa > overlaped.total.pe.fa`;
}##overlap_fq

sub primercut{

    ##This function is to cut primer
    ##There are three parameters 
    #1. primer with one row one primer seq
    #2. the overlaped pe file
    #3. The output prefix 
    my @ar= @_;

    open FILE1, $ar[0] || die "Can not open FILE1: $!";
    open FILE2, $ar[1] || die "Can not open FILE2: $!";

    my $outfile1 = $ar[2].".primercut.fa";
    my $outfile3 = $ar[2].".distribution_primer_error.number";
    my $outfile5 = $ar[2].".barcodeprimer";

    open OUTFILE1, ">$outfile1" || die "Can not open OUTFILE1: $!";
    open OUTFILE3, ">$outfile3" || die "Can not open OUTFILE3: $!";
    open BP,">$outfile5" || die "can't create the barcode primer file\n";
#open R,">$ar[3]" or die "can'te create the bool file\n";

    my %hash=(
            "R"=>"A,G",
            "Y"=>"C,T",
            "M"=>"A,C",
            "K"=>"G,T",
            "S"=>"G,C",
            "W"=>"A,T",
            "H"=>"A,T,C",
            "B"=>"G,T,C",
            "V"=>"G,A,C",
            "D"=>"G,A,T",
            "N"=>"A,T,G,C",
            );

    my $frontseq = <FILE1>;
    chomp($frontseq);

    die "the forward primer contains merger bases that the program doesn't recognize\n" if ($frontseq=~/![ATGCatgcRYMKSWHBVDN]/);

    my @frontseq;
    $frontseq[0]=$frontseq;
    while ($frontseq[0]=~/[RYMKSWHBVDN]/) #if the unfolding is not completed, continue this circle
#foreach (1..3)
    {
        my @temp=();
        foreach my $temp(@frontseq)
        {
            my $base_p=0;
            my @temp2=();
            my @temp_base=split(//,$temp);
            my @base_un;
            my $bool2=0;
            my %base=();#$base{0} are the bases before the merger base and $base{1} are the bases afer the merger base
                foreach my $temp_base(@temp_base)
                {
                    if ($temp_base=~/([RYMKSWHBVDN])/ && $bool2==0)
                    {
                        $base_p++;
                        @base_un=split(/,/,$hash{$1});
                        $bool2=1;
                    }else{
                        $base{$base_p}.=$temp_base;
                    }
                }
            foreach (@base_un)
            {
                push @temp,$base{0}.$_.$base{1};
            }

        }
        @frontseq=@temp;
    }
    my $fp_length=length $frontseq;

    my $backseq = <FILE1>;
    chomp($backseq);

    die "the backward primer contains merger bases that the program doesn't recognize\n" if ($backseq=~/![ATGCatgcRYMKSWHBVDN]/);
    my @backseq;
    $backseq[0]=$backseq;
    while ($backseq[0]=~/[RYMKSWHBVDN]/)
    {
        my @temp=();
        foreach my $temp(@backseq)
        {
            my $base_p=0;
            my @temp2=();
            my @temp_base=split(//,$temp);
            my @base_un;
            my $bool2=0;
            my %base=();
            foreach my $temp_base(@temp_base)
            {
                if ($temp_base=~/([RYMKSWHBVDN])/ && $bool2==0)
                {
                    $base_p++;
                    @base_un=split(/,/,$hash{$1});
                    $bool2=1;
                }else{
                    $base{$base_p}.=$temp_base;
                }
            }
            foreach (@base_un)
            {
                my $primer_temp=$base{0}.$_.$base{1};
                $primer_temp=reverse $primer_temp;
                $primer_temp=~tr/ATGCatgc/TACGTACG/;
                push @temp,$primer_temp;
            }

        }
        @backseq=@temp;
    }
    my $bp_length=length $backseq;

    close FILE1;

    print OUTFILE3 "Primer error number:\n";
    while(<FILE2>)
    {
        my %hash=();
        my $s1;
        my $s2;
        my $p1;
        my $p2;

        my $iden1=0;
        my $iden2=0;
        my $error_num = 0;
        my $myname0;
        my $bool=0;
        chomp($s1 = $_);
        chomp($s2 = <FILE2>);
        my $length=length $s2;
        $s2=~tr/atgc/ATGC/;
C:foreach my $c(0..1)
  {

      if ($c==1)
      {
          $s2=reverse $s2;
          $s2=~tr/ATGCatgc/TACGTACG/;
          #print $s2,"\n";
      }

      foreach my $fp(@frontseq)
      {
          if ($s2=~/.+($fp)/)
          {
              $iden1=1;
              $p1=$+[1];
              last;
          }
      }
      foreach my $bp(@backseq)
      {
          if ($s2=~/($bp)/)
          {
              $iden2=1;
              $p2=$-[1]+1;
              last;
          }
      }

      $myname0 = substr $s1,1;
      if ($iden1 == 1 && $iden2 == 1)
      {             
          #print "yes\n";
          my $seq0 = substr $s2,$p1,($p2-$p1-1);
          print OUTFILE1 $s1,"\n",$seq0,"\n";
          my $barpri1 = substr $s2,0,$p1;
          my $barpri2 = substr $s2,-($length-$p2+1);
          $barpri2=reverse $barpri2;
          $barpri2=~tr/ATGCatgc/TACGTACG/;
          if ((length $barpri1 > 32) or (length $barpri2 > 32))#if the length of barpri is too long, it may be located in the wrong place
          {
              next;
          }
          print BP $s1," ",$barpri1," ",$barpri2,"\n";
          print OUTFILE3 $myname0,":   ";
          printf OUTFILE3 "(%d)\n",$error_num;
          $bool=1;
          last C;
      }
  }
  next if $bool == 1;

  $error_num=20;

  print OUTFILE3 $myname0,":   ";
  printf OUTFILE3 "(%d)\n",$error_num;
    }


    close FILE2;
    close OUTFILE1;
    close OUTFILE3;

}#primercut

sub primer_0_error_and_lt2_mismatch {
    ##1. error.number
    ##2. position.mismatch
    ##3. outprefix
    my @ar = @_;

open FILE1,"$ar[0]" || die "can not open FILE1:$! ";
open FILE2,"$ar[1]" || die "can not open FILE2I:$! ";
<FILE1>;
<FILE2>;
my $out=$ar[2].".primer_0_error_and_40_70_mismatch.names";
open FILE3,">$out" || die "can not open FILE3:$!";
print FILE3 "seqID (this file contains the name of sequences which have zero error in primer ,zero or one mismatch in 40-70 site)","\n";

my %primere;
while(<FILE1>)
{
chomp;
if(/\s*([A-Za-z0-9_#:-]*):\s*\(\s*(\d+)\s*\)/)
{
my $seqn = $1;
my $num = $2;
if($num == 0)
{
$num++;
$primere{$seqn} = $num;  
}
## print "$seqn\t$num\n";
}
}

my @n;
while(<FILE2>)
{
chomp;
if(/\s*([A-Za-z0-9_#:-]*):\s*\(\s*(\d+)\s*\)/)
{
if((exists $primere{$1})&& $2<=1)
{
print FILE3 "$1\n";
}
}
}
##print "$#n\n";
close FILE1;
close FILE2;
close FILE3;

}##primer_0_error_and_lt2_mismatch


sub find_fa_new {
    #*.primer_0_error_and_40_70_mismatch.names
    #*.primercut.total.fa
    #*oprefix
    my @ar = @_;
    
    open FILE0, $ar[0] || die "Can not open FILE0: $!";
    open FILE1, $ar[1] || die "Can not open FILE1: $!";
    <FILE0>;
    unless(-d "$ar[2]")
    {
        mkdir "$ar[2]",0755, or warn "Cannot make directory: $!";
    }

    my @dir = split/\//,$ar[2];
    my $outfile = "./$ar[2]/".$dir[-1].".unique.tag.clean.fa";


    open OUTFILE1, ">$outfile" || die "Can not open OUTFILE1: $!";

    my %seq;
    while(<FILE0>)
    {
        chomp;
        $seq{$_}=0;   
    }

    while(<FILE1>)
    {
        if(/>/)
        {
            chomp;
            my $id = $_;
            my $na = substr $id,1;
            if (exists $seq{$na})
            {
                my $temp1 = <FILE1>;
                my $temp = uc ($temp1);     
                print OUTFILE1 "$id\n$temp"; 
            }
        }
    }
    close FILE0;
    close FILE1;
    close OUTFILE1;
}#find_fa_new

sub primer_cut{

    #There are two parameters for this function
    #1. the primer file
    #2. the total overlap file
    my @ar = @_;
    #unique total overlap file
    `mothur "#unique.seqs(fasta=$ar[1])"`;
    my $uq = $ar[1];
    $uq =~ s/\.fa$//;
    $uq = "$uq.unique.fa";

    split_file($uq,$ar[2]);

    unless (-d "primercut")
    {
        `mkdir primercut`;
    }else{
        `rm -rf primercut; mkdir primercut`;
    }
    
    ##make use of multi thread method to paralell overlap
    my @thre;
    my @record;
    for(my $i = 1; $i <= $ar[2]; $i++)
    {
        unless(-d "primercut/$i")
        {
             `mkdir primercut/$i; mv $uq.$i.split primercut/$i/;`;    
        }
        my $oprefix = "primercut/$i/oprefix.$i";
        $thre[$i] =  threads->create(\&primercut, "$ar[0]", "primercut/$i/$uq.$i.split","$oprefix");
    }
    ##join all the threads 
    for(my $i = 1; $i <= $ar[2]; $i++)
    {
        $record[$i] = $thre[$i]->join();   
    }

    ##Statistic all the qualified files
    unless(-d "qualified_sequences")
    {
        `mkdir qualified_sequences`;
    }

    ##Files 
    my $op = "qualified_sequences/qualified";
    `cat overlap/*/*40_70_position.mismatch >$op.40_70_position.mismatch`;
    `cat primercut/*/*distribution_primer_error.number >$op.distribution_primer_error.number`;
    `cat primercut/*/*primercut.fa >$op.primercut.total.fa`;
    `cat primercut/*/*.barcodeprimer > $op.barcodeprimer`;
   
    primer_0_error_and_lt2_mismatch("$op.distribution_primer_error.number","$op.40_70_position.mismatch",$op);
    
    find_fa_new("$op.primer_0_error_and_40_70_mismatch.names", "$op.primercut.total.fa","$op");
    
#    unless(-d "sample_seq")
#    {
#    `mkdir sample_seq`; 
#    }

    my $uqna = $ar[1];
    $uqna =~ s/\.fa$//;
    $uqna = "$uqna.names";
    #print "$ar[3]\n$uqna\n$op/qualified.unique.tag.clean.fa\n$op.barcodeprimer\n";
    barcode_select_samples($ar[3], $uqna, "$op/qualified.unique.tag.clean.fa","$op.barcodeprimer");

}#primer_cut

sub barcode_select_samples{

#1. xls table primer sample id information
#2. mothur unique .name 
#3. qualified unique .tag
#4. barcodefile from primercut

my @ar = @_;
open FILE1,"$ar[0]" || die "can not open FILE1:$!";
open FILE2,"$ar[1]" || die "can not open FILE2:$!";
open FILE3,"$ar[2]" || die "can not open FILE3:$!";

my %seqname;
while(<FILE2>)
{
chomp;
my @otu = split/\s+/;
$seqname{$otu[0]} = $otu[1];
}
close FILE2;

my %fasta;
while(<FILE3>)
{
chomp;
if(/^>/)
{
my $name = substr $_,1;
my $seq = <FILE3>;
$fasta{$name} = $seq;
}
}
close FILE3;

my @dir = split/\//,$ar[0];
pop @dir;
my $prefix;
if ($#dir >=0)
{
$prefix = join "/",@dir;
##print "$prefix\n";
}
else
{
$prefix = 1;
}
while(<FILE1>)
{
## print "$_\n";die;
chomp;
my @m = split/\./;
pop @m;
my $dirnew = join ".",@m;
##print "$dirnew\n";die;
system ("mkdir $dirnew");
if (!($prefix eq 1 ))
{
my $file = "$prefix\/$_";
open II,"$file" || die "can not open II:$!";
##print "$file\n";
}
else
{
    open II,"$_" || die "can not open II:$!";
##print "cwd:$_\n";
}
<II>;
while (<II>)
{
## print "$_\n";
    chomp;
    my @column = split/\t/;
    if ($#column == 5)
    {
        my $dir = "$dirnew\/$column[-1]";
        system ("mkdir $dir");
        my $out = "$dir\/$column[-1].tag.clean.fa";
##print "$out\n";die;
        open III,">$out"||die "can not open III:$!";
        my $flength = (length $column[2]) - 18;
        my $blength = (length $column[4]) - 14;  
        my $forbarcode = substr $column[2],0,$flength;
        my $b1 = substr $column[4],0,$blength;
        my $backbarcode = uc($b1);     
#   my $forbarcode = $column[2];
#   my $b1 = $column[4];
#   my $backbarcode = uc($b1);
#       print $forbarcode,"\n";
##  print "flength:$flength\tforebarcode:$forbarcode\tblength:$blength\tbackbarcode:$backbarcode\n";
#       print $column[2],"\n",$column[4],"\n";  
        open FILE4,"$ar[3]" || die "can not open FILE4:$!";
        while(<FILE4>)
        {
## print "$_";die;
            chomp;
            my @primer = split/\s+/;
# print "$primer[1]\n$primer[2]\n";
            my $fbarlocal = index $primer[1],$forbarcode;
            my $backbarlocal = index $primer[2],$backbarcode;
            my $seqn = substr $primer[0],1;
##   print "967F:$primer[1]\tfbarlocal:$fbarlocal\t1046R:$primer[2]\tbackbarlocal:$backbarlocal\n";die;  
#   print $primer[1],"\n",$column[2],"\n";
            if (($fbarlocal == 1 || $fbarlocal ==0 ) && ($backbarlocal ==1 || $backbarlocal ==0))
# && (length $primer[1] == length $column[2]) && (length $backbarcode == length $column[4]))
            {
#       print $primer[1],"\n",$primer[2],"\n";
                if ((length $primer[1] == length $column[2]) && (length $primer[2] == length $column[4]))
                {
                    if (exists $seqname{$seqn})
                    {
                        my @temp = split/\,/,$seqname{$seqn};
                        if (exists $fasta{$seqn})
                        {
                            for my $i(@temp)
                            {
                                print III ">$i\n$fasta{$seqn}";   
                            }
                        }
                    }
                }
            } 
 #               print "here\n";
        }
        close FILE4;
    } 
}
close FILE2;
}

close FILE1;
}#barcode_select_samples

1;
__END__
