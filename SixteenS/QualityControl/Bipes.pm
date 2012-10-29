###
package SixteenS::QualityContorl::Bipes;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(overlap_fq);

use strict;
BEGIN{
    unshift @INC, "../../";
}

use SixteenS::Tools::SequenceProcess;
use threads;

sub overlap_onef{

    my @ar = @_;

    open FILE1, $ar[0] || die "Can not open FILE1: $!";
    open FILE2, $ar[1] || die "Can not open FILE2: $!";

    my $minscore = $ar[3];
    my $maxscore = $ar[4];
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

        open TEMP1, ">temp1.seq" || die "Can not open TEMP1: $!";
        open TEMP2, ">temp2.seq" || die "Can not open TEMP2: $!";
        print TEMP1 ">", $name1, "\n", $seq1, "\n";
        print TEMP2 ">", $name1, "\n", $seq2, "\n";
        close TEMP1;
        close TEMP2;

        system "merger temp1.seq temp2.seq -sreverse2 -outseq temp3.seq -outfile temp4 -datafile EDNAMAT -auto";

        open TEMP4, "temp4" || die "Can not open TEMP4: $!";
        my $lengthseq;
        my $identity;
        my $similarity;
        my $gaps;
        my $score;
        my $str1;
        my $str2;
        my $str3;
        my @array;

        open TEMP3, "temp3.seq" || die "Can not open TEMP3: $!";
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
    unlink "temp3.seq","temp1.seq","temp2.seq","temp4";
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
    `mkdir split;cd split;link -n ../$ar[0]; link -n ../$ar[1]`;
    split_file($ar[0], $ar[2]);
    split_file($ar[1], $ar[2]);
    
    ##make use of multi thread method to paralell overlap
    my @thre;
    my @record;
    for(my $i = 1; $i <= $ar[3]; $i++)
    {
         $thre[$i] =  threads->create(\&overlap_onef, "$ar[0].$i.split", "$ar[1].$i.split","", "-100"."1000");
    }
    ##join all the threads 
    for(my $i = 1; $i <= $ar[3]; $i++)
    {
         $record[$i] = $thre[$i]->join();   
    }
    
    ##merge all the split results into the target one     
}##overlap_fq

sub remove_primerwrong{


}##remove_primerwrong

sub part_seqs2sample{



}##part_seqs2sample

1;
__END__
