package  Sixteen::Tools::Matrix;
use strict;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(normalized);

#This package is mainly to process matrix
#such as normalize matrix by column
sub normalized{
    #F P  is matrix file

    my @tm = @_;
    
    die "This File does not exits $tm[0] $!\n" unless open(I, "$tm[0]") ;
    my $head = <I>; chomp($head);
    my @nam = split(/\s+/,$head);
    my $num = $#nam;
    my %sum;
    while(<I>)
    {
        my @tem = split(/\s+/,$_);
        for(my $i =1; $i <= $num; $i++)
        {
            $sum{$i}+=$tem[$i];
        }
    }
    close I;

    die "Normalized $!\n" unless open(I,"$tm[0]");
    my $out = "Normaliezed_1.$tm[0]";
    die "OutPut file $out\n" unless open(T,">$out");
    print T "$head\n";
    <I>;
    while(<I>)
    {
        my @tem = split(/\s+/,$_);
        print "$tem[0]";
        for(my $i =1; $i <= $num; $i++)
        {
            my $on = $tem[$i] / $sum{$i};
            print T "\t$on";
        }
        print T "\n";
    }
    close I;
    close T;
}#sub function
