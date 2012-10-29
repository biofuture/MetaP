package SixteenS::Tools::SequenceProcess;
use strict;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(fetch_fastainfo split_file);

sub fetch_fastainfo{
    
    my ($fa, $idseqref, $indexfaref) = @_;
    ##the $idseqref and $indexfaref are hash reference 
    ##we could modify these hash in this subroution
    
    die "Not a fasta file\n" unless open(I,"$fa");
    my $index=0;
    $/ = ">";
    <I>;
    while(<I>)
    {
        chomp;
        my @tem = split /\n/;
        my $fname = shift @tem;
        $fname = (split(/\s+/,$fname))[0];
        my $seq = join("",@tem);
        ${$idseqref}{$fname} = $seq;
        ${$indexfaref}{$index} = $fname;
        $index++;
    }
    close I;
    $/ = "\n";
}


sub split_file{
    ##split_file split file into n parts based on line
    ##the output file were in the same dir with the $file
    my ($file, $num) = @_;
    die "Cant not find $file $!\n" unless open(I,"$file");
    my $totallines = 0;
    my $tem = `wc -l $file`;
    chomp($tem);
    $totallines = (split(/\s+/,$tem))[1];
    my $line = int( $totallines / $num ); 
    #my $restline = $totallines - $num * $line;
    
    for (my $i = 1; $i < $num; $i++)
    {
        my $outf = join(".",$file,$i,"split");
        die "Can not write to $outf $!\n" unless open(T,">$outf");
        my $flag = 1;
        my $seq = "";
        while($flag <= $line)
        {
            $seq = <I>;
            print T "$seq";
            $flag++
        }
        close T;
    }
    
    my $loutf = "$file.$num.split";
    die "Can not write to $loutf $!\n" unless open(T,">$loutf");
    while(<I>)
    {
        print T "$_";
    }
    close I;
}##split_file
1;
__END__
