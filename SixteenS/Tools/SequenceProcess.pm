package SixteenS::Tools::SequenceProcess;
use strict;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(fetch_fastainfo split_file rename_tags);

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
    ##split_file split fastq file into n parts based on line
    ##the output file were in the same dir with the $file
    my ($file, $num) = @_;
    die "Cant not find $file $!\n" unless open(I,"$file");
    my $totallines = 0;
    my $tem = `wc -l $file`;
    chomp($tem);
    $totallines = (split(/\s+/,$tem))[1];

    my $line = int( $totallines / ($num * 4 ) ); 
    #my $restline = $totallines - $num * $line;
    
    for (my $i = 1; $i < $num; $i++)
    {
        my $outf = join(".",$file,$i,"split");
        die "Can not write to $outf $!\n" unless open(T,">$outf");
        my $flag = 1;
        my $seq = "";
        while($flag <= ($line*4))
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

sub rename_tags{
    #This funciton is mainly for rename all the clean 
    #nochimera tags for convenient of the following analysis
    #then generate the name.list file and the total tags file
    #The 1st input is the list of all the nochimera tags
    #the 2ed input paramerer is the output dir
    my @ar = @_;
    unless(-d "$ar[1]")
    {
        `mkdir $ar[1]`;
    }
    my $onamelist = "$ar[1]/total.name.list";
    die "$!\n" unless open(NAL,">$onamelist");
    my $ototalfa = "$ar[1]/total.nochimera.fa";
    die "$!\n" unless open(TA, ">$ototalfa");

    die "can not open $ar[0] $!\n" unless open(REN, "$ar[0]");
    while(<REN>)
    {
        chomp;
        die "$! \n"  unless open(I, "$_");
        my @tm = split(/\//, $_);
        unless(-d "$ar[1]/$tm[-2]"){
            `mkdir $ar[1]/$tm[-2]`;
        }
        my $out = "$ar[1]/$tm[-2]/$tm[-2].nochimera.fa";
        die "$! \n" unless open(T,">$out");
        my $index = 1;
        while(<I>)
        {   
            my $nam = "$tm[-2]\_$index";
            my $seq = <I>;
            print T ">$nam\n$seq";
            print NAL "$tm[-2]\t$nam\n";
            print TA ">$nam\n$seq";
            $index++;
        }
        close I;
        close T;
    }
    close REN;
    close TA;
    close NAL;

}#reanme_tags

1;
__END__
