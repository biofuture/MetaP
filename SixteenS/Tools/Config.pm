package SixteenS::Tools::Config;
use strict;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_config);

##This package supply a function to go through the configure file
##to extract some basic information for MetaP install and running

sub read_config{
    ##The input is the Config file
    ##the parameter is a hash to store all the variables

    die "Config file not exist $!\n" unless open(I,"../Config");
    while(<I>)
    {
        next if(/^#/);
        my @tem = split("=",$_);
        $$varible_hash{$tem[0]} = $tem[1]};
    }
    close I;
}#read_config
1;
__END__
