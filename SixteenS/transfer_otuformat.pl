#!/usr/bin/perl -w
BEGIN{
    unshift (@INC,"/Users/jiangxiaotao/Documents/MetaP");
}
use strict;    
use SixteenS::Tools::OtuFormatTransfer;
#die "perl $0";
my @stotu;
transfer_tscotu($ARGV[0],$ARGV[1],$ARGV[2],\@stotu);
print "MetaP Version 1.0\nOTU_ID\tRepresentTag\tRepresentSeq\tUniqueNumber\tTotalNumber\tComposition\n";
for my $key (@stotu){
    print "$key\n"; 
}
