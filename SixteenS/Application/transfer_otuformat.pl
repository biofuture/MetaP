#!/usr/bin/perl -w
BEGIN{
    unshift (@INC,"/Users/jiangxiaotao/Documents/MetaP");
}
use strict;    
use SixteenS::Tools::OtuFormatTransfer;
if($#ARGV < 0){
    die "perl $0 <.fa> <.name> <tsc> <out>\n";
}
#die "perl $0";
#my @stotu;
#transfer_tscotu($ARGV[0],$ARGV[1],$ARGV[2],\@stotu);
transfer_tscotu($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
=head
print "MetaP Version 1.0\nOTU_ID\tRepresentTag\tRepresentSeq\tUniqueNumber\tTotalNumber\tComposition\n";
for my $key (@stotu){
    print "$key\n"; 
}
=cut
