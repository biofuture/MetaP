#!/usr/bin/perl -w
use strict;

use Getopt::Std;
#######################################
##BIPES quality control step for V6 data
#######################################




########################################
##CHIMERA remove
########################################


########################################
##Taxonomic assignment
########################################


########################################
##OTU picking steps
########################################

########################################
##Generate Otu table and taxon table
########################################


########################################
##-v verbose usage information for this script
########################################
my  $usage = <<USE;
    Author Xiao-tao JIANG
    Email: biofuture.jiang\@gmail.com
    perl $0 -a <1.fq> -b <2.fq> -c <primer.txt> -d <out dir> -q [quality control method] 






USE

print $usage;
#__END__
