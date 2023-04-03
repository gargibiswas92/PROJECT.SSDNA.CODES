#!/usr/bin/perl
use strict;
#################################
#will run bzip2 -v on all input files
#################################
my $user = `whoami`;
chomp $user;
foreach my $arg (@ARGV){
    my ($crdPath,$crdName)= $arg =~ /(.*\/)*(.*)/ ;
    `mv $arg /scratch/$user/`;
    `bzip2 -v /scratch/$user/$crdName`;
    `mv /scratch/$user/$crdName.bz2 $crdPath$crdName.bz2`;
}