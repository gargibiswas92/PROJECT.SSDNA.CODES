#!/usr/bin/perl
use strict;
#################################
#will run bzip2 -v on all input files
#################################
my $user = `whoami`;
chomp $user;
foreach my $arg (@ARGV){
    my ($crdPath,$crdName)=  $arg  =~ /(.*\/)*(.*)/ ;
    `mv  $arg /scratch/$user/`;
    `bunzip2 -v /scratch/$user/$crdName`;
    $crdName =~ s/\.bz2//g;
    `mv /scratch/$user/$crdName $crdPath$crdName`;
}