#!/usr/bin/perl -w
use strict;

die "usage filter.pl filterValue filterFile fileToFilter output\n" unless (scalar(@ARGV) == 4) ;

my $filterValue = $ARGV[0];
my $zipped1 =0;
if ($ARGV[1] =~ /.bz2$/){
    print `/home/ohad/scripts/bunzip2l.pl  $ARGV[1]`;
    $ARGV[1] =~ s/\.bz2$//g;
    $zipped1 =1;
}
my $zipped2 =0;
if ($ARGV[2] =~ /.bz2$/){
    print `/home/ohad/scripts/bunzip2l.pl  $ARGV[2]`;
    $ARGV[2] =~ s/\.bz2$//g;
    $zipped2 =1;
}

open(FILTER,"$ARGV[1]")|| die("Cannot Open File $ARGV[1] for reading");
open(TO_FILTER,"$ARGV[2]")|| die("Cannot Open File $ARGV[2] for reading");
open(OUT,">$ARGV[3]")|| die("Cannot Open File $ARGV[3] for reading");

while (my $filter = <FILTER>){
    chomp($filter);
    my $to_filter = <TO_FILTER>;
    if ($filter == $filterValue){
        print OUT $to_filter;
    }                              
}
if ($zipped1==1){print `/home/ohad/scripts/bzip2l.pl  $ARGV[1]`;}
if ($zipped2==1){print `/home/ohad/scripts/bzip2l.pl  $ARGV[2]`;}


