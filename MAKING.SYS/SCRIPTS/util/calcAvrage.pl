#!/usr/bin/perl
use strict;
###############################
#prints avrages for all columns
#
################################
die "usage calcAvrage.pl File [max/false min/false colNum] \n" unless (scalar(@ARGV) == 1||scalar(@ARGV) == 4) ;
my $max = "false";
my $min = "false";
my $colNum = -1;
if (scalar(@ARGV) == 4){
    $max = $ARGV[1];
    $min = $ARGV[2];
    $colNum = $ARGV[3];
}


my $File = $ARGV[0];
open(FILE_HANDLE,$File) or die "Cannot open $File, $!";
my $line;
my @newLine;
my @sumLine;
my $lineNum = 0;
while($line= <FILE_HANDLE>){
    @newLine =split(/\s+/,$line);
    my $i= 0;
 
    foreach my $val( @newLine ){
        if ($colNum == -1 ||$colNum == $i){
            if((lc($max) eq "false" || $val<=$max) && (lc($min) eq "false" || $val>=$min) ){
                $sumLine[$i]+=$val;
            }
        }
        $i++;
    }
    if ($i>0){
    $lineNum++;}
}
print "$File\t";
if($colNum == -1){
    foreach my $val( @sumLine ){
            print (($val/$lineNum)."\t");
        }
}
else {
    print (($sumLine[$colNum]/$lineNum)."\t");
}
print "\n";