#!/usr/bin/perl
use strict;
###############################
#prints the value which was most populated in the input file
#if a outputHistogramFile is inserted as input the hiltogram will be printed to it
################################
die "usage calcAvrage.pl File colNum interval [outputHistogramFile summingBy]  \n adding summingBy will sum column summingBy bined by colNum" unless (scalar(@ARGV) == 3||scalar(@ARGV) == 4||scalar(@ARGV) == 5) ;
my $File = $ARGV[0];
my $colNum = $ARGV[1];
my $interval = $ARGV[2];
my $biningBy=-1;
my $outputHistogramFile;
if (scalar(@ARGV) >= 4){
    $outputHistogramFile= $ARGV[3];
}
if (scalar(@ARGV) == 5){
    $biningBy= $ARGV[4];
}
open(FILE_HANDLE,$File) or die "Cannot open $File, $!";
my $line;
my @newLine;
my @histoUp;
my @histoDown;
my @histoUpBB;
my @histoDownBB;
my $currentCell;
my $val;
#my $lineNum = 0;
my $startVal=0;
$line= <FILE_HANDLE>;
@newLine =split(/\s+/,$line);
$startVal = $newLine[$colNum];
$startVal = int($startVal/$interval)*$interval;
$histoUp[0]++;
#make histogram
do{
    @newLine =split(/\s+/,$line);
    $val = $newLine[$colNum];
    #I add 0.00000001 so that presision mistakes will be rouned to the full value
    $currentCell = (($val - $startVal)/$interval +0.00000001);
    if($currentCell>=0){
        $histoUp[int($currentCell)]++;
        if ($biningBy>= 0){
            $histoUpBB[int($currentCell)]+= $newLine[$biningBy];
        }
    }
    else{
        $histoDown[int(-$currentCell)]++;
        if ($biningBy>= 0){
            $histoDownBB[int(-$currentCell)]+= $newLine[$biningBy];}
    }
#    $lineNum++;
}while($line= <FILE_HANDLE>);


# find most probable cell
my $maxProb = 0;
my $maxProbCell = 0;
for (my $i = 0;$i<scalar(@histoUp);$i++){
    if($histoUp[$i]>$maxProb){
        $maxProbCell= $i+1;
        $maxProb= $histoUp[$i];
    }
}
for (my $i = 0;$i<scalar(@histoDown);$i++){
    if($histoDown[$i]>$maxProb){
        $maxProbCell= -$i;
        $maxProb= $histoDown[$i];
    }
}
#print "$File\t";
print ($maxProbCell*$interval+$startVal-$interval);
print "\n";
if ($outputHistogramFile ne ""){
    open(OUT,">$outputHistogramFile") or die "Cannot open $outputHistogramFile for writing, $!";
 
    #print histoDown
    for (my $j = scalar(@histoDown)-1;$j>=0;$j--){
        if ($biningBy< 0){
        print OUT ((-$j)*$interval+$startVal-$interval)."\t".$histoDown[$j];
        unless (defined $histoDown[$j]){
            print OUT 0;
        }
        }
        else {
            print OUT ((-$j)*$interval+$startVal-$interval)."\t";
            if(defined $histoDown[$j]){
                print OUT ($histoDownBB[$j]/$histoDown[$j]);
            }
            else {
                print OUT "undefined";
            }
        }
        print OUT "\n";  
    }
    
    #print histoUp
    for (my $j = 0;$j<scalar(@histoUp);$j++){
        if ($biningBy< 0){
        print OUT (($j+1)*$interval+$startVal-$interval)."\t".$histoUp[$j];
            unless (defined $histoUp[$j]){
                print OUT 0;
            }
        }
        else {
            print OUT (($j+1)*$interval+$startVal-$interval)."\t";
            if(defined $histoUp[$j]){
                print OUT ($histoUpBB[$j]/$histoUp[$j]);
            }
            else {
                print OUT "undefined";
            }
        }
        print OUT "\n";
    }
}