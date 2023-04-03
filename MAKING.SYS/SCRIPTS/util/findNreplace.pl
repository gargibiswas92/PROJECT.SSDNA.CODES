#!/usr/bin/perl -w

use strict;
###############################################
#this will replace any replace word with its replacewith word
#the output will be in the outFile
#
#programer:Ohad Givaty givatyo@pob.huji.ac.il
################################################


if(scalar(@ARGV) < 4){
    die "usage: perl findNrplace.pl inFile outFile replaceWord replacewith [replaceWord2 replacewith2 replaceWord3 replacewith3 ...] \n";
}
##print " $ARGV[2] $ARGV[3]\n"; ##@@

my $infile= $ARGV[0];
my $outfile= $ARGV[1];

open(IN,"$infile")|| die("Cannot Open File $infile  for reading");


if($infile eq $outfile){
    my @lines  =<IN>;
    close(IN);

open(OUT,">$outfile") || die("Cannot Open File $outfile for writing");
my $first = 0;
foreach my $line (@lines ){
    for(my $i = 2; $i + 1 < scalar(@ARGV) ;$i = $i +2)
    {
	#for some strange reason in the replace with argument is \n we insert \n into the file and not newline !! the following condition solves this problem
	if($first ==0){	    
	    if($ARGV[$i+1] eq "\\n")
	    {$ARGV[$i+1]  ="\n";}	          
	}
	my $replaceWord = $ARGV[$i];
	my $replacewith = $ARGV[$i+1];	

	$line =~ s/$replaceWord/$replacewith/g;
    }
    $first++;
    print OUT $line;
}

close(OUT);
}
else {
open(OUT,">$outfile") || die("Cannot Open File $outfile for writing");
my $first = 0;
while ( my $line = <IN>){
    for(my $i = 2; $i + 1 < scalar(@ARGV) ;$i = $i +2)
    {
	#for some strange reason in the replace with argument is \n we insert \n into the file and not newline !! the following condition solves this problem
	if($first ==0){	    
	    if($ARGV[$i+1] eq "\\n")
	    {$ARGV[$i+1]  ="\n";}	          
	}
	my $replaceWord = $ARGV[$i];
	my $replacewith = $ARGV[$i+1];	    
	$line =~ s/$replaceWord/$replacewith/g;
    }
    $first++;
    print OUT $line;
}

close(OUT);
}
