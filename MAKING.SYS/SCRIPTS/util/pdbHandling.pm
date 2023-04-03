sub readPdb {

    my ($pdbCode)= @_;
 
    
    unless (open(PDB,$pdbCode)) {print("cant open file\n"); die "cant open $pdbCode, $!";}

    my @lines=<PDB>;
    close PDB;
    return make2dinfo(@lines);
}
sub readNMRPdb {

    my ($pdbCode,$structNum)= @_;
    unless (open(PDB,$pdbCode)) {print("cant open file\n"); die "cant open $pdbCode, $!";}
    my $dropedStructs = 0;
    my $line = 1;
    while ($dropedStructs<$structNum-1 ){
      if( $line = <PDB>){
        if ($line=~/^ENDMDL/){
          $dropedStructs++;}}
      else{
        last;}
    }
    
    my @lines=<PDB>;
    close PDB;
    return make2dinfo(@lines);
}
sub make2dinfo{
    my @twoDinfo; 
    my @lines = @_;
    my @info ="";
    my $line;
    foreach $line (@lines)
      {#for each line
        #in order not to count residues in another structure (i.e. different structures in NMR)
        chomp($line);
        if ($line=~/^ENDMDL/){last;}
        else{
          if($line=~ /^ATOM/){
            if(length($line)<=46){
              print("The pdb file you iserted is invalid line too short:<br>\n$line\n<br> " );
              exit 0;
            }
            
            @info = split(/ +/,$line);
            
            $info[1]=substr($line,6,5);#atom NUMBER
            if(($info[1] =~/^\s*-?\d+.?\d*\s*$/) == 0){
              print("The pdb file you iserted is invalid $info[1] an invalid atom NUMBER" );
              exit 0;
            }
            $info[2]=substr($line,12,4);#atom type
            $info[3]=substr($line,16,1);#alternative location indicator
            #print("$info[2]\t");
            $info[4]=substr($line,17,3);#AA name
            #print("$info[4]\t");
            $info[5]=substr($line,21,1);#chain ID
            $info[6]=substr($line,22,4);#AA NUMBER
            #print("$info[6]\t");
            if(($info[6] =~/^\s*-?\d+.?\d*\s*$/) == 0){
              print("The pdb file you iserted is invalid $info[6] an invalid AA NUMBER" );
              exit 0;
            }
            $info[7]=substr($line,26,1);#code for insertion of residue
            $info[8]=substr($line,30,8);#x cordinate
            if(($info[8] =~/^\s*-?\d+.?\d*\s*$/) == 0){
              print("The pdb file you iserted is invalid $info[8] an invalid x cordinate" );
              exit 0;
            }
            $info[9]=substr($line,38,8);#y cordinate
            if(($info[9] =~/^\s*-?\d+.?\d*\s*$/) == 0){
              print("The pdb file you iserted is invalid $info[9] an invalid y cordinate" );
              exit 0;
            }
            $info[10]=substr($line,46,8);#z cordinate
            if(($info[10] =~/^\s*-?\d+.?\d*\s*$/) == 0){
              print("The pdb file you iserted is invalid $info[10] an invalid z cordinate" );
              exit 0;
            }
            if(length($line)>=66){
              $info[12]=substr($line,60,6)
            }
            #there are 4 more fields I didn't code for
            $info[16]= $line;
            #print("$info[$pdbField{xCord}]\t$info[$pdbField{yCord}]\t$info[$pdbField{zCord}]\n");
            push (@twoDinfo ,[@info]);
          }
        }
      }
    return  @twoDinfo;
}

##################################################################################################
sub writePDB{
    my ($twoDinfoRef)= @_;
    my @twoDinfo = @$twoDinfoRef;
    
    foreach my $atom (@twoDinfo){
	my $line ="ATOM  ".fillLeft($$atom[1],5)." ".fillRight($$atom[2],4).fillRight($$atom[3],1).$$atom[4]." ".$$atom[5].fillLeft($$atom[6],4).$$atom[7]."   ".fillRight($$atom[8],8).fillRight($$atom[9],8).fillRight($$atom[10],8)."\n";
	print($line);
    }
}
sub writeNormPDB{
    my ($twoDinfoRef)= @_;
    my @twoDinfo = @$twoDinfoRef;
    my $AAnum = 0;
    my $AtomNum = 1;
    my $lastAA = "";
    foreach my $atom (@twoDinfo){
	if($lastAA ne $$atom[6]){
		$AAnum++;
	}
	my $line ="ATOM  ".fillLeft($AtomNum,5)." ".fillRight($$atom[2],4).fillRight($$atom[3],1).$$atom[4]." ".$$atom[5].fillLeft($AAnum,4).$$atom[7]."   ".fillRight($$atom[8],8).fillRight($$atom[9],8).fillRight($$atom[10],8)."\n";
	print($line);
	$AtomNum++;
	$lastAA = $$atom[6]
    }
}

sub fillLeft{
    my ($inputString,$wantedLength)=@_;
    while(length($inputString)< $wantedLength){
	$inputString = " ".$inputString;
    }
    return substr($inputString,0,$wantedLength ) ;
}
sub fillRight{
    my ($inputString,$wantedLength)=@_;
    while(length($inputString)< $wantedLength){
	$inputString = "$inputString ";
    }
    return  substr($inputString,0,$wantedLength ) ;
}
sub getPdbMap{
return (
                atomNum => 1,
                atomType => 2,
                chainID => 5,
                AAType => 4,
                AANum => 6,
                xCord => 8,
                yCord => 9,
                zCord => 10,
                tempFactor => 12,#score for stability
                pdbLine => 16,#the original pdb line
);
}

sub printCrd{
  my $numInLine =10;
  my ($twoDinfoRef,$name,$atomNum)= @_;
  my @twoDinfo = @$twoDinfoRef;
  my %pdbField = getPdbMap();
  my $inLine = 0;
  #print "$name\n ";
  foreach my $atom (@twoDinfo){
    if ($$atom[$pdbField{atomType}] =~ /CA/ ){
     print "$$atom[$pdbField{xCord}]";
     $inLine++;
      if($inLine ==$numInLine){ 
        print "\n"; 
       $inLine =0;
      }
     print "$$atom[$pdbField{yCord}]";
     $inLine++;
     if($inLine ==$numInLine){
       print "\n";
        $inLine =0;
     }
      print "$$atom[$pdbField{zCord}]"; 
      $inLine++;
     if($inLine ==$numInLine){
        print "\n"; 
        $inLine =0;
      }
    }
  }
  print "\n";
}

sub printPrmCrd{
  my ($twoDinfoRef,$name,$atomNum)= @_;
  my @twoDinfo = @$twoDinfoRef;
  my %pdbField = getPdbMap();
  my $inLine = 0;
  print " $name\n $atomNum\n ";
  foreach my $atom (@twoDinfo){
 #   if ($$atom[$pdbField{AAType}] =~ "CA" ){
      print " $$atom[$pdbField{xCord}] $$atom[$pdbField{yCord}] $$atom[$pdbField{yCord}]";
      $inLine +=3;
      if($inLine ==6){
        print "\n ";
	$inLine =0;
      }
  #  }
  }
  print "\n";
}
#pdb file must include only dna and protein as ATOM
#all phospates must have 4 oxigens (one from previus sugar)
#currently we assume that the atom numbebering and the AA numbering is the same at the begining since proteins are represented by CA only.
sub printGoPaulDna{
    my ($twoDinfoRef,$currentAAnum)= @_;
    my @pdbfile = @$twoDinfoRef;
    
    my %pdbField = getPdbMap();
  # not neede for paul  my @dnaNames = getDnaName();
    
    my $AAType;
   # my $atomType;
    my $origAAnum = -100;    
    my $currentAtomNum = $currentAAnum; #in the protein each beed is a AA
    $currentAAnum++;
    my $readyFlag = 0;
    my $dnaBeed = "non";
    my @phospate;
    my @sugar;
    my @base;
    my $firstInChain = 1;
    my $currentChain;
    my $lastChain= "OOOOO";
    my $curentAA;
    foreach my $atom (@pdbfile)  {
      $currentChain =$$atom[$pdbField{chainID}];
      $curentAA=$$atom[$pdbField{AANum}];   
      if($currentChain eq $lastChain){
         $firstInChain= 0;}
      else{
        $firstInChain= 1;
        $lastChain = $currentChain;
        }
        #SKIP protein atoms and hidrogen atoms
        if($$atom[$pdbField{AAType}] !~ / / ){
           next;
        }
        if($$atom[$pdbField{atomType}] =~ /^ *(O3(\*|\')|O5(\*|\')|O1P|O2P|P) *$/){
            push(@phospate,$atom);
        }
        else{
            if($$atom[$pdbField{atomType}] =~ /^ *(C3(\*|\')|C5(\*|\')|C4(\*|\')|O4(\*|\')|C1(\*|\')|C2(\*|\')) *$/){
                push(@sugar,$atom);
            }
        
            else{
                push(@base,$atom);
            }
        }
      #the phospate is full if it has all 5 atoms OR if the first O3 is missing it needs only 4
      if(scalar(@phospate) == 5|| (scalar(@phospate) == 4 && $phospate[0][$pdbField{atomType}] !~ /O3(\*|\')/  )){
            my ($x,$y,$z) =calcCenter(@phospate);
            $currentAtomNum++;

            $AAType = $$atom[$pdbField{AAType}]; 
            my $atomType = "P ";


            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($$atom[3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4)  .$$atom[7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
            print  (fillLeft($currentAtomNum,5).fillLeft($currentAtomNum,4).fillLeft($atomType,3).fillLeft($AAType,4).fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)." 99999.9\n");
#	    print "for atom  AA :$phospate[0][$pdbField{AANum}] atoms:  ";#@@@
#            foreach my $atom1 (@phospate){#@@@
#		print " $$atom1[$pdbField{atomType}]";#@@atomNum@
#	    }#@@@
#	    print "\n";#@@@
	    @phospate=();
        }
        if(scalar(@sugar) == 6){
            my ($x,$y,$z) =calcCenter(@sugar);
            $currentAtomNum++;
            my $atomType= "S ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($$atom[3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$$atom[7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
            print  (fillLeft($currentAtomNum,5).fillLeft($currentAtomNum,4).fillLeft($atomType,3).fillLeft($AAType,4).fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)." 99999.9\n");
#	    print "for atom  AA :$sugar[0][$pdbField{AANum}] atoms:  ";#@@@
#            foreach my $atom1 (@sugar){#@@@
#		print " $$atom1[$pdbField{atomType}]";#@@atomNum@
#	    }#@@@
#	    print "\n";#@@@
	    @sugar =();
        }
        #if it is a new nucleotide
        if(scalar(@base) > 5 && $curentAA ne $base[0][$pdbField{AANum}]){
            my ($x,$y,$z) =calcCenter(@base);
            $currentAtomNum++;
            my $atomType= "B ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($base[0][3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$base[0][7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
            print  (fillLeft($currentAtomNum,5).fillLeft($currentAtomNum,4).fillLeft($atomType,3).fillLeft($AAType,4).fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)." 99999.9\n");
            $currentAAnum++;
#	    print "for atom  AA :$base[0][$pdbField{AANum}] atoms:  ";#@@@
#            foreach my $atom1 (@base){#@@@
#		print " $$atom1[$pdbField{atomType}]";#@@atomNum@
#	    }#@@@
#	    print "\n";#@@@
	    @base = ();
        }
  }#end foreach 
      if(scalar(@base) > 0){
      my ($x,$y,$z) =calcCenter(@base);
            $currentAtomNum++;
            my $atomType= "B ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($base[0][3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$base[0][7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
            print  (fillLeft($currentAtomNum,5).fillLeft($currentAtomNum,4).fillLeft($atomType,3).fillLeft($AAType,4).fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)." 99999.9\n");
            #|   59  59 CA SER  21.880   8.303 -18.148   1.000|
  }
      return ($currentAtomNum,$currentAAnum);
      
}#end printGoPaulDna

###############################################
#prints only P,B,S for phosphate Base and sugar
#in a standard pdb format
###############################################
sub printGoNormDna{
    my ($twoDinfoRef,$currentAAnum)= @_;
    my @pdbfile = @$twoDinfoRef;
    
    my %pdbField = getPdbMap();
  # not neede for paul  my @dnaNames = getDnaName();
    
    my $AAType;
   # my $atomType;
    my $origAAnum = -100;    
    my $currentAtomNum = $currentAAnum; #in the protein each beed is a AA
    $currentAAnum++;
    my $readyFlag = 0;
    my $dnaBeed = "non";
    my @phospate;
    my @sugar;
    my @base;
    my $firstInChain = 1;
    my $currentChain;
    my $lastChain= "OOOOO";
    my $curentAA;
    foreach my $atom (@pdbfile)  {
      $currentChain =$$atom[$pdbField{chainID}];
      $curentAA=$$atom[$pdbField{AANum}];   
      if($currentChain eq $lastChain){
         $firstInChain= 0;}
      else{
        $firstInChain= 1;
        $lastChain = $currentChain;
        }
        #SKIP protein atoms and hidrogen atoms
        if($$atom[$pdbField{AAType}] !~ / / ){
           next;
        }
        if($$atom[$pdbField{atomType}] =~ /^ *(O3(\*|\')|O5(\*|\')|O1P|O2P|P) *$/){
            push(@phospate,$atom);
        }
        else{
            if($$atom[$pdbField{atomType}] =~ /^ *(C3(\*|\')|C5(\*|\')|C4(\*|\')|O4(\*|\')|C1(\*|\')|C2(\*|\')) *$/){
                push(@sugar,$atom);
            }
        
            else{
                push(@base,$atom);
            }
        }
      #the phospate is full if it has all 5 atoms OR if the first O3 is missing it needs only 4
      if(scalar(@phospate) == 5|| (scalar(@phospate) == 4 && $phospate[0][$pdbField{atomType}] !~ /O3(\*|\')/  )){
            my ($x,$y,$z) =calcCenter(@phospate);
            $currentAtomNum++;

            $AAType = $$atom[$pdbField{AAType}]; 
            my $atomType = "P ";

           print "ATOM  ".fillLeft($currentAtomNum,5)." ".fillRight($atomType,4).fillRight($$atom[3],1).$AAType." ".$$atom[5].fillLeft($currentAAnum,4).$$atom[7]."   ".fillRight($x,8).fillRight($y,8).fillRight($z,8)."\n";

    
#	    print "for atom  AA :$phospate[0][$pdbField{AANum}] atoms:  ";#@@@
#            foreach my $atom1 (@phospate){#@@@
#		print " $$atom1[$pdbField{atomType}]";#@@atomNum@
#	    }#@@@
#	    print "\n";#@@@
	    @phospate=();
        }
        if(scalar(@sugar) == 6){
            my ($x,$y,$z) =calcCenter(@sugar);
            $currentAtomNum++;
            my $atomType= "S ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($$atom[3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$$atom[7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
          print "ATOM  ".fillLeft($currentAtomNum,5)." ".fillRight($atomType,4).fillRight($$atom[3],1).$AAType." ".$$atom[5].fillLeft($currentAAnum,4).$$atom[7]."   ".fillRight($x,8).fillRight($y,8).fillRight($z,8)."\n";
#	    print "for atom  AA :$sugar[0][$pdbField{AANum}] atoms:  ";#@@@
#            foreach my $atom1 (@sugar){#@@@
#		print " $$atom1[$pdbField{atomType}]";#@@atomNum@
#	    }#@@@
#	    print "\n";#@@@
	    @sugar =();
        }
        #if it is a new nucleotide
        if(scalar(@base) > 5 && $curentAA ne $base[0][$pdbField{AANum}]){
            my ($x,$y,$z) =calcCenter(@base);
            $currentAtomNum++;
            my $atomType= "B ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($base[0][3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$base[0][7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
          print "ATOM  ".fillLeft($currentAtomNum,5)." ".fillRight($atomType,4).fillRight($$atom[3],1).$AAType." ".$$atom[5].fillLeft($currentAAnum,4).$$atom[7]."   ".fillRight($x,8).fillRight($y,8).fillRight($z,8)."\n";
            $currentAAnum++;
#	    print "for atom  AA :$base[0][$pdbField{AANum}] atoms:  ";#@@@
#            foreach my $atom1 (@base){#@@@
#		print " $$atom1[$pdbField{atomType}]";#@@atomNum@
#	    }#@@@
#	    print "\n";#@@@
	    @base = ();
        }
  }#end foreach 
      if(scalar(@base) > 0){
          my ($x,$y,$z) =calcCenter(@base);
          my $atom = $pdbfile[scalar(@pdbfile)-1];
          $currentAtomNum++;
          my $atomType= "B ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($base[0][3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$base[0][7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
          print "ATOM  ".fillLeft($currentAtomNum,5)." ".fillRight($atomType,4).fillRight($$atom[3],1).$AAType." ".$$atom[5].fillLeft($currentAAnum,4).$$atom[7]."   ".fillRight($x,8).fillRight($y,8).fillRight($z,8)."\n";
            #|   59  59 CA SER  21.880   8.303 -18.148   1.000|
  }
      return ($currentAtomNum,$currentAAnum);
      
}#end printGoNormDna

1;
