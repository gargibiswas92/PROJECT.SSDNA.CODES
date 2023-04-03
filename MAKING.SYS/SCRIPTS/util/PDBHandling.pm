package PDBHandling;
use strict;
use warnings;
use Math::Trig;
use lib ("/home_b/gargi/SCRIPTS/MD/util/");
use MDMath;

use base 'Exporter';
our @EXPORT_OK = qw(readPDBFile readTrajPDBFrame writeChainAsPDBFile writeAllChainsAsPDBFile calculateChirals calculateImpropers calculateAngles calculateBeadDistances calculateSuccessiveDihedrals calculateVectors calculateNormalizedVectors);

my $atomMassByType = { "H" => 1,
                       "C" => 12,
                       "N" => 14,
                       "O" => 16,
                       "P" => 31,
                       "S" => 32};

sub readTrajPDBFrame
{
  my $fileHandleRef = $_[0];
  my $pdbAtomsRef = generateAtomsArray($_[0]);
  my ($pdbResiduesRef,$pdbResiduesByIDRef) = generateResiduesArray($pdbAtomsRef);
  my ($pdbChainsRef,$pdbChainsByIDRef) = generateChainsArray($pdbAtomsRef,$pdbResiduesRef);
  
  my $pdbDataRef = { "ATOMS" => $pdbAtomsRef,
                     "RESIDUES" => $pdbResiduesRef,
		     "RESIDUES_BY_ID" => $pdbResiduesByIDRef,
                     "CHAINS" => $pdbChainsRef,
		     "CHAINS_BY_ID" => $pdbChainsByIDRef};
  
  return $pdbDataRef;

}

sub readPDBFile
{
  my $pdbFile = $_[0];
  open (PDB_FILE_HANDLE,$pdbFile) or die "Cannot open PDB file $pdbFile: $!.\n";
  my $pdbAtomsRef = generateAtomsArray(\*PDB_FILE_HANDLE);
  close (PDB_FILE_HANDLE);  
  my ($pdbResiduesRef,$pdbResiduesByIDRef) = generateResiduesArray($pdbAtomsRef);
  my ($pdbChainsRef,$pdbChainsByIDRef) = generateChainsArray($pdbAtomsRef,$pdbResiduesRef);
  
  my $pdbDataRef = { "ATOMS" => $pdbAtomsRef,
                     "RESIDUES" => $pdbResiduesRef,
		     "RESIDUES_BY_ID" => $pdbResiduesByIDRef,
                     "CHAINS" => $pdbChainsRef,
		     "CHAINS_BY_ID" => $pdbChainsByIDRef};
  
  return $pdbDataRef;
}

sub generateChainsArray
{
  my ($pdbAtomsRef,$pdbResiduesRef) = @_;
  
  my $chainsRef = [];
  
  my $currentChainID = "@";
  my $chainIndex = -1;
  
  for(my $pdbResidueIter = 0; $pdbResidueIter < scalar(@$pdbResiduesRef); $pdbResidueIter++)
  {
    if ($pdbResiduesRef->[$pdbResidueIter]->{"CHAIN_ID"} ne $currentChainID)
    {
      $chainIndex++;
      $currentChainID = $pdbResiduesRef->[$pdbResidueIter]->{"CHAIN_ID"};
      $chainsRef->[$chainIndex]->{"CHAIN_ID"} = $pdbResiduesRef->[$pdbResidueIter]->{"CHAIN_ID"};
    }
    push(@{$chainsRef->[$chainIndex]->{"RESIDUES"}},$pdbResiduesRef->[$pdbResidueIter]);
  }

  $currentChainID = "@";
  $chainIndex = -1;
  
  for(my $pdbAtomIter = 0; $pdbAtomIter < scalar(@$pdbAtomsRef); $pdbAtomIter++)
  {
    if ($pdbAtomsRef->[$pdbAtomIter]->{"CHAIN_ID"} ne $currentChainID)
    {
      $chainIndex++;
      $currentChainID = $pdbAtomsRef->[$pdbAtomIter]->{"CHAIN_ID"};
      $chainsRef->[$chainIndex]->{"CHAIN_ID"} = $pdbAtomsRef->[$pdbAtomIter]->{"CHAIN_ID"};
    }
    push(@{$chainsRef->[$chainIndex]->{"ATOMS"}},$pdbAtomsRef->[$pdbAtomIter]);
  }


  my $chainsByIDRef = {};

  for(my $pdbResidueIter = 0; $pdbResidueIter < scalar(@$pdbResiduesRef); $pdbResidueIter++)
  {
    my $chainID = $pdbResiduesRef->[$pdbResidueIter]->{"CHAIN_ID"};
    push(@{$chainsByIDRef->{$chainID}->{"RESIDUES"}},$pdbResiduesRef->[$pdbResidueIter]);
  }

  for(my $pdbAtomIter = 0; $pdbAtomIter < scalar(@$pdbAtomsRef); $pdbAtomIter++)
  {
    my $chainID = $pdbAtomsRef->[$pdbAtomIter]->{"CHAIN_ID"};
    push(@{$chainsByIDRef->{$chainID}->{"ATOMS"}},$pdbAtomsRef->[$pdbAtomIter]);
  }
  
  return ($chainsRef,$chainsByIDRef);
}

sub generateResiduesArray
{
  my $pdbAtomsRef = $_[0];
  
  my $residuesRef = [];
  
  my $currentResidueID = "@";
  my $residueIndex = -1;
  
  for(my $pdbAtomIter = 0; $pdbAtomIter < scalar(@$pdbAtomsRef); $pdbAtomIter++)
  {
    my $newResidueID = $pdbAtomsRef->[$pdbAtomIter]->{"CHAIN_ID"}.".".$pdbAtomsRef->[$pdbAtomIter]->{"RESIDUE_ID"};
    if ($newResidueID ne $currentResidueID)
    {
      $residueIndex++;
      $currentResidueID = $newResidueID;
      $residuesRef->[$residueIndex]->{"CHAIN_ID"} = $pdbAtomsRef->[$pdbAtomIter]->{"CHAIN_ID"};
      $residuesRef->[$residueIndex]->{"RESIDUE_ID"} = $pdbAtomsRef->[$pdbAtomIter]->{"RESIDUE_ID"};
      $residuesRef->[$residueIndex]->{"RESIDUE_TYPE"} = $pdbAtomsRef->[$pdbAtomIter]->{"RESIDUE_TYPE"};
    }
    my ($atomType) = $pdbAtomsRef->[$pdbAtomIter]->{"TYPE"} =~ m/(\w+)/;
    $residuesRef->[$residueIndex]->{"ATOMS_BY_TYPE"}->{$atomType} = $pdbAtomsRef->[$pdbAtomIter];
    push(@{$residuesRef->[$residueIndex]->{"ATOMS"}},$pdbAtomsRef->[$pdbAtomIter]);
    
  }

  my $residuesByIDRef = {};
  for(my $pdbAtomIter = 0; $pdbAtomIter < scalar(@$pdbAtomsRef); $pdbAtomIter++)
  {
    my $currentResidueID = $pdbAtomsRef->[$pdbAtomIter]->{"RESIDUE_ID"};
    $residuesByIDRef->{$currentResidueID}->{"CHAIN_ID"} = $pdbAtomsRef->[$pdbAtomIter]->{"CHAIN_ID"};
    $residuesByIDRef->{$currentResidueID}->{"RESIDUE_ID"} = $pdbAtomsRef->[$pdbAtomIter]->{"RESIDUE_ID"};
    $residuesByIDRef->{$currentResidueID}->{"RESIDUE_TYPE"} = $pdbAtomsRef->[$pdbAtomIter]->{"RESIDUE_TYPE"};
    my ($atomType) = $pdbAtomsRef->[$pdbAtomIter]->{"TYPE"} =~ m/(\w+)/;
    $residuesByIDRef->{$currentResidueID}->{"ATOMS_BY_TYPE"}->{$atomType} = $pdbAtomsRef->[$pdbAtomIter];
    push(@{$residuesByIDRef->{$currentResidueID}->{"ATOMS"}},$pdbAtomsRef->[$pdbAtomIter]);
  }
  
  return ($residuesRef,$residuesByIDRef);
}

sub fixAtomType
{
  my ($atomType) = @_;
  
  if (substr($atomType,0,1) =~ m/[A-Z]/)
  {
    $atomType = substr($atomType,3,1).substr($atomType,0,3);
  }

  return $atomType;  

}

sub revertAtomType
{
  my ($atomType) = @_;
  
  $atomType = substr($atomType,1,3).substr($atomType,0,1);
  
  return $atomType;  

}

sub generateAtomsArray
{
  my $pdbFileHandleRef = $_[0];
  my $atomsRef = [];

  my $pdbLine;
  for (my $pdbLineIter = 1;$pdbLine = <$pdbFileHandleRef>; $pdbLineIter++)
  {
    chomp($pdbLine);
    last if ($pdbLine=~/^ENDMDL/);
    
    if($pdbLine =~ /^ATOM/ and (substr($pdbLine,72,3) !~ m/WAT/))
    {
      die "PDB line $pdbLineIter is too short:\n$pdbLine\n" unless (length($pdbLine) > 46);
      my $atomRef = {};
      
      $atomRef->{"ID"} = substr($pdbLine,6,5);
      $atomRef->{"ID"} =~ s/\s+//;
      die "PDB line $pdbLineIter contains invalid atom number .".$atomRef->{"ID"}.".\n"
      unless ($atomRef->{"ID"} =~/^\d+$/);
      
      $atomRef->{"TYPE"}= fixAtomType(substr($pdbLine,12,4));

      $atomRef->{"ALTERNATIVE_LOCATION"}= substr($pdbLine,16,1);
      $atomRef->{"ALTERNATIVE_LOCATION"} =~ s/\s+//;
      
      $atomRef->{"RESIDUE_TYPE"}= substr($pdbLine,17,3);
      $atomRef->{"RESIDUE_TYPE"} =~ s/\s+//;
      
      $atomRef->{"CHAIN_ID"}= substr($pdbLine,21,1);
      if ($atomRef->{"CHAIN_ID"} =~ m/^\s{1}$/)
      {
        if (length($pdbLine) >= 73) 
        {
          $atomRef->{"CHAIN_ID"}= substr($pdbLine,72,1);
	}
      }
      $atomRef->{"CHAIN_ID"} =~ s/\s+//;
      
      $atomRef->{"RESIDUE_ID"}= substr($pdbLine,22,4);
      $atomRef->{"RESIDUE_ID"} =~ s/\s+//;
      die "PDB line $pdbLineIter contains invalid residue number ".$atomRef->{"RESIDUE_ID"}.".\n"
      unless ($atomRef->{"RESIDUE_ID"} =~/^\d+$/);
      
      $atomRef->{"INSERTION_CODE"}= substr($pdbLine,26,1);
      $atomRef->{"INSERTION_CODE"} =~ s/\s+//;
      
      $atomRef->{"X"}= substr($pdbLine,30,8);
      $atomRef->{"X"} =~ s/\s+//;
      die "PDB line $pdbLineIter contains invalid x coordinate .".$atomRef->{"X"}.".\n"
      unless ($atomRef->{"X"} =~/^-?\d+(\.\d+)?$/);
      $atomRef->{"Y"}= substr($pdbLine,38,8);
      $atomRef->{"Y"} =~ s/\s+//;
      die "PDB line $pdbLineIter contains invalid y coordinate .".$atomRef->{"Y"}.".\n"
      unless ($atomRef->{"Y"} =~/^-?\d+(\.\d+)?$/);
      $atomRef->{"Z"}= substr($pdbLine,46,8);
      $atomRef->{"Z"} =~ s/\s+//;
      die "PDB line $pdbLineIter contains invalid z coordinate .".$atomRef->{"Z"}.".\n"
      unless ($atomRef->{"Z"} =~/^-?\d+(\.\d+)?$/);
      
      my $atomType = uc(substr($pdbLine,13,1));
      $atomType =~ s/\s+//;
      $atomRef->{"MASS"} = $atomMassByType->{$atomType};
      push(@$atomsRef,$atomRef);
    }
  }
  
  return $atomsRef;
}

sub fillLeft
{
    my ($inputString,$wantedLength)=@_;
    while(length($inputString)< $wantedLength)
    {
	$inputString = " ".$inputString;
    }
    return substr($inputString,0,$wantedLength);
}

sub fillRight
{
    my ($inputString,$wantedLength)=@_;
    while(length($inputString)< $wantedLength)
    {
	$inputString = "$inputString ";
    }
    return substr($inputString,0,$wantedLength );
}

sub writeChainAsPDBFile
{
  my ($chainRef,$outputFile) = @_;

  my $residueAtomTypesWritten;
  
  open(OUTPUT_FILE_HANDLE,">",$outputFile) or die "Cannot open temp file $outputFile:$!.\n";
  for (my $residueIter =0 ; $residueIter < scalar(@{$chainRef->{"RESIDUES"}}); $residueIter++)
  {
    my $residueRef = $chainRef->{"RESIDUES"}->[$residueIter];
    $residueAtomTypesWritten = {};
    for (my $atomIter =0; $atomIter < scalar(@{$residueRef->{"ATOMS"}}); $atomIter++)
    {
      my $atomRef = $residueRef->{"ATOMS"}->[$atomIter];
      my $type = $atomRef->{"TYPE"};
      $type =~ s/\s+//;
      next if (exists($residueAtomTypesWritten->{$type}));
      my $line ="ATOM  ".fillLeft($atomRef->{"ID"},5)." ".fillRight($atomRef->{"TYPE"},4)." ".
                         fillRight($atomRef->{"RESIDUE_TYPE"},3)."  ".fillLeft($atomRef->{"RESIDUE_ID"},4).
                         fillLeft($atomRef->{"X"},12).fillLeft($atomRef->{"Y"},8).fillLeft($atomRef->{"Z"},8)."\n";
      print OUTPUT_FILE_HANDLE $line;
      $residueAtomTypesWritten->{$type}++;
    }
  }
  print OUTPUT_FILE_HANDLE "END\n";
  close(OUTPUT_FILE_HANDLE);
}

sub writeAllChainsAsSingleChainPDBFile
{
  my ($pdbResiduesRef,$outputFile) = @_;
  my $residueAtomTypesWritten;
  open(OUTPUT_FILE_HANDLE,">",$outputFile) or die "Cannot open temp file $outputFile :$!.\n";
  for (my $residueIter =0 ; $residueIter < scalar(@$pdbResiduesRef); $residueIter++)
  {
    $residueAtomTypesWritten = {};
    my $residueRef = $pdbResiduesRef->[$residueIter];
    for (my $atomIter =0; $atomIter < scalar(@{$residueRef->{"ATOMS"}}); $atomIter++)
    {
 
      my $atomRef = $residueRef->{"ATOMS"}->[$atomIter];
      my $type = $atomRef->{"TYPE"};
      $type =~ s/\s+//;
      next if (exists($residueAtomTypesWritten->{$type}));
      my $line ="ATOM  ".fillLeft($atomRef->{"ID"},5)." ".fillRight($atomRef->{"TYPE"},4)." ".
                         fillRight($atomRef->{"RESIDUE_TYPE"},3)."  ".fillLeft($residueIter+1,4).
                         fillLeft($atomRef->{"X"},12).fillLeft($atomRef->{"Y"},8).fillLeft($atomRef->{"Z"},8)."\n";
      print OUTPUT_FILE_HANDLE $line;
      $residueAtomTypesWritten->{$type}++;
    }
  }
  print OUTPUT_FILE_HANDLE "END\n";
  close(OUTPUT_FILE_HANDLE);
}

sub writeAllChainsAsPDBFile
{
  my ($pdbResiduesRef,$outputFile) = @_;
  my $residueAtomTypesWritten;
  open(OUTPUT_FILE_HANDLE,">",$outputFile) or die "Cannot open temp file $outputFile :$!.\n";
  for (my $residueIter =0 ; $residueIter < scalar(@$pdbResiduesRef); $residueIter++)
  {
    $residueAtomTypesWritten = {};
    my $residueRef = $pdbResiduesRef->[$residueIter];
    for (my $atomIter =0; $atomIter < scalar(@{$residueRef->{"ATOMS"}}); $atomIter++)
    {
 
      my $atomRef = $residueRef->{"ATOMS"}->[$atomIter];
      my $type = $atomRef->{"TYPE"};
      $type =~ s/\s+//;
      next if (exists($residueAtomTypesWritten->{$type}));
      my $line ="ATOM  ".fillLeft($atomRef->{"ID"},5)." ".fillRight($atomRef->{"TYPE"},4)." ".
                         fillRight($atomRef->{"RESIDUE_TYPE"},3)." ".$atomRef->{"CHAIN_ID"}.fillLeft($residueIter+1,4).
                         fillLeft($atomRef->{"X"},12).fillLeft($atomRef->{"Y"},8).fillLeft($atomRef->{"Z"},8)."\n";
      print OUTPUT_FILE_HANDLE $line;
      $residueAtomTypesWritten->{$type}++;
    }
  }
  print OUTPUT_FILE_HANDLE "END\n";
  close(OUTPUT_FILE_HANDLE);
}



sub calculateImpropers
{
  my ($pdbResiduesRef) = @_;
  my $dihedralsRef = [];
  my $dihedralIter = -1;
  if ($pdbResiduesRef->[0]->{"RESIDUE_TYPE"} ne "GLY")
  {
    $dihedralIter++;
    $dihedralsRef->[$dihedralIter] = MDMath::calculateDihedral($pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                 $pdbResiduesRef->[1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                 $pdbResiduesRef->[2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                 $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CB"});
  }
  for (my $residueIter = 1; $residueIter < scalar(@$pdbResiduesRef)-1; $residueIter++)
  {
    next if ($pdbResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY");
    $dihedralIter++;
    $dihedralsRef->[$dihedralIter] = MDMath::calculateDihedral($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"});    
  }

  if ($pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"RESIDUE_TYPE"} ne "GLY")
  {
    $dihedralsRef->[$dihedralIter] = MDMath::calculateDihedral($pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-3]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CB"});
  }

  return $dihedralsRef;
}



sub calculateSuccessiveDihedrals
{
  my ($pdbResiduesRef) = @_;
  my $succesiveDihedralsRef = [];
  my $successiveDihedralIter = -1;
  for (my $residueIter = 0; $residueIter < scalar(@$pdbResiduesRef)-1; $residueIter++)
  {
    next if (($pdbResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY") or
             ($pdbResiduesRef->[$residueIter+1]->{"RESIDUE_TYPE"} eq "GLY"));
    $successiveDihedralIter++;
    $succesiveDihedralsRef->[$successiveDihedralIter] = MDMath::calculateDihedral($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"},
                                       $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CB"});    
  }
  return $succesiveDihedralsRef;
}

sub calculateDihedrals
{
  my ($pdbResiduesRef) = @_;
  my $dihedralsRef = [];
  my $dihedralIter = -1;
  for (my $residueIter = 0; $residueIter < scalar(@$pdbResiduesRef)-3; $residueIter++)
  {
    $dihedralIter++;
    $dihedralsRef->[$dihedralIter] = MDMath::calculateDihedral($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+3]->{"ATOMS_BY_TYPE"}->{"CA"});    
  }
  return $dihedralsRef;
}


sub calculateChirals
{
  my ($pdbResiduesRef) = @_;
  my $chiralsRef = [];
  my $coschiralsRef = [];
  my $tripleProductsRef = [];
  my $chiralIter = -1;
  if ($pdbResiduesRef->[0]->{"RESIDUE_TYPE"} ne "GLY")
  {
    $chiralIter++;
    ($tripleProductsRef->[$chiralIter],
     $coschiralsRef->[$chiralIter],
     $chiralsRef->[$chiralIter]) = MDMath::calculateChiral($pdbResiduesRef->[1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                   $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                   $pdbResiduesRef->[2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                   $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CB"});
  }
  for (my $residueIter = 1; $residueIter < scalar(@$pdbResiduesRef)-1; $residueIter++)
  {
    next if ($pdbResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY");
    $chiralIter++;
    ($tripleProductsRef->[$chiralIter],
     $coschiralsRef->[$chiralIter],
     $chiralsRef->[$chiralIter]) = MDMath::calculateChiral($pdbResiduesRef->[$residueIter-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"});    
  }

  if ($pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"RESIDUE_TYPE"} ne "GLY")
  {
    ($tripleProductsRef->[$chiralIter],
     $coschiralsRef->[$chiralIter],
     $chiralsRef->[$chiralIter]) = MDMath::calculateChiral($pdbResiduesRef->[scalar(@$pdbResiduesRef)-2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-3]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CB"});
  }

  return ($tripleProductsRef,$coschiralsRef,$chiralsRef);
}

sub calculateAngles
{
  my ($pdbResiduesRef) = @_;
  my $CAanglesRef = [];
  my $CBanglesRef = [];
  my $angleIter = -1;
  if ($pdbResiduesRef->[0]->{"RESIDUE_TYPE"} ne "GLY")
  {
    $angleIter++;
    $CAanglesRef->[$angleIter] = MDMath::calculateAngle($pdbResiduesRef->[1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                 $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                 $pdbResiduesRef->[2]->{"ATOMS_BY_TYPE"}->{"CA"});
    $CBanglesRef->[$angleIter] = MDMath::calculateAngle($pdbResiduesRef->[1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                 $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CA"},
                                                 $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CB"});
  }
  for (my $residueIter = 1; $residueIter < scalar(@$pdbResiduesRef)-1; $residueIter++)
  {
    next if ($pdbResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY");
    $angleIter++;
    $CAanglesRef->[$angleIter] = MDMath::calculateAngle($pdbResiduesRef->[$residueIter-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CA"});    
    $CBanglesRef->[$angleIter] = MDMath::calculateAngle($pdbResiduesRef->[$residueIter-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"});    
  }

  if ($pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"RESIDUE_TYPE"} ne "GLY")
  {
    $CAanglesRef->[$angleIter] = MDMath::calculateAngle($pdbResiduesRef->[scalar(@$pdbResiduesRef)-2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-3]->{"ATOMS_BY_TYPE"}->{"CA"});
    $CBanglesRef->[$angleIter] = MDMath::calculateAngle($pdbResiduesRef->[scalar(@$pdbResiduesRef)-2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                       $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CB"});

  }

  return ($CAanglesRef,$CBanglesRef);
}


sub calculateBeadDistances
{
  my ($pdbResiduesRef) = @_;
  my $beadDistancesRef = [];
  my $beadDistanceIter = -1;
  
  for (my $residueIter = 0; $residueIter < scalar(@$pdbResiduesRef); $residueIter++)
  {
    for (my $secondResidueIter = $residueIter+1; $secondResidueIter < scalar(@$pdbResiduesRef); $secondResidueIter++)
    {
      $beadDistanceIter++;
      $beadDistancesRef->[$beadDistanceIter]->{"TYPE1"} = "CA";
      $beadDistancesRef->[$beadDistanceIter]->{"ID1"} = $residueIter+1;
      $beadDistancesRef->[$beadDistanceIter]->{"TYPE2"} = "CA";
      $beadDistancesRef->[$beadDistanceIter]->{"ID2"} = $secondResidueIter+1;
      $beadDistancesRef->[$beadDistanceIter]->{"DISTANCE"} =
      MDMath::calculateDistance($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                            $pdbResiduesRef->[$secondResidueIter]->{"ATOMS_BY_TYPE"}->{"CA"});      
  
      my $CB1exists = exists($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"});
      my $CB2exists = exists($pdbResiduesRef->[$secondResidueIter]->{"ATOMS_BY_TYPE"}->{"CB"});
      if ($CB1exists)
      {
        $beadDistanceIter++;
        $beadDistancesRef->[$beadDistanceIter]->{"TYPE1"} = "CB";
        $beadDistancesRef->[$beadDistanceIter]->{"ID1"} = $residueIter+1;
        $beadDistancesRef->[$beadDistanceIter]->{"TYPE2"} = "CA";
        $beadDistancesRef->[$beadDistanceIter]->{"ID2"} = $secondResidueIter+1;
        $beadDistancesRef->[$beadDistanceIter]->{"DISTANCE"} =
        MDMath::calculateDistance($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"},
                              $pdbResiduesRef->[$secondResidueIter]->{"ATOMS_BY_TYPE"}->{"CA"});      
      }
      if ($CB2exists)
      {
        $beadDistanceIter++;
        $beadDistancesRef->[$beadDistanceIter]->{"TYPE1"} = "CA";
        $beadDistancesRef->[$beadDistanceIter]->{"ID1"} = $residueIter+1;
        $beadDistancesRef->[$beadDistanceIter]->{"TYPE2"} = "CB";
        $beadDistancesRef->[$beadDistanceIter]->{"ID2"} = $secondResidueIter+1;
        $beadDistancesRef->[$beadDistanceIter]->{"DISTANCE"} =
        MDMath::calculateDistance($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                              $pdbResiduesRef->[$secondResidueIter]->{"ATOMS_BY_TYPE"}->{"CB"});      
        }
      if ($CB1exists and $CB2exists)
      {
        $beadDistanceIter++;
        $beadDistancesRef->[$beadDistanceIter]->{"TYPE1"} = "CB";
        $beadDistancesRef->[$beadDistanceIter]->{"ID1"} = $residueIter+1;
        $beadDistancesRef->[$beadDistanceIter]->{"TYPE2"} = "CB";
        $beadDistancesRef->[$beadDistanceIter]->{"ID2"} = $secondResidueIter+1;
        $beadDistancesRef->[$beadDistanceIter]->{"DISTANCE"} =
        MDMath::calculateDistance($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"},
                              $pdbResiduesRef->[$secondResidueIter]->{"ATOMS_BY_TYPE"}->{"CB"});      
      }

    }      
  }

  return $beadDistancesRef;
}

sub calculateVectors
{
  my ($pdbResiduesRef) = @_;
  my $caVectorsRef = [];
  my $cbVectorsRef = [];
  my $vectorIter = -1;
  
  for (my $residueIter = 0; $residueIter < scalar(@$pdbResiduesRef)-1; $residueIter++)
  {
      $vectorIter++;
      $caVectorsRef->[$vectorIter]=MDMath::calculateVector($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                            $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CA"});      
      next if ($pdbResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY");
      $cbVectorsRef->[$vectorIter] = MDMath::calculateVector($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                            $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"})
  }
  return ($caVectorsRef,$cbVectorsRef);
}

sub calculateNormalizedVectors
{
  my ($pdbResiduesRef) = @_;
  my $caVectorsRef = [];
  my $cbVectorsRef = [];
  my $vectorIter = -1;
  
  for (my $residueIter = 3; $residueIter < scalar(@$pdbResiduesRef); $residueIter++)
  {
      $vectorIter++;
      my $residuePlaneTransitionMatrixRef = MDMath::calculatePlaneTransitionMatrix($pdbResiduesRef->[$residueIter-3]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[$residueIter-2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[$residueIter-1]->{"ATOMS_BY_TYPE"}->{"CA"});
      my $CAvector = MDMath::calculateVector($pdbResiduesRef->[$residueIter-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                            $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"});
      $caVectorsRef->[$vectorIter] = MDMath::multiplyMatrixByVector($residuePlaneTransitionMatrixRef,$CAvector);
      MDMath::scaleVector($caVectorsRef->[$vectorIter]);
  }
  
  $vectorIter = -1;
  if ($pdbResiduesRef->[0]->{"RESIDUE_TYPE"} ne "GLY")
  {
    $vectorIter++;
      my $residuePlaneTransitionMatrixRef = MDMath::calculatePlaneTransitionMatrix($pdbResiduesRef->[1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[2]->{"ATOMS_BY_TYPE"}->{"CA"});
      my $CBvector = MDMath::calculateVector($pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CA"},
                            $pdbResiduesRef->[0]->{"ATOMS_BY_TYPE"}->{"CB"});
      $cbVectorsRef->[$vectorIter] = MDMath::multiplyMatrixByVector($residuePlaneTransitionMatrixRef,$CBvector);
      MDMath::scaleVector($cbVectorsRef->[$vectorIter]);
  }
  for (my $residueIter = 1; $residueIter < scalar(@$pdbResiduesRef)-1; $residueIter++)
  {
      next if ($pdbResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY");
      $vectorIter++;
      my $residuePlaneTransitionMatrixRef = MDMath::calculatePlaneTransitionMatrix($pdbResiduesRef->[$residueIter-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[$residueIter+1]->{"ATOMS_BY_TYPE"}->{"CA"});
      my $CBvector = MDMath::calculateVector($pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CA"},
                            $pdbResiduesRef->[$residueIter]->{"ATOMS_BY_TYPE"}->{"CB"});
      $cbVectorsRef->[$vectorIter] = MDMath::multiplyMatrixByVector($residuePlaneTransitionMatrixRef,$CBvector);
      MDMath::scaleVector($cbVectorsRef->[$vectorIter]);
  }  
  if ($pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"RESIDUE_TYPE"} ne "GLY")
  {
    $vectorIter++;
      my $residuePlaneTransitionMatrixRef = MDMath::calculatePlaneTransitionMatrix($pdbResiduesRef->[scalar(@$pdbResiduesRef)-2]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                                     $pdbResiduesRef->[scalar(@$pdbResiduesRef)-3]->{"ATOMS_BY_TYPE"}->{"CA"});
      my $CBvector = MDMath::calculateVector($pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CA"},
                            $pdbResiduesRef->[scalar(@$pdbResiduesRef)-1]->{"ATOMS_BY_TYPE"}->{"CB"});
      $cbVectorsRef->[$vectorIter] = MDMath::multiplyMatrixByVector($residuePlaneTransitionMatrixRef,$CBvector);
      MDMath::scaleVector($cbVectorsRef->[$vectorIter]);
  } 
  return ($caVectorsRef,$cbVectorsRef);
}

1;
