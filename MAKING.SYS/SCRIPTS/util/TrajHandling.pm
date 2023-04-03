package TrajHandling;
use strict;
use warnings;
use Math::Trig;

use base 'Exporter';
our @EXPORT_OK = qw(readBeadsData readResiduesAndBeadsData generatePDBResiduesArray readFrame skipFrame);

sub readBeadsData
{
    my ($trajFileHandleRef) = @_;
    
    my $currentLine = <$trajFileHandleRef>;
    my ($numberOfChains)  = $currentLine =~ m/(\d+)/;
    my $numberOfBeads = 0;
    for (my $chainIter = 0; $chainIter < $numberOfChains; $chainIter++)
    {
        $currentLine = <$trajFileHandleRef>;
        my ($numberOfBeadsInChain)  = $currentLine =~ m/(\d+)/;
        $numberOfBeads += $numberOfBeadsInChain;
    }
    my $beadsRef = [];
    my $currentResidueType = "@@@";
    my $residueIter = 0;
    my $beadTypesRef = {};
    for (my $beadIter = 0; $beadIter < $numberOfBeads; $beadIter++)
    {
      my $beadRef = {};
      $currentLine = <$trajFileHandleRef>;
      my ($beadType, $residueType) = $currentLine =~ m/\s*\d+\s+(\w+)\s+(\w+)/;
      if (exists($beadTypesRef->{$beadType}) or ($currentResidueType ne $residueType))
      {
        foreach my $beadTypeKey (keys(%$beadTypesRef))
        {
          delete $beadTypesRef->{$beadTypeKey} ;
        }
                              
        $residueIter++;
        $currentResidueType = $residueType;
      }
      
      $beadRef->{"TYPE"} = $beadType;
      $beadRef->{"RESIDUE_TYPE"} = $residueType;
      $beadRef->{"RESIDUE_ID"} = $residueIter;
      $beadTypesRef->{$beadType}++;
      $beadsRef->[$beadIter] = $beadRef;
    }
    
    return $beadsRef;
}

sub readResiduesAndBeadsData
{
    my ($trajFileHandleRef) = @_;
    
    my $currentLine = <$trajFileHandleRef>;
    my ($numberOfChains)  = $currentLine =~ m/(\d+)/;
    my $numberOfBeads = 0;
    for (my $chainIter = 0; $chainIter < $numberOfChains; $chainIter++)
    {
        $currentLine = <$trajFileHandleRef>;
        my ($numberOfBeadsInChain)  = $currentLine =~ m/(\d+)/;
        $numberOfBeads += $numberOfBeadsInChain;
    }
    
    my $beadsRef = [];
    my $residuesRef = [];
    my $currentResidueType = "@@@";
    my $residueIter = -1;
    my $beadTypesRef = {};
    my $residueRef;
    for (my $beadIter = 0; $beadIter < $numberOfBeads; $beadIter++)
    {
      $currentLine = <$trajFileHandleRef>;
      my ($beadType, $residueType) = $currentLine =~ m/\s+\d+\s+(\w+)\s+(\w+)/;
      if (exists($beadTypesRef->{$beadType}) or ($currentResidueType ne $residueType))
      {
        $residueIter++;

        foreach my $beadTypeKey (keys(%$beadTypesRef))
        {
          delete $beadTypesRef->{$beadTypeKey} ;
        }
                              
        $currentResidueType = $residueType;
        $residueRef = {};
      }
      
      my $beadRef = {};
      $beadRef->{"TYPE"} = $beadType;
      $beadRef->{"RESIDUE_TYPE"} = $residueType;
      $beadRef->{"RESIDUE_ID"} = $residueIter;
      $beadRef->{"ID"} = $beadIter;
      $beadTypesRef->{$beadType}++;
      $beadsRef->[$beadIter] = $beadRef;
      $residueRef->{"BEADS"}->{$beadType} = $beadRef;
      $residueRef->{"TYPE"} = $residueType;
      $residueRef->{"ID"} = $residueIter;
      $residuesRef->[$residueIter] = $residueRef;
    }
    
    return ($residuesRef,$beadsRef);
}

sub readFrame
{
    my ($trajFileHandleRef,$beadsDataRef) = @_;
    my $lastFrameFound  = 0;
    
    my $currentLine = <$trajFileHandleRef>;
    if (!$currentLine)
    {
      $lastFrameFound = 1;
      return ($lastFrameFound,[]);	
    }
    if ($currentLine =~ m/continue/i or $currentLine =~ m/\d+/i)
    {
      $currentLine = <$trajFileHandleRef>;

      my $beadCoordinatesRef = [];
      for (my $beadIter = 0; $beadIter < scalar(@$beadsDataRef); $beadIter++)
      {
        $currentLine = <$trajFileHandleRef>;
        $currentLine =~ m/\s*(\-?\d+(\.\d+)?)\s*(\-?\d+(\.\d+)?)\s*(\-?\d+(\.\d+)?)/;
        my ($x, $y, $z) = ($1, $3, $5);
        $beadCoordinatesRef->[$beadIter]->{"X"} = $x;
        $beadCoordinatesRef->[$beadIter]->{"Y"} = $y;
        $beadCoordinatesRef->[$beadIter]->{"Z"} = $z;
      }
      return ($lastFrameFound,$beadCoordinatesRef);

    }
    else
    {
      $lastFrameFound = 1;
      return ($lastFrameFound,[]);	
    }
    
}

sub readFrameFaster
{
    my ($trajFileHandleRef,$beadsDataRef) = @_;
    my $lastFrameFound  = 0;
    
    my $currentLine = <$trajFileHandleRef>;
    if (!$currentLine)
    {
      $lastFrameFound = 1;
      return ($lastFrameFound,[]);	
    }
    if ($currentLine =~ m/continue/i or $currentLine =~ m/\d+/i)
    {
      $currentLine = <$trajFileHandleRef>;

      my $beadCoordinatesRef = [];
      for (my $beadIter = 0; $beadIter < scalar(@$beadsDataRef); $beadIter++)
      {
        $currentLine = <$trajFileHandleRef>;
        my @split = split(/\s+/,$currentLine);
        ($beadCoordinatesRef->[$beadIter*3],$beadCoordinatesRef->[$beadIter*3+1],$beadCoordinatesRef->[$beadIter*3+2]) =
        ($split[1],$split[2],$split[3]);
      }
      return ($lastFrameFound,$beadCoordinatesRef);

    }
    else
    {
      $lastFrameFound = 1;
      return ($lastFrameFound,[]);	
    }
    
}


sub generatePDBResiduesArray
{
  my ($beadsDataRef,$beadCoordinatesRef) = @_;
  my $pdbResiduesRef = [];
  
  for (my $beadIter =0; $beadIter < scalar(@$beadsDataRef); $beadIter++)
  {
    my $residueType = $beadsDataRef->[$beadIter]->{"RESIDUE_TYPE"};
    my $residueID = $beadsDataRef->[$beadIter]->{"RESIDUE_ID"};
    my $beadType = $beadsDataRef->[$beadIter]->{"TYPE"};
    my $x = $beadCoordinatesRef->[$beadIter]->{"X"};
    my $y = $beadCoordinatesRef->[$beadIter]->{"Y"};
    my $z = $beadCoordinatesRef->[$beadIter]->{"Z"};
    
    $pdbResiduesRef->[$residueID-1]->{"RESIDUE_TYPE"} = $residueType;
    $pdbResiduesRef->[$residueID-1]->{"ATOMS_BY_TYPE"}->{$beadType}->{"X"} = $x;
    $pdbResiduesRef->[$residueID-1]->{"ATOMS_BY_TYPE"}->{$beadType}->{"Y"} = $y;
    $pdbResiduesRef->[$residueID-1]->{"ATOMS_BY_TYPE"}->{$beadType}->{"Z"} = $z;
  }

  return $pdbResiduesRef;
}

sub writeCoordinatesToPDBFile
{
  my ($trajFileHandleRef,$beadsDataRef,$executionPreferencesRef) = @_;
    
  my ($trajFilePrefix) = $executionPreferencesRef->{"TRAJ_FILE"} =~ m/(.*)\.dat/g;
  open(PDB_FILE_HANDLE, ">",$trajFilePrefix.".pdb") or die "Error opening pdb file for writing: $!\n";

  my $frameIter = 0;
  while (my $currentLine = <$trajFileHandleRef>)
  {
    $frameIter++;
    my $frameIsValid = 0;
    $frameIsValid = 1 if (($frameIter >= $executionPreferencesRef->{"FIRST_VALID_FRAME"})
    and ($frameIter <= $executionPreferencesRef->{"LAST_VALID_FRAME"})
    and ($frameIter % $executionPreferencesRef->{"SKIP_EVERY_X_FRAME"} == 0));
    printf PDB_FILE_HANDLE ("MODEL%8d\n",$frameIter) if $frameIsValid;
    for (my $beadIter =0; $beadIter < scalar(@$beadsDataRef); $beadIter++)
    {
        $currentLine = <$trajFileHandleRef>;
        chomp($currentLine);
        printf PDB_FILE_HANDLE ("ATOM%7d%4s%5s%6d%28s%6.2f%6.2f\n",
                               $beadIter+1,
                               $beadsDataRef->[$beadIter]->{"TYPE"},
                               $beadsDataRef->[$beadIter]->{"RESIDUE_TYPE"},
                               $beadsDataRef->[$beadIter]->{"RESIDUE_ID"},
                               $currentLine, '1', '0') if $frameIsValid;
    }
    print PDB_FILE_HANDLE ("ENDMDL\n") if $frameIsValid;
    $currentLine = <$trajFileHandleRef>;
  }
  
  close(PDB_FILE_HANDLE);
}

sub skipFrame
{
    my ($trajFileHandleRef,$beadsDataRef) = @_;
    my $lastFrameFound  = 0;
    
    my $currentLine = <$trajFileHandleRef>;
    if (!$currentLine)
    {
      $lastFrameFound = 1;
      return ($lastFrameFound,[]);	
    }
    if ($currentLine =~ m/continue/i or $currentLine =~ m/\d+/i)
    {
      $currentLine = <$trajFileHandleRef>;

      my $beadCoordinatesRef = [];
      for (my $beadIter = 0; $beadIter < scalar(@$beadsDataRef); $beadIter++)
      {
        <$trajFileHandleRef>;
      }
      return ($lastFrameFound,[]);
    }
    else
    {
      $lastFrameFound = 1;
      return ($lastFrameFound);	
    }
    
}
1;
