package CA;
use strict;
use warnings;
use Math::Trig;
$| = 1;

#######################3
#
#
#########################
#correction 1.0 ohad givaty adding static DNA
use base 'Exporter';
our @EXPORT_OK = qw(createMDInputFile);

#constants
my $utilDirectory = "/home/arielaz/scripts/util/";
require $utilDirectory."MDMath.pm";

my $caBeadAtoms = { "C" => 1,
                    "CA" => 1,
                    "H" => 1,
                    "HA" => 1,
                    "N" => 1,
                    "O" => 1 } ;

my $residueCharges = { "LYS" => 1,
                       "ARG" => 1,
                       "GLU" => -1,
                       "ASP" => -1 };

my $beadRepulsionRadius = { "CA" => 2,
                            "CB" => 1.5} ;



my $minimumRepulsionDistance = { "CA_CA" => 4,
                                 "CA_CB" => 2,
                                 "CB_CA" => 2,
                                 "CB_CB" => 1} ;


sub createMDInputFile
{
    my ($executionPreferencesRef,$pdbDataRef) = @_;
    my $mdDataRef = {};
    determineBeads($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineDNABeads($executionPreferencesRef,$pdbDataRef,$mdDataRef);## 1.0
    determineBonds($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineAngles($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineDihedrals($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineContacts($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineRepulsions($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    #determineDNARepulsions($executionPreferencesRef,$pdbDataRef,$mdDataRef);## 1.0
    determineElectrostatics($executionPreferencesRef,$pdbDataRef,$mdDataRef) if ($executionPreferencesRef->{"ELECTROSTATICS"} =~ m/YES/i);
    #determineDNAElectrostatics($executionPreferencesRef,$pdbDataRef,$mdDataRef) if ($executionPreferencesRef->{"ELECTROSTATICS"} =~ m/YES/i);## 1.0
    
    printMDInputFile($executionPreferencesRef,$pdbDataRef,$mdDataRef);## 1.0 needs to be changed
}
sub determineDNABeads
{
    #my @validResTypes = ["C","A","G","T"];
    my ($executionPreferencesRef,$pdbDataRef,$mdDataRef)=@_;
    my $beadsDataRef = [];
    my $residuesDataRef = $mdDataRef->{"RESIDUES"};
    my $pdbfile = $pdbDataRef->{"ATOMS"} ;
    my @pdbfile = @$pdbfile;
    my $AAType;
    my $chainsDataRef = $mdDataRef->{"CHAINS"};
   # my $atomType;
    my $origAAnum = -100;
    my $currentRESNum = $mdDataRef->{"BEADS"}->[scalar(@{$mdDataRef->{"BEADS"}})-1]->{"ATOM_ID"};
    $currentRESNum++;
    my $currentBeadNum = $currentRESNum-1;
    my $readyFlag = 0;
    my $dnaBeed = "non";
    my @phospate;
    my @sugar;
    my @base;
    my $currentChain;
    my $lastChain= "OOOOO";
    my $curentRes;
    my $lastRes = -1000;
    my $currentResidueRef = {};
    my $currentChainRef ={};
    my $firstBeadInCain =$currentBeadNum+1;
    foreach my $atom (@pdbfile)  {
        
      
        $curentRes=$$atom{"RESIDUE_ID"};
        if(    $curentRes ne $lastRes ){
            if ($lastRes ne -1000){
                my $lastAtom =  $base[0];
                my($phospateBead,$baseBead,$sugarBead,$currentBeadNum1)
                =makeDNARes(\@phospate  ,\@base,\@sugar,$currentBeadNum,$lastAtom,$currentChain,$executionPreferencesRef);
                $currentBeadNum =$currentBeadNum1;
                $currentResidueRef->{"RESIDUE_ID"} =  $lastAtom->{"RESIDUE_ID"};
                $currentResidueRef->{"RESIDUE_TYPE"} = $lastAtom->{"RESIDUE_TYPE"};
                
                $currentResidueRef->{"BEADS_DNA"}->{"P"} = $phospateBead;
                push(@$beadsDataRef, $phospateBead);
                push(@{$currentChainRef->{"BEADS_DNA"}},$phospateBead);
                
                $currentResidueRef->{"BEADS_DNA"}->{"B"} = $baseBead;
                push(@{$currentChainRef->{"BEADS_DNA"}},$baseBead);
                push(@$beadsDataRef, ,$baseBead);
                
                $currentResidueRef->{"BEADS_DNA"}->{"S"} = $sugarBead;
                push(@{$currentChainRef->{"BEADS_DNA"}},$sugarBead);
                push(@$beadsDataRef, $sugarBead);
                
                push(@$residuesDataRef,$currentResidueRef);
                $currentChainRef->{"RESIDUES"}->{$lastAtom->{"RESIDUE_ID"}} = $currentResidueRef;
                $currentRESNum++;
                $currentResidueRef = {};
                
            }
            $lastRes = $curentRes;
        }
        $currentChain =$$atom{"CHAIN_ID"};
        
        
        if($currentChain ne $lastChain){
            $currentChainRef->{"SIZE"} = $currentBeadNum - $firstBeadInCain;
            if ($currentChainRef->{"SIZE"} >0){
                $chainsDataRef->{$lastChain} = $currentChainRef;
            }
            $currentChainRef =$chainsDataRef->{$currentChain};
            $lastChain = $currentChain;
            $firstBeadInCain = $currentBeadNum+1;
        }
        #SKIP protein atoms and hidrogen atoms
        if(length ($$atom{"RESIDUE_TYPE"}) == 3){  #alternative $$atom{"RESIDUE_TYPE"} !~ / / ){
           next;
        }
        if($$atom{"TYPE"} =~ /^ *(O3(\*|\')|O5(\*|\')|O1P|O2P|P) *$/){
            push(@phospate,$atom);
        }
        else{
            if($$atom{"TYPE"} =~ /^ *(C3(\*|\')|C5(\*|\')|C4(\*|\')|O4(\*|\')|C1(\*|\')|C2(\*|\')) *$/){
                push(@sugar,$atom);
            }
        
            else{
                push(@base,$atom);
            }
        }
        
  }#end foreach
    
    my $lastAtom =  $base[0];
            my($phospateBead,$baseBead,$sugarBead,$currentBeadNum1)=makeDNARes(\@phospate  ,\@base,\@sugar,$currentBeadNum,$lastAtom,$currentChain,$executionPreferencesRef);
            $currentBeadNum =$currentBeadNum1;
            $currentResidueRef->{"RESIDUE_ID"} =  $lastAtom->{"RESIDUE_ID"};
            $currentResidueRef->{"RESIDUE_TYPE"} = $lastAtom->{"RESIDUE_TYPE"};
                
            if (keys(%$phospateBead) != 0)
            {
                $currentResidueRef->{"BEADS_DNA"}->{"P"} = $phospateBead;
                push(@$beadsDataRef, $phospateBead);
                push(@{$currentChainRef->{"BEADS_DNA"}},$phospateBead);
            }
            if (keys(%$baseBead) != 0)
            {
                $currentResidueRef->{"BEADS_DNA"}->{"B"} = $baseBead;
                push(@{$currentChainRef->{"BEADS_DNA"}},$baseBead);
                push(@$beadsDataRef, ,$baseBead);
            }
            if (keys(%$sugarBead) != 0)
            {
                $currentResidueRef->{"BEADS_DNA"}->{"S"} = $sugarBead;
                push(@{$currentChainRef->{"BEADS_DNA"}},$sugarBead);
                push(@$beadsDataRef, $sugarBead);
            }    
                push(@$residuesDataRef,$currentResidueRef);
                $currentChainRef->{"RESIDUES"}->{$lastAtom->{"RESIDUE_ID"}} = $currentResidueRef;
                $currentRESNum++;
                $currentResidueRef = {};
                   
    if (exists($currentChainRef->{"BEADS_DNA"}) && scalar(@{$currentChainRef->{"BEADS_DNA"}}) > 0){
        $chainsDataRef->{$currentChain} = $currentChainRef;
    }
    $mdDataRef->{"BEADS_DNA"} = $beadsDataRef;
    
      
}

sub makeDNARes{
    
    my ( $phospate  ,$base,$sugar,$currentBeadNum,$atom,$currentChain,$executionPreferencesRef) =@_;
    my @phospate  =@$phospate;
    my @base =@$base;
    my @sugar  = @$sugar;
    my $sugarBead={};
    my $baseBead={};
    my $phospateBead={};
        
        
        
#phosphate
    if (scalar(@phospate) != 0)
    {
            $currentBeadNum++;

            #$AAType = $$atom{"TYPE"}; 
            my $atomType = "P ";
             my ($x,$y,$z) =calcCenter(@phospate);

            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($$atom[3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4)  .$$atom[7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
           ####################################3

            if ($executionPreferencesRef->{"NORMALIZE"} =~ m/YES/i)
            {
                $phospateBead->{"ATOM_ID"} = $currentBeadNum;
            }
            else
            {
                $phospateBead->{"ATOM_ID"} = $currentBeadNum;
            }
              $phospateBead->{"ATOM_TYPE"} = $atomType;
              $phospateBead->{"CHAIN_ID"} = $currentChain;
              $phospateBead->{"RESIDUE_ID"} = $atom->{"RESIDUE_ID"};
              $phospateBead->{"RESIDUE_TYPE"} = $atom->{"RESIDUE_TYPE"};
              $phospateBead->{"X"} = $x;
              $phospateBead->{"Y"} = $y;
              $phospateBead->{"Z"} = $z;
              $phospateBead->{"MASS"} = 999.0;
            #############################################
    }
#  sugar
    if (scalar(@phospate) != 0)
    {
            my ($x,$y,$z) =calcCenter(@sugar);
            $currentBeadNum++;
            my $atomType= "S ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($$atom[3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$$atom[7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
            
           ####################################3
            if ($executionPreferencesRef->{"NORMALIZE"} =~ m/YES/i)
            {
                $sugarBead->{"ATOM_ID"} = $currentBeadNum;
            }
            else
            {
                $sugarBead->{"ATOM_ID"} = $currentBeadNum;
            }
              $sugarBead->{"ATOM_TYPE"} = $atomType;
              $sugarBead->{"CHAIN_ID"} = $currentChain;
              $sugarBead->{"RESIDUE_ID"} = $atom->{"RESIDUE_ID"};
              $sugarBead->{"RESIDUE_TYPE"} = $atom->{"RESIDUE_TYPE"};
              $sugarBead->{"X"} = $x;
              $sugarBead->{"Y"} = $y;
              $sugarBead->{"Z"} = $z;
              $sugarBead->{"MASS"} = 999.0;

            ############################################
    }
#   base
    if (scalar(@phospate) != 0)
    {
        #if it is a new nucleotide
            my ($x,$y,$z) =calcCenter(@base);
            $currentBeadNum++;
            my $atomType= "B ";
            
            #print  ("ATOM  ".fillLeft($currentAtomNum,5)."  ".fillRight($atomType,3).fillRight($base[0][3],1).fillRight($AAType."-",3)." "." ".fillRight($currentAAnum,4).$base[0][7]."   ".fillLeft($x,8).fillLeft($y,8).fillLeft($z,8)."\n");
           ####################################3
            if ($executionPreferencesRef->{"NORMALIZE"} =~ m/YES/i)
            {
                $baseBead->{"ATOM_ID"} = $currentBeadNum;
            }
            else
            {
                $baseBead->{"ATOM_ID"} = $currentBeadNum;
            }
              $baseBead->{"ATOM_TYPE"} = $atomType;
              $baseBead->{"CHAIN_ID"} = $currentChain;
              $baseBead->{"RESIDUE_ID"} = $atom->{"RESIDUE_ID"};
              $baseBead->{"RESIDUE_TYPE"} = $atom->{"RESIDUE_TYPE"};
              $baseBead->{"X"} = $x;
              $baseBead->{"Y"} = $y;
              $baseBead->{"Z"} = $z;
              $baseBead->{"MASS"} = 999.0;
            ############################################
    }
            return ($phospateBead,$baseBead,$sugarBead,$currentBeadNum);
    
}







sub determineBeads
{
    print "Generating Beads.\n";
    my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
    
    my $beadsDataRef = [];
    my $residuesDataRef = [];
    my $residuesByIDDataRef = {};
    my $chainsDataRef = {};

    my $firstAtomID = ( exists($executionPreferencesRef->{"FIRST_VALID_ATOM"}) ?
                        $executionPreferencesRef->{"FIRST_VALID_ATOM"} : 
                        $pdbDataRef->{"ATOMS"}->[0]->{"ID"} ) ;
    my $lastAtomID = ( exists($executionPreferencesRef->{"LAST_VALID_ATOM"}) ?
                        $executionPreferencesRef->{"LAST_VALID_ATOM"} : 
                        $pdbDataRef->{"ATOMS"}->[scalar(@{$pdbDataRef->{"ATOMS"}})-1]->{"ID"} ) ;


    my $beadIter = -1;
    my $residueIter = -1;
    for (my $pdbChainIter = 0; $pdbChainIter < scalar(@{$pdbDataRef->{"CHAINS"}}) ; $pdbChainIter++)
    {
      my $currentChainRef = {};
      my $pdbChainRef = $pdbDataRef->{"CHAINS"}->[$pdbChainIter];
      my $currentChainID = $pdbChainRef->{"CHAIN_ID"};
      my $firstBeadIndex = $beadIter;
      for (my $pdbResidueIter = 0; $pdbResidueIter < scalar(@{$pdbChainRef->{"RESIDUES"}}) ; $pdbResidueIter++
)
      {
        $residueIter++;        
        my $pdbResidueRef = $pdbChainRef->{"RESIDUES"}->[$pdbResidueIter];
        my $currentResidueID = $pdbResidueRef->{"RESIDUE_ID"};
        my $pdbAtomsRef = $pdbResidueRef->{"ATOMS"};
        my $currentResidueRef = {};

        for (my $pdbAtomIter = 0; $pdbAtomIter < scalar(@{$pdbResidueRef->{"ATOMS"}}); $pdbAtomIter++)
        {
          my $pdbAtomRef = $pdbResidueRef->{"ATOMS"}->[$pdbAtomIter];
          my $currentAtomID = $pdbAtomRef->{"ID"};
          if (($currentAtomID < $firstAtomID) or ($currentAtomID > $lastAtomID))
          {
            $residueIter--;
            next;
          }


          if (my ($currentAtomType) = $pdbAtomRef->{"TYPE"} =~ m/(CA)/)
          {
              next if (exists($currentResidueRef->{"BEADS"}->{"CA"}));
              $beadIter++;
              my $beadRef = {};
              if ($executionPreferencesRef->{"NORMALIZE"} =~ m/YES/i)
              {
                $beadRef->{"ATOM_ID"} = $beadIter+1;
              }
              else
              {
                $beadRef->{"ATOM_ID"} = $pdbAtomRef->{"ID"};
              }
              $beadRef->{"ATOM_TYPE"} = $currentAtomType;
              $beadRef->{"CHAIN_ID"} = $currentChainID;
              $beadRef->{"RESIDUE_ID"} = $currentResidueID;
              $beadRef->{"RESIDUE_TYPE"} = $pdbAtomRef->{"RESIDUE_TYPE"};
              $beadRef->{"X"} = $pdbAtomRef->{"X"};
              $beadRef->{"Y"} = $pdbAtomRef->{"Y"};
              $beadRef->{"Z"} = $pdbAtomRef->{"Z"};
              $beadRef->{"MASS"} = 1.0;
              $beadsDataRef->[$beadIter] = $beadRef;
              $currentResidueRef->{"BEADS"}->{$currentAtomType} = $beadRef;
              $currentResidueRef->{"RESIDUE_ID"} = $currentResidueID;
              $currentResidueRef->{"RESIDUE_TYPE"} = $pdbAtomRef->{"RESIDUE_TYPE"};
              push(@{$currentChainRef->{"BEADS"}},$beadRef);
           }
         }
        if (scalar(keys(%{$currentResidueRef->{"BEADS"}})) > 0)
        {
          $residuesDataRef->[$residueIter] = $currentResidueRef;
          $residuesByIDDataRef->{$currentResidueID} = $currentResidueRef;
          $currentChainRef->{"RESIDUES"}->{$currentResidueID} = $currentResidueRef;            
        }
        else
        {
          $residueIter--;
        }
      }
      if (exists($currentChainRef->{"BEADS"}) &&  scalar(@{$currentChainRef->{"BEADS"}}) > 0){
        $chainsDataRef->{$currentChainID} = $currentChainRef;
      }
    }
    $mdDataRef->{"BEADS"} = $beadsDataRef;
    $mdDataRef->{"RESIDUES"} = $residuesDataRef;
    $mdDataRef->{"RESIDUES_BY_ID"} = $residuesByIDDataRef;
    $mdDataRef->{"CHAINS"} = $chainsDataRef;
}


sub determineBonds
{
  print "Calculating Bonds.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  my $bondsRef = [];

  my $mdResiduesRef = $mdDataRef->{"RESIDUES"};
  
  my $bondIter = -1;
  for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef)-1; $residueIter++)
  {  
    $bondIter++;
    my $CACAbondRef = {};
    my $currentCABeadRef = $mdResiduesRef->[$residueIter] ->{"BEADS"}->{"CA"};
    my $nextCABeadRef = $mdResiduesRef->[$residueIter+1] ->{"BEADS"}->{"CA"};
    $CACAbondRef->{"FIRST_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
    $CACAbondRef->{"SECOND_BEAD_ID"} = $nextCABeadRef->{"ATOM_ID"};
    $CACAbondRef->{"DISTANCE"} = MDMath::calculateDistance($currentCABeadRef,$nextCABeadRef);
    if ($currentCABeadRef->{"CHAIN_ID"} eq $nextCABeadRef->{"CHAIN_ID"})
    {
       $CACAbondRef->{"COEFFICIENT"} = 100.0;
    }
    else
    {        
       $CACAbondRef->{"COEFFICIENT"} = 0.0;
    }
    $bondsRef->[$bondIter]=$CACAbondRef;    
  }

  $mdDataRef->{"BONDS"} = $bondsRef;
}

sub determineAngles
{
  print "Calculating Angles.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  my $anglesRef = [];

  my $mdResiduesRef = $mdDataRef->{"RESIDUES"};
  
  my $angleIter = -1;
  for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef)-2; $residueIter++)
  {
    $angleIter++;
    my $CACACAangleRef = {};
    my $prevCABeadRef = $mdResiduesRef->[$residueIter]->{"BEADS"}->{"CA"};
    my $currentCABeadRef = $mdResiduesRef->[$residueIter+1]->{"BEADS"}->{"CA"};
    my $nextCABeadRef = $mdResiduesRef->[$residueIter+2]->{"BEADS"}->{"CA"};
    $CACACAangleRef->{"FIRST_BEAD_ID"} = $prevCABeadRef->{"ATOM_ID"};
    $CACACAangleRef->{"SECOND_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
    $CACACAangleRef->{"THIRD_BEAD_ID"} = $nextCABeadRef->{"ATOM_ID"};
    $CACACAangleRef->{"ANGLE"} = MDMath::calculateAngle($prevCABeadRef,$currentCABeadRef,$nextCABeadRef);
    if (($prevCABeadRef->{"CHAIN_ID"} eq $currentCABeadRef->{"CHAIN_ID"}) and ($prevCABeadRef->{"CHAIN_ID"} eq $nextCABeadRef->{"CHAIN_ID"}))
    {
       $CACACAangleRef->{"COEFFICIENT"} = 20.0;
    }
    else
    {        
       $CACACAangleRef->{"COEFFICIENT"} = 0.0;
    }

    $anglesRef->[$angleIter]=$CACACAangleRef;
  }
    
  $mdDataRef->{"ANGLES"} = $anglesRef;
}


sub determineDihedrals
{
  print "Calculating Dihedrals.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  my $dihedralsRef = [];
  my $dihedralIter = -1;

  my $mdResiduesRef = $mdDataRef->{"RESIDUES"};
  
  for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef)-3; $residueIter++)
  {
    $dihedralIter++;
    my $CACACACADihedralRef = {};
    my $firstCABeadRef = $mdResiduesRef->[$residueIter] ->{"BEADS"}->{"CA"};
    my $secondCABeadRef = $mdResiduesRef->[$residueIter+1] ->{"BEADS"}->{"CA"};
    my $thirdCABeadRef = $mdResiduesRef->[$residueIter+2] ->{"BEADS"}->{"CA"};
    my $fourthCABeadRef = $mdResiduesRef->[$residueIter+3] ->{"BEADS"}->{"CA"};
    $CACACACADihedralRef->{"FIRST_BEAD_ID"} = $firstCABeadRef->{"ATOM_ID"};
    $CACACACADihedralRef->{"SECOND_BEAD_ID"} = $secondCABeadRef->{"ATOM_ID"};
    $CACACACADihedralRef->{"THIRD_BEAD_ID"} = $thirdCABeadRef->{"ATOM_ID"};
    $CACACACADihedralRef->{"FOURTH_BEAD_ID"} = $fourthCABeadRef->{"ATOM_ID"};
    $CACACACADihedralRef->{"DIHEDRAL"} = MDMath::calculateDihedral($firstCABeadRef,$secondCABeadRef,$thirdCABeadRef,$fourthCABeadRef);
    if (($firstCABeadRef->{"CHAIN_ID"} eq $secondCABeadRef->{"CHAIN_ID"})
    and ($firstCABeadRef->{"CHAIN_ID"} eq $thirdCABeadRef->{"CHAIN_ID"})
    and ($firstCABeadRef->{"CHAIN_ID"} eq $fourthCABeadRef->{"CHAIN_ID"}))
    {
       $CACACACADihedralRef->{"COEFFICIENT1"} = 1.0;
       $CACACACADihedralRef->{"COEFFICIENT2"} = 0.0;
       $CACACACADihedralRef->{"COEFFICIENT3"} = 0.0;
    }
    else
    {        
       $CACACACADihedralRef->{"COEFFICIENT1"} = 0.0;
       $CACACACADihedralRef->{"COEFFICIENT2"} = 0.0;
       $CACACACADihedralRef->{"COEFFICIENT3"} = 0.0;
    }

    $dihedralsRef->[$dihedralIter]=$CACACACADihedralRef;
  }
  $mdDataRef->{"DIHEDRALS"} = $dihedralsRef;
}

sub determineBead
{
      return "CA";
}

sub contactDefined
{
  my ($contactsDataRef,$firstResidueNumber,$firstResidueAtomBead,$secondResidueNumber,$secondResidueAtomBead,$distance) = @_;
  return 0 if (!defined($contactsDataRef->{$firstResidueNumber.".".$firstResidueAtomBead}));
  return (defined($contactsDataRef->{$firstResidueNumber.".".$firstResidueAtomBead}->{$secondResidueNumber.".".$secondResidueAtomBead}));
}

sub addContact
{
  my ($contactsDataRef,$firstResidueNumber,$firstResidueAtomBead,$secondResidueNumber,$secondResidueAtomBead,$distance) = @_;
    $contactsDataRef->{$firstResidueNumber.".".$firstResidueAtomBead}->{$secondResidueNumber.".".$secondResidueAtomBead}++;

}

sub _determineContacts
{
  my ($csuOutputFile, $currentResidueID,$contactsDataRef) = @_;
  
  
  open (CSU_OUTPUT_HANDLE, $csuOutputFile) or die "Error reading the csu output file $csuOutputFile: $!\n";
  my $currentLine;
  my $numberOfContacts = 0;
  
  while ($currentLine = <CSU_OUTPUT_HANDLE>)
  {
      last if ($currentLine =~ m/Full list of atomic contacts formed by/);
  }

  $currentLine = <CSU_OUTPUT_HANDLE>;
  die "Cannot find number of contacts in CSU output file" unless $currentLine =~ m/Total number of contacts is\s+\d+/;
  ($numberOfContacts) = $currentLine =~ m/Total number of contacts is\s+(\d+)/;
  
  # skipping header lines
  $currentLine = <CSU_OUTPUT_HANDLE>;
  $currentLine = <CSU_OUTPUT_HANDLE>;
  $currentLine = <CSU_OUTPUT_HANDLE>;
  $currentLine = <CSU_OUTPUT_HANDLE>;
  $currentLine = <CSU_OUTPUT_HANDLE>;
  
  my $firstResidueNumber = $currentResidueID;
  for (my $contactIter = 0; $contactIter < $numberOfContacts; $contactIter++)
  {
    $currentLine = <CSU_OUTPUT_HANDLE>;
    my ($firstResidueAtomType,$secondResidueNumber,$secondResidueAtomType,$distance) =  $currentLine =~ m/^\s*(\w+)\s+.+\s+\w+\s+(\d+)\s+(\w+)\s+.+\s+(\d+.\d+)/;
    $distance = $distance**2;
    if ((abs($secondResidueNumber-$firstResidueNumber) >= 4))
    {
      my $firstResidueAtomBead = determineBead($firstResidueAtomType);
      my $secondResidueAtomBead = determineBead($secondResidueAtomType);
      if ($secondResidueNumber > $firstResidueNumber)
      {
        addContact($contactsDataRef,$firstResidueNumber,$firstResidueAtomBead,$secondResidueNumber,$secondResidueAtomBead,$distance);
      }
      else
      {
        if (!contactDefined($contactsDataRef,$secondResidueNumber,$secondResidueAtomBead,$firstResidueNumber,$firstResidueAtomBead))
        {
          addContact($contactsDataRef,$secondResidueNumber,$secondResidueAtomBead,$firstResidueNumber,$firstResidueAtomBead,$distance);
        }
      }
    }
  }
  close CSU_OUTPUT_HANDLE;
}

sub translateChainContactsToBeads
{
  my ($currentChainID,$chainContactsDataRef,$contactsDataRef,$mdDataRef,$executionPreferencesRef) = @_;
  print "Translating CSU contacts to bead contacts for chain $currentChainID\n";
  my $mdChainRef = $mdDataRef->{"CHAINS"}->{$currentChainID};
  foreach my $chainContactFirstBead (keys(%$chainContactsDataRef))
  {
    my ($firstResidueID,$firstBeadType) = split (/\./,$chainContactFirstBead);
    my $firstResidueRef = $mdChainRef->{"RESIDUES"}->{$firstResidueID};
    my $firstBeadRef;
    if ($firstResidueRef->{"RESIDUE_TYPE"} eq "GLY")
    {
      $firstBeadRef = $mdChainRef->{"RESIDUES"}->{$firstResidueID}->{"BEADS"}->{"CA"};
    }
    else
    {
      $firstBeadRef = $mdChainRef->{"RESIDUES"}->{$firstResidueID}->{"BEADS"}->{$firstBeadType};
    }
    foreach my $chainContactSecondBead (keys(%{$chainContactsDataRef->{$chainContactFirstBead}}))
    {
      my ($secondResidueID,$secondBeadType) = split (/\./,$chainContactSecondBead);
      my $secondResidueRef = $mdChainRef->{"RESIDUES"}->{$secondResidueID};
      my $secondBeadRef;
      if ($secondResidueRef->{"RESIDUE_TYPE"} eq "GLY")
      {
        $secondBeadRef = $mdChainRef->{"RESIDUES"}->{$secondResidueID}->{"BEADS"}->{"CA"};
      }
      else
      {
        $secondBeadRef = $mdChainRef->{"RESIDUES"}->{$secondResidueID}->{"BEADS"}->{$secondBeadType};
      }

      my $CBdistance = MDMath::calculateDistance($firstBeadRef,$secondBeadRef);
      $CBdistance = $CBdistance**2;
      $contactsDataRef->{$firstBeadRef->{"ATOM_ID"}.".".$secondBeadRef->{"ATOM_ID"}}->{"DISTANCE"} = $CBdistance;
      $contactsDataRef->{$firstBeadRef->{"ATOM_ID"}.".".$secondBeadRef->{"ATOM_ID"}}->{"COEFFICIENT"} = 1.0;
    }
  }
}

sub translateContactsToBeads
{
  print "Translating CSU contacts to bead contacts.\n";
  my ($CSUContactsDataRef,$mdDataRef,$executionPreferencesRef,$contactsDataRef) = @_;
  foreach my $contactFirstBead (keys(%$CSUContactsDataRef))
  {
    my ($firstResidueIter,$firstBeadType) = split (/\./,$contactFirstBead);
    $firstResidueIter--;
    my $firstResidueRef = $mdDataRef->{"RESIDUES"}->[$firstResidueIter];
    my $firstBeadRef;
    if ($firstResidueRef->{"RESIDUE_TYPE"} eq "GLY")
    {
      $firstBeadRef = $mdDataRef->{"RESIDUES"}->[$firstResidueIter]->{"BEADS"}->{"CA"};
    }
    else
    {
      $firstBeadRef = $mdDataRef->{"RESIDUES"}->[$firstResidueIter]->{"BEADS"}->{$firstBeadType};
    }
    foreach my $contactSecondBead (keys(%{$CSUContactsDataRef->{$contactFirstBead}}))
    {
      my ($secondResidueIter,$secondBeadType) = split (/\./,$contactSecondBead);
      $secondResidueIter--;
      my $secondResidueRef = $mdDataRef->{"RESIDUES"}->[$secondResidueIter];
      my $secondBeadRef;
      if ($secondResidueRef->{"RESIDUE_TYPE"} eq "GLY")
      {
        $secondBeadRef = $mdDataRef->{"RESIDUES"}->[$secondResidueIter]->{"BEADS"}->{"CA"};
      }
      else
      {
        $secondBeadRef = $mdDataRef->{"RESIDUES"}->[$secondResidueIter]->{"BEADS"}->{$secondBeadType};
      }

      my $CBdistance = MDMath::calculateDistance($firstBeadRef,$secondBeadRef);
      $CBdistance = $CBdistance**2;
      $contactsDataRef->{$firstBeadRef->{"ATOM_ID"}.".".$secondBeadRef->{"ATOM_ID"}}->{"DISTANCE"} = $CBdistance;
      $contactsDataRef->{$firstBeadRef->{"ATOM_ID"}.".".$secondBeadRef->{"ATOM_ID"}}->{"COEFFICIENT"} = 1.0;
    }
  }
}


sub determineContacts
{

  print "Calculating CSU contacts.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;

  my @dnaBeads =[];
    my @ignoreRange = [-1,-10];#a range that dose not exist
    if( exists($mdDataRef->{"BEADS_DNA"}) && scalar (@{$mdDataRef->{"BEADS_DNA"}})>0 ){
        my $firstDNAbead =$mdDataRef->{"BEADS_DNA"}->[0];
        my $lastDNAbead=$mdDataRef->{"BEADS_DNA"}->[scalar(@{$mdDataRef->{"BEADS_DNA"}})-1];
        my $firstDNAres = getPdbDataAtom($firstDNAbead,$pdbDataRef);
        my $lastDNAres = getPdbDataAtom($lastDNAbead,$pdbDataRef);
        print "firstDnaAtom ".$firstDNAres->{"ATOMS"}->[0]->{"ID"}." lastDnaAtom ".
            $lastDNAres->{"ATOMS"}->[scalar(@{$lastDNAres->{"ATOMS"}})-1]->{"ID"}."\n";
        @ignoreRange =[$firstDNAres->{"ATOMS"}->[0]->{"ID"},
                       $lastDNAres->{"ATOMS"}->[scalar(@{$lastDNAres->{"ATOMS"}})-1]->{"ID"}];
    }
  my $scriptDirectory = $executionPreferencesRef->{"SCRIPT_DIRECTORY"};
# ohad  require $utilDirectory."PDBHandling.pm";
  require "/home/ohad/scripts/fromAriel/PDBHandling.pm"; 
  my $contactsDataRef = {};
  
  if ($executionPreferencesRef->{"INTER_CHAIN_CONTACTS"} =~ m/NO/)
  {
    foreach my $pdbChainRef (@{$pdbDataRef->{"CHAINS"}})
    {
      my $chainContactsDataRef = {};
      my $currentChainID = $pdbChainRef->{"CHAIN_ID"};
      PDBHandling::writeChainAsPDBFile($pdbChainRef,"csu.input.".$currentChainID,@ignoreRange);
      
      # this foreach loop was replaced by the for loop below (NOT TESTED !!!!!)
      #
      #foreach my $pdbResidueRef (@{$pdbChainRef->{"RESIDUES"}})
      #{
      #  
      #  my $currentResidueID = $pdbResidueRef->{"RESIDUE_ID"};
      #  `$scriptDirectory/resc csu.input.$currentChainID $currentResidueID >> csu.output.$currentChainID.$currentResidueID`;
      #  _determineContacts("csu.output.".$currentChainID.".".$currentResidueID,$currentResidueID,$chainContactsDataRef);
      #  unlink("csu.output.".$currentChainID.".".$currentResidueID);
      #}
    my $Dyresnum = 
    $executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/? $executionPreferencesRef->{"LAST_DYNAMIC_ATOM"}  : 
    scalar(@{$pdbDataRef->{"RESIDUES"}});
    for (my $residueIter = 1; $residueIter <= $Dyresnum; $residueIter++){
        my $pdbResidueRef = $pdbChainRef->{"RESIDUES"}->[$residueIter];
        my $currentResidueID = $pdbResidueRef->{"RESIDUE_ID"};
        `$scriptDirectory/resc csu.input.$currentChainID $currentResidueID >> csu.output.$currentChainID.$currentResidueID`;
        _determineContacts("csu.output.".$currentChainID.".".$currentResidueID,$currentResidueID,$chainContactsDataRef);
        unlink("csu.output.".$currentChainID.".".$currentResidueID);
      }
      unlink("csu.input.".$currentChainID);
      translateChainContactsToBeads($currentChainID,$chainContactsDataRef,$contactsDataRef,$mdDataRef,$executionPreferencesRef);  
    }
  }
  else
  {
    my $CSUContactsDataRef = {};
    PDBHandling::writeAllChainsAsSingleChainPDBFile($pdbDataRef->{"RESIDUES"},"csu.input.A",@ignoreRange);
    #this will have to change whene we insert renges of static atoms
    if ($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/) {
        my $DyRangeList = $executionPreferencesRef->{"DYNAMIC_ATOM_RANGE"};
        my @DyRange = split /,/, $DyRangeList;
        my $DyRangeLong = @DyRange;
        my $i = 0;
#        print "*** The array contains: $DyRangeLong values ***\n";
#        print "*** The first value is: $DyRange[0] ***\n";
#        print "*** The second value is: $DyRange[1] ***\n";
#        print "*** The 3th value is: $DyRange[2] ***\n";
#        print "*** The 4th value is: $DyRange[3] ***\n";
        while ($i < $DyRangeLong) {
            for (my $residueIter = $DyRange[$i]; $residueIter <= $DyRange[$i+1]; $residueIter++)
            {
                `$scriptDirectory/resc csu.input.A $residueIter >> csu.output.A.$residueIter`;
                _determineContacts("csu.output.A.".$residueIter,$residueIter,$CSUContactsDataRef);
                unlink("csu.output.A.".$residueIter);
            }
            $i = $i + 2;
        }
    }
    else {
        my $Dyresnum = scalar(@{$pdbDataRef->{"RESIDUES"}});
        for (my $residueIter = 1; $residueIter <= $Dyresnum; $residueIter++)
        {
          `$scriptDirectory/resc csu.input.A $residueIter >> csu.output.A.$residueIter`;
          _determineContacts("csu.output.A.".$residueIter,$residueIter,$CSUContactsDataRef);
          unlink("csu.output.A.".$residueIter);
        }
    }
#    my $Dyresnum = 
#    $executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/ ? $executionPreferencesRef->{"DYNAMIC_ATOM_RANGE"}  : 
#    scalar(@{$pdbDataRef->{"RESIDUES"}});
#    for (my $residueIter = 1; $residueIter <= $Dyresnum; $residueIter++)
#    {
#      `$scriptDirectory/resc csu.input.A $residueIter >> csu.output.A.$residueIter`;
#      _determineContacts("csu.output.A.".$residueIter,$residueIter,$CSUContactsDataRef);
#      unlink("csu.output.A.".$residueIter);
#    }
    unlink("csu.input.A");
    translateContactsToBeads($CSUContactsDataRef,$mdDataRef,$executionPreferencesRef,$contactsDataRef);  
  }
  $mdDataRef->{"CONTACTS"} = $contactsDataRef;
}

#############################################
#input a bead from $mdDataRef and $pdbDataRef
#output the residue from $pdbDataRef
####################################
sub getPdbDataAtom {
    my ($mdDataBead,$pdbDataRef) = @_;
    my $chain = $mdDataBead->{"CHAIN_ID"};
    my $residue = $mdDataBead->{"RESIDUE_ID"};
    my $pdbRes;
    foreach my $pdbRes1 (@{$pdbDataRef->{"RESIDUES"}}){
        $pdbRes = $pdbRes1;
        last if ($pdbRes1->{"CHAIN_ID"} eq $chain &&
              $pdbRes1->{"RESIDUE_ID"} eq $residue );
    }
    die "requsted residue $chain $residue not found in getPdbDataAtom\n"
        unless ($pdbRes->{"CHAIN_ID"} eq $chain &&
                $pdbRes->{"RESIDUE_ID"} eq $residue);
    return $pdbRes; 
}


sub calculateRepulsionDistance
{
  my ($firstBeadType,$secondBeadType) = @_;
  $firstBeadType =~ s/\s+//;
  $secondBeadType=~ s/\s+//;
  return ($beadRepulsionRadius->{$firstBeadType} + $beadRepulsionRadius->{$secondBeadType});
  
}
sub beadsShouldRepulse
{
    my ($firstChainID,$firstResidueID,$firstBeadType,$secondChainID,$secondResidueID,$secondBeadType) = @_;
    return (((($secondResidueID-$firstResidueID) >= 4)) or ($firstChainID ne $secondChainID));
}

sub determineRepulsions
{
  print "Calculating Repulsions.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  $beadRepulsionRadius->{"P"} = $executionPreferencesRef->{"DNA_REPULSION_DISTANCE"};
  $beadRepulsionRadius->{"S"} = $executionPreferencesRef->{"DNA_REPULSION_DISTANCE"};
  $beadRepulsionRadius->{"B"} = $executionPreferencesRef->{"DNA_REPULSION_DISTANCE"};
  
  my $contactsDataRef = $mdDataRef->{"CONTACTS"};
  my $repulsionsDataRef = {};
  my @allBeads;
  if (exists($mdDataRef->{"BEADS"})){
    push(@allBeads,@{$mdDataRef->{"BEADS"}});
  }
  if (exists($mdDataRef->{"BEADS_DNA"})){
    push(@allBeads,@{$mdDataRef->{"BEADS_DNA"}});
  }
  #For dynamic DNA this must change
#  my $lastDynamicBead=
#    $executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/? $executionPreferencesRef->{"LAST_DYNAMIC_ATOM"}  : 
#    scalar(@allBeads);
  my @DynamicBeadRange;
  if ($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/)
  {
      my $DynamicBeadList = $executionPreferencesRef->{"DYNAMIC_ATOM_RANGE"};
      @DynamicBeadRange = split /,/, $DynamicBeadList;
  } else {
      @DynamicBeadRange = (1,scalar(@allBeads));
  }
  my $i = 0;
  my $DyBeadRangeLong = @DynamicBeadRange;
#  print "***222*** The array contains: $DyBeadRangeLong values ***\n";
#  print "***222*** The first value is: $DynamicBeadRange[0] ***\n";
#  print "***222*** The second value is: $DynamicBeadRange[1] ***\n";
  while ($i < $DyBeadRangeLong)
  {
#  for (my $firstBeadIter = 0; $firstBeadIter < $lastDynamicBead; $firstBeadIter++)
    for (my $firstBeadIter = $DynamicBeadRange[$i]-1; $firstBeadIter < $DynamicBeadRange[$i+1]; $firstBeadIter++)
      {
        my $firstBeadID = $allBeads[$firstBeadIter]->{"ATOM_ID"};
        my $firstChainID = $allBeads[$firstBeadIter]->{"CHAIN_ID"};
        my $firstResidueID = $allBeads[$firstBeadIter]->{"RESIDUE_ID"};
        my $firstBeadType = $allBeads[$firstBeadIter]->{"ATOM_TYPE"};
        for (my $secondBeadIter = $firstBeadIter+1; $secondBeadIter < scalar(@allBeads); $secondBeadIter++)
        {
          my $secondBeadID = $allBeads[$secondBeadIter]->{"ATOM_ID"};
          my $secondChainID = $allBeads[$secondBeadIter]->{"CHAIN_ID"};
          my $secondResidueID = $allBeads[$secondBeadIter]->{"RESIDUE_ID"};
          my $secondBeadType = $allBeads[$secondBeadIter]->{"ATOM_TYPE"};
          
          if (beadsShouldRepulse($firstChainID,$firstResidueID,$firstBeadType,$secondChainID,$secondResidueID,$secondBeadType)
          and (!defined($mdDataRef->{"CONTACTS"}->{$firstBeadID.".".$secondBeadID}))
          and ($firstBeadID <=$DynamicBeadRange[$i+1] or  $secondBeadID <=$DynamicBeadRange[$i+1] ))
          {
            $repulsionsDataRef->{$firstBeadID.".".$secondBeadID}->{"DISTANCE"} = calculateRepulsionDistance($firstBeadType,$secondBeadType)**2;
            $repulsionsDataRef->{$firstBeadID.".".$secondBeadID}->{"COEFFICIENT"} = 1.0;
          }
        }
      }
    $i = $i + 2;
  }
  
  $mdDataRef->{"REPULSIONS"} = $repulsionsDataRef;
}


sub determineCharge
{
    my $currentResidueType = $_[0];
    
    return (exists($residueCharges->{$currentResidueType}) ? $residueCharges->{$currentResidueType} : 0);
}

sub determineElectrostatics
{
    print "Calculating electrostatic residues.\n";
    my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
    my $electrostaticsRef = {};
    if (exists($mdDataRef->{"BEADS"})){
    for (my $beadIter = 0; $beadIter < scalar(@{$mdDataRef->{"BEADS"}}); $beadIter++)
    {
        my $currentCABeadID = $mdDataRef->{"BEADS"}->[$beadIter] ->{"ATOM_ID"};
        my $currentCharge = determineCharge($mdDataRef->{"BEADS"}->[$beadIter]->{"RESIDUE_TYPE"});
        if ($currentCharge != 0){
            $electrostaticsRef->{$currentCABeadID} = $currentCharge;
        }
    }
    }
    if (exists($mdDataRef->{"BEADS_DNA"})){
    for (my $beadIter = 0; $beadIter < scalar(@{$mdDataRef->{"BEADS_DNA"}}); $beadIter++)
    {
        my $currentDNABeadID = $mdDataRef->{"BEADS_DNA"}->[$beadIter] ->{"ATOM_ID"};
        my $currentCharge = $mdDataRef->{"BEADS_DNA"}->[$beadIter]->{"ATOM_TYPE"} =~ / *P */? -1 :0;
        if ($currentCharge != 0){
            $electrostaticsRef->{$currentDNABeadID} = $currentCharge;
        }
    }
    }
    $mdDataRef->{"ELECTROSTATICS"} = $electrostaticsRef;
}

#sub determineElectrostatics
#{
#  print "Calculating electrostatic residues.\n";
#  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
#  my $electrostaticsRef = {};
#
#  for (my $residueIter = 0; $residueIter < scalar(@{$mdDataRef->{"RESIDUES"}}); $residueIter++)
#  {
#    
#    my $currentResidueType = $mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}->{"CA"}->{"RESIDUE_TYPE"};
#    my $currentCharge = determineCharge($currentResidueType);
#    if ($currentCharge != 0)
#    {
#        my $currentCABeadID =  $mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}->{"CA"}->{"ATOM_ID"};
#        $electrostaticsRef->{$currentCABeadID} = $currentCharge;
#    }
#  }
#
#  $mdDataRef->{"ELECTROSTATICS"} = $electrostaticsRef;
#
#}

sub printMDInputFile
{
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;

  my $outputFile = $executionPreferencesRef->{"MD_INPUT_FILE"};

  print "Printing to output file $outputFile.\n";
  
  open(MD_INPUT_FILE_HANDLE,">",$outputFile);
  
  my $mdBondsRef = $mdDataRef->{"BONDS"};
  printf MD_INPUT_FILE_HANDLE ("%12d Hookean bonded pairs.\n",scalar(@$mdBondsRef));
  for (my $bondIter = 0; $bondIter < scalar(@$mdBondsRef); $bondIter++)
  {
    printf MD_INPUT_FILE_HANDLE ("%5d%5d%5d%8.3f%8.3f\n",
                                 $bondIter+1,
                                 $mdBondsRef->[$bondIter]->{"FIRST_BEAD_ID"},
                                 $mdBondsRef->[$bondIter]->{"SECOND_BEAD_ID"},
                                 $mdBondsRef->[$bondIter]->{"DISTANCE"},
                                 $mdBondsRef->[$bondIter]->{"COEFFICIENT"});
  }

  my $mdAnglesRef = $mdDataRef->{"ANGLES"};
  printf MD_INPUT_FILE_HANDLE ("%12d bond angle trios.\n",scalar(@$mdAnglesRef));
  for (my $angleIter = 0; $angleIter < scalar(@$mdAnglesRef); $angleIter++)
  {
    printf MD_INPUT_FILE_HANDLE ("%5d%5d%5d%5d%8.3f%8.3f\n",
                                 $angleIter+1,
                                 $mdAnglesRef->[$angleIter]->{"FIRST_BEAD_ID"},
                                 $mdAnglesRef->[$angleIter]->{"SECOND_BEAD_ID"},
                                 $mdAnglesRef->[$angleIter]->{"THIRD_BEAD_ID"},
                                 $mdAnglesRef->[$angleIter]->{"ANGLE"},
                                 $mdAnglesRef->[$angleIter]->{"COEFFICIENT"});
  }

  my $mdDihedralsRef = $mdDataRef->{"DIHEDRALS"};
  printf MD_INPUT_FILE_HANDLE ("%12d dihedral quartets.\n",scalar(@$mdDihedralsRef));
  for (my $dihedralIter = 0; $dihedralIter < scalar(@$mdDihedralsRef); $dihedralIter++)
  {
    printf MD_INPUT_FILE_HANDLE ("%5d%5d%5d%5d%5d%8.3f%8.3f%8.3f%8.3f\n",
                                 $dihedralIter+1,
                                 $mdDihedralsRef->[$dihedralIter]->{"FIRST_BEAD_ID"},
                                 $mdDihedralsRef->[$dihedralIter]->{"SECOND_BEAD_ID"},
                                 $mdDihedralsRef->[$dihedralIter]->{"THIRD_BEAD_ID"},
                                 $mdDihedralsRef->[$dihedralIter]->{"FOURTH_BEAD_ID"},
                                 $mdDihedralsRef->[$dihedralIter]->{"DIHEDRAL"},
                                 $mdDihedralsRef->[$dihedralIter]->{"COEFFICIENT1"},
                                 $mdDihedralsRef->[$dihedralIter]->{"COEFFICIENT2"},
                                 $mdDihedralsRef->[$dihedralIter]->{"COEFFICIENT3"});
  }
  my $mdContactsRef = $mdDataRef->{"CONTACTS"};
  printf MD_INPUT_FILE_HANDLE ("%12d contacts.\n",scalar(keys(%$mdContactsRef)));
  my $contactIter = 0;
  foreach my $contactKey (sort contactCompare keys(%$mdContactsRef))
  {
    $contactIter++;
    my ($firstBeadID,$secondBeadID) = split (/\./,$contactKey);
    printf MD_INPUT_FILE_HANDLE ("%5d%5d%5d%10.3f%9.6f\n",
                                 $contactIter,
                                 $firstBeadID,
                                 $secondBeadID,
                                 $mdContactsRef->{$contactKey}->{"DISTANCE"},
                                 $mdContactsRef->{$contactKey}->{"COEFFICIENT"});
  }

  my $mdRepulsionsRef = $mdDataRef->{"REPULSIONS"};
  printf MD_INPUT_FILE_HANDLE ("%12d repulsive pairs.\n",scalar(keys(%$mdRepulsionsRef)));
  my $repulsionIter = 0;
  foreach my $repulsionKey (sort contactCompare keys(%$mdRepulsionsRef))
  {
    $repulsionIter++;
    my ($firstBeadID,$secondBeadID) = split (/\./,$repulsionKey);
    printf MD_INPUT_FILE_HANDLE ("%8d%5d%5d%10.3f%10.3f\n",
                                 $repulsionIter,
                                 $firstBeadID,
                                 $secondBeadID,
                                 $mdRepulsionsRef->{$repulsionKey}->{"DISTANCE"},
                                 $mdRepulsionsRef->{$repulsionKey}->{"COEFFICIENT"},
                                 );
  }

  if ($executionPreferencesRef->{"ELECTROSTATICS"} =~ m/YES/i)
  {
    my $mdElectrostaticsRef = $mdDataRef->{"ELECTROSTATICS"};
    printf MD_INPUT_FILE_HANDLE ("%12d electrostatic residues.\n",scalar(keys(%$mdElectrostaticsRef)));
    my $erIter = 0;
    foreach my $beadID (sort { $a <=> $b } keys(%$mdElectrostaticsRef))
    {
      $erIter++;
      printf MD_INPUT_FILE_HANDLE ("%8d%5d%8.3f\n",
                                   $erIter,
                                   $beadID,
                                   $mdElectrostaticsRef->{$beadID});
    }
  }
  
    #printing atom possitions 
  my $atomNumber = (exists (  $mdDataRef->{"BEADS_DNA"}))
    ? scalar(@{$mdDataRef->{"BEADS"}})+scalar(@{$mdDataRef->{"BEADS_DNA"}})
    : scalar(@{$mdDataRef->{"BEADS"}});
  printf MD_INPUT_FILE_HANDLE ("%12d atom positions.\n",$atomNumber);
  printf MD_INPUT_FILE_HANDLE ("%12d chains.\n",scalar(keys(%{$mdDataRef->{"CHAINS"}})));

  my $chainIter = 0;

  my @chains = sort(
    {
     chainsBeads($mdDataRef->{"CHAINS"}->{$a})->[0]->{"ATOM_ID"} <=>
     chainsBeads($mdDataRef->{"CHAINS"}->{$b})->[0]->{"ATOM_ID"}
     }
    keys (%{$mdDataRef->{"CHAINS"}})  );
  foreach my $mdChainID (@chains)
   
  {
    $chainIter++;
    my $chainLength =0;
    if (exists($mdDataRef->{"CHAINS"}->{$mdChainID}->{"BEADS"})){
        $chainLength+=scalar(@{$mdDataRef->{"CHAINS"}->{$mdChainID}->{"BEADS"}});
    }
    if (exists($mdDataRef->{"CHAINS"}->{$mdChainID}->{"BEADS_DNA"})){
        $chainLength+=scalar(@{$mdDataRef->{"CHAINS"}->{$mdChainID}->{"BEADS_DNA"}});
    }

    printf MD_INPUT_FILE_HANDLE ("%12d is the size of chain%12d.\n",
                                 $chainLength,
                                 $chainIter);
  }
  foreach my $bead (@{$mdDataRef->{"BEADS"}}){
    printf MD_INPUT_FILE_HANDLE ("%5d%4d%3s%4s%8.3f%8.3f%8.3f%8.3f\n",
                                     $bead->{"ATOM_ID"},
                                     $bead->{"ATOM_ID"},
                                     $bead->{"ATOM_TYPE"},
                                     $bead->{"RESIDUE_TYPE"},
                                     $bead->{"X"},
                                     $bead->{"Y"},
                                     $bead->{"Z"},
                                     $bead->{"MASS"});
  }
  #if you need to make sure this prints P -> S -> B
  # you must go over all residues and print the BEADS_DNA in each residue
  #according to the desired order
  foreach my $bead (@{$mdDataRef->{"BEADS_DNA"}}){
    printf MD_INPUT_FILE_HANDLE ("%5d%4d%3s%4s%8.3f%8.3f%8.3f%8.3f\n",
                                     $bead->{"ATOM_ID"},
                                     $bead->{"ATOM_ID"},
                                     $bead->{"ATOM_TYPE"},
                                     $bead->{"RESIDUE_TYPE"},
                                     $bead->{"X"},
                                     $bead->{"Y"},
                                     $bead->{"Z"},
                                     $bead->{"MASS"});
  }
  
        

  #foreach my $mdChainID (@chains)
  #{
  #  my $mdResiduesRef = $mdDataRef->{"CHAINS"}->{$mdChainID}->{"RESIDUES"};
  #  my $beadIter = 0;
  #  foreach my $mdResidueID (sort {$a <=> $b} keys(%$mdResiduesRef))
  #  {
  #      my @beadNames =(sort
  #          {$mdResiduesRef->{$mdResidueID}->{"BEADS"}->{$a}->{"ATOM_ID"} <=>
  #          $mdResiduesRef->{$mdResidueID}->{"BEADS"}->{$b}->{"ATOM_ID"} }
  #          (keys (%{$mdResiduesRef->{$mdResidueID}->{"BEADS"}}))); 
  #      foreach my $beadName (@beadNames){
  #      $beadIter++;
  #      my $mdCABeadRef = $mdResiduesRef->{$mdResidueID}->{"BEADS"}->{$beadName};
  #      printf MD_INPUT_FILE_HANDLE ("%5d%4d%3s%4s%8.3f%8.3f%8.3f%8.3f\n",
  #                                   $beadIter,
  #                                   $mdCABeadRef->{"ATOM_ID"},
  #                                   $beadName,
  #                                   $mdCABeadRef->{"RESIDUE_TYPE"},
  #                                   $mdCABeadRef->{"X"},
  #                                   $mdCABeadRef->{"Y"},
  #                                   $mdCABeadRef->{"Z"},
  #                                   $mdCABeadRef->{"MASS"});
  #      }
  #      
  #      next if ($mdResiduesRef->{$mdResidueID}->{"RESIDUE_TYPE"} eq "GLY");
  #  }
  #}  
  ##foreach my $mdChainID (sort keys (%{$mdDataRef->{"CHAINS"}}))
  ##{
  ##  my $mdResiduesRef = $mdDataRef->{"CHAINS"}->{$mdChainID}->{"RESIDUES"};
  ##  my $beadIter = 0;
  ##  foreach my $mdResidueID (sort {$a <=> $b} keys(%$mdResiduesRef))
  ##  {
  ##      $beadIter++;
  ##      my $mdCABeadRef = $mdResiduesRef->{$mdResidueID}->{"BEADS"}->{"CA"};
  ##      printf MD_INPUT_FILE_HANDLE ("%5d%4d CA %3s%8.3f%8.3f%8.3f%8.3f\n",
  ##                                   $beadIter,
  ##                                   $mdCABeadRef->{"ATOM_ID"},
  ##                                   $mdCABeadRef->{"RESIDUE_TYPE"},
  ##                                   $mdCABeadRef->{"X"},
  ##                                   $mdCABeadRef->{"Y"},
  ##                                   $mdCABeadRef->{"Z"},
  ##                                   $mdCABeadRef->{"MASS"});
  ##      
  ##      next if ($mdResiduesRef->{$mdResidueID}->{"RESIDUE_TYPE"} eq "GLY");
  ##  }
  ##}

  close(MD_INPUT_FILE_HANDLE);
}

#returns all bead in the chain first protein the DNA
sub chainsBeads{
    my ($chain) =@_;
    my @output ;
    if (exists ($chain->{"BEADS"})){
        push(@output,@{$chain->{"BEADS"}});
    }
    if (exists ($chain->{"BEADS_DNA"})){
        push(@output,@{$chain->{"BEADS_DNA"}});
    }
    return \@output;
}

sub calcCenter{
    my @atoms = @_ ;
    my $y=0;my $z=0;my $x = 0;
    foreach my $atom (@atoms){
        $x+=$atom->{"X"}; 
        $y+=$atom->{"Y"};
        $z+=$atom->{"Z"};
    }
    $x/=scalar(@atoms);
    $y/=scalar(@atoms);
    $z/=scalar(@atoms);
    return(sprintf("%.3f",$x),sprintf("%.3f",$y),sprintf("%.3f",$z)); 
}

sub contactCompare
{ 
  my @a  = split (/\./,$a);
  my @b  = split (/\./,$b);
  return 1  if ($a[0] > $b[0]);
  return -1  if ($a[0] < $b[0]);
  return $a[1] <=> $b[1];
}


#####################
#$mdDataRef
# {ANGLES}
# {BEADS} ->array[BEAD]
# {BEADS_DNA}->array[BEAD](containing dna beads)
# {CHAINS}->hash{CHAIN}
#  ...
# {RESIDUES} -> array[RESIDUE]
#
#
##CHAIN
#   {BEADS} ->array[BEAD]
#   {BEADS_DNA}->array[BEAD](containing dna beads)
#   {RESIDUES} -> array[RESIDUE]
#   size    scalar
#
#RESIDUE
# {BEADS}
# {BEADS_DNA}
# {RESIDUE_ID}
# {RESIDUE_TYPE}
#
#BEAD
#  {ATOM_ID}  (= bead id)
#  {ATOM_TYPE}
#  {CHAIN_ID}
#  {MASS}
#  {RESIDUE_ID}
#  {RESIDUE_TYPE}
#  {X}
#  {Y}
#  {Z}
#
#####################