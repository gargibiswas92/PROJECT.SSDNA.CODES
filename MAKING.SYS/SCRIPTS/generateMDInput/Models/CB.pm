package CB;
use strict;
use warnings;
use Math::Trig;

use base 'Exporter';
our @EXPORT_OK = qw(createMDInputFile);
#see explanation regarding data structures at end of file

#constants
#my $utilDirectory = "/home_c/arielaz/scripts/util/";
# my $utilDirectory = "/home_d/arumay/scripts/util/";
my $utilDirectory = $scriptDirectory."../../util/";
require $utilDirectory."MDMath.pm";
#################################################################
#correction 1.0 ohad givaty adding static DNA
#1.1 AGNES adding protein-DNA native contacts and contact ranges 
################################################################
my @DynamicRange;
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

my $beadRepulsionRadius = { "CA" => 1.9,
                            "CB" => 1.5 };

my $minimumRepulsionDistance_short = { "CA_CA" => 4,
                                       "CA_CB" => 2,
                                       "CB_CA" => 2,
                                       "CB_CB" => 1} ;

my $minimumRepulsionDistance_long = { "CA_CA" => 4,
                                      "CA_CB" => 2,
                                      "CB_CA" => 2,
                                      "CB_CB" => 2} ;

my $minimumContactDistance = { "CA_CA" => 4,
                               "CA_CB" => 3,
                               "CB_CA" => 3,
                               "CB_CB" => 2} ;

#change this to arrays from wich you must take center mass
my $farthestHeavyAtomByResidueType = {"MET" => "CG",
                                      "ASP" => "OD1",
                                      "SER" => "OG",
                                      "ALA" => "CB",
                                      "ILE" => "CD1",
                                      "THR" => "OG1",
                                      "LEU" => "CD1",
                                      "TRP" => "CZ2",
                                      "GLN" => "NE2",
                                      "PHE" => "CZ",
                                      "LYS" => "NZ",
                                      "PRO" => "CG",
                                      "ASN" => "OD1",
                                      "HIS" => "NE2",
                                      "CYS" => "SG",
                                      "GLU" => "OE1",
                                      "VAL" => "CG1",
                                      "ARG" => "NH1",
                                      "TYR" => "OH"};
sub createMDInputFile
{
    my ($executionPreferencesRef,$pdbDataRef,@ChainOrder) = @_;
    my $mdDataRef = {};
    determineBeads($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineDNABeads($executionPreferencesRef,$pdbDataRef,$mdDataRef);## 1.0
    setDynamicRange($executionPreferencesRef,$mdDataRef);#1.0 @DynamicRange = 
    determineBonds($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineAngles($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineDihedrals($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    determineChirals($executionPreferencesRef,$pdbDataRef,$mdDataRef) if ($executionPreferencesRef->{"CHIRALITY"} =~ m/YES/i);
    determineContacts($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    if ($executionPreferencesRef->{"HAS_PROTEIN_DNA_CONTACTS"} =~ m/YES/i)
    {
     determineDNAContacts($executionPreferencesRef,$pdbDataRef,$mdDataRef);
    }
    determineRepulsions($executionPreferencesRef,$mdDataRef);
    if ($executionPreferencesRef->{"CB_POSITION"}  =~ m/FHA/i)
    {
      determineEllipsoidRepulsions($executionPreferencesRef,$mdDataRef)  if ($executionPreferencesRef->{"REPULSION_ELLIPSOIDS"}  =~ m/YES/i);
    }
    determineElectrostatics($executionPreferencesRef,$pdbDataRef,$mdDataRef) if ($executionPreferencesRef->{"ELECTROSTATICS"} =~ m/YES/i);
    printMDInputFile($executionPreferencesRef,$pdbDataRef,$mdDataRef,@ChainOrder);
}



sub determineDNABeads
{
    #my @validResTypes = ["C","A","G","T"];
    my ($executionPreferencesRef,$pdbDataRef,$mdDataRef)=@_;
    print "Generating DNA Beads.\n";
    my $beadsDataRef = [];
    $mdDataRef->{"BASES"}=[];
    my $dnaResiduesDataRef = $mdDataRef->{"BASES"};
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
    
        #take only DNA residues
        if( length ($$atom{"RESIDUE_TYPE"}) < 3){
            $curentRes=$$atom{"RESIDUE_ID"};
        }
        else{
            next;
        }

        # handle last residue
        # (when we get to the next residue $lastRes ne $curentRes )
        if(    $curentRes ne $lastRes  ){
            #not first iteration
            if ($lastRes ne -1000){
                
                #handle last residue
                my $lastAtom =  $base[0];
                my($phospateBead,$baseBead,$sugarBead,$currentBeadNum1)
                =makeDNARes(\@phospate  ,\@base,\@sugar,$currentBeadNum,$lastAtom,$currentChain,$executionPreferencesRef);
                $currentBeadNum =$currentBeadNum1;
                
                $currentResidueRef->{"RESIDUE_ID"} =  $lastAtom->{"RESIDUE_ID"};
                $currentResidueRef->{"RESIDUE_TYPE"} = $lastAtom->{"RESIDUE_TYPE"};
                if (keys(%$phospateBead) != 0)
                {
                    $currentResidueRef->{"BEADS_DNA"}->{"P"} = $phospateBead;
                    push(@$beadsDataRef, $phospateBead);
                    push(@{$currentChainRef->{"BEADS_DNA"}},$phospateBead);
                }
                if (keys(%$sugarBead) != 0)
                {
                    $currentResidueRef->{"BEADS_DNA"}->{"S"} = $sugarBead;
                    push(@{$currentChainRef->{"BEADS_DNA"}},$sugarBead);
                    push(@$beadsDataRef, $sugarBead);
                }
                if (keys(%$baseBead) != 0)
                {
                    $currentResidueRef->{"BEADS_DNA"}->{"B"} = $baseBead;
                    push(@{$currentChainRef->{"BEADS_DNA"}},$baseBead);
                    push(@$beadsDataRef, ,$baseBead);
                }
                
                if (keys(%$currentResidueRef) != 0)
                {
                    push(@$dnaResiduesDataRef,$currentResidueRef);
                    $currentChainRef->{"BASES"}->{$lastAtom->{"RESIDUE_ID"}} = $currentResidueRef;
                    $currentRESNum++;
                    $currentResidueRef = {};
                }
                @phospate = ();
                @sugar = ();
                @base = ();
            }
            $lastRes = $curentRes;

        }#end handling last residue


        #handling last chain 
        
        $currentChain =$$atom{"CHAIN_ID"};
      
        
        if($currentChain ne $lastChain ){
            if (exists($currentChainRef->{"BEADS_DNA"}) && scalar(@{$currentChainRef->{"BEADS_DNA"}}) > 0){
            $currentChainRef->{"SIZE"} = scalar(@{$currentChainRef->{"BEADS_DNA"}});  #$currentBeadNum - $firstBeadInCain;
                if ($currentChainRef->{"SIZE"} >0){
                    $chainsDataRef->{$lastChain} = $currentChainRef;
                }
                $currentChainRef =$chainsDataRef->{$currentChain};
            }
            $lastChain = $currentChain;
            $firstBeadInCain = $currentBeadNum+1;
             
        }
        
        #handle the new atom
        
        
        #SKIP protein atoms and hidrogen atoms
        if(length ($$atom{"RESIDUE_TYPE"}) == 3){  
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
    
    #insert last base
    if (scalar(@base) != 0) {
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
            if (keys(%$sugarBead) != 0)
            {
                $currentResidueRef->{"BEADS_DNA"}->{"S"} = $sugarBead;
                push(@{$currentChainRef->{"BEADS_DNA"}},$sugarBead);
                push(@$beadsDataRef, $sugarBead);
            }
            if (keys(%$baseBead) != 0)
            {
                $currentResidueRef->{"BEADS_DNA"}->{"B"} = $baseBead;
                push(@{$currentChainRef->{"BEADS_DNA"}},$baseBead);
                push(@$beadsDataRef, ,$baseBead);
            }
            
            if (keys(%$currentResidueRef) != 0)
            {
                push(@$dnaResiduesDataRef,$currentResidueRef);
                $currentChainRef->{"BASES"}->{$lastAtom->{"RESIDUE_ID"}} = $currentResidueRef;
                $currentRESNum++;
                $currentResidueRef = {};
            }
    }

               
    #insert last chain
    if (exists($currentChainRef->{"BEADS_DNA"}) && scalar(@{$currentChainRef->{"BEADS_DNA"}}) > 0){
        $currentChainRef->{"SIZE"} = scalar(@{$currentChainRef->{"BEADS_DNA"}}); 
        $chainsDataRef->{$currentChain} = $currentChainRef;
        #set chain size
    }
    #inset data to $mdDataRef
    $mdDataRef->{"BEADS_DNA"} = $beadsDataRef;
}  



sub makeDNARes
{
    
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
              $phospateBead->{"RADIUS"} = $executionPreferencesRef->{"DNA_REPULSION_DISTANCE"} ;
            #############################################
    }
#  sugar
    if (scalar(@sugar) != 0)
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
              $sugarBead->{"RADIUS"} = $executionPreferencesRef->{"DNA_REPULSION_DISTANCE"} ;
            ############################################
    }
#   base
    if (scalar(@base) != 0)
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
              $baseBead->{"RADIUS"} = $executionPreferencesRef->{"DNA_REPULSION_DISTANCE"} ;
            ############################################
    }

    return ($phospateBead,$baseBead,$sugarBead,$currentBeadNum);
    
}


sub calcCenter
{
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


sub determineBeadPosition
{
    my ($executionPreferencesRef,$beadRef,$pdbResidueRef,$pdbAtomRef) = @_;
    if ($pdbAtomRef->{"TYPE"} =~ m/CA/i)
    {
        $beadRef->{"X"} = $pdbAtomRef->{"X"};
        $beadRef->{"Y"} = $pdbAtomRef->{"Y"};
        $beadRef->{"Z"} = $pdbAtomRef->{"Z"};
        $beadRef->{"RADIUS"} = $beadRepulsionRadius->{"CA"};
        return;
    }
    
    if ($pdbAtomRef->{"TYPE"} =~ m/CB/i)
    {
        if ($executionPreferencesRef->{"CB_POSITION"} =~ m/CB/i)
        {
          $beadRef->{"X"} = $pdbAtomRef->{"X"};
          $beadRef->{"Y"} = $pdbAtomRef->{"Y"};
          $beadRef->{"Z"} = $pdbAtomRef->{"Z"};
          $beadRef->{"RADIUS"} = $beadRepulsionRadius->{"CB"};
          return;
        }
        if ($executionPreferencesRef->{"CB_POSITION"} =~ m/FHA/i)
        {
          my $farthestAtomType = $farthestHeavyAtomByResidueType->{$pdbAtomRef->{"RESIDUE_TYPE"}};
          
          ##### die if you do not know where to put the CB since FHA dose not exist # added by ohad  ###
          die "farthest havy Atom ".$pdbResidueRef->{"ATOMS_BY_TYPE"}->{$farthestAtomType}." dose not exist for bead
               ".$beadRef->{"RESIDUE_TYPE"}." ".$beadRef->{"RESIDUE_ID"}." ".$beadRef->{"ATOM_TYPE"}." : $farthestAtomType not found\n"
            unless defined($pdbResidueRef->{"ATOMS_BY_TYPE"}->{$farthestAtomType});
          
          my $fhaAtomRef = $pdbResidueRef->{"ATOMS_BY_TYPE"}->{$farthestAtomType};
          my $caAtomRef = $pdbResidueRef->{"ATOMS_BY_TYPE"}->{"CA"};
          my $fhaRadius = MDMath::calculateDistance($caAtomRef,$fhaAtomRef) * $executionPreferencesRef->{"FHA_RADIUS_FACTOR"};
          $beadRef->{"X"} = $fhaAtomRef->{"X"};
          $beadRef->{"Y"} = $fhaAtomRef->{"Y"};
          $beadRef->{"Z"} = $fhaAtomRef->{"Z"};
          $beadRef->{"RADIUS"} = $fhaRadius;          
          return;
        }
        
    }
    
}


sub determineBeads
{
    print "Generating Beads.\n";
    my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
    
    my $beadsDataRef = [];
    my $residuesDataRef = [];
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
      for (my $pdbResidueIter = 0; $pdbResidueIter < scalar(@{$pdbChainRef->{"RESIDUES"}}) ; $pdbResidueIter++)
      {
        $residueIter++;        
        my $currentResidueRef = {};
        my $pdbResidueRef = $pdbChainRef->{"RESIDUES"}->[$pdbResidueIter];
        my $currentResidueID;
        $currentResidueID = $pdbResidueRef->{"RESIDUE_ID"};

        my $pdbAtomsRef = $pdbResidueRef->{"ATOMS"};
        for (my $pdbAtomIter = 0; $pdbAtomIter < scalar(@{$pdbResidueRef->{"ATOMS"}}); $pdbAtomIter++)
        {
          my $pdbAtomRef = $pdbResidueRef->{"ATOMS"}->[$pdbAtomIter];
          my $currentAtomID = $pdbAtomRef->{"ID"};
          #skip if not in range or this is a DNA residue
          if (($currentAtomID < $firstAtomID) or ($currentAtomID > $lastAtomID)
              || length($pdbAtomRef->{"RESIDUE_TYPE"})<3)
          {
            $residueIter--;
            next;
          }

          if (my ($currentAtomType) = $pdbAtomRef->{"TYPE"} =~ m/(CA|CB)/)
          {
              next if (exists($currentResidueRef->{"BEADS"}->{$currentAtomType}));
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
              $beadRef->{"RESIDUE_INDEX"} = $residueIter;
              $beadRef->{"RESIDUE_TYPE"} = $pdbAtomRef->{"RESIDUE_TYPE"};
              determineBeadPosition($executionPreferencesRef,$beadRef,$pdbResidueRef,$pdbAtomRef);
              $beadRef->{"MASS"} = 1.0;
              $beadsDataRef->[$beadIter] = $beadRef;
              $currentResidueRef->{"BEADS"}->{$currentAtomType} = $beadRef;
              $currentResidueRef->{"RESIDUE_ID"} = $currentResidueID;
              $currentResidueRef->{"RESIDUE_INDEX"} = $residueIter;
              $currentResidueRef->{"RESIDUE_TYPE"} = $pdbAtomRef->{"RESIDUE_TYPE"};
              push(@{$currentChainRef->{"BEADS"}},$beadRef);
           }
        }
        if(scalar(keys(%$currentResidueRef))==4){ #residue was not skipped
            $residuesDataRef->[$residueIter] = $currentResidueRef;
            $currentChainRef->{"RESIDUES"}->{$currentResidueID} = $currentResidueRef;
        }
        else
        {
            $residueIter--;
        }
      }
      $currentChainRef->{"SIZE"} = $beadIter-$firstBeadIndex;    
      $chainsDataRef->{$currentChainID} = $currentChainRef;
    }
    $mdDataRef->{"BEADS"} = $beadsDataRef;
    $mdDataRef->{"RESIDUES"} = $residuesDataRef;
    $mdDataRef->{"CHAINS"} = $chainsDataRef;
}

sub determineBonds
{
  print "Calculating Bonds.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  my $bondsRef = [];

  my $mdResiduesRef = $mdDataRef->{"RESIDUES"};
  
  my $bondIter = -1;
  #make CA CB bond ???
  for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef); $residueIter++)
  {
    #skip GLY or non dynamic residues OR missing CA/CB bead
    next if ($mdResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY"
             || isDynamic($mdResiduesRef->[$residueIter]->{"RESIDUE_ID"}) == 0 #not dynamic
             || !exists($mdResiduesRef->[$residueIter]->{"BEADS"}->{"CA"})
             || !exists($mdResiduesRef->[$residueIter]->{"BEADS"}->{"CB"}));
    $bondIter++;
    my $CACBbondRef = {};
    my $currentCABeadRef = $mdResiduesRef->[$residueIter]->{"BEADS"}->{"CA"};
    my $currentCBBeadRef = $mdResiduesRef->[$residueIter]->{"BEADS"}->{"CB"};
    $CACBbondRef->{"FIRST_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
    $CACBbondRef->{"SECOND_BEAD_ID"} = $currentCBBeadRef->{"ATOM_ID"};
    $CACBbondRef->{"DISTANCE"} = MDMath::calculateDistance($currentCABeadRef,$currentCBBeadRef);
    $CACBbondRef->{"COEFFICIENT"} = 100.0;
    $bondsRef->[$bondIter]=$CACBbondRef;
  }
  #make CA CA bond ???
  for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef)-1; $residueIter++)
  {
    #skip if both residues are static
    next if (isDynamic($mdResiduesRef->[$residueIter+1]->{"RESIDUE_ID"})==0
             && isDynamic($mdResiduesRef->[$residueIter]->{"RESIDUE_ID"})==0 );
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
  for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef)-2; $residueIter++){
     #skip DNA residues
     if(length ($mdResiduesRef->[$residueIter]->{"RESIDUE_TYPE"}) < 3){  #alternative $$atom{"RESIDUE_TYPE"} !~ / / ){
           next;
     }
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
  
    for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef)-1; $residueIter++)    {  
         #skip DNA residues
     if(length ($mdResiduesRef->[$residueIter]->{"RESIDUE_TYPE"}) < 3){  #alternative $$atom{"RESIDUE_TYPE"} !~ / / ){
           next;
     }
      if ($mdResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} ne "GLY"
           &&  exists($mdResiduesRef->[$residueIter]->{"BEADS"}->{"CA"})
           && exists($mdResiduesRef->[$residueIter]->{"BEADS"}->{"CB"}))
      {
        $angleIter++;
        my $CACACBangleRef = {};
        my $currentCABeadRef = $mdResiduesRef->[$residueIter] ->{"BEADS"}->{"CA"};
        my $nextCABeadRef = $mdResiduesRef->[$residueIter+1] ->{"BEADS"}->{"CA"};
        my $currentCBBeadRef = $mdResiduesRef->[$residueIter] ->{"BEADS"}->{"CB"};
        $CACACBangleRef->{"FIRST_BEAD_ID"} = $currentCBBeadRef->{"ATOM_ID"};
        $CACACBangleRef->{"SECOND_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
        $CACACBangleRef->{"THIRD_BEAD_ID"} = $nextCABeadRef->{"ATOM_ID"};
        $CACACBangleRef->{"ANGLE"} = MDMath::calculateAngle($currentCBBeadRef,$currentCABeadRef,$nextCABeadRef);
        if ($currentCABeadRef->{"CHAIN_ID"} eq $nextCABeadRef->{"CHAIN_ID"})
        {
           $CACACBangleRef->{"COEFFICIENT"} = 20.0;
        }
        else
        {        
          $CACACBangleRef->{"COEFFICIENT"} = 0.0;
        }
        $anglesRef->[$angleIter]=$CACACBangleRef;
      }
      
      if ($mdResiduesRef->[$residueIter+1]->{"RESIDUE_TYPE"} ne "GLY"
          && exists($mdResiduesRef->[$residueIter+1]->{"BEADS"}->{"CB"}))
      {
        $angleIter++;
        my $CBCACAangleRef = {};
        my $currentCABeadRef = $mdResiduesRef->[$residueIter] ->{"BEADS"}->{"CA"};
        my $nextCABeadRef = $mdResiduesRef->[$residueIter+1] ->{"BEADS"}->{"CA"};
        my $nextCBBeadRef = $mdResiduesRef->[$residueIter+1] ->{"BEADS"}->{"CB"};
        $CBCACAangleRef->{"FIRST_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
        $CBCACAangleRef->{"SECOND_BEAD_ID"} = $nextCABeadRef->{"ATOM_ID"};
        $CBCACAangleRef->{"THIRD_BEAD_ID"} = $nextCBBeadRef->{"ATOM_ID"};
        $CBCACAangleRef->{"ANGLE"} = MDMath::calculateAngle($currentCABeadRef,$nextCABeadRef,$nextCBBeadRef);
        if (($currentCABeadRef->{"CHAIN_ID"} eq $nextCABeadRef->{"CHAIN_ID"}) and ($currentCABeadRef->{"CHAIN_ID"} eq $nextCBBeadRef->{"CHAIN_ID"}))
        {
           $CBCACAangleRef->{"COEFFICIENT"} = 20.0
        }
        else
        {        
          $CBCACAangleRef->{"COEFFICIENT"} = 0.0;
        }
        $anglesRef->[$angleIter]=$CBCACAangleRef;
      }

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
  
  for (my $residueIter = 0; $residueIter < scalar(@$mdResiduesRef)-3; $residueIter++){
    if(length ($mdResiduesRef->[$residueIter]->{"RESIDUE_TYPE"}) < 3){  #skip DNA residues
        next;
     }
   
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

sub determineChirals
{
  print "Calculating chirals.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  my $chiralsRef = [];
  my $chiralsIter = -1;

  my $mdResiduesRef = $mdDataRef->{"RESIDUES"};
  
    if ($mdResiduesRef->[0]->{"RESIDUE_TYPE"} ne "GLY"
           && exists($mdResiduesRef->[0]->{"BEADS"}->{"CB"}))
    {
        $chiralsIter++;
        my $CACACACBDihedralRef = {};
        my $currentCABeadRef = $mdResiduesRef->[0] ->{"BEADS"}->{"CA"};
        my $previousCABeadRef = $mdResiduesRef->[1] ->{"BEADS"}->{"CA"};
        my $nextCABeadRef = $mdResiduesRef->[2] ->{"BEADS"}->{"CA"};
        my $currentCBBeadRef = $mdResiduesRef->[0] ->{"BEADS"}->{"CB"};
        $CACACACBDihedralRef->{"FIRST_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"SECOND_BEAD_ID"} = $previousCABeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"THIRD_BEAD_ID"} = $nextCABeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"FOURTH_BEAD_ID"} = $currentCBBeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"CHIRAL"} = MDMath::calculateTripleProduct($currentCABeadRef,$previousCABeadRef,$nextCABeadRef,$currentCBBeadRef);
        if (($currentCABeadRef->{"CHAIN_ID"} eq $nextCABeadRef->{"CHAIN_ID"})
        and ($currentCABeadRef->{"CHAIN_ID"} eq $previousCABeadRef->{"CHAIN_ID"}))
        {
          $CACACACBDihedralRef->{"COEFFICIENT"} = $executionPreferencesRef->{"CHIRALITY_COEFFICIENT"};
        }
        else
        {           
          $CACACACBDihedralRef->{"COEFFICIENT"} = 0.0;
        }
        $chiralsRef->[$chiralsIter]=$CACACACBDihedralRef;
    }
    for (my $residueIter = 1; $residueIter < scalar(@$mdResiduesRef)-1; $residueIter++)
    {
        #skip for GLY and DNA
      next if ($mdResiduesRef->[$residueIter]->{"RESIDUE_TYPE"} eq "GLY"
                ||length ($mdResiduesRef->[$residueIter]->{"RESIDUE_TYPE"}) < 3
                ||! exists($mdResiduesRef->[$residueIter]->{"BEADS"}->{"CB"}));
      $chiralsIter++;
      my $CACACACBDihedralRef = {};
      my $currentCABeadRef = $mdResiduesRef->[$residueIter] ->{"BEADS"}->{"CA"};
      my $previousCABeadRef = $mdResiduesRef->[$residueIter-1] ->{"BEADS"}->{"CA"};
      my $nextCABeadRef = $mdResiduesRef->[$residueIter+1] ->{"BEADS"}->{"CA"};
      my $currentCBBeadRef = $mdResiduesRef->[$residueIter] ->{"BEADS"}->{"CB"};
      $CACACACBDihedralRef->{"FIRST_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
      $CACACACBDihedralRef->{"SECOND_BEAD_ID"} = $previousCABeadRef->{"ATOM_ID"};
      $CACACACBDihedralRef->{"THIRD_BEAD_ID"} = $nextCABeadRef->{"ATOM_ID"};
      $CACACACBDihedralRef->{"FOURTH_BEAD_ID"} = $currentCBBeadRef->{"ATOM_ID"};
      $CACACACBDihedralRef->{"CHIRAL"} = MDMath::calculateTripleProduct($currentCABeadRef,$previousCABeadRef,$nextCABeadRef,$currentCBBeadRef);
      if (($currentCABeadRef->{"CHAIN_ID"} eq $nextCABeadRef->{"CHAIN_ID"})
      and ($currentCABeadRef->{"CHAIN_ID"} eq $previousCABeadRef->{"CHAIN_ID"}))
      {
        $CACACACBDihedralRef->{"COEFFICIENT"} = $executionPreferencesRef->{"CHIRALITY_COEFFICIENT"};
      }
      else
      {         
        $CACACACBDihedralRef->{"COEFFICIENT"} = 0.0;
      }
      $chiralsRef->[$chiralsIter]=$CACACACBDihedralRef;
    }

    if ($mdResiduesRef->[scalar(@$mdResiduesRef)-1]->{"RESIDUE_TYPE"} ne "GLY"
        && exists($mdResiduesRef->[scalar(@$mdResiduesRef)-1]->{"BEADS"}->{"CB"}))
    {
        $chiralsIter++;
        my $CACACACBDihedralRef = {};
        my $currentCABeadRef = $mdResiduesRef->[scalar(@$mdResiduesRef)-1] ->{"BEADS"}->{"CA"};
        my $previousCABeadRef = $mdResiduesRef->[scalar(@$mdResiduesRef)-3] ->{"BEADS"}->{"CA"};
        my $nextCABeadRef = $mdResiduesRef->[scalar(@$mdResiduesRef)-2] ->{"BEADS"}->{"CA"};
        my $currentCBBeadRef = $mdResiduesRef->[scalar(@$mdResiduesRef)-1] ->{"BEADS"}->{"CB"};
        $CACACACBDihedralRef->{"FIRST_BEAD_ID"} = $currentCABeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"SECOND_BEAD_ID"} = $previousCABeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"THIRD_BEAD_ID"} = $nextCABeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"FOURTH_BEAD_ID"} = $currentCBBeadRef->{"ATOM_ID"};
        $CACACACBDihedralRef->{"CHIRAL"} = MDMath::calculateTripleProduct($currentCABeadRef,$previousCABeadRef,$nextCABeadRef,$currentCBBeadRef);
        if (($currentCABeadRef->{"CHAIN_ID"} eq $nextCABeadRef->{"CHAIN_ID"})
        and ($currentCABeadRef->{"CHAIN_ID"} eq $previousCABeadRef->{"CHAIN_ID"}))
        {
          $CACACACBDihedralRef->{"COEFFICIENT"} = $executionPreferencesRef->{"CHIRALITY_COEFFICIENT"};
        }
        else
        {           
          $CACACACBDihedralRef->{"COEFFICIENT"} = 0.0;
        }
        $chiralsRef->[$chiralsIter]=$CACACACBDihedralRef;
    }
  
  $mdDataRef->{"CHIRALS"} = $chiralsRef;
  
}

sub isHydrogen
{
  my $atomType = $_[0];
  return ($atomType =~ m/H/);

}

sub determineBead
{
    my $atomType = $_[0];
    if (defined($caBeadAtoms->{$atomType}))
    {
      return "CA";
    }
    else
    {
      return "CB";
    }
    
}

sub shouldContact
{
    my ($firstResidueIndex, $firstResidueAtomType, $firstResidueAtomBead,  $secondResidueIndex,$secondResidueAtomType,$secondResidueAtomBead) = @_;
    return 0 if (isHydrogen($firstResidueAtomType) or isHydrogen($secondResidueAtomType));
   return (abs($firstResidueIndex-$secondResidueIndex) >= $minimumContactDistance->{$firstResidueAtomBead."_".$secondResidueAtomBead});
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
    my $firstResidueAtomBead = determineBead($firstResidueAtomType);
    my $secondResidueAtomBead = determineBead($secondResidueAtomType);
#    next if (isHydrogen($firstResidueAtomType) or isHydrogen($secondResidueAtomType));
#   if ((abs($secondResidueNumber-$firstResidueNumber) >= 4))
    if (shouldContact($firstResidueNumber,$firstResidueAtomType,$firstResidueAtomBead,$secondResidueNumber,$secondResidueAtomType,$secondResidueAtomBead))
    {
      if ($secondResidueNumber > $firstResidueNumber &&
            (isDynamic($secondResidueNumber)==1
             || isDynamic($firstResidueNumber)==1 ))
      {
        addContact($contactsDataRef,$firstResidueNumber,$firstResidueAtomBead,$secondResidueNumber,$secondResidueAtomBead,$distance);
      }
      else
      {
        if (!contactDefined($contactsDataRef,$secondResidueNumber,$secondResidueAtomBead,$firstResidueNumber,$firstResidueAtomBead)&&
            (isDynamic($secondResidueNumber)==1
             || isDynamic($firstResidueNumber)==1 ))
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
  print "Translating CSU contacts to bead contacts for chain ".$currentChainID.".\n";
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
      $contactsDataRef->{$firstBeadRef->{"ATOM_ID"}.".".$secondBeadRef->{"ATOM_ID"}}->{"COEFFICIENT"} = $executionPreferencesRef->{"CONTACTS_COEFFICIENT"};
    }
  }
}

sub determineContacts
{
  print "Calculating CSU contacts.\n";
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  my $scriptDirectory = $executionPreferencesRef->{"SCRIPT_DIRECTORY"};
  require $utilDirectory."PDBHandling.pm";
  my $contactsDataRef = {};

  if ($executionPreferencesRef->{"INTER_CHAIN_CONTACTS"} =~ m/NO/)
  {
    foreach my $pdbChainRef (@{$pdbDataRef->{"CHAINS"}})
    {
      my $chainContactsDataRef = {};
      my $currentChainID = $pdbChainRef->{"CHAIN_ID"};
      PDBHandling::writeChainAsPDBFile($pdbChainRef,"csu.input.".$currentChainID);
      foreach my $pdbResidueRef (@{$pdbChainRef->{"RESIDUES"}})
      {
        next if(length ($pdbResidueRef->{"RESIDUE_TYPE"}) < 3); #skip DNA residues
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
 
    PDBHandling::writeAllChainsAsSingleChainPDBFile($pdbDataRef->{"RESIDUES"},"csu.input.A");
    for (my $residueIter = 1; $residueIter <= scalar(@{$pdbDataRef->{"RESIDUES"}}); $residueIter++)
    {
      #skip if residue is static
      #print "$residueIter\n";
      #sleep(0.5);
      next if (isDynamic($pdbDataRef->{RESIDUES}[$residueIter-1]->{"RESIDUE_ID"})==0);
      print isDynamic($pdbDataRef->{RESIDUES}[$residueIter]->{"RESIDUE_ID"})==0;
      `$scriptDirectory/resc csu.input.A $residueIter >> csu.output.A.$residueIter`;
      _determineContacts("csu.output.A.".$residueIter,$residueIter,$CSUContactsDataRef);
      unlink("csu.output.A.".$residueIter);
    }
    
    unlink("csu.input.A");
    translateContactsToBeads($CSUContactsDataRef,$mdDataRef,$executionPreferencesRef,$contactsDataRef);
  }
  $mdDataRef->{"CONTACTS"} = $contactsDataRef;
}

sub determineDNAContacts {
  my($executionPreferencesRef,$pdbDataRef, $mdDataRef) = @_;
  my $DNAContRegionNum=0;
  my @DNAContactRegions;
  if ($executionPreferencesRef->{HAS_DNA_CONTACT_REGION}=~ /YES/)
  {
    @DNAContactRegions=split(/\,/,$executionPreferencesRef->{"DNA_CONTACT_RANGE"});
    $DNAContRegionNum=scalar(@DNAContactRegions);
  }
  my @DNAChainID;
  my $counter=1;
  $DNAChainID[0]=$mdDataRef->{"BEADS_DNA"}->[0]->{"CHAIN_ID"};
  foreach my $i (@{$mdDataRef->{"BEADS_DNA"}})
  {
    if ($i->{"CHAIN_ID"} ne $DNAChainID[$counter-1])
    {
    $DNAChainID[$counter]=$i->{"CHAIN_ID"};
    $counter++;
    }
  }
  
  $counter=0;

  foreach my $c (@{$pdbDataRef->{"CHAINS"}})
  {
    my $BeadType="CA";
    my $control=0;
    foreach my $dnachain (@DNAChainID)
    {
     if ($c->{"CHAIN_ID"} eq $dnachain) {$control=1;}
    }
    if ($control == 0)
    {
      foreach my $dnachain (@DNAChainID)
      {
      my $firstresID=0;
      my ($firstbead,$secondbead,$firstbeadID,$secondbeadID);
      my $secondresID=0;
      my $DNABead='D';
      my $ChainLength=@{$c->{"ATOMS"}};
      my $DNALength=@{$pdbDataRef->{"CHAINS_BY_ID"}->{$dnachain}->{"ATOMS"}};
      for (my $i = 0; $i < $ChainLength; $i++)
      {
        my $region_control=0;
        for (my $rc=0; $rc < $DNAContRegionNum; $rc=$rc+2)
        {
            if (($c->{"ATOMS"}->[$i]->{"RESIDUE_ID"} < $DNAContactRegions[$rc]) || ($c->{"ATOMS"}->[$i]->{"RESIDUE_ID"} > $DNAContactRegions[$rc+1]))
            { $region_control=1;}
        }
        if ($region_control==0)
        {
        
          my $dnaatom=$pdbDataRef->{"CHAINS_BY_ID"}->{$dnachain};
          for (my $j = 0; $j < $DNALength; $j++)
          {
            if ($c->{"ATOMS"}->[$i]->{"TYPE"} =~ m/O|N|C|S/)
            {
             my $dist = MDMath::calculateDistance($c->{"ATOMS"}->[$i],$dnaatom->{"ATOMS"}->[$j]);
              if ($dist < 4.0)
              {
                if ($c->{"ATOMS"}->[$i]->{"TYPE"} =~ m/ N | O /)
                {$BeadType="CA";}
                else
                {$BeadType="CB";}
                 #Determining DNA bead type
                  if ($dnaatom->{"ATOMS"}->[$j]->{"TYPE"} eq "P")
                  {
                    $DNABead='P';
                  }
                  elsif ($dnaatom->{"ATOMS"}->[$j]->{"TYPE"} =~ m/\*/)
                  {
                    $DNABead= 'S';
                  }
                  else
                  {
                    $DNABead= 'B';
                  }
                  
                 $firstresID=$c->{"ATOMS"}->[$i]->{"RESIDUE_ID"};
                 $firstbead=$mdDataRef->{"CHAINS"}->{$c->{"CHAIN_ID"}}->{"RESIDUES"}->{$firstresID}->{"BEADS"}->{"$BeadType"};
                 $firstbeadID=$firstbead->{"ATOM_ID"};
                 $secondresID=$dnaatom->{"ATOMS"}->[$j]->{"RESIDUE_ID"};
                 $secondbead=$mdDataRef->{"CHAINS"}->{$dnachain}->{"BASES"}->{$secondresID}->{"BEADS_DNA"}->{$DNABead};
                 $secondbeadID=$secondbead->{"ATOM_ID"};
                 
                 
                if (!defined($mdDataRef->{"DNA.CONTACTS"}->{$firstbeadID.".".$secondbeadID}))
                 {
                  $mdDataRef->{"DNA_CONTACTS"}->[$counter]->{"DISTANCE"}=$dist; 
                  $mdDataRef->{"DNA_CONTACTS"}->[$counter]->{"DNARES"}=$secondbeadID;
                  $mdDataRef->{"DNA_CONTACTS"}->[$counter]->{"PROTRES"}=$firstbeadID;
                  $mdDataRef->{"DNA_CONTACTS"}->[$counter]->{"DNA_BEAD_TYPE"}=$DNABead;
                  $mdDataRef->{"DNA.CONTACTS"}->{$firstbeadID.".".$secondbeadID}=1;
                  $mdDataRef->{"DNA_CONTACTS"}->[$counter]->{"BEAD_DISTANCE"}= (MDMath::calculateDistance($firstbead,$secondbead))**2;
                  $counter++;
                  }
               }
            }
         }
        }
      }
      }
    }
  }
}

sub calculateRepulsionDistance
{
  my ($executionPreferencesRef,$firstBeadRef,$secondBeadRef) = @_;
  my $repulsionDistance = (($executionPreferencesRef->{"REPULSION_RADIUS_FACTOR"}*($firstBeadRef->{"RADIUS"} + $secondBeadRef->{"RADIUS"}))**2);
  return $repulsionDistance;
}

sub calculateEllipsoidRepulsionRadius
{
  my ($executionPreferencesRef,$caBeadRef,$cbBeadRef,$intruderBeadRef) = @_;
  my $caBeadRadius = $caBeadRef->{"RADIUS"}*$executionPreferencesRef->{"REPULSION_RADIUS_FACTOR"};
  my $cbBeadRadius = $cbBeadRef->{"RADIUS"}*$executionPreferencesRef->{"REPULSION_RADIUS_FACTOR"};
  my $intruderBeadRadius = $intruderBeadRef->{"RADIUS"}*$executionPreferencesRef->{"REPULSION_RADIUS_FACTOR"};
  my $minimalRadius = ($caBeadRadius < $cbBeadRadius ?
                       ($caBeadRadius+$intruderBeadRadius)*$executionPreferencesRef->{"REPULSION_ELLIPSOIDS_FACTOR"} :
                       ($cbBeadRadius+$intruderBeadRadius)*$executionPreferencesRef->{"REPULSION_ELLIPSOIDS_FACTOR"} );

  my $repulsionRadius = sqrt($minimalRadius**2 + $cbBeadRef->{"RADIUS"}**2)*2;
  return $repulsionRadius;
}

sub addRepulsion
{
    my ($executionPreferencesRef,$repulsionsDataRef,$firstBeadRef,$secondBeadRef) = @_;
    my $newRepulsionRef = {};
    $newRepulsionRef->{"FIRST_BEAD_ID"} = $firstBeadRef->{"ATOM_ID"};
    $newRepulsionRef->{"SECOND_BEAD_ID"} = $secondBeadRef->{"ATOM_ID"};
    $newRepulsionRef->{"DISTANCE"} = calculateRepulsionDistance($executionPreferencesRef,$firstBeadRef,$secondBeadRef);
    
    $newRepulsionRef->{"COEFFICIENT"} = 1.0;
    $repulsionsDataRef->{10000*$firstBeadRef->{"ATOM_ID"}+$secondBeadRef->{"ATOM_ID"}} = $newRepulsionRef;
}

sub addRepulsionForDNA
{
    my ($executionPreferencesRef,$mdDataRef,$repulsionsDataRef,$firstBeadRef,$secondBeadRef) = @_;
    my $newRepulsionRef = {};
    if (!defined($mdDataRef->{"DNA.CONTACTS"}->{$firstBeadRef->{"ATOM_ID"}.".".$secondBeadRef->{"ATOM_ID"}}))
    {
    $newRepulsionRef->{"FIRST_BEAD_ID"} = $firstBeadRef->{"ATOM_ID"};
    $newRepulsionRef->{"SECOND_BEAD_ID"} = $secondBeadRef->{"ATOM_ID"};
    #if (($executionPreferencesRef->{"CB_POSITION"} =~ m/FHA/i) & ($firstBeadRef->{"ATOM_TYPE"} =~ m/CB/i))
    if ($executionPreferencesRef->{"CB_POSITION"} =~ m/FHA/i)
    {
    $newRepulsionRef->{"DISTANCE"} = (2*($executionPreferencesRef->{"DNA_REPULSION_DISTANCE"})*($executionPreferencesRef->{"DNA_FHA_REPULSION_FACTOR"}))**2;
    }
    else
    {
    $newRepulsionRef->{"DISTANCE"} = (($executionPreferencesRef->{"DNA_REPULSION_DISTANCE"}+$firstBeadRef->{"RADIUS"})**2);
    }
    $newRepulsionRef->{"COEFFICIENT"} = 1.0;
    $repulsionsDataRef->{10000*$firstBeadRef->{"ATOM_ID"}+$secondBeadRef->{"ATOM_ID"}} = $newRepulsionRef;
    }
}

sub addEllipsoidRepulsion
{
    my ($executionPreferencesRef,$ellipsoidRepulsionsDataRef,$caBeadRef,$cbBeadRef,$intruderBeadRef) = @_;
    my $newRepulsionRef = {};
    $newRepulsionRef->{"CA_BEAD_ID"} = $caBeadRef->{"ATOM_ID"};
    $newRepulsionRef->{"CB_BEAD_ID"} = $cbBeadRef->{"ATOM_ID"};
    $newRepulsionRef->{"INTRUDER_BEAD_ID"} = $intruderBeadRef->{"ATOM_ID"};
    $newRepulsionRef->{"DISTANCE"} = calculateEllipsoidRepulsionRadius($executionPreferencesRef,$caBeadRef,$cbBeadRef,$intruderBeadRef);
    $newRepulsionRef->{"COEFFICIENT"} = $executionPreferencesRef->{"REPULSION_ELLIPSOIDS_COEFFICIENT"};
    $ellipsoidRepulsionsDataRef->{100000000*$caBeadRef->{"ATOM_ID"}+10000*$cbBeadRef->{"ATOM_ID"}+$intruderBeadRef->{"ATOM_ID"}} = $newRepulsionRef;
}


sub beadsShouldRepulse
{
    my ($executionPreferencesRef,$firstBeadRef,$secondBeadRef) = @_;
    my $firstResidueID = $firstBeadRef->{"RESIDUE_ID"};
    my $firstBeadType = $firstBeadRef->{"ATOM_TYPE"};
    my $secondResidueID = $secondBeadRef->{"RESIDUE_ID"};
    my $secondBeadType = $secondBeadRef->{"ATOM_TYPE"};

#      print $firstResidueID,$secondResidueID,$firstBeadType,$secondBeadType, $minimumRepulsionDistance_short->{$firstBeadType."_".$secondBeadType},"ok \n";
#      sleep(1);


    if ($executionPreferencesRef->{"REPULSION_NEIGHBORS"} =~ m/NORMAL/i )
    {
      return (($secondResidueID-$firstResidueID) >= $minimumRepulsionDistance_short->{$firstBeadType."_".$secondBeadType});
    }
    if ($executionPreferencesRef->{"REPULSION_NEIGHBORS"} =~ m/FAR/i )
    {
      return (($secondResidueID-$firstResidueID) >= $minimumRepulsionDistance_long->{$firstBeadType."_".$secondBeadType});
    }
}

sub ellipsoidShouldRepulse
{
    my ($executionPreferencesRef,$caBeadRef,$cbBeadRef,$intruderBeadRef) = @_;
    my $caChainID = $caBeadRef->{"CHAIN_ID"};
    my $residueID = $caBeadRef->{"RESIDUE_ID"};
    my $intruderChainID = $intruderBeadRef->{"RESIDUE_ID"};
    my $intruderResidueID = $intruderBeadRef->{"RESIDUE_ID"};
    my $intruderBeadType = $intruderBeadRef->{"ATOM_TYPE"};
    if ($executionPreferencesRef->{"REPULSION_NEIGHBORS"} =~ m/NORMAL/i )
    {
      my $caRepulsionDistance = $minimumRepulsionDistance_short->{"CA_".$intruderBeadType};
      my $cbRepulsionDistance = $minimumRepulsionDistance_short->{"CB_".$intruderBeadType};
      my $minimumRepulsionDistance = ( $caRepulsionDistance < $cbRepulsionDistance ? $caRepulsionDistance : $cbRepulsionDistance);
      return ((abs($intruderResidueID-$residueID) >= $minimumRepulsionDistance) or ($caChainID ne $intruderChainID));

    }
    if ($executionPreferencesRef->{"REPULSION_NEIGHBORS"} =~ m/FAR/i )
    {
      my $caRepulsionDistance = $minimumRepulsionDistance_long->{"CA_".$intruderBeadType};
      my $cbRepulsionDistance = $minimumRepulsionDistance_long->{"CB_".$intruderBeadType};
      my $minimumRepulsionDistance = ( $caRepulsionDistance < $cbRepulsionDistance ? $caRepulsionDistance : $cbRepulsionDistance);
      return ((abs($intruderResidueID-$residueID) >= $minimumRepulsionDistance) or ($caChainID ne $intruderChainID));
    }
}


sub beadContactDefined
{
  my ($contactsRef,$firstBeadRef,$secondBeadRef) = @_;
  return (defined($contactsRef->{$firstBeadRef->{"ATOM_ID"}.".".$secondBeadRef->{"ATOM_ID"}}));
}


sub determineRepulsionsForBead
{
  my ($executionPreferencesRef,$repulsionsDataRef,$mdDataRef,$firstBeadRef) = @_;
  #repultion from protein
  for (my $residueIter = $firstBeadRef->{"RESIDUE_INDEX"}+1; $residueIter < scalar(@{$mdDataRef->{"RESIDUES"}}); $residueIter++)
  {
    
    foreach my $secondBeadRef (values(%{$mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}}))
    {
	
      if (beadsShouldRepulse($executionPreferencesRef,$firstBeadRef,$secondBeadRef)
      and (!beadContactDefined($mdDataRef->{"CONTACTS"},$firstBeadRef,$secondBeadRef)))
      {
#	print $firstBeadRef->{"RESIDUE_ID"},$secondBeadRef->{"RESIDUE_ID"},"ok1\n";
        addRepulsion($executionPreferencesRef,$repulsionsDataRef,$firstBeadRef,$secondBeadRef);
      }
    }
  }
  
  #repultion from DNA
  # this assumes DNA is static and no DNA -DNA repultions exist !!!
  for (my $secondBead = 0; $secondBead < scalar(@{$mdDataRef->{"BEADS_DNA"}}); $secondBead++)
  {
         my $secondBeadRef = $mdDataRef->{"BEADS_DNA"}->[$secondBead];
       if (!beadContactDefined($mdDataRef->{"CONTACTS"},$firstBeadRef,$secondBeadRef) and (!defined($mdDataRef->{"DNA.CONTACTS"}->{$firstBeadRef->{"RESIDUE_ID"}.".".$secondBeadRef->{"RESIDUE_ID"}})))
        {
            addRepulsionForDNA($executionPreferencesRef,$mdDataRef,$repulsionsDataRef,$firstBeadRef,$secondBeadRef);
        }
  
    }
  
}

sub determineRepulsions
{
  print "Calculating Repulsions.\n";
  my ($executionPreferencesRef, $mdDataRef) = @_;
  my $repulsionsDataRef = {};
  
  for (my $firstBeadIter = 0; $firstBeadIter < scalar(@{$mdDataRef->{"BEADS"}}); $firstBeadIter++)
  {
    my $firstBeadRef = $mdDataRef->{"BEADS"}->[$firstBeadIter];
    determineRepulsionsForBead($executionPreferencesRef,$repulsionsDataRef,$mdDataRef,$firstBeadRef);
  }
    
    
  $mdDataRef->{"REPULSIONS"} = $repulsionsDataRef;
}

sub determineEllipsoidRepulsionsForResidue
{
  my ($executionPreferencesRef,$ellipsoidRepulsionsDataRef,$mdDataRef,$residueRef,$currentResidueIter) = @_;
  my $caBeadRef = $residueRef->{"BEADS"}->{"CA"};
  my $cbBeadRef = $residueRef->{"BEADS"}->{"CB"};

  for (my $residueIter = 0; $residueIter < scalar(@{$mdDataRef->{"RESIDUES"}}); $residueIter++)
  {
    next if ($currentResidueIter == $residueIter);
    foreach my $intruderBeadRef (values(%{$mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}}))
    {
      if (ellipsoidShouldRepulse($executionPreferencesRef,$caBeadRef,$cbBeadRef,$intruderBeadRef))
      {
        addEllipsoidRepulsion($executionPreferencesRef,$ellipsoidRepulsionsDataRef,$caBeadRef,$cbBeadRef,$intruderBeadRef);
      }
    }
  }
}

sub determineEllipsoidRepulsions
{
  print "Calculating Ellipsoid Repulsions.\n";
  my ($executionPreferencesRef, $mdDataRef) = @_;
  my $ellipsoidRepulsionsDataRef = {};
  
  for (my $firstResidueIter = 0; $firstResidueIter < scalar(@{$mdDataRef->{"RESIDUES"}}); $firstResidueIter++)
  {
    if ($mdDataRef->{"RESIDUES"}->[$firstResidueIter]->{"RESIDUE_TYPE"} ne "GLY"
        && exists($mdDataRef->{"RESIDUES"}->[$firstResidueIter]->{"BEADS"}->{"CB"}))
    {
      determineEllipsoidRepulsionsForResidue($executionPreferencesRef,
                                             $ellipsoidRepulsionsDataRef,
                                             $mdDataRef,
                                             $mdDataRef->{"RESIDUES"}->[$firstResidueIter],
                                             $firstResidueIter);
    }
  }    
  
  $mdDataRef->{"ELLIPSOID_REPULSIONS"} = $ellipsoidRepulsionsDataRef;
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

  for (my $residueIter = 0; $residueIter < scalar(@{$mdDataRef->{"RESIDUES"}}); $residueIter++)
  {
    my $currentResidueType = $mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}->{"CA"}->{"RESIDUE_TYPE"};
    my $currentCharge = determineCharge($currentResidueType);
    if ($currentCharge != 0)
    {
        my $currentCBBeadID ;
        if (exists  ($mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}->{"CB"})){
            $currentCBBeadID =  $mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}->{"CB"}->{"ATOM_ID"};
        }
        else{
            $currentCBBeadID =  $mdDataRef->{"RESIDUES"}->[$residueIter]->{"BEADS"}->{"CA"}->{"ATOM_ID"};
        }
        $electrostaticsRef->{$currentCBBeadID} = $currentCharge;
        
        
    }
  }
  #add electorstatic for DNA phospates
    for (my $baseIter = 0; $baseIter < scalar(@{$mdDataRef->{"BASES"}}); $baseIter++){
        if (exists  ($mdDataRef->{"BASES"}->[$baseIter]->{"BEADS_DNA"}->{"P"})){
            my $currentPid =  $mdDataRef->{"BASES"}->[$baseIter]->{"BEADS_DNA"}->{"P"}->{"ATOM_ID"};
            $electrostaticsRef->{$currentPid} = -1.0;
        }
        
    }
  $mdDataRef->{"ELECTROSTATICS"} = $electrostaticsRef;

}

sub printMDInputFile
{
  my ($executionPreferencesRef,$pdbDataRef, $mdDataRef,@ChainOrder) = @_;
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

  if ($executionPreferencesRef->{"CHIRALITY"} =~ m/YES/i)
  {
    my $mdChiralsRef = $mdDataRef->{"CHIRALS"};
    printf MD_INPUT_FILE_HANDLE ("%12d chiral constraints.\n",scalar(@$mdChiralsRef));
    for (my $chiralIter = 0; $chiralIter < scalar(@$mdChiralsRef); $chiralIter++)
    {
      printf MD_INPUT_FILE_HANDLE ("%5d%5d%5d%5d%5d%8.3f%8.3f\n",
                                 $chiralIter+1,
                                 $mdChiralsRef->[$chiralIter]->{"FIRST_BEAD_ID"},
                                 $mdChiralsRef->[$chiralIter]->{"SECOND_BEAD_ID"},
                                 $mdChiralsRef->[$chiralIter]->{"THIRD_BEAD_ID"},
                                 $mdChiralsRef->[$chiralIter]->{"FOURTH_BEAD_ID"},
                                 $mdChiralsRef->[$chiralIter]->{"CHIRAL"},
                                 $mdChiralsRef->[$chiralIter]->{"COEFFICIENT"});
    }
  }

  my $mdContactsRef = $mdDataRef->{"CONTACTS"};
  my $contacts=scalar(keys(%$mdContactsRef));
  if ($executionPreferencesRef->{"HAS_PROTEIN_DNA_CONTACTS"} =~ m/YES/i)
    {
    $contacts = $contacts + @{$mdDataRef->{"DNA_CONTACTS"}};
    }
  printf MD_INPUT_FILE_HANDLE ("%12d contacts.\n",$contacts);
  my $contactIter = 0;
  foreach my $contactKey (sort { 1000000*$a <=> 1000000*$b } keys(%$mdContactsRef))
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
  #Print DNA-PROT contacts
  foreach my $dnacontact (@{$mdDataRef->{"DNA_CONTACTS"}})
  {
    $contactIter++;
    printf MD_INPUT_FILE_HANDLE ("%5d%5d%5d%10.3f%9.6f\n",
                                 $contactIter,
                                 $dnacontact->{"PROTRES"},
                                 $dnacontact->{"DNARES"},
                                 $dnacontact->{"BEAD_DISTANCE"},
                                 $executionPreferencesRef->{"CONTACTS_COEFFICIENT_DNA"});
  }

  my $mdRepulsionsRef = $mdDataRef->{"REPULSIONS"};
  printf MD_INPUT_FILE_HANDLE ("%12d repulsive pairs.\n",scalar(keys(%$mdRepulsionsRef)));
  my $repulsionIter = 0;
#  foreach my $repulsionKey ( sort contactCompare keys(%$mdRepulsionsRef))
   foreach my $repulsionKey ( sort { 10000*$a <=> 10000*$b } keys(%$mdRepulsionsRef))

  {
    $repulsionIter++;
    printf MD_INPUT_FILE_HANDLE ("%8d%5d%5d%10.3f%10.3f\n",
                                 $repulsionIter,
                                 $mdRepulsionsRef->{$repulsionKey}->{"FIRST_BEAD_ID"},
                                 $mdRepulsionsRef->{$repulsionKey}->{"SECOND_BEAD_ID"},
                                 $mdRepulsionsRef->{$repulsionKey}->{"DISTANCE"},
                                 $mdRepulsionsRef->{$repulsionKey}->{"COEFFICIENT"});
  }

  if ($executionPreferencesRef->{"CB_POSITION"}  =~ m/FHA/i)
  {
    if ($executionPreferencesRef->{"REPULSION_ELLIPSOIDS"} =~ m/YES/i)
    {
      my $mdEllipsoidRepulsionsRef = $mdDataRef->{"ELLIPSOID_REPULSIONS"};
      printf MD_INPUT_FILE_HANDLE ("%12d ellipsoid repulsive trios.\n",scalar(keys(%$mdEllipsoidRepulsionsRef)));
      my $ellipsoidRepulsionIter = 0;
      foreach my $repulsionKey ( sort { 10000*$a <=> 10000*$b } keys(%$mdEllipsoidRepulsionsRef))
      {
        $ellipsoidRepulsionIter++;
        printf MD_INPUT_FILE_HANDLE ("%8d%5d%5d%5d%10.3f%10.3f\n",
                                   $ellipsoidRepulsionIter,
                                   $mdEllipsoidRepulsionsRef->{$repulsionKey}->{"CA_BEAD_ID"},
                                   $mdEllipsoidRepulsionsRef->{$repulsionKey}->{"CB_BEAD_ID"},
                                   $mdEllipsoidRepulsionsRef->{$repulsionKey}->{"INTRUDER_BEAD_ID"},
                                   $mdEllipsoidRepulsionsRef->{$repulsionKey}->{"DISTANCE"},
                                   $mdEllipsoidRepulsionsRef->{$repulsionKey}->{"COEFFICIENT"});
      }
    }
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

  printf MD_INPUT_FILE_HANDLE ("%12d atom positions.\n",(scalar(@{$mdDataRef->{"BEADS"}})+ scalar(@{$mdDataRef->{"BEADS_DNA"}})) );
  printf MD_INPUT_FILE_HANDLE ("%12d chains.\n",scalar(keys(%{$mdDataRef->{"CHAINS"}})));

    
  my $chainIter = 0;
 
 # foreach my $mdChainID (sort keys (%{$mdDataRef->{"CHAINS"}}))
   
   foreach my $mdChainID (@ChainOrder)
   {
    $chainIter++;
    printf MD_INPUT_FILE_HANDLE ("%12d is the size of chain%12d.\n",
                                 $mdDataRef->{"CHAINS"}->{$mdChainID}->{"SIZE"},
                                 $chainIter);

   }
  
  $chainIter = 0;
 
 # foreach my $mdChainID (sort keys (%{$mdDataRef->{"CHAINS"}}))
  foreach my $mdChainID (@ChainOrder)
  {
    $chainIter++;
    my $mdResiduesRef = $mdDataRef->{"CHAINS"}->{$mdChainID}->{"RESIDUES"};
    my $mdBasesRef = $mdDataRef->{"CHAINS"}->{$mdChainID}->{"BASES"};

    #printf MD_INPUT_FILE_HANDLE ("%12d is the size of chain%12d.\n",
    #                             $mdDataRef->{"CHAINS"}->{$mdChainID}->{"SIZE"},
    #                             $chainIter);
    my $beadIter = 0;
    my @mdResidues = (sort {$a <=> $b} keys(%$mdResiduesRef));
    my @mdBases = (sort {$a <=> $b} keys(%$mdBasesRef));
    my @mdResBase = (@mdResidues,@mdBases);
   # foreach my $mdResidueID (sort {$a <=> $b} keys(%$mdResiduesRef))
   foreach my $mdResidueID (@mdResBase)
    {
        my @printOrder =("CA","CB","P","S","B");
        foreach my $bead (@printOrder)
           {
            my $mdBeadRef;
            if (exists($mdResiduesRef->{$mdResidueID}->{"BEADS"}->{$bead})){
                $beadIter++;
                $mdBeadRef = $mdResiduesRef->{$mdResidueID}->{"BEADS"}->{$bead};
            }
            elsif (exists($mdBasesRef->{$mdResidueID}->{"BEADS_DNA"}->{$bead})){
                $beadIter++;
                $mdBeadRef = $mdBasesRef->{$mdResidueID}->{"BEADS_DNA"}->{$bead};
            }
            else {next;}
            printf MD_INPUT_FILE_HANDLE ("%5d%4d %2s %3s%8.3f%8.3f%8.3f%8.3f\n",
                                     #$beadIter,
                                     $mdBeadRef->{"ATOM_ID"},
                                     $mdBeadRef->{"RESIDUE_ID"},
                                     $bead,
                                     $mdBeadRef->{"RESIDUE_TYPE"},
                                     $mdBeadRef->{"X"},
                                     $mdBeadRef->{"Y"},
                                     $mdBeadRef->{"Z"},
                                     $mdBeadRef->{"MASS"});
            }
     }
  }
  close(MD_INPUT_FILE_HANDLE);
}


sub setDynamicRange
{
    my ($executionPreferencesRef,$mdDataRef) = @_;
    my $mdResiduesRef = $mdDataRef->{"RESIDUES"};
    if ($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/)
    {
        my $DynamicRangeList = $executionPreferencesRef->{"DYNAMIC_ATOM_RANGE"};
        @DynamicRange = split /,/, $DynamicRangeList;
        my $DynamicRangeLong = scalar(@DynamicRange);
        my $counter = 0;
        while ( $counter < $DynamicRangeLong)
        {
            if ($DynamicRange[$counter] != 1) {
                $DynamicRange[$counter] = $DynamicRange[$counter] -1;
            }
            if ($DynamicRange[$counter+1] != scalar(@$mdResiduesRef)) {
                $DynamicRange[$counter+1] = $DynamicRange[$counter+1] + 1;
            }
        $counter = $counter +2;
        }
  }
    else {
    $DynamicRange[0] = 0;
    $DynamicRange[1] = scalar(@$mdResiduesRef);
  }
    return @DynamicRange;
}
sub isDynamic
{
    my ($resID) = @_;
    for(my $i=0;($i+1)<scalar(@DynamicRange);$i=$i+2){
        if ($DynamicRange[$i]<=$resID && $DynamicRange[$i+1]>=$resID){
        return 1;
        }
    }
    return 0;
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
# {BASES} -> array[dnaRESIDUE]
# {RESIDUES} -> array[RESIDUE]
#
#
#CHAINS
#   {A}->CHAIN
#   {B}->CHAIN
#
##CHAIN
#   {BEADS} ->array[BEAD]
#   {BEADS_DNA}->array[BEAD](containing dna beads)
#   {RESIDUES} -> array[RESIDUE]
#   {BASES}-> array[RESIDUE](containing dna RESIDUES)
#   size    scalar
#
#RESIDUE
# {BEADS}->array[BEAD]
# {BEADS_DNA}->array[BEAD]
# {RESIDUE_ID}
# {RESIDUE_TYPE}
#
##{BASE} ->
#      {BEADS_DNA}->array[BEAD]
#      {RESIDUE_ID}
#      {RESIDUE_TYPE}
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
