package MDInputHandling;
use strict;
use warnings;
use Math::Trig;

use base 'Exporter';
our @EXPORT_OK = qw(readContacts readESContacts);

sub readContacts
{
  my ($beadsDataRef,$executionPreferencesRef) = @_;
  my $mdInputFileName = $executionPreferencesRef->{"MD_INPUT_FILE"};
  my $numberOfContacts;
  my $contactsRef = [];
  my $contactMapFound = 0;
  
  open (MD_INPUT_FILE_HANDLE, $mdInputFileName) or die "Error reading from ".$mdInputFileName.": $!\n";

  my $currentLine;
  while (!$contactMapFound)
  {
    $currentLine = <MD_INPUT_FILE_HANDLE>;
    if ($currentLine =~ m/\s*(\d+)\s+contacts\./)
    {
      ($numberOfContacts) = $currentLine =~ m/\s*(\d+)\s+contacts\./;
      $contactMapFound = 1;
    }
  }
  for (my $contactIter = 0 ; $contactIter < $numberOfContacts; $contactIter++)
  {
    $currentLine = <MD_INPUT_FILE_HANDLE>;
    my ($firstBeadID, $secondBeadID, $contactDistance)
    = $currentLine =~ m/^\s*\d+\s+(\d+)\s+(\d+)\s+(\d+\.+\d+)/;
    $contactsRef->[$contactIter]->{"ID1"} = $firstBeadID-1;
    $contactsRef->[$contactIter]->{"TYPE1"} = $beadsDataRef->[$firstBeadID-1]->{"TYPE"};    
    $contactsRef->[$contactIter]->{"ID2"} = $secondBeadID-1;
    $contactsRef->[$contactIter]->{"TYPE2"} = $beadsDataRef->[$secondBeadID-1]->{"TYPE"};
    $contactsRef->[$contactIter]->{"DISTANCE"} = sqrt($contactDistance);    
  }
  close (MD_INPUT_FILE_HANDLE);
  
  return $contactsRef;
}

sub readESContacts
{
  my ($beadsDataRef,$executionPreferencesRef) = @_;
  my $mdInputFileName = $executionPreferencesRef->{"MD_INPUT_FILE"};
  my $numberOfESbeads;
  my $esBeadsRef = [];
  my $contactsRef = [];
  my $ESFound = 0;
  
  open (MD_INPUT_FILE_HANDLE, $mdInputFileName) or die "Error reading from ".$mdInputFileName.": $!\n";

  my $currentLine;
  while (!$ESFound)
  {
    $currentLine = <MD_INPUT_FILE_HANDLE>;
    if ($currentLine =~ m/\s*(\d+)\s+electrostatic/)
    {
      ($numberOfESbeads) = $currentLine =~ m/\s*(\d+)\s+electrostatic/;
      $ESFound = 1;
    }
  }
  for (my $esBeadIter = 0 ; $esBeadIter < $numberOfESbeads; $esBeadIter++)
  {
    $currentLine = <MD_INPUT_FILE_HANDLE>;
    my ($beadIndex,$beadCharge) = $currentLine =~ m/^\s*\d+\s+(\d+)\s+(\-?\d+\.+\d+)/;
    $esBeadsRef->[$esBeadIter]->{"ID"} = $beadIndex-1;
    $esBeadsRef->[$esBeadIter]->{"CHARGE"} = $beadCharge;    
  }
  close (MD_INPUT_FILE_HANDLE);

  my $contactIter = -1;
  for (my $firstEsBeadIter = 0 ; $firstEsBeadIter < $numberOfESbeads; $firstEsBeadIter++)
  {
    for (my $secondEsBeadIter = $firstEsBeadIter+1 ; $secondEsBeadIter < $numberOfESbeads; $secondEsBeadIter++)
    {
      if ( ($esBeadsRef->[$firstEsBeadIter]->{"CHARGE"} * $esBeadsRef->[$secondEsBeadIter]->{"CHARGE"}) < 0)
      {
        $contactIter++;
        my $firstBeadID = $esBeadsRef->[$firstEsBeadIter]->{"ID"};
        $contactsRef->[$contactIter]->{"ID1"} = $firstBeadID;
        $contactsRef->[$contactIter]->{"TYPE1"} = $beadsDataRef->[$firstBeadID]->{"TYPE"};    
        my $secondBeadID = $esBeadsRef->[$secondEsBeadIter]->{"ID"};
        $contactsRef->[$contactIter]->{"ID2"} = $secondBeadID;
        $contactsRef->[$contactIter]->{"TYPE2"} = $beadsDataRef->[$secondBeadID]->{"TYPE"};
        $contactsRef->[$contactIter]->{"DISTANCE"} = $executionPreferencesRef->{"ES_CONTACT_DISTANCE"};    

      }
    }
  }
  close (MD_INPUT_FILE_HANDLE);

  return $contactsRef;
}

1;