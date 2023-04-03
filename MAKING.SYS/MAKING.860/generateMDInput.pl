#!/usr/bin/perl -w
use strict;
use warnings;
$| = 1;

### Directory from which 'generateMDInput.pl' is running
my ($scriptDirectory) = $0 =~ m/(.*\/)/g;

#my $utilDirectory = "/home_d/netaly/scripts/util/";
#my $utilDirectory = "/home_d/arumay/scripts/util/";
#my $utilDirectory = "/home_c/lavi/scripts/util/";
my $utilDirectory = "/home_c/garima/scripts/util/";

### "generateMDInput.prefs.schema" should be inside Script Directory
my $schemaFile = $scriptDirectory."/generateMDInput.prefs.schema";

### "generateMDInput.prefs" should be inside PWD
my $executionPreferencesFile = "./generateMDInput.prefs";

### "/home_d/arumay/scripts/util/" should have "Preferences.pm"
require $utilDirectory."Preferences.pm";

# Accessing the main subroutine of 'Preferences.pm' for getting preferences
my $executionPreferencesRef = Preferences::getPreferences($schemaFile,$executionPreferencesFile);

#require $utilDirectory."PDBHandling.pm";
require $utilDirectory."PDBHandling.pm";
 
#if has static atoms then file must be normalized !!!\
if($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/){
    my $origPdb = $executionPreferencesRef->{"PDB_FILE"};
    `/home_d/netaly/scripts/normalizePDBResidueIndexes.pl $origPdb > $origPdb.norm`;
#    `/home_c/lavi/scripts/normalizePDBResidueIndexes.pl $origPdb > $origPdb.norm`;
    $executionPreferencesRef->{"PDB_FILE"} = "$origPdb.norm";
}


my $pdbDataRef = PDBHandling::readPDBFile($executionPreferencesRef->{"PDB_FILE"});


#store chain order
my $i=0;
my @ChainOrder;
for (my $i = 0; $i < scalar(@{$pdbDataRef->{"CHAINS"}}); $i++)
  {
    $ChainOrder[$i] = $pdbDataRef->{"CHAINS"}->[$i]->{"CHAIN_ID"};
   }
   
require $scriptDirectory."/Models/".$executionPreferencesRef->{"MODEL_TYPE"}.".pm";
$executionPreferencesRef->{"SCRIPT_DIRECTORY"} = $scriptDirectory;
no strict 'refs';
&{$executionPreferencesRef->{"MODEL_TYPE"}."::createMDInputFile"}($executionPreferencesRef,$pdbDataRef,@ChainOrder);
