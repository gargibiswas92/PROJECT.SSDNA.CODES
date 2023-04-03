#!/usr/bin/perl -w
use strict;
use warnings;

my ($scriptDirectory) = $0 =~ m/(.*\/)/g;
my $utilDirectory = "/home_a/golbin/scripts/util/";

my $schemaFile = $scriptDirectory."/generateMDInput.prefs.schema";
my $executionPreferencesFile = "./generateMDInput.prefs";
require $utilDirectory."Preferences.pm";
my $executionPreferencesRef = Preferences::getPreferences($schemaFile,$executionPreferencesFile);
#require $utilDirectory."PDBHandling.pm";
require $utilDirectory."PDBHandling.pm";

#if has static atoms then file must be normalized !!!\
if($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ /YES/){
    my $origPdb = $executionPreferencesRef->{"PDB_FILE"};
    `~golbin/scripts/normalizePDBResidueIndexes.pl $origPdb > $origPdb.norm`;
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
