package Preferences;
use strict;
use warnings;
use base 'Exporter';
our @EXPORT_OK = qw(getPreferences);

my $validationFunctions = {};

# Validation functions

sub checkOptionalOrDefault
{
    my ($preferenceName,$userInput,$formatRef) = @_;
    my $userInputIsValid = 0;
    if (exists($formatRef->{"OPTIONAL"}))
    {
        
        $userInputIsValid = 1 if (($userInput eq "") or (!defined($userInput)));
    }
    if (exists($formatRef->{"DEFAULT"}))
    {
        if (($userInput eq "") or (!defined($userInput)))
        {
          warn "Using default value ".$formatRef->{"DEFAULT"}." for $preferenceName\n";
          $userInput = $formatRef->{"DEFAULT"};
          $userInputIsValid = 1;
        }
    }

    return ($userInputIsValid,$userInput);
}

sub isValidArray
{
    my ($preferenceName,$userInput,$formatRef) = @_;
    my $arrayType = "STRING";
    $arrayType = $formatRef->{"TYPE"} if exists($formatRef->{"TYPE"});
    
    my $deserializeArray = 0;
    $deserializeArray = 1 if (exists($formatRef->{"DESERIALIZE"}) and ($formatRef->{"DESERIALIZE"} =~ m/ARRAY/));
    my $userInputIsValid = 1;
    my @userInputSegements = split(/\,/,$userInput);
    my $deserializedUserInput;
    $deserializedUserInput = [] if $deserializeArray;
    
    foreach my $userInputSegment (@userInputSegements)
    {
      my ($isValid, $userInputSegment) = &{$validationFunctions->{$arrayType}}($preferenceName,$userInputSegment,$formatRef);
      return (0,$userInput) if ($isValid == 0);
      push(@$deserializedUserInput,$userInputSegment) if $deserializeArray;
    }
    return ( $deserializeArray ? ($userInputIsValid,$deserializedUserInput) : ($userInputIsValid,$userInput));
}

sub isValidBinary
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = -B $userInput;
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidDirectory
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  if (exists($formatRef->{"NEW"}))
  {
    mkdir $userInput,0755 unless -d $userInput;
  }
  
  my ($userInputIsValid) = -d $userInput;
  if (exists($formatRef->{"WRITABLE"}))
  {
    $userInputIsValid = ($userInputIsValid and -w $userInput);
  }
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}


sub isValidFlag
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  my $flagOptions = $formatRef->{"OPTIONS"};
  my ($userInputIsValid) = $userInput =~ m/$flagOptions/i;
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,uc($userInput));
}

sub isValidFloat
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = ($userInput =~ m/^-?\d+(.\d+)?(E[+-]\d+)?$/i);
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidFloatRange
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = ($userInput =~ m/(-?\d+(.\d+)?(E[+-]\d+))?(\<\>-?\d+(.\d+)?(E[+-]\d+)?)?/i);
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  if (exists($formatRef->{"DESERIALIZE"}) and ($formatRef->{"DESERIALIZE"} =~ m/RANGE/))
  {
    my @range = split (/\<\>/,$userInput);
    return ($userInputIsValid,\@range);
  }
  else
  {
    return ($userInputIsValid,$userInput);
  }
}

sub isValidInteger
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = 0;
  #$userInputIsValid = ($userInput =~ m/^(-?\d+)$/);
  #is a number && is an integer
  $userInputIsValid = ($userInput =~ m/^-?\d+(.\d+)?(E[+-]\d+)?$/i) && ($userInput - int($userInput) ==0);
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidIntegerRange
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = ($userInput =~ m/^-?\d+(.\d+)?(E[+-]\d+)?(\<\>-?\d+(.\d+)?(E[+-]\d+)?)?/i);
  my @range;
  if ($userInputIsValid)
  {
    @range = split (/\<\>/,$userInput);
    if (scalar(@range) == 2)
    {
      $userInputIsValid = (($range[0] - int($range[0]) ==0) && ($range[1] - int($range[1]) ==0));
    }
    if (scalar(@range) == 1)
    {
      $userInputIsValid = ($range[0] - int($range[0]) ==0);
    }
    
  }
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  if (exists($formatRef->{"DESERIALIZE"}) and ($formatRef->{"DESERIALIZE"} =~ m/RANGE/))
  {
    return ($userInputIsValid,\@range);
  }
  else
  {
    return ($userInputIsValid,$userInput);
  }
}


sub isValidString
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = (($userInput ne "") or (!defined($userInput)));
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidText
{
  my ($preferenceName,$userInput,$formatRef) = @_;
  if (exists($formatRef->{"NEW"}))
  {
    my ($newTextDirectory) = $userInput =~ m/(.*\/)/g;
    $newTextDirectory = $ENV{"PWD"} if (!defined($newTextDirectory));
    my ($userInputIsValid) = -w $newTextDirectory;
    ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
    return ($userInputIsValid,$userInput);
    
  }
  else
  {
    my ($userInputIsValid) = (-T $userInput) or (-T $userInput.".bz2");
    ($userInputIsValid,$userInput) = checkOptionalOrDefault($preferenceName,$userInput,$formatRef) if (!$userInputIsValid);
    return ($userInputIsValid,$userInput);
  }
}


sub loadValidationFunctions
{
    $validationFunctions->{"ARRAY"}=\&isValidArray;
    $validationFunctions->{"BINARY"}=\&isValidBinary;
    $validationFunctions->{"DIRECTORY"}=\&isValidDirectory;
    $validationFunctions->{"FLAG"}=\&isValidFlag;
    $validationFunctions->{"FLOAT"}=\&isValidFloat;
    $validationFunctions->{"FLOAT_RANGE"}=\&isValidFloatRange;
    $validationFunctions->{"INTEGER"}=\&isValidInteger;
    $validationFunctions->{"INTEGER_RANGE"}=\&isValidIntegerRange;
    $validationFunctions->{"STRING"}=\&isValidString;
    $validationFunctions->{"TEXT"}=\&isValidText;    
}


sub loadPreferenceSchemaLine
{
  my $preferencesIndexedMetadataRef = $_[0];
  my $preferencesMetadataRef = $_[1];
  my $currentSchemaLine = $_[2];

  my @currentSchemaLineSegments = split(/\=/,$currentSchemaLine);
  my $currentPreferenceName = shift(@currentSchemaLineSegments);
  my $currentPreferenceMetadata = $currentSchemaLineSegments[0];
  my @currentSchemaFormatSegments = split(/\s+/,$currentPreferenceMetadata);
  my $currentPreferenceType = shift(@currentSchemaFormatSegments);
  my $currentPreferenceFormatRef = {};
  foreach my $formatItem (@currentSchemaFormatSegments)
  {
    my ($formatKey,$formatValue) = split(/\:/,$formatItem);
    $currentPreferenceFormatRef->{$formatKey} = $formatValue;
  }
  my $preferenceIndexesMetadataRef = { "NAME" => uc($currentPreferenceName),
                                         "TYPE" => uc($currentPreferenceType),
                                         "FORMAT" => $currentPreferenceFormatRef,
                                         "METADATA" => $currentPreferenceMetadata};
  push(@$preferencesIndexedMetadataRef,$preferenceIndexesMetadataRef);
  my $preferenceMetadataRef = { "TYPE" => uc($currentPreferenceType),
                                         "FORMAT" => $currentPreferenceFormatRef,
                                         "METADATA" => $currentPreferenceMetadata };
  $preferencesMetadataRef->{uc($currentPreferenceName)} = $preferenceMetadataRef;
}

sub loadPreferenceFlagDependency
{
  my $preferencesFlagDependenciesRef = $_[0];
  my $currentDependencyLine = $_[1];

  my @currentDependencyLineSegments = split(/\:/,$currentDependencyLine);
  my @currentDependenciesSegments = split(/\,/,$currentDependencyLineSegments[1]);
  foreach my $currentDependencySegment (@currentDependenciesSegments)
  {
    push(@{$preferencesFlagDependenciesRef->{$currentDependencySegment}},$currentDependencyLineSegments[0]);
  }
}

sub loadPreferencesSchema
{
  my $preferencesSchemaFile  = $_[0];
  
  open(PREFS_SCHEMA_FILE_HANDLE,$preferencesSchemaFile)
  or die "Could not open schema file $preferencesSchemaFile.\n";
    
  my $currentSchemaLine;
  my $readPreferencesState = -1;  
  
  my $preferencesIndexedMetadataRef = [];
  my $preferencesMetadataRef = {};
  my $preferencesFlagDependenciesRef = {};
  
  while ($currentSchemaLine = <PREFS_SCHEMA_FILE_HANDLE>)
  {
    if ($readPreferencesState == -1)
    {
      if ($currentSchemaLine =~ m/\[PREFERENCES\]/i)
      {
        $readPreferencesState = 0;
      }
      next;
    }

    if ($readPreferencesState == 0)
    {
      if ($currentSchemaLine =~ m/\[DEPENDENCIES\]/i)
      {
        $readPreferencesState = 1;
        next;
      }
      else
      {
        chomp($currentSchemaLine);
        loadPreferenceSchemaLine($preferencesIndexedMetadataRef,
                                          $preferencesMetadataRef,
                                          $currentSchemaLine);
      }
    }

    if ($readPreferencesState == 1)
    {
      chomp($currentSchemaLine);
      loadPreferenceFlagDependency($preferencesFlagDependenciesRef,$currentSchemaLine);
    }

  }

  close(PREFS_SCHEMA_FILE_HANDLE);

  return ($preferencesIndexedMetadataRef,$preferencesMetadataRef,$preferencesFlagDependenciesRef);
} 

sub savePreferences
{
  my ($preferencesFile, $preferencesSchemaFile,$preferencesRef,$preferencesIndexedMetadataRef) = @_;
  
  print "Storing preferences in $preferencesFile.\n";

  open (PREFERENCES_FILE_HANDLE,">",$preferencesFile);
  for (my $preferenceIter = 0; $preferenceIter < scalar(@$preferencesIndexedMetadataRef); $preferenceIter++)
  {
    my $preferenceName = $preferencesIndexedMetadataRef->[$preferenceIter]->{"NAME"};
    if ((defined($preferencesRef->{$preferenceName})) and ($preferencesRef->{$preferenceName} ne ""))
    {
      print PREFERENCES_FILE_HANDLE $preferenceName."=".$preferencesRef->{$preferenceName}."\n";
    }
  }
  close(PREFERENCES_FILE_HANDLE);
}

# main subroutine for getting preferences
sub getPreferences
{
  my $preferencesRef = {};
  
  my $preferencesSchemaFile = $_[0];
  die "Preferences schema file $preferencesSchemaFile was not found.\n" unless -T $preferencesSchemaFile;

  my ($preferencesIndexedMetadataRef,$preferencesMetadataRef,$preferencesFlagDependenciesRef) = loadPreferencesSchema($preferencesSchemaFile);
 
  my $preferencesFile = $_[1];
  if (!-T $preferencesFile)
  {
    #Preference file was not found
    print "Preferences file was not found\n";
    $preferencesRef = getUserPreferences($preferencesIndexedMetadataRef,$preferencesMetadataRef,$preferencesFlagDependenciesRef);
    
    my $isValidAnswer = 0;
    while ($isValidAnswer == 0)
    {
      
      print "Would you like to store the preferences? (YES|NO):";
      chomp (my $userInput = <STDIN>);
      $userInput = "YES" if ($userInput eq "");
      if ($userInput !~ m/YES|NO/i)
      {
        warn "Invalid answer $userInput.\n";
      }
      else
      {
        $isValidAnswer = 1;
        if ($userInput =~ m/YES/i)
        {
          my ($preferencesFile) = $0 =~ m/\/([\w|\.]+)\./g;
          $preferencesFile .= ".prefs";
          print "Enter name of new preferences file (default is $preferencesFile):";
          chomp (my $userInputPreferenceFile = <STDIN>);
          $userInputPreferenceFile = $preferencesFile if ($userInputPreferenceFile eq "");
          savePreferences($preferencesFile, $preferencesSchemaFile,$preferencesRef,$preferencesIndexedMetadataRef);
        }
      }
    }
  }
  else
  {
    # preference file was found
    print "Preferences file was found\n";
    $preferencesRef = loadPreferences($preferencesMetadataRef,$preferencesIndexedMetadataRef,$preferencesFlagDependenciesRef,$preferencesFile);
  }
  
  return $preferencesRef;
}

sub validatePreferenceDependencies
{
  my $preferenceName = $_[0];
  my $preferencesRef = $_[1];
  my $preferencesFlagDependenciesRef = $_[2];
  my $preferencesMetadataRef = $_[3];
  
  my $allRequiredFlagsAreSet = 1;
  
  my $requiredFlagsRef = $preferencesFlagDependenciesRef->{$preferenceName};
  foreach $requiredFlagsRef (@$requiredFlagsRef)
  {
    my @requiredFlagSegments = split(/\=/,$requiredFlagsRef);
    updateDefaultValue($preferencesRef->{$requiredFlagSegments[0]},$preferencesRef,$preferencesMetadataRef)
    if (!exists($preferencesRef->{$requiredFlagSegments[0]}));
    
    if (exists($preferencesRef->{$requiredFlagSegments[0]}))
    {
      my @possibleFlagValues = split(/\,/,$requiredFlagSegments[1]);
      my $flagIsRaised = 0;
      foreach my $possibleFlagValue (@possibleFlagValues)
      {
        $flagIsRaised++ if ($preferencesRef->{$requiredFlagSegments[0]} eq $possibleFlagValue);
      }
      $allRequiredFlagsAreSet = 0 if ($flagIsRaised == 0);
    }
    else
    {
      $allRequiredFlagsAreSet = 0;
    }
  }
  
  return $allRequiredFlagsAreSet;
}

sub getUserPreference
{
    my $preferencesRef = $_[0];
    my $preferenceMetadataIter = $_[1];
    my $preferencesIndexedMetadataRef = $_[2];
    my $preferenceName = $preferencesIndexedMetadataRef->[$preferenceMetadataIter]->{"NAME"};
    my $preferenceMetadata = $preferencesIndexedMetadataRef->[$preferenceMetadataIter]->{"METADATA"};
    my $preferenceType = $preferencesIndexedMetadataRef->[$preferenceMetadataIter]->{"TYPE"};
    my $preferenceFormat = $preferencesIndexedMetadataRef->[$preferenceMetadataIter]->{"FORMAT"};
    my $preferencesFlagDependenciesRef = $_[3];
    my $preferencesMetadataRef = $_[4];
    return if (!validatePreferenceDependencies($preferenceName,$preferencesRef,$preferencesFlagDependenciesRef,$preferencesMetadataRef)); 
    
    my ($userInput,$processedUserInput);
    my $isValidAnswer = 0;
    while ($isValidAnswer == 0)
    {
      
      print "Enter $preferenceName ($preferenceMetadata):";
      chomp ($userInput = <STDIN>);
      ($isValidAnswer,$processedUserInput) = &{$validationFunctions->{$preferenceType}}($preferenceName,$userInput,$preferenceFormat);
      if (!$isValidAnswer)
      {
        warn "Invalid value $userInput for $preferenceName, Should be $preferenceType.\n";
      }
      else
      {
        $isValidAnswer = 1;
      }
    }
    
    if (defined($processedUserInput) and ($processedUserInput ne ""))
    {
      $preferencesRef->{$preferenceName} = $processedUserInput;
    }
}

sub getUserPreferences
{
    my $preferencesRef = {};
    
    my $preferencesIndexedMetadataRef = $_[0];
    my $preferencesMetadataRef = $_[1];
    my $preferencesFlagDependenciesRef = $_[2];

    
    for ( my $preferenceMetadataIter = 0; $preferenceMetadataIter < scalar(@$preferencesIndexedMetadataRef); $preferenceMetadataIter++)
    {
      my $preferenceName = $preferencesIndexedMetadataRef->[$preferenceMetadataIter]->{"NAME"};
      if (validatePreferenceDependencies($preferenceName,$preferencesRef,$preferencesFlagDependenciesRef,$preferencesMetadataRef))
      {
        getUserPreference($preferencesRef,
                          $preferenceMetadataIter,
                          $preferencesIndexedMetadataRef,
                          $preferencesFlagDependenciesRef,
                          $preferencesMetadataRef);
      }
    }
    return $preferencesRef;
}


sub readPreferenceLine
{
    my $preferencesRef = $_[0];
    my $preferenceLine = $_[1];
    my $preferencesMetadataRef = $_[2];
    my $preferencesFlagDependenciesRef = $_[3];

    chomp($preferenceLine);
    return if ($preferenceLine =~ m/^\#/);
    my @preferenceLineSegments = split(/\=/,$preferenceLine);
    my $preferenceName = shift(@preferenceLineSegments);
    my $preferenceValue = shift(@preferenceLineSegments);
    warn "Invalid preference $preferenceName found in preferences file.\n"
    unless (exists($preferencesMetadataRef->{$preferenceName}));
    my $preferenceType = $preferencesMetadataRef->{$preferenceName}->{"TYPE"};
    my $preferenceFormat = $preferencesMetadataRef->{$preferenceName}->{"FORMAT"};
    return if (!validatePreferenceDependencies($preferenceName,$preferencesRef,$preferencesFlagDependenciesRef,$preferencesMetadataRef)); 
    
    my ($userInputisValid,$processedPreferenceValue) = &{$validationFunctions->{$preferenceType}}($preferenceName,$preferenceValue,$preferenceFormat);
    die "Invalid $preferenceType for $preferenceName: $preferenceValue.\n" unless ($userInputisValid);
    
    $preferencesRef->{$preferenceName} = $processedPreferenceValue;
}

sub updateDefaultValue
{
    my ($preferenceName,$preferencesRef,$preferencesMetadataRef) = @_;
    return if (!defined($preferenceName) or $preferenceName eq "");
    if (defined($preferencesMetadataRef->{$preferenceName}->{"FORMAT"}->{"DEFAULT"}))
    {
      warn "Using default value ".$preferencesMetadataRef->{$preferenceName}->{"FORMAT"}->{"DEFAULT"}." for $preferenceName\n";
      $preferencesRef->{$preferenceName} = $preferencesMetadataRef->{$preferenceName}->{"FORMAT"}->{"DEFAULT"};
    }
    
}

sub updateDefaultPreferences
{
  my ($preferencesRef,$preferencesMetadataRef,$preferencesIndexedMetadataRef,$preferencesFlagDependenciesRef) = @_;
  for (my $preferenceIter = 0; $preferenceIter < scalar(@$preferencesIndexedMetadataRef); $preferenceIter++)
  {
    my $preferenceName = $preferencesIndexedMetadataRef->[$preferenceIter]->{"NAME"};
    if (!defined($preferencesRef->{$preferenceName}) and
        validatePreferenceDependencies($preferenceName,$preferencesRef,$preferencesFlagDependenciesRef,$preferencesMetadataRef))
    {
      updateDefaultValue($preferenceName,$preferencesRef,$preferencesMetadataRef)
    }
  }
}

sub loadPreferences
{
  my ($preferencesMetadataRef,
      $preferencesIndexedMetadataRef,
      $preferencesFlagDependenciesRef,
      $preferencesFile) = @_;

  my $preferencesRef = {};
    
  open(PREFERENCES_FILE_HANDLE,$preferencesFile)
  or die "Could not open preferences file $preferencesFile.\n";
    
  my $currentPreferencesLine;
  while ($currentPreferencesLine = <PREFERENCES_FILE_HANDLE>)
  {
    readPreferenceLine($preferencesRef,
                                $currentPreferencesLine,
                                $preferencesMetadataRef,
                                $preferencesFlagDependenciesRef);
  }
  updateDefaultPreferences($preferencesRef,$preferencesMetadataRef,$preferencesIndexedMetadataRef,$preferencesFlagDependenciesRef);

  close(PREFERENCES_FILE_HANDLE);

 return $preferencesRef;
  
}


loadValidationFunctions();
1;
