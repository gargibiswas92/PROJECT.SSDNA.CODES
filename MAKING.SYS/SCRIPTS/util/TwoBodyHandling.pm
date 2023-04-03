package TwoBodyHandling;
use strict;
use warnings;
use Math::Trig;

use base 'Exporter';
our @EXPORT_OK = qw(findFramesWithinQRange);

sub findFramesWithinQRange
{
  my ($twoBodyFileName, $minimalQ, $maximalQ) = @_;
  my $frameRangesRef = [];
  my $validFramesFound = 0;
  
  open (TWO_BODY_FILE_HANDLE, $twoBodyFileName) or die "Error reading from ".$twoBodyFileName.": $!\n";

  my $currentTimeStep = 0;
  my $currentLine;
  my $currentRangeIter = -1;
  my $currentRangeStatus = 0;
  
  while ($currentLine = <TWO_BODY_FILE_HANDLE>)
  {
    $currentTimeStep++;
    my ($currentQ) = $currentLine =~ m/(\d+)/;
    if ($currentRangeStatus == 0)
    {
      if (($currentQ >= $minimalQ) and ($currentQ <= $maximalQ))
      {
        $validFramesFound++;
        $currentRangeIter++;
        $frameRangesRef->[$currentRangeIter]->{"START"} = $currentTimeStep;
        $currentRangeStatus = 1;
        next;
      }
    }
    if ($currentRangeStatus == 1)
    {
      if (($currentQ < $minimalQ) or ($currentQ > $maximalQ))
      {
        $frameRangesRef->[$currentRangeIter]->{"END"} = $currentTimeStep-1;
        $currentRangeStatus = 0;
      }
      else
      {
        $validFramesFound++;
      }
    }
  }
  
  if ($currentRangeStatus == 1)
  {
    $frameRangesRef->[$currentRangeIter]->{"END"} = $currentTimeStep;
  }
  return ($validFramesFound,$frameRangesRef);
}

sub frameIsInRange
{
  my ($frameRangesRef,$currentRangeIter,$frameIter) = @_;
  return (0,$currentRangeIter) if ($currentRangeIter >= scalar(@$frameRangesRef));
  if ($frameIter > $frameRangesRef->[$currentRangeIter]->{"END"})
  {
    $currentRangeIter++;
    return (0,$currentRangeIter) if ($currentRangeIter >= scalar(@$frameRangesRef));
  }
  
  my $frameIsInRange = ( ($frameIter >= $frameRangesRef->[$currentRangeIter]->{"START"}) and
                         ($frameIter <= $frameRangesRef->[$currentRangeIter]->{"END"}) );
  return ($frameIsInRange,$currentRangeIter);
}

1;