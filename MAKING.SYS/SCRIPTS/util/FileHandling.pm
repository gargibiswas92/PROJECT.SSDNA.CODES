package FileHandling;
use strict;
use warnings;
use Math::Trig;
use base 'Exporter';
chomp(my $id = `whoami`);
my $randomSeed = time() % 100000;
my $tempFolder = "";

our @EXPORT_OK = qw(unzipIfNeeded zipIfNeeded toScratch fromScratch getTempFolder deleteTempFolder);

sub getTempFolder
{
  if ($tempFolder eq "")
  {
    my $randomNumber = int rand $randomSeed;
    $tempFolder = "/scratch/".$id."/".$randomNumber."/";
    mkdir $tempFolder,0755 or die "Cannot create temporary folder $tempFolder\n";
  }

  return $tempFolder;
}

sub deleteTempFolder
{
  my ($tempFolder) = @_;
  `rm -rf $tempFolder`;
}


sub toScratch
{
  my ($filePath) = @_;
  my $tempFolder = getTempFolder();
  my ($originalFolder,$fileName) = $filePath =~ m/(.*\/)?(.*)/g;
  if (!defined($originalFolder))
  {
    $originalFolder = "./";
    $filePath = "./".$filePath;
  }

  if (!-f $filePath)
  {
    my $bz2FilePath = $filePath.".bz2";
    if (-f $bz2FilePath)
    {
      `cp $bz2FilePath $tempFolder`;
      my $tempFilePath =  $tempFolder.$fileName.".bz2";
      `bunzip2 $tempFilePath`;
    }
    else
    {
      die "Could not find file $filePath or zipped file $bz2FilePath\n";
    }
  }
  else
  {
    `cp $filePath $tempFolder`;
    my $tempFilePath =  $tempFolder.$fileName;
    if ($fileName =~ m /\.bz2$/)
    {
      `bunzip2 $tempFilePath`;      
      $fileName =~ s/.bz2//;
    }
  }
  return ($tempFolder,$fileName,$originalFolder);
}

sub fromScratch
{
  my ($tempFolder,$fileName,$deleteTempFolder,$originalFolder) = @_;
  chomp ($originalFolder = `pwd`) if (!defined($originalFolder));
  $deleteTempFolder = 0 if (!defined($deleteTempFolder));
  my $tempFilePath = $tempFolder.$fileName;
  `cp $tempFilePath $originalFolder`;
  `rm -rf $tempFolder` if $deleteTempFolder;
}

sub unzipIfNeeded
{
  my ($filePath) = @_;

  if ($filePath =~ m/\.bz2$/)
  {
    `bunzip2l $filePath`;
     $filePath =~ s/\.bz2//;
     return (1,$filePath);
  }
  
    return (0,$filePath);
}

sub zipIfNeeded
{
  my ($filePath,$needed) = @_;

  if ($needed)
  {
    `bzip2l $filePath`;
  }
}

1;