#!/usr/bin/perl -w

# MD wrapper perl script
# based on run_GO scripts written by Amit Mor
# created by: Ariel Azia
# first creation date: 14/5/2007
#2.1 corrections for DNA
#2.2 adding support for debye huckel feature of MD
#2.3 adding chirality
$| = 1;

use strict;
use warnings;
use Term::ANSIColor qw(:constants);
my $minimumEsDistanceByModelType = { "CA_NO" => 4,					
                                     "CA_YES" => 3,
                                     "CB_NO" => 3,
                                     "CB_YES" => 2
                                   };

my ($scriptDirectory) = $0 =~ m/(.*\/)/g;                 # gets the directory of the execution
my $utilDirectory = "/home/garima/Scripts/util/";        # path to ariel util directory 
my $mdWrapperSchemaFileName = $scriptDirectory."/MDWrapper.prefs.schema";    # full path to ariel schema file 
my $mdWrapperPrefsFileName = "./MDWrapper.prefs";
require $utilDirectory."Preferences.pm";

my $executionPreferencesRef = Preferences::getPreferences($mdWrapperSchemaFileName,$mdWrapperPrefsFileName);
processExecutionPreferences($executionPreferencesRef);
createOutputFolders();
launchExecutions();

sub processExecutionPreferences
{
    my $executionPreferencesRef = $_[0];
    chomp(my $id = `whoami`);
    $executionPreferencesRef->{"ID"} = $id;
    $executionPreferencesRef->{"RANDOM_SEED"} = time() % 100000;

  ($executionPreferencesRef->{"INPUT_FILE_NAME"}) = $executionPreferencesRef->{"INPUT_FILE_PATH"} =~ m/.*\/(.*)/;
  ($executionPreferencesRef->{"INPUT_FILE_PREFIX"}) = $executionPreferencesRef->{"INPUT_FILE_NAME"} =~ m/(.*)\./;
  if ($executionPreferencesRef->{"USE_INITIAL_CONDITIONS"} != 3)
  {
    ($executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_NAME"}) = $executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_PATH"} =~ m/.*\/(.*)/;
  }
#  $executionPreferencesRef->{"TEMPERATURES_BUFFER"} = $executionPreferencesRef->{"TEMPERATURES"};
#  my @temperatures = split(/\,/,$executionPreferencesRef->{"TEMPERATURES_BUFFER"});
#  $executionPreferencesRef->{"TEMPERATURES"} = \@temperatures;
 
  $executionPreferencesRef->{"FORCES_BUFFER"} = $executionPreferencesRef->{"FORCES"};
  my @forces = split(/\,/,$executionPreferencesRef->{"FORCES_BUFFER"});
  $executionPreferencesRef->{"FORCES"} = \@forces;
 
  use File::Spec;
  $executionPreferencesRef->{"BASE_PATH"} = File::Spec->rel2abs( $executionPreferencesRef->{"BASE_PATH"} );
}

sub createOutputFolders
{
    my $path;
    
    if ($executionPreferencesRef->{"EXECUTION_LOCATION"} =~ m/Cluster/i)
    {
	if (($executionPreferencesRef->{"OUTPUT_PATH"} eq "/cluster_data/"))
	{
	    $path = $executionPreferencesRef->{"OUTPUT_PATH"}."/".$executionPreferencesRef->{"ID"}."/".$executionPreferencesRef->{"RANDOM_SEED"};
	}
	else
	{
	    $path = $executionPreferencesRef->{"OUTPUT_PATH"};
	}
    }
    else
    {
       if ((!defined($executionPreferencesRef->{"EXECUTION_DIRECTORY"})) or ($executionPreferencesRef->{"EXECUTION_DIRECTORY"} eq ""))
       {
            $path = $executionPreferencesRef->{"BASE_PATH"}."/".$executionPreferencesRef->{"INPUT_FILE_PREFIX"}."_".$executionPreferencesRef->{"RANDOM_SEED"};
            $path .= "_hsa" if ($executionPreferencesRef->{"HAS_STATIC_ATOMS"});
            if ($executionPreferencesRef->{"CONFINE_IN_BOX"} =~ m/YES/i)
            {
	       my $boxXdim = $executionPreferencesRef->{"BOX_X_MAX"} - $executionPreferencesRef->{"BOX_X_MIN"};
	       my $boxYdim = $executionPreferencesRef->{"BOX_Y_MAX"} - $executionPreferencesRef->{"BOX_Y_MIN"};
	       my $boxZdim = $executionPreferencesRef->{"BOX_Z_MAX"} - $executionPreferencesRef->{"BOX_Z_MIN"};
               $path .= "_cib_".$boxXdim."_".$boxYdim."_".$boxZdim;
	    }
	    if ($executionPreferencesRef->{"USE_ELECTROSTATICS"} =~ m/YES/i)
	    {
	      $path .= "_es_".$executionPreferencesRef->{"DIELECTRIC_CONSTANT"};
	      if ($executionPreferencesRef->{"USE_DEBYE_HUCKEL"} =~ m/YES/i)
	      {
              $path .= "_dh_".$executionPreferencesRef->{"IONIC_STRENGTH"}."_";
              $path .= $executionPreferencesRef->{"IONIC_RADIUS"}."_";
              $path .= $executionPreferencesRef->{"SOLVENT_DENSITY"};
	      }
	    }
	}
	else
	{
	    $path = $executionPreferencesRef->{"BASE_PATH"}."/".$executionPreferencesRef->{"EXECUTION_DIRECTORY"};
	}
    }

    $executionPreferencesRef->{"PATH"} = $path;
    
    mkdir($path,0755) unless -d $path;
    -d $path or die "Could not create execution directory\n";

    mkdir($path.'/input',0755) unless -d $path.'/input';
    -d $path.'/input' or die "Could not create input directory\n";
    mkdir($path.'/output',0755) unless -d $path.'/output';
    -d $path.'/output' or die "Could not create output directory\n";
    mkdir($path.'/output/Temperature',0755) unless -d $path.'/output/Temperature';
    -d $path.'/output/Temperature' or die "Could not create output directory: temperature\n";
    mkdir($path.'/output/Distance',0755) unless -d $path.'/output/Distance';
    -d $path.'/output/Distance' or die "Could not create output directory: Distance\n";
    mkdir($path.'/output/2Body',0755) unless -d $path.'/output/2Body';
    -d $path.'/output/2Body' or die "Could not create output directory: 2Body\n";
    mkdir($path.'/output/3Body',0755) unless -d $path.'/output/3Body';
    -d $path.'/output/3Body' or die "Could not create output directory: 3Body\n";
    mkdir($path.'/output/Etotal',0755) unless -d $path.'/output/Etotal';
    -d $path.'/output/Etotal' or die "Could not create output directory: Etotal\n";
    mkdir($path.'/output/EbyType',0755) unless -d $path.'/output/EbyType';
    -d $path.'/output/EbyType' or die "Could not create output directory:EbyType\n";
    mkdir($path.'/output/Force_Spring',0755) unless -d $path.'/output/Force_Spring';
    -d $path.'/output/Force_Spring' or die "Could not create output directory:Force_Spring\n";
    mkdir($path.'/output/Log',0755) unless -d $path.'/output/Log';
    -d $path.'/output/Log' or die "Could not create output directory:Log\n";
    mkdir($path.'/output/LastCord',0755) unless -d $path.'/output/LastCord';
    -d $path.'/output/LastCord' or die "Could not create output directory:LastCord\n";
    if ($executionPreferencesRef->{"STORE_TRAJECTORIES"} =~ m/YES/i)
    {
      mkdir($path.'/output/Traj',0755) unless -d $path.'/output/Traj';
      -d $path.'/output/Traj' or die "Could not create output directory:Traj\n";
      if ($executionPreferencesRef->{"STORE_TRAJECTORY_DISTANCES"} =~ m/YES/i)
      {
        mkdir($path.'/output/TrajDist',0755) unless -d $path.'/output/TrajDist';
        -d $path.'/output/TrajDist' or die "Could not create output directory:TrajDist\n";
      }
    }

#    if ($executionPreferencesRef->{"WRITE_ALL_CONTACTS"} =~ m/YES/i)
#    {
#      mkdir($path.'/output/AllContacts',0755) unless -d $path.'/output/AllContacts';
#      -d $path.'/output/AllContacts' or die "Could not create output directory:AllContacts\n";
#    }

#    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
#    {
#      mkdir($path.'/output/ContactRanges',0755) unless -d $path.'/output/ContactRanges';
#      -d $path.'/output/ContactRanges' or die "Could not create output directory:ContactRanges\n";
#    }

}

sub createSettingsFile
{
  my ($currentTemperature,$currentForce, $randomNumber) = @_;
  my $currentTemperatureRange;
  if ($executionPreferencesRef->{"TEMPERATURE_ANNEALING"} =~ m/YES/i)
  {
    $currentTemperatureRange = $currentTemperature;
    $currentTemperature = $currentTemperatureRange->[0];
  }	
	
  my $settingsFileName = 'settings_'.$randomNumber.'.dat';
  my $executionUID = $executionPreferencesRef->{"INPUT_FILE_PREFIX"}."_".$currentTemperature."_".$currentForce."_".$randomNumber;
#  print "$executionUID";
  open(SETTINGS_FILE_HANDLE, ">",$executionPreferencesRef->{"PATH"}."/input/$settingsFileName") or die "could not create settings file in input subdirectory ".$executionPreferencesRef->{"PATH"}."\/input\/\n";
  print SETTINGS_FILE_HANDLE "This is a settings file to be used for running MD simulations\n\n";
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_INITIAL_CONDITIONS"}." : whether to use initial conditions\n";
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"INPUT_FILE_NAME"}." : name of input file for execution\n";
  if ($executionPreferencesRef->{"USE_INITIAL_CONDITIONS"} == 3)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"INPUT_FILE_NAME"};
  }
  else
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_NAME"}
  }
  print SETTINGS_FILE_HANDLE ", Final_".$executionUID.".dat";

# print ", Final_".$executionUID.".dat";
#  sleep(5);


  print SETTINGS_FILE_HANDLE " : initial coordinates/velocities, final conformation\n";
  print SETTINGS_FILE_HANDLE $executionUID.".log : log file for current execution\n";
  if ($executionPreferencesRef->{"STORE_TRAJECTORIES"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE "Traj_".$executionUID.".dat, ".$executionPreferencesRef->{"TRAJECTORY_OUTPUT_FREQUENCY"};
    print SETTINGS_FILE_HANDLE " : Trajectory file name, Trajectory output frequency(in integer number of time steps)\n";
    if ($executionPreferencesRef->{"STORE_TRAJECTORY_DISTANCES"} =~ m/YES/i)
    {
      print SETTINGS_FILE_HANDLE "Traj_Dist_".$executionUID.".dat ";
      print SETTINGS_FILE_HANDLE " : Trajectory distances file name\n";
    }
    else
    {
      print SETTINGS_FILE_HANDLE "NO ";
      print SETTINGS_FILE_HANDLE " : Trajectory distances file name\n";
    }
  }
  else
  {
    print SETTINGS_FILE_HANDLE "NO, 1";
    print SETTINGS_FILE_HANDLE " : Trajectory file name, Trajectory output frequency(in integer number of time steps)\n";
  }
  print SETTINGS_FILE_HANDLE "Etotal_".$executionUID.".dat, EbyType_".$executionUID.".dat";
  print SETTINGS_FILE_HANDLE " : Total potential energy, Energy by term\n";
  print SETTINGS_FILE_HANDLE "2body_".$executionUID.".dat, 3body_".$executionUID.".dat";
  print SETTINGS_FILE_HANDLE " : 2 body output file, 3 body output file\n";
  print SETTINGS_FILE_HANDLE "DistanceFile_".$executionUID.".dat";
  print SETTINGS_FILE_HANDLE " : distance output file\n";
  print SETTINGS_FILE_HANDLE "TemperatureFile_".$executionUID.".dat";
  print SETTINGS_FILE_HANDLE " : Temperature\n";
   
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"EXECUTION_TYPE"};
  print SETTINGS_FILE_HANDLE " : type of execution (MD / LD)\n";
  if ($executionPreferencesRef->{"EXECUTION_TYPE"} =~ m/LD/i)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"GAMMA"};
    print SETTINGS_FILE_HANDLE " : gamma drag coefficient for bath coupling\n";
  }
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"EXECUTION_STEPS"}.", ".$executionPreferencesRef->{"OUTPUT_FREQUENCY"};
  print SETTINGS_FILE_HANDLE " : total steps for this execution, writing to output frequency\n";
  printf SETTINGS_FILE_HANDLE "%10.5f%10.5f%s\n",$executionPreferencesRef->{"TAU"}, 0.2," : tau, Rescale By Tau";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"TEMPERATURE_ANNEALING"};
  print SETTINGS_FILE_HANDLE " : whether to apply temperature annealing\n";
  
  if ($executionPreferencesRef->{"TEMPERATURE_ANNEALING"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE $currentTemperatureRange->[0]." : start temperature for this execution\n";
    print SETTINGS_FILE_HANDLE $currentTemperatureRange->[1]." : end temperature for this execution\n";
    
  }
  else
  {
    print SETTINGS_FILE_HANDLE $currentTemperature." : temperature for this execution\n";
  }

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"PULL_SPRINGS"};
  print SETTINGS_FILE_HANDLE " : whether to pull with spring\n";

  if ($executionPreferencesRef->{"PULL_SPRINGS"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{SPRING_COEFF}." : KSpring\n";
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{SPRING_LENGTH}." : RSpring\n";
#    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{PULLING_VELOCITY}." : Vp the rate of pulling the spring\n";
    print SETTINGS_FILE_HANDLE $currentForce." : Vp the rate of pulling the spring\n";
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{PULLING_RESIDUE1}." ".$executionPreferencesRef->{PULLING_RESIDUE2};
    print SETTINGS_FILE_HANDLE "  : the attachment residues2 of the springs\n";
    print SETTINGS_FILE_HANDLE "FstaticSpring_".$executionUID.".dat, FmovingSpring_".$executionUID.".dat \n";
  }
  else
  {
    print SETTINGS_FILE_HANDLE $currentForce." : force for this execution\n";
  }

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"HAS_STATIC_ATOMS"};
  print SETTINGS_FILE_HANDLE " : whether input file has static atoms\n";
  if ($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"LAST_DYNAMIC_ATOM"};
    print SETTINGS_FILE_HANDLE " : last dynamic atom\n";
  }

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"CONFINE_IN_BOX"};
  print SETTINGS_FILE_HANDLE " : whether to confine model in box\n";
  if ($executionPreferencesRef->{"CONFINE_IN_BOX"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"BOX_X_MIN"}.", ".$executionPreferencesRef->{"BOX_Y_MIN"}.", ".$executionPreferencesRef->{"BOX_Z_MIN"};
    print SETTINGS_FILE_HANDLE " : box minimal dimensions\n";
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"BOX_X_MAX"}.", ".$executionPreferencesRef->{"BOX_Y_MAX"}.", ".$executionPreferencesRef->{"BOX_Z_MAX"};
    print SETTINGS_FILE_HANDLE " : box maximal dimensions\n";
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"BOX_COEFFICIENT"};
    print SETTINGS_FILE_HANDLE " : box force coefficient\n";
  }

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_ELECTROSTATICS"};
  print SETTINGS_FILE_HANDLE " : whether to use electrostatics in execution\n";
  if ($executionPreferencesRef->{"USE_ELECTROSTATICS"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"DIELECTRIC_CONSTANT"};
    print SETTINGS_FILE_HANDLE " : dielectric constant\n";
    my $minimumEsDistance = $minimumEsDistanceByModelType->{uc($executionPreferencesRef->{"MODEL_TYPE"})."_".uc($executionPreferencesRef->{"APPLY_ELECTROSTATICS_BETWEEN_NEIGHBORS"})};
    print SETTINGS_FILE_HANDLE $minimumEsDistance;
    print SETTINGS_FILE_HANDLE " : minimum index difference between beads to apply electrostatics\n";
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"ELECTROSTATICS_CUTOFF_TYPE"};
    print SETTINGS_FILE_HANDLE " : type of cutoff for electrostatic interaction\n";

    if ($executionPreferencesRef->{"ELECTROSTATICS_CUTOFF_TYPE"} =~ m/ENERGY/i)
    {
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"ELECTROSTATICS_ENERGY_CUTOFF"};
      print SETTINGS_FILE_HANDLE " : maximum fraction of energy for electrostatic interaction\n";
    }

    if ($executionPreferencesRef->{"ELECTROSTATICS_CUTOFF_TYPE"} =~ m/DISTANCE/i)
    {
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"ELECTROSTATICS_DISTANCE_CUTOFF"};
      print SETTINGS_FILE_HANDLE " : maximum distance for electrostatic interaction\n";
    }

    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_DEBYE_HUCKEL"};
    print SETTINGS_FILE_HANDLE " : whether to use debye huckel in execution\n";
  
    if ($executionPreferencesRef->{"USE_DEBYE_HUCKEL"} =~ m/YES/i)
    {
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"IONIC_STRENGTH"};
      print SETTINGS_FILE_HANDLE " : ionic strength of solution\n";
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"IONIC_RADIUS"};
      print SETTINGS_FILE_HANDLE " : ionic radius of ions in solution (1.4 for NaCl) \n";
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"SOLVENT_DENSITY"};
      print SETTINGS_FILE_HANDLE " : solventDensity (1 for water) \n";
    }

  }
  
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_CHIRALS"};
  print SETTINGS_FILE_HANDLE " : whether to use chirals in execution\n";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_ELLIPSOID_REPULSIONS"};
  print SETTINGS_FILE_HANDLE " : whether to use ellipsoid repulsions in execution\n";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"HAS_SINGLE_STRAND_DNA"};
  print SETTINGS_FILE_HANDLE " : whether to simulate dynamics of the single stranded DNA\n";

  if ($executionPreferencesRef->{"HAS_SINGLE_STRAND_DNA"} =~ m/YES/i)
  {
     print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_CENTERED_DNA_STRAND"};
     print SETTINGS_FILE_HANDLE " : whether to center the single stranded DNA\n";
     
     if ($executionPreferencesRef->{"USE_CENTERED_DNA_STRAND"} =~ m/YES/i)
     {
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"FIRST_DNA_BEAD"};
      print SETTINGS_FILE_HANDLE " : the first DNA bead\n";
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"LAST_DNA_BEAD"};
      print SETTINGS_FILE_HANDLE " : the last DNA bead\n";
     }  
  }

  close(SETTINGS_FILE_HANDLE);

  return $settingsFileName;
}

sub createRandomFile
{
    my ($seed, $randomNumber) = @_;
    my $i = int rand $seed;
    my $j = int rand $seed;
    my $k = int rand $seed;
    my $randomFileName = 'random_' . $randomNumber .'.dat';
    open (RANDOM_FILE_HANDLE, '>',$executionPreferencesRef->{"PATH"}.'/input/'.$randomFileName) or die "could not create random input file ".$executionPreferencesRef->{"PATH"}."/input/".$randomFileName."\n";
    
    print RANDOM_FILE_HANDLE "$i $j $k";
    
    close RANDOM_FILE_HANDLE;
    open (MDW_LOG_FILE_HANDLE, ">>",$executionPreferencesRef->{"LOG_FILE"}) or die "can't open log file :$!\n";
    print MDW_LOG_FILE_HANDLE localtime(time),":Random values: $i $j $k\n";
    close MDW_LOG_FILE_HANDLE;
    return $randomFileName;
}

sub reportExecution
{
    my ($randomNumber) = @_;
    open (MDW_LOG_FILE_HANDLE, ">>",$executionPreferencesRef->{"LOG_FILE"}) or die "can't open log file:$!\n";
    print MDW_LOG_FILE_HANDLE localtime(time),":Number of steps in simulation: ",$executionPreferencesRef->{"EXECUTION_STEPS"},"\n";
    print MDW_LOG_FILE_HANDLE localtime(time),":Temperatures submitted: ",$executionPreferencesRef->{"TEMPERATURES_BUFFER"},"\n";
    print MDW_LOG_FILE_HANDLE localtime(time),":Forces submitted: ",$executionPreferencesRef->{"FORCES_BUFFER"},"\n";
    print MDW_LOG_FILE_HANDLE localtime(time),":Base seed for simulation:", $executionPreferencesRef->{"RANDOM_SEED"},"\n";
    close(MDW_LOG_FILE_HANDLE);
}

sub launchExecutions
{

  
  foreach (1..$executionPreferencesRef->{"ITERATIONS_PER_TEMPERATURE"})
  {
    if ($executionPreferencesRef->{"TEMPERATURE_ANNEALING"} =~ m/YES/i)
    {
      foreach my $currentTemprature (@{$executionPreferencesRef->{"TEMPERATURE_RANGES"}})
      {
        foreach my $currentForce (@{$executionPreferencesRef->{"FORCES"}})
        {
        my $randomNumber = int rand $executionPreferencesRef->{"RANDOM_SEED"};
        my $settingsFileName = createSettingsFile($currentTemprature,$currentForce, $randomNumber);
        my $randomFileName = createRandomFile($executionPreferencesRef->{"RANDOM_SEED"},$randomNumber);
        launchExecution($currentTemprature->[0],$currentForce, $settingsFileName, $randomNumber, $randomFileName);
        reportExecution($randomNumber);
        }
      }
    }
    else
    {
    foreach my $currentTemprature (@{$executionPreferencesRef->{"TEMPERATURES"}})
      {
	 foreach my $currentForce (@{$executionPreferencesRef->{"FORCES"}})
        {
        my $randomNumber = int rand $executionPreferencesRef->{"RANDOM_SEED"};
        my $settingsFileName = createSettingsFile($currentTemprature,$currentForce,$randomNumber);
        my $randomFileName = createRandomFile($executionPreferencesRef->{"RANDOM_SEED"},$randomNumber);
        launchExecution($currentTemprature,$currentForce,$settingsFileName, $randomNumber, $randomFileName);
        reportExecution($randomNumber);
      }
    }
  }
}
}
sub launchExecution
{
    my ($currentTemprature,$currentForce, $settingFileName, $randomNumber, $randomFileName) = @_;
    print "$settingFileName,$randomFileName,$randomFileName\n";
    my $settingsFilePath = $executionPreferencesRef->{"PATH"}.'/input/'.$settingFileName;
    my $randomFilePath = $executionPreferencesRef->{"PATH"}.'/input/'.$randomFileName;
    my $executionUID = $executionPreferencesRef->{"INPUT_FILE_PREFIX"}."_".$currentTemprature."_".$currentForce."_".$randomNumber;
    my $cshScript = '_'.$executionUID.'.csh';
    my $path = $executionPreferencesRef->{"PATH"};
    my $id = $executionPreferencesRef->{"ID"};
    open (EXECUTION_SCRIPT_FILE_HANDLE, ">",$cshScript) or die "Couldn't create the running script: $!\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#!/bin/csh -f\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -cwd\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -j y\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -S /bin/csh\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -M myemail\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -e error_file\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -q all.q\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -o ",$executionPreferencesRef->{"LOG_FILE"},"\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "set out_dir = ",$path,"/output\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "set scratch = /scratch/",$id,"/",$executionUID,"\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "echo \'#########################\'\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "echo Running on host `hostname`\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "echo Time is `date`\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "echo Directory is ",$path,"\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "echo temp Directory is \$scratch\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "mkdir \$scratch\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$settingsFilePath," \$scratch/settings.dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$randomFilePath," \$scratch/random.dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"INPUT_FILE_PATH"}," ",$path,"/input/\n";
    if ($executionPreferencesRef->{"USE_INITIAL_CONDITIONS"} != 3)
    { 
      print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_PATH"}," ",$path,"/input/\n";
    }
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"GO_PATH"}," \$scratch\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"INPUT_FILE_PATH"}," \$scratch\n";
    if ($executionPreferencesRef->{"USE_INITIAL_CONDITIONS"} != 3)
    {
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_PATH"}," \$scratch\n";
    }
    
    if ($executionPreferencesRef->{"RUN_POST_EXECUTION_ANALYSIS"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"POST_EXECUTION_ANALYSIS_PROG"}," \$scratch\n";
      print EXECUTION_SCRIPT_FILE_HANDLE "chmod 0755 ",$executionPreferencesRef->{"POST_EXECUTION_ANALYSIS_PROG"}," \$scratch\n";
      print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$mdWrapperPrefsFileName," \$scratch\n";
    }
    
    print EXECUTION_SCRIPT_FILE_HANDLE "cd \$scratch\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "./MD.exe\n";

    if ($executionPreferencesRef->{"RUN_POST_EXECUTION_ANALYSIS"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE $executionPreferencesRef->{"POST_EXECUTION_ANALYSIS_PROG"}," $executionUID $path\n";
    }
	
    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 DistanceFile_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp DistanceFile_",$executionUID,".dat.bz2 \$out_dir/Distance\n";

    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 TemperatureFile_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp TemperatureFile_",$executionUID,".dat.bz2 \$out_dir/Temperature\n";

    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 2body_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp 2body_",$executionUID,".dat.bz2 \$out_dir/2Body\n";

    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 3body_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp 3body_",$executionUID,".dat.bz2 \$out_dir/3Body\n";

    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 Etotal_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp Etotal_",$executionUID,".dat.bz2 \$out_dir/Etotal\n";

    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 EbyType_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp EbyType_",$executionUID,".dat.bz2 \$out_dir/EbyType\n";

     print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 ",$executionUID,".log\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionUID,".log.bz2 \$out_dir/Log\n";


     if ($executionPreferencesRef->{"PULL_SPRINGS"} =~ m/YES/i)
     {
     print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 FmovingSpring_",$executionUID,".dat\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp FmovingSpring_",$executionUID,".dat.bz2 \$out_dir/Force_Spring\n";
     
     print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 FstaticSpring_",$executionUID,".dat\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp FstaticSpring_",$executionUID,".dat.bz2 \$out_dir/Force_Spring\n";
     }
     print EXECUTION_SCRIPT_FILE_HANDLE "cp Final_",$executionUID,".dat \$out_dir/LastCord\n";

    if ($executionPreferencesRef->{"STORE_TRAJECTORIES"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 Traj_",$executionUID,".dat\n";
      print EXECUTION_SCRIPT_FILE_HANDLE "cp Traj_",$executionUID,".dat.bz2 \$out_dir/Traj\n";

      if ($executionPreferencesRef->{"STORE_TRAJECTORY_DISTANCES"} =~ m/YES/i)
      {
        print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 Traj_Dist_",$executionUID,".dat\n";
        print EXECUTION_SCRIPT_FILE_HANDLE "cp Traj_Dist_",$executionUID,".dat.bz2 \$out_dir/TrajDist\n";
      }
    }
    print EXECUTION_SCRIPT_FILE_HANDLE "rm -r \$scratch\n";

    close(EXECUTION_SCRIPT_FILE_HANDLE);
    chmod(0755,$cshScript);
    
    if ($executionPreferencesRef->{"EXECUTION_LOCATION"} =~ m/Cluster/i)
    {
      chomp(my $totalJobsWaitingInLine = `qstat | grep " qw " | wc -l`);
      if (($totalJobsWaitingInLine != 0) or ($executionPreferencesRef->{"USE_EMPTY_SLOTS"} =~ m/NO/i))
      {
        #enforce max processors allocation
        chomp(my $myJobs = `qstat | grep " $id " | wc -l`);
        while ($myJobs >= $executionPreferencesRef->{"MAX_PROCESSORS_ALLOCATION"})
        {
            #if there are queued jobs and you use at least $maxProcessorsAllocation - sleep!
            sleep (60);
            #if there are no queued jobs leave the loop #### 3.1 ###
            chomp($totalJobsWaitingInLine = `qstat | grep " qw " | wc -l`);#### 3.1 ###
            last unless ($totalJobsWaitingInLine > 0);#### 3.1 ###
            chomp($myJobs = `qstat | grep " $id " | wc -l`);#### 3.1 ###
        }
      }
      open(MDW_LOCK_FILE_HANDLE,"/home/arielaz/scripts/MDWrapper/MDW.lock") or die "can't open the lock file: $!";
      flock(MDW_LOCK_FILE_HANDLE, 2) or die "can't lock the MDW lock file: $!";
      `cp $cshScript  $path/input/$cshScript`;
      open (PIPE, "|qsub $cshScript") or die "Couldn't execute the csh script: $!\n";
      sleep (10);
     flock(MDW_LOCK_FILE_HANDLE, 8) or die "can't release the MDW lock file: $!";
      close(MDW_LOCK_FILE_HANDLE);
    }
    else
    {
      `cp $cshScript  $path/input/$cshScript`;
      open (PIPE, "|$cshScript") or die "Couldn't execute the csh script: $!\n";
    }
    close PIPE;
    unlink $cshScript;
}
