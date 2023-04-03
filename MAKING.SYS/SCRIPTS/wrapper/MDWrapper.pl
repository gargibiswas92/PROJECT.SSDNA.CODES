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

# print "debug"
my ($scriptDirectory) = $0 =~ m/(.*\/)/g;
my $utilDirectory = $scriptDirectory."../util/";

my $mdWrapperSchemaFileName = $scriptDirectory."/MDWrapper.prefs.schema";
my $mdWrapperPrefsFileName = "./MDWrapper.prefs";
require $utilDirectory."ExecutionPreferences.pm";
my $executionPreferencesRef = ExecutionPreferences::getExecutionPreferences($mdWrapperSchemaFileName,$mdWrapperPrefsFileName);
processExecutionPreferences($executionPreferencesRef);
createOutputFolders();
my $laterCopyScript = $executionPreferencesRef->{"PATH"}."/".$executionPreferencesRef->{"RANDOM_SEED"}.".cpy";
launchExecutions();
`cp $executionPreferencesRef->{"LOG_FILE"} $executionPreferencesRef->{"PATH"}`; 

sub processExecutionPreferences
{
    my $executionPreferencesRef = $_[0];
    chomp(my $id = `whoami`);
    $executionPreferencesRef->{"ID"} = $id;
    $executionPreferencesRef->{"RANDOM_SEED"} = time() % 100000;

  ($executionPreferencesRef->{"INPUT_FILE_NAME"}) = $executionPreferencesRef->{"INPUT_FILE_PATH"} =~ m/.*\/(.*)/;
  ($executionPreferencesRef->{"INPUT_FILE_PREFIX"}) = $executionPreferencesRef->{"INPUT_FILE_NAME"} =~ m/(.*)\./;
  
    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
  {
    ($executionPreferencesRef->{"CONTACT_RANGES_FILE_NAME"}) = $executionPreferencesRef->{"CONTACT_RANGES_FILE_PATH"} =~ m/.*\/(.*)/;
  }

  if ($executionPreferencesRef->{"USE_INITIAL_CONDITIONS"} != 3)
  {
    ($executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_NAME"}) = $executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_PATH"} =~ m/.*\/(.*)/;
  }
  $executionPreferencesRef->{"TEMPERATURES_BUFFER"} = $executionPreferencesRef->{"TEMPERATURES"};
  my @temperatures = split(/\,/,$executionPreferencesRef->{"TEMPERATURES_BUFFER"});
  $executionPreferencesRef->{"TEMPERATURES"} = \@temperatures;
  
  use File::Spec;
  $executionPreferencesRef->{"BASE_PATH"} = File::Spec->rel2abs( $executionPreferencesRef->{"BASE_PATH"} );

  
  if ($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ m/YES/i)
  {
    my @dynamicAtomRange = split(/\,/,$executionPreferencesRef->{"DYNAMIC_ATOMS_RANGE"});
    $executionPreferencesRef->{"DYNAMIC_ATOMS_RANGES_SIZE"} = scalar(@dynamicAtomRange);
    die "number of values in dynamic atom range is not even" 
    unless $executionPreferencesRef->{"DYNAMIC_ATOMS_RANGES_SIZE"} % 2 == 0;
    $executionPreferencesRef->{"DYNAMIC_ATOMS_RANGE"} = "\"".$executionPreferencesRef->{"DYNAMIC_ATOMS_RANGE"}."\"";
  }
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
    mkdir($path.'/output/2Body',0755) unless -d $path.'/output/2Body';
    -d $path.'/output/2Body' or die "Could not create output directory: 2Body\n";
    mkdir($path.'/output/3Body',0755) unless -d $path.'/output/3Body';
    -d $path.'/output/3Body' or die "Could not create output directory: 3Body\n";
    mkdir($path.'/output/Etotal',0755) unless -d $path.'/output/Etotal';
    -d $path.'/output/Etotal' or die "Could not create output directory: Etotal\n";
    mkdir($path.'/output/EbyType',0755) unless -d $path.'/output/EbyType';
    -d $path.'/output/EbyType' or die "Could not create output directory:EbyType\n";
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

    if ($executionPreferencesRef->{"WRITE_ALL_CONTACTS"} =~ m/YES/i)
    {
      mkdir($path.'/output/AllContacts',0755) unless -d $path.'/output/AllContacts';
      -d $path.'/output/AllContacts' or die "Could not create output directory:AllContacts\n";
    }

    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
    {
      mkdir($path.'/output/ContactRanges',0755) unless -d $path.'/output/ContactRanges';
      -d $path.'/output/ContactRanges' or die "Could not create output directory:ContactRanges\n";
    }

}

sub createSettingsFile
{
  my ($currentTemperature, $randomNumber) = @_;
  my $settingsFileName = 'settings_'.$randomNumber.'.dat';
  my $executionUID = $executionPreferencesRef->{"INPUT_FILE_PREFIX"}."_".$currentTemperature."_".$randomNumber;
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
  print SETTINGS_FILE_HANDLE " : initial coordinates/velocities, final conformation\n";
  print SETTINGS_FILE_HANDLE $executionUID.".log : log file for current execution\n";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_DUAL_GAUSSIAN_CONTACTS"};
  print SETTINGS_FILE_HANDLE " : dual Gausian potential contacts\n";
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_DUAL_GAUSSIAN_DIHEDRALS"};
  print SETTINGS_FILE_HANDLE " : dual Gaussian potential dihedral angles\n";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"OUTPUT_MINIMALIST_TRAJ"};
  print SETTINGS_FILE_HANDLE " : Whether to output minimalists Traj file (xyz only)\n";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"WRITE_3BODY"};
  print SETTINGS_FILE_HANDLE " : whether to calculate and write 3body interactions\n";

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
    print SETTINGS_FILE_HANDLE "NO, 1000000";
    print SETTINGS_FILE_HANDLE " : Trajectory file name, Trajectory output frequency(in integer number of time steps)\n";
  }
  print SETTINGS_FILE_HANDLE "Etotal_".$executionUID.".dat, EbyType_".$executionUID.".dat";
  print SETTINGS_FILE_HANDLE " : Total potential energy, Energy by term\n";
  print SETTINGS_FILE_HANDLE "2body_".$executionUID.".dat, 3body_".$executionUID.".dat";
  print SETTINGS_FILE_HANDLE " : 2 body output file, 3 body output file\n";
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"WRITE_ALL_CONTACTS"};
  print SETTINGS_FILE_HANDLE " : whether to write explicit list of all contacts\n";
  if ($executionPreferencesRef->{"WRITE_ALL_CONTACTS"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE "allContacts_".$executionUID.".dat";
    print SETTINGS_FILE_HANDLE " : allContacts file\n";
  }
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"CONTACT_RANGES_DEFINED"};
  print SETTINGS_FILE_HANDLE " : whether to write contact ranges files\n";
   
  if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"CONTACT_RANGES_FILE_NAME"};
    print SETTINGS_FILE_HANDLE " : contact ranges file\n";
    print SETTINGS_FILE_HANDLE "2body_ranges_".$executionUID.".dat";
    print SETTINGS_FILE_HANDLE " : 2body ranges file\n";
    print SETTINGS_FILE_HANDLE "Energy_ranges_".$executionUID.".dat";
    print SETTINGS_FILE_HANDLE " : Energy ranges file\n";
  }
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

  print SETTINGS_FILE_HANDLE $currentTemperature." : temperature for this execution\n";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"HAS_STATIC_ATOMS"};
  print SETTINGS_FILE_HANDLE " : whether there are static atoms in the simulation\n";
  if ($executionPreferencesRef->{"HAS_STATIC_ATOMS"} =~ m/YES/i)
  {
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"DYNAMIC_ATOMS_RANGES_SIZE"};
    print SETTINGS_FILE_HANDLE " : number of values in dynamic atoms range array\n";
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"DYNAMIC_ATOMS_RANGE"};
    print SETTINGS_FILE_HANDLE " : dynamic atoms range array\n";
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"EXCLUDE_STATIC_ATOMS_FROM_TRAJ"};
    print SETTINGS_FILE_HANDLE " : whether to exclude static atoms coordinates from the trajectories \n";
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

    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"COMPENSATE_ELECTROSTATIC_CONTACTS"};
      print SETTINGS_FILE_HANDLE " : whether to fix potential of electrostaic contacts\n";

    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_DEBYE_HUCKEL"};
    print SETTINGS_FILE_HANDLE " : whether to use debye huckel in execution\n";
  
    if ($executionPreferencesRef->{"USE_DEBYE_HUCKEL"} =~ m/YES/i)
    {
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_DH_ENERGY_TABLE"};
      print SETTINGS_FILE_HANDLE " : whether to use precalculated energy and force values\n";
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"IONIC_STRENGTH"};
      print SETTINGS_FILE_HANDLE " : ionic strength of solution\n";
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"IONIC_RADIUS"};
      print SETTINGS_FILE_HANDLE " : ionic radius of ions in solution (1.4 for NaCl) \n";
      print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"SOLVENT_DENSITY"};
      print SETTINGS_FILE_HANDLE " : solventDensity (1 for water) \n";
      # print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_ELECTROSTATICS_GROUPS"};
      # print SETTINGS_FILE_HANDLE " : whether to use ES groups \n";
    }
  }
  
  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_CHIRALS"};
  print SETTINGS_FILE_HANDLE " : whether to use chirals in execution\n";

  print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_ELLIPSOID_REPULSIONS"};
  print SETTINGS_FILE_HANDLE " : whether to use ellipsoid repulsions in execution\n";

  # print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"USE_SLIDING_OBSTACLES"};
  # print SETTINGS_FILE_HANDLE " : whether to add sliding obstacles for protein sliding on DNA\n";
  # if ($executionPreferencesRef->{"USE_SLIDING_OBSTACLES"} =~ m/YES/i)
  #   {
  #     print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"NUMBER_OF_OBSTACLES"};
  #     print SETTINGS_FILE_HANDLE " : Number of Obstacles on DNA\n";
  #     print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"DNA_BASE_PAIR_LENGTH"};
  #     print SETTINGS_FILE_HANDLE " : Base pairs in static DNA\n";
  #     print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"OBSTACLES_BP_SIZE"};
  #     print SETTINGS_FILE_HANDLE " : Obstacle size on DNA\n";
  #     if ($executionPreferencesRef->{"OBSTACLES_BP_SIZE"}>$executionPreferencesRef->{"DNA_BASE_PAIR_LENGTH"}){
  #   die "\n\nFailed!! Obstacle Size is Larger than DNA Size!\n\n";
  #     }
    print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"SIGMA_ELECTROSTATIC_REPULSION_OF_OBSTACLE_PHOSPHATES"};
    print SETTINGS_FILE_HANDLE " : Sigma, Electrostatic repulsion strength by Obstacle Phosphates\n";
    # print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"DIFFUSIVE_OBSTACLES"};
    # print SETTINGS_FILE_HANDLE " : whether to move sliding obstacles on DNA during the run\n";
  #   if ($executionPreferencesRef->{"DIFFUSIVE_OBSTACLES"} =~ m/YES/i)
  #   {
	# print SETTINGS_FILE_HANDLE $executionPreferencesRef->{"DIFFUSION_TIMESTEP"};
  #       print SETTINGS_FILE_HANDLE " : Timestep window for Obstacle move\n";
  #   }
  # }

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
    print MDW_LOG_FILE_HANDLE localtime(time),":Base seed for simulation:", $executionPreferencesRef->{"RANDOM_SEED"},"\n";
    close(MDW_LOG_FILE_HANDLE);
}

sub launchExecutions
{
  foreach (1..$executionPreferencesRef->{"ITERATIONS_PER_TEMPERATURE"})
  {
    foreach my $currentTemprature (@{$executionPreferencesRef->{"TEMPERATURES"}})
    {
        my $randomNumber = int rand $executionPreferencesRef->{"RANDOM_SEED"};
        my $settingsFileName = createSettingsFile($currentTemprature, $randomNumber);
        my $randomFileName = createRandomFile($executionPreferencesRef->{"RANDOM_SEED"},$randomNumber);
        launchExecution($currentTemprature, $settingsFileName, $randomNumber, $randomFileName);
        reportExecution($randomNumber);
    }
  }
}

sub launchExecution
{
    my ($currentTemprature, $settingFileName, $randomNumber, $randomFileName) = @_;
    my $settingsFilePath = $executionPreferencesRef->{"PATH"}.'/input/'.$settingFileName;
    my $randomFilePath = $executionPreferencesRef->{"PATH"}.'/input/'.$randomFileName;
    my $executionUID = $executionPreferencesRef->{"INPUT_FILE_PREFIX"}."_".$currentTemprature."_".$randomNumber;
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
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -o ",$executionPreferencesRef->{"LOG_FILE"},"\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "#\$ -q all.q\n";

    print EXECUTION_SCRIPT_FILE_HANDLE "set hostname = `hostname`\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "set out_dir = ",$path,"/output\n";
    my $perlOutDir = $path."/output";
    my $perlScratch = "/scratch/".$id."/".$executionUID;
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
    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"CONTACT_RANGES_FILE_PATH"}," ",$path,"/input/\n";
    }
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"GO_PATH"}," \$scratch\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"INPUT_FILE_PATH"}," \$scratch\n";
    if ($executionPreferencesRef->{"USE_INITIAL_CONDITIONS"} != 3)
    {
    print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"INITIAL_CONDITIONS_FILE_PATH"}," \$scratch\n";
    }
    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionPreferencesRef->{"CONTACT_RANGES_FILE_PATH"}," \$scratch\n";
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

    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 TemperatureFile_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 2body_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 3body_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 Etotal_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 EbyType_",$executionUID,".dat\n";
    print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 ",$executionUID,".log\n";

    if ($executionPreferencesRef->{"STORE_TRAJECTORIES"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 Traj_",$executionUID,".dat\n";
      if ($executionPreferencesRef->{"STORE_TRAJECTORY_DISTANCES"} =~ m/YES/i)
      {
        print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 Traj_Dist_",$executionUID,".dat\n";
      }
    }
    if ($executionPreferencesRef->{"WRITE_ALL_CONTACTS"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 allContacts_",$executionUID,".dat\n";
    }
    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 2body_ranges_",$executionUID,".dat\n";
      print EXECUTION_SCRIPT_FILE_HANDLE "bzip2 Energy_ranges_",$executionUID,".dat\n";
    }


   if ($executionPreferencesRef->{"USE_COPY_QUEUE"} =~ m/YES/i)
   {

     print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/2body_$executionUID.dat.bz2 $perlOutDir/2Body\" >> $executionUID.cpy\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/3body_$executionUID.dat.bz2 $perlOutDir/3Body\" >> $executionUID.cpy\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/Etotal_$executionUID.dat.bz2 $perlOutDir/Etotal\" >> $executionUID.cpy\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/EbyType_$executionUID.dat.bz2 $perlOutDir/EbyType\" >> $executionUID.cpy\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/$executionUID.log.bz2 $perlOutDir/Log\" >> $executionUID.cpy\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/Final_$executionUID.dat $perlOutDir/LastCord\" >> $executionUID.cpy\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/TemperatureFile_$executionUID.dat.bz2 $perlOutDir/Temperature\" >> $executionUID.cpy\n";

    if ($executionPreferencesRef->{"STORE_TRAJECTORIES"} =~ m/YES/i)
    {
    print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/Traj_$executionUID.dat.bz2 $perlOutDir/Traj\" >> $executionUID.cpy\n";

      if ($executionPreferencesRef->{"STORE_TRAJECTORY_DISTANCES"} =~ m/YES/i)
      {
        print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/Traj_Dist_$executionUID.dat.bz2 $perlOutDir/TrajDist\" >> $executionUID.cpy\n";
      }
    }
    if ($executionPreferencesRef->{"WRITE_ALL_CONTACTS"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/allContacts_$executionUID.dat.bz2 $perlOutDir/AllContacts\" >> $executionUID.cpy\n";
    }
    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/2body_ranges_$executionUID.dat.bz2 $perlOutDir/ContactRanges\" >> $executionUID.cpy\n";
      print EXECUTION_SCRIPT_FILE_HANDLE "echo \"scp \$hostname\\:$perlScratch/Energy_ranges_$executionUID.dat.bz2 $perlOutDir/ContactRanges\" >> $executionUID.cpy\n";
    }
      print EXECUTION_SCRIPT_FILE_HANDLE "echo \"ssh \$hostname rm -rf $perlScratch\" >> $executionUID.del\n";
      print EXECUTION_SCRIPT_FILE_HANDLE "cat $executionUID.cpy >> $laterCopyScript \n";
      print EXECUTION_SCRIPT_FILE_HANDLE "cat $executionUID.del >> $laterCopyScript.del \n";
      print EXECUTION_SCRIPT_FILE_HANDLE "echo 'END' >> $laterCopyScript.end \n";

   }
   else
   {
     print EXECUTION_SCRIPT_FILE_HANDLE "cp Final_",$executionUID,".dat \$out_dir/LastCord\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp ",$executionUID,".log.bz2 \$out_dir/Log\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp TemperatureFile_",$executionUID,".dat.bz2 \$out_dir/Temperature\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp 2body_",$executionUID,".dat.bz2 \$out_dir/2Body\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp 3body_",$executionUID,".dat.bz2 \$out_dir/3Body\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp Etotal_",$executionUID,".dat.bz2 \$out_dir/Etotal\n";
     print EXECUTION_SCRIPT_FILE_HANDLE "cp EbyType_",$executionUID,".dat.bz2 \$out_dir/EbyType\n";
     if ($executionPreferencesRef->{"USE_SLIDING_OBSTACLES"} =~ m/YES/i){
        print EXECUTION_SCRIPT_FILE_HANDLE "cp AlteredPhosphates.dat \$out_dir/\n";
        }

     if ($executionPreferencesRef->{"STORE_TRAJECTORIES"} =~ m/YES/i)
     {
       print EXECUTION_SCRIPT_FILE_HANDLE "cp Traj_",$executionUID,".dat.bz2 \$out_dir/Traj\n";

       if ($executionPreferencesRef->{"STORE_TRAJECTORY_DISTANCES"} =~ m/YES/i)
       {
        print EXECUTION_SCRIPT_FILE_HANDLE "cp Traj_Dist_",$executionUID,".dat.bz2 \$out_dir/TrajDist\n";
      }
     }
    if ($executionPreferencesRef->{"WRITE_ALL_CONTACTS"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "cp allContacts_",$executionUID,".dat.bz2 \$out_dir/AllContacts\n";
    }
    if ($executionPreferencesRef->{"CONTACT_RANGES_DEFINED"} =~ m/YES/i)
    {
      print EXECUTION_SCRIPT_FILE_HANDLE "cp 2body_ranges_",$executionUID,".dat.bz2 \$out_dir/ContactRanges\n";
      print EXECUTION_SCRIPT_FILE_HANDLE "cp Energy_ranges_",$executionUID,".dat.bz2 \$out_dir/ContactRanges\n";
    }
     print EXECUTION_SCRIPT_FILE_HANDLE "rm -r \$scratch\n";
   }
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
      open(MDW_LOCK_FILE_HANDLE,$scriptDirectory."/MDW.lock") or die "cannot open the lock file: $!";
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
