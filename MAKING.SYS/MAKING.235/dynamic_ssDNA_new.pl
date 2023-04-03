#!/usr/bin/perl -w
use strict;
use Math::Trig;
################################################################################
#   This script will receive a dat file of a protein + single strand DNA       #
#   (& maybe double strand soon..) and will create a dynamic data for the DNA  #
#   as well. it will create bonds between the DNA beads,bond angles            #
#   and repulsions according to the user prefrences                            # 
################################################################################
my $dat_file = shift @ARGV or die "Didn't recieve any arguments\n" ;
my ($basename) = $dat_file =~ /(.*?)\.dat/; exit unless defined $basename;   #extracts the name of the file without the .dat suffix
my $new_basename='dynamic_'.$basename;
my $new_name = $new_basename.'.dat';        # new_name holds the name of the output file
my @arr_of_DNA_beads;
my @rep_arr;
my $static=1000000;                         # default value for the phosphate bead
my %dat_file_map=qw(1 Hookean
                  2 angle_trios
                  3 dihedral
                  4 repulsive
                  5 contacts
                  6 other);         # this hash will indicate our location in the dat file
my $key=1;
my $original_repulsive;
my $coordinate_flag=0;

###########   required input data: ##########################################
open (SOURCE,"$dat_file") or die "can't open input dat file for reading $!\n";
print "\n\nENTER BEAD NUMBER OF THE FIRST PHOSPHATE IN THE ssDNA (INTEGER):";
chomp (my $first_phosphate=<STDIN>);       # this indicates the first line to read from
print "ENTER BEAD NUMBER OF THE LAST BASE IN THE ssDNA (INTEGER):";
chomp (my $last_base=<STDIN>);
### Next input added by Amir at March 2013:
print "\n\nENTER NUMBER OF PROTEIN SUBUNITS:";
chomp (my $n_prot_chains=<STDIN>);
print "\n\nENTER SPRING CONSTANT OF THE BONDS BETWEEN THE BEADS (INTEGER, DEFAULT 100):";
chomp (my $bonds_spring=<STDIN>);
if ($bonds_spring eq ""){
    $bonds_spring=100;
}
print "ENTER SPRING CONSTANT OF THE ANGLE BETWEEN DNA BEADS TRIOS (INTEGER, DEFAULT 20):";
chomp (my $angle_spring=<STDIN>);
if ($angle_spring eq ""){
    $angle_spring=20;
}
print "ENTER SPRING CONSTANT OF THE REPULSION BETWEEN THE BEADS (INTEGER, DEFAULT 1):";
chomp (my $repulsion_spring=<STDIN>);
if ($repulsion_spring eq ""){
    $repulsion_spring=1;
}
print "ENTER EPSILON VALUE FOR BASE-PAIRING LJ INTERACTIONS (INTEGER, DEFAULT 0.5):";
chomp (my $epsilon=<STDIN>);
if ($epsilon eq ""){
    $epsilon=0.5;
}
#print "APPLY Polyanion (P) or Single strand (S) DNA dynamics? (P/S):";   # if you want to model the strand as a polymer or DNA
my $poly_type="S";
# print "HAS STATIC ATOM IN CHAIN? (y/n):";   # if you want to hold your ssDNA chain static in one spot
#chomp (my $has_static=<STDIN>);
#if ($has_static eq 'y'){
#    print "ENTER NUMBER OF STATIC DNA BEAD (INTEGER):";
#    chomp ($static=<STDIN>);
#}
my $has_static="n";

#now we will read the input dat file until the coordinate section
my $trash;
while ($trash = <SOURCE>){
    last if $trash =~/ chains/;  
}
my ($chains_num) = $trash =~ /\s*(\d+)\s+\w+/;
for (my $i=0 ; $i<$first_phosphate+$chains_num ; $i++){  #this skips lines until the first phosphate
    $trash=<SOURCE>;}
for (my $i=$first_phosphate ; $i<=$last_base ; $i++){
    $trash =~ /\s*\d+\s*(\d+)\s*(\w)\s*(\w)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\d+\.\d+)\s*(\w)/;
    $arr_of_DNA_beads[$i]->[1]=$1;   # captures the index of the bead
    $arr_of_DNA_beads[$i]->[2]=$2;   # captures the type of the bead (P,S or B)
    $arr_of_DNA_beads[$i]->[3]=$3;   # captures the nucleotide type (C,G,A,T)
    $arr_of_DNA_beads[$i]->[4]=$4;   # captures the x coordinate of the bead
    $arr_of_DNA_beads[$i]->[5]=$5;   # captures the y coordinate of the bead
    $arr_of_DNA_beads[$i]->[6]=$6;   # captures the z coordinate of the bead
    if ($i==$static){
        $arr_of_DNA_beads[$i]->[7]=$7;}   # captures the mass of the bead
    else {$arr_of_DNA_beads[$i]->[7]=1.000;}
    $arr_of_DNA_beads[$i]->[8]=$8;   #chain id


    $trash=<SOURCE>;
}
close SOURCE;
open (SOURCE,"$dat_file") or die "can't open input dat file for reading $!\n";
open (TARGET,">>","$new_name") or die "can't open output dat file for writting $!\n";

my $hookean_bonded_pairs=$last_base-1;
#########################################################################

generate_stack($first_phosphate,$last_base,@arr_of_DNA_beads);
sub generate_stack
{
my ($fp,$lb,@arr) = @_;
my $oldstack=0;
my $check1;
my $check2;
my $stength_aromatic;
 
#my @trp_resi=(185,266,281,430,487,578,608);
#my @tyr_resi=(16,171,204,231,254,262,315,324,345,369,395,488,500,541,550,557,558,587,592,594,650,669,693);
#my @phe_resi=(60,74,82,106,110,143,155,214,236,291,307,314,338,349,371,437,451,459,467,497,542,612,638,652,661,679);
#my @his_resi=(28,89,138,224,237,243,245,280,322,482,494,574,628,643);

my @trp_resi=(356,515,545,836,943,1112,1171);
my @tyr_resi=(31,330,394,446,492,507,609,627,669,717,768,945,968,1039,1057,1071,1073,1130,1140,1144,1252,1290,1337);
my @phe_resi=(115,141,157,204,212,276,298,413,456,563,593,607,655,677,721,849,876,890,904,962,1041,1178,1228,1256,1274,1310);
my @his_resi=(53,170,266,432,458,470,474,543,623,934,956,1104,1208,1238);


my @stack_prot_resi= @trp_resi;
push(@stack_prot_resi,@tyr_resi,);
push(@stack_prot_resi,@phe_resi);
push(@stack_prot_resi,@his_resi);

my $size_arom=scalar(@stack_prot_resi);

#print "$size_arom , @stack_prot_resi";

my $no_prot_DNA=$size_arom*(($lb+1-$fp)/3);
my $number_nucleotides=($lb+1-$fp)/3;
my $number_stack=$number_nucleotides-1+$no_prot_DNA;
#my $number_stack=$number_nucleotides-1;
print TARGET "         $number_stack stacks.\n";


for(my $i=$fp+2 ; $i<=$lb-3 ;$i=$i+3){
         my $j=$i+3;
	 $check1=$arr[$i]->[3];
         $check2=$arr[$j]->[3];

         $oldstack=$oldstack+1;
	 my $sigma=3;
	 my $GAMMA_s=1.5;
	 my $strength=0.1;
         printf TARGET ("%5d%5d%5d%8.3f%8.3f\n",$oldstack,$i,$j,13.000,$strength);
}


for(my $i=$fp+2 ; $i<=$lb ;$i=$i+3){
	foreach my $j(@stack_prot_resi){
	$oldstack=$oldstack+1;
        foreach my $k(@trp_resi){
	if($j eq $k)
        {$stength_aromatic=3.0;}
        else{$stength_aromatic=2.0}
	last if ($j eq $k)	 
}
	printf TARGET ("%5d%5d%5d%8.3f%8.3f\n",$oldstack,$i,$j,10.000, $stength_aromatic);
}
  }


}


#########################################################################
######## now we will construct a new dat file with the data from the array
######## we will add bonds, bond angles (90 degrees) and repulsion
while (my $str=<SOURCE>){  #a
    
    if ($dat_file_map{$key} eq "Hookean"){ #b
        print TARGET "         $hookean_bonded_pairs Hookean bonded pairs.\n";
        for (my $i=1 ; $i<=$first_phosphate-2 ; $i++){ #c
            $str=<SOURCE>;
            print TARGET "$str";
        }  #c
        generate_hookean($first_phosphate,$last_base,$bonds_spring,@arr_of_DNA_beads);
        $str=<SOURCE>;
        $key=2;
    } #b
    
    if ($dat_file_map{$key} eq "angle_trios"){ #d
        my ($old)= $str =~ /\s*(\d+)/;     #old number of bond angle trios
        my $new=$old+4*($last_base+1-$first_phosphate)/3-3;    #this is the new number of bond angles
        print TARGET "         $new bond angle trios.\n";
        for (my $i=1 ; $i<=$old ; $i++){  #e
            $str=<SOURCE>;
            print TARGET "$str";
        }  #e
        generate_angle_trio($old,$new,$first_phosphate,$last_base,$angle_spring,@arr_of_DNA_beads);
        $str=<SOURCE>;
        if ($poly_type eq "S"){
            $key=3;
        }
        else {$key=6;}       
    }  #d
    
    if ($dat_file_map{$key} eq "dihedral"){
       my ($old)= $str =~ /\s*(\d+)/;     #old number of dihedral angle quartets
       my $new=$old+($last_base+1-$first_phosphate)/3-1;    #new number of dihedrals
       $new=$new+($last_base+1-$first_phosphate)/3-3;


       print TARGET "         $new dihedral quartets.\n";
       for (my $i=1 ; $i<=$old ; $i++){  #e
            $str=<SOURCE>;
            print TARGET "$str";
        }
       $old++;   #from now on this will count the index of the lines
       for (my $i=$first_phosphate+2 ; $i<=$last_base-3 ; $i=$i+3){
        my $dihedral=calculateDihedral($i,$i-1,$i+2,$i+3);
        printf TARGET ("%5d%5d%5d%5d%5d%8.3f%8.3f%8.3f%8.3f\n",$old,$i,$i-1,$i+2,$i+3,$dihedral,0.5,0,0);
        $old++;
       }

	for (my $i=$first_phosphate ; $i<=$last_base-11 ; $i=$i+3){
#       print $i,$i+3,$i+6,$i+9,$arr_of_DNA_beads[$i]->[2], $arr_of_DNA_beads[$i+3]->[2],$arr_of_DNA_beads[$i+6]->[2],$arr_of_DNA_beads[$i+9]->[2],"\n";
#       my $dihedral=calculateDihedral($i,$i+3,$i+6,$i+9);
        printf TARGET ("%5d%5d%5d%5d%5d%8.3f%8.3f%8.3f%8.3f\n",$old,$i,$i+3,$i+6,$i+9,3.140,0.7,0,0);
        $old++;
        }




       $str=<SOURCE>;
       $key=5;
    }
    
    if ($dat_file_map{$key} eq "contacts"){

my @check_cn;
my $counter_cn=1;

       my ($old)= $str =~ /\s*(\d+)/;     #old number of contacts
       my $number_of_nucleotides=($last_base+1-$first_phosphate)/3;
       my $new=$old;
       if ($number_of_nucleotides>3){    #we are counting only contacts between the i'th and i+3'th base
        for (my $i=$first_phosphate; $i<=$last_base-9 ; $i=$i+3){
            my $j=$i+9;
            while ($j<=$last_base){
                $new++;
                $j=$j+3;
            }
        }
       }
       else {$new=$old;}
#       print TARGET "         $new contacts.\n";
       print TARGET "         $old contacts.\n";
#       for (my $i=1 ; $i<=$old ; $i++){  #e
#            $str=<SOURCE>;
#            print TARGET "$str";
#        }



	    open(my $fh, '>', 'abc');


	    while ($str = <SOURCE> and $counter_cn <= $old ){
	    $str =~ /s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)/;
            $check_cn[$counter_cn]->[1]=$1;
	    $check_cn[$counter_cn]->[2]=$2;
            $check_cn[$counter_cn]->[3]=$3;
	    $check_cn[$counter_cn]->[4]=$4;
            $check_cn[$counter_cn]->[5]=$5;
            
        if(($check_cn[$counter_cn]->[2] >=160 and $check_cn[$counter_cn]->[2] <= 170)
 or ($check_cn[$counter_cn]->[3] >=160 and $check_cn[$counter_cn]->[3] <= 170))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }

        if(($check_cn[$counter_cn]->[2] >=217 and $check_cn[$counter_cn]->[2] <= 230)
 or ($check_cn[$counter_cn]->[3] >=217 and $check_cn[$counter_cn]->[3] <= 230))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }

        if(($check_cn[$counter_cn]->[2] >=357 and $check_cn[$counter_cn]->[2] <= 376)
 or ($check_cn[$counter_cn]->[3] >=357 and $check_cn[$counter_cn]->[3] <= 376))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }

        if(($check_cn[$counter_cn]->[2] >=544 and $check_cn[$counter_cn]->[2] <= 556)
 or ($check_cn[$counter_cn]->[3] >=544 and $check_cn[$counter_cn]->[3] <= 556))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }

        if(($check_cn[$counter_cn]->[2] >=692 and $check_cn[$counter_cn]->[2] <= 707)
 or ($check_cn[$counter_cn]->[3] >=692 and $check_cn[$counter_cn]->[3] <= 707))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }


        if(($check_cn[$counter_cn]->[2] >=776 and $check_cn[$counter_cn]->[2] <= 795)
 or ($check_cn[$counter_cn]->[3] >=776 and $check_cn[$counter_cn]->[3] <= 795))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }

        if(($check_cn[$counter_cn]->[2] >=852 and $check_cn[$counter_cn]->[2] <= 870)
 or ($check_cn[$counter_cn]->[3] >=852 and $check_cn[$counter_cn]->[3] <= 870))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }

        if(($check_cn[$counter_cn]->[2] >=919 and $check_cn[$counter_cn]->[2] <= 986)
 or ($check_cn[$counter_cn]->[3] >=919 and $check_cn[$counter_cn]->[3] <= 986))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }


        if(($check_cn[$counter_cn]->[2] >=1269 and $check_cn[$counter_cn]->[2] <= 1280)
 or ($check_cn[$counter_cn]->[3] >=1269 and $check_cn[$counter_cn]->[3] <= 1280))
	{
	print $fh "$check_cn[$counter_cn]->[2] $check_cn[$counter_cn]->[3]\n";
	$check_cn[$counter_cn]->[5]=0.0
        }

	printf TARGET ("%5d%5d%5d%10.3f%9.6f\n",$check_cn[$counter_cn]->[1],$check_cn[$counter_cn]->[2], $check_cn[$counter_cn]->[3],
$check_cn[$counter_cn]->[4], $check_cn[$counter_cn]->[5]);
	 
	last if ($counter_cn eq $old);
 	
        $counter_cn=$counter_cn+1;
        }


       $old++;      #from now on this will count the index of the lines
       for (my $i=$first_phosphate+2; $i<=$last_base-9 ; $i=$i+3){  #this loop jump on the bases
        my $j=$i+9;
        while ($j<=$last_base){
#            printf TARGET ("%5d%5d%5d%10.3f%9.6f\n",$old,$i,$j,25,$epsilon);
            $j=$j+3;
            $old++;
        }
       }
       $str=<SOURCE>;
       $key=4;
    }
       
    if ($dat_file_map{$key} eq "repulsive"){  #f

my @check_rep;
my $counter_rep=1;


#my @trp_resi=(185,266,281,430,487,578,608);
#my @tyr_resi=(16,171,204,231,254,262,315,324,345,369,395,488,500,541,550,557,558,587,592,594,650,669,693);
#my @phe_resi=(60,74,82,106,110,143,155,214,236,291,307,314,338,349,371,437,451,459,467,497,542,612,638,652,661,679);
#my @his_resi=(28,89,138,224,237,243,245,280,322,482,494,574,628,643);

my @trp_resi=(356,515,545,836,943,1112,1171);
my @tyr_resi=(31,330,394,446,492,507,609,627,669,717,768,945,968,1039,1057,1071,1073,1130,1140,1144,1252,1290,1337);
my @phe_resi=(115,141,157,204,212,276,298,413,456,563,593,607,655,677,721,849,876,890,904,962,1041,1178,1228,1256,1274,1310);
my @his_resi=(53,170,266,432,458,470,474,543,623,934,956,1104,1208,1238);


my @stack_prot_resi= @trp_resi;
push(@stack_prot_resi,@tyr_resi,);
push(@stack_prot_resi,@phe_resi);
push(@stack_prot_resi,@his_resi);


        if ($poly_type eq "S"){
            ($original_repulsive) = $str =~ /s*(\d+)/;
        }
	
        my $new_reuplsion_number=generate_repulsion($first_phosphate,$last_base,$original_repulsive,@arr_of_DNA_beads);
        print TARGET "       $new_reuplsion_number repulsive pairs.\n";  
#        for (my $i=1 ; $i<=$original_repulsive ; $i++){   #g
#            $str=<SOURCE>;
#            print TARGET "$str";
#        }  #g
	    while ($str = <SOURCE> and $counter_rep <= $original_repulsive ){
	    $str =~ /s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)/;
            $check_rep[$counter_rep]->[1]=$1;
	    $check_rep[$counter_rep]->[2]=$2;
            $check_rep[$counter_rep]->[3]=$3;
	    $check_rep[$counter_rep]->[4]=$4;
            $check_rep[$counter_rep]->[5]=$5;
   
	   my $DNA_check=$check_rep[$counter_rep]->[3];
	
	   
	  
	  if($DNA_check  >= $first_phosphate and $arr_of_DNA_beads[$DNA_check]->[2] eq 'B'){
		
	  my $aromatic_resi=$check_rep[$counter_rep]->[2];
		
	 foreach my $j(@stack_prot_resi)
	 { 
	 if ($j ne $aromatic_resi){$check_rep[$counter_rep]->[5]=1;}
	 else {$check_rep[$counter_rep]->[5]=0;}
		
         last if ($j eq $aromatic_resi)
	 }

         }

#	my $protein_check=$check_rep[$counter_rep]->[2];
#	if($protein_check  >= 776 and $protein_check <= 799 and $protein_check ne 777 and $protein_check ne 785 and $protein_check ne 793 and $DNA_check >= $first_phosphate)
#	{
#print "$protein_check,$DNA_check,$check_rep[$counter_rep]->[5] \n";
#	 $check_rep[$counter_rep]->[5]=0.0;
#print "$protein_check,$DNA_check,$check_rep[$counter_rep]->[5] \n";
#	sleep(0.1);
#	}


	     printf TARGET ("%8d%5d%5d%10.3f%10.3f\n",$check_rep[$counter_rep]->[1],$check_rep[$counter_rep]->[2], $check_rep[$counter_rep]->[3],$check_rep[$counter_rep]->[4], $check_rep[$counter_rep]->[5]);

	  last if ($counter_rep eq $original_repulsive);

	   $counter_rep=$counter_rep+1;
	    
          }
       
	 for (my $i=$original_repulsive+1 ; $i<=$new_reuplsion_number ; $i++){  #h
            printf TARGET ("%8d%5d%5d%10.3f%10.3f\n",$i,$rep_arr[$i]->[1],$rep_arr[$i]->[2],$rep_arr[$i]->[3],$repulsion_spring);
        }  #h
        $key=6;
        $str=<SOURCE>;
    }   #f
    
    if ($dat_file_map{$key} eq "other"){  #i 
        if ($str =~/atom positions./){   #j
            for (my $i=1 ; $i<=($n_prot_chains+3) ; $i++){  #k
                print TARGET "$str";
                $str=<SOURCE>;
            }   #k
            $coordinate_flag=1;   #this flag indicates that we are in the last part of the dat file
        }   #j
        if (($str =~/repulsive/) and ($poly_type ne "S")){  #l
            $key=4;
            ($original_repulsive) = $str =~ /s*(\d+)/;
        }   #l
        else {  #m
            if ((!$coordinate_flag) or ($has_static eq 'n')){  #n
                print TARGET "$str";
            }   #n
        }  #m
        if ($coordinate_flag){   #o
            if ($has_static eq 'y'){   #p
                for (my $i=1 ; $i<$first_phosphate ;$i++){    #copy the protein coordinates
                    print TARGET "$str";
                    $str=<SOURCE>;
                }  #copy the protein coordinates
                for (my $i=$first_phosphate ; $i<=$last_base ; $i++){  #q
                    if ($i==$static){  #r
                        printf TARGET ("%5d%4d%2s%5s%8.3f%8.3f%8.3f%8.3f\n",$i,$i,$arr_of_DNA_beads[$i]->[2],$arr_of_DNA_beads[$i]->[3],$arr_of_DNA_beads[$i]->[4],$arr_of_DNA_beads[$i]->[5],$arr_of_DNA_beads[$i]->[6],999);
                        $str=<SOURCE>;
                    }   #r
                    else{   #s
                        printf TARGET ("%5d%4d%2s%5s%8.3f%8.3f%8.3f%8.3f\n",$i,$i,$arr_of_DNA_beads[$i]->[2],$arr_of_DNA_beads[$i]->[3],$arr_of_DNA_beads[$i]->[4],$arr_of_DNA_beads[$i]->[5],$arr_of_DNA_beads[$i]->[6],1);
                        $str=<SOURCE>;
                    }   #s
                }   #q
            }  #p
            else {   # no static atoms
                for (my $i=2 ; $i<$first_phosphate ;$i++){    #copy the protein coordinates
                    $str=<SOURCE>;
                    print TARGET "$str";
                }  #copy the protein coordinates
                for (my $i=$first_phosphate ; $i<=$last_base ; $i++){
                    printf TARGET ("%5d%4d%2s%5s%8.3f%8.3f%8.3f%8.3f%3s\n",$i,$i,$arr_of_DNA_beads[$i]->[2],$arr_of_DNA_beads[$i]->[3],$arr_of_DNA_beads[$i]->[4],$arr_of_DNA_beads[$i]->[5],$arr_of_DNA_beads[$i]->[6],1,
$arr_of_DNA_beads[$i]->[8]);
                    $str=<SOURCE>;
                }
            }
        }   #o
    }   #i
} #a
sub generate_hookean
{
    my ($fp,$lb,$spring,@arr) = @_;
    printf TARGET ("%5d%5d%5d%8.3f%8.3f\n",$fp-1,$fp-1,$fp,1.0,0.0);
    for (my $i=$fp ; $i<=$lb-2 ; $i=$i+3)
    {
        my $dist_P_S=($arr[$i]->[4]-$arr[$i+1]->[4])**2+($arr[$i]->[5]-$arr[$i+1]->[5])**2+($arr[$i]->[6]-$arr[$i+1]->[6])**2;
        $dist_P_S=sqrt($dist_P_S);
        printf TARGET ("%5d%5d%5d%8.3f%8.3f\n",$i,$i,$i+1,$dist_P_S,$spring); 
        my $dist_S_B=($arr[$i+1]->[4]-$arr[$i+2]->[4])**2+($arr[$i+1]->[5]-$arr[$i+2]->[5])**2+($arr[$i+1]->[6]-$arr[$i+2]->[6])**2;
        $dist_S_B=sqrt($dist_S_B);
        printf TARGET ("%5d%5d%5d%8.3f%8.3f\n",$i+1,$i+1,$i+2,$dist_S_B,$spring);
        if (($lb-$i)>2)
        {
            my $dist_S_P=($arr[$i+1]->[4]-$arr[$i+3]->[4])**2+($arr[$i+1]->[5]-$arr[$i+3]->[5])**2+($arr[$i+1]->[6]-$arr[$i+3]->[6])**2;
            $dist_S_P=sqrt($dist_S_P);
            printf TARGET ("%5d%5d%5d%8.3f%8.3f\n",$i+2,$i+1,$i+3,$dist_S_P,$spring);
        }
    }
}

sub generate_angle_trio
{
    my ($old,$new,$fp,$lb,$spring,@arr) = @_;
    my $bb_index=$fp+1;   #this is a backbone index which pass from sugar to phophate on the backbone
    my $cos_angle1;        #we will calculate the angles between vectors with vector multiplication
    my $angle1;
    my $cos_angle2;        #we will calculate the angles between vectors with vector multiplication
    my $angle2;
    my $l1;               #we will store here the length of a vector
    my $l2;
    my $l3;
    my $formed_PiSiBi_angle=0;
    my $i=$old+1;
    #for (my $i=$old+1 ; $i<=$new ; $i++){
    while ($i<=$new){    
        if (($arr[$bb_index]->[2] eq "S") and (!$formed_PiSiBi_angle)){    # we are on a sugar and need to create Pi-Si-Bi angle
            #we will calculate the vector from this sugar to its phophate (in the same nucleotide_
            my @vector1=($arr[$bb_index]->[4]-$arr[$bb_index-1]->[4],$arr[$bb_index]->[5]-$arr[$bb_index-1]->[5],$arr[$bb_index]->[6]-$arr[$bb_index-1]->[6]);   
            #and the vector from this sugar to the preceding base
            my @vector2=($arr[$bb_index]->[4]-$arr[$bb_index+1]->[4],$arr[$bb_index]->[5]-$arr[$bb_index+1]->[5],$arr[$bb_index]->[6]-$arr[$bb_index+1]->[6]);
            #we will calculate the length of the vectors
            $l1=length_vect(@vector1);
            $l2=length_vect(@vector2);
            $cos_angle1=($vector1[0]*$vector2[0]+$vector1[1]*$vector2[1]+$vector1[2]*$vector2[2])/($l1*$l2);
            $angle1=acos($cos_angle1);   # this angle is Pi-Si-Bi
            printf TARGET ("%5d%5d%5d%5d%8.3f%8.3f\n",$i,$bb_index-1,$bb_index,$bb_index+1,$angle1,$spring);
            $i++;
            $formed_PiSiBi_angle=1;
        }
        if (($arr[$bb_index]->[2] eq "S") and ($formed_PiSiBi_angle)){    # we are on a sugar and need to create Bi-Si-Pi+1 and Pi-Si-Pi+1 angles (if exists Pi+1)
            if ($bb_index<$lb-1){     #if we are in the last sugar in the backbone we dont need to do this calculations
            #we will calculate the vector from this sugar to the next phophate (in the next nucleotide)   Si-Pi+1
            my @vector1=($arr[$bb_index]->[4]-$arr[$bb_index+2]->[4],$arr[$bb_index]->[5]-$arr[$bb_index+2]->[5],$arr[$bb_index]->[6]-$arr[$bb_index+2]->[6]);   
            #and the vector from this sugar to the preceding base                                         Si-Bi
            my @vector2=($arr[$bb_index]->[4]-$arr[$bb_index+1]->[4],$arr[$bb_index]->[5]-$arr[$bb_index+1]->[5],$arr[$bb_index]->[6]-$arr[$bb_index+1]->[6]);
            #and the vector from this sugar to the its phophate (in the same nucleotide)                  Si-Pi
            my @vector3=($arr[$bb_index]->[4]-$arr[$bb_index-1]->[4],$arr[$bb_index]->[5]-$arr[$bb_index-1]->[5],$arr[$bb_index]->[6]-$arr[$bb_index-1]->[6]);
            #we will calculate the length of the vectors
            $l1=length_vect(@vector1);
            $l2=length_vect(@vector2);
            $l3=length_vect(@vector3);
            $cos_angle1=($vector1[0]*$vector2[0]+$vector1[1]*$vector2[1]+$vector1[2]*$vector2[2])/($l1*$l2);
            $angle1=acos($cos_angle1);   # this angle is Bi-Si-Pi+1
            $cos_angle2=($vector1[0]*$vector3[0]+$vector1[1]*$vector3[1]+$vector1[2]*$vector3[2])/($l1*$l3);
            $angle2=acos($cos_angle2);   # this angle is Pi-Si-Pi+1
            printf TARGET ("%5d%5d%5d%5d%8.3f%8.3f\n",$i,$bb_index+1,$bb_index,$bb_index+2,$angle1,$spring);
            $i++;
            printf TARGET ("%5d%5d%5d%5d%8.3f%8.3f\n",$i,$bb_index-1,$bb_index,$bb_index+2,$angle2,$spring);
            $i++;
            $formed_PiSiBi_angle=0;
            $bb_index +=2;
            }
        }
        if ($arr[$bb_index]->[2] eq "P"){    # we are on a phophate and need to create the angle Si-1-Pi-Si
            #we will calculate the vector from this phophate to the previous sugar (in the previous nucleotide)
            my @vector1=($arr[$bb_index]->[4]-$arr[$bb_index-2]->[4],$arr[$bb_index]->[5]-$arr[$bb_index-2]->[5],$arr[$bb_index]->[6]-$arr[$bb_index-2]->[6]);   
            #and the vector from this sugar to the preceding base
            my @vector2=($arr[$bb_index]->[4]-$arr[$bb_index+1]->[4],$arr[$bb_index]->[5]-$arr[$bb_index+1]->[5],$arr[$bb_index]->[6]-$arr[$bb_index+1]->[6]);
            #we will calculate the length of the vectors
            $l1=length_vect(@vector1);
            $l2=length_vect(@vector2);
            $cos_angle1=($vector1[0]*$vector2[0]+$vector1[1]*$vector2[1]+$vector1[2]*$vector2[2])/($l1*$l2);
            $angle1=acos($cos_angle1);   # this angle is Pi-Si-Bi
            printf TARGET ("%5d%5d%5d%5d%8.3f%8.3f\n",$i,$bb_index-2,$bb_index,$bb_index+1,$angle1,$spring);
            $i++;
            $bb_index++;
        }
    }
}
sub generate_repulsion
{   # we create a repulsion for each 2 beads which are at least 4 bonds apart
    my ($fp,$lb,$old,@arr) = @_;    #we will insert the repulsion data to array of the beads
    my $counter=$old+1;
    for (my $i=$fp ; $i<=$lb-3 ; $i++){
        if ($arr[$i]->[2] eq "P"){
            if ($i<=$lb-5){
                for (my $j=$i+5 ; $j<=$lb ; $j++){
                    $rep_arr[$counter]->[1]=$i;
                    $rep_arr[$counter]->[2]=$j;
                    if ($arr[$j]->[2] eq "B"){       #if the J'th bead is a phosphate we use different repulsion distance
                        $rep_arr[$counter]->[3]=27.04;
                    }
                    else {
                        $rep_arr[$counter]->[3]=54.76;
                    }
                    $counter++;
                }
            }
        }
        if ($arr[$i]->[2] eq "S"){
            if ($i<=$lb-6){
                for (my $j=$i+6 ; $j<=$lb ; $j++){
                    $rep_arr[$counter]->[1]=$i;
                    $rep_arr[$counter]->[2]=$j;
                    if ($arr[$j]->[2] eq "B"){       #if the J'th bead is a phosphate we use different repulsion distance
                        $rep_arr[$counter]->[3]=27.04;
                    }
                    else {
                        $rep_arr[$counter]->[3]=54.76;
                    }
                    $counter++;
                }
            }           
        }
        if ($arr[$i]->[2] eq "B"){
            if ($i<=$lb-3){
                for (my $j=$i+3 ; $j<=$lb ; $j++){
		    if($j == $i+3){
		    next ;
		    }
                    do {
                    $rep_arr[$counter]->[1]=$i;
                    $rep_arr[$counter]->[2]=$j;
                    my $bond_gap=$j-$i;
                    if ($bond_gap%3==0){    # this means that i'th and j'th beads are bases
                        $rep_arr[$counter]->[3]=9.00;    # base-base repulsion distance
                    }
                    else {
                        $rep_arr[$counter]->[3]=27.04;   # base-other repulsion distance
                    }   
                    $counter++;
                    }
#unless (($arr[$j]->[2] eq "B") and ($j>=$i+9))   #there is a contact term between bases that are more than 3 nucleic acid away
                }
            }           
        }
    }
    return $counter-1;
}
           
sub length_vect      #this subroutine recieves a vector and returns its length    
{
    my (@arr) = @_;
    my $sum=0;
    foreach my $ind (@arr){
        $sum += $ind**2;
    }
    return sqrt($sum);
}

sub calculateDihedral{
    my ($l1,$l2,$l3,$l4) = @_;
    
    #calculation of the vector from bead 2 to bead 1 (sugar to its base)
    my $X_21 = $arr_of_DNA_beads[$l1]->[4]-$arr_of_DNA_beads[$l2]->[4];
    my $Y_21 = $arr_of_DNA_beads[$l1]->[5]-$arr_of_DNA_beads[$l2]->[5];
    my $Z_21 = $arr_of_DNA_beads[$l1]->[6]-$arr_of_DNA_beads[$l2]->[6];
    
    #calculation of the vector from bead 2 to bead 3 (sugar i'th to sugar i+1'th)
    my $X_23 = $arr_of_DNA_beads[$l3]->[4]-$arr_of_DNA_beads[$l2]->[4];
    my $Y_23 = $arr_of_DNA_beads[$l3]->[5]-$arr_of_DNA_beads[$l2]->[5];
    my $Z_23 = $arr_of_DNA_beads[$l3]->[6]-$arr_of_DNA_beads[$l2]->[6];
    
    #calculation of the vector from bead 4 to bead 3 (base i+1'th to sugar i+1'th)
    my $X_43 = $arr_of_DNA_beads[$l3]->[4]-$arr_of_DNA_beads[$l4]->[4];
    my $Y_43 = $arr_of_DNA_beads[$l3]->[5]-$arr_of_DNA_beads[$l4]->[5];
    my $Z_43 = $arr_of_DNA_beads[$l3]->[6]-$arr_of_DNA_beads[$l4]->[6];

    #vector multiplication to obtain the normal to the surface formed by beads 1,2,3
    my $X_normal_21_23 = $Y_21*$Z_23 - $Z_21*$Y_23;
    my $Y_normal_21_23 = $Z_21*$X_23 - $X_21*$Z_23;
    my $Z_normal_21_23 = $X_21*$Y_23 - $Y_21*$X_23;

    #vector multiplication to obtain the normal to the surface formed by beads 2,3,4
    my $X_normal_43_23 = $Y_43*$Z_23 - $Z_43*$Y_23;
    my $Y_normal_43_23 = $Z_43*$X_23 - $X_43*$Z_23;
    my $Z_normal_43_23 = $X_43*$Y_23 - $Y_43*$X_23;

    #size of the 2 normals
    my $size_normal_21_23 = sqrt($X_normal_21_23**2 + $Y_normal_21_23**2 + $Z_normal_21_23**2 + 1e-24);
    my $size_normal_43_23 = sqrt($X_normal_43_23**2 + $Y_normal_43_23**2 + $Z_normal_43_23**2 + 1e-24);
    my $mul_normals = $X_normal_21_23*$X_normal_43_23 + $Y_normal_21_23*$Y_normal_43_23 + $Z_normal_21_23*$Z_normal_43_23;

    my $mul_normals_factor = 1 / ($size_normal_21_23 * $size_normal_43_23);
    $mul_normals_factor = 0 if (($size_normal_21_23 < 1.0e-3) or ($size_normal_43_23 < 1.0e-3));

    my $cos_dihedral = $mul_normals*$mul_normals_factor;
    $cos_dihedral = ( $cos_dihedral < -1 ? -1 : $cos_dihedral);
    $cos_dihedral = ( $cos_dihedral > 1 ? 1 : $cos_dihedral);
 
    my $mul_23_normal_normals = $X_23*($Z_normal_21_23*$Y_normal_43_23 - $Y_normal_21_23*$Z_normal_43_23)
                              + $Y_23*($X_normal_21_23*$Z_normal_43_23 - $Z_normal_21_23*$X_normal_43_23)
                              + $Z_23*($Y_normal_21_23*$X_normal_43_23 - $X_normal_21_23*$Y_normal_43_23);
    my $dihedral = acos($cos_dihedral);
    $dihedral = ( $mul_23_normal_normals < 0 ? pi + $dihedral : pi - $dihedral);
 
    return $dihedral;
}



