#!/usr/bin/perl
use strict;
#############################################################
#any script calling ParseCRD must includ a function runAnalysis($crd_str ,$number_of_beads,$numOf_dynamic_Protein_beads)
#
#the idea is that each of the 3 sublutins here parse a crd file and for each frame they read they call the subrutin
#runAnalysis (wich must be defined in the file calling the subrutine)with the $crd_str which they created
#
#$crd_str
#     hush{bead_id}->{x}
#                  ->{y}
#                  ->{z}
##################################################################


#any script calling ParseCRD must includ a function runAnalysis($crd_str ,$number_of_aa,$numOfProteinAAs)
sub ParseCRD {
	### the centre subroutine of the code
    
    my ( $crd, $number_of_aa,$numOfProteinAAs) = @_; ##@@ , $csu_str, $cont_in_step,$send_out, $hydro) = @_;
    
    
    open (CRD, "$crd") or die "can't open the crd file: $!\n";
    
    my $crd_str = {};my $cont_str = {};
    
    ### to calculate the amount of lines in a single trajectory:
    ### version 1.1:
    ### my $number_of_aa = (34 * $repeats) + 15;
    my $pdb = <CRD>;
    my $line_number = int (($number_of_aa *3) / 10);
    my $modulus = $number_of_aa % 10;
    ###
    $line_number -=1 if ($modulus == 0);
    ###
    my $traject_number = 0;
    my $crd_line = "";
    my $resi = "";
    my @check = ();
    my $time = localtime time;
    print "--->$time: working on file $crd <---\n";

    
    LOOP: while (<CRD>) {

        $crd_line = "";
        $crd_line .= $_;
        my $count_line = 1;
        my $test_arr = [];
        
        ### for normal lines
        while ($count_line < $line_number) {

            my $line = <CRD>;
            #$crd_line =~ s/\-/ -/g;
          #  print "$crd_line\n";
          @{$test_arr} = split /(-*\d{1,3}\.\d{1,3})/, $line;
          my @clean;
          foreach (@{$test_arr}) {
            if (/(-*\d{1,3}\.\d{1,3})/) {
                push @clean, $_;
            }
        }
          @{$test_arr} = @clean;
            #@{$test_arr} = split /\s+/, $line;            
            #shift @{$test_arr};


            if (scalar @{$test_arr} == 10) {
                ### join the line to $crd_line
                $crd_line .= $line;
                $count_line++;
            }
            
            elsif (scalar @{$test_arr} < 10) {
                ### it is a FUCKED up traject step
                ### need to start the top while again
                ### and do nothing with the collected data

                next LOOP;
            }
        }
        
        ### for last line in NORMAL step
        if ($count_line == $line_number) {
            my $line = <CRD>;
            ### the last line in the traect step
            ### join this one as well
            $crd_line .= $line;
            $traject_number++;

        }
  
        $crd_line =~ s/\n/ /g;
        
        my @crd = split /(-*\d{0,3}\.\d{0,3})/, $crd_line;
        my $crd_new =  [];
        
        ### bugfix
        foreach (@crd) {
            if (/(-*\d{0,3}\.\d{0,3})/) {
                push @{$crd_new}, $_;
            }
        }
        
        #### debug 1.0.4:
        #print scalar @crd, "\n";

        my $resi = 1;
        while (@{$crd_new}) {
            my $x = shift @{$crd_new};
            my $y = shift @{$crd_new};
            my $z = shift @{$crd_new};
            ### debug:
            if ((defined $z) and (defined $y) and (defined $x)) {
            ###
            $crd_str->{"$resi"} = {
                                             'x' => "$x",
                                             'y' => "$y",
                                             'z' => "$z",
                                            };
            ###debug :$resi++;
            ### debug
            }
            ### debug
            else {
                print "---> $resi is fucked in $crd<---\n";
                }
            ###
          #@@  print " lalalala $crd_str->{$resi}->{'x'} $crd_str->{$resi}->{'x'} $crd_str->{$resi}->{'x'}\n";
            $resi++; 
        }
        #print "read trej $traject_number\n";    
        ### From here you can call whatever you'd like that
        ### can use the parsed crd trajectory step
        
        ### Call the subroutine defined in the main that
        ### count the inter and intra contacts within each
        ### conformer along the trajectory
        
        runAnalysis($crd_str ,$number_of_aa,$numOfProteinAAs);
        
        ##@@  $cont_str = $cont_in_step->($native_str, $crd_str, $csu_str, $traject_number)
        ###@@  if defined $cont_in_step;
        
        ### Calculate Hydrodynamic Radius - Rh
        ##@@ $hydro->($crd_str, $resi, $cont_str, $traject_number)
        ###@@@ if defined $hydro;
        
        ### Subroutine defined in main that prints the
        ### contact analysis results to a file
        ###@@ $send_out->($cont_str, $traject_number, $crd)
        ##@@@ if defined $send_out;

        ### somewhere at the end of working with the trajectory:
        ### version 1.0.5: create an escape sub
        
        #if (scalar @check < 10) {
        #        
        #        ($traject_number, $cont_str, $crd_line, $crd_str,
        #        $parse_counter, $resi) =
        #        $skip->('within', $traject_number, $cont_str, $crd_line, $crd_str,
        #        $parse_counter, $resi);
        #        next Loop;
        #};
        
        #($traject_number, $cont_str, $crd_line, $crd_str,
        #$parse_counter, $resi) =
        #$skip->('new', $traject_number, $cont_str, $crd_line, $crd_str,
        #$parse_counter, $resi);
        #
        ### versio 1.0.5:
        #$traject_number++;
        $cont_str = {};
        #$crd_line = "";
        $crd_str = {};
        #$parse_counter = 0;
        $resi = 1;

    }
    #@@ $native_str = {};
}

sub ParsePrmCRD {
	### the centre subroutine of the code
    
    my ( $crd, $number_of_aa,$numOfProteinAAs) = @_; ##@@ , $csu_str, $cont_in_step,$send_out, $hydro) = @_;
    
    
    open (CRD, "$crd") or die "can't open the crd file: $!\n";
    
    my $crd_str = {};my $cont_str = {};
    
    ### to calculate the amount of lines in a single trajectory:
    ### version 1.1:
    ### my $number_of_aa = (34 * $repeats) + 15;
    my $pdb = <CRD>;
    my $atomNum = <CRD>;
    my $line_number = int (($number_of_aa *3) / 6);
    my $modulus = $number_of_aa % 6;
    ###
    $line_number -=1 if ($modulus == 0);
    ###
    my $traject_number = 0;
    my $crd_line = "";
    my $resi = "";
    my @check = ();
    my $time = localtime time;
    print "--->$time: working on file $crd <---\n";

    
    LOOP: while (<CRD>) {

        $crd_line = "";
        $crd_line .= $_;
        my $count_line = 1;
        my $test_arr = [];
        
        ### for normal lines
        while ($count_line < $line_number) {

            my $line = <CRD>;
            #$crd_line =~ s/\-/ -/g;
          #  print "$crd_line\n";
          @{$test_arr} = split /(-*\d{1,3}\.\d{1,3})/, $line;
          my @clean;
          foreach (@{$test_arr}) {
            if (/(-*\d{1,3}\.\d{1,3})/) {
                push @clean, $_;
            }
        }
          @{$test_arr} = @clean;
            #@{$test_arr} = split /\s+/, $line;            
            #shift @{$test_arr};


            if (scalar @{$test_arr} == 6) {
                ### join the line to $crd_line
                $crd_line .= $line;
                $count_line++;
            }
            
            elsif (scalar @{$test_arr} < 6) {
                ### it is a FUCKED up traject step
                ### need to start the top while again
                ### and do nothing with the collected data

                next LOOP;
            }
        }
        
        ### for last line in NORMAL step
        if ($count_line == $line_number) {
            my $line = <CRD>;
            ### the last line in the traect step
            ### join this one as well
            $crd_line .= $line;
            $traject_number++;

        }
  
        $crd_line =~ s/\n/ /g;
        
        my @crd = split /(-*\d{0,3}\.\d{0,3})/, $crd_line;
        my $crd_new =  [];
        
        ### bugfix
        foreach (@crd) {
            if (/(-*\d{0,3}\.\d{0,3})/) {
                push @{$crd_new}, $_;
            }
        }
        
        #### debug 1.0.4:
        #print scalar @crd, "\n";

        my $resi = 1;
        while (@{$crd_new}) {
            my $x = shift @{$crd_new};
            my $y = shift @{$crd_new};
            my $z = shift @{$crd_new};
            ### debug:
            if ((defined $z) and (defined $y) and (defined $x)) {
            ###
            $crd_str->{"$resi"} = {
                                             'x' => "$x",
                                             'y' => "$y",
                                             'z' => "$z",
                                            };
            ###debug :$resi++;
            ### debug
            }
            ### debug
            else {
                print "---> $resi is fucked in $crd<---\n";
                }
            ###
          #@@  print " lalalala $crd_str->{$resi}->{'x'} $crd_str->{$resi}->{'x'} $crd_str->{$resi}->{'x'}\n";
            $resi++; 
        }
        
        runAnalysis($crd_str ,$number_of_aa,$numOfProteinAAs);
 
        $cont_str = {};
        #$crd_line = "";
        $crd_str = {};
        #$parse_counter = 0;
        $resi = 1;

    }
    #@@ $native_str = {};
}

sub parseTrj{
	use DB_File;
	use Fcntl;
	
	my $flag = ''; ### (gets to be either 'continue' or 'end');
	my $rg_str = {};
	my ($traj, $number_of_aa ,$numOfProteinAAs) = @_;
	#my ($basename) = $traj =~ /(.*?)\.dat/; exit unless defined $basename;
        my $line;
                
        open(TRJ,"$traj")or die "cant open input file $traj\n";
        
	### Create an object named 'tie' to hold the file:
	# my $tie = tie(@lines, DB_File, $traj, O_RDWR, 0666, $DB_RECNO) or die "can't openfile $traj: $!\n";
	### Get line number 2 that holds the number of atoms:
        chomp (my $chainNumber  = <TRJ> );
	#chomp(my $number_atoms_ch1 = <TRJ>);
        my $number_atoms;
        for (my $i= 0;$i<$chainNumber;$i++){
            $number_atoms+=<TRJ>;
        }
	#$number_atoms_ch1 =~ s/\s+//;
        die "wrong number of atoms at /home/ohad/scripts/crdHandling.pm -> parseTrj number_of_aa $number_of_aa number_atoms_ch1 $number_atoms\n" unless ($number_of_aa == $number_atoms);
	
	
	my $crd_str = {};
	my $resi = -1;
	my $traject_number =-1;
	
        while ($line = <TRJ>){
            if (my ($model) = $line =~ /^\s+(\d+)$/) {
                # ignor the number on atoms at the end of the atom names section
                if ($traject_number<0){
                    $traject_number++;
                    next;
                }

    		$resi = 1;
                $traject_number++;
                #$traject_number = $model / 1000;
                if($traject_number>1){
                    runAnalysis($crd_str ,$number_of_aa,$numOfProteinAAs); 
                }
		next;
	    }
            ### Build the data structure for a single configuration:
            if (my ($x, $y, $z) = $line =~ /\s*(-*\d+\.\d{3})\s*(-*\d+\.\d{3})\s*(-*\d+\.\d{3})/) {
			$crd_str->{"$resi"} = {
                                'x' => "$x",
                                'y' => "$y",
                                'z' => "$z",              
                                };
			$resi++;
			next;
		}
                
        }
        
	
}
1;