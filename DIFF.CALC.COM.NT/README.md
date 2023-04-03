# Workflow to calculate the diffusion coefficient: (All yoy need is the .dat obtained from the trajectory)  
1. First run the code calc_NT_com1.f (You will obviously run the executable, calc_NT_com1). It will calculate the position of the center of mass of the protein, and will calculate the distance of the each of the nucleotides from the protein center of mass. It will create two output files, one having the x-, y-, z- coordinates of the protein center of mass at each 1000 MD steps. And in the other output file it will write the nucleotide index of the DNA closest to the protein COM.  
  You will see some magic numbers here. You might have to change it for different systems.
  * The bead number of the last DNA bead.  
  * The total number of MD time steps.  
  * The total number of simulations ran. In other words the number of trajectory files present in the directory.  
  * The bead number of the last protein bead.  
 2. Then you will run the executable of NT_avgbin.f, this calculates the average nucleotide index in 100 consecutive MD steps. This just smoothens the graphs, if we average.  
 3. Then you will run the executable version of NT_MSD_com.f, this will essentially calculates the mean square displacement (MSD) at different time steps.  
 4. Then you run new_graph.py program, to plot the MSD versus time step plot and calculate the diffusion coefficient from the slope. Remember, I have some specific values in both the filenames inside this code. you might need to change the filenames while you run.  
 5. After getting the value of the diffusion coefficient you might want to get the number of intersegmental transfers obtained from the trajectory. For that will will need to run NT_int_trans code. This will give you the difference of nucleotide index between 1000 consecutive MD steps.  
 6. Then you might run the script cat_file to concatenate all the intersegmental transfer files obtained from running the previous step.  
 7. Also, you will want to run the script file_t2p2 which will run  TrajToPDB3.pl to generate .pdb files from the .dat file.
