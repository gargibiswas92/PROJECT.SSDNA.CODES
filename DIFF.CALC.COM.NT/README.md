# Workflow to calculate the diffusion coefficient:   
1. First run the code calc_NT_com1.f. It will calculate the position of the center of mass of the protein, and will calculate the distance of the each of the nucleotides from the protein center of mass. It will create two output files, one having the x-, y-, z- coordinates of the protein center of mass at each 1000 MD steps. And in the other output file it will write the nucleotide index of the DNA closest to the protein COM.  
  You will see some magic numbers here. You might have to change it for different systems.
  * The bead number of the last DNA bead.  
  * The total number of MD time steps.  
  * The total number of simulations ran. In other words the number of trajectory files present in the directory.  
  * The bead number of the last protein bead.
