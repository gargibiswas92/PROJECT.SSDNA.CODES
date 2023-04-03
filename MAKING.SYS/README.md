# Necessary codes for ssDNA project
1. We started from a PDB file with 60 nt DNA and the RPA protein.  
2. Then I modified the PDB file to add more nucleotides at the two ends of the 60 nt DNA, to make its size 235 nt, 460 nt, 860 nt and 1400 nt respectively.  
3. Then I generated the CG models for each of the cases.  
4. Then I applied some pulling forces at two ends of the newly generated DNA, to get desired conformations.  
(sometimes I applied some more forces to other parts of the DNA to get the desired conformation specially for long DNA)  
5. After getting the desired conformations I modified some contacts and some dihedrals.  
6. I modified the x-, y- and z- axis of the structure by combined rotation and translation matrices.
7. After I obtained the final structure I also modified the .dat file (dat file contains all information regarding co-ordinates, bonds, angles, dihedrals, contacts, repulsions and electrostatics) to incorporate the new co-ordinates into it.   
8. After that, in case of 235 nt DNA I applied another round of pulling at the two ends of the DNA, to obtain structure representing a force-strained DNA.  9. I got four more structures like this. In each cases the end-to-end distances are different. The lower force means lower end-to-end distances and for higher end-to-end distance represents higher force on the DNA.  
10. With all these structures I ran simulations using MD codes to understand the effect of force on DNA.  
11. I ran simulations for 5x10^7 MD time steps for each cases.  
12. Then I calculated the diffusion coefficient using the codes in the /RESULTS/ directory.  


