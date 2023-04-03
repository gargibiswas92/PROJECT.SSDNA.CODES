!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Distances calculate the distance between any two beads
!***********************************************************************
      subroutine extforce
      include 'MD.com'
        
          
          integer IN_B2, IN_B4, Lbead
!	  real r1_x, r2_x, r1_y, r2_y, r1_z, r2_z, r1_mod, r2_mod

c          write(87,*)Fext,FX(firstDNAbead),FX(LastDNAbead-2)

c          Lbead=((LastDNAbead+1)/2)+1

	   IN_B2=1340
	   	
	   IN_B4=2717

	
	   FX(IN_B2)=FX(IN_B2)-Fext*0.989
	   FY(IN_B2)=FY(IN_B2)+Fext*0.023
	   FZ(IN_B2)=FZ(IN_B2)-Fext*0.143

	
	   FX(IN_B4)=FX(IN_B4)+Fext*0.989
	   FY(IN_B4)=FY(IN_B4)-Fext*0.023
	   FZ(IN_B4)=FZ(IN_B4)+Fext*0.143
	   
          
c          if(X(Lbead).ge.X(firstDNAbead))then
c          FX(Lbead) = FX(Lbead)+Fext
c          FX(firstDNAbead)=FX(firstDNAbead)-Fext
c          else
c          FX(Lbead) = FX(Lbead)-Fext
c          FX(firstDNAbead)=FX(firstDNAbead)+Fext
c          endif

         end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF DISTANCES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
