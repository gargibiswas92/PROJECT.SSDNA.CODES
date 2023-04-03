!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Distances calculate the distance between any two beads
!***********************************************************************
      subroutine extforce
      include 'MD.com'
        
          
          integer IN_B2, IN_B4, IN_B22, IN_B44
!	  real r1_x, r2_x, r1_y, r2_y, r1_z, r2_z, r1_mod, r2_mod

c          write(87,*)Fext,FX(firstDNAbead),FX(LastDNAbead-2)

c          Lbead=((LastDNAbead+1)/2)+1

	   IN_B2=1913
	   IN_B22=3359	
	   IN_B4=2411
           IN_B44=2855
	
	   FX(IN_B2)=FX(IN_B2)+Fext*0.669
	   FY(IN_B2)=FY(IN_B2)-Fext*0.078
	   FZ(IN_B2)=FZ(IN_B2)-Fext*0.738

	
	   FX(IN_B22)=FX(IN_B22)+Fext*0.476
	   FY(IN_B22)=FY(IN_B22)-Fext*0.603
	   FZ(IN_B22)=FZ(IN_B22)-Fext*0.639
	   
          
	   FX(IN_B4)=FX(IN_B4)+Fext*0.628
	   FY(IN_B4)=FY(IN_B4)-Fext*0.169
	   FZ(IN_B4)=FZ(IN_B4)-Fext*0.758
          

	   FX(IN_B44)=FX(IN_B44)+Fext*0.429
	   FY(IN_B44)=FY(IN_B44)-Fext*0.671
	   FZ(IN_B44)=FZ(IN_B44)-Fext*0.604

c           if(X(Lbead).ge.X(firstDNAbead))then
c          FX(Lbead) = FX(Lbead)+Fext
c          FX(firstDNAbead)=FX(firstDNAbead)-Fext
c          else
c          FX(Lbead) = FX(Lbead)-Fext
c          FX(firstDNAbead)=FX(firstDNAbead)+Fext
c          endif

         end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF DISTANCES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
