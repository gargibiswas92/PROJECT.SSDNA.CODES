!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Distances calculate the distance between any two beads
!***********************************************************************
      subroutine extforce
      include 'MD.com'
        
          
          integer Lbead, IN_B1, IN_B2, IN_B3, IN_B4
	  real r1_x, r2_x, r1_y, r2_y, r1_z, r2_z, r1_mod, r2_mod

c          write(87,*)Fext,FX(firstDNAbead),FX(LastDNAbead-2)

c          Lbead=((LastDNAbead+1)/2)+1

	   IN_B1=1352
	   IN_B2=1340
	   	
           IN_B3=1505
	   IN_B4=1517


	   r1_x= x(IN_B1)-x(IN_B2)
	   r1_y= y(IN_B1)-y(IN_B2)
	   r1_z= z(IN_B1)-z(IN_B2)
	  
	   r1_mod=sqrt((r1_x)**2+(r1_y)**2+(r1_z)**2)

	   r2_x= x(IN_B3)-x(IN_B4)
	   r2_y= y(IN_B3)-y(IN_B4)
	   r2_z= z(IN_B3)-z(IN_B4)

	   r2_mod=sqrt((r2_x)**2+(r2_y)**2+(r2_z)**2)

	
	   FX(IN_B2)=FX(IN_B2)+Fext*(r1_x/r1_mod)
	   FY(IN_B2)=FY(IN_B2)+Fext*(r1_y/r1_mod)
	   FZ(IN_B2)=FZ(IN_B2)+Fext*(r1_z/r1_mod)

	
	   FX(IN_B4)=FX(IN_B4)+Fext*(r2_x/r2_mod)
	   FY(IN_B4)=FY(IN_B4)+Fext*(r2_y/r2_mod)
	   FZ(IN_B4)=FZ(IN_B4)+Fext*(r2_z/r2_mod)
	   
          
c          if(X(Lbead).ge.X(firstDNAbead))then
          FX(Lbead) = FX(Lbead)+Fext
c          FX(firstDNAbead)=FX(firstDNAbead)-Fext
c          else
c          FX(Lbead) = FX(Lbead)-Fext
c          FX(firstDNAbead)=FX(firstDNAbead)+Fext
c          endif

         end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF DISTANCES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
