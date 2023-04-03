!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Distances calculate the distance between any two beads
!***********************************************************************
      subroutine extforce
      include 'MD.com'
        
          
          integer IN_B2, IN_B4, IN_B22, IN_B44
!	  real r1_x, r2_x, r1_y, r2_y, r1_z, r2_z, r1_mod, r2_mod

c          write(87,*)Fext,FX(firstDNAbead),FX(LastDNAbead-2)

c          Lbead=((LastDNAbead+1)/2)+1

	   IN_B2=2732
	   IN_B22=4163	
	   IN_B4=1796
           IN_B44=5078
	
	   FX(IN_B2)=FX(IN_B2)+Fext*0.605
	   FY(IN_B2)=FY(IN_B2)-Fext*0.203
	   FZ(IN_B2)=FZ(IN_B2)-Fext*0.769

	
	   FX(IN_B22)=FX(IN_B22)+Fext*0.536
	   FY(IN_B22)=FY(IN_B22)-Fext*0.466
	   FZ(IN_B22)=FZ(IN_B22)-Fext*0.703
	   
          
	   FX(IN_B4)=FX(IN_B4)+Fext*0.547
	   FY(IN_B4)=FY(IN_B4)-Fext*0.017
	   FZ(IN_B4)=FZ(IN_B4)-Fext*0.836
          

	   FX(IN_B44)=FX(IN_B44)+Fext*0.585
	   FY(IN_B44)=FY(IN_B44)-Fext*0.524
	   FZ(IN_B44)=FZ(IN_B44)-Fext*0.618

c           if(X(Lbead).ge.X(firstDNAbead))then
c          FX(Lbead) = FX(Lbead)+Fext
c          FX(firstDNAbead)=FX(firstDNAbead)-Fext
c          else
c          FX(Lbead) = FX(Lbead)-Fext
c          FX(firstDNAbead)=FX(firstDNAbead)+Fext
c          endif

         end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF DISTANCES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
