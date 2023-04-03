!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Distances calculate the distance between any two beads
!***********************************************************************
      subroutine extforce
      include 'MD.com'
        
          

c          write(87,*)Fext,FX(firstDNAbead),FX(LastDNAbead-2)


          if(X(LastDNAbead-2).ge.X(firstDNAbead))then
          FX(LastDNAbead-2) = FX(LastDNAbead-2)+Fext
          FX(firstDNAbead)=FX(firstDNAbead)-Fext
          else
          FX(LastDNAbead-2) = FX(LastDNAbead-2)-Fext
          FX(firstDNAbead)=FX(firstDNAbead)+Fext
          endif

         end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF DISTANCES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
