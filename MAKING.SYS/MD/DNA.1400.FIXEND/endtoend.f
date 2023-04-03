!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Distances calculate the distance between any two beads
!***********************************************************************
      subroutine enddis(disP)
c      subroutine enddis(disP,disS,disB)
      include 'MD.com'
        
      integer sugar,base,phosphate
      real disB,disP,disS
c-------------------------------------------------------------------------------------------------------          

          phosphate=firstDNAbead
!          dx = X(firstDNAbead+6) - X(firstDNAbead)
!          dy = Y(firstDNAbead+6) - Y(firstDNAbead)
!          dz = Z(firstDNAbead+6) - Z(firstDNAbead)

          dx = X(firstDNAbead) - X(LastDNAbead-2)
          dy = Y(firstDNAbead) - Y(LastDNAbead-2)
          dz = Z(firstDNAbead) - Z(LastDNAbead-2)


          disP=sqrt(dx**2+dy**2+dz**2)
          
          sugar=firstDNAbead+1
        
!          dx = X(sugar) - X(sugar+6)
!          dy = Y(sugar) - Y(sugar+6)
!          dz = Z(sugar) - Z(sugar+6)
 
!          dx = X(sugar) - X(LastDNAbead-1)
!          dy = Y(sugar) - Y(LastDNAbead-1)
!          dz = Z(sugar) - Z(LastDNAbead-1)
!          disS=sqrt(dx**2+dy**2+dz**2)
          
        
          base=firstDNAbead+2       

!          dx = X(base) - X(base+6)
!          dy = Y(base) - Y(base+6)
!          dz = Z(base) - Z(base+6)
          
!          dx = X(base) - X(LastDNAbead)
!          dy = Y(base) - Y(LastDNAbead)
!          dz = Z(base) - Z(LastDNAbead)
!          disB=sqrt(dx**2+dy**2+dz**2)
c----------------------------------------------------------------------------------------------------------
      end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF DISTANCES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
