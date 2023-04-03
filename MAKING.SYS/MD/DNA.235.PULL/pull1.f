!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Bonds  computes the hookean force between chosen atoms              *
!***********************************************************************
! pull the last bead with a spring
! the virtual end of the spring is moved by Vpdt
! this subroutine calculates the force in the spring that is connected to the first bead the static spring
! XYZspring1 are fixed from the begining of the simulation
      subroutine pull1 (E) !, Eprot, Ecarbo )
      include 'MD.com'
      integer I2, J2,  outE
      real r2, f, RT,r1
c     dimension RT(nBA)
c      real kspring, !Rspring,Vp
      real Fspring1 
      E = 0.0
c      kspring=0.12
c      kspring=60.
c     Rspring=5.
      
      
      dx = Xspring1 - X(Spring1Res)
      dy = Yspring1 - Y(Spring1Res)
      dz = Zspring1 - Z(Spring1Res)
      
      r2 = dx**2 + dy**2 + dz**2
      r1 = sqrt(r2)
      write(105,*)r1
!     energy calculation
c     E = E + kspring*(r1-Rspring)**2
      E = E + kspring*(r1-Rspring)**2
      
      
!     End energy calculation
      
      
!     f_over_r is the force over the magnitude of r so there is no need to resolve
!     the dx, dy and dz into unit vectors
      
      
c     f = kspring*Rspring - kspring*r1
      f = (kspring*Rspring)/r1 - kspring
      
                                ! now add the force  
      Fx(Spring1Res) = Fx(Spring1Res) - f * dx
      Fy(Spring1Res) = Fy(Spring1Res) - f * dy
      Fz(Spring1Res) = Fz(Spring1Res) - f * dz

      
      Fspring1=sqrt ((f*dx)**2 + (f*dy)**2 + (f*dz)**2)
      
!     To add a writing condition
      write(101,*)Fspring1
      
      
      E = E/2.0
      
      END
      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF BONDS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


