!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Bonds  computes the hookean force between chosen atoms              *
!***********************************************************************
! pull the last bead with a spring
! the virtual end of the spring is moved by Vpdt
! the moving spring


      subroutine pull2 (E)      !, Eprot, Ecarbo )
      include 'MD.com'
      integer I2, J2,  outE
      real r2, f, RT,r1
c     dimension RT(nBA)
c     real kspring, Vp ,Fspring2!Rspring,Vp
      real Fspring2
       E = 0.0
c       kspring=0.12
c       kspring=60.
c       Rspring=5.
c       Vp=0.005 ! in units of distance (A) per time step (tau)
c       write(*,*) Xspring2
        Xspring2= Xspring2 + ( Ue2eX*Vp)! without *tau because the Vp is in units of distance per tau
	Yspring2= Yspring2 + ( Ue2eY*Vp)
	Zspring2= Zspring2 + ( Ue2eZ*Vp)


        write(85,*)Ue2eX,Ue2eY,Ue2eZ,Xspring2,Yspring2,Zspring2



c        write(107,"(3f8.3)")Xspring2,Yspring2,Zspring2
c        write(*,*) Xspring2 , Ue2eX*Vp*tau
	dx = Xspring2 - X(Spring2Res)
	dy = Yspring2 - Y(Spring2Res)
	dz = Zspring2 - Z(Spring2Res)
        
        r2 = dx**2 + dy**2 + dz**2
        r1 = sqrt(r2)
c        write(106,*)r1
! energy calculation
c        E = E + kspring*(r1-Rspring)**2
         E = E + kspring*(r1-Rspring)**2



! End energy calculation



! f_over_r is the force over the magnitude of r so there is no need to resolve
! the dx, dy and dz into unit vectors


c            f = kspring*Rspring - kspring*r1
            f = (kspring*Rspring)/r1 - kspring

            ! now add the force  
	      Fx(Spring2Res) = Fx(Spring2Res) - f * dx
	      Fy(Spring2Res) = Fy(Spring2Res) - f * dy
	      Fz(Spring2Res) = Fz(Spring2Res) - f * dz

              Fspring2=sqrt ((f*dx)**2 + (f*dy)**2 + (f*dz)**2)

!To add a writing condition


                dx = X(Spring1Res) - X(Spring2Res)
                dy = Y(Spring1Res) - Y(Spring2Res)
                dz = Z(Spring1Res) - Z(Spring2Res)

               r2 = dx**2 + dy**2 + dz**2
               r1 = sqrt(r2)


               write(102,"(2F9.3)")Fspring2,r1
c              write(106,*)r1
              

             E = E/2.0

      END
      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF BONDS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


