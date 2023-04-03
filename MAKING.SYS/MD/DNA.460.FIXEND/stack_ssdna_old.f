!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Bonds  computes the hookean force between chosen atoms              *
!***********************************************************************

      subroutine stacking(E)
      include 'MD.com'
      integer I2, J2
      real r2,r1, st, stf,Edna,Estack,f


      E = 0.0
      Edna = 0.0  !this will calculate the bonds energy in the DNA

	do 1 i=1, nSTK
           I2 = IS1(i)
           J2 = IS2(i)

	dx = X(I2) - X(J2)
	dy = Y(I2) - Y(J2)
	dz = Z(I2) - Z(J2)
       
	  r2 = dx**2 + dy**2 + dz**2
          r1 = sqrt(r2)

! energy calculation

             st=(r1/sigma_STK(i))
             st=20*(st-GAMMA_STK(i))
             st=st+1
             stf=1/exp(st)
             Estack = delta_STK(i)*stf 
             E = E - Estack
         



! End energy calculation

! f_over_r is the force over the magnitude of r so there is no need to resolve
! the dx, dy and dz into unit vectors

! the index i indicates the interaction between particle i and i+1

            f = -Estack*(20/sigma_STK(i))   
            f=f/r1

            ! now add the force  
	      Fx(I2) = Fx(I2) + f * dx
	      Fy(I2) = Fy(I2) + f * dy
	      Fz(I2) = Fz(I2) + f * dz
c! the negative sign is due to the computation of dx, dy and dz
	      Fx(J2) = Fx(J2) - f * dx
	      Fy(J2) = Fy(J2) - f * dy
	      Fz(J2) = Fz(J2) - f * dz

1         continue
cc             write(87,*)E,Edna
c             E = E/2.0
c             Edna = Edna/2.0
cc            E = E - Edna     !changed 29/01/08
c             E = Edna     !changed 29/01/08
        write(87,*)E

      END
      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF BONDS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

