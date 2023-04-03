!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Bonds  computes the hookean force between chosen atoms              *
!***********************************************************************

      subroutine stacking(E)
      include 'MD.com'
      integer I2, J2,  outE
      real r2,r1, Edna


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

         write(*,*)pause
c! energy calculation
c             E = E + bk(i)*(r1-Rb(i))**2
c
c
c      ! changed 29/1/08 
c      if ((I2 .ge. firstDNABead) .and. (J2 .ge. firstDNABead))then
c        Edna = Edna + bk(i)*(r1-Rb(i))**2
c      endif
c      ! end changed 29/1/08
c
c	if( bk(i)*(r1-Rb(i))**2 .gt. 50.0)then
c
c!	write(46,*) 
c!	write(46,*) bk(i)*(r1-Rb(i))**2
c!	write(46,*) I2,J2
c!	write(46,*) r1, Rb(i), bk(i)
c	endif
c! End energy calculation
c
c! f_over_r is the force over the magnitude of r so there is no need to resolve
c! the dx, dy and dz into unit vectors
c
c! the index i indicates the interaction between particle i and i+1
c
c            f = RBC(i)/r1 - bK(i)
c
c            ! now add the force  
c	      Fx(I2) = Fx(I2) + f * dx
c	      Fy(I2) = Fy(I2) + f * dy
c	      Fz(I2) = Fz(I2) + f * dz
c! the negative sign is due to the computation of dx, dy and dz
c	      Fx(J2) = Fx(J2) - f * dx
c	      Fy(J2) = Fy(J2) - f * dy
c	      Fz(J2) = Fz(J2) - f * dz
c
1         continue
cc             write(87,*)E,Edna
c             E = E/2.0
c             Edna = Edna/2.0
cc            E = E - Edna     !changed 29/01/08
c             E = Edna     !changed 29/01/08

      END
      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF BONDS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* BONDSP does the same as ANGLP but for the bond lengths              *
!* computes RBC for each bond (this computation is done only once)     *
!* to compute the force: f = RBC(i)/r1 - bK(i)
!* where r1 is the delta from the optimal length 
!***********************************************************************

      subroutine bondsp
      include 'MD.com'

      do i=1, nBA

      RBC(i) = Rb(i)*bK(i)

      end do

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF BONDSP^^^^^^^^^^^^^^^^^^^^^^^^^^^^
