!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* CONTACTS: computes the force on all atoms due to contacts via a     *
!* 10-12 potential                                                     *
!***********************************************************************

      subroutine contacts(E, Conts, count)
      include 'MD.com'

      integer C1, C2, ConfID,ContN, Q, SC1, SC2, Cf1, cf2, count(NC),
     Q summc
      real  r2, rm2, rm10, f_over_r, dsig, deps, Energy, Res, ResCon,
     Q s1, s2, ep1, ep2, r1, rc,r, summm , Edna
      dimension Energy(2), ContN(2), ResCon(AN), Res(AN)



      Q = 0
      E = 0.0
      Edna = 0.0
      Conts = 0

      do i=1, NC
         count(i) = 0
      end do

	do i=1, NC
	
           C1 = IC(i)
           C2 = JC(i)
!        write(*,*) C1, C2




	dx = X(C1) - X(C2)

 	dy = Y(C1) - Y(C2)

	dz = Z(C1) - Z(C2)

	  r2 = dx**2 + dy**2 + dz**2

	      rm2 = 1.0/r2
              rm2 = rm2*sigma(i)
	      rm10 = rm2**5

         E=E+epsC(i)*rm2**5*(5*rm2-6)

         f_over_r = epsC(i)*rm2**5*(rm2-1)*60/r2

         ! added 29/01/08
         if ((C1 .ge. firstDNABead) .and. 
     Q   (C2 .ge. firstDNABead))then
         Edna=Edna+epsC(i)*rm2**5*(5*rm2-6)
         endif
        ! end change

       ! we will only count native protein contacts
       ! and will ignore intra DNA contacts
       ! changed by amir (4-3-2008)
         if((r2 .le. sigma(i)*ConCut**2) .and.
     Q   (C1 .lt. firstDNABead) .and.
     Q   (C2 .lt. firstDNABead))then
            conts = conts+1
            count(conts) = i
         endif


! now add the acceleration 
	      FX(C1) = Fx(C1) + f_over_r * dx
	      FY(C1) = Fy(C1) + f_over_r * dy
	      FZ(C1) = Fz(C1) + f_over_r * dz

 	      Fx(C2) =  Fx(C2) - f_over_r * dx
	      Fy(C2) =  Fy(C2) - f_over_r * dy
	      Fz(C2) =  Fz(C2) - f_over_r * dz

              enddo


        E=E-Edna
c        E=Edna


      end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of Contacts^^^^^^^^^^^^^^^^^^^^^^^^^^^


