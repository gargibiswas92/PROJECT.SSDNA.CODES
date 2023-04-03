!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* CONTACTS: computes the force on all atoms due to stack bases via a     *
!* 10-12 potential                                                     *
!***********************************************************************

      subroutine stacking(E,Edna,Edna_prot,En1,Arom_resi_Ene,idna)
      include 'MD.com'

      integer S1, S2,arom_resi,iprot,idna
      real  r2, rm2, rm10, f_over_r,Edna,Echeck,Edna_prot,En1,Arom_resi_Ene,
     & Emin
      dimension En1(Nmax),arom_resi(100),Arom_resi_Ene(Nmax),
     &idna(Nmax),Emin(Nmax)


      arom_resi(1)=413
      arom_resi(2)=356
      arom_resi(3)=432
      arom_resi(4)=1057
      arom_resi(5)=1171
      arom_resi(6)=1274
      arom_resi(7)=1178
      arom_resi(8)=1073
      arom_resi(9)=1039
      arom_resi(10)=890
      arom_resi(11)=836
      arom_resi(12)=904
      arom_resi(13)=655
      arom_resi(14)=545
      arom_resi(15)=593


      E = 0.0
      Edna = 0.0
      Edna_prot=0.0

      do i=1,15
      Emin(i)=10000.0
      enddo
      

	do i=1, NSTK
	
           S1 = IS(i)
           S2 = JS(i)
!        write(*,*) C1, C2


	dx = X(S1) - X(S2)

 	dy = Y(S1) - Y(S2)

	dz = Z(S1) - Z(S2)

	  r2 = dx**2 + dy**2 + dz**2

	      rm2 = 1.0/r2
              rm2 = rm2*sigma_stk(i)
	      rm10 = rm2**5

         Echeck=eps_stk(i)*rm2**5*(5*rm2-6)
         
         E=E+Echeck
        
        if((S1.ge.firstDNAbead).and.(S2.ge.firstDNAbead))then
        Edna=Edna+Echeck

        else
        Edna_prot=Edna_prot+Echeck



        if(S1.gt.S2)then
c       write(87,*)S1,S2,Echeck
        En1(S1)=En1(S1)+Echeck
        endif

        if(S1.lt.S2)then
c        write(87,*)S1,S2,Echeck,'ok'
        En1(S2)=En1(S2)+Echeck
        endif


c        write(87,*)S1,S2,Echeck    



        do iprot=1,15

        if(S1.eq.arom_resi(iprot).or.S2.eq.arom_resi(iprot))then


        if(S1.gt.S2.and.Echeck.lt.Emin(iprot))then
        idna(iprot)=S1
        Emin(iprot)=Echeck
        endif


        if(S1.lt.S2.and.Echeck.lt.Emin(iprot))then
        idna(iprot)=S2
        Emin(iprot)=Echeck
        endif

        
        goto 111

        endif

        enddo

111    continue  

        
c        do iprot=1,15
c        if(S1.eq.arom_resi(iprot).or.S2.eq.arom_resi(iprot))then
c        Arom_resi_Ene(iprot)=Arom_resi_Ene(iprot)+Echeck
c        endif
c        enddo

        endif



         f_over_r = eps_stk(i)*rm2**5*(rm2-1)*60/r2


! now add the acceleration 
	      FX(S1) = Fx(S1) + f_over_r * dx
	      FY(S1) = Fy(S1) + f_over_r * dy
	      FZ(S1) = Fz(S1) + f_over_r * dz

 	      Fx(S2) =  Fx(S2) - f_over_r * dx
	      Fy(S2) =  Fy(S2) - f_over_r * dy
	      Fz(S2) =  Fz(S2) - f_over_r * dz

              enddo

c         do n=1,15
c         write(87,*)arom_resi(n),idna(n)
c         enddo



      end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of Contacts^^^^^^^^^^^^^^^^^^^^^^^^^^^


