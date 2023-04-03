!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* IntLD() : integrates over one time step for all particles           *
!***********************************************************************

      subroutine intLD(E,ET, Conts, outE, count )
! E is the total potential energy
! ET is the total potential and kinetic energy 
      include 'MD.com'
      integer outE, Qi, SC1, SC2, count,esPairs 
      Real DumE, DumKE, Res,FR,SD,distance1,distance2,distance3
      dimension DumE(17), Res(AN), count(NC)
      real RsBy,sum1,sum2

      E = 0.0


      do i=1, 17
         DumE(i)=0.0
      end do


       do i=1,17
       index_dna(i)=0
       enddo


c###########################################
        do i=1, lastDynamicAtom
        Ene(i)=0
        EneStck(i)=0
        EneStck_resi(i)=0
        EneElec_resi(i)=0
        enddo


        do i=1,lastDynamicAtom
        idna_elec(i)=0
        enddo


c###########################################
        sum1=0
        sum2=0


! these routines compute the forces
      call Fstart 

      call Bonds( DumE(1))

      call ANGL(DumE(2))
 
      call dihedral(DumE(3))
       

      if (useChirals .eq. 'YES') then
       call chiral(DumE(8))
      end if


      call contacts(DumE(4), Conts, count)
      
      call stacking(DumE(10),DumE(11),DumE(12),EneStck,EneStck_resi,index_dna)



c     call base_pairing(DumE(10))

      call NonContacts( DumE(5))

      if (useEllipsoidRepulsions .eq. 'YES') then
        call ellipsoidRepulsions(DumE(9))
      endif

! calculate the box force
! if there is a ssDNA it will not be confined by  the box
      if (confineInBox .eq. 'YES')then
         if (hasSSDNAdynamics .eq. 'YES') then
             call box(X,Y,Z,Fx,Fy,Fz ,LastDNABead,DumE(6),
     Q       boxMin,boxMax,boxCoeff)
         endif
         else
        call box(X,Y,Z,Fx,Fy,Fz ,lastDynamicAtom,DumE(6),boxMin,
     Q       boxMax,boxCoeff)
      endif

! calculate electrostatic force
      if (useElectrostatics .eq. 'YES') then
      ! E1,2 - residue indexes
      ! Q102 - charge of bond
      ! X,Y,Z - residue position
      !  FX, FY, FZ - residue applied force
      ! NE - number of electrostatic contacts
        if (useDebyeHuckel .eq. 'YES') then
          call debyehuckel(esFirstAtomIndex,esSecondAtomIndex,esCharge,
     Q  X,Y,Z,Fx,Fy,Fz,Ene,EneElec_resi,idna_elec,esPairsNum,
     &esDistanceCutoff,DumE(7),DumE(13), DumE(14),DumE(15),
     Q  DebyeHuckelPotentials, DebyeHuckelForces , firstDNABead)
	else

          call coulomb(esFirstAtomIndex,esSecondAtomIndex,
     Q  esCharge,X,Y,Z,Fx,Fy,Fz,esPairsNum,DumE(7), deConstant)
	endif
      end if	


c      call enddis(distance1,distance2,distance3)
c      call enddis(distance1)

      if(PullWithSpring .eq. 'YES')then 
c     call Pull1( DumE(16))   
      call unitvec 
      call Pull2( DumE(17))
      else
      call extforce
      endif
 
      Call Findtemp(DumKE)
      ET = E + DumKE

      if(outE .eq. 0)then
      E = DumE(1) + DumE(2) + DumE(3) +  DumE(4) + 
     Q    DumE(5) + DumE(6) + DumE(7) + DumE(8) + DumE(9)+DumE(10)
        if(EnergyTerm .ne. 'NO')then
      write(70,"(15F15.6)") DumE(1),DumE(2),DumE(3),DumE(4),DumE(5),
     Q DumE(6),DumE(8), DumE(9),DumE(10),DumE(11),DumE(12),DumE(7),
     QDumE(13),DumE(14),DumE(15), DumE(16), DumE(17)
        
c#########################################################################
        do i= firstDNABead,lastDynamicAtom-2,3
        sum1=sum1+Ene(i)
        enddo
        write(12,"(67F15.6)")(Ene(i),i=firstDNABead,
     &lastDynamicAtom-2,3),sum1,DumE(15)


        do i= firstDNABead+2,lastDynamicAtom,3
        sum2=sum2+EneStck(i)
        enddo
        write(13,"(67F15.6)")(EneStck(i),i=firstDNABead+2,
     &lastDynamicAtom,3),sum2,DumE(12)



c      write(14,"(29F15.6)")(EneStck_resi(i),i=1,15),(EneElec_resi(j),j=1,14)
      write(14,"(29I10)")(index_dna(i),i=1,15),(idna_elec(j),j=1,14)

c###########################################################################	
        endif
c      write(14,*)distance1
      endif


      do i=1, lastDynamicAtom
	SD = sqrt(2*ms(i)*gamma*T/tau)

	call gauss(0., SD, fr)
      Vx(i) = (Vx(i)*c_e + (Fx(i)+fr)*tau/ms(i))*c_i
        call gauss(0., SD, fr)
      Vy(i) = (Vy(i)*c_e + (Fy(i)+fr)*tau/ms(i))*c_i
        call gauss(0., SD, fr)
      Vz(i) = (Vz(i)*c_e + (Fz(i)+fr)*tau/ms(i))*c_i

c        if(i.eq.lastDynamicAtom)then
c	Vx(i)=0.0
c	Vy(i)=0.0
c	Vz(i)=0.0
c        endif

        enddo


!      do i=1, AN

!      Fx(i) = X(i)
!      Fy(i) = Y(i)
!      Fz(i) = Z(i)

!      enddo

      do i=1, AN
      if(i.eq.firstDNABead.or.i.eq.2717)goto 742
      X(i) = X(i) + Vx(i)*tau
      Y(i) = Y(i) + Vy(i)*tau
      Z(i) = Z(i) + Vz(i)*tau
742   continue
      enddo

!      do i=1, AN
!      X(i) = X(i) + tau*Vx(i)
!      Y(i) = Y(i) + tau*Vy(i)
!      Z(i) = Z(i) + tau*Vz(i)
!      enddo

!      do i=1, AN

!      Vx(i) = (X(i)

!      enddo


      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of IntSymp^^^^^^^^^^^^^^^^^^^^^^^^^^^^
