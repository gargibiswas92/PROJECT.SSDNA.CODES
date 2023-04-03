!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* IntLD() : integrates over one time step for all particles           *
!***********************************************************************

      subroutine intLD(E,ET, Conts, outE, count )
! E is the total potential energy
! ET is the total potential and kinetic energy 
      include 'MD.com'
      integer outE, Qi, SC1, SC2, count 
      Real DumE, DumKE, Res,FR,SD,distance1,distance2,distance3
      dimension DumE(10), Res(AN), count(NC)
      real RsBy

      E = 0.0


      do i=1, 10
         DumE(i)=0.0
      end do
! these routines compute the forces
      call Fstart 

      call Bonds( DumE(1))

      call ANGL(DumE(2))
 
      call dihedral(DumE(3))
       

      if (useChirals .eq. 'YES') then
       call chiral(DumE(8))
      end if


      call contacts(DumE(4), Conts, count)
      
       call stacking(DumE(10))
c      call base_pairing(DumE(10))

      call NonContacts( DumE(5))

      if (useEllipsoidRepulsions .eq. 'YES') then
        call ellipsoidRepulsions(DumE(9))
      endif

! calculate the box force
! if there is a ssDNA it will not be confined by  the box
      if (confineInBox .eq. 'YES')then
         if (hasSSDNAdynamics .eq. 'YES') then
             call box(X,Y,Z,Fx,Fy,Fz ,firstDNABead-1,DumE(6),
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
c          write(87,*)esPairsNum
        if (useDebyeHuckel .eq. 'YES') then
          call debyehuckel(esFirstAtomIndex,esSecondAtomIndex,esCharge,
     Q  X,Y,Z,Fx,Fy,Fz,esPairsNum,esDistanceCutoff,DumE(7),
     Q  DebyeHuckelPotentials, DebyeHuckelForces , firstDNABead)
	else
          call coulomb(esFirstAtomIndex,esSecondAtomIndex,
     Q  esCharge,X,Y,Z,Fx,Fy,Fz,esPairsNum,DumE(7), deConstant)
	endif
      end if	

c      call enddis(distance1,distance2,distance3)
      call enddis(distance1)
    

      call extforce
 
      Call Findtemp(DumKE)
      ET = E + DumKE

      if(outE .eq. 0)then
      E = DumE(1) + DumE(2) + DumE(3) +  DumE(4) + 
     Q    DumE(5) + DumE(6) + DumE(7) + DumE(8) + DumE(9)+DumE(10)
        if(EnergyTerm .ne. 'NO')then
      write(70,"(10F15.6)") DumE(1),DumE(2),DumE(3),DumE(4),DumE(5),
     Q DumE(6),DumE(7),DumE(8), DumE(9),DumE(10)
	endif
      write(14,*)distance1
      endif


      do i=1, lastDynamicAtom
	SD = sqrt(2*ms(i)*gamma*T/tau)

	call gauss(0., SD, fr)
      Vx(i) = (Vx(i)*c_e + (Fx(i)+fr)*tau/ms(i))*c_i
        call gauss(0., SD, fr)
      Vy(i) = (Vy(i)*c_e + (Fy(i)+fr)*tau/ms(i))*c_i
        call gauss(0., SD, fr)
      Vz(i) = (Vz(i)*c_e + (Fz(i)+fr)*tau/ms(i))*c_i


      enddo


!      do i=1, AN

!      Fx(i) = X(i)
!      Fy(i) = Y(i)
!      Fz(i) = Z(i)

!      enddo

      do i=1, AN
      X(i) = X(i) + Vx(i)*tau
      Y(i) = Y(i) + Vy(i)*tau
      Z(i) = Z(i) + Vz(i)*tau
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
