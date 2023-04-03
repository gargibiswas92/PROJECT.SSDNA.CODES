
!* contains subrutines  start, stophere

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* start creates reads the run settings from settings.dat, reads in the*
!* random variables for the random subroutine and creates output files.*
!***********************************************************************

      subroutine start

      include 'MD.com'

      open(60, file = 'settings.dat', access ='sequential',
     Q status = 'old')

      read(60,*)
      read(60,*)
      read(60,*) startt ! indicates the format of the initial conditions
!      read(60,*) AN ! number of atoms in simulation
      read(60,*) Conf ! Conf is the input file generated by Read.f
      read(60,*) initval, finalpx1 ! initial coordinates,restart file
      read(60,*) output ! output information about simulation
      read(60,*) Trajectory, WOT ! Trajectory is the trajectory file name

	if(Trajectory .ne. 'NO')then
          open(6, FILE= Trajectory, status = 'unknown')
          read(60,*) TrajDist ! whether to store distances
	  if(TrajDist .ne. 'NO')then
            open(93, FILE= TrajDist, access = 'SEQUENTIAL',
     Q      status = 'unknown')
          end if
	endif

      read(60,*) EnergyTot, EnergyTerm ! These are the energy files
	if(EnergyTot .ne. 'NO')then
        open(99,file=EnergyTot,status='unknown',access='sequential')
	endif
        if(EnergyTerm .ne. 'NO')then
        open(70,file=EnergyTerm,status='unknown',access='sequential')
        endif

      read(60,*) ContactFile, ThreeBodyFile ! this stores the 
                                            ! contacts and 3 body contacts

        if(ContactFile .ne. 'NO')then
        open(12,file=ContactFile,status='unknown',access='sequential')
        endif
        if(ThreeBodyFile .ne. 'NO')then
        open(13,file=ThreeBodyFile,status='unknown',access='sequential')
        endif

        read(60,*) ENDdisFile ! this is the temperature file
        if(ENDdisFile .ne. 'NO')then
        open(14,file=ENDdisFile,
     Q  status='unknown',access='sequential')
        endif
        
        read(60,*) TemperatureVtime ! this is the temperature file

        if(TemperatureVtime .ne. 'NO')then
        open(58,file=TemperatureVtime,
     Q  status='unknown',access='sequential')
        endif

      read(60,*) symtype ! MD for Molecular Dynamics, LD for Langevin D.
      if(symtype .eq. 'LD')then
      read(60,*) gamma ! gamma is the drag coefficient
      endif
      read(60,*) stepstop, WO ! number of integration steps, write out
                                   ! measurements every WO steps, write out
                                   
      read(60,"(2F10.5)") tau, RsTau ! time step, coupling constant for Berendsen 
                            ! thermostat

       read(60,*) tempAnnealing
      if(tempAnnealing .ne. 'NO')then
        read(60,*) Tstart ! temperature to start the simulation, in reducedunits
        read(60,*) Tend ! temperature to end the simulation, in reduced units
        Tdiff = Tstart - Tend
        T = Tstart
      else
        read(60,*) T ! temperature of the simulation
      endif

       read(60,*) PullWithSpring 
       if(PullWithSpring .eq. 'YES')then
       read(60,*) Kspring !pulling spring constant
       read (60,*)Rspring ! Equlibriuum spring length
       read(60,*)Vp !the velocity of the spring in units of distance (A) pertime step (tau)
       read(60,*)Spring1Res,Spring2Res
       read (60,*)F1Vtime,F2Vtime ! the files of forces with time
       open (101,file=F1Vtime)
       open (102,file=F2Vtime)
        else
        
        read(60,*) Fext ! Force, in reduced units

       endif

      read(60,*) hasStaticAtoms ! whether to keep some atoms static
      if(hasStaticAtoms .eq. 'YES')then
      read(60,*) lastDynamicAtom !number of protein atoms (the rest are static atoms) !ADDED BY OHAD
      endif

      read(60,*) confineInBox ! whether to confine the molecule in a box
      if(confineInBox .eq. 'YES')then
       read(60,*) boxMin(1),boxMin(2),boxMin(3) !minimum values for box
       read(60,*) boxMax(1),boxMax(2),boxMax(3) !maximum values for box
       read(60,*) boxCoeff !box force cofficient
      endif


      read(60,*) useElectrostatics !whether to apply coulombic interactions
      if(useElectrostatics .eq. 'YES')then
        read(60,*) deConstant !dielectric constant in epsilon 0 units
        read(60,*) esMinBeadDistance !dielectric constant in epsilon 0 units
        read(60,*) esCutoffType ! minimum energy to apply electrostatic force
        if(esCutoffType .eq. 'ENERGY')then
          read(60,*) esEnergyCutoff ! minimum energy to apply electrostatic force
	endif
        if(esCutoffType .eq. 'DISTANCE')then
          read(60,*) esDistanceCutoff ! minimum energy to apply electrostatic force
	endif
        read(60,*) useDebyeHuckel ! whether to use debye huckel screening factor
        if(useDebyeHuckel .eq. 'YES')then
          read(60,*) ionicStrength ! ionic strength of the solution
          read(60,*) ionicRadius ! the average radius of the ions in the solution (for example 1.4 for NaCl)
          read(60,*) solventDensity ! the specific density of the solvent (for example 1 for water)
          ! this supposedly transforms simulation temperature to real temperatures
        
	call initES
        call dhenergytable(screeningFactor,saltCoefficient,
     Q deConstant,esDistanceCutoff,DebyeHuckelPotentials,
     Q DebyeHuckelForces)
        endif
      endif
       read(60,*) useChirals
       read(60,*) useEllipsoidRepulsions
       read(60,*) hasSSDNAdynamics
       if(hasSSDNAdynamics .eq. 'YES')then   !wheter to simulate dynamics of single strand DNA
         read(60,*) DnaCentered              !wheter to keep the DNA strand at the center of mass
         if(DnaCentered .eq. 'YES')then      
          read(60,*) firstDNABead ! the first phosphate in the strand
          read(60,*) LastDNABead !  the last base in the strand 
         endif
         
       endif
 

      open(63, file = 'random.dat', access = 
     Q 'SEQUENTIAL',  status = 'UNKNOWN')
      read(63,*) xrandom,yrandom,zrandom
      close (63, status = 'KEEP')


	if(symtype  .eq. 'LD')then
	c_e = 1.0 - gamma*tau/2.0
	c_i = 1.0/(1+gamma*tau/2.0)
	endif

	writecount = 0
      end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^end of start^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* stophere closes all files that are open.  It also stores            *
!* the new xrandom, yrandom and zrandom in file random.dat for use of  *
!* the random subroutine, next time this is run.                       *
!***********************************************************************

      subroutine stophere

      include 'MD.com'

! Close the files that holds the temperature

      close(71)

! this writes the new xrandom, yrandom and zrandom 
      open(63, file = 'random.dat', access = 
     Q 'SEQUENTIAL',  status = 'UNKNOWN')
      write(63,*) xrandom,yrandom,zrandom
      close (63, status = 'KEEP')

      return

      end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of stophere^^^^^^^^^^^^^^^^^^^^^^^^

 
      subroutine initES()
      include 'MD.com'
      
       real screeningFactorConstant,realTemperature, esCutoff

       real esEnergyError, esMinDistance, esMaxDistance, 
     Q     esCurrentDistance, esEnergy 

      parameter (screeningFactorConstant = 2529.11892) 

       if(useDebyeHuckel .eq. 'YES')then
! this supposedly transforms simulation temperature to real temperatures
         realTemperature = 295*T
c          realTemperature = 350-((1-T)*100)

! calculates the screening factor kappa
         screeningFactor = SQRT((screeningFactorConstant*
     Q   ionicStrength*solventDensity)
     Q   /(deConstant*realTemperature))
	 saltCoefficient = (exp(screeningFactor*ionicRadius)/
     Q   (1+screeningFactor*ionicRadius))
       endif
       if (esCutoffType .eq. 'ENERGY')then
         esMinDistance = 16.0
         esMaxDistance = 100000000.0

           if(useDebyeHuckel .eq. 'YES')then
	     call debyehuckelfactor(esMinDistance,deConstant,
     Q            screeningFactor,saltCoefficient,esEnergy)
           else
	     call coulombfactor(esMinDistance,deConstant,esEnergy)
           endif
         esCutoff = esEnergy * esEnergyCutoff
       
         esEnergyError = esCutoff / 100.0

         do i=1,5000
           esCurrentDistance = (esMaxDistance+esMinDistance)/2
           if(useDebyeHuckel .eq. 'YES')then
	     call debyehuckelfactor(esCurrentDistance,deConstant,
     Q            screeningFactor,saltCoefficient, esEnergy)
           else
	     call coulombfactor(esCurrentDistance,deConstant,
     Q                          esEnergy)
           endif
           	
	   ! esCutoff distance found
           if (abs(esCutoff-esEnergy) .le. esEnergyError)then
             esDistanceCutoff =  esCurrentDistance
             if (esDistanceCutoff .ge. esCutoffMax) then
               write(*,*) 
     Q         'distance Cutoff value exceeds maximum allowed'
               call abort
             end if             
             return 
           endif
           if(esEnergy .le. esCutoff)then
             esMaxDistance = esCurrentDistance
	   else
             esMinDistance = esCurrentDistance
           endif
         end do

         write(*,*) 'es cutoff distance not found'
         call abort
       end if

       if (esCutoffType .eq. 'DISTANCE')then
         if(esDistanceCutoff**2 .gt. esCutoffMax)then
           write(*,*) 
     Q     'distance Cutoff value exceeds maximum allowed'
         call abort
	 end if
       end if
      end