!***********************************************************************
! Common parameters and variables for DualGo.f                         *
!***********************************************************************

        IMPLICIT NONE

! These are settings
      	integer startt, N, stepstop, thermtime, ThermInt, neartime, WO,
     Q WOT
	real percent1, tau,T, Fext, Tstart,Tend,Tdiff,pchange, RsTau,gamma
	character (LEN=60) initval, output, finalpx1,  
     Q Conf, symtype


	character(LEN=20) PDB, atom
	parameter(PDB="(A6,I5,A4,TR1,A3,I6,3F8.3)")
	parameter(atom='ATOM  ')

	character(LEN=45)FMTX
        parameter(FMTX="(I5,I4,A4,A3,4F8.3)")

! LD parameters
	real c_e, c_i

! Pi is Pi
       real Pi
       parameter (pi = 3.141592653589793)
       integer Nmax, i, j, k, l,JN,MC,
     Q KEtimes

       integer writecount,index_dna,idna_elec
	character(LEN = 15) list,list1, list2
	parameter(list="(3F8.3,F6.2)")
	parameter(list1="(I5,4F8.3)")
	parameter(list2="(I5,3F8.3)")
       real  X, Y, Z, ms, Vx, Vy, Vz,
     Q Fx, Fy, Fz,Ene,EneStck,EneElec_resi,EneStck_resi

! This limits the size of the simulation.
       parameter (Nmax=20000)
       dimension X(Nmax), Y(Nmax), Z(Nmax), ms(Nmax),
     Q Vx(Nmax), Vy(Nmax), Vz(Nmax),
     Q Fx(Nmax), Fy(Nmax), Fz(Nmax),Ene(Nmax),EneStck(Nmax),
     &EneElec_resi(Nmax),EneStck_resi(Nmax),index_dna(Nmax),
     &idna_elec(Nmax)


       real temp, rand, temprand, xrandom, yrandom, zrandom,
     Q KEaveall, KEdev2


! These store information about the chain
	integer ResSeq,MDT,BeadIndex,ResidueIndex,ChainIndex
	character(LEN=4) AtType
	character(LEN=3) ResID
	dimension ResSeq(Nmax)
	dimension BeadIndex(Nmax)
	dimension ResidueIndex(Nmax)
	dimension ChainIndex(Nmax)
	dimension AtType(Nmax)
	dimension ResID(Nmax)
	integer CLmax
	parameter (clmax = 10)
	integer ChainLength
	dimension ChainLength(ClMax)


! XT, YT, and ZT are temp arrays to be used in one routine at a time
       real XT, YT, ZT
       dimension XT(Nmax), YT(Nmax), ZT(Nmax)
       integer AN, ANo

! these are dummy variables used in several routines
      real XIJ,YIJ,ZIJ,XKJ,YKJ, ZKJ, XKL,YKL,ZKL,DX,DY,
     + DZ, GX,GY,GZ,CT,CPHI,SPHI,Z1, Z2,FXI,FYI,FZI,
     + FXJ,FYJ,FZJ, FXK,FYK,FZK,FXL,FYL,FZL,DF, CT0, CT1,CT2
! end of dummy variables.

! These varaibles are used for indexing the atom types
	integer nAT, ATn, AT, ATbyType
	character (LEN = 2) ATindex
	parameter (nAT = 1)
	dimension ATindex(nAT), ATn(nAT), AT(Nmax),
     Q ATbyType(nAT, Nmax)
	parameter (ATindex =(/'CA'/))

! End of index varaibles

! Energy Variables
	real E, ET, P1, KE

! End of Energy Variables


! Variables used for contacts LJ
	integer NC, trip, tripi
	real PeakH
	parameter (peakH = 0.5)
	integer maxCon, CO, PairNum,
     Q noPairNum, N1, N2, IC, JC, Conts
	real sigma, epsC, ConCut, Shift
	parameter (ConCut = 1.2)
	parameter (maxCon =11*Nmax)
	dimension sigma(maxCon), epsC(Maxcon)
	dimension IC(maxcon), JC(maxcon), trip(3,maxcon)

! these are the non native contact pairs' variables
	integer NNoPairNum, NPairNum, INC,JNC, NNC, NNCmax,NNCt,NCset
	parameter (NNCmax = (Nmax-4)**2/2)
	real NCSigma, NNCsigma
	dimension  INC(NNCmax),JNC(NNCmax), 
     Q NCSigma(NNCmax), NNCsigma(NNCmax),NCset(NNCmax)
! End of non native

! these are the ellipsoid repulsions pairs' variables
	integer IEllipsoid, JEllipsoid, KEllipsoid,
     Q	        ellipsoidRepulsionsNum,ellipsoidRepulsionsMax,
     Q          currentEllipsoidRepulsions, 
     Q          currentEllipsoidRepulsionsNum
	parameter (ellipsoidRepulsionsMax = (Nmax-4)**2/2)
	real ellipsoidSigma, ellipsoidCoeff
	dimension  IEllipsoid(ellipsoidRepulsionsMax),
     Q             JEllipsoid(ellipsoidRepulsionsMax),
     Q             KEllipsoid(ellipsoidRepulsionsMax),
     Q		   ellipsoidSigma(ellipsoidRepulsionsMax),
     Q             ellipsoidCoeff(ellipsoidRepulsionsMax),
     Q             currentEllipsoidRepulsions(ellipsoidRepulsionsMax)
! End of ellipsoid repulsions variables

! these are the electrostatics variables
	integer esAtomsNum, maxEsAtoms, esFirstAtomIndex, esSecondAtomIndex,
     Q          esPairsNum, esMinBeadDistance, esMinNeighbor
	real esCharge, esEnergyCutoff, ionicStrength,ionicRadius,
     Q       solventDensity, esDistanceCutoff, DebyeHuckelPotentials,
     Q       DebyeHuckelForces, esCutoffMax
        parameter (esCutoffMax = 40000.0)
	parameter (maxEsAtoms = Nmax)
	dimension esFirstAtomIndex(maxEsAtoms**2)
	dimension esSecondAtomIndex(maxEsAtoms**2)
        dimension esCharge(maxEsAtoms**2)
        dimension DebyeHuckelPotentials(8000000)
        dimension DebyeHuckelForces(8000000)
! End of electrostatics  	


!Variables for stacked bases
	Integer NSTKmax
        parameter (NSTKmax=Nmax*3)
        integer IS, JS, nSTK
        dimension IS(NSTKmax), JS(NSTKmax)
	real sigma_STK, GAMMA_STK, delta_STK, eps_STK
	dimension sigma_STK(NSTKmax),Gamma_STK(NSTKmax),
     &delta_STK(NSTKmax), eps_STK(NSTKmax)
!End of varaible for stackes bases


!Variables for stacked bases
        Integer NBPSmax
        parameter (NBPSmax=Nmax*2)
        integer nBPS
        real sigmaB, epsB
        dimension sigmaB(NBPSmax), epsB(NBPSmax)
!End of varaible for stackes bases




! Variables for bonded pairs
	integer NBmax
	parameter (NBmax=Nmax*2)
	integer Ib1, Ib2, nBA, intB, nPROT
	dimension Ib1(NBmax), Ib2(NBmax)
	real Rb, RbT, bK, RBC
	dimension  Rb(NBmax), RbT(NBmax),
     Q  bK(NBmax), RBC(NBmax)
! End of varaible for bonded pairs

! Variables for single strand dynamics   (added 1/1/2008) by amir
        integer firstDNABead, LastDNABead

! Variables for bond angles
	integer NTmax
	parameter (NTmax =Nmax*2)
	real ANTT, TK, intT, ANTC
	dimension ANTT(Ntmax), TK(Ntmax), ANTC(Ntmax)
	integer IT, JT, KT, nTA
	dimension IT(Ntmax), JT(Ntmax), KT(Ntmax)
! End of varaibles for Bond Angles

! Variables for Phi angles
	integer npmax
	parameter (npmax =NMAX*3)
	integer nPA, IP, JP, KP, LP, intP
	real APT, PK, AP0, AP1, Z10, Z20, Z11, Z22, Z12, Dums, 
     Q DFLIM, DF1, DF0, DR1, DR2, DR3, DR4, DR5, DR6, DRX, DRY, DRZ
	dimension APT(Npmax), IP(Npmax), JP(Npmax), KP(Npmax),
     Q LP(Npmax), PK(Npmax)

	real GAMC1, GAMC3, GAMS1, GAMS3, DihAng
	dimension GAMC1(Npmax), GAMC3(Npmax), GAMS1(Npmax),
     Q DihAng(Npmax), GAMS3(Npmax)

      real GMUL, TM24,TM06,tenm3, zero,one,two,four,six,twelve,ellipsoidRepulsionsNum ftem,
     Q S, COSNP, SINNP, COSNP3, SINNP3, DC1, DC2, DC3,
     Q DC4, DC5, DC6, EPW
	real cosarray, sinarray
	integer refinephi
	parameter (refinephi=10000)
	dimension cosarray(3*refinephi), sinarray(3*refinephi)

! end of varaibles for Phi Angles.
	
! Variables for Chiral angles
	integer nChiralMax
	parameter (nChiralMax =NMAX*3)
	integer nChirals, Ichiral, Jchiral, Kchiral, Lchiral
        real chiralCoeff, chiralValue
	dimension chiralValue(Npmax), Ichiral(Npmax), Jchiral(Npmax),
     Q             Kchiral(Npmax), Lchiral(Npmax), chiralCoeff(Npmax)
! end of variables for chiral angles

	character(LEN=200) Trajectory, EnergyTot,EnergyTerm,ContactFile,
     Q TemperatureVtime, ThreeBodyFile, TrajDist, ENDdisFile
     &,F1Vtime,F2Vtime
! conditional flags for additional execution settings 
        character(LEN=25) hasStaticAtoms, confineInBox,
     Q                    useElectrostatics, useDebyeHuckel, useChirals,
     Q                    useEllipsoidRepulsions, esCutoffType,
     Q                    hasSSDNAdynamics, DnaCentered,tempAnnealing


! addtional parameters, used only if conditional flag is set 

        integer lastDynamicAtom,Spring1Res ,Spring2Res
        real boxMin,boxMax, boxCoeff
	dimension boxMin(3),boxMax(3)
        real deConstant, screeningFactor, saltCoefficient

! Variables
        character (LEN=5) PullWithSpring
        real Re2e, Ue2eX, Ue2eY,Ue2eZ
        real Xspring1,Yspring1,Zspring1
        real Xspring2,Yspring2,Zspring2
        real Kspring,Rspring,Vp

       COMMON /char/ initval, output, finalpx1, 
     Q conf, AtType, ResID,Trajectory, EnergyTot,EnergyTerm,ContactFile,
     Q TemperatureVtime, ThreeBodyFile, symtype, hasStaticAtoms,
     Q confineInBox, useElectrostatics, useDebyeHuckel, useChirals,
     Q useEllipsoidRepulsions, TrajDist, esCutoffType, hasSSDNAdynamics,
     Q DnaCentered,F1Vtime,F2Vtime,PullWithSpring

       COMMON /real/ X, Y, Z, ms, Vx,
     Q Vy, Vz, Fx, Fy, Fz,rand, temprand, xrandom, yrandom, zrandom,
     Q KEaveall, KEdev2, tau, rstau,gamma,c_i, c_e,
     Q T, Fext,pchange, temp, Rb, RbT, RBC, bK, ANTC, ANTT, APT, 
     Q GAMS1, GAMS3, GAMC1, GAMC3, PK, sigma,sigmaB, shift,
     Q epsC, epsB, NCSigma,NNCsigma, esCharge, DihAng,
     Q xt, yt, zt, cosarray, sinarray,
     Q boxMin,boxMax,boxCoeff,
     Q deConstant, screeningFactor, saltCoefficient,
     Q chiralCoeff, chiralValue, ellipsoidSigma, ellipsoidCoeff,
     Q esEnergyCutoff, ionicStrength, ionicRadius, solventDensity,
     Q esDistanceCutoff, DebyeHuckelPotentials, DebyeHuckelForces,
     & sigma_STK, GAMMA_STK, delta_STK,eps_STK,Tstart,Tend,Tdiff, 
     & Re2e, Ue2eX, Ue2eY,Ue2eZ,
     Q Xspring1,Yspring1,Zspring1,
     Q Xspring2,Yspring2,Zspring2,Kspring,Rspring,Vp

 
     
       COMMON /int/ KEtimes, writecount,WO,WOT,NNC,NNCt,NCset,
     Q startt, N, stepstop, thermtime, ThermInt, neartime,
     Q AN, ANo, ATn, AT, ATbyType, Ib1, Ib2, trip, tripi,
     Q nBA, IT, JT, KT, nTA, TK, intT, intB, intP, IP, JP, KP,
     Q  LP, nPA, IC, JC,IS,JS, N1, N2, PairNum, NoPairNum, INC,JNC,
     Q NPairNum, NNoPairNum, ResSeq, chainlength, MDT, NC, nPROT,
     Q lastDynamicAtom, esAtomsNum,esFirstAtomIndex,esSecondAtomIndex,
     Q esPairsNum,nChirals, Ichiral, Jchiral, Kchiral, Lchiral,
     Q IEllipsoid, JEllipsoid, KEllipsoid,ellipsoidRepulsionsNum,
     Q currentEllipsoidRepulsions, currentEllipsoidRepulsionsNum,
     Q esMinBeadDistance, BeadIndex,ResidueIndex,ChainIndex,
     Q firstDNABead,LastDNABead, nSTK,Spring1Res ,Spring2Res
