!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     This part is suitable for a system that contains a dynamic       *
!     single stranded DNA.                                             *
!     it will enable dynamics of the strand and will keep the center   *
!     of mass of DNA at the center of the box                          *
!***********************************************************************
      subroutine centerDNA

      include 'MD.com'

	real  Mass, CMX, CMY, CMZ


        if ((hasSSDNAdynamics .eq. 'YES') 
     Q      .and. (DnaCentered .eq. 'YES')) then


        CMX = 0.0
        CMY = 0.0
        CMZ = 0.0

        do i=firstDNABead , LastDNABead
c        do i=1,firstDNABead-1
           CMX = CMX + X(i)*ms(i)
           CMY = CMY + Y(i)*ms(i)
           CMZ = CMZ + Z(i)*ms(i)
           Mass=Mass+ms(i)
        end do

        CMX=CMX/Mass
        CMY=CMY/Mass
        CMZ=CMZ/Mass

c        write(87,*)CMX,CMY,CMZ,firstDNABead-1    


c        do i=firstDNABead , LastDNABead
c           X(i) = X(i)-CMX
c           Y(i) = Y(i)-CMY
c           Z(i) = Z(i)-CMZ
c        end do
       end if
       end
        
         
        
