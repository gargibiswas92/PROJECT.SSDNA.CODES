        implicit real*8(a-h,o-z)
        integer stacktable
        dimension stacktable(4,4)     

        do i=1,4
        do j=i,4
        if(j.eq.i)then
        if (i.eq.1)stacktable(i,j)=3
        if (i.eq.2)stacktable(i,j)=1
        if (i.eq.3)stacktable(i,j)=3
        if (i.eq.4)stacktable(i,j)=4
        else if(j.eq.i+1)then       
        if (i.eq.1)stacktable(i,j)=2
        if (i.eq.2)stacktable(i,j)=1
        if (i.eq.3)stacktable(i,j)=3
        else if(j.eq.i+2)then        
        if (i.eq.1)stacktable(i,j)=2
        if (i.eq.2)stacktable(i,j)=2
        else if(j.eq.i+3)then     
        if (i.eq.1)stacktable(i,j)=3
        else
        endif

        enddo
        enddo

        do k=2,4
        do l=1,k-1
        stacktable(k,l)=stacktable(l,k)
        enddo
        enddo

        do i=1,4
        write(*,*)stacktable(i,1),stacktable(i,2),stacktable(i,3),
     &stacktable(i,4)
        enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        open(1,file='sequence.dat',ACCESS= 'sequential',status='old')
        read(25, "(A4,I7,A6,A4,A2,I3,4x,3F8.3,2x,F4.2,F6.2,11x,A3)")
     &BeadIndex,Index, AtType, Residuename,chainid,noresidue,x,y,z,
     &s1,s2,last

       
       if(Residuename.eq.'A')then
       inedxstack=1
       elseif(Residuename.eq.'T')then
       inedxstack=2
       elseif(Residuename.eq.'C')then
       inedxstack=3
       elseif(Residuename.eq.'G')then
       inedxstack=4
       else
       endif
 





        END

c        open(1,file='sequence.dat',ACCESS= 'sequential',status='old')
c        read(25, "(A4,I7,A6,A4,A2,I3,4x,3F8.3,2x,F4.2,F6.2,11x,A3)")
c     &BeadIndex,Index, AtType, Residuename,chainid,noresidue,x,y,z,
c     &s1,s2,last
        
