      subroutine base_pairing(E)
      include 'MD.com'

      integer C1, C2, INUM, IBP
      real r_hb, r2, rm2, rm10, f_over_r, Edna
  
      E = 0.0
      Edna = 0.0
      INUM=0.0
      r_hb=4.0

	do i=firstDNAbead+2,lastDNAbead-9, 3
        IBP=0
        do j=i+9,lastDNAbead,3
        INUM=INUM+1
       
        C1=I
        C2=J
 
        if(IBP.eq.0)then
	dx = X(C1) - X(C2)
 	dy = Y(C1) - Y(C2)
	dz = Z(C1) - Z(C2)

         r2 = dx**2 + dy**2 + dz**2

         if(r2.le.r_hb)then
c        write(87,*)INUM,i,j,sigmaB(INUM),epsB(INUM),r2,IBP

	      rm2 = 1.0/r2
              rm2 = rm2*sigmaB(INUM)
	      rm10 = rm2**5

         E=E+epsB(INUM)*rm2**5*(5*rm2-6)

         f_over_r = epsB(INUM)*rm2**5*(rm2-1)*60/r2

! now add the acceleration 
	      FX(C1) = Fx(C1) + f_over_r * dx
	      FY(C1) = Fy(C1) + f_over_r * dy
	      FZ(C1) = Fz(C1) + f_over_r * dz

 	      Fx(C2) =  Fx(C2) - f_over_r * dx
	      Fy(C2) =  Fy(C2) - f_over_r * dy
	      Fz(C2) =  Fz(C2) - f_over_r * dz

             IBP=1

             else
        
             IBP=0

             endif

             endif
               
c            if(IBP.eq.1)goto 200    
                

             enddo
c200          continue
      
             enddo


             E=Edna


             end



