! This subroutine calculates the initial end-to-end vector and its unit vector in X,Y,Z directions
! later I calculate the XYZ components of the pulling velocity with the U vectors. 

      subroutine endtoend ! (Ue2eX,Ue2eY,Ue2eZ)
      include 'MD.com'

c      real Re2e, Ue2eX, Ue2eY, Ue2eZ


c       dx=X(Spring2Res)-X(Spring1Res)
c       dy=Y(Spring2Res)-Y(Spring1Res)
c       dz=Z(Spring2Res)-Z(Spring1Res)

c--------------1nd end pulling direction based on protein residues-------------- 
c      dx=X(356)-X(1146)
c      dy=Y(356)-Y(1146)
c      dz=Z(356)-Z(1146)
c-----------------------------------

c--------------2nd end pulling direction based on protein residues-------------- 
      dx=X(675)-X(801)
      dy=Y(675)-Y(801)
      dz=Z(675)-Z(801)
c-----------------------------------

      Re2e=sqrt(dx**2+dy**2+dz**2) ! end to end vector of the initial conformation
c     write(*,*)'Main Re2e=',Re2e

      Ue2eX=dx/Re2e ! unit vector in x direction
      Ue2eY=dy/Re2e
      Ue2eZ=dz/Re2e

c      write(84,*)Ue2eX,Ue2eY,Ue2eZ,"ok"

      Xspring1= X(Spring1Res)- Ue2eX*Rspring
      Yspring1= Y(Spring1Res)- Ue2eY*Rspring
      Zspring1= Z(Spring1Res)- Ue2eZ*Rspring

c      write(84,*)'XYZ spring1', Xspring1, Yspring1,Zspring1

      Xspring2= X(Spring2Res)+ Ue2eX*Rspring
      Yspring2= Y(Spring2Res)+ Ue2eY*Rspring
      Zspring2= Z(Spring2Res)+ Ue2eZ*Rspring
c      write(84,*)'XYZ spring2', Xspring2, Yspring2,Zspring2
c      write(84,*)Rspring
      end
