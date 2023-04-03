! This subroutine calculates the initial end-to-end vector and its unit vector in X,Y,Z directions
! later I calculate the XYZ components of the pulling velocity with the U vectors. 

      subroutine unitvec ! (Ve2eX,Ve2eY,Ve2eZ)
      include 'MD.com'

      real Ve2eX,Ve2eY,Ve2eZ

c      dx=X(1340)-X(1349)
c      dy=Y(1340)-Y(1349)
c      dz=Z(1340)-Z(1349)


      dx=X(Spring2Res)-X(Spring1Res)
      dy=Y(Spring2Res)-Y(Spring1Res)
      dz=Z(Spring2Res)-Z(Spring1Res)


      Re2e=sqrt(dx**2+dy**2+dz**2) ! end to end vector of the initial conformation
c      write(*,*)'Main Re2e=',Re2e

      Ve2eX=dx/Re2e ! unit vector in x direction
      Ve2eY=dy/Re2e
      Ve2eZ=dz/Re2e

      write(84,*)Ve2eX,Ve2eY,Ve2eZ,"ok"


      end
