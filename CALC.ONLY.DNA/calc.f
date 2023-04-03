        implicit real*8(a-h,o-z)
        character*8 bead*3,residue*4,filename*1000
        real x,y,z,CMX,CMY,CMZ,CMX_prot,CMY_prot,CMZ_prot,d_min
	character*8 comm*8
        integer index,idnaend, idnastart,nt_in
        dimension x(1000000),y(1000000),z(1000000)
        character filename1*1000
        character*30 myfile1, myfile2


        idnaend=541	

        Call System("ls Traj_*.dat > input.txt")

        open(30, FILE='input.txt', ACCESS='sequential',status='old')

        do ifile=1,30

        read(30,*)filename1

        open(1,file=filename1,status="old")  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        read(1,*)no_chain
        do i=1,no_chain
        read(1,*)nc
        enddo
 
        do i=1,idnaend
        read(1,*)index,bead,residue
        enddo

        read(1,*)ntot
	write(*,*)ntot


	do mstep=1,99999


        read(1,*)kb
	write(*,*)kb
        do ij=1,idnaend

        read(1,'(a24)')line
        if (line(1:8)=='********')then
           continue
        else
           read(line,'(3f8.3)')x(ij),y(ij),z(ij)      
c        read(1,"(3f8.3)")x(ij),y(ij),z(ij)
c	write(*,*)x(ij),y(ij),z(ij)
        end if
        enddo
        
        ik=541
        il=2
        dx=x(ik)-x(il)
        dy=y(ik)-y(il)
        dz=z(ik)-z(il)
        r=sqrt(dx**2+dy**2+dz**2)

        write(myfile2,'(a,a,i3.3,a)')'time_dis','_',ifile,'.txt'

        open(14, file=myfile2, position="append")

        write(14,*)kb, r

        read(1,*)comm

	enddo
        close(1)
        close(14)
        enddo

        close(30)
c###############################################
        stop


        end
