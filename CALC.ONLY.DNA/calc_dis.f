        implicit real*8(a-h,o-z)
        character*8 bead*3,residue*4,filename*1000
        real x,y,z,x1,y1,z1,x2,y2,z2,r, rt
	character*8 comm*8
        integer index,idnaend,d, d1,q,p
        dimension x(1000000),y(1000000),z(1000000),x1(5),y1(5)
        dimension z1(5),x2(5),y2(5),z2(5),d(5),d1(5),rt(5)
        character filename1*1000
        character*30 myfile1, myfile2, line


        idnaend=541	

        Call System("ls Traj_*.dat > input.txt")

        open(30, FILE='input.txt', ACCESS='sequential',status='old')

        do ifile=1,30

        read(30,*)filename1

        open(1,file=filename1,status="old")  
        write(myfile2,'(a,a,i3.3,a)')'new_dis','_',ifile,'.txt'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        d(1)=4
        d(2)=7
        d(3)=10
        d(4)=13
        d(5)=16


        d1(1)=541
        d1(2)=538
        d1(3)=535
        d1(4)=532
        d1(5)=529

        read(1,*)no_chain
        do i=1,no_chain
        read(1,*)nc
        enddo
 
        do i=1,idnaend
        read(1,*)index,bead,residue
        enddo

        read(1,*)ntot
c	write(*,*)ntot


	do mstep=1,99999


        read(1,*)kb
c	write(*,*)kb
        do ij=1,idnaend
        read(1,'(a24)')line
        do p=1,5
        if(ij == d(p))then
        read(line,'(3f8.3)')x(ij),y(ij),z(ij)
            x1(p) = x(ij)
            y1(p) = y(ij)
            z1(p) = z(ij)
        else if(ij == d1(p))then
        read(line,'(3f8.3)')x(ij),y(ij),z(ij)
            x2(p) = x(ij)
            y2(p) = y(ij)
            z2(p) = z(ij)
        else
            continue
        end if
        enddo
        enddo

        do q=1,5
        dx=x1(q)-x2(q)
        dy=y1(q)-y2(q)
        dz=z1(q)-z2(q)
        r=sqrt(dx**2+dy**2+dz**2)
        rt(q)=r
        end do

        open(14, file=myfile2, position="append")

        write(14,'(i10, (5f8.3))')kb, (rt(q), q=1,5)

        read(1,*)comm

	enddo
        close(1)
        close(14)
        enddo

        close(30)
c###############################################
        stop


        end
