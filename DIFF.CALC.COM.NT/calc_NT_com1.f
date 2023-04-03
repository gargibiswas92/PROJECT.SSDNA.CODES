        implicit real*8(a-h,o-z)
        character*8 bead*3,residue*4,filename*1000
        real x,y,z,CMX,CMY,CMZ,CMX_prot,CMY_prot,CMZ_prot,d_min
	character*8 comm*8
        integer index,idnaend, idnastart,nt_in
        dimension x(100000),y(100000),z(100000)
        character filename1*1000
        character*30 myfile1, myfile2, line


        idnastart= 1340
        idnaend=2719	
        iprot=1339	
	Prot_bead_totat=1339

        Call System("ls Traj_*.dat > input.txt")

        open(30, FILE='input.txt', ACCESS='sequential',status='old')

        do ifile=1,10

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


	do mstep=1,49999


        read(1,*)kb
	write(*,*)kb
        do ij=1,idnaend
        read(1,'(a24)')line
        if (line(1:8)=='********')then
        continue
        else
        read(line,'(3f8.3)')x(ij),y(ij),z(ij)
c            x(ij)=line(1:8)
c            y(ij)=line(9:16)
c            z(ij)=line(17:24)
c        read(1,"(3f8.3)")x(ij),y(ij),z(ij)
	write(*,*)x(ij),y(ij),z(ij)
        end if
        enddo
        

	CMX=0.0
	CMY=0.0
	CMZ=0.0
       
        bead_no = 0

	do ik=1, iprot
	CMX=CMX+x(ik)
	CMY=CMY+y(ik)
	CMZ=CMZ+z(ik)
        bead_no = bead_no + 1
	enddo

	CMX_prot=CMX/bead_no
	CMY_prot=CMY/bead_no
	CMZ_prot=CMZ/bead_no



	nt_in=0
	dmin= 1000000.0

        do j=idnastart,idnaend,3
        dx=CMX_prot-x(j)
        dy=CMY_prot-y(j)
        dz=CMZ_prot-z(j)
        r=sqrt(dx**2+dy**2+dz**2)

	if(r.lt.dmin)then 
	dmin=r
	nt_in=j	
	endif
       
        enddo
       

        write(myfile1,'(a,a,i3.3,a)')'prot_COM','_',ifile,'.txt'
        write(myfile2,'(a,a,i3.3,a)')'NT_ind','_',ifile,'.txt'

        open(13, file=myfile1, position="append")
        open(14, file=myfile2, position="append")


	write(13,*)Kb,CMX_prot, CMY_prot, CMZ_prot

        write(14,*)kb, nt_in, dmin

        read(1,*)comm

	enddo
        close(1)
        close(13)
        close(14)
        enddo

        close(30)
c###############################################
        stop


        end
