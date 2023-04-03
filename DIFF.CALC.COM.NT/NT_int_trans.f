        implicit real*8(a-h,o-z)
        real dist,dis
        integer kb, idna
        dimension kb(100000),idna(100000),dist(100000),dis(100000)

        character*30 filename1, myfile1

        Call System("ls NT_ind*.txt > NT_ind1.input")

        open(30, FILE='NT_ind.input', ACCESS='sequential',status='old')

        do ifile=1,10

        read(30,*)filename1
        open(1,file=filename1,status="old")  

        write(myfile1,'(a,a,i3.3,a)')'int_trans','_',ifile,'.txt'
        open(13, file=myfile1, position="append")

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        do ij=1,49999

        read(1,*)kb(ij),idna(ij),dist(ij)
	write(*,*)kb(ij),idna(ij),dist(ij)

        enddo

        do i=1,49998
        dis(i)=abs((idna(i+1)-idna(i))/3)
        write(13,15)i,dis(i)
15      format(i10,f12.3)
        enddo

        enddo
c##############################################################################
        stop
        end
