        implicit real*8(a-h,o-z)
        character*8 filename*1000
        real kb, d
        dimension d(30,1000000),kb(1000000)
        character filename1*1000
        character*30 myfile1, myfile2, line
	

        Call System("ls dis_*.txt > input_dis.txt")

        open(30, FILE='input_dis.txt', ACCESS='sequential',status='old')

        do ifile=1,20

        read(30,*)filename1

        open(1,file=filename1,status="old")  
        

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        

	do ij=1,99999


        read(1,*)kb(ij), d(ifile,ij)
	
        enddo
        enddo

        open(14, file='all_dis.txt', position="append")

        do i=1,99999

        write(14, '(i10, 20f10.3)')i, (d(j,i), j=1,20)

        enddo

        close(1)
        close(14)
        close(30)
c###############################################
        stop


        end
