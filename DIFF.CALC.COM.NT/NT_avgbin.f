        implicit real*8(a-h,o-z)
        real dist,sum,ssum,avgbin,idna1
        integer kb, idna
        dimension kb(100000),idna(100000),dist(100000),avgbin(100000),
     &sum(100000),ssum(100000),idna1(100000)

        character*30 filename1, myfile1

        Call System("ls NT_ind*.txt > NT_ind.input")

        open(30, FILE='NT_ind.input', ACCESS='sequential',status='old')

        do ifile=1,10

        read(30,*)filename1
        open(1,file=filename1,status="old")  

        write(myfile1,'(a,a,i3.3,a)')'avg_dis','_',ifile,'.txt'
        open(13, file=myfile1, position="append")

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        do ij=1,49999

        read(1,*)kb(ij),idna(ij),dist(ij)
c	write(*,*)kb(ij),idna(ij),dist(ij)

        enddo

        do i=1,49999
        avgbin(i)=0
        sum(i)=0
        enddo

        do j=1,49999
        idna1(j)=(idna(j)-1340)/3
        idna1(j)=idna1(j)+1
        enddo

        do il=1,49999-100
        bin=0.0
           do im=il,il+100
           bin=bin+im
           sum(il)=sum(il)+idna1(im)
           enddo
           ssum(il)=sum(il)/101.0
           avgbin(il)=bin/101.0

           write(13,20)avgbin(il),ssum(il)
20         format(f8.1, f12.3)
        enddo
        enddo
c###############################################
        stop
        end
