        implicit real*8(a-h,o-z)
        character*8 filename*1000
        real d
        integer flag, kb, count
        dimension d(50,1000000),flag(50,1000000), kb(1000000)
        character filename1*1000
        character*30 myfile1, myfile2, line
	

        open(1,file='all_dis.txt',status="old")  
        open(2,file='binary.txt',status="new")  
        open(3,file='sum_binary.txt',status="new")  
        

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do i=1,99999
        
        read(1, '(i10, 49f10.3)')kb(i), (d(j,i), j=1,49)
        end do

        do k=1,99999
        do l=1,49
        if (d(l,k) > 5)then
                flag(l,k) = 1
        else if (d(l,k) <= 5)then
                flag(l,k) = 0
        end if
        enddo
        enddo

        do ij=1,99999

        write(2, '(i10, 49i5)')ij, (flag(jk,ij), jk=1,49)   
        enddo


        do m=1,99999

        count = 0
        do n=1,49

        if (flag(n,m).eq.1)then
                count = count+1
        else
                continue
        end if
        end do
        write(3,'(i10,i10)')m, count
        end do

        close(1)
        close(2)
c        close(3)
c###############################################
        stop


        end
