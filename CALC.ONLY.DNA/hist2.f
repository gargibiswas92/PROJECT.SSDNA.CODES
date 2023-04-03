        implicit real*8(a-h,o-z)
        character*8 filename*1000
        real d, bin
        integer time, hist, t, count1
        dimension time(1000000),hist(1000000),t(1000000),bin(1000000)
        character filename1*1000

        open(3,file='sum_binary.txt',status="old")  
        open(4,file='binary_bin.txt',status="new")  
        

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        do m=1,99999

        read(3,'(i10,i10)')time(m),hist(m)
c        write(*,'(i10,i10)')time(m),hist(m)
        end do


        do j=1,98999
        count1 = 0
        sum=0
         do k=j,j+999
            sum=sum+hist(k)
            count1 = count1+ 1
         end do
         bin(j)=sum/count1
         t(j)=j+(count1/2)
         write(4,'(i10,f8.2)')t(j), bin(j)
         end do
        close(3)
        close(4)
c###############################################
        stop


        end
