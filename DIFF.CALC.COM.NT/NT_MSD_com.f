        implicit real*8(a-h,o-z)
        real msd_nt,tt,dis
        integer num
        dimension tt(100000),dis(100000)
        character*30 filename1, myfile1

        Call System("ls avg_dis*txt > avg_dis.input")
        open(30, FILE='avg_dis.input',ACCESS='sequential',status='old')

        do ifile=1,10

        read(30,*)filename1

        open(1,file=filename1,status="old")

        write(myfile1,'(a,a,i3.3,a)')'MSD_NT','_',ifile,'.txt'
        open(13, file=myfile1, position="append")  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       

        do ij=1,49898

        read(1,*)tt(ij),dis(ij)
	write(*,*)tt(ij),dis(ij)

        enddo

         
        do ns=1,1000
        msd_nt=0.0
        num=0
           do il=1,49898-ns
              im=il+ns
              num= num+1
              dis_nt=(dis(im)-dis(il))**2
              msd_nt=msd_nt+dis_nt

           enddo
        write(13,20)ns,msd_nt/num
20      format(i8, f12.3)

         enddo 
         enddo       
c###############################################
        stop
        close(1)
        close(13)
        end
