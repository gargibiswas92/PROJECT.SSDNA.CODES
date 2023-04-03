        implicit real*8(a-h,o-z)
        CHARACTER bead*4 ,residue*4, pfresidue*3, ck*3
        character filename1*50, filename2*50
        dimension ck(10000)
        

        write(*,*)"file reads beads, sequence file, total_no_residue,
     &old_repulsive pair, new_added_pair"


        read(*,*)filename1, filename2, num_of_total_resi, irep, nadded
        write(*,*)filename1, filename2, num_of_total_resi, irep, nadded


c        num_of_total_resi=1414

        d1=4.410
        d2=5.664
        d3=7.076
        ep=1.000

c        irep=990869

        open(27, FILE='added_repulsion.dat', ACCESS= 'sequential',
     &  status='unknown')
        open(25, FILE=filename2, ACCESS= 'sequential',
     &  status='old')
        do j=1,num_of_total_resi
        read(25,70)index1,index2, bead, residue
70      format(I4,I4,A3,A4)
        ck(index1)=bead
c        write(*,*)index1,index2,bead,residue
        enddo
        close(25)


        open(26, FILE=filename1, ACCESS= 'sequential',
     &  status='old')
        do i=1,nadded
        read(26,*)i1,i2
        write(*,*)i1,i2,ck(i1),ck(i2) 
 
        if(ck(i1).eq." CB" .and. ck(i2) .eq. " CB")then
        irep=irep+1
        write(27,80)irep,i1,i2,d1,ep
        endif
        if((ck(i1).eq." CA" .and. ck(i2) .eq. " CB").or.(ck(i1).eq." CB"
     &.and. ck(i2) .eq. " CA"))then
        irep=irep+1
        write(27,80)irep,i1,i2,d2,ep
        endif     
   
        if(ck(i1).eq." CA" .and. ck(i2) .eq. " CA")then
        irep=irep+1
        write(27,80)irep,i1,i2,d3,ep
        endif

80      format(I8,I5,I5,f10.3,f10.3)
        enddo









        stop
        end
