C This file contains some auxiliary Fortran routines which can not
C be simulated in C. They are 
C    i)  file input/output routines
C    ii) Fortran functions which return complex numbers. 
C        ( Because authors are not sure that there is an unique standard 
C         for writing corresponding bridge in C) 
                                                         

C    This auxillary routine opens file 'fname' and writes its content
C    into   preliminary opened  by Fortran tools file channel Nch 
 
      subroutine fortreread(Nch,fname)
      integer Nch,Nrd,i
      character*(*) fname
      character*200 buff
      character*10 fff
      logical open

      
      Nrd=12
4     Nrd=Nrd+1 
      INQUIRE(UNIT=Nrd,OPENED=open)
      if(open) goto 4       

      open(Nrd,FILE=fname,STATUS='OLD')
9     continue
        read(Nrd,fmt='(a200)',end=11) buff
        do i=200,1,-1
          if(buff(i:i) .ne. ' ') goto 10
        enddo
10      continue
        if(i.eq.0) i=1
        if(i.lt.10) then
           write(fff,fmt='(A2,I1,A1)')  '(a',i,')'
        else 
          if(i.lt.100) then 
            write(fff,fmt='(A2,I2,A1)') '(a',i,')' 
          else 
            write(fff,fmt='(A2,I3,A1)') '(a',i,')'
          endif
        endif    
        write(Nch,fmt=fff) buff
      goto 9
11    close(Nrd)
      end  

C  Read line routine for Fortran version of readSLHA  
C 

      integer function fortranreadline(N,buff)
      integer N
      character *(*) buff
      character *20 forma
C      write(forma,'(A2,I3,A1)')'(A',LEN(BUFF),')'
      write(forma,*) '(A',LEN(BUFF),')'
      read(N,fmt=forma,end=10,err=10) buff 
      fortranreadline=0 
      return
10    fortranreadline=-1
      return
      end       

C Fortran routines  with complex return value.
       
      Complex*16 function cMixMatrix(id,i,j)
      integer id,i,j
      real*8  re,im,reMixMatrix,imMixMatrix
            re=reMixMatrix(id,i,j)
            im=imMixMatrix(id,i,j)
            cMixMatrix=COMPLEX(re,im)
      end      

      Complex*16 function cMixMatrixU(id,i,j)
      integer id,i,j 
      real*8  re,im,reMixMatrixU,imMixMatrixU
            re=reMixMatrixU(id,i,j)
            im=imMixMatrixU(id,i,j)
            cMixMatrixU=COMPLEX(re,im)
      end

