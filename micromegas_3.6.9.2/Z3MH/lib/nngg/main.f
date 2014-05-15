      real*8 function findValW(name)
      implicit none
      character*(*) name
      findValW=0
      end

      implicit none

      character*100  fInput 
      character*100  fOutput
      character*20 rdBuff
      real *8 v,vcsnngg,vcsnngz,vcsgg,vcsgz
      integer nArgs,err ,slharead,mode
      
      OPEN(UNIT=78,FILE='nngg.out',STATUS='UNKNOWN')
      write(78,fmt='("BLOCK lGamma  # AZ and AA cross sections")')
              
      
      err=slhaRead('nngg.in',0) 
      if(err.ne.0)  stop 'Problem in reading  input file'                                                    
      call ModelConstIni
	  
c      if(err.ne.0) goto 3
      v=0.001
      vcsgg=vcsnngg(v)
c      vcsgz=vcsnngz(v)      
 
c      write(78,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(78,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(78)
      return 
3     close(78)
      return   
      end 
