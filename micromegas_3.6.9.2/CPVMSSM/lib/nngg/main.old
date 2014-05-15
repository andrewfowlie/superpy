
      implicit none

      character*100  fInput 
      character*100  fOutput
      real *8 v,vcsnngg,vcsnngz,vcsgg,vcsgz
      integer nArgs,err ,slharead
      
      nArgs=iargc()
      if(nArgs.ne.2) stop 'Two arguments expected: in/output file names'
       
      call getarg(1,fInput) 
      call getarg(2,fOutput)
      
      err=slharead(fInput,1)
      if(err.ne.0)  stop 'Problem in reading  input file'
      v=0.01
      vcsgg=vcsnngg(v)
      vcsgz=vcsnngz(v)      
 
      OPEN(UNIT=78,FILE=fOutput,STATUS='UNKNOWN')   
      write(78,*) 'BLOCK lGamma  # AZ and AA cross sections'
      write(78,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(78,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(78)
      end 
