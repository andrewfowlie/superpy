
      implicit none

      character*100  fInput 
      character*100  fOutput
      real *8 v,vcsnngg,vcsnngz,vcsgg,vcsgz,ModelConstIni,zero
      integer nArgs,err ,slharead
            
      nArgs=iargc()
      if(nArgs.lt.2) stop 
     >'2 or 3 arguments  expected for lGamma.exe: in/out files, GI key'       
      call getarg(1,fInput) 
      call getarg(2,fOutput)
      
      err=slharead(fInput,1)
      if(err.ne.0)  stop 'Problem in reading  input file'
      zero=ModelConstIni(nArgs -1 )
      v=0.001
      vcsgg=vcsnngg(v)
      vcsgz=vcsnngz(v)      
 
      OPEN(UNIT=78,FILE=fOutput,STATUS='UNKNOWN')   
      write(78,fmt='("BLOCK lGamma  # AZ and AA cross sections")')
      write(78,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(78,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(78)
      end 
