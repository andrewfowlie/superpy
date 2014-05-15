      real*8 function findValW(name)
      implicit none
      character*(*) name
      findValW=0
      end

      implicit none

      character*100  fInput 
      character*100  fOutput
      character*20   rdBuff
      real *8 v,vcsnngg,vcsnngz,vcsgg,vcsgz
      integer nArgs,err ,slharead,mode
      
      
      

      nArgs=iargc()

      if(nArgs.eq.3) then
           mode=1
      else 
          if(nArgs.eq.2) then
               mode=0
          else
            return
          endif
      endif       
 
      call getarg(1,fInput)
      call getarg(2,fOutput)   
          
      if(slharead(fInput,0).ne.0) return
                              
      call ModelConstIni(mode,err)
	  
      if(err.ne.0) return
      
      v=0.01
      vcsgg=vcsnngg(v)
      vcsgz=vcsnngz(v)      
      OPEN(UNIT=78,FILE=fOutput,STATUS='UNKNOWN')
      write(78,fmt='("BLOCK lGamma  # AZ and AA cross sections")')
 
      write(78,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(78,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(78)
      return 
      end 
