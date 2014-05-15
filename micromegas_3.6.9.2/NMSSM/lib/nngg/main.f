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
      real*8 Wh(5), slhaWidth
          
      external slhaWidth
              
      if(slharead('spectr',0).ne.0.or.slharead('decay',1).ne.0) goto 3
      nArgs=iargc()
      if(nArgs.eq.0) then
           mode=0
      else
           mode=1       
      endif

      Wh(1)=slhaWidth(25) 
      Wh(2)=slhaWidth(35)
      Wh(3)=slhaWidth(45)
      Wh(4)=slhaWidth(36)
      Wh(5)=slhaWidth(46)
                                                     
      call ModelConstIni(mode,Wh,err)
      if(err.ne.0) goto 3
      v=0.01
      vcsgg=vcsnngg(v)
      vcsgz=vcsnngz(v)      

      OPEN(UNIT=78,FILE='nngg.out',STATUS='UNKNOWN')
      write(78,fmt='("BLOCK lGamma  # AZ and AA cross sections")')
      write(78,fmt='(A6, I2,A20)')      '  0   ', mode,'#  GI flag'
      write(78,fmt='(A6, 1PE10.4,A20)') '  1   ', vcsgz,'# ~o1,~o1->A,Z [pb]'
      write(78,fmt='(A6, 1PE10.4,A20)') '  2   ', vcsgg,'# ~o1,~o1->A,A [pb]'
      close(78)
      return 
3     close(78)
      return   
      end 
