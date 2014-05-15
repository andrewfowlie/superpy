*     LanHEP output produced at Wed Apr 10 22:10:25 2013
*     Model named 'MSSM&CPV'

#include "DgetNEWv0.inc"
#include "CgetNEWv0.inc"

      subroutine ModelDefaults
      implicit none
c      real*8 slhaVal0,slhaVal1,slhaVal2
c      integer slhaValExists0,slhaValExists1,slhaValExists2
#include "SLHAplus.fh"
      integer err
	  
	  real*8 Q 

#include "model.h"
	  Q=100

      EE = 0.31223D0
      SW = 0.481D0
      s12 = 0.221D0
      s23 = 0.041D0
      s13 = 0.0035D0
      GG = 1.238D0
      MZ = 91.1884D0
      Me = 0.000511D0
      Mm = 0.1057D0
      Mu = 0.046D0
      Md = 0.046D0
      Mc = 1.42D0
      Ms = 0.2D0
      wt = 1.7524D0
      wZ = 2.4944D0
      wW = 2.08895D0
      Ml = 1.777D0
      Mt = slhaVal1('SMINPUTS',Q,6)
      Mb = 2.7D0
	  
        tb= slhaVal1('EXTPAR',Q,25) 
	  
      MG3 = slhaVal1('EXTPAR',Q,3)
      MG1r=slhaVal1('EXTPAR',Q,1)
      MG2r=slhaVal1('EXTPAR',Q,2)
      MG1i=slhaVal1('IMEXTPAR',Q,1)
      MG2i=slhaVal1('IMEXTPAR',Q,2)
      MG3i=slhaVal1('IMEXTPAR',Q,3)
      mur=slhaVal1('EXTPAR', Q,23)
      mui=slhaVal1('IMEXTPAR', Q,23)
	  
      Alr= slhaVal1('EXTPAR', Q,13)
	  Ali= slhaVal1('IMEXTPAR', Q,13)
	  Abr= slhaVal1('EXTPAR', Q,12)
	  Atr= slhaVal1('EXTPAR', Q,11)
	  Abi= slhaVal1('IMEXTPAR', Q,12)
	  Ati= slhaVal1('IMEXTPAR', Q,11)
	  
	  
	  aMG1=sqrt(MG1r**2+MG1i**2)
	  aMG2=sqrt(MG2r**2+MG2i**2)
	  amu=sqrt(mur**2+mui**2)
c	  aAt=sqrt(Atr**2+Ati**2)
c	  aAb=sqrt(Abr**2+Abi**2)
c	  aAl=sqrt(Alr**2+Ali**2)
	  
      fiMG1=atan2(MG1i,MG1r)
      fiMG2=atan2(MG2i,MG2r)
      fimu=atan2(mui,mur)
c      fiAt=atan2(Ati,Atr)
c      fiAb=atan2(Abi,Abr)
c      fiAl=atan2(Ali,Alr)
	  
      Ml2=slhaVal1('EXTPAR',Q,32)
	  Ml1=slhaVal1('EXTPAR',Q,31)
      Ml3=slhaVal1('EXTPAR',Q,33)
      Mr2=slhaVal1('EXTPAR',Q,35)
	  Mr1=slhaVal1('EXTPAR',Q,34)
      Mr3=slhaVal1('EXTPAR',Q,36)
	  
      Mq2=slhaVal1('EXTPAR',Q,42) 
	  Mq1=slhaVal1('EXTPAR',Q,41) 
      Mq3=slhaVal1('EXTPAR',Q,43)
      Mu2=slhaVal1('EXTPAR',Q,45)  
	  Mu1=slhaVal1('EXTPAR',Q,44)  
      Mu3=slhaVal1('EXTPAR',Q,46)      
      Md2=slhaVal1('EXTPAR',Q,48)
	  Md1=slhaVal1('EXTPAR',Q,47)
      Md3=slhaVal1('EXTPAR',Q,49)
	  
      Ae = 0D0
      Am = 0D0
      Au = 0D0
      Ad = 0D0
      As = 0D0
      Ac = 0D0
	  
	  
      MHc = slhaVal1('MASS', Q,37)
	  
      wHc = 0D0
      wNE1 = 0D0
      wNE2 = 0D0
      wNE3 = 0D0
      wNE4 = 0D0
      wC1 = 0D0
      wC2 = 0D0
      wSG = 0D0
      wh1 = 0D0
      wh2 = 0D0
      wh3 = 0D0
      Zl33i = 0D0
      wSe1 = 0D0
      wSe2 = 0D0
      wSm1 = 0D0
      wSm2 = 0D0
      wSl1 = 0D0
      wSl2 = 0D0
      wSne = 0D0
      wSnm = 0D0
      wSnl = 0D0
      wSu1 = 0D0
      wSu2 = 0D0
      wSd1 = 0D0
      wSd2 = 0D0
      wSc1 = 0D0
      wSc2 = 0D0
      wSs1 = 0D0
      wSs2 = 0D0
      wSt1 = 0D0
      wSt2 = 0D0
      wSb1 = 0D0
      wSb2 = 0D0
	  
	  hisW(1)=slhaWidth(25)
	  hisW(2)=slhaWidth(35)
	  hisW(3)=slhaWidth(36)
	  wh1=hisW(1)
	  wh2=hisW(2)
	  wh3=hisW(3)
	  
      end


      subroutine ModelConstIni(mode, err)
      implicit none
      integer fail

#include "model.h"

c      double precision initQCD,McRun,MbRun,MtRun,alphaQCD,McEff,MbEff,
c     - MtEff,slhaVal,System,aPrintF,openAppend,slhaRead,
c     - rDiagonal3,rDiagonalA,MassArray,MixMatrix,MixMatrixU,HggF,HggS,
c     - HggV,HggA,neutcp1,neutcp2,chacp1,chacp2

      double precision neutcp1,neutcp2,chacp1,chacp2, Q
	  integer inidi, htk, mode, err, i, j

#include "SLHAplus.fh"

      call ModelDefaults
	  Q=100
	  err=0
      CW = sqrt(1D0-SW**2)
      S2W = 2D0*SW*CW
      C2W = CW**2-SW**2
      MW = MZ*CW
      c12 = sqrt(1D0-s12**2)
      c23 = sqrt(1D0-s23**2)
      c13 = sqrt(1D0-s13**2)
      sb = tb/sqrt(1D0+tb**2)
      cb = sqrt(1D0-sb**2)
      t2b = 2D0*tb/(1D0-tb**2)
      c2b = cb**2-sb**2
      s2b = 4D0*sb*cb/2D0
      ls1 = Me*Ae/Sqrt2*EE/SW/MW/cb
      ls2 = Mm*Am/Sqrt2*EE/SW/MW/cb
      MSne = sqrt(c2b*MW**2/2D0+CW**2*Ml1**2)/CW
      MSnm = sqrt(c2b*MW**2/2D0+CW**2*Ml2**2)/CW
      MSnl = sqrt(c2b*MW**2/2D0+CW**2*Ml3**2)/CW
      MSeLL = Ml1**2+Me**2-MW**2/CW**2*(1D0/2D0-SW**2)*c2b
      MSeRR = Mr1**2+Me**2-MW**2/CW**2*SW**2*c2b
      MSeLR = -Me*amu*tb
      MSe1 = sqrt((MSeLL+MSeRR-sqrt((MSeLL-MSeRR)**2+4D0
     &      *MSeLR**2))/2D0)
      MSe2 = sqrt((MSeLL+MSeRR+sqrt((MSeLL-MSeRR)**2+4D0
     &      *MSeLR**2))/2D0)
      MSeth = atan2(-2D0*MSeLR, -MSeLL+MSeRR)/2D0
      Zl11 = cos(MSeth)
      Zl14 = -sin(MSeth)
      Zl41 = -Zl14
      Zl44 = Zl11
      MSmuLL = Ml2**2+Mm**2-MW**2/CW**2*(1D0/2D0-SW**2)*c2b
      MSmuLR = MW*SW*Sqrt2*cb*ls2/EE-Mm*amu*tb
      MSmuRR = Mr2**2+Mm**2-MW**2/CW**2*SW**2*c2b
      MSm1 = sqrt((MSmuLL+MSmuRR-sqrt((MSmuLL-MSmuRR)**2
     &      +4D0*MSmuLR**2))/2D0)
      MSm2 = sqrt((MSmuLL+MSmuRR+sqrt((MSmuLL-MSmuRR)**2
     &      +4D0*MSmuLR**2))/2D0)
      MSmuth = atan2(-2D0*MSmuLR, -MSmuLL+MSmuRR)/2D0
      Zl22 = cos(MSmuth)
      Zl25 = -sin(MSmuth)
      Zl52 = -Zl25
      Zl55 = Zl22
      MtopS = Mt
      us1 = -Mu*Au*EE/SW/MW/Sqrt2/sb
      ds1 = Md*Ad*EE/SW/MW/Sqrt2/cb
      ds2 = Ms*As*EE/SW/MW/Sqrt2/cb
      us2 = -Mc*Ac*EE/SW/MW/Sqrt2/sb
      MSuLL = Mq1**2+Mu**2+MW**2/CW**2*(1D0/2D0-2D0/3D0*SW
     &      **2)*c2b
      MSuLR = -MW*SW*Sqrt2*sb*us1/EE-Mu*amu/tb
      MSuRR = Mu1**2+Mu**2+2D0/3D0*MW**2/CW**2*SW**2*c2b
      MSu1 = sqrt((MSuLL+MSuRR-sqrt((MSuLL-MSuRR)**2+4D0
     &      *MSuLR**2))/2D0)
      MSu2 = sqrt((MSuLL+MSuRR+sqrt((MSuLL-MSuRR)**2+4D0
     &      *MSuLR**2))/2D0)
      MSuth = atan2(-2D0*MSuLR, -MSuLL+MSuRR)/2D0
      Zu11 = cos(MSuth)
      Zu14 = -sin(MSuth)
      Zu41 = -Zu14
      Zu44 = Zu11
      MSdLL = Mq1**2+Md**2-MW**2/CW**2*(1D0/2D0-1D0/3D0*SW
     &      **2)*c2b
      MSdLR = MW*SW*Sqrt2*cb*ds1/EE-Md*amu*tb
      MSdRR = Md1**2+Md**2-1D0/3D0*MW**2/CW**2*SW**2*c2b
      MSd1 = sqrt((MSdLL+MSdRR-sqrt((MSdLL-MSdRR)**2+4D0
     &      *MSdLR**2))/2D0)
      MSd2 = sqrt((MSdLL+MSdRR+sqrt((MSdLL-MSdRR)**2+4D0
     &      *MSdLR**2))/2D0)
      MSdth = atan2(-2D0*MSdLR, -MSdLL+MSdRR)/2D0
      Zd11 = cos(MSdth)
      Zd14 = -sin(MSdth)
      Zd41 = -Zd14
      Zd44 = Zd11
      MScLL = Mq2**2+Mc**2+MW**2/CW**2*(1D0/2D0-2D0/3D0*SW
     &      **2)*c2b
      MScLR = -MW*SW*Sqrt2*sb*us2/EE-Mc*amu/tb
      MScRR = Mu2**2+Mc**2+2D0/3D0*MW**2/CW**2*SW**2*c2b
      MSc1 = sqrt((MScLL+MScRR-sqrt((MScLL-MScRR)**2+4D0
     &      *MScLR**2))/2D0)
      MSc2 = sqrt((MScLL+MScRR+sqrt((MScLL-MScRR)**2+4D0
     &      *MScLR**2))/2D0)
      MScth = atan2(-2D0*MScLR, -MScLL+MScRR)/2D0
      Zu22 = cos(MScth)
      Zu25 = -sin(MScth)
      Zu52 = -Zu25
      Zu55 = Zu22
      MSsLL = Mq2**2+Ms**2-MW**2/CW**2*(1D0/2D0-1D0/3D0*SW
     &      **2)*c2b
      MSsLR = MW*SW*Sqrt2*cb*ds2/EE-Ms*amu*tb
      MSsRR = Md2**2+Ms**2-1D0/3D0*MW**2/CW**2*SW**2*c2b
      MSs1 = sqrt((MSsLL+MSsRR-sqrt((MSsLL-MSsRR)**2+4D0
     &      *MSsLR**2))/2D0)
      MSs2 = sqrt((MSsLL+MSsRR+sqrt((MSsLL-MSsRR)**2+4D0
     &      *MSsLR**2))/2D0)
      MSsth = atan2(-2D0*MSsLR, -MSsLL+MSsRR)/2D0
      Zd22 = cos(MSsth)
      Zd25 = -sin(MSsth)
      Zd52 = -Zd25
      Zd55 = Zd22
c      MG1r = aMG1*cos(fiMG1*pi/180D0)
c      MG1i = aMG1*sin(fiMG1*pi/180D0)
c      MG2r = aMG2*cos(fiMG2*pi/180D0)
c      MG2i = aMG2*sin(fiMG2*pi/180D0)
c      mur = amu*cos(fimu/180D0*pi)
c      mui = amu*sin(fimu/180D0*pi)
      fMG1 = fiMG1*pi/180D0
      fMG2 = fiMG2*pi/180D0
      fmu = fimu*pi/180D0
	  
	  if(mode.eq.1) then
      neutct = neutcp1(MZ, MW, SW, aMG1, aMG2, amu, tb, fiMG1
     &      , fiMG2, fimu)
      Zn11r = neutcp2(neutct, 1D0, 1D0, 1D0)
      Zn11i = -neutcp2(neutct, 2D0, 1D0, 1D0)
      Zn12r = neutcp2(neutct, 1D0, 1D0, 2D0)
      Zn12i = -neutcp2(neutct, 2D0, 1D0, 2D0)
      Zn13r = neutcp2(neutct, 1D0, 1D0, 3D0)
      Zn13i = -neutcp2(neutct, 2D0, 1D0, 3D0)
      Zn14r = neutcp2(neutct, 1D0, 1D0, 4D0)
      Zn14i = -neutcp2(neutct, 2D0, 1D0, 4D0)
      Zn21r = neutcp2(neutct, 1D0, 2D0, 1D0)
      Zn21i = -neutcp2(neutct, 2D0, 2D0, 1D0)
      Zn22r = neutcp2(neutct, 1D0, 2D0, 2D0)
      Zn22i = -neutcp2(neutct, 2D0, 2D0, 2D0)
      Zn23r = neutcp2(neutct, 1D0, 2D0, 3D0)
      Zn23i = -neutcp2(neutct, 2D0, 2D0, 3D0)
      Zn24r = neutcp2(neutct, 1D0, 2D0, 4D0)
      Zn24i = -neutcp2(neutct, 2D0, 2D0, 4D0)
      Zn31r = neutcp2(neutct, 1D0, 3D0, 1D0)
      Zn31i = -neutcp2(neutct, 2D0, 3D0, 1D0)
      Zn32r = neutcp2(neutct, 1D0, 3D0, 2D0)
      Zn32i = -neutcp2(neutct, 2D0, 3D0, 2D0)
      Zn33r = neutcp2(neutct, 1D0, 3D0, 3D0)
      Zn33i = -neutcp2(neutct, 2D0, 3D0, 3D0)
      Zn34r = neutcp2(neutct, 1D0, 3D0, 4D0)
      Zn34i = -neutcp2(neutct, 2D0, 3D0, 4D0)
      Zn41r = neutcp2(neutct, 1D0, 4D0, 1D0)
      Zn41i = -neutcp2(neutct, 2D0, 4D0, 1D0)
      Zn42r = neutcp2(neutct, 1D0, 4D0, 2D0)
      Zn42i = -neutcp2(neutct, 2D0, 4D0, 2D0)
      Zn43r = neutcp2(neutct, 1D0, 4D0, 3D0)
      Zn43i = -neutcp2(neutct, 2D0, 4D0, 3D0)
      Zn44r = neutcp2(neutct, 1D0, 4D0, 4D0)
      Zn44i = -neutcp2(neutct, 2D0, 4D0, 4D0)
      MNE1 = neutcp2(neutct, 3D0, 1D0, 1D0)
      MNE2 = neutcp2(neutct, 3D0, 2D0, 2D0)
      MNE3 = neutcp2(neutct, 3D0, 3D0, 3D0)
      MNE4 = neutcp2(neutct, 3D0, 4D0, 4D0)
	  else
	  MNE1=slhaVal1('MASS',Q,1000022)
	  MNE2=slhaVal1('MASS',Q,1000023)
	  MNE3=slhaVal1('MASS',Q,1000025)
	  MNE4=slhaVal1('MASS',Q,1000035)
	  Zn11r=slhaVal2('NMIX',Q,1,1)
	  Zn12r=slhaVal2('NMIX',Q,2,1)
	  Zn13r=slhaVal2('NMIX',Q,3,1)
	  Zn14r=slhaVal2('NMIX',Q,4,1)
	  Zn21r=slhaVal2('NMIX',Q,1,2)
	  Zn22r=slhaVal2('NMIX',Q,2,2)
	  Zn23r=slhaVal2('NMIX',Q,3,2)
	  Zn24r=slhaVal2('NMIX',Q,4,2)
	  Zn31r=slhaVal2('NMIX',Q,1,3)
	  Zn32r=slhaVal2('NMIX',Q,2,3)
	  Zn33r=slhaVal2('NMIX',Q,3,3)
	  Zn34r=slhaVal2('NMIX',Q,4,3)
	  Zn41r=slhaVal2('NMIX',Q,1,4)
	  Zn42r=slhaVal2('NMIX',Q,2,4)
	  Zn43r=slhaVal2('NMIX',Q,3,4)
	  Zn44r=slhaVal2('NMIX',Q,4,4)
	  Zn11i=-slhaVal2('IMNMIX',Q,1,1)
	  Zn12i=-slhaVal2('IMNMIX',Q,2,1)
	  Zn13i=-slhaVal2('IMNMIX',Q,3,1)
	  Zn14i=-slhaVal2('IMNMIX',Q,4,1)
	  Zn21i=-slhaVal2('IMNMIX',Q,1,2)
	  Zn22i=-slhaVal2('IMNMIX',Q,2,2)
	  Zn23i=-slhaVal2('IMNMIX',Q,3,2)
	  Zn24i=-slhaVal2('IMNMIX',Q,4,2)
	  Zn31i=-slhaVal2('IMNMIX',Q,1,3)
	  Zn32i=-slhaVal2('IMNMIX',Q,2,3)
	  Zn33i=-slhaVal2('IMNMIX',Q,3,3)
	  Zn34i=-slhaVal2('IMNMIX',Q,4,3)
	  Zn41i=-slhaVal2('IMNMIX',Q,1,4)
	  Zn42i=-slhaVal2('IMNMIX',Q,2,4)
	  Zn43i=-slhaVal2('IMNMIX',Q,3,4)
	  Zn44i=-slhaVal2('IMNMIX',Q,4,4)
	  endif
	  
	  if(mode.eq.1) then
      chatkn = chacp1(MZ, MW, SW, aMG2, amu, tb, fiMG2, fimu)
      Zm11r = chacp2(chatkn, 1D0, 1D0, 1D0)
      Zp11r = chacp2(chatkn, 3D0, 1D0, 1D0)
      Zm11i = chacp2(chatkn, 2D0, 1D0, 1D0)
      Zp11i = -chacp2(chatkn, 4D0, 1D0, 1D0)
      Zm12r = chacp2(chatkn, 1D0, 1D0, 2D0)
      Zp12r = chacp2(chatkn, 3D0, 1D0, 2D0)
      Zm12i = chacp2(chatkn, 2D0, 1D0, 2D0)
      Zp12i = -chacp2(chatkn, 4D0, 1D0, 2D0)
      Zm21r = chacp2(chatkn, 1D0, 2D0, 1D0)
      Zp21r = chacp2(chatkn, 3D0, 2D0, 1D0)
      Zm21i = chacp2(chatkn, 2D0, 2D0, 1D0)
      Zp21i = -chacp2(chatkn, 4D0, 2D0, 1D0)
      Zm22r = chacp2(chatkn, 1D0, 2D0, 2D0)
      Zp22r = chacp2(chatkn, 3D0, 2D0, 2D0)
      Zm22i = chacp2(chatkn, 2D0, 2D0, 2D0)
      Zp22i = -chacp2(chatkn, 4D0, 2D0, 2D0)
      MC1 = chacp2(chatkn, 5D0, 1D0, 1D0)
      MC2 = chacp2(chatkn, 5D0, 2D0, 2D0)
	  else
	  MC1=slhaVal1('MASS',Q,1000024)
	  MC2=slhaVal1('MASS',Q,1000037)
	  Zp11r=slhaVal2('VMIX',Q,1,1)
	  Zp12r=slhaVal2('VMIX',Q,2,1)
	  Zp21r=slhaVal2('VMIX',Q,1,2)
	  Zp22r=slhaVal2('VMIX',Q,2,2)
	  Zp11i=slhaVal2('IMVMIX',Q,1,1)
	  Zp12i=slhaVal2('IMVMIX',Q,2,1)
	  Zp21i=slhaVal2('IMVMIX',Q,1,2)
	  Zp22i=slhaVal2('IMVMIX',Q,2,2)
	  Zm11r=slhaVal2('UMIX',Q,1,1)
	  Zm12r=slhaVal2('UMIX',Q,2,1)
	  Zm21r=slhaVal2('UMIX',Q,1,2)
	  Zm22r=slhaVal2('UMIX',Q,2,2)
	  Zm11i=slhaVal2('IMUMIX',Q,1,1)
	  Zm12i=slhaVal2('IMUMIX',Q,2,1)
	  Zm21i=slhaVal2('IMUMIX',Q,1,2)
	  Zm22i=slhaVal2('IMUMIX',Q,2,2)
	  endif

      MSG = MG3
	  
	  if(mode.eq.1) then
	  inidi=initdiagonal()
      hM11 = MZ**2-2D0*MZ**2*sb**2+MZ**2*SW**2*sb**2+MHc
     &      **2*sb**2
      hM12 = -MZ**2*SW**2*cb*sb-MHc**2*cb*sb
      hM22 = -MZ**2+MZ**2*SW**2+2D0*MZ**2*sb**2-MZ**2*SW
     &      **2*sb**2+MHc**2*cb**2
      hM33 = -MW**2+MHc**2
      htk = rDiagonal3(hM11, hM12, 0D0, hM22, 0D0, hM33
     &      )
      Mh1 = sqrt(MassArray(htk, 1))
      Mh2 = sqrt(MassArray(htk, 2))
      Mh3 = sqrt(MassArray(htk, 3))
      Zh11 = MixMatrix(htk, 1, 1)
      Zh12 = MixMatrix(htk, 2, 1)
      Zh13 = MixMatrix(htk, 3, 1)
      Zh21 = MixMatrix(htk, 1, 2)
      Zh22 = MixMatrix(htk, 2, 2)
      Zh23 = MixMatrix(htk, 3, 2)
      Zh31 = MixMatrix(htk, 1, 3)
      Zh32 = MixMatrix(htk, 2, 3)
      Zh33 = MixMatrix(htk, 3, 3)
	  else
	  Mh1=slhaVal1('MASS',Q,25)
	  Mh2=slhaVal1('MASS',Q,35)
	  Mh3=slhaVal1('MASS',Q,36)
      Zh11=slhaVal2('CVHMIX',Q,1,1)
      Zh12=slhaVal2('CVHMIX',Q,2,1)
      Zh13=slhaVal2('CVHMIX',Q,3,1)
      Zh21=slhaVal2('CVHMIX',Q,1,2)
      Zh22=slhaVal2('CVHMIX',Q,2,2)
      Zh23=slhaVal2('CVHMIX',Q,3,2)
      Zh31=slhaVal2('CVHMIX',Q,1,3)/sb
      Zh32=slhaVal2('CVHMIX',Q,2,3)/sb
      Zh33=slhaVal2('CVHMIX',Q,3,3)/sb
	  endif
	  
      MSbLL = Mq3**2+Mb**2-MW**2/CW**2*(1D0/2D0-1D0/3D0*SW
     &      **2)*c2b
      MSbLRr = Mb*(Abr-mur*tb)
      MSbLRi = Mb*(Abi+mui*tb)
      MSbRR = Md3**2+Mb**2-1D0/3D0*MW**2/CW**2*SW**2*c2b
      MSb1 = sqrt((MSbLL+MSbRR-sqrt((MSbLL-MSbRR)**2+4D0
     &      *(MSbLRr**2+MSbLRi**2)))/2D0)
      MSb2 = sqrt((MSbLL+MSbRR+sqrt((MSbLL-MSbRR)**2+4D0
     &      *(MSbLRr**2+MSbLRi**2)))/2D0)
      MSbth = atan2(-2D0*sqrt(MSbLRr**2+MSbLRi**2), -MSbLL
     &      +MSbRR)*MSbLRr/sqrt(MSbLRr**2)/2D0
      MSbfi = atan(-MSbLRi/MSbLRr)
      Zd33r = cos(MSbth)
      Zd36r = -sin(MSbth)*cos(MSbfi)
      Zd63r = -Zd36r
      Zd66r = Zd33r
      Zd33i = sin(0D0)
      Zd36i = -sin(MSbth)*sin(MSbfi)
      Zd63i = Zd36i
      Zd66i = Zd33i
      MStLL = Mq3**2+MtopS**2+MW**2/CW**2*(1D0/2D0-2D0/3D0
     &      *SW**2)*c2b
      MStLRr = MtopS*(Atr-mur/tb)
      MStLRi = MtopS*(-Ati-mui/tb)
      MStRR = Mu3**2+MtopS**2+2D0/3D0*MW**2/CW**2*SW**2*c2b
      MSt1 = sqrt((MStLL+MStRR-sqrt((MStLL-MStRR)**2+4D0
     &      *(MStLRr**2+MStLRi**2)))/2D0)
      MSt2 = sqrt((MStLL+MStRR+sqrt((MStLL-MStRR)**2+4D0
     &      *(MStLRr**2+MStLRi**2)))/2D0)
      MStth = atan2(-2D0*sqrt(MStLRr**2+MStLRi**2), -MStLL
     &      +MStRR)*MStLRr/sqrt(MStLRr**2)/2D0
      MStfi = atan(MStLRi/MStLRr)
      Zu33r = cos(MStth)
      Zu36r = -sin(MStth)*cos(MStfi)
      Zu63r = -Zu36r
      Zu66r = Zu33r
      Zu33i = sin(0D0)
      Zu36i = -sin(MStth)*sin(MStfi)
      Zu63i = Zu36i
      Zu66i = Zu33i
      MSlLL = Ml3**2+Ml**2-MW**2/CW**2*(1D0/2D0-SW**2)*c2b
      MSlLRr = Ml*(Alr-mur*tb)
      MSlLRi = Ml*(Ali+mui*tb)
      MSlRR = Mr3**2+Ml**2-MW**2/CW**2*SW**2*c2b
      MSl1 = sqrt((MSlLL+MSlRR-sqrt((MSlLL-MSlRR)**2+4D0
     &      *(MSlLRr**2+MSlLRi**2)))/2D0)
      MSl2 = sqrt((MSlLL+MSlRR+sqrt((MSlLL-MSlRR)**2+4D0
     &      *(MSlLRr**2+MSlLRi**2)))/2D0)
      MSlth = atan2(-2D0*sqrt(MSlLRr**2+MSlLRi**2), -MSlLL
     &      +MSlRR)*MSlLRr/sqrt(MSlLRr**2)/2D0
      MSlfi = atan(-MSlLRi/MSlLRr)
      Zl33r = cos(MSlth)
      Zl36r = -sin(MSlth)*cos(MSlfi)
      Zl63r = -Zl36r
      Zl66r = Zl33r
      Zl36i = -sin(MSlth)*sin(MSlfi)
      Zl63i = Zl36i
      Zl66i = Zl33i
	  
	  if(mode.eq.0) then
	  MSl1=slhaVal1('MASS',Q,1000015)
	  MSl2=slhaVal1('MASS',Q,2000015)
	  Zl33r=slhaVal2('STAUMIX',Q,1,1)
	  Zl36r=slhaVal2('STAUMIX',Q,2,1)
	  Zl63r=slhaVal2('STAUMIX',Q,1,2)
	  Zl66r=slhaVal2('STAUMIX',Q,2,2)
	  Zl33i=-slhaVal2('IMSTAUMIX',Q,1,1)
	  Zl36i=-slhaVal2('IMSTAUMIX',Q,2,1)
	  Zl63i=-slhaVal2('IMSTAUMIX',Q,1,2)
	  Zl66i=-slhaVal2('IMSTAUMIX',Q,2,2)
	  
	  MSnl=slhaVal1('MASS',Q,1000016)
	  
	  MSt1=slhaVal1('MASS',Q,1000006)
	  MSt2=slhaVal1('MASS',Q,2000006)
	  Zu33r=slhaVal2('STOPMIX',Q,1,1)
	  Zu36r=slhaVal2('STOPMIX',Q,2,1)
	  Zu63r=slhaVal2('STOPMIX',Q,1,2)
	  Zu66r=slhaVal2('STOPMIX',Q,2,2)
	  Zu33i=-slhaVal2('IMSTOPMIX',Q,1,1)
	  Zu36i=-slhaVal2('IMSTOPMIX',Q,2,1)
	  Zu63i=-slhaVal2('IMSTOPMIX',Q,1,2)
	  Zu66i=-slhaVal2('IMSTOPMIX',Q,2,2)
	  	  
	  MSb1=slhaVal1('MASS',Q,1000005)
	  MSb2=slhaVal1('MASS',Q,2000005)
	  Zd33r=slhaVal2('SBOTMIX',Q,1,1)
	  Zd36r=slhaVal2('SBOTMIX',Q,2,1)
	  Zd63r=slhaVal2('SBOTMIX',Q,1,2)
	  Zd66r=slhaVal2('SBOTMIX',Q,2,2)
	  Zd33i=-slhaVal2('IMSBOTMIX',Q,1,1)
	  Zd36i=-slhaVal2('IMSBOTMIX',Q,2,1)
	  Zd63i=-slhaVal2('IMSBOTMIX',Q,1,2)
	  Zd66i=-slhaVal2('IMSBOTMIX',Q,2,2)
	  
	  
	  endif
	  
	  
      Mh1sq = Mh1**2
      Mh2sq = Mh2**2
      Mh3sq = Mh3**2
      Mesq = Me**2
      Mmsq = Mm**2
      Mlsq = Ml**2
      Musq = Mu**2
      Mdsq = Md**2
      Mcsq = Mc**2
      Mssq = Ms**2
      Mbsq = Mb**2
      Mtsq = Mt**2
      MSe1sq = MSe1**2
      MSm1sq = MSm1**2
      MSl1sq = MSl1**2
      MSu1sq = MSu1**2
      MSd1sq = MSd1**2
      MSc1sq = MSc1**2
      MSs1sq = MSs1**2
      MSb1sq = MSb1**2
      MSt1sq = MSt1**2
      MSe2sq = MSe2**2
      MSm2sq = MSm2**2
      MSl2sq = MSl2**2
      MSu2sq = MSu2**2
      MSd2sq = MSd2**2
      MSc2sq = MSc2**2
      MSs2sq = MSs2**2
      MSb2sq = MSb2**2
      MSt2sq = MSt2**2
      MSnesq = MSne**2
      MSnmsq = MSnm**2
      MSnlsq = MSnl**2
	  
	  call ModelDigest
	  

      call aaini01
      call aaini02
      call mtrini
      end

      subroutine aaini01
      implicit none
#include "model.h"

      AAABR(1) = MW*(1D0/CW-CW-SW**2/CW)
      AAABR(2) = CW**2+SW**2
      AAABR(3) = 1D0-CW**2-SW**2
      AAABR(4) = MW**2*(CW**2+2D0*SW**2+SW**4/CW**2)
      AAABR(5) = EE*MW
      AAABR(6) = EE*MW/SW
      AAABR(7) = EE*MW/SW*(Zh12*cb+Zh22*sb)
      AAABR(8) = EE*MW/SW*(Zh13*cb+Zh23*sb)
      AAABR(9) = EE*MW/SW*(Zh11*cb+Zh21*sb)
      AAABR(10) = EE*MW*(SW/CW-CW/SW)
      AAABR(11) = EE/CW*MW/SW
      AAABR(12) = EE/CW**2*MW/SW*(Zh11*cb+Zh21*sb)
      AAABR(13) = EE/CW**2*MW/SW*(Zh13*cb+Zh23*sb)
      AAABR(14) = EE/CW**2*MW/SW*(Zh12*cb+Zh22*sb)
      AAABR(15) = CW*EE/SW
      AAABR(16) = EE*MW*(Zh12/SW*cb+3D0*Zh22/SW*sb-2D0*Zh22
     &      /SW*sb**3-SW/CW**2*Zh12*cb+SW/CW**2*Zh22*sb-2D0
     &      *SW/CW**2*Zh22*sb**3+2D0*Zh12/SW*cb*sb**2+2D0
     &      *SW/CW**2*Zh12*cb*sb**2)
      AAABR(17) = EE*MW*(Zh13/SW*cb+3D0*Zh23/SW*sb-2D0*Zh23
     &      /SW*sb**3-SW/CW**2*Zh13*cb+SW/CW**2*Zh23*sb-2D0
     &      *SW/CW**2*Zh23*sb**3+2D0*Zh13/SW*cb*sb**2+2D0
     &      *SW/CW**2*Zh13*cb*sb**2)
      AAABR(18) = EE*MW*(Zh11/SW*cb+3D0*Zh21/SW*sb-2D0*Zh21
     &      /SW*sb**3-SW/CW**2*Zh11*cb+SW/CW**2*Zh21*sb-2D0
     &      *SW/CW**2*Zh21*sb**3+2D0*Zh11/SW*cb*sb**2+2D0
     &      *SW/CW**2*Zh11*cb*sb**2)
      AAABR(19) = 2D0*Zl41*ls1*sb+EE/MW*Me/SW*Sqrt2*Zl41
     &      *mur-2D0*EE*MW/SW*Sqrt2*Zl11/cb*sb+2D0*EE*MW
     &      /SW*Sqrt2*Zl11/cb*sb**3+EE/MW*Me**2/SW*Sqrt2
     &      *Zl11/cb*sb
      AAABR(20) = EE/MW*Me/SW*Sqrt2*Zl41*mui
      AAABR(21) = 2D0*Zl52*ls2*sb+EE/MW*Mm/SW*Sqrt2*Zl52
     &      *mur-2D0*EE*MW/SW*Sqrt2*Zl22/cb*sb+2D0*EE*MW
     &      /SW*Sqrt2*Zl22/cb*sb**3+EE/MW*Mm**2/SW*Sqrt2
     &      *Zl22/cb*sb
      AAABR(22) = EE/MW*Mm/SW*Sqrt2*Zl52*mui
      AAABR(23) = EE/SW*Sqrt2*(Ml/MW*Zl63i*mui-Ml/MW*Zl63r
     &      *mur+2D0*MW*Zl33r/cb*sb-2D0*MW*Zl33r/cb*sb**3
     &      -Ml**2/MW*Zl33r/cb*sb-Ali/MW*Ml*Zl63i/cb*sb-Alr
     &      /MW*Ml*Zl63r/cb*sb)
      AAABR(24) = EE/SW*Sqrt2*(Ml/MW*Zl63i*mur+Ml/MW*Zl63r
     &      *mui-2D0*MW*Zl33i/cb*sb+2D0*MW*Zl33i/cb*sb**3
     &      +Ml**2/MW*Zl33i/cb*sb-Ali/MW*Ml*Zl63r/cb*sb+Alr
     &      /MW*Ml*Zl63i/cb*sb)
      AAABR(25) = 2D0*Zl55*ls2*sb+EE/MW*Mm/SW*Sqrt2*Zl55
     &      *mur-2D0*EE*MW/SW*Sqrt2*Zl25/cb*sb+2D0*EE*MW
     &      /SW*Sqrt2*Zl25/cb*sb**3+EE/MW*Mm**2/SW*Sqrt2
     &      *Zl25/cb*sb
      AAABR(26) = EE/MW*Mm/SW*Sqrt2*Zl55*mui
      AAABR(27) = EE/SW*Sqrt2*(Ml/MW*Zl66i*mui-Ml/MW*Zl66r
     &      *mur+2D0*MW*Zl36r/cb*sb-2D0*MW*Zl36r/cb*sb**3
     &      -Ml**2/MW*Zl36r/cb*sb-Ali/MW*Ml*Zl66i/cb*sb-Alr
     &      /MW*Ml*Zl66r/cb*sb)
      AAABR(28) = EE/SW*Sqrt2*(Ml/MW*Zl66i*mur+Ml/MW*Zl66r
     &      *mui-2D0*MW*Zl36i/cb*sb+2D0*MW*Zl36i/cb*sb**3
     &      +Ml**2/MW*Zl36i/cb*sb-Ali/MW*Ml*Zl66r/cb*sb+Alr
     &      /MW*Ml*Zl66i/cb*sb)
      AAABR(29) = 2D0*Zl44*ls1*sb+EE/MW*Me/SW*Sqrt2*Zl44
     &      *mur-2D0*EE*MW/SW*Sqrt2*Zl14/cb*sb+2D0*EE*MW
     &      /SW*Sqrt2*Zl14/cb*sb**3+EE/MW*Me**2/SW*Sqrt2
     &      *Zl14/cb*sb
      AAABR(30) = EE/MW*Me/SW*Sqrt2*Zl44*mui
      AAABR(31) = Zd11*Zu41*cb*us1-Zd41*Zu11*ds1*sb+EE*MW
     &      /SW*Sqrt2*Zd11*Zu11*cb*sb
      AAABR(32) = EE/SW*Sqrt2*(Mb/MW*Zd63i*Zu33i*mur-Mb/MW
     &      *Zd63i*Zu33r*mui+Mb/MW*Zd63r*Zu33i*mui+Mb/MW
     &      *Zd63r*Zu33r*mur+MtopS/MW*Zd33i*Zu63i*mur+MtopS
     &      /MW*Zd33i*Zu63r*mui-MtopS/MW*Zd33r*Zu63i*mui
     &      +MtopS/MW*Zd33r*Zu63r*mur-2D0*MW*Zd33i*Zu33i
     &      /cb*sb+2D0*MW*Zd33i*Zu33i/cb*sb**3-2D0*MW*Zd33r
     &      *Zu33r/cb*sb+2D0*MW*Zd33r*Zu33r/cb*sb**3+Mb**2
     &      /MW*Zd33i*Zu33i/cb*sb+Mb**2/MW*Zd33r*Zu33r/cb
     &      *sb+MtopS**2/MW*Zd33i*Zu33i/cb/sb-MtopS**2/MW
     &      *Zd33i*Zu33i/cb*sb+MtopS**2/MW*Zd33r*Zu33r/cb
     &      /sb-MtopS**2/MW*Zd33r*Zu33r/cb*sb+Abi/MW*Mb*Zd63i
     &      *Zu33r/cb*sb-Abi/MW*Mb*Zd63r*Zu33i/cb*sb+Abr
     &      /MW*Mb*Zd63i*Zu33i/cb*sb+Abr/MW*Mb*Zd63r*Zu33r
     &      /cb*sb-Ati/MW*MtopS*Zd33i*Zu63r/cb/sb+Ati/MW
     &      *MtopS*Zd33i*Zu63r/cb*sb+Ati/MW*MtopS*Zd33r*Zu63i
     &      /cb/sb-Ati/MW*MtopS*Zd33r*Zu63i/cb*sb+Atr/MW
     &      *MtopS*Zd33i*Zu63i/cb/sb-Atr/MW*MtopS*Zd33i*Zu63i
     &      /cb*sb+Atr/MW*MtopS*Zd33r*Zu63r/cb/sb-Atr/MW
     &      *MtopS*Zd33r*Zu63r/cb*sb+Mb/MW*MtopS*Zd63i*Zu63i
     &      /cb/sb+Mb/MW*MtopS*Zd63r*Zu63r/cb/sb)
      AAABR(33) = EE/SW*Sqrt2*(Mb/MW*Zd63i*Zu33i*mui+Mb/MW
     &      *Zd63i*Zu33r*mur-Mb/MW*Zd63r*Zu33i*mur+Mb/MW
     &      *Zd63r*Zu33r*mui-MtopS/MW*Zd33i*Zu63i*mui+MtopS
     &      /MW*Zd33i*Zu63r*mur-MtopS/MW*Zd33r*Zu63i*mur
     &      -MtopS/MW*Zd33r*Zu63r*mui-2D0*MW*Zd33i*Zu33r
     &      /cb*sb+2D0*MW*Zd33i*Zu33r/cb*sb**3+2D0*MW*Zd33r
     &      *Zu33i/cb*sb-2D0*MW*Zd33r*Zu33i/cb*sb**3+Mb**2
     &      /MW*Zd33i*Zu33r/cb*sb-Mb**2/MW*Zd33r*Zu33i/cb
     &      *sb+MtopS**2/MW*Zd33i*Zu33r/cb/sb-MtopS**2/MW
     &      *Zd33i*Zu33r/cb*sb-MtopS**2/MW*Zd33r*Zu33i/cb
     &      /sb+MtopS**2/MW*Zd33r*Zu33i/cb*sb-Abi/MW*Mb*Zd63i
     &      *Zu33i/cb*sb-Abi/MW*Mb*Zd63r*Zu33r/cb*sb+Abr
     &      /MW*Mb*Zd63i*Zu33r/cb*sb-Abr/MW*Mb*Zd63r*Zu33i
     &      /cb*sb+Ati/MW*MtopS*Zd33i*Zu63i/cb/sb-Ati/MW
     &      *MtopS*Zd33i*Zu63i/cb*sb+Ati/MW*MtopS*Zd33r*Zu63r
     &      /cb/sb-Ati/MW*MtopS*Zd33r*Zu63r/cb*sb+Atr/MW
     &      *MtopS*Zd33i*Zu63r/cb/sb-Atr/MW*MtopS*Zd33i*Zu63r
     &      /cb*sb-Atr/MW*MtopS*Zd33r*Zu63i/cb/sb+Atr/MW
     &      *MtopS*Zd33r*Zu63i/cb*sb+Mb/MW*MtopS*Zd63i*Zu63r
     &      /cb/sb-Mb/MW*MtopS*Zd63r*Zu63i/cb/sb)
      AAABR(34) = 2D0*Zd22*Zu52/cb*us2-2D0*Zd52*Zu22*ds2
     &      *sb-2D0*Zd22*Zu52/cb*sb**2*us2-EE/MW*Mc/SW*Sqrt2
     &      *Zd22*Zu52*mur-EE/MW*Ms/SW*Sqrt2*Zd52*Zu22*mur
     &      +2D0*EE*MW/SW*Sqrt2*Zd22*Zu22/cb*sb-2D0*EE*MW
     &      /SW*Sqrt2*Zd22*Zu22/cb*sb**3-EE/MW*Mc**2/SW*Sqrt2
     &      *Zd22*Zu22/cb/sb+EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu22
     &      /cb*sb-EE/MW*Ms**2/SW*Sqrt2*Zd22*Zu22/cb*sb-EE
     &      /MW*Mc*Ms/SW*Sqrt2*Zd52*Zu52/cb/sb
      AAABR(35) = EE/MW/SW*Sqrt2*mui*(Mc*Zd22*Zu52-Ms*Zd52
     &      *Zu22)
      AAABR(36) = EE/SW*Sqrt2*(Mb/MW*Zd66i*Zu33i*mur-Mb/MW
     &      *Zd66i*Zu33r*mui+Mb/MW*Zd66r*Zu33i*mui+Mb/MW
     &      *Zd66r*Zu33r*mur+MtopS/MW*Zd36i*Zu63i*mur+MtopS
     &      /MW*Zd36i*Zu63r*mui-MtopS/MW*Zd36r*Zu63i*mui
     &      +MtopS/MW*Zd36r*Zu63r*mur-2D0*MW*Zd36i*Zu33i
     &      /cb*sb+2D0*MW*Zd36i*Zu33i/cb*sb**3-2D0*MW*Zd36r
     &      *Zu33r/cb*sb+2D0*MW*Zd36r*Zu33r/cb*sb**3+Mb**2
     &      /MW*Zd36i*Zu33i/cb*sb+Mb**2/MW*Zd36r*Zu33r/cb
     &      *sb+MtopS**2/MW*Zd36i*Zu33i/cb/sb-MtopS**2/MW
     &      *Zd36i*Zu33i/cb*sb+MtopS**2/MW*Zd36r*Zu33r/cb
     &      /sb-MtopS**2/MW*Zd36r*Zu33r/cb*sb+Abi/MW*Mb*Zd66i
     &      *Zu33r/cb*sb-Abi/MW*Mb*Zd66r*Zu33i/cb*sb+Abr
     &      /MW*Mb*Zd66i*Zu33i/cb*sb+Abr/MW*Mb*Zd66r*Zu33r
     &      /cb*sb-Ati/MW*MtopS*Zd36i*Zu63r/cb/sb+Ati/MW
     &      *MtopS*Zd36i*Zu63r/cb*sb+Ati/MW*MtopS*Zd36r*Zu63i
     &      /cb/sb-Ati/MW*MtopS*Zd36r*Zu63i/cb*sb+Atr/MW
     &      *MtopS*Zd36i*Zu63i/cb/sb-Atr/MW*MtopS*Zd36i*Zu63i
     &      /cb*sb+Atr/MW*MtopS*Zd36r*Zu63r/cb/sb-Atr/MW
     &      *MtopS*Zd36r*Zu63r/cb*sb+Mb/MW*MtopS*Zd66i*Zu63i
     &      /cb/sb+Mb/MW*MtopS*Zd66r*Zu63r/cb/sb)
      AAABR(37) = EE/SW*Sqrt2*(Mb/MW*Zd66i*Zu33i*mui+Mb/MW
     &      *Zd66i*Zu33r*mur-Mb/MW*Zd66r*Zu33i*mur+Mb/MW
     &      *Zd66r*Zu33r*mui-MtopS/MW*Zd36i*Zu63i*mui+MtopS
     &      /MW*Zd36i*Zu63r*mur-MtopS/MW*Zd36r*Zu63i*mur
     &      -MtopS/MW*Zd36r*Zu63r*mui-2D0*MW*Zd36i*Zu33r
     &      /cb*sb+2D0*MW*Zd36i*Zu33r/cb*sb**3+2D0*MW*Zd36r
     &      *Zu33i/cb*sb-2D0*MW*Zd36r*Zu33i/cb*sb**3+Mb**2
     &      /MW*Zd36i*Zu33r/cb*sb-Mb**2/MW*Zd36r*Zu33i/cb
     &      *sb+MtopS**2/MW*Zd36i*Zu33r/cb/sb-MtopS**2/MW
     &      *Zd36i*Zu33r/cb*sb-MtopS**2/MW*Zd36r*Zu33i/cb
     &      /sb+MtopS**2/MW*Zd36r*Zu33i/cb*sb-Abi/MW*Mb*Zd66i
     &      *Zu33i/cb*sb-Abi/MW*Mb*Zd66r*Zu33r/cb*sb+Abr
     &      /MW*Mb*Zd66i*Zu33r/cb*sb-Abr/MW*Mb*Zd66r*Zu33i
     &      /cb*sb+Ati/MW*MtopS*Zd36i*Zu63i/cb/sb-Ati/MW
     &      *MtopS*Zd36i*Zu63i/cb*sb+Ati/MW*MtopS*Zd36r*Zu63r
     &      /cb/sb-Ati/MW*MtopS*Zd36r*Zu63r/cb*sb+Atr/MW
     &      *MtopS*Zd36i*Zu63r/cb/sb-Atr/MW*MtopS*Zd36i*Zu63r
     &      /cb*sb-Atr/MW*MtopS*Zd36r*Zu63i/cb/sb+Atr/MW
     &      *MtopS*Zd36r*Zu63i/cb*sb+Mb/MW*MtopS*Zd66i*Zu63r
     &      /cb/sb-Mb/MW*MtopS*Zd66r*Zu63i/cb/sb)
      AAABR(38) = 2D0*Zd25*Zu52/cb*us2-2D0*Zd55*Zu22*ds2
     &      *sb-2D0*Zd25*Zu52/cb*sb**2*us2-EE/MW*Mc/SW*Sqrt2
     &      *Zd25*Zu52*mur-EE/MW*Ms/SW*Sqrt2*Zd55*Zu22*mur
     &      +2D0*EE*MW/SW*Sqrt2*Zd25*Zu22/cb*sb-2D0*EE*MW
     &      /SW*Sqrt2*Zd25*Zu22/cb*sb**3-EE/MW*Mc**2/SW*Sqrt2
     &      *Zd25*Zu22/cb/sb+EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu22
     &      /cb*sb-EE/MW*Ms**2/SW*Sqrt2*Zd25*Zu22/cb*sb-EE
     &      /MW*Mc*Ms/SW*Sqrt2*Zd55*Zu52/cb/sb
      AAABR(39) = EE/MW/SW*Sqrt2*mui*(Mc*Zd25*Zu52-Ms*Zd55
     &      *Zu22)
      AAABR(40) = Zd14*Zu41*cb*us1-Zd44*Zu11*ds1*sb+EE*MW
     &      /SW*Sqrt2*Zd14*Zu11*cb*sb
      AAABR(41) = EE/SW*Sqrt2*(Mb/MW*Zd63i*Zu36i*mur-Mb/MW
     &      *Zd63i*Zu36r*mui+Mb/MW*Zd63r*Zu36i*mui+Mb/MW
     &      *Zd63r*Zu36r*mur+MtopS/MW*Zd33i*Zu66i*mur+MtopS
     &      /MW*Zd33i*Zu66r*mui-MtopS/MW*Zd33r*Zu66i*mui
     &      +MtopS/MW*Zd33r*Zu66r*mur-2D0*MW*Zd33i*Zu36i
     &      /cb*sb+2D0*MW*Zd33i*Zu36i/cb*sb**3-2D0*MW*Zd33r
     &      *Zu36r/cb*sb+2D0*MW*Zd33r*Zu36r/cb*sb**3+Mb**2
     &      /MW*Zd33i*Zu36i/cb*sb+Mb**2/MW*Zd33r*Zu36r/cb
     &      *sb+MtopS**2/MW*Zd33i*Zu36i/cb/sb-MtopS**2/MW
     &      *Zd33i*Zu36i/cb*sb+MtopS**2/MW*Zd33r*Zu36r/cb
     &      /sb-MtopS**2/MW*Zd33r*Zu36r/cb*sb+Abi/MW*Mb*Zd63i
     &      *Zu36r/cb*sb-Abi/MW*Mb*Zd63r*Zu36i/cb*sb+Abr
     &      /MW*Mb*Zd63i*Zu36i/cb*sb+Abr/MW*Mb*Zd63r*Zu36r
     &      /cb*sb-Ati/MW*MtopS*Zd33i*Zu66r/cb/sb+Ati/MW
     &      *MtopS*Zd33i*Zu66r/cb*sb+Ati/MW*MtopS*Zd33r*Zu66i
     &      /cb/sb-Ati/MW*MtopS*Zd33r*Zu66i/cb*sb+Atr/MW
     &      *MtopS*Zd33i*Zu66i/cb/sb-Atr/MW*MtopS*Zd33i*Zu66i
     &      /cb*sb+Atr/MW*MtopS*Zd33r*Zu66r/cb/sb-Atr/MW
     &      *MtopS*Zd33r*Zu66r/cb*sb+Mb/MW*MtopS*Zd63i*Zu66i
     &      /cb/sb+Mb/MW*MtopS*Zd63r*Zu66r/cb/sb)
      AAABR(42) = EE/SW*Sqrt2*(Mb/MW*Zd63i*Zu36i*mui+Mb/MW
     &      *Zd63i*Zu36r*mur-Mb/MW*Zd63r*Zu36i*mur+Mb/MW
     &      *Zd63r*Zu36r*mui-MtopS/MW*Zd33i*Zu66i*mui+MtopS
     &      /MW*Zd33i*Zu66r*mur-MtopS/MW*Zd33r*Zu66i*mur
     &      -MtopS/MW*Zd33r*Zu66r*mui-2D0*MW*Zd33i*Zu36r
     &      /cb*sb+2D0*MW*Zd33i*Zu36r/cb*sb**3+2D0*MW*Zd33r
     &      *Zu36i/cb*sb-2D0*MW*Zd33r*Zu36i/cb*sb**3+Mb**2
     &      /MW*Zd33i*Zu36r/cb*sb-Mb**2/MW*Zd33r*Zu36i/cb
     &      *sb+MtopS**2/MW*Zd33i*Zu36r/cb/sb-MtopS**2/MW
     &      *Zd33i*Zu36r/cb*sb-MtopS**2/MW*Zd33r*Zu36i/cb
     &      /sb+MtopS**2/MW*Zd33r*Zu36i/cb*sb-Abi/MW*Mb*Zd63i
     &      *Zu36i/cb*sb-Abi/MW*Mb*Zd63r*Zu36r/cb*sb+Abr
     &      /MW*Mb*Zd63i*Zu36r/cb*sb-Abr/MW*Mb*Zd63r*Zu36i
     &      /cb*sb+Ati/MW*MtopS*Zd33i*Zu66i/cb/sb-Ati/MW
     &      *MtopS*Zd33i*Zu66i/cb*sb+Ati/MW*MtopS*Zd33r*Zu66r
     &      /cb/sb-Ati/MW*MtopS*Zd33r*Zu66r/cb*sb+Atr/MW
     &      *MtopS*Zd33i*Zu66r/cb/sb-Atr/MW*MtopS*Zd33i*Zu66r
     &      /cb*sb-Atr/MW*MtopS*Zd33r*Zu66i/cb/sb+Atr/MW
     &      *MtopS*Zd33r*Zu66i/cb*sb+Mb/MW*MtopS*Zd63i*Zu66r
     &      /cb/sb-Mb/MW*MtopS*Zd63r*Zu66i/cb/sb)
      AAABR(43) = 2D0*Zd22*Zu55/cb*us2-2D0*Zd52*Zu25*ds2
     &      *sb-2D0*Zd22*Zu55/cb*sb**2*us2-EE/MW*Mc/SW*Sqrt2
     &      *Zd22*Zu55*mur-EE/MW*Ms/SW*Sqrt2*Zd52*Zu25*mur
     &      +2D0*EE*MW/SW*Sqrt2*Zd22*Zu25/cb*sb-2D0*EE*MW
     &      /SW*Sqrt2*Zd22*Zu25/cb*sb**3-EE/MW*Mc**2/SW*Sqrt2
     &      *Zd22*Zu25/cb/sb+EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu25
     &      /cb*sb-EE/MW*Ms**2/SW*Sqrt2*Zd22*Zu25/cb*sb-EE
     &      /MW*Mc*Ms/SW*Sqrt2*Zd52*Zu55/cb/sb
      AAABR(44) = EE/MW/SW*Sqrt2*mui*(Mc*Zd22*Zu55-Ms*Zd52
     &      *Zu25)
      AAABR(45) = Zd11*Zu44*cb*us1-Zd41*Zu14*ds1*sb+EE*MW
     &      /SW*Sqrt2*Zd11*Zu14*cb*sb
      AAABR(46) = EE/SW*Sqrt2*(Mb/MW*Zd66i*Zu36i*mur-Mb/MW
     &      *Zd66i*Zu36r*mui+Mb/MW*Zd66r*Zu36i*mui+Mb/MW
     &      *Zd66r*Zu36r*mur+MtopS/MW*Zd36i*Zu66i*mur+MtopS
     &      /MW*Zd36i*Zu66r*mui-MtopS/MW*Zd36r*Zu66i*mui
     &      +MtopS/MW*Zd36r*Zu66r*mur-2D0*MW*Zd36i*Zu36i
     &      /cb*sb+2D0*MW*Zd36i*Zu36i/cb*sb**3-2D0*MW*Zd36r
     &      *Zu36r/cb*sb+2D0*MW*Zd36r*Zu36r/cb*sb**3+Mb**2
     &      /MW*Zd36i*Zu36i/cb*sb+Mb**2/MW*Zd36r*Zu36r/cb
     &      *sb+MtopS**2/MW*Zd36i*Zu36i/cb/sb-MtopS**2/MW
     &      *Zd36i*Zu36i/cb*sb+MtopS**2/MW*Zd36r*Zu36r/cb
     &      /sb-MtopS**2/MW*Zd36r*Zu36r/cb*sb+Abi/MW*Mb*Zd66i
     &      *Zu36r/cb*sb-Abi/MW*Mb*Zd66r*Zu36i/cb*sb+Abr
     &      /MW*Mb*Zd66i*Zu36i/cb*sb+Abr/MW*Mb*Zd66r*Zu36r
     &      /cb*sb-Ati/MW*MtopS*Zd36i*Zu66r/cb/sb+Ati/MW
     &      *MtopS*Zd36i*Zu66r/cb*sb+Ati/MW*MtopS*Zd36r*Zu66i
     &      /cb/sb-Ati/MW*MtopS*Zd36r*Zu66i/cb*sb+Atr/MW
     &      *MtopS*Zd36i*Zu66i/cb/sb-Atr/MW*MtopS*Zd36i*Zu66i
     &      /cb*sb+Atr/MW*MtopS*Zd36r*Zu66r/cb/sb-Atr/MW
     &      *MtopS*Zd36r*Zu66r/cb*sb+Mb/MW*MtopS*Zd66i*Zu66i
     &      /cb/sb+Mb/MW*MtopS*Zd66r*Zu66r/cb/sb)
      AAABR(47) = EE/SW*Sqrt2*(Mb/MW*Zd66i*Zu36i*mui+Mb/MW
     &      *Zd66i*Zu36r*mur-Mb/MW*Zd66r*Zu36i*mur+Mb/MW
     &      *Zd66r*Zu36r*mui-MtopS/MW*Zd36i*Zu66i*mui+MtopS
     &      /MW*Zd36i*Zu66r*mur-MtopS/MW*Zd36r*Zu66i*mur
     &      -MtopS/MW*Zd36r*Zu66r*mui-2D0*MW*Zd36i*Zu36r
     &      /cb*sb+2D0*MW*Zd36i*Zu36r/cb*sb**3+2D0*MW*Zd36r
     &      *Zu36i/cb*sb-2D0*MW*Zd36r*Zu36i/cb*sb**3+Mb**2
     &      /MW*Zd36i*Zu36r/cb*sb-Mb**2/MW*Zd36r*Zu36i/cb
     &      *sb+MtopS**2/MW*Zd36i*Zu36r/cb/sb-MtopS**2/MW
     &      *Zd36i*Zu36r/cb*sb-MtopS**2/MW*Zd36r*Zu36i/cb
     &      /sb+MtopS**2/MW*Zd36r*Zu36i/cb*sb-Abi/MW*Mb*Zd66i
     &      *Zu36i/cb*sb-Abi/MW*Mb*Zd66r*Zu36r/cb*sb+Abr
     &      /MW*Mb*Zd66i*Zu36r/cb*sb-Abr/MW*Mb*Zd66r*Zu36i
     &      /cb*sb+Ati/MW*MtopS*Zd36i*Zu66i/cb/sb-Ati/MW
     &      *MtopS*Zd36i*Zu66i/cb*sb+Ati/MW*MtopS*Zd36r*Zu66r
     &      /cb/sb-Ati/MW*MtopS*Zd36r*Zu66r/cb*sb+Atr/MW
     &      *MtopS*Zd36i*Zu66r/cb/sb-Atr/MW*MtopS*Zd36i*Zu66r
     &      /cb*sb-Atr/MW*MtopS*Zd36r*Zu66i/cb/sb+Atr/MW
     &      *MtopS*Zd36r*Zu66i/cb*sb+Mb/MW*MtopS*Zd66i*Zu66r
     &      /cb/sb-Mb/MW*MtopS*Zd66r*Zu66i/cb/sb)
      AAABR(48) = Zd14*Zu44*cb*us1-Zd44*Zu14*ds1*sb+EE*MW
     &      /SW*Sqrt2*Zd14*Zu14*cb*sb
      AAABR(49) = 2D0*Zd25*Zu55/cb*us2-2D0*Zd55*Zu25*ds2
     &      *sb-2D0*Zd25*Zu55/cb*sb**2*us2-EE/MW*Mc/SW*Sqrt2
     &      *Zd25*Zu55*mur-EE/MW*Ms/SW*Sqrt2*Zd55*Zu25*mur
     &      +2D0*EE*MW/SW*Sqrt2*Zd25*Zu25/cb*sb-2D0*EE*MW
     &      /SW*Sqrt2*Zd25*Zu25/cb*sb**3-EE/MW*Mc**2/SW*Sqrt2
     &      *Zd25*Zu25/cb/sb+EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu25
     &      /cb*sb-EE/MW*Ms**2/SW*Sqrt2*Zd25*Zu25/cb*sb-EE
     &      /MW*Mc*Ms/SW*Sqrt2*Zd55*Zu55/cb/sb
      AAABR(50) = EE/MW/SW*Sqrt2*mui*(Mc*Zd25*Zu55-Ms*Zd55
     &      *Zu25)
      AAABR(51) = EE*MW*(Zh11/SW*sb-2D0*Zh11/SW*sb**3+Zh21
     &      /SW*cb+2D0*SW/CW**2*Zh11*sb-2D0*SW/CW**2*Zh11
     &      *sb**3-2D0*Zh21/SW*cb*sb**2-2D0*SW/CW**2*Zh21
     &      *cb*sb**2)
      AAABR(52) = EE*MW/SW*Zh31
      AAABR(53) = EE*MW*(Zh12/SW*sb-2D0*Zh12/SW*sb**3+Zh22
     &      /SW*cb+2D0*SW/CW**2*Zh12*sb-2D0*SW/CW**2*Zh12
     &      *sb**3-2D0*Zh22/SW*cb*sb**2-2D0*SW/CW**2*Zh22
     &      *cb*sb**2)
      AAABR(54) = EE*MW/SW*Zh32
      AAABR(55) = EE*MW*(Zh13/SW*sb-2D0*Zh13/SW*sb**3+Zh23
     &      /SW*cb+2D0*SW/CW**2*Zh13*sb-2D0*SW/CW**2*Zh13
     &      *sb**3-2D0*Zh23/SW*cb*sb**2-2D0*SW/CW**2*Zh23
     &      *cb*sb**2)
      AAABR(56) = EE*MW/SW*Zh33
      AAABR(57) = 2D0*Zl52/cb*ls2-2D0*Zl52/cb*ls2*sb**2-EE
     &      *MW/SW*Sqrt2*Zl22+EE/MW*Mm**2/SW*Sqrt2*Zl22+2D0
     &      *EE*MW/SW*Sqrt2*Zl22*sb**2-EE/MW*Mm/SW*Sqrt2
     &      *Zl52/cb*mur*sb
      AAABR(58) = EE/MW*Mm/SW*Sqrt2*Zl52/cb*mui*sb
      AAABR(59) = 2D0*Zl41/cb*ls1-2D0*Zl41/cb*ls1*sb**2-EE
     &      *MW/SW*Sqrt2*Zl11+EE/MW*Me**2/SW*Sqrt2*Zl11+2D0
     &      *EE*MW/SW*Sqrt2*Zl11*sb**2-EE/MW*Me/SW*Sqrt2
     &      *Zl41/cb*mur*sb
      AAABR(60) = EE/MW*Me/SW*Sqrt2*Zl41/cb*mui*sb
      AAABR(61) = EE/SW*Sqrt2*(MW*Zl33r-Ml**2/MW*Zl33r-2D0
     &      *MW*Zl33r*sb**2-Ali/MW*Ml*Zl63i-Alr/MW*Ml*Zl63r
     &      -Ml/MW*Zl63i/cb*mui*sb+Ml/MW*Zl63r/cb*mur*sb)
      AAABR(62) = EE/SW*Sqrt2*(MW*Zl33i-Ml**2/MW*Zl33i-2D0
     &      *MW*Zl33i*sb**2+Ali/MW*Ml*Zl63r-Alr/MW*Ml*Zl63i
     &      +Ml/MW*Zl63i/cb*mur*sb+Ml/MW*Zl63r/cb*mui*sb)
      AAABR(63) = EE/MW*Me/SW*Zl11*Zl41*mui*tb
      AAABR(64) = EE/MW*Ml/SW*(Ali*Zl33i*Zl63i+Ali*Zl33r
     &      *Zl63r+Alr*Zl33i*Zl63r-Alr*Zl33r*Zl63i+Zl33i
     &      *Zl63i/cb*mui*sb-Zl33i*Zl63r/cb*mur*sb+Zl33r
     &      *Zl63i/cb*mur*sb+Zl33r*Zl63r/cb*mui*sb)
      AAABR(65) = EE/MW*Mm/SW*Zl22*Zl52*mui*tb
      AAABR(66) = 1D0/cb*(Sqrt2*Zl22*Zl55*ls2-Sqrt2*Zl25
     &      *Zl52*ls2-Sqrt2*Zl22*Zl55*ls2*sb**2+Sqrt2*Zl25
     &      *Zl52*ls2*sb**2-EE/MW*Mm/SW*Zl22*Zl55*mur*sb
     &      +EE/MW*Mm/SW*Zl25*Zl52*mur*sb)
      AAABR(67) = EE/MW*Mm/SW/cb*mui*sb*(Zl22*Zl55+Zl25*Zl52
     &      )
      AAABR(68) = 1D0/cb*(Sqrt2*Zl11*Zl44*ls1-Sqrt2*Zl14
     &      *Zl41*ls1-Sqrt2*Zl11*Zl44*ls1*sb**2+Sqrt2*Zl14
     &      *Zl41*ls1*sb**2-EE/MW*Me/SW*Zl11*Zl44*mur*sb
     &      +EE/MW*Me/SW*Zl14*Zl41*mur*sb)
      AAABR(69) = EE/MW*Me/SW/cb*mui*sb*(Zl11*Zl44+Zl14*Zl41
     &      )
      AAABR(70) = EE/MW*Ml/SW*(Ali*Zl33i*Zl66r-Ali*Zl33r
     &      *Zl66i-Ali*Zl36i*Zl63r+Ali*Zl36r*Zl63i-Alr*Zl33i
     &      *Zl66i-Alr*Zl33r*Zl66r+Alr*Zl36i*Zl63i+Alr*Zl36r
     &      *Zl63r+Zl33i*Zl66i/cb*mur*sb+Zl33i*Zl66r/cb*mui
     &      *sb-Zl33r*Zl66i/cb*mui*sb+Zl33r*Zl66r/cb*mur
     &      *sb-Zl36i*Zl63i/cb*mur*sb-Zl36i*Zl63r/cb*mui
     &      *sb+Zl36r*Zl63i/cb*mui*sb-Zl36r*Zl63r/cb*mur*sb)
      AAABR(71) = EE/MW*Ml/SW*(Ali*Zl33i*Zl66i+Ali*Zl33r
     &      *Zl66r+Ali*Zl36i*Zl63i+Ali*Zl36r*Zl63r+Alr*Zl33i
     &      *Zl66r-Alr*Zl33r*Zl66i+Alr*Zl36i*Zl63r-Alr*Zl36r
     &      *Zl63i+Zl33i*Zl66i/cb*mui*sb-Zl33i*Zl66r/cb*mur
     &      *sb+Zl33r*Zl66i/cb*mur*sb+Zl33r*Zl66r/cb*mui
     &      *sb+Zl36i*Zl63i/cb*mui*sb-Zl36i*Zl63r/cb*mur
     &      *sb+Zl36r*Zl63i/cb*mur*sb+Zl36r*Zl63r/cb*mui*sb)
      AAABR(72) = 2D0*Sqrt2*Zh11*Zl11*Zl41*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl11**2/cb+EE
     &      *MW/SW*Zh21*Zl11**2*sb+EE/CW**2*MW*SW*Zh11*Zl11
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl41**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl11**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl41**2*sb+EE*MW/SW*Zh11*Zl11**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl11**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl41**2/cb*sb**2+2D0*EE
     &      /MW*Me/SW*Zh31*Zl11*Zl41*mui-2D0*EE/MW*Me/SW
     &      *Zh21*Zl11*Zl41/cb*mur
      AAABR(73) = EE*(2D0*Ml**2/MW/SW*Zh11/cb-MW/SW*Zh11
     &      *Zl33i**2/cb-MW/SW*Zh11*Zl33r**2/cb+MW/SW*Zh21
     &      *Zl33i**2*sb+MW/SW*Zh21*Zl33r**2*sb+MW/CW**2
     &      *SW*Zh11*Zl33i**2/cb+MW/CW**2*SW*Zh11*Zl33r**2
     &      /cb-2D0*MW/CW**2*SW*Zh11*Zl63i**2/cb-2D0*MW/CW
     &      **2*SW*Zh11*Zl63r**2/cb-MW/CW**2*SW*Zh21*Zl33i
     &      **2*sb-MW/CW**2*SW*Zh21*Zl33r**2*sb+2D0*MW/CW
     &      **2*SW*Zh21*Zl63i**2*sb+2D0*MW/CW**2*SW*Zh21
     &      *Zl63r**2*sb+MW/SW*Zh11*Zl33i**2/cb*sb**2+MW
     &      /SW*Zh11*Zl33r**2/cb*sb**2-MW/CW**2*SW*Zh11*Zl33i
     &      **2/cb*sb**2-MW/CW**2*SW*Zh11*Zl33r**2/cb*sb
     &      **2+2D0*MW/CW**2*SW*Zh11*Zl63i**2/cb*sb**2+2D0
     &      *MW/CW**2*SW*Zh11*Zl63r**2/cb*sb**2+2D0*Ml/MW
     &      /SW*Zh31*Zl33i*Zl63i*mui-2D0*Ml/MW/SW*Zh31*Zl33i
     &      *Zl63r*mur+2D0*Ml/MW/SW*Zh31*Zl33r*Zl63i*mur
     &      +2D0*Ml/MW/SW*Zh31*Zl33r*Zl63r*mui-2D0*Ali/MW
     &      *Ml/SW*Zh11*Zl33i*Zl63r/cb+2D0*Ali/MW*Ml/SW*Zh11
     &      *Zl33r*Zl63i/cb+2D0*Alr/MW*Ml/SW*Zh11*Zl33i*Zl63i
     &      /cb+2D0*Alr/MW*Ml/SW*Zh11*Zl33r*Zl63r/cb-2D0
     &      *Ml/MW/SW*Zh21*Zl33i*Zl63i/cb*mur-2D0*Ml/MW/SW
     &      *Zh21*Zl33i*Zl63r/cb*mui+2D0*Ml/MW/SW*Zh21*Zl33r
     &      *Zl63i/cb*mui-2D0*Ml/MW/SW*Zh21*Zl33r*Zl63r/cb
     &      *mur-2D0*Ali/MW*Ml/SW*Zh31*Zl33i*Zl63i/cb*sb
     &      -2D0*Ali/MW*Ml/SW*Zh31*Zl33r*Zl63r/cb*sb-2D0
     &      *Alr/MW*Ml/SW*Zh31*Zl33i*Zl63r/cb*sb+2D0*Alr
     &      /MW*Ml/SW*Zh31*Zl33r*Zl63i/cb*sb)
      AAABR(74) = 2D0*Sqrt2*Zh12*Zl11*Zl41*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl11**2/cb+EE
     &      *MW/SW*Zh22*Zl11**2*sb+EE/CW**2*MW*SW*Zh12*Zl11
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl41**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl11**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl41**2*sb+EE*MW/SW*Zh12*Zl11**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl11**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl41**2/cb*sb**2+2D0*EE
     &      /MW*Me/SW*Zh32*Zl11*Zl41*mui-2D0*EE/MW*Me/SW
     &      *Zh22*Zl11*Zl41/cb*mur
      AAABR(75) = EE*(2D0*Ml**2/MW/SW*Zh12/cb-MW/SW*Zh12
     &      *Zl33i**2/cb-MW/SW*Zh12*Zl33r**2/cb+MW/SW*Zh22
     &      *Zl33i**2*sb+MW/SW*Zh22*Zl33r**2*sb+MW/CW**2
     &      *SW*Zh12*Zl33i**2/cb+MW/CW**2*SW*Zh12*Zl33r**2
     &      /cb-2D0*MW/CW**2*SW*Zh12*Zl63i**2/cb-2D0*MW/CW
     &      **2*SW*Zh12*Zl63r**2/cb-MW/CW**2*SW*Zh22*Zl33i
     &      **2*sb-MW/CW**2*SW*Zh22*Zl33r**2*sb+2D0*MW/CW
     &      **2*SW*Zh22*Zl63i**2*sb+2D0*MW/CW**2*SW*Zh22
     &      *Zl63r**2*sb+MW/SW*Zh12*Zl33i**2/cb*sb**2+MW
     &      /SW*Zh12*Zl33r**2/cb*sb**2-MW/CW**2*SW*Zh12*Zl33i
     &      **2/cb*sb**2-MW/CW**2*SW*Zh12*Zl33r**2/cb*sb
     &      **2+2D0*MW/CW**2*SW*Zh12*Zl63i**2/cb*sb**2+2D0
     &      *MW/CW**2*SW*Zh12*Zl63r**2/cb*sb**2+2D0*Ml/MW
     &      /SW*Zh32*Zl33i*Zl63i*mui-2D0*Ml/MW/SW*Zh32*Zl33i
     &      *Zl63r*mur+2D0*Ml/MW/SW*Zh32*Zl33r*Zl63i*mur
     &      +2D0*Ml/MW/SW*Zh32*Zl33r*Zl63r*mui-2D0*Ali/MW
     &      *Ml/SW*Zh12*Zl33i*Zl63r/cb+2D0*Ali/MW*Ml/SW*Zh12
     &      *Zl33r*Zl63i/cb+2D0*Alr/MW*Ml/SW*Zh12*Zl33i*Zl63i
     &      /cb+2D0*Alr/MW*Ml/SW*Zh12*Zl33r*Zl63r/cb-2D0
     &      *Ml/MW/SW*Zh22*Zl33i*Zl63i/cb*mur-2D0*Ml/MW/SW
     &      *Zh22*Zl33i*Zl63r/cb*mui+2D0*Ml/MW/SW*Zh22*Zl33r
     &      *Zl63i/cb*mui-2D0*Ml/MW/SW*Zh22*Zl33r*Zl63r/cb
     &      *mur-2D0*Ali/MW*Ml/SW*Zh32*Zl33i*Zl63i/cb*sb
     &      -2D0*Ali/MW*Ml/SW*Zh32*Zl33r*Zl63r/cb*sb-2D0
     &      *Alr/MW*Ml/SW*Zh32*Zl33i*Zl63r/cb*sb+2D0*Alr
     &      /MW*Ml/SW*Zh32*Zl33r*Zl63i/cb*sb)
      AAABR(76) = 2D0*Sqrt2*Zh11*Zl22*Zl52*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl22**2/cb+EE
     &      *MW/SW*Zh21*Zl22**2*sb+EE/CW**2*MW*SW*Zh11*Zl22
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl52**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl22**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl52**2*sb+EE*MW/SW*Zh11*Zl22**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl22**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl52**2/cb*sb**2+2D0*EE
     &      /MW*Mm/SW*Zh31*Zl22*Zl52*mui-2D0*EE/MW*Mm/SW
     &      *Zh21*Zl22*Zl52/cb*mur
      AAABR(77) = 2D0*Sqrt2*Zh13*Zl22*Zl52*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl22**2/cb+EE
     &      *MW/SW*Zh23*Zl22**2*sb+EE/CW**2*MW*SW*Zh13*Zl22
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl52**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl22**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl52**2*sb+EE*MW/SW*Zh13*Zl22**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl22**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl52**2/cb*sb**2+2D0*EE
     &      /MW*Mm/SW*Zh33*Zl22*Zl52*mui-2D0*EE/MW*Mm/SW
     &      *Zh23*Zl22*Zl52/cb*mur
      AAABR(78) = 2D0*Sqrt2*Zh13*Zl11*Zl41*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl11**2/cb+EE
     &      *MW/SW*Zh23*Zl11**2*sb+EE/CW**2*MW*SW*Zh13*Zl11
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl41**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl11**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl41**2*sb+EE*MW/SW*Zh13*Zl11**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl11**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl41**2/cb*sb**2+2D0*EE
     &      /MW*Me/SW*Zh33*Zl11*Zl41*mui-2D0*EE/MW*Me/SW
     &      *Zh23*Zl11*Zl41/cb*mur
      AAABR(79) = 2D0*Sqrt2*Zh12*Zl22*Zl52*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl22**2/cb+EE
     &      *MW/SW*Zh22*Zl22**2*sb+EE/CW**2*MW*SW*Zh12*Zl22
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl52**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl22**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl52**2*sb+EE*MW/SW*Zh12*Zl22**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl22**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl52**2/cb*sb**2+2D0*EE
     &      /MW*Mm/SW*Zh32*Zl22*Zl52*mui-2D0*EE/MW*Mm/SW
     &      *Zh22*Zl22*Zl52/cb*mur
      AAABR(80) = EE*(2D0*Ml**2/MW/SW*Zh13/cb-MW/SW*Zh13
     &      *Zl33i**2/cb-MW/SW*Zh13*Zl33r**2/cb+MW/SW*Zh23
     &      *Zl33i**2*sb+MW/SW*Zh23*Zl33r**2*sb+MW/CW**2
     &      *SW*Zh13*Zl33i**2/cb+MW/CW**2*SW*Zh13*Zl33r**2
     &      /cb-2D0*MW/CW**2*SW*Zh13*Zl63i**2/cb-2D0*MW/CW
     &      **2*SW*Zh13*Zl63r**2/cb-MW/CW**2*SW*Zh23*Zl33i
     &      **2*sb-MW/CW**2*SW*Zh23*Zl33r**2*sb+2D0*MW/CW
     &      **2*SW*Zh23*Zl63i**2*sb+2D0*MW/CW**2*SW*Zh23
     &      *Zl63r**2*sb+MW/SW*Zh13*Zl33i**2/cb*sb**2+MW
     &      /SW*Zh13*Zl33r**2/cb*sb**2-MW/CW**2*SW*Zh13*Zl33i
     &      **2/cb*sb**2-MW/CW**2*SW*Zh13*Zl33r**2/cb*sb
     &      **2+2D0*MW/CW**2*SW*Zh13*Zl63i**2/cb*sb**2+2D0
     &      *MW/CW**2*SW*Zh13*Zl63r**2/cb*sb**2+2D0*Ml/MW
     &      /SW*Zh33*Zl33i*Zl63i*mui-2D0*Ml/MW/SW*Zh33*Zl33i
     &      *Zl63r*mur+2D0*Ml/MW/SW*Zh33*Zl33r*Zl63i*mur
     &      +2D0*Ml/MW/SW*Zh33*Zl33r*Zl63r*mui-2D0*Ali/MW
     &      *Ml/SW*Zh13*Zl33i*Zl63r/cb+2D0*Ali/MW*Ml/SW*Zh13
     &      *Zl33r*Zl63i/cb+2D0*Alr/MW*Ml/SW*Zh13*Zl33i*Zl63i
     &      /cb+2D0*Alr/MW*Ml/SW*Zh13*Zl33r*Zl63r/cb-2D0
     &      *Ml/MW/SW*Zh23*Zl33i*Zl63i/cb*mur-2D0*Ml/MW/SW
     &      *Zh23*Zl33i*Zl63r/cb*mui+2D0*Ml/MW/SW*Zh23*Zl33r
     &      *Zl63i/cb*mui-2D0*Ml/MW/SW*Zh23*Zl33r*Zl63r/cb
     &      *mur-2D0*Ali/MW*Ml/SW*Zh33*Zl33i*Zl63i/cb*sb
     &      -2D0*Ali/MW*Ml/SW*Zh33*Zl33r*Zl63r/cb*sb-2D0
     &      *Alr/MW*Ml/SW*Zh33*Zl33i*Zl63r/cb*sb+2D0*Alr
     &      /MW*Ml/SW*Zh33*Zl33r*Zl63i/cb*sb)
      AAABR(81) = Sqrt2*Zh11*Zl22*Zl55*ls2+Sqrt2*Zh11*Zl25
     &      *Zl52*ls2-EE*MW/SW*Zh11*Zl22*Zl25/cb+EE*MW/SW
     &      *Zh21*Zl22*Zl25*sb+EE/CW**2*MW*SW*Zh11*Zl22*Zl25
     &      /cb-2D0*EE/CW**2*MW*SW*Zh11*Zl52*Zl55/cb-EE/CW
     &      **2*MW*SW*Zh21*Zl22*Zl25*sb+2D0*EE/CW**2*MW*SW
     &      *Zh21*Zl52*Zl55*sb+EE/MW*Mm/SW*Zh31*Zl22*Zl55
     &      *mui+EE/MW*Mm/SW*Zh31*Zl25*Zl52*mui+EE*MW/SW
     &      *Zh11*Zl22*Zl25/cb*sb**2-EE/CW**2*MW*SW*Zh11
     &      *Zl22*Zl25/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl52
     &      *Zl55/cb*sb**2-EE/MW*Mm/SW*Zh21*Zl22*Zl55/cb
     &      *mur-EE/MW*Mm/SW*Zh21*Zl25*Zl52/cb*mur
      AAABR(82) = Sqrt2*Zh31*Zl22*Zl55*ls2*sb-Sqrt2*Zh31
     &      *Zl25*Zl52*ls2*sb+EE/MW*Mm/SW*Zh31*Zl22*Zl55
     &      *mur-EE/MW*Mm/SW*Zh31*Zl25*Zl52*mur+EE/MW*Mm
     &      /SW*Zh21*Zl22*Zl55/cb*mui-EE/MW*Mm/SW*Zh21*Zl25
     &      *Zl52/cb*mui
      AAABR(83) = EE*(MW/SW*Zh11*Zl33i*Zl36i/cb+MW/SW*Zh11
     &      *Zl33r*Zl36r/cb-MW/SW*Zh21*Zl33i*Zl36i*sb-MW
     &      /SW*Zh21*Zl33r*Zl36r*sb-MW/CW**2*SW*Zh11*Zl33i
     &      *Zl36i/cb-MW/CW**2*SW*Zh11*Zl33r*Zl36r/cb+2D0
     &      *MW/CW**2*SW*Zh11*Zl63i*Zl66i/cb+2D0*MW/CW**2
     &      *SW*Zh11*Zl63r*Zl66r/cb+MW/CW**2*SW*Zh21*Zl33i
     &      *Zl36i*sb+MW/CW**2*SW*Zh21*Zl33r*Zl36r*sb-2D0
     &      *MW/CW**2*SW*Zh21*Zl63i*Zl66i*sb-2D0*MW/CW**2
     &      *SW*Zh21*Zl63r*Zl66r*sb-Ml/MW/SW*Zh31*Zl33i*Zl66i
     &      *mui+Ml/MW/SW*Zh31*Zl33i*Zl66r*mur-Ml/MW/SW*Zh31
     &      *Zl33r*Zl66i*mur-Ml/MW/SW*Zh31*Zl33r*Zl66r*mui
     &      -Ml/MW/SW*Zh31*Zl36i*Zl63i*mui+Ml/MW/SW*Zh31
     &      *Zl36i*Zl63r*mur-Ml/MW/SW*Zh31*Zl36r*Zl63i*mur
     &      -Ml/MW/SW*Zh31*Zl36r*Zl63r*mui-MW/SW*Zh11*Zl33i
     &      *Zl36i/cb*sb**2-MW/SW*Zh11*Zl33r*Zl36r/cb*sb
     &      **2+Ali/MW*Ml/SW*Zh11*Zl33i*Zl66r/cb-Ali/MW*Ml
     &      /SW*Zh11*Zl33r*Zl66i/cb+Ali/MW*Ml/SW*Zh11*Zl36i
     &      *Zl63r/cb-Ali/MW*Ml/SW*Zh11*Zl36r*Zl63i/cb-Alr
     &      /MW*Ml/SW*Zh11*Zl33i*Zl66i/cb-Alr/MW*Ml/SW*Zh11
     &      *Zl33r*Zl66r/cb-Alr/MW*Ml/SW*Zh11*Zl36i*Zl63i
     &      /cb-Alr/MW*Ml/SW*Zh11*Zl36r*Zl63r/cb+MW/CW**2
     &      *SW*Zh11*Zl33i*Zl36i/cb*sb**2+MW/CW**2*SW*Zh11
     &      *Zl33r*Zl36r/cb*sb**2-2D0*MW/CW**2*SW*Zh11*Zl63i
     &      *Zl66i/cb*sb**2-2D0*MW/CW**2*SW*Zh11*Zl63r*Zl66r
     &      /cb*sb**2+Ml/MW/SW*Zh21*Zl33i*Zl66i/cb*mur+Ml
     &      /MW/SW*Zh21*Zl33i*Zl66r/cb*mui-Ml/MW/SW*Zh21
     &      *Zl33r*Zl66i/cb*mui+Ml/MW/SW*Zh21*Zl33r*Zl66r
     &      /cb*mur+Ml/MW/SW*Zh21*Zl36i*Zl63i/cb*mur+Ml/MW
     &      /SW*Zh21*Zl36i*Zl63r/cb*mui-Ml/MW/SW*Zh21*Zl36r
     &      *Zl63i/cb*mui+Ml/MW/SW*Zh21*Zl36r*Zl63r/cb*mur
     &      +Ali/MW*Ml/SW*Zh31*Zl33i*Zl66i/cb*sb+Ali/MW*Ml
     &      /SW*Zh31*Zl33r*Zl66r/cb*sb+Ali/MW*Ml/SW*Zh31
     &      *Zl36i*Zl63i/cb*sb+Ali/MW*Ml/SW*Zh31*Zl36r*Zl63r
     &      /cb*sb+Alr/MW*Ml/SW*Zh31*Zl33i*Zl66r/cb*sb-Alr
     &      /MW*Ml/SW*Zh31*Zl33r*Zl66i/cb*sb+Alr/MW*Ml/SW
     &      *Zh31*Zl36i*Zl63r/cb*sb-Alr/MW*Ml/SW*Zh31*Zl36r
     &      *Zl63i/cb*sb)
      AAABR(84) = EE*(MW/SW*Zh11*Zl33i*Zl36r/cb-MW/SW*Zh11
     &      *Zl33r*Zl36i/cb-MW/SW*Zh21*Zl33i*Zl36r*sb+MW
     &      /SW*Zh21*Zl33r*Zl36i*sb-MW/CW**2*SW*Zh11*Zl33i
     &      *Zl36r/cb+MW/CW**2*SW*Zh11*Zl33r*Zl36i/cb+2D0
     &      *MW/CW**2*SW*Zh11*Zl63i*Zl66r/cb-2D0*MW/CW**2
     &      *SW*Zh11*Zl63r*Zl66i/cb+MW/CW**2*SW*Zh21*Zl33i
     &      *Zl36r*sb-MW/CW**2*SW*Zh21*Zl33r*Zl36i*sb-2D0
     &      *MW/CW**2*SW*Zh21*Zl63i*Zl66r*sb+2D0*MW/CW**2
     &      *SW*Zh21*Zl63r*Zl66i*sb-Ml/MW/SW*Zh31*Zl33i*Zl66i
     &      *mur-Ml/MW/SW*Zh31*Zl33i*Zl66r*mui+Ml/MW/SW*Zh31
     &      *Zl33r*Zl66i*mui-Ml/MW/SW*Zh31*Zl33r*Zl66r*mur
     &      +Ml/MW/SW*Zh31*Zl36i*Zl63i*mur+Ml/MW/SW*Zh31
     &      *Zl36i*Zl63r*mui-Ml/MW/SW*Zh31*Zl36r*Zl63i*mui
     &      +Ml/MW/SW*Zh31*Zl36r*Zl63r*mur-MW/SW*Zh11*Zl33i
     &      *Zl36r/cb*sb**2+MW/SW*Zh11*Zl33r*Zl36i/cb*sb
     &      **2-Ali/MW*Ml/SW*Zh11*Zl33i*Zl66i/cb-Ali/MW*Ml
     &      /SW*Zh11*Zl33r*Zl66r/cb+Ali/MW*Ml/SW*Zh11*Zl36i
     &      *Zl63i/cb+Ali/MW*Ml/SW*Zh11*Zl36r*Zl63r/cb-Alr
     &      /MW*Ml/SW*Zh11*Zl33i*Zl66r/cb+Alr/MW*Ml/SW*Zh11
     &      *Zl33r*Zl66i/cb+Alr/MW*Ml/SW*Zh11*Zl36i*Zl63r
     &      /cb-Alr/MW*Ml/SW*Zh11*Zl36r*Zl63i/cb+MW/CW**2
     &      *SW*Zh11*Zl33i*Zl36r/cb*sb**2-MW/CW**2*SW*Zh11
     &      *Zl33r*Zl36i/cb*sb**2-2D0*MW/CW**2*SW*Zh11*Zl63i
     &      *Zl66r/cb*sb**2+2D0*MW/CW**2*SW*Zh11*Zl63r*Zl66i
     &      /cb*sb**2-Ml/MW/SW*Zh21*Zl33i*Zl66i/cb*mui+Ml
     &      /MW/SW*Zh21*Zl33i*Zl66r/cb*mur-Ml/MW/SW*Zh21
     &      *Zl33r*Zl66i/cb*mur-Ml/MW/SW*Zh21*Zl33r*Zl66r
     &      /cb*mui+Ml/MW/SW*Zh21*Zl36i*Zl63i/cb*mui-Ml/MW
     &      /SW*Zh21*Zl36i*Zl63r/cb*mur+Ml/MW/SW*Zh21*Zl36r
     &      *Zl63i/cb*mur+Ml/MW/SW*Zh21*Zl36r*Zl63r/cb*mui
     &      +Ali/MW*Ml/SW*Zh31*Zl33i*Zl66r/cb*sb-Ali/MW*Ml
     &      /SW*Zh31*Zl33r*Zl66i/cb*sb-Ali/MW*Ml/SW*Zh31
     &      *Zl36i*Zl63r/cb*sb+Ali/MW*Ml/SW*Zh31*Zl36r*Zl63i
     &      /cb*sb-Alr/MW*Ml/SW*Zh31*Zl33i*Zl66i/cb*sb-Alr
     &      /MW*Ml/SW*Zh31*Zl33r*Zl66r/cb*sb+Alr/MW*Ml/SW
     &      *Zh31*Zl36i*Zl63i/cb*sb+Alr/MW*Ml/SW*Zh31*Zl36r
     &      *Zl63r/cb*sb)
      AAABR(85) = EE*(MW/SW*Zh13*Zl33i*Zl36i/cb+MW/SW*Zh13
     &      *Zl33r*Zl36r/cb-MW/SW*Zh23*Zl33i*Zl36i*sb-MW
     &      /SW*Zh23*Zl33r*Zl36r*sb-MW/CW**2*SW*Zh13*Zl33i
     &      *Zl36i/cb-MW/CW**2*SW*Zh13*Zl33r*Zl36r/cb+2D0
     &      *MW/CW**2*SW*Zh13*Zl63i*Zl66i/cb+2D0*MW/CW**2
     &      *SW*Zh13*Zl63r*Zl66r/cb+MW/CW**2*SW*Zh23*Zl33i
     &      *Zl36i*sb+MW/CW**2*SW*Zh23*Zl33r*Zl36r*sb-2D0
     &      *MW/CW**2*SW*Zh23*Zl63i*Zl66i*sb-2D0*MW/CW**2
     &      *SW*Zh23*Zl63r*Zl66r*sb-Ml/MW/SW*Zh33*Zl33i*Zl66i
     &      *mui+Ml/MW/SW*Zh33*Zl33i*Zl66r*mur-Ml/MW/SW*Zh33
     &      *Zl33r*Zl66i*mur-Ml/MW/SW*Zh33*Zl33r*Zl66r*mui
     &      -Ml/MW/SW*Zh33*Zl36i*Zl63i*mui+Ml/MW/SW*Zh33
     &      *Zl36i*Zl63r*mur-Ml/MW/SW*Zh33*Zl36r*Zl63i*mur
     &      -Ml/MW/SW*Zh33*Zl36r*Zl63r*mui-MW/SW*Zh13*Zl33i
     &      *Zl36i/cb*sb**2-MW/SW*Zh13*Zl33r*Zl36r/cb*sb
     &      **2+Ali/MW*Ml/SW*Zh13*Zl33i*Zl66r/cb-Ali/MW*Ml
     &      /SW*Zh13*Zl33r*Zl66i/cb+Ali/MW*Ml/SW*Zh13*Zl36i
     &      *Zl63r/cb-Ali/MW*Ml/SW*Zh13*Zl36r*Zl63i/cb-Alr
     &      /MW*Ml/SW*Zh13*Zl33i*Zl66i/cb-Alr/MW*Ml/SW*Zh13
     &      *Zl33r*Zl66r/cb-Alr/MW*Ml/SW*Zh13*Zl36i*Zl63i
     &      /cb-Alr/MW*Ml/SW*Zh13*Zl36r*Zl63r/cb+MW/CW**2
     &      *SW*Zh13*Zl33i*Zl36i/cb*sb**2+MW/CW**2*SW*Zh13
     &      *Zl33r*Zl36r/cb*sb**2-2D0*MW/CW**2*SW*Zh13*Zl63i
     &      *Zl66i/cb*sb**2-2D0*MW/CW**2*SW*Zh13*Zl63r*Zl66r
     &      /cb*sb**2+Ml/MW/SW*Zh23*Zl33i*Zl66i/cb*mur+Ml
     &      /MW/SW*Zh23*Zl33i*Zl66r/cb*mui-Ml/MW/SW*Zh23
     &      *Zl33r*Zl66i/cb*mui+Ml/MW/SW*Zh23*Zl33r*Zl66r
     &      /cb*mur+Ml/MW/SW*Zh23*Zl36i*Zl63i/cb*mur+Ml/MW
     &      /SW*Zh23*Zl36i*Zl63r/cb*mui-Ml/MW/SW*Zh23*Zl36r
     &      *Zl63i/cb*mui+Ml/MW/SW*Zh23*Zl36r*Zl63r/cb*mur
     &      +Ali/MW*Ml/SW*Zh33*Zl33i*Zl66i/cb*sb+Ali/MW*Ml
     &      /SW*Zh33*Zl33r*Zl66r/cb*sb+Ali/MW*Ml/SW*Zh33
     &      *Zl36i*Zl63i/cb*sb+Ali/MW*Ml/SW*Zh33*Zl36r*Zl63r
     &      /cb*sb+Alr/MW*Ml/SW*Zh33*Zl33i*Zl66r/cb*sb-Alr
     &      /MW*Ml/SW*Zh33*Zl33r*Zl66i/cb*sb+Alr/MW*Ml/SW
     &      *Zh33*Zl36i*Zl63r/cb*sb-Alr/MW*Ml/SW*Zh33*Zl36r
     &      *Zl63i/cb*sb)
      AAABR(86) = EE*(MW/SW*Zh13*Zl33i*Zl36r/cb-MW/SW*Zh13
     &      *Zl33r*Zl36i/cb-MW/SW*Zh23*Zl33i*Zl36r*sb+MW
     &      /SW*Zh23*Zl33r*Zl36i*sb-MW/CW**2*SW*Zh13*Zl33i
     &      *Zl36r/cb+MW/CW**2*SW*Zh13*Zl33r*Zl36i/cb+2D0
     &      *MW/CW**2*SW*Zh13*Zl63i*Zl66r/cb-2D0*MW/CW**2
     &      *SW*Zh13*Zl63r*Zl66i/cb+MW/CW**2*SW*Zh23*Zl33i
     &      *Zl36r*sb-MW/CW**2*SW*Zh23*Zl33r*Zl36i*sb-2D0
     &      *MW/CW**2*SW*Zh23*Zl63i*Zl66r*sb+2D0*MW/CW**2
     &      *SW*Zh23*Zl63r*Zl66i*sb-Ml/MW/SW*Zh33*Zl33i*Zl66i
     &      *mur-Ml/MW/SW*Zh33*Zl33i*Zl66r*mui+Ml/MW/SW*Zh33
     &      *Zl33r*Zl66i*mui-Ml/MW/SW*Zh33*Zl33r*Zl66r*mur
     &      +Ml/MW/SW*Zh33*Zl36i*Zl63i*mur+Ml/MW/SW*Zh33
     &      *Zl36i*Zl63r*mui-Ml/MW/SW*Zh33*Zl36r*Zl63i*mui
     &      +Ml/MW/SW*Zh33*Zl36r*Zl63r*mur-MW/SW*Zh13*Zl33i
     &      *Zl36r/cb*sb**2+MW/SW*Zh13*Zl33r*Zl36i/cb*sb
     &      **2-Ali/MW*Ml/SW*Zh13*Zl33i*Zl66i/cb-Ali/MW*Ml
     &      /SW*Zh13*Zl33r*Zl66r/cb+Ali/MW*Ml/SW*Zh13*Zl36i
     &      *Zl63i/cb+Ali/MW*Ml/SW*Zh13*Zl36r*Zl63r/cb-Alr
     &      /MW*Ml/SW*Zh13*Zl33i*Zl66r/cb+Alr/MW*Ml/SW*Zh13
     &      *Zl33r*Zl66i/cb+Alr/MW*Ml/SW*Zh13*Zl36i*Zl63r
     &      /cb-Alr/MW*Ml/SW*Zh13*Zl36r*Zl63i/cb+MW/CW**2
     &      *SW*Zh13*Zl33i*Zl36r/cb*sb**2-MW/CW**2*SW*Zh13
     &      *Zl33r*Zl36i/cb*sb**2-2D0*MW/CW**2*SW*Zh13*Zl63i
     &      *Zl66r/cb*sb**2+2D0*MW/CW**2*SW*Zh13*Zl63r*Zl66i
     &      /cb*sb**2-Ml/MW/SW*Zh23*Zl33i*Zl66i/cb*mui+Ml
     &      /MW/SW*Zh23*Zl33i*Zl66r/cb*mur-Ml/MW/SW*Zh23
     &      *Zl33r*Zl66i/cb*mur-Ml/MW/SW*Zh23*Zl33r*Zl66r
     &      /cb*mui+Ml/MW/SW*Zh23*Zl36i*Zl63i/cb*mui-Ml/MW
     &      /SW*Zh23*Zl36i*Zl63r/cb*mur+Ml/MW/SW*Zh23*Zl36r
     &      *Zl63i/cb*mur+Ml/MW/SW*Zh23*Zl36r*Zl63r/cb*mui
     &      +Ali/MW*Ml/SW*Zh33*Zl33i*Zl66r/cb*sb-Ali/MW*Ml
     &      /SW*Zh33*Zl33r*Zl66i/cb*sb-Ali/MW*Ml/SW*Zh33
     &      *Zl36i*Zl63r/cb*sb+Ali/MW*Ml/SW*Zh33*Zl36r*Zl63i
     &      /cb*sb-Alr/MW*Ml/SW*Zh33*Zl33i*Zl66i/cb*sb-Alr
     &      /MW*Ml/SW*Zh33*Zl33r*Zl66r/cb*sb+Alr/MW*Ml/SW
     &      *Zh33*Zl36i*Zl63i/cb*sb+Alr/MW*Ml/SW*Zh33*Zl36r
     &      *Zl63r/cb*sb)
      AAABR(87) = Sqrt2*Zh13*Zl22*Zl55*ls2+Sqrt2*Zh13*Zl25
     &      *Zl52*ls2-EE*MW/SW*Zh13*Zl22*Zl25/cb+EE*MW/SW
     &      *Zh23*Zl22*Zl25*sb+EE/CW**2*MW*SW*Zh13*Zl22*Zl25
     &      /cb-2D0*EE/CW**2*MW*SW*Zh13*Zl52*Zl55/cb-EE/CW
     &      **2*MW*SW*Zh23*Zl22*Zl25*sb+2D0*EE/CW**2*MW*SW
     &      *Zh23*Zl52*Zl55*sb+EE/MW*Mm/SW*Zh33*Zl22*Zl55
     &      *mui+EE/MW*Mm/SW*Zh33*Zl25*Zl52*mui+EE*MW/SW
     &      *Zh13*Zl22*Zl25/cb*sb**2-EE/CW**2*MW*SW*Zh13
     &      *Zl22*Zl25/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh13*Zl52
     &      *Zl55/cb*sb**2-EE/MW*Mm/SW*Zh23*Zl22*Zl55/cb
     &      *mur-EE/MW*Mm/SW*Zh23*Zl25*Zl52/cb*mur
      AAABR(88) = Sqrt2*Zh33*Zl22*Zl55*ls2*sb-Sqrt2*Zh33
     &      *Zl25*Zl52*ls2*sb+EE/MW*Mm/SW*Zh33*Zl22*Zl55
     &      *mur-EE/MW*Mm/SW*Zh33*Zl25*Zl52*mur+EE/MW*Mm
     &      /SW*Zh23*Zl22*Zl55/cb*mui-EE/MW*Mm/SW*Zh23*Zl25
     &      *Zl52/cb*mui
      AAABR(89) = Sqrt2*Zh13*Zl11*Zl44*ls1+Sqrt2*Zh13*Zl14
     &      *Zl41*ls1-EE*MW/SW*Zh13*Zl11*Zl14/cb+EE*MW/SW
     &      *Zh23*Zl11*Zl14*sb+EE/CW**2*MW*SW*Zh13*Zl11*Zl14
     &      /cb-2D0*EE/CW**2*MW*SW*Zh13*Zl41*Zl44/cb-EE/CW
     &      **2*MW*SW*Zh23*Zl11*Zl14*sb+2D0*EE/CW**2*MW*SW
     &      *Zh23*Zl41*Zl44*sb+EE/MW*Me/SW*Zh33*Zl11*Zl44
     &      *mui+EE/MW*Me/SW*Zh33*Zl14*Zl41*mui+EE*MW/SW
     &      *Zh13*Zl11*Zl14/cb*sb**2-EE/CW**2*MW*SW*Zh13
     &      *Zl11*Zl14/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh13*Zl41
     &      *Zl44/cb*sb**2-EE/MW*Me/SW*Zh23*Zl11*Zl44/cb
     &      *mur-EE/MW*Me/SW*Zh23*Zl14*Zl41/cb*mur
      AAABR(90) = Sqrt2*Zh33*Zl11*Zl44*ls1*sb-Sqrt2*Zh33
     &      *Zl14*Zl41*ls1*sb+EE/MW*Me/SW*Zh33*Zl11*Zl44
     &      *mur-EE/MW*Me/SW*Zh33*Zl14*Zl41*mur+EE/MW*Me
     &      /SW*Zh23*Zl11*Zl44/cb*mui-EE/MW*Me/SW*Zh23*Zl14
     &      *Zl41/cb*mui
      AAABR(91) = Sqrt2*Zh11*Zl11*Zl44*ls1+Sqrt2*Zh11*Zl14
     &      *Zl41*ls1-EE*MW/SW*Zh11*Zl11*Zl14/cb+EE*MW/SW
     &      *Zh21*Zl11*Zl14*sb+EE/CW**2*MW*SW*Zh11*Zl11*Zl14
     &      /cb-2D0*EE/CW**2*MW*SW*Zh11*Zl41*Zl44/cb-EE/CW
     &      **2*MW*SW*Zh21*Zl11*Zl14*sb+2D0*EE/CW**2*MW*SW
     &      *Zh21*Zl41*Zl44*sb+EE/MW*Me/SW*Zh31*Zl11*Zl44
     &      *mui+EE/MW*Me/SW*Zh31*Zl14*Zl41*mui+EE*MW/SW
     &      *Zh11*Zl11*Zl14/cb*sb**2-EE/CW**2*MW*SW*Zh11
     &      *Zl11*Zl14/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl41
     &      *Zl44/cb*sb**2-EE/MW*Me/SW*Zh21*Zl11*Zl44/cb
     &      *mur-EE/MW*Me/SW*Zh21*Zl14*Zl41/cb*mur
      AAABR(92) = Sqrt2*Zh31*Zl11*Zl44*ls1*sb-Sqrt2*Zh31
     &      *Zl14*Zl41*ls1*sb+EE/MW*Me/SW*Zh31*Zl11*Zl44
     &      *mur-EE/MW*Me/SW*Zh31*Zl14*Zl41*mur+EE/MW*Me
     &      /SW*Zh21*Zl11*Zl44/cb*mui-EE/MW*Me/SW*Zh21*Zl14
     &      *Zl41/cb*mui
      AAABR(93) = Sqrt2*Zh12*Zl22*Zl55*ls2+Sqrt2*Zh12*Zl25
     &      *Zl52*ls2-EE*MW/SW*Zh12*Zl22*Zl25/cb+EE*MW/SW
     &      *Zh22*Zl22*Zl25*sb+EE/CW**2*MW*SW*Zh12*Zl22*Zl25
     &      /cb-2D0*EE/CW**2*MW*SW*Zh12*Zl52*Zl55/cb-EE/CW
     &      **2*MW*SW*Zh22*Zl22*Zl25*sb+2D0*EE/CW**2*MW*SW
     &      *Zh22*Zl52*Zl55*sb+EE/MW*Mm/SW*Zh32*Zl22*Zl55
     &      *mui+EE/MW*Mm/SW*Zh32*Zl25*Zl52*mui+EE*MW/SW
     &      *Zh12*Zl22*Zl25/cb*sb**2-EE/CW**2*MW*SW*Zh12
     &      *Zl22*Zl25/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh12*Zl52
     &      *Zl55/cb*sb**2-EE/MW*Mm/SW*Zh22*Zl22*Zl55/cb
     &      *mur-EE/MW*Mm/SW*Zh22*Zl25*Zl52/cb*mur
      AAABR(94) = Sqrt2*Zh32*Zl22*Zl55*ls2*sb-Sqrt2*Zh32
     &      *Zl25*Zl52*ls2*sb+EE/MW*Mm/SW*Zh32*Zl22*Zl55
     &      *mur-EE/MW*Mm/SW*Zh32*Zl25*Zl52*mur+EE/MW*Mm
     &      /SW*Zh22*Zl22*Zl55/cb*mui-EE/MW*Mm/SW*Zh22*Zl25
     &      *Zl52/cb*mui
      AAABR(95) = Sqrt2*Zh12*Zl11*Zl44*ls1+Sqrt2*Zh12*Zl14
     &      *Zl41*ls1-EE*MW/SW*Zh12*Zl11*Zl14/cb+EE*MW/SW
     &      *Zh22*Zl11*Zl14*sb+EE/CW**2*MW*SW*Zh12*Zl11*Zl14
     &      /cb-2D0*EE/CW**2*MW*SW*Zh12*Zl41*Zl44/cb-EE/CW
     &      **2*MW*SW*Zh22*Zl11*Zl14*sb+2D0*EE/CW**2*MW*SW
     &      *Zh22*Zl41*Zl44*sb+EE/MW*Me/SW*Zh32*Zl11*Zl44
     &      *mui+EE/MW*Me/SW*Zh32*Zl14*Zl41*mui+EE*MW/SW
     &      *Zh12*Zl11*Zl14/cb*sb**2-EE/CW**2*MW*SW*Zh12
     &      *Zl11*Zl14/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh12*Zl41
     &      *Zl44/cb*sb**2-EE/MW*Me/SW*Zh22*Zl11*Zl44/cb
     &      *mur-EE/MW*Me/SW*Zh22*Zl14*Zl41/cb*mur
      AAABR(96) = Sqrt2*Zh32*Zl11*Zl44*ls1*sb-Sqrt2*Zh32
     &      *Zl14*Zl41*ls1*sb+EE/MW*Me/SW*Zh32*Zl11*Zl44
     &      *mur-EE/MW*Me/SW*Zh32*Zl14*Zl41*mur+EE/MW*Me
     &      /SW*Zh22*Zl11*Zl44/cb*mui-EE/MW*Me/SW*Zh22*Zl14
     &      *Zl41/cb*mui
      AAABR(97) = EE*(MW/SW*Zh12*Zl33i*Zl36i/cb+MW/SW*Zh12
     &      *Zl33r*Zl36r/cb-MW/SW*Zh22*Zl33i*Zl36i*sb-MW
     &      /SW*Zh22*Zl33r*Zl36r*sb-MW/CW**2*SW*Zh12*Zl33i
     &      *Zl36i/cb-MW/CW**2*SW*Zh12*Zl33r*Zl36r/cb+2D0
     &      *MW/CW**2*SW*Zh12*Zl63i*Zl66i/cb+2D0*MW/CW**2
     &      *SW*Zh12*Zl63r*Zl66r/cb+MW/CW**2*SW*Zh22*Zl33i
     &      *Zl36i*sb+MW/CW**2*SW*Zh22*Zl33r*Zl36r*sb-2D0
     &      *MW/CW**2*SW*Zh22*Zl63i*Zl66i*sb-2D0*MW/CW**2
     &      *SW*Zh22*Zl63r*Zl66r*sb-Ml/MW/SW*Zh32*Zl33i*Zl66i
     &      *mui+Ml/MW/SW*Zh32*Zl33i*Zl66r*mur-Ml/MW/SW*Zh32
     &      *Zl33r*Zl66i*mur-Ml/MW/SW*Zh32*Zl33r*Zl66r*mui
     &      -Ml/MW/SW*Zh32*Zl36i*Zl63i*mui+Ml/MW/SW*Zh32
     &      *Zl36i*Zl63r*mur-Ml/MW/SW*Zh32*Zl36r*Zl63i*mur
     &      -Ml/MW/SW*Zh32*Zl36r*Zl63r*mui-MW/SW*Zh12*Zl33i
     &      *Zl36i/cb*sb**2-MW/SW*Zh12*Zl33r*Zl36r/cb*sb
     &      **2+Ali/MW*Ml/SW*Zh12*Zl33i*Zl66r/cb-Ali/MW*Ml
     &      /SW*Zh12*Zl33r*Zl66i/cb+Ali/MW*Ml/SW*Zh12*Zl36i
     &      *Zl63r/cb-Ali/MW*Ml/SW*Zh12*Zl36r*Zl63i/cb-Alr
     &      /MW*Ml/SW*Zh12*Zl33i*Zl66i/cb-Alr/MW*Ml/SW*Zh12
     &      *Zl33r*Zl66r/cb-Alr/MW*Ml/SW*Zh12*Zl36i*Zl63i
     &      /cb-Alr/MW*Ml/SW*Zh12*Zl36r*Zl63r/cb+MW/CW**2
     &      *SW*Zh12*Zl33i*Zl36i/cb*sb**2+MW/CW**2*SW*Zh12
     &      *Zl33r*Zl36r/cb*sb**2-2D0*MW/CW**2*SW*Zh12*Zl63i
     &      *Zl66i/cb*sb**2-2D0*MW/CW**2*SW*Zh12*Zl63r*Zl66r
     &      /cb*sb**2+Ml/MW/SW*Zh22*Zl33i*Zl66i/cb*mur+Ml
     &      /MW/SW*Zh22*Zl33i*Zl66r/cb*mui-Ml/MW/SW*Zh22
     &      *Zl33r*Zl66i/cb*mui+Ml/MW/SW*Zh22*Zl33r*Zl66r
     &      /cb*mur+Ml/MW/SW*Zh22*Zl36i*Zl63i/cb*mur+Ml/MW
     &      /SW*Zh22*Zl36i*Zl63r/cb*mui-Ml/MW/SW*Zh22*Zl36r
     &      *Zl63i/cb*mui+Ml/MW/SW*Zh22*Zl36r*Zl63r/cb*mur
     &      +Ali/MW*Ml/SW*Zh32*Zl33i*Zl66i/cb*sb+Ali/MW*Ml
     &      /SW*Zh32*Zl33r*Zl66r/cb*sb+Ali/MW*Ml/SW*Zh32
     &      *Zl36i*Zl63i/cb*sb+Ali/MW*Ml/SW*Zh32*Zl36r*Zl63r
     &      /cb*sb+Alr/MW*Ml/SW*Zh32*Zl33i*Zl66r/cb*sb-Alr
     &      /MW*Ml/SW*Zh32*Zl33r*Zl66i/cb*sb+Alr/MW*Ml/SW
     &      *Zh32*Zl36i*Zl63r/cb*sb-Alr/MW*Ml/SW*Zh32*Zl36r
     &      *Zl63i/cb*sb)
      AAABR(98) = EE*(MW/SW*Zh12*Zl33i*Zl36r/cb-MW/SW*Zh12
     &      *Zl33r*Zl36i/cb-MW/SW*Zh22*Zl33i*Zl36r*sb+MW
     &      /SW*Zh22*Zl33r*Zl36i*sb-MW/CW**2*SW*Zh12*Zl33i
     &      *Zl36r/cb+MW/CW**2*SW*Zh12*Zl33r*Zl36i/cb+2D0
     &      *MW/CW**2*SW*Zh12*Zl63i*Zl66r/cb-2D0*MW/CW**2
     &      *SW*Zh12*Zl63r*Zl66i/cb+MW/CW**2*SW*Zh22*Zl33i
     &      *Zl36r*sb-MW/CW**2*SW*Zh22*Zl33r*Zl36i*sb-2D0
     &      *MW/CW**2*SW*Zh22*Zl63i*Zl66r*sb+2D0*MW/CW**2
     &      *SW*Zh22*Zl63r*Zl66i*sb-Ml/MW/SW*Zh32*Zl33i*Zl66i
     &      *mur-Ml/MW/SW*Zh32*Zl33i*Zl66r*mui+Ml/MW/SW*Zh32
     &      *Zl33r*Zl66i*mui-Ml/MW/SW*Zh32*Zl33r*Zl66r*mur
     &      +Ml/MW/SW*Zh32*Zl36i*Zl63i*mur+Ml/MW/SW*Zh32
     &      *Zl36i*Zl63r*mui-Ml/MW/SW*Zh32*Zl36r*Zl63i*mui
     &      +Ml/MW/SW*Zh32*Zl36r*Zl63r*mur-MW/SW*Zh12*Zl33i
     &      *Zl36r/cb*sb**2+MW/SW*Zh12*Zl33r*Zl36i/cb*sb
     &      **2-Ali/MW*Ml/SW*Zh12*Zl33i*Zl66i/cb-Ali/MW*Ml
     &      /SW*Zh12*Zl33r*Zl66r/cb+Ali/MW*Ml/SW*Zh12*Zl36i
     &      *Zl63i/cb+Ali/MW*Ml/SW*Zh12*Zl36r*Zl63r/cb-Alr
     &      /MW*Ml/SW*Zh12*Zl33i*Zl66r/cb+Alr/MW*Ml/SW*Zh12
     &      *Zl33r*Zl66i/cb+Alr/MW*Ml/SW*Zh12*Zl36i*Zl63r
     &      /cb-Alr/MW*Ml/SW*Zh12*Zl36r*Zl63i/cb+MW/CW**2
     &      *SW*Zh12*Zl33i*Zl36r/cb*sb**2-MW/CW**2*SW*Zh12
     &      *Zl33r*Zl36i/cb*sb**2-2D0*MW/CW**2*SW*Zh12*Zl63i
     &      *Zl66r/cb*sb**2+2D0*MW/CW**2*SW*Zh12*Zl63r*Zl66i
     &      /cb*sb**2-Ml/MW/SW*Zh22*Zl33i*Zl66i/cb*mui+Ml
     &      /MW/SW*Zh22*Zl33i*Zl66r/cb*mur-Ml/MW/SW*Zh22
     &      *Zl33r*Zl66i/cb*mur-Ml/MW/SW*Zh22*Zl33r*Zl66r
     &      /cb*mui+Ml/MW/SW*Zh22*Zl36i*Zl63i/cb*mui-Ml/MW
     &      /SW*Zh22*Zl36i*Zl63r/cb*mur+Ml/MW/SW*Zh22*Zl36r
     &      *Zl63i/cb*mur+Ml/MW/SW*Zh22*Zl36r*Zl63r/cb*mui
     &      +Ali/MW*Ml/SW*Zh32*Zl33i*Zl66r/cb*sb-Ali/MW*Ml
     &      /SW*Zh32*Zl33r*Zl66i/cb*sb-Ali/MW*Ml/SW*Zh32
     &      *Zl36i*Zl63r/cb*sb+Ali/MW*Ml/SW*Zh32*Zl36r*Zl63i
     &      /cb*sb-Alr/MW*Ml/SW*Zh32*Zl33i*Zl66i/cb*sb-Alr
     &      /MW*Ml/SW*Zh32*Zl33r*Zl66r/cb*sb+Alr/MW*Ml/SW
     &      *Zh32*Zl36i*Zl63i/cb*sb+Alr/MW*Ml/SW*Zh32*Zl36r
     &      *Zl63r/cb*sb)
      AAABR(99) = 2D0*Zl44/cb*ls1-2D0*Zl44/cb*ls1*sb**2-EE
     &      *MW/SW*Sqrt2*Zl14+EE/MW*Me**2/SW*Sqrt2*Zl14+2D0
     &      *EE*MW/SW*Sqrt2*Zl14*sb**2-EE/MW*Me/SW*Sqrt2
     &      *Zl44/cb*mur*sb
      AAABR(100) = EE/MW*Me/SW*Sqrt2*Zl44/cb*mui*sb
      AAABR(101) = 2D0*Zl55/cb*ls2-2D0*Zl55/cb*ls2*sb**2
     &      -EE*MW/SW*Sqrt2*Zl25+EE/MW*Mm**2/SW*Sqrt2*Zl25
     &      +2D0*EE*MW/SW*Sqrt2*Zl25*sb**2-EE/MW*Mm/SW*Sqrt2
     &      *Zl55/cb*mur*sb
      AAABR(102) = EE/MW*Mm/SW*Sqrt2*Zl55/cb*mui*sb
      AAABR(103) = EE/SW*Sqrt2*(MW*Zl36r-Ml**2/MW*Zl36r-2D0
     &      *MW*Zl36r*sb**2-Ali/MW*Ml*Zl66i-Alr/MW*Ml*Zl66r
     &      -Ml/MW*Zl66i/cb*mui*sb+Ml/MW*Zl66r/cb*mur*sb)
      AAABR(104) = EE/SW*Sqrt2*(MW*Zl36i-Ml**2/MW*Zl36i-2D0
     &      *MW*Zl36i*sb**2+Ali/MW*Ml*Zl66r-Alr/MW*Ml*Zl66i
     &      +Ml/MW*Zl66i/cb*mur*sb+Ml/MW*Zl66r/cb*mui*sb)
      AAABR(105) = EE/MW*Me/SW*Zl14*Zl44*mui*tb
      AAABR(106) = EE/MW*Ml/SW*(Ali*Zl36i*Zl66i+Ali*Zl36r
     &      *Zl66r+Alr*Zl36i*Zl66r-Alr*Zl36r*Zl66i+Zl36i
     &      *Zl66i/cb*mui*sb-Zl36i*Zl66r/cb*mur*sb+Zl36r
     &      *Zl66i/cb*mur*sb+Zl36r*Zl66r/cb*mui*sb)
      AAABR(107) = EE/MW*Mm/SW*Zl25*Zl55*mui*tb
      AAABR(108) = EE*(2D0*Ml**2/MW/SW*Zh12/cb-MW/SW*Zh12
     &      *Zl36i**2/cb-MW/SW*Zh12*Zl36r**2/cb+MW/SW*Zh22
     &      *Zl36i**2*sb+MW/SW*Zh22*Zl36r**2*sb+MW/CW**2
     &      *SW*Zh12*Zl36i**2/cb+MW/CW**2*SW*Zh12*Zl36r**2
     &      /cb-2D0*MW/CW**2*SW*Zh12*Zl66i**2/cb-2D0*MW/CW
     &      **2*SW*Zh12*Zl66r**2/cb-MW/CW**2*SW*Zh22*Zl36i
     &      **2*sb-MW/CW**2*SW*Zh22*Zl36r**2*sb+2D0*MW/CW
     &      **2*SW*Zh22*Zl66i**2*sb+2D0*MW/CW**2*SW*Zh22
     &      *Zl66r**2*sb+MW/SW*Zh12*Zl36i**2/cb*sb**2+MW
     &      /SW*Zh12*Zl36r**2/cb*sb**2-MW/CW**2*SW*Zh12*Zl36i
     &      **2/cb*sb**2-MW/CW**2*SW*Zh12*Zl36r**2/cb*sb
     &      **2+2D0*MW/CW**2*SW*Zh12*Zl66i**2/cb*sb**2+2D0
     &      *MW/CW**2*SW*Zh12*Zl66r**2/cb*sb**2+2D0*Ml/MW
     &      /SW*Zh32*Zl36i*Zl66i*mui-2D0*Ml/MW/SW*Zh32*Zl36i
     &      *Zl66r*mur+2D0*Ml/MW/SW*Zh32*Zl36r*Zl66i*mur
     &      +2D0*Ml/MW/SW*Zh32*Zl36r*Zl66r*mui-2D0*Ali/MW
     &      *Ml/SW*Zh12*Zl36i*Zl66r/cb+2D0*Ali/MW*Ml/SW*Zh12
     &      *Zl36r*Zl66i/cb+2D0*Alr/MW*Ml/SW*Zh12*Zl36i*Zl66i
     &      /cb+2D0*Alr/MW*Ml/SW*Zh12*Zl36r*Zl66r/cb-2D0
     &      *Ml/MW/SW*Zh22*Zl36i*Zl66i/cb*mur-2D0*Ml/MW/SW
     &      *Zh22*Zl36i*Zl66r/cb*mui+2D0*Ml/MW/SW*Zh22*Zl36r
     &      *Zl66i/cb*mui-2D0*Ml/MW/SW*Zh22*Zl36r*Zl66r/cb
     &      *mur-2D0*Ali/MW*Ml/SW*Zh32*Zl36i*Zl66i/cb*sb
     &      -2D0*Ali/MW*Ml/SW*Zh32*Zl36r*Zl66r/cb*sb-2D0
     &      *Alr/MW*Ml/SW*Zh32*Zl36i*Zl66r/cb*sb+2D0*Alr
     &      /MW*Ml/SW*Zh32*Zl36r*Zl66i/cb*sb)
      AAABR(109) = 2D0*Sqrt2*Zh11*Zl25*Zl55*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl25**2/cb+EE
     &      *MW/SW*Zh21*Zl25**2*sb+EE/CW**2*MW*SW*Zh11*Zl25
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl55**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl25**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl55**2*sb+EE*MW/SW*Zh11*Zl25**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl25**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl55**2/cb*sb**2+2D0*EE
     &      /MW*Mm/SW*Zh31*Zl25*Zl55*mui-2D0*EE/MW*Mm/SW
     &      *Zh21*Zl25*Zl55/cb*mur
      AAABR(110) = EE*(2D0*Ml**2/MW/SW*Zh13/cb-MW/SW*Zh13
     &      *Zl36i**2/cb-MW/SW*Zh13*Zl36r**2/cb+MW/SW*Zh23
     &      *Zl36i**2*sb+MW/SW*Zh23*Zl36r**2*sb+MW/CW**2
     &      *SW*Zh13*Zl36i**2/cb+MW/CW**2*SW*Zh13*Zl36r**2
     &      /cb-2D0*MW/CW**2*SW*Zh13*Zl66i**2/cb-2D0*MW/CW
     &      **2*SW*Zh13*Zl66r**2/cb-MW/CW**2*SW*Zh23*Zl36i
     &      **2*sb-MW/CW**2*SW*Zh23*Zl36r**2*sb+2D0*MW/CW
     &      **2*SW*Zh23*Zl66i**2*sb+2D0*MW/CW**2*SW*Zh23
     &      *Zl66r**2*sb+MW/SW*Zh13*Zl36i**2/cb*sb**2+MW
     &      /SW*Zh13*Zl36r**2/cb*sb**2-MW/CW**2*SW*Zh13*Zl36i
     &      **2/cb*sb**2-MW/CW**2*SW*Zh13*Zl36r**2/cb*sb
     &      **2+2D0*MW/CW**2*SW*Zh13*Zl66i**2/cb*sb**2+2D0
     &      *MW/CW**2*SW*Zh13*Zl66r**2/cb*sb**2+2D0*Ml/MW
     &      /SW*Zh33*Zl36i*Zl66i*mui-2D0*Ml/MW/SW*Zh33*Zl36i
     &      *Zl66r*mur+2D0*Ml/MW/SW*Zh33*Zl36r*Zl66i*mur
     &      +2D0*Ml/MW/SW*Zh33*Zl36r*Zl66r*mui-2D0*Ali/MW
     &      *Ml/SW*Zh13*Zl36i*Zl66r/cb+2D0*Ali/MW*Ml/SW*Zh13
     &      *Zl36r*Zl66i/cb+2D0*Alr/MW*Ml/SW*Zh13*Zl36i*Zl66i
     &      /cb+2D0*Alr/MW*Ml/SW*Zh13*Zl36r*Zl66r/cb-2D0
     &      *Ml/MW/SW*Zh23*Zl36i*Zl66i/cb*mur-2D0*Ml/MW/SW
     &      *Zh23*Zl36i*Zl66r/cb*mui+2D0*Ml/MW/SW*Zh23*Zl36r
     &      *Zl66i/cb*mui-2D0*Ml/MW/SW*Zh23*Zl36r*Zl66r/cb
     &      *mur-2D0*Ali/MW*Ml/SW*Zh33*Zl36i*Zl66i/cb*sb
     &      -2D0*Ali/MW*Ml/SW*Zh33*Zl36r*Zl66r/cb*sb-2D0
     &      *Alr/MW*Ml/SW*Zh33*Zl36i*Zl66r/cb*sb+2D0*Alr
     &      /MW*Ml/SW*Zh33*Zl36r*Zl66i/cb*sb)
      AAABR(111) = 2D0*Sqrt2*Zh13*Zl14*Zl44*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl14**2/cb+EE
     &      *MW/SW*Zh23*Zl14**2*sb+EE/CW**2*MW*SW*Zh13*Zl14
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl44**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl14**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl44**2*sb+EE*MW/SW*Zh13*Zl14**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl14**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl44**2/cb*sb**2+2D0*EE
     &      /MW*Me/SW*Zh33*Zl14*Zl44*mui-2D0*EE/MW*Me/SW
     &      *Zh23*Zl14*Zl44/cb*mur
      AAABR(112) = 2D0*Sqrt2*Zh13*Zl25*Zl55*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl25**2/cb+EE
     &      *MW/SW*Zh23*Zl25**2*sb+EE/CW**2*MW*SW*Zh13*Zl25
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl55**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl25**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl55**2*sb+EE*MW/SW*Zh13*Zl25**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl25**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl55**2/cb*sb**2+2D0*EE
     &      /MW*Mm/SW*Zh33*Zl25*Zl55*mui-2D0*EE/MW*Mm/SW
     &      *Zh23*Zl25*Zl55/cb*mur
      AAABR(113) = 2D0*Sqrt2*Zh12*Zl14*Zl44*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl14**2/cb+EE
     &      *MW/SW*Zh22*Zl14**2*sb+EE/CW**2*MW*SW*Zh12*Zl14
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl44**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl14**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl44**2*sb+EE*MW/SW*Zh12*Zl14**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl14**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl44**2/cb*sb**2+2D0*EE
     &      /MW*Me/SW*Zh32*Zl14*Zl44*mui-2D0*EE/MW*Me/SW
     &      *Zh22*Zl14*Zl44/cb*mur
      AAABR(114) = EE*(2D0*Ml**2/MW/SW*Zh11/cb-MW/SW*Zh11
     &      *Zl36i**2/cb-MW/SW*Zh11*Zl36r**2/cb+MW/SW*Zh21
     &      *Zl36i**2*sb+MW/SW*Zh21*Zl36r**2*sb+MW/CW**2
     &      *SW*Zh11*Zl36i**2/cb+MW/CW**2*SW*Zh11*Zl36r**2
     &      /cb-2D0*MW/CW**2*SW*Zh11*Zl66i**2/cb-2D0*MW/CW
     &      **2*SW*Zh11*Zl66r**2/cb-MW/CW**2*SW*Zh21*Zl36i
     &      **2*sb-MW/CW**2*SW*Zh21*Zl36r**2*sb+2D0*MW/CW
     &      **2*SW*Zh21*Zl66i**2*sb+2D0*MW/CW**2*SW*Zh21
     &      *Zl66r**2*sb+MW/SW*Zh11*Zl36i**2/cb*sb**2+MW
     &      /SW*Zh11*Zl36r**2/cb*sb**2-MW/CW**2*SW*Zh11*Zl36i
     &      **2/cb*sb**2-MW/CW**2*SW*Zh11*Zl36r**2/cb*sb
     &      **2+2D0*MW/CW**2*SW*Zh11*Zl66i**2/cb*sb**2+2D0
     &      *MW/CW**2*SW*Zh11*Zl66r**2/cb*sb**2+2D0*Ml/MW
     &      /SW*Zh31*Zl36i*Zl66i*mui-2D0*Ml/MW/SW*Zh31*Zl36i
     &      *Zl66r*mur+2D0*Ml/MW/SW*Zh31*Zl36r*Zl66i*mur
     &      +2D0*Ml/MW/SW*Zh31*Zl36r*Zl66r*mui-2D0*Ali/MW
     &      *Ml/SW*Zh11*Zl36i*Zl66r/cb+2D0*Ali/MW*Ml/SW*Zh11
     &      *Zl36r*Zl66i/cb+2D0*Alr/MW*Ml/SW*Zh11*Zl36i*Zl66i
     &      /cb+2D0*Alr/MW*Ml/SW*Zh11*Zl36r*Zl66r/cb-2D0
     &      *Ml/MW/SW*Zh21*Zl36i*Zl66i/cb*mur-2D0*Ml/MW/SW
     &      *Zh21*Zl36i*Zl66r/cb*mui+2D0*Ml/MW/SW*Zh21*Zl36r
     &      *Zl66i/cb*mui-2D0*Ml/MW/SW*Zh21*Zl36r*Zl66r/cb
     &      *mur-2D0*Ali/MW*Ml/SW*Zh31*Zl36i*Zl66i/cb*sb
     &      -2D0*Ali/MW*Ml/SW*Zh31*Zl36r*Zl66r/cb*sb-2D0
     &      *Alr/MW*Ml/SW*Zh31*Zl36i*Zl66r/cb*sb+2D0*Alr
     &      /MW*Ml/SW*Zh31*Zl36r*Zl66i/cb*sb)
      AAABR(115) = 2D0*Sqrt2*Zh12*Zl25*Zl55*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl25**2/cb+EE
     &      *MW/SW*Zh22*Zl25**2*sb+EE/CW**2*MW*SW*Zh12*Zl25
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl55**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl25**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl55**2*sb+EE*MW/SW*Zh12*Zl25**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl25**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl55**2/cb*sb**2+2D0*EE
     &      /MW*Mm/SW*Zh32*Zl25*Zl55*mui-2D0*EE/MW*Mm/SW
     &      *Zh22*Zl25*Zl55/cb*mur
      AAABR(116) = 2D0*Sqrt2*Zh11*Zl14*Zl44*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl14**2/cb+EE
     &      *MW/SW*Zh21*Zl14**2*sb+EE/CW**2*MW*SW*Zh11*Zl14
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl44**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl14**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl44**2*sb+EE*MW/SW*Zh11*Zl14**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl14**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl44**2/cb*sb**2+2D0*EE
     &      /MW*Me/SW*Zh31*Zl14*Zl44*mui-2D0*EE/MW*Me/SW
     &      *Zh21*Zl14*Zl44/cb*mur
      AAABR(117) = EE*MW*(Zh13/SW*cb-Zh23/SW*sb+SW/CW**2
     &      *Zh13*cb-SW/CW**2*Zh23*sb)
      AAABR(118) = EE*MW*(Zh11/SW*cb-Zh21/SW*sb+SW/CW**2
     &      *Zh11*cb-SW/CW**2*Zh21*sb)
      AAABR(119) = EE*MW*(Zh12/SW*cb-Zh22/SW*sb+SW/CW**2
     &      *Zh12*cb-SW/CW**2*Zh22*sb)
      AAABR(120) = 2D0*Zd22*Zu52*sb*us2+2D0*Zd52*Zu22/cb
     &      *ds2-2D0*Zd52*Zu22/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd22*Zu22-EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu22+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd22*Zu22+2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu22*sb**2+EE/MW*Mc/SW*Sqrt2*Zd22*Zu52
     &      /cb*mur/sb-EE/MW*Mc/SW*Sqrt2*Zd22*Zu52/cb*mur
     &      *sb-EE/MW*Ms/SW*Sqrt2*Zd52*Zu22/cb*mur*sb
      AAABR(121) = EE/MW/SW*Sqrt2/cb*mui*(Mc*Zd22*Zu52/sb
     &      -Mc*Zd22*Zu52*sb+Ms*Zd52*Zu22*sb)
      AAABR(122) = 2D0*Zd11*Zu41*sb*us1+2D0*Zd41*Zu11*cb
     &      *ds1-EE*MW/SW*Sqrt2*Zd11*Zu11+2D0*EE*MW/SW*Sqrt2
     &      *Zd11*Zu11*sb**2
      AAABR(123) = EE/SW*Sqrt2*(MW*Zd33i*Zu33i+MW*Zd33r*Zu33r
     &      -Mb**2/MW*Zd33i*Zu33i-Mb**2/MW*Zd33r*Zu33r+MtopS
     &      **2/MW*Zd33i*Zu33i+MtopS**2/MW*Zd33r*Zu33r-2D0
     &      *MW*Zd33i*Zu33i*sb**2-2D0*MW*Zd33r*Zu33r*sb**2
     &      -Abi/MW*Mb*Zd63i*Zu33r+Abi/MW*Mb*Zd63r*Zu33i
     &      -Abr/MW*Mb*Zd63i*Zu33i-Abr/MW*Mb*Zd63r*Zu33r
     &      -Ati/MW*MtopS*Zd33i*Zu63r+Ati/MW*MtopS*Zd33r
     &      *Zu63i+Atr/MW*MtopS*Zd33i*Zu63i+Atr/MW*MtopS
     &      *Zd33r*Zu63r+Mb/MW*Zd63i*Zu33i/cb*mur*sb-Mb/MW
     &      *Zd63i*Zu33r/cb*mui*sb+Mb/MW*Zd63r*Zu33i/cb*mui
     &      *sb+Mb/MW*Zd63r*Zu33r/cb*mur*sb-MtopS/MW*Zd33i
     &      *Zu63i/cb*mur/sb+MtopS/MW*Zd33i*Zu63i/cb*mur
     &      *sb-MtopS/MW*Zd33i*Zu63r/cb*mui/sb+MtopS/MW*Zd33i
     &      *Zu63r/cb*mui*sb+MtopS/MW*Zd33r*Zu63i/cb*mui
     &      /sb-MtopS/MW*Zd33r*Zu63i/cb*mui*sb-MtopS/MW*Zd33r
     &      *Zu63r/cb*mur/sb+MtopS/MW*Zd33r*Zu63r/cb*mur*sb)
      AAABR(124) = EE/SW*Sqrt2*(MW*Zd33i*Zu33r-MW*Zd33r*Zu33i
     &      -Mb**2/MW*Zd33i*Zu33r+Mb**2/MW*Zd33r*Zu33i+MtopS
     &      **2/MW*Zd33i*Zu33r-MtopS**2/MW*Zd33r*Zu33i-2D0
     &      *MW*Zd33i*Zu33r*sb**2+2D0*MW*Zd33r*Zu33i*sb**2
     &      +Abi/MW*Mb*Zd63i*Zu33i+Abi/MW*Mb*Zd63r*Zu33r
     &      -Abr/MW*Mb*Zd63i*Zu33r+Abr/MW*Mb*Zd63r*Zu33i
     &      +Ati/MW*MtopS*Zd33i*Zu63i+Ati/MW*MtopS*Zd33r
     &      *Zu63r+Atr/MW*MtopS*Zd33i*Zu63r-Atr/MW*MtopS
     &      *Zd33r*Zu63i+Mb/MW*Zd63i*Zu33i/cb*mui*sb+Mb/MW
     &      *Zd63i*Zu33r/cb*mur*sb-Mb/MW*Zd63r*Zu33i/cb*mur
     &      *sb+Mb/MW*Zd63r*Zu33r/cb*mui*sb+MtopS/MW*Zd33i
     &      *Zu63i/cb*mui/sb-MtopS/MW*Zd33i*Zu63i/cb*mui
     &      *sb-MtopS/MW*Zd33i*Zu63r/cb*mur/sb+MtopS/MW*Zd33i
     &      *Zu63r/cb*mur*sb+MtopS/MW*Zd33r*Zu63i/cb*mur
     &      /sb-MtopS/MW*Zd33r*Zu63i/cb*mur*sb+MtopS/MW*Zd33r
     &      *Zu63r/cb*mui/sb-MtopS/MW*Zd33r*Zu63r/cb*mui*sb)
      AAABR(125) = EE/SW*Sqrt2*(MW*Zd33i*Zu36i+MW*Zd33r*Zu36r
     &      -Mb**2/MW*Zd33i*Zu36i-Mb**2/MW*Zd33r*Zu36r+MtopS
     &      **2/MW*Zd33i*Zu36i+MtopS**2/MW*Zd33r*Zu36r-2D0
     &      *MW*Zd33i*Zu36i*sb**2-2D0*MW*Zd33r*Zu36r*sb**2
     &      -Abi/MW*Mb*Zd63i*Zu36r+Abi/MW*Mb*Zd63r*Zu36i
     &      -Abr/MW*Mb*Zd63i*Zu36i-Abr/MW*Mb*Zd63r*Zu36r
     &      -Ati/MW*MtopS*Zd33i*Zu66r+Ati/MW*MtopS*Zd33r
     &      *Zu66i+Atr/MW*MtopS*Zd33i*Zu66i+Atr/MW*MtopS
     &      *Zd33r*Zu66r+Mb/MW*Zd63i*Zu36i/cb*mur*sb-Mb/MW
     &      *Zd63i*Zu36r/cb*mui*sb+Mb/MW*Zd63r*Zu36i/cb*mui
     &      *sb+Mb/MW*Zd63r*Zu36r/cb*mur*sb-MtopS/MW*Zd33i
     &      *Zu66i/cb*mur/sb+MtopS/MW*Zd33i*Zu66i/cb*mur
     &      *sb-MtopS/MW*Zd33i*Zu66r/cb*mui/sb+MtopS/MW*Zd33i
     &      *Zu66r/cb*mui*sb+MtopS/MW*Zd33r*Zu66i/cb*mui
     &      /sb-MtopS/MW*Zd33r*Zu66i/cb*mui*sb-MtopS/MW*Zd33r
     &      *Zu66r/cb*mur/sb+MtopS/MW*Zd33r*Zu66r/cb*mur*sb)
      AAABR(126) = EE/SW*Sqrt2*(MW*Zd33i*Zu36r-MW*Zd33r*Zu36i
     &      -Mb**2/MW*Zd33i*Zu36r+Mb**2/MW*Zd33r*Zu36i+MtopS
     &      **2/MW*Zd33i*Zu36r-MtopS**2/MW*Zd33r*Zu36i-2D0
     &      *MW*Zd33i*Zu36r*sb**2+2D0*MW*Zd33r*Zu36i*sb**2
     &      +Abi/MW*Mb*Zd63i*Zu36i+Abi/MW*Mb*Zd63r*Zu36r
     &      -Abr/MW*Mb*Zd63i*Zu36r+Abr/MW*Mb*Zd63r*Zu36i
     &      +Ati/MW*MtopS*Zd33i*Zu66i+Ati/MW*MtopS*Zd33r
     &      *Zu66r+Atr/MW*MtopS*Zd33i*Zu66r-Atr/MW*MtopS
     &      *Zd33r*Zu66i+Mb/MW*Zd63i*Zu36i/cb*mui*sb+Mb/MW
     &      *Zd63i*Zu36r/cb*mur*sb-Mb/MW*Zd63r*Zu36i/cb*mur
     &      *sb+Mb/MW*Zd63r*Zu36r/cb*mui*sb+MtopS/MW*Zd33i
     &      *Zu66i/cb*mui/sb-MtopS/MW*Zd33i*Zu66i/cb*mui
     &      *sb-MtopS/MW*Zd33i*Zu66r/cb*mur/sb+MtopS/MW*Zd33i
     &      *Zu66r/cb*mur*sb+MtopS/MW*Zd33r*Zu66i/cb*mur
     &      /sb-MtopS/MW*Zd33r*Zu66i/cb*mur*sb+MtopS/MW*Zd33r
     &      *Zu66r/cb*mui/sb-MtopS/MW*Zd33r*Zu66r/cb*mui*sb)
      AAABR(127) = 2D0*Zd11*Zu44*sb*us1+2D0*Zd41*Zu14*cb
     &      *ds1-EE*MW/SW*Sqrt2*Zd11*Zu14+2D0*EE*MW/SW*Sqrt2
     &      *Zd11*Zu14*sb**2
      AAABR(128) = 2D0*Zd22*Zu55*sb*us2+2D0*Zd52*Zu25/cb
     &      *ds2-2D0*Zd52*Zu25/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd22*Zu25-EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu25+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd22*Zu25+2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu25*sb**2+EE/MW*Mc/SW*Sqrt2*Zd22*Zu55
     &      /cb*mur/sb-EE/MW*Mc/SW*Sqrt2*Zd22*Zu55/cb*mur
     &      *sb-EE/MW*Ms/SW*Sqrt2*Zd52*Zu25/cb*mur*sb
      AAABR(129) = EE/MW/SW*Sqrt2/cb*mui*(Mc*Zd22*Zu55/sb
     &      -Mc*Zd22*Zu55*sb+Ms*Zd52*Zu25*sb)
      AAABR(130) = EE/MW*Ms/SW*Zd22*Zd52*mui*tb
      AAABR(131) = EE/MW*Mb/SW*(Abi*Zd33i*Zd63i+Abi*Zd33r
     &      *Zd63r+Abr*Zd33i*Zd63r-Abr*Zd33r*Zd63i+Zd33i
     &      *Zd63i/cb*mui*sb-Zd33i*Zd63r/cb*mur*sb+Zd33r
     &      *Zd63i/cb*mur*sb+Zd33r*Zd63r/cb*mui*sb)
      AAABR(132) = 1D0/cb*(Sqrt2*Zd22*Zd55*ds2-Sqrt2*Zd25
     &      *Zd52*ds2-Sqrt2*Zd22*Zd55*ds2*sb**2+Sqrt2*Zd25
     &      *Zd52*ds2*sb**2-EE/MW*Ms/SW*Zd22*Zd55*mur*sb
     &      +EE/MW*Ms/SW*Zd25*Zd52*mur*sb)
      AAABR(133) = EE/MW*Ms/SW/cb*mui*sb*(Zd22*Zd55+Zd25
     &      *Zd52)
      AAABR(134) = EE/MW*Mb/SW*(Abi*Zd33i*Zd66r-Abi*Zd33r
     &      *Zd66i-Abi*Zd36i*Zd63r+Abi*Zd36r*Zd63i-Abr*Zd33i
     &      *Zd66i-Abr*Zd33r*Zd66r+Abr*Zd36i*Zd63i+Abr*Zd36r
     &      *Zd63r+Zd33i*Zd66i/cb*mur*sb+Zd33i*Zd66r/cb*mui
     &      *sb-Zd33r*Zd66i/cb*mui*sb+Zd33r*Zd66r/cb*mur
     &      *sb-Zd36i*Zd63i/cb*mur*sb-Zd36i*Zd63r/cb*mui
     &      *sb+Zd36r*Zd63i/cb*mui*sb-Zd36r*Zd63r/cb*mur*sb)
      AAABR(135) = EE/MW*Mb/SW*(Abi*Zd33i*Zd66i+Abi*Zd33r
     &      *Zd66r+Abi*Zd36i*Zd63i+Abi*Zd36r*Zd63r+Abr*Zd33i
     &      *Zd66r-Abr*Zd33r*Zd66i+Abr*Zd36i*Zd63r-Abr*Zd36r
     &      *Zd63i+Zd33i*Zd66i/cb*mui*sb-Zd33i*Zd66r/cb*mur
     &      *sb+Zd33r*Zd66i/cb*mur*sb+Zd33r*Zd66r/cb*mui
     &      *sb+Zd36i*Zd63i/cb*mui*sb-Zd36i*Zd63r/cb*mur
     &      *sb+Zd36r*Zd63i/cb*mur*sb+Zd36r*Zd63r/cb*mui*sb)
      AAABR(136) = Sqrt2*cb*ds1*(Zd11*Zd44-Zd14*Zd41)
      AAABR(137) = 6D0*Sqrt2*Zd22*Zd52*Zh11*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh11/cb-3D0*EE*MW/SW*Zd22**2*Zh11/cb
     &      +3D0*EE*MW/SW*Zd22**2*Zh21*sb-EE/CW**2*MW*SW
     &      *Zd22**2*Zh11/cb+EE/CW**2*MW*SW*Zd22**2*Zh21
     &      *sb-2D0*EE/CW**2*MW*SW*Zd52**2*Zh11/cb+2D0*EE
     &      /CW**2*MW*SW*Zd52**2*Zh21*sb+3D0*EE*MW/SW*Zd22
     &      **2*Zh11/cb*sb**2+EE/CW**2*MW*SW*Zd22**2*Zh11
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd52**2*Zh11/cb
     &      *sb**2+6D0*EE/MW*Ms/SW*Zd22*Zd52*Zh31*mui-6D0
     &      *EE/MW*Ms/SW*Zd22*Zd52*Zh21/cb*mur
      AAABR(138) = 6D0*Sqrt2*Zd11*Zd41*Zh12*ds1-3D0*EE*MW
     &      /SW*Zd11**2*Zh12*cb+3D0*EE*MW/SW*Zd11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zd11**2*Zh12*cb+EE/CW**2*MW
     &      *SW*Zd11**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zd41**2
     &      *Zh12*cb+2D0*EE/CW**2*MW*SW*Zd41**2*Zh22*sb
      AAABR(139) = 6D0*Sqrt2*Zd22*Zd52*Zh13*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh13/cb-3D0*EE*MW/SW*Zd22**2*Zh13/cb
     &      +3D0*EE*MW/SW*Zd22**2*Zh23*sb-EE/CW**2*MW*SW
     &      *Zd22**2*Zh13/cb+EE/CW**2*MW*SW*Zd22**2*Zh23
     &      *sb-2D0*EE/CW**2*MW*SW*Zd52**2*Zh13/cb+2D0*EE
     &      /CW**2*MW*SW*Zd52**2*Zh23*sb+3D0*EE*MW/SW*Zd22
     &      **2*Zh13/cb*sb**2+EE/CW**2*MW*SW*Zd22**2*Zh13
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd52**2*Zh13/cb
     &      *sb**2+6D0*EE/MW*Ms/SW*Zd22*Zd52*Zh33*mui-6D0
     &      *EE/MW*Ms/SW*Zd22*Zd52*Zh23/cb*mur
      AAABR(140) = EE*(6D0*Mb**2/MW/SW*Zh11/cb-3D0*MW/SW
     &      *Zd33i**2*Zh11/cb+3D0*MW/SW*Zd33i**2*Zh21*sb
     &      -3D0*MW/SW*Zd33r**2*Zh11/cb+3D0*MW/SW*Zd33r**2
     &      *Zh21*sb-MW/CW**2*SW*Zd33i**2*Zh11/cb+MW/CW**2
     &      *SW*Zd33i**2*Zh21*sb-MW/CW**2*SW*Zd33r**2*Zh11
     &      /cb+MW/CW**2*SW*Zd33r**2*Zh21*sb-2D0*MW/CW**2
     &      *SW*Zd63i**2*Zh11/cb+2D0*MW/CW**2*SW*Zd63i**2
     &      *Zh21*sb-2D0*MW/CW**2*SW*Zd63r**2*Zh11/cb+2D0
     &      *MW/CW**2*SW*Zd63r**2*Zh21*sb+3D0*MW/SW*Zd33i
     &      **2*Zh11/cb*sb**2+3D0*MW/SW*Zd33r**2*Zh11/cb
     &      *sb**2+MW/CW**2*SW*Zd33i**2*Zh11/cb*sb**2+MW
     &      /CW**2*SW*Zd33r**2*Zh11/cb*sb**2+2D0*MW/CW**2
     &      *SW*Zd63i**2*Zh11/cb*sb**2+2D0*MW/CW**2*SW*Zd63r
     &      **2*Zh11/cb*sb**2+6D0*Mb/MW/SW*Zd33i*Zd63i*Zh31
     &      *mui-6D0*Mb/MW/SW*Zd33i*Zd63r*Zh31*mur+6D0*Mb
     &      /MW/SW*Zd33r*Zd63i*Zh31*mur+6D0*Mb/MW/SW*Zd33r
     &      *Zd63r*Zh31*mui-6D0*Abi/MW*Mb/SW*Zd33i*Zd63r
     &      *Zh11/cb+6D0*Abi/MW*Mb/SW*Zd33r*Zd63i*Zh11/cb
     &      +6D0*Abr/MW*Mb/SW*Zd33i*Zd63i*Zh11/cb+6D0*Abr
     &      /MW*Mb/SW*Zd33r*Zd63r*Zh11/cb-6D0*Mb/MW/SW*Zd33i
     &      *Zd63i*Zh21/cb*mur-6D0*Mb/MW/SW*Zd33i*Zd63r*Zh21
     &      /cb*mui+6D0*Mb/MW/SW*Zd33r*Zd63i*Zh21/cb*mui
     &      -6D0*Mb/MW/SW*Zd33r*Zd63r*Zh21/cb*mur-6D0*Abi
     &      /MW*Mb/SW*Zd33i*Zd63i*Zh31/cb*sb-6D0*Abi/MW*Mb
     &      /SW*Zd33r*Zd63r*Zh31/cb*sb-6D0*Abr/MW*Mb/SW*Zd33i
     &      *Zd63r*Zh31/cb*sb+6D0*Abr/MW*Mb/SW*Zd33r*Zd63i
     &      *Zh31/cb*sb)
      AAABR(141) = 6D0*Sqrt2*Zd11*Zd41*Zh13*ds1-3D0*EE*MW
     &      /SW*Zd11**2*Zh13*cb+3D0*EE*MW/SW*Zd11**2*Zh23
     &      *sb-EE/CW**2*MW*SW*Zd11**2*Zh13*cb+EE/CW**2*MW
     &      *SW*Zd11**2*Zh23*sb-2D0*EE/CW**2*MW*SW*Zd41**2
     &      *Zh13*cb+2D0*EE/CW**2*MW*SW*Zd41**2*Zh23*sb
      AAABR(142) = EE*(6D0*Mb**2/MW/SW*Zh13/cb-3D0*MW/SW
     &      *Zd33i**2*Zh13/cb+3D0*MW/SW*Zd33i**2*Zh23*sb
     &      -3D0*MW/SW*Zd33r**2*Zh13/cb+3D0*MW/SW*Zd33r**2
     &      *Zh23*sb-MW/CW**2*SW*Zd33i**2*Zh13/cb+MW/CW**2
     &      *SW*Zd33i**2*Zh23*sb-MW/CW**2*SW*Zd33r**2*Zh13
     &      /cb+MW/CW**2*SW*Zd33r**2*Zh23*sb-2D0*MW/CW**2
     &      *SW*Zd63i**2*Zh13/cb+2D0*MW/CW**2*SW*Zd63i**2
     &      *Zh23*sb-2D0*MW/CW**2*SW*Zd63r**2*Zh13/cb+2D0
     &      *MW/CW**2*SW*Zd63r**2*Zh23*sb+3D0*MW/SW*Zd33i
     &      **2*Zh13/cb*sb**2+3D0*MW/SW*Zd33r**2*Zh13/cb
     &      *sb**2+MW/CW**2*SW*Zd33i**2*Zh13/cb*sb**2+MW
     &      /CW**2*SW*Zd33r**2*Zh13/cb*sb**2+2D0*MW/CW**2
     &      *SW*Zd63i**2*Zh13/cb*sb**2+2D0*MW/CW**2*SW*Zd63r
     &      **2*Zh13/cb*sb**2+6D0*Mb/MW/SW*Zd33i*Zd63i*Zh33
     &      *mui-6D0*Mb/MW/SW*Zd33i*Zd63r*Zh33*mur+6D0*Mb
     &      /MW/SW*Zd33r*Zd63i*Zh33*mur+6D0*Mb/MW/SW*Zd33r
     &      *Zd63r*Zh33*mui-6D0*Abi/MW*Mb/SW*Zd33i*Zd63r
     &      *Zh13/cb+6D0*Abi/MW*Mb/SW*Zd33r*Zd63i*Zh13/cb
     &      +6D0*Abr/MW*Mb/SW*Zd33i*Zd63i*Zh13/cb+6D0*Abr
     &      /MW*Mb/SW*Zd33r*Zd63r*Zh13/cb-6D0*Mb/MW/SW*Zd33i
     &      *Zd63i*Zh23/cb*mur-6D0*Mb/MW/SW*Zd33i*Zd63r*Zh23
     &      /cb*mui+6D0*Mb/MW/SW*Zd33r*Zd63i*Zh23/cb*mui
     &      -6D0*Mb/MW/SW*Zd33r*Zd63r*Zh23/cb*mur-6D0*Abi
     &      /MW*Mb/SW*Zd33i*Zd63i*Zh33/cb*sb-6D0*Abi/MW*Mb
     &      /SW*Zd33r*Zd63r*Zh33/cb*sb-6D0*Abr/MW*Mb/SW*Zd33i
     &      *Zd63r*Zh33/cb*sb+6D0*Abr/MW*Mb/SW*Zd33r*Zd63i
     &      *Zh33/cb*sb)
      AAABR(143) = EE*(6D0*Mb**2/MW/SW*Zh12/cb-3D0*MW/SW
     &      *Zd33i**2*Zh12/cb+3D0*MW/SW*Zd33i**2*Zh22*sb
     &      -3D0*MW/SW*Zd33r**2*Zh12/cb+3D0*MW/SW*Zd33r**2
     &      *Zh22*sb-MW/CW**2*SW*Zd33i**2*Zh12/cb+MW/CW**2
     &      *SW*Zd33i**2*Zh22*sb-MW/CW**2*SW*Zd33r**2*Zh12
     &      /cb+MW/CW**2*SW*Zd33r**2*Zh22*sb-2D0*MW/CW**2
     &      *SW*Zd63i**2*Zh12/cb+2D0*MW/CW**2*SW*Zd63i**2
     &      *Zh22*sb-2D0*MW/CW**2*SW*Zd63r**2*Zh12/cb+2D0
     &      *MW/CW**2*SW*Zd63r**2*Zh22*sb+3D0*MW/SW*Zd33i
     &      **2*Zh12/cb*sb**2+3D0*MW/SW*Zd33r**2*Zh12/cb
     &      *sb**2+MW/CW**2*SW*Zd33i**2*Zh12/cb*sb**2+MW
     &      /CW**2*SW*Zd33r**2*Zh12/cb*sb**2+2D0*MW/CW**2
     &      *SW*Zd63i**2*Zh12/cb*sb**2+2D0*MW/CW**2*SW*Zd63r
     &      **2*Zh12/cb*sb**2+6D0*Mb/MW/SW*Zd33i*Zd63i*Zh32
     &      *mui-6D0*Mb/MW/SW*Zd33i*Zd63r*Zh32*mur+6D0*Mb
     &      /MW/SW*Zd33r*Zd63i*Zh32*mur+6D0*Mb/MW/SW*Zd33r
     &      *Zd63r*Zh32*mui-6D0*Abi/MW*Mb/SW*Zd33i*Zd63r
     &      *Zh12/cb+6D0*Abi/MW*Mb/SW*Zd33r*Zd63i*Zh12/cb
     &      +6D0*Abr/MW*Mb/SW*Zd33i*Zd63i*Zh12/cb+6D0*Abr
     &      /MW*Mb/SW*Zd33r*Zd63r*Zh12/cb-6D0*Mb/MW/SW*Zd33i
     &      *Zd63i*Zh22/cb*mur-6D0*Mb/MW/SW*Zd33i*Zd63r*Zh22
     &      /cb*mui+6D0*Mb/MW/SW*Zd33r*Zd63i*Zh22/cb*mui
     &      -6D0*Mb/MW/SW*Zd33r*Zd63r*Zh22/cb*mur-6D0*Abi
     &      /MW*Mb/SW*Zd33i*Zd63i*Zh32/cb*sb-6D0*Abi/MW*Mb
     &      /SW*Zd33r*Zd63r*Zh32/cb*sb-6D0*Abr/MW*Mb/SW*Zd33i
     &      *Zd63r*Zh32/cb*sb+6D0*Abr/MW*Mb/SW*Zd33r*Zd63i
     &      *Zh32/cb*sb)
      AAABR(144) = 6D0*Sqrt2*Zd22*Zd52*Zh12*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh12/cb-3D0*EE*MW/SW*Zd22**2*Zh12/cb
     &      +3D0*EE*MW/SW*Zd22**2*Zh22*sb-EE/CW**2*MW*SW
     &      *Zd22**2*Zh12/cb+EE/CW**2*MW*SW*Zd22**2*Zh22
     &      *sb-2D0*EE/CW**2*MW*SW*Zd52**2*Zh12/cb+2D0*EE
     &      /CW**2*MW*SW*Zd52**2*Zh22*sb+3D0*EE*MW/SW*Zd22
     &      **2*Zh12/cb*sb**2+EE/CW**2*MW*SW*Zd22**2*Zh12
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd52**2*Zh12/cb
     &      *sb**2+6D0*EE/MW*Ms/SW*Zd22*Zd52*Zh32*mui-6D0
     &      *EE/MW*Ms/SW*Zd22*Zd52*Zh22/cb*mur
      AAABR(145) = 6D0*Sqrt2*Zd11*Zd41*Zh11*ds1-3D0*EE*MW
     &      /SW*Zd11**2*Zh11*cb+3D0*EE*MW/SW*Zd11**2*Zh21
     &      *sb-EE/CW**2*MW*SW*Zd11**2*Zh11*cb+EE/CW**2*MW
     &      *SW*Zd11**2*Zh21*sb-2D0*EE/CW**2*MW*SW*Zd41**2
     &      *Zh11*cb+2D0*EE/CW**2*MW*SW*Zd41**2*Zh21*sb
      AAABR(146) = EE*(3D0*MW/SW*Zd33i*Zd36i*Zh11/cb-3D0
     &      *MW/SW*Zd33i*Zd36i*Zh21*sb+3D0*MW/SW*Zd33r*Zd36r
     &      *Zh11/cb-3D0*MW/SW*Zd33r*Zd36r*Zh21*sb+MW/CW
     &      **2*SW*Zd33i*Zd36i*Zh11/cb-MW/CW**2*SW*Zd33i
     &      *Zd36i*Zh21*sb+MW/CW**2*SW*Zd33r*Zd36r*Zh11/cb
     &      -MW/CW**2*SW*Zd33r*Zd36r*Zh21*sb+2D0*MW/CW**2
     &      *SW*Zd63i*Zd66i*Zh11/cb-2D0*MW/CW**2*SW*Zd63i
     &      *Zd66i*Zh21*sb+2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh11
     &      /cb-2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh21*sb-3D0*Mb
     &      /MW/SW*Zd33i*Zd66i*Zh31*mui+3D0*Mb/MW/SW*Zd33i
     &      *Zd66r*Zh31*mur-3D0*Mb/MW/SW*Zd33r*Zd66i*Zh31
     &      *mur-3D0*Mb/MW/SW*Zd33r*Zd66r*Zh31*mui-3D0*Mb
     &      /MW/SW*Zd36i*Zd63i*Zh31*mui+3D0*Mb/MW/SW*Zd36i
     &      *Zd63r*Zh31*mur-3D0*Mb/MW/SW*Zd36r*Zd63i*Zh31
     &      *mur-3D0*Mb/MW/SW*Zd36r*Zd63r*Zh31*mui-3D0*MW
     &      /SW*Zd33i*Zd36i*Zh11/cb*sb**2-3D0*MW/SW*Zd33r
     &      *Zd36r*Zh11/cb*sb**2+3D0*Abi/MW*Mb/SW*Zd33i*Zd66r
     &      *Zh11/cb-3D0*Abi/MW*Mb/SW*Zd33r*Zd66i*Zh11/cb
     &      +3D0*Abi/MW*Mb/SW*Zd36i*Zd63r*Zh11/cb-3D0*Abi
     &      /MW*Mb/SW*Zd36r*Zd63i*Zh11/cb-3D0*Abr/MW*Mb/SW
     &      *Zd33i*Zd66i*Zh11/cb-3D0*Abr/MW*Mb/SW*Zd33r*Zd66r
     &      *Zh11/cb-3D0*Abr/MW*Mb/SW*Zd36i*Zd63i*Zh11/cb
     &      -3D0*Abr/MW*Mb/SW*Zd36r*Zd63r*Zh11/cb-MW/CW**2
     &      *SW*Zd33i*Zd36i*Zh11/cb*sb**2-MW/CW**2*SW*Zd33r
     &      *Zd36r*Zh11/cb*sb**2-2D0*MW/CW**2*SW*Zd63i*Zd66i
     &      *Zh11/cb*sb**2-2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh11
     &      /cb*sb**2+3D0*Mb/MW/SW*Zd33i*Zd66i*Zh21/cb*mur
     &      +3D0*Mb/MW/SW*Zd33i*Zd66r*Zh21/cb*mui-3D0*Mb
     &      /MW/SW*Zd33r*Zd66i*Zh21/cb*mui+3D0*Mb/MW/SW*Zd33r
     &      *Zd66r*Zh21/cb*mur+3D0*Mb/MW/SW*Zd36i*Zd63i*Zh21
     &      /cb*mur+3D0*Mb/MW/SW*Zd36i*Zd63r*Zh21/cb*mui
     &      -3D0*Mb/MW/SW*Zd36r*Zd63i*Zh21/cb*mui+3D0*Mb
     &      /MW/SW*Zd36r*Zd63r*Zh21/cb*mur+3D0*Abi/MW*Mb
     &      /SW*Zd33i*Zd66i*Zh31/cb*sb+3D0*Abi/MW*Mb/SW*Zd33r
     &      *Zd66r*Zh31/cb*sb+3D0*Abi/MW*Mb/SW*Zd36i*Zd63i
     &      *Zh31/cb*sb+3D0*Abi/MW*Mb/SW*Zd36r*Zd63r*Zh31
     &      /cb*sb+3D0*Abr/MW*Mb/SW*Zd33i*Zd66r*Zh31/cb*sb
     &      -3D0*Abr/MW*Mb/SW*Zd33r*Zd66i*Zh31/cb*sb+3D0
     &      *Abr/MW*Mb/SW*Zd36i*Zd63r*Zh31/cb*sb-3D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd63i*Zh31/cb*sb)
      AAABR(147) = EE*(3D0*MW/SW*Zd33i*Zd36r*Zh11/cb-3D0
     &      *MW/SW*Zd33i*Zd36r*Zh21*sb-3D0*MW/SW*Zd33r*Zd36i
     &      *Zh11/cb+3D0*MW/SW*Zd33r*Zd36i*Zh21*sb+MW/CW
     &      **2*SW*Zd33i*Zd36r*Zh11/cb-MW/CW**2*SW*Zd33i
     &      *Zd36r*Zh21*sb-MW/CW**2*SW*Zd33r*Zd36i*Zh11/cb
     &      +MW/CW**2*SW*Zd33r*Zd36i*Zh21*sb+2D0*MW/CW**2
     &      *SW*Zd63i*Zd66r*Zh11/cb-2D0*MW/CW**2*SW*Zd63i
     &      *Zd66r*Zh21*sb-2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh11
     &      /cb+2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh21*sb-3D0*Mb
     &      /MW/SW*Zd33i*Zd66i*Zh31*mur-3D0*Mb/MW/SW*Zd33i
     &      *Zd66r*Zh31*mui+3D0*Mb/MW/SW*Zd33r*Zd66i*Zh31
     &      *mui-3D0*Mb/MW/SW*Zd33r*Zd66r*Zh31*mur+3D0*Mb
     &      /MW/SW*Zd36i*Zd63i*Zh31*mur+3D0*Mb/MW/SW*Zd36i
     &      *Zd63r*Zh31*mui-3D0*Mb/MW/SW*Zd36r*Zd63i*Zh31
     &      *mui+3D0*Mb/MW/SW*Zd36r*Zd63r*Zh31*mur-3D0*MW
     &      /SW*Zd33i*Zd36r*Zh11/cb*sb**2+3D0*MW/SW*Zd33r
     &      *Zd36i*Zh11/cb*sb**2-3D0*Abi/MW*Mb/SW*Zd33i*Zd66i
     &      *Zh11/cb-3D0*Abi/MW*Mb/SW*Zd33r*Zd66r*Zh11/cb
     &      +3D0*Abi/MW*Mb/SW*Zd36i*Zd63i*Zh11/cb+3D0*Abi
     &      /MW*Mb/SW*Zd36r*Zd63r*Zh11/cb-3D0*Abr/MW*Mb/SW
     &      *Zd33i*Zd66r*Zh11/cb+3D0*Abr/MW*Mb/SW*Zd33r*Zd66i
     &      *Zh11/cb+3D0*Abr/MW*Mb/SW*Zd36i*Zd63r*Zh11/cb
     &      -3D0*Abr/MW*Mb/SW*Zd36r*Zd63i*Zh11/cb-MW/CW**2
     &      *SW*Zd33i*Zd36r*Zh11/cb*sb**2+MW/CW**2*SW*Zd33r
     &      *Zd36i*Zh11/cb*sb**2-2D0*MW/CW**2*SW*Zd63i*Zd66r
     &      *Zh11/cb*sb**2+2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh11
     &      /cb*sb**2-3D0*Mb/MW/SW*Zd33i*Zd66i*Zh21/cb*mui
     &      +3D0*Mb/MW/SW*Zd33i*Zd66r*Zh21/cb*mur-3D0*Mb
     &      /MW/SW*Zd33r*Zd66i*Zh21/cb*mur-3D0*Mb/MW/SW*Zd33r
     &      *Zd66r*Zh21/cb*mui+3D0*Mb/MW/SW*Zd36i*Zd63i*Zh21
     &      /cb*mui-3D0*Mb/MW/SW*Zd36i*Zd63r*Zh21/cb*mur
     &      +3D0*Mb/MW/SW*Zd36r*Zd63i*Zh21/cb*mur+3D0*Mb
     &      /MW/SW*Zd36r*Zd63r*Zh21/cb*mui+3D0*Abi/MW*Mb
     &      /SW*Zd33i*Zd66r*Zh31/cb*sb-3D0*Abi/MW*Mb/SW*Zd33r
     &      *Zd66i*Zh31/cb*sb-3D0*Abi/MW*Mb/SW*Zd36i*Zd63r
     &      *Zh31/cb*sb+3D0*Abi/MW*Mb/SW*Zd36r*Zd63i*Zh31
     &      /cb*sb-3D0*Abr/MW*Mb/SW*Zd33i*Zd66i*Zh31/cb*sb
     &      -3D0*Abr/MW*Mb/SW*Zd33r*Zd66r*Zh31/cb*sb+3D0
     &      *Abr/MW*Mb/SW*Zd36i*Zd63i*Zh31/cb*sb+3D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd63r*Zh31/cb*sb)
      AAABR(148) = 3D0*Sqrt2*Zd22*Zd55*Zh11*ds2+3D0*Sqrt2
     &      *Zd25*Zd52*Zh11*ds2-3D0*EE*MW/SW*Zd22*Zd25*Zh11
     &      /cb+3D0*EE*MW/SW*Zd22*Zd25*Zh21*sb-EE/CW**2*MW
     &      *SW*Zd22*Zd25*Zh11/cb+EE/CW**2*MW*SW*Zd22*Zd25
     &      *Zh21*sb-2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh11/cb
     &      +2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh21*sb+3D0*EE
     &      /MW*Ms/SW*Zd22*Zd55*Zh31*mui+3D0*EE/MW*Ms/SW
     &      *Zd25*Zd52*Zh31*mui+3D0*EE*MW/SW*Zd22*Zd25*Zh11
     &      /cb*sb**2+EE/CW**2*MW*SW*Zd22*Zd25*Zh11/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh11/cb*sb**2
     &      -3D0*EE/MW*Ms/SW*Zd22*Zd55*Zh21/cb*mur-3D0*EE
     &      /MW*Ms/SW*Zd25*Zd52*Zh21/cb*mur
      AAABR(149) = Sqrt2*Zd22*Zd55*Zh31*ds2*sb-Sqrt2*Zd25
     &      *Zd52*Zh31*ds2*sb+EE/MW*Ms/SW*Zd22*Zd55*Zh31
     &      *mur-EE/MW*Ms/SW*Zd25*Zd52*Zh31*mur+EE/MW*Ms
     &      /SW*Zd22*Zd55*Zh21/cb*mui-EE/MW*Ms/SW*Zd25*Zd52
     &      *Zh21/cb*mui
      AAABR(150) = 3D0*Sqrt2*Zd22*Zd55*Zh13*ds2+3D0*Sqrt2
     &      *Zd25*Zd52*Zh13*ds2-3D0*EE*MW/SW*Zd22*Zd25*Zh13
     &      /cb+3D0*EE*MW/SW*Zd22*Zd25*Zh23*sb-EE/CW**2*MW
     &      *SW*Zd22*Zd25*Zh13/cb+EE/CW**2*MW*SW*Zd22*Zd25
     &      *Zh23*sb-2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh13/cb
     &      +2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh23*sb+3D0*EE
     &      /MW*Ms/SW*Zd22*Zd55*Zh33*mui+3D0*EE/MW*Ms/SW
     &      *Zd25*Zd52*Zh33*mui+3D0*EE*MW/SW*Zd22*Zd25*Zh13
     &      /cb*sb**2+EE/CW**2*MW*SW*Zd22*Zd25*Zh13/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh13/cb*sb**2
     &      -3D0*EE/MW*Ms/SW*Zd22*Zd55*Zh23/cb*mur-3D0*EE
     &      /MW*Ms/SW*Zd25*Zd52*Zh23/cb*mur
      AAABR(151) = Sqrt2*Zd22*Zd55*Zh33*ds2*sb-Sqrt2*Zd25
     &      *Zd52*Zh33*ds2*sb+EE/MW*Ms/SW*Zd22*Zd55*Zh33
     &      *mur-EE/MW*Ms/SW*Zd25*Zd52*Zh33*mur+EE/MW*Ms
     &      /SW*Zd22*Zd55*Zh23/cb*mui-EE/MW*Ms/SW*Zd25*Zd52
     &      *Zh23/cb*mui
      AAABR(152) = EE*(3D0*MW/SW*Zd33i*Zd36i*Zh13/cb-3D0
     &      *MW/SW*Zd33i*Zd36i*Zh23*sb+3D0*MW/SW*Zd33r*Zd36r
     &      *Zh13/cb-3D0*MW/SW*Zd33r*Zd36r*Zh23*sb+MW/CW
     &      **2*SW*Zd33i*Zd36i*Zh13/cb-MW/CW**2*SW*Zd33i
     &      *Zd36i*Zh23*sb+MW/CW**2*SW*Zd33r*Zd36r*Zh13/cb
     &      -MW/CW**2*SW*Zd33r*Zd36r*Zh23*sb+2D0*MW/CW**2
     &      *SW*Zd63i*Zd66i*Zh13/cb-2D0*MW/CW**2*SW*Zd63i
     &      *Zd66i*Zh23*sb+2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh13
     &      /cb-2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh23*sb-3D0*Mb
     &      /MW/SW*Zd33i*Zd66i*Zh33*mui+3D0*Mb/MW/SW*Zd33i
     &      *Zd66r*Zh33*mur-3D0*Mb/MW/SW*Zd33r*Zd66i*Zh33
     &      *mur-3D0*Mb/MW/SW*Zd33r*Zd66r*Zh33*mui-3D0*Mb
     &      /MW/SW*Zd36i*Zd63i*Zh33*mui+3D0*Mb/MW/SW*Zd36i
     &      *Zd63r*Zh33*mur-3D0*Mb/MW/SW*Zd36r*Zd63i*Zh33
     &      *mur-3D0*Mb/MW/SW*Zd36r*Zd63r*Zh33*mui-3D0*MW
     &      /SW*Zd33i*Zd36i*Zh13/cb*sb**2-3D0*MW/SW*Zd33r
     &      *Zd36r*Zh13/cb*sb**2+3D0*Abi/MW*Mb/SW*Zd33i*Zd66r
     &      *Zh13/cb-3D0*Abi/MW*Mb/SW*Zd33r*Zd66i*Zh13/cb
     &      +3D0*Abi/MW*Mb/SW*Zd36i*Zd63r*Zh13/cb-3D0*Abi
     &      /MW*Mb/SW*Zd36r*Zd63i*Zh13/cb-3D0*Abr/MW*Mb/SW
     &      *Zd33i*Zd66i*Zh13/cb-3D0*Abr/MW*Mb/SW*Zd33r*Zd66r
     &      *Zh13/cb-3D0*Abr/MW*Mb/SW*Zd36i*Zd63i*Zh13/cb
     &      -3D0*Abr/MW*Mb/SW*Zd36r*Zd63r*Zh13/cb-MW/CW**2
     &      *SW*Zd33i*Zd36i*Zh13/cb*sb**2-MW/CW**2*SW*Zd33r
     &      *Zd36r*Zh13/cb*sb**2-2D0*MW/CW**2*SW*Zd63i*Zd66i
     &      *Zh13/cb*sb**2-2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh13
     &      /cb*sb**2+3D0*Mb/MW/SW*Zd33i*Zd66i*Zh23/cb*mur
     &      +3D0*Mb/MW/SW*Zd33i*Zd66r*Zh23/cb*mui-3D0*Mb
     &      /MW/SW*Zd33r*Zd66i*Zh23/cb*mui+3D0*Mb/MW/SW*Zd33r
     &      *Zd66r*Zh23/cb*mur+3D0*Mb/MW/SW*Zd36i*Zd63i*Zh23
     &      /cb*mur+3D0*Mb/MW/SW*Zd36i*Zd63r*Zh23/cb*mui
     &      -3D0*Mb/MW/SW*Zd36r*Zd63i*Zh23/cb*mui+3D0*Mb
     &      /MW/SW*Zd36r*Zd63r*Zh23/cb*mur+3D0*Abi/MW*Mb
     &      /SW*Zd33i*Zd66i*Zh33/cb*sb+3D0*Abi/MW*Mb/SW*Zd33r
     &      *Zd66r*Zh33/cb*sb+3D0*Abi/MW*Mb/SW*Zd36i*Zd63i
     &      *Zh33/cb*sb+3D0*Abi/MW*Mb/SW*Zd36r*Zd63r*Zh33
     &      /cb*sb+3D0*Abr/MW*Mb/SW*Zd33i*Zd66r*Zh33/cb*sb
     &      -3D0*Abr/MW*Mb/SW*Zd33r*Zd66i*Zh33/cb*sb+3D0
     &      *Abr/MW*Mb/SW*Zd36i*Zd63r*Zh33/cb*sb-3D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd63i*Zh33/cb*sb)
      AAABR(153) = EE*(3D0*MW/SW*Zd33i*Zd36r*Zh13/cb-3D0
     &      *MW/SW*Zd33i*Zd36r*Zh23*sb-3D0*MW/SW*Zd33r*Zd36i
     &      *Zh13/cb+3D0*MW/SW*Zd33r*Zd36i*Zh23*sb+MW/CW
     &      **2*SW*Zd33i*Zd36r*Zh13/cb-MW/CW**2*SW*Zd33i
     &      *Zd36r*Zh23*sb-MW/CW**2*SW*Zd33r*Zd36i*Zh13/cb
     &      +MW/CW**2*SW*Zd33r*Zd36i*Zh23*sb+2D0*MW/CW**2
     &      *SW*Zd63i*Zd66r*Zh13/cb-2D0*MW/CW**2*SW*Zd63i
     &      *Zd66r*Zh23*sb-2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh13
     &      /cb+2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh23*sb-3D0*Mb
     &      /MW/SW*Zd33i*Zd66i*Zh33*mur-3D0*Mb/MW/SW*Zd33i
     &      *Zd66r*Zh33*mui+3D0*Mb/MW/SW*Zd33r*Zd66i*Zh33
     &      *mui-3D0*Mb/MW/SW*Zd33r*Zd66r*Zh33*mur+3D0*Mb
     &      /MW/SW*Zd36i*Zd63i*Zh33*mur+3D0*Mb/MW/SW*Zd36i
     &      *Zd63r*Zh33*mui-3D0*Mb/MW/SW*Zd36r*Zd63i*Zh33
     &      *mui+3D0*Mb/MW/SW*Zd36r*Zd63r*Zh33*mur-3D0*MW
     &      /SW*Zd33i*Zd36r*Zh13/cb*sb**2+3D0*MW/SW*Zd33r
     &      *Zd36i*Zh13/cb*sb**2-3D0*Abi/MW*Mb/SW*Zd33i*Zd66i
     &      *Zh13/cb-3D0*Abi/MW*Mb/SW*Zd33r*Zd66r*Zh13/cb
     &      +3D0*Abi/MW*Mb/SW*Zd36i*Zd63i*Zh13/cb+3D0*Abi
     &      /MW*Mb/SW*Zd36r*Zd63r*Zh13/cb-3D0*Abr/MW*Mb/SW
     &      *Zd33i*Zd66r*Zh13/cb+3D0*Abr/MW*Mb/SW*Zd33r*Zd66i
     &      *Zh13/cb+3D0*Abr/MW*Mb/SW*Zd36i*Zd63r*Zh13/cb
     &      -3D0*Abr/MW*Mb/SW*Zd36r*Zd63i*Zh13/cb-MW/CW**2
     &      *SW*Zd33i*Zd36r*Zh13/cb*sb**2+MW/CW**2*SW*Zd33r
     &      *Zd36i*Zh13/cb*sb**2-2D0*MW/CW**2*SW*Zd63i*Zd66r
     &      *Zh13/cb*sb**2+2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh13
     &      /cb*sb**2-3D0*Mb/MW/SW*Zd33i*Zd66i*Zh23/cb*mui
     &      +3D0*Mb/MW/SW*Zd33i*Zd66r*Zh23/cb*mur-3D0*Mb
     &      /MW/SW*Zd33r*Zd66i*Zh23/cb*mur-3D0*Mb/MW/SW*Zd33r
     &      *Zd66r*Zh23/cb*mui+3D0*Mb/MW/SW*Zd36i*Zd63i*Zh23
     &      /cb*mui-3D0*Mb/MW/SW*Zd36i*Zd63r*Zh23/cb*mur
     &      +3D0*Mb/MW/SW*Zd36r*Zd63i*Zh23/cb*mur+3D0*Mb
     &      /MW/SW*Zd36r*Zd63r*Zh23/cb*mui+3D0*Abi/MW*Mb
     &      /SW*Zd33i*Zd66r*Zh33/cb*sb-3D0*Abi/MW*Mb/SW*Zd33r
     &      *Zd66i*Zh33/cb*sb-3D0*Abi/MW*Mb/SW*Zd36i*Zd63r
     &      *Zh33/cb*sb+3D0*Abi/MW*Mb/SW*Zd36r*Zd63i*Zh33
     &      /cb*sb-3D0*Abr/MW*Mb/SW*Zd33i*Zd66i*Zh33/cb*sb
     &      -3D0*Abr/MW*Mb/SW*Zd33r*Zd66r*Zh33/cb*sb+3D0
     &      *Abr/MW*Mb/SW*Zd36i*Zd63i*Zh33/cb*sb+3D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd63r*Zh33/cb*sb)
      AAABR(154) = 3D0*Sqrt2*Zd11*Zd44*Zh13*ds1+3D0*Sqrt2
     &      *Zd14*Zd41*Zh13*ds1-3D0*EE*MW/SW*Zd11*Zd14*Zh13
     &      *cb+3D0*EE*MW/SW*Zd11*Zd14*Zh23*sb-EE/CW**2*MW
     &      *SW*Zd11*Zd14*Zh13*cb+EE/CW**2*MW*SW*Zd11*Zd14
     &      *Zh23*sb-2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh13*cb
     &      +2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh23*sb
      AAABR(155) = Sqrt2*Zh33*ds1*sb*(Zd11*Zd44-Zd14*Zd41)
      AAABR(156) = 3D0*Sqrt2*Zd22*Zd55*Zh12*ds2+3D0*Sqrt2
     &      *Zd25*Zd52*Zh12*ds2-3D0*EE*MW/SW*Zd22*Zd25*Zh12
     &      /cb+3D0*EE*MW/SW*Zd22*Zd25*Zh22*sb-EE/CW**2*MW
     &      *SW*Zd22*Zd25*Zh12/cb+EE/CW**2*MW*SW*Zd22*Zd25
     &      *Zh22*sb-2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh12/cb
     &      +2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh22*sb+3D0*EE
     &      /MW*Ms/SW*Zd22*Zd55*Zh32*mui+3D0*EE/MW*Ms/SW
     &      *Zd25*Zd52*Zh32*mui+3D0*EE*MW/SW*Zd22*Zd25*Zh12
     &      /cb*sb**2+EE/CW**2*MW*SW*Zd22*Zd25*Zh12/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh12/cb*sb**2
     &      -3D0*EE/MW*Ms/SW*Zd22*Zd55*Zh22/cb*mur-3D0*EE
     &      /MW*Ms/SW*Zd25*Zd52*Zh22/cb*mur
      AAABR(157) = Sqrt2*Zd22*Zd55*Zh32*ds2*sb-Sqrt2*Zd25
     &      *Zd52*Zh32*ds2*sb+EE/MW*Ms/SW*Zd22*Zd55*Zh32
     &      *mur-EE/MW*Ms/SW*Zd25*Zd52*Zh32*mur+EE/MW*Ms
     &      /SW*Zd22*Zd55*Zh22/cb*mui-EE/MW*Ms/SW*Zd25*Zd52
     &      *Zh22/cb*mui
      AAABR(158) = EE*(3D0*MW/SW*Zd33i*Zd36i*Zh12/cb-3D0
     &      *MW/SW*Zd33i*Zd36i*Zh22*sb+3D0*MW/SW*Zd33r*Zd36r
     &      *Zh12/cb-3D0*MW/SW*Zd33r*Zd36r*Zh22*sb+MW/CW
     &      **2*SW*Zd33i*Zd36i*Zh12/cb-MW/CW**2*SW*Zd33i
     &      *Zd36i*Zh22*sb+MW/CW**2*SW*Zd33r*Zd36r*Zh12/cb
     &      -MW/CW**2*SW*Zd33r*Zd36r*Zh22*sb+2D0*MW/CW**2
     &      *SW*Zd63i*Zd66i*Zh12/cb-2D0*MW/CW**2*SW*Zd63i
     &      *Zd66i*Zh22*sb+2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh12
     &      /cb-2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh22*sb-3D0*Mb
     &      /MW/SW*Zd33i*Zd66i*Zh32*mui+3D0*Mb/MW/SW*Zd33i
     &      *Zd66r*Zh32*mur-3D0*Mb/MW/SW*Zd33r*Zd66i*Zh32
     &      *mur-3D0*Mb/MW/SW*Zd33r*Zd66r*Zh32*mui-3D0*Mb
     &      /MW/SW*Zd36i*Zd63i*Zh32*mui+3D0*Mb/MW/SW*Zd36i
     &      *Zd63r*Zh32*mur-3D0*Mb/MW/SW*Zd36r*Zd63i*Zh32
     &      *mur-3D0*Mb/MW/SW*Zd36r*Zd63r*Zh32*mui-3D0*MW
     &      /SW*Zd33i*Zd36i*Zh12/cb*sb**2-3D0*MW/SW*Zd33r
     &      *Zd36r*Zh12/cb*sb**2+3D0*Abi/MW*Mb/SW*Zd33i*Zd66r
     &      *Zh12/cb-3D0*Abi/MW*Mb/SW*Zd33r*Zd66i*Zh12/cb
     &      +3D0*Abi/MW*Mb/SW*Zd36i*Zd63r*Zh12/cb-3D0*Abi
     &      /MW*Mb/SW*Zd36r*Zd63i*Zh12/cb-3D0*Abr/MW*Mb/SW
     &      *Zd33i*Zd66i*Zh12/cb-3D0*Abr/MW*Mb/SW*Zd33r*Zd66r
     &      *Zh12/cb-3D0*Abr/MW*Mb/SW*Zd36i*Zd63i*Zh12/cb
     &      -3D0*Abr/MW*Mb/SW*Zd36r*Zd63r*Zh12/cb-MW/CW**2
     &      *SW*Zd33i*Zd36i*Zh12/cb*sb**2-MW/CW**2*SW*Zd33r
     &      *Zd36r*Zh12/cb*sb**2-2D0*MW/CW**2*SW*Zd63i*Zd66i
     &      *Zh12/cb*sb**2-2D0*MW/CW**2*SW*Zd63r*Zd66r*Zh12
     &      /cb*sb**2+3D0*Mb/MW/SW*Zd33i*Zd66i*Zh22/cb*mur
     &      +3D0*Mb/MW/SW*Zd33i*Zd66r*Zh22/cb*mui-3D0*Mb
     &      /MW/SW*Zd33r*Zd66i*Zh22/cb*mui+3D0*Mb/MW/SW*Zd33r
     &      *Zd66r*Zh22/cb*mur+3D0*Mb/MW/SW*Zd36i*Zd63i*Zh22
     &      /cb*mur+3D0*Mb/MW/SW*Zd36i*Zd63r*Zh22/cb*mui
     &      -3D0*Mb/MW/SW*Zd36r*Zd63i*Zh22/cb*mui+3D0*Mb
     &      /MW/SW*Zd36r*Zd63r*Zh22/cb*mur+3D0*Abi/MW*Mb
     &      /SW*Zd33i*Zd66i*Zh32/cb*sb+3D0*Abi/MW*Mb/SW*Zd33r
     &      *Zd66r*Zh32/cb*sb+3D0*Abi/MW*Mb/SW*Zd36i*Zd63i
     &      *Zh32/cb*sb+3D0*Abi/MW*Mb/SW*Zd36r*Zd63r*Zh32
     &      /cb*sb+3D0*Abr/MW*Mb/SW*Zd33i*Zd66r*Zh32/cb*sb
     &      -3D0*Abr/MW*Mb/SW*Zd33r*Zd66i*Zh32/cb*sb+3D0
     &      *Abr/MW*Mb/SW*Zd36i*Zd63r*Zh32/cb*sb-3D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd63i*Zh32/cb*sb)
      AAABR(159) = EE*(3D0*MW/SW*Zd33i*Zd36r*Zh12/cb-3D0
     &      *MW/SW*Zd33i*Zd36r*Zh22*sb-3D0*MW/SW*Zd33r*Zd36i
     &      *Zh12/cb+3D0*MW/SW*Zd33r*Zd36i*Zh22*sb+MW/CW
     &      **2*SW*Zd33i*Zd36r*Zh12/cb-MW/CW**2*SW*Zd33i
     &      *Zd36r*Zh22*sb-MW/CW**2*SW*Zd33r*Zd36i*Zh12/cb
     &      +MW/CW**2*SW*Zd33r*Zd36i*Zh22*sb+2D0*MW/CW**2
     &      *SW*Zd63i*Zd66r*Zh12/cb-2D0*MW/CW**2*SW*Zd63i
     &      *Zd66r*Zh22*sb-2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh12
     &      /cb+2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh22*sb-3D0*Mb
     &      /MW/SW*Zd33i*Zd66i*Zh32*mur-3D0*Mb/MW/SW*Zd33i
     &      *Zd66r*Zh32*mui+3D0*Mb/MW/SW*Zd33r*Zd66i*Zh32
     &      *mui-3D0*Mb/MW/SW*Zd33r*Zd66r*Zh32*mur+3D0*Mb
     &      /MW/SW*Zd36i*Zd63i*Zh32*mur+3D0*Mb/MW/SW*Zd36i
     &      *Zd63r*Zh32*mui-3D0*Mb/MW/SW*Zd36r*Zd63i*Zh32
     &      *mui+3D0*Mb/MW/SW*Zd36r*Zd63r*Zh32*mur-3D0*MW
     &      /SW*Zd33i*Zd36r*Zh12/cb*sb**2+3D0*MW/SW*Zd33r
     &      *Zd36i*Zh12/cb*sb**2-3D0*Abi/MW*Mb/SW*Zd33i*Zd66i
     &      *Zh12/cb-3D0*Abi/MW*Mb/SW*Zd33r*Zd66r*Zh12/cb
     &      +3D0*Abi/MW*Mb/SW*Zd36i*Zd63i*Zh12/cb+3D0*Abi
     &      /MW*Mb/SW*Zd36r*Zd63r*Zh12/cb-3D0*Abr/MW*Mb/SW
     &      *Zd33i*Zd66r*Zh12/cb+3D0*Abr/MW*Mb/SW*Zd33r*Zd66i
     &      *Zh12/cb+3D0*Abr/MW*Mb/SW*Zd36i*Zd63r*Zh12/cb
     &      -3D0*Abr/MW*Mb/SW*Zd36r*Zd63i*Zh12/cb-MW/CW**2
     &      *SW*Zd33i*Zd36r*Zh12/cb*sb**2+MW/CW**2*SW*Zd33r
     &      *Zd36i*Zh12/cb*sb**2-2D0*MW/CW**2*SW*Zd63i*Zd66r
     &      *Zh12/cb*sb**2+2D0*MW/CW**2*SW*Zd63r*Zd66i*Zh12
     &      /cb*sb**2-3D0*Mb/MW/SW*Zd33i*Zd66i*Zh22/cb*mui
     &      +3D0*Mb/MW/SW*Zd33i*Zd66r*Zh22/cb*mur-3D0*Mb
     &      /MW/SW*Zd33r*Zd66i*Zh22/cb*mur-3D0*Mb/MW/SW*Zd33r
     &      *Zd66r*Zh22/cb*mui+3D0*Mb/MW/SW*Zd36i*Zd63i*Zh22
     &      /cb*mui-3D0*Mb/MW/SW*Zd36i*Zd63r*Zh22/cb*mur
     &      +3D0*Mb/MW/SW*Zd36r*Zd63i*Zh22/cb*mur+3D0*Mb
     &      /MW/SW*Zd36r*Zd63r*Zh22/cb*mui+3D0*Abi/MW*Mb
     &      /SW*Zd33i*Zd66r*Zh32/cb*sb-3D0*Abi/MW*Mb/SW*Zd33r
     &      *Zd66i*Zh32/cb*sb-3D0*Abi/MW*Mb/SW*Zd36i*Zd63r
     &      *Zh32/cb*sb+3D0*Abi/MW*Mb/SW*Zd36r*Zd63i*Zh32
     &      /cb*sb-3D0*Abr/MW*Mb/SW*Zd33i*Zd66i*Zh32/cb*sb
     &      -3D0*Abr/MW*Mb/SW*Zd33r*Zd66r*Zh32/cb*sb+3D0
     &      *Abr/MW*Mb/SW*Zd36i*Zd63i*Zh32/cb*sb+3D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd63r*Zh32/cb*sb)
      AAABR(160) = 3D0*Sqrt2*Zd11*Zd44*Zh12*ds1+3D0*Sqrt2
     &      *Zd14*Zd41*Zh12*ds1-3D0*EE*MW/SW*Zd11*Zd14*Zh12
     &      *cb+3D0*EE*MW/SW*Zd11*Zd14*Zh22*sb-EE/CW**2*MW
     &      *SW*Zd11*Zd14*Zh12*cb+EE/CW**2*MW*SW*Zd11*Zd14
     &      *Zh22*sb-2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh12*cb
     &      +2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh22*sb
      AAABR(161) = Sqrt2*Zh32*ds1*sb*(Zd11*Zd44-Zd14*Zd41)
      AAABR(162) = 3D0*Sqrt2*Zd11*Zd44*Zh11*ds1+3D0*Sqrt2
     &      *Zd14*Zd41*Zh11*ds1-3D0*EE*MW/SW*Zd11*Zd14*Zh11
     &      *cb+3D0*EE*MW/SW*Zd11*Zd14*Zh21*sb-EE/CW**2*MW
     &      *SW*Zd11*Zd14*Zh11*cb+EE/CW**2*MW*SW*Zd11*Zd14
     &      *Zh21*sb-2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh11*cb
     &      +2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh21*sb
      AAABR(163) = Sqrt2*Zh31*ds1*sb*(Zd11*Zd44-Zd14*Zd41)
      AAABR(164) = 2D0*Zd14*Zu41*sb*us1+2D0*Zd44*Zu11*cb
     &      *ds1-EE*MW/SW*Sqrt2*Zd14*Zu11+2D0*EE*MW/SW*Sqrt2
     &      *Zd14*Zu11*sb**2
      AAABR(165) = 2D0*Zd25*Zu52*sb*us2+2D0*Zd55*Zu22/cb
     &      *ds2-2D0*Zd55*Zu22/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd25*Zu22-EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu22+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd25*Zu22+2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu22*sb**2+EE/MW*Mc/SW*Sqrt2*Zd25*Zu52
     &      /cb*mur/sb-EE/MW*Mc/SW*Sqrt2*Zd25*Zu52/cb*mur
     &      *sb-EE/MW*Ms/SW*Sqrt2*Zd55*Zu22/cb*mur*sb
      AAABR(166) = EE/MW/SW*Sqrt2/cb*mui*(Mc*Zd25*Zu52/sb
     &      -Mc*Zd25*Zu52*sb+Ms*Zd55*Zu22*sb)
      AAABR(167) = EE/SW*Sqrt2*(MW*Zd36i*Zu33i+MW*Zd36r*Zu33r
     &      -Mb**2/MW*Zd36i*Zu33i-Mb**2/MW*Zd36r*Zu33r+MtopS
     &      **2/MW*Zd36i*Zu33i+MtopS**2/MW*Zd36r*Zu33r-2D0
     &      *MW*Zd36i*Zu33i*sb**2-2D0*MW*Zd36r*Zu33r*sb**2
     &      -Abi/MW*Mb*Zd66i*Zu33r+Abi/MW*Mb*Zd66r*Zu33i
     &      -Abr/MW*Mb*Zd66i*Zu33i-Abr/MW*Mb*Zd66r*Zu33r
     &      -Ati/MW*MtopS*Zd36i*Zu63r+Ati/MW*MtopS*Zd36r
     &      *Zu63i+Atr/MW*MtopS*Zd36i*Zu63i+Atr/MW*MtopS
     &      *Zd36r*Zu63r+Mb/MW*Zd66i*Zu33i/cb*mur*sb-Mb/MW
     &      *Zd66i*Zu33r/cb*mui*sb+Mb/MW*Zd66r*Zu33i/cb*mui
     &      *sb+Mb/MW*Zd66r*Zu33r/cb*mur*sb-MtopS/MW*Zd36i
     &      *Zu63i/cb*mur/sb+MtopS/MW*Zd36i*Zu63i/cb*mur
     &      *sb-MtopS/MW*Zd36i*Zu63r/cb*mui/sb+MtopS/MW*Zd36i
     &      *Zu63r/cb*mui*sb+MtopS/MW*Zd36r*Zu63i/cb*mui
     &      /sb-MtopS/MW*Zd36r*Zu63i/cb*mui*sb-MtopS/MW*Zd36r
     &      *Zu63r/cb*mur/sb+MtopS/MW*Zd36r*Zu63r/cb*mur*sb)
      AAABR(168) = EE/SW*Sqrt2*(MW*Zd36i*Zu33r-MW*Zd36r*Zu33i
     &      -Mb**2/MW*Zd36i*Zu33r+Mb**2/MW*Zd36r*Zu33i+MtopS
     &      **2/MW*Zd36i*Zu33r-MtopS**2/MW*Zd36r*Zu33i-2D0
     &      *MW*Zd36i*Zu33r*sb**2+2D0*MW*Zd36r*Zu33i*sb**2
     &      +Abi/MW*Mb*Zd66i*Zu33i+Abi/MW*Mb*Zd66r*Zu33r
     &      -Abr/MW*Mb*Zd66i*Zu33r+Abr/MW*Mb*Zd66r*Zu33i
     &      +Ati/MW*MtopS*Zd36i*Zu63i+Ati/MW*MtopS*Zd36r
     &      *Zu63r+Atr/MW*MtopS*Zd36i*Zu63r-Atr/MW*MtopS
     &      *Zd36r*Zu63i+Mb/MW*Zd66i*Zu33i/cb*mui*sb+Mb/MW
     &      *Zd66i*Zu33r/cb*mur*sb-Mb/MW*Zd66r*Zu33i/cb*mur
     &      *sb+Mb/MW*Zd66r*Zu33r/cb*mui*sb+MtopS/MW*Zd36i
     &      *Zu63i/cb*mui/sb-MtopS/MW*Zd36i*Zu63i/cb*mui
     &      *sb-MtopS/MW*Zd36i*Zu63r/cb*mur/sb+MtopS/MW*Zd36i
     &      *Zu63r/cb*mur*sb+MtopS/MW*Zd36r*Zu63i/cb*mur
     &      /sb-MtopS/MW*Zd36r*Zu63i/cb*mur*sb+MtopS/MW*Zd36r
     &      *Zu63r/cb*mui/sb-MtopS/MW*Zd36r*Zu63r/cb*mui*sb)
      AAABR(169) = EE/SW*Sqrt2*(MW*Zd36i*Zu36i+MW*Zd36r*Zu36r
     &      -Mb**2/MW*Zd36i*Zu36i-Mb**2/MW*Zd36r*Zu36r+MtopS
     &      **2/MW*Zd36i*Zu36i+MtopS**2/MW*Zd36r*Zu36r-2D0
     &      *MW*Zd36i*Zu36i*sb**2-2D0*MW*Zd36r*Zu36r*sb**2
     &      -Abi/MW*Mb*Zd66i*Zu36r+Abi/MW*Mb*Zd66r*Zu36i
     &      -Abr/MW*Mb*Zd66i*Zu36i-Abr/MW*Mb*Zd66r*Zu36r
     &      -Ati/MW*MtopS*Zd36i*Zu66r+Ati/MW*MtopS*Zd36r
     &      *Zu66i+Atr/MW*MtopS*Zd36i*Zu66i+Atr/MW*MtopS
     &      *Zd36r*Zu66r+Mb/MW*Zd66i*Zu36i/cb*mur*sb-Mb/MW
     &      *Zd66i*Zu36r/cb*mui*sb+Mb/MW*Zd66r*Zu36i/cb*mui
     &      *sb+Mb/MW*Zd66r*Zu36r/cb*mur*sb-MtopS/MW*Zd36i
     &      *Zu66i/cb*mur/sb+MtopS/MW*Zd36i*Zu66i/cb*mur
     &      *sb-MtopS/MW*Zd36i*Zu66r/cb*mui/sb+MtopS/MW*Zd36i
     &      *Zu66r/cb*mui*sb+MtopS/MW*Zd36r*Zu66i/cb*mui
     &      /sb-MtopS/MW*Zd36r*Zu66i/cb*mui*sb-MtopS/MW*Zd36r
     &      *Zu66r/cb*mur/sb+MtopS/MW*Zd36r*Zu66r/cb*mur*sb)
      AAABR(170) = EE/SW*Sqrt2*(MW*Zd36i*Zu36r-MW*Zd36r*Zu36i
     &      -Mb**2/MW*Zd36i*Zu36r+Mb**2/MW*Zd36r*Zu36i+MtopS
     &      **2/MW*Zd36i*Zu36r-MtopS**2/MW*Zd36r*Zu36i-2D0
     &      *MW*Zd36i*Zu36r*sb**2+2D0*MW*Zd36r*Zu36i*sb**2
     &      +Abi/MW*Mb*Zd66i*Zu36i+Abi/MW*Mb*Zd66r*Zu36r
     &      -Abr/MW*Mb*Zd66i*Zu36r+Abr/MW*Mb*Zd66r*Zu36i
     &      +Ati/MW*MtopS*Zd36i*Zu66i+Ati/MW*MtopS*Zd36r
     &      *Zu66r+Atr/MW*MtopS*Zd36i*Zu66r-Atr/MW*MtopS
     &      *Zd36r*Zu66i+Mb/MW*Zd66i*Zu36i/cb*mui*sb+Mb/MW
     &      *Zd66i*Zu36r/cb*mur*sb-Mb/MW*Zd66r*Zu36i/cb*mur
     &      *sb+Mb/MW*Zd66r*Zu36r/cb*mui*sb+MtopS/MW*Zd36i
     &      *Zu66i/cb*mui/sb-MtopS/MW*Zd36i*Zu66i/cb*mui
     &      *sb-MtopS/MW*Zd36i*Zu66r/cb*mur/sb+MtopS/MW*Zd36i
     &      *Zu66r/cb*mur*sb+MtopS/MW*Zd36r*Zu66i/cb*mur
     &      /sb-MtopS/MW*Zd36r*Zu66i/cb*mur*sb+MtopS/MW*Zd36r
     &      *Zu66r/cb*mui/sb-MtopS/MW*Zd36r*Zu66r/cb*mui*sb)
      AAABR(171) = 2D0*Zd14*Zu44*sb*us1+2D0*Zd44*Zu14*cb
     &      *ds1-EE*MW/SW*Sqrt2*Zd14*Zu14+2D0*EE*MW/SW*Sqrt2
     &      *Zd14*Zu14*sb**2
      AAABR(172) = 2D0*Zd25*Zu55*sb*us2+2D0*Zd55*Zu25/cb
     &      *ds2-2D0*Zd55*Zu25/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd25*Zu25-EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu25+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd25*Zu25+2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu25*sb**2+EE/MW*Mc/SW*Sqrt2*Zd25*Zu55
     &      /cb*mur/sb-EE/MW*Mc/SW*Sqrt2*Zd25*Zu55/cb*mur
     &      *sb-EE/MW*Ms/SW*Sqrt2*Zd55*Zu25/cb*mur*sb
      AAABR(173) = EE/MW/SW*Sqrt2/cb*mui*(Mc*Zd25*Zu55/sb
     &      -Mc*Zd25*Zu55*sb+Ms*Zd55*Zu25*sb)
      AAABR(174) = EE/MW*Mb/SW*(Abi*Zd36i*Zd66i+Abi*Zd36r
     &      *Zd66r+Abr*Zd36i*Zd66r-Abr*Zd36r*Zd66i+Zd36i
     &      *Zd66i/cb*mui*sb-Zd36i*Zd66r/cb*mur*sb+Zd36r
     &      *Zd66i/cb*mur*sb+Zd36r*Zd66r/cb*mui*sb)
      AAABR(175) = EE/MW*Ms/SW*Zd25*Zd55*mui*tb
      AAABR(176) = EE*(6D0*Mb**2/MW/SW*Zh13/cb-3D0*MW/SW
     &      *Zd36i**2*Zh13/cb+3D0*MW/SW*Zd36i**2*Zh23*sb
     &      -3D0*MW/SW*Zd36r**2*Zh13/cb+3D0*MW/SW*Zd36r**2
     &      *Zh23*sb-MW/CW**2*SW*Zd36i**2*Zh13/cb+MW/CW**2
     &      *SW*Zd36i**2*Zh23*sb-MW/CW**2*SW*Zd36r**2*Zh13
     &      /cb+MW/CW**2*SW*Zd36r**2*Zh23*sb-2D0*MW/CW**2
     &      *SW*Zd66i**2*Zh13/cb+2D0*MW/CW**2*SW*Zd66i**2
     &      *Zh23*sb-2D0*MW/CW**2*SW*Zd66r**2*Zh13/cb+2D0
     &      *MW/CW**2*SW*Zd66r**2*Zh23*sb+3D0*MW/SW*Zd36i
     &      **2*Zh13/cb*sb**2+3D0*MW/SW*Zd36r**2*Zh13/cb
     &      *sb**2+MW/CW**2*SW*Zd36i**2*Zh13/cb*sb**2+MW
     &      /CW**2*SW*Zd36r**2*Zh13/cb*sb**2+2D0*MW/CW**2
     &      *SW*Zd66i**2*Zh13/cb*sb**2+2D0*MW/CW**2*SW*Zd66r
     &      **2*Zh13/cb*sb**2+6D0*Mb/MW/SW*Zd36i*Zd66i*Zh33
     &      *mui-6D0*Mb/MW/SW*Zd36i*Zd66r*Zh33*mur+6D0*Mb
     &      /MW/SW*Zd36r*Zd66i*Zh33*mur+6D0*Mb/MW/SW*Zd36r
     &      *Zd66r*Zh33*mui-6D0*Abi/MW*Mb/SW*Zd36i*Zd66r
     &      *Zh13/cb+6D0*Abi/MW*Mb/SW*Zd36r*Zd66i*Zh13/cb
     &      +6D0*Abr/MW*Mb/SW*Zd36i*Zd66i*Zh13/cb+6D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd66r*Zh13/cb-6D0*Mb/MW/SW*Zd36i
     &      *Zd66i*Zh23/cb*mur-6D0*Mb/MW/SW*Zd36i*Zd66r*Zh23
     &      /cb*mui+6D0*Mb/MW/SW*Zd36r*Zd66i*Zh23/cb*mui
     &      -6D0*Mb/MW/SW*Zd36r*Zd66r*Zh23/cb*mur-6D0*Abi
     &      /MW*Mb/SW*Zd36i*Zd66i*Zh33/cb*sb-6D0*Abi/MW*Mb
     &      /SW*Zd36r*Zd66r*Zh33/cb*sb-6D0*Abr/MW*Mb/SW*Zd36i
     &      *Zd66r*Zh33/cb*sb+6D0*Abr/MW*Mb/SW*Zd36r*Zd66i
     &      *Zh33/cb*sb)
      AAABR(177) = 6D0*Sqrt2*Zd25*Zd55*Zh11*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh11/cb-3D0*EE*MW/SW*Zd25**2*Zh11/cb
     &      +3D0*EE*MW/SW*Zd25**2*Zh21*sb-EE/CW**2*MW*SW
     &      *Zd25**2*Zh11/cb+EE/CW**2*MW*SW*Zd25**2*Zh21
     &      *sb-2D0*EE/CW**2*MW*SW*Zd55**2*Zh11/cb+2D0*EE
     &      /CW**2*MW*SW*Zd55**2*Zh21*sb+3D0*EE*MW/SW*Zd25
     &      **2*Zh11/cb*sb**2+EE/CW**2*MW*SW*Zd25**2*Zh11
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd55**2*Zh11/cb
     &      *sb**2+6D0*EE/MW*Ms/SW*Zd25*Zd55*Zh31*mui-6D0
     &      *EE/MW*Ms/SW*Zd25*Zd55*Zh21/cb*mur
      AAABR(178) = 6D0*Sqrt2*Zd14*Zd44*Zh13*ds1-3D0*EE*MW
     &      /SW*Zd14**2*Zh13*cb+3D0*EE*MW/SW*Zd14**2*Zh23
     &      *sb-EE/CW**2*MW*SW*Zd14**2*Zh13*cb+EE/CW**2*MW
     &      *SW*Zd14**2*Zh23*sb-2D0*EE/CW**2*MW*SW*Zd44**2
     &      *Zh13*cb+2D0*EE/CW**2*MW*SW*Zd44**2*Zh23*sb
      AAABR(179) = 6D0*Sqrt2*Zd14*Zd44*Zh12*ds1-3D0*EE*MW
     &      /SW*Zd14**2*Zh12*cb+3D0*EE*MW/SW*Zd14**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zd14**2*Zh12*cb+EE/CW**2*MW
     &      *SW*Zd14**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zd44**2
     &      *Zh12*cb+2D0*EE/CW**2*MW*SW*Zd44**2*Zh22*sb
      AAABR(180) = 6D0*Sqrt2*Zd14*Zd44*Zh11*ds1-3D0*EE*MW
     &      /SW*Zd14**2*Zh11*cb+3D0*EE*MW/SW*Zd14**2*Zh21
     &      *sb-EE/CW**2*MW*SW*Zd14**2*Zh11*cb+EE/CW**2*MW
     &      *SW*Zd14**2*Zh21*sb-2D0*EE/CW**2*MW*SW*Zd44**2
     &      *Zh11*cb+2D0*EE/CW**2*MW*SW*Zd44**2*Zh21*sb
      AAABR(181) = 6D0*Sqrt2*Zd25*Zd55*Zh12*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh12/cb-3D0*EE*MW/SW*Zd25**2*Zh12/cb
     &      +3D0*EE*MW/SW*Zd25**2*Zh22*sb-EE/CW**2*MW*SW
     &      *Zd25**2*Zh12/cb+EE/CW**2*MW*SW*Zd25**2*Zh22
     &      *sb-2D0*EE/CW**2*MW*SW*Zd55**2*Zh12/cb+2D0*EE
     &      /CW**2*MW*SW*Zd55**2*Zh22*sb+3D0*EE*MW/SW*Zd25
     &      **2*Zh12/cb*sb**2+EE/CW**2*MW*SW*Zd25**2*Zh12
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd55**2*Zh12/cb
     &      *sb**2+6D0*EE/MW*Ms/SW*Zd25*Zd55*Zh32*mui-6D0
     &      *EE/MW*Ms/SW*Zd25*Zd55*Zh22/cb*mur
      AAABR(182) = 6D0*Sqrt2*Zd25*Zd55*Zh13*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh13/cb-3D0*EE*MW/SW*Zd25**2*Zh13/cb
     &      +3D0*EE*MW/SW*Zd25**2*Zh23*sb-EE/CW**2*MW*SW
     &      *Zd25**2*Zh13/cb+EE/CW**2*MW*SW*Zd25**2*Zh23
     &      *sb-2D0*EE/CW**2*MW*SW*Zd55**2*Zh13/cb+2D0*EE
     &      /CW**2*MW*SW*Zd55**2*Zh23*sb+3D0*EE*MW/SW*Zd25
     &      **2*Zh13/cb*sb**2+EE/CW**2*MW*SW*Zd25**2*Zh13
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd55**2*Zh13/cb
     &      *sb**2+6D0*EE/MW*Ms/SW*Zd25*Zd55*Zh33*mui-6D0
     &      *EE/MW*Ms/SW*Zd25*Zd55*Zh23/cb*mur
      AAABR(183) = EE*(6D0*Mb**2/MW/SW*Zh12/cb-3D0*MW/SW
     &      *Zd36i**2*Zh12/cb+3D0*MW/SW*Zd36i**2*Zh22*sb
     &      -3D0*MW/SW*Zd36r**2*Zh12/cb+3D0*MW/SW*Zd36r**2
     &      *Zh22*sb-MW/CW**2*SW*Zd36i**2*Zh12/cb+MW/CW**2
     &      *SW*Zd36i**2*Zh22*sb-MW/CW**2*SW*Zd36r**2*Zh12
     &      /cb+MW/CW**2*SW*Zd36r**2*Zh22*sb-2D0*MW/CW**2
     &      *SW*Zd66i**2*Zh12/cb+2D0*MW/CW**2*SW*Zd66i**2
     &      *Zh22*sb-2D0*MW/CW**2*SW*Zd66r**2*Zh12/cb+2D0
     &      *MW/CW**2*SW*Zd66r**2*Zh22*sb+3D0*MW/SW*Zd36i
     &      **2*Zh12/cb*sb**2+3D0*MW/SW*Zd36r**2*Zh12/cb
     &      *sb**2+MW/CW**2*SW*Zd36i**2*Zh12/cb*sb**2+MW
     &      /CW**2*SW*Zd36r**2*Zh12/cb*sb**2+2D0*MW/CW**2
     &      *SW*Zd66i**2*Zh12/cb*sb**2+2D0*MW/CW**2*SW*Zd66r
     &      **2*Zh12/cb*sb**2+6D0*Mb/MW/SW*Zd36i*Zd66i*Zh32
     &      *mui-6D0*Mb/MW/SW*Zd36i*Zd66r*Zh32*mur+6D0*Mb
     &      /MW/SW*Zd36r*Zd66i*Zh32*mur+6D0*Mb/MW/SW*Zd36r
     &      *Zd66r*Zh32*mui-6D0*Abi/MW*Mb/SW*Zd36i*Zd66r
     &      *Zh12/cb+6D0*Abi/MW*Mb/SW*Zd36r*Zd66i*Zh12/cb
     &      +6D0*Abr/MW*Mb/SW*Zd36i*Zd66i*Zh12/cb+6D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd66r*Zh12/cb-6D0*Mb/MW/SW*Zd36i
     &      *Zd66i*Zh22/cb*mur-6D0*Mb/MW/SW*Zd36i*Zd66r*Zh22
     &      /cb*mui+6D0*Mb/MW/SW*Zd36r*Zd66i*Zh22/cb*mui
     &      -6D0*Mb/MW/SW*Zd36r*Zd66r*Zh22/cb*mur-6D0*Abi
     &      /MW*Mb/SW*Zd36i*Zd66i*Zh32/cb*sb-6D0*Abi/MW*Mb
     &      /SW*Zd36r*Zd66r*Zh32/cb*sb-6D0*Abr/MW*Mb/SW*Zd36i
     &      *Zd66r*Zh32/cb*sb+6D0*Abr/MW*Mb/SW*Zd36r*Zd66i
     &      *Zh32/cb*sb)
      AAABR(184) = EE*(6D0*Mb**2/MW/SW*Zh11/cb-3D0*MW/SW
     &      *Zd36i**2*Zh11/cb+3D0*MW/SW*Zd36i**2*Zh21*sb
     &      -3D0*MW/SW*Zd36r**2*Zh11/cb+3D0*MW/SW*Zd36r**2
     &      *Zh21*sb-MW/CW**2*SW*Zd36i**2*Zh11/cb+MW/CW**2
     &      *SW*Zd36i**2*Zh21*sb-MW/CW**2*SW*Zd36r**2*Zh11
     &      /cb+MW/CW**2*SW*Zd36r**2*Zh21*sb-2D0*MW/CW**2
     &      *SW*Zd66i**2*Zh11/cb+2D0*MW/CW**2*SW*Zd66i**2
     &      *Zh21*sb-2D0*MW/CW**2*SW*Zd66r**2*Zh11/cb+2D0
     &      *MW/CW**2*SW*Zd66r**2*Zh21*sb+3D0*MW/SW*Zd36i
     &      **2*Zh11/cb*sb**2+3D0*MW/SW*Zd36r**2*Zh11/cb
     &      *sb**2+MW/CW**2*SW*Zd36i**2*Zh11/cb*sb**2+MW
     &      /CW**2*SW*Zd36r**2*Zh11/cb*sb**2+2D0*MW/CW**2
     &      *SW*Zd66i**2*Zh11/cb*sb**2+2D0*MW/CW**2*SW*Zd66r
     &      **2*Zh11/cb*sb**2+6D0*Mb/MW/SW*Zd36i*Zd66i*Zh31
     &      *mui-6D0*Mb/MW/SW*Zd36i*Zd66r*Zh31*mur+6D0*Mb
     &      /MW/SW*Zd36r*Zd66i*Zh31*mur+6D0*Mb/MW/SW*Zd36r
     &      *Zd66r*Zh31*mui-6D0*Abi/MW*Mb/SW*Zd36i*Zd66r
     &      *Zh11/cb+6D0*Abi/MW*Mb/SW*Zd36r*Zd66i*Zh11/cb
     &      +6D0*Abr/MW*Mb/SW*Zd36i*Zd66i*Zh11/cb+6D0*Abr
     &      /MW*Mb/SW*Zd36r*Zd66r*Zh11/cb-6D0*Mb/MW/SW*Zd36i
     &      *Zd66i*Zh21/cb*mur-6D0*Mb/MW/SW*Zd36i*Zd66r*Zh21
     &      /cb*mui+6D0*Mb/MW/SW*Zd36r*Zd66i*Zh21/cb*mui
     &      -6D0*Mb/MW/SW*Zd36r*Zd66r*Zh21/cb*mur-6D0*Abi
     &      /MW*Mb/SW*Zd36i*Zd66i*Zh31/cb*sb-6D0*Abi/MW*Mb
     &      /SW*Zd36r*Zd66r*Zh31/cb*sb-6D0*Abr/MW*Mb/SW*Zd36i
     &      *Zd66r*Zh31/cb*sb+6D0*Abr/MW*Mb/SW*Zd36r*Zd66i
     &      *Zh31/cb*sb)
      AAABR(185) = EE/MW*Mc/SW*Zu22*Zu52*mui/tb
      AAABR(186) = EE/MW*MtopS/SW*(Ati*Zu33i*Zu63i+Ati*Zu33r
     &      *Zu63r+Atr*Zu33i*Zu63r-Atr*Zu33r*Zu63i+Zu33i
     &      *Zu63i*cb*mui/sb-Zu33i*Zu63r*cb*mur/sb+Zu33r
     &      *Zu63i*cb*mur/sb+Zu33r*Zu63r*cb*mui/sb)
      AAABR(187) = Sqrt2*sb*us1*(Zu11*Zu44-Zu14*Zu41)
      AAABR(188) = Sqrt2*Zu22*Zu55*sb*us2-Sqrt2*Zu25*Zu52
     &      *sb*us2+EE/MW*Mc/SW*Zu22*Zu55*cb*mur/sb-EE/MW
     &      *Mc/SW*Zu25*Zu52*cb*mur/sb
      AAABR(189) = EE/MW*Mc/SW*cb*mui/sb*(Zu22*Zu55+Zu25
     &      *Zu52)
      AAABR(190) = EE/MW*MtopS/SW*(Ati*Zu33i*Zu66r-Ati*Zu33r
     &      *Zu66i-Ati*Zu36i*Zu63r+Ati*Zu36r*Zu63i-Atr*Zu33i
     &      *Zu66i-Atr*Zu33r*Zu66r+Atr*Zu36i*Zu63i+Atr*Zu36r
     &      *Zu63r+Zu33i*Zu66i*cb*mur/sb+Zu33i*Zu66r*cb*mui
     &      /sb-Zu33r*Zu66i*cb*mui/sb+Zu33r*Zu66r*cb*mur
     &      /sb-Zu36i*Zu63i*cb*mur/sb-Zu36i*Zu63r*cb*mui
     &      /sb+Zu36r*Zu63i*cb*mui/sb-Zu36r*Zu63r*cb*mur/sb)
      AAABR(191) = EE/MW*MtopS/SW*(Ati*Zu33i*Zu66i+Ati*Zu33r
     &      *Zu66r+Ati*Zu36i*Zu63i+Ati*Zu36r*Zu63r+Atr*Zu33i
     &      *Zu66r-Atr*Zu33r*Zu66i+Atr*Zu36i*Zu63r-Atr*Zu36r
     &      *Zu63i+Zu33i*Zu66i*cb*mui/sb-Zu33i*Zu66r*cb*mur
     &      /sb+Zu33r*Zu66i*cb*mur/sb+Zu33r*Zu66r*cb*mui
     &      /sb+Zu36i*Zu63i*cb*mui/sb-Zu36i*Zu63r*cb*mur
     &      /sb+Zu36r*Zu63i*cb*mur/sb+Zu36r*Zu63r*cb*mui/sb)
      AAABR(192) = 6D0*Sqrt2*Zh21*Zu11*Zu41*us1-3D0*EE*MW
     &      /SW*Zh11*Zu11**2*cb+3D0*EE*MW/SW*Zh21*Zu11**2
     &      *sb+EE/CW**2*MW*SW*Zh11*Zu11**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh11*Zu41**2*cb-EE/CW**2*MW*SW*Zh21
     &      *Zu11**2*sb+4D0*EE/CW**2*MW*SW*Zh21*Zu41**2*sb
      AAABR(193) = EE*(6D0*MtopS**2/MW/SW*Zh21/sb+3D0*MW
     &      /SW*Zh11*Zu33i**2*cb+3D0*MW/SW*Zh11*Zu33r**2
     &      *cb-3D0*MW/SW*Zh21*Zu33i**2*sb-3D0*MW/SW*Zh21
     &      *Zu33r**2*sb-MW/CW**2*SW*Zh11*Zu33i**2*cb-MW
     &      /CW**2*SW*Zh11*Zu33r**2*cb+4D0*MW/CW**2*SW*Zh11
     &      *Zu63i**2*cb+4D0*MW/CW**2*SW*Zh11*Zu63r**2*cb
     &      +MW/CW**2*SW*Zh21*Zu33i**2*sb+MW/CW**2*SW*Zh21
     &      *Zu33r**2*sb-4D0*MW/CW**2*SW*Zh21*Zu63i**2*sb
     &      -4D0*MW/CW**2*SW*Zh21*Zu63r**2*sb+6D0*MtopS/MW
     &      /SW*Zh31*Zu33i*Zu63i*mui-6D0*MtopS/MW/SW*Zh31
     &      *Zu33i*Zu63r*mur+6D0*MtopS/MW/SW*Zh31*Zu33r*Zu63i
     &      *mur+6D0*MtopS/MW/SW*Zh31*Zu33r*Zu63r*mui-6D0
     &      *Ati/MW*MtopS/SW*Zh21*Zu33i*Zu63r/sb+6D0*Ati
     &      /MW*MtopS/SW*Zh21*Zu33r*Zu63i/sb+6D0*Atr/MW*MtopS
     &      /SW*Zh21*Zu33i*Zu63i/sb+6D0*Atr/MW*MtopS/SW*Zh21
     &      *Zu33r*Zu63r/sb-6D0*MtopS/MW/SW*Zh11*Zu33i*Zu63i
     &      *mur/sb-6D0*MtopS/MW/SW*Zh11*Zu33i*Zu63r*mui
     &      /sb+6D0*MtopS/MW/SW*Zh11*Zu33r*Zu63i*mui/sb-6D0
     &      *MtopS/MW/SW*Zh11*Zu33r*Zu63r*mur/sb-6D0*Ati
     &      /MW*MtopS/SW*Zh31*Zu33i*Zu63i*cb/sb-6D0*Ati/MW
     &      *MtopS/SW*Zh31*Zu33r*Zu63r*cb/sb-6D0*Atr/MW*MtopS
     &      /SW*Zh31*Zu33i*Zu63r*cb/sb+6D0*Atr/MW*MtopS/SW
     &      *Zh31*Zu33r*Zu63i*cb/sb)
      AAABR(194) = EE*(6D0*MtopS**2/MW/SW*Zh23/sb+3D0*MW
     &      /SW*Zh13*Zu33i**2*cb+3D0*MW/SW*Zh13*Zu33r**2
     &      *cb-3D0*MW/SW*Zh23*Zu33i**2*sb-3D0*MW/SW*Zh23
     &      *Zu33r**2*sb-MW/CW**2*SW*Zh13*Zu33i**2*cb-MW
     &      /CW**2*SW*Zh13*Zu33r**2*cb+4D0*MW/CW**2*SW*Zh13
     &      *Zu63i**2*cb+4D0*MW/CW**2*SW*Zh13*Zu63r**2*cb
     &      +MW/CW**2*SW*Zh23*Zu33i**2*sb+MW/CW**2*SW*Zh23
     &      *Zu33r**2*sb-4D0*MW/CW**2*SW*Zh23*Zu63i**2*sb
     &      -4D0*MW/CW**2*SW*Zh23*Zu63r**2*sb+6D0*MtopS/MW
     &      /SW*Zh33*Zu33i*Zu63i*mui-6D0*MtopS/MW/SW*Zh33
     &      *Zu33i*Zu63r*mur+6D0*MtopS/MW/SW*Zh33*Zu33r*Zu63i
     &      *mur+6D0*MtopS/MW/SW*Zh33*Zu33r*Zu63r*mui-6D0
     &      *Ati/MW*MtopS/SW*Zh23*Zu33i*Zu63r/sb+6D0*Ati
     &      /MW*MtopS/SW*Zh23*Zu33r*Zu63i/sb+6D0*Atr/MW*MtopS
     &      /SW*Zh23*Zu33i*Zu63i/sb+6D0*Atr/MW*MtopS/SW*Zh23
     &      *Zu33r*Zu63r/sb-6D0*MtopS/MW/SW*Zh13*Zu33i*Zu63i
     &      *mur/sb-6D0*MtopS/MW/SW*Zh13*Zu33i*Zu63r*mui
     &      /sb+6D0*MtopS/MW/SW*Zh13*Zu33r*Zu63i*mui/sb-6D0
     &      *MtopS/MW/SW*Zh13*Zu33r*Zu63r*mur/sb-6D0*Ati
     &      /MW*MtopS/SW*Zh33*Zu33i*Zu63i*cb/sb-6D0*Ati/MW
     &      *MtopS/SW*Zh33*Zu33r*Zu63r*cb/sb-6D0*Atr/MW*MtopS
     &      /SW*Zh33*Zu33i*Zu63r*cb/sb+6D0*Atr/MW*MtopS/SW
     &      *Zh33*Zu33r*Zu63i*cb/sb)
      AAABR(195) = 6D0*Sqrt2*Zh22*Zu11*Zu41*us1-3D0*EE*MW
     &      /SW*Zh12*Zu11**2*cb+3D0*EE*MW/SW*Zh22*Zu11**2
     &      *sb+EE/CW**2*MW*SW*Zh12*Zu11**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh12*Zu41**2*cb-EE/CW**2*MW*SW*Zh22
     &      *Zu11**2*sb+4D0*EE/CW**2*MW*SW*Zh22*Zu41**2*sb
      AAABR(196) = EE*(6D0*MtopS**2/MW/SW*Zh22/sb+3D0*MW
     &      /SW*Zh12*Zu33i**2*cb+3D0*MW/SW*Zh12*Zu33r**2
     &      *cb-3D0*MW/SW*Zh22*Zu33i**2*sb-3D0*MW/SW*Zh22
     &      *Zu33r**2*sb-MW/CW**2*SW*Zh12*Zu33i**2*cb-MW
     &      /CW**2*SW*Zh12*Zu33r**2*cb+4D0*MW/CW**2*SW*Zh12
     &      *Zu63i**2*cb+4D0*MW/CW**2*SW*Zh12*Zu63r**2*cb
     &      +MW/CW**2*SW*Zh22*Zu33i**2*sb+MW/CW**2*SW*Zh22
     &      *Zu33r**2*sb-4D0*MW/CW**2*SW*Zh22*Zu63i**2*sb
     &      -4D0*MW/CW**2*SW*Zh22*Zu63r**2*sb+6D0*MtopS/MW
     &      /SW*Zh32*Zu33i*Zu63i*mui-6D0*MtopS/MW/SW*Zh32
     &      *Zu33i*Zu63r*mur+6D0*MtopS/MW/SW*Zh32*Zu33r*Zu63i
     &      *mur+6D0*MtopS/MW/SW*Zh32*Zu33r*Zu63r*mui-6D0
     &      *Ati/MW*MtopS/SW*Zh22*Zu33i*Zu63r/sb+6D0*Ati
     &      /MW*MtopS/SW*Zh22*Zu33r*Zu63i/sb+6D0*Atr/MW*MtopS
     &      /SW*Zh22*Zu33i*Zu63i/sb+6D0*Atr/MW*MtopS/SW*Zh22
     &      *Zu33r*Zu63r/sb-6D0*MtopS/MW/SW*Zh12*Zu33i*Zu63i
     &      *mur/sb-6D0*MtopS/MW/SW*Zh12*Zu33i*Zu63r*mui
     &      /sb+6D0*MtopS/MW/SW*Zh12*Zu33r*Zu63i*mui/sb-6D0
     &      *MtopS/MW/SW*Zh12*Zu33r*Zu63r*mur/sb-6D0*Ati
     &      /MW*MtopS/SW*Zh32*Zu33i*Zu63i*cb/sb-6D0*Ati/MW
     &      *MtopS/SW*Zh32*Zu33r*Zu63r*cb/sb-6D0*Atr/MW*MtopS
     &      /SW*Zh32*Zu33i*Zu63r*cb/sb+6D0*Atr/MW*MtopS/SW
     &      *Zh32*Zu33r*Zu63i*cb/sb)
      AAABR(197) = 6D0*Sqrt2*Zh23*Zu22*Zu52*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh23/sb-3D0*EE*MW/SW*Zh13*Zu22**2*cb
     &      +3D0*EE*MW/SW*Zh23*Zu22**2*sb+EE/CW**2*MW*SW
     &      *Zh13*Zu22**2*cb-4D0*EE/CW**2*MW*SW*Zh13*Zu52
     &      **2*cb-EE/CW**2*MW*SW*Zh23*Zu22**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh23*Zu52**2*sb-6D0*EE/MW*Mc/SW
     &      *Zh33*Zu22*Zu52*mui+6D0*EE/MW*Mc/SW*Zh13*Zu22
     &      *Zu52*mur/sb
      AAABR(198) = 6D0*Sqrt2*Zh21*Zu22*Zu52*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh21/sb-3D0*EE*MW/SW*Zh11*Zu22**2*cb
     &      +3D0*EE*MW/SW*Zh21*Zu22**2*sb+EE/CW**2*MW*SW
     &      *Zh11*Zu22**2*cb-4D0*EE/CW**2*MW*SW*Zh11*Zu52
     &      **2*cb-EE/CW**2*MW*SW*Zh21*Zu22**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh21*Zu52**2*sb-6D0*EE/MW*Mc/SW
     &      *Zh31*Zu22*Zu52*mui+6D0*EE/MW*Mc/SW*Zh11*Zu22
     &      *Zu52*mur/sb
      AAABR(199) = 6D0*Sqrt2*Zh23*Zu11*Zu41*us1-3D0*EE*MW
     &      /SW*Zh13*Zu11**2*cb+3D0*EE*MW/SW*Zh23*Zu11**2
     &      *sb+EE/CW**2*MW*SW*Zh13*Zu11**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh13*Zu41**2*cb-EE/CW**2*MW*SW*Zh23
     &      *Zu11**2*sb+4D0*EE/CW**2*MW*SW*Zh23*Zu41**2*sb
      AAABR(200) = 6D0*Sqrt2*Zh22*Zu22*Zu52*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh22/sb-3D0*EE*MW/SW*Zh12*Zu22**2*cb
     &      +3D0*EE*MW/SW*Zh22*Zu22**2*sb+EE/CW**2*MW*SW
     &      *Zh12*Zu22**2*cb-4D0*EE/CW**2*MW*SW*Zh12*Zu52
     &      **2*cb-EE/CW**2*MW*SW*Zh22*Zu22**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh22*Zu52**2*sb-6D0*EE/MW*Mc/SW
     &      *Zh32*Zu22*Zu52*mui+6D0*EE/MW*Mc/SW*Zh12*Zu22
     &      *Zu52*mur/sb
      AAABR(201) = EE*(3D0*MW/SW*Zh13*Zu33i*Zu36i*cb+3D0
     &      *MW/SW*Zh13*Zu33r*Zu36r*cb-3D0*MW/SW*Zh23*Zu33i
     &      *Zu36i*sb-3D0*MW/SW*Zh23*Zu33r*Zu36r*sb-MW/CW
     &      **2*SW*Zh13*Zu33i*Zu36i*cb-MW/CW**2*SW*Zh13*Zu33r
     &      *Zu36r*cb+4D0*MW/CW**2*SW*Zh13*Zu63i*Zu66i*cb
     &      +4D0*MW/CW**2*SW*Zh13*Zu63r*Zu66r*cb+MW/CW**2
     &      *SW*Zh23*Zu33i*Zu36i*sb+MW/CW**2*SW*Zh23*Zu33r
     &      *Zu36r*sb-4D0*MW/CW**2*SW*Zh23*Zu63i*Zu66i*sb
     &      -4D0*MW/CW**2*SW*Zh23*Zu63r*Zu66r*sb+3D0*MtopS
     &      /MW/SW*Zh33*Zu33i*Zu66i*mui-3D0*MtopS/MW/SW*Zh33
     &      *Zu33i*Zu66r*mur+3D0*MtopS/MW/SW*Zh33*Zu33r*Zu66i
     &      *mur+3D0*MtopS/MW/SW*Zh33*Zu33r*Zu66r*mui+3D0
     &      *MtopS/MW/SW*Zh33*Zu36i*Zu63i*mui-3D0*MtopS/MW
     &      /SW*Zh33*Zu36i*Zu63r*mur+3D0*MtopS/MW/SW*Zh33
     &      *Zu36r*Zu63i*mur+3D0*MtopS/MW/SW*Zh33*Zu36r*Zu63r
     &      *mui-3D0*Ati/MW*MtopS/SW*Zh23*Zu33i*Zu66r/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh23*Zu33r*Zu66i/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh23*Zu36i*Zu63r/sb+3D0*Ati
     &      /MW*MtopS/SW*Zh23*Zu36r*Zu63i/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh23*Zu33i*Zu66i/sb+3D0*Atr/MW*MtopS/SW*Zh23
     &      *Zu33r*Zu66r/sb+3D0*Atr/MW*MtopS/SW*Zh23*Zu36i
     &      *Zu63i/sb+3D0*Atr/MW*MtopS/SW*Zh23*Zu36r*Zu63r
     &      /sb-3D0*MtopS/MW/SW*Zh13*Zu33i*Zu66i*mur/sb-3D0
     &      *MtopS/MW/SW*Zh13*Zu33i*Zu66r*mui/sb+3D0*MtopS
     &      /MW/SW*Zh13*Zu33r*Zu66i*mui/sb-3D0*MtopS/MW/SW
     &      *Zh13*Zu33r*Zu66r*mur/sb-3D0*MtopS/MW/SW*Zh13
     &      *Zu36i*Zu63i*mur/sb-3D0*MtopS/MW/SW*Zh13*Zu36i
     &      *Zu63r*mui/sb+3D0*MtopS/MW/SW*Zh13*Zu36r*Zu63i
     &      *mui/sb-3D0*MtopS/MW/SW*Zh13*Zu36r*Zu63r*mur
     &      /sb-3D0*Ati/MW*MtopS/SW*Zh33*Zu33i*Zu66i*cb/sb
     &      -3D0*Ati/MW*MtopS/SW*Zh33*Zu33r*Zu66r*cb/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh33*Zu36i*Zu63i*cb/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh33*Zu36r*Zu63r*cb/sb-3D0*Atr/MW
     &      *MtopS/SW*Zh33*Zu33i*Zu66r*cb/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh33*Zu33r*Zu66i*cb/sb-3D0*Atr/MW*MtopS/SW
     &      *Zh33*Zu36i*Zu63r*cb/sb+3D0*Atr/MW*MtopS/SW*Zh33
     &      *Zu36r*Zu63i*cb/sb)
      AAABR(202) = EE*(3D0*MW/SW*Zh13*Zu33i*Zu36r*cb-3D0
     &      *MW/SW*Zh13*Zu33r*Zu36i*cb-3D0*MW/SW*Zh23*Zu33i
     &      *Zu36r*sb+3D0*MW/SW*Zh23*Zu33r*Zu36i*sb-MW/CW
     &      **2*SW*Zh13*Zu33i*Zu36r*cb+MW/CW**2*SW*Zh13*Zu33r
     &      *Zu36i*cb+4D0*MW/CW**2*SW*Zh13*Zu63i*Zu66r*cb
     &      -4D0*MW/CW**2*SW*Zh13*Zu63r*Zu66i*cb+MW/CW**2
     &      *SW*Zh23*Zu33i*Zu36r*sb-MW/CW**2*SW*Zh23*Zu33r
     &      *Zu36i*sb-4D0*MW/CW**2*SW*Zh23*Zu63i*Zu66r*sb
     &      +4D0*MW/CW**2*SW*Zh23*Zu63r*Zu66i*sb+3D0*MtopS
     &      /MW/SW*Zh33*Zu33i*Zu66i*mur+3D0*MtopS/MW/SW*Zh33
     &      *Zu33i*Zu66r*mui-3D0*MtopS/MW/SW*Zh33*Zu33r*Zu66i
     &      *mui+3D0*MtopS/MW/SW*Zh33*Zu33r*Zu66r*mur-3D0
     &      *MtopS/MW/SW*Zh33*Zu36i*Zu63i*mur-3D0*MtopS/MW
     &      /SW*Zh33*Zu36i*Zu63r*mui+3D0*MtopS/MW/SW*Zh33
     &      *Zu36r*Zu63i*mui-3D0*MtopS/MW/SW*Zh33*Zu36r*Zu63r
     &      *mur+3D0*Ati/MW*MtopS/SW*Zh23*Zu33i*Zu66i/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh23*Zu33r*Zu66r/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh23*Zu36i*Zu63i/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh23*Zu36r*Zu63r/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh23*Zu33i*Zu66r/sb-3D0*Atr/MW*MtopS/SW*Zh23
     &      *Zu33r*Zu66i/sb-3D0*Atr/MW*MtopS/SW*Zh23*Zu36i
     &      *Zu63r/sb+3D0*Atr/MW*MtopS/SW*Zh23*Zu36r*Zu63i
     &      /sb+3D0*MtopS/MW/SW*Zh13*Zu33i*Zu66i*mui/sb-3D0
     &      *MtopS/MW/SW*Zh13*Zu33i*Zu66r*mur/sb+3D0*MtopS
     &      /MW/SW*Zh13*Zu33r*Zu66i*mur/sb+3D0*MtopS/MW/SW
     &      *Zh13*Zu33r*Zu66r*mui/sb-3D0*MtopS/MW/SW*Zh13
     &      *Zu36i*Zu63i*mui/sb+3D0*MtopS/MW/SW*Zh13*Zu36i
     &      *Zu63r*mur/sb-3D0*MtopS/MW/SW*Zh13*Zu36r*Zu63i
     &      *mur/sb-3D0*MtopS/MW/SW*Zh13*Zu36r*Zu63r*mui
     &      /sb-3D0*Ati/MW*MtopS/SW*Zh33*Zu33i*Zu66r*cb/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh33*Zu33r*Zu66i*cb/sb+3D0
     &      *Ati/MW*MtopS/SW*Zh33*Zu36i*Zu63r*cb/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh33*Zu36r*Zu63i*cb/sb+3D0*Atr/MW
     &      *MtopS/SW*Zh33*Zu33i*Zu66i*cb/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh33*Zu33r*Zu66r*cb/sb-3D0*Atr/MW*MtopS/SW
     &      *Zh33*Zu36i*Zu63i*cb/sb-3D0*Atr/MW*MtopS/SW*Zh33
     &      *Zu36r*Zu63r*cb/sb)
      AAABR(203) = 3D0*Sqrt2*Zh21*Zu11*Zu44*us1+3D0*Sqrt2
     &      *Zh21*Zu14*Zu41*us1-3D0*EE*MW/SW*Zh11*Zu11*Zu14
     &      *cb+3D0*EE*MW/SW*Zh21*Zu11*Zu14*sb+EE/CW**2*MW
     &      *SW*Zh11*Zu11*Zu14*cb-4D0*EE/CW**2*MW*SW*Zh11
     &      *Zu41*Zu44*cb-EE/CW**2*MW*SW*Zh21*Zu11*Zu14*sb
     &      +4D0*EE/CW**2*MW*SW*Zh21*Zu41*Zu44*sb
      AAABR(204) = Sqrt2*Zh31*cb*us1*(Zu11*Zu44-Zu14*Zu41)
      AAABR(205) = 3D0*Sqrt2*Zh23*Zu11*Zu44*us1+3D0*Sqrt2
     &      *Zh23*Zu14*Zu41*us1-3D0*EE*MW/SW*Zh13*Zu11*Zu14
     &      *cb+3D0*EE*MW/SW*Zh23*Zu11*Zu14*sb+EE/CW**2*MW
     &      *SW*Zh13*Zu11*Zu14*cb-4D0*EE/CW**2*MW*SW*Zh13
     &      *Zu41*Zu44*cb-EE/CW**2*MW*SW*Zh23*Zu11*Zu14*sb
     &      +4D0*EE/CW**2*MW*SW*Zh23*Zu41*Zu44*sb
      AAABR(206) = Sqrt2*Zh33*cb*us1*(Zu11*Zu44-Zu14*Zu41)
      AAABR(207) = 3D0*Sqrt2*Zh22*Zu22*Zu55*us2+3D0*Sqrt2
     &      *Zh22*Zu25*Zu52*us2-3D0*EE*MW/SW*Zh12*Zu22*Zu25
     &      *cb+3D0*EE*MW/SW*Zh22*Zu22*Zu25*sb+EE/CW**2*MW
     &      *SW*Zh12*Zu22*Zu25*cb-4D0*EE/CW**2*MW*SW*Zh12
     &      *Zu52*Zu55*cb-EE/CW**2*MW*SW*Zh22*Zu22*Zu25*sb
     &      +4D0*EE/CW**2*MW*SW*Zh22*Zu52*Zu55*sb-3D0*EE
     &      /MW*Mc/SW*Zh32*Zu22*Zu55*mui-3D0*EE/MW*Mc/SW
     &      *Zh32*Zu25*Zu52*mui+3D0*EE/MW*Mc/SW*Zh12*Zu22
     &      *Zu55*mur/sb+3D0*EE/MW*Mc/SW*Zh12*Zu25*Zu52*mur
     &      /sb
      AAABR(208) = Sqrt2*Zh32*Zu22*Zu55*cb*us2-Sqrt2*Zh32
     &      *Zu25*Zu52*cb*us2-EE/MW*Mc/SW*Zh32*Zu22*Zu55
     &      *mur+EE/MW*Mc/SW*Zh32*Zu25*Zu52*mur-EE/MW*Mc
     &      /SW*Zh12*Zu22*Zu55*mui/sb+EE/MW*Mc/SW*Zh12*Zu25
     &      *Zu52*mui/sb
      AAABR(209) = 3D0*Sqrt2*Zh22*Zu11*Zu44*us1+3D0*Sqrt2
     &      *Zh22*Zu14*Zu41*us1-3D0*EE*MW/SW*Zh12*Zu11*Zu14
     &      *cb+3D0*EE*MW/SW*Zh22*Zu11*Zu14*sb+EE/CW**2*MW
     &      *SW*Zh12*Zu11*Zu14*cb-4D0*EE/CW**2*MW*SW*Zh12
     &      *Zu41*Zu44*cb-EE/CW**2*MW*SW*Zh22*Zu11*Zu14*sb
     &      +4D0*EE/CW**2*MW*SW*Zh22*Zu41*Zu44*sb
      AAABR(210) = Sqrt2*Zh32*cb*us1*(Zu11*Zu44-Zu14*Zu41)
      AAABR(211) = EE*(3D0*MW/SW*Zh11*Zu33i*Zu36i*cb+3D0
     &      *MW/SW*Zh11*Zu33r*Zu36r*cb-3D0*MW/SW*Zh21*Zu33i
     &      *Zu36i*sb-3D0*MW/SW*Zh21*Zu33r*Zu36r*sb-MW/CW
     &      **2*SW*Zh11*Zu33i*Zu36i*cb-MW/CW**2*SW*Zh11*Zu33r
     &      *Zu36r*cb+4D0*MW/CW**2*SW*Zh11*Zu63i*Zu66i*cb
     &      +4D0*MW/CW**2*SW*Zh11*Zu63r*Zu66r*cb+MW/CW**2
     &      *SW*Zh21*Zu33i*Zu36i*sb+MW/CW**2*SW*Zh21*Zu33r
     &      *Zu36r*sb-4D0*MW/CW**2*SW*Zh21*Zu63i*Zu66i*sb
     &      -4D0*MW/CW**2*SW*Zh21*Zu63r*Zu66r*sb+3D0*MtopS
     &      /MW/SW*Zh31*Zu33i*Zu66i*mui-3D0*MtopS/MW/SW*Zh31
     &      *Zu33i*Zu66r*mur+3D0*MtopS/MW/SW*Zh31*Zu33r*Zu66i
     &      *mur+3D0*MtopS/MW/SW*Zh31*Zu33r*Zu66r*mui+3D0
     &      *MtopS/MW/SW*Zh31*Zu36i*Zu63i*mui-3D0*MtopS/MW
     &      /SW*Zh31*Zu36i*Zu63r*mur+3D0*MtopS/MW/SW*Zh31
     &      *Zu36r*Zu63i*mur+3D0*MtopS/MW/SW*Zh31*Zu36r*Zu63r
     &      *mui-3D0*Ati/MW*MtopS/SW*Zh21*Zu33i*Zu66r/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh21*Zu33r*Zu66i/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh21*Zu36i*Zu63r/sb+3D0*Ati
     &      /MW*MtopS/SW*Zh21*Zu36r*Zu63i/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh21*Zu33i*Zu66i/sb+3D0*Atr/MW*MtopS/SW*Zh21
     &      *Zu33r*Zu66r/sb+3D0*Atr/MW*MtopS/SW*Zh21*Zu36i
     &      *Zu63i/sb+3D0*Atr/MW*MtopS/SW*Zh21*Zu36r*Zu63r
     &      /sb-3D0*MtopS/MW/SW*Zh11*Zu33i*Zu66i*mur/sb-3D0
     &      *MtopS/MW/SW*Zh11*Zu33i*Zu66r*mui/sb+3D0*MtopS
     &      /MW/SW*Zh11*Zu33r*Zu66i*mui/sb-3D0*MtopS/MW/SW
     &      *Zh11*Zu33r*Zu66r*mur/sb-3D0*MtopS/MW/SW*Zh11
     &      *Zu36i*Zu63i*mur/sb-3D0*MtopS/MW/SW*Zh11*Zu36i
     &      *Zu63r*mui/sb+3D0*MtopS/MW/SW*Zh11*Zu36r*Zu63i
     &      *mui/sb-3D0*MtopS/MW/SW*Zh11*Zu36r*Zu63r*mur
     &      /sb-3D0*Ati/MW*MtopS/SW*Zh31*Zu33i*Zu66i*cb/sb
     &      -3D0*Ati/MW*MtopS/SW*Zh31*Zu33r*Zu66r*cb/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh31*Zu36i*Zu63i*cb/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh31*Zu36r*Zu63r*cb/sb-3D0*Atr/MW
     &      *MtopS/SW*Zh31*Zu33i*Zu66r*cb/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh31*Zu33r*Zu66i*cb/sb-3D0*Atr/MW*MtopS/SW
     &      *Zh31*Zu36i*Zu63r*cb/sb+3D0*Atr/MW*MtopS/SW*Zh31
     &      *Zu36r*Zu63i*cb/sb)
      AAABR(212) = EE*(3D0*MW/SW*Zh11*Zu33i*Zu36r*cb-3D0
     &      *MW/SW*Zh11*Zu33r*Zu36i*cb-3D0*MW/SW*Zh21*Zu33i
     &      *Zu36r*sb+3D0*MW/SW*Zh21*Zu33r*Zu36i*sb-MW/CW
     &      **2*SW*Zh11*Zu33i*Zu36r*cb+MW/CW**2*SW*Zh11*Zu33r
     &      *Zu36i*cb+4D0*MW/CW**2*SW*Zh11*Zu63i*Zu66r*cb
     &      -4D0*MW/CW**2*SW*Zh11*Zu63r*Zu66i*cb+MW/CW**2
     &      *SW*Zh21*Zu33i*Zu36r*sb-MW/CW**2*SW*Zh21*Zu33r
     &      *Zu36i*sb-4D0*MW/CW**2*SW*Zh21*Zu63i*Zu66r*sb
     &      +4D0*MW/CW**2*SW*Zh21*Zu63r*Zu66i*sb+3D0*MtopS
     &      /MW/SW*Zh31*Zu33i*Zu66i*mur+3D0*MtopS/MW/SW*Zh31
     &      *Zu33i*Zu66r*mui-3D0*MtopS/MW/SW*Zh31*Zu33r*Zu66i
     &      *mui+3D0*MtopS/MW/SW*Zh31*Zu33r*Zu66r*mur-3D0
     &      *MtopS/MW/SW*Zh31*Zu36i*Zu63i*mur-3D0*MtopS/MW
     &      /SW*Zh31*Zu36i*Zu63r*mui+3D0*MtopS/MW/SW*Zh31
     &      *Zu36r*Zu63i*mui-3D0*MtopS/MW/SW*Zh31*Zu36r*Zu63r
     &      *mur+3D0*Ati/MW*MtopS/SW*Zh21*Zu33i*Zu66i/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh21*Zu33r*Zu66r/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh21*Zu36i*Zu63i/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh21*Zu36r*Zu63r/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh21*Zu33i*Zu66r/sb-3D0*Atr/MW*MtopS/SW*Zh21
     &      *Zu33r*Zu66i/sb-3D0*Atr/MW*MtopS/SW*Zh21*Zu36i
     &      *Zu63r/sb+3D0*Atr/MW*MtopS/SW*Zh21*Zu36r*Zu63i
     &      /sb+3D0*MtopS/MW/SW*Zh11*Zu33i*Zu66i*mui/sb-3D0
     &      *MtopS/MW/SW*Zh11*Zu33i*Zu66r*mur/sb+3D0*MtopS
     &      /MW/SW*Zh11*Zu33r*Zu66i*mur/sb+3D0*MtopS/MW/SW
     &      *Zh11*Zu33r*Zu66r*mui/sb-3D0*MtopS/MW/SW*Zh11
     &      *Zu36i*Zu63i*mui/sb+3D0*MtopS/MW/SW*Zh11*Zu36i
     &      *Zu63r*mur/sb-3D0*MtopS/MW/SW*Zh11*Zu36r*Zu63i
     &      *mur/sb-3D0*MtopS/MW/SW*Zh11*Zu36r*Zu63r*mui
     &      /sb-3D0*Ati/MW*MtopS/SW*Zh31*Zu33i*Zu66r*cb/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh31*Zu33r*Zu66i*cb/sb+3D0
     &      *Ati/MW*MtopS/SW*Zh31*Zu36i*Zu63r*cb/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh31*Zu36r*Zu63i*cb/sb+3D0*Atr/MW
     &      *MtopS/SW*Zh31*Zu33i*Zu66i*cb/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh31*Zu33r*Zu66r*cb/sb-3D0*Atr/MW*MtopS/SW
     &      *Zh31*Zu36i*Zu63i*cb/sb-3D0*Atr/MW*MtopS/SW*Zh31
     &      *Zu36r*Zu63r*cb/sb)
      AAABR(213) = 3D0*Sqrt2*Zh21*Zu22*Zu55*us2+3D0*Sqrt2
     &      *Zh21*Zu25*Zu52*us2-3D0*EE*MW/SW*Zh11*Zu22*Zu25
     &      *cb+3D0*EE*MW/SW*Zh21*Zu22*Zu25*sb+EE/CW**2*MW
     &      *SW*Zh11*Zu22*Zu25*cb-4D0*EE/CW**2*MW*SW*Zh11
     &      *Zu52*Zu55*cb-EE/CW**2*MW*SW*Zh21*Zu22*Zu25*sb
     &      +4D0*EE/CW**2*MW*SW*Zh21*Zu52*Zu55*sb-3D0*EE
     &      /MW*Mc/SW*Zh31*Zu22*Zu55*mui-3D0*EE/MW*Mc/SW
     &      *Zh31*Zu25*Zu52*mui+3D0*EE/MW*Mc/SW*Zh11*Zu22
     &      *Zu55*mur/sb+3D0*EE/MW*Mc/SW*Zh11*Zu25*Zu52*mur
     &      /sb
      AAABR(214) = Sqrt2*Zh31*Zu22*Zu55*cb*us2-Sqrt2*Zh31
     &      *Zu25*Zu52*cb*us2-EE/MW*Mc/SW*Zh31*Zu22*Zu55
     &      *mur+EE/MW*Mc/SW*Zh31*Zu25*Zu52*mur-EE/MW*Mc
     &      /SW*Zh11*Zu22*Zu55*mui/sb+EE/MW*Mc/SW*Zh11*Zu25
     &      *Zu52*mui/sb
      AAABR(215) = EE*(3D0*MW/SW*Zh12*Zu33i*Zu36i*cb+3D0
     &      *MW/SW*Zh12*Zu33r*Zu36r*cb-3D0*MW/SW*Zh22*Zu33i
     &      *Zu36i*sb-3D0*MW/SW*Zh22*Zu33r*Zu36r*sb-MW/CW
     &      **2*SW*Zh12*Zu33i*Zu36i*cb-MW/CW**2*SW*Zh12*Zu33r
     &      *Zu36r*cb+4D0*MW/CW**2*SW*Zh12*Zu63i*Zu66i*cb
     &      +4D0*MW/CW**2*SW*Zh12*Zu63r*Zu66r*cb+MW/CW**2
     &      *SW*Zh22*Zu33i*Zu36i*sb+MW/CW**2*SW*Zh22*Zu33r
     &      *Zu36r*sb-4D0*MW/CW**2*SW*Zh22*Zu63i*Zu66i*sb
     &      -4D0*MW/CW**2*SW*Zh22*Zu63r*Zu66r*sb+3D0*MtopS
     &      /MW/SW*Zh32*Zu33i*Zu66i*mui-3D0*MtopS/MW/SW*Zh32
     &      *Zu33i*Zu66r*mur+3D0*MtopS/MW/SW*Zh32*Zu33r*Zu66i
     &      *mur+3D0*MtopS/MW/SW*Zh32*Zu33r*Zu66r*mui+3D0
     &      *MtopS/MW/SW*Zh32*Zu36i*Zu63i*mui-3D0*MtopS/MW
     &      /SW*Zh32*Zu36i*Zu63r*mur+3D0*MtopS/MW/SW*Zh32
     &      *Zu36r*Zu63i*mur+3D0*MtopS/MW/SW*Zh32*Zu36r*Zu63r
     &      *mui-3D0*Ati/MW*MtopS/SW*Zh22*Zu33i*Zu66r/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh22*Zu33r*Zu66i/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh22*Zu36i*Zu63r/sb+3D0*Ati
     &      /MW*MtopS/SW*Zh22*Zu36r*Zu63i/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh22*Zu33i*Zu66i/sb+3D0*Atr/MW*MtopS/SW*Zh22
     &      *Zu33r*Zu66r/sb+3D0*Atr/MW*MtopS/SW*Zh22*Zu36i
     &      *Zu63i/sb+3D0*Atr/MW*MtopS/SW*Zh22*Zu36r*Zu63r
     &      /sb-3D0*MtopS/MW/SW*Zh12*Zu33i*Zu66i*mur/sb-3D0
     &      *MtopS/MW/SW*Zh12*Zu33i*Zu66r*mui/sb+3D0*MtopS
     &      /MW/SW*Zh12*Zu33r*Zu66i*mui/sb-3D0*MtopS/MW/SW
     &      *Zh12*Zu33r*Zu66r*mur/sb-3D0*MtopS/MW/SW*Zh12
     &      *Zu36i*Zu63i*mur/sb-3D0*MtopS/MW/SW*Zh12*Zu36i
     &      *Zu63r*mui/sb+3D0*MtopS/MW/SW*Zh12*Zu36r*Zu63i
     &      *mui/sb-3D0*MtopS/MW/SW*Zh12*Zu36r*Zu63r*mur
     &      /sb-3D0*Ati/MW*MtopS/SW*Zh32*Zu33i*Zu66i*cb/sb
     &      -3D0*Ati/MW*MtopS/SW*Zh32*Zu33r*Zu66r*cb/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh32*Zu36i*Zu63i*cb/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh32*Zu36r*Zu63r*cb/sb-3D0*Atr/MW
     &      *MtopS/SW*Zh32*Zu33i*Zu66r*cb/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh32*Zu33r*Zu66i*cb/sb-3D0*Atr/MW*MtopS/SW
     &      *Zh32*Zu36i*Zu63r*cb/sb+3D0*Atr/MW*MtopS/SW*Zh32
     &      *Zu36r*Zu63i*cb/sb)
      AAABR(216) = EE*(3D0*MW/SW*Zh12*Zu33i*Zu36r*cb-3D0
     &      *MW/SW*Zh12*Zu33r*Zu36i*cb-3D0*MW/SW*Zh22*Zu33i
     &      *Zu36r*sb+3D0*MW/SW*Zh22*Zu33r*Zu36i*sb-MW/CW
     &      **2*SW*Zh12*Zu33i*Zu36r*cb+MW/CW**2*SW*Zh12*Zu33r
     &      *Zu36i*cb+4D0*MW/CW**2*SW*Zh12*Zu63i*Zu66r*cb
     &      -4D0*MW/CW**2*SW*Zh12*Zu63r*Zu66i*cb+MW/CW**2
     &      *SW*Zh22*Zu33i*Zu36r*sb-MW/CW**2*SW*Zh22*Zu33r
     &      *Zu36i*sb-4D0*MW/CW**2*SW*Zh22*Zu63i*Zu66r*sb
     &      +4D0*MW/CW**2*SW*Zh22*Zu63r*Zu66i*sb+3D0*MtopS
     &      /MW/SW*Zh32*Zu33i*Zu66i*mur+3D0*MtopS/MW/SW*Zh32
     &      *Zu33i*Zu66r*mui-3D0*MtopS/MW/SW*Zh32*Zu33r*Zu66i
     &      *mui+3D0*MtopS/MW/SW*Zh32*Zu33r*Zu66r*mur-3D0
     &      *MtopS/MW/SW*Zh32*Zu36i*Zu63i*mur-3D0*MtopS/MW
     &      /SW*Zh32*Zu36i*Zu63r*mui+3D0*MtopS/MW/SW*Zh32
     &      *Zu36r*Zu63i*mui-3D0*MtopS/MW/SW*Zh32*Zu36r*Zu63r
     &      *mur+3D0*Ati/MW*MtopS/SW*Zh22*Zu33i*Zu66i/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh22*Zu33r*Zu66r/sb-3D0
     &      *Ati/MW*MtopS/SW*Zh22*Zu36i*Zu63i/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh22*Zu36r*Zu63r/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh22*Zu33i*Zu66r/sb-3D0*Atr/MW*MtopS/SW*Zh22
     &      *Zu33r*Zu66i/sb-3D0*Atr/MW*MtopS/SW*Zh22*Zu36i
     &      *Zu63r/sb+3D0*Atr/MW*MtopS/SW*Zh22*Zu36r*Zu63i
     &      /sb+3D0*MtopS/MW/SW*Zh12*Zu33i*Zu66i*mui/sb-3D0
     &      *MtopS/MW/SW*Zh12*Zu33i*Zu66r*mur/sb+3D0*MtopS
     &      /MW/SW*Zh12*Zu33r*Zu66i*mur/sb+3D0*MtopS/MW/SW
     &      *Zh12*Zu33r*Zu66r*mui/sb-3D0*MtopS/MW/SW*Zh12
     &      *Zu36i*Zu63i*mui/sb+3D0*MtopS/MW/SW*Zh12*Zu36i
     &      *Zu63r*mur/sb-3D0*MtopS/MW/SW*Zh12*Zu36r*Zu63i
     &      *mur/sb-3D0*MtopS/MW/SW*Zh12*Zu36r*Zu63r*mui
     &      /sb-3D0*Ati/MW*MtopS/SW*Zh32*Zu33i*Zu66r*cb/sb
     &      +3D0*Ati/MW*MtopS/SW*Zh32*Zu33r*Zu66i*cb/sb+3D0
     &      *Ati/MW*MtopS/SW*Zh32*Zu36i*Zu63r*cb/sb-3D0*Ati
     &      /MW*MtopS/SW*Zh32*Zu36r*Zu63i*cb/sb+3D0*Atr/MW
     &      *MtopS/SW*Zh32*Zu33i*Zu66i*cb/sb+3D0*Atr/MW*MtopS
     &      /SW*Zh32*Zu33r*Zu66r*cb/sb-3D0*Atr/MW*MtopS/SW
     &      *Zh32*Zu36i*Zu63i*cb/sb-3D0*Atr/MW*MtopS/SW*Zh32
     &      *Zu36r*Zu63r*cb/sb)
      AAABR(217) = 3D0*Sqrt2*Zh23*Zu22*Zu55*us2+3D0*Sqrt2
     &      *Zh23*Zu25*Zu52*us2-3D0*EE*MW/SW*Zh13*Zu22*Zu25
     &      *cb+3D0*EE*MW/SW*Zh23*Zu22*Zu25*sb+EE/CW**2*MW
     &      *SW*Zh13*Zu22*Zu25*cb-4D0*EE/CW**2*MW*SW*Zh13
     &      *Zu52*Zu55*cb-EE/CW**2*MW*SW*Zh23*Zu22*Zu25*sb
     &      +4D0*EE/CW**2*MW*SW*Zh23*Zu52*Zu55*sb-3D0*EE
     &      /MW*Mc/SW*Zh33*Zu22*Zu55*mui-3D0*EE/MW*Mc/SW
     &      *Zh33*Zu25*Zu52*mui+3D0*EE/MW*Mc/SW*Zh13*Zu22
     &      *Zu55*mur/sb+3D0*EE/MW*Mc/SW*Zh13*Zu25*Zu52*mur
     &      /sb
      AAABR(218) = Sqrt2*Zh33*Zu22*Zu55*cb*us2-Sqrt2*Zh33
     &      *Zu25*Zu52*cb*us2-EE/MW*Mc/SW*Zh33*Zu22*Zu55
     &      *mur+EE/MW*Mc/SW*Zh33*Zu25*Zu52*mur-EE/MW*Mc
     &      /SW*Zh13*Zu22*Zu55*mui/sb+EE/MW*Mc/SW*Zh13*Zu25
     &      *Zu52*mui/sb
      AAABR(219) = EE/MW*MtopS/SW*(Ati*Zu36i*Zu66i+Ati*Zu36r
     &      *Zu66r+Atr*Zu36i*Zu66r-Atr*Zu36r*Zu66i+Zu36i
     &      *Zu66i*cb*mui/sb-Zu36i*Zu66r*cb*mur/sb+Zu36r
     &      *Zu66i*cb*mur/sb+Zu36r*Zu66r*cb*mui/sb)
      AAABR(220) = EE/MW*Mc/SW*Zu25*Zu55*mui/tb
      AAABR(221) = 6D0*Sqrt2*Zh22*Zu14*Zu44*us1-3D0*EE*MW
     &      /SW*Zh12*Zu14**2*cb+3D0*EE*MW/SW*Zh22*Zu14**2
     &      *sb+EE/CW**2*MW*SW*Zh12*Zu14**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh12*Zu44**2*cb-EE/CW**2*MW*SW*Zh22
     &      *Zu14**2*sb+4D0*EE/CW**2*MW*SW*Zh22*Zu44**2*sb
      AAABR(222) = 6D0*Sqrt2*Zh21*Zu25*Zu55*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh21/sb-3D0*EE*MW/SW*Zh11*Zu25**2*cb
     &      +3D0*EE*MW/SW*Zh21*Zu25**2*sb+EE/CW**2*MW*SW
     &      *Zh11*Zu25**2*cb-4D0*EE/CW**2*MW*SW*Zh11*Zu55
     &      **2*cb-EE/CW**2*MW*SW*Zh21*Zu25**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh21*Zu55**2*sb-6D0*EE/MW*Mc/SW
     &      *Zh31*Zu25*Zu55*mui+6D0*EE/MW*Mc/SW*Zh11*Zu25
     &      *Zu55*mur/sb
      AAABR(223) = 6D0*Sqrt2*Zh21*Zu14*Zu44*us1-3D0*EE*MW
     &      /SW*Zh11*Zu14**2*cb+3D0*EE*MW/SW*Zh21*Zu14**2
     &      *sb+EE/CW**2*MW*SW*Zh11*Zu14**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh11*Zu44**2*cb-EE/CW**2*MW*SW*Zh21
     &      *Zu14**2*sb+4D0*EE/CW**2*MW*SW*Zh21*Zu44**2*sb
      AAABR(224) = 6D0*Sqrt2*Zh23*Zu25*Zu55*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh23/sb-3D0*EE*MW/SW*Zh13*Zu25**2*cb
     &      +3D0*EE*MW/SW*Zh23*Zu25**2*sb+EE/CW**2*MW*SW
     &      *Zh13*Zu25**2*cb-4D0*EE/CW**2*MW*SW*Zh13*Zu55
     &      **2*cb-EE/CW**2*MW*SW*Zh23*Zu25**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh23*Zu55**2*sb-6D0*EE/MW*Mc/SW
     &      *Zh33*Zu25*Zu55*mui+6D0*EE/MW*Mc/SW*Zh13*Zu25
     &      *Zu55*mur/sb
      AAABR(225) = 6D0*Sqrt2*Zh23*Zu14*Zu44*us1-3D0*EE*MW
     &      /SW*Zh13*Zu14**2*cb+3D0*EE*MW/SW*Zh23*Zu14**2
     &      *sb+EE/CW**2*MW*SW*Zh13*Zu14**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh13*Zu44**2*cb-EE/CW**2*MW*SW*Zh23
     &      *Zu14**2*sb+4D0*EE/CW**2*MW*SW*Zh23*Zu44**2*sb
      AAABR(226) = EE*(6D0*MtopS**2/MW/SW*Zh22/sb+3D0*MW
     &      /SW*Zh12*Zu36i**2*cb+3D0*MW/SW*Zh12*Zu36r**2
     &      *cb-3D0*MW/SW*Zh22*Zu36i**2*sb-3D0*MW/SW*Zh22
     &      *Zu36r**2*sb-MW/CW**2*SW*Zh12*Zu36i**2*cb-MW
     &      /CW**2*SW*Zh12*Zu36r**2*cb+4D0*MW/CW**2*SW*Zh12
     &      *Zu66i**2*cb+4D0*MW/CW**2*SW*Zh12*Zu66r**2*cb
     &      +MW/CW**2*SW*Zh22*Zu36i**2*sb+MW/CW**2*SW*Zh22
     &      *Zu36r**2*sb-4D0*MW/CW**2*SW*Zh22*Zu66i**2*sb
     &      -4D0*MW/CW**2*SW*Zh22*Zu66r**2*sb+6D0*MtopS/MW
     &      /SW*Zh32*Zu36i*Zu66i*mui-6D0*MtopS/MW/SW*Zh32
     &      *Zu36i*Zu66r*mur+6D0*MtopS/MW/SW*Zh32*Zu36r*Zu66i
     &      *mur+6D0*MtopS/MW/SW*Zh32*Zu36r*Zu66r*mui-6D0
     &      *Ati/MW*MtopS/SW*Zh22*Zu36i*Zu66r/sb+6D0*Ati
     &      /MW*MtopS/SW*Zh22*Zu36r*Zu66i/sb+6D0*Atr/MW*MtopS
     &      /SW*Zh22*Zu36i*Zu66i/sb+6D0*Atr/MW*MtopS/SW*Zh22
     &      *Zu36r*Zu66r/sb-6D0*MtopS/MW/SW*Zh12*Zu36i*Zu66i
     &      *mur/sb-6D0*MtopS/MW/SW*Zh12*Zu36i*Zu66r*mui
     &      /sb+6D0*MtopS/MW/SW*Zh12*Zu36r*Zu66i*mui/sb-6D0
     &      *MtopS/MW/SW*Zh12*Zu36r*Zu66r*mur/sb-6D0*Ati
     &      /MW*MtopS/SW*Zh32*Zu36i*Zu66i*cb/sb-6D0*Ati/MW
     &      *MtopS/SW*Zh32*Zu36r*Zu66r*cb/sb-6D0*Atr/MW*MtopS
     &      /SW*Zh32*Zu36i*Zu66r*cb/sb+6D0*Atr/MW*MtopS/SW
     &      *Zh32*Zu36r*Zu66i*cb/sb)
      AAABR(227) = EE*(6D0*MtopS**2/MW/SW*Zh21/sb+3D0*MW
     &      /SW*Zh11*Zu36i**2*cb+3D0*MW/SW*Zh11*Zu36r**2
     &      *cb-3D0*MW/SW*Zh21*Zu36i**2*sb-3D0*MW/SW*Zh21
     &      *Zu36r**2*sb-MW/CW**2*SW*Zh11*Zu36i**2*cb-MW
     &      /CW**2*SW*Zh11*Zu36r**2*cb+4D0*MW/CW**2*SW*Zh11
     &      *Zu66i**2*cb+4D0*MW/CW**2*SW*Zh11*Zu66r**2*cb
     &      +MW/CW**2*SW*Zh21*Zu36i**2*sb+MW/CW**2*SW*Zh21
     &      *Zu36r**2*sb-4D0*MW/CW**2*SW*Zh21*Zu66i**2*sb
     &      -4D0*MW/CW**2*SW*Zh21*Zu66r**2*sb+6D0*MtopS/MW
     &      /SW*Zh31*Zu36i*Zu66i*mui-6D0*MtopS/MW/SW*Zh31
     &      *Zu36i*Zu66r*mur+6D0*MtopS/MW/SW*Zh31*Zu36r*Zu66i
     &      *mur+6D0*MtopS/MW/SW*Zh31*Zu36r*Zu66r*mui-6D0
     &      *Ati/MW*MtopS/SW*Zh21*Zu36i*Zu66r/sb+6D0*Ati
     &      /MW*MtopS/SW*Zh21*Zu36r*Zu66i/sb+6D0*Atr/MW*MtopS
     &      /SW*Zh21*Zu36i*Zu66i/sb+6D0*Atr/MW*MtopS/SW*Zh21
     &      *Zu36r*Zu66r/sb-6D0*MtopS/MW/SW*Zh11*Zu36i*Zu66i
     &      *mur/sb-6D0*MtopS/MW/SW*Zh11*Zu36i*Zu66r*mui
     &      /sb+6D0*MtopS/MW/SW*Zh11*Zu36r*Zu66i*mui/sb-6D0
     &      *MtopS/MW/SW*Zh11*Zu36r*Zu66r*mur/sb-6D0*Ati
     &      /MW*MtopS/SW*Zh31*Zu36i*Zu66i*cb/sb-6D0*Ati/MW
     &      *MtopS/SW*Zh31*Zu36r*Zu66r*cb/sb-6D0*Atr/MW*MtopS
     &      /SW*Zh31*Zu36i*Zu66r*cb/sb+6D0*Atr/MW*MtopS/SW
     &      *Zh31*Zu36r*Zu66i*cb/sb)
      AAABR(228) = 6D0*Sqrt2*Zh22*Zu25*Zu55*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh22/sb-3D0*EE*MW/SW*Zh12*Zu25**2*cb
     &      +3D0*EE*MW/SW*Zh22*Zu25**2*sb+EE/CW**2*MW*SW
     &      *Zh12*Zu25**2*cb-4D0*EE/CW**2*MW*SW*Zh12*Zu55
     &      **2*cb-EE/CW**2*MW*SW*Zh22*Zu25**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh22*Zu55**2*sb-6D0*EE/MW*Mc/SW
     &      *Zh32*Zu25*Zu55*mui+6D0*EE/MW*Mc/SW*Zh12*Zu25
     &      *Zu55*mur/sb
      AAABR(229) = EE*(6D0*MtopS**2/MW/SW*Zh23/sb+3D0*MW
     &      /SW*Zh13*Zu36i**2*cb+3D0*MW/SW*Zh13*Zu36r**2
     &      *cb-3D0*MW/SW*Zh23*Zu36i**2*sb-3D0*MW/SW*Zh23
     &      *Zu36r**2*sb-MW/CW**2*SW*Zh13*Zu36i**2*cb-MW
     &      /CW**2*SW*Zh13*Zu36r**2*cb+4D0*MW/CW**2*SW*Zh13
     &      *Zu66i**2*cb+4D0*MW/CW**2*SW*Zh13*Zu66r**2*cb
     &      +MW/CW**2*SW*Zh23*Zu36i**2*sb+MW/CW**2*SW*Zh23
     &      *Zu36r**2*sb-4D0*MW/CW**2*SW*Zh23*Zu66i**2*sb
     &      -4D0*MW/CW**2*SW*Zh23*Zu66r**2*sb+6D0*MtopS/MW
     &      /SW*Zh33*Zu36i*Zu66i*mui-6D0*MtopS/MW/SW*Zh33
     &      *Zu36i*Zu66r*mur+6D0*MtopS/MW/SW*Zh33*Zu36r*Zu66i
     &      *mur+6D0*MtopS/MW/SW*Zh33*Zu36r*Zu66r*mui-6D0
     &      *Ati/MW*MtopS/SW*Zh23*Zu36i*Zu66r/sb+6D0*Ati
     &      /MW*MtopS/SW*Zh23*Zu36r*Zu66i/sb+6D0*Atr/MW*MtopS
     &      /SW*Zh23*Zu36i*Zu66i/sb+6D0*Atr/MW*MtopS/SW*Zh23
     &      *Zu36r*Zu66r/sb-6D0*MtopS/MW/SW*Zh13*Zu36i*Zu66i
     &      *mur/sb-6D0*MtopS/MW/SW*Zh13*Zu36i*Zu66r*mui
     &      /sb+6D0*MtopS/MW/SW*Zh13*Zu36r*Zu66i*mui/sb-6D0
     &      *MtopS/MW/SW*Zh13*Zu36r*Zu66r*mur/sb-6D0*Ati
     &      /MW*MtopS/SW*Zh33*Zu36i*Zu66i*cb/sb-6D0*Ati/MW
     &      *MtopS/SW*Zh33*Zu36r*Zu66r*cb/sb-6D0*Atr/MW*MtopS
     &      /SW*Zh33*Zu36i*Zu66r*cb/sb+6D0*Atr/MW*MtopS/SW
     &      *Zh33*Zu36r*Zu66i*cb/sb)
      AAABR(230) = EE*MW*(Zh12/SW*cb-Zh22/SW*sb+2D0*Zh22
     &      /SW*sb**3+SW/CW**2*Zh12*cb-SW/CW**2*Zh22*sb+2D0
     &      *SW/CW**2*Zh22*sb**3-2D0*Zh12/SW*cb*sb**2-2D0
     &      *SW/CW**2*Zh12*cb*sb**2)
      AAABR(231) = EE*MW*(Zh13/SW*cb-Zh23/SW*sb+2D0*Zh23
     &      /SW*sb**3+SW/CW**2*Zh13*cb-SW/CW**2*Zh23*sb+2D0
     &      *SW/CW**2*Zh23*sb**3-2D0*Zh13/SW*cb*sb**2-2D0
     &      *SW/CW**2*Zh13*cb*sb**2)
      AAABR(232) = EE*MW*(Zh11/SW*cb-Zh21/SW*sb+2D0*Zh21
     &      /SW*sb**3+SW/CW**2*Zh11*cb-SW/CW**2*Zh21*sb+2D0
     &      *SW/CW**2*Zh21*sb**3-2D0*Zh11/SW*cb*sb**2-2D0
     &      *SW/CW**2*Zh11*cb*sb**2)
      AAABR(233) = EE*MW*(Zh12/SW*Zh33*cb**2*sb+Zh13/SW*Zh32
     &      *cb**2*sb-Zh22/SW*Zh33*cb*sb**2-Zh23/SW*Zh32
     &      *cb*sb**2+SW/CW**2*Zh12*Zh33*cb**2*sb+SW/CW**2
     &      *Zh13*Zh32*cb**2*sb-SW/CW**2*Zh22*Zh33*cb*sb
     &      **2-SW/CW**2*Zh23*Zh32*cb*sb**2)
      AAABR(234) = EE*MW*(Zh11/SW*Zh32*cb**2*sb+Zh12/SW*Zh31
     &      *cb**2*sb-Zh21/SW*Zh32*cb*sb**2-Zh22/SW*Zh31
     &      *cb*sb**2+SW/CW**2*Zh11*Zh32*cb**2*sb+SW/CW**2
     &      *Zh12*Zh31*cb**2*sb-SW/CW**2*Zh21*Zh32*cb*sb
     &      **2-SW/CW**2*Zh22*Zh31*cb*sb**2)
      AAABR(235) = EE*MW*Zh31*(Zh11/SW*cb**2*sb-Zh21/SW*cb
     &      *sb**2+SW/CW**2*Zh11*cb**2*sb-SW/CW**2*Zh21*cb
     &      *sb**2)
      AAABR(236) = EE*MW*Zh33*(Zh13/SW*cb**2*sb-Zh23/SW*cb
     &      *sb**2+SW/CW**2*Zh13*cb**2*sb-SW/CW**2*Zh23*cb
     &      *sb**2)
      AAABR(237) = EE*MW*(Zh11/SW*Zh33*cb**2*sb+Zh13/SW*Zh31
     &      *cb**2*sb-Zh21/SW*Zh33*cb*sb**2-Zh23/SW*Zh31
     &      *cb*sb**2+SW/CW**2*Zh11*Zh33*cb**2*sb+SW/CW**2
     &      *Zh13*Zh31*cb**2*sb-SW/CW**2*Zh21*Zh33*cb*sb
     &      **2-SW/CW**2*Zh23*Zh31*cb*sb**2)
      AAABR(238) = EE*MW*Zh32*(Zh12/SW*cb**2*sb-Zh22/SW*cb
     &      *sb**2+SW/CW**2*Zh12*cb**2*sb-SW/CW**2*Zh22*cb
     &      *sb**2)
      AAABR(239) = EE*MW*(3D0*Zh11**2/SW*Zh12*cb-Zh11**2
     &      /SW*Zh22*sb-Zh12/SW*Zh21**2*cb-Zh12/SW*Zh31**2
     &      *cb+3D0*Zh21**2/SW*Zh22*sb+Zh22/SW*Zh31**2*sb
     &      -2D0*Zh22/SW*Zh31**2*sb**3+3D0*SW/CW**2*Zh11
     &      **2*Zh12*cb-SW/CW**2*Zh11**2*Zh22*sb-SW/CW**2
     &      *Zh12*Zh21**2*cb-SW/CW**2*Zh12*Zh31**2*cb+3D0
     &      *SW/CW**2*Zh21**2*Zh22*sb+SW/CW**2*Zh22*Zh31
     &      **2*sb-2D0*SW/CW**2*Zh22*Zh31**2*sb**3-2D0*Zh11
     &      /SW*Zh12*Zh21*sb-2D0*Zh11/SW*Zh21*Zh22*cb-2D0
     &      *Zh11/SW*Zh31*Zh32*cb+2D0*Zh12/SW*Zh31**2*cb
     &      *sb**2+2D0*Zh21/SW*Zh31*Zh32*sb-4D0*Zh21/SW*Zh31
     &      *Zh32*sb**3-2D0*SW/CW**2*Zh11*Zh12*Zh21*sb-2D0
     &      *SW/CW**2*Zh11*Zh21*Zh22*cb-2D0*SW/CW**2*Zh11
     &      *Zh31*Zh32*cb+2D0*SW/CW**2*Zh12*Zh31**2*cb*sb
     &      **2+2D0*SW/CW**2*Zh21*Zh31*Zh32*sb-4D0*SW/CW
     &      **2*Zh21*Zh31*Zh32*sb**3+4D0*Zh11/SW*Zh31*Zh32
     &      *cb*sb**2+4D0*SW/CW**2*Zh11*Zh31*Zh32*cb*sb**2)
      AAABR(240) = EE*MW*(3D0*Zh11/SW*Zh12**2*cb-Zh11/SW
     &      *Zh22**2*cb-Zh11/SW*Zh32**2*cb-Zh12**2/SW*Zh21
     &      *sb+3D0*Zh21/SW*Zh22**2*sb+Zh21/SW*Zh32**2*sb
     &      -2D0*Zh21/SW*Zh32**2*sb**3+3D0*SW/CW**2*Zh11
     &      *Zh12**2*cb-SW/CW**2*Zh11*Zh22**2*cb-SW/CW**2
     &      *Zh11*Zh32**2*cb-SW/CW**2*Zh12**2*Zh21*sb+3D0
     &      *SW/CW**2*Zh21*Zh22**2*sb+SW/CW**2*Zh21*Zh32
     &      **2*sb-2D0*SW/CW**2*Zh21*Zh32**2*sb**3-2D0*Zh11
     &      /SW*Zh12*Zh22*sb+2D0*Zh11/SW*Zh32**2*cb*sb**2
     &      -2D0*Zh12/SW*Zh21*Zh22*cb-2D0*Zh12/SW*Zh31*Zh32
     &      *cb+2D0*Zh22/SW*Zh31*Zh32*sb-4D0*Zh22/SW*Zh31
     &      *Zh32*sb**3-2D0*SW/CW**2*Zh11*Zh12*Zh22*sb+2D0
     &      *SW/CW**2*Zh11*Zh32**2*cb*sb**2-2D0*SW/CW**2
     &      *Zh12*Zh21*Zh22*cb-2D0*SW/CW**2*Zh12*Zh31*Zh32
     &      *cb+2D0*SW/CW**2*Zh22*Zh31*Zh32*sb-4D0*SW/CW
     &      **2*Zh22*Zh31*Zh32*sb**3+4D0*Zh12/SW*Zh31*Zh32
     &      *cb*sb**2+4D0*SW/CW**2*Zh12*Zh31*Zh32*cb*sb**2)
      AAABR(241) = EE*MW*(3D0*Zh11/SW*Zh12*Zh13*cb-Zh11/SW
     &      *Zh12*Zh23*sb-Zh11/SW*Zh13*Zh22*sb-Zh11/SW*Zh22
     &      *Zh23*cb-Zh11/SW*Zh32*Zh33*cb-Zh12/SW*Zh13*Zh21
     &      *sb-Zh12/SW*Zh21*Zh23*cb-Zh12/SW*Zh31*Zh33*cb
     &      -Zh13/SW*Zh21*Zh22*cb-Zh13/SW*Zh31*Zh32*cb+3D0
     &      *Zh21/SW*Zh22*Zh23*sb+Zh21/SW*Zh32*Zh33*sb-2D0
     &      *Zh21/SW*Zh32*Zh33*sb**3+Zh22/SW*Zh31*Zh33*sb
     &      -2D0*Zh22/SW*Zh31*Zh33*sb**3+Zh23/SW*Zh31*Zh32
     &      *sb-2D0*Zh23/SW*Zh31*Zh32*sb**3+3D0*SW/CW**2
     &      *Zh11*Zh12*Zh13*cb-SW/CW**2*Zh11*Zh12*Zh23*sb
     &      -SW/CW**2*Zh11*Zh13*Zh22*sb-SW/CW**2*Zh11*Zh22
     &      *Zh23*cb-SW/CW**2*Zh11*Zh32*Zh33*cb-SW/CW**2
     &      *Zh12*Zh13*Zh21*sb-SW/CW**2*Zh12*Zh21*Zh23*cb
     &      -SW/CW**2*Zh12*Zh31*Zh33*cb-SW/CW**2*Zh13*Zh21
     &      *Zh22*cb-SW/CW**2*Zh13*Zh31*Zh32*cb+3D0*SW/CW
     &      **2*Zh21*Zh22*Zh23*sb+SW/CW**2*Zh21*Zh32*Zh33
     &      *sb-2D0*SW/CW**2*Zh21*Zh32*Zh33*sb**3+SW/CW**2
     &      *Zh22*Zh31*Zh33*sb-2D0*SW/CW**2*Zh22*Zh31*Zh33
     &      *sb**3+SW/CW**2*Zh23*Zh31*Zh32*sb-2D0*SW/CW**2
     &      *Zh23*Zh31*Zh32*sb**3+2D0*Zh11/SW*Zh32*Zh33*cb
     &      *sb**2+2D0*Zh12/SW*Zh31*Zh33*cb*sb**2+2D0*Zh13
     &      /SW*Zh31*Zh32*cb*sb**2+2D0*SW/CW**2*Zh11*Zh32
     &      *Zh33*cb*sb**2+2D0*SW/CW**2*Zh12*Zh31*Zh33*cb
     &      *sb**2+2D0*SW/CW**2*Zh13*Zh31*Zh32*cb*sb**2)
      AAABR(242) = EE*MW*(3D0*Zh11/SW*Zh13**2*cb-Zh11/SW
     &      *Zh23**2*cb-Zh11/SW*Zh33**2*cb-Zh13**2/SW*Zh21
     &      *sb+3D0*Zh21/SW*Zh23**2*sb+Zh21/SW*Zh33**2*sb
     &      -2D0*Zh21/SW*Zh33**2*sb**3+3D0*SW/CW**2*Zh11
     &      *Zh13**2*cb-SW/CW**2*Zh11*Zh23**2*cb-SW/CW**2
     &      *Zh11*Zh33**2*cb-SW/CW**2*Zh13**2*Zh21*sb+3D0
     &      *SW/CW**2*Zh21*Zh23**2*sb+SW/CW**2*Zh21*Zh33
     &      **2*sb-2D0*SW/CW**2*Zh21*Zh33**2*sb**3-2D0*Zh11
     &      /SW*Zh13*Zh23*sb+2D0*Zh11/SW*Zh33**2*cb*sb**2
     &      -2D0*Zh13/SW*Zh21*Zh23*cb-2D0*Zh13/SW*Zh31*Zh33
     &      *cb+2D0*Zh23/SW*Zh31*Zh33*sb-4D0*Zh23/SW*Zh31
     &      *Zh33*sb**3-2D0*SW/CW**2*Zh11*Zh13*Zh23*sb+2D0
     &      *SW/CW**2*Zh11*Zh33**2*cb*sb**2-2D0*SW/CW**2
     &      *Zh13*Zh21*Zh23*cb-2D0*SW/CW**2*Zh13*Zh31*Zh33
     &      *cb+2D0*SW/CW**2*Zh23*Zh31*Zh33*sb-4D0*SW/CW
     &      **2*Zh23*Zh31*Zh33*sb**3+4D0*Zh13/SW*Zh31*Zh33
     &      *cb*sb**2+4D0*SW/CW**2*Zh13*Zh31*Zh33*cb*sb**2)
      AAABR(243) = EE*MW*(3D0*Zh12/SW*Zh13**2*cb-Zh12/SW
     &      *Zh23**2*cb-Zh12/SW*Zh33**2*cb-Zh13**2/SW*Zh22
     &      *sb+3D0*Zh22/SW*Zh23**2*sb+Zh22/SW*Zh33**2*sb
     &      -2D0*Zh22/SW*Zh33**2*sb**3+3D0*SW/CW**2*Zh12
     &      *Zh13**2*cb-SW/CW**2*Zh12*Zh23**2*cb-SW/CW**2
     &      *Zh12*Zh33**2*cb-SW/CW**2*Zh13**2*Zh22*sb+3D0
     &      *SW/CW**2*Zh22*Zh23**2*sb+SW/CW**2*Zh22*Zh33
     &      **2*sb-2D0*SW/CW**2*Zh22*Zh33**2*sb**3-2D0*Zh12
     &      /SW*Zh13*Zh23*sb+2D0*Zh12/SW*Zh33**2*cb*sb**2
     &      -2D0*Zh13/SW*Zh22*Zh23*cb-2D0*Zh13/SW*Zh32*Zh33
     &      *cb+2D0*Zh23/SW*Zh32*Zh33*sb-4D0*Zh23/SW*Zh32
     &      *Zh33*sb**3-2D0*SW/CW**2*Zh12*Zh13*Zh23*sb+2D0
     &      *SW/CW**2*Zh12*Zh33**2*cb*sb**2-2D0*SW/CW**2
     &      *Zh13*Zh22*Zh23*cb-2D0*SW/CW**2*Zh13*Zh32*Zh33
     &      *cb+2D0*SW/CW**2*Zh23*Zh32*Zh33*sb-4D0*SW/CW
     &      **2*Zh23*Zh32*Zh33*sb**3+4D0*Zh13/SW*Zh32*Zh33
     &      *cb*sb**2+4D0*SW/CW**2*Zh13*Zh32*Zh33*cb*sb**2)
      AAABR(244) = EE*MW*(Zh13**3/SW*cb+Zh23**3/SW*sb+SW
     &      /CW**2*Zh13**3*cb+SW/CW**2*Zh23**3*sb-Zh13/SW
     &      *Zh23**2*cb-Zh13/SW*Zh33**2*cb-Zh13**2/SW*Zh23
     &      *sb+Zh23/SW*Zh33**2*sb-2D0*Zh23/SW*Zh33**2*sb
     &      **3-SW/CW**2*Zh13*Zh23**2*cb-SW/CW**2*Zh13*Zh33
     &      **2*cb-SW/CW**2*Zh13**2*Zh23*sb+SW/CW**2*Zh23
     &      *Zh33**2*sb-2D0*SW/CW**2*Zh23*Zh33**2*sb**3+2D0
     &      *Zh13/SW*Zh33**2*cb*sb**2+2D0*SW/CW**2*Zh13*Zh33
     &      **2*cb*sb**2)
      AAABR(245) = EE*MW*(Zh11**3/SW*cb+Zh21**3/SW*sb+SW
     &      /CW**2*Zh11**3*cb+SW/CW**2*Zh21**3*sb-Zh11/SW
     &      *Zh21**2*cb-Zh11/SW*Zh31**2*cb-Zh11**2/SW*Zh21
     &      *sb+Zh21/SW*Zh31**2*sb-2D0*Zh21/SW*Zh31**2*sb
     &      **3-SW/CW**2*Zh11*Zh21**2*cb-SW/CW**2*Zh11*Zh31
     &      **2*cb-SW/CW**2*Zh11**2*Zh21*sb+SW/CW**2*Zh21
     &      *Zh31**2*sb-2D0*SW/CW**2*Zh21*Zh31**2*sb**3+2D0
     &      *Zh11/SW*Zh31**2*cb*sb**2+2D0*SW/CW**2*Zh11*Zh31
     &      **2*cb*sb**2)
      AAABR(246) = EE*MW*(3D0*Zh12**2/SW*Zh13*cb-Zh12**2
     &      /SW*Zh23*sb-Zh13/SW*Zh22**2*cb-Zh13/SW*Zh32**2
     &      *cb+3D0*Zh22**2/SW*Zh23*sb+Zh23/SW*Zh32**2*sb
     &      -2D0*Zh23/SW*Zh32**2*sb**3+3D0*SW/CW**2*Zh12
     &      **2*Zh13*cb-SW/CW**2*Zh12**2*Zh23*sb-SW/CW**2
     &      *Zh13*Zh22**2*cb-SW/CW**2*Zh13*Zh32**2*cb+3D0
     &      *SW/CW**2*Zh22**2*Zh23*sb+SW/CW**2*Zh23*Zh32
     &      **2*sb-2D0*SW/CW**2*Zh23*Zh32**2*sb**3-2D0*Zh12
     &      /SW*Zh13*Zh22*sb-2D0*Zh12/SW*Zh22*Zh23*cb-2D0
     &      *Zh12/SW*Zh32*Zh33*cb+2D0*Zh13/SW*Zh32**2*cb
     &      *sb**2+2D0*Zh22/SW*Zh32*Zh33*sb-4D0*Zh22/SW*Zh32
     &      *Zh33*sb**3-2D0*SW/CW**2*Zh12*Zh13*Zh22*sb-2D0
     &      *SW/CW**2*Zh12*Zh22*Zh23*cb-2D0*SW/CW**2*Zh12
     &      *Zh32*Zh33*cb+2D0*SW/CW**2*Zh13*Zh32**2*cb*sb
     &      **2+2D0*SW/CW**2*Zh22*Zh32*Zh33*sb-4D0*SW/CW
     &      **2*Zh22*Zh32*Zh33*sb**3+4D0*Zh12/SW*Zh32*Zh33
     &      *cb*sb**2+4D0*SW/CW**2*Zh12*Zh32*Zh33*cb*sb**2)
      AAABR(247) = EE*MW*(3D0*Zh11**2/SW*Zh13*cb-Zh11**2
     &      /SW*Zh23*sb-Zh13/SW*Zh21**2*cb-Zh13/SW*Zh31**2
     &      *cb+3D0*Zh21**2/SW*Zh23*sb+Zh23/SW*Zh31**2*sb
     &      -2D0*Zh23/SW*Zh31**2*sb**3+3D0*SW/CW**2*Zh11
     &      **2*Zh13*cb-SW/CW**2*Zh11**2*Zh23*sb-SW/CW**2
     &      *Zh13*Zh21**2*cb-SW/CW**2*Zh13*Zh31**2*cb+3D0
     &      *SW/CW**2*Zh21**2*Zh23*sb+SW/CW**2*Zh23*Zh31
     &      **2*sb-2D0*SW/CW**2*Zh23*Zh31**2*sb**3-2D0*Zh11
     &      /SW*Zh13*Zh21*sb-2D0*Zh11/SW*Zh21*Zh23*cb-2D0
     &      *Zh11/SW*Zh31*Zh33*cb+2D0*Zh13/SW*Zh31**2*cb
     &      *sb**2+2D0*Zh21/SW*Zh31*Zh33*sb-4D0*Zh21/SW*Zh31
     &      *Zh33*sb**3-2D0*SW/CW**2*Zh11*Zh13*Zh21*sb-2D0
     &      *SW/CW**2*Zh11*Zh21*Zh23*cb-2D0*SW/CW**2*Zh11
     &      *Zh31*Zh33*cb+2D0*SW/CW**2*Zh13*Zh31**2*cb*sb
     &      **2+2D0*SW/CW**2*Zh21*Zh31*Zh33*sb-4D0*SW/CW
     &      **2*Zh21*Zh31*Zh33*sb**3+4D0*Zh11/SW*Zh31*Zh33
     &      *cb*sb**2+4D0*SW/CW**2*Zh11*Zh31*Zh33*cb*sb**2)
      AAABR(248) = EE*MW*(Zh12**3/SW*cb+Zh22**3/SW*sb+SW
     &      /CW**2*Zh12**3*cb+SW/CW**2*Zh22**3*sb-Zh12/SW
     &      *Zh22**2*cb-Zh12/SW*Zh32**2*cb-Zh12**2/SW*Zh22
     &      *sb+Zh22/SW*Zh32**2*sb-2D0*Zh22/SW*Zh32**2*sb
     &      **3-SW/CW**2*Zh12*Zh22**2*cb-SW/CW**2*Zh12*Zh32
     &      **2*cb-SW/CW**2*Zh12**2*Zh22*sb+SW/CW**2*Zh22
     &      *Zh32**2*sb-2D0*SW/CW**2*Zh22*Zh32**2*sb**3+2D0
     &      *Zh12/SW*Zh32**2*cb*sb**2+2D0*SW/CW**2*Zh12*Zh32
     &      **2*cb*sb**2)
      AAABR(249) = EE*(SW/CW-CW/SW)
      AAABR(250) = EE/SW*(Zh13*sb-Zh23*cb)
      AAABR(251) = EE/SW*Zh33
      AAABR(252) = EE/SW*(Zh12*sb-Zh22*cb)
      AAABR(253) = EE/SW*Zh32
      AAABR(254) = EE/SW*(Zh11*sb-Zh21*cb)
      AAABR(255) = EE/SW*Zh31
      AAABR(256) = EE*(SW/CW*Zl33i**2+SW/CW*Zl33r**2+2D0
     &      *SW/CW*Zl63i**2+2D0*SW/CW*Zl63r**2-CW/SW*Zl33i
     &      **2-CW/SW*Zl33r**2)
      AAABR(257) = EE*(SW/CW*Zl11**2+2D0*SW/CW*Zl41**2-CW
     &      /SW*Zl11**2)
      AAABR(258) = EE*(SW/CW*Zl22**2+2D0*SW/CW*Zl52**2-CW
     &      /SW*Zl22**2)
      AAABR(259) = EE*(SW/CW*Zl22*Zl25+2D0*SW/CW*Zl52*Zl55
     &      -CW/SW*Zl22*Zl25)
      AAABR(260) = EE*(SW/CW*Zl11*Zl14+2D0*SW/CW*Zl41*Zl44
     &      -CW/SW*Zl11*Zl14)
      AAABR(261) = EE*(SW/CW*Zl33i*Zl36i+SW/CW*Zl33r*Zl36r
     &      +2D0*SW/CW*Zl63i*Zl66i+2D0*SW/CW*Zl63r*Zl66r
     &      -CW/SW*Zl33i*Zl36i-CW/SW*Zl33r*Zl36r)
      AAABR(262) = EE*(SW/CW*Zl33i*Zl36r-SW/CW*Zl33r*Zl36i
     &      +2D0*SW/CW*Zl63i*Zl66r-2D0*SW/CW*Zl63r*Zl66i
     &      -CW/SW*Zl33i*Zl36r+CW/SW*Zl33r*Zl36i)
      AAABR(263) = EE/SW*Sqrt2*Zl33r
      AAABR(264) = EE/SW*Sqrt2*Zl33i
      AAABR(265) = EE/SW*Sqrt2*Zl11
      AAABR(266) = EE/SW*Sqrt2*Zl22
      AAABR(267) = EE*(SW/CW*Zl14**2+2D0*SW/CW*Zl44**2-CW
     &      /SW*Zl14**2)
      AAABR(268) = EE*(SW/CW*Zl36i**2+SW/CW*Zl36r**2+2D0
     &      *SW/CW*Zl66i**2+2D0*SW/CW*Zl66r**2-CW/SW*Zl36i
     &      **2-CW/SW*Zl36r**2)
      AAABR(269) = EE*(SW/CW*Zl25**2+2D0*SW/CW*Zl55**2-CW
     &      /SW*Zl25**2)
      AAABR(270) = EE/SW*Sqrt2*Zl36r
      AAABR(271) = EE/SW*Sqrt2*Zl36i
      AAABR(272) = EE/SW*Sqrt2*Zl25
      AAABR(273) = EE/SW*Sqrt2*Zl14
      AAABR(274) = EE*(SW/CW+CW/SW)
      AAABR(275) = EE*(SW/CW*Zd33i**2+SW/CW*Zd33r**2-2D0
     &      *SW/CW*Zd63i**2-2D0*SW/CW*Zd63r**2+3D0*CW/SW
     &      *Zd33i**2+3D0*CW/SW*Zd33r**2)
      AAABR(276) = EE*(SW/CW*Zd22**2-2D0*SW/CW*Zd52**2+3D0
     &      *CW/SW*Zd22**2)
      AAABR(277) = EE*(SW/CW*Zd11**2-2D0*SW/CW*Zd41**2+3D0
     &      *CW/SW*Zd11**2)
      AAABR(278) = EE*(SW/CW*Zd22*Zd25-2D0*SW/CW*Zd52*Zd55
     &      +3D0*CW/SW*Zd22*Zd25)
      AAABR(279) = EE*(SW/CW*Zd33i*Zd36i+SW/CW*Zd33r*Zd36r
     &      -2D0*SW/CW*Zd63i*Zd66i-2D0*SW/CW*Zd63r*Zd66r
     &      +3D0*CW/SW*Zd33i*Zd36i+3D0*CW/SW*Zd33r*Zd36r)
      AAABR(280) = EE*(SW/CW*Zd33i*Zd36r-SW/CW*Zd33r*Zd36i
     &      -2D0*SW/CW*Zd63i*Zd66r+2D0*SW/CW*Zd63r*Zd66i
     &      +3D0*CW/SW*Zd33i*Zd36r-3D0*CW/SW*Zd33r*Zd36i)
      AAABR(281) = EE*(SW/CW*Zd11*Zd14-2D0*SW/CW*Zd41*Zd44
     &      +3D0*CW/SW*Zd11*Zd14)
      AAABR(282) = EE/SW*Sqrt2*Zd22*Zu22
      AAABR(283) = EE/SW*Sqrt2*(Zd33i*Zu33i+Zd33r*Zu33r)
      AAABR(284) = EE/SW*Sqrt2*(Zd33i*Zu33r-Zd33r*Zu33i)
      AAABR(285) = EE/SW*Sqrt2*Zd11*Zu11
      AAABR(286) = EE/SW*Sqrt2*Zd22*Zu25
      AAABR(287) = EE/SW*Sqrt2*Zd11*Zu14
      AAABR(288) = EE/SW*Sqrt2*(Zd33i*Zu36i+Zd33r*Zu36r)
      AAABR(289) = EE/SW*Sqrt2*(Zd33i*Zu36r-Zd33r*Zu36i)
      AAABR(290) = EE*(SW/CW*Zd25**2-2D0*SW/CW*Zd55**2+3D0
     &      *CW/SW*Zd25**2)
      AAABR(291) = EE*(SW/CW*Zd36i**2+SW/CW*Zd36r**2-2D0
     &      *SW/CW*Zd66i**2-2D0*SW/CW*Zd66r**2+3D0*CW/SW
     &      *Zd36i**2+3D0*CW/SW*Zd36r**2)
      AAABR(292) = EE*(SW/CW*Zd14**2-2D0*SW/CW*Zd44**2+3D0
     &      *CW/SW*Zd14**2)
      AAABR(293) = EE/SW*Sqrt2*Zd25*Zu22
      AAABR(294) = EE/SW*Sqrt2*Zd14*Zu11
      AAABR(295) = EE/SW*Sqrt2*(Zd36i*Zu33i+Zd36r*Zu33r)
      AAABR(296) = EE/SW*Sqrt2*(Zd36i*Zu33r-Zd36r*Zu33i)
      AAABR(297) = EE/SW*Sqrt2*Zd25*Zu25
      AAABR(298) = EE/SW*Sqrt2*(Zd36i*Zu36i+Zd36r*Zu36r)
      AAABR(299) = EE/SW*Sqrt2*(Zd36i*Zu36r-Zd36r*Zu36i)
      AAABR(300) = EE/SW*Sqrt2*Zd14*Zu14
      AAABR(301) = EE*(SW/CW*Zu11**2+4D0*SW/CW*Zu41**2-3D0
     &      *CW/SW*Zu11**2)
      AAABR(302) = EE*(SW/CW*Zu33i**2+SW/CW*Zu33r**2+4D0
     &      *SW/CW*Zu63i**2+4D0*SW/CW*Zu63r**2-3D0*CW/SW
     &      *Zu33i**2-3D0*CW/SW*Zu33r**2)
      AAABR(303) = EE*(SW/CW*Zu22**2+4D0*SW/CW*Zu52**2-3D0
     &      *CW/SW*Zu22**2)
      AAABR(304) = EE*(SW/CW*Zu11*Zu14+4D0*SW/CW*Zu41*Zu44
     &      -3D0*CW/SW*Zu11*Zu14)
      AAABR(305) = EE*(SW/CW*Zu22*Zu25+4D0*SW/CW*Zu52*Zu55
     &      -3D0*CW/SW*Zu22*Zu25)
      AAABR(306) = EE*(SW/CW*Zu33i*Zu36i+SW/CW*Zu33r*Zu36r
     &      +4D0*SW/CW*Zu63i*Zu66i+4D0*SW/CW*Zu63r*Zu66r
     &      -3D0*CW/SW*Zu33i*Zu36i-3D0*CW/SW*Zu33r*Zu36r)
      AAABR(307) = EE*(SW/CW*Zu33i*Zu36r-SW/CW*Zu33r*Zu36i
     &      +4D0*SW/CW*Zu63i*Zu66r-4D0*SW/CW*Zu63r*Zu66i
     &      -3D0*CW/SW*Zu33i*Zu36r+3D0*CW/SW*Zu33r*Zu36i)
      AAABR(308) = EE*(SW/CW*Zu36i**2+SW/CW*Zu36r**2+4D0
     &      *SW/CW*Zu66i**2+4D0*SW/CW*Zu66r**2-3D0*CW/SW
     &      *Zu36i**2-3D0*CW/SW*Zu36r**2)
      AAABR(309) = EE*(SW/CW*Zu14**2+4D0*SW/CW*Zu44**2-3D0
     &      *CW/SW*Zu14**2)
      AAABR(310) = EE*(SW/CW*Zu25**2+4D0*SW/CW*Zu55**2-3D0
     &      *CW/SW*Zu25**2)
      AAABR(311) = EE/SW
      AAABR(312) = EE/SW*(Zh11*cb+Zh21*sb)
      AAABR(313) = EE/SW*(Zh12*cb+Zh22*sb)
      AAABR(314) = EE/SW*(Zh13*cb+Zh23*sb)
      AAABR(315) = EE*(SW/CW*Zh12*cb+SW/CW*Zh22*sb+CW/SW
     &      *Zh12*cb+CW/SW*Zh22*sb)
      AAABR(316) = EE*(SW/CW*Zh13*cb+SW/CW*Zh23*sb+CW/SW
     &      *Zh13*cb+CW/SW*Zh23*sb)
      AAABR(317) = EE*(SW/CW*Zh11*cb+SW/CW*Zh21*sb+CW/SW
     &      *Zh11*cb+CW/SW*Zh21*sb)
      AAABR(318) = EE*(SW/CW*Zh12*Zh33*sb-SW/CW*Zh13*Zh32
     &      *sb-SW/CW*Zh22*Zh33*cb+SW/CW*Zh23*Zh32*cb+CW
     &      /SW*Zh12*Zh33*sb-CW/SW*Zh13*Zh32*sb-CW/SW*Zh22
     &      *Zh33*cb+CW/SW*Zh23*Zh32*cb)
      AAABR(319) = EE*(SW/CW*Zh11*Zh33*sb-SW/CW*Zh13*Zh31
     &      *sb-SW/CW*Zh21*Zh33*cb+SW/CW*Zh23*Zh31*cb+CW
     &      /SW*Zh11*Zh33*sb-CW/SW*Zh13*Zh31*sb-CW/SW*Zh21
     &      *Zh33*cb+CW/SW*Zh23*Zh31*cb)
      AAABR(320) = EE*(SW/CW*Zh11*Zh32*sb-SW/CW*Zh12*Zh31
     &      *sb-SW/CW*Zh21*Zh32*cb+SW/CW*Zh22*Zh31*cb+CW
     &      /SW*Zh11*Zh32*sb-CW/SW*Zh12*Zh31*sb-CW/SW*Zh21
     &      *Zh32*cb+CW/SW*Zh22*Zh31*cb)
      AAABR(321) = EE/CW*MW*SW
      AAABR(322) = EE*MW*(2D0*SW*Zh11*cb+2D0*SW*Zh21*sb+SW
     &      **3/CW**2*Zh11*cb+SW**3/CW**2*Zh21*sb+CW**2/SW
     &      *Zh11*cb+CW**2/SW*Zh21*sb)
      AAABR(323) = EE*MW*(2D0*SW*Zh12*cb+2D0*SW*Zh22*sb+SW
     &      **3/CW**2*Zh12*cb+SW**3/CW**2*Zh22*sb+CW**2/SW
     &      *Zh12*cb+CW**2/SW*Zh22*sb)
      AAABR(324) = EE*MW*(2D0*SW*Zh13*cb+2D0*SW*Zh23*sb+SW
     &      **3/CW**2*Zh13*cb+SW**3/CW**2*Zh23*sb+CW**2/SW
     &      *Zh13*cb+CW**2/SW*Zh23*sb)
      AAABR(325) = EE/SW*Sqrt2*(Zm12i*Zp22r*sb-Zm12r*Zp22i
     &      *sb-Zm22i*Zp12r*cb+Zm22r*Zp12i*cb)
      AAABR(326) = EE/SW*Sqrt2*(Zm12i*Zp22i*sb+Zm12r*Zp22r
     &      *sb-Zm22i*Zp12i*cb-Zm22r*Zp12r*cb)
      AAABR(327) = EE/SW*Sqrt2*(Zm11i*Zp22i*sb+Zm11r*Zp22r
     &      *sb-Zm21i*Zp12i*cb-Zm21r*Zp12r*cb)
      AAABR(328) = EE/SW*Sqrt2*(Zm11i*Zp22r*sb-Zm11r*Zp22i
     &      *sb-Zm21i*Zp12r*cb+Zm21r*Zp12i*cb)
      AAABR(329) = EE/SW*Sqrt2*(Zm12i*Zp21i*sb+Zm12r*Zp21r
     &      *sb-Zm22i*Zp11i*cb-Zm22r*Zp11r*cb)
      AAABR(330) = EE/SW*Sqrt2*(Zm12i*Zp21r*sb-Zm12r*Zp21i
     &      *sb-Zm22i*Zp11r*cb+Zm22r*Zp11i*cb)
      AAABR(331) = EE/SW*Sqrt2*(Zm11i*Zp21r*sb-Zm11r*Zp21i
     &      *sb-Zm21i*Zp11r*cb+Zm21r*Zp11i*cb)
      AAABR(332) = EE/SW*Sqrt2*(Zm11i*Zp21i*sb+Zm11r*Zp21r
     &      *sb-Zm21i*Zp11i*cb-Zm21r*Zp11r*cb)
      AAABR(333) = EE/SW*Sqrt2*(Zh11*Zm22i*Zp12i+Zh11*Zm22r
     &      *Zp12r+Zh21*Zm12i*Zp22i+Zh21*Zm12r*Zp22r-Zh31
     &      *Zm12i*Zp22r*cb+Zh31*Zm12r*Zp22i*cb-Zh31*Zm22i
     &      *Zp12r*sb+Zh31*Zm22r*Zp12i*sb)
      AAABR(334) = EE/SW*Sqrt2*(Zh11*Zm22i*Zp12r-Zh11*Zm22r
     &      *Zp12i+Zh21*Zm12i*Zp22r-Zh21*Zm12r*Zp22i+Zh31
     &      *Zm12i*Zp22i*cb+Zh31*Zm12r*Zp22r*cb+Zh31*Zm22i
     &      *Zp12i*sb+Zh31*Zm22r*Zp12r*sb)
      AAABR(335) = EE/SW*Sqrt2*(Zh12*Zm22i*Zp12i+Zh12*Zm22r
     &      *Zp12r+Zh22*Zm12i*Zp22i+Zh22*Zm12r*Zp22r-Zh32
     &      *Zm12i*Zp22r*cb+Zh32*Zm12r*Zp22i*cb-Zh32*Zm22i
     &      *Zp12r*sb+Zh32*Zm22r*Zp12i*sb)
      AAABR(336) = EE/SW*Sqrt2*(Zh12*Zm22i*Zp12r-Zh12*Zm22r
     &      *Zp12i+Zh22*Zm12i*Zp22r-Zh22*Zm12r*Zp22i+Zh32
     &      *Zm12i*Zp22i*cb+Zh32*Zm12r*Zp22r*cb+Zh32*Zm22i
     &      *Zp12i*sb+Zh32*Zm22r*Zp12r*sb)
      AAABR(337) = EE/SW*Sqrt2*(Zh11*Zm21i*Zp11i+Zh11*Zm21r
     &      *Zp11r+Zh21*Zm11i*Zp21i+Zh21*Zm11r*Zp21r-Zh31
     &      *Zm11i*Zp21r*cb+Zh31*Zm11r*Zp21i*cb-Zh31*Zm21i
     &      *Zp11r*sb+Zh31*Zm21r*Zp11i*sb)
      AAABR(338) = EE/SW*Sqrt2*(Zh11*Zm21i*Zp11r-Zh11*Zm21r
     &      *Zp11i+Zh21*Zm11i*Zp21r-Zh21*Zm11r*Zp21i+Zh31
     &      *Zm11i*Zp21i*cb+Zh31*Zm11r*Zp21r*cb+Zh31*Zm21i
     &      *Zp11i*sb+Zh31*Zm21r*Zp11r*sb)
      AAABR(339) = EE/SW*Sqrt2*(Zh12*Zm22i*Zp11i+Zh12*Zm22r
     &      *Zp11r+Zh22*Zm12i*Zp21i+Zh22*Zm12r*Zp21r-Zh32
     &      *Zm12i*Zp21r*cb+Zh32*Zm12r*Zp21i*cb-Zh32*Zm22i
     &      *Zp11r*sb+Zh32*Zm22r*Zp11i*sb)
      AAABR(340) = EE/SW*Sqrt2*(Zh12*Zm22i*Zp11r-Zh12*Zm22r
     &      *Zp11i+Zh22*Zm12i*Zp21r-Zh22*Zm12r*Zp21i+Zh32
     &      *Zm12i*Zp21i*cb+Zh32*Zm12r*Zp21r*cb+Zh32*Zm22i
     &      *Zp11i*sb+Zh32*Zm22r*Zp11r*sb)
      AAABR(341) = EE/SW*Sqrt2*(Zh12*Zm21i*Zp12i+Zh12*Zm21r
     &      *Zp12r+Zh22*Zm11i*Zp22i+Zh22*Zm11r*Zp22r-Zh32
     &      *Zm11i*Zp22r*cb+Zh32*Zm11r*Zp22i*cb-Zh32*Zm21i
     &      *Zp12r*sb+Zh32*Zm21r*Zp12i*sb)
      AAABR(342) = EE/SW*Sqrt2*(Zh12*Zm21i*Zp12r-Zh12*Zm21r
     &      *Zp12i+Zh22*Zm11i*Zp22r-Zh22*Zm11r*Zp22i+Zh32
     &      *Zm11i*Zp22i*cb+Zh32*Zm11r*Zp22r*cb+Zh32*Zm21i
     &      *Zp12i*sb+Zh32*Zm21r*Zp12r*sb)
      AAABR(343) = EE/SW*Sqrt2*(Zh13*Zm22i*Zp12i+Zh13*Zm22r
     &      *Zp12r+Zh23*Zm12i*Zp22i+Zh23*Zm12r*Zp22r-Zh33
     &      *Zm12i*Zp22r*cb+Zh33*Zm12r*Zp22i*cb-Zh33*Zm22i
     &      *Zp12r*sb+Zh33*Zm22r*Zp12i*sb)
      AAABR(344) = EE/SW*Sqrt2*(Zh13*Zm22i*Zp12r-Zh13*Zm22r
     &      *Zp12i+Zh23*Zm12i*Zp22r-Zh23*Zm12r*Zp22i+Zh33
     &      *Zm12i*Zp22i*cb+Zh33*Zm12r*Zp22r*cb+Zh33*Zm22i
     &      *Zp12i*sb+Zh33*Zm22r*Zp12r*sb)
      AAABR(345) = EE/SW*Sqrt2*(Zh12*Zm21i*Zp11i+Zh12*Zm21r
     &      *Zp11r+Zh22*Zm11i*Zp21i+Zh22*Zm11r*Zp21r-Zh32
     &      *Zm11i*Zp21r*cb+Zh32*Zm11r*Zp21i*cb-Zh32*Zm21i
     &      *Zp11r*sb+Zh32*Zm21r*Zp11i*sb)
      AAABR(346) = EE/SW*Sqrt2*(Zh12*Zm21i*Zp11r-Zh12*Zm21r
     &      *Zp11i+Zh22*Zm11i*Zp21r-Zh22*Zm11r*Zp21i+Zh32
     &      *Zm11i*Zp21i*cb+Zh32*Zm11r*Zp21r*cb+Zh32*Zm21i
     &      *Zp11i*sb+Zh32*Zm21r*Zp11r*sb)
      AAABR(347) = EE/SW*Sqrt2*(Zh11*Zm22i*Zp11i+Zh11*Zm22r
     &      *Zp11r+Zh21*Zm12i*Zp21i+Zh21*Zm12r*Zp21r-Zh31
     &      *Zm12i*Zp21r*cb+Zh31*Zm12r*Zp21i*cb-Zh31*Zm22i
     &      *Zp11r*sb+Zh31*Zm22r*Zp11i*sb)
      AAABR(348) = EE/SW*Sqrt2*(Zh11*Zm22i*Zp11r-Zh11*Zm22r
     &      *Zp11i+Zh21*Zm12i*Zp21r-Zh21*Zm12r*Zp21i+Zh31
     &      *Zm12i*Zp21i*cb+Zh31*Zm12r*Zp21r*cb+Zh31*Zm22i
     &      *Zp11i*sb+Zh31*Zm22r*Zp11r*sb)
      AAABR(349) = EE/SW*Sqrt2*(Zh11*Zm21i*Zp12i+Zh11*Zm21r
     &      *Zp12r+Zh21*Zm11i*Zp22i+Zh21*Zm11r*Zp22r-Zh31
     &      *Zm11i*Zp22r*cb+Zh31*Zm11r*Zp22i*cb-Zh31*Zm21i
     &      *Zp12r*sb+Zh31*Zm21r*Zp12i*sb)
      AAABR(350) = EE/SW*Sqrt2*(Zh11*Zm21i*Zp12r-Zh11*Zm21r
     &      *Zp12i+Zh21*Zm11i*Zp22r-Zh21*Zm11r*Zp22i+Zh31
     &      *Zm11i*Zp22i*cb+Zh31*Zm11r*Zp22r*cb+Zh31*Zm21i
     &      *Zp12i*sb+Zh31*Zm21r*Zp12r*sb)
      AAABR(351) = EE/SW*Sqrt2*(Zh13*Zm22i*Zp11i+Zh13*Zm22r
     &      *Zp11r+Zh23*Zm12i*Zp21i+Zh23*Zm12r*Zp21r-Zh33
     &      *Zm12i*Zp21r*cb+Zh33*Zm12r*Zp21i*cb-Zh33*Zm22i
     &      *Zp11r*sb+Zh33*Zm22r*Zp11i*sb)
      AAABR(352) = EE/SW*Sqrt2*(Zh13*Zm22i*Zp11r-Zh13*Zm22r
     &      *Zp11i+Zh23*Zm12i*Zp21r-Zh23*Zm12r*Zp21i+Zh33
     &      *Zm12i*Zp21i*cb+Zh33*Zm12r*Zp21r*cb+Zh33*Zm22i
     &      *Zp11i*sb+Zh33*Zm22r*Zp11r*sb)
      AAABR(353) = EE/SW*Sqrt2*(Zh13*Zm21i*Zp12i+Zh13*Zm21r
     &      *Zp12r+Zh23*Zm11i*Zp22i+Zh23*Zm11r*Zp22r-Zh33
     &      *Zm11i*Zp22r*cb+Zh33*Zm11r*Zp22i*cb-Zh33*Zm21i
     &      *Zp12r*sb+Zh33*Zm21r*Zp12i*sb)
      AAABR(354) = EE/SW*Sqrt2*(Zh13*Zm21i*Zp12r-Zh13*Zm21r
     &      *Zp12i+Zh23*Zm11i*Zp22r-Zh23*Zm11r*Zp22i+Zh33
     &      *Zm11i*Zp22i*cb+Zh33*Zm11r*Zp22r*cb+Zh33*Zm21i
     &      *Zp12i*sb+Zh33*Zm21r*Zp12r*sb)
      AAABR(355) = EE/SW*Sqrt2*(Zh13*Zm21i*Zp11i+Zh13*Zm21r
     &      *Zp11r+Zh23*Zm11i*Zp21i+Zh23*Zm11r*Zp21r-Zh33
     &      *Zm11i*Zp21r*cb+Zh33*Zm11r*Zp21i*cb-Zh33*Zm21i
     &      *Zp11r*sb+Zh33*Zm21r*Zp11i*sb)
      AAABR(356) = EE/SW*Sqrt2*(Zh13*Zm21i*Zp11r-Zh13*Zm21r
     &      *Zp11i+Zh23*Zm11i*Zp21r-Zh23*Zm11r*Zp21i+Zh33
     &      *Zm11i*Zp21i*cb+Zh33*Zm11r*Zp21r*cb+Zh33*Zm21i
     &      *Zp11i*sb+Zh33*Zm21r*Zp11r*sb)
      AAABR(357) = EE/SW*(2D0*Zl22*Zm12r-Mm/MW*Sqrt2*Zl52
     &      *Zm22r/cb)
      AAABR(358) = EE/SW*(2D0*Zl22*Zm12i-Mm/MW*Sqrt2*Zl52
     &      *Zm22i/cb)
      AAABR(359) = EE/SW*(2D0*Zl11*Zm12r-Me/MW*Sqrt2*Zl41
     &      *Zm22r/cb)
      AAABR(360) = EE/SW*(2D0*Zl11*Zm12i-Me/MW*Sqrt2*Zl41
     &      *Zm22i/cb)
      AAABR(361) = EE/SW*(2D0*Zl11*Zm11r-Me/MW*Sqrt2*Zl41
     &      *Zm21r/cb)
      AAABR(362) = EE/SW*(2D0*Zl11*Zm11i-Me/MW*Sqrt2*Zl41
     &      *Zm21i/cb)
      AAABR(363) = EE/SW*(2D0*Zl22*Zm11r-Mm/MW*Sqrt2*Zl52
     &      *Zm21r/cb)
      AAABR(364) = EE/SW*(2D0*Zl22*Zm11i-Mm/MW*Sqrt2*Zl52
     &      *Zm21i/cb)
      AAABR(365) = EE/SW*(2D0*Zl33i*Zm11i-2D0*Zl33r*Zm11r
     &      -Ml/MW*Sqrt2*Zl63i*Zm21i/cb+Ml/MW*Sqrt2*Zl63r
     &      *Zm21r/cb)
      AAABR(366) = EE/SW*(2D0*Zl33i*Zm11r+2D0*Zl33r*Zm11i
     &      -Ml/MW*Sqrt2*Zl63i*Zm21r/cb-Ml/MW*Sqrt2*Zl63r
     &      *Zm21i/cb)
      AAABR(367) = EE/SW*(2D0*Zl33i*Zm12i-2D0*Zl33r*Zm12r
     &      -Ml/MW*Sqrt2*Zl63i*Zm22i/cb+Ml/MW*Sqrt2*Zl63r
     &      *Zm22r/cb)
      AAABR(368) = EE/SW*(2D0*Zl33i*Zm12r+2D0*Zl33r*Zm12i
     &      -Ml/MW*Sqrt2*Zl63i*Zm22r/cb-Ml/MW*Sqrt2*Zl63r
     &      *Zm22i/cb)
      AAABR(369) = EE/SW*(2D0*Zl25*Zm12r-Mm/MW*Sqrt2*Zl55
     &      *Zm22r/cb)
      AAABR(370) = EE/SW*(2D0*Zl25*Zm12i-Mm/MW*Sqrt2*Zl55
     &      *Zm22i/cb)
      AAABR(371) = EE/SW*(2D0*Zl25*Zm11r-Mm/MW*Sqrt2*Zl55
     &      *Zm21r/cb)
      AAABR(372) = EE/SW*(2D0*Zl25*Zm11i-Mm/MW*Sqrt2*Zl55
     &      *Zm21i/cb)
      AAABR(373) = EE/SW*(2D0*Zl14*Zm12r-Me/MW*Sqrt2*Zl44
     &      *Zm22r/cb)
      AAABR(374) = EE/SW*(2D0*Zl14*Zm12i-Me/MW*Sqrt2*Zl44
     &      *Zm22i/cb)
      AAABR(375) = EE/SW*(2D0*Zl36i*Zm11i-2D0*Zl36r*Zm11r
     &      -Ml/MW*Sqrt2*Zl66i*Zm21i/cb+Ml/MW*Sqrt2*Zl66r
     &      *Zm21r/cb)
      AAABR(376) = EE/SW*(2D0*Zl36i*Zm11r+2D0*Zl36r*Zm11i
     &      -Ml/MW*Sqrt2*Zl66i*Zm21r/cb-Ml/MW*Sqrt2*Zl66r
     &      *Zm21i/cb)
      AAABR(377) = EE/SW*(2D0*Zl14*Zm11r-Me/MW*Sqrt2*Zl44
     &      *Zm21r/cb)
      AAABR(378) = EE/SW*(2D0*Zl14*Zm11i-Me/MW*Sqrt2*Zl44
     &      *Zm21i/cb)
      AAABR(379) = EE/SW*(2D0*Zl36i*Zm12i-2D0*Zl36r*Zm12r
     &      -Ml/MW*Sqrt2*Zl66i*Zm22i/cb+Ml/MW*Sqrt2*Zl66r
     &      *Zm22r/cb)
      AAABR(380) = EE/SW*(2D0*Zl36i*Zm12r+2D0*Zl36r*Zm12i
     &      -Ml/MW*Sqrt2*Zl66i*Zm22r/cb-Ml/MW*Sqrt2*Zl66r
     &      *Zm22i/cb)
      AAABR(381) = EE*sb*(2D0*Zm11i/SW*Zn32i+2D0*Zm11r/SW
     &      *Zn32r-Sqrt2/CW*Zm21i*Zn12i-Sqrt2/CW*Zm21r*Zn12r
     &      -Sqrt2/SW*Zm21i*Zn22i-Sqrt2/SW*Zm21r*Zn22r)
      AAABR(382) = EE*sb*(2D0*Zm11i/SW*Zn32r-2D0*Zm11r/SW
     &      *Zn32i-Sqrt2/CW*Zm21i*Zn12r+Sqrt2/CW*Zm21r*Zn12i
     &      -Sqrt2/SW*Zm21i*Zn22r+Sqrt2/SW*Zm21r*Zn22i)
      AAABR(383) = EE*cb*(2D0*Zn42i/SW*Zp11i-2D0*Zn42r/SW
     &      *Zp11r+Sqrt2/CW*Zn12i*Zp21i-Sqrt2/CW*Zn12r*Zp21r
     &      +Sqrt2/SW*Zn22i*Zp21i-Sqrt2/SW*Zn22r*Zp21r)
      AAABR(384) = EE*cb*(2D0*Zn42i/SW*Zp11r+2D0*Zn42r/SW
     &      *Zp11i+Sqrt2/CW*Zn12i*Zp21r+Sqrt2/CW*Zn12r*Zp21i
     &      +Sqrt2/SW*Zn22i*Zp21r+Sqrt2/SW*Zn22r*Zp21i)
      AAABR(385) = EE*sb*(2D0*Zm11i/SW*Zn34i+2D0*Zm11r/SW
     &      *Zn34r-Sqrt2/CW*Zm21i*Zn14i-Sqrt2/CW*Zm21r*Zn14r
     &      -Sqrt2/SW*Zm21i*Zn24i-Sqrt2/SW*Zm21r*Zn24r)
      AAABR(386) = EE*sb*(2D0*Zm11i/SW*Zn34r-2D0*Zm11r/SW
     &      *Zn34i-Sqrt2/CW*Zm21i*Zn14r+Sqrt2/CW*Zm21r*Zn14i
     &      -Sqrt2/SW*Zm21i*Zn24r+Sqrt2/SW*Zm21r*Zn24i)
      AAABR(387) = EE*cb*(2D0*Zn44i/SW*Zp11i-2D0*Zn44r/SW
     &      *Zp11r+Sqrt2/CW*Zn14i*Zp21i-Sqrt2/CW*Zn14r*Zp21r
     &      +Sqrt2/SW*Zn24i*Zp21i-Sqrt2/SW*Zn24r*Zp21r)
      AAABR(388) = EE*cb*(2D0*Zn44i/SW*Zp11r+2D0*Zn44r/SW
     &      *Zp11i+Sqrt2/CW*Zn14i*Zp21r+Sqrt2/CW*Zn14r*Zp21i
     &      +Sqrt2/SW*Zn24i*Zp21r+Sqrt2/SW*Zn24r*Zp21i)
      AAABR(389) = EE*sb*(2D0*Zm11i/SW*Zn33i+2D0*Zm11r/SW
     &      *Zn33r-Sqrt2/CW*Zm21i*Zn13i-Sqrt2/CW*Zm21r*Zn13r
     &      -Sqrt2/SW*Zm21i*Zn23i-Sqrt2/SW*Zm21r*Zn23r)
      AAABR(390) = EE*sb*(2D0*Zm11i/SW*Zn33r-2D0*Zm11r/SW
     &      *Zn33i-Sqrt2/CW*Zm21i*Zn13r+Sqrt2/CW*Zm21r*Zn13i
     &      -Sqrt2/SW*Zm21i*Zn23r+Sqrt2/SW*Zm21r*Zn23i)
      AAABR(391) = EE*cb*(2D0*Zn43i/SW*Zp11i-2D0*Zn43r/SW
     &      *Zp11r+Sqrt2/CW*Zn13i*Zp21i-Sqrt2/CW*Zn13r*Zp21r
     &      +Sqrt2/SW*Zn23i*Zp21i-Sqrt2/SW*Zn23r*Zp21r)
      AAABR(392) = EE*cb*(2D0*Zn43i/SW*Zp11r+2D0*Zn43r/SW
     &      *Zp11i+Sqrt2/CW*Zn13i*Zp21r+Sqrt2/CW*Zn13r*Zp21i
     &      +Sqrt2/SW*Zn23i*Zp21r+Sqrt2/SW*Zn23r*Zp21i)
      AAABR(393) = EE*sb*(2D0*Zm11i/SW*Zn31i+2D0*Zm11r/SW
     &      *Zn31r-Sqrt2/CW*Zm21i*Zn11i-Sqrt2/CW*Zm21r*Zn11r
     &      -Sqrt2/SW*Zm21i*Zn21i-Sqrt2/SW*Zm21r*Zn21r)
      AAABR(394) = EE*sb*(2D0*Zm11i/SW*Zn31r-2D0*Zm11r/SW
     &      *Zn31i-Sqrt2/CW*Zm21i*Zn11r+Sqrt2/CW*Zm21r*Zn11i
     &      -Sqrt2/SW*Zm21i*Zn21r+Sqrt2/SW*Zm21r*Zn21i)
      AAABR(395) = EE*cb*(2D0*Zn41i/SW*Zp11i-2D0*Zn41r/SW
     &      *Zp11r+Sqrt2/CW*Zn11i*Zp21i-Sqrt2/CW*Zn11r*Zp21r
     &      +Sqrt2/SW*Zn21i*Zp21i-Sqrt2/SW*Zn21r*Zp21r)
      AAABR(396) = EE*cb*(2D0*Zn41i/SW*Zp11r+2D0*Zn41r/SW
     &      *Zp11i+Sqrt2/CW*Zn11i*Zp21r+Sqrt2/CW*Zn11r*Zp21i
     &      +Sqrt2/SW*Zn21i*Zp21r+Sqrt2/SW*Zn21r*Zp21i)
      AAABR(397) = EE*sb*(2D0*Zm12i/SW*Zn33i+2D0*Zm12r/SW
     &      *Zn33r-Sqrt2/CW*Zm22i*Zn13i-Sqrt2/CW*Zm22r*Zn13r
     &      -Sqrt2/SW*Zm22i*Zn23i-Sqrt2/SW*Zm22r*Zn23r)
      AAABR(398) = EE*sb*(2D0*Zm12i/SW*Zn33r-2D0*Zm12r/SW
     &      *Zn33i-Sqrt2/CW*Zm22i*Zn13r+Sqrt2/CW*Zm22r*Zn13i
     &      -Sqrt2/SW*Zm22i*Zn23r+Sqrt2/SW*Zm22r*Zn23i)
      AAABR(399) = EE*cb*(2D0*Zn43i/SW*Zp12i-2D0*Zn43r/SW
     &      *Zp12r+Sqrt2/CW*Zn13i*Zp22i-Sqrt2/CW*Zn13r*Zp22r
     &      +Sqrt2/SW*Zn23i*Zp22i-Sqrt2/SW*Zn23r*Zp22r)
      AAABR(400) = EE*cb*(2D0*Zn43i/SW*Zp12r+2D0*Zn43r/SW
     &      *Zp12i+Sqrt2/CW*Zn13i*Zp22r+Sqrt2/CW*Zn13r*Zp22i
     &      +Sqrt2/SW*Zn23i*Zp22r+Sqrt2/SW*Zn23r*Zp22i)
      AAABR(401) = EE*sb*(2D0*Zm12i/SW*Zn32i+2D0*Zm12r/SW
     &      *Zn32r-Sqrt2/CW*Zm22i*Zn12i-Sqrt2/CW*Zm22r*Zn12r
     &      -Sqrt2/SW*Zm22i*Zn22i-Sqrt2/SW*Zm22r*Zn22r)
      AAABR(402) = EE*sb*(2D0*Zm12i/SW*Zn32r-2D0*Zm12r/SW
     &      *Zn32i-Sqrt2/CW*Zm22i*Zn12r+Sqrt2/CW*Zm22r*Zn12i
     &      -Sqrt2/SW*Zm22i*Zn22r+Sqrt2/SW*Zm22r*Zn22i)
      AAABR(403) = EE*cb*(2D0*Zn42i/SW*Zp12i-2D0*Zn42r/SW
     &      *Zp12r+Sqrt2/CW*Zn12i*Zp22i-Sqrt2/CW*Zn12r*Zp22r
     &      +Sqrt2/SW*Zn22i*Zp22i-Sqrt2/SW*Zn22r*Zp22r)
      AAABR(404) = EE*cb*(2D0*Zn42i/SW*Zp12r+2D0*Zn42r/SW
     &      *Zp12i+Sqrt2/CW*Zn12i*Zp22r+Sqrt2/CW*Zn12r*Zp22i
     &      +Sqrt2/SW*Zn22i*Zp22r+Sqrt2/SW*Zn22r*Zp22i)
      AAABR(405) = EE*sb*(2D0*Zm12i/SW*Zn31i+2D0*Zm12r/SW
     &      *Zn31r-Sqrt2/CW*Zm22i*Zn11i-Sqrt2/CW*Zm22r*Zn11r
     &      -Sqrt2/SW*Zm22i*Zn21i-Sqrt2/SW*Zm22r*Zn21r)
      AAABR(406) = EE*sb*(2D0*Zm12i/SW*Zn31r-2D0*Zm12r/SW
     &      *Zn31i-Sqrt2/CW*Zm22i*Zn11r+Sqrt2/CW*Zm22r*Zn11i
     &      -Sqrt2/SW*Zm22i*Zn21r+Sqrt2/SW*Zm22r*Zn21i)
      AAABR(407) = EE*cb*(2D0*Zn41i/SW*Zp12i-2D0*Zn41r/SW
     &      *Zp12r+Sqrt2/CW*Zn11i*Zp22i-Sqrt2/CW*Zn11r*Zp22r
     &      +Sqrt2/SW*Zn21i*Zp22i-Sqrt2/SW*Zn21r*Zp22r)
      AAABR(408) = EE*cb*(2D0*Zn41i/SW*Zp12r+2D0*Zn41r/SW
     &      *Zp12i+Sqrt2/CW*Zn11i*Zp22r+Sqrt2/CW*Zn11r*Zp22i
     &      +Sqrt2/SW*Zn21i*Zp22r+Sqrt2/SW*Zn21r*Zp22i)
      AAABR(409) = EE*sb*(2D0*Zm12i/SW*Zn34i+2D0*Zm12r/SW
     &      *Zn34r-Sqrt2/CW*Zm22i*Zn14i-Sqrt2/CW*Zm22r*Zn14r
     &      -Sqrt2/SW*Zm22i*Zn24i-Sqrt2/SW*Zm22r*Zn24r)
      AAABR(410) = EE*sb*(2D0*Zm12i/SW*Zn34r-2D0*Zm12r/SW
     &      *Zn34i-Sqrt2/CW*Zm22i*Zn14r+Sqrt2/CW*Zm22r*Zn14i
     &      -Sqrt2/SW*Zm22i*Zn24r+Sqrt2/SW*Zm22r*Zn24i)
      AAABR(411) = EE*cb*(2D0*Zn44i/SW*Zp12i-2D0*Zn44r/SW
     &      *Zp12r+Sqrt2/CW*Zn14i*Zp22i-Sqrt2/CW*Zn14r*Zp22r
     &      +Sqrt2/SW*Zn24i*Zp22i-Sqrt2/SW*Zn24r*Zp22r)
      AAABR(412) = EE*cb*(2D0*Zn44i/SW*Zp12r+2D0*Zn44r/SW
     &      *Zp12i+Sqrt2/CW*Zn14i*Zp22r+Sqrt2/CW*Zn14r*Zp22i
     &      +Sqrt2/SW*Zn24i*Zp22r+Sqrt2/SW*Zn24r*Zp22i)
      AAABR(413) = EE*cb*(2D0*Zm12i/SW*Zn33i+2D0*Zm12r/SW
     &      *Zn33r-Sqrt2/CW*Zm22i*Zn13i-Sqrt2/CW*Zm22r*Zn13r
     &      -Sqrt2/SW*Zm22i*Zn23i-Sqrt2/SW*Zm22r*Zn23r)
      AAABR(414) = EE*cb*(2D0*Zm12i/SW*Zn33r-2D0*Zm12r/SW
     &      *Zn33i-Sqrt2/CW*Zm22i*Zn13r+Sqrt2/CW*Zm22r*Zn13i
     &      -Sqrt2/SW*Zm22i*Zn23r+Sqrt2/SW*Zm22r*Zn23i)
      AAABR(415) = EE*sb*(2D0*Zn43i/SW*Zp12i-2D0*Zn43r/SW
     &      *Zp12r+Sqrt2/CW*Zn13i*Zp22i-Sqrt2/CW*Zn13r*Zp22r
     &      +Sqrt2/SW*Zn23i*Zp22i-Sqrt2/SW*Zn23r*Zp22r)
      AAABR(416) = EE*sb*(2D0*Zn43i/SW*Zp12r+2D0*Zn43r/SW
     &      *Zp12i+Sqrt2/CW*Zn13i*Zp22r+Sqrt2/CW*Zn13r*Zp22i
     &      +Sqrt2/SW*Zn23i*Zp22r+Sqrt2/SW*Zn23r*Zp22i)
      AAABR(417) = EE*cb*(2D0*Zm12i/SW*Zn32i+2D0*Zm12r/SW
     &      *Zn32r-Sqrt2/CW*Zm22i*Zn12i-Sqrt2/CW*Zm22r*Zn12r
     &      -Sqrt2/SW*Zm22i*Zn22i-Sqrt2/SW*Zm22r*Zn22r)
      AAABR(418) = EE*cb*(2D0*Zm12i/SW*Zn32r-2D0*Zm12r/SW
     &      *Zn32i-Sqrt2/CW*Zm22i*Zn12r+Sqrt2/CW*Zm22r*Zn12i
     &      -Sqrt2/SW*Zm22i*Zn22r+Sqrt2/SW*Zm22r*Zn22i)
      AAABR(419) = EE*sb*(2D0*Zn42i/SW*Zp12i-2D0*Zn42r/SW
     &      *Zp12r+Sqrt2/CW*Zn12i*Zp22i-Sqrt2/CW*Zn12r*Zp22r
     &      +Sqrt2/SW*Zn22i*Zp22i-Sqrt2/SW*Zn22r*Zp22r)
      AAABR(420) = EE*sb*(2D0*Zn42i/SW*Zp12r+2D0*Zn42r/SW
     &      *Zp12i+Sqrt2/CW*Zn12i*Zp22r+Sqrt2/CW*Zn12r*Zp22i
     &      +Sqrt2/SW*Zn22i*Zp22r+Sqrt2/SW*Zn22r*Zp22i)
      AAABR(421) = EE*cb*(2D0*Zm11i/SW*Zn34i+2D0*Zm11r/SW
     &      *Zn34r-Sqrt2/CW*Zm21i*Zn14i-Sqrt2/CW*Zm21r*Zn14r
     &      -Sqrt2/SW*Zm21i*Zn24i-Sqrt2/SW*Zm21r*Zn24r)
      AAABR(422) = EE*cb*(2D0*Zm11i/SW*Zn34r-2D0*Zm11r/SW
     &      *Zn34i-Sqrt2/CW*Zm21i*Zn14r+Sqrt2/CW*Zm21r*Zn14i
     &      -Sqrt2/SW*Zm21i*Zn24r+Sqrt2/SW*Zm21r*Zn24i)
      AAABR(423) = EE*sb*(2D0*Zn44i/SW*Zp11i-2D0*Zn44r/SW
     &      *Zp11r+Sqrt2/CW*Zn14i*Zp21i-Sqrt2/CW*Zn14r*Zp21r
     &      +Sqrt2/SW*Zn24i*Zp21i-Sqrt2/SW*Zn24r*Zp21r)
      AAABR(424) = EE*sb*(2D0*Zn44i/SW*Zp11r+2D0*Zn44r/SW
     &      *Zp11i+Sqrt2/CW*Zn14i*Zp21r+Sqrt2/CW*Zn14r*Zp21i
     &      +Sqrt2/SW*Zn24i*Zp21r+Sqrt2/SW*Zn24r*Zp21i)
      AAABR(425) = EE*cb*(2D0*Zm11i/SW*Zn33i+2D0*Zm11r/SW
     &      *Zn33r-Sqrt2/CW*Zm21i*Zn13i-Sqrt2/CW*Zm21r*Zn13r
     &      -Sqrt2/SW*Zm21i*Zn23i-Sqrt2/SW*Zm21r*Zn23r)
      AAABR(426) = EE*cb*(2D0*Zm11i/SW*Zn33r-2D0*Zm11r/SW
     &      *Zn33i-Sqrt2/CW*Zm21i*Zn13r+Sqrt2/CW*Zm21r*Zn13i
     &      -Sqrt2/SW*Zm21i*Zn23r+Sqrt2/SW*Zm21r*Zn23i)
      AAABR(427) = EE*sb*(2D0*Zn43i/SW*Zp11i-2D0*Zn43r/SW
     &      *Zp11r+Sqrt2/CW*Zn13i*Zp21i-Sqrt2/CW*Zn13r*Zp21r
     &      +Sqrt2/SW*Zn23i*Zp21i-Sqrt2/SW*Zn23r*Zp21r)
      AAABR(428) = EE*sb*(2D0*Zn43i/SW*Zp11r+2D0*Zn43r/SW
     &      *Zp11i+Sqrt2/CW*Zn13i*Zp21r+Sqrt2/CW*Zn13r*Zp21i
     &      +Sqrt2/SW*Zn23i*Zp21r+Sqrt2/SW*Zn23r*Zp21i)
      AAABR(429) = EE*cb*(2D0*Zm12i/SW*Zn31i+2D0*Zm12r/SW
     &      *Zn31r-Sqrt2/CW*Zm22i*Zn11i-Sqrt2/CW*Zm22r*Zn11r
     &      -Sqrt2/SW*Zm22i*Zn21i-Sqrt2/SW*Zm22r*Zn21r)
      AAABR(430) = EE*cb*(2D0*Zm12i/SW*Zn31r-2D0*Zm12r/SW
     &      *Zn31i-Sqrt2/CW*Zm22i*Zn11r+Sqrt2/CW*Zm22r*Zn11i
     &      -Sqrt2/SW*Zm22i*Zn21r+Sqrt2/SW*Zm22r*Zn21i)
      AAABR(431) = EE*sb*(2D0*Zn41i/SW*Zp12i-2D0*Zn41r/SW
     &      *Zp12r+Sqrt2/CW*Zn11i*Zp22i-Sqrt2/CW*Zn11r*Zp22r
     &      +Sqrt2/SW*Zn21i*Zp22i-Sqrt2/SW*Zn21r*Zp22r)
      AAABR(432) = EE*sb*(2D0*Zn41i/SW*Zp12r+2D0*Zn41r/SW
     &      *Zp12i+Sqrt2/CW*Zn11i*Zp22r+Sqrt2/CW*Zn11r*Zp22i
     &      +Sqrt2/SW*Zn21i*Zp22r+Sqrt2/SW*Zn21r*Zp22i)
      AAABR(433) = EE*cb*(2D0*Zm11i/SW*Zn31i+2D0*Zm11r/SW
     &      *Zn31r-Sqrt2/CW*Zm21i*Zn11i-Sqrt2/CW*Zm21r*Zn11r
     &      -Sqrt2/SW*Zm21i*Zn21i-Sqrt2/SW*Zm21r*Zn21r)
      AAABR(434) = EE*cb*(2D0*Zm11i/SW*Zn31r-2D0*Zm11r/SW
     &      *Zn31i-Sqrt2/CW*Zm21i*Zn11r+Sqrt2/CW*Zm21r*Zn11i
     &      -Sqrt2/SW*Zm21i*Zn21r+Sqrt2/SW*Zm21r*Zn21i)
      AAABR(435) = EE*sb*(2D0*Zn41i/SW*Zp11i-2D0*Zn41r/SW
     &      *Zp11r+Sqrt2/CW*Zn11i*Zp21i-Sqrt2/CW*Zn11r*Zp21r
     &      +Sqrt2/SW*Zn21i*Zp21i-Sqrt2/SW*Zn21r*Zp21r)
      AAABR(436) = EE*sb*(2D0*Zn41i/SW*Zp11r+2D0*Zn41r/SW
     &      *Zp11i+Sqrt2/CW*Zn11i*Zp21r+Sqrt2/CW*Zn11r*Zp21i
     &      +Sqrt2/SW*Zn21i*Zp21r+Sqrt2/SW*Zn21r*Zp21i)
      AAABR(437) = EE*cb*(2D0*Zm11i/SW*Zn32i+2D0*Zm11r/SW
     &      *Zn32r-Sqrt2/CW*Zm21i*Zn12i-Sqrt2/CW*Zm21r*Zn12r
     &      -Sqrt2/SW*Zm21i*Zn22i-Sqrt2/SW*Zm21r*Zn22r)
      AAABR(438) = EE*cb*(2D0*Zm11i/SW*Zn32r-2D0*Zm11r/SW
     &      *Zn32i-Sqrt2/CW*Zm21i*Zn12r+Sqrt2/CW*Zm21r*Zn12i
     &      -Sqrt2/SW*Zm21i*Zn22r+Sqrt2/SW*Zm21r*Zn22i)
      AAABR(439) = EE*sb*(2D0*Zn42i/SW*Zp11i-2D0*Zn42r/SW
     &      *Zp11r+Sqrt2/CW*Zn12i*Zp21i-Sqrt2/CW*Zn12r*Zp21r
     &      +Sqrt2/SW*Zn22i*Zp21i-Sqrt2/SW*Zn22r*Zp21r)
      AAABR(440) = EE*sb*(2D0*Zn42i/SW*Zp11r+2D0*Zn42r/SW
     &      *Zp11i+Sqrt2/CW*Zn12i*Zp21r+Sqrt2/CW*Zn12r*Zp21i
     &      +Sqrt2/SW*Zn22i*Zp21r+Sqrt2/SW*Zn22r*Zp21i)
      AAABR(441) = EE*cb*(2D0*Zm12i/SW*Zn34i+2D0*Zm12r/SW
     &      *Zn34r-Sqrt2/CW*Zm22i*Zn14i-Sqrt2/CW*Zm22r*Zn14r
     &      -Sqrt2/SW*Zm22i*Zn24i-Sqrt2/SW*Zm22r*Zn24r)
      AAABR(442) = EE*cb*(2D0*Zm12i/SW*Zn34r-2D0*Zm12r/SW
     &      *Zn34i-Sqrt2/CW*Zm22i*Zn14r+Sqrt2/CW*Zm22r*Zn14i
     &      -Sqrt2/SW*Zm22i*Zn24r+Sqrt2/SW*Zm22r*Zn24i)
      AAABR(443) = EE*sb*(2D0*Zn44i/SW*Zp12i-2D0*Zn44r/SW
     &      *Zp12r+Sqrt2/CW*Zn14i*Zp22i-Sqrt2/CW*Zn14r*Zp22r
     &      +Sqrt2/SW*Zn24i*Zp22i-Sqrt2/SW*Zn24r*Zp22r)
      AAABR(444) = EE*sb*(2D0*Zn44i/SW*Zp12r+2D0*Zn44r/SW
     &      *Zp12i+Sqrt2/CW*Zn14i*Zp22r+Sqrt2/CW*Zn14r*Zp22i
     &      +Sqrt2/SW*Zn24i*Zp22r+Sqrt2/SW*Zn24r*Zp22i)
      AAABR(445) = EE/SW*(2D0*Zd33i*Zm11i-2D0*Zd33r*Zm11r
     &      -Mb/MW*Sqrt2*Zd63i*Zm21i/cb+Mb/MW*Sqrt2*Zd63r
     &      *Zm21r/cb)
      AAABR(446) = EE/SW*(2D0*Zd33i*Zm11r+2D0*Zd33r*Zm11i
     &      -Mb/MW*Sqrt2*Zd63i*Zm21r/cb-Mb/MW*Sqrt2*Zd63r
     &      *Zm21i/cb)
      AAABR(447) = EE/MW*Mt/SW*Sqrt2/sb*(Zd33i*Zp21i-Zd33r
     &      *Zp21r)
      AAABR(448) = EE/MW*Mt/SW*Sqrt2/sb*(Zd33i*Zp21r+Zd33r
     &      *Zp21i)
      AAABR(449) = EE/SW*Zd11*Zm12r
      AAABR(450) = EE/SW*Zd11*Zm12i
      AAABR(451) = EE/SW*(2D0*Zd22*Zm11r-Ms/MW*Sqrt2*Zd52
     &      *Zm21r/cb)
      AAABR(452) = EE/SW*(2D0*Zd22*Zm11i-Ms/MW*Sqrt2*Zd52
     &      *Zm21i/cb)
      AAABR(453) = EE/MW*Mc/SW*Sqrt2*Zd22*Zp21r/sb
      AAABR(454) = EE/MW*Mc/SW*Sqrt2*Zd22*Zp21i/sb
      AAABR(455) = EE/SW*Zd11*Zm11r
      AAABR(456) = EE/SW*Zd11*Zm11i
      AAABR(457) = EE/SW*(2D0*Zd22*Zm12r-Ms/MW*Sqrt2*Zd52
     &      *Zm22r/cb)
      AAABR(458) = EE/SW*(2D0*Zd22*Zm12i-Ms/MW*Sqrt2*Zd52
     &      *Zm22i/cb)
      AAABR(459) = EE/MW*Mc/SW*Sqrt2*Zd22*Zp22r/sb
      AAABR(460) = EE/MW*Mc/SW*Sqrt2*Zd22*Zp22i/sb
      AAABR(461) = EE/SW*(2D0*Zd33i*Zm12i-2D0*Zd33r*Zm12r
     &      -Mb/MW*Sqrt2*Zd63i*Zm22i/cb+Mb/MW*Sqrt2*Zd63r
     &      *Zm22r/cb)
      AAABR(462) = EE/SW*(2D0*Zd33i*Zm12r+2D0*Zd33r*Zm12i
     &      -Mb/MW*Sqrt2*Zd63i*Zm22r/cb-Mb/MW*Sqrt2*Zd63r
     &      *Zm22i/cb)
      AAABR(463) = EE/MW*Mt/SW*Sqrt2/sb*(Zd33i*Zp22i-Zd33r
     &      *Zp22r)
      AAABR(464) = EE/MW*Mt/SW*Sqrt2/sb*(Zd33i*Zp22r+Zd33r
     &      *Zp22i)
      AAABR(465) = EE/SW*Zd14*Zm12r
      AAABR(466) = EE/SW*Zd14*Zm12i
      AAABR(467) = EE/SW*(2D0*Zd36i*Zm12i-2D0*Zd36r*Zm12r
     &      -Mb/MW*Sqrt2*Zd66i*Zm22i/cb+Mb/MW*Sqrt2*Zd66r
     &      *Zm22r/cb)
      AAABR(468) = EE/SW*(2D0*Zd36i*Zm12r+2D0*Zd36r*Zm12i
     &      -Mb/MW*Sqrt2*Zd66i*Zm22r/cb-Mb/MW*Sqrt2*Zd66r
     &      *Zm22i/cb)
      AAABR(469) = EE/MW*Mt/SW*Sqrt2/sb*(Zd36i*Zp22i-Zd36r
     &      *Zp22r)
      AAABR(470) = EE/MW*Mt/SW*Sqrt2/sb*(Zd36i*Zp22r+Zd36r
     &      *Zp22i)
      AAABR(471) = EE/SW*(2D0*Zd36i*Zm11i-2D0*Zd36r*Zm11r
     &      -Mb/MW*Sqrt2*Zd66i*Zm21i/cb+Mb/MW*Sqrt2*Zd66r
     &      *Zm21r/cb)
      AAABR(472) = EE/SW*(2D0*Zd36i*Zm11r+2D0*Zd36r*Zm11i
     &      -Mb/MW*Sqrt2*Zd66i*Zm21r/cb-Mb/MW*Sqrt2*Zd66r
     &      *Zm21i/cb)
      AAABR(473) = EE/MW*Mt/SW*Sqrt2/sb*(Zd36i*Zp21i-Zd36r
     &      *Zp21r)
      AAABR(474) = EE/MW*Mt/SW*Sqrt2/sb*(Zd36i*Zp21r+Zd36r
     &      *Zp21i)
      AAABR(475) = EE/SW*Zd14*Zm11r
      AAABR(476) = EE/SW*Zd14*Zm11i
      AAABR(477) = EE/SW*(2D0*Zd25*Zm11r-Ms/MW*Sqrt2*Zd55
     &      *Zm21r/cb)
      AAABR(478) = EE/SW*(2D0*Zd25*Zm11i-Ms/MW*Sqrt2*Zd55
     &      *Zm21i/cb)
      AAABR(479) = EE/MW*Mc/SW*Sqrt2*Zd25*Zp21r/sb
      AAABR(480) = EE/MW*Mc/SW*Sqrt2*Zd25*Zp21i/sb
      AAABR(481) = EE/SW*(2D0*Zd25*Zm12r-Ms/MW*Sqrt2*Zd55
     &      *Zm22r/cb)
      AAABR(482) = EE/SW*(2D0*Zd25*Zm12i-Ms/MW*Sqrt2*Zd55
     &      *Zm22i/cb)
      AAABR(483) = EE/MW*Mc/SW*Sqrt2*Zd25*Zp22r/sb
      AAABR(484) = EE/MW*Mc/SW*Sqrt2*Zd25*Zp22i/sb
      AAABR(485) = EE/MW*Mm/SW
      AAABR(486) = EE/MW*Me/SW
      AAABR(487) = EE/MW*Ml/SW
      AAABR(488) = EE/MW*Mm/SW*Zh11/cb
      AAABR(489) = EE/MW*Mm/SW*Zh31/cb*sb
      AAABR(490) = EE/MW*Ml/SW*Zh13/cb
      AAABR(491) = EE/MW*Ml/SW*Zh33/cb*sb
      AAABR(492) = EE/MW*Mm/SW*Zh13/cb
      AAABR(493) = EE/MW*Mm/SW*Zh33/cb*sb
      AAABR(494) = EE/MW*Me/SW*Zh13/cb
      AAABR(495) = EE/MW*Me/SW*Zh33/cb*sb
      AAABR(496) = EE/MW*Mm/SW*Zh12/cb
      AAABR(497) = EE/MW*Mm/SW*Zh32/cb*sb
      AAABR(498) = EE/MW*Me/SW*Zh11/cb
      AAABR(499) = EE/MW*Me/SW*Zh31/cb*sb
      AAABR(500) = EE/MW*Me/SW*Zh12/cb
      AAABR(501) = EE/MW*Me/SW*Zh32/cb*sb
      AAABR(502) = EE/MW*Ml/SW*Zh11/cb
      AAABR(503) = EE/MW*Ml/SW*Zh31/cb*sb
      AAABR(504) = EE/MW*Ml/SW*Zh12/cb
      AAABR(505) = EE/MW*Ml/SW*Zh32/cb*sb
      AAABR(506) = EE/MW*Me/SW*Sqrt2*tb
      AAABR(507) = EE/MW*Mm/SW*Sqrt2*tb
      AAABR(508) = EE/MW*Ml/SW*Sqrt2*tb
      AAABR(509) = EE/MW*Ml/SW*Sqrt2
      AAABR(510) = EE/MW*Me/SW*Sqrt2
      AAABR(511) = EE/MW*Mm/SW*Sqrt2
      AAABR(512) = EE*Sqrt2*(Zl11/CW*Zn12r+Zl11/SW*Zn22r
     &      -Me/MW/SW*Zl41*Zn32r/cb)
      AAABR(513) = EE*Sqrt2*(Zl11/CW*Zn12i+Zl11/SW*Zn22i
     &      -Me/MW/SW*Zl41*Zn32i/cb)
      AAABR(514) = EE*Sqrt2*(2D0*Zl41/CW*Zn12r+Me/MW/SW*Zl11
     &      *Zn32r/cb)
      AAABR(515) = EE*Sqrt2*(2D0*Zl41/CW*Zn12i+Me/MW/SW*Zl11
     &      *Zn32i/cb)
      AAABR(516) = EE*Sqrt2*(Zl11/CW*Zn11r+Zl11/SW*Zn21r
     &      -Me/MW/SW*Zl41*Zn31r/cb)
      AAABR(517) = EE*Sqrt2*(Zl11/CW*Zn11i+Zl11/SW*Zn21i
     &      -Me/MW/SW*Zl41*Zn31i/cb)
      AAABR(518) = EE*Sqrt2*(2D0*Zl41/CW*Zn11r+Me/MW/SW*Zl11
     &      *Zn31r/cb)
      AAABR(519) = EE*Sqrt2*(2D0*Zl41/CW*Zn11i+Me/MW/SW*Zl11
     &      *Zn31i/cb)
      AAABR(520) = EE*Sqrt2*(Zl22/CW*Zn14r+Zl22/SW*Zn24r
     &      -Mm/MW/SW*Zl52*Zn34r/cb)
      AAABR(521) = EE*Sqrt2*(Zl22/CW*Zn14i+Zl22/SW*Zn24i
     &      -Mm/MW/SW*Zl52*Zn34i/cb)
      AAABR(522) = EE*Sqrt2*(2D0*Zl52/CW*Zn14r+Mm/MW/SW*Zl22
     &      *Zn34r/cb)
      AAABR(523) = EE*Sqrt2*(2D0*Zl52/CW*Zn14i+Mm/MW/SW*Zl22
     &      *Zn34i/cb)
      AAABR(524) = EE*Sqrt2*(Zl11/CW*Zn13r+Zl11/SW*Zn23r
     &      -Me/MW/SW*Zl41*Zn33r/cb)
      AAABR(525) = EE*Sqrt2*(Zl11/CW*Zn13i+Zl11/SW*Zn23i
     &      -Me/MW/SW*Zl41*Zn33i/cb)
      AAABR(526) = EE*Sqrt2*(2D0*Zl41/CW*Zn13r+Me/MW/SW*Zl11
     &      *Zn33r/cb)
      AAABR(527) = EE*Sqrt2*(2D0*Zl41/CW*Zn13i+Me/MW/SW*Zl11
     &      *Zn33i/cb)
      AAABR(528) = EE*Sqrt2*(Zl33i/CW*Zn14i+Zl33r/CW*Zn14r
     &      +Zl33i/SW*Zn24i+Zl33r/SW*Zn24r-Ml/MW/SW*Zl63i
     &      *Zn34i/cb-Ml/MW/SW*Zl63r*Zn34r/cb)
      AAABR(529) = EE*Sqrt2*(Zl33i/CW*Zn14r-Zl33r/CW*Zn14i
     &      +Zl33i/SW*Zn24r-Zl33r/SW*Zn24i-Ml/MW/SW*Zl63i
     &      *Zn34r/cb+Ml/MW/SW*Zl63r*Zn34i/cb)
      AAABR(530) = EE*Sqrt2*(2D0*Zl63i/CW*Zn14i-2D0*Zl63r
     &      /CW*Zn14r+Ml/MW/SW*Zl33i*Zn34i/cb-Ml/MW/SW*Zl33r
     &      *Zn34r/cb)
      AAABR(531) = EE*Sqrt2*(2D0*Zl63i/CW*Zn14r+2D0*Zl63r
     &      /CW*Zn14i+Ml/MW/SW*Zl33i*Zn34r/cb+Ml/MW/SW*Zl33r
     &      *Zn34i/cb)
      AAABR(532) = EE*Sqrt2*(Zl33i/CW*Zn11i+Zl33r/CW*Zn11r
     &      +Zl33i/SW*Zn21i+Zl33r/SW*Zn21r-Ml/MW/SW*Zl63i
     &      *Zn31i/cb-Ml/MW/SW*Zl63r*Zn31r/cb)
      AAABR(533) = EE*Sqrt2*(Zl33i/CW*Zn11r-Zl33r/CW*Zn11i
     &      +Zl33i/SW*Zn21r-Zl33r/SW*Zn21i-Ml/MW/SW*Zl63i
     &      *Zn31r/cb+Ml/MW/SW*Zl63r*Zn31i/cb)
      AAABR(534) = EE*Sqrt2*(2D0*Zl63i/CW*Zn11i-2D0*Zl63r
     &      /CW*Zn11r+Ml/MW/SW*Zl33i*Zn31i/cb-Ml/MW/SW*Zl33r
     &      *Zn31r/cb)
      AAABR(535) = EE*Sqrt2*(2D0*Zl63i/CW*Zn11r+2D0*Zl63r
     &      /CW*Zn11i+Ml/MW/SW*Zl33i*Zn31r/cb+Ml/MW/SW*Zl33r
     &      *Zn31i/cb)
      AAABR(536) = EE*Sqrt2*(Zl22/CW*Zn13r+Zl22/SW*Zn23r
     &      -Mm/MW/SW*Zl52*Zn33r/cb)
      AAABR(537) = EE*Sqrt2*(Zl22/CW*Zn13i+Zl22/SW*Zn23i
     &      -Mm/MW/SW*Zl52*Zn33i/cb)
      AAABR(538) = EE*Sqrt2*(2D0*Zl52/CW*Zn13r+Mm/MW/SW*Zl22
     &      *Zn33r/cb)
      AAABR(539) = EE*Sqrt2*(2D0*Zl52/CW*Zn13i+Mm/MW/SW*Zl22
     &      *Zn33i/cb)
      AAABR(540) = EE*Sqrt2*(Zl33i/CW*Zn12i+Zl33r/CW*Zn12r
     &      +Zl33i/SW*Zn22i+Zl33r/SW*Zn22r-Ml/MW/SW*Zl63i
     &      *Zn32i/cb-Ml/MW/SW*Zl63r*Zn32r/cb)
      AAABR(541) = EE*Sqrt2*(Zl33i/CW*Zn12r-Zl33r/CW*Zn12i
     &      +Zl33i/SW*Zn22r-Zl33r/SW*Zn22i-Ml/MW/SW*Zl63i
     &      *Zn32r/cb+Ml/MW/SW*Zl63r*Zn32i/cb)
      AAABR(542) = EE*Sqrt2*(2D0*Zl63i/CW*Zn12i-2D0*Zl63r
     &      /CW*Zn12r+Ml/MW/SW*Zl33i*Zn32i/cb-Ml/MW/SW*Zl33r
     &      *Zn32r/cb)
      AAABR(543) = EE*Sqrt2*(2D0*Zl63i/CW*Zn12r+2D0*Zl63r
     &      /CW*Zn12i+Ml/MW/SW*Zl33i*Zn32r/cb+Ml/MW/SW*Zl33r
     &      *Zn32i/cb)
      AAABR(544) = EE*Sqrt2*(Zl33i/CW*Zn13i+Zl33r/CW*Zn13r
     &      +Zl33i/SW*Zn23i+Zl33r/SW*Zn23r-Ml/MW/SW*Zl63i
     &      *Zn33i/cb-Ml/MW/SW*Zl63r*Zn33r/cb)
      AAABR(545) = EE*Sqrt2*(Zl33i/CW*Zn13r-Zl33r/CW*Zn13i
     &      +Zl33i/SW*Zn23r-Zl33r/SW*Zn23i-Ml/MW/SW*Zl63i
     &      *Zn33r/cb+Ml/MW/SW*Zl63r*Zn33i/cb)
      AAABR(546) = EE*Sqrt2*(2D0*Zl63i/CW*Zn13i-2D0*Zl63r
     &      /CW*Zn13r+Ml/MW/SW*Zl33i*Zn33i/cb-Ml/MW/SW*Zl33r
     &      *Zn33r/cb)
      AAABR(547) = EE*Sqrt2*(2D0*Zl63i/CW*Zn13r+2D0*Zl63r
     &      /CW*Zn13i+Ml/MW/SW*Zl33i*Zn33r/cb+Ml/MW/SW*Zl33r
     &      *Zn33i/cb)
      AAABR(548) = EE*Sqrt2*(Zl11/CW*Zn14r+Zl11/SW*Zn24r
     &      -Me/MW/SW*Zl41*Zn34r/cb)
      AAABR(549) = EE*Sqrt2*(Zl11/CW*Zn14i+Zl11/SW*Zn24i
     &      -Me/MW/SW*Zl41*Zn34i/cb)
      AAABR(550) = EE*Sqrt2*(2D0*Zl41/CW*Zn14r+Me/MW/SW*Zl11
     &      *Zn34r/cb)
      AAABR(551) = EE*Sqrt2*(2D0*Zl41/CW*Zn14i+Me/MW/SW*Zl11
     &      *Zn34i/cb)
      AAABR(552) = EE*Sqrt2*(Zl22/CW*Zn12r+Zl22/SW*Zn22r
     &      -Mm/MW/SW*Zl52*Zn32r/cb)
      AAABR(553) = EE*Sqrt2*(Zl22/CW*Zn12i+Zl22/SW*Zn22i
     &      -Mm/MW/SW*Zl52*Zn32i/cb)
      AAABR(554) = EE*Sqrt2*(2D0*Zl52/CW*Zn12r+Mm/MW/SW*Zl22
     &      *Zn32r/cb)
      AAABR(555) = EE*Sqrt2*(2D0*Zl52/CW*Zn12i+Mm/MW/SW*Zl22
     &      *Zn32i/cb)
      AAABR(556) = EE*Sqrt2*(Zl22/CW*Zn11r+Zl22/SW*Zn21r
     &      -Mm/MW/SW*Zl52*Zn31r/cb)
      AAABR(557) = EE*Sqrt2*(Zl22/CW*Zn11i+Zl22/SW*Zn21i
     &      -Mm/MW/SW*Zl52*Zn31i/cb)
      AAABR(558) = EE*Sqrt2*(2D0*Zl52/CW*Zn11r+Mm/MW/SW*Zl22
     &      *Zn31r/cb)
      AAABR(559) = EE*Sqrt2*(2D0*Zl52/CW*Zn11i+Mm/MW/SW*Zl22
     &      *Zn31i/cb)
      AAABR(560) = EE*Sqrt2*(Zl14/CW*Zn13r+Zl14/SW*Zn23r
     &      -Me/MW/SW*Zl44*Zn33r/cb)
      AAABR(561) = EE*Sqrt2*(Zl14/CW*Zn13i+Zl14/SW*Zn23i
     &      -Me/MW/SW*Zl44*Zn33i/cb)
      AAABR(562) = EE*Sqrt2*(2D0*Zl44/CW*Zn13r+Me/MW/SW*Zl14
     &      *Zn33r/cb)
      AAABR(563) = EE*Sqrt2*(2D0*Zl44/CW*Zn13i+Me/MW/SW*Zl14
     &      *Zn33i/cb)
      AAABR(564) = EE*Sqrt2*(Zl36i/CW*Zn14i+Zl36r/CW*Zn14r
     &      +Zl36i/SW*Zn24i+Zl36r/SW*Zn24r-Ml/MW/SW*Zl66i
     &      *Zn34i/cb-Ml/MW/SW*Zl66r*Zn34r/cb)
      AAABR(565) = EE*Sqrt2*(Zl36i/CW*Zn14r-Zl36r/CW*Zn14i
     &      +Zl36i/SW*Zn24r-Zl36r/SW*Zn24i-Ml/MW/SW*Zl66i
     &      *Zn34r/cb+Ml/MW/SW*Zl66r*Zn34i/cb)
      AAABR(566) = EE*Sqrt2*(2D0*Zl66i/CW*Zn14i-2D0*Zl66r
     &      /CW*Zn14r+Ml/MW/SW*Zl36i*Zn34i/cb-Ml/MW/SW*Zl36r
     &      *Zn34r/cb)
      AAABR(567) = EE*Sqrt2*(2D0*Zl66i/CW*Zn14r+2D0*Zl66r
     &      /CW*Zn14i+Ml/MW/SW*Zl36i*Zn34r/cb+Ml/MW/SW*Zl36r
     &      *Zn34i/cb)
      AAABR(568) = EE*Sqrt2*(Zl36i/CW*Zn13i+Zl36r/CW*Zn13r
     &      +Zl36i/SW*Zn23i+Zl36r/SW*Zn23r-Ml/MW/SW*Zl66i
     &      *Zn33i/cb-Ml/MW/SW*Zl66r*Zn33r/cb)
      AAABR(569) = EE*Sqrt2*(Zl36i/CW*Zn13r-Zl36r/CW*Zn13i
     &      +Zl36i/SW*Zn23r-Zl36r/SW*Zn23i-Ml/MW/SW*Zl66i
     &      *Zn33r/cb+Ml/MW/SW*Zl66r*Zn33i/cb)
      AAABR(570) = EE*Sqrt2*(2D0*Zl66i/CW*Zn13i-2D0*Zl66r
     &      /CW*Zn13r+Ml/MW/SW*Zl36i*Zn33i/cb-Ml/MW/SW*Zl36r
     &      *Zn33r/cb)
      AAABR(571) = EE*Sqrt2*(2D0*Zl66i/CW*Zn13r+2D0*Zl66r
     &      /CW*Zn13i+Ml/MW/SW*Zl36i*Zn33r/cb+Ml/MW/SW*Zl36r
     &      *Zn33i/cb)
      AAABR(572) = EE*Sqrt2*(Zl14/CW*Zn14r+Zl14/SW*Zn24r
     &      -Me/MW/SW*Zl44*Zn34r/cb)
      AAABR(573) = EE*Sqrt2*(Zl14/CW*Zn14i+Zl14/SW*Zn24i
     &      -Me/MW/SW*Zl44*Zn34i/cb)
      AAABR(574) = EE*Sqrt2*(2D0*Zl44/CW*Zn14r+Me/MW/SW*Zl14
     &      *Zn34r/cb)
      AAABR(575) = EE*Sqrt2*(2D0*Zl44/CW*Zn14i+Me/MW/SW*Zl14
     &      *Zn34i/cb)
      AAABR(576) = EE*Sqrt2*(Zl25/CW*Zn11r+Zl25/SW*Zn21r
     &      -Mm/MW/SW*Zl55*Zn31r/cb)
      AAABR(577) = EE*Sqrt2*(Zl25/CW*Zn11i+Zl25/SW*Zn21i
     &      -Mm/MW/SW*Zl55*Zn31i/cb)
      AAABR(578) = EE*Sqrt2*(2D0*Zl55/CW*Zn11r+Mm/MW/SW*Zl25
     &      *Zn31r/cb)
      AAABR(579) = EE*Sqrt2*(2D0*Zl55/CW*Zn11i+Mm/MW/SW*Zl25
     &      *Zn31i/cb)
      AAABR(580) = EE*Sqrt2*(Zl25/CW*Zn14r+Zl25/SW*Zn24r
     &      -Mm/MW/SW*Zl55*Zn34r/cb)
      AAABR(581) = EE*Sqrt2*(Zl25/CW*Zn14i+Zl25/SW*Zn24i
     &      -Mm/MW/SW*Zl55*Zn34i/cb)
      AAABR(582) = EE*Sqrt2*(2D0*Zl55/CW*Zn14r+Mm/MW/SW*Zl25
     &      *Zn34r/cb)
      AAABR(583) = EE*Sqrt2*(2D0*Zl55/CW*Zn14i+Mm/MW/SW*Zl25
     &      *Zn34i/cb)
      AAABR(584) = EE*Sqrt2*(Zl25/CW*Zn13r+Zl25/SW*Zn23r
     &      -Mm/MW/SW*Zl55*Zn33r/cb)
      AAABR(585) = EE*Sqrt2*(Zl25/CW*Zn13i+Zl25/SW*Zn23i
     &      -Mm/MW/SW*Zl55*Zn33i/cb)
      AAABR(586) = EE*Sqrt2*(2D0*Zl55/CW*Zn13r+Mm/MW/SW*Zl25
     &      *Zn33r/cb)
      AAABR(587) = EE*Sqrt2*(2D0*Zl55/CW*Zn13i+Mm/MW/SW*Zl25
     &      *Zn33i/cb)
      AAABR(588) = EE*Sqrt2*(Zl25/CW*Zn12r+Zl25/SW*Zn22r
     &      -Mm/MW/SW*Zl55*Zn32r/cb)
      AAABR(589) = EE*Sqrt2*(Zl25/CW*Zn12i+Zl25/SW*Zn22i
     &      -Mm/MW/SW*Zl55*Zn32i/cb)
      AAABR(590) = EE*Sqrt2*(2D0*Zl55/CW*Zn12r+Mm/MW/SW*Zl25
     &      *Zn32r/cb)
      AAABR(591) = EE*Sqrt2*(2D0*Zl55/CW*Zn12i+Mm/MW/SW*Zl25
     &      *Zn32i/cb)
      AAABR(592) = EE*Sqrt2*(Zl14/CW*Zn12r+Zl14/SW*Zn22r
     &      -Me/MW/SW*Zl44*Zn32r/cb)
      AAABR(593) = EE*Sqrt2*(Zl14/CW*Zn12i+Zl14/SW*Zn22i
     &      -Me/MW/SW*Zl44*Zn32i/cb)
      AAABR(594) = EE*Sqrt2*(2D0*Zl44/CW*Zn12r+Me/MW/SW*Zl14
     &      *Zn32r/cb)
      AAABR(595) = EE*Sqrt2*(2D0*Zl44/CW*Zn12i+Me/MW/SW*Zl14
     &      *Zn32i/cb)
      AAABR(596) = EE*Sqrt2*(Zl36i/CW*Zn12i+Zl36r/CW*Zn12r
     &      +Zl36i/SW*Zn22i+Zl36r/SW*Zn22r-Ml/MW/SW*Zl66i
     &      *Zn32i/cb-Ml/MW/SW*Zl66r*Zn32r/cb)
      AAABR(597) = EE*Sqrt2*(Zl36i/CW*Zn12r-Zl36r/CW*Zn12i
     &      +Zl36i/SW*Zn22r-Zl36r/SW*Zn22i-Ml/MW/SW*Zl66i
     &      *Zn32r/cb+Ml/MW/SW*Zl66r*Zn32i/cb)
      AAABR(598) = EE*Sqrt2*(2D0*Zl66i/CW*Zn12i-2D0*Zl66r
     &      /CW*Zn12r+Ml/MW/SW*Zl36i*Zn32i/cb-Ml/MW/SW*Zl36r
     &      *Zn32r/cb)
      AAABR(599) = EE*Sqrt2*(2D0*Zl66i/CW*Zn12r+2D0*Zl66r
     &      /CW*Zn12i+Ml/MW/SW*Zl36i*Zn32r/cb+Ml/MW/SW*Zl36r
     &      *Zn32i/cb)
      AAABR(600) = EE*Sqrt2*(Zl14/CW*Zn11r+Zl14/SW*Zn21r
     &      -Me/MW/SW*Zl44*Zn31r/cb)
      AAABR(601) = EE*Sqrt2*(Zl14/CW*Zn11i+Zl14/SW*Zn21i
     &      -Me/MW/SW*Zl44*Zn31i/cb)
      AAABR(602) = EE*Sqrt2*(2D0*Zl44/CW*Zn11r+Me/MW/SW*Zl14
     &      *Zn31r/cb)
      AAABR(603) = EE*Sqrt2*(2D0*Zl44/CW*Zn11i+Me/MW/SW*Zl14
     &      *Zn31i/cb)
      AAABR(604) = EE*Sqrt2*(Zl36i/CW*Zn11i+Zl36r/CW*Zn11r
     &      +Zl36i/SW*Zn21i+Zl36r/SW*Zn21r-Ml/MW/SW*Zl66i
     &      *Zn31i/cb-Ml/MW/SW*Zl66r*Zn31r/cb)
      AAABR(605) = EE*Sqrt2*(Zl36i/CW*Zn11r-Zl36r/CW*Zn11i
     &      +Zl36i/SW*Zn21r-Zl36r/SW*Zn21i-Ml/MW/SW*Zl66i
     &      *Zn31r/cb+Ml/MW/SW*Zl66r*Zn31i/cb)
      AAABR(606) = EE*Sqrt2*(2D0*Zl66i/CW*Zn11i-2D0*Zl66r
     &      /CW*Zn11r+Ml/MW/SW*Zl36i*Zn31i/cb-Ml/MW/SW*Zl36r
     &      *Zn31r/cb)
      AAABR(607) = EE*Sqrt2*(2D0*Zl66i/CW*Zn11r+2D0*Zl66r
     &      /CW*Zn11i+Ml/MW/SW*Zl36i*Zn31r/cb+Ml/MW/SW*Zl36r
     &      *Zn31i/cb)
      AAABR(608) = EE/SW*Zp11r
      AAABR(609) = EE/SW*Zp11i
      AAABR(610) = EE/MW*Me/SW*Sqrt2*Zm21r/cb
      AAABR(611) = EE/MW*Me/SW*Sqrt2*Zm21i/cb
      AAABR(612) = EE/MW*Ml/SW*Sqrt2*Zm21r/cb
      AAABR(613) = EE/MW*Ml/SW*Sqrt2*Zm21i/cb
      AAABR(614) = EE/MW*Mm/SW*Sqrt2*Zm21r/cb
      AAABR(615) = EE/MW*Mm/SW*Sqrt2*Zm21i/cb
      AAABR(616) = EE/SW*Zp12r
      AAABR(617) = EE/SW*Zp12i
      AAABR(618) = EE/MW*Ml/SW*Sqrt2*Zm22r/cb
      AAABR(619) = EE/MW*Ml/SW*Sqrt2*Zm22i/cb
      AAABR(620) = EE/MW*Me/SW*Sqrt2*Zm22r/cb
      AAABR(621) = EE/MW*Me/SW*Sqrt2*Zm22i/cb
      AAABR(622) = EE/MW*Mm/SW*Sqrt2*Zm22r/cb
      AAABR(623) = EE/MW*Mm/SW*Sqrt2*Zm22i/cb
      AAABR(624) = EE*Sqrt2*(Zn14r/CW-Zn24r/SW)
      AAABR(625) = EE*Sqrt2*(Zn14i/CW-Zn24i/SW)
      AAABR(626) = EE*Sqrt2*(Zn13r/CW-Zn23r/SW)
      AAABR(627) = EE*Sqrt2*(Zn13i/CW-Zn23i/SW)
      AAABR(628) = EE*Sqrt2*(Zn12r/CW-Zn22r/SW)
      AAABR(629) = EE*Sqrt2*(Zn12i/CW-Zn22i/SW)
      AAABR(630) = EE*Sqrt2*(Zn11r/CW-Zn21r/SW)
      AAABR(631) = EE*Sqrt2*(Zn11i/CW-Zn21i/SW)
      AAABR(632) = EE*Sqrt2*(Zd33i/CW*Zn11i+Zd33r/CW*Zn11r
     &      -3D0*Zd33i/SW*Zn21i-3D0*Zd33r/SW*Zn21r+3D0*Mb
     &      /MW/SW*Zd63i*Zn31i/cb+3D0*Mb/MW/SW*Zd63r*Zn31r
     &      /cb)
      AAABR(633) = EE*Sqrt2*(Zd33i/CW*Zn11r-Zd33r/CW*Zn11i
     &      -3D0*Zd33i/SW*Zn21r+3D0*Zd33r/SW*Zn21i+3D0*Mb
     &      /MW/SW*Zd63i*Zn31r/cb-3D0*Mb/MW/SW*Zd63r*Zn31i
     &      /cb)
      AAABR(634) = EE*Sqrt2*(2D0*Zd63i/CW*Zn11i-2D0*Zd63r
     &      /CW*Zn11r+3D0*Mb/MW/SW*Zd33i*Zn31i/cb-3D0*Mb
     &      /MW/SW*Zd33r*Zn31r/cb)
      AAABR(635) = EE*Sqrt2*(2D0*Zd63i/CW*Zn11r+2D0*Zd63r
     &      /CW*Zn11i+3D0*Mb/MW/SW*Zd33i*Zn31r/cb+3D0*Mb
     &      /MW/SW*Zd33r*Zn31i/cb)
      AAABR(636) = EE*Sqrt2*Zd11*(Zn12r/CW-3D0*Zn22r/SW)
      AAABR(637) = EE*Sqrt2*Zd11*(Zn12i/CW-3D0*Zn22i/SW)
      AAABR(638) = EE/CW*Sqrt2*Zd41*Zn12r
      AAABR(639) = EE/CW*Sqrt2*Zd41*Zn12i
      AAABR(640) = EE*Sqrt2*Zd11*(Zn13r/CW-3D0*Zn23r/SW)
      AAABR(641) = EE*Sqrt2*Zd11*(Zn13i/CW-3D0*Zn23i/SW)
      AAABR(642) = EE/CW*Sqrt2*Zd41*Zn13r
      AAABR(643) = EE/CW*Sqrt2*Zd41*Zn13i
      AAABR(644) = EE*Sqrt2*(Zd33i/CW*Zn14i+Zd33r/CW*Zn14r
     &      -3D0*Zd33i/SW*Zn24i-3D0*Zd33r/SW*Zn24r+3D0*Mb
     &      /MW/SW*Zd63i*Zn34i/cb+3D0*Mb/MW/SW*Zd63r*Zn34r
     &      /cb)
      AAABR(645) = EE*Sqrt2*(Zd33i/CW*Zn14r-Zd33r/CW*Zn14i
     &      -3D0*Zd33i/SW*Zn24r+3D0*Zd33r/SW*Zn24i+3D0*Mb
     &      /MW/SW*Zd63i*Zn34r/cb-3D0*Mb/MW/SW*Zd63r*Zn34i
     &      /cb)
      AAABR(646) = EE*Sqrt2*(2D0*Zd63i/CW*Zn14i-2D0*Zd63r
     &      /CW*Zn14r+3D0*Mb/MW/SW*Zd33i*Zn34i/cb-3D0*Mb
     &      /MW/SW*Zd33r*Zn34r/cb)
      AAABR(647) = EE*Sqrt2*(2D0*Zd63i/CW*Zn14r+2D0*Zd63r
     &      /CW*Zn14i+3D0*Mb/MW/SW*Zd33i*Zn34r/cb+3D0*Mb
     &      /MW/SW*Zd33r*Zn34i/cb)
      AAABR(648) = EE*Sqrt2*Zd11*(Zn11r/CW-3D0*Zn21r/SW)
      AAABR(649) = EE*Sqrt2*Zd11*(Zn11i/CW-3D0*Zn21i/SW)
      AAABR(650) = EE/CW*Sqrt2*Zd41*Zn11r
      AAABR(651) = EE/CW*Sqrt2*Zd41*Zn11i
      AAABR(652) = EE*Sqrt2*Zd11*(Zn14r/CW-3D0*Zn24r/SW)
      AAABR(653) = EE*Sqrt2*Zd11*(Zn14i/CW-3D0*Zn24i/SW)
      AAABR(654) = EE/CW*Sqrt2*Zd41*Zn14r
      AAABR(655) = EE/CW*Sqrt2*Zd41*Zn14i
      AAABR(656) = EE*Sqrt2*(Zd22/CW*Zn13r-3D0*Zd22/SW*Zn23r
     &      +3D0*Ms/MW/SW*Zd52*Zn33r/cb)
      AAABR(657) = EE*Sqrt2*(Zd22/CW*Zn13i-3D0*Zd22/SW*Zn23i
     &      +3D0*Ms/MW/SW*Zd52*Zn33i/cb)
      AAABR(658) = EE*Sqrt2*(2D0*Zd52/CW*Zn13r+3D0*Ms/MW
     &      /SW*Zd22*Zn33r/cb)
      AAABR(659) = EE*Sqrt2*(2D0*Zd52/CW*Zn13i+3D0*Ms/MW
     &      /SW*Zd22*Zn33i/cb)
      AAABR(660) = EE*Sqrt2*(Zd33i/CW*Zn13i+Zd33r/CW*Zn13r
     &      -3D0*Zd33i/SW*Zn23i-3D0*Zd33r/SW*Zn23r+3D0*Mb
     &      /MW/SW*Zd63i*Zn33i/cb+3D0*Mb/MW/SW*Zd63r*Zn33r
     &      /cb)
      AAABR(661) = EE*Sqrt2*(Zd33i/CW*Zn13r-Zd33r/CW*Zn13i
     &      -3D0*Zd33i/SW*Zn23r+3D0*Zd33r/SW*Zn23i+3D0*Mb
     &      /MW/SW*Zd63i*Zn33r/cb-3D0*Mb/MW/SW*Zd63r*Zn33i
     &      /cb)
      AAABR(662) = EE*Sqrt2*(2D0*Zd63i/CW*Zn13i-2D0*Zd63r
     &      /CW*Zn13r+3D0*Mb/MW/SW*Zd33i*Zn33i/cb-3D0*Mb
     &      /MW/SW*Zd33r*Zn33r/cb)
      AAABR(663) = EE*Sqrt2*(2D0*Zd63i/CW*Zn13r+2D0*Zd63r
     &      /CW*Zn13i+3D0*Mb/MW/SW*Zd33i*Zn33r/cb+3D0*Mb
     &      /MW/SW*Zd33r*Zn33i/cb)
      AAABR(664) = EE*Sqrt2*(Zd22/CW*Zn11r-3D0*Zd22/SW*Zn21r
     &      +3D0*Ms/MW/SW*Zd52*Zn31r/cb)
      AAABR(665) = EE*Sqrt2*(Zd22/CW*Zn11i-3D0*Zd22/SW*Zn21i
     &      +3D0*Ms/MW/SW*Zd52*Zn31i/cb)
      AAABR(666) = EE*Sqrt2*(2D0*Zd52/CW*Zn11r+3D0*Ms/MW
     &      /SW*Zd22*Zn31r/cb)
      AAABR(667) = EE*Sqrt2*(2D0*Zd52/CW*Zn11i+3D0*Ms/MW
     &      /SW*Zd22*Zn31i/cb)
      AAABR(668) = EE*Sqrt2*(Zd33i/CW*Zn12i+Zd33r/CW*Zn12r
     &      -3D0*Zd33i/SW*Zn22i-3D0*Zd33r/SW*Zn22r+3D0*Mb
     &      /MW/SW*Zd63i*Zn32i/cb+3D0*Mb/MW/SW*Zd63r*Zn32r
     &      /cb)
      AAABR(669) = EE*Sqrt2*(Zd33i/CW*Zn12r-Zd33r/CW*Zn12i
     &      -3D0*Zd33i/SW*Zn22r+3D0*Zd33r/SW*Zn22i+3D0*Mb
     &      /MW/SW*Zd63i*Zn32r/cb-3D0*Mb/MW/SW*Zd63r*Zn32i
     &      /cb)
      AAABR(670) = EE*Sqrt2*(2D0*Zd63i/CW*Zn12i-2D0*Zd63r
     &      /CW*Zn12r+3D0*Mb/MW/SW*Zd33i*Zn32i/cb-3D0*Mb
     &      /MW/SW*Zd33r*Zn32r/cb)
      AAABR(671) = EE*Sqrt2*(2D0*Zd63i/CW*Zn12r+2D0*Zd63r
     &      /CW*Zn12i+3D0*Mb/MW/SW*Zd33i*Zn32r/cb+3D0*Mb
     &      /MW/SW*Zd33r*Zn32i/cb)
      AAABR(672) = EE*Sqrt2*(Zd22/CW*Zn14r-3D0*Zd22/SW*Zn24r
     &      +3D0*Ms/MW/SW*Zd52*Zn34r/cb)
      AAABR(673) = EE*Sqrt2*(Zd22/CW*Zn14i-3D0*Zd22/SW*Zn24i
     &      +3D0*Ms/MW/SW*Zd52*Zn34i/cb)
      AAABR(674) = EE*Sqrt2*(2D0*Zd52/CW*Zn14r+3D0*Ms/MW
     &      /SW*Zd22*Zn34r/cb)
      AAABR(675) = EE*Sqrt2*(2D0*Zd52/CW*Zn14i+3D0*Ms/MW
     &      /SW*Zd22*Zn34i/cb)
      AAABR(676) = EE*Sqrt2*(Zd22/CW*Zn12r-3D0*Zd22/SW*Zn22r
     &      +3D0*Ms/MW/SW*Zd52*Zn32r/cb)
      AAABR(677) = EE*Sqrt2*(Zd22/CW*Zn12i-3D0*Zd22/SW*Zn22i
     &      +3D0*Ms/MW/SW*Zd52*Zn32i/cb)
      AAABR(678) = EE*Sqrt2*(2D0*Zd52/CW*Zn12r+3D0*Ms/MW
     &      /SW*Zd22*Zn32r/cb)
      AAABR(679) = EE*Sqrt2*(2D0*Zd52/CW*Zn12i+3D0*Ms/MW
     &      /SW*Zd22*Zn32i/cb)
      AAABR(680) = EE*Sqrt2*(Zd25/CW*Zn12r-3D0*Zd25/SW*Zn22r
     &      +3D0*Ms/MW/SW*Zd55*Zn32r/cb)
      AAABR(681) = EE*Sqrt2*(Zd25/CW*Zn12i-3D0*Zd25/SW*Zn22i
     &      +3D0*Ms/MW/SW*Zd55*Zn32i/cb)
      AAABR(682) = EE*Sqrt2*(2D0*Zd55/CW*Zn12r+3D0*Ms/MW
     &      /SW*Zd25*Zn32r/cb)
      AAABR(683) = EE*Sqrt2*(2D0*Zd55/CW*Zn12i+3D0*Ms/MW
     &      /SW*Zd25*Zn32i/cb)
      AAABR(684) = EE*Sqrt2*(Zd36i/CW*Zn14i+Zd36r/CW*Zn14r
     &      -3D0*Zd36i/SW*Zn24i-3D0*Zd36r/SW*Zn24r+3D0*Mb
     &      /MW/SW*Zd66i*Zn34i/cb+3D0*Mb/MW/SW*Zd66r*Zn34r
     &      /cb)
      AAABR(685) = EE*Sqrt2*(Zd36i/CW*Zn14r-Zd36r/CW*Zn14i
     &      -3D0*Zd36i/SW*Zn24r+3D0*Zd36r/SW*Zn24i+3D0*Mb
     &      /MW/SW*Zd66i*Zn34r/cb-3D0*Mb/MW/SW*Zd66r*Zn34i
     &      /cb)
      AAABR(686) = EE*Sqrt2*(2D0*Zd66i/CW*Zn14i-2D0*Zd66r
     &      /CW*Zn14r+3D0*Mb/MW/SW*Zd36i*Zn34i/cb-3D0*Mb
     &      /MW/SW*Zd36r*Zn34r/cb)
      AAABR(687) = EE*Sqrt2*(2D0*Zd66i/CW*Zn14r+2D0*Zd66r
     &      /CW*Zn14i+3D0*Mb/MW/SW*Zd36i*Zn34r/cb+3D0*Mb
     &      /MW/SW*Zd36r*Zn34i/cb)
      AAABR(688) = EE*Sqrt2*Zd14*(Zn13r/CW-3D0*Zn23r/SW)
      AAABR(689) = EE*Sqrt2*Zd14*(Zn13i/CW-3D0*Zn23i/SW)
      AAABR(690) = EE/CW*Sqrt2*Zd44*Zn13r
      AAABR(691) = EE/CW*Sqrt2*Zd44*Zn13i
      AAABR(692) = EE*Sqrt2*(Zd25/CW*Zn14r-3D0*Zd25/SW*Zn24r
     &      +3D0*Ms/MW/SW*Zd55*Zn34r/cb)
      AAABR(693) = EE*Sqrt2*(Zd25/CW*Zn14i-3D0*Zd25/SW*Zn24i
     &      +3D0*Ms/MW/SW*Zd55*Zn34i/cb)
      AAABR(694) = EE*Sqrt2*(2D0*Zd55/CW*Zn14r+3D0*Ms/MW
     &      /SW*Zd25*Zn34r/cb)
      AAABR(695) = EE*Sqrt2*(2D0*Zd55/CW*Zn14i+3D0*Ms/MW
     &      /SW*Zd25*Zn34i/cb)
      AAABR(696) = EE*Sqrt2*(Zd36i/CW*Zn11i+Zd36r/CW*Zn11r
     &      -3D0*Zd36i/SW*Zn21i-3D0*Zd36r/SW*Zn21r+3D0*Mb
     &      /MW/SW*Zd66i*Zn31i/cb+3D0*Mb/MW/SW*Zd66r*Zn31r
     &      /cb)
      AAABR(697) = EE*Sqrt2*(Zd36i/CW*Zn11r-Zd36r/CW*Zn11i
     &      -3D0*Zd36i/SW*Zn21r+3D0*Zd36r/SW*Zn21i+3D0*Mb
     &      /MW/SW*Zd66i*Zn31r/cb-3D0*Mb/MW/SW*Zd66r*Zn31i
     &      /cb)
      AAABR(698) = EE*Sqrt2*(2D0*Zd66i/CW*Zn11i-2D0*Zd66r
     &      /CW*Zn11r+3D0*Mb/MW/SW*Zd36i*Zn31i/cb-3D0*Mb
     &      /MW/SW*Zd36r*Zn31r/cb)
      AAABR(699) = EE*Sqrt2*(2D0*Zd66i/CW*Zn11r+2D0*Zd66r
     &      /CW*Zn11i+3D0*Mb/MW/SW*Zd36i*Zn31r/cb+3D0*Mb
     &      /MW/SW*Zd36r*Zn31i/cb)
      AAABR(700) = EE*Sqrt2*(Zd36i/CW*Zn12i+Zd36r/CW*Zn12r
     &      -3D0*Zd36i/SW*Zn22i-3D0*Zd36r/SW*Zn22r+3D0*Mb
     &      /MW/SW*Zd66i*Zn32i/cb+3D0*Mb/MW/SW*Zd66r*Zn32r
     &      /cb)
      AAABR(701) = EE*Sqrt2*(Zd36i/CW*Zn12r-Zd36r/CW*Zn12i
     &      -3D0*Zd36i/SW*Zn22r+3D0*Zd36r/SW*Zn22i+3D0*Mb
     &      /MW/SW*Zd66i*Zn32r/cb-3D0*Mb/MW/SW*Zd66r*Zn32i
     &      /cb)
      AAABR(702) = EE*Sqrt2*(2D0*Zd66i/CW*Zn12i-2D0*Zd66r
     &      /CW*Zn12r+3D0*Mb/MW/SW*Zd36i*Zn32i/cb-3D0*Mb
     &      /MW/SW*Zd36r*Zn32r/cb)
      AAABR(703) = EE*Sqrt2*(2D0*Zd66i/CW*Zn12r+2D0*Zd66r
     &      /CW*Zn12i+3D0*Mb/MW/SW*Zd36i*Zn32r/cb+3D0*Mb
     &      /MW/SW*Zd36r*Zn32i/cb)
      AAABR(704) = EE*Sqrt2*Zd14*(Zn11r/CW-3D0*Zn21r/SW)
      AAABR(705) = EE*Sqrt2*Zd14*(Zn11i/CW-3D0*Zn21i/SW)
      AAABR(706) = EE/CW*Sqrt2*Zd44*Zn11r
      AAABR(707) = EE/CW*Sqrt2*Zd44*Zn11i
      AAABR(708) = EE*Sqrt2*Zd14*(Zn12r/CW-3D0*Zn22r/SW)
      AAABR(709) = EE*Sqrt2*Zd14*(Zn12i/CW-3D0*Zn22i/SW)
      AAABR(710) = EE/CW*Sqrt2*Zd44*Zn12r
      AAABR(711) = EE/CW*Sqrt2*Zd44*Zn12i
      AAABR(712) = EE*Sqrt2*(Zd25/CW*Zn11r-3D0*Zd25/SW*Zn21r
     &      +3D0*Ms/MW/SW*Zd55*Zn31r/cb)
      AAABR(713) = EE*Sqrt2*(Zd25/CW*Zn11i-3D0*Zd25/SW*Zn21i
     &      +3D0*Ms/MW/SW*Zd55*Zn31i/cb)
      AAABR(714) = EE*Sqrt2*(2D0*Zd55/CW*Zn11r+3D0*Ms/MW
     &      /SW*Zd25*Zn31r/cb)
      AAABR(715) = EE*Sqrt2*(2D0*Zd55/CW*Zn11i+3D0*Ms/MW
     &      /SW*Zd25*Zn31i/cb)
      AAABR(716) = EE*Sqrt2*(Zd25/CW*Zn13r-3D0*Zd25/SW*Zn23r
     &      +3D0*Ms/MW/SW*Zd55*Zn33r/cb)
      AAABR(717) = EE*Sqrt2*(Zd25/CW*Zn13i-3D0*Zd25/SW*Zn23i
     &      +3D0*Ms/MW/SW*Zd55*Zn33i/cb)
      AAABR(718) = EE*Sqrt2*(2D0*Zd55/CW*Zn13r+3D0*Ms/MW
     &      /SW*Zd25*Zn33r/cb)
      AAABR(719) = EE*Sqrt2*(2D0*Zd55/CW*Zn13i+3D0*Ms/MW
     &      /SW*Zd25*Zn33i/cb)
      AAABR(720) = EE*Sqrt2*(Zd36i/CW*Zn13i+Zd36r/CW*Zn13r
     &      -3D0*Zd36i/SW*Zn23i-3D0*Zd36r/SW*Zn23r+3D0*Mb
     &      /MW/SW*Zd66i*Zn33i/cb+3D0*Mb/MW/SW*Zd66r*Zn33r
     &      /cb)
      AAABR(721) = EE*Sqrt2*(Zd36i/CW*Zn13r-Zd36r/CW*Zn13i
     &      -3D0*Zd36i/SW*Zn23r+3D0*Zd36r/SW*Zn23i+3D0*Mb
     &      /MW/SW*Zd66i*Zn33r/cb-3D0*Mb/MW/SW*Zd66r*Zn33i
     &      /cb)
      AAABR(722) = EE*Sqrt2*(2D0*Zd66i/CW*Zn13i-2D0*Zd66r
     &      /CW*Zn13r+3D0*Mb/MW/SW*Zd36i*Zn33i/cb-3D0*Mb
     &      /MW/SW*Zd36r*Zn33r/cb)
      AAABR(723) = EE*Sqrt2*(2D0*Zd66i/CW*Zn13r+2D0*Zd66r
     &      /CW*Zn13i+3D0*Mb/MW/SW*Zd36i*Zn33r/cb+3D0*Mb
     &      /MW/SW*Zd36r*Zn33i/cb)
      AAABR(724) = EE*Sqrt2*Zd14*(Zn14r/CW-3D0*Zn24r/SW)
      AAABR(725) = EE*Sqrt2*Zd14*(Zn14i/CW-3D0*Zn24i/SW)
      AAABR(726) = EE/CW*Sqrt2*Zd44*Zn14r
      AAABR(727) = EE/CW*Sqrt2*Zd44*Zn14i
      AAABR(728) = EE/MW*Ms/SW
      AAABR(729) = EE/MW*Mb/SW
      AAABR(730) = EE/MW*Ms/SW*Zh11/cb
      AAABR(731) = EE/MW*Ms/SW*Zh31/cb*sb
      AAABR(732) = EE/MW*Ms/SW*Zh13/cb
      AAABR(733) = EE/MW*Ms/SW*Zh33/cb*sb
      AAABR(734) = EE/MW*Ms/SW*Zh12/cb
      AAABR(735) = EE/MW*Ms/SW*Zh32/cb*sb
      AAABR(736) = EE/MW*Mb/SW*Zh11/cb
      AAABR(737) = EE/MW*Mb/SW*Zh31/cb*sb
      AAABR(738) = EE/MW*Mb/SW*Zh13/cb
      AAABR(739) = EE/MW*Mb/SW*Zh33/cb*sb
      AAABR(740) = EE/MW*Mb/SW*Zh12/cb
      AAABR(741) = EE/MW*Mb/SW*Zh32/cb*sb
      AAABR(742) = EE/MW*Mc/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(743) = EE/MW*Ms/SW*Sqrt2/cb*sb
      AAABR(744) = EE/MW*Mt/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(745) = EE/MW*Mb/SW*Sqrt2/cb*sb
      AAABR(746) = EE/MW*Mc/SW*Sqrt2
      AAABR(747) = EE/MW*Ms/SW*Sqrt2
      AAABR(748) = EE/MW*Mt/SW*Sqrt2
      AAABR(749) = EE/MW*Mb/SW*Sqrt2
      AAABR(750) = EE/SW*(2D0*Zp11r*Zu22-Mc/MW*Sqrt2*Zp21r
     &      *Zu52/sb)
      AAABR(751) = EE/SW*(2D0*Zp11i*Zu22-Mc/MW*Sqrt2*Zp21i
     &      *Zu52/sb)
      AAABR(752) = EE/MW*Ms/SW*Sqrt2*Zm21r*Zu22/cb
      AAABR(753) = EE/MW*Ms/SW*Sqrt2*Zm21i*Zu22/cb
      AAABR(754) = EE/SW*(2D0*Zp11i*Zu33i+2D0*Zp11r*Zu33r
     &      -Mt/MW*Sqrt2*Zp21i*Zu63i/sb-Mt/MW*Sqrt2*Zp21r
     &      *Zu63r/sb)
      AAABR(755) = EE/SW*(2D0*Zp11i*Zu33r-2D0*Zp11r*Zu33i
     &      -Mt/MW*Sqrt2*Zp21i*Zu63r/sb+Mt/MW*Sqrt2*Zp21r
     &      *Zu63i/sb)
      AAABR(756) = EE/MW*Mb/SW*Sqrt2/cb*(Zm21i*Zu33i+Zm21r
     &      *Zu33r)
      AAABR(757) = EE/MW*Mb/SW*Sqrt2/cb*(Zm21i*Zu33r-Zm21r
     &      *Zu33i)
      AAABR(758) = EE/SW*Zp11r*Zu11
      AAABR(759) = EE/SW*Zp11i*Zu11
      AAABR(760) = EE/SW*(2D0*Zp11i*Zu36i+2D0*Zp11r*Zu36r
     &      -Mt/MW*Sqrt2*Zp21i*Zu66i/sb-Mt/MW*Sqrt2*Zp21r
     &      *Zu66r/sb)
      AAABR(761) = EE/SW*(2D0*Zp11i*Zu36r-2D0*Zp11r*Zu36i
     &      -Mt/MW*Sqrt2*Zp21i*Zu66r/sb+Mt/MW*Sqrt2*Zp21r
     &      *Zu66i/sb)
      AAABR(762) = EE/MW*Mb/SW*Sqrt2/cb*(Zm21i*Zu36i+Zm21r
     &      *Zu36r)
      AAABR(763) = EE/MW*Mb/SW*Sqrt2/cb*(Zm21i*Zu36r-Zm21r
     &      *Zu36i)
      AAABR(764) = EE/SW*(2D0*Zp11r*Zu25-Mc/MW*Sqrt2*Zp21r
     &      *Zu55/sb)
      AAABR(765) = EE/SW*(2D0*Zp11i*Zu25-Mc/MW*Sqrt2*Zp21i
     &      *Zu55/sb)
      AAABR(766) = EE/MW*Ms/SW*Sqrt2*Zm21r*Zu25/cb
      AAABR(767) = EE/MW*Ms/SW*Sqrt2*Zm21i*Zu25/cb
      AAABR(768) = EE/SW*Zp11r*Zu14
      AAABR(769) = EE/SW*Zp11i*Zu14
      AAABR(770) = EE/SW*(2D0*Zp12i*Zu33i+2D0*Zp12r*Zu33r
     &      -Mt/MW*Sqrt2*Zp22i*Zu63i/sb-Mt/MW*Sqrt2*Zp22r
     &      *Zu63r/sb)
      AAABR(771) = EE/SW*(2D0*Zp12i*Zu33r-2D0*Zp12r*Zu33i
     &      -Mt/MW*Sqrt2*Zp22i*Zu63r/sb+Mt/MW*Sqrt2*Zp22r
     &      *Zu63i/sb)
      AAABR(772) = EE/MW*Mb/SW*Sqrt2/cb*(Zm22i*Zu33i+Zm22r
     &      *Zu33r)
      AAABR(773) = EE/MW*Mb/SW*Sqrt2/cb*(Zm22i*Zu33r-Zm22r
     &      *Zu33i)
      AAABR(774) = EE/SW*Zp12r*Zu11
      AAABR(775) = EE/SW*Zp12i*Zu11
      AAABR(776) = EE/SW*(2D0*Zp12r*Zu22-Mc/MW*Sqrt2*Zp22r
     &      *Zu52/sb)
      AAABR(777) = EE/SW*(2D0*Zp12i*Zu22-Mc/MW*Sqrt2*Zp22i
     &      *Zu52/sb)
      AAABR(778) = EE/MW*Ms/SW*Sqrt2*Zm22r*Zu22/cb
      AAABR(779) = EE/MW*Ms/SW*Sqrt2*Zm22i*Zu22/cb
      AAABR(780) = EE/SW*(2D0*Zp12i*Zu36i+2D0*Zp12r*Zu36r
     &      -Mt/MW*Sqrt2*Zp22i*Zu66i/sb-Mt/MW*Sqrt2*Zp22r
     &      *Zu66r/sb)
      AAABR(781) = EE/SW*(2D0*Zp12i*Zu36r-2D0*Zp12r*Zu36i
     &      -Mt/MW*Sqrt2*Zp22i*Zu66r/sb+Mt/MW*Sqrt2*Zp22r
     &      *Zu66i/sb)
      AAABR(782) = EE/MW*Mb/SW*Sqrt2/cb*(Zm22i*Zu36i+Zm22r
     &      *Zu36r)
      AAABR(783) = EE/MW*Mb/SW*Sqrt2/cb*(Zm22i*Zu36r-Zm22r
     &      *Zu36i)
      AAABR(784) = EE/SW*(2D0*Zp12r*Zu25-Mc/MW*Sqrt2*Zp22r
     &      *Zu55/sb)
      AAABR(785) = EE/SW*(2D0*Zp12i*Zu25-Mc/MW*Sqrt2*Zp22i
     &      *Zu55/sb)
      AAABR(786) = EE/MW*Ms/SW*Sqrt2*Zm22r*Zu25/cb
      AAABR(787) = EE/MW*Ms/SW*Sqrt2*Zm22i*Zu25/cb
      AAABR(788) = EE/SW*Zp12r*Zu14
      AAABR(789) = EE/SW*Zp12i*Zu14
      AAABR(790) = GG*Sqrt2*Zd11
      AAABR(791) = GG*Sqrt2*Zd41
      AAABR(792) = GG*Sqrt2*Zd22
      AAABR(793) = GG*Sqrt2*Zd52
      AAABR(794) = GG*Sqrt2*Zd33r
      AAABR(795) = GG*Sqrt2*Zd33i
      AAABR(796) = GG*Sqrt2*Zd63r
      AAABR(797) = GG*Sqrt2*Zd63i
      AAABR(798) = GG*Sqrt2*Zd14
      AAABR(799) = GG*Sqrt2*Zd44
      AAABR(800) = GG*Sqrt2*Zd36r
      AAABR(801) = GG*Sqrt2*Zd36i
      AAABR(802) = GG*Sqrt2*Zd66r
      AAABR(803) = GG*Sqrt2*Zd66i
      AAABR(804) = GG*Sqrt2*Zd25
      AAABR(805) = GG*Sqrt2*Zd55
      AAABR(806) = EE*Sqrt2*(Zn11i/CW*Zu33i+Zn11r/CW*Zu33r
     &      +3D0*Zn21i/SW*Zu33i+3D0*Zn21r/SW*Zu33r+3D0*Mt
     &      /MW/SW*Zn41i*Zu63i/sb+3D0*Mt/MW/SW*Zn41r*Zu63r
     &      /sb)
      AAABR(807) = EE*Sqrt2*(Zn11i/CW*Zu33r-Zn11r/CW*Zu33i
     &      +3D0*Zn21i/SW*Zu33r-3D0*Zn21r/SW*Zu33i+3D0*Mt
     &      /MW/SW*Zn41i*Zu63r/sb-3D0*Mt/MW/SW*Zn41r*Zu63i
     &      /sb)
      AAABR(808) = EE*Sqrt2*(4D0*Zn11i/CW*Zu63i-4D0*Zn11r
     &      /CW*Zu63r-3D0*Mt/MW/SW*Zn41i*Zu33i/sb+3D0*Mt
     &      /MW/SW*Zn41r*Zu33r/sb)
      AAABR(809) = EE*Sqrt2*(4D0*Zn11i/CW*Zu63r+4D0*Zn11r
     &      /CW*Zu63i-3D0*Mt/MW/SW*Zn41i*Zu33r/sb-3D0*Mt
     &      /MW/SW*Zn41r*Zu33i/sb)
      AAABR(810) = EE*Sqrt2*(Zn11r/CW*Zu22+3D0*Zn21r/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn41r*Zu52/sb)
      AAABR(811) = EE*Sqrt2*(Zn11i/CW*Zu22+3D0*Zn21i/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn41i*Zu52/sb)
      AAABR(812) = EE*Sqrt2*(4D0*Zn11r/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn41r*Zu22/sb)
      AAABR(813) = EE*Sqrt2*(4D0*Zn11i/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn41i*Zu22/sb)
      AAABR(814) = EE*Sqrt2*(Zn12i/CW*Zu33i+Zn12r/CW*Zu33r
     &      +3D0*Zn22i/SW*Zu33i+3D0*Zn22r/SW*Zu33r+3D0*Mt
     &      /MW/SW*Zn42i*Zu63i/sb+3D0*Mt/MW/SW*Zn42r*Zu63r
     &      /sb)
      AAABR(815) = EE*Sqrt2*(Zn12i/CW*Zu33r-Zn12r/CW*Zu33i
     &      +3D0*Zn22i/SW*Zu33r-3D0*Zn22r/SW*Zu33i+3D0*Mt
     &      /MW/SW*Zn42i*Zu63r/sb-3D0*Mt/MW/SW*Zn42r*Zu63i
     &      /sb)
      AAABR(816) = EE*Sqrt2*(4D0*Zn12i/CW*Zu63i-4D0*Zn12r
     &      /CW*Zu63r-3D0*Mt/MW/SW*Zn42i*Zu33i/sb+3D0*Mt
     &      /MW/SW*Zn42r*Zu33r/sb)
      AAABR(817) = EE*Sqrt2*(4D0*Zn12i/CW*Zu63r+4D0*Zn12r
     &      /CW*Zu63i-3D0*Mt/MW/SW*Zn42i*Zu33r/sb-3D0*Mt
     &      /MW/SW*Zn42r*Zu33i/sb)
      AAABR(818) = EE*Sqrt2*Zu11*(Zn14r/CW+3D0*Zn24r/SW)
      AAABR(819) = EE*Sqrt2*Zu11*(Zn14i/CW+3D0*Zn24i/SW)
      AAABR(820) = EE/CW*Sqrt2*Zn14r*Zu41
      AAABR(821) = EE/CW*Sqrt2*Zn14i*Zu41
      AAABR(822) = EE*Sqrt2*Zu11*(Zn13r/CW+3D0*Zn23r/SW)
      AAABR(823) = EE*Sqrt2*Zu11*(Zn13i/CW+3D0*Zn23i/SW)
      AAABR(824) = EE/CW*Sqrt2*Zn13r*Zu41
      AAABR(825) = EE/CW*Sqrt2*Zn13i*Zu41
      AAABR(826) = EE*Sqrt2*Zu11*(Zn11r/CW+3D0*Zn21r/SW)
      AAABR(827) = EE*Sqrt2*Zu11*(Zn11i/CW+3D0*Zn21i/SW)
      AAABR(828) = EE/CW*Sqrt2*Zn11r*Zu41
      AAABR(829) = EE/CW*Sqrt2*Zn11i*Zu41
      AAABR(830) = EE*Sqrt2*(Zn14r/CW*Zu22+3D0*Zn24r/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn44r*Zu52/sb)
      AAABR(831) = EE*Sqrt2*(Zn14i/CW*Zu22+3D0*Zn24i/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn44i*Zu52/sb)
      AAABR(832) = EE*Sqrt2*(4D0*Zn14r/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn44r*Zu22/sb)
      AAABR(833) = EE*Sqrt2*(4D0*Zn14i/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn44i*Zu22/sb)
      AAABR(834) = EE*Sqrt2*Zu11*(Zn12r/CW+3D0*Zn22r/SW)
      AAABR(835) = EE*Sqrt2*Zu11*(Zn12i/CW+3D0*Zn22i/SW)
      AAABR(836) = EE/CW*Sqrt2*Zn12r*Zu41
      AAABR(837) = EE/CW*Sqrt2*Zn12i*Zu41
      AAABR(838) = EE*Sqrt2*(Zn12r/CW*Zu22+3D0*Zn22r/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn42r*Zu52/sb)
      AAABR(839) = EE*Sqrt2*(Zn12i/CW*Zu22+3D0*Zn22i/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn42i*Zu52/sb)
      AAABR(840) = EE*Sqrt2*(4D0*Zn12r/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn42r*Zu22/sb)
      AAABR(841) = EE*Sqrt2*(4D0*Zn12i/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn42i*Zu22/sb)
      AAABR(842) = EE*Sqrt2*(Zn13i/CW*Zu33i+Zn13r/CW*Zu33r
     &      +3D0*Zn23i/SW*Zu33i+3D0*Zn23r/SW*Zu33r+3D0*Mt
     &      /MW/SW*Zn43i*Zu63i/sb+3D0*Mt/MW/SW*Zn43r*Zu63r
     &      /sb)
      AAABR(843) = EE*Sqrt2*(Zn13i/CW*Zu33r-Zn13r/CW*Zu33i
     &      +3D0*Zn23i/SW*Zu33r-3D0*Zn23r/SW*Zu33i+3D0*Mt
     &      /MW/SW*Zn43i*Zu63r/sb-3D0*Mt/MW/SW*Zn43r*Zu63i
     &      /sb)
      AAABR(844) = EE*Sqrt2*(4D0*Zn13i/CW*Zu63i-4D0*Zn13r
     &      /CW*Zu63r-3D0*Mt/MW/SW*Zn43i*Zu33i/sb+3D0*Mt
     &      /MW/SW*Zn43r*Zu33r/sb)
      AAABR(845) = EE*Sqrt2*(4D0*Zn13i/CW*Zu63r+4D0*Zn13r
     &      /CW*Zu63i-3D0*Mt/MW/SW*Zn43i*Zu33r/sb-3D0*Mt
     &      /MW/SW*Zn43r*Zu33i/sb)
      AAABR(846) = EE*Sqrt2*(Zn14i/CW*Zu33i+Zn14r/CW*Zu33r
     &      +3D0*Zn24i/SW*Zu33i+3D0*Zn24r/SW*Zu33r+3D0*Mt
     &      /MW/SW*Zn44i*Zu63i/sb+3D0*Mt/MW/SW*Zn44r*Zu63r
     &      /sb)
      AAABR(847) = EE*Sqrt2*(Zn14i/CW*Zu33r-Zn14r/CW*Zu33i
     &      +3D0*Zn24i/SW*Zu33r-3D0*Zn24r/SW*Zu33i+3D0*Mt
     &      /MW/SW*Zn44i*Zu63r/sb-3D0*Mt/MW/SW*Zn44r*Zu63i
     &      /sb)
      AAABR(848) = EE*Sqrt2*(4D0*Zn14i/CW*Zu63i-4D0*Zn14r
     &      /CW*Zu63r-3D0*Mt/MW/SW*Zn44i*Zu33i/sb+3D0*Mt
     &      /MW/SW*Zn44r*Zu33r/sb)
      AAABR(849) = EE*Sqrt2*(4D0*Zn14i/CW*Zu63r+4D0*Zn14r
     &      /CW*Zu63i-3D0*Mt/MW/SW*Zn44i*Zu33r/sb-3D0*Mt
     &      /MW/SW*Zn44r*Zu33i/sb)
      AAABR(850) = EE*Sqrt2*(Zn13r/CW*Zu22+3D0*Zn23r/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn43r*Zu52/sb)
      AAABR(851) = EE*Sqrt2*(Zn13i/CW*Zu22+3D0*Zn23i/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn43i*Zu52/sb)
      AAABR(852) = EE*Sqrt2*(4D0*Zn13r/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn43r*Zu22/sb)
      AAABR(853) = EE*Sqrt2*(4D0*Zn13i/CW*Zu52-3D0*Mc/MW
     &      /SW*Zn43i*Zu22/sb)
      AAABR(854) = EE*Sqrt2*Zu14*(Zn13r/CW+3D0*Zn23r/SW)
      AAABR(855) = EE*Sqrt2*Zu14*(Zn13i/CW+3D0*Zn23i/SW)
      AAABR(856) = EE/CW*Sqrt2*Zn13r*Zu44
      AAABR(857) = EE/CW*Sqrt2*Zn13i*Zu44
      AAABR(858) = EE*Sqrt2*(Zn13r/CW*Zu25+3D0*Zn23r/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn43r*Zu55/sb)
      AAABR(859) = EE*Sqrt2*(Zn13i/CW*Zu25+3D0*Zn23i/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn43i*Zu55/sb)
      AAABR(860) = EE*Sqrt2*(4D0*Zn13r/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn43r*Zu25/sb)
      AAABR(861) = EE*Sqrt2*(4D0*Zn13i/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn43i*Zu25/sb)
      AAABR(862) = EE*Sqrt2*Zu14*(Zn12r/CW+3D0*Zn22r/SW)
      AAABR(863) = EE*Sqrt2*Zu14*(Zn12i/CW+3D0*Zn22i/SW)
      AAABR(864) = EE/CW*Sqrt2*Zn12r*Zu44
      AAABR(865) = EE/CW*Sqrt2*Zn12i*Zu44
      AAABR(866) = EE*Sqrt2*(Zn14r/CW*Zu25+3D0*Zn24r/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn44r*Zu55/sb)
      AAABR(867) = EE*Sqrt2*(Zn14i/CW*Zu25+3D0*Zn24i/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn44i*Zu55/sb)
      AAABR(868) = EE*Sqrt2*(4D0*Zn14r/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn44r*Zu25/sb)
      AAABR(869) = EE*Sqrt2*(4D0*Zn14i/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn44i*Zu25/sb)
      AAABR(870) = EE*Sqrt2*(Zn13i/CW*Zu36i+Zn13r/CW*Zu36r
     &      +3D0*Zn23i/SW*Zu36i+3D0*Zn23r/SW*Zu36r+3D0*Mt
     &      /MW/SW*Zn43i*Zu66i/sb+3D0*Mt/MW/SW*Zn43r*Zu66r
     &      /sb)
      AAABR(871) = EE*Sqrt2*(Zn13i/CW*Zu36r-Zn13r/CW*Zu36i
     &      +3D0*Zn23i/SW*Zu36r-3D0*Zn23r/SW*Zu36i+3D0*Mt
     &      /MW/SW*Zn43i*Zu66r/sb-3D0*Mt/MW/SW*Zn43r*Zu66i
     &      /sb)
      AAABR(872) = EE*Sqrt2*(4D0*Zn13i/CW*Zu66i-4D0*Zn13r
     &      /CW*Zu66r-3D0*Mt/MW/SW*Zn43i*Zu36i/sb+3D0*Mt
     &      /MW/SW*Zn43r*Zu36r/sb)
      AAABR(873) = EE*Sqrt2*(4D0*Zn13i/CW*Zu66r+4D0*Zn13r
     &      /CW*Zu66i-3D0*Mt/MW/SW*Zn43i*Zu36r/sb-3D0*Mt
     &      /MW/SW*Zn43r*Zu36i/sb)
      AAABR(874) = EE*Sqrt2*Zu14*(Zn14r/CW+3D0*Zn24r/SW)
      AAABR(875) = EE*Sqrt2*Zu14*(Zn14i/CW+3D0*Zn24i/SW)
      AAABR(876) = EE/CW*Sqrt2*Zn14r*Zu44
      AAABR(877) = EE/CW*Sqrt2*Zn14i*Zu44
      AAABR(878) = EE*Sqrt2*(Zn11i/CW*Zu36i+Zn11r/CW*Zu36r
     &      +3D0*Zn21i/SW*Zu36i+3D0*Zn21r/SW*Zu36r+3D0*Mt
     &      /MW/SW*Zn41i*Zu66i/sb+3D0*Mt/MW/SW*Zn41r*Zu66r
     &      /sb)
      AAABR(879) = EE*Sqrt2*(Zn11i/CW*Zu36r-Zn11r/CW*Zu36i
     &      +3D0*Zn21i/SW*Zu36r-3D0*Zn21r/SW*Zu36i+3D0*Mt
     &      /MW/SW*Zn41i*Zu66r/sb-3D0*Mt/MW/SW*Zn41r*Zu66i
     &      /sb)
      AAABR(880) = EE*Sqrt2*(4D0*Zn11i/CW*Zu66i-4D0*Zn11r
     &      /CW*Zu66r-3D0*Mt/MW/SW*Zn41i*Zu36i/sb+3D0*Mt
     &      /MW/SW*Zn41r*Zu36r/sb)
      AAABR(881) = EE*Sqrt2*(4D0*Zn11i/CW*Zu66r+4D0*Zn11r
     &      /CW*Zu66i-3D0*Mt/MW/SW*Zn41i*Zu36r/sb-3D0*Mt
     &      /MW/SW*Zn41r*Zu36i/sb)
      AAABR(882) = EE*Sqrt2*(Zn12i/CW*Zu36i+Zn12r/CW*Zu36r
     &      +3D0*Zn22i/SW*Zu36i+3D0*Zn22r/SW*Zu36r+3D0*Mt
     &      /MW/SW*Zn42i*Zu66i/sb+3D0*Mt/MW/SW*Zn42r*Zu66r
     &      /sb)
      AAABR(883) = EE*Sqrt2*(Zn12i/CW*Zu36r-Zn12r/CW*Zu36i
     &      +3D0*Zn22i/SW*Zu36r-3D0*Zn22r/SW*Zu36i+3D0*Mt
     &      /MW/SW*Zn42i*Zu66r/sb-3D0*Mt/MW/SW*Zn42r*Zu66i
     &      /sb)
      AAABR(884) = EE*Sqrt2*(4D0*Zn12i/CW*Zu66i-4D0*Zn12r
     &      /CW*Zu66r-3D0*Mt/MW/SW*Zn42i*Zu36i/sb+3D0*Mt
     &      /MW/SW*Zn42r*Zu36r/sb)
      AAABR(885) = EE*Sqrt2*(4D0*Zn12i/CW*Zu66r+4D0*Zn12r
     &      /CW*Zu66i-3D0*Mt/MW/SW*Zn42i*Zu36r/sb-3D0*Mt
     &      /MW/SW*Zn42r*Zu36i/sb)
      AAABR(886) = EE*Sqrt2*(Zn11r/CW*Zu25+3D0*Zn21r/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn41r*Zu55/sb)
      AAABR(887) = EE*Sqrt2*(Zn11i/CW*Zu25+3D0*Zn21i/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn41i*Zu55/sb)
      AAABR(888) = EE*Sqrt2*(4D0*Zn11r/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn41r*Zu25/sb)
      AAABR(889) = EE*Sqrt2*(4D0*Zn11i/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn41i*Zu25/sb)
      AAABR(890) = EE*Sqrt2*Zu14*(Zn11r/CW+3D0*Zn21r/SW)
      AAABR(891) = EE*Sqrt2*Zu14*(Zn11i/CW+3D0*Zn21i/SW)
      AAABR(892) = EE/CW*Sqrt2*Zn11r*Zu44
      AAABR(893) = EE/CW*Sqrt2*Zn11i*Zu44
      AAABR(894) = EE*Sqrt2*(Zn12r/CW*Zu25+3D0*Zn22r/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn42r*Zu55/sb)
      AAABR(895) = EE*Sqrt2*(Zn12i/CW*Zu25+3D0*Zn22i/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn42i*Zu55/sb)
      AAABR(896) = EE*Sqrt2*(4D0*Zn12r/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn42r*Zu25/sb)
      AAABR(897) = EE*Sqrt2*(4D0*Zn12i/CW*Zu55-3D0*Mc/MW
     &      /SW*Zn42i*Zu25/sb)
      AAABR(898) = EE*Sqrt2*(Zn14i/CW*Zu36i+Zn14r/CW*Zu36r
     &      +3D0*Zn24i/SW*Zu36i+3D0*Zn24r/SW*Zu36r+3D0*Mt
     &      /MW/SW*Zn44i*Zu66i/sb+3D0*Mt/MW/SW*Zn44r*Zu66r
     &      /sb)
      AAABR(899) = EE*Sqrt2*(Zn14i/CW*Zu36r-Zn14r/CW*Zu36i
     &      +3D0*Zn24i/SW*Zu36r-3D0*Zn24r/SW*Zu36i+3D0*Mt
     &      /MW/SW*Zn44i*Zu66r/sb-3D0*Mt/MW/SW*Zn44r*Zu66i
     &      /sb)
      AAABR(900) = EE*Sqrt2*(4D0*Zn14i/CW*Zu66i-4D0*Zn14r
     &      /CW*Zu66r-3D0*Mt/MW/SW*Zn44i*Zu36i/sb+3D0*Mt
     &      /MW/SW*Zn44r*Zu36r/sb)
      AAABR(901) = EE*Sqrt2*(4D0*Zn14i/CW*Zu66r+4D0*Zn14r
     &      /CW*Zu66i-3D0*Mt/MW/SW*Zn44i*Zu36r/sb-3D0*Mt
     &      /MW/SW*Zn44r*Zu36i/sb)
      AAABR(902) = EE/MW*Mt/SW
      AAABR(903) = EE/MW*Mc/SW
      AAABR(904) = EE/MW*Mt/SW*Zh23/sb
      AAABR(905) = EE/MW*Mt/SW*Zh33*cb/sb
      AAABR(906) = EE/MW*Mc/SW*Zh21/sb
      AAABR(907) = EE/MW*Mc/SW*Zh31*cb/sb
      AAABR(908) = EE/MW*Mc/SW*Zh22/sb
      AAABR(909) = EE/MW*Mc/SW*Zh32*cb/sb
      AAABR(910) = EE/MW*Mt/SW*Zh22/sb
      AAABR(911) = EE/MW*Mt/SW*Zh32*cb/sb
      AAABR(912) = EE/MW*Mc/SW*Zh23/sb
      AAABR(913) = EE/MW*Mc/SW*Zh33*cb/sb
      AAABR(914) = EE/MW*Mt/SW*Zh21/sb
      AAABR(915) = EE/MW*Mt/SW*Zh31*cb/sb
      AAABR(916) = GG*Sqrt2*Zu11
      AAABR(917) = GG*Sqrt2*Zu41
      AAABR(918) = GG*Sqrt2*Zu33r
      AAABR(919) = GG*Sqrt2*Zu33i
      AAABR(920) = GG*Sqrt2*Zu63r
      AAABR(921) = GG*Sqrt2*Zu63i
      AAABR(922) = GG*Sqrt2*Zu22
      AAABR(923) = GG*Sqrt2*Zu52
      AAABR(924) = GG*Sqrt2*Zu36r
      AAABR(925) = GG*Sqrt2*Zu36i
      AAABR(926) = GG*Sqrt2*Zu66r
      AAABR(927) = GG*Sqrt2*Zu66i
      AAABR(928) = GG*Sqrt2*Zu25
      AAABR(929) = GG*Sqrt2*Zu55
      AAABR(930) = GG*Sqrt2*Zu14
      AAABR(931) = GG*Sqrt2*Zu44
      AAABR(932) = EE*(Zn11i/CW*Zn31r*cb+Zn11i/CW*Zn41r*sb
     &      +Zn11r/CW*Zn31i*cb+Zn11r/CW*Zn41i*sb-Zn21i/SW
     &      *Zn31r*cb-Zn21i/SW*Zn41r*sb-Zn21r/SW*Zn31i*cb
     &      -Zn21r/SW*Zn41i*sb)
      AAABR(933) = EE*(Zn11i/CW*Zn31i*cb+Zn11i/CW*Zn41i*sb
     &      -Zn11r/CW*Zn31r*cb-Zn11r/CW*Zn41r*sb-Zn21i/SW
     &      *Zn31i*cb-Zn21i/SW*Zn41i*sb+Zn21r/SW*Zn31r*cb
     &      +Zn21r/SW*Zn41r*sb)
      AAABR(934) = EE*(Zn11i/CW*Zn34r*cb+Zn11i/CW*Zn44r*sb
     &      +Zn11r/CW*Zn34i*cb+Zn11r/CW*Zn44i*sb+Zn14i/CW
     &      *Zn31r*cb+Zn14i/CW*Zn41r*sb+Zn14r/CW*Zn31i*cb
     &      +Zn14r/CW*Zn41i*sb-Zn21i/SW*Zn34r*cb-Zn21i/SW
     &      *Zn44r*sb-Zn21r/SW*Zn34i*cb-Zn21r/SW*Zn44i*sb
     &      -Zn24i/SW*Zn31r*cb-Zn24i/SW*Zn41r*sb-Zn24r/SW
     &      *Zn31i*cb-Zn24r/SW*Zn41i*sb)
      AAABR(935) = EE*(Zn11i/CW*Zn34i*cb+Zn11i/CW*Zn44i*sb
     &      -Zn11r/CW*Zn34r*cb-Zn11r/CW*Zn44r*sb+Zn14i/CW
     &      *Zn31i*cb+Zn14i/CW*Zn41i*sb-Zn14r/CW*Zn31r*cb
     &      -Zn14r/CW*Zn41r*sb-Zn21i/SW*Zn34i*cb-Zn21i/SW
     &      *Zn44i*sb+Zn21r/SW*Zn34r*cb+Zn21r/SW*Zn44r*sb
     &      -Zn24i/SW*Zn31i*cb-Zn24i/SW*Zn41i*sb+Zn24r/SW
     &      *Zn31r*cb+Zn24r/SW*Zn41r*sb)
      AAABR(936) = EE*(Zn11i/CW*Zn32r*cb+Zn11i/CW*Zn42r*sb
     &      +Zn11r/CW*Zn32i*cb+Zn11r/CW*Zn42i*sb+Zn12i/CW
     &      *Zn31r*cb+Zn12i/CW*Zn41r*sb+Zn12r/CW*Zn31i*cb
     &      +Zn12r/CW*Zn41i*sb-Zn21i/SW*Zn32r*cb-Zn21i/SW
     &      *Zn42r*sb-Zn21r/SW*Zn32i*cb-Zn21r/SW*Zn42i*sb
     &      -Zn22i/SW*Zn31r*cb-Zn22i/SW*Zn41r*sb-Zn22r/SW
     &      *Zn31i*cb-Zn22r/SW*Zn41i*sb)
      AAABR(937) = EE*(Zn11i/CW*Zn32i*cb+Zn11i/CW*Zn42i*sb
     &      -Zn11r/CW*Zn32r*cb-Zn11r/CW*Zn42r*sb+Zn12i/CW
     &      *Zn31i*cb+Zn12i/CW*Zn41i*sb-Zn12r/CW*Zn31r*cb
     &      -Zn12r/CW*Zn41r*sb-Zn21i/SW*Zn32i*cb-Zn21i/SW
     &      *Zn42i*sb+Zn21r/SW*Zn32r*cb+Zn21r/SW*Zn42r*sb
     &      -Zn22i/SW*Zn31i*cb-Zn22i/SW*Zn41i*sb+Zn22r/SW
     &      *Zn31r*cb+Zn22r/SW*Zn41r*sb)
      AAABR(938) = EE*(Zn11i/CW*Zn33r*cb+Zn11i/CW*Zn43r*sb
     &      +Zn11r/CW*Zn33i*cb+Zn11r/CW*Zn43i*sb+Zn13i/CW
     &      *Zn31r*cb+Zn13i/CW*Zn41r*sb+Zn13r/CW*Zn31i*cb
     &      +Zn13r/CW*Zn41i*sb-Zn21i/SW*Zn33r*cb-Zn21i/SW
     &      *Zn43r*sb-Zn21r/SW*Zn33i*cb-Zn21r/SW*Zn43i*sb
     &      -Zn23i/SW*Zn31r*cb-Zn23i/SW*Zn41r*sb-Zn23r/SW
     &      *Zn31i*cb-Zn23r/SW*Zn41i*sb)
      AAABR(939) = EE*(Zn11i/CW*Zn33i*cb+Zn11i/CW*Zn43i*sb
     &      -Zn11r/CW*Zn33r*cb-Zn11r/CW*Zn43r*sb+Zn13i/CW
     &      *Zn31i*cb+Zn13i/CW*Zn41i*sb-Zn13r/CW*Zn31r*cb
     &      -Zn13r/CW*Zn41r*sb-Zn21i/SW*Zn33i*cb-Zn21i/SW
     &      *Zn43i*sb+Zn21r/SW*Zn33r*cb+Zn21r/SW*Zn43r*sb
     &      -Zn23i/SW*Zn31i*cb-Zn23i/SW*Zn41i*sb+Zn23r/SW
     &      *Zn31r*cb+Zn23r/SW*Zn41r*sb)
      AAABR(940) = EE*(Zh13/CW*Zn11i*Zn33i-Zh13/CW*Zn11r
     &      *Zn33r+Zh13/CW*Zn13i*Zn31i-Zh13/CW*Zn13r*Zn31r
     &      -Zh23/CW*Zn11i*Zn43i+Zh23/CW*Zn11r*Zn43r-Zh23
     &      /CW*Zn13i*Zn41i+Zh23/CW*Zn13r*Zn41r-Zh13/SW*Zn21i
     &      *Zn33i+Zh13/SW*Zn21r*Zn33r-Zh13/SW*Zn23i*Zn31i
     &      +Zh13/SW*Zn23r*Zn31r+Zh23/SW*Zn21i*Zn43i-Zh23
     &      /SW*Zn21r*Zn43r+Zh23/SW*Zn23i*Zn41i-Zh23/SW*Zn23r
     &      *Zn41r-Zh33/CW*Zn11i*Zn33r*sb+Zh33/CW*Zn11i*Zn43r
     &      *cb-Zh33/CW*Zn11r*Zn33i*sb+Zh33/CW*Zn11r*Zn43i
     &      *cb-Zh33/CW*Zn13i*Zn31r*sb+Zh33/CW*Zn13i*Zn41r
     &      *cb-Zh33/CW*Zn13r*Zn31i*sb+Zh33/CW*Zn13r*Zn41i
     &      *cb+Zh33/SW*Zn21i*Zn33r*sb-Zh33/SW*Zn21i*Zn43r
     &      *cb+Zh33/SW*Zn21r*Zn33i*sb-Zh33/SW*Zn21r*Zn43i
     &      *cb+Zh33/SW*Zn23i*Zn31r*sb-Zh33/SW*Zn23i*Zn41r
     &      *cb+Zh33/SW*Zn23r*Zn31i*sb-Zh33/SW*Zn23r*Zn41i
     &      *cb)
      AAABR(941) = EE*(Zh13/CW*Zn11i*Zn33r+Zh13/CW*Zn11r
     &      *Zn33i+Zh13/CW*Zn13i*Zn31r+Zh13/CW*Zn13r*Zn31i
     &      -Zh23/CW*Zn11i*Zn43r-Zh23/CW*Zn11r*Zn43i-Zh23
     &      /CW*Zn13i*Zn41r-Zh23/CW*Zn13r*Zn41i-Zh13/SW*Zn21i
     &      *Zn33r-Zh13/SW*Zn21r*Zn33i-Zh13/SW*Zn23i*Zn31r
     &      -Zh13/SW*Zn23r*Zn31i+Zh23/SW*Zn21i*Zn43r+Zh23
     &      /SW*Zn21r*Zn43i+Zh23/SW*Zn23i*Zn41r+Zh23/SW*Zn23r
     &      *Zn41i+Zh33/CW*Zn11i*Zn33i*sb-Zh33/CW*Zn11i*Zn43i
     &      *cb-Zh33/CW*Zn11r*Zn33r*sb+Zh33/CW*Zn11r*Zn43r
     &      *cb+Zh33/CW*Zn13i*Zn31i*sb-Zh33/CW*Zn13i*Zn41i
     &      *cb-Zh33/CW*Zn13r*Zn31r*sb+Zh33/CW*Zn13r*Zn41r
     &      *cb-Zh33/SW*Zn21i*Zn33i*sb+Zh33/SW*Zn21i*Zn43i
     &      *cb+Zh33/SW*Zn21r*Zn33r*sb-Zh33/SW*Zn21r*Zn43r
     &      *cb-Zh33/SW*Zn23i*Zn31i*sb+Zh33/SW*Zn23i*Zn41i
     &      *cb+Zh33/SW*Zn23r*Zn31r*sb-Zh33/SW*Zn23r*Zn41r
     &      *cb)
      AAABR(942) = EE*(Zh11/CW*Zn11i*Zn33i-Zh11/CW*Zn11r
     &      *Zn33r+Zh11/CW*Zn13i*Zn31i-Zh11/CW*Zn13r*Zn31r
     &      -Zh21/CW*Zn11i*Zn43i+Zh21/CW*Zn11r*Zn43r-Zh21
     &      /CW*Zn13i*Zn41i+Zh21/CW*Zn13r*Zn41r-Zh11/SW*Zn21i
     &      *Zn33i+Zh11/SW*Zn21r*Zn33r-Zh11/SW*Zn23i*Zn31i
     &      +Zh11/SW*Zn23r*Zn31r+Zh21/SW*Zn21i*Zn43i-Zh21
     &      /SW*Zn21r*Zn43r+Zh21/SW*Zn23i*Zn41i-Zh21/SW*Zn23r
     &      *Zn41r-Zh31/CW*Zn11i*Zn33r*sb+Zh31/CW*Zn11i*Zn43r
     &      *cb-Zh31/CW*Zn11r*Zn33i*sb+Zh31/CW*Zn11r*Zn43i
     &      *cb-Zh31/CW*Zn13i*Zn31r*sb+Zh31/CW*Zn13i*Zn41r
     &      *cb-Zh31/CW*Zn13r*Zn31i*sb+Zh31/CW*Zn13r*Zn41i
     &      *cb+Zh31/SW*Zn21i*Zn33r*sb-Zh31/SW*Zn21i*Zn43r
     &      *cb+Zh31/SW*Zn21r*Zn33i*sb-Zh31/SW*Zn21r*Zn43i
     &      *cb+Zh31/SW*Zn23i*Zn31r*sb-Zh31/SW*Zn23i*Zn41r
     &      *cb+Zh31/SW*Zn23r*Zn31i*sb-Zh31/SW*Zn23r*Zn41i
     &      *cb)
      AAABR(943) = EE*(Zh11/CW*Zn11i*Zn33r+Zh11/CW*Zn11r
     &      *Zn33i+Zh11/CW*Zn13i*Zn31r+Zh11/CW*Zn13r*Zn31i
     &      -Zh21/CW*Zn11i*Zn43r-Zh21/CW*Zn11r*Zn43i-Zh21
     &      /CW*Zn13i*Zn41r-Zh21/CW*Zn13r*Zn41i-Zh11/SW*Zn21i
     &      *Zn33r-Zh11/SW*Zn21r*Zn33i-Zh11/SW*Zn23i*Zn31r
     &      -Zh11/SW*Zn23r*Zn31i+Zh21/SW*Zn21i*Zn43r+Zh21
     &      /SW*Zn21r*Zn43i+Zh21/SW*Zn23i*Zn41r+Zh21/SW*Zn23r
     &      *Zn41i+Zh31/CW*Zn11i*Zn33i*sb-Zh31/CW*Zn11i*Zn43i
     &      *cb-Zh31/CW*Zn11r*Zn33r*sb+Zh31/CW*Zn11r*Zn43r
     &      *cb+Zh31/CW*Zn13i*Zn31i*sb-Zh31/CW*Zn13i*Zn41i
     &      *cb-Zh31/CW*Zn13r*Zn31r*sb+Zh31/CW*Zn13r*Zn41r
     &      *cb-Zh31/SW*Zn21i*Zn33i*sb+Zh31/SW*Zn21i*Zn43i
     &      *cb+Zh31/SW*Zn21r*Zn33r*sb-Zh31/SW*Zn21r*Zn43r
     &      *cb-Zh31/SW*Zn23i*Zn31i*sb+Zh31/SW*Zn23i*Zn41i
     &      *cb+Zh31/SW*Zn23r*Zn31r*sb-Zh31/SW*Zn23r*Zn41r
     &      *cb)
      AAABR(944) = EE*(Zh12/CW*Zn11i*Zn34i-Zh12/CW*Zn11r
     &      *Zn34r+Zh12/CW*Zn14i*Zn31i-Zh12/CW*Zn14r*Zn31r
     &      -Zh22/CW*Zn11i*Zn44i+Zh22/CW*Zn11r*Zn44r-Zh22
     &      /CW*Zn14i*Zn41i+Zh22/CW*Zn14r*Zn41r-Zh12/SW*Zn21i
     &      *Zn34i+Zh12/SW*Zn21r*Zn34r-Zh12/SW*Zn24i*Zn31i
     &      +Zh12/SW*Zn24r*Zn31r+Zh22/SW*Zn21i*Zn44i-Zh22
     &      /SW*Zn21r*Zn44r+Zh22/SW*Zn24i*Zn41i-Zh22/SW*Zn24r
     &      *Zn41r-Zh32/CW*Zn11i*Zn34r*sb+Zh32/CW*Zn11i*Zn44r
     &      *cb-Zh32/CW*Zn11r*Zn34i*sb+Zh32/CW*Zn11r*Zn44i
     &      *cb-Zh32/CW*Zn14i*Zn31r*sb+Zh32/CW*Zn14i*Zn41r
     &      *cb-Zh32/CW*Zn14r*Zn31i*sb+Zh32/CW*Zn14r*Zn41i
     &      *cb+Zh32/SW*Zn21i*Zn34r*sb-Zh32/SW*Zn21i*Zn44r
     &      *cb+Zh32/SW*Zn21r*Zn34i*sb-Zh32/SW*Zn21r*Zn44i
     &      *cb+Zh32/SW*Zn24i*Zn31r*sb-Zh32/SW*Zn24i*Zn41r
     &      *cb+Zh32/SW*Zn24r*Zn31i*sb-Zh32/SW*Zn24r*Zn41i
     &      *cb)
      AAABR(945) = EE*(Zh12/CW*Zn11i*Zn34r+Zh12/CW*Zn11r
     &      *Zn34i+Zh12/CW*Zn14i*Zn31r+Zh12/CW*Zn14r*Zn31i
     &      -Zh22/CW*Zn11i*Zn44r-Zh22/CW*Zn11r*Zn44i-Zh22
     &      /CW*Zn14i*Zn41r-Zh22/CW*Zn14r*Zn41i-Zh12/SW*Zn21i
     &      *Zn34r-Zh12/SW*Zn21r*Zn34i-Zh12/SW*Zn24i*Zn31r
     &      -Zh12/SW*Zn24r*Zn31i+Zh22/SW*Zn21i*Zn44r+Zh22
     &      /SW*Zn21r*Zn44i+Zh22/SW*Zn24i*Zn41r+Zh22/SW*Zn24r
     &      *Zn41i+Zh32/CW*Zn11i*Zn34i*sb-Zh32/CW*Zn11i*Zn44i
     &      *cb-Zh32/CW*Zn11r*Zn34r*sb+Zh32/CW*Zn11r*Zn44r
     &      *cb+Zh32/CW*Zn14i*Zn31i*sb-Zh32/CW*Zn14i*Zn41i
     &      *cb-Zh32/CW*Zn14r*Zn31r*sb+Zh32/CW*Zn14r*Zn41r
     &      *cb-Zh32/SW*Zn21i*Zn34i*sb+Zh32/SW*Zn21i*Zn44i
     &      *cb+Zh32/SW*Zn21r*Zn34r*sb-Zh32/SW*Zn21r*Zn44r
     &      *cb-Zh32/SW*Zn24i*Zn31i*sb+Zh32/SW*Zn24i*Zn41i
     &      *cb+Zh32/SW*Zn24r*Zn31r*sb-Zh32/SW*Zn24r*Zn41r
     &      *cb)
      AAABR(946) = EE*(Zh12/CW*Zn11i*Zn32i-Zh12/CW*Zn11r
     &      *Zn32r+Zh12/CW*Zn12i*Zn31i-Zh12/CW*Zn12r*Zn31r
     &      -Zh22/CW*Zn11i*Zn42i+Zh22/CW*Zn11r*Zn42r-Zh22
     &      /CW*Zn12i*Zn41i+Zh22/CW*Zn12r*Zn41r-Zh12/SW*Zn21i
     &      *Zn32i+Zh12/SW*Zn21r*Zn32r-Zh12/SW*Zn22i*Zn31i
     &      +Zh12/SW*Zn22r*Zn31r+Zh22/SW*Zn21i*Zn42i-Zh22
     &      /SW*Zn21r*Zn42r+Zh22/SW*Zn22i*Zn41i-Zh22/SW*Zn22r
     &      *Zn41r-Zh32/CW*Zn11i*Zn32r*sb+Zh32/CW*Zn11i*Zn42r
     &      *cb-Zh32/CW*Zn11r*Zn32i*sb+Zh32/CW*Zn11r*Zn42i
     &      *cb-Zh32/CW*Zn12i*Zn31r*sb+Zh32/CW*Zn12i*Zn41r
     &      *cb-Zh32/CW*Zn12r*Zn31i*sb+Zh32/CW*Zn12r*Zn41i
     &      *cb+Zh32/SW*Zn21i*Zn32r*sb-Zh32/SW*Zn21i*Zn42r
     &      *cb+Zh32/SW*Zn21r*Zn32i*sb-Zh32/SW*Zn21r*Zn42i
     &      *cb+Zh32/SW*Zn22i*Zn31r*sb-Zh32/SW*Zn22i*Zn41r
     &      *cb+Zh32/SW*Zn22r*Zn31i*sb-Zh32/SW*Zn22r*Zn41i
     &      *cb)
      AAABR(947) = EE*(Zh12/CW*Zn11i*Zn32r+Zh12/CW*Zn11r
     &      *Zn32i+Zh12/CW*Zn12i*Zn31r+Zh12/CW*Zn12r*Zn31i
     &      -Zh22/CW*Zn11i*Zn42r-Zh22/CW*Zn11r*Zn42i-Zh22
     &      /CW*Zn12i*Zn41r-Zh22/CW*Zn12r*Zn41i-Zh12/SW*Zn21i
     &      *Zn32r-Zh12/SW*Zn21r*Zn32i-Zh12/SW*Zn22i*Zn31r
     &      -Zh12/SW*Zn22r*Zn31i+Zh22/SW*Zn21i*Zn42r+Zh22
     &      /SW*Zn21r*Zn42i+Zh22/SW*Zn22i*Zn41r+Zh22/SW*Zn22r
     &      *Zn41i+Zh32/CW*Zn11i*Zn32i*sb-Zh32/CW*Zn11i*Zn42i
     &      *cb-Zh32/CW*Zn11r*Zn32r*sb+Zh32/CW*Zn11r*Zn42r
     &      *cb+Zh32/CW*Zn12i*Zn31i*sb-Zh32/CW*Zn12i*Zn41i
     &      *cb-Zh32/CW*Zn12r*Zn31r*sb+Zh32/CW*Zn12r*Zn41r
     &      *cb-Zh32/SW*Zn21i*Zn32i*sb+Zh32/SW*Zn21i*Zn42i
     &      *cb+Zh32/SW*Zn21r*Zn32r*sb-Zh32/SW*Zn21r*Zn42r
     &      *cb-Zh32/SW*Zn22i*Zn31i*sb+Zh32/SW*Zn22i*Zn41i
     &      *cb+Zh32/SW*Zn22r*Zn31r*sb-Zh32/SW*Zn22r*Zn41r
     &      *cb)
      AAABR(948) = EE*(Zh13/CW*Zn11i*Zn34i-Zh13/CW*Zn11r
     &      *Zn34r+Zh13/CW*Zn14i*Zn31i-Zh13/CW*Zn14r*Zn31r
     &      -Zh23/CW*Zn11i*Zn44i+Zh23/CW*Zn11r*Zn44r-Zh23
     &      /CW*Zn14i*Zn41i+Zh23/CW*Zn14r*Zn41r-Zh13/SW*Zn21i
     &      *Zn34i+Zh13/SW*Zn21r*Zn34r-Zh13/SW*Zn24i*Zn31i
     &      +Zh13/SW*Zn24r*Zn31r+Zh23/SW*Zn21i*Zn44i-Zh23
     &      /SW*Zn21r*Zn44r+Zh23/SW*Zn24i*Zn41i-Zh23/SW*Zn24r
     &      *Zn41r-Zh33/CW*Zn11i*Zn34r*sb+Zh33/CW*Zn11i*Zn44r
     &      *cb-Zh33/CW*Zn11r*Zn34i*sb+Zh33/CW*Zn11r*Zn44i
     &      *cb-Zh33/CW*Zn14i*Zn31r*sb+Zh33/CW*Zn14i*Zn41r
     &      *cb-Zh33/CW*Zn14r*Zn31i*sb+Zh33/CW*Zn14r*Zn41i
     &      *cb+Zh33/SW*Zn21i*Zn34r*sb-Zh33/SW*Zn21i*Zn44r
     &      *cb+Zh33/SW*Zn21r*Zn34i*sb-Zh33/SW*Zn21r*Zn44i
     &      *cb+Zh33/SW*Zn24i*Zn31r*sb-Zh33/SW*Zn24i*Zn41r
     &      *cb+Zh33/SW*Zn24r*Zn31i*sb-Zh33/SW*Zn24r*Zn41i
     &      *cb)
      AAABR(949) = EE*(Zh13/CW*Zn11i*Zn34r+Zh13/CW*Zn11r
     &      *Zn34i+Zh13/CW*Zn14i*Zn31r+Zh13/CW*Zn14r*Zn31i
     &      -Zh23/CW*Zn11i*Zn44r-Zh23/CW*Zn11r*Zn44i-Zh23
     &      /CW*Zn14i*Zn41r-Zh23/CW*Zn14r*Zn41i-Zh13/SW*Zn21i
     &      *Zn34r-Zh13/SW*Zn21r*Zn34i-Zh13/SW*Zn24i*Zn31r
     &      -Zh13/SW*Zn24r*Zn31i+Zh23/SW*Zn21i*Zn44r+Zh23
     &      /SW*Zn21r*Zn44i+Zh23/SW*Zn24i*Zn41r+Zh23/SW*Zn24r
     &      *Zn41i+Zh33/CW*Zn11i*Zn34i*sb-Zh33/CW*Zn11i*Zn44i
     &      *cb-Zh33/CW*Zn11r*Zn34r*sb+Zh33/CW*Zn11r*Zn44r
     &      *cb+Zh33/CW*Zn14i*Zn31i*sb-Zh33/CW*Zn14i*Zn41i
     &      *cb-Zh33/CW*Zn14r*Zn31r*sb+Zh33/CW*Zn14r*Zn41r
     &      *cb-Zh33/SW*Zn21i*Zn34i*sb+Zh33/SW*Zn21i*Zn44i
     &      *cb+Zh33/SW*Zn21r*Zn34r*sb-Zh33/SW*Zn21r*Zn44r
     &      *cb-Zh33/SW*Zn24i*Zn31i*sb+Zh33/SW*Zn24i*Zn41i
     &      *cb+Zh33/SW*Zn24r*Zn31r*sb-Zh33/SW*Zn24r*Zn41r
     &      *cb)
      AAABR(950) = EE*(Zh11/CW*Zn11i*Zn32i-Zh11/CW*Zn11r
     &      *Zn32r+Zh11/CW*Zn12i*Zn31i-Zh11/CW*Zn12r*Zn31r
     &      -Zh21/CW*Zn11i*Zn42i+Zh21/CW*Zn11r*Zn42r-Zh21
     &      /CW*Zn12i*Zn41i+Zh21/CW*Zn12r*Zn41r-Zh11/SW*Zn21i
     &      *Zn32i+Zh11/SW*Zn21r*Zn32r-Zh11/SW*Zn22i*Zn31i
     &      +Zh11/SW*Zn22r*Zn31r+Zh21/SW*Zn21i*Zn42i-Zh21
     &      /SW*Zn21r*Zn42r+Zh21/SW*Zn22i*Zn41i-Zh21/SW*Zn22r
     &      *Zn41r-Zh31/CW*Zn11i*Zn32r*sb+Zh31/CW*Zn11i*Zn42r
     &      *cb-Zh31/CW*Zn11r*Zn32i*sb+Zh31/CW*Zn11r*Zn42i
     &      *cb-Zh31/CW*Zn12i*Zn31r*sb+Zh31/CW*Zn12i*Zn41r
     &      *cb-Zh31/CW*Zn12r*Zn31i*sb+Zh31/CW*Zn12r*Zn41i
     &      *cb+Zh31/SW*Zn21i*Zn32r*sb-Zh31/SW*Zn21i*Zn42r
     &      *cb+Zh31/SW*Zn21r*Zn32i*sb-Zh31/SW*Zn21r*Zn42i
     &      *cb+Zh31/SW*Zn22i*Zn31r*sb-Zh31/SW*Zn22i*Zn41r
     &      *cb+Zh31/SW*Zn22r*Zn31i*sb-Zh31/SW*Zn22r*Zn41i
     &      *cb)
      AAABR(951) = EE*(Zh11/CW*Zn11i*Zn32r+Zh11/CW*Zn11r
     &      *Zn32i+Zh11/CW*Zn12i*Zn31r+Zh11/CW*Zn12r*Zn31i
     &      -Zh21/CW*Zn11i*Zn42r-Zh21/CW*Zn11r*Zn42i-Zh21
     &      /CW*Zn12i*Zn41r-Zh21/CW*Zn12r*Zn41i-Zh11/SW*Zn21i
     &      *Zn32r-Zh11/SW*Zn21r*Zn32i-Zh11/SW*Zn22i*Zn31r
     &      -Zh11/SW*Zn22r*Zn31i+Zh21/SW*Zn21i*Zn42r+Zh21
     &      /SW*Zn21r*Zn42i+Zh21/SW*Zn22i*Zn41r+Zh21/SW*Zn22r
     &      *Zn41i+Zh31/CW*Zn11i*Zn32i*sb-Zh31/CW*Zn11i*Zn42i
     &      *cb-Zh31/CW*Zn11r*Zn32r*sb+Zh31/CW*Zn11r*Zn42r
     &      *cb+Zh31/CW*Zn12i*Zn31i*sb-Zh31/CW*Zn12i*Zn41i
     &      *cb-Zh31/CW*Zn12r*Zn31r*sb+Zh31/CW*Zn12r*Zn41r
     &      *cb-Zh31/SW*Zn21i*Zn32i*sb+Zh31/SW*Zn21i*Zn42i
     &      *cb+Zh31/SW*Zn21r*Zn32r*sb-Zh31/SW*Zn21r*Zn42r
     &      *cb-Zh31/SW*Zn22i*Zn31i*sb+Zh31/SW*Zn22i*Zn41i
     &      *cb+Zh31/SW*Zn22r*Zn31r*sb-Zh31/SW*Zn22r*Zn41r
     &      *cb)
      AAABR(952) = EE*(Zh13/CW*Zn11i*Zn31i-Zh13/CW*Zn11r
     &      *Zn31r-Zh23/CW*Zn11i*Zn41i+Zh23/CW*Zn11r*Zn41r
     &      -Zh13/SW*Zn21i*Zn31i+Zh13/SW*Zn21r*Zn31r+Zh23
     &      /SW*Zn21i*Zn41i-Zh23/SW*Zn21r*Zn41r-Zh33/CW*Zn11i
     &      *Zn31r*sb+Zh33/CW*Zn11i*Zn41r*cb-Zh33/CW*Zn11r
     &      *Zn31i*sb+Zh33/CW*Zn11r*Zn41i*cb+Zh33/SW*Zn21i
     &      *Zn31r*sb-Zh33/SW*Zn21i*Zn41r*cb+Zh33/SW*Zn21r
     &      *Zn31i*sb-Zh33/SW*Zn21r*Zn41i*cb)
      AAABR(953) = EE*(Zh13/CW*Zn11i*Zn31r+Zh13/CW*Zn11r
     &      *Zn31i-Zh23/CW*Zn11i*Zn41r-Zh23/CW*Zn11r*Zn41i
     &      -Zh13/SW*Zn21i*Zn31r-Zh13/SW*Zn21r*Zn31i+Zh23
     &      /SW*Zn21i*Zn41r+Zh23/SW*Zn21r*Zn41i+Zh33/CW*Zn11i
     &      *Zn31i*sb-Zh33/CW*Zn11i*Zn41i*cb-Zh33/CW*Zn11r
     &      *Zn31r*sb+Zh33/CW*Zn11r*Zn41r*cb-Zh33/SW*Zn21i
     &      *Zn31i*sb+Zh33/SW*Zn21i*Zn41i*cb+Zh33/SW*Zn21r
     &      *Zn31r*sb-Zh33/SW*Zn21r*Zn41r*cb)
      AAABR(954) = EE*(Zh13/CW*Zn11i*Zn32i-Zh13/CW*Zn11r
     &      *Zn32r+Zh13/CW*Zn12i*Zn31i-Zh13/CW*Zn12r*Zn31r
     &      -Zh23/CW*Zn11i*Zn42i+Zh23/CW*Zn11r*Zn42r-Zh23
     &      /CW*Zn12i*Zn41i+Zh23/CW*Zn12r*Zn41r-Zh13/SW*Zn21i
     &      *Zn32i+Zh13/SW*Zn21r*Zn32r-Zh13/SW*Zn22i*Zn31i
     &      +Zh13/SW*Zn22r*Zn31r+Zh23/SW*Zn21i*Zn42i-Zh23
     &      /SW*Zn21r*Zn42r+Zh23/SW*Zn22i*Zn41i-Zh23/SW*Zn22r
     &      *Zn41r-Zh33/CW*Zn11i*Zn32r*sb+Zh33/CW*Zn11i*Zn42r
     &      *cb-Zh33/CW*Zn11r*Zn32i*sb+Zh33/CW*Zn11r*Zn42i
     &      *cb-Zh33/CW*Zn12i*Zn31r*sb+Zh33/CW*Zn12i*Zn41r
     &      *cb-Zh33/CW*Zn12r*Zn31i*sb+Zh33/CW*Zn12r*Zn41i
     &      *cb+Zh33/SW*Zn21i*Zn32r*sb-Zh33/SW*Zn21i*Zn42r
     &      *cb+Zh33/SW*Zn21r*Zn32i*sb-Zh33/SW*Zn21r*Zn42i
     &      *cb+Zh33/SW*Zn22i*Zn31r*sb-Zh33/SW*Zn22i*Zn41r
     &      *cb+Zh33/SW*Zn22r*Zn31i*sb-Zh33/SW*Zn22r*Zn41i
     &      *cb)
      AAABR(955) = EE*(Zh13/CW*Zn11i*Zn32r+Zh13/CW*Zn11r
     &      *Zn32i+Zh13/CW*Zn12i*Zn31r+Zh13/CW*Zn12r*Zn31i
     &      -Zh23/CW*Zn11i*Zn42r-Zh23/CW*Zn11r*Zn42i-Zh23
     &      /CW*Zn12i*Zn41r-Zh23/CW*Zn12r*Zn41i-Zh13/SW*Zn21i
     &      *Zn32r-Zh13/SW*Zn21r*Zn32i-Zh13/SW*Zn22i*Zn31r
     &      -Zh13/SW*Zn22r*Zn31i+Zh23/SW*Zn21i*Zn42r+Zh23
     &      /SW*Zn21r*Zn42i+Zh23/SW*Zn22i*Zn41r+Zh23/SW*Zn22r
     &      *Zn41i+Zh33/CW*Zn11i*Zn32i*sb-Zh33/CW*Zn11i*Zn42i
     &      *cb-Zh33/CW*Zn11r*Zn32r*sb+Zh33/CW*Zn11r*Zn42r
     &      *cb+Zh33/CW*Zn12i*Zn31i*sb-Zh33/CW*Zn12i*Zn41i
     &      *cb-Zh33/CW*Zn12r*Zn31r*sb+Zh33/CW*Zn12r*Zn41r
     &      *cb-Zh33/SW*Zn21i*Zn32i*sb+Zh33/SW*Zn21i*Zn42i
     &      *cb+Zh33/SW*Zn21r*Zn32r*sb-Zh33/SW*Zn21r*Zn42r
     &      *cb-Zh33/SW*Zn22i*Zn31i*sb+Zh33/SW*Zn22i*Zn41i
     &      *cb+Zh33/SW*Zn22r*Zn31r*sb-Zh33/SW*Zn22r*Zn41r
     &      *cb)
      AAABR(956) = EE*(Zh12/CW*Zn11i*Zn33i-Zh12/CW*Zn11r
     &      *Zn33r+Zh12/CW*Zn13i*Zn31i-Zh12/CW*Zn13r*Zn31r
     &      -Zh22/CW*Zn11i*Zn43i+Zh22/CW*Zn11r*Zn43r-Zh22
     &      /CW*Zn13i*Zn41i+Zh22/CW*Zn13r*Zn41r-Zh12/SW*Zn21i
     &      *Zn33i+Zh12/SW*Zn21r*Zn33r-Zh12/SW*Zn23i*Zn31i
     &      +Zh12/SW*Zn23r*Zn31r+Zh22/SW*Zn21i*Zn43i-Zh22
     &      /SW*Zn21r*Zn43r+Zh22/SW*Zn23i*Zn41i-Zh22/SW*Zn23r
     &      *Zn41r-Zh32/CW*Zn11i*Zn33r*sb+Zh32/CW*Zn11i*Zn43r
     &      *cb-Zh32/CW*Zn11r*Zn33i*sb+Zh32/CW*Zn11r*Zn43i
     &      *cb-Zh32/CW*Zn13i*Zn31r*sb+Zh32/CW*Zn13i*Zn41r
     &      *cb-Zh32/CW*Zn13r*Zn31i*sb+Zh32/CW*Zn13r*Zn41i
     &      *cb+Zh32/SW*Zn21i*Zn33r*sb-Zh32/SW*Zn21i*Zn43r
     &      *cb+Zh32/SW*Zn21r*Zn33i*sb-Zh32/SW*Zn21r*Zn43i
     &      *cb+Zh32/SW*Zn23i*Zn31r*sb-Zh32/SW*Zn23i*Zn41r
     &      *cb+Zh32/SW*Zn23r*Zn31i*sb-Zh32/SW*Zn23r*Zn41i
     &      *cb)
      AAABR(957) = EE*(Zh12/CW*Zn11i*Zn33r+Zh12/CW*Zn11r
     &      *Zn33i+Zh12/CW*Zn13i*Zn31r+Zh12/CW*Zn13r*Zn31i
     &      -Zh22/CW*Zn11i*Zn43r-Zh22/CW*Zn11r*Zn43i-Zh22
     &      /CW*Zn13i*Zn41r-Zh22/CW*Zn13r*Zn41i-Zh12/SW*Zn21i
     &      *Zn33r-Zh12/SW*Zn21r*Zn33i-Zh12/SW*Zn23i*Zn31r
     &      -Zh12/SW*Zn23r*Zn31i+Zh22/SW*Zn21i*Zn43r+Zh22
     &      /SW*Zn21r*Zn43i+Zh22/SW*Zn23i*Zn41r+Zh22/SW*Zn23r
     &      *Zn41i+Zh32/CW*Zn11i*Zn33i*sb-Zh32/CW*Zn11i*Zn43i
     &      *cb-Zh32/CW*Zn11r*Zn33r*sb+Zh32/CW*Zn11r*Zn43r
     &      *cb+Zh32/CW*Zn13i*Zn31i*sb-Zh32/CW*Zn13i*Zn41i
     &      *cb-Zh32/CW*Zn13r*Zn31r*sb+Zh32/CW*Zn13r*Zn41r
     &      *cb-Zh32/SW*Zn21i*Zn33i*sb+Zh32/SW*Zn21i*Zn43i
     &      *cb+Zh32/SW*Zn21r*Zn33r*sb-Zh32/SW*Zn21r*Zn43r
     &      *cb-Zh32/SW*Zn23i*Zn31i*sb+Zh32/SW*Zn23i*Zn41i
     &      *cb+Zh32/SW*Zn23r*Zn31r*sb-Zh32/SW*Zn23r*Zn41r
     &      *cb)
      AAABR(958) = EE*(Zh11/CW*Zn11i*Zn31i-Zh11/CW*Zn11r
     &      *Zn31r-Zh21/CW*Zn11i*Zn41i+Zh21/CW*Zn11r*Zn41r
     &      -Zh11/SW*Zn21i*Zn31i+Zh11/SW*Zn21r*Zn31r+Zh21
     &      /SW*Zn21i*Zn41i-Zh21/SW*Zn21r*Zn41r-Zh31/CW*Zn11i
     &      *Zn31r*sb+Zh31/CW*Zn11i*Zn41r*cb-Zh31/CW*Zn11r
     &      *Zn31i*sb+Zh31/CW*Zn11r*Zn41i*cb+Zh31/SW*Zn21i
     &      *Zn31r*sb-Zh31/SW*Zn21i*Zn41r*cb+Zh31/SW*Zn21r
     &      *Zn31i*sb-Zh31/SW*Zn21r*Zn41i*cb)
      AAABR(959) = EE*(Zh11/CW*Zn11i*Zn31r+Zh11/CW*Zn11r
     &      *Zn31i-Zh21/CW*Zn11i*Zn41r-Zh21/CW*Zn11r*Zn41i
     &      -Zh11/SW*Zn21i*Zn31r-Zh11/SW*Zn21r*Zn31i+Zh21
     &      /SW*Zn21i*Zn41r+Zh21/SW*Zn21r*Zn41i+Zh31/CW*Zn11i
     &      *Zn31i*sb-Zh31/CW*Zn11i*Zn41i*cb-Zh31/CW*Zn11r
     &      *Zn31r*sb+Zh31/CW*Zn11r*Zn41r*cb-Zh31/SW*Zn21i
     &      *Zn31i*sb+Zh31/SW*Zn21i*Zn41i*cb+Zh31/SW*Zn21r
     &      *Zn31r*sb-Zh31/SW*Zn21r*Zn41r*cb)
      AAABR(960) = EE*(Zh11/CW*Zn11i*Zn34i-Zh11/CW*Zn11r
     &      *Zn34r+Zh11/CW*Zn14i*Zn31i-Zh11/CW*Zn14r*Zn31r
     &      -Zh21/CW*Zn11i*Zn44i+Zh21/CW*Zn11r*Zn44r-Zh21
     &      /CW*Zn14i*Zn41i+Zh21/CW*Zn14r*Zn41r-Zh11/SW*Zn21i
     &      *Zn34i+Zh11/SW*Zn21r*Zn34r-Zh11/SW*Zn24i*Zn31i
     &      +Zh11/SW*Zn24r*Zn31r+Zh21/SW*Zn21i*Zn44i-Zh21
     &      /SW*Zn21r*Zn44r+Zh21/SW*Zn24i*Zn41i-Zh21/SW*Zn24r
     &      *Zn41r-Zh31/CW*Zn11i*Zn34r*sb+Zh31/CW*Zn11i*Zn44r
     &      *cb-Zh31/CW*Zn11r*Zn34i*sb+Zh31/CW*Zn11r*Zn44i
     &      *cb-Zh31/CW*Zn14i*Zn31r*sb+Zh31/CW*Zn14i*Zn41r
     &      *cb-Zh31/CW*Zn14r*Zn31i*sb+Zh31/CW*Zn14r*Zn41i
     &      *cb+Zh31/SW*Zn21i*Zn34r*sb-Zh31/SW*Zn21i*Zn44r
     &      *cb+Zh31/SW*Zn21r*Zn34i*sb-Zh31/SW*Zn21r*Zn44i
     &      *cb+Zh31/SW*Zn24i*Zn31r*sb-Zh31/SW*Zn24i*Zn41r
     &      *cb+Zh31/SW*Zn24r*Zn31i*sb-Zh31/SW*Zn24r*Zn41i
     &      *cb)
      AAABR(961) = EE*(Zh11/CW*Zn11i*Zn34r+Zh11/CW*Zn11r
     &      *Zn34i+Zh11/CW*Zn14i*Zn31r+Zh11/CW*Zn14r*Zn31i
     &      -Zh21/CW*Zn11i*Zn44r-Zh21/CW*Zn11r*Zn44i-Zh21
     &      /CW*Zn14i*Zn41r-Zh21/CW*Zn14r*Zn41i-Zh11/SW*Zn21i
     &      *Zn34r-Zh11/SW*Zn21r*Zn34i-Zh11/SW*Zn24i*Zn31r
     &      -Zh11/SW*Zn24r*Zn31i+Zh21/SW*Zn21i*Zn44r+Zh21
     &      /SW*Zn21r*Zn44i+Zh21/SW*Zn24i*Zn41r+Zh21/SW*Zn24r
     &      *Zn41i+Zh31/CW*Zn11i*Zn34i*sb-Zh31/CW*Zn11i*Zn44i
     &      *cb-Zh31/CW*Zn11r*Zn34r*sb+Zh31/CW*Zn11r*Zn44r
     &      *cb+Zh31/CW*Zn14i*Zn31i*sb-Zh31/CW*Zn14i*Zn41i
     &      *cb-Zh31/CW*Zn14r*Zn31r*sb+Zh31/CW*Zn14r*Zn41r
     &      *cb-Zh31/SW*Zn21i*Zn34i*sb+Zh31/SW*Zn21i*Zn44i
     &      *cb+Zh31/SW*Zn21r*Zn34r*sb-Zh31/SW*Zn21r*Zn44r
     &      *cb-Zh31/SW*Zn24i*Zn31i*sb+Zh31/SW*Zn24i*Zn41i
     &      *cb+Zh31/SW*Zn24r*Zn31r*sb-Zh31/SW*Zn24r*Zn41r
     &      *cb)
      AAABR(962) = EE*(Zh12/CW*Zn11i*Zn31i-Zh12/CW*Zn11r
     &      *Zn31r-Zh22/CW*Zn11i*Zn41i+Zh22/CW*Zn11r*Zn41r
     &      -Zh12/SW*Zn21i*Zn31i+Zh12/SW*Zn21r*Zn31r+Zh22
     &      /SW*Zn21i*Zn41i-Zh22/SW*Zn21r*Zn41r-Zh32/CW*Zn11i
     &      *Zn31r*sb+Zh32/CW*Zn11i*Zn41r*cb-Zh32/CW*Zn11r
     &      *Zn31i*sb+Zh32/CW*Zn11r*Zn41i*cb+Zh32/SW*Zn21i
     &      *Zn31r*sb-Zh32/SW*Zn21i*Zn41r*cb+Zh32/SW*Zn21r
     &      *Zn31i*sb-Zh32/SW*Zn21r*Zn41i*cb)
      AAABR(963) = EE*(Zh12/CW*Zn11i*Zn31r+Zh12/CW*Zn11r
     &      *Zn31i-Zh22/CW*Zn11i*Zn41r-Zh22/CW*Zn11r*Zn41i
     &      -Zh12/SW*Zn21i*Zn31r-Zh12/SW*Zn21r*Zn31i+Zh22
     &      /SW*Zn21i*Zn41r+Zh22/SW*Zn21r*Zn41i+Zh32/CW*Zn11i
     &      *Zn31i*sb-Zh32/CW*Zn11i*Zn41i*cb-Zh32/CW*Zn11r
     &      *Zn31r*sb+Zh32/CW*Zn11r*Zn41r*cb-Zh32/SW*Zn21i
     &      *Zn31i*sb+Zh32/SW*Zn21i*Zn41i*cb+Zh32/SW*Zn21r
     &      *Zn31r*sb-Zh32/SW*Zn21r*Zn41r*cb)
      AAABR(964) = EE*(Zn12i/CW*Zn33r*cb+Zn12i/CW*Zn43r*sb
     &      +Zn12r/CW*Zn33i*cb+Zn12r/CW*Zn43i*sb+Zn13i/CW
     &      *Zn32r*cb+Zn13i/CW*Zn42r*sb+Zn13r/CW*Zn32i*cb
     &      +Zn13r/CW*Zn42i*sb-Zn22i/SW*Zn33r*cb-Zn22i/SW
     &      *Zn43r*sb-Zn22r/SW*Zn33i*cb-Zn22r/SW*Zn43i*sb
     &      -Zn23i/SW*Zn32r*cb-Zn23i/SW*Zn42r*sb-Zn23r/SW
     &      *Zn32i*cb-Zn23r/SW*Zn42i*sb)
      AAABR(965) = EE*(Zn12i/CW*Zn33i*cb+Zn12i/CW*Zn43i*sb
     &      -Zn12r/CW*Zn33r*cb-Zn12r/CW*Zn43r*sb+Zn13i/CW
     &      *Zn32i*cb+Zn13i/CW*Zn42i*sb-Zn13r/CW*Zn32r*cb
     &      -Zn13r/CW*Zn42r*sb-Zn22i/SW*Zn33i*cb-Zn22i/SW
     &      *Zn43i*sb+Zn22r/SW*Zn33r*cb+Zn22r/SW*Zn43r*sb
     &      -Zn23i/SW*Zn32i*cb-Zn23i/SW*Zn42i*sb+Zn23r/SW
     &      *Zn32r*cb+Zn23r/SW*Zn42r*sb)
      AAABR(966) = EE*(Zn12i/CW*Zn34r*cb+Zn12i/CW*Zn44r*sb
     &      +Zn12r/CW*Zn34i*cb+Zn12r/CW*Zn44i*sb+Zn14i/CW
     &      *Zn32r*cb+Zn14i/CW*Zn42r*sb+Zn14r/CW*Zn32i*cb
     &      +Zn14r/CW*Zn42i*sb-Zn22i/SW*Zn34r*cb-Zn22i/SW
     &      *Zn44r*sb-Zn22r/SW*Zn34i*cb-Zn22r/SW*Zn44i*sb
     &      -Zn24i/SW*Zn32r*cb-Zn24i/SW*Zn42r*sb-Zn24r/SW
     &      *Zn32i*cb-Zn24r/SW*Zn42i*sb)
      AAABR(967) = EE*(Zn12i/CW*Zn34i*cb+Zn12i/CW*Zn44i*sb
     &      -Zn12r/CW*Zn34r*cb-Zn12r/CW*Zn44r*sb+Zn14i/CW
     &      *Zn32i*cb+Zn14i/CW*Zn42i*sb-Zn14r/CW*Zn32r*cb
     &      -Zn14r/CW*Zn42r*sb-Zn22i/SW*Zn34i*cb-Zn22i/SW
     &      *Zn44i*sb+Zn22r/SW*Zn34r*cb+Zn22r/SW*Zn44r*sb
     &      -Zn24i/SW*Zn32i*cb-Zn24i/SW*Zn42i*sb+Zn24r/SW
     &      *Zn32r*cb+Zn24r/SW*Zn42r*sb)
      AAABR(968) = EE*(Zn12i/CW*Zn32r*cb+Zn12i/CW*Zn42r*sb
     &      +Zn12r/CW*Zn32i*cb+Zn12r/CW*Zn42i*sb-Zn22i/SW
     &      *Zn32r*cb-Zn22i/SW*Zn42r*sb-Zn22r/SW*Zn32i*cb
     &      -Zn22r/SW*Zn42i*sb)
      AAABR(969) = EE*(Zn12i/CW*Zn32i*cb+Zn12i/CW*Zn42i*sb
     &      -Zn12r/CW*Zn32r*cb-Zn12r/CW*Zn42r*sb-Zn22i/SW
     &      *Zn32i*cb-Zn22i/SW*Zn42i*sb+Zn22r/SW*Zn32r*cb
     &      +Zn22r/SW*Zn42r*sb)
      AAABR(970) = EE*(Zh13/CW*Zn12i*Zn34i-Zh13/CW*Zn12r
     &      *Zn34r+Zh13/CW*Zn14i*Zn32i-Zh13/CW*Zn14r*Zn32r
     &      -Zh23/CW*Zn12i*Zn44i+Zh23/CW*Zn12r*Zn44r-Zh23
     &      /CW*Zn14i*Zn42i+Zh23/CW*Zn14r*Zn42r-Zh13/SW*Zn22i
     &      *Zn34i+Zh13/SW*Zn22r*Zn34r-Zh13/SW*Zn24i*Zn32i
     &      +Zh13/SW*Zn24r*Zn32r+Zh23/SW*Zn22i*Zn44i-Zh23
     &      /SW*Zn22r*Zn44r+Zh23/SW*Zn24i*Zn42i-Zh23/SW*Zn24r
     &      *Zn42r-Zh33/CW*Zn12i*Zn34r*sb+Zh33/CW*Zn12i*Zn44r
     &      *cb-Zh33/CW*Zn12r*Zn34i*sb+Zh33/CW*Zn12r*Zn44i
     &      *cb-Zh33/CW*Zn14i*Zn32r*sb+Zh33/CW*Zn14i*Zn42r
     &      *cb-Zh33/CW*Zn14r*Zn32i*sb+Zh33/CW*Zn14r*Zn42i
     &      *cb+Zh33/SW*Zn22i*Zn34r*sb-Zh33/SW*Zn22i*Zn44r
     &      *cb+Zh33/SW*Zn22r*Zn34i*sb-Zh33/SW*Zn22r*Zn44i
     &      *cb+Zh33/SW*Zn24i*Zn32r*sb-Zh33/SW*Zn24i*Zn42r
     &      *cb+Zh33/SW*Zn24r*Zn32i*sb-Zh33/SW*Zn24r*Zn42i
     &      *cb)
      AAABR(971) = EE*(Zh13/CW*Zn12i*Zn34r+Zh13/CW*Zn12r
     &      *Zn34i+Zh13/CW*Zn14i*Zn32r+Zh13/CW*Zn14r*Zn32i
     &      -Zh23/CW*Zn12i*Zn44r-Zh23/CW*Zn12r*Zn44i-Zh23
     &      /CW*Zn14i*Zn42r-Zh23/CW*Zn14r*Zn42i-Zh13/SW*Zn22i
     &      *Zn34r-Zh13/SW*Zn22r*Zn34i-Zh13/SW*Zn24i*Zn32r
     &      -Zh13/SW*Zn24r*Zn32i+Zh23/SW*Zn22i*Zn44r+Zh23
     &      /SW*Zn22r*Zn44i+Zh23/SW*Zn24i*Zn42r+Zh23/SW*Zn24r
     &      *Zn42i+Zh33/CW*Zn12i*Zn34i*sb-Zh33/CW*Zn12i*Zn44i
     &      *cb-Zh33/CW*Zn12r*Zn34r*sb+Zh33/CW*Zn12r*Zn44r
     &      *cb+Zh33/CW*Zn14i*Zn32i*sb-Zh33/CW*Zn14i*Zn42i
     &      *cb-Zh33/CW*Zn14r*Zn32r*sb+Zh33/CW*Zn14r*Zn42r
     &      *cb-Zh33/SW*Zn22i*Zn34i*sb+Zh33/SW*Zn22i*Zn44i
     &      *cb+Zh33/SW*Zn22r*Zn34r*sb-Zh33/SW*Zn22r*Zn44r
     &      *cb-Zh33/SW*Zn24i*Zn32i*sb+Zh33/SW*Zn24i*Zn42i
     &      *cb+Zh33/SW*Zn24r*Zn32r*sb-Zh33/SW*Zn24r*Zn42r
     &      *cb)
      AAABR(972) = EE*(Zh12/CW*Zn12i*Zn33i-Zh12/CW*Zn12r
     &      *Zn33r+Zh12/CW*Zn13i*Zn32i-Zh12/CW*Zn13r*Zn32r
     &      -Zh22/CW*Zn12i*Zn43i+Zh22/CW*Zn12r*Zn43r-Zh22
     &      /CW*Zn13i*Zn42i+Zh22/CW*Zn13r*Zn42r-Zh12/SW*Zn22i
     &      *Zn33i+Zh12/SW*Zn22r*Zn33r-Zh12/SW*Zn23i*Zn32i
     &      +Zh12/SW*Zn23r*Zn32r+Zh22/SW*Zn22i*Zn43i-Zh22
     &      /SW*Zn22r*Zn43r+Zh22/SW*Zn23i*Zn42i-Zh22/SW*Zn23r
     &      *Zn42r-Zh32/CW*Zn12i*Zn33r*sb+Zh32/CW*Zn12i*Zn43r
     &      *cb-Zh32/CW*Zn12r*Zn33i*sb+Zh32/CW*Zn12r*Zn43i
     &      *cb-Zh32/CW*Zn13i*Zn32r*sb+Zh32/CW*Zn13i*Zn42r
     &      *cb-Zh32/CW*Zn13r*Zn32i*sb+Zh32/CW*Zn13r*Zn42i
     &      *cb+Zh32/SW*Zn22i*Zn33r*sb-Zh32/SW*Zn22i*Zn43r
     &      *cb+Zh32/SW*Zn22r*Zn33i*sb-Zh32/SW*Zn22r*Zn43i
     &      *cb+Zh32/SW*Zn23i*Zn32r*sb-Zh32/SW*Zn23i*Zn42r
     &      *cb+Zh32/SW*Zn23r*Zn32i*sb-Zh32/SW*Zn23r*Zn42i
     &      *cb)
      AAABR(973) = EE*(Zh12/CW*Zn12i*Zn33r+Zh12/CW*Zn12r
     &      *Zn33i+Zh12/CW*Zn13i*Zn32r+Zh12/CW*Zn13r*Zn32i
     &      -Zh22/CW*Zn12i*Zn43r-Zh22/CW*Zn12r*Zn43i-Zh22
     &      /CW*Zn13i*Zn42r-Zh22/CW*Zn13r*Zn42i-Zh12/SW*Zn22i
     &      *Zn33r-Zh12/SW*Zn22r*Zn33i-Zh12/SW*Zn23i*Zn32r
     &      -Zh12/SW*Zn23r*Zn32i+Zh22/SW*Zn22i*Zn43r+Zh22
     &      /SW*Zn22r*Zn43i+Zh22/SW*Zn23i*Zn42r+Zh22/SW*Zn23r
     &      *Zn42i+Zh32/CW*Zn12i*Zn33i*sb-Zh32/CW*Zn12i*Zn43i
     &      *cb-Zh32/CW*Zn12r*Zn33r*sb+Zh32/CW*Zn12r*Zn43r
     &      *cb+Zh32/CW*Zn13i*Zn32i*sb-Zh32/CW*Zn13i*Zn42i
     &      *cb-Zh32/CW*Zn13r*Zn32r*sb+Zh32/CW*Zn13r*Zn42r
     &      *cb-Zh32/SW*Zn22i*Zn33i*sb+Zh32/SW*Zn22i*Zn43i
     &      *cb+Zh32/SW*Zn22r*Zn33r*sb-Zh32/SW*Zn22r*Zn43r
     &      *cb-Zh32/SW*Zn23i*Zn32i*sb+Zh32/SW*Zn23i*Zn42i
     &      *cb+Zh32/SW*Zn23r*Zn32r*sb-Zh32/SW*Zn23r*Zn42r
     &      *cb)
      AAABR(974) = EE*(Zh11/CW*Zn12i*Zn34i-Zh11/CW*Zn12r
     &      *Zn34r+Zh11/CW*Zn14i*Zn32i-Zh11/CW*Zn14r*Zn32r
     &      -Zh21/CW*Zn12i*Zn44i+Zh21/CW*Zn12r*Zn44r-Zh21
     &      /CW*Zn14i*Zn42i+Zh21/CW*Zn14r*Zn42r-Zh11/SW*Zn22i
     &      *Zn34i+Zh11/SW*Zn22r*Zn34r-Zh11/SW*Zn24i*Zn32i
     &      +Zh11/SW*Zn24r*Zn32r+Zh21/SW*Zn22i*Zn44i-Zh21
     &      /SW*Zn22r*Zn44r+Zh21/SW*Zn24i*Zn42i-Zh21/SW*Zn24r
     &      *Zn42r-Zh31/CW*Zn12i*Zn34r*sb+Zh31/CW*Zn12i*Zn44r
     &      *cb-Zh31/CW*Zn12r*Zn34i*sb+Zh31/CW*Zn12r*Zn44i
     &      *cb-Zh31/CW*Zn14i*Zn32r*sb+Zh31/CW*Zn14i*Zn42r
     &      *cb-Zh31/CW*Zn14r*Zn32i*sb+Zh31/CW*Zn14r*Zn42i
     &      *cb+Zh31/SW*Zn22i*Zn34r*sb-Zh31/SW*Zn22i*Zn44r
     &      *cb+Zh31/SW*Zn22r*Zn34i*sb-Zh31/SW*Zn22r*Zn44i
     &      *cb+Zh31/SW*Zn24i*Zn32r*sb-Zh31/SW*Zn24i*Zn42r
     &      *cb+Zh31/SW*Zn24r*Zn32i*sb-Zh31/SW*Zn24r*Zn42i
     &      *cb)
      AAABR(975) = EE*(Zh11/CW*Zn12i*Zn34r+Zh11/CW*Zn12r
     &      *Zn34i+Zh11/CW*Zn14i*Zn32r+Zh11/CW*Zn14r*Zn32i
     &      -Zh21/CW*Zn12i*Zn44r-Zh21/CW*Zn12r*Zn44i-Zh21
     &      /CW*Zn14i*Zn42r-Zh21/CW*Zn14r*Zn42i-Zh11/SW*Zn22i
     &      *Zn34r-Zh11/SW*Zn22r*Zn34i-Zh11/SW*Zn24i*Zn32r
     &      -Zh11/SW*Zn24r*Zn32i+Zh21/SW*Zn22i*Zn44r+Zh21
     &      /SW*Zn22r*Zn44i+Zh21/SW*Zn24i*Zn42r+Zh21/SW*Zn24r
     &      *Zn42i+Zh31/CW*Zn12i*Zn34i*sb-Zh31/CW*Zn12i*Zn44i
     &      *cb-Zh31/CW*Zn12r*Zn34r*sb+Zh31/CW*Zn12r*Zn44r
     &      *cb+Zh31/CW*Zn14i*Zn32i*sb-Zh31/CW*Zn14i*Zn42i
     &      *cb-Zh31/CW*Zn14r*Zn32r*sb+Zh31/CW*Zn14r*Zn42r
     &      *cb-Zh31/SW*Zn22i*Zn34i*sb+Zh31/SW*Zn22i*Zn44i
     &      *cb+Zh31/SW*Zn22r*Zn34r*sb-Zh31/SW*Zn22r*Zn44r
     &      *cb-Zh31/SW*Zn24i*Zn32i*sb+Zh31/SW*Zn24i*Zn42i
     &      *cb+Zh31/SW*Zn24r*Zn32r*sb-Zh31/SW*Zn24r*Zn42r
     &      *cb)
      AAABR(976) = EE*(Zh12/CW*Zn12i*Zn34i-Zh12/CW*Zn12r
     &      *Zn34r+Zh12/CW*Zn14i*Zn32i-Zh12/CW*Zn14r*Zn32r
     &      -Zh22/CW*Zn12i*Zn44i+Zh22/CW*Zn12r*Zn44r-Zh22
     &      /CW*Zn14i*Zn42i+Zh22/CW*Zn14r*Zn42r-Zh12/SW*Zn22i
     &      *Zn34i+Zh12/SW*Zn22r*Zn34r-Zh12/SW*Zn24i*Zn32i
     &      +Zh12/SW*Zn24r*Zn32r+Zh22/SW*Zn22i*Zn44i-Zh22
     &      /SW*Zn22r*Zn44r+Zh22/SW*Zn24i*Zn42i-Zh22/SW*Zn24r
     &      *Zn42r-Zh32/CW*Zn12i*Zn34r*sb+Zh32/CW*Zn12i*Zn44r
     &      *cb-Zh32/CW*Zn12r*Zn34i*sb+Zh32/CW*Zn12r*Zn44i
     &      *cb-Zh32/CW*Zn14i*Zn32r*sb+Zh32/CW*Zn14i*Zn42r
     &      *cb-Zh32/CW*Zn14r*Zn32i*sb+Zh32/CW*Zn14r*Zn42i
     &      *cb+Zh32/SW*Zn22i*Zn34r*sb-Zh32/SW*Zn22i*Zn44r
     &      *cb+Zh32/SW*Zn22r*Zn34i*sb-Zh32/SW*Zn22r*Zn44i
     &      *cb+Zh32/SW*Zn24i*Zn32r*sb-Zh32/SW*Zn24i*Zn42r
     &      *cb+Zh32/SW*Zn24r*Zn32i*sb-Zh32/SW*Zn24r*Zn42i
     &      *cb)
      AAABR(977) = EE*(Zh12/CW*Zn12i*Zn34r+Zh12/CW*Zn12r
     &      *Zn34i+Zh12/CW*Zn14i*Zn32r+Zh12/CW*Zn14r*Zn32i
     &      -Zh22/CW*Zn12i*Zn44r-Zh22/CW*Zn12r*Zn44i-Zh22
     &      /CW*Zn14i*Zn42r-Zh22/CW*Zn14r*Zn42i-Zh12/SW*Zn22i
     &      *Zn34r-Zh12/SW*Zn22r*Zn34i-Zh12/SW*Zn24i*Zn32r
     &      -Zh12/SW*Zn24r*Zn32i+Zh22/SW*Zn22i*Zn44r+Zh22
     &      /SW*Zn22r*Zn44i+Zh22/SW*Zn24i*Zn42r+Zh22/SW*Zn24r
     &      *Zn42i+Zh32/CW*Zn12i*Zn34i*sb-Zh32/CW*Zn12i*Zn44i
     &      *cb-Zh32/CW*Zn12r*Zn34r*sb+Zh32/CW*Zn12r*Zn44r
     &      *cb+Zh32/CW*Zn14i*Zn32i*sb-Zh32/CW*Zn14i*Zn42i
     &      *cb-Zh32/CW*Zn14r*Zn32r*sb+Zh32/CW*Zn14r*Zn42r
     &      *cb-Zh32/SW*Zn22i*Zn34i*sb+Zh32/SW*Zn22i*Zn44i
     &      *cb+Zh32/SW*Zn22r*Zn34r*sb-Zh32/SW*Zn22r*Zn44r
     &      *cb-Zh32/SW*Zn24i*Zn32i*sb+Zh32/SW*Zn24i*Zn42i
     &      *cb+Zh32/SW*Zn24r*Zn32r*sb-Zh32/SW*Zn24r*Zn42r
     &      *cb)
      AAABR(978) = EE*(Zh11/CW*Zn12i*Zn33i-Zh11/CW*Zn12r
     &      *Zn33r+Zh11/CW*Zn13i*Zn32i-Zh11/CW*Zn13r*Zn32r
     &      -Zh21/CW*Zn12i*Zn43i+Zh21/CW*Zn12r*Zn43r-Zh21
     &      /CW*Zn13i*Zn42i+Zh21/CW*Zn13r*Zn42r-Zh11/SW*Zn22i
     &      *Zn33i+Zh11/SW*Zn22r*Zn33r-Zh11/SW*Zn23i*Zn32i
     &      +Zh11/SW*Zn23r*Zn32r+Zh21/SW*Zn22i*Zn43i-Zh21
     &      /SW*Zn22r*Zn43r+Zh21/SW*Zn23i*Zn42i-Zh21/SW*Zn23r
     &      *Zn42r-Zh31/CW*Zn12i*Zn33r*sb+Zh31/CW*Zn12i*Zn43r
     &      *cb-Zh31/CW*Zn12r*Zn33i*sb+Zh31/CW*Zn12r*Zn43i
     &      *cb-Zh31/CW*Zn13i*Zn32r*sb+Zh31/CW*Zn13i*Zn42r
     &      *cb-Zh31/CW*Zn13r*Zn32i*sb+Zh31/CW*Zn13r*Zn42i
     &      *cb+Zh31/SW*Zn22i*Zn33r*sb-Zh31/SW*Zn22i*Zn43r
     &      *cb+Zh31/SW*Zn22r*Zn33i*sb-Zh31/SW*Zn22r*Zn43i
     &      *cb+Zh31/SW*Zn23i*Zn32r*sb-Zh31/SW*Zn23i*Zn42r
     &      *cb+Zh31/SW*Zn23r*Zn32i*sb-Zh31/SW*Zn23r*Zn42i
     &      *cb)
      AAABR(979) = EE*(Zh11/CW*Zn12i*Zn33r+Zh11/CW*Zn12r
     &      *Zn33i+Zh11/CW*Zn13i*Zn32r+Zh11/CW*Zn13r*Zn32i
     &      -Zh21/CW*Zn12i*Zn43r-Zh21/CW*Zn12r*Zn43i-Zh21
     &      /CW*Zn13i*Zn42r-Zh21/CW*Zn13r*Zn42i-Zh11/SW*Zn22i
     &      *Zn33r-Zh11/SW*Zn22r*Zn33i-Zh11/SW*Zn23i*Zn32r
     &      -Zh11/SW*Zn23r*Zn32i+Zh21/SW*Zn22i*Zn43r+Zh21
     &      /SW*Zn22r*Zn43i+Zh21/SW*Zn23i*Zn42r+Zh21/SW*Zn23r
     &      *Zn42i+Zh31/CW*Zn12i*Zn33i*sb-Zh31/CW*Zn12i*Zn43i
     &      *cb-Zh31/CW*Zn12r*Zn33r*sb+Zh31/CW*Zn12r*Zn43r
     &      *cb+Zh31/CW*Zn13i*Zn32i*sb-Zh31/CW*Zn13i*Zn42i
     &      *cb-Zh31/CW*Zn13r*Zn32r*sb+Zh31/CW*Zn13r*Zn42r
     &      *cb-Zh31/SW*Zn22i*Zn33i*sb+Zh31/SW*Zn22i*Zn43i
     &      *cb+Zh31/SW*Zn22r*Zn33r*sb-Zh31/SW*Zn22r*Zn43r
     &      *cb-Zh31/SW*Zn23i*Zn32i*sb+Zh31/SW*Zn23i*Zn42i
     &      *cb+Zh31/SW*Zn23r*Zn32r*sb-Zh31/SW*Zn23r*Zn42r
     &      *cb)
      AAABR(980) = EE*(Zh13/CW*Zn12i*Zn32i-Zh13/CW*Zn12r
     &      *Zn32r-Zh23/CW*Zn12i*Zn42i+Zh23/CW*Zn12r*Zn42r
     &      -Zh13/SW*Zn22i*Zn32i+Zh13/SW*Zn22r*Zn32r+Zh23
     &      /SW*Zn22i*Zn42i-Zh23/SW*Zn22r*Zn42r-Zh33/CW*Zn12i
     &      *Zn32r*sb+Zh33/CW*Zn12i*Zn42r*cb-Zh33/CW*Zn12r
     &      *Zn32i*sb+Zh33/CW*Zn12r*Zn42i*cb+Zh33/SW*Zn22i
     &      *Zn32r*sb-Zh33/SW*Zn22i*Zn42r*cb+Zh33/SW*Zn22r
     &      *Zn32i*sb-Zh33/SW*Zn22r*Zn42i*cb)
      AAABR(981) = EE*(Zh13/CW*Zn12i*Zn32r+Zh13/CW*Zn12r
     &      *Zn32i-Zh23/CW*Zn12i*Zn42r-Zh23/CW*Zn12r*Zn42i
     &      -Zh13/SW*Zn22i*Zn32r-Zh13/SW*Zn22r*Zn32i+Zh23
     &      /SW*Zn22i*Zn42r+Zh23/SW*Zn22r*Zn42i+Zh33/CW*Zn12i
     &      *Zn32i*sb-Zh33/CW*Zn12i*Zn42i*cb-Zh33/CW*Zn12r
     &      *Zn32r*sb+Zh33/CW*Zn12r*Zn42r*cb-Zh33/SW*Zn22i
     &      *Zn32i*sb+Zh33/SW*Zn22i*Zn42i*cb+Zh33/SW*Zn22r
     &      *Zn32r*sb-Zh33/SW*Zn22r*Zn42r*cb)
      AAABR(982) = EE*(Zh12/CW*Zn12i*Zn32i-Zh12/CW*Zn12r
     &      *Zn32r-Zh22/CW*Zn12i*Zn42i+Zh22/CW*Zn12r*Zn42r
     &      -Zh12/SW*Zn22i*Zn32i+Zh12/SW*Zn22r*Zn32r+Zh22
     &      /SW*Zn22i*Zn42i-Zh22/SW*Zn22r*Zn42r-Zh32/CW*Zn12i
     &      *Zn32r*sb+Zh32/CW*Zn12i*Zn42r*cb-Zh32/CW*Zn12r
     &      *Zn32i*sb+Zh32/CW*Zn12r*Zn42i*cb+Zh32/SW*Zn22i
     &      *Zn32r*sb-Zh32/SW*Zn22i*Zn42r*cb+Zh32/SW*Zn22r
     &      *Zn32i*sb-Zh32/SW*Zn22r*Zn42i*cb)
      AAABR(983) = EE*(Zh12/CW*Zn12i*Zn32r+Zh12/CW*Zn12r
     &      *Zn32i-Zh22/CW*Zn12i*Zn42r-Zh22/CW*Zn12r*Zn42i
     &      -Zh12/SW*Zn22i*Zn32r-Zh12/SW*Zn22r*Zn32i+Zh22
     &      /SW*Zn22i*Zn42r+Zh22/SW*Zn22r*Zn42i+Zh32/CW*Zn12i
     &      *Zn32i*sb-Zh32/CW*Zn12i*Zn42i*cb-Zh32/CW*Zn12r
     &      *Zn32r*sb+Zh32/CW*Zn12r*Zn42r*cb-Zh32/SW*Zn22i
     &      *Zn32i*sb+Zh32/SW*Zn22i*Zn42i*cb+Zh32/SW*Zn22r
     &      *Zn32r*sb-Zh32/SW*Zn22r*Zn42r*cb)
      AAABR(984) = EE*(Zh11/CW*Zn12i*Zn32i-Zh11/CW*Zn12r
     &      *Zn32r-Zh21/CW*Zn12i*Zn42i+Zh21/CW*Zn12r*Zn42r
     &      -Zh11/SW*Zn22i*Zn32i+Zh11/SW*Zn22r*Zn32r+Zh21
     &      /SW*Zn22i*Zn42i-Zh21/SW*Zn22r*Zn42r-Zh31/CW*Zn12i
     &      *Zn32r*sb+Zh31/CW*Zn12i*Zn42r*cb-Zh31/CW*Zn12r
     &      *Zn32i*sb+Zh31/CW*Zn12r*Zn42i*cb+Zh31/SW*Zn22i
     &      *Zn32r*sb-Zh31/SW*Zn22i*Zn42r*cb+Zh31/SW*Zn22r
     &      *Zn32i*sb-Zh31/SW*Zn22r*Zn42i*cb)
      AAABR(985) = EE*(Zh11/CW*Zn12i*Zn32r+Zh11/CW*Zn12r
     &      *Zn32i-Zh21/CW*Zn12i*Zn42r-Zh21/CW*Zn12r*Zn42i
     &      -Zh11/SW*Zn22i*Zn32r-Zh11/SW*Zn22r*Zn32i+Zh21
     &      /SW*Zn22i*Zn42r+Zh21/SW*Zn22r*Zn42i+Zh31/CW*Zn12i
     &      *Zn32i*sb-Zh31/CW*Zn12i*Zn42i*cb-Zh31/CW*Zn12r
     &      *Zn32r*sb+Zh31/CW*Zn12r*Zn42r*cb-Zh31/SW*Zn22i
     &      *Zn32i*sb+Zh31/SW*Zn22i*Zn42i*cb+Zh31/SW*Zn22r
     &      *Zn32r*sb-Zh31/SW*Zn22r*Zn42r*cb)
      AAABR(986) = EE*(Zh13/CW*Zn12i*Zn33i-Zh13/CW*Zn12r
     &      *Zn33r+Zh13/CW*Zn13i*Zn32i-Zh13/CW*Zn13r*Zn32r
     &      -Zh23/CW*Zn12i*Zn43i+Zh23/CW*Zn12r*Zn43r-Zh23
     &      /CW*Zn13i*Zn42i+Zh23/CW*Zn13r*Zn42r-Zh13/SW*Zn22i
     &      *Zn33i+Zh13/SW*Zn22r*Zn33r-Zh13/SW*Zn23i*Zn32i
     &      +Zh13/SW*Zn23r*Zn32r+Zh23/SW*Zn22i*Zn43i-Zh23
     &      /SW*Zn22r*Zn43r+Zh23/SW*Zn23i*Zn42i-Zh23/SW*Zn23r
     &      *Zn42r-Zh33/CW*Zn12i*Zn33r*sb+Zh33/CW*Zn12i*Zn43r
     &      *cb-Zh33/CW*Zn12r*Zn33i*sb+Zh33/CW*Zn12r*Zn43i
     &      *cb-Zh33/CW*Zn13i*Zn32r*sb+Zh33/CW*Zn13i*Zn42r
     &      *cb-Zh33/CW*Zn13r*Zn32i*sb+Zh33/CW*Zn13r*Zn42i
     &      *cb+Zh33/SW*Zn22i*Zn33r*sb-Zh33/SW*Zn22i*Zn43r
     &      *cb+Zh33/SW*Zn22r*Zn33i*sb-Zh33/SW*Zn22r*Zn43i
     &      *cb+Zh33/SW*Zn23i*Zn32r*sb-Zh33/SW*Zn23i*Zn42r
     &      *cb+Zh33/SW*Zn23r*Zn32i*sb-Zh33/SW*Zn23r*Zn42i
     &      *cb)
      AAABR(987) = EE*(Zh13/CW*Zn12i*Zn33r+Zh13/CW*Zn12r
     &      *Zn33i+Zh13/CW*Zn13i*Zn32r+Zh13/CW*Zn13r*Zn32i
     &      -Zh23/CW*Zn12i*Zn43r-Zh23/CW*Zn12r*Zn43i-Zh23
     &      /CW*Zn13i*Zn42r-Zh23/CW*Zn13r*Zn42i-Zh13/SW*Zn22i
     &      *Zn33r-Zh13/SW*Zn22r*Zn33i-Zh13/SW*Zn23i*Zn32r
     &      -Zh13/SW*Zn23r*Zn32i+Zh23/SW*Zn22i*Zn43r+Zh23
     &      /SW*Zn22r*Zn43i+Zh23/SW*Zn23i*Zn42r+Zh23/SW*Zn23r
     &      *Zn42i+Zh33/CW*Zn12i*Zn33i*sb-Zh33/CW*Zn12i*Zn43i
     &      *cb-Zh33/CW*Zn12r*Zn33r*sb+Zh33/CW*Zn12r*Zn43r
     &      *cb+Zh33/CW*Zn13i*Zn32i*sb-Zh33/CW*Zn13i*Zn42i
     &      *cb-Zh33/CW*Zn13r*Zn32r*sb+Zh33/CW*Zn13r*Zn42r
     &      *cb-Zh33/SW*Zn22i*Zn33i*sb+Zh33/SW*Zn22i*Zn43i
     &      *cb+Zh33/SW*Zn22r*Zn33r*sb-Zh33/SW*Zn22r*Zn43r
     &      *cb-Zh33/SW*Zn23i*Zn32i*sb+Zh33/SW*Zn23i*Zn42i
     &      *cb+Zh33/SW*Zn23r*Zn32r*sb-Zh33/SW*Zn23r*Zn42r
     &      *cb)
      AAABR(988) = EE*(Zn13i/CW*Zn33r*cb+Zn13i/CW*Zn43r*sb
     &      +Zn13r/CW*Zn33i*cb+Zn13r/CW*Zn43i*sb-Zn23i/SW
     &      *Zn33r*cb-Zn23i/SW*Zn43r*sb-Zn23r/SW*Zn33i*cb
     &      -Zn23r/SW*Zn43i*sb)
      AAABR(989) = EE*(Zn13i/CW*Zn33i*cb+Zn13i/CW*Zn43i*sb
     &      -Zn13r/CW*Zn33r*cb-Zn13r/CW*Zn43r*sb-Zn23i/SW
     &      *Zn33i*cb-Zn23i/SW*Zn43i*sb+Zn23r/SW*Zn33r*cb
     &      +Zn23r/SW*Zn43r*sb)
      AAABR(990) = EE*(Zn13i/CW*Zn34r*cb+Zn13i/CW*Zn44r*sb
     &      +Zn13r/CW*Zn34i*cb+Zn13r/CW*Zn44i*sb+Zn14i/CW
     &      *Zn33r*cb+Zn14i/CW*Zn43r*sb+Zn14r/CW*Zn33i*cb
     &      +Zn14r/CW*Zn43i*sb-Zn23i/SW*Zn34r*cb-Zn23i/SW
     &      *Zn44r*sb-Zn23r/SW*Zn34i*cb-Zn23r/SW*Zn44i*sb
     &      -Zn24i/SW*Zn33r*cb-Zn24i/SW*Zn43r*sb-Zn24r/SW
     &      *Zn33i*cb-Zn24r/SW*Zn43i*sb)
      AAABR(991) = EE*(Zn13i/CW*Zn34i*cb+Zn13i/CW*Zn44i*sb
     &      -Zn13r/CW*Zn34r*cb-Zn13r/CW*Zn44r*sb+Zn14i/CW
     &      *Zn33i*cb+Zn14i/CW*Zn43i*sb-Zn14r/CW*Zn33r*cb
     &      -Zn14r/CW*Zn43r*sb-Zn23i/SW*Zn34i*cb-Zn23i/SW
     &      *Zn44i*sb+Zn23r/SW*Zn34r*cb+Zn23r/SW*Zn44r*sb
     &      -Zn24i/SW*Zn33i*cb-Zn24i/SW*Zn43i*sb+Zn24r/SW
     &      *Zn33r*cb+Zn24r/SW*Zn43r*sb)
      AAABR(992) = EE*(Zh12/CW*Zn13i*Zn33i-Zh12/CW*Zn13r
     &      *Zn33r-Zh22/CW*Zn13i*Zn43i+Zh22/CW*Zn13r*Zn43r
     &      -Zh12/SW*Zn23i*Zn33i+Zh12/SW*Zn23r*Zn33r+Zh22
     &      /SW*Zn23i*Zn43i-Zh22/SW*Zn23r*Zn43r-Zh32/CW*Zn13i
     &      *Zn33r*sb+Zh32/CW*Zn13i*Zn43r*cb-Zh32/CW*Zn13r
     &      *Zn33i*sb+Zh32/CW*Zn13r*Zn43i*cb+Zh32/SW*Zn23i
     &      *Zn33r*sb-Zh32/SW*Zn23i*Zn43r*cb+Zh32/SW*Zn23r
     &      *Zn33i*sb-Zh32/SW*Zn23r*Zn43i*cb)
      AAABR(993) = EE*(Zh12/CW*Zn13i*Zn33r+Zh12/CW*Zn13r
     &      *Zn33i-Zh22/CW*Zn13i*Zn43r-Zh22/CW*Zn13r*Zn43i
     &      -Zh12/SW*Zn23i*Zn33r-Zh12/SW*Zn23r*Zn33i+Zh22
     &      /SW*Zn23i*Zn43r+Zh22/SW*Zn23r*Zn43i+Zh32/CW*Zn13i
     &      *Zn33i*sb-Zh32/CW*Zn13i*Zn43i*cb-Zh32/CW*Zn13r
     &      *Zn33r*sb+Zh32/CW*Zn13r*Zn43r*cb-Zh32/SW*Zn23i
     &      *Zn33i*sb+Zh32/SW*Zn23i*Zn43i*cb+Zh32/SW*Zn23r
     &      *Zn33r*sb-Zh32/SW*Zn23r*Zn43r*cb)
      AAABR(994) = EE*(Zh11/CW*Zn13i*Zn33i-Zh11/CW*Zn13r
     &      *Zn33r-Zh21/CW*Zn13i*Zn43i+Zh21/CW*Zn13r*Zn43r
     &      -Zh11/SW*Zn23i*Zn33i+Zh11/SW*Zn23r*Zn33r+Zh21
     &      /SW*Zn23i*Zn43i-Zh21/SW*Zn23r*Zn43r-Zh31/CW*Zn13i
     &      *Zn33r*sb+Zh31/CW*Zn13i*Zn43r*cb-Zh31/CW*Zn13r
     &      *Zn33i*sb+Zh31/CW*Zn13r*Zn43i*cb+Zh31/SW*Zn23i
     &      *Zn33r*sb-Zh31/SW*Zn23i*Zn43r*cb+Zh31/SW*Zn23r
     &      *Zn33i*sb-Zh31/SW*Zn23r*Zn43i*cb)
      AAABR(995) = EE*(Zh11/CW*Zn13i*Zn33r+Zh11/CW*Zn13r
     &      *Zn33i-Zh21/CW*Zn13i*Zn43r-Zh21/CW*Zn13r*Zn43i
     &      -Zh11/SW*Zn23i*Zn33r-Zh11/SW*Zn23r*Zn33i+Zh21
     &      /SW*Zn23i*Zn43r+Zh21/SW*Zn23r*Zn43i+Zh31/CW*Zn13i
     &      *Zn33i*sb-Zh31/CW*Zn13i*Zn43i*cb-Zh31/CW*Zn13r
     &      *Zn33r*sb+Zh31/CW*Zn13r*Zn43r*cb-Zh31/SW*Zn23i
     &      *Zn33i*sb+Zh31/SW*Zn23i*Zn43i*cb+Zh31/SW*Zn23r
     &      *Zn33r*sb-Zh31/SW*Zn23r*Zn43r*cb)
      AAABR(996) = EE*(Zh13/CW*Zn13i*Zn34i-Zh13/CW*Zn13r
     &      *Zn34r+Zh13/CW*Zn14i*Zn33i-Zh13/CW*Zn14r*Zn33r
     &      -Zh23/CW*Zn13i*Zn44i+Zh23/CW*Zn13r*Zn44r-Zh23
     &      /CW*Zn14i*Zn43i+Zh23/CW*Zn14r*Zn43r-Zh13/SW*Zn23i
     &      *Zn34i+Zh13/SW*Zn23r*Zn34r-Zh13/SW*Zn24i*Zn33i
     &      +Zh13/SW*Zn24r*Zn33r+Zh23/SW*Zn23i*Zn44i-Zh23
     &      /SW*Zn23r*Zn44r+Zh23/SW*Zn24i*Zn43i-Zh23/SW*Zn24r
     &      *Zn43r-Zh33/CW*Zn13i*Zn34r*sb+Zh33/CW*Zn13i*Zn44r
     &      *cb-Zh33/CW*Zn13r*Zn34i*sb+Zh33/CW*Zn13r*Zn44i
     &      *cb-Zh33/CW*Zn14i*Zn33r*sb+Zh33/CW*Zn14i*Zn43r
     &      *cb-Zh33/CW*Zn14r*Zn33i*sb+Zh33/CW*Zn14r*Zn43i
     &      *cb+Zh33/SW*Zn23i*Zn34r*sb-Zh33/SW*Zn23i*Zn44r
     &      *cb+Zh33/SW*Zn23r*Zn34i*sb-Zh33/SW*Zn23r*Zn44i
     &      *cb+Zh33/SW*Zn24i*Zn33r*sb-Zh33/SW*Zn24i*Zn43r
     &      *cb+Zh33/SW*Zn24r*Zn33i*sb-Zh33/SW*Zn24r*Zn43i
     &      *cb)
      AAABR(997) = EE*(Zh13/CW*Zn13i*Zn34r+Zh13/CW*Zn13r
     &      *Zn34i+Zh13/CW*Zn14i*Zn33r+Zh13/CW*Zn14r*Zn33i
     &      -Zh23/CW*Zn13i*Zn44r-Zh23/CW*Zn13r*Zn44i-Zh23
     &      /CW*Zn14i*Zn43r-Zh23/CW*Zn14r*Zn43i-Zh13/SW*Zn23i
     &      *Zn34r-Zh13/SW*Zn23r*Zn34i-Zh13/SW*Zn24i*Zn33r
     &      -Zh13/SW*Zn24r*Zn33i+Zh23/SW*Zn23i*Zn44r+Zh23
     &      /SW*Zn23r*Zn44i+Zh23/SW*Zn24i*Zn43r+Zh23/SW*Zn24r
     &      *Zn43i+Zh33/CW*Zn13i*Zn34i*sb-Zh33/CW*Zn13i*Zn44i
     &      *cb-Zh33/CW*Zn13r*Zn34r*sb+Zh33/CW*Zn13r*Zn44r
     &      *cb+Zh33/CW*Zn14i*Zn33i*sb-Zh33/CW*Zn14i*Zn43i
     &      *cb-Zh33/CW*Zn14r*Zn33r*sb+Zh33/CW*Zn14r*Zn43r
     &      *cb-Zh33/SW*Zn23i*Zn34i*sb+Zh33/SW*Zn23i*Zn44i
     &      *cb+Zh33/SW*Zn23r*Zn34r*sb-Zh33/SW*Zn23r*Zn44r
     &      *cb-Zh33/SW*Zn24i*Zn33i*sb+Zh33/SW*Zn24i*Zn43i
     &      *cb+Zh33/SW*Zn24r*Zn33r*sb-Zh33/SW*Zn24r*Zn43r
     &      *cb)
      AAABR(998) = EE*(Zh11/CW*Zn13i*Zn34i-Zh11/CW*Zn13r
     &      *Zn34r+Zh11/CW*Zn14i*Zn33i-Zh11/CW*Zn14r*Zn33r
     &      -Zh21/CW*Zn13i*Zn44i+Zh21/CW*Zn13r*Zn44r-Zh21
     &      /CW*Zn14i*Zn43i+Zh21/CW*Zn14r*Zn43r-Zh11/SW*Zn23i
     &      *Zn34i+Zh11/SW*Zn23r*Zn34r-Zh11/SW*Zn24i*Zn33i
     &      +Zh11/SW*Zn24r*Zn33r+Zh21/SW*Zn23i*Zn44i-Zh21
     &      /SW*Zn23r*Zn44r+Zh21/SW*Zn24i*Zn43i-Zh21/SW*Zn24r
     &      *Zn43r-Zh31/CW*Zn13i*Zn34r*sb+Zh31/CW*Zn13i*Zn44r
     &      *cb-Zh31/CW*Zn13r*Zn34i*sb+Zh31/CW*Zn13r*Zn44i
     &      *cb-Zh31/CW*Zn14i*Zn33r*sb+Zh31/CW*Zn14i*Zn43r
     &      *cb-Zh31/CW*Zn14r*Zn33i*sb+Zh31/CW*Zn14r*Zn43i
     &      *cb+Zh31/SW*Zn23i*Zn34r*sb-Zh31/SW*Zn23i*Zn44r
     &      *cb+Zh31/SW*Zn23r*Zn34i*sb-Zh31/SW*Zn23r*Zn44i
     &      *cb+Zh31/SW*Zn24i*Zn33r*sb-Zh31/SW*Zn24i*Zn43r
     &      *cb+Zh31/SW*Zn24r*Zn33i*sb-Zh31/SW*Zn24r*Zn43i
     &      *cb)
      AAABR(999) = EE*(Zh11/CW*Zn13i*Zn34r+Zh11/CW*Zn13r
     &      *Zn34i+Zh11/CW*Zn14i*Zn33r+Zh11/CW*Zn14r*Zn33i
     &      -Zh21/CW*Zn13i*Zn44r-Zh21/CW*Zn13r*Zn44i-Zh21
     &      /CW*Zn14i*Zn43r-Zh21/CW*Zn14r*Zn43i-Zh11/SW*Zn23i
     &      *Zn34r-Zh11/SW*Zn23r*Zn34i-Zh11/SW*Zn24i*Zn33r
     &      -Zh11/SW*Zn24r*Zn33i+Zh21/SW*Zn23i*Zn44r+Zh21
     &      /SW*Zn23r*Zn44i+Zh21/SW*Zn24i*Zn43r+Zh21/SW*Zn24r
     &      *Zn43i+Zh31/CW*Zn13i*Zn34i*sb-Zh31/CW*Zn13i*Zn44i
     &      *cb-Zh31/CW*Zn13r*Zn34r*sb+Zh31/CW*Zn13r*Zn44r
     &      *cb+Zh31/CW*Zn14i*Zn33i*sb-Zh31/CW*Zn14i*Zn43i
     &      *cb-Zh31/CW*Zn14r*Zn33r*sb+Zh31/CW*Zn14r*Zn43r
     &      *cb-Zh31/SW*Zn23i*Zn34i*sb+Zh31/SW*Zn23i*Zn44i
     &      *cb+Zh31/SW*Zn23r*Zn34r*sb-Zh31/SW*Zn23r*Zn44r
     &      *cb-Zh31/SW*Zn24i*Zn33i*sb+Zh31/SW*Zn24i*Zn43i
     &      *cb+Zh31/SW*Zn24r*Zn33r*sb-Zh31/SW*Zn24r*Zn43r
     &      *cb)
      AAABR(1000) = EE*(Zh13/CW*Zn13i*Zn33i-Zh13/CW*Zn13r
     &      *Zn33r-Zh23/CW*Zn13i*Zn43i+Zh23/CW*Zn13r*Zn43r
     &      -Zh13/SW*Zn23i*Zn33i+Zh13/SW*Zn23r*Zn33r+Zh23
     &      /SW*Zn23i*Zn43i-Zh23/SW*Zn23r*Zn43r-Zh33/CW*Zn13i
     &      *Zn33r*sb+Zh33/CW*Zn13i*Zn43r*cb-Zh33/CW*Zn13r
     &      *Zn33i*sb+Zh33/CW*Zn13r*Zn43i*cb+Zh33/SW*Zn23i
     &      *Zn33r*sb-Zh33/SW*Zn23i*Zn43r*cb+Zh33/SW*Zn23r
     &      *Zn33i*sb-Zh33/SW*Zn23r*Zn43i*cb)
      end

      subroutine aaini02
      implicit none
#include "model.h"

      AAABR(1001) = EE*(Zh13/CW*Zn13i*Zn33r+Zh13/CW*Zn13r
     &      *Zn33i-Zh23/CW*Zn13i*Zn43r-Zh23/CW*Zn13r*Zn43i
     &      -Zh13/SW*Zn23i*Zn33r-Zh13/SW*Zn23r*Zn33i+Zh23
     &      /SW*Zn23i*Zn43r+Zh23/SW*Zn23r*Zn43i+Zh33/CW*Zn13i
     &      *Zn33i*sb-Zh33/CW*Zn13i*Zn43i*cb-Zh33/CW*Zn13r
     &      *Zn33r*sb+Zh33/CW*Zn13r*Zn43r*cb-Zh33/SW*Zn23i
     &      *Zn33i*sb+Zh33/SW*Zn23i*Zn43i*cb+Zh33/SW*Zn23r
     &      *Zn33r*sb-Zh33/SW*Zn23r*Zn43r*cb)
      AAABR(1002) = EE*(Zh12/CW*Zn13i*Zn34i-Zh12/CW*Zn13r
     &      *Zn34r+Zh12/CW*Zn14i*Zn33i-Zh12/CW*Zn14r*Zn33r
     &      -Zh22/CW*Zn13i*Zn44i+Zh22/CW*Zn13r*Zn44r-Zh22
     &      /CW*Zn14i*Zn43i+Zh22/CW*Zn14r*Zn43r-Zh12/SW*Zn23i
     &      *Zn34i+Zh12/SW*Zn23r*Zn34r-Zh12/SW*Zn24i*Zn33i
     &      +Zh12/SW*Zn24r*Zn33r+Zh22/SW*Zn23i*Zn44i-Zh22
     &      /SW*Zn23r*Zn44r+Zh22/SW*Zn24i*Zn43i-Zh22/SW*Zn24r
     &      *Zn43r-Zh32/CW*Zn13i*Zn34r*sb+Zh32/CW*Zn13i*Zn44r
     &      *cb-Zh32/CW*Zn13r*Zn34i*sb+Zh32/CW*Zn13r*Zn44i
     &      *cb-Zh32/CW*Zn14i*Zn33r*sb+Zh32/CW*Zn14i*Zn43r
     &      *cb-Zh32/CW*Zn14r*Zn33i*sb+Zh32/CW*Zn14r*Zn43i
     &      *cb+Zh32/SW*Zn23i*Zn34r*sb-Zh32/SW*Zn23i*Zn44r
     &      *cb+Zh32/SW*Zn23r*Zn34i*sb-Zh32/SW*Zn23r*Zn44i
     &      *cb+Zh32/SW*Zn24i*Zn33r*sb-Zh32/SW*Zn24i*Zn43r
     &      *cb+Zh32/SW*Zn24r*Zn33i*sb-Zh32/SW*Zn24r*Zn43i
     &      *cb)
      AAABR(1003) = EE*(Zh12/CW*Zn13i*Zn34r+Zh12/CW*Zn13r
     &      *Zn34i+Zh12/CW*Zn14i*Zn33r+Zh12/CW*Zn14r*Zn33i
     &      -Zh22/CW*Zn13i*Zn44r-Zh22/CW*Zn13r*Zn44i-Zh22
     &      /CW*Zn14i*Zn43r-Zh22/CW*Zn14r*Zn43i-Zh12/SW*Zn23i
     &      *Zn34r-Zh12/SW*Zn23r*Zn34i-Zh12/SW*Zn24i*Zn33r
     &      -Zh12/SW*Zn24r*Zn33i+Zh22/SW*Zn23i*Zn44r+Zh22
     &      /SW*Zn23r*Zn44i+Zh22/SW*Zn24i*Zn43r+Zh22/SW*Zn24r
     &      *Zn43i+Zh32/CW*Zn13i*Zn34i*sb-Zh32/CW*Zn13i*Zn44i
     &      *cb-Zh32/CW*Zn13r*Zn34r*sb+Zh32/CW*Zn13r*Zn44r
     &      *cb+Zh32/CW*Zn14i*Zn33i*sb-Zh32/CW*Zn14i*Zn43i
     &      *cb-Zh32/CW*Zn14r*Zn33r*sb+Zh32/CW*Zn14r*Zn43r
     &      *cb-Zh32/SW*Zn23i*Zn34i*sb+Zh32/SW*Zn23i*Zn44i
     &      *cb+Zh32/SW*Zn23r*Zn34r*sb-Zh32/SW*Zn23r*Zn44r
     &      *cb-Zh32/SW*Zn24i*Zn33i*sb+Zh32/SW*Zn24i*Zn43i
     &      *cb+Zh32/SW*Zn24r*Zn33r*sb-Zh32/SW*Zn24r*Zn43r
     &      *cb)
      AAABR(1004) = EE*(Zn14i/CW*Zn34r*cb+Zn14i/CW*Zn44r
     &      *sb+Zn14r/CW*Zn34i*cb+Zn14r/CW*Zn44i*sb-Zn24i
     &      /SW*Zn34r*cb-Zn24i/SW*Zn44r*sb-Zn24r/SW*Zn34i
     &      *cb-Zn24r/SW*Zn44i*sb)
      AAABR(1005) = EE*(Zn14i/CW*Zn34i*cb+Zn14i/CW*Zn44i
     &      *sb-Zn14r/CW*Zn34r*cb-Zn14r/CW*Zn44r*sb-Zn24i
     &      /SW*Zn34i*cb-Zn24i/SW*Zn44i*sb+Zn24r/SW*Zn34r
     &      *cb+Zn24r/SW*Zn44r*sb)
      AAABR(1006) = EE*(Zh12/CW*Zn14i*Zn34i-Zh12/CW*Zn14r
     &      *Zn34r-Zh22/CW*Zn14i*Zn44i+Zh22/CW*Zn14r*Zn44r
     &      -Zh12/SW*Zn24i*Zn34i+Zh12/SW*Zn24r*Zn34r+Zh22
     &      /SW*Zn24i*Zn44i-Zh22/SW*Zn24r*Zn44r-Zh32/CW*Zn14i
     &      *Zn34r*sb+Zh32/CW*Zn14i*Zn44r*cb-Zh32/CW*Zn14r
     &      *Zn34i*sb+Zh32/CW*Zn14r*Zn44i*cb+Zh32/SW*Zn24i
     &      *Zn34r*sb-Zh32/SW*Zn24i*Zn44r*cb+Zh32/SW*Zn24r
     &      *Zn34i*sb-Zh32/SW*Zn24r*Zn44i*cb)
      AAABR(1007) = EE*(Zh12/CW*Zn14i*Zn34r+Zh12/CW*Zn14r
     &      *Zn34i-Zh22/CW*Zn14i*Zn44r-Zh22/CW*Zn14r*Zn44i
     &      -Zh12/SW*Zn24i*Zn34r-Zh12/SW*Zn24r*Zn34i+Zh22
     &      /SW*Zn24i*Zn44r+Zh22/SW*Zn24r*Zn44i+Zh32/CW*Zn14i
     &      *Zn34i*sb-Zh32/CW*Zn14i*Zn44i*cb-Zh32/CW*Zn14r
     &      *Zn34r*sb+Zh32/CW*Zn14r*Zn44r*cb-Zh32/SW*Zn24i
     &      *Zn34i*sb+Zh32/SW*Zn24i*Zn44i*cb+Zh32/SW*Zn24r
     &      *Zn34r*sb-Zh32/SW*Zn24r*Zn44r*cb)
      AAABR(1008) = EE*(Zh11/CW*Zn14i*Zn34i-Zh11/CW*Zn14r
     &      *Zn34r-Zh21/CW*Zn14i*Zn44i+Zh21/CW*Zn14r*Zn44r
     &      -Zh11/SW*Zn24i*Zn34i+Zh11/SW*Zn24r*Zn34r+Zh21
     &      /SW*Zn24i*Zn44i-Zh21/SW*Zn24r*Zn44r-Zh31/CW*Zn14i
     &      *Zn34r*sb+Zh31/CW*Zn14i*Zn44r*cb-Zh31/CW*Zn14r
     &      *Zn34i*sb+Zh31/CW*Zn14r*Zn44i*cb+Zh31/SW*Zn24i
     &      *Zn34r*sb-Zh31/SW*Zn24i*Zn44r*cb+Zh31/SW*Zn24r
     &      *Zn34i*sb-Zh31/SW*Zn24r*Zn44i*cb)
      AAABR(1009) = EE*(Zh11/CW*Zn14i*Zn34r+Zh11/CW*Zn14r
     &      *Zn34i-Zh21/CW*Zn14i*Zn44r-Zh21/CW*Zn14r*Zn44i
     &      -Zh11/SW*Zn24i*Zn34r-Zh11/SW*Zn24r*Zn34i+Zh21
     &      /SW*Zn24i*Zn44r+Zh21/SW*Zn24r*Zn44i+Zh31/CW*Zn14i
     &      *Zn34i*sb-Zh31/CW*Zn14i*Zn44i*cb-Zh31/CW*Zn14r
     &      *Zn34r*sb+Zh31/CW*Zn14r*Zn44r*cb-Zh31/SW*Zn24i
     &      *Zn34i*sb+Zh31/SW*Zn24i*Zn44i*cb+Zh31/SW*Zn24r
     &      *Zn34r*sb-Zh31/SW*Zn24r*Zn44r*cb)
      AAABR(1010) = EE*(Zh13/CW*Zn14i*Zn34i-Zh13/CW*Zn14r
     &      *Zn34r-Zh23/CW*Zn14i*Zn44i+Zh23/CW*Zn14r*Zn44r
     &      -Zh13/SW*Zn24i*Zn34i+Zh13/SW*Zn24r*Zn34r+Zh23
     &      /SW*Zn24i*Zn44i-Zh23/SW*Zn24r*Zn44r-Zh33/CW*Zn14i
     &      *Zn34r*sb+Zh33/CW*Zn14i*Zn44r*cb-Zh33/CW*Zn14r
     &      *Zn34i*sb+Zh33/CW*Zn14r*Zn44i*cb+Zh33/SW*Zn24i
     &      *Zn34r*sb-Zh33/SW*Zn24i*Zn44r*cb+Zh33/SW*Zn24r
     &      *Zn34i*sb-Zh33/SW*Zn24r*Zn44i*cb)
      AAABR(1011) = EE*(Zh13/CW*Zn14i*Zn34r+Zh13/CW*Zn14r
     &      *Zn34i-Zh23/CW*Zn14i*Zn44r-Zh23/CW*Zn14r*Zn44i
     &      -Zh13/SW*Zn24i*Zn34r-Zh13/SW*Zn24r*Zn34i+Zh23
     &      /SW*Zn24i*Zn44r+Zh23/SW*Zn24r*Zn44i+Zh33/CW*Zn14i
     &      *Zn34i*sb-Zh33/CW*Zn14i*Zn44i*cb-Zh33/CW*Zn14r
     &      *Zn34r*sb+Zh33/CW*Zn14r*Zn44r*cb-Zh33/SW*Zn24i
     &      *Zn34i*sb+Zh33/SW*Zn24i*Zn44i*cb+Zh33/SW*Zn24r
     &      *Zn34r*sb-Zh33/SW*Zn24r*Zn44r*cb)
      AAABR(1012) = EE*(SW/CW*Zm21i*Zm22i+SW/CW*Zm21r*Zm22r
     &      -2D0*CW/SW*Zm11i*Zm12i-2D0*CW/SW*Zm11r*Zm12r
     &      -CW/SW*Zm21i*Zm22i-CW/SW*Zm21r*Zm22r)
      AAABR(1013) = EE*(SW/CW*Zm21i*Zm22r-SW/CW*Zm21r*Zm22i
     &      -2D0*CW/SW*Zm11i*Zm12r+2D0*CW/SW*Zm11r*Zm12i
     &      -CW/SW*Zm21i*Zm22r+CW/SW*Zm21r*Zm22i)
      AAABR(1014) = EE*(SW/CW*Zp21i*Zp22i+SW/CW*Zp21r*Zp22r
     &      -2D0*CW/SW*Zp11i*Zp12i-2D0*CW/SW*Zp11r*Zp12r
     &      -CW/SW*Zp21i*Zp22i-CW/SW*Zp21r*Zp22r)
      AAABR(1015) = EE*(SW/CW*Zp21i*Zp22r-SW/CW*Zp21r*Zp22i
     &      -2D0*CW/SW*Zp11i*Zp12r+2D0*CW/SW*Zp11r*Zp12i
     &      -CW/SW*Zp21i*Zp22r+CW/SW*Zp21r*Zp22i)
      AAABR(1016) = EE*(SW/CW*Zm21i**2+SW/CW*Zm21r**2-2D0
     &      *CW/SW*Zm11i**2-2D0*CW/SW*Zm11r**2-CW/SW*Zm21i
     &      **2-CW/SW*Zm21r**2)
      AAABR(1017) = EE*(SW/CW*Zp21i**2+SW/CW*Zp21r**2-2D0
     &      *CW/SW*Zp11i**2-2D0*CW/SW*Zp11r**2-CW/SW*Zp21i
     &      **2-CW/SW*Zp21r**2)
      AAABR(1018) = EE*(SW/CW*Zm22i**2+SW/CW*Zm22r**2-2D0
     &      *CW/SW*Zm12i**2-2D0*CW/SW*Zm12r**2-CW/SW*Zm22i
     &      **2-CW/SW*Zm22r**2)
      AAABR(1019) = EE*(SW/CW*Zp22i**2+SW/CW*Zp22r**2-2D0
     &      *CW/SW*Zp12i**2-2D0*CW/SW*Zp12r**2-CW/SW*Zp22i
     &      **2-CW/SW*Zp22r**2)
      AAABR(1020) = EE/SW*(2D0*Zm12i*Zn24i-2D0*Zm12r*Zn24r
     &      +Sqrt2*Zm22i*Zn34i-Sqrt2*Zm22r*Zn34r)
      AAABR(1021) = EE/SW*(2D0*Zm12i*Zn24r+2D0*Zm12r*Zn24i
     &      +Sqrt2*Zm22i*Zn34r+Sqrt2*Zm22r*Zn34i)
      AAABR(1022) = EE/SW*(2D0*Zn24i*Zp12i+2D0*Zn24r*Zp12r
     &      -Sqrt2*Zn44i*Zp22i-Sqrt2*Zn44r*Zp22r)
      AAABR(1023) = EE/SW*(2D0*Zn24i*Zp12r-2D0*Zn24r*Zp12i
     &      -Sqrt2*Zn44i*Zp22r+Sqrt2*Zn44r*Zp22i)
      AAABR(1024) = EE/SW*(2D0*Zm12i*Zn23i-2D0*Zm12r*Zn23r
     &      +Sqrt2*Zm22i*Zn33i-Sqrt2*Zm22r*Zn33r)
      AAABR(1025) = EE/SW*(2D0*Zm12i*Zn23r+2D0*Zm12r*Zn23i
     &      +Sqrt2*Zm22i*Zn33r+Sqrt2*Zm22r*Zn33i)
      AAABR(1026) = EE/SW*(2D0*Zn23i*Zp12i+2D0*Zn23r*Zp12r
     &      -Sqrt2*Zn43i*Zp22i-Sqrt2*Zn43r*Zp22r)
      AAABR(1027) = EE/SW*(2D0*Zn23i*Zp12r-2D0*Zn23r*Zp12i
     &      -Sqrt2*Zn43i*Zp22r+Sqrt2*Zn43r*Zp22i)
      AAABR(1028) = EE/SW*(2D0*Zm12i*Zn21i-2D0*Zm12r*Zn21r
     &      +Sqrt2*Zm22i*Zn31i-Sqrt2*Zm22r*Zn31r)
      AAABR(1029) = EE/SW*(2D0*Zm12i*Zn21r+2D0*Zm12r*Zn21i
     &      +Sqrt2*Zm22i*Zn31r+Sqrt2*Zm22r*Zn31i)
      AAABR(1030) = EE/SW*(2D0*Zn21i*Zp12i+2D0*Zn21r*Zp12r
     &      -Sqrt2*Zn41i*Zp22i-Sqrt2*Zn41r*Zp22r)
      AAABR(1031) = EE/SW*(2D0*Zn21i*Zp12r-2D0*Zn21r*Zp12i
     &      -Sqrt2*Zn41i*Zp22r+Sqrt2*Zn41r*Zp22i)
      AAABR(1032) = EE/SW*(2D0*Zm11i*Zn22i-2D0*Zm11r*Zn22r
     &      +Sqrt2*Zm21i*Zn32i-Sqrt2*Zm21r*Zn32r)
      AAABR(1033) = EE/SW*(2D0*Zm11i*Zn22r+2D0*Zm11r*Zn22i
     &      +Sqrt2*Zm21i*Zn32r+Sqrt2*Zm21r*Zn32i)
      AAABR(1034) = EE/SW*(2D0*Zn22i*Zp11i+2D0*Zn22r*Zp11r
     &      -Sqrt2*Zn42i*Zp21i-Sqrt2*Zn42r*Zp21r)
      AAABR(1035) = EE/SW*(2D0*Zn22i*Zp11r-2D0*Zn22r*Zp11i
     &      -Sqrt2*Zn42i*Zp21r+Sqrt2*Zn42r*Zp21i)
      AAABR(1036) = EE/SW*(2D0*Zm12i*Zn22i-2D0*Zm12r*Zn22r
     &      +Sqrt2*Zm22i*Zn32i-Sqrt2*Zm22r*Zn32r)
      AAABR(1037) = EE/SW*(2D0*Zm12i*Zn22r+2D0*Zm12r*Zn22i
     &      +Sqrt2*Zm22i*Zn32r+Sqrt2*Zm22r*Zn32i)
      AAABR(1038) = EE/SW*(2D0*Zn22i*Zp12i+2D0*Zn22r*Zp12r
     &      -Sqrt2*Zn42i*Zp22i-Sqrt2*Zn42r*Zp22r)
      AAABR(1039) = EE/SW*(2D0*Zn22i*Zp12r-2D0*Zn22r*Zp12i
     &      -Sqrt2*Zn42i*Zp22r+Sqrt2*Zn42r*Zp22i)
      AAABR(1040) = EE/SW*(2D0*Zm11i*Zn24i-2D0*Zm11r*Zn24r
     &      +Sqrt2*Zm21i*Zn34i-Sqrt2*Zm21r*Zn34r)
      AAABR(1041) = EE/SW*(2D0*Zm11i*Zn24r+2D0*Zm11r*Zn24i
     &      +Sqrt2*Zm21i*Zn34r+Sqrt2*Zm21r*Zn34i)
      AAABR(1042) = EE/SW*(2D0*Zn24i*Zp11i+2D0*Zn24r*Zp11r
     &      -Sqrt2*Zn44i*Zp21i-Sqrt2*Zn44r*Zp21r)
      AAABR(1043) = EE/SW*(2D0*Zn24i*Zp11r-2D0*Zn24r*Zp11i
     &      -Sqrt2*Zn44i*Zp21r+Sqrt2*Zn44r*Zp21i)
      AAABR(1044) = EE/SW*(2D0*Zm11i*Zn23i-2D0*Zm11r*Zn23r
     &      +Sqrt2*Zm21i*Zn33i-Sqrt2*Zm21r*Zn33r)
      AAABR(1045) = EE/SW*(2D0*Zm11i*Zn23r+2D0*Zm11r*Zn23i
     &      +Sqrt2*Zm21i*Zn33r+Sqrt2*Zm21r*Zn33i)
      AAABR(1046) = EE/SW*(2D0*Zn23i*Zp11i+2D0*Zn23r*Zp11r
     &      -Sqrt2*Zn43i*Zp21i-Sqrt2*Zn43r*Zp21r)
      AAABR(1047) = EE/SW*(2D0*Zn23i*Zp11r-2D0*Zn23r*Zp11i
     &      -Sqrt2*Zn43i*Zp21r+Sqrt2*Zn43r*Zp21i)
      AAABR(1048) = EE/SW*(2D0*Zm11i*Zn21i-2D0*Zm11r*Zn21r
     &      +Sqrt2*Zm21i*Zn31i-Sqrt2*Zm21r*Zn31r)
      AAABR(1049) = EE/SW*(2D0*Zm11i*Zn21r+2D0*Zm11r*Zn21i
     &      +Sqrt2*Zm21i*Zn31r+Sqrt2*Zm21r*Zn31i)
      AAABR(1050) = EE/SW*(2D0*Zn21i*Zp11i+2D0*Zn21r*Zp11r
     &      -Sqrt2*Zn41i*Zp21i-Sqrt2*Zn41r*Zp21r)
      AAABR(1051) = EE/SW*(2D0*Zn21i*Zp11r-2D0*Zn21r*Zp11i
     &      -Sqrt2*Zn41i*Zp21r+Sqrt2*Zn41r*Zp21i)
      AAABR(1052) = EE/CW*SW
      AAABR(1053) = EE/SW*Sqrt2
      AAABR(1054) = EE*(SW/CW+3D0*CW/SW)
      AAABR(1055) = EE*(SW/CW-3D0*CW/SW)
      AAABR(1056) = EE*(SW/CW*Zn31i*Zn32r-SW/CW*Zn31r*Zn32i
     &      -SW/CW*Zn41i*Zn42r+SW/CW*Zn41r*Zn42i+CW/SW*Zn31i
     &      *Zn32r-CW/SW*Zn31r*Zn32i-CW/SW*Zn41i*Zn42r+CW
     &      /SW*Zn41r*Zn42i)
      AAABR(1057) = EE*(SW/CW*Zn31i*Zn32i+SW/CW*Zn31r*Zn32r
     &      -SW/CW*Zn41i*Zn42i-SW/CW*Zn41r*Zn42r+CW/SW*Zn31i
     &      *Zn32i+CW/SW*Zn31r*Zn32r-CW/SW*Zn41i*Zn42i-CW
     &      /SW*Zn41r*Zn42r)
      AAABR(1058) = EE*(SW/CW*Zn31i*Zn33r-SW/CW*Zn31r*Zn33i
     &      -SW/CW*Zn41i*Zn43r+SW/CW*Zn41r*Zn43i+CW/SW*Zn31i
     &      *Zn33r-CW/SW*Zn31r*Zn33i-CW/SW*Zn41i*Zn43r+CW
     &      /SW*Zn41r*Zn43i)
      AAABR(1059) = EE*(SW/CW*Zn31i*Zn33i+SW/CW*Zn31r*Zn33r
     &      -SW/CW*Zn41i*Zn43i-SW/CW*Zn41r*Zn43r+CW/SW*Zn31i
     &      *Zn33i+CW/SW*Zn31r*Zn33r-CW/SW*Zn41i*Zn43i-CW
     &      /SW*Zn41r*Zn43r)
      AAABR(1060) = EE*(SW/CW*Zn31i*Zn34r-SW/CW*Zn31r*Zn34i
     &      -SW/CW*Zn41i*Zn44r+SW/CW*Zn41r*Zn44i+CW/SW*Zn31i
     &      *Zn34r-CW/SW*Zn31r*Zn34i-CW/SW*Zn41i*Zn44r+CW
     &      /SW*Zn41r*Zn44i)
      AAABR(1061) = EE*(SW/CW*Zn31i*Zn34i+SW/CW*Zn31r*Zn34r
     &      -SW/CW*Zn41i*Zn44i-SW/CW*Zn41r*Zn44r+CW/SW*Zn31i
     &      *Zn34i+CW/SW*Zn31r*Zn34r-CW/SW*Zn41i*Zn44i-CW
     &      /SW*Zn41r*Zn44r)
      AAABR(1062) = EE*(SW/CW*Zn31i**2+SW/CW*Zn31r**2-SW
     &      /CW*Zn41i**2-SW/CW*Zn41r**2+CW/SW*Zn31i**2+CW
     &      /SW*Zn31r**2-CW/SW*Zn41i**2-CW/SW*Zn41r**2)
      AAABR(1063) = EE*(SW/CW*Zn32i**2+SW/CW*Zn32r**2-SW
     &      /CW*Zn42i**2-SW/CW*Zn42r**2+CW/SW*Zn32i**2+CW
     &      /SW*Zn32r**2-CW/SW*Zn42i**2-CW/SW*Zn42r**2)
      AAABR(1064) = EE*(SW/CW*Zn32i*Zn34r-SW/CW*Zn32r*Zn34i
     &      -SW/CW*Zn42i*Zn44r+SW/CW*Zn42r*Zn44i+CW/SW*Zn32i
     &      *Zn34r-CW/SW*Zn32r*Zn34i-CW/SW*Zn42i*Zn44r+CW
     &      /SW*Zn42r*Zn44i)
      AAABR(1065) = EE*(SW/CW*Zn32i*Zn34i+SW/CW*Zn32r*Zn34r
     &      -SW/CW*Zn42i*Zn44i-SW/CW*Zn42r*Zn44r+CW/SW*Zn32i
     &      *Zn34i+CW/SW*Zn32r*Zn34r-CW/SW*Zn42i*Zn44i-CW
     &      /SW*Zn42r*Zn44r)
      AAABR(1066) = EE*(SW/CW*Zn32i*Zn33r-SW/CW*Zn32r*Zn33i
     &      -SW/CW*Zn42i*Zn43r+SW/CW*Zn42r*Zn43i+CW/SW*Zn32i
     &      *Zn33r-CW/SW*Zn32r*Zn33i-CW/SW*Zn42i*Zn43r+CW
     &      /SW*Zn42r*Zn43i)
      AAABR(1067) = EE*(SW/CW*Zn32i*Zn33i+SW/CW*Zn32r*Zn33r
     &      -SW/CW*Zn42i*Zn43i-SW/CW*Zn42r*Zn43r+CW/SW*Zn32i
     &      *Zn33i+CW/SW*Zn32r*Zn33r-CW/SW*Zn42i*Zn43i-CW
     &      /SW*Zn42r*Zn43r)
      AAABR(1068) = EE*(SW/CW*Zn33i*Zn34r-SW/CW*Zn33r*Zn34i
     &      -SW/CW*Zn43i*Zn44r+SW/CW*Zn43r*Zn44i+CW/SW*Zn33i
     &      *Zn34r-CW/SW*Zn33r*Zn34i-CW/SW*Zn43i*Zn44r+CW
     &      /SW*Zn43r*Zn44i)
      AAABR(1069) = EE*(SW/CW*Zn33i*Zn34i+SW/CW*Zn33r*Zn34r
     &      -SW/CW*Zn43i*Zn44i-SW/CW*Zn43r*Zn44r+CW/SW*Zn33i
     &      *Zn34i+CW/SW*Zn33r*Zn34r-CW/SW*Zn43i*Zn44i-CW
     &      /SW*Zn43r*Zn44r)
      AAABR(1070) = EE*(SW/CW*Zn33i**2+SW/CW*Zn33r**2-SW
     &      /CW*Zn43i**2-SW/CW*Zn43r**2+CW/SW*Zn33i**2+CW
     &      /SW*Zn33r**2-CW/SW*Zn43i**2-CW/SW*Zn43r**2)
      AAABR(1071) = EE*(SW/CW*Zn34i**2+SW/CW*Zn34r**2-SW
     &      /CW*Zn44i**2-SW/CW*Zn44r**2+CW/SW*Zn34i**2+CW
     &      /SW*Zn34r**2-CW/SW*Zn44i**2-CW/SW*Zn44r**2)
      AAABR(1072) = EE**2
      AAABR(1073) = EE**2*(SW/CW-CW/SW)
      AAABR(1074) = EE**2*(2D0-SW**2/CW**2-CW**2/SW**2)
      AAABR(1075) = EE**2/SW**2
      AAABR(1076) = EE**2/CW*(Zh12*sb-Zh22*cb)
      AAABR(1077) = EE**2/SW*(Zh12*sb-Zh22*cb)
      AAABR(1078) = EE**2/CW*Zh32
      AAABR(1079) = EE**2/SW*Zh32
      AAABR(1080) = EE**2/CW*(Zh11*sb-Zh21*cb)
      AAABR(1081) = EE**2/SW*(Zh11*sb-Zh21*cb)
      AAABR(1082) = EE**2/CW*Zh31
      AAABR(1083) = EE**2/SW*Zh31
      AAABR(1084) = EE**2/CW*(Zh13*sb-Zh23*cb)
      AAABR(1085) = EE**2/SW*(Zh13*sb-Zh23*cb)
      AAABR(1086) = EE**2/CW*Zh33
      AAABR(1087) = EE**2/SW*Zh33
      AAABR(1088) = EE**2*(SW/CW*Zl11**2+2D0*SW/CW*Zl41**2
     &      -CW/SW*Zl11**2)
      AAABR(1089) = EE**2*(SW/CW*Zl33i**2+SW/CW*Zl33r**2
     &      +2D0*SW/CW*Zl63i**2+2D0*SW/CW*Zl63r**2-CW/SW
     &      *Zl33i**2-CW/SW*Zl33r**2)
      AAABR(1090) = EE**2*(SW/CW*Zl22**2+2D0*SW/CW*Zl52**2
     &      -CW/SW*Zl22**2)
      AAABR(1091) = EE**2*(2D0*Zl22**2-SW**2/CW**2*Zl22**2
     &      -4D0*SW**2/CW**2*Zl52**2-CW**2/SW**2*Zl22**2)
      AAABR(1092) = EE**2*(2D0*Zl33i**2+2D0*Zl33r**2-SW**2
     &      /CW**2*Zl33i**2-SW**2/CW**2*Zl33r**2-4D0*SW**2
     &      /CW**2*Zl63i**2-4D0*SW**2/CW**2*Zl63r**2-CW**2
     &      /SW**2*Zl33i**2-CW**2/SW**2*Zl33r**2)
      AAABR(1093) = EE**2*(2D0*Zl11**2-SW**2/CW**2*Zl11**2
     &      -4D0*SW**2/CW**2*Zl41**2-CW**2/SW**2*Zl11**2)
      AAABR(1094) = EE**2/SW**2*(Zl33i**2+Zl33r**2)
      AAABR(1095) = EE**2/SW**2*Zl22**2
      AAABR(1096) = EE**2/SW**2*Zl11**2
      AAABR(1097) = EE**2*(SW/CW*Zl11*Zl14+2D0*SW/CW*Zl41
     &      *Zl44-CW/SW*Zl11*Zl14)
      AAABR(1098) = EE**2*(SW/CW*Zl33i*Zl36i+SW/CW*Zl33r
     &      *Zl36r+2D0*SW/CW*Zl63i*Zl66i+2D0*SW/CW*Zl63r
     &      *Zl66r-CW/SW*Zl33i*Zl36i-CW/SW*Zl33r*Zl36r)
      AAABR(1099) = EE**2*(SW/CW*Zl33i*Zl36r-SW/CW*Zl33r
     &      *Zl36i+2D0*SW/CW*Zl63i*Zl66r-2D0*SW/CW*Zl63r
     &      *Zl66i-CW/SW*Zl33i*Zl36r+CW/SW*Zl33r*Zl36i)
      AAABR(1100) = EE**2*(SW/CW*Zl22*Zl25+2D0*SW/CW*Zl52
     &      *Zl55-CW/SW*Zl22*Zl25)
      AAABR(1101) = EE**2*(2D0*Zl22*Zl25-SW**2/CW**2*Zl22
     &      *Zl25-4D0*SW**2/CW**2*Zl52*Zl55-CW**2/SW**2*Zl22
     &      *Zl25)
      AAABR(1102) = EE**2*(2D0*Zl33i*Zl36i+2D0*Zl33r*Zl36r
     &      -SW**2/CW**2*Zl33i*Zl36i-SW**2/CW**2*Zl33r*Zl36r
     &      -4D0*SW**2/CW**2*Zl63i*Zl66i-4D0*SW**2/CW**2
     &      *Zl63r*Zl66r-CW**2/SW**2*Zl33i*Zl36i-CW**2/SW
     &      **2*Zl33r*Zl36r)
      AAABR(1103) = EE**2*(2D0*Zl33i*Zl36r-2D0*Zl33r*Zl36i
     &      -SW**2/CW**2*Zl33i*Zl36r+SW**2/CW**2*Zl33r*Zl36i
     &      -4D0*SW**2/CW**2*Zl63i*Zl66r+4D0*SW**2/CW**2
     &      *Zl63r*Zl66i-CW**2/SW**2*Zl33i*Zl36r+CW**2/SW
     &      **2*Zl33r*Zl36i)
      AAABR(1104) = EE**2*(2D0*Zl11*Zl14-SW**2/CW**2*Zl11
     &      *Zl14-4D0*SW**2/CW**2*Zl41*Zl44-CW**2/SW**2*Zl11
     &      *Zl14)
      AAABR(1105) = EE**2/SW**2*Zl11*Zl14
      AAABR(1106) = EE**2/SW**2*Zl22*Zl25
      AAABR(1107) = EE**2/SW**2*(Zl33i*Zl36i+Zl33r*Zl36r)
      AAABR(1108) = EE**2/SW**2*(Zl33i*Zl36r-Zl33r*Zl36i)
      AAABR(1109) = EE**2/CW*Sqrt2*Zl22
      AAABR(1110) = EE**2/SW*Sqrt2*Zl22
      AAABR(1111) = EE**2/CW*Sqrt2*Zl33r
      AAABR(1112) = EE**2/CW*Sqrt2*Zl33i
      AAABR(1113) = EE**2/SW*Sqrt2*Zl33r
      AAABR(1114) = EE**2/SW*Sqrt2*Zl33i
      AAABR(1115) = EE**2/CW*Sqrt2*Zl11
      AAABR(1116) = EE**2/SW*Sqrt2*Zl11
      AAABR(1117) = EE**2*(SW/CW*Zl25**2+2D0*SW/CW*Zl55**2
     &      -CW/SW*Zl25**2)
      AAABR(1118) = EE**2*(SW/CW*Zl36i**2+SW/CW*Zl36r**2
     &      +2D0*SW/CW*Zl66i**2+2D0*SW/CW*Zl66r**2-CW/SW
     &      *Zl36i**2-CW/SW*Zl36r**2)
      AAABR(1119) = EE**2*(SW/CW*Zl14**2+2D0*SW/CW*Zl44**2
     &      -CW/SW*Zl14**2)
      AAABR(1120) = EE**2*(2D0*Zl36i**2+2D0*Zl36r**2-SW**2
     &      /CW**2*Zl36i**2-SW**2/CW**2*Zl36r**2-4D0*SW**2
     &      /CW**2*Zl66i**2-4D0*SW**2/CW**2*Zl66r**2-CW**2
     &      /SW**2*Zl36i**2-CW**2/SW**2*Zl36r**2)
      AAABR(1121) = EE**2*(2D0*Zl14**2-SW**2/CW**2*Zl14**2
     &      -4D0*SW**2/CW**2*Zl44**2-CW**2/SW**2*Zl14**2)
      AAABR(1122) = EE**2*(2D0*Zl25**2-SW**2/CW**2*Zl25**2
     &      -4D0*SW**2/CW**2*Zl55**2-CW**2/SW**2*Zl25**2)
      AAABR(1123) = EE**2/SW**2*(Zl36i**2+Zl36r**2)
      AAABR(1124) = EE**2/SW**2*Zl14**2
      AAABR(1125) = EE**2/SW**2*Zl25**2
      AAABR(1126) = EE**2/CW*Sqrt2*Zl14
      AAABR(1127) = EE**2/SW*Sqrt2*Zl14
      AAABR(1128) = EE**2/CW*Sqrt2*Zl25
      AAABR(1129) = EE**2/SW*Sqrt2*Zl25
      AAABR(1130) = EE**2/CW*Sqrt2*Zl36r
      AAABR(1131) = EE**2/CW*Sqrt2*Zl36i
      AAABR(1132) = EE**2/SW*Sqrt2*Zl36r
      AAABR(1133) = EE**2/SW*Sqrt2*Zl36i
      AAABR(1134) = EE**2*(2D0+SW**2/CW**2+CW**2/SW**2)
      AAABR(1135) = EE**2*(SW/CW*Zd11**2-2D0*SW/CW*Zd41**2
     &      +3D0*CW/SW*Zd11**2)
      AAABR(1136) = EE**2*(SW/CW*Zd33i**2+SW/CW*Zd33r**2
     &      -2D0*SW/CW*Zd63i**2-2D0*SW/CW*Zd63r**2+3D0*CW
     &      /SW*Zd33i**2+3D0*CW/SW*Zd33r**2)
      AAABR(1137) = EE**2*(SW/CW*Zd22**2-2D0*SW/CW*Zd52**2
     &      +3D0*CW/SW*Zd22**2)
      AAABR(1138) = EE*GG*(SW/CW*Zd33i**2+SW/CW*Zd33r**2
     &      -2D0*SW/CW*Zd63i**2-2D0*SW/CW*Zd63r**2+3D0*CW
     &      /SW*Zd33i**2+3D0*CW/SW*Zd33r**2)
      AAABR(1139) = EE*GG
      AAABR(1140) = EE*GG*(SW/CW*Zd11**2-2D0*SW/CW*Zd41**2
     &      +3D0*CW/SW*Zd11**2)
      AAABR(1141) = EE*GG*(SW/CW*Zd22**2-2D0*SW/CW*Zd52**2
     &      +3D0*CW/SW*Zd22**2)
      AAABR(1142) = EE**2*(6D0*Zd22**2+SW**2/CW**2*Zd22**2
     &      +4D0*SW**2/CW**2*Zd52**2+9D0*CW**2/SW**2*Zd22
     &      **2)
      AAABR(1143) = EE**2*(6D0*Zd11**2+SW**2/CW**2*Zd11**2
     &      +4D0*SW**2/CW**2*Zd41**2+9D0*CW**2/SW**2*Zd11
     &      **2)
      AAABR(1144) = EE**2*(6D0*Zd33i**2+6D0*Zd33r**2+SW**2
     &      /CW**2*Zd33i**2+SW**2/CW**2*Zd33r**2+4D0*SW**2
     &      /CW**2*Zd63i**2+4D0*SW**2/CW**2*Zd63r**2+9D0
     &      *CW**2/SW**2*Zd33i**2+9D0*CW**2/SW**2*Zd33r**2)
      AAABR(1145) = GG**2
      AAABR(1146) = EE**2/SW**2*(Zd33i**2+Zd33r**2)
      AAABR(1147) = EE**2/SW**2*Zd22**2
      AAABR(1148) = EE**2/SW**2*Zd11**2
      AAABR(1149) = EE**2*(SW/CW*Zd33i*Zd36i+SW/CW*Zd33r
     &      *Zd36r-2D0*SW/CW*Zd63i*Zd66i-2D0*SW/CW*Zd63r
     &      *Zd66r+3D0*CW/SW*Zd33i*Zd36i+3D0*CW/SW*Zd33r
     &      *Zd36r)
      AAABR(1150) = EE**2*(SW/CW*Zd33i*Zd36r-SW/CW*Zd33r
     &      *Zd36i-2D0*SW/CW*Zd63i*Zd66r+2D0*SW/CW*Zd63r
     &      *Zd66i+3D0*CW/SW*Zd33i*Zd36r-3D0*CW/SW*Zd33r
     &      *Zd36i)
      AAABR(1151) = EE**2*(SW/CW*Zd11*Zd14-2D0*SW/CW*Zd41
     &      *Zd44+3D0*CW/SW*Zd11*Zd14)
      AAABR(1152) = EE**2*(SW/CW*Zd22*Zd25-2D0*SW/CW*Zd52
     &      *Zd55+3D0*CW/SW*Zd22*Zd25)
      AAABR(1153) = EE*GG*(SW/CW*Zd33i*Zd36i+SW/CW*Zd33r
     &      *Zd36r-2D0*SW/CW*Zd63i*Zd66i-2D0*SW/CW*Zd63r
     &      *Zd66r+3D0*CW/SW*Zd33i*Zd36i+3D0*CW/SW*Zd33r
     &      *Zd36r)
      AAABR(1154) = EE*GG*(SW/CW*Zd33i*Zd36r-SW/CW*Zd33r
     &      *Zd36i-2D0*SW/CW*Zd63i*Zd66r+2D0*SW/CW*Zd63r
     &      *Zd66i+3D0*CW/SW*Zd33i*Zd36r-3D0*CW/SW*Zd33r
     &      *Zd36i)
      AAABR(1155) = EE*GG*(SW/CW*Zd22*Zd25-2D0*SW/CW*Zd52
     &      *Zd55+3D0*CW/SW*Zd22*Zd25)
      AAABR(1156) = EE*GG*(SW/CW*Zd11*Zd14-2D0*SW/CW*Zd41
     &      *Zd44+3D0*CW/SW*Zd11*Zd14)
      AAABR(1157) = EE**2*(6D0*Zd33i*Zd36i+6D0*Zd33r*Zd36r
     &      +SW**2/CW**2*Zd33i*Zd36i+SW**2/CW**2*Zd33r*Zd36r
     &      +4D0*SW**2/CW**2*Zd63i*Zd66i+4D0*SW**2/CW**2
     &      *Zd63r*Zd66r+9D0*CW**2/SW**2*Zd33i*Zd36i+9D0
     &      *CW**2/SW**2*Zd33r*Zd36r)
      AAABR(1158) = EE**2*(6D0*Zd33i*Zd36r-6D0*Zd33r*Zd36i
     &      +SW**2/CW**2*Zd33i*Zd36r-SW**2/CW**2*Zd33r*Zd36i
     &      +4D0*SW**2/CW**2*Zd63i*Zd66r-4D0*SW**2/CW**2
     &      *Zd63r*Zd66i+9D0*CW**2/SW**2*Zd33i*Zd36r-9D0
     &      *CW**2/SW**2*Zd33r*Zd36i)
      AAABR(1159) = EE**2*(6D0*Zd22*Zd25+SW**2/CW**2*Zd22
     &      *Zd25+4D0*SW**2/CW**2*Zd52*Zd55+9D0*CW**2/SW
     &      **2*Zd22*Zd25)
      AAABR(1160) = EE**2*(6D0*Zd11*Zd14+SW**2/CW**2*Zd11
     &      *Zd14+4D0*SW**2/CW**2*Zd41*Zd44+9D0*CW**2/SW
     &      **2*Zd11*Zd14)
      AAABR(1161) = EE**2/SW**2*(Zd33i*Zd36i+Zd33r*Zd36r)
      AAABR(1162) = EE**2/SW**2*(Zd33i*Zd36r-Zd33r*Zd36i)
      AAABR(1163) = EE**2/SW**2*Zd11*Zd14
      AAABR(1164) = EE**2/SW**2*Zd22*Zd25
      AAABR(1165) = EE**2/CW*Sqrt2*Zd11*Zu11
      AAABR(1166) = EE**2/SW*Sqrt2*Zd11*Zu11
      AAABR(1167) = EE**2/CW*Sqrt2*Zd22*Zu22
      AAABR(1168) = EE**2/SW*Sqrt2*Zd22*Zu22
      AAABR(1169) = EE**2/CW*Sqrt2*(Zd33i*Zu33i+Zd33r*Zu33r)
      AAABR(1170) = EE**2/CW*Sqrt2*(Zd33i*Zu33r-Zd33r*Zu33i)
      AAABR(1171) = EE**2/SW*Sqrt2*(Zd33i*Zu33i+Zd33r*Zu33r)
      AAABR(1172) = EE**2/SW*Sqrt2*(Zd33i*Zu33r-Zd33r*Zu33i)
      AAABR(1173) = EE*GG/SW*Sqrt2*Zd22*Zu22
      AAABR(1174) = EE*GG/SW*Sqrt2*(Zd33i*Zu33i+Zd33r*Zu33r)
      AAABR(1175) = EE*GG/SW*Sqrt2*(Zd33i*Zu33r-Zd33r*Zu33i)
      AAABR(1176) = EE*GG/SW*Sqrt2*Zd11*Zu11
      AAABR(1177) = EE**2/CW*Sqrt2*Zd11*Zu14
      AAABR(1178) = EE**2/SW*Sqrt2*Zd11*Zu14
      AAABR(1179) = EE**2/CW*Sqrt2*(Zd33i*Zu36i+Zd33r*Zu36r)
      AAABR(1180) = EE**2/CW*Sqrt2*(Zd33i*Zu36r-Zd33r*Zu36i)
      AAABR(1181) = EE**2/SW*Sqrt2*(Zd33i*Zu36i+Zd33r*Zu36r)
      AAABR(1182) = EE**2/SW*Sqrt2*(Zd33i*Zu36r-Zd33r*Zu36i)
      AAABR(1183) = EE**2/CW*Sqrt2*Zd22*Zu25
      AAABR(1184) = EE**2/SW*Sqrt2*Zd22*Zu25
      AAABR(1185) = EE*GG/SW*Sqrt2*Zd11*Zu14
      AAABR(1186) = EE*GG/SW*Sqrt2*(Zd33i*Zu36i+Zd33r*Zu36r)
      AAABR(1187) = EE*GG/SW*Sqrt2*(Zd33i*Zu36r-Zd33r*Zu36i)
      AAABR(1188) = EE*GG/SW*Sqrt2*Zd22*Zu25
      AAABR(1189) = EE**2*(SW/CW*Zd14**2-2D0*SW/CW*Zd44**2
     &      +3D0*CW/SW*Zd14**2)
      AAABR(1190) = EE**2*(SW/CW*Zd36i**2+SW/CW*Zd36r**2
     &      -2D0*SW/CW*Zd66i**2-2D0*SW/CW*Zd66r**2+3D0*CW
     &      /SW*Zd36i**2+3D0*CW/SW*Zd36r**2)
      AAABR(1191) = EE**2*(SW/CW*Zd25**2-2D0*SW/CW*Zd55**2
     &      +3D0*CW/SW*Zd25**2)
      AAABR(1192) = EE*GG*(SW/CW*Zd14**2-2D0*SW/CW*Zd44**2
     &      +3D0*CW/SW*Zd14**2)
      AAABR(1193) = EE*GG*(SW/CW*Zd25**2-2D0*SW/CW*Zd55**2
     &      +3D0*CW/SW*Zd25**2)
      AAABR(1194) = EE*GG*(SW/CW*Zd36i**2+SW/CW*Zd36r**2
     &      -2D0*SW/CW*Zd66i**2-2D0*SW/CW*Zd66r**2+3D0*CW
     &      /SW*Zd36i**2+3D0*CW/SW*Zd36r**2)
      AAABR(1195) = EE**2*(6D0*Zd14**2+SW**2/CW**2*Zd14**2
     &      +4D0*SW**2/CW**2*Zd44**2+9D0*CW**2/SW**2*Zd14
     &      **2)
      AAABR(1196) = EE**2*(6D0*Zd36i**2+6D0*Zd36r**2+SW**2
     &      /CW**2*Zd36i**2+SW**2/CW**2*Zd36r**2+4D0*SW**2
     &      /CW**2*Zd66i**2+4D0*SW**2/CW**2*Zd66r**2+9D0
     &      *CW**2/SW**2*Zd36i**2+9D0*CW**2/SW**2*Zd36r**2)
      AAABR(1197) = EE**2*(6D0*Zd25**2+SW**2/CW**2*Zd25**2
     &      +4D0*SW**2/CW**2*Zd55**2+9D0*CW**2/SW**2*Zd25
     &      **2)
      AAABR(1198) = EE**2/SW**2*Zd25**2
      AAABR(1199) = EE**2/SW**2*Zd14**2
      AAABR(1200) = EE**2/SW**2*(Zd36i**2+Zd36r**2)
      AAABR(1201) = EE**2/CW*Sqrt2*(Zd36i*Zu33i+Zd36r*Zu33r)
      AAABR(1202) = EE**2/CW*Sqrt2*(Zd36i*Zu33r-Zd36r*Zu33i)
      AAABR(1203) = EE**2/SW*Sqrt2*(Zd36i*Zu33i+Zd36r*Zu33r)
      AAABR(1204) = EE**2/SW*Sqrt2*(Zd36i*Zu33r-Zd36r*Zu33i)
      AAABR(1205) = EE**2/CW*Sqrt2*Zd14*Zu11
      AAABR(1206) = EE**2/SW*Sqrt2*Zd14*Zu11
      AAABR(1207) = EE**2/CW*Sqrt2*Zd25*Zu22
      AAABR(1208) = EE**2/SW*Sqrt2*Zd25*Zu22
      AAABR(1209) = EE*GG/SW*Sqrt2*(Zd36i*Zu33i+Zd36r*Zu33r)
      AAABR(1210) = EE*GG/SW*Sqrt2*(Zd36i*Zu33r-Zd36r*Zu33i)
      AAABR(1211) = EE*GG/SW*Sqrt2*Zd25*Zu22
      AAABR(1212) = EE*GG/SW*Sqrt2*Zd14*Zu11
      AAABR(1213) = EE**2/CW*Sqrt2*(Zd36i*Zu36i+Zd36r*Zu36r)
      AAABR(1214) = EE**2/CW*Sqrt2*(Zd36i*Zu36r-Zd36r*Zu36i)
      AAABR(1215) = EE**2/SW*Sqrt2*(Zd36i*Zu36i+Zd36r*Zu36r)
      AAABR(1216) = EE**2/SW*Sqrt2*(Zd36i*Zu36r-Zd36r*Zu36i)
      AAABR(1217) = EE**2/CW*Sqrt2*Zd25*Zu25
      AAABR(1218) = EE**2/SW*Sqrt2*Zd25*Zu25
      AAABR(1219) = EE**2/CW*Sqrt2*Zd14*Zu14
      AAABR(1220) = EE**2/SW*Sqrt2*Zd14*Zu14
      AAABR(1221) = EE*GG/SW*Sqrt2*Zd14*Zu14
      AAABR(1222) = EE*GG/SW*Sqrt2*Zd25*Zu25
      AAABR(1223) = EE*GG/SW*Sqrt2*(Zd36i*Zu36i+Zd36r*Zu36r)
      AAABR(1224) = EE*GG/SW*Sqrt2*(Zd36i*Zu36r-Zd36r*Zu36i)
      AAABR(1225) = EE**2*(SW/CW*Zu11**2+4D0*SW/CW*Zu41**2
     &      -3D0*CW/SW*Zu11**2)
      AAABR(1226) = EE**2*(SW/CW*Zu33i**2+SW/CW*Zu33r**2
     &      +4D0*SW/CW*Zu63i**2+4D0*SW/CW*Zu63r**2-3D0*CW
     &      /SW*Zu33i**2-3D0*CW/SW*Zu33r**2)
      AAABR(1227) = EE**2*(SW/CW*Zu22**2+4D0*SW/CW*Zu52**2
     &      -3D0*CW/SW*Zu22**2)
      AAABR(1228) = EE*GG*(SW/CW*Zu11**2+4D0*SW/CW*Zu41**2
     &      -3D0*CW/SW*Zu11**2)
      AAABR(1229) = EE*GG*(SW/CW*Zu33i**2+SW/CW*Zu33r**2
     &      +4D0*SW/CW*Zu63i**2+4D0*SW/CW*Zu63r**2-3D0*CW
     &      /SW*Zu33i**2-3D0*CW/SW*Zu33r**2)
      AAABR(1230) = EE*GG*(SW/CW*Zu22**2+4D0*SW/CW*Zu52**2
     &      -3D0*CW/SW*Zu22**2)
      AAABR(1231) = EE**2*(6D0*Zu22**2-SW**2/CW**2*Zu22**2
     &      -16D0*SW**2/CW**2*Zu52**2-9D0*CW**2/SW**2*Zu22
     &      **2)
      AAABR(1232) = EE**2*(6D0*Zu11**2-SW**2/CW**2*Zu11**2
     &      -16D0*SW**2/CW**2*Zu41**2-9D0*CW**2/SW**2*Zu11
     &      **2)
      AAABR(1233) = EE**2*(6D0*Zu33i**2+6D0*Zu33r**2-SW**2
     &      /CW**2*Zu33i**2-SW**2/CW**2*Zu33r**2-16D0*SW
     &      **2/CW**2*Zu63i**2-16D0*SW**2/CW**2*Zu63r**2
     &      -9D0*CW**2/SW**2*Zu33i**2-9D0*CW**2/SW**2*Zu33r
     &      **2)
      AAABR(1234) = EE**2/SW**2*(Zu33i**2+Zu33r**2)
      AAABR(1235) = EE**2/SW**2*Zu11**2
      AAABR(1236) = EE**2/SW**2*Zu22**2
      AAABR(1237) = EE**2*(SW/CW*Zu33i*Zu36i+SW/CW*Zu33r
     &      *Zu36r+4D0*SW/CW*Zu63i*Zu66i+4D0*SW/CW*Zu63r
     &      *Zu66r-3D0*CW/SW*Zu33i*Zu36i-3D0*CW/SW*Zu33r
     &      *Zu36r)
      AAABR(1238) = EE**2*(SW/CW*Zu33i*Zu36r-SW/CW*Zu33r
     &      *Zu36i+4D0*SW/CW*Zu63i*Zu66r-4D0*SW/CW*Zu63r
     &      *Zu66i-3D0*CW/SW*Zu33i*Zu36r+3D0*CW/SW*Zu33r
     &      *Zu36i)
      AAABR(1239) = EE**2*(SW/CW*Zu11*Zu14+4D0*SW/CW*Zu41
     &      *Zu44-3D0*CW/SW*Zu11*Zu14)
      AAABR(1240) = EE**2*(SW/CW*Zu22*Zu25+4D0*SW/CW*Zu52
     &      *Zu55-3D0*CW/SW*Zu22*Zu25)
      AAABR(1241) = EE*GG*(SW/CW*Zu33i*Zu36i+SW/CW*Zu33r
     &      *Zu36r+4D0*SW/CW*Zu63i*Zu66i+4D0*SW/CW*Zu63r
     &      *Zu66r-3D0*CW/SW*Zu33i*Zu36i-3D0*CW/SW*Zu33r
     &      *Zu36r)
      AAABR(1242) = EE*GG*(SW/CW*Zu33i*Zu36r-SW/CW*Zu33r
     &      *Zu36i+4D0*SW/CW*Zu63i*Zu66r-4D0*SW/CW*Zu63r
     &      *Zu66i-3D0*CW/SW*Zu33i*Zu36r+3D0*CW/SW*Zu33r
     &      *Zu36i)
      AAABR(1243) = EE*GG*(SW/CW*Zu22*Zu25+4D0*SW/CW*Zu52
     &      *Zu55-3D0*CW/SW*Zu22*Zu25)
      AAABR(1244) = EE*GG*(SW/CW*Zu11*Zu14+4D0*SW/CW*Zu41
     &      *Zu44-3D0*CW/SW*Zu11*Zu14)
      AAABR(1245) = EE**2*(6D0*Zu22*Zu25-SW**2/CW**2*Zu22
     &      *Zu25-16D0*SW**2/CW**2*Zu52*Zu55-9D0*CW**2/SW
     &      **2*Zu22*Zu25)
      AAABR(1246) = EE**2*(6D0*Zu11*Zu14-SW**2/CW**2*Zu11
     &      *Zu14-16D0*SW**2/CW**2*Zu41*Zu44-9D0*CW**2/SW
     &      **2*Zu11*Zu14)
      AAABR(1247) = EE**2*(6D0*Zu33i*Zu36i+6D0*Zu33r*Zu36r
     &      -SW**2/CW**2*Zu33i*Zu36i-SW**2/CW**2*Zu33r*Zu36r
     &      -16D0*SW**2/CW**2*Zu63i*Zu66i-16D0*SW**2/CW**2
     &      *Zu63r*Zu66r-9D0*CW**2/SW**2*Zu33i*Zu36i-9D0
     &      *CW**2/SW**2*Zu33r*Zu36r)
      AAABR(1248) = EE**2*(6D0*Zu33i*Zu36r-6D0*Zu33r*Zu36i
     &      -SW**2/CW**2*Zu33i*Zu36r+SW**2/CW**2*Zu33r*Zu36i
     &      -16D0*SW**2/CW**2*Zu63i*Zu66r+16D0*SW**2/CW**2
     &      *Zu63r*Zu66i-9D0*CW**2/SW**2*Zu33i*Zu36r+9D0
     &      *CW**2/SW**2*Zu33r*Zu36i)
      AAABR(1249) = EE**2/SW**2*Zu11*Zu14
      AAABR(1250) = EE**2/SW**2*Zu22*Zu25
      AAABR(1251) = EE**2/SW**2*(Zu33i*Zu36i+Zu33r*Zu36r)
      AAABR(1252) = EE**2/SW**2*(Zu33i*Zu36r-Zu33r*Zu36i)
      AAABR(1253) = EE**2*(SW/CW*Zu25**2+4D0*SW/CW*Zu55**2
     &      -3D0*CW/SW*Zu25**2)
      AAABR(1254) = EE**2*(SW/CW*Zu14**2+4D0*SW/CW*Zu44**2
     &      -3D0*CW/SW*Zu14**2)
      AAABR(1255) = EE**2*(SW/CW*Zu36i**2+SW/CW*Zu36r**2
     &      +4D0*SW/CW*Zu66i**2+4D0*SW/CW*Zu66r**2-3D0*CW
     &      /SW*Zu36i**2-3D0*CW/SW*Zu36r**2)
      AAABR(1256) = EE*GG*(SW/CW*Zu25**2+4D0*SW/CW*Zu55**2
     &      -3D0*CW/SW*Zu25**2)
      AAABR(1257) = EE*GG*(SW/CW*Zu36i**2+SW/CW*Zu36r**2
     &      +4D0*SW/CW*Zu66i**2+4D0*SW/CW*Zu66r**2-3D0*CW
     &      /SW*Zu36i**2-3D0*CW/SW*Zu36r**2)
      AAABR(1258) = EE*GG*(SW/CW*Zu14**2+4D0*SW/CW*Zu44**2
     &      -3D0*CW/SW*Zu14**2)
      AAABR(1259) = EE**2*(6D0*Zu36i**2+6D0*Zu36r**2-SW**2
     &      /CW**2*Zu36i**2-SW**2/CW**2*Zu36r**2-16D0*SW
     &      **2/CW**2*Zu66i**2-16D0*SW**2/CW**2*Zu66r**2
     &      -9D0*CW**2/SW**2*Zu36i**2-9D0*CW**2/SW**2*Zu36r
     &      **2)
      AAABR(1260) = EE**2*(6D0*Zu14**2-SW**2/CW**2*Zu14**2
     &      -16D0*SW**2/CW**2*Zu44**2-9D0*CW**2/SW**2*Zu14
     &      **2)
      AAABR(1261) = EE**2*(6D0*Zu25**2-SW**2/CW**2*Zu25**2
     &      -16D0*SW**2/CW**2*Zu55**2-9D0*CW**2/SW**2*Zu25
     &      **2)
      AAABR(1262) = EE**2/SW**2*Zu14**2
      AAABR(1263) = EE**2/SW**2*(Zu36i**2+Zu36r**2)
      AAABR(1264) = EE**2/SW**2*Zu25**2
      AAABR(1265) = EE**2/CW
      AAABR(1266) = EE**2/SW
      AAABR(1267) = EE**2/CW*(Zh13*cb+Zh23*sb)
      AAABR(1268) = EE**2/SW*(Zh13*cb+Zh23*sb)
      AAABR(1269) = EE**2/CW*(Zh11*cb+Zh21*sb)
      AAABR(1270) = EE**2/SW*(Zh11*cb+Zh21*sb)
      AAABR(1271) = EE**2/CW*(Zh12*cb+Zh22*sb)
      AAABR(1272) = EE**2/SW*(Zh12*cb+Zh22*sb)
      AAABR(1273) = CW*EE**2/SW
      AAABR(1274) = CW**2*EE**2/SW**2
      end

      subroutine mtrini
      implicit none
#include "model.h"

      integer m1,m2,m3,m4

      do m1 = 1,3
      MTR033(m1)=0D0
      MTR034(m1)=0D0
      MTR037(m1)=0D0
      MTR038(m1)=0D0
      MTR039(m1)=0D0
      MTR040(m1)=0D0
      MTR043(m1)=0D0
      MTR044(m1)=0D0
      MTR045(m1)=0D0
      MTR046(m1)=0D0
      MTR049(m1)=0D0
      MTR050(m1)=0D0
      MTR051(m1)=0D0
      MTR052(m1)=0D0
      MTR055(m1)=0D0
      MTR058(m1)=0D0
      MTR059(m1)=0D0
      MTR060(m1)=0D0
      MTR061(m1)=0D0
      MTR062(m1)=0D0
      MTR063(m1)=0D0
      MTR064(m1)=0D0
      MTR065(m1)=0D0
      MTR066(m1)=0D0
      MTR067(m1)=0D0
      MTR068(m1)=0D0
      MTR069(m1)=0D0
      MTR070(m1)=0D0
      MTR071(m1)=0D0
      MTR072(m1)=0D0
      MTR073(m1)=0D0
      MTR074(m1)=0D0
      MTR075(m1)=0D0
      MTR076(m1)=0D0
      MTR077(m1)=0D0
      MTR078(m1)=0D0
      MTR079(m1)=0D0
      MTR080(m1)=0D0
      MTR081(m1)=0D0
      MTR082(m1)=0D0
      MTR083(m1)=0D0
      MTR084(m1)=0D0
      MTR085(m1)=0D0
      MTR087(m1)=0D0
      MTR102(m1)=0D0
      MTR105(m1)=0D0
      MTR106(m1)=0D0
      MTR120(m1)=0D0
      MTR123(m1)=0D0
      MTR124(m1)=0D0
      MTR125(m1)=0D0
      MTR126(m1)=0D0
      MTR131(m1)=0D0
      MTR132(m1)=0D0
      MTR133(m1)=0D0
      MTR134(m1)=0D0
      MTR143(m1)=0D0
      MTR146(m1)=0D0
      MTR147(m1)=0D0
      MTR148(m1)=0D0
      MTR149(m1)=0D0
      MTR158(m1)=0D0
      MTR159(m1)=0D0
      MTR160(m1)=0D0
      MTR161(m1)=0D0
      MTR166(m1)=0D0
      MTR167(m1)=0D0
      MTR168(m1)=0D0
      MTR169(m1)=0D0
      MTR193(m1)=0D0
      MTR194(m1)=0D0
      MTR195(m1)=0D0
      MTR196(m1)=0D0
      MTR197(m1)=0D0
      MTR198(m1)=0D0
      MTR199(m1)=0D0
      MTR200(m1)=0D0
      MTR201(m1)=0D0
      MTR202(m1)=0D0
      MTR203(m1)=0D0
      MTR204(m1)=0D0
      MTR205(m1)=0D0
      MTR206(m1)=0D0
      MTR207(m1)=0D0
      MTR208(m1)=0D0
      MTR209(m1)=0D0
      MTR210(m1)=0D0
      MTR211(m1)=0D0
      MTR212(m1)=0D0
      MTR213(m1)=0D0
      MTR214(m1)=0D0
      MTR215(m1)=0D0
      MTR216(m1)=0D0
      MTR217(m1)=0D0
      MTR218(m1)=0D0
      MTR219(m1)=0D0
      MTR220(m1)=0D0
      MTR221(m1)=0D0
      MTR222(m1)=0D0
      MTR223(m1)=0D0
      MTR224(m1)=0D0
      MTR225(m1)=0D0
      MTR226(m1)=0D0
      MTR227(m1)=0D0
      MTR228(m1)=0D0
      MTR229(m1)=0D0
      MTR230(m1)=0D0
      MTR231(m1)=0D0
      MTR232(m1)=0D0
      MTR233(m1)=0D0
      MTR234(m1)=0D0
      MTR235(m1)=0D0
      MTR236(m1)=0D0
      MTR237(m1)=0D0
      MTR238(m1)=0D0
      MTR239(m1)=0D0
      MTR240(m1)=0D0
      MTR241(m1)=0D0
      MTR242(m1)=0D0
      MTR243(m1)=0D0
      MTR244(m1)=0D0
      MTR245(m1)=0D0
      MTR246(m1)=0D0
      MTR247(m1)=0D0
      MTR248(m1)=0D0
      MTR249(m1)=0D0
      MTR250(m1)=0D0
      MTR251(m1)=0D0
      MTR252(m1)=0D0
      MTR253(m1)=0D0
      MTR254(m1)=0D0
      MTR255(m1)=0D0
      MTR256(m1)=0D0
      MTR257(m1)=0D0
      MTR258(m1)=0D0
      MTR259(m1)=0D0
      MTR260(m1)=0D0
      MTR261(m1)=0D0
      MTR262(m1)=0D0
      MTR263(m1)=0D0
      MTR264(m1)=0D0
      MTR265(m1)=0D0
      MTR266(m1)=0D0
      MTR267(m1)=0D0
      MTR268(m1)=0D0
      MTR269(m1)=0D0
      MTR270(m1)=0D0
      MTR271(m1)=0D0
      MTR272(m1)=0D0
      MTR273(m1)=0D0
      MTR274(m1)=0D0
      enddo
      do m1 = 1,3
      do m2 = 1,3
      MTR086(m1,m2)=0D0
      MTR103(m1,m2)=0D0
      MTR104(m1,m2)=0D0
      MTR121(m1,m2)=0D0
      MTR122(m1,m2)=0D0
      MTR144(m1,m2)=0D0
      MTR145(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,2
      do m2 = 1,3
      MTR099(m1,m2)=0D0
      MTR100(m1,m2)=0D0
      MTR101(m1,m2)=0D0
      MTR111(m1,m2)=0D0
      MTR112(m1,m2)=0D0
      MTR113(m1,m2)=0D0
      MTR114(m1,m2)=0D0
      MTR127(m1,m2)=0D0
      MTR128(m1,m2)=0D0
      MTR129(m1,m2)=0D0
      MTR130(m1,m2)=0D0
      MTR135(m1,m2)=0D0
      MTR136(m1,m2)=0D0
      MTR137(m1,m2)=0D0
      MTR138(m1,m2)=0D0
      MTR150(m1,m2)=0D0
      MTR151(m1,m2)=0D0
      MTR152(m1,m2)=0D0
      MTR153(m1,m2)=0D0
      MTR170(m1,m2)=0D0
      MTR171(m1,m2)=0D0
      enddo
      enddo
      quuMass(1)=Mu
      quuMass(2)=Mc
      quuMass(3)=Mt
      qudMass(1)=Md
      qudMass(2)=Ms
      qudMass(3)=Mb
      lpdMass(1)=Me
      lpdMass(2)=Mm
      lpdMass(3)=Ml
      neuMass(1)=MNE1
      neuMass(2)=MNE2
      neuMass(3)=MNE3
      neuMass(4)=MNE4
      chaMass(1)=MC1
      chaMass(2)=MC2
      sluMass(1)=MSne
      sluMass(2)=MSnm
      sluMass(3)=MSnl
      sldMass(1)=MSe1
      sldMass(2)=MSm1
      sldMass(3)=MSl1
      sleMass(1)=MSe2
      sleMass(2)=MSm2
      sleMass(3)=MSl2
      squMass(1)=MSu1
      squMass(2)=MSc1
      squMass(3)=MSt1
      sqvMass(1)=MSu2
      sqvMass(2)=MSc2
      sqvMass(3)=MSt2
      sqdMass(1)=MSd1
      sqdMass(2)=MSs1
      sqdMass(3)=MSb1
      sqeMass(1)=MSd2
      sqeMass(2)=MSs2
      sqeMass(3)=MSb2
      hisMass(1)=Mh1
      hisMass(2)=Mh2
      hisMass(3)=Mh3
      MTR001(1)=AAABR(9)
      MTR001(2)=AAABR(7)
      MTR001(3)=AAABR(8)
      MTR002(1)=AAABR(12)
      MTR002(2)=AAABR(14)
      MTR002(3)=AAABR(13)
      MTR003(1)=AAABR(18)
      MTR003(2)=AAABR(16)
      MTR003(3)=AAABR(17)
      MTR004(1)=AAABR(19)+cI*AAABR(20)
      MTR004(2)=AAABR(21)+cI*AAABR(22)
      MTR004(3)=-AAABR(23)+cI*AAABR(24)
      MTR005(1)=AAABR(29)+cI*AAABR(30)
      MTR005(2)=AAABR(25)+cI*AAABR(26)
      MTR005(3)=-AAABR(27)+cI*AAABR(28)
      MTR006(1)=2D0*AAABR(31)
      MTR006(2)=AAABR(34)+cI*AAABR(35)
      MTR006(3)=-AAABR(32)-cI*AAABR(33)
      MTR007(1)=2D0*AAABR(40)
      MTR007(2)=AAABR(38)+cI*AAABR(39)
      MTR007(3)=-AAABR(36)-cI*AAABR(37)
      MTR008(1)=2D0*AAABR(45)
      MTR008(2)=AAABR(43)+cI*AAABR(44)
      MTR008(3)=-AAABR(41)-cI*AAABR(42)
      MTR009(1)=2D0*AAABR(48)
      MTR009(2)=AAABR(49)+cI*AAABR(50)
      MTR009(3)=-AAABR(46)-cI*AAABR(47)
      MTR010(1)=AAABR(51)-cI*AAABR(52)
      MTR010(2)=AAABR(53)-cI*AAABR(54)
      MTR010(3)=AAABR(55)-cI*AAABR(56)
      MTR011(1)=AAABR(19)-cI*AAABR(20)
      MTR011(2)=AAABR(21)-cI*AAABR(22)
      MTR011(3)=-AAABR(23)-cI*AAABR(24)
      MTR012(1)=AAABR(29)-cI*AAABR(30)
      MTR012(2)=AAABR(25)-cI*AAABR(26)
      MTR012(3)=-AAABR(27)-cI*AAABR(28)
      MTR013(1)=2D0*AAABR(31)
      MTR013(2)=AAABR(34)-cI*AAABR(35)
      MTR013(3)=-AAABR(32)+cI*AAABR(33)
      MTR014(1)=2D0*AAABR(45)
      MTR014(2)=AAABR(43)-cI*AAABR(44)
      MTR014(3)=-AAABR(41)+cI*AAABR(42)
      MTR015(1)=2D0*AAABR(40)
      MTR015(2)=AAABR(38)-cI*AAABR(39)
      MTR015(3)=-AAABR(36)+cI*AAABR(37)
      MTR016(1)=2D0*AAABR(48)
      MTR016(2)=AAABR(49)-cI*AAABR(50)
      MTR016(3)=-AAABR(46)+cI*AAABR(47)
      MTR017(1)=AAABR(51)+cI*AAABR(52)
      MTR017(2)=AAABR(53)+cI*AAABR(54)
      MTR017(3)=AAABR(55)+cI*AAABR(56)
      MTR018(1)=AAABR(59)+cI*AAABR(60)
      MTR018(2)=AAABR(57)+cI*AAABR(58)
      MTR018(3)=-AAABR(61)+cI*AAABR(62)
      MTR019(1)=AAABR(63)
      MTR019(2)=AAABR(65)
      MTR019(3)=AAABR(64)
      MTR020(1)=cI*AAABR(68)+AAABR(69)
      MTR020(2)=cI*AAABR(66)+AAABR(67)
      MTR020(3)=-cI*AAABR(70)+AAABR(71)
      MTR021(1,1)=AAABR(72)
      MTR021(1,2)=AAABR(76)
      MTR021(1,3)=AAABR(73)
      MTR021(2,1)=AAABR(74)
      MTR021(2,2)=AAABR(79)
      MTR021(2,3)=AAABR(75)
      MTR021(3,1)=AAABR(78)
      MTR021(3,2)=AAABR(77)
      MTR021(3,3)=AAABR(80)
      MTR022(1,1)=AAABR(91)-cI*AAABR(92)
      MTR022(1,2)=AAABR(81)-cI*AAABR(82)
      MTR022(1,3)=-AAABR(83)+cI*AAABR(84)
      MTR022(2,1)=AAABR(95)-cI*AAABR(96)
      MTR022(2,2)=AAABR(93)-cI*AAABR(94)
      MTR022(2,3)=-AAABR(97)+cI*AAABR(98)
      MTR022(3,1)=AAABR(89)-cI*AAABR(90)
      MTR022(3,2)=AAABR(87)-cI*AAABR(88)
      MTR022(3,3)=-AAABR(85)+cI*AAABR(86)
      MTR023(1)=AAABR(99)+cI*AAABR(100)
      MTR023(2)=AAABR(101)+cI*AAABR(102)
      MTR023(3)=-AAABR(103)+cI*AAABR(104)
      MTR024(1)=cI*AAABR(68)-AAABR(69)
      MTR024(2)=cI*AAABR(66)-AAABR(67)
      MTR024(3)=-cI*AAABR(70)-AAABR(71)
      MTR025(1)=AAABR(105)
      MTR025(2)=AAABR(107)
      MTR025(3)=AAABR(106)
      MTR026(1,1)=AAABR(91)+cI*AAABR(92)
      MTR026(1,2)=AAABR(81)+cI*AAABR(82)
      MTR026(1,3)=-AAABR(83)-cI*AAABR(84)
      MTR026(2,1)=AAABR(95)+cI*AAABR(96)
      MTR026(2,2)=AAABR(93)+cI*AAABR(94)
      MTR026(2,3)=-AAABR(97)-cI*AAABR(98)
      MTR026(3,1)=AAABR(89)+cI*AAABR(90)
      MTR026(3,2)=AAABR(87)+cI*AAABR(88)
      MTR026(3,3)=-AAABR(85)-cI*AAABR(86)
      MTR027(1,1)=AAABR(116)
      MTR027(1,2)=AAABR(109)
      MTR027(1,3)=AAABR(114)
      MTR027(2,1)=AAABR(113)
      MTR027(2,2)=AAABR(115)
      MTR027(2,3)=AAABR(108)
      MTR027(3,1)=AAABR(111)
      MTR027(3,2)=AAABR(112)
      MTR027(3,3)=AAABR(110)
      MTR028(1)=AAABR(59)-cI*AAABR(60)
      MTR028(2)=AAABR(57)-cI*AAABR(58)
      MTR028(3)=-AAABR(61)-cI*AAABR(62)
      MTR029(1)=AAABR(99)-cI*AAABR(100)
      MTR029(2)=AAABR(101)-cI*AAABR(102)
      MTR029(3)=-AAABR(103)-cI*AAABR(104)
      MTR030(1,1)=AAABR(118)
      MTR030(1,2)=AAABR(118)
      MTR030(1,3)=AAABR(118)
      MTR030(2,1)=AAABR(119)
      MTR030(2,2)=AAABR(119)
      MTR030(2,3)=AAABR(119)
      MTR030(3,1)=AAABR(117)
      MTR030(3,2)=AAABR(117)
      MTR030(3,3)=AAABR(117)
      MTR031(1)=AAABR(122)
      MTR031(2)=AAABR(120)+cI*AAABR(121)
      MTR031(3)=-AAABR(123)+cI*AAABR(124)
      MTR032(1)=AAABR(127)
      MTR032(2)=AAABR(128)+cI*AAABR(129)
      MTR032(3)=-AAABR(125)+cI*AAABR(126)
      MTR033(2)=AAABR(130)
      MTR033(3)=AAABR(131)
      MTR034(1)=AAABR(136)
      MTR034(2)=-cI**2*AAABR(132)-cI*AAABR(133)
      MTR034(3)=cI**2*AAABR(134)-cI*AAABR(135)
      MTR035(1,1)=AAABR(145)
      MTR035(1,2)=AAABR(137)
      MTR035(1,3)=AAABR(140)
      MTR035(2,1)=AAABR(138)
      MTR035(2,2)=AAABR(144)
      MTR035(2,3)=AAABR(143)
      MTR035(3,1)=AAABR(141)
      MTR035(3,2)=AAABR(139)
      MTR035(3,3)=AAABR(142)
      MTR036(1,1)=AAABR(162)-3D0*cI*AAABR(163)
      MTR036(1,2)=AAABR(148)-3D0*cI*AAABR(149)
      MTR036(1,3)=-AAABR(146)+cI*AAABR(147)
      MTR036(2,1)=AAABR(160)-3D0*cI*AAABR(161)
      MTR036(2,2)=AAABR(156)-3D0*cI*AAABR(157)
      MTR036(2,3)=-AAABR(158)+cI*AAABR(159)
      MTR036(3,1)=AAABR(154)-3D0*cI*AAABR(155)
      MTR036(3,2)=AAABR(150)-3D0*cI*AAABR(151)
      MTR036(3,3)=-AAABR(152)+cI*AAABR(153)
      MTR037(1)=AAABR(164)
      MTR037(2)=AAABR(165)+cI*AAABR(166)
      MTR037(3)=-AAABR(167)+cI*AAABR(168)
      MTR038(1)=AAABR(171)
      MTR038(2)=AAABR(172)+cI*AAABR(173)
      MTR038(3)=-AAABR(169)+cI*AAABR(170)
      MTR039(1)=AAABR(136)
      MTR039(2)=-cI**2*AAABR(132)+cI*AAABR(133)
      MTR039(3)=cI**2*AAABR(134)+cI*AAABR(135)
      MTR040(2)=AAABR(175)
      MTR040(3)=AAABR(174)
      MTR041(1,1)=AAABR(162)+3D0*cI*AAABR(163)
      MTR041(1,2)=AAABR(148)+3D0*cI*AAABR(149)
      MTR041(1,3)=-AAABR(146)-cI*AAABR(147)
      MTR041(2,1)=AAABR(160)+3D0*cI*AAABR(161)
      MTR041(2,2)=AAABR(156)+3D0*cI*AAABR(157)
      MTR041(2,3)=-AAABR(158)-cI*AAABR(159)
      MTR041(3,1)=AAABR(154)+3D0*cI*AAABR(155)
      MTR041(3,2)=AAABR(150)+3D0*cI*AAABR(151)
      MTR041(3,3)=-AAABR(152)-cI*AAABR(153)
      MTR042(1,1)=AAABR(180)
      MTR042(1,2)=AAABR(177)
      MTR042(1,3)=AAABR(184)
      MTR042(2,1)=AAABR(179)
      MTR042(2,2)=AAABR(181)
      MTR042(2,3)=AAABR(183)
      MTR042(3,1)=AAABR(178)
      MTR042(3,2)=AAABR(182)
      MTR042(3,3)=AAABR(176)
      MTR043(1)=AAABR(122)
      MTR043(2)=AAABR(120)-cI*AAABR(121)
      MTR043(3)=-AAABR(123)-cI*AAABR(124)
      MTR044(1)=AAABR(164)
      MTR044(2)=AAABR(165)-cI*AAABR(166)
      MTR044(3)=-AAABR(167)-cI*AAABR(168)
      MTR045(2)=AAABR(185)
      MTR045(3)=AAABR(186)
      MTR046(1)=AAABR(187)
      MTR046(2)=-cI**2*AAABR(188)+cI*AAABR(189)
      MTR046(3)=-cI**2*AAABR(190)+cI*AAABR(191)
      MTR047(1,1)=AAABR(192)
      MTR047(1,2)=AAABR(198)
      MTR047(1,3)=-AAABR(193)
      MTR047(2,1)=AAABR(195)
      MTR047(2,2)=AAABR(200)
      MTR047(2,3)=-AAABR(196)
      MTR047(3,1)=AAABR(199)
      MTR047(3,2)=AAABR(197)
      MTR047(3,3)=-AAABR(194)
      MTR048(1,1)=AAABR(203)-3D0*cI*AAABR(204)
      MTR048(1,2)=AAABR(213)-3D0*cI*AAABR(214)
      MTR048(1,3)=-AAABR(211)+cI*AAABR(212)
      MTR048(2,1)=AAABR(209)-3D0*cI*AAABR(210)
      MTR048(2,2)=AAABR(207)-3D0*cI*AAABR(208)
      MTR048(2,3)=-AAABR(215)+cI*AAABR(216)
      MTR048(3,1)=AAABR(205)-3D0*cI*AAABR(206)
      MTR048(3,2)=AAABR(217)-3D0*cI*AAABR(218)
      MTR048(3,3)=-AAABR(201)+cI*AAABR(202)
      MTR049(1)=AAABR(127)
      MTR049(2)=AAABR(128)-cI*AAABR(129)
      MTR049(3)=-AAABR(125)-cI*AAABR(126)
      MTR050(1)=AAABR(171)
      MTR050(2)=AAABR(172)-cI*AAABR(173)
      MTR050(3)=-AAABR(169)-cI*AAABR(170)
      MTR051(1)=AAABR(187)
      MTR051(2)=-cI**2*AAABR(188)-cI*AAABR(189)
      MTR051(3)=-cI**2*AAABR(190)-cI*AAABR(191)
      MTR052(2)=AAABR(220)
      MTR052(3)=AAABR(219)
      MTR053(1,1)=AAABR(203)+3D0*cI*AAABR(204)
      MTR053(1,2)=AAABR(213)+3D0*cI*AAABR(214)
      MTR053(1,3)=-AAABR(211)-cI*AAABR(212)
      MTR053(2,1)=AAABR(209)+3D0*cI*AAABR(210)
      MTR053(2,2)=AAABR(207)+3D0*cI*AAABR(208)
      MTR053(2,3)=-AAABR(215)-cI*AAABR(216)
      MTR053(3,1)=AAABR(205)+3D0*cI*AAABR(206)
      MTR053(3,2)=AAABR(217)+3D0*cI*AAABR(218)
      MTR053(3,3)=-AAABR(201)-cI*AAABR(202)
      MTR054(1,1)=AAABR(223)
      MTR054(1,2)=AAABR(222)
      MTR054(1,3)=-AAABR(227)
      MTR054(2,1)=AAABR(221)
      MTR054(2,2)=AAABR(228)
      MTR054(2,3)=-AAABR(226)
      MTR054(3,1)=AAABR(225)
      MTR054(3,2)=AAABR(224)
      MTR054(3,3)=-AAABR(229)
      MTR055(1)=AAABR(232)
      MTR055(2)=AAABR(230)
      MTR055(3)=AAABR(231)
      MTR056(1,1)=2D0*AAABR(235)
      MTR056(1,2)=AAABR(234)
      MTR056(1,3)=AAABR(237)
      MTR056(2,1)=-MTR056(1,2)
      MTR056(2,2)=2D0*AAABR(238)
      MTR056(2,3)=AAABR(233)
      MTR056(3,1)=-MTR056(1,3)
      MTR056(3,2)=-MTR056(2,3)
      MTR056(3,3)=2D0*AAABR(236)
      MTR057(1,1,1)=3D0*AAABR(245)
      MTR057(1,1,2)=AAABR(239)
      MTR057(1,1,3)=AAABR(247)
      MTR057(1,2,1)=MTR057(1,1,2)
      MTR057(1,2,2)=AAABR(240)
      MTR057(1,2,3)=AAABR(241)
      MTR057(1,3,1)=MTR057(1,1,3)
      MTR057(1,3,2)=MTR057(1,2,3)
      MTR057(1,3,3)=AAABR(242)
      MTR057(2,1,1)=MTR057(1,1,2)
      MTR057(2,1,2)=MTR057(1,2,2)
      MTR057(2,1,3)=MTR057(1,2,3)
      MTR057(2,2,1)=MTR057(1,2,2)
      MTR057(2,2,2)=3D0*AAABR(248)
      MTR057(2,2,3)=AAABR(246)
      MTR057(2,3,1)=MTR057(1,2,3)
      MTR057(2,3,2)=MTR057(2,2,3)
      MTR057(2,3,3)=AAABR(243)
      MTR057(3,1,1)=MTR057(1,1,3)
      MTR057(3,1,2)=MTR057(1,2,3)
      MTR057(3,1,3)=MTR057(1,3,3)
      MTR057(3,2,1)=MTR057(1,2,3)
      MTR057(3,2,2)=MTR057(2,2,3)
      MTR057(3,2,3)=MTR057(2,3,3)
      MTR057(3,3,1)=MTR057(1,3,3)
      MTR057(3,3,2)=MTR057(2,3,3)
      MTR057(3,3,3)=3D0*AAABR(244)
      MTR058(1)=AAABR(254)+cI*AAABR(255)
      MTR058(2)=AAABR(252)+cI*AAABR(253)
      MTR058(3)=AAABR(250)+cI*AAABR(251)
      MTR059(1)=AAABR(254)-cI*AAABR(255)
      MTR059(2)=AAABR(252)-cI*AAABR(253)
      MTR059(3)=AAABR(250)-cI*AAABR(251)
      MTR060(1)=AAABR(257)
      MTR060(2)=AAABR(258)
      MTR060(3)=AAABR(256)
      MTR061(1)=AAABR(260)
      MTR061(2)=AAABR(259)
      MTR061(3)=AAABR(261)-cI*AAABR(262)
      MTR062(1)=AAABR(265)
      MTR062(2)=AAABR(266)
      MTR062(3)=AAABR(263)-cI*AAABR(264)
      MTR063(1)=AAABR(260)
      MTR063(2)=AAABR(259)
      MTR063(3)=AAABR(261)+cI*AAABR(262)
      MTR064(1)=AAABR(267)
      MTR064(2)=AAABR(269)
      MTR064(3)=AAABR(268)
      MTR065(1)=AAABR(273)
      MTR065(2)=AAABR(272)
      MTR065(3)=AAABR(270)-cI*AAABR(271)
      MTR066(1)=AAABR(265)
      MTR066(2)=AAABR(266)
      MTR066(3)=AAABR(263)+cI*AAABR(264)
      MTR067(1)=AAABR(273)
      MTR067(2)=AAABR(272)
      MTR067(3)=AAABR(270)+cI*AAABR(271)
      MTR068(1)=AAABR(277)
      MTR068(2)=AAABR(276)
      MTR068(3)=AAABR(275)
      MTR069(1)=AAABR(281)
      MTR069(2)=AAABR(278)
      MTR069(3)=AAABR(279)-cI*AAABR(280)
      MTR070(1)=AAABR(285)
      MTR070(2)=AAABR(282)
      MTR070(3)=AAABR(283)-cI*AAABR(284)
      MTR071(1)=AAABR(287)
      MTR071(2)=AAABR(286)
      MTR071(3)=AAABR(288)-cI*AAABR(289)
      MTR072(1)=AAABR(281)
      MTR072(2)=AAABR(278)
      MTR072(3)=AAABR(279)+cI*AAABR(280)
      MTR073(1)=AAABR(292)
      MTR073(2)=AAABR(290)
      MTR073(3)=AAABR(291)
      MTR074(1)=AAABR(294)
      MTR074(2)=AAABR(293)
      MTR074(3)=AAABR(295)-cI*AAABR(296)
      MTR075(1)=AAABR(300)
      MTR075(2)=AAABR(297)
      MTR075(3)=AAABR(298)-cI*AAABR(299)
      MTR076(1)=AAABR(285)
      MTR076(2)=AAABR(282)
      MTR076(3)=AAABR(283)+cI*AAABR(284)
      MTR077(1)=AAABR(294)
      MTR077(2)=AAABR(293)
      MTR077(3)=AAABR(295)+cI*AAABR(296)
      MTR078(1)=AAABR(301)
      MTR078(2)=AAABR(303)
      MTR078(3)=AAABR(302)
      MTR079(1)=AAABR(304)
      MTR079(2)=AAABR(305)
      MTR079(3)=AAABR(306)-cI*AAABR(307)
      MTR080(1)=AAABR(287)
      MTR080(2)=AAABR(286)
      MTR080(3)=AAABR(288)+cI*AAABR(289)
      MTR081(1)=AAABR(300)
      MTR081(2)=AAABR(297)
      MTR081(3)=AAABR(298)+cI*AAABR(299)
      MTR082(1)=AAABR(304)
      MTR082(2)=AAABR(305)
      MTR082(3)=AAABR(306)+cI*AAABR(307)
      MTR083(1)=AAABR(309)
      MTR083(2)=AAABR(310)
      MTR083(3)=AAABR(308)
      MTR084(1)=AAABR(312)
      MTR084(2)=AAABR(313)
      MTR084(3)=AAABR(314)
      MTR085(1)=AAABR(317)
      MTR085(2)=AAABR(315)
      MTR085(3)=AAABR(316)
      MTR086(1,2)=AAABR(320)
      MTR086(1,3)=AAABR(319)
      MTR086(2,1)=-MTR086(1,2)
      MTR086(2,3)=AAABR(318)
      MTR086(3,1)=-MTR086(1,3)
      MTR086(3,2)=-MTR086(2,3)
      MTR087(1)=AAABR(322)
      MTR087(2)=AAABR(323)
      MTR087(3)=AAABR(324)
      MTR088(1,1)=AAABR(331)-cI*AAABR(332)
      MTR088(1,2)=-cI*AAABR(329)+AAABR(330)
      MTR088(2,1)=-cI*AAABR(327)+AAABR(328)
      MTR088(2,2)=AAABR(325)-cI*AAABR(326)
      MTR089(1,1)=cI*AAABR(332)+AAABR(331)
      MTR089(1,2)=cI*AAABR(327)+AAABR(328)
      MTR089(2,1)=cI*AAABR(329)+AAABR(330)
      MTR089(2,2)=cI*AAABR(326)+AAABR(325)
      MTR090(1,1,1)=AAABR(337)+cI*AAABR(338)
      MTR090(1,1,2)=AAABR(345)+cI*AAABR(346)
      MTR090(1,1,3)=AAABR(355)+cI*AAABR(356)
      MTR090(1,2,1)=AAABR(347)+cI*AAABR(348)
      MTR090(1,2,2)=AAABR(339)+cI*AAABR(340)
      MTR090(1,2,3)=AAABR(351)+cI*AAABR(352)
      MTR090(2,1,1)=AAABR(349)+cI*AAABR(350)
      MTR090(2,1,2)=AAABR(341)+cI*AAABR(342)
      MTR090(2,1,3)=AAABR(353)+cI*AAABR(354)
      MTR090(2,2,1)=AAABR(333)+cI*AAABR(334)
      MTR090(2,2,2)=AAABR(335)+cI*AAABR(336)
      MTR090(2,2,3)=AAABR(343)+cI*AAABR(344)
      MTR091(1,1,1)=cI*AAABR(338)-AAABR(337)
      MTR091(1,1,2)=cI*AAABR(346)-AAABR(345)
      MTR091(1,1,3)=cI*AAABR(356)-AAABR(355)
      MTR091(1,2,1)=-AAABR(349)+cI*AAABR(350)
      MTR091(1,2,2)=-AAABR(341)+cI*AAABR(342)
      MTR091(1,2,3)=-AAABR(353)+cI*AAABR(354)
      MTR091(2,1,1)=-AAABR(347)+cI*AAABR(348)
      MTR091(2,1,2)=-AAABR(339)+cI*AAABR(340)
      MTR091(2,1,3)=-AAABR(351)+cI*AAABR(352)
      MTR091(2,2,1)=cI*AAABR(334)-AAABR(333)
      MTR091(2,2,2)=cI*AAABR(336)-AAABR(335)
      MTR091(2,2,3)=cI*AAABR(344)-AAABR(343)
      MTR092(1,1)=AAABR(361)-cI*AAABR(362)
      MTR092(1,2)=AAABR(363)-cI*AAABR(364)
      MTR092(1,3)=-AAABR(365)-cI*AAABR(366)
      MTR092(2,1)=AAABR(359)-cI*AAABR(360)
      MTR092(2,2)=AAABR(357)-cI*AAABR(358)
      MTR092(2,3)=-AAABR(367)-cI*AAABR(368)
      MTR093(1,1)=AAABR(377)-cI*AAABR(378)
      MTR093(1,2)=AAABR(371)-cI*AAABR(372)
      MTR093(1,3)=-AAABR(375)-cI*AAABR(376)
      MTR093(2,1)=AAABR(373)-cI*AAABR(374)
      MTR093(2,2)=AAABR(369)-cI*AAABR(370)
      MTR093(2,3)=-AAABR(379)-cI*AAABR(380)
      MTR094(1,1)=AAABR(393)-cI*AAABR(394)
      MTR094(1,2)=AAABR(381)-cI*AAABR(382)
      MTR094(1,3)=AAABR(389)-cI*AAABR(390)
      MTR094(1,4)=AAABR(385)-cI*AAABR(386)
      MTR094(2,1)=AAABR(405)-cI*AAABR(406)
      MTR094(2,2)=AAABR(401)-cI*AAABR(402)
      MTR094(2,3)=AAABR(397)-cI*AAABR(398)
      MTR094(2,4)=AAABR(409)-cI*AAABR(410)
      MTR095(1,1)=AAABR(395)+cI*AAABR(396)
      MTR095(1,2)=AAABR(383)+cI*AAABR(384)
      MTR095(1,3)=AAABR(391)+cI*AAABR(392)
      MTR095(1,4)=AAABR(387)+cI*AAABR(388)
      MTR095(2,1)=AAABR(407)+cI*AAABR(408)
      MTR095(2,2)=AAABR(403)+cI*AAABR(404)
      MTR095(2,3)=AAABR(399)+cI*AAABR(400)
      MTR095(2,4)=AAABR(411)+cI*AAABR(412)
      MTR096(1,1)=AAABR(433)-cI*AAABR(434)
      MTR096(1,2)=AAABR(437)-cI*AAABR(438)
      MTR096(1,3)=AAABR(425)-cI*AAABR(426)
      MTR096(1,4)=AAABR(421)-cI*AAABR(422)
      MTR096(2,1)=AAABR(429)-cI*AAABR(430)
      MTR096(2,2)=AAABR(417)-cI*AAABR(418)
      MTR096(2,3)=AAABR(413)-cI*AAABR(414)
      MTR096(2,4)=AAABR(441)-cI*AAABR(442)
      MTR097(1,1)=AAABR(435)+cI*AAABR(436)
      MTR097(1,2)=AAABR(439)+cI*AAABR(440)
      MTR097(1,3)=AAABR(427)+cI*AAABR(428)
      MTR097(1,4)=AAABR(423)+cI*AAABR(424)
      MTR097(2,1)=AAABR(431)+cI*AAABR(432)
      MTR097(2,2)=AAABR(419)+cI*AAABR(420)
      MTR097(2,3)=AAABR(415)+cI*AAABR(416)
      MTR097(2,4)=AAABR(443)+cI*AAABR(444)
      MTR098(1,1)=2D0*AAABR(455)-2D0*cI*AAABR(456)
      MTR098(1,2)=AAABR(451)-cI*AAABR(452)
      MTR098(1,3)=-AAABR(445)-cI*AAABR(446)
      MTR098(2,1)=2D0*AAABR(449)-2D0*cI*AAABR(450)
      MTR098(2,2)=AAABR(457)-cI*AAABR(458)
      MTR098(2,3)=-AAABR(461)-cI*AAABR(462)
      MTR099(1,2)=AAABR(453)-cI*AAABR(454)
      MTR099(1,3)=-AAABR(447)-cI*AAABR(448)
      MTR099(2,2)=AAABR(459)-cI*AAABR(460)
      MTR099(2,3)=-AAABR(463)-cI*AAABR(464)
      MTR100(1,1)=2D0*AAABR(475)-2D0*cI*AAABR(476)
      MTR100(1,2)=AAABR(477)-cI*AAABR(478)
      MTR100(1,3)=-AAABR(471)-cI*AAABR(472)
      MTR100(2,1)=2D0*AAABR(465)-2D0*cI*AAABR(466)
      MTR100(2,2)=AAABR(481)-cI*AAABR(482)
      MTR100(2,3)=-AAABR(467)-cI*AAABR(468)
      MTR101(1,2)=AAABR(479)-cI*AAABR(480)
      MTR101(1,3)=-AAABR(473)-cI*AAABR(474)
      MTR101(2,2)=AAABR(483)-cI*AAABR(484)
      MTR101(2,3)=-AAABR(469)-cI*AAABR(470)
      MTR102(1)=AAABR(486)
      MTR102(2)=AAABR(485)
      MTR102(3)=AAABR(487)
      MTR103(1,1)=AAABR(498)-cI*AAABR(499)
      MTR103(1,2)=AAABR(500)-cI*AAABR(501)
      MTR103(1,3)=AAABR(494)-cI*AAABR(495)
      MTR103(2,1)=AAABR(488)-cI*AAABR(489)
      MTR103(2,2)=AAABR(496)-cI*AAABR(497)
      MTR103(2,3)=AAABR(492)-cI*AAABR(493)
      MTR103(3,1)=AAABR(502)-cI*AAABR(503)
      MTR103(3,2)=AAABR(504)-cI*AAABR(505)
      MTR103(3,3)=AAABR(490)-cI*AAABR(491)
      MTR104(1,1)=cI*AAABR(499)+AAABR(498)
      MTR104(1,2)=cI*AAABR(501)+AAABR(500)
      MTR104(1,3)=cI*AAABR(495)+AAABR(494)
      MTR104(2,1)=cI*AAABR(489)+AAABR(488)
      MTR104(2,2)=cI*AAABR(497)+AAABR(496)
      MTR104(2,3)=cI*AAABR(493)+AAABR(492)
      MTR104(3,1)=cI*AAABR(503)+AAABR(502)
      MTR104(3,2)=cI*AAABR(505)+AAABR(504)
      MTR104(3,3)=cI*AAABR(491)+AAABR(490)
      MTR105(1)=AAABR(506)
      MTR105(2)=AAABR(507)
      MTR105(3)=AAABR(508)
      MTR106(1)=AAABR(510)
      MTR106(2)=AAABR(511)
      MTR106(3)=AAABR(509)
      MTR107(1,1)=AAABR(516)-cI*AAABR(517)
      MTR107(1,2)=AAABR(556)-cI*AAABR(557)
      MTR107(1,3)=AAABR(532)+cI*AAABR(533)
      MTR107(2,1)=AAABR(512)-cI*AAABR(513)
      MTR107(2,2)=AAABR(552)-cI*AAABR(553)
      MTR107(2,3)=AAABR(540)+cI*AAABR(541)
      MTR107(3,1)=AAABR(524)-cI*AAABR(525)
      MTR107(3,2)=AAABR(536)-cI*AAABR(537)
      MTR107(3,3)=AAABR(544)+cI*AAABR(545)
      MTR107(4,1)=AAABR(548)-cI*AAABR(549)
      MTR107(4,2)=AAABR(520)-cI*AAABR(521)
      MTR107(4,3)=AAABR(528)+cI*AAABR(529)
      MTR108(1,1)=AAABR(518)+cI*AAABR(519)
      MTR108(1,2)=AAABR(558)+cI*AAABR(559)
      MTR108(1,3)=-AAABR(534)+cI*AAABR(535)
      MTR108(2,1)=AAABR(514)+cI*AAABR(515)
      MTR108(2,2)=AAABR(554)+cI*AAABR(555)
      MTR108(2,3)=-AAABR(542)+cI*AAABR(543)
      MTR108(3,1)=AAABR(526)+cI*AAABR(527)
      MTR108(3,2)=AAABR(538)+cI*AAABR(539)
      MTR108(3,3)=-AAABR(546)+cI*AAABR(547)
      MTR108(4,1)=AAABR(550)+cI*AAABR(551)
      MTR108(4,2)=AAABR(522)+cI*AAABR(523)
      MTR108(4,3)=-AAABR(530)+cI*AAABR(531)
      MTR109(1,1)=AAABR(600)-cI*AAABR(601)
      MTR109(1,2)=AAABR(576)-cI*AAABR(577)
      MTR109(1,3)=AAABR(604)+cI*AAABR(605)
      MTR109(2,1)=AAABR(592)-cI*AAABR(593)
      MTR109(2,2)=AAABR(588)-cI*AAABR(589)
      MTR109(2,3)=AAABR(596)+cI*AAABR(597)
      MTR109(3,1)=AAABR(560)-cI*AAABR(561)
      MTR109(3,2)=AAABR(584)-cI*AAABR(585)
      MTR109(3,3)=AAABR(568)+cI*AAABR(569)
      MTR109(4,1)=AAABR(572)-cI*AAABR(573)
      MTR109(4,2)=AAABR(580)-cI*AAABR(581)
      MTR109(4,3)=AAABR(564)+cI*AAABR(565)
      MTR110(1,1)=AAABR(602)+cI*AAABR(603)
      MTR110(1,2)=AAABR(578)+cI*AAABR(579)
      MTR110(1,3)=-AAABR(606)+cI*AAABR(607)
      MTR110(2,1)=AAABR(594)+cI*AAABR(595)
      MTR110(2,2)=AAABR(590)+cI*AAABR(591)
      MTR110(2,3)=-AAABR(598)+cI*AAABR(599)
      MTR110(3,1)=AAABR(562)+cI*AAABR(563)
      MTR110(3,2)=AAABR(586)+cI*AAABR(587)
      MTR110(3,3)=-AAABR(570)+cI*AAABR(571)
      MTR110(4,1)=AAABR(574)+cI*AAABR(575)
      MTR110(4,2)=AAABR(582)+cI*AAABR(583)
      MTR110(4,3)=-AAABR(566)+cI*AAABR(567)
      MTR111(1,1)=AAABR(608)-cI*AAABR(609)
      MTR111(1,2)=AAABR(608)-cI*AAABR(609)
      MTR111(1,3)=AAABR(608)-cI*AAABR(609)
      MTR111(2,1)=AAABR(616)-cI*AAABR(617)
      MTR111(2,2)=AAABR(616)-cI*AAABR(617)
      MTR111(2,3)=AAABR(616)-cI*AAABR(617)
      MTR112(1,1)=AAABR(610)-cI*AAABR(611)
      MTR112(1,2)=AAABR(614)-cI*AAABR(615)
      MTR112(1,3)=AAABR(612)-cI*AAABR(613)
      MTR112(2,1)=AAABR(620)-cI*AAABR(621)
      MTR112(2,2)=AAABR(622)-cI*AAABR(623)
      MTR112(2,3)=AAABR(618)-cI*AAABR(619)
      MTR113(1,1)=AAABR(361)+cI*AAABR(362)
      MTR113(1,2)=AAABR(363)+cI*AAABR(364)
      MTR113(1,3)=-AAABR(365)+cI*AAABR(366)
      MTR113(2,1)=AAABR(359)+cI*AAABR(360)
      MTR113(2,2)=AAABR(357)+cI*AAABR(358)
      MTR113(2,3)=-AAABR(367)+cI*AAABR(368)
      MTR114(1,1)=AAABR(377)+cI*AAABR(378)
      MTR114(1,2)=AAABR(371)+cI*AAABR(372)
      MTR114(1,3)=-AAABR(375)+cI*AAABR(376)
      MTR114(2,1)=AAABR(373)+cI*AAABR(374)
      MTR114(2,2)=AAABR(369)+cI*AAABR(370)
      MTR114(2,3)=-AAABR(379)+cI*AAABR(380)
      MTR115(1,1)=AAABR(630)-cI*AAABR(631)
      MTR115(1,2)=AAABR(630)-cI*AAABR(631)
      MTR115(1,3)=AAABR(630)-cI*AAABR(631)
      MTR115(2,1)=AAABR(628)-cI*AAABR(629)
      MTR115(2,2)=AAABR(628)-cI*AAABR(629)
      MTR115(2,3)=AAABR(628)-cI*AAABR(629)
      MTR115(3,1)=AAABR(626)-cI*AAABR(627)
      MTR115(3,2)=AAABR(626)-cI*AAABR(627)
      MTR115(3,3)=AAABR(626)-cI*AAABR(627)
      MTR115(4,1)=AAABR(624)-cI*AAABR(625)
      MTR115(4,2)=AAABR(624)-cI*AAABR(625)
      MTR115(4,3)=AAABR(624)-cI*AAABR(625)
      MTR116(1,1)=AAABR(648)-cI*AAABR(649)
      MTR116(1,2)=AAABR(664)-cI*AAABR(665)
      MTR116(1,3)=AAABR(632)+cI*AAABR(633)
      MTR116(2,1)=AAABR(636)-cI*AAABR(637)
      MTR116(2,2)=AAABR(676)-cI*AAABR(677)
      MTR116(2,3)=AAABR(668)+cI*AAABR(669)
      MTR116(3,1)=AAABR(640)-cI*AAABR(641)
      MTR116(3,2)=AAABR(656)-cI*AAABR(657)
      MTR116(3,3)=AAABR(660)+cI*AAABR(661)
      MTR116(4,1)=AAABR(652)-cI*AAABR(653)
      MTR116(4,2)=AAABR(672)-cI*AAABR(673)
      MTR116(4,3)=AAABR(644)+cI*AAABR(645)
      MTR117(1,1)=2D0*AAABR(650)+2D0*cI*AAABR(651)
      MTR117(1,2)=AAABR(666)+cI*AAABR(667)
      MTR117(1,3)=-AAABR(634)+cI*AAABR(635)
      MTR117(2,1)=2D0*AAABR(638)+2D0*cI*AAABR(639)
      MTR117(2,2)=AAABR(678)+cI*AAABR(679)
      MTR117(2,3)=-AAABR(670)+cI*AAABR(671)
      MTR117(3,1)=2D0*AAABR(642)+2D0*cI*AAABR(643)
      MTR117(3,2)=AAABR(658)+cI*AAABR(659)
      MTR117(3,3)=-AAABR(662)+cI*AAABR(663)
      MTR117(4,1)=2D0*AAABR(654)+2D0*cI*AAABR(655)
      MTR117(4,2)=AAABR(674)+cI*AAABR(675)
      MTR117(4,3)=-AAABR(646)+cI*AAABR(647)
      MTR118(1,1)=AAABR(704)-cI*AAABR(705)
      MTR118(1,2)=AAABR(712)-cI*AAABR(713)
      MTR118(1,3)=AAABR(696)+cI*AAABR(697)
      MTR118(2,1)=AAABR(708)-cI*AAABR(709)
      MTR118(2,2)=AAABR(680)-cI*AAABR(681)
      MTR118(2,3)=AAABR(700)+cI*AAABR(701)
      MTR118(3,1)=AAABR(688)-cI*AAABR(689)
      MTR118(3,2)=AAABR(716)-cI*AAABR(717)
      MTR118(3,3)=AAABR(720)+cI*AAABR(721)
      MTR118(4,1)=AAABR(724)-cI*AAABR(725)
      MTR118(4,2)=AAABR(692)-cI*AAABR(693)
      MTR118(4,3)=AAABR(684)+cI*AAABR(685)
      MTR119(1,1)=2D0*AAABR(706)+2D0*cI*AAABR(707)
      MTR119(1,2)=AAABR(714)+cI*AAABR(715)
      MTR119(1,3)=-AAABR(698)+cI*AAABR(699)
      MTR119(2,1)=2D0*AAABR(710)+2D0*cI*AAABR(711)
      MTR119(2,2)=AAABR(682)+cI*AAABR(683)
      MTR119(2,3)=-AAABR(702)+cI*AAABR(703)
      MTR119(3,1)=2D0*AAABR(690)+2D0*cI*AAABR(691)
      MTR119(3,2)=AAABR(718)+cI*AAABR(719)
      MTR119(3,3)=-AAABR(722)+cI*AAABR(723)
      MTR119(4,1)=2D0*AAABR(726)+2D0*cI*AAABR(727)
      MTR119(4,2)=AAABR(694)+cI*AAABR(695)
      MTR119(4,3)=-AAABR(686)+cI*AAABR(687)
      MTR120(2)=AAABR(728)
      MTR120(3)=AAABR(729)
      MTR121(2,1)=AAABR(730)-cI*AAABR(731)
      MTR121(2,2)=AAABR(734)-cI*AAABR(735)
      MTR121(2,3)=AAABR(732)-cI*AAABR(733)
      MTR121(3,1)=AAABR(736)-cI*AAABR(737)
      MTR121(3,2)=AAABR(740)-cI*AAABR(741)
      MTR121(3,3)=AAABR(738)-cI*AAABR(739)
      MTR122(2,1)=cI*AAABR(731)+AAABR(730)
      MTR122(2,2)=cI*AAABR(735)+AAABR(734)
      MTR122(2,3)=cI*AAABR(733)+AAABR(732)
      MTR122(3,1)=cI*AAABR(737)+AAABR(736)
      MTR122(3,2)=cI*AAABR(741)+AAABR(740)
      MTR122(3,3)=cI*AAABR(739)+AAABR(738)
      MTR123(2)=AAABR(742)
      MTR123(3)=AAABR(744)
      MTR124(2)=AAABR(743)
      MTR124(3)=AAABR(745)
      MTR125(2)=AAABR(746)
      MTR125(3)=AAABR(748)
      MTR126(2)=AAABR(747)
      MTR126(3)=AAABR(749)
      MTR127(1,1)=2D0*AAABR(758)-2D0*cI*AAABR(759)
      MTR127(1,2)=AAABR(750)-cI*AAABR(751)
      MTR127(1,3)=AAABR(754)-cI*AAABR(755)
      MTR127(2,1)=2D0*AAABR(774)-2D0*cI*AAABR(775)
      MTR127(2,2)=AAABR(776)-cI*AAABR(777)
      MTR127(2,3)=AAABR(770)-cI*AAABR(771)
      MTR128(1,2)=AAABR(752)-cI*AAABR(753)
      MTR128(1,3)=AAABR(756)-cI*AAABR(757)
      MTR128(2,2)=AAABR(778)-cI*AAABR(779)
      MTR128(2,3)=AAABR(772)-cI*AAABR(773)
      MTR129(1,1)=2D0*AAABR(768)-2D0*cI*AAABR(769)
      MTR129(1,2)=AAABR(764)-cI*AAABR(765)
      MTR129(1,3)=AAABR(760)-cI*AAABR(761)
      MTR129(2,1)=2D0*AAABR(788)-2D0*cI*AAABR(789)
      MTR129(2,2)=AAABR(784)-cI*AAABR(785)
      MTR129(2,3)=AAABR(780)-cI*AAABR(781)
      MTR130(1,2)=AAABR(766)-cI*AAABR(767)
      MTR130(1,3)=AAABR(762)-cI*AAABR(763)
      MTR130(2,2)=AAABR(786)-cI*AAABR(787)
      MTR130(2,3)=AAABR(782)-cI*AAABR(783)
      MTR131(1)=AAABR(790)
      MTR131(2)=AAABR(792)
      MTR131(3)=AAABR(794)+cI*AAABR(795)
      MTR132(1)=AAABR(791)
      MTR132(2)=AAABR(793)
      MTR132(3)=AAABR(796)+cI*AAABR(797)
      MTR133(1)=AAABR(798)
      MTR133(2)=AAABR(804)
      MTR133(3)=AAABR(800)+cI*AAABR(801)
      MTR134(1)=AAABR(799)
      MTR134(2)=AAABR(805)
      MTR134(3)=AAABR(802)+cI*AAABR(803)
      MTR135(1,1)=2D0*AAABR(455)+2D0*cI*AAABR(456)
      MTR135(1,2)=AAABR(451)+cI*AAABR(452)
      MTR135(1,3)=-AAABR(445)+cI*AAABR(446)
      MTR135(2,1)=2D0*AAABR(449)+2D0*cI*AAABR(450)
      MTR135(2,2)=AAABR(457)+cI*AAABR(458)
      MTR135(2,3)=-AAABR(461)+cI*AAABR(462)
      MTR136(1,2)=AAABR(453)+cI*AAABR(454)
      MTR136(1,3)=-AAABR(447)+cI*AAABR(448)
      MTR136(2,2)=AAABR(459)+cI*AAABR(460)
      MTR136(2,3)=-AAABR(463)+cI*AAABR(464)
      MTR137(1,1)=2D0*AAABR(475)+2D0*cI*AAABR(476)
      MTR137(1,2)=AAABR(477)+cI*AAABR(478)
      MTR137(1,3)=-AAABR(471)+cI*AAABR(472)
      MTR137(2,1)=2D0*AAABR(465)+2D0*cI*AAABR(466)
      MTR137(2,2)=AAABR(481)+cI*AAABR(482)
      MTR137(2,3)=-AAABR(467)+cI*AAABR(468)
      MTR138(1,2)=AAABR(479)+cI*AAABR(480)
      MTR138(1,3)=-AAABR(473)+cI*AAABR(474)
      MTR138(2,2)=AAABR(483)+cI*AAABR(484)
      MTR138(2,3)=-AAABR(469)+cI*AAABR(470)
      MTR139(1,1)=AAABR(826)-cI*AAABR(827)
      MTR139(1,2)=AAABR(810)-cI*AAABR(811)
      MTR139(1,3)=AAABR(806)-cI*AAABR(807)
      MTR139(2,1)=AAABR(834)-cI*AAABR(835)
      MTR139(2,2)=AAABR(838)-cI*AAABR(839)
      MTR139(2,3)=AAABR(814)-cI*AAABR(815)
      MTR139(3,1)=AAABR(822)-cI*AAABR(823)
      MTR139(3,2)=AAABR(850)-cI*AAABR(851)
      MTR139(3,3)=AAABR(842)-cI*AAABR(843)
      MTR139(4,1)=AAABR(818)-cI*AAABR(819)
      MTR139(4,2)=AAABR(830)-cI*AAABR(831)
      MTR139(4,3)=AAABR(846)-cI*AAABR(847)
      MTR140(1,1)=4D0*AAABR(828)+4D0*cI*AAABR(829)
      MTR140(1,2)=AAABR(812)+cI*AAABR(813)
      MTR140(1,3)=-AAABR(808)+cI*AAABR(809)
      MTR140(2,1)=4D0*AAABR(836)+4D0*cI*AAABR(837)
      MTR140(2,2)=AAABR(840)+cI*AAABR(841)
      MTR140(2,3)=-AAABR(816)+cI*AAABR(817)
      MTR140(3,1)=4D0*AAABR(824)+4D0*cI*AAABR(825)
      MTR140(3,2)=AAABR(852)+cI*AAABR(853)
      MTR140(3,3)=-AAABR(844)+cI*AAABR(845)
      MTR140(4,1)=4D0*AAABR(820)+4D0*cI*AAABR(821)
      MTR140(4,2)=AAABR(832)+cI*AAABR(833)
      MTR140(4,3)=-AAABR(848)+cI*AAABR(849)
      MTR141(1,1)=AAABR(890)-cI*AAABR(891)
      MTR141(1,2)=AAABR(886)-cI*AAABR(887)
      MTR141(1,3)=AAABR(878)-cI*AAABR(879)
      MTR141(2,1)=AAABR(862)-cI*AAABR(863)
      MTR141(2,2)=AAABR(894)-cI*AAABR(895)
      MTR141(2,3)=AAABR(882)-cI*AAABR(883)
      MTR141(3,1)=AAABR(854)-cI*AAABR(855)
      MTR141(3,2)=AAABR(858)-cI*AAABR(859)
      MTR141(3,3)=AAABR(870)-cI*AAABR(871)
      MTR141(4,1)=AAABR(874)-cI*AAABR(875)
      MTR141(4,2)=AAABR(866)-cI*AAABR(867)
      MTR141(4,3)=AAABR(898)-cI*AAABR(899)
      MTR142(1,1)=4D0*AAABR(892)+4D0*cI*AAABR(893)
      MTR142(1,2)=AAABR(888)+cI*AAABR(889)
      MTR142(1,3)=-AAABR(880)+cI*AAABR(881)
      MTR142(2,1)=4D0*AAABR(864)+4D0*cI*AAABR(865)
      MTR142(2,2)=AAABR(896)+cI*AAABR(897)
      MTR142(2,3)=-AAABR(884)+cI*AAABR(885)
      MTR142(3,1)=4D0*AAABR(856)+4D0*cI*AAABR(857)
      MTR142(3,2)=AAABR(860)+cI*AAABR(861)
      MTR142(3,3)=-AAABR(872)+cI*AAABR(873)
      MTR142(4,1)=4D0*AAABR(876)+4D0*cI*AAABR(877)
      MTR142(4,2)=AAABR(868)+cI*AAABR(869)
      MTR142(4,3)=-AAABR(900)+cI*AAABR(901)
      MTR143(2)=AAABR(903)
      MTR143(3)=AAABR(902)
      MTR144(2,1)=AAABR(906)-cI*AAABR(907)
      MTR144(2,2)=AAABR(908)-cI*AAABR(909)
      MTR144(2,3)=AAABR(912)-cI*AAABR(913)
      MTR144(3,1)=AAABR(914)-cI*AAABR(915)
      MTR144(3,2)=AAABR(910)-cI*AAABR(911)
      MTR144(3,3)=AAABR(904)-cI*AAABR(905)
      MTR145(2,1)=cI*AAABR(907)+AAABR(906)
      MTR145(2,2)=cI*AAABR(909)+AAABR(908)
      MTR145(2,3)=cI*AAABR(913)+AAABR(912)
      MTR145(3,1)=cI*AAABR(915)+AAABR(914)
      MTR145(3,2)=cI*AAABR(911)+AAABR(910)
      MTR145(3,3)=cI*AAABR(905)+AAABR(904)
      MTR146(1)=AAABR(916)
      MTR146(2)=AAABR(922)
      MTR146(3)=AAABR(918)+cI*AAABR(919)
      MTR147(1)=AAABR(917)
      MTR147(2)=AAABR(923)
      MTR147(3)=AAABR(920)+cI*AAABR(921)
      MTR148(1)=AAABR(930)
      MTR148(2)=AAABR(928)
      MTR148(3)=AAABR(924)+cI*AAABR(925)
      MTR149(1)=AAABR(931)
      MTR149(2)=AAABR(929)
      MTR149(3)=AAABR(926)+cI*AAABR(927)
      MTR150(1,1)=2D0*AAABR(758)+2D0*cI*AAABR(759)
      MTR150(1,2)=AAABR(750)+cI*AAABR(751)
      MTR150(1,3)=AAABR(754)+cI*AAABR(755)
      MTR150(2,1)=2D0*AAABR(774)+2D0*cI*AAABR(775)
      MTR150(2,2)=AAABR(776)+cI*AAABR(777)
      MTR150(2,3)=AAABR(770)+cI*AAABR(771)
      MTR151(1,2)=AAABR(752)+cI*AAABR(753)
      MTR151(1,3)=AAABR(756)+cI*AAABR(757)
      MTR151(2,2)=AAABR(778)+cI*AAABR(779)
      MTR151(2,3)=AAABR(772)+cI*AAABR(773)
      MTR152(1,1)=2D0*AAABR(768)+2D0*cI*AAABR(769)
      MTR152(1,2)=AAABR(764)+cI*AAABR(765)
      MTR152(1,3)=AAABR(760)+cI*AAABR(761)
      MTR152(2,1)=2D0*AAABR(788)+2D0*cI*AAABR(789)
      MTR152(2,2)=AAABR(784)+cI*AAABR(785)
      MTR152(2,3)=AAABR(780)+cI*AAABR(781)
      MTR153(1,2)=AAABR(766)+cI*AAABR(767)
      MTR153(1,3)=AAABR(762)+cI*AAABR(763)
      MTR153(2,2)=AAABR(786)+cI*AAABR(787)
      MTR153(2,3)=AAABR(782)+cI*AAABR(783)
      MTR154(1,1)=AAABR(648)+cI*AAABR(649)
      MTR154(1,2)=AAABR(664)+cI*AAABR(665)
      MTR154(1,3)=AAABR(632)-cI*AAABR(633)
      MTR154(2,1)=AAABR(636)+cI*AAABR(637)
      MTR154(2,2)=AAABR(676)+cI*AAABR(677)
      MTR154(2,3)=AAABR(668)-cI*AAABR(669)
      MTR154(3,1)=AAABR(640)+cI*AAABR(641)
      MTR154(3,2)=AAABR(656)+cI*AAABR(657)
      MTR154(3,3)=AAABR(660)-cI*AAABR(661)
      MTR154(4,1)=AAABR(652)+cI*AAABR(653)
      MTR154(4,2)=AAABR(672)+cI*AAABR(673)
      MTR154(4,3)=AAABR(644)-cI*AAABR(645)
      MTR155(1,1)=2D0*AAABR(650)-2D0*cI*AAABR(651)
      MTR155(1,2)=AAABR(666)-cI*AAABR(667)
      MTR155(1,3)=-AAABR(634)-cI*AAABR(635)
      MTR155(2,1)=2D0*AAABR(638)-2D0*cI*AAABR(639)
      MTR155(2,2)=AAABR(678)-cI*AAABR(679)
      MTR155(2,3)=-AAABR(670)-cI*AAABR(671)
      MTR155(3,1)=2D0*AAABR(642)-2D0*cI*AAABR(643)
      MTR155(3,2)=AAABR(658)-cI*AAABR(659)
      MTR155(3,3)=-AAABR(662)-cI*AAABR(663)
      MTR155(4,1)=2D0*AAABR(654)-2D0*cI*AAABR(655)
      MTR155(4,2)=AAABR(674)-cI*AAABR(675)
      MTR155(4,3)=-AAABR(646)-cI*AAABR(647)
      MTR156(1,1)=AAABR(704)+cI*AAABR(705)
      MTR156(1,2)=AAABR(712)+cI*AAABR(713)
      MTR156(1,3)=AAABR(696)-cI*AAABR(697)
      MTR156(2,1)=AAABR(708)+cI*AAABR(709)
      MTR156(2,2)=AAABR(680)+cI*AAABR(681)
      MTR156(2,3)=AAABR(700)-cI*AAABR(701)
      MTR156(3,1)=AAABR(688)+cI*AAABR(689)
      MTR156(3,2)=AAABR(716)+cI*AAABR(717)
      MTR156(3,3)=AAABR(720)-cI*AAABR(721)
      MTR156(4,1)=AAABR(724)+cI*AAABR(725)
      MTR156(4,2)=AAABR(692)+cI*AAABR(693)
      MTR156(4,3)=AAABR(684)-cI*AAABR(685)
      MTR157(1,1)=2D0*AAABR(706)-2D0*cI*AAABR(707)
      MTR157(1,2)=AAABR(714)-cI*AAABR(715)
      MTR157(1,3)=-AAABR(698)-cI*AAABR(699)
      MTR157(2,1)=2D0*AAABR(710)-2D0*cI*AAABR(711)
      MTR157(2,2)=AAABR(682)-cI*AAABR(683)
      MTR157(2,3)=-AAABR(702)-cI*AAABR(703)
      MTR157(3,1)=2D0*AAABR(690)-2D0*cI*AAABR(691)
      MTR157(3,2)=AAABR(718)-cI*AAABR(719)
      MTR157(3,3)=-AAABR(722)-cI*AAABR(723)
      MTR157(4,1)=2D0*AAABR(726)-2D0*cI*AAABR(727)
      MTR157(4,2)=AAABR(694)-cI*AAABR(695)
      MTR157(4,3)=-AAABR(686)-cI*AAABR(687)
      MTR158(1)=AAABR(790)
      MTR158(2)=AAABR(792)
      MTR158(3)=AAABR(794)-cI*AAABR(795)
      MTR159(1)=AAABR(791)
      MTR159(2)=AAABR(793)
      MTR159(3)=AAABR(796)-cI*AAABR(797)
      MTR160(1)=AAABR(798)
      MTR160(2)=AAABR(804)
      MTR160(3)=AAABR(800)-cI*AAABR(801)
      MTR161(1)=AAABR(799)
      MTR161(2)=AAABR(805)
      MTR161(3)=AAABR(802)-cI*AAABR(803)
      MTR162(1,1)=AAABR(826)+cI*AAABR(827)
      MTR162(1,2)=AAABR(810)+cI*AAABR(811)
      MTR162(1,3)=AAABR(806)+cI*AAABR(807)
      MTR162(2,1)=AAABR(834)+cI*AAABR(835)
      MTR162(2,2)=AAABR(838)+cI*AAABR(839)
      MTR162(2,3)=AAABR(814)+cI*AAABR(815)
      MTR162(3,1)=AAABR(822)+cI*AAABR(823)
      MTR162(3,2)=AAABR(850)+cI*AAABR(851)
      MTR162(3,3)=AAABR(842)+cI*AAABR(843)
      MTR162(4,1)=AAABR(818)+cI*AAABR(819)
      MTR162(4,2)=AAABR(830)+cI*AAABR(831)
      MTR162(4,3)=AAABR(846)+cI*AAABR(847)
      MTR163(1,1)=4D0*AAABR(828)-4D0*cI*AAABR(829)
      MTR163(1,2)=AAABR(812)-cI*AAABR(813)
      MTR163(1,3)=-AAABR(808)-cI*AAABR(809)
      MTR163(2,1)=4D0*AAABR(836)-4D0*cI*AAABR(837)
      MTR163(2,2)=AAABR(840)-cI*AAABR(841)
      MTR163(2,3)=-AAABR(816)-cI*AAABR(817)
      MTR163(3,1)=4D0*AAABR(824)-4D0*cI*AAABR(825)
      MTR163(3,2)=AAABR(852)-cI*AAABR(853)
      MTR163(3,3)=-AAABR(844)-cI*AAABR(845)
      MTR163(4,1)=4D0*AAABR(820)-4D0*cI*AAABR(821)
      MTR163(4,2)=AAABR(832)-cI*AAABR(833)
      MTR163(4,3)=-AAABR(848)-cI*AAABR(849)
      MTR164(1,1)=AAABR(890)+cI*AAABR(891)
      MTR164(1,2)=AAABR(886)+cI*AAABR(887)
      MTR164(1,3)=AAABR(878)+cI*AAABR(879)
      MTR164(2,1)=AAABR(862)+cI*AAABR(863)
      MTR164(2,2)=AAABR(894)+cI*AAABR(895)
      MTR164(2,3)=AAABR(882)+cI*AAABR(883)
      MTR164(3,1)=AAABR(854)+cI*AAABR(855)
      MTR164(3,2)=AAABR(858)+cI*AAABR(859)
      MTR164(3,3)=AAABR(870)+cI*AAABR(871)
      MTR164(4,1)=AAABR(874)+cI*AAABR(875)
      MTR164(4,2)=AAABR(866)+cI*AAABR(867)
      MTR164(4,3)=AAABR(898)+cI*AAABR(899)
      MTR165(1,1)=4D0*AAABR(892)-4D0*cI*AAABR(893)
      MTR165(1,2)=AAABR(888)-cI*AAABR(889)
      MTR165(1,3)=-AAABR(880)-cI*AAABR(881)
      MTR165(2,1)=4D0*AAABR(864)-4D0*cI*AAABR(865)
      MTR165(2,2)=AAABR(896)-cI*AAABR(897)
      MTR165(2,3)=-AAABR(884)-cI*AAABR(885)
      MTR165(3,1)=4D0*AAABR(856)-4D0*cI*AAABR(857)
      MTR165(3,2)=AAABR(860)-cI*AAABR(861)
      MTR165(3,3)=-AAABR(872)-cI*AAABR(873)
      MTR165(4,1)=4D0*AAABR(876)-4D0*cI*AAABR(877)
      MTR165(4,2)=AAABR(868)-cI*AAABR(869)
      MTR165(4,3)=-AAABR(900)-cI*AAABR(901)
      MTR166(1)=AAABR(916)
      MTR166(2)=AAABR(922)
      MTR166(3)=AAABR(918)-cI*AAABR(919)
      MTR167(1)=AAABR(917)
      MTR167(2)=AAABR(923)
      MTR167(3)=AAABR(920)-cI*AAABR(921)
      MTR168(1)=AAABR(930)
      MTR168(2)=AAABR(928)
      MTR168(3)=AAABR(924)-cI*AAABR(925)
      MTR169(1)=AAABR(931)
      MTR169(2)=AAABR(929)
      MTR169(3)=AAABR(926)-cI*AAABR(927)
      MTR170(1,1)=AAABR(608)+cI*AAABR(609)
      MTR170(1,2)=AAABR(608)+cI*AAABR(609)
      MTR170(1,3)=AAABR(608)+cI*AAABR(609)
      MTR170(2,1)=AAABR(616)+cI*AAABR(617)
      MTR170(2,2)=AAABR(616)+cI*AAABR(617)
      MTR170(2,3)=AAABR(616)+cI*AAABR(617)
      MTR171(1,1)=AAABR(610)+cI*AAABR(611)
      MTR171(1,2)=AAABR(614)+cI*AAABR(615)
      MTR171(1,3)=AAABR(612)+cI*AAABR(613)
      MTR171(2,1)=AAABR(620)+cI*AAABR(621)
      MTR171(2,2)=AAABR(622)+cI*AAABR(623)
      MTR171(2,3)=AAABR(618)+cI*AAABR(619)
      MTR172(1,1)=AAABR(516)+cI*AAABR(517)
      MTR172(1,2)=AAABR(556)+cI*AAABR(557)
      MTR172(1,3)=AAABR(532)-cI*AAABR(533)
      MTR172(2,1)=AAABR(512)+cI*AAABR(513)
      MTR172(2,2)=AAABR(552)+cI*AAABR(553)
      MTR172(2,3)=AAABR(540)-cI*AAABR(541)
      MTR172(3,1)=AAABR(524)+cI*AAABR(525)
      MTR172(3,2)=AAABR(536)+cI*AAABR(537)
      MTR172(3,3)=AAABR(544)-cI*AAABR(545)
      MTR172(4,1)=AAABR(548)+cI*AAABR(549)
      MTR172(4,2)=AAABR(520)+cI*AAABR(521)
      MTR172(4,3)=AAABR(528)-cI*AAABR(529)
      MTR173(1,1)=AAABR(518)-cI*AAABR(519)
      MTR173(1,2)=AAABR(558)-cI*AAABR(559)
      MTR173(1,3)=-AAABR(534)-cI*AAABR(535)
      MTR173(2,1)=AAABR(514)-cI*AAABR(515)
      MTR173(2,2)=AAABR(554)-cI*AAABR(555)
      MTR173(2,3)=-AAABR(542)-cI*AAABR(543)
      MTR173(3,1)=AAABR(526)-cI*AAABR(527)
      MTR173(3,2)=AAABR(538)-cI*AAABR(539)
      MTR173(3,3)=-AAABR(546)-cI*AAABR(547)
      MTR173(4,1)=AAABR(550)-cI*AAABR(551)
      MTR173(4,2)=AAABR(522)-cI*AAABR(523)
      MTR173(4,3)=-AAABR(530)-cI*AAABR(531)
      MTR174(1,1)=AAABR(600)+cI*AAABR(601)
      MTR174(1,2)=AAABR(576)+cI*AAABR(577)
      MTR174(1,3)=AAABR(604)-cI*AAABR(605)
      MTR174(2,1)=AAABR(592)+cI*AAABR(593)
      MTR174(2,2)=AAABR(588)+cI*AAABR(589)
      MTR174(2,3)=AAABR(596)-cI*AAABR(597)
      MTR174(3,1)=AAABR(560)+cI*AAABR(561)
      MTR174(3,2)=AAABR(584)+cI*AAABR(585)
      MTR174(3,3)=AAABR(568)-cI*AAABR(569)
      MTR174(4,1)=AAABR(572)+cI*AAABR(573)
      MTR174(4,2)=AAABR(580)+cI*AAABR(581)
      MTR174(4,3)=AAABR(564)-cI*AAABR(565)
      MTR175(1,1)=AAABR(602)-cI*AAABR(603)
      MTR175(1,2)=AAABR(578)-cI*AAABR(579)
      MTR175(1,3)=-AAABR(606)-cI*AAABR(607)
      MTR175(2,1)=AAABR(594)-cI*AAABR(595)
      MTR175(2,2)=AAABR(590)-cI*AAABR(591)
      MTR175(2,3)=-AAABR(598)-cI*AAABR(599)
      MTR175(3,1)=AAABR(562)-cI*AAABR(563)
      MTR175(3,2)=AAABR(586)-cI*AAABR(587)
      MTR175(3,3)=-AAABR(570)-cI*AAABR(571)
      MTR175(4,1)=AAABR(574)-cI*AAABR(575)
      MTR175(4,2)=AAABR(582)-cI*AAABR(583)
      MTR175(4,3)=-AAABR(566)-cI*AAABR(567)
      MTR176(1,1)=AAABR(630)+cI*AAABR(631)
      MTR176(1,2)=AAABR(630)+cI*AAABR(631)
      MTR176(1,3)=AAABR(630)+cI*AAABR(631)
      MTR176(2,1)=AAABR(628)+cI*AAABR(629)
      MTR176(2,2)=AAABR(628)+cI*AAABR(629)
      MTR176(2,3)=AAABR(628)+cI*AAABR(629)
      MTR176(3,1)=AAABR(626)+cI*AAABR(627)
      MTR176(3,2)=AAABR(626)+cI*AAABR(627)
      MTR176(3,3)=AAABR(626)+cI*AAABR(627)
      MTR176(4,1)=AAABR(624)+cI*AAABR(625)
      MTR176(4,2)=AAABR(624)+cI*AAABR(625)
      MTR176(4,3)=AAABR(624)+cI*AAABR(625)
      MTR177(1,1)=AAABR(393)+cI*AAABR(394)
      MTR177(1,2)=AAABR(381)+cI*AAABR(382)
      MTR177(1,3)=AAABR(389)+cI*AAABR(390)
      MTR177(1,4)=AAABR(385)+cI*AAABR(386)
      MTR177(2,1)=AAABR(405)+cI*AAABR(406)
      MTR177(2,2)=AAABR(401)+cI*AAABR(402)
      MTR177(2,3)=AAABR(397)+cI*AAABR(398)
      MTR177(2,4)=AAABR(409)+cI*AAABR(410)
      MTR178(1,1)=AAABR(395)-cI*AAABR(396)
      MTR178(1,2)=AAABR(383)-cI*AAABR(384)
      MTR178(1,3)=AAABR(391)-cI*AAABR(392)
      MTR178(1,4)=AAABR(387)-cI*AAABR(388)
      MTR178(2,1)=AAABR(407)-cI*AAABR(408)
      MTR178(2,2)=AAABR(403)-cI*AAABR(404)
      MTR178(2,3)=AAABR(399)-cI*AAABR(400)
      MTR178(2,4)=AAABR(411)-cI*AAABR(412)
      MTR179(1,1)=AAABR(433)+cI*AAABR(434)
      MTR179(1,2)=AAABR(437)+cI*AAABR(438)
      MTR179(1,3)=AAABR(425)+cI*AAABR(426)
      MTR179(1,4)=AAABR(421)+cI*AAABR(422)
      MTR179(2,1)=AAABR(429)+cI*AAABR(430)
      MTR179(2,2)=AAABR(417)+cI*AAABR(418)
      MTR179(2,3)=AAABR(413)+cI*AAABR(414)
      MTR179(2,4)=AAABR(441)+cI*AAABR(442)
      MTR180(1,1)=AAABR(435)-cI*AAABR(436)
      MTR180(1,2)=AAABR(439)-cI*AAABR(440)
      MTR180(1,3)=AAABR(427)-cI*AAABR(428)
      MTR180(1,4)=AAABR(423)-cI*AAABR(424)
      MTR180(2,1)=AAABR(431)-cI*AAABR(432)
      MTR180(2,2)=AAABR(419)-cI*AAABR(420)
      MTR180(2,3)=AAABR(415)-cI*AAABR(416)
      MTR180(2,4)=AAABR(443)-cI*AAABR(444)
      MTR181(1,1)=2D0*AAABR(932)-2D0*cI*AAABR(933)
      MTR181(1,2)=AAABR(936)-cI*AAABR(937)
      MTR181(1,3)=AAABR(938)-cI*AAABR(939)
      MTR181(1,4)=AAABR(934)-cI*AAABR(935)
      MTR181(2,1)=AAABR(936)-cI*AAABR(937)
      MTR181(2,2)=2D0*AAABR(968)-2D0*cI*AAABR(969)
      MTR181(2,3)=AAABR(964)-cI*AAABR(965)
      MTR181(2,4)=AAABR(966)-cI*AAABR(967)
      MTR181(3,1)=AAABR(938)-cI*AAABR(939)
      MTR181(3,2)=AAABR(964)-cI*AAABR(965)
      MTR181(3,3)=2D0*AAABR(988)-2D0*cI*AAABR(989)
      MTR181(3,4)=AAABR(990)-cI*AAABR(991)
      MTR181(4,1)=AAABR(934)-cI*AAABR(935)
      MTR181(4,2)=AAABR(966)-cI*AAABR(967)
      MTR181(4,3)=AAABR(990)-cI*AAABR(991)
      MTR181(4,4)=2D0*AAABR(1004)-2D0*cI*AAABR(1005)
      MTR182(1,1)=2D0*cI*AAABR(933)+2D0*AAABR(932)
      MTR182(1,2)=cI*AAABR(937)+AAABR(936)
      MTR182(1,3)=cI*AAABR(939)+AAABR(938)
      MTR182(1,4)=cI*AAABR(935)+AAABR(934)
      MTR182(2,1)=cI*AAABR(937)+AAABR(936)
      MTR182(2,2)=2D0*cI*AAABR(969)+2D0*AAABR(968)
      MTR182(2,3)=cI*AAABR(965)+AAABR(964)
      MTR182(2,4)=cI*AAABR(967)+AAABR(966)
      MTR182(3,1)=cI*AAABR(939)+AAABR(938)
      MTR182(3,2)=cI*AAABR(965)+AAABR(964)
      MTR182(3,3)=2D0*cI*AAABR(989)+2D0*AAABR(988)
      MTR182(3,4)=cI*AAABR(991)+AAABR(990)
      MTR182(4,1)=cI*AAABR(935)+AAABR(934)
      MTR182(4,2)=cI*AAABR(967)+AAABR(966)
      MTR182(4,3)=cI*AAABR(991)+AAABR(990)
      MTR182(4,4)=2D0*cI*AAABR(1005)+2D0*AAABR(1004)
      MTR183(1,1,1)=2D0*AAABR(958)+2D0*cI*AAABR(959)
      MTR183(1,1,2)=2D0*AAABR(962)+2D0*cI*AAABR(963)
      MTR183(1,1,3)=2D0*AAABR(952)+2D0*cI*AAABR(953)
      MTR183(1,2,1)=AAABR(950)+cI*AAABR(951)
      MTR183(1,2,2)=AAABR(946)+cI*AAABR(947)
      MTR183(1,2,3)=AAABR(954)+cI*AAABR(955)
      MTR183(1,3,1)=AAABR(942)+cI*AAABR(943)
      MTR183(1,3,2)=AAABR(956)+cI*AAABR(957)
      MTR183(1,3,3)=AAABR(940)+cI*AAABR(941)
      MTR183(1,4,1)=AAABR(960)+cI*AAABR(961)
      MTR183(1,4,2)=AAABR(944)+cI*AAABR(945)
      MTR183(1,4,3)=AAABR(948)+cI*AAABR(949)
      MTR183(2,1,1)=AAABR(950)+cI*AAABR(951)
      MTR183(2,1,2)=AAABR(946)+cI*AAABR(947)
      MTR183(2,1,3)=AAABR(954)+cI*AAABR(955)
      MTR183(2,2,1)=2D0*AAABR(984)+2D0*cI*AAABR(985)
      MTR183(2,2,2)=2D0*AAABR(982)+2D0*cI*AAABR(983)
      MTR183(2,2,3)=2D0*AAABR(980)+2D0*cI*AAABR(981)
      MTR183(2,3,1)=AAABR(978)+cI*AAABR(979)
      MTR183(2,3,2)=AAABR(972)+cI*AAABR(973)
      MTR183(2,3,3)=AAABR(986)+cI*AAABR(987)
      MTR183(2,4,1)=AAABR(974)+cI*AAABR(975)
      MTR183(2,4,2)=AAABR(976)+cI*AAABR(977)
      MTR183(2,4,3)=AAABR(970)+cI*AAABR(971)
      MTR183(3,1,1)=AAABR(942)+cI*AAABR(943)
      MTR183(3,1,2)=AAABR(956)+cI*AAABR(957)
      MTR183(3,1,3)=AAABR(940)+cI*AAABR(941)
      MTR183(3,2,1)=AAABR(978)+cI*AAABR(979)
      MTR183(3,2,2)=AAABR(972)+cI*AAABR(973)
      MTR183(3,2,3)=AAABR(986)+cI*AAABR(987)
      MTR183(3,3,1)=2D0*AAABR(994)+2D0*cI*AAABR(995)
      MTR183(3,3,2)=2D0*AAABR(992)+2D0*cI*AAABR(993)
      MTR183(3,3,3)=2D0*AAABR(1000)+2D0*cI*AAABR(1001)
      MTR183(3,4,1)=AAABR(998)+cI*AAABR(999)
      MTR183(3,4,2)=AAABR(1002)+cI*AAABR(1003)
      MTR183(3,4,3)=AAABR(996)+cI*AAABR(997)
      MTR183(4,1,1)=AAABR(960)+cI*AAABR(961)
      MTR183(4,1,2)=AAABR(944)+cI*AAABR(945)
      MTR183(4,1,3)=AAABR(948)+cI*AAABR(949)
      MTR183(4,2,1)=AAABR(974)+cI*AAABR(975)
      MTR183(4,2,2)=AAABR(976)+cI*AAABR(977)
      MTR183(4,2,3)=AAABR(970)+cI*AAABR(971)
      MTR183(4,3,1)=AAABR(998)+cI*AAABR(999)
      MTR183(4,3,2)=AAABR(1002)+cI*AAABR(1003)
      MTR183(4,3,3)=AAABR(996)+cI*AAABR(997)
      MTR183(4,4,1)=2D0*AAABR(1008)+2D0*cI*AAABR(1009)
      MTR183(4,4,2)=2D0*AAABR(1006)+2D0*cI*AAABR(1007)
      MTR183(4,4,3)=2D0*AAABR(1010)+2D0*cI*AAABR(1011)
      MTR184(1,1,1)=2D0*cI*AAABR(959)-2D0*AAABR(958)
      MTR184(1,1,2)=2D0*cI*AAABR(963)-2D0*AAABR(962)
      MTR184(1,1,3)=2D0*cI*AAABR(953)-2D0*AAABR(952)
      MTR184(1,2,1)=cI*AAABR(951)-AAABR(950)
      MTR184(1,2,2)=cI*AAABR(947)-AAABR(946)
      MTR184(1,2,3)=cI*AAABR(955)-AAABR(954)
      MTR184(1,3,1)=cI*AAABR(943)-AAABR(942)
      MTR184(1,3,2)=cI*AAABR(957)-AAABR(956)
      MTR184(1,3,3)=cI*AAABR(941)-AAABR(940)
      MTR184(1,4,1)=cI*AAABR(961)-AAABR(960)
      MTR184(1,4,2)=cI*AAABR(945)-AAABR(944)
      MTR184(1,4,3)=cI*AAABR(949)-AAABR(948)
      MTR184(2,1,1)=cI*AAABR(951)-AAABR(950)
      MTR184(2,1,2)=cI*AAABR(947)-AAABR(946)
      MTR184(2,1,3)=cI*AAABR(955)-AAABR(954)
      MTR184(2,2,1)=2D0*cI*AAABR(985)-2D0*AAABR(984)
      MTR184(2,2,2)=2D0*cI*AAABR(983)-2D0*AAABR(982)
      MTR184(2,2,3)=2D0*cI*AAABR(981)-2D0*AAABR(980)
      MTR184(2,3,1)=cI*AAABR(979)-AAABR(978)
      MTR184(2,3,2)=cI*AAABR(973)-AAABR(972)
      MTR184(2,3,3)=cI*AAABR(987)-AAABR(986)
      MTR184(2,4,1)=cI*AAABR(975)-AAABR(974)
      MTR184(2,4,2)=cI*AAABR(977)-AAABR(976)
      MTR184(2,4,3)=cI*AAABR(971)-AAABR(970)
      MTR184(3,1,1)=cI*AAABR(943)-AAABR(942)
      MTR184(3,1,2)=cI*AAABR(957)-AAABR(956)
      MTR184(3,1,3)=cI*AAABR(941)-AAABR(940)
      MTR184(3,2,1)=cI*AAABR(979)-AAABR(978)
      MTR184(3,2,2)=cI*AAABR(973)-AAABR(972)
      MTR184(3,2,3)=cI*AAABR(987)-AAABR(986)
      MTR184(3,3,1)=2D0*cI*AAABR(995)-2D0*AAABR(994)
      MTR184(3,3,2)=2D0*cI*AAABR(993)-2D0*AAABR(992)
      MTR184(3,3,3)=2D0*cI*AAABR(1001)-2D0*AAABR(1000)
      MTR184(3,4,1)=cI*AAABR(999)-AAABR(998)
      MTR184(3,4,2)=cI*AAABR(1003)-AAABR(1002)
      MTR184(3,4,3)=cI*AAABR(997)-AAABR(996)
      MTR184(4,1,1)=cI*AAABR(961)-AAABR(960)
      MTR184(4,1,2)=cI*AAABR(945)-AAABR(944)
      MTR184(4,1,3)=cI*AAABR(949)-AAABR(948)
      MTR184(4,2,1)=cI*AAABR(975)-AAABR(974)
      MTR184(4,2,2)=cI*AAABR(977)-AAABR(976)
      MTR184(4,2,3)=cI*AAABR(971)-AAABR(970)
      MTR184(4,3,1)=cI*AAABR(999)-AAABR(998)
      MTR184(4,3,2)=cI*AAABR(1003)-AAABR(1002)
      MTR184(4,3,3)=cI*AAABR(997)-AAABR(996)
      MTR184(4,4,1)=2D0*cI*AAABR(1009)-2D0*AAABR(1008)
      MTR184(4,4,2)=2D0*cI*AAABR(1007)-2D0*AAABR(1006)
      MTR184(4,4,3)=2D0*cI*AAABR(1011)-2D0*AAABR(1010)
      MTR185(1,1)=AAABR(1016)
      MTR185(1,2)=AAABR(1012)-cI*AAABR(1013)
      MTR185(2,1)=AAABR(1012)+cI*AAABR(1013)
      MTR185(2,2)=AAABR(1018)
      MTR186(1,1)=AAABR(1017)
      MTR186(1,2)=AAABR(1014)-cI*AAABR(1015)
      MTR186(2,1)=AAABR(1014)+cI*AAABR(1015)
      MTR186(2,2)=AAABR(1019)
      MTR187(1,1)=AAABR(1048)+cI*AAABR(1049)
      MTR187(1,2)=AAABR(1032)+cI*AAABR(1033)
      MTR187(1,3)=AAABR(1044)+cI*AAABR(1045)
      MTR187(1,4)=AAABR(1040)+cI*AAABR(1041)
      MTR187(2,1)=AAABR(1028)+cI*AAABR(1029)
      MTR187(2,2)=AAABR(1036)+cI*AAABR(1037)
      MTR187(2,3)=AAABR(1024)+cI*AAABR(1025)
      MTR187(2,4)=AAABR(1020)+cI*AAABR(1021)
      MTR188(1,1)=AAABR(1050)+cI*AAABR(1051)
      MTR188(1,2)=AAABR(1034)+cI*AAABR(1035)
      MTR188(1,3)=AAABR(1046)+cI*AAABR(1047)
      MTR188(1,4)=AAABR(1042)+cI*AAABR(1043)
      MTR188(2,1)=AAABR(1030)+cI*AAABR(1031)
      MTR188(2,2)=AAABR(1038)+cI*AAABR(1039)
      MTR188(2,3)=AAABR(1026)+cI*AAABR(1027)
      MTR188(2,4)=AAABR(1022)+cI*AAABR(1023)
      MTR189(1,1)=AAABR(1048)-cI*AAABR(1049)
      MTR189(1,2)=AAABR(1032)-cI*AAABR(1033)
      MTR189(1,3)=AAABR(1044)-cI*AAABR(1045)
      MTR189(1,4)=AAABR(1040)-cI*AAABR(1041)
      MTR189(2,1)=AAABR(1028)-cI*AAABR(1029)
      MTR189(2,2)=AAABR(1036)-cI*AAABR(1037)
      MTR189(2,3)=AAABR(1024)-cI*AAABR(1025)
      MTR189(2,4)=AAABR(1020)-cI*AAABR(1021)
      MTR190(1,1)=AAABR(1050)-cI*AAABR(1051)
      MTR190(1,2)=AAABR(1034)-cI*AAABR(1035)
      MTR190(1,3)=AAABR(1046)-cI*AAABR(1047)
      MTR190(1,4)=AAABR(1042)-cI*AAABR(1043)
      MTR190(2,1)=AAABR(1030)-cI*AAABR(1031)
      MTR190(2,2)=AAABR(1038)-cI*AAABR(1039)
      MTR190(2,3)=AAABR(1026)-cI*AAABR(1027)
      MTR190(2,4)=AAABR(1022)-cI*AAABR(1023)
      MTR191(1,1)=AAABR(1062)
      MTR191(1,2)=AAABR(1057)-cI*AAABR(1056)
      MTR191(1,3)=AAABR(1059)-cI*AAABR(1058)
      MTR191(1,4)=AAABR(1061)-cI*AAABR(1060)
      MTR191(2,1)=cI*AAABR(1056)+AAABR(1057)
      MTR191(2,2)=AAABR(1063)
      MTR191(2,3)=AAABR(1067)-cI*AAABR(1066)
      MTR191(2,4)=AAABR(1065)-cI*AAABR(1064)
      MTR191(3,1)=cI*AAABR(1058)+AAABR(1059)
      MTR191(3,2)=cI*AAABR(1066)+AAABR(1067)
      MTR191(3,3)=AAABR(1070)
      MTR191(3,4)=AAABR(1069)-cI*AAABR(1068)
      MTR191(4,1)=cI*AAABR(1060)+AAABR(1061)
      MTR191(4,2)=cI*AAABR(1064)+AAABR(1065)
      MTR191(4,3)=cI*AAABR(1068)+AAABR(1069)
      MTR191(4,4)=AAABR(1071)
      MTR192(1,1)=AAABR(1062)
      MTR192(1,2)=cI*AAABR(1056)+AAABR(1057)
      MTR192(1,3)=cI*AAABR(1058)+AAABR(1059)
      MTR192(1,4)=cI*AAABR(1060)+AAABR(1061)
      MTR192(2,1)=AAABR(1057)-cI*AAABR(1056)
      MTR192(2,2)=AAABR(1063)
      MTR192(2,3)=cI*AAABR(1066)+AAABR(1067)
      MTR192(2,4)=cI*AAABR(1064)+AAABR(1065)
      MTR192(3,1)=AAABR(1059)-cI*AAABR(1058)
      MTR192(3,2)=AAABR(1067)-cI*AAABR(1066)
      MTR192(3,3)=AAABR(1070)
      MTR192(3,4)=cI*AAABR(1068)+AAABR(1069)
      MTR192(4,1)=AAABR(1061)-cI*AAABR(1060)
      MTR192(4,2)=AAABR(1065)-cI*AAABR(1064)
      MTR192(4,3)=AAABR(1069)-cI*AAABR(1068)
      MTR192(4,4)=AAABR(1071)
      MTR193(1)=AAABR(1081)+cI*AAABR(1083)
      MTR193(2)=AAABR(1077)+cI*AAABR(1079)
      MTR193(3)=AAABR(1085)+cI*AAABR(1087)
      MTR194(1)=AAABR(1080)+cI*AAABR(1082)
      MTR194(2)=AAABR(1076)+cI*AAABR(1078)
      MTR194(3)=AAABR(1084)+cI*AAABR(1086)
      MTR195(1)=AAABR(1081)-cI*AAABR(1083)
      MTR195(2)=AAABR(1077)-cI*AAABR(1079)
      MTR195(3)=AAABR(1085)-cI*AAABR(1087)
      MTR196(1)=AAABR(1080)-cI*AAABR(1082)
      MTR196(2)=AAABR(1076)-cI*AAABR(1078)
      MTR196(3)=AAABR(1084)-cI*AAABR(1086)
      MTR197(1)=AAABR(1088)
      MTR197(2)=AAABR(1090)
      MTR197(3)=AAABR(1089)
      MTR198(1)=AAABR(1096)
      MTR198(2)=AAABR(1095)
      MTR198(3)=AAABR(1094)
      MTR199(1)=AAABR(1093)
      MTR199(2)=AAABR(1091)
      MTR199(3)=AAABR(1092)
      MTR200(1)=AAABR(1097)
      MTR200(2)=AAABR(1100)
      MTR200(3)=AAABR(1098)-cI*AAABR(1099)
      MTR201(1)=AAABR(1105)
      MTR201(2)=AAABR(1106)
      MTR201(3)=AAABR(1107)-cI*AAABR(1108)
      MTR202(1)=AAABR(1104)
      MTR202(2)=AAABR(1101)
      MTR202(3)=AAABR(1102)-cI*AAABR(1103)
      MTR203(1)=AAABR(1116)
      MTR203(2)=AAABR(1110)
      MTR203(3)=AAABR(1113)-cI*AAABR(1114)
      MTR204(1)=AAABR(1115)
      MTR204(2)=AAABR(1109)
      MTR204(3)=AAABR(1111)-cI*AAABR(1112)
      MTR205(1)=AAABR(1097)
      MTR205(2)=AAABR(1100)
      MTR205(3)=AAABR(1098)+cI*AAABR(1099)
      MTR206(1)=AAABR(1105)
      MTR206(2)=AAABR(1106)
      MTR206(3)=AAABR(1107)+cI*AAABR(1108)
      MTR207(1)=AAABR(1104)
      MTR207(2)=AAABR(1101)
      MTR207(3)=AAABR(1102)+cI*AAABR(1103)
      MTR208(1)=AAABR(1119)
      MTR208(2)=AAABR(1117)
      MTR208(3)=AAABR(1118)
      MTR209(1)=AAABR(1124)
      MTR209(2)=AAABR(1125)
      MTR209(3)=AAABR(1123)
      MTR210(1)=AAABR(1121)
      MTR210(2)=AAABR(1122)
      MTR210(3)=AAABR(1120)
      MTR211(1)=AAABR(1127)
      MTR211(2)=AAABR(1129)
      MTR211(3)=AAABR(1132)-cI*AAABR(1133)
      MTR212(1)=AAABR(1126)
      MTR212(2)=AAABR(1128)
      MTR212(3)=AAABR(1130)-cI*AAABR(1131)
      MTR213(1)=AAABR(1116)
      MTR213(2)=AAABR(1110)
      MTR213(3)=AAABR(1113)+cI*AAABR(1114)
      MTR214(1)=AAABR(1115)
      MTR214(2)=AAABR(1109)
      MTR214(3)=AAABR(1111)+cI*AAABR(1112)
      MTR215(1)=AAABR(1127)
      MTR215(2)=AAABR(1129)
      MTR215(3)=AAABR(1132)+cI*AAABR(1133)
      MTR216(1)=AAABR(1126)
      MTR216(2)=AAABR(1128)
      MTR216(3)=AAABR(1130)+cI*AAABR(1131)
      MTR217(1)=AAABR(1135)
      MTR217(2)=AAABR(1137)
      MTR217(3)=AAABR(1136)
      MTR218(1)=AAABR(1140)
      MTR218(2)=AAABR(1141)
      MTR218(3)=AAABR(1138)
      MTR219(1)=AAABR(1148)
      MTR219(2)=AAABR(1147)
      MTR219(3)=AAABR(1146)
      MTR220(1)=AAABR(1143)
      MTR220(2)=AAABR(1142)
      MTR220(3)=AAABR(1144)
      MTR221(1)=AAABR(1151)
      MTR221(2)=AAABR(1152)
      MTR221(3)=AAABR(1149)-cI*AAABR(1150)
      MTR222(1)=AAABR(1156)
      MTR222(2)=AAABR(1155)
      MTR222(3)=AAABR(1153)-cI*AAABR(1154)
      MTR223(1)=AAABR(1163)
      MTR223(2)=AAABR(1164)
      MTR223(3)=AAABR(1161)-cI*AAABR(1162)
      MTR224(1)=AAABR(1160)
      MTR224(2)=AAABR(1159)
      MTR224(3)=AAABR(1157)-cI*AAABR(1158)
      MTR225(1)=AAABR(1166)
      MTR225(2)=AAABR(1168)
      MTR225(3)=AAABR(1171)-cI*AAABR(1172)
      MTR226(1)=AAABR(1176)
      MTR226(2)=AAABR(1173)
      MTR226(3)=AAABR(1174)-cI*AAABR(1175)
      MTR227(1)=AAABR(1165)
      MTR227(2)=AAABR(1167)
      MTR227(3)=AAABR(1169)-cI*AAABR(1170)
      MTR228(1)=AAABR(1178)
      MTR228(2)=AAABR(1184)
      MTR228(3)=AAABR(1181)-cI*AAABR(1182)
      MTR229(1)=AAABR(1185)
      MTR229(2)=AAABR(1188)
      MTR229(3)=AAABR(1186)-cI*AAABR(1187)
      MTR230(1)=AAABR(1177)
      MTR230(2)=AAABR(1183)
      MTR230(3)=AAABR(1179)-cI*AAABR(1180)
      MTR231(1)=AAABR(1151)
      MTR231(2)=AAABR(1152)
      MTR231(3)=AAABR(1149)+cI*AAABR(1150)
      MTR232(1)=AAABR(1156)
      MTR232(2)=AAABR(1155)
      MTR232(3)=AAABR(1153)+cI*AAABR(1154)
      MTR233(1)=AAABR(1163)
      MTR233(2)=AAABR(1164)
      MTR233(3)=AAABR(1161)+cI*AAABR(1162)
      MTR234(1)=AAABR(1160)
      MTR234(2)=AAABR(1159)
      MTR234(3)=AAABR(1157)+cI*AAABR(1158)
      MTR235(1)=AAABR(1189)
      MTR235(2)=AAABR(1191)
      MTR235(3)=AAABR(1190)
      MTR236(1)=AAABR(1192)
      MTR236(2)=AAABR(1193)
      MTR236(3)=AAABR(1194)
      MTR237(1)=AAABR(1199)
      MTR237(2)=AAABR(1198)
      MTR237(3)=AAABR(1200)
      MTR238(1)=AAABR(1195)
      MTR238(2)=AAABR(1197)
      MTR238(3)=AAABR(1196)
      MTR239(1)=AAABR(1206)
      MTR239(2)=AAABR(1208)
      MTR239(3)=AAABR(1203)-cI*AAABR(1204)
      MTR240(1)=AAABR(1212)
      MTR240(2)=AAABR(1211)
      MTR240(3)=AAABR(1209)-cI*AAABR(1210)
      MTR241(1)=AAABR(1205)
      MTR241(2)=AAABR(1207)
      MTR241(3)=AAABR(1201)-cI*AAABR(1202)
      MTR242(1)=AAABR(1220)
      MTR242(2)=AAABR(1218)
      MTR242(3)=AAABR(1215)-cI*AAABR(1216)
      MTR243(1)=AAABR(1221)
      MTR243(2)=AAABR(1222)
      MTR243(3)=AAABR(1223)-cI*AAABR(1224)
      MTR244(1)=AAABR(1219)
      MTR244(2)=AAABR(1217)
      MTR244(3)=AAABR(1213)-cI*AAABR(1214)
      MTR245(1)=AAABR(1166)
      MTR245(2)=AAABR(1168)
      MTR245(3)=AAABR(1171)+cI*AAABR(1172)
      MTR246(1)=AAABR(1176)
      MTR246(2)=AAABR(1173)
      MTR246(3)=AAABR(1174)+cI*AAABR(1175)
      MTR247(1)=AAABR(1165)
      MTR247(2)=AAABR(1167)
      MTR247(3)=AAABR(1169)+cI*AAABR(1170)
      MTR248(1)=AAABR(1206)
      MTR248(2)=AAABR(1208)
      MTR248(3)=AAABR(1203)+cI*AAABR(1204)
      MTR249(1)=AAABR(1212)
      MTR249(2)=AAABR(1211)
      MTR249(3)=AAABR(1209)+cI*AAABR(1210)
      MTR250(1)=AAABR(1205)
      MTR250(2)=AAABR(1207)
      MTR250(3)=AAABR(1201)+cI*AAABR(1202)
      MTR251(1)=AAABR(1225)
      MTR251(2)=AAABR(1227)
      MTR251(3)=AAABR(1226)
      MTR252(1)=AAABR(1228)
      MTR252(2)=AAABR(1230)
      MTR252(3)=AAABR(1229)
      MTR253(1)=AAABR(1235)
      MTR253(2)=AAABR(1236)
      MTR253(3)=AAABR(1234)
      MTR254(1)=AAABR(1232)
      MTR254(2)=AAABR(1231)
      MTR254(3)=AAABR(1233)
      MTR255(1)=AAABR(1239)
      MTR255(2)=AAABR(1240)
      MTR255(3)=AAABR(1237)-cI*AAABR(1238)
      MTR256(1)=AAABR(1244)
      MTR256(2)=AAABR(1243)
      MTR256(3)=AAABR(1241)-cI*AAABR(1242)
      MTR257(1)=AAABR(1249)
      MTR257(2)=AAABR(1250)
      MTR257(3)=AAABR(1251)-cI*AAABR(1252)
      MTR258(1)=AAABR(1246)
      MTR258(2)=AAABR(1245)
      MTR258(3)=AAABR(1247)-cI*AAABR(1248)
      MTR259(1)=AAABR(1178)
      MTR259(2)=AAABR(1184)
      MTR259(3)=AAABR(1181)+cI*AAABR(1182)
      MTR260(1)=AAABR(1185)
      MTR260(2)=AAABR(1188)
      MTR260(3)=AAABR(1186)+cI*AAABR(1187)
      MTR261(1)=AAABR(1177)
      MTR261(2)=AAABR(1183)
      MTR261(3)=AAABR(1179)+cI*AAABR(1180)
      MTR262(1)=AAABR(1220)
      MTR262(2)=AAABR(1218)
      MTR262(3)=AAABR(1215)+cI*AAABR(1216)
      MTR263(1)=AAABR(1221)
      MTR263(2)=AAABR(1222)
      MTR263(3)=AAABR(1223)+cI*AAABR(1224)
      MTR264(1)=AAABR(1219)
      MTR264(2)=AAABR(1217)
      MTR264(3)=AAABR(1213)+cI*AAABR(1214)
      MTR265(1)=AAABR(1239)
      MTR265(2)=AAABR(1240)
      MTR265(3)=AAABR(1237)+cI*AAABR(1238)
      MTR266(1)=AAABR(1244)
      MTR266(2)=AAABR(1243)
      MTR266(3)=AAABR(1241)+cI*AAABR(1242)
      MTR267(1)=AAABR(1249)
      MTR267(2)=AAABR(1250)
      MTR267(3)=AAABR(1251)+cI*AAABR(1252)
      MTR268(1)=AAABR(1246)
      MTR268(2)=AAABR(1245)
      MTR268(3)=AAABR(1247)+cI*AAABR(1248)
      MTR269(1)=AAABR(1254)
      MTR269(2)=AAABR(1253)
      MTR269(3)=AAABR(1255)
      MTR270(1)=AAABR(1258)
      MTR270(2)=AAABR(1256)
      MTR270(3)=AAABR(1257)
      MTR271(1)=AAABR(1262)
      MTR271(2)=AAABR(1264)
      MTR271(3)=AAABR(1263)
      MTR272(1)=AAABR(1260)
      MTR272(2)=AAABR(1261)
      MTR272(3)=AAABR(1259)
      MTR273(1)=AAABR(1270)
      MTR273(2)=AAABR(1272)
      MTR273(3)=AAABR(1268)
      MTR274(1)=AAABR(1269)
      MTR274(2)=AAABR(1271)
      MTR274(3)=AAABR(1267)

      end

***********************************************

      subroutine ModelVarIni(fail, sqrtS)
      implicit none
      double precision sqrtS
      integer fail
      double precision Alfas

#include "model.h"

c      double precision ALPHAS2
c      external ALPHAS2

c      Alfas = ALPHAS2(sqrtS)
c      GG = sqrt(4*pi*Alfas)
      
      end

************************************************

      subroutine ModelDigest
      implicit none

#include "model.h"
#ifdef TEST
      write(16,*) 'pi=',pi
      write(16,*) 'EE=',EE
      write(16,*) 'SW=',SW
      write(16,*) 's12=',s12
      write(16,*) 's23=',s23
      write(16,*) 's13=',s13
      write(16,*) 'GG=',GG
      write(16,*) 'CW=',CW
      write(16,*) 'S2W=',S2W
      write(16,*) 'C2W=',C2W
      write(16,*) 'MZ=',MZ
      write(16,*) 'MW=',MW
      write(16,*) 'c12=',c12
      write(16,*) 'c23=',c23
      write(16,*) 'c13=',c13
      write(16,*) 'Me=',Me
      write(16,*) 'Mm=',Mm
      write(16,*) 'Mu=',Mu
      write(16,*) 'Md=',Md
      write(16,*) 'Mc=',Mc
      write(16,*) 'Ms=',Ms
      write(16,*) 'wt=',wt
      write(16,*) 'wZ=',wZ
      write(16,*) 'wW=',wW
      write(16,*) 'tb=',tb
      write(16,*) 'Ml=',Ml
      write(16,*) 'Mt=',Mt
      write(16,*) 'Mb=',Mb
      write(16,*) 'sb=',sb
      write(16,*) 'cb=',cb
      write(16,*) 't2b=',t2b
      write(16,*) 'c2b=',c2b
      write(16,*) 's2b=',s2b
      write(16,*) 'MG3=',MG3
      write(16,*) 'Ml1=',Ml1
      write(16,*) 'Ml2=',Ml2
      write(16,*) 'Ml3=',Ml3
      write(16,*) 'Mr1=',Mr1
      write(16,*) 'Mr2=',Mr2
      write(16,*) 'Mr3=',Mr3
      write(16,*) 'Ae=',Ae
      write(16,*) 'ls1=',ls1
      write(16,*) 'Am=',Am
      write(16,*) 'ls2=',ls2
      write(16,*) 'MSne=',MSne
      write(16,*) 'MSnm=',MSnm
      write(16,*) 'MSnl=',MSnl
      write(16,*) 'MSeLL=',MSeLL
      write(16,*) 'MSeRR=',MSeRR
      write(16,*) 'MSeLR=',MSeLR
      write(16,*) 'MSe1=',MSe1
      write(16,*) 'MSe2=',MSe2
      write(16,*) 'MSeth=',MSeth
      write(16,*) 'Zl11=',Zl11
      write(16,*) 'Zl14=',Zl14
      write(16,*) 'Zl41=',Zl41
      write(16,*) 'Zl44=',Zl44
      write(16,*) 'MSmuLL=',MSmuLL
      write(16,*) 'MSmuLR=',MSmuLR
      write(16,*) 'MSmuRR=',MSmuRR
      write(16,*) 'MSm1=',MSm1
      write(16,*) 'MSm2=',MSm2
      write(16,*) 'MSmuth=',MSmuth
      write(16,*) 'Zl22=',Zl22
      write(16,*) 'Zl25=',Zl25
      write(16,*) 'Zl52=',Zl52
      write(16,*) 'Zl55=',Zl55
      write(16,*) 'Mq1=',Mq1
      write(16,*) 'Mq2=',Mq2
      write(16,*) 'Mq3=',Mq3
      write(16,*) 'Mu1=',Mu1
      write(16,*) 'Mu2=',Mu2
      write(16,*) 'Mu3=',Mu3
      write(16,*) 'Md1=',Md1
      write(16,*) 'Md2=',Md2
      write(16,*) 'Md3=',Md3
      write(16,*) 'Zual4=',Zual4
      write(16,*) 'Zuas4=',Zuas4
      write(16,*) 'Zur0=',Zur0
      write(16,*) 'ZuQMt=',ZuQMt
      write(16,*) 'Zur1=',Zur1
      write(16,*) 'Zulam5=',Zulam5
      write(16,*) 'Zur2=',Zur2
      write(16,*) 'Zur3=',Zur3
      write(16,*) 'Zulam6=',Zulam6
      write(16,*) 'Zur4=',Zur4
      write(16,*) 'MtopS=',MtopS
      write(16,*) 'Au=',Au
      write(16,*) 'Ad=',Ad
      write(16,*) 'us1=',us1
      write(16,*) 'ds1=',ds1
      write(16,*) 'As=',As
      write(16,*) 'ds2=',ds2
      write(16,*) 'Ac=',Ac
      write(16,*) 'us2=',us2
      write(16,*) 'MSuLL=',MSuLL
      write(16,*) 'MSuLR=',MSuLR
      write(16,*) 'MSuRR=',MSuRR
      write(16,*) 'MSu1=',MSu1
      write(16,*) 'MSu2=',MSu2
      write(16,*) 'MSuth=',MSuth
      write(16,*) 'Zu11=',Zu11
      write(16,*) 'Zu14=',Zu14
      write(16,*) 'Zu41=',Zu41
      write(16,*) 'Zu44=',Zu44
      write(16,*) 'MSdLL=',MSdLL
      write(16,*) 'MSdLR=',MSdLR
      write(16,*) 'MSdRR=',MSdRR
      write(16,*) 'MSd1=',MSd1
      write(16,*) 'MSd2=',MSd2
      write(16,*) 'MSdth=',MSdth
      write(16,*) 'Zd11=',Zd11
      write(16,*) 'Zd14=',Zd14
      write(16,*) 'Zd41=',Zd41
      write(16,*) 'Zd44=',Zd44
      write(16,*) 'MScLL=',MScLL
      write(16,*) 'MScLR=',MScLR
      write(16,*) 'MScRR=',MScRR
      write(16,*) 'MSc1=',MSc1
      write(16,*) 'MSc2=',MSc2
      write(16,*) 'MScth=',MScth
      write(16,*) 'Zu22=',Zu22
      write(16,*) 'Zu25=',Zu25
      write(16,*) 'Zu52=',Zu52
      write(16,*) 'Zu55=',Zu55
      write(16,*) 'MSsLL=',MSsLL
      write(16,*) 'MSsLR=',MSsLR
      write(16,*) 'MSsRR=',MSsRR
      write(16,*) 'MSs1=',MSs1
      write(16,*) 'MSs2=',MSs2
      write(16,*) 'MSsth=',MSsth
      write(16,*) 'Zd22=',Zd22
      write(16,*) 'Zd25=',Zd25
      write(16,*) 'Zd52=',Zd52
      write(16,*) 'Zd55=',Zd55
      write(16,*) 'aMG1=',aMG1
      write(16,*) 'fiMG1=',fiMG1
      write(16,*) 'aMG2=',aMG2
      write(16,*) 'fiMG2=',fiMG2
      write(16,*) 'MG1r=',MG1r
      write(16,*) 'MG1i=',MG1i
      write(16,*) 'MG2r=',MG2r
      write(16,*) 'MG2i=',MG2i
      write(16,*) 'MG3i=',MG3i
      write(16,*) 'amu=',amu
      write(16,*) 'fimu=',fimu
      write(16,*) 'mur=',mur
      write(16,*) 'mui=',mui
      write(16,*) 'Alr=',Alr
      write(16,*) 'Ali=',Ali
      write(16,*) 'Atr=',Atr
      write(16,*) 'Ati=',Ati
      write(16,*) 'Abr=',Abr
      write(16,*) 'Abi=',Abi
      write(16,*) 'MHc=',MHc
      write(16,*) 'wHc=',wHc
      write(16,*) 'fMG1=',fMG1
      write(16,*) 'fMG2=',fMG2
      write(16,*) 'fmu=',fmu
      write(16,*) 'neutct=',neutct
      write(16,*) 'Zn11r=',Zn11r
      write(16,*) 'Zn11i=',Zn11i
      write(16,*) 'Zn12r=',Zn12r
      write(16,*) 'Zn12i=',Zn12i
      write(16,*) 'Zn13r=',Zn13r
      write(16,*) 'Zn13i=',Zn13i
      write(16,*) 'Zn14r=',Zn14r
      write(16,*) 'Zn14i=',Zn14i
      write(16,*) 'Zn21r=',Zn21r
      write(16,*) 'Zn21i=',Zn21i
      write(16,*) 'Zn22r=',Zn22r
      write(16,*) 'Zn22i=',Zn22i
      write(16,*) 'Zn23r=',Zn23r
      write(16,*) 'Zn23i=',Zn23i
      write(16,*) 'Zn24r=',Zn24r
      write(16,*) 'Zn24i=',Zn24i
      write(16,*) 'Zn31r=',Zn31r
      write(16,*) 'Zn31i=',Zn31i
      write(16,*) 'Zn32r=',Zn32r
      write(16,*) 'Zn32i=',Zn32i
      write(16,*) 'Zn33r=',Zn33r
      write(16,*) 'Zn33i=',Zn33i
      write(16,*) 'Zn34r=',Zn34r
      write(16,*) 'Zn34i=',Zn34i
      write(16,*) 'Zn41r=',Zn41r
      write(16,*) 'Zn41i=',Zn41i
      write(16,*) 'Zn42r=',Zn42r
      write(16,*) 'Zn42i=',Zn42i
      write(16,*) 'Zn43r=',Zn43r
      write(16,*) 'Zn43i=',Zn43i
      write(16,*) 'Zn44r=',Zn44r
      write(16,*) 'Zn44i=',Zn44i
      write(16,*) 'MNE1=',MNE1
      write(16,*) 'wNE1=',wNE1
      write(16,*) 'MNE2=',MNE2
      write(16,*) 'wNE2=',wNE2
      write(16,*) 'MNE3=',MNE3
      write(16,*) 'wNE3=',wNE3
      write(16,*) 'MNE4=',MNE4
      write(16,*) 'wNE4=',wNE4
      write(16,*) 'chatkn=',chatkn
      write(16,*) 'Zm11r=',Zm11r
      write(16,*) 'Zp11r=',Zp11r
      write(16,*) 'Zm11i=',Zm11i
      write(16,*) 'Zp11i=',Zp11i
      write(16,*) 'Zm12r=',Zm12r
      write(16,*) 'Zp12r=',Zp12r
      write(16,*) 'Zm12i=',Zm12i
      write(16,*) 'Zp12i=',Zp12i
      write(16,*) 'Zm21r=',Zm21r
      write(16,*) 'Zp21r=',Zp21r
      write(16,*) 'Zm21i=',Zm21i
      write(16,*) 'Zp21i=',Zp21i
      write(16,*) 'Zm22r=',Zm22r
      write(16,*) 'Zp22r=',Zp22r
      write(16,*) 'Zm22i=',Zm22i
      write(16,*) 'Zp22i=',Zp22i
      write(16,*) 'MC1=',MC1
      write(16,*) 'wC1=',wC1
      write(16,*) 'MC2=',MC2
      write(16,*) 'wC2=',wC2
      write(16,*) 'MSG=',MSG
      write(16,*) 'wSG=',wSG
      write(16,*) 'hM11=',hM11
      write(16,*) 'hM12=',hM12
      write(16,*) 'hM22=',hM22
      write(16,*) 'hM33=',hM33
      write(16,*) 'htkn=',htkn
      write(16,*) 'Mh1=',Mh1
      write(16,*) 'wh1=',wh1
      write(16,*) 'Mh2=',Mh2
      write(16,*) 'wh2=',wh2
      write(16,*) 'Mh3=',Mh3
      write(16,*) 'wh3=',wh3
      write(16,*) 'Zh11=',Zh11
      write(16,*) 'Zh12=',Zh12
      write(16,*) 'Zh13=',Zh13
      write(16,*) 'Zh21=',Zh21
      write(16,*) 'Zh22=',Zh22
      write(16,*) 'Zh23=',Zh23
      write(16,*) 'Zh31=',Zh31
      write(16,*) 'Zh32=',Zh32
      write(16,*) 'Zh33=',Zh33
      write(16,*) 'MSbLL=',MSbLL
      write(16,*) 'MSbLRr=',MSbLRr
      write(16,*) 'MSbLRi=',MSbLRi
      write(16,*) 'MSbRR=',MSbRR
      write(16,*) 'MSb1=',MSb1
      write(16,*) 'MSb2=',MSb2
      write(16,*) 'MSbth=',MSbth
      write(16,*) 'MSbfi=',MSbfi
      write(16,*) 'Zd33r=',Zd33r
      write(16,*) 'Zd36r=',Zd36r
      write(16,*) 'Zd63r=',Zd63r
      write(16,*) 'Zd66r=',Zd66r
      write(16,*) 'Zd33i=',Zd33i
      write(16,*) 'Zd36i=',Zd36i
      write(16,*) 'Zd63i=',Zd63i
      write(16,*) 'Zd66i=',Zd66i
      write(16,*) 'MStLL=',MStLL
      write(16,*) 'MStLRr=',MStLRr
      write(16,*) 'MStLRi=',MStLRi
      write(16,*) 'MStRR=',MStRR
      write(16,*) 'MSt1=',MSt1
      write(16,*) 'MSt2=',MSt2
      write(16,*) 'MStth=',MStth
      write(16,*) 'MStfi=',MStfi
      write(16,*) 'Zu33r=',Zu33r
      write(16,*) 'Zu36r=',Zu36r
      write(16,*) 'Zu63r=',Zu63r
      write(16,*) 'Zu66r=',Zu66r
      write(16,*) 'Zu33i=',Zu33i
      write(16,*) 'Zu36i=',Zu36i
      write(16,*) 'Zu63i=',Zu63i
      write(16,*) 'Zu66i=',Zu66i
      write(16,*) 'MSlLL=',MSlLL
      write(16,*) 'MSlLRr=',MSlLRr
      write(16,*) 'MSlLRi=',MSlLRi
      write(16,*) 'MSlRR=',MSlRR
      write(16,*) 'MSl1=',MSl1
      write(16,*) 'MSl2=',MSl2
      write(16,*) 'MSlth=',MSlth
      write(16,*) 'MSlfi=',MSlfi
      write(16,*) 'Zl33r=',Zl33r
      write(16,*) 'Zl36r=',Zl36r
      write(16,*) 'Zl63r=',Zl63r
      write(16,*) 'Zl66r=',Zl66r
      write(16,*) 'Zl33i=',Zl33i
      write(16,*) 'Zl36i=',Zl36i
      write(16,*) 'Zl63i=',Zl63i
      write(16,*) 'Zl66i=',Zl66i
      write(16,*) 'wSe1=',wSe1
      write(16,*) 'wSe2=',wSe2
      write(16,*) 'wSm1=',wSm1
      write(16,*) 'wSm2=',wSm2
      write(16,*) 'wSl1=',wSl1
      write(16,*) 'wSl2=',wSl2
      write(16,*) 'wSne=',wSne
      write(16,*) 'wSnm=',wSnm
      write(16,*) 'wSnl=',wSnl
      write(16,*) 'wSu1=',wSu1
      write(16,*) 'wSu2=',wSu2
      write(16,*) 'wSd1=',wSd1
      write(16,*) 'wSd2=',wSd2
      write(16,*) 'wSc1=',wSc1
      write(16,*) 'wSc2=',wSc2
      write(16,*) 'wSs1=',wSs1
      write(16,*) 'wSs2=',wSs2
      write(16,*) 'wSt1=',wSt1
      write(16,*) 'wSt2=',wSt2
      write(16,*) 'wSb1=',wSb1
      write(16,*) 'wSb2=',wSb2
      write(16,*) 'Mh1sq=',Mh1sq
      write(16,*) 'Mh2sq=',Mh2sq
      write(16,*) 'Mh3sq=',Mh3sq
      write(16,*) 'Mesq=',Mesq
      write(16,*) 'Mmsq=',Mmsq
      write(16,*) 'Mlsq=',Mlsq
      write(16,*) 'Musq=',Musq
      write(16,*) 'Mdsq=',Mdsq
      write(16,*) 'Mcsq=',Mcsq
      write(16,*) 'Mssq=',Mssq
      write(16,*) 'Mbsq=',Mbsq
      write(16,*) 'Mtsq=',Mtsq
      write(16,*) 'MSe1sq=',MSe1sq
      write(16,*) 'MSm1sq=',MSm1sq
      write(16,*) 'MSl1sq=',MSl1sq
      write(16,*) 'MSu1sq=',MSu1sq
      write(16,*) 'MSd1sq=',MSd1sq
      write(16,*) 'MSc1sq=',MSc1sq
      write(16,*) 'MSs1sq=',MSs1sq
      write(16,*) 'MSb1sq=',MSb1sq
      write(16,*) 'MSt1sq=',MSt1sq
      write(16,*) 'MSe2sq=',MSe2sq
      write(16,*) 'MSm2sq=',MSm2sq
      write(16,*) 'MSl2sq=',MSl2sq
      write(16,*) 'MSu2sq=',MSu2sq
      write(16,*) 'MSd2sq=',MSd2sq
      write(16,*) 'MSc2sq=',MSc2sq
      write(16,*) 'MSs2sq=',MSs2sq
      write(16,*) 'MSb2sq=',MSb2sq
      write(16,*) 'MSt2sq=',MSt2sq
      write(16,*) 'MSnesq=',MSnesq
      write(16,*) 'MSnmsq=',MSnmsq
      write(16,*) 'MSnlsq=',MSnlsq
#endif

      end

#include "cpvf.inc"

