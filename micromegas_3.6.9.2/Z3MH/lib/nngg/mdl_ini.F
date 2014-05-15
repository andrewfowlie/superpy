*     LanHEP output produced at Tue Apr 23 15:14:21 2013
*     Model named 'Z3 Inert Dublet'

#include "DgetNEWv0.inc"
#include "CgetNEWv0.inc"

      subroutine ModelDefaults
      implicit none
#include "model.h"

      EE = 0.31223D0
      GG = 1.117D0
      SW = 0.481D0
      MZ = 91.187D0
      Q = 100D0
      wZ = 2.502D0
      wW = 2.094D0
      Mm = 0.1057D0
      Mt = 1.777D0
      Mu = 0.01D0
      Md = 0.01D0
      Mc = 1.3D0
      Ms = 0.2D0
      Mtop = 170D0
      wtop = 1.442D0
      Mh = 200D0
      la3 = 0.2D0
      la2 = 0.01D0
      la4 = 0.01D0
      laS = 0.1D0
      laS1 = 0.1D0
      laS2 = 0.1D0
      laS21 = 0.1D0
      Mdm1 = 50D0
      Mdm2 = 200D0
      sinDm = 0.1D0
      muS = 100D0

      end


      subroutine ModelConstIni
      implicit none
      integer fail

#include "model.h"
      real*8 zero/0.d0/ 
      double precision MbEff
            
      fail=0
      call ModelDefaults
      EE   = slhaValFormat('Z3_Parameters',zero,'EE    %lf')
      SW   = slhaValFormat('Z3_Parameters',zero,'SW    %lf')
      MZ   = slhaValFormat('Z3_Parameters',zero,'MZ    %lf')
      Q    = slhaValFormat('Z3_Parameters',zero,'Q     %lf')
      wZ   = slhaValFormat('Z3_Parameters',zero,'wZ    %lf')
      wW   = slhaValFormat('Z3_Parameters',zero,'wW    %lf')
      Mm   = slhaValFormat('Z3_Parameters',zero,'Mm    %lf')
      Mt   = slhaValFormat('Z3_Parameters',zero,'Ml    %lf')
      Mu   = slhaValFormat('Z3_Parameters',zero,'Mu    %lf')
      Md   = slhaValFormat('Z3_Parameters',zero,'Md    %lf')
      Mc   = slhaValFormat('Z3_Parameters',zero,'Mc    %lf')
      Ms   = slhaValFormat('Z3_Parameters',zero,'Ms    %lf')
      Mtop = slhaValFormat('Z3_Parameters',zero,'Mtop  %lf')
      wtop = slhaValFormat('Z3_Parameters',zero,'wtop  %lf')
      Mh   = slhaValFormat('Z3_Parameters',zero,'Mh    %lf')
      la3  = slhaValFormat('Z3_Parameters',zero,'la3   %lf')
      la2  = slhaValFormat('Z3_Parameters',zero,'la2   %lf')
      la4  = slhaValFormat('Z3_Parameters',zero,'la4   %lf')
      laS  = slhaValFormat('Z3_Parameters',zero,'laS   %lf')
      laS1 = slhaValFormat('Z3_Parameters',zero,'laS1  %lf')
      laS2 = slhaValFormat('Z3_Parameters',zero,'laS2  %lf')
      laS21= slhaValFormat('Z3_Parameters',zero,'laS21 %lf')
      Mdm1 = slhaValFormat('Z3_Parameters',zero,'Mdm1  %lf')
      Mdm2 = slhaValFormat('Z3_Parameters',zero,'Mdm2  %lf')
      sinDm= slhaValFormat('Z3_Parameters',zero,'sinDm %lf')
      muS  = slhaValFormat('Z3_Parameters',zero,'muppS %lf')

#ifdef TEST

        write(*,*) ' EE   ',EE
        write(*,*) ' GG   ',GG
        write(*,*) ' SW   ',SW
        write(*,*) ' MZ   ',MZ
        write(*,*) ' Q    ',Q
        write(*,*) ' wZ   ',wZ
        write(*,*) ' wW   ',wW
        write(*,*) ' Mm   ',Mm
        write(*,*) ' Mt   ',Mt
        write(*,*) ' Mu   ',Mu
        write(*,*) ' Md   ',Md
        write(*,*) ' Mc   ',Mc
        write(*,*) ' Ms   ',Ms
        write(*,*) ' Mtop ',Mtop
        write(*,*) ' wtop ',wtop
        write(*,*) ' Mh   ',Mh
        write(*,*) ' la3  ',la3
        write(*,*) ' la2  ',la2
        write(*,*) ' la4  ',la4
        write(*,*) ' laS  ',laS
        write(*,*) ' laS1 ',laS1
        write(*,*) ' laS2 ',laS2
        write(*,*) ' laS21',laS21
        write(*,*) ' Mdm1 ',Mdm1
        write(*,*) ' Mdm2 ',Mdm2
        write(*,*) ' sinDm',sinDm
        write(*,*) ' muS  ',muS

#endif

      CW = sqrt(1D0-SW**2)
      MW = MZ*CW
      Mb = MbEff(Q)
      V = 2D0*MW/EE*SW
      cosDm = sqrt(1D0-sinDm**2)
      muSq = Mdm2**2*cosDm**2+Mdm1**2*sinDm**2-V**2*laS1/2D0
      mush = -4D0*(Mdm2**2-Mdm1**2)*sinDm*cosDm/V/Sqrt2
      mu2q = -V**2*(la4+la3)/2D0+cosDm**2*Mdm1**2+sinDm**2
     &      *Mdm2**2
      MHC = sqrt(mu2q+V**2*la3/2D0)

      call aaini01
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
      AAABR(7) = EE/CW*MW*(1D0/SW-2D0*SW)
      AAABR(8) = EE/CW*MW/SW
      AAABR(9) = EE/CW**2*MW/SW
      AAABR(10) = CW*EE/SW
      AAABR(11) = EE/MW*Mh**2/SW
      AAABR(12) = Sqrt2*mush*sinDm+4D0*MW/EE*SW*cosDm*la4
      AAABR(13) = Sqrt2*cosDm*mush-4D0*MW/EE*SW*la4*sinDm
      AAABR(14) = Sqrt2*mush*(cosDm**2+sinDm**2)
      AAABR(15) = Sqrt2*cosDm*mush*sinDm+4D0*MW/EE*SW*cosDm
     &      **2*la3+4D0*MW/EE*SW*cosDm**2*la4+4D0*MW/EE*SW
     &      *laS1*sinDm**2
      AAABR(16) = Sqrt2*cosDm**2*mush-Sqrt2*mush*sinDm**2
     &      -8D0*MW/EE*SW*cosDm*la3*sinDm-8D0*MW/EE*SW*cosDm
     &      *la4*sinDm+8D0*MW/EE*SW*cosDm*laS1*sinDm
      AAABR(17) = Sqrt2*cosDm*mush*sinDm-4D0*MW/EE*SW*cosDm
     &      **2*laS1-4D0*MW/EE*SW*la3*sinDm**2-4D0*MW/EE
     &      *SW*la4*sinDm**2
      AAABR(18) = MW/EE*SW*la3
      AAABR(19) = Sqrt2*muS*sinDm**3+2D0*MW/EE*SW*cosDm*laS21
     &      *sinDm**2
      AAABR(20) = 3D0*Sqrt2*cosDm*muS*sinDm**2-2D0*MW/EE
     &      *SW*laS21*sinDm**3+4D0*MW/EE*SW*cosDm**2*laS21
     &      *sinDm
      AAABR(21) = 3D0*Sqrt2*cosDm**2*muS*sinDm+2D0*MW/EE
     &      *SW*cosDm**3*laS21-4D0*MW/EE*SW*cosDm*laS21*sinDm
     &      **2
      AAABR(22) = Sqrt2*cosDm**3*muS-2D0*MW/EE*SW*cosDm**2
     &      *laS21*sinDm
      AAABR(23) = EE*(SW/CW-CW/SW)
      AAABR(24) = EE/SW
      AAABR(25) = EE*(SW/CW+CW/SW)
      AAABR(26) = EE/SW*cosDm
      AAABR(27) = EE*cosDm**2*(SW/CW+CW/SW)
      AAABR(28) = EE*cosDm*sinDm*(SW/CW+CW/SW)
      AAABR(29) = EE/SW*sinDm
      AAABR(30) = EE*sinDm**2*(SW/CW+CW/SW)
      AAABR(31) = EE/CW*MW*SW
      AAABR(32) = EE*MW*(2D0*SW+SW**3/CW**2+CW**2/SW)
      AAABR(33) = EE/MW*Mb/SW
      AAABR(34) = EE/MW*Mb/SW*Sqrt2
      AAABR(35) = EE/MW*Mtop/SW*Sqrt2
      AAABR(36) = EE/MW*Mc/SW
      AAABR(37) = EE/MW*Ms/SW*Sqrt2
      AAABR(38) = EE/MW*Mc/SW*Sqrt2
      AAABR(39) = EE/MW*Md/SW
      AAABR(40) = EE/MW*Md/SW*Sqrt2
      AAABR(41) = EE/MW*Mu/SW*Sqrt2
      AAABR(42) = EE/MW*Mt/SW
      AAABR(43) = EE/MW*Mt/SW*Sqrt2
      AAABR(44) = EE/MW*Mm/SW
      AAABR(45) = EE/MW*Mm/SW*Sqrt2
      AAABR(46) = EE/MW*Ms/SW
      AAABR(47) = EE/MW*Mtop/SW
      AAABR(48) = EE/MW*Mu/SW
      AAABR(49) = EE*(SW/CW+3D0*CW/SW)
      AAABR(50) = EE/CW*SW
      AAABR(51) = EE/SW*Sqrt2
      AAABR(52) = EE*(SW/CW-3D0*CW/SW)
      AAABR(53) = EE**2/MW**2*Mh**2/SW**2
      AAABR(54) = cosDm**2*la3+laS1*sinDm**2
      AAABR(55) = cosDm*sinDm*(la3-laS1)
      AAABR(56) = cosDm**2*laS1+la3*sinDm**2
      AAABR(57) = la3+la4
      AAABR(58) = cosDm*la4
      AAABR(59) = la4*sinDm
      AAABR(60) = laS21*sinDm**2
      AAABR(61) = cosDm*laS21*sinDm
      AAABR(62) = cosDm**2*laS21
      AAABR(63) = cosDm**2*la3+cosDm**2*la4+laS1*sinDm**2
      AAABR(64) = cosDm*sinDm*(la3+la4-laS1)
      AAABR(65) = cosDm**2*laS1+la3*sinDm**2+la4*sinDm**2
      AAABR(66) = cosDm*laS21*sinDm**2
      AAABR(67) = laS21*(sinDm**3-2D0*cosDm**2*sinDm)
      AAABR(68) = laS21*(cosDm**3-2D0*cosDm*sinDm**2)
      AAABR(69) = cosDm**2*laS21*sinDm
      AAABR(70) = cosDm**4*la2+laS*sinDm**4+cosDm**2*laS2
     &      *sinDm**2
      AAABR(71) = 2D0*cosDm*laS*sinDm**3-cosDm*laS2*sinDm
     &      **3-2D0*cosDm**3*la2*sinDm+cosDm**3*laS2*sinDm
      AAABR(72) = cosDm**4*laS2+laS2*sinDm**4+6D0*cosDm**2
     &      *la2*sinDm**2+6D0*cosDm**2*laS*sinDm**2-4D0*cosDm
     &      **2*laS2*sinDm**2
      AAABR(73) = 2D0*cosDm**2*la2+laS2*sinDm**2
      AAABR(74) = cosDm**4*laS2+laS2*sinDm**4+2D0*cosDm**2
     &      *la2*sinDm**2+2D0*cosDm**2*laS*sinDm**2
      AAABR(75) = 2D0*cosDm*la2*sinDm**3-cosDm*laS2*sinDm
     &      **3-2D0*cosDm**3*laS*sinDm+cosDm**3*laS2*sinDm
      AAABR(76) = cosDm*sinDm*(2D0*la2-laS2)
      AAABR(77) = cosDm**2*sinDm**2*(la2+laS-laS2)
      AAABR(78) = cosDm**4*laS+la2*sinDm**4+cosDm**2*laS2
     &      *sinDm**2
      AAABR(79) = cosDm**2*laS2+2D0*la2*sinDm**2
      AAABR(80) = EE**2*(SW/CW-CW/SW)
      AAABR(81) = EE**2
      AAABR(82) = EE**2*(2D0-SW**2/CW**2-CW**2/SW**2)
      AAABR(83) = EE**2/SW**2
      AAABR(84) = EE**2/CW
      AAABR(85) = EE**2/SW
      AAABR(86) = EE**2*(2D0+SW**2/CW**2+CW**2/SW**2)
      AAABR(87) = EE**2*cosDm**2*(2D0+SW**2/CW**2+CW**2/SW
     &      **2)
      AAABR(88) = EE**2/SW**2*cosDm**2
      AAABR(89) = EE**2*cosDm*sinDm*(2D0+SW**2/CW**2+CW**2
     &      /SW**2)
      AAABR(90) = EE**2/SW**2*cosDm*sinDm
      AAABR(91) = EE**2/CW*cosDm
      AAABR(92) = EE**2/SW*cosDm
      AAABR(93) = EE**2*sinDm**2*(2D0+SW**2/CW**2+CW**2/SW
     &      **2)
      AAABR(94) = EE**2/SW**2*sinDm**2
      AAABR(95) = EE**2/CW*sinDm
      AAABR(96) = EE**2/SW*sinDm
      AAABR(97) = CW*EE**2/SW
      AAABR(98) = CW**2*EE**2/SW**2
      AAABR(99) = GG**2
      end

      subroutine mtrini
      implicit none
#include "model.h"

      integer m1,m2,m3,m4


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
      fail=0
      end

************************************************

      subroutine ModelDigest
      implicit none

#include "model.h"

      write(16,*) 'EE=',EE
      write(16,*) 'GG=',GG
      write(16,*) 'SW=',SW
      write(16,*) 'CW=',CW
      write(16,*) 'MZ=',MZ
      write(16,*) 'MW=',MW
      write(16,*) 'Q=',Q
      write(16,*) 'Mb=',Mb
      write(16,*) 'V=',V
      write(16,*) 'wZ=',wZ
      write(16,*) 'wW=',wW
      write(16,*) 'Mm=',Mm
      write(16,*) 'Mt=',Mt
      write(16,*) 'Mu=',Mu
      write(16,*) 'Md=',Md
      write(16,*) 'Mc=',Mc
      write(16,*) 'Ms=',Ms
      write(16,*) 'Mtop=',Mtop
      write(16,*) 'wtop=',wtop
      write(16,*) 'Mh=',Mh
      write(16,*) 'la3=',la3
      write(16,*) 'la2=',la2
      write(16,*) 'la4=',la4
      write(16,*) 'laS=',laS
      write(16,*) 'laS1=',laS1
      write(16,*) 'laS2=',laS2
      write(16,*) 'laS21=',laS21
      write(16,*) 'Mdm1=',Mdm1
      write(16,*) 'Mdm2=',Mdm2
      write(16,*) 'sinDm=',sinDm
      write(16,*) 'muS=',muS
      write(16,*) 'cosDm=',cosDm
      write(16,*) 'muSq=',muSq
      write(16,*) 'mush=',mush
      write(16,*) 'mu2q=',mu2q
      write(16,*) 'MHC=',MHC

      end

