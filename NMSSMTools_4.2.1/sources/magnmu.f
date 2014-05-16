*      Computation of the Muon Anomalous Moment
*
*       - Literature:
*      [1]:  F. Jegerlehner, "Essentials of the Muon g-2.",
*     Acta Phys.Polon.B38:3021,2007, hep-ph/0703125
*      [2]:  J. Bijnens and J. Prades, "The hadronic light-by-light
*      contribution to the muon anomalous magnetic moment: Where
*      do we stand?", Mod. Phys. Lett. A 22 (2007) 767,
*      arXiv:hep-ph/0702170
*      [3]:  A. Czarnecki, W. J. Marciano, A. Vainshtein," Refinements in
*     electroweak contributions to the muon anomalous magnetic moment."
*     Phys.Rev.D67:073006,2003, Erratum-ibid.D73:119901,2006, hep-ph/0212229
*      [4]:  S. P. Martin, J. D. Wells
*     "Muon anomalous magnetic dipole moment in supersymmetric theories."
*     Phys.Rev.D64:035003,2001, hep-ph/0103067
*      [5]:  J. P. Leveille, "The Second Order Weak Correction to (G-2) of the
*     Muon in Arbitrary Gauge Models."
*     Nucl.Phys.B137:63,1978
*      [6]:  J. F. Gunion, D. Hooper, B. McElrath
*     "Light neutralino dark matter in the NMSSM"
*     Phys.Rev.D73:015011,2006, hep-ph/0509024
*      [7]:  G. Degrassi, G.F. Giudice, "QED logarithms in the electroweak
*     corrections to the muon anomalous magnetic moment."
*     Phys.Rev.D58:053007,1998, hep-ph/9803384
*      [8]:  S. Heinemeyer, D. Stockinger, G. Weiglein,
*     "Electroweak and supersymmetric two-loop corrections to (g-2)(mu)."
*     Nucl.Phys.B699:103-123,2004, hep-ph/0405255
*      [9]:  K. Cheung, C.-H. Chou, O.C.W. Kong, "Muon anomalous magnetic
*     moment, two Higgs doublet model, and supersymmetry."
*     Phys.Rev.D64:111301,2001, hep-ph/0103183
*      [10]: A. Arhrib, S. Baek, "Two loop Barr-Zee type
*     contributions to (g-2)(muon) in the MSSM."
*     Phys.Rev.D65:075002,2002, hep-ph/0104225
*      [11]: D. Stockinger, "The Muon Magnetic Moment and Supersymmetry."
*     J.Phys.G34:R45-R92,2007, hep-ph/0609168


      SUBROUTINE MAGNMU(PAR,PROB)
      IMPLICIT NONE

      INTEGER I,K
      DOUBLE PRECISION PAR(*),PROB(*)

      DOUBLE PRECISION MSMU(2),USMU(2,2),UST(2,2),USB(2,2),USL(2,2),
     .   MST(2),MSB(2)
      DOUBLE PRECISION aux,Ymu,Pi,tanb,cosb,sinb,ALEM0,QNP
      DOUBLE PRECISION fc1,fc2,fn1,fn2,fhs,fhp,fhc,runmb,asf
      DOUBLE PRECISION At,Ab,Atau,lambda,Yl,fS,fPS,fSF,mtq,mbq
      DOUBLE PRECISION MSTAU(2),lambdHT(3,2),lambdHB(3,2),lambdHL(3,2)
      DOUBLE PRECISION amu2Lbos,amu2LHf,amu2LCh,amu2LSF,amu2L,amuerr,
     .   amuChar,amuNeutr,amuHiggs
      DOUBLE PRECISION amuexp,amuqed,damuqed,amuhadlo,amuhadnlo
      DOUBLE PRECISION amulbl,amuweaklo,amufermnlo,amubosnlo,damusm
      DOUBLE PRECISION errexp,errqed,errhadlo,errhadnlo
      DOUBLE PRECISION errlbl,errfermnlo,errbosnlo,errsm
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2)
      DOUBLE PRECISION PCOMP(2,2),CMASS
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION amu1L
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      
      COMMON/ALEM0/ALEM0
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/STSBSCALE/QSTSB
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin

      pi=4d0*datan(1d0)

      TANB=PAR(3)
      sinb=tanb/dsqrt(1d0+tanb**2)
      cosb=sinb/tanb

********************************************************************
*   Evaluation of the Discrepancy between Experiment and
*  SM Contributions
*  (except for 2-loop bosonic: (-2.2 +/- 0.2)d-10: see below)

*       -> Experimental Measurement BNL:
      amuexp=11659208.0d-10
      errexp=(6.3d-10)**2
      
*       -> QED Contributions (4 loops):
      amuqed=11658471.8113d-10
      errqed=(0.0162d-10)**2
      
*      -> Difference Exp.-QED:
      damuqed=amuexp-amuqed
      
*      -> LO Hadronic contr. from [1]:
      amuhadlo=692.1d-10
      errhadlo=(5.6d-10)**2
      
*      -> NLO Hadronic contr. from [1]:
      amuhadnlo=-10.03d-10
      errhadnlo=(.22d-10)**2
      
*      -> Light-by-Light contr. from [2]:
      amulbl=11.0D-10
      errlbl=(4.0D-10)**2
      
*      -> Weak SM contr., 1 Loop (error negl.)
      amuweaklo=19.482d-10
      
*      -> Weak SM fermionic 2 loop contr. from [3]
      amufermnlo=-1.512d-10
      errfermnlo=(.1d-10)**2
      
*      -> Weak SM bosonic 2 loop contr., e.g. from [8]
*       (Serves just to estimate pure SM contribution.
*      It will be deduced from the 2-loop bosonic contr. below
*      in order to avoid double counting.)
      amubosnlo=GF*MMUON**2*ALEMMZ/(72d0*dsqrt(2d0)*pi**3)
     .      *(107d0+23d0*(1d0-4d0*s2TW)**2)*dlog(MMUON/MW)
      errbosnlo=(0.2D-10)**2
      
*      -> Difference Exp.-SM:
      damusm=damuqed-amuhadlo-amuhadnlo-amulbl-amuweaklo-amufermnlo
     .    -amubosnlo

*      -> 1-sigma error of this difference:
      errsm=dsqrt(errexp+errqed+errhadlo+errhadnlo+errlbl+errfermnlo
     .     +errbosnlo)

*      -> 2-sigma lower bound to the required SUSY contr.:
      damumin=damusm-2d0*errsm

*      -> 2-sigma upper bound to the required SUSY contr.:
      damumax=damusm+2d0*errsm
      
********************************************************************

*       - SMuon masses and mixing matrix
      MSMU(1)=MSMU1
      MSMU(2)=MSMU2

      USMU(1,1)=CSMU
      USMU(1,2)=dsqrt(1d0-CSMU**2)
      USMU(2,1)=-USMU(1,2)
      USMU(2,2)=USMU(1,1)

*       - 1-Loop Chargino/SNeutrino Contribution (cf. [4]):
      Ymu=MMUON/H2Q
      aux=0d0
      do k=1,2
      aux=aux+MMUON/(12d0*MSMUNT**2)
     .     *(g2q*V(k,1)**2+(Ymu*U(k,2))**2)*Fc1((MCH(k)/MSMUNT)**2)
     .     +2d0*MCH(k)/(3d0*MSMUNT**2)
     .       *(-dsqrt(g2q)*V(k,1))*(Ymu*U(k,2))*Fc2((MCH(k)/MSMUNT)**2)
      enddo
      amuChar=MMUON/(16d0*pi**2)*aux

*       - 1-Loop Neutralino/SMuon Contribution (cf. [4]):
      aux=0d0
      do i=1,5
       do k=1,2
       aux=aux-MMUON/(12d0*MSMU(k)**2)*((1/dsqrt(2d0)
     .  *(dsqrt(g1q)*N(i,1)+dsqrt(g2q)*N(i,2))*USMU(k,1)
     .     -Ymu*N(i,4)*USMU(k,2))**2
     .    +(dsqrt(2d0*g1q)*N(i,1)*USMU(k,2)
     .     +Ymu*N(i,4)*USMU(k,1))**2)*fn1((MNEU(i)/MSMU(k))**2)
     .    +MNEU(i)/(3d0*MSMU(k)**2)*(1/dsqrt(2d0)
     .      *(dsqrt(g1q)*N(i,1)+dsqrt(g2q)*N(i,2))*USMU(k,1)
     .    -Ymu*N(i,4)*USMU(k,2))
     .      *(dsqrt(2d0*g1q)*N(i,1)*USMU(k,2)
     .    +Ymu*N(i,4)*USMU(k,1))*fn2((MNEU(i)/MSMU(k))**2)
       enddo
      enddo
      amuNeutr=MMUON/(16d0*pi**2)*aux

*       - 1-Loop Higgs/Muon Contribution (cf. [5,6]):
      aux=0d0
      do i=1,3
      aux=aux+SCOMP(i,2)**2*fhs(SMASS(i)/MMUON)
      enddo
      do i=1,2
      aux=aux-PCOMP(i,1)**2*sinb**2*fhp(PMASS(i)/MMUON)
      enddo
      aux=aux+sinb**2*fhc(CMASS/MMUON)
      amuHiggs=(Ymu/(4d0*pi))**2*aux
      
      amu1L=amuChar+amuNeutr+amuHiggs
      
*       - Large Logarithms: 1-loop-like 2-loop contributions (cf. [7]):
      QNP=Max(MCH(2),MSMU(2))
      amu1L=amu1L*(1d0-4d0*ALEM0/pi*dlog(QNP/MMUON))
      
*       - 2-loop bosonic contributions (Leading Log) (cf. [8]):
      aux=0d0
      do i=1,3
      aux=aux+SCOMP(i,2)*(SCOMP(i,2)*cosb-SCOMP(i,1)*sinb)/SMASS(i)**2
      enddo
      aux=(cosb**2-sinb**2)*MZ**2/cosb*aux
      aux=(98d0+9d0*aux+23d0*(1d0-4d0*s2TW)**2)/30d0
      amu2Lbos=5d0/(24d0*dsqrt(2d0)*pi**3)*GF*MMUON**2*ALEMMZ
     .   *(aux*dlog(MMUON**2/MW**2))
      
*       - 2-loop Higgs/Fermion contribution (cf. [9]):
      mtq=MT/(1d0+4d0/(3d0*pi)*asf(MT)+11d0/pi**2*asf(MT)**2)
      mbq=runmb(MB)
      aux=0d0
      do i=1,3
      aux=aux+4d0/3d0*SCOMP(i,1)*SCOMP(i,2)/sinb
     .                  *FS((MTQ/SMASS(i))**2)
     .       +1d0/3d0*SCOMP(i,2)**2/cosb*FS((MBQ/SMASS(i))**2)
     .       +SCOMP(i,2)**2/cosb*FS((mtau/SMASS(i))**2)
      enddo
      aux=aux/cosb
      do i=1,2
      aux=aux+PCOMP(i,1)**2*(4d0/3d0*FPS((MTQ/PMASS(i))**2)
     .      +1d0/3d0*FPS((MBQ/PMASS(i))**2)*tanb**2
     .      +FPS((mtau/PMASS(i))**2)*tanb**2)
      enddo
      amu2LHf=GF*MMUON**2*ALEM0/(4d0*dsqrt(2d0)*pi**3)*aux
      
*       - 2-loop Sfermion contribution (cf. [10]):
      AT=PAR(12)
      AB=PAR(13)
      Atau=PAR(14)
      lambda=PAR(1)
      Yl=MTAU/H2Q
      
      MST(1)=MST1
      MST(2)=MST2
      UST(1,1)=CST
      UST(1,2)=dsqrt(1d0-CST**2)
      UST(2,2)=UST(1,1)
      UST(2,1)=-UST(1,2)
      
      MSB(1)=MSB1
      MSB(2)=MSB2
      USB(1,1)=CSB
      USB(1,2)=dsqrt(1d0-CSB**2)
      USB(2,2)=USB(1,1)
      USB(2,1)=-USB(1,2)
      
      MSTAU(1)=MSL1
      MSTAU(2)=MSL2
      USL(1,1)=CSL
      USL(1,2)=dsqrt(1d0-CSL**2)
      USL(2,2)=USL(1,1)
      USL(2,1)=-USL(1,2)
      
      aux=0d0
      do i=1,3
      do k=1,2
      lambdHT(i,k)=HTQ*(At*SCOMP(i,1)-muq*SCOMP(i,2)
     .       -lambda*H2Q*SCOMP(i,3))*UST(k,2)*UST(k,1)
     .       +(HTQ**2*H1Q*SCOMP(i,1)-g1q/3d0
     .  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*UST(k,2)**2
     .       +(HTQ**2*H1Q*SCOMP(i,1)-(3d0*g2q-g1q)/12d0
     .  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*UST(k,1)**2
      lambdHB(i,k)=HBQ*(Ab*SCOMP(i,2)-muq*SCOMP(i,1)
     .       -lambda*H1Q*SCOMP(i,3))*USB(k,2)*USB(k,1)
     .       +(HBQ**2*H2Q*SCOMP(i,2)+g1q/6d0
     .  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USB(k,2)**2
     .       +(HBQ**2*H2Q*SCOMP(i,2)+(3d0*g2q+g1q)/12d0
     .  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USB(k,1)**2
      lambdHL(i,k)=Yl*(Atau*SCOMP(i,2)-muq*SCOMP(i,1)
     .       -lambda*H1Q*SCOMP(i,3))*USL(k,1)*USL(k,2)
     .       +(Yl**2*H2Q*SCOMP(i,2)+g1q/2d0
     .  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USL(k,2)**2
     .       +(Yl**2*H2Q*SCOMP(i,2)+(g2q-g1q)/4d0
     .  *(H1Q*SCOMP(i,1)-H2Q*SCOMP(i,2)))*USL(k,1)**2
      
      lambdHT(i,k)=lambdHT(i,k)*dsqrt(2d0)*2d0*MW/dsqrt(g2)
      lambdHB(i,k)=lambdHB(i,k)*dsqrt(2d0)*2d0*MW/dsqrt(g2)
      lambdHL(i,k)=lambdHL(i,k)*dsqrt(2d0)*2d0*MW/dsqrt(g2)
      
      aux=aux+SCOMP(i,2)*(4d0/3d0*lambdHT(i,k)
     .      *FSF((MST(k)/SMASS(i))**2)/MST(k)**2
     .     +1d0/3d0*lambdHB(i,k)
     .      *FSF((MSB(k)/SMASS(i))**2)/MSB(k)**2
     .     +lambdHL(i,k)
     .      *FSF((MSTAU(k)/SMASS(i))**2)/MSTAU(k)**2)
      enddo
      enddo
      amu2LSF=GF*MMUON**2*ALEMMZ/(4d0*dsqrt(2d0)*pi**3*cosb)*aux
      
*       - 2-loop Chargino contribution (cf. [11]):
      aux=0d0
      do k=1,2
      do i=1,3
      aux=aux+SCOMP(i,2)/cosb*dsqrt(2d0)*MW/MCH(k)
     .  *((U(k,1)*V(k,2)*SCOMP(i,1)+U(k,2)*V(k,1)*SCOMP(i,2))
     .    +lambda/dsqrt(g2)*U(k,2)*V(k,2)*SCOMP(i,3))
     .  *fS((MCH(k)/SMASS(i))**2)
      enddo
      do i=1,2
      aux=aux-PCOMP(i,1)*tanb*dsqrt(2d0)*MW/MCH(k)
     .  *((U(k,1)*V(k,2)*cosb+U(k,2)*V(k,1)*sinb)*PCOMP(i,1)
     .    -lambda/dsqrt(g2)*PCOMP(i,2)*U(k,2)*V(k,2))
     .  *fPS((MCH(k)/PMASS(i))**2)
      enddo
      enddo
      amu2LCh=GF*MMUON**2*dsqrt(2d0)*ALEMMZ/(8d0*pi**3)*aux
      
*       - Conclusion and errors  (cf. [11])
*      -> Pure SUSY contribution to the discrepancy (G-2)mu:
*       (Here the bosonic SM 2-loop contr. is deduced)
      amu2L=amu1L+amu2Lbos+amu2LHf+amu2LSF+amu2LCh-amubosnlo
*      -> Theoretical errors (added linearly):
      amuerr=2.8d-10+.02d0*dabs(amu1L)+.3d0*dabs(amu2L-amu1L)
*       -> Theoretical Central Value:
      delmagmu=amu2L
*       -> Theoretical Error Bounds:
      amuthmax=amu2L+amuerr
      amuthmin=amu2L-amuerr
      
*       - Comparison with the required SUSY contr.:
      IF(damumin-amuthmax.ge.0d0)THEN
      PROB(37)=1d0-amuthmax/damumin
      ENDIF
      IF(damuMax-amuthmin.le.0d0)THEN
      PROB(37)=1d0-amuthmin/damuMax
      ENDIF

      return
      END
      
      
**********************************************************************      
**********************************************************************
**********************************************************************

      double precision function fn1(x)

      implicit none
      double precision x
      if(dabs(x).gt.1.d-3)then
      if(dabs(x-1d0).gt.1.d-3)then
      fn1=2d0/(1d0-x)**4
     .       *(1d0-6d0*x+3d0*x**2+2d0*x**3-6d0*x**2*dlog(x))
      else
      fn1=1d0
      endif
      else
      fn1=2d0
      endif
      return
      end

*********************************************************************

      double precision function fn2(x)

      implicit none
      double precision x
      if(dabs(x).gt.1.d-3)then
      if(dabs(x-1d0).gt.1.d-3)then
      fn2=3d0/(1d0-x)**3*(1d0-x**2+2d0*x*dlog(x))
      else
      fn2=1d0
      endif
      else
      fn2=3d0
      endif
      return
      end

*********************************************************************

      double precision function fc1(x)

      implicit none
      double precision x
      if(dabs(x).gt.1.d-3)then
      if(dabs(x-1d0).gt.1.d-3)then
      fc1=2d0/(1d0-x)**4
     .       *(2d0+3d0*x-6d0*x**2+x**3+6d0*x*dlog(x))
      else
      fc1=1d0
      endif
      else
      fc1=4d0
      endif
      return
      end

*********************************************************************

      double precision function fc2(x)

      implicit none
      double precision x
      if(dabs(x-1d0).gt.1.d-3)then
      fc2=3d0/(1d0-x)**3*(-3d0+4d0*x-x**2-2d0*dlog(x))/2d0
      else
      fc2=1d0
      endif
      return
      end

*********************************************************************

      double precision function fhs(x)

      implicit none
      double precision x,aux
      if(abs(x).gt.2.001d0)then
       if(abs(x).gt.10d0)then
       aux=(3d0+x**2)*dlog(x**2)/x**4-7d0/(6d0*x**2)
       else
       aux=3d0/2d0-x**2+x**2/2d0*(x**2-3d0)*dlog(x**2)
     .       +(5d0*x**4-x**6-4d0*x**2)/(2d0*dsqrt(x**2*(x**2-4d0)))
     .       *dlog((2d0-x**2-dsqrt(x**2*(x**2-4)))
     .      *(x**2-dsqrt(x**2*(x**2-4d0)))
     .      /((2d0-x**2+dsqrt(x**2*(x**2-4d0)))
     .      *(x**2+dsqrt(x**2*(x**2-4d0)))))
       endif
      else
      aux=-5d0/2d0+4d0*dlog(2d0)
      endif
      fhs=aux
      return
      end

*********************************************************************

      double precision function fhp(x)

      implicit none
      double precision x,aux
      if(abs(x).gt.2.001d0)then
       if(abs(x).gt.15d0)then
       aux=(5d0+x**2)*dlog(x**2)/x**4-11d0/(6d0*x**2)
       else
       aux=1d0/2d0+x**2+x**2/2d0*(1d0-x**2)*dlog(x**2)
     .       +x**4*(x**2-3d0)/(2d0*dsqrt(x**2*(x**2-4d0)))
     .       *dlog((2d0-x**2-dsqrt(x**2*(x**2-4)))
     .      *(x**2-dsqrt(x**2*(x**2-4d0)))
     .      /((2d0-x**2+dsqrt(x**2*(x**2-4d0)))
     .      *(x**2+dsqrt(x**2*(x**2-4d0)))))
       endif
      else
      aux=17d0/2d0-12d0*dlog(2d0)
      endif
      fhp=aux
      return
      end

*********************************************************************

      double precision function fhc(x)

      implicit none
      double precision x,aux
      if(abs(x).gt.2.001d0)then
       if(abs(x).gt.10d0)then
       aux=-(1d0+2d0*x**2)/(12d0*x**4)
       else
       aux=1d0/2d0-x**2+x**2*(x**2-1d0)*dlog(x**2/(x**2-1d0))
       endif
      else
      aux=-7d0/2d0+24d0*dlog(2d0)-12d0*dlog(3d0)
      endif
      fhc=aux
      return
      end

*********************************************************************

      double precision function fPS(x)

      implicit none
      double precision x,aux
      IF(x.lt.0.002d0)THEN
      aux=0.4995501108d0*dsqrt(x)+1133.437039d0*x
     .      -371.1353082d0*x*dsqrt(x)-1087.028512d0*dlog(1d0+x)
      ELSE
       IF(x.lt.0.02d0)THEN
       aux=-1.018546244d-3+7.096260514d-1*dsqrt(x)
     .     +38.60712022d0*x-331.8418712d0*x*dsqrt(x)
     .     +1522.657713d0*x**2-2874.048399d0*x**2*dsqrt(x)
       ELSE
        IF(x.lt.0.1d0)THEN
        aux=4.306498122d-2+32.56225304d0*dlog(1d0+x)
     .      +925.7514262d0*dlog(1d0+x)**2
     .      -148.1323591d0*dsqrt(x)*dlog(1d0+x)
     .      -648.4904538d0*x*dlog(1d0+x)
        ELSE
         IF(x.lt.2.5d0)THEN
         aux=-1.365013496d-1-1.858455623d0*dlog(1d0+x)
     .       -5.996763746d-1*dlog(1d0+x)**2
     .       +4.390843985d-1*dsqrt(x)*dlog(1d0+x)
     .       -1.444359743d-1*x*dlog(1d0+x)
     .       +3.852425143d0*dsqrt(x)
         ELSE
          IF(x.lt.100d0)THEN
          aux=4.304425955d-1+6.766323794d-2*dlog(1d0+x)
     .  -1.584446296d-1*dlog(1d0+x)**2
     .  -2.787080541d-1*dsqrt(x)*dlog(1d0+x)
     .  +1.557845370d-3*x*dlog(1d0+x)
     .  +2.139180566d0*dsqrt(x)
          ELSE
           IF(x.lt.10000d0)THEN
           aux=2.025445594d0+9.960866255d-1*dlog(x)
     .   +1.122896720d-4*dsqrt(x)
           ELSE
            aux=2.000835136d0+9.9992369d-1*dlog(x)
     .    +2.327105016d-7*dsqrt(x)
           ENDIF
          ENDIF
         ENDIF
        ENDIF
       ENDIF
      ENDIF
      fPS=aux
      return
      end

*********************************************************************

      double precision function fS(x)

      implicit none
      double precision x,fPS
      fS=(2d0*x-1d0)*fPS(x)-2d0*x*(2d0+dlog(x))
      return
      end
            
*********************************************************************

      double precision function fSF(x)

      implicit none
      double precision x,fPS
      fSF=x/2d0*(2d0+dlog(x)-fPS(x))
      return
      end
