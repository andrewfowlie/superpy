      SUBROUTINE RGESGM(PAR,IFAIL)

*   Subroutine to integrate the RG equations for the gauge and Yukawa
*   couplings up to the messenger scale MMESS, through a CALL of the
*   subroutine ODEINTGM that is part of the file integgm.f
*
*   Below Q2 all sparticle/heavy Higgs thresholds are taken into
*   account in the naive step function approximation.
*   Above the Susy scale Q2 the two loop beta functions are used.
*   (Note: The sparticle thresholds are consistent even if a
*   sparticle mass is above Q2: then the threshold effect between
*   MT and Q2 "anticipates" the threshold effect between Q2 and MMESS)
*
***********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,NN
      PARAMETER (NN=8)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN)
      DOUBLE PRECISION HBOT,HTOP,HTAU,PI,COEF
      DOUBLE PRECISION TANB,h1,h2,sb2,cb2
      DOUBLE PRECISION g1z,g2z,g3z,g1t,g2t,g3t
      DOUBLE PRECISION MA2,Q2,DELMB,RUNMB,QSTSB
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,SW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION LS,KS,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     .      HBOTMES,HTAUMES
      DOUBLE PRECISION MUS,NUS,SIGNKAPPA
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION ALSMT,ALSMA,ALSQ,DLA,DLQA,F1,F2,HTMA,LMAMT

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     .      HBOTMES,HTAUMES      
      COMMON/DELMB/DELMB
      COMMON/STSBSCALE/QSTSB      
      COMMON/MESCAL/MSUSYEFF,MMESS,N5

      EXTERNAL DERIVS,RKQSGM

      EPS=1.D-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)
      SW=.2312d0
      TANB=PAR(3)
      cb2=1d0/(1d0+tanb**2)
      sb2=tanb**2*cb2

* Parameters at Q2=M_SUSY**2
      LS=PAR(1)
      KS=PAR(2)
      IF(KS.NE.0d0)THEN
       SIGNKAPPA=DABS(KS)/KS
      ELSE
       SIGNKAPPA=1d0
      ENDIF
      MUS=PAR(4)
      NUS=MUS*KS/LS
      MA2=PAR(23)**2

* At MZ
      h2=1d0/DSQRT(2d0*DSQRT(2d0)*(1d0+TANB**2)*GF)
      h1=h2*TANB

* Running strong coupling
* Input: ALSMZ = Alphas_s(MZ,MS_bar)

* First: g_3**2 at MZ in the DR_bar scheme:
      g3z=4d0*PI*ALSMZ/(1d0-ALSMZ/(4d0*PI))

* Next: g_3**2 at M_top (for the running Yukawas):
      g3t=g3z/(1d0+g3z*COEF*DLOG((MT/MZ)**2)*23d0/3d0)

* Finally: g_3**2 at Q2=M_SUSY**2
* including the top and sparticle thresholds:
      g3s=g3z/(1d0+g3z*COEF*(DLOG(Q2/MZ**2)*23d0/3d0
     .    -DLOG(Q2/MT**2)*2d0/3d0
     .    -DLOG(Q2/MAX(PAR(15),MZ**2))*2d0/3d0
     .    -DLOG(Q2/MAX(PAR(16),MZ**2))/3d0
     .    -DLOG(Q2/MAX(PAR(17),MZ**2))/3d0
     .    -DLOG(Q2/MAX(PAR(7),MZ**2))/3d0
     .    -DLOG(Q2/MAX(PAR(8),MZ**2))/6d0
     .    -DLOG(Q2/MAX(PAR(9),MZ**2))/6d0
     .    -DLOG(Q2/MAX(PAR(22)**2,MZ**2))*2d0))

* Running SU(2) coupling
* Use: the above value for SW=sin^2_theta and
* ALEMMZ=alpha_em(MZ,MS_bar)

* First: g_2**2 at MZ in the DR_bar scheme:
      g2z=4d0*PI*ALEMMZ/(SW-ALEMMZ/(6d0*PI))

* Next: g_2**2 at M_top (for the running Yukawas):
      g2t=g2z/(1d0+g2z*COEF*DLOG((MT/MZ)**2)*19d0/6d0)

* Finally: g_2**2 at Q2=M_SUSY**2
* including the Higgs and sparticle thresholds:
      g2s=g2z/(1d0+g2z*COEF*(DLOG(Q2/MZ**2)*19d0/6d0
     .    -DLOG(Q2/MAX(MA2,MZ**2))/6d0
     .    -DLOG(Q2/MAX(MUS**2,MZ**2))*2d0/3d0
     .    -DLOG(Q2/MAX(PAR(7),MZ**2))/2d0
     .    -DLOG(Q2/MAX(PAR(10),MZ**2))/6d0
     .    -DLOG(Q2/MAX(PAR(15),MZ**2))
     .    -DLOG(Q2/MAX(PAR(18),MZ**2))/3d0
     .    -DLOG(Q2/MAX(PAR(21)**2,MZ**2))*4d0/3d0))

* Running U(1) coupling
* Use: the above value for SW=sin^2_theta and
* ALEMMZ=alpha_em(MZ,MS_bar)

* First: g_1**2 at MZ in the DR_bar=MS_bar scheme:
      g1z=4d0*PI*ALEMMZ/(1d0-SW)

* Next: g_1**2 at M_top (for the running Yukawas):
      g1t=g1z/(1d0+g1z*COEF*DLOG((MT/MZ)**2)*53d0/9d0)

* Finally: g_1**2 at Q2=M_SUSY**2
*        including the top, Higgs and sparticle thresholds:
      g1s=g1z/(1d0-g1z*COEF*(DLOG(Q2/MZ**2)*53d0/9d0
     .    +DLOG(Q2/MT**2)*17d0/18d0
     .    +DLOG(Q2/MAX(MA2,MZ**2))/6d0
     .    +DLOG(Q2/MAX(MUS**2,MZ**2))*2d0/3d0
     .    +DLOG(Q2/MAX(PAR(7),MZ**2))/18d0
     .    +DLOG(Q2/MAX(PAR(8),MZ**2))*4d0/9d0
     .    +DLOG(Q2/MAX(PAR(9),MZ**2))/9d0
     .    +DLOG(Q2/MAX(PAR(10),MZ**2))/6d0
     .    +DLOG(Q2/MAX(PAR(11),MZ**2))/3d0
     .    +DLOG(Q2/MAX(PAR(15),MZ**2))/9d0
     .    +DLOG(Q2/MAX(PAR(16),MZ**2))*8d0/9d0
     .    +DLOG(Q2/MAX(PAR(17),MZ**2))*2d0/9d0
     .    +DLOG(Q2/MAX(PAR(18),MZ**2))/3d0
     .    +DLOG(Q2/MAX(PAR(19),MZ**2))*2d0/3d0))

* Running Yukawa couplings:

* First: HTOP at MT, input: MT=top pole mass
      HTOP=MT/(1d0+g3t*COEF*16d0/3d0+176d0*(g3t*COEF)**2)/h1

* Conversion to DR_bar:
      HTOP=HTOP*(1d0-g3t*COEF*4d0/3d0+g2t*COEF*3d0/8d0)

* Second: HBOT at MT, input: MB(MT,MS_bar) from RUNMB(MT)
      HBOT=RUNMB(MT)/H2/(1d0+DELMB)

* Conversion to DR_bar:
      HBOT=HBOT*(1d0-g3t*COEF*4d0/3d0+g2t*COEF*3d0/8d0)

* Third: HTAU at MZ; use: MTAU(MZ)~1.775
      HTAU=1.775d0/h2

* Conversion to DR_bar:
      HTAU=HTAU*(1d0+g2t*COEF*3d0/8d0)

* Aux. quantities for the resummation of logs ~ht^2*LQT
* and      ~ht^2*LMAMT:

      ALSMT=G3T/(4d0*PI)
      ALSQ=G3S/(4d0*PI)
      LMAMT=DLOG(MAX(MA2,MT**2)/MT**2)
      ALSMA=ALSMT/(1d0+ALSMT/(4d0*PI)*(7d0*LMAMT-2d0*
     .    DLOG(MAX(MA2,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))
      DLA=(ALSMA/ALSMT)**(1d0/7d0)
      DLQA=(ALSQ/ALSMA)**(1d0/7d0)
      F1=1d0-9d0*SB2*HTOP**2*(1d0-DLA)/(8d0*PI*ALSMT)
      HTMA=HTOP*DLA**4/DSQRT(DABS(F1))
      F2=1d0-9d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA)

* HTOP at Q2=M_SUSY**2:
      HTOPS=HTOP*(1d0+7d0*COEF*G3T*DLOG(Q2/MT**2))**(-4d0/7d0)
     .    /DSQRT(F1*F2)
     .    *(1d0+COEF/4d0*((-26d0/9d0*g1t-6d0*g2t
     .    +2d0*HBOT**2+2d0*LS**2)*DLOG(Q2/MT**2)
     .    -(HBOT**2*(1d0-3d0*CB2)-2d0*HTAU**2*CB2)*LMAMT
     .    -8d0*g3t/3d0*(DLOG(MAX(PAR(7),PAR(22)**2,MT**2)/Q2)
     .    +DLOG(MAX(PAR(8),PAR(22)**2,MT**2)/Q2))
     .    -HTOP**2*(2d0*DLOG(MAX(PAR(7),MUS**2,MT**2)/Q2)
     .    +DLOG(MAX(PAR(8),MUS**2,MT**2)/Q2))
     .    -HBOT**2*DLOG(MAX(PAR(9),MUS**2,MT**2)/MT**2)
     .    -2d0*LS**2*DLOG(MAX(4d0*NUS**2,MUS**2,MT**2)/MT**2)
     .    +g1t*(-DLOG(MAX(PAR(20)**2,MUS**2,MT**2)/MT**2)
     .    -14d0/9d0*DLOG(MAX(PAR(8),PAR(20)**2,MT**2)/MT**2)
     .    +47d0/18d0*DLOG(MAX(PAR(7),PAR(20)**2,MT**2)/MT**2))
     .    +3d0*g2t*(-DLOG(MAX(PAR(21)**2,MUS**2,MT**2)/MT**2)
     .    +3d0/2d0*DLOG(MAX(PAR(7),PAR(21)**2,MT**2)/MT**2))))

* HBOT at Q2=M_SUSY**2:
      HBOTS=HBOT*(1d0+7d0*COEF*G3T*DLOG(Q2/MT**2))**(-4d0/7d0)
     .    *F1**(-1d0/6d0)
     .    *(1d0-3d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))**(-1d0/6d0)
     .    *(1d0+COEF/4d0*((-14d0/9d0*g1t-6d0*g2t+12d0*HBOT**2
     .    +2d0*HTAU**2+2d0*LS**2)*DLOG(Q2/MT**2)
     .    -(9d0*HBOT**2*SB2+2d0*HTAU**2*SB2)*LMAMT
     .    -8d0*g3t/3d0*(DLOG(MAX(PAR(7),PAR(22)**2,MT**2)/Q2)
     .    +DLOG(MAX(PAR(9),PAR(22)**2,MT**2)/Q2))
     .    -HBOT**2*(2d0*DLOG(MAX(PAR(7),MUS**2,MT**2)/MT**2)
     .    +DLOG(MAX(PAR(9),MUS**2,MT**2)/MT**2))
     .    -HTOP**2*DLOG(MAX(PAR(8),MUS**2,MT**2)/Q2)
     .    -2d0*LS**2*DLOG(MAX(4d0*NUS**2,MUS**2,MT**2)/MT**2)
     .    +g1t*(-DLOG(MAX(PAR(20)**2,MUS**2,MT**2)/MT**2)
     .    -8d0/9d0*DLOG(MAX(PAR(9),PAR(20)**2,MT**2)/MT**2)
     .    +47d0/18d0*DLOG(MAX(PAR(7),PAR(20)**2,MT**2)/MT**2))
     .    +3d0*g2t*(-DLOG(MAX(PAR(21)**2,MUS**2,MT**2)/MT**2)
     .    +3d0/2d0*DLOG(MAX(PAR(7),PAR(21)**2,MT**2)/MT**2))))

* HTAU at Q2=M_SUSY**2; assume 2 Higgs doublet beta function between
* MZ and M_SUSY
      HTAUS=HTAU*(1d0+COEF/2d0*(-15d0/4d0*g1z-9d0/4d0*g2z
     .    +3d0*HBOT**2+5d0/2d0*HTAU**2)*DLOG(Q2/MZ**2))

* Definition of the couplings squared Y(I) at M_SUSY

      Y(1)=g1s
      Y(2)=g2s
      Y(3)=g3s
      Y(4)=LS**2
      Y(5)=KS**2
      Y(6)=HTOPS**2
      Y(7)=HBOTS**2
      Y(8)=HTAUS**2

      X1=0d0
      X2=COEF*DLOG(MMESS**2/Q2)
      
      !WRITE(0,*)"CALL RGESGM"
      !WRITE(0,*)""
      !WRITE(0,*)"MSUSY =",DSQRT(Q2)
      !WRITE(0,*)"G1 =",Y(1)
      !WRITE(0,*)"G2 =",Y(2)
      !WRITE(0,*)"G3 =",Y(3)
      !WRITE(0,*)"L2 =",Y(4)
      !WRITE(0,*)"K2 =",Y(5)
      !WRITE(0,*)"HT2 =",Y(6)
      !WRITE(0,*)"HB2 =",Y(7)
      !WRITE(0,*)"HL2 =",Y(8)
      !WRITE(0,*)""

      CALL ODEINTGM(Y,NN,X1,X2,EPS,DERIVS,RKQSGM,IFAIL)
      
      !WRITE(0,*)"MMESS =",MMESS
      !WRITE(0,*)"G1 =",Y(1)
      !WRITE(0,*)"G2 =",Y(2)
      !WRITE(0,*)"G3 =",Y(3)
      !WRITE(0,*)"L2 =",Y(4)
      !WRITE(0,*)"K2 =",Y(5)
      !WRITE(0,*)"HT2 =",Y(6)
      !WRITE(0,*)"HB2 =",Y(7)
      !WRITE(0,*)"HL2 =",Y(8)
      !WRITE(0,*)""

      IF(IFAIL.GT.0)THEN
       IFAIL=11
      ELSE
       IFAIL=0
      ENDIF

* Couplings at the messenger scale

      G1MES=Y(1)
      G2MES=Y(2)
      G3MES=Y(3)
      LMES=Y(4)
      KMES=DSQRT(Y(5))*SIGNKAPPA
      HTOPMES=Y(6)
      HBOTMES=Y(7)
      HTAUMES=Y(8)

      END
