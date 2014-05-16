      SUBROUTINE RGESOFT(PAR,IFAIL)

*   Subroutine to integrate the RGEs for all 21 soft terms
*   from the SUSY scale Q2 up to MGUT, through a CALL of the
*   subroutine ODEINTS that is part of the file integs.f
*   Q2 is either computed in terms of the first generation squarks,
*   or put in by the user.
*
*   MGUT and the gauge/Yukawa couplings at Q2 are read in from
*   COMMON/SUSYCOUP.
*
*   The soft terms at Q2 are read in from PAR(*), the soft Higgs
*   masses squared at the scale QSTSB from COMMON/QMHIGGS.
*   They still have to be integrated to the SUSY scale Q2.
*   All sparticle threshold effects are taken into account.
*
*   CAUTION: Near an infrared quasi fixed point, very small changes in
*   A_top, M3 or M_HU^2 at low energy can cause a very large changes
*   of M_HU^2 and the stop masses at MGUT!
*
***********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,NN
      PARAMETER (NN=35)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN),PI,COEF
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION LS2,KS2,HTOPS2,HBOTS2,HTAUS2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MH1S,MH2S,MSS,Q2,QSTSB,ANOMQSTSB
      DOUBLE PRECISION LQSTSB,LM1QSTSB,LM2QSTSB
      DOUBLE PRECISION LM3QSTSB,LMQ3QSTSB,LMU3QSTSB
      DOUBLE PRECISION LMD3QSTSB,LMQQSTSB,LMUQSTSB,LMDQSTSB
      DOUBLE PRECISION LML3QSTSB,LME3QSTSB,LMLQSTSB,LMEQSTSB
      DOUBLE PRECISION RTOP,RBOT,RTAU,RL,RK,RG
      DOUBLE PRECISION M1,M2,M3,ALS,AKS,AT,AB,ATAU,AMUON
      DOUBLE PRECISION MH1Q,MH2Q,MSQ,MQ3,MU3,MD3,MQ
      DOUBLE PRECISION MU,MD,ML3,ME3,ML,ME
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
      DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT
      DOUBLE PRECISION MH1GUT,MH2GUT,MSGUT,MQ3GUT
      DOUBLE PRECISION MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT
      DOUBLE PRECISION MLGUT,MEGUT,XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY

      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/QMHIGGS/MH1Q,MH2Q,MSQ
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     .      HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     .      ATAUGUT,AMUGUT,MH1GUT,MH2GUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     .      MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      
      EXTERNAL DERIVSS,RKQSS

      IF(IFAIL.NE.0)RETURN

      !WRITE(0,*),"CALL RGESOFT"
      !WRITE(0,*),""

      EPS=1.D-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)
      
      LS2=PAR(1)**2
      KS2=PAR(2)**2
      HTOPS2=HTOPS**2
      HBOTS2=HBOTS**2
      HTAUS2=HTAUS**2
      ALS=PAR(5)
      AKS=PAR(6)
      AT=PAR(12)
      AB=PAR(13)
      ATAU=PAR(14)
      AMUON=PAR(25)
      M1=PAR(20)
      M2=PAR(21)
      M3=PAR(22)
      MQ3=PAR(7)
      MU3=PAR(8)
      MD3=PAR(9)
      ML3=PAR(10)
      ME3=PAR(11)
      MQ=PAR(15)
      MU=PAR(16)
      MD=PAR(17)
      ML=PAR(18)
      ME=PAR(19)

*  Useful logarithms:

      LQSTSB=DLOG(Q2/QSTSB)
      LM1QSTSB=DLOG(Q2/MAX(M1**2,QSTSB))
      LM2QSTSB=DLOG(Q2/MAX(M2**2,QSTSB))
      LM3QSTSB=DLOG(Q2/MAX(M3**2,QSTSB))
      LMQ3QSTSB=DLOG(Q2/MAX(MQ3,QSTSB))
      LMU3QSTSB=DLOG(Q2/MAX(MU3,QSTSB))
      LMD3QSTSB=DLOG(Q2/MAX(MD3,QSTSB))
      LMQQSTSB=DLOG(Q2/MAX(MQ,QSTSB))
      LMUQSTSB=DLOG(Q2/MAX(MU,QSTSB))
      LMDQSTSB=DLOG(Q2/MAX(MD,QSTSB))
      LML3QSTSB=DLOG(Q2/MAX(ML3,QSTSB))
      LME3QSTSB=DLOG(Q2/MAX(ME3,QSTSB))
      LMLQSTSB=DLOG(Q2/MAX(ML,QSTSB))
      LMEQSTSB=DLOG(Q2/MAX(ME,QSTSB))

* The Higgs masses squared:

      ANOMQSTSB=G1S*(-MH2Q*LQSTSB+MQ3*LMQ3QSTSB
     .  -2d0*MU3*LMU3QSTSB
     .  +MD3*LMD3QSTSB+2d0*(MQ*LMQQSTSB
     .  -2d0*MU*LMUQSTSB+MD*LMDQSTSB)
     .  +ME3*LME3QSTSB-ML3*LML3QSTSB
     .  +2d0*(ME*LMEQSTSB-ML*LMLQSTSB))

      RTOP=HTOPS2*(MQ3*LMQ3QSTSB+MU3*LMU3QSTSB
     .       +AT**2*LQSTSB)

      RBOT=HBOTS2*(MH2Q*LQSTSB+MQ3*LMQ3QSTSB+MD3*LMD3QSTSB
     .       +AB**2*LQSTSB)

      RTAU=HTAUS2*(MH2Q*LQSTSB+ML3*LML3QSTSB+ME3*LME3QSTSB
     .       +ATAU**2*LME3QSTSB)

      RL=LS2*(MH2Q+MSQ+ALS**2)*LQSTSB

      RK=KS2*(3d0*MSQ+AKS**2)*LQSTSB

      RG=G1S*M1**2*LM1QSTSB+3d0*G2S*M2**2*LM2QSTSB

      MH1S=(MH1Q+COEF*(RL+3d0*RTOP-RG+ANOMQSTSB/2d0))
     .      *(Q2/QSTSB)**(COEF*(LS2+3d0*HTOPS2+G1S/2d0))

      MH2S=MH2Q+COEF*(RL+3d0*RBOT+RTAU-RG-ANOMQSTSB/2d0
     .      +(LS2-G1S/2d0)*MH1Q*LQSTSB)

      MSS=MSQ+COEF*(2d0*RK+2d0*RL+2d0*LS2*MH1Q*LQSTSB)

* Definition of the couplings squared Y(I) at Q2=M_SUSY

      Y(1)=g1s
      Y(2)=g2s
      Y(3)=g3s
      Y(4)=LS2
* NOTE: Y(5)=K, NOT K**2
      Y(5)=PAR(2)
      Y(6)=HTOPS2
      Y(7)=HBOTS2
      Y(8)=HTAUS2

* Definition of the soft terms Y(I) at Q2=M_SUSY

      Y(9)=M1
      Y(10)=M2
      Y(11)=M3
      Y(12)=ALS
      Y(13)=AKS
      Y(14)=AT
      Y(15)=AB
      Y(16)=ATAU
      Y(17)=MH1S
      Y(18)=MH2S
      Y(19)=MSS
      Y(20)=MQ3
      Y(21)=MU3
      Y(22)=MD3
      Y(23)=MQ
      Y(24)=MU
      Y(25)=MD
      Y(26)=ML3
      Y(27)=ME3
      Y(28)=ML
      Y(29)=ME
      Y(30)=XIFSUSY
      Y(31)=XISSUSY
      Y(32)=MUPSUSY
      Y(33)=MSPSUSY
      Y(34)=M3HSUSY
      Y(35)=AMUON
      
      X1=0d0
      X2=COEF*DLOG(MGUT**2/Q2)

      !WRITE(0,*),"MSUSY =",DSQRT(Q2)
      !WRITE(0,*),"G1 =",5d0/3d0*Y(1)
      !WRITE(0,*),"G2 =",Y(2)
      !WRITE(0,*),"G3 =",Y(3)
      !WRITE(0,*),"L =",Y(4)
      !WRITE(0,*),"K =",Y(5)
      !WRITE(0,*),"HT =",Y(6)
      !WRITE(0,*),"HB =",Y(7)
      !WRITE(0,*),"HL =",Y(8)
      !WRITE(0,*),"M1 =",Y(9)
      !WRITE(0,*),"M2 =",Y(10)
      !WRITE(0,*),"M3 =",Y(11)
      !WRITE(0,*),"AL =",Y(12)
      !WRITE(0,*),"AK =",Y(13)
      !WRITE(0,*),"ATOP =",Y(14)
      !WRITE(0,*),"ABOT =",Y(15)
      !WRITE(0,*),"ATAU =",Y(16)
      !WRITE(0,*),"AMUON =",Y(35)
      !WRITE(0,*),"MH1 =",Y(17)
      !WRITE(0,*),"MH2 =",Y(18)
      !WRITE(0,*),"MS =",Y(19)
      !WRITE(0,*),"MQ3 =",Y(20)
      !WRITE(0,*),"MU3 =",Y(21)
      !WRITE(0,*),"MD3 =",Y(22)
      !WRITE(0,*),"MQ =",Y(23)
      !WRITE(0,*),"MU =",Y(24)
      !WRITE(0,*),"MD =",Y(25)
      !WRITE(0,*),"ML3 =",Y(26)
      !WRITE(0,*),"ME3 =",Y(27)
      !WRITE(0,*),"ML =",Y(28)
      !WRITE(0,*),"ME =",Y(29)
      !WRITE(0,*),"XIF =",Y(30)
      !WRITE(0,*),"XIS =",Y(31)
      !WRITE(0,*),"MUP =",Y(32)
      !WRITE(0,*),"MSP =",Y(33)
      !WRITE(0,*),"M3H =",Y(34)
      !WRITE(0,*),""

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)
      
      !WRITE(0,*),"MGUT =",MGUT
      !WRITE(0,*),"G1GUT =",5d0/3d0*Y(1)
      !WRITE(0,*),"G2GUT =",Y(2)
      !WRITE(0,*),"G3GUT =",Y(3)
      !WRITE(0,*),"LGUT =",Y(4)
      !WRITE(0,*),"KGUT =",Y(5)
      !WRITE(0,*),"HTGUT =",Y(6)
      !WRITE(0,*),"HBGUT =",Y(7)
      !WRITE(0,*),"HLGUT =",Y(8)
      !WRITE(0,*),"M1GUT =",Y(9)
      !WRITE(0,*),"M2GUT =",Y(10)
      !WRITE(0,*),"M3GUT =",Y(11)
      !WRITE(0,*),"ALGUT =",Y(12)
      !WRITE(0,*),"AKGUT =",Y(13)
      !WRITE(0,*),"ATOPGUT =",Y(14)
      !WRITE(0,*),"ABOTGUT =",Y(15)
      !WRITE(0,*),"ATAUGUT =",Y(16)
      !WRITE(0,*),"AMUGUT =",Y(35)
      !WRITE(0,*),"MH1GUT =",Y(17)
      !WRITE(0,*),"MH2GUT =",Y(18)
      !WRITE(0,*),"MSGUT =",Y(19)
      !WRITE(0,*),"MQ3GUT =",Y(20)
      !WRITE(0,*),"MU3GUT =",Y(21)
      !WRITE(0,*),"MD3GUT =",Y(22)
      !WRITE(0,*),"MQGUT =",Y(23)
      !WRITE(0,*),"MUGUT =",Y(24)
      !WRITE(0,*),"MDGUT =",Y(25)
      !WRITE(0,*),"ML3GUT =",Y(26)
      !WRITE(0,*),"ME3GUT =",Y(27)
      !WRITE(0,*),"MLGUT =",Y(28)
      !WRITE(0,*),"MEGUT =",Y(29)
      !WRITE(0,*),"XIFGUT =",Y(30)
      !WRITE(0,*),"XISGUT =",Y(31)
      !WRITE(0,*),"MUPGUT =",Y(32)
      !WRITE(0,*),"MSPGUT =",Y(33)
      !WRITE(0,*),"M3HGUT =",Y(34)
      !WRITE(0,*),""

      IF(IFAIL.NE.0)THEN
       !WRITE(0,*),"IFAIL =",IFAIL
       !WRITE(0,*),""
       !WRITE(0,*),""
       IFAIL=12
       RETURN
      ENDIF
      !WRITE(0,*),""
      
* Couplings at the GUT scale

      G1GUT=Y(1)
      G2GUT=Y(2)
      G3GUT=Y(3)
      LGUT=Y(4)
      KGUT=Y(5)
      HTOPGUT=Y(6)
      HBOTGUT=Y(7)
      HTAUGUT=Y(8)
      
* Soft terms at the GUT scale

      M1GUT=Y(9)
      M2GUT=Y(10)
      M3GUT=Y(11)
      ALGUT=Y(12)
      AKGUT=Y(13)
      ATGUT=Y(14)
      ABGUT=Y(15)
      ATAUGUT=Y(16)
      AMUGUT=Y(35)
      MH1GUT=Y(17)
      MH2GUT=Y(18)
      MSGUT=Y(19)
      MQ3GUT=Y(20)
      MU3GUT=Y(21)
      MD3GUT=Y(22)
      MQGUT=Y(23)
      MUGUT=Y(24)
      MDGUT=Y(25)
      ML3GUT=Y(26)
      ME3GUT=Y(27)
      MLGUT=Y(28)
      MEGUT=Y(29)
      XIFGUT=Y(30)
      XISGUT=Y(31)
      MUPGUT=Y(32)
      MSPGUT=Y(33)
      M3HGUT=Y(34)
      
      END
