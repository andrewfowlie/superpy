      SUBROUTINE MINIMIZE(PAR,CHECK)

**********************************************************************      
* Subroutine to compute
*   MAFLAG=-5 :
*     - MHU, MHD, MS
*   MAFLAG=/=-5 :
*     - KAPPA (if MAFLAG=-1 or -2) or XIF (if MAFLAG=-3 or -4)
*     - MS (if MAFLAG=-1 or -3) or XIS (if MAFLAG=-2 or -4)
*     - MU 
*  at the scales QSTSB and Q2
*
* MU at QSTSB is stored in COMMON/QPAR, MU at Q2 in PAR(4)
* KAPPA at QSTSB is stored in COMMON/QPAR, KAPPA at Q2 in PAR(2)
* MS at QSTSB is stored in COMMON/QMHIGGS, MS at Q2 in OMMON/SUSYMH
* XIF at QSTSB is stored in COMMON/QEXT, XIF at Q2 in OMMON/SUSYEXT
* XIS at QSTSB is stored in COMMON/QEXT, XIS at Q2 in OMMON/SUSYEXT
* MA at QSTSB is stored in PAR(23)
* MP at QSTSB is stored in PAR(24)
*
**********************************************************************

      IMPLICIT NONE

      INTEGER OMGFLAG,MAFLAG

      DOUBLE PRECISION PAR(*),CHECK,SIGMU,pi,At,Ab,COEF
      DOUBLE PRECISION B,MUOLD,KOLD,MSSOLD,XIFOLD,XISOLD
      DOUBLE PRECISION mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
      DOUBLE PRECISION ct,fmt1,fmt2,fmt,gmt
      DOUBLE PRECISION cb,fmb1,fmb2,fmb,gmb
      DOUBLE PRECISION MH1S,MH2S,MSS,QSTSB
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION RT1,RT2,RB1,RB2,AM,BM,CM,DETM
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION LS2,KS2,HTOPS2,HBOTS2,HTAUS2
      DOUBLE PRECISION M1,M2,M3,ATAU,MQ3,MU3,MD3
      DOUBLE PRECISION MQ,MU,MD,ML3,ME3,ML,ME,AL
      DOUBLE PRECISION Q2,MH1Q,MH2Q,MSQ
      DOUBLE PRECISION LQSTSB,LM1QSTSB,LM2QSTSB
      DOUBLE PRECISION LM3QSTSB,LMQ3QSTSB,LMU3QSTSB
      DOUBLE PRECISION LMD3QSTSB,LMQQSTSB,LMUQSTSB,LMDQSTSB
      DOUBLE PRECISION LML3QSTSB,LME3QSTSB,LMLQSTSB,LMEQSTSB
      DOUBLE PRECISION RL,RTOP,RBOT,RTAU,RG,ANOMQSTSB,MUFAIL
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/STSBSCALE/QSTSB
      COMMON/RADCOR/mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/SIGMU/SIGMU
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/RENSCALE/Q2
      COMMON/QMHIGGS/MH1Q,MH2Q,MSQ
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/DETM/DETM
      COMMON/MUFAIL/MUFAIL

      pi=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

!      WRITE(0,*)"CALL MINIMIZE"
!      WRITE(0,*)""

* (S)top/(S)/bottom loop corrections to the minimization equations:

      At=PAR(12)
      ct=3d0*htq**2*COEF
      fmt1=mst1*(DLOG(mst1/QSTSB)-1d0)
      fmt2=mst2*(DLOG(mst2/QSTSB)-1d0)
      fmt=mtopq**2*(DLOG(mtopq**2/QSTSB)-1d0)
      IF(mst1-mst2.NE.0d0)THEN
       gmt=(fmt2-fmt1)/(mst2-mst1)
      ELSE
       gmt=DLOG(mst1/QSTSB)
      ENDIF

      Ab=PAR(13)
      cb=3d0*hbq**2*COEF
      fmb1=msb1*(DLOG(msb1/QSTSB)-1d0)
      fmb2=msb2*(DLOG(msb2/QSTSB)-1d0)
      fmb=mbotq**2*(DLOG(mbotq**2/QSTSB)-1d0)
      IF(msb1-msb2.NE.0d0)THEN
       gmb=(fmb2-fmb1)/(msb2-msb1)
      ELSE
       gmb=DLOG(msb1/QSTSB)
      ENDIF

      RT1=CT*(FMT1+FMT2-2d0*FMT+AT*XT*GMT)
      RT2=CT*XT*GMT
      RB1=CB*(FMB1+FMB2-2d0*FMB+AB*XB*GMB)
      RB2=CB*XB*GMB

* If MU is an input parameters MH1Q, MH2Q, MSQ are computed

      IF(MAFLAG.EQ.-5)THEN
       B=ALQ+NUQ
       MH1Q=-LQ**2*H2Q**2 - MUQ**2 + MUQ*B/TANBQ
     .      +GQ/2d0*(H2Q**2-H1Q**2)
     .      -RT1+RB2*MUQ/TANBQ
       MH2Q=-LQ**2*H1Q**2 - MUQ**2 + MUQ*B*TANBQ
     .      +GQ/2d0*(H1Q**2-H2Q**2)
     .      -RB1+RT2*MUQ*TANBQ
       MSQ=-LQ**2*(H1Q**2+H2Q**2) - 2d0*NUQ**2
     .     +LQ**2*H1Q*H2Q/MUQ*(ALQ+2d0*NUQ+MUPQ) - NUQ*AKQ
     .     -XIFQ*(2d0*KQ+LQ*MUPQ/MUQ)-MUPQ**2-3d0*MUPQ*NUQ
     .     -MSPQ-LQ*XISQ/MUQ
     .     +RT2*LQ**2*H1Q*H2Q/MUQ
     .     +RB2*LQ**2*H1Q*H2Q/MUQ
      ENDIF

* MH1, MH2 are integrated between Q2 and QSTSB:

      LS2=PAR(1)**2
      KS2=PAR(2)**2
      HTOPS2=HTOPS**2
      HBOTS2=HBOTS**2
      HTAUS2=HTAUS**2

      M1=PAR(20)
      M2=PAR(21)
      M3=PAR(22)
      ATAU=PAR(14)
      AL=PAR(5)

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

* Useful logarithms:

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

      ANOMQSTSB=G1S*(-MH2S*LQSTSB+MQ3*LMQ3QSTSB
     .  -2d0*MU3*LMU3QSTSB
     .  +MD3*LMD3QSTSB+2d0*(MQ*LMQQSTSB
     .  -2d0*MU*LMUQSTSB+MD*LMDQSTSB)
     .  +ME3*LME3QSTSB-ML3*LML3QSTSB
     .  +2d0*(ME*LMEQSTSB-ML*LMLQSTSB))

      RTOP=HTOPS2*(MQ3*LMQ3QSTSB+MU3*LMU3QSTSB
     .       +AT**2*LQSTSB)

      RBOT=HBOTS2*(MH2S*LQSTSB+MQ3*LMQ3QSTSB+MD3*LMD3QSTSB
     .       +AB**2*LQSTSB)

      RTAU=HTAUS2*(MH2S*LQSTSB+ML3*LML3QSTSB+ME3*LME3QSTSB
     .       +ATAU**2*LME3QSTSB)

      RL=LS2*(MH2S+MSS+AL**2)*LQSTSB

      RG=G1S*M1**2*LM1QSTSB+3d0*G2S*M2**2*LM2QSTSB

* Running MH1, MH2
*   from QSTSB to Q2:

      IF(MAFLAG.EQ.-5)THEN

      MH1S=(MH1Q+COEF*(RL+3d0*RTOP-RG+ANOMQSTSB/2d0))
     .     *(Q2/QSTSB)**(COEF*(LS2+3d0*HTOPS2+G1S/2d0))
      MH2S=MH2Q+COEF*(RL+3d0*RBOT+RTAU-RG-ANOMQSTSB/2d0
     .     +(LS2-G1S/2d0)*MH1S*LQSTSB)
      MSS=MSQ+COEF*(2d0*LQ**2*(MH1Q+MH2Q+MSQ+ALQ**2)
     .    +2d0*KQ**2*(3d0*MSQ+AKQ**2))*LQSTSB

!      WRITE(0,*)"MH1S =",MH1S
!      WRITE(0,*)"MH2S =",MH2S
!      WRITE(0,*)"MSS =",MSS

      CHECK=0d0

*   or from Q2 to QSTSB:

      ELSE

      MH1Q=MH1S*(QSTSB/Q2)**(COEF*(LS2+3d0*HTOPS2+G1S/2d0))
     .     -COEF*(RL+3d0*RTOP-RG+ANOMQSTSB/2d0)
      MH2Q=MH2S-COEF*(RL+3d0*RBOT+RTAU-RG-ANOMQSTSB/2d0
     .     +(LS2-G1S/2d0)*MH1S*LQSTSB)

* Store the previous value of MU, (K or XIF) and (MS or XIS) for CHECK:

      MUOLD=PAR(4)
      KOLD=PAR(2)
      MSSOLD=MSS
      XIFOLD=XIF
      XISOLD=XIS

* Prepare the computation of MU:

      am=tanbq-1d0/tanbq
      bm=rt2-rb2
      cm=gq/2d0*(H1Q**2+H2Q**2)*(tanbq-1d0/tanbq)
     .     +(rt1+MH1Q)*tanbq-(rb1+MH2Q)/tanbq
      detm=bm**2-4d0*am*cm

* MU at QSTSB is computed here:

      IF(DETM.GE.0d0)THEN
       MUQ=(-BM+SIGMU*DSQRT(DETM))/(2d0*AM)
      ELSE
!       WRITE(0,*)"DETM < 0"
!       WRITE(0,*)""
       MUQ=MUFAIL
      ENDIF

      IF(DABS(MUQ).LT.DABS(MUFAIL)) THEN
!       WRITE(0,*)"DETM = 0"
!       WRITE(0,*)""
       DETM=0d0
       MUQ=MUFAIL
      ENDIF

* PAR(4) = MU at Q2, MU at QSTSB = MUQ is stored in COMMON/QPAR     

      PAR(4)=MUQ*(1d0+COEF/2d0*(2d0*LQ**2+3d0*(HTQ**2+HBQ**2)
     .       +(MTAU/H2Q)**2-G1Q-3d0*G2Q)*LQSTSB)

!      WRITE(0,*)"MUS =",MUOLD,PAR(4)

* Effective B parameter:

      B=TANBQ/((1d0+TANBQ**2)*MUQ)*(MH1Q+MH2Q+2d0*MUQ**2
     .  +LQ**2*(H1Q**2+H2Q**2)+RT1+RB1-MUQ*(RT2*TANBQ+RB2/TANBQ))

* If MAFLAG=-1,-2:
*  Compute KQ = Kappa at QSTSB, then PAR(2) = KAPPA at Q2

      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN

       KQ=LQ*(B-ALQ-MUPQ-(LQ*XIFQ+M3HQ)/MUQ)/MUQ
       PAR(2)=KQ*(1d0+3d0*COEF*(LQ**2+KQ**2)*LQSTSB)
       NUQ=MUQ*KQ/LQ
!      WRITE(0,*)"KS =",KOLD,PAR(2)

*  Else compute XIF

      ELSE

       XIFQ=(MUQ*(B-ALQ-NUQ-MUPQ)-M3HQ)/LQ
       XIF=XIFQ*(1d0+COEF*(LQ**2+KQ**2)*LQSTSB)
!      WRITE(0,*)"XIF =",XIFOLD,XIF

      ENDIF

* If MAFLAG=-1,-3:
*  Compute MSSX = MS at QSTSB, then MSS = MS at Q2

      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN

       MSQ=-LQ**2*(H1Q**2+H2Q**2) - 2d0*NUQ**2
     .    +LQ**2*H1Q*H2Q/MUQ*(ALQ+2d0*NUQ+MUPQ) - NUQ*AKQ
     .    -XIFQ*(2d0*KQ+LQ*MUPQ/MUQ)-MUPQ**2-3d0*MUPQ*NUQ
     .    -MSPQ-LQ*XISQ/MUQ
     .    +ct*LQ**2*H1Q*H2Q/MUQ*Xt*gmt
     .    +cb*LQ**2*H1Q*H2Q/MUQ*Xb*gmb
       MSS=MSQ+COEF*(2d0*LQ**2*(MH1Q+MH2Q+MSQ+ALQ**2)
     .      +2d0*KQ**2*(3d0*MSQ+AKQ**2))*LQSTSB

!       WRITE(0,*)"MSS =",MSSOLD,MSS

* Else compute XIS

      ELSE

       MSQ=MSS-COEF*(2d0*LQ**2*(MH1Q+MH2Q+MSS+ALQ**2)
     .    +2d0*KQ**2*(3d0*MSS+AKQ**2))*LQSTSB
       XISQ=LQ*(H1Q*H2Q*(ALQ+2d0*NUQ+MUPQ)-MUQ*(H1Q**2+H2Q**2))
     .     -MUQ/LQ*(MSQ+2d0*NUQ**2+NUQ*AKQ
     .     +XIFQ*(2d0*KQ+LQ*MUPQ/MUQ)
     .     +MUPQ**2+3d0*MUPQ*NUQ+MSPQ
     .     -ct*LQ**2*H1Q*H2Q/MUQ*Xt*gmt
     .     -cb*LQ**2*H1Q*H2Q/MUQ*Xb*gmb)
       XIS=XISQ+COEF*(LQ**2*(XISQ+2d0*ALQ*XIF)
     .    +KQ**2*(XISQ+2d0*AKQ*XIFQ)
     .    +2d0*LQ*M3HQ*(ALQ+MUPQ)
     .    +KQ*MSPQ*(AKQ+MUPQ))*LQSTSB

!       WRITE(0,*)"XIS =",XISOLD,XIS

      ENDIF

*  CHECK of MU, K or XIF, MS or XIS against the previous value

      CHECK=(MUOLD-PAR(4))**2/(MUOLD**2+1.D2)
      CHECK=CHECK+(KOLD-PAR(2))**2/(KOLD**2+1.D-12)
      CHECK=CHECK+(XIFOLD-XIF)**2/(XIFOLD**2+1.D4)
      CHECK=CHECK+(MSSOLD-MSS)**2/(MSSOLD**2+1.D4)
      CHECK=CHECK+(XISOLD-XIS)**2/(XISOLD**2+1.D6)

      ENDIF

* Approximate value for the MSSM-like CP odd Higgs mass MA:

      PAR(23)=DSQRT(MAX((MUQ*B+M3HQ+MUQ*MUPQ+LQ*XIFQ)
     .       *(tanbQ+1d0/tanbQ),1d0))

* Approximate value for the singlet-like CP odd Higgs mass MP:

      PAR(24)=DSQRT(MAX(LQ**2*(B+3d0*NUQ+MUPQ)*H1Q*H2Q/MUQ
     .       -3d0*AKQ*NUQ-MUPQ*NUQ-2d0*MSPQ-4d0*KQ*XIFQ-LQ/MUQ
     .       *(XIFQ*MUPQ+XISQ),1d0))

!      WRITE(0,*)"MA =",PAR(23)
!      WRITE(0,*)"MP =",PAR(24)
!      WRITE(0,*)""
!      WRITE(0,*)"CHECK =",CHECK
!      WRITE(0,*)""
!      WRITE(0,*)""

      END
