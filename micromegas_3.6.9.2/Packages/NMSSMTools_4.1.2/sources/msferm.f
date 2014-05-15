      SUBROUTINE MSFERM(PAR,IFAIL,FLAG)
      
*******************************************************************
* Subroutine to compute the squark pole masses/mixing angles and
* slepton masses/mixing angles (with thanks to S. Kraml).
*
* The (large tan(beta)) SUSY correction DELMB is computed here.
*
* Negative squark masses squared give IFAIL = 8.
*
* The running masses of the 3rd generation of squarks are first
* evaluated at a scale QSTSB taking all possible thresholds
* (i.e. possibly large logs) into account. Pole masses are computed
* including 1 loop (S)QCD rad. corrs.
* No electroweak rad. corrs. are included at present.
*
* Concerning the 1st and 2nd generation squarks, it is assumed that
* the input parameters PAR(15)=MQ, PAR(16)=MU, PAR(17)=MD are the
* running squark masses at the SUSY scale Q2 ~ MQ^2 ~ MU^2 ~ MD^2.
*
* The conventions for the cosines/sines CSF/SSF (F=T, B, L=Stau) of the
* 3rd generation sfermion mixing angle are
* F1 = CSF*FL + SSF*FR, F2 = CSF*FR - SSF*FL
* where FR, FL are the left/right handed weak eigenstates, and
* F1, F2 are the physical states ordered in mass (M_F1 < M_F2)
*******************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,OMGFLAG,MAFLAG,FLAG

      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION cor,Q2,QSTSB,pi
      DOUBLE PRECISION Wt,At,mstL,mstR
      DOUBLE PRECISION Wb,Ab,msbL,msbR
      DOUBLE PRECISION Wl,Atau,AMUON,Xl,mslL,mslR
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW,X
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION MGL,nen,fac,DMSQUARK
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION M1,M2,M3,MQ3,MU3,MD3,ML3,ME3,MQ,MU,MD,ML,ME
      DOUBLE PRECISION MQ3P,MU3P,MD3P,MH1S,MH2S,MSS,ANOMQSTSB,B
      DOUBLE PRECISION LQSTSB,LM1QSTSB,LM2QSTSB
      DOUBLE PRECISION LM3QSTSB,LMQ3QSTSB,LMU3QSTSB
      DOUBLE PRECISION LMD3QSTSB,LMQQSTSB,LMUQSTSB,LMDQSTSB
      DOUBLE PRECISION LML3QSTSB,LME3QSTSB,LMLQSTSB,LMEQSTSB
      DOUBLE PRECISION RTOP,RBOT,COEF,ATP,ABP
      DOUBLE PRECISION INTEG,DELMB,C2T,C2B
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTQ,MBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/RADCOR2/MQ3P,MU3P,MD3P,ATP,ABP
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/DELMB/DELMB
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/SUSYMH/MH1S,MH2S,MSS

      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

!      WRITE(0,*)"CALL MSFERM"
!      WRITE(0,*)""

      AT=PAR(12)
      AB=PAR(13)
      ATAU=PAR(14)
      AMUON=PAR(25)

*   In order to run the squark and slepton masses from Q2 to QSTSB
*   including all thresholds, tree level expressions for the soft
*   Higgs masses and various logarithms are needed:

      IF(MAFLAG.GE.0)THEN
       B=ALQ+NUQ
       MH1S= -LQ**2*H2Q**2 - MUQ**2 + MUQ*B/TANBQ
     .      + GQ/2d0*(H2Q**2-H1Q**2)
       MH2S= -LQ**2*H2Q**2 - MUQ**2 + MUQ*B*TANBQ
     .      + GQ/2d0*(H1Q**2-H2Q**2)
      ENDIF

      M1=PAR(20)
      M2=PAR(21)
      M3=PAR(22)
      MGL=PAR(22)
      
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
      
      ANOMQSTSB=G1Q*(MH1S*LQSTSB-MH2S*LQSTSB+MQ3*LMQ3QSTSB
     .  -2d0*MU3*LMU3QSTSB
     .  +MD3*LMD3QSTSB+2d0*(MQ*LMQQSTSB
     .  -2d0*MU*LMUQSTSB+MD*LMDQSTSB)
     .  +ME3*LME3QSTSB-ML3*LML3QSTSB
     .  +2d0*(ME*LMEQSTSB-ML*LMLQSTSB))

*   Running stop/sbottom masses squared and mixings:
*   On input, they are defined at the scale Q2.
*   First, they have to be evaluated at the scale QSTSB
*   (as used in the calculation of the pole masses)
*   taking all possible thresholds into account
*   (with the exception of the gluino threshold that
*   is taken into account in the calculation of the pole masses).

      RTOP=HTQ**2*(MH1S*LQSTSB+MQ3*LMQ3QSTSB+MU3*LMU3QSTSB
     .       +AT**2*LQSTSB)

      RBOT=HBQ**2*(MH2S*LQSTSB+MQ3*LMQ3QSTSB+MD3*LMD3QSTSB
     .       +AB**2*LQSTSB)

      MQ3P=MQ3-COEF*(RTOP+RBOT
     .       -G1Q*M1**2*LM1QSTSB/9d0-3d0*G2Q*M2**2*LM2QSTSB
     .       +ANOMQSTSB/3d0)

      MU3P=MU3-COEF*(2d0*RTOP
     .       -16d0/9d0*G1Q*M1**2*LM1QSTSB-2d0/3d0*ANOMQSTSB)

      MD3P=MD3-COEF*(2d0*RBOT
     .       -4d0/9d0*G1Q*M1**2*LM1QSTSB+1d0/3d0*ANOMQSTSB)

*  Integrate the trilinears from Q2 to QSTSB:

      ATP=AT-COEF*((6d0*HTQ**2*AT+HBQ**2*AB+LQ**2*ALQ)*LQSTSB
     .    +13d0*G1Q/3d0*M1*LM1QSTSB+3d0*G2Q*M2*LM2QSTSB
     .    +64d0*PI*ALSQ/3d0*M3*LM3QSTSB)

      ABP=AB-COEF*((6d0*HBQ**2*AB+HTQ**2*AT+LQ**2*ALQ)*LQSTSB
     .    +7d0*G1Q/9d0*M1*LM1QSTSB+3d0*G2Q*M2*LM2QSTSB
     .    +64d0*PI*ALSQ/3d0*M3*LM3QSTSB)

*   Running stop masses squared and mixings

      mstL= MQ3P + mtq**2 + (gQ/2d0-g1Q/3d0)*(h2q**2-h1q**2)
      mstR= MU3P + mtq**2 + g1Q/3d0*(h2q**2-h1q**2)
      Xt= ATP-muq/tanbq
      Wt= DSQRT( (mstL-mstR)**2 + 4d0*(Xt*mtq)**2)
      RMST1= 0.5d0*(mstL+mstR-Wt)
      RMST2= 0.5d0*(mstL+mstR+Wt)

      IF(RMST1.LE.1d1)THEN
       IFAIL=8
       nen=0d0
       IF(FLAG.EQ.0)RMST1=1d1
!      WRITE(0,*)"MST1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ELSE
       nen= DSQRT((mstL-RMST1)**2 + (Xt*mtq)**2)
      ENDIF

      IF(RMST2.LE.1d1)THEN
       IFAIL=8
       IF(FLAG.EQ.0)RMST2=1d1
!      WRITE(0,*)"MST2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(nen.EQ.0d0)THEN
       IF(mstL.LE.mstR)THEN
        CST= 1d0
       ELSE
        CST= 0d0
       ENDIF
      ELSE
       CST= -mtq*Xt/nen
      ENDIF
      S2T=2d0*CST*DSQRT(1d0-CST**2)

*   Running sbottom masses squared and mixings

      msbL= MQ3P + mbq**2 - (gQ/2d0-g1Q/6d0)*(h2q**2-h1q**2)
      msbR= MD3P + mbq**2 - g1Q/6d0*(h2q**2-h1q**2)
      Xb= ABP-muq*tanbq
      Wb= DSQRT((msbL-msbR)**2 + 4d0*(Xb*mbq)**2)
      RMSB1= 0.5d0*(msbL+msbR-Wb)
      RMSB2= 0.5d0*(msbL+msbR+Wb)

      IF(RMSB1.LE.1d1)THEN
       IFAIL=8
       nen=0d0
       IF(FLAG.EQ.0)RMSB1=1d1
!      WRITE(0,*)"MSB1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ELSE
       nen= DSQRT((msbL-RMSB1)**2 + (Xb*mbq)**2)
      ENDIF

      IF(RMSB2.LE.1d1)THEN
       IFAIL=8
       IF(FLAG.EQ.0)RMSB2=1d1
!      WRITE(0,*)"MSB2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(nen.EQ.0d0)THEN
       IF(msbL.LE.msbR)THEN
        CSB= 1d0
       ELSE
        CSB= 0d0
       ENDIF
      ELSE
       CSB= -mbq*Xb/nen
      ENDIF
      S2B=2*CSB*DSQRT(1d0-CSB**2)

*   Pole Squark Masses squared

      fac=ALSQ/(3d0*PI)
      MST1=RMST1
      MST2=RMST2
      MSB1=RMSB1
      MSB2=RMSB2

      IF(MST1.GT.1d1)THEN
       MST1=RMST1-fac*DMSQUARK(1,RMST1,RMST2,S2T,mtq,MGL,QSTSB)
       MST2=RMST2-fac*DMSQUARK(2,RMST1,RMST2,S2T,mtq,MGL,QSTSB)
      ENDIF
      IF(MSB1.GT.1d1)THEN
       MSB1=RMSB1-fac*DMSQUARK(1,RMSB1,RMSB2,S2B,mbq,MGL,QSTSB)
       MSB2=RMSB2-fac*DMSQUARK(2,RMSB1,RMSB2,S2B,mbq,MGL,QSTSB)
      ENDIF

      IF(MST1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MST1=1d2
!      WRITE(0,*)"MST1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MST2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MST2=1d2
!      WRITE(0,*)"MST2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSB1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSB1=1d2
!      WRITE(0,*)"MSB1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSB2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSB2=1d2
!      WRITE(0,*)"MSB2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

*   First generation squark masses squared

      MUR=MU+g1Q/3d0*(h2q**2-h1q**2)
      MUL=MQ+(gQ/2d0-g1Q/3d0)*(h2q**2-h1q**2)
      MDR=MD-g1Q/6d0*(h2q**2-h1q**2)
      MDL=MQ+(-gQ/2d0+g1Q/6d0)*(h2q**2-h1q**2)

      IF(MUR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUR=1d2
!      WRITE(0,*)"MUR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MUL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUL=1d2
!      WRITE(0,*)"MUL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDR=1d2
!      WRITE(0,*)"MDR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDL=1d2
!      WRITE(0,*)"MDL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

*   Pole Squark Masses squared
*   On input they are defined at the scale Q2 that is assumed to not
*   too far from their actual masses
*   (Threshold effects are neglected)

      IF(MUR.GT.1d2)THEN
       X=M3**2/MUR
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MUR)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MUR)
       ENDIF
       MUR=MUR*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUL.GT.1d2)THEN
       X=M3**2/MUL
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MUL)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MUL)
       ENDIF
       MUL=MUL*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUL.GT.1d2)THEN
       X=M3**2/MDR
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MDR)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MDR)
       ENDIF
       MDR=MDR*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUL.GT.1d2)THEN
       X=M3**2/MDL
       IF(X.NE.1d0)THEN
        cor=1d0+3d0*X+1d0/2d0*(1d0-X)**2*DLOG((1d0-X)**2)
     .     -X**2*DLOG(X)+2d0*X*DLOG(Q2/MDL)
       ELSE
        cor=4d0+2d0*DLOG(Q2/MDL)
       ENDIF
       MDL=MDL*(1d0+2d0*fac*cor)
      ENDIF

      IF(MUR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUR=1d2
!      WRITE(0,*)"MUR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MUL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MUL=1d2
!      WRITE(0,*)"MUL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDR=1d2
!      WRITE(0,*)"MDR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MDL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MDL=1d2
!      WRITE(0,*)"MDL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MST1.GE.0d0)MST1=DSQRT(MST1)
      IF(MST2.GE.0d0)MST2=DSQRT(MST2)
      IF(MSB1.GE.0d0)MSB1=DSQRT(MSB1)
      IF(MSB2.GE.0d0)MSB2=DSQRT(MSB2)
      IF(MUL.GE.0d0)MUL=DSQRT(MUL)
      IF(MUR.GE.0d0)MUR=DSQRT(MUR)
      IF(MDL.GE.0d0)MDL=DSQRT(MDL)
      IF(MDR.GE.0d0)MDR=DSQRT(MDR)

* Calculation of the SUSY corrections to h_bot, DELMB, as in
* Carena et al., hep-ph/9912516
      
      C2T=1d0-S2T
      C2B=1d0-S2B

      DELMB=MUQ*TANBQ*(2d0/(3d0*PI)*ALSMZ*M3*INTEG(MSB1,MSB2,M3)
     .   +COEF*HTQ**2*AT*INTEG(MST1,MST2,MUQ)
     .   -COEF*G2Q*M2*(CST**2*INTEG(MST1,M2,MUQ)
     .   +(1d0-CST**2)*INTEG(MST2,M2,MUQ)
     .   +1d0/2d0*(CSB**2*INTEG(MSB1,M2,MUQ)
     .   +(1d0-CSB**2)*INTEG(MSB2,M2,MUQ))))

!      WRITE(0,*)"DELMB =",DELMB
!      WRITE(0,*)""

*   Stau masses squared and mixings

      mslL= ML3 + MTAU**2 - (gQ/2d0-g1Q/2d0)*(h2q**2-h1q**2)
      mslR= ME3 + MTAU**2 - g1Q/2d0*(h2q**2-h1q**2)
      Xl= Atau-MUQ*tanbq
      Wl= DSQRT((mslL-mslR)**2 + 4d0*(Xl*MTAU)**2)
      msl1= 0.5d0*(mslL+mslR-Wl)
      msl2= 0.5d0*(mslL+mslR+Wl)

      nen= DSQRT((mslL-msl1)**2 + (Xl*MTAU)**2)
      IF(nen.EQ.0d0)THEN
       IF(mslL.LE.mslR)THEN
        CSL= 1d0
       ELSE
        CSL= 0d0
       ENDIF
      ELSE
       CSL= -MTAU*Xl/nen
      ENDIF

      MSNT= ML3 + gQ/2d0*(h2q**2-h1q**2)

      IF(MSL1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSL1=1d2
!      WRITE(0,*)"MSL1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSL2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSL2=1d2
!      WRITE(0,*)"MSL2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSNT.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSNT=1d2
!      WRITE(0,*)"MSNT^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSL1.GE.0d0)MSL1=DSQRT(MSL1)
      IF(MSL2.GE.0d0)MSL2=DSQRT(MSL2)
      IF(MSNT.GE.0d0)MSNT=DSQRT(MSNT)

*   Smuon masses squared and mixings

      mslL= ML + MMUON**2 - (gQ/2d0-g1Q/2d0)*(h2q**2-h1q**2)
      mslR= ME + MMUON**2 - g1Q/2d0*(h2q**2-h1q**2)
      Xl= AMUON-MUQ*tanbq
      Wl= DSQRT((mslL-mslR)**2 + 4d0*(Xl*MMUON)**2)
      MSMU1= 0.5d0*(mslL+mslR-Wl)
      MSMU2= 0.5d0*(mslL+mslR+Wl)

      nen= DSQRT((mslL-MSMU1)**2 + (Xl*MMUON)**2)
      IF(nen.EQ.0d0)THEN
       IF(mslL.LE.mslR)THEN
        CSMU= 1d0
       ELSE
        CSMU= 0d0
       ENDIF
      ELSE
       CSMU= -MMUON*Xl/nen
      ENDIF

      MSMUNT= ML + gQ/2d0*(h2q**2-h1q**2)

      IF(MSMU1.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSMU1=1d2
!      WRITE(0,*)"MSMU1^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSMU2.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSMU2=1d2
!      WRITE(0,*)"MSMU2^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSMUNT.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MSMUNT=1d2
!      WRITE(0,*)"MSMUNT^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MSMU1.GE.0d0)MSMU1=DSQRT(MSMU1)
      IF(MSMU1.GE.0d0)MSMU2=DSQRT(MSMU2)
      IF(MSMUNT.GE.0d0)MSMUNT=DSQRT(MSMUNT)

*   Selectron masses squared

      MLR=ME - g1Q/2d0*(h2q**2-h1q**2)
      MLL=ML + (-gQ/2d0+g1Q/2d0)*(h2q**2-h1q**2)
      MNL=ML + gQ/2d0*(h2q**2-h1q**2)

      IF(MLR.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MLR=1d2
!      WRITE(0,*)"MLR^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MLL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MLL=1d2
!      WRITE(0,*)"MLL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MNL.LE.1d2)THEN
       IFAIL=8
       IF(FLAG.EQ.0)MNL=1d2
!      WRITE(0,*)"MNL^2 < 0"
!      WRITE(0,*)""
!      WRITE(0,*)""
      ENDIF

      IF(MLR.GE.0d0)MLR=DSQRT(MLR)
      IF(MLL.GE.0d0)MLL=DSQRT(MLL)
      IF(MNL.GE.0d0)MNL=DSQRT(MNL)

      END


      DOUBLE PRECISION FUNCTION DMSQUARK(i,msq1,msq2,s2t,mq,mglu,q2)

*    with thanks to S. Kraml
*    msq1, msq2 are the squark masses squared

      implicit none
      double precision msq1,msq2,s2t,mq,mglu,msq,msqp
      double precision mq2,mglu2,sgn,dmg,dmsg,dmsq,q2
      double precision NMA0,NMB0,NMB1
      integer i

      mglu=dabs(mglu)
      mq2 = mq*mq
      mglu2 = mglu*mglu
      if (i.eq.1) then
        msq=msq1
        msqp=msq2
        sgn = -1d0
      else
        msq=msq2
        msqp=msq1
        sgn = 1d0
      endif

      dmg = -2d0*msq*( 2d0*NMB0(msq,0d0,msq,q2)-NMB1(msq,0d0,msq,q2) )
      dmsg = -4d0*( NMA0(mq2,q2) - msq*NMB1(msq,mglu2,mq2,q2) +
     .    (mglu2 + sgn*mglu*mq*s2t)*NMB0(msq,mglu2,mq2,q2) )
      dmsq = (1d0-s2t**2)*NMA0(msq,q2) + s2t**2*NMA0(msqp,q2)
      DMSQUARK = dmg + dmsg + dmsq

      end


      DOUBLE PRECISION FUNCTION INTEG(X,Y,Z)

* Function for DELMB:

      IMPLICIT NONE

      DOUBLE PRECISION X,Y,Z

       IF(DABS(X).EQ.DABS(Y) .AND. DABS(X).EQ.DABS(Z))THEN
        INTEG=.5d0/X
      ELSEIF(DABS(X).EQ.DABS(Y))THEN
        INTEG=(X**2-Z**2+Z**2*DLOG(Z**2/X**2))/(X**2-Z**2)**2
      ELSEIF(DABS(Y).EQ.DABS(Z))THEN
        INTEG=(Y**2-X**2+X**2*DLOG(X**2/Y**2))/(Y**2-X**2)**2
      ELSEIF(DABS(X).EQ.DABS(Z))THEN
        INTEG=(X**2-Y**2+Y**2*DLOG(Y**2/X**2))/(X**2-Y**2)**2
      ELSE
        INTEG=(X**2*Y**2*DLOG(X**2/Y**2)
     .   +Y**2*Z**2*DLOG(Y**2/Z**2)+Z**2*X**2*DLOG(Z**2/X**2))/
     .   ((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
      ENDIF

      END
