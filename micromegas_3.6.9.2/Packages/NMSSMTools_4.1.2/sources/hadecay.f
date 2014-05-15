      SUBROUTINE HDECAY(MH,BRJJ,BRMM,BRLL,BRSS,BRCC,
     .        BRBB,BRTT,BRWW,BRZZ,BRGG,BRZG)

      IMPLICIT NONE

      INTEGER N0,NF,NFEXT,NFGG

      DOUBLE PRECISION MH,CJ,CG,PI,HIGTOP,ASG,ASH,AS3,AS4,ASMT
      DOUBLE PRECISION SQR2,EPS,FQCD,XFAC,X,Y,RATCOUP,RAT
      DOUBLE PRECISION HJJ,HMM,HLL,HSS,HCC,HBB,HTT,HWW,HZZ,HGG,HZG
      DOUBLE PRECISION HS1,HS2,HC1,HC2,HB1,HB2,HT1,HT2,DCC,DBB
      DOUBLE PRECISION DLU,DLD,XM1,XM2,CWW,CZZ,XX(4),YY(4)
      DOUBLE PRECISION WIDTH,BRJJ,BRMM,BRLL,BRSS,BRCC
      DOUBLE PRECISION BRBB,BRTT,BRWW,BRZZ,BRGG,BRZG
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0,RMS,RMC,RMB,RMT
      DOUBLE PRECISION C2TW,T2TW,ALEM0,RMTTOP,FT,FB,RUNMB
      DOUBLE PRECISION HVV,HV,HFF,QCD0,HQCDM,HQCD,QCDH,TQCDH,HGGQCD
      DOUBLE PRECISION BETA,SP,ALPHAS,RUNM,QQINT,FINT,ACOUP
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      DOUBLE COMPLEX CTT,CTB,CTC,CTL,CTW,CXT,CXB,CXC,CXL,CXW
      DOUBLE COMPLEX CLT,CLB,CLC,CLW,CXTZ,CXBZ,CXCZ,CXWZ
      DOUBLE COMPLEX CI1,CI2,CGZ,CF,CA,CB

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0

      QQINT(RAT,X,Y)= RAT**2*X+(1d0-RAT**2)*Y
      BETA(X)= DSQRT(1d0-4d0*X)
      CF(CA)= -CDLOG(-(1d0+CDSQRT(1d0-CA))
     .       / (1d0-CDSQRT(1d0-CA)))**2/4d0
      CGZ(CA)= CDSQRT(1d0-CA)/2d0*CDLOG(-(1d0+CDSQRT(1d0-CA))
     .       / (1d0-CDSQRT(1d0-CA)))
      CI1(CA,CB)= CA*CB/2d0/(CA-CB)
     .       + CA**2*CB**2/2/(CA-CB)**2*(CF(CA)-CF(CB))
     .       + CA**2*CB/(CA-CB)**2*(CGZ(CA)-CGZ(CB))
      CI2(CA,CB)= -CA*CB/2d0/(CA-CB)*(CF(CA)-CF(CB))
      HV(X)= 3d0*(1d0-8d0*X+20d0*X**2)/DSQRT((4d0*X-1d0))
     .       * DACOS((3d0*X-1d0)/2d0/DSQRT(X**3))
     .       - (1d0-X)*(47d0/2d0*X-13d0/2d0+1d0/X)
     .       - 3d0/2d0*(1d0-6d0*X+4d0*X**2)*DLOG(X)
      HVV(X,Y)= GF/(4d0*PI*SQR2)*X**3/2d0*BETA(Y)
     .       * (1d0-4d0*Y+12d0*Y**2)
      HFF(X,Y)= GF/(4d0*PI*SQR2)*X**3*Y*(BETA(Y))**3
      QCD0(X)= (1d0+X**2)*(4d0*SP((1d0-X)/(1d0+X))
     .       + 2d0*SP((X-1d0)/(X+1d0))
     .       - 3d0*DLOG((1d0+X)/(1d0-X))*DLOG(2d0/(1d0+X))
     .       - 2d0*DLOG((1d0+X)/(1d0-X))*DLOG(X))
     .       - 3d0*X*DLOG(4d0/(1d0-X**2))-4d0*X*DLOG(X)
      HQCDM(X)= QCD0(X)/X+(3d0+34d0*X**2-13d0*X**4)/16d0/X**3
     .       * DLOG((1d0+X)/(1d0-X))+3d0/8d0/X**2*(7d0*X**2-1d0)
      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
     .   + (29.14671d0 + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
     .   + DLOG(X)**2/9d0))*(ASH/PI)**2
     .   +(164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
      QCDH(X)= 1d0+HQCD(X)
      TQCDH(X)= 1d0+4d0/3d0*HQCDM(BETA(X))*ASH/PI
      HGGQCD(ASG,NF)= 1d0+ASG/PI*(95d0/4d0-NF*7d0/6d0)

      EPS= 1.D-8
      PI= 4d0*DATAN(1d0)
      SQR2= DSQRT(2d0)

*   Number of light flavours included in the gluonic decays
*   Higgs -> gg* -> gqq (see hdecay): NFGG = 3
      NFGG= 3

*   Alpha_EM(0)
      ALEM0= 1d0/137.04d0

*   Weak angle theta_W (S2TW = sin(theta_W)):
      C2TW= 1d0-S2TW
      T2TW= S2TW/C2TW

*   Alpha_s at the top pole mass scales, used for the running
*   Yukawa coupling ht and running quark masses RMT below
*   NOTE: MT = top pole mass
      ASMT= ALPHAS(MT,2)
      
*   MT = Top pole mass; RMTTOP = running mass at Mtop (MS_bar):
      RMTTOP= MT/(1d0+4d0*ASMT/(3d0*PI)+11d0*(ASMT/PI)**2)

      HIGTOP= (MH/MT)**2
      MT0= 3.D8
      ASH= ALPHAS(MH,2)
      MC0= 1.D8
      MB0= 2.D8
      AS3= ALPHAS(MH,2)
      MC0= MC
      AS4= ALPHAS(MH,2)
      MB0= MBP
      MT0= MT

*  Running quark masses at MH

      RMS= RUNM(MH,3)

      RMC= RUNM(MH,4)

      RMB= RUNMB(MH)

      IF(MH.GE.MT)THEN
       RMT= RMTTOP
     .  *(1d0+7d0/(4d0*PI)*ASMT*DLOG(MH**2/MT**2))
     .  **(-4d0/7d0)
      ELSE
        RMT= RMTTOP
     .  *(1d0+23d0/(12d0*PI)*ASMT*DLOG(MH**2/MT**2))
     .  **(-12d0/23d0)
      ENDIF

*  Radiative couplings

      CTT= 4d0*(MT/MH)**2*DCMPLX(1d0,-EPS)
      CTB= 4d0*(MBP/MH)**2*DCMPLX(1d0,-EPS)
      CTC= 4d0*(MC/MH)**2*DCMPLX(1d0,-EPS)
      CTL= 4d0*(MTAU/MH)**2*DCMPLX(1d0,-EPS)
      CTW= 4d0*(MW/MH)**2*DCMPLX(1d0,-EPS)
      CXT= 2d0*CTT*(1d0+(1d0-CTT)*CF(CTT))
      CXB= 2d0*CTB*(1d0+(1d0-CTB)*CF(CTB))
      CXC= 2d0*CTC*(1d0+(1d0-CTC)*CF(CTC))
      CXL= 2d0*CTL*(1d0+(1d0-CTL)*CF(CTL))
      CXW= -(2d0+3d0*CTW+3d0*CTW*(2d0-CTW)*CF(CTW))
      CJ= CDABS(CXT+CXC+CXB)
      CG= CDABS(4d0/3d0*(CXT+CXC)+CXB/3d0+CXL+CXW)

*  Partial widths

*   h -> gg

      NFEXT= 3
      ASG= AS3
      FQCD= HGGQCD(ASG,NFEXT)
      XFAC= CJ**2*FQCD
      HJJ= HVV(MH,0d0)*(ASG/PI)**2*XFAC/8d0

*   h -> gg* -> gcc to be added to h -> cc

      NFEXT= 4
      ASG= AS4
      FQCD= HGGQCD(ASG,NFEXT)
      XFAC= CJ**2*FQCD
      DCC= HVV(MH,0d0)*(ASG/PI)**2*XFAC/8d0-HJJ

*   h -> gg* -> gbb to be added to h -> bb

      NFEXT= 5
      ASG= ASH
      FQCD= HGGQCD(ASG,NFEXT)
      XFAC= CJ**2*FQCD
      DBB= HVV(MH,0d0)*(ASG/PI)**2*XFAC/8d0-HJJ-DCC

      IF(NFGG.EQ.5)THEN
       HJJ= HJJ+DBB+DCC
       DBB= 0d0
       DCC= 0d0
      ELSEIF(NFGG.EQ.4)THEN
       HJJ= HJJ+DCC
       DCC= 0d0
      ENDIF

*   h -> mumu

      IF(MH.LE.2d0*MMUON)THEN
       HMM= 0d0
      ELSE
       HMM= HFF(MH,(MMUON/MH)**2)
      ENDIF

*   h -> tautau

      IF(MH.LE.2d0*MTAU)THEN
       HLL= 0d0
      ELSE
       HLL= HFF(MH,(MTAU/MH)**2)
      ENDIF

*   h -> ss

      RATCOUP= 1d0
      IF(MH.LE.2d0*MS)THEN
       HSS= 0d0
      ELSE
       HS1= 3d0*HFF(MH,(MS/MH)**2)
     .    * TQCDH((MS/MH)**2)
       HS2= 3d0*HFF(MH,(RMS/MH)**2)
     .    * QCDH((RMS/MH)**2)
       IF(HS2.LT.0d0) HS2=0d0
       RAT= 2d0*MS/MH
       HSS= QQINT(RAT,HS1,HS2)
      ENDIF

*   h -> cc

      RATCOUP= 1d0
      IF(MH.LE.2d0*MC)THEN
       HCC= 0d0
      ELSE
       HC1= 3d0*HFF(MH,(MC/MH)**2)
     .    * TQCDH((MC/MH)**2)
       HC2= 3d0*HFF(MH,(RMC/MH)**2)
     .   * QCDH((RMC/MH)**2)
     .   + DCC
       IF(HC2.LT.0d0) HC2=0d0
       RAT= 2d0*MC/MH
       HCC= QQINT(RAT,HC1,HC2)
      ENDIF

*   h -> bb

      RATCOUP= 1d0
      IF(MH.LE.2d0*MBP)THEN
       HBB= 0d0
      ELSE
       HB1= 3d0*HFF(MH,(MBP/MH)**2)
     .    * TQCDH((MBP/MH)**2)
       HB2= 3d0*HFF(MH,(RMB/MH)**2)
     .    * QCDH((RMB/MH)**2)
     .    + DBB
       IF(HB2.LT.0d0) HB2=0d0
       RAT= 2d0*MBP/MH
       HBB= QQINT(RAT,HB1,HB2)
      ENDIF

*   h -> tt

      RATCOUP= 0d0
      IF (MH.LE.2d0*MT)THEN
       HTT= 0d0
      ELSE
       RMT= RUNM(MH,6)
       HT1= 3d0*HFF(MH,(MT/MH)**2)
     .    * TQCDH((MT/MH)**2)
       HT2= 3d0*HFF(MH,(RMT/MH)**2)
     .    * QCDH((RMT/MH)**2)
       IF(HT2.LT.0d0) HT2=0d0
       RAT= 2d0*MT/MH
       HTT= QQINT(RAT,HT1,HT2)
      ENDIF

*   h -> WW

      DLD= 2d0
      DLU= 2d0
      XM1= 2d0*MW-DLD
      XM2= 2d0*MW+DLU
      IF(MH.LE.MW)THEN
       HWW= 0d0
      ELSEIF(MH.LE.XM1)THEN
       CWW= 3d0*GF**2*MW**4/16d0/PI**3
       HWW= HV((MW/MH)**2)*CWW*MH
      ELSEIF(MH.LT.XM2)THEN
       CWW= 3d0*GF**2*MW**4/16d0/PI**3
       XX(1)= XM1-1d0
       XX(2)= XM1
       XX(3)= XM2
       XX(4)= XM2+1d0
       YY(1)= HV((MW/XX(1))**2)*CWW*XX(1)
       YY(2)= HV((MW/XX(2))**2)*CWW*XX(2)
       YY(3)= HVV(XX(3),(MW/XX(3))**2)
       YY(4)= HVV(XX(4),(MW/XX(4))**2)
       HWW= FINT(MH,XX,YY)
      ELSE
       HWW= HVV(MH,(MW/MH)**2)
      ENDIF

*   h -> ZZ

      DLD= 2d0
      DLU= 2d0
      XM1= 2d0*MZ-DLD
      XM2= 2d0*MZ+DLU
      IF(MH.LE.MZ)THEN
       HZZ= 0d0
      ELSEIF(MH.LE.XM1)THEN
       CZZ= 3d0*GF**2*MZ**4/192d0/PI**3
     .    * (7d0-40d0/3d0*S2TW+160d0/9d0*S2TW**2)
       HZZ= HV((MZ/MH)**2)*CZZ*MH
      ELSEIF(MH.LT.XM2)THEN
       CZZ= 3d0*GF**2*MZ**4/192d0/PI**3
     .    * (7d0-40d0/3d0*S2TW+160d0/9d0*S2TW**2)
       XX(1)= XM1-1d0
       XX(2)= XM1
       XX(3)= XM2
       XX(4)= XM2+1d0
       YY(1)= HV((MZ/XX(1))**2)*CZZ*XX(1)
       YY(2)= HV((MZ/XX(2))**2)*CZZ*XX(2)
       YY(3)= HVV(XX(3),(MZ/XX(3))**2)/2d0
       YY(4)= HVV(XX(4),(MZ/XX(4))**2)/2d0
       HZZ= FINT(MH,XX,YY)
      ELSE
       HZZ= HVV(MH,(MZ/MH)**2)/2d0
      ENDIF

*   h -> gamma gamma

      XFAC= CG**2
      HGG= HVV(MH,0d0)*(ALEM0/PI)**2/16d0*XFAC

*  h -> Z gamma

      IF(MH.LE.MZ)THEN
       HZG= 0d0
      ELSE
       FT= -2d0*(1d0-8d0/3d0*S2TW)/DSQRT(S2TW*C2TW)
       FB= (-1d0+4d0/3d0*S2TW)/DSQRT(S2TW*C2TW)
       CLT= 4d0*(MT/MZ)**2*DCMPLX(1d0,-EPS)
       CLB= 4d0*(MBP/MZ)**2*DCMPLX(1d0,-EPS)
       CLC= 4d0*(MC/MZ)**2*DCMPLX(1d0,-EPS)
       CLW= 4d0*(MW/MZ)**2*DCMPLX(1d0,-EPS)
       CXTZ= FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CXBZ= FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CXCZ= FT*(CI1(CTC,CLC) - CI2(CTC,CLC))
       CXWZ= -1d0/DSQRT(T2TW)*(4d0*(3d0-T2TW)*CI2(CTW,CLW)
     .     + ((1d0+2d0/CTW)*T2TW - (5d0+2d0/CTW))*CI1(CTW,CLW))
       XFAC= CDABS(CXTZ+CXBZ+CXCZ+CXWZ)**2
       ACOUP= SQR2*GF*MZ**2*S2TW*C2TW/PI**2
       HZG= GF/(4d0*PI*SQR2)*MH**3*(ALEM0/PI)*ACOUP/16d0
     .    * XFAC*(1d0-(MZ/MH)**2)**3
      ENDIF

*  Branching ratios

      WIDTH=HJJ+HMM+HLL+HSS+HCC+HBB+HTT+HWW+HZZ+HGG+HZG
      BRJJ= HJJ/WIDTH
      BRMM= HMM/WIDTH
      BRLL= HLL/WIDTH
      BRSS= HSS/WIDTH
      BRCC= HCC/WIDTH
      BRBB= HBB/WIDTH
      BRTT= HTT/WIDTH
      BRWW= HWW/WIDTH
      BRZZ= HZZ/WIDTH
      BRGG= HGG/WIDTH
      BRZG= HZG/WIDTH

      END


      SUBROUTINE ADECAY(MA,BRJJ,BRMM,BRLL,BRSS,BRCC,
     .        BRBB,BRTT,BRGG,BRZG)

      IMPLICIT NONE

      INTEGER N0,NF,NFEXT,NFGG

      DOUBLE PRECISION MA,CJ,CG,PI,HIGTOP,ASG,ASH,AS3,AS4,ASMT
      DOUBLE PRECISION SQR2,EPS,FQCD,XFAC,X,Y,RATCOUP,RAT
      DOUBLE PRECISION HJJ,HMM,HLL,HSS,HCC,HBB,HTT,HGG,HZG
      DOUBLE PRECISION HS1,HS2,HC1,HC2,HB1,HB2,HT1,HT2,DCC,DBB
      DOUBLE PRECISION WIDTH,BRJJ,BRMM,BRLL,BRSS,BRCC
      DOUBLE PRECISION BRBB,BRTT,BRGG,BRZG
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0,RMS,RMC,RMB,RMT
      DOUBLE PRECISION C2TW,T2TW,ALEM0,RMTTOP,FT,FB,RUNMB
      DOUBLE PRECISION AFF,QCD0,AQCDM,AQCD,QCDA,TQCDA,AGGQCD
      DOUBLE PRECISION BETA,SP,ALPHAS,RUNM,QQINT,ACOUP
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      DOUBLE COMPLEX CTT,CTB,CTC,CTL,CXT,CXB,CXC,CXL
      DOUBLE COMPLEX CLT,CLB,CLC
      DOUBLE COMPLEX CI1,CI2,CGZ,CF,CA,CB

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0

      QQINT(RAT,X,Y)= RAT**2*X+(1d0-RAT**2)*Y
      BETA(X)= DSQRT(1d0-4d0*X)
      CF(CA)= -CDLOG(-(1d0+CDSQRT(1d0-CA))
     .       / (1d0-CDSQRT(1d0-CA)))**2/4d0
      CGZ(CA)= CDSQRT(1d0-CA)/2d0*CDLOG(-(1d0+CDSQRT(1d0-CA))
     .       / (1d0-CDSQRT(1d0-CA)))
      CI1(CA,CB)= CA*CB/2d0/(CA-CB)
     .       + CA**2*CB**2/2/(CA-CB)**2*(CF(CA)-CF(CB))
     .       + CA**2*CB/(CA-CB)**2*(CGZ(CA)-CGZ(CB))
      CI2(CA,CB)= -CA*CB/2d0/(CA-CB)*(CF(CA)-CF(CB))
      AFF(X,Y)= GF/(4d0*PI*SQR2)*X**3*Y*(BETA(Y))
      QCD0(X)= (1d0+X**2)*(4d0*SP((1d0-X)/(1d0+X))
     .       + 2d0*SP((X-1d0)/(X+1d0))
     .       - 3d0*DLOG((1d0+X)/(1d0-X))*DLOG(2d0/(1d0+X))
     .       - 2d0*DLOG((1d0+X)/(1d0-X))*DLOG(X))
     .       - 3d0*X*DLOG(4d0/(1d0-X**2))-4d0*X*DLOG(X)
      AQCDM(X)= QCD0(X)/X+(19d0+2d0*X**2+3d0*X**4)/16d0/X
     .       * DLOG((1d0+X)/(1d0-X))+3d0/8d0*(7d0-X**2)
      AQCD(X)=(4d0/3d0*AQCDM(BETA(X))
     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-6d0*X)/(1d0-4d0*X))*ASH/PI
     .   + (29.14671d0 + RATCOUP*(23d0/6d0 - DLOG(HIGTOP)
     .   + DLOG(X)**2/6d0))*(ASH/PI)**2
     .   + (164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
      QCDA(X)= 1d0+AQCD(X)
      TQCDA(X)= 1d0+4d0/3d0*AQCDM(BETA(X))*ASH/PI
      AGGQCD(ASG,NF)= 1d0+ASG/PI*(97d0/4d0-NF*7d0/6d0)

      EPS= 1.D-8
      PI= 4d0*DATAN(1d0)
      SQR2= DSQRT(2d0)

*   Number of light flavours included in the gluonic decays
*   Higgs -> gg* -> gqq (see hdecay): NFGG = 3
      NFGG= 3

*   Alpha_EM(0)
      ALEM0= 1d0/137.04d0

*   Weak angle theta_W (S2TW = sin(theta_W)):
      C2TW= 1d0-S2TW
      T2TW= S2TW/C2TW

*   Alpha_s at the top pole mass scales, used for the running
*   Yukawa coupling ht and running quark masses RMT below
*   NOTE: MT = top pole mass
      ASMT= ALPHAS(MT,2)
      
*   MT = Top pole mass; RMTTOP = running mass at Mtop (MS_bar):
      RMTTOP= MT/(1d0+4d0*ASMT/(3d0*PI)+11d0*(ASMT/PI)**2)

      HIGTOP= (MA/MT)**2
      MT0= 3.D8
      ASH= ALPHAS(MA,2)
      MC0= 1.D8
      MB0= 2.D8
      AS3= ALPHAS(MA,2)
      MC0= MC
      AS4= ALPHAS(MA,2)
      MB0= MBP
      MT0= MT

*  Running quark masses at MA

      RMS= RUNM(MA,3)

      RMC= RUNM(MA,4)

      RMB= RUNMB(MA)

      IF(MA.GE.MT)THEN
       RMT= RMTTOP
     .  *(1d0+7d0/(4d0*PI)*ASMT*DLOG(MA**2/MT**2))
     .  **(-4d0/7d0)
      ELSE
        RMT= RMTTOP
     .  *(1d0+23d0/(12d0*PI)*ASMT*DLOG(MA**2/MT**2))
     .  **(-12d0/23d0)
      ENDIF

*  Radiative couplings

      CTT= 4d0*(MT/MA)**2*DCMPLX(1d0,-EPS)
      CTB= 4d0*(MBP/MA)**2*DCMPLX(1d0,-EPS)
      CTC= 4d0*(MC/MA)**2*DCMPLX(1d0,-EPS)
      CTL= 4d0*(MTAU/MA)**2*DCMPLX(1d0,-EPS)
      CXT= CTT*CF(CTT)
      CXB= CTB*CF(CTB)
      CXC= CTC*CF(CTC)
      CXL= CTL*CF(CTL)
      CJ= CDABS(CXT+CXC+CXB)
      CG= CDABS(4d0/3d0*(CXT+CXC)+CXB/3d0+CXL)

*  Partial widths

*   a -> gg

      NFEXT= 3
      ASG= AS3
      FQCD= AGGQCD(ASG,NFEXT)
      XFAC= CJ**2*FQCD
      HJJ= GF/(16d0*PI*SQR2)*MA**3*(ASG/PI)**2*XFAC

*   a -> gg* -> gcc to be added to a -> cc

      NFEXT= 4
      ASG= AS4
      FQCD= AGGQCD(ASG,NFEXT)
      XFAC= CJ**2*FQCD
      DCC= GF/(16d0*PI*SQR2)*MA**3*(ASG/PI)**2*XFAC-HJJ

*   a -> gg* -> gbb to be added to a -> bb

      NFEXT= 5
      ASG= ASH
      FQCD= AGGQCD(ASG,NFEXT)
      XFAC= CJ**2*FQCD
      DBB= GF/(16d0*PI*SQR2)*MA**3*(ASG/PI)**2*XFAC-HJJ-DCC

      IF(NFGG.EQ.5)THEN
       HJJ= HJJ+DBB+DCC
       DBB= 0d0
       DCC= 0d0
      ELSEIF(NFGG.EQ.4)THEN
       HJJ= HJJ+DCC
       DCC= 0d0
      ENDIF

*   a -> mumu

      IF(MA.LE.2*MMUON)THEN
       HMM= 0d0
      ELSE
       HMM= AFF(MA,(MMUON/MA)**2)
      ENDIF

*   a -> tautau

      IF(MA.LE.2*MTAU)THEN
       HLL= 0d0
      ELSE
       HLL= AFF(MA,(MTAU/MA)**2)
      ENDIF

*   a -> ss

      RATCOUP= 1d0
      IF(MA.LE.2d0*MS)THEN
       HSS= 0d0
      ELSE
       HS1= 3d0*AFF(MA,(MS/MA)**2)
     .    * TQCDA((MS/MA)**2)
       HS2= 3d0*AFF(MA,(RMS/MA)**2)
     .    * QCDA((RMS/MA)**2)
       IF(HS2.LT.0d0) HS2=0d0
       RAT= 2d0*MS/MA
       HSS= QQINT(RAT,HS1,HS2)
      ENDIF

*   a -> cc

      RATCOUP= 1d0
      IF(MA.LE.2*MC)THEN
       HCC= 0d0
      ELSE
       HC1= 3d0*AFF(MA,(MC/MA)**2)
     .    * TQCDA((MC/MA)**2)
       HC2= 3d0*AFF(MA,(RMC/MA)**2)
     .    * QCDA((RMC/MA)**2)
     .     + DCC
       IF(HC2.LT.0d0) HC2=0d0
       RAT= 2d0*MC/MA
       HCC= QQINT(RAT,HC1,HC2)
      ENDIF

*   a -> bb

      RATCOUP= 1d0
      IF(MA.LE.2*MBP) THEN
       HBB= 0d0
      ELSE
       HB1= 3d0*AFF(MA,(MBP/MA)**2)
     .    * TQCDA((MBP/MA)**2)
       HB2= 3d0*AFF(MA,(RMB/MA)**2)
     .    * QCDA((RMB/MA)**2)
     .     + DBB
       IF(HB2.LT.0d0) HB2=0d0
       RAT= 2d0*MBP/MA
       HBB= QQINT(RAT,HB1,HB2)
      ENDIF

*   a -> tt

      RATCOUP= 0d0
      IF (MA.LE.2d0*MT) THEN
       HTT= 0d0
      ELSE
       HT1= 3d0*AFF(MA,(MT/MA)**2)
     .    * TQCDA((MT/MA)**2)
       HT2= 3d0*AFF(MA,(RMT/MA)**2)
     .    * QCDA((RMT/MA)**2)
       IF(HT2.LT.0d0) HT2=0d0
       RAT= 2d0*MT/MA
       HTT= QQINT(RAT,HT1,HT2)
      ENDIF

*   a -> gamma gamma

      XFAC= CG**2
      HGG= GF/(32d0*PI*SQR2)*MA**3*(ALEM0/PI)**2*XFAC

*   a -> Z gamma

      IF(MA.LE.MZ)THEN
       HZG= 0d0
      ELSE
       FT= -2d0*(1d0-8d0/3d0*S2TW)/DSQRT(S2TW*C2TW)
       FB= (-1d0+4d0/3d0*S2TW)/DSQRT(S2TW*C2TW)
       CLT= 4d0*(MT/MZ)**2*DCMPLX(1d0,-EPS)
       CLB= 4d0*(MBP/MZ)**2*DCMPLX(1d0,-EPS)
       CLC= 4d0*(MC/MZ)**2*DCMPLX(1d0,-EPS)
       CXT= FT*(-CI2(CTT,CLT))
       CXB= FB*(-CI2(CTB,CLB))
       CXC= FT*(-CI2(CTC,CLC))
       XFAC= CDABS(CXT+CXB+CXC)**2
       ACOUP= SQR2*GF*MZ**2*S2TW*C2TW/PI**2
       HZG= GF/(4d0*PI*SQR2)*MA**3*(ALEM0/PI)*ACOUP/16d0
     .    * XFAC*(1d0-(MZ/MA)**2)**3
      ENDIF

*   Branching ratios

       WIDTH=HJJ+HMM+HLL+HSS+HCC+HBB+HTT+HGG+HZG
       BRJJ= HJJ/WIDTH
       BRMM= HMM/WIDTH
       BRLL= HLL/WIDTH
       BRSS= HSS/WIDTH
       BRCC= HCC/WIDTH
       BRBB= HBB/WIDTH
       BRTT= HTT/WIDTH
       BRGG= HGG/WIDTH
       BRZG= HZG/WIDTH

      END
