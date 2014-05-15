      SUBROUTINE RUNPAR(PAR)

*******************************************************************
* Subroutine to evolve the parameters from the scale Q2 to QSTSB
*
* The mass scale QSTSB of the 3rd generation squarks is
* computed here as QSTSB = MQ3*MU3 (QSTSB has the dimension mass^2)
*
* Several parameters at the scale QSTSB (relevant for the
* calculation of the Higgs masses) are computed here and stored in
* COMMON/QGAUGE, /QHIGGS, /QQUARK, /QPAR, /QEXT:
* - The electroweak and strong gauge couplings (the electroweak
*   couplings appear in the tree level Higgs mass matrix in MHIGGS,
*   alpha_s only in the two loop Higgs mass corrections)
* - the Higgs wave function normalization constants ZHU, ZHD, ZS
*   and the Higgs vevs H1Q, H2Q and TANBQ
* - the top/bottom Yukawa couplings HTQ/HBQ and masses MTOPQ/MBOTQ
* - the NMSSM parameters LQ, KQ, MUQ, NUQ, ALQ, AKQ, XIFQ, XISQ,
*   MUPQ, MSPQ, M3HQ at QSTSB
* - if MOD(MAFLAG,3) = 0: AL at Q2 is given in PAR(5), XIF at Q2 is
*   given, MA at QSTSB is computed here and stored in PAR(23)
* - if MOD(MAFLAG,3) = 1: MA at QSTSB is given in PAR(23), XIF at Q2
*   is given, AL at Q2 is computed here and stored in PAR(5)
* - if MOD(MAFLAG,3) = 2: MA at QSTSB is given in PAR(23), AL at Q2
*   is given in PAR(5), XIF at Q2 is computed here
* - if INT(MAFLAG/3) = 0: AK at Q2 is given in PAR(6), XIS at Q2 is
*   given, MP at QSTSB is computed here and stored in PAR(24)
* - if INT(MAFLAG/3) = 1: MP at QSTSB is given in PAR(24), XIS at Q2
*   is given, AK at Q2 is computed here and stored in PAR(6)
* - if INT(MAFLAG/3) = 2: MP at QSTSB is given in PAR(24), AK at Q2
*   is given in PAR(6), XIS at Q2 is computed here
* - if MAFLAG < 0: this is the case NMSPEC and NMGMSB
*   then MUQ, (KQ+NUQ or XIFQ), (MSQ or XISQ), MA and MP are
*   subsequently recomputed in the subroutine MINIMIZE
*
* The SUSY scale Q2, where the soft terms are
* defined on input, is possibly much larger than QSTSB.
* Unless Q2 is defined by the user, Q2 is defined here in terms
* of the first generation squark masses as
* MSUSY**2 == Q2 = MAX((2*MQ**2+MU**2+MD**2)/4,Q2MIN).
*
*******************************************************************

      IMPLICIT NONE

      INTEGER Q2FIX,OMGFLAG,MAFLAG

      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION PI,COEF
      DOUBLE PRECISION tanb,SB2,CB2,h1,h2,M1,M2,HTAU
      DOUBLE PRECISION L,K,AL,AK,MU,NU,RUNMB,LQT,ALSMT,HT,HB      
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION Q2MIN,Q2,QSTSB,MA2,MP2
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION Lmu,Lnu,LM1mu,LM2MU,Lmunu,LQ2,LMAMT
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION ALSMA,DLA,DLQA,F1,F2,HTMA

      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/FLAGS/OMGFLAG,MAFLAG

      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

!      WRITE(0,*)"CALL RUNPAR"
!      WRITE(0,*)""

*   Definition of the SUSY scale Q2, unless defined initially
      IF(Q2FIX.EQ.0)THEN
       Q2=MAX((2d0*PAR(15)+PAR(16)+PAR(17))/4d0,Q2MIN)
      ENDIF
!      WRITE(0,*)"QSUSY =",DSQRT(Q2)

*   Definition of the scale QSTSB
      QSTSB=DSQRT(MAX(PAR(7)*PAR(8),Q2MIN**2))
*      QSTSB=MAX((2d0*PAR(7)+PAR(8)+PAR(9))/4,Q2MIN)
      LQ2=DLOG(QSTSB/Q2)
      LQT=DLOG(MAX(QSTSB,MT**2)/MT**2)
!      WRITE(0,*)"QSTSB =",DSQRT(QSTSB)

*   NMSSM parameters at Q2
      L=PAR(1)
      K=PAR(2)
      tanb=PAR(3)      
      SB2=tanb**2/(1d0+tanb**2)
      CB2=1d0-SB2
      h1=DSQRT(SB2/(2d0*DSQRT(2d0)*GF))
      h2=h1/tanb
      MU=PAR(4)
      NU=K/L*MU
      M1=PAR(20)
      M2=PAR(21)

      IF(MAFLAG.LT.0)THEN
       AL=PAR(5)
       AK=PAR(6)
       MA2=PAR(23)**2
      ELSE
       IF(MOD(MAFLAG,3).EQ.0)THEN
        AL=PAR(5)
        MA2=MAX(((AL+NU)*MU+M3H+MU*MUP+L*XIF)*(tanb+1d0/tanb),1d0)
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        MA2=PAR(23)**2
        AL=(MA2*TANB/(1d0+TANB**2)-M3H-L*XIF)/MU-MUP-NU
       ELSE
        AL=PAR(5)
        MA2=PAR(23)**2
        XIF=(MA2*TANB/(1d0+TANB**2)-(AL+NU)*MU-M3H-MU*MUP)/L
       ENDIF
       IF(MAFLAG/3.EQ.0)THEN
        AK=PAR(6)
       ELSEIF(MAFLAG/3.EQ.1)THEN
        MP2=PAR(24)**2
        IF(K.EQ.0d0)THEN
         AK=0d0
        ELSE
         AK=(L**2*(AL+4d0*NU+MUP)*H1*H2/MU-L*(XIF*MUP+XIS)/MU
     .      -MUP*NU-4d0*K*XIF-2d0*MSP-MP2)/(3d0*NU)
        ENDIF
       ELSE
        MP2=PAR(24)**2
        AK=PAR(6)
        XIS=L*(AL+4d0*NU+MUP)*H1*H2-(3d0*AK*NU+MUP*NU
     .      +4d0*K*XIF+2d0*MSP+MP2)*MU/L-XIF*MUP
       ENDIF
      ENDIF

*   Electroweak gauge couplings at QSTSB:

*   g_2**2 including the Higgs and sparticle thresholds:
      g2q=g2/(1d0+g2*COEF*(DLOG(QSTSB/MZ**2)*19d0/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2d0/3d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/2d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2)))
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3d0
     .    -DLOG(QSTSB/MIN(QSTSB,MAX(M2**2,MZ**2)))*4d0/3d0))

*   g_1**2 including the top, Higgs and sparticle thresholds:
      g1q=g1/(1d0-g1*COEF*(DLOG(QSTSB/MZ**2)*53d0/9d0
     .    +DLOG(QSTSB/MT**2)*17d0/18d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(MA2,MZ**2)))/6d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(MU**2,MZ**2)))*2d0/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(7),MZ**2)))/18d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(8),MZ**2)))*4d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(9),MZ**2)))/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(10),MZ**2)))/6d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(11),MZ**2)))/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(15),MZ**2)))/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(16),MZ**2)))*8d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(17),MZ**2)))*2d0/9d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(18),MZ**2)))/3d0
     .    +DLOG(QSTSB/MIN(QSTSB,MAX(PAR(19),MZ**2)))*2d0/3d0))

      gq=(g1q+g2q)/2d0

*   Alphas at MT and QSTSB
      ALSMT=ALSMZ/(1d0+23d0/(12d0*PI)*ALSMZ*DLOG(MT**2/MZ**2))
      ALSQ=ALSMT/(1d0+ALSMT/(4d0*PI)*(7d0*LQT-2d0*
     .    DLOG(MAX(QSTSB,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))

*   Yukawas at MT (ht, hb: running MS_bar couplings)
      HT=MT/(1d0+4d0*ALSMT/(3d0*PI)+11d0*(ALSMT/PI)**2)/H1
      HB=RUNMB(MT)/H2
      HTAU=MTAU/H2

*   Logs for the Wave Function Renormalization Constants
      Lmu = DLOG(MIN(MAX(mu**2,MZ**2),QSTSB)/QSTSB)
      Lnu = DLOG(MIN(MAX(4d0*nu**2,MZ**2),QSTSB)/QSTSB)
      LM1mu = DLOG(MIN(MAX(M1**2,mu**2,MZ**2),QSTSB)/QSTSB)
      LM2mu = DLOG(MIN(MAX(M2**2,mu**2,MZ**2),QSTSB)/QSTSB)
      Lmunu = DLOG(MIN(MAX(mu**2,4d0*nu**2,MZ**2),QSTSB)/QSTSB)
      LMAMT = DLOG(MIN(MAX(MA2,MT**2),QSTSB)/MT**2)

* Aux. quantities for the resummation of logs ~ht^2*LQT
* and ~ht^2*LMAMT:

      ALSMA=ALSMT/(1d0+ALSMT/(4d0*PI)*(7d0*LMAMT-2d0*
     .    DLOG(MAX(MA2,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))
      DLA=(ALSMA/ALSMT)**(1d0/7d0)
      DLQA=(ALSQ/ALSMA)**(1d0/7d0)
      F1=DABS(1d0-9d0*SB2*HT**2*(1d0-DLA)/(8d0*PI*ALSMT))
      HTMA=HT*DLA**4/DSQRT(F1)
      F2=DABS(1d0-9d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))
      
      ZHU=(F1*F2)**(-2d0/3d0)*(1d0+COEF*(
     .    +CB2*(3d0*HB**2+(MTAU/H2)**2)*LMAMT
     .    -G1Q/2d0*LM1mu-3d0*G2Q/2d0*LM2mu
     .    -3d0/4d0*(G1Q+3d0*G2Q)*DLOG(QSTSB/MZ**2)
     .    -L**2*LMUNU))

      ZHD=F1**(-2d0/3d0)*(1d0+COEF*
     .     (3d0*hb**2*LQT+(MTAU/H2)**2*DLOG(QSTSB/MZ**2)
     .    +SB2*(-3d0*HB**2-(MTAU/H2)**2)*LMAMT
     .    -G1Q/2d0*LM1mu-3d0*G2Q/2d0*LM2mu
     .    -3d0/4d0*(G1Q+3d0*G2Q)*DLOG(QSTSB/MZ**2)
     .    -L**2*LMUNU))

      ZS=1d0-2d0*COEF*(L**2*Lmu+K**2*Lnu)

*   Higgs Vevs at QSTSB
      H1Q=H1/DSQRT(ZHU)
      H2Q=H2/DSQRT(ZHD)
      TANBQ=H1Q/H2Q

*   Top/Bottom Yukawas at QSTSB
*   (Note: RUNMB(Q) includes QCD corrections only)
*   including electroweak contributions, Logs ht^2*LQT resummed:
      HTQ=HT*(1d0+7d0/(4d0*PI)*ALSMT*LQT)**(-4d0/7d0)
     .   /DSQRT(F1*F2)
     .   *(1d0+COEF/4d0*((-17d0/6d0*g1q-9d0/2d0*g2q+hb**2)*LQT
     .   +((3d0*CB2-1d0)*HB**2+2d0*HTAU**2*CB2)*LMAMT
     .   -2d0*L**2*Lmunu-G1Q*LM1mu-3d0*G2Q*LM2mu))

      HBQ=RUNMB(DSQRT(QSTSB))/H2Q*F1**(-1d0/6d0)
     .   *(1d0-3d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))**(-1d0/6d0)
     .   *(1d0+COEF/4d0*((-5d0/6d0*g1q-9d0/2d0*g2q
     .   +9d0*hb**2+2d0*HTAU**2)*LQT
     .   +(-9d0*SB2*HB**2-2d0*HTAU**2*SB2)*LMAMT
     .   -2d0*L**2*Lmunu-G1Q*LM1mu-3d0*G2Q*LM2mu))

*   Conversion to DR_bar:
      HTQ=HTQ*(1d0-ALSQ/(3d0*PI)+g2q*COEF*3d0/8d0)
      HBQ=HBQ*(1d0-ALSQ/(3d0*PI)+g2q*COEF*3d0/8d0)
      
*   Running Top and Bottom Quark Masses
      MTOPQ=HTQ*H1Q
      MBOTQ=HBQ*H2Q

*   NMSSM Parameters at QSTSB, stored in COMMON/QPAR and QEXT
      LQ=L*(1d0+COEF/2d0*(-G1Q-3d0*G2Q+4d0*L**2+2d0*K**2
     .  +3d0*HTQ**2+3d0*HBQ**2+HTAU**2)*LQ2)
      KQ=K*(1d0+3d0*COEF*(L**2+K**2)*LQ2)
      MUQ=MU*(1d0+COEF/2d0*(3d0*HTQ**2+3d0*HBQ**2+HTAU**2
     .   +2d0*L**2-G1Q-3d0*G2Q)*LQ2)
      NUQ=MUQ*KQ/LQ
      MUPQ=MUP*(1d0+2d0*COEF*(L**2+K**2)*LQ2)
      MSPQ=MSP+COEF*(2d0*L**2*(MSP+2d0*MUP*AL)
     .    +4d0*K**2*(MSP+2d0*MUP*AK)+4d0*L*K*M3H)*LQ2
      M3HQ=M3H+COEF/2d0*((3d0*HTQ**2+3d0*HBQ**2+HTAU**2+6d0*L**2
     .    -G1Q-3d0*G2Q)*M3H+2d0*L*K*MSP)*LQ2

      IF(MAFLAG.LT.0)THEN
       ALQ=AL+COEF*(4d0*L**2*AL+2d0*K**2*AK+3d0*HTQ**2*PAR(12)
     .    +3d0*HBQ**2*PAR(13)+HTAU**2*PAR(14)+G1Q*M1+3d0*G2Q*M2)*LQ2
       AKQ=AK+6d0*COEF*(L**2*AL+K**2*AK)*LQ2
       XIFQ=XIF*(1d0+COEF*(L**2+K**2)*LQ2)
       XISQ=XIS+COEF*(L**2*(XIS+2d0*AL*XIF)+K**2*(XIS+2d0*AK*XIF)
     .    +2d0*L*M3H*(AL+MUP)+K*MSP*(AK+MUP))*LQ2
      ELSE
       IF(MOD(MAFLAG,3).EQ.0)THEN
        ALQ=AL+COEF*(4d0*L**2*AL+2d0*K**2*AK+3d0*HTQ**2*PAR(12)
     .    +3d0*HBQ**2*PAR(13)+HTAU**2*PAR(14)+G1Q*M1+3d0*G2Q*M2)*LQ2
        XIFQ=XIF*(1d0+COEF*(L**2+K**2)*LQ2)
        MA2=((ALQ+NUQ+MUPQ)*MUQ+M3HQ+LQ*XIFQ)*(TANBQ+1d0/TANBQ)
        PAR(23)=DSQRT(MAX(MA2,1d0))
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        XIFQ=XIF*(1d0+COEF*(L**2+K**2)*LQ2)
        ALQ=(MA2*TANBQ/(1d0+TANBQ**2)-M3HQ-LQ*XIFQ)/MUQ-MUPQ-NUQ
        PAR(5)=ALQ-COEF*(4d0*L**2*AL+2d0*K**2*AK+3d0*HTQ**2*PAR(12)
     .    +3d0*HBQ**2*PAR(13)+HTAU**2*PAR(14)+G1Q*M1+3d0*G2Q*M2)*LQ2
       ELSE
        ALQ=AL+COEF*(4d0*L**2*AL+2d0*K**2*AK+3d0*HTQ**2*PAR(12)
     .    +3d0*HBQ**2*PAR(13)+HTAU**2*PAR(14)+G1Q*M1+3d0*G2Q*M2)*LQ2
        XIFQ=(MA2*TANBQ/(1d0+TANBQ**2)-(ALQ+NUQ)*MUQ-M3H-MUQ*MUP)/LQ
        XIF=XIFQ*(1d0-COEF*(L**2+K**2)*LQ2)
       ENDIF
       IF(MAFLAG/3.EQ.0)THEN
        AKQ=AK+6d0*COEF*(L**2*AL+K**2*AK)*LQ2
        XISQ=XIS+COEF*(L**2*(XIS+2d0*AL*XIF)+K**2*(XIS+2d0*AK*XIF)
     .    +2d0*L*M3H*(AL+MUP)+K*MSP*(AK+MUP))*LQ2
        MP2=LQ**2*(ALQ+4d0*NUQ+MUPQ)*H1Q*H2Q/MUQ-3d0*AKQ*NUQ
     .    -LQ*(XIFQ*MUPQ+XISQ)/MUQ-MUPQ*NUQ-4d0*KQ*XIFQ-2d0*MSPQ
        PAR(24)=DSQRT(MAX(MP2,1d0))
       ELSEIF(MAFLAG/3.EQ.1)THEN
        XISQ=XIS+COEF*(L**2*(XIS+2d0*AL*XIF)+K**2*(XIS+2d0*AK*XIF)
     .    +2d0*L*M3H*(AL+MUP)+K*MSP*(AK+MUP))*LQ2
        IF(K.EQ.0d0)THEN
         AKQ=0d0
         PAR(6)=0d0
        ELSE
         AKQ=(LQ**2*(ALQ+4d0*NUQ+MUPQ)*H1Q*H2Q/MUQ
     .      -LQ*(XIFQ*MUPQ+XISQ)/MUQ-MUPQ*NUQ
     .      -4d0*KQ*XIFQ-2d0*MSPQ-MP2)/(3d0*NUQ)
         PAR(6)=AKQ-6d0*COEF*(L**2*AL+K**2*AK)*LQ2
        ENDIF
       ELSE
        AKQ=AK+6d0*COEF*(L**2*AL+K**2*AK)*LQ2
        XISQ=LQ*(ALQ+4d0*NUQ+MUPQ)*H1Q*H2Q-(3d0*AKQ*NUQ+MUPQ*NUQ
     .    +4d0*KQ*XIFQ+2d0*MSPQ+MP2)*MUQ/LQ-XIFQ*MUPQ
        XIS=XISQ-COEF*(L**2*(XIS+2d0*AL*XIF)+K**2*(XIS+2d0*AK*XIF)
     .    +2d0*L*M3H*(AL+MUP)+K*MSP*(AK+MUP))*LQ2
       ENDIF
      ENDIF

!      WRITE(0,*)"LQ =",LQ
!      WRITE(0,*)"KQ =",KQ
!      WRITE(0,*)"MUQ =",MUQ
!      WRITE(0,*)"NUQ =",NUQ
!      WRITE(0,*)"ALQ =",ALQ
!      WRITE(0,*)"AKQ =",AKQ
!      WRITE(0,*)"XIFQ =",XIFQ
!      WRITE(0,*)"XISQ =",XISQ
!      WRITE(0,*)"MUPQ =",MUPQ
!      WRITE(0,*)"MSPQ =",MSPQ
!      WRITE(0,*)"M3HQ =",M3HQ
!      WRITE(0,*)""
!      WRITE(0,*)""

      END
