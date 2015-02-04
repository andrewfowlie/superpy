      PROGRAM MAIN
      
      IMPLICIT NONE

*     May need these for problem point test.
      INTEGER NPROB,NPAR,ERR
      PARAMETER (NPROB=51,NPAR=25)
      DOUBLE PRECISION PAR(NPAR),PROB(NPROB)
     
      INTEGER IFAIL,I
      DOUBLE PRECISION INTEG,DELMB,C2T,C2B,S2T,S2B
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ 
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C     MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C     CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION M3,COEF,HTQ,HBQ,MTQ,MBQ  
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ,AT,M2,PI
      DOUBLE PRECISION MHUS,MHDS,MHSS

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ 
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ !NUQ = kappa * s
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .     MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .     CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ 
      COMMON/DELMB/DELMB        
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/SUSYMH/MHUS,MHDS,MHSS
      
*     I/O files

      OPEN(15,FILE='inp', STATUS= 'UNKNOWN')
      OPEN(17,FILE='spectr', STATUS= 'UNKNOWN') 
      OPEN(18,FILE='decay', STATUS= 'UNKNOWN') 
      OPEN(19,FILE='omega', STATUS= 'UNKNOWN')

*   Initialization

      CALL INITIALIZE()

*     Reading of SLHA input file
      
      CALL SLHAINPUT(PAR,NPAR)

      CALL NONSLHAINPUT
      DO I=1,NPROB
         PROB(I)=0d0
      ENDDO      

      CALL RUNPAR(PAR)

* Calculate DELMB
* Calculation of the SUSY corrections to h_bot, DELMB, as in
* Carena et al., hep-ph/9912516
      M3 = PAR(22)
      M2 = PAR(21)
      AT = PAR(12)
      C2T=1d0-S2T
      C2B=1d0-S2B
     
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

      DELMB=MUQ*TANBQ*(2d0/(3d0*PI)*ALSMZ*M3*INTEG(MSB1,MSB2,M3) !SQCD 
     .   +COEF*HTQ**2*AT*INTEG(MST1,MST2,MUQ) !SEW
     .   -COEF*G2Q*M2*(CST**2*INTEG(MST1,M2,MUQ)
     .   +(1d0-CST**2)*INTEG(MST2,M2,MUQ)
     .   +1d0/2d0*(CSB**2*INTEG(MSB1,M2,MUQ)
     .   +(1d0-CSB**2)*INTEG(MSB2,M2,MUQ))))

*   Computation of Higgs + top branching ratios

      CALL DECAY(PAR)
      CALL TDECAY(PAR)

*   Exp. constraints

      CALL SUBEXP(PAR,PROB)
      
*   b -> s gamma + B physics

      CALL BSG(PAR,PROB)
      
*   Anom. magn. moment of the Muon

      CALL MAGNMU(PAR,PROB)
      
*   Landau Pole?

      CALL RGES(PAR,PROB,IFAIL)
      
*   Get effective coouplings needed for RELDEN
      
      CALL EFFECTIVECOUPLINGS(PAR)

*   Relic density

      CALL RELDEN(PAR,PROB)

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)IFAIL=10
      ENDDO

*   Sparticle decays

      CALL NMSDECAY(PAR)

*   Recording of the results
      
 11   CALL OUTPUT(PAR,PROB,IFAIL)
      
      END


*******************************************************************
*   This subroutine reads an SLHA file to get the
*   mass spectrum, mixing, DR-bar parameters and entries of PAR.
*******************************************************************
  
      SUBROUTINE SLHAINPUT(PAR,NPAR)

*Fill PAR from SLHA file as

*      PAR(1) = lambda
*      PAR(2) = kappa
*      PAR(3) = tan(beta)
*      PAR(4) = mu (effective mu term = lambda*s)
*      PAR(5) = Alambda
*      PAR(6) = Akappa
*      PAR(7) = mQ3**2
*      PAR(8) = mU3**2
*      PAR(9) = mD3**2
*      PAR(10) = mL3**2
*      PAR(11) = mE3**2
*      PAR(12) = AU3
*      PAR(13) = AD3
*      PAR(14) = AE3
*      PAR(15) = mQ2**2
*      PAR(16) = mU2**2
*      PAR(17) = mD2**2
*      PAR(18) = mL2**2
*      PAR(19) = mE2**2
*      PAR(20) = M1
*      PAR(21) = M2
*      PAR(22) = M3
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element)
*      PAR(25) = AE2
*

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120
      INTEGER I,J,NLINE,INL,ICH,IX,IVAL,Q2FIX
      INTEGER N0,NLOOP,NBER,NPAR,ERR,NTOT,GMUFLAG,HFLAG
      INTEGER OMGFLAG,MAFLAG,PFLAG,NMSFLAG,VFLAG
      
      DOUBLE PRECISION PAR(*),VAL,PI
      !Parameters to be filled
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,MA2,Q2MIN
      DOUBLE PRECISION SMASS(3),S(3,3),SCOMP(3,3)
      DOUBLE PRECISION PMASS(2),P(2,3),P2(2,2), CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),
     C     N(5,5),NEU(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C     MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C     CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,
     C     MSMUL,MSMUR
      DOUBLE PRECISION HUQ,HDQ,MTQ,MBQ
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION QSTSB,Q2
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
*      Added for testing - second gensfermion masses
      DOUBLE PRECISION MDL2,MUL2,MUR2,MDR2 
      DOUBLE PRECISION MHUS,MHDS,MHSS
      DOUBLE PRECISION MQ3P,MU3P,MD3P,ATP,ABP
           
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB  
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,P2,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/DELMB/DELMB   !See line 425 in msferm.f
      COMMON/STSBSCALE/QSTSB
      COMMON/QQUARK/HUQ,HDQ,MTQ,MBQ !HUQ,HDQ are DRbar Yukawa couplings.
      COMMON/QEXT/XIF,XIS,MUP,MSP,M3H 
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/VFLAG/VFLAG
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/SUSYMH/MHUS,MHDS,MHSS
      COMMON/RENSCALE/Q2
      COMMON/RADCOR2/MQ3P,MU3P,MD3P,ATP,ABP

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=0d0
      ENDDO    
      DO I=1,3
       SMASS(I)=1d99
      ENDDO
      DO I=1,2
       PMASS(I)=1d99
      ENDDO
      DO I=1,3
       DO J=1,3
        S(I,J)=1d99
      ENDDO
      ENDDO
      DO I=1,2
       DO J=1,3
        P(I,J)=1d99
      ENDDO
      ENDDO
      CMASS=1d99
      MGL=1d99
      MCHA(1)=1d99
      MCHA(2)=1d99
      DO I=1,2
       DO J=1,2
        U(I,J)=1d99
        V(I,J)=1d99
      ENDDO
      ENDDO
      DO I=1,5
       MNEU(I)=1d99
      ENDDO
      DO I=1,5
       DO J=1,5
        N(I,J)=1d99
      ENDDO
      ENDDO
      MUR=1d99
      MUL=1d99
      MDR=1d99
      MDL=1d99
      MLR=1d99
      MLL=1d99 
      MNL=1d99
      MST1=1d99
      MST2=1d99
      MSB2=1d99
      MSB1=1d99
      MSL1=1d99
      MSL2=1d99
      MSNT=1d99
      CST=1d99
      CSB=1d99
      CSL=1d99
      MSMU1=1d99
      MSMU2=1d99
      MSMUNT=1d99      

*   DEFAULT VALUE FOR FLAGS
      OMGFLAG=0
      MAFLAG=-1
      PFLAG=0
      NMSFLAG=0
      HFLAG=0
      VFLAG=0

*   Set default for Q2

      Q2=1d99

*   INITIALIZE READ LOOP
      NLINE=0
      CHBLCK=' '

*     START TO READ NEW LINE INTO CHINL
 21   CHINL=' '
     
*     LINE NUMBER
      NLINE=NLINE+1
      READ(15,'(A120)',END=29,ERR=999) CHINL
*     CHECK FOR EMPTY OR COMMENT LINES
      IF(CHINL.EQ.' '.OR.CHINL(1:1).EQ.'#'
     .  .OR.CHINL(1:1).EQ.'*') GOTO 21
*     FORCE UPPER CASE LETTERS IN CHINL (AS REQUIRED BELOW)
      INL=0
 22   INL=INL+1
    
      IF(CHINL(INL:INL).NE.'#')THEN
       DO ICH=97,122
        IF(CHINL(INL:INL).EQ.CHAR(ICH)) CHINL(INL:INL)=CHAR(ICH-32)
       ENDDO
       IF(INL.LT.120) GOTO 22
      ENDIF

*   CHECK FOR BLOCK STATEMENT
      IF(CHINL(1:1).EQ.'B')THEN
       READ(CHINL,'(A6,A)',ERR=999) CHDUM,CHBLCK
       GOTO 21
      ENDIF

*   CHECK FOR NMSSM MODEL IN MODSEL
*   IF THE RELIC DENSITY SHOULD BE COMPUTED
*   THE BLOCK MODSEL MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.8) PFLAG=NINT(VAL)
       IF(IX.EQ.9) OMGFLAG=NINT(VAL)  !flag sets whether micromegas is called or not
       IF(IX.EQ.11) GMUFLAG=NINT(VAL) !flag sets whether (g-2) routine called or not
       IF(IX.EQ.12) Q2=VAL**2
       IF(IX.EQ.13) NMSFLAG=NINT(VAL) !flag sets whether NMSDECAY is called or not
       IF(IX.EQ.14) VFLAG=NINT(VAL)   !flag sets whether H->V*V* decays are included
      
*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL

* From PAR array 
* For DECAY I need   PAR(3) = tan(beta)
*                    PAR(12) = AU3
*                    PAR(13) = AD3
*                    PAR(14) = AE3                        
* Additionally 
* For TDECAY I need   PAR(4) =  mueff = lambda <S> 
* For BSG I need      PAR(21) = M2
* For MAGNMU I need   PAR(1) = lambda
* For RELDEN I need   PAR(1) = kappa
*                     PAR(5) = Alambda
*                     PAR(6) = Akappa
*                     PAR(20) = M1
*                     PAR(21) = M2
*                     PAR(22) = M3
*                     PAR(18) = mL2**2
*                     PAR(19) = mE2**2
*                     PAR(10) = mL3**2
*                     PAR(11) = mE3**2   
*                     PAR(15) = mQ2**2
*                     PAR(16) = mU2**2
*                     PAR(17) = mD2**2      
*                     PAR(7) = mQ3**2
*                     PAR(8) = mU3**2
*                     PAR(9) = mD3**2
     
*   READ Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.3) PAR(3)=VAL  ! fills DRbar tan(beta) at MZ
   
*     READ NMSSMRUN
      ELSEIF(CHBLCK(1:8).EQ.'NMSSMRUN')THEN
       READ(CHINL,*,ERR=999) IX,VAL   
       IF(IX.EQ.1)THEN 
          PAR(1)=VAL            ! lambda(SQRT(Q2))
       ENDIF
       IF(IX.EQ.2)THEN 
          PAR(2)=VAL            ! kappa(SQRT(Q2))
       ENDIF
       IF(IX.EQ.3)THEN
          PAR(5)=VAL            ! Alambda(SQRT(Q2))
       ENDIF
       IF(IX.EQ.4)THEN
          PAR(6)=VAL            ! Akappa(SQRT(Q2))
       ENDIF
       IF(IX.EQ.5)THEN
          PAR(4)=VAL            ! mueff(SQRT(Q2))
       ENDIF
       IF(IX.EQ.6)XIF=VAL
       IF(IX.EQ.7)XIS=VAL
       IF(IX.EQ.8)MUP=VAL
       IF(IX.EQ.9)MSP=VAL
       IF(IX.EQ.10)MHSS=VAL
       IF(IX.EQ.12)M3H=VAL
      ELSEIF(CHBLCK(1:4).EQ.'HMIX')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       if(IX.eq.4)MA2=VAL
     
*     READ au
      ELSEIF(CHBLCK(1:2).EQ.'AU')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL
       IF(IX1.EQ.3 .AND. IX2.EQ.3) PAR(12)=VAL  !At = AU3
       IF(IX1.EQ.3 .AND. IX2.EQ.3)ATP=VAL  !At = AU3

*     READ ad
      ELSEIF(CHBLCK(1:2).EQ.'AD')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL
       IF(IX1.EQ.3 .AND. IX2.EQ.3) PAR(13)=VAL!Ab = AD3
      IF(IX1.EQ.3 .AND. IX2.EQ.3) ABP=VAL!Ab = AD3
     
*     READ ae
      ELSEIF(CHBLCK(1:2).EQ.'AE')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL
       IF(IX1.EQ.3 .AND. IX2.EQ.3) PAR(14)=VAL  !Atau = AE3
       IF(IX1.EQ.2 .AND. IX2.EQ.2) PAR(25)=VAL  !Amu = AE3

*     READ msoft
      ELSEIF(CHBLCK(1:5).EQ.'MSOFT')THEN   
       READ(CHINL,*,ERR=999) IX,VAL 
       IF(IX.EQ.43) PAR(7)=VAL**2 !mqL3**2 = mQ3**2  
       IF(IX.EQ.46) PAR(8)=VAL**2 !mtR**2 = mU3**2
       IF(IX.EQ.49) PAR(9)=VAL**2 !mbR**2 = mD3**2
       IF(IX.EQ.33) PAR(10)=VAL**2 !mtauL**2 = mL3**2
       IF(IX.EQ.36) PAR(11)=VAL**2 !mtauR**2 = mE3**2
       IF(IX.EQ.42) PAR(15)=VAL**2 !mqL2**2 = mQ2**2  
       IF(IX.EQ.45) PAR(16)=VAL**2 !mcR**2 = mU2**2
       IF(IX.EQ.48) PAR(17)=VAL**2 !msR**2 = mD2**2
       IF(IX.EQ.32) PAR(18)=VAL**2 !mmuL**2 = mL2**2
       IF(IX.EQ.35) PAR(19)=VAL**2 !mmuR**2 = mE2**2
       IF(IX.EQ.1) PAR(20)=VAL  !M1
       IF(IX.EQ.2) PAR(21)=VAL  !M2
       IF(IX.EQ.3) PAR(22)=VAL  !M3
       IF(IX.EQ.21)MHDS=VAL !mH1^2
       IF(IX.EQ.22)MHUS=VAL !mH2^2
       IF(IX.EQ.43) MQ3P=VAL**2 !mqL3**2 = mQ3**2  
       IF(IX.EQ.46) MU3P=VAL**2 !mtR**2 = mU3**2
       IF(IX.EQ.49) MD3P=VAL**2 !mbR**2 = mD3**2
c$$$       IF(IX.EQ.21)MHSS=VAL !mH1^2      

*     PAR is filled.  Now we fill the masses and mixings.

*     READ MASS
      ELSEIF(CHBLCK(1:4).EQ.'MASS')THEN   
       READ(CHINL,*,ERR=999) IX,VAL  
       IF(IX.EQ.24) MW=VAL  ! Taking softsusy predicted MW.  
       IF(IX.EQ.25) SMASS(1)=VAL !mh0(1).  
       IF(IX.EQ.35) SMASS(2)=VAL !mh0(2).
       IF(IX.EQ.45) SMASS(3)=VAL !mh0(3).  
       IF(IX.EQ.36) PMASS(1)=VAL !mA0(1).
       IF(IX.EQ.46) PMASS(2)=VAL !mA0(2).  
       IF(IX.EQ.37) CMASS=VAL
       IF(IX.EQ.1000021) MGL=VAL
       IF(IX.EQ.1000022) MNEU(1)=VAL
       IF(IX.EQ.1000023) MNEU(2)=VAL
       IF(IX.EQ.1000024) MCHA(1)=VAL
       IF(IX.EQ.1000025) MNEU(3)=VAL
       IF(IX.EQ.1000035) MNEU(4)=VAL
       IF(IX.EQ.1000037) MCHA(2)=VAL
       IF(IX.EQ.1000045) MNEU(5)=VAL
       IF(IX.EQ.1000001) MDL=VAL
       IF(IX.EQ.1000002) MUL=VAL
       !Add warning if the 1st two generations are different by some amount
       IF(IX.EQ.1000003) MDL2=VAL
       IF(IX.EQ.1000004) MUL2=VAL
       IF(IX.EQ.1000005) MSB1=VAL
       IF(IX.EQ.1000006) MST1=VAL
       IF(IX.EQ.1000011) MLL=VAL
       IF(IX.EQ.1000012) MNL=VAL
       IF(IX.EQ.1000013) MSMUL=VAL
       IF(IX.EQ.1000014) MSMUNT=VAL
       IF(IX.EQ.1000015) MSL1=VAL
       IF(IX.EQ.1000016) MSNT=VAL
       IF(IX.EQ.2000001) MDR=VAL
       IF(IX.EQ.2000002) MUR=VAL
       IF(IX.EQ.2000003) MDR2=VAL
       IF(IX.EQ.2000004) MUR2=VAL
       IF(IX.EQ.2000005) MSB2=VAL
       IF(IX.EQ.2000006) MST2=VAL
       IF(IX.EQ.2000011) MLR=VAL
       IF(IX.EQ.2000013) MSMUR=VAL
       IF(IX.EQ.2000015) MSL2=VAL
  
*     READ NMHmix
      ELSEIF(CHBLCK(1:6).EQ.'NMHMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       S(IX1,IX2)=VAL
     
*     READ NMAmix
      ELSEIF(CHBLCK(1:6).EQ.'NMAMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       P(IX1,IX2)=VAL   

*     READ NMNmix
      ELSEIF(CHBLCK(1:6).EQ.'NMNMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       N(IX1,IX2)=VAL


*     READ Umix
      ELSEIF(CHBLCK(1:4).EQ.'UMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       U(IX1,IX2)=VAL


*     READ NVmix
      ELSEIF(CHBLCK(1:4).EQ.'VMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       V(IX1,IX2)=VAL

*     READ STOPMIX
      ELSEIF(CHBLCK(1:7).EQ.'STOPMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       IF(IX1.EQ.1)THEN
          IF(IX2.EQ.1)CST=VAL
        !  IF(IX2.EQ.2)SST=VAL
       ENDIF

*     READ SBOTMIX
      ELSEIF(CHBLCK(1:7).EQ.'SBOTMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       IF(IX1.EQ.1)THEN
          IF(IX2.EQ.1)CSB=VAL
       ENDIF

*     READ STAUMIX
      ELSEIF(CHBLCK(1:7).EQ.'STAUMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       IF(IX1.EQ.1)THEN
          IF(IX2.EQ.1)CSL=VAL
       ENDIF

*     READ SMUMIX
      ELSEIF(CHBLCK(1:6).EQ.'SMUMIX')THEN   
       READ(CHINL,*,ERR=999) IX1,IX2,VAL        
       IF(IX1.EQ.1)THEN
          IF(IX2.EQ.1)CSMU=VAL
       ENDIF
      
      ENDIF 

       GOTO 21

*   END OF READING FROM INPUT FILE
    
 29   ERR=0
      
*   g1,g2  and sin(theta)^2 in the on-shell scheme in terms of
*   GF, MZ(pole) and MW(pole)

      g2=4d0*DSQRT(2d0)*GF*MW**2
      g1=4d0*DSQRT(2d0)*GF*(MZ**2-MW**2)
      S2TW=1d0-(MW/MZ)**2

      IF(Q2.EQ.1D99) THEN
        Q2=DSQRT(MAX(PAR(7)*PAR(8),100d0**2))
      ENDIF
!      QSTSB=Q2

*   Set Q2MIN, Q2FIX:      

      Q2MIN=100d0**2
      Q2FIX=1
      IF(Q2.LE.Q2MIN)THEN
       Q2FIX=0
      ENDIF

      IF(MA2.LT.0d0)THEN
         WRITE(0,1)"MA2 MUST BE GIVEN WITH POSITIVE VALUE" 
         ERR=1
      ENDIF
      PAR(23)=DSQRT(MA2)
     
      P2(1,1) = P(1,1)/ SIN(ATAN(PAR(3)))
      P2(2,1) = P(2,1)/ SIN(ATAN(PAR(3)))
      P2(1,2) = P(1,3)
      P2(2,2) = P(2,3)

       
!     Conventions for CP EVEN are not SLHA
       SCOMP(1,2)=S(1,1)
       SCOMP(1,1)=S(1,2)
       SCOMP(1,3)=S(1,3)
       SCOMP(2,2)=S(2,1)
       SCOMP(2,1)=S(2,2)
       SCOMP(2,3)=S(2,3)
       SCOMP(3,2)=S(3,1)
       SCOMP(3,1)=S(3,2)
       SCOMP(3,3)=S(3,3)

!     Conventions for neutralino mixing are not SLHA
       NEU=N
       NEU(1,3)=N(1,4)
       NEU(1,4)=N(1,3)
       NEU(2,3)=N(2,4)
       NEU(2,4)=N(2,3)
       NEU(3,3)=N(3,4)
       NEU(3,4)=N(3,3)
       NEU(4,3)=N(4,4)
       NEU(4,4)=N(4,3)
       NEU(5,3)=N(5,4)
       NEU(5,4)=N(5,3)
   
*     MASS ORDER SMUONS
      IF(MSMUL.LT.MSMUR)THEN
         MSMU1=MSMUL
         MSMU2=MSMUR
      ELSE
         MSMU1=MSMUR
         MSMU2=MSMUL
      ENDIF

*   Check for errors

      IF(SMASS(1).EQ.1d99)THEN
         WRITE(0,1)"mh0(1) MUST BE GIVEN IN BLOCK MASS" 
         ERR=1
      ENDIF
      IF(SMASS(2).EQ.1d99)THEN
         WRITE(0,1)"mh0(2) MUST BE GIVEN IN BLOCK MASS" 
         ERR=1
      ENDIF
      IF(SMASS(3).EQ.1d99)THEN
         WRITE(0,1)"mh0(3) MUST BE GIVEN IN BLOCK MASS" 
         ERR=1
      ENDIF
        IF(PMASS(1).EQ.1d99)THEN
         WRITE(0,1)"mA0(1) MUST BE GIVEN IN BLOCK MASS" 
         ERR=1
      ENDIF
      IF(PMASS(2).EQ.1d99)THEN
         WRITE(0,1)"mA0(2) MUST BE GIVEN IN BLOCK MASS" 
         ERR=1
      ENDIF
      IF(CMASS.EQ.1d99)THEN
       WRITE(0,1)"MHPM MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      DO I=1,3
       DO J=1,3
        IF(S(I,J).EQ.1d99)THEN
         WRITE(0,1)"S MUST BE FILLED IN BLOCK NMHMIX" 
         ERR=1
         ENDIF
      ENDDO
      ENDDO
      DO I=1,2
       DO J=1,3
        IF(P(I,J).EQ.1d99)THEN
         WRITE(0,1)"P MUST BE FILLED IN BLOCK NMAMIX" 
         ERR=1
         ENDIF
      ENDDO
      ENDDO
      IF(MCHA(1).EQ.1d99)THEN
         WRITE(0,1)"MCH(1) MUST BE GIVEN IN BLOCK MASS" 
         ERR=1
      ENDIF
      IF(MCHA(2).EQ.1d99)THEN
       WRITE(0,1)"MCH(2) MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      DO I=1,5
       IF(MNEU(I).EQ.1d99)THEN
        WRITE(0,1)"MNEU MUST BE GIVEN IN BLOCK MASS" 
        ERR=1
      ENDIF
      ENDDO
       DO I=1,2
       DO J=1,2
        IF(V(I,J).EQ.1d99)THEN
         WRITE(0,1)"CHARGINO V MATRIX MUST BE FILLED IN BLOCK VMIX" 
         ERR=1
         ENDIF
      ENDDO
      ENDDO
       DO I=1,2
       DO J=1,2
        IF(U(I,J).EQ.1d99)THEN
         WRITE(0,1)"CHARGINO U MATRIX MUST BE FILLED IN BLOCK UMIX" 
         ERR=1
      ENDIF
      ENDDO
      ENDDO
       DO I=1,5
       DO J=1,5
        IF(N(I,J).EQ.1d99)THEN
         WRITE(0,1)"NEUTRALINO N MATRIX MUST BE FILLED IN BLOCK NMNMIX" 
         ERR=1
      ENDIF
      ENDDO
      ENDDO
      IF(MUR.EQ.1d99)THEN
       WRITE(0,1)"MUR MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MUL.EQ.1d99)THEN
       WRITE(0,1)"MUL MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MDR.EQ.1d99)THEN
       WRITE(0,1)"MDR MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MLR.EQ.1d99)THEN
       WRITE(0,1)"MLR MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MLL.EQ.1d99)THEN
       WRITE(0,1)"MLL MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MNL.EQ.1d99)THEN
       WRITE(0,1)"MNL MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MST1.EQ.1d99)THEN
       WRITE(0,1)"MST1 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MST2.EQ.1d99)THEN
       WRITE(0,1)"MST2 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MSB1.EQ.1d99)THEN
       WRITE(0,1)"MSB1 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MSB2.EQ.1d99)THEN
       WRITE(0,1)"MSB2 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MSL1.EQ.1d99)THEN
       WRITE(0,1)"MSL1 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MSL2.EQ.1d99)THEN
       WRITE(0,1)"MSL2 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MSNT.EQ.1d99)THEN
       WRITE(0,1)"MSNT MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MSMU1.EQ.1d99)THEN
       WRITE(0,1)"MSMU1 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(MSMU2.EQ.1d99)THEN
       WRITE(0,1)"MSMU2 MUST BE GIVEN IN BLOCK MASS" 
       ERR=1
      ENDIF
      IF(CST.EQ.1d99)THEN
       WRITE(0,1)"STOP MIXING MUST BE GIVEN IN BLOCK STOPMIX" 
       ERR=1
      ENDIF
      IF(CSB.EQ.1d99)THEN
       WRITE(0,1)"SBOTTOM MIXING MUST BE GIVEN IN BLOCK SBOTMIX" 
       ERR=1
      ENDIF
      IF(CSL.EQ.1d99)THEN
       WRITE(0,1)"STAU MIXING MUST BE GIVEN IN BLOCK STAUMIX" 
       ERR=1
      ENDIF      
      
*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF
      


      RETURN
      
 999  WRITE(0,1)"READ ERROR ON LINE:", NLINE
!     WRITE(0,*)CHINL(1:80)
      STOP 1
      
 1    FORMAT(A)
      
      END

      SUBROUTINE NONSLHAINPUT
      INTEGER N0
      INTEGER NLOOP,NBER
      DOUBLE PRECISION ACC,XITLA
      DOUBLE PRECISION MSB,MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW

*     Initialization for ALPHAS and RUNM (as in hdecay)
*     The bottom quark pole mass MBP is set in INIT and can be changed
*     only there (changing its running mass MB above has no effect
*     on MBP, since one would have to compute alpha_s(MB) first)
      
      MSB=MS
      MC0=MC
      MB0=MBP
      MT0=MT
      N0=5
      NLOOP=2
      NBER=18
      ACC=1.D-8
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      CALL ALSINI(ACC)
      CALL BERNINI(NBER)

      RETURN
      
      END
      

      SUBROUTINE OUTPUT(PAR,PROB,IFAIL)

*********************************************************************      
*   Subroutine writing all the results in the the output files.
*********************************************************************      
 
      IMPLICIT NONE

      INTEGER I,NBIN,IFAIL
      INTEGER NMSFLAG,OMGFLAG,MAFLAG,PFLAG
      INTEGER NSUSY,NGUT,NMES
      PARAMETER (NSUSY=14,NGUT=21,NMES=15)

      DOUBLE PRECISION PAR(*),PROB(*),SIG(3,8)
      DOUBLE PRECISION SMASS(3),AMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,TANB,SINB,COSB
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION SST,SSB,SSL,Q2
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX,LOAMASS
      DOUBLE PRECISION Xf,sigmaV,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION DELMB,PX,PA(6),PB(2),PL(7),PK(8),MH(3),MMH(3)
      DOUBLE PRECISION DMH(3),MA(2),MMA(2),DMA(2),MHC,MMHC,DMHC
      DOUBLE PRECISION MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      DOUBLE PRECISION LUX,PRINTCHANNELS,omg_

      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,AMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFHIGM/MH,MMH,DMH,MA,MMA,DMA,MHC,MMHC,DMHC
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax,chi2max,chi2gam,chi2bb,chi2zz
      
      TANB=PAR(3)
      COSB=1d0/DSQRT(1d0+TANB**2)
      SINB=TANB*COSB

      WRITE(17,899) "# NMSSMTools OUTPUT IN SLHA FORMAT"
      WRITE(17,899) "# Info about spectrum calculator"
      WRITE(17,899) "BLOCK SPINFO   # Program information"
      WRITE(17,900) 1,"NMSSMTools # Spectrum calculator"
      WRITE(17,900) 2,"4.1.1      # Version number"

      IF(PROB(1).NE.0d0)
     . WRITE(17,900) 3,"# Chargino too light"
      IF(PROB(2).NE.0d0)
     . WRITE(17,900) 3,"# Neutralinos too light"
      IF(PROB(3).NE.0d0)
     . WRITE(17,900) 3,"# Charged Higgs too light"
      IF(PROB(4).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, ind. of h decay"
      IF(PROB(5).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> bb"
      IF(PROB(6).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> tautau"
      IF(PROB(7).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> invisible"
      IF(PROB(8).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> 2jets"
      IF(PROB(9).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> 2photons"
      IF(PROB(10).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> AA -> 4bs"
      IF(PROB(11).NE.0d0 .OR. PROB(41).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> AA -> 4taus"
      IF(PROB(12).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> AA -> 2bs 2taus"
      IF(PROB(19).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hZ, h -> AA,A -> light pair"
      IF(PROB(13).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> Z -> hA (Z width)"
      IF(PROB(14).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hA -> 4bs"
      IF(PROB(15).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hA -> 4taus"
      IF(PROB(16).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hA -> 2bs 2taus"
      IF(PROB(17).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hA -> AAA -> 6bs"
      IF(PROB(18).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by ee -> hA -> AAA -> 6taus"
      IF(PROB(20).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by stop -> b l sneutrino"
      IF(PROB(21).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by stop -> neutralino c"
      IF(PROB(22).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by sbottom -> neutralino b"
      IF(PROB(23).NE.0d0)
     . WRITE(17,900) 3,"# Squark/gluino too light"
      IF(PROB(24).NE.0d0)
     . WRITE(17,900) 3,"# Selectron/smuon too light"
      IF(PROB(25).NE.0d0)
     . WRITE(17,900) 3,"# Stau too light"
      IF(PROB(26).GT.0d0)
     . WRITE(17,900) 3,"# Lightest neutralino is not the LSP"
      IF(PROB(26).LT.0d0)
     . WRITE(17,900) 3,"# Mass of the lightest neutralino < 511 keV"
      IF(PROB(27).NE.0d0)
     . WRITE(17,900) 3,"# Landau Pole below MGUT"
      IF(PROB(28).NE.0d0)
     . WRITE(17,900) 3,"# Unphysical global minimum"
      IF(PROB(29).NE.0d0)
     . WRITE(17,900) 3,"# Higgs soft masses >> Msusy"
      IF(PROB(30).GT.0d0)WRITE(17,900) 3,
     . "# Relic density too large (WMAP)"
      IF(PROB(30).LT.0d0.AND.PROB(30).NE.-1d0)
     . WRITE(17,900) 3,"# Relic density too small (WMAP)"
      IF(PROB(30).EQ.-1d0)
     . WRITE(17,900) 3,"# Problem in micrOMEGAs"
      IF(PROB(31).NE.0d0)WRITE(17,900) 3,
     . "# Excluded by LUX"
      IF(PROB(32).NE.0d0)
     . WRITE(17,900) 3,"# b -> s gamma more than 2 sigma away"
      IF(PROB(33).NE.0d0)
     . WRITE(17,900) 3,"# Delta M_s more than 2 sigma away"
      IF(PROB(34).NE.0d0)
     . WRITE(17,900) 3,"# Delta M_d more than 2 sigma away"
      IF(PROB(35).NE.0d0)
     . WRITE(17,900) 3,"# B_s -> mu+ mu- more than 2 sigma away"
      IF(PROB(36).NE.0d0)
     . WRITE(17,900) 3,"# B+ -> tau nu_tau more than 2 sigma away"
      IF(PROB(37).NE.0d0)
     . WRITE(17,900) 3,"# Muon magn. mom. more than 2 sigma away"
      IF(PROB(38).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by Upsilon(1S) -> A gamma (CLEO)"
      IF(PROB(38).LT.0d0)
     . WRITE(17,900) 3,"# (but A width> 10 MeV)"
      IF(PROB(39).NE.0d0)
     . WRITE(17,900) 3,
     . "# Excluded etab(1S) mass difference (BABAR - theory)"
       IF(PROB(40).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by BR(B -> X_s mu +mu-)"
       IF(PROB(42).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by top -> b H+, H+ -> c s"
       IF(PROB(43).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by top -> b H+, H+ -> tau nu_tau"
       IF(PROB(44).NE.0d0)
     . WRITE(17,900) 3,
     . "# Excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus"
       IF(PROB(45).NE.0d0)
     . WRITE(17,900) 3,"# Excluded by t -> bH+ (LHC)"
       IF(PROB(46).NE.0d0)
     . WRITE(17,918) 3,"# No Higgs in the",MHMIN,MHMAX," GeV mass range"
       IF(PROB(47).NE.0d0)
     . WRITE(17,919) 3,"# chi2(H->gg) > ",chi2MAX
       IF(PROB(48).NE.0d0)
     . WRITE(17,919) 3,"# chi2(H->bb) > ",chi2MAX
       IF(PROB(49).NE.0d0)
     . WRITE(17,919) 3,"# chi2(H->ZZ) > ",chi2MAX

      IF(IFAIL.EQ.1.OR.IFAIL.EQ.3.OR.IFAIL.EQ.5.OR.IFAIL.EQ.7)
     . WRITE(17,900) 4,"# M_H1^2<1"
      IF(IFAIL.EQ.2.OR.IFAIL.EQ.3.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)
     . WRITE(17,900) 4,"# M_A1^2<1"
      IF(IFAIL.EQ.4.OR.IFAIL.EQ.5.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)
     . WRITE(17,900) 4,"# M_HC^2<1"
      IF(IFAIL.EQ.8)
     . WRITE(17,900) 4,"# Negative sfermion mass squared"
      IF(IFAIL.EQ.9)
     . WRITE(17,900) 4,"# Disallowed parameters: lambda or tan(beta)=0"
      IF(IFAIL.EQ.11)
     . WRITE(17,900) 4,"# Integration problem in RGES"
      IF(IFAIL.EQ.12)
     . WRITE(17,900) 4,"# Integration problem in RGESUNI"
      IF(IFAIL.EQ.13)
     . WRITE(17,900) 4,"# Integration problem in RGESINV"
      IF(IFAIL.EQ.14.OR.IFAIL.EQ.15)
     . WRITE(17,900) 4,"# Convergence Problem"
      IF(IFAIL.EQ.16)
     . WRITE(17,900) 4,"# No electroweak symmetry breaking"

      WRITE(17,899) "# Input parameters"
      WRITE(17,899) "BLOCK MODSEL"
      WRITE(17,921) 3,1,"NMSSM particle content"
      WRITE(17,921) 1,1,"IMOD"
      WRITE(17,921) 10,0,"ISCAN"
      WRITE(17,921) 9,OMGFLAG,"Call micrOmegas"
      WRITE(17,921) 8,PFLAG,"Precision for Higgs masses"
      WRITE(17,921) 13,NMSFLAG,"Sparticle decays via NMSDECAY"

      WRITE(17,899) "BLOCK SMINPUTS"
      WRITE(17,901) 1,1d0/ALEMMZ,"ALPHA_EM^-1(MZ)"
      WRITE(17,901) 2,GF,"GF"
      WRITE(17,901) 3,ALSMZ,"ALPHA_S(MZ)"
      WRITE(17,901) 4,MZ,"MZ"
      WRITE(17,901) 5,MB,"MB(MB)"
      WRITE(17,901) 6,MT,"MTOP (POLE MASS)"
      WRITE(17,901) 7,MTAU,"MTAU"
      WRITE(17,899) "# SMINPUTS Beyond SLHA:"
      WRITE(17,906) "MW:",MW
      WRITE(17,906) "MS:",MS
      WRITE(17,906) "MC:",MC
      WRITE(17,906) "VUS:",VUS
      WRITE(17,906) "VCB:",VCB
      WRITE(17,906) "VUB:",VUB
      
      WRITE(17,899) "BLOCK MINPAR"
      WRITE(17,901) 0,DSQRT(Q2),"REN. SCALE"
      WRITE(17,901) 3,TANB,"TANBETA(MZ)"

      WRITE(17,899) "BLOCK EXTPAR"       
      WRITE(17,901) 61,PAR(1),"LAMBDA AT THE SUSY SCALE"
      WRITE(17,901) 62,PAR(2),"KAPPA AT THE SUSY SCALE"
      WRITE(17,901) 65,PAR(4),"MUEFF AT THE SUSY SCALE"
      
      WRITE(17,899) "# "
      WRITE(17,899) "# NMSSM SPECIFIC PARAMETERS THE SUSY SCALE"
      WRITE(17,907) "BLOCK NMSSMRUN Q=",DSQRT(Q2),
     .   " # (INPUTS AT THE SUSY SCALE)"
      WRITE(17,901) 1,PAR(1),"LAMBDA(Q,DR_bar)"
      WRITE(17,901) 2,PAR(2),"KAPPA(Q,DR_bar)"
      WRITE(17,901) 3,PAR(5),"ALAMBDA"
      WRITE(17,901) 4,PAR(6),"AKAPPA"
      WRITE(17,901) 5,PAR(4),"MUEFF"
      WRITE(17,901) 6,XIFSUSY,"XIF"
      WRITE(17,901) 7,XISSUSY,"XIS"
      WRITE(17,901) 8,MUPSUSY,"MU'"
      WRITE(17,901) 9,MSPSUSY,"MS'^2"
      WRITE(17,901) 10,MSS,"MS^2"
      WRITE(17,901) 12,M3HSUSY,"M3H^2"
     
      WRITE(17,899) "# "
      WRITE(17,907) "BLOCK HMIX Q=",DSQRT(QSTSB),
     .    " # (COMPUTED AT THE SCALE OF STOP/SBOTTOM MASSES)"
      WRITE(17,901) 1,MUQ,"MUEFF"
      WRITE(17,901) 2,TANBQ,"TAN(BETA)"
      WRITE(17,901) 3,DSQRT(2d0*(H1Q**2+H2Q**2)),"V(Q)"
      WRITE(17,901) 4,PAR(23)**2,"MA^2"            

      IF(IFAIL.NE.0.AND.IFAIL.NE.10) GOTO 1
      
      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK MASS   # Mass spectrum "
      WRITE(17,899) "#  PDG Code     mass             particle "
      WRITE(17,902) 25,SMASS(1),"lightest neutral scalar"
      WRITE(17,902) 35,SMASS(2),"second neutral scalar"
      WRITE(17,902) 45,SMASS(3),"third neutral scalar"
      WRITE(17,902) 36,AMASS(1),"lightest pseudoscalar"
      WRITE(17,902) 46,AMASS(2),"second pseudoscalar"
      WRITE(17,902) 37,CMASS,"charged Higgs"
      WRITE(17,902) 1000001,MDL," ~d_L"
      WRITE(17,902) 2000001,MDR," ~d_R"
      WRITE(17,902) 1000002,MUL," ~u_L"
      WRITE(17,902) 2000002,MUR," ~u_R"       
      WRITE(17,902) 1000003,MDL," ~s_L"
      WRITE(17,902) 2000003,MDR," ~s_R"       
      WRITE(17,902) 1000004,MUL," ~c_L"
      WRITE(17,902) 2000004,MUR," ~c_R"       
      WRITE(17,902) 1000005,MSB1," ~b_1"
      WRITE(17,902) 2000005,MSB2," ~b_2"       
      WRITE(17,902) 1000006,MST1," ~t_1"
      WRITE(17,902) 2000006,MST2," ~t_2"       
      WRITE(17,902) 1000011,MLL," ~e_L"
      WRITE(17,902) 2000011,MLR," ~e_R"
      WRITE(17,902) 1000012,MNL," ~nue_L"
      WRITE(17,902) 1000013,MLL," ~mu_L"
      WRITE(17,902) 2000013,MLR," ~mu_R"
      WRITE(17,902) 1000014,MSMUNT," ~numu_L"
      WRITE(17,902) 1000015,MSL1," ~tau_1"
      WRITE(17,902) 2000015,MSL2," ~tau_2"
      WRITE(17,902) 1000016,MSNT," ~nutau_L"
      WRITE(17,902) 1000021,MGL," ~g"
      WRITE(17,902) 1000022,MNEU(1),"neutralino(1)"
      WRITE(17,902) 1000023,MNEU(2),"neutralino(2)"
      WRITE(17,902) 1000025,MNEU(3),"neutralino(3)"
      WRITE(17,902) 1000035,MNEU(4),"neutralino(4)"
      WRITE(17,902) 1000045,MNEU(5),"neutralino(5)"
      WRITE(17,902) 1000024,MCHA(1),"chargino(1)"
      WRITE(17,902) 1000037,MCHA(2),"chargino(2)"
      
      WRITE(17,899) "# " 
      WRITE(17,899) "# Low energy observables"
      WRITE(17,899) "BLOCK LOWEN"
      WRITE(17,899)
     .   "# Exp. 2 Sigma: 3.04E-4 < BR(b -> s gamma) < 4.06E-4:"
      WRITE(17,901) 1,BRSG,"BR(b -> s gamma)"
      WRITE(17,901) 11,BRSGMAX,"(BR(b -> s gamma)+Theor.Err.)"
      WRITE(17,901) 12,BRSGMIN,"(BR(b -> s gamma)-Theor.Err.)"
      WRITE(17,899) "# Exp. 2 Sigma: 4.99E-1 < Delta M_d < 5.15E-1:"
      WRITE(17,901) 2,DMD,"Delta M_d in ps^-1"
      WRITE(17,901) 21,DMdmax,"Delta M_d +Theor.Err."
      WRITE(17,901) 22,DMdmin,"Delta M_d -Theor.Err."
      WRITE(17,899) 
     .   "# Exp. 2 Sigma: 1.7633E+1 < Delta Ms < 1.7805E+1:"
      WRITE(17,901) 3,DMS,"Delta M_s in ps^-1"
      WRITE(17,901) 31,DMsmax,"Delta M_s +Theor.Err."
      WRITE(17,901) 32,DMsmin,"Delta M_s -Theor.Err."
      WRITE(17,899) "# Exp. 2 Sigma: 2.0E-9 < BR(Bs->mu+mu-) < 4.7E-9:"
      WRITE(17,901) 4,BRBMUMU,"BR(Bs -> mu+mu-)"
      WRITE(17,901) 41,BRBMUMUmax,"BR(Bs -> mu+mu-)+Theor.Err."
      WRITE(17,901) 42,BRBMUMUmin,"BR(Bs -> mu+mu-)-Theor.Err."
      WRITE(17,899) 
     .   "# Exp. 2 Sigma: 1.07E-4 < BR(B+ > tau+ + nu_tau) < 2.27E-4:"
      WRITE(17,901) 5,BRBtaunu,"BR(B+ -> tau+ + nu_tau)"
      WRITE(17,901) 51,BRBtaunumax,
     .   "BR(B+ -> tau+ + nu_tau) + Theor.Err."
      WRITE(17,901) 52,BRBtaunumin,
     .   "BR(B+ -> tau+ + nu_tau) - Theor.Err."
      WRITE(17,899) "# " 
      WRITE(17,899) "# BSM contr. to the muon anomalous magn. moment:"
      WRITE(17,901) 6,delmagmu,"Del_a_mu"
      WRITE(17,901) 61,amuthmax,"Del_a_mu + Theor.Err."
      WRITE(17,901) 62,amuthmin,"Del_a_mu - Theor.Err."
      WRITE(17,907) "# Minimal Exp.-SM (2 sigma):",damumin
      WRITE(17,907) "# Maximal Exp.-SM (2 sigma):",damumax
     
      IF(OMGFLAG.NE.0)THEN
	WRITE(17,899) "# " 
	WRITE(17,911)
     .   "# Omega h^2 (allowed:",OMGMIN," < Omega h^2 <",OMGMAX,"):"
	IF(OMG.EQ.0d0)THEN
	  WRITE(17,899) "# Cannot compute Omega h^2 (mLSP < 1 GeV)"
	ELSEIF(OMG.EQ.-1d0)THEN
	  WRITE(17,899) 
     .      "# Lightest neutralino is not the LSP in micrOMEGAs"
	ELSEIF(OMG.EQ.-2d0)THEN
	  WRITE(17,899) "# Problem in micrOMEGAs"
	ELSE
	  WRITE(17,901) 10,OMG,"Omega h^2"
	  omg_=printChannels(Xf,1.D-3,1.D-4,1,17)
	ENDIF
      ENDIF
      IF(OMGFLAG.EQ.2 .OR. OMGFLAG.EQ.4)THEN
	WRITE(17,907)"# sigma(p)_SI (allowed: sigma_p^SI < ",
     .   LUX(DABS(MNEU(1))),"):"
	WRITE(17,901) 20,CSPSI,"sigma_p^SI"
	WRITE(17,915)"# values used for sigma_piN,sigma_S",
     .  " (strange content of the proton)"
	WRITE(17,901) 30,sigmapiN,"sigma_piN"
	WRITE(17,901) 40,sigmaS,"sigma_S"
      ENDIF
!      IF(OMGFLAG.NE.0)CALL printRelDen(PROB,17)
      
      WRITE(17,899) "# "
      WRITE(17,899) "# 3*3 Higgs mixing"
      WRITE(17,899) "BLOCK NMHMIX"
      WRITE(17,903) 1,1,SCOMP(1,2),"S_(1,1)"
      WRITE(17,903) 1,2,SCOMP(1,1),"S_(1,2)"
      WRITE(17,903) 1,3,SCOMP(1,3),"S_(1,3)"
      WRITE(17,903) 2,1,SCOMP(2,2),"S_(2,1)"
      WRITE(17,903) 2,2,SCOMP(2,1),"S_(2,2)"
      WRITE(17,903) 2,3,SCOMP(2,3),"S_(2,3)"
      WRITE(17,903) 3,1,SCOMP(3,2),"S_(3,1)"
      WRITE(17,903) 3,2,SCOMP(3,1),"S_(3,2)"
      WRITE(17,903) 3,3,SCOMP(3,3),"S_(3,3)"

      WRITE(17,899) "# "      
      WRITE(17,899) "# 3*3 Pseudoscalar Higgs mixing"
      WRITE(17,899) "BLOCK NMAMIX"
      WRITE(17,903) 1,1,SINB*PCOMP(1,1),"P_(1,1)"
      WRITE(17,903) 1,2,COSB*PCOMP(1,1),"P_(1,2)"
      WRITE(17,903) 1,3,PCOMP(1,2),"P_(1,3)"
      WRITE(17,903) 2,1,SINB*PCOMP(2,1),"P_(2,1)"
      WRITE(17,903) 2,2,COSB*PCOMP(2,1),"P_(2,2)"
      WRITE(17,903) 2,3,PCOMP(2,2),"P_(2,3)"
           
      SST=DSQRT(1-CST**2)
      SSB=DSQRT(1-CSB**2)
      SSL=DSQRT(1-CSL**2)

      WRITE(17,899) "# "
      WRITE(17,899) "# 3rd generation sfermion mixing"
      WRITE(17,899) "BLOCK STOPMIX  # Stop mixing matrix"
      WRITE(17,903) 1,1,CST,"Rst_(1,1)"
      WRITE(17,903) 1,2,SST,"Rst_(1,2)"
      WRITE(17,903) 2,1,-SST,"Rst_(2,1)"
      WRITE(17,903) 2,2,CST,"Rst_(2,2)"
      WRITE(17,899) "BLOCK SBOTMIX  # Sbottom mixing matrix"
      WRITE(17,903) 1,1,CSB,"Rsb_(1,1)"
      WRITE(17,903) 1,2,SSB,"Rsb_(1,2)"
      WRITE(17,903) 2,1,-SSB,"Rsb_(2,1)"
      WRITE(17,903) 2,2,CSB,"Rsb_(2,2)"
      WRITE(17,899) "BLOCK STAUMIX  # Stau mixing matrix"
      WRITE(17,903) 1,1,CSL,"Rsl_(1,1)"
      WRITE(17,903) 1,2,SSL,"Rsl_(1,2)"
      WRITE(17,903) 2,1,-SSL,"Rsl_(2,1)"
      WRITE(17,903) 2,2,CSL,"Rsl_(2,2)"

      WRITE(17,899) "# "
      WRITE(17,899) "# Gaugino-Higgsino mixing"
      WRITE(17,899) "BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix"
      WRITE(17,903) 1,1,NEU(1,1),"N_(1,1)"
      WRITE(17,903) 1,2,NEU(1,2),"N_(1,2)"
      WRITE(17,903) 1,3,NEU(1,4),"N_(1,3)"
      WRITE(17,903) 1,4,NEU(1,3),"N_(1,4)"
      WRITE(17,903) 1,5,NEU(1,5),"N_(1,5)"
      WRITE(17,903) 2,1,NEU(2,1),"N_(2,1)"
      WRITE(17,903) 2,2,NEU(2,2),"N_(2,2)"
      WRITE(17,903) 2,3,NEU(2,4),"N_(2,3)"
      WRITE(17,903) 2,4,NEU(2,3),"N_(2,4)"
      WRITE(17,903) 2,5,NEU(2,5),"N_(2,5)"
      WRITE(17,903) 3,1,NEU(3,1),"N_(3,1)"
      WRITE(17,903) 3,2,NEU(3,2),"N_(3,2)"
      WRITE(17,903) 3,3,NEU(3,4),"N_(3,3)"
      WRITE(17,903) 3,4,NEU(3,3),"N_(3,4)"
      WRITE(17,903) 3,5,NEU(3,5),"N_(3,5)"
      WRITE(17,903) 4,1,NEU(4,1),"N_(4,1)"
      WRITE(17,903) 4,2,NEU(4,2),"N_(4,2)"
      WRITE(17,903) 4,3,NEU(4,4),"N_(4,3)"
      WRITE(17,903) 4,4,NEU(4,3),"N_(4,4)"
      WRITE(17,903) 4,5,NEU(4,5),"N_(4,5)"
      WRITE(17,903) 5,1,NEU(5,1),"N_(5,1)"
      WRITE(17,903) 5,2,NEU(5,2),"N_(5,2)"
      WRITE(17,903) 5,3,NEU(5,4),"N_(5,3)"
      WRITE(17,903) 5,4,NEU(5,3),"N_(5,4)"
      WRITE(17,903) 5,5,NEU(5,5),"N_(5,5)"

      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK UMIX  # Chargino U Mixing Matrix"
      WRITE(17,903) 1,1,U(1,1),"U_(1,1)"
      WRITE(17,903) 1,2,U(1,2),"U_(1,2)"
      WRITE(17,903) 2,1,U(2,1),"U_(2,1)"
      WRITE(17,903) 2,2,U(2,2),"U_(2,2)"

      WRITE(17,899) "# "
      WRITE(17,899) "BLOCK VMIX  # Chargino V Mixing Matrix"
      WRITE(17,903) 1,1,V(1,1),"V_(1,1)"
      WRITE(17,903) 1,2,V(1,2),"V_(1,2)"
      WRITE(17,903) 2,1,V(2,1),"V_(2,1)"
      WRITE(17,903) 2,2,V(2,2),"V_(2,2)"

      WRITE(17,899) "# "
      WRITE(17,899) "# Higgs reduced couplings"
      WRITE(17,899) "# (as compared to a SM Higgs with same mass)"
      WRITE(17,899) "BLOCK REDCOUP"
      WRITE(17,899) "# H1"
      WRITE(17,903) 1,1,CU(1),"U-type fermions"
      WRITE(17,903) 1,2,CD(1),"D-type fermions"
      WRITE(17,903) 1,3,CV(1),"W,Z bosons"
      WRITE(17,903) 1,4,CJ(1),"Gluons"
      WRITE(17,903) 1,5,CG(1),"Photons"
      WRITE(17,899) "# H2"
      WRITE(17,903) 2,1,CU(2),"U-type fermions"
      WRITE(17,903) 2,2,CD(2),"D-type fermions"
      WRITE(17,903) 2,3,CV(2),"W,Z bosons"
      WRITE(17,903) 2,4,CJ(2),"Gluons"
      WRITE(17,903) 2,5,CG(2),"Photons"
      WRITE(17,899) "# H3"
      WRITE(17,903) 3,1,CU(3),"U-type fermions"
      WRITE(17,903) 3,2,CD(3),"D-type fermions"
      WRITE(17,903) 3,3,CV(3),"W,Z bosons"
      WRITE(17,903) 3,4,CJ(3),"Gluons"
      WRITE(17,903) 3,5,CG(3),"Photons"
      WRITE(17,899) "# A1"
      WRITE(17,903) 4,1,CU(4),"U-type fermions"
      WRITE(17,903) 4,2,CD(4),"D-type fermions"
      WRITE(17,903) 4,3,0.,"W,Z bosons"
      WRITE(17,903) 4,4,CJ(4),"Gluons"
      WRITE(17,903) 4,5,CG(4),"Photons"
      WRITE(17,899) "# A2"
      WRITE(17,903) 5,1,CU(5),"U-type fermions"
      WRITE(17,903) 5,2,CD(5),"D-type fermions"
      WRITE(17,903) 5,3,0.,"W,Z bosons"
      WRITE(17,903) 5,4,CJ(5),"Gluons"
      WRITE(17,903) 5,5,CG(5),"Photons"

      WRITE(17,899) "# "
      WRITE(17,899) "# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE"
      WRITE(17,907) "BLOCK GAUGE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,901) 1,DSQRT(G1S),"g1(Q,DR_bar)"
      WRITE(17,901) 2,DSQRT(G2S),"g2(Q,DR_bar)"
      WRITE(17,901) 3,DSQRT(G3S),"g3(Q,DR_bar)"
      
      WRITE(17,907) "BLOCK YU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,HTOPS,"HTOP(Q,DR_bar)"
      WRITE(17,907) "BLOCK YD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,HBOTS,"HBOT(Q,DR_bar)"
      WRITE(17,907) "BLOCK YE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,HTAUS,"HTAU(Q,DR_bar)"

      WRITE(17,899) "# "
      WRITE(17,899) "# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE"
      WRITE(17,907) "BLOCK AU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,PAR(12),"ATOP"
      WRITE(17,907) "BLOCK AD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 3,3,PAR(13),"ABOT"
      WRITE(17,907) "BLOCK AE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,903) 2,2,PAR(25),"AMUON"
      WRITE(17,903) 3,3,PAR(14),"ATAU"

      WRITE(17,899) "# "
      WRITE(17,899) "# SOFT MASSES AT THE SUSY SCALE"
      WRITE(17,907) "BLOCK MSOFT Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(17,901) 1,PAR(20),"M1"
      WRITE(17,901) 2,PAR(21),"M2"
      WRITE(17,901) 3,PAR(22),"M3"
      WRITE(17,901) 21,MHDS,"M_HD^2"
      WRITE(17,901) 22,MHUS,"M_HU^2"
      WRITE(17,901) 31,PAR(18)/DSQRT(DABS(PAR(18))),"M_eL"
      WRITE(17,901) 32,PAR(18)/DSQRT(DABS(PAR(18))),"M_muL"
      WRITE(17,901) 33,PAR(10)/DSQRT(DABS(PAR(10))),"M_tauL"
      WRITE(17,901) 34,PAR(19)/DSQRT(DABS(PAR(19))),"M_eR"
      WRITE(17,901) 35,PAR(19)/DSQRT(DABS(PAR(19))),"M_muR"
      WRITE(17,901) 36,PAR(11)/DSQRT(DABS(PAR(11))),"M_tauR"
      WRITE(17,901) 41,PAR(15)/DSQRT(DABS(PAR(15))),"M_q1L"
      WRITE(17,901) 42,PAR(15)/DSQRT(DABS(PAR(15))),"M_q2L"
      WRITE(17,901) 43,PAR(7)/DSQRT(DABS(PAR(7))),"M_q3L"
      WRITE(17,901) 44,PAR(16)/DSQRT(DABS(PAR(16))),"M_uR"
      WRITE(17,901) 45,PAR(16)/DSQRT(DABS(PAR(16))),"M_cR"
      WRITE(17,901) 46,PAR(8)/DSQRT(DABS(PAR(8))),"M_tR"
      WRITE(17,901) 47,PAR(17)/DSQRT(DABS(PAR(17))),"M_dR"
      WRITE(17,901) 48,PAR(17)/DSQRT(DABS(PAR(17))),"M_sR"
      WRITE(17,901) 49,PAR(9)/DSQRT(DABS(PAR(9))),"M_bR"

      WRITE(17,899) "# "
      WRITE(17,899) "# REDUCED CROSS SECTIONS AT LHC"
      WRITE(17,899) "BLOCK LHCCROSSSECTIONS"
      WRITE(17,901) 11,SIG(1,1),"VBF/VH -> H1 -> tautau"
      WRITE(17,901) 12,SIG(1,2),"ggF -> H1 -> tautau"
      WRITE(17,901) 13,SIG(1,3),"VBF/VH -> H1 -> bb"
      WRITE(17,901) 14,SIG(1,4),"ttH -> H1 -> bb"
      WRITE(17,901) 15,SIG(1,5),"VBF/VH -> H1 -> ZZ/WW"
      WRITE(17,901) 16,SIG(1,6),"ggF -> H1 -> ZZ/WW"
      WRITE(17,901) 17,SIG(1,7),"VBF/VH -> H1 -> gammagamma"
      WRITE(17,901) 18,SIG(1,8),"ggF -> H1 -> gammagamma"
      WRITE(17,901) 21,SIG(2,1),"VBF/VH -> H2 -> tautau"
      WRITE(17,901) 22,SIG(2,2),"ggF -> H2 -> tautau"
      WRITE(17,901) 23,SIG(2,3),"VBF/VH -> H2 -> bb"
      WRITE(17,901) 24,SIG(2,4),"ttH -> H2 -> bb"
      WRITE(17,901) 25,SIG(2,5),"VBF/VH -> H2 -> ZZ/WW"
      WRITE(17,901) 26,SIG(2,6),"ggF -> H2 -> ZZ/WW"
      WRITE(17,901) 27,SIG(2,7),"VBF/VH -> H2 -> gammagamma"
      WRITE(17,901) 28,SIG(2,8),"ggF -> H2 -> gammagamma"
      WRITE(17,901) 31,SIG(3,1),"VBF/VH -> H3 -> tautau"
      WRITE(17,901) 32,SIG(3,2),"ggF -> H3 -> tautau"
      WRITE(17,901) 33,SIG(3,3),"VBF/VH -> H3 -> bb"
      WRITE(17,901) 34,SIG(3,4),"ttH -> H3 -> bb"
      WRITE(17,901) 35,SIG(3,5),"VBF/VH -> H3 -> ZZ/WW"
      WRITE(17,901) 36,SIG(3,6),"ggF -> H3 -> ZZ/WW"
      WRITE(17,901) 37,SIG(3,7),"VBF/VH -> H3 -> gammagamma"
      WRITE(17,901) 38,SIG(3,8),"ggF -> H3 -> gammagamma"

      WRITE(17,899) "# "
      WRITE(17,899) "# CHI^2 FOR HIGGS COUPLINGS"
      WRITE(17,899) "BLOCK LHCFIT"
      WRITE(17,901) 1,chi2gam,"Hgammagamma"
      WRITE(17,901) 2,chi2bb,"Hff"
      WRITE(17,901) 3,chi2zz,"HVV"

      WRITE(17,899) "# "
      WRITE(17,899)

 1    WRITE(18,899) "# HIGGS + TOP BRANCHING RATIOS IN SLHA FORMAT"
      WRITE(18,899) "# Info about decay package"
      WRITE(18,899) "BLOCK DCINFO   # Program information"
      WRITE(18,900) 1,"NMSSMTools # Decay package"
      WRITE(18,900) 2,"4.1.1      # Version number"

      IF(IFAIL.NE.0.AND.IFAIL.NE.10) GOTO 2

      WRITE(18,899) "#           PDG          Width"
      WRITE(18,904) 25,WIDTH(1),"Lightest neutral Higgs scalar"
      IF(BRJJ(1).GT.0d0)
     .  WRITE(18,905) BRJJ(1),2,21,21,"BR(H_1 -> gluon gluon)"
      IF(BRMM(1).GT.0d0)
     .  WRITE(18,905) BRMM(1),2,13,-13,"BR(H_1 -> muon muon)"
      IF(BRLL(1).GT.0d0)
     .  WRITE(18,905) BRLL(1),2,15,-15,"BR(H_1 -> tau tau)"
      IF(BRSS(1).GT.0d0)
     .  WRITE(18,905) BRSS(1),2,3,-3,"BR(H_1 -> s sbar)"
      IF(BRCC(1).GT.0d0)
     .  WRITE(18,905) BRCC(1),2,4,-4,"BR(H_1 -> c cbar)"
      IF(BRBB(1).GT.0d0)
     .  WRITE(18,905) BRBB(1),2,5,-5,"BR(H_1 -> b bbar)"
      IF(BRTT(1).GT.0d0)
     .  WRITE(18,905) BRTT(1),2,6,-6,"BR(H_1 -> t tbar)"
      IF(BRWW(1).GT.0d0)
     .  WRITE(18,905) BRWW(1),2,24,-24,"BR(H_1 -> W+ W-)"
      IF(BRZZ(1).GT.0d0)
     .  WRITE(18,905) BRZZ(1),2,23,23,"BR(H_1 -> Z Z)"
      IF(BRGG(1).GT.0d0)
     .  WRITE(18,905) BRGG(1),2,22,22,"BR(H_1 -> gamma gamma)"
      IF(BRZG(1).GT.0d0)
     .  WRITE(18,905) BRZG(1),2,23,22,"BR(H_1 -> Z gamma)"
      IF(BRHAA(1,1).GT.0d0)
     .  WRITE(18,905) BRHAA(1,1),2,36,36,"BR(H_1 -> A_1 A_1)"
      IF(BRHAA(1,2).GT.0d0)
     .  WRITE(18,905) BRHAA(1,2),2,36,46,"BR(H_1 -> A_1 A_2)"
      IF(BRHAA(1,3).GT.0d0)
     .  WRITE(18,905) BRHAA(1,3),2,46,46,"BR(H_1 -> A_2 A_2)"
      IF(BRHAZ(1,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(1,1),2,23,36,"BR(H_1 -> A_1 Z)"
      IF(BRNEU(1,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,1),2,1000022,1000022,
     .    "BR(H_1 -> neu_1 neu_1)"
      IF(BRNEU(1,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,2),2,1000022,1000023,
     .    "BR(H_1 -> neu_1 neu_2)"
      IF(BRNEU(1,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,3),2,1000022,1000025,
     .    "BR(H_1 -> neu_1 neu_3)"
      IF(BRNEU(1,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,4),2,1000022,1000035,
     .    "BR(H_1 -> neu_1 neu_4)"
      IF(BRNEU(1,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,5),2,1000022,1000045,
     .    "BR(H_1 -> neu_1 neu_5)"
      IF(BRNEU(1,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,2),2,1000023,1000023,
     .    "BR(H_1 -> neu_2 neu_2)"
      IF(BRNEU(1,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,3),2,1000023,1000025,
     .    "BR(H_1 -> neu_2 neu_3)"
      IF(BRNEU(1,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,4),2,1000023,1000035,
     .    "BR(H_1 -> neu_2 neu_4)"
      IF(BRNEU(1,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,5),2,1000023,1000045,
     .    "BR(H_1 -> neu_2 neu_5)"
      IF(BRNEU(1,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,3),2,1000025,1000025,
     .    "BR(H_1 -> neu_3 neu_3)"
      IF(BRNEU(1,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,4),2,1000025,1000035,
     .    "BR(H_1 -> neu_3 neu_4)"
      IF(BRNEU(1,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,5),2,1000025,1000045,
     .    "BR(H_1 -> neu_3 neu_5)"
      IF(BRNEU(1,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,4,4),2,1000035,1000035,
     .    "BR(H_1 -> neu_4 neu_4)"
      IF(BRNEU(1,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,4,5),2,1000035,1000045,
     .    "BR(H_1 -> neu_4 neu_5)"
      IF(BRNEU(1,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,5,5),2,1000045,1000045,
     .    "BR(H_1 -> neu_5 neu_5)"
      IF(BRCHA(1,1).GT.0d0)
     .  WRITE(18,905) BRCHA(1,1),2,1000024,-1000024,
     .    "BR(H_1 -> cha_1 cha_1bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(18,905) BRCHA(1,2),2,1000024,-1000037,
     .    "BR(H_1 -> cha_1 cha_2bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(18,905) BRCHA(1,2),2,1000037,-1000024,
     .    "BR(H_1 -> cha_2 cha_1bar)"
      IF(BRCHA(1,3).GT.0d0)
     .  WRITE(18,905) BRCHA(1,3),2,1000037,-1000037,
     .    "BR(H_1 -> cha_2 cha_2bar)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,1),2,1000002,-1000002,
     .    "BR(H_1 -> ~u_L ~ubar_L)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,1),2,1000004,-1000004,
     .    "BR(H_1 -> ~c_L ~cbar_L)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,2),2,2000002,-2000002,
     .    "BR(H_1 -> ~u_R ~ubar_R)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,2),2,2000004,-2000004,
     .    "BR(H_1 -> ~c_R ~cbar_R)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,3),2,1000001,-1000001,
     .    "BR(H_1 -> ~d_L ~dbar_L)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,3),2,1000003,-1000003,
     .    "BR(H_1 -> ~s_L ~sbar_L)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,4),2,2000001,-2000001,
     .    "BR(H_1 -> ~d_R ~dbar_R)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,4),2,2000003,-2000003,
     .    "BR(H_1 -> ~s_R ~sbar_R)"
      IF(BRHSQ(1,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,5),2,1000006,-1000006,
     .    "BR(H_1 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(1,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,6),2,2000006,-2000006,
     .    "BR(H_1 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,7),2,1000006,-2000006,
     .    "BR(H_1 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,7),2,2000006,-1000006,
     .    "BR(H_1 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(1,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,8),2,1000005,-1000005,
     .    "BR(H_1 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(1,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,9),2,2000005,-2000005,
     .    "BR(H_1 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,10),2,1000005,-2000005,
     .    "BR(H_1 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,10),2,2000005,-1000005,
     .    "BR(H_1 -> ~b_2 ~bbar_1)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(18,905) BRHSL(1,1),2,1000011,-1000011,
     .    "BR(H_1 -> ~e_L ~ebar_L)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(18,905) BRHSL(1,1),2,1000013,-1000013,
     .    "BR(H_1 -> ~mu_L ~mubar_L)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(18,905) BRHSL(1,2),2,2000011,-2000011,
     .    "BR(H_1 -> ~e_R ~ebar_R)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(18,905) BRHSL(1,2),2,2000013,-2000013,
     .    "BR(H_1 -> ~mu_R ~mubarRL)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(18,905) BRHSL(1,3),2,1000012,-1000012,
     .    "BR(H_1 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(18,905) BRHSL(1,3),2,1000014,-1000014,
     .    "BR(H_1 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(1,4).GT.0d0)
     .  WRITE(18,905) BRHSL(1,4),2,1000015,-1000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(1,5).GT.0d0)
     .  WRITE(18,905) BRHSL(1,5),2,2000015,-2000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(18,905) BRHSL(1,6),2,1000015,-2000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(18,905) BRHSL(1,6),2,2000015,-1000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(1,7).GT.0d0)
     .  WRITE(18,905) BRHSL(1,7),2,1000016,-1000016,
     .    "BR(H_1 -> ~nu_tau_L ~nu_taubar_L)"
      
      WRITE(18,904) 35,WIDTH(2),"2nd neutral Higgs scalar"
      IF(BRJJ(2).GT.0d0)
     .  WRITE(18,905) BRJJ(2),2,21,21,"BR(H_2 -> gluon gluon)"
      IF(BRMM(2).GT.0d0)
     .  WRITE(18,905) BRMM(2),2,13,-13,"BR(H_2 -> muon muon)"
      IF(BRLL(2).GT.0d0)
     .  WRITE(18,905) BRLL(2),2,15,-15,"BR(H_2 -> tau tau)"
      IF(BRSS(2).GT.0d0)
     .  WRITE(18,905) BRSS(2),2,3,-3,"BR(H_2 -> s sbar)"
      IF(BRCC(2).GT.0d0)
     .  WRITE(18,905) BRCC(2),2,4,-4,"BR(H_2 -> c cbar)"
      IF(BRBB(2).GT.0d0)
     .  WRITE(18,905) BRBB(2),2,5,-5,"BR(H_2 -> b bbar)"
      IF(BRTT(2).GT.0d0)
     .  WRITE(18,905) BRTT(2),2,6,-6,"BR(H_2 -> t tbar)"
      IF(BRWW(2).GT.0d0)
     .  WRITE(18,905) BRWW(2),2,24,-24,"BR(H_2 -> W+ W-)"
      IF(BRZZ(2).GT.0d0)
     .  WRITE(18,905) BRZZ(2),2,23,23,"BR(H_2 -> Z Z)"
      IF(BRGG(2).GT.0d0)
     .  WRITE(18,905) BRGG(2),2,22,22,"BR(H_2 -> gamma gamma)"
      IF(BRZG(2).GT.0d0)
     .  WRITE(18,905) BRZG(2),2,23,22,"BR(H_2 -> Z gamma)"
      IF(BRHHH(1).GT.0d0)
     .  WRITE(18,905) BRHHH(1),2,25,25,"BR(H_2 -> H_1 H_1)"
      IF(BRHAA(2,1).GT.0d0)
     .  WRITE(18,905) BRHAA(2,1),2,36,36,"BR(H_2 -> A_1 A_1)"
      IF(BRHAA(2,2).GT.0d0)
     .  WRITE(18,905) BRHAA(2,2),2,36,46,"BR(H_2 -> A_1 A_2)"
      IF(BRHAA(2,3).GT.0d0)
     .  WRITE(18,905) BRHAA(2,3),2,46,46,"BR(H_2 -> A_2 A_2)"
      IF(BRHAZ(2,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(2,1),2,23,36,"BR(H_2 -> A_1 Z)"
      IF(BRNEU(2,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,1),2,1000022,1000022,
     .    "BR(H_2 -> neu_1 neu_1)"
      IF(BRNEU(2,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,2),2,1000022,1000023,
     .    "BR(H_2 -> neu_1 neu_2)"
      IF(BRNEU(2,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,3),2,1000022,1000025,
     .    "BR(H_2 -> neu_1 neu_3)"
      IF(BRNEU(2,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,4),2,1000022,1000035,
     .    "BR(H_2 -> neu_1 neu_4)"
      IF(BRNEU(2,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,5),2,1000022,1000045,
     .    "BR(H_2 -> neu_1 neu_5)"
      IF(BRNEU(2,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,2),2,1000023,1000023,
     .    "BR(H_2 -> neu_2 neu_2)"
      IF(BRNEU(2,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,3),2,1000023,1000025,
     .    "BR(H_2 -> neu_2 neu_3)"
      IF(BRNEU(2,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,4),2,1000023,1000035,
     .    "BR(H_2 -> neu_2 neu_4)"
      IF(BRNEU(2,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,5),2,1000023,1000045,
     .    "BR(H_2 -> neu_2 neu_5)"
      IF(BRNEU(2,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,3),2,1000025,1000025,
     .    "BR(H_2 -> neu_3 neu_3)"
      IF(BRNEU(2,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,4),2,1000025,1000035,
     .    "BR(H_2 -> neu_3 neu_4)"
      IF(BRNEU(2,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,5),2,1000025,1000045,
     .    "BR(H_2 -> neu_3 neu_5)"
      IF(BRNEU(2,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,4,4),2,1000035,1000035,
     .    "BR(H_2 -> neu_4 neu_4)"
      IF(BRNEU(2,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,4,5),2,1000035,1000045,
     .    "BR(H_2 -> neu_4 neu_5)"
      IF(BRNEU(2,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,5,5),2,1000045,1000045,
     .    "BR(H_2 -> neu_5 neu_5)"
      IF(BRCHA(2,1).GT.0d0)
     .  WRITE(18,905) BRCHA(2,1),2,1000024,-1000024,
     .    "BR(H_2 -> cha_1 cha_1bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(18,905) BRCHA(2,2),2,1000024,-1000037,
     .    "BR(H_2 -> cha_1 cha_2bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(18,905) BRCHA(2,2),2,1000037,-1000024,
     .    "BR(H_2 -> cha_2 cha_1bar)"
      IF(BRCHA(2,3).GT.0d0)
     .  WRITE(18,905) BRCHA(2,3),2,1000037,-1000037,
     .    "BR(H_2 -> cha_2 cha_2bar)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,1),2,1000002,-1000002,
     .    "BR(H_2 -> ~u_L ~ubar_L)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,1),2,1000004,-1000004,
     .    "BR(H_2 -> ~c_L ~cbar_L)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,2),2,2000002,-2000002,
     .    "BR(H_2 -> ~u_R ~ubar_R)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,2),2,2000004,-2000004,
     .    "BR(H_2 -> ~c_R ~cbar_R)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,3),2,1000001,-1000001,
     .    "BR(H_2 -> ~d_L ~dbar_L)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,3),2,1000003,-1000003,
     .    "BR(H_2 -> ~s_L ~sbar_L)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,4),2,2000001,-2000001,
     .    "BR(H_2 -> ~d_R ~dbar_R)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,4),2,2000003,-2000003,
     .    "BR(H_2 -> ~s_R ~sbar_R)"
      IF(BRHSQ(2,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,5),2,1000006,-1000006,
     .    "BR(H_2 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(2,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,6),2,2000006,-2000006,
     .    "BR(H_2 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,7),2,1000006,-2000006,
     .    "BR(H_2 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,7),2,2000006,-1000006,
     .    "BR(H_2 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(2,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,8),2,1000005,-1000005,
     .    "BR(H_2 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(2,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,9),2,2000005,-2000005,
     .    "BR(H_2 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,10),2,1000005,-2000005,
     .    "BR(H_2 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,10),2,2000005,-1000005,
     .    "BR(H_2 -> ~b_2 ~bbar_1)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(18,905) BRHSL(2,1),2,1000011,-1000011,
     .    "BR(H_2 -> ~e_L ~ebar_L)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(18,905) BRHSL(2,1),2,1000013,-1000013,
     .    "BR(H_2 -> ~mu_L ~mubar_L)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(18,905) BRHSL(2,2),2,2000011,-2000011,
     .    "BR(H_2 -> ~e_R ~ebar_R)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(18,905) BRHSL(2,2),2,2000013,-2000013,
     .    "BR(H_2 -> ~mu_R ~mubarRL)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(18,905) BRHSL(2,3),2,1000012,-1000012,
     .    "BR(H_2 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(18,905) BRHSL(2,3),2,1000014,-1000014,
     .    "BR(H_2 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(2,4).GT.0d0)
     .  WRITE(18,905) BRHSL(2,4),2,1000015,-1000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(2,5).GT.0d0)
     .  WRITE(18,905) BRHSL(2,5),2,2000015,-2000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(18,905) BRHSL(2,6),2,1000015,-2000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(18,905) BRHSL(2,6),2,2000015,-1000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(2,7).GT.0d0)
     .  WRITE(18,905) BRHSL(2,7),2,1000016,-1000016,
     .    "BR(H_2 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 45,WIDTH(3),"3rd neutral Higgs scalar"
      IF(BRJJ(3).GT.0d0)
     .  WRITE(18,905) BRJJ(3),2,21,21,"BR(H_3 -> gluon gluon)"
      IF(BRMM(3).GT.0d0)
     .  WRITE(18,905) BRMM(3),2,13,-13,"BR(H_3 -> muon muon)"
      IF(BRLL(3).GT.0d0)
     .  WRITE(18,905) BRLL(3),2,15,-15,"BR(H_3 -> tau tau)"
      IF(BRSS(3).GT.0d0)
     .  WRITE(18,905) BRSS(3),2,3,-3,"BR(H_3 -> s sbar)"
      IF(BRCC(3).GT.0d0)
     .  WRITE(18,905) BRCC(3),2,4,-4,"BR(H_3 -> c cbar)"
      IF(BRBB(3).GT.0d0)
     .  WRITE(18,905) BRBB(3),2,5,-5,"BR(H_3 -> b bbar)"
      IF(BRTT(3).GT.0d0)
     .  WRITE(18,905) BRTT(3),2,6,-6,"BR(H_3 -> t tbar)"
      IF(BRWW(3).GT.0d0)
     .  WRITE(18,905) BRWW(3),2,24,-24,"BR(H_3 -> W+ W-)"
      IF(BRZZ(3).GT.0d0)
     .  WRITE(18,905) BRZZ(3),2,23,23,"BR(H_3 -> Z Z)"
      IF(BRGG(3).GT.0d0)
     .  WRITE(18,905) BRGG(3),2,22,22,"BR(H_3 -> gamma gamma)"
      IF(BRZG(3).GT.0d0)
     .  WRITE(18,905) BRZG(3),2,23,22,"BR(H_3 -> Z gamma)"
      IF(BRHHH(2).GT.0d0)
     .  WRITE(18,905) BRHHH(2),2,25,25,"BR(H_3 -> H_1 H_1)"
      IF(BRHHH(3).GT.0d0)
     .  WRITE(18,905) BRHHH(3),2,25,35,"BR(H_3 -> H_1 H_2)"
      IF(BRHHH(4).GT.0d0)
     .  WRITE(18,905) BRHHH(4),2,35,35,"BR(H_3 -> H_2 H_2)"
      IF(BRHAA(3,1).GT.0d0)
     .  WRITE(18,905) BRHAA(3,1),2,36,36,"BR(H_3 -> A_1 A_1)"
      IF(BRHAA(3,2).GT.0d0)
     .  WRITE(18,905) BRHAA(3,2),2,36,46,"BR(H_3 -> A_1 A_2)"
      IF(BRHAA(3,3).GT.0d0)
     .  WRITE(18,905) BRHAA(3,3),2,46,46,"BR(H_3 -> A_2 A_2)"
      IF(BRHAZ(3,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(3,1),2,23,36,"BR(H_3 -> A_1 Z)"
      IF(BRHAZ(3,2).GT.0d0)
     .  WRITE(18,905) BRHAZ(3,2),2,23,46,"BR(H_3 -> A_2 Z)"
      IF(BRHCHC(3).GT.0d0)
     .  WRITE(18,905) BRHCHC(3),2,37,-37,"BR(H_3 -> H+ H-)"
      IF(BRHCW(3).GT.0d0)
     .  WRITE(18,905) BRHCW(3),2,24,-37,"BR(H_3 -> W+ H-)"
      IF(BRNEU(3,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,1),2,1000022,1000022,
     .    "BR(H_3 -> neu_1 neu_1)"
      IF(BRNEU(3,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,2),2,1000022,1000023,
     .    "BR(H_3 -> neu_1 neu_2)"
      IF(BRNEU(3,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,3),2,1000022,1000025,
     .    "BR(H_3 -> neu_1 neu_3)"
      IF(BRNEU(3,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,4),2,1000022,1000035,
     .    "BR(H_3 -> neu_1 neu_4)"
      IF(BRNEU(3,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,5),2,1000022,1000045,
     .    "BR(H_3 -> neu_1 neu_5)"
      IF(BRNEU(3,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,2),2,1000023,1000023,
     .    "BR(H_3 -> neu_2 neu_2)"
      IF(BRNEU(3,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,3),2,1000023,1000025,
     .    "BR(H_3 -> neu_2 neu_3)"
      IF(BRNEU(3,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,4),2,1000023,1000035,
     .    "BR(H_3 -> neu_2 neu_4)"
      IF(BRNEU(3,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,5),2,1000023,1000045,
     .    "BR(H_3 -> neu_2 neu_5)"
      IF(BRNEU(3,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,3),2,1000025,1000025,
     .    "BR(H_3 -> neu_3 neu_3)"
      IF(BRNEU(3,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,4),2,1000025,1000035,
     .    "BR(H_3 -> neu_3 neu_4)"
      IF(BRNEU(3,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,5),2,1000025,1000045,
     .    "BR(H_3 -> neu_3 neu_5)"
      IF(BRNEU(3,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,4,4),2,1000035,1000035,
     .    "BR(H_3 -> neu_4 neu_4)"
      IF(BRNEU(3,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,4,5),2,1000035,1000045,
     .    "BR(H_3 -> neu_4 neu_5)"
      IF(BRNEU(3,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,5,5),2,1000045,1000045,
     .    "BR(H_3 -> neu_5 neu_5)"
      IF(BRCHA(3,1).GT.0d0)
     .  WRITE(18,905) BRCHA(3,1),2,1000024,-1000024,
     .    "BR(H_3 -> cha_1 cha_1bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(18,905) BRCHA(3,2),2,1000024,-1000037,
     .    "BR(H_3 -> cha_1 cha_2bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(18,905) BRCHA(3,2),2,1000037,-1000024,
     .    "BR(H_3 -> cha_2 cha_1bar)"
      IF(BRCHA(3,3).GT.0d0)
     .  WRITE(18,905) BRCHA(3,3),2,1000037,-1000037,
     .    "BR(H_3 -> cha_2 cha_2bar)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,1),2,1000002,-1000002,
     .    "BR(H_3 -> ~u_L ~ubar_L)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,1),2,1000004,-1000004,
     .    "BR(H_3 -> ~c_L ~cbar_L)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,2),2,2000002,-2000002,
     .    "BR(H_3 -> ~u_R ~ubar_R)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,2),2,2000004,-2000004,
     .    "BR(H_3 -> ~c_R ~cbar_R)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,3),2,1000001,-1000001,
     .    "BR(H_3 -> ~d_L ~dbar_L)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,3),2,1000003,-1000003,
     .    "BR(H_3 -> ~s_L ~sbar_L)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,4),2,2000001,-2000001,
     .    "BR(H_3 -> ~d_R ~dbar_R)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,4),2,2000003,-2000003,
     .    "BR(H_3 -> ~s_R ~sbar_R)"
      IF(BRHSQ(3,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,5),2,1000006,-1000006,
     .    "BR(H_3 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(3,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,6),2,2000006,-2000006,
     .    "BR(H_3 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,7),2,1000006,-2000006,
     .    "BR(H_3 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,7),2,2000006,-1000006,
     .    "BR(H_3 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(3,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,8),2,1000005,-1000005,
     .    "BR(H_3 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(3,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,9),2,2000005,-2000005,
     .    "BR(H_3 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,10),2,1000005,-2000005,
     .    "BR(H_3 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,10),2,2000005,-1000005,
     .    "BR(H_3 -> ~b_2 ~bbar_1)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(18,905) BRHSL(3,1),2,1000011,-1000011,
     .    "BR(H_3 -> ~e_L ~ebar_L)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(18,905) BRHSL(3,1),2,1000013,-1000013,
     .    "BR(H_3 -> ~mu_L ~mubar_L)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(18,905) BRHSL(3,2),2,2000011,-2000011,
     .    "BR(H_3 -> ~e_R ~ebar_R)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(18,905) BRHSL(3,2),2,2000013,-2000013,
     .    "BR(H_3 -> ~mu_R ~mubarRL)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(18,905) BRHSL(3,3),2,1000012,-1000012,
     .    "BR(H_3 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(18,905) BRHSL(3,3),2,1000014,-1000014,
     .    "BR(H_3 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(3,4).GT.0d0)
     .  WRITE(18,905) BRHSL(3,4),2,1000015,-1000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(3,5).GT.0d0)
     .  WRITE(18,905) BRHSL(3,5),2,2000015,-2000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(18,905) BRHSL(3,6),2,1000015,-2000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(18,905) BRHSL(3,6),2,2000015,-1000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(3,7).GT.0d0)
     .  WRITE(18,905) BRHSL(3,7),2,1000016,-1000016,
     .    "BR(H_3 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 36,WIDTH(4),"Lightest pseudoscalar"
      IF(BRJJ(4).GT.0d0)
     .  WRITE(18,905) BRJJ(4),2,21,21,"BR(A_1 -> gluon gluon)"
      IF(BRMM(4).GT.0d0)
     .  WRITE(18,905) BRMM(4),2,13,-13,"BR(A_1 -> muon muon)"
      IF(BRLL(4).GT.0d0)
     .  WRITE(18,905) BRLL(4),2,15,-15,"BR(A_1 -> tau tau)"
      IF(BRSS(4).GT.0d0)
     .  WRITE(18,905) BRSS(4),2,3,-3,"BR(A_1 -> s sbar)"
      IF(BRCC(4).GT.0d0)
     .  WRITE(18,905) BRCC(4),2,4,-4,"BR(A_1 -> c cbar)"
      IF(BRBB(4).GT.0d0)
     .  WRITE(18,905) BRBB(4),2,5,-5,"BR(A_1 -> b bbar)"
      IF(BRTT(4).GT.0d0)
     .  WRITE(18,905) BRTT(4),2,6,-6,"BR(A_1 -> t tbar)"
      IF(BRGG(4).GT.0d0)
     .  WRITE(18,905) BRGG(4),2,22,22,"BR(A_1 -> gamma gamma)"
      IF(BRZG(4).GT.0d0)
     .  WRITE(18,905) BRZG(4),2,23,22,"BR(A_1 -> Z gamma)"
      IF(BRAHZ(1,1).GT.0d0)
     .  WRITE(18,905) BRAHZ(1,1),2,23,25,"BR(A_1 -> Z H_1)"
      IF(BRAHZ(1,2).GT.0d0)
     .  WRITE(18,905) BRAHZ(1,2),2,23,35,"BR(A_1 -> Z H_2)"
      IF(BRNEU(4,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,1),2,1000022,1000022,
     .    "BR(A_1 -> neu_1 neu_1)"
      IF(BRNEU(4,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,2),2,1000022,1000023,
     .    "BR(A_1 -> neu_1 neu_2)"
      IF(BRNEU(4,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,3),2,1000022,1000025,
     .    "BR(A_1 -> neu_1 neu_3)"
      IF(BRNEU(4,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,4),2,1000022,1000035,
     .    "BR(A_1 -> neu_1 neu_4)"
      IF(BRNEU(4,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,5),2,1000022,1000045,
     .    "BR(A_1 -> neu_1 neu_5)"
      IF(BRNEU(4,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,2),2,1000023,1000023,
     .    "BR(A_1 -> neu_2 neu_2)"
      IF(BRNEU(4,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,3),2,1000023,1000025,
     .    "BR(A_1 -> neu_2 neu_3)"
      IF(BRNEU(4,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,4),2,1000023,1000035,
     .    "BR(A_1 -> neu_2 neu_4)"
      IF(BRNEU(4,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,5),2,1000023,1000045,
     .    "BR(A_1 -> neu_2 neu_5)"
      IF(BRNEU(4,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,3),2,1000025,1000025,
     .    "BR(A_1 -> neu_3 neu_3)"
      IF(BRNEU(4,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,4),2,1000025,1000035,
     .    "BR(A_1 -> neu_3 neu_4)"
      IF(BRNEU(4,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,5),2,1000025,1000045,
     .    "BR(A_1 -> neu_3 neu_5)"
      IF(BRNEU(4,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,4,4),2,1000035,1000035,
     .    "BR(A_1 -> neu_4 neu_4)"
      IF(BRNEU(4,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,4,5),2,1000035,1000045,
     .    "BR(A_1 -> neu_4 neu_5)"
      IF(BRNEU(4,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,5,5),2,1000045,1000045,
     .    "BR(A_1 -> neu_5 neu_5)"
      IF(BRCHA(4,1).GT.0d0)
     .  WRITE(18,905) BRCHA(4,1),2,1000024,-1000024,
     .    "BR(A_1 -> cha_1 cha_1bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(18,905) BRCHA(4,2),2,1000024,-1000037,
     .    "BR(A_1 -> cha_1 cha_2bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(18,905) BRCHA(4,2),2,1000037,-1000024,
     .    "BR(A_1 -> cha_2 cha_1bar)"
      IF(BRCHA(4,3).GT.0d0)
     .  WRITE(18,905) BRCHA(4,3),2,1000037,-1000037,
     .    "BR(A_1 -> cha_2 cha_2bar)"
      IF(BRASQ(1,1).GT.0d0)
     .  WRITE(18,905) BRASQ(1,1),2,1000006,-2000006,
     .    "BR(A_1 -> ~t_1 ~tbar_2)"
      IF(BRASQ(1,1).GT.0d0)
     .  WRITE(18,905) BRASQ(1,1),2,2000006,-1000006,
     .    "BR(A_1 -> ~t_2 ~tbar_1)"
      IF(BRASQ(1,2).GT.0d0)
     .  WRITE(18,905) BRASQ(1,2),2,1000005,-2000005,
     .    "BR(A_1 -> ~b_1 ~bbar_2)"
      IF(BRASQ(1,2).GT.0d0)
     .  WRITE(18,905) BRASQ(1,2),2,2000005,-1000005,
     .    "BR(A_1 -> ~b_2 ~bbar_1)"
      IF(BRASL(1).GT.0d0)
     .  WRITE(18,905) BRASL(1),2,1000015,-2000015,
     .    "BR(A_1 -> ~tau_1 ~taubar_2)"
      IF(BRASL(1).GT.0d0)
     .  WRITE(18,905) BRASL(1),2,2000015,-1000015,
     .    "BR(A_1 -> ~tau_2 ~taubar_1)"

      WRITE(18,904) 46,WIDTH(5),"2nd pseudoscalar"
      IF(BRJJ(5).GT.0d0)
     .  WRITE(18,905) BRJJ(5),2,21,21,"BR(A_2 -> gluon gluon)"
      IF(BRMM(5).GT.0d0)
     .  WRITE(18,905) BRMM(5),2,13,-13,"BR(A_2 -> muon muon)"
      IF(BRLL(5).GT.0d0)
     .  WRITE(18,905) BRLL(5),2,15,-15,"BR(A_2 -> tau tau)"
      IF(BRSS(5).GT.0d0)
     .  WRITE(18,905) BRSS(5),2,3,-3,"BR(A_2 -> s sbar)"
      IF(BRCC(5).GT.0d0)
     .  WRITE(18,905) BRCC(5),2,4,-4,"BR(A_2 -> c cbar)"
      IF(BRBB(5).GT.0d0)
     .  WRITE(18,905) BRBB(5),2,5,-5,"BR(A_2 -> b bbar)"
      IF(BRTT(5).GT.0d0)
     .  WRITE(18,905) BRTT(5),2,6,-6,"BR(A_2 -> t tbar)"
      IF(BRGG(5).GT.0d0)
     .  WRITE(18,905) BRGG(5),2,22,22,"BR(A_2 -> gamma gamma)"
      IF(BRZG(5).GT.0d0)
     .  WRITE(18,905) BRZG(5),2,23,22,"BR(A_2 -> Z gamma)"
      IF(BRAHA(1).GT.0d0)
     .  WRITE(18,905) BRAHA(1),2,36,25,"BR(A_2 -> A_1 H_1)"
      IF(BRAHA(2).GT.0d0)
     .  WRITE(18,905) BRAHA(2),2,36,35,"BR(A_2 -> A_1 H_2)"
      IF(BRAHA(3).GT.0d0)
     .  WRITE(18,905) BRAHA(3),2,36,45,"BR(A_2 -> A_1 H_3)"
      IF(BRAHZ(2,1).GT.0d0)
     .  WRITE(18,905) BRAHZ(2,1),2,23,25,"BR(A_2 -> Z H_1)"
      IF(BRAHZ(2,2).GT.0d0)
     .  WRITE(18,905) BRAHZ(2,2),2,23,35,"BR(A_2 -> Z H_2)"
      IF(BRHCW(5).GT.0d0)
     .  WRITE(18,905) BRHCW(5),2,24,-37,"BR(A_2 -> W+ H-)"
      IF(BRNEU(5,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,1),2,1000022,1000022,
     .    "BR(A_2 -> neu_1 neu_1)"
      IF(BRNEU(5,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,2),2,1000022,1000023,
     .    "BR(A_2 -> neu_1 neu_2)"
      IF(BRNEU(5,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,3),2,1000022,1000025,
     .    "BR(A_2 -> neu_1 neu_3)"
      IF(BRNEU(5,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,4),2,1000022,1000035,
     .    "BR(A_2 -> neu_1 neu_4)"
      IF(BRNEU(5,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,5),2,1000022,1000045,
     .    "BR(A_2 -> neu_1 neu_5)"
      IF(BRNEU(5,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,2),2,1000023,1000023,
     .    "BR(A_2 -> neu_2 neu_2)"
      IF(BRNEU(5,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,3),2,1000023,1000025,
     .    "BR(A_2 -> neu_2 neu_3)"
      IF(BRNEU(5,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,4),2,1000023,1000035,
     .    "BR(A_2 -> neu_2 neu_4)"
      IF(BRNEU(5,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,5),2,1000023,1000045,
     .    "BR(A_2 -> neu_2 neu_5)"
      IF(BRNEU(5,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,3),2,1000025,1000025,
     .    "BR(A_2 -> neu_3 neu_3)"
      IF(BRNEU(5,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,4),2,1000025,1000035,
     .    "BR(A_2 -> neu_3 neu_4)"
      IF(BRNEU(5,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,5),2,1000025,1000045,
     .    "BR(A_2 -> neu_3 neu_5)"
      IF(BRNEU(5,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,4,4),2,1000035,1000035,
     .    "BR(A_2 -> neu_4 neu_4)"
      IF(BRNEU(5,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,4,5),2,1000035,1000045,
     .    "BR(A_2 -> neu_4 neu_5)"
      IF(BRNEU(5,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,5,5),2,1000045,1000045,
     .    "BR(A_2 -> neu_5 neu_5)"
      IF(BRCHA(5,1).GT.0d0)
     .  WRITE(18,905) BRCHA(5,1),2,1000024,-1000024,
     .    "BR(A_2 -> cha_1 cha_1bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(18,905) BRCHA(5,2),2,1000024,-1000037,
     .    "BR(A_2 -> cha_1 cha_2bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(18,905) BRCHA(5,2),2,1000037,-1000024,
     .    "BR(A_2 -> cha_2 cha_1bar)"
      IF(BRCHA(5,3).GT.0d0)
     .  WRITE(18,905) BRCHA(5,3),2,1000037,-1000037,
     .    "BR(A_2 -> cha_2 cha_2bar)"
      IF(BRASQ(2,1).GT.0d0)
     .  WRITE(18,905) BRASQ(2,1),2,1000006,-2000006,
     .    "BR(A_2 -> ~t_1 ~tbar_2)"
      IF(BRASQ(2,1).GT.0d0)
     .  WRITE(18,905) BRASQ(2,1),2,2000006,-1000006,
     .    "BR(A_2 -> ~t_2 ~tbar_1)"
      IF(BRASQ(2,2).GT.0d0)
     .  WRITE(18,905) BRASQ(2,2),2,1000005,-2000005,
     .    "BR(A_2 -> ~b_1 ~bbar_2)"
      IF(BRASQ(2,2).GT.0d0)
     .  WRITE(18,905) BRASQ(2,2),2,2000005,-1000005,
     .    "BR(A_2 -> ~b_2 ~bbar_1)"
      IF(BRASL(2).GT.0d0)
     .  WRITE(18,905) BRASL(2),2,1000015,-2000015,
     .    "BR(A_2 -> ~tau_1 ~taubar_2)"
      IF(BRASL(2).GT.0d0)
     .  WRITE(18,905) BRASL(2),2,2000015,-1000015,
     .    "BR(A_2 -> ~tau_2 ~taubar_1)"

      WRITE(18,904) 37,HCWIDTH,"Charged Higgs"
      IF(HCBRM.GT.0d0)
     .  WRITE(18,905) HCBRM,2,-13,14,"BR(H+ -> muon nu_muon)"
      IF(HCBRL.GT.0d0)
     .  WRITE(18,905) HCBRL,2,-15,16,"BR(H+ -> tau nu_tau)"
      IF(HCBRSU.GT.0d0)
     .  WRITE(18,905) HCBRSU,2,2,-3,"BR(H+ -> u sbar)"
      IF(HCBRSC.GT.0d0)
     .  WRITE(18,905) HCBRSC,2,4,-3,"BR(H+ -> c sbar)"
      IF(HCBRBU.GT.0d0)
     .  WRITE(18,905) HCBRBU,2,2,-5,"BR(H+ -> u bbar)"
      IF(HCBRBC.GT.0d0)
     .  WRITE(18,905) HCBRBC,2,4,-5,"BR(H+ -> c bbar)"
      IF(HCBRBT.GT.0d0)
     .  WRITE(18,905) HCBRBT,2,6,-5,"BR(H+ -> t bbar)"
      IF(HCBRWH(1).GT.0d0)
     .  WRITE(18,905) HCBRWH(1),2,24,25,"BR(H+ -> W+ H_1)"
      IF(HCBRWH(2).GT.0d0)
     .  WRITE(18,905) HCBRWH(2),2,24,35,"BR(H+ -> W+ H_2)"
      IF(HCBRWH(4).GT.0d0)
     .  WRITE(18,905) HCBRWH(4),2,24,36,"BR(H+ -> W+ A_1)"
      IF(HCBRNC(1,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(1,1),2,1000024,1000022,
     .    "BR(H+ -> cha_1 neu_1)"
      IF(HCBRNC(2,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(2,1),2,1000024,1000023,
     .    "BR(H+ -> cha_1 neu_2)"
      IF(HCBRNC(3,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(3,1),2,1000024,1000025,
     .    "BR(H+ -> cha_1 neu_3)"
      IF(HCBRNC(4,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(4,1),2,1000024,1000035,
     .    "BR(H+ -> cha_1 neu_4)"
      IF(HCBRNC(5,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(5,1),2,1000024,1000045,
     .    "BR(H+ -> cha_1 neu_5)"
      IF(HCBRNC(1,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(1,2),2,1000037,1000022,
     .    "BR(H+ -> cha_2 neu_1)"
      IF(HCBRNC(2,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(2,2),2,1000037,1000023,
     .    "BR(H+ -> cha_2 neu_2)"
      IF(HCBRNC(3,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(3,2),2,1000037,1000025,
     .    "BR(H+ -> cha_2 neu_3)"
      IF(HCBRNC(4,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(4,2),2,1000037,1000035,
     .    "BR(H+ -> cha_2 neu_4)"
      IF(HCBRNC(5,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(5,2),2,1000037,1000045,
     .    "BR(H+ -> cha_2 neu_5)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(18,905) HCBRSQ(1),2,1000002,-1000001,
     .    "BR(H+ -> ~u_L ~dbar_L)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(18,905) HCBRSQ(1),2,1000004,-1000003,
     .    "BR(H+ -> ~c_L ~sbar_L)"
      IF(HCBRSQ(2).GT.0d0)
     .  WRITE(18,905) HCBRSQ(2),2,1000006,-1000005,
     .    "BR(H+ -> ~t_1 ~bbar_1)"
      IF(HCBRSQ(3).GT.0d0)
     .  WRITE(18,905) HCBRSQ(3),2,1000006,-2000005,
     .    "BR(H+ -> ~t_1 ~bbar_2)"
      IF(HCBRSQ(4).GT.0d0)
     .  WRITE(18,905) HCBRSQ(4),2,2000006,-1000005,
     .    "BR(H+ -> ~t_2 ~bbar_1)"
      IF(HCBRSQ(5).GT.0d0)
     .  WRITE(18,905) HCBRSQ(5),2,2000006,-2000005,
     .    "BR(H+ -> ~t_2 ~bbar_2)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(18,905) HCBRSL(1),2,1000012,-1000011,
     .    "BR(H+ -> ~nu_e_L ~ebar_L)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(18,905) HCBRSL(1),2,1000014,-1000013,
     .    "BR(H+ -> ~nu_mu_L ~mubar_L)"
      IF(HCBRSL(2).GT.0d0)
     .  WRITE(18,905) HCBRSL(2),2,1000016,-1000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_1)"
      IF(HCBRSL(3).GT.0d0)
     .  WRITE(18,905) HCBRSL(3),2,1000016,-2000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_2)"

      WRITE(18,904) 6,toptot,'Top Quark'
      IF(brtopbw.ne.0.D0)
     .  WRITE(18,905) brtopbw,2,5,24,'BR(t ->  b    W+)'
      IF(brtopbh.ne.0.D0)
     .  WRITE(18,905) brtopbh,2,5,37,'BR(t ->  b    H+)'
      IF(brtopneutrstop(1,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(1,1),2,1000006,1000022,
     . 'BR(t -> ~t_1 ~chi_10)'
      IF(brtopneutrstop(2,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(2,1),2,1000006,1000023,
     . 'BR(t -> ~t_1 ~chi_20)'
      IF(brtopneutrstop(3,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(3,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_30)'
      IF(brtopneutrstop(4,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(4,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_40)'
      IF(brtopneutrstop(5,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(5,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_50)'
      IF(brtopneutrstop(1,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(1,2),2,2000006,1000022,
     . 'BR(t -> ~t_2 ~chi_10)'
      IF(brtopneutrstop(2,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(2,2),2,2000006,1000023,
     . 'BR(t -> ~t_2 ~chi_20)'
      IF(brtopneutrstop(3,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(3,2),2,2000006,1000025,
     .'BR(t -> ~t_2 ~chi_30)'
      IF(brtopneutrstop(4,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(4,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_40)'
      IF(brtopneutrstop(5,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(5,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_50)'

      IF(NMSFLAG.NE.0)CALL NS_OUTPUT

 2    IF(OMGFLAG.EQ.0) RETURN
      WRITE(19,899) "# RELIC DENSITY CALCULATED BY MICROMEGAS"
      WRITE(19,899) "#"
      WRITE(19,899) "BLOCK RDINFO   # Program information"
      WRITE(19,900) 1,"MicrOmegas # Dark matter package"
      WRITE(19,900) 2,"4.1.1      # Version number"
      IF(PROB(30).GT.0d0)WRITE(19,900) 3,
     . "# Relic density too large (WMAP)"
      IF(PROB(30).LT.0d0.AND.PROB(30).NE.-1d0.AND.PROB(30).NE.-2d0)
     . WRITE(19,900) 3,"# Relic density too small (WMAP)"
      IF(PROB(31).NE.0d0)WRITE(19,900) 3,
     . "# Excluded by LUX"
      IF(IFAIL.EQ.0.OR.IFAIL.GE.10)THEN
        IF(OMGFLAG.NE.0)CALL printRelDen(PROB,19) 
      ELSE
        WRITE(19,900) 4,"# Cannot compute Omega h^2 (0<IFAIL<10)"
      ENDIF

 899  FORMAT(A)
 900  FORMAT(1X,I5,3X,A)
 901  FORMAT(1X,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 902  FORMAT(1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 903  FORMAT(1X,I2,1X,I2,3X,1P,E16.8,0P,3X,'#',1X,A)
 904  FORMAT('DECAY',1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 905  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A)
 906  FORMAT('#',1X,A,3X,E16.8)
 907  FORMAT(A,1P,E16.8,A)
 908  FORMAT(2E16.8)
 909  FORMAT(8X,A,12X,A)
 910  FORMAT(E16.8,3X,A)
 911  FORMAT(A,F6.3,A,F6.3,A)
 913  FORMAT(2X,I4,I4,10X,A)
 914  FORMAT(1X,I5,3X,1P,I16,0P,3X,'#',1X,A)
 915  FORMAT(A,A)
 916  FORMAT(A,I1,4E16.8)
 917  FORMAT(A,4E16.8)
 918  FORMAT(1X,I5,3X,A,F6.1,'-',F5.1,A)
 919  FORMAT(1X,I5,3X,A,F6.2)
 920  FORMAT('#',0P,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 921  FORMAT(1X,I2,1X,I2,3X,'#',1X,A)

      END


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

      SUBROUTINE EFFECTIVECOUPLINGS(PAR)

      DOUBLE PRECISION tanb,SB2,CB2,h1,h2,M1,M2,HTAU
      DOUBLE PRECISION L,K,AL,AK,MU,NU,RUNMB,LQT,ALSMT,HT,HB
      DOUBLE PRECISION ALSMA,DLA,DLQA,F1,F2,HTMA
      DOUBLE PRECISION MA2,PI,COEF,T,ATB,fmt,fmb,rt,rb
      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION sb,cb,s2,Lmu_mh,Lnu_mh
      DOUBLE PRECISION LM2,Lmu,Lnu,Lmunu,Lmunu_mh,LM2mu_mh
      DOUBLE PRECISION LQZ,LA,LP,LS,LPP,P1,P2,LMAMT,LM1mu,LM2mu
      DOUBLE PRECISION SUBDET,BOS,GAUGE,SFERM,MP2,MS2,M12,Lmax1,MGAU


      DOUBLE PRECISION PX,PA(6),PB(2),PL(7),PK(8)
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION Q2,QSTSB
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION MQ3P,MU3P,MD3P,ATP,ABP !ATP,ABP trilinears at QSTSB
     
       DOUBLE PRECISION HTQ,HBQ,MTQ,MBQ
       DOUBLE PRECISION mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
       DOUBLE PRECISION mstL,mstR,Wt,msbL,msbR,Wb


      COMMON/QEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/RADCOR2/MQ3P,MU3P,MD3P,ATP,ABP
      COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ
      !These are the DRbar stop and sbottom (masses)**2
      COMMON/RADCOR/mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
   
*     For loops
      LQT=DLOG(MAX(QSTSB,MT**2)/MT**2)
      T= DLOG(QSTSB/MT**2)
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)
            
*   Trig. functions of beta
      cb=1d0/DSQRT(1d0+tanbq**2)
      sb=tanbq*cb
      s2=2d0*sb*cb

*   VEVS
      VD=1d0/DSQRT(2d0*DSQRT(2d0)*GF*(1d0+PAR(3)**2))
      VU=VD*PAR(3)
      V2=VU**2+VD**2
      !     VS=MUQ/LQ*DSQRT(ZS)

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
    
      !MA2=MAX(((PAR(5)+NU)*MU+M3H+MU*MUP+L*XIF)*(tanb+1d0/tanb),1d0)
      MA2= PAR(23)**2
 
*   Could alternatively pull from slha2 file if we make defs match
*   Approximate value for the Singlet-like CP odd Higgs mass squared:   
      MP2=-3d0*NUQ*AKQ-XIF*(4d0*KQ+LQ*MUP/MUQ)
     .     -2d0*MSP-MUP*NUQ-LQ*XIS/MUQ
*   Approximate value for the Singlet-like CP even Higgs mass squared:
      MS2=MAX(NUQ*(AKQ+4d0*NUQ+3d0*MUP)-LQ*(XIS+XIF*MUP)/MUQ,MZ**2)
*   Approximate value for the off-diag. CP odd mass matrix element:
      M12=LQ*DSQRT(h1q**2+h2q**2)*(ALQ-2d0*NUQ-MUP)


      MGAU= (G1q*M1+3d0*G2q*M2)

      bos=0d0

      P2= MAX(MA2+MP2,MZ**2)
      P1= MAX((MA2*MP2-M12**2)/P2,MZ**2)
      LA= DLOG(MAX(MA2,MZ**2)/MZ**2)
      LS= DLOG(MS2/MZ**2)
      LP= DLOG(P2/MZ**2)
      LPP= DLOG(P2/P1)
      LQZ= DLOG(QSTSB/MZ**2)
      IF(2d0*NUQ+MUP.EQ.0d0)THEN
         Lmax1=0d0
      ELSE
         Lmax1= DLOG(MAX(MUQ**2/(2d0*NUQ+MUP)**2,1d0))
      ENDIF
    
*    loop corrections (can we replace these?)
  
      bos= COEF*MZ**2/GQ*((GQ**2*(2d0*S2TW**2
     .     -2d0*S2TW*(1d0+s2**2)-11d0/4d0*s2**4+5d0*s2**2+3d0/4d0)
     .     +GQ*LQ**2*(2d0*S2TW*s2**2+11d0/2d0*s2**4-15d0/2d0*s2**2
     .     -1d0)+LQ**4*(-11d0/4d0*s2**4+5d0/2d0*s2**2+1d0))*LA
     .     +(LQ**2*(LQ-KQ*s2)**2+3d0*LQ**2/MS2*(GQ+(LQ**2-GQ)*s2**2)
     .     *(2d0*MUQ-s2*(ALQ+2d0*NUQ+MUP))**2
     .     -LQ**4/MS2**2*(2d0*MUQ-s2*(ALQ+2d0*NUQ+MUP))**4)*LS
     .     +(GQ**2/4d0*(1d0-s2**4)+GQ*LQ**2*(1d0/2d0*s2**4
     .     +1d0/2d0*s2**2-1d0)+LQ**4*(-1d0/4d0*s2**4-1d0/2d0*s2**2
     .     +1d0)+LQ**2*(LQ+KQ*s2)**2)*LP
     .     -((GQ-LQ**2)**2/2d0*MP2/P2*s2**2*(1d0-s2**2)
     .     +(LQ*MA2*(LQ+KQ*s2)-LQ**2*(ALQ-2d0*NUQ-MUP)**2
     .     -MP2/2d0*(GQ*(1d0-s2**2)
     .     -LQ**2*(2d0-s2**2)))**2/P2**2-LQ*MA2*MP2*(LQ+KQ*s2)
     .     *(GQ*(1d0-s2**2)-LQ**2*(2d0-s2**2))/P2**2)*LPP
     .     +(GQ**2*(-4d0+S2**2+2d0*S2TW*(1d0+S2**2)-2d0*S2TW**2)
     .     +GQ*LQ**2*(2d0+S2**2-2d0*S2**2*S2TW)
     .     -LQ**4*(4d0+2d0*S2**2)-2d0*LQ**2*KQ**2*S2**2)*LQZ)

    
      GAUGE= MZ**2*GQ*COEF*(-9d0+12d0*S2TW-6d0*S2TW**2)*LQZ
      

      sferm= COEF/2d0*(h1q**2+h2q**2)*
     .   ((g1q**2/6d0+3d0/2d0*g2q**2**2)*DLOG(MAX(PAR(7),MZ**2)/QSTSB)
     .   +4d0/3d0*g1q**2*DLOG(MAX(PAR(8),MZ**2)/QSTSB)
     .   +1d0/3d0*g1q**2*DLOG(MAX(PAR(9),MZ**2)/QSTSB)
     .   +(g1q**2/3d0+3d0*g2q**2)*DLOG(MAX(PAR(15),MZ**2)/QSTSB)
     .   +8d0/3d0*g1q**2*DLOG(MAX(PAR(16),MZ**2)/QSTSB)
     .   +2d0/3d0*g1q**2*DLOG(MAX(PAR(17),MZ**2)/QSTSB)
     .   +(g1q**2/2d0+g2q**2/2d0)*DLOG(MAX(PAR(10),MZ**2)/QSTSB)
     .   +g1q**2*DLOG(MAX(PAR(11),MZ**2)/QSTSB)
     .   +(g1q**2+g2q**2)*DLOG(MAX(PAR(18),MZ**2)/QSTSB)
     .      +2d0*g1q**2*DLOG(MAX(PAR(19),MZ**2)/QSTSB))

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
    
*    PROBLEM:  mhiggs.f has a different definition of Lmunu etc from runpar.f
      Lmunu_mh= DLOG(MAX(MUQ**2,(2d0*NUQ+MUP)**2,MZ**2)/QSTSB)
      LM2mu_mh = DLOG(MAX(M2**2,MUQ**2,MZ**2)/QSTSB)
      Lmu_mh= DLOG(MAX(MUQ**2,MZ**2)/QSTSB)
      Lnu_mh= DLOG(MAX((2d0*NUQ+MUP)**2,MZ**2)/QSTSB)

*     MIXING
      XT=2d0*ATP**2/QSTSB*(1d0-ATP**2/QSTSB/12d0)
      XB=2d0*ABP**2/QSTSB*(1d0-ABP**2/QSTSB/12d0)
      ATB=-MUQ**2/QSTSB-(MUQ**2-ATP*ABP)**2/QSTSB**2/6d0
     .   +(ATP+ABP)**2/QSTSB/2d0
      ALSMA=ALSMT/(1d0+ALSMT/(4d0*PI)*(7d0*LMAMT-2d0*
     .    DLOG(MAX(MA2,PAR(22)**2)/MAX(PAR(22)**2,MT**2))))
      DLA=(ALSMA/ALSMT)**(1d0/7d0)
      DLQA=(ALSQ/ALSMA)**(1d0/7d0)
      F1=DABS(1d0-9d0*SB2*HT**2*(1d0-DLA)/(8d0*PI*ALSMT))
      HTMA=HT*DLA**4/DSQRT(F1)
      F2=DABS(1d0-9d0*HTMA**2*(1d0-DLQA)/(8d0*PI*ALSMA))

*    We need DRbar running stop and sbottom masses      

      mstL= MQ3P + mtq**2 + (gQ/2d0-g1Q/3d0)*(h2q**2-h1q**2)
      mstR= MU3P + mtq**2 + g1Q/3d0*(h2q**2-h1q**2)
      Wt= DSQRT( (mstL-mstR)**2 + 4d0*((ATP-muq/tanbq)*mtq)**2)
      MST1= 0.5d0*(mstL+mstR-Wt)
      MST2= 0.5d0*(mstL+mstR+Wt)
*   Running sbottom masses squared and mixings

      msbL= MQ3P + mbq**2 - (gQ/2d0-g1Q/6d0)*(h2q**2-h1q**2)
      msbR= MD3P + mbq**2 - g1Q/6d0*(h2q**2-h1q**2)
      Wb= DSQRT((msbL-msbR)**2 + 4d0*((ABP-muq*tanbq)*mbq)**2)
      MSB1= 0.5d0*(msbL+msbR-Wb)
      MSB2= 0.5d0*(msbL+msbR+Wb)


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
      
 
      
      IF(mst1.NE.mst2)THEN
         fmt= (mst2*DLOG(mst2/QSTSB)-mst1*DLOG(mst1/QSTSB))/
     .        (mst2-mst1)-1d0
      ELSE
         fmt= DLOG(mst1/QSTSB)
      ENDIF

     
      IF(msb1.NE.msb2)THEN
       fmb= (msb2*DLOG(msb2/QSTSB)-msb1*DLOG(msb1/QSTSB))/
     .    (msb2-msb1)-1d0
      ELSE
         fmb= DLOG(msb1/QSTSB)
      ENDIF
      
    
      rt= 3d0/2d0*COEF*htq**2
      rb= 3d0/2d0*COEF*hbq**2
   
*   Now we fill the effective couplings
  
      PX=(MUP*XIF+XIS
     .     -COEF*2d0*KQ*MUP**3*LNU
     .     )/DSQRT(ZS)

      PA(1)=(0d0
     .  -COEF*KQ*(4d0*LQ**2*MUP*LMUNU_mh+8d0*KQ**2*Lmax1
     .  ))/DSQRT(ZHD**2*ZS)

      PA(2)=(0d0
     .  -COEF*KQ*(4d0*LQ**2*MUP*LMUNU_mh+8d0*KQ**2*Lmax1
     .  ))/DSQRT(ZHU**2*ZS)

      PA(3)=(KQ*MUP
     .  -COEF*8d0*KQ**3*MUP*LNU
     .  )/DSQRT(ZS**3)

      PA(4)=(KQ/3d0*AKQ
     .  )/DSQRT(ZS**3)

      PA(5)=(LQ*(ALQ+2d0*rt*Atp*fmt+2d0*rb*Abp*fmb)
     .  -COEF*8d0*LQ*KQ**2*MUP*Lmax1
     .  )/DSQRT(ZS*ZHU*ZHD)

      PA(6)=(LQ*MUP
     .  +COEF*(LQ*MGAU*LM2MU_mh-2d0*LQ**3*MUP*LMUNU_mh
     .  -8d0*LQ*KQ**2*MUP*Lmax1))/DSQRT(ZS*ZHU*ZHD)

    

      PB(1)=(M3H+LQ*XIF
     .  -COEF*4d0*LQ*KQ*MUP**2*Lmax1
     .  )/DSQRT(ZHU*ZHD)

      PB(2)=(MSP/2d0+KQ*XIF
     .  -COEF*2d0*KQ**2*MUP**2*Lnu
     .  )/ZS

      PL(1)=2d0*(GQ/4d0
     .  +COEF*(3d0*hbQ**4*(T+XB/2d0)-htQ**4/4d0*MUQ**4/QSTSB**2)
     .  +3d0/4d0*COEF*GQ*(htq**2*MUQ**2-hbq**2*ABP**2)/QSTSB
     .  +3d0/2d0*COEF**2*hbq**4*
     .   (T**2*(64d0*PI*ALSQ-2d0/3d0*g1q-3d0*hbq**2*cb**2
     .  +3d0*htq**2*sb**2)+((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
     .  -(dlog(MAX(MA2,MT**2)/MT**2))**2)*
     .  (-3d0*hbq**2*sb**2-htq**2*(3d0*sb**2+1d0)))
     .  +(bos+gauge+sferm)/(4d0*V2)
     .  -COEF*((GQ**2+G2Q**2)*LM2MU_mh+LQ**4*LMUNU_mh-8d0*KQ**4*Lmax1)
     .   )/ZHD**2



      PL(2)=2d0*(GQ/4d0
     .  +COEF*(3d0*htQ**4*(T+XT/2d0)-hbQ**4/4d0*MUQ**4/QSTSB**2)
     .  +3d0/4d0*COEF*GQ*(hbq**2*MUQ**2-htq**2*ATP**2)/QSTSB
     .  +3d0/2d0*COEF**2*htq**4*
     .  (T**2*(64d0*PI*ALSQ+4d0/3d0*g1q-3d0*htq**2*sb**2
     .  +3d0*hbq**2*cb**2)+((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
     .  -(dlog(MAX(MA2,MT**2)/MT**2))**2)*
     .  (-3d0*htq**2*cb**2-hbq**2*(3d0*cb**2+1d0)))
     .  +(bos+gauge+sferm)/(4d0*V2)
     .  -COEF*((GQ**2+G2Q**2)*LM2MU_mh+LQ**4*LMUNU_mh-8d0*KQ**4*Lmax1)
     .   )/ZHU**2

      PL(3)=((G2Q-G1Q)/4d0
     .  +6d0*COEF*htQ**2*hbQ**2*(T+ATB/2d0)
     .  +COEF/2d0**htQ**4*MUQ**2/QSTSB*(3d0-ATP**2/QSTSB)
     .  +COEF/2d0**hbQ**4*MUQ**2/QSTSB*(3d0-ABP**2/QSTSB)
     .  +3d0/8d0*COEF*(G1Q-G2Q)*(htq**2*(ATP**2-MUQ**2)
     .  +hbq**2*(ABP**2-MUQ**2))/QSTSB
     .  +6d0*COEF**2*htQ**2*hbQ**2*(32d0*PI*ALSQ-htq**2-hbq**2)
     .  *((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
     .  -(dlog(MAX(MA2,MT**2)/MT**2))**2)
     .  +(bos+gauge-sferm)/(2d0*V2)
     .  +COEF*G2Q**2*(3d0/4d0*DLOG(MAX(PAR(7),MZ**2)/QSTSB)
     .  +3d0/2d0*DLOG(MAX(PAR(15),MZ**2)/QSTSB)
     .  +1d0/4d0*DLOG(MAX(PAR(10),MZ**2)/QSTSB)
     .  +1d0/2d0*DLOG(MAX(PAR(18),MZ**2)/QSTSB))
     .  +COEF*(-(2d0*GQ**2-4d0*GQ*G2Q+4d0*G2Q**2)*LM2mu_mh
     .  -2d0*LQ**4*Lmunu_mh+16d0*KQ**4*Lmax1)
     .   )/(ZHU*ZHD)

      PL(4)=(LQ**2-G2Q/2d0
     .  -6d0*COEF*htQ**2*hbQ**2*(T+ATB/2d0)
     .  +COEF/2d0**htQ**4*MUQ**2/QSTSB*(3d0-ATP**2/QSTSB)
     .  +COEF/2d0**hbQ**4*MUQ**2/QSTSB*(3d0-ABP**2/QSTSB)
     .  +3d0/4d0*COEF*G2Q*(htq**2*(ATP**2-MUQ**2)
     .  +hbq**2*(ABP**2-MUQ**2))/QSTSB
     .  -6d0*COEF**2*htQ**2*hbQ**2*(32d0*PI*ALSQ-htq**2-hbq**2)
     .  *((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
     .  -(dlog(MAX(MA2,MT**2)/MT**2))**2)
     .  -COEF*G2Q**2*(3d0/4d0*DLOG(MAX(PAR(7),MZ**2)/QSTSB)
     .  +3d0/2d0*DLOG(MAX(PAR(15),MZ**2)/QSTSB)
     .  +1d0/4d0*DLOG(MAX(PAR(10),MZ**2)/QSTSB)
     .  +1d0/2d0*DLOG(MAX(PAR(18),MZ**2)/QSTSB))
     .  +COEF*(G2Q*(G2Q-7d0*G1Q)/4d0*LQZ+2d0*G2Q*(G2Q-G1Q)*LM2mu_mh
     .  +16d0*LQ**2*KQ**2*Lmax1)
     .   )/(ZHU*ZHD)

      PL(5)=2d0*(0d0
     .  -COEF/4d0*htQ**4*MUQ**2*ATP**2/QSTSB**2
     .  -COEF/4d0*hbQ**4*MUQ**2*ABP**2/QSTSB**2
     .  +COEF*8d0*LQ**2*KQ**2*Lmax1
     .  )/(ZHU*ZHD)

      PL(6)=(0d0
     .  -COEF/2d0*htQ**4*MUQ**3*ATP/QSTSB**2
     .  +3d0*COEF*hbQ**4*ABP*MUQ/QSTSB*(1d0-ABP**2/6d0/QSTSB)
     .  +3d0/4d0*COEF*GQ*MUQ*(htq**2*ATP-hbq**2*ABP)/QSTSB
     .  +COEF*16d0*LQ*KQ**3*Lmax1
     .   )/DSQRT(ZHU*ZHD**3)

      PL(7)=(0d0
     .  -COEF/2d0*hbQ**4*MUQ**3*ABP/QSTSB**2
     .  +3d0*COEF*htQ**4*ATP*MUQ/QSTSB*(1d0-ATP**2/6d0/QSTSB)
     .  +3d0/4d0*COEF*GQ*MUQ*(hbq**2*ABP-htq**2*ATP)/QSTSB
     .  +COEF*16d0*LQ*KQ**3*Lmax1
     .   )/DSQRT(ZHU**3*ZHD)

      PK(1)=(LQ**2
     .  +COEF*(-LQ**2*(G1Q+3d0*G2Q)*LM2MU_mh
     .  -(2d0*LQ**4+8d0*LQ**2*KQ**2)*LMUNU_mh
     .  -16d0*KQ**4*Lmax1
     .  ))/(ZS*ZHD)

      PK(2)=(LQ**2
     .  +COEF*(-LQ**2*(G1Q+3d0*G2Q)*LM2MU_mh
     .  -(2d0*LQ**4+8d0*LQ**2*KQ**2)*LMUNU_mh
     .  -16d0*KQ**4*Lmax1
     .  ))/(ZS*ZHU)

      PK(3)=(KQ**2
     .  +COEF*(-2d0*LQ**4*LMU_mh-8d0*KQ**4*LNU_mh)
     .  )/ZS**2

      PK(4)=(0d0
     .  -COEF*16d0*LQ*KQ**3*Lmax1
     .  )/DSQRT(ZS**2*ZHU*ZHD)

      PK(5)=(0d0
     .  )/DSQRT(ZS**2*ZHU*ZHD)

      PK(6)=(LQ*KQ
     .  -COEF*4d0*LQ**3*KQ*LMUNU_mh
     .  )/DSQRT(ZS**2*ZHU*ZHD)

      PK(7)=(0d0
     .  )/ZS**2

      PK(8)=(0d0
     .     )/ZS**2

      
      END
      
