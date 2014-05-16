      SUBROUTINE RGESUNIGM(PAR,IFAIL,MESTEST)

*   Subroutine to integrate the 2-loop RGEs for all soft terms
*   from the SUSY scale Q2 up to MMESS, through a CALL of the
*   subroutine ODEINTS that is part of the file integs.f
*
*   The gauge/Yukawa couplings at Q2 are read in from
*   COMMON/SUSYCOUP, initialized in RGESGM.
*
*   The soft terms at Q2 are read in from PAR(*),
*   COMMON/SUSYMH and COMMON/SUSYEXT
*
*   PURPOSE: Once the gauge/Yukawa couplings at Q2 have been computed
*   in RGESGM (including SUSY thresholds), the agreement of the soft
*   terms at MMESS with the inputs is tested via MESTEST.
*
*   If MAFLAG=-1,-3 MS**2 at MMESS is computed here
*   If MAFLAG=-2,-4 XIS at MMESS is computed here
*   If MAFLAG=-3,-4 XIF at MMESS is computed here
*   They are stored in COMMON/MESEXT
*
*   NOTE: Y(5)=kappa=PAR(2), NOT kappa^2, since the sign of kappa
*   matters in the RGEs (lambda is assumed to be positive)
*
***********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,OMGFLAG,MAFLAG,NN
      PARAMETER (NN=35)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN),PI,COEF
      DOUBLE PRECISION g1s,g2s,g3s,HTOPS,HBOTS,HTAUS,Q2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION M1INP,M2INP,M3INP,AKINP,ATINP,ABINP,
     .      ATAUINP,AMUINP,MH1INP,MH2INP,MQ3INP,MU3INP,MD3INP,
     .      MQINP,MUINP,MDINP,ML3INP,ME3INP,MLINP,MEINP
      DOUBLE PRECISION M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     .      ATAUMES,AMUMES,MH1MES,MH2MES,MQ3MES,MU3MES,MD3MES,
     .      MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      DOUBLE PRECISION MH1S,MH2S,MSS,MESTEST
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     .      HBOTMES,HTAUMES
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP,DELHINP
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION MSUSYEFF,MMESS,N5

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/SOFTINP/M1INP,M2INP,M3INP,AKINP,ATINP,ABINP,
     .      ATAUINP,AMUINP,MH1INP,MH2INP,MQ3INP,MU3INP,MD3INP,
     .      MQINP,MUINP,MDINP,ML3INP,ME3INP,MLINP,MEINP
      COMMON/SOFTMES/M1MES,M2MES,M3MES,ALMES,AKMES,ATMES,ABMES,
     .      ATAUMES,AMUMES,MH1MES,MH2MES,MQ3MES,MU3MES,MD3MES,
     .      MQMES,MUMES,MDMES,ML3MES,ME3MES,MLMES,MEMES
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     .      HBOTMES,HTAUMES      
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP,DELHINP
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES

      EXTERNAL DERIVSS,RKQSS

      EPS=1.D-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

* Definition of the couplings squared Y(I) at Q2=M_SUSY

      Y(1)=g1s
      Y(2)=g2s
      Y(3)=g3s
      Y(4)=PAR(1)**2
* NOTE: Y(5)=K, NOT K**2
      Y(5)=PAR(2)
      Y(6)=HTOPS**2
      Y(7)=HBOTS**2
      Y(8)=HTAUS**2

* Definition of the soft terms Y(I) at Q2=M_SUSY

      Y(9)=PAR(20)
      Y(10)=PAR(21)
      Y(11)=PAR(22)
      Y(12)=PAR(5)
      Y(13)=PAR(6)
      Y(14)=PAR(12)
      Y(15)=PAR(13)
      Y(16)=PAR(14)
      Y(17)=MH1S
      Y(18)=MH2S
      Y(19)=MSS
      Y(20)=PAR(7)
      Y(21)=PAR(8)
      Y(22)=PAR(9)
      Y(23)=PAR(15)
      Y(24)=PAR(16)
      Y(25)=PAR(17)
      Y(26)=PAR(10)
      Y(27)=PAR(11)
      Y(28)=PAR(18)
      Y(29)=PAR(19)
      Y(30)=XIFSUSY
      Y(31)=XISSUSY
      Y(32)=MUPSUSY
      Y(33)=MSPSUSY
      Y(34)=M3HSUSY
      Y(35)=PAR(25)
      
      X1=0d0
      X2=COEF*DLOG(MMESS**2/Q2)

      !WRITE(0,*)"CALL RGESUNI"
      !WRITE(0,*)""
      !WRITE(0,*)"MSUSY =",DSQRT(Q2)
      !WRITE(0,*)"G1 =",5d0/3d0*Y(1)
      !WRITE(0,*)"G2 =",Y(2)
      !WRITE(0,*)"G3 =",Y(3)
      !WRITE(0,*)"L =",DSQRT(Y(4))
      !WRITE(0,*)"K =",Y(5)
      !WRITE(0,*)"HT =",Y(6)
      !WRITE(0,*)"HB =",Y(7)
      !WRITE(0,*)"HL =",Y(8)
      !WRITE(0,*)"M1 =",Y(9)
      !WRITE(0,*)"M2 =",Y(10)
      !WRITE(0,*)"M3 =",Y(11)
      !WRITE(0,*)"AL =",Y(12)
      !WRITE(0,*)"AK =",Y(13)
      !WRITE(0,*)"ATOP =",Y(14)
      !WRITE(0,*)"ABOT =",Y(15)
      !WRITE(0,*)"ATAU =",Y(16)
      !WRITE(0,*)"AMUON =",Y(35)
      !WRITE(0,*)"MH1 =",Y(17)
      !WRITE(0,*)"MH2 =",Y(18)
      !WRITE(0,*)"MS =",Y(19)
      !WRITE(0,*)"MQ3 =",Y(20)
      !WRITE(0,*)"MU3 =",Y(21)
      !WRITE(0,*)"MD3 =",Y(22)
      !WRITE(0,*)"MQ =",Y(23)
      !WRITE(0,*)"MU =",Y(24)
      !WRITE(0,*)"MD =",Y(25)
      !WRITE(0,*)"ML3 =",Y(26)
      !WRITE(0,*)"ME3 =",Y(27)
      !WRITE(0,*)"ML =",Y(28)
      !WRITE(0,*)"ME =",Y(29)
      !WRITE(0,*)"XIF =",Y(30)
      !WRITE(0,*)"XIS =",Y(31)
      !WRITE(0,*)"MUP =",Y(32)
      !WRITE(0,*)"MSP =",Y(33)
      !WRITE(0,*)"M3H =",Y(34)
      !WRITE(0,*)""

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)

      !WRITE(0,*)"MMESS =",MMESS
      !WRITE(0,*)"G1MES =",Y(1)
      !WRITE(0,*)"G2MES =",Y(2)
      !WRITE(0,*)"G3MES =",Y(3)
      !WRITE(0,*)"LMES =",DSQRT(Y(4))
      !WRITE(0,*)"KMES =",Y(5)
      !WRITE(0,*)"HTMES =",Y(6)
      !WRITE(0,*)"HBMES =",Y(7)
      !WRITE(0,*)"HLMES =",Y(8)
      !WRITE(0,*)"M1MES =",Y(9)," M1INP =",M1INP," ERR =",
!     .      (Y(9)-M1INP)**2/(1.D2+M1INP**2)
      !WRITE(0,*)"M2MES =",Y(10)," M2INP =",M2INP," ERR =",
!     .      (Y(10)-M2INP)**2/(1.D2+M2INP**2)
      !WRITE(0,*)"M3MES =",Y(11)," M3INP =",M3INP," ERR =",
!     .      (Y(11)-M3INP)**2/(1.D2+M3INP**2)
      !WRITE(0,*)"ALMES =",Y(12)," ALINP =",ALINP," ERR =",
!     .      (Y(12)-ALINP)**2/(1.D2+ALINP**2)
      !WRITE(0,*)"AKMES =",Y(13)," AKINP =",AKINP," ERR =",
!     .      (Y(13)-AKINP)**2/(1.D2+AKINP**2)
      !WRITE(0,*)"ATOPMES =",Y(14)," ATINP =",ATINP," ERR =",
!     .      (Y(14)-ATINP)**2/(1.D2+ATINP**2)
      !WRITE(0,*)"ABOTMES =",Y(15)," ABINP =",ABINP," ERR =",
!     .      (Y(15)-ABINP)**2/(1.D2+ABINP**2)
      !WRITE(0,*)"ATAUMES =",Y(16)," ATAUINP =",ATAUINP," ERR =",
!     .      (Y(16)-ATAUINP)**2/(1.D2+ATAUINP**2)
      !WRITE(0,*)"AMUMES =",Y(35)," AMUINP =",AMUINP," ERR =",
!     .      (Y(35)-AMUINP)**2/(1.D2+AMUINP**2)
      !WRITE(0,*)"MH1MES =",Y(17)," MH1INP=",MH1INP," ERR =",
!     .      (Y(17)-MH1INP)**2/(1.D4+MH1INP**2)
      !WRITE(0,*)"MH2MES =",Y(18)," MH2INP=",MH2INP," ERR =",
!     .      (Y(18)-MH2INP)**2/(1.D4+MH2INP**2)
!      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
       !WRITE(0,*)"MSMES =",Y(19)," MSINP=",MSINP," ERR =",
!     .      (Y(19)-MSINP)**2/(1.D4+MSINP**2)
!      ELSE
       !WRITE(0,*)"MSMES =",Y(19)
!      ENDIF
      !WRITE(0,*)"MQ3MES =",Y(20)," MQ3INP =",MQ3INP," ERR =",
!     .      (Y(20)-MQ3INP)**2/(1.D4+MQ3INP**2)
      !WRITE(0,*)"MU3MES =",Y(21)," MU3INP =",MU3INP," ERR =",
!     .      (Y(21)-MU3INP)**2/(1.D4+MU3INP**2)
      !WRITE(0,*)"MD3MES =",Y(22)," MD3INP =",MD3INP," ERR =",
!     .      (Y(22)-MD3INP)**2/(1.D4+MD3INP**2)
      !WRITE(0,*)"MQMES =",Y(23)," MQINP =",MQINP," ERR =",
!     .      (Y(23)-MQINP)**2/(1.D4+MQINP**2)
      !WRITE(0,*)"MUMES =",Y(24)," MUINP =",MUINP," ERR =",
!     .      (Y(24)-MUINP)**2/(1.D4+MUINP**2)
      !WRITE(0,*)"MDMES =",Y(25)," MDINP =",MDINP," ERR =",
!     .      (Y(25)-MDINP)**2/(1.D4+MDINP**2)
      !WRITE(0,*)"ML3MES =",Y(26)," ML3INP =",ML3INP," ERR =",
!     .      (Y(26)-ML3INP)**2/(1.D4+ML3INP)
      !WRITE(0,*)"ME3MES =",Y(27)," ME3INP =",ME3INP," ERR =",
!     .      (Y(27)-ME3INP)**2/(1.D4+ME3INP**2)
      !WRITE(0,*)"MLMES =",Y(28)," MLINP =",MLINP," ERR =",
!     .      (Y(28)-MLINP)**2/(1.D4+MLINP**2)
      !WRITE(0,*)"MEMES =",Y(29)," MEINP =",MEINP," ERR =",
!     .      (Y(29)-MEINP)**2/(1.D4+MEINP**2)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
       !WRITE(0,*)"XIFMES =",Y(30)," XIFINP=",XIFINP," ERR =",
!     .      (Y(30)-XIFINP)**2/(1.D4+XIFINP**2)
!      ELSE
       !WRITE(0,*)"XIFMES =",Y(30)
!      ENDIF
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       !WRITE(0,*)"XISMES =",Y(31)," XISINP=",XISINP," ERR =",
!     .      (Y(31)-XISINP)**2/(1.D6+XISINP**2)
!      ELSE
       !WRITE(0,*)"XISMES =",Y(31)
!      ENDIF
      !WRITE(0,*)"MUPMES =",Y(32)," MUPINP  =",MUPINP," ERR =",
!     .      (Y(32)-MUPINP)**2/(1.D2+MUPINP**2)
      !WRITE(0,*)"MSPMES =",Y(33)," MSPINP =",MSPINP," ERR =",
!     .      (Y(33)-MSPINP)**2/(1.D4+MSPINP**2)
      !WRITE(0,*)""

      IF(IFAIL.NE.0)THEN
       IFAIL=12
       RETURN
      ENDIF

* MESTEST:

      MESTEST=0d0
      MESTEST=MESTEST+(Y(9)-M1INP)**2/(1.D2+M1INP**2)
      MESTEST=MESTEST+(Y(10)-M2INP)**2/(1.D2+M2INP**2)
      MESTEST=MESTEST+(Y(11)-M3INP)**2/(1.D2+M3INP**2)
      MESTEST=MESTEST+(Y(12)-ALINP)**2/(1.D2+ALINP**2)
      MESTEST=MESTEST+(Y(13)-AKINP)**2/(1.D2+AKINP**2)
      MESTEST=MESTEST+(Y(14)-ATINP)**2/(1.D2+ATINP**2)
      MESTEST=MESTEST+(Y(15)-ABINP)**2/(1.D2+ABINP**2)
      MESTEST=MESTEST+(Y(16)-ATAUINP)**2/(1.D2+ATAUINP**2)
      MESTEST=MESTEST+(Y(17)-MH1INP)**2/(1.D4+MH1INP**2)
      MESTEST=MESTEST+(Y(18)-MH2INP)**2/(1.D4+MH2INP**2)
      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
       MESTEST=MESTEST+(Y(19)-MSINP)**2/(1.D4+MSINP**2)
      ENDIF
      MESTEST=MESTEST+(Y(20)-MQ3INP)**2/(1.D4+MQ3INP**2)
      MESTEST=MESTEST+(Y(21)-MU3INP)**2/(1.D4+MU3INP**2)
      MESTEST=MESTEST+(Y(22)-MD3INP)**2/(1.D4+MD3INP**2)
      MESTEST=MESTEST+(Y(23)-MQINP)**2/(1.D4+MQINP**2)
      MESTEST=MESTEST+(Y(24)-MUINP)**2/(1.D4+MUINP**2)
      MESTEST=MESTEST+(Y(25)-MDINP)**2/(1.D4+MDINP**2)
      MESTEST=MESTEST+(Y(26)-ML3INP)**2/(1.D4+ML3INP**2)
      MESTEST=MESTEST+(Y(27)-ME3INP)**2/(1.D4+ME3INP**2)
      MESTEST=MESTEST+(Y(28)-MLINP)**2/(1.D4+MLINP**2)
      MESTEST=MESTEST+(Y(29)-MEINP)**2/(1.D4+MEINP**2)
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
       MESTEST=MESTEST+(Y(30)-XIFINP)**2/(1.D4+XIFINP**2)
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       MESTEST=MESTEST+(Y(31)-XISINP)**2/(1.D6+XISINP**2)
      ENDIF
      MESTEST=MESTEST+(Y(32)-MUPINP)**2/(1.D2+MUPINP**2)
      MESTEST=MESTEST+(Y(33)-MSPINP)**2/(1.D4+MSPINP**2)
      MESTEST=MESTEST+(Y(34))**2/1.D4
      MESTEST=MESTEST+(Y(35)-AMUINP)**2/(1.D2+AMUINP**2)
      
* Couplings at the messenger scale

      G1MES=Y(1)
      G2MES=Y(2)
      G3MES=Y(3)
      LMES=Y(4)
      KMES=Y(5)
      HTOPMES=Y(6)
      HBOTMES=Y(7)
      HTAUMES=Y(8)
      
* Soft terms at the messenger scale

      M1MES=Y(9)
      M2MES=Y(10)
      M3MES=Y(11)
      ALMES=Y(12)
      AKMES=Y(13)
      ATMES=Y(14)
      ABMES=Y(15)
      ATAUMES=Y(16)
      MH1MES=Y(17)
      MH2MES=Y(18)
      MSMES=Y(19)
      MQ3MES=Y(20)
      MU3MES=Y(21)
      MD3MES=Y(22)
      MQMES=Y(23)
      MUMES=Y(24)
      MDMES=Y(25)
      ML3MES=Y(26)
      ME3MES=Y(27)
      MLMES=Y(28)
      MEMES=Y(29)
      XIFMES=Y(30)
      XISMES=Y(31)
      MUPMES=Y(32)
      MSPMES=Y(33)
      M3HMES=Y(34)
      AMUMES=Y(35)
      
      END
