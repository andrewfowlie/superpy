      SUBROUTINE RGESUNI(PAR,IFAIL,GUTEST)

*   Subroutine to integrate the 2-loop RGEs for all soft terms
*   from the SUSY scale Q2 up to MGUT, through a CALL of the
*   subroutine ODEINTS that is part of the file integ.f
*
*   MGUT and the gauge/Yukawa couplings at Q2 are read in from
*   COMMON/MGUT and COMMON/SUSYCOUP, initialized in RGES.
*
*   The soft terms at Q2 are read in from PAR(*),
*   COMMON/SUSYMH and COMMON/SUSYEXT
*
*   PURPOSE: Once the gauge/Yukawa couplings at Q2 have been computed
*   in RGES (including SUSY thresholds), the agreement of the soft
*   terms at MGUT with the inputs is tested via GUTEST.
*
*   If MAFLAG=-5 MHU**2, MHD**2 at MGUT are computed here
*   If MAFLAG=-1,-3,-5 MS**2 at MGUT is computed here
*   If MAFLAG=-2,-4 XIS at MGUT is computed here
*   If MAFLAG=-3,-4 XIF at MGUT is computed here
*   They are stored in COMMON/GUTEXT and COMMON/GUTPAR
*
*   NOTE: Y(5)=kappa=PAR(2), NOT kappa^2, since the sign of kappa
*   matters in the RGEs (lambda is assumed to be positive)
*
***********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,NN,I,IM,OMGFLAG,MAFLAG
      PARAMETER (NN=35)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN),PI,COEF
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS,Q2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW,GUTEST
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT
      DOUBLE PRECISION HTOPGUT,HBOTGUT,HTAUGUT
      DOUBLE PRECISION M1INP,M2INP,M3INP,MH2INP,MH1INP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT
      DOUBLE PRECISION ATAUGUT,AMUGUT,MH1GUT,MH2GUT,MSGUT,MQ3GUT
      DOUBLE PRECISION MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT
      DOUBLE PRECISION ME3GUT,MLGUT,MEGUT,MH1S,MH2S,MSS
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION M0,M12,A0,MUM,MDM,MSM,MUT,MDT,MST

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     .      HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     .      ATAUGUT,AMUGUT,MH1GUT,MH2GUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     .      MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/INPPAR/M1INP,M2INP,M3INP,MH2INP,MH1INP,ALINP,AKINP,
     .      XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      COMMON/MSAVE/MUM,MDM,MSM,MUT,MDT,MST,IM

      EXTERNAL DERIVSS,RKQSS

      EPS=1.D-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)
      IM=IM+1
      IF(IM.EQ.21)THEN
       MUT=0d0
       MDT=0d0
       MST=0d0
       IM=1
      ENDIF
      
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
      X2=COEF*DLOG(MGUT**2/Q2)

!      WRITE(0,*)"CALL RGESUNI"
!      WRITE(0,*)""
!      WRITE(0,*)"MSUSY =",DSQRT(Q2)
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L =",DSQRT(Y(4))
!      WRITE(0,*)"K =",Y(5)
!      WRITE(0,*)"HT =",Y(6)
!      WRITE(0,*)"HB =",Y(7)
!      WRITE(0,*)"HL =",Y(8)
!      WRITE(0,*)"M1 =",Y(9)
!      WRITE(0,*)"M2 =",Y(10)
!      WRITE(0,*)"M3 =",Y(11)
!      WRITE(0,*)"AL =",Y(12)
!      WRITE(0,*)"AK =",Y(13)
!      WRITE(0,*)"ATOP =",Y(14)
!      WRITE(0,*)"ABOT =",Y(15)
!      WRITE(0,*)"ATAU =",Y(16)
!      WRITE(0,*)"AMUON =",Y(35)
!      WRITE(0,*)"MH1 =",Y(17)
!      WRITE(0,*)"MH2 =",Y(18)
!      WRITE(0,*)"MS =",Y(19)
!      WRITE(0,*)"MQ3 =",Y(20)
!      WRITE(0,*)"MU3 =",Y(21)
!      WRITE(0,*)"MD3 =",Y(22)
!      WRITE(0,*)"MQ =",Y(23)
!      WRITE(0,*)"MU =",Y(24)
!      WRITE(0,*)"MD =",Y(25)
!      WRITE(0,*)"ML3 =",Y(26)
!      WRITE(0,*)"ME3 =",Y(27)
!      WRITE(0,*)"ML =",Y(28)
!      WRITE(0,*)"ME =",Y(29)
!      WRITE(0,*)"XIF =",Y(30)
!      WRITE(0,*)"XIS =",Y(31)
!      WRITE(0,*)"MUP =",Y(32)
!      WRITE(0,*)"MSP =",Y(33)
!      WRITE(0,*)"M3H =",Y(34)
!      WRITE(0,*)""

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)

!      WRITE(0,*)"MGUT =",MGUT
!      WRITE(0,*)"G1GUT =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2GUT =",Y(2)
!      WRITE(0,*)"G3GUT =",Y(3)
!      WRITE(0,*)"LGUT =",DSQRT(Y(4))
!      WRITE(0,*)"KGUT =",Y(5)
!      WRITE(0,*)"HTGUT =",Y(6)
!      WRITE(0,*)"HBGUT =",Y(7)
!      WRITE(0,*)"HLGUT =",Y(8)
!      WRITE(0,*)"M1GUT =",Y(9)," M1INP =",M1INP," ERR =",
!     .      (Y(9)-M1INP)**2/(1.D2+M1INP**2)
!      WRITE(0,*)"M2GUT =",Y(10)," M2INP =",M2INP," ERR =",
!     .      (Y(10)-M2INP)**2/(1.D2+M2INP**2)
!      WRITE(0,*)"M3GUT =",Y(11)," M3INP =",M3INP," ERR =",
!     .      (Y(11)-M3INP)**2/(1.D2+M3INP**2)
!      WRITE(0,*)"ALGUT =",Y(12)," ALINP =",ALINP," ERR =",
!     .      (Y(12)-ALINP)**2/(1.D2+ALINP**2)
!      WRITE(0,*)"AKGUT =",Y(13)," AKINP =",AKINP," ERR =",
!     .      (Y(13)-AKINP)**2/(1.D2+AKINP**2)
!      WRITE(0,*)"ATOPGUT =",Y(14)," A0 =",A0," ERR =",
!     .      (Y(14)-A0)**2/(1.D2+A0**2)
!      WRITE(0,*)"ABOTGUT =",Y(15)," A0 =",A0," ERR =",
!     .      (Y(15)-A0)**2/(1.D2+A0**2)
!      WRITE(0,*)"ATAUGUT =",Y(16)," A0 =",A0," ERR =",
!     .      (Y(16)-A0)**2/(1.D2+A0**2)
!      WRITE(0,*)"AMUGUT =",Y(35)," A0 =",A0," ERR =",
!     .      (Y(35)-A0)**2/(1.D2+A0**2)
!      IF(MAFLAG.NE.-5)THEN
!       WRITE(0,*)"MH1GUT =",Y(17)," MH1INP=",MH1INP," ERR =",
!     .      (Y(17)-MH1INP)**2/(1.D4+MH1INP**2)
!      ELSE
!       WRITE(0,*)"MH1GUT =",Y(17)
!      ENDIF
!      IF(MAFLAG.NE.-5)THEN
!       WRITE(0,*)"MH2GUT =",Y(18)," MH2INP=",MH2INP," ERR =",
!     .      (Y(18)-MH2INP)**2/(1.D4+MH2INP**2)
!      ELSE
!       WRITE(0,*)"MH2GUT =",Y(18)
!      ENDIF
!      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
!       WRITE(0,*)"MSGUT =",Y(19)," MSINP=",MSINP," ERR =",
!     .      (Y(19)-MSINP)**2/(1.D4+MSINP**2)
!      ELSE
!       WRITE(0,*)"MSGUT =",Y(19)
!      ENDIF
!      WRITE(0,*)"MQ3GUT =",Y(20)," M0 =",M0**2," ERR =",
!     .      (Y(20)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"MU3GUT =",Y(21)," M0 =",M0**2," ERR =",
!     .      (Y(21)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"MD3GUT =",Y(22)," M0 =",M0**2," ERR =",
!     .      (Y(22)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"MQGUT =",Y(23)," M0 =",M0**2," ERR =",
!     .      (Y(23)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"MUGUT =",Y(24)," M0 =",M0**2," ERR =",
!     .      (Y(24)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"MDGUT =",Y(25)," M0 =",M0**2," ERR =",
!     .      (Y(25)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"ML3GUT =",Y(26)," M0 =",M0**2," ERR =",
!     .      (Y(26)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"ME3GUT =",Y(27)," M0 =",M0**2," ERR =",
!     .      (Y(27)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"MLGUT =",Y(28)," M0 =",M0**2," ERR =",
!     .      (Y(28)-M0**2)**2/(1.D4+M0**4)
!      WRITE(0,*)"MEGUT =",Y(29)," M0 =",M0**2," ERR =",
!     .      (Y(29)-M0**2)**2/(1.D4+M0**4)
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
!       WRITE(0,*)"XIFGUT =",Y(30)," XIFINP=",XIFINP," ERR =",
!     .      (Y(30)-XIFINP)**2/(1.D4+XIFINP**2)
!      ELSE
!       WRITE(0,*)"XIFGUT =",Y(30)
!      ENDIF
!      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
!       WRITE(0,*)"XISGUT =",Y(31)," XISINP=",XISINP," ERR =",
!     .      (Y(31)-XISINP)**2/(1.D6+XISINP**2)
!      ELSE
!       WRITE(0,*)"XISGUT =",Y(31)
!      ENDIF
!      WRITE(0,*)"MUPGUT =",Y(32)," MUPINP  =",MUPINP," ERR =",
!     .      (Y(32)-MUPINP)**2/(1.D2+MUPINP**2)
!      WRITE(0,*)"MSPGUT =",Y(33)," MSPINP =",MSPINP," ERR =",
!     .      (Y(33)-MSPINP)**2/(1.D4+MSPINP**2)
!      WRITE(0,*)"M3HGUT =",Y(34)," MSPINP =",M3HINP," ERR =",
!     .      (Y(34)-M3HINP)**2/(1.D4+M3HINP**2)
!      WRITE(0,*)""

      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"IFAIL =",IFAIL
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=12
       RETURN
      ENDIF

* GUTEST:

      GUTEST=0d0
      GUTEST=GUTEST+(Y(9)-M1INP)**2/(1.D2+M1INP**2)
      GUTEST=GUTEST+(Y(10)-M2INP)**2/(1.D2+M2INP**2)
      GUTEST=GUTEST+(Y(11)-M3INP)**2/(1.D2+M3INP**2)
      GUTEST=GUTEST+(Y(12)-ALINP)**2/(1.D2+ALINP**2)
      GUTEST=GUTEST+(Y(13)-AKINP)**2/(1.D2+AKINP**2)
      DO I=14,16
       GUTEST=GUTEST+(Y(I)-A0)**2/(1.D2+A0**2)
      ENDDO
      IF(MAFLAG.NE.-5)THEN
       GUTEST=GUTEST+(Y(17)-MH1INP)**2/(1.D4+MH1INP**2)
       GUTEST=GUTEST+(Y(18)-MH2INP)**2/(1.D4+MH2INP**2)
      ENDIF
      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
       GUTEST=GUTEST+(Y(19)-MSINP)**2/(1.D4+MSINP**2)
      ENDIF
      DO I=20,29
       GUTEST=GUTEST+(Y(I)-M0**2)**2/(1.D4+M0**4)
      ENDDO
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
       GUTEST=GUTEST+(Y(30)-XIFINP)**2/(1.D4+XIFINP**2)
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
       GUTEST=GUTEST+(Y(31)-XISINP)**2/(1.D6+XISINP**2)
      ENDIF
      GUTEST=GUTEST+(Y(32)-MUPINP)**2/(1.D2+MUPINP**2)
      GUTEST=GUTEST+(Y(33)-MSPINP)**2/(1.D4+MSPINP**2)
      GUTEST=GUTEST+(Y(34)-M3HINP)**2/(1.D4+M3HINP**2)
      GUTEST=GUTEST+(Y(35)-A0)**2/(1.D2+A0**2)
!      WRITE(0,*)"GUTEST =",GUTEST
!      WRITE(0,*)""
!      WRITE(0,*)""

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
      MH1GUT=Y(17)
      MH2GUT=Y(18)
      MSGUT=Y(19)
      IF(MAFLAG.EQ.-5)THEN
       MUT=MUT+Y(17)
       MUM=MUT/IM
       MDT=MDT+Y(18)
       MDM=MDT/IM
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-5)THEN
       MST=MST+Y(19)
       MSM=MST/IM
      ENDIF
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
      AMUGUT=Y(35)

      END
