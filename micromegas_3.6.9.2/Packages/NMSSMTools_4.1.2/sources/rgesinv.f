      SUBROUTINE RGESINV(PAR,IFAIL)

*   Subroutine to integrate the 2-loop RGEs for all soft terms
*   from MGUT down to Q2 = M_SUSY, through a CALL of the
*   subroutine ODEINTS that is part of the file integs.f
*
*   At MGUT:
*   MHD^2 = MDM (MAFLAG=-5), MHDINP (MAFLAG=/=-5)
*   MHU^2 = MUM (MAFLAG=-5), MHUINP (MAFLAG=/=-5)
*   MS^2 = MSM (MAFLAG=-1,-3,-5), MSINP (MAFLAG=-2,-4)
*   XIS = XISGUT (MAFLAG=-2,-4), XISINP (MAFLAG=-1,-3)
*   XIF = XIFGUT (MAFLAG=-3,-4), XIFINP (MAFLAG=-1,-2)
*
*   It uses COMMON/INPPAR, /GUTPAR and /GUTEXT for the soft terms,
*   COMMON/GUTCOUP for the gauge/Yukawa couplings at MGUT.
*
*   The output (soft terms at Q2) is written into PAR(*),
*   COMMON/SUSYMH and /SUSYEXT
*
***********************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,NN,I,IM,OMGFLAG,MAFLAG
      PARAMETER (NN=35)

      DOUBLE PRECISION PAR(*),EPS,X1,X2,Y(NN),PI
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION Q2,COEF,M0,M12,A0,MH1S,MH2S,MSS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT
      DOUBLE PRECISION HTOPGUT,HBOTGUT,HTAUGUT,MGUT
      DOUBLE PRECISION M1INP,M2INP,M3INP,MH2INP,MH1INP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT
      DOUBLE PRECISION ATAUGUT,AMUGUT,MH1GUT,MH2GUT,MSGUT,MQ3GUT
      DOUBLE PRECISION MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT
      DOUBLE PRECISION ME3GUT,MLGUT,MEGUT,MUM,MDM,MSM,MUT,MDT,MST
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/RENSCALE/Q2
      COMMON/MGUT/MGUT
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/INPPAR/M1INP,M2INP,M3INP,MH2INP,MH1INP,ALINP,AKINP,
     .      XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     .      HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     .      ATAUGUT,AMUGUT,MH1GUT,MH2GUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     .      MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/SUSYMH/MH1S,MH2S,MSS
      COMMON/MSAVE/MUM,MDM,MSM,MUT,MDT,MST,IM

      EXTERNAL DERIVSS,RKQSS

      EPS=1.D-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

* Definition of the couplings squared Y(I) at MGUT

      Y(1)=g1GUT
      Y(2)=g2GUT
      Y(3)=g3GUT
      Y(4)=LGUT
* NOTE: KGUT = KAPPA, NOT KAPPA**2
      Y(5)=KGUT
      Y(6)=HTOPGUT
      Y(7)=HBOTGUT
      Y(8)=HTAUGUT
      
* Input values for the soft terms at MGUT:
      
      Y(9)=M1INP
      Y(10)=M2INP
      Y(11)=M3INP
      Y(12)=ALINP
      Y(13)=AKINP
      Y(14)=A0
      Y(15)=A0
      Y(16)=A0
      IF(MAFLAG.EQ.-5)THEN
       Y(17)=MUM
       Y(18)=MDM
      ELSE
       Y(17)=MH1INP
       Y(18)=MH2INP
      ENDIF
      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-5)THEN
       Y(19)=MSM
      ELSE
       Y(19)=MSINP
      ENDIF
      DO I=20,29
        Y(I)=M0**2
      ENDDO
      IF(MAFLAG.EQ.-3 .OR. MAFLAG.EQ.-4)THEN
       Y(30)=XIFGUT
      ELSE
       Y(30)=XIFINP
      ENDIF
      IF(MAFLAG.EQ.-2 .OR. MAFLAG.EQ.-4)THEN
       Y(31)=XISGUT
      ELSE
       Y(31)=XISINP
      ENDIF
      Y(32)=MUPINP
      Y(33)=MSPINP
      Y(34)=M3HINP
      Y(35)=A0

      X1=COEF*DLOG(MGUT**2/Q2)
      X2=0d0

!      WRITE(0,*)"CALL RGESINV"
!      WRITE(0,*)""
!      WRITE(0,*)"MGUT =",MGUT
!      WRITE(0,*)"G1GUT =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2GUT =",Y(2)
!      WRITE(0,*)"G3GUT =",Y(3)
!      WRITE(0,*)"LGUT =",DSQRT(Y(4))
!      WRITE(0,*)"KGUT =",Y(5)
!      WRITE(0,*)"HTGUT =",Y(6)
!      WRITE(0,*)"HBGUT =",Y(7)
!      WRITE(0,*)"HLGUT =",Y(8)
!      WRITE(0,*)"M1GUT =",Y(9)
!      WRITE(0,*)"M2GUT =",Y(10)
!      WRITE(0,*)"M3GUT =",Y(11)
!      WRITE(0,*)"ALGUT =",Y(12)
!      WRITE(0,*)"AKGUT =",Y(13)
!      WRITE(0,*)"ATOPGUT =",Y(14)
!      WRITE(0,*)"ABOTGUT =",Y(15)
!      WRITE(0,*)"ATAUGUT =",Y(16)
!      WRITE(0,*)"AMUGUT =",Y(35)
!      WRITE(0,*)"MH1GUT =",Y(17)
!      WRITE(0,*)"MH2GUT =",Y(18)
!      WRITE(0,*)"MSGUT =",Y(19)
!      WRITE(0,*)"MQ3GUT =",Y(20)
!      WRITE(0,*)"MU3GUT =",Y(21)
!      WRITE(0,*)"MD3GUT =",Y(22)
!      WRITE(0,*)"MQGUT =",Y(23)
!      WRITE(0,*)"MUGUT =",Y(24)
!      WRITE(0,*)"MDGUT =",Y(25)
!      WRITE(0,*)"ML3GUT =",Y(26)
!      WRITE(0,*)"ME3GUT =",Y(27)
!      WRITE(0,*)"MLGUT =",Y(28)
!      WRITE(0,*)"MEGUT =",Y(29)
!      WRITE(0,*)"XIFGUT =",Y(30)
!      WRITE(0,*)"XISGUT =",Y(31)
!      WRITE(0,*)"MUPGUT =",Y(32)
!      WRITE(0,*)"MSPGUT =",Y(33)
!      WRITE(0,*)"M3HGUT =",Y(34)
!      WRITE(0,*)""

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)

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

      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"IFAIL =",IFAIL
!       WRITE(0,*)""
!       WRITE(0,*)""
       IFAIL=13
       RETURN
      ENDIF
!      WRITE(0,*)""

* Y(5) = KAPPA, NOT KAPPA**2

      IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)PAR(2)=Y(5)

* SOFT TERMS AT THE SUSY SCALE

      PAR(5)=Y(12)
      PAR(6)=Y(13)
      PAR(7)=Y(20)
      PAR(8)=Y(21)
      PAR(9)=Y(22)
      PAR(10)=Y(26)
      PAR(11)=Y(27)
      PAR(12)=Y(14)
      PAR(13)=Y(15)
      PAR(14)=Y(16)
      PAR(15)=Y(23)
      PAR(16)=Y(24)
      PAR(17)=Y(25)
      PAR(18)=Y(28)
      PAR(19)=Y(29)
      PAR(20)=Y(9)
      PAR(21)=Y(10)
      PAR(22)=Y(11)
      PAR(25)=Y(35)

* MH1S, MH2S AND MSS at Q2 are stored in COMMON/SUSYMH:
      
      MH1S=Y(17)
      MH2S=Y(18)
      MSS=Y(19)

* EXT parameters at Q2, stored in COMMON/SUSYEXT:

      XIFSUSY=Y(30)
      XISSUSY=Y(31)
      MUPSUSY=Y(32)
      MSPSUSY=Y(33)
      M3HSUSY=Y(34)
      
      END
