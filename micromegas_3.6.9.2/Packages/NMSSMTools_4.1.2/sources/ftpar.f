      SUBROUTINE FTPAR(PAR,FLAG) 

************************************************************************   
*   Subroutine extracting the fine-tuning
*
*   At the SUSY scale FTSUSY(i)=d(ln MZ^2)/d(ln PSi^2)
*
*     FTSUSY(1) --> PS=m_Hu
*     FTSUSY(2) --> PS=m_Hd
*     FTSUSY(3) --> PS=m_S
*     FTSUSY(4) --> PS=A_lambda
*     FTSUSY(5) --> PS=A_kappa
*     FTSUSY(6) --> PS=XiF
*     FTSUSY(7) --> PS=XiS
*     FTSUSY(8) --> PS=MUP
*     FTSUSY(9) --> PS=MSP
*     FTSUSY(10) --> PS=M_3H
*     FTSUSY(11) --> PS=lambda
*     FTSUSY(12) --> PS=kappa
*     FTSUSY(13) --> PS=htop
*     FTSUSY(14) --> PS=g
*     FTSUSY(15) --> MAX(FTSUSY(i))
*     FTSUSY(16) --> i FOR FTMAX
*
*   Based on the 3 minimization equations:
*
*     (1)=mhu+mu**2+mz**2*(lh+(tb**2-1)/2)/(1+tb**2)
*        -(mu*b+m3+mu*mup+xfh)/tb-ct*mz**2*tb**2/(1+tb**2)
*     (2)=mhd+mu**2+mz**2*(lh*tb**2+(1-tb**2)/2)/(1+tb**2)
*        -(mu*b+m3+mu*mup+xfh)*tb
*     (3)=ms+msp+mup**2+2*kh*xfh+kh*mu*(ak+2*kh*mu+3*mup)
*        +mz**2*lh*(1-(2*kh*mu+al+mup)/mu*tb/(1+tb**2))
*        +(xsh+xfh*mup)/mu
*
*     where lh=l**2/g**2, kh=k/l, xfh=l*xif, xsh = l*xis, b=al+kh*mu
*     ct0 = (3*ht**4)/(8*pi**2*g**2)*ln(mt**2/QSTSB)
* 
*   At the GUT scale FTGUT(i)=d(ln MZ^2)/d(ln PGi^2)
*
*     FTGUT(1) --> PG=m_Hu
*     FTGUT(2) --> PG=m_Hd
*     FTGUT(3) --> PG=m_S
*     FTGUT(4) --> PG=m0
*     FTGUT(5) --> PG=M1
*     FTGUT(6) --> PG=M2
*     FTGUT(7) --> PG=M3
*     FTGUT(8) --> PG=M12
*     FTGUT(9) --> PG=A_lambda
*     FTGUT(10) --> PG=A_kappa
*     FTGUT(11) --> PG=A0
*     FTGUT(12) --> PG=XiF
*     FTGUT(13) --> PG=XiS
*     FTGUT(14) --> PG=MUP
*     FTGUT(15) --> PG=MSP
*     FTGUT(16) --> PG=M3H
*     FTGUT(17) --> PG=lambda
*     FTGUT(18) --> PG=kappa
*     FTGUT(19) --> PG=htop
*     FTGUT(20) --> PG=g
*     FTGUT(21) --> PG=MGUT
*     FTGUT(22) --> MAX(FTGUT(i))
*     FTGUT(23) --> i FOR FTMAX
*
*   JACG(j,i) is the Jacobian defined as d(ln PSj^2)/d(ln PGi^2)
* 
*   At the MESS scale FTMES(i)=d(ln MZ^2)/d(ln PMi^2)
*
*     FTMES(1) --> PM=M_susyeff
*     FTMES(2) --> PM=m_S
*     FTMES(3) --> PM=Delta_H
*     FTMES(4) --> PM=A_lambda
*     FTMES(5) --> PM=XiF
*     FTMES(6) --> PM=XiS
*     FTMES(7) --> PM=MUP
*     FTMES(8) --> PM=MSP
*     FTMES(9) --> PM=lambda
*     FTMES(10) --> PM=kappa
*     FTMES(11) --> PM=htop
*     FTMES(12) --> PM=g1
*     FTMES(13) --> PM=g2
*     FTMES(14) --> PM=g3
*     FTMES(15) --> PM=MMESS
*     FTMES(16) --> MAX(FTGUT(i))
*     FTMES(17) --> i FOR FTMAX
*
*   JACM(j,i) is the Jacobian defined as d(ln PSj^2)/d(ln PMi^2)
* 
*
************************************************************************
 
      IMPLICIT NONE 

      INTEGER I,J,FLAG,M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER MSFLAG,AKFLAG,ALFLAG,NSUSY,NGUT,NMES
      PARAMETER (NSUSY=14,NGUT=21,NMES=15)

      DOUBLE PRECISION PAR(*),FTSUSY(NSUSY+2),FTGUT(NGUT+2)
      DOUBLE PRECISION PS(NSUSY),PG(NGUT),JACG(NSUSY,NGUT)
      DOUBLE PRECISION PM(NMES),JACM(NSUSY,NMES),FTMES(NMES+2)
      DOUBLE PRECISION L,K,TB,TB2,MU,AL,AK,PI,ALS,B
      DOUBLE PRECISION LH,KH,XFH,XSH,CT
      DOUBLE PRECISION A1,A2,A3,D1Z,D2Z,D3Z
      DOUBLE PRECISION D1T,D2T,D3T,D1M,D2M,D3M,D
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW      
      DOUBLE PRECISION G,G1,G2,G3,HT,HB,HL
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION MHU2,MHD2,MS2,M0,M12,A0
      DOUBLE PRECISION QSTSB,HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT
      DOUBLE PRECISION HTGUT,HBGUT,HLGUT,MGUT
      DOUBLE PRECISION M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT
      DOUBLE PRECISION ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT
      DOUBLE PRECISION MU3GUT,MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT
      DOUBLE PRECISION ME3GUT,MLGUT,MEGUT
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP
      DOUBLE PRECISION MSPINP,DELHINP,MSUSYEFF,MMESS,N5
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES
      DOUBLE PRECISION HTMES,HBMES,HLMES,XIFMES,XISMES
      DOUBLE PRECISION MSMES,MUPMES,MSPMES,M3HMES

      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYCOUP/g1,g2,g3,HT,HB,HL
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/SUSYMH/MHU2,MHD2,MS2
      COMMON/SOFTGUT/M0,M12,A0
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,
     .      HBGUT,HLGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     .      ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     .      MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     .      MSFLAG,AKFLAG,ALFLAG
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/ALSHIFT/ALS,B
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP,DELHINP
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES

      PI=4d0*DATAN(1d0)

      DO I=1,NSUSY+2
       FTSUSY(I)=0d0
      ENDDO

      L=PAR(1)
      K=PAR(2)
      TB=PAR(3)
      MU=PAR(4)
      AL=PAR(5)
      AK=PAR(6)
      G=(G1+G2)/2d0

      LH=L**2/G
      KH=K/L
      XFH=L*XIF
      XSH=L*XIS
      TB2=TB**2
      CT=3d0*HT**4/(8d0*PI**2*G)*DLOG(MTOPQ**2/QSTSB)
    
*   DiT = d(i)/d(tb)

      D1T = 2d0*MZ**2*TB/(1d0+TB2)**2*(1d0-LH-CT)
     .    + (MU*B+M3H+MU*MUP+XFH)/TB2
      D2T = 2d0*MZ**2*TB/(1d0+TB2)**2*(LH-1d0)
     .    - (MU*B+M3H+MU*MUP+XFH)
      D3T = LH*MZ**2*(TB2-1d0)/(1d0+TB2)**2
     .     *(B+KH*MU+MUP)/MU

*   DiM = d(i)/d(mu)

      D1M = 2d0*MU-(B+KH*MU+MUP)/TB
      D2M = 2d0*MU-(B+KH*MU+MUP)*TB
      D3M = KH*AK+4d0*KH**2*MU+3d0*KH*MUP
     .    + LH*MZ**2*TB/(1d0+TB2)*(AL+MUP)/MU**2
     .    - (XSH+XFH*MUP)/MU**2

*   DiZ = d(i)/d(MZ^2)

      D1Z = (LH+(TB2-1d0)/2d0-CT*TB2)/(1d0+TB2)
      D2Z = (LH*TB2+(1d0-TB2)/2d0)/(1d0+TB2)
      D3Z = LH*(1d0-TB/(1d0+TB2)*(B+KH*MU+MUP)/MU)

*   Ai=d(j)/d(tb)*d(k)/d(mu)-d(k)/d(tb)*d(j)/d(mu)

      A1=D2T*D3M-D2M*D3T
      A2=D3T*D1M-D3M*D1T
      A3=D1T*D2M-D1M*D2T

*   Determinant of the linear system

      D=A1*D1Z+A2*D2Z+A3*D3Z

*   Fine-tuning measures at the SUSY scale

      FTSUSY(1)=MHU2/MZ**2*A1/D
      FTSUSY(2)=MHD2/MZ**2*A2/D
      FTSUSY(3)=MS2/MZ**2*A3/D
      FTSUSY(4)=AL/2d0/MZ**2*(-MU/TB*A1-MU*TB*A2
     .          -MZ**2*LH*TB/(1d0+TB2)/MU*A3)/D
      FTSUSY(5)=AK/2d0/MZ**2*KH*MU*A3/D
      FTSUSY(6)=XFH/MZ**2*(-A1/TB-A2*TB+(2d0*KH+MUP/MU)*A3)/D
      FTSUSY(7)=3d0/2d0*XSH/MZ**2*A3/MU/D
      FTSUSY(8)=MUP/2d0/MZ**2*(-MU/TB*A1-MU*TB*A2+(2d0*MUP
     .          +3d0*KH*MU+XFH/MU-LH*MZ**2/MU*TB/(1d0+TB2)))/D
      FTSUSY(9)=MSP/MZ**2*A3/D
      FTSUSY(10)=M3H/MZ**2*(-A1/TB-A2*TB)/D
      FTSUSY(11)=(LH*(A1+A2*TB2)/(1d0+TB2)+D3Z*A3)/D
      FTSUSY(12)=KH/2d0/MZ**2*(-MU**2/TB*A1-MU**2*TB*A2
     .          +(2d0*XFH+AK*MU+4d0*KH*MU**2+3d0*MU*MUP
     .          -2d0*LH*MZ**2*TB/(1d0+TB2))*A3)/D
      FTSUSY(13)=-2d0*CT*TB2/(1d0+TB2)*A1/D
      FTSUSY(14)=-FTSUSY(11)+CT*TB2/(1d0+TB2)*A1/D
      FTSUSY(11)=FTSUSY(11)-FTSUSY(12)

      FTSUSY(NSUSY+1)=DABS(FTSUSY(1))
      FTSUSY(NSUSY+2)=1d0
      DO I=2,NSUSY
       IF(DABS(FTSUSY(I)).GE.FTSUSY(NSUSY+1))THEN
        FTSUSY(NSUSY+1)=DABS(FTSUSY(I))
        FTSUSY(NSUSY+2)=DFLOAT(I)
       ENDIF
      ENDDO

*   mSUGRA boundary conditions

      IF(FLAG.EQ.1)THEN

       DO I=1,NGUT+2
        FTGUT(I)=0d0
       ENDDO

*   GUT scale parameters

       IF(MHUFLAG.EQ.1)THEN
        PG(1)=MHUGUT
       ELSE
        PG(1)=0d0
       ENDIF
       IF(MHDFLAG.EQ.1)THEN
        PG(2)=MHDGUT
       ELSE
        PG(2)=0d0
       ENDIF
       IF(MSFLAG.EQ.1)THEN
        PG(3)=MSGUT
       ELSE
        PG(3)=0d0
       ENDIF
       PG(4)=M0**2
       IF(M1FLAG.EQ.1)THEN
        PG(5)=M1GUT
       ELSE
        PG(5)=0d0
       ENDIF
       IF(M2FLAG.EQ.1)THEN
        PG(6)=M2GUT
       ELSE
        PG(6)=0d0
       ENDIF
       IF(M3FLAG.EQ.1)THEN
        PG(7)=M3GUT
       ELSE
        PG(7)=0d0
       ENDIF
       IF(M1FLAG*M2FLAG*M3FLAG.EQ.0)THEN
        PG(8)=M12
       ELSE
        PG(8)=0d0
       ENDIF
       IF(ALFLAG.EQ.1)THEN
        PG(9)=ALGUT
       ELSE
        PG(9)=0d0
       ENDIF
       IF(AKFLAG.EQ.1)THEN
        PG(10)=AKGUT
       ELSE
        PG(10)=0d0
       ENDIF
       PG(11)=A0
       PG(12)=XIFGUT
       PG(13)=XISGUT
       PG(14)=MUPGUT
       PG(15)=MSPGUT
       PG(16)=M3HGUT
       PG(17)=LGUT
       PG(18)=KGUT
       PG(19)=HTGUT
       PG(20)=G2GUT
       PG(21)=MGUT**2

       CALL RGESVAR(NGUT,NSUSY,PG,PS)

       IF(MHU2.NE.0d0)
     .  FTSUSY(1)=PS(1)/MHU2*FTSUSY(1)
       IF(MHD2.NE.0d0)
     .  FTSUSY(2)=PS(2)/MHD2*FTSUSY(2)
       IF(MS2.NE.0d0)
     .  FTSUSY(3)=PS(3)/MS2*FTSUSY(3)
       IF(AL.NE.0d0)
     .  FTSUSY(4)=PS(4)/AL*FTSUSY(4)
       IF(AK.NE.0d0)
     .  FTSUSY(5)=PS(5)/AK*FTSUSY(5)
       IF(XIF.NE.0d0)
     .  FTSUSY(6)=PS(6)/XIF*FTSUSY(6)
       IF(XIS.NE.0d0)
     .  FTSUSY(7)=PS(7)/XIS*FTSUSY(7)
       IF(MUP.NE.0d0)
     .  FTSUSY(8)=PS(8)/MUP*FTSUSY(8)
       IF(MSP.NE.0d0)
     .  FTSUSY(9)=PS(9)/MSP*FTSUSY(9)
       IF(M3H.NE.0d0)
     .  FTSUSY(10)=PS(10)/M3H*FTSUSY(10)
       FTSUSY(11)=PS(11)/L**2*FTSUSY(11)
       IF(K.NE.0d0)
     .  FTSUSY(12)=PS(12)/K*FTSUSY(12)
       FTSUSY(13)=PS(13)/HT**2*FTSUSY(13)
       FTSUSY(14)=PS(14)/G*FTSUSY(14)

       CALL JACOBIAN(NGUT,NSUSY,PG,PS,JACG)

       DO I=1,NGUT
        FTGUT(I)=0d0
        DO J=1,NSUSY
         FTGUT(I)=FTGUT(I)+FTSUSY(J)*JACG(J,I)
        ENDDO
       ENDDO

       FTGUT(NGUT+1)=DABS(FTGUT(1))
       FTGUT(NGUT+2)=1d0
       DO I=2,NGUT-2
        IF(DABS(FTGUT(I)).GE.FTGUT(NGUT+1))THEN
         FTGUT(NGUT+1)=DABS(FTGUT(I))
         FTGUT(NGUT+2)=DFLOAT(I)
        ENDIF
       ENDDO

      ENDIF

*   GMSB boundary conditions

      IF(FLAG.EQ.2)THEN

       DO I=1,NMES+2
        FTMES(I)=0d0
       ENDDO

*   MES scale parameters

       PM(1)=MSUSYEFF
       PM(2)=MSMES
       PM(3)=DELHINP
       PM(4)=ALINP
       PM(5)=XIFMES
       PM(6)=XISMES
       PM(7)=MUPMES
       PM(8)=MSPMES
       PM(9)=LMES
       PM(10)=KMES
       PM(11)=HTMES
       PM(12)=G1MES
       PM(13)=G2MES
       PM(14)=G3MES
       PM(15)=MMESS**2

       CALL RGESVARGM(NMES,NSUSY,PM,PS,N5)

       IF(MHU2.NE.0d0)
     .  FTSUSY(1)=PS(1)/MHU2*FTSUSY(1)
       IF(MHD2.NE.0d0)
     .  FTSUSY(2)=PS(2)/MHD2*FTSUSY(2)
       IF(MS2.NE.0d0)
     .  FTSUSY(3)=PS(3)/MS2*FTSUSY(3)
       IF(AL.NE.0d0)
     .  FTSUSY(4)=PS(4)/AL*FTSUSY(4)
       IF(AK.NE.0d0)
     .  FTSUSY(5)=PS(5)/AK*FTSUSY(5)
       IF(XIF.NE.0d0)
     .  FTSUSY(6)=PS(6)/XIF*FTSUSY(6)
       IF(XIS.NE.0d0)
     .  FTSUSY(7)=PS(7)/XIS*FTSUSY(7)
       IF(MUP.NE.0d0)
     .  FTSUSY(8)=PS(8)/MUP*FTSUSY(8)
       IF(MSP.NE.0d0)
     .  FTSUSY(9)=PS(9)/MSP*FTSUSY(9)
       IF(M3H.NE.0d0)
     .  FTSUSY(10)=PS(10)/M3H*FTSUSY(10)
       FTSUSY(11)=PS(11)/L**2*FTSUSY(11)
       IF(K.NE.0d0)
     .  FTSUSY(12)=PS(12)/K*FTSUSY(12)
       FTSUSY(13)=PS(13)/HT**2*FTSUSY(13)
       FTSUSY(14)=PS(14)/G*FTSUSY(14)

       CALL JACOBIANGM(NMES,NSUSY,PM,PS,N5,JACM)

       DO I=1,NMES
        FTMES(I)=0d0
        DO J=1,NSUSY
         FTMES(I)=FTMES(I)+FTSUSY(J)*JACM(J,I)
        ENDDO
       ENDDO

       FTMES(NMES+1)=DABS(FTMES(1))
       FTMES(NMES+2)=1d0
       DO I=2,NMES
        IF(DABS(FTMES(I)).GE.FTMES(NMES+1))THEN
         FTMES(NMES+1)=DABS(FTMES(I))
         FTMES(NMES+2)=DFLOAT(I)
        ENDIF
       ENDDO

      ENDIF

      END


      SUBROUTINE JACOBIAN(M,N,PG,PS,JACG)

      IMPLICIT NONE 

      INTEGER N,M,I,J
      DOUBLE PRECISION JACG(N,M),PG(M),PS(N),PP(N),T,H,EPS

      EPS=1.d-4

      DO J=1,M
       T=PG(J)
       IF(T.EQ.0d0)THEN
        DO I=1,N
         JACG(I,J)=0d0
        ENDDO
       ELSE
        H=EPS*DABS(T)
        PG(J)=T+H
        CALL RGESVAR(M,N,PG,PP)
        PG(J)=T
        DO I=1,N
         IF(PS(I).EQ.0d0)THEN
          JACG(I,J)=0d0
         ELSE
          JACG(I,J)=PG(J)/PS(I)*PP(I)/(2d0*H)
         ENDIF
        ENDDO
        PG(J)=T-H
        CALL RGESVAR(M,N,PG,PP)
        PG(J)=T
        DO I=1,N
         IF(PS(I).EQ.0d0)THEN
          JACG(I,J)=0d0
         ELSE
          JACG(I,J)=JACG(I,J)-PG(J)/PS(I)*PP(I)/(2d0*H)
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      DO J=1,M
       JACG(4,J)=JACG(4,J)*2d0
       JACG(5,J)=JACG(5,J)*2d0
       JACG(7,J)=JACG(7,J)*2d0/3d0
       JACG(8,J)=JACG(8,J)*2d0
       JACG(12,J)=JACG(12,J)*2d0
      ENDDO
      DO I=1,N
       DO J=5,11
        JACG(I,J)=JACG(I,J)/2d0
       ENDDO
       JACG(I,13)=JACG(I,14)*3d0/2d0
       JACG(I,14)=JACG(I,14)/2d0
       JACG(I,18)=JACG(I,18)/2d0
      ENDDO

      END


      SUBROUTINE RGESVAR(M,N,PG,PS)

      IMPLICIT NONE

      INTEGER I,IFAIL,M,N,M1,M2,M3,MHD,MHU,MS,AK,AL,NN
      PARAMETER (NN=35)

      DOUBLE PRECISION PG(M),PS(N),EPS,X1,X2,Y(NN),PI
      DOUBLE PRECISION G1,G2,G3,L,K,HT,HB,HL,Q2,COEF

      COMMON/GUTCOUP/G1,G2,G3,L,K,HT,HB,HL
      COMMON/RENSCALE/Q2
      COMMON/SCANFLAGS/M1,M2,M3,MHD,MHU,MS,AK,AL
 
      EXTERNAL DERIVSS,RKQSS

      EPS=1.d-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

      DO I=1,N
       PS(I)=0d0
      ENDDO

      Y(1)=3d0/5d0*PG(20)
      Y(2)=PG(20)
      Y(3)=G3-G2+PG(20)
      Y(4)=PG(17)
      Y(5)=PG(18)
      Y(6)=PG(19)
      Y(7)=HB
      Y(8)=HL
      IF(M1.EQ.1)THEN
       Y(9)=PG(5)
      ELSE
       Y(9)=PG(8)
      ENDIF
      IF(M2.EQ.1)THEN
       Y(10)=PG(6)
      ELSE
       Y(10)=PG(8)
      ENDIF
      IF(M3.EQ.1)THEN
       Y(11)=PG(7)
      ELSE
       Y(11)=PG(8)
      ENDIF
      IF(AL.EQ.1)THEN
       Y(12)=PG(9)
      ELSE
       Y(12)=PG(11)
      ENDIF
      IF(AK.EQ.1)THEN
       Y(13)=PG(10)
      ELSE
       Y(13)=PG(11)
      ENDIF
      Y(14)=PG(11)
      Y(15)=PG(11)
      Y(16)=PG(11)
      IF(MHU.EQ.1)THEN
       Y(17)=PG(1)
      ELSE
       Y(17)=PG(4)
      ENDIF
      IF(MHD.EQ.1)THEN
       Y(18)=PG(2)
      ELSE
       Y(18)=PG(4)
      ENDIF
      IF(MS.EQ.1)THEN
       Y(19)=PG(3)
      ELSE
       Y(19)=PG(4)
      ENDIF
      DO I=20,29
        Y(I)=PG(4)
      ENDDO
      Y(30)=PG(12)
      Y(31)=PG(13)
      Y(32)=PG(14)
      Y(33)=PG(15)
      Y(34)=PG(16)
      Y(35)=PG(11)

      X1=COEF*DLOG(PG(21)/Q2)
      X2=0d0

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)

      IF(IFAIL.NE.0)THEN
       RETURN
      ENDIF

      PS(1)=Y(17)
      PS(2)=Y(18)
      PS(3)=Y(19)
      PS(4)=Y(12)
      PS(5)=Y(13)
      PS(6)=Y(30)
      PS(7)=Y(31)
      PS(8)=Y(32)
      PS(9)=Y(33)
      PS(10)=Y(34)
      PS(11)=Y(4)
      PS(12)=Y(5)
      PS(13)=Y(6)
      PS(14)=(Y(1)+Y(2))/2d0
      
      END


      SUBROUTINE JACOBIANGM(M,N,PM,PS,N5,JACM)

      IMPLICIT NONE 

      INTEGER N,M,I,J
      DOUBLE PRECISION JACM(N,M),PM(M),PS(N),PP(N),N5,T,H,EPS

      EPS=1.d-4

      DO J=1,M
       T=PM(J)
       IF(T.EQ.0d0)THEN
        DO I=1,N
         JACM(I,J)=0d0
        ENDDO
       ELSE
        H=EPS*DABS(T)
        PM(J)=T+H
        CALL RGESVARGM(M,N,PM,PP,N5)
        PM(J)=T
        DO I=1,N
         IF(PS(I).EQ.0d0)THEN
          JACM(I,J)=0d0
         ELSE
          JACM(I,J)=PM(J)/PS(I)*PP(I)/(2d0*H)
         ENDIF
        ENDDO
        PM(J)=T-H
        CALL RGESVARGM(M,N,PM,PP,N5)
        PM(J)=T
        DO I=1,N
         IF(PS(I).EQ.0d0)THEN
          JACM(I,J)=0d0
         ELSE
          JACM(I,J)=JACM(I,J)-PM(J)/PS(I)*PP(I)/(2d0*H)
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      DO J=1,M
       JACM(4,J)=JACM(4,J)*2d0
       JACM(5,J)=JACM(5,J)*2d0
       JACM(7,J)=JACM(7,J)*2d0/3d0
       JACM(8,J)=JACM(8,J)*2d0
       JACM(12,J)=JACM(12,J)*2d0
      ENDDO
      DO I=1,N
       JACM(I,1)=JACM(I,1)/2d0
       JACM(I,4)=JACM(I,4)/2d0
       JACM(I,6)=JACM(I,6)*3d0/2d0
       JACM(I,7)=JACM(I,7)/2d0
       JACM(I,10)=JACM(I,10)/2d0
       JACM(I,15)=JACM(I,15)/2d0
      ENDDO

      END


      SUBROUTINE RGESVARGM(M,N,PM,PS,N5)

      IMPLICIT NONE

      INTEGER I,IFAIL,M,N,NN
      PARAMETER (NN=35)

      DOUBLE PRECISION PM(M),PS(N),N5,EPS,X1,X2,Y(NN)
      DOUBLE PRECISION G1,G2,G3,L,K,HT,HB,HL,Q2,COEF
      DOUBLE PRECISION X,F1,F2,ALP1,ALP2,ALP3,SP2,PI

      COMMON/MESCOUP/G1,G2,G3,L,K,HT,HB,HL
      COMMON/RENSCALE/Q2
 
      EXTERNAL DERIVSS,RKQSS

      EPS=1.d-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

      DO I=1,N
       PS(I)=0d0
      ENDDO

      Y(1)=PM(12)
      Y(2)=PM(13)
      Y(3)=PM(14)
      Y(4)=PM(9)
      Y(5)=PM(10)
      Y(6)=PM(11)
      Y(7)=HB
      Y(8)=HL

      X=PM(1)/PM(15)
      IF(X.GE.1.D-7) THEN
        IF(X.LT.1d0) THEN
          F1=((1d0+X)*DLOG(1d0+X)+(1d0-X)*DLOG(1d0-X))/X**2
          F2=(1d0+X)/X**2*(DLOG(1d0+X)-2d0*SP2(X/(1d0+X))
     .       +SP2(2d0*X/(1d0+X))/2d0)
     .       +(1d0-X)/X**2*(DLOG(1d0-X)-2d0*SP2(-X/(1d0-X))
     .       +SP2(-2d0*X/(1d0-X))/2d0)
        ELSE
          F1=0d0
          F2=0d0
        ENDIF
      ELSE
        F1=1d0
        F2=1d0
      ENDIF

      ALP1=PM(12)/(4d0*PI)
      ALP2=PM(13)/(4d0*PI)
      ALP3=PM(14)/(4d0*PI)
      
      Y(9)=5d0/3d0*ALP1*N5*PM(1)*F1/(4d0*PI)
      Y(10)=ALP2*N5*PM(1)*F1/(4d0*PI)
      Y(11)=ALP3*N5*PM(1)*F1/(4d0*PI)
      Y(12)=PM(4)
      IF(PM(10).NE.0d0)THEN
       Y(13)=PM(4)*3d0
      ELSE
       Y(13)=0d0
      ENDIF
      Y(14)=0d0
      Y(15)=0d0
      Y(16)=0d0
      Y(17)=COEF*((5d0/6d0*ALP1**2+3d0/2d0*ALP2**2)*N5*F2
     .     -COEF*PM(9)*PM(3))*PM(1)**2
      Y(18)=COEF*((5d0/6d0*ALP1**2+3d0/2d0*ALP2**2)*N5*F2
     .     -COEF*PM(9)*PM(3))*PM(1)**2
      Y(19)=PM(2)
      Y(20)=COEF*(5d0/54d0*ALP1**2+3d0/2d0*ALP2**2
     .     +8d0/3d0*ALP3**2)*F2*N5*PM(1)**2
      Y(21)=COEF*(40d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*PM(1)**2
      Y(22)=COEF*(10d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*PM(1)**2
      Y(23)=COEF*(5d0/54d0*ALP1**2+3d0/2d0*ALP2**2
     .     +8d0/3d0*ALP3**2)*F2*N5*PM(1)**2
      Y(24)=COEF*(40d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*PM(1)**2
      Y(25)=COEF*(10d0/27d0*ALP1**2
     .     +8d0/3d0*ALP3**2)*F2*N5*PM(1)**2
      Y(26)=COEF*(5d0/6d0*ALP1**2
     .     +3d0/2d0*ALP2**2)*F2*N5*PM(1)**2
      Y(27)=COEF*(10d0/3d0*ALP1**2)*F2*N5*PM(1)**2
      Y(28)=COEF*(5d0/6d0*ALP1**2
     .     +3d0/2d0*ALP2**2)*F2*N5*PM(1)**2
      Y(29)=COEF*(10d0/3d0*ALP1**2)*F2*N5*PM(1)**2
      Y(30)=PM(5)
      Y(31)=PM(6)
      Y(32)=PM(7)
      Y(33)=PM(8)
      Y(34)=0d0
      Y(35)=0d0

      X1=COEF*DLOG(PM(15)/Q2)
      X2=0d0

      CALL ODEINTS(Y,NN,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)

      IF(IFAIL.NE.0)THEN
       RETURN
      ENDIF

      PS(1)=Y(17)
      PS(2)=Y(18)
      PS(3)=Y(19)
      PS(4)=Y(12)
      PS(5)=Y(13)
      PS(6)=Y(30)
      PS(7)=Y(31)
      PS(8)=Y(32)
      PS(9)=Y(33)
      PS(10)=Y(34)
      PS(11)=Y(4)
      PS(12)=Y(5)
      PS(13)=Y(6)
      PS(14)=(Y(1)+Y(2))/2d0

      END
