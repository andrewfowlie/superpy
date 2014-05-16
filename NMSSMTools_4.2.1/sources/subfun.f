      SUBROUTINE DIAGN(N,A,D,V,EPS)

*   Diagonalization of a real symmetric NxN matrix
*   Using the Jacobi method (from Numerical Recipes)

      IMPLICIT NONE

      INTEGER N,I,J,IP,IQ,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION A(N,N),D(N),V(N,N),B(N),Z(N)
      DOUBLE PRECISION EPS,SM,THR,G,H,C,S,T,THETA,TAU

      DO IP=1,N
       DO IQ=1,N
        V(IP,IQ)=0d0
       ENDDO
       V(IP,IP)=1d0
       B(IP)=A(IP,IP)
       D(IP)=B(IP)
       Z(IP)=0d0
      ENDDO

      DO I=1,NMAX

       SM=0d0
       DO IP=1,N-1
        DO IQ=IP+1,N
         SM=SM+DABS(A(IP,IQ))
        ENDDO
       ENDDO
       IF(SM.LE.EPS)RETURN

       IF(I.LT.4)THEN
        THR=.2d0*SM/N**2
       ELSE
        THR=0d0
       ENDIF

       DO IP=1,N-1
       DO IQ=IP+1,N
        IF(DABS(A(IP,IQ)).LT.EPS*MIN(DABS(D(IP)),DABS(D(IQ)))
     .     .AND. I.GT.4)THEN
         A(IP,IQ)=0d0
        ELSEIF(DABS(A(IP,IQ)).GT.THR)THEN
         H=D(IQ)-D(IP)
         IF(DABS(A(IP,IQ)).LT.EPS*DABS(H))THEN
          T=A(IP,IQ)/H
         ELSE
          THETA=.5d0*H/A(IP,IQ)
          T=1d0/(DABS(THETA)+DSQRT(1d0+THETA**2))
          IF(THETA.LT.0d0)T=-T
         ENDIF
         C=1d0/DSQRT(1d0+T**2)
         S=T*C
         TAU=S/(1d0+C)
         H=T*A(IP,IQ)
         Z(IP)=Z(IP)-H
         Z(IQ)=Z(IQ)+H
         D(IP)=D(IP)-H
         D(IQ)=D(IQ)+H
         A(IP,IQ)=0d0
         DO J=1,IP-1
          G=A(J,IP)
          H=A(J,IQ)
          A(J,IP)=G-S*(H+G*TAU)
          A(J,IQ)=H+S*(G-H*TAU)
         ENDDO
         DO J=IP+1,IQ-1
          G=A(IP,J)
          H=A(J,IQ)
          A(IP,J)=G-S*(H+G*TAU)
          A(J,IQ)=H+S*(G-H*TAU)
         ENDDO
         DO J=IQ+1,N
          G=A(IP,J)
          H=A(IQ,J)
          A(IP,J)=G-S*(H+G*TAU)
          A(IQ,J)=H+S*(G-H*TAU)
         ENDDO
         DO J=1,N
          G=V(J,IP)
          H=V(J,IQ)
          V(J,IP)=G-S*(H+G*TAU)
          V(J,IQ)=H+S*(G-H*TAU)
         ENDDO
        ENDIF
       ENDDO
       ENDDO

       DO IP=1,N
        B(IP)=B(IP)+Z(IP)
        D(IP)=B(IP)
        Z(IP)=0d0
       ENDDO

      ENDDO

      RETURN
      END


      SUBROUTINE SORTN(N,D,V)

*   Reordering of the eigenvalues D(I), I=1..N
*   and corresponding eigenvectors V(J,I), J=1..N in ascending order
*   (from Numerical Recipes)

      IMPLICIT NONE

      INTEGER N,I,J,K
      DOUBLE PRECISION D(N),V(N,N),P

      DO I=1,N-1
       K=I
       P=D(I)
       DO J=I+1,N
        IF(D(J).LT.P)THEN
         K=J
         P=D(J)
        ENDIF
       ENDDO
       IF(K.NE.I)THEN
        D(K)=D(I)
        D(I)=P
        DO J=1,N
         P=V(J,I)
         V(J,I)=V(J,K)
         V(J,K)=P
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE SORTNA(N,D,V)

*   Reordering of the absolute value of the eigenvalues D(I), I=1..N
*   and corresponding eigenvectors V(J,I), J=1..N in ascending order
*   (from Numerical Recipes)

      IMPLICIT NONE

      INTEGER N,I,J,K
      DOUBLE PRECISION D(N),V(N,N),P

      DO I=1,N-1
       K=I
       P=D(I)
       DO J=I+1,N
        IF(DABS(D(J)).LT.DABS(P))THEN
         K=J
         P=D(J)
        ENDIF
       ENDDO
       IF(K.NE.I)THEN
        D(K)=D(I)
        D(I)=P
        DO J=1,N
         P=V(J,I)
         V(J,I)=V(J,K)
         V(J,K)=P
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END


      DOUBLE PRECISION FUNCTION RUNM(Q,NF)

*   Subroutine to calculate the quark running masses

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594d0)

      DIMENSION AM(NN),YMSB(NN)

      COMMON/ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,gg1,gg2,S2TW
      COMMON/SMSPEC/AMS,AMC,AMB,AMBP,AMT,AMTAU,AMMUON,AMZ,AMW

      B0(NF)= (33d0-2d0*NF)/12d0
      B1(NF)= (102d0-38d0/3d0*NF)/16d0
      B2(NF)= (2857d0/2d0-5033d0/18d0*NF+325d0/54d0*NF**2)/64d0
      G0(NF)= 1d0
      G1(NF)= (202d0/3d0-20d0/9d0*NF)/16d0
      G2(NF)= (1249d0-(2216d0/27d0+160d0/3d0*ZETA3)*NF
     .      - 140d0/81d0*NF**2)/64d0
      C1(NF)= G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF)= ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .      + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .      - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2d0
      TRAN(X,XK)= 1d0+4d0/3d0*ALPHAS(X,2)/PI+XK*(ALPHAS(X,2)/PI)**2
      CQ(X,NF)= (2d0*B0(NF)*X)**(G0(NF)/B0(NF))
     .  * (1d0+C1(NF)*X+C2(NF)*X**2)

      PI= 4d0*DATAN(1d0)
      ACC= 1.D-8
      AM(1)= 0
      AM(2)= 0
      AM(3)= AMS
      AM(4)= AMC
      AM(5)= AMBP
      AM(6)= AMT
      XK= 16.11d0
      DO 1 I=1,NF-1
       XK= XK - 1.04d0*(1d0-AM(I)/AM(NF))
1     CONTINUE
      IF(NF.GE.4)THEN
       XMSB= AM(NF)/TRAN(AM(NF),0d0)
       XMHAT= XMSB/CQ(ALPHAS(AM(NF),2)/PI,NF)
      ELSE
       XMSB= 0
       XMHAT= 0
      ENDIF
      YMSB(3)= AMS
      IF(NF.EQ.3)THEN
       YMSB(4)= YMSB(3)*CQ(ALPHAS(AM(4),2)/PI,3)/
     .            CQ(ALPHAS(1d0,2)/PI,3)
       YMSB(5)= YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .            CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .            CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4)= XMSB
       YMSB(5)= YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .            CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .            CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5)= XMSB
       YMSB(4)= YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .            CQ(ALPHAS(AM(5),2)/PI,4)
       YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .            CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6)= XMSB
       YMSB(5)= YMSB(6)*CQ(ALPHAS(AM(5),2)/PI,5)/
     .            CQ(ALPHAS(AM(6),2)/PI,5)
       YMSB(4)= YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .            CQ(ALPHAS(AM(5),2)/PI,4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0= 3
       Q0= 1d0
      ELSEIF(Q.LE.AMBP)THEN
       N0= 4
       Q0= AMC
      ELSEIF(Q.LE.AMT)THEN
       N0= 5
       Q0= AMBP
      ELSE
       N0= 6
       Q0= AMT
      ENDIF
      IF(NF.GT.3)THEN
       XKFAC= TRAN(AM(NF),0d0)/TRAN(AM(NF),XK)
      ELSE
       XKFAC= 1d0
      ENDIF
      RUNM= YMSB(N0)*CQ(ALPHAS(Q,2)/PI,N0)/
     .         CQ(ALPHAS(Q0,2)/PI,N0)
     .       * XKFAC

      RETURN
      END


*  Running alpha_s and aux. subroutines/functions as in HDECAY

      DOUBLE PRECISION FUNCTION ALPHAS(Q,N)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION XLB(6)

      COMMON/ALSLAM/XLB1(6),XLB2(6)
      COMMON/ALS/XLAMBDA,AMC,AMBP,AMT,N0

      B0(NF)=33d0-2d0*NF
      B1(NF)=6d0*(153d0-19d0*NF)/B0(NF)**2
      ALS1(NF,X)=12d0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12d0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .    *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .     /DLOG(X**2/XLB(NF)**2))


      PI=4d0*DATAN(1d0)

      IF(N.EQ.1)THEN
       DO 1 I=1,6
        XLB(I)=XLB1(I)
1     CONTINUE
      ELSE
       DO 2 I=1,6
        XLB(I)=XLB2(I)
2     CONTINUE
      ENDIF
      IF(Q.LT.AMC)THEN
       NF=3
      ELSEIF(Q.LE.AMBP)THEN
       NF=4
      ELSEIF(Q.LE.AMT)THEN
       NF=5
      ELSE
       NF=6
      ENDIF
      IF(N.EQ.1)THEN
        ALPHAS=ALS1(NF,Q)
      ELSE
        ALPHAS=ALS2(NF,Q)
      ENDIF

      RETURN
      END


      SUBROUTINE ALSINI(ACC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION XLB(6)

      COMMON/ALSLAM/XLB1(6),XLB2(6)
      COMMON/ALS/XLAMBDA,AMC,AMBP,AMT,N0

      PI=4d0*DATAN(1d0)
      XLB1(1)=0d0
      XLB1(2)=0d0
      XLB2(1)=0d0
      XLB2(2)=0d0
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2d0/25d0)
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2d0/23d0)
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
      ENDIF
      DO 1 I=1,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2d0/25d0)
     .       *(2d0*DLOG(AMC/XLB(3)))**(-107d0/1875d0)
       XLB(4)=XITER(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
     .       *(2d0*DLOG(AMBP/XLB(4)))**(-963d0/13225d0)
       XLB(5)=XITER(AMBP,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .      *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
     .       *(2d0*DLOG(AMBP/XLB(4)))**(-963d0/13225d0)
       XLB(5)=XITER(AMBP,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .       *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .      *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
     .       *(2d0*DLOG(AMBP/XLB(5)))**(963d0/14375d0)
       XLB(4)=XITER(AMBP,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .       *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .      *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2d0/23d0)
     .      *(2d0*DLOG(AMT/XLB(6)))**(321d0/3703d0)
       XLB(5)=XITER(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
     .       *(2d0*DLOG(AMBP/XLB(5)))**(963d0/14375d0)
       XLB(4)=XITER(AMBP,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .       *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 2 I=1,6
       XLB2(I)=XLB(I)
2     CONTINUE

      RETURN
      END


      DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      B0(NF)=33d0-2d0*NF
      B1(NF)=6d0*(153d0-19d0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12d0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .        *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .        /DLOG(X**2/XLB**2))
      AA(NF)=12d0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2d0*(1d0+DSQRT(1d0-4d0*B*DLOG(X)))
      PI=4d0*DATAN(1d0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2d0)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
      XITER=XLB2

      RETURN
      END


      DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)

*  Iteration routine to determine improved Lambda's

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,gg1,gg2,S2TW
      COMMON/SMSPEC/AMS,AMC,AMB,AMBP,AMT,AMTAU,AMMUON,AMZ,AMW

      B0(NF)=33d0-2d0*NF
      B1(NF)=6d0*(153d0-19d0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12d0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .        *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .        /DLOG(X**2/XLB**2))
      AA(NF)=12d0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2d0*(1d0+DSQRT(1d0-4d0*B*DLOG(X)))
      PI=4d0*DATAN(1d0)
      NF=5
      Q=AMZ
      XLB=Q*DEXP(-AA(NF)/ALP/2d0)
      IF(NO.EQ.1)GOTO 111
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB**2)
      A=AA(NF)/ALP
      B=BB(NF)*ALP
      XX=XIT(A,B,X)
      XLB=Q*DEXP(-XX/2d0)
      Y1=ALP
      Y2=ALS2(NF,Q,XLB)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
111   XITLA=XLB

      RETURN
      END


      DOUBLE PRECISION FUNCTION FINT(Z,XX,YY)

*  One-dimensional cubic interpolation
*  Z  = wanted point
*  XX = array of 4 discrete x-values around Z
*  YY = array of 4 discrete function-values around Z

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION XX(4),YY(4)

      X = DLOG(Z)
      X0=DLOG(XX(1))
      X1=DLOG(XX(2))
      X2=DLOG(XX(3))
      X3=DLOG(XX(4))
      Y0=DLOG(YY(1))
      Y1=DLOG(YY(2))
      Y2=DLOG(YY(3))
      Y3=DLOG(YY(4))
      A0=(X-X1)*(X-X2)*(X-X3)/(X0-X1)/(X0-X2)/(X0-X3)
      A1=(X-X0)*(X-X2)*(X-X3)/(X1-X0)/(X1-X2)/(X1-X3)
      A2=(X-X0)*(X-X1)*(X-X3)/(X2-X0)/(X2-X1)/(X2-X3)
      A3=(X-X0)*(X-X1)*(X-X2)/(X3-X0)/(X3-X1)/(X3-X2)
      FINT=DEXP(A0*Y0+A1*Y1+A2*Y2+A3*Y3)

      RETURN
      END


*   Spence function and auxiliary functions as in HDECAY

      DOUBLE PRECISION FUNCTION SP(X)

*  REAL dilogarithm (Spence-function)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX CX,LI2

      CX = DCMPLX(X,0d0)
      SP = DREAL(LI2(CX))

      RETURN
      END


      DOUBLE COMPLEX FUNCTION LI2(X)

*  COMPLEX dilogarithm (Spence-function)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX X,Y,CLI2

      COMMON/CONST/ZETA2,ZETA3

      ZERO=1.D-16
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      LI2=0
      IF(R2.LE.ZERO)THEN
        LI2=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1d0.AND.XI.EQ.0d0)THEN
        IF(XR.EQ.1d0)THEN
          LI2=DCMPLX(ZETA2)
        ELSE
          LI2=-DCMPLX(ZETA2/2d0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1d0.AND.RR.GT.0.5d0)THEN
        Y=(X-1d0)/X
        LI2=CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1d0-X)+0.5d0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1d0.AND.RR.LE.0.5d0)THEN
        Y=1d0/X
        LI2=-CLI2(Y)-ZETA2-0.5d0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1d0.AND.XR.GT.0.5d0)THEN
        Y=1d0-X
        LI2=-CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1d0-X)
       RETURN
      ELSEIF(R2.LE.1d0.AND.XR.LE.0.5d0)THEN
        Y=X
        LI2=CLI2(Y)
        RETURN

      ENDIF
      END


      DOUBLE COMPLEX FUNCTION CLI2(X)

*  Taylor-expansion for complex dilogarithm (Spence-function)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX X,Z

      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/POLY/NBER

      N=NBER-1
      Z=-CDLOG(1d0-X)
      CLI2=B2(NBER)
      DO 111 I=N,1,-1
        CLI2=Z*CLI2+B2(I)
111   CONTINUE
      CLI2=Z**2*CLI2+Z

      RETURN
      END


      DOUBLE PRECISION FUNCTION FACULT(N)

*  DOUBLE PRECISION version of FACULTY

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      FACULT=1d0
      IF(N.EQ.0)RETURN
      DO 999 I=1,N
        FACULT=FACULT*DFLOAT(I)
999   CONTINUE

      RETURN
      END


      SUBROUTINE BERNINI(N)

*  Initialization of coefficients for polylogarithms

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(18),PB(19)

      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/CONST/ZETA2,ZETA3
      COMMON/POLY/NBER

      NBER=N
      PI=4d0*DATAN(1d0)

      B(1)=-1d0/2d0
      B(2)=1d0/6d0
      B(3)=0d0
      B(4)=-1d0/30d0
      B(5)=0d0
      B(6)=1d0/42d0
      B(7)=0d0
      B(8)=-1d0/30d0
      B(9)=0d0
      B(10)=5d0/66d0
      B(11)=0d0
      B(12)=-691d0/2730d0
      B(13)=0d0
      B(14)=7d0/6d0
      B(15)=0d0
      B(16)=-3617d0/510d0
      B(17)=0d0
      B(18)=43867d0/798d0
      ZETA2=PI**2/6d0
      ZETA3=1.202056903159594d0

      DO 995 I=1,18
        B2(I)=B(I)/FACULT(I+1)
        B12(I)=DFLOAT(I+1)/FACULT(I+2)*B(I)/2d0
        PB(I+1)=B(I)
        B3(I)=0d0
995   CONTINUE
      PB(1)=1d0
      DO 996 I=1,18
      DO 996 J=0,I
       B3(I)=B3(I)+PB(J+1)*PB(I-J+1)/FACULT(I-J)/FACULT(J+1)
     .       /DFLOAT(I+1)
996   CONTINUE

      RETURN
      END


*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*   Passarino-Veltman one- and two-points functions A0, B0 and B1
*   orig from LoopTools, http://www.feynarts.de/looptools/
*   taken from Suspect2.3, modified by S. Kraml, 7 March 2005
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function NMA0(m2,q)

      implicit none
      double precision m2,q
      if(m2.ne.0d0) then
       NMA0 = m2 * (1d0-dlog( m2/q ))
      else
       NMA0 = 0d0
      endif
      end

      double precision function NMB0(p,m1,m2,q)

*     note: all input is quadratical: p=p^2, m1=m1^2, m2=m2^2, q=q^2

      implicit none
      double precision p, m1, m2
      double precision mudim2, divergence, lambda2, q
      double precision acc, eps, minacc
      double complex x1, x2, y1, y2, r, be0
      double complex Ieps, onePeps, oneMeps
      COMMON/cutoff/mudim2, divergence, lambda2
      parameter (acc = 1.D-12)
      parameter (eps = 1.D-20)
      parameter (Ieps = (0d0,1d0)*eps)
      parameter (onePeps = 1d0 + Ieps)
      parameter (oneMeps = 1d0 - Ieps)

      double complex fpv, xlogx
      external fpv, xlogx

      divergence = 0d0
      lambda2 = 0d0
      mudim2 = q
      minacc = acc*(m1 + m2)

* general case
      if(abs(p) .gt. minacc) then
  
      CALL roots(p, m1, m2, x1, x2, y1, y2, r)
        if(abs(y1) .gt. .5d0 .and. abs(y2) .gt. .5d0) then
          be0 = -log(m2/mudim2) -
     +      fpv(1, x1, y1) - fpv(1, x2, y2)
        else if(abs(x1) .lt. 10d0 .and. abs(x2) .lt. 10d0) then
          be0 = 2 - log(p*oneMeps/mudim2) +
     +      xlogx(-x1) + xlogx(-x2) - xlogx(y1) - xlogx(y2)
        else if(abs(x1) .gt. .5d0 .and. abs(x2) .gt. .5d0) then
          be0 = -log(m1/mudim2) -
     +      fpv(1, y1, x1) - fpv(1, y2, x2)
        else
          be0 = 1.D100
        endif

* zero momentum
      else if(abs(m1 - m2) .gt. minacc) then
        x2 = oneMeps*m1/(m1 - m2)
        y2 = oneMeps*m2/(m2 - m1)
        if(abs(y2) .gt. .5d0) then
          be0 = -log(m2/mudim2) - fpv(1, x2, y2)
        else
          be0 = -log(m1/mudim2) - fpv(1, y2, x2)
        endif
      else
        be0 = -log(m2/mudim2)
      endif

      NMB0 = dble(be0 + divergence)

      end


      double precision function NMB1(s,mi,mj,q)

* note: all input is quadratical: s=p^2, mi=m1^2, mj=m2^2, q=q^2

      implicit none
      double precision s,mi,mj,NMB0,NMA0,q

      if(mi.eq.mj) then
       NMB1 = NMB0(s,mi,mj,q)/2d0
      else
       NMB1= (NMA0(mj,q) - NMA0(mi,q)
     .   + (s+mi-mj)*NMB0(s,mi,mj,q))/(2d0*s)
      endif
      end


*---------------------------------------------------------------------
* auxiliary functions used by the B0,B1 two-point functions
* from Looptools http://www.feynarts.de/looptools/
*---------------------------------------------------------------------

      subroutine roots(p, m1, m2, x1, x2, y1, y2, r)

      implicit none
      double precision p, m1, m2
      double complex x1, x2, y1, y2, r
      double precision mudim2, divergence, lambda2
      COMMON/cutoff/mudim2, divergence, lambda2
      double precision acc, eps
      double complex Ieps, onePeps, oneMeps
      parameter (acc = 1D-12)
      parameter (eps = 1D-20)
      parameter (Ieps = (0d0,1d0)*eps)
      parameter (onePeps = 1d0 + Ieps)
      parameter (oneMeps = 1d0 - Ieps)
      double precision q

      r = sqrt(dcmplx(p*(p - 2*(m1 + m2)) + (m1 - m2)**2))
      q = p + m1 - m2
      x1 = (q + r)/2d0/p
      x2 = (q - r)/2d0/p
      if(abs(x2) .gt. abs(x1)) then
        x1 = m1/p/x2
      else if(abs(x1) .gt. abs(x2)) then
        x2 = m1/p/x1
      endif
      x1 = x1 + abs(p*x1)/p*Ieps
      x2 = x2 - abs(p*x2)/p*Ieps
      q = p - m1 + m2
      y2 = (q + r)/2d0/p
      y1 = (q - r)/2d0/p
      if(abs(y2) .gt. abs(y1)) then
        y1 = m2/p/y2
      else if(abs(y1) .gt. abs(y2)) then
        y2 = m2/p/y1
      endif
      y1 = y1 - abs(p*y1)/p*Ieps
      y2 = y2 + abs(p*y2)/p*Ieps
      end


      double complex function fpv(n, x, y)

      implicit none
      integer n
      double complex x, y
      double precision mudim2, divergence, lambda2
      COMMON/cutoff/mudim2, divergence, lambda2
      double precision acc, eps
      double complex Ieps, onePeps, oneMeps
      parameter (acc = 1D-12)
      parameter (eps = 1D-20)
      parameter (Ieps = (0,1)*eps)
      parameter (onePeps = 1 + Ieps)
      parameter (oneMeps = 1 - Ieps)
      integer m
      double complex xm
      if(abs(x) .lt. 10d0) then
        if(n .eq. 0) then
          fpv = -log(-y/x)
        else if(abs(x) .lt. acc) then
          fpv = -1d0/n
        else
          fpv = 0
          xm = 1
          do m = 0, n - 1
            fpv = fpv - xm/(n - m)
            xm = xm*x
          enddo
          fpv = fpv - xm*log(-y/x)
        endif
      else
        fpv = 0
        xm = 1
        do m = 1, 30
          xm = xm/x
          fpv = fpv + xm/(m + n)
          if(abs(xm/fpv) .lt. acc**2) return
        enddo
      endif
      end


      double complex function yfpv(n, x, y)

      implicit none
      integer n
      double complex x, y
      double complex fpv
      external fpv
      if(abs(y) .eq. 0d0) then
        yfpv = 0
      else
        yfpv = y*fpv(n, x, y)
      endif
      end


      double complex function xlogx(x)

      implicit none
      double complex x
      if(abs(x) .eq. 0d0) then
        xlogx = 0
      else
        xlogx = x*log(x)
      endif
      end


      DOUBLE PRECISION FUNCTION RUNMB(Q)
      
*   Subroutine to calculate the running b quark mass for Q > MB

      IMPLICIT NONE
      DOUBLE PRECISION Q
      DOUBLE PRECISION PI,ALPHAS,ALMB,ALMT,ALQ,U5MTMB,U6QMT
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      PI=4d0*DATAN(1d0)

      ALQ=ALPHAS(Q,2)
      ALMB=ALPHAS(MB,2)

      IF(Q.LE.MT) THEN

       RUNMB=MB*(ALQ/ALMB)**(12d0/23d0)*(1d0+7462d0*(ALQ-ALMB)/
     .      (4d0*PI*1587d0))

      ELSE

       ALMT=ALPHAS(MT,2)
       U5MTMB=(ALMT/ALMB)**(12d0/23d0)*(1d0+7462d0*(ALMT-ALMB)/
     .      (4d0*PI*1587d0))
        U6QMT=(ALQ/ALMT)**(4d0/7d0)*(1d0+7398d0*(ALQ-ALMT)/
     .       (4d0*PI*1323d0))
       RUNMB=MB*U6QMT*U5MTMB

      ENDIF

      END


      DOUBLE PRECISION FUNCTION RAN2(IDUM)

      IMPLICIT NONE

      INTEGER IDUM
      DOUBLE PRECISION DA,DB,DC
      PARAMETER(DA=16807.d0,DB=2147483647.d0,DC=2147483648.d0)

      IDUM=INT(ABS(MOD(DA*IDUM,DB)+0.5d0))
      RAN2=DFLOAT(IDUM)/DC

      END


      DOUBLE PRECISION FUNCTION GAU(IDUM)

      IMPLICIT NONE

      INTEGER IDUM,ISET
      DOUBLE PRECISION F,SET,R,V1,V2,RAN2
      SAVE ISET,SET
      DATA ISET/0/
      IF(ISET.EQ.0)THEN
 5      V1=2d0*RAN2(IDUM)-1d0
        V2=2d0*RAN2(IDUM)-1d0
        R=V1**2+V2**2
        IF(R.GE.1d0.OR.R.EQ.0d0)GOTO 5
        F=DSQRT(-2d0*DLOG(R)/R)
        SET=V1*F
        GAU=V2*F
        ISET=1
      ELSE
        GAU=SET
        ISET=0
      ENDIF
      END


