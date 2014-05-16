*   Subroutines to integrate the RGEs for the soft terms


      SUBROUTINE ODEINTS(YSTART,NVAR,X1,X2,EPS,DERIVSS,RKQSS,IFAIL)

*   Driver subroutine to integrate Ordinary Differential Equations
*   using the Runge-Kutta 5th order method with adaptative step size
*   (from Numerical Recipes)
*   IFAIL=1 stepsize smaller than minimum in ODEINTS
*   IFAIL=2 too many steps in ODEINTS
*   IFAIL=3 stepsize underflow in RKQSS
*   IFAIL=4 max(|dy(i)/dx|*(1+x)/(1+|y(i)|)) > 1/eps^2
*   NOTE: It is assumed that YSTART(5)=Kappa, NOT Kappa**2

      IMPLICIT NONE

      INTEGER NVAR,NMAX,IFAIL,MAXSTP,NSTP,I
      PARAMETER(MAXSTP=10000,NMAX=500)

      DOUBLE PRECISION YSTART(NVAR),Y(NMAX),YSCAL(NMAX)
      DOUBLE PRECISION DYDX(NMAX),X,X1,X2,EPS,TINY
      DOUBLE PRECISION H1,HMIN,H,HDID,HNEXT

      EXTERNAL DERIVSS,RKQSS

      IFAIL=0
      TINY=EPS**4
      H1=DSQRT(EPS)*DABS(X2-X1)
      HMIN=EPS**2*DABS(X2-X1)

      X=X1
      H=DSIGN(H1,X2-X1)
      DO I=1,NVAR
       Y(I)=YSTART(I)
      ENDDO

      DO NSTP=1,MAXSTP

       CALL DERIVSS(NVAR,X,Y,DYDX)
      
       DO I=1,NVAR
        IF(DABS(DYDX(I))*(1d0+X)/(1d0+DABS(Y(I))).GT.EPS**(-2))IFAIL=4
       ENDDO
       IF(IFAIL.NE.0)RETURN

       DO I=1,NVAR
        YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
       ENDDO

       IF((X+H-X2)*(X+H-X1).GT.0d0) H=X2-X

       CALL
     .  RKQSS(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVSS,IFAIL)

       IF(IFAIL.NE.0)RETURN

       IF((X-X2)*(X2-X1).GE.0d0)THEN
        DO I=1,NVAR
         YSTART(I)=Y(I)
        ENDDO

        RETURN
       ENDIF

       IF(ABS(HNEXT).LT.HMIN)THEN
        IFAIL=1
        RETURN
       ENDIF

       H=HNEXT

      ENDDO

      IFAIL=2

      RETURN
      END


      SUBROUTINE
     .       RKQSS(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVSS,IFAIL)

*   Stepper subroutine for ODEINTS

      IMPLICIT NONE

      INTEGER IFAIL,I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,DYDX(N),Y(N),YSCAL(N)
      DOUBLE PRECISION ERRMAX,H,HTEMP,XNEW,YERR(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION SAFETY,PGROW,PSHRNK,ERRCON

      EXTERNAL DERIVSS

      SAFETY=.9d0
      PGROW=-.2d0
      PSHRNK=-.25d0
      ERRCON=(5d0/SAFETY)**(1d0/PGROW)
      H=HTRY

1     CALL RKCKS(Y,DYDX,N,X,H,YTEMP,YERR,DERIVSS)

      ERRMAX=0d0
      DO I=1,N
       ERRMAX=MAX(ERRMAX,DABS(YERR(I)/YSCAL(I)))
      ENDDO
      ERRMAX=ERRMAX/EPS

      IF(ERRMAX.GT.1d0)THEN

       HTEMP=SAFETY*H*(ERRMAX**PSHRNK)
       H=SIGN(MAX(DABS(HTEMP),.1d0*DABS(H)),H)
       XNEW=X+H
       IF(XNEW.EQ.X)THEN
        IFAIL=3
        RETURN
       ENDIF
       GOTO 1

      ELSE

       IF(ERRMAX.GT.ERRCON)THEN
        HNEXT=SAFETY*H*(ERRMAX**PGROW)
       ELSE
        HNEXT=5d0*H
       ENDIF
       HDID=H
       X=X+H
       DO I=1,N
        Y(I)=YTEMP(I)
       ENDDO
       RETURN

      ENDIF

      END


      SUBROUTINE RKCKS(Y,DYDX,N,X,H,YOUT,YERR,DERIVSS)

*   Algorithm subroutine for ODEINT

      IMPLICIT NONE

      INTEGER I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION H,X,DYDX(N),Y(N),YERR(N),YOUT(N)
      DOUBLE PRECISION AK2(NMAX),AK3(NMAX),AK4(NMAX)
      DOUBLE PRECISION AK5(NMAX),AK6(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION A2,A3,A4,A5,A6,C1,C3,C4,C6
      DOUBLE PRECISION B21,B31,B32,B41,B42,B43,B51,B52,B53,B54
      DOUBLE PRECISION B61,B62,B63,B64,B65,DC1,DC3,DC4,DC5,DC6

      PARAMETER(A2=.2d0, A3=.3d0, A4=.6d0, A5=1d0, A6=.875d0,
     .       B21=.2d0, B31=3d0/40d0, B32=9d0/40d0, B41=.3d0, B42=-.9d0,
     .       B43=1.2d0,B51=-11d0/54d0, B52=2.5d0, B53=-70d0/27d0,
     .       B54=35d0/27d0, B61=1631d0/55296d0, B62=175d0/512d0,
     .       B63=575d0/13824d0, B64=44275d0/110592d0,
     .       B65=253d0/4096d0, C1=37d0/378d0, C3=250d0/621d0,
     .       C4=125d0/594d0, C6=512d0/1771d0, DC1=C1-2825d0/27648d0,
     .       DC3=C3-18575d0/48384d0, DC4=C4-13525d0/55296d0,
     .       DC5=-277d0/14336d0, DC6=C6-.25d0)

      EXTERNAL DERIVSS

      DO I=1,N
       YTEMP(I)=Y(I)+B21*H*DYDX(I)
      ENDDO
      CALL DERIVSS(N,X+A2*H,YTEMP,AK2)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B31*DYDX(I)+B32*AK2(I))
      ENDDO
      CALL DERIVSS(N,X+A3*H,YTEMP,AK3)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B41*DYDX(I)+B42*AK2(I)+B43*AK3(I))
      ENDDO
      CALL DERIVSS(N,X+A4*H,YTEMP,AK4)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B51*DYDX(I)+B52*AK2(I)+B53*AK3(I)+B54*AK4(I))
      ENDDO
      CALL DERIVSS(N,X+A5*H,YTEMP,AK5)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B61*DYDX(I)+B62*AK2(I)+B63*AK3(I)+B64*AK4(I)+
     .  B65*AK5(I))
      ENDDO
      CALL DERIVSS(N,X+A6*H,YTEMP,AK6)
      DO I=1,N
       YOUT(I)=Y(I)+H*(C1*DYDX(I)+C3*AK3(I)+C4*AK4(I)+C6*AK6(I))
      ENDDO
      DO I=1,N
       YERR(I)=H*(DC1*DYDX(I)+DC3*AK3(I)+DC4*AK4(I)+DC5*AK5(I)+DC6*
     .  AK6(I))
      ENDDO

      RETURN
      END


      SUBROUTINE DERIVSS(N,X,Y,F)

*   2-loop Renormalization group equations for g1, g2, g3,
*   lambda, kappa, htop, hbot, htau and for all soft terms
*   to be integrated by ODEINTS

      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION X,Y(N),F(N),PI,c2,S
      DOUBLE PRECISION G1,G2,G3,L2,K2,HT2,HB2,HTAU2
      DOUBLE PRECISION M1,M2,M3,AL,AK,AT,AB,ATAU,AMUON
      DOUBLE PRECISION MH1,MH2,MS,MQ3,MU3,MD3,MQ,MU,MD
      DOUBLE PRECISION ML3,ME3,ML,ME,SP,SIG1,SIG2,SIG3
      DOUBLE PRECISION TMQ,TMU,TMD,TML,TME
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H,LL,KK

      PI=4d0*DATAN(1d0)
      c2=1d0/(16d0*PI**2)
      
      G1=Y(1)
      G2=Y(2)
      G3=Y(3)
      L2=Y(4)
      LL=DSQRT(L2)
* NOTE: Y(5)=K, NOT K**2
      K2=Y(5)**2
      KK=Y(5)
      HT2=Y(6)
      HB2=Y(7)
      HTAU2=Y(8)
      M1=Y(9)
      M2=Y(10)
      M3=Y(11)
      AL=Y(12)
      AK=Y(13)
      AT=Y(14)
      AB=Y(15)
      ATAU=Y(16)
      MH1=Y(17)
      MH2=Y(18)
      MS=Y(19)
      MQ3=Y(20)
      MU3=Y(21)
      MD3=Y(22)
      MQ=Y(23)
      MU=Y(24)
      MD=Y(25)
      ML3=Y(26)
      ME3=Y(27)
      ML=Y(28)
      ME=Y(29)
      XIF=Y(30)
      XIS=Y(31)
      MUP=Y(32)
      MSP=Y(33)
      M3H=Y(34)
      AMUON=Y(35)

      TMQ=MQ3+2d0*MQ
      TMU=MU3+2d0*MU
      TMD=MD3+2d0*MD
      TML=ML3+2d0*ML
      TME=ME3+2d0*ME

      S= g1*(MH1 - MH2 + TMQ - 2d0*TMU + TMD + TME - TML)

      SP= HT2*(-3d0*MH1 - MQ3 + 4d0*MU3)
     .  + HB2*(3d0*MH2 - MQ3 - 2d0*MD3)
     .  + HTAU2*(MH2 + ML3 - 2d0*ME3) + L2*(MH2 - MH1)
     .  + (G1/18d0 + 3d0/2d0*G2 + 8d0/3d0*G3)*TMQ
     .  - (16d0/9d0*G1 + 16d0/3d0*G3)*TMU
     .  + (2d0/9d0*G1 + 8d0/3d0*G3)*TMD
     .  + (G1/2d0 + 3d0/2d0*G2)*(MH1-MH2-TML) + 2d0*G1*TME

      SIG1= G1*(MH1 + MH2 + TMQ/3d0 + 8d0/3d0*TMU + 2d0/3d0*TMD
     .    + TML + 2d0*TME)

      SIG2=G2*(MH1 + MH2 + 3d0*TMQ + TML)

      SIG3=G3*(2d0*TMQ + TMU + TMD)

      F(1)= 11d0*g1**2
     .    + c2*g1**2*(199d0/9d0*g1 + 9d0*g2 + 88d0/3d0*g3
     .    - 2d0*L2 - 26d0/3d0*HT2 - 14d0/3d0*HB2 - 6d0*HTAU2)

      F(2)= g2**2
     .    + c2*g2**2*(3d0*g1 + 25d0*g2 + 24d0*g3
     .    - 2d0*L2 - 6d0*HT2 - 6d0*HB2 - 2d0*HTAU2)

      F(3)= -3d0*g3**2
     .    + c2*g3**2*(11d0/3d0*g1 + 9d0*g2 + 14d0*g3
     .    - 4 d0*HT2 - 4d0*HB2)

      F(4)= L2*(-g1 - 3d0*g2
     .    + 4d0*L2 + 2d0*K2 + 3d0*HT2 + 3d0*HB2 + HTAU2)
     .    + c2*L2*(L2*(2d0*g1 + 6d0*g2
     .    - 10d0*L2 - 12d0*K2 - 9d0*HT2 - 9d0*HB2 - 3d0*HTAU2)
     .    + 23d0/2d0*g1**2 + 15d0/2d0*g2**2 + 3d0*g1*g2
     .    + 4d0/3d0*g1*HT2 - 2d0/3d0*g1*HB2 + 2d0*g1*HTAU2
     .    + 16d0*g3*HT2 + 16d0*g3*HB2 - 3d0*HTAU2**2
     .    - 9d0*HT2**2 - 9d0*HB2**2 - 6d0*HT2*HB2 - 8d0*K2**2)

* NOTE: KK=Kappa, K2=Kappa**2

      F(5)= 3d0*KK*(L2 + K2)
     .    + c2*3d0*KK*(-4d0*K2**2 + L2*(g1 + 3d0*g2
     .    - 2d0*L2 - 4d0*K2 - 3d0*HT2 - 3d0*HB2 - HTAU2))

      F(6)= HT2*(-13d0/9d0*g1 - 3d0*g2 - 16d0/3d0*g3
     .    + L2 + 6d0*HT2 + HB2)
     .    + c2*HT2*(HT2*(2d0*g1 + 6d0*g2 + 16d0*g3
     .    - 3d0*L2 - 22d0*HT2 - 5d0*HB2)
     .    + 2d0/3d0*g1*HB2 - 3d0*L2**2 - 2d0*L2*K2
     .    - 4d0*L2*HB2 - 5d0*HB2**2 - HB2*HTAU2 - L2*HTAU2
     .    + 2743d0/162d0*g1**2 + 15d0/2d0*g2**2- 16d0/9d0*g3**2
     .    + 5d0/3d0*g1*g2 + 136d0/27d0*g1*g3 + 8d0*g2*g3)

      F(7)= HB2*(-7d0/9d0*g1 - 3d0*g2 - 16d0/3d0*g3
     .    + L2 + HT2 + 6d0*HB2 + HTAU2)
     .    + c2*HB2*(HB2*(2d0/3d0*g1 + 6d0*g2 + 16d0*g3
     .    - 3d0*L2 - 5d0*HT2 - 22d0*HB2 - 3d0*HTAU2)
     .    + 4d0/3d0*g1*HT2 - 3d0*L2**2 - 2d0*L2*K2
     .    - 4d0*L2*HT2 - 5d0*HT2**2 - 3d0*HTAU2**2 + 2d0*g1*HTAU2
     .    + 1435d0/162d0*g1**2 + 15d0/2d0*g2**2 - 16d0/9d0*g3**2
     .    + 5d0/3d0*g1*g2 + 40d0/27d0*g1*g3 + 8d0*g2*g3)

      F(8)= HTAU2*(-3d0*g1 - 3d0*g2
     .    + L2 + 3d0*HB2 + 4d0*HTAU2)
     .    + c2*HTAU2*(-10d0*HTAU2**2 - 9d0*HTAU2*HB2 - 9d0*HB2**2
     .    - 3d0*HB2*HT2 + HTAU2*(2d0*G1 + 6d0*G2)
     .    + HB2*(-2d0/3d0*G1 + 16*G3)
     .    - L2*(3d0*HTAU2 + 3d0*HT2 + 3d0*L2 + 2d0*K2)
     .    + 75d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2)

      F(9)= 11d0*g1*M1
     .    + c2*g1*(398d0/9d0*g1*M1 + 9d0*g2*(M1+M2)
     .    + 88d0/3d0*g3*(M1+M3) + 26d0/3d0*HT2*(AT-M1)
     .    + 14d0/3d0*HB2*(AB-M1) + 6d0*HTAU2*(ATAU-M1)
     .    + 2d0*L2*(AL-M1))

      F(10)= g2*M2
     .     + c2*g2*(3d0*g1*(M1+M2) + 50d0*g2*M2
     .     + 24d0*g3*(M2+M3) + 6d0*HT2*(AT-M2)
     .     + 6d0*HB2*(AB-M2) + 2d0*HTAU2*(ATAU-M2)
     .     + 2d0*L2*(AL-M2))

      F(11)= -3d0*g3*M3
     .     + c2*g3*(11d0/3d0*g1*(M1+M3) + 9d0*g2*(M2+M3)
     .     + 28d0*g3*M3 + 4d0*HT2*(AT-M3)
     .     + 4d0*HB2*(AB-M3))

      F(12)= 4d0*L2*AL + 2d0*K2*AK + 3d0*HT2*AT
     .     + 3d0*HB2*AB + HTAU2*ATAU + g1*M1 + 3d0*g2*M2
     .     + c2*(-20d0*L2**2*AL - 9d0*L2*HT2*(AL+AT)
     .     - 9d0*L2*HB2*(AL+AB) - 3d0*L2*HTAU2*(AL+ATAU)
     .     - 12d0*L2*K2*(AL+AK) - 18d0*HT2**2*AT
     .     - 18d0*HB2**2*AB - 6d0*HT2*HB2*(AT+AB)
     .     - 6d0*HTAU2**2*ATAU - 16d0*K2**2*AK
     .     + 4d0/3d0*g1*HT2*(AT-M1) - 2d0/3d0*g1*HB2*(AB-M1)
     .     + 2d0*g1*HTAU2*(ATAU-M1) + 16d0*g3*HT2*(AT-M3)
     .     + 16d0*g3*HB2*(AB-M3) + 2d0*g1*L2*(AL-M1)
     .     + 6d0*g2*L2*(AL-M2)
     .     - 23d0*g1**2*M1 - 3d0*g1*g2*(M1+M2) - 15d0*g2**2*M2)

      F(13)= 6d0*(L2*AL + K2*AK)
     .     + c2*(-48d0*K2**2*AK - 24d0*L2*K2*(AL+AK)
     .     - 24d0*L2**2*AL - 18d0*L2*HT2*(AL+AT)
     .     - 18d0*L2*HB2*(AL+AB) - 6d0*L2*HTAU2*(AL+ATAU)
     .     + 6d0*g1*L2*(AL-M1) + 18d0*g2*L2*(AL-M2))

      F(14)= 6d0*HT2*AT + HB2*AB + L2*AL
     .     + 13d0/9d0*g1*M1 + 3d0*g2*M2 + 16d0/3d0*g3*M3
     .     + c2*(-44d0*HT2**2*AT - 5d0*HT2*HB2*(AT+AB)
     .     - 3d0*HT2*L2*(AT+AL) - 10d0*HB2**2*AB
     .     - HB2*HTAU2*(AB+ATAU) - 4d0*HB2*L2*(AB+AL)
     .     - 6d0*L2**2*AL - L2*HTAU2*(AL+ATAU)
     .     - 2d0*L2*K2*(AL+AK)
     .     + 2d0*g1*HT2*(AT-M1) + 6d0*g2*HT2*(AT-M2)
     .     + 16d0*g3*HT2*(AT-M3) + 2d0/3d0*g1*HB2*(AB-M1)
     .     - 2743d0/81d0*g1**2*M1 - 5d0/3d0*g1*g2*(M1+M2)
     .     - 136d0/27d0*g1*g3*(M1+M3) - 15d0*g2**2*M2
     .     - 8d0*g2*g3*(M2+M3) + 32d0/9d0*g3**2*M3)

      F(15)= 6d0*HB2*AB + HT2*AT + HTAU2*ATAU + L2*AL
     .     + 7d0/9d0*g1*M1 + 3d0*g2*M2 + 16d0/3d0*g3*M3
     .     + c2*(-44d0*HB2**2*AB - 5d0*HT2*HB2*(AT+AB)
     .     - 3d0*HB2*HTAU2*(AB+ATAU) - 3d0*HB2*L2*(AB+AL)
     .     - 10d0*HT2**2*AT - 4d0*HT2*L2*(AT+AL)
     .     - 6d0*HTAU2**2*ATAU - 6d0*L2**2*AL - 2d0*L2*K2*(AL+AK)
     .     + 2d0/3d0*g1*HB2*(AB-M1) + 6d0*g2*HB2*(AB-M2)
     .     + 16d0*g3*HB2*(AB-M3) + 4d0/3d0*g1*HT2*(AT-M1)
     .     + 2d0*HTAU2*g1*(ATAU-M1)
     .     - 1435d0/81d0*g1**2*M1 - 5d0/3d0*g1*g2*(M1+M2)
     .     - 40d0/27d0*g1*g3*(M1+M3) - 15d0*g2**2*M2
     .     - 8d0*g2*g3*(M2+M3) + 32d0/9d0*g3**2*M3)

      F(16)= 4d0*HTAU2*ATAU + 3d0*HB2*AB + L2*AL
     .     + 3d0*g1*M1 + 3d0*g2*M2
     .     + c2*(-20d0*HTAU2**2*ATAU - 9d0*HTAU2*HB2*(ATAU+AB)
     .     - 3d0*HTAU2*L2*(ATAU+AL) - 18d0*HB2**2*AB
     .     - 3d0*HT2*HB2*(AT+AB) - 6d0*L2**2*AL - 2d0*L2*K2*(AL+AK)
     .     + 2d0*HTAU2*g1*(ATAU-M1) + 6d0*HTAU2*g2*(ATAU-M2)
     .     - 2d0/3d0*HB2*g1*(AB-M1) + 16d0*HB2*g3*(AB-M3)
     .     - 75d0*g1**2*M1 - 3d0*g1*g2*(M1+M2) - 15d0*g2**2*M2)

      F(17)= L2*(MH1+MH2+MS+AL**2) + 3d0*HT2*(MH1+MQ3+MU3+AT**2)
     .       - g1*M1**2 - 3d0*g2*M2**2 + S/2d0
     .       + C2*(-18d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 3d0*HT2*HB2*(MH1+MH2+2d0*MQ3+MU3+MD3+(AT+AB)**2)
     .       - 6d0*L2**2*(MH1+MH2+MS+2d0*AL**2)
     .       - 3d0*L2*HB2*(MH1+2d0*MH2+MS+MQ3+MD3+(AL+AB)**2)
     .       - L2*HTAU2*(MH1+2d0*MH2+MS+ML3+ME3+(AL+ATAU)**2)
     .       - 2d0*L2*K2*(MH1+MH2+4d0*MS+(AL+AK)**2)
     .       + 4d0/3d0*G1*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1*(M1-AT))
     .       + 16d0*G3*HT2*(MH1+MQ3+MU3+AT**2+2d0*M3*(M3-AT))
     .       + 69d0/2d0*G1**2*M1**2 + 33d0/2d0*G2**2*M2**2
     .     + 3d0*G1*G2*(M1**2+M2**2+M1*M2)
     .       + G1*SP + G1/2d0*SIG1 + 3d0/2d0*G2*SIG2)

      F(18)= L2*(MH1+MH2+MS+AL**2) + 3d0*HB2*(MH2+MQ3+MD3+AB**2)
     .       + HTAU2*(MH2+ML3+ME3+ATAU**2)
     .       - g1*M1**2 - 3d0*g2*M2**2 - S/2d0
     .       + C2*(-18d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 6d0*HTAU2**2*(MH2+ML3+ME3+2d0*ATAU**2)
     .       - 3d0*HT2*HB2*(MH1+MH2+2d0*MQ3+MU3+MD3+(AT+AB)**2)
     .       - 6d0*L2**2*(MH1+MH2+MS+2d0*AL**2)
     .       - 3d0*L2*HT2*(2d0*MH1+MH2+MS+MQ3+MU3+(AL+AT)**2)
     .       - 2d0*L2*K2*(MH1+MH2+4d0*MS+(AL+AK)**2)
     .       - 2d0/3d0*G1*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1*(M1-AB))
     .       + 2d0*G1*HTAU2*(MH2+ML3+ME3+ATAU**2+2d0*M1*(M1-ATAU))
     .       + 16d0*G3*HB2*(MH2+MQ3+MD3+AB**2+2d0*M3*(M3-AB))
     .       + 69d0/2d0*G1**2*M1**2 + 33d0/2d0*G2**2*M2**2
     .     + 3d0*G1*G2*(M1**2+M2**2+M1*M2)
     .       - G1*SP + G1/2d0*SIG1 + 3d0/2d0*G2*SIG2)

      F(19)= 2d0*L2*(MH1+MH2+MS+AL**2) + 2d0*K2*(3d0*MS+AK**2)
     .     + C2*(-16d0*K2**2*(3d0*MS+2d0*AK**2)
     .     - 8d0*L2**2*(MH1+MH2+MS+2d0*AL**2)
     .       - 8d0*L2*K2*(MH1+MH2+4d0*MS+(AL+AK)**2)
     .       - 6d0*L2*HT2*(2d0*MH1+MH2+MS+MQ3+MU3+(AL+AT)**2)
     .       - 6d0*L2*HB2*(MH1+2d0*MH2+MS+MQ3+MD3+(AL+AB)**2)
     .       - 2d0*L2*HTAU2*(MH1+2d0*MH2+MS+ML3+ME3+(AL+ATAU)**2)
     .       + 2d0*G1*L2*(MH1+MH2+MS+AL**2+2d0*M1*(M1-AL))
     .       + 6d0*G2*L2*(MH1+MH2+MS+AL**2+2d0*M2*(M2-AL)))

      F(20)= HT2*(MH1+MQ3+MU3+AT**2) + HB2*(MH2+MQ3+MD3+AB**2)
     .       - g1*M1**2/9d0 - 3d0*g2*M2**2
     .       - 16d0/3d0*g3*M3**2 + S/6d0
     .       + C2*(-10d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 10d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+(AB+ATAU)**2)
     .       - L2*HT2*(2d0*MH1+MH2+MS+MQ3+MU3+(AL+AT)**2)
     .       - L2*HB2*(MH1+2d0*MH2+MS+MQ3+MD3+(AL+AB)**2)
     .       + 4d0/3d0*G1*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1*(M1-AT))
     .       + 2d0/3d0*G1*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1*(M1-AB))
     .       + 199d0/54d0*G1**2*M1**2 + 33d0/2d0*G2**2*M2**2
     .     - 64d0/3d0*G3**2*M3**2
     .     + 1d0/3d0*G1*G2*(M1**2+M2**2+M1*M2)
     .       + 16d0/27d0*G1*G3*(M1**2+M3**2+M1*M3)
     .     + 16d0*G2*G3*(M2**2+M3**2+M2*M3)
     .       + 1d0/3d0*G1*SP + 1d0/18d0*G1*SIG1
     .       + 3d0/2d0*G2*SIG2 + 8d0/3d0*G3*SIG3)

      F(21)= 2d0*HT2*(MH1+MQ3+MU3+AT**2)
     .       - 16d0/9d0*g1*M1**2 - 16d0/3d0*g3*M3**2 - 2d0*S/3d0
     .       + C2*(-16d0*HT2**2*(MH1+MQ3+MU3+2d0*AT**2)
     .       - 2d0*HT2*HB2*(MH1+MH2+2d0*MQ3+MU3+MD3+(AT+AB)**2)
     .       - 2d0*HT2*L2*(2d0*MH1+MH2+MS+MQ3+MU3+(AT+AL)**2)
     .       - 2d0/3d0*G1*HT2*(MH1+MQ3+MU3+AT**2+2d0*M1*(M1-AT))
     .       + 6d0*G2*HT2*(MH1+MQ3+MU3+AT**2+2d0*M2*(M2-AT))
     .       + 1712d0/27d0*G1**2*M1**2 - 64d0/3d0*G3**2*M3**2
     .       + 256d0/27d0*G1*G3*(M1**2+M3**2+M1*M3)
     .     - 4d0/3d0*G1*SP + 8d0/9d0*G1*SIG1 + 8d0/3d0*G3*SIG3)

      F(22)= 2d0*HB2*(MH2+MQ3+MD3+AB**2)
     .       - 4d0/9d0*g1*M1**2 - 16d0/3d0*g3*M3**2 + S/3d0
     .       + C2*(-16d0*HB2**2*(MH2+MQ3+MD3+2d0*AB**2)
     .       - 2d0*HB2*HT2*(MH1+MH2+2d0*MQ3+MU3+MD3+(AB+AT)**2)
     .       - 2d0*HB2*HTAU2*(2d0*MH2+MQ3+MD3+ML3+ME3+(AB+ATAU)**2)
     .       - 2d0*HB2*L2*(MH1+2d0*MH2+MS+MQ3+MD3+(AB+AL)**2)
     .       + 2d0/3d0*G1*HB2*(MH2+MQ3+MD3+AB**2+2d0*M1*(M1-AB))
     .       + 6d0*G2*HB2*(MH2+MQ3+MD3+AB**2+2d0*M2*(M2-AB))
     .       + 404d0/27d0*G1**2*M1**2 - 64d0/3d0*G3**2*M3**2
     .       + 64d0/27d0*G1*G3*(M1**2+M3**2+M1*M3)
     .     + 2d0/3d0*G1*SP + 2d0/9d0*G1*SIG1 + 8d0/3d0*G3*SIG3)

      F(23)= -g1*M1**2/9d0 - 3d0*g2*M2**2
     .       - 16d0/3d0*g3*M3**2 + S/6d0
     .       + C2*(199d0/54d0*G1**2*M1**2 + 33d0/2d0*G2**2*M2**2
     .     - 64d0/3d0*G3**2*M3**2
     .     + 1d0/3d0*G1*G2*(M1**2+M2**2+M1*M2)
     .       + 16d0/27d0*G1*G3*(M1**2+M3**2+M1*M3)
     .     + 16d0*G2*G3*(M2**2+M3**2+M2*M3)
     .       + 1d0/3d0*G1*SP + 1d0/18d0*G1*SIG1
     .       + 3d0/2d0*G2*SIG2 + 8d0/3d0*G3*SIG3)

      F(24)= -16d0/9d0*g1*M1**2 - 16d0/3d0*g3*M3**2 - 2d0*S/3d0
     .       + C2*(1712d0/27d0*G1**2*M1**2 - 64d0/3d0*G3**2*M3**2
     .       + 256d0/27d0*G1*G3*(M1**2+M3**2+M1*M3)
     .     - 4d0/3d0*G1*SP + 8d0/9d0*G1*SIG1 + 8d0/3d0*G3*SIG3)

      F(25)= -4d0/9d0*g1*M1**2 - 16d0/3d0*g3*M3**2 + S/3d0
     .       + C2*(404d0/27d0*G1**2*M1**2 - 64d0/3d0*G3**2*M3**2
     .       + 64d0/27d0*G1*G3*(M1**2+M3**2+M1*M3)
     .     + 2d0/3d0*G1*SP + 2d0/9d0*G1*SIG1 + 8d0/3d0*G3*SIG3)

      F(26)= HTAU2*(MH2+ML3+ME3+ATAU**2)
     .       - g1*M1**2 - 3d0*g2*M2**2 - S/2d0
     .       + C2*(-6d0*HTAU2**2*(MH2+ML3+ME3+2d0*ATAU**2)
     .       - 3d0*HTAU2*HB2*(2d0*MH2+ML3+ME3+MQ3+MD3+(ATAU+AB)**2)
     .       - HTAU2*L2*(MH1+2d0*MH2+MS+ML3+ME3+(ATAU+AL)**2)
     .       + 2d0*G1*HTAU2*(MH2+ML3+ME3+ATAU**2+2d0*M1*(M1-ATAU))
     .       + 69d0/2d0*G1**2*M1**2 + 33d0/2d0*G2**2*M2**2
     .     + 3d0*G1*G2*(M1**2+M2**2+M1*M2)
     .       - G1*SP + G1/2d0*SIG1 + 3d0/2d0*G2*SIG2)

      F(27)= 2d0*HTAU2*(MH2+ML3+ME3+ATAU**2) - 4d0*g1*M1**2 + S
     .       + C2*(-8d0*HTAU2**2*(MH2+ML3+ME3+2d0*ATAU**2)
     .       - 6d0*HTAU2*HB2*(2d0*MH2+ML3+ME3+MQ3+MD3+(ATAU+AB)**2)
     .       - 2d0*HTAU2*L2*(MH1+2d0*MH2+MS+ML3+ME3+(ATAU+AL)**2)
     .       - 2d0*G1*HTAU2*(MH2+ML3+ME3+ATAU**2+2d0*M1*(M1-ATAU))
     .       + 6d0*G2*HTAU2*(MH2+ME3+ML3+ATAU**2+2d0*M2*(M2-ATAU))
     .       + 156d0*G1**2*M1**2 + 2d0*G1*SP + 2d0*G1*SIG1)

      F(28)= -g1*M1**2 - 3d0*g2*M2**2 - S/2d0
     .       + C2*(69d0/2d0*G1**2*M1**2 + 33d0/2d0*G2**2*M2**2
     .     + 3d0*G1*G2*(M1**2+M2**2+M1*M2)
     .       - G1*SP + G1/2d0*SIG1 + 3d0/2d0*G2*SIG2)

      F(29)= -4d0*g1*M1**2 + S
     .       + C2*(156d0*G1**2*M1**2 + 2d0*G1*SP + 2d0*G1*SIG1)

      F(30)= XIF*(L2+K2)
     .     + c2*XIF*(L2*(g1 + 3d0*g2 - 3d0*HT2 - 3d0*HB2
     .     - HTAU2- 2d0*L2 - 4d0*K2) - 4d0*K2**2)
      
      F(31)= L2*(XIS + 2d0*AL*XIF) + K2*(XIS + 2d0*AK*XIF)
     .     + 2d0*LL*M3H*(AL + MUP)
     .     + KK*(MSP*(AK + MUP) + 2d0*MUP*MS)
     .     - c2*(4d0*K2**2*(XIS + 4d0*AK*XIF)
     .     + 2d0*L2**2*(XIS + 4d0*AL*XIF)
     .     + 4d0*L2*K2*(XIS + 2d0*(AL+AK)*XIF)
     .     + 3d0*L2*HT2*(XIS + 2d0*(AL+AT)*XIF)
     .     + 3d0*L2*HB2*(XIS + 2d0*(AL+AB)*XIF)
     .     + L2*HTAU2*(XIS + 2d0*(AL+ATAU)*XIF)
     .     + 6d0*LL*HT2*M3H*(MUP + AL + AT)
     .     + 6d0*LL*HB2*M3H*(MUP + AL + AB)
     .     + 2d0*LL*HTAU2*M3H*(MUP + AL + ATAU)
     .     + 4d0*LL*L2*M3H*(MUP + 2d0*AL)
     .     + 4d0*KK*L2*((MSP + MUP*AL)*(MUP + AL + AK)
     .     + (3d0*MS + MH1 + MH2)*MUP)
     .     + 4d0*KK*K2*((MSP + MUP*AK)*(MUP + 2d0*AK)
     .     + 5d0*MS*MUP)
     .     + LL*G1*(2D0*M3H*(AL+MUP-M1)+LL*(2D0*XIF*(AL-M1)+XIS))
     .     + 3D0*LL*G2*(2D0*M3H*(AL+MUP-M2)+LL*(2D0*XIF*(AL-M2)+XIS)))

      F(32)= MUP*2d0*(L2+K2)
     .     + c2*MUP*2d0*(L2*(g1 + 3d0*g2 - 3d0*HT2 - 3d0*HB2
     .     - HTAU2- 2d0*L2 - 4d0*K2) - 4d0*K2**2)
      
      F(33)= 2d0*MSP*(L2+2d0*K2) + 4d0*MUP*(L2*AL+K2*AK)
     .     + 4d0*LL*KK*M3H + c2*(
     .     - 4d0*L2**2*(MSP + 4d0*MUP*AL)
     .     - 8d0*K2**2*(2d0*MSP + 5d0*MUP*AK)
     .     - 8d0*L2*K2*(2d0*MSP + MUP*(3d0*AL + 2d0*AK))
     .     - 6d0*L2*HT2*(MSP + 2d0*MUP*(AL + AT))
     .     - 6d0*L2*HB2*(MSP + 2d0*MUP*(AL + AB))
     .     - 2d0*L2*HTAU2*(MSP + 2d0*MUP*(AL + ATAU))
     .     - 4d0*LL*KK*(3d0*HT2 + 3d0*HB2 + HTAU2 + 2d0*L2
     .     - G1 - 3d0*G2)*M3H + 2d0*L2*G1*(MSP + 2d0*MUP*(AL-M1))
     .     + 6d0*L2*G2*(MSP + 2d0*MUP*(AL-M2)))

      F(34)= M3H/2d0*(3d0*HT2 + 3d0*HB2 + HTAU2 + 6d0*L2
     .     - G1 - 3d0*G2) + LL*KK*MSP
     .     + c2*(M3H/2d0*(-9d0*HT2**2 - 9d0*HB2**2 - 3d0*HTAU2**2
     .     - 6d0*HT2*HB2 - 14d0*L2**2 - 15d0*L2*HT2 - 15d0*L2*HB2
     .     - 5d0*L2*HTAU2 - 4d0*L2*K2 + 4d0*L2*G1 + 12d0*L2*G2
     .     + 4d0/3d0*HT2*G1 + 16d0*HT2*G3 - 2d0/3d0*HB2*G1
     .     + 16d0*HB2*G3 + 2d0*HTAU2*G1 + 23d0/2d0*G1**2
     .     + 3d0*G1*G2 + 15d0/2d0*G2**2)
     .     - 4d0*LL*KK*(L2+K2)*MSP
     .     - 4d0*LL*KK*(L2*AL+K2*AK)*MUP)

      F(35)= 3d0*HB2*AB + HTAU2*ATAU + L2*AL
     .     + 3d0*g1*M1 + 3d0*g2*M2
     .     + c2*(-6d0*HTAU2**2*ATAU - 3d0*HT2*HB2*(AT+AB)
     .     - 18d0*AB*HB2**2 - 6d0*L2**2*AL
     .     - 3d0*L2*HT2*(AL+AT) - 2d0*L2*K2*(AL+AK)
     .     + 2d0*HTAU2*g1*(ATAU-M1) - 2d0/3d0*HB2*g1*(AB-M1)
     .     + 16d0*HB2*g3*(AB-M3)
     .     - 75d0*g1**2*M1 - 3d0*g1*g2*(M1+M2) - 15d0*g2**2*M2)

      RETURN
      END
