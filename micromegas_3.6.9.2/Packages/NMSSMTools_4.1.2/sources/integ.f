*   Subroutines to integrate the RGEs for the gauge and Yukawa couplings


      SUBROUTINE
     .       ODEINT(YSTART,NVAR,X1,X2,EPS,DERIVS,RKQS,IFAIL)

*   Driver subroutine to integrate Ordinary Differential Equations
*   using the Runge-Kutta 5th order method with adaptative step size
*   (from Numerical Recipes)
*   IFAIL=1 stepsize smaller than minimum in ODEINT
*   IFAIL=2 too many steps in ODEINT
*   IFAIL=3 stepsize underflow in RKQS
*   IFAIL=4 max(|dy(i)/dx|*(1+x)/(1+|y(i)|)) > 1/eps^2

      IMPLICIT NONE

      INTEGER NVAR,NMAX,IFAIL,MAXSTP,NSTP,I
      PARAMETER(MAXSTP=10000,NMAX=500)

      DOUBLE PRECISION YSTART(NVAR),Y(NMAX),YSCAL(NMAX)
      DOUBLE PRECISION DYDX(NMAX),X,X1,X2,EPS,TINY
      DOUBLE PRECISION H1,HMIN,H,HDID,HNEXT

      EXTERNAL DERIVS,RKQS

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

       CALL DERIVS(NVAR,X,Y,DYDX)

       DO I=1,NVAR
        IF(DABS(DYDX(I))*(1d0+X)/(1d0+DABS(Y(I))).GT.EPS**(-2))IFAIL=4
       ENDDO
       IF(IFAIL.GT.0)RETURN

       DO I=1,NVAR
        YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
       ENDDO

       CALL RKQS(Y,DYDX,NVAR,X,X2,H,EPS,YSCAL,HDID,HNEXT,DERIVS,IFAIL)
       IF(IFAIL.GT.0)RETURN

       IF(DABS(Y(2)-5d0/3d0*Y(1)).LT.EPS)THEN
        X2=X
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


      SUBROUTINE RKQS(Y,DYDX,N,X,X2,HTRY,EPS,YSCAL,HDID,
     .       HNEXT,DERIVS,IFAIL)

*   Stepper subroutine for ODEINT

      IMPLICIT NONE

      INTEGER IFAIL,I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,X2,DYDX(N),Y(N),YSCAL(N)
      DOUBLE PRECISION ERRMAX,H,HTEMP,XNEW,YERR(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION SAFETY,PGROW,PSHRNK,ERRCON,G,GTEMP

      EXTERNAL DERIVS

      SAFETY=.9d0
      PGROW=-.2d0
      PSHRNK=-.25d0
      ERRCON=(5d0/SAFETY)**(1d0/PGROW)
      H=HTRY
      G=Y(2)-5d0/3d0*Y(1)

1     CALL RKCK(Y,DYDX,N,X,H,YTEMP,YERR,DERIVS)

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
      ENDIF

      GTEMP=YTEMP(2)-5d0/3d0*YTEMP(1)
      IF(GTEMP.LT.-EPS)THEN
       IFAIL=-1
       X2=X+H
       H=G/(G-GTEMP)*H
       GOTO 1
      ENDIF

      HDID=H
      X=X+H
      DO I=1,N
       Y(I)=YTEMP(I)
      ENDDO

      IF(ERRMAX.GT.ERRCON)THEN
       HNEXT=SAFETY*H*(ERRMAX**PGROW)
      ELSE
        HNEXT=5d0*H
      ENDIF
      IF(X+HNEXT.GT.X2 .AND. IFAIL.EQ.-1) HNEXT=X2-X

      RETURN

      END


      SUBROUTINE RKCK(Y,DYDX,N,X,H,YOUT,YERR,DERIVS)

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
     .       B43=1.2d0, B51=-11d0/54d0, B52=2.5d0, B53=-70d0/27d0,
     .       B54=35d0/27d0, B61=1631d0/55296d0, B62=175d0/512d0,
     .       B63=575d0/13824d0, B64=44275d0/110592d0,
     .       B65=253d0/4096d0, C1=37d0/378d0, C3=250d0/621d0,
     .       C4=125d0/594d0, C6=512d0/1771d0, DC1=C1-2825d0/27648d0,
     .       DC3=C3-18575d0/48384d0, DC4=C4-13525d0/55296d0,
     .       DC5=-277d0/14336d0, DC6=C6-.25d0)

      EXTERNAL DERIVS

      DO I=1,N
       YTEMP(I)=Y(I)+B21*H*DYDX(I)
      ENDDO
      CALL DERIVS(N,X+A2*H,YTEMP,AK2)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B31*DYDX(I)+B32*AK2(I))
      ENDDO
      CALL DERIVS(N,X+A3*H,YTEMP,AK3)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B41*DYDX(I)+B42*AK2(I)+B43*AK3(I))
      ENDDO
      CALL DERIVS(N,X+A4*H,YTEMP,AK4)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B51*DYDX(I)+B52*AK2(I)+B53*AK3(I)+B54*AK4(I))
      ENDDO
      CALL DERIVS(N,X+A5*H,YTEMP,AK5)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B61*DYDX(I)+B62*AK2(I)+B63*AK3(I)+B64*AK4(I)+
     .  B65*AK5(I))
      ENDDO
      CALL DERIVS(N,X+A6*H,YTEMP,AK6)
      DO I=1,N
       YOUT(I)=Y(I)+H*(C1*DYDX(I)+C3*AK3(I)+C4*AK4(I)+C6*AK6(I))
      ENDDO
      DO I=1,N
       YERR(I)=H*(DC1*DYDX(I)+DC3*AK3(I)+DC4*AK4(I)+DC5*AK5(I)+DC6*
     .  AK6(I))
      ENDDO

      RETURN
      END

      SUBROUTINE DERIVS(N,X,Y,F)

*   2-loop Renormalization group equations for g1, g2, g3,
*   lambda, kappa, htop, hbot, htau to be integrated by ODEINT

      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION X,Y(N),F(N),PI,c2
      DOUBLE PRECISION G1,G2,G3,L2,K2,HT2,HB2,HTAU2

      PI=4d0*DATAN(1d0)
      c2=1d0/(16d0*PI**2)

      G1=Y(1)
      G2=Y(2)
      G3=Y(3)
      L2=Y(4)
      K2=Y(5)
      HT2=Y(6)
      HB2=Y(7)
      HTAU2=Y(8)

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

      F(5)= 6d0*K2*(L2 + K2)
     .    + c2*6d0*K2*(-4d0*K2**2 + L2*(g1 + 3d0*g2
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

      RETURN
      END
