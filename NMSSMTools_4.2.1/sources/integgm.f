*   Subroutines to integrate the RGEs for the gauge and Yukawa couplings


      SUBROUTINE
     .       ODEINTGM(YSTART,NVAR,X1,X2,EPS,DERIVS,RKQSGM,IFAIL)

*   Driver subroutine to integrate Ordinary Differential Equations
*   using the Runge-Kutta 5th order method with adaptative step size
*   (from Numerical Recipes)
*   IFAIL=1 stepsize smaller than minimum in ODEINTGM
*   IFAIL=2 too many steps in ODEINTGM
*   IFAIL=3 stepsize underflow in RKQSGM
*   IFAIL=4 max(|dy(i)/dx|*(1+x)/(1+|y(i)|)) > 1/eps^2

      IMPLICIT NONE

      INTEGER NVAR,NMAX,IFAIL,MAXSTP,NSTP,I
      PARAMETER(MAXSTP=10000,NMAX=500)
      
      DOUBLE PRECISION YSTART(NVAR),Y(NMAX),YSCAL(NMAX)
      DOUBLE PRECISION DYDX(NMAX),X,X1,X2,EPS,TINY
      DOUBLE PRECISION H1,HMIN,H,HDID,HNEXT

      EXTERNAL DERIVS,RKQSGM

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
       IF(IFAIL.NE.0)RETURN

       DO I=1,NVAR
        YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
       ENDDO

       IF((X+H-X2)*(X+H-X1).GT.0d0) H=X2-X
      
       CALL
     .    RKQSGM(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS,IFAIL)
      
        IF(IFAIL.GT.0)RETURN
      
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


      SUBROUTINE RKQSGM(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,
     .       HNEXT,DERIVS,IFAIL)

*   Stepper subroutine for ODEINTGM

      IMPLICIT NONE

      INTEGER IFAIL,I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,DYDX(N),Y(N),YSCAL(N)
      DOUBLE PRECISION ERRMAX,H,HTEMP,XNEW,YERR(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION SAFETY,PGROW,PSHRNK,ERRCON

      EXTERNAL DERIVS

      SAFETY=.9d0
      PGROW=-.2d0
      PSHRNK=-.25d0
      ERRCON=(5d0/SAFETY)**(1d0/PGROW)
      H=HTRY

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
