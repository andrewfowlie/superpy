      REAL*8 FUNCTION neut1(MZ, MW, SW, MG1, MG2, vv, tb,hk,hl,xv)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8  MZ, MW, SW, MG1, MG2, mue, tb,hl,hk,xv, nm(5)
      COMPLEX*16 neum(5,5)
      common /neumres/ neum, nm
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN
      
      DO 1 I = 1, 5
      DO 1 J = 1, 5
 1    Neum(I,J) = (0D0, 0D0)

      DO 2 I=1,5
 2    nm(I)=0D0
 
      S2W_CN=SW**2
      MW_CN=MW
      MZ_CN=MZ
      
      call  neutdiag(MG1,MG2,vv,tb,hl,hk,xv,nm, neum)
      
      neut1=0D0
      
      end
      
      REAL*8 FUNCTION neutd(tk, o1, o2, o3)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 nm(5)
      COMPLEX*16 neum(5,5)
      common /neumres/ neum, nm
      integer i1,i2,i3
	  i1=o1+0.1
	  i2=o2+0.1
	  i3=o3+0.1
	  
      neutd=-10000D0
      
      if(i1.eq.1) neutd=DREAL(neum(i3,i2))
      if(i1.eq.2) neutd=DIMAG(neum(i3,i2))
      if(i1.eq.3) neutd=nm(i2)
      
      return 
      end

      SUBROUTINE NEUTDIAG(M1,M2,vv,TANB,hl,hk,xv,MN,N)
C********************************************************************
C* Diagonalizes the complex, symmetric MSSM neutralino mass matrix. *
C* Gaugino mass unification is assumed.                             *
C********************************************************************

      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 Y(5,5),N(5,5),MDIAG(5,5),XI
      DIMENSION AUX(10,10),EV(10),H(10),MN(5)
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN

      MW   = MW_CN
      MZ   = MZ_CN
      CW   = DSQRT(1.D0-S2W_CN)
      SW   = DSQRT(S2W_CN)

c     M1 = 5.D0*S2W_CN*M2/(3.D0*CW*CW)     !U(1) gaugino mass

*      print*,m1,m2,mu,tanb,th1,th2,thmu

C  *** Define complex neutralino mass matrix  ***

      XI = (0.D0,1.D0)


      BETA = DATAN(TANB)
      CB   = DCOS(BETA)
      SB   = DSIN(BETA)

      Y(1,1) = M1
      Y(1,2) = (0.D0,0.D0)
      Y(1,3) =-MZ*SW*CB
      Y(1,4) = MZ*SW*SB
	  Y(1,5) = 0D0
      Y(2,2) = M2
      Y(2,3) = MW*CB
      Y(2,4) =-MW*SB
	  Y(2,5) = 0D0
      Y(3,3) = 0.D0
      Y(3,4) =-xv*hl
	  Y(3,5) =-vv*hl*sb
      Y(4,4) = 0.D0
	  Y(4,5) =-vv*hl*cb
	  Y(5,5) = 2D0*hk*xv

      DO 1 I = 2, 5
      DO 1 J = 1, I-1
 1    Y(I,J) = Y(J,I)


C  *** Define auxiliary real, symmetric 8x8 matrix AUX ***

      DO 10 I = 1,5
      DO 10 J = 1,5
      AUX(I,J)     = DREAL(Y(I,J))
      AUX(I,J+5)   = DIMAG(Y(I,J))
      AUX(I+5,J)   = DIMAG(Y(I,J))
 10   AUX(I+5,J+5) =-DREAL(Y(I,J))


C  *** Diagonalize AUX; eigenvalues in EV, eigenvectors in AUX  ***

      CALL DIAGRS(10,10,AUX,EV,H,IERR)

c      WRITE(*,11) (EV(K),K=1,8)
c 11   FORMAT(' Eigenvalues of AUX:',/2(4(2x,e11.4),/))
c      WRITE(*,*) '                   Eigenvectors of AUX:'
c      DO 12 I = 1,8
c 12   WRITE(*,13) (AUX(K,I),K=1,8)
c 13   FORMAT(8(1X,F7.4))

      DO 20 I = 1,5
 20   MN(I) = EV(I+5)        !First 4 eigenvalues are negative!


C  *** Define N. Recall that the eigenvectors are the COLUMNS of AUX! ***

      DO 30 I = 1,5
      DO 30 J = 1,5

 30   N(I,J) = AUX(J,5+I) + XI*AUX(J+5,5+I)


C  *** Check diagonalization  ***

      DO 40 I = 1,5
      DO 40 L = 1,5
      MDIAG(I,L) = (0.D0,0.D0)

      DO 41 J = 1,5
      DO 41 K = 1,5
 41   MDIAG(I,L) = MDIAG(I,L)
     &           + DCONJG(N(I,J))*Y(J,K)*DCONJG(N(L,K))
 40   CONTINUE


c      DO 50 I=1,4
c	  DO 50 J=1,4
c	  write(*,*) I,J,MDIAG(I,J)
c  50   continue

      RETURN
      END


      SUBROUTINE DIAGRS(NM,N,Z,D,E,IERR)
C*************************************************************
C*  NM: MAXIMAL DIMENSION OF THE MATRIX Z                    *
C*  N : ACTUAL DIMENSION IN THE CALLING PROGRAM              *
C*  E : AUXILIARY VECTOR                                     *
C*  D : VECTOR CONTAINING THE EIGENVALUES                    *
C*  AFTER THE DIAGONALIZATION, THE COLUMNS -                 *
C*  NOT THE ROWS!! - OF THE MATRIX                           *
C*  Z ARE THE EIGENVECTORS                                   *
C*  IERR: ERROR PARAMETER (IERR=0: EVERYTHING OK!)           *
C*************************************************************

C**** EISPACK TRED2
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(NM,N),E(N),D(N)
      IF(N.EQ.1) GOTO 320
      DO 300 II=2,N
      I=N+2-II
      L=I-1
      H=0.D0
      SCALE=0.D0
      IF(L.LT.2) GOTO 130
      DO 120 K=1,L
  120 SCALE=SCALE+DABS(Z(I,K))
      IF(SCALE.NE.0.D0) GOTO 140
  130 E(I)=Z(I,L)
      GOTO 290
  140 DO 150 K=1,L
      Z(I,K)=Z(I,K)/SCALE
      H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
      F=Z(I,L)
      G=-DSIGN(DSQRT(H),F)
      E(I)=SCALE*G
      H=H-F*G
      Z(I,L)=F-G
      F=0.D0
      DO 240 J=1,L
      Z(J,I)=Z(I,J)/H
      G=0.D0
      DO 180 K=1,J
  180 G=G+Z(J,K)*Z(I,K)
      JP1=J+1
      IF(L.LT.JP1) GOTO 220
      DO 200 K=JP1,L
  200 G=G+Z(K,J)*Z(I,K)
  220 E(J)=G/H
      F=F+E(J)*Z(I,J)
  240 CONTINUE
      HH=F/(H+H)
      DO 260 J=1,L
      F=Z(I,J)
      G=E(J)-HH*F
      E(J)=G
      DO 260 K=1,J
      Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
  290 D(I)=H
  300 CONTINUE
  320 D(1)=0.D0
      E(1)=0.D0
      DO 500 I=1,N
      L=I-1
      IF(D(I).EQ.0.D0) GOTO 380
      DO 360 J=1,L
      G=0.D0
      DO 340 K=1,L
  340 G=G+Z(I,K)*Z(K,J)
      DO 360 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
  380 D(I)=Z(I,I)
      Z(I,I)=1.D0
      IF(L.LT.1) GOTO 500
      DO 400 J=1,L
      Z(I,J)=0.D0
      Z(J,I)=0.D0
  400 CONTINUE
  500 CONTINUE

C**** EISPACK IMTQL2

      GENAU = 2.D0  **(-40.D0)
      IERR = 0
      IF (N .EQ. 1) GO TO 5001
      DO 5100 I = 2, N
 5100 E(I-1) = E(I)
      E(N) = 0.D0
      DO 5240 L = 1, N
         J = 0
 5105    DO 5110 M = L, N
            IF (M .EQ. N) GO TO 5120
            IF ( DABS(E(M)) .LE. GENAU * ( DABS(D(M)) +  DABS(D(M+1))))
     X         GO TO 5120
 5110    CONTINUE
 5120    P = D(L)
         IF (M .EQ. L) GO TO 5240
         IF (J .EQ. 30) GO TO 5000
         J = J + 1
         G = (D(L+1) - P) / (2.D0   * E(L))
         R =  DSQRT(G*G+1.D0  )
         G = D(M) - P + E(L) / (G +  DSIGN(R,G))
         S = 1.D0
         C = 1.D0
         P = 0.D0
         MML = M - L
         DO 5200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF ( DABS(F) .LT.  DABS(G)) GO TO 5150
            C = G / F
            R =  DSQRT(C*C+1.D0  )
            E(I+1) = F * R
            S = 1.D0   / R
            C = C * S
            GO TO 5160
 5150       S = F / G
            R =  DSQRT(S*S+1.D0  )
            E(I+1) = G * R
            C = 1.D0   / R
            S = S * C
 5160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.D0   * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
            DO 5180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
 5180       CONTINUE
 5200    CONTINUE
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.D0
         GO TO 5105
 5240 CONTINUE
      DO 5300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 5260 J = II, N
            IF (D(J) .GE. P) GO TO 5260
            K = J
            P = D(J)
 5260    CONTINUE
         IF (K .EQ. I) GO TO 5300
         D(K) = D(I)
         D(I) = P
         DO 5280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
 5280    CONTINUE
 5300 CONTINUE
      GO TO 5001
 5000 IERR = L
 5001 RETURN
      END

