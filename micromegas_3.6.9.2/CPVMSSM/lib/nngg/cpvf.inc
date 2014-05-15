      REAL*8 FUNCTION neutcp1(MZ, MW, SW, MG1, MG2, mue, tb, 
     c   th1, th2, thmu)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8  MZ, MW, SW, MG1, MG2, mue, tb, th1, th2, thmu, nm(4)
      COMPLEX*16 neum(4,4)
      common /neumres/ neum, nm
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN
      
      DO 1 I = 1, 4
      DO 1 J = 1, 4
 1    Neum(I,J) = (0D0, 0D0)

      DO 2 I=1,4
 2    nm(I)=0D0
 

      S2W_CN=SW**2
      MW_CN=MW
      MZ_CN=MZ
	  
      call  neutdiag(MG1, MG2, mue, tb, th1, th2, thmu, nm, neum)
      
      neutcp1=0D0
      
      end
      
      REAL*8 FUNCTION neutcp2(tk, o1, o2, o3)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 nm(4)
      COMPLEX*16 neum(4,4)
      common /neumres/ neum, nm
      integer i1,i2,i3
	  i1=o1+0.1
	  i2=o2+0.1
	  i3=o3+0.1
	  
      neutcp2=-10000D0
      
      if(i1.eq.1) neutcp2=DREAL(neum(i3,i2))
      if(i1.eq.2) neutcp2=DIMAG(neum(i3,i2))
      if(i1.eq.3) neutcp2=nm(i2)
      
      return 
      end

      SUBROUTINE NEUTDIAG(M1,M2,MU,TANB,TH1,TH2,THMU,MN,N)
C********************************************************************
C* Diagonalizes the complex, symmetric MSSM neutralino mass matrix. *
C* Gaugino mass unification is assumed.                             *
C********************************************************************

      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 Y(4,4),N(4,4),MDIAG(4,4),XI
      DIMENSION AUX(8,8),EV(8),H(8),MN(4)
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN

      MW   = MW_CN
      MZ   = MZ_CN
      CW   = DSQRT(1.D0-S2W_CN)
      SW   = DSQRT(S2W_CN)

c     M1 = 5.D0*S2W_CN*M2/(3.D0*CW*CW)     !U(1) gaugino mass

*      print*,m1,m2,mu,tanb,th1,th2,thmu

C  *** Define complex neutralino mass matrix  ***

      XI = (0.D0,1.D0)

      C1 = DCOS(TH1)
      S1 = DSIN(TH1)

      C2 = DCOS(TH2)
      S2 = DSIN(TH2)

      CMU = DCOS(THMU)
      SMU = DSIN(THMU)

      BETA = DATAN(TANB)
      CB   = DCOS(BETA)
      SB   = DSIN(BETA)

      Y(1,1) = M1*(C1+XI*S1)
      Y(1,2) = (0.D0,0.D0)
      Y(1,3) =-MZ*SW*CB
      Y(1,4) = MZ*SW*SB
      Y(2,2) = M2*(C2+XI*S2)
      Y(2,3) = MW*CB
      Y(2,4) =-MW*SB
      Y(3,3) = 0.D0
      Y(3,4) =-MU*(CMU+XI*SMU)
      Y(4,4) = 0.D0

      DO 1 I = 2, 4
      DO 1 J = 1, I-1
 1    Y(I,J) = Y(J,I)


C  *** Define auxiliary real, symmetric 8x8 matrix AUX ***

      DO 10 I = 1,4
      DO 10 J = 1,4
      AUX(I,J)     = DREAL(Y(I,J))
      AUX(I,J+4)   = DIMAG(Y(I,J))
      AUX(I+4,J)   = DIMAG(Y(I,J))
 10   AUX(I+4,J+4) =-DREAL(Y(I,J))


C  *** Diagonalize AUX; eigenvalues in EV, eigenvectors in AUX  ***

      CALL DIAGRS(8,8,AUX,EV,H,IERR)

c      WRITE(*,11) (EV(K),K=1,8)
c 11   FORMAT(' Eigenvalues of AUX:',/2(4(2x,e11.4),/))
c      WRITE(*,*) '                   Eigenvectors of AUX:'
c      DO 12 I = 1,8
c 12   WRITE(*,13) (AUX(K,I),K=1,8)
c 13   FORMAT(8(1X,F7.4))

      DO 20 I = 1,4
 20   MN(I) = EV(I+4)        !First 4 eigenvalues are negative!


C  *** Define N. Recall that the eigenvectors are the COLUMNS of AUX! ***

      DO 30 I = 1,4
      DO 30 J = 1,4

 30   N(I,J) = AUX(J,4+I) + XI*AUX(J+4,4+I)


C  *** Check diagonalization  ***

      DO 40 I = 1,4
      DO 40 L = 1,4
      MDIAG(I,L) = (0.D0,0.D0)

      DO 41 J = 1,4
      DO 41 K = 1,4
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

      REAL*8 FUNCTION chacp1(MZ, MW, SW, MG2, mue, tb, th2, thmu)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8  MZ, MW, SW, MG2, mue, tb, th2, thmu, cm(2)
      COMPLEX*16 vm(2,2), um(2,2)
      common /chamres/ cm, vm, um
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN
      
      DO 1 I = 1, 2
      DO 1 J = 1, 2
	  vm(I,J) = (0D0, 0D0)
 1    um(I,J) = (0D0, 0D0)

      DO 2 I=1,2
 2    cm(I)=0D0
 

      S2W_CN=SW**2
      MW_CN=MW
      MZ_CN=MZ
	  
      call  chardiag(MG2, mue, tb, th2, thmu, cm, vm, um)
      
      chacp1=0D0
      
      end
      
      REAL*8 FUNCTION chacp2(tk, o1, o2, o3)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 cm(2)
      COMPLEX*16 vm(2,2), um(2,2)
      common /chamres/ cm, vm, um
      integer i1,i2,i3
	  i1=o1+0.1
	  i2=o2+0.1
	  i3=o3+0.1
	  
      chacp2=-10000D0
      
      if(i1.eq.1) chacp2=DREAL(vm(i3,i2))
      if(i1.eq.2) chacp2=DIMAG(vm(i3,i2))
      if(i1.eq.3) chacp2=DREAL(um(i3,i2))
      if(i1.eq.4) chacp2=DIMAG(um(i3,i2))
      if(i1.eq.5) chacp2=cm(i2)
      
      return 
      end

      SUBROUTINE CHARDIAG(M2,MU,TANB,TH2,THMU,MC,U,V)
C******************************************************************
C* Diagonalizes the chargino mass matrix of the MSSM in Gunion    *
C* and Haber notation, allowing the gaugino mass M_2 and higgsino *
C* mass mu to be complex. M2 and MU are the absolute values of    *
C* these quantities, and TH2 and THMU their phases. TANB is the   *
C* ratio of vevs tan(beta). U and V are the diagonalization ma-   *
C* trices, and MC(1) and MC(2) the two eigenvalues, 1 standing    *
C* for the lighter one.                                           *
C******************************************************************


      IMPLICIT REAL*8(A-H,M,O-Z)
      COMPLEX*16 U(2,2), V(2,2), CHECK3(2,2), X(2,2), XI
      DIMENSION MC(2)
      COMMON /WEINBERG/ S2W_CN,MW_CN,MZ_CN

C  *** Eigenvalues ***

      MW   = MW_CN
      MZ   = MZ_CN

      M2SQ = M2*M2
      MUSQ = MU*MU
      MWSQ = MW*MW

      THETA = TH2 + THMU
      BETA = DATAN(TANB)
      CB = DCOS(BETA)
      SB = DSIN(BETA)

      TERM1 = M2SQ + MUSQ + 2.D0*MWSQ
      TERM2 = DSQRT( (M2SQ - MUSQ - 2.D0*MWSQ*DCOS(2.D0*BETA))**2
     &             + 8.D0*MWSQ*( (M2*CB)**2 + (MU*SB)**2
     &                         + M2*MU*DSIN(2.D0*BETA)*DCOS(THETA) ) )

      MC1SQ = .5D0*(TERM1 - TERM2)
      MC2SQ = .5D0*(TERM1 + TERM2)

      IF(MC1SQ.LT.0.D0) THEN
         WRITE(*,*) ' Squared chargino mass is negative!!'
         STOP 77
      ENDIF

      MC(1) = DSQRT(MC1SQ)
      MC(2) = DSQRT(MC2SQ)

C  *** Phases Gamma_L, Gamma_R  ***

      C2 = DCOS(TH2)
      S2 = DSIN(TH2)
      CMU = DCOS(THMU)
      SMU = DSIN(THMU)

      XR = CB*C2*M2 + SB*CMU*MU
      IF(DABS(XR).LT.1.D-5) XR = 1.D-5
      GAMR = DATAN( (CB*S2*M2 - SB*SMU*MU) / XR )

      XL = SB*C2*M2 + CB*CMU*MU
      IF(DABS(XL).LT.1.D-5) XL = 1.D-5
      GAML = DATAN( (SB*S2*M2 - CB*SMU*MU) / XL )

C  *** Mixing angle theta_L, theta_R  ***

      RT2 = DSQRT(2.D0)
      YR = M2*CB*DCOS(GAMR-TH2) + MU*SB*DCOS(GAMR+THMU)
      IF(DABS(YR).LT.1.D-5) YR = 1.D-5
      THR = DATAN( ( M2SQ + 2.D0*(MW*SB)**2 - MC1SQ )
     &           / (RT2*MW*YR)                       )

      YL = M2*SB*DCOS(GAML-TH2) + MU*CB*DCOS(GAML+THMU)
      IF(DABS(YL).LT.1.D-5) YL = 1.D-5
      THL = DATAN( ( M2SQ + 2.D0*(MW*CB)**2 - MC1SQ )
     &           / (RT2*MW*YL)                       )


C  *** Phases gamma_1, gamma_2     ***

      CTHL = DCOS(THL)
      CTHR = DCOS(THR)
      STHL = DSIN(THL)
      STHR = DSIN(THR)
      IF(DABS(CTHL).LT.1.D-7) CTHL = 1.D-7
      CGR = DCOS(GAMR)
      CGL = DCOS(GAML)
      SGR = DSIN(GAMR)
      SGL = DSIN(GAML)

*JSL:06/Jun/06: CG1 and CG2 are in one line (Thanks to Pukhov)
      CG1 = ( M2*C2*CTHR - RT2*MW*CB*CGR*STHR) / (MC(1)*CTHL)
      IF(DABS(CG1).GT.1.0001D0) THEN
        WRITE(*,*) ' |cos(theta_1)| > 1!!'
        GAM1 = 0.D0
      ELSEIF(CG1.GT.1.D0) THEN
        GAM1 = 0.D0
      ELSEIF(CG1.LT.-1.D0) THEN
        GAM1 = 3.141592654D0
      ELSEIF(DABS(1.D0-CG1).LT.1.D-12) THEN
        GAM1 = 0.D0
      ELSE
        GAM1 = DACOS(CG1)
      ENDIF
      SG1 = ( M2*S2*CTHR - RT2*MW*CB*SGR*STHR ) / (MC(1)*CTHL)
      IF(SG1.LT.0.D0) GAM1 = -GAM1

      CG2 = ( MU*CMU*CTHR + RT2*MW*SB*CGR*STHR) / (MC(2)*CTHL)
      IF(DABS(CG2).GT.1.0001D0) THEN
        WRITE(*,*) ' |cos(gamma_2)| > 1!!'
        GAM2 = 0.D0
      ELSEIF(CG2.GT.1.D0) THEN
        GAM2 = 0.D0
      ELSEIF(CG2.LT.-1.D0) THEN
        GAM2 = 3.141592654D0
      ELSEIF(DABS(1.D0-CG2).LT.1.D-12) THEN
        GAM2 = 0.D0
      ELSE
        GAM2 = DACOS(CG2)
      ENDIF
      SG2 = ( MU*SMU*CTHR - RT2*MW*SB*SGR*STHR ) / (MC(2)*CTHL)
      IF(SG2.LT.0.D0) GAM2 = -GAM2

C  *** Define diagonalizing matrices U, V  ***

      XI = (0.D0,1.D0)
      
      U(1,1) = CTHR
      U(1,2) = -STHR*( CGR - XI*SGR )
      U(2,1) = STHR*( CGR + XI*SGR )
      U(2,2) = CTHR

      V(1,1) = CTHL*( CG1 + XI*SG1 )
      V(1,2) = -STHL*( DCOS(GAM1-GAML) + XI*DSIN(GAM1-GAML) )
      V(2,1) = STHL*( DCOS(GAML+GAM2) + XI*DSIN(GAML+GAM2) )
      V(2,2) = CTHL*( CG2 + XI*SG2 )
*      print*,'JSLEE:1-CG1,1-CG2',1.D0-CG1,1.D0-CG2
*      print*,'JSLEE:GAMR,GAML,GAM1,GAM2',GAMR,GAML,GAM1,GAM2
*      print*,'JSLEE:V(1,2)',V(1,2)

C  *** Checks  ***

      CHECK1 = CG1*CG1 + SG1*SG1
      IF(DABS(CHECK1-1.D0).GT.1.D-5) THEN
        WRITE(*,*) ' gamma_1 not computed correctly!'
        WRITE(1,22) GAML,GAMR,THL,THR,GAM1,GAM2
      ENDIF

      CHECK2 = CG2*CG2 + SG2*SG2
      IF(DABS(CHECK2-1.D0).GT.1.D-5) THEN
        WRITE(*,*) ' gamma_2 not computed correctly!'
        WRITE(1,22) GAML,GAMR,THL,THR,GAM1,GAM2
      ENDIF
 22   FORMAT(' Gamma_L, Gamma_R = ',2(e11.4,2x),/,
     &       ' theta_L, theta_R = ',2(e11.4,2x),/,
     &       ' gamma_1, gamma_2 = ',2(e11.4,2x))

      X(1,1) = M2*( C2 + XI*S2 )
      X(2,2) = MU*( CMU + XI*SMU )
      X(1,2) = RT2*MW*SB
      X(2,1) = RT2*MW*CB

      DO 10 I = 1,2
      DO 10 L = 1,2
      CHECK3(I,L) = (0.D0,0.D0)

      DO 20 J = 1,2
      DO 20 K = 1,2
 20   CHECK3(I,L) = CHECK3(I,L) 
     &            + DCONJG(U(I,J))*X(J,K)*DCONJG(V(L,K))
 10   CONTINUE

      IF(ABS(CHECK3(1,2)).GT.1.D-3.OR.
     &   ABS(CHECK3(2,1)).GT.1.D-3.OR.
     &   DABS(DIMAG(CHECK3(1,1))).GT.1.D-3.OR.
     &   DABS(DIMAG(CHECK3(2,2))).GT.1.D-3)
     & WRITE(*,*) ' Diagonalization of chargino mass matrix failed!!'

c      WRITE(1,*) '       Diagonalized chargino mass matrix:'
c      DO 30 I = 1,2
c 30   WRITE(1,31) CHECK3(I,1), CHECK3(I,2)
c 31   FORMAT(2(2X,E11.4,' + i*',e11.4))

      RETURN
      END


