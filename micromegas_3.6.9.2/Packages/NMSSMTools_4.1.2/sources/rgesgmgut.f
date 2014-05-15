      SUBROUTINE RGESGMGUT(PROB,IFAIL)

*   Subroutine to integrate the RG equations for the gauge and Yukawa
*   couplings from MMESS up to the GUT scale (which is determined here)
*   through a CALL of the subroutine ODEINTGMGUT that is part of
*   the file integgmgut.f
*
*   It checks whether there is a Landau Pole below M_GUT
*   for the couplings lambda, kappa, htop and hbot
*   If yes: PROB(27) =/= 0
*
***********************************************************************

      IMPLICIT NONE

      INTEGER I,IFAIL,NN
      PARAMETER (NN=8)

      DOUBLE PRECISION PROB(*),EPS,X1,X2,Y(NN)
      DOUBLE PRECISION PI,COEF,YMAX
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     .      HBOTMES,HTAUMES
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     .      HBOTGUT,HTAUGUT
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,MGUT

      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTOPMES,
     .      HBOTMES,HTAUMES
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     .      HBOTGUT,HTAUGUT
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MGUT/MGUT

      EXTERNAL DERIVSGMGUT,RKQSGMGUT

      EPS=1.D-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

* Definition of the couplings squared Y(I) at MMESS

      Y(1)=G1MES
      Y(2)=G2MES
      Y(3)=G3MES
      Y(4)=LMES
      Y(5)=KMES
      Y(6)=HTOPMES
      Y(7)=HBOTMES
      Y(8)=HTAUMES

      X1=0d0
      X2=(3d0/G1MES-5d0/G2MES)/28d0


      CALL ODEINTGMGUT(Y,NN,X1,X2,EPS,DERIVSGMGUT,RKQSGMGUT,IFAIL)

* The GUT scale in GeV:

      MGUT=MMESS*DEXP(8d0*PI**2*X2)

      YMAX=0.
      DO I=1,8
       YMAX=MAX(YMAX,Y(I))
      ENDDO

      PROB(27)=DDIM(YMAX/(4d0*PI),1d0)
      
      IF(IFAIL.GT.0)THEN
       IFAIL=0
*       IFAIL=11
      ELSE
       IFAIL=0
      ENDIF

* Couplings at the GUT scale

      G1GUT=Y(1)
      G2GUT=Y(2)
      G3GUT=Y(3)
      LGUT=Y(4)
      KGUT=Y(5)
      HTOPGUT=Y(6)
      HBOTGUT=Y(7)
      HTAUGUT=Y(8)

      END
