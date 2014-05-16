      SUBROUTINE CHECKMIN(PROB)

**********************************************************************      
* Subroutine to check whether the physical minimum of the effective
* potential (<h1>, <h2> and <s> =/= 0) is deeper than minima with
* <h1>, <h2> or <s> = 0
*
* If not: PROB(28) =/= 0
*
* The effective potential includes 1 loop contributions (large logs
* + finite) from (s)top and (s)bottom loops
*
* The soft masses squared mh1q, mh2q and msq (at the scale QSTSB)
* are computed here and stored in COMMON/QMHIGGS (in NMHDECAY)
* or directly taken from COMMON/QMHIGGS (in NMSPEC and NMGMSB)
* If mh1q, mh2q >> Q2 then PROB(29) =/= 0
*
**********************************************************************

      IMPLICIT NONE

      INTEGER OMGFLAG,MAFLAG

      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION QSTSB,pi,cc,mQ3,mU3,mD3,At,Ab
      DOUBLE PRECISION mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
      DOUBLE PRECISION Mstop,dMst,Wt
      DOUBLE PRECISION Msbot,dMsb,Wb
      DOUBLE PRECISION ct,fmt1,fmt2,fmt,gmt
      DOUBLE PRECISION cb,fmb1,fmb2,fmb,gmb
      DOUBLE PRECISION V,V0,V1,V2,V3,D,mh1q,mh2q,msq
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION Bq,EPS
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ,Q2
      DOUBLE PRECISION S1,S2,S3,S4,S5,S6,S7
      DOUBLE PRECISION XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      DOUBLE PRECISION A,B,C,P,Q,DET,AUX,XM
      DOUBLE PRECISION V31,V32,PHI,SIGQ,R

      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/STSBSCALE/QSTSB
      COMMON/QMHIGGS/MH1Q,MH2Q,MSQ
      COMMON/RADCOR/mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
      COMMON/RADCOR2/MQ3,MU3,MD3,AT,AB
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QEXT/XIFQ,XISQ,MUPQ,MSPQ,M3HQ
      COMMON/RENSCALE/Q2

      pi=4d0*DATAN(1d0)
      cc=3d0/(32d0*pi**2)
      EPS=-1.D-2
      
* Store previous values:
      S1=MST1
      S2=MST2
      S3=MSB1
      S4=MSB2
      S5=MTOPQ
      S6=MBOTQ
      S7=NUQ
      
*   Parameters for (S)top/(S)bottom rad. corrs. in wrong minima:

      Mstop=(mQ3+mU3)/2d0
      dMst=(mQ3-mU3)/2d0
      Msbot=(mQ3+mD3)/2d0
      dMsb=(mQ3-mD3)/2d0
      
      IF(MIN(mst1,msb1).LE.0d0)RETURN

      ct=3d0*htq**2/(16d0*pi**2)
      fmt1=mst1*(DLOG(mst1/QSTSB)-1d0)
      fmt2=mst2*(DLOG(mst2/QSTSB)-1d0)
      fmt=mtopq**2*(DLOG(mtopq**2/QSTSB)-1d0)
      IF(mst1-mst2.NE.0d0)THEN
       gmt=(fmt2-fmt1)/(mst2-mst1)
      ELSE
       gmt=DLOG(mst1/QSTSB)
      ENDIF

      cb=3d0*hbq**2/(16d0*pi**2)
      fmb1=msb1*(DLOG(msb1/QSTSB)-1d0)
      fmb2=msb2*(DLOG(msb2/QSTSB)-1d0)
      fmb=mbotq**2*(DLOG(mbotq**2/QSTSB)-1d0)
      IF(msb1-msb2.NE.0d0)THEN
       gmb=(fmb2-fmb1)/(msb2-msb1)
      ELSE
       gmb=DLOG(msb1/QSTSB)
      ENDIF

*   Soft masses

      IF(MAFLAG.GE.0)THEN
       Bq=Alq+nuq
       MH1Q=-lq**2*h2q**2 - muq**2 + muq*Bq/tanbq
     .      +gq/2d0*(h2q**2-h1q**2)
     .      -ct*(fmt1+fmt2-2d0*fmt+At*Xt*gmt)
     .      +cb*muq/tanbq*Xb*gmb
       MH2Q=-lq**2*h1q**2 - muq**2 + muq*Bq*tanbq
     .      +gq/2d0*(h1q**2-h2q**2)
     .      +ct*muq*tanbq*Xt*gmt
     .      -cb*(fmb1+fmb2-2d0*fmb+Ab*Xb*gmb)
       MSQ=-LQ**2*(H1Q**2+H2Q**2) - 2d0*NUQ**2
     .     +LQ**2*H1Q*H2Q/MUQ*(ALQ+2d0*NUQ+MUPQ) - NUQ*AKQ
     .     -XIFQ*(2d0*KQ+LQ*MUPQ/MUQ)-MUPQ**2-3d0*MUPQ*NUQ
     .     -MSPQ-LQ*XISQ/MUQ
     .     +ct*LQ**2*H1Q*H2Q/MUQ*Xt*gmt
     .     +cb*LQ**2*H1Q*H2Q/MUQ*Xb*gmb
      ENDIF

*   Physical minimum

      V= muq**2*(h1q**2+h2q**2) + lq**2*h1q**2*h2q**2
     .  + nuq**2*muq**2/lq**2 + mupq**2*muq**2/lq**2
     .  - 2d0*nuq*muq*h1q*h2q - 2d0*lq*xifq*h1q*h2q
     .  - 2d0*mupq*muq*h1q*h2q + 2d0*xifq*kq*muq**2/lq**2
     .  + 2d0*mupq*kq*muq**3/lq**3 + 2d0*xifq*mupq*muq/lq
     .  - 2d0*(Alq*muq+m3hq)*h1q*h2q
     .  + 2d0/3d0*Akq*kq*muq**3/lq**3
     .  + 2d0*xisq*muq/lq + mspq*muq**2/lq**2
     .  + gq/4d0*(h1q**2-h2q**2)**2
     .  + mh1q*h1q**2 + mh2q*h2q**2 + msq*muq**2/lq**2
     .  + cc*
     .  ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .  + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .  - 2d0*mtopq**4*(DLOG(mtopq**2/QSTSB)-1.5d0)
     .  + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .  + msb2**2*(DLOG(msb2/QSTSB)-1.5d0)
     .  - 2d0*mbotq**4*(DLOG(mbotq**2/QSTSB)-1.5d0))

*   Minimum with h1=h2=s=0

      mst1=mU3
      mst2=mQ3
      msb1=mD3
      msb2=mQ3
      V0= cc*
     .  ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .  + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .  + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .  + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
      IF(V.NE.0d0)THEN
       PROB(28)=DDIM(EPS,(V0-V)/DABS(V))
      ELSE
       PROB(28)=DDIM(EPS,V0)
      ENDIF

*   Minimum with h2=s=0

      V1=0d0
      IF(mh1q.LT.0d0)THEN
       mtopq=htq*DSQRT(-2d0*mh1q/gq)
       Wt=DSQRT(dMst**2+mtopq**2*At**2)
       mst1=Mstop+mtopq**2-Wt
       mst2=Mstop+mtopq**2+Wt
       msb1=mD3
       msb2=mQ3
       IF(mst1.GT.0d0)THEN
        V1= -mh1q**2/gq
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    - 2d0*mtopq**4*(DLOG(mtopq**2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
        IF(V.NE.0d0)THEN
         PROB(28)=PROB(28)+DDIM(EPS,(V1-V)/DABS(V))
        ELSE
         PROB(28)=PROB(28)+DDIM(EPS,V1)
        ENDIF
       ENDIF
      ENDIF

*   Minimum with h1=s=0

      V2=0d0
      IF(mh2q.LT.0d0)THEN
       mst1=mU3
       mst2=mQ3
       mbotq=hbq*DSQRT(-2d0*mh2q/gq)
       Wb=DSQRT(dMsb**2+mbotq**2*Ab**2)
       msb1=Msbot+mbotq**2-Wb
       msb2=Msbot+mbotq**2+Wb
       IF(msb1.GT.0d0)THEN
        V2= -mh2q**2/gq
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0)
     .    - 2d0*mbotq**4*(DLOG(mbotq**2/QSTSB)-1.5d0))
        IF(V.NE.0d0)THEN
         PROB(28)=PROB(28)+DDIM(EPS,(V2-V)/DABS(V))
        ELSE
         PROB(28)=PROB(28)+DDIM(EPS,V2)
        ENDIF
       ENDIF
      ENDIF


*   Minimum with h1=h2=0

      V3=0d0
      mst1=mU3
      mst2=mQ3
      msb1=mD3
      msb2=mQ3
      A=4d0*KQ**2
      B=KQ*(6d0*MUPQ+2.*AKQ)
      C=2d0*MUPQ**2+4d0*KQ*XIFQ+2d0*MSQ+2d0*MSPQ
      D=2d0*MUPQ*XIFQ+2d0*XISQ
      IF(A.EQ.0d0)THEN
       IF(D.NE.0d0)THEN
        XM=-C/D
        V3=C*XM**2/2d0+D*XM
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
       ENDIF
      ELSE
       P=C/(3d0*A)-B**2/(9d0*A**2)
       Q=B**3/(27d0*A**3)-B*C/(6d0*A**2)+D/(2d0*A)
       IF(P.EQ.0d0)THEN
        IF(Q.EQ.0d0)THEN
         XM=-B/(3d0*A)
        ELSE
         SIGQ=Q/DABS(Q)
         AUX=(DABS(Q))**(1d0/3d0)
         XM=SIGQ*AUX-B/(3d0*A)
        ENDIF
        V3=A*XM**4/4d0+B*XM**3/3d0+C*XM**2/2d0+D*XM
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
       ELSEIF(Q.EQ.0d0)THEN
        XM=-B/(3d0*A)
        V3=A*XM**4/4d0+B*XM**3/3d0+C*XM**2/2d0+D*XM
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
        IF(P.LT.0d0)THEN
         XM=DSQRT(-P)-B/(3d0*A)
         V31=A*XM**4/4d0+B*XM**3/3d0+C*XM**2/2d0+D*XM
     .      + cc*
     .      ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .      + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .      + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .      + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
         XM=-DSQRT(-P)-B/(3d0*A)
         V32=A*XM**4/4d0+B*XM**3/3d0+C*XM**2/2d0+D*XM
     .      + cc*
     .      ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .      + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .      + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .      + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
         V3=MIN(V3,V31,V32)
        ENDIF
       ELSE
        DET=Q**2+P**3
        SIGQ=Q/DABS(Q)
        R=SIGQ*DSQRT(DABS(P))
        IF(DET.GT.0d0) THEN
         AUX=(DABS(Q)+DSQRT(DET))**(1d0/3d0)
         XM=SIGQ*(P/AUX-AUX)-B/(3d0*A)
         V3=A*XM**4/4d0+B*XM**3/3d0+C*XM**2/2d0+D*XM
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
            ELSE
          PHI=DACOS(DABS(Q)/(DABS(P))**(3d0/2d0))
          XM=-2d0*R*DCOS(PHI/3d0)-B/(3d0*A)
          V31=A*XM**4/4d0+B*XM**3/3d0+C*XM**2/2d0+D*XM
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
          XM=2d0*R*DCOS((PHI-PI)/3d0)-B/(3d0*A)
          V32=A*XM**4/4d0+B*XM**3/3d0+C*XM**2/2d0+D*XM
     .    + cc*
     .    ( mst1**2*(DLOG(mst1/QSTSB)-1.5d0)
     .    + mst2**2*(DLOG(mst2/QSTSB)-1.5d0)
     .    + msb1**2*(DLOG(msb1/QSTSB)-1.5d0)
     .    + msb2**2*(DLOG(msb2/QSTSB)-1.5d0))
             V3=MIN(V31,V32)
         ENDIF
        ENDIF
        ENDIF
          IF(V.NE.0d0)THEN
        PROB(28)=PROB(28)+DDIM(EPS,(V3-V)/DABS(V))
       ELSE
        PROB(28)=PROB(28)+DDIM(EPS,V3)
       ENDIF

      PROB(29)=DDIM(MAX(DABS(MH1Q),DABS(MH2Q))/Q2,10d0)
      
* RETURN previous values:
      MST1=S1
      MST2=S2
      MSB1=S3
      MSB2=S4
      MTOPQ=S5
      MBOTQ=S6
      NUQ=S7

      !WRITE(0,*)"CALL CHECKMIN"
      !WRITE(0,*)""
      !WRITE(0,*)"V =",V
      !WRITE(0,*)"V0 =",V0
      !WRITE(0,*)"V1 =",V1
      !WRITE(0,*)"V2 =",V2
      !WRITE(0,*)"V3 =",V3
      !WRITE(0,*)""
      !WRITE(0,*)""

      END
