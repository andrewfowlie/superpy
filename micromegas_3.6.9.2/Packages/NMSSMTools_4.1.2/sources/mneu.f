      SUBROUTINE NEUTRALINO(PAR)

* Subroutine to compute the neutralino masses MNEU(i) (i=1..5,
* ordered in mass) and the neutralino mixing matrix NEU(i,j)
* 1 loop rad. corrs. are taken into account (assuming degenerate
* squarks/sleptons) as in BMPZ (Bagger et al., hep-ph/9606211).
* It is assumed that the input parameters M1, M2 and mu
* are the running DRbar masses at the scale Q2.
* The auxiliary functions NMB0 and NMB1 are defined in subfun.f

      IMPLICIT NONE

      INTEGER I,J
      DOUBLE PRECISION PAR(*),PI
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION VALP(5),VECP(5,5),EPS
      DOUBLE PRECISION l,k,TanB,mu,M1,M2,v1,v2
      DOUBLE PRECISION ALSMZ,ALP0,GF,g1,g2,S2TW,Q2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MSQ,MSL,SIN2B,M22,MU2,MZ2,MW2,M12
      DOUBLE PRECISION QSTSB,NMB0,NMB1
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTQ,MBQ,MA2
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY

      COMMON/GAUGE/ALSMZ,ALP0,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,NEU
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY

      EPS=1.D-8
      PI=4d0*DATAN(1d0)

      l=par(1)
      k=par(2)
      TanB=PAR(3)
      M1=PAR(20)
      M2=PAR(21)
      MA2=PAR(23)**2
      v2=1d0/DSQRT(2d0*DSQRT(2d0)*(1d0+TanB**2)*GF)
      v1=v2*TanB
      SIN2B=2d0*TANB/(1d0+TANB**2)
      M22=M2**2
      M12=M1**2
      MZ2=MZ**2
      MW2=MW**2

* For 1 loop rad. corrs.:
* Average squark mass squared:
      MSQ=(2d0*PAR(7)+PAR(8)+PAR(9)
     .       +4d0*PAR(15)+2d0*PAR(16)+2d0*PAR(17))/12d0

* Average slepton mass squared:
      MSL=(2d0*PAR(10)+PAR(11)+4d0*PAR(18)
     .       +2d0*PAR(19))/9d0

* First: MU at the scale Q2:

      MU=PAR(4)
      MU2=MU**2

* MU incl. rad. corrs.:

      MU=MU*(1d0-1d0/(64d0*PI**2)*(
     .       12d0*(HTQ**2+HBQ**2)*NMB1(MU2,0d0,QSTSB,Q2)
     .       +(3d0*G2+G1)*(NMB1(MU2,M22,MA2,Q2)+NMB1(MU2,M22,MZ2,Q2)
     .       +2d0*NMB1(MU2,MU2,MZ2,Q2)-4d0*NMB0(MU2,MU2,MZ2,Q2))))

* M1 incl. rad. corrs.:

      M1=M1*(1d0-G1/(16d0*PI**2)*(11d0*NMB1(M12,0d0,MSQ,Q2)
     .     +9d0* NMB1(M12,0d0,MSL,Q2)
     .     +SIN2B*MU/M1*(NMB0(M12,MU2,MA2,Q2)-NMB0(M12,MU2,MZ2,Q2))
     .     +NMB1(M12,MU2,MA2,Q2)+NMB1(M12,MU2,MZ2,Q2)))

* M2 incl. rad. corrs.:

      M2=M2*(1d0-G2/(16d0*PI**2)*(9d0*NMB1(M22,0d0,MSQ,Q2)
     .     +3d0* NMB1(M22,0d0,MSL,Q2)
     .     +SIN2B*MU/M2*(NMB0(M22,MU2,MA2,Q2)-NMB0(M22,MU2,MZ2,Q2))
     .     +NMB1(M22,MU2,MA2,Q2)+NMB1(M22,MU2,MZ2,Q2)
     .     -8d0*NMB0(M22,M22,MW2,Q2)+4d0*NMB1(M22,M22,MW2,Q2)))

      NEU(1,1)=M1
      NEU(1,2)=0d0
      NEU(1,3)=DSQRT(g1/2d0)*v1
      NEU(1,4)=-DSQRT(g1/2d0)*v2
      NEU(1,5)=0d0
      NEU(2,2)=M2
      NEU(2,3)=-DSQRT(g2/2d0)*v1
      NEU(2,4)=DSQRT(g2/2d0)*v2
      NEU(2,5)=0d0
      NEU(3,3)=0d0
      NEU(3,4)=-mu
      NEU(3,5)=-l*v2
      NEU(4,4)=0d0
      NEU(4,5)=-l*v1
      NEU(5,5)=2d0*k/l*mu+MUPSUSY

      CALL DIAGN(5,NEU,VALP,VECP,EPS)
      CALL SORTNA(5,VALP,VECP)
      DO I=1,5
       MNEU(I)=VALP(I)
       DO J=1,5
        NEU(I,J)=VECP(J,I)
       ENDDO
      ENDDO

      END
