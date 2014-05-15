      SUBROUTINE CHARGINO(PAR)

* Subroutine to compute the chargino masses MCH(i) (i=1,2, ordered in
* mass) and the chargino mixing matrices U(i,j), V(i,j)
* 1 loop rad. corrs. are taken into account (assuming degenerate
* squarks/sleptons) as in BMPZ (Bagger et al., hep-ph/9606211)
* It is assumed that the input parameters M2 and mu are the running
* DRbar masses at the scale Q2
* The auxiliary functions NMB0 and NMB1 are defined in subfun.f

      IMPLICIT NONE

      INTEGER I

      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION TanB,mu,M2,g2v1,g2v2,TrX2,detX,D,tU,tV
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION Q2,v1,v2,PI,CU,SU,CV,SV
      DOUBLE PRECISION MSQ,MSL,SIN2B,M22,MU2,MZ2,MW2
      DOUBLE PRECISION QSTSB,NMB0,NMB1,HTQ,HBQ,MTQ,MBQ,MA2
      
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/QQUARK/HTQ,HBQ,MTQ,MBQ
      
      PI=4d0*DATAN(1d0)

      TanB=PAR(3)
      M2=PAR(21)
      MA2=PAR(23)**2
      g2v2=DSQRT(2d0/(1d0+TanB**2))*MW
      g2v1=g2v2*TanB
      v2=1d0/DSQRT(2d0*DSQRT(2d0)*(1d0+TanB**2)*GF)
      v1=v2*TanB
      SIN2B=2d0*TANB/(1d0+TANB**2)
      M22=M2**2
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

* M2 incl. rad. corrs.:

      M2=M2*(1d0-G2/(16d0*PI**2)*(9d0*NMB1(M22,0d0,MSQ,Q2)
     .     +3d0* NMB1(M22,0d0,MSL,Q2)
     .     +SIN2B*MU/M2*(NMB0(M22,MU2,MA2,Q2)-NMB0(M22,MU2,MZ2,Q2))
     .     +NMB1(M22,MU2,MA2,Q2)+NMB1(M22,MU2,MZ2,Q2)
     .     -8d0*NMB0(M22,M22,MW2,Q2)+4d0*NMB1(M22,M22,MW2,Q2)))

      TrX2=2d0*MW**2+M2**2+mu**2
      detX=mu*M2-2d0*MW**2*g2v1*g2v2/(g2v1**2+g2v2**2)
      D=TrX2**2-4d0*detX**2

      tU=(g2v2**2-g2v1**2-M2**2+mu**2-DSQRT(D))/2d0/(M2*g2v2+mu*g2v1)
      tV=(g2v1**2-g2v2**2-M2**2+mu**2-DSQRT(D))/2d0/(M2*g2v1+mu*g2v2)

      CU=DCOS(DATAN(tU))
      SU=DSIN(DATAN(tU))
      CV=DCOS(DATAN(tV))
      SV=DSIN(DATAN(tV))

      U(1,1)=CU
      U(1,2)=SU
      U(2,1)=-SU
      U(2,2)=CU

      V(1,1)=CV
      V(1,2)=SV
      V(2,1)=-SV
      V(2,2)=CV

      DO I=1,2
       MCH(I)=g2v1*U(I,1)*V(I,2)+g2v2*U(I,2)*V(I,1)
     .  +M2*U(I,1)*V(I,1)+mu*U(I,2)*V(I,2)
      ENDDO

      END
