      SUBROUTINE GLUINO(PAR)

* Subroutine to compute the gluino pole mass MGL, taking into account
* gluon/gluino, quark/squark loops as in BMPZ (Bagger et al.,
* hep-ph/9606211), with thanks to S. Kraml.
* It is assumed that the input parameter PAR(22) = M3 is the running
* DRbar gluino mass at the scale Q2.
* The auxiliary functions NMB0 and NMB1 are defined in subfun.f

      IMPLICIT NONE

      DOUBLE PRECISION PAR(*)
      DOUBLE PRECISION MGL,MCH(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,PI,Q2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      DOUBLE PRECISION S2T,S2B,ALSMT,ALSQ
      DOUBLE PRECISION MGL2,DMG,DMQ,R,NMB0,NMB1

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SUSYSPEC/MGL,MCH,U,V,MNEU,N
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RENSCALE/Q2
      
      PI=4d0*DATAN(1d0)
      MGL=PAR(22)
      MGL2=MGL**2
      ALSMT=ALSMZ/(1d0+23d0/(12d0*PI)*ALSMZ*DLOG(MT**2/MZ**2))
      ALSQ=ALSMT/(1d0+7d0*ALSMT/(4d0*pi)*DLOG(Q2/MT**2))
      S2T=2d0*CST*DSQRT(1d0-CST**2)
      S2B=2d0*CSB*DSQRT(1d0-CSB**2)
      R=0d0 ! 1=MSbar, 0=DRbar scheme

*  gluon/gluino correction, eq.(22) of BMPZ,
*  plus possible shift to MSbar

      DMG= 3d0*(2d0*NMB0(MGL2,MGL2,0d0,Q2)-NMB1(MGL2,MGL2,0d0,Q2)
     .     -R/2d0)

*  quark-squark loops, eq.(D.44) of BMPZ

      DMQ= 2d0*(NMB1(MGL2,0d0,MUL**2,Q2) + NMB1(MGL2,0d0,MUR**2,Q2)
     .       + NMB1(MGL2,0d0,MDL**2,Q2) + NMB1(MGL2,0d0,MDR**2,Q2))
     .       + NMB1(MGL2,MT**2,MST1**2,Q2) + NMB1(MGL2,MT**2,MST2**2,Q2)
     .       + NMB1(MGL2,MB**2,MSB1**2,Q2) + NMB1(MGL2,MB**2,MSB2**2,Q2)
     .       + MT/MGL*S2T*(NMB0(MGL2,MT**2,MST1**2,Q2)
     .       + NMB0(MGL2,MT**2,MST2**2,Q2))
     .       + MB/MGL*S2B*(NMB0(MGL2,MB**2,MSB1**2,Q2)
     .       + NMB0(MGL2,MB**2,MSB2**2,Q2))

      MGL=MGL/(1d0-ALSQ/(2d0*PI)*(DMG-DMQ/2d0))

      END
