* DeltaR1L_susy.F
* generated 24 Jun 2010 12:04
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"


	subroutine DeltaR1L_susy(dr)
	implicit none
	double complex dr

#include "FH.h"
#include "looptools.h"

	double complex dup1, dup10, dup11, dup12, dup13, dup14, dup15
	double complex dup16, dup17, dup18, dup19, dup2, dup20, dup21
	double complex dup22, dup23, dup24, dup25, dup26, dup27, dup28
	double complex dup29, dup3, dup30, dup31, dup32, dup33, dup34
	double complex dup35, dup36, dup37, dup38, dup39, dup4, dup40
	double complex dup5, dup6, dup7, dup8, dup9, tmp1, tmp10
	double complex tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp2
	double complex tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
	integer Cha3, Cha4, Cha5, Gen3, Neu3, Neu4, Neu5, Sfe3, Sfe4
	integer Sfe5, Sfe6

	dup1 = Mh02*(Mh02 - 4*MW2) + 12*MW2**2

	dup2 = MHH2*(MHH2 - 4*MW2) + 12*MW2**2

        tmp16 = (((CW2*(4*
     -                (CBA2*
     -                   (3*(MA02 + Mh02 + MHH2)*MW2 - 
     -                     (-MW2 + 
     -                       CW2*(3*(Mh02 + MHH2 + MHp2) + MW2))*
     -                      MZ2) + 
     -                  2*CW2**2*
     -                   (-12*MZ2**2 + 
     -                     MW2*(3*MHp2 + 27*MW2 + 5*MZ2)) + 
     -                  MW2*(3*(MA02 + Mh02 + MHH2) + MZ2)*SBA2 - 
     -                  CW2*MZ2*
     -                   (37*MW2 + 3*(MA02 + MHp2 + MZ2) + 
     -                     (3*(Mh02 + MHH2 + MHp2) + MW2)*SBA2)) + 
     -               (-48*CW2*MW2*(MHp2 + MW2) + 
     -                  (3*(MA02 + Mh02 + MHH2) + 6*MHp2 + 
     -                     (170 - 200*CW2)*MW2)*MZ2 + 
     -                  3*(1 + 8*CW2)*MZ2**2)*SW2) + 
     -            MW2*(4*CW2*(6*MHp2 + 78*MW2 + 143*MZ2)*SW2**2 + 
     -               36*(8*MW2 + 3*MZ2)*SW2**3))/(MW2*MZ2) + 
     -         (6*(-(CW2*(MW2 - MZ2)*
     -                 (CW2*
     -                    (-3*MW2*MZ2**2 + MZ2**3 + 
     -                      CW2*(-24*MW2*MZ2**2 + 8*MZ2**3)) + 
     -                   MW2**2*
     -                    (CBA2*(MHH2 - 2*MZ2) + 
     -                      MZ2*(3 - 2*SBA2) + Mh02*SBA2))) + 
     -              MZ2*(((-1 - 8*CW2)*CW2*MZ2**3 + 
     -                    MW2*
     -                     (CW2*(4 + 44*CW2)*MZ2**2 + 
     -                       MW2*(-30*MW2 + CW2*(72 + 12*CW2)*MZ2))
     -                    )*SW2 + 
     -                 12*MW2**2*(5*MW2 - 6*CW2*MZ2)*SW2**2 - 
     -                 24*MW2**3*SW2**3))*A0(MZ2))/
     -          (MW2**2*(MW2 - MZ2)*MZ2**2))/CW2 + 
     -      (6*((((MA02 - MHp2)*
     -                 (MZ2*
     -                    (CW2*(-MA02 + MHp2)*MZ2 + 
     -                      MW2*(-3*MW2 + 2*CW2*MZ2)) + 
     -                   MW2**2*
     -                    (CBA2*(MA02 - Mh02 + MZ2) + 
     -                      (MA02 - MHH2 + MZ2)*SBA2)) + 
     -                (MHp2*(MHp2 - MW2) + 
     -                   MA02*(MA02 - 2*(MHp2 + MW2)))*MZ2**2*SW2)*
     -              A0(MA02))/((MA02 - MHp2)*MZ2**2) + 
     -           CBA2*(MHp2*(MHp2 - 2*MW2) + MW2**2 + 
     -              Mh02*(Mh02 - 2*(MHp2 + MW2)))*(CW2 - SW2)*
     -            B0(MW2,Mh02,MHp2)))/MW2**2)/SW2**2 + 
     -   6*((4/MHp2 + (-((MA02 + CBA2*Mh02 + MHH2*SBA2)/MW2**2) + 
     -            (-4 + 3*(MA02/(MA02 - MHp2) + 
     -                  (CBA2*Mh02)/(Mh02 - MHp2) + 
     -                  (MHH2*SBA2)/(MHH2 - MHp2)))/MW2)/SW2 + 
     -         ((-4*(CW2 - SW2)**2)/MZ2 + 
     -            (CW2*(MA02 + CBA2*Mh02 + 4*MW2 + MHH2*SBA2) + 
     -               2*MHp2*(-CW2 + SW2))/MW2**2)/SW2**2)*A0(MHp2)+
     -        ((CW2 - SW2)*
     -         ((MHp2*(MHp2 - 2*MW2) + MW2**2 + 
     -              MA02*(MA02 - 2*(MHp2 + MW2)))*B0(MW2,MA02,MHp2)
     -             + (MHp2*(MHp2 - 2*MW2) + MW2**2 + 
     -              MHH2*(MHH2 - 2*(MHp2 + MW2)))*SBA2*
     -            B0(MW2,MHH2,MHp2)))/(MW2**2*SW2**2))

        tmp16 = tmp16 + (6*
     -      (((((CBA2*(Mh02 - MW2)*
     -                    ((Mh02 - MHp2)*
     -                       (MA02*MW2**2 + 
     -                       CW2*(Mh02 - MHp2)*MZ2**2 + 
     -                       MW2*(MW2*(-Mh02 - MZ2) + CW2*MZ2**2))-
     -                        (MHp2*(MHp2 - MW2) + 
     -                       Mh02*(Mh02 - 2*(MHp2 + MW2)))*MZ2**2*
     -                       SW2) + 
     -                   (Mh02 - MHp2)*
     -                    (-((Mh02 - MW2)*
     -                       (-3*MW2*MZ2*(MW2 - CW2*MZ2) + 
     -                       Mh02*(MW2**2 - CW2*MZ2**2)*SBA2)) - 
     -                      dup1*MZ2**2*SBA2*SW2))*A0(Mh02))/
     -               ((Mh02 - MHp2)*MZ2**2) + 
     -              (((Mh02 - MW2)*
     -                    (CW2**2*(MHH2 - MW2)*
     -                       ((MW2 - MZ2)*MZ2*
     -                       (CBA2*(MHH2 - 2*MW2) + MZ2 + 
     -                       MW2*(40 - 2*SBA2) + Mh02*SBA2) + 
     -                       8*
     -                       (MZ2**3 + 
     -                       MW2*(MW2*(MW2 - 13*MZ2) + 5*MZ2**2))*
     -                       SW2) + 
     -                      (-MHH2 + MW2)*
     -                       (4*CW2**3*(9*MW2 - 2*MZ2)*(MW2 - MZ2)*
     -                       (MW2 + MZ2) + 
     -                       6*MW2**2*SW2*
     -                       (8*MW2*SW2**2 + 
     -                       MZ2*(-5 + (10 - 12*SW2)*SW2)))) + 
     -                   CW2*SW2*
     -                    (CBA2*dup2*(-Mh02 + MW2)*(MW2 - MZ2)*
     -                      MZ2 + 
     -                      (MHH2 - MW2)*
     -                       (dup1*MZ2*(-MW2 + MZ2)*SBA2 + 
     -                       (-Mh02 + MW2)*
     -                       (-MZ2**3 + 
     -                       MW2*
     -                       (MZ2**2*(-10 - 80*SW2) + 
     -                       MW2*(MZ2*(86 - 44*SW2) + 52*MW2*SW2)))
     -                       )))*A0(MW2))/
     -               (CW2*(-MHH2 + MW2)*(MW2 - MZ2)*MZ2))/
     -            (-Mh02 + MW2) + 
     -           (CW2 - SW2)*
     -            (dup1*SBA2*B0(MW2,Mh02,MW2) + 
     -              CBA2*dup2*B0(MW2,MHH2,MW2) - 
     -              ((CW2*(4*MW2*MZ2 - MZ2**2 + 
     -                      CW2*(-8*MZ2**2 + MW2*(60*MW2 + 44*MZ2))
     -                      ) - 12*MW2**2*SW2**2)*B0(MW2,MW2,MZ2))/
     -               CW2))/MW2**2 - 
     -        (SBA2*((MHH2*(MHH2 - 2*MZ2) + MZ2**2 + 
     -                 MA02*(MA02 - 2*(MHH2 + MZ2)))*
     -               B0(MZ2,MA02,MHH2) + 
     -              ((CW2*Mh02*(Mh02 - 4*MZ2) + 12*MW2*MZ2)*
     -                 B0(MZ2,Mh02,MZ2))/CW2) + 
     -           CBA2*((Mh02*(Mh02 - 2*MZ2) + MZ2**2 + 
     -                 MA02*(MA02 - 2*(Mh02 + MZ2)))*
     -               B0(MZ2,MA02,Mh02) + 
     -              ((CW2*MHH2*(MHH2 - 4*MZ2) + 12*MW2*MZ2)*
     -                 B0(MZ2,MHH2,MZ2))/CW2))/MZ2**2 + 
     -        ((4*MHp2 - MZ2)*(CW2 - SW2)**2*B0(MZ2,MHp2,MHp2) + 
     -           ((-20*MW2 - MZ2)*SW2**2 + 
     -              CW2*(CW2*(60*MW2 + 39*MZ2) + 
     -                 (-8*MW2 + 2*MZ2)*SW2))*B0(MZ2,MW2,MW2))/MZ2)
     -      )/SW2**2

        tmp16 = tmp16 - (6*
     -      (CBA2*(MHH2 - MHp2)*
     -         (MHH2*(MHH2 - MW2)*(MW2**2 - CW2*MZ2**2) + 
     -           dup2*MZ2**2*SW2) + 
     -        (MHH2 - MW2)*
     -         ((MHH2 - MHp2)*
     -            (MW2**2*
     -               (MZ2*(-3 + SBA2) + (-MA02 + MHH2)*SBA2) + 
     -              CW2*MZ2**2*
     -               (-(MW2*(-3 + SBA2)) + (-MHH2 + MHp2)*SBA2)) + 
     -           (MHp2*(MHp2 - MW2) + 
     -              MHH2*(MHH2 - 2*(MHp2 + MW2)))*MZ2**2*SBA2*SW2))
     -       *A0(MHH2))/
     -    ((MHH2 - MHp2)*MW2**2*(-MHH2 + MW2)*MZ2**2*SW2**2)

	dr = dr + 1/(288.D0*Pi)*(Alfa0*tmp16)

	do Cha3 = 1,2

        dr = dr + 1/(48.D0*Pi)*
     -    (Alfa0*(16*SW2*A0(MCha2(Cha3)) - 
     -         MCha2(Cha3)*
     -          (16*SW2 + (3*
     -                (B0(0.D0,MCha2(Cha3),MSf2(1,1,1)) + 
     -                  B0(0.D0,MCha2(Cha3),MSf2(1,1,2))) + 
     -               DB0(0.D0,MCha2(Cha3),MSf2(1,1,1))*
     -                (3*MCha2(Cha3) - 3*MSf2(1,1,1)) + 
     -               DB0(0.D0,MCha2(Cha3),MSf2(1,1,2))*
     -                (3*MCha2(Cha3) - 3*MSf2(1,1,2)))*
     -             VCha(Cha3,1)*VChaC(Cha3,1))))/(SW2*MCha2(Cha3))

	enddo

	do Gen3 = 1,3

        dr = dr - 1/(144.D0*Pi)*
     -    (Alfa0*(6*(2*MW2 - 3*CW2*MZ2)*A0(MSf2(1,1,Gen3)) + 
     -         MW2*(2*(MZ2 - 6*MSf2(1,1,Gen3)) + 
     -            3*B0(MZ2,MSf2(1,1,Gen3),MSf2(1,1,Gen3))*
     -             (MZ2 - 4*MSf2(1,1,Gen3)))))/(MW2*MZ2*SW2**2)

	enddo

	do Neu3 = 1,4

        dr = dr + 1/(32.D0*Pi)*
     -    (Alfa0*(B0(0.D0,MNeu2(Neu3),MSf2(1,1,1)) + 
     -         B0(0.D0,MNeu2(Neu3),MSf2(1,1,2)) + 
     -         DB0(0.D0,MNeu2(Neu3),MSf2(1,1,1))*
     -          (MNeu2(Neu3) - MSf2(1,1,1)) + 
     -         DB0(0.D0,MNeu2(Neu3),MSf2(1,1,2))*
     -          (MNeu2(Neu3) - MSf2(1,1,2)))*
     -       (SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2))*
     -       (-(SW*ZNeuC(Neu3,1)) + CW*ZNeuC(Neu3,2)))/(CW2*SW2)

	enddo

	do Cha4 = 1,2

	tmp1 = A0(MCha2(Cha4))

	do Cha3 = 1,2

        dup3 = -(1/2.D0*(UCha(Cha4,2)*UChaC(Cha3,2))) + 
     -   SW2*IndexDelta(Cha3,Cha4) - UCha(Cha4,1)*UChaC(Cha3,1)

        dup4 = -(1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2))) + 
     -   SW2*IndexDelta(Cha3,Cha4) - UCha(Cha3,1)*UChaC(Cha4,1)

        dup5 = -(1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2))) + 
     -   SW2*IndexDelta(Cha3,Cha4) - VCha(Cha4,1)*VChaC(Cha3,1)

        dup6 = -(1/2.D0*(VCha(Cha3,2)*VChaC(Cha4,2))) + 
     -   SW2*IndexDelta(Cha3,Cha4) - VCha(Cha3,1)*VChaC(Cha4,1)

        dr = dr + 1/(36.D0*Pi)*
     -    (Alfa0*(-3*B0(MZ2,MCha2(Cha3),MCha2(Cha4))*
     -          (6*(dup3*dup5 + dup4*dup6)*MZ2*MCha(Cha3)*
     -             MCha(Cha4) + 
     -            (dup3*dup4 + dup5*dup6)*
     -             (2*MZ2**2 - MCha2(Cha3)**2 - 
     -               MCha2(Cha3)*(MZ2 - 2*MCha2(Cha4)) - 
     -               MZ2*MCha2(Cha4) - MCha2(Cha4)**2)) + 
     -         (dup3*dup4 + dup5*dup6)*
     -          (3*B0(0.D0,MCha2(Cha3),MCha2(Cha4))*
     -             (MCha2(Cha3) - MCha2(Cha4))*
     -             (2*MZ2 - MCha2(Cha3) + MCha2(Cha4)) + 
     -            MZ2*(12*tmp1 + 
     -               2*(MZ2 - 3*(MCha2(Cha3) + MCha2(Cha4)))))))/
     -     (MZ2**2*SW2**2)

	enddo

	enddo

	do Cha3 = 1,2

	tmp2 = A0(MCha2(Cha3))

	do Neu3 = 1,4

        dup7 = sqrt2*UChaC(Cha3,1)*ZNeu(Neu3,2) + 
     -   UChaC(Cha3,2)*ZNeu(Neu3,3)

        dup8 = (sqrt2*UChaC(Cha3,1)*ZNeu(Neu3,2) + 
     -     UChaC(Cha3,2)*ZNeu(Neu3,3))/sqrt2

        dup9 = 2*UChaC(Cha3,1)*ZNeu(Neu3,2) + 
     -   sqrt2*UChaC(Cha3,2)*ZNeu(Neu3,3)

        dup10 = -(sqrt2*VChaC(Cha3,1)*ZNeu(Neu3,2)) + 
     -   VChaC(Cha3,2)*ZNeu(Neu3,4)

        dup11 = (sqrt2*VChaC(Cha3,1)*ZNeu(Neu3,2) - 
     -     VChaC(Cha3,2)*ZNeu(Neu3,4))/sqrt2

        dup12 = 2*VChaC(Cha3,1)*ZNeu(Neu3,2) - 
     -   sqrt2*VChaC(Cha3,2)*ZNeu(Neu3,4)

        dup13 = sqrt2*UCha(Cha3,1)*ZNeuC(Neu3,2) + 
     -   UCha(Cha3,2)*ZNeuC(Neu3,3)

        dup14 = (sqrt2*UCha(Cha3,1)*ZNeuC(Neu3,2) + 
     -     UCha(Cha3,2)*ZNeuC(Neu3,3))/sqrt2

        dup15 = 2*UCha(Cha3,1)*ZNeuC(Neu3,2) + 
     -   sqrt2*UCha(Cha3,2)*ZNeuC(Neu3,3)

        dup16 = (sqrt2*VCha(Cha3,1)*ZNeuC(Neu3,2) - 
     -     VCha(Cha3,2)*ZNeuC(Neu3,4))/sqrt2

        dr = dr + 1/(72.D0*Pi)*
     -    (Alfa0*(MW2*((dup11*dup16 + dup14*dup8)*
     -             (SW2*(4*MW2 - 3*(MCha2(Cha3) + MNeu2(Neu3))) + 
     -               CW2*(-4*MW2 - 36*tmp2 + 
     -                  12*(MCha2(Cha3) + MNeu2(Neu3)))) + 
     -            3*(2*CW2 + SW2)*A0(MNeu2(Neu3))*
     -             ((dup9*UCha(Cha3,1) + dup12*VCha(Cha3,1))*
     -                ZNeuC(Neu3,2) + 
     -               dup7*UCha(Cha3,2)*ZNeuC(Neu3,3) + 
     -               dup10*VCha(Cha3,2)*ZNeuC(Neu3,4))) + 
     -         3*(CW2 - SW2)*B0(MW2,MCha2(Cha3),MNeu2(Neu3))*
     -          ((2*MW2**2 - MCha2(Cha3)**2 - 
     -               MCha2(Cha3)*(MW2 - 2*MNeu2(Neu3)) - 
     -               MW2*MNeu2(Neu3) - MNeu2(Neu3)**2)*
     -             ((dup9*UCha(Cha3,1) + dup12*VCha(Cha3,1))*
     -                ZNeuC(Neu3,2) + 
     -               dup7*UCha(Cha3,2)*ZNeuC(Neu3,3) + 
     -               dup10*VCha(Cha3,2)*ZNeuC(Neu3,4)) + 
     -            6*MW2*MCha(Cha3)*MNeu(Neu3)*
     -             (dup9*VChaC(Cha3,1)*ZNeu(Neu3,2) - 
     -               dup7*VChaC(Cha3,2)*ZNeu(Neu3,4) + 
     -               dup15*VCha(Cha3,1)*ZNeuC(Neu3,2) - 
     -               dup13*VCha(Cha3,2)*ZNeuC(Neu3,4))) + 
     -         6*B0(0.D0,MCha2(Cha3),MNeu2(Neu3))*
     -          (1/2.D0*((MCha2(Cha3)*
     -                  (-(SW2*(MW2 + MCha2(Cha3))) + 
     -                    CW2*(4*MW2 + MCha2(Cha3))) - 
     -                 2*(SW2*(MW2 - MCha2(Cha3)) + 
     -                    CW2*(2*MW2 + MCha2(Cha3)))*MNeu2(Neu3) + 
     -                 (CW2 - SW2)*MNeu2(Neu3)**2)*
     -               ((dup9*UCha(Cha3,1) + dup12*VCha(Cha3,1))*
     -                  ZNeuC(Neu3,2) + 
     -                 dup7*UCha(Cha3,2)*ZNeuC(Neu3,3) + 
     -                 dup10*VCha(Cha3,2)*ZNeuC(Neu3,4))) + 
     -            3*MW2*SW2*MCha(Cha3)*MNeu(Neu3)*
     -             (dup9*VChaC(Cha3,1)*ZNeu(Neu3,2) - 
     -               dup7*VChaC(Cha3,2)*ZNeu(Neu3,4) + 
     -               dup15*VCha(Cha3,1)*ZNeuC(Neu3,2) - 
     -               dup13*VCha(Cha3,2)*ZNeuC(Neu3,4)))))/
     -     (MW2**2*SW2**2)

	enddo

	enddo

	do Sfe3 = 1,2
	do Cha3 = 1,2

        dr = dr - 1/(16.D0*Pi)*
     -    (Alfa0*UCha(Cha3,1)*UChaC(Cha3,1)*
     -       ((B0(0.D0,MCha2(Cha3),MSf2(Sfe3,2,1)) + 
     -            DB0(0.D0,MCha2(Cha3),MSf2(Sfe3,2,1))*
     -             (MCha2(Cha3) - MSf2(Sfe3,2,1)))*USf(Sfe3,1,2,1)*
     -          USfC(Sfe3,1,2,1) + 
     -         (B0(0.D0,MCha2(Cha3),MSf2(Sfe3,2,2)) + 
     -            DB0(0.D0,MCha2(Cha3),MSf2(Sfe3,2,2))*
     -             (MCha2(Cha3) - MSf2(Sfe3,2,2)))*USf(Sfe3,1,2,2)*
     -          USfC(Sfe3,1,2,2)))/SW2

	enddo
	enddo

	do Cha5 = 1,2

	tmp3 = A0(MCha2(Cha5))

	do Neu5 = 1,4

	dup17 = MSf2(1,1,1) - MSf2(1,1,2)

	dup18 = SW*ZNeu(Neu5,1) - CW*ZNeu(Neu5,2)

        dup19 = (sqrt2*UChaC(Cha5,1)*ZNeu(Neu5,2) + 
     -     UChaC(Cha5,2)*ZNeu(Neu5,3))/sqrt2

        dup20 = (sqrt2*VChaC(Cha5,1)*ZNeu(Neu5,2) - 
     -     VChaC(Cha5,2)*ZNeu(Neu5,4))/sqrt2

	dup21 = SW*ZNeuC(Neu5,1) - CW*ZNeuC(Neu5,2)

        dup22 = (sqrt2*UCha(Cha5,1)*ZNeuC(Neu5,2) + 
     -     UCha(Cha5,2)*ZNeuC(Neu5,3))/sqrt2

        dup23 = (sqrt2*VCha(Cha5,1)*ZNeuC(Neu5,2) - 
     -     VCha(Cha5,2)*ZNeuC(Neu5,4))/sqrt2

        dr = dr + 1/(32.D0*Pi)*
     -    (Alfa0*(VCha(Cha5,1)*
     -          (dup21*sqrt2*VChaC(Cha5,2)*ZNeu(Neu5,4) - 
     -            2*VChaC(Cha5,1)*
     -             (SW*ZNeu(Neu5,1)*ZNeuC(Neu5,2) + 
     -               ZNeu(Neu5,2)*
     -                (SW*ZNeuC(Neu5,1) - 2*CW*ZNeuC(Neu5,2)))) + 
     -         dup18*sqrt2*VCha(Cha5,2)*VChaC(Cha5,1)*
     -          ZNeuC(Neu5,4) + 
     -         (2*(((tmp3*(2*CW2*dup21*
     -                       (dup20*MCha2(Cha5) - 
     -                       2*dup22*MCha(Cha5)*MNeu(Neu5))*
     -                       (MCha2(Cha5) - MSf2(1,1,2))*
     -                       VCha(Cha5,1) + 
     -                       dup18*
     -                       (2*CW2*
     -                       (dup23*MCha2(Cha5) - 
     -                       2*dup19*MCha(Cha5)*MNeu(Neu5))*
     -                       (MCha2(Cha5) - MSf2(1,1,1)) - 
     -                       CW*dup21*MW2*MCha2(Cha5)*VCha(Cha5,1))
     -                       *VChaC(Cha5,1)))/
     -                   ((MCha2(Cha5) - MSf2(1,1,1))*
     -                     (MCha2(Cha5) - MSf2(1,1,2))) + 
     -                  (A0(MNeu2(Neu5))*
     -                     (4*CW2*MCha(Cha5)*MNeu(Neu5)*
     -                       (dup21*dup22*
     -                       (MNeu2(Neu5) - MSf2(1,1,2))*
     -                       VCha(Cha5,1) + 
     -                       dup18*dup19*
     -                       (MNeu2(Neu5) - MSf2(1,1,1))*
     -                       VChaC(Cha5,1)) + 
     -                       MNeu2(Neu5)*
     -                       (CW2*
     -                       (2*dup18*dup23*MSf2(1,1,1)*
     -                       VChaC(Cha5,1) - 
     -                       2*MNeu2(Neu5)*
     -                       (dup20*dup21*VCha(Cha5,1) + 
     -                       dup18*dup23*VChaC(Cha5,1))) + 
     -                       VCha(Cha5,1)*
     -                       (2*CW2*dup20*dup21*MSf2(1,1,2) + 
     -                       dup18*MW2*VChaC(Cha5,1)*
     -                       (CW*SW*ZNeuC(Neu5,1) - 
     -                       CW2*ZNeuC(Neu5,2))))))/
     -                   ((MNeu2(Neu5) - MSf2(1,1,1))*
     -                     (MNeu2(Neu5) - MSf2(1,1,2))))/
     -                (MCha2(Cha5) - MNeu2(Neu5)) + 
     -               (A0(MSf2(1,1,1))*VCha(Cha5,1)*
     -                  (-(MSf2(1,1,1)*VChaC(Cha5,1)*
     -                       (CW2*dup21*
     -                       (MW2 + 2*MSf2(1,1,1) - 2*MSf2(1,1,2))*
     -                       ZNeu(Neu5,2) + 
     -                       MW2*ZNeu(Neu5,1)*
     -                       (-(CW*SW2*ZNeuC(Neu5,1)) + 
     -                       CW2*SW*ZNeuC(Neu5,2)))) + 
     -                    CW2*dup17*dup21*
     -                     (sqrt2*MSf2(1,1,1)*VChaC(Cha5,2)*
     -                       ZNeu(Neu5,4) + 
     -                       2*MCha(Cha5)*MNeu(Neu5)*
     -                       (2*UCha(Cha5,1)*ZNeuC(Neu5,2) + 
     -                       sqrt2*UCha(Cha5,2)*ZNeuC(Neu5,3)))))/
     -                (dup17*(MCha2(Cha5) - MSf2(1,1,1))*
     -                  (-MNeu2(Neu5) + MSf2(1,1,1)))) - 
     -            (2*A0(MSf2(1,1,2))*VChaC(Cha5,1)*
     -               (2*CW2*dup17*dup18*MCha(Cha5)*MNeu(Neu5)*
     -                  (2*UChaC(Cha5,1)*ZNeu(Neu5,2) + 
     -                    sqrt2*UChaC(Cha5,2)*ZNeu(Neu5,3)) + 
     -                 MSf2(1,1,2)*
     -                  (VCha(Cha5,1)*
     -                     (MW2*
     -                       (-(CW*SW2*ZNeu(Neu5,1)) + 
     -                       CW2*SW*ZNeu(Neu5,2))*ZNeuC(Neu5,1) + 
     -                       CW2*
     -                       (MW2 - 2*MSf2(1,1,1) + 2*MSf2(1,1,2))*
     -                       (SW*ZNeu(Neu5,1) - CW*ZNeu(Neu5,2))*
     -                       ZNeuC(Neu5,2)) + 
     -                    CW2*dup17*dup18*sqrt2*VCha(Cha5,2)*
     -                     ZNeuC(Neu5,4))))/
     -             ((MCha2(Cha5) - MSf2(1,1,2))*
     -               (-MNeu2(Neu5) + MSf2(1,1,2))*
     -               (-MSf2(1,1,1) + MSf2(1,1,2))))/CW2))/(CW*SW2)

	enddo

	enddo

	do Sfe3 = 1,2
	do Gen3 = 1,3

        dr = dr + 1/(144.D0*Pi)*
     -    (Alfa0*((6*(CW2 - SW2)*
     -             B0(MW2,MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3))*
     -             (MSf2(1,1,Gen3)**2 + 
     -               (MW2 - MSf2(Sfe3,2,Gen3))**2 - 
     -               2*MSf2(1,1,Gen3)*(MW2 + MSf2(Sfe3,2,Gen3))) - 
     -            6*B0(0.D0,MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3))*
     -             (MSf2(1,1,Gen3)*
     -                (SW2*(2*MW2 - MSf2(1,1,Gen3)) + 
     -                  CW2*(MW2 + MSf2(1,1,Gen3))) + 
     -               (CW2 - SW2)*
     -                (-((MW2 + 2*MSf2(1,1,Gen3))*
     -                     MSf2(Sfe3,2,Gen3)) + 
     -                  MSf2(Sfe3,2,Gen3)**2)))*USf(Sfe3,1,2,Gen3)*
     -          USfC(Sfe3,1,2,Gen3) + 
     -         MW2*(-32*MW2*SW2**2 + 
     -            (CW2*(4*MW2 - 
     -                  12*(MSf2(1,1,Gen3) + MSf2(Sfe3,2,Gen3))) + 
     -               SW2*(-4*MW2 + 
     -                  3*(MSf2(1,1,Gen3) + MSf2(Sfe3,2,Gen3))))*
     -             USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     -            ((6*A0(MSf2(Sfe3,2,Gen3))*
     -                  (2*MW2*MZ2*SW2**2 - 
     -                    MSf2(Sfe3,2,Gen3)*
     -                     ((3*MW2*(1 - 2*SW2)**2 + 
     -                       MZ2*(-CW2 + SW2))*USf(Sfe3,1,2,Gen3)*
     -                       USfC(Sfe3,1,2,Gen3) + 
     -                       12*MW2*SW2**2*USf(Sfe3,2,2,Gen3)*
     -                       USfC(Sfe3,2,2,Gen3))))/
     -                MSf2(Sfe3,2,Gen3) + 
     -               2*((A0(MSf2(Sfe3,3,Gen3))*
     -                     (8*MW2*MZ2*SW2**2 - 
     -                       3*MSf2(Sfe3,3,Gen3)*
     -                       ((-9*CW2*MZ2 + MW2*(3 - 4*SW2)**2)*
     -                       USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3)
     -                        + 16*MW2*SW2**2*USf(Sfe3,2,3,Gen3)*
     -                       USfC(Sfe3,2,3,Gen3))))/
     -                   MSf2(Sfe3,3,Gen3) + 
     -                  (A0(MSf2(Sfe3,4,Gen3))*
     -                     (2*MW2*MZ2*SW2**2 - 
     -                       3*MSf2(Sfe3,4,Gen3)*
     -                       ((-9*CW2*MZ2 + MW2*(3 - 2*SW2)**2)*
     -                       USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3)
     -                        + 4*MW2*SW2**2*USf(Sfe3,2,4,Gen3)*
     -                       USfC(Sfe3,2,4,Gen3))))/
     -                   MSf2(Sfe3,4,Gen3)))/MZ2)))/(MW2**2*SW2**2)

	enddo
	enddo

	do Neu4 = 1,4

	tmp4 = A0(MNeu2(Neu4))

	do Neu3 = 1,4

	dup24 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - ZNeu(Neu4,4)*ZNeuC(Neu3,4)

	dup25 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - ZNeu(Neu3,4)*ZNeuC(Neu4,4)

        dr = dr + 1/(144.D0*Pi)*
     -    (Alfa0*(dup24*dup25*
     -          (2*MZ2**2 + 
     -            6*MNeu2(Neu3)*
     -             (MZ2*(-1 + B0(0.D0,MNeu2(Neu3),MNeu2(Neu4))) + 
     -               B0(0.D0,MNeu2(Neu3),MNeu2(Neu4))*MNeu2(Neu4)) + 
     -            MZ2*(12*tmp4 - 
     -               6*(1 + B0(0.D0,MNeu2(Neu3),MNeu2(Neu4)))*
     -                MNeu2(Neu4)) - 
     -            3*B0(0.D0,MNeu2(Neu3),MNeu2(Neu4))*
     -             (MNeu2(Neu3)**2 + MNeu2(Neu4)**2)) + 
     -         3*B0(MZ2,MNeu2(Neu3),MNeu2(Neu4))*
     -          (3*(dup24**2 + dup25**2)*MZ2*MNeu(Neu3)*
     -             MNeu(Neu4) + 
     -            dup24*dup25*
     -             (-2*MZ2**2 + MNeu2(Neu3)**2 + 
     -               MNeu2(Neu3)*(MZ2 - 2*MNeu2(Neu4)) + 
     -               MZ2*MNeu2(Neu4) + MNeu2(Neu4)**2))))/
     -     (MZ2**2*SW2**2)

	enddo

	enddo

	do Sfe3 = 1,2
	do Neu3 = 1,4

        dr = dr - 1/(32.D0*Pi)*
     -    (Alfa0*((B0(0.D0,MNeu2(Neu3),MSf2(Sfe3,2,1)) + 
     -            DB0(0.D0,MNeu2(Neu3),MSf2(Sfe3,2,1))*
     -             (MNeu2(Neu3) - MSf2(Sfe3,2,1)))*USf(Sfe3,1,2,1)*
     -          USfC(Sfe3,1,2,1) + 
     -         (B0(0.D0,MNeu2(Neu3),MSf2(Sfe3,2,2)) + 
     -            DB0(0.D0,MNeu2(Neu3),MSf2(Sfe3,2,2))*
     -             (MNeu2(Neu3) - MSf2(Sfe3,2,2)))*USf(Sfe3,1,2,2)*
     -          USfC(Sfe3,1,2,2))*
     -       (SW*ZNeu(Neu3,1) + CW*ZNeu(Neu3,2))*
     -       (SW*ZNeuC(Neu3,1) + CW*ZNeuC(Neu3,2)))/(CW2*SW2)

	enddo
	enddo

	do Sfe5 = 1,2

	tmp5 = A0(MSf2(Sfe5,2,1))

	tmp6 = A0(MSf2(Sfe5,2,2))

	do Neu5 = 1,4

	dup26 = SW*ZNeu(Neu5,1) - CW*ZNeu(Neu5,2)

	dup27 = SW*ZNeu(Neu5,1) + CW*ZNeu(Neu5,2)

	dup28 = SW*ZNeuC(Neu5,1) - CW*ZNeuC(Neu5,2)

	dup29 = SW*ZNeuC(Neu5,1) + CW*ZNeuC(Neu5,2)

        dr = dr + 1/(32.D0*Pi)*
     -    (Alfa0*(2*(dup27*dup28*
     -             ((A0(MSf2(1,1,1))*MSf2(1,1,1))/
     -                ((MNeu2(Neu5) - MSf2(1,1,1))*
     -                  (MSf2(1,1,1) - MSf2(Sfe5,2,1))) + 
     -               (tmp5*MSf2(Sfe5,2,1))/
     -                ((MNeu2(Neu5) - MSf2(Sfe5,2,1))*
     -                  (-MSf2(1,1,1) + MSf2(Sfe5,2,1))))*
     -             USf(Sfe5,1,2,1)*USfC(Sfe5,1,2,1) + 
     -            (dup26*dup29*tmp6*MSf2(Sfe5,2,2)*USf(Sfe5,1,2,2)*
     -               USfC(Sfe5,1,2,2))/
     -             ((MNeu2(Neu5) - MSf2(Sfe5,2,2))*
     -               (-MSf2(1,1,2) + MSf2(Sfe5,2,2)))) - 
     -         2*((dup26*dup29*A0(MSf2(1,1,2))*MSf2(1,1,2)*
     -               USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2))/
     -             ((-MNeu2(Neu5) + MSf2(1,1,2))*
     -               (MSf2(1,1,2) - MSf2(Sfe5,2,2))) + 
     -            A0(MNeu2(Neu5))*MNeu2(Neu5)*
     -             ((dup27*dup28*USf(Sfe5,1,2,1)*USfC(Sfe5,1,2,1))/
     -                ((MNeu2(Neu5) - MSf2(1,1,1))*
     -                  (MNeu2(Neu5) - MSf2(Sfe5,2,1))) + 
     -               (dup26*dup29*USf(Sfe5,1,2,2)*
     -                  USfC(Sfe5,1,2,2))/
     -                ((MNeu2(Neu5) - MSf2(1,1,2))*
     -                  (MNeu2(Neu5) - MSf2(Sfe5,2,2))))) + 
     -         USf(Sfe5,1,2,1)*USfC(Sfe5,1,2,1)*
     -          (-((SW2*ZNeu(Neu5,1) + CW*SW*ZNeu(Neu5,2))*
     -               ZNeuC(Neu5,1)) + 
     -            (CW*SW*ZNeu(Neu5,1) + CW2*ZNeu(Neu5,2))*
     -             ZNeuC(Neu5,2)) + 
     -         USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     -          (ZNeu(Neu5,2)*
     -             (CW*SW*ZNeuC(Neu5,1) + CW2*ZNeuC(Neu5,2)) - 
     -            ZNeu(Neu5,1)*
     -             (SW2*ZNeuC(Neu5,1) + CW*SW*ZNeuC(Neu5,2)))))/
     -     (CW2*SW2)

	enddo

	enddo

	do Sfe5 = 1,2

	tmp8 = A0(MSf2(Sfe5,2,1))

	tmp9 = A0(MSf2(Sfe5,2,2))

	do Cha5 = 1,2

	tmp7 = A0(MCha2(Cha5))

	do Neu5 = 1,4

	dup30 = SW*ZNeu(Neu5,1) + CW*ZNeu(Neu5,2)

        dup31 = 2*UChaC(Cha5,1)*ZNeu(Neu5,2) + 
     -   sqrt2*UChaC(Cha5,2)*ZNeu(Neu5,3)

	dup32 = SW*ZNeuC(Neu5,1) + CW*ZNeuC(Neu5,2)

        dup33 = (sqrt2*UCha(Cha5,1)*ZNeuC(Neu5,2) + 
     -     UCha(Cha5,2)*ZNeuC(Neu5,3))/sqrt2

        dup34 = 2*VCha(Cha5,1)*ZNeuC(Neu5,2) - 
     -   sqrt2*VCha(Cha5,2)*ZNeuC(Neu5,4)

        dr = dr + 1/(32.D0*Pi)*
     -    (Alfa0*((-4*(A0(MNeu2(Neu5))*
     -               (MNeu2(Neu5)*
     -                  ((dup30*dup33*UChaC(Cha5,1)*
     -                       USf(Sfe5,1,2,1)*USfC(Sfe5,1,2,1))/
     -                     (MNeu2(Neu5) - MSf2(Sfe5,2,1)) + 
     -                    (dup31*dup32*UCha(Cha5,1)*
     -                       USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2))/
     -                     (2*MNeu2(Neu5) - 2*MSf2(Sfe5,2,2))) + 
     -                 MCha(Cha5)*MNeu(Neu5)*
     -                  ((dup32*dup34*UCha(Cha5,1)*USf(Sfe5,1,2,2)*
     -                       USfC(Sfe5,1,2,2))/
     -                     (-MNeu2(Neu5) + MSf2(Sfe5,2,2)) + 
     -                    (2*dup30*UChaC(Cha5,1)*USf(Sfe5,1,2,1)*
     -                       USfC(Sfe5,1,2,1)*
     -                       (-(sqrt2*VChaC(Cha5,1)*
     -                      ZNeu(Neu5,2)) + 
     -                       VChaC(Cha5,2)*ZNeu(Neu5,4)))/
     -                     (sqrt2*(MNeu2(Neu5) - MSf2(Sfe5,2,1)))))
     -                + tmp7*
     -               (MCha2(Cha5)*
     -                  (-(1/2.D0*
     -                       (dup31*dup32*UCha(Cha5,1)*
     -                       USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2))/
     -                       (MCha2(Cha5) - MSf2(Sfe5,2,2))) + 
     -                    (dup30*dup33*UChaC(Cha5,1)*
     -                       USf(Sfe5,1,2,1)*USfC(Sfe5,1,2,1))/
     -                     (-MCha2(Cha5) + MSf2(Sfe5,2,1))) + 
     -                 MCha(Cha5)*MNeu(Neu5)*
     -                  ((dup32*dup34*UCha(Cha5,1)*USf(Sfe5,1,2,2)*
     -                       USfC(Sfe5,1,2,2))/
     -                     (MCha2(Cha5) - MSf2(Sfe5,2,2)) + 
     -                    (dup30*UChaC(Cha5,1)*USf(Sfe5,1,2,1)*
     -                       USfC(Sfe5,1,2,1)*
     -                       (-2*VChaC(Cha5,1)*ZNeu(Neu5,2) + 
     -                       sqrt2*VChaC(Cha5,2)*ZNeu(Neu5,4)))/
     -                     (-MCha2(Cha5) + MSf2(Sfe5,2,1))))))/
     -          (-MCha2(Cha5) + MNeu2(Neu5)) + 
     -         dup30*UChaC(Cha5,1)*USf(Sfe5,1,2,1)*
     -          USfC(Sfe5,1,2,1)*
     -          (sqrt2*UCha(Cha5,2)*ZNeuC(Neu5,3) + 
     -            (2*tmp8*(2*MCha(Cha5)*MNeu(Neu5)*
     -                  (-2*VChaC(Cha5,1)*ZNeu(Neu5,2) + 
     -                    sqrt2*VChaC(Cha5,2)*ZNeu(Neu5,4)) + 
     -                 MSf2(Sfe5,2,1)*
     -                  (2*UCha(Cha5,1)*ZNeuC(Neu5,2) + 
     -                    sqrt2*UCha(Cha5,2)*ZNeuC(Neu5,3))))/
     -             ((MCha2(Cha5) - MSf2(Sfe5,2,1))*
     -               (-MNeu2(Neu5) + MSf2(Sfe5,2,1)))) + 
     -         UCha(Cha5,1)*
     -          (2*dup30*UChaC(Cha5,1)*USf(Sfe5,1,2,1)*
     -             USfC(Sfe5,1,2,1)*ZNeuC(Neu5,2) + 
     -            dup32*USf(Sfe5,1,2,2)*USfC(Sfe5,1,2,2)*
     -             (sqrt2*UChaC(Cha5,2)*ZNeu(Neu5,3) + 
     -               2*(UChaC(Cha5,1)*ZNeu(Neu5,2) + 
     -                  (tmp9*
     -                     (dup31*MSf2(Sfe5,2,2) + 
     -                       2*MCha(Cha5)*MNeu(Neu5)*
     -                       (-2*VCha(Cha5,1)*ZNeuC(Neu5,2) + 
     -                       sqrt2*VCha(Cha5,2)*ZNeuC(Neu5,4))))/
     -                   ((MCha2(Cha5) - MSf2(Sfe5,2,2))*
     -                     (-MNeu2(Neu5) + MSf2(Sfe5,2,2))))))))/
     -     (CW*SW2)

	enddo

	enddo

	enddo

	do Sfe4 = 1,2
	do Gen3 = 1,3

	tmp10 = A0(MSf2(Sfe4,2,Gen3))

	tmp11 = A0(MSf2(Sfe4,3,Gen3))

	tmp12 = A0(MSf2(Sfe4,4,Gen3))

	do Sfe3 = 1,2

        dup35 = -(USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     -   2*SW2*(USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     -      USf(Sfe4,2,2,Gen3)*USfC(Sfe3,2,2,Gen3))

        dup36 = -3*USf(Sfe4,1,3,Gen3)*USfC(Sfe3,1,3,Gen3) + 
     -   4*SW2*(USf(Sfe4,1,3,Gen3)*USfC(Sfe3,1,3,Gen3) + 
     -      USf(Sfe4,2,3,Gen3)*USfC(Sfe3,2,3,Gen3))

        dup37 = -3*USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,4,Gen3) + 
     -   2*SW2*(USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,4,Gen3) + 
     -      USf(Sfe4,2,4,Gen3)*USfC(Sfe3,2,4,Gen3))

        dup38 = -(USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3)) + 
     -   2*SW2*(USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     -      USf(Sfe3,2,2,Gen3)*USfC(Sfe4,2,2,Gen3))

        dup39 = -3*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,1,3,Gen3) + 
     -   4*SW2*(USf(Sfe3,1,3,Gen3)*USfC(Sfe4,1,3,Gen3) + 
     -      USf(Sfe3,2,3,Gen3)*USfC(Sfe4,2,3,Gen3))

        dup40 = -3*USf(Sfe3,1,4,Gen3)*USfC(Sfe4,1,4,Gen3) + 
     -   2*SW2*(USf(Sfe3,1,4,Gen3)*USfC(Sfe4,1,4,Gen3) + 
     -      USf(Sfe3,2,4,Gen3)*USfC(Sfe4,2,4,Gen3))

        dr = dr + 1/(432.D0*Pi)*
     -    (Alfa0*((dup35*dup38*
     -             (9*B0(0.D0,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     -                (MSf2(Sfe3,2,Gen3) - MSf2(Sfe4,2,Gen3))*
     -                (MZ2 + MSf2(Sfe3,2,Gen3) - MSf2(Sfe4,2,Gen3))
     -                 - 9*
     -                B0(MZ2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     -                (MSf2(Sfe3,2,Gen3)**2 + 
     -                  (MZ2 - MSf2(Sfe4,2,Gen3))**2 - 
     -                  2*MSf2(Sfe3,2,Gen3)*
     -                   (MZ2 + MSf2(Sfe4,2,Gen3)))) + 
     -            dup36*dup39*
     -             (3*B0(0.D0,MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3))*
     -                (MSf2(Sfe3,3,Gen3) - MSf2(Sfe4,3,Gen3))*
     -                (MZ2 + MSf2(Sfe3,3,Gen3) - MSf2(Sfe4,3,Gen3))
     -                 - 3*
     -                B0(MZ2,MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3))*
     -                (MSf2(Sfe3,3,Gen3)**2 + 
     -                  (MZ2 - MSf2(Sfe4,3,Gen3))**2 - 
     -                  2*MSf2(Sfe3,3,Gen3)*
     -                   (MZ2 + MSf2(Sfe4,3,Gen3)))) + 
     -            dup37*dup40*
     -             (3*B0(0.D0,MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3))*
     -                (MSf2(Sfe3,4,Gen3) - MSf2(Sfe4,4,Gen3))*
     -                (MZ2 + MSf2(Sfe3,4,Gen3) - MSf2(Sfe4,4,Gen3))
     -                 - 3*
     -                B0(MZ2,MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3))*
     -                (MSf2(Sfe3,4,Gen3)**2 + 
     -                  (MZ2 - MSf2(Sfe4,4,Gen3))**2 - 
     -                  2*MSf2(Sfe3,4,Gen3)*
     -                   (MZ2 + MSf2(Sfe4,4,Gen3)))))/MZ2**2 + 
     -         ((54*(CW2 - SW2)*
     -               B0(MW2,MSf2(Sfe3,3,Gen3),MSf2(Sfe4,4,Gen3))*
     -               (MSf2(Sfe3,3,Gen3)**2 + 
     -                 (MW2 - MSf2(Sfe4,4,Gen3))**2 - 
     -                 2*MSf2(Sfe3,3,Gen3)*
     -                  (MW2 + MSf2(Sfe4,4,Gen3))) - 
     -              54*B0(0.D0,MSf2(Sfe3,3,Gen3),MSf2(Sfe4,4,Gen3))*
     -               (MSf2(Sfe3,3,Gen3)*
     -                  (SW2*(2*MW2 - MSf2(Sfe3,3,Gen3)) + 
     -                    CW2*(MW2 + MSf2(Sfe3,3,Gen3))) + 
     -                 (CW2 - SW2)*
     -                  (-((MW2 + 2*MSf2(Sfe3,3,Gen3))*
     -                       MSf2(Sfe4,4,Gen3)) + 
     -                    MSf2(Sfe4,4,Gen3)**2)))*
     -            USf(Sfe3,1,3,Gen3)*USf(Sfe4,1,4,Gen3)*
     -            USfC(Sfe3,1,3,Gen3)*USfC(Sfe4,1,4,Gen3))/MW2**2+
     -           (18*dup35*dup38*tmp10 + 
     -            (2*MW2*((-(dup36*dup39) - dup37*dup40)*MZ2 + 
     -                  dup35*dup38*
     -                   (-3*MZ2 + 
     -                     9*
     -                      (MSf2(Sfe3,2,Gen3) + MSf2(Sfe4,2,Gen3))
     -                     ) + 
     -                  3*(dup36*dup39*
     -                      (MSf2(Sfe3,3,Gen3) + MSf2(Sfe4,3,Gen3))
     -                       + dup37*dup40*
     -                      (MSf2(Sfe3,4,Gen3) + MSf2(Sfe4,4,Gen3))
     -                     )) + 
     -               9*MZ2*
     -                (CW2*
     -                   (4*MW2 - 
     -                     12*
     -                      (MSf2(Sfe3,3,Gen3) + MSf2(Sfe4,4,Gen3))
     -                     ) + 
     -                  SW2*
     -                   (-4*MW2 + 
     -                     3*
     -                      (MSf2(Sfe3,3,Gen3) + MSf2(Sfe4,4,Gen3))
     -                     ))*USf(Sfe3,1,3,Gen3)*
     -                USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,3,Gen3)*
     -                USfC(Sfe4,1,4,Gen3))/MW2 + 
     -            6*(dup36*dup39*tmp11 + 
     -               (tmp12*
     -                  (dup37*dup40*MW2 + 
     -                    MZ2*(-18*CW2 - 9*SW2)*USf(Sfe3,1,3,Gen3)*
     -                     USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,3,Gen3)*
     -                     USfC(Sfe4,1,4,Gen3)))/MW2))/MZ2))/SW2**2

	enddo

	enddo
	enddo

	do Cha5 = 1,2

	tmp15 = A0(MCha2(Cha5))

	do Neu5 = 1,4

	tmp14 = A0(MNeu2(Neu5))

	do Sfe5 = 1,2

	tmp13 = A0(MSf2(Sfe5,2,2))

	do Sfe6 = 1,2

        dr = dr + 1/(16.D0*Pi)*
     -    (Alfa0*MW2*((-((tmp15*MCha2(Cha5))/
     -               ((MCha2(Cha5) - MSf2(Sfe5,2,2))*
     -                 (MCha2(Cha5) - MSf2(Sfe6,2,1)))) + 
     -            (tmp14*MNeu2(Neu5))/
     -             ((MNeu2(Neu5) - MSf2(Sfe5,2,2))*
     -               (MNeu2(Neu5) - MSf2(Sfe6,2,1))))/
     -          (MCha2(Cha5) - MNeu2(Neu5)) + 
     -         ((tmp13*MSf2(Sfe5,2,2))/
     -             ((MCha2(Cha5) - MSf2(Sfe5,2,2))*
     -               (-MNeu2(Neu5) + MSf2(Sfe5,2,2))) + 
     -            (A0(MSf2(Sfe6,2,1))*MSf2(Sfe6,2,1))/
     -             ((MCha2(Cha5) - MSf2(Sfe6,2,1))*
     -               (MNeu2(Neu5) - MSf2(Sfe6,2,1))))/
     -          (MSf2(Sfe5,2,2) - MSf2(Sfe6,2,1)))*UCha(Cha5,1)*
     -       UChaC(Cha5,1)*USf(Sfe5,1,2,2)*USf(Sfe6,1,2,1)*
     -       USfC(Sfe5,1,2,2)*USfC(Sfe6,1,2,1)*
     -       (SW*ZNeu(Neu5,1) + CW*ZNeu(Neu5,2))*
     -       (SW*ZNeuC(Neu5,1) + CW*ZNeuC(Neu5,2)))/(CW2*SW2)

	enddo

	enddo

	enddo

	enddo

#ifdef DETAILED_DEBUG
	DPOPE 'dr =', dr ENDL
#endif

	end


