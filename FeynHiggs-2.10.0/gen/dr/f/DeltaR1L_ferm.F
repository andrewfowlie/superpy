* DeltaR1L_ferm.F
* generated 24 Jun 2010 11:55
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"


	subroutine DeltaR1L_ferm(dr)
	implicit none
	double complex dr

#include "FH.h"
#include "looptools.h"

	integer Gen3

        dr = dr - 1/(216.D0*Pi)*
     -    (Alfa0*(-96*(2*MT2 + MZ2)*SW2**2*A0(MT2) + 
     -         MT2*(192*MT2*SW2**2*(1 + B0(MZ2,MT2,MT2)) + 
     -            MZ2*(-9 + 27*B0(MZ2,0.D0,0.D0) + 
     -               SW2**2*(64 + 96*B0(MZ2,MT2,MT2))))))/
     -     (MT2*MZ2*SW2**2)

	do Gen3 = 1,3

        dr = dr + 1/(144.D0*Pi)*
     -    (Alfa0*(6*(((CW2 - SW2)*(2*MW2 + Mf2(2,Gen3))*
     -               (B0(MW2,0.D0,Mf2(2,Gen3))*(MW2 - Mf2(2,Gen3)) + 
     -                 B0(0.D0,0.D0,Mf2(2,Gen3))*Mf2(2,Gen3)))/MW2**2+
     -              (B0(MZ2,Mf2(3,Gen3),Mf2(3,Gen3))*
     -               (MZ2*(-3 + 8*SW2) + (3 + 16*SW2)*Mf2(3,Gen3)))
     -              /MZ2) + 
     -         (12*(3 - 8*SW2)*A0(Mf2(3,Gen3)) + 
     -            6*(((MW2*(2 - 8*SW2) + MZ2*(-4*CW2 + SW2))*
     -                   A0(Mf2(2,Gen3)) + 
     -                  (MW2*(6 - 8*SW2) + MZ2*(-12*CW2 + 3*SW2))*
     -                   A0(Mf2(4,Gen3)))/MW2 + 
     -               B0(MZ2,Mf2(2,Gen3),Mf2(2,Gen3))*
     -                (MZ2*(-1 + 4*SW2) + (1 + 8*SW2)*Mf2(2,Gen3))+
     -                 B0(MZ2,Mf2(4,Gen3),Mf2(4,Gen3))*
     -                (MZ2*(-3 + 4*SW2) + (3 + 8*SW2)*Mf2(4,Gen3)))
     -              + (-2*MW2*MZ2*(-7 + 8*(CW2 + SW2)) + 
     -               (MZ2*(36*CW2 - 9*SW2) + MW2*(-36 + 96*SW2))*
     -                Mf2(3,Gen3) + 
     -               3*((MZ2*(4*CW2 - SW2) + MW2*(-4 + 16*SW2))*
     -                   Mf2(2,Gen3) + 
     -                  (MZ2*(12*CW2 - 3*SW2) + 
     -                     MW2*(-12 + 16*SW2))*Mf2(4,Gen3)))/MW2)/
     -          MZ2 + (-18*(CW2 - SW2)*
     -             B0(MW2,Mf2(3,Gen3),Mf2(4,Gen3))*
     -             (-2*MW2**2 + Mf2(3,Gen3)**2 + 
     -               Mf2(3,Gen3)*(MW2 - 2*Mf2(4,Gen3)) + 
     -               Mf2(4,Gen3)*(MW2 + Mf2(4,Gen3))) + 
     -            18*B0(0.D0,Mf2(3,Gen3),Mf2(4,Gen3))*
     -             (-(Mf2(3,Gen3)*
     -                  (SW2*(MW2 - 2*Mf2(4,Gen3)) + 
     -                    2*CW2*(MW2 + Mf2(4,Gen3)))) + 
     -               (CW2 - SW2)*
     -                (Mf2(3,Gen3)**2 + 
     -                  Mf2(4,Gen3)*(2*MW2 + Mf2(4,Gen3)))))/MW2**2
     -         ))/SW2**2

	enddo

#ifdef DETAILED_DEBUG
	DPOPE 'dr =', dr ENDL
#endif

	end


