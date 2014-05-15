* DeltaR1L_sm.F
* generated 24 Jun 2010 11:55
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"


	subroutine DeltaR1L_sm(dr)
	implicit none
	double complex dr

#include "FH.h"
#include "looptools.h"

	double complex dup1

	dup1 = MH2*(MH2 - 4*MW2) + 12*MW2**2

        dr = dr + 1/(288.D0*Pi)*
     -    (Alfa0*(((6*((((MH2 - MW2)*
     -                       (-3*MW2*MZ2*(MW2 - CW2*MZ2) + 
     -                       MH2*(MW2**2 - CW2*MZ2**2)) + 
     -                      dup1*MZ2**2*SW2)*A0(MH2))/(MH2 - MW2)+
     -                   ((-(CW2*(MW2 - MZ2)*
     -                       (MH2*MW2**2 + 
     -                       MZ2*
     -                       (CW2*(1 + 8*CW2)*MZ2**2 + 
     -                       MW2*(MW2 + (-3 - 24*CW2)*CW2*MZ2))))+
     -                        MZ2*
     -                       (((-1 - 8*CW2)*CW2*MZ2**3 + 
     -                       MW2*
     -                       (CW2*(4 + 44*CW2)*MZ2**2 + 
     -                       MW2*(-30*MW2 + CW2*(72 + 12*CW2)*MZ2))
     -                       )*SW2 + 
     -                       12*MW2**2*(5*MW2 - 6*CW2*MZ2)*
     -                      SW2**2 - 24*MW2**3*SW2**3))*A0(MZ2))/
     -                  (CW2*(MW2 - MZ2))))/MZ2**2 - 
     -            (6*(CW2 - SW2)*
     -               (CW2*(4*MW2*MZ2 - MZ2**2 + 
     -                    CW2*(-8*MZ2**2 + MW2*(60*MW2 + 44*MZ2)))-
     -                   12*MW2**2*SW2**2)*B0(MW2,MW2,MZ2))/CW2)/
     -          MW2**2 + (((CW2*
     -                   (4*
     -                      (54*CW2**2*MW2**2 + 
     -                       (2 + CW2*(-40 + 11*CW2))*MW2*MZ2 - 
     -                       3*CW2*(1 + 8*CW2)*MZ2**2 + 
     -                       3*MH2*(MW2 - CW2*MZ2)) + 
     -                     (MZ2*(178*MW2 + 3*(MH2 + MZ2)) - 
     -                       8*CW2*
     -                       (-3*MZ2**2 + MW2*(6*MW2 + 26*MZ2)))*
     -                      SW2) + 
     -                  MW2*
     -                   (24*CW2*(13*MW2 + 25*MZ2)*SW2**2 + 
     -                     36*(8*MW2 + 3*MZ2)*SW2**3))/MW2 + 
     -               (6*(CW2*SW2*
     -                     (MZ2*
     -                       (MW2*
     -                       (-74*MW2**2 + MH2*(MH2 + 82*MW2)) - 
     -                       (2*MW2**2 + MH2*(MH2 + 6*MW2))*MZ2 + 
     -                       (-MH2 + MW2)*MZ2**2) + 
     -                       4*(MH2 - MW2)*MW2*
     -                       (MW2*(13*MW2 - 11*MZ2) - 20*MZ2**2)*
     -                       SW2) + 
     -                    CW2**2*(-MH2 + MW2)*
     -                     ((MW2 - MZ2)*MZ2*(MH2 + 38*MW2 + MZ2) + 
     -                       8*
     -                       (MZ2**3 + 
     -                       MW2*(MW2*(MW2 - 13*MZ2) + 5*MZ2**2))*
     -                       SW2) + 
     -                    (MH2 - MW2)*
     -                     (4*CW2**3*(9*MW2 - 2*MZ2)*(MW2 - MZ2)*
     -                       (MW2 + MZ2) + 
     -                       6*MW2**2*SW2*
     -                       (8*MW2*SW2**2 + 
     -                       MZ2*(-5 + (10 - 12*SW2)*SW2))))*
     -                  A0(MW2))/(MW2**2*(-MH2 + MW2)*(MW2 - MZ2)))
     -              /MZ2 - 
     -            (6*(CW2*MH2*(MH2 - 4*MZ2) + 12*MW2*MZ2)*
     -               B0(MZ2,MH2,MZ2))/MZ2**2)/CW2 + 
     -         6*((dup1*(CW2 - SW2)*B0(MW2,MH2,MW2))/MW2**2 + 
     -            (((-20*MW2 - MZ2)*SW2**2 + 
     -                 CW2*
     -                  (CW2*(60*MW2 + 39*MZ2) + 
     -                    (-8*MW2 + 2*MZ2)*SW2))*B0(MZ2,MW2,MW2))/
     -             MZ2)))/SW2**2

#ifdef DETAILED_DEBUG
	DPOPE 'dr =', dr ENDL
#endif

	end


