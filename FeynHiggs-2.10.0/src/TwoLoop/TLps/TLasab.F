* TLasab.F
* Two-loop O(a_s a_b) corrections to the CP-even Higgs mass matrix.
* Routine written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
* Based on A. Brignole, G. Degrassi, P. Slavich, and F. Zwirner,
* Nucl. Phys. B643 (2002) 79 [hep-ph/0206101].
* this file is part of FeynHiggs
* last modified 4 Dec 13 th

* Revision history:
* 14/07/2004: rewritten completely by TH for FeynHiggs
* 03/03/2004: mg is given as input instead of mg^2;
*             value of pi corrected (10th digit);
*             unused variables cleaned up.
* 27/03/2003: problem fixed for s2b=0.

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine TLasab(dPhi)
	implicit none
	RealType dPhi(3)

#include "TLps.h"

	RealType x1, x3, x4, x6, x7, x5, x8, x9, GG, FF, FFo
	RealType dMSb12, dMSb22, ds2tb, dhbhb, dAb, c2, k

	x1 = Li2omx(MSb12/pMGl2)

	x3 = Li2omx(MSb22/pMGl2)

	x4 = Log(pMGl2)

	x6 = Log(MSb12)

	x7 = Log(MSb22)

	x5 = Log(Q)

	x8 = Log(Abs(1 - MSb12/pMGl2))

	x9 = Log(Abs(1 - MSb22/pMGl2))

        GG = (-4*(4*MSb2diff - 2*pMGl2*x1 + MSb2sum*x1 + 
     -       2*pMGl2*x3 - MSb2sum*x3 - 2*MSb12*x5 + 2*MSb22*x5 + 
     -       2*MSb12*x6 + MSb2sum*x6 - MSb2sum*x4*x6 + 
     -       MSb2sum*x5*x6 - 2*MSb22*x7 - MSb2sum*x7 + 
     -       MSb2sum*x4*x7 - MSb2sum*x5*x7))/MSb2diff

        FF = (MSb12**3*MSb2diff*s2tb2*(x5 - x6) + 
     -     MSb2diff*(4*pMGl2*MSb2diff*MSb2sum*(1 - x4 + x5) - 
     -        MSb22**2*s2tb2*(MSb2diff + MSb22*(x5 - x7))) - 
     -     MSb12**2*(MSb2diff**2*s2tb2 + 
     -        (3 + 13*c2tb2)*MSb22*MSb2diff*(x5 - x6) + 
     -        6*c2tb2*MSb22*MSb2sum*(x5 - x6)*(x6 - x7)) + 
     -     MSb12*MSb22*(2*MSb2diff**2*
     -         (5 + 7*c2tb2 + s2tb2*(2*x5 - x6 - x7)) + 
     -        c2tb2*MSb2sum*
     -         (6*MSb22*(x5 - x7) + MSb2sum*(x6 - x7))*(x6 - x7) + 
     -        MSb2diff*(8*pMGl2*
     -            (x1 - x3 - (-1 + x4 - x5)*(x6 - x7)) + 
     -           (3 + 13*c2tb2)*MSb22*(x5 - x7) + 
     -           MSb2sum*(-4*x1 + 4*x3 - 
     -              (x6 - x7)*
     -               (-3 - 9*c2tb2 - 4*x4 - 2*x5 + 3*x6 + 3*x7)))))
     -    /(MSb12*MSb22*MSb2diff**2)

* COUNTER TERMS:

	if( onshell ) then
	  FFo = 2 + (MSb2sum*(x6 - x7))/MSb2diff

	  dMSb12 = s2tb2*(MSb2diff + MSb22*x5 - MSb22*x7) - 
     -     (2*pMGl2**2*x8)/MSb12 - 
     -     MSb12*(2 + 2*x4 + s2tb2*x5 - 3*x6 + c2tb2*x6 + 2*x8) + 
     -     pMGl2*(-6 + 4*x4 - 4*x5 + 4*x8)

	  dMSb22 = (-2*pMGl2**2*x9 + 
     -     2*pMGl2*MSb22*(-3 + 2*x4 - 2*x5 + 2*x9) - 
     -     MSb22*(s2tb2*(MSb2diff - MSb12*x5 + MSb12*x6) + 
     -        MSb22*(2 + 2*x4 + s2tb2*x5 - 3*x7 + c2tb2*x7 + 2*x9))
     -     )/MSb22

	  ds2tb = (-2*c2tb2*s2tb*
     -     (MSb2diff + MSb12*(-x5 + x6) + MSb22*(x5 - x7)))/
     -     MSb2diff

	  dhbhb = -4 + (4*MSb12*(x5 - x6) + 4*MSb22*(-x5 + x7) + 
     -      (2*(pMGl2 - MSb12)**2*x8)/MSb12)/MSb2diff - 
     -     2*(x4 - x5 + ((pMGl2 - MSb22)**2*x9)/(MSb22*MSb2diff))

	  c2 = 4 - 2*x4 + 2*x5 - x8 + (pMGl2*x8)/MSb12 - x9 + 
     -     (pMGl2*x9)/MSb22

	  FF = FF + (2*MSb12*MSb22/MSb2diff**2*(x6 - x7) +
     &      MSb2sum/MSb2diff)*(dMSb12/MSb12 - dMSb22/MSb22)

	  dAb = Abr*dhbhb + 2*pMGl*c2

	  dPhi(se11) = s2tb*Abr*(Abr*ds2tb + s2tb*dAb)*FFo

	  dPhi(se12) = -s2tb*MUEr*(Abr*ds2tb +
     &      .5D0*s2tb*(dAb + Abr*dhbhb))*FFo

	  dPhi(se22) = s2tb*MUEr2*(ds2tb + s2tb*dhbhb)*FFo
	else
	  dPhi(se11) = 0
	  dPhi(se12) = 0
	  dPhi(se22) = 0
	endif

	k = C_F*Nc/(64*pi**3)*AlfasMT*hb2

	dPhi(se11) = k*(dPhi(se11) + .5D0*s2tb2*Abr*(Abr*FF + 2*pMGl*GG))
	dPhi(se12) = k*(dPhi(se12) - .5D0*s2tb2*MUEr*(Abr*FF + pMGl*GG))
	dPhi(se22) = k*(dPhi(se22) + .5D0*s2tb2*MUEr2*FF)

	if( debuglevel .lt. 6 ) return

	DHIGGS "asab: GG         =", GG		ENDL
	DHIGGS "asab: FF         =", FF		ENDL
	DHIGGS "asab:   se11     =", dPhi(se11)	ENDL
	DHIGGS "asab:   se12     =", dPhi(se12)	ENDL
	DHIGGS "asab:   se22     =", dPhi(se22)	ENDL
	end
