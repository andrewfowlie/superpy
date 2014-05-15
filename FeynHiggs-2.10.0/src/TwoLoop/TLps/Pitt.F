* the top self-energy
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function Pitt(p)
	implicit none
	RealType p

#include "TLps.h"

        Pitt = 1/2.D0*(hb*ht*MUEr*s2tb*
     -      (myB0(p,MUEr2,MSb12,Q) - myB0(p,MUEr2,MSb22,Q)) +
     -     pMT*(hb2*myB1(p,0.D0,MA02,Q) +
     -        (ctb2*ht2 + hb2*stb2)*myB1(p,MUEr2,MSb12,Q) +
     -        (ctb2*hb2 + ht2*stb2)*myB1(p,MUEr2,MSb22,Q) +
     -        ht2*(myB1(p,0.D0,0.D0,Q) + 2*myB1(p,pMT2,0.D0,Q) +
     -           myB1(p,MUEr2,MSt12,Q) + myB1(p,MUEr2,MSt22,Q))))
	end
