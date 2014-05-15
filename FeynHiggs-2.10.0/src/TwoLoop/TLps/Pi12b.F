* the sbottom1-sbottom2 self-energy
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function Pi12b(p)
	implicit none
	RealType p

#include "TLps.h"

	Pi12b = 1/2.D0*(s2tb*
     -       ((ctt2*hb2 - ht2*stt2)*myAA(MSt12,Q) + 
     -         (-(ctt2*ht2) + hb2*stt2)*myAA(MSt22,Q) + 
     -         (-hb2 + ht2)*(myAA(pMT2,Q) + myAA(MUEr2,Q)) + 
     -         hb2*(myAA(MA02,Q) + 
     -            c2tb*((1 + Nc)*(myAA(MSb12,Q) - myAA(MSb22,Q)) + 
     -               MUEr2*
     -                (myB0(p,MSb12,0D0,Q) - myB0(p,MSb22,0D0,Q)) + 
     -               Abr**2*
     -                (myB0(p,MSb12,MA02,Q) - myB0(p,MSb22,MA02,Q))
     -               ))) + 
     -      (c2tb*hb*ht*MUEr*(2*ctt2*pMT + Atr*s2tt) + 
     -         s2tb*(ctt2*hb2*MUEr2 - 
     -            ht2*(ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2))))*
     -       myB0(p,MSt12,0D0,Q) + 
     -      (Abr**2*ctt2*hb2*s2tb + 
     -         Abr*(c2tb*hb*ht*MUEr + hb2*pMT*s2tb)*s2tt + 
     -         (2*c2tb*hb*ht*pMT*MUEr + hb2*pMT2*s2tb - 
     -            ht2*MUEr2*s2tb)*stt2)*myB0(p,MSt12,MA02,Q) + 
     -      (c2tb*hb*ht*MUEr*(-(Atr*s2tt) + 2*pMT*stt2) + 
     -         s2tb*(-(Atr**2*ctt2*ht2) + Atr*ht2*pMT*s2tt - 
     -            ht2*pMT2*stt2 + hb2*MUEr2*stt2))*myB0(p,MSt22,0D0,Q)+
     -        (c2tb*hb*ht*MUEr*(2*ctt2*pMT - Abr*s2tt) + 
     -         s2tb*(ctt2*(hb2*pMT2 - ht2*MUEr2) + 
     -            Abr*hb2*(-(pMT*s2tt) + Abr*stt2)))*
     -       myB0(p,MSt22,MA02,Q)) + 
     -   (-2*c2tb*hb*ht*pMT*MUEr + 
     -      1/2.D0*((-hb2 + ht2)*(pMT2 + MUEr2 - p)*s2tb))*
     -    myB0(p,pMT2,MUEr2,Q)
	end
