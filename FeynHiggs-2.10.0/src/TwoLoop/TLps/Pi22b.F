* the sbottom2-sbottom2 self-energy
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function Pi22b(p)
	implicit none
	RealType p

#include "TLps.h"

        Pi22b = 1/2.D0*(hb2*
     -      (-2*myAA(0D0,Q) + 2*(1 + ctb2)*myAA(MA02,Q) - 
     -        (-2 + s2tb2 + Nc*s2tb2)*myAA(MSb12,Q) + 
     -        2*(-MUEr2 + p)*myB0(p,0D0,MUEr2,Q) + 
     -        (1 + c2tb2)*(MUEr2*myB0(p,MSb12,0D0,Q) + 
     -           Abr**2*myB0(p,MSb12,MA02,Q)) + 
     -        s2tb2*((1 + Nc)*myAA(MSb22,Q) + 
     -           MUEr2*myB0(p,MSb22,0D0,Q) + 
     -           Abr**2*myB0(p,MSb22,MA02,Q)))) + 
     -   (ctb2*ctt2*hb2 + ht2*stb2*stt2)*myAA(MSt12,Q) + 
     -   (ctt2*ht2*stb2 + ctb2*hb2*stt2)*myAA(MSt22,Q) - 
     -   (ctb2*hb2 + ht2*stb2)*myAA(pMT2,Q) - 
     -   (hb2 + ctb2*hb2 + ht2*stb2)*myAA(MUEr2,Q) + 
     -   (ctb2*ctt2*hb2*MUEr2 + 
     -      ht2*stb2*(ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2)) - 
     -      1/2.D0*(hb*ht*MUEr*s2tb*(2*ctt2*pMT + Atr*s2tt)))*
     -    myB0(p,MSt12,0D0,Q) + 
     -   (ht2*MUEr2*stb2*stt2 + 
     -      ctb2*hb2*(Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2) - 
     -      1/2.D0*(hb*ht*MUEr*s2tb*(Abr*s2tt + 2*pMT*stt2)))*
     -    myB0(p,MSt12,MA02,Q) + 
     -   (ctb2*hb2*MUEr2*stt2 + 
     -      ht2*stb2*(Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2) + 
     -      1/2.D0*(hb*ht*MUEr*s2tb*(Atr*s2tt - 2*pMT*stt2)))*
     -    myB0(p,MSt22,0D0,Q) + 
     -   (ctt2*ht2*MUEr2*stb2 + 
     -      ctb2*hb2*(ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)) + 
     -      1/2.D0*(hb*ht*MUEr*s2tb*(-2*ctt2*pMT + Abr*s2tt)))*
     -    myB0(p,MSt22,MA02,Q) + 
     -   (2*hb*ht*pMT*MUEr*s2tb - 
     -      (pMT2 + MUEr2 - p)*(ctb2*hb2 + ht2*stb2))*
     -    myB0(p,pMT2,MUEr2,Q)
	end

