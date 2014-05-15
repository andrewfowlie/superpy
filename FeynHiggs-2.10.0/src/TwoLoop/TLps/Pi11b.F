* the sbottom1-sbottom1 self-energy
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function Pi11b(p)
	implicit none
	RealType p

#include "TLps.h"

        Pi11b = 1/2.D0*(hb2*
     -      (-2*myAA(0D0,Q) + 2*(1 + stb2)*myAA(MA02,Q) - 
     -        (-2 + s2tb2 + Nc*s2tb2)*myAA(MSb22,Q) + 
     -        2*(-MUEr2 + p)*myB0(p,0D0,MUEr2,Q) + 
     -        s2tb2*((1 + Nc)*myAA(MSb12,Q) + 
     -           MUEr2*myB0(p,MSb12,0D0,Q) + 
     -           Abr**2*myB0(p,MSb12,MA02,Q)) + 
     -        (1 + c2tb2)*(MUEr2*myB0(p,MSb22,0D0,Q) + 
     -           Abr**2*myB0(p,MSb22,MA02,Q)))) + 
     -   (ctt2*hb2*stb2 + ctb2*ht2*stt2)*myAA(MSt12,Q) + 
     -   (ctb2*ctt2*ht2 + hb2*stb2*stt2)*myAA(MSt22,Q) - 
     -   (ctb2*ht2 + hb2*stb2)*myAA(pMT2,Q) - 
     -   (hb2 + ctb2*ht2 + hb2*stb2)*myAA(MUEr2,Q) + 
     -   (ctt2*hb2*MUEr2*stb2 + 
     -      ctb2*ht2*(ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2)) + 
     -      hb*ht*MUEr*s2tb*(ctt2*pMT + 1/2.D0*(Atr*s2tt)))*
     -    myB0(p,MSt12,0D0,Q) + 
     -   (Abr**2*ctt2*hb2*stb2 + Abr*hb2*pMT*s2tt*stb2 + 
     -      ctb2*ht2*MUEr2*stt2 + hb*ht*pMT*MUEr*s2tb*stt2 + 
     -      hb2*pMT2*stb2*stt2 + 1/2.D0*(Abr*hb*ht*MUEr*s2tb*s2tt))*
     -    myB0(p,MSt12,MA02,Q) + 
     -   (hb2*MUEr2*stb2*stt2 + 
     -      ctb2*ht2*(Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2) + 
     -      hb*ht*MUEr*s2tb*(pMT*stt2 - 1/2.D0*(Atr*s2tt)))*
     -    myB0(p,MSt22,0D0,Q) + 
     -   (ctt2*(ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb + 
     -         hb2*pMT2*stb2) + Abr**2*hb2*stb2*stt2 - 
     -      1/2.D0*(Abr*s2tt*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2)))*
     -    myB0(p,MSt22,MA02,Q) + 
     -   (-2*hb*ht*pMT*MUEr*s2tb - 
     -      (pMT2 + MUEr2 - p)*(ctb2*ht2 + hb2*stb2))*
     -    myB0(p,pMT2,MUEr2,Q)
	end

