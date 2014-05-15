* the stop1-stop1 self-energy
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function Pi11t(p)
	implicit none
	RealType p

#include "TLps.h"

	Pi11t = 1/2.D0*(ht2*
     -      (-((-2 + s2tt2 + Nc*s2tt2)*myAA(MSt22,Q)) - 
     -        2*myAA(pMT2,Q) + 
     -        (2*pMT + Atr*s2tt)**2*myB0(p,MSt12,0D0,Q) + 
     -        s2tt2*((1 + Nc)*myAA(MSt12,Q) + 
     -           MUEr2*myB0(p,MSt12,MA02,Q)) + 
     -        (1 + c2tt2)*(Atr**2*myB0(p,MSt22,0D0,Q) + 
     -           MUEr2*myB0(p,MSt22,MA02,Q)) - 
     -        2*(pMT2 + MUEr2 - p)*myB0(p,pMT2,MUEr2,Q))) - 
     -   (ctt2*hb2 + ht2*stt2)*myAA(0D0,Q) + ctt2*hb2*myAA(MA02,Q) + 
     -   (ctt2*hb2*stb2 + ctb2*ht2*stt2)*myAA(MSb12,Q) + 
     -   (ctb2*ctt2*hb2 + ht2*stb2*stt2)*myAA(MSb22,Q) - 
     -   (ctt2*hb2 + ht2 + ht2*stt2)*myAA(MUEr2,Q) + 
     -   (-MUEr2 + p)*(ctt2*hb2 + ht2*stt2)*myB0(p,0D0,MUEr2,Q) + 
     -   (ctt2*hb2*MUEr2*stb2 + 
     -      ctb2*ht2*(ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2)) + 
     -      hb*ht*MUEr*s2tb*(ctt2*pMT + 1/2.D0*(Atr*s2tt)))*
     -    myB0(p,MSb12,0D0,Q) + 
     -   (Abr**2*ctt2*hb2*stb2 + Abr*hb2*pMT*s2tt*stb2 + 
     -      ctb2*ht2*MUEr2*stt2 + hb*ht*pMT*MUEr*s2tb*stt2 + 
     -      hb2*pMT2*stb2*stt2 + 1/2.D0*(Abr*hb*ht*MUEr*s2tb*s2tt))*
     -    myB0(p,MSb12,MA02,Q) + 
     -   (ctb2*ctt2*hb2*MUEr2 + 
     -      ht2*stb2*(ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2)) - 
     -      1/2.D0*(hb*ht*MUEr*s2tb*(2*ctt2*pMT + Atr*s2tt)))*
     -    myB0(p,MSb22,0D0,Q) + 
     -   (ht2*MUEr2*stb2*stt2 + 
     -      ctb2*hb2*(Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2) - 
     -      1/2.D0*(hb*ht*MUEr*s2tb*(Abr*s2tt + 2*pMT*stt2)))*
     -    myB0(p,MSb22,MA02,Q)
	end

