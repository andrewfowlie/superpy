* DRbartoOS.F
* this file is part of FeynHiggs
* last modified 30 Jul 13 th

* Code by Pietro Slavich
* 10/04/2013: different Ab conversion according to hep-ph/0206101
* 28/05/2004: error in charg12 corrected
* 18/06/2004: rewritten by TH for FeynHiggs

#include "externals.h"
#include "types.h"
#include "debug.h"

#define __SUBROUTINE__ "DRbartoOS"


	subroutine DRbartoOS(t, i, ospara)
	implicit none
	integer t, i
	RealType ospara(3)

#include "FH.h"
#include "looptools.h"

	integer try, b, ty, by
	RealType Q, dtt
	RealType stt, stt2, ctt, ctt2, s2tt, c2tt
	RealType stb, stb2, ctb, ctb2, s2tb, c2tb
	RealType MT_, MT2_, MB_, MB2_
	RealType TB_, CB_, SB_, CB2_, SB2_
	RealType ht, hb, Xt, Xb, Yt, Yb, MUE2
	RealType MSt12, MSt22, MSb12, MSb22
	RealType dMSt12, dMSt12_s, dMSt12_y
	RealType dMSt22, dMSt22_s, dMSt22_y
	RealType astt, dtt_s, dtt_y

	integer Nc
	parameter (Nc = 3)

	RealType Pitt11_s, Pitt11_y, Pitt12_s, Pitt12_y
	RealType Charge11, Charge12, p

* Eq. (B3) of DSZ
	Pitt11_s(p, MSt12, MSt22, s2tt) = AlfasMT/(3*pi)*(
     &    MSt12*((3 - c2tt**2)*(log(MSt12/Q) - 1) - 4) -
     &    MSt22*s2tt**2*(log(MSt22/Q) - 1) +
     &    2*(MGl2*(log(MGl2/Q) - 1) +
     &       MT2_*(log(MT2_/Q) - 1) +
     &       (p - MGl2 - MT2_ + s2tt*MT_*Re(M_3 + M_3C))*
     &         myB0(p, MT2_, MGl2, Q)) )

* Eq. (B7) of DSZ
	Pitt12_s(p) = AlfasMT/(3*pi)*(
     &    4*MT_*Re(M_3*ctt2 - M_3C*stt2)*myB0(p, MT2_, MGl2, Q) +
     &    c2tt*s2tt*(MSt22*(1 - log(MSt22/Q)) -
     &               MSt12*(1 - log(MSt12/Q))) )

	Charge11(s2tt, stt2, ctt2, s2tb, c2tb, stb2, ctb2,
     &           SB_, CB_, Xt, Xb) =
     &    (ht*SB_)**2*(ctt2*ctb2*MT2_ + stt2*stb2*MB2_ +
     &      .5D0*s2tt*s2tb*MB_*MT_ +
     &      Xt*(ctb2*s2tt*MT_ + stt2*s2tb*MB_ + stt2*ctb2*Xt)) +
     &    (hb*CB_)**2*(ctb2*ctt2*MB2_ + stb2*stt2*MT2_ +
     &      .5D0*s2tb*s2tt*MB_*MT_ +
     &      Xb*(ctt2*s2tb*MB_ + stb2*s2tt*MT_ + stb2*ctt2*Xb)) -
     &    ht*hb*SB_*CB_*(2*(ctt2*ctb2 + stt2*stb2)*MT_*MB_ +
     &      .5D0*s2tt*s2tb*(MT2_ + MB2_ + Xt*Xb) +
     &      s2tb*MT_*(ctt2*Xb + stt2*Xt) +
     &      s2tt*MB_*(ctb2*Xt + stb2*Xb))

	Pitt11_y(p, MSt12, MSt22, s2tt, c2tt, stt2, ctt2) =
     &  1/(16*pi**2)*(
     &    ht**2*(myG(p, MT2_, MUE2, Q) +
     &      stt2*(ctb2*myAA(MSb12, Q) + stb2*myAA(MSb22, Q)) +
     &      (c2tt**2 - .5D0*(Nc - 1)*s2tt**2)*myAA(MSt22, Q) +
     &      .5D0*(Nc + 1)*s2tt**2*myAA(MSt12, Q) +
     &      CB2_*(1 + stt2)*myAA(MA02, q) +
     &      .5D0*(SB2_*(2*MT_ + s2tt*Xt)**2*myB0(p, MSt12, 0D0, Q) +
     &         CB2_*(2*MT_ + s2tt*Yt)**2*myB0(p, MSt12, MA02, Q) +
     &         (1 + c2tt**2)*(SB2_*Xt**2*myB0(p, MSt22, 0D0, Q) +
     &                        CB2_*Yt**2*myB0(p, MSt22, MA02, Q)))) +
     &    (ht**2*stt2 + hb**2*ctt2)*myG(p, MB2_, MUE2, Q) -
     &    2*ht*hb*MB_*reMUE*s2tt*myB0(p, MB2_, MUE2, Q) +
     &    hb**2*ctt2*(SB2_*myAA(MA02, Q) +
     &      stb2*myAA(MSb12, Q) + ctb2*myAA(MSb22, Q)) +
     &    Charge11(s2tt, stt2, ctt2, s2tb, c2tb, stb2, ctb2,
     &      SB_, CB_, Xt, Xb)*myB0(p, MSb12, 0D0, Q) +
     &    Charge11(s2tt, stt2, ctt2, -s2tb, -c2tb, ctb2, stb2,
     &      SB_, CB_, Xt, Xb)*myB0(p, MSb22, 0D0, Q) +
     &    Charge11(s2tt, stt2, ctt2, s2tb, c2tb, stb2, ctb2,
     &      CB_, -SB_, Yt, Yb)*myB0(p, MSb12, MA02, Q) +
     &    Charge11(s2tt, stt2, ctt2, -s2tb, -c2tb, ctb2, stb2,
     &      CB_, -SB_, Yt, Yb)*myB0(p, MSb22, MA02, Q) )

	Charge12(s2tb, c2tb, stb2, ctb2, SB_, CB_, Xt, Xb) =
     &    .5D0*(ht*SB_)**2*(s2tt*(stb2*MB2_ - ctb2*MT2_) + 
     &      Xt*(2*c2tt*ctb2*MT_ + s2tt*(s2tb*MB_ + ctb2*Xt)) +
     &      c2tt*s2tb*MT_*MB_) +
     &    .5D0*(hb*CB_)**2*(s2tt*(stb2*MT2_ - ctb2*MB2_) + 
     &      Xb*(2*c2tt*ctb2*MT_ - s2tt*(s2tb*MB_ + stb2*Xb)) +
     &      c2tt*s2tb*MT_*MB_) -
     &    .5D0*ht*hb*SB_*CB_*(c2tt*s2tb*(MT2_ + MB2_ + Xt*Xb) +
     &      s2tt*s2tb*MT_*(Xt - Xb) +
     &      2*MB_*(c2tt*(ctb2*Xt + stb2*Xb) - s2tt*c2tb*MT_))

	Pitt12_y(p) = 1/(16*pi**2)*(
     &    .5D0*ht**2*(
     &      s2tt*(CB2_*myAA(MA02, Q) +
     &        (Nc + 1)*c2tt*(myAA(MSt12, Q) - myAA(MSt22, Q)) +
     &        ctb2*myAA(MSb12, Q) + stb2*myAA(MSb22, Q)) +
     &      c2tt*(SB2_*Xt*(2*MT_ + s2tt*Xt)*myB0(p, MSt12, 0D0, Q) +
     &        CB2_*Yt*(2*MT_ + s2tt*Yt)*myB0(p, MSt12, MA02, Q) +
     &        SB2_*Xt*(2*MT_ - s2tt*Xt)*myB0(p, MSt22, 0D0, Q) +
     &        CB2_*Yt*(2*MT_ - s2tt*Yt)*myB0(p, MSt22, MA02, Q))) +
     &    .5D0*(ht**2 - hb**2)*s2tt*myG(p, MB2_, MUE2, Q) -
     &    2*ht*hb*MB_*reMUE*c2tt*myB0(p, MB2_, MUE2, Q) -
     &    .5D0*hb**2*s2tt*(SB2_*myAA(MA02, Q) +
     &      stb2*myAA(MSb12, Q) + ctb2*myAA(MSb22, Q)) +
     &    Charge12(s2tb, c2tb, stb2, ctb2, SB_, CB_, Xt, Xb)*
     &      myB0(p, MSb12, 0D0, Q) +
     &    Charge12(-s2tb, -c2tb, ctb2, stb2, SB_, CB_, Xt, Xb)*
     &      myB0(p, MSb22, 0D0, Q) +
     &    Charge12(s2tb, c2tb, stb2, ctb2, CB_, -SB_, Yt, Yb)*
     &      myB0(p, MSb12, MA02, Q) +
     &    Charge12(-s2tb, -c2tb, ctb2, stb2, CB_, -SB_, Yt, Yb)*
     &      myB0(p, MSb22, MA02, Q) )

	b = 7 - t
* the caller is assumed to have set the DRbar scale with SetDecayScale:
	ty = t + (tH - 3)*isQ(t)
	by = b + (bH - 4)*isQ(t)

	MSt12 = MSf2(1,ty,3)
	MSt22 = MSf2(2,ty,3)
	stt = Re(USf(1,2,ty,3))
	astt = asin(stt)
	Q = QSf(t)**2

	if( debuglevel .ge. 2 ) then
	  DPARA "MSf("//Digit(t)//")_DRbar =", sqrt(MSt12), sqrt(MSt22) ENDL
	  DPARA "ssf("//Digit(t)//")_DRbar =", stt ENDL
	endif

	MSb12 = MSf2(1,by,3)
	MSb22 = MSf2(2,by,3)
	stb = Re(USf(1,2,by,3))

	stb2 = stb**2
	ctb2 = (1 - stb)*(1 + stb)
	ctb = sqrt(ctb2)
	s2tb = 2*ctb*stb
	c2tb = (ctb - stb)*(ctb + stb)

	MT_ = Mf(t,3)
	MB_ = Mf(b,3)
c future extension (hr+kw)
c	MT_ = Mf(ty,3)
c	MB_ = Mf(by,3)
	MT2_ = MT_**2
	MB2_ = MB_**2

	CB_ = SCB(b)
	CB2_ = CB_**2
	SB_ = SCB(t)
	SB2_ = SB_**2
	TB_ = SB_/CB_

	ht = sqrt2/vev*MT_/SB_
	hb = sqrt2/vev*MB_/CB_

	Xt = Re(Af0(t,3) + MUETB(t))
	Yt = Xt + reMUE/(SB_*CB_)
	Xb = Re(Af0(b,3) + MUETB(b))
	Yb = Xb + reMUE/(SB_*CB_)

	MUE2 = reMUE**2

	dMSt12 = 0
	dMSt22 = 0
	dtt = 0

c	do try = 1, 20
	do try = 1, 1
	  stt2 = stt**2
	  ctt2 = (1 - stt)*(1 + stt)
	  ctt = sqrt(ctt2)
	  s2tt = 2*ctt*stt
	  c2tt = (ctt - stt)*(ctt + stt)

* strong shifts

	  dMSt12_s = Pitt11_s(MSt12, MSt12, MSt22, s2tt)

	  dMSt22_s = Pitt11_s(MSt22, MSt22, MSt12, -s2tt)

	  dtt_s = .5D0/(MSt12 - MSt22)*
     &      (Pitt12_s(MSt12) + Pitt12_s(MSt22))

* Yukawa shifts

	  dMSt12_y = Pitt11_y(MSt12,
     &      MSt12, MSt22, s2tt, c2tt, stt2, ctt2)

	  dMSt22_y = Pitt11_y(MSt22,
     &      MSt22, MSt12, -s2tt, -c2tt, ctt2, stt2)

	  dtt_y = .5D0/(MSt12 - MSt22)*
     &      (Pitt12_y(MSt12) + Pitt12_y(MSt22))

	  if( debuglevel .ge. 3 ) then
	    DPARA "MSf21    =", MSt12		ENDL
	    DPARA "dMSf12_s =", dMSt12_s	ENDL
	    DPARA "dMSf12_y =", dMSt12_y	ENDL
	    DPARA "MSf22    =", MSt22		ENDL
	    DPARA "dMSf22_s =", dMSt22_s	ENDL
	    DPARA "dMSf22_y =", dMSt22_y	ENDL
	    DPARA "tf       =", astt		ENDL
	    DPARA "dtf_s    =", dtt_s		ENDL
	    DPARA "dtf_y    =", dtt_y		ENDL
	  endif

	  dMSt12_y = dMSt12_y + dMSt12_s
	  dMSt22_y = dMSt22_y + dMSt22_s
	  dtt_y = dtt_y + dtt_s
	  if( MSt12 .lt. dMSt12_y .or.
     &        MSt22 .lt. dMSt22_y ) then
	    if( MSt12 .lt. dMSt12_s .or.
     &          MSt22 .lt. dMSt22_s ) return
	    Warning("Using DRbartoOS without Yukawa corrections")
	    dMSt12_y = dMSt12_s
	    dMSt22_y = dMSt22_s
	    dtt_y = dtt_s
	  endif

	  dMSt12 = dMSt12 + dMSt12_y
	  dMSt22 = dMSt22 + dMSt22_y
	  dtt = dtt + dtt_y

	  MSt12 = MSf2(1,ty,3) - dMSt12
	  MSt22 = MSf2(2,ty,3) - dMSt22
	  stt = sin(astt - dtt)

	  if( abs(dMSt12_y) + abs(dMSt22_y) + abs(dtt_y) .lt. 1D-8 )
     &      goto 1
	enddo

1	continue

* From X = USf^+ diag(MSf) USf we get
*   DSS2 = X11 - X22
*   DSS2 = X22 - Mf^2
*   Xf = Conjugate(X12)/Mf

	ctt = sqrt((1 - stt)*(1 + stt))
	ospara(1) = (MSt12 - MSt22)*(ctt - stt)*(ctt + stt)
	ospara(2) = MSt12*stt**2 + MSt22*ctt**2 - Mf2(t,3)
	if( t .eq. 4 ) then
* hep-ph/0206101 Eq. (19)
	  ospara(3) = Xt + C_F/(2*pi)*AlfasMT*MGl*(4 - 2*log(MGl2/Q) -
     &      (1 - MGl2/MSt12)*log(abs(1 - MSt12/MGl2)) -
     &      (1 - MGl2/MSt22)*log(abs(1 - MSt22/MGl2)))
	else
	  ospara(3) = (MSt12 - MSt22)*ctt*stt/Mf(t,3)
	endif
	i = ibset(i, t)

	if( debuglevel .ge. 2 ) then
	  DPARA "MSf("//Digit(t)//")_OS    =", sqrt(MSt12), sqrt(MSt22) ENDL
	  DPARA "ssf("//Digit(t)//")_OS    =", stt ENDL
	endif
	end

