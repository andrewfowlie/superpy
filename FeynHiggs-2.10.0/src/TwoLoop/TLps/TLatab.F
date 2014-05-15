* TLatab.F
* Two-loop O(at ab + ab^2) corrections to the Higgs masses
* in the on-shell renormalization scheme
* Written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
* Based on A. Dedes, G. Degrassi, and P. Slavich, hep-ph/0305127.
* this file is part of FeynHiggs
* last modified 10 Dec 13 th

* Revision history:
* 14/07/2004: rewritten completely by TH for FeynHiggs
* 28/05/2004: error in fPi12t corrected
* 25/09/2003: error in F1b corrected
*             numerical value of Pi corrected (10th digit)
* 20/07/2003: routines optimized, bug for mu<0 corrected
*
* Notice: we assume that the 1-loop part is computed in terms of
*         OS parameters, defined as described in the paper.  In
*         particular, the parameters in the bottom/sbottom sector
*         should be computed in terms of the renormalized coupling hb.
*         Also, mL_sbot differs from mL_stop at the one-loop level.

#include "externals.h"
#include "types.h"
#include "debug.h"


* COMPLETE_REN adds also the renormalization for MUEr and the
* wave function ren. which are usually zero in DRbar
c#define COMPLETE_REN


	subroutine TLatab(dPhi)
	implicit none
	RealType dPhi(3)

#include "TLps.h"

	RealType F1t, s2ttF2t, s2tt2F3t, s2ttF4t, F1b, s2tbF2b
	RealType s2tb2F3b, s2tbF4b, s2tts2tbF5, F6
	RealType k, F1tol, F2tol, F2bol, F3tol, F3bol
	RealType dMSt12, dMSt22, dMSb12, dMSb22, pi1, pi2
	RealType dpMTpMT, dvev2vev2
	RealType ds2tt, ds2tb, dhtht, s2tbdhbhb
	RealType s2ttdAtrAtr, s2ttdAtrF2t, s2ttdAtrF3t
	RealType dAbthb, s2tbdAb, s2tbAbdAbF3b

#ifdef COMPLETE_REN
	RealType dMUErMUEr, dZ1, dZ2, M11ol, M12ol, M22ol
#endif
	RealType dMUErF2t, dMUErF3t, dMUErF3b
	RealType WF11, WF12, WF22

	RealType x44, x48, x11, x12, x46, x13, x14, x16, x15
	RealType x17, x23, x25, x26, x27, x50, x33, x34, x35
	RealType x36, x43, x8, x51, x52, x6, x7, x9, x21, x22
	RealType x18, x19, x20, x24, x10, x3, x4, x5, x1, x45
	RealType x49, x28, x29, x47, x30, x31, x38, x39, x40
	RealType x32, x37, x41, x42

	RealType tmp1, tmp12, tmp13, tmp14
	RealType tmp15, tmp16, tmp17, tmp18, tmp19, tmp2
	RealType tmp20, tmp21, tmp22, tmp23, tmp24, tmp25
	RealType tmp26, tmp27, tmp28, tmp29, tmp3, tmp30
	RealType tmp31, tmp32, tmp4, tmp5, tmp6, tmp7, tmp8
	RealType tmp9
	RealType lc11, lc12, lc13, lc14, lc19, lc20, lc21
	RealType lc22, lc23, lc24, lc25, lc3, lc4
	RealType lc5, lc6, lc7, lc8, lc9

	RealType Pi11t, Pi22t, Pi12t
	RealType Pi11b, Pi22b, Pi12b, Pitt
	external Pi11t, Pi22t, Pi12t, Pi11b, Pi22b, Pi12b, Pitt

	RealType p, Lambhb, PiAAb

	RealType x, y, fW
	fW(x, y) = .25D0*(x + y) - .5D0*x*y/(x - y)*log(x/y)

        Lambhb(p) = -(1/2.D0*
     &     (2*Abr*hb2*MSt2diff*Nc*x1 + 
     &        2*Abr*hb2*MSt2diff*myB0(p,MSb12,MA02,Q) + 
     &        2*Abr*hb2*MSt2diff*myB0(p,MSb22,MA02,Q) - 
     &        2*Atr*ht2*MUEr2*myB0(p,MSt12,0.D0,Q) + 
     &        ht2*MSt2diff*pMT*s2tt*myB0(p,MSt12,0.D0,Q) + 
     &        2*Atr*ht2*MSt2diff*stt2*myB0(p,MSt12,0.D0,Q) + 
     &        2*Abr*ctt2*hb2*MSt2diff*myB0(p,MSt12,MA02,Q) - 
     &        2*Abr*ht2*MUEr2*myB0(p,MSt12,MA02,Q) + 
     &        hb2*MSt2diff*pMT*s2tt*myB0(p,MSt12,MA02,Q) + 
     &        2*Atr*ctt2*ht2*MSt2diff*myB0(p,MSt22,0.D0,Q) + 
     &        2*Atr*ht2*MUEr2*myB0(p,MSt22,0.D0,Q) - 
     &        ht2*MSt2diff*pMT*s2tt*myB0(p,MSt22,0.D0,Q) + 
     &        2*Abr*ht2*MUEr2*myB0(p,MSt22,MA02,Q) - 
     &        hb2*MSt2diff*pMT*s2tt*myB0(p,MSt22,MA02,Q) + 
     &        2*Abr*hb2*MSt2diff*stt2*myB0(p,MSt22,MA02,Q))/
     &      (MSt2diff*Sqrt(2.D0)))

        PiAAb(p) = 1/2.D0*(Nc*
     &     (hb2*(-4*x3 + 2*(x4 + x5 + p*myB0(p,0.D0,0.D0,Q)) + 
     &          Abr**2*myB0(p,MSb12,MSb22,Q)) + 
     &       ht2*MUEr2*myB0(p,MSt12,MSt22,Q)))

	x44 = InvTri(MA02,MSb12,MSb12)

	x48 = InvTri(MA02,MSb12,MSb22)

	x11 = InvTri(MA02,MSb12,MSt12)

	x12 = InvTri(MA02,MSb12,MSt22)

	x46 = InvTri(MA02,MSb22,MSb22)

	x13 = InvTri(MA02,MSb22,MSt12)

	x14 = InvTri(MA02,MSb22,MSt22)

	x16 = InvTri(MA02,MSt12,MSt12)

	x15 = InvTri(MA02,MSt12,MSt22)

	x17 = InvTri(MA02,MSt22,MSt22)

	x23 = InvTri(MSb12,pMT2,MUEr2)

	x25 = InvTri(MSb22,pMT2,MUEr2)

	x26 = InvTri(MSt12,pMT2,MUEr2)

	x27 = InvTri(MSt22,pMT2,MUEr2)

	x50 = Li2omx(MSb12/MSb22)

	x33 = Li2omx(MSt12/MSb12)

	x34 = Li2omx(MSt12/MSb22)

	x35 = Li2omx(MSt22/MSb12)

	x36 = Li2omx(MSt22/MSb22)

	x43 = Li2omx(MSt12/MSt22)

	x8 = Li2omx(pMT2/MA02)

	x51 = Li2omx(MSb12/MUEr2)

	x52 = Li2omx(MSb22/MUEr2)

	x6 = Li2omx(MSt12/MUEr2)

	x7 = Li2omx(MSt22/MUEr2)

	x9 = Log(MA02)

	x21 = Log(MSb12)

	x22 = Log(MSb22)

	x18 = Log(MSt12)

	x19 = Log(MSt22)

	x20 = Log(pMT2)

	x24 = Log(MUEr2)

	x10 = Log(Q)

	x3 = myAA(0.D0,Q)

	x4 = myAA(MSb12,Q)

	x5 = myAA(MSb22,Q)

	x1 = myB0(0.D0,MSb12,MSb22,Q)

	x45 = Phi(MA02,MSb12,MSb12)

	x49 = Phi(MA02,MSb12,MSb22)

	x28 = Phi(MA02,MSb12,MSt12)

	x29 = Phi(MA02,MSb12,MSt22)

	x47 = Phi(MA02,MSb22,MSb22)

	x30 = Phi(MA02,MSb22,MSt12)

	x31 = Phi(MA02,MSb22,MSt22)

	x38 = Phi(MA02,MSt12,MSt12)

	x39 = Phi(MA02,MSt12,MSt22)

	x40 = Phi(MA02,MSt22,MSt22)

	x32 = Phi(MSb12,pMT2,MUEr2)

	x37 = Phi(MSb22,pMT2,MUEr2)

	x41 = Phi(MUEr2,pMT2,MSt12)

	x42 = Phi(MUEr2,pMT2,MSt22)

        tmp1 = 1/4.D0*((-2*ctb2*ht2*(2*MSb12 + MUEr2)*stt2)/
     &      (MSb12 - MSt12) + 
     &     ctb2*ht2*((2*(ctt2*(2*MSb12 + MUEr2) + 
     &             (MSb12*(c2tt*pMT2 + 2*Atr*pMT*s2tt))/MSt22))/
     &         (-MSb12 + MSt22) + 
     &        (2*MSb12*(2*ctt2*pMT2 + 2*Atr*pMT*s2tt + 2*Atr**2*stt2 + 
     &             (2*MSb12 + MUEr2)*stt2))/((MSb12 - MSt12)*MSt12))
     &       + ((hb*ht*MUEr*s2tb*
     &           (MSb12**2*
     &              (-4*MSt12*MSt22 + c2tt*MSt2diff*pMT2 + 
     &                3*MSt2sum*pMT2 + Abr*MSt2diff*pMT*s2tt + 
     &                2*Atr*MSt2diff*pMT*s2tt) + 
     &             MSb12*(-((3 + c2tt)*MSt22**2*pMT2) + 
     &                MSt12**2*(4*MSt22 + (-3 + c2tt)*pMT2) - 
     &                Abr*MSt22*MSt2diff*pMT*s2tt - 
     &                2*Atr*MSt2diff*MSt2sum*pMT*s2tt + 
     &                MSt12*
     &                 (4*MSt22**2 - 2*MSt22*pMT2 - 
     &                   Abr*MSt2diff*pMT*s2tt)) + 
     &             MSt12*MSt22*
     &              (-4*MSt12*MSt22 + 2*ctt2*MSt12*pMT2 + 
     &                Abr*MSt2diff*pMT*s2tt + 2*MSt22*pMT2*stt2)))/pMT+
     &          (hb2*stb2*
     &           (MSt12*MSt22*
     &              ((-(c2tt*MSt2diff) + MSt2sum)*pMT2**2 - 
     &                2*MSt12*MSt22*(pMT2 + MUEr2)) + 
     &             2*(MSb12**3*
     &                 (-(MSt12*MSt22) + 
     &                   (c2tt*MSt2diff + MSt2sum)*pMT2) + 
     &                Abr*(MSb12 - MSt12)*(MSb12 - MSt22)*pMT2*
     &                 (Abr*ctt2*MSt22 + MSt2diff*pMT*s2tt + 
     &                   Abr*MSt12*stt2)) + 
     &             MSb12*(pMT2*
     &                 ((c2tt*MSt2diff - MSt2sum)*MSt2sum*pMT2 - 
     &                   4*ctt2*MSt22**2*MUEr2) + 
     &                2*MSt12*MSt22*
     &                 (2*ctt2*MSt22*pMT2 + MSt2sum*(pMT2 + MUEr2)) - 
     &                2*MSt12**2*
     &                 (MSt22**2 - 2*MSt22*pMT2*stt2 + 
     &                   2*pMT2*MUEr2*stt2)) + 
     &             MSb12**2*
     &              (2*MSt12*MSt22*(MSt2sum - 3*pMT2 - MUEr2) + 
     &                pMT2*((-(c2tt*MSt2diff) + MSt2sum)*pMT2 + 
     &                   2*(c2tt*MSt2diff + MSt2sum)*MUEr2 - 
     &                   4*(ctt2*MSt22**2 + MSt12**2*stt2)))))/pMT2)
     &       /((MSb12 - MSt12)*MSt12*(MSb12 - MSt22)*MSt22) + 
     &     2*(-(ctb2*ht2) - (ctb2*ht2*(MSb12 + MUEr2))/pMT2 + 
     &        (ctb2*ht2*MSb12*
     &           (pMT2 + ctt2*(2*(Atr**2 + MSb12) + MUEr2)))/
     &         ((MSb12 - MSt22)*MSt22) + 
     &        (((1 + ctt2)*hb*ht*MUEr*s2tb + 
     &              2*ctb2*ht2*(2*ctt2*pMT + Atr*s2tt) + 
     &              hb2*stb2*(Abr*s2tt + 2*pMT*stt2))*(-x10 + x18) + 
     &           (-2*Atr*ctb2*ht2*s2tt + 2*ctt2*hb2*pMT*stb2 - 
     &              Abr*hb2*s2tt*stb2 + 4*ctb2*ht2*pMT*stt2 + 
     &              hb*ht*MUEr*s2tb*(1 + stt2))*(-x10 + x19))/pMT) + 
     &     (4*(ctb2*ht2*pMT + hb*ht*MUEr*s2tb + hb2*pMT*stb2)*
     &        (x10 - x20))/pMT)

        tmp1 = tmp1 + (ctb2*ht2 + (hb*ht*MUEr*s2tb)/pMT + 
     &      hb2*stb2)*(-x10 + x24) + 
     &   1/2.D0*((-MSb12**2 + (pMT2 - MUEr2)**2)*
     &       (ctb2*ht2*pMT*(-MSb12 + pMT2 + MUEr2) + 
     &         2*hb*ht*pMT2*MUEr*s2tb + 
     &         hb2*pMT*(-MSb12 + pMT2 + MUEr2)*stb2)*x23)/(pMT*pMT2)

        tmp2 = 1/8.D0*(ht2*
     &      ((MSb22*(2*Atr**2*(1 + c2tb*c2tt) + 
     &             2*(1 + c2tb*c2tt)*MSb22 + 4*ctt2*pMT2 + MUEr2 + 
     &             c2tb*c2tt*MUEr2 + 4*Atr*pMT*s2tt))/
     &         ((MSb22 - MSt12)*MSt12) - 
     &        (2*(2*Atr**2*c2tb*ctt2*MSb22 + 2*ctt2*MSb22*MSt22 + 
     &             c2tt*MSb22*pMT2 + ctt2*MSt22*MUEr2 + 
     &             c2tb*(MSb22*pMT2 + 
     &                ctt2*(MSb22 - MSt22)*(2*MSb22 + MUEr2)) + 
     &             2*Atr*MSb22*pMT*s2tt))/((MSb22 - MSt22)*MSt22) + 
     &        MSb22*((2*(pMT2 + ctt2*(2*(Atr**2 + MSb22) + MUEr2) + 
     &                c2tb*(c2tt*pMT2 + 2*Atr*pMT*s2tt)))/
     &            ((MSb22 - MSt22)*MSt22) + 
     &           (2*(c2tb + c2tt)*(Atr**2 + MSb22) + c2tt*MUEr2 + 
     &              c2tb*(4*ctt2*pMT2 + MUEr2 + 4*Atr*pMT*s2tt))/
     &            (MSt12*(-MSb22 + MSt12))) - 
     &        (4*(MSb22 + pMT2 + MUEr2)*stb2)/pMT2 - 
     &        (4*(2*MSb22 + MUEr2)*stb2*stt2)/(MSb22 - MSt12)) + 
     &     (2*(-((hb*ht*MUEr*s2tb*
     &               (MSb22**2*
     &                  (-4*MSt12*MSt22 + c2tt*MSt2diff*pMT2 + 
     &                    3*MSt2sum*pMT2 + Abr*MSt2diff*pMT*s2tt + 
     &                    2*Atr*MSt2diff*pMT*s2tt) + 
     &                 MSb22*
     &                  (-((3 + c2tt)*MSt22**2*pMT2) + 
     &                    MSt12**2*(4*MSt22 + (-3 + c2tt)*pMT2) - 
     &                    Abr*MSt22*MSt2diff*pMT*s2tt - 
     &                    2*Atr*MSt2diff*MSt2sum*pMT*s2tt + 
     &                    MSt12*
     &                     (4*MSt22**2 - 2*MSt22*pMT2 - 
     &                       Abr*MSt2diff*pMT*s2tt)) + 
     &                 MSt12*MSt22*
     &                  (-4*MSt12*MSt22 + 2*ctt2*MSt12*pMT2 + 
     &                    Abr*MSt2diff*pMT*s2tt + 2*MSt22*pMT2*stt2)))
     &              /pMT) + 
     &          (ctb2*hb2*(MSt12*MSt22*
     &                ((-(c2tt*MSt2diff) + MSt2sum)*pMT2**2 - 
     &                  2*MSt12*MSt22*(pMT2 + MUEr2)) + 
     &               2*(MSb22**3*
     &                   (-(MSt12*MSt22) + 
     &                     (c2tt*MSt2diff + MSt2sum)*pMT2) + 
     &                  Abr*(MSb22 - MSt12)*(MSb22 - MSt22)*pMT2*
     &                   (Abr*ctt2*MSt22 + MSt2diff*pMT*s2tt + 
     &                     Abr*MSt12*stt2)) + 
     &               MSb22*
     &                (pMT2*
     &                   ((c2tt*MSt2diff - MSt2sum)*MSt2sum*pMT2 - 
     &                     4*ctt2*MSt22**2*MUEr2) + 
     &                  2*MSt12*MSt22*
     &                   (2*ctt2*MSt22*pMT2 + MSt2sum*(pMT2 + MUEr2))-
     &                    2*MSt12**2*
     &                   (MSt22**2 - 2*MSt22*pMT2*stt2 + 
     &                     2*pMT2*MUEr2*stt2)) + 
     &               MSb22**2*
     &                (2*MSt12*MSt22*(MSt2sum - 3*pMT2 - MUEr2) + 
     &                  pMT2*
     &                   ((-(c2tt*MSt2diff) + MSt2sum)*pMT2 + 
     &                     2*(c2tt*MSt2diff + MSt2sum)*MUEr2 - 
     &                     4*(ctt2*MSt22**2 + MSt12**2*stt2)))))/
     &           pMT2))/
     &      ((MSb22 - MSt12)*MSt12*(MSb22 - MSt22)*MSt22) + 
     &     (4*((-((1 + ctt2)*hb*ht*MUEr*s2tb) + Abr*ctb2*hb2*s2tt + 
     &             4*ctt2*ht2*pMT*stb2 + 2*Atr*ht2*s2tt*stb2 + 
     &             2*ctb2*hb2*pMT*stt2)*(-x10 + x18) + 
     &          (ctb2*hb2*(2*ctt2*pMT - Abr*s2tt) - 
     &             hb*ht*MUEr*s2tb*(1 + stt2) + 
     &             ht2*stb2*(-2*Atr*s2tt + 4*pMT*stt2))*(-x10 + x19))
     &        )/pMT + (8*(ctb2*hb2*pMT - hb*ht*MUEr*s2tb + 
     &          ht2*pMT*stb2)*(x10 - x20))/pMT)

        tmp2 = tmp2 + (ctb2*hb2 - (hb*ht*MUEr*s2tb)/pMT + 
     &      ht2*stb2)*(-x10 + x24) + 
     &   1/2.D0*((-MSb22**2 + (pMT2 - MUEr2)**2)*
     &       (ctb2*hb2*pMT*(-MSb22 + pMT2 + MUEr2) - 
     &         2*hb*ht*pMT2*MUEr*s2tb + 
     &         ht2*pMT*(-MSb22 + pMT2 + MUEr2)*stb2)*x25)/(pMT*pMT2)

        tmp3 = 2*(ctt2*hb2 + ht2*stt2)*(-x10 + x24) + 
     &   1/4.D0*(2*hb*ht*MSb2diff*MUEr*s2tb*(2*ctt2*pMT + Atr*s2tt) + 
     &       (hb2*(4*Abr*(MSb12 - MSt12)*(MSb22 - MSt12)*s2tt + 
     &            pMT*(MSb22*
     &                (4*ctt2*MUEr2*stb2 - 4*MSt12*(1 + stt2)) + 
     &               4*(MSt12*(MSt12 - ctt2*MUEr2 + MSt12*stt2) + 
     &                  MSb12*
     &                   (ctb2*ctt2*MUEr2 + 
     &                     (MSb22 - MSt12)*(1 + stt2))))))/pMT)/
     &     ((MSb12 - MSt12)*(-MSb22 + MSt12))

        tmp3 = tmp3 + ht2*
     &    ((-2 + s2tt2 + Nc*s2tt2)*(x10 - x19) + 2*(x10 - x20) + 
     &      1/4.D0*((MA02*(MSb12 - MSt12)*(-MSb22 + MSt12)*pMT2*MUEr2*
     &              ((1 + c2tt2)*MSt12**2 + MSt22**2*s2tt2) + 
     &             MSt12*(MSt12*
     &                 (MSt12*
     &                    (MUEr2*
     &                       (MSt12*
     &                       (-6*MSt22**2 + 
     &                       (1 + c2tt2)*(3*MSt22 + MSt2sum)*pMT2)+
     &                        2*MSt22**2*pMT2*s2tt2) + 
     &                      2*MSt22*
     &                       (MSt12**2*
     &                       (MSt22 + pMT2 + c2tt2*pMT2 - 
     &                       pMT2*Nc*s2tt2) + 
     &                       MSt12*
     &                       (Atr**2*(1 + c2tt2)*pMT2 - 
     &                       6*Atr*MSt22*pMT*s2tt + 
     &                       MSt22*pMT2*
     &                       (-14 - c2tt + s2tt2 + Nc*s2tt2)) + 
     &                       MSt22*pMT2*
     &                       ((5 + c2tt)*pMT2 + 
     &                       Atr*(6*pMT*s2tt + Atr*s2tt2 + 2*Atr*stt2))
     &                       )) - 
     &                   MSb22*
     &                    (MUEr2*
     &                       (MSt12*
     &                       (-6*MSt22**2 + 
     &                       (1 + c2tt2)*(3*MSt22 + MSt2sum)*pMT2)+
     &                        2*MSt22**2*pMT2*s2tt2) + 
     &                      MSt22*
     &                       (2*MSt12**2*
     &                       (MSt22 + pMT2 + c2tt2*pMT2 - 
     &                       pMT2*Nc*s2tt2) + 
     &                       2*MSt12*
     &                       (Atr**2*(1 + c2tt2)*pMT2 - 
     &                       6*Atr*MSt22*pMT*s2tt + 
     &                       MSt22*pMT2*
     &                       (-14 - c2tt + s2tt2 + Nc*s2tt2)) + 
     &                       MSt22*pMT2*
     &                       ((9 + c2tb + 2*c2tt*ctb2)*pMT2 + 
     &                       4*Atr*
     &                       ((2 + ctb2)*pMT*s2tt + 
     &                       Atr*(ctb2 + 2*ctt2)*stt2))))) + 
     &                MSb12*
     &                 (MSb22*
     &                    (2*MSt12**2*MSt22*(MSt22 + 2*pMT2) + 
     &                      2*MSt22*pMT2*
     &                       (Atr**2*
     &                      (-(c2tt2*MSt2diff) + MSt2sum) + 
     &                       4*Atr*MSt22*pMT*s2tt + 
     &                       MSt22*(4*pMT2 + MUEr2*s2tt2)) + 
     &                      MSt12*
     &                       ((1 + c2tt2)*MSt2sum*pMT2*MUEr2 - 
     &                       2*MSt22**2*
     &                       (13*pMT2 + 2*ctt2*pMT2 + 3*MUEr2 + 
     &                       6*Atr*pMT*s2tt) + 
     &                       MSt22*pMT2*
     &                       (3*(1 + c2tt2)*MUEr2 + 
     &                       2*MSt2diff*(1 + Nc)*s2tt2))) - 
     &                   MSt12*
     &                    (MUEr2*
     &                       (MSt12*
     &                       (-6*MSt22**2 + 
     &                       (1 + c2tt2)*(3*MSt22 + MSt2sum)*pMT2)+
     &                        2*MSt22**2*pMT2*s2tt2) + 
     &                      MSt22*
     &                       (2*MSt12**2*
     &                       (MSt22 + pMT2 + c2tt2*pMT2 - 
     &                       pMT2*Nc*s2tt2) + 
     &                       2*MSt12*
     &                       (Atr**2*(1 + c2tt2)*pMT2 - 
     &                       6*Atr*MSt22*pMT*s2tt + 
     &                       MSt22*pMT2*
     &                       (-14 - c2tt + s2tt2 + Nc*s2tt2)) + 
     &                       MSt22*pMT2*
     &                       (pMT2*(9 - c2tb + 2*c2tt*stb2) + 
     &                       4*Atr*
     &                       (pMT*s2tt*(2 + stb2) + 
     &                       Atr*(2*ctt2 + stb2)*stt2)))))))/
     &           (MSt12**2*MSt22**2) - 
     &          2*(MSb12 - MSt12)*(MSb22 - MSt12)*
     &           (MSt12 - pMT2 - MUEr2)*
     &           ((MSt12 - pMT2)**2 - 4*MSt12*MUEr2 + 3*MUEr2**2)*x26)
     &         /((-MSb12 + MSt12)*(-MSb22 + MSt12)*pMT2))

        F1t = hb2 + tmp3*(-x10 + x18) + tmp1*(-x10 + x21) + 
     &   tmp2*(-x10 + x22) - (hb2 + ht2)*(x10 - x24)**2 + 
     &   (-x10 + x20)*(3*hb2 + (hb2 + ht2)*(x10 - x24) + 
     &      ht2*(10 - 1/2.D0*(MSt2sum*(3*pMT2 + MUEr2))/(MSt12*MSt22)))
     &     + (x10 - x22)**2*
     &    (-(ht2*stb2) + 1/2.D0*(hb*ht*MUEr*s2tb)/pMT) - 
     &   1/2.D0*((2*ctb2*ht2*pMT + hb*ht*MUEr*s2tb)*(x10 - x21)**2)/pMT

        F1t = F1t + (-x10 + x19)*
     &    (2*(ctt2*ht2 + hb2*stt2)*(-x10 + x24) + 
     &      ((hb2*((MSb12 - MSt22)*(MSb22 - MSt22)*
     &               (pMT + ctt2*pMT - Abr*s2tt) + 
     &              pMT*MUEr2*(ctb2*MSb12 - MSt22 + MSb22*stb2)*stt2)
     &            )/pMT - 1/2.D0*
     &          (hb*ht*MSb2diff*MUEr*s2tb*(Atr*s2tt - 2*pMT*stt2)))/
     &       ((MSb12 - MSt22)*(-MSb22 + MSt22)) + 
     &      1/4.D0*(ht2*(4*ctt2 - 
     &           ((MSt12 - MSt22)*(MUEr2 + c2tt2*(4*MSt22 + MUEr2)))/
     &            (MSt12*MSt2diff) + 
     &           (2*Atr*(6*MSt22**2 + c2tb*MSb2diff*pMT2 - 
     &                MSb2sum*pMT2 - 6*MSt22*pMT2)*s2tt)/
     &            (MSt22**2*pMT) - 
     &           (2*((-5 + c2tt)*pMT2 + 2*c2tt2*MUEr2 - 
     &                Atr**2*(2*ctt2 + s2tt2) + 
     &                (2*ctb2*MSb12**2*
     &                   (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))/
     &                 (MSb12 - MSt22)**2 + 
     &                (2*MSb22**2*stb2*
     &                   (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))/
     &                 (MSb22 - MSt22)**2))/MSt22 + 
     &           8*(x10 - x20) + 
     &           2*(-15 + s2tt2 + Nc*s2tt2 + 
     &              ((MSt12 - MSt22)*
     &                 (-(Atr**2*(1 + c2tt2)) + 
     &                   MSt22*(-1 + Nc)*s2tt2))/(MSt12*MSt2diff)-
     &                (4*MSb22**2*stb2*
     &                 (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))/
     &               ((MSb22 - MSt22)*MSt22**2) + 
     &              (MSt22 - 3*MUEr2 + 
     &                 (-MSt22 + pMT2 + MUEr2)*
     &                  ((MSt22 - pMT2)**2 - 4*MSt22*MUEr2 + 
     &                    3*MUEr2**2)*x27)/pMT2) + 
     &           (4*(-((ctb2*MSb12**2*(MSb12 - 2*MSt22)*
     &                     (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))/
     &                   (MSb12 - MSt22)**2) + 
     &                (MSb22**3*stb2*
     &                   (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))/
     &                 (MSb22 - MSt22)**2 + 
     &                1/4.D0*
     &                 (-(MA02*MUEr2*s2tt2) - 
     &                   2*(c2tb*MSb2diff - MSb2sum)*
     &                    (Atr**2*ctt2 + pMT2*stt2))))/MSt22**2)))

        F1t = F1t + 1/4.D0*
     &    (((MA02 + MSb12 - MSt12)*x11*
     &         ((2*Abr**2*ctt2*hb2*stb2 + 
     &              Abr*s2tt*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2) + 
     &              2*(ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb + 
     &                 hb2*pMT2*stb2)*stt2)*
     &            (MA02*(x10 - x21) + MSb12*(-x10 + x21) + 
     &              MSt12*(x10 - 2*x18 + x21)) - 
     &           ((MA02 + MSb12 - MSt12)*
     &              (Abr*(hb*ht*pMT*MUEr*s2tb*s2tt + 
     &                   2*hb2*(Abr*ctt2*pMT + pMT2*s2tt)*stb2) + 
     &                2*(ctb2*ht2*pMT*MUEr2 + hb*ht*pMT2*MUEr*s2tb + 
     &                   hb2*pMT*pMT2*stb2)*stt2)*x28)/pMT))/MSt12 + 
     &      (MA02 + MSb12 - MSt22)*x12*
     &       (-2*(2*ctt2*(ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb + 
     &               hb2*pMT2*stb2) + 
     &            Abr*(-(hb*ht*MUEr*s2tb*s2tt) + 
     &               2*hb2*stb2*(-(pMT*s2tt) + Abr*stt2)))*
     &          (-x10 + x19) + 
     &         ((MA02 - MSb12 - MSt22)*
     &            (2*ctb2*ctt2*ht2*MUEr2 + 
     &              2*ctt2*(hb*ht*pMT*MUEr*s2tb + hb2*pMT2*stb2) + 
     &              Abr*(-(hb*ht*MUEr*s2tb*s2tt) - 
     &                 2*hb2*pMT*s2tt*stb2 + 2*Abr*hb2*stb2*stt2))*
     &            (x10 - x21))/MSt22 + 
     &         ((MA02 + MSb12 - MSt22)*
     &            (-2*ctb2*ctt2*ht2*pMT*MUEr2 - 
     &              2*ctt2*pMT2*(hb*ht*MUEr*s2tb + hb2*pMT*stb2) + 
     &              Abr*(hb*ht*pMT*MUEr*s2tb*s2tt + 
     &                 2*hb2*pMT2*s2tt*stb2 - 2*Abr*hb2*pMT*stb2*stt2)
     &              )*x29)/(MSt22*pMT)))

        F1t = F1t + 1/4.D0*
     &    ((2*ht2*MSt22*(-1 + Nc)*s2tt2 - 
     &         (c2tb*MSb2diff*
     &             (c2tt*(hb2 + ht2)*MSt2diff + 
     &               (hb2 - ht2)*MSt2sum) + 
     &            MSb2sum*(c2tt*(hb2 - ht2)*MSt2diff + 
     &               (hb2 + ht2)*MSt2sum) + 
     &            2*(c2tt*hb2*MA02*MSt2diff + hb2*MA02*MSt2sum + 
     &               Abr**2*hb2*(c2tt*MSt2diff + MSt2sum) + 
     &               Atr**2*
     &                (-(c2tt*ht2*MSt2diff) + 3*ht2*MSt2sum) - 
     &               c2tt*hb2*MSt2diff*pMT2 + 
     &               c2tt*ht2*MSt2diff*pMT2 + hb2*MSt2sum*pMT2 + 
     &               ht2*MSt2sum*pMT2 - c2tt*hb2*MSt2diff*MUEr2 + 
     &               c2tt*ht2*MSt2diff*MUEr2 - hb2*MSt2sum*MUEr2 - 
     &               3*ht2*MSt2sum*MUEr2) + 
     &            4*(Abr*hb2 + 3*Atr*ht2)*MSt2diff*pMT*s2tt - 
     &            ((-2*ht2*
     &                  (pMT2*
     &                     (MSt2sum*pMT2 + (5 - 2*c2tt)*MSt22*MUEr2)+
     &                      MSt12*
     &                     ((5 + 2*c2tt)*pMT2*MUEr2 - 
     &                       MSt22*
     &                       (ctb2*MSb12 - MSt2sum + 3*pMT2 + 
     &                       7*MUEr2 + MSb22*stb2))) + 
     &                 hb2*
     &                  (MSt22*
     &                     (MSt12*
     &                       (c2tb*MSb2diff + MSb2sum - 2*pMT2) + 
     &                       2*(MSt12 - 4*ctt2*pMT2)*MUEr2) - 
     &                    8*MSt12*pMT2*MUEr2*stt2))*(-x10 + x24))/pMT2
     &            )/MSt22 + 
     &         ((2*Abr**2*ctt2*hb2*pMT*stb2 + 
     &              Abr*s2tt*
     &               (hb*ht*pMT*MUEr*s2tb + 
     &                 2*hb2*(-MA02 - MSb12 + MSt12 + pMT2)*stb2) + 
     &              2*(ctb2*ht2*pMT*MUEr2 + 
     &                 hb*ht*(-MA02 - MSb12 + MSt12 + pMT2)*MUEr*
     &                  s2tb + 
     &                 hb2*pMT*(-2*MA02 - 2*MSb12 + 2*MSt12 + pMT2)*
     &                  stb2)*stt2)*x28)/pMT + 
     &         ((2*Abr**2*ctb2*ctt2*hb2*pMT - 
     &              Abr*(2*ctb2*hb2*(MA02 + MSb22 - MSt12 - pMT2) + 
     &                 hb*ht*pMT*MUEr*s2tb)*s2tt + 
     &              2*(ctb2*hb2*pMT*
     &                  (-2*MA02 - 2*MSb22 + 2*MSt12 + pMT2) + 
     &                 hb*ht*(MA02 + MSb22 - MSt12 - pMT2)*MUEr*
     &                  s2tb + ht2*pMT*MUEr2*stb2)*stt2)*x30)/pMT)/
     &       MSt12 + (MA02 + MSb22 - MSt12)*x13*
     &       (2*(2*Abr**2*ctb2*ctt2*hb2 + 
     &            Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt + 
     &            2*(ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + 
     &               ht2*MUEr2*stb2)*stt2)*(x10 - x18) + 
     &         ((MA02 - MSb22 - MSt12)*
     &             (hb*ht*MUEr*s2tb*(Abr*s2tt + 2*pMT*stt2) - 
     &               2*(Abr**2*ctb2*ctt2*hb2 + 
     &                  Abr*ctb2*hb2*pMT*s2tt + ctb2*hb2*pMT2*stt2 + 
     &                  ht2*MUEr2*stb2*stt2))*(-x10 + x22) - 
     &            ((MA02 + MSb22 - MSt12)*
     &               (Abr*(-(hb*ht*pMT*MUEr*s2tb*s2tt) + 
     &                    2*ctb2*hb2*(Abr*ctt2*pMT + pMT2*s2tt)) + 
     &                 2*(ctb2*hb2*pMT*pMT2 - hb*ht*pMT2*MUEr*s2tb + 
     &                    ht2*pMT*MUEr2*stb2)*stt2)*x30)/pMT)/MSt12))

        F1t = F1t + 1/4.D0*
     &    ((MA02 + MSb22 - MSt22)*x14*
     &      (-2*(hb*ht*MUEr*s2tb*(-2*ctt2*pMT + Abr*s2tt) + 
     &           2*(ctt2*ht2*MUEr2*stb2 + 
     &              ctb2*hb2*(ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))
     &              ))*(-x10 + x19) + 
     &        ((MA02 - MSb22 - MSt22)*
     &           (-2*ctt2*hb*ht*pMT*MUEr*s2tb + 
     &             Abr*hb*ht*MUEr*s2tb*s2tt + 2*ctt2*ht2*MUEr2*stb2 + 
     &             2*ctb2*hb2*
     &              (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*
     &           (x10 - x22))/MSt22 - 
     &        ((MA02 + MSb22 - MSt22)*
     &           (-2*ctt2*hb*ht*pMT2*MUEr*s2tb + 
     &             Abr*hb*ht*pMT*MUEr*s2tb*s2tt + 
     &             2*ctt2*ht2*pMT*MUEr2*stb2 + 
     &             2*ctb2*hb2*
     &              (ctt2*pMT*pMT2 + Abr*(-(pMT2*s2tt) + Abr*pMT*stt2)))*
     &           x31)/(MSt22*pMT))) + 
     &   1/4.D0*(Abr*(2*hb2*stb2*
     &           (MA02*s2tt + MSb12*s2tt - MSt22*s2tt - pMT2*s2tt + 
     &             Abr*pMT*stt2)*x29 + 
     &          hb*ht*pMT*MUEr*s2tb*s2tt*(-x29 + x31)) + 
     &       2*ctb2*(Abr*hb2*
     &           (MA02*s2tt + MSb22*s2tt - MSt22*s2tt - pMT2*s2tt + 
     &             Abr*pMT*stt2)*x31 + 
     &          ctt2*pMT*(ht2*MUEr2*x29 + 
     &             hb2*(-2*MA02 - 2*MSb22 + 2*MSt22 + pMT2)*x31)) - 
     &       2*ctt2*(hb*ht*MUEr*s2tb*
     &           (MSb12*x29 - MSt22*x29 - pMT2*x29 + 
     &             MA02*(x29 - x31) - MSb22*x31 + MSt22*x31 + 
     &             pMT2*x31) + 
     &          pMT*stb2*(hb2*(2*MA02 + 2*MSb12 - 2*MSt22 - pMT2)*
     &              x29 - ht2*MUEr2*x31)))/(MSt22*pMT) + 
     &   x23*(-((MSb12 - pMT2 + MUEr2)*
     &         (-2*hb*ht*pMT*MUEr*s2tb + 
     &           (MSb12 - pMT2 - MUEr2)*(ctb2*ht2 + hb2*stb2))*
     &         (-x10 + x20)) + 
     &      1/2.D0*((((MSb12 - pMT2)**2 - MUEr2**2)*
     &             (2*hb*ht*pMT2*MUEr*s2tb + 
     &               pMT*(-MSb12 + pMT2 + MUEr2)*(ctb2*ht2 + hb2*stb2)
     &               )*(-x10 + x24))/pMT2 + 
     &          ((MSb12 - pMT2 + MUEr2)**2*
     &             (2*hb*ht*pMT2*MUEr2*s2tb + 
     &               pMT*MUEr*(-MSb12 + pMT2 + MUEr2)*
     &                (ctb2*ht2 + hb2*stb2))*x32)/(MUEr*MUEr2))/pMT)+
     &     (-((ctt2*hb*ht*MUEr*s2tb + 
     &           ctb2*ht2*(2*ctt2*pMT + Atr*s2tt))*x33) + 
     &      (ctt2*hb*ht*MUEr*s2tb - ht2*(2*ctt2*pMT + Atr*s2tt)*stb2)*
     &       x34 + (Atr*ctb2*ht2*s2tt - 
     &         (2*ctb2*ht2*pMT + hb*ht*MUEr*s2tb)*stt2)*x35 + 
     &      (Atr*ht2*s2tt*stb2 + hb*ht*MUEr*s2tb*stt2 - 
     &         2*ht2*pMT*stb2*stt2)*x36 + 
     &      1/2.D0*((2*hb*ht*MUEr2*(MSb12 - 2*pMT2 + MUEr2)*s2tb + 
     &             pMT*MUEr*(3*MSb12 - 3*pMT2 + MUEr2)*
     &              (ctb2*ht2 + hb2*stb2))*x32 + 
     &          (-2*hb*ht*MUEr2*(MSb22 - 2*pMT2 + MUEr2)*s2tb + 
     &             pMT*MUEr*(3*MSb22 - 3*pMT2 + MUEr2)*
     &              (ctb2*hb2 + ht2*stb2))*x37)/(MUEr*MUEr2))/pMT

        F1t = F1t + x25*(-((MSb22 - pMT2 + MUEr2)*
     &         (2*hb*ht*pMT*MUEr*s2tb + 
     &           (MSb22 - pMT2 - MUEr2)*(ctb2*hb2 + ht2*stb2))*
     &         (-x10 + x20)) + 
     &      1/2.D0*((((MSb22 - pMT2)**2 - MUEr2**2)*
     &             (-2*hb*ht*pMT2*MUEr*s2tb + 
     &               pMT*(-MSb22 + pMT2 + MUEr2)*(ctb2*hb2 + ht2*stb2)
     &               )*(-x10 + x24))/pMT2 + 
     &          ((MSb22 - pMT2 + MUEr2)**2*
     &             (-2*hb*ht*pMT2*MUEr2*s2tb + 
     &               pMT*MUEr*(-MSb22 + pMT2 + MUEr2)*
     &                (ctb2*hb2 + ht2*stb2))*x37)/(MUEr*MUEr2))/pMT)

        F1t = F1t + 1/24.D0*
     &    (ht2*(72 - (24*c2tt2*MSt12)/MSt22 - 
     &        (24*c2tt2*MSt22)/MSt12 - 8*Pi**2 - 24*s2tt2 + 
     &        (12*MSt12*(-1 + Nc)*s2tt2)/MSt22 - 24*Nc*s2tt2 - 
     &        72*(x10 - x20)**2 + 
     &        (6*MUEr2*s2tt2*x38)/(MSt12**3*x16) + 
     &        (6*(1 + c2tt2)*
     &           (-3*MA02**2 + 6*MA02*MSt12 + 8*MA02*MSt22 + 
     &             3*MSt12*MSt2diff - MSt22*MSt2diff)*MUEr2*x39)/
     &         MSt22**3 + 12*
     &         ((4 + (2*Atr*s2tt)/pMT + s2tt2 + Nc*s2tt2)*
     &            (x10 - x18)**2 + 
     &           (4 - (2*Atr*s2tt)/pMT + s2tt2 + Nc*s2tt2)*
     &            (x10 - x19)**2 + 
     &           ((1 + c2tt2)*MUEr2*x39)/(MSt22**3*x15)) + 
     &        (6*(1 + c2tt2)*MUEr2*x15*
     &           (-((MSt22**2*
     &                  (-((MSt12 - MSt22)*
     &                       (MSt22**3 + 
     &                       MSt12**2*(3*MSt22 - 2*MSt2diff) - 
     &                       4*MSt12*MSt22*(MSt22 - MSt2diff))) - 
     &                    4*MA02*MSt12*(MSt12 + 2*MSt22)*
     &                     MSt2diff + 
     &                    MA02**2*
     &                     (-MSt22**2 + MSt12*(MSt22 + 2*MSt2diff))
     &                    )*(x10 - x19))/(MSt12*MSt2diff)) + 
     &             MA02**4*x39 - 
     &             MSt2diff**2*
     &              (MSt22**2*(-x10 + x18) + 
     &                MSt22*MSt2sum*(-x10 + x18) + 
     &                MSt2diff*MSt2sum*x39) + 
     &             MA02**2*
     &              (3*MSt22**2*(x10 - x18) + 
     &                6*MSt12*MSt2sum*x39 + 
     &                MSt22*MSt2sum*(3*x10 - 3*x18 + 4*x39)) - 
     &             MA02**3*
     &              (4*MSt12*x39 + MSt22*(x10 - x18 + 6*x39)) + 
     &             MA02*(4*MSt12**2*MSt2diff*x39 - 
     &                2*MSt22**2*MSt2diff*x39 + 
     &                MSt22*
     &                 (-3*MSt2sum**2*(x10 - x18) + 
     &                   6*MSt12*MSt2diff*x39))))/MSt22**3 + 
     &        (6*MUEr2*s2tt2*x40)/(MSt22**3*x17) + 
     &        3*MA02*MUEr2*s2tt2*
     &         ((-3*(MA02 - 4*MSt12)*x38 + 
     &              MA02*(MA02 - 6*MSt12)*x16*
     &               (MA02*x38 - 2*MSt12*(x10 - x18 + x38)))/
     &            MSt12**3 + 
     &           (-3*(MA02 - 4*MSt22)*x40 + 
     &              MA02*(MA02 - 6*MSt22)*x17*
     &               (MA02*x40 - 2*MSt22*(x10 - x19 + x40)))/
     &            MSt22**3) + 
     &        (12*x26*((MSt12 - pMT2 - MUEr2)*
     &              ((MSt12 - pMT2)**2 - 4*MSt12*MUEr2 - MUEr2**2)*
     &              (-x10 + x20) - 
     &             ((-MSt12 + pMT2 + MUEr2)**2*
     &                (-(MSt12 - pMT2)**2 + (3*MSt12 + pMT2)*MUEr2)*
     &                (-x10 + x24))/pMT2 + 
     &             4*MUEr2**2*(-MSt12 + pMT2 + MUEr2)*x41))/MSt12 + 
     &        (12*x27*((MSt22 - pMT2 - MUEr2)*
     &              ((MSt22 - pMT2)**2 - 4*MSt22*MUEr2 - MUEr2**2)*
     &              (-x10 + x20) - 
     &             ((-MSt22 + pMT2 + MUEr2)**2*
     &                (-(MSt22 - pMT2)**2 + (3*MSt22 + pMT2)*MUEr2)*
     &                (-x10 + x24))/pMT2 + 
     &             4*MUEr2**2*(-MSt22 + pMT2 + MUEr2)*x42))/MSt22))

        F1t = F1t - 2*(ctt2*hb2*x6 + ht2*stt2*x6 + ctt2*ht2*x7 + 
     &      hb2*stt2*x7) + hb2*(2*x8 + (x10 - x9)**2) + 
     &   1/4.D0*((-x10 + x9)*
     &      ((-((-2 + c2tt - c2tt2)*ht2*MUEr2) + 
     &           (MA02 + MSb12 - MSt12)*(MA02 - MSb12 + MSt12)*
     &            (hb*ht*MUEr*s2tb*(Abr*s2tt + 2*pMT*stt2) + 
     &              2*(ctb2*ht2*MUEr2*stt2 + 
     &                 hb2*stb2*
     &                  (Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2)))*x11)
     &          /MSt12 + ((MA02 + MSb12 - MSt22)*
     &           (MA02 - MSb12 + MSt22)*
     &           (hb*ht*MUEr*s2tb*(2*ctt2*pMT - Abr*s2tt) + 
     &             2*(ctb2*ctt2*ht2*MUEr2 + 
     &                hb2*stb2*
     &                 (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))))*x12
     &           )/MSt22 + 
     &        ((MA02**2 - (MSb22 - MSt12)**2)*
     &           (2*Abr**2*ctb2*ctt2*hb2 + 
     &             Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt + 
     &             2*(ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + 
     &                ht2*MUEr2*stb2)*stt2)*x13)/MSt12 + 
     &        ((MA02**2 - (MSb22 - MSt22)**2)*
     &           (-2*ctt2*hb*ht*pMT*MUEr*s2tb + 
     &             Abr*hb*ht*MUEr*s2tb*s2tt + 2*ctt2*ht2*MUEr2*stb2 + 
     &             2*ctb2*hb2*
     &              (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*x14)/
     &         MSt22 + ht2*MUEr2*
     &         ((s2tt2*(MA02 - MA02**2*(MA02 - 6*MSt12)*x16))/
     &            MSt12**2 - 
     &           ((1 + c2tt2)*MSt12 + 
     &              MA02**2*(MA02 - 6*MSt22)*s2tt2*x17)/MSt22**2)+
     &          (2*(ht2*(MA02 + ctt2*MSt22)*pMT*MUEr2 + 
     &             hb2*MSt22**2*
     &              (2*pMT*stt2*(-x10 + x18) + 
     &                Abr*s2tt*(x18 - x19) + 2*ctt2*pMT*(-x10 + x19))
     &             ))/(MSt22**2*pMT) + 
     &        (hb2*(4*ctt2*MA02*MSt22 - c2tt*MSt2diff*pMT2 + 
     &             MSt2sum*pMT2 + 2*Abr*MSt2diff*pMT*s2tt + 
     &             4*MA02*MSt12*stt2 + 
     &             2*Abr**2*(ctt2*MSt22 + MSt12*stt2) + 
     &             8*MSt12*MSt22*x10 - 8*MSt12*MSt22*x20))/
     &         (MSt12*MSt22) + 
     &        (2*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2)*(x10 - x21))/
     &         pMT + (2*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*
     &           (x10 - x22))/pMT + 
     &        (1 + c2tt2)*ht2*MUEr2*x15*
     &         (5*MA02 - (MA02 - MSt12)**3/MSt22**2 + 2*MSt22 + 
     &           ((MA02 - MSt22)*(MA02 + MSt22))/MSt12 + 
     &           (4*(MA02**2 - 1/2.D0*(MSt12*(MA02 + MSt12))))/MSt22)
     &        ))

        tmp4 = s2tt*(ht2*
     &       (stt2 + ((ctb2*MSb12**2*
     &               (Atr**2 + 2*ctt2*pMT2 + 2*Atr*pMT*s2tt))/
     &             (-MSb12 + MSt12) + 
     &            (MSb22**3*stb2*
     &               (ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2)))/
     &             (MSb22 - MSt12)**2)/MSt12**2) + 
     &      (hb2 + ht2 + (c2tt*(hb2 - ht2)*(MSt2sum - 2*MUEr2))/
     &          MSt2diff)*(-x10 + x24) + 
     &      1/2.D0*((hb*ht*MSb2diff*MUEr*s2tb*(2*ctt2*pMT + Atr*s2tt))/
     &          ((MSb12 - MSt12)*(-MSb22 + MSt12)) + 
     &         (ht2*(-(Atr**2*(-3 + c2tt)) + (5 + c2tt)*pMT2 + 
     &              MUEr2*s2tt2))/MSt12)) + 
     &   (ht2*(c2tt2*(6*Atr*pMT - 2*MSt22*s2tt) + 
     &         s2tt*((c2tt + c2tt2)*MSt12 + 
     &            c2tt*(Atr**2 - pMT2 - MUEr2) + 
     &            2*c2tt2*(Atr**2 + MSt2sum + MUEr2 + MSt2sum*Nc)*
     &             (-x10 + x19))) - 
     &      1/4.D0*(hb2*(2*(MSb12 - MSt12)*(MSb22 - MSt12)*
     &             (2*Abr**2*c2tt*pMT*s2tt + 2*c2tt*MA02*pMT*s2tt + 
     &               Abr*(-4*c2tt2*pMT2 + (2*MSt12 + MSt2sum)*s2tt2))
     &              + pMT*s2tt*
     &             (c2tt*(-4*ctb2*MSb22**2*MSt12 + 
     &                  4*MSb12**2*(MSb22 - MSt12)*stb2) + 
     &               MSb22*
     &                ((-4*c2tt*MSt12**2 + 4*ctt2*MSt22*MUEr2)*
     &                   stb2 + 
     &                  MSt12*
     &                   (-4*MSt2diff + 
     &                     2*c2tt*(2*pMT2 + MUEr2 + ctb2*MUEr2) - 
     &                     2*MUEr2*stb2)) + 
     &               MSb12*
     &                (4*ctb2*(c2tt*MSb22**2 + ctt2*MSt22*MUEr2) + 
     &                  MSb22*(4*MSt2diff - 4*c2tt*(pMT2 + MUEr2)) + 
     &                  MSt12*
     &                   (-4*MSt22 - 2*ctb2*MUEr2 + 
     &                     2*c2tt*(2*pMT2 + MUEr2 + MUEr2*stb2)) + 
     &                  MSt12**2*(2 - 2*c2tb*c2tt + 4*stt2)) + 
     &               4*MSt12*
     &                (-(ctt2*MSt22*MUEr2) + 
     &                  MSt12*
     &                   (MSt2diff + c2tt*(MSt12 - pMT2) + 
     &                     MUEr2*stt2)))))/
     &        ((MSb12 - MSt12)*(MSb22 - MSt12)*pMT))/MSt2diff

        tmp4 = tmp4 + ht2*
     &    (Atr*(((-3/2.D0 - (6*MSt12)/MSt2diff + (3*pMT2)/MSt12)*
     &            s2tt2)/pMT + 
     &         1/2.D0*(Atr*(1 + c2tt2)*MSt22*s2tt + 
     &             (-(c2tb*MSb2diff) + MSb2sum)*pMT*s2tt2)/MSt12**2)
     &        + 1/4.D0*(s2tt*
     &         (-20 + (2*Atr**2*
     &              ((2*c2tt*ctb2*MSb12**2)/(MSb12 - MSt12) + 
     &                ((1 + c2tt2)*MSt22**2*
     &                   (-MSt12 + MSt22 - 2*MSt2diff))/MSt2diff**2
     &                ))/MSt12**2 - 8*c2tt2*Nc + 
     &           (4*ctb2*MSb12**2*
     &              (ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2)))/
     &            ((MSb12 - MSt12)*MSt12**2) - 
     &           (4*MSb22**2*stb2*
     &              (ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2)))/
     &            ((MSb22 - MSt12)**2*MSt12) - 
     &           (MA02*MUEr2*s2tt2 + 
     &              2*(c2tb*MSb2diff - MSb2sum)*
     &               (ctt2*pMT2 + Atr**2*stt2))/MSt12**2 + 
     &           (2*(-(MSb22*MSt12**2*
     &                   (MSt12 + c2tt*pMT2 - MUEr2 + pMT2*s2tt2 + 
     &                     pMT2*Nc*s2tt2)) + 
     &                MSt12**3*
     &                 (MSt12 + c2tt*pMT2 - MUEr2 + pMT2*s2tt2 + 
     &                   pMT2*Nc*s2tt2) + 
     &                4*MSb22**2*pMT2*stb2*
     &                 (ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2))))/
     &            (MSt12**2*(-MSb22 + MSt12)*pMT2) + 
     &           (2*(-MSt12 + pMT2 + MUEr2)**3*x26)/pMT2 + 
     &           (4*((c2tt2*MSt12**2 + 
     &                   1/4.D0*((1 + c2tt2)*(MA02 + MSt12)*MUEr2))/
     &                 MSt22 + 
     &                1/4.D0*
     &                 ((1 + c2tt2)*MSt12*(-MA02 + MSt12)*MUEr2)/
     &                  MSt22**2 + 
     &                1/2.D0*
     &                 (c2tt*(-(c2tb*MSb2diff) + MSb2sum) - 
     &                   (1 + c2tt2)*MUEr2 - 
     &                   (MSt12*
     &                      (-(Atr**2*(1 + c2tt2)) + 
     &                       (MSt12 - MSt22)*(-1 + Nc)*s2tt2))/
     &                    MSt22)))/MSt2diff)))

        tmp5 = s2tt*((-hb2 - ht2 - 
     &         (c2tt*(hb2 - ht2)*(MSt2sum - 2*MUEr2))/MSt2diff)*
     &       (-x10 + x24) + 
     &      1/2.D0*(ht2*(c2tt - (1 + Nc)*s2tt2 + 
     &            (-(Atr**2*(1 + c2tt2)) + c2tb*c2tt*MSb2diff - 
     &               c2tt*MSb2sum + MSt22*(-1 + Nc)*s2tt2)/MSt2diff
     &            ) + (hb*ht*MSb2diff*MUEr*s2tb*
     &            (Atr*s2tt - 2*pMT*stt2))/
     &          ((MSb12 - MSt22)*(-MSb22 + MSt22)))) + 
     &   1/4.D0*((ht2*s2tt*(2*Atr**2*ctt2*
     &            (2*ctb2*MSb12**3 + 
     &              2*MSb12*(-(c2tb*MSb2diff) + MSb2sum)*MSt22 + 
     &              (c2tb*MSb2diff - MSb2sum)*MSt22**2 - 
     &              MSb12**2*
     &               (-(c2tb*MSb2diff) + MSb2sum + 2*MSt22 + 
     &                 2*c2tb*MSt22)) - 
     &           4*Atr*ctb2*MSb12**2*(MSb12 - 2*MSt22)*pMT*s2tt + 
     &           MA02*(MSb12 - MSt22)**2*MUEr2*s2tt2 + 
     &           2*(2*ctb2*MSb12**2*(MSb12 - 2*MSt22) + 
     &              (c2tb*MSb2diff - MSb2sum)*(MSb12 - MSt22)**2)*
     &            pMT2*stt2))/((MSb12 - MSt22)**2*MSt22**2) + 
     &      (hb2*((MSb12 - MSt22)*
     &            (4*c2tt*(Abr**2 + MA02)*(MSb22 - MSt22)*pMT*s2tt + 
     &              2*Abr*(-MSb22 + MSt22)*
     &               (4*c2tt2*pMT2 + (-4*MSt22 + MSt2diff)*s2tt2))+
     &             pMT*s2tt*
     &            (c2tt*(-4*ctb2*MSb22**2*MSt22 + 
     &                 4*MSb12**2*(MSb22 - MSt22)*stb2) + 
     &              4*MSt22*
     &               (MSt22*
     &                  (MSt2diff + c2tt*(MSt22 - pMT2) - ctt2*MUEr2)
     &                   + MSt12*MUEr2*stt2) + 
     &              MSb22*(MSt22*
     &                  (4*MSt12 + 4*c2tt*pMT2 + MUEr2 + 
     &                    3*c2tt*MUEr2 - 4*MSt22*(1 + c2tt*stb2)) - 
     &                 2*MUEr2*(c2tb*MSt22 + 2*MSt12*stb2)*stt2) + 
     &              MSb12*(-2*MSt22*
     &                  (MSt22 + 2*ctt2*MSt22 - ctb2*MUEr2) + 
     &                 2*c2tt*
     &                  (2*ctb2*MSb22**2 - 2*MSb22*(pMT2 + MUEr2) + 
     &                    MSt22*
     &                     (-(c2tb*MSt22) + 2*pMT2 + MUEr2 + 
     &                       MUEr2*stb2)) + 
     &                 4*(MSb22*MSt2diff + 
     &                    MSt12*(MSt22 - ctb2*MUEr2*stt2))))))/
     &       ((-MSb12 + MSt22)*(-MSb22 + MSt22)*MSt2diff*pMT))

        tmp5 = tmp5 + ht2*
     &    ((Atr*s2tt2*(-3/2.D0 + (6*MSt22)/MSt2diff + 
     &           1/2.D0*((-(c2tb*MSb2diff) + MSb2sum + 6*MSt22)*pMT2)/
     &             MSt22**2))/pMT + 
     &      (-6*Atr*c2tt2*pMT + 
     &         s2tt*(c2tt2*(-2*MSt12 + MSt22) + 
     &            c2tt*(-Atr**2 - MSt22 + pMT2 + MUEr2) - 
     &            1/4.D0*((1 + c2tt2)*(MSt12 - MSt22)*MUEr2)/MSt12 + 
     &            1/2.D0*(MSt22*
     &                (Atr**2*(1 + c2tt2) - MSt22*(-1 + Nc)*s2tt2))/
     &              MSt12))/MSt2diff + 
     &      s2tt*(5 - ctt2 + 
     &         c2tt2*(MSt22**2/(MSt12*MSt2diff) + 2*Nc) + 
     &         (MSb22**2*stb2*
     &            (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))/
     &          ((MSb22 - MSt22)*MSt22**2) + 
     &         1/2.D0*(-MSt22 + MUEr2)/pMT2 + 
     &         ((ctb2*MSb12**2*
     &               (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))/
     &             (MSb12 - MSt22)**2 + 
     &            1/2.D0*((-5 + c2tt)*pMT2 + 2*c2tt2*MUEr2 - 
     &               Atr**2*(2*ctt2 + s2tt2)))/MSt22 - 
     &         1/2.D0*((-MSt22 + pMT2 + MUEr2)**3*x27)/pMT2))

        tmp6 = x33*(-(c2tt*hb2*MUEr2*s2tt*stb2) + 
     &      (-(1/4.D0*(hb*ht*MUEr*s2tb*
     &              (-4*Atr*c2tt2*pMT + 
     &                (MSt2diff + 
     &                   c2tt*(-2*MSb12 + MSt2sum + 4*pMT2))*s2tt)))
     &           - 1/2.D0*(ctb2*ht2*
     &            (pMT*(MSt2diff + 
     &                 c2tt*(-2*(Atr**2 + MSb12) + MSt2sum + 2*pMT2))
     &                *s2tt + 
     &              Atr*(-4*c2tt2*pMT2 + (-2*MSb12 + MSt2sum)*s2tt2))
     &            ))/pMT) + 
     &   x35*(c2tt*hb2*MUEr2*s2tt*stb2 + 
     &      1/4.D0*(hb*ht*MUEr*s2tb*
     &           (-4*Atr*c2tt2*pMT + 
     &             (MSt2diff + c2tt*(-2*MSb12 + MSt2sum + 4*pMT2))*
     &              s2tt) + 
     &          2*ctb2*ht2*
     &           (pMT*(MSt2diff + 
     &                c2tt*(-2*(Atr**2 + MSb12) + MSt2sum + 2*pMT2))*
     &              s2tt + 
     &             Atr*(-4*c2tt2*pMT2 + (-2*MSb12 + MSt2sum)*s2tt2)))
     &         /pMT) + x34*(-(c2tt*ctb2*hb2*MUEr2*s2tt) + 
     &      1/4.D0*(hb*ht*MUEr*s2tb*
     &           (-4*Atr*c2tt2*pMT + 
     &             (MSt2diff + c2tt*(-2*MSb22 + MSt2sum + 4*pMT2))*
     &              s2tt) - 
     &          2*ht2*(pMT*(MSt2diff + 
     &                c2tt*(-2*(Atr**2 + MSb22) + MSt2sum + 2*pMT2))*
     &              s2tt + 
     &             Atr*(-4*c2tt2*pMT2 + (-2*MSb22 + MSt2sum)*s2tt2))*
     &           stb2)/pMT) - 
     &   1/4.D0*((Abr**2*hb2*(MSt2diff + 
     &             c2tt*(-2*MA02 - 2*MSb12 + MSt2sum))*pMT*s2tt*stb2
     &            + Abr*(MSt2diff*s2tt2*
     &              (-(hb*ht*pMT*MUEr*s2tb) + 
     &                hb2*(MA02 + MSb12 - MSt22 - 2*pMT2)*stb2) + 
     &             2*c2tt2*(MA02 + MSb12 - MSt22)*
     &              (hb*ht*pMT*MUEr*s2tb + 2*hb2*pMT2*stb2)) + 
     &          s2tt*(c2tt*(2*MA02 + 2*MSb12 - MSt2sum)*
     &              (ctb2*ht2*pMT*MUEr2 + hb*ht*pMT2*MUEr*s2tb + 
     &                hb2*pMT*pMT2*stb2) + 
     &             MSt2diff*
     &              (ctb2*ht2*pMT*MUEr2 + 
     &                pMT2*(hb*ht*MUEr*s2tb + hb2*pMT*stb2) - 
     &                ctt2*(MA02 + MSb12 - MSt22)*
     &                 (hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2))))*x29 + 
     &       (Abr**2*ctb2*hb2*
     &           (MSt2diff + c2tt*(-2*MA02 - 2*MSb22 + MSt2sum))*
     &           pMT*s2tt + 
     &          Abr*(2*c2tt2*(MA02 + MSb22 - MSt22)*
     &              (2*ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb) + 
     &             MSt2diff*
     &              (ctb2*hb2*(MA02 + MSb22 - MSt22 - 2*pMT2) + 
     &                hb*ht*pMT*MUEr*s2tb)*s2tt2) + 
     &          s2tt*(ctb2*hb2*pMT*
     &              (-2*ctt2*(MA02 + MSb22 - MSt22)*MSt2diff + 
     &                (2*c2tt*MA02 + 2*c2tt*MSb22 + MSt2diff - 
     &                   c2tt*MSt2sum)*pMT2) + 
     &             ctt2*hb*ht*(MA02 + MSb22 - MSt22)*MSt2diff*MUEr*
     &              s2tb - 
     &             (MSt2diff + c2tt*(2*MA02 + 2*MSb22 - MSt2sum))*
     &              (hb*ht*pMT2*MUEr*s2tb - ht2*pMT*MUEr2*stb2)))*x31)/
     &     (MSt22*pMT) + 1/4.D0*
     &    (2*Abr*hb2*s2tt2*(MSt22*stb2*x14*x28 - 
     &          ctb2*MSt12*x11*x31) - 
     &       c2tt*s2tt*(hb*ht*MUEr*s2tb*
     &           (MSt22*x14*x28 + MSt12*x11*x31) + 
     &          2*hb2*pMT*(MSt22*stb2*x14*x28 - ctb2*MSt12*x11*x31))
     &       )/(MSt12*MSt22*pMT*x11*x14)

        tmp6 = tmp6 + x36*
     &    (c2tt*ctb2*hb2*MUEr2*s2tt + 
     &      (-(1/4.D0*(hb*ht*MUEr*s2tb*
     &              (-4*Atr*c2tt2*pMT + 
     &                (MSt2diff + 
     &                   c2tt*(-2*MSb22 + MSt2sum + 4*pMT2))*s2tt)))
     &           + 1/2.D0*(ht2*
     &            (pMT*(MSt2diff + 
     &                 c2tt*(-2*(Atr**2 + MSb22) + MSt2sum + 2*pMT2))
     &                *s2tt + 
     &              Atr*(-4*c2tt2*pMT2 + (-2*MSb22 + MSt2sum)*s2tt2))
     &             *stb2))/pMT) + 
     &   s2tt*(((hb2 + ht2)*MSt2diff + 
     &         c2tt*(hb2 - ht2)*(MSt2sum - 2*MUEr2))*(-x6 + x7) + 
     &      1/8.D0*(ht2*MUEr2*
     &         ((-2*(MA02*(3*(1 + c2tt2)*(MA02 - 2*MSt12)*MSt12 + 
     &                   (-((3 + 7*c2tt2)*MA02) + 8*c2tt2*MSt12)*
     &                    MSt22 + 2*(3 + 11*c2tt2)*MSt22**2) + 
     &                (3*(1 + c2tt2)*MSt12 - (3 + 7*c2tt2)*MSt22)*
     &                 MSt2diff**2)*x39)/MSt22**3 + 
     &           MA02*(-(((4*c2tt2*(MA02 - 6*MSt12)*MSt12 + 
     &                     3*(MA02 - 4*MSt12)*MSt2diff*s2tt2)*x38)/
     &                 MSt12**3) + 
     &              ((-4*c2tt2*(MA02 - 6*MSt22)*MSt22 + 
     &                   3*(MA02 - 4*MSt22)*MSt2diff*s2tt2)*x40)/
     &               MSt22**3))))

        tmp7 = 1/4.D0*(s2tt*
     &     ((hb2*(4*ctt2*MA02*MSt22 + 
     &            Abr**2*(MSt2diff + c2tt*MSt2sum) + MSt2diff*pMT2 - 
     &            c2tt*MSt2sum*pMT2 + 2*Abr*MSt2sum*pMT*s2tt - 
     &            4*MA02*MSt12*stt2))/(MSt12*MSt22) + 
     &       ((MA02 + MSb12 - MSt12)*(MA02 - MSb12 + MSt12)*
     &          (hb*ht*MUEr*s2tb*(Abr*s2tt + 2*pMT*stt2) + 
     &            2*(ctb2*ht2*MUEr2*stt2 + 
     &               hb2*stb2*(Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2))
     &            )*x11)/MSt12 - 
     &       ((MA02 + MSb12 - MSt22)*(MA02 - MSb12 + MSt22)*
     &          (2*ctb2*ctt2*ht2*MUEr2 + 
     &            2*ctt2*(hb*ht*pMT*MUEr*s2tb + hb2*pMT2*stb2) + 
     &            Abr*(-(hb*ht*MUEr*s2tb*s2tt) - 
     &               2*hb2*pMT*s2tt*stb2 + 2*Abr*hb2*stb2*stt2))*x12)
     &         /MSt22 + ((MA02**2 - (MSb22 - MSt12)**2)*
     &          (2*Abr**2*ctb2*ctt2*hb2 + 
     &            Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt + 
     &            2*(ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + 
     &               ht2*MUEr2*stb2)*stt2)*x13)/MSt12 + 
     &       ((-MA02**2 + (MSb22 - MSt22)**2)*
     &          (-2*ctt2*hb*ht*pMT*MUEr*s2tb + 
     &            Abr*hb*ht*MUEr*s2tb*s2tt + 2*ctt2*ht2*MUEr2*stb2 + 
     &            2*ctb2*hb2*(ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))
     &            )*x14)/MSt22 + 
     &       ht2*MUEr2*((s2tt2*
     &             (MA02 - MA02**2*(MA02 - 6*MSt12)*x16))/MSt12**2+
     &            ((1 + c2tt2)*MSt12 + 
     &             MA02**2*(MA02 - 6*MSt22)*s2tt2*x17)/MSt22**2) + 
     &       ht2*MUEr2*((2 - c2tt + c2tt2)/MSt12 - 
     &          (3 + c2tt)/MSt22 - 
     &          (2*(MA02 + c2tt2*MSt22))/MSt22**2 + 
     &          4*x15*(MSt22 + c2tt2*MSt22 + 
     &             1/4.D0*((1 + c2tt2)*
     &                (-6*MSt12 + (MA02 - MSt12)**3/MSt22**2 + 
     &                  (4*MSt12**2)/MSt22 + 
     &                  ((MA02 - MSt22)*(MA02 + MSt22))/MSt12 - 
     &                  (MA02*(2*MA02 + 2*MSt12 + MSt22))/MSt22))))
     &       ))

        tmp7 = tmp7 + 1/4.D0*
     &    ((2*Abr**2*c2tt*hb2*pMT*s2tt + 
     &         Abr*hb2*(-4*c2tt2*pMT2 + 
     &            s2tt2*(2*MA02 - 2*ctb2*MSb22 + MSt2sum - 
     &               2*MSb12*stb2)) + 
     &         s2tt*(hb2*MSt2diff*pMT + 
     &            c2tt*(-2*ht2*pMT*MUEr2 - hb*ht*MSb2diff*MUEr*s2tb + 
     &               hb2*pMT*
     &                (2*MA02 + 2*ctb2*MSb22 - MSt2sum - 2*pMT2 + 
     &                  2*MSb12*stb2))))*(x18 - x19))/(MSt2diff*pMT)

        s2ttF2t = tmp6/MSt2diff + tmp4*(-x10 + x18) + 
     &   tmp5*(-x10 + x19) + 
     &   1/8.D0*(((2*s2tt*(-(hb2*stb2*
     &                (-2*Abr*(MSb12 - MSt12)*(MSb12 - MSt22)*
     &                   MSt2sum*pMT*s2tt - 
     &                  4*MSb12**3*(ctt2*MSt22 - MSt12*stt2) + 
     &                  2*Abr**2*(MSb12 - MSt12)*(MSb12 - MSt22)*
     &                   (-(ctt2*MSt22) + MSt12*stt2) + 
     &                  2*MSt12*MSt22*pMT2*
     &                   (ctt2*MSt12 - MSt22*stt2) + 
     &                  MSb12**2*
     &                   (MSt2diff*(2*MSt2sum - pMT2) + 
     &                     c2tt*MSt2sum*(2*MSt2sum + pMT2) - 
     &                     4*ctt2*MSt22*MUEr2 + 4*MSt12*MUEr2*stt2)-
     &                    MSb12*
     &                   (2*MSt12*MSt22*
     &                      (MSt2diff + c2tt*MSt2sum) + 
     &                     MSt2sum*(-MSt2diff + c2tt*MSt2sum)*
     &                      pMT2 - 4*ctt2*MSt22**2*MUEr2 + 
     &                     4*MSt12**2*MUEr2*stt2))) + 
     &             ctb2*ht2*
     &              (4*Atr*MSb12*
     &                 (-MSt12**2 - MSt22**2 + MSb12*MSt2sum)*pMT*
     &                 s2tt - 
     &                4*MSb12**3*(ctt2*MSt12 - MSt22*stt2) - 
     &                2*MSt12*MSt22*MUEr2*
     &                 (ctt2*MSt12 - MSt22*stt2) - 
     &                4*Atr**2*MSb12*
     &                 (ctt2*(MSb12 - MSt12)*MSt12 + 
     &                   MSt22*(-MSb12 + MSt22)*stt2) + 
     &                2*MSb12**2*
     &                 (MSt2diff*(-MSt2sum + pMT2) + 
     &                   c2tt*MSt2sum*(MSt2sum + pMT2) - 
     &                   ctt2*MSt12*MUEr2 + MSt22*MUEr2*stt2) - 
     &                MSb12*
     &                 (4*ctt2*MSt22*(MSt12**2 + MSt22*pMT2) - 
     &                   2*ctt2*MSt12*MSt2sum*MUEr2 + 
     &                   2*
     &                    (-2*MSt12*(MSt22**2 + MSt12*pMT2) + 
     &                      MSt22*MSt2sum*MUEr2)*stt2)) - 
     &             hb*ht*MUEr*s2tb*
     &              (-(Abr*(MSb12 - MSt12)*(MSb12 - MSt22)*MSt2sum*
     &                   s2tt) - 
     &                2*Atr*MSb12*
     &                 (-MSt12**2 - MSt22**2 + MSb12*MSt2sum)*s2tt+
     &                  pMT*
     &                 (-2*MSb12**2*
     &                    (ctt2*MSt22 + MSt2diff - MSt12*stt2) + 
     &                   2*MSb12*MSt2diff*
     &                    (ctt2*MSt22 + MSt2sum + MSt12*stt2) + 
     &                   2*MSt12*MSt22*(ctt2*MSt12 - MSt22*stt2))))
     &           )/((MSb12 - MSt12)*MSt12*(MSb12 - MSt22)*MSt22) + 
     &        ((4*Abr**2*c2tt*hb2*pMT*s2tt*stb2 + 
     &              Abr*(2*hb2*(-2*MA02 + 2*MSb12 + MSt2sum)*s2tt2*
     &                  stb2 - 
     &                 4*c2tt2*(hb*ht*pMT*MUEr*s2tb + 2*hb2*pMT2*stb2)
     &                 ) + 
     &              s2tt*(MSt2diff*
     &                  (4*ctb2*ht2*pMT + 3*hb*ht*MUEr*s2tb + 
     &                    2*hb2*pMT*stb2) + 
     &                 c2tt*
     &                  (4*ctb2*ht2*pMT*
     &                     (-2*MSb12 + MSt2sum + 2*pMT2 - MUEr2) + 
     &                    hb*ht*
     &                     (2*MA02 - 2*MSb12 + MSt2sum + 4*pMT2)*
     &                     MUEr*s2tb + 
     &                    2*hb2*pMT*
     &                     (2*MA02 + 2*MSb12 - MSt2sum - 2*pMT2 + 
     &                       4*MUEr2)*stb2)))*(x18 - x19) + 
     &           8*Atr**2*c2tt*ctb2*ht2*pMT*s2tt*(-x18 + x19) + 
     &           4*Atr*(2*c2tt2*hb*ht*pMT*MUEr*s2tb*(-x18 + x19) + 
     &              ctb2*ht2*
     &               (-(MSt2sum*s2tt2*x10) - 4*c2tt2*pMT2*x18 + 
     &                 MSt2sum*s2tt2*x18 + 
     &                 2*MSt22*s2tt2*(x10 - x19) + 
     &                 4*c2tt2*pMT2*x19 + 
     &                 MSt2diff*s2tt2*(-x10 + x19))))/(MSt2diff*pMT)
     &        )*(-x10 + x21))

        s2ttF2t = s2ttF2t + 
     &   1/8.D0*(((2*s2tt*(-(ctb2*hb2*
     &                (-2*Abr*(MSb22 - MSt12)*(MSb22 - MSt22)*
     &                   MSt2sum*pMT*s2tt - 
     &                  4*MSb22**3*(ctt2*MSt22 - MSt12*stt2) + 
     &                  2*Abr**2*(MSb22 - MSt12)*(MSb22 - MSt22)*
     &                   (-(ctt2*MSt22) + MSt12*stt2) + 
     &                  2*MSt12*MSt22*pMT2*
     &                   (ctt2*MSt12 - MSt22*stt2) + 
     &                  MSb22**2*
     &                   (MSt2diff*(2*MSt2sum - pMT2) + 
     &                     c2tt*MSt2sum*(2*MSt2sum + pMT2) - 
     &                     4*ctt2*MSt22*MUEr2 + 4*MSt12*MUEr2*stt2)-
     &                    MSb22*
     &                   (2*MSt12*MSt22*
     &                      (MSt2diff + c2tt*MSt2sum) + 
     &                     MSt2sum*(-MSt2diff + c2tt*MSt2sum)*
     &                      pMT2 - 4*ctt2*MSt22**2*MUEr2 + 
     &                     4*MSt12**2*MUEr2*stt2))) + 
     &             ht2*stb2*
     &              (4*Atr*MSb22*
     &                 (-MSt12**2 - MSt22**2 + MSb22*MSt2sum)*pMT*
     &                 s2tt - 
     &                4*MSb22**3*(ctt2*MSt12 - MSt22*stt2) - 
     &                2*MSt12*MSt22*MUEr2*
     &                 (ctt2*MSt12 - MSt22*stt2) - 
     &                4*Atr**2*MSb22*
     &                 (ctt2*(MSb22 - MSt12)*MSt12 + 
     &                   MSt22*(-MSb22 + MSt22)*stt2) + 
     &                2*MSb22**2*
     &                 (MSt2diff*(-MSt2sum + pMT2) + 
     &                   c2tt*MSt2sum*(MSt2sum + pMT2) - 
     &                   ctt2*MSt12*MUEr2 + MSt22*MUEr2*stt2) - 
     &                MSb22*
     &                 (4*ctt2*MSt22*(MSt12**2 + MSt22*pMT2) - 
     &                   2*ctt2*MSt12*MSt2sum*MUEr2 + 
     &                   2*
     &                    (-2*MSt12*(MSt22**2 + MSt12*pMT2) + 
     &                      MSt22*MSt2sum*MUEr2)*stt2)) - 
     &             hb*ht*MUEr*s2tb*
     &              (Abr*(MSb22 - MSt12)*(MSb22 - MSt22)*MSt2sum*
     &                 s2tt + 
     &                2*Atr*MSb22*
     &                 (-MSt12**2 - MSt22**2 + MSb22*MSt2sum)*s2tt+
     &                  pMT*
     &                 (2*MSb22**2*
     &                    (ctt2*MSt22 + MSt2diff - MSt12*stt2) - 
     &                   2*MSb22*MSt2diff*
     &                    (ctt2*MSt22 + MSt2sum + MSt12*stt2) + 
     &                   2*MSt12*MSt22*(-(ctt2*MSt12) + MSt22*stt2)
     &                   ))))/
     &         ((MSb22 - MSt12)*MSt12*(MSb22 - MSt22)*MSt22) + 
     &        ((4*Abr**2*c2tt*ctb2*hb2*pMT*s2tt + 
     &              Abr*(c2tt2*
     &                  (-8*ctb2*hb2*pMT2 + 4*hb*ht*pMT*MUEr*s2tb) + 
     &                 2*ctb2*hb2*(-2*MA02 + 2*MSb22 + MSt2sum)*
     &                  s2tt2) + 
     &              s2tt*(MSt2diff*
     &                  (2*ctb2*hb2*pMT - 3*hb*ht*MUEr*s2tb + 
     &                    4*ht2*pMT*stb2) + 
     &                 c2tt*
     &                  (2*ctb2*hb2*pMT*
     &                     (2*MA02 + 2*MSb22 - MSt2sum - 2*pMT2 + 
     &                       4*MUEr2) - 
     &                    hb*ht*
     &                     (2*MA02 - 2*MSb22 + MSt2sum + 4*pMT2)*
     &                     MUEr*s2tb + 
     &                    4*ht2*pMT*
     &                     (-2*MSb22 + MSt2sum + 2*pMT2 - MUEr2)*stb2
     &                    )))*(x18 - x19) + 
     &           8*Atr**2*c2tt*ht2*pMT*s2tt*stb2*(-x18 + x19) + 
     &           4*Atr*(2*c2tt2*hb*ht*pMT*MUEr*s2tb*(x18 - x19) + 
     &              ht2*stb2*
     &               (-(MSt2sum*s2tt2*x10) - 4*c2tt2*pMT2*x18 + 
     &                 MSt2sum*s2tt2*x18 + 
     &                 2*MSt22*s2tt2*(x10 - x19) + 
     &                 4*c2tt2*pMT2*x19 + 
     &                 MSt2diff*s2tt2*(-x10 + x19))))/(MSt2diff*pMT)
     &        )*(-x10 + x22))

        s2ttF2t = s2ttF2t + 
     &   1/4.D0*((2*(-(ht2*(4*Atr*c2tt2*pMT2 + 
     &                 2*pMT*
     &                  (-MSt2diff + 
     &                    c2tt2*(Atr**2 + MUEr2 + 2*MSt12*(1 + Nc)))*
     &                  s2tt - 
     &                 (Atr*MSt2sum + MSt2diff*pMT*(1 + Nc)*s2tt)*
     &                  s2tt2)*(x10 - x18)**2) - 
     &            ht2*(-4*Atr*c2tt2*pMT2 + 
     &               2*pMT*(MSt2diff + 
     &                  c2tt2*(Atr**2 + MUEr2 + 2*MSt22*(1 + Nc)))*
     &                s2tt + 
     &               (Atr*(2*MSt22 - MSt2diff) + 
     &                  MSt2diff*pMT*(1 + Nc)*s2tt)*s2tt2)*
     &             (x10 - x19)**2))/MSt2diff + 
     &       ((-(hb2*(4*Abr*(2*MSt12*MSt22 + MSt2sum*pMT2)*s2tt2 + 
     &                  2*pMT*s2tt*
     &                   (MSt2diff*pMT2 + 
     &                     c2tt*(2*MSt12*MSt22 - MSt2sum*pMT2) + 
     &                     2*(Abr**2 + MA02)*
     &                      (ctt2*MSt22 - MSt12*stt2) + 
     &                     2*(ctb2*MSb22 - MUEr2 + MSb12*stb2)*
     &                      (ctt2*MSt22 - MSt12*stt2)))) + 
     &             ht2*(-12*Atr*(2*MSt12*MSt22 + MSt2sum*pMT2)*
     &                 s2tt2 + 
     &                2*pMT*s2tt*
     &                 (-((MSt2diff + c2tt*MSt2sum)*pMT2) - 
     &                   2*MSt12*
     &                    (-(c2tt*MSt22) + MUEr2 + ctt2*MUEr2) - 
     &                   MSt12**2*(-2 + s2tt2 + Nc*s2tt2) + 
     &                   MSt22**2*(-2 + s2tt2 + Nc*s2tt2) + 
     &                   2*MSt22*MUEr2*(1 + stt2) + 
     &                   2*(ctb2*MSb12 + MSb22*stb2)*
     &                    (ctt2*MSt12 - MSt22*stt2) + 
     &                   Atr**2*
     &                    (2*(1 + ctt2)*MSt12 - 2*MSt22*(1 + stt2))
     &                   )))/MSt12 + 
     &          ((-2*Abr*hb2*s2tt2*stb2 + 
     &               c2tt*s2tt*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2))*
     &             x29)/(MSt2diff*x12))/MSt22)/pMT

        s2ttF2t = s2ttF2t + 
     &   1/4.D0*(s2tt*(((MA02 + MSb12 - MSt12)*
     &           (2*Abr**2*ctt2*hb2*stb2 + 
     &             Abr*s2tt*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2) + 
     &             2*(ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb + 
     &                hb2*pMT2*stb2)*stt2)*x11*
     &           (MA02*(x10 - x21 - x28) - 
     &             MSb12*(x10 - x21 + x28) + 
     &             MSt12*(x10 - 2*x18 + x21 + x28)))/MSt12 - 
     &        ((MA02 + MSb12 - MSt22)*
     &           (2*ctb2*ctt2*ht2*MUEr2 + 
     &             2*ctt2*(hb*ht*pMT*MUEr*s2tb + hb2*pMT2*stb2) + 
     &             Abr*(-(hb*ht*MUEr*s2tb*s2tt) - 
     &                2*hb2*pMT*s2tt*stb2 + 2*Abr*hb2*stb2*stt2))*
     &           x12*(MA02*(x10 - x21 - x29) - 
     &             MSb12*(x10 - x21 + x29) + 
     &             MSt22*(x10 - 2*x19 + x21 + x29)))/MSt22 + 
     &        ((MA02 + MSb22 - MSt12)*
     &           (2*Abr**2*ctb2*ctt2*hb2 + 
     &             Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt + 
     &             2*(ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + 
     &                ht2*MUEr2*stb2)*stt2)*x13*
     &           (MA02*(x10 - x22 - x30) - 
     &             MSb22*(x10 - x22 + x30) + 
     &             MSt12*(x10 - 2*x18 + x22 + x30)))/MSt12 - 
     &        ((MA02 + MSb22 - MSt22)*
     &           (-2*ctt2*hb*ht*pMT*MUEr*s2tb + 
     &             Abr*hb*ht*MUEr*s2tb*s2tt + 2*ctt2*ht2*MUEr2*stb2 + 
     &             2*ctb2*hb2*
     &              (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*x14*
     &           (MA02*(x10 - x22 - x31) - 
     &             MSb22*(x10 - x22 + x31) + 
     &             MSt22*(x10 - 2*x19 + x22 + x31)))/MSt22))

        s2ttF2t = s2ttF2t + 
     &   (s2tt*(-(ht2*MUEr2*x41) + 
     &         1/2.D0*(ht2*MSt2diff*(3*pMT2 + MUEr2)*(x10 - x20) + 
     &             ((ht2*(-(pMT2*
     &                       (MSt2diff*pMT2 + 
     &                       (5 - 2*c2tt)*MSt22*MUEr2)) + 
     &                     MSt12*
     &                      (MSt22*MSt2diff + 
     &                       (5 + 2*c2tt)*pMT2*MUEr2)) + 
     &                  4*hb2*pMT2*MUEr2*(-(ctt2*MSt22) + MSt12*stt2)
     &                  )*(-x10 + x24))/pMT2)/MSt22) + 
     &      1/8.D0*((2*(Abr*(MSt2diff*s2tt2*
     &                   (hb*ht*pMT*MUEr*s2tb - 
     &                     hb2*(MA02 + MSb12 - MSt12 - 2*pMT2)*stb2)
     &                    + 2*c2tt2*(MA02 + MSb12 - MSt12)*
     &                   (hb*ht*pMT*MUEr*s2tb + 2*hb2*pMT2*stb2)) - 
     &               2*Abr**2*hb2*pMT*s2tt*stb2*
     &                (c2tt*(MA02 + MSb12) - ctt2*MSt22 + 
     &                  MSt12*stt2) + 
     &               s2tt*(c2tt*(2*MA02 + 2*MSb12 - MSt2sum)*
     &                   (ctb2*ht2*pMT*MUEr2 + hb*ht*pMT2*MUEr*s2tb + 
     &                     hb2*pMT*pMT2*stb2) + 
     &                  MSt2diff*
     &                   (ctb2*ht2*pMT*MUEr2 + 
     &                     hb2*pMT*stb2*
     &                      (pMT2 - 2*(MA02 + MSb12 - MSt12)*stt2)+
     &                       hb*ht*MUEr*s2tb*
     &                      (pMT2 - (MA02 + MSb12 - MSt12)*stt2))))*
     &             x28)/pMT + 
     &          (2*ht2*(MSt2diff + 
     &                c2tt*(2*MA02 + 2*MSb22 - MSt2sum))*MUEr2*s2tt*
     &              stb2 + 
     &             (hb*ht*MUEr*s2tb*
     &                 (-2*Abr*pMT*
     &                    (2*c2tt2*(MA02 + MSb22 - MSt12) + 
     &                      MSt2diff*s2tt2) + 
     &                   2*s2tt*
     &                    ((-MSt2diff + 
     &                       c2tt*(-2*(MA02 + MSb22) + MSt2sum))*
     &                       pMT2 + 
     &                      (MA02 + MSb22 - MSt12)*MSt2diff*stt2))+
     &                  2*
     &                 (ctb2*hb2*
     &                    (Abr*
     &                       (4*c2tt2*(MA02 + MSb22 - MSt12)*pMT2 - 
     &                       MSt2diff*
     &                       (MA02 + MSb22 - MSt12 - 2*pMT2)*s2tt2)+
     &                        pMT*s2tt*
     &                       ((MSt2diff + 
     &                       c2tt*(2*(MA02 + MSb22) - MSt2sum))*pMT2
     &                        - 2*
     &                       ((MA02 + MSb22 - MSt12)*MSt2diff*
     &                      stt2 + 
     &                       Abr**2*
     &                       (c2tt*(MA02 + MSb22) - ctt2*MSt22 + 
     &                       MSt12*stt2)))) + 
     &                   (-2*c2tt*ctb2*hb2*pMT*s2tt + 
     &                      c2tt*hb*ht*MUEr*s2tb*s2tt + 
     &                      2*Abr*ctb2*hb2*s2tt2)/x13))/pMT)*x30)/
     &        MSt2diff)/MSt12

        s2ttF2t = s2ttF2t + tmp7*(-x10 + x9) + 
     &   ht2*(1/4.D0*((1 + c2tt2)*MUEr2*s2tt*x15*
     &         ((((MSt12 - MSt22)*
     &                 (MA02**3 - MSt12**3 + 2*MSt12**2*MSt22 + 
     &                   MSt12*MSt22**2 - 2*MSt22**3 - 
     &                   MA02**2*(3*MSt12 + 4*MSt22) + 
     &                   MA02*
     &                    (3*MSt12**2 + 2*MSt12*MSt22 + 5*MSt22**2)
     &                   ))/(MSt22**2*MSt2diff) - 
     &              2*(MA02 + MSt2diff))*(-x10 + x18) - 
     &           ((MA02**2*
     &                 (-MSt22**2 + MSt12*(MSt22 - 2*MSt2diff)) + 
     &                4*MA02*MSt12*(MSt12 + 2*MSt22)*MSt2diff + 
     &                (MSt12 - MSt22)*
     &                 (-MSt22**3 + 
     &                   MSt12**2*(MSt22 - 2*MSt2diff) + 
     &                   4*MSt12*MSt22*MSt2diff))*(x10 - x19))/
     &            (MSt12*MSt22*MSt2diff) - 
     &           ((MA02**4 + 
     &                MA02**2*
     &                 (6*MSt12**2 + 4*MSt12*MSt22 + 6*MSt22**2) + 
     &                4*MA02*MSt12**2*MSt2diff + MSt2diff**4 - 
     &                4*MA02**3*MSt2sum)*x39)/MSt22**3)) + 
     &      s2tt*(((-MSt12 + pMT2 + MUEr2)**2*x26*
     &            (MUEr2*x41 + 
     &              1/2.D0*((MSt12 - pMT2 + MUEr2)*(-x10 + x20)) - 
     &              1/2.D0*((-(MSt12 - pMT2)**2 + (MSt12 + pMT2)*MUEr2)*
     &                  (-x10 + x24))/pMT2))/MSt12 + 
     &         MUEr2*(x42/MSt22 + 
     &            (1/4.D0*((2*c2tt2*MSt12 + MSt2diff*s2tt2)*x38)/
     &                 (MSt2diff*x16) + 
     &               1/8.D0*(MA02**2*(MA02 - 6*MSt12)*s2tt2*x16*
     &                  (MA02*x38 - 2*MSt12*(x10 - x18 + x38))))/
     &             MSt12**3 + 
     &            (1/4.D0*((-2*
     &                      (MSt2diff + 
     &                      c2tt2*(2*MSt22 + MSt2diff))*x39)/x15 + 
     &                   ((2*c2tt2*MSt22 - MSt2diff*s2tt2)*x40)/x17
     &                   )/MSt2diff - 
     &               1/8.D0*(MA02**2*(MA02 - 6*MSt22)*s2tt2*x17*
     &                  (MA02*x40 - 2*MSt22*(x10 - x19 + x40))))/
     &             MSt22**3) + 
     &         1/2.D0*((-MSt22 + pMT2 + MUEr2)**2*x27*
     &             ((MSt22 - pMT2 + MUEr2)*(x10 - x20) + 
     &               ((-(MSt22 - pMT2)**2 + (MSt22 + pMT2)*MUEr2)*
     &                  (-x10 + x24))/pMT2 - 2*MUEr2*x42))/MSt22))

        tmp8 = 1/4.D0*(s2tt2*
     &      (hb*ht*MUEr*s2tb*
     &         (MSt2diff*(Abr*(MSb12 - MSt12)*(MSb12 - MSt22) + 
     &              2*Atr*MSb12*(MSb12 - MSt2sum))*s2tt + 
     &           pMT*(MSt12*MSt22*(-(c2tt*MSt2diff) + MSt2sum) + 
     &              2*MSb12**2*
     &               (ctt2*MSt22 + MSt2sum + MSt12*stt2) - 
     &              2*MSb12*
     &               (MSt12*MSt22 + (1 + ctt2)*MSt22**2 + 
     &                 MSt12**2*(1 + stt2)))) + 
     &        (-(hb2*stb2*(-(MSb12**2*
     &                   (MSt2diff*MSt2sum*pMT2 - 
     &                     c2tt*
     &                      (20*MSt12*MSt22*MSt2sum + 
     &                       MSt2diff**2*pMT2) + 
     &                     4*ctt2*MSt22*
     &                      (-(MSt22*MSt2sum) + MSt2diff*MUEr2) + 
     &                     4*MSt12*(MSt12*MSt2sum + MSt2diff*MUEr2)*
     &                      stt2)) + 
     &                MSt2diff*
     &                 (MSt12*MSt22*(c2tt*MSt2diff - MSt2sum)*
     &                    pMT2 - 
     &                   2*Abr**2*(MSb12 - MSt12)*(MSb12 - MSt22)*
     &                    (ctt2*MSt22 + MSt12*stt2)) - 
     &                2*(Abr*(MSb12 - MSt12)*(MSb12 - MSt22)*
     &                    MSt2diff**2*pMT*s2tt + 
     &                   2*MSb12**3*
     &                    (5*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &                      MSt12**2*stt2)) - 
     &                MSb12*
     &                 (MSt2diff*(c2tt*MSt2diff - MSt2sum)*MSt2sum*
     &                    pMT2 + 
     &                   4*MSt12*MSt22*
     &                    (5*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &                      MSt12**2*stt2) - 
     &                   4*MSt2diff*MUEr2*
     &                    (ctt2*MSt22**2 + MSt12**2*stt2)))) + 
     &           ctb2*ht2*(MSb12**2*
     &               (20*c2tt*MSt12*MSt22*MSt2sum + 
     &                 4*ctt2*
     &                  (MSt12**2*MSt2sum + MSt22*MSt2diff*pMT2) - 
     &                 c2tt*MSt2diff**2*MUEr2 + 
     &                 MSt2diff*MSt2sum*MUEr2 - 
     &                 4*MSt22**2*MSt2sum*stt2 + 
     &                 4*MSt12*MSt2diff*pMT2*stt2) - 
     &              4*MSb12**3*
     &               (ctt2*MSt12**2 + 
     &                 MSt22*(5*c2tt*MSt12 - MSt22*stt2)) + 
     &              MSt2diff*
     &               (MSt12*MSt22*(-(c2tt*MSt2diff) + MSt2sum)*
     &                  MUEr2 - 
     &                 2*Atr**2*MSb12*
     &                  (2*ctt2*MSt12**2 + c2tt*MSb12*MSt2diff - 
     &                    MSb12*MSt2sum + 2*MSt22**2*stt2)) + 
     &              MSb12*(-4*ctt2*MSt12**3*MSt22 - 
     &                 MSt22**2*
     &                  (4*ctt2*MSt2diff*pMT2 + MSt2sum*MUEr2) + 
     &                 MSt2diff**2*
     &                  (c2tt*MSt2sum*MUEr2 + 
     &                    4*Atr*(MSb12 - MSt2sum)*pMT*s2tt) + 
     &                 4*MSt12*MSt22**3*stt2 - 
     &                 MSt12**2*
     &                  (20*c2tt*MSt22**2 - MSt2sum*MUEr2 + 
     &                    4*MSt2diff*pMT2*stt2))))/MSt2diff))/
     &    ((MSb12 - MSt12)*MSt12*(MSb12 - MSt22)*MSt22)

        tmp8 = tmp8 + 1/2.D0*
     &    (s2tt*(2*Atr**2*ctb2*ht2*(MSt2diff - 3*c2tt*MSt2sum)*
     &          s2tt + 2*Atr*MSt2sum*
     &          (2*ctb2*ht2*pMT + hb*ht*MUEr*s2tb)*(-2 + 3*s2tt2) + 
     &         Abr**2*hb2*(MSt2diff + 
     &            c2tt*(-6*MA02 + 6*MSb12 + 3*MSt2sum))*s2tt*stb2-
     &           Abr*(2*MA02 - 2*MSb12 - MSt2sum)*(-2 + 3*s2tt2)*
     &          (hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2) + 
     &         s2tt*(ctb2*ht2*
     &             (2*MSt2diff*pMT2 + 6*c2tt*MSt2sum*pMT2 + 
     &               6*c2tt*MA02*MUEr2 + MSt2diff*MUEr2 - 
     &               3*c2tt*MSt2sum*MUEr2 + 
     &               2*MSb12*(MSt2diff - 3*c2tt*(MSt2sum + MUEr2)))+
     &              MSt2diff*
     &             (3*hb*ht*pMT*MUEr*s2tb + 
     &               hb2*(2*MSb12 + pMT2 + 2*MUEr2)*stb2) + 
     &            3*c2tt*(hb*ht*(2*MA02 - 2*MSb12 + MSt2sum)*pMT*
     &                MUEr*s2tb + 
     &               hb2*(2*MSb12*MSt2sum + 2*MA02*pMT2 - 
     &                  2*MSb12*pMT2 - MSt2sum*pMT2 + 2*MSt2sum*MUEr2)
     &                 *stb2)))*(x18 - x19))/MSt2diff**2

        tmp9 = 1/4.D0*(s2tt2*
     &      (-(hb*ht*MUEr*s2tb*
     &           (MSt2diff*
     &              (Abr*(MSb22 - MSt12)*(MSb22 - MSt22) + 
     &                2*Atr*MSb22*(MSb22 - MSt2sum))*s2tt + 
     &             pMT*(MSt12*MSt22*(-(c2tt*MSt2diff) + MSt2sum) + 
     &                2*MSb22**2*
     &                 (ctt2*MSt22 + MSt2sum + MSt12*stt2) - 
     &                2*MSb22*
     &                 (MSt12*MSt22 + (1 + ctt2)*MSt22**2 + 
     &                   MSt12**2*(1 + stt2))))) + 
     &        (-(ctb2*hb2*(-(MSb22**2*
     &                   (MSt2diff*MSt2sum*pMT2 - 
     &                     c2tt*
     &                      (20*MSt12*MSt22*MSt2sum + 
     &                       MSt2diff**2*pMT2) + 
     &                     4*ctt2*MSt22*
     &                      (-(MSt22*MSt2sum) + MSt2diff*MUEr2) + 
     &                     4*MSt12*(MSt12*MSt2sum + MSt2diff*MUEr2)*
     &                      stt2)) + 
     &                MSt2diff*
     &                 (MSt12*MSt22*(c2tt*MSt2diff - MSt2sum)*
     &                    pMT2 - 
     &                   2*Abr**2*(MSb22 - MSt12)*(MSb22 - MSt22)*
     &                    (ctt2*MSt22 + MSt12*stt2)) - 
     &                2*(Abr*(MSb22 - MSt12)*(MSb22 - MSt22)*
     &                    MSt2diff**2*pMT*s2tt + 
     &                   2*MSb22**3*
     &                    (5*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &                      MSt12**2*stt2)) - 
     &                MSb22*
     &                 (MSt2diff*(c2tt*MSt2diff - MSt2sum)*MSt2sum*
     &                    pMT2 + 
     &                   4*MSt12*MSt22*
     &                    (5*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &                      MSt12**2*stt2) - 
     &                   4*MSt2diff*MUEr2*
     &                    (ctt2*MSt22**2 + MSt12**2*stt2)))) + 
     &           ht2*stb2*(MSb22**2*
     &               (20*c2tt*MSt12*MSt22*MSt2sum + 
     &                 4*ctt2*
     &                  (MSt12**2*MSt2sum + MSt22*MSt2diff*pMT2) - 
     &                 c2tt*MSt2diff**2*MUEr2 + 
     &                 MSt2diff*MSt2sum*MUEr2 - 
     &                 4*MSt22**2*MSt2sum*stt2 + 
     &                 4*MSt12*MSt2diff*pMT2*stt2) - 
     &              4*MSb22**3*
     &               (ctt2*MSt12**2 + 
     &                 MSt22*(5*c2tt*MSt12 - MSt22*stt2)) + 
     &              MSt2diff*
     &               (MSt12*MSt22*(-(c2tt*MSt2diff) + MSt2sum)*
     &                  MUEr2 - 
     &                 2*Atr**2*MSb22*
     &                  (2*ctt2*MSt12**2 + c2tt*MSb22*MSt2diff - 
     &                    MSb22*MSt2sum + 2*MSt22**2*stt2)) + 
     &              MSb22*(-4*ctt2*MSt12**3*MSt22 - 
     &                 MSt22**2*
     &                  (4*ctt2*MSt2diff*pMT2 + MSt2sum*MUEr2) + 
     &                 MSt2diff**2*
     &                  (c2tt*MSt2sum*MUEr2 + 
     &                    4*Atr*(MSb22 - MSt2sum)*pMT*s2tt) + 
     &                 4*MSt12*MSt22**3*stt2 - 
     &                 MSt12**2*
     &                  (20*c2tt*MSt22**2 - MSt2sum*MUEr2 + 
     &                    4*MSt2diff*pMT2*stt2))))/MSt2diff))/
     &    ((MSb22 - MSt12)*MSt12*(MSb22 - MSt22)*MSt22)

        tmp9 = tmp9 + 1/2.D0*
     &    ((Abr**2*ctb2*hb2*
     &          (MSt2diff + c2tt*(-6*MA02 + 6*MSb22 + 3*MSt2sum))*
     &          s2tt2 - Abr*(2*MA02 - 2*MSb22 - MSt2sum)*
     &          (2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt*
     &          (-2 + 3*s2tt2) + 
     &         2*Atr**2*ht2*(MSt2diff - 3*c2tt*MSt2sum)*s2tt2*
     &          stb2 - 2*Atr*MSt2sum*s2tt*(-2 + 3*s2tt2)*
     &          (hb*ht*MUEr*s2tb - 2*ht2*pMT*stb2) + 
     &         s2tt2*(ctb2*hb2*
     &             (6*c2tt*MA02*pMT2 + MSt2diff*pMT2 - 
     &               3*c2tt*MSt2sum*pMT2 + 
     &               2*MSb22*
     &                (MSt2diff + 3*c2tt*MSt2sum - 3*c2tt*pMT2) + 
     &               2*MSt2diff*MUEr2 + 6*c2tt*MSt2sum*MUEr2) + 
     &            MSt2diff*
     &             (-3*hb*ht*pMT*MUEr*s2tb + 
     &               ht2*(2*MSb22 + 2*pMT2 + MUEr2)*stb2) - 
     &            3*c2tt*(hb*ht*(2*MA02 - 2*MSb22 + MSt2sum)*pMT*
     &                MUEr*s2tb + 
     &               ht2*(-2*MSt2sum*pMT2 - 2*MA02*MUEr2 + 
     &                  MSt2sum*MUEr2 + 2*MSb22*(MSt2sum + MUEr2))*
     &                stb2)))*(x18 - x19))/MSt2diff**2

        tmp12 = (-3*Abr**2*c2tt*hb2*s2tt2*
     &       (stb2*x14*x29 + ctb2*x12*x31) - 
     &      Abr*s2tt*(-2 + 3*s2tt2)*
     &       (hb*ht*MUEr*s2tb*(x14*x29 - x12*x31) + 
     &         2*hb2*pMT*(stb2*x14*x29 + ctb2*x12*x31)) + 
     &      3*c2tt*s2tt2*(hb*ht*pMT*MUEr*s2tb*(x14*x29 - x12*x31) + 
     &         ctb2*(ht2*MUEr2*x14*x29 + hb2*pMT2*x12*x31) + 
     &         stb2*(hb2*pMT2*x14*x29 + ht2*MUEr2*x12*x31)))/
     &    (MSt22*x12*x14) + 
     &   (-(Atr**2*ctb2*ht2*
     &         (6*c2tt*MSb12 + MSt2diff - 3*c2tt*MSt2sum)*s2tt2) + 
     &      Atr*(2*MSb12 - MSt2sum)*
     &       (2*ctb2*ht2*pMT + hb*ht*MUEr*s2tb)*s2tt*(-2 + 3*s2tt2)+
     &        (6*c2tt*MSb12 - MSt2diff - 3*c2tt*MSt2sum)*s2tt2*
     &       (ctb2*ht2*pMT2 + hb*ht*pMT*MUEr*s2tb + hb2*MUEr2*stb2))*
     &    x33 + (-(Atr**2*ht2*
     &         (6*c2tt*MSb22 + MSt2diff - 3*c2tt*MSt2sum)*s2tt2*
     &         stb2) - Atr*(2*MSb22 - MSt2sum)*s2tt*(-2 + 3*s2tt2)*
     &       (hb*ht*MUEr*s2tb - 2*ht2*pMT*stb2) + 
     &      (6*c2tt*MSb22 - MSt2diff - 3*c2tt*MSt2sum)*s2tt2*
     &       (ctb2*hb2*MUEr2 - hb*ht*pMT*MUEr*s2tb + ht2*pMT2*stb2))*
     &    x34 + (Atr**2*ctb2*ht2*
     &       (6*c2tt*MSb12 + MSt2diff - 3*c2tt*MSt2sum)*s2tt2 - 
     &      Atr*(2*MSb12 - MSt2sum)*
     &       (2*ctb2*ht2*pMT + hb*ht*MUEr*s2tb)*s2tt*(-2 + 3*s2tt2)-
     &        (6*c2tt*MSb12 - MSt2diff - 3*c2tt*MSt2sum)*s2tt2*
     &       (ctb2*ht2*pMT2 + hb*ht*pMT*MUEr*s2tb + hb2*MUEr2*stb2))*
     &    x35 + (Atr**2*ht2*
     &       (6*c2tt*MSb22 + MSt2diff - 3*c2tt*MSt2sum)*s2tt2*stb2+
     &        Atr*(2*MSb22 - MSt2sum)*s2tt*(-2 + 3*s2tt2)*
     &       (hb*ht*MUEr*s2tb - 2*ht2*pMT*stb2) - 
     &      (6*c2tt*MSb22 - MSt2diff - 3*c2tt*MSt2sum)*s2tt2*
     &       (ctb2*hb2*MUEr2 - hb*ht*pMT*MUEr*s2tb + ht2*pMT2*stb2))*
     &    x36

        tmp12 = tmp12 + 1/2.D0*
     &    (ht2*(s2tt*((4*Atr*MSt2sum*pMT*(-2 + 3*s2tt2) + 
     &              s2tt*(c2tt2*
     &                  (8*MA02*MUEr2 - 4*MSt22*(Atr**2 + MUEr2) + 
     &                    4*MSt12*
     &                     (-3*(Atr**2 + MUEr2) - 2*MSt2sum*(1 + Nc))
     &                    ) + 
     &                 MSt2diff*
     &                  (4*pMT2 + 
     &                    (Atr**2 + MSt2sum + MUEr2 + MSt2sum*Nc)*
     &                     s2tt2)))*(x10 - x18)**2 - 
     &           (-4*Atr*pMT*
     &               (2*c2tt2*MSt12 + MSt22*(2 - 4*s2tt2) + 
     &                 MSt2diff*s2tt2) + 
     &              s2tt*(4*c2tt2*
     &                  ((-2*MA02 + MSt12)*MUEr2 + 
     &                    MSt22*(3*MUEr2 + 2*MSt2sum*(1 + Nc))) + 
     &                 Atr**2*
     &                  (8*MSt2sum - 11*MSt12*s2tt2 - 
     &                    5*MSt22*s2tt2) + 
     &                 MSt2diff*
     &                  (4*pMT2 + 
     &                    (MSt2sum + MUEr2 + MSt2sum*Nc)*s2tt2)))*
     &            (x10 - x19)**2) + 
     &        (((1 + c2tt2)*MSt12**2 - (3 + 7*c2tt2)*MSt12*MSt22 + 
     &             (2 + 22*c2tt2)*MSt22**2)*MUEr2*s2tt2*x39)/
     &         (MSt22**3*x15))) - 
     &   1/4.D0*(s2tt2*(4*c2tt2*ht2*MSt12*(4*MSt12 - MSt2diff)*MUEr2*
     &          x38 + 4*hb2*MSt12**3*
     &          (MSt2diff*(MSt2sum - 2*MUEr2) + 
     &            c2tt*(MSt12**2 + 4*MSt12*MSt22 + MSt22**2 + 
     &               6*MUEr2*(-MSt2sum + MUEr2)))*x16*(x6 - x7) + 
     &         ht2*(MSt12*MSt2diff*MUEr2*s2tt2*x38 - 
     &            MSt2diff**2*MUEr2*s2tt2*x38 - 
     &            16*c2tt*MSt12**4*MSt22*x16*(x6 - x7) - 
     &            8*MSt12**3*
     &             (MSt2diff*MUEr2 + 
     &               3*c2tt*MUEr2*(-MSt2sum + MUEr2) - MSt22**2*stt2)
     &              *x16*(x6 - x7) + 8*ctt2*MSt12**5*x16*(-x6 + x7)
     &            )))/(MSt12**3*x16)

        tmp13 = 1/4.D0*((MA02**2 - (MSb22 - MSt12)**2)*s2tt2*
     &      (2*Abr**2*ctb2*ctt2*hb2 + 
     &        Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt + 
     &        2*(ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + ht2*MUEr2*stb2)*
     &         stt2)*x13)/MSt12

        tmp13 = tmp13 + 1/4.D0*
     &    ((MA02**2 - (MSb22 - MSt22)**2)*s2tt2*
     &       (-2*ctt2*hb*ht*pMT*MUEr*s2tb + Abr*hb*ht*MUEr*s2tb*s2tt + 
     &         2*ctt2*ht2*MUEr2*stb2 + 
     &         2*ctb2*hb2*(ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*
     &       x14)/MSt22

        tmp13 = tmp13 + 1/4.D0*
     &    (s2tt2*((24*c2tt*hb2*MA02)/MSt2diff + 
     &         (2*ht2*MA02*MUEr2)/MSt22**2 + 
     &         2*((2*ctt2*hb2*MA02 + ht2*MUEr2 + Abr*hb2*pMT*s2tt)/
     &             MSt12 + 
     &            (hb2*(-(Abr*pMT*s2tt) + 2*MA02*stt2))/MSt22) + 
     &         ((-c2tt + c2tt2)*ht2*MUEr2 + 
     &            2*hb2*(Abr**2*ctt2 + pMT2*stt2) + 
     &            (MA02 + MSb12 - MSt12)*(MA02 - MSb12 + MSt12)*
     &             (hb*ht*MUEr*s2tb*(Abr*s2tt + 2*pMT*stt2) + 
     &               2*(ctb2*ht2*MUEr2*stt2 + 
     &                  hb2*stb2*
     &                   (Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2)))*x11
     &            )/MSt12 + 
     &         (2*ctt2*hb2*pMT2 + (5 + c2tt)*ht2*MUEr2 + 
     &            4*c2tt2*ht2*MUEr2 + 2*Abr**2*hb2*stt2 + 
     &            (MA02 + MSb12 - MSt22)*(MA02 - MSb12 + MSt22)*
     &             (hb*ht*MUEr*s2tb*(2*ctt2*pMT - Abr*s2tt) + 
     &               2*(ctb2*ctt2*ht2*MUEr2 + 
     &                  hb2*stb2*
     &                   (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))))*
     &             x12)/MSt22 + 
     &         (1 + c2tt2)*ht2*
     &          (-3*MA02 + 8*MSt12 - (MA02 - MSt12)**3/MSt22**2 + 
     &            (6*(MA02 - MSt12)*MSt12)/MSt22 - 2*MSt22 + 
     &            ((MA02 - MSt22)*(MA02 + MSt22))/MSt12)*MUEr2*x15)+
     &        ht2*MUEr2*s2tt2*
     &       ((s2tt2*(MA02 - MA02**2*(MA02 - 6*MSt12)*x16))/
     &          MSt12**2 - 
     &         ((1 + c2tt2)*MSt12 + 
     &            MA02**2*(MA02 - 6*MSt22)*s2tt2*x17)/MSt22**2) + 
     &      (2*(Abr**2*hb2*s2tt2*
     &            (MSt2diff + 
     &              3*c2tt*
     &               (2*MA02 - 2*ctb2*MSb22 + MSt2sum - 
     &                 2*MSb12*stb2)) + 
     &           2*Abr*s2tt*(-2 + 3*s2tt2)*
     &            (hb*ht*MSb2diff*MUEr*s2tb + 
     &              hb2*pMT*
     &               (2*MA02 - 2*ctb2*MSb22 + MSt2sum - 
     &                 2*MSb12*stb2)) + 
     &           s2tt2*(3*(1 - 2*c2tt2)*ht2*MSt2diff*MUEr2 + 
     &              hb2*(2*MA02*
     &                  (MSt2diff + 3*c2tt*MSt2sum - 3*c2tt*pMT2) + 
     &                 pMT2*
     &                  (6*c2tt*ctb2*MSb22 + MSt2diff - 
     &                    3*c2tt*MSt2sum + 6*c2tt*MSb12*stb2)) - 
     &              3*c2tt*
     &               (2*hb*ht*MSb2diff*pMT*MUEr*s2tb + 
     &                 ht2*MUEr2*
     &                  (2*MA02 - 2*ctb2*MSb12 + MSt2sum - 
     &                    2*MSb22*stb2))))*(x18 - x19))/MSt2diff**2
     &      )

        lc3 = (ht2*(Atr**2*
     &         (c2tt2*(11*MSt12 + 5*MSt22) + MSt2diff) - 
     &        c2tt2*(8*MA02*MUEr2 - 
     &           MSt2sum*(8*MUEr2 + MSt2sum*(4 + 4*Nc))))*s2tt2)/
     &    MSt2diff**2

        lc23 = 1/4.D0*((hb2*
     &          (2*MSt2diff**3 - c2tt*MSt2diff**2*MSt2sum + 
     &            8*MSt2diff*MSt2sum**2 + 15*c2tt*MSt2sum**3 - 
     &            2*((1 + c2tb)*MSb22**2*
     &                (MSt2diff + 3*c2tt*MSt2sum) + 
     &               (c2tt*MSt2diff**2 + MSt2diff*MSt2sum)*MUEr2)+
     &              (4*MSt2diff*MSt2sum - 30*c2tt*MSt2sum**2)*
     &             pMT2 - c2tt*
     &             (12*MSt2sum**2*MUEr2 + 6*MSt2diff**2*pMT2) + 
     &            2*(MSb22*
     &                (c2tt*
     &                   (2*MSt2diff**2 - 
     &                     3*MSt2sum*
     &                      ((2 - c2tb)*MSt2sum - 2*MUEr2 - 6*pMT2)
     &                     ) - 
     &                  MSt2diff*
     &                   (4*MSt2sum - MUEr2 - 
     &                     c2tb*(MSt2sum + MUEr2) + 2*pMT2)) + 
     &               MSb12*
     &                (c2tt*
     &                   (2*MSt2diff**2 - 
     &                     3*MSt2sum*
     &                      ((2 + c2tb)*MSt2sum - 2*MUEr2 - 6*pMT2)
     &                     ) - 
     &                  MSt2diff*
     &                   (4*MSt2sum - MUEr2 + 
     &                     c2tb*(MSt2sum + MUEr2) + 2*pMT2) + 
     &                  2*MSb22*(MSt2diff - 6*c2tt*pMT2)))) - 
     &         4*hb*ht*MSb2diff*MSt2diff*MUEr*pMT*s2tb)*s2tt2)/
     &     MSt2diff

        lc25 = -(1/4.D0*(ht2*
     &         (MSb22**2*(Atr**2*(1 + c2tt - c2tb*(1 + c2tt)) + 
     &              pMT2 - c2tt*pMT2 - c2tb*(pMT2 - c2tt*pMT2)) + 
     &           Atr*(2*c2tb*MSb22**2 - 
     &              MSb22*(2*MSb22 - 12*MSt22) - 12*MSt22**2)*pMT*
     &            s2tt)*s2tt2)/MSt22**2)

        lc21 = -(1/8.D0*((1 + c2tb)*ht2*MSb12**2*
     &         ((MSt2diff**2 + 2*c2tt*MSt2diff*MSt2sum + 
     &              MSt2sum**2)*pMT2 + 
     &           Atr*(Atr*(MSt2diff**2 - 2*c2tt*MSt2diff*MSt2sum + 
     &                 MSt2sum**2) + 4*MSt2diff*MSt2sum*pMT*s2tt))*
     &         s2tt2)/(MSt12**2*MSt22**2))

        lc22 = -(1/32.D0*
     &      (ht2*(Atr**2*(4*(1 - c2tb)*MSb22**2*
     &               (MSt2diff**2 - 2*c2tt*MSt2diff*MSt2sum + 
     &                 MSt2sum**2) - 
     &              (3 - c2tt)*
     &               (16*MSb22*MSt12*MSt22**2 - 
     &                 (MSt2diff**2 - MSt2sum**2)**2)) + 
     &           4*(1 - c2tb)*MSb22**2*
     &            (MSt2diff**2 + 2*c2tt*MSt2diff*MSt2sum + 
     &              MSt2sum**2)*pMT2 + 
     &           2*Atr*(48*MSb22*MSt12**2*MSt22 + 
     &              8*(1 - c2tb)*MSb22**2*MSt2diff*MSt2sum - 
     &              3*(MSt2diff**2 - MSt2sum**2)**2)*pMT*s2tt)*
     &         s2tt2)/(MSt12**2*MSt22**2))

        lc11 = -(s2tt2*
     &      (1/2.D0*(hb*ht*MSb2diff*MUEr*s2tb*
     &           (pMT - c2tt*pMT - Atr*s2tt)) + 
     &        (hb2*(1/2.D0*(MSt22*
     &                ((1 + c2tb)*MSb22**2*
     &                   (MSt2diff + 3*c2tt*MSt2sum) + 
     &                  (1 - c2tt)*MSt12**2*(MSt22 - MUEr2) - 
     &                  MSt22*
     &                   ((5 + 7*c2tt)*MSt22**2 - 
     &                     MSt12*
     &                      (4*MSt22 - 
     &                       c2tt*(16*MSt22 - 4*MUEr2 - 6*pMT2)) + 
     &                     2*MSt2diff*pMT2 - 
     &                     MSt22*(MUEr2 + c2tt*(7*MUEr2 + 18*pMT2))
     &                     ))) + 
     &             1/4.D0*(-(MSb22*
     &                   (2*
     &                      (MSt12*MSt22*
     &                       ((3 - c2tb*(1 - 3*c2tt) - 13*c2tt)*
     &                       MSt22 - MUEr2 - 
     &                       c2tb*(1 - c2tt)*MUEr2 + 5*c2tt*MUEr2)+
     &                        (1 - c2tt)*MSt12**2*
     &                       (MSt22 + MUEr2*(-1/2.D0 + 1/2.D0*c2tb)))-
     &                       MSt22*
     &                      (8*(1 + c2tt)*MSt22**2 + 
     &                       (4*MSt2diff - 12*c2tt*MSt2sum)*pMT2 - 
     &                       MSt22*
     &                       ((3 + 13*c2tt)*MUEr2 + 24*c2tt*pMT2 + 
     &                       2*c2tb*
     &                       (MSt22 + 
     &                       c2tt*(3*MSt22 - 1/2.D0*MUEr2) + 
     &                       1/2.D0*MUEr2))))) - 
     &                MSb12*
     &                 (-(MSt22*
     &                      (2*(4 + 4*c2tt + c2tb*(1 + 3*c2tt))*
     &                       MSt22**2 + 4*MSt2diff*pMT2 - 
     &                       2*MSt12*
     &                       ((3 + c2tb*(1 - 3*c2tt) - 13*c2tt)*
     &                       MSt22 - 
     &                       (1 - c2tb*(1 - c2tt) - 5*c2tt)*
     &                      MUEr2 + 6*c2tt*pMT2) - 
     &                       MSt22*
     &                       ((3 - c2tb*(1 - c2tt) + 13*c2tt)*
     &                       MUEr2 + 36*c2tt*pMT2))) - 
     &                   4*MSb22*
     &                    ((MSt2diff + 3*c2tt*MSt2sum)*MUEr2 + 
     &                      MSt12*
     &                       (MSt22 - c2tt*(5*MSt22 - 3*pMT2)) + 
     &                      (9*c2tt*MSt22 - MSt2diff)*pMT2 - 
     &                      MSt22**2*(3/2.D0 + 1/2.D0*c2tt) + 
     &                      1/2.D0*((1 - c2tt)*MSt12**2)) + 
     &                   2*
     &                    ((1 + c2tb)*MSb22**2*
     &                       (MSt2diff + 3*c2tt*MSt2sum) + 
     &                      (1 - c2tt)*MSt12**2*
     &                       (MSt22 + 1/2.D0*(-MUEr2 - c2tb*MUEr2))))
     &                )))/MSt2diff**2))

        lc19 = -(1/8.D0*((hb2*
     &            (MSb22*((3 - c2tb - 4*c2tt)*MSt2diff**3 + 
     &                 3*c2tt*MSt2sum**2*
     &                  ((3 - c2tb)*MSt2sum - 4*MUEr2 - 8*pMT2) + 
     &                 MSt2diff**2*
     &                  ((8 - c2tt - c2tb*(2 + 3*c2tt))*MSt2sum - 
     &                    2*
     &                     (MUEr2 + c2tb*(MUEr2 - c2tt*MUEr2) - 
     &                       2*pMT2 + c2tt*(MUEr2 + 6*pMT2))) + 
     &                 MSt2diff*MSt2sum*
     &                  ((5 + 12*c2tt - c2tb*(1 + 6*c2tt))*
     &                     MSt2sum - 
     &                    4*
     &                     (MUEr2 - pMT2 + c2tt*(3*MUEr2 + 9*pMT2))
     &                    )) + 
     &              MSb12*((3 + c2tb - 4*c2tt)*MSt2diff**3 - 
     &                 4*(1 + c2tb)*MSb22**2*
     &                  (MSt2diff + 3*c2tt*MSt2sum) + 
     &                 3*c2tt*MSt2sum**2*
     &                  ((3 + c2tb)*MSt2sum - 4*MUEr2 - 8*pMT2) + 
     &                 MSt2diff**2*
     &                  ((8 - c2tt + c2tb*(2 + 3*c2tt))*MSt2sum - 
     &                    2*
     &                     ((1 - c2tb*(1 - c2tt) + c2tt)*MUEr2 - 
     &                       2*(1 - 3*c2tt)*pMT2)) - 
     &                 4*MSb22*
     &                  ((1 - 2*c2tt)*MSt2diff**2 + 
     &                    3*c2tt*MSt2sum*
     &                     (MSt2sum - 2*MUEr2 - 4*pMT2) + 
     &                    2*MSt2diff*
     &                     (MSt2sum - MUEr2 + pMT2 - 3*c2tt*pMT2))+
     &                   MSt2diff*MSt2sum*
     &                  ((5 + 12*c2tt + c2tb*(1 + 6*c2tt))*
     &                     MSt2sum - 
     &                    4*
     &                     (MUEr2 - pMT2 + c2tt*(3*MUEr2 + 9*pMT2))
     &                    )) + 
     &              MSt22*(4*(1 + c2tb)*MSb22**2*
     &                  (MSt2diff + 3*c2tt*MSt2sum) - 
     &                 2*(2*(1 - c2tt)*MSt2diff**3 + 
     &                    6*c2tt*MSt2sum**2*
     &                     (MSt2sum - MUEr2 - 2*pMT2) + 
     &                    MSt2diff**2*
     &                     (5*MSt2sum + 
     &                       c2tt*(MSt2sum - 2*MUEr2 - 6*pMT2) + 
     &                       2*pMT2) + 
     &                    MSt2diff*MSt2sum*
     &                     (3*(1 + 3*c2tt)*MSt2sum - 
     &                       2*
     &                       (MUEr2 - pMT2 + 
     &                       c2tt*(3*MUEr2 + 9*pMT2)))))) + 
     &           4*hb*ht*MSb2diff*MSt2diff**2*MUEr*s2tb*
     &            (pMT - c2tt*pMT - Atr*s2tt))*s2tt2)/
     &       (MSb2diff*MSt2diff**2))

        lc24 = -(1/8.D0*(MSb12**2*
     &         (2*Atr**2*(1 + c2tb)*(1 + c2tt)*ht2*MSt2diff**2 - 
     &           4*(1 - c2tb)*hb2*MSt22**2*
     &            (MSt2diff + 3*c2tt*MSt2sum) + 
     &           (1 + c2tb)*ht2*MSt2diff**2*
     &            ((2 - 2*c2tt)*pMT2 - 4*Atr*pMT*s2tt))*s2tt2)/
     &       (MSt22**2*MSt2diff**2))

        lc13 = -((ht2*s2tt2*
     &        (Atr*(c2tb*MSb22**2 - MSb22*(MSb22 - 3*MSt22) - 
     &             3*MSt22**2)*pMT*s2tt + 
     &          1/2.D0*((1 - c2tb)*MSb22**2*
     &             (Atr**2*(1 + c2tt) + pMT2 - c2tt*pMT2)) - 
     &          1/2.D0*((1 - c2tb)*MSb22**2*
     &             (-(Atr*pMT*s2tt) + 
     &               1/2.D0*(Atr**2*(1 + c2tt) + pMT2 - c2tt*pMT2))))
     &        )/MSt22**2)

        lc12 = MSb12**2*s2tt2*
     &    (1/2.D0*((1 - c2tb)*hb2*(MSt2diff + 3*c2tt*MSt2sum))/
     &        MSt2diff**2 - 
     &      ((1 + c2tb)*ht2*
     &         (1/2.D0*(Atr**2*(1 + c2tt) + pMT2) + 
     &           1/2.D0*(-(c2tt*pMT2) - Atr*pMT*s2tt + 
     &              1/2.D0*(-(Atr**2*(1 + c2tt)) - pMT2 + c2tt*pMT2))
     &           ))/MSt22**2)

        lc4 = -(s2tt2*((ht2*(MSt2sum - 2*(MUEr2 + pMT2)))/
     &         MSt2diff - (ht2*
     &            ((2 + c2tt)*MSt12**2 - 
     &              MSt22*(2*MSt22 - c2tt*(MSt22 - 6*MUEr2) - 
     &                 4*MUEr2) + 
     &              2*MSt12*(c2tt*(2*MSt22 - 3*MUEr2) - 2*MUEr2) + 
     &              2*MSt2diff*pMT2) + 
     &           2*hb2*(MSt2diff*MUEr2 + 1/2.D0*MSt12**2 - 
     &              1/2.D0*MSt22**2 - 
     &              c2tt*(-3*MSt2sum*MUEr2 + 
     &                 MSt12*(2*MSt22 + 1/2.D0*MSt12) + 1/2.D0*MSt22**2
     &                 )))/MSt2diff**2))

        lc14 = 1/4.D0*(hb2*s2tt*
     &        (4*(MA02*(MSt2diff + 3*c2tt*MSt2sum) + 
     &             Abr**2*(MSt2diff + c2tt*(6*MSt22 + 3*MSt2sum)))*
     &           s2tt - 8*Abr*pMT*
     &           (2*c2tt2*MSt12 + MSt22*(6 - 10*s2tt2) + 
     &             MSt2diff*s2tt2)))/MSt2diff**2 - 
     &    ht2*((6*Atr*(1 - 3*c2tt2)*pMT*s2tt)/MSt2diff + 
     &       1/4.D0*((4*MSt12*MSt22**2*Nc + 
     &             MSt2diff*(MA02*MUEr2 - MSt22**2*(2 + 2*Nc)))*
     &           s2tt2**2)/(MSt22**2*MSt2diff) + 
     &       s2tt2*(3 + c2tt + 5*c2tt2 - 
     &          (MSt12 + 9*(MSt22 + c2tt*MSt22))/MSt2diff + 
     &          4*c2tt2*Nc - 
     &          (pMT2 + c2tt*(3*Atr**2 - 3*(MUEr2 + pMT2)))/
     &           MSt2diff + 1/2.D0*(-1 - c2tt) + 
     &          ((c2tb*MSb2diff - MSb2sum)*
     &             (1/4.D0*(Atr**2*(1 + c2tt) + pMT2 - c2tt*pMT2) - 
     &               1/2.D0*(Atr*pMT*s2tt)))/MSt22**2 + 
     &          1/2.D0*(2*c2tt2*MUEr2 + (-5 + c2tt)*pMT2 - 
     &              Atr**2*(1 + c2tt + s2tt2))/MSt22) - 
     &       MSt22*((24*Atr*(1 - 3*c2tt2)*pMT*s2tt - 
     &             (16*c2tt2*MSt12 + 
     &                4*(3*c2tt*MSt22 - 
     &                   4*c2tt2*(MSt22 + MSt2diff*Nc)) + 
     &                3*c2tt*
     &                 (4*Atr**2 - c2tb*MSb2diff + MSb2sum - 
     &                   2*MUEr2 - 4*pMT2))*s2tt2)/MSt2diff**2 + 
     &          1/2.D0*(s2tt2*
     &              (Atr**2*(1 + c2tt2) + MSt22*(s2tt2 - Nc*s2tt2))
     &              )/(MSt12*MSt2diff))) - 
     &    s2tt2*(-(ht2*(MSt22 - MUEr2 - pMT2)**2*x27) - 
     &       (x24*(ht2*((2 + c2tt)*MSt12**2 - 
     &               MSt22*
     &                (2*MSt22 - c2tt*(MSt22 - 6*MUEr2) - 4*MUEr2)+
     &                 2*(MSt12*
     &                   (c2tt*(2*MSt22 - 3*MUEr2) - 2*MUEr2) + 
     &                  MSt2diff*pMT2)) + 
     &            2*hb2*(MSt2diff*MUEr2 + 1/2.D0*MSt12**2 - 
     &               1/2.D0*MSt22**2 - 
     &               c2tt*(-3*MSt2sum*MUEr2 + 
     &                  MSt12*(2*MSt22 + 1/2.D0*MSt12) + 
     &                  1/2.D0*MSt22**2))))/MSt2diff**2 + 
     &       (ht2*((MSt2sum - 2*(MUEr2 + pMT2))*x20 + 
     &            1/4.D0*((13 + c2tt2)*MUEr2 - 
     &               (MSt22*(MUEr2 + c2tt2*(4*MSt22 + MUEr2)))/
     &                MSt12 - 
     &               2*(Atr**2*(9 + 3*c2tt2) - 
     &                  c2tb*(MSb2diff + 3*c2tt*MSb2diff) + 
     &                  MSb2sum + 3*c2tt*MSb2sum + 
     &                  MSt22*(s2tt2 + 3*Nc*s2tt2)))))/MSt2diff)

        lc20 = -(1/64.D0*
     &      (Abr*hb2*(MSt2diff**2 - MSt2sum**2)**2*
     &          (pMT*s2tt*(32 - 48*s2tt2) - 24*Abr*c2tt*s2tt2) - 
     &         ht2*(4*Atr**2*
     &             (-((2*c2tb*MSb2diff - 2*MSb2sum)*MSt2diff*
     &                  (MSt2diff**2 - 2*c2tt*MSt2diff*MSt2sum + 
     &                    MSt2sum**2)) + 
     &               8*MSt12**2*MSt22*
     &                (3*MSt2diff - c2tt*(5*MSt2diff + 6*MSt2sum)))
     &              *s2tt2 + 
     &            8*Atr*pMT*s2tt*
     &             (6*(1 - 3*c2tt2)*MSt2sum**4 + 
     &               3*MSt2diff**4*(3 - 7*c2tt2 - 2*s2tt2) - 
     &               MSt2diff*
     &                (3*MSt2diff**2*MSt2sum - 3*MSt2sum**3)*s2tt2+
     &                 MSt2diff**2*MSt2sum*
     &                ((-4*c2tb*MSb2diff + 4*MSb2sum)*s2tt2 + 
     &                  MSt2sum*(-18 + 42*c2tt2 + 9*s2tt2))) - 
     &            s2tt2*(12*c2tt*MSt2sum**4*(MSt2sum - 2*pMT2) + 
     &               2*(MSt2diff**4*
     &                   (6*c2tt*MSt2sum + 
     &                     4*(MUEr2 + c2tt2*MUEr2 - 2*c2tt*pMT2) - 
     &                     MSt2sum**2*
     &                      (2*MSt2sum - 4*(MUEr2 + pMT2))*
     &                      (x26 - x27)) + 
     &                  MSt2diff**6*(MSt2sum - 2*(MUEr2 + pMT2))*
     &                   (x26 - x27)) - 
     &               2*MSt2diff**2*MSt2sum*
     &                (4*((1 + c2tt2)*MA02*MUEr2 + 
     &                     c2tt*
     &                      (3*MSt2sum**2 - 
     &                       (2*c2tb*MSb2diff - 2*MSb2sum + 
     &                       5*MSt2sum)*pMT2)) - 
     &                  MSt2sum**3*(MSt2sum - 2*(MUEr2 + pMT2))*
     &                   (x26 - x27)) - MSt2diff**7*(x26 + x27) - 
     &               MSt2diff**5*
     &                (14 + 
     &                  2*(1 + c2tt2*(8 + 16*Nc) + 8*Nc*s2tt2) + 
     &                  (4*pMT2**2 + MUEr2*(4*MUEr2 + 8*pMT2) - 
     &                     MSt2sum*(MSt2sum + 4*(MUEr2 + pMT2)))*
     &                   (x26 + x27)) - 
     &               MSt2diff*MSt2sum**2*
     &                (8*MSt2sum*(5*pMT2 + MUEr2*s2tt2) - 
     &                  4*((2*c2tb*MSb2diff - 2*MSb2sum)*pMT2 + 
     &                     MA02*MUEr2*(2 + s2tt2)) + 
     &                  MSt2sum**3*(MSt2sum - 4*(MUEr2 + pMT2))*
     &                   (x26 + x27) + 
     &                  2*MSt2sum**2*
     &                   (7 + c2tt2*(17 + 16*Nc) + 
     &                     (5 + 4*Nc)*s2tt2 + 
     &                     (2*pMT2**2 + MUEr2*(2*MUEr2 + 4*pMT2))*
     &                      (x26 + x27))) + 
     &               MSt2diff**3*
     &                (8*(c2tb*MSb2diff - MSb2sum)*pMT2 - 
     &                  4*MSt2sum*
     &                   (3*c2tt2*MUEr2 - 10*pMT2 + 
     &                     MUEr2*(1 - s2tt2)) + 
     &                  4*MA02*MUEr2*(2 + s2tt2) + 
     &                  MSt2sum**3*(MSt2sum - 8*(MUEr2 + pMT2))*
     &                   (x26 + x27) + 
     &                  4*MSt2sum**2*
     &                   (7 + c2tt2*(13 + 16*Nc) + 
     &                     (3 + 6*Nc)*s2tt2 + 
     &                     (2*pMT2**2 + MUEr2*(2*MUEr2 + 4*pMT2))*
     &                      (x26 + x27))))))/
     &       (MSt12**2*MSt22**2*MSt2diff))

        s2tt2F3t = (lc14 + lc12/(MSb12 - MSt22) + 
     &       (lc13 + lc11/(MSb12 - MSt22))/(MSb22 - MSt22))*
     &     Log(MSt22/MSt12) + 
     &    Log(MSt12/Q)*(lc20 + lc21/(MSb12 - MSt12) + 
     &       (lc22 + lc23/(MSb12 - MSt12))/(MSb22 - MSt12) + 
     &       ((lc19 + lc24)/((MSb12 - MSt12)*(MSb12 - MSt22)) - 
     &          (lc19 - lc25)/((MSb22 - MSt12)*(MSb22 - MSt22)))*
     &        MSt2diff + lc3*Log(MSt22/Q)) + 
     &    lc4*Log(MSt12/MSt22)*Log(Q)

        s2tt2F3t = s2tt2F3t +
     &    tmp12/MSt2diff**2 + tmp8*(x21 - x10) + tmp9*(x22 - x10)

        s2tt2F3t = s2tt2F3t + 1/8D0*(ht2*s2tt2*
     &       (MA02*MSt2diff*MUEr2*x17*
     &          (2*MSt12**2*
     &             (2*MSt12*
     &                ((3 + 11*c2tt2)*MSt12*MSt22 - 
     &                  (1 + c2tt2)*MSt12**2*
     &                   (3 + (2*MSt12 - 5*MSt22)*MSt2diff*x15) + 
     &                  MSt22**2*
     &                   (4 - 3*MSt22*MSt2diff*x15 + 
     &                     c2tt2*(12 - 3*MSt22*MSt2diff*x15)))*x39+
     &                 MA02*
     &                (6*MSt22**3*MSt2diff*s2tt2*x16*x38 + 
     &                  (3*(1 + c2tt2)*MSt12**2*
     &                      (1 + 2*MSt12*MSt2diff*x15) - 
     &                     MSt12*MSt22*
     &                      (5 + 2*MSt12*MSt2diff*x15 + 
     &                       c2tt2*(13 + 2*MSt12*MSt2diff*x15)))*
     &                   x39 + 
     &                  2*MSt12*MSt22**2*MSt2diff*
     &                   (2*x15*(x39 + c2tt2*x39) + 
     &                     3*s2tt2*x17*x40))) + 
     &            MA02**2*MSt2diff*
     &             (MA02*(MSt22**3*s2tt2*x16*x38 + 
     &                  MSt12**3*
     &                   (2*x15*(x39 + c2tt2*x39) + s2tt2*x17*x40))
     &                 - 4*MSt12*
     &                (2*MSt22**3*s2tt2*x16*x38 + 
     &                  MSt12**2*
     &                   ((MSt12 + MSt2sum + 
     &                       c2tt2*(MSt12 + MSt2sum))*x15*x39 + 
     &                     2*MSt22*s2tt2*x17*x40)))) - 
     &         2*MSt12**2*(MSt2diff*x17*
     &             (2*MSt22**3*(MUEr2 + pMT2)*
     &                (2*MSt22 - 2*(MUEr2 + pMT2) + 
     &                  MSt2diff*
     &                   (1 - 
     &                     (pMT2**2 + MUEr2*(MUEr2 + 2*pMT2))*x26))
     &                 *x41 + 
     &               MSt12**3*MSt2diff*
     &                ((1 + c2tt2)*MUEr2*(3 - MSt2diff**2*x15)*
     &                   x39 + 2*MSt22**3*x26*x41)) - 
     &            MSt12**2*
     &             (MUEr2*(2*MSt22*MSt2diff**2*
     &                   (5 - 2*MSt2diff**2*x15 + 
     &                     c2tt2*(9 - 2*MSt2diff**2*x15))*x17*x39-
     &                    (2*MSt22 + MSt2diff)*s2tt2*x40) + 
     &               2*MSt22**2*MSt2diff*x17*
     &                (MSt22*
     &                   (3*MSt2diff*(MUEr2 + pMT2)*x26*x41 - 
     &                     2*x42) + 2*(MUEr2 + pMT2)*x42)) + 
     &            MSt12*MSt22*
     &             (MSt22*(MUEr2*
     &                   (MSt2diff**2*
     &                      (7 + MSt2diff**2*x15 + 
     &                       c2tt2*(31 + MSt2diff**2*x15))*x17*x39-
     &                       2*(1 - 9*c2tt2)*x40) + 
     &                  2*MSt2diff*(MUEr2 + pMT2)*
     &                   (MSt2diff + 2*(MUEr2 + pMT2))*x17*x42) + 
     &               MSt2diff*
     &                (4*c2tt2*MUEr2*x40 - 
     &                  2*MSt22**2*x17*
     &                   ((2*MSt22 - 2*(MUEr2 + pMT2) + 
     &                       MSt2diff*
     &                       (1 + 
     &                       (-3*pMT2**2 - 
     &                       MUEr2*(3*MUEr2 + 6*pMT2))*x26))*x41 + 
     &                     (MSt2diff + 2*(MUEr2 + pMT2))*x42))))))/
     &     (MSt12**3*MSt22**3*MSt2diff**2*x17)

        s2tt2F3t = s2tt2F3t + 1/8D0*(-2*(2*MSt22*
     &         ((2*Abr*hb2 + 6*Atr*ht2)*pMT*s2tt*
     &            (MSt22**2*s2tt2 - 
     &              MSt12*(MSt22*(16 - 22*s2tt2) - MSt12*s2tt2)) + 
     &           s2tt2*(2*hb2*
     &               (MSt12*MSt22*
     &                  (3*MSt22 + 
     &                    c2tt*
     &                     (5*MA02 + 6*MSt2sum - 5*MUEr2 - 
     &                       11*pMT2 + 5*MSb12*stb2)) + 
     &                 Abr**2*
     &                  (11*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &                    MSt12**2*stt2) + 
     &                 MSt22**2*
     &                  (ctt2*(MA02 - MUEr2 + MSb12*stb2) + 
     &                    pMT2*stt2) - 
     &                 MSt12**2*
     &                  (3*MSt22 + ctt2*pMT2 + 
     &                    (MA02 - MUEr2 + MSb12*stb2)*stt2)) - 
     &              ht2*(MSt12**3*(2 - s2tt2 - Nc*s2tt2) - 
     &                 MSt12*MSt22*
     &                  (MSt22*
     &                     (32 + Nc*(16 - 13*s2tt2) - 13*s2tt2) - 
     &                    c2tt*
     &                     (12*MSt22 - 10*MUEr2 - 22*pMT2 + 
     &                       10*MSb22*stb2)) + 
     &                 MSt12**2*
     &                  (12*c2tt*MSt22 - 2*MUEr2 + 
     &                    MSt22*
     &                     (32 + Nc*(16 - 13*s2tt2) - 13*s2tt2) - 
     &                    ctt2*(2*MUEr2 - 2*MSb22*stb2) + 
     &                    2*pMT2*stt2) - 
     &                 MSt22**2*
     &                  (2*ctt2*pMT2 + 
     &                    MSt22*(2 - s2tt2 - Nc*s2tt2) + 
     &                    2*MSb22*stb2*stt2 - MUEr2*(2 + 2*stt2)))-
     &                2*(Atr**2*ht2*
     &                  (MSt12**2 + ctt2*MSt12**2 + 
     &                    11*c2tt*MSt12*MSt22 - MSt22**2*(1 + stt2)
     &                    ) - 
     &                 ctb2*
     &                  (hb2*MSb22*
     &                     (5*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &                       MSt12**2*stt2) - 
     &                    ht2*MSb12*
     &                     (ctt2*MSt12**2 + 5*c2tt*MSt12*MSt22 - 
     &                       MSt22**2*stt2))))) - 
     &        s2tt2*x10*(8*hb2*MSt22*MUEr2*
     &            (2*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &              MSt12**2*stt2) + 
     &           ht2*(2*MSt22*
     &               (5*MSt22**2*MUEr2 - 
     &                 c2tt*
     &                  (2*MSt22**2 + MSt12*(2*MSt12 + 8*MSt22))*
     &                  MUEr2 - 
     &                 MSt12*(4*MSt22*MSt2diff + 5*MSt12*MUEr2) + 
     &                 MSt2diff*MSt2sum*(MUEr2 + 4*pMT2)) - 
     &              (1 + c2tt2)*MSt2diff*
     &               (MA02**2*
     &                  (MA02*MSt12 - 3*MSt12**2 - MSt22**2) + 
     &                 (3*MA02*MSt12 - MSt12*(MSt12 - 4*MSt22) + 
     &                    MSt22**2)*MSt2diff**2)*MUEr2*x15))))/
     &    (MSt12*MSt22**2*MSt2diff)

        s2tt2F3t = s2tt2F3t + 1/8D0*(4*
     &       ((MSt2diff*s2tt2*
     &            ((ht2*((5*MSt12**2 - 5*MSt22**2 + 
     &                       c2tt*
     &                       (2*MSt22**2 + 
     &                       MSt12*(2*MSt12 + 8*MSt22)))*MUEr2 - 
     &                    MSt2diff*MSt2sum*pMT2) - 
     &                 4*hb2*MUEr2*
     &                  (2*c2tt*MSt12*MSt22 + ctt2*MSt22**2 - 
     &                    MSt12**2*stt2))*x24 + 
     &              ht2*MSt22*MSt2diff*(MSt12 - MUEr2 - pMT2)**2*
     &               (MSt12*(2*x10 - x20 - x24) + 
     &                 (MUEr2 - pMT2)*(x20 - x24))*x26))/MSt22 + 
     &         (6*Abr**2*c2tt*hb2*s2tt2*
     &             (stb2*x13*x28 + ctb2*x11*x30) - 
     &            2*Abr*s2tt*(2 - 3*s2tt2)*
     &             (hb*ht*MUEr*s2tb*(x13*x28 - x11*x30) + 
     &               2*hb2*pMT*(stb2*x13*x28 + ctb2*x11*x30)) - 
     &            6*c2tt*s2tt2*
     &             (hb*ht*MUEr*pMT*s2tb*(x13*x28 - x11*x30) + 
     &               stb2*(hb2*pMT2*x13*x28 + ht2*MUEr2*x11*x30) + 
     &               ctb2*(ht2*MUEr2*x13*x28 + hb2*pMT2*x11*x30)))/
     &          (x11*x13)))/(MSt12*MSt2diff**2) - 
     &    ht2*s2tt2*(-x24 + 
     &       1/4.D0*(MUEr2*(MA02**3*
     &              (MSt22**2*s2tt2*x16*(x10 - x18) - 
     &                MSt12**2*
     &                 (x15*(x18 + c2tt2*x18) - 
     &                   s2tt2*x17*(x10 - x19))) + 
     &             MA02**2*MSt12*
     &              (3*(1 + c2tt2)*MSt12**2*x15*x18 + 
     &                2*MSt12*MSt22*
     &                 (-(s2tt2*x17*(3*x10 - 3*x19)) + 
     &                   x15*(x18 + c2tt2*(x18 - x19) - x19)) - 
     &                MSt22**2*
     &                 (s2tt2*x16*(6*x10 - 6*x18) - 
     &                   x15*(x19 + c2tt2*x19))) + 
     &             (1 + c2tt2)*MSt12*x15*
     &              (MSt12**3*(MSt12 - 4*MSt22)*x18 - 
     &                MSt22**4*x19 + 
     &                MSt12**2*MSt22*
     &                 (-(MSt22*(x18 - 5*x19)) + 2*MSt2diff*x19) + 
     &                MSt12*
     &                 (2*MSt22**2*
     &                    (MSt2diff*(x18 - 2*x19) + 
     &                      MSt22*(2*x18 - 2*x19)) - 
     &                   MA02*
     &                    (3*MSt22**2*x18 + 
     &                      MSt12*
     &                       (3*MSt12*x18 - MSt22*(2*x18 + 4*x19)))
     &                   ))))/(MSt12**2*MSt22**2) - 
     &       1/2.D0*((2*MSt12*MSt22 - MSt2sum*(MUEr2 + 3*pMT2))*x20)/
     &         (MSt12*MSt22))

        s2tt2F3t = s2tt2F3t + 
     &   1/4.D0*(((2*(MA02 + 3*c2tt*MA02 + MSb12 + 3*c2tt*MSb12 - 
     &                2*c2tt*MSt22 - ctt2*MSt2sum)*s2tt2*
     &              (ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb + 
     &                hb2*pMT2*stb2) + 
     &             Abr*s2tt*
     &              (hb*ht*MUEr*s2tb*
     &                 (4*MSb12 - 4*MSt22 + MA02*(4 - 6*s2tt2) - 
     &                   6*MSb12*s2tt2 + 4*MSt22*s2tt2 + 
     &                   MSt2sum*s2tt2) + 
     &                2*hb2*pMT*
     &                 (4*MSb12 - 4*c2tt2*MSt22 + 
     &                   MA02*(4 - 6*s2tt2) - 6*MSb12*s2tt2 + 
     &                   MSt2sum*s2tt2)*stb2) + 
     &             2*Abr**2*hb2*s2tt2*stb2*
     &              (MA02 - 3*c2tt*MA02 + MSb12 - 3*c2tt*MSb12 + 
     &                2*c2tt*MSt22 - MSt2sum*stt2))*x29 + 
     &          (-(Abr*s2tt*
     &                (hb*ht*MUEr*s2tb*
     &                   (4*MSb22 - 4*MSt22 + MA02*(4 - 6*s2tt2) - 
     &                     6*MSb22*s2tt2 + 4*MSt22*s2tt2 + 
     &                     MSt2sum*s2tt2) + 
     &                  2*ctb2*hb2*pMT*
     &                   (-4*MSb22 + 4*c2tt2*MSt22 + 
     &                     6*MSb22*s2tt2 - MSt2sum*s2tt2 + 
     &                     MA02*(-4 + 6*s2tt2)))) + 
     &             2*(MA02 + 3*c2tt*MA02 + MSb22 + 3*c2tt*MSb22 - 
     &                2*c2tt*MSt22 - ctt2*MSt2sum)*s2tt2*
     &              (ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + 
     &                ht2*MUEr2*stb2) + 
     &             2*Abr**2*ctb2*hb2*s2tt2*
     &              (MA02 - 3*c2tt*MA02 + MSb22 - 3*c2tt*MSb22 + 
     &                2*c2tt*MSt22 - MSt2sum*stt2))*x31)/MSt2diff+
     &         2*ht2*(-MSt22 + pMT2 + MUEr2)**2*s2tt2*x27*
     &        (MSt22*(2*x10 - x20 - x24 - x42) + 
     &          MUEr2*(x20 - x24 + x42) + pMT2*(-x20 + x24 + x42)))/
     &     MSt22

        s2tt2F3t = s2tt2F3t + 
     &   (((ht2*s2tt2*x41)/x26 + 
     &         1/4.D0*((ctb2*ht2*MUEr2*s2tt2*
     &                ((-2 + 6*c2tt)*(MA02 + MSb12) - 
     &                  4*c2tt*MSt12 + 2*MSt2sum*stt2) + 
     &               s2tt*(-2*Abr**2*hb2*
     &                   (MA02 + 3*c2tt*MA02 + MSb12 + 
     &                     3*c2tt*MSb12 - 2*c2tt*MSt12 - 
     &                     ctt2*MSt2sum)*s2tt*stb2 + 
     &                  Abr*
     &                   (hb*ht*MUEr*s2tb*
     &                      (4*MSb12 - 4*MSt12 + 
     &                       MA02*(4 - 6*s2tt2) - 6*MSb12*s2tt2 + 
     &                       4*MSt12*s2tt2 + MSt2sum*s2tt2) + 
     &                     2*hb2*pMT*
     &                      (4*MSb12 - 4*c2tt2*MSt12 + 
     &                       MA02*(4 - 6*s2tt2) - 6*MSb12*s2tt2 + 
     &                       MSt2sum*s2tt2)*stb2) + 
     &                  2*s2tt*(hb*ht*pMT*MUEr*s2tb + hb2*pMT2*stb2)*
     &                   ((-1 + 3*c2tt)*MA02 + 
     &                     (-1 + 3*c2tt)*MSb12 - 2*c2tt*MSt12 + 
     &                     MSt2sum*stt2)))*x28 + 
     &            (ht2*MUEr2*s2tt2*stb2*
     &                ((-2 + 6*c2tt)*(MA02 + MSb22) - 
     &                  4*c2tt*MSt12 + 2*MSt2sum*stt2) + 
     &               s2tt*(hb*ht*MUEr*s2tb*
     &                   (Abr*
     &                      (-4*(MA02 + MSb22 - MSt12) + 
     &                       (6*(MA02 + MSb22) - 4*MSt12 - 
     &                      MSt2sum)*s2tt2) + 
     &                     2*pMT*s2tt*
     &                      (MA02 + MSb22 + 
     &                       c2tt*(-3*(MA02 + MSb22) + 2*MSt12) - 
     &                       MSt2sum*stt2)) - 
     &                  2*ctb2*hb2*
     &                   (Abr*pMT*
     &                      (-4*(MA02 + MSb22 - c2tt2*MSt12) + 
     &                       (6*(MA02 + MSb22) - MSt2sum)*s2tt2) + 
     &                     s2tt*
     &                      (Abr**2*
     &                       (MA02 + MSb22 + 
     &                       c2tt*(3*(MA02 + MSb22) - 2*MSt12) - 
     &                       ctt2*MSt2sum) + 
     &                       pMT2*
     &                       (MA02 + MSb22 + 
     &                       c2tt*(-3*(MA02 + MSb22) + 2*MSt12) - 
     &                       MSt2sum*stt2)))))*x30))/MSt12 + 
     &      ht2*s2tt2*(-(x42/(MSt22*x27)) + 
     &         2*Atr**2*(-1 + 3*c2tt2)*x43 + 
     &         1/8.D0*(MA02*MUEr2*
     &            (((-8*c2tt2*MA02*MSt12 - 12*MSt12**2 + 
     &                   60*c2tt2*MSt12**2 + 5*MA02*MSt12*s2tt2 - 
     &                   3*MA02*MSt22*s2tt2 + 
     &                   12*MSt12*MSt2diff*s2tt2)*x38)/MSt12**3 + 
     &              ((8*c2tt2*(MA02 - 6*MSt22)*MSt22 + 
     &                   (MA02*(3*MSt12 - 5*MSt22) + 
     &                      12*MSt22*(MSt22 + MSt2diff))*s2tt2)*x40
     &                 )/MSt22**3))))/MSt2diff

        s2tt2F3t = s2tt2F3t + tmp13*(-x10 + x9) +
     &   1/4.D0*(s2tt2*(((MA02 + MSb12 - MSt12)*
     &           (2*Abr**2*ctt2*hb2*stb2 + 
     &             Abr*s2tt*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2) + 
     &             2*(ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb + 
     &                hb2*pMT2*stb2)*stt2)*x11*
     &           (MA02*(x10 - x21 - x28) - 
     &             MSb12*(x10 - x21 + x28) + 
     &             MSt12*(x10 - 2*x18 + x21 + x28)))/MSt12 + 
     &        ((MA02 + MSb12 - MSt22)*
     &           (2*ctb2*ctt2*ht2*MUEr2 + 
     &             2*ctt2*(hb*ht*pMT*MUEr*s2tb + hb2*pMT2*stb2) + 
     &             Abr*(-(hb*ht*MUEr*s2tb*s2tt) - 
     &                2*hb2*pMT*s2tt*stb2 + 2*Abr*hb2*stb2*stt2))*
     &           x12*(MA02*(x10 - x21 - x29) - 
     &             MSb12*(x10 - x21 + x29) + 
     &             MSt22*(x10 - 2*x19 + x21 + x29)))/MSt22 + 
     &        ((MA02 + MSb22 - MSt12)*
     &           (2*Abr**2*ctb2*ctt2*hb2 + 
     &             Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt + 
     &             2*(ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + 
     &                ht2*MUEr2*stb2)*stt2)*x13*
     &           (MA02*(x10 - x22 - x30) - 
     &             MSb22*(x10 - x22 + x30) + 
     &             MSt12*(x10 - 2*x18 + x22 + x30)))/MSt12 + 
     &        ((MA02 + MSb22 - MSt22)*
     &           (-2*ctt2*hb*ht*pMT*MUEr*s2tb + 
     &             Abr*hb*ht*MUEr*s2tb*s2tt + 2*ctt2*ht2*MUEr2*stb2 + 
     &             2*ctb2*hb2*
     &              (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*x14*
     &           (MA02*(x10 - x22 - x31) - 
     &             MSb22*(x10 - x22 + x31) + 
     &             MSt22*(x10 - 2*x19 + x22 + x31)))/MSt22))

        s2ttF4t = 1/8.D0*(((-((MA02 + MSb12 - MSt12)*
     &              (2*Abr*ctt2*hb2*s2tb*s2tt + 
     &                2*ctb2*hb*ht*MUEr*s2tt2 + hb2*pMT*s2tb*s2tt2))+
     &             (2*(Abr*c2tt*hb2*s2tb*s2tt + 
     &                2*ctb2*hb*ht*MUEr*s2tt2 + hb2*pMT*s2tb*s2tt2))/
     &            (MSt2diff*x11))*x28 + 
     &        ((MA02 + MSb22 - MSt12)*
     &            (2*Abr*ctt2*hb2*s2tb*s2tt + hb2*pMT*s2tb*s2tt2 - 
     &              2*hb*ht*MUEr*s2tt2*stb2) + 
     &           (-2*hb2*s2tb*(Abr*c2tt*s2tt + pMT*s2tt2) + 
     &              4*hb*ht*MUEr*s2tt2*stb2)/(MSt2diff*x13))*x30)/
     &      MSt12 + (-(((MA02 + MSb12 - MSt22)*
     &              (2*ctb2*hb*ht*MUEr*s2tt2 + 
     &                hb2*s2tb*(pMT*s2tt2 - 2*Abr*s2tt*stt2)) + 
     &             (2*(Abr*c2tt*hb2*s2tb*s2tt + 
     &                  2*ctb2*hb*ht*MUEr*s2tt2 + hb2*pMT*s2tb*s2tt2)
     &                )/(MSt2diff*x12))*x29) + 
     &        ((s2tt2*(hb2*pMT*s2tb - 2*hb*ht*MUEr*stb2)*
     &              (2 + MA02*MSt2diff*x14 + MSb22*MSt2diff*x14 - 
     &                MSt22*MSt2diff*x14) + 
     &             2*Abr*hb2*s2tb*s2tt*
     &              (c2tt - 
     &                (MA02 + MSb22 - MSt22)*MSt2diff*stt2*x14))*
     &           x31)/(MSt2diff*x14))/MSt22)

        s2ttF4t = s2ttF4t + 
     &   (hb*ht*MUEr*s2tt2*((2*MSt22 - MSt2diff)*(x10 - x19)*
     &          (x10 - x24) + (MSt2sum - 2*MUEr2)*x6 + 
     &         (-2*MSt22 + MSt2diff + 2*MUEr2)*x7) + 
     &      1/8.D0*((-x10 + x19)*
     &         ((2*hb*ht*MUEr*s2tt2*
     &               (2*ctb2*(MA02 - MSb12) + 
     &                 (-2*MSt22 + MSt2diff)*(1 + stb2)) + 
     &              hb2*s2tb*s2tt*
     &               ((2*MA02 - 2*(MSb12 + MSt22) + MSt2diff)*pMT*
     &                  s2tt + 
     &                 2*Abr*
     &                  (c2tt*(MA02 - MSb12) - ctt2*MSt22 + 
     &                    MSt12*stt2)) + 
     &              2*ht2*s2tb*s2tt*
     &               ((-2*MSt22 + MSt2diff)*pMT*s2tt + 
     &                 2*Atr*(ctt2*MSt12 - MSt22*stt2)))*
     &            (-x10 + x21) + 
     &           (s2tt2*((2*ht2*(-2*MSt22 + MSt2diff) + 
     &                    hb2*
     &                     (2*MA02 - 2*MSb22 - 2*MSt22 + MSt2diff))
     &                   *pMT*s2tb + 
     &                 2*hb*ht*MUEr*
     &                  (2*(1 + ctb2)*MSt22 - 
     &                    (1 + ctb2)*MSt2diff + 
     &                    2*(-MA02 + MSb22)*stb2)) + 
     &              2*Abr*hb2*s2tb*s2tt*
     &               (c2tt*(MA02 - MSb22) - ctt2*MSt22 + 
     &                 MSt12*stt2) + 
     &              4*Atr*ht2*s2tb*s2tt*(ctt2*MSt12 - MSt22*stt2))*
     &            (x10 - x22))) + 
     &      (-x10 + x18)*(hb*ht*MSt2sum*MUEr*s2tt2*(x10 - x24) + 
     &         1/8.D0*((s2tt2*
     &               (4*ctb2*hb*ht*(-MA02 + MSb12)*MUEr + 
     &                 hb2*(-2*MA02 + 2*MSb12 + MSt2sum)*pMT*s2tb + 
     &                 2*hb*ht*MSt2sum*MUEr*(1 + stb2)) + 
     &              2*Abr*hb2*s2tb*s2tt*
     &               (c2tt*(-MA02 + MSb12) + ctt2*MSt22 - 
     &                 MSt12*stt2) + 
     &              2*ht2*s2tb*s2tt*
     &               (-2*Atr*ctt2*MSt12 + MSt2sum*pMT*s2tt + 
     &                 2*Atr*MSt22*stt2))*(-x10 + x21)) - 
     &         1/8.D0*((s2tt2*
     &               ((2*hb2*MA02 - 2*hb2*MSb22 - hb2*MSt2sum - 
     &                    2*ht2*MSt2sum)*pMT*s2tb + 
     &                 2*hb*ht*MUEr*
     &                  (MSt2sum + ctb2*MSt2sum - 2*MA02*stb2 + 
     &                    2*MSb22*stb2)) + 
     &              2*Abr*hb2*s2tb*s2tt*
     &               (c2tt*(MA02 - MSb22) - ctt2*MSt22 + 
     &                 MSt12*stt2) + 
     &              4*Atr*ht2*s2tb*s2tt*(ctt2*MSt12 - MSt22*stt2))*
     &            (x10 - x22))) + 
     &      1/4.D0*(((2*MSb12 - MSt2sum)*s2tt2*
     &             (ht2*pMT*s2tb + 2*hb*ht*MUEr*stb2) - 
     &            2*Atr*ht2*s2tb*s2tt*
     &             (c2tt*MSb12 - ctt2*MSt12 + MSt22*stt2))*x33 + 
     &         ((2*MSb22 - MSt2sum)*
     &             (2*ctb2*hb*ht*MUEr - ht2*pMT*s2tb)*s2tt2 + 
     &            2*Atr*ht2*s2tb*s2tt*
     &             (c2tt*MSb22 - ctt2*MSt12 + MSt22*stt2))*x34 + 
     &         (-((2*MSb12 - 2*MSt22 + MSt2diff)*s2tt2*
     &               (ht2*pMT*s2tb + 2*hb*ht*MUEr*stb2)) + 
     &            2*Atr*ht2*s2tb*s2tt*
     &             (c2tt*MSb12 - ctt2*MSt12 + MSt22*stt2))*x35 - 
     &         (Atr*ht2*(2*c2tt*MSb22 + MSt2diff - c2tt*MSt2sum)*
     &             s2tb*s2tt + 
     &            (2*MSb22 - 2*MSt22 + MSt2diff)*
     &             (2*ctb2*hb*ht*MUEr - ht2*pMT*s2tb)*s2tt2)*x36 + 
     &         (Abr*c2tt*hb2*MSb2diff*s2tb*s2tt + 
     &            s2tt2*(hb2*MSb2diff*pMT*s2tb + 
     &               hb*ht*MUEr*
     &                (2*MA02 - 2*ctb2*MSb12 + MSt2sum - 
     &                  2*MSb22*stb2)))*(x18 - x19)*(-x10 + x9)))/
     &    MSt2diff

	F1b = 0

        lc3 = 1/4.D0*((2*Abr*hb2*(MA02 - MSb2sum) - 
     &          Atr*ht2*(MSb2sum - MSt12 - MSt22 - c2tt*MSt2diff))*
     &        s2tb2 + MSt2diff*
     &        (c2tb*hb*ht*MUEr*s2tb - ht2*pMT*s2tb2)*s2tt)/MSb2diff

        lc4 = 1/4.D0*(Abr*hb2*
     &        (6*MA02 + MSb2sum + c2tt*MSt2diff - MSt2sum)*s2tb2 + 
     &       MSt2diff*(c2tb*hb*ht*MUEr*s2tb + hb2*pMT*s2tb2)*s2tt)/
     &     MSb2diff

        lc5 = -(1/8.D0*(hb*ht*
     &          (MSb2diff - c2tb*(2*MA02 - MSb2sum + 2*MSt12))*
     &          MUEr*s2tb*s2tt - 
     &         s2tb2*(Abr*hb2*
     &             (2*MA02 - MSb2sum + 
     &               c2tt*(2*MA02 - MSb2sum - 2*MSt12) - 2*MSt12)+
     &              hb2*(2*MA02 - MSb2sum - 2*MSt12)*pMT*s2tt + 
     &            ht2*MSt12*(-(Atr*(4 - 4*c2tt)) - 4*pMT*s2tt)))/
     &       MSb2diff)

        lc6 = -(1/8.D0*((-(Abr*hb2*
     &               (2*MA02 - MSb2sum - 
     &                 c2tt*(2*MA02 - MSb2sum - 2*MSt22) - 2*MSt22)
     &               ) + 4*Atr*ht2*(MSt22 + c2tt*MSt22))*s2tb2 - 
     &         (hb*ht*(MSb2diff - 
     &               c2tb*(2*MA02 - MSb2sum + 2*MSt22))*MUEr*s2tb+
     &              (-(hb2*(2*MA02 - MSb2sum - 2*MSt22)) + 
     &               4*ht2*MSt22)*pMT*s2tb2)*s2tt)/MSb2diff)

        lc7 = 1/2.D0*((3*Abr*hb2 + Atr*ht2)*(2*MSb12 + MSb2sum)*
     &       s2tb2)/MSb2diff

        lc8 = 1/2.D0*((3*Abr*hb2 + Atr*ht2)*(2*MSb22 + MSb2sum)*
     &       s2tb2)/MSb2diff

        lc9 = 1/8.D0*((((MA02 - MSb12 + MSt12)*
     &             (Abr*(1 + c2tt)*hb2*s2tb2 + 
     &               (hb*ht*MUEr*(s2tb + c2tb*s2tb) + 
     &                  hb2*pMT*s2tb2)*s2tt) - 
     &            (2*(Abr*hb2*(s2tb2 + c2tt*s2tb2) + 
     &                 (c2tb*hb*ht*MUEr*s2tb + hb2*pMT*s2tb2)*s2tt)
     &               )/(MSb2diff*x11))*x28 + 
     &         ((MA02 - MSb22 + MSt12)*
     &             (Abr*(1 + c2tt)*hb2*s2tb2 - 
     &               (hb*ht*MUEr*(s2tb - c2tb*s2tb) - 
     &                  hb2*pMT*s2tb2)*s2tt) + 
     &            (2*(Abr*hb2*(s2tb2 + c2tt*s2tb2) + 
     &                 (c2tb*hb*ht*MUEr*s2tb + hb2*pMT*s2tb2)*s2tt)
     &               )/(MSb2diff*x13))*x30)/MSt12 - 
     &      (((MA02 - MSb12 + MSt22)*
     &             ((1 + c2tb)*hb*ht*MUEr*s2tb*s2tt - 
     &               hb2*s2tb2*(Abr*(1 - c2tt) - pMT*s2tt)) + 
     &            (2*(Abr*hb2*(s2tb2 - c2tt*s2tb2) - 
     &                 (c2tb*hb*ht*MUEr*s2tb + hb2*pMT*s2tb2)*s2tt)
     &               )/(MSb2diff*x12))*x29 - 
     &         ((MA02 - MSb22 + MSt22)*
     &             (Abr*(1 - c2tt)*hb2*s2tb2 + 
     &               (hb*ht*MUEr*(s2tb - c2tb*s2tb) - 
     &                  hb2*pMT*s2tb2)*s2tt) + 
     &            (2*(Abr*hb2*(s2tb2 - c2tt*s2tb2) - 
     &                 (c2tb*hb*ht*MUEr*s2tb + hb2*pMT*s2tb2)*s2tt)
     &               )/(MSb2diff*x14))*x31)/MSt22 - 
     &      (2*((Atr*(1 - c2tt)*ht2*(MSb2sum - 2*MSt12)*s2tb2 + 
     &              (hb*ht*(MSb2diff - c2tb*(MSb2sum - 2*MSt12))*
     &                  MUEr*s2tb + 
     &                 ht2*(MSb2sum - 2*MSt12)*pMT*s2tb2)*s2tt)*
     &            (x33 - x34) + 
     &           (Atr*(1 + c2tt)*ht2*(MSb2sum - 2*MSt22)*s2tb2 - 
     &              (hb*ht*(MSb2diff - c2tb*(MSb2sum - 2*MSt22))*
     &                  MUEr*s2tb + 
     &                 ht2*(MSb2sum - 2*MSt22)*pMT*s2tb2)*s2tt)*
     &            (x35 - x36)))/MSb2diff + 
     &      s2tb2*((2*Abr*hb2*
     &            (2*MSb22 + MSb2diff - 
     &              MA02*(MA02 - 6*MSb22)*MSb2diff*x46)*x47)/
     &          (MSb22**2*MSb2diff*x46) + 
     &         8*(6*Abr*hb2 + 2*Atr*ht2 - 
     &            1/4.D0*(Abr*hb2*
     &                (2*MSb12 + 
     &                  MSb2diff*(-1 + MA02*(MA02 - 6*MSb12)*x44))*
     &                x45)/(MSb12**2*MSb2diff*x44))))

        s2tbF2b = -lc9 - lc7*Log(MSb12/Q) + 
     &    lc8*Log(MSb22/Q) - 
     &    Log(MSb22/MSb12)*
     &     (lc4*Log(MA02/Q) - lc3*(Log(MSb12/Q) + Log(MSb22/Q)) - 
     &       lc5*Log(MSt12/Q) - lc6*Log(MSt22/Q))

        tmp14 = (-2*hb*ht*pMT*MUEr*(MSb2sum + 2*pMT2 - 2*MUEr2)*s2tb*
     &       (-2 + 3*s2tb2) + 
     &      s2tb2*(ctb2*(-(hb2*MSb12**2) + ht2*MSb22**2) - 
     &         hb2*MSb2diff*pMT2 - ht2*MSb2diff*pMT2 - 
     &         hb2*MSb2diff*MUEr2 - ht2*MSb2diff*MUEr2 - 
     &         c2tb*(hb2 - ht2)*
     &          (2*MSb12*MSb22 - 
     &            3*(MSb2sum + pMT2 - MUEr2)*(pMT2 + MUEr2)) - 
     &         ht2*MSb12**2*stb2 + hb2*MSb22**2*stb2))*(-x10 + x20)
     &     + hb2*s2tb2*(MSb12*
     &       (3*c2tb*(2*MA02 + 4*(Abr**2 + MSb12) - c2tt*MSt2diff + 
     &            MSt2sum - 2*(pMT2 + MUEr2)) - 
     &         16*c2tb2*MSb2diff*(1 + Nc)) + 
     &      (MSb2diff*MUEr2 + 
     &         c2tb2*(8*Abr**2*(-MA02 + MSb2sum) + 
     &            (11*MSb12 + 5*MSb22)*MUEr2 + 4*MSb2sum**2*(1 + Nc)
     &            ))*(-x10 + x22)) + 
     &   1/4.D0*((-2*hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &         (-2*MSt12*(2*(1 + ctt2)*pMT + (Abr + 2*Atr)*s2tt) + 
     &           (2*MA02 - MSb2sum)*(Abr*s2tt + 2*pMT*stt2)) - 
     &        hb2*s2tb2*(-2*MSb22*MSt12 + 6*c2tb*MSb22*MSt12 - 
     &           2*c2tt*MSb22*MSt12 + 6*c2tb*c2tt*MSb22*MSt12 - 
     &           2*Abr**2*ctt2*
     &            (MSb2diff + 
     &              c2tb*(6*MA02 - 3*(MSb2sum + 2*MSt12))) - 
     &           6*c2tb*MA02*pMT2 + 6*c2tb*c2tt*MA02*pMT2 - 
     &           MSb22*pMT2 + 3*c2tb*MSb22*pMT2 + c2tt*MSb22*pMT2 - 
     &           3*c2tb*c2tt*MSb22*pMT2 + 6*c2tb*MSt12*pMT2 - 
     &           6*c2tb*c2tt*MSt12*pMT2 + 12*c2tb*MSt12*MUEr2 + 
     &           12*c2tb*c2tt*MSt12*MUEr2 - 
     &           2*Abr*(MSb2diff + 
     &              c2tb*(6*MA02 - 3*(MSb2sum + 2*MSt12)))*pMT*s2tt+
     &             2*(1 + 3*c2tb)*MSb12*(2*ctt2*MSt12 + pMT2*stt2))
     &          - ht2*s2tb2*
     &         (-2*MSb2diff*(2*MSt12 + MUEr2)*stt2 - 
     &           6*c2tb*(4*ctt2*MSt12*pMT2 + 4*Atr*MSt12*pMT*s2tt + 
     &              (2*(2*Atr**2 + MSb2sum)*MSt12 + 
     &                 (-2*MA02 + MSb2sum + 2*MSt12)*MUEr2)*stt2)))*
     &      (-x10 + x18)) + 
     &   1/2.D0*((ctt2*(-((12*Atr**2*c2tb*ht2*MSt22 + 
     &                hb2*(-MSb12 + MSb22 + 
     &                   c2tb*
     &                    (6*MA02 - 3*(MSb12 + MSb22 + 2*MSt22)))*
     &                 pMT2 + 
     &                ht2*(MSb2diff*(2*MSt22 + MUEr2) + 
     &                   3*c2tb*
     &                    (2*MSb2sum*MSt22 - 2*MA02*MUEr2 + 
     &                      MSb2sum*MUEr2 + 2*MSt22*MUEr2)))*s2tb2)+
     &             2*hb*ht*(2*MA02 - MSb2sum)*pMT*MUEr*s2tb*
     &            (-2 + 3*s2tb2)) + 
     &        s2tb2*(Abr*hb2*
     &            (MSb2diff + c2tb*(6*MA02 - 3*MSb2sum - 6*MSt22))*
     &            pMT*s2tt - 
     &           Abr**2*hb2*
     &            (MSb2diff + c2tb*(6*MA02 - 3*MSb2sum - 6*MSt22))*
     &            stt2 + 2*MSt22*
     &            (6*Atr*c2tb*ht2*pMT*s2tt + 
     &              (-6*c2tb*ht2*pMT2 + 
     &                 hb2*
     &                  (MSb12 + 3*c2tb*MSb12 - MSb22 + 
     &                    3*c2tb*MSb22 + 6*c2tb*MUEr2))*stt2)) - 
     &        hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &         (Abr*(2*MA02 - MSb2sum)*s2tt + 
     &           MSt22*(-2*(Abr + 2*Atr)*s2tt + 4*pMT*(1 + stt2))))*
     &      (x10 - x19))

        tmp14 = tmp14 + (-2*hb*ht*pMT*MUEr*
     &       (MSb2sum - 2*pMT2 + 2*MUEr2)*s2tb*(-2 + 3*s2tb2) + 
     &      s2tb2*(ht2*(ctb2*MSb22**2 - MSb2diff*(pMT2 + MUEr2) + 
     &            c2tb*(2*MSb12*MSb22 - 
     &               3*(MSb2sum - pMT2 + MUEr2)*(pMT2 + MUEr2)) - 
     &            MSb12**2*stb2) + 
     &         hb2*(-((2 + ctb2)*MSb12**2) + 
     &            MSb22*((-1 + 3*c2tb)*pMT2 + (-5 + 3*c2tb)*MUEr2) + 
     &            3*c2tb*(-pMT2**2 + MUEr2**2) + 
     &            MSb12*(pMT2 + 5*MUEr2 + 
     &               c2tb*(-2*MSb22 + 3*(pMT2 + MUEr2))) + 
     &            MSb22**2*(2 + stb2))))*(-x10 + x24)

        tmp15 = hb2*(1/2.D0*((1 + Nc)*s2tb2**2) + 
     &     s2tb2*(-2 + c2tb - 4*c2tb2*Nc + 
     &        (ctt2*MSt12**2*MUEr2*stb2)/
     &         (MSb12**2*(-MSb12 + MSt12)) + 
     &        1/4.D0*((1 + c2tb2)*
     &            (Abr**2*(MA02 - MSb12)*MSb12**3*MSb2diff + 
     &              2*MSb22**3*
     &               (-(MSb12*MSb22) + MSb22**2 + MSb2diff**2)*MUEr2
     &              ))/(MSb12**2*MSb22**2*MSb2diff**2) + 
     &        (MUEr2*(1 + stb2 - 
     &              (MSt22**2*stb2*stt2)/(MSb12 - MSt22)**2) + 
     &           1/2.D0*(Abr**2*s2tb2))/MSb12 + 
     &        1/4.D0*((2 - 6*c2tb)*c2tt*MSt2diff + 
     &            (-2 + 6*c2tb)*MSt2sum + (4 - 12*c2tb)*pMT2 + 
     &            (8 - 12*c2tb)*MUEr2 - 
     &            (4*(1 + c2tb2)*MSb22**2*MUEr2)/MSb12**2 - 
     &            4*c2tb2*(6*MSb22 + MUEr2) + 
     &            4*MSb22*(-1 + Nc)*s2tb2 + 
     &            MSb12*(36*c2tb + 20*c2tb2 - 
     &               2*(18 + s2tb2 + 3*Nc*s2tb2)) + 
     &            (Abr**2*(1 + c2tb2)*(-MA02 + MSb12) + 
     &               2*MSb12*
     &                (-((1 + c2tb2)*MUEr2) + 
     &                  MSb12*(-2 + s2tb2 + Nc*s2tb2)))/MSb22 - 
     &            24*(Abr**2 + MA02)*stb2)/MSb2diff + 
     &        (((2*MSb12 - MSt22)*MSt22**2*MUEr2*stb2*stt2)/
     &            (MSb12 - MSt22)**2 + 
     &           1/4.D0*(-(Abr**2*MA02*s2tb2) + 
     &              2*(-(c2tt*MSt2diff) + MSt2sum)*MUEr2*stb2))/
     &         MSb12**2))

        tmp15 = tmp14/MSb2diff**2 + tmp15 + 
     &   (MSb12 - pMT2 - MUEr2)*s2tb2*
     &    (ctb2*ht2*(MSb12 - pMT2 - MUEr2) - 2*hb*ht*pMT*MUEr*s2tb + 
     &      hb2*(MSb12 - pMT2 - MUEr2)*stb2)*x23 + 
     &   (s2tb2*((2*ctt2*hb2*MSt12**2*MUEr2*stb2)/MSb12**2 + 
     &        1/4.D0*(2*hb*ht*MSb2diff**2*
     &             (2*MSb12 - c2tt*MSt2diff - MSt2sum)*pMT*MUEr*s2tb-
     &              2*Atr*MSb2diff**2*MSt2diff*
     &             (2*ctb2*ht2*pMT + hb*ht*MUEr*s2tb)*s2tt + 
     &            2*Atr**2*ht2*
     &             ((3 - 17*c2tb)*MSb12**3 + 
     &               MSb22*
     &                (-2*MSt12*(MSt22 + 3*c2tb*MSt22) + 
     &                  ctb2*MSb22*(c2tt*MSt2diff - MSt2sum)) + 
     &               2*MSb12*
     &                (ctt2*MSb22*MSt12 + MSt12*MSt22 + 
     &                  c2tb*
     &                   (-9*MSt12*MSt22 + 
     &                     MSb22*((3 + ctt2)*MSt12 + 3*MSt22)) + 
     &                  MSb22*MSt2sum + 
     &                  ctb2*MSb22*(MSb22 + 2*MSt22*stt2)) + 
     &               MSb12**2*
     &                (-4*(1 + 2*c2tb)*MSb22 + 
     &                  (-3 + ctt2 + c2tb*(17 + ctt2))*MSt22 + 
     &                  MSt12*(-2 - ctt2 + c2tb*(17 + stt2)))) - 
     &            ht2*(8*(-1 + 2*c2tb)*MSb12**4 + 
     &               2*MSb12**3*
     &                (2*(1 + 7*c2tb)*MSb22 + 
     &                  (-1 + 3*c2tb)*c2tt*MSt2diff + 3*MSt2sum - 
     &                  5*c2tb*MSt2sum + pMT2 - 7*c2tb*pMT2 + 
     &                  2*MUEr2 - 6*c2tb*MUEr2) + 
     &               2*MSb22*
     &                (ctb2*MSb22*
     &                   (4*MSt12*MSt22 + 
     &                     (c2tt*MSt2diff + MSt2sum)*pMT2) + 
     &                  (1 + 3*c2tb)*MSt12*MSt22*
     &                   (c2tt*MSt2diff + MSt2sum - 2*(pMT2 + MUEr2))
     &                  ) + 
     &               MSb12**2*
     &                (8*ctb2*MSb22**2 + 
     &                  4*(1 - 3*c2tb)*ctt2*MSt22**2 + 
     &                  2*c2tt*MSb22*MSt2diff + 
     &                  6*c2tb*c2tt*MSb22*MSt2diff - 
     &                  2*MSb22*MSt2sum - 22*c2tb*MSb22*MSt2sum - 
     &                  8*c2tb*MSb22*pMT2 + c2tt*MSt2diff*pMT2 + 
     &                  c2tb*c2tt*MSt2diff*pMT2 - 3*MSt2sum*pMT2 + 
     &                  13*c2tb*MSt2sum*pMT2 - 4*MSb22*MUEr2 - 
     &                  12*c2tb*MSb22*MUEr2 - 4*MSt2sum*MUEr2 + 
     &                  12*c2tb*MSt2sum*MUEr2 - 
     &                  8*MSt12*MSt22*stb2 + 4*MSt12**2*stt2 - 
     &                  12*c2tb*MSt12**2*stt2) - 
     &               2*MSb12*
     &                (2*ctb2*MSb22**2*(2*MSt2sum + pMT2) - 
     &                  (-1 + 3*c2tb)*MSt12*MSt22*
     &                   (c2tt*MSt2diff + MSt2sum - 2*(pMT2 + MUEr2))
     &                    + MSb22*
     &                   (2*
     &                      (ctt2*(MSt22**2 - MSt12*pMT2) - 
     &                       MSt2sum*MUEr2 + 
     &                       (MSt12**2 - MSt22*pMT2)*stt2) + 
     &                     c2tb*
     &                      (6*ctt2*MSt22**2 - 
     &                       2*MSt12*(4*MSt22 + (2 + ctt2)*pMT2) - 
     &                       6*MSt2sum*MUEr2 + 6*MSt12**2*stt2 - 
     &                       2*MSt22*pMT2*(2 + stt2))))))/
     &          (MSb2diff**2*(MSb12 - MSt22))))/(MSb12 - MSt12)

        tmp16 = 1/8.D0*(((2*hb*ht*MSb2diff**2*MUEr*s2tb*
     &            ((MSb22 - MSt22)*MSt22*
     &               (-((Abr + 2*Atr)*s2tt) + 2*pMT*(1 + stt2)) + 
     &              MSb12*(-2*ctt2*MSb22*pMT + Abr*MSb22*s2tt + 
     &                 MSt22*
     &                  (-((Abr + 2*Atr)*s2tt) + 2*pMT*(1 + stt2))))-
     &             ht2*(-2*
     &               (ctb2*MSb22**2*(MSb22 - MSt22)*MSt22*
     &                  (2*ctt2*(2*(Atr**2 + MSt22) + MUEr2) - 
     &                    4*Atr*pMT*s2tt + 4*pMT2*stt2) + 
     &                 MSb12**3*stb2*
     &                  (-4*Atr**2*ctt2*MSt22 + 
     &                    2*ctt2*(MSb22 - MSt22)*
     &                     (2*MSt22 + MUEr2) + 4*Atr*MSt22*pMT*s2tt - 
     &                    4*MSt22*pMT2*stt2)) + 
     &              MSb12**2*
     &               (-4*c2tb*ctt2*MSb22*
     &                  (11*MSt22**2 + MSb22*(-10*MSt22 + MUEr2)) + 
     &                 4*MSt22*
     &                  (ctt2*
     &                     (MSb22*(MSt22 + ctb2*MUEr2) - 
     &                       (2*Atr**2*(MSb22 + MSt22) + 
     &                       MSt22*(2*MSt22 + MUEr2))*stb2) + 
     &                    2*(MSb22 + MSt22)*stb2*
     &                     (Atr*pMT*s2tt - pMT2*stt2))) + 
     &              MSb12*MSb22*
     &               (4*ctb2*ctt2*MSb22**2*(2*MSt22 + MUEr2) + 
     &                 2*c2tb*MSt22**2*
     &                  (-2*ctt2*(2*Atr**2 - 10*MSt22 + MUEr2) + 
     &                    4*Atr*pMT*s2tt - 4*pMT2*stt2) - 
     &                 4*MSb22*MSt22*
     &                  (ctt2*
     &                     (MSt22 + 11*c2tb*MSt22 + MUEr2*stb2) - 
     &                    2*ctb2*
     &                     (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))))+
     &             hb2*(2*ctb2*MSb12**3*
     &               (2*ctt2*(MSb22 - MSt22)*pMT2 - 
     &                 4*MSt22*(-MSb22 + MSt22 + MUEr2)*stt2) + 
     &              MSb12*MSb22*
     &               (-4*ctt2*pMT2*
     &                  (-(ctb2*MSb22*MSt22) + c2tb*MSt22**2 + 
     &                    MSb22**2*stb2) + 
     &                 4*MSt22*
     &                  (MSt22*
     &                     (MSb22 - 11*c2tb*MSb22 + 
     &                       10*c2tb*MSt22 - 2*c2tb*MUEr2) - 
     &                    2*MSb22*(MSb22 + MUEr2)*stb2)*stt2) + 
     &              (MSb22 - MSt22)*
     &               (2*Abr*MSb2diff*(-(c2tb*MSb2diff) + MSb2sum)*
     &                  (MSb12 - MSt22)*(pMT*s2tt - Abr*stt2) + 
     &                 4*MSb22**2*MSt22*stb2*
     &                  (ctt2*pMT2 + 2*(MSt22 + MUEr2)*stt2)) + 
     &              MSb12**2*
     &               (-4*c2tb*MSb22*
     &                  (ctt2*MSb22*pMT2 + 
     &                    MSt22*(-10*MSb22 + 11*MSt22)*stt2) + 
     &                 4*MSt22*
     &                  (-(MSb22*(ctt2*pMT2*stb2 + MSt22*stt2)) + 
     &                    ctb2*
     &                     (ctt2*MSt22*pMT2 + 
     &                       2*(MSt22**2 + (MSb22 + MSt22)*MUEr2)*
     &                       stt2)))))*(-x10 + x19))/
     &       ((MSb12 - MSt22)*(-MSb22 + MSt22)) + 
     &      4*(2*hb*ht*MSb2diff**2*pMT*MUEr*s2tb + 
     &         hb2*(-(MSb12*MSb22*
     &               (MSb22 + 9*c2tb*pMT2 - c2tb*MUEr2)) + 
     &            MSb12**2*(MSb22 - ctb2*(3*pMT2 + MUEr2)) + 
     &            MSb22**2*(3*pMT2 + MUEr2)*stb2) + 
     &         ht2*(ctb2*MSb22**2*(3*pMT2 + MUEr2) - 
     &            MSb12*MSb22*(MSb22 + c2tb*(-9*pMT2 + MUEr2)) + 
     &            MSb12**2*(MSb22 - (3*pMT2 + MUEr2)*stb2)))*
     &       (x10 - x20))/MSb2diff

        tmp16 = tmp16 + (-x10 + x24)*
     &    (-(hb*ht*MSb2diff*pMT*MUEr*s2tb) + 
     &      1/4.D0*(hb2*(MSb12**2*
     &              (-2*MSb22 + pMT2 + c2tb*pMT2 + 11*MUEr2 + 
     &                3*c2tb*MUEr2) + 
     &             2*MSb12*MSb22*
     &              (MSb22 - c2tb*pMT2 + 9*c2tb*MUEr2) + 
     &             MSb22**2*(-11*MUEr2 + 3*c2tb*MUEr2 - 2*pMT2*stb2))+
     &            2*ht2*(-(ctb2*MSb22**2*(pMT2 + 3*MUEr2)) + 
     &             MSb12*MSb22*(MSb22 + c2tb*pMT2 - 9*c2tb*MUEr2) + 
     &             MSb12**2*(-MSb22 + pMT2*stb2 + 3*MUEr2*stb2)))/
     &        MSb2diff)

        tmp17 = hb2*MSb22*
     &    (-3*c2tb*(2*MA02 + 4*(Abr**2 + MSb22) - c2tt*MSt2diff + 
     &         MSt2sum - 2*(pMT2 + MUEr2)) + 
     &      16*c2tb2*MSb2diff*(1 + Nc))*s2tb2 + 
     &   (2*hb*ht*pMT*MUEr*(MSb2sum + 2*pMT2 - 2*MUEr2)*s2tb*
     &       (-2 + 3*s2tb2) + 
     &      s2tb2*(ctb2*(hb2*MSb12**2 - ht2*MSb22**2) + 
     &         hb2*MSb2diff*pMT2 + ht2*MSb2diff*pMT2 + 
     &         hb2*MSb2diff*MUEr2 + ht2*MSb2diff*MUEr2 + 
     &         c2tb*(hb2 - ht2)*
     &          (2*MSb12*MSb22 - 
     &            3*(MSb2sum + pMT2 - MUEr2)*(pMT2 + MUEr2)) + 
     &         ht2*MSb12**2*stb2 - hb2*MSb22**2*stb2))*(-x10 + x20)
     &     + 1/4.D0*((2*hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &         (-2*MSt12*(2*(1 + ctt2)*pMT + (Abr + 2*Atr)*s2tt) + 
     &           (2*MA02 - MSb2sum)*(Abr*s2tt + 2*pMT*stt2)) + 
     &        s2tb2*(ht2*(-2*MSb2diff*(2*MSt12 + MUEr2)*stt2 - 
     &              6*c2tb*
     &               (4*MSt12*(ctt2*pMT2 + Atr*pMT*s2tt) + 
     &                 (2*(2*Atr**2 + MSb2sum)*MSt12 + 
     &                    (-2*MA02 + MSb2sum + 2*MSt12)*MUEr2)*stt2)
     &              ) + hb2*
     &            ((c2tt*MSb22 - 
     &                 3*c2tb*
     &                  (2*MA02 + c2tt*MSb22 - 2*MSt12 + 
     &                    2*c2tt*MSt12))*pMT2 + 
     &              (-1 + 3*c2tb)*MSb22*(2*MSt12 + pMT2) + 
     &              12*c2tb*MSt12*MUEr2 + 
     &              2*c2tt*
     &               ((-1 + 3*c2tb)*MSb22*MSt12 + 
     &                 3*c2tb*(MA02*pMT2 + 2*MSt12*MUEr2)) + 
     &              2*(Abr*(-MSb2diff + 
     &                    3*c2tb*(-2*MA02 + MSb2sum + 2*MSt12))*
     &                  (Abr*ctt2 + pMT*s2tt) + 
     &                 (1 + 3*c2tb)*MSb12*(2*ctt2*MSt12 + pMT2*stt2)
     &                 ))))*(-x10 + x18)) + 
     &   1/2.D0*((ctt2*((12*Atr**2*c2tb*ht2*MSt22 + 
     &              hb2*(-MSb12 + MSb22 + 
     &                 c2tb*(6*MA02 - 3*(MSb12 + MSb22 + 2*MSt22)))
     &                *pMT2 + 
     &              ht2*(MSb2diff*(2*MSt22 + MUEr2) + 
     &                 3*c2tb*
     &                  (2*MSb2sum*MSt22 - 2*MA02*MUEr2 + 
     &                    MSb2sum*MUEr2 + 2*MSt22*MUEr2)))*s2tb2 - 
     &           2*hb*ht*(2*MA02 - MSb2sum)*pMT*MUEr*s2tb*
     &            (-2 + 3*s2tb2)) + 
     &        s2tb2*(-(Abr*hb2*
     &              (MSb2diff + 
     &                c2tb*(6*MA02 - 3*MSb2sum - 6*MSt22))*pMT*s2tt)
     &             + Abr**2*hb2*
     &            (MSb2diff + c2tb*(6*MA02 - 3*MSb2sum - 6*MSt22))*
     &            stt2 - 2*MSt22*
     &            (6*Atr*c2tb*ht2*pMT*s2tt + 
     &              (-6*c2tb*ht2*pMT2 + 
     &                 hb2*
     &                  (MSb12 + 3*c2tb*MSb12 - MSb22 + 
     &                    3*c2tb*MSb22 + 6*c2tb*MUEr2))*stt2)) + 
     &        hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &         (Abr*(2*MA02 - MSb2sum)*s2tt + 
     &           MSt22*(-2*(Abr + 2*Atr)*s2tt + 4*pMT*(1 + stt2))))*
     &      (x10 - x19))

        tmp17 = tmp17 + (2*hb*ht*pMT*MUEr*
     &       (MSb2sum - 2*pMT2 + 2*MUEr2)*s2tb*(-2 + 3*s2tb2) + 
     &      s2tb2*(-(ht2*(ctb2*MSb22**2 - MSb2diff*(pMT2 + MUEr2) + 
     &              c2tb*(2*MSb12*MSb22 - 
     &                 3*(MSb2sum - pMT2 + MUEr2)*(pMT2 + MUEr2)) - 
     &              MSb12**2*stb2)) + 
     &         hb2*((2 + ctb2)*MSb12**2 + 
     &            3*c2tb*(pMT2 - MUEr2)*(pMT2 + MUEr2) + 
     &            MSb22*(pMT2 - 3*c2tb*pMT2 + (5 - 3*c2tb)*MUEr2) - 
     &            MSb12*(pMT2 + 5*MUEr2 + 
     &               c2tb*(-2*MSb22 + 3*(pMT2 + MUEr2))) - 
     &            MSb22**2*(2 + stb2))))*(-x10 + x24)

        tmp18 = hb2*(1/2.D0*((1 + Nc)*s2tb2**2) + 
     &     s2tb2*(-2 - c2tb - 4*c2tb2*Nc + 
     &        (ctb2*MUEr2*((ctt2*(MSb22 - 2*MSt12)*MSt12**2)/
     &              (MSb22 - MSt12)**2 - 
     &             (MSb22*MSt22**2*stt2)/(MSb22 - MSt22)**2))/
     &         MSb22**2 + 1/4.D0*
     &         ((-4*Abr**2*c2tb2 + 2*MUEr2*(2*ctb2 + s2tb2))/MSb22 + 
     &           (Abr**2*(11 + 12*c2tb - c2tb2) + 24*ctb2*MA02 + 
     &              72*ctb2*MSb22 - 2*(1 + 3*c2tb)*c2tt*MSt2diff + 
     &              (2 + 6*c2tb)*MSt2sum - 4*(1 + 3*c2tb)*pMT2 - 
     &              6*(1 + 2*c2tb)*MUEr2 + 
     &              c2tb2*(24*MSb12 - 20*MSb22 + 6*MUEr2) + 
     &              (-4*MSb12*(-1 + Nc) + MSb22*(2 + 6*Nc))*
     &               s2tb2 + 
     &              (MSb22*
     &                 ((1 + c2tb2)*(Abr**2 + 2*MUEr2) - 
     &                   2*MSb22*(-2 + s2tb2 + Nc*s2tb2)))/MSb12)/
     &            MSb2diff) + 
     &        (-(1/4.D0*(Abr**2*MA02*s2tb2)) + 
     &           1/2.D0*(ctb2*MUEr2*
     &              ((2*ctt2*MSt12**3)/(MSb22 - MSt12)**2 - 
     &                c2tt*MSt2diff + MSt2sum + 
     &                (4*MSt22**2*stt2)/(MSb22 - MSt22))))/MSb22**2
     &        ))

        tmp18 = tmp17/MSb2diff**2 + tmp18 + 
     &   (MSb22 - pMT2 - MUEr2)*s2tb2*
     &    (ctb2*hb2*(MSb22 - pMT2 - MUEr2) + 2*hb*ht*pMT*MUEr*s2tb + 
     &      ht2*(MSb22 - pMT2 - MUEr2)*stb2)*x25 + 
     &   (s2tb2*((ctb2*hb2*MSt22**3*MUEr2*stt2)/MSb22**2 + 
     &        1/2.D0*((MSb22 - MSt22)*
     &            (hb*ht*MSb2diff**2*
     &               (-2*MSb22 + c2tt*MSt2diff + MSt2sum)*pMT*MUEr*
     &               s2tb + 
     &              Atr*MSb2diff**2*MSt2diff*s2tt*
     &               (hb*ht*MUEr*s2tb - 2*ht2*pMT*stb2) + 
     &              Atr**2*ht2*
     &               (MSb12**2*(2*MSb22 + c2tt*MSt2diff - MSt2sum)*
     &                  stb2 + 
     &                 2*MSb12*
     &                  ((-2 + 4*c2tb)*MSb22**2 + 
     &                    (-1 + 3*c2tb)*MSt12*MSt22 + 
     &                    MSb22*
     &                     (ctt2*MSt12 - 
     &                       c2tb*((3 + ctt2)*MSt12 + 3*MSt22) + 
     &                       MSt2sum + 2*MSt22*stb2*stt2)) + 
     &                 MSb22*
     &                  ((3 + 17*c2tb)*MSb22**2 + 
     &                    2*(1 + 9*c2tb)*MSt12*MSt22 - 
     &                    MSb22*
     &                     (MSt22*(3 + 17*c2tb - 2*ctt2*stb2) + 
     &                       MSt12*(3 + 17*c2tb - 2*stb2*stt2))))+
     &                ht2*
     &               (-(MSb12**2*
     &                    (4*MSb22**2 + 4*MSt12*MSt22 + 
     &                      (c2tt*MSt2diff + MSt2sum)*pMT2 - 
     &                      2*MSb22*(2*MSt2sum + pMT2))*stb2) + 
     &                 MSb22*
     &                  ((4 + 8*c2tb)*MSb22**3 + 
     &                    MSb22**2*
     &                     (c2tt*(MSt2diff + 3*c2tb*MSt2diff) - 
     &                       (3 + 5*c2tb)*MSt2sum - pMT2 - 
     &                       7*c2tb*pMT2 - 2*MUEr2 - 6*c2tb*MUEr2) + 
     &                    (1 + 3*c2tb)*MSt12*MSt22*
     &                     (c2tt*MSt2diff + MSt2sum - 
     &                       2*(pMT2 + MUEr2)) + 
     &                    2*MSb22*
     &                     (2*ctb2*MSt12*MSt22 + MSt12*pMT2 + 
     &                       3*c2tb*MSt12*pMT2 + MSt22*pMT2 + 
     &                       3*c2tb*MSt22*pMT2 + MSt2sum*MUEr2 + 
     &                       3*c2tb*MSt2sum*MUEr2 - 
     &                       ctt2*MSt22*
     &                       (MSt22 + 3*c2tb*MSt22 + pMT2*stb2) - 
     &                       MSt12**2*stt2 - 
     &                       3*c2tb*MSt12**2*stt2 - 
     &                       MSt12*pMT2*stb2*stt2)) + 
     &                 MSb12*
     &                  (2*(-1 + 7*c2tb)*MSb22**3 + 
     &                    MSb22**2*
     &                     ((-1 + 3*c2tb)*c2tt*MSt2diff + 
     &                       MSt2sum - 11*c2tb*MSt2sum - 
     &                       4*c2tb*pMT2 + 2*MUEr2 - 6*c2tb*MUEr2) + 
     &                    (-1 + 3*c2tb)*MSt12*MSt22*
     &                     (c2tt*MSt2diff + MSt2sum - 
     &                       2*(pMT2 + MUEr2)) + 
     &                    MSb22*
     &                     (c2tb*
     &                       (8*MSt12*MSt22 - 6*ctt2*MSt22**2 + 
     &                       4*MSt2sum*pMT2 + 6*MSt2sum*MUEr2 - 
     &                       6*MSt12**2*stt2) + 
     &                       2*
     &                       (-(MSt2sum*MUEr2) + 
     &                       ctt2*(MSt22**2 - 2*MSt12*pMT2*stb2) + 
     &                       MSt12**2*stt2 - 2*MSt22*pMT2*stb2*stt2)
     &                       )))))/(MSb2diff**2*(MSb22 - MSt12))))/
     &    (MSb22 - MSt22)**2

        tmp19 = ((3*Abr**2*c2tb*ctt2*hb2*s2tb2 + 
     &        Abr*(hb*ht*MUEr*s2tb*(2 - 3*s2tb2) + 
     &           3*c2tb*hb2*pMT*s2tb2)*s2tt + 
     &        (2*hb*ht*pMT*MUEr*s2tb*(2 - 3*s2tb2) + 
     &           3*c2tb*(hb2*pMT2 - ht2*MUEr2)*s2tb2)*stt2)*
     &      (-(x13*x28) + x11*x30))/(MSt12*x11*x13) + 
     &   ((ctt2*(3*c2tb*(-(hb2*pMT2) + ht2*MUEr2)*s2tb2 + 
     &           2*hb*ht*pMT*MUEr*s2tb*(-2 + 3*s2tb2)) + 
     &        Abr*(hb*ht*MUEr*s2tb*(2 - 3*s2tb2)*s2tt + 
     &           3*c2tb*hb2*s2tb2*(pMT*s2tt - Abr*stt2)))*
     &      (x14*x29 - x12*x31))/(MSt22*x12*x14) + 
     &   (hb*ht*(MSb2sum - 2*MSt12)*MUEr*s2tb*(-2 + 3*s2tb2)*
     &       (2*ctt2*pMT + Atr*s2tt) + 
     &      s2tb2*(ctt2*hb2*
     &          (MSb2diff - 3*c2tb*MSb2sum + 6*c2tb*MSt12)*MUEr2 + 
     &         ht2*(MSb2diff + 3*c2tb*(MSb2sum - 2*MSt12))*
     &          (ctt2*pMT2 + Atr*(pMT*s2tt + Atr*stt2))))*(x33 - x34)+
     &     (x10 - x22)**2*
     &    (-(hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &         ((MSb2sum + c2tt*MSt2diff - MSt2sum)*pMT + 
     &           Atr*MSt2diff*s2tt)) + 
     &      1/2.D0*(s2tb2*(-(hb2*
     &              ((5*c2tb2*MSb22 + 4*MSb2diff - 
     &                   3*c2tb*(MSb2sum + c2tt*MSt2diff - MSt2sum)
     &                   )*MUEr2 + 
     &                c2tb2*
     &                 (-8*Abr**2*MA02 + 11*MSb12*MUEr2 + 
     &                   8*MSb22*(MSb12 + MSb22)*(1 + Nc)) - 
     &                (MSb12**2 - MSb22**2)*(1 + Nc)*s2tb2 + 
     &                Abr**2*
     &                 (4*c2tb2*(2*MSb22 + MSb2sum) + 
     &                   MSb2diff*s2tb2))) + 
     &           ht2*(-(Atr**2*
     &                 (MSb2diff - 
     &                   3*c2tb*
     &                    (-MSb2sum + c2tt*MSt2diff + MSt2sum))) - 
     &              (MSb2diff + 
     &                 3*c2tb*(MSb2sum + c2tt*MSt2diff - MSt2sum))*
     &               pMT2 - 6*Atr*c2tb*MSt2diff*pMT*s2tt)))) + 
     &   (x10 - x21)**2*(hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &       ((MSb2sum + c2tt*MSt2diff - MSt2sum)*pMT + 
     &         Atr*MSt2diff*s2tt) + 
     &      1/2.D0*(s2tb2*(-(hb2*
     &              (c2tb2*
     &                 (-8*Abr**2*MA02 + 
     &                   4*MSb12*
     &                    (3*MUEr2 + 2*MSb12*(1 + Nc) + 
     &                      2*MSb22*(1 + Nc))) + 
     &                (MSb12**2 - MSb22**2)*(1 + Nc)*s2tb2 + 
     &                Abr**2*
     &                 (4*c2tb2*(2*MSb12 + MSb2sum) - 
     &                   MSb2diff*s2tb2) + 
     &                MUEr2*
     &                 (4*c2tb2*MSb22 + 
     &                   3*c2tb*
     &                    (MSb2sum + c2tt*MSt2diff - MSt2sum) - 
     &                   MSb2diff*(1 + s2tb2)))) + 
     &           ht2*((MSb2diff + 
     &                 3*c2tb*(MSb2sum + c2tt*MSt2diff - MSt2sum))*
     &               pMT2 + 6*Atr*c2tb*MSt2diff*pMT*s2tt + 
     &              Atr**2*(MSb2diff + 
     &                 3*c2tb*
     &                  (MSb2sum - 2*ctt2*MSt22 - 2*MSt12*stt2)))))
     &      )

        tmp19 = tmp19 + (hb*ht*(MSb2sum - 2*MSt22)*MUEr*s2tb*
     &       (-2 + 3*s2tb2)*(-(Atr*s2tt) + 2*pMT*stt2) + 
     &      s2tb2*(hb2*(MSb2diff - 3*c2tb*MSb2sum + 6*c2tb*MSt22)*
     &          MUEr2*stt2 + 
     &         ht2*(MSb2diff + 3*c2tb*(MSb2sum - 2*MSt22))*
     &          (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2)))*x35 + 
     &   (-(hb2*(MSb2diff - 3*c2tb*MSb2sum + 6*c2tb*MSt22)*MUEr2*
     &         s2tb2*stt2) + 
     &      hb*ht*(MSb2sum - 2*MSt22)*MUEr*s2tb*(-2 + 3*s2tb2)*
     &       (Atr*s2tt - 2*pMT*stt2) - 
     &      ht2*(MSb2diff + 3*c2tb*(MSb2sum - 2*MSt22))*s2tb2*
     &       (Atr**2*ctt2 - Atr*pMT*s2tt + pMT2*stt2))*x36 + 
     &   1/2.D0*((MUEr*(-((hb2 + ht2)*MSb2diff) + 
     &            c2tb*(hb2 - ht2)*(3*MSb2sum - 6*(pMT2 + MUEr2)))*
     &          s2tb2 + 8*hb*ht*pMT*MUEr2*s2tb*(-2 + 3*s2tb2))*x37)/
     &     (MUEr*MUEr2*x25) + 
     &   1/8.D0*(Abr**2*hb2*s2tb2*
     &      (((8*c2tb2*MSb12*(-4*MSb12 + MSb2diff) + 
     &             2*MSb2diff*(-MSb12 + MSb2diff)*s2tb2)*x45)/
     &         (MSb12**3*x44) + 
     &        (2*(2*MSb22*(2*c2tb2*MSb12 + MSb22 - 
     &                11*c2tb2*MSb22) - 
     &             MSb12*(2*MSb22 + MSb2diff)*s2tb2)*x47)/
     &         (MSb22**3*x46)))

        tmp20 = x32*((hb*ht*pMT*s2tb*
     &         (4*c2tb2*MSb12 - MSb2sum*s2tb2 + 
     &           2*(pMT2 + MUEr2)*(-2 + 3*s2tb2)))/MUEr + 
     &      1/4.D0*((MSb12 - pMT2 - MUEr2)*s2tb2*
     &          (ht2*((-1 + 3*c2tb)*MSb12 + 6*ctb2*MSb22 - 
     &               2*(1 + 3*c2tb)*(pMT2 + MUEr2)) - 
     &            hb2*(MSb12 + 3*c2tb*MSb12 - 
     &               2*(-1 + 3*c2tb)*(pMT2 + MUEr2) - 6*MSb22*stb2)))
     &         /MUEr2) + (x28*
     &      (hb*ht*(MA02 - MSb12 + MSt12)*MUEr*s2tb*
     &         (Abr*s2tt + 2*pMT*stt2) + 
     &        1/4.D0*(s2tb2*(-(ht2*
     &                ((2 + 6*c2tb)*MA02 - MSb2sum - 
     &                  c2tb*(4*MSb12 + MSb2sum - 6*MSt12) + 
     &                  2*MSt12)*MUEr2*stt2) - 
     &             hb*ht*(6*MA02 - 4*MSb12 - MSb2sum + 6*MSt12)*
     &              MUEr*s2tb*(Abr*s2tt + 2*pMT*stt2) + 
     &             hb2*((-2 + 6*c2tb)*MA02 + MSb2sum - 
     &                c2tb*(4*MSb12 + MSb2sum - 6*MSt12) - 2*MSt12)
     &               *(Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2)))))/
     &    MSt12 + (x29*(hb*ht*(MA02 - MSb12 + MSt22)*MUEr*s2tb*
     &         (2*ctt2*pMT - Abr*s2tt) + 
     &        1/4.D0*(s2tb2*(-(ctt2*ht2*
     &                ((2 + 6*c2tb)*MA02 - MSb2sum - 
     &                  c2tb*(4*MSb12 + MSb2sum - 6*MSt22) + 
     &                  2*MSt22)*MUEr2) - 
     &             hb*ht*(6*MA02 - 4*MSb12 - MSb2sum + 6*MSt22)*
     &              MUEr*s2tb*(2*ctt2*pMT - Abr*s2tt) + 
     &             hb2*((-2 + 6*c2tb)*MA02 + MSb2sum - 
     &                c2tb*(4*MSb12 + MSb2sum - 6*MSt22) - 2*MSt22)
     &               *(ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))))))/
     &    MSt22 + 1/4.D0*((-4*hb*ht*(MSb22 - MSt12)*MUEr*s2tb*
     &          (Abr*s2tt + 2*pMT*stt2) + 
     &         hb*ht*MA02*MUEr*s2tb*(4 - 6*s2tb2)*
     &          (Abr*s2tt + 2*pMT*stt2) + 
     &         s2tb2*(ht2*MUEr2*
     &             (4*c2tb*MSb22 + (2 - 6*c2tb)*(MA02 + MSt12) - 
     &               2*MSb2sum*stb2)*stt2 + 
     &            hb*ht*(4*MSb22 + MSb2sum - 6*MSt12)*MUEr*s2tb*
     &             (Abr*s2tt + 2*pMT*stt2) - 
     &            hb2*(MSb2sum + 
     &               c2tb*(4*MSb22 + MSb2sum - 6*MSt12) - 
     &               2*(MA02 + 3*c2tb*MA02 + MSt12))*
     &             (Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2)))*x30)/
     &     MSt12 + 1/4.D0*((-4*hb*ht*(MSb22 - MSt22)*MUEr*s2tb*
     &          (2*ctt2*pMT - Abr*s2tt) + 
     &         hb*ht*MUEr*s2tb*
     &          (MA02*(4 - 6*s2tb2) + 
     &            (4*MSb22 + MSb2sum - 6*MSt22)*s2tb2)*
     &          (2*ctt2*pMT - Abr*s2tt) + 
     &         s2tb2*(ctt2*ht2*MUEr2*
     &             (4*c2tb*MSb22 + (2 - 6*c2tb)*(MA02 + MSt22) - 
     &               2*MSb2sum*stb2) - 
     &            hb2*(MSb2sum + 
     &               c2tb*(4*MSb22 + MSb2sum - 6*MSt22) - 
     &               2*(MA02 + 3*c2tb*MA02 + MSt22))*
     &             (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))))*x31)/
     &     MSt22

        tmp20 = tmp20 + x37*
     &    ((hb*ht*pMT*s2tb*(MSb22*(4 - 5*s2tb2) - MSb12*s2tb2 + 
     &           2*(pMT2 + MUEr2)*(-2 + 3*s2tb2)))/MUEr + 
     &      1/4.D0*((MSb22 - pMT2 - MUEr2)*s2tb2*
     &          (hb2*(-6*ctb2*MSb12 + MSb22 - 3*c2tb*MSb22 + 
     &               2*(1 + 3*c2tb)*(pMT2 + MUEr2)) + 
     &            ht2*(MSb22 + 3*c2tb*MSb22 - 
     &               2*(-1 + 3*c2tb)*(pMT2 + MUEr2) - 6*MSb12*stb2)))
     &         /MUEr2) + s2tb2*
     &    (-(1/2.D0*(ht2*(ctb2*MSb22**2*
     &               (c2tt*MSt2diff + MSt2sum - 2*(pMT2 + MUEr2)) + 
     &              MSb12*MSb22*
     &               (6*(1 + 2*c2tb)*MSb22 + 
     &                 5*c2tb*
     &                  (c2tt*MSt2diff + MSt2sum - 2*(pMT2 + MUEr2)))
     &                + 2*Atr**2*
     &               (11*c2tb*MSb12*MSb22 + ctb2*MSb22**2 - 
     &                 MSb12**2*stb2) + 
     &              MSb12**2*
     &               (6*(-1 + 2*c2tb)*MSb22 - 
     &                 (c2tt*MSt2diff + MSt2sum - 2*(pMT2 + MUEr2))*
     &                  stb2)))/(MSb12*MSb22)) + 
     &      hb2*(c2tb*(11*Abr**2 + 5*MA02 + 6*MSb2sum - 
     &            5*(pMT2 + MUEr2)) + (-2 + 6*c2tb2)*MUEr2*x50 - 
     &         2*(MSb2sum - 2*MUEr2)*(x51 - x52) + 
     &         1/4.D0*(Abr**2*
     &             (-((5 + 13*c2tb2)*MA02**2*MSb22) - 
     &               3*(1 + c2tb2)*MSb12**2*(2*MA02 + MSb2diff) + 
     &               MSb22**2*
     &                (8*(1 + 3*c2tb2)*MA02 - 
     &                  (7 + 31*c2tb2)*MSb2diff) + 
     &               MSb12*
     &                (MA02*
     &                   (3*(1 + c2tb2)*MA02 + 
     &                     (6 + 22*c2tb2)*MSb22) + 
     &                  2*(5 + 9*c2tb2)*MSb22*MSb2diff))*x49)/
     &           MSb22**3))

        tmp21 = -(1/8.D0*(s2tb2*
     &        (4*Abr**2*hb2*stb2*
     &           (ctt2*(-MA02**2 + (MSb12 - MSt12)**2)*x11 + 
     &             (-MA02**2 + (MSb12 - MSt22)**2)*stt2*x12) - 
     &          2*Abr*s2tt*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2)*
     &           (-(MSt12**2*x11) + MA02**2*(x11 - x12) + 
     &             MSt22**2*x12 + MSb12**2*(-x11 + x12) + 
     &             2*MSb12*(MSt12*x11 - MSt22*x12)) + 
     &          2*hb2*pMT2*stb2*
     &           (2*MSt12**2*stt2*x11 + MSt22**2*x12 + 
     &             c2tt*MSt22**2*x12 - 
     &             2*MA02**2*(stt2*x11 + ctt2*x12) + 
     &             2*MSb12**2*(stt2*x11 + ctt2*x12) - 
     &             4*MSb12*(MSt12*stt2*x11 + ctt2*MSt22*x12)) + 
     &          4*(ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb)*
     &           (MSt12**2*stt2*x11 + ctt2*MSt22**2*x12 - 
     &             MA02**2*(stt2*x11 + ctt2*x12) + 
     &             MSb12**2*(stt2*x11 + ctt2*x12) - 
     &             2*MSb12*(MSt12*stt2*x11 + ctt2*MSt22*x12))))/
     &      MSb12)

        tmp21 = tmp21 + 1/8.D0*
     &    (s2tb2*(-2*Abr**2*ctb2*hb2*
     &          (2*ctt2*(-MA02**2 + (MSb22 - MSt12)**2)*x13 + 
     &            2*(-MA02**2 + (MSb22 - MSt22)**2)*stt2*x14) + 
     &         2*Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt*
     &          (-(MSt12**2*x13) + MA02**2*(x13 - x14) + 
     &            MSt22**2*x14 + MSb22**2*(-x13 + x14) + 
     &            2*MSb22*(MSt12*x13 - MSt22*x14)) + 
     &         4*(-(ctb2*hb2*pMT2) + hb*ht*pMT*MUEr*s2tb - 
     &            ht2*MUEr2*stb2)*
     &          (MSt12**2*stt2*x13 + ctt2*MSt22**2*x14 - 
     &            MA02**2*(stt2*x13 + ctt2*x14) + 
     &            MSb22**2*(stt2*x13 + ctt2*x14) - 
     &            2*MSb22*(MSt12*stt2*x13 + ctt2*MSt22*x14))))/
     &     MSb22

        tmp21 = tmp21 + ((-x10 + x21)*
     &       (hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &          ((2*MA02 + MSb2sum - c2tt*MSt2diff - MSt2sum)*pMT + 
     &            Abr*MSt2diff*s2tt) + 
     &         1/2.D0*(s2tb2*
     &            (hb2*((MSb2diff + 
     &                    3*c2tb*
     &                     (-MSb2sum + c2tt*MSt2diff + MSt2sum))*
     &                  pMT2 - 
     &                 3*Abr**2*
     &                  (2*ctb2*MSb12 + 2*c2tb2*MSb2diff + 
     &                    c2tb*c2tt*MSt2diff - c2tb*MSt2sum - 
     &                    2*MSb22*stb2) - 
     &                 6*(Abr*c2tb*MSt2diff*pMT*s2tt + 
     &                    MA02*
     &                     (2*ctb2*MSb12 + c2tb*(Abr**2 + pMT2) - 
     &                       2*MSb22*stb2))) + 
     &              ht2*MUEr2*
     &               (MSb2diff + 
     &                 3*c2tb*
     &                  (2*MA02 + MSb2sum - 2*ctt2*MSt22 - 
     &                    2*MSt12*stt2))))) + 
     &      (-x10 + x22)*(-(hb*ht*MUEr*s2tb*(-2 + 3*s2tb2)*
     &            ((2*MA02 + MSb2sum - c2tt*MSt2diff - MSt2sum)*
     &               pMT + Abr*MSt2diff*s2tt)) + 
     &         1/2.D0*(s2tb2*
     &            (ht2*(-MSb2diff + 
     &                 3*c2tb*
     &                  (-2*MA02 - MSb2sum + c2tt*MSt2diff + 
     &                    MSt2sum))*MUEr2 + 
     &              hb2*(3*Abr**2*
     &                  (2*ctb2*MSb12 + 2*c2tb2*MSb2diff + 
     &                    c2tb*c2tt*MSt2diff - c2tb*MSt2sum - 
     &                    2*MSb22*stb2) + 
     &                 6*(Abr*c2tb*MSt2diff*pMT*s2tt + 
     &                    MA02*
     &                     (2*ctb2*MSb12 + c2tb*(Abr**2 + pMT2) - 
     &                       2*MSb22*stb2)) - 
     &                 pMT2*
     &                  (MSb2diff - 3*c2tb*MSb2sum + 
     &                    6*c2tb*ctt2*MSt22 + 6*c2tb*MSt12*stt2))))
     &         ))/MSb2diff**2 - 
     &   1/4.D0*(Abr**2*hb2*s2tb2*
     &       ((1 + c2tb2)*MSb12**3 - 
     &         6*MA02**2*MSb12*MSb22**2*s2tb2*x44 + 
     &         MA02*MSb22**2*s2tb2*(-1 + MA02**2*x44) + 
     &         MA02**2*MSb12**2*(MA02 - 6*MSb22)*s2tb2*x46))/
     &     (MSb12**2*MSb22**2) + 
     &   s2tb2*(1/2.D0*(hb2*(Abr**2*MSb22 + 
     &             MA02*((3 + c2tb)*MSb12 - (-3 + c2tb)*MSb22)) + 
     &          ctb2*ht2*MSb22*MUEr2 + hb*ht*MSb2diff*pMT*MUEr*s2tb + 
     &          ht2*MSb12*MUEr2*stb2)/(MSb12*MSb22) + 
     &      hb2*((-6*c2tb*MA02)/MSb2diff + 
     &         1/4.D0*(Abr**2*
     &              ((5 + c2tb)*MSb12 + (-c2tb + c2tb2)*MSb22) + 
     &             2*pMT2*(ctb2*MSb12 + MSb22*stb2))/(MSb12*MSb22)+
     &           (Abr**2*(c2tb2*MSb22 + 1/2.D0*MA02 - 
     &              1/4.D0*((1 + c2tb2)*
     &                  (MA02**2*
     &                     ((MA02 - 3*MSb12)*MSb12 - MSb22**2) + 
     &                    (3*MA02*MSb12 - MSb12**2 + 
     &                       4*MSb12*MSb22 + MSb22**2)*MSb2diff**2)
     &                   *x48)/MSb12))/MSb22**2))

        s2tb2F3b = tmp19/MSb2diff**2 + tmp20/MSb2diff + 
     &   tmp15*(-x10 + x21) + tmp18*(-x10 + x22) + 
     &   (s2tb2*tmp16 + 1/8.D0*
     &       (s2tb2*(-2*hb*ht*MSb2diff**2*MUEr*s2tb*
     &             ((MSb22 - MSt12)*
     &                (Abr*(MSb12 - MSt12)*s2tt - 
     &                  2*MSt12*(pMT + ctt2*pMT + Atr*s2tt)) + 
     &               MSb12*
     &                (-2*MSt12*(pMT + ctt2*pMT + Atr*s2tt) + 
     &                  2*MSb22*pMT*stt2)) - 
     &            hb2*(2*Abr**2*ctt2*MSb2diff*
     &                (-(c2tb*MSb2diff) + MSb2sum)*(MSb12 - MSt12)*
     &                (MSb22 - MSt12) + 
     &               2*Abr*MSb2diff*(-(c2tb*MSb2diff) + MSb2sum)*
     &                (MSb12 - MSt12)*(MSb22 - MSt12)*pMT*s2tt - 
     &               2*MSb22**2*(MSb22 - MSt12)*MSt12*stb2*
     &                (4*ctt2*(MSt12 + MUEr2) + 2*pMT2*stt2) + 
     &               2*ctb2*MSb12**3*
     &                (4*ctt2*MSt12*(-MSb22 + MSt12 + MUEr2) + 
     &                  2*(-MSb22 + MSt12)*pMT2*stt2) - 
     &               MSb12*MSb22*
     &                (4*ctt2*MSt12*
     &                   (MSt12*
     &                      (MSb22 - 11*c2tb*MSb22 + 
     &                       10*c2tb*MSt12 - 2*c2tb*MUEr2) - 
     &                     2*MSb22*(MSb22 + MUEr2)*stb2) - 
     &                  4*pMT2*
     &                   (-(ctb2*MSb22*MSt12) + c2tb*MSt12**2 + 
     &                     MSb22**2*stb2)*stt2) - 
     &               MSb12**2*
     &                (4*ctt2*MSt12*
     &                   (c2tb*MSb22*(10*MSb22 - 11*MSt12) + 
     &                     2*ctb2*MSt12*(MSt12 + MUEr2) - 
     &                     MSb22*(MSt12 - 2*ctb2*MUEr2)) - 
     &                  4*pMT2*
     &                   (c2tb*MSb22**2 + 
     &                     MSt12*(-(ctb2*MSt12) + MSb22*stb2))*stt2
     &                  )) - 
     &            ht2*(2*MSb12**3*stb2*
     &                (4*ctt2*MSt12*pMT2 + 4*Atr*MSt12*pMT*s2tt + 
     &                  4*Atr**2*MSt12*stt2 - 
     &                  2*(MSb22 - MSt12)*(2*MSt12 + MUEr2)*stt2) - 
     &               4*ctb2*MSb22**2*(MSb22 - MSt12)*MSt12*
     &                (2*ctt2*pMT2 + 2*Atr*pMT*s2tt + 
     &                  (2*(Atr**2 + MSt12) + MUEr2)*stt2) + 
     &               MSb12**2*
     &                (-8*MSt12*(MSb22 + MSt12)*
     &                   (ctt2*pMT2 + Atr*pMT*s2tt)*stb2 - 
     &                  4*(-(MSb22*MSt12*(MSt12 + ctb2*MUEr2)) + 
     &                     c2tb*MSb22*
     &                      (11*MSt12**2 + 
     &                       MSb22*(-10*MSt12 + MUEr2)) + 
     &                     MSt12*
     &                      (2*Atr**2*(MSb22 + MSt12) + 
     &                       MSt12*(2*MSt12 + MUEr2))*stb2)*stt2) + 
     &               MSb12*MSb22*
     &                (4*ctb2*MSb22*
     &                   (2*ctt2*MSt12*pMT2 + 2*Atr*MSt12*pMT*s2tt + 
     &                     2*Atr**2*MSt12*stt2 + 
     &                     MSb22*(2*MSt12 + MUEr2)*stt2) - 
     &                  4*MSt12*
     &                   (MSb22*(MSt12 + MUEr2*stb2)*stt2 + 
     &                     c2tb*MSt12*
     &                      (2*ctt2*pMT2 + 2*Atr*pMT*s2tt + 
     &                       (2*Atr**2 + 11*MSb22 - 10*MSt12 + 
     &                      MUEr2)*stt2)))))*(-x10 + x18))/
     &        (MSb2diff*(MSb12 - MSt12)*(-MSb22 + MSt12)))/
     &    (MSb12*MSb22)

        s2tb2F3b = s2tb2F3b + 
     &   1/4.D0*(s2tb2*(-(((MA02 - MSb12 + MSt12)*
     &              (2*Abr**2*ctt2*hb2*stb2 + 
     &                Abr*s2tt*(hb*ht*MUEr*s2tb + 2*hb2*pMT*stb2) + 
     &                2*(ctb2*ht2*MUEr2 + hb*ht*pMT*MUEr*s2tb + 
     &                   hb2*pMT2*stb2)*stt2)*x11*
     &              (MSt12**2*(x10 - x18) - 
     &                MSb12*MSt12*(x10 + x18 - 2*x21 - x28) - 
     &                MSb12**2*x28 + 
     &                MA02*(MSt12*(-x10 + x18) + MSb12*x28)))/MSt12
     &            ) + ((MA02 - MSb12 + MSt22)*
     &            (2*ctb2*ctt2*ht2*MUEr2 + 
     &              2*ctt2*(hb*ht*pMT*MUEr*s2tb + hb2*pMT2*stb2) + 
     &              Abr*(-(hb*ht*MUEr*s2tb*s2tt) - 
     &                 2*hb2*pMT*s2tt*stb2 + 2*Abr*hb2*stb2*stt2))*
     &            x12*(MSt22**2*(-x10 + x19) + 
     &              MSb12*MSt22*(x10 + x19 - 2*x21 - x29) + 
     &              MSb12**2*x29 + 
     &              MA02*(MSt22*(x10 - x19) - MSb12*x29)))/MSt22))/
     &     MSb12

        s2tb2F3b = s2tb2F3b - 
     &   1/4.D0*((MA02 - MSb22 + MSt12)*s2tb2*
     &       (2*Abr**2*ctb2*ctt2*hb2 + 
     &         Abr*(2*ctb2*hb2*pMT - hb*ht*MUEr*s2tb)*s2tt + 
     &         2*(ctb2*hb2*pMT2 - hb*ht*pMT*MUEr*s2tb + 
     &            ht2*MUEr2*stb2)*stt2)*x13*
     &       (MSt12**2*(x10 - x18) - 
     &         MSb22*MSt12*(x10 + x18 - 2*x22 - x30) - 
     &         MSb22**2*x30 + MA02*(MSt12*(-x10 + x18) + MSb22*x30)
     &         ))/(MSb22*MSt12)

        s2tb2F3b = s2tb2F3b - 
     &   1/4.D0*((MA02 - MSb22 + MSt22)*s2tb2*
     &       (-2*ctt2*hb*ht*pMT*MUEr*s2tb + Abr*hb*ht*MUEr*s2tb*s2tt + 
     &         2*ctt2*ht2*MUEr2*stb2 + 
     &         2*ctb2*hb2*(ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*
     &       x14*(MSt22**2*(x10 - x19) - 
     &         MSb22*MSt22*(x10 + x19 - 2*x22 - x31) - 
     &         MSb22**2*x31 + MA02*(MSt22*(-x10 + x19) + MSb22*x31)
     &         ))/(MSb22*MSt22)

        s2tb2F3b = s2tb2F3b + 
     &   1/2.D0*((hb2*s2tb2*((MSb12*
     &              (ctb2*(-(c2tt*MSt2diff) + MSt2sum - 2*pMT2) + 
     &                (3 + c2tb)*(Abr**2 + MA02 - MUEr2) - 
     &                MSb12*(-1 + Nc)*s2tb2))/MSb22 + 
     &           (MSb22*((-3 + c2tb)*(Abr**2 + MA02 - MUEr2) + 
     &                MSb22*(-1 + Nc)*s2tb2 + 
     &                (c2tt*MSt2diff - MSt2sum + 2*pMT2)*stb2))/
     &            MSb12))/MSb2diff + 
     &      x23*(((MSb12 - pMT2 - MUEr2)*s2tb2*
     &            (-2*hb*ht*pMT*MUEr*s2tb + 
     &              (MSb12 - pMT2 - MUEr2)*(ctb2*ht2 + hb2*stb2))*
     &            (MSb12*(2*x10 - x20 - x24) - 
     &              (pMT2 - MUEr2)*(x20 - x24)))/MSb12 + 
     &         (-MSb12 + pMT2 + MUEr2)**2*s2tb2*
     &          ((2*hb*ht*pMT*s2tb)/MUEr + 
     &            ((-MSb12 + pMT2 + MUEr2)*(ctb2*ht2 + hb2*stb2))/
     &             MUEr2)*x32) + 
     &      x25*(((MSb22 - pMT2 - MUEr2)*s2tb2*
     &            (2*hb*ht*pMT*MUEr*s2tb + 
     &              (MSb22 - pMT2 - MUEr2)*(ctb2*hb2 + ht2*stb2))*
     &            (MSb22*(2*x10 - x20 - x24) - 
     &              (pMT2 - MUEr2)*(x20 - x24)))/MSb22 + 
     &         (-MSb22 + pMT2 + MUEr2)**2*s2tb2*
     &          ((-2*hb*ht*pMT*s2tb)/MUEr + 
     &            ((-MSb22 + pMT2 + MUEr2)*(ctb2*hb2 + ht2*stb2))/
     &             MUEr2)*x37) + 
     &      (((MUEr*((hb2 + ht2)*MSb2diff - 
     &                 3*c2tb*(hb2 - ht2)*
     &                  (MSb2sum - 2*(pMT2 + MUEr2)))*s2tb2 - 
     &              8*hb*ht*pMT*MUEr2*s2tb*(-2 + 3*s2tb2))*x32)/
     &          (MUEr*MUEr2*x23) + 
     &         (Abr**2*hb2*((1 + c2tb2)*MSb12**2 - 
     &              (3 + 7*c2tb2)*MSb12*MSb22 + 
     &              (2 + 22*c2tb2)*MSb22**2)*s2tb2*x49)/
     &          (MSb22**3*x48))/MSb2diff**2)

        s2tb2F3b = s2tb2F3b + tmp21*(-x10 + x9) + 
     &   1/8.D0*(hb2*s2tb2*(-72 + 
     &        8*c2tb2*(-7 + 
     &           (MSb12**3 - MSb22**3)/(MSb12*MSb22*MSb2diff) - 
     &           8*Nc) + (20*c2tb*(-(c2tt*MSt2diff) + MSt2sum) - 
     &           (Abr**2*MA02*
     &              (8*c2tb2*(MA02 - 6*MSb12)*MSb12 + 
     &                (MA02*(-5*MSb12 + 3*MSb22) + 
     &                   12*MSb12*(MSb12 - MSb2diff))*s2tb2)*x45)/
     &            MSb12**3 + 
     &           (Abr**2*MA02*
     &              (8*c2tb2*(MA02 - 6*MSb22)*MSb22 + 
     &                (MA02*(3*MSb12 - 5*MSb22) + 
     &                   12*MSb22*(MSb22 + MSb2diff))*s2tb2)*x47)/
     &            MSb22**3)/MSb2diff + 
     &        s2tb2*(-4*(1 + 3*Nc) + 
     &           Abr**2*MA02**2*
     &            (((MA02 - 6*MSb12)*x44*
     &                 (MA02*x45 - 2*MSb12*(x10 - x21 + x45)))/
     &               MSb12**3 + 
     &              ((MA02 - 6*MSb22)*x46*
     &                 (MA02*x47 - 2*MSb22*(x10 - x22 + x47)))/
     &               MSb22**3)) + 
     &        (2*Abr**2*(1 + c2tb2)*x48*
     &           (MSb22*((MA02 - MSb12)**3 - 
     &                2*(MA02 - MSb12)*(MA02 + 2*MSb12)*MSb22 + 
     &                3*(MA02 + MSb12)*MSb22**2 - 6*MSb22**3)*
     &              (-x10 + x21) + 
     &             (MSb22**2*
     &                (MA02*(MA02 - 4*MSb12)*MSb12 + 
     &                  (-MA02**2 - 2*MSb12**2 + 9*MSb12*MSb22 + 
     &                     MSb22**2)*MSb2diff)*(-x10 + x22))/MSb12+
     &               (MA02**4 + 
     &                MA02**2*
     &                 (6*MSb12**2 - 2*MSb12*MSb22 + 4*MSb22**2) + 
     &                (MSb12**2 - 4*MSb12*MSb22 - MSb22**2)*
     &                 MSb2diff**2 - 
     &                2*MA02*
     &                 (2*MSb12**3 - 5*MSb12**2*MSb22 + 
     &                   3*MSb22**3 + MA02**2*(MSb12 + MSb2sum)))*
     &              x49))/MSb22**3))

        tmp22 = 1/8.D0*(((-2*hb*ht*MUEr*s2tb2*
     &           (Abr*MSb2sum*(MSb12 - MSt12)*(MSb22 - MSt12)*
     &              s2tt + 
     &             2*MSb12*MSb22*(MSb2sum - 2*MSt12)*
     &              (pMT + ctt2*pMT + Atr*s2tt) - 
     &             2*MSt12*(MSb12**2 + MSb22**2 - MSb2sum*MSt12)*
     &              pMT*stt2))/(-MSb22 + MSt12) + 
     &        (s2tb*(hb2*(4*Abr*(MSb12 - MSt12)*(-MSb22 + MSt12)*
     &                 (Abr*ctt2 + pMT*s2tt)*
     &                 (ctb2*MSb12 - MSb22*stb2) - 
     &                4*MSb22*(MSb22 - MSt12)*MSt12*pMT2*stb2*
     &                 stt2 - 
     &                MSb12**2*
     &                 (8*ctt2*MSb22*
     &                    (c2tb*(MSb22 - MSt12) + ctb2*MUEr2) + 
     &                   4*ctb2*(MSb22 - MSt12)*pMT2*stt2) - 
     &                MSb12*
     &                 (4*ctb2*MSt12**2*pMT2*stt2 + 
     &                   4*c2tb*MSb22*MSt12*
     &                    (-2*ctt2*(MSb22 - MSt12 + MUEr2) - 
     &                      pMT2*stt2) - 
     &                   4*MSb22**2*stb2*(2*ctt2*MUEr2 + pMT2*stt2)))
     &               + ht2*
     &              (MSb12**2*
     &                 (-8*ctt2*MSb22*pMT2*stb2 + 
     &                   4*(MSb22 - MSt12)*
     &                    (2*c2tb*MSb22 - MUEr2*stb2)*stt2) + 
     &                MSb22*
     &                 (-4*ctb2*(MSb22 - MSt12)*MSt12*MUEr2*stt2 - 
     &                   8*Atr*MSb12*
     &                    (-(ctb2*MSb22) + c2tb*MSt12 + 
     &                      MSb12*stb2)*(pMT*s2tt + Atr*stt2)) + 
     &                MSb12*
     &                 (4*ctb2*MSb22**2*(2*ctt2*pMT2 + MUEr2*stt2) - 
     &                   4*MSt12*
     &                    (MSt12*MUEr2*stb2*stt2 + 
     &                      c2tb*MSb22*
     &                       (2*ctt2*pMT2 + 
     &                       (2*MSb22 - 2*MSt12 + MUEr2)*stt2))))))/
     &         (MSb22 - MSt12))*(-x10 + x18))/(MSb12 - MSt12)

        tmp22 = tmp22 + 1/8.D0*
     &    (((s2tb*(-(hb2*(MSb12**2*
     &                    (4*ctb2*ctt2*(MSb22 - MSt22)*pMT2 + 
     &                      8*MSb22*
     &                       (c2tb*(MSb22 - MSt22) + ctb2*MUEr2)*
     &                       stt2) + 
     &                   MSb12*
     &                    (-4*ctt2*pMT2*
     &                       (c2tb*MSb22*MSt22 - ctb2*MSt22**2 + 
     &                       MSb22**2*stb2) - 
     &                      8*MSb22*
     &                       (c2tb*MSt22*(MSb22 - MSt22 + MUEr2) + 
     &                       MSb22*MUEr2*stb2)*stt2) + 
     &                   (MSb22 - MSt22)*
     &                    (4*ctt2*MSb22*MSt22*pMT2*stb2 + 
     &                      4*Abr*(MSb12 - MSt22)*
     &                       (ctb2*MSb12 - MSb22*stb2)*
     &                       (-(pMT*s2tt) + Abr*stt2)))) + 
     &              ht2*(ctt2*MSb22*
     &                  (4*ctb2*MSt22*(-MSb22 + MSt22)*MUEr2 - 
     &                    8*Atr**2*MSb12*
     &                     (-(ctb2*MSb22) + c2tb*MSt22 + 
     &                       MSb12*stb2)) + 
     &                 MSb12**2*
     &                  (4*ctt2*(MSb22 - MSt22)*
     &                     (2*c2tb*MSb22 - MUEr2*stb2) - 
     &                    8*MSb22*pMT2*stb2*stt2) + 
     &                 MSb12*
     &                  (-4*ctt2*MSt22**2*MUEr2*stb2 + 
     &                    8*Atr*MSb12*MSb22*pMT*s2tt*stb2 + 
     &                    4*ctb2*MSb22**2*
     &                     (ctt2*MUEr2 - 2*Atr*pMT*s2tt + 2*pMT2*stt2)-
     &                      4*c2tb*MSb22*MSt22*
     &                     (ctt2*(2*MSb22 - 2*MSt22 + MUEr2) - 
     &                       2*Atr*pMT*s2tt + 2*pMT2*stt2)))))/
     &          (MSb22 - MSt22) - 
     &         (2*hb*ht*MUEr*s2tb2*
     &            ((-2*Atr*MSb12*MSb22*(MSb2sum - 2*MSt22) + 
     &                 Abr*MSb2sum*(MSb12 - MSt22)*(-MSb22 + MSt22))
     &                *s2tt + 
     &              2*pMT*(ctt2*MSb22*MSt22*(-MSb22 + MSt22) + 
     &                 MSb12**2*(MSb22 - ctt2*MSt22) + 
     &                 MSb12*
     &                  (MSb22**2 + ctt2*MSt22**2 + 
     &                    MSb22*(MSb2sum*stt2 - 2*MSt22*(1 + stt2))
     &                    ))))/(-MSb22 + MSt22))*(-x10 + x19))/
     &     (MSb12 - MSt22)

        tmp23 = 1/4.D0*((Abr**2*ctt2*hb2*pMT*s2tb*
     &         (MA02*MSb2diff + 
     &           c2tb*(3*MA02*MSb12 + 2*MSb12**2 - MA02*MSb22 - 
     &              3*MSb12*MSt12) + MSt12*(MSb12 - 2*MSb22*stb2))+
     &          Abr*s2tt*(hb*ht*pMT*MUEr*
     &            (2*c2tb2*MSb12*(MA02 + MSb12 - MSt12) + 
     &              MSb2diff*(MA02 - MSt12)*s2tb2) + 
     &           hb2*s2tb*(c2tb*MSb12*(2*MSb12 - 3*MSt12)*pMT2 + 
     &              MSb12*MSt12*pMT2 + MA02**2*MSb2diff*stb2 + 
     &              2*MSb12**2*MSb2diff*stb2 - 
     &              3*MSb12*MSb2diff*MSt12*stb2 + 
     &              MSb2diff*MSt12**2*stb2 - 
     &              2*MSb22*MSt12*pMT2*stb2 - 
     &              MA02*(2*(MSb2diff*MSt12 - MSb22*pMT2)*stb2 + 
     &                 MSb12*(pMT2 - 3*c2tb*pMT2 + 3*MSb2diff*stb2)))
     &           ) + (hb*ht*MUEr*
     &            (4*c2tb2*MSb12**2*pMT2 - 4*MSb12*MSt12*pMT2 - 
     &              2*MSb12**3*s2tb2 + 2*MSb12**2*MSb22*s2tb2 + 
     &              MA02**2*MSb2diff*s2tb2 + 
     &              3*MSb12**2*MSt12*s2tb2 - 
     &              3*MSb12*MSb22*MSt12*s2tb2 - 
     &              MSb12*MSt12**2*s2tb2 + MSb22*MSt12**2*s2tb2 + 
     &              6*MSb12*MSt12*pMT2*s2tb2 - 
     &              2*MSb22*MSt12*pMT2*s2tb2 + 
     &              MA02*(-2*MSb2diff*MSt12*s2tb2 + 
     &                 2*MSb22*pMT2*s2tb2 + 
     &                 MSb12*
     &                  (4*pMT2 - 3*MSb2diff*s2tb2 - 6*pMT2*s2tb2)))+
     &             pMT*s2tb*
     &            (ht2*(MA02*MSb2diff + 
     &                 (MSb12 - 2*ctb2*MSb22)*MSt12 + 
     &                 c2tb*
     &                  (MA02*(-2*MSb12 + MSb2diff) + 
     &                    MSb12*(-2*MSb12 + 3*MSt12)))*MUEr2 + 
     &              hb2*(c2tb*MSb12*(2*MSb12 - 3*MSt12)*pMT2 + 
     &                 MSb12*MSt12*pMT2 + 2*MA02**2*MSb2diff*stb2 + 
     &                 4*MSb12**2*MSb2diff*stb2 - 
     &                 6*MSb12*MSb2diff*MSt12*stb2 + 
     &                 2*MSb2diff*MSt12**2*stb2 - 
     &                 2*MSb22*MSt12*pMT2*stb2 + 
     &                 MA02*
     &                  (6*MSb12**2*stb2 + 
     &                    2*MSb22*(-2*MSt12 + pMT2)*stb2 + 
     &                    MSb12*
     &                     (-pMT2 + 3*c2tb*pMT2 - 6*MSb22*stb2 + 
     &                       4*MSt12*stb2)))))*stt2)*x28)/
     &    (MSt12*pMT)

        tmp23 = tmp23 + 1/8.D0*
     &    ((2*ctt2*ht2*(MA02*MSb2diff + 
     &            (MSb12 - 2*ctb2*MSb22)*MSt22 + 
     &            c2tb*(MA02*(-2*MSb12 + MSb2diff) + 
     &               MSb12*(-2*MSb12 + 3*MSt22)))*MUEr2*s2tb + 
     &         (2*hb*ht*MUEr*
     &             (ctt2*(4*MSb12*(MA02 + c2tb2*MSb12 - MSt22)*
     &                   pMT2 + 
     &                  (MA02**2*MSb2diff + 
     &                     MSb12*
     &                      (2*MSb12*MSb2diff - 
     &                       MSt22*(3*MSb2diff + MSt22 - 6*pMT2)) - 
     &                     MA02*
     &                      (2*MSb2diff*MSt22 - 2*MSb22*pMT2 + 
     &                       3*MSb12*(MSb2diff + 2*pMT2)))*s2tb2) + 
     &               Abr*MSb12*pMT*
     &                (-2*c2tb2*MSb12 + MSt22*(2 - 3*s2tb2) + 
     &                  MA02*(-2 + 3*s2tb2))*s2tt + 
     &               MSb22*s2tb2*
     &                (ctt2*MSt22*(MSt22 - 2*pMT2) + 
     &                  Abr*(-MA02 + MSt22)*pMT*s2tt)) - 
     &            hb2*s2tb*
     &             (2*(MA02**2*MSb2diff*(-2*ctt2*pMT + Abr*s2tt)*
     &                   stb2 + 
     &                  Abr*s2tt*
     &                   (MSb12*
     &                      (2*c2tb*MSb12 + MSt22 - 3*c2tb*MSt22)*
     &                      pMT2 + 
     &                     (MSb2diff*
     &                       (2*MSb12**2 - 3*MSb12*MSt22 + 
     &                       MSt22**2) - 2*MSb22*MSt22*pMT2)*stb2))-
     &                 2*pMT*
     &                (ctt2*
     &                   (4*MSb12**2*MSb2diff*stb2 + 
     &                     2*MSt22*(MSb2diff*MSt22 - MSb22*pMT2)*
     &                      stb2 + 
     &                     MSb12*MSt22*(pMT2 - 6*MSb2diff*stb2)) + 
     &                  Abr**2*MSt22*(MSb12 - 2*MSb22*stb2)*stt2 + 
     &                  c2tb*MSb12*(2*MSb12 - 3*MSt22)*
     &                   (ctt2*pMT2 + Abr**2*stt2)) + 
     &               MA02*(-2*Abr*s2tt*
     &                   (2*(MSb2diff*MSt22 - MSb22*pMT2)*stb2 + 
     &                     MSb12*
     &                      (pMT2 - 3*c2tb*pMT2 + 3*MSb2diff*stb2))+
     &                    pMT*
     &                   (-2*ctt2*
     &                      (6*MSb12**2*stb2 + 
     &                       2*MSb22*(-2*MSt22 + pMT2)*stb2 + 
     &                       MSb12*
     &                       (-pMT2 + 3*c2tb*pMT2 - 6*MSb22*stb2 + 
     &                       4*MSt22*stb2)) - 
     &                     2*Abr**2*
     &                      (3*c2tb*MSb12 - c2tb*MSb22 + MSb2diff)*
     &                      stt2))))/pMT)*x29)/MSt22

        tmp24 = -(1/4.D0*((Abr**2*ctt2*hb2*
     &           (MA02*MSb2diff + 
     &             c2tb*(MA02*(2*MSb22 + MSb2diff) + 
     &                MSb22*(2*MSb22 - 3*MSt12)) + 
     &             2*ctb2*MSb12*MSt12 - MSb22*MSt12)*pMT*s2tb + 
     &          Abr*(hb2*(MSb22*
     &                 (MA02 + 3*c2tb*MA02 + 2*c2tb*MSb22 - 
     &                   MSt12 - 3*c2tb*MSt12)*pMT2 + 
     &                ctb2*
     &                 (MA02**2*MSb2diff + 2*MSb22**2*MSb2diff - 
     &                   3*MSb22*MSb2diff*MSt12 + 
     &                   MSb2diff*MSt12**2 + 2*MSb12*MSt12*pMT2 - 
     &                   MA02*
     &                    (3*MSb22*MSb2diff + 2*MSb2diff*MSt12 + 
     &                      2*MSb12*pMT2)))*s2tb + 
     &             hb*ht*pMT*MUEr*
     &              (2*c2tb2*MSb22**2 - 2*MSb22*MSt12 - 
     &                MSb12*MSt12*s2tb2 + 3*MSb22*MSt12*s2tb2 + 
     &                MA02*(2*MSb22 + MSb12*s2tb2 - 3*MSb22*s2tb2))
     &             )*s2tt + 
     &          (hb*ht*MUEr*
     &              (4*c2tb2*MSb22**2*pMT2 - 4*MSb22*MSt12*pMT2 + 
     &                2*MSb12*MSb22**2*s2tb2 - 2*MSb22**3*s2tb2 - 
     &                MA02**2*MSb2diff*s2tb2 - 
     &                3*MSb12*MSb22*MSt12*s2tb2 + 
     &                3*MSb22**2*MSt12*s2tb2 + 
     &                MSb12*MSt12**2*s2tb2 - 
     &                MSb22*MSt12**2*s2tb2 - 
     &                2*MSb12*MSt12*pMT2*s2tb2 + 
     &                6*MSb22*MSt12*pMT2*s2tb2 + 
     &                MA02*
     &                 (2*(MSb2diff*MSt12 + MSb12*pMT2)*s2tb2 + 
     &                   MSb22*
     &                    (4*pMT2 + 3*MSb2diff*s2tb2 - 6*pMT2*s2tb2))
     &                ) + pMT*s2tb*
     &              (hb2*MSb22*
     &                 (MA02 + 3*c2tb*MA02 + 2*c2tb*MSb22 - 
     &                   MSt12 - 3*c2tb*MSt12)*pMT2 + 
     &                2*ctb2*hb2*
     &                 (MA02**2*MSb2diff + 2*MSb22**2*MSb2diff - 
     &                   3*MSb22*MSb2diff*MSt12 + 
     &                   MSb2diff*MSt12**2 + MSb12*MSt12*pMT2 - 
     &                   MA02*
     &                    (3*MSb22*MSb2diff + 2*MSb2diff*MSt12 + 
     &                      MSb12*pMT2)) + 
     &                ht2*MUEr2*
     &                 (-(MSb22*MSt12) + 
     &                   c2tb*MSb22*
     &                    (-2*MA02 - 2*MSb22 + 3*MSt12) + 
     &                   2*MA02*MSb2diff*stb2 + 2*MSb12*MSt12*stb2)
     &                ))*stt2)*x30)/(MSt12*pMT))

        tmp24 = tmp24 + 1/4.D0*
     &    ((ctt2*(hb*ht*MUEr*
     &             (-4*c2tb2*MSb22**2*pMT2 + 4*MSb22*MSt22*pMT2 + 
     &               MA02**2*MSb2diff*s2tb2 + 
     &               2*MSb22**2*MSb2diff*s2tb2 - 
     &               3*MSb22*MSb2diff*MSt22*s2tb2 - 
     &               MSb12*MSt22**2*s2tb2 + MSb22*MSt22**2*s2tb2 + 
     &               2*MSb12*MSt22*pMT2*s2tb2 - 
     &               6*MSb22*MSt22*pMT2*s2tb2 - 
     &               MA02*(2*(MSb2diff*MSt22 + MSb12*pMT2)*s2tb2 + 
     &                  MSb22*
     &                   (4*pMT2 + 3*MSb2diff*s2tb2 - 6*pMT2*s2tb2)))
     &              + pMT*s2tb*
     &             (hb2*MSb22*
     &                (-((1 + 3*c2tb)*MA02) - 2*c2tb*MSb22 + 
     &                  MSt22 + 3*c2tb*MSt22)*pMT2 - 
     &               2*ctb2*hb2*
     &                (MA02**2*MSb2diff + 2*MSb22**2*MSb2diff - 
     &                  3*MSb22*MSb2diff*MSt22 + 
     &                  MSt22*(MSb2diff*MSt22 + MSb12*pMT2) - 
     &                  MA02*
     &                   (3*MSb22*MSb2diff + 2*MSb2diff*MSt22 + 
     &                     MSb12*pMT2)) + 
     &               ht2*MUEr2*
     &                (c2tb*MSb22*(2*MA02 + 2*MSb22 - 3*MSt22) + 
     &                  MSb22*MSt22 - 
     &                  2*(MA02*MSb2diff + MSb12*MSt22)*stb2))) + 
     &         Abr*(hb*ht*pMT*MUEr*
     &             (2*c2tb2*MSb22**2 + 
     &               MA02*(MSb22*(2 - 3*s2tb2) + MSb12*s2tb2) + 
     &               MSt22*(-2*MSb22 - MSb12*s2tb2 + 3*MSb22*s2tb2)
     &               )*s2tt + 
     &            hb2*s2tb*
     &             (MSb22*((2*c2tb*MSb22 - MSt22 - 3*c2tb*MSt22)*
     &                   (pMT2*s2tt - Abr*pMT*stt2) + 
     &                  MA02*
     &                   (pMT2*(s2tt + 3*c2tb*s2tt) - 
     &                     2*Abr*c2tb*pMT*stt2)) + 
     &               ctb2*(MA02**2*MSb2diff*s2tt + 
     &                  2*MSb22**2*MSb2diff*s2tt - 
     &                  3*MSb22*MSb2diff*MSt22*s2tt + 
     &                  MSt22*
     &                   (MSb2diff*MSt22*s2tt + 2*MSb12*pMT2*s2tt - 
     &                     2*Abr*MSb12*pMT*stt2) - 
     &                  MA02*
     &                   (3*MSb22*MSb2diff*s2tt + 
     &                     2*
     &                      (MSb2diff*MSt22*s2tt + 
     &                       MSb12*pMT2*s2tt + Abr*MSb2diff*pMT*stt2))
     &                  ))))*x31)/(MSt22*pMT)

        tmp25 = tmp24/MSb22 + 
     &   (x10 - x21)**2*(1/4.D0*
     &       (hb*ht*MUEr*(4*c2tb2*pMT2 + 
     &            (MSb2sum + c2tt*MSt2diff - MSt2sum)*s2tb2))/pMT + 
     &      s2tb*(1/2.D0*(c2tb*hb2*MUEr2) - 
     &         1/4.D0*(ht2*(MSb12 - 2*ctb2*MSb22 + 
     &              c2tb*(-MSb12 - c2tt*MSt2diff + MSt2sum + 
     &                 2*(Atr**2 + pMT2) - (Atr*MSt2diff*s2tt)/pMT)))))
     &     + x32*((hb*ht*(-2*pMT2*(MSb12 + c2tb2*(-pMT2 + MUEr2)) + 
     &           1/2.D0*((3*MSb22*pMT2 + MSb12*(-MSb2diff + pMT2) + 
     &                MSb2diff*MUEr2)*s2tb2)))/(pMT*MUEr) + 
     &      1/2.D0*(s2tb*(hb2*
     &             (c2tb*pMT2**2 - c2tb*MUEr2**2 + 
     &               3*MSb22*pMT2*stb2 + MSb2diff*MUEr2*stb2 + 
     &               MSb12**2*(1 + stb2) - 
     &               MSb12*(pMT2 + ctb2*pMT2 + 3*MSb22*stb2)) + 
     &            ht2*((1 + ctb2)*MSb12**2 + 3*ctb2*MSb22*pMT2 - 
     &               c2tb*pMT2**2 + ctb2*MSb2diff*MUEr2 + 
     &               c2tb*MUEr2**2 - 
     &               MSb12*(3*ctb2*MSb22 + pMT2 + pMT2*stb2))))/MUEr2)
     &     + x35*(1/2.D0*(hb*ht*MUEr*
     &          (-2*Atr*c2tb2*pMT*s2tt + 
     &            (4*c2tb2*pMT2 + (MSb2sum - 2*MSt22)*s2tb2)*stt2))/
     &        pMT + s2tb*(c2tb*hb2*MUEr2*stt2 - 
     &         1/4.D0*(ht2*(-2*Atr*s2tt*
     &                (-(ctb2*MSb22) + c2tb*(MSt22 + 2*pMT2) + 
     &                  MSb12*stb2) + 
     &               4*pMT*(Atr**2*c2tb*ctt2 + 
     &                  (-(ctb2*MSb22) + c2tb*(MSt22 + pMT2) + 
     &                     MSb12*stb2)*stt2)))/pMT)) + 
     &   x33*(c2tb*ctt2*hb2*MUEr2*s2tb + 
     &      1/2.D0*(hb*ht*MUEr*
     &           (ctt2*(MSb2sum - 2*MSt12)*s2tb2 + 
     &             2*c2tb2*(2*ctt2*pMT2 + Atr*pMT*s2tt)) - 
     &          ht2*s2tb*(-(ctb2*MSb22*(2*ctt2*pMT + Atr*s2tt)) + 
     &             MSb12*(2*ctt2*pMT + Atr*s2tt)*stb2 + 
     &             c2tb*(2*ctt2*pMT*(MSt12 + pMT2) + 
     &                Atr*(MSt12*s2tt + 2*pMT2*s2tt + 2*Atr*pMT*stt2)))
     &          )/pMT) + 1/4.D0*
     &    ((-(hb*ht*MUEr*(MSt12 - MSt2sum + c2tb2*(MSt22 + 4*pMT2) + 
     &              MSb2sum*s2tb2 - 2*ctt2*MSt12*s2tb2 + 
     &              c2tt*MSt22*s2tb2)) + 
     &         s2tb*(2*Atr**2*c2tb*ht2*pMT - 
     &            Atr*c2tb*ht2*MSt2diff*s2tt + 
     &            pMT*(-2*c2tb*hb2*MUEr2 + 
     &               ht2*(-2*ctb2*MSb22 + 
     &                  c2tb*
     &                   (-(c2tt*MSt2diff) + MSt2sum + 2*pMT2) + 
     &                  2*MSb12*stb2))))*(x10 - x22)**2)/pMT + 
     &   1/2.D0*((-(hb*ht*MUEr*
     &            (ctt2*(MSb2sum - 2*MSt12)*s2tb2 + 
     &              2*c2tb2*(2*ctt2*pMT2 + Atr*pMT*s2tt))) + 
     &         s2tb*(-(ctb2*ht2*MSb22*(2*ctt2*pMT + Atr*s2tt)) + 
     &            ht2*MSb12*(2*ctt2*pMT + Atr*s2tt)*stb2 + 
     &            c2tb*(2*ctt2*pMT*(ht2*(MSt12 + pMT2) - hb2*MUEr2) + 
     &               Atr*ht2*
     &                (MSt12*s2tt + 2*pMT2*s2tt + 2*Atr*pMT*stt2))))*
     &       x34)/pMT

        tmp25 = tmp25 + x37*
     &    ((hb*ht*(2*pMT2*(MSb22 + c2tb2*(-pMT2 + MUEr2)) - 
     &           1/2.D0*((3*MSb12*pMT2 + MSb22*(MSb2diff + pMT2) - 
     &                MSb2diff*MUEr2)*s2tb2)))/(pMT*MUEr) + 
     &      1/2.D0*(s2tb*(-(hb2*
     &               (-MSb22**2 + c2tb*(pMT2 - MUEr2)*(pMT2 + MUEr2) + 
     &                 ctb2*
     &                  (-MSb22**2 + 3*MSb12*(MSb22 - pMT2) + 
     &                    MSb2diff*MUEr2) + MSb22*pMT2*(1 + stb2)))+
     &              ht2*(c2tb*pMT2**2 - c2tb*MUEr2**2 + 
     &               3*MSb12*pMT2*stb2 - MSb2diff*MUEr2*stb2 + 
     &               MSb22**2*(1 + stb2) - 
     &               MSb22*(pMT2 + ctb2*pMT2 + 3*MSb12*stb2))))/MUEr2)
     &     + 1/2.D0*((hb*ht*MUEr*
     &          (2*Atr*c2tb2*pMT*s2tt - 
     &            (4*c2tb2*pMT2 + (MSb2sum - 2*MSt22)*s2tb2)*stt2)+
     &           s2tb*(2*Atr**2*c2tb*ctt2*ht2*pMT + 
     &            Atr*ht2*s2tt*
     &             (ctb2*MSb22 - c2tb*(MSt22 + 2*pMT2) - MSb12*stb2)
     &              + 2*pMT*
     &             (-(ctb2*ht2*MSb22) + 
     &               c2tb*(ht2*MSt22 + ht2*pMT2 - hb2*MUEr2) + 
     &               ht2*MSb12*stb2)*stt2))*x36)/pMT

        s2tbF4b = (tmp22/MSb22 + tmp23/MSb2diff)/MSb12 + 
     &   tmp25/MSb2diff + c2tb*(hb2 - ht2)*s2tb*(-x10 + x20)

        s2tbF4b = s2tbF4b + 
     &   1/8.D0*((-x10 + x21)*
     &      (((hb*ht*MUEr*(-4*c2tb2*
     &                  (2*(1 + ctt2)*pMT2 + (Abr + 2*Atr)*pMT*s2tt) + 
     &                 2*s2tb2*
     &                  (2*(1 + ctt2)*MSt12 + 
     &                    (-2*MA02 + MSb2sum)*stt2)) - 
     &              2*ht2*s2tb*
     &               (4*pMT*(-(ctb2*MSb22) + MSb12*stb2)*stt2 - 
     &                 2*c2tb*
     &                  (2*ctt2*pMT*(MSt12 + pMT2) + 
     &                    Atr*(MSt12 + 2*pMT2)*s2tt + 
     &                    2*Atr**2*pMT*stt2 + pMT*MUEr2*stt2)) - 
     &              hb2*s2tb*
     &               (2*Abr*s2tt*
     &                  (ctb2*MSb12 + 
     &                    c2tb*(-MA02 + MSt12 + 2*pMT2) - MSb22*stb2
     &                    ) + 
     &                 4*pMT*
     &                  (ctb2*(1 + ctt2)*MSb12 - 
     &                    (1 + ctt2)*MSb22*stb2 + 
     &                    c2tb*
     &                     (Abr**2*ctt2 + 2*ctt2*MUEr2 + 
     &                       (-MA02 + MSt12 + pMT2)*stt2))))*
     &            (-x10 + x18) + 
     &           (-(hb2*s2tb*
     &                 (c2tb*
     &                    (4*ctt2*pMT*(-MA02 + MSt22 + pMT2) + 
     &                      2*Abr*(MA02 - MSt22 - 2*pMT2)*s2tt + 
     &                      4*pMT*(Abr**2 + 2*MUEr2)*stt2) + 
     &                   2*(ctb2*MSb12 - MSb22*stb2)*
     &                    (-(Abr*s2tt) + 2*pMT*(1 + stt2)))) - 
     &              hb*ht*MUEr*
     &               (-2*s2tb2*
     &                  (ctt2*(-2*MA02 + MSb2sum) + 
     &                    2*MSt22*(1 + stt2)) + 
     &                 c2tb2*
     &                  (-4*(Abr + 2*Atr)*pMT*s2tt + 8*pMT2*(1 + stt2))
     &                 ) - 
     &              2*ht2*s2tb*
     &               (-2*ctt2*pMT*
     &                  (2*ctb2*MSb22 + c2tb*MUEr2 - 2*MSb12*stb2)+
     &                   c2tb*
     &                  (2*Atr*(MSt22 + 2*pMT2)*s2tt - 
     &                    4*pMT*(Atr**2*ctt2 + (MSt22 + pMT2)*stt2))))
     &             *(-x10 + x19))/(MSb2diff*pMT) + 
     &        4*(((hb*ht*MUEr*s2tb2*
     &                 (-2*(MSb12**2 + MSt12*MSt22) + 
     &                   2*MSb12*(MSt2sum - pMT2) + 
     &                   (c2tt*MSt2diff + MSt2sum)*pMT2 + 
     &                   Atr*MSt2diff*pMT*s2tt))/pMT + 
     &              (s2tb*(ctb2*ht2*
     &                    (MSb12**3 - MSt12*MSt22*(pMT2 + MUEr2) - 
     &                      MSb12**2*(MSt2sum + pMT2 + MUEr2) + 
     &                      MSb12*
     &                       (MSt12*MSt22 - 2*Atr**2*pMT2 + 
     &                       MSt2sum*pMT2 - 2*pMT2**2 + MSt2sum*MUEr2)
     &                        + pMT2*
     &                       (Atr**2*(-(c2tt*MSt2diff) + MSt2sum) + 
     &                       (c2tt*MSt2diff + MSt2sum)*pMT2 + 
     &                       2*Atr*MSt2diff*pMT*s2tt)) + 
     &                   hb2*
     &                    (MSb12**3 + 
     &                      (c2tt*MSt2diff + MSt2sum)*pMT2*MUEr2 - 
     &                      MSt12*MSt22*(pMT2 + MUEr2) - 
     &                      MSb12**2*(MSt2sum + pMT2 + MUEr2) + 
     &                      MSb12*
     &                       (MSt12*MSt22 - 2*pMT2*MUEr2 + 
     &                       MSt2sum*(pMT2 + MUEr2)))*stb2))/pMT2)/
     &            ((MSb12 - MSt12)*(MSb12 - MSt22)) + 
     &           (MSb12 + pMT2 - MUEr2)*(-MSb12 + pMT2 + MUEr2)*
     &            ((-2*hb*ht*MUEr*s2tb2)/pMT - 
     &              ((-MSb12 + pMT2 + MUEr2)*s2tb*
     &                 (ctb2*ht2 + hb2*stb2))/pMT2)*x23 + 
     &           ((hb*ht*MUEr*
     &                  (4*pMT2 - (MSb2sum + 6*pMT2 - 2*MUEr2)*s2tb2)+
     &                   pMT*s2tb*
     &                  (hb2*
     &                     (-MSb2diff + c2tb*(MSb2sum + 2*pMT2)) - 
     &                    2*ht2*
     &                     (ctb2*MSb22 + c2tb*pMT2 - MSb12*stb2)))*
     &               (-x10 + x20) + 
     &              (-(hb*ht*MUEr*
     &                    (pMT2*(4 - 6*s2tb2) + 
     &                      (MSb2sum + 2*MUEr2)*s2tb2)) + 
     &                 pMT*s2tb*
     &                  (-2*ctb2*ht2*MSb22 - hb2*MSb2diff + 
     &                    c2tb*hb2*MSb2sum - 2*c2tb*hb2*pMT2 + 
     &                    2*c2tb*ht2*pMT2 + 2*ht2*MSb12*stb2))*
     &               (-x10 + x24))/(MSb2diff*pMT))))

        s2tbF4b = s2tbF4b + 
     &   1/8.D0*((-x10 + x22)*
     &      (((-2*ht2*s2tb*
     &               (c2tb*
     &                  (4*ctt2*pMT*(MSt12 + pMT2) + 
     &                    2*Atr*(MSt12 + 2*pMT2)*s2tt) + 
     &                 2*pMT*
     &                  (2*Atr**2*c2tb + 2*ctb2*MSb22 + c2tb*MUEr2 - 
     &                    2*MSb12*stb2)*stt2) - 
     &              hb*ht*MUEr*
     &               (-4*c2tb2*
     &                  (2*(1 + ctt2)*pMT2 + (Abr + 2*Atr)*pMT*s2tt) + 
     &                 2*s2tb2*
     &                  (2*(1 + ctt2)*MSt12 + 
     &                    (-2*MA02 + MSb2sum)*stt2)) - 
     &              hb2*s2tb*
     &               (-4*Abr**2*c2tb*ctt2*pMT + 
     &                 2*Abr*s2tt*
     &                  (-(ctb2*MSb12) + 
     &                    c2tb*(MA02 - MSt12 - 2*pMT2) + MSb22*stb2)
     &                   + pMT*
     &                  (-4*ctb2*(1 + ctt2)*MSb12 + 
     &                    4*(1 + ctt2)*MSb22*stb2 - 
     &                    4*c2tb*
     &                     (2*ctt2*MUEr2 + 
     &                       (-MA02 + MSt12 + pMT2)*stt2))))*
     &            (-x10 + x18) + 
     &           (-2*ht2*s2tb*
     &               (2*ctt2*pMT*
     &                  (2*ctb2*MSb22 + c2tb*MUEr2 - 2*MSb12*stb2)+
     &                   c2tb*
     &                  (4*Atr**2*ctt2*pMT - 
     &                    2*Atr*(MSt22 + 2*pMT2)*s2tt + 
     &                    4*pMT*(MSt22 + pMT2)*stt2)) + 
     &              hb*ht*MUEr*
     &               (-2*s2tb2*
     &                  (ctt2*(-2*MA02 + MSb2sum) + 
     &                    2*MSt22*(1 + stt2)) + 
     &                 c2tb2*
     &                  (-4*(Abr + 2*Atr)*pMT*s2tt + 8*pMT2*(1 + stt2))
     &                 ) - 
     &              hb2*s2tb*
     &               (2*(ctb2*MSb12 - MSb22*stb2)*
     &                  (Abr*s2tt - 2*pMT*(1 + stt2)) - 
     &                 2*c2tb*
     &                  (Abr*(MA02 - MSt22 - 2*pMT2)*s2tt + 
     &                    2*pMT*
     &                     (ctt2*(-MA02 + MSt22 + pMT2) + 
     &                       (Abr**2 + 2*MUEr2)*stt2))))*(-x10 + x19)
     &           )/(MSb2diff*pMT) + 
     &        4*(((hb*ht*MUEr*s2tb2*
     &                 (-2*(MSb22**2 + MSt12*MSt22) + 
     &                   2*MSb22*(MSt2sum - pMT2) + 
     &                   (c2tt*MSt2diff + MSt2sum)*pMT2 + 
     &                   Atr*MSt2diff*pMT*s2tt))/pMT + 
     &              (s2tb*(ctb2*hb2*
     &                    ((MSb22**2 + MSt12*MSt22 - 
     &                      MSb22*MSt2sum)*(-MSb22 + pMT2) + 
     &                      (MSb22**2 + MSt12*MSt22 - 
     &                       MSb22*(MSt2sum - 2*pMT2) - 
     &                       (c2tt*MSt2diff + MSt2sum)*pMT2)*MUEr2)+
     &                     ht2*
     &                    (-MSb22**3 - MSt2sum*pMT2*(Atr**2 + pMT2) + 
     &                      MSt12*MSt22*(pMT2 + MUEr2) + 
     &                      MSb22**2*(MSt2sum + pMT2 + MUEr2) - 
     &                      MSb22*
     &                       (MSt12*MSt22 - 2*pMT2*(Atr**2 + pMT2) + 
     &                       MSt2sum*(pMT2 + MUEr2)) - 
     &                      MSt2diff*pMT2*
     &                       (-(Atr**2*c2tt) + c2tt*pMT2 + 
     &                       2*Atr*pMT*s2tt))*stb2))/pMT2)/
     &            ((MSb22 - MSt12)*(MSb22 - MSt22)) + 
     &           ((hb*ht*MUEr*
     &                  (-4*pMT2 + (MSb2sum + 6*pMT2 - 2*MUEr2)*s2tb2)
     &                   + pMT*s2tb*
     &                  (hb2*(MSb2diff - c2tb*(MSb2sum + 2*pMT2)) + 
     &                    2*ht2*
     &                     (ctb2*MSb22 + c2tb*pMT2 - MSb12*stb2)))*
     &               (-x10 + x20) + 
     &              (hb*ht*MUEr*
     &                  (pMT2*(4 - 6*s2tb2) + 
     &                    (MSb2sum + 2*MUEr2)*s2tb2) + 
     &                 pMT*s2tb*
     &                  (2*ctb2*ht2*MSb22 + 
     &                    hb2*
     &                     (MSb2diff - c2tb*MSb2sum + 2*c2tb*pMT2)-
     &                      2*ht2*(c2tb*pMT2 + MSb12*stb2)))*
     &               (-x10 + x24))/(MSb2diff*pMT) + 
     &           (MSb22 + pMT2 - MUEr2)*(-MSb22 + pMT2 + MUEr2)*
     &            ((-2*hb*ht*MUEr*s2tb2)/pMT + 
     &              ((-MSb22 + pMT2 + MUEr2)*s2tb*
     &                 (ctb2*hb2 + ht2*stb2))/pMT2)*x25)))

        s2tbF4b = s2tbF4b + 
     &   (MA02 + MSb12 - MSt12)*x11*
     &    ((((MA02**2 + MSt12*(-MSb12 + MSt12) - 
     &              MA02*(MSb12 + 2*MSt12))*x28*
     &            (-(1/4.D0*(hb*ht*MUEr*s2tb2*
     &                    (Abr*pMT*s2tt + 2*pMT2*stt2))/pMT) - 
     &              1/2.D0*(s2tb*
     &                 (Abr**2*ctt2*hb2*stb2 + 
     &                   Abr*hb2*pMT*s2tt*stb2 + 
     &                   ctb2*ht2*MUEr2*stt2 + hb2*pMT2*stb2*stt2))))
     &           /MSt12 + 1/4.D0*
     &          ((MA02 - MSb12 - MSt12)*
     &            (2*Abr**2*ctt2*hb2*s2tb*stb2 + 
     &              Abr*s2tt*
     &               (hb*ht*MUEr*s2tb2 + 2*hb2*pMT*s2tb*stb2) + 
     &              2*(ctb2*ht2*MUEr2*s2tb + hb*ht*pMT*MUEr*s2tb2 + 
     &                 hb2*pMT2*s2tb*stb2)*stt2)*(x10 - x18)))/MSb12
     &        - 1/2.D0*((Abr*(hb*ht*MUEr*s2tb2*s2tt + 
     &              2*hb2*s2tb*(Abr*ctt2 + pMT*s2tt)*stb2) + 
     &           2*(ctb2*ht2*MUEr2*s2tb + hb*ht*pMT*MUEr*s2tb2 + 
     &              hb2*pMT2*s2tb*stb2)*stt2)*(-x10 + x21))) + 
     &   1/4.D0*((MA02 + MSb12 - MSt22)*x12*
     &      (2*(hb*ht*MUEr*s2tb2*(-2*ctt2*pMT + Abr*s2tt) - 
     &           2*s2tb*(ctb2*ctt2*ht2*MUEr2 + 
     &              hb2*stb2*(ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))
     &              ))*(-x10 + x21) + 
     &        ((MA02 - MSb12 - MSt22)*
     &            (2*ctb2*ctt2*ht2*MUEr2*s2tb + 
     &              2*ctt2*
     &               (hb*ht*pMT*MUEr*s2tb2 + hb2*pMT2*s2tb*stb2) + 
     &              Abr*(-(hb*ht*MUEr*s2tb2*s2tt) + 
     &                 2*hb2*s2tb*stb2*(-(pMT*s2tt) + Abr*stt2)))*
     &            (x10 - x19) + 
     &           ((MA02**2 + MSt22*(-MSb12 + MSt22) - 
     &                MA02*(MSb12 + 2*MSt22))*
     &              ((hb*ht*MUEr*s2tb2*(-2*ctt2*pMT2 + Abr*pMT*s2tt))/
     &                 pMT - 
     &                2*s2tb*
     &                 (ctb2*ctt2*ht2*MUEr2 + 
     &                   hb2*stb2*
     &                    (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))))*
     &              x29)/MSt22)/MSb12)) - 
     &   1/4.D0*((MA02 + MSb22 - MSt12)*
     &       (2*Abr**2*ctb2*ctt2*hb2*s2tb + 
     &         Abr*(2*ctb2*hb2*pMT*s2tb - hb*ht*MUEr*s2tb2)*s2tt + 
     &         2*(ctb2*hb2*pMT2*s2tb - hb*ht*pMT*MUEr*s2tb2 + 
     &            ht2*MUEr2*s2tb*stb2)*stt2)*x13*
     &       (-(MA02**2*x30) + 
     &         MSt12*(-(MSt12*(x10 - x18 + x30)) + 
     &            MSb22*(x10 + x18 - 2*x22 + x30)) + 
     &         MA02*(MSb22*x30 + MSt12*(x10 - x18 + 2*x30))))/
     &     (MSb22*MSt12)

        s2tbF4b = s2tbF4b - 
     &   1/4.D0*((MA02 + MSb22 - MSt22)*
     &       (-2*ctt2*hb*ht*pMT*MUEr*s2tb2 + 
     &         Abr*hb*ht*MUEr*s2tb2*s2tt + 
     &         2*ctt2*ht2*MUEr2*s2tb*stb2 + 
     &         2*ctb2*hb2*s2tb*
     &          (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*x14*
     &       (-(MA02**2*x31) + 
     &         MSt22*(-(MSt22*(x10 - x19 + x31)) + 
     &            MSb22*(x10 + x19 - 2*x22 + x31)) + 
     &         MA02*(MSb22*x31 + MSt22*(x10 - x19 + 2*x31))))/
     &     (MSb22*MSt22) + 
     &   x23*(1/2.D0*((-MSb12 + pMT2 + MUEr2)**2*
     &         ((-2*hb*ht*MUEr*s2tb2)/pMT - 
     &           ((-MSb12 + pMT2 + MUEr2)*s2tb*
     &              (ctb2*ht2 + hb2*stb2))/pMT2)*(-x10 + x24)) + 
     &      (-MSb12 + pMT2 + MUEr2)*
     &       ((2*hb*ht*pMT*MUEr*s2tb2 - 
     &            (MSb12 - pMT2 - MUEr2)*s2tb*(ctb2*ht2 + hb2*stb2))*
     &          (-x10 + x20) + 
     &         1/2.D0*((MSb12 - pMT2 + MUEr2)*
     &             ((2*hb*ht*pMT2*MUEr2*s2tb2)/(pMT*MUEr) + 
     &               (-MSb12 + pMT2 + MUEr2)*s2tb*
     &                (ctb2*ht2 + hb2*stb2))*x32)/MUEr2)) + 
     &   x25*(1/2.D0*((-MSb22 + pMT2 + MUEr2)**2*
     &         ((-2*hb*ht*MUEr*s2tb2)/pMT + 
     &           ((-MSb22 + pMT2 + MUEr2)*s2tb*
     &              (ctb2*hb2 + ht2*stb2))/pMT2)*(-x10 + x24)) + 
     &      (-MSb22 + pMT2 + MUEr2)*
     &       ((2*hb*ht*pMT*MUEr*s2tb2 + 
     &            (MSb22 - pMT2 - MUEr2)*s2tb*(ctb2*hb2 + ht2*stb2))*
     &          (-x10 + x20) - 
     &         1/2.D0*((MSb22 - pMT2 + MUEr2)*
     &             (-2*hb*ht*pMT2*MUEr2*s2tb2 + 
     &               pMT*MUEr*(-MSb22 + pMT2 + MUEr2)*s2tb*
     &                (ctb2*hb2 + ht2*stb2))*x37)/(pMT*MUEr*MUEr2)))

        s2tbF4b = s2tbF4b + 
     &   1/8.D0*(((-2*(2*c2tb*hb2*MSb12*pMT*s2tb + 
     &                hb*ht*(-2*MSb12 + MSb2diff)*MUEr*s2tb2 + 
     &                2*hb2*MSb2diff*pMT*s2tb*stb2)*
     &              (MSt22*stt2*x12*x28 + ctt2*MSt12*x11*x29) + 
     &             Abr*hb2*s2tb*s2tt*
     &              (-3*c2tb*MSb12*MSt22*x12*x28 + 
     &                c2tb*MSb22*MSt22*x12*x28 - 
     &                MSb2diff*MSt22*x12*x28 + 
     &                2*c2tb*MSb12*MSt12*x11*x29 + 
     &                2*MSb2diff*MSt12*stb2*x11*x29))/
     &           (MSb12*MSt12*MSt22*x11*x12) + 
     &          (((Abr*hb2*(MSb2diff + c2tb*(2*MSb22 + MSb2diff))*
     &                   s2tb*s2tt + 
     &                  2*(hb2*
     &                      (MSb2diff + c2tb*(2*MSb22 + MSb2diff))*
     &                      pMT*s2tb - 
     &                     hb*ht*(2*MSb22 + MSb2diff)*MUEr*s2tb2)*
     &                   stt2)*x30)/(MSt12*x13) - 
     &             (2*(ctt2*
     &                   (-2*hb2*(c2tb*MSb22 + ctb2*MSb2diff)*pMT*
     &                      s2tb + 
     &                     hb*ht*(2*MSb22 + MSb2diff)*MUEr*s2tb2) + 
     &                  Abr*hb2*(c2tb*MSb22 + ctb2*MSb2diff)*s2tb*
     &                   s2tt)*x31)/(MSt22*x14))/MSb22)/MSb2diff + 
     &       4*(((-(c2tb*(hb2 - ht2)*pMT*(pMT2 + MUEr2)*s2tb) + 
     &               4*hb*ht*pMT2*MUEr*s2tb2 - 
     &               MSb12*pMT*s2tb*(ctb2*ht2 + hb2*stb2) + 
     &               MSb22*pMT*s2tb*(ctb2*hb2 + ht2*stb2))*
     &             (-x10 + x24))/pMT2 + 
     &          ((c2tb*(hb2 - ht2)*pMT*MUEr*s2tb - 
     &               2*hb*ht*MUEr2*s2tb2)*(x25*x32 - x23*x37))/
     &           (MSb2diff*MUEr*MUEr2*x23*x25)))/pMT

        s2tbF4b = s2tbF4b + 
     &   1/4.D0*((2*(-((s2tb*
     &                (ctb2*ht2*MUEr2 + hb2*(Abr**2 + pMT2)*stb2))/
     &              MSb12) + 
     &           (-((hb*ht*MSb2sum*pMT*MUEr*s2tb2)/MSb12) + 
     &              s2tb*(ctb2*hb2*(Abr**2 + pMT2) + ht2*MUEr2*stb2))/
     &            MSb22) + 
     &        ((MA02 + MSb12 - MSt12)**2*
     &            (hb*ht*MUEr*s2tb2*(Abr*s2tt + 2*pMT*stt2) + 
     &              2*s2tb*
     &               (ctb2*ht2*MUEr2*stt2 + 
     &                 hb2*stb2*
     &                  (Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2)))*x11+
     &             (MA02 + MSb12 - MSt22)**2*
     &            (2*ctt2*(ctb2*ht2*MUEr2*s2tb + 
     &                 hb*ht*pMT*MUEr*s2tb2 + hb2*pMT2*s2tb*stb2) + 
     &              Abr*(-(hb*ht*MUEr*s2tb2*s2tt) + 
     &                 2*hb2*s2tb*stb2*(-(pMT*s2tt) + Abr*stt2)))*x12
     &           )/MSb12 + 
     &        ((MA02 + MSb22 - MSt12)**2*
     &            (hb*ht*MUEr*s2tb2*(Abr*s2tt + 2*pMT*stt2) - 
     &              2*s2tb*
     &               (ht2*MUEr2*stb2*stt2 + 
     &                 ctb2*hb2*
     &                  (Abr**2*ctt2 + Abr*pMT*s2tt + pMT2*stt2)))*x13-
     &             (MA02 + MSb22 - MSt22)**2*
     &            (hb*ht*MUEr*s2tb2*(-2*ctt2*pMT + Abr*s2tt) + 
     &              2*s2tb*
     &               (ctt2*ht2*MUEr2*stb2 + 
     &                 ctb2*hb2*
     &                  (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))))*
     &            x14)/MSb22 + 
     &        ((4*c2tb2*hb*ht*pMT2*MUEr + 2*Abr**2*c2tb*hb2*pMT*s2tb - 
     &             2*c2tb*hb2*MA02*pMT*s2tb - 
     &             2*ctb2*hb2*MSb12*pMT*s2tb + 
     &             c2tb*c2tt*hb2*MSt2diff*pMT*s2tb + 
     &             c2tb*hb2*MSt2sum*pMT*s2tb + 
     &             2*c2tb*hb2*pMT*pMT2*s2tb - 
     &             2*c2tb*ht2*pMT*MUEr2*s2tb + 
     &             2*hb*ht*MA02*MUEr*s2tb2 + 
     &             hb*ht*MSb2sum*MUEr*s2tb2 - 
     &             c2tt*hb*ht*MSt2diff*MUEr*s2tb2 - 
     &             hb*ht*MSt2sum*MUEr*s2tb2 - 
     &             Abr*c2tb*hb2*MSt2diff*s2tb*s2tt + 
     &             2*hb2*MSb22*pMT*s2tb*stb2)*(x21 - x22))/
     &         (MSb2diff*pMT))*(-x10 + x9))

        tmp26 = 1/4.D0*(-2*hb*ht*MUEr*s2tb2*s2tt*
     &       (-((MSt2diff + c2tt*(-2*MSb12 + MSt2sum))*pMT) + 
     &         Atr*(2*MSb12 - MSt2sum)*s2tt) + 
     &      (s2tb*s2tt*(-(hb2*
     &              (4*c2tt*
     &                 (ctb2*MSb12*
     &                    (MSb12**2 + MSt12*MSt22 - MSb12*MSt2sum)+
     &                     (-(MSb22*
     &                       (MSb12**2 + MSt12*MSt22 - 
     &                       MSb12*MSt2sum)) + MSb12*MSb2diff*MUEr2)
     &                     *stb2) + 
     &                4*MSb2diff*MUEr2*stb2*
     &                 (-(ctt2*MSt22) + MSt12*stt2))) - 
     &           ht2*(-2*Atr**2*ctb2*MSb2diff*
     &               (2*c2tt*MSb12 + MSt2diff - c2tt*MSt2sum) + 
     &              4*Atr*ctb2*MSb2diff*(2*MSb12 - MSt2sum)*pMT*
     &               s2tt - 4*c2tt*MSb12**3*stb2 + 
     &              MSb12*(2*ctb2*(MSt2diff + c2tt*MSt2sum)*pMT2 + 
     &                 4*c2tt*ctb2*MSb22*(-MSt2sum + pMT2) - 
     &                 4*c2tt*MSt12*MSt22*stb2) + 
     &              4*c2tt*MSb12**2*
     &               (ctb2*(MSb22 - pMT2) + MSt2sum*stb2) + 
     &              4*ctb2*MSb22*
     &               (c2tt*MSt12*MSt22 - ctt2*MSt22*pMT2 + 
     &                 MSt12*pMT2*stt2))))/MSb2diff)/
     &    ((MSb12 - MSt12)*(MSb12 - MSt22))

        tmp26 = tmp26 - 1/4.D0*
     &    (Abr**2*hb2*s2tb*s2tt*
     &        (c2tb*(MSt2diff + c2tt*(-2*MA02 + MSt2sum)) + 
     &          2*c2tt*(ctb2*MSb12 - MSb22*stb2))*(x18 - x19) + 
     &       s2tt*(2*c2tb2*hb*ht*
     &           (3*MSt2diff + c2tt*(2*MA02 + MSt2sum))*pMT*MUEr + 
     &          MSt2diff*s2tb*
     &           (2*ctb2*hb2*MSb12 - 2*ctb2*ht2*MSb22 + 
     &             c2tb*hb2*pMT2 - 2*c2tb*ht2*pMT2 + 
     &             2*c2tb*hb2*MUEr2 - c2tb*ht2*MUEr2 + 
     &             2*ht2*MSb12*stb2 - 2*hb2*MSb22*stb2) + 
     &          c2tt*(2*hb*ht*pMT*MUEr*(-2*MSb12 + MSb2sum*s2tb2) + 
     &             s2tb*(2*ctb2*
     &                 (hb2*MSb12*(MSt2sum - pMT2) + 
     &                   ht2*MSb22*(MSt2sum + MUEr2)) + 
     &                c2tb*
     &                 (ht2*
     &                    (-2*MSt2sum*pMT2 - 2*MA02*MUEr2 + 
     &                      MSt2sum*MUEr2) + 
     &                   hb2*
     &                    (2*MA02*pMT2 - MSt2sum*pMT2 + 
     &                      2*MSt2sum*MUEr2)) - 
     &                2*(hb2*MSb22*(MSt2sum - pMT2) + 
     &                   ht2*MSb12*(MSt2sum + MUEr2))*stb2)))*
     &        (x18 - x19) + 
     &       4*Atr**2*c2tb*ht2*s2tb*s2tt*(ctt2*MSt12 - MSt22*stt2)*
     &        (x18 - x19) + 
     &       2*Abr*(hb*ht*MUEr*
     &           (2*MSt12*s2tb2*(-x10 + x18) + 
     &             2*MSt22*s2tb2*(x10 - x19) + 
     &             (-2*MA02 + 2*c2tb2*c2tt2*MA02 + 
     &                c2tt2*MSb2sum*s2tb2 + 2*MSb12*s2tt2 + 
     &                c2tb2*MSt2sum*s2tt2)*(x18 - x19)) - 
     &          2*c2tt2*hb2*pMT*s2tb*(ctb2*MSb12 - MSb22*stb2)*
     &           (x18 - x19) + 
     &          c2tb*hb2*pMT*s2tb*
     &           (MSt2diff*s2tt2*(-2*x10 + x18 + x19) + 
     &             2*c2tt2*
     &              (-(MSt22*x10) + MSt12*(x10 - x18) + MA02*x18 - 
     &                MA02*x19 + MSt22*x19))) - 
     &       4*Atr*(hb*ht*MUEr*
     &           (2*MSt12*s2tb2*(x10 - x18) + 
     &             2*MSt22*s2tb2*(-x10 + x19) + 
     &             c2tb2*MSt2sum*s2tt2*(-x18 + x19)) + 
     &          c2tb*ht2*pMT*s2tb*
     &           (MSt2diff*s2tt2*(-2*x10 + x18 + x19) + 
     &             2*c2tt2*(MSt12*(x10 - x18) + MSt22*(-x10 + x19))
     &             )))/(MSb2diff*MSt2diff)

        tmp27 = 1/4.D0*(s2tt*
     &      (2*hb*ht*MUEr*s2tb2*
     &         (MSt2diff*pMT + c2tt*(-2*MSb22 + MSt2sum)*pMT + 
     &           Atr*(-2*MSb22 + MSt2sum)*s2tt) - 
     &        (2*s2tb*(ctb2*hb2*MSb2diff*MSt2diff*MUEr2 + 
     &             ht2*(Atr**2*MSb2diff*MSt2diff + 
     &                (-MSb12 + MSb22)*MSt2diff*pMT2 + 
     &                2*Atr*MSb2diff*(-2*MSb22 + MSt2sum)*pMT*s2tt)*
     &              stb2 + 
     &             c2tt*(ctb2*
     &                 (-2*ht2*MSb22*
     &                    (MSb22**2 + MSt12*MSt22 - MSb22*MSt2sum)+
     &                     hb2*
     &                    (-2*MSb12*
     &                       (MSb22**2 + MSt12*MSt22 - 
     &                       MSb22*MSt2sum) + 
     &                      MSb2diff*(-2*MSb22 + MSt2sum)*MUEr2)) + 
     &                (2*hb2*MSb22*
     &                    (MSb22**2 + MSt12*MSt22 - MSb22*MSt2sum)+
     &                     ht2*
     &                    ((2*MSb22 - MSt2sum)*
     &                       (Atr**2*MSb2diff - MSb22*pMT2) + 
     &                      MSb12*
     &                       (2*MSb22**2 + 2*MSt12*MSt22 - 
     &                       2*MSb22*MSt2sum + 2*MSb22*pMT2 - 
     &                       MSt2sum*pMT2)))*stb2)))/MSb2diff))/
     &    ((MSb22 - MSt12)*(MSb22 - MSt22))

        tmp27 = tmp27 + 1/4.D0*
     &    (Abr**2*hb2*s2tb*s2tt*
     &        (c2tb*(MSt2diff + c2tt*(-2*MA02 + MSt2sum)) + 
     &          2*c2tt*(ctb2*MSb12 - MSb22*stb2))*(x18 - x19) + 
     &       s2tt*(2*c2tb2*hb*ht*
     &           (3*MSt2diff + c2tt*(2*MA02 + MSt2sum))*pMT*MUEr + 
     &          MSt2diff*s2tb*
     &           (2*ctb2*hb2*MSb12 - 2*ctb2*ht2*MSb22 + 
     &             c2tb*hb2*pMT2 - 2*c2tb*ht2*pMT2 + 
     &             2*c2tb*hb2*MUEr2 - c2tb*ht2*MUEr2 + 
     &             2*ht2*MSb12*stb2 - 2*hb2*MSb22*stb2) + 
     &          c2tt*(2*hb*ht*pMT*MUEr*(-2*MSb22 + MSb2sum*s2tb2) + 
     &             s2tb*(2*ctb2*
     &                 (hb2*MSb12*(MSt2sum - pMT2) + 
     &                   ht2*MSb22*(MSt2sum + MUEr2)) + 
     &                c2tb*
     &                 (ht2*
     &                    (-2*MSt2sum*pMT2 - 2*MA02*MUEr2 + 
     &                      MSt2sum*MUEr2) + 
     &                   hb2*
     &                    (2*MA02*pMT2 - MSt2sum*pMT2 + 
     &                      2*MSt2sum*MUEr2)) - 
     &                2*(hb2*MSb22*(MSt2sum - pMT2) + 
     &                   ht2*MSb12*(MSt2sum + MUEr2))*stb2)))*
     &        (x18 - x19) + 
     &       4*Atr**2*c2tb*ht2*s2tb*s2tt*(ctt2*MSt12 - MSt22*stt2)*
     &        (x18 - x19) + 
     &       2*Abr*(hb*ht*MUEr*
     &           (2*MSt12*s2tb2*(-x10 + x18) + 
     &             2*MSt22*s2tb2*(x10 - x19) + 
     &             (-2*MA02 + 2*c2tb2*c2tt2*MA02 + 
     &                c2tt2*MSb2sum*s2tb2 + 2*MSb22*s2tt2 + 
     &                c2tb2*MSt2sum*s2tt2)*(x18 - x19)) - 
     &          2*c2tt2*hb2*pMT*s2tb*(ctb2*MSb12 - MSb22*stb2)*
     &           (x18 - x19) + 
     &          c2tb*hb2*pMT*s2tb*
     &           (MSt2diff*s2tt2*(-2*x10 + x18 + x19) + 
     &             2*c2tt2*
     &              (-(MSt22*x10) + MSt12*(x10 - x18) + MA02*x18 - 
     &                MA02*x19 + MSt22*x19))) - 
     &       4*Atr*(hb*ht*MUEr*
     &           (2*MSt12*s2tb2*(x10 - x18) + 
     &             2*MSt22*s2tb2*(-x10 + x19) + 
     &             c2tb2*MSt2sum*s2tt2*(-x18 + x19)) + 
     &          c2tb*ht2*pMT*s2tb*
     &           (MSt2diff*s2tt2*(-2*x10 + x18 + x19) + 
     &             2*c2tt2*(MSt12*(x10 - x18) + MSt22*(-x10 + x19))
     &             )))/(MSb2diff*MSt2diff)

        tmp28 = 1/4.D0*(hb*ht*MUEr*
     &      (Abr*(-2*c2tt2*MSb2diff*
     &            (MA02**2 - 3*MA02*MSb12 + 2*MSb12**2 - 
     &              2*MA02*MSt12 - 3*MSb12*MSt12 + MSt12**2)*s2tb2+
     &             MSt2diff*
     &            (2*c2tb2*MSb12*(MA02 + MSb12 - MSt12) + 
     &              MSb2diff*(MA02 - MSt12)*s2tb2)*s2tt2) + 
     &        pMT*s2tt*(c2tt*MSb2diff*
     &            (-2*MA02**2 - 4*MSb12**2 + 6*MSb12*MSt12 + 
     &              MA02*(6*MSb12 + 5*MSt12) + 
     &              MSt12*(-3*MSt12 + MSt22))*s2tb2 - 
     &           MSt2diff*(MSb2diff*MSt12*s2tb2 + 
     &              4*c2tb2*MSb12*(-MSb12 + MSt12)*stt2) + 
     &           MA02*(4*c2tb2*MSb12*MSt2diff*stt2 - 
     &              MSb2diff*s2tb2*(MSt12 - 2*MSt22*stt2)))))/
     &    (MSb12*MSb2diff*MSt12*MSt2diff)

        tmp28 = tmp28 + 1/8.D0*
     &    (s2tb*(ht2*MUEr2*s2tt*
     &           (c2tt*ctb2*(8*MSb12**3 - 4*MA02**2*MSb2diff) - 
     &             2*(ctb2*MSb22*MSt12*
     &                 (3*c2tt*MSt12 - c2tt*MSt22 + MSt2diff) + 
     &                MSb12**2*
     &                 (c2tt*
     &                    (4*ctb2*MSb22 + 3*MSt12 + 4*c2tb*MSt12 - 
     &                      c2tb*MSt22) + c2tb*MSt2diff)) + 
     &             MSb12*MSt12*
     &              (c2tt*(12*ctb2*MSb22 + (3 + 5*c2tb)*MSt12) - 
     &                (1 + 3*c2tb)*(MSt12 - 2*MSt22*stt2)) + 
     &             MA02*(-(c2tt*
     &                   ((5 + 7*c2tb)*MSb12*MSt12 + 
     &                     2*ctb2*
     &                      (6*MSb12**2 + 
     &                       MSb22*(-6*MSb12 - 5*MSt12 + MSt22))))+
     &                  2*ctb2*MSb22*MSt2diff + 
     &                (1 + 3*c2tb)*MSb12*(MSt12 - 2*MSt22*stt2)))+
     &            hb2*(-4*c2tt*MA02**2*MSb2diff*
     &              (2*Abr*c2tt*pMT - Abr**2*s2tt + pMT2*s2tt)*stb2 - 
     &             2*Abr*pMT*
     &              (-(MSb12*(2*c2tb*MSb12 + MSt12 - 3*c2tb*MSt12)*
     &                   MSt2diff*s2tt2) + 
     &                2*(2*c2tt2*MSb2diff*
     &                    (2*MSb12**2 - 3*MSb12*MSt12 + MSt12**2)+
     &                     MSb22*MSt12*MSt2diff*s2tt2)*stb2) + 
     &             s2tt*(Abr**2*
     &                 ((-8*c2tt*MSb12**3 + 
     &                      2*MSb22*MSt12*
     &                       (MSt12 + 3*c2tt*MSt12 - 2*ctt2*MSt22))
     &                     *stb2 + 
     &                   MSb12*MSt12*
     &                    ((-1 - 3*c2tt + c2tb*(3 + 5*c2tt))*
     &                      MSt12 + (2 - 6*c2tb)*ctt2*MSt22 - 
     &                      12*c2tt*MSb22*stb2) + 
     &                   MSb12**2*
     &                    (-2*c2tb*
     &                       (MSt12 + 4*c2tt*MSt12 - 2*ctt2*MSt22)+
     &                        c2tt*(6*MSt12 + 8*MSb22*stb2))) - 
     &                pMT2*((-8*c2tt*MSb12**3 + 
     &                      2*MSb22*MSt12*
     &                       (3*c2tt*MSt12 - c2tt*MSt22 + MSt2diff)
     &                      )*stb2 + 
     &                   MSb12**2*
     &                    (-2*c2tb*MSt2diff + 
     &                      c2tt*
     &                       ((6 - 8*c2tb)*MSt12 + 2*c2tb*MSt22 + 
     &                       8*MSb22*stb2)) + 
     &                   MSb12*MSt12*
     &                    (MSt12 - 3*(c2tb + c2tt)*MSt12 + 
     &                      c2tt*(5*c2tb*MSt12 - 12*MSb22*stb2) + 
     &                      (-2 + 6*c2tb)*MSt22*stt2))) + 
     &             MA02*(2*Abr*pMT*
     &                 ((-1 + 3*c2tb)*MSb12*MSt2diff*s2tt2 + 
     &                   2*
     &                    (6*c2tt2*MSb12*MSb2diff + 
     &                      4*c2tt2*MSb2diff*MSt12 + 
     &                      MSb22*MSt2diff*s2tt2)*stb2) + 
     &                s2tt*
     &                 (Abr**2*
     &                    (-2*MSb22*
     &                       (MSt12 + 5*c2tt*MSt12 - 2*ctt2*MSt22)*
     &                       stb2 + 
     &                      MSb12*
     &                       ((1 + 5*c2tt - c2tb*(3 + 7*c2tt))*
     &                       MSt12 + (-2 + 6*c2tb)*ctt2*MSt22 - 
     &                       12*c2tt*MSb2diff*stb2)) + 
     &                   pMT2*
     &                    (-12*c2tt*MSb12**2*stb2 + 
     &                      2*MSb22*
     &                       (5*c2tt*MSt12 - c2tt*MSt22 + 
     &                      MSt2diff)*stb2 + 
     &                      MSb12*
     &                       ((1 - 3*c2tb - 5*c2tt + 7*c2tb*c2tt)*
     &                       MSt12 + 12*c2tt*MSb22*stb2 - 
     &                       2*MSt22*stt2 + 6*c2tb*MSt22*stt2))))))
     &         + (2*(2*hb*ht*MUEr*
     &             (c2tt*pMT*(MSb12*(2 - 3*s2tb2) + MSb22*s2tb2)*
     &                s2tt + 
     &               Abr*(c2tt2*(-2*MSb12 + MSb2diff)*s2tb2 - 
     &                  2*MSb12*s2tt2)) + 
     &            c2tt*s2tb*
     &             (ht2*(MSb2diff + c2tb*(-2*MSb12 + MSb2diff))*
     &                MUEr2*s2tt + 
     &               hb2*(2*Abr*c2tt*pMT - Abr**2*s2tt + pMT2*s2tt)*
     &                ((-1 + 3*c2tb)*MSb12 + 2*MSb22*stb2))))/x11)/
     &     (MSb12*MSb2diff*MSt12*MSt2diff)

        tmp29 = -(1/4.D0*(hb*ht*MUEr*
     &        (Abr*(2*c2tt2*MSb2diff*
     &              (MA02**2 - 3*MA02*MSb22 + 2*MSb22**2 - 
     &                2*MA02*MSt12 - 3*MSb22*MSt12 + MSt12**2)*
     &              s2tb2 + 
     &             MSt2diff*
     &              (2*c2tb2*MSb22*(MA02 + MSb22 - MSt12) + 
     &                MSb2diff*(-MA02 + MSt12)*s2tb2)*s2tt2) + 
     &          pMT*s2tt*(c2tt*MSb2diff*
     &              (2*MA02**2 - 6*MA02*MSb22 + 4*MSb22**2 - 
     &                5*MA02*MSt12 - 6*MSb22*MSt12 + 3*MSt12**2 - 
     &                MSt12*MSt22)*s2tb2 + 
     &             MSt2diff*
     &              (MSb2diff*MSt12*s2tb2 + 
     &                4*c2tb2*MSb22*(MSb22 - MSt12)*stt2) + 
     &             MA02*(4*c2tb2*MSb22*MSt2diff*stt2 + 
     &                MSb2diff*s2tb2*(MSt12 - 2*MSt22*stt2)))))/
     &      (MSb22*MSb2diff*MSt12*MSt2diff))

        tmp29 = tmp29 - 1/8.D0*
     &    (ht2*MUEr2*s2tb*s2tt*
     &        (-4*c2tt*MA02**2*MSb2diff*stb2 + 
     &          2*MSb12*(4*c2tt*MSb22**2 - 6*c2tt*MSb22*MSt12 + 
     &             MSt12*(3*c2tt*MSt12 - c2tt*MSt22 + MSt2diff))*
     &           stb2 + MSb22*
     &           ((1 - 3*c2tt + c2tb*(-3 + 5*c2tt))*MSt12**2 + 
     &             MSb22*((6 - 8*c2tb)*c2tt*MSt12 + 
     &                2*c2tb*c2tt*MSt22 - 2*c2tb*MSt2diff) - 
     &             8*c2tt*MSb22**2*stb2 + 
     &             2*(-1 + 3*c2tb)*MSt12*MSt22*stt2) + 
     &          MA02*(-2*MSb12*
     &              (c2tt*(6*MSb22 + 4*MSt12 - MSt2diff) + 
     &                MSt2diff)*stb2 + 
     &             MSb22*((-1 + c2tb*(3 - 7*c2tt) + 5*c2tt)*
     &                 MSt12 + 12*c2tt*MSb22*stb2 + 2*MSt22*stt2 - 
     &                6*c2tb*MSt22*stt2))) + 
     &       hb2*s2tb*(Abr**2*
     &           (c2tt*(-(MSb22*MSt12*
     &                   ((5 + 7*c2tb)*MA02 + (6 + 8*c2tb)*MSb22 - 
     &                     (3 + 5*c2tb)*MSt12)) + 
     &                2*ctb2*
     &                 (2*(MA02 - 2*MSb22)*(MA02 - MSb22)*
     &                    MSb2diff + 
     &                   MSb12*(5*MA02 + 6*MSb22)*MSt12 - 
     &                   3*MSb12*MSt12**2)) + 
     &             (2*ctb2*MSb12*(MA02 - MSt12) + 
     &                MSb22*
     &                 (-((1 + 3*c2tb)*MA02) - 2*c2tb*MSb22 + 
     &                   MSt12 + 3*c2tb*MSt12))*
     &              (MSt12 - 2*ctt2*MSt22))*s2tt - 
     &          2*Abr*pMT*(4*c2tt2*ctb2*
     &              (2*MSb22**3 - 3*MA02*MSb22*MSb2diff + 
     &                MSb2diff*(MA02 - MSt12)**2 - 
     &                3*MSb22**2*MSt12 + 
     &                MSb12*MSb22*(-2*MSb22 + 3*MSt12)) + 
     &             (2*ctb2*MSb12*(MA02 - MSt12) + 
     &                MSb22*
     &                 (-((1 + 3*c2tb)*MA02) - 2*c2tb*MSb22 + 
     &                   MSt12 + 3*c2tb*MSt12))*MSt2diff*s2tt2) + 
     &          pMT2*s2tt*(2*ctb2*MSb12*(-MA02 + MSt12)*MSt2diff + 
     &             c2tt*(MSb22*
     &                 (MSt12*
     &                    ((5 + 7*c2tb)*MA02 + 
     &                      (6 + 8*c2tb)*MSb22 - (3 + 5*c2tb)*MSt12
     &                      ) - 2*c2tb*MSb22*MSt22) - 
     &                2*ctb2*
     &                 (4*MSb22**3 + 2*MA02**2*MSb2diff + 
     &                   MA02*
     &                    (-6*MSb22**2 + 
     &                      MSb12*(6*MSb22 + 4*MSt12 - MSt2diff))+
     &                     MSb12*
     &                    (-4*MSb22**2 + 6*MSb22*MSt12 + 
     &                      MSt12*(-2*MSt12 + MSt2diff)))) + 
     &             MSb22*(-(MA02*MSt12) - 3*c2tb*MA02*MSt12 + 
     &                MSt12**2 + 3*c2tb*MSt12**2 + 
     &                2*c2tb*MSb22*MSt2diff + 
     &                2*(1 + 3*c2tb)*(MA02 - MSt12)*MSt22*stt2)))-
     &         (2*(2*(c2tt*hb2*(c2tb*MSb22 + ctb2*MSb2diff)*s2tb*
     &                (-2*Abr*c2tt*pMT + Abr**2*s2tt - pMT2*s2tt) + 
     &               hb*ht*MUEr*
     &                (c2tt*pMT*
     &                   (-(MSb12*s2tb2) + MSb22*(-2 + 3*s2tb2))*
     &                   s2tt + 
     &                  Abr*
     &                   (c2tt2*(2*MSb22 + MSb2diff)*s2tb2 + 
     &                     2*MSb22*s2tt2))) + 
     &            c2tt*ht2*MUEr2*s2tb*s2tt*
     &             ((-1 + 3*c2tb)*MSb22 + 2*MSb12*stb2)))/x13)/
     &     (MSb22*MSb2diff*MSt12*MSt2diff)

        tmp30 = -(c2tb*c2tt*(hb2 + ht2)*s2tb) + 
     &   1/8.D0*(((-2*hb*ht*MUEr*s2tb2*
     &           (Abr*MSb2sum*(MSb12 - MSt12)*(MSb22 - MSt12)*
     &              s2tt + 
     &             2*MSb12*MSb22*(MSb2sum - 2*MSt12)*
     &              (pMT + ctt2*pMT + Atr*s2tt) - 
     &             2*MSt12*(MSb12**2 + MSb22**2 - MSb2sum*MSt12)*
     &              pMT*stt2))/
     &         (MSb12*MSb22*(MSb12 - MSt12)*(-MSb22 + MSt12)) + 
     &        s2tb*(8*c2tt*
     &            (-((c2tb*(hb2 + ht2)*MSt12)/MSt2diff) + 
     &              (ctb2*hb2*pMT2)/(4*MSb22 - 4*MSt12)) - 
     &           (hb2*(2*Abr**2*c2tb*ctt2*(MSb12 - MSt12) + 
     &                c2tt*(MSb12 - MSt12)*pMT2 + 
     &                c2tb*
     &                 (MSb12*pMT2 - MSt12*pMT2 + 4*ctt2*MSb12*MUEr2)+
     &                  2*Abr*c2tb*(MSb12 - MSt12)*pMT*s2tt))/
     &            (MSb12*(MSb12 - MSt12)) + 
     &           2*ctb2*MSt12*
     &            (c2tt*((hb2*pMT2)/(MSb22*(-MSb22 + MSt12)) + 
     &                 (ht2*MUEr2)/(MSb12**2 - MSb12*MSt12)) + 
     &              (hb2*(pMT2 + 2*Abr*(Abr*ctt2 + pMT*s2tt)))/
     &               (MSb22*(MSb22 - MSt12))) + 
     &           hb2*((-2*ctb2*
     &                 (2*Abr**2*ctt2 + pMT2 + 4*ctt2*MUEr2 + 
     &                   2*Abr*pMT*s2tt))/(MSb22 - MSt12) + 
     &              (2*Abr**2*ctt2 + pMT2 + c2tb*c2tt*pMT2 + 
     &                 4*ctt2*MUEr2 + 2*Abr*pMT*s2tt)/(MSb12 - MSt12)+
     &                (MSt12*
     &                 (pMT2 + c2tb*c2tt*pMT2 + 
     &                   2*Abr*(Abr*ctt2 + pMT*s2tt)))/
     &               (MSb12*(-MSb12 + MSt12))) + 
     &           8*c2tb*(-(ctt2*hb2) + ht2*stt2) - 
     &           (4*ht2*stb2*
     &              (2*ctt2*pMT2 + 2*Atr*pMT*s2tt + 2*Atr**2*stt2 + 
     &                MUEr2*stt2))/(MSb22 - MSt12) + 
     &           (2*ht2*(2*MSb12*(MSb12 - MSt12)*MSt12*MUEr2*stb2*
     &                 stt2 + 
     &                ctb2*MSb22*(MSb22 - MSt12)*
     &                 (4*ctt2*MSb12*pMT2 - MSt12*MUEr2 + 
     &                   2*MSb12*
     &                    (2*Atr*pMT*s2tt + 2*Atr**2*stt2 + MUEr2*stt2)
     &                   )))/
     &            (MSb12*MSb22*(MSb12 - MSt12)*(MSb22 - MSt12))))*
     &      (-x10 + x18))

        tmp30 = tmp30 + 1/8.D0*
     &    ((2*(MA02 + MSb12 - MSt12)*
     &         (2*Abr**2*ctt2*hb2*s2tb*stb2 + 
     &           Abr*s2tt*(hb*ht*MUEr*s2tb2 + 2*hb2*pMT*s2tb*stb2) + 
     &           2*(ctb2*ht2*MUEr2*s2tb + hb*ht*pMT*MUEr*s2tb2 + 
     &              hb2*pMT2*s2tb*stb2)*stt2)*x11*
     &         (-(MA02**2*x28) + 
     &           MSt12*(-(MSt12*(x10 - x18 + x28)) + 
     &              MSb12*(x10 + x18 - 2*x21 + x28)) + 
     &           MA02*(MSb12*x28 + MSt12*(x10 - x18 + 2*x28))))/
     &       (MSb12*MSt12) + 
     &      (-x10 + x19)*((-2*hb*ht*MUEr*s2tb2*
     &            ((2*Atr*MSb12*MSb22*(MSb2sum - 2*MSt22) + 
     &                 Abr*MSb2sum*(MSb12 - MSt22)*(MSb22 - MSt22))*
     &               s2tt - 
     &              2*pMT*(ctt2*MSb22*MSt22*(-MSb22 + MSt22) + 
     &                 MSb12**2*(MSb22 - ctt2*MSt22) + 
     &                 MSb12*
     &                  (MSb22**2 + ctt2*MSt22**2 + 
     &                    MSb22*(MSb2sum*stt2 - 2*MSt22*(1 + stt2))
     &                    ))))/
     &          (MSb12*MSb22*(MSb12 - MSt22)*(-MSb22 + MSt22)) + 
     &         s2tb*(8*c2tb*(-(ctt2*ht2) + hb2*stt2) - 
     &            (4*ctb2*ht2*
     &               (2*Atr**2*ctt2 + ctt2*MUEr2 - 2*Atr*pMT*s2tt + 
     &                 2*pMT2*stt2))/(MSb12 - MSt22) + 
     &            (4*(ht2*stb2*
     &                  (2*Atr**2*ctt2*MSb22 + 
     &                    ctt2*(MSb22 - MSt22)*MUEr2 - 
     &                    2*Atr*MSb22*pMT*s2tt + 2*MSb22*pMT2*stt2) + 
     &                 ctb2*hb2*
     &                  (ctt2*(MSb22 - MSt22)*pMT2 + 
     &                    Abr*(-MSb22 + MSt22)*pMT*s2tt + 
     &                    Abr**2*(MSb22 - MSt22)*stt2 + 
     &                    2*MSb22*MUEr2*stt2)))/
     &             (MSb22*(MSb22 - MSt22)) - 
     &            (hb2*(2*c2tb*ctt2*MSt22*pMT2 + 
     &                 2*Abr*pMT*s2tt*(MSt22 - 2*MSb12*stb2) + 
     &                 4*MSb12*stb2*(ctt2*pMT2 + 2*MUEr2*stt2) + 
     &                 Abr**2*
     &                  ((c2tb + c2tt)*MSt22 + 4*MSb12*stb2*stt2)))
     &              /(MSb12*(MSb12 - MSt22)) + 
     &            8*MSt22*((c2tb*c2tt*(hb2 + ht2))/MSt2diff + 
     &               1/8.D0*(4*ctb2*ctt2*ht2*MUEr2 + 
     &                   hb2*
     &                    (2*ctt2*pMT2 + 
     &                      Abr*(Abr + Abr*c2tb*c2tt + 2*c2tb*pMT*s2tt)
     &                      ))/(MSb12*(MSb12 - MSt22))))))

        tmp30 = tmp30 + 1/4.D0*
     &    (-(((MA02 + MSb12 - MSt22)*
     &           (2*ctb2*ctt2*ht2*MUEr2*s2tb + 
     &             2*ctt2*(hb*ht*pMT*MUEr*s2tb2 + 
     &                hb2*pMT2*s2tb*stb2) + 
     &             Abr*(-(hb*ht*MUEr*s2tb2*s2tt) + 
     &                2*hb2*s2tb*stb2*(-(pMT*s2tt) + Abr*stt2)))*x12*
     &           (-(MA02**2*x29) + 
     &             MSt22*(-(MSt22*(x10 - x19 + x29)) + 
     &                MSb12*(x10 + x19 - 2*x21 + x29)) + 
     &             MA02*(MSb12*x29 + MSt22*(x10 - x19 + 2*x29))))/
     &         (MSb12*MSt22)) - 
     &      ((MA02 + MSb22 - MSt12)*
     &         (2*Abr**2*ctb2*ctt2*hb2*s2tb + 
     &           Abr*(2*ctb2*hb2*pMT*s2tb - hb*ht*MUEr*s2tb2)*s2tt + 
     &           2*(ctb2*hb2*pMT2*s2tb - hb*ht*pMT*MUEr*s2tb2 + 
     &              ht2*MUEr2*s2tb*stb2)*stt2)*x13*
     &         (-(MA02**2*x30) + 
     &           MSt12*(-(MSt12*(x10 - x18 + x30)) + 
     &              MSb22*(x10 + x18 - 2*x22 + x30)) + 
     &           MA02*(MSb22*x30 + MSt12*(x10 - x18 + 2*x30))))/
     &       (MSb22*MSt12) + 
     &      ((MA02 + MSb22 - MSt22)*
     &         (-2*ctt2*hb*ht*pMT*MUEr*s2tb2 + 
     &           Abr*hb*ht*MUEr*s2tb2*s2tt + 
     &           2*ctt2*ht2*MUEr2*s2tb*stb2 + 
     &           2*ctb2*hb2*s2tb*
     &            (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*x14*
     &         (-(MA02**2*x31) + 
     &           MSt22*(-(MSt22*(x10 - x19 + x31)) + 
     &              MSb22*(x10 + x19 - 2*x22 + x31)) + 
     &           MA02*(MSb22*x31 + MSt22*(x10 - x19 + 2*x31))))/
     &       (MSb22*MSt22))

        tmp31 = 1/4.D0*(hb*ht*MUEr*
     &     (pMT*(-4*c2tb2*ctt2*MSb12*(MA02 + MSb12 - MSt22)*
     &           MSt2diff + 
     &          MSb2diff*(MSt22**2 + 
     &             c2tt*(2*MA02**2 - 6*MA02*MSb12 + 4*MSb12**2 - 
     &                4*MA02*MSt22 - 6*MSb12*MSt22 + 3*MSt22**2) - 
     &             2*ctt2*(MSt12*MSt22 + MA02*MSt2diff))*s2tb2)*
     &        s2tt + Abr*(2*c2tt2*MSb2diff*
     &           (MA02**2 - 3*MA02*MSb12 + 2*MSb12**2 - 
     &             2*MA02*MSt22 - 3*MSb12*MSt22 + MSt22**2)*s2tb2+
     &            MSt2diff*
     &           (2*c2tb2*MSb12*(MA02 + MSb12 - MSt22) + 
     &             MSb2diff*(MA02 - MSt22)*s2tb2)*s2tt2)))

        tmp31 = tmp31 + 1/8.D0*
     &    (s2tb*(ht2*(c2tt*ctb2*
     &             (-8*MSb12**3 + 4*MA02**2*MSb2diff) + 
     &            2*MSb12**2*
     &             (c2tb*(-2*ctt2*MSt12 + MSt22) + 
     &               c2tt*(4*ctb2*MSb22 + (3 + 4*c2tb)*MSt22)) + 
     &            MSt22*((MSb12 + 3*c2tb*MSb12 - 2*ctb2*MSb22)*
     &                (2*ctt2*MSt12 - MSt22) - 
     &               c2tt*(6*ctb2*MSb22*(2*MSb12 - MSt22) + 
     &                  (3 + 5*c2tb)*MSb12*MSt22)) + 
     &            MA02*((1 + 3*c2tb)*MSb12*
     &                (-2*ctt2*MSt12 + MSt22) + 
     &               c2tt*((5 + 7*c2tb)*MSb12*MSt22 - 
     &                  2*ctb2*
     &                   (6*MSb12*MSb2diff - MSb22*MSt12 + 
     &                     5*MSb22*MSt22)) - 2*ctb2*MSb22*MSt2diff)
     &            )*MUEr2*s2tt + 
     &         hb2*(2*Abr*pMT*
     &             (2*(2*c2tt2*
     &                   (-2*MSb12**3 + MA02**2*MSb2diff - 
     &                     2*MA02*MSb22*MSt22 + MSb22*MSt22**2) + 
     &                  MSb22*(MA02 - MSt22)*MSt2diff*s2tt2)*stb2+
     &                 MSb12*
     &                ((-1 + 3*c2tb)*(MA02 - MSt22)*MSt2diff*
     &                   s2tt2 - 
     &                  4*c2tt2*
     &                   (3*MA02*MSb22 - 2*MA02*MSt22 + 
     &                     3*MSb22*MSt22 + MSt22**2)*stb2) + 
     &               2*MSb12**2*
     &                (c2tb*MSt2diff*s2tt2 + 
     &                  2*c2tt2*(2*MSb22 + 3*(MA02 + MSt22))*stb2))
     &              + s2tt*
     &             (pMT2*(-2*MSb22*(MA02 - MSt22)*MSt2diff*stb2 + 
     &                  MSb12*
     &                   (2*(-1 + 3*c2tb)*ctt2*MSt12*
     &                      (MA02 - MSt22) - 
     &                     (1 + 3*c2tt)*MSt22**2 + 
     &                     c2tb*
     &                      (-((3 + 7*c2tt)*MA02*MSt22) + 
     &                       (3 + 5*c2tt)*MSt22**2 + 
     &                       2*MSb12*
     &                       (c2tt*(MSt12 - 4*MSt22) - MSt2diff))+
     &                       MA02*(MSt22 - 12*c2tt*MSb22*stb2)) + 
     &                  c2tt*
     &                   ((-8*MSb12**3 + 
     &                       4*MSb12**2*(3*MA02 + 2*MSb22) + 
     &                       4*MA02**2*MSb2diff + 
     &                       2*MSb22*
     &                       (MA02*MSt12 + 
     &                       MSt22*(-6*MSb12 - MSt12 + 3*MSt22)))*
     &                      stb2 + 
     &                     MSt22*
     &                      (6*MSb12**2 + 
     &                       5*MA02*(MSb12 - 2*MSb22*stb2)))) + 
     &               Abr**2*
     &                (c2tt*
     &                   (3*MSb12*MSt22**2 + 8*MSb12**3*stb2 - 
     &                     4*MA02**2*MSb2diff*stb2 - 
     &                     2*MSb12**2*(3*MSt22 + 4*MSb22*stb2)) + 
     &                  MSb22*MSt22*stb2*
     &                   ((2 - 6*c2tt)*MSt22 - 4*MSt12*stt2) + 
     &                  MA02*
     &                   (-12*c2tt*MSb12**2*stb2 + 
     &                     MSb12*
     &                      (MSt22 - 5*c2tt*MSt22 + 
     &                       12*c2tt*MSb22*stb2 - 2*MSt12*stt2) + 
     &                     2*MSb22*stb2*
     &                      ((-1 + 5*c2tt)*MSt22 + 2*MSt12*stt2))+
     &                    MSb12*
     &                   (-MSt22**2 + 
     &                     2*MSt22*
     &                      (6*c2tt*MSb22*stb2 + MSt12*stt2) + 
     &                     c2tb*
     &                      (MSt22*
     &                       ((-3 + 7*c2tt)*MA02 + 
     &                       (-2 + 8*c2tt)*MSb12 + 
     &                       (3 - 5*c2tt)*MSt22) + 
     &                       2*MSt12*(3*MA02 + 2*MSb12 - 3*MSt22)*
     &                       stt2)))))) + 
     &      (2*(2*hb*ht*MUEr*
     &            (c2tt*pMT*(-(MSb22*s2tb2) + MSb12*(-2 + 3*s2tb2))*
     &               s2tt + 
     &              Abr*(c2tt2*(2*MSb12 - MSb2diff)*s2tb2 + 
     &                 2*MSb12*s2tt2)) + 
     &           c2tt*s2tb*
     &            (ht2*(MSb12 + 3*c2tb*MSb12 - 2*ctb2*MSb22)*MUEr2*
     &               s2tt + 
     &              hb2*(2*Abr*c2tt*pMT - Abr**2*s2tt + pMT2*s2tt)*
     &               (MSb12 - 3*c2tb*MSb12 - 2*MSb22*stb2))))/x12)

        tmp32 = 1/8.D0*(2*hb*ht*MUEr*
     &      (pMT*(4*c2tb2*ctt2*MSb22*(MA02 + MSb22 - MSt22)*
     &            MSt2diff + 
     &           MSb2diff*(MSt22**2 + 
     &              c2tt*(2*MA02**2 - 6*MA02*MSb22 + 4*MSb22**2 - 
     &                 4*MA02*MSt22 - 6*MSb22*MSt22 + 3*MSt22**2)-
     &                2*ctt2*(MSt12*MSt22 + MA02*MSt2diff))*s2tb2)
     &          *s2tt + Abr*
     &         (2*c2tt2*MSb2diff*
     &            (MA02**2 - 3*MA02*MSb22 + 2*MSb22**2 - 
     &              2*MA02*MSt22 - 3*MSb22*MSt22 + MSt22**2)*s2tb2+
     &             MSt2diff*
     &            (-2*c2tb2*MSb22*(MA02 + MSb22 - MSt22) + 
     &              MSb2diff*(MA02 - MSt22)*s2tb2)*s2tt2)) - 
     &     ht2*MUEr2*s2tb*s2tt*
     &      (4*c2tt*MA02**2*MSb2diff*stb2 - 
     &        2*MSb12*(4*c2tt*MSb22**2 - 6*c2tt*MSb22*MSt22 + 
     &           MSt22*(-2*ctt2*MSt12 + MSt22 + 3*c2tt*MSt22))*stb2
     &          + MSb22*(MSt22*
     &            (2*(-1 + 3*c2tb)*ctt2*MSt12 + 
     &              (1 + 3*c2tt - c2tb*(3 + 5*c2tt))*MSt22) + 
     &           2*MSb22*(-3*c2tt*MSt22 + 
     &              c2tb*(-(c2tt*MSt12) + 4*c2tt*MSt22 + MSt2diff))
     &             + 8*c2tt*MSb22**2*stb2) + 
     &        MA02*(2*(1 + 5*c2tt)*MSb12*MSt22*stb2 + 
     &           2*ctt2*MSt12*
     &            (MSb22 - 3*c2tb*MSb22 - 2*MSb12*stb2) + 
     &           MSb22*((-1 + 3*c2tb - 5*c2tt + 7*c2tb*c2tt)*
     &               MSt22 - 12*c2tt*MSb2diff*stb2))) - 
     &     hb2*s2tb*((2*ctb2*MSb12*(MA02 - MSt22)*
     &            (-2*ctt2*MSt12 + MSt22) + 
     &           c2tt*(MSb22*
     &               (2*c2tb*MSb22*MSt12 - 
     &                 (5*MA02 + 7*c2tb*MA02 + 6*MSb22 + 
     &                    8*c2tb*MSb22)*MSt22 + 
     &                 (3 + 5*c2tb)*MSt22**2) + 
     &              2*ctb2*
     &               (2*(MA02 - 2*MSb22)*(MA02 - MSb22)*MSb2diff + 
     &                 MSb12*(5*MA02 + 6*MSb22)*MSt22 - 
     &                 3*MSb12*MSt22**2)) + 
     &           MSb22*(-((1 + 3*c2tb)*(MA02 - MSt22)*
     &                 (-2*ctt2*MSt12 + MSt22)) - 
     &              2*c2tb*MSb22*MSt2diff))*pMT2*s2tt + 
     &        2*Abr*pMT*(2*ctb2*
     &            (2*c2tt2*
     &               (2*MSb22**3 - 3*MA02*MSb22*MSb2diff + 
     &                 MSb2diff*(MA02 - MSt22)**2 + 
     &                 MSb12*MSb22*(-2*MSb22 + 3*MSt22)) + 
     &              MSb12*(-MA02 + MSt22)*MSt2diff*s2tt2) + 
     &           MSb22*((1 + 3*c2tb)*(MA02 - MSt22)*MSt2diff*
     &               s2tt2 + 
     &              2*MSb22*
     &               (-6*c2tt2*ctb2*MSt22 + c2tb*MSt2diff*s2tt2)))+
     &          Abr**2*s2tt*
     &         (2*ctb2*MSb12*(MA02 - MSt22)*MSt2diff + 
     &           c2tt*(MSb22*MSt22*
     &               ((5 + 7*c2tb)*MA02 + (6 + 8*c2tb)*MSb22 - 
     &                 (3 + 5*c2tb)*MSt22) - 
     &              2*ctb2*
     &               (4*MSb22**3 + 2*MA02**2*MSb2diff - 
     &                 MSb12*
     &                  (4*MSb22**2 - 6*MSb22*MSt22 + 
     &                    MSt22*(2*MSt22 + MSt2diff)) + 
     &                 MA02*
     &                  (-6*MSb22**2 + 
     &                    MSb12*(6*MSb22 + 4*MSt22 + MSt2diff))))-
     &             MSb22*(MA02 + 3*c2tb*MA02 + 2*c2tb*MSb22 - 
     &              MSt22 - 3*c2tb*MSt22)*(MSt22 - 2*MSt12*stt2))))

        tmp32 = tmp32 + 1/4.D0*
     &    (2*(c2tt*hb2*(c2tb*MSb22 + ctb2*MSb2diff)*s2tb*
     &           (2*Abr*c2tt*pMT - Abr**2*s2tt + pMT2*s2tt) + 
     &          hb*ht*MUEr*(c2tt*pMT*
     &              (MSb22*(2 - 3*s2tb2) + MSb12*s2tb2)*s2tt - 
     &             Abr*(c2tt2*(2*MSb22 + MSb2diff)*s2tb2 + 
     &                2*MSb22*s2tt2))) + 
     &       c2tt*ht2*MUEr2*s2tb*s2tt*
     &        (MSb22 - 3*c2tb*MSb22 - 2*MSb12*stb2))/x14

        s2tts2tbF5 = s2tt*tmp30 + tmp26*(-x10 + x21) + 
     &   tmp27*(-x10 + x22) + tmp28*x28 + tmp29*x30 + 
     &   (Atr*(c2tb*ht2*pMT*s2tb + hb*ht*MUEr*s2tb2)*
     &       (2*x10 - x21 - x22)*(x21 - x22) + 
     &      ((tmp31*x29)/(MSb12*MSt22) + 
     &         MSb12*(2*Atr*
     &             (hb*ht*MUEr*s2tt2*(-x33 + x35) + 
     &               c2tt2*ht2*pMT*s2tb*stb2*(x33 - x34 - x35 + x36)
     &               ) + c2tt*s2tt*
     &             (2*hb*ht*pMT*MUEr*(-x33 + x35) - 
     &               s2tb*(ctb2*hb2*MUEr2 + 
     &                  ht2*(-Atr**2 + pMT2)*stb2)*
     &                (x33 - x34 - x35 + x36))) + 
     &         1/2.D0*(2*tmp32*x31 + 
     &             MSb22*MSt22*
     &              (-(Atr**2*ht2*
     &                   (2*c2tt*ctb2*MSb22 + c2tb*MSt2diff - 
     &                     c2tb*c2tt*MSt2sum)*s2tb*s2tt*
     &                   (x33 - x34 - x35 + x36)) - 
     &                2*Atr*
     &                 (c2tb*ht2*MSt2diff*pMT*s2tb*s2tt2*
     &                    (x33 - x34 + x35 - x36) + 
     &                   hb*ht*MUEr*
     &                    (-2*MSb22*s2tt2*x34 + 
     &                      2*MSt12*s2tb2*(-x33 + x34) + 
     &                      2*MSt22*s2tb2*x35 + 
     &                      c2tb2*MSt2sum*s2tt2*
     &                       (-x33 + x34 + x35 - x36) - 
     &                      2*MSt22*s2tb2*x36 + 2*MSb22*s2tt2*x36)+
     &                     c2tt2*
     &                    (2*ctb2*ht2*MSb22*pMT*s2tb*
     &                       (x33 - x34 - x35 + x36) + 
     &                      hb*ht*MSb2sum*MUEr*s2tb2*
     &                       (x33 - x34 - x35 + x36) - 
     &                      2*c2tb*ht2*pMT*s2tb*
     &                       (MSt12*x33 - MSt12*x34 - MSt22*x35 + 
     &                       MSt22*x36))) + 
     &                s2tt*
     &                 (2*c2tb2*hb*ht*(MSt2diff + c2tt*MSt2sum)*pMT*
     &                    MUEr*(x33 - x34 - x35 + x36) - 
     &                   c2tb*MSt2diff*(ht2*pMT2 - hb2*MUEr2)*s2tb*
     &                    (x33 - x34 - x35 + x36) + 
     &                   c2tt*
     &                    (2*ctb2*ht2*MSb22*pMT2*s2tb*
     &                       (x33 - x34 - x35 + x36) - 
     &                      c2tb*MSt2sum*(ht2*pMT2 - hb2*MUEr2)*s2tb*
     &                       (x33 - x34 - x35 + x36) + 
     &                      2*
     &                       (hb2*MSb22*MUEr2*s2tb*stb2*
     &                       (x33 - x34 - x35 + x36) + 
     &                       hb*ht*pMT*MUEr*
     &                       (2*MSb22*(x34 - x36) + 
     &                       MSb2sum*s2tb2*(x33 - x34 - x35 + x36))
     &                       )))))/(MSb22*MSt22))/MSt2diff)/
     &    MSb2diff

        s2tts2tbF5 = s2tts2tbF5 + 
     &   (-x10 + x9)*((hb*ht*MUEr*(c2tt*pMT*s2tt - Abr*s2tt2)*
     &         (x18 - x19))/MSt2diff + 
     &      (Abr*(c2tb*hb2*pMT*s2tb - hb*ht*MUEr*s2tb2)*(x21 - x22))/
     &       MSb2diff + 1/4.D0*
     &       ((2*(hb*ht*MSb2sum*MUEr*s2tb2*
     &                (c2tt*pMT*s2tt - Abr*s2tt2) + 
     &               hb2*s2tb*s2tt*
     &                (Abr**2*c2tt - c2tt*pMT2 + 2*Abr*pMT*s2tt)*
     &                (ctb2*MSb12 - MSb22*stb2)))/MSb22 + 
     &          s2tt*((MA02 + MSb12 - MSt12)**2*
     &              (Abr*(hb*ht*MUEr*s2tb2*s2tt + 
     &                   2*hb2*s2tb*(Abr*ctt2 + pMT*s2tt)*stb2) + 
     &                2*(ctb2*ht2*MUEr2*s2tb + hb*ht*pMT*MUEr*s2tb2 + 
     &                   hb2*pMT2*s2tb*stb2)*stt2)*x11 + 
     &             (MA02 + MSb12 - MSt22)**2*
     &              (hb*ht*MUEr*s2tb2*(-2*ctt2*pMT + Abr*s2tt) - 
     &                2*s2tb*
     &                 (ctb2*ctt2*ht2*MUEr2 + 
     &                   hb2*stb2*
     &                    (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2))))*
     &              x12))/MSb12 + 
     &      1/4.D0*(s2tt*((2*c2tt*ht2*MUEr2*s2tb*
     &               (ctb2*MSb22 - MSb12*stb2))/MSb12 + 
     &            (MA02 + MSb22 - MSt12)**2*
     &             (hb*ht*MUEr*s2tb2*(Abr*s2tt + 2*pMT*stt2) - 
     &               2*s2tb*
     &                (Abr*ctb2*hb2*(Abr*ctt2 + pMT*s2tt) + 
     &                  (ctb2*hb2*pMT2 + ht2*MUEr2*stb2)*stt2))*x13+
     &              (MA02 + MSb22 - MSt22)**2*
     &             (-2*ctt2*hb*ht*pMT*MUEr*s2tb2 + 
     &               Abr*hb*ht*MUEr*s2tb2*s2tt + 
     &               2*ctt2*ht2*MUEr2*s2tb*stb2 + 
     &               2*ctb2*hb2*s2tb*
     &                (ctt2*pMT2 + Abr*(-(pMT*s2tt) + Abr*stt2)))*x14))
     &         /MSb22)

        F6 = 1/8.D0*(((-((MA02 + MSb12 - MSt12)*
     &              (2*Abr*ctt2*hb2*s2tb + 2*ctb2*hb*ht*MUEr*s2tt + 
     &                hb2*pMT*s2tb*s2tt)) + (hb2*s2tb*s2tt)/(pMT*x11)
     &           )*x28 + ((MA02 + MSb22 - MSt12)*
     &            (2*Abr*ctt2*hb2*s2tb + hb2*pMT*s2tb*s2tt - 
     &              2*hb*ht*MUEr*s2tt*stb2) - 
     &           (hb2*s2tb*s2tt)/(pMT*x13))*x30)/MSt12 + 
     &     (((MA02 + MSb12 - MSt22)*
     &            (2*ctb2*hb*ht*MUEr*s2tt + 
     &              hb2*s2tb*(pMT*s2tt - 2*Abr*stt2)) - 
     &           (hb2*s2tb*s2tt)/(pMT*x12))*x29 + 
     &        ((MA02 + MSb22 - MSt22)*
     &            (-(hb2*pMT*s2tb*s2tt) + 2*hb*ht*MUEr*s2tt*stb2 + 
     &              2*Abr*hb2*s2tb*stt2) + (hb2*s2tb*s2tt)/(pMT*x14))
     &          *x31)/MSt22 + 
     &     (2*(ht2*s2tb*(-(pMT2*s2tt*x33) - 2*Atr*pMT*stt2*x33 - 
     &             MSb22*s2tt*x34 + pMT2*s2tt*x34 + 
     &             2*Atr*pMT*stt2*x34 + MSt12*s2tt*(-x33 + x34) + 
     &             MSb12*s2tt*(x33 - x35) - 2*Atr*ctt2*pMT*x35 + 
     &             MSt22*s2tt*x35 + pMT2*s2tt*x35 + 
     &             2*Atr*ctt2*pMT*x36 + MSb22*s2tt*x36 - 
     &             MSt22*s2tt*x36 - pMT2*s2tt*x36) + 
     &          2*hb*ht*pMT*MUEr*s2tt*
     &           (stb2*(-x33 + x35) + ctb2*(-x34 + x36))))/pMT + 
     &     8*hb*ht*MUEr*s2tt*
     &      (x10*(x18 - x19) - x18*x24 + x19*x24 + x6 - x7) + 
     &     (((s2tt*((2*ht2*(MSt12 + pMT2) + 
     &                    hb2*(-MA02 + MSb12 + MSt12 + pMT2))*s2tb+
     &                   2*hb*ht*pMT*MUEr*(1 + stb2)) + 
     &              2*pMT*s2tb*(Abr*ctt2*hb2 + 2*Atr*ht2*stt2))*
     &            (-x10 + x18) + 
     &           (-(s2tt*((2*ht2*(MSt22 + pMT2) + 
     &                      hb2*(-MA02 + MSb12 + MSt22 + pMT2))*s2tb
     &                     + 2*hb*ht*pMT*MUEr*(1 + stb2))) + 
     &              2*pMT*s2tb*(2*Atr*ctt2*ht2 + Abr*hb2*stt2))*
     &            (-x10 + x19))*(-x10 + x21) + 
     &        ht2*s2tb*(2*Atr*pMT - MSt2diff*s2tt)*
     &         (2*x10 - x21 - x22)*(x21 - x22) + 
     &        (((2*(1 + ctb2)*hb*ht*pMT*MUEr + 
     &                 (-2*ht2*(MSt12 + pMT2) - 
     &                    hb2*(-MA02 + MSb22 + MSt12 + pMT2))*s2tb)*
     &               s2tt - 2*pMT*s2tb*(Abr*ctt2*hb2 + 2*Atr*ht2*stt2)
     &              )*(-x10 + x18) + 
     &           (pMT*(-4*Atr*ctt2*ht2*s2tb - 
     &                 2*(1 + ctb2)*hb*ht*MUEr*s2tt) + 
     &              s2tb*(2*ht2*(MSt22 + pMT2)*s2tt + 
     &                 hb2*
     &                  ((-MA02 + MSb22 + MSt22 + pMT2)*s2tt - 
     &                    2*Abr*pMT*stt2)))*(-x10 + x19))*
     &         (-x10 + x22) + 
     &        ((2*hb*ht*pMT*MUEr + hb2*MSb2diff*s2tb)*s2tt*
     &            (x18 - x19) + 
     &           hb2*s2tb*(2*Abr*pMT - MSt2diff*s2tt)*(-x21 + x22))*
     &         (-x10 + x9))/pMT)

* COUNTER TERMS:

	dMUErF2t = 0
	dMUErF3t = 0
	dMUErF3b = 0
	s2ttdAtrF2t = 0
	s2ttdAtrF3t = 0
	s2tbAbdAbF3b = 0
	WF11 = 1
	WF12 = 1
	WF22 = 1

	if( onshell ) then
	  F1tol = x18 + x19 - 2*x20
	  F2tol = x18 - x19
	  F3tol = 2 + MSt2sum/MSt2diff*F2tol
	  F2bol = x21 - x22
	  F3bol = 2 + MSb2sum/MSb2diff*F2bol

	  dvev2vev2 = -(1/2.D0*
     &     (ht2*Nc*(pMT2 + 2*pMT2*x10 - 2*pMT2*x20 + 
     &          4*ctb2*ctt2*fW(MSt12,MSb12) + 
     &          4*ctt2*stb2*fW(MSt12,MSb22) + 
     &          4*ctb2*stt2*fW(MSt22,MSb12) + 
     &          4*stb2*stt2*fW(MSt22,MSb22)))/pMT2)

	  dpMTpMT = Pitt(pMT2)/pMT
	  dMSt12 = Pi11t(MSt12)
	  dMSt22 = Pi22t(MSt22)
	  pi1 = Pi12t(MSt12)
	  pi2 = Pi12t(MSt22)
	  ds2tt = -c2tt*(pi1 + pi2)/MSt2diff

	  dhtht = dpMTpMT - .5D0*dvev2vev2

	  dMSb12 = Pi11b(MSb12)
	  dMSb22 = Pi22b(MSb22)
	  pi1 = Pi12b(MSb12)
	  pi2 = Pi12b(MSb22)
	  ds2tb = -c2tb*(pi1 + pi2)/MSb2diff

	  s2tbdhbhb =
     &      (s2tb*(dMSb22 - dMSb12) - c2tb*(pi2 + pi1))/MSb2diff -
     &      .5D0*s2tb*dvev2vev2

	  s2ttdAtrAtr = ds2tt + s2tt*((dMSt22 - dMSt12)/MSt2diff - dpMTpMT)
	  s2ttdAtrF2t = s2ttdAtrAtr*F2tol
	  s2ttdAtrF3t = s2ttdAtrAtr*F3tol

	  dAbthb = -.5D0*Abr*(PiAAb(MSb12) - PiAAb(MSb22) +
     &      dMSb12 + Pi22b(MSb12) - Pi11b(MSb22) - dMSb22)/MSb2diff -
     &      1/sqrt2*(Lambhb(MSb12) + Lambhb(MSb22))
	  s2tbdAb = s2tb*dAbthb - Abr*s2tbdhbhb
	  s2tbAbdAbF3b = s2tbdAb*F3bol

	  k = s2tt*(dMSt12/MSt12 - dMSt22/MSt22)

	  F1t = F1t + dMSt12/MSt12 + dMSt22/MSt22 - 4*dpMTpMT +
     &      2*(dhtht + dpMTpMT)*F1tol

	  s2ttF2t = s2ttF2t + k +
     &      (ds2tt + s2tt*(2*dhtht + dpMTpMT))*F2tol

	  s2tt2F3t = s2tt2F3t + s2tt*(
     &      k*(2*MSt12*MSt22/MSt2diff*F2tol + MSt2sum)/MSt2diff +
     &      2*(s2tt*dhtht + ds2tt)*F3tol )

	  k = s2tb*(dMSb12/MSb12 - dMSb22/MSb22)

	  s2tb2F3b = s2tb2F3b + s2tb*(
     &      k*(2*MSb12*MSb22/MSb2diff*F2bol + MSb2sum)/MSb2diff +
     &      2*(s2tbdhbhb + ds2tb)*F3bol )

#ifdef COMPLETE_REN
	  dMUErMUEr = 0
	  dMUErF2t = dMUErMUEr*F2tol
	  dMUErF3t = dMUErMUEr*F3tol
	  dMUErF3b = dMUErMUEr*F3bol

	  M11ol = .5D0*(ht2*s2tt2*MUEr2*F3tol + hb2*s2tb2*Abr**2*F3bol)

	  M12ol = -ht2*s2tt*MUEr*(pMT*F2tol + .5D0*s2tt*Atr*F3tol) -
     &      .5D0*hb2*s2tb2*MUEr*Abr*F3bol

	  M22ol =
     &      ht2*(2*pMT2*F1tol + 2*s2tt*pMT*Atr*F2tol +
     &        .5D0*s2tt2*Atr**2*F3tol) +
     &      .5D0*hb2*s2tb2*MUEr2*F3bol

	  dZ1 = 0
	  dZ2 = 0
	  WF11 = 1 + 2*dZ1*M11ol
	  WF12 = 1 + (dZ1 + dZ2)*M12ol
	  WF22 = 1 + 2*dZ2*M22ol
#endif
	endif

	k = Nc/(16*pi**2)**2

	dPhi(se11) = k*WF11*(
     &    ht2 * .5D0*MUEr2*(s2tt2F3t + 2*s2tt2*dMUErF3t) -
     &    ht*hb * MUEr*(2*s2ttF4t + Abr*s2tts2tbF5) +
     &    hb2* (2*F1b + 2*Abr*s2tbF2b +
     &      Abr*(.5D0*Abr*s2tb2F3b + s2tb*s2tbAbdAbF3b)) )

	dPhi(se12) = -k*WF12*(
     &    ht2 * MUEr*(pMT*(s2ttF2t + s2tt*dMUErF2t)
     &     + .5D0*Atr*(s2tt2F3t + s2tt*(s2ttdAtrF3t + s2tt*dMUErF3t))) -
     &    ht*hb * (Atr*s2ttF4t + pMT*Abr*s2tbF4b
     &     + .5D0*(Atr*Abr + MUEr2)*s2tts2tbF5 + 2*pMT*F6) +
     &    hb2 * MUEr*(s2tbF2b + .5D0*(Abr*s2tb2F3b +
     &      s2tb*(s2tbAbdAbF3b + s2tb*Abr*dMUErF3b))) )

	dPhi(se22) = k*WF22*(
     &    ht2 * (2*pMT2*F1t
     &     + 2*Atr*pMT*(s2ttF2t + s2ttdAtrF2t)
     &     + .5D0*Atr**2*(s2tt2F3t + 2*s2tt*s2ttdAtrF3t)) -
     &    ht*hb * MUEr*(2*pMT*s2tbF4b + Atr*s2tts2tbF5) +
     &    hb2 * .5D0*MUEr2*(s2tb2F3b + 2*s2tb2*dMUErF3b) )

	if( debuglevel .lt. 6 ) return

	DHIGGS "atab: F1t        =", F1t	ENDL
	DHIGGS "atab: s2ttF2t    =", s2ttF2t	ENDL
	DHIGGS "atab: s2tt2F3t   =", s2tt2F3t	ENDL
	DHIGGS "atab: s2ttF4t    =", s2ttF4t	ENDL
	DHIGGS "atab: s2tts2tbF5 =", s2tts2tbF5	ENDL
	DHIGGS "atab: s2tbF2b    =", s2tbF2b	ENDL
	DHIGGS "atab: s2tb2F3b   =", s2tb2F3b	ENDL
	DHIGGS "atab: s2tbF4b    =", s2tbF4b	ENDL
	DHIGGS "atab: s2tts2tbF5 =", s2tts2tbF5	ENDL
	DHIGGS "atab: F6         =", F6		ENDL
	DHIGGS "atab:   se11     =", dPhi(se11)	ENDL
	DHIGGS "atab:   se12     =", dPhi(se12)	ENDL
	DHIGGS "atab:   se22     =", dPhi(se22)	ENDL
	end

