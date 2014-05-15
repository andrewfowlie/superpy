* TLatat.F
* Two-loop O(a_t^2) corrections to the CP-even Higgs mass matrix.
* Routine written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
* Based on A. Brignole, G. Degrassi, P. Slavich and F. Zwirner,
* hep-ph/0112177.
* this file is part of FeynHiggs
* last modified 10 Dec 13 th

* Revision history:
* 14/07/2004: rewritten completely by TH for FeynHiggs
* 22/10/2002: first written

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine TLatat(dPhi)
	implicit none
	RealType dPhi(3)

#include "TLps.h"

	RealType F1, s2ttF2, s2ttF3, F1o, F2o, F3o
	RealType dMUErF2, dMUErF3
	RealType s2ttAtrdAtrF2, s2ttAtrdAtrF3, s2ttdAtr
	RealType dvev2vev2, dpMTpMT, dMSt12, dMSt22, ds2tt, k
	RealType MbSL2

	RealType x20, x19, x21, x22, x6, x13, x14, x1, x3, x31
	RealType x32, x17, x4, x5, x18, x12, x9, x11, x7, x8
	RealType x10, x44, x41, x46, x38, x39, x40, x42, x43
	RealType x45, x47, x49, x48, x50, x51, x52, x53, x54
	RealType x55, x56, x57, x58, x59, x34, x33, x35, x36
	RealType x37, x15, x16, x25, x28, x26, x29, x23, x24
	RealType x27, x30
	RealType tmp1, tmp10, tmp11, tmp12, tmp13, tmp14
	RealType tmp15, tmp16, tmp17, tmp18, tmp19, tmp2
	RealType tmp20, tmp21, tmp22, tmp23, tmp3, tmp4, tmp5
	RealType tmp6, tmp8, tmp9
	RealType lc10, lc11, lc12, lc13, lc14, lc15, lc16
	RealType lc17, lc18, lc19, lc20, lc21, lc3, lc4, lc5, lc6
	RealType lc7, lc8, lc9

	MbSL2 = BSS2(2) + BSS2(1) - DSf(1,4)

	x20 = InvTri(MA02,MbSL2,MSt12)

	x19 = InvTri(MA02,MbSL2,MSt22)

	x21 = InvTri(MA02,MSt12,MSt22)

	x22 = InvTri(MA02,MSt22,MSt12)

	x6 = InvTri(MbSL2,pMT2,MUEr2)

	x13 = InvTri(MSt12,pMT2,MUEr2)

	x14 = InvTri(MSt22,pMT2,MUEr2)

	x1 = Li2omx(MbSL2/MSt12)

	x3 = Li2omx(MbSL2/MSt22)

	x31 = Li2omx(MSt22/MSt12)

	x32 = Li2omx(MSt12/MSt22)

	x17 = Li2omx(MA02/pMT2)

	x4 = Li2omx(MUEr2/MSt12)

	x5 = Li2omx(MUEr2/MSt22)

	x18 = Log(MA02)

	x12 = Log(MbSL2)

	x9 = Log(MSt12)

	x11 = Log(MSt22)

	x7 = Log(pMT2)

	x8 = Log(MUEr2)

	x10 = Log(Q)

	x44 = myAA(MA02,Q)

	x41 = myAA(MbSL2,Q)

	x46 = myAA(MSt12,Q)

	x38 = myAA(MSt22,Q)

	x39 = myAA(pMT2,Q)

	x40 = myAA(MUEr2,Q)

	x42 = myB0(MSt12,0.D0,MUEr2,Q)

	x43 = myB0(MSt12,MbSL2,0.D0,Q)

	x45 = myB0(MSt12,MbSL2,MA02,Q)

	x47 = myB0(MSt12,MSt12,0.D0,Q)

	x49 = myB0(MSt12,MSt12,MA02,Q)

	x48 = myB0(MSt12,MSt22,0.D0,Q)

	x50 = myB0(MSt12,MSt22,MA02,Q)

	x51 = myB0(MSt12,pMT2,MUEr2,Q)

	x52 = myB0(MSt22,0.D0,MUEr2,Q)

	x53 = myB0(MSt22,MbSL2,0.D0,Q)

	x54 = myB0(MSt22,MbSL2,MA02,Q)

	x55 = myB0(MSt22,MSt12,0.D0,Q)

	x56 = myB0(MSt22,MSt12,MA02,Q)

	x57 = myB0(MSt22,MSt22,0.D0,Q)

	x58 = myB0(MSt22,MSt22,MA02,Q)

	x59 = myB0(MSt22,pMT2,MUEr2,Q)

	x34 = myB0(pMT2,0.D0,MA02,Q)

	x33 = myB0(pMT2,pMT2,MA02,Q)

	x35 = myB0(pMT2,MUEr2,MbSL2,Q)

	x36 = myB0(pMT2,MUEr2,MSt12,Q)

	x37 = myB0(pMT2,MUEr2,MSt22,Q)

	x15 = Phi(MA02,MbSL2,MSt12)

	x16 = Phi(MA02,MbSL2,MSt22)

	x25 = Phi(MA02,MSt12,MSt12)

	x28 = Phi(MA02,MSt12,MSt22)

	x26 = Phi(MA02,MSt22,MSt12)

	x29 = Phi(MA02,MSt22,MSt22)

	x23 = Phi(MA02,pMT2,pMT2)

	x24 = Phi(MbSL2,MUEr2,pMT2)

	x27 = Phi(MUEr2,pMT2,MSt12)

	x30 = Phi(MUEr2,pMT2,MSt22)

	tmp1 = 2*MSt12**2 + 5*MSt12*MSt22 + MSt22**2

	tmp2 = MSt12**2 + 5*MSt12*MSt22 + 2*MSt22**2

	tmp3 = 2*stt2 - (s2tt*Ytr)/pMT

	tmp4 = pSB2*Xtr**2 + pCB2*Ytr**2

        tmp8 = (MSt12 - pMT2)**3 + 2*MSt12*(-MSt12 + pMT2)*MUEr2 + 
     &   (MSt12 + pMT2)*MUEr2**2

	tmp9 = (MA02 - MbSL2 + MSt12)*(x10 - x12) - 2*MA02*x15

	tmp10 = (MA02 - MbSL2 + MSt22)*(x10 - x12) - 2*MA02*x16

	tmp11 = pMT2*stt2 + Xtr*(-(pMT*s2tt) + ctt2*Xtr)

	tmp13 = ctt2*pMT2 + Xtr*(pMT*s2tt + stt2*Xtr)

	tmp14 = pMT2*stt2 + Ytr*(-(pMT*s2tt) + ctt2*Ytr)

	tmp16 = ctt2*pMT2 + Ytr*(pMT*s2tt + stt2*Ytr)

	tmp19 = ctt2 + (pSB2*tmp11)/(MbSL2 - MSt22)

	tmp20 = stt2 + (pSB2*tmp13)/(MbSL2 - MSt12)

        tmp22 = (1 + c2tt2)*pCB2*
     &   ((MA02**2*MSt12 + MSt2diff**2*MSt2sum - MA02*tmp1)*
     &      (-x10 + x11) + 
     &     (MA02**2*MSt22 + MSt2diff**2*MSt2sum - MA02*tmp2)*
     &      (-x10 + x9))*Ytr**2

	tmp23 = tmp3*(-x10 + x11) + (-x10 + x9)*(2*ctt2 + (s2tt*Ytr)/pMT)

        F1 = 3 - s2tt2 + 2*(stt2*x4 + ctt2*x5) + 
     &   (-3 + Nc)*(x10 - x7)**2 + 
     &   (-x10 + x7)*(-(MSt2sum/pMT2) - pMT2/MSt12 - pMT2/MSt22 + 
     &      (4*pMT2*(-MbSL2 + 3*(4*pMT2 + MUEr2)) + 
     &         MA02*(MbSL2 - 3*MUEr2 - 4*pMT2*(2 + pSB2)))/
     &       (pMT2*(-MA02 + 4*pMT2)) + x10 - x8) + 
     &   MUEr2*(-((1 + 2*ctt2)/MSt22) + 3/pMT2 - (1 + 2*stt2)/MSt12)*
     &    (-x10 + x8) + ((MbSL2 - pMT2 + MUEr2)*x24 + 
     &      (-MbSL2 + pMT2 + MUEr2)*x6*
     &       ((MbSL2**2 - 3*MbSL2*pMT2 + 2*pMT2**2 - 
     &            2*MbSL2*MUEr2 - 3*pMT2*MUEr2 + MUEr2**2)*(x10 - x7)
     &           + MUEr2*(2*MbSL2*x24 + 
     &            (MbSL2 + pMT2 - MUEr2)*(-x10 + x8))))/pMT2 + 
     &   (pSB2*(2*ctt2*pMT*x1 + 2*pMT*stt2*x3 + s2tt*(x1 - x3)*Xtr))/
     &    pMT - Pi**2/3.D0*pSB2 + 
     &   (-x10 + x12)*(MbSL2*
     &       (tmp19/MSt22 + tmp20/MSt12 + 
     &         (-1 + (-MbSL2 + pMT2 + MUEr2)**2*x6)/pMT2) - x7 + 
     &      x8 + 1/2.D0*(pCB2*tmp23)) - 
     &   (ctt2*pMT2 + MbSL2*stt2 + pCB2*MA02*(1 + stt2) - 
     &      2*MUEr2*(1 + stt2) - 
     &      x13*(2*MUEr2**2*(-MSt12 + pMT2 + MUEr2)*x27 + 
     &         ((MSt12 - pMT2 - MUEr2)*
     &             ((MSt12 - pMT2)**2*(MSt12 + pMT2) - 
     &               (2*MSt12**2 + 5*MSt12*pMT2 + pMT2**2)*MUEr2 + 
     &               MSt12*MUEr2**2)*(-x10 + x7) + 
     &            MUEr2*tmp8*(x10 - x8))/pMT2) + 
     &      pCB2*MA02*x25*(2 + (s2tt*Ytr)/pMT + 
     &         1/2.D0*(2*pMT + s2tt*Ytr)**2/(MA02 - 4*MSt12)) - 
     &      1/4.D0*(x22*(tmp22 - 
     &            2*(1 + c2tt2)*pCB2*MA02**2*MSt22*x26*Ytr**2))/MSt22
     &      )/MSt12 - (pCB2*(1 + ctt2)*MA02 + ctt2*MbSL2 - 
     &      2*(1 + ctt2)*MUEr2 + pMT2*stt2 - 
     &      x14*(2*MUEr2**2*(-MSt22 + pMT2 + MUEr2)*x30 + 
     &         ((MSt22 - pMT2 - MUEr2)*
     &             ((MSt22 - pMT2)**2*(MSt22 + pMT2) - 
     &               (2*MSt22**2 + 5*MSt22*pMT2 + pMT2**2)*MUEr2 + 
     &               MSt22*MUEr2**2)*(-x10 + x7) - 
     &            MUEr2*((MSt22 - pMT2)**3 + 
     &               2*MSt22*(-MSt22 + pMT2)*MUEr2 + 
     &               (MSt22 + pMT2)*MUEr2**2)*(-x10 + x8))/pMT2) + 
     &      pCB2*MA02*x29*(2 - (s2tt*Ytr)/pMT + 
     &         1/2.D0*(-2*pMT + s2tt*Ytr)**2/(MA02 - 4*MSt22)) - 
     &      1/4.D0*(x21*(tmp22 - 
     &            2*(1 + c2tt2)*pCB2*MA02**2*MSt12*x28*Ytr**2))/MSt12
     &      )/MSt22

        F1 = F1 + (-x10 + x9)*
     &    ((-2 + s2tt2 + Nc*s2tt2)*(x10 - x11) + 
     &      (-MSt12 + ((MSt12 - pMT2 - MUEr2)*
     &            (MSt12**3 + pMT2*(pMT2 - MUEr2)**2 - 
     &              MSt12**2*(pMT2 + MUEr2) - 
     &              MSt12*pMT2*(pMT2 + 5*MUEr2))*x13)/MSt12)/pMT2 + 
     &      2*(x10 - x7) + 
     &      1/4.D0*((MSt22 + MSt2diff)*Nc*s2tt2 + 
     &          MSt12*(4 - (2 + Nc)*s2tt2) + 
     &          2*(1 + c2tt2)*pSB2*Xtr**2)/MSt22 + 
     &      (MUEr2 + pMT2*(-1 + (-1 + pCB2)*ctt2 + 4*pSB2) + 
     &         ((MbSL2 - 2*MSt12)*pSB2*tmp13)/(MbSL2 - MSt12) + 
     &         stt2*(2*MSt12 + c2tt*pSB2*Xtr**2) + 
     &         (s2tt*((-3*MSt12 + pMT2)*pSB2*Xtr - 
     &              3*pCB2*(MSt12 + pMT2)*Ytr))/pMT + 
     &         pCB2*(Ytr*(pMT*s2tt + stt2*Ytr) + 
     &            1/2.D0*((MA02 - 8*MSt12)*(2*pMT + s2tt*Ytr)**2)/
     &              (MA02 - 4*MSt12)))/MSt12 + 
     &      1/2.D0*(-13 - c2tt2 - 2*stt2 + 
     &         ((-3 + c2tt)*pCB2*Ytr**2)/MSt12)) + 
     &   (-x10 + x11)*((-MSt22 + 
     &         ((MSt22 - pMT2 - MUEr2)*
     &            (MSt22**3 + pMT2*(pMT2 - MUEr2)**2 - 
     &              MSt22**2*(pMT2 + MUEr2) - 
     &              MSt22*pMT2*(pMT2 + 5*MUEr2))*x14)/MSt22)/pMT2 + 
     &      2*(x10 - x7) + 
     &      (MUEr2 + ((MbSL2 - 2*MSt22)*pSB2*tmp11)/
     &          (MbSL2 - MSt22) + 
     &         (s2tt*((3*MSt22 - pMT2)*pSB2*Xtr + 
     &              3*pCB2*(MSt22 + pMT2)*Ytr))/pMT + 
     &         pCB2*(pMT2*stt2 - pMT*s2tt*Ytr) + 
     &         ctt2*(2*MSt22 - c2tt*pSB2*Xtr**2 + pCB2*Ytr**2))/MSt22+
     &        1/4.D0*(4*MSt22 - 2*(MSt22 + MSt2diff*Nc)*s2tt2 + 
     &          2*(1 + c2tt2)*pSB2*Xtr**2)/MSt12 + 
     &      1/2.D0*(-13 - c2tt2 - 2*ctt2 + 
     &         (2*pMT2*(-1 + 4*pSB2 - stt2) - (3 + c2tt)*pCB2*Ytr**2 + 
     &            (pCB2*(MA02 - 8*MSt22)*(-2*pMT + s2tt*Ytr)**2)/
     &             (MA02 - 4*MSt22))/MSt22)) + 
     &   1/2.D0*(((x10 - x11)**2*
     &          (pMT*(6 + c2tt*pCB2 - c2tt2*(1 + Nc) + pSB2) - 
     &            s2tt*(6*pSB2*Xtr + 5*pCB2*Ytr)) + 
     &         (x10 - x9)**2*
     &          (pMT*(5 - c2tt*pCB2 - Nc + s2tt2 + Nc*s2tt2 + pSB2) + 
     &            s2tt*(6*pSB2*Xtr + 5*pCB2*Ytr)) - 
     &         (3*(MSt22 - 2*pMT2)*s2tt*(pSB2*Xtr + pCB2*Ytr) + 
     &            pCB2*(MA02 + MbSL2 - MSt22)*x16*
     &             (2*pMT*stt2 - s2tt*Ytr))/MSt22)/pMT + 
     &      ((-((1 + c2tt2)*(MSt12**2 + MSt22**2)) + 
     &            MSt2diff**2*Nc*s2tt2 - 
     &            2*(ctt2*MSt12 + MSt2sum + MSt22*stt2)*tmp4)/MSt22
     &           + (3*(MSt12 - 2*pMT2)*s2tt*(pSB2*Xtr + pCB2*Ytr) - 
     &            pCB2*(MA02 + MbSL2 - MSt12)*x15*
     &             (2*ctt2*pMT + s2tt*Ytr))/pMT)/MSt12)

        F1 = F1 + pCB2*(-2*x17 + 
     &      tmp14*((MbSL2*tmp10)/MSt22 + 
     &         (MA02 + MbSL2 - MSt22)*(x10 - x11))*x19 + 
     &      (2*MA02*(MA02 - 6*pMT2)*x23)/((MA02 - 4*pMT2)*pMT2) + 
     &      tmp16*x20*((MbSL2*tmp9)/MSt12 + 
     &         (MA02 + MbSL2 - MSt12)*(x10 - x9)) + 
     &      (-x10 + x18)*(x10 - x12 + 1/2.D0*tmp23 + 
     &         MA02*(4/(MA02 - 4*pMT2) + 
     &            (1 + ctt2 + (MA02 - MbSL2 - MSt22)*tmp14*x19)/
     &             MSt22 + 
     &            (1 + stt2 + (MA02 - MbSL2 - MSt12)*tmp16*x20 - 
     &               1/4.D0*((1 + c2tt2)*
     &                   (MSt2diff**2 - MA02*MSt2sum)*(x21 + x22)*
     &                   Ytr**2)/MSt22)/MSt12 + 
     &            1/2.D0*((-2*pMT + s2tt*Ytr)**2/
     &                ((MA02 - 4*MSt22)*MSt22) + 
     &               (2*pMT + s2tt*Ytr)**2/((MA02 - 4*MSt12)*MSt12)))
     &         ))

	tmp5 = MSt2sum + MSt2sum*Nc + pSB2*Xtr**2 + pCB2*Ytr**2

	tmp6 = (MA02 - MbSL2)**2 - (MA02 + MbSL2)*MSt22

	tmp17 = MbSL2 - pMT2 - 2*MUEr2 + pSB2*Xtr**2 + pCB2*(MA02 + Ytr**2)

        tmp18 = -3*(MA02 + MbSL2) + 
     &   (MA02 - MbSL2)**2/MSt12 + 2*MSt12

        s2ttF2 = 1/8.D0*((8*
     &        (-3*Atr*pMT*s2tt2 + 
     &          s2tt*x14*(-4*pMT2*MUEr2**2*x30 - 
     &             ((MSt22 - pMT2)**3 - 
     &                2*MSt22*(MSt22 - pMT2)*MUEr2 + 
     &                (MSt22 + pMT2)*MUEr2**2)*(-x10 + x7) + 
     &             MUEr2*(-(MSt22 - pMT2)**2 + MUEr2**2)*(-x10 + x8)))
     &        )/MSt22 + 8*MbSL2*(-x10 + x12)*
     &      (s2tt*((-ctt2 + (pSB2*tmp11)/(-MbSL2 + MSt22))/
     &            MSt22 + tmp20/MSt12) + 
     &        (Atr*s2tt2*(-x11 + x9))/(MSt2diff*pMT)) + 
     &     ((x10 - x9)**2*(-4*Atr*
     &           (6*c2tt2*pMT2 + 
     &             (MbSL2 - 3*(2*MSt12 + MSt2diff))*s2tt2) + 
     &          c2tt*pMT*s2tt*
     &           (-4*pCB2*MA02 - 4*MbSL2 - 7*MSt12 + MSt22 - 
     &             4*MSt2diff + MSt2sum + 4*pMT2 + 8*MUEr2 - 4*tmp4)+
     &            2*(MSt2diff*pMT*s2tt*
     &              (5 - 2*Nc + 4*s2tt2 + 4*Nc*s2tt2) + 
     &             c2tt2*pMT*s2tt*
     &              (MSt22 + MSt2sum + MSt22*Nc + MSt2sum*Nc - 
     &                11*MSt12*(1 + Nc) - 4*tmp4) - 
     &             4*pCB2*MA02*s2tt2*Ytr)))/(MSt2diff*pMT) + 
     &     ((x10 - x11)**2*
     &        (4*Atr*(6*c2tt2*pMT2 + 
     &             (MbSL2 - 6*MSt22 + 3*MSt2diff)*s2tt2) + 
     &          c2tt*pMT*s2tt*
     &           (4*pCB2*MA02 + 4*MbSL2 - MSt12 + 7*MSt22 - 
     &             4*MSt2diff - MSt2sum - 4*pMT2 - 8*MUEr2 + 4*tmp4)+
     &            2*(-(MSt2diff*pMT*s2tt*
     &                (5 + 4*s2tt2 + Nc*(-2 + 4*s2tt2))) + 
     &             c2tt2*pMT*s2tt*
     &              (MSt2sum + MSt2sum*Nc + MSt12*(1 + Nc) - 
     &                11*MSt22*(1 + Nc) - 4*tmp4) + 
     &             4*pCB2*MA02*s2tt2*Ytr)))/(MSt2diff*pMT) + 
     &     (8*(-3*Atr*pMT*s2tt2 + 
     &          (MSt2diff*pMT2*s2tt*(x10 - x7) - 
     &             MUEr2*s2tt*
     &              (-2*ctt2*MSt12 + MSt2diff + 2*MSt22*stt2)*
     &              (-x10 + x8))/MSt22 + 
     &          s2tt*(x13*(4*pMT2*MUEr2**2*x27 + tmp8*(-x10 + x7) + 
     &                MUEr2*((MSt12 - pMT2)**2 - MUEr2**2)*(-x10 + x8)
     &                ) - (c2tt2*pCB2*MA02*x26*Ytr**2)/MSt2diff)))/
     &      MSt12 + 4*(c2tt*s2tt + 
     &        (pSB2*(x1*((-4*c2tt2*pMT2 + 
     &                    (-2*MbSL2 + MSt2sum)*s2tt2)*Xtr + 
     &                 2*pMT*s2tt*
     &                  (ctt2*MSt2diff + 
     &                    c2tt*(-MbSL2 + MSt12 + pMT2 - Xtr**2)))+
     &                x3*((4*c2tt2*pMT2 + 
     &                    (2*MbSL2 - 2*MSt22 + MSt2diff)*s2tt2)*
     &                  Xtr + 
     &                 pMT*s2tt*
     &                  (-2*MSt2diff*stt2 + 
     &                    2*c2tt*(MbSL2 - MSt22 - pMT2 + Xtr**2))))
     &             + (pCB2*x15*(c2tt*pMT*s2tt + s2tt2*Ytr))/
     &            (MSt12*x20))/(MSt2diff*pMT)))

        s2ttF2 = s2ttF2 + 
     &   (-x10 + x9)*((s2tt*
     &          (MUEr2 + pMT2*(3 - 4*pCB2 + ctt2*(-1 + pSB2)) + 
     &            (MSt12 - pMT2 - MUEr2)*
     &             (MSt12**2 + (pMT2 - MUEr2)**2 - 
     &               2*MSt12*(pMT2 + 2*MUEr2))*x13) + 
     &         pSB2*Xtr*(2*pMT*s2tt2 + s2tt*stt2*Xtr) - 
     &         3*pCB2*pMT*s2tt2*Ytr)/MSt12 + 
     &      (6*Atr*c2tt2*pMT + 
     &         s2tt*(-(c2tt2*
     &               (MSt22 + MSt22*Nc - 3*MSt12*(1 + Nc) + 
     &                 2*tmp5*x10 - 2*tmp5*x11)) + 
     &            c2tt*(MbSL2 + MSt12 - pMT2 - 2*MUEr2 + 
     &               2*ctt2*pSB2*Xtr**2 + pCB2*(MA02 + 2*ctt2*Ytr**2)))
     &         )/MSt2diff - 
     &      3/2.D0*(Atr*(6*MSt12 + MSt2diff)*s2tt2)/(MSt2diff*pMT) + 
     &      1/4.D0*(s2tt*(-20 + (4*pSB2*tmp13)/(-MbSL2 + MSt12) + 
     &           ((-3 + 2*c2tt - c2tt2)*tmp4)/MSt12 - 
     &           (8*pCB2*(2*pMT + s2tt*Ytr)**2)/(MA02 - 4*MSt12) + 
     &           (2*(s2tt2*(MSt12 + MSt12*Nc + pSB2*Xtr**2) + 
     &                4*pCB2*(pMT2 + pMT*s2tt*Ytr)))/MSt12 + 
     &           (2*MSt12**2*MSt2diff*(-2 + s2tt2 + Nc*s2tt2) - 
     &              4*c2tt2*MSt12*MSt22*tmp5 - 
     &              (1 + c2tt2)*MSt2diff*
     &               (MSt22*pSB2*Xtr**2 - 2*MSt2diff*pSB2*Xtr**2 + 
     &                 pCB2*MSt22*Ytr**2))/(MSt12*MSt22*MSt2diff))))+
     &     1/4.D0*(s2tt*((6*pCB2*MA02 + 2*MbSL2 + 4*MSt12 + 
     &           c2tt*MSt22 + 2*pMT2 - 12*MUEr2 - 2*MSt12*s2tt2 - 
     &           2*MSt12*Nc*s2tt2 + 2*c2tt*tmp17 + 6*tmp4)/MSt22 + 
     &        8*((c2tt*(-MSt12 + MUEr2))/MSt2diff + stt2)*x4 - 
     &        8*(ctt2 + (c2tt*(-MSt22 + MUEr2))/MSt2diff)*x5 - 
     &        (4*c2tt2*pCB2*MA02*x28*Ytr**2)/(MSt22*MSt2diff) + 
     &        (-2*(MbSL2 + MSt22 + pMT2 - 6*MUEr2) + 
     &           2*MSt22*(-c2tt2 + Nc*s2tt2) + 
     &           c2tt*(MSt12 + 2*tmp17) - 
     &           6*(pSB2*Xtr**2 + pCB2*(MA02 + Ytr**2)))/MSt12)) + 
     &   (-9*Atr*s2tt2 + 1/4.D0*
     &       (pCB2*((-2*x16*(c2tt*pMT*s2tt + s2tt2*Ytr))/
     &             (MSt22*x19) + 
     &            (x12 - x9)*(x18 - x9)*
     &             ((4*c2tt2*pMT2 + 
     &                  (2*(MA02 + MbSL2) - MSt2sum)*s2tt2)*Ytr + 
     &               pMT*s2tt*
     &                (-2*ctt2*MSt2diff + 
     &                  2*c2tt*
     &                   (MA02 + MbSL2 - MSt12 - pMT2 + Ytr**2)))+
     &              (-x11 + x12)*(-x11 + x18)*
     &             ((-4*c2tt2*pMT2 + 
     &                  (-2*(MA02 + MbSL2) + MSt2sum)*s2tt2)*Ytr+
     &                 pMT*s2tt*
     &                (2*MSt2diff*stt2 - 
     &                  2*c2tt*
     &                   (MA02 + MbSL2 - MSt22 - pMT2 + Ytr**2)))))
     &         /MSt2diff)/pMT

        s2ttF2 = s2ttF2 + 
     &   pCB2*(s2tt*(tmp14*((2*(MA02 + MbSL2 - MSt22) - 
     &               tmp6/MSt22)*(-x10 + x11) + 
     &            (MbSL2*
     &               ((MA02 - MbSL2 + MSt22)*(-x10 + x12) + 
     &                 2*MA02*x16))/MSt22)*x19 + 
     &         tmp16*x20*((MbSL2*tmp9)/MSt12 + 
     &            tmp18*(-x10 + x9)) + 
     &         1/4.D0*(((x21*
     &                 ((-1 - c2tt2)*
     &                    (((MA02*
     &                       (MSt12*(MA02 - 2*MSt12 - 3*MSt22) + 
     &                       MSt22**2) - MSt2diff**3)*(-x10 + x11))
     &                       /MSt12 + 2*MA02*MSt2diff*x28) - 
     &                   ((1 + c2tt2)*
     &                      (MSt2diff**3 + 
     &                       MA02**2*(MSt12 + MSt2diff) + 
     &                       MA02*
     &                       (2*MSt12**2 - MSt22**2 - 
     &                       MSt2diff*MSt2sum - 
     &                       MSt12*(MSt22 + 2*MSt2sum)))*(x10 - x9)
     &                      )/MSt12))/MSt22 + 
     &              ((1 + c2tt2)*x22*
     &                 (((MA02**2*(MSt22 - MSt2diff) - 
     &                       MSt2diff**3 - 
     &                       MA02*
     &                       (MSt12**2 + MSt12*MSt22 - 
     &                       2*MSt22**2 + 2*MSt22*MSt2sum - 
     &                       MSt2diff*MSt2sum))*(x10 - x11))/MSt22-
     &                     2*MA02*MSt2diff*x26 + 
     &                   ((MA02*
     &                       (MSt12**2 + 
     &                       (MA02 - 3*MSt12 - 2*MSt22)*MSt22) + 
     &                       MSt2diff**3)*(-x10 + x9))/MSt22))/
     &               MSt12)*Ytr**2)) + 
     &      1/4.D0*((MA02 + MbSL2 - MSt22)*x16*
     &          (2*MSt2diff*pMT*s2tt*stt2 - 4*c2tt2*pMT2*Ytr - 
     &            MSt2diff*s2tt2*Ytr + 2*c2tt*pMT*s2tt*(pMT2 - Ytr**2))
     &          )/(MSt22*MSt2diff*pMT) - 
     &      1/4.D0*((MA02 + MbSL2 - MSt12)*x15*
     &          (2*ctt2*MSt2diff*pMT*s2tt - 4*c2tt2*pMT2*Ytr + 
     &            MSt2diff*s2tt2*Ytr + 2*c2tt*pMT*s2tt*(pMT2 - Ytr**2))
     &          )/(MSt12*MSt2diff*pMT) + 
     &      MA02*((x29*(s2tt + 
     &              (Ytr*(c2tt2*(-2*pMT + s2tt*Ytr) - 
     &                   1/2.D0*
     &                    ((MA02 - 4*MSt22 + MSt2diff)*s2tt2)/pMT))/
     &               MSt2diff + 
     &              1/2.D0*(s2tt*(-2*pMT + s2tt*Ytr)**2)/
     &                (MA02 - 4*MSt22)))/MSt22 + 
     &         (x25*((Ytr*(c2tt2*(2*pMT + s2tt*Ytr) + 
     &                   1/2.D0*
     &                    ((MA02 - 4*MSt12 - MSt2diff)*s2tt2)/pMT))/
     &               MSt2diff + 
     &              s2tt*(-1 - 
     &                 1/2.D0*(2*pMT + s2tt*Ytr)**2/(MA02 - 4*MSt12))))
     &           /MSt12 + (-x10 + x18)*
     &          ((-2*s2tt2*(x11 - x9)*Ytr)/(MSt2diff*pMT) + 
     &            s2tt*(-(((MA02 - MbSL2 - MSt22)*tmp14*x19)/
     &                  MSt22) + 
     &               ((MA02 - MbSL2 - MSt12)*tmp16*x20 + 
     &                  1/4.D0*
     &                   ((1 + c2tt2)*MSt2diff*(MA02 - MSt2sum)*
     &                      (x21 + x22)*Ytr**2)/MSt22)/MSt12 + 
     &               1/2.D0*(-((2*(1 + ctt2) + 
     &                       (-2*pMT + s2tt*Ytr)**2/(MA02 - 4*MSt22))
     &                      /MSt22) + 
     &                  (2*(1 + stt2) + 
     &                     (2*pMT + s2tt*Ytr)**2/(MA02 - 4*MSt12))/
     &                   MSt12)))))

        s2ttF2 = s2ttF2 + 
     &   (-x10 + x11)*((Atr*(-3/2.D0 + (9*MSt22)/MSt2diff)*s2tt2)/
     &       pMT + (-6*Atr*c2tt2*pMT + 
     &         s2tt*(c2tt2*(2*MSt22 + MSt2diff)*(1 + Nc) - 
     &            c2tt*(pCB2*MA02 + MbSL2 + MSt22 - pMT2 - 
     &               2*MUEr2 + 2*pSB2*stt2*Xtr**2 + 2*pCB2*stt2*Ytr**2))
     &         )/MSt2diff + 
     &      (pMT*s2tt2*(2*pSB2*Xtr - 3*pCB2*Ytr) + 
     &         s2tt*(-2*pCB2*pMT2 - MUEr2 + 
     &            pMT2*(-3 + 4*pCB2 + stt2 - pSB2*stt2) - 
     &            (MSt22 - pMT2 - MUEr2)*
     &             (MSt22**2 + (pMT2 - MUEr2)**2 - 
     &               2*MSt22*(pMT2 + 2*MUEr2))*x14 - 
     &            ctt2*pSB2*Xtr**2 + 2*pCB2*pMT*s2tt*Ytr + 
     &            1/4.D0*((3 + 2*c2tt + c2tt2)*tmp4) - 
     &            1/2.D0*(s2tt2*pSB2*Xtr**2)))/MSt22 + 
     &      s2tt*(5 + (pSB2*tmp11)/(MbSL2 - MSt22) + 
     &         (2*pCB2*(-2*pMT + s2tt*Ytr)**2)/(MA02 - 4*MSt22) - 
     &         1/2.D0*((1 + Nc)*s2tt2) + 
     &         1/4.D0*(-2*MSt22**2*MSt2diff*
     &              (-2 + s2tt2 + Nc*s2tt2) - 
     &             4*c2tt2*MSt12*MSt22*tmp5 + 
     &             (1 + c2tt2)*MSt2diff*
     &              (MSt12*pSB2*Xtr**2 + 2*MSt2diff*pSB2*Xtr**2 + 
     &                pCB2*MSt12*Ytr**2))/(MSt12*MSt22*MSt2diff)))

	lc10 = (-6*c2tt*MUEr2**2*s2tt)/MSt2diff**2

	lc11 = (MSt2sum*(2 + Nc)*s2tt)/MSt2diff

        lc12 = (4*MSt2sum*c2tt2*s2tt*
     &      (MSt2sum*(1 + Nc) + 2*(pSB2*Xtr**2 + pCB2*Ytr**2)))/
     &    MSt2diff**2

        lc13 = -(MA02*pCB2*s2tt*
     &      ((6*c2tt)/MSt2diff - 
     &        1/4.D0*(6 - 2*c2tt - 
     &            ((1 + c2tt2)*
     &               (MSt2sum*(2*MSt2sum - MA02) - MSt2diff**2)*
     &               (x21 + x22)*Ytr**2)/MSt22 + 
     &            (2*(2*pMT + s2tt*Ytr)**2)/(MA02 - 4*MSt12) + 
     &            (2*MA02 - 2*(MbSL2 + MSt12))*x20*
     &             (pMT2 + 2*pMT*s2tt*Ytr + Ytr**2 + 
     &               c2tt*(pMT2 - Ytr**2)))/MSt12 - 
     &        1/2.D0*(3 + c2tt + 
     &            (2*pMT - s2tt*Ytr)**2/(MA02 - 4*MSt22) + 
     &            1/2.D0*((2*MA02 - 2*(MbSL2 + MSt22))*x19*
     &               (pMT2 - 2*pMT*s2tt*Ytr + Ytr**2 - 
     &                 c2tt*(pMT2 - Ytr**2))))/MSt22))

        lc14 = -(1/4.D0*(MbSL2*s2tt*
     &        ((24*c2tt)/MSt2diff - 
     &          ((1 - c2tt)*(2*MbSL2 - 2*MSt12) + 
     &             2*pMT2*(pSB2 + c2tt*pSB2) + 
     &             pSB2*(4*pMT*s2tt*Xtr + (2 - 2*c2tt)*Xtr**2))/
     &           ((MbSL2 - MSt12)*MSt12) - 
     &          (2*((1 + c2tt)*(MbSL2 - MSt22) + 
     &               pSB2*(pMT2 - 2*pMT*s2tt*Xtr + Xtr**2 - 
     &                  c2tt*(pMT2 - Xtr**2))))/
     &           ((MbSL2 - MSt22)*MSt22) + 
     &          pCB2*(((2*MA02 - 2*MbSL2 + 2*MSt22)*x19*
     &                (pMT2 - 2*pMT*s2tt*Ytr + Ytr**2 - 
     &                  c2tt*(pMT2 - Ytr**2)))/MSt22 + 
     &             ((2*MA02 - 2*MbSL2 + 2*MSt12)*x20*
     &                (pMT2 + 2*pMT*s2tt*Ytr + Ytr**2 + 
     &                  c2tt*(pMT2 - Ytr**2)))/MSt12))))

        lc15 = -(1/4.D0*(12*MSt2sum*Atr*pMT*(2 - 3*s2tt2) + 
     &         s2tt*(2*(MbSL2*(3*MSt2sum*c2tt - MSt2diff) + 
     &               MSt2diff*
     &                (6*MUEr2 - 
     &                  3*(MSt2sum*(2 + Nc + c2tt2*(1 + Nc)) + 
     &                     MA02*pCB2) - pMT2 - 
     &                  3*(pSB2*Xtr**2 + pCB2*Ytr**2)) + 
     &               4*MSt2sum*c2tt2*
     &                (MSt2sum*(1 + Nc) + 
     &                  2*(pSB2*Xtr**2 + pCB2*Ytr**2))) - 
     &            c2tt*(MSt2diff**2 - 
     &               3*MSt2sum*
     &                (MSt2sum - 4*MUEr2 + 2*MA02*pCB2 - 2*pMT2 + 
     &                  2*(pSB2*Xtr**2 + pCB2*Ytr**2)))))/
     &       MSt2diff**2)

        lc16 = -((6*Atr*(2*MSt2sum - MSt2diff)*pMT*
     &          (2 - 3*s2tt2) + 
     &         s2tt*(MbSL2*(3*MSt2sum*c2tt - MSt2diff) - 
     &            MSt2diff*(3*MA02*pCB2 + pMT2) - 
     &            3/2.D0*(c2tt*
     &               (-(MSt2sum*
     &                    (2*MSt2sum - 6*MUEr2 + 2*MA02*pCB2 - 
     &                      4*pMT2)) + 
     &                 MSt2diff*(2*MSt22 - 6*MUEr2 - 2*pMT2)))))/
     &       MSt2diff**2) - 
     &    (-3*pCB2*pMT*s2tt2*Ytr + 
     &       pSB2*Xtr*(2*pMT*s2tt2 + 
     &          s2tt*Xtr*(1/2.D0 + s2tt2 - 1/2.D0*c2tt)) + 
     &       s2tt*(MUEr2 + 1/2.D0*((1 + c2tt)*pMT2*pSB2) + 
     &          1/4.D0*((2*MSt12 - 2*(MUEr2 + pMT2))*
     &             (2*(MUEr2 - pMT2)**2 + 
     &               MSt2sum*(MSt2sum - 3*(MUEr2 + pMT2)) - 
     &               MSt2diff*
     &                (2*MSt2sum - MSt2diff - 3*(MUEr2 + pMT2)))*
     &             x13)) + 
     &       1/2.D0*((5 - c2tt)*pMT2*s2tt - 
     &          pCB2*(4*pMT2*s2tt + 
     &             s2tt2*(-4*pMT*Ytr - s2tt*Ytr**2))))/MSt12 + 
     &    1/4.D0*(s2tt*(8 + 4*Nc - 
     &         (MSt12*(MSt2diff*
     &                (MSt2sum*(1 - Nc) - 
     &                  MSt2diff*(1 - Nc + c2tt2*(2 + 2*Nc))) + 
     &               2*MSt22*
     &                (12*(1 - c2tt)*MUEr2 + 
     &                  MSt2diff*
     &                   (7 + c2tt + 7*Nc + c2tt2*(6 + 6*Nc)) - 
     &                  MSt2sum*(15 + 4*Nc + c2tt2*(16 + 16*Nc))))-
     &              2*(1 + c2tt2)*(MSt2sum**2 - 2*MSt22*MSt2diff)*
     &             pSB2*Xtr**2)/(MSt12*MSt22*MSt2diff) + 
     &         (2*pSB2*(pMT2 + 2*pMT*s2tt*Xtr + Xtr**2 + 
     &              c2tt*(pMT2 - Xtr**2)))/(MbSL2 - MSt12) + 
     &         (8*pCB2*(2*pMT + s2tt*Ytr)**2)/(MA02 - 4*MSt12) + 
     &         (2*(3*MSt2sum*MSt2diff - 
     &              c2tt*(MSt2sum*(6*MSt2sum - 9*MSt2diff) + 
     &                 4*MSt2diff**2))*(pSB2*Xtr**2 + pCB2*Ytr**2))
     &           /(MSt12*MSt2diff**2) + 
     &         1/2.D0*(pCB2*(((1 + c2tt2)*
     &                 (2*MA02*(MSt2sum**2 + MSt22*MSt2diff) - 
     &                   2*
     &                    (MA02**2*MSt22 - 
     &                      MSt2diff*
     &                       (MSt2sum**2 - 2*MSt22*MSt2diff)))*
     &                 (x21 + x22)*Ytr**2)/(MSt12*MSt22) + 
     &              4*(3*(MA02 + MbSL2) - 
     &                 (MA02 - MbSL2)**2/MSt12 - 2*MSt12)*x20*
     &               (pMT2 + 2*pMT*s2tt*Ytr + Ytr**2 + 
     &                 c2tt*(pMT2 - Ytr**2))))))

        lc17 = 1/4.D0*(12*MSt2sum*Atr*pMT*(2 - 3*s2tt2) + 
     &       s2tt*(2*(MbSL2*(3*MSt2sum*c2tt - MSt2diff) + 
     &             MSt2diff*
     &              (6*MUEr2 - 
     &                3*(MSt2sum*(2 + Nc + c2tt2*(1 + Nc)) + 
     &                   MA02*pCB2) - pMT2 - 
     &                3*(pSB2*Xtr**2 + pCB2*Ytr**2)) - 
     &             4*MSt2sum*c2tt2*
     &              (MSt2sum*(1 + Nc) + 
     &                2*(pSB2*Xtr**2 + pCB2*Ytr**2))) - 
     &          c2tt*(MSt2diff**2 - 
     &             3*MSt2sum*
     &              (MSt2sum - 4*MUEr2 + 2*MA02*pCB2 - 2*pMT2 + 
     &                2*(pSB2*Xtr**2 + pCB2*Ytr**2)))))/MSt2diff**2

        lc18 = 1/2.D0*(((2*
     &             (MbSL2*(3*MSt2sum*c2tt - MSt2diff) + 
     &               MSt2diff*(-3*MA02*pCB2 - pMT2)) + 
     &            3*c2tt*(2*MSt12*MSt2diff + 
     &               MSt2sum*
     &                (2*MSt2sum - 6*MUEr2 + 2*MA02*pCB2 - 4*pMT2)-
     &                 MSt2diff*(6*MUEr2 + 2*pMT2)))*s2tt - 
     &         12*Atr*pMT*(2*c2tt2*MSt2diff - 
     &            MSt2sum*(4 - 6*s2tt2) - MSt2diff*(4 - 5*s2tt2)))/
     &       MSt2diff**2 - 
     &      (pSB2*(pMT2*(s2tt - c2tt*s2tt) - 
     &            Xtr*(4*pMT*s2tt2 - 
     &               s2tt*(Xtr + c2tt*Xtr + 2*s2tt2*Xtr))) + 
     &         pCB2*s2tt2*Ytr*(2*pMT + s2tt*Ytr) + 
     &         s2tt*((5 + c2tt - 4*pCB2)*pMT2 + 
     &            2*(MUEr2 + 
     &               1/2.D0*((2*MSt22 - 2*(MUEr2 + pMT2))*
     &                  (MUEr2*(MUEr2 - 2*pMT2) + pMT2**2 + 
     &                    MSt22*(2*MSt22 - 3*(MUEr2 + pMT2)))*x14))
     &            ))/MSt22 + 
     &      s2tt*(11 + c2tt + 2*Nc - 
     &         (MSt22*(1 - Nc) + c2tt2*MSt2diff*(1 + Nc))/MSt12 - 
     &         s2tt2 - Nc*s2tt2 + 
     &         (MSt12*MSt22*
     &             (12*(1 + c2tt)*MUEr2 - 
     &               MSt2diff*(13 + 6*Nc + c2tt2*(7 + 7*Nc)) - 
     &               MSt2sum*(15 + 4*Nc + c2tt2*(16 + 16*Nc))) - 
     &            (1 + c2tt2)*(MSt2sum**2 + 2*MSt12*MSt2diff)*pSB2*
     &             Xtr**2)/(MSt12*MSt22*MSt2diff) + 
     &         (pSB2*(pMT2 - 2*pMT*s2tt*Xtr + Xtr**2 - 
     &              c2tt*(pMT2 - Xtr**2)))/(MbSL2 - MSt22) + 
     &         (4*pCB2*(2*pMT - s2tt*Ytr)**2)/(MA02 - 4*MSt22) - 
     &         1/2.D0*((6*MSt2sum*MSt2diff - 
     &               2*c2tt*
     &                (4*MSt2diff**2 + 
     &                  MSt2sum*(6*MSt2sum + 9*MSt2diff)))*
     &             (pSB2*Xtr**2 + pCB2*Ytr**2))/(MSt22*MSt2diff**2)
     &           + 1/4.D0*(pCB2*
     &            (((1 + c2tt2)*
     &                 (MA02*(2*MSt2sum**2 - 2*MSt12*MSt2diff) - 
     &                   2*
     &                    (MA02**2*MSt12 + MSt2sum**2*MSt2diff + 
     &                      2*MSt12*MSt2diff**2))*(x21 + x22)*
     &                 Ytr**2)/(MSt12*MSt22) + 
     &              4*(2*(MA02 + MbSL2) - 
     &                 (MA02*(MA02 - 2*MbSL2 - MSt22) + 
     &                    MbSL2*(MbSL2 - MSt22))/MSt22 - 2*MSt22)*
     &               x19*(pMT2 - 2*pMT*s2tt*Ytr + Ytr**2 - 
     &                 c2tt*(pMT2 - Ytr**2))))))

        lc19 = -(1/8.D0*(MUEr2*s2tt*
     &         (4*c2tt*(3*MSt2sum**2 - MSt2diff**2) - 
     &           MSt2diff*(MSt2sum**3*(x13 + x14) + 
     &              MSt2diff*
     &               (4*(MUEr2**2 - pMT2**2)*(x13 - x14) + 
     &                 MSt2diff*
     &                  (4*MUEr2*x13 + MSt2diff*(x13 - x14) + 
     &                    4*MUEr2*x14)) - 
     &              MSt2sum*
     &               (MSt2diff**2*(x13 + x14) + 
     &                 MSt2sum*
     &                  (4*MUEr2*x13 + MSt2diff*(x13 - x14) + 
     &                    4*MUEr2*x14) - 
     &                 4*(4 + (MUEr2**2 - pMT2**2)*(x13 + x14))))))
     &        /(MSt12*MSt22*MSt2diff))

        lc20 = 1/8.D0*(pMT2*s2tt*
     &       (8*MSt2sum + 2*MSt22*
     &          (4*MSt12**2 - 4*(MUEr2**2 - pMT2**2))*x13 + 
     &         8*MSt12*((MSt22**2 - MUEr2**2 + pMT2**2)*x14 - 
     &            2*MSt22*pMT2*(x13 + x14))))/(MSt12*MSt22)

        lc21 = 1/2.D0*(-((2*pSB2*(x1 - x3)*
     &            (MSt2diff*s2tt*(pMT2 + Xtr**2) - 
     &              (MSt2sum - 2*MbSL2)*
     &               (2*pMT*(2 - 3*s2tt2)*Xtr - 
     &                 3*c2tt*s2tt*(pMT2 - Xtr**2))) + 
     &           s2tt*((c2tt*
     &                  (MSt2diff**2 - 3*(MSt2sum - 2*MUEr2)**2) + 
     &                 2*MSt2diff*(MSt2sum - 2*MUEr2))*x4 - 
     &              (MSt2diff*
     &                  ((4 + 
     &                       x14*
     &                       (8*MUEr2*pMT2 + 
     &                       MSt2diff*
     &                       (2*pMT2 + 
     &                       MUEr2*
     &                       (2 - 4*(MUEr2*pMT2 + pMT2**2)*x14))))*
     &                     x30 + 
     &                    (3 + c2tt)*MSt2diff*pSB2*x14*Xtr**2) - 
     &                 MSt22*x14*
     &                  (3*c2tt*(MSt2sum - 2*MUEr2)**2*x5 + 
     &                    MSt2diff*
     &                     (MSt2diff*
     &                       (8 + 4*Nc + 
     &                       (2 - 4*MUEr2*pMT2*x14)*x30 - c2tt*x5)-
     &                        2*
     &                       ((MSt2sum - 2*MUEr2)*x5 + 
     &                       pSB2*
     &                       (-x31 + c2tt2*(3*x31 - 3*x32) + x32)*
     &                       Xtr**2))))/(MSt22*x14) - 
     &              ((MSt2diff**2 + 
     &                   c2tt2*
     &                    (16*MSt2sum*MA02 - 8*MA02**2 - 
     &                      3*MSt2diff**2))*pCB2*x26*Ytr**2)/MSt12)
     &           )/MSt2diff**2) - 
     &      1/2.D0*((MSt2diff*
     &             (MSt2diff*(15 + 4*Nc + c2tt2*(8 + 8*Nc)) - 
     &               3*MA02*pCB2) - 
     &            MSt2sum*(2*MbSL2 - 12*MUEr2 + 
     &               MSt2sum*(15 + 4*Nc + c2tt2*(8 + 8*Nc)) + 
     &               3*MA02*pCB2 + 2*pMT2) + 
     &            c2tt*(2*(MbSL2*MSt2diff + MA02*MSt22*pCB2) - 
     &               2*MSt2diff*pMT2))*s2tt)/(MSt12*MSt22) + 
     &      (2*(6*(Atr*(7 - 15*c2tt2)*pMT - MSt2sum*c2tt*s2tt) - 
     &           1/2.D0*(s2tt*
     &               (8*MUEr2*pMT2 + 
     &                 MSt2diff*(2*MSt12 - 2*(MUEr2 + pMT2))*
     &                  (1 - 2*MUEr2*pMT2*x13))*x27)/MSt12))/
     &       MSt2diff + (2*
     &         (-2*s2tt*(x27/(MSt12*x13) - 
     &              c2tt*(2*Atr**2 - 3*MbSL2 + 9*MUEr2 - 
     &                 3*MA02*pCB2 + 6*pMT2 - 
     &                 8*(pSB2*Xtr**2 + pCB2*Ytr**2))) + 
     &           (MSt2sum*c2tt*MUEr2*s2tt - 
     &              3/2.D0*(Atr*(3*MSt2sum + 5*MSt2diff)*pMT*
     &                 s2tt2) - 
     &              1/4.D0*((2*(MA02 + MbSL2) - 2*MSt22)*pCB2*x16*
     &                 ((1 - 3*c2tt)*pMT2*s2tt + 
     &                   Ytr*
     &                    (pMT*(4*c2tt2 - 2*s2tt2) + 
     &                      s2tt*(Ytr + 3*c2tt*Ytr)))))/MSt22))/
     &       MSt2diff)

        lc21 = lc21 + 1/2.D0*
     &     ((s2tt*(MSt2diff**4*(2 - s2tt2 - Nc*s2tt2) + 
     &            8*c2tt2*(2*MSt2sum - MA02)*MA02*MSt12*pCB2*x28*
     &             Ytr**2 + 
     &            MSt2diff**2*
     &             ((3 - c2tt)*MSt22*(pSB2*Xtr**2 + pCB2*Ytr**2) + 
     &               MA02*pCB2*
     &                (2*MSt22*
     &                   (MSt2sum*(1 + c2tt2)*x22*x26*Ytr**2 + 
     &                     MbSL2*x15*x20*
     &                      (pMT2 + 2*pMT*s2tt*Ytr + Ytr**2 + 
     &                       c2tt*(pMT2 - Ytr**2))) + 
     &                  MSt12*
     &                   (3 + 
     &                     c2tt*
     &                      (1 - 2*MbSL2*x16*x19*(pMT2 - Ytr**2))+
     &                       2*
     &                      (MSt2sum*(1 + c2tt2)*x21*x28*Ytr**2 + 
     &                       MbSL2*x16*x19*
     &                       (pMT2 - 2*pMT*s2tt*Ytr + Ytr**2)))) + 
     &               1/2.D0*(pCB2*
     &                  (2*MSt12*(3 + c2tt + x28 - 3*c2tt2*x28) - 
     &                    (1 + c2tt2)*MA02**2*
     &                     (MSt2diff*(x22*x26 - x21*x28) + 
     &                       MSt2sum*(x22*x26 + x21*x28)))*Ytr**2))
     &            ))/(MSt12*MSt22*MSt2diff**2) + 
     &       2*((MA02*pCB2*x25*
     &              ((2*Ytr*
     &                   (pMT*
     &                      (-(MSt12*(8 - 12*s2tt2)) + 
     &                       MA02*(2 - 3*s2tt2)) + 
     &                     c2tt2*
     &                      (-4*MSt2sum + 2*MA02 + 3*MSt2diff)*
     &                      s2tt*Ytr))/MSt2diff**2 - 
     &                (4*c2tt2*pMT*Ytr - 
     &                   1/2.D0*(s2tt*(2*pMT + s2tt*Ytr)**2))/
     &                 MSt2diff + 
     &                1/2.D0*
     &                 (s2tt*(2*pMT + s2tt*Ytr)**2)/
     &                  (MA02 - 4*MSt12)) + 
     &             (MSt2sum*c2tt*MUEr2*s2tt - 
     &                3/2.D0*
     &                 (Atr*(3*MSt2sum - 5*MSt2diff)*pMT*s2tt2) + 
     &                1/4.D0*
     &                 ((2*(MA02 + MbSL2) - 2*MSt12)*pCB2*x15*
     &                   (pMT2*(s2tt + 3*c2tt*s2tt) - 
     &                     Ytr*
     &                      (pMT*(4*c2tt2 - 2*s2tt2) + 
     &                       (-1 + 3*c2tt)*s2tt*Ytr))))/MSt2diff)/
     &           MSt12 + (pCB2*
     &             (MA02*x29*
     &                ((-2*Ytr*
     &                     (pMT*
     &                       (-(MSt22*(8 - 12*s2tt2)) + 
     &                       MA02*(2 - 3*s2tt2)) + 
     &                       c2tt2*
     &                       (4*MSt2sum - 2*MA02 + 3*MSt2diff)*
     &                       s2tt*Ytr))/MSt2diff**2 - 
     &                  (4*c2tt2*pMT*Ytr + 
     &                     1/2.D0*(s2tt*(2*pMT - s2tt*Ytr)**2))/
     &                   MSt2diff + 
     &                  1/2.D0*
     &                   (s2tt*(2*pMT - s2tt*Ytr)**2)/
     &                    (MA02 - 4*MSt22)) + 
     &               1/2.D0*((2*MSt22*x15*x19 - 2*MSt12*x16*x20)*
     &                   (2*pMT*(2 - 3*s2tt2)*Ytr - 
     &                     3*c2tt*s2tt*(pMT2 - Ytr**2)))/
     &                 (MSt12*MSt2diff**2*x19*x20)))/MSt22))

        lc3 = -(1/2.D0*(pCB2*
     &         (MSt2diff*s2tt*(pMT2 + Ytr**2) - 
     &           (MSt2sum - 2*(MA02 + MbSL2))*
     &            (2*pMT*(2 - 3*s2tt2)*Ytr - 
     &              3*c2tt*s2tt*(pMT2 - Ytr**2))))/MSt2diff**2)

        lc4 = -(1/2.D0*((MSt2diff - c2tt2*(8*MA02 + 3*MSt2diff))*
     &         pCB2*s2tt*Ytr**2)/MSt2diff**2)

        lc5 = 1/2.D0*((c2tt2*(8*MA02 - 3*MSt2diff) + MSt2diff)*
     &       pCB2*s2tt*Ytr**2)/MSt2diff**2

        lc6 = (MA02*pCB2*
     &      (3*MSt2diff*s2tt - 8*pMT*(2 - 3*s2tt2)*Ytr - 
     &        3*c2tt*s2tt*(MSt2sum - 2*pMT2 + 2*Ytr**2)))/
     &    MSt2diff**2

        lc7 = (MbSL2*(4*Atr*pMT*(2 - 3*s2tt2) - 
     &        s2tt*(MSt2diff - 
     &           3*c2tt*(MSt2sum - 2*pMT2 + 
     &              2*(pSB2*Xtr**2 + pCB2*Ytr**2)))))/MSt2diff**2

        lc8 = -(1/2.D0*((MSt2diff*
     &             (MSt2sum*(2 + Nc) + 3*MA02*pCB2) + 
     &            3*c2tt*(2*MUEr2**2 + 
     &               MA02*pCB2*(-MSt2sum + 2*pMT2)))*s2tt - 
     &         4*pMT*(2 - 3*s2tt2)*(Atr*MbSL2 + 2*MA02*pCB2*Ytr) - 
     &         s2tt*(6*c2tt*MA02*pCB2*Ytr**2 - 
     &            MbSL2*(MSt2diff - 
     &               3*c2tt*
     &                (MSt2sum - 2*pMT2 + 
     &                  2*(pSB2*Xtr**2 + pCB2*Ytr**2)))))/
     &       MSt2diff**2)

	lc9 = ((MSt2sum - 2*(MUEr2 + pMT2))*s2tt)/MSt2diff

        s2ttF3 = -lc21 + lc13*Log(MA02/Q) + 
     &    lc14*Log(MbSL2/Q) + lc15*Log(MSt12/Q)**2 - 
     &    lc18*Log(MSt22/Q) + lc17*Log(MSt22/Q)**2 - 
     &    Log(MSt12/Q)*(lc16 - lc12*Log(MSt22/Q)) - 
     &    lc19*Log(MUEr2/Q) - lc20*Log(pMT2/Q) - 
     &    Log(MSt22/MSt12)*
     &     (lc6*Log(MA02) - lc7*Log(MbSL2) - lc4*Log(MA02/MSt12) + 
     &       lc5*Log(MA02/MSt22) - 
     &       lc3*Log((MA02*MbSL2)/(MSt12*MSt22)) + 
     &       lc8*Log(MSt12*MSt22) - lc10*Log(MUEr2) - 
     &       lc9*Log((MSt12*MSt22)/(MUEr2*pMT2)) + 
     &       lc11*(1 + Log(Q)))

	tmp21 = -2*x10*x11 + x11**2 + 2*x10*x9 - x9**2

	dMUErF2 = 1/4.D0*(Nc*tmp21)

        dMUErF3 = -(1/4.D0*(Nc*(-2 - 2*x10 + x11 + x9)*
     &        (2*MSt2diff + MSt2sum*(-x11 + x9)))/MSt2diff)

	s2ttAtrdAtrF2 = 1/2.D0*(Atr*(3 + Nc)*s2tt*tmp21)

        s2ttAtrdAtrF3 = 1/2.D0*
     &   (Atr*(3 + Nc)*s2tt*(2 + 2*x10 - x11 - x9)*
     &      (2*MSt2diff + MSt2sum*(-x11 + x9)))/MSt2diff

        dPhi(se12) = 1/8.D0*(MUEr*Nc*
     &      (2*MSt2diff*pMT*s2tt*
     &         (-2*x10*x11 + x11**2 + 2*x10*x9 - x9**2) + 
     &        Atr*s2tt2*(2 + 2*x10 - x11 - x9)*
     &         (2*MSt2diff + MSt2sum*(-x11 + x9))))/MSt2diff

        dPhi(se22) = Nc*((pMT2 - Atr*pMT*s2tt)*(x10 - x11)**2 - 
     &     2*pMT2*(x10 - x7)**2 + 
     &     (pMT2 + Atr*pMT*s2tt)*(x10 - x9)**2 + 
     &     1/4.D0*(Atr**2*s2tt2*(-2 - 2*x10 + x11 + x9)*
     &         (2*MSt2diff + MSt2sum*(-x11 + x9)))/MSt2diff)

* COUNTER TERMS:

	if( onshell ) then
	  F1o = x11 - 2*x7 + x9
	  F2o = -x11 + x9
	  F3o = 2 + (MSt2sum*(-x11 + x9))/MSt2diff

	  dvev2vev2 = 1/2.D0*
     &     (Nc*pSB2*(-((MbSL2 + ctt2*MSt12 + pMT2 + MSt22*stt2)/
     &          pMT2) + 2*(-x10 + x7 + 
     &          (MbSL2*((MSt22*stt2*(x11 - x12))/
     &                (-MbSL2 + MSt22) + 
     &               (ctt2*MSt12*(x12 - x9))/(MbSL2 - MSt12)))/
     &           pMT2)))

	  dpMTpMT = 1/4.D0*(MSt2sum - 2*pMT2 - 3*MUEr2 + 10*pMT2*pSB2 + 
     &      MSt12*x10 + MSt22*x10 - 2*pMT2*x10 - 3*MUEr2*x10 + 
     &      5*pMT2*pSB2*x10 - MSt22*x11 + 
     &      pCB2*(MA02*(3 + 3*x10 - 3*x18 - 2*x33 - x34) + 
     &         pMT2*(4*x33 + x34)) + 
     &      MbSL2*(1 + x10 - x12 - x35) + pMT2*x35 + MUEr2*x35 - 
     &      MSt12*x36 + pMT2*x36 + MUEr2*x36 - MSt22*x37 + pMT2*x37 + 
     &      MUEr2*x37 + 2*pMT2*x7 - 5*pMT2*pSB2*x7 + 3*MUEr2*x8 - 
     &      MSt12*x9)/pMT2

	  dMSt12 = -x39 - (1 + stt2)*x40 + 
     &      stt2*(x41 + (MSt12 - MUEr2)*x42) + pSB2*tmp13*x43 + 
     &      pCB2*(x44 + stt2*x44 + tmp16*x45) + 
     &      (MSt12 - pMT2 - MUEr2)*x51 - 
     &      1/2.D0*((-2 + s2tt2 + Nc*s2tt2)*x38) + 
     &      1/2.D0*((1 + Nc)*s2tt2*x46 + 
     &        pSB2*((1 + c2tt2)*x48*Xtr**2 + 
     &           x47*(2*pMT + s2tt*Xtr)**2) + 
     &        pCB2*((1 + c2tt2)*x50*Ytr**2 + x49*(2*pMT + s2tt*Ytr)**2))

	  dMSt22 = -x39 - (1 + ctt2)*x40 + 
     &      ctt2*(x41 + (MSt22 - MUEr2)*x52) + pSB2*tmp11*x53 + 
     &      pCB2*(x44 + ctt2*x44 + tmp14*x54) + 
     &      (MSt22 - pMT2 - MUEr2)*x59 - 
     &      1/2.D0*((-2 + s2tt2 + Nc*s2tt2)*x46) + 
     &      1/2.D0*((1 + Nc)*s2tt2*x38 + pSB2*x57*(-2*pMT + s2tt*Xtr)**2 + 
     &        pCB2*x58*(-2*pMT + s2tt*Ytr)**2 + 
     &        (1 + c2tt2)*(pSB2*x55*Xtr**2 + pCB2*x56*Ytr**2))

	  tmp12 = pMT2*s2tt - Xtr*(2*c2tt*pMT + s2tt*Xtr)

	  tmp15 = pMT2*s2tt - Ytr*(2*c2tt*pMT + s2tt*Ytr)

	  ds2tt = 1/2.D0*(c2tt*
     &       (s2tt*(2*x40 - 2*x41 - MSt12*x42 + MUEr2*x42 - 
     &            2*pCB2*x44 - MSt22*x52 + MUEr2*x52) + 
     &         pSB2*tmp12*(x43 + x53) + pCB2*tmp15*(x45 + x54)) + 
     &      c2tt2*(-2*pMT*(pSB2*(x47 + x48 + x55 + x57)*Xtr + 
     &            pCB2*(x49 + x50 + x56 + x58)*Ytr) + 
     &         s2tt*(2*(1 + Nc)*x38 - 2*(1 + Nc)*x46 - 
     &            pSB2*Xtr**2*(x47 - x48 + x55 - x57) -
     &            pCB2*Ytr**2*(x49 - x50 + x56 - x58))))/MSt2diff

	  k = s2tt*(dMSt12/MSt12 - dMSt22/MSt22)

	  F1 = F1 + dMSt12/MSt12 + dMSt22/MSt22 -
     &      4*dpMTpMT*(1 - F1o) - dvev2vev2*F1o

	  s2ttF2 = s2ttF2 + k +
     &      (ds2tt + s2tt*(3*dpMTpMT - dvev2vev2))*F2o

          s2ttF3 = s2ttF3 +
     &      (s2tt*(2*dpMTpMT - dvev2vev2) + 2*ds2tt)*F3o +
     &      k*(2*MSt12*MSt22/MSt2diff*F2o + MSt2sum)/MSt2diff

	  s2ttdAtr = Xtr*(ds2tt +
     &      s2tt*((dMSt22 - dMSt12)/MSt2diff - dpMTpMT))
	  s2ttAtrdAtrF2 = s2ttAtrdAtrF2 + s2ttdAtr*F2o
	  s2ttAtrdAtrF3 = s2ttAtrdAtrF3 + s2ttdAtr*F3o
	endif

	k = Nc/(16*pi**2)**2*ht2**2

	dPhi(se11) = k*.5D0*MUEr2*s2tt*(s2ttF3 + 2*s2tt*dMUErF3)

	dPhi(se12) = k*(dPhi(se12) - MUEr*
     &    (.5D0*s2tt*(Atr*(s2ttF3 + s2tt*dMUErF3) + s2ttAtrdAtrF3) +
     &     pMT*(s2ttF2 + s2tt*dMUErF2)))

	dPhi(se22) = k*(dPhi(se22) +
     &    s2tt*Atr*(.5D0*Atr*s2ttF3 + s2ttAtrdAtrF3) +
     &    2*pMT*(Atr*s2ttF2 + s2ttAtrdAtrF2) +
     &    2*pMT2*F1)

	if( debuglevel .lt. 6 ) return

	DHIGGS "atat: F1         =", F1		ENDL
	DHIGGS "atat: s2ttF2     =", s2ttF2	ENDL
	DHIGGS "atat: s2ttF3     =", s2ttF3	ENDL
	DHIGGS "atat:   se11     =", dPhi(se11)	ENDL
	DHIGGS "atat:   se12     =", dPhi(se12)	ENDL
	DHIGGS "atat:   se22     =", dPhi(se22)	ENDL
	end
