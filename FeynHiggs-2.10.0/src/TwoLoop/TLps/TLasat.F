* TLasat.F
* Two-loop O(a_s a_t) corrections to the CP-even Higgs mass matrix.
* Routine written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
* Based on G. Degrassi, P. Slavich, and F. Zwirner,
* Nucl. Phys. B611 (2001) 403 [hep-ph/0105096].
* this file is part of FeynHiggs
* last modified 4 Dec 13 th

* Revision history:
* 14/07/2004: rewritten completely by TH for FeynHiggs
* 24/02/2004: mg is given as input instead of mg^2;
*             value of pi corrected (10th digit);
*             unused variables cleaned up.
* 22/10/2002: gs is given as input.

#include "externals.h"
#include "types.h"
#include "debug.h"


* dPhi(1..3) = corrections to the Higgs mass matrix
* where 1 = 11, 2 = 12, 3 = 22

	subroutine TLasat(dPhi)
	implicit none
	RealType dPhi(3)

#include "TLps.h"

	RealType F1, s2ttF2, s2ttF3, F1o, F2o, F3o
	RealType s2ttAtrdF2, s2ttAtrdF3, s2ttdAtr
	RealType dpMTpMT, dMSt12, dMSt22, ds2tt, k

	RealType x1, x6, x3, x4, x5, x7, x8, x9, x10, x13, x14
	RealType x11, x12
	RealType tmp1, tmp2, tmp3, tmp4, tmp5

	x1 = InvTri(pMGl2,pMT2,MSt12)

	x6 = InvTri(pMGl2,pMT2,MSt22)

	x3 = Log(pMGl2)

	x4 = Log(MSt12)

	x5 = Log(MSt22)

	x7 = Log(pMT2)

	x8 = Log(Q)

	x9 = Phi(pMT2,MSt12,pMGl2)

	x10 = Phi(pMT2,MSt22,pMGl2)

	x13 = myB0(MSt12,pMT2,pMGl2,Q)

	x14 = myB0(MSt22,pMT2,pMGl2,Q)

	x11 = myB0(pMT2,pMGl2,MSt12,Q)

	x12 = myB0(pMT2,pMGl2,MSt22,Q)

	tmp1 = pMGl2 - MSt12 + pMT2 - 2*pMGl*pMT*s2tt

	tmp2 = pMGl2 - MSt22 + pMT2 + 2*pMGl*pMT*s2tt

        F1 = -6 + (4*(pMGl2 + pMT2 - pMGl*pMT*s2tt) - MSt22*s2tt2)/
     -    MSt12 + (-2 + (4*pMGl2*
     -         (2*pMGl2*pMT - pMGl*(pMGl2 - MSt12 + pMT2)*s2tt)*x1)/pMT)*
     -    (-x3 + x4) + c2tt2*(-2 + (x4 - x5)**2) + 
     -   (-x3 + x5)*(-2 + 4*pMGl2*
     -       (2*pMGl2 + (pMGl*(pMGl2 - MSt22 + pMT2)*s2tt)/pMT)*x6) + 
     -   (5 + 2*x5 - 2*x7)*(x4 - x7) + (x4 - x7)**2 + 
     -   5*(x5 - x7) + (x5 - x7)**2 + 
     -   x6*(4*pMGl2*tmp2*x10 + 
     -      ((2*pMGl2**2*MSt22*pMT + 
     -           4*pMGl2*pMT2*
     -            (pMT*(5*MSt22 + pMT2) + pMGl*(3*MSt22 + pMT2)*s2tt)-
     -             2*(MSt22 - pMT2)**2*
     -            (MSt22*pMT + 2*pMT2*(pMT + pMGl*s2tt)))*(-x3 + x7))/
     -       (MSt22*pMT)) + 
     -   (4*(pMGl2 + pMT2 + pMGl*pMT*s2tt) - MSt12*s2tt2 + 
     -      ((-4*MSt2sum*(pMGl2 + pMT2) + 4*pMGl*MSt2diff*pMT*s2tt)*
     -         (x3 - x8))/MSt12)/MSt22 + 
     -   8*(-x3 + x7 + (x7 - x8)**2) + 
     -   (1 - (MSt2diff*s2tt2)/MSt22 + (4*pMGl*s2tt*(x7 - x8))/pMT)*
     -    (x4 - x8) + (-4 + (2*pMGl*s2tt)/pMT)*(x4 - x8)**2 + 
     -   (1 + (MSt2diff*s2tt2)/MSt12 - (4*pMGl*s2tt*(x7 - x8))/pMT)*
     -    (x5 - x8) + (-4 - (2*pMGl*s2tt)/pMT)*(x5 - x8)**2 + 
     -   6*(x7 - x8) + (4*pMGl*s2tt*(x10 - x9))/pMT + 
     -   x1*(((2*pMGl2**2*MSt12*pMT + 
     -           4*pMGl2*pMT2*
     -            (pMT*(5*MSt12 + pMT2) - pMGl*(3*MSt12 + pMT2)*s2tt)-
     -             2*(MSt12 - pMT2)**2*
     -            (MSt12*pMT + 2*pMT2*(pMT - pMGl*s2tt)))*(-x3 + x7))/
     -       (MSt12*pMT) + 4*pMGl2*tmp1*x9)

        s2ttF2 = (-2*s2tt + 
     -      4*pMGl2*((pMGl2 - MSt12 + pMT2)*s2tt - 2*pMGl*pMT*s2tt2)*x1)
     -     *(-x3 + x4) + (-x3 + x5)*
     -    (2*s2tt - 4*pMGl2*
     -       ((pMGl2 - MSt22 + pMT2)*s2tt + 2*pMGl*pMT*s2tt2)*x6) + 
     -   x6*(-(((8*pMGl*pMGl2*pMT*pMT2*s2tt + 
     -             (pMGl2 - MSt22 + pMT2)**3*s2tt2)*x10)/(pMGl*pMT)) - 
     -      (4*(-((MSt22 - pMT2)**2*(pMT2*s2tt + pMGl*pMT*s2tt2)) + 
     -           pMGl2*(pMT2*(3*MSt22 + pMT2)*s2tt + 
     -              pMGl*pMT*(MSt22 + pMT2)*s2tt2))*(-x3 + x7))/MSt22)
     -     + ((4*MSt2diff*pMT*(pMGl2 + pMT2)*s2tt - 
     -         (4*pMGl*(3*MSt12*MSt22 + MSt2sum*pMT2) + 
     -            MSt2diff*MSt2sum*pMT*s2tt)*s2tt2 + 
     -         ((4*MSt2diff*
     -               (-(MSt2diff*pMT*(pMGl2 + pMT2)*s2tt) + 
     -                 pMGl*MSt2sum*pMT2*s2tt2) + 
     -              pMGl*MSt12*MSt22*(2*pMGl2 + MSt2sum - 2*pMT2)*
     -               s2tt2*(x4 - x5))*(x3 - x8))/MSt2diff)/MSt22 + 
     -      (((2*pMGl*MSt12*(4*MSt22 + MSt2sum) + 
     -              MSt22*MSt2diff*pMT*s2tt)*s2tt2 + 
     -           pMGl*MSt12*
     -            (-8*c2tt2*pMT2 + 
     -              (8*pMT2 - (-2*pMGl2 + MSt2sum + 10*pMT2)*s2tt2)*
     -               (x7 - x8)))*(x5 - x8))/MSt2diff)/(MSt12*pMT) + 
     -   (-3*s2tt + (pMGl*MSt2sum*s2tt2)/(MSt2diff*pMT))*
     -    (x4 - x8)**2 + (3*s2tt - 
     -      (pMGl*MSt2sum*s2tt2)/(MSt2diff*pMT))*(x5 - x8)**2 + 
     -   x1*((4*((MSt12 - pMT2)**2*(-(pMT2*s2tt) + pMGl*pMT*s2tt2) + 
     -           pMGl2*(pMT2*(3*MSt12 + pMT2)*s2tt - 
     -              pMGl*pMT*(MSt12 + pMT2)*s2tt2))*(-x3 + x7))/MSt12+
     -        ((8*pMGl*pMGl2*pMT*pMT2*s2tt - 
     -           (pMGl2 - MSt12 + pMT2)**3*s2tt2)*x9)/(pMGl*pMT)) + 
     -   (s2tt*(x4 - x5)*(5*(1 + c2tt2)*MSt2diff + 
     -         c2tt2*(MSt2sum*x4 - 2*MSt22*x5 - MSt2sum*x5 + 
     -            2*MSt12*(x4 - x8) + 2*MSt22*x8)) + 
     -      (((-((2*pMGl*MSt22*(4*MSt12 + MSt2sum) + 
     -                   MSt12*MSt2diff*pMT*s2tt)*s2tt2) + 
     -              pMGl*MSt22*
     -               (8*c2tt2*pMT2 + 
     -                 (-8*pMT2 + 
     -                    (-2*pMGl2 + MSt2sum + 10*pMT2)*s2tt2)*
     -                  (x7 - x8)))*(x4 - x8))/MSt22 + 
     -         (8*c2tt2*pMGl2*pMT2*(-x10 + x9) + 
     -            s2tt2*((-2*x10)/x6 + (2*x9)/x1))/pMGl)/pMT)/
     -    MSt2diff

	tmp5 = -((x4 - x5)*(x4 + x5 - 2*x8))

	s2ttAtrdF2 = 2*pMGl*s2tt*tmp5

        tmp3 = MSt2diff*(-2*MSt2sum + 4*(pMGl2 + pMT2))*s2tt + 
     -   4*pMGl*MSt2sum*pMT*(2 - 3*s2tt2)

        tmp4 = MSt2diff*(MSt2sum - 2*(pMGl2 + pMT2))*s2tt + 
     -   2*pMGl*MSt2sum*pMT*(-2 + 3*s2tt2)

        s2ttF3 = s2tt*(10 + 18*c2tt2 - 
     -      4*(c2tt2 + (MSt12*(pMGl2 + pMT2))/(MSt22*MSt2diff)) + 
     -      ((4*MSt22*(pMGl2 + pMT2))/MSt2diff - 
     -         ((MSt12**2 + MSt22**2)*s2tt2)/MSt22)/MSt12) + 
     -   ((-2 + (2*(MSt22 - pMT2))/pMGl2 - (8*pMT2)/MSt2diff)*s2tt - 
     -      (4*pMT*(pMGl2 - MSt22 + pMT2)*(-2 + 3*s2tt2))/
     -       (pMGl*MSt2diff))*x10 + 
     -   (pMGl*pMT*(64 + (-88 - (4*MSt12)/MSt22 - (4*MSt22)/MSt12)*
     -          s2tt2) - 3*MSt2sum*s2tt*(x4 - x5)*(x4 + x5 - 2*x8))
     -     /MSt2diff + x6*(4*pMT2*
     -       (s2tt*tmp2*x10 + 
     -         (((-(MSt22 - pMT2)**2 + pMGl2*(MSt22 + pMT2))*s2tt + 
     -              pMGl*pMT*(pMGl2 + MSt22 - pMT2)*s2tt2)*(x7 - x8))/
     -          MSt22) + (8*pMGl2*pMT2*s2tt + 
     -         4*pMGl*pMT*(pMGl2 - MSt22 + pMT2)*s2tt2)*(x5 - x8)) + 
     -   (x3 - x8)*(pMGl2*(((-4*MSt2sum*s2tt)/MSt22 + 
     -            4*(-((pMGl2 + MSt12 - pMT2)*pMT2*s2tt) - 
     -               pMGl*pMT*(-pMGl2 + MSt12 + pMT2)*s2tt2)*x1)/MSt12+
     -           (4*(-((pMGl2 + MSt22 - pMT2)*pMT2*s2tt) + 
     -              pMGl*pMT*(-pMGl2 + MSt22 + pMT2)*s2tt2)*x6)/MSt22)+
     -        (2*tmp4*x4 + tmp3*x5 - tmp3*x8 - 2*tmp4*x8)/
     -       MSt2diff**2) + 
     -   ((-2 + (2*(MSt12 - pMT2))/pMGl2 + (8*pMT2)/MSt2diff)*s2tt + 
     -      (pMT*(pMGl2 - MSt12 + pMT2)*(8 - 12*s2tt2))/(pMGl*MSt2diff)
     -      )*x9 + x1*((8*pMGl2*pMT2*s2tt - 
     -         4*pMGl*pMT*(pMGl2 - MSt12 + pMT2)*s2tt2)*(x4 - x8) + 
     -      4*pMT2*((((-(MSt12 - pMT2)**2 + pMGl2*(MSt12 + pMT2))*
     -               s2tt - pMGl*pMT*(pMGl2 + MSt12 - pMT2)*s2tt2)*
     -            (x7 - x8))/MSt12 + s2tt*tmp1*x9)) + 
     -   ((x5 - x8)*(-8*pMGl2*MSt2diff*s2tt - 
     -         3*MSt22*MSt2diff*s2tt - 
     -         13*c2tt2*MSt22*MSt2diff*s2tt - 
     -         8*MSt2diff*pMT2*s2tt + 
     -         (MSt22**2*MSt2diff*s2tt*s2tt2)/MSt12 - 
     -         2*MSt2diff**2*s2tt*s2tt2 + 
     -         8*pMGl*pMT*(-2 + 3*s2tt2)*
     -          (2*MSt22 + MSt2sum - (pMGl2 - pMT2)*(x3 - x7)) + 
     -         tmp3*x7 - tmp3*x8) + 
     -      (x4 - x8)*(8*pMGl2*MSt2diff*s2tt + 
     -         3*MSt12*MSt2diff*s2tt + 
     -         13*c2tt2*MSt12*MSt2diff*s2tt + 
     -         8*MSt2diff*pMT2*s2tt - 
     -         (MSt12**2*MSt2diff*s2tt*s2tt2)/MSt22 - 
     -         2*MSt2diff**2*s2tt*s2tt2 - 
     -         8*pMGl*pMT*(-2 + 3*s2tt2)*
     -          (2*MSt12 + MSt2sum - (pMGl2 - pMT2)*(x3 - x7)) + 
     -         2*tmp4*x7 - 2*tmp4*x8) + 
     -      MSt2sum*s2tt*(x4 - x5)*
     -       ((3 + 9*c2tt2)*MSt2diff + 
     -         c2tt2*(MSt2sum*x4 - 6*MSt22*x5 - MSt2sum*x5 + 
     -            6*MSt12*(x4 - x8) + 6*MSt22*x8)) + 
     -      ((pMGl*MSt2diff*s2tt + 2*pMGl2*pMT*(-2 + 3*s2tt2))*
     -         ((-4*x10)/x6 + (4*x9)/x1))/(pMGl*pMGl2))/MSt2diff**2

        s2ttAtrdF3 = (2*pMGl*s2tt*
     -     (4*MSt2diff + MSt2sum*tmp5 + 4*MSt12*x4 - 4*MSt22*x5 - 
     -       4*MSt12*x8 + 4*MSt22*x8))/MSt2diff

* COUNTER TERMS:

	if( onshell ) then
	  F1o = x4 + x5 - 2*x7
	  F2o = x4 - x5
	  F3o = 2 + (MSt2sum*(x4 - x5))/MSt2diff

	  dpMTpMT = 1/2.D0*(MSt2sum + 2*MSDR*pMT2 + tmp1*x11 + 
     -      tmp2*x12 - MSt12*x4 - MSt22*x5 + 6*pMT2*x7 + 
     -      2*pMGl2*(-1 + x3 - x8) + MSt12*x8 + MSt22*x8 - 6*pMT2*x8)
     -     /pMT2

	  dMSt12 = -2*(pMGl2 + 3*MSt12 + pMT2 + tmp1*x13) + 
     -     2*(pMGl2*(x3 - x8) + pMT2*(x7 - x8)) - 
     -     (-3 + c2tt2)*MSt12*(x4 - x8) + 
     -     s2tt2*(MSt2diff + MSt22*(-x5 + x8))

	  dMSt22 = -(MSt2diff*s2tt2) - 
     -     2*(pMGl2 + 3*MSt22 + pMT2 + tmp2*x14) + 2*pMGl2*(x3 - x8) - 
     -     MSt12*s2tt2*(x4 - x8) - (-3 + c2tt2)*MSt22*(x5 - x8) + 
     -     2*pMT2*(x7 - x8)

	  ds2tt = (-2*c2tt2*
     -     (2*pMGl*pMT*(x13 + x14) + 
     -       s2tt*(MSt2diff + MSt12*(x4 - x8) - MSt22*(x5 - x8))))
     -      /MSt2diff

	  k = s2tt*(dMSt12/MSt12 - dMSt22/MSt22)

* Eq. (37)
	  F1 = F1 + dMSt12/MSt12 + dMSt22/MSt22 - 4*dpMTpMT*(1 - F1o)

* Eq. (38)
	  s2ttF2 = s2ttF2 + k + (3*s2tt*dpMTpMT + ds2tt)*F2o

* Eq. (39)
	  s2ttF3 = s2ttF3 + 2*(s2tt*dpMTpMT + ds2tt)*F3o +
     &      k*(2*MSt12*MSt22/MSt2diff*F2o + MSt2sum)/MSt2diff

* Eq. (B8)
	  s2ttdAtr = Xtr*(ds2tt +
     &      s2tt*((dMSt22 - dMSt12)/MSt2diff - dpMTpMT))
	  s2ttAtrdF2 = s2ttAtrdF2 + s2ttdAtr*F2o
	  s2ttAtrdF3 = s2ttAtrdF3 + s2ttdAtr*F3o
	endif

	k = C_F*Nc/(64*pi**3)*AlfasMT*ht2

* Eq. (25)
	dPhi(se11) = k*.5D0*MUEr2*s2tt*s2ttF3

* Eq. (26)
	dPhi(se12) = -k*MUEr*(.5D0*s2tt*(Atr*s2ttF3 + s2ttAtrdF3) +
     &    pMT*s2ttF2)

* Eq. (27)
	dPhi(se22) = k*(s2tt*Atr*(.5D0*Atr*s2ttF3 + s2ttAtrdF3) +
     &    2*pMT*(Atr*s2ttF2 + s2ttAtrdF2) + 2*pMT2*F1)

	if( debuglevel .lt. 6 ) return

	DHIGGS "asat: F1         =", F1		ENDL
	DHIGGS "asat: s2ttF2     =", s2ttF2	ENDL
	DHIGGS "asat: s2ttF3     =", s2ttF3	ENDL
	DHIGGS "asat:   se11     =", dPhi(se11)	ENDL
	DHIGGS "asat:   se12     =", dPhi(se12)	ENDL
	DHIGGS "asat:   se22     =", dPhi(se22)	ENDL
	end
