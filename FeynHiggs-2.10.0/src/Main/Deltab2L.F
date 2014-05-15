* Deltab2L.F
* 2-loop parts of Delta b
* code adapted from HDecay
* this file is part of FeynHiggs
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	RealType function Deltab2L_ew(scale)
	implicit none
	RealType scale

#include "FH.h"

	ComplexType Ieps
	parameter (Ieps = (1D0, 1D-15))

	RealType MSu12, MSu22, MSd12, MSd22
	RealType MSt12, MSt22, MSb12, MSb22
	equivalence (MSu12, MSf2(1,3,2)), (MSu22, MSf2(2,3,2))
	equivalence (MSt12, MSf2(1,3,3)), (MSt22, MSf2(2,3,3))
	equivalence (MSd12, MSdL2(1,2)), (MSd22, MSdL2(2,2))
	equivalence (MSb12, MSdL2(1,3)), (MSb22, MSdL2(2,3))

	RealType MUE2
	RealType fnorm, anomalous, finscale
	ComplexType d, yt1, yt2, XGl2, XT2
	ComplexType XSt12, XSt22, XSb12, XSb22
	ComplexType R22, R42, R52, R72, R82, R92
	ComplexType Rat2, Rlt2, RMSt12, RMSt22
	RealType Rew

	RealType T134s, T134p, TLdb
	external T134s, T134p, TLdb

	MUE2 = Sq(MUE)

	fnorm = 1/MUE2
	d = fnorm*Ieps
	XT2 = MT2*d
	XGl2 = MGl2*d
	XSt12 = MSt12*d
	XSt22 = MSt22*d
	XSb12 = MSb12*d
	XSb22 = MSb22*d

	yt1 = sqrt(XGl2*(XGl2 - 2*(XT2 + XSt12)) + (XT2 - XSt12)**2)
	yt2 = sqrt(XGl2*(XGl2 - 2*(XT2 + XSt22)) + (XT2 - XSt22)**2)

        R22 = -(1/2.D0*(2*(XSt12 - XSt22)*(-1 + XSt12*XSt22)*
     &         (7 + log(XGl2)*(6 + 2*log(XGl2))) + 
     &        (-1 + XSt22)**2*
     &         (XSt12*(2*XSt12 + 4*(1 + log(XGl2)) - log(XSt12))*
     &            log(XSt12) + 
     &           XGl2*(1 + XSt12)*T134p(MUE2,MSt12,MGl2)) + 
     &        (-1 + XSt12)**2*
     &         (XSt22*log(XSt22)*
     &            (-2*(XSt22 + 2*(1 + log(XGl2))) + log(XSt22)) - 
     &           XGl2*(1 + XSt22)*T134p(MUE2,MSt22,MGl2)))/
     &      ((-1 + XSt12)**2*(XSt12 - XSt22)*(-1 + XSt22)**2))

        R42 = 1/2.D0*(4*(XSt12 - XSt22)*(-1 + XSt12*XSt22)*
     &       (7 + log(XGl2)*(6 + 2*log(XGl2))) + 
     &      (-1 + XSt22)**2*
     &       (XSt12*(7 + 5*XSt12 + 8*log(XGl2) - 2*log(XSt12))*
     &          log(XSt12) + 
     &         2*XGl2*(1 + XSt12)*T134p(MSt12,MUE2,MGl2)) + 
     &      (-1 + XSt12)**2*
     &       (-(XSt22*(7 + 5*XSt22 + 8*log(XGl2) - 2*log(XSt22))*
     &            log(XSt22)) - 
     &         2*XGl2*(1 + XSt22)*T134p(MSt22,MUE2,MGl2)))/
     &    ((-1 + XSt12)**2*(XSt12 - XSt22)*(-1 + XSt22)**2)

        R52 = 1/4.D0*((-1 + XSt22)*
     &       (XSt12*((4*XSb12*XSb22 + 
     &               XGl2*(2*XGl2 - 3*(XSb12 + XSb22)))*log(XGl2)+
     &              XSb12*(XGl2 - XSb22)*log(XSb12) - 
     &            (XGl2 - XSb12)*
     &             (XSb22*(4 - log(XSb22) - log(XSt12)) + 
     &               XGl2*(-4 + log(XSt12))))*log(XSt12) + 
     &         XGl2*(-((2*XGl2 - XSb12 - XSb22)*
     &               (XGl2 - XSt12 + XT2)*
     &               T134s(MSt12,MGl2,MT2,MGl2)) + 
     &            (XGl2 - XSb22)*(XSb12 - XSt12 + XT2)*
     &             T134s(MSt12,MSb12,MT2,MGl2) + 
     &            (XGl2 - XSb12)*(XSb22 - XSt12 + XT2)*
     &             T134s(MSt12,MSb22,MT2,MGl2))) + 
     &      (-1 + XSt12)*(XSt22*
     &          (((-4*XSb12*XSb22 + 
     &                  XGl2*(-2*XGl2 + 3*(XSb12 + XSb22)))*
     &                log(XGl2) + 
     &               (XGl2 - XSb22)*
     &                (-4*(XGl2 - XSb12) - XSb12*log(XSb12)) - 
     &               (XGl2 - XSb12)*XSb22*log(XSb22))*log(XSt22) + 
     &            (XGl2 - XSb12)*(XGl2 - XSb22)*log(XSt22)**2) + 
     &         XGl2*((2*XGl2 - XSb12 - XSb22)*(XGl2 - XSt22 + XT2)*
     &             T134s(MSt22,MGl2,MT2,MGl2) - 
     &            (XGl2 - XSb22)*(XSb12 - XSt22 + XT2)*
     &             T134s(MSt22,MSb12,MT2,MGl2) - 
     &            (XGl2 - XSb12)*(XSb22 - XSt22 + XT2)*
     &             T134s(MSt22,MSb22,MT2,MGl2))) + 
     &      XGl2*(XSt12 - XSt22)*
     &       (-((2*XGl2 - XSb12 - XSb22)*(-1 + XGl2 + XT2)*
     &            T134s(MUE2,MGl2,MT2,MGl2)) + 
     &         (XGl2 - XSb22)*(-1 + XSb12 + XT2)*
     &          T134s(MUE2,MSb12,MT2,MGl2) + 
     &         (XGl2 - XSb12)*(-1 + XSb22 + XT2)*
     &          T134s(MUE2,MSb22,MT2,MGl2)))/
     &    ((XGl2 - XSb12)*(XGl2 - XSb22)*(-1 + XSt12)*
     &      (XSt12 - XSt22)*(-1 + XSt22))

        R72 = 1/8.D0*((-1 + XSt22)**2*
     &       (4*XGl2*T134p(MSt12,MSt12,MGl2) + 
     &         XSt12*((5 - 8*XSt12)*log(XSt12)**2 + 
     &            4*XGl2*T134p(MSt12,MSt12,MGl2))) + 
     &      (-1 + XSt12)**2*
     &       (-4*XGl2*T134p(MSt22,MSt22,MGl2) + 
     &         XSt22*((-5 + 8*XSt22)*log(XSt22)**2 - 
     &            4*XGl2*T134p(MSt22,MSt22,MGl2))) - 
     &      4*((-1 + XSt12)**2*XSt22*
     &          (-2 - 3*log(XGl2) + XSt22*(3 + 4*log(XGl2)))*
     &          log(XSt22) + 
     &         XGl2*(1 + XSt12)*(-1 + XSt22)**2*
     &          T134p(MUE2,MSt12,MGl2)) + 
     &      4*((-1 + XSt12)*(XSt12 - XSt22)*(-1 + XSt22)*
     &          (2 + log(XGl2)*(3 + 2*log(XGl2))) + 
     &         XSt12*(-1 + XSt22)**2*
     &          (-2 - 3*log(XGl2) + XSt12*(3 + 4*log(XGl2)))*
     &          log(XSt12) + 
     &         XGl2*(-1 + XSt12)**2*(1 + XSt22)*
     &          T134p(MUE2,MSt22,MGl2)))/
     &    ((-1 + XSt12)**2*(XSt12 - XSt22)*(-1 + XSt22)**2)

        R82 = 1/4.D0*(-((-1 + XSt12)**2*yt1**2*
     &         (4*(XGl2**3 + 
     &              (XSt22 - XT2)**2*((1 - XSt22)*XSt22 + XT2) - 
     &              XGl2**2*(XSt22*(1 + XSt22) + XT2) - 
     &              XGl2*(XT2**2 + XSt22**2*(1 - 2*XSt22 + 4*XT2)))
     &             + (-(XGl2**3*(-2 + XSt22)) - 
     &              (-2 + XSt22)*(XSt22 - XT2)**2*XT2 + 
     &              XGl2**2*(-2 + XSt22)*(2*XSt22 + XT2) + 
     &              XGl2*(-2*XT2**2 + 
     &                 XSt22*((2 - XSt22)*XSt22 + (-4 + XT2)*XT2)))
     &             *(log(XGl2) - log(XSt22)))*
     &         (log(XGl2) - log(XSt22))) + 
     &      2*XT2*((XSt12 - XSt22)*
     &          (2*(-1 + XSt12)*(-1 + XSt22)*
     &             (XGl2**4 + (XSt12 - XT2)**2*(XSt22 - XT2)**2 - 
     &               2*XGl2**3*(XSt12 + XSt22 - XT2) + 
     &               XGl2**2*
     &                (XSt22*(XSt22 - 10*XT2) + 
     &                  XSt12*(XSt12 + 4*XSt22 - 10*XT2) + 
     &                  (12 - 6*XT2)*XT2) + 
     &               2*XGl2*
     &                (-(XSt12*XSt22*(XSt12 + XSt22)) + 
     &                  (XSt22*(-3 + 2*XSt22) + 
     &                     XSt12*(-3 + 2*XSt12 + 7*XSt22))*XT2 + 
     &                  (6 - 5*XSt12 - 5*XSt22)*XT2**2 + XT2**3))+
     &              (-2 + XSt12 + XSt22)*yt1**2*yt2**2*log(XGl2))-
     &           (-1 + XSt22)**2*
     &          (-(XGl2**2*XSt12) + 
     &            (-2 + XSt12)*(XSt12 - XT2)**2 + 
     &            2*XGl2*(XSt12 + XT2))*yt2**2*
     &          (log(XGl2) - log(XSt12)) + 
     &         (-1 + XSt12)**2*
     &          (-(XGl2**2*XSt22) + 
     &            (-2 + XSt22)*(XSt22 - XT2)**2 + 
     &            2*XGl2*(XSt22 + XT2))*yt1**2*
     &          (log(XGl2) - log(XSt22)))*(log(XGl2) - log(XT2)) - 
     &      (-1 + XSt12)*(XSt12 - XSt22)*(-1 + XSt22)*XT2*
     &       (-XGl2**4 + 2*XGl2**3*(XSt12 + XSt22) - 
     &         (XSt12 - XT2)**2*(XSt22 - XT2)**2 - 
     &         XGl2**2*(XSt22*(XSt22 - 6*XT2) + 
     &            XSt12*(XSt12 + 4*XSt22 - 6*XT2) + (8 - 2*XT2)*XT2
     &            ) + 2*XGl2*
     &          (XSt12*XSt22*(XSt12 + XSt22) - 
     &            ((-2 + XSt22)*XSt22 + 
     &               XSt12*(-2 + XSt12 + 6*XSt22))*XT2 + 
     &            (-4 + 3*XSt12 + 3*XSt22)*XT2**2))*
     &       (log(XGl2) - log(XT2))**2 + 
     &      (-1 + XSt22)**2*yt2**2*
     &       ((4*(XGl2**3 + 
     &               (XSt12 - XT2)**2*((1 - XSt12)*XSt12 + XT2) - 
     &               XGl2**2*(XSt12*(1 + XSt12) + XT2) - 
     &               XGl2*(XT2**2 + XSt12**2*(1 - 2*XSt12 + 4*XT2))
     &               ) + (-(XGl2**3*(-2 + XSt12)) - 
     &               (-2 + XSt12)*(XSt12 - XT2)**2*XT2 + 
     &               XGl2**2*(-2 + XSt12)*(2*XSt12 + XT2) + 
     &               XGl2*(-2*XT2**2 + 
     &                  XSt12*((2 - XSt12)*XSt12 + (-4 + XT2)*XT2))
     &               )*(log(XGl2) - log(XSt12)))*
     &          (log(XGl2) - log(XSt12)) - 
     &         2*XGl2*(-((XGl2 - XSt12)**2*(-2 + XGl2 + XSt12)) + 
     &            ((-4 + XSt12)*XSt12 + XGl2*(XGl2 + 2*XSt12))*
     &             XT2 + (2 + XGl2 + XSt12)*XT2**2 - XT2**3)*
     &          T134s(MSt12,MT2,MGl2,MGl2)) + 
     &      yt1**2*((XSt12 - XSt22)*(-2 + XSt12 + XSt22)*
     &          (XGl2 + XT2)*yt2**2*log(XGl2)*(4 + log(XGl2)) + 
     &         2*XGl2*(-1 + XSt12)**2*
     &          (-((XGl2 - XSt22)**2*(-2 + XGl2 + XSt22)) + 
     &            ((-4 + XSt22)*XSt22 + XGl2*(XGl2 + 2*XSt22))*
     &             XT2 + (2 + XGl2 + XSt22)*XT2**2 - XT2**3)*
     &          T134s(MSt22,MT2,MGl2,MGl2)) - 
     &      2*(XSt12 - XSt22)*
     &       (-((-1 + XSt12)*(-1 + XSt22)*
     &            ((3 + 2*XGl2)*(XGl2 - XSt12)**2*
     &               (XGl2 - XSt22)**2 + 
     &              2*(XGl2*
     &                  (-3*XSt12**2 + 
     &                    XGl2*
     &                     (XGl2*(15 + 4*XGl2 - 14*XSt12) + 
     &                       XSt12*(-4 + 6*XSt12))) + 
     &                 (-3*XSt12**2 + 
     &                    XGl2*
     &                     (XSt12*(5 + 3*XSt12) + 
     &                       XGl2*(-4 - 14*XGl2 + 5*XSt12)))*XSt22+
     &                   ((-3 + XSt12)*XSt12 + 
     &                    XGl2*(-3 + 6*XGl2 + 3*XSt12))*XSt22**2)*
     &               XT2 + 
     &              (-10*XGl2**3 + 
     &                 XGl2**2*(46 - 48*XSt12 - 48*XSt22) + 
     &                 XSt12**2*(3 - 4*XSt22) - 
     &                 4*XSt12*(-3 + XSt22)*XSt22 + 3*XSt22**2 + 
     &                 2*XGl2*
     &                  (XSt12*(-4 + 6*XSt12 + 5*XSt22) + 
     &                    XSt22*(-4 + 6*XSt22)))*XT2**2 + 
     &              2*(XGl2*(15 - 5*XGl2 - 14*XSt12 - 14*XSt22) + 
     &                 (-3 + XSt22)*XSt22 + 
     &                 XSt12*(-3 + XSt12 + 4*XSt22))*XT2**3 + 
     &              (3 + 8*XGl2 - 4*XSt12 - 4*XSt22)*XT2**4 + 
     &              2*XT2**5)) - 
     &         XGl2*(-2 + XSt12 + XSt22)*(-1 + XGl2 + XT2)*yt1**2*
     &          yt2**2*T134s(MUE2,MT2,MGl2,MGl2)))/
     &    ((-1 + XSt12)**2*(XSt12 - XSt22)*(-1 + XSt22)**2*yt1**2*
     &      yt2**2)

        R92 = 1/8.D0*(XSt12*(-1 + XSt22)**2*log(XSt12)*
     &       (4*(-2 + 3*XSt12)*(XSt12 - XSt22) + 
     &         (XSt22 + XSt12*(3 - 6*XSt12 + 2*XSt22))*log(XSt12))+
     &        (-1 + XSt12)**2*XSt22*
     &       ((3 - 6*XSt22)*XSt22 + XSt12*(1 + 2*XSt22))*
     &       log(XSt22)**2 + 
     &      4*((-1 + XSt22)*
     &          ((-1 + XSt12)*(XSt12 - XSt22)**2 + 
     &            (XSt12 - XSt22)*log(XGl2)*
     &             (XSt22 + 
     &               XSt12*
     &                (-1 + XSt22*(-1 - log(XSt12)) + log(XSt12) + 
     &                  XSt12*(1 + (-2 + 2*XSt22)*log(XSt12))))) - 
     &         (-1 + XSt12)*XSt22*
     &          ((-1 + XSt12)*(XSt12 - XSt22)*
     &             (-2 - log(XGl2) + XSt22*(3 + 2*log(XGl2))) - 
     &            XSt12*(-1 + XSt22)*(-2 + XSt12 + XSt22)*
     &             log(XSt12))*log(XSt22)))/
     &    ((-1 + XSt12)**2*(XSt12 - XSt22)**2*(-1 + XSt22)**2)

        Rat2 = -(1/16.D0*((-(XSt22*log(XSt22)) + 
     &          XSt12*((1 - XSt22)*log(XSt12) + XSt22*log(XSt22)))*
     &        (XSt22*((XGl2 - XSt12 + XT2)*
     &              (-(yt1*
     &                   (-(XSt12*(XSt12 - XSt22)*
     &                       (log(-XGl2 + XSt12 - XT2 - yt1) - 
     &                       log(-XGl2 + XSt12 - XT2 + yt1))) + 
     &                     4*XT2**2*
     &                      (-log(-XGl2 - XSt12 + XT2 - yt1) + 
     &                       log(-XGl2 - XSt12 + XT2 + yt1)) + 
     &                     (XSt12*(XSt12 - XSt22) + 4*XT2**2)*
     &                      (log(-XGl2 + XSt12 + XT2 - yt1) - 
     &                       log(-XGl2 + XSt12 + XT2 + yt1)))) + 
     &                (XGl2 - XSt12 - XT2)*
     &                 (4*XT2**2*
     &                    (log(-XGl2 - XSt12 + XT2 - yt1) + 
     &                      log(-XGl2 - XSt12 + XT2 + yt1)) + 
     &                   (XSt12*(-XSt12 + XSt22) - 4*XT2**2)*
     &                    (log(-XGl2 + XSt12 + XT2 - yt1) + 
     &                      log(-XGl2 + XSt12 + XT2 + yt1)))) + 
     &             XSt12*(XSt12 - XSt22)*
     &              ((XGl2 - XSt12 - XT2)*(XGl2 - XSt12 + XT2)*
     &                 (log(-XGl2 + XSt12 - XT2 - yt1) + 
     &                   log(-XGl2 + XSt12 - XT2 + yt1)) + 
     &                (XGl2 - XSt22 - XT2)*(XGl2 - XSt22 + XT2)*
     &                 (log(-XGl2 + XSt22 - XT2 - yt2) + 
     &                   log(-XGl2 + XSt22 - XT2 + yt2)))) + 
     &          XSt12*(2*XSt22*XT2*
     &              ((XSt12 - XSt22)*
     &                 ((XSt12 + XSt22)*(1 - log(XGl2)) - 
     &                   2*(XGl2 + XT2*(1 + log(XGl2)))) + 
     &                XSt12*(XSt12 - XSt22 + 8*XT2)*log(XSt12) + 
     &                XSt22*(XSt12 - XSt22 - 8*XT2)*log(XSt22) + 
     &                6*(-XSt12 + XSt22)*XT2*log(XT2)) + 
     &             (XGl2 - XSt22 + XT2)*
     &              (-(yt2*
     &                   (-((XSt12 - XSt22)*XSt22*
     &                       (log(-XGl2 + XSt22 - XT2 - yt2) - 
     &                       log(-XGl2 + XSt22 - XT2 + yt2))) + 
     &                     4*XT2**2*
     &                      (log(-XGl2 - XSt22 + XT2 - yt2) - 
     &                       log(-XGl2 - XSt22 + XT2 + yt2)) + 
     &                     (XSt12*XSt22 - XSt22**2 - 4*XT2**2)*
     &                      (log(-XGl2 + XSt22 + XT2 - yt2) - 
     &                       log(-XGl2 + XSt22 + XT2 + yt2)))) + 
     &                (XGl2 - XSt22 - XT2)*
     &                 (-4*XT2**2*
     &                    (log(-XGl2 - XSt22 + XT2 - yt2) + 
     &                      log(-XGl2 - XSt22 + XT2 + yt2)) + 
     &                   (-(XSt12*XSt22) + XSt22**2 + 4*XT2**2)*
     &                    (log(-XGl2 + XSt22 + XT2 - yt2) + 
     &                      log(-XGl2 + XSt22 + XT2 + yt2)))))))/
     &      ((-1 + XSt12)*XSt12*(XSt12 - XSt22)**2*(-1 + XSt22)*
     &        XSt22*XT2**2))

        Rlt2 = -(1/8.D0*(2*XT2*
     &         (2*(-XSt12 + XSt22)*XT2*log(XGl2)**2 + 
     &           XSt12*(-1 + XSt22)*(XSt12 - 2*XT2)*
     &            (log(XGl2) - log(XSt12))**2 - 
     &           (-1 + XSt12)*XSt22*(XSt22 - 2*XT2)*
     &            (log(XGl2) - log(XSt22))**2) - 
     &        XSt22*(log(XGl2) - log(XSt22))*
     &         (2*(-XSt12 + XSt22)*XT2*log(XGl2) + 
     &           (-1 + XSt12)*
     &            (-2*(-2*XGl2 + XSt12 + XSt22)*XT2 + 
     &              (XGl2 - XSt12 + XT2)*(-XGl2 + XSt12 + XT2)*
     &               (log(-XGl2 + XSt12 - XT2 - yt1) - 
     &                 log(-XGl2 + XSt12 + XT2 - yt1) + 
     &                 log(-XGl2 + XSt12 - XT2 + yt1) - 
     &                 log(-XGl2 + XSt12 + XT2 + yt1)) + 
     &              (XGl2 - XSt22 + XT2)*(-XGl2 + XSt22 + XT2)*
     &               (log(-XGl2 + XSt22 - XT2 - yt2) - 
     &                 log(-XGl2 + XSt22 + XT2 - yt2) + 
     &                 log(-XGl2 + XSt22 - XT2 + yt2) - 
     &                 log(-XGl2 + XSt22 + XT2 + yt2)))) + 
     &        XSt12*(log(XGl2) - log(XSt12))*
     &         (-2*(XSt12 - XSt22)*XT2*
     &            ((-1 + XSt22)*log(XGl2) - XSt22*log(XSt22)) + 
     &           (-1 + XSt22)*
     &            (-2*(-2*XGl2 + XSt12 + XSt22)*XT2 + 
     &              (XGl2 - XSt12 + XT2)*(-XGl2 + XSt12 + XT2)*
     &               (log(-XGl2 + XSt12 - XT2 - yt1) - 
     &                 log(-XGl2 + XSt12 + XT2 - yt1) + 
     &                 log(-XGl2 + XSt12 - XT2 + yt1) - 
     &                 log(-XGl2 + XSt12 + XT2 + yt1)) + 
     &              (XGl2 - XSt22 + XT2)*(-XGl2 + XSt22 + XT2)*
     &               (log(-XGl2 + XSt22 - XT2 - yt2) - 
     &                 log(-XGl2 + XSt22 + XT2 - yt2) + 
     &                 log(-XGl2 + XSt22 - XT2 + yt2) - 
     &                 log(-XGl2 + XSt22 + XT2 + yt2)))) + 
     &        (XSt12 - XSt22)*log(XGl2)*
     &         (-2*(-2*XGl2 + XSt12 + XSt22)*XT2 + 
     &           (XGl2 - XSt12 + XT2)*(-XGl2 + XSt12 + XT2)*
     &            (log(-XGl2 + XSt12 - XT2 - yt1) - 
     &              log(-XGl2 + XSt12 + XT2 - yt1) + 
     &              log(-XGl2 + XSt12 - XT2 + yt1) - 
     &              log(-XGl2 + XSt12 + XT2 + yt1)) + 
     &           (XGl2 - XSt22 + XT2)*(-XGl2 + XSt22 + XT2)*
     &            (log(-XGl2 + XSt22 - XT2 - yt2) - 
     &              log(-XGl2 + XSt22 + XT2 - yt2) + 
     &              log(-XGl2 + XSt22 - XT2 + yt2) - 
     &              log(-XGl2 + XSt22 + XT2 + yt2))) + 
     &        (-(XSt22*log(XSt22)) + 
     &           XSt12*((1 - XSt22)*log(XSt12) + XSt22*log(XSt22)))
     &          *((XGl2 - XSt12 + XT2)*yt1*
     &            (-log(-XGl2 + XSt12 - XT2 - yt1) + 
     &              log(-XGl2 + XSt12 + XT2 - yt1) + 
     &              log(-XGl2 + XSt12 - XT2 + yt1) - 
     &              log(-XGl2 + XSt12 + XT2 + yt1)) + 
     &           (XGl2 - XSt22 + XT2)*yt2*
     &            (-log(-XGl2 + XSt22 - XT2 - yt2) + 
     &              log(-XGl2 + XSt22 + XT2 - yt2) + 
     &              log(-XGl2 + XSt22 - XT2 + yt2) - 
     &              log(-XGl2 + XSt22 + XT2 + yt2))))/
     &      ((-1 + XSt12)*(XSt12 - XSt22)*(-1 + XSt22)*XT2**2))

        RMSt12 = 1/4.D0*(2*XSt12*
     &       ((-1 + XSt22)*(-XSt12**2 + XSt22)*
     &          (XGl2 + XSt12 + XT2)*log(XSt12)**2 + 
     &         (2*XGl2 + XSt12 + 3*XT2)*log(XGl2)*
     &          ((XSt12**2 - XSt22)*(-1 + XSt22)*log(XSt12) - 
     &            (-1 + XSt12)*
     &             (XSt22*(1 - XSt22 - log(XSt22)) + 
     &               XSt12*(-1 + XSt22*(1 + log(XSt22)))))) + 
     &      (-1 + XSt12)*((XSt12 - XSt22)*(-1 + XSt22)*
     &          (XSt12*(-2*(XSt12 + 5*(XGl2 + XT2)) + 
     &               2*XT2*log(XT2)) + 
     &            (XGl2 - XSt12 + XT2)*
     &             ((XGl2 - XSt12 - XT2 + yt1)*
     &                (log(-XGl2 - XSt12 + XT2 - yt1) - 
     &                  log(-XGl2 + XSt12 + XT2 - yt1)) + 
     &               (XGl2 - XSt12 - XT2 - yt1)*
     &                log(-XGl2 - XSt12 + XT2 + yt1) + 
     &               (-XGl2 + XSt12 + XT2 + yt1)*
     &                log(-XGl2 + XSt12 + XT2 + yt1))) + 
     &         (-1 + XSt12)*XSt22*log(XSt22)*
     &          (XSt12*(-2*(XSt12 + 5*(XGl2 + XT2)) + 
     &               2*((XGl2 + XT2)*log(XSt22) + XT2*log(XT2))) + 
     &            (XGl2 - XSt12 + XT2)*
     &             ((XGl2 - XSt12 - XT2 + yt1)*
     &                (log(-XGl2 - XSt12 + XT2 - yt1) - 
     &                  log(-XGl2 + XSt12 + XT2 - yt1)) + 
     &               (XGl2 - XSt12 - XT2 - yt1)*
     &                log(-XGl2 - XSt12 + XT2 + yt1) + 
     &               (-XGl2 + XSt12 + XT2 + yt1)*
     &                log(-XGl2 + XSt12 + XT2 + yt1)))) + 
     &      log(XSt12)*(2*(-1 + XSt12)**2*XSt12**2*XSt22*
     &          log(XSt22) - 
     &         (-1 + XSt22)*
     &          (2*XSt12*(-2*XSt12**3 + 
     &               XGl2*(3*XSt22 + 
     &                  XSt12*(-7*XSt12 + 2*(1 + XSt22))) + 
     &               XSt12**2*(1 + XSt22 - 7*XT2) + 
     &               (3*XSt22 + 2*XSt12*(1 + XSt22))*XT2) + 
     &            (XSt12**2 - XSt22)*
     &             (2*XSt12*XT2*log(XT2) + 
     &               (XGl2 - XSt12 + XT2)*
     &                ((XGl2 - XSt12 - XT2 + yt1)*
     &                   (log(-XGl2 - XSt12 + XT2 - yt1) - 
     &                     log(-XGl2 + XSt12 + XT2 - yt1)) + 
     &                  (XGl2 - XSt12 - XT2 - yt1)*
     &                   log(-XGl2 - XSt12 + XT2 + yt1) + 
     &                  (-XGl2 + XSt12 + XT2 + yt1)*
     &                   log(-XGl2 + XSt12 + XT2 + yt1))))))/
     &    ((-1 + XSt12)**2*XSt12*(XSt12 - XSt22)**2*(-1 + XSt22))

        RMSt22 = 1/4.D0*(2*XSt22*
     &       ((-1 + XSt12)*(XSt12 - XSt22)*(-1 + XSt22)*
     &          (XSt22 + 5*(XGl2 + XT2)) + 
     &         XSt12*(-1 + XSt22)**2*(XGl2 + XT2)*log(XSt12)**2 + 
     &         (2*XGl2 + XSt22 + 3*XT2)*log(XGl2)*
     &          (-(XSt12*(-1 + XSt22)**2*log(XSt12)) + 
     &            (-1 + XSt12)*
     &             (XSt22*(1 + XSt22*(-1 + log(XSt22))) + 
     &               XSt12*(-1 + XSt22 - log(XSt22))))) + 
     &      XSt12*(-1 + XSt22)**2*log(XSt12)*
     &       (2*XSt22**2*log(XSt22) + 
     &         XSt22*(-2*(XSt22 + 5*(XGl2 + XT2)) + 
     &            2*XT2*log(XT2)) + 
     &         (XGl2 - XSt22 + XT2)*
     &          ((XGl2 - XSt22 - XT2 + yt2)*
     &             (log(-XGl2 - XSt22 + XT2 - yt2) - 
     &               log(-XGl2 + XSt22 + XT2 - yt2)) + 
     &            (XGl2 - XSt22 - XT2 - yt2)*
     &             log(-XGl2 - XSt22 + XT2 + yt2) + 
     &            (-XGl2 + XSt22 + XT2 + yt2)*
     &             log(-XGl2 + XSt22 + XT2 + yt2))) - 
     &      (-1 + XSt12)*(2*XSt22*(-XSt12 + XSt22**2)*
     &          (XGl2 + XSt22 + XT2)*log(XSt22)**2 + 
     &         (XSt12 - XSt22)*(-1 + XSt22)*
     &          (2*XSt22*XT2*log(XT2) + 
     &            (XGl2 - XSt22 + XT2)*
     &             ((XGl2 - XSt22 - XT2 + yt2)*
     &                (log(-XGl2 - XSt22 + XT2 - yt2) - 
     &                  log(-XGl2 + XSt22 + XT2 - yt2)) + 
     &               (XGl2 - XSt22 - XT2 - yt2)*
     &                log(-XGl2 - XSt22 + XT2 + yt2) + 
     &               (-XGl2 + XSt22 + XT2 + yt2)*
     &                log(-XGl2 + XSt22 + XT2 + yt2))) + 
     &         log(XSt22)*(2*XSt22*
     &             (-2*XSt22**3 + 
     &               XGl2*((2 - 7*XSt22)*XSt22 + 
     &                  XSt12*(3 + 2*XSt22)) + 
     &               XSt22**2*(1 + XSt12 - 7*XT2) + 
     &               (3*XSt12 + 2*(1 + XSt12)*XSt22)*XT2) - 
     &            (XSt12 - XSt22**2)*
     &             (2*XSt22*XT2*log(XT2) + 
     &               (XGl2 - XSt22 + XT2)*
     &                ((XGl2 - XSt22 - XT2 + yt2)*
     &                   (log(-XGl2 - XSt22 + XT2 - yt2) - 
     &                     log(-XGl2 + XSt22 + XT2 - yt2)) + 
     &                  (XGl2 - XSt22 - XT2 - yt2)*
     &                   log(-XGl2 - XSt22 + XT2 + yt2) + 
     &                  (-XGl2 + XSt22 + XT2 + yt2)*
     &                   log(-XGl2 + XSt22 + XT2 + yt2))))))/
     &    ((-1 + XSt12)*(XSt12 - XSt22)**2*(-1 + XSt22)**2*XSt22)

	Rew = C_F*Re(R22 + R42 + R52 + R72 + R82 + R92 +
     &    Rat2 + Rlt2 + RMSt12 + RMSt22)/TLdb(MSt12, MSt22, MUE2)
	anomalous = -C_F
	finscale = 2*log(scale**2/MGl2)
c??? should be MUE2???

	Deltab2L_ew = Rew*fnorm + anomalous + finscale
	end


************************************************************************

	RealType function Deltab2L_qcd(scale)
	implicit none
	RealType scale

#include "FH.h"

	integer Nu, Nd, Nf
	parameter (Nu = 2, Nd = 2, Nf = Nu + Nd + 1)

	ComplexType Ieps
	parameter (Ieps = (1D0, 1D-15))

	RealType MSu12, MSu22, MSd12, MSd22
	RealType MSt12, MSt22, MSb12, MSb22
	equivalence (MSu12, MSf2(1,3,2)), (MSu22, MSf2(2,3,2))
	equivalence (MSt12, MSf2(1,3,3)), (MSt22, MSf2(2,3,3))
	equivalence (MSd12, MSdL2(1,2)), (MSd22, MSdL2(2,2))
	equivalence (MSb12, MSdL2(1,3)), (MSb22, MSdL2(2,3))

	RealType fnorm, bo, anomalous, finscale
	ComplexType d, r, XGl2, XT2
	ComplexType XSt12, XSt22, XSb12, XSb22
	ComplexType XSu12, XSu22, XSd12, XSd22
	ComplexType XpSt12, XmSt12, XpSt22, XmSt22
	ComplexType R12, R22, R32, R42, R52, R62, R72, R82, R92
	ComplexType RCtCA, RAlsCA2, RMGlCA2, RCA
	ComplexType RCtCF, RAlsCF2, RMSb12, RMSb22, RCF
	ComplexType Ru102, RuCtTF, RuAlsTF2, RuMGlTF2
	ComplexType Rd102, RdCtTF, RdAlsTF2, RdMGlTF2
	ComplexType Rb102, RbCtTF, RbAlsTF2, RbMGlTF2
	ComplexType Rt102, RtCtTF, RtAlsTF2, RtMGlTF2
	ComplexType RTF, RTFb, RTFt
	RealType Rqcd

	RealType T134s, T134p, TLdb
	external T134s, T134p, TLdb

#define Q2 MT2
#define XQ2 XT2

	fnorm = 4/MGl2
	d = Ieps/MGl2
	XGl2  = Ieps
	XT2   = MT2*d
	XSt12 = MSt12*d
	XSt22 = MSt22*d
	XSb12 = MSb12*d
	XSb22 = MSb22*d
	XSu12 = MSu12*d
	XSu22 = MSu22*d
	XSd12 = MSd12*d
	XSd22 = MSd22*d

	d = (MGl2 + MT2 - MSt12)/MGl2*(.5D0/Ieps)
	r = sqrt(d**2 - MT2/MGl2/Ieps)
	XpSt12 = d + r
	XmSt12 = d - r
	d = (MGl2 + MT2 - MSt22)/MGl2*(.5D0/Ieps)
	r = sqrt(d**2 - MT2/MGl2/Ieps)
	XpSt22 = d + r
	XmSt22 = d - r


        R12 = 1/16.D0*(-(XSb12*(-1 + XSb22)**2*log(XSb12)*
     &         (4 - log(XSb12) + XSb12*(-4 + 2*log(XSb12)))) + 
     &      (-1 + XSb12)**2*XSb22*log(XSb22)*
     &       (4 - log(XSb22) + XSb22*(-4 + 2*log(XSb22))) + 
     &      4*(-T134p(MSb12,MGl2,MGl2) + T134p(MSb22,MGl2,MGl2) + 
     &         XSb22*((-2 + XSb22)*T134p(MGl2,MGl2,MGl2) + 
     &            (2 - XSb22)*T134p(MSb12,MGl2,MGl2) + 
     &            T134p(MSb22,MGl2,MGl2) - T134p(MSb22,MSb22,MGl2))
     &           + XSb12*(2*T134p(MGl2,MGl2,MGl2) - 
     &            T134p(MSb12,MGl2,MGl2) + 
     &            T134p(MSb12,MSb12,MGl2) - 
     &            2*T134p(MSb22,MGl2,MGl2) + 
     &            XSb12*(-T134p(MGl2,MGl2,MGl2) + 
     &               T134p(MSb22,MGl2,MGl2) + 
     &               XSb22*
     &                (T134p(MSb22,MGl2,MGl2) - 
     &                  T134p(MSb22,MSb22,MGl2))) + 
     &            XSb22*(XSb22*
     &                (-T134p(MSb12,MGl2,MGl2) + 
     &                  T134p(MSb12,MSb12,MGl2)) - 
     &               2*(T134p(MSb12,MSb12,MGl2) + 
     &                  T134p(MSb22,MGl2,MGl2)) + 
     &               2*(T134p(MSb12,MGl2,MGl2) + 
     &                  T134p(MSb22,MSb22,MGl2))))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        R22 = -(1/2.D0*(2*(XSb12 - XSb22)*(-1 + XSb12*XSb22)*
     &         (7 + log(XGl2)*(6 + 2*log(XGl2))) + 
     &        (-1 + XSb22)**2*
     &         (XSb12*(2*XSb12 + 4*(1 + log(XGl2)) - log(XSb12))*
     &            log(XSb12) + 
     &           XGl2*(1 + XSb12)*T134p(MGl2,MSb12,MGl2)) + 
     &        (-1 + XSb12)**2*
     &         (XSb22*log(XSb22)*
     &            (-2*(XSb22 + 2*(1 + log(XGl2))) + log(XSb22)) - 
     &           XGl2*(1 + XSb22)*T134p(MGl2,MSb22,MGl2)))/
     &      ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2))

        R32 = 1/8.D0*(14*(1 + XSb12**2)*XSb22 + 
     &      3*(-1 + XSb22)*(XSb12**2 + XSb22)*
     &       T134p(MGl2,MGl2,MGl2) + 
     &      (-2 + XSb22*(4 - 8*XSb12 + (-2 + 4*XSb12)*XSb22))*
     &       T134p(MSb12,MGl2,MGl2) + 
     &      XSb12*(-14*(1 + XSb22**2) + 
     &         (-1 + XSb22)**2*(2 + 4*XSb12 - log(XSb12))*
     &          log(XSb12) - 
     &         3*(-1 + XSb22**2)*T134p(MGl2,MGl2,MGl2) + 
     &         4*T134p(MSb12,MGl2,MGl2)) + 
     &      (-1 + XSb12)**2*
     &       (XSb22*(log(XSb22)*(-2 - 4*XSb22 + log(XSb22)) - 
     &            4*T134p(MSb22,MGl2,MGl2)) + 
     &         2*T134p(MSb22,MGl2,MGl2)))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        R42 = 1/2.D0*(4*(XSb12 - XSb22)*(-1 + XSb12*XSb22)*
     &       (7 + log(XGl2)*(6 + 2*log(XGl2))) + 
     &      (-1 + XSb22)**2*
     &       (XSb12*(7 + 5*XSb12 + 8*log(XGl2) - 2*log(XSb12))*
     &          log(XSb12) + 
     &         2*XGl2*(1 + XSb12)*T134p(MSb12,MGl2,MGl2)) + 
     &      (-1 + XSb12)**2*
     &       (-(XSb22*(7 + 5*XSb22 + 8*log(XGl2) - 2*log(XSb22))*
     &            log(XSb22)) - 
     &         2*XGl2*(1 + XSb22)*T134p(MSb22,MGl2,MGl2)))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        R52 = 1/16.D0*(XSb12*(-1 + XSb22)**2*log(XSb12)*
     &       (4 - log(XSb12) + XSb12*(-4 + 2*log(XSb12))) + 
     &      (-1 + XSb12)*(-((-1 + XSb12)*XSb22*log(XSb22)*
     &            (4 - log(XSb22) + XSb22*(-4 + 2*log(XSb22)))) + 
     &         2*(-1 + XSb22)*
     &          ((1 - XSb12)*T134p(MSb12,MGl2,MGl2) + 
     &            (XSb12 - XSb22)*T134p(MSb12,MSb22,MGl2) + 
     &            (-1 + XSb22)*T134p(MSb22,MGl2,MGl2))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        R62 = 1/8.D0*(XSb12*(-1 + XSb22)**2*(-2 + log(XSb12))*
     &       log(XSb12) - (-1 + XSb12)**2*XSb22*(-2 + log(XSb22))*
     &       log(XSb22) + 2*
     &       (T134p(MSb12,MGl2,MGl2) + 
     &         XSb22*((3 - XSb22)*T134p(MGl2,MGl2,MGl2) + 
     &            (-2 + XSb22)*T134p(MSb12,MGl2,MGl2) - 
     &            T134p(MSb22,MGl2,MGl2)) - 
     &         T134p(MSb22,MGl2,MGl2) + 
     &         XSb12*((-3 - XSb22**2)*T134p(MGl2,MGl2,MGl2) + 
     &            T134p(MSb12,MGl2,MGl2) + 
     &            2*T134p(MSb22,MGl2,MGl2) + 
     &            XSb22*((-2 + XSb22)*T134p(MSb12,MGl2,MGl2) + 
     &               2*T134p(MSb22,MGl2,MGl2)) + 
     &            XSb12*((1 + XSb22)*T134p(MGl2,MGl2,MGl2) + 
     &               (-1 - XSb22)*T134p(MSb22,MGl2,MGl2)))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        R72 = 1/8.D0*(4*((-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22)*
     &          (2 + log(XGl2)*(3 + 2*log(XGl2))) + 
     &         XSb12*(-1 + XSb22)**2*
     &          (-2 - 3*log(XGl2) + XSb12*(3 + 4*log(XGl2)))*
     &          log(XSb12) + 
     &         XGl2*(-1 + XSb12)**2*(1 + XSb22)*
     &          T134p(MGl2,MSb22,MGl2)) + 
     &      (-1 + XSb22)**2*
     &       ((5 - 8*XSb12)*XSb12*log(XSb12)**2 - 
     &         4*XGl2*(1 + XSb12)*
     &          (T134p(MGl2,MSb12,MGl2) - T134p(MSb12,MSb12,MGl2)))
     &        + (-1 + XSb12)**2*
     &       (XSb22*(4*(2 + XSb22*(-3 - 4*log(XGl2)) + 
     &               3*log(XGl2))*log(XSb22) + 
     &            (-5 + 8*XSb22)*log(XSb22)**2) - 
     &         4*XGl2*(1 + XSb22)*T134p(MSb22,MSb22,MGl2)))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        R82 = 1/8.D0*((-1 + XSb22)**2*log(XSb12)*
     &       (4 - 2*log(XSb12) + XSb12*(4 - 4*XSb12 + log(XSb12)))+
     &        (-1 + XSb12)*
     &       ((-1 + XSb12)*log(XSb22)*
     &          (4*(-1 + (-1 + XSb22)*XSb22) - 
     &            (-2 + XSb22)*log(XSb22)) + 
     &         2*(-1 + XSb22)*
     &          (XSb22*(5 - T134p(MSb12,MGl2,MGl2)) + 
     &            T134p(MSb12,MGl2,MGl2) + 
     &            XSb12*(-5 + T134p(MSb22,MGl2,MGl2)) - 
     &            T134p(MSb22,MGl2,MGl2))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        R92 = 1/8.D0*(XSb12*(-1 + XSb22)**2*log(XSb12)*
     &       (4*(-2 + 3*XSb12)*(XSb12 - XSb22) + 
     &         (XSb22 + XSb12*(3 - 6*XSb12 + 2*XSb22))*log(XSb12))+
     &        (-1 + XSb12)**2*XSb22*
     &       ((3 - 6*XSb22)*XSb22 + XSb12*(1 + 2*XSb22))*
     &       log(XSb22)**2 + 
     &      4*((-1 + XSb22)*
     &          ((-1 + XSb12)*(XSb12 - XSb22)**2 + 
     &            (XSb12 - XSb22)*log(XGl2)*
     &             (XSb22 + 
     &               XSb12*
     &                (-1 + XSb22*(-1 - log(XSb12)) + log(XSb12) + 
     &                  XSb12*(1 + (-2 + 2*XSb22)*log(XSb12))))) - 
     &         (-1 + XSb12)*XSb22*
     &          ((-1 + XSb12)*(XSb12 - XSb22)*
     &             (-2 - log(XGl2) + XSb22*(3 + 2*log(XGl2))) - 
     &            XSb12*(-1 + XSb22)*(-2 + XSb12 + XSb22)*
     &             log(XSb12))*log(XSb22)))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)**2*(-1 + XSb22)**2)

        RAlsCA2 = 1/48.D0*(-(XSb12*(-1 + XSb22)*log(XSb12)*
     &         (-10 + 9*log(XSb12))) + 
     &      (-1 + XSb12)*XSb22*log(XSb22)*(-10 + 9*log(XSb22)))/
     &    ((-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22))

        RMGlCA2 = 1/16.D0*(28*(-1 + XSb12)*(XSb12 - XSb22)*
     &       (-1 + XSb22) - 
     &      XSb12*(1 + XSb12)*(-1 + XSb22)**2*log(XSb12)*
     &       (-14 + 3*log(XSb12)) + 
     &      (-1 + XSb12)**2*XSb22*(1 + XSb22)*log(XSb22)*
     &       (-14 + 3*log(XSb22)))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        RAlsCF2 = 1/8.D0*(XSb22*log(XSb22) + 
     &      XSb12*((-1 + XSb22)*log(XSb12) - XSb22*log(XSb22)))/
     &    ((-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22))

        RMSb12 = 1/4.D0*(-(log(XSb12)*
     &         (-(XSb12*(-1 + XSb22)*
     &              (-2*XSb12*(-1 + XSb12*(3 + XSb12)) + 
     &                (3 + XSb12*(2 + XSb12))*XSb22)) + 
     &           (-1 + XSb12)**2*
     &            ((XSb12**2 - XSb22)*(-1 + XSb22)*
     &               log(1 - XSb12) + XSb12**2*XSb22*log(XSb22))))+
     &        (-1 + XSb12)*
     &       (XSb12*(5 + XSb12) + (-1 + XSb12)**2*log(1 - XSb12))*
     &       (XSb22*(1 - XSb22 - log(XSb22)) + 
     &         XSb12*(-1 + XSb22*(1 + log(XSb22)))) + 
     &      XSb12*(XSb22*((1 - XSb22)*log(XSb12)**2 - 
     &            log(XSb22)**2) + 
     &         XSb12*(XSb22*
     &             ((1 - XSb22)*log(XSb12)**2 + 2*log(XSb22)**2) + 
     &            XSb12*((-1 + XSb12*(-1 + XSb22) + XSb22)*
     &                log(XSb12)**2 - XSb22*log(XSb22)**2))))/
     &    ((-1 + XSb12)**2*XSb12*(XSb12 - XSb22)**2*(-1 + XSb22))

        RMSb22 = -(1/4.D0*(XSb12*(-1 + XSb22)**2*log(XSb12)*
     &         (-((-1 + XSb22)**2*log(1 - XSb22)) + 
     &           XSb22*(-5 + log(XSb12) + XSb22*(-1 + log(XSb22))))
     &          + (-1 + XSb12)*
     &         ((-1 + XSb22)**2*log(1 - XSb22)*
     &            (XSb22*(1 + XSb22*(-1 + log(XSb22))) + 
     &              XSb12*(-1 + XSb22 - log(XSb22))) + 
     &           XSb22*((XSb12 - XSb22)*(-1 + XSb22)*(5 + XSb22) + 
     &              log(XSb22)*
     &               (-(XSb12*(3 + XSb22*(2 + XSb22))) + 
     &                 2*XSb22*(-1 + XSb22*(3 + XSb22)) + 
     &                 (1 + XSb22)*(XSb12 - XSb22**2)*log(XSb22))))
     &        )/
     &      ((-1 + XSb12)*(XSb12 - XSb22)**2*(-1 + XSb22)**2*XSb22)
     &     )

        Ru102 = 1/8.D0*(XSb12*(-1 + XSb22)**2*log(XSb12)*
     &       (-((XSu12 + XSu22)*(-4 + log(XSb12))) - 
     &         2*(XSu12*log(XSu12) + XSu22*log(XSu22))) + 
     &      (-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22)*
     &       (12 + XSu12*(12 + (-8 + log(XSu12))*log(XSu12)) + 
     &         XSu22*(12 + (-8 + log(XSu22))*log(XSu22))) - 
     &      2*((XSb12 - XSb22)*
     &          ((XSb22*(-1 - XSu12) + 
     &               XSb12*(-1 + 2*XSb22 - XSu12) + 2*XSu12)*
     &             T134p(MGl2,MSu12,MGl2) + 
     &            (XSb22*(-1 - XSu22) + 
     &               XSb12*(-1 + 2*XSb22 - XSu22) + 2*XSu22)*
     &             T134p(MGl2,MSu22,MGl2)) + 
     &         (-1 + XSb22)**2*
     &          (-(XSu12*T134p(MSb12,MSu12,MGl2)) - 
     &            XSu22*T134p(MSb12,MSu22,MGl2) + 
     &            XSb12*(T134p(MSb12,MSu12,MGl2) + 
     &               T134p(MSb12,MSu22,MGl2)))) + 
     &      (-1 + XSb12)**2*
     &       (XSb22*log(XSb22)*
     &          ((XSu12 + XSu22)*(-4 + log(XSb22)) + 
     &            2*(XSu12*log(XSu12) + XSu22*log(XSu22))) + 
     &         2*(-(XSu12*T134p(MSb22,MSu12,MGl2)) - 
     &            XSu22*T134p(MSb22,MSu22,MGl2) + 
     &            XSb22*(T134p(MSb22,MSu12,MGl2) + 
     &               T134p(MSb22,MSu22,MGl2)))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        Rd102 = 1/8.D0*(XSb12*(-1 + XSb22)**2*log(XSb12)*
     &       (-((XSd12 + XSd22)*(-4 + log(XSb12))) - 
     &         2*(XSd12*log(XSd12) + XSd22*log(XSd22))) + 
     &      (-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22)*
     &       (12 + XSd12*(12 + (-8 + log(XSd12))*log(XSd12)) + 
     &         XSd22*(12 + (-8 + log(XSd22))*log(XSd22))) - 
     &      2*((XSb12 - XSb22)*
     &          ((XSb22*(-1 - XSd12) + 
     &               XSb12*(-1 + 2*XSb22 - XSd12) + 2*XSd12)*
     &             T134p(MGl2,MSd12,MGl2) + 
     &            (XSb22*(-1 - XSd22) + 
     &               XSb12*(-1 + 2*XSb22 - XSd22) + 2*XSd22)*
     &             T134p(MGl2,MSd22,MGl2)) + 
     &         (-1 + XSb22)**2*
     &          (-(XSd12*T134p(MSb12,MSd12,MGl2)) - 
     &            XSd22*T134p(MSb12,MSd22,MGl2) + 
     &            XSb12*(T134p(MSb12,MSd12,MGl2) + 
     &               T134p(MSb12,MSd22,MGl2)))) + 
     &      (-1 + XSb12)**2*
     &       (XSb22*log(XSb22)*
     &          ((XSd12 + XSd22)*(-4 + log(XSb22)) + 
     &            2*(XSd12*log(XSd12) + XSd22*log(XSd22))) + 
     &         2*(-(XSd12*T134p(MSb22,MSd12,MGl2)) - 
     &            XSd22*T134p(MSb22,MSd22,MGl2) + 
     &            XSb22*(T134p(MSb22,MSd12,MGl2) + 
     &               T134p(MSb22,MSd22,MGl2)))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        Rb102 = 1/8.D0*(XSb12*log(XSb12)*
     &       ((-1 + XSb22)*
     &          (4*((-3 + XSb22)*XSb22 + 
     &               XSb12*(1 - 2*XSb12 + 3*XSb22)) + 
     &            (XSb12*(2 + XSb12 - 4*XSb22) + 
     &               (2 - XSb22)*XSb22)*log(XSb12)) + 
     &         2*(XSb12 - XSb22)*XSb22*(-2 + XSb12 + XSb22)*
     &          log(XSb22)) + 
     &      (-1 + XSb12)*XSb22*log(XSb22)*
     &       (-4*((1 - 2*XSb22)*XSb22 + 
     &            XSb12*(-3 + XSb12 + 3*XSb22)) + 
     &         ((-2 - XSb22)*XSb22 + XSb12*(-2 + XSb12 + 4*XSb22))*
     &          log(XSb22)) + 
     &      2*(XSb12 - XSb22)*
     &       ((-1 + XSb12)*(XSb12 - XSb22)*
     &          T134p(MSb12,MGl2,MGl2) - 
     &         (2 + (-2 + XSb12)*XSb12 + (-2 + XSb22)*XSb22)*
     &          T134p(MSb12,MSb22,MGl2) + 
     &         (-1 + XSb22)*
     &          (-6 + XSb12*
     &             (6*(XSb12 + XSb22) - T134p(MSb22,MGl2,MGl2)) + 
     &            XSb22*(-6 + T134p(MSb22,MGl2,MGl2)))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        Rt102 = 1/8.D0*((-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22)*
     &       (-2*XQ2*(XQ2**4 - (-1 + XSt12)**2*(-1 + XSt22)**2 - 
     &            XQ2**3*(2 + XSt12 + XSt22) + 
     &            XQ2*(2 + (-1 + XSt22)*XSt22 + 
     &               XSt12**2*(1 + XSt22) + 
     &               XSt12*(-1 + (-4 + XSt22)*XSt22)))*log(XQ2)**2+
     &           XSt12*(1 + (-2 + XSt22)*XSt22 + 
     &            XQ2*(XQ2 - 2*(1 + XSt22)))*log(XSt12)*
     &          (4*(-2 + (4 - 2*XSt12)*XSt12 + 
     &               XQ2*(1 + XQ2 + XSt12)) + 
     &            (1 - XQ2**2 + (-2 + XSt12)*XSt12)*log(XSt12)) + 
     &         4*(-6*XQ2**5 + 
     &            3*(-1 + XSt12)**2*(-1 + XSt22)**2*
     &             (1 + XSt12 + XSt22) + 
     &            XQ2**4*(13 + 8*(XSt12 + XSt22)) - 
     &            XQ2**3*(6 + XSt22*(2 + 5*XSt22) + 
     &               XSt12*(2 + 5*XSt12 + 8*XSt22)) + 
     &            XQ2**2*((5 + 6*XSt12**2)*XSt22 + 
     &               3*(XSt12**3 + XSt22**3) + 
     &               XSt12*(5 + XSt22*(4 + 6*XSt22))) + 
     &            XQ2*(2*(-2 - 3*XSt12**3*(1 + XSt22) + 
     &                  XSt22*(-4 + (4 - 3*XSt22)*XSt22) + 
     &                  XSt12**2*(4 + (5 - 2*XSt22)*XSt22) + 
     &                  XSt12*
     &                   (-4 + XSt22*(6 + (5 - 3*XSt22)*XSt22))) + 
     &               log(XQ2)*
     &                (4*XQ2**4 - 
     &                  2*(-1 + XSt12)**2*(-1 + XSt22)**2 - 
     &                  5*XQ2**3*(2 + XSt12 + XSt22) + 
     &                  XQ2**2*
     &                   (6 + XSt22*(2 + XSt22) + 
     &                     XSt12*(2 + XSt12 + 4*XSt22)) + 
     &                  XQ2*
     &                   (2 + (-1 + XSt22)*XSt22 + 
     &                     XSt12**2*(1 + XSt22) + 
     &                     XSt12*(-1 + (-4 + XSt22)*XSt22)) + 
     &                  XSt12*(-1 - XQ2 + XSt12)*
     &                   (1 + (-2 + XSt22)*XSt22 + 
     &                     XQ2*(XQ2 - 2*(1 + XSt22)))*log(XSt12) + 
     &                  (1 + (-2 + XSt12)*XSt12 + 
     &                     XQ2*(XQ2 - 2*(1 + XSt12)))*XSt22*
     &                   (-1 - XQ2 + XSt22)*log(XSt22))))) + 
     &      (-1 + (2 - XSt22)*XSt22 + XQ2*(-XQ2 + 2*(1 + XSt22)))*
     &       (-2*(-1 + (2 - XSt12)*XSt12 + 
     &            XQ2*(-XQ2 + 2*(1 + XSt12)))*
     &          (XSb22*log(XSb22) + 
     &            XSb12*((-1 + (2 - XSb22)*XSb22)*log(XSb12) + 
     &               (-2 + XSb12)*XSb22*log(XSb22)))*
     &          (2*(XSt12 + XSt22 + XQ2*(-2 + log(XQ2))) - 
     &            XSt12*log(XSt12) - XSt22*log(XSt22)) + 
     &         2*(XSb12 - XSb22)*
     &          (XQ2**2*(-2 + XSb12 + XSb22)*
     &             (-1 + XQ2 - 3*XSt12) - 
     &            (-1 + XSt12)**2*
     &             (-2*XSt12 + XSb22*(1 + XSt12) + 
     &               XSb12*(1 - 2*XSb22 + XSt12)) + 
     &            XQ2*((2 - 6*XSt12)*XSt12 + 
     &               XSb22*(1 + 3*XSt12**2) + 
     &               XSb12*(1 + 3*XSt12**2 - 2*XSb22*(1 + XSt12))))
     &           *T134s(MGl2,Q2,MSt12,MGl2)) - 
     &      (1 + (-2 + XSt12)*XSt12 + XQ2*(XQ2 - 2*(1 + XSt12)))*
     &       ((2*XQ2 - XSt12 - XSt22)*
     &          (1 + (-2 + XSt22)*XSt22 + 
     &            XQ2*(XQ2 - 2*(1 + XSt22)))*
     &          (-(XSb12*(-1 + XSb22)**2*log(XSb12)**2) + 
     &            (-1 + XSb12)**2*XSb22*log(XSb22)**2) + 
     &         (-1 + XSb22)*
     &          ((-1 + XSb12)*(XSb12 - XSb22)*XSt22*
     &             (-4*(-2 + (4 - 2*XSt22)*XSt22 + 
     &                  XQ2*(1 + XQ2 + XSt22))*log(XSt22) + 
     &               (1 + XQ2 - XSt22)*(-1 + XQ2 + XSt22)*
     &                log(XSt22)**2) + 
     &            2*(-1 + XSb22)*
     &             (1 + (-2 + XSt22)*XSt22 + 
     &               XQ2*(XQ2 - 2*(1 + XSt22)))*
     &             ((XSb12 - XSt12)*T134s(MSb12,Q2,MSt12,MGl2) + 
     &               (XSb12 - XSt22)*T134s(MSb12,Q2,MSt22,MGl2) + 
     &               XQ2*(T134s(MSb12,Q2,MSt12,MGl2) + 
     &                  T134s(MSb12,Q2,MSt22,MGl2))))) + 
     &      2*(-1 + (2 - XSt12)*XSt12 + 
     &         XQ2*(-XQ2 + 2*(1 + XSt12)))*
     &       ((XSb12 - XSb22)*
     &          (XQ2**2*(-2 + XSb12 + XSb22)*
     &             (-1 + XQ2 - 3*XSt22) - 
     &            (-1 + XSt22)**2*
     &             (-2*XSt22 + XSb22*(1 + XSt22) + 
     &               XSb12*(1 - 2*XSb22 + XSt22)) + 
     &            XQ2*((2 - 6*XSt22)*XSt22 + 
     &               XSb22*(1 + 3*XSt22**2) + 
     &               XSb12*(1 + 3*XSt22**2 - 2*XSb22*(1 + XSt22))))
     &           *T134s(MGl2,Q2,MSt22,MGl2) - 
     &         (-1 + XSb12)**2*
     &          (1 + (-2 + XSt22)*XSt22 + 
     &            XQ2*(XQ2 - 2*(1 + XSt22)))*
     &          ((XSb22 - XSt12)*T134s(MSb22,Q2,MSt12,MGl2) + 
     &            (XSb22 - XSt22)*T134s(MSb22,Q2,MSt22,MGl2) + 
     &            XQ2*(T134s(MSb22,Q2,MSt12,MGl2) + 
     &               T134s(MSb22,Q2,MSt22,MGl2)))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2*
     &      (1 + (-2 + XSt12)*XSt12 + XQ2*(XQ2 - 2*(1 + XSt12)))*
     &      (1 + (-2 + XSt22)*XSt22 + XQ2*(XQ2 - 2*(1 + XSt22))))

        RuAlsTF2 = 1/24.D0*
     &   (XSb12*(-1 + XSb22)*log(XSb12)*
     &       (-6 + 3*log(XSb12) + log(XSu12) + log(XSu22)) - 
     &      (-1 + XSb12)*XSb22*log(XSb22)*
     &       (-6 + 3*log(XSb22) + log(XSu12) + log(XSu22)))/
     &    ((-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22))

        RuMGlTF2 = 1/8.D0*(XSb12*(1 + XSb12)*(-1 + XSb22)**2*
     &       log(XSb12)**2 - 
     &      (-1 + XSb12)**2*XSb22*(1 + XSb22)*log(XSb22)**2 + 
     &      (XSb12*(1 + XSb12)*(-1 + XSb22)**2*log(XSb12) + 
     &         (-1 + XSb12)*
     &          (XSb22*(2 + XSb22*(-2 + log(XSb22)) + 
     &               log(XSb22)) + 
     &            XSb12*(-2 + XSb22*(2 + (-1 - XSb22)*log(XSb22))))
     &         )*(-6 + XSu12 + XSu22 + 
     &         (-1 + XSu12)**2*(log(-1 + XSu12) - log(XSu12)) + 
     &         log(XSu12) + 
     &         (-1 + XSu22)**2*(log(-1 + XSu22) - log(XSu22)) + 
     &         log(XSu22)))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        RdAlsTF2 = 1/24.D0*
     &   (XSb12*(-1 + XSb22)*log(XSb12)*
     &       (-6 + 3*log(XSb12) + log(XSd12) + log(XSd22)) - 
     &      (-1 + XSb12)*XSb22*log(XSb22)*
     &       (-6 + 3*log(XSb22) + log(XSd12) + log(XSd22)))/
     &    ((-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22))

        RdMGlTF2 = 1/8.D0*(XSb12*(1 + XSb12)*(-1 + XSb22)**2*
     &       log(XSb12)**2 - 
     &      (-1 + XSb12)**2*XSb22*(1 + XSb22)*log(XSb22)**2 + 
     &      (XSb12*(1 + XSb12)*(-1 + XSb22)**2*log(XSb12) + 
     &         (-1 + XSb12)*
     &          (XSb22*(2 + XSb22*(-2 + log(XSb22)) + 
     &               log(XSb22)) + 
     &            XSb12*(-2 + XSb22*(2 + (-1 - XSb22)*log(XSb22))))
     &         )*(-6 + XSd12 + XSd22 + 
     &         (-1 + XSd12)**2*(log(-1 + XSd12) - log(XSd12)) + 
     &         log(XSd12) + 
     &         (-1 + XSd22)**2*(log(-1 + XSd22) - log(XSd22)) + 
     &         log(XSd22)))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        RbAlsTF2 = 1/24.D0*
     &   (-4*XSb12*(-1 + XSb22)*log(XSb12)**2 + 
     &      2*(-1 + XSb12)*XSb22*log(XSb22)*(-3 + 2*log(XSb22)) + 
     &      log(XSb12)*(-(XSb22*log(XSb22)) + 
     &         XSb12*(-6 + 6*XSb22 + log(XSb22))))/
     &    ((-1 + XSb12)*(-1 + XSb22)*(-XSb12 + XSb22))

        RbMGlTF2 = 1/8.D0*(XSb12*(1 + XSb12)*(-1 + XSb22)**2*
     &       log(XSb12)**2 - 
     &      (-1 + XSb12)**2*XSb22*(1 + XSb22)*log(XSb22)**2 + 
     &      (-6 + XSb12 + XSb22 + 
     &         (-1 + XSb12)**2*(log(-1 + XSb12) - log(XSb12)) + 
     &         log(XSb12) + 
     &         (-1 + XSb22)**2*(log(-1 + XSb22) - log(XSb22)) + 
     &         log(XSb22))*
     &       (XSb12*(1 + XSb12)*(-1 + XSb22)**2*log(XSb12) + 
     &         (-1 + XSb12)*
     &          (XSb22*(2 + XSb22*(-2 + log(XSb22)) + 
     &               log(XSb22)) + 
     &            XSb12*(-2 + XSb22*(2 + (-1 - XSb22)*log(XSb22))))
     &         ))/((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

        RtAlsTF2 = 1/24.D0*
     &   (XSb12*(-1 + XSb22)*log(XSb12)*
     &       (-6 + 3*log(XSb12) + log(XSt12) + log(XSt22)) - 
     &      (-1 + XSb12)*XSb22*log(XSb22)*
     &       (-6 + 3*log(XSb22) + log(XSt12) + log(XSt22)))/
     &    ((-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22))

        RtMGlTF2 = 1/8.D0*(XSb12*(1 + XSb12)*(-1 + XSb22)**2*
     &       log(XSb12)**2 - 
     &      (-1 + XSb12)**2*XSb22*(1 + XSb22)*log(XSb22)**2 + 
     &      2*(-1 + XSb12)*(XSb12 - XSb22)*(-1 + XSb22)*
     &       (-6 + XSt12 + XSt22 - 
     &         (1 + XQ2 - XSt12)*
     &          (XmSt12*log(-1 + XmSt12) - XmSt12*log(XmSt12) + 
     &            XpSt12*(log(-1 + XpSt12) - log(XpSt12))) - 
     &         (1 + XQ2 - XSt22)*
     &          (XmSt22*log(-1 + XmSt22) - XmSt22*log(XmSt22) + 
     &            XpSt22*(log(-1 + XpSt22) - log(XpSt22))) - 
     &         2*XQ2*(1 + log(XQ2)) + 
     &         (1 + XQ2)*(log(XSt12) + log(XSt22))) + 
     &      ((-1 - XSb22)*XSb22*log(XSb22) + 
     &         XSb12*((1 + (-2 + XSb22)*XSb22)*log(XSb12) + 
     &            2*XSb22*(1 + XSb22)*log(XSb22) + 
     &            XSb12*((1 + (-2 + XSb22)*XSb22)*log(XSb12) + 
     &               (-1 - XSb22)*XSb22*log(XSb22))))*
     &       (-6 + XSt12 + XSt22 - 
     &         XmSt12*(1 + XQ2 - XSt12)*
     &          (log(-1 + XmSt12) - log(XmSt12)) - 
     &         XmSt22*(1 + XQ2 - XSt22)*
     &          (log(-1 + XmSt22) - log(XmSt22)) - 
     &         XpSt12*(1 + XQ2 - XSt12)*
     &          (log(-1 + XpSt12) - log(XpSt12)) - 
     &         XpSt22*(1 + XQ2 - XSt22)*
     &          (log(-1 + XpSt22) - log(XpSt22)) - 
     &         2*XQ2*(1 + log(XQ2)) + 
     &         (1 + XQ2)*(log(XSt12) + log(XSt22))))/
     &    ((-1 + XSb12)**2*(XSb12 - XSb22)*(-1 + XSb22)**2)

	RCtCA = RAlsCA2 + RMGlCA2
	RCtCF = RAlsCF2 + RMSb12 + RMSb22
	RuCtTF = RuAlsTF2 + RuMGlTF2
	RdCtTF = RdAlsTF2 + RdMGlTF2
	RbCtTF = RbalsTF2 + RbMGlTF2
	RtCtTF = RtAlsTF2 + RtMGlTF2

	bo = -2*TLdb(MSb12, MSb22, MGl2)

	RCA = R12 + .25D0*R22 + R32 + R52 + R62 + RCtCA
	RCF = -.5D0*(R22 + R42 + R72 + R92) - 2*R52 + R82 + RCtCF
	RTF = Nu*(Ru102 + RuCtTF) + Nd*(Rd102 + RdCtTF)
	RTFb = Rb102 + RbCtTF
	RTFt = Rt102 + RtCtTF + bo*log(MGl2/MT2)/(3*fnorm)

	Rqcd = Re(C_A*RCA + C_F*RCF + T_F*(RTF + RTFb + RTFt))/bo
	anomalous = -C_F/4D0
	finscale = (11*C_A - 4*T_F*Nf)/12D0*log(scale**2/MGl2)

	Deltab2L_qcd = Rqcd*fnorm + anomalous + finscale
	end
 

************************************************************************

	RealType function TLdb(x, y, z)
	implicit none
	RealType x, y, z

	RealType n, d

	n = 0
	d = 1

	if( x .ne. y ) then
	  if( x*y .ne. 0 ) n = n + x*y*log(x/y)
	  d = d*(x - y)
	endif
	if( x .ne. z ) then
	  if( x*z .ne. 0 ) n = n + x*z*log(z/x)
	  d = d*(x - z)
	endif
	if( y .ne. z ) then
	  if( y*z .ne. 0 ) n = n + y*z*log(y/z)
	  d = d*(y - z)
	endif

	TLdb = n/d
	end

