c#define HHHTEST 1
c#define HFFTEST 1

* Couplings.F
* evaluation of Higgs couplings
* this file is part of FeynHiggs
* last modified 18 Dec 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"
#include "FHCouplings.h"


	subroutine FHCouplings(error,
     &    couplings_, couplingsms_, gammas_, gammasms_, fast)
	implicit none
	integer error, fast
	ComplexType couplings_(ncouplings)
	ComplexType couplingsms_(ncouplingsms)
	RealType gammas_(ngammas)
	RealType gammasms_(ngammasms)

#include "FH.h"
#define __SUBROUTINE__ "FHCouplings"

	integer i

	call CheckCouplings(error, fast)

	if( error .gt. 0 ) return

	do i = 1, ncouplings
	  couplings_(i) = couplings(i)
	enddo
	do i = 1, ncouplingsms
	  couplingsms_(i) = couplingsms(i)
	enddo
	do i = 1, ngammas
	  gammas_(i) = gammas(i)
	enddo
	do i = 1, ngammasms
	  gammasms_(i) = gammasms(i)
	enddo
	end


************************************************************************

	subroutine CheckCouplings(error, fast)
	implicit none
	integer error, fast

#include "Decay.h"

	integer i, f, h, h1, h2, hv
	integer t, tc, bc, ty, by, tk, dteff, dbReff
	integer g, g1, g2, s1, s2
	integer c1, c2, n1, n2, uz
	integer Nc
	RealType fac, sqme(2)
	RealType kw, kt, kx, x, y, rtb(3)
	RealType r, rr, r2, logrr, logr2
	RealType mhi, mlo, m, gam, qcd
	ComplexType DD(3,3,0:2)
	ComplexType cL(0:2), cR(0:2)
	ComplexType chffL(3,0:1), chhv(3,0:2), chhh(3,4,4,0:2)
	ComplexType chss(2,2), chssy(2,2)
	ComplexType u1, u2, u3, u4, c
	RealType ga(4), gasm(4), gatot_top

	RealType Kallen, KallenLR, HXVWidth, Li2
	external Kallen, KallenLR, HXVWidth, Li2
	external hee_SquaredME, heeSM_SquaredME
	external huu_SquaredME, huuSM_SquaredME
	external hdd_SquaredME, hddSM_SquaredME

	character*2 higgs(4)
	data higgs /"h0", "HH", "A0", "Hp"/

	ComplexType xhffL, xhffR
	xhffL(h,i,t,g,uz) =
     &    XHiggs(h,1,uz)*chffL(1,i)*HffDb(0,1,t,g) +
     &    XHiggs(h,2,uz)*chffL(2,i)*HffDb(0,2,t,g) +
     &    XHiggs(h,3,uz)*chffL(3,i)*HffDb(0,3,t,g)
	xhffR(h,i,t,g,uz) =
     &    XHiggs(h,1,uz)*chffL(1,i)*Conjugate(HffDb(0,1,t,g)) +
     &    XHiggs(h,2,uz)*chffL(2,i)*Conjugate(HffDb(0,2,t,g)) -
     &    XHiggs(h,3,uz)*chffL(3,i)*Conjugate(HffDb(0,3,t,g))

	error = 0
	if( coup_valid .eq. ieor(valid, fast) ) return

	call CheckHiggsCorr(error)
	if( error .gt. 0 ) return

	do i = 1, ncouplings
	  couplings(i) = 0
	enddo
	do i = 1, ncouplingsms
	  couplingsms(i) = 0
	enddo
	do i = 1, ngammas
	  gammas(i) = 0
	enddo
	do i = 1, ngammasms
	  gammasms(i) = 0
	enddo
	do i = 1, H0FF(3,4,3,3)
	  ratios(i) = 0
	enddo

* We first compute DD(1,2) = (U,Z)Higgs.Y, where
* Y = {{CA, -SA, 0}, {SA, CA, 0}, {0, 0, 1}}
* Note that SA & CA are the sine and cosine of alpha_*tree*.

* for checking with the MSSM.mod couplings
#if 0
#define YHiggs(i,j,uz) Delta(i,j)
#else
#define YHiggs(i,j,uz) XHiggs(i,j,uz)
#endif

#if defined(HHHTEST) || defined(HFFTEST)
	do uz = 1, 2
	  XHiggs(1,1,uz) = 1
	  XHiggs(2,2,uz) = 1
	  XHiggs(3,3,uz) = 1
	  XHiggs(2,1,uz) = 0
	  XHiggs(3,1,uz) = 0
	  XHiggs(1,2,uz) = 0
	  XHiggs(3,2,uz) = 0
	  XHiggs(1,3,uz) = 0
	  XHiggs(2,3,uz) = 0
	enddo
#endif

	do uz = 0, 2
	  DD(1,1,uz) = CA*YHiggs(1,1,uz) + SA*YHiggs(1,2,uz)
	  DD(2,1,uz) = CA*YHiggs(2,1,uz) + SA*YHiggs(2,2,uz)
	  DD(3,1,uz) = CA*YHiggs(3,1,uz) + SA*YHiggs(3,2,uz)

	  DD(1,2,uz) = CA*YHiggs(1,2,uz) - SA*YHiggs(1,1,uz)
	  DD(2,2,uz) = CA*YHiggs(2,2,uz) - SA*YHiggs(2,1,uz)
	  DD(3,2,uz) = CA*YHiggs(3,2,uz) - SA*YHiggs(3,1,uz)

	  DD(1,3,uz) = YHiggs(1,3,uz)
	  DD(2,3,uz) = YHiggs(2,3,uz)
	  DD(3,3,uz) = YHiggs(3,3,uz)
	enddo


*** COUPLINGS AND DECAY WIDTHS ***

c	call setdelta(100D0)
	Divergence = getdelta()
	call clearcache

* For couplings where the order is essential, the following
* conventions apply:
* - A charged Higgs at pos. 1 is always the particle, i.e. H-,
*   this fixes all other particles by charge conservation.
* - For the H0ChaCha coupling, the first Chargino is the
*   particle, chi-, and the second one the antiparticle, chi+.

	chffL(1,0) = CA/SB
	chffL(2,0) = SA/SB
	chffL(3,0) = cI/TB

	chffL(1,1) = -SA/CB
	chffL(2,1) = CA/CB
	chffL(3,1) = cI*TB

* rtb = Xtt^*/Xbb where
* htt = CA/SB;	hbb = -SA/CB
* Htt = SA/SB;	Hbb = CA/CB
* Att = I/TB;	Abb = I TB
	rtb(1) = -CA/(SA*TB)
	rtb(2) = SA/(CA*TB)
	rtb(3) = -1/TB2

	do g = 1, 3
* 1103.1335 Eq. (34)
	  x = Im(Deltaf(bTR,g))/(1 + Re(Deltaf(bTR,g)))
	  y = Re(Deltaf(bTR,g)) + x*Im(Deltaf(bTR,g))
	  do h = 1, 3
	    HffDb(0,h,2,g) = 1
	    HffDb(1,h,2,g) = 0
	    HffDb(0,h,3,g) = 1
	    HffDb(1,h,3,g) = 0

* Hbb resummation factors a la hep-ph/9912516 Eq. (21)
c	    HffDb(0,h,4,g) = (1 + rtb(h)*Deltaf(bTR,g))/abs(1 + Deltaf(bTR,g))
* Hbb resummation factors a la 1103.1335 Eqs. (113)-(116)
	    HffDb(0,h,4,g) = (1 + rtb(h)*y - cI*(1 - rtb(h))*x)/(1 + y)
	    HffDb(1,h,4,g) = Deltaf(4,g)*(1 - rtb(h))
	  enddo
	enddo

	do i = 1, 3*3
	  hggU_flat(i) = 1
	  hggDRe_flat(i) = 1
	  hggDIm_flat(i) = 1
	enddo

	dbReff = Delta(mfeff, 1)*(bBR - 4)
	dteff = Delta(mfeff, 2)*(tH - 3)

	f = -ibits(ior(fast, -fast), 31, 1)

* H0 -> Vector Vector

	do h = 1, 3
	  m = MHiggs(h)
	  call SetDecayScale(error, m, -1D0, 1)
	  if( error .ne. 0 ) return

* leading correction factors in O(alpha_s)
	  hgagaQ = 1 - AlfasMH/pi
	  hgagaSq = 1 + 8/(3*pi)*AlfasMH

	  hggU(1,3) = hggSq(1)
	  hggU(2,3) = hggSq(2)
	  hggU(3,3) = hggSq(3)

	  hggDRe(1,3) = hggSq(1)
	  hggDRe(2,3) = hggSq(2)
	  hggDRe(3,3) = hggSq(3)

	  hggDIm(1,3) = hggSq(1)
	  hggDIm(2,3) = hggSq(2)
	  hggDIm(3,3) = hggSq(3)

* H0 -> gamma gamma

	  call Kin1to2(fac, h, 0D0, 0D0, .5D0)

	  call hgagaSM_SquaredME(sqme)
	  GammaSM(H0VV(h,1)) = fac*sqme(2)
	  GammaSMTot(h) = GammaSMTot(h) + GammaSM(H0VV(h,1))

	  call hgaga_SquaredME(sqme)
	  Gamma(H0VV(h,1)) = fac*sqme(2)
	  GammaTot(h) = GammaTot(h) + Gamma(H0VV(h,1))

* H0 -> gamma Z

#if 1
	  call Kin1to2(fac, h, 0D0, MZ, 1D0)

	  call hgaZSM_CalcRenConst
	  call hgaZSM_SquaredME(sqme)
	  GammaSM(H0VV(h,2)) = fac*sqme(2)
	  GammaSMTot(h) = GammaSMTot(h) + GammaSM(H0VV(h,2))

	  call hgaZ_CalcRenConst
	  call hgaZ_SquaredME(sqme)
	  Gamma(H0VV(h,2)) = fac*sqme(2)
	  GammaTot(h) = GammaTot(h) + Gamma(H0VV(h,2))
#endif

* H0 -> Z Z (3), W W (4)

#ifndef HVV
#define HVV NLO
#endif
#define LO 0
#define NLO 1
#define OLD 2

	  do uz = 0, 2
	    cL(uz) = SB*DD(h,1,uz) + CB*DD(h,2,uz)
	  enddo

	  CouplingSM(H0VV(h,3)) = cI*EL0*MW/(SW*CW2)
	  Coupling(H0VV(h,3)) = CouplingSM(H0VV(h,3))*cL(uzint)

	  CouplingSM(H0VV(h,4)) = cI*EL0*MW/SW
	  Coupling(H0VV(h,4)) = CouplingSM(H0VV(h,4))*cL(uzint)

	  ratios(H0VV(h,3)) = Sq(cL(uzext))
	  ratios(H0VV(h,4)) = Sq(cL(uzext))

#if HVV == LO
#include "NHiggsDecayFitsLO.f"
#elif HVV == NLO
#include "NHiggsDecayFitsNLO.f"
#else
	  GammaSM(H0VV(h,3)) = .5D0*
     &      HXVWidth(m, MZ, MZ, GammaZ, .8D0, 1D0, 1)
	  GammaSM(H0VV(h,4)) =
     &      HXVWidth(m, MW, MW, GammaW, 1D0, 1D0, 1)
#endif
	  GammaSMTot(h) = GammaSMTot(h) + GammaSM(H0VV(h,3))
	  GammaSMTot(h) = GammaSMTot(h) + GammaSM(H0VV(h,4))

	  Gamma(H0VV(h,3)) = GammaSM(H0VV(h,3))*ratios(H0VV(h,3))
	  GammaTot(h) = GammaTot(h) + Gamma(H0VV(h,3))

	  Gamma(H0VV(h,4)) = GammaSM(H0VV(h,4))*ratios(H0VV(h,4))
	  GammaTot(h) = GammaTot(h) + Gamma(H0VV(h,4))

* H0 -> gluon gluon

	  call Kin1to2(fac, h, 0D0, 0D0, .5D0)

*#include "kfactor-gghEW.F"

	  call hggSM_SquaredME(sqme)
	  GammaSM(H0VV(h,5)) = fac*sqme(2)
	  GammaSMTot(h) = GammaSMTot(h) + GammaSM(H0VV(h,5))

	  call hgg_SquaredME(sqme)
	  Gamma(H0VV(h,5)) = fac*sqme(2)
	  GammaTot(h) = GammaTot(h) + Gamma(H0VV(h,5))

* H0 -> Fermion Fermion

	  do t = 2, 4
	    i = iand(t - 1, 1)
* note: for mfeff = 0, t = 4, we take HffDb(tc = 3) (which is 1)
* to omit the Db resummation
	    tc = t - ifB(t, Delta(mfeff, 0))
	    ty = t + ifQ(t, dteff)
	    do g = 1, 3
	      c = -.5D0*EL0/(SW*MW)*Mf(ty,g)*cI
	      LCouplingSM(H0FF(h,t,g,g)) = c
	      RCouplingSM(H0FF(h,t,g,g)) = c
	      LCoupling(H0FF(h,t,g,g)) = c*xhffL(h,i,tc,g, uzint)
	      RCoupling(H0FF(h,t,g,g)) = c*xhffR(h,i,tc,g, uzint)

	      cL(1) = xhffL(h,i,t,g, uzext)
	      cR(1) = xhffR(h,i,t,g, uzext)
	      ratios(H0FF(h,t,g,g)) = .5D0*(Sq(cL(1)) + Sq(cR(1)))
	    enddo
	  enddo

	  call hee_CalcRenConst
	  call heeSM_CalcRenConst
	  call huu_CalcRenConst
	  call huuSM_CalcRenConst
	  call hdd_CalcRenConst
	  call hddSM_CalcRenConst

	  do g2 = 1, 3
	    do g1 = 1 + iand(f, g2 - 1), 3 + iand(f, g2 - 3)
	      call HffDecays(h, 2, g1, g2,
     &          hee_SquaredME, heeSM_SquaredME)
	      call HffDecays(h, 3, g1, g2,
     &          huu_SquaredME, huuSM_SquaredME)
	      call HffDecays(h, 4, g1, g2,
     &          hdd_SquaredME, hddSM_SquaredME)
	    enddo
	  enddo

* H0 -> Sfermion Sfermion

	  do t = 1, 4
	    Nc = Ncolor(t)
	    tk = t + ifB(t, bBR - 4)
	    tc = t + ifQ(t, dteff) + ifB(t, dbReff)
	    ty = t + ifQ(t, tT - 3) + ifB(t, bHR - bTR)
	    do g = 1, 3
	      call H0SfSfCouplings(chss, X5(h,tc,t,g,uzint), DD)
	      call H0SfSfCouplings(chssy, X5(h,ty,t,g,uzext), DD)
	      do s2 = 1, 2
	        do s1 = 1, 2
	          Coupling(H0SfSf(h,s1,s2,t,g)) = chss(s1,s2)
	          Gamma(H0SfSf(h,s1,s2,t,g)) =
     &              Nc/(16*pi*MHiggs(h))*
     &              Kallen(MHiggs(h), MSf(s1,tk,g), MSf(s2,tk,g))*
     &              Sq(chssy(s1,s2))
	          GammaTot(h) = GammaTot(h) + Gamma(H0SfSf(h,s1,s2,t,g))
	        enddo
	      enddo
	    enddo
	    if( t .eq. 3 )
     &        chSt1St1(H0SfSf(h,1,1,3,3)) = Sq(chssy(1,1))
	  enddo

* enddo for h = 1..3
	enddo

	m = MHiggs(4)
	call SetDecayScale(error, m, -1D0, 1)
	if( error .ne. 0 ) return

	c = cI/sqrt2*EL0/(SW*MW)

	mlo = MT - .5D0
	mhi = MT + MB + .8D0

* Hp -> Fermion Fermion

	do t = 2, 4, 2
	  f = -ibits(t, 1, 1)
	  Nc = Ncolor(t)
c	  i = -iand(isB(t), isMFV(t))
	  i = -isB(t)
	  tc = t + ifQ(t, dteff)
	  bc = tc + iand(i, dbReff)
	  ty = t + ifB(t, bH - 4)
	  by = ty + iand(i, bHR - bH)
	  do g2 = 1, 3
	    do g1 = 1 + iand(f, g2 - 1), 3 + iand(f, g2 - 3)
	      cL(1) = c*TB
	      cR(1) = c/TB
	      qcd = 1
	      if( t .eq. 4 ) then
	        cL(1) = cL(1)*Conjugate(CKMin(g1,g2))
	        cR(1) = cR(1)*Conjugate(CKMin(g1,g2))
	        if( g1 .eq. g2 ) then
	          qcd = 1 + AlfasMH/pi*(C_F*9/4D0 + 8/3D0)
	          if( g1 .eq. 3 ) m = max(m, mhi)
	        endif
	      endif

	      LCoupling(HpFF(t/2,g1,g2)) = Mf(bc,g2)*cL(1)
	      RCoupling(HpFF(t/2,g1,g2)) = Mf(tc-1,g1)*cR(1)

	      Gamma(HpFF(t/2,g1,g2)) = Nc/(8*pi)*m*qcd*
     &          KallenLR(m, Mf(t-1,g1), Mf(t,g2),
     &            Mf(by,g2)*cL(1), Mf(ty-1,g1)*cR(1))
	      GammaTot(4) = GammaTot(4) + Gamma(HpFF(t/2,g1,g2))
	    enddo
	  enddo
	enddo

	m = MHiggs(4)
	if( m .lt. mhi ) then
	  GammaTot(4) = GammaTot(4) - Gamma(HpFF(2,3,3))
	  gam = 0
	  if( m .gt. MW + 2*MB ) then
	    m = min(m, mlo)
* off-shell decay to top* bottom -> W bb
* taken from hep-ph/9511342, Eq. (65)
	    kw = MW2/m**2
	    kt = MT2/m**2
	    kx = (kt - 1)/(kt - kw)
	    gam = 3/(64*pi**3)*(GF*MT2/TB)**2*m*
     &        ( -2.5D0 + (4 - 1.5D0*kw)*kw +
     &          (1 - kw)/kt**2*(2*kw**2*(2 - kt) + kt*(3*kt**2 - kw)) +
     &          (1 - 4*kt + 3*(kt - kw)*(kt + kw))*log(kx) +
     &          kw**2/kt**3*(3*kt + 4*kw*(kt - 1))*log(kx*kw) )*
     &        (1 - MB/MT)**3
*             ^^^^^^^^^^^^^^ fudge factor for non-zero bottom mass
* interpolate in the region MT - .5 GeV < MHp < MT + MB + .8 GeV
	    if( MHiggs(4) .gt. mlo ) gam = gam +
     &        (MHiggs(4) - mlo)/(mhi - mlo)*
     &        (Gamma(HpFF(2,3,3)) - gam)
	  endif
	  Gamma(HpFF(2,3,3)) = gam
	  GammaTot(4) = GammaTot(4) + Gamma(HpFF(2,3,3))
	endif

* Higgs -> Higgs Vector

	do uz = 0, 2
	  do hv = 1, 3
	    chhv(hv,uz) = CB*DD(hv,1,uz) - SB*DD(hv,2,uz)
	  enddo
	enddo

* Higgs -> Higgs Z

	do hv = 1, 3
	  do h = 1, 3
	    do uz = 0, 2
	      cL(uz) = DD(h,3,uz)*chhv(hv,uz) - DD(hv,3,uz)*chhv(h,uz)
	    enddo

	    Coupling(H0HV(h,hv)) = -EL0/(2*SW*CW)*cL(uzint)

	    Gamma(H0HV(h,hv)) =
     &        HXVWidth(MHiggs(h), MHiggs(hv), MZ, GammaZ,
     &          1.8D0, 2.4D0, 0)*
     &        Sq(cL(uzext))
	    GammaTot(h) = GammaTot(h) + Gamma(H0HV(h,hv))
	  enddo
	enddo

* Hp -> H0 W

	do hv = 1, 3
	  do uz = 0, 2
	    cL(uz) = chhv(hv,uz) + cI*DD(hv,3,uz)
	  enddo

	  Coupling(HpHV(hv)) = -cI*EL0/(2*SW)*cL(uzint)

	  Gamma(HpHV(hv)) =
     &      HXVWidth(MHiggs(4), MHiggs(hv), MW, GammaW,
     &        1.7D0, 2.2D0, 0)*
     &      Sq(cL(uzext))
	  GammaTot(4) = GammaTot(4) + Gamma(HpHV(hv))
	enddo

* H0 -> Higgs Higgs

	call H0HHCouplings(chhh, DD)
	call hhh_CalcRenConst

	do h2 = 1, 3
	  do h1 = 1, 3
	    do h = 1, 3
	      Coupling(H0HH(h,h1,h2)) = chhh(h,h1,h2,uzint)
	      call HhhDecays(h, h1, h2)
	    enddo
	  enddo
	enddo

	do h = 1, 3
	  Coupling(H0HH(h,4,4)) = chhh(h,4,4,uzint)
	  Gamma(H0HH(h,4,4)) = 1/(32*pi*MHiggs(h))*
     &      Kallen(MHiggs(h), MHiggs(4), MHiggs(4))*
     &      Sq(chhh(h,4,4,uzext))
	  GammaTot(h) = GammaTot(h) + Gamma(H0HH(h,4,4))
	enddo

* Hp -> Sfermion Sfermion

	do t = 2, 4, 2
	  f = -ibits(t, 1, 1)
	  Nc = Ncolor(t)
c	  i = -iand(isB(t), isMFV(t))
	  i = -isB(t)
	  tc = t + ifQ(t, dteff) + iand(i, dbReff)
	  ty = t + ifB(t, bH - 4)
	  by = ty + iand(i, bHR - bH)
	  tk = t + iand(i, bBR - 4)
	  do g2 = 1, 3
	    do g1 = 1 + iand(f, g2 - 1), 3 + iand(f, g2 - 3)
	      call HpSfSfCouplings(chss, X5(tc,t,t,g1,g2))
	      call HpSfSfCouplings(chssy, X5(by,ty,t,g1,g2))
	      do s2 = 1, 2
	        do s1 = 1, 2
	          Coupling(HpSfSf(s1,s2,t/2,g1,g2)) = chss(s1,s2)
	          Gamma(HpSfSf(s1,s2,t/2,g1,g2)) =
     &              Nc/(16*pi*MHiggs(4))*
     &              Kallen(MHiggs(4), MSf(s1,t-1,g1), MSf(s2,tk,g2))*
     &              Sq(chssy(s1,s2))
	          GammaTot(4) = GammaTot(4) +
     &              Gamma(HpSfSf(s1,s2,t/2,g1,g2))
	        enddo
	      enddo
	    enddo
	  enddo
	enddo

* H0 -> Chargino Chargino

	do c2 = 1, 2
          do c1 = 1, 2
	    u1 = UChaC(c1,2)*VChaC(c2,1)
	    u2 = UChaC(c1,1)*VChaC(c2,2)
	    u3 = VCha(c1,1)*UCha(c2,2)
	    u4 = VCha(c1,2)*UCha(c2,1)
	    do h = 1, 3
	      do uz = 0, 2
	        cL(uz) = 1/sqrt2*EL0/SW*(
     &            DD(h,3,uz)*(-u1*SB - u2*CB) +
     &            DD(h,2,uz)*(-cI*u1) +
     &            DD(h,1,uz)*(-cI*u2) )
	        cR(uz) = 1/sqrt2*EL0/SW*(
     &            DD(h,3,uz)*(u3*SB + u4*CB) +
     &            DD(h,2,uz)*(-cI*u3) +
     &            DD(h,1,uz)*(-cI*u4) )
	      enddo

	      LCoupling(H0ChaCha(h,c1,c2)) = cL(uzint)
	      RCoupling(H0ChaCha(h,c1,c2)) = cR(uzint)

	      Gamma(H0ChaCha(h,c1,c2)) = MHiggs(h)/(8*pi)*
     &          KallenLR(MHiggs(h), MCha(c1), MCha(c2),
     &            cL(uzext), cR(uzext))
	      GammaTot(h) = GammaTot(h) + Gamma(H0ChaCha(h,c1,c2))
	    enddo
	  enddo
	enddo

* H0 -> Neutralino Neutralino

	do n2 = 1, 4
          do n1 = 1, 4
	    u1 = ZNeu(n1,3)*(CW*ZNeu(n2,2) - SW*ZNeu(n2,1)) +
     &           ZNeu(n2,3)*(CW*ZNeu(n1,2) - SW*ZNeu(n1,1))
	    u2 = ZNeu(n1,4)*(CW*ZNeu(n2,2) - SW*ZNeu(n2,1)) +
     &           ZNeu(n2,4)*(CW*ZNeu(n1,2) - SW*ZNeu(n1,1))
	    do h = 1, 3
	      do uz = 0, 2
	        cL(uz) = .5D0*EL0/(CW*SW)*(
     &            DD(h,3,uz)*(-Conjugate(u1)*SB + Conjugate(u2)*CB) +
     &            DD(h,2,uz)*(-cI*Conjugate(u1)) +
     &            DD(h,1,uz)*( cI*Conjugate(u2)) )
	        cR(uz) = .5D0*EL0/(CW*SW)*(
     &            DD(h,3,uz)*( u1*SB - u2*CB) +
     &            DD(h,2,uz)*(-cI*u1) +
     &            DD(h,1,uz)*( cI*u2) )
	      enddo

	      LCoupling(H0NeuNeu(h,n1,n2)) = cL(uzint)
	      RCoupling(H0NeuNeu(h,n1,n2)) = cR(uzint)

	      Gamma(H0NeuNeu(h,n1,n2)) = MHiggs(h)/(16*pi)*
     &          KallenLR(MHiggs(h), MNeu(n1), MNeu(n2),
     &            cL(uzext), cR(uzext))
	      GammaTot(h) = GammaTot(h) + Gamma(H0NeuNeu(h,n1,n2))
	      if( n1 .ne. n2 )
     &          Gamma(H0NeuNeu(h,n1,n2)) = 2*Gamma(H0NeuNeu(h,n1,n2))
	    enddo
	  enddo
	enddo

* Higgs -> Neutralino Chargino

	do c2 = 1, 2
	  do n1 = 1, 4
	    cL(1) = -cI*EL0/SW*CB*(ZNeuC(n1,4)*VChaC(c2,1) +
     &        1/sqrt2*(ZNeuC(n1,2) + ZNeuC(n1,1)*SW/CW)*VChaC(c2,2))
	    cR(1) = -cI*EL0/SW*SB*(ZNeu(n1,3)*UCha(c2,1) -
     &        1/sqrt2*(ZNeu(n1,2) + ZNeu(n1,1)*SW/CW)*UCha(c2,2))

	    LCoupling(HpNeuCha(n1,c2)) = cL(1)
	    RCoupling(HpNeuCha(n1,c2)) = cR(1)

	    Gamma(HpNeuCha(n1,c2)) = MHiggs(4)/(8*pi)*
     &        KallenLR(MHiggs(4), MNeu(n1), MCha(c2),
     &          cL(1), cR(1))
	    GammaTot(4) = GammaTot(4) + Gamma(HpNeuCha(n1,c2))
	  enddo
	enddo

* top -> W bottom

	r = MW/MT
	rr = (1 - r)*(1 + r)
	r2 = 1 + 2*r**2
	logrr = log(rr)
	logr2 = log(r**2)

* hep-ph/0408158, Eqs 6 + 26:
	Gamma(tBF(1)) = GF*MT*MT2/(8*pi*sqrt2)*rr**2*r2*(
     &    1 - C_F/(2*pi)*AlfasMT*(
     &      2/3D0*pi**2 + 1/2D0 + 2*(logr2 + 1)*(logrr - 1) -
     &      4/(3*rr) + (22 - 34*r**2)/(9*rr**2)*logr2 +
     &      (3 + 27*logrr - 4*logr2)/(9*r2) + 4*Li2(r**2) ) )

* top -> H+ bottom

	if( MT .gt. MHiggs(4) + MB ) then
	  r = 1 - MHiggs2(4)/MT2

* these is the H+ -> t b coupling with MB(MT)/|1 + Db|, MT(MT)
	  c = cI/sqrt2*EL0/(SW*MW)
	  cL(1) = Mf(bTR,3)*TB*c
	  cR(1) = Mf(tT,3)/TB*c

* Corresponds to hep-ph/0211098, Eq. (7a), but with Delta_b
* corrections included (via the couplings).
* Here only a factor of 1/MT is used due to different definition
* of the Kallen function.  Additional factor 1/2 from mismatch of
* 'c' in  Couplings.F and \tilde c.
* Finite O(alpha_s) corrections from hep-ph/9912516, Eq. (29).
	  Gamma(tBF(2)) = 1/(32*pi*MT) *
     &      Kallen(MT, MHiggs(4), MB) * (
     &        ((MT - MHiggs(4))*(MT + MHiggs(4)) + MB2)*
     &          (Sq(cL(1)) + Sq(cR(1))) +
     &        4*MT*MB*abs(cL(1)*cR(1)) ) *
     &      ( 1 + AlfasMT/pi*(7 - 8/9D0*pi**2 +
     &          (2/3D0*r**2 - 2)*log(r) + r*(2 + 4/9D0*r)) )
	endif

*** BRANCHING RATIOS ***

	do h = 1, 4
	  if( GammaTot(h) .gt. 0 ) then
	    ga(h) = 1/GammaTot(h)
	  else
	    ga(h) = 0
	  endif
	  if( GammaSMTot(h) .gt. 0 ) then
	    gasm(h) = 1/GammaSMTot(h)
	  else
	    gasm(h) = 0
	  endif
	enddo

	h = 0

*#define Next(h) mod(h, 3) + 1
#define Next(h) (iand(h + 1, 3) + ishft(h + 1, -2))

	do i = H0VV(1,1), H0VV(3,5)
	  h = Next(h)
	  BR(i) = Gamma(i)*ga(h)
	  BRSM(i) = GammaSM(i)*gasm(h)
	enddo

	do i = H0FF(1,2,1,1), H0FF(3,4,3,3)
	  h = Next(h)
	  BR(i) = Gamma(i)*ga(h)
	  BRSM(i) = GammaSM(i)*gasm(h)
	enddo

	do i = H0ChaCha(1,1,1), H0ChaCha(3,2,2)
	  h = Next(h)
	  BR(i) = Gamma(i)*ga(h)
	enddo

	do i = H0NeuNeu(1,1,1), H0NeuNeu(3,4,4)
	  h = Next(h)
	  BR(i) = Gamma(i)*ga(h)
	enddo

	do i = H0HV(1,1), H0HV(3,3)
	  h = Next(h)
	  BR(i) = Gamma(i)*ga(h)
	enddo

	do i = H0HH(1,1,1), H0HH(3,4,4)
	  h = Next(h)
	  BR(i) = Gamma(i)*ga(h)
	enddo

	do i = H0SfSf(1,1,1,1,1), H0SfSf(3,2,2,4,3)
	  h = Next(h)
	  BR(i) = Gamma(i)*ga(h)
	enddo

	do i = HpFF(1,1,1), HpFF(2,3,3)
	  BR(i) = Gamma(i)*ga(4)
	enddo

	do i = HpNeuCha(1,1), HpNeuCha(4,2)
	  BR(i) = Gamma(i)*ga(4)
	enddo

	do i = HpHV(1), HpHV(3)
	  BR(i) = Gamma(i)*ga(4)
	enddo

	do i = HpSfSf(1,1,1,1,1), HpSfSf(2,2,2,3,3)
	  BR(i) = Gamma(i)*ga(4)
	enddo

	gatot_top = Gamma(tBF(1)) + Gamma(tBF(2))
	if( gatot_top .ne. 0 ) then
	  BR(tBF(1)) = Gamma(tBF(1))/gatot_top
	  BR(tBF(2)) = Gamma(tBF(2))/gatot_top
	endif

	coup_valid = ieor(valid, fast)

	if( debuglevel .lt. 1 ) return

#define Header(s) s, " decays    Width                    WidthSM"
#define HiggsHeader(h) Header(higgs(h))

	DCOUP "============================" ENDL
	do h = 1, 3
	  DHEAD higgs(h), " total width =", GammaTot(h) ENDL
	  DHEAD HiggsHeader(h) ENDL
	  DCOUP "gamma-gamma", Gamma(H0VV(h,1)), GammaSM(H0VV(h,1)) ENDL
	  DCOUP "gamma-Z    ", Gamma(H0VV(h,2)), GammaSM(H0VV(h,2)) ENDL
	  DCOUP "Z-Z        ", Gamma(H0VV(h,3)), GammaSM(H0VV(h,3)) ENDL
	  DCOUP "W-W        ", Gamma(H0VV(h,4)), GammaSM(H0VV(h,4)) ENDL
	  DCOUP "gluon-gluon", Gamma(H0VV(h,5)), GammaSM(H0VV(h,5)) ENDL
	  DCOUP " " ENDL
	  DCOUP "h0-Z       ", Gamma(H0HV(h,1)) ENDL
	  DCOUP "HH-Z       ", Gamma(H0HV(h,2)) ENDL
	  DCOUP "A0-Z       ", Gamma(H0HV(h,3)) ENDL
	  DCOUP " " ENDL
	  DCOUP "h0-h0      ", Gamma(H0HH(h,1,1)) ENDL
	  DCOUP "Hp-Hp      ", Gamma(H0HH(h,4,4)) ENDL
	  DCOUP " " ENDL
	  DCOUP "f-f    e   ", Gamma(H0FF(h,2,1,1)),
     &                       GammaSM(H0FF(h,2,1,1)) ENDL
	  DCOUP "      mu   ", Gamma(H0FF(h,2,2,2)),
     &                       GammaSM(H0FF(h,2,2,2)) ENDL
	  DCOUP "     tau   ", Gamma(H0FF(h,2,3,3)),
     &                       GammaSM(H0FF(h,2,3,3)) ENDL
	  DCOUP "       u   ", Gamma(H0FF(h,3,1,1)),
     &                       GammaSM(H0FF(h,3,1,1)) ENDL
	  DCOUP "       c   ", Gamma(H0FF(h,3,2,2)),
     &                       GammaSM(H0FF(h,3,2,2)) ENDL
	  DCOUP "       t   ", Gamma(H0FF(h,3,3,3)),
     &                       GammaSM(H0FF(h,3,3,3)) ENDL
	  DCOUP "       d   ", Gamma(H0FF(h,4,1,1)),
     &                       GammaSM(H0FF(h,4,1,1)) ENDL
	  DCOUP "       s   ", Gamma(H0FF(h,4,2,2)),
     &                       GammaSM(H0FF(h,4,2,2)) ENDL
	  DCOUP "       b   ", Gamma(H0FF(h,4,3,3)),
     &                       GammaSM(H0FF(h,4,3,3)) ENDL
	  DCOUP " " ENDL
	  DCOUP "Cha1-Cha1  ", Gamma(H0ChaCha(h,1,1)) ENDL
	  DCOUP "Cha1-Cha2  ", Gamma(H0ChaCha(h,1,2)) ENDL
	  DCOUP "Cha2-Cha1  ", Gamma(H0ChaCha(h,2,1)) ENDL
	  DCOUP "Cha2-Cha2  ", Gamma(H0ChaCha(h,2,1)) ENDL
	  DCOUP " " ENDL
	  DCOUP "Neu1-Neu1  ", Gamma(H0NeuNeu(h,1,1)) ENDL
	  DCOUP "Neu1-Neu2  ", Gamma(H0NeuNeu(h,1,2)) ENDL
	  DCOUP "Neu1-Neu3  ", Gamma(H0NeuNeu(h,1,3)) ENDL
	  DCOUP "Neu1-Neu4  ", Gamma(H0NeuNeu(h,1,4)) ENDL
	  DCOUP "Neu2-Neu2  ", Gamma(H0NeuNeu(h,2,2)) ENDL
	  DCOUP "Neu2-Neu3  ", Gamma(H0NeuNeu(h,2,3)) ENDL
	  DCOUP "Neu2-Neu4  ", Gamma(H0NeuNeu(h,2,4)) ENDL
	  DCOUP "Neu3-Neu3  ", Gamma(H0NeuNeu(h,3,3)) ENDL
	  DCOUP "Neu3-Neu4  ", Gamma(H0NeuNeu(h,3,4)) ENDL
	  DCOUP "Neu4-Neu4  ", Gamma(H0NeuNeu(h,4,4)) ENDL
	  DCOUP " " ENDL
	  DCOUP "sf-sf ne 11", Gamma(H0SfSf(h,1,1,1,1)) ENDL
	  DCOUP "     nmu 11", Gamma(H0SfSf(h,1,1,1,2)) ENDL
	  DCOUP "    ntau 11", Gamma(H0SfSf(h,1,1,1,3)) ENDL
	  DCOUP "       e 11", Gamma(H0SfSf(h,1,1,2,1)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,2,1)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,2,1)) ENDL
	  DCOUP "      mu 11", Gamma(H0SfSf(h,1,1,2,2)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,2,2)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,2,2)) ENDL
	  DCOUP "     tau 11", Gamma(H0SfSf(h,1,1,2,3)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,2,3)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,2,3)) ENDL
	  DCOUP "       u 11", Gamma(H0SfSf(h,1,1,3,1)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,3,1)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,3,1)) ENDL
	  DCOUP "       c 11", Gamma(H0SfSf(h,1,1,3,2)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,3,2)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,3,2)) ENDL
	  DCOUP "       t 11", Gamma(H0SfSf(h,1,1,3,3)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,3,3)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,3,3)) ENDL
	  DCOUP "       d 11", Gamma(H0SfSf(h,1,1,4,1)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,4,1)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,4,1)) ENDL
	  DCOUP "       s 11", Gamma(H0SfSf(h,1,1,4,2)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,4,2)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,4,2)) ENDL
	  DCOUP "       b 11", Gamma(H0SfSf(h,1,1,4,3)) ENDL
	  DCOUP "         12", Gamma(H0SfSf(h,1,2,4,3)) ENDL
	  DCOUP "         22", Gamma(H0SfSf(h,2,2,4,3)) ENDL
	  DCOUP " " ENDL
	enddo

	DHEAD higgs(4), " total width =", GammaTot(4) ENDL
	DHEAD HiggsHeader(4) ENDL
	DCOUP "h0-W       ", Gamma(HpHV(1)) ENDL
	DCOUP "HH-W       ", Gamma(HpHV(2)) ENDL
	DCOUP "A0-W       ", Gamma(HpHV(3)) ENDL
	DCOUP " " ENDL
	DCOUP "f-f'   e   ", Gamma(HpFF(1,1,1)) ENDL
	DCOUP "      mu   ", Gamma(HpFF(1,2,2)) ENDL
	DCOUP "     tau   ", Gamma(HpFF(1,3,3)) ENDL
	DCOUP "     u-d   ", Gamma(HpFF(2,1,1)) ENDL
	DCOUP "     u-s   ", Gamma(HpFF(2,1,2)) ENDL
	DCOUP "     u-b   ", Gamma(HpFF(2,1,3)) ENDL
	DCOUP "     c-d   ", Gamma(HpFF(2,2,1)) ENDL
	DCOUP "     c-s   ", Gamma(HpFF(2,2,2)) ENDL
	DCOUP "     c-b   ", Gamma(HpFF(2,2,3)) ENDL
	DCOUP "     t-d   ", Gamma(HpFF(2,3,1)) ENDL
	DCOUP "     t-s   ", Gamma(HpFF(2,3,2)) ENDL
	DCOUP "     t-b   ", Gamma(HpFF(2,3,3)) ENDL
	DCOUP " " ENDL
	DCOUP "Neu1-Cha1  ", Gamma(HpNeuCha(1,1)) ENDL
	DCOUP "Neu1-Cha2  ", Gamma(HpNeuCha(1,2)) ENDL
	DCOUP "Neu2-Cha1  ", Gamma(HpNeuCha(2,1)) ENDL
	DCOUP "Neu2-Cha2  ", Gamma(HpNeuCha(2,2)) ENDL
	DCOUP "Neu3-Cha1  ", Gamma(HpNeuCha(3,1)) ENDL
	DCOUP "Neu3-Cha2  ", Gamma(HpNeuCha(3,2)) ENDL
	DCOUP "Neu4-Cha1  ", Gamma(HpNeuCha(4,1)) ENDL
	DCOUP "Neu4-Cha2  ", Gamma(HpNeuCha(4,2)) ENDL
	DCOUP " " ENDL
	DCOUP "sf-sf' e 11", Gamma(HpSfSf(1,1,1,1,1)) ENDL
	DCOUP "         12", Gamma(HpSfSf(1,2,1,1,1)) ENDL
	DCOUP "         22", Gamma(HpSfSf(2,2,1,1,1)) ENDL
	DCOUP "      mu 11", Gamma(HpSfSf(1,1,1,2,2)) ENDL
	DCOUP "         12", Gamma(HpSfSf(1,2,1,2,2)) ENDL
	DCOUP "         22", Gamma(HpSfSf(2,2,1,2,2)) ENDL
	DCOUP "     tau 11", Gamma(HpSfSf(1,1,1,3,3)) ENDL
	DCOUP "         12", Gamma(HpSfSf(1,2,1,3,3)) ENDL
	DCOUP "         22", Gamma(HpSfSf(2,2,1,3,3)) ENDL
	DCOUP "     u-d 11", Gamma(HpSfSf(1,1,2,1,1)) ENDL
	DCOUP "         12", Gamma(HpSfSf(1,2,2,1,1)) ENDL
	DCOUP "         22", Gamma(HpSfSf(2,2,2,1,1)) ENDL
	DCOUP "     c-s 11", Gamma(HpSfSf(1,1,2,2,2)) ENDL
	DCOUP "         12", Gamma(HpSfSf(1,2,2,2,2)) ENDL
	DCOUP "         22", Gamma(HpSfSf(2,2,2,2,2)) ENDL
	DCOUP "     t-b 11", Gamma(HpSfSf(1,1,2,3,3)) ENDL
	DCOUP "         12", Gamma(HpSfSf(1,2,2,3,3)) ENDL
	DCOUP "         22", Gamma(HpSfSf(2,2,2,3,3)) ENDL
	DCOUP " " ENDL

	DHEAD "top total width =", gatot_top ENDL
	DHEAD Header("top") ENDL
	DCOUP "W-b        ", Gamma(tBF(1)) ENDL
	DCOUP "Hp-b       ", Gamma(tBF(2)) ENDL
	DCOUP "============================" ENDL
	end


************************************************************************

	subroutine H0HHCouplings(ghhh, DD)
	implicit none
	ComplexType ghhh(3,4,4,0:2), DD(3,3,0:2)

#include "FH.h"

	integer h, h1, h2, uz
	ComplexType omega_t
	ComplexType hc(6,6,6)
* 1 = Phi2, 2 = Phi1, 3 = A, 4 = H+-, 5 = G, 6 = G+-

	call CalcSelf(Mh02, cpeven, 0)
c	call CalcSelf(0D0, cpeven, 0)
	omega_t = CA2*Re(seR(h0h0)) + SA2*Re(seR(HHHH)) +
     &    S2A*Re(seR(h0HH))

	hc(1,1,1) = -3*SB + 3/MZ2/SB*(omega_t +
     &    GF*sqrt2/(pi**2*SB2)*MT2**2)		! hep-ph/0108245, eq26
	hc(1,1,2) = CB 
	hc(1,2,2) = SB 
	hc(2,2,2) = -3*CB
	hc(1,3,3) = -SB*C2B + SB*omega_t/(MZ2*TB2)
	hc(2,3,3) = CB*C2B
#if HHHTEST
c	hc(1,1,1) = -3*SB
c	hc(1,3,3) = -SB*C2B
#endif

	do uz = 0, 2
	  do h2 = 1, 3
	    do h1 = 1, 3
	      do h = 1, 3
	        ghhh(h,h1,h2,uz) = cI*EL0*MZ/(2*CW*SW)*(
     &            hc(1,1,1)* DD(h,1,uz)*DD(h1,1,uz)*DD(h2,1,uz) +
     &            hc(1,1,2)*(DD(h,1,uz)*DD(h1,1,uz)*DD(h2,2,uz) +
     &                       DD(h,1,uz)*DD(h1,2,uz)*DD(h2,1,uz) +
     &                       DD(h,2,uz)*DD(h1,1,uz)*DD(h2,1,uz)) +
     &            hc(1,2,2)*(DD(h,1,uz)*DD(h1,2,uz)*DD(h2,2,uz) +
     &                       DD(h,2,uz)*DD(h1,1,uz)*DD(h2,2,uz) +
     &                       DD(h,2,uz)*DD(h1,2,uz)*DD(h2,1,uz)) +
     &            hc(2,2,2)* DD(h,2,uz)*DD(h1,2,uz)*DD(h2,2,uz) +
     &            hc(1,3,3)*(DD(h,1,uz)*DD(h1,3,uz)*DD(h2,3,uz) +
     &                       DD(h,3,uz)*DD(h1,1,uz)*DD(h2,3,uz) +
     &                       DD(h,3,uz)*DD(h1,3,uz)*DD(h2,1,uz)) +
     &            hc(2,3,3)*(DD(h,2,uz)*DD(h1,3,uz)*DD(h2,3,uz) +
     &                       DD(h,3,uz)*DD(h1,2,uz)*DD(h2,3,uz) +
     &                       DD(h,3,uz)*DD(h1,3,uz)*DD(h2,2,uz)) )
	      enddo
	    enddo
	  enddo
	enddo

* hep-ph/9612362, eq16
	hc(1,4,4) = SB*(C2B - 2*CW2 - omega_t/(MZ2*TB2))
	hc(2,4,4) = CB*(C2B - 2*CW2)
	hc(3,4,4) = 0

	do uz = 0, 2
	  do h = 1, 3
	    ghhh(h,4,1,uz) = 0
	    ghhh(h,4,2,uz) = 0
	    ghhh(h,4,3,uz) = 0
	    ghhh(h,1,4,uz) = 0
	    ghhh(h,2,4,uz) = 0
	    ghhh(h,3,4,uz) = 0
	    ghhh(h,4,4,uz) = cI*EL0*MZ/(2*CW*SW)*(
     &        hc(1,4,4)*DD(h,1,uz) +
     &        hc(2,4,4)*DD(h,2,uz) +
     &        hc(3,4,4)*DD(h,3,uz) )
	  enddo
	enddo

c	hc(2,6,6) = -C2B*CB
c	hc(1,6,6) =  C2B*SB
c	hc(2,5,5) = -C2B*CB
c	hc(1,5,5) =  C2B*SB
c	hc(2,3,6) =  S2B*CB
c	hc(1,3,6) = -S2B*SB
c	hc(2,4,6) = -CW2*SB + S2B*CB
c	hc(1,4,6) =  CW2*CB - S2B*SB
c	hc(3,4,6) = CW2		! (-I) is missing for (AH+G-), (+I) for AH-G+
	end


************************************************************************

	subroutine H0SfSfCouplings(ghss, h_ty_t_g_uz, DD)
	implicit none
	ComplexType ghss(2,2), DD(3,3,0:2)
	integer h_ty_t_g_uz

#include "FH.h"

	integer h, ty, t, g, uz
	ComplexType u1111, u1112, u1121, u1122, u1212
	ComplexType u1222, u2112, u2121, u2122, u2222
	RealType d1, d2, m11a, m22a, m11b, m22b, m, mx
	ComplexType m12a, m12b, m3, c
	ComplexType hc(3,3)

	h = nib4(h_ty_t_g_uz)
	ty = nib3(h_ty_t_g_uz)
	t = nib2(h_ty_t_g_uz)
	g = nib1(h_ty_t_g_uz)
	uz = nib0(h_ty_t_g_uz)

	u1111 = USf(1,1,ty,g)*USfC(1,1,ty,g)
	u1112 = USf(1,1,ty,g)*USfC(1,2,ty,g)
	u1121 = USf(1,1,ty,g)*USfC(2,1,ty,g)
	u1122 = USf(1,1,ty,g)*USfC(2,2,ty,g)
	u1212 = USf(1,2,ty,g)*USfC(1,2,ty,g)
	u1222 = USf(1,2,ty,g)*USfC(2,2,ty,g)
	u2112 = USf(2,1,ty,g)*USfC(1,2,ty,g)
	u2121 = USf(2,1,ty,g)*USfC(2,1,ty,g)
	u2122 = USf(2,1,ty,g)*USfC(2,2,ty,g)
	u2222 = USf(2,2,ty,g)*USfC(2,2,ty,g)

	d2 = Qf(t)*SW2*MZ2
	d1 = (iand(t, 1) - .5D0)*MZ2 - d2
	m11a = SB*d1
	m22a = SB*d2
	m11b = CB*d1
	m22b = CB*d2
	m = Mf(ty,g)

	if( iand(t, 1) .eq. 1 ) then
	  mx = -m/SB
	  m11a = m11a + mx*m
	  m22a = m22a + mx*m
	  m12a = mx*Conjugate(Af(t,g))
	  m12b = mx*MUE
	  m3 = -m12a*CB - m12b*SB
	else
	  mx = m/CB
	  m11b = m11b + mx*m
	  m22b = m22b + mx*m
	  m12a = mx*MUE
	  m12b = mx*Conjugate(Af(t,g))
	  m3 = m12a*CB + m12b*SB
	endif

	hc(1,1) = m11a*u1111 + m22a*u1212 + Re(m12a*u1112)
	hc(2,1) = m11b*u1111 + m22b*u1212 + Re(m12b*u1112)
	hc(3,1) = Im(m3*u1112)

	hc(1,2) = m11a*u2121 + m22a*u2222 + Re(m12a*u2122)
	hc(2,2) = m11b*u2121 + m22b*u2222 + Re(m12b*u2122)
	hc(3,2) = Im(m3*u2122)

	hc(1,3) = m11a*u1121 + m22a*u1222 +
     &    .5D0*(m12a*u1122 + Conjugate(m12a*u2112))
	hc(2,3) = m11b*u1121 + m22b*u1222 +
     &    .5D0*(m12b*u1122 + Conjugate(m12b*u2112))
	hc(3,3) = -.5D0*cI*(m3*u1122 - Conjugate(m3*u2112))

	c = cI*EL0/(MW*SW)
	ghss(1,1) = c*(DD(h,1,uz)*hc(1,1) -
     &                 DD(h,2,uz)*hc(2,1) -
     &                 DD(h,3,uz)*hc(3,1))
	ghss(2,2) = c*(DD(h,1,uz)*hc(1,2) -
     &                 DD(h,2,uz)*hc(2,2) -
     &                 DD(h,3,uz)*hc(3,2))
	ghss(2,1) = c*(DD(h,1,uz)*hc(1,3) -
     &                 DD(h,2,uz)*hc(2,3) -
     &                 DD(h,3,uz)*hc(3,3))
	ghss(1,2) = -Conjugate(ghss(2,1))

	if( DEBUGLEVEL .lt. 5 ) return

	DCOUP "H0SfSfCouplings h,t,ty,g = ",
     &    Digit(h), Digit(t), Hex(ty), Digit(g) ENDL
	DCOUP "H0SfSfCouplings Mf =", m ENDL
	end


************************************************************************

	subroutine HpSfSfCouplings(ghss, by_ty_t_g1_g2)
	implicit none
	ComplexType ghss(2,2)
	integer by_ty_t_g1_g2

#include "FH.h"

	integer by, ty, t, g1, g2, s1, s2
	RealType m11, m22
	ComplexType m12, m21, c

	by = nib4(by_ty_t_g1_g2)
	ty = nib3(by_ty_t_g1_g2)
	t = nib2(by_ty_t_g1_g2)
	g1 = nib1(by_ty_t_g1_g2)
	g2 = nib0(by_ty_t_g1_g2)

	m11 = Mf2(ty-1,g1)/TB + Mf2(by,g2)*TB - MW2*S2B
	m22 = Mf(ty-1,g1)*Mf(by,g2)*(1/TB + TB)
	m12 = Mf(by,g2)*(Conjugate(MUE) + Af(t,g2)*TB)
	m21 = Mf(ty-1,g1)*(MUE + Conjugate(Af(t-1,g1))/TB)

	c = cI/sqrt2*EL0/(MW*SW)
	if( t .eq. 4 ) then
	  c = c*Conjugate(CKMin(g1,g2))
	else
	  c = c*Delta(g1, g2)
	endif

	do s1 = 1, 2
	  do s2 = 1, 2
	    ghss(s1,s2) = c*(
     &        USfC(s1,1,ty-1,g1)*
     &          (m11*USf(s2,1,by,g2) + m12*USf(s2,2,by,g2)) +
     &        USfC(s1,2,ty-1,g1)*
     &          (m21*USf(s2,1,by,g2) + m22*USf(s2,2,by,g2)) )
	  enddo
	enddo
	end


************************************************************************
	
	RealType function Kallen(m, m1, m2)
	implicit none
	RealType m, m1, m2

	RealType km, kp

	kp = (m1 + m2)/m
	kp = (1 - kp)*(1 + kp)
	if( kp .gt. 0 ) then
	  km = (m1 - m2)/m
	  km = (1 - km)*(1 + km)
	  Kallen = sqrt(km*kp)
	else
	  Kallen = 0
	endif
	end


************************************************************************
	
	RealType function KallenLR(m, m1, m2, cL, cR)
	implicit none
	RealType m, m1, m2
	ComplexType cL, cR

	RealType km, kp

	kp = (m1 + m2)/m
	kp = (1 - kp)*(1 + kp)
	if( kp .gt. 0 ) then
	  km = (m1 - m2)/m
	  km = (1 - km)*(1 + km)
	  KallenLR = .25D0*sqrt(km*kp)*
     &      (km*Sq(cL - cR) + kp*Sq(cL + cR))
	else
	  KallenLR = 0
	endif
	end


************************************************************************

	RealType function HXVWidth(mh, mx, mv, gammav,
     &    delta1, delta2, isvector)
	implicit none
	RealType mh, mx, mv, gammav, delta1, delta2
	integer isvector

#include "FH.h"

	RealType m, m1, g, c
	common /HXVVars/ m, m1, g, c

	RealType Kallen
	external Kallen, HXVIntegrand

	RealType lambda, d1, d2, w1, w2

	d1 = mx + mv - delta1
	d2 = mx + mv + delta2

	if( mh .gt. d2 ) then
* real decay above threshold
	  lambda = Kallen(mh, mx, mv)
	  HXVWidth = lambda*(lambda**2 + isvector*12*(mv/mh)**4)
	else if( mh .gt. d1 ) then
* interpolate inside [threshold - delta1, threshold + delta2]
	  lambda = Kallen(d2, mx, mv)
	  w2 = lambda*(lambda**2 + isvector*12*(mv/d2)**4)
	  m = d1/mv
	  m1 = mx/mv
	  g = gammav/mv
	  c = isvector*12
	  call Gauss(1, w1, HXVIntegrand, 0D0, m - m1, 16)
	  w1 = (isvector + 1)*w1
	  HXVWidth = w1 + (w2 - w1)*(mh - d1)/(delta1 + delta2)
	else if( mh .gt. mx + 1 ) then
* effects from finite width of gauge boson
	  m = mh/mv
	  m1 = mx/mv
	  g = gammav/mv
	  c = isvector*12
	  call Gauss(1, w1, HXVIntegrand, 0D0, m - m1, 16)
	  HXVWidth = (isvector + 1)*w1
	else
	  HXVWidth = 0
	endif

	HXVWidth = GF/(8*sqrt2*pi)*mh**3*HXVWidth
	end


************************************************************************

	subroutine HXVIntegrand(res, x)
	implicit none
	RealType res, x

#include "const.h"

	RealType Kallen
	external Kallen

	RealType m, m1, g, c
	common /HXVVars/ m, m1, g, c

	RealType lambda

	lambda = Kallen(m, m1, x)

	res = 2/pi*x*g*lambda*(lambda**2 + c*(x/m**2)**2)/
     &    (((x - 1)*(x + 1))**2 + g**2)
	end


************************************************************************

	subroutine SetDecayScale(fail, mh, mtr, sort)
	implicit none
	integer fail, sort
	RealType mh, mtr

#include "Decay.h"

	RealType qrun
	integer nf

	RealType Run
	external Run

	call AlphaS(AlfasMH, nf, mh)

	qrun = Run(2D0, mh)

	Mf(tH,1) = MU*qrun
	Mf(tH,2) = MC*Run(MC, mh)
	Mf(tH,3) = mtr
	if( Mf(tH,3) .lt. 0 ) Mf(tH,3) = Mf(tT,3)*Run(MT, mh)
	if( debuglevel .ge. 5 ) then
	  DCOUP "SetDecayScale Q =", mh ENDL
	  DCOUP "Alfas(Q) =", AlfasMH	ENDL
	  DCOUP "Mf(tH,1) =", Mf(tH,1)	ENDL
	  DCOUP "  (tH,2) =", Mf(tH,2)	ENDL
	  DCOUP "  (tH,3) =", Mf(tH,3)	ENDL
	endif
	call SfUpdate(fail, X3(tH,3,3), sort, 1)

	Mf(bH,1) = MD*qrun
	Mf(bH,2) = MS*qrun
	Mf(bH,3) = MB*Run(MB, mh)
	if( debuglevel .ge. 5 ) then
	  DCOUP "Mf(bH,1) =", Mf(bH,1)	ENDL
	  DCOUP "  (bH,2) =", Mf(bH,2)	ENDL
	  DCOUP "  (bH,3) =", Mf(bH,3)	ENDL
	endif
	call SfUpdate(fail, X3(bH,4,4), sort, 1)

	Mf(bHR,1) = Mf(bH,1)/abs(1 + Deltaf(bTR,1))
	Mf(bHR,2) = Mf(bH,2)/abs(1 + Deltaf(bTR,2))
	Mf(bHR,3) = Mf(bH,3)/abs(1 + Deltaf(bTR,3))
	call SfUpdate(fail, X3(bHR,4,4), sort, 1)

	hggSq(1) = sqrt(1 + AlfasMH/pi*(95/4D0 - 7/6D0*nf))
	hggSq(2) = hggSq(1)
	hggSq(3) = sqrt(1 + AlfasMH/pi*(97/4D0 - 7/6D0*nf))
	end


************************************************************************

	subroutine Kin1to2(fac, h, m2_, m3_, symm)
	implicit none
	integer h
	RealType fac, m2_, m3_, symm

#include "Decay.h"

	RealType sqrtS, Pout

	sqrtS = MHiggs(h)
	if( sqrtS .le. m2_ + m3_ ) then
	  fac = 0
	  return
	endif

	Pout = .5D0*(sqrtS - (m2_ - m3_)*(m2_ + m3_)/sqrtS)
	Pout = (Pout - m3_)*(Pout + m3_)
	Pout = sqrt(Pout)

	call VecSet(1, sqrtS, 0D0, 0D0, 0D0, 1D0)
	call VecSet(2, m2_, Pout, 1D0, 0D0, 0D0)
	call VecSet(3, m3_, Pout, -1D0, 0D0, 0D0)

* 4 pi from the integration, symmetry factors, and flux
	fac = 4*pi/(2*(2*pi)**2*4)*Pout/sqrtS**2*symm

	hno = h

	if( debuglevel .lt. 5 ) return

	DCOUP "Kin1to2: h =", h	ENDL
	DCOUP "  sqrtS =", sqrtS ENDL
	DCOUP "  m2    =", m2_	ENDL
	DCOUP "  m3    =", m3_	ENDL
	end


************************************************************************

	subroutine HffDecays(h, t, g1, g2, sq, sqSM)
	implicit none
	integer h, t, g1, g2
	external sq, sqSM

#include "Decay.h"

	RealType fac, sqmeSM(3), sqme(3)
	RealType m, b2, b, lb, d, ld, irH, ir
	character*8 id
	integer ty

	RealType Li2
	external Li2

	ty = t + ifB(t, bH - 4)
	call Kin1to2(fac, h, Mf(ty,g1), Mf(ty,g2), 1D0)
	if( fac .eq. 0 ) return

	gno1 = g1
	gno2 = g2

	ir = 0
	if( g1 .eq. g2 ) then
* IR factor from Dabelstein, hep-ph/9503443, Eq. (4.4)
	  m = (m2 + m3)/m1
	  b2 = (1 - m)*(1 + m)
	  b = sqrt(b2)
	  lb = log(b)
	  d = (1 - b)/(1 + b)
	  ld = log(d)
	  irH = (1 + b2)/b*
     &        (4*Li2(d) + 2*Li2(-d) + (3*log(2/(1 + b)) + 2*lb)*ld) +
     &      6*log(.5D0*m) - 4*lb -
     &      (3/16D0 + b2*(34/16D0 - 13/16D0*b2))/(b*b2)*ld +
     &      21/8D0 - 3/8D0/b2
* we ignore the extra term for the pseudoscalar, Eq. (4.3),
* as the contribution is tiny and is straightforward to apply
* only in the 2x2 mixing case: irH = irH + m**4
	  ir = Alfa0/pi*Qf(t)**2*irH
	  if( t .gt. 2 )
     &      ir = ir + AlfasMH/pi*C_F*(irH + 3*log(m1/m2) + 2)
* 2 (= 8/4) is due to the MSbar shift, K Williams PhD Eq. (7.5)
	  fac = fac*(1 + ir)
	endif

	call sqSM(sqmeSM, O'155')
	sqmeSM(3) = sqmeSM(1) + sqmeSM(2)
	GammaSM(H0FF(h,t,g1,g2)) = fac*sqmeSM(3)
c	GammaSM(H0FF(h,t,g1,g2)) = max(GammaSM(H0FF(h,t,g1,g2)), 0D0)
	GammaSMTot(h) = GammaSMTot(h) + GammaSM(H0FF(h,t,g1,g2))

	call sq(sqme, O'155')
	sqme(3) = sqme(1) + sqme(2)
	Gamma(H0FF(h,t,g1,g2)) = fac*sqme(3)
c	Gamma(H0FF(h,t,g1,g2)) = max(Gamma(H0FF(h,t,g1,g2)), 0D0)
	GammaTot(h) = GammaTot(h) + Gamma(H0FF(h,t,g1,g2))

	if( debuglevel .lt. 3 ) return

	id = "Hff"//Digit(h)//Digit(t)//Digit(g1)//Digit(g2)//" "

	DCOUP id, "GammaSM(tree) =", fac*sqmeSM(1)	ENDL
	DCOUP id, "GammaSM(loop) =", fac*sqmeSM(2)	ENDL
	DCOUP id, "GammaSM       =", fac*sqmeSM(3)	ENDL
	DCOUP id, "Gamma(tree)   =", fac*sqme(1)	ENDL
	DCOUP id, "Gamma(loop)   =", fac*sqme(2)	ENDL
	DCOUP id, "Gamma         =", fac*sqme(3)	ENDL
	DCOUP id, "  Mh        =", m1			ENDL
	DCOUP id, "  Mf1       =", m2			ENDL
	DCOUP id, "  Mf2       =", m3			ENDL
	DCOUP id, "  prefactor =", fac			ENDL
	DCOUP id, "  IR factor =", ir			ENDL
	end

************************************************************************

	subroutine HhhDecays(h, h1, h2)
	implicit none
	integer h, h1, h2

#include "Decay.h"

	RealType fac, sqme(3)
	character*7 id

	call Kin1to2(fac, h, MHiggs(h1), MHiggs(h2),
     &    .5D0*(2 - Delta(h1, h2)))
	if( fac .eq. 0 ) return

	hno1 = h1
	hno2 = h2
	call hhh_SquaredME(sqme)

	sqme(3) = sqme(1) + sqme(2)
	Gamma(H0HH(h,h1,h2)) = fac*sqme(3)
	GammaTot(h) = GammaTot(h) + Gamma(H0HH(h,h1,h2))

	if( debuglevel .lt. 3 ) return

	id = "Hhh"//Digit(h)//Digit(h1)//Digit(h2)//" "

	DCOUP id, "Gamma(tree) =", fac*sqme(1)	ENDL
	DCOUP id, "Gamma(loop) =", fac*sqme(2)	ENDL
	DCOUP id, "Gamma       =", fac*sqme(3)	ENDL
	DCOUP id, "  Mh        =", m1		ENDL
	DCOUP id, "  Mh1       =", m2		ENDL
	DCOUP id, "  Mh2       =", m3		ENDL
	DCOUP id, "  prefactor =", fac		ENDL
	end

