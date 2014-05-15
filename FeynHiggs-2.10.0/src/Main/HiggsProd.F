* HiggsProd.F
* approximate calculation of the Higgs production cross-sections
* neutral-Higgs cross-section data are from Fabio Maltoni:
*   http://maltoni.home.cern.ch/maltoni/TeV4LHC
* charged-Higgs cross-section data are from Tilman Plehn
*   http://www.ph.ed.ac.uk/~tplehn/charged_higgs
* this file is part of FeynHiggs
* last modified 18 Dec 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"
#include "FHCouplings.h"

#define cTopNLO ck(1)
#define cTopNNLO ck(2)
#define cBotReNLO ck(3)
#define cBotImNLO ck(4)
#define kTopNLO ck(5)
#define kTopNNLO ck(6)
#define kBotNLO ck(7)
#define kTopBotNLO ck(8)
#define nck 8


	subroutine FHHiggsProd(error, sqrts, prodxs)
	implicit none
	integer error
	RealType sqrts, prodxs(nprodxs)

#include "FH.h"
#define __SUBROUTINE__ "FHHiggsProd"

	integer i, h
	character*8 s

	logical blurb
	character*2 higgs(4)
	save blurb, higgs
	data blurb /.TRUE./
	data higgs /"h0", "HH", "A0", "Hp"/

	do i = 1, nprodxs
	  prodxs(i) = -999
	enddo

	call CheckCouplings(error, iand(ieor(coup_valid, valid), 1))
	if( error .gt. 0 ) return

	if( blurb ) then
	  print *, "FHHiggsProd contains code by"
	  print *, "- SM XS for VBF, WH, ZH, ttH are taken"
	  print *, "  from the LHC Higgs Cross Section Working Group,"
	  print *, "  https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CrossSections"
	  print *, "  see there for an extensive list of references"
	  print *, "- SM bbH XS: Harlander et al. hep-ph/0304035"
	  print *, "- SM ggH XS: Grazzini et al."
	  print *, "  http://theory.fi.infn.it/grazzini/hcalculators.html"
	  print *, "  see there for an extensive list of references"
	  print *, "- 2HDM charged Higgs XS: Plehn et al."
	  print *, "- heavy charged Higgs XS:"
	  print *, "  Dittmaier et al. arXiv:0906.2648"
	  print *, "  Flechl et al. arXiv:1307.1347"
	  blurb = .FALSE.
	endif

	if( int(sqrts + .05D0) .eq. 2 ) then
	  call HiggsProdTev(error, prodxs)
	  s = "Tevatron"
	else if( sqrts .gt. 2.5D0 ) then
	  call HiggsProdLHC(error, sqrts, prodxs)
	  s = "LHC"
	else
	  Error(error, "Bad energy")
	  return
	endif

	if( debuglevel .lt. 2 ) return

	DPROD s, " production cross-sections (MSSM, SM):" ENDL
	DPROD "sqrt(s) =", sqrts ENDL

	do h = 1, 3
	  DPROD "b-b-"//higgs(h)//"    =", bbh(h), bbhSM(h)	ENDL
	  DPROD "btag-b-"//higgs(h)//" =", btagbh(h), btagbhSM(h) ENDL
	  DPROD "g-g-"//higgs(h)//",c  =", ggh(h), gghSM(h)	ENDL
	  DPROD "g-g-"//higgs(h)//",k  =", ggh2(h)		ENDL
	  DPROD "q-q-"//higgs(h)//"    =", qqh(h), qqhSM(h)	ENDL
	  DPROD "t-t-"//higgs(h)//"    =", tth(h), tthSM(h)	ENDL
	  DPROD "W-"//higgs(h)//"      =", Wh(h), WhSM(h)	ENDL
	  DPROD "Z-"//higgs(h)//"      =", Zh(h), ZhSM(h)	ENDL
	  DPROD "St-St-"//higgs(h)//"  =", StSth(h)		ENDL
	enddo
	DPROD "t-Hm      =", tHm ENDL
	end


************************************************************************

	subroutine HiggsProdTev(error, prodxs)
	implicit none
	integer error
	RealType prodxs(nprodxs)

#include "FH.h"

	integer h
	RealType m, sqrtm, logm, c, ck(nck)
	RealType gghSMNLO(3)

	do h = 1, 3
	  m = MHiggs(h)
	  if( m .gt. 80 ) then
	    sqrtm = sqrt(m)
	    logm = log(m)

#include "NHiggsProdFits-Tev.f"
#include "gghcTev.f"
#include "gghkTev.f"

	    c = GammaSM(H0FF(h,4,3,3))
	    if( c .eq. 0 ) then
	      bbh(h) = 0
	      bbhSM(h) = 0
	      btagbh(h) = 0
	      btagbhSM(h) = 0
	    else
	      c = Gamma(H0FF(h,4,3,3))/c
	      bbh(h) = c*bbhSM(h)
	      btagbh(h) = c*btagbhSM(h)
	    endif

	    c = ratios(H0FF(h,3,3,3))
	    tth(h) = c*tthSM(h)

	    if( m .ge. 90 .and. m .le. 300 )
     &        call gghxs(error, h, ggh(h), ggh2(h), gghSMNLO(h), ck)

	    c = ratios(H0VV(h,4))
	    qqh(h) = c*qqhSM(h)
	    Wh(h) = c*WhSM(h)

	    c = ratios(H0VV(h,3))
	    Zh(h) = c*ZhSM(h)
	  endif
	enddo
	end


************************************************************************

	subroutine HiggsProdLHC(error, sqrts, prodxs)
	implicit none
	integer error
	RealType sqrts, prodxs(nprodxs)

#include "FH.h"

	integer h
	RealType logsqrts, m, sqrtm, logm, mst1
	RealType c, TBeff, ck(nck)
	RealType gghSMNLO(3)

	logsqrts = log(sqrts)

	do h = 1, 3
	  m = MHiggs(h)
	  if( m .gt. 80 ) then
	    sqrtm = sqrt(m)
	    logm = log(m)

#include "NHiggsProdFits-LHC.f"
#include "gghcLHC.f"
#include "gghkLHC.f"

	    if( m .gt. 90 ) then
	      call gghxs(error, h, ggh(h), ggh2(h), gghSMNLO(h), ck)
	      if( m .gt. 400 ) ggh2(h) = -999
	    endif

	    c = GammaSM(H0FF(h,4,3,3))
	    if( c .eq. 0 ) then
	      bbh(h) = 0
	      bbhSM(h) = 0
	    else
	      c = Gamma(H0FF(h,4,3,3))/c
	      bbh(h) = c*bbhSM(h)
	    endif

	    if( abs(sqrts - 7) .lt. 1D-5 ) then

#include "NHiggsProdFits-LHC7.f"

	    else if( abs(sqrts - 8) .lt. 1D-5 ) then

#include "NHiggsProdFits-LHC8.f"

	    else if( abs(sqrts - 14) .lt. 1D-5 ) then

	      mst1 = MSf(1,3,3)

#include "NHiggsProdFits-LHC14.f"

	      if( c .eq. 0 ) then
	        btagbh(h) = 0
	        btagbhSM(h) = 0
	      else
	        btagbh(h) = c*btagbhSM(h)
	      endif

	      StSth(h) = StSth(h)*chSt1St1(H0SfSf(h,1,1,3,3))

	    else
	      goto 99
	    endif

	    c = ratios(H0FF(h,3,3,3))
	    tth(h) = c*tthSM(h)

	    c = ratios(H0VV(h,4))
	    qqh(h) = c*qqhSM(h)
	    Wh(h) = c*WhSM(h)

	    c = ratios(H0VV(h,3))
	    Zh(h) = c*ZhSM(h)

99	  continue
	  endif
	enddo

	if( MHiggs(4) .gt. 80 ) then
	  m = MHiggs(4)
	  sqrtm = sqrt(m)
	  TBeff = TB/sqrt(abs(1 + Deltaf(bTR,3)))

	  if( abs(sqrts - 14) .lt. 1D-5 ) then

#include "CHiggsProdFits-LHC14.f"

	    tHm = tHm/abs(1 + Deltaf(bTR,3))
	  endif

	  if( abs(sqrts - 8) .lt. 1D-5 .and. TBeff .ge. 30 ) then

#include "CHiggsProdFits-LHC8.f"

	  endif
	endif
	end


************************************************************************

	subroutine SetKFactors(k0, kt, kbre, kbim)
	implicit none
	RealType k0, kt, kbre, kbim

#include "FH.h"

	hggU(1,1) = k0
	hggU(2,1) = k0
	hggU(3,1) = k0
	hggU(1,2) = k0
	hggU(2,2) = k0
	hggU(3,2) = k0
	hggU(1,3) = kt
	hggU(2,3) = kt
	hggU(3,3) = kt

	hggDRe(1,1) = k0
	hggDRe(2,1) = k0
	hggDRe(3,1) = k0
	hggDRe(1,2) = k0
	hggDRe(2,2) = k0
	hggDRe(3,2) = k0
	hggDRe(1,3) = kbre
	hggDRe(2,3) = kbre
	hggDRe(3,3) = kbre

	hggDIm(1,1) = k0
	hggDIm(2,1) = k0
	hggDIm(3,1) = k0
	hggDIm(1,2) = k0
	hggDIm(2,2) = k0
	hggDIm(3,2) = k0
	hggDIm(1,3) = kbim
	hggDIm(2,3) = kbim
	hggDIm(3,3) = kbim
	end


************************************************************************

	subroutine gghxs(error, h, ggh_c, ggh_k, ggh_SM, ck)
	implicit none
	integer error, h
	RealType ggh_c, ggh_k, ggh_SM, ck(nck)

#include "FH.h"

	RealType c, csm, fac, sqme(2)

	if( debuglevel .gt. 4 ) then
	  DPROD "ggh factors"			ENDL
	  DPROD "cTopNLO    =", cTopNLO		ENDL
	  DPROD "cTopNNLO   =", cTopNNLO	ENDL
	  DPROD "cBotReNLO  =", cBotReNLO	ENDL
	  DPROD "cBotImNLO  =", cBotImNLO	ENDL
	  DPROD "kTopNLO    =", kTopNLO		ENDL
	  DPROD "kTopNNLO   =", kTopNNLO	ENDL
	  DPROD "kBotNLO    =", kBotNLO		ENDL
	  DPROD "kTopBotNLO =", kTopBotNLO	ENDL
	  DPROD "ggh_SM     =", ggh_SM		ENDL
	endif

* 172 hard-coded to match the Vicini data
	call SetDecayScale(error, MHiggs(h), 172D0, 0)
	call Kin1to2(fac, h, 0D0, 0D0, .5D0)

* ggh version 1 w/amplitude correction factors
	call SetKFactors(0D0, cTopNLO, cBotReNLO, cBotImNLO)
	call hggSM_SquaredME(sqme, O'155')
	csm = sqme(2)
	if( csm .eq. 0 ) then
	  ggh_c = 0
	else
	  call SetKFactors(1D0, cTopNNLO, cBotReNLO, cBotImNLO)
	  call hgg_SquaredME(sqme, O'155')
	  ggh_c = ggh_SM*sqme(2)/csm
	endif

* ggh version 2 w/genuine k-factors
	call SetKFactors(0D0, 1D0, 0D0, 0D0)
	call hggSM_SquaredME(sqme, O'155')
	csm = (kTopNLO - kTopBotNLO)*sqme(2)
	call hgg_SquaredME(sqme, O'155')
	c = (kTopNNLO - kTopBotNLO)*sqme(2)

	call SetKFactors(0D0, 0D0, 1D0, 1D0)
	call hggSM_SquaredME(sqme, O'155')
	csm = csm + (kBotNLO - kTopBotNLO)*sqme(2)
	call hgg_SquaredME(sqme, O'155')
	c = c + (kBotNLO - kTopBotNLO)*sqme(2)

	call SetKFactors(0D0, 1D0, 1D0, 1D0)
	call hggSM_SquaredME(sqme, O'155')
	csm = csm + kTopBotNLO*sqme(2)
	call hgg_SquaredME(sqme, O'155')
	c = c + (kTopBotNLO - 1)*sqme(2)

	if( csm .eq. 0 ) then
	  ggh_k = 0
	else
	  call SetKFactors(1D0, 1D0, 1D0, 1D0)
	  call hgg_SquaredME(sqme, O'155')
	  ggh_k = ggh_SM*(sqme(2) + c)/csm
	endif
	end
