* CalcSelf.F
* Calculation of the renormalized Higgs self-energies
* this file is part of FeynHiggs
* last modified 12 Nov 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine CalcSelf(k2, key, dkey)
	implicit none
	RealType k2
	integer key, dkey

#include "FH.h"
#include "OneLoop.h"

#define SelfEnergy(i,f,bX) if( ibits(key,i-1,1) .eq. 1 ) \
seR(i) = f(p2,bX)+dZ(i)*(p2-Msq(i))-dMsq(i)+se2R(i)+seAdd(i)

#define DSelfEnergy(i,f,bX) if( ibits(dkey,i-1,1) .eq. 1 ) \
dseR(i) = f(p2,bX) + dZ(i)

	RealType p2
	logical onshell
	integer i

* p2approx decides at which p^2 the self-energies are to be evaluated
	p2 = k2
	onshell = p2approx .eq. 2 .or.
     &    ((p2approx .eq. 0 .or. p2approx .eq. 3) .and. p2 .eq. -1)
	if( p2approx .eq. 1 .or. p2 .eq. -1 ) p2 = 0

	do i = 1, semax
	  seR(i) = 0
	  dseR(i) = 0
	enddo

* 1) the CP-even bosons h & H + the mixing h-H

	if( iand(key, cpeven) .ne. 0 .or.
     &      iand(dkey, cpeven) .ne. 0 ) then
	  if( mssmpart .eq. 0 ) then
	    call CPevenSEs_MTapprox
	  else
	    if( onshell ) p2 = Msq(h0h0)
	    SelfEnergy(h0h0, se_h0h0, bM0)
	    DSelfEnergy(h0h0, dse_h0h0, bM0)

	    if( onshell ) p2 = Msq(HHHH)
	    SelfEnergy(HHHH, se_HHHH, bM0)
	    DSelfEnergy(HHHH, dse_HHHH, bM0)

	    if( onshell ) p2 = Msq(h0HH)
	    SelfEnergy(h0HH, se_h0HH, bM0)
	    DSelfEnergy(h0HH, dse_h0HH, bM0)
	  endif
	endif

* 2) the CP-odd boson A

	if( onshell ) p2 = Msq(A0A0)
	SelfEnergy(A0A0, se_A0A0, bM)
	DSelfEnergy(A0A0, dse_A0A0, bM)

* 3) the CP-violating h-A & H-A mixing

	if( onshell ) p2 = Msq(h0A0)
	SelfEnergy(h0A0, se_h0A0, bM)
	DSelfEnergy(h0A0, dse_h0A0, bM)

	if( onshell ) p2 = Msq(HHA0)
	SelfEnergy(HHA0, se_HHA0, bM)
	DSelfEnergy(HHA0, dse_HHA0, bM)

* 4) the charged boson H+

#if 1
	if( iand(key, Key(HmHp)) .ne. 0 .or.
     &      iand(dkey, Key(HmHp)) .ne. 0 ) then
	  if( mssmpart .eq. 0 ) then
	    call HpSEs_MTapprox
	  else
	    if( onshell ) p2 = Msq(HmHp)
	    SelfEnergy(HmHp, se_HmHp, bM)
	    DSelfEnergy(HmHp, dse_HmHp, bM)
	  endif
	endif
#else
	if( onshell ) p2 = Msq(HmHp)
	SelfEnergy(HmHp, se_HmHp, bM)
	DSelfEnergy(HmHp, dse_HmHp, bM)
#endif

* 5) the Goldstones

	if( onshell ) p2 = Msq(G0G0)
	SelfEnergy(G0G0, se_G0G0, bM)
	DSelfEnergy(G0G0, dse_G0G0, bM)

	if( onshell ) p2 = Msq(h0G0)
	SelfEnergy(h0G0, se_h0G0, bM)
	DSelfEnergy(h0G0, dse_h0G0, bM)

	if( onshell ) p2 = Msq(HHG0)
	SelfEnergy(HHG0, se_HHG0, bM)
	DSelfEnergy(HHG0, dse_HHG0, bM)

	if( onshell ) p2 = Msq(A0G0)
	SelfEnergy(A0G0, se_A0G0, bM)
	DSelfEnergy(A0G0, dse_A0G0, bM)

	if( onshell ) p2 = Msq(GmGp)
	SelfEnergy(GmGp, se_GmGp, bM)
	DSelfEnergy(GmGp, dse_GmGp, bM)

	if( onshell ) p2 = Msq(HmGp)
	SelfEnergy(HmGp, se_HmGp, bM)
	DSelfEnergy(HmGp, dse_HmGp, bM)


	if( p2approx .eq. 3 ) then
	  do i = 1, semax
	    seR(i) = Re(seR(i))
	    dseR(i) = Re(dseR(i))
	  enddo
	endif

	if( debuglevel .lt. 4 ) return

	if( onshell ) then
	  DHIGGS "k2 = onshell" ENDL
	else
	  DHIGGS "k2 =", p2 ENDL
	endif
	DHIGGS "seR(h0h0) =", seR(h0h0) ENDL
	DHIGGS "seR(HHHH) =", seR(HHHH) ENDL
	DHIGGS "seR(A0A0) =", seR(A0A0) ENDL
	DHIGGS "seR(HmHp) =", seR(HmHp) ENDL
	DHIGGS "seR(h0HH) =", seR(h0HH) ENDL
	DHIGGS "seR(h0A0) =", seR(h0A0) ENDL
	DHIGGS "seR(HHA0) =", seR(HHA0) ENDL

	if( debuglevel .lt. 5 ) return

	DHIGGS "seR(G0G0) =", seR(G0G0) ENDL
	DHIGGS "seR(h0G0) =", seR(h0G0) ENDL
	DHIGGS "seR(HHG0) =", seR(HHG0) ENDL
	DHIGGS "seR(A0G0) =", seR(A0G0) ENDL
	DHIGGS "seR(GmGp) =", seR(GmGp) ENDL
	DHIGGS "seR(HmGp) =", seR(HmGp) ENDL
	end


************************************************************************
* CPevenSEs_MTapprox: the one-loop CP-even renormalized self-energies
* in the MT^4-approximation as specified in Sven Heinemeyer's PhD
* thesis, p. 74.

	subroutine CPevenSEs_MTapprox
	implicit none

#include "FH.h"

	RealType mlog, x, x2, c
	ComplexType sigma_t, lam_t, omega_t

	if( MSf2(3,tM1,3) .lt. 1D-12 ) then
	  x = Re(Xf(3,3))/MSf2(1,tM1,3)
	  x2 = -1/12D0*x**2
	else
	  mlog = log(MSf2(2,tM1,3)/MSf2(1,tM1,3))
	  x = Re(Xf(3,3))/MSf2(3,tM1,3)
	  x2 = x**2*(1 - .5D0*MSf2(4,tM1,3)/MSf2(3,tM1,3)*mlog)
	  x = x*mlog
	endif

	c = 3*GF*Mf2(tM1,3)**2/(sqrt2*pi**2*SB2)

	sigma_t = c*MUE**2*x2

	lam_t = c*MUE*(-.5D0*x - Af(3,3)*x2)

	omega_t = c*(log(MSf(1,tM1,3)*MSf(2,tM1,3)/Mf2(tM1,3)) +
     &    Af(3,3)*(x + Af(3,3)*x2))

	seR(h0h0) = -(SA2*sigma_t + CA2*omega_t - 2*SA*CA*lam_t)
	seR(HHHH) = -(CA2*sigma_t + SA2*omega_t + 2*SA*CA*lam_t)
	seR(h0HH) = -((CA2 - SA2)*lam_t + SA*CA*(omega_t - sigma_t))
	end


************************************************************************
* HpSEs_MTapprox: the one-loop renormalized charged-Higgs self-energy
* in the MT^4-approximation

#define __SUBROUTINE__ "HpSEs_MTapprox"

	subroutine HpSEs_MTapprox
	implicit none

#include "FH.h"

	RealType DSStop2(2), MStop2(4), MStop(2)
	ComplexType UStop(2,2)
	RealType d, MSbot12, mlog1, dmlog
	integer fail

	RealType ctt, stt
	equivalence (UStop(1,1), ctt)
	equivalence (UStop(1,2), stt)

* use squark masses w/o D-terms
	DSStop2(1) = DSS2(1,3,3) - (DSf(1,3) - DSf(2,3))
	DSStop2(2) = DSS2(2,3,3) - DSf(2,3)
	if( DSStop2(1)**2 + Sq(Xf(3,3)) .lt. .1D0 )
     &    DSStop2(1) = -.2D0*sqrt(DSStop2(2))
	fail = 0
	call EigenSf(fail, 1, MStop, MStop2, UStop,
     &    DSStop2, Mf(tM1,3), Xf(3,3), 0)
	if( fail .ne. 0 ) Warning("negative stop mass squares")

	d = Mf2(tM1,3) - stt**2*MStop2(3)
	MSbot12 = MStop2(1) - d

	dmlog = log(MStop2(1))
	mlog1 = MStop2(1)*(dmlog - 1)
	if( abs(MStop2(3)) .gt. 1D-12 )
     &    dmlog = (MStop2(2)*(log(MStop2(2)) - 1) - mlog1)/MStop2(3)

	seR(HmHp) = 3/(8.D0*pi)*Mf2(tM1,3)*
     &    AlfaGF/(MW2*SW2*TB2*(Mf2(tM1,3) + ctt**2*MStop2(3)))*
     &    (Mf(tM1,3)*(Af(3,3) + MUE*TB) + ctt*stt*MStop2(3))**2*
     &    (dmlog + (MSbot12*(log(MSbot12) - 1) - mlog1)/d)
	end


************************************************************************
* CalcCTs: calculate the relevant counterterms in the Higgs sector.

	subroutine CalcCTs
	implicit none

#include "FH.h"
#include "OneLoop.h"

	RealType d1, d2, c
	RealType dZH(2,4:5), dTB(4:5)
	RealType dMZ2, dMA2, dZA0G0, dTHH, dTh0, dTA0

	Msq(h0h0) = Mh02
	Msq(HHHH) = MHH2
	Msq(A0A0) = MA02
	Msq(HmHp) = MHp2
	Msq(h0HH) = .5D0*(Msq(h0h0) + Msq(HHHH))
	Msq(h0A0) = .5D0*(Msq(h0h0) + Msq(A0A0))
	Msq(HHA0) = .5D0*(Msq(HHHH) + Msq(A0A0))
	Msq(G0G0) = 0
	Msq(h0G0) = .5D0*(Msq(h0h0) + Msq(G0G0))
	Msq(HHG0) = .5D0*(Msq(HHHH) + Msq(G0G0))
	Msq(A0G0) = .5D0*(Msq(A0A0) + Msq(G0G0))
	Msq(GmGp) = 0
	Msq(HmGp) = .5D0*(Msq(HmHp) + Msq(GmGp))

* choose the field renormalization:
	if( fieldren .eq. 4 ) then
* on-shell field renormalization IV goes like this...
	  dZ(h0h0) = -Re(dse_h0h0(0D0, bM0))
	  dZ(HHHH) = -Re(dse_HHHH(0D0, bM0))
	  dZ(h0HH) = -Re(dse_h0HH(0D0, bM0))

	  dZ(A0A0) = -Re(dse_A0A0(0D0, bM))
	  dZ(HmHp) = -Re(dse_HmHp(0D0, bM))
	  dZ(h0A0) = -Re(dse_h0A0(0D0, bM))
	  dZ(HHA0) = -Re(dse_HHA0(0D0, bM))
	  dZ(G0G0) = -Re(dse_G0G0(0D0, bM))
	  dZ(h0G0) = -Re(dse_h0G0(0D0, bM))
	  dZ(HHG0) = -Re(dse_HHG0(0D0, bM))
	  dZ(A0G0) = -Re(dse_A0G0(0D0, bM))
	  dZA0G0 = -Re(dse_A0G0(0D0, bM0))
	  dZ(GmGp) = -Re(dse_GmGp(0D0, bM))
	  dZ(HmGp) = -Re(dse_HmGp(0D0, bM))
	else
* ... all others use these relations:
	  if( fieldren .eq. 0 ) then
	    dZH(1,4) = 0
	    dZH(2,4) = 0
	    dZH(1,5) = 0
	    dZH(2,5) = 0
	  else if( fieldren .eq. 1 ) then
	    d1 = Re(dse_A0A0(MA02, bM))
	    d2 = Re(se_A0Z(MA02, bM))
	    dZH(1,4) = -d1 - d2/(MZ*TB)
	    dZH(2,4) = -d1 + d2*TB/MZ
	    d1 = Re(dse_A0A0(MA02, bM0))
	    d2 = Re(se_A0Z(MA02, bM0))
	    dZH(1,5) = -d1 - d2/(MZ*TB)
	    dZH(2,5) = -d1 + d2*TB/MZ
	  else if( fieldren .eq. 2 ) then
	    d1 = Re(dse_G0G0(MZ2, bM))
	    d2 = Re(dse_A0A0(MA02, bM))
	    dZH(1,4) = -( CB2*d1 - SB2*d2)/C2B
	    dZH(2,4) = -(-SB2*d1 + CB2*d2)/C2B
	    d1 = Re(dse_G0G0(MZ2, bM0))
	    d2 = Re(dse_A0A0(MA02, bM0))
	    dZH(1,5) = -( CB2*d1 - SB2*d2)/C2B
	    dZH(2,5) = -(-SB2*d1 + CB2*d2)/C2B
	  else
	    d1 = Re(dse_HHHH(MHH2, bM))
	    d2 = Re(dse_h0h0(Mh02, bM))
	    dZH(1,4) = -( CA2*d1 - SA2*d2)/C2A
	    dZH(2,4) = -(-SA2*d1 + CA2*d2)/C2A
	    d1 = Re(dse_HHHH(MHH2, bM0))
	    d2 = Re(dse_h0h0(Mh02, bM0))
	    dZH(1,5) = -( CA2*d1 - SA2*d2)/C2A
	    dZH(2,5) = -(-SA2*d1 + CA2*d2)/C2A
	  endif

	  dZ(h0h0) = SA2*dZH(1,5) + CA2*dZH(2,5)
	  dZ(HHHH) = CA2*dZH(1,5) + SA2*dZH(2,5)
	  dZ(h0HH) = SA*CA*(dZH(2,5) - dZH(1,5))

	  dZ(A0A0) = SB2*dZH(1,4) + CB2*dZH(2,4)
	  dZ(HmHp) = dZ(A0A0)
	  dZ(h0A0) = 0
	  dZ(HHA0) = 0
	  dZ(G0G0) = CB2*dZH(1,4) + SB2*dZH(2,4)
	  dZ(GmGp) = dZ(G0G0)
	  dZ(h0G0) = 0
	  dZ(HHG0) = 0
	  dZ(A0G0) = SB*CB*(dZH(2,4) - dZH(1,4))
	  dZA0G0 = SB*CB*(dZH(2,5) - dZH(1,5))
	  dZ(HmGp) = dZ(A0G0)
	endif

* renormalization of TB:
	if( tanbren .eq. 0 ) then
	  dTB(4) = 0
	  dTB(5) = 0
	else if( tanbren .eq. 1 ) then
	  dTB(4) = .5D0*(dZH(2,4) - dZH(1,4))
	  dTB(5) = .5D0*(dZH(2,5) - dZH(1,5))
	else
	  dTB(4) = (Re(se_A0Z(MA02, bM)) + .5D0*dZ(A0G0)*MZ)/(SB*CB*MZ)
	  dTB(5) = (Re(se_A0Z(MA02, bM0)) + .5D0*dZA0G0*MZ)/(SB*CB*MZ)
	endif

* depending on `inputmass', either MA0 or MHp is defined on-shell:
	if( inputmass .eq. 3 ) then
	  dMsq(A0A0) = Re(se_A0A0(MA02, bM))
	  dMsq(HmHp) = dMsq(A0A0) + Re(se_WW(MW2, bM))
	  dMA2 = Re(se_A0A0(MA02, bM0))
	else
	  dMsq(HmHp) = Re(se_HmHp(MHp2, bM))
	  dMsq(A0A0) = dMsq(HmHp) - Re(se_WW(MW2, bM))
	  dMA2 = Re(se_HmHp(MHp2, bM0)) - Re(se_WW(MW2, bM0))
	endif

	c = ELGF/(2*MZ*SW*CW)

* and finally the Higgs mass counterterms:
	dMZ2 = Re(se_ZZ(MZ2, bM0))
	dTh0 = -Re(td_h0(bM0))
	dTHH = -Re(td_HH(bM0))
	dTA0 = -Re(td_A0(bM0))

	dMsq(h0h0) = dMA2*CBA2 + dMZ2*SAB**2 +
     &    c*SBA*(dTHH*CBA*SBA - dTh0*(1 + CBA2)) -
     &    dTB(5)*S2B*(MA02*SBA*CBA - MZ2*SAB*CAB)

	dMsq(HHHH) = dMA2*SBA2 + dMZ2*CAB**2 -
     &    c*CBA*(dTHH*(1 + SBA2) - dTh0*SBA*CBA) +
     &    dTB(5)*S2B*(MA02*SBA*CBA - MZ2*SAB*CAB)

	dMsq(h0HH) = -dMA2*SBA*CBA - dMZ2*SAB*CAB -
     &    c*(dTHH*SBA**3 + dTh0*CBA**3) -
     &    dTB(5)*SB*CB*(MA02*(CBA - SBA)*(CBA + SBA) +
     &      MZ2*(CAB - SAB)*(CAB + SAB))

	dTh0 = -Re(td_h0(bM))
	dTHH = -Re(td_HH(bM))
	dTA0 = -Re(td_A0(bM))

	dMsq(h0A0) = -c*dTA0*SBA
	dMsq(HHG0) = dMsq(h0A0)
	dMsq(h0G0) = c*dTA0*CBA
	dMsq(HHA0) = -dMsq(h0G0)

	dMsq(G0G0) = -c*(dTHH*CBA + dTh0*SBA)
	dMsq(GmGp) = dMsq(G0G0)
	dMsq(A0G0) = c*(dTHH*SBA - dTh0*CBA) - dTB(4)*SB*CB*MA02
	dMsq(HmGp) = c*(dTHH*SBA - dTh0*CBA - cI*dTA0) -
     &    dTB(4)*SB*CB*MHp2

	if( debuglevel .lt. 2 ) return

	DSELF "========================"	ENDL
	DSELF "counter terms by CalcCTs"	ENDL
	DSELF "------------------------"	ENDL
	DSELF "dMsq(h0h0) =", dMsq(h0h0)	ENDL
	DSELF "dZ(h0h0)   =", dZ(h0h0)		ENDL
	DSELF "dMsq(HHHH) =", dMsq(HHHH)	ENDL
	DSELF "dZ(HHHH)   =", dZ(HHHH)		ENDL
	DSELF "dMsq(A0A0) =", dMsq(A0A0)	ENDL
	DSELF "dZ(A0A0)   =", dZ(A0A0)		ENDL
	DSELF "dMsq(HmHp) =", dMsq(HmHp)	ENDL
	DSELF "dZ(HmHp)   =", dZ(HmHp)		ENDL
	DSELF "dMsq(h0HH) =", dMsq(h0HH)	ENDL
	DSELF "dZ(h0HH)   =", dZ(h0HH)		ENDL
	DSELF "dMsq(h0A0) =", dMsq(h0A0)	ENDL
	DSELF "dZ(h0A0)   =", dZ(h0A0)		ENDL
	DSELF "dMsq(HHA0) =", dMsq(HHA0)	ENDL
	DSELF "dZ(HHA0)   =", dZ(HHA0)		ENDL
	DSELF "dMsq(G0G0) =", dMsq(G0G0)	ENDL
	DSELF "dZ(G0G0)   =", dZ(G0G0)		ENDL
	DSELF "dMsq(h0G0) =", dMsq(h0G0)	ENDL
	DSELF "dZ(h0G0)   =", dZ(h0G0)		ENDL
	DSELF "dMsq(HHG0) =", dMsq(HHG0)	ENDL
	DSELF "dZ(HHG0)   =", dZ(HHG0)		ENDL
	DSELF "dMsq(A0G0) =", dMsq(A0G0)	ENDL
	DSELF "dZ(A0G0)   =", dZ(A0G0)		ENDL
	DSELF "dMsq(GmGp) =", dMsq(GmGp)	ENDL
	DSELF "dZ(GmGp)   =", dZ(GmGp)		ENDL
	DSELF "dMsq(HmGp) =", dMsq(HmGp)	ENDL
	DSELF "dZ(HmGp)   =", dZ(HmGp)		ENDL
	end

