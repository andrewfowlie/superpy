* TLShifts.F
* computes the two-loop corrections
* this file is part of FeynHiggs
* last modified 18 Dec 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"


************************************************************************
* calculate the two-loop self-energies
* Note: this subroutine is actually called from Sfermions.F since it
* also computes the Deltab compatible with the two-loop corrections.

	subroutine CalcRenSETL(error)
	implicit none
	integer error

#include "TLps.h"
#define __SUBROUTINE__ "CalcRenSETL"

	RealType dMTs, dMT, dPhi22
	RealType axax(3,5), dPhi(6,MaxSlots)
	integer i, fail

	RealType asat(3), asab(3), atat(3)
	RealType atab(3), atat_largeTB(3)
	equivalence (axax(1,1), asat)
	equivalence (axax(1,2), asab)
	equivalence (axax(1,3), atat)
	equivalence (axax(1,4), atab)
	equivalence (axax(1,5), atat_largeTB)

	logical blurb
	save blurb
	data blurb /.TRUE./

	call TLShifts(error)
	if( error .ne. 0 ) return
	call SfUpdate(error, X3(bTR0,bTR0,4), SORT_SF, 1)
	if( error .ne. 0 ) return

	Q = getmudim()

	do i = 1, semax
	  se2R(i) = 0
	enddo

	do i = 1, 5
	  axax(1,i) = 0
	  axax(2,i) = 0
	  axax(3,i) = 0
	enddo

	do i = 1, 6
	  dPhi(i,1) = 0
	enddo

	dMTs = 4/(3*pi)*Alfas2L*Mf(tM2,3)

	if( tlCplxApprox .eq. 1 ) goto 1

	if( blurb ) then
	  print *, "FHHiggsCorr contains code by:"
	  print *, "P. Slavich et al. (2-loop rMSSM Higgs self-energies)"
	  blurb = .FALSE.
	endif

	dMT = -1/(2*pi)*Alfat*Mf(tM2,3)
	if( tlCplxApprox .eq. 0 ) dMT = dMT + dMTs

	do vs = 1, ibset(0, nvars)
	  if( debuglevel .ge. 3 ) DSELF "point ", vs ENDL

	  dPhi(se11,vs) = bogus
	  dPhi(se12,vs) = bogus
	  dPhi(se22,vs) = bogus
	  dPhi(se1A,vs) = 0
	  dPhi(se2A,vs) = 0
	  dPhi(seAA,vs) = 0
	  fail = 0

	  call TLNormalTB(fail, MBy)
	  if( fail .ne. 0 ) goto 999

	  if( tlCplxApprox .eq. 0 ) call TLasat(asat)

	  call TLatat(atat)

	  if( Ab_bad .eq. 0 ) then
	    call TLLargeTB(fail)
	    if( fail .ne. 0 ) goto 999
	    call TLasab(asab)
	    call TLatab(atab)
* the O(at^2) part at large TB will be subtracted:
	    call TLatat(atat_largeTB)
	  endif

	  if( debuglevel .ge. 4 ) then
	    DSELF "Sigma11_asat =", asat(se11)		ENDL
	    DSELF "        atat =", atat(se11)		ENDL
	    DSELF "        asab =", asab(se11)		ENDL
	    DSELF "        atab =", atab(se11)		ENDL
	    DSELF "atat_largeTB =", atat_largeTB(se11)	ENDL
	    DSELF "Sigma12_asat =", asat(se12)		ENDL
	    DSELF "        atat =", atat(se12)		ENDL
	    DSELF "        asab =", asab(se12)		ENDL
	    DSELF "        atab =", atab(se12)		ENDL
	    DSELF "atat_largeTB =", atat_largeTB(se12)	ENDL
	    DSELF "Sigma22_asat =", asat(se22)		ENDL
	    DSELF "        atat =", atat(se22)		ENDL
	    DSELF "        asab =", asab(se22)		ENDL
	    DSELF "        atab =", atab(se22)		ENDL
	    DSELF "atat_largeTB =", atat_largeTB(se22)	ENDL
	  endif

	  do i = 1, 3
	    dPhi(i,vs) = atat_largeTB(i) - atab(i) -
     &        atat(i) - asat(i) - asab(i)
	  enddo

	  if( runningMT .eq. 1 .or. tM1 .ne. tM2 ) then
	    call TLMSf_gaugeless(fail, MUEc, Xtc)
	    if( fail .ne. 0 ) goto 999
	    call MTrepara(3, dMT, dPhi(1,vs))
	  endif

	  if( looplevel .ge. 3 ) then 
	    call Resum4H(dPhi22)
	    dPhi(se22,vs) = dPhi(se22,vs) - dPhi22
	  endif

999	  continue
	enddo

	do i = 1, 3*(1 + runningMT)
	  if( debuglevel .ge. 3 )
     &      DSELF "interpolating self-energy ", Digit(i) ENDL
	  call TLInterpolate(error, dPhi(i,1), dPhi(i,1), 6)
	enddo

1	if( tlCplxApprox .ne. 0 ) then
	  call TLhrSE(error, Q)
	  if( runningMT .eq. 1 )
     &      call MTrepara(inputmass, dMTs, dPhi(1,1))
	endif

	se2R(h0h0) = se2R(h0h0) +
     &    SA2*dPhi(se11,1) + CA2*dPhi(se22,1) - S2A*dPhi(se12,1)
	se2R(HHHH) = se2R(HHHH) +
     &    CA2*dPhi(se11,1) + SA2*dPhi(se22,1) + S2A*dPhi(se12,1)
	se2R(h0HH) = se2R(h0HH) +
     &    C2A*dPhi(se12,1) + .5D0*S2A*(dPhi(se22,1) - dPhi(se11,1))
	se2R(h0A0) = se2R(h0A0) +
     &    CA*dPhi(se2A,1) - SA*dPhi(se1A,1)
	se2R(HHA0) = se2R(HHA0) +
     &    CA*dPhi(se1A,1) + SA*dPhi(se2A,1)
	se2R(A0A0) = se2R(A0A0) +
     &    dPhi(seAA,1)

	if( error .ne. 0 ) tl_valid = valid

	if( debuglevel .lt. 2 ) return

	DSELF "========================"	ENDL
	DSELF "two-loop self-energies"		ENDL
	DSELF "------------------------"	ENDL
	DSELF "Sigma_11 =", dPhi(se11,1)	ENDL
	DSELF "Sigma_12 =", dPhi(se12,1)	ENDL
	DSELF "Sigma_22 =", dPhi(se22,1)	ENDL
	DSELF "se2R(h0h0) =", se2R(h0h0)	ENDL
	DSELF "se2R(HHHH) =", se2R(HHHH)	ENDL
	DSELF "se2R(h0HH) =", se2R(h0HH)	ENDL
	DSELF "se2R(A0A0) =", se2R(A0A0)	ENDL
	DSELF "se2R(h0A0) =", se2R(h0A0)	ENDL
	DSELF "se2R(HHA0) =", se2R(HHA0)	ENDL
	DSELF "se2R(HmHp) =", se2R(HmHp)	ENDL
	end

#undef __SUBROUTINE__

************************************************************************

	subroutine TLShifts(error)
	implicit none
	integer error

#include "TLps.h"

	ComplexType var(MaxVars)
	RealType phase, db
	integer i, j, n, fail, bit(MaxVars)

	character*3 name(4)
	data name /"M_3", "MUE", "Xt", "Xb"/

	if( tlCplxApprox .eq. 1 ) return

	var(1) = M_3
	var(2) = MUE
	var(3) = Xf(3,3)
	if( tlCplxApprox .eq. 3 .or. tlCplxApprox .eq. 4 )
     &    var(3) = Af(3,3)
	var(4) = Xf(4,3)
	if( tlCplxApprox .eq. 3 .or. tlCplxApprox .eq. 5 )
     &    var(4) = Af(4,3)

	nvars = 0
	if( tlCplxApprox .ge. 3 ) then
	  monomial(1) = 1
	  vdmb(1) = 1
	  do i = 1, MaxVars
	    bit(i) = MaxVars
	    if( Im(var(i)) .ne. 0 ) then
	      phase = atan2(abs(Im(var(i))), Re(var(i)))/pi
	      if( debuglevel .ge. 3 )
     &          DSELF "phase ", name(i), phase*pi ENDL
	      n = ibset(0, nvars)
	      do j = 1, n
	        monomial(j+n) = monomial(j)*phase
	        vdmb(j+n) = vdmb(j) + ishft(vdmb(j), n)
	      enddo
	      bit(i) = nvars
	      nvars = nvars + 1
	      var(i) = abs(var(i))
	    endif
	  enddo
	endif

#define SignedVar(i) (1-2*ibits(vs-1,bit(i),1))*var(i)

	do vs = 1, ibset(0, nvars)
	  pMGl = Re(SignedVar(1))
	  pMGl2 = max(pMGl**2, 1D-14)

* note that {MUEr,Atr,Abr} \equiv Re({MUEc,Atc,Abc}) (see TLps.m)

	  MUEc = SignedVar(2)
	  if( abs(MUEc) .lt. 1D-7 ) MUEc = 1D-7
	  MUEr2 = MUEr**2

	  if( tlCplxApprox .eq. 3 .or. tlCplxApprox .eq. 4 ) then
	    Atc = SignedVar(3)
	    Xtc = Atc - MUEc/TB
	  else
	    Xtc = SignedVar(3)
	    Atc = Xtc + MUEc/TB
	  endif
	  if( tlCplxApprox .eq. 3 .or. tlCplxApprox .eq. 5 ) then
	    Abc = SignedVar(4)
	    Xbc = Abc - MUEc*TB
	  else
	    Xbc = SignedVar(4)
	    Abc = Xbc + MUEc*TB
	  endif

	  if( debuglevel .ge. 3 ) then
	    DSELF "point ", vs	ENDL
	    DSELF name(1), pMGl	ENDL
	    DSELF name(2), MUEc	ENDL
	    DSELF name(3), Xtc	ENDL
	    DSELF name(4), Xbc	ENDL
	  endif

	  BSS2(1) = DSS2(1,4,3)
	  BSS2(2) = DSS2(2,4,3)
	  Deltab = 0
	  MBy = Mf(bM,3)
	  if( bM .ne. 4 ) call TLDeltab(fail)

	  if( debuglevel .gt. 3 ) then
	    DSELF "BSS2(1) =", BSS2(1)
	    DSELF "Deltab  =", Deltab
	    DSELF "MBy     =", MBy
	  endif
	enddo

	if( bM .eq. 4 ) return

	vs = 1

	if( debuglevel .ge. 3 ) DSELF "interpolating BSS" ENDL
	call TLInterpolate(error, DSS2(1,bTR0,3), BSS2(1), Nvr)
c	DSS2(2,bTR0,3) = DSS2(2,4,3)

	if( debuglevel .ge. 3 ) DSELF "interpolating Deltab" ENDL
	call TLInterpolate(error, db, Deltab, Nvr)
	Deltaf(bTR0,3) = db

	if( debuglevel .ge. 3 ) DSELF "interpolating MBy" ENDL
	call TLInterpolate(error, Mf(bTR0,3), MBy, Nvr)
	Mf2(bTR0,3) = Mf(bTR0,3)**2

	if( debuglevel .gt. 3 ) then
	  DPARA "============================"		ENDL
	  DPARA "DSS2(1,bTR0,3) =", DSS2(1,bTR0,3)	ENDL
	  DPARA "Deltaf(bTR0,3) =", Deltaf(bTR0,3)	ENDL
	  DPARA "Mf(bTR0,3) =", Mf(bTR0,3)		ENDL
	endif
	end


************************************************************************

	subroutine TLInterpolate(error, res, data, lddata)
	implicit none
	integer error, lddata
	RealType res, data(lddata,*)

#include "TLps.h"
#define __SUBROUTINE__ "TLInterpolate"

	RealType val(MaxSlots), c
	integer n, nn, i, j, k, p, ivdm(2,MaxSlots,MaxSlots)

	if( nvars .eq. 0 ) then
	  res = data(1,1)
	  if( res .eq. bogus ) Error(error, "Value not defined")
	  return
	endif

	n = 0
	nn = ibset(0, nvars)
	do j = 1, nn
	  if( data(1,j) .ne. bogus ) then
	    n = n + 1
	    val(n) = data(1,j)
	    do i = 1, nn
	      ivdm(1,i,n) = ibits(vdmb(j), i - 1, 1)
	      ivdm(2,i,n) = Delta(i, n)
	    enddo
	  endif
	enddo

	if( debuglevel .ge. 3 ) then
	  do i = 1, nn
	    DSELF "  data =", data(1,i) ENDL
	  enddo
	  DSELF "  n =", n ENDL
	endif

	if( n .lt. nn ) then
	  if( n .eq. 0 ) then
	    res = bogus
	    Error(error, "No points for interpolation")
	    return
	  endif
	  Warning("Missing points for interpolation")
	endif

	do i = n - 1, 1, -1
	  do j = n, i + 1, -1
	    p = ivdm(1,j,i)
	    if( p .ne. 0 ) then
	      do k = 1, n
	        ivdm(1,k,i) = ivdm(1,k,i) - p*ivdm(1,k,j)
	        ivdm(2,k,i) = ivdm(2,k,i) - p*ivdm(2,k,j)
	      enddo
	    endif
	  enddo
	  if( ivdm(1,i,i) .ne. 1 ) then
	    do k = 1, n
	      ivdm(1,k,i) = -ivdm(1,k,i)
	      ivdm(2,k,i) = -ivdm(2,k,i)
	    enddo
	  endif
	enddo

	res = 0

	do i = 1, n
	  do j = 1, i - 1
	    p = ivdm(1,j,i)
	    if( p .ne. 0 ) then
	      do k = 1, n
	        ivdm(1,k,i) = ivdm(1,k,i) - p*ivdm(1,k,j)
	        ivdm(2,k,i) = ivdm(2,k,i) - p*ivdm(2,k,j)
	      enddo
	    endif
	  enddo
	  if( ivdm(1,i,i) .ne. 1 ) then
	    do k = 1, n
	      ivdm(1,k,i) = -ivdm(1,k,i)
	      ivdm(2,k,i) = -ivdm(2,k,i)
	    enddo
	  endif

	  c = 0
	  do j = 1, n
	    c = c + ivdm(2,j,i)*val(j)
	  enddo
	  res = res + c*monomial(i)
	enddo

	if( debuglevel .ge. 3 ) DSELF "  result =", res ENDL
	end

#undef __SUBROUTINE__


************************************************************************

	RealType function mld(m1, m2)
	implicit none
	RealType m1, m2

	if( abs(m1 - m2) .lt. 1D-12 ) then
	  mld = 1
	else
	  mld = m1/(m1 - m2)*log(m1/m2)
	endif
	end

************************************************************************

	subroutine TLDeltab(error)
	implicit none
	integer error

#include "TLps.h"

	RealType db(2), eb(2), dhbhb(2), dMSQ2(2), dvevvev

	RealType mld
	external mld

	error = 0

	call TLNormalTB(error, TLMBrun)
	if( error .ne. 0 ) goto 999

	Q = TLMBscale**2

	call TLStrongShifts(db(1), dhbhb(1), dMSQ2(1))

	eb(1) = C_F/(4*pi)*AlfasMT*2*MUEr*pMGl*TB*
     &    (mld(MSb22, pMGl2) - mld(MSb12, pMGl2))/MSb2diff

	eb(2) = 1/(16*pi**2)*ht2*Atr*MUEr*TB*
     &    (mld(MSt22, MUEr2) - mld(MSt12, MUEr2))/MSt2diff

	call TLLargeTB(error)
	if( error .ne. 0 ) goto 999

	call TLYukawaShifts(db(2), dhbhb(2), dMSQ2(2), dvevvev)

	BSS2(1) = BSS2(1) + dMSQ2(1) + dMSQ2(2)

	if( btest(botResum, 8) ) then
* for FHUncertainties: expanded to first order
	  MBy = TLMBrun*(1 + db(1) + db(2) -
     &      dvevvev - eb(1) - eb(2) - dhbhb(1) - dhbhb(2))
	else
	  MBy = TLMBrun*((1 + db(1) + db(2))/
     &        ((1 + dvevvev)*abs(1 + eb(1) + eb(2))) -
     &      dhbhb(1) - dhbhb(2))
	endif
	Deltab = TLMBrun/MBy - 1

	if( debuglevel .lt. 6 ) return
	DHIGGS "dMSQ2_s   =", dMSQ2(1)	ENDL
	DHIGGS "dMSQ2_yuk =", dMSQ2(2)	ENDL
	DHIGGS "dhbhb_s   =", dhbhb(1)	ENDL
	DHIGGS "dhbhb_yuk =", dhbhb(2)	ENDL
	DHIGGS "db_s      =", db(1)	ENDL
	DHIGGS "db_yuk    =", db(2)	ENDL
	DHIGGS "eb_s      =", eb(1)	ENDL
	DHIGGS "eb_yuk    =", eb(2)	ENDL
	return

999	BSS2(1) = bogus
	Deltab = bogus
	MBy = bogus
	end


************************************************************************

	subroutine TLNormalTB(error, MB_)
	implicit none
	integer error
	RealType MB_

#include "TLps.h"
#define __SUBROUTINE__ "TLNormalTB"

	pCB2 = CB2
	pSB2 = SB2

	pMT = Mf(tM2,3)
	pMT2 = pMT**2

	ht = sqrt2/vev*pMT/SB
	ht2 = ht**2

	call TLMSf(error, StopPara(1),
     &    DSS2(1,3,3), Mf(tM1,3), Xtc,
     &    "normal TB: Mstop^2 negative")

	Xtr = -.5D0*MSt2diff*s2tt/Mf(tM1,3)
	Ytr = Atr + MUEr*TB

	pMB = MB_
	pMB2 = pMB**2

	hb = sqrt2/vev*pMB/CB
	hb2 = hb**2

	call TLMSf(error, SbotPara(1),
     &    BSS2(1), pMB, Xbc,
     &    "normal TB: Msbot^2 negative")

	if( debuglevel .ge. 4 ) then
	  DSELF "normalTB:"			ENDL
	  DSELF "MT     =", pMT			ENDL
	  DSELF "MStop1 =", SignSqrt(MSt12)	ENDL
	  DSELF "MStop2 =", SignSqrt(MSt22)	ENDL
	  DSELF "ctt    =", ctt			ENDL
	  DSELF "stt    =", stt			ENDL
	  DSELF "MB     =", pMB			ENDL
	  DSELF "MSbot1 =", SignSqrt(MSb12)	ENDL
	  DSELF "MSbot2 =", SignSqrt(MSb22)	ENDL
	  DSELF "ctb    =", ctb			ENDL
	  DSELF "stb    =", stb			ENDL
	endif
	end

#undef __SUBROUTINE__


************************************************************************

	subroutine TLLargeTB(error)
	implicit none
	integer error

#include "TLps.h"
#define __SUBROUTINE__ "TLLargeTB"

	RealType TBlarge
	parameter (TBlarge = 1000)

	RealType DSS2large(2)

	pCB2 = 1/(1 + TBlarge**2)
	pSB2 = 1

* Note: for convenience, SB, CB, and C2B keep their old
* (normal TB) values -- they are not used in the TLa?a? terms.

	Xtr = Atr
	Ytr = Atr + MUEr*TBlarge

	ht = sqrt2/vev*pMT
	ht2 = ht**2

c	DSS2large(1) = DSS2(1,3,3) - MZ2*C2B*(.5D0 - 2*Qu*SW2)
c	DSS2large(2) = DSS2(2,3,3) - MZ2*C2B*Qu*SW2
	DSS2large(2) = Re(MSS2(3,3,tU(3)))
	DSS2large(1) = Re(MSS2(3,3,tQ(3))) - DSS2large(2)

	call TLMSf(error, StopPara(1),
     &    DSS2large, Mf(tM1,3), Atc,
     &    "large TB: Mstop^2 negative")

	DSS2large(1) = BSS2(1) - MZ2*C2B*(-.5D0 - 2*Qd*SW2)
	DSS2large(2) = BSS2(2) - MZ2*C2B*Qd*SW2 - pMB2

	call TLMSf(error, SbotPara(1),
     &    DSS2large, pMB, -MUEc/CB,
     &    "large TB: Msbot^2 negative")

	if( debuglevel .ge. 4 ) then
	  DSELF "largeTB:"			ENDL
	  DSELF "MStop1 =", SignSqrt(MSt12)	ENDL
	  DSELF "MStop2 =", SignSqrt(MSt22)	ENDL
	  DSELF "ctt    =", ctt			ENDL
	  DSELF "stt    =", stt			ENDL
	  DSELF "MSbot1 =", SignSqrt(MSb12)	ENDL
	  DSELF "MSbot2 =", SignSqrt(MSb22)	ENDL
	  DSELF "ctb    =", ctb			ENDL
	  DSELF "stb    =", stb			ENDL
	endif
	end

#undef __SUBROUTINE__


************************************************************************

	subroutine TLStrongShifts(db, dhbhb, dMSQ2)
	implicit none
	RealType db, dhbhb, dMSQ2

#include "TLps.h"

	RealType logMT2, logMGl2
	RealType logMSt12, logMSt22, logMSb12, logMSb22
	RealType b0MSt12, b0MSt22
	RealType stop11, stop12, stop22
	RealType sbot11, sbot12, sbot22, top

	logMT2 = log(pMT2/Q)
	logMGl2 = log(pMGl2/Q)
	logMSt12 = log(MSt12/Q)
	logMSt22 = log(MSt22/Q)
	logMSb12 = log(MSb12/Q)
	logMSb22 = log(MSb22/Q)
	b0MSt12 = myB0(MSt12, pMT2, pMGl2, Q)
	b0MSt22 = myB0(MSt22, pMT2, pMGl2, Q)

	db = C_F/(4*pi)*AlfasMT*(3/2D0 - logMGl2 + .5D0*(
     &    MSb12/(pMGl2 - MSb12)*(1 - (logMGl2 - logMSb12)*
     &      (1 + pMGl2/(pMGl2 - MSb12) + 4*pMGl*Abr/MSb2diff)) +
     &    MSb22/(pMGl2 - MSb22)*(1 - (logMGl2 - logMSb22)*
     &      (1 + pMGl2/(pMGl2 - MSb22) - 4*pMGl*Abr/MSb2diff)) ))

	dhbhb = C_F/(4*pi)*AlfasMT*(-4 - 2*logMGl2 - 2/MSb2diff*(
     &    MSb12*(2*logMSb12 -
     &      (1 - pMGl2/MSb12)**2*log(abs(1 - MSb12/pMGl2))) -
     &    MSb22*(2*logMSb22 -
     &      (1 - pMGl2/MSb22)**2*log(abs(1 - MSb22/pMGl2))) ))

	top = pMT*(3*logMT2 + MSDR) +
     &    .5D0/pMT*(2*pMGl2*(logMGl2 - 1) -
     &      MSt12*(logMSt12 - 1) - MSt22*(logMSt22 - 1) +
     &      (pMGl2 + pMT2 - MSt12 - 2*s2tt*pMT*pMGl)*
     &        myB0(pMT2, pMGl2, MSt12, Q) +
     &      (pMGl2 + pMT2 - MSt22 + 2*s2tt*pMT*pMGl)*
     &        myB0(pMT2, pMGl2, MSt22, Q))

	stop11 = MSt12*(3*logMSt12 - 7 - c2tt2*(logMSt12 - 1)) -
     &    s2tt2*MSt22*(logMSt22 - 1) +
     &    2*(pMGl2*(logMGl2 - 1) + pMT2*(logMT2 - 1) +
     &      (MSt12 - pMGl2 - pMT2 + 2*s2tt*pMT*pMGl)*b0MSt12)

	stop22 = MSt22*(3*logMSt22 - 7 - c2tt2*(logMSt22 - 1)) -
     &    s2tt2*MSt12*(logMSt12 - 1) +
     &    2*(pMGl2*(logMGl2 - 1) + pMT2*(logMT2 - 1) +
     &      (MSt22 - pMGl2 - pMT2 - 2*s2tt*pMT*pMGl)*b0MSt22)

	stop12 =
     &    c2tt*s2tt*(MSt22*(1 - logMSt22) - MSt12*(1 - logMSt12)) +
     &    2*pMT*pMGl*c2tt*(b0MSt12 + b0MSt22)

	sbot11 = MSb12*(3*logMSb12 - 3 - c2tb2*(logMSb12 - 1)) -
     &    s2tb2*MSb22*(logMSb22 - 1) - 6*pMGl2 -
     &    2*(MSb12 - 2*pMGl2)*logMGl2 -
     &    2*(MSb12 - pMGl2)**2/MSb12*log(abs(1 - MSb12/pMGl2))

	sbot22 = MSb22*(3*logMSb22 - 3 - c2tb2*(logMSb22 - 1)) -
     &    s2tb2*MSb12*(logMSb12 - 1) - 6*pMGl2 -
     &    2*(MSb22 - 2*pMGl2)*logMGl2 -
     &    2*(MSb22 - pMGl2)**2/MSb22*log(abs(1 - MSb22/pMGl2))

	sbot12 =
     &    c2tb*s2tb*(MSb22*(1 - logMSb22) - MSb12*(1 - logMSb12))

	dMSQ2 = C_F/(4*pi)*AlfasMT*(
     &    ctt2*stop11 + stt2*stop22 - 2*stt*ctt*stop12 -
     &    ctb2*sbot11 - stb2*sbot22 + 2*stb*ctb*sbot12 -
     &    2*pMT*top )
	end


************************************************************************

	subroutine TLYukawaShifts(db, dhbhb, dMSQ2, dvevvev)
	implicit none
	RealType db, dhbhb, dMSQ2, dvevvev

#include "TLps.h"

	RealType logMT2, logMA02, top, deltaM, t1, t2
	RealType stop11, stop12, stop22
	RealType sbot11, sbot12, sbot22

	RealType Pi11t, Pi22t, Pi12t
	RealType Pi11b, Pi22b, Pi12b, Pitt
	external Pi11t, Pi22t, Pi12t, Pi11b, Pi22b, Pi12b, Pitt

	RealType x, y, fW
	fW(x, y) = .25D0*(x + y) - .5D0*x*y/(x - y)*log(x/y)

	logMT2 = log(pMT2/Q)
	logMA02 = log(MA02/Q)

	top = Pitt(pMT2)

	stop11 = Pi11t(MSt12)
	stop22 = Pi22t(MSt22)
	stop12 = .5D0*(Pi12t(MSt12) + Pi12t(MSt22))

	sbot11 = Pi11b(MSb12)
	sbot22 = Pi22b(MSb22)
	sbot12 = .5D0*(Pi12b(MSb12) + Pi12b(MSb22))

	deltaM = pMT2 - MA02
	if( abs(deltaM) .lt. 1D-13 ) then
	  t1 = -logMT2
	  t2 = -.5D0
	else
	  t1 = 1 - (pMT2*logMT2 - MA02*logMA02)/deltaM
	  t2 = pMT2*(deltaM - pMT2*(logMT2 - logMA02))/deltaM**2
	endif

	db = 1/(32*pi**2)*(
     &    ht2*(6/4D0*logMT2 - 5/4D0 + 2*t1 +
     &         2*MUEr2/MSt2diff*(myB0(0D0, MUEr2, MSt12, Q) -
     &                           myB0(0D0, MUEr2, MSt22, Q))) +
     &    hb2*(.5D0*t2 + 3/4D0 - 3/2D0*logMA02 +
     &         myB1(0D0, MUEr2, MSb12, Q) +
     &         myB1(0D0, MUEr2, MSb22, Q)) +
     &    (hb2*ctt2 + ht2*stt2)*myB1(0D0, MUEr2, MSt12, Q) +
     &    (hb2*stt2 + ht2*ctt2)*myB1(0D0, MUEr2, MSt22, Q) )

	dvevvev = Nc/(16*pi**2)*ht2*(.5D0*logMT2 - .25D0 -
     &    (ctt2*ctb2*fW(MSt12, MSb12) + 
     &     ctt2*stb2*fW(MSt12, MSb22) +
     &     stt2*ctb2*fW(MSt22, MSb12) +
     &     stt2*stb2*fW(MSt22, MSb22))/pMT2)

	dhbhb = (sbot22 - sbot11 - 2*c2tb/s2tb*sbot12)/
     &      (16*pi**2*MSb2diff) -
     &    dvevvev

	dMSQ2 = 1/(16*pi**2)*(
     &    ctt2*stop11 + stt2*stop22 - s2tt*stop12 -
     &    ctb2*sbot11 - stb2*sbot22 + s2tb*sbot12 -
     &    2*pMT*top )
	end


************************************************************************

#undef StopPara

	subroutine TLMSf(error, StopPara, M11, Mf_, Xf_, msg)
	implicit none
	integer error
	RealType StopPara(*), M11(2), Mf_
	ComplexType Xf_
	character*(*) msg

#include "FH.h"
#define __SUBROUTINE__ "TLMSf"

	integer fail
	RealType M12, M12sq, delta, h

	M12sq = Mf_**2*Sq(Xf_)
	h = .5D0*(M11(1) + sign(sqrt(M11(1)**2 + 4*M12sq), M11(1)))
	if( h .ne. 0 ) h = 1/h
	delta = h*M12sq
	stt = Mf_*Re(Xf_)*h/sqrt(delta*h + 1)
	stt2 = stt**2
	ctt2 = (1 - stt)*(1 + stt)
	ctt = sqrt(ctt2)
	c2tt = (ctt - stt)*(ctt + stt)
	c2tt2 = c2tt**2
	s2tt = 2*stt*ctt
	s2tt2 = s2tt**2

	MSt2diff = -M11(1) - 2*delta + 1D-15
	MSt22 = M11(2) + Mf_**2
	MSt12 = MSt22 + M11(1) + delta
	MSt2sum = 2*MSt22 + M11(1)
	MSt22 = MSt22 - delta

	if( MSt12 .lt. 0 .or. MSt22 .lt. 0 ) Error(error, msg)

	if( debuglevel .lt. 4 ) return

	DHIGGS msg(1:index(msg, "^")-1), sqrt(MSt12), sqrt(MSt22) ENDL
	DHIGGS "  MSq2diff =", MSt2diff			ENDL
	DHIGGS "  M11 =", M11(1) + M11(2)		ENDL
	DHIGGS "  M22 =", M11(2)			ENDL
	DHIGGS "  Xf  =", Re(Xf_)			ENDL
	DHIGGS "  Mf  =", Mf_				ENDL
	end
