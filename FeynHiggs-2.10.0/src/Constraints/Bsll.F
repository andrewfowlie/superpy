* Bsll.F
* computes the decay rate of B_s -> l+ l-
* the calculation follows closely hep-ph/0112305
* this file is part of FeynHiggs
* last modified 16 Jun 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine Bsll(bsllMSSM, bsllSM)
	implicit none
	RealType bsllMSSM, bsllSM

#include "FH.h"
#include "looptools.h"

	ComplexType CAL, CAR, CALSM, CALHp, CARHp
	ComplexType CALCha, CARCha, CALGlu, CARGlu
	ComplexType CALNeu, CARNeu, CSLNeu, CPLNeu, CSRNeu, CPRNeu
	ComplexType FS, FP, FA
	RealType beta2, fac, asMT
	integer nfdummy

	integer Cha5, Cha6, Neu5, Neu6
	integer Sfe5, All5, All6, Ind1, Ind2, Ind3
	ComplexType dup1, dup2, dup3, dup4, dup5, dup6, dup7
	ComplexType dup8, dup9, dup10, dup11, dup12, dup13
	ComplexType dup14, dup15, dup16, dup17, dup18, dup19
	ComplexType dup20, dup21, dup22, dup23, dup24
	ComplexType dup25, dup26, dup27, dup28, dup29
	ComplexType dup30, dup31, dup32, dup33, dup34
	ComplexType dup35, dup36, dup37, dup38, dup39
	ComplexType dup40, dup41, dup42, dup43, dup44
	ComplexType dup45, dup46
	ComplexType tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
	ComplexType tmp7, tmp8, tmp9

	integer l, q
	parameter (l = 2)	! muon
	parameter (q = 2)	! strange

	RealType MBs, fBs, tauBs, Mshat
	parameter (MBs = 5367.5D-3)
	parameter (fBs = 240D-3)
	parameter (tauBs = 2.227D12)
	parameter (Mshat = 90/4400D0)

	call AlphaS(asMT, nfdummy, Mf(tT,3))

#include "bsll.h"

	CAL = CALSM + CALHp + CALCha + CALNeu + CALGlu
	CAR =         CARHp + CARCha + CARNeu + CARGlu
#define CSL CSLNeu
#define CSR CSRNeu
#define CPL CPLNeu
#define CPR CPRNeu

	if( debuglevel .ge. 2 ) then
	  DCONST "Bsll Wilson coefficients before running" ENDL
	  DCONST "CALSM =", CALSM	ENDL
	  DCONST "CAL =", CAL		ENDL
	  DCONST "CAR =", CAR		ENDL
	  DCONST "CSL =", CSL		ENDL
	  DCONST "CSR =", CSR		ENDL
	  DCONST "CPL =", CPL		ENDL
	  DCONST "CPR =", CPR		ENDL
	endif

	beta2 = 1 - 4/MBs**2*Mf2(2,l)
#if 0
	fac = 1/(16*pi**3)*MBs*fBs**2*tauBs*(GF*Alfa0/SW2)**2*
     &    Sq(CKM(3,3)*CKMC(3,q))*
     &    sqrt(beta2)
#else
	fac = 1/(16*pi**3)*MBs*fBs**2*tauBs*(GF*Alfa0/SW2)**2*
     &    Sq(CKM(3,3)*CKMC(3,q))*
     &    sqrt(beta2)
c     &    2*MW2*SW2**2/Alfa0**2*
c	fac = 1/(16*pi**3)*MBs*fBs**2*tauBs/GF**2*
c     &    sqrt(beta2)
c	PRINT *, "fac=", fac
#endif

	bsllSM = fac*Sq(Mf(2,l)*CALSM)

	FS = .5D0*MBs**2/(1 + Mshat)*(CSR - CSL)/Mf(4,3)
	FP = .5D0*MBs**2/(1 + Mshat)*(CPR - CPL)/Mf(4,3)
	FA = .5D0*(CAR - CAL)
	bsllMSSM = fac*(beta2*Sq(FS) + Sq(FP + 2*Mf(2,l)*FA))
	end

