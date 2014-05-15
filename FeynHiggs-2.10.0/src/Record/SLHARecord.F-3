* SLHARecord.F
* initialize FeynHiggs record from an SLHA data structure
* this file is part of FeynHiggs
* last modified 28 Jun 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHSLHARecord(error, record, slhadata)
	implicit none
	integer error
	ComplexType slhadata(*)

#include "FH.h"
#include "FHRecord.h"
#include "SLHA.h"
#include "looptools.h"
#define __SUBROUTINE__ "FHSLHAtoRecord"

	RecordDecl(record)

#define Origin(x) if(debuglevel.ge.1) DSLHA x ENDL
#define Origin2(x,y) if(debuglevel.ge.1) DSLHA x,y ENDL

	integer t, nQ, nU, qua, qua2
	RealType extpar_Q

	character*12 qsf_from
	data qsf_from /"QSf(t) from "/

#include "DeltaAlfaTop.h"

	error = 0

	call FHClearRecord(record)

	extpar_Q = Re(ExtPar_Q)
* defaults to the GUT scale if not given
c	if( extpar_Q .eq. invalid ) extpar_Q = 2D16
	if( extpar_Q .eq. invalid ) then
	  extpar_Q = -1
	  Warning("ExtPar_Q taken to be sqrt(M3SQ M3SU)")
	endif

#define R(i) record(i,iLower)

	R(iinvAlfaMZ) = Re(SMInputs_invAlfaMZ) + 1/Alfa0*
     &    DeltaAlfaTop(Re(SMInputs_Mt)**2, Re(SMInputs_MZ)**2)
	R(iAlfasMZ) = Re(SMInputs_AlfasMZ)
	R(iGF) = Re(SMInputs_GF)

	R(iME) = Re(SMInputs_Me)
	R(iMU) = Re(SMInputs_Mu)
	R(iMD) = Re(SMInputs_Md)
	R(iMM) = Re(SMInputs_Mmu)
	R(iMC) = Re(SMInputs_Mc)
	R(iMS) = Re(SMInputs_Ms)
	R(iML) = Re(SMInputs_Mtau)
	R(iMT) = Re(SMInputs_Mt)
	R(iMB) = Re(SMInputs_Mb)

	R(iMW) = -1
	R(iMZ) = Re(SMInputs_MZ)

	if( SLHAValid(BlockVCKMIn(1), LengthVCKMIn) .ne. 0 ) then
	  Origin("CKM from VCKMIn")
	  R(iCKMlambda) = Re(VCKMIn_lambda)
	  R(iCKMA) = Re(VCKMIn_A)
	  R(iCKMrhobar) = Re(VCKMIn_rhobar)
	  R(iCKMetabar) = Re(VCKMIn_etabar)
	endif

	call SMDefaults(R(iinvAlfaMZ), R(iAlfasMZ), R(iGF),
     &    R(iME), R(iMU), R(iMD),
     &    R(iMM), R(iMC), R(iMS),
     &    R(iML),         R(iMB),
     &    R(iMW), R(iMZ),
     &    R(iCKMlambda), R(iCKMA), R(iCKMrhobar), R(iCKMetabar))

	if( int(ModSel_Content) .gt. 0 )
     &    Warning("SLHA record indicates model other than MSSM")
	if( int(ModSel_RPV) .gt. 0 )
     &    Warning("SLHA record indicates RPV")

	if( DMass_Q .ne. invalid ) then
	  Origin("scalefactor from DMass_Q")
	  R(iscalefactor) = Re(DMass_Q)/Re(SMInputs_Mt)
	endif

	if( MinPar_TB .ne. invalid ) then
	  Origin("TB from MinPar_TB")
	  R(iTB) = Re(MinPar_TB)
	else if( HMix_TB .ne. invalid ) then
	  Origin("TB from HMix_TB")
	  R(iTB) = Re(HMix_TB)
	else if( ExtPar_TB .ne. invalid ) then
	  Origin("TB from ExtPar_TB")
	  R(iTB) = Re(ExtPar_TB)
	endif

	if( Mass_MA0 .ne. invalid ) then
	  Origin("MA0 from Mass_MA0")
	  R(iMA0) = Re(Mass_MA0)
c	  R(iQH) = 0
	else if( ExtPar_MA02 .ne. invalid ) then
	  Origin("MA0 from ExtPar_MA02")
	  R(iMA0) = sqrt(Re(ExtPar_MA02))
c	  R(iQH) = extpar_Q
	else if( ExtPar_MA0 .ne. invalid ) then
	  Origin("MA0 from ExtPar_MA0")
	  R(iMA0) = Re(ExtPar_MA0)
c	  R(iQH) = extpar_Q
	endif

	if( Mass_MHp .ne. invalid ) then
	  Origin("MHp from Mass_MHp")
	  R(iMHp) = Re(Mass_MHp)
c	  R(iQH) = 0
	else if( ExtPar_MHp .ne. invalid ) then
	  Origin("MHp from ExtPar_MHp")
	  R(iMHp) = Re(ExtPar_MHp)
c	  R(iQH) = extpar_Q
	endif

	if( HMix_MUE .ne. invalid ) then
	  Origin("MUE from HMix_MUE")
	  R(iRe(iMUE)) = Re(HMix_MUE)
	  R(iIm(iMUE)) = Im(HMix_MUE)
	else if( ExtPar_MUE .ne. invalid ) then
	  Origin("MUE from ExtPar_MUE")
	  R(iRe(iMUE)) = Re(ExtPar_MUE)
	  R(iIm(iMUE)) = Im(ExtPar_MUE)
	endif

	if( MSoft_M1 .ne. invalid ) then
	  Origin("M_1 from MSoft_M1")
	  R(iRe(iM1)) = Re(MSoft_M1)
	  R(iIm(iM1)) = Im(MSoft_M1)
	else if( ExtPar_M1 .ne. invalid ) then
	  Origin("M_1 from ExtPar_M1")
	  R(iRe(iM1)) = Re(ExtPar_M1)
	  R(iIm(iM1)) = Im(ExtPar_M1)
	endif

	if( MSoft_M2 .ne. invalid ) then
	  Origin("M_2 from MSoft_M2")
	  R(iRe(iM2)) = Re(MSoft_M2)
	  R(iIm(iM2)) = Im(MSoft_M2)
	else if( ExtPar_M2 .ne. invalid ) then
	  Origin("M_2 from ExtPar_M2")
	  R(iRe(iM2)) = Re(ExtPar_M2)
	  R(iIm(iM2)) = Im(ExtPar_M2)
	endif

	if( MSoft_M3 .ne. invalid ) then
	  Origin("M_3 from MSoft_M3")
	  R(iRe(iM3)) = Re(MSoft_M3)
	  R(iIm(iM3)) = Im(MSoft_M3)
	else if( Mass_MGl .ne. invalid ) then
	  Origin("M_3 from Mass_MGl")
	  R(iRe(iM3)) = Re(Mass_MGl)
	  R(iIm(iM3)) = Im(Mass_MGl)
	else if( ExtPar_M3 .ne. invalid ) then
	  Origin("M_3 from ExtPar_M3")
	  R(iRe(iM3)) = Re(ExtPar_M3)
	  R(iIm(iM3)) = Im(ExtPar_M3)
	endif

	do t = 2, 4
	  nQ = tQ(t)
	  nU = tU(t)
	  qsf_from(5:5) = Digit(t)

	  qua = 0

	  call FVRecord(record, qua, t,
     &      MSS2_MSS2(1,1,nQ),
     &      MSS2_MSS2(1,1,nU),
     &      Tf_Tf(1,1,t), "MSS2_MSS2", "Tf_Tf")
	  if( qua .ne. 0 ) then
	    Origin(qsf_from//"MSS2_Q")
	    R(iQSf(t)) = Re(MSS2_Q(nQ))
	  endif
	  if( qua .eq. 7 ) goto 1
	  call FVRecord(record, qua2, t,
     &      MSS2In_MSS2(1,1,nQ),
     &      MSS2In_MSS2(1,1,nU),
     &      TfIn_Tf(1,1,t), "MSS2In_MSS2", "TfIn_Tf")
	  if( qua2 .eq. 7 ) goto 1
	  qua = ior(qua, qua2)

	  if( btest(qua, 0) ) then
	  else if( SLHAValid(MSoft_MSS(1,nQ), 3) .ne. 0 ) then
	    Origin("MS"//cMSS(nQ)//" from MSoft_MSS2")
	    Origin(qsf_from//"MSoft_Q")
	    R(iMSS(nQ,1)) = Re(MSoft_MSS(1,nQ))
	    R(iMSS(nQ,2)) = Re(MSoft_MSS(2,nQ))
	    R(iMSS(nQ,3)) = Re(MSoft_MSS(3,nQ))
	    R(iQSf(t)) = Re(MSoft_Q)
	  else if( SLHAValid(ExtPar_MSS(1,nQ), 3) .ne. 0 ) then
	    Origin("MS"//cMSS(nQ)//" from ExtPar_MSS")
	    Origin(qsf_from//"ExtPar_Q")
	    R(iMSS(nQ,1)) = Re(ExtPar_MSS(1,nQ))
	    R(iMSS(nQ,2)) = Re(ExtPar_MSS(2,nQ))
	    R(iMSS(nQ,3)) = Re(ExtPar_MSS(3,nQ))
	    R(iQSf(t)) = extpar_Q
	  endif

	  if( btest(qua, 1) ) then
	  else if( SLHAValid(MSoft_MSS(1,nU), 3) .ne. 0 ) then
	    Origin("MS"//cMSS(nU)//" from MSoft_MSS2")
	    Origin(qsf_from//"MSoft_Q")
	    R(iMSS(nU,1)) = Re(MSoft_MSS(1,nU))
	    R(iMSS(nU,2)) = Re(MSoft_MSS(2,nU))
	    R(iMSS(nU,3)) = Re(MSoft_MSS(3,nU))
	    R(iQSf(t)) = Re(MSoft_Q)
	  else if( SLHAValid(ExtPar_MSS(1,nU), 3) .ne. 0 ) then
	    Origin("MS"//cMSS(nU)//" from ExtPar_MSS")
	    Origin(qsf_from//"ExtPar_Q")
	    R(iMSS(nU,1)) = Re(ExtPar_MSS(1,nU))
	    R(iMSS(nU,2)) = Re(ExtPar_MSS(2,nU))
	    R(iMSS(nU,3)) = Re(ExtPar_MSS(3,nU))
	    R(iQSf(t)) = extpar_Q
	  endif

	  if( btest(qua, 2) ) then
	  else if( Af_Af(3,3,t) .ne. invalid ) then
	    Origin("A"//cAf(t)//" from Af_Af")
	    Origin(qsf_from//"Af_Q")
	    R(iRe(iAf(t,1))) = Re(Af_Af(1,1,t))
	    R(iIm(iAf(t,1))) = Im(Af_Af(1,1,t))
	    R(iRe(iAf(t,2))) = Re(Af_Af(2,2,t))
	    R(iIm(iAf(t,2))) = Im(Af_Af(2,2,t))
	    R(iRe(iAf(t,3))) = Re(Af_Af(3,3,t))
	    R(iIm(iAf(t,3))) = Im(Af_Af(3,3,t))
	    R(iQSf(t)) = Re(Af_Q(t))
	  else if( SLHAValid(ExtPar_Af(1), 3) .ne. 0 ) then
	    Origin("A"//cAf(t)//" from ExtPar_Af")
	    Origin(qsf_from//"ExtPar_Q")
	    R(iRe(iAf(t,1))) = Re(ExtPar_Af(t))
	    R(iIm(iAf(t,1))) = Im(ExtPar_Af(t))
	    R(iRe(iAf(t,2))) = Re(ExtPar_Af(t))
	    R(iIm(iAf(t,2))) = Im(ExtPar_Af(t))
	    R(iRe(iAf(t,3))) = Re(ExtPar_Af(t))
	    R(iIm(iAf(t,3))) = Im(ExtPar_Af(t))
	    R(iQSf(t)) = extpar_Q
	  endif

1	  continue
	enddo
	end


************************************************************************

	subroutine FVRecord(record, qua, t, MSQ2_, MSU2_, TU_,
     &    Mname, Tname)
	implicit none
	integer qua, t
	ComplexType MSQ2_(3,3), MSU2_(3,3), TU_(3,3)
	character*(*) Mname, Tname

#include "FH.h"
#include "FHRecord.h"

	RecordDecl(record)

	integer nQ, nU, g, j, o, i
	RealType x
	ComplexType c
	ComplexType MSL2_(3,3), MSR2_(3,3), Kf_(3,3)
	RealType MSL_(3), MSR_(3)

	nQ = tQ(t)
	nU = tU(t)

	if( debuglevel .ge. 4 ) then
	  DSLHA "============================"		ENDL
	  DSLHA "FVRecord input from ", Mname, ", ", Tname, ":" ENDL
	  DSLHA "----------------------------"		ENDL
	  DSLHA "MS"//cMSS(nQ)//"2(1,1) =", MSQ2_(1,1)	ENDL
	  DSLHA "    (2,2) =", MSQ2_(2,2)		ENDL
	  DSLHA "    (3,3) =", MSQ2_(3,3)		ENDL
	  DSLHA "    (1,2) =", MSQ2_(1,2)		ENDL
	  DSLHA "    (1,3) =", MSQ2_(1,3)		ENDL
	  DSLHA "    (2,3) =", MSQ2_(2,3)		ENDL
	  DSLHA "MS"//cMSS(nU)//"2(1,1) =", MSU2_(1,1)	ENDL
	  DSLHA "    (2,2) =", MSU2_(2,2)		ENDL
	  DSLHA "    (3,3) =", MSU2_(3,3)		ENDL
	  DSLHA "    (1,2) =", MSU2_(1,2)		ENDL
	  DSLHA "    (1,3) =", MSU2_(1,3)		ENDL
	  DSLHA "    (2,3) =", MSU2_(2,3)		ENDL
	  DSLHA "  T"//cAf(t)//"(1,1) =", TU_(1,1)	ENDL
	  DSLHA "    (2,2) =", TU_(2,2)			ENDL
	  DSLHA "    (3,3) =", TU_(3,3)			ENDL
	  DSLHA "    (1,2) =", TU_(1,2)			ENDL
	  DSLHA "    (1,3) =", TU_(1,3)			ENDL
	  DSLHA "    (2,1) =", TU_(2,1)			ENDL
	  DSLHA "    (2,3) =", TU_(2,3)			ENDL
	  DSLHA "    (3,1) =", TU_(3,1)			ENDL
	  DSLHA "    (3,2) =", TU_(3,2)			ENDL
	endif

#define diag(j,g) iand(ibset(0, j-1), ibset(0, g-1))

	qua = 0

	o = 0
	do g = 1, 3
	  do j = 1, g
	    MSL2_(j,g) = 0
	    if( MSQ2_(j,g) .ne. invalid ) then
	      MSL2_(j,g) = MSQ2_(j,g)
	      o = ior(o, diag(j,g))
	    endif
	  enddo
	  MSL_(g) = SignSqrt(MSL2_(g,g))
	enddo

	if( o .eq. 7 ) then
	  Origin2("MS"//cMSS(nQ)//" from ",Mname)
	  R(iMSS(nQ,1)) = MSL_(1)
	  R(iMSS(nQ,2)) = MSL_(2)
	  R(iMSS(nQ,3)) = MSL_(3)
	  qua = ibset(qua, 0)
	endif

	o = 0
	do g = 1, 3
	  do j = 1, g
	    MSR2_(j,g) = 0
	    if( MSU2_(j,g) .ne. invalid ) then
	      MSR2_(j,g) = MSU2_(j,g)
	      o = ior(o, diag(j,g))
	    endif
	  enddo
	  MSR_(g) = SignSqrt(MSR2_(g,g))
	enddo

	if( o .eq. 7 ) then
	  Origin2("MS"//cMSS(nU)//" from ",Mname)
	  R(iMSS(nU,1)) = MSR_(1)
	  R(iMSS(nU,2)) = MSR_(2)
	  R(iMSS(nU,3)) = MSR_(3)
	  qua = ibset(qua, 1)
	endif

	o = 0

c	vev = 1/sqrt(sqrt2*R(iGF))
c	CB = 1/sqrt(1 + R(iTB)**2)
c	x = vev/sqrt2*CB
	x = 1/sqrt(2*sqrt2*R(iGF)*(1 + R(iTB)**2))
	if( iand(t, 1) .eq. 0 ) x = x*R(iTB)
	do g = 1, 3
	  do j = 1, 3
	    Kf_(j,g) = 0
	    if( TU_(j,g) .ne. invalid ) then
	      Kf_(j,g) = x*Conjugate(TU_(g,j))
	      o = ior(o, diag(j,g))
	    endif
	  enddo
	enddo

	if( o .eq. 7 ) then
	  Origin2("K"//cAf(t)//" from ",Tname)
	  c = Kf_(1,1)/R(iMf(t,1))
	  R(iRe(iAf(t,1))) = Re(c)
	  R(iIm(iAf(t,1))) = Im(c)
	  c = Kf_(2,2)/R(iMf(t,2))
	  R(iRe(iAf(t,2))) = Re(c)
	  R(iIm(iAf(t,2))) = Im(c)
	  c = Kf_(3,3)/R(iMf(t,3))
	  R(iRe(iAf(t,3))) = Re(c)
	  R(iIm(iAf(t,3))) = Im(c)
	  qua = ibset(qua, 2)
	endif

	if( debuglevel .ge. 3 ) DSLHA "qua =", qua ENDL

	if( qua .ne. 7 ) return

	Origin2("deltaLL,RR("//Digit(t)//") from ",Mname)
	Origin2("deltaLR,RL("//Digit(t)//") from ",Tname)

* note: the following assignments assume that 12,23,13 are
* consecutive in the RecordIndex.h enumeration

	i = ideltaLLL12 + ifQ(t, ideltaQLL12 - ideltaLLL12)
	c = MSL2_(1,2)/(MSL_(1)*MSL_(2))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)
	i = i + ideltaLLL23 - ideltaLLL12
	c = MSL2_(2,3)/(MSL_(2)*MSL_(3))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)
	i = i + ideltaLLL13 - ideltaLLL23
	c = MSL2_(1,3)/(MSL_(1)*MSL_(3))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)

	i = ideltaELR12 +
     &    ifQ(t, ideltaULR12 - ideltaELR12) +
     &    ifB(t, ideltaDLR12 - ideltaULR12)
	c = Kf_(1,2)/(MSL_(1)*MSR_(2))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)
	i = i + ideltaELR23 - ideltaELR12
	c = Kf_(2,3)/(MSL_(2)*MSR_(3))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)
	i = i + ideltaELR13 - ideltaELR23
	c = Kf_(1,3)/(MSL_(1)*MSR_(3))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)

	i = i + ideltaERL12 - ideltaELR13
	c = Conjugate(Kf_(2,1))/(MSR_(1)*MSL_(2))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)
	i = i + ideltaERL23 - ideltaERL12
	c = Conjugate(Kf_(3,2))/(MSR_(2)*MSL_(3))
	R(iRe(i+1)) = Re(c)
	R(iIm(i+1)) = Im(c)
	i = i + ideltaERL13 - ideltaERL23
	c = Conjugate(Kf_(3,1))/(MSR_(1)*MSL_(3))
	R(iRe(i+2)) = Re(c)
	R(iIm(i+2)) = Im(c)

	i = i + ideltaERR12 - ideltaERL13
	c = MSR2_(1,2)/(MSR_(1)*MSR_(2))
	R(iRe(i)) = Re(c)
	R(iIm(i)) = Im(c)
	i = i + ideltaERR23 - ideltaERR12
	c = MSR2_(2,3)/(MSR_(2)*MSR_(3))
	R(iRe(i+1)) = Re(c)
	R(iIm(i+1)) = Im(c)
	i = i + ideltaERR13 - ideltaERR23
	c = MSR2_(1,3)/(MSR_(1)*MSR_(3))
	R(iRe(i+2)) = Re(c)
	R(iIm(i+2)) = Im(c)
	end

