* Para.F
* initialization routines for the various parameter blocks
* this file is part of FeynHiggs
* last modified 22 Nov 13 th

* based on mssm_ini.F from FormCalc:
* Christian Schappacher, 17 Feb 99
* version of 15 Jun 01 th
* (almost) completely rewritten to be independent
* of any external libraries and
* adapted to Higgs mass routines by Markus Frank.
* 01.10.2002 sh: inclusion of O(als alb) resummation

#include "externals.h"
#include "types.h"
#include "debug.h"

#define __SUBROUTINE__ "FHSetPara"


	subroutine BeginPara(error, scalefactor_, MT_, TB_, M_3_)
	implicit none
	integer error
	RealType scalefactor_, MT_, TB_
	ComplexType M_3_

#include "FH.h"
#include "FHRecord.h"
#include "looptools.h"

	external recordnamesini

	integer i, t, ty
	RealType i3part, qpart, tmp

	RealType Run
	external Run

	para_valid = 0
	sf_valid = 0
	tl_valid = 0
	higgs_valid = 0
	coup_valid = 0

	error = 0
	CheckFlags(error)

	if( debuglevel .ge. 2 ) then
	  DPARA "============================"	ENDL
	  DPARA "BeginPara input parameters:"	ENDL
	  DPARA "----------------------------"	ENDL
	  DTAG(iscalefactor, scalefactor_)
	  DTAG(iMT, MT_)
	  DTAG(iTB, TB_)
	  DTAGre(iM3, M_3_)
	  DTAGimz(iM3, M_3_)
	endif

	fv = 0

	do i = 1, semax
	  se2R(i) = 0
	enddo
	do i = 1, semax
	  seAdd(i) = 0
	enddo

	Qf(1) = 0
	Qf(2) = Qe
	Qf(3) = Qu
	Qf(4) = Qd

	QSf(2) = 0
	QSf(3) = 0
	QSf(4) = 0

	MT = MT_
	if( MT .lt. 0 ) Error(error, "invalid MT")
	MT2 = MT**2

	do ty = 5, SfSlots
	  t = int(ibits(SfIni, 3*(SfSlots - ty), 3))
	  Mf(ty,1) = Mf(t,1)
	  Mf(ty,2) = Mf(t,2)
	  Mf(ty,3) = Mf(t,3)
	  Mf2(ty,1) = Mf2(t,1)
	  Mf2(ty,2) = Mf2(t,2)
	  Mf2(ty,3) = Mf2(t,3)
	enddo

	call ltini
	scalefactor = scalefactor_
	if( scalefactor .lt. 0 ) scalefactor = 1
	call setmudim((scalefactor*MT)**2)

	AlfasMT = AlfasMZ/(1 + (11 - 10/3D0)/(4*pi)*AlfasMZ*log(MT2/MZ2))
c	call AlphaS(AlfasMT, nf, MT)
	GSMT = sqrt(4*pi*AlfasMT)

	if( tlCplxApprox .eq. 1 ) then
	  Mf(tT,3) = MT/(1 + 4/(3*pi)*AlfasMT)
	  Alfat = (Mf(tT,3)/(vev/sqrt2))**2/(4*pi)
	else
	  tmp = 0
	  do while( abs(Mf(tT,3) - tmp) .gt. 1D-5 )
	    tmp = Mf(tT,3)
	    Alfat = (Mf(tT,3)/(vev/sqrt2))**2/(4*pi)
	    Mf(tT,3) = MT/(1 - Alfat/(2*pi) + 4/(3*pi)*AlfasMT)
	  enddo
	endif
	Mf2(tT,3) = Mf(tT,3)**2

	MB_MT = MB*Run(MB, MT)

	M_3 = M_3_
	MGl = abs(M_3)
	MGl2 = MGl**2

	TB = TB_
	if( TB .lt. 0 ) Error(error, "invalid TB")
	TB2 = TB**2
	CB2 = 1/(1 + TB2)
	SB2 = CB2*TB2
	CB = sqrt(CB2)
	SB = CB*TB
	S2B = 2*TB*CB2
	C2B = (CB - SB)*(CB + SB)

	SCB(2) = CB
	SCB(3) = SB
	SCB(4) = CB

	i3part = .5D0*MZ2*C2B
	qpart = SW2*MZ2*C2B
	DSf(2,1) = 0
	DSf(1,1) = +i3part

	DSf(2,2) = Qe*qpart
	DSf(1,2) = -i3part - DSf(2,2)

	DSf(2,3) = Qu*qpart
	DSf(1,3) = +i3part - DSf(2,3)

	DSf(2,4) = Qd*qpart
	DSf(2,5) = DSf(2,4)
	DSf(1,4) = -i3part - DSf(2,4)
	DSf(1,5) = DSf(1,4)
	end


************************************************************************

	subroutine HiggsMasses(error, MA0_, MHp_)
	implicit none
	integer error
	RealType MA0_, MHp_

#include "FH.h"
#include "FHRecord.h"

	RealType avgMH2, deltaMH2

	if( debuglevel .ge. 2 ) then
	  DPARA "============================"	ENDL
	  DPARA "HiggsMasses input:"		ENDL
	  DPARA "----------------------------"	ENDL
	  DTAGm(iMA0, MA0_)
	  DTAGm(iMHp, MHp_)
	endif

	if( MA0_ .ge. 0 ) then
	  MA0 = max(MA0_, 1D-3)
	  MA02 = MA0**2
	  MHp2 = MA02 + MW2
	  MHp = sqrt(MHp2)
	  inputmass = 3
	else if( MHp_ .gt. MW ) then
	  MHp = MHp_
	  MHp2 = MHp**2
	  MA02 = (MHp - MW)*(MHp + MW)
	  MA0 = sqrt(MA02)
	  inputmass = 4
	else
	  Error(error, "neither MA0 nor MHp given")
	  return
	endif

* these are the tree-level Higgs masses:
	avgMH2 = .5D0*(MA02 + MZ2)
	deltaMH2 = MA0*MZ*C2B
	deltaMH2 = sqrt((avgMH2 - deltaMH2)*(avgMH2 + deltaMH2))
* the 1D-6 is so we don't get killed by e.g. B0(0, Mh02, Mh02) later
	Mh02 = max(avgMH2 - deltaMH2, 1D-6)
	MHH2 = avgMH2 + deltaMH2
	Mh0 = -1
	if( Mh02 .ge. 0 ) Mh0 = sqrt(Mh02)
	MHH = -1
	if( MHH2 .ge. 0 ) MHH = sqrt(MHH2)

	CA2 = .5D0*(SB2*MA02 + CB2*MZ2 - Mh02)/deltaMH2
	CA  = sqrt(CA2)
	SA  = -SB*CB*avgMH2/(CA*deltaMH2)
	SA2 = SA**2
	S2A = 2*SA*CA
	C2A = (CA - SA)*(CA + SA)
	SAB = SA*CB + CA*SB
	CAB = CA*CB - SA*SB
	SBA = SB*CA - CB*SA
	SBA2 = SBA**2
	CBA = CB*CA + SB*SA
	CBA2 = CBA**2
	end


************************************************************************

	subroutine CharginoMasses(M_2_, MUE_)
	implicit none
	ComplexType M_2_, MUE_

#include "FH.h"
#include "FHRecord.h"

	ComplexType X(2,2)

	if( debuglevel .ge. 2 ) then
	  DPARA "============================"	ENDL
	  DPARA "CharginoMasses input:"		ENDL
	  DPARA "----------------------------"	ENDL
	  DTAGre(iM2, M_2_)
	  DTAGimz(iM2, M_2_)
	  DTAGre(iMUE, MUE_)
	  DTAGimz(iMUE, MUE_)
	endif

	M_2 = M_2_
	MUE = MUE_

	X(1,1) = M_2
	X(1,2) = sqrt2*MW*SB
	X(2,1) = sqrt2*MW*CB
	X(2,2) = MUE

	call SVD(2, 2, X,2, MCha, UCha,2, VCha,2, 1)

	MCha2(1) = MCha(1)**2
	MCha2(2) = MCha(2)**2

	MUETB(3) = -Conjugate(MUE)/TB
	MUETB(2) = -Conjugate(MUE)*TB
	MUETB(4) = MUETB(2)
	end


************************************************************************

	subroutine NeutralinoMasses(M_1_)
	implicit none
	ComplexType M_1_

#include "FH.h"
#include "FHRecord.h"

	ComplexType Y(4,4)

	if( debuglevel .ge. 2 ) then
	  DPARA "============================"	ENDL
	  DPARA "NeutralinoMasses input:"	ENDL
	  DPARA "----------------------------"	ENDL
	  DTAGrez(iM1, M_1_)
	  DTAGimz(iM1, M_1_)
	endif

	M_1 = M_1_
* We use the SUSY-GUT relation here:
	if( abs(M_1) .lt. 1D-5 ) M_1 = 5/3D0*SW2/CW2*M_2

* Calculation of Neutralino masses and mixing matrix:
	Y(1,1) = M_1
	Y(2,2) = M_2
	Y(3,3) = 0
	Y(4,4) = 0
	Y(1,2) = 0
	Y(1,3) = -MZ*SW*CB
	Y(1,4) = MZ*SW*SB
	Y(2,3) = MZ*CW*CB
	Y(2,4) = -MZ*CW*SB
	Y(3,4) = -MUE

	call TakagiFactor(4, Y,4, MNeu, ZNeu,4, 1)

	MNeu2(1) = MNeu(1)**2
	MNeu2(2) = MNeu(2)**2
	MNeu2(3) = MNeu(3)**2
	MNeu2(4) = MNeu(4)**2
	end


************************************************************************

	subroutine SfermionMasses(t,
     &    M1SL_, M2SL_, M3SL_,
     &    M1SR_, M2SR_, M3SR_,
     &    Af1_, Af2_, Af3_, Q_)
	implicit none
	integer t
	RealType M1SL_, M2SL_, M3SL_
	RealType M1SR_, M2SR_, M3SR_, Q_
	ComplexType Af1_, Af2_, Af3_

#include "FH.h"
#include "FHRecord.h"

	integer nQ, nU

	nQ = tQ(t)
	nU = tU(t)

	if( debuglevel .ge. 2 ) then
	  DPARA "============================"	ENDL
	  DPARA "SfermionMasses input:"		ENDL
	  DPARA "----------------------------"	ENDL
	  DTAG(iMSS(nQ,1), M1SL_)
	  DTAG(iMSS(nQ,2), M2SL_)
	  DTAG(iMSS(nQ,3), M3SL_)
	  DTAG(iMSS(nU,1), M1SR_)
	  DTAG(iMSS(nU,2), M2SR_)
	  DTAG(iMSS(nU,3), M3SR_)
	  DTAGre(iAf(t,1), Af1_)
	  DTAGimz(iAf(t,1), Af1_)
	  DTAGre(iAf(t,2), Af2_)
	  DTAGimz(iAf(t,2), Af2_)
	  DTAGre(iAf(t,3), Af3_)
	  DTAGimz(iAf(t,3), Af3_)
	  DTAGz(iQSf(t), Q_)
	endif

	MSS0(nQ,1) = M1SL_
	MSS(nQ,1) = MSS0(nQ,1)
	MSS0(nQ,2) = M2SL_
	MSS(nQ,2) = MSS0(nQ,2)
	MSS0(nQ,3) = M3SL_
	MSS(nQ,3) = MSS0(nQ,3)

	MSS0(nU,1) = M1SR_
	MSS(nU,1) = MSS0(nU,1)
	MSS0(nU,2) = M2SR_
	MSS(nU,2) = MSS0(nU,2)
	MSS0(nU,3) = M3SR_
	MSS(nU,3) = MSS0(nU,3)

	Af0(t,1) = Af1_
	Af(t,1) = Af0(t,1)
	Af0(t,2) = Af2_
	Af(t,2) = Af0(t,2)
	Af0(t,3) = Af3_
	Af(t,3) = Af0(t,3)

	QSf(t) = Q_
	end


************************************************************************

	subroutine EndPara(error)
	implicit none
	integer error

#include "FH.h"
#include "looptools.h"

	RealType m1, m2, Qcur, Qdefault, ospara(3,3:4)
	integer t, g, fail, drbar, nQ, nU

	RealType Run
	external Run

	character*(*) cplxmsg
	parameter (cplxmsg =
     &    " in the presence of complex parameters")

	do g = 1, 3*3*5
	  MSS2_flat(g) = 0
	enddo
	do g = 1, 3
	  MSS2(g,g,1) = SignSq(MSS(1,g))
	  MSS2(g,g,2) = SignSq(MSS(2,g))
	  MSS2(g,g,3) = SignSq(MSS(3,g))
	  MSS2(g,g,4) = SignSq(MSS(4,g))
	  MSS2(g,g,5) = SignSq(MSS(5,g))

	  do t = 2, 4
	    m1 = Re(MSS2(g,g,tQ(t)))
	    m2 = Re(MSS2(g,g,tU(t)))
	    DSS2(1,t,g) = (m1 - m2) + (DSf(1,t) - DSf(2,t))
	    DSS2(2,t,g) = m2 + DSf(2,t)
	    Xf(t,g) = Af(t,g) + MUETB(t)
	  enddo
	enddo

* DRbar -> OS transition
	Qdefault = sqrt(abs(MSS(3,3)*MSS(4,3)))
	Qcur = -1
	drbar = 0
	fail = 0
	do t = 3, 4
	  if( QSf(t) .lt. 0 ) QSf(t) = Qdefault
	  if( QSf(t) .gt. 0 ) then
	    if( abs(QSf(t) - Qcur) .gt. 1D-10 ) then
	      fail = 0
	      Qcur = QSf(t)
	      call SetDecayScale(fail, Qcur, -1D0, SORT_SF)
	    endif
	    if( fail .eq. 0 ) call DRbartoOS(t, drbar, ospara(1,t))
	    if( .not. btest(drbar, t) )
     &        Error(error, "DRbartoOS transition failed")
	  endif
	enddo

	Ab_bad = 0
	do t = 3, 4
	  if( btest(drbar, t) ) then
	    if( t .eq. 4 .and.
     &          abs(ospara(3,4) - Xf(4,3)) .gt. .3D0*abs(Af(4,3)) ) then
	      Ab_bad = 1
	      if( debuglevel .ge. 2 ) then
	        DPARA "Ab_DRbar     =", Af(4,3) ENDL
	        DPARA "Ab_shift     =", ospara(3,4) - Xf(4,3) ENDL
	      endif
	      Warning("Questionable Ab from DRbar-to-OS conversion")
	    else
	      nQ = tQ(t)
	      nU = tU(t)

	      if( debuglevel .ge. 2 ) then 
	        DPARA "MS"//cMSS(nQ)//"_DRbar    =", MSS(nQ,3) ENDL
	        DPARA "MS"//cMSS(nU)//"_DRbar    =", MSS(nU,3) ENDL
	        DPARA "A"//cAf(t)//"_DRbar     =", Af(t,3) ENDL
	      endif

	      DSS2(1,t,3) = ospara(1,t)
	      DSS2(2,t,3) = ospara(2,t)
	      MSS2(3,3,nQ) = DSS2(1,t,3) + DSS2(2,t,3) - DSf(1,t)
	      MSS(nQ,3) = sqrt(Re(MSS2(3,3,nQ)))
	      MSS2(3,3,nU) = DSS2(2,t,3) - DSf(2,t)
	      MSS(nU,3) = sqrt(Re(MSS2(3,3,nU)))

	      Xf(t,3) = ospara(3,t)
	      Af(t,3) = Xf(t,3) - MUETB(t)

	      if( debuglevel .ge. 2 ) then 
	        DPARA "MS"//cMSS(nQ)//"_OS       =", MSS(nQ,3) ENDL
	        DPARA "MS"//cMSS(nU)//"_OS       =", MSS(nU,3) ENDL
	        DPARA "A"//cAf(t)//"_OS        =", Af(t,3) ENDL
	      endif
	    endif
	  endif
	enddo

	if( abs(Im(Af(2,1))) + abs(Im(Af(3,1))) + abs(Im(Af(4,1))) +
     &      abs(Im(Af(2,2))) + abs(Im(Af(3,2))) + abs(Im(Af(4,2))) +
     &      abs(Im(Af(2,3))) + abs(Im(Af(3,3))) + abs(Im(Af(4,3))) +
     &      abs(Im(MUE)) +
     &      abs(Im(M_1)) + abs(Im(M_2)) + abs(Im(M_3)) .ne. 0 ) then
	  if( higgsmix .lt. 3 )
     &      Warning(Digit(higgsmix)//"x"//Digit(higgsmix)//" mixing"//cplxmsg)
	  if( inputmass .eq. 3 )
     &      Warning("Input parameter MA0"//cplxmsg)
	endif

	if( error .eq. 0 ) para_valid = valid

	if( debuglevel .lt. 1 ) return

	DPARA "============================"	ENDL
	DPARA "parameters set by FHSetPara"	ENDL
	DPARA "----------------------------"	ENDL
	DPARA "mudim =", sqrt(getmudim())	ENDL
	DPARA "MT =", MT			ENDL
	DPARA "MTrun =", Mf(tT,3)		ENDL
	DPARA "1/Alfa(MZ) =", 1/AlfaMZ		ENDL
	DPARA "AlfasMT =", AlfasMT		ENDL
	DPARA "Alfat =", Alfat			ENDL
	DPARA "TB =", TB			ENDL
	DPARA "M_1 =", M_1			ENDL
	DPARA "M_2 =", M_2			ENDL
	DPARA "M_3 =", M_3			ENDL
	DPARA "MUE =", MUE			ENDL
	DPARA "----------------------------"	ENDL
	DPARA "M1SL =", MSS(1,1)		ENDL
	DPARA "M1SE =", MSS(2,1)		ENDL
	DPARA "M1SQ =", MSS(3,1)		ENDL
	DPARA "M1SU =", MSS(4,1)		ENDL
	DPARA "M1SD =", MSS(5,1)		ENDL
	DPARA "Ae   =", Af(2,1)			ENDL
	DPARA "Au   =", Af(3,1)			ENDL
	DPARA "Ad   =", Af(4,1)			ENDL
	DPARA "----------------------------"	ENDL
	DPARA "M2SL =", MSS(1,2)		ENDL
	DPARA "M2SE =", MSS(2,2)		ENDL
	DPARA "M2SQ =", MSS(3,2)		ENDL
	DPARA "M2SU =", MSS(4,2)		ENDL
	DPARA "M2SD =", MSS(5,2)		ENDL
	DPARA "Amu  =", Af(2,2)			ENDL
	DPARA "Ac   =", Af(3,2)			ENDL
	DPARA "As   =", Af(4,2)			ENDL
	DPARA "----------------------------"	ENDL
	DPARA "M3SL =", MSS(1,3)		ENDL
	DPARA "M3SE =", MSS(2,3)		ENDL
	DPARA "M3SQ =", MSS(3,3)		ENDL
	DPARA "M3SU =", MSS(4,3)		ENDL
	DPARA "M3SD =", MSS(5,3)		ENDL
	DPARA "Atau =", Af(2,3)			ENDL
	DPARA "At   =", Af(3,3)			ENDL
	DPARA "Ab   =", Af(4,3)			ENDL
	DPARA "----------------------------"	ENDL
	DPARA "Higgs masses (tree-level):"	ENDL
	DPARA "Mh0 =", Mh0			ENDL
	DPARA "MHH =", MHH			ENDL
	DPARA "MA0 =", MA0			ENDL
	DPARA "MHp =", MHp			ENDL
	DPARA "SA =", SA			ENDL
	DPARA "----------------------------"	ENDL
	DPARA "Chargino masses:"		ENDL
	DPARA "MCha(1) =", MCha(1)		ENDL
	DPARA "MCha(2) =", MCha(2)		ENDL
	DPARA "----------------------------"	ENDL
	DPARA "Neutralino masses:"		ENDL
	DPARA "MNeu(1) =", MNeu(1)		ENDL
	DPARA "MNeu(2) =", MNeu(2)		ENDL
	DPARA "MNeu(3) =", MNeu(3)		ENDL
	DPARA "MNeu(4) =", MNeu(4)		ENDL

	call flush(debugunit)
	end


************************************************************************

	block data recordnamesini

#include "FH.h"
#include "FHRecord.h"

	data cMSS /"L", "E", "Q", "U", "D"/
	data cAf /"n", "e", "u", "d"/

	data FHName(iinvAlfaMZ) /"invAlfaMZ"/
	data FHName(iAlfasMZ) /"AlfasMZ"/
	data FHName(iGF) /"GF"/

	data FHName(iME) /"ME"/
	data FHName(iMU) /"MU"/
	data FHName(iMD) /"MD"/
	data FHName(iMM) /"MM"/
	data FHName(iMC) /"MC"/
	data FHName(iMS) /"MS"/
	data FHName(iML) /"ML"/
	data FHName(iMT) /"MT"/
	data FHName(iMB) /"MB"/

	data FHName(iMW) /"MW"/
	data FHName(iMZ) /"MZ"/

	data FHName(iCKMlambda) /"CKMlambda"/
	data FHName(iCKMA) /"CKMA"/
	data FHName(iCKMrhobar) /"CKMrhobar"/
	data FHName(iCKMetabar) /"CKMetabar"/

	data FHName(iTB) /"TB"/
	data FHName(iMA0) /"MA0"/
	data FHName(iMHp) /"MHp"/

	data FHName(iMSusy) /"MSusy"/
	data FHName(iM1SL) /"M1SL"/
	data FHName(iM2SL) /"M2SL"/
	data FHName(iM3SL) /"M3SL"/
	data FHName(iM1SE) /"M1SE"/
	data FHName(iM2SE) /"M2SE"/
	data FHName(iM3SE) /"M3SE"/
	data FHName(iM1SQ) /"M1SQ"/
	data FHName(iM2SQ) /"M2SQ"/
	data FHName(iM3SQ) /"M3SQ"/
	data FHName(iM1SU) /"M1SU"/
	data FHName(iM2SU) /"M2SU"/
	data FHName(iM3SU) /"M3SU"/
	data FHName(iM1SD) /"M1SD"/
	data FHName(iM2SD) /"M2SD"/
	data FHName(iM3SD) /"M3SD"/
	data FHName(iQtau) /"Qtau"/
	data FHName(iQt) /"Qt"/
	data FHName(iQb) /"Qb"/
	data FHName(iscalefactor) /"scalefactor"/
	data FHName(iprodSqrts) /"prodSqrts"/

	data FHName(iRe(iAe)) /"Re(Ae)"/
	data FHName(iIm(iAe)) /"Im(Ae)"/
	data FHName(iAbs(iAe)) /"Abs(Ae)"/
	data FHName(iArg(iAe)) /"Arg(Ae)"/

	data FHName(iRe(iAu)) /"Re(Au)"/
	data FHName(iIm(iAu)) /"Im(Au)"/
	data FHName(iAbs(iAu)) /"Abs(Au)"/
	data FHName(iArg(iAu)) /"Arg(Au)"/

	data FHName(iRe(iAd)) /"Re(Ad)"/
	data FHName(iIm(iAd)) /"Im(Ad)"/
	data FHName(iAbs(iAd)) /"Abs(Ad)"/
	data FHName(iArg(iAd)) /"Arg(Ad)"/

	data FHName(iRe(iAmu)) /"Re(Amu)"/
	data FHName(iIm(iAmu)) /"Im(Amu)"/
	data FHName(iAbs(iAmu)) /"Abs(Amu)"/
	data FHName(iArg(iAmu)) /"Arg(Amu)"/

	data FHName(iRe(iAc)) /"Re(Ac)"/
	data FHName(iIm(iAc)) /"Im(Ac)"/
	data FHName(iAbs(iAc)) /"Abs(Ac)"/
	data FHName(iArg(iAc)) /"Arg(Ac)"/

	data FHName(iRe(iAs)) /"Re(As)"/
	data FHName(iIm(iAs)) /"Im(As)"/
	data FHName(iAbs(iAs)) /"Abs(As)"/
	data FHName(iArg(iAs)) /"Arg(As)"/

	data FHName(iRe(iAtau)) /"Re(Atau)"/
	data FHName(iIm(iAtau)) /"Im(Atau)"/
	data FHName(iAbs(iAtau)) /"Abs(Atau)"/
	data FHName(iArg(iAtau)) /"Arg(Atau)"/

	data FHName(iRe(iAt)) /"Re(At)"/
	data FHName(iIm(iAt)) /"Im(At)"/
	data FHName(iAbs(iAt)) /"Abs(At)"/
	data FHName(iArg(iAt)) /"Arg(At)"/

	data FHName(iRe(iAb)) /"Re(Ab)"/
	data FHName(iIm(iAb)) /"Im(Ab)"/
	data FHName(iAbs(iAb)) /"Abs(Ab)"/
	data FHName(iArg(iAb)) /"Arg(Ab)"/

	data FHName(iRe(iXtau)) /"Re(Xtau)"/
	data FHName(iIm(iXtau)) /"Im(Xtau)"/
	data FHName(iAbs(iXtau)) /"Abs(Xtau)"/
	data FHName(iArg(iXtau)) /"Arg(Xtau)"/

	data FHName(iRe(iXt)) /"Re(Xt)"/
	data FHName(iIm(iXt)) /"Im(Xt)"/
	data FHName(iAbs(iXt)) /"Abs(Xt)"/
	data FHName(iArg(iXt)) /"Arg(Xt)"/

	data FHName(iRe(iXb)) /"Re(Xb)"/
	data FHName(iIm(iXb)) /"Im(Xb)"/
	data FHName(iAbs(iXb)) /"Abs(Xb)"/
	data FHName(iArg(iXb)) /"Arg(Xb)"/

	data FHName(iRe(iMUE)) /"Re(MUE)"/
	data FHName(iIm(iMUE)) /"Im(MUE)"/
	data FHName(iAbs(iMUE)) /"Abs(MUE)"/
	data FHName(iArg(iMUE)) /"Arg(MUE)"/

	data FHName(iRe(iM1)) /"Re(M_1)"/
	data FHName(iIm(iM1)) /"Im(M_1)"/
	data FHName(iAbs(iM1)) /"Abs(M_1)"/
	data FHName(iArg(iM1)) /"Arg(M_1)"/

	data FHName(iRe(iM2)) /"Re(M_2)"/
	data FHName(iIm(iM2)) /"Im(M_2)"/
	data FHName(iAbs(iM2)) /"Abs(M_2)"/
	data FHName(iArg(iM2)) /"Arg(M_2)"/

	data FHName(iRe(iM3)) /"Re(M_3)"/
	data FHName(iIm(iM3)) /"Im(M_3)"/
	data FHName(iAbs(iM3)) /"Abs(M_3)"/
	data FHName(iArg(iM3)) /"Arg(M_3)"/

	data FHName(iRe(ideltaLLL12)) /"Re(deltaLLL12)"/
	data FHName(iIm(ideltaLLL12)) /"Im(deltaLLL12)"/
	data FHName(iAbs(ideltaLLL12)) /"Abs(deltaLLL12)"/
	data FHName(iArg(ideltaLLL12)) /"Arg(deltaLLL12)"/

	data FHName(iRe(ideltaLLL23)) /"Re(deltaLLL23)"/
	data FHName(iIm(ideltaLLL23)) /"Im(deltaLLL23)"/
	data FHName(iAbs(ideltaLLL23)) /"Abs(deltaLLL23)"/
	data FHName(iArg(ideltaLLL23)) /"Arg(deltaLLL23)"/

	data FHName(iRe(ideltaLLL13)) /"Re(deltaLLL13)"/
	data FHName(iIm(ideltaLLL13)) /"Im(deltaLLL13)"/
	data FHName(iAbs(ideltaLLL13)) /"Abs(deltaLLL13)"/
	data FHName(iArg(ideltaLLL13)) /"Arg(deltaLLL13)"/

	data FHName(iRe(ideltaELR12)) /"Re(deltaELR12)"/
	data FHName(iIm(ideltaELR12)) /"Im(deltaELR12)"/
	data FHName(iAbs(ideltaELR12)) /"Abs(deltaELR12)"/
	data FHName(iArg(ideltaELR12)) /"Arg(deltaELR12)"/

	data FHName(iRe(ideltaELR23)) /"Re(deltaELR23)"/
	data FHName(iIm(ideltaELR23)) /"Im(deltaELR23)"/
	data FHName(iAbs(ideltaELR23)) /"Abs(deltaELR23)"/
	data FHName(iArg(ideltaELR23)) /"Arg(deltaELR23)"/

	data FHName(iRe(ideltaELR13)) /"Re(deltaELR13)"/
	data FHName(iIm(ideltaELR13)) /"Im(deltaELR13)"/
	data FHName(iAbs(ideltaELR13)) /"Abs(deltaELR13)"/
	data FHName(iArg(ideltaELR13)) /"Arg(deltaELR13)"/

	data FHName(iRe(ideltaERL12)) /"Re(deltaERL12)"/
	data FHName(iIm(ideltaERL12)) /"Im(deltaERL12)"/
	data FHName(iAbs(ideltaERL12)) /"Abs(deltaERL12)"/
	data FHName(iArg(ideltaERL12)) /"Arg(deltaERL12)"/

	data FHName(iRe(ideltaERL23)) /"Re(deltaERL23)"/
	data FHName(iIm(ideltaERL23)) /"Im(deltaERL23)"/
	data FHName(iAbs(ideltaERL23)) /"Abs(deltaERL23)"/
	data FHName(iArg(ideltaERL23)) /"Arg(deltaERL23)"/

	data FHName(iRe(ideltaERL13)) /"Re(deltaERL13)"/
	data FHName(iIm(ideltaERL13)) /"Im(deltaERL13)"/
	data FHName(iAbs(ideltaERL13)) /"Abs(deltaERL13)"/
	data FHName(iArg(ideltaERL13)) /"Arg(deltaERL13)"/

	data FHName(iRe(ideltaERR12)) /"Re(deltaERR12)"/
	data FHName(iIm(ideltaERR12)) /"Im(deltaERR12)"/
	data FHName(iAbs(ideltaERR12)) /"Abs(deltaERR12)"/
	data FHName(iArg(ideltaERR12)) /"Arg(deltaERR12)"/

	data FHName(iRe(ideltaERR23)) /"Re(deltaERR23)"/
	data FHName(iIm(ideltaERR23)) /"Im(deltaERR23)"/
	data FHName(iAbs(ideltaERR23)) /"Abs(deltaERR23)"/
	data FHName(iArg(ideltaERR23)) /"Arg(deltaERR23)"/

	data FHName(iRe(ideltaERR13)) /"Re(deltaERR13)"/
	data FHName(iIm(ideltaERR13)) /"Im(deltaERR13)"/
	data FHName(iAbs(ideltaERR13)) /"Abs(deltaERR13)"/
	data FHName(iArg(ideltaERR13)) /"Arg(deltaERR13)"/

	data FHName(iRe(ideltaQLL12)) /"Re(deltaQLL12)"/
	data FHName(iIm(ideltaQLL12)) /"Im(deltaQLL12)"/
	data FHName(iAbs(ideltaQLL12)) /"Abs(deltaQLL12)"/
	data FHName(iArg(ideltaQLL12)) /"Arg(deltaQLL12)"/

	data FHName(iRe(ideltaQLL23)) /"Re(deltaQLL23)"/
	data FHName(iIm(ideltaQLL23)) /"Im(deltaQLL23)"/
	data FHName(iAbs(ideltaQLL23)) /"Abs(deltaQLL23)"/
	data FHName(iArg(ideltaQLL23)) /"Arg(deltaQLL23)"/

	data FHName(iRe(ideltaQLL13)) /"Re(deltaQLL13)"/
	data FHName(iIm(ideltaQLL13)) /"Im(deltaQLL13)"/
	data FHName(iAbs(ideltaQLL13)) /"Abs(deltaQLL13)"/
	data FHName(iArg(ideltaQLL13)) /"Arg(deltaQLL13)"/

	data FHName(iRe(ideltaULR12)) /"Re(deltaULR12)"/
	data FHName(iIm(ideltaULR12)) /"Im(deltaULR12)"/
	data FHName(iAbs(ideltaULR12)) /"Abs(deltaULR12)"/
	data FHName(iArg(ideltaULR12)) /"Arg(deltaULR12)"/

	data FHName(iRe(ideltaULR23)) /"Re(deltaULR23)"/
	data FHName(iIm(ideltaULR23)) /"Im(deltaULR23)"/
	data FHName(iAbs(ideltaULR23)) /"Abs(deltaULR23)"/
	data FHName(iArg(ideltaULR23)) /"Arg(deltaULR23)"/

	data FHName(iRe(ideltaULR13)) /"Re(deltaULR13)"/
	data FHName(iIm(ideltaULR13)) /"Im(deltaULR13)"/
	data FHName(iAbs(ideltaULR13)) /"Abs(deltaULR13)"/
	data FHName(iArg(ideltaULR13)) /"Arg(deltaULR13)"/

	data FHName(iRe(ideltaURL12)) /"Re(deltaURL12)"/
	data FHName(iIm(ideltaURL12)) /"Im(deltaURL12)"/
	data FHName(iAbs(ideltaURL12)) /"Abs(deltaURL12)"/
	data FHName(iArg(ideltaURL12)) /"Arg(deltaURL12)"/

	data FHName(iRe(ideltaURL23)) /"Re(deltaURL23)"/
	data FHName(iIm(ideltaURL23)) /"Im(deltaURL23)"/
	data FHName(iAbs(ideltaURL23)) /"Abs(deltaURL23)"/
	data FHName(iArg(ideltaURL23)) /"Arg(deltaURL23)"/

	data FHName(iRe(ideltaURL13)) /"Re(deltaURL13)"/
	data FHName(iIm(ideltaURL13)) /"Im(deltaURL13)"/
	data FHName(iAbs(ideltaURL13)) /"Abs(deltaURL13)"/
	data FHName(iArg(ideltaURL13)) /"Arg(deltaURL13)"/

	data FHName(iRe(ideltaURR12)) /"Re(deltaURR12)"/
	data FHName(iIm(ideltaURR12)) /"Im(deltaURR12)"/
	data FHName(iAbs(ideltaURR12)) /"Abs(deltaURR12)"/
	data FHName(iArg(ideltaURR12)) /"Arg(deltaURR12)"/

	data FHName(iRe(ideltaURR23)) /"Re(deltaURR23)"/
	data FHName(iIm(ideltaURR23)) /"Im(deltaURR23)"/
	data FHName(iAbs(ideltaURR23)) /"Abs(deltaURR23)"/
	data FHName(iArg(ideltaURR23)) /"Arg(deltaURR23)"/

	data FHName(iRe(ideltaURR13)) /"Re(deltaURR13)"/
	data FHName(iIm(ideltaURR13)) /"Im(deltaURR13)"/
	data FHName(iAbs(ideltaURR13)) /"Abs(deltaURR13)"/
	data FHName(iArg(ideltaURR13)) /"Arg(deltaURR13)"/

	data FHName(iRe(ideltaDLR12)) /"Re(deltaDLR12)"/
	data FHName(iIm(ideltaDLR12)) /"Im(deltaDLR12)"/
	data FHName(iAbs(ideltaDLR12)) /"Abs(deltaDLR12)"/
	data FHName(iArg(ideltaDLR12)) /"Arg(deltaDLR12)"/

	data FHName(iRe(ideltaDLR23)) /"Re(deltaDLR23)"/
	data FHName(iIm(ideltaDLR23)) /"Im(deltaDLR23)"/
	data FHName(iAbs(ideltaDLR23)) /"Abs(deltaDLR23)"/
	data FHName(iArg(ideltaDLR23)) /"Arg(deltaDLR23)"/

	data FHName(iRe(ideltaDLR13)) /"Re(deltaDLR13)"/
	data FHName(iIm(ideltaDLR13)) /"Im(deltaDLR13)"/
	data FHName(iAbs(ideltaDLR13)) /"Abs(deltaDLR13)"/
	data FHName(iArg(ideltaDLR13)) /"Arg(deltaDLR13)"/

	data FHName(iRe(ideltaDRL12)) /"Re(deltaDRL12)"/
	data FHName(iIm(ideltaDRL12)) /"Im(deltaDRL12)"/
	data FHName(iAbs(ideltaDRL12)) /"Abs(deltaDRL12)"/
	data FHName(iArg(ideltaDRL12)) /"Arg(deltaDRL12)"/

	data FHName(iRe(ideltaDRL23)) /"Re(deltaDRL23)"/
	data FHName(iIm(ideltaDRL23)) /"Im(deltaDRL23)"/
	data FHName(iAbs(ideltaDRL23)) /"Abs(deltaDRL23)"/
	data FHName(iArg(ideltaDRL23)) /"Arg(deltaDRL23)"/

	data FHName(iRe(ideltaDRL13)) /"Re(deltaDRL13)"/
	data FHName(iIm(ideltaDRL13)) /"Im(deltaDRL13)"/
	data FHName(iAbs(ideltaDRL13)) /"Abs(deltaDRL13)"/
	data FHName(iArg(ideltaDRL13)) /"Arg(deltaDRL13)"/

	data FHName(iRe(ideltaDRR12)) /"Re(deltaDRR12)"/
	data FHName(iIm(ideltaDRR12)) /"Im(deltaDRR12)"/
	data FHName(iAbs(ideltaDRR12)) /"Abs(deltaDRR12)"/
	data FHName(iArg(ideltaDRR12)) /"Arg(deltaDRR12)"/

	data FHName(iRe(ideltaDRR23)) /"Re(deltaDRR23)"/
	data FHName(iIm(ideltaDRR23)) /"Im(deltaDRR23)"/
	data FHName(iAbs(ideltaDRR23)) /"Abs(deltaDRR23)"/
	data FHName(iArg(ideltaDRR23)) /"Arg(deltaDRR23)"/

	data FHName(iRe(ideltaDRR13)) /"Re(deltaDRR13)"/
	data FHName(iIm(ideltaDRR13)) /"Im(deltaDRR13)"/
	data FHName(iAbs(ideltaDRR13)) /"Abs(deltaDRR13)"/
	data FHName(iArg(ideltaDRR13)) /"Arg(deltaDRR13)"/
	end

