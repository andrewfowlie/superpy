:Evaluate: BeginPackage["FeynHiggs`"]

:Evaluate: FHSetFlags::usage =
	"FHSetFlags sets the FeynHiggs flags."

:Evaluate: FHSetFlagsString::usage =
	"FHSetFlagsString sets the FeynHiggs flags."

:Evaluate: FHRetrieveFlags::usage =
	"FHRetrieveFlags retrieves the FeynHiggs flags."

:Evaluate: FHRetrieveFlagsString::usage =
	"FHRetrieveFlagsString retrieves the FeynHiggs flags."

:Evaluate: FHSetSMPara::usage =
	"FHSetSMPara sets the FeynHiggs SM parameters."

:Evaluate: FHRetrieveSMPara::usage =
	"FHRetrieveSMPara retrieves the FeynHiggs SM parameters."

:Evaluate: FHGetSMPara::usage =
	"FHGetSMPara returns the parameters computed by FHSetSMPara."

:Evaluate: FHSetPara::usage =
	"FHSetPara sets the FeynHiggs input parameters."

:Evaluate: FHRetrievePara::usage =
	"FHRetrievePara retrieves the FeynHiggs input parameters."

:Evaluate: FHRetrieveOSPara::usage =
	"FHRetrieveOSPara retrieves the on-shell FeynHiggs input parameters."

:Evaluate: FHSetSLHA::usage =
	"FHSetSLHA sets the FeynHiggs parameters from an SLHA file."

:Evaluate: FHSetLFV::usage =
	"FHSetLFV sets the lepton-flavour-violating parameters."

:Evaluate: FHSetNMFV::usage =
	"FHSetNMFV sets the non-minimal flavour-violating parameters."

:Evaluate: FHRetrieveLFV::usage =
	"FHRetrieveLFV retrieves the lepton-flavour-violating parameters."

:Evaluate: FHRetrieveNMFV::usage =
	"FHRetrieveNMFV retrieves the non-minimal flavour-violating parameters."

:Evaluate: FHGetPara::usage =
	"FHGetPara returns the parameters computed by FHSetPara."

:Evaluate: FHGetTLPara::usage =
	"FHGetTLPara returns parameters used in the internal computation of the neutral Higgs masses in FeynHiggs."

:Evaluate: FHGetNMFV::usage =
	"FHGetNMFV returns the parameters computed by FHSetNMFV."

:Evaluate: FHHiggsCorr::usage =
	"FHHiggsCorr computes the Higgs masses and mixings."

:Evaluate: FHUncertainties::usage =
	"FHUncertainties computes error estimates for the Higgs masses and mixings."

:Evaluate: FHCouplings::usage =
	"FHCouplings computes the Higgs couplings, widths, and branching ratios."

:Evaluate: FHConstraints::usage =
	"FHConstraints evaluates electroweak precision observables as further constraints on the MSSM parameter space."

:Evaluate: FHFlavour::usage =
	"FHFlavour evaluates flavour observables as further constraints on the MSSM parameter space."

:Evaluate: FHHiggsProd::usage =
	"FHHiggsProd computes (approximate) Higgs production cross-sections."

:Evaluate: FHGetSelf::usage =
	"FHGetSelf computes various Higgs self-energies plus their derivatives."

:Evaluate: FHAddSelf::usage =
	"FHAddSelf registers user-defined shifts for the Higgs self-energies."

:Evaluate: FHOutput::usage =
	"FHOutput writes the FeynHiggs input and outputs to a file."

:Evaluate: FHOutputSLHA::usage =
	"FHOutputSLHA writes the FeynHiggs input and outputs to an SLHA file."

:Evaluate: FHRecord::usage =
	"FHRecord contains the parameters of a FeynHiggs Record."

:Evaluate: FHRecordIndex::usage =
	"FHRecordIndex looks up the FeynHiggs Record index of a parameter name."

:Evaluate: FHClearRecord::usage =
	"FHClearRecord returns an empty FeynHiggs Record."

:Evaluate: FHReadRecord::usage =
	"FHReadRecord reads a FeynHiggs Record from a file."

:Evaluate: FHSLHARecord::usage =
	"FHSLHARecord reads a FeynHiggs Record from an SLHA file."

:Evaluate: FHLoopRecord::usage =
	"FHLoopRecord advances the parameter values in a FeynHiggs Record and returns the new record if the loop continues and False if it stops."

:Evaluate: FHSetRecord::usage =
	"FHSetRecord sets the input parameters from a FeynHiggs Record."

:Evaluate: FHRetrieveRecord::usage =
	"FHRetrieveRecord fills a FeynHiggs Record from the parameters currently set."

:Evaluate: FHLoadTable::usage =
	"FHLoadTable loads a parameter table into internal storage."

:Evaluate: FHTableRecord::usage =
	"FHTableRecord associates a FeynHiggs Record with the internal table."

:Evaluate: FHSetDebug::usage =
	"FHSetDebug sets the FeynHiggs debug level."

:Evaluate: FHSelectUZ::usage =
	"FHSelectUZ chooses which of UHiggs (= 1) or ZHiggs (= 2) to use for internal and external Higgs bosons, i.e. in the couplings and the decays, respectively, and whether to use resummed masses in the couplings."

:Evaluate: FHError::usage =
	"FHError is an error message returned by FeynHiggs."

:Evaluate: FHAbort[f_Symbol] := (Message[f::badsyntax]; Abort[])

:Evaluate: General::badsyntax =
	"Probably not all arguments have numerical values."

:Evaluate: MapIndexed[(Key[#] = 2^(#2[[1]] - 1))&,
	SelfID = {h0h0, HHHH, A0A0, HmHp,
	  h0HH, h0A0, HHA0,
	  G0G0, h0G0, HHG0, A0G0,
	  GmGp, HmGp}]

:Evaluate: Module[ {offset = 1, indexdef},
	Attributes[indexdef] = {HoldAll, Listable};
	indexdef[stride_, i_] :=
	  (i =.; ToString[i] -> (i = (offset += stride) - stride));
	FHRecordIndices = Flatten[{
	  indexdef[1, {iVar, iLower, iUpper, iStep}],
	  offset = 1;
	  indexdef[1, iAdmin],
	  indexdef[0, FHRecordR],
	  indexdef[1, {iinvAlfaMZ, iAlfasMZ, iGF,
	    iME, iMU, iMD,
	    iMM, iMC, iMS,
	    iML, iMT, iMB,
	    iMW, iMZ,
	    iCKMlambda, iCKMA, iCKMrhobar, iCKMetabar,
	    iTB, iMA0, iMHp,
	    iMSusy,
	    iM1SL, iM1SE, iM1SQ, iM1SU, iM1SD,
	    iM2SL, iM2SE, iM2SQ, iM2SU, iM2SD,
	    iM3SL, iM3SE, iM3SQ, iM3SU, iM3SD,
	    iQtau, iQt, iQb, iscalefactor, iprodSqrts}],
	  indexdef[0, FHRecordC],
	  indexdef[4, {iAe, iAu, iAd,
	    iAmu, iAc, iAs,
	    iAtau, iAt, iAb,
	    iXtau, iXt, iXb,
	    iMUE, iM1, iM2, iM3,
	    ideltaLLL12, ideltaLLL23, ideltaLLL13,
	    ideltaELR12, ideltaELR23, ideltaELR13,
	    ideltaERL12, ideltaERL23, ideltaERL13,
	    ideltaERR12, ideltaERR23, ideltaERR13,
	    ideltaQLL12, ideltaQLL23, ideltaQLL13,
	    ideltaULR12, ideltaULR23, ideltaULR13,
	    ideltaURL12, ideltaURL23, ideltaURL13,
	    ideltaURR12, ideltaURR23, ideltaURR13,
	    ideltaDLR12, ideltaDLR23, ideltaDLR13,
	    ideltaDRL12, ideltaDRL23, ideltaDRL13,
	    ideltaDRR12, ideltaDRR23, ideltaDRR13}],
	  indexdef[-1, FHRecordE],
	  indexdef[0, FHRecordN] }] ];
	iRe[v_] := v;
	iIm[v_] := v + 1;
	iAbs[v_] := v + 2;
	iArg[v_] := v + 3;
	FHWriteIndex[] := WriteString["RecordIndices.h",
	  "#ifndef RECORDINDICES_H\n" <>
	  "#define RECORDINDICES_H\n\n" <>
	  ({"#define ", #1, " ", ToString[#2], "\n"}&@@@ FHRecordIndices) <>
	  "\n#endif\n"]

:Evaluate: Begin["`Private`"]

:Begin:
:Function: mFHSetFlags
:Pattern:
  FHSetFlags[mssmpart_, fieldren_, tanbren_,
    higgsmix_, p2approx_, looplevel_,
    runningMT_, botResum_, tlCplxApprox_]
:Arguments: {
  mssmpart, fieldren, tanbren,
  higgsmix, p2approx, looplevel,
  runningMT, botResum, tlCplxApprox }
:ArgumentTypes: {
  Integer, Integer, Integer,
  Integer, Integer, Integer,
  Integer, Integer, Integer }
:ReturnType: Manual
:End:

:Evaluate: FHSetFlags[s_String] := FHSetFlagsString[s]

:Evaluate: _FHSetFlags := FHAbort[FHSetFlags]

:Begin:
:Function: mFHSetFlagsString
:Pattern: FHSetFlagsString[flags_]
:Arguments: {flags}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Evaluate: _FHSetFlagsString := FHAbort[FHSetFlagsString]

:Begin:
:Function: mFHRetrieveFlags
:Pattern: FHRetrieveFlags[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Evaluate: _FHRetrieveFlags := FHAbort[FHRetrieveFlags]

:Begin:
:Function: mFHRetrieveFlagsString
:Pattern: FHRetrieveFlagsString[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Evaluate: _FHRetrieveFlagsString := FHAbort[FHRetrieveFlagsString]

:Begin:
:Function: mFHSetSMPara
:Pattern:
  FHSetSMPara[invAlfa_, AlfasMZ_, GF_,
    ME_, MM_, ML_, MU_, MC_, MD_, MS_, MB_,
    MW_, MZ_,
    CKMlambda_, CKMA_, CKMrhobar_, CKMetabar_]
:Arguments: {
  N[invAlfa], N[AlfasMZ], N[GF],
  N[ME], N[MM], N[ML], N[MU], N[MC], N[MD], N[MS], N[MB],
  N[MW], N[MZ],
  N[CKMlambda], N[CKMA], N[CKMrhobar], N[CKMetabar]}
:ArgumentTypes: {
  Real, Real, Real,
  Real, Real, Real, Real, Real, Real, Real, Real,
  Real, Real,
  Real, Real, Real, Real }
:ReturnType: Manual
:End:

:Evaluate: _FHSetSMPara := FHAbort[FHSetSMPara]

:Begin:
:Function: mFHRetrieveSMPara
:Pattern: FHRetrieveSMPara[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Evaluate: _FHRetrieveSMPara := FHAbort[FHRetrieveSMPara]

:Begin:
:Function: mFHGetSMPara
:Pattern: FHGetSMPara[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Evaluate: _FHGetSMPara := FHAbort[FHGetSMPara]

:Begin:
:Function: mFHSetPara
:Pattern:
  FHSetPara[scalefactor_,
    MT_, TB_, MA0_, MHp_,
    M3SL_, M3SE_, M3SQ_, M3SU_, M3SD_,
    M2SL_, M2SE_, M2SQ_, M2SU_, M2SD_,
    M1SL_, M1SE_, M1SQ_, M1SU_, M1SD_,
    MUE_,
    Atau_, At_, Ab_,
    Amu_, Ac_, As_,
    Ae_, Au_, Ad_,
    M1_, M2_, M3_,
    Qtau_, Qt_, Qb_]
:Arguments: {
  N[scalefactor],
  N[MT], N[TB], N[MA0], N[MHp], 
  N[M3SL], N[M3SE], N[M3SQ], N[M3SU], N[M3SD],
  N[M2SL], N[M2SE], N[M2SQ], N[M2SU], N[M2SD],
  N[M1SL], N[M1SE], N[M1SQ], N[M1SU], N[M1SD],
  N[Re[MUE]], N[Im[MUE]], 
  N[Re[Atau]], N[Im[Atau]], N[Re[At]], N[Im[At]], N[Re[Ab]], N[Im[Ab]],
  N[Re[Amu]], N[Im[Amu]], N[Re[Ac]], N[Im[Ac]], N[Re[As]], N[Im[As]],
  N[Re[Ae]], N[Im[Ae]], N[Re[Au]], N[Im[Au]], N[Re[Ad]], N[Im[Ad]],
  N[Re[M1]], N[Im[M1]], N[Re[M2]], N[Im[M2]], N[Re[M3]], N[Im[M3]], 
  N[Qtau], N[Qt], N[Qb] }
:ArgumentTypes: {
  Real,
  Real, Real, Real, Real,
  Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real,
  Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real }
:ReturnType: Manual
:End:

:Begin:
:Function: mFHRetrievePara
:Pattern: FHRetrievePara[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHRetrieveOSPara
:Pattern: FHRetrieveOSPara[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHSetSLHA
:Pattern: FHSetSLHA[file_]
:Arguments: {file}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHSetLFV
:Pattern: FHSetLFV[
  deltaLLL12_, deltaLLL23_, deltaLLL13_,
  deltaELR12_, deltaELR23_, deltaELR13_,
  deltaERL12_, deltaERL23_, deltaERL13_,
  deltaERR12_, deltaERR23_, deltaERR13_ ]
:Arguments: {
  N[Re[deltaLLL12]], N[Im[deltaLLL12]],
  N[Re[deltaLLL23]], N[Im[deltaLLL23]],
  N[Re[deltaLLL13]], N[Im[deltaLLL13]],
  N[Re[deltaELR12]], N[Im[deltaELR12]],
  N[Re[deltaELR23]], N[Im[deltaELR23]],
  N[Re[deltaELR13]], N[Im[deltaELR13]],
  N[Re[deltaERL12]], N[Im[deltaERL12]],
  N[Re[deltaERL23]], N[Im[deltaERL23]],
  N[Re[deltaERL13]], N[Im[deltaERL13]],
  N[Re[deltaERR12]], N[Im[deltaERR12]],
  N[Re[deltaERR23]], N[Im[deltaERR23]],
  N[Re[deltaERR13]], N[Im[deltaERR13]] }
:ArgumentTypes: {
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHSetNMFV
:Pattern: FHSetNMFV[
  deltaQLL12_, deltaQLL23_, deltaQLL13_,
  deltaULR12_, deltaULR23_, deltaULR13_,
  deltaURL12_, deltaURL23_, deltaURL13_,
  deltaURR12_, deltaURR23_, deltaURR13_,
  deltaDLR12_, deltaDLR23_, deltaDLR13_,
  deltaDRL12_, deltaDRL23_, deltaDRL13_,
  deltaDRR12_, deltaDRR23_, deltaDRR13_ ]
:Arguments: {
  N[Re[deltaQLL12]], N[Im[deltaQLL12]],
  N[Re[deltaQLL23]], N[Im[deltaQLL23]],
  N[Re[deltaQLL13]], N[Im[deltaQLL13]],
  N[Re[deltaULR12]], N[Im[deltaULR12]],
  N[Re[deltaULR23]], N[Im[deltaULR23]],
  N[Re[deltaULR13]], N[Im[deltaULR13]],
  N[Re[deltaURL12]], N[Im[deltaURL12]],
  N[Re[deltaURL23]], N[Im[deltaURL23]],
  N[Re[deltaURL13]], N[Im[deltaURL13]],
  N[Re[deltaURR12]], N[Im[deltaURR12]],
  N[Re[deltaURR23]], N[Im[deltaURR23]],
  N[Re[deltaURR13]], N[Im[deltaURR13]],
  N[Re[deltaDLR12]], N[Im[deltaDLR12]],
  N[Re[deltaDLR23]], N[Im[deltaDLR23]],
  N[Re[deltaDLR13]], N[Im[deltaDLR13]],
  N[Re[deltaDRL12]], N[Im[deltaDRL12]],
  N[Re[deltaDRL23]], N[Im[deltaDRL23]],
  N[Re[deltaDRL13]], N[Im[deltaDRL13]],
  N[Re[deltaDRR12]], N[Im[deltaDRR12]],
  N[Re[deltaDRR23]], N[Im[deltaDRR23]],
  N[Re[deltaDRR13]], N[Im[deltaDRR13]] }
:ArgumentTypes: {
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real,
  Real, Real, Real, Real, Real, Real}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHRetrieveLFV
:Pattern: FHRetrieveLFV[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHRetrieveNMFV
:Pattern: FHRetrieveNMFV[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHGetPara
:Pattern: FHGetPara[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHGetTLPara
:Pattern: FHGetTLPara[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHGetNMFV
:Pattern: FHGetNMFV[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHHiggsCorr
:Pattern: FHHiggsCorr[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHUncertainties
:Pattern: FHUncertainties[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHCouplings
:Pattern: FHCouplings[fast_:1]
:Arguments: {fast}
:ArgumentTypes: {Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHConstraints
:Pattern: FHConstraints[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHFlavour
:Pattern: FHFlavour[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHHiggsProd
:Pattern: FHHiggsProd[sqrts_]
:Arguments: {sqrts}
:ArgumentTypes: {Real}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHGetSelf
:Pattern: FHGetSelf[p2_, key_, dkey_]
:Arguments: {N[p2], key, dkey}
:ArgumentTypes: {Real, Integer, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHAddSelf
:Pattern: FHAddSelf[sig_List, rotate_]
:Arguments: {Flatten[Transpose[{Re[sig], Im[sig]}]], rotate}
:ArgumentTypes: {RealList, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHOutput
:Pattern: FHOutput[file_, key_, sqrts_:0]
:Arguments: {file, key, sqrts}
:ArgumentTypes: {String, Integer, Real}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHOutputSLHA
:Pattern: FHOutputSLHA[file_, key_]
:Arguments: {file, key}
:ArgumentTypes: {String, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHClearRecord
:Pattern: FHClearRecord[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHLoopRecord
:Pattern: FHLoopRecord[FHRecord[para__List]]
:Arguments: {N[Flatten[Transpose[{para}]]]}
:ArgumentTypes: {RealList}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHSetRecord
:Pattern: FHSetRecord[FHRecord[para__List]]
:Arguments: {N[Flatten[Transpose[{para}]]]}
:ArgumentTypes: {RealList}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHRetrieveRecord
:Pattern: FHRetrieveRecord[FHRecord[para__List], iX_]
:Arguments: {N[Flatten[Transpose[{para}]]], iX}
:ArgumentTypes: {RealList, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHRecordIndex
:Pattern: FHRecordIndex[para_]
:Arguments: {para}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHReadRecord
:Pattern: FHReadRecord[file_]
:Arguments: {file}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHSLHARecord
:Pattern: FHSLHARecord[file_]
:Arguments: {file}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHLoadTable
:Pattern: FHLoadTable[file_]
:Arguments: {file}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHTableRecord
:Pattern: FHTableRecord[FHRecord[para__List], i1_, i2_]
:Arguments: {N[Flatten[Transpose[{para}]]], i1, i2}
:ArgumentTypes: {RealList, Integer, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHSetDebug
:Pattern: FHSetDebug[debuglevel_]
:Arguments: {debuglevel}
:ArgumentTypes: {Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: mFHSelectUZ
:Pattern: FHSelectUZ[uzint_, uzext_, mfeff_]
:Arguments: {uzint, uzext, mfeff}
:ArgumentTypes: {Integer, Integer, Integer}
:ReturnType: Manual
:End:

:Evaluate: CTensor[a_, dims_] :=
  RTensor[Apply[Complex, Partition[a, 2], 1], dims]

:Evaluate: RTensor[a_, dims_] :=
  Transpose[Fold[Partition, Chop[a], dims][[1]], Range[Length[dims], 1, -1]]

:Evaluate: Format[_FHRecord] := "-FHRecord-"

:Evaluate: ToRecord[para_List] :=
	FHRecord@@ Transpose[Partition[para, FHRecordN]]

:Evaluate: Cases[ LinkPatterns[$CurrentLink],
	_[s_[___]] :> (_s := FHAbort[s]) ]

:Evaluate: End[]

:Evaluate: EndPackage[]



/*
	MFeynHiggs.tm
		the Mathematica frontend for FeynHiggs
		this file is part of FeynHiggs
		last modified 7 Oct 13 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <pthread.h>

#include "mathlink.h"
#ifndef MLCONST
#define MLCONST
#endif

#include "CFeynHiggs.h"
#include "CSLHA.h"

typedef MLCONST char cchar;
typedef unsigned char byte;
typedef const int cint;

#if QUAD
#define MLPutReal MLPutReal128
#define MLPutRealList MLPutReal128List
typedef const int len_t;
#else
typedef const long len_t;
#endif

#define _c_ ,
#define _s_ ;

#define _Vr_(v,d) v
#define _Vc_(v,d) ToComplex2(re_##v, im_##v)
#define _Va_(f) f(_c_, _Vr_, _Vc_, _Vr_, _Vc_)

#define _Rx_(v,d) &v
#define _Rax_(v,d) v
#define _Ra_(f) f(_c_, _Rx_, _Rx_, _Rax_, _Rax_)

#define _Mr_(v,d) cRealType v
#define _Mc_(v,d) cRealType re_##v, cRealType im_##v
#define _Mar_(v,d) cRealType v d
#define _Mac_(v,d) cRealType re_##v d, cRealType im_##v d
#define _Ma_(f) f(_c_, _Mr_, _Mc_, _Mar_, _Mac_)

#define _Lr_(v,d) RealType v
#define _Lar_(v,d) RealType v d
#define _Lc_(v,d) ComplexType v
#define _Lac_(v,d) ComplexType v d
#define _La_(f) f(_s_, _Lr_, _Lc_, _Lar_, _Lac_)

extern void FORTRAN(fortranflush)();

/******************************************************************/

static int forcestderr = 0;
static int stdoutorig;
static int stdoutpipe[2];
static pthread_t stdouttid;

static void *MLstdout(void *fd)
{
  static byte *buf = NULL;
  static long size = 0;
  enum { unit = 10240 };
  long len = 0, n = 0;

  do {
    len += n;
    if( size - len < 128 ) buf = realloc(buf, size += unit);
    n = read(*(int *)fd, buf + len, size - len);
  } while( n > 0 );

  if( len ) {
    MLPutFunction(stdlink, "EvaluatePacket", 1);
    MLPutFunction(stdlink, "WriteString", 2);
    MLPutString(stdlink, "stdout");
    MLPutByteString(stdlink, buf, len);
    MLEndPacket(stdlink);

    MLNextPacket(stdlink);
    MLNewPacket(stdlink);
  }

  return NULL;
}

/******************************************************************/

static inline void BeginRedirect()
{
  if( forcestderr == 0 &&
      pipe(stdoutpipe) == 0 &&
      pthread_create(&stdouttid, NULL, MLstdout, stdoutpipe) == 0 ) {
    dup2(stdoutpipe[1], 1);
    close(stdoutpipe[1]);
  }
  else dup2(2, 1);
}

/******************************************************************/

static void EndRedirect()
{
  void *ret;

  FORTRAN(fortranflush)();
  fflush(stdout);
  dup2(stdoutorig, 1);
  if( stdouttid ) pthread_join(stdouttid, &ret);
}

/******************************************************************/

static void MLPutStatus(MLINK mlp, int error)
{
  if( error ) {
    MLPutFunction(mlp, "FHError", 1);
    MLPutInteger(mlp, error);
  }
  else MLPutSymbol(mlp, "True");
}

/******************************************************************/

/*#define Context "FeynHiggs`"*/
#define Context

#define MLPutFHSymbol(mlp,s) \
  MLPutSymbol(mlp, Context s)

#define MLPutRule(mlp,s) \
  MLPutFunction(mlp, "Rule", 2); \
  MLPutFHSymbol(mlp, #s)

#define MLPutRules(mlp,s,n) \
  MLPutFunction(mlp, "Rule", 2); \
  MLPutFunction(mlp, Context #s, n)

#define MLPutIRule(mlp,v) \
  MLPutRule(mlp, v); \
  MLPutInteger(mlp, v)

#define MLPutRRule(mlp,v) \
  MLPutRule(mlp, v); \
  MLPutReal(mlp, v)

#define MLPutCRule(mlp,v) \
  MLPutRule(mlp, v); \
  MLPutComplex(mlp, v)

#define MLPutRLRule(mlp,v,n) \
  MLPutRule(mlp, v); \
  MLPutRealList(mlp, v, n)

/******************************************************************/

static void MLPutComplex(MLINK mlp, cComplexType c)
{
  if( Im(c) == 0 ) MLPutReal(mlp, Re(c));
  else {
    MLPutFunction(mlp, "Complex", 2);
    MLPutReal(mlp, Re(c));
    MLPutReal(mlp, Im(c));
  }
}

/******************************************************************/

static void MLPutRealTensor(MLINK mlp, RealType *a, len_t len,
  int *dims, len_t depth)
{
  MLPutFunction(mlp, "FeynHiggs`Private`RTensor", 2);
  MLPutRealList(mlp, a, len);
  MLPutIntegerList(mlp, dims, depth);
}

/******************************************************************/

static void MLPutComplexTensor(MLINK mlp, ComplexType *a, cint len,
  int *dims, cint depth)
{
  MLPutFunction(mlp, "FeynHiggs`Private`CTensor", 2);
  MLPutRealList(mlp, (RealType *)a, 2*len);
  MLPutIntegerList(mlp, dims, depth);
}

/******************************************************************/

static void mFHSetFlags(cint mssmpart, cint fieldren, cint tanbren,
  cint higgsmix, cint p2approx, cint looplevel,
  cint runningMT, cint botResum, cint tlCplxApprox)
{
  int error;

  BeginRedirect();

  FHSetFlags(&error, mssmpart, fieldren, tanbren,
    higgsmix, p2approx, looplevel,
    runningMT, botResum, tlCplxApprox);

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetFlagsString(cchar *flags)
{
  int error;

  BeginRedirect();

  FHSetFlagsString(&error, flags);

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrieveFlags(void)
{
  int error, mssmpart, fieldren, tanbren;
  int higgsmix, p2approx, looplevel;
  int runningMT, botResum, tlCplxApprox;

  BeginRedirect();

  FHRetrieveFlags(&error, &mssmpart, &fieldren, &tanbren,
    &higgsmix, &p2approx, &looplevel,
    &runningMT, &botResum, &tlCplxApprox);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 9);

    MLPutIRule(stdlink, mssmpart);
    MLPutIRule(stdlink, fieldren);
    MLPutIRule(stdlink, tanbren);
    MLPutIRule(stdlink, higgsmix);
    MLPutIRule(stdlink, p2approx);
    MLPutIRule(stdlink, looplevel);
    MLPutIRule(stdlink, runningMT);
    MLPutIRule(stdlink, botResum);
    MLPutIRule(stdlink, tlCplxApprox);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrieveFlagsString(void)
{
  int error;
  char flags[10];

  BeginRedirect();

  FHRetrieveFlagsString(&error, flags);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else MLPutByteString(stdlink, (unsigned char *)flags, 9);

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetSMPara(_Ma_(argsSetSMPara))
{
  int error;

  BeginRedirect();

  FHSetSMPara(&error, _Va_(argsSetSMPara));

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrieveSMPara(void)
{
  int error;
  _La_(argsSetSMPara);

  BeginRedirect();

  FHRetrieveSMPara(&error, _Ra_(argsSetSMPara));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 17);

    MLPutRRule(stdlink, invAlfa);
    MLPutRRule(stdlink, AlfasMZ);
    MLPutRRule(stdlink, GF);

    MLPutRRule(stdlink, ME);
    MLPutRRule(stdlink, MM);
    MLPutRRule(stdlink, ML);
    MLPutRRule(stdlink, MU);
    MLPutRRule(stdlink, MC);
    MLPutRRule(stdlink, MD);
    MLPutRRule(stdlink, MS);
    MLPutRRule(stdlink, MB);

    MLPutRRule(stdlink, MW);
    MLPutRRule(stdlink, MZ);

    MLPutRRule(stdlink, CKMlambda);
    MLPutRRule(stdlink, CKMA);
    MLPutRRule(stdlink, CKMrhobar);
    MLPutRRule(stdlink, CKMetabar);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHGetSMPara(void)
{
  int error;
  _La_(argsGetSMPara);

  BeginRedirect();

  FHGetSMPara(&error, _Ra_(argsGetSMPara));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 1);

    MLPutRule(stdlink, CKM);
    MLPutComplexTensor(stdlink, (ComplexType *)CKM, 3*3, (int[]){3, 3}, 2);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetPara(_Ma_(argsOSPara), _Ma_(argsQPara))
{
  int error;

  BeginRedirect();

  FHSetPara(&error, _Va_(argsOSPara), _Va_(argsQPara));

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrievePara(void)
{
  int error;
  _La_(argsOSPara);
  _La_(argsQPara);

  BeginRedirect();

  FHRetrievePara(&error, _Ra_(argsOSPara), _Ra_(argsQPara));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 36);

    MLPutRRule(stdlink, scalefactor);

    MLPutRRule(stdlink, MT);
    MLPutRRule(stdlink, TB);
    MLPutRRule(stdlink, MA0);
    MLPutRRule(stdlink, MHp);

    MLPutRRule(stdlink, M3SL);
    MLPutRRule(stdlink, M3SE);
    MLPutRRule(stdlink, M3SQ);
    MLPutRRule(stdlink, M3SU);
    MLPutRRule(stdlink, M3SD);

    MLPutRRule(stdlink, M2SL);
    MLPutRRule(stdlink, M2SE);
    MLPutRRule(stdlink, M2SQ);
    MLPutRRule(stdlink, M2SU);
    MLPutRRule(stdlink, M2SD);

    MLPutRRule(stdlink, M1SL);
    MLPutRRule(stdlink, M1SE);
    MLPutRRule(stdlink, M1SQ);
    MLPutRRule(stdlink, M1SU);
    MLPutRRule(stdlink, M1SD);

    MLPutCRule(stdlink, MUE);

    MLPutCRule(stdlink, Atau);
    MLPutCRule(stdlink, At);
    MLPutCRule(stdlink, Ab);

    MLPutCRule(stdlink, Amu);
    MLPutCRule(stdlink, Ac);
    MLPutCRule(stdlink, As);

    MLPutCRule(stdlink, Ae);
    MLPutCRule(stdlink, Au);
    MLPutCRule(stdlink, Ad);

    MLPutCRule(stdlink, M1);
    MLPutCRule(stdlink, M2);
    MLPutCRule(stdlink, M3);

    MLPutRRule(stdlink, Qtau);
    MLPutRRule(stdlink, Qt);
    MLPutRRule(stdlink, Qb);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrieveOSPara(void)
{
  int error;
  _La_(argsOSPara);

  BeginRedirect();

  FHRetrieveOSPara(&error, _Ra_(argsOSPara));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 33);

    MLPutRRule(stdlink, scalefactor);

    MLPutRRule(stdlink, MT);
    MLPutRRule(stdlink, TB);
    MLPutRRule(stdlink, MA0);
    MLPutRRule(stdlink, MHp);

    MLPutRRule(stdlink, M3SL);
    MLPutRRule(stdlink, M3SE);
    MLPutRRule(stdlink, M3SQ);
    MLPutRRule(stdlink, M3SU);
    MLPutRRule(stdlink, M3SD);

    MLPutRRule(stdlink, M2SL);
    MLPutRRule(stdlink, M2SE);
    MLPutRRule(stdlink, M2SQ);
    MLPutRRule(stdlink, M2SU);
    MLPutRRule(stdlink, M2SD);

    MLPutRRule(stdlink, M1SL);
    MLPutRRule(stdlink, M1SE);
    MLPutRRule(stdlink, M1SQ);
    MLPutRRule(stdlink, M1SU);
    MLPutRRule(stdlink, M1SD);

    MLPutCRule(stdlink, MUE);

    MLPutCRule(stdlink, Atau);
    MLPutCRule(stdlink, At);
    MLPutCRule(stdlink, Ab);

    MLPutCRule(stdlink, Amu);
    MLPutCRule(stdlink, Ac);
    MLPutCRule(stdlink, As);

    MLPutCRule(stdlink, Ae);
    MLPutCRule(stdlink, Au);
    MLPutCRule(stdlink, Ad);

    MLPutCRule(stdlink, M3);
    MLPutCRule(stdlink, M2);
    MLPutCRule(stdlink, M1);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetSLHA(cchar *file)
{
  int error;
  COMPLEX slhadata[nslhadata];

  BeginRedirect();

  SLHARead(&error, slhadata, file, 0);
  if( error == 0 ) FHSetSLHA(&error, slhadata);

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetLFV(_Ma_(argsSetLFV))
{
  int error;

  BeginRedirect();

  FHSetLFV(&error, _Va_(argsSetLFV));

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrieveLFV(void)
{
  int error;
  _La_(argsSetLFV);

  BeginRedirect();

  FHRetrieveLFV(&error, _Ra_(argsSetLFV));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 12);

    MLPutCRule(stdlink, deltaLLL12);
    MLPutCRule(stdlink, deltaLLL23);
    MLPutCRule(stdlink, deltaLLL13);

    MLPutCRule(stdlink, deltaELR12);
    MLPutCRule(stdlink, deltaELR23);
    MLPutCRule(stdlink, deltaELR13);

    MLPutCRule(stdlink, deltaERL12);
    MLPutCRule(stdlink, deltaERL23);
    MLPutCRule(stdlink, deltaERL13);

    MLPutCRule(stdlink, deltaERR12);
    MLPutCRule(stdlink, deltaERR23);
    MLPutCRule(stdlink, deltaERR13);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetNMFV(_Ma_(argsSetNMFV))
{
  int error;

  BeginRedirect();

  FHSetNMFV(&error, _Va_(argsSetNMFV));

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrieveNMFV(void)
{
  int error;
  _La_(argsSetNMFV);

  BeginRedirect();

  FHRetrieveNMFV(&error, _Ra_(argsSetNMFV));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 21);

    MLPutCRule(stdlink, deltaQLL12);
    MLPutCRule(stdlink, deltaQLL23);
    MLPutCRule(stdlink, deltaQLL13);

    MLPutCRule(stdlink, deltaULR12);
    MLPutCRule(stdlink, deltaULR23);
    MLPutCRule(stdlink, deltaULR13);

    MLPutCRule(stdlink, deltaURL12);
    MLPutCRule(stdlink, deltaURL23);
    MLPutCRule(stdlink, deltaURL13);

    MLPutCRule(stdlink, deltaURR12);
    MLPutCRule(stdlink, deltaURR23);
    MLPutCRule(stdlink, deltaURR13);

    MLPutCRule(stdlink, deltaDLR12);
    MLPutCRule(stdlink, deltaDLR23);
    MLPutCRule(stdlink, deltaDLR13);

    MLPutCRule(stdlink, deltaDRL12);
    MLPutCRule(stdlink, deltaDRL23);
    MLPutCRule(stdlink, deltaDRL13);

    MLPutCRule(stdlink, deltaDRR12);
    MLPutCRule(stdlink, deltaDRR23);
    MLPutCRule(stdlink, deltaDRR13);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHGetPara(void)
{
  int error, nmfv, t, g;
  _La_(argsGetPara);

  BeginRedirect();

  FHGetPara(&error, &nmfv, _Ra_(argsGetPara));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 5*(3*2 + 2) + 9);

    for( t = 0; t < 5; ++t ) {
      for( g = 0; g < 3; ++g ) {
        MLPutRules(stdlink, MSf, 2);
        MLPutInteger(stdlink, t + 1);
        MLPutInteger(stdlink, g + 1);
        MLPutRealList(stdlink, MSf[g][t], 2);

        MLPutRules(stdlink, USf, 2);
        MLPutInteger(stdlink, t + 1);
        MLPutInteger(stdlink, g + 1);
        MLPutComplexTensor(stdlink, (ComplexType *)USf[g][t], 2*2, (int[]){2, 2}, 2);
      }

      MLPutRules(stdlink, MASf, 1);
      MLPutInteger(stdlink, t + 1);
      MLPutRealList(stdlink, MASf[t], 6);

      MLPutRules(stdlink, UASf, 1);
      MLPutInteger(stdlink, t + 1);
      MLPutComplexTensor(stdlink, (ComplexType *)UASf[t],
        6*6, (int[]){6, 6}, 2);
    }

    MLPutRLRule(stdlink, MCha, 2);

    MLPutRule(stdlink, UCha);
    MLPutComplexTensor(stdlink, (ComplexType *)UCha, 2*2, (int[]){2, 2}, 2);

    MLPutRule(stdlink, VCha);
    MLPutComplexTensor(stdlink, (ComplexType *)VCha, 2*2, (int[]){2, 2}, 2);

    MLPutRLRule(stdlink, MNeu, 4);

    MLPutRule(stdlink, ZNeu);
    MLPutComplexTensor(stdlink, (ComplexType *)ZNeu, 4*4, (int[]){4, 4}, 2);

    MLPutCRule(stdlink, Deltab);

    MLPutRRule(stdlink, MGl);

    MLPutRLRule(stdlink, MHtree, 4);
    MLPutRRule(stdlink, SAtree);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHGetTLPara(void)
{
  int error;
  _La_(argsGetTLPara);

  BeginRedirect();

  FHGetTLPara(&error, _Ra_(argsGetTLPara));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 4);

    MLPutRule(stdlink, MSb);
    MLPutRealList(stdlink, MSb, 2);

    MLPutRule(stdlink, USb);
    MLPutComplexTensor(stdlink, (ComplexType *)USb, 2*2, (int[]){2, 2}, 2);

    MLPutRule(stdlink, MbSL2);
    MLPutReal(stdlink, MbSL2);

    MLPutRule(stdlink, Deltab);
    MLPutComplex(stdlink, Deltab);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHGetNMFV(void)
{
  int error, t, n;
  _La_(argsGetNMFV);

  BeginRedirect();

  FHGetNMFV(&error, _Ra_(argsGetNMFV));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 8);

    for( n = 0; n < 5; ++n ) {
      MLPutRules(stdlink, MSS2, 1);
      MLPutInteger(stdlink, n + 1);
      MLPutComplexTensor(stdlink, (ComplexType *)MSS2[n],
        3*3, (int[]){3, 3}, 2);
    }

    for( t = 0; t < 3; ++t ) {
      MLPutRules(stdlink, Kf, 1);
      MLPutInteger(stdlink, t + 2);
      MLPutComplexTensor(stdlink, (ComplexType *)Kf[t],
        3*3, (int[]){3, 3}, 2);
    }
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHHiggsCorr(void)
{
  int error;
  _La_(argsHiggsCorr);

  BeginRedirect();

  FHHiggsCorr(&error, _Ra_(argsHiggsCorr));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 4);

    MLPutRLRule(stdlink, MHiggs, 4);

    MLPutCRule(stdlink, SAeff);

    MLPutRule(stdlink, UHiggs);
    MLPutComplexTensor(stdlink, (ComplexType *)UHiggs, 3*3, (int[]){3, 3}, 2);

    MLPutRule(stdlink, ZHiggs);
    MLPutComplexTensor(stdlink, (ComplexType *)ZHiggs, 3*3, (int[]){3, 3}, 2);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHUncertainties(void)
{
  int error;
  _La_(argsUncertainties);

  BeginRedirect();

  FHUncertainties(&error, _Ra_(argsUncertainties));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 4);

    MLPutRLRule(stdlink, DeltaMHiggs, 4);

    MLPutCRule(stdlink, DeltaSAeff);

    MLPutRule(stdlink, DeltaUHiggs);
    MLPutComplexTensor(stdlink, (ComplexType *)DeltaUHiggs, 3*3, (int[]){3, 3}, 2);

    MLPutRule(stdlink, DeltaZHiggs);
    MLPutComplexTensor(stdlink, (ComplexType *)DeltaZHiggs, 3*3, (int[]){3, 3}, 2);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHCouplings(cint fast)
{
  int error;
  _La_(argsCouplings);

  BeginRedirect();

  FHCouplings(&error, _Ra_(argsCouplings), fast);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 49);

#define MLPutLHS(array, channel) \
  MLPutRules(stdlink, array, 1); \
  MLPutSymbol(stdlink, Context #channel)

#define MLPutArray1(array, channel, i) \
  MLPutLHS(array, channel); \
  MLPutTensor(stdlink, &array(channel(1)), i, (int[]){i}, 1)

#define MLPutArray2(array, channel, i,j) \
  MLPutLHS(array, channel); \
  MLPutTensor(stdlink, &array(channel(1,1)), i*j, (int[]){i,j}, 2)

#define MLPutArray3(array, channel, i,j,k) \
  MLPutLHS(array, channel); \
  MLPutTensor(stdlink, &array(channel(1,1,1)), i*j*k, (int[]){i,j,k}, 3)

#define MLPutArray4(array, channel, i,j,k,l) \
  MLPutLHS(array, channel); \
  MLPutTensor(stdlink, &array(channel(1,1,1,1)), i*j*k*l, (int[]){i,j,k,l}, 4)

#define MLPutArray5(array, channel, i,j,k,l,m) \
  MLPutLHS(array, channel); \
  MLPutTensor(stdlink, &array(channel(1,1,1,1,1)), i*j*k*l*m, (int[]){i,j,k,l,m}, 5)

/* COUPLINGS */

#define MLPutTensor MLPutComplexTensor

    MLPutArray2(Coupling,  H0VV,     3,5);
    MLPutArray4(LCoupling, H0FF,     3,4,3,3);
    MLPutArray4(RCoupling, H0FF,     3,4,3,3);
    MLPutArray3(LCoupling, HpFF,     2,3,3);
    MLPutArray3(RCoupling, HpFF,     2,3,3);
    MLPutArray3(LCoupling, H0ChaCha, 3,2,2);
    MLPutArray3(RCoupling, H0ChaCha, 3,2,2);
    MLPutArray3(LCoupling, H0NeuNeu, 3,4,4);
    MLPutArray3(RCoupling, H0NeuNeu, 3,4,4);
    MLPutArray2(LCoupling, HpNeuCha, 4,3);
    MLPutArray2(RCoupling, HpNeuCha, 4,3);
    MLPutArray2(Coupling,  H0HV,     3,3);
    MLPutArray1(Coupling,  HpHV,     3);
    MLPutArray3(Coupling,  H0HH,     3,4,4);
    MLPutArray5(Coupling,  H0SfSf,   3,2,2,4,3);
    MLPutArray5(Coupling,  HpSfSf,   2,2,2,3,3);

    MLPutArray2(CouplingSM,  H0VV,   3,5);
    MLPutArray4(LCouplingSM, H0FF,   3,4,3,3);
    MLPutArray4(RCouplingSM, H0FF,   3,4,3,3);

#undef MLPutTensor

/* DECAY WIDTHS */

#define MLPutTensor MLPutRealTensor

    MLPutRule(stdlink, GammaTot);
    MLPutRealList(stdlink, &GammaTot(1), 4);

    MLPutArray2(Gamma, H0VV,     3,5);
    MLPutArray4(Gamma, H0FF,     3,4,3,3);
    MLPutArray3(Gamma, HpFF,     2,3,3);
    MLPutArray3(Gamma, H0ChaCha, 3,2,2);
    MLPutArray3(Gamma, H0NeuNeu, 3,4,4);
    MLPutArray2(Gamma, HpNeuCha, 4,3);
    MLPutArray2(Gamma, H0HV,     3,3);
    MLPutArray1(Gamma, HpHV,     3);
    MLPutArray3(Gamma, H0HH,     3,4,4);
    MLPutArray5(Gamma, H0SfSf,   3,2,2,4,3);
    MLPutArray5(Gamma, HpSfSf,   2,2,2,3,3);
    MLPutArray1(Gamma, tBF,      2);

    MLPutArray2(BR, H0VV,     3,5);
    MLPutArray4(BR, H0FF,     3,4,3,3);
    MLPutArray3(BR, HpFF,     2,3,3);
    MLPutArray3(BR, H0ChaCha, 3,2,2);
    MLPutArray3(BR, H0NeuNeu, 3,4,4);
    MLPutArray2(BR, HpNeuCha, 4,3);
    MLPutArray2(BR, H0HV,     3,3);
    MLPutArray1(BR, HpHV,     3);
    MLPutArray3(BR, H0HH,     3,4,4);
    MLPutArray5(BR, H0SfSf,   3,2,2,4,3);
    MLPutArray5(BR, HpSfSf,   2,2,2,3,3);
    MLPutArray1(BR, tBF,      2);

    MLPutRule(stdlink, GammaSMTot);
    MLPutRealList(stdlink, &GammaSMTot(1), 3);

    MLPutArray2(GammaSM, H0VV, 3,5);
    MLPutArray4(GammaSM, H0FF, 3,4,3,3);

    MLPutArray2(BRSM,    H0VV, 3,5);
    MLPutArray4(BRSM,    H0FF, 3,4,3,3);

#undef MLPutTensor
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHConstraints(void)
{
  int error, ccb;
  _La_(argsConstraints);

  BeginRedirect();

  FHConstraints(&error, _Ra_(argsConstraints), &ccb);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 9);

    MLPutRRule(stdlink, gm2);
    MLPutRRule(stdlink, DeltaRho);
    MLPutRRule(stdlink, MWMSSM);
    MLPutRRule(stdlink, MWSM);
    MLPutRRule(stdlink, SW2MSSM);
    MLPutRRule(stdlink, SW2SM);
    MLPutRRule(stdlink, EDMeTh);
    MLPutRRule(stdlink, EDMn);
    MLPutRRule(stdlink, EDMHg);
    MLPutIRule(stdlink, ccb);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHFlavour(void)
{
  int error;
  _La_(argsFlavour);

  BeginRedirect();

  FHFlavour(&error, _Ra_(argsFlavour));

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 4);

    MLPutRRule(stdlink, BsgammaMSSM);
    MLPutRRule(stdlink, BsgammaSM);
    MLPutRRule(stdlink, DeltaMsMSSM);
    MLPutRRule(stdlink, DeltaMsSM);
    MLPutRRule(stdlink, BsmumuMSSM);
    MLPutRRule(stdlink, BsmumuSM);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHHiggsProd(cRealType sqrts)
{
  int error;
  RealType prodxs[nprodxs];

  BeginRedirect();

  FHHiggsProd(&error, sqrts, prodxs);

  EndRedirect();

#define ProdXS(channel) \
  MLPutRule(stdlink, channel); \
  MLPutRealList(stdlink, &channel(1), 3)

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "List", 16);

    ProdXS(bbh);
    ProdXS(bbhSM);
    ProdXS(btagbh);
    ProdXS(btagbhSM);
    ProdXS(ggh);
    ProdXS(gghSM);
    ProdXS(qqh);
    ProdXS(qqhSM);
    ProdXS(tth);
    ProdXS(tthSM);
    ProdXS(Wh);
    ProdXS(WhSM);
    ProdXS(Zh);
    ProdXS(ZhSM);
    ProdXS(StSth);

/* cannot use MLPutRule here because tHmLHC expands immediately */
    MLPutFunction(stdlink, "Rule", 2);
    MLPutFHSymbol(stdlink, "tHm");
    MLPutReal(stdlink, tHm);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHGetSelf(cRealType p2, cint key, cint dkey)
{
  int error, i;
  ComplexType sig[nsig], dsig[nsig];

  BeginRedirect();

  FHGetSelf(&error, p2, key, sig, dkey, dsig);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    int n = 0;
    for( i = 0; i < nsig; ++i )
      n += ((key >> i) & 1) + ((dkey >> i) & 1);
    MLPutFunction(stdlink, "List", n);

    for( i = 0; i < nsig; ++i ) {
      if( (key >> i) & 1 ) {
        MLPutRules(stdlink, Sigma, 1);
        MLPutFunction(stdlink, "Part", 2);
        MLPutFHSymbol(stdlink, "SelfID");
        MLPutInteger(stdlink, i + 1);
        MLPutComplex(stdlink, sig[i]);
      }
      if( (dkey >> i) & 1 ) {
        MLPutRules(stdlink, DSigma, 1);
        MLPutFunction(stdlink, "Part", 2);
        MLPutFHSymbol(stdlink, "SelfID");
        MLPutInteger(stdlink, i + 1);
        MLPutComplex(stdlink, dsig[i]);
      }
    }
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHAddSelf(RealType *sig, len_t sig_len, cint rotate)
{
  int error = 999;

  if( sig_len == 2*nsig ) {
    BeginRedirect();
    FHAddSelf(&error, (cComplexType *)sig, rotate);
    EndRedirect();
  }

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHOutput(cchar *file, cint key, cRealType sqrts)
{
  int error;

  BeginRedirect();

  FHOutput(&error, file, key, sqrts);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else MLPutString(stdlink, file);

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHOutputSLHA(cchar *file, cint key)
{
  int error;
  COMPLEX slhadata[nslhadata];

  BeginRedirect();

  SLHAClear(slhadata);
  FHOutputSLHA(&error, slhadata, key);
  if( error == 0 ) SLHAWrite(&error, slhadata, file);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else MLPutString(stdlink, file);

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRecordIndex(cchar *para)
{
  int ind;

  BeginRedirect();

  FHRecordIndex(&ind, para);

  EndRedirect();

  MLPutInteger(stdlink, ind);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHClearRecord(void)
{
  RealType record[nrecord];

  BeginRedirect();

  FHClearRecord(record);

  EndRedirect();

  MLPutFunction(stdlink, "FeynHiggs`Private`ToRecord", 1);
  MLPutRealList(stdlink, record, nrecord);
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHReadRecord(cchar *file)
{
  int error;
  RealType record[nrecord];
  COMPLEX slhadata[nslhadata];

  BeginRedirect();

  FHReadRecord(&error, record, slhadata, file);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "FeynHiggs`Private`ToRecord", 1);
    MLPutRealList(stdlink, record, nrecord);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSLHARecord(cchar *file)
{
  int error;
  RealType record[nrecord];
  COMPLEX slhadata[nslhadata];

  BeginRedirect();

  SLHARead(&error, slhadata, file, 0);
  if( error == 0 ) FHSLHARecord(&error, record, slhadata);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "FeynHiggs`Private`ToRecord", 1);
    MLPutRealList(stdlink, record, nrecord);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHLoopRecord(RealType *record, len_t record_len)
{
  int error = 999;

  if( record_len == nrecord ) {
    BeginRedirect();
    FHLoopRecord(&error, record);
    EndRedirect();
  }

  if( error > 0 ) MLPutStatus(stdlink, error);
  else if( error < 0 ) MLPutSymbol(stdlink, "False");
  else {
    MLPutFunction(stdlink, "FeynHiggs`Private`ToRecord", 1);
    MLPutRealList(stdlink, record, nrecord);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetRecord(RealType *record, len_t record_len)
{
  int error = 999;

  if( record_len == nrecord ) {
    BeginRedirect();
    FHSetRecord(&error, record);
    EndRedirect();
  }

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "FeynHiggs`Private`ToRecord", 1);
    MLPutRealList(stdlink, record, nrecord);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHRetrieveRecord(RealType *record, len_t record_len, cint iX)
{
  int error = 999;

  if( record_len == nrecord ) {
    BeginRedirect();
    FHRetrieveRecord(&error, record, iX);
    EndRedirect();
  }

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "FeynHiggs`Private`ToRecord", 1);
    MLPutRealList(stdlink, record, nrecord);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHLoadTable(cchar *file)
{
  int error;

  BeginRedirect();

  FHLoadTable(&error, file);

  EndRedirect();

  if( error ) MLPutStatus(stdlink, error);
  else MLPutSymbol(stdlink, "True");

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHTableRecord(RealType *record, len_t record_len,
  cint i1, cint i2)
{
  int error = 999;

  if( record_len == nrecord ) {
    BeginRedirect();
    FHTableRecord(&error, record, i1, i2);
    EndRedirect();
  }

  if( error ) MLPutStatus(stdlink, error);
  else {
    MLPutFunction(stdlink, "FeynHiggs`Private`ToRecord", 1);
    MLPutRealList(stdlink, record, nrecord);
  }

  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSetDebug(cint debuglevel)
{
  BeginRedirect();

  FHSetDebug(debuglevel);

  EndRedirect();

  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}

/******************************************************************/

static void mFHSelectUZ(cint uzint, cint uzext, cint mfeff)
{
  int error;

  BeginRedirect();

  FHSelectUZ(&error, uzint, uzext, mfeff);

  EndRedirect();

  MLPutStatus(stdlink, error);
  MLEndPacket(stdlink);
}

/******************************************************************/

int main(int argc, char **argv)
{
  int fd;

	/* make sure a pipe will not overlap with 0, 1, 2 */
  do { fd = open("/dev/null", O_WRONLY); } while( fd <= 2 );
  close(fd);

  if( getenv("FHFORCESTDERR") ) forcestderr = 1;
  stdoutorig = dup(1);

  return MLMain(argc, argv);
}

