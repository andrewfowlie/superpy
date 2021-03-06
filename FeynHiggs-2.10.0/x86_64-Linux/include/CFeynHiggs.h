#ifndef FTYPES_H
#define FTYPES_H

#if 0
#define FORTRAN(s) s
#else
#define FORTRAN(s) s##_
#endif

#if QUAD

#define RealType long double

#pragma pack(push, 1)
typedef union {
  long double r10;
  struct {
    unsigned long long frac;
    unsigned short exp;
  } i10;
  struct {
    char zero[6];
    unsigned long long frac;
    unsigned short exp;
  } i16;
  unsigned long long i8[2];
unsigned char b[16];
} REAL;
#pragma pack(pop)

static inline REAL ToREAL(const RealType r) {
  REAL new;
  new.i8[0] = 0;
  new.i16.frac = ((REAL *)&r)->i10.frac << 1;
  new.i16.exp = ((REAL *)&r)->i10.exp;
  return new;
}

static inline RealType ToReal(const REAL r) {
  REAL new;
  const long long z = r.i16.frac | (r.i16.exp & 0x7fff);
  new.i10.frac = (r.i16.frac >> 1) | ((z | -z) & 0x8000000000000000LL);
  new.i10.exp = r.i16.exp;
  return new.r10;
}

static inline void ToRealArray(RealType *out, const REAL *in, const int n) {
  int i;
  for( i = 0; i < n; ++i ) out[i] = ToReal(in[i]);
}

static inline void ToREALArray(REAL *out, const RealType *in, const int n) {
  int i;
  for( i = 0; i < n; ++i ) out[i] = ToREAL(in[i]);
}

#else

#define RealType double
typedef double REAL;

#define ToReal(r) (r)
#define ToREAL(r) (r)

#endif

typedef int INTEGER;
typedef const INTEGER CINTEGER;
typedef const REAL CREAL;
typedef struct { REAL re, im; } COMPLEX;
typedef const COMPLEX CCOMPLEX;
typedef char CHARACTER;
typedef const CHARACTER CCHARACTER;

#ifdef __cplusplus

#include <complex>
typedef std::complex<RealType> ComplexType;
#define ToComplex(c) ComplexType(ToReal((c).re), ToReal((c).im))
#define ToComplex2(r,i) ComplexType(r, i)
#define Re(x) std::real(x)
#define Im(x) std::imag(x)

#elif __STDC_VERSION__ >= 199901L

#include <complex.h>
typedef RealType complex ComplexType;
#define ToComplex(c) (ToReal((c).re) + I*ToReal((c).im))
#define ToComplex2(r,i) (r + I*(i))
#define Re(x) creal(x)
#define Im(x) cimag(x)

#else

typedef struct { RealType re, im; } ComplexType;
#define ToComplex(c) (ComplexType){ToReal((c).re), ToReal((c).im)}
#define ToComplex2(r,i) (ComplexType){r, i}
#define Re(x) (x).re
#define Im(x) (x).im

#endif

typedef const RealType cRealType;
typedef const ComplexType cComplexType;

#endif

#ifndef RECORDINDICES_H
#define RECORDINDICES_H

#define iVar 1
#define iLower 2
#define iUpper 3
#define iStep 4
#define iAdmin 1
#define FHRecordR 2
#define iinvAlfaMZ 2
#define iAlfasMZ 3
#define iGF 4
#define iME 5
#define iMU 6
#define iMD 7
#define iMM 8
#define iMC 9
#define iMS 10
#define iML 11
#define iMT 12
#define iMB 13
#define iMW 14
#define iMZ 15
#define iCKMlambda 16
#define iCKMA 17
#define iCKMrhobar 18
#define iCKMetabar 19
#define iTB 20
#define iMA0 21
#define iMHp 22
#define iMSusy 23
#define iM1SL 24
#define iM1SE 25
#define iM1SQ 26
#define iM1SU 27
#define iM1SD 28
#define iM2SL 29
#define iM2SE 30
#define iM2SQ 31
#define iM2SU 32
#define iM2SD 33
#define iM3SL 34
#define iM3SE 35
#define iM3SQ 36
#define iM3SU 37
#define iM3SD 38
#define iQtau 39
#define iQt 40
#define iQb 41
#define iscalefactor 42
#define iprodSqrts 43
#define FHRecordC 44
#define iAe 44
#define iAu 48
#define iAd 52
#define iAmu 56
#define iAc 60
#define iAs 64
#define iAtau 68
#define iAt 72
#define iAb 76
#define iXtau 80
#define iXt 84
#define iXb 88
#define iMUE 92
#define iM1 96
#define iM2 100
#define iM3 104
#define ideltaLLL12 108
#define ideltaLLL23 112
#define ideltaLLL13 116
#define ideltaELR12 120
#define ideltaELR23 124
#define ideltaELR13 128
#define ideltaERL12 132
#define ideltaERL23 136
#define ideltaERL13 140
#define ideltaERR12 144
#define ideltaERR23 148
#define ideltaERR13 152
#define ideltaQLL12 156
#define ideltaQLL23 160
#define ideltaQLL13 164
#define ideltaULR12 168
#define ideltaULR23 172
#define ideltaULR13 176
#define ideltaURL12 180
#define ideltaURL23 184
#define ideltaURL13 188
#define ideltaURR12 192
#define ideltaURR23 196
#define ideltaURR13 200
#define ideltaDLR12 204
#define ideltaDLR23 208
#define ideltaDLR13 212
#define ideltaDRL12 216
#define ideltaDRL23 220
#define ideltaDRL13 224
#define ideltaDRR12 228
#define ideltaDRR23 232
#define ideltaDRR13 236
#define FHRecordE 240
#define FHRecordN 239

#endif
/*
	CFeynHiggs.h
		C/C++ wrapper functions for the FH subroutines
		this file is part of FeynHiggs
		last modified 31 Jul 13 th
*/


#ifndef CFEYNHIGGS_H
#define CFEYNHIGGS_H

#include <string.h>
#include "FHCouplings.h"

enum { nsig = 13 };
enum { nrecord = 4*FHRecordN };

#define couplingS(c) couplings[c-1]
#define couplingsmS(c) couplingsms[c-1]
#define gammaS(c) gammas[c-1]
#define gammasmS(c) gammasms[c-1]
#define prodXS(c) prodxs[c-1]

#define _fh_c_ ,
#define _fh_s_ ;

#define _fh_Fcr_(v,d) CREAL *v
#define _fh_Fcc_(v,d) CCOMPLEX *v
#define _fh_Fca_(f) f(_fh_c_, _fh_Fcr_, _fh_Fcc_, _fh_Fcr_, _fh_Fcc_)

#define _fh_Fr_(v,d) REAL *v
#define _fh_Fc_(v,d) COMPLEX *v
#define _fh_Fa_(f) f(_fh_c_, _fh_Fr_, _fh_Fc_, _fh_Fr_, _fh_Fc_)

#define _fh_Ccr_(v,d) cRealType v
#define _fh_Ccar_(v,d) cRealType v d
#define _fh_Ccc_(v,d) cComplexType v
#define _fh_Ccac_(v,d) cComplexType v d
#define _fh_Cca_(f) f(_fh_c_, _fh_Ccr_, _fh_Ccc_, _fh_Ccar_, _fh_Ccac_)

#define _fh_Cr_(v,d) RealType *v
#define _fh_Car_(v,d) RealType v d
#define _fh_Cc_(v,d) ComplexType *v
#define _fh_Cac_(v,d) ComplexType v d
#define _fh_Ca_(f) f(_fh_c_, _fh_Cr_, _fh_Cc_, _fh_Car_, _fh_Cac_)

#if QUAD

#define _fh_Pcr_(v,d) (CREAL *)v##_
#define _fh_Pcc_(v,d) (CCOMPLEX *)v##_
#define _fh_Pr_(v,d) (REAL *)v##_
#define _fh_Pc_(v,d) (COMPLEX *)v##_
#define _fh_Pcac_(v,d) (CCOMPLEX *)v##_
#define _fh_Dr_(v,d) REAL v##_ d
#define _fh_Dc_(v,d) COMPLEX v##_ d
#define _fh_Ir_(v,d) ToREALArray((REAL *)v##_, (cRealType *)&v, sizeof v##_/sizeof(REAL))
#define _fh_Ic_(v,d) ToREALArray((REAL *)v##_, (cRealType *)&v, sizeof v##_/sizeof(REAL))
#define _fh_Or_(v,d) ToRealArray((RealType *)v, (CREAL *)v##_, sizeof v##_/sizeof(REAL))
#define _fh_Oc_(v,d) ToRealArray((RealType *)v, (CREAL *)v##_, sizeof v##_/sizeof(REAL))

#else

#define _fh_Pcr_(v,d) (CREAL *)&v
#define _fh_Pcc_(v,d) (CCOMPLEX *)&v
#define _fh_Pr_(v,d) (REAL *)v
#define _fh_Pc_(v,d) (COMPLEX *)v
#define _fh_Pcac_(v,d) (CCOMPLEX *)v
#define _fh_Dr_(v,d)
#define _fh_Dc_(v,d)
#define _fh_Or_(v,d)
#define _fh_Oc_(v,d)
#define _fh_Ir_(v,d)
#define _fh_Ic_(v,d)

#endif

#define _fh_Pca_(f) f(_fh_c_, _fh_Pcr_, _fh_Pcc_, _fh_Pr_, _fh_Pc_)
#define _fh_Pa_(f) f(_fh_c_, _fh_Pr_, _fh_Pc_, _fh_Pr_, _fh_Pc_)
#define _fh_Da_(f) f(_fh_s_, _fh_Dr_, _fh_Dc_, _fh_Dr_, _fh_Dc_)
#define _fh_Ia_(f) f(_fh_s_, _fh_Ir_, _fh_Ic_, _fh_Ir_, _fh_Ic_)
#define _fh_Oa_(f) f(_fh_s_, _fh_Or_, _fh_Oc_, _fh_Or_, _fh_Oc_)

/****************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

extern void FORTRAN(fhsetflags)(INTEGER *error,
  CINTEGER *mssmpart, CINTEGER *fieldren, CINTEGER *tanbren,
  CINTEGER *higgsmix, CINTEGER *p2approx, CINTEGER *looplevel,
  CINTEGER *tl_running_mt, CINTEGER *tl_bot_resum,
  CINTEGER *tl_cplx_approx);

extern void FORTRAN(fhsetflagsstring)(INTEGER *error,
  const char *flags, const int flags_len);

extern void FORTRAN(fhretrieveflags)(INTEGER *error,
  INTEGER *mssmpart, INTEGER *fieldren, INTEGER *tanbren,
  INTEGER *higgsmix, INTEGER *p2approx, INTEGER *looplevel,
  INTEGER *tl_running_mt, INTEGER *tl_bot_resum,
  INTEGER *tl_cplx_approx);

extern void FORTRAN(fhretrieveflagsstring)(INTEGER *error,
  char *flags, const int flags_len);


#define argsSetSMPara(x,r,c,ra,ca) \
  r(invAlfa,[1]) x r(AlfasMZ,[1]) x r(GF,[1]) x \
  r(ME,[1]) x r(MU,[1]) x r(MD,[1]) x \
  r(MM,[1]) x r(MC,[1]) x r(MS,[1]) x \
  r(ML,[1]) x             r(MB,[1]) x \
  r(MW,[1]) x r(MZ,[1]) x \
  r(CKMlambda,[1]) x r(CKMA,[1]) x r(CKMrhobar,[1]) x r(CKMetabar,[1])

extern void FORTRAN(fhsetsmpara)(INTEGER *error,
  _fh_Fca_(argsSetSMPara));

extern void FORTRAN(fhretrievesmpara)(INTEGER *error,
  _fh_Fa_(argsSetSMPara));

#define argsGetSMPara(x,r,c,ra,ca) \
  ca(CKM,[3][3])

extern void FORTRAN(fhgetsmpara)(INTEGER *error,
  _fh_Fa_(argsGetSMPara));


#define argsOSPara(x,r,c,ra,ca) \
  r(scalefactor,[1]) x \
  r(MT,[1]) x r(TB,[1]) x \
  r(MA0,[1]) x r(MHp,[1]) x \
  r(M3SL,[1]) x r(M3SE,[1]) x r(M3SQ,[1]) x r(M3SU,[1]) x r(M3SD,[1]) x \
  r(M2SL,[1]) x r(M2SE,[1]) x r(M2SQ,[1]) x r(M2SU,[1]) x r(M2SD,[1]) x \
  r(M1SL,[1]) x r(M1SE,[1]) x r(M1SQ,[1]) x r(M1SU,[1]) x r(M1SD,[1]) x \
  c(MUE,[1]) x \
  c(Atau,[1]) x c(At,[1]) x c(Ab,[1]) x \
  c(Amu,[1]) x c(Ac,[1]) x c(As,[1]) x \
  c(Ae,[1]) x c(Au,[1]) x c(Ad,[1]) x \
  c(M1,[1]) x c(M2,[1]) x c(M3,[1])

#define argsQPara(x,r,c,ra,ca) \
  r(Qtau,[1]) x r(Qt,[1]) x r(Qb,[1])

extern void FORTRAN(fhsetpara)(INTEGER *error,
  _fh_Fca_(argsOSPara), _fh_Fca_(argsQPara));

extern void FORTRAN(fhretrievepara)(INTEGER *error,
  _fh_Fa_(argsOSPara), _fh_Fa_(argsQPara));

extern void FORTRAN(fhretrieveospara)(INTEGER *error,
  _fh_Fa_(argsOSPara));

extern void FORTRAN(fhsetslha)(INTEGER *error,
  CCOMPLEX *slhadata);

#define argsGetPara(x,r,c,ra,ca) \
  ra(MSf,[3][5][2]) x ca(USf,[3][5][2][2]) x \
  ra(MASf,[5][6]) x ca(UASf,[5][6][6]) x \
  ra(MCha,[2]) x ca(UCha,[2][2]) x ca(VCha,[2][2]) x \
  ra(MNeu,[4]) x ca(ZNeu,[4][4]) x \
  c(Deltab,[1]) x r(MGl,[1]) x ra(MHtree,[4]) x r(SAtree,[1])

extern void FORTRAN(fhgetpara)(INTEGER *error,
  INTEGER *nmfv, _fh_Fa_(argsGetPara));

#define argsGetTLPara(x,r,c,ra,ca) \
  ra(MSb,[2]) x ca(USb,[2][2]) x \
  r(MbSL2,[1]) x c(Deltab,[1])

extern void FORTRAN(fhgettlpara)(INTEGER *error,
  _fh_Fa_(argsGetTLPara));


#define argsSetLFV(x,r,c,ra,ca) \
  c(deltaLLL12,[1]) x c(deltaLLL23,[1]) x c(deltaLLL13,[1]) x \
  c(deltaELR12,[1]) x c(deltaELR23,[1]) x c(deltaELR13,[1]) x \
  c(deltaERL12,[1]) x c(deltaERL23,[1]) x c(deltaERL13,[1]) x \
  c(deltaERR12,[1]) x c(deltaERR23,[1]) x c(deltaERR13,[1])

extern void FORTRAN(fhsetlfv)(INTEGER *error,
  _fh_Fca_(argsSetLFV));

extern void FORTRAN(fhretrievelfv)(INTEGER *error,
  _fh_Fa_(argsSetLFV));

#define argsSetNMFV(x,r,c,ra,ca) \
  c(deltaQLL12,[1]) x c(deltaQLL23,[1]) x c(deltaQLL13,[1]) x \
  c(deltaULR12,[1]) x c(deltaULR23,[1]) x c(deltaULR13,[1]) x \
  c(deltaURL12,[1]) x c(deltaURL23,[1]) x c(deltaURL13,[1]) x \
  c(deltaURR12,[1]) x c(deltaURR23,[1]) x c(deltaURR13,[1]) x \
  c(deltaDLR12,[1]) x c(deltaDLR23,[1]) x c(deltaDLR13,[1]) x \
  c(deltaDRL12,[1]) x c(deltaDRL23,[1]) x c(deltaDRL13,[1]) x \
  c(deltaDRR12,[1]) x c(deltaDRR23,[1]) x c(deltaDRR13,[1])

extern void FORTRAN(fhsetnmfv)(INTEGER *error,
  _fh_Fca_(argsSetNMFV));

extern void FORTRAN(fhretrievenmfv)(INTEGER *error,
  _fh_Fa_(argsSetNMFV));

#define argsGetNMFV(x,r,c,ra,ca) \
  ca(MSS2,[5][3][3]) x ca(Kf,[3][3][3])

extern void FORTRAN(fhgetnmfv)(INTEGER *error,
  _fh_Fa_(argsGetNMFV));


#define argsHiggsCorr(x,r,c,ra,ca) \
  ra(MHiggs,[4]) x c(SAeff,[1]) x \
  ca(UHiggs,[3][3]) x ca(ZHiggs,[3][3])

extern void FORTRAN(fhhiggscorr)(INTEGER *error,
  _fh_Fa_(argsHiggsCorr));

#define argsUncertainties(x,r,c,ra,ca) \
  ra(DeltaMHiggs,[4]) x c(DeltaSAeff,[1]) x \
  ca(DeltaUHiggs,[3][3]) x ca(DeltaZHiggs,[3][3])

extern void FORTRAN(fhuncertainties)(INTEGER *error,
  _fh_Fa_(argsUncertainties));


#define argsCouplings(x,r,c,ra,ca) \
  ca(couplings,[ncouplings]) x ca(couplingsms,[ncouplingsms]) x \
  ra(gammas,[ngammas]) x ra(gammasms,[ngammasms])

extern void FORTRAN(fhcouplings)(INTEGER *error,
  _fh_Fa_(argsCouplings), CINTEGER *fast);


#define argsConstraints(x,r,c,ra,ca) \
  r(gm2,[1]) x r(DeltaRho,[1]) x \
  r(MWMSSM,[1]) x r(MWSM,[1]) x \
  r(SW2MSSM,[1]) x r(SW2SM,[1]) x \
  r(EDMeTh,[1]) x r(EDMn,[1]) x r(EDMHg,[1])

extern void FORTRAN(fhconstraints)(INTEGER *error,
  _fh_Fa_(argsConstraints), INTEGER *ccb);


#define argsFlavour(x,r,c,ra,ca) \
  r(BsgammaMSSM,[1]) x r(BsgammaSM,[1]) x \
  r(DeltaMsMSSM,[1]) x r(DeltaMsSM,[1]) x \
  r(BsmumuMSSM,[1]) x r(BsmumuSM,[1])

extern void FORTRAN(fhflavour)(INTEGER *error,
  _fh_Fa_(argsFlavour));


extern void FORTRAN(fhhiggsprod)(INTEGER *error,
  _fh_Fcr_(sqrts,[1]), _fh_Fr_(prodxs,[nprodxs]));

extern void FORTRAN(fhgetself)(INTEGER *error, _fh_Fcr_(k2,[1]),
  CINTEGER *key, _fh_Fc_(sig,[nsig]),
  CINTEGER *dkey, _fh_Fc_(dsig,[nsig]));

extern void FORTRAN(fhaddself)(INTEGER *error,
  _fh_Fcc_(sig,[nsig]), CINTEGER *rotate);

extern void FORTRAN(fhoutput)(INTEGER *error,
  const char *filename, CINTEGER *key, _fh_Fcr_(sqrts,[1]),
  const int filename_len);

extern void FORTRAN(fhoutputslha)(INTEGER *error,
  COMPLEX *slhadata, CINTEGER *key);

extern void FORTRAN(fhrecordindex)(INTEGER *index,
  CCHARACTER *para, INTEGER para_len);

#define argsRecord(x,r,c,ra,ca) \
  ra(record,[nrecord])

extern void FORTRAN(fhclearrecord)(_fh_Fa_(argsRecord));

extern void FORTRAN(fhreadrecord)(INTEGER *error,
  _fh_Fa_(argsRecord), COMPLEX *slhadata,
  CCHARACTER *inputfile, INTEGER inputfile_len);

extern void FORTRAN(fhslharecord)(INTEGER *error,
  _fh_Fa_(argsRecord), CCOMPLEX *slhadata);

extern void FORTRAN(fhlooprecord)(INTEGER *error,
  _fh_Fa_(argsRecord));

extern void FORTRAN(fhsetrecord)(INTEGER *error,
  _fh_Fca_(argsRecord));

extern void FORTRAN(fhretrieverecord)(INTEGER *error,
  _fh_Fa_(argsRecord), CINTEGER *iX);

extern void FORTRAN(fhloadtable)(INTEGER *error,
  CCHARACTER *inputfile, CINTEGER *inputunit,
  INTEGER inputfile_len);

extern void FORTRAN(fhtablerecord)(INTEGER *error,
  _fh_Fa_(argsRecord), CINTEGER *ind1, CINTEGER *ind2);


extern void FORTRAN(fhsetdebug)(CINTEGER *debuglevel);

extern void FORTRAN(fhselectuz)(INTEGER *error,
  CINTEGER *uzint, CINTEGER *uzext, CINTEGER *mfeff);

#ifdef __cplusplus
}
#endif

/****************************************************************/

static inline void FHSetFlags(int *error,
  const int mssmpart, const int fieldren, const int tanbren,
  const int higgsmix, const int p2approx, const int looplevel,
  const int tl_running_mt, const int tl_bot_resum,
  const int tl_cplx_approx)
{
  FORTRAN(fhsetflags)(error, &mssmpart, &fieldren, &tanbren,
    &higgsmix, &p2approx, &looplevel,
    &tl_running_mt, &tl_bot_resum, &tl_cplx_approx);
}


static inline void FHSetFlagsString(int *error, const char *flags)
{
  FORTRAN(fhsetflagsstring)(error, flags, strlen(flags));
}


static inline void FHRetrieveFlags(int *error,
  int *mssmpart, int *fieldren, int *tanbren,
  int *higgsmix, int *p2approx, int *looplevel,
  int *runningMT, int *botResum, int *tlCplxApprox)
{
  FORTRAN(fhretrieveflags)(error, mssmpart, fieldren, tanbren,
    higgsmix, p2approx, looplevel,
    runningMT, botResum, tlCplxApprox);
}


#define FHRetrieveFlagsString(error, flags) { \
  FORTRAN(fhretrieveflagsstring)(error, flags, sizeof(flags) - 1); \
  if( *(error) == 0 ) (flags)[9] = 0; \
}


static inline void FHSetSMPara(int *error,
  _fh_Cca_(argsSetSMPara))
{
  _fh_Da_(argsSetSMPara);
  _fh_Ia_(argsSetSMPara);
  FORTRAN(fhsetsmpara)(error, _fh_Pca_(argsSetSMPara));
}


static inline void FHRetrieveSMPara(int *error,
  _fh_Ca_(argsSetSMPara))
{
  _fh_Da_(argsSetSMPara);
  FORTRAN(fhretrievesmpara)(error, _fh_Pa_(argsSetSMPara));
  _fh_Oa_(argsSetSMPara);
}


static inline void FHGetSMPara(int *error,
  _fh_Ca_(argsGetSMPara))
{
  _fh_Da_(argsGetSMPara);
  FORTRAN(fhgetsmpara)(error, _fh_Pa_(argsGetSMPara));
  _fh_Oa_(argsGetSMPara);
}


static inline void FHSetPara(int *error,
  _fh_Cca_(argsOSPara), _fh_Cca_(argsQPara))
{
  _fh_Da_(argsOSPara);
  _fh_Da_(argsQPara);
  _fh_Ia_(argsOSPara);
  _fh_Ia_(argsQPara);
  FORTRAN(fhsetpara)(error, _fh_Pca_(argsOSPara), _fh_Pca_(argsQPara));
}


static inline void FHRetrievePara(int *error,
  _fh_Ca_(argsOSPara), _fh_Ca_(argsQPara))
{
  _fh_Da_(argsOSPara);
  _fh_Da_(argsQPara);
  FORTRAN(fhretrievepara)(error, _fh_Pa_(argsOSPara), _fh_Pa_(argsQPara));
  _fh_Oa_(argsOSPara);
  _fh_Oa_(argsQPara);
}


static inline void FHRetrieveOSPara(int *error,
  _fh_Ca_(argsOSPara))
{
  _fh_Da_(argsOSPara);
  FORTRAN(fhretrieveospara)(error, _fh_Pa_(argsOSPara));
  _fh_Oa_(argsOSPara);
}


static inline void FHSetSLHA(int *error,
  CCOMPLEX *slhadata)
{
  FORTRAN(fhsetslha)(error, slhadata);
}


static inline void FHSetLFV(int *error,
  _fh_Cca_(argsSetLFV))
{
  _fh_Da_(argsSetLFV);
  _fh_Ia_(argsSetLFV);
  FORTRAN(fhsetlfv)(error, _fh_Pca_(argsSetLFV));
}

static inline void FHRetrieveLFV(int *error,
  _fh_Ca_(argsSetLFV))
{
  _fh_Da_(argsSetLFV);
  FORTRAN(fhretrievelfv)(error, _fh_Pa_(argsSetLFV));
  _fh_Oa_(argsSetLFV);
}

static inline void FHSetNMFV(int *error,
  _fh_Cca_(argsSetNMFV))
{
  _fh_Da_(argsSetNMFV);
  _fh_Ia_(argsSetNMFV);
  FORTRAN(fhsetnmfv)(error, _fh_Pca_(argsSetNMFV));
}

static inline void FHRetrieveNMFV(int *error,
  _fh_Ca_(argsSetNMFV))
{
  _fh_Da_(argsSetNMFV);
  FORTRAN(fhretrievenmfv)(error, _fh_Pa_(argsSetNMFV));
  _fh_Oa_(argsSetNMFV);
}


static inline void FHGetPara(int *error,
  int *nmfv, _fh_Ca_(argsGetPara))
{
  _fh_Da_(argsGetPara);
  FORTRAN(fhgetpara)(error, nmfv, _fh_Pa_(argsGetPara));
  _fh_Oa_(argsGetPara);
}


static inline void FHGetTLPara(int *error,
  _fh_Ca_(argsGetTLPara))
{
  _fh_Da_(argsGetTLPara);
  FORTRAN(fhgettlpara)(error, _fh_Pa_(argsGetTLPara));
  _fh_Oa_(argsGetTLPara);
}


static inline void FHGetNMFV(int *error,
  _fh_Ca_(argsGetNMFV))
{
  _fh_Da_(argsGetNMFV);
  FORTRAN(fhgetnmfv)(error, _fh_Pa_(argsGetNMFV));
  _fh_Oa_(argsGetNMFV);
}


static inline void FHHiggsCorr(int *error,
  _fh_Ca_(argsHiggsCorr))
{
  _fh_Da_(argsHiggsCorr);
  FORTRAN(fhhiggscorr)(error, _fh_Pa_(argsHiggsCorr));
  _fh_Oa_(argsHiggsCorr);
}


static inline void FHUncertainties(int *error,
  _fh_Ca_(argsUncertainties))
{
  _fh_Da_(argsUncertainties);
  FORTRAN(fhuncertainties)(error, _fh_Pa_(argsUncertainties));
  _fh_Oa_(argsUncertainties);
}


static inline void FHCouplings(int *error,
  _fh_Ca_(argsCouplings), const int fast)
{
  _fh_Da_(argsCouplings);
  FORTRAN(fhcouplings)(error, _fh_Pa_(argsCouplings), &fast);
  _fh_Oa_(argsCouplings);
}


static inline void FHConstraints(int *error,
  _fh_Ca_(argsConstraints), int *ccb)
{
  _fh_Da_(argsConstraints);
  FORTRAN(fhconstraints)(error, _fh_Pa_(argsConstraints), ccb);
  _fh_Oa_(argsConstraints);
}


static inline void FHFlavour(int *error,
  _fh_Ca_(argsFlavour))
{
  _fh_Da_(argsFlavour);
  FORTRAN(fhflavour)(error, _fh_Pa_(argsFlavour));
  _fh_Oa_(argsFlavour);
}


static inline void FHHiggsProd(int *error,
  _fh_Ccr_(sqrts,[1]), _fh_Car_(prodxs,[nprodxs]))
{
  _fh_Dr_(sqrts,[1]);
  _fh_Dr_(prodxs,[nprodxs]);
  _fh_Ir_(sqrts,[1]);
  FORTRAN(fhhiggsprod)(error, _fh_Pcr_(sqrts,[1]), _fh_Pr_(prodxs,[nprodxs]));
  _fh_Or_(prodxs,[nprodxs]);
}


static inline void FHGetSelf(int *error, _fh_Ccr_(k2,[1]),
  const int key, _fh_Cac_(sig,[nsig]),
  const int dkey, _fh_Cac_(dsig,[nsig]))
{
  _fh_Dr_(k2,[1]);
  _fh_Dc_(sig,[nsig]);
  _fh_Dc_(dsig,[nsig]);
  _fh_Ir_(k2,[1]);
  FORTRAN(fhgetself)(error, _fh_Pcr_(k2,[1]),
    &key, _fh_Pc_(sig,[nsig]), &dkey, _fh_Pc_(dsig,[nsig]));
  _fh_Oc_(sig,[nsig]);
  _fh_Oc_(dsig,[nsig]);
}

#define Key(se) (1 << (se-1))


static inline void FHAddSelf(int *error,
  _fh_Ccac_(sig,[nsig]), const int rotate)
{
  _fh_Dc_(sig,[nsig]);
  _fh_Ic_(sig,[nsig]);
  FORTRAN(fhaddself)(error, _fh_Pcac_(sig,[nsig]), &rotate);
}


static inline void FHOutput(int *error, const char *filename,
  const int key, _fh_Ccr_(sqrts,[1]))
{
  _fh_Dr_(sqrts,[1]);
  _fh_Ir_(sqrts,[1]);
  FORTRAN(fhoutput)(error, filename, &key, _fh_Pcr_(sqrts,[1]),
    strlen(filename));
}


static inline void FHOutputSLHA(int *error,
  COMPLEX *slhadata, const int key)
{
  FORTRAN(fhoutputslha)(error, slhadata, &key);
}


static inline void FHRecordIndex(int *index, const char *para)
{
  FORTRAN(fhrecordindex)(index, para, strlen(para));
}


static inline void FHClearRecord(_fh_Ca_(argsRecord))
{
  _fh_Da_(argsRecord);
  FORTRAN(fhclearrecord)(_fh_Pa_(argsRecord));
  _fh_Oa_(argsRecord);
}


static inline void FHReadRecord(int *error,
  _fh_Ca_(argsRecord), COMPLEX *slhadata, const char *inputfile)
{
  _fh_Da_(argsRecord);
  FORTRAN(fhreadrecord)(error, _fh_Pa_(argsRecord), slhadata,
    inputfile, strlen(inputfile));
  _fh_Oa_(argsRecord);
}


static inline void FHSLHARecord(int *error,
  _fh_Ca_(argsRecord), CCOMPLEX *slhadata)
{
  _fh_Da_(argsRecord);
  FORTRAN(fhslharecord)(error, _fh_Pa_(argsRecord), slhadata);
  _fh_Oa_(argsRecord);
}


static inline void FHLoopRecord(int *error, _fh_Ca_(argsRecord))
{
  _fh_Da_(argsRecord);
  _fh_Ia_(argsRecord);
  FORTRAN(fhlooprecord)(error, _fh_Pa_(argsRecord));
  _fh_Oa_(argsRecord);
}


static inline void FHSetRecord(int *error, _fh_Cca_(argsRecord))
{
  _fh_Da_(argsRecord);
  _fh_Ia_(argsRecord);
  FORTRAN(fhsetrecord)(error, _fh_Pca_(argsRecord));
}


static inline void FHRetrieveRecord(int *error,
  _fh_Ca_(argsRecord), const int iX)
{
  _fh_Da_(argsRecord);
  FORTRAN(fhretrieverecord)(error, _fh_Pa_(argsRecord), &iX);
  _fh_Oa_(argsRecord);
}


static inline void FHLoadTable(int *error, const char *inputfile)
{
  int inputunit = 5;
  FORTRAN(fhloadtable)(error, inputfile, &inputunit, strlen(inputfile));
}


static inline void FHTableRecord(int *error,
  _fh_Ca_(argsRecord), const int ind1, const int ind2)
{
  _fh_Da_(argsRecord);
  _fh_Ia_(argsRecord);
  FORTRAN(fhtablerecord)(error, _fh_Pa_(argsRecord), &ind1, &ind2);
  _fh_Oa_(argsRecord);
}


static inline void FHSetDebug(const int debuglevel)
{
  FORTRAN(fhsetdebug)(&debuglevel);
}


static inline void FHSelectUZ(int *error,
  const int uzint, const int uzext, const int mfeff)
{
  FORTRAN(fhselectuz)(error, &uzint, &uzext, &mfeff);
}

#endif

