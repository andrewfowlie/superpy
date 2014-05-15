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

/*
	CSLHA.h
		C/C++ wrapper functions for the SLHALib subroutines
		last modified 20 Sep 12 th
*/


#ifndef CSLHA_H
#define CSLHA_H

#include <string.h>
#include <stdarg.h>
#include "SLHADefs.h"
#include "PDG.h"

#define Slhadata(i) ToREAL(slhadata[i-1])
#define SlhaData(i) ToREAL(slhadata[i-1])


#ifdef __cplusplus
extern "C" {
#endif

extern void FORTRAN(slhaclear)(COMPLEX *slhadata);

extern void FORTRAN(slharead)(INTEGER *error,
  COMPLEX *slhadata, const char *filename,
  CINTEGER *abort, const int filename_len);

extern void FORTRAN(slhawrite)(INTEGER *error,
  CCOMPLEX *slhadata, const char *filename,
  const int filename_len);

extern void FORTRAN(slhaputinfo)(COMPLEX *slhablock,
  CINTEGER *code, const char *text, const int text_len);

extern void FORTRAN(slhagetinfo)(COMPLEX *slhablock,
  CINTEGER *i, char *text, const int text_len);

extern INTEGER FORTRAN(slhanewdecay)(COMPLEX *slhadata,
  CREAL *width, CINTEGER *parent_id);

extern INTEGER FORTRAN(slhafinddecay)(COMPLEX *slhadata,
  CINTEGER *parent_id);

extern void FORTRAN(slhaadddecay)(COMPLEX *slhadata,
  CREAL *br, CINTEGER *decay, CINTEGER *nchildren,
  CINTEGER *child1_id, CINTEGER *child2_id,
  CINTEGER *child3_id, CINTEGER *child4_id);

extern REAL FORTRAN(slhagetdecay)(
  CCOMPLEX *slhadata,
  CINTEGER *parent_id, CINTEGER *nchildren,
  CINTEGER *child1_id, CINTEGER *child2_id,
  CINTEGER *child3_id, CINTEGER *child4_id);

extern INTEGER FORTRAN(slhadecaytable)(CCOMPLEX *slhadata,
  CINTEGER *parent_id, REAL *width, INTEGER *id,
  CINTEGER *maxparticles, CINTEGER *maxchannels);

extern INTEGER FORTRAN(slhaexist)(CCOMPLEX *slhablock,
  CINTEGER *length);

extern INTEGER FORTRAN(slhavalid)(CCOMPLEX *slhablock,
  CINTEGER *length);

extern void FORTRAN(slhapdgname)(CINTEGER *code,
  const char *name, const int name_len);

#ifdef __cplusplus
}
#endif


static inline void SLHAClear(COMPLEX *slhadata)
{
  FORTRAN(slhaclear)(slhadata);
}


static inline void SLHARead(int *error,
  COMPLEX *slhadata, const char *filename, const int abort)
{
  FORTRAN(slharead)(error, slhadata,
    filename, &abort, strlen(filename));
}


static inline void SLHAWrite(int *error,
  CCOMPLEX *slhadata, const char *filename)
{
  FORTRAN(slhawrite)(error, slhadata,
    filename, strlen(filename));
}


static inline void SLHAPutInfo(COMPLEX *slhablock,
  const int code, const char *text)
{
  FORTRAN(slhaputinfo)(slhablock, &code,
    text, strlen(text));
}


static inline void SLHAGetInfo(ComplexType *slhablock,
  const int i, char *text, const int textsize)
{
  int off;

  FORTRAN(slhagetinfo)((COMPLEX *)slhablock, &i,
    text, textsize);

  off = textsize - 1;
  while( off && text[off - 1] == ' ' ) --off;
  text[off] = 0;
}


static inline int SLHANewDecay(COMPLEX *slhadata,
  cRealType width, const int parent_id)
{
  CREAL width_ = ToREAL(width);
  return FORTRAN(slhanewdecay)(slhadata, &width_, &parent_id);
}


static inline int SLHAFindDecay(COMPLEX *slhadata,
  const int parent_id)
{
  return FORTRAN(slhafinddecay)(slhadata, &parent_id);
}


static inline void SLHAAddDecay(COMPLEX *slhadata,
  cRealType br, const int decay,
  const int nchildren, const int child1_id, ...)
{
  va_list child_ids;
  int i, child_id[4];
  CREAL br_ = ToREAL(br);

  va_start(child_ids, child1_id);
  for( i = 0; i < nchildren; ++i )
    child_id[i] = va_arg(child_ids, int);
  va_end(child_ids);

  FORTRAN(slhaadddecay)(slhadata,
    &br_, &decay, &nchildren,
    &child_id[0], &child_id[1], &child_id[2], &child_id[3]);
}


static inline void SLHAGetDecay(COMPLEX *slhadata,
  const int parent_id,
  const int nchildren, const int child1_id, ...)
{
  va_list child_ids;
  int i, child_id[4];

  va_start(child_ids, child1_id);
  for( i = 0; i < nchildren; ++i )
    child_id[i] = va_arg(child_ids, int);
  va_end(child_ids);

  FORTRAN(slhagetdecay)(slhadata,
    &parent_id, &nchildren,
    &child_id[0], &child_id[1], &child_id[2], &child_id[3]);
}


static inline int SLHADecayTable(CCOMPLEX *slhadata,
  const int parent_id, double *width, int *id,
  const int maxparticles, const int maxchannels)
{
  REAL width_;
  const int ret = FORTRAN(slhadecaytable)(slhadata,
    &parent_id, &width_, id, &maxparticles, &maxchannels);
  *width = ToReal(width_);
  return ret;
}


static inline int SLHAExist(cComplexType *slhablock,
  const int length)
{
  return FORTRAN(slhaexist)((CCOMPLEX *)slhablock, &length);
}


static inline int SLHAValid(cComplexType *slhablock,
  const int length)
{
  return FORTRAN(slhavalid)((CCOMPLEX *)slhablock, &length);
}


static inline void SLHAPDGName(const int code,
  char name[PDGLen+1])
{
  char *n;

  FORTRAN(slhapdgname)(&code, name, PDGLen);
  n = (char *)memchr(name, ' ', PDGLen);
  *(n ? n : name + PDGLen) = 0;
}

#endif

