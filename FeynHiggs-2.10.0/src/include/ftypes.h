#ifndef FTYPES_H
#define FTYPES_H

#if NOUNDERSCORE
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

