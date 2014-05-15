#ifndef __COMPLEX_H__
#define __COMPLEX_H__

extern fcomplex Cadd(fcomplex a, fcomplex b);
extern fcomplex Csub(fcomplex a, fcomplex b);
extern fcomplex Cmul(fcomplex a, fcomplex b);
extern fcomplex Conjg(fcomplex z);
extern fcomplex COmplex(float re, float im);
extern fcomplex Cdiv(fcomplex a, fcomplex b);
extern float	Cabs(fcomplex z);
extern fcomplex Csqrt(fcomplex z);
extern fcomplex RCmul(float x, fcomplex a);
#endif
