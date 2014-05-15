/*
	vegas-f.c
		Fortran interface for Vegas
		this file is part of Vegas
		last modified 16 Jul 04 th
*/


#include "decl.h"

#ifdef HAVE_UNDERSCORE
#define vegas vegas_
#endif

void Vegas(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount nstart, ccount nincrease,
  count *pneval, int *pfail,
  real *integral, real *error, real *prob);


void vegas(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, ccount *pmineval, ccount *pmaxeval,
  ccount *pnstart, ccount *pnincrease, 
  count *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  Vegas(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnstart, *pnincrease,
    pneval, pfail,
    integral, error, prob);
}

