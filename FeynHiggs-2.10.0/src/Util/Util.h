* Util.h
* limits and user preferences for the utility routines
* this file is part of FeynHiggs
* last modified 23 Apr 10 th


* for Eigensystem, SingularValues, PseudoEigensystem:
* the maximum size of the matrix

#define MAXDIM 8
* i.e. 8x8


* A matrix is considered diagonal if the sum of the squares
* of the off-diagonal elements is less than EPS.  SYM_EPS is
* half of EPS since only the upper triangle is counted for
* symmetric matrices.
* 52 bits is the mantissa length for IEEE double precision.

#define EPS 2D0**(-102)

#define SYM_EPS 2D0**(-103)

#define DBL_EPS 2D0**(-52)


* for Gauss:
* the maximum number of components of the integrand vector

#define MAXCOMP 8

