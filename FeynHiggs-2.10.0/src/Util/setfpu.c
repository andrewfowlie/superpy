/*
	setfpu.c
		set FPU to throw an exception (terminate the
		program) rather than deliver NaNs
		necessary only for g77/gfortran
		this file is part of FeynHiggs
		last modified 20 Jul 11 th
*/


#include <stdlib.h>

#ifdef __GLIBC_PREREQ
#if __GLIBC_PREREQ(2,2)
#define _GNU_SOURCE
#define __USE_GNU
#include <fenv.h>
#define HAVE_FENV
#endif
#endif

#include "externals.h"

#if NOUNDERSCORE
#define setfpu_ setfpu
#endif

void setfpu_(void)
{
#ifdef HAVE_FENV
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
}

