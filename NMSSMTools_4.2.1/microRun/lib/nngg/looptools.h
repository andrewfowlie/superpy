* looptools.h
* the header file for Fortran with all definitions for LoopTools
* this file is part of LoopTools
* last modified 16 Dec 99 th


#ifndef _LT_GLOBAL_DEFS_
#define _LT_GLOBAL_DEFS_

#define cc0 1
#define cc1 2
#define cc2 3
#define cc00 4
#define cc11 5
#define cc12 6
#define cc22 7
#define cc001 8
#define cc002 9
#define cc111 10
#define cc112 11
#define cc122 12
#define cc222 13

#define dd0 1
#define dd1 2
#define dd2 3
#define dd3 4
#define dd00 5
#define dd11 6
#define dd12 7
#define dd13 8
#define dd22 9
#define dd23 10
#define dd33 11
#define dd001 12
#define dd002 13
#define dd003 14
#define dd111 15
#define dd112 16
#define dd113 17
#define dd122 18
#define dd123 19
#define dd133 20
#define dd222 21
#define dd223 22
#define dd233 23
#define dd333 24
#define dd0000 25
#define dd0011 26
#define dd0012 27
#define dd0013 28
#define dd0022 29
#define dd0023 30
#define dd0033 31
#define dd1111 32
#define dd1112 33
#define dd1113 34
#define dd1122 35
#define dd1123 36
#define dd1133 37
#define dd1222 38
#define dd1223 39
#define dd1233 40
#define dd1333 41
#define dd2222 42
#define dd2223 43
#define dd2233 44
#define dd2333 45
#define dd3333 46

* for compatibility:

#define Cval(id, pos) Ccache(pos + id)
#define Dval(id, pos) Dcache(pos + id)

#define bcaini ffini
#define bcaexi ffexi

#endif


	double complex Ccache(1)
	common /Cbase/ Ccache

	double complex Dcache(1)
	common /Dbase/ Dcache

	double complex A0, B0, DB0, C0, D0
	double complex B1, DB1, B00, DB00, B11, DB11
	double complex C0i, D0i
	integer*8 Cget, Dget, getcachelast
	double precision getmudim, getdelta, getlambda

	external A0, B0, DB0, C0, D0
	external B1, DB1, B00, DB00, B11, DB11
	external C0i, D0i
	external Cget, Dget, getcachelast
	external getmudim, getdelta, getlambda

