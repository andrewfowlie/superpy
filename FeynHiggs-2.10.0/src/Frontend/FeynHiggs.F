* FeynHiggs.F
* the command-line front end for FeynHiggs
* this file is part of FeynHiggs
* last modified 18 Dec 13 th

#include "FHCouplings.h"
#include "types.h"
#include "debug.h"
#include "version.h"


	program FeynHiggs
	implicit none

#include "FHRecord.h"
#include "SLHA.h"

	character*128 argv
	integer argc, type, error, loop, keyslha, i
	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, runningMT, botResum, tlCplxApprox
	RecordDecl(record)
	ComplexType slhadata(nslhadata)

#if U77EXT
#define len_trim lnblnk
	integer iargc, lnblnk
	external iargc, lnblnk
#endif

	character*16 flags

	character*55 dmssmpart(-1:5), dfieldren(-1:5),
     &    dtanbren(-1:3), dhiggsmix(0:5), dp2approx(-1:5),
     &    dlooplevel(-1:4), drunningMT(-1:2),
     &    dbotResum(-1:3), dtlCplxApprox(-1:7)

	data flags /"400242110"/

	data dmssmpart /
     &    "???",
     &    "MT^4 approximation (for real parameters only)",
     &    "top/stop sector",
     &    "top/stop + bottom/sbottom sector",
     &    "full (s)quark/(s)lepton sector",
     &    "full MSSM [recommended]",
     &    "???" /

	data dfieldren /
     &    "???",
     &    "DRbar field ren. [strongly recommended]",
     &    "on-shell field ren., Dabelstein convention",
     &    "on-shell field ren., `Goldstone-pole' version",
     &    "on-shell field ren., MF I",
     &    "on-shell field ren., MF II",
     &    "???" /

	data dtanbren /
     &    "???",
     &    "DRbar TB-ren. [strongly recommended]",
     &    "field renormalization part only for TB (Dabelstein)",
     &    "mixed field/on-shell renormalization for TB",
     &    "???" /

	data dhiggsmix /
     &    "???",
     &    "no mixing in the neutral Higgs sector",
     &    "2x2 (h0-HH) mixing in the neutral Higgs sector",
     &    "3x3 (h0-HH-A0) mixing in the neutral Higgs sector",
     &    "4x4 (G0-h0-HH-A0) mixing [experimental]",
     &    "???" /

	data dp2approx /
     &    "???",
     &    "no approximation [recommended]",
     &    "p^2 = 0 approximation",
     &    "self-energies are calculated at p^2 = M_tree^2",
     &    "imaginary parts of self-energies are discarded",
     &    "UHiggs is evaluated at p^2 = 0",
     &    "???" /

	data dlooplevel /
     &    "???",
     &    "tree level only",
     &    "include 1-loop corrections",
     &    "include 2-loop corrections",
     &    "include resummation of large logs",
     &    "???" /

	data drunningMT /
     &    "???",
     &    "use pole top mass",
     &    "use running top mass",
     &    "???" /

	data dbotResum /
     &    "???",
     &    "non-resummed MB in 1L/2L corrections",
     &    "resummed 1L-MB in 1L/2L corrections [recommended]",
     &    "resummed 2L-MB in 1L/2L corrections",
     &    "???" /

	data dtlCplxApprox /
     &    "???",
     &    "all corrections evaluated in the rMSSM",
     &    "a_s a_t in cMSSM only",
     &    "a_s a_t in cMSSM, rest in rMSSM",
     &    "a_s a_t in cMSSM, rest interpolated in At,Ab,M_3,MUE",
     &    "a_s a_t in cMSSM, rest interpolated in At,Xb,M_3,MUE",
     &    "a_s a_t in cMSSM, rest interpolated in Xt,Ab,M_3,MUE",
     &    "a_s a_t in cMSSM, rest interpolated in Xt,Xb,M_3,MUE",
     &    "???" /

	call pageron

	argc = iargc()
	if( argc .lt. 1 .or. argc .gt. 3 ) goto 999

	if( argc .ge. 2 ) call getarg(2, flags)
	call FHSetFlagsString(error, flags)

	call FHRetrieveFlags(error,
     &    mssmpart, fieldren, tanbren, higgsmix, p2approx,
     &    looplevel, runningMT, botResum, tlCplxApprox)

#define ShowFlag(t,f,lo,hi) t, f, d f(min(max(f,lo-1),hi+1))

101	format(A, " = ", I1, T20, "# ", A)
	print 101,
     &    ShowFlag("mssmpart",     mssmpart,     0, 4),
     &    ShowFlag("fieldren",     fieldren,     0, 4),
     &    ShowFlag("tanbren",      tanbren,      0, 2),
     &    ShowFlag("higgsmix",     higgsmix,     1, 4),
     &    ShowFlag("p2approx",     p2approx,     0, 4),
     &    ShowFlag("looplevel",    looplevel,    0, 3),
     &    ShowFlag("runningMT",    runningMT,    0, 1),
     &    ShowFlag("botResum",     botResum,     0, 2),
     &    ShowFlag("tlCplxApprox", tlCplxApprox, 0, 6)

	if( error .ne. 0 ) stop

	call getarg(1, argv)
	i = index(argv, "#")
	type = -1
	if( i .gt. 0 ) then
	  type = 2
	  if( argv(i+1:) .eq. "SLHA" ) type = 0
	  argv(i:) = " "
	endif

	call FHReadRecord(error, record, slhadata, argv)
* 0 = SLHA, 2 = FH file, all others are true errors
	if( ibclr(error, 1) .ne. 0 ) goto 1

#if SLHAPARA
	keyslha = -2
#else
	keyslha = -4
#endif
	if( type .eq. 0 .and. error .ne. 0 ) then
	  call SLHAClear(slhadata)
	  keyslha = -1
	endif

	if( type .eq. -1 ) type = error

	if( type .eq. 0 ) then
	  i = 0
	  if( argv .ne. "-" ) i = len_trim(argv) + 1
	endif

* turn on verbose mode:
	record(iAdmin,iStep) = 1

	loop = 0

2	call FHLoopRecord(error, record)
	if( error .eq. -2 ) goto 2
	if( error .ne. 0 ) goto 1

	loop = loop + 1

	call FHSetRecord(error, record)
	if( error .ne. 0 ) goto 2

	if( type .eq. 0 ) then
	  call FHOutputSLHA(error, slhadata, keyslha)
	  if( error .ne. 0 ) print *, "Error in FHOutputSLHA ", error
	  if( i .gt. 0 ) write(argv(i:),'(".fh-",I3.3)') loop
	  call SLHAWrite(error, slhadata, argv)
	else
	  call FHOutput(error, "stdout", -1, record(iprodSqrts,iVar))
	  if( error .ne. 0 ) print *, "Error in FHOutput ", error
	  print 200
	endif

	goto 2


#define DESC_FMT (/I5,": ",A)

200	format(/"-------------------- END OF OUTPUT ---------------------")

100	format ("FeynHiggs ", FHVERSION
     &    /"built on ", FHBUILDDATE
     &    //"USAGE:  ", A, " <inputfile> [<flags>]"
     &    //"<inputfile> is a file containing the parameters"
     &    //"<flags> is a string of 9 digits with"
     &    /"  #1 = mssmpart = scope of the calculation",
     &    5 DESC_FMT
     &    /"  #2 = fieldren = 1-loop field renormalization",
     &    5 DESC_FMT
     &    /"  #3 = tanbren = 1-loop tan(beta) renormalization",
     &    3 DESC_FMT
     &    /"  #4 = higgsmix = mixing in the neutral Higgs sector",
     &    4 DESC_FMT
     &    /"  #5 = p2approx = approximation for the 1-loop result",
     &    5 DESC_FMT
     &    /"  #6 = looplevel, inclusion of higher-order corrections",
     &    4 DESC_FMT
     &    /"  #7 = runningMT, pole vs. running MT",
     &    2 DESC_FMT
     &    /"  #8 = botResum, MB resummation",
     &    3 DESC_FMT
     &    /"  #9 = tlCplxApprox, complex parameters in the 2-loop corrections",
     &    7 DESC_FMT
     &    /)

999	call getarg(0, argv)
	print 100, Strip(argv),
     &    (i, dmssmpart(i),     i = 0, 4),
     &    (i, dfieldren(i),     i = 0, 4),
     &    (i, dtanbren(i),      i = 0, 2),
     &    (i, dhiggsmix(i),     i = 1, 4),
     &    (i, dp2approx(i),     i = 0, 4),
     &    (i, dlooplevel(i),    i = 0, 3),
     &    (i, drunningMT(i),    i = 0, 1),
     &    (i, dbotResum(i),     i = 0, 2),
     &    (i, dtlCplxApprox(i), i = 0, 6)

1	call flush(6)
	call pageroff
	end

