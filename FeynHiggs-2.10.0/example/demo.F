* demo.F
* demonstration program for calling FeynHiggs from Fortran
* this file is part of FeynHiggs
* last modified 31 Jul 13 th

* compile this file with something like
* f77 -Ibuild example/demo.F -Lbuild -lFH


	program FHdemo
	implicit none

	call setFlags

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
	call setPara
c	call setSLHA("myfile.slha")

	call getPara

	call higgsCorr
	end


************************************************************************

	subroutine setFlags
	implicit none

	integer mssmpart, fieldren, tanbren, higgsmix, p2approx
	integer looplevel, runningMT, botResum, tlCplxApprox
	parameter (mssmpart = 4)
	parameter (fieldren = 0)
	parameter (tanbren = 0)
	parameter (higgsmix = 2)
	parameter (p2approx = 0)
	parameter (looplevel = 2)
	parameter (runningMT = 1)
	parameter (botResum = 1)
	parameter (tlCplxApprox = 3)

	integer error

	call FHSetFlags(error, mssmpart, fieldren, tanbren,
     &    higgsmix, p2approx, looplevel,
     &    runningMT, botResum, tlCplxApprox)
	if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setPara
	implicit none

	double precision invAlfa, AlfasMZ, GF
	double precision ME, MU, MD, MM, MC, MS, ML, MB, MZ, MW
	double precision CKMlambda, CKMA, CKMrhobar, CKMetabar
	parameter (invAlfa = -1)
	parameter (AlfasMZ = -1)
	parameter (GF = -1)
	parameter (ME = -1)
	parameter (MU = -1)
	parameter (MD = -1)
	parameter (MM = -1)
	parameter (MC = -1)
	parameter (MS = -1)
	parameter (ML = -1)
	parameter (MB = -1)
	parameter (MW = -1)
	parameter (MZ = -1)
	parameter (CKMlambda = -1)
	parameter (CKMA = -1)
	parameter (CKMrhobar = -1)
	parameter (CKMetabar = -1)

	double precision MT, TB, MA0, MHp
	parameter (MT = 172)
	parameter (TB = 5)
	parameter (MA0 = 250)
	parameter (MHp = -1)

	double precision MSusy
	double precision M3SL, M2SL, M1SL
	double precision M3SE, M2SE, M1SE
	double precision M3SQ, M2SQ, M1SQ
	double precision M3SU, M2SU, M1SU
	double precision M3SD, M2SD, M1SD
	parameter (MSusy = 1000)
	parameter (M3SL = MSusy)
	parameter (M2SL = M3SL)
	parameter (M1SL = M2SL)
	parameter (M3SE = MSusy)
	parameter (M2SE = M3SE)
	parameter (M1SE = M2SE)
	parameter (M3SQ = MSusy)
	parameter (M2SQ = M3SQ)
	parameter (M1SQ = M2SQ)
	parameter (M3SU = MSusy)
	parameter (M2SU = M3SU)
	parameter (M1SU = M2SU)
	parameter (M3SD = MSusy)
	parameter (M2SD = M3SD)
	parameter (M1SD = M2SD)

	double complex Atau, At, Ab
	double complex Amu, Ac, As
	double complex Ae, Au, Ad
	parameter (At = 2000)
	parameter (Ab = At)
	parameter (Atau = At)
	parameter (Ac = At)
	parameter (As = Ab)
	parameter (Amu = Atau)
	parameter (Au = Ac)
	parameter (Ad = As)
	parameter (Ae = Amu)

	double complex MUE, M_1, M_2, M_3
	parameter (MUE = 200)
	parameter (M_1 = 0)
	parameter (M_2 = 500)
	parameter (M_3 = 800)

	double precision Qtau, Qt, Qb
	parameter (Qtau = 0)
	parameter (Qt = 0)
	parameter (Qb = 0)

	double precision scalefactor
	parameter (scalefactor = 1)

	integer error

	call FHSetSMPara(error,
     &    invAlfa, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MB,
     &    MW, MZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)
	if( error .ne. 0 ) stop

	call FHSetPara(error, scalefactor,
     &    MT, TB, MA0, MHp,
     &    M3SL, M3SE, M3SQ, M3SU, M3SD,
     &    M2SL, M2SE, M2SQ, M2SU, M2SD,
     &    M1SL, M1SE, M1SQ, M1SU, M1SD,
     &    MUE,
     &    Atau, At, Ab,
     &    Amu, Ac, As,
     &    Ae, Au, Ad,
     &    M_1, M_2, M_3,
     &    Qtau, Qt, Qb)
        if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine setSLHA(filename)
	implicit none
	character*(*) filename

#include "SLHA.h"

	integer error
	double complex slhadata(nslhadata)

	call SLHARead(error, slhadata, filename, 1)
	if( error .ne. 0 ) stop

	call FHSetSLHA(error, slhadata)
	if( error .ne. 0 ) stop
	end

************************************************************************

	subroutine getPara
	implicit none

	integer nmfv
	double precision MSf(2,5,3), MASf(6,5), MCha(2), MNeu(4)
	double complex USf(2,2,5,3), UASf(6,6,5)
	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex DeltaMB
	double precision MGl
	double precision MHtree(4), SAtree

	integer error

	call FHGetPara(error, nmfv, MSf, USf, MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     &    MHtree, SAtree)
	if( error .ne. 0 ) stop

* print some sample output:
	print *, "MCha =", MCha
	print *, "MNeu =", MNeu
	end


************************************************************************

	subroutine higgsCorr
	implicit none

	integer error
	double precision MHiggs(4)
	double complex SAeff, UHiggs(3,3), ZHiggs(3,3)

	call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
	if( error .ne. 0 ) stop

* print some sample output:
	print *, "Mh1 =", MHiggs(1)
	print *, "Mh2 =", MHiggs(2)
	print *, "Mh3 =", MHiggs(3)
	print *, "MHp =", MHiggs(4)
	print *, "sin alpha_eff =", SAeff
	end

