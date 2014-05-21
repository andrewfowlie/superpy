c This file is based on demo.F.
c SuperPy changes:
c Made program read name of SLHA file from command line.
c Made program calculate low-energy and flavour constraints only.
c Made program rewrite SLHA file with new Higgs masses.

* demo.F
* demonstration program for calling FeynHiggs from Fortran
* this file is part of FeynHiggs
* last modified 13 May 11 th

* compile this file with something like
* f77 -Ibuild example/demo.F -Lbuild -lFH


	program FHdemo
	implicit none
c	Added by SuperPy.
	character file*50
	call setFlags

c	Added by SuperPy.
	CALL GETARG(1, file)

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
c	call setPara
	call setSLHA(file)
c	Added by SuperPy.
	call flavour
c	Added by SuperPy.
	call couple
c	Added by SuperPy.
	call constraints
	call higgsCorr
c	Re-write SLHA file with new Higgs masses.
	call writeSLHA(file)
	end


************************************************************************s

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
c FHConstraints - evaluate electroweak precision observables.

	subroutine constraints

	implicit none

	integer error
	double precision gm2
	double precision Deltarho, MWMSSM, MWSM, SW2MSSM, SW2SM
	double precision edmeTh, edmn, edmHg
	double precision ccb

	call FHConstraints(error, gm2, Deltarho, MWMSSM, MWSM, SW2MSSM,
     &   SW2SM, edmeTh, edmn, edmHg, ccb)

	if( error .ne. 0 ) stop

* print some sample output:
	print *, "gm2=", gm2
	print *, "MWMSSM=", MWMSSM
	print *, "SW2MSSM=", SW2MSSM

	end

************************************************************************
c FHFlavour - evaluate flavour observables

	subroutine flavour

	implicit none

	integer error
	double precision bsgammaMSSM, bsgammaSM
	double precision deltaMsMSSM, deltaMsSM
	double precision bsmumuMSSM, bsmumuSM

	call FHFlavour(error, bsgammaMSSM, bsgammaSM, deltaMsMSSM,
     &	deltaMsSM, bsmumuMSSM, bsmumuSM)

	if( error .ne. 0 ) stop

* print some sample output:
	print *, "deltaMsMSSM=", deltaMsMSSM
	print *, "bsgammaMSSM=", bsgammaMSSM
	print *, "bsmumuMSSM=", bsmumuMSSM

	end

************************************************************************
c FHCouplings - compute Higgs couplings

	subroutine couple

#include "FHCouplings.h"

	integer error
c For FHCouplins
	double complex couplings(ncouplings), couplingsms(ncouplingsms)
	double precision gammas(ngammas), gammasms(ngammasms)
	integer fast
c For FHHiggsProd & 7 TeV and with fast turned on + neglect off-diagonals.
	double precision sqrts, prodxs(nprodxs)
	sqrts=7d0
	fast=1

	call FHCouplings(error,
     & couplings, couplingsms,
     & gammas, gammasms, fast)
	if( error .ne. 0 ) stop

	call FHHiggsProd(error, sqrts, prodxs)
	if( error .ne. 0 ) stop

c NB that this uses a clever pre-processor macro - see manual.
* print some sample output:
	print *, "g_ZZh in MSSM=", couplings(H0VV(1,3))
	print *, "g_ZZh in SM=", couplingsms(H0VV(1,3))
	print *, "R_ZZh MSSM/SM=", abs(couplings(H0VV(1,3)))/
     & abs(couplingsms(H0VV(1,3)))
	print *, "BR(h->gamma gamma)=", gammas(H0VV(1,1))
	print *, "inclusive sigma at 7TeV (fb)=", sum(prodxs)
	print *, "sigma * BR(h->gamma gamma)=", sum(prodxs)*
     & gammas(H0VV(1,1))
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

************************************************************************

	subroutine writeSLHA(filename)
	implicit none

#include "SLHA.h"

	integer error
	double complex slhadata(nslhadata)
	integer key
	character*(*) filename

	call SLHARead(error, slhadata, filename, 1)
	if( error .ne. 0 ) stop

	call FHSetSLHA(error, slhadata)
	if( error .ne. 0 ) stop

c   Add Higgs masses to SLHA file - key = 4.
	key = 4
	call FHOutputSLHA(error, slhadata, key)
	if( error .ne. 0 ) stop

	call SLHAWrite(error, slhadata, filename)
	if( error .ne. 0 ) stop

	end
