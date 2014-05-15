* Sfermions.F
* initialization routines for the sfermions
* this file is part of FeynHiggs
* last modified 3 Jun 13 th

#include "externals.h"
#include "types.h"
#include "debug.h"

#define MX2_LL(i,j) MX2(i,j)
#define MX2_LR(i,j) MX2(i,j+3)
#define MX2_RR(i,j) MX2(i+3,j+3)


	subroutine Sfermions(error)
	implicit none
	integer error

#include "FH.h"
#define __SUBROUTINE__ "Sfermions"

	integer i, j, g, fail, nfdummy, pass
	RealType Qew, Qqcd, asew, db2Lqcd, db2Lew, adb(3,3,2)
	RealType Mdin(3), Mdbar(3)
	ComplexType MX2(3,3), db, Mcorr(3,3), se(3,3)
	ComplexType DL(3,3), DR(3,3)

	RealType Deltab2L_ew, Deltab2L_qcd
	external Deltab2L_ew, Deltab2L_qcd

	error = 0
	CheckPara(error)
	if( error .ne. 0 ) return

	do i = 1, 3*3
	  CKM_flat(i) = CKMin_flat(i)
	enddo

* Sneutrino masses:
	do g = 1, 3
	  USf(1,1,1,g) = 1
	  USf(2,2,1,g) = 1
	  USf(2,1,1,g) = 0
	  USf(1,2,1,g) = 0
	enddo
	do i = 1, 6*6
	  UASf_flat(i,1) = 0
	enddo
	do i = 1, 6
	  UASf(i,i,1) = 1
	enddo

	MSf2(1,1,1) = Re(MSS2(1,1,1)) + DSf(1,1)
	MASf2(1,1) = MSf2(1,1,1)
	MSf(1,1,1) = sqrt(MSf2(1,1,1))
	MASf(1,1) = MSf(1,1,1)

	MSf2(1,1,2) = Re(MSS2(2,2,1)) + DSf(1,1)
	MASf2(2,1) = MSf2(1,1,2)
	MSf(1,1,2) = sqrt(MSf2(1,1,2))
	MASf(2,1) = MSf(1,1,2)

	MSf2(1,1,3) = Re(MSS2(3,3,1)) + DSf(1,1)
	MASf2(3,1) = MSf2(1,1,3)
	MSf(1,1,3) = sqrt(MSf2(1,1,3))
	MASf(3,1) = MSf(1,1,3)

	if( btest(fv, 2) ) then
	  MX2_LL(1,1) = MSf2(1,1,1)
	  MX2_LL(2,2) = MSf2(1,1,2)
	  MX2_LL(3,3) = MSf2(1,1,3)
	  MX2_LL(1,2) = MSS2(1,2,1)
	  MX2_LL(2,3) = MSS2(2,3,1)
	  MX2_LL(1,3) = MSS2(1,3,1)
	  call SfEigensystem(3, MX2,3, MASf2(1,1), UASf(1,1,1),6, 0, 3)
	endif

	do g = 1, 3
	  if( MSf2(1,1,g) .lt. 0 ) then
	    Error(error, "negative sneutrino mass squares")
	    MSf2(1,1,g) = 0
	  endif
	  MSf2(2,1,g) = 1D123
	  MSf(2,1,g) = 1D123
	  MASf2(g+3,1) = 1D123
	  MASf(g+3,1) = 1D123
	enddo

	if( debuglevel .ge. 1 ) then
	  DPARA "----------------------------"	ENDL
	  DPARA "MSf(1,1) =", MSf(1,1,1)	ENDL
	  DPARA "MSf(1,2) =", MSf(1,1,2)	ENDL
	  DPARA "MSf(1,3) =", MSf(1,1,3)	ENDL
	endif

	call SfUpdate(error, X3(2,2,2), 1, 1)
	call SfUpdate(error, X3(3,3,3), 1, 1)
	call SfUpdate(error, X3(4,4,4), 1, 1)
	call SfUpdate(error, X3(tT,3,3), SORT_SF, 1)
	if( error .ne. 0 ) return

	do i = 1, NSf
	  Sf(i,bTR) = Sf(i,4)
	enddo

	Mdin(1) = Mf(bTR,1)
	Mdin(2) = Mf(bTR,2)
	Mf(bTR,3) = MB_MT
	Mdin(3) = Mf(bTR,3)

	call GauginoMassesLargeTB
	fail = 0
	call SfUpdateLargeTB(fail)

	Qew = (MSf(1,3,3) + MSf(2,3,3) + abs(MUE))/3D0
	Qqcd = (MSdL(1,3) + MSdL(2,3) + MGl)/3D0
* for comparison with 2.7.4, for the moment:
	if( botResum .eq. 1 ) then
	  Qew = sqrt(MSdL(1,3)*MSdL(2,3))
	  Qqcd = Qew
	endif
	call AlphaS(asew, nfdummy, Qew)
	call AlphaS(AlfasDb, nfdummy, Qqcd)

	if( DEBUGLEVEL .gt. 3 ) then
	  DPARA "Delta b  Q(QCD)  =", Qqcd	ENDL
	  DPARA "Delta b  as(QCD) =", AlfasDb	ENDL
	  DPARA "Delta b  Q(ew)   =", Qew	ENDL
	  DPARA "Delta b  as(ew)  =", asew	ENDL
	endif

	Qqcd = Qqcd**2
	Qew = Qew**2

	adb(1,1,1) = AlfasDb
	adb(2,1,1) = AlfasDb
	adb(3,1,1) = AlfasDb
	adb(1,2,1) = AlfasDb
	adb(2,2,1) = AlfasDb
	adb(3,2,1) = AlfasDb    
	adb(1,3,1) = AlfasDb
	adb(2,3,1) = AlfasDb
	adb(3,3,1) = AlfasDb

	adb(1,1,2) = AlfaGF
	adb(2,1,2) = AlfaGF
	adb(3,1,2) = AlfaGF
	adb(1,2,2) = AlfaGF
	adb(2,2,2) = AlfaGF
	adb(3,2,2) = AlfaGF
	adb(1,3,2) = AlfaGF
	adb(2,3,2) = AlfaGF
	adb(3,3,2) = AlfaGF

	if( botResum .eq. 2 ) then
	  db2Lqcd = AlfasDb/pi*Deltab2L_qcd(Qqcd)
	  adb(3,3,1) = AlfasDb*(1 + db2Lqcd)
	  db2Lew = asew/pi*Deltab2L_ew(Qew)
	  adb(3,3,2) = AlfaGF*(1 + db2Lew)
	  if( debuglevel .ge. 2 ) then
	    DPARA "Delta b  f^2L_QCD =", db2Lqcd ENDL
	    DPARA "Delta b  f^2L_ew  =", db2Lew  ENDL
	  endif
	endif

	do g = 1, 3
	  call MDcorr(se, g, g, Qqcd, adb(g,g,1), Qew, adb(g,g,2))
	  Deltaf(bTR,g) = se(1,1)/Mf(bTR,g)
	  Mf(bTR,g) = Re(Mdin(g)/abs(1 + Deltaf(bTR,g))*
     &      (1 - .5D0*(se(2,1) + se(3,1))))
	enddo

	do pass = 1, 20
	  db = Deltaf(bTR,3)
	  if( debuglevel .ge. 2 ) then
	    DPARA "============================"	ENDL
	    DPARA "pass =", pass			ENDL
	    DPARA "Deltaf(bTR,3) =", Deltaf(bTR,3)	ENDL
	    DPARA "Mf(bTR,3) =", Mf(bTR,3)		ENDL
	  endif

	  call SfUpdateLargeTB(fail)
	  if( botResum .eq. 2 ) then
	    db2Lqcd = AlfasDb/pi*Deltab2L_qcd(Qqcd)
	    adb(3,3,1) = AlfasDb*(1 + db2Lqcd)
	    db2Lew = asew/pi*Deltab2L_ew(Qew)
	    adb(3,3,2) = AlfaGF*(1 + db2Lew)
	    if( debuglevel .ge. 2 ) then
	      DPARA "Delta b  f^2L_QCD =", db2Lqcd ENDL
	      DPARA "Delta b  f^2L_ew  =", db2Lew  ENDL
	    endif
	  endif

	  do g = 1, 3
	    do j = 1, 3
	      call MDcorr(se(1,j), j, g,
     &          Qqcd, adb(j,g,1), Qew, adb(j,g,2))
	      Mcorr(j,g) = se(1,j) +
     &          Mf(bTR,j)*se(2,j) + Mf(bTR,g)*se(3,j)
	    enddo
	    Deltaf(bTR,g) = se(1,g)/Mf(bTR,g)
	    Mf(bTR,g) = Re(Mdin(g)/abs(1 + Deltaf(bTR,g))*
     &        (1 - .5D0*(se(2,g) + se(3,g))))
	    Mcorr(g,g) = Mcorr(g,g) + Mf(bTR,g)
	  enddo

	  call SVD(3, 3, Mcorr,3, Mdbar, DR,3, DL,3, 0)

	  if( debuglevel .ge. 3 ) then
	    DPARA "MDc1 =", Mdbar(1), Mdin(1) ENDL
	    DPARA "MDc2 =", Mdbar(2), Mdin(2) ENDL
	    DPARA "MDc3 =", Mdbar(3), Mdin(3) ENDL
	  endif

	  do g = 1, 3
	    do j = 1, 3
	      CKM(j,g) = CKMin(j,1)*DL(1,g) +
     &                   CKMin(j,2)*DL(2,g) + 
     &                   CKMin(j,3)*DL(3,g)
	    enddo
c	    DL(g,g) = DL(g,g) - 1
c	    DR(g,g) = DR(g,g) - 1
	  enddo

	  if( abs(Deltaf(bTR,3) - db) .lt. 1D-6 ) goto 99
	enddo

	Warning("resummation of MB did not converge")

99	continue

* Delta_b with just the 1L corrections for subtraction in Hff:
	Deltaf(4,1) = Deltaf(bTR,1)
	Deltaf(4,2) = Deltaf(bTR,2)
	call MDcorr(se, 3, 3, Qqcd, AlfasDb, Qew, AlfaGF)
	Deltaf(4,3) = se(1,1)/Mf(bTR,3)

	if( fail .ne. 0 )
     &    Error(error, "negative sbottom mass squares in Delta b")

	call SfUpdate(error, X3(bTR,4,4), SORT_SF, 1)
c	call SfUpdate(error, X3(bTR,4,4), 1, 1)

	do i = 1, NSf
	  Sf(i,bTR0) = Sf(i,bTR)
	enddo

	Mf(bBR,1) = Mf(4,1)/abs(1 + Deltaf(bTR,1))
	Mf(bBR,2) = Mf(4,2)/abs(1 + Deltaf(bTR,2))
	Mf(bBR,3) = Mf(4,3)/abs(1 + Deltaf(bTR,3))
	call SfUpdate(error, X3(bBR,4,4), 1, 1)

	if( error .eq. 0 ) sf_valid = valid
	end

************************************************************************

	subroutine GauginoMassesLargeTB
	implicit none

#include "FH.h"

	ComplexType X(2,2), Y(4,4)

* chargino masses in large-TB (CB -> 0) limit
        X(1,1) = M_2
        X(1,2) = sqrt2*MW*SB
        X(2,1) = 0
        X(2,2) = MUE
        call SVD(2, 2, X,2, MChaL, UChaL,2, VChaL,2, 0)

* neutralino masses in large-TB (CB -> 0) limit
	Y(1,1) = M_1
	Y(2,2) = M_2
	Y(3,3) = 0
	Y(4,4) = 0
	Y(1,2) = 0
	Y(1,3) = 0
	Y(1,4) = MZ*SW*SB
	Y(2,3) = 0
	Y(2,4) = -MZ*CW*SB
	Y(3,4) = -MUE
	call TakagiFactor(4, Y,4, MNeuL, ZNeuL,4, 0)
	end

************************************************************************

	subroutine SfUpdateLargeTB(fail)
	implicit none
	integer fail

#include "FH.h"

	call EigenSf(fail, 1,
     &    MSdL(1,1), MSdL2(1,1), USdL(1,1,1),
     &    DSS2(1,4,1), Mf(bTR,1), MUETB(4), 0)
	call EigenSf(fail, 1,
     &    MSdL(1,2), MSdL2(1,2), USdL(1,1,2),
     &    DSS2(1,4,2), Mf(bTR,2), MUETB(4), 0)
	call EigenSf(fail, 1,
     &    MSdL(1,3), MSdL2(1,3), USdL(1,1,3),
     &    DSS2(1,4,3), Mf(bTR,3), MUETB(4), 0)

	if( debuglevel .lt. 2 ) return

	DPARA "----------------------------" ENDL
	DPARA "MSdL(1) =", MSdL(1,1), MSdL(2,1) ENDL
	DPARA "MSdL(2) =", MSdL(1,2), MSdL(2,2) ENDL
	DPARA "MSdL(3) =", MSdL(1,3), MSdL(2,3) ENDL
	end

************************************************************************

	subroutine SfUpdate(error, ty_ts_t, sort, show)
	implicit none
	integer error, ty_ts_t, sort, show

#include "FH.h"

	integer g, t, ty, fail
	character h
	character*32 msg(2:4)
	data msg /
     &    "negative stau mass squares",
     &    "negative stop mass squares",
     &    "negative sbottom mass squares" /

	fail = 0
	ty = nib2(ty_ts_t)
	t = nib0(ty_ts_t)

	Mf2(ty,1) = Mf(ty,1)**2
	Mf2(ty,2) = Mf(ty,2)**2
	Mf2(ty,3) = Mf(ty,3)**2

	Kf(1,1,ty) = Mf(ty,1)*Af(t,1)
	Kf(2,2,ty) = Mf(ty,2)*Af(t,2)
	Kf(3,3,ty) = Mf(ty,3)*Af(t,3)

	call SfUpdateMFV(fail, ty_ts_t, sort)
	do g = 1, 3
	  MASf(g,ty) = MSf(1,ty,g)
	  MASf(g+3,ty) = MSf(2,ty,g)
	  MASf2(g,ty) = MSf2(1,ty,g)
	  MASf2(g+3,ty) = MSf2(2,ty,g)
	  UASf(g,g,ty) = USf(1,1,ty,g)
	  USf2(1,1,ty,g) = Sq(USf(1,1,ty,g))
	  UASf(g+3,g,ty) = USf(2,1,ty,g)
	  USf2(2,1,ty,g) = Sq(USf(2,1,ty,g))
	  UASf(g,g+3,ty) = USf(1,2,ty,g)
	  USf2(1,2,ty,g) = Sq(USf(1,2,ty,g))
	  UASf(g+3,g+3,ty) = USf(2,2,ty,g)
	  USf2(2,2,ty,g) = Sq(USf(2,2,ty,g))
	enddo

	if( btest(fv, t) ) call SfUpdateNMFV(fail, ty_ts_t, sort)

	error = ior(error, fail)
	if( iand(fail, -show) .ne. 0 ) Error(error, msg(t))

	if( debuglevel .lt. 1 ) return

	h = Hex(ty)

	DPARA "----------------------------" ENDL
	DPARA "Mf("//h//",3)  =", Mf(ty,3) ENDL
	DPARA "MSf("//h//",1) =", MSf(1,ty,1), MSf(2,ty,1) ENDL
	DPARA "MSf("//h//",2) =", MSf(1,ty,2), MSf(2,ty,2) ENDL
	DPARA "MSf("//h//",3) =", MSf(1,ty,3), MSf(2,ty,3) ENDL

	if( debuglevel .lt. 2 ) return

	DPARA "USf(1,1,"//h//",3) =", USf(1,1,ty,3) ENDL
	DPARA "   (1,2)     =", USf(1,2,ty,3) ENDL
	DPARA "   (2,1)     =", USf(2,1,ty,3) ENDL
	DPARA "   (2,2)     =", USf(2,2,ty,3) ENDL
	DPARA "   11^2+12^2 =", USf2(1,1,ty,3) + USf2(1,2,ty,3) ENDL

	if( debuglevel .lt. 3 ) return

	DPARA "USf(1,1,"//h//",2) =", USf(1,1,ty,2) ENDL
	DPARA "   (1,2)     =", USf(1,2,ty,2) ENDL
	DPARA "   (2,1)     =", USf(2,1,ty,2) ENDL
	DPARA "   (2,2)     =", USf(2,2,ty,2) ENDL
	DPARA "   11^2+12^2 =", USf2(1,1,ty,2) + USf2(1,2,ty,2) ENDL

	DPARA "USf(1,1,"//h//",1) =", USf(1,1,ty,1) ENDL
	DPARA "   (1,2)     =", USf(1,2,ty,1) ENDL
	DPARA "   (2,1)     =", USf(2,1,ty,1) ENDL
	DPARA "   (2,2)     =", USf(2,2,ty,1) ENDL
	DPARA "   11^2+12^2 =", USf2(1,1,ty,1) + USf2(1,2,ty,1) ENDL
	end

************************************************************************

	subroutine SfUpdateMFV(fail, ty_ts_t, sort)
	implicit none
	integer fail, ty_ts_t, sort

#include "FH.h"

	integer ty, ts, t, i

	ty = nib2(ty_ts_t)
	ts = nib1(ty_ts_t)
	t = nib0(ty_ts_t)

#define CSf(i,t) i
	do i = UASf(1,1,*), UASf(6,6,*)
#undef CSf
	  CSf(i,ty) = 0
	enddo

	call EigenSf(fail, 1,
     &    MSf(1,ty,1), MSf2(1,ty,1), USf(1,1,ty,1),
     &    DSS2(1,ts,1), Mf(ty,1), Xf(t,1), sort)
	call EigenSf(fail, 1,
     &    MSf(1,ty,2), MSf2(1,ty,2), USf(1,1,ty,2),
     &    DSS2(1,ts,2), Mf(ty,2), Xf(t,2), sort)
	call EigenSf(fail, 1,
     &    MSf(1,ty,3), MSf2(1,ty,3), USf(1,1,ty,3),
     &    DSS2(1,ts,3), Mf(ty,3), Xf(t,3), sort)
	end

************************************************************************

	subroutine SfUpdateNMFV(fail, ty_ts_t, sort)
	implicit none
	integer fail, ty_ts_t, sort

#include "FH.h"

	integer t, ts, ty, nQ, nU, as, g, j
	ComplexType MX2(6,6), m1, m2, m3
	character h

	ty = nib2(ty_ts_t)
	ts = nib1(ty_ts_t)
	t = nib0(ty_ts_t)

	nQ = tQ(t)
	nU = tU(t)

	MX2_LL(1,1) = DSS2(1,ts,1) + DSS2(2,ts,1) + Mf2(ty,1)
	MX2_LL(2,2) = DSS2(1,ts,2) + DSS2(2,ts,2) + Mf2(ty,2)
	MX2_LL(3,3) = DSS2(1,ts,3) + DSS2(2,ts,3) + Mf2(ty,3)
	MX2_LL(1,2) = MSS2(1,2,nQ)
	MX2_LL(2,3) = MSS2(2,3,nQ)
	MX2_LL(1,3) = MSS2(1,3,nQ)

	MX2_RR(1,1) = DSS2(2,ts,1) + Mf2(ty,1)
	MX2_RR(2,2) = DSS2(2,ts,2) + Mf2(ty,2)
	MX2_RR(3,3) = DSS2(2,ts,3) + Mf2(ty,3)
	MX2_RR(1,2) = MSS2(1,2,nU)
	MX2_RR(2,3) = MSS2(2,3,nU)
	MX2_RR(1,3) = MSS2(1,3,nU)

	MX2_LR(1,1) = Mf(ty,1)*Xf(t,1)
	MX2_LR(2,2) = Mf(ty,2)*Xf(t,2)
	MX2_LR(3,3) = Mf(ty,3)*Xf(t,3)
	MX2_LR(1,2) = Kf(1,2,ty)
	MX2_LR(2,3) = Kf(2,3,ty)
	MX2_LR(1,3) = Kf(1,3,ty)
	MX2_LR(2,1) = Kf(2,1,ty)
	MX2_LR(3,2) = Kf(3,2,ty)
	MX2_LR(3,1) = Kf(3,1,ty)

	if( t .eq. 3 ) then
	  do g = 1, 3
	    m1 = MSS2(1,1,nQ)*CKMinC(g,1) +
     &           MSS2(1,2,nQ)*CKMinC(g,2) +
     &           MSS2(1,3,nQ)*CKMinC(g,3)
	    m2 = Conjugate(MSS2(1,2,nQ))*CKMinC(g,1) +
     &           MSS2(2,2,nQ)*CKMinC(g,2) +
     &           MSS2(2,3,nQ)*CKMinC(g,3)
	    m3 = Conjugate(MSS2(1,3,nQ))*CKMinC(g,1) +
     &           Conjugate(MSS2(2,3,nQ))*CKMinC(g,2) +
     &           MSS2(3,3,nQ)*CKMinC(g,3)
	    do j = 1, g
	      MX2_LL(j,g) = MX2_LL(j,g) +
     &          (CKMin(j,1)*m1 + CKMin(j,2)*m2 + CKMin(j,3)*m3 -
     &            MSS2(j,g,nQ))
	    enddo
	  enddo
	endif

	if( debuglevel .ge. 3 ) then
	  DPARA "============================"	ENDL
	  DPARA "MX2_LL(1,1) =", MX2_LL(1,1)	ENDL
	  DPARA "      (2,2) =", MX2_LL(2,2)	ENDL
	  DPARA "______(3,3) =", MX2_LL(3,3)	ENDL
	  DPARA "      (1,2) =", MX2_LL(1,2)	ENDL
	  DPARA "      (1,3) =", MX2_LL(1,3)	ENDL
	  DPARA "      (2,3) =", MX2_LL(2,3)	ENDL
	  DPARA "----------------------------"	ENDL
	  DPARA "MX2_LR(1,1) =", MX2_LR(1,1)	ENDL
	  DPARA "      (2,2) =", MX2_LR(2,2)	ENDL
	  DPARA "______(3,3) =", MX2_LR(3,3)	ENDL
	  DPARA "      (1,2) =", MX2_LR(1,2)	ENDL
	  DPARA "      (1,3) =", MX2_LR(1,3)	ENDL
	  DPARA "      (2,3) =", MX2_LR(2,3)	ENDL
	  DPARA "----------------------------"	ENDL
	  DPARA "MX2_RL(1,2) =", Conjugate(MX2_LR(2,1)) ENDL
	  DPARA "      (1,3) =", Conjugate(MX2_LR(3,1)) ENDL
	  DPARA "      (2,3) =", Conjugate(MX2_LR(3,2)) ENDL
	  DPARA "----------------------------"	ENDL
	  DPARA "MX2_RR(1,1) =", MX2_RR(1,1)	ENDL
	  DPARA "      (2,2) =", MX2_RR(2,2)	ENDL
	  DPARA "______(3,3) =", MX2_RR(3,3)	ENDL
	  DPARA "      (1,2) =", MX2_RR(1,2)	ENDL
	  DPARA "      (1,3) =", MX2_RR(1,3)	ENDL
	  DPARA "      (2,3) =", MX2_RR(2,3)	ENDL
        endif

	call SfEigensystem(6, MX2,6, MASf2(1,ty), UASf(1,1,ty),6, sort, 3)

	if( min(MASf2(1,ty), MASf2(2,ty), MASf2(3,ty),
     &          MASf2(4,ty), MASf2(5,ty), MASf2(6,ty)) .lt. 0 ) then
	  fail = 1
	  return
	endif

	do as = 1, 6
	  MASf(as,ty) = sqrt(MASf2(as,ty))
	enddo

	if( debuglevel .lt. 1 ) return

	h = Hex(ty)

	DPARA "----------------------------"	ENDL
	DPARA "MASf("//h//") =", MASf(1,ty), MASf(2,ty) ENDL
	DPARA "         ", MASf(3,ty), MASf(4,ty) ENDL
	DPARA "         ", MASf(5,ty), MASf(6,ty) ENDL

	if( debuglevel .lt. 2 ) return

	DPARA "UASf(1,1,"//h//") =", UASf(1,1,ty) ENDL
	DPARA "    (1,2)   =", UASf(1,2,ty)	ENDL
	DPARA "    (1,3)   =", UASf(1,3,ty)	ENDL
	DPARA "    (1,4)   =", UASf(1,4,ty)	ENDL
	DPARA "    (1,5)   =", UASf(1,5,ty)	ENDL
	DPARA "    (1,6)   =", UASf(1,6,ty)	ENDL
	DPARA "    (2,1)   =", UASf(2,1,ty)	ENDL
	DPARA "    (2,2)   =", UASf(2,2,ty)	ENDL
	DPARA "    (2,3)   =", UASf(2,3,ty)	ENDL
	DPARA "    (2,4)   =", UASf(2,4,ty)	ENDL
	DPARA "    (2,5)   =", UASf(2,5,ty)	ENDL
	DPARA "    (2,6)   =", UASf(2,6,ty)	ENDL
	DPARA "    (3,1)   =", UASf(3,1,ty)	ENDL
	DPARA "    (3,2)   =", UASf(3,2,ty)	ENDL
	DPARA "    (3,3)   =", UASf(3,3,ty)	ENDL
	DPARA "    (3,4)   =", UASf(3,4,ty)	ENDL
	DPARA "    (3,5)   =", UASf(3,5,ty)	ENDL
	DPARA "    (3,6)   =", UASf(3,6,ty)	ENDL
	DPARA "    (4,1)   =", UASf(4,1,ty)	ENDL
	DPARA "    (4,2)   =", UASf(4,2,ty)	ENDL
	DPARA "    (4,3)   =", UASf(4,3,ty)	ENDL
	DPARA "    (4,4)   =", UASf(4,4,ty)	ENDL
	DPARA "    (4,5)   =", UASf(4,5,ty)	ENDL
	DPARA "    (4,6)   =", UASf(4,6,ty)	ENDL
	DPARA "    (5,1)   =", UASf(5,1,ty)	ENDL
	DPARA "    (5,2)   =", UASf(5,2,ty)	ENDL
	DPARA "    (5,3)   =", UASf(5,3,ty)	ENDL
	DPARA "    (5,4)   =", UASf(5,4,ty)	ENDL
	DPARA "    (5,5)   =", UASf(5,5,ty)	ENDL
	DPARA "    (5,6)   =", UASf(5,6,ty)	ENDL
	DPARA "    (6,1)   =", UASf(6,1,ty)	ENDL
	DPARA "    (6,2)   =", UASf(6,2,ty)	ENDL
	DPARA "    (6,3)   =", UASf(6,3,ty)	ENDL
	DPARA "    (6,4)   =", UASf(6,4,ty)	ENDL
	DPARA "    (6,5)   =", UASf(6,5,ty)	ENDL
	DPARA "    (6,6)   =", UASf(6,6,ty)	ENDL
	end

************************************************************************

	RealType function Tdmb(x, y, z)
	implicit none
	RealType x, y, z

	RealType n, d

	n = 0
	d = 1

	if( x .ne. y ) then
	  if( x*y .ne. 0 ) n = n + x*y*log(x/y)
	  d = d*(x - y)
	endif
	if( x .ne. z ) then
	  if( x*z .ne. 0 ) n = n + x*z*log(z/x)
	  d = d*(x - z)
	endif
	if( y .ne. z ) then
	  if( y*z .ne. 0 ) n = n + y*z*log(y/z)
	  d = d*(y - z)
	endif

	Tdmb = n/d
	end

