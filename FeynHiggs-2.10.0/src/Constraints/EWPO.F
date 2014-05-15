* EWPO.F
* calculate electroweak precision observables
* this file is part of FeynHiggs
* last modified 26 Sep 12 th

#include "externals.h"
#include "types.h"


	subroutine EWPO(deltarho, MWMSSM, MWSM, SW2MSSM, SW2SM)
	implicit none
	RealType deltarho, MWMSSM, MWSM, SW2MSSM, SW2SM

#include "FH.h"

	RealType MH, gH, gX, dh, logdh, dt, dZZ, dAlfa, dAlfas

* MWSM from hep-ph/0311148, Eqs. 6-8
	RealType MH0_MW, MT0_MW, MZ0_MW, AlfasMZ0_MW
	parameter (MH0_MW = 100)
	parameter (MT0_MW = 174.3D0)
	parameter (MZ0_MW = 91.1875D0)
	parameter (AlfasMZ0_MW = .119D0)

	RealType MW0
	RealType c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11
	parameter (MW0 = 80.3799D0)
	parameter (c1  = -.05429D0)
	parameter (c2  = -.008939D0)
	parameter (c3  = .0000890D0)
	parameter (c4  = .000161D0)
	parameter (c5  = -1.070D0)
	parameter (c6  = .5256D0)
	parameter (c7  = -.0678D0)
	parameter (c8  = -.00179D0)
	parameter (c9  = .0000659D0)
	parameter (c10 = -.0737D0)
	parameter (c11 = 114.9D0)

* SW2SM from hep-ph/0407317, Eqs. 3-5
	RealType MT0_SW, MZ0_SW, AlfasMZ0_SW
	parameter (MT0_SW = 178)
	parameter (MZ0_SW = 91.1876D0)
	parameter (AlfasMZ0_SW = .117D0)

	RealType SW20
	RealType d1, d2, d3, d4, d5, d6, d7, d8, d9, d10
	parameter (SW20 = .2312527D0)
	parameter (d1   = 4.729D-4)
	parameter (d2   = 2.07D-5)
	parameter (d3   = 3.85D-6)
	parameter (d4   = -1.85D-6)
	parameter (d5   = .0207D0)
	parameter (d6   = -.002851D0)
	parameter (d7   = 1.82D-4)
	parameter (d8   = -9.74D-6)
	parameter (d9   = 3.98D-4)
	parameter (d10  = -.655D0)

	RealType DeltaRho_1L, DeltaRho_2L_gluon
	external DeltaRho_1L, DeltaRho_2L_gluon

* which Higgs is considered SM-like?

#if 0
* this choice valid for MA > 150 GeV
	MH = MHiggs(1)
#else
	MH = MHiggs(1)
	gH = abs(Coupling(H0VV(1,4)))
	gX = abs(Coupling(H0VV(2,4)))
	if( gX .gt. gH ) then
	  gH = gX
	  MH = MHiggs(2)
	endif
	if( abs(Coupling(H0VV(3,4))) .gt. gH ) MH = MHiggs(3)
#endif

* a) Delta Rho

	deltarho = DeltaRho_1L() + DeltaRho_2L_gluon()

* b) MW

	dh = MH/MH0_MW
	logdh = log(dh)
	dt = MT/MT0_MW
	dt = (dt - 1)*(dt + 1)
* Delta alpha_had is not input (yet)
	dAlfa = 0
	dAlfas = AlfasMZ/AlfasMZ0_MW - 1
	dZZ = MZ/MZ0_MW - 1

	MWSM = MW0 + logdh*(c1 + logdh*(c2 + logdh**2*c3)) +
     &    c4*(dh - 1)*(dh + 1) +
     &    c5*dAlfa +
     &    (c6 + c7*dt + c8*logdh + c9*dh**2)*dt +
     &    c10*dAlfas +
     &    c11*dZZ
	MWMSSM = MWSM*(1 + .5D0*CW2/(CW2 - SW2)*deltarho)

* c) SW2

	dt = MT/MT0_SW
	dt = (dt - 1)*(dt + 1)
	dAlfas = AlfasMZ/AlfasMZ0_SW - 1
	dZZ = MZ/MZ0_SW - 1

	SW2SM = SW20 + logdh*(d1 + logdh*(d2 + logdh**2*d3)) +
     &    d4*(dh - 1)*(dh + 1) +
     &    d5*dAlfa +
     &    (d6 + d7*dt + d8*(dh - 1))*dt +
     &    d9*dAlfas +
     &    d10*dZZ
	SW2MSSM = SW2SM - CW2*SW2/(CW2 - SW2)*deltarho
	end

