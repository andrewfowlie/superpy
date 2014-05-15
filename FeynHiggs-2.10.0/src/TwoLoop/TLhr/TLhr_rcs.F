* TLhr_rcs.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_rcs
	implicit none

#include "TLhr.h"


	dMSfsq1133eps(-1) = 1/(3.D0*Pi)*(Alfas2L*(Cr(3) - Cr(6)))

	dMSfsq1133eps(0) = 1/(3.D0*Pi)*(Alfas2L*(Cr(12) - Cr(19)))

	dMSfsq1133eps(1) = -(1/(3.D0*Pi)*(Alfas2L*(Cr(21) + Cr(25))))

	dMSfsq1233eps(-1) = 1/(3.D0*Pi)*(Alfas2L*(Cr(29) - 2*Cr(31)))

	dMSfsq1233eps(0) = 1/(3.D0*Pi)*(Alfas2L*(Cr(35) - Cr(39)))

	dMSfsq1233eps(1) = 1/(3.D0*Pi)*(Alfas2L*(-Cr(40) + Cr(42)))

	dMSfsq1143eps(-1) = -(2/(3.D0*Pi)*(Alfas2L*Cr(64)))

	dMSfsq1143eps(0) = 1/(3.D0*Pi)*Cr(78)

	dMSfsq1143eps(1) = 1/(6.D0*Pi)*(Alfas2L*Cr(79))

	dMf133eps(-1) = -(4/(3.D0*Pi)*(Alfas2L*MTy))

        dMf133eps(0) = 1/(6.D0*Pi)*
     &   (-(Alfas2L*Cr(82)) + 6*Pi*(dMTfin*MTy))/MTy

	dMf133eps(1) = -(1/(12.D0*Pi)*(Alfas2L*Cr(83))/MTy)

        dY33eps(-1) = 4/(3.D0*Pi)*
     &   (Alfas2L*(2*MTy*Cr(138) - Cr(75)*MStop2(3)*UCStop(3,1)))

        dY33eps(0) = 2/(3.D0*Pi)*
     &   (Alfas2L*(MTy*Cr(67)*Cr(138) - 
     &       Cr(75)*Cr(139)*UCStop(3,1)))

	dY33eps(1) = 1/(3.D0*Pi)*(Alfas2L*Cr(75)*Cr(132)*UCStop(3,1))

        dAf133eps(-1) = (-(((-MUEC + At*TB)*dMf133eps(-1))/TB) + 
     &     (Conjugate(dMSfsq1133eps(-1)) - 
     &        Conjugate(dMSfsq1233eps(-1)))*UStop(1,1)*UStopC(1,2)+
     &       dY33eps(-1)*UStop(2,1)*UStopC(1,2) + 
     &     Conjugate(dY33eps(-1))*UStop(1,1)*UStopC(2,2))/MTy

        dAf133eps(0) = (-(((-MUEC + At*TB)*dMf133eps(0))/TB) + 
     &     (Conjugate(dMSfsq1133eps(0)) - 
     &        Conjugate(dMSfsq1233eps(0)))*UStop(1,1)*UStopC(1,2)+
     &       dY33eps(0)*UStop(2,1)*UStopC(1,2) + 
     &     Conjugate(dY33eps(0))*UStop(1,1)*UStopC(2,2))/MTy

        dAf133eps(1) = (-(((-MUEC + At*TB)*dMf133eps(1))/TB) + 
     &     (Conjugate(dMSfsq1133eps(1)) - 
     &        Conjugate(dMSfsq1233eps(1)))*UStop(1,1)*UStopC(1,2)+
     &       dY33eps(1)*UStop(2,1)*UStopC(1,2) + 
     &     Conjugate(dY33eps(1))*UStop(1,1)*UStopC(2,2))/MTy

	end
