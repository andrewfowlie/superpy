* DeltaRho_1L.F
* generated 28-Feb-2013 11:34
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"
#include "types.h"
#include "debug.h"


	RealType function DeltaRho_1L()
	implicit none

#include "FH.h"
#include "looptools.h"

	ComplexType drho
	ComplexType tmp1, tmp2, tmp3, tmp4, tmp5
	integer All3, All4, Gen3, Ind1, Ind2, Sfe3, Sfe4

	drho = 0

	LOOP(Gen3, 1,3,1)

        drho = drho - 1/(2.D0*Pi**2*sqrt2)*
     &     (GF*B00(0.D0,MSf2(1,1,Gen3),MSf2(1,1,Gen3)))

	ENDLOOP(Gen3)

	LOOP(All3, 1,6,1)

	tmp1 = A0(MASf2(All3,3))

	tmp2 = A0(MASf2(All3,bBR))

	LOOP(Ind1, 1,3,1)

        drho = drho + 1/(3.D0*Pi**2*sqrt2)*
     &     (GF*SW2*(sqrt2**2*
     &          (tmp1*(-3*UASf(All3,Ind1,3)*UASfC(All3,Ind1,3) + 
     &              2*SW2*(UASf(All3,Ind1,3)*UASfC(All3,Ind1,3) + 
     &                 UASf(All3,3 + Ind1,3)*UASfC(All3,3 + Ind1,3)
     &                 ))) + 
     &         tmp2*(-3*UASf(All3,Ind1,bBR)*UASfC(All3,Ind1,bBR) + 
     &            SW2*(UASf(All3,Ind1,bBR)*UASfC(All3,Ind1,bBR) + 
     &               UASf(All3,3 + Ind1,bBR)*
     &                UASfC(All3,3 + Ind1,bBR)))))

	ENDLOOP(Ind1)

	ENDLOOP(All3)

	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        drho = drho + 1/(Pi**2*sqrt2)*
     &     (GF*(B00(0.D0,MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3))*
     &          USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     &         SW2*A0(MSf2(Sfe3,2,Gen3))*
     &          (-(CW2*USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     &            SW2*USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3))))

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        drho = drho - 1/(2.D0*Pi**2*sqrt2)*
     &     (GF*B00(0.D0,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &       (-(USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     &         2*SW2*(USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     &            USf(Sfe4,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)))*
     &       (-(USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3)) + 
     &         2*SW2*(USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     &            USf(Sfe3,2,2,Gen3)*USfC(Sfe4,2,2,Gen3))))

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

	LOOP(All4, 1,6,1)
	LOOP(All3, 1,6,1)

	tmp3 = B00(0.D0,MASf2(All3,3),MASf2(All4,3))

	tmp4 = B00(0.D0,MASf2(All3,3),MASf2(All4,bBR))

	tmp5 = B00(0.D0,MASf2(All3,bBR),MASf2(All4,bBR))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        drho = drho + 1/(Pi**2*sqrt2)*
     &     (GF*(1/6.D0*(-(tmp3*
     &               (-3*UASf(All4,Ind2,3)*UASfC(All3,Ind2,3) + 
     &                 SW2*
     &                  (4*UASf(All4,Ind2,3)*UASfC(All3,Ind2,3) + 
     &                    4*UASf(All4,3 + Ind2,3)*
     &                     UASfC(All3,3 + Ind2,3)))*
     &               (-3*UASf(All3,Ind1,3)*UASfC(All4,Ind1,3) + 
     &                 SW2*
     &                  (4*UASf(All3,Ind1,3)*UASfC(All4,Ind1,3) + 
     &                    4*UASf(All3,3 + Ind1,3)*
     &                     UASfC(All4,3 + Ind1,3)))) - 
     &            tmp5*(-3*UASf(All4,Ind2,bBR)*
     &                UASfC(All3,Ind2,bBR) + 
     &               SW2*(2*UASf(All4,Ind2,bBR)*
     &                   UASfC(All3,Ind2,bBR) + 
     &                  2*UASf(All4,3 + Ind2,bBR)*
     &                   UASfC(All3,3 + Ind2,bBR)))*
     &             (-3*UASf(All3,Ind1,bBR)*UASfC(All4,Ind1,bBR) + 
     &               SW2*(2*UASf(All3,Ind1,bBR)*
     &                   UASfC(All4,Ind1,bBR) + 
     &                  2*UASf(All3,3 + Ind1,bBR)*
     &                   UASfC(All4,3 + Ind1,bBR)))) + 
     &         3*tmp4*UASf(All3,Ind1,3)*UASf(All4,Ind2,bBR)*
     &          UASfC(All3,Ind2,3)*UASfC(All4,Ind1,bBR)))

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All3)
	ENDLOOP(All4)

	DeltaRho_1L = Re(drho)

	if( debuglevel .lt. 3 ) return

	DCONST 'DeltaRho_1L =', DeltaRho_1L ENDL
	end

