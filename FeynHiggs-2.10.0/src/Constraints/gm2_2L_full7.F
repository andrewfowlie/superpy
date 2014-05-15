* gm2_2L_full7.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full7
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp30, tmp31, tmp32, tmp33, tmp34, tmp35, tmp36
	ComplexType tmp37, tmp38, tmp39, tmp40, tmp41, tmp42, tmp43
	ComplexType tmp44, tmp45, tmp46, tmp47, tmp48, tmp49, tmp50
	ComplexType tmp51

	Sf8619 = MStau2(2)*UStau22**2*(6 - 6*UStau22**2)

	Sf435 = -9*C2B*MW2*Sf407 + 22*MSl2Diff(1,1)**7 + 
     &   MW2*(-(Sf396*MSl2Diff(1,1)**5) + 
     &      MW2*(3*Sf398*MSl2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf397) + 3*Sf399) - 
     &               2*Sf401) - 2*Sf400*MSl2Diff(1,1))))

	Sf393 = 9*C2B*MW2*Sf363 - 22*MSl2Diff(1,1)**7 + 
     &   MW2*(Sf352*MSl2Diff(1,1)**5 + 
     &      MW2*(-3*Sf354*MSl2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf353) + 3*Sf355) - 
     &               2*Sf357) + 2*Sf356*MSl2Diff(1,1))))

	Sf1358 = -9*C2B*MW2*Sf1344 + 22*MSl2Diff(2,1)**7 + 
     &   MW2*(-(Sf1333*MSl2Diff(2,1)**5) + 
     &      MW2*(3*Sf1335*MSl2Diff(2,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf1334) + 3*Sf1336) - 
     &               2*Sf1338) - 2*Sf1337*MSl2Diff(2,1))))

	Sf1330 = 9*C2B*MW2*Sf1314 - 22*MSl2Diff(2,1)**7 + 
     &   MW2*(Sf1303*MSl2Diff(2,1)**5 + 
     &      MW2*(-3*Sf1305*MSl2Diff(2,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf1304) + 3*Sf1306) - 
     &               2*Sf1308) + 2*Sf1307*MSl2Diff(2,1))))

	Sf495 = 4*C2B*MW2**4*Sf483 - 3*MSl2Diff(1,1)**7 + 
     &   MW2*(3*Sf471*MSl2Diff(1,1)**5 + 
     &      MW2*(MW2*(-2*Sf485 + 
     &            MW2*(MW2*(MW2*Sf492 - 3*Sf493) + Sf494)) - 
     &         6*Sf472*MSl2Diff(1,1)**3))

	Sf1412 = 4*C2B*MW2**4*Sf1400 - 3*MSl2Diff(2,1)**7 + 
     &   MW2*(3*Sf1388*MSl2Diff(2,1)**5 + 
     &      MW2*(MW2*(-2*Sf1402 + 
     &            MW2*(MW2*(MW2*Sf1409 - 3*Sf1410) + Sf1411)) - 
     &         6*Sf1389*MSl2Diff(2,1)**3))

	Sf5415 = 180*MBot2*MTop2*MSbot2(2)**2 + 
     &   MSbot2(2)**2*((-9*Sf2881*Ss30 - 1412*MSbot2(2))*
     &       MSbot2(2) + Sf3012*MStop2(2)) + 
     &   MStop2(2)**2*(-9*Sf3013 + 
     &      MStop2(2)*(-Sf3014 + 1412*MStop2(2)))

	Sf3015 = 180*MBot2*MTop2*MSbot2(2)**2 + 
     &   MSbot2(2)**2*((-9*Sf2881*Ss30 - 1412*MSbot2(2))*
     &       MSbot2(2) + Sf3012*MStop2(1)) + 
     &   MStop2(1)**2*(-9*Sf3013 + 
     &      MStop2(1)*(-Sf3014 + 1412*MStop2(1)))

	Sf5309 = 180*MBot2*MTop2*MSbot2(1)**2 + 
     &   MSbot2(1)**2*((-9*Sf2881*Ss30 - 1412*MSbot2(1))*
     &       MSbot2(1) + Sf2883*MStop2(2)) + 
     &   MStop2(2)**2*(-9*Sf2885 + 
     &      MStop2(2)*(-Sf2887 + 1412*MStop2(2)))

	Sf2888 = 180*MBot2*MTop2*MSbot2(1)**2 + 
     &   MSbot2(1)**2*((-9*Sf2881*Ss30 - 1412*MSbot2(1))*
     &       MSbot2(1) + Sf2883*MStop2(1)) + 
     &   MStop2(1)**2*(-9*Sf2885 + 
     &      MStop2(1)*(-Sf2887 + 1412*MStop2(1)))

	tmp48 = MW2*(32*MW2 - 7*MSbot2(2)) - 22*MSbot2(2)**2

	tmp51 = -128*MSbot2(2)**4 + MW2*(MW2*tmp48 + 53*MSbot2(2)**3)

	Sf6288 = MW2*tmp51 + 72*MSbot2(2)**5 + Sf6287*MStop2(2)**2

	Sf3980 = MW2*tmp51 + 72*MSbot2(2)**5 + Sf3979*MStop2(1)**2

	tmp47 = MW2*(32*MW2 - 7*MSbot2(1)) - 22*MSbot2(1)**2

	tmp50 = -128*MSbot2(1)**4 + MW2*(MW2*tmp47 + 53*MSbot2(1)**3)

	Sf6274 = MW2*tmp50 + 72*MSbot2(1)**5 + Sf6273*MStop2(2)**2

	Sf3883 = MW2*tmp50 + 72*MSbot2(1)**5 + Sf3882*MStop2(1)**2

	Sf5429 = 20*MTop2*MSbot2(2)**2 + 
     &   MSbot2(2)**2*(-(Sf2881*MSbot2(2)) + 
     &      27*Sf2882*MStop2(2)) + 
     &   MStop2(2)**2*(-Sf3032 - Sf2886*MStop2(2))

	Sf3033 = 20*MTop2*MSbot2(2)**2 + 
     &   MSbot2(2)**2*(-(Sf2881*MSbot2(2)) + 
     &      27*Sf2882*MStop2(1)) + 
     &   MStop2(1)**2*(-Sf3032 - Sf2886*MStop2(1))

	Sf5323 = 20*MTop2*MSbot2(1)**2 + 
     &   MSbot2(1)**2*(-(Sf2881*MSbot2(1)) + 
     &      27*Sf2882*MStop2(2)) + 
     &   MStop2(2)**2*(-Sf2905 - Sf2886*MStop2(2))

	Sf2906 = 20*MTop2*MSbot2(1)**2 + 
     &   MSbot2(1)**2*(-(Sf2881*MSbot2(1)) + 
     &      27*Sf2882*MStop2(1)) + 
     &   MStop2(1)**2*(-Sf2905 - Sf2886*MStop2(1))

	Sf1814 = 3*Sf1813 - Sf1810*MSneu2(1) + 
     &   MSneu2(1)**2*(-3*Sf1809 + 
     &      (30*Sf1594 - 15*MSneu2(1))*MSneu2(1))

	Sf1003 = 3*Sf1002 - Sf999*MSneu2(1) + 
     &   MSneu2(1)**2*(-3*Sf997 + 
     &      (30*Sf734 - 15*MSneu2(1))*MSneu2(1))

	Sf10034 = MZ2**4*MHiggs2(2) + Sf10033*MStau2(2)**2 + 
     &   MZ2*(-2*Sf9937*MStau2(2)**3 + 
     &      MZ2**2*(-MZ2**2 - 36*Sf9024*MStau2(2)))

	Sf9940 = MZ2**4*MHiggs2(2) + Sf9939*MStau2(1)**2 + 
     &   MZ2*(-2*Sf9937*MStau2(1)**3 + 
     &      MZ2**2*(-MZ2**2 - 36*Sf9024*MStau2(1)))

	Sf2136 = Sf2135*MSneu2(1) + 
     &   MW2*(-4*MW2**2 + Sf733*MStau2(2)) + 
     &   MHiggs2(4)*(7*MW2**2 + 3*MStau2(2)**2)

	Sf2056 = Sf2055*MSneu2(1) + 
     &   MW2*(-4*MW2**2 + Sf733*MStau2(1)) + 
     &   MHiggs2(4)*(7*MW2**2 + 3*MStau2(1)**2)

	Sf736 = -(Sf735*MSneu2(1)) + 
     &   MW2*(-MW2**2 + Sf733*MStau2(1)) + 
     &   MHiggs2(4)*(2*MW2**2 - 3*MStau2(1)**2)

	Sf1596 = -(Sf1595*MSneu2(1)) + 
     &   MW2*(-MW2**2 + Sf733*MStau2(2)) + 
     &   MHiggs2(4)*(2*MW2**2 - 3*MStau2(2)**2)

	tmp43 = 2*MW2**2 - 3*MSbot2(2)**2

	Sf6844 = tmp43*MHiggs2(4) + 
     &   MW2*(-MW2**2 + Sf733*MSbot2(2)) - Sf6729*MStop2(2)

	Sf6730 = tmp43*MHiggs2(4) + 
     &   MW2*(-MW2**2 + Sf733*MSbot2(2)) - Sf6729*MStop2(1)

	tmp42 = 2*MW2**2 - 3*MSbot2(1)**2

	Sf6589 = tmp42*MHiggs2(4) + 
     &   MW2*(-MW2**2 + Sf733*MSbot2(1)) - Sf6462*MStop2(2)

	Sf6463 = tmp42*MHiggs2(4) + 
     &   MW2*(-MW2**2 + Sf733*MSbot2(1)) - Sf6462*MStop2(1)

	tmp45 = -177*MZ2**2 + MHiggs2(2)*(355*MZ2 + 56*MHiggs2(2))

	Sf9723 = -83*MZ2**3 + tmp45*MHiggs2(2) - Sf9722*MStau2(2)

	Sf9466 = -83*MZ2**3 + tmp45*MHiggs2(2) - Sf9465*MStau2(1)

	tmp44 = -9*MZ2**2 + MHiggs2(2)*(5*MZ2 + 16*MHiggs2(2))

	Sf9752 = -22*MZ2**3 + tmp44*MHiggs2(2) + Sf9751*MStau2(2)

	Sf9525 = -22*MZ2**3 + tmp44*MHiggs2(2) + Sf9524*MStau2(1)

	Sf5692 = MTop2*(6*MTop2 + MW2*Sf5689) - MW2**2*Sf5690 + 
     &   MBot2*(6*MBot2 - Sf5691)

	Sf5667 = MTop2*(3*MTop2 + 3*MW2*Sf5664) - MW2**2*Sf5665 + 
     &   MBot2*(3*MBot2 - 3*Sf5666)

	Sf3866 = MTop2*(3*MTop2 - MW2*Sf3863) - MW2**2*Sf3864 + 
     &   MBot2*(3*MBot2 - Sf3865)

	Sf3370 = MBot2*(3*MBot2 - Sf3369) + 
     &   MTop2*(3*MTop2 - MW2*Sf504) - C2B*MW2**2*Sf601

	Sf5682 = MTop2*(3*MTop2 + 2*MW2*Sf169) - MW2**2*Sf5680 + 
     &   MBot2*(3*MBot2 - 2*Sf5681)

	Sf3350 = MTop2*(2*MTop2 + 4*MW2*Sf174) + 
     &   MBot2*(2*MBot2 - 4*Sf2238) + MW2**2*Sf3349

	Sf3853 = C2B*MW2**2*Sf174 + MBot2*(MBot2 - Sf3852) + 
     &   MTop2*(MTop2 + MW2*Sf566)

	Sf10593 = Sf10592*MStop2(1) + 
     &   MStop2(1)**2*(Sf10591 + MHiggs2(2)*MStop2(1)) - 
     &   Sf10245*Sf10246*Sf10472*MStop2(2)

	Sf10332 = Sf10331*MStop2(1) + 
     &   MStop2(1)**2*(Sf10329 + MHiggs2(1)*MStop2(1)) - 
     &   Sf10245*Sf10246*Sf10330*MStop2(2)

	tmp49 = 2*MZ2**2*MStau2(1) + 
     &   MZ2*(3*MStau2(1) + MStau2(2))*MStau2(3) - MStau2(3)**3

	Sf9869 = tmp49*MHiggs2(2) + 
     &   MZ2**2*(MStau2(3)**2 - MZ2*MStau2(4))

	Sf10218 = tmp49*MHiggs2(1) + 
     &   MZ2**2*(MStau2(3)**2 - MZ2*MStau2(4))

	Sf9682 = 64*MStau2(1)**3 + 
     &   MZ2*(37*MZ2*MStau2(1) - 12*MStau2(1)**2) + 
     &   Sf9681*MStau2(2)

	Sf9652 = 16*MStau2(1)**3 + 
     &   MZ2*(7*MZ2*MStau2(1) + 6*MStau2(1)**2) + Sf9651*MStau2(2)

	Sf8412 = Sf8395*MSq2Diff(2,2)**2 + Sf8402*MStop2(2) + 
     &   Sf4906*(Sf8388*MSq2Diff(2,2)**2 + Sf8383*MStop2(2))

	Sf8316 = Sf8299*MSq2Diff(1,2)**2 + Sf8306*MStop2(2) + 
     &   Sf4820*(Sf8292*MSq2Diff(1,2)**2 + Sf8287*MStop2(2))

	Sf6895 = Sf4696*MSbot2(2)**2 + 3*Sf4697*MStop2(2) + 
     &   (Sf4698 - 260*MStop2(2))*MStop2(2)**2

	Sf4699 = Sf4696*MSbot2(2)**2 + 3*Sf4697*MStop2(1) + 
     &   (Sf4698 - 260*MStop2(1))*MStop2(1)**2

	Sf5428 = 10*MTop2*MSbot2(2)**2 + Sf5427*MStop2(2)**2 + 
     &   MSbot2(2)*(-123*MSbot2(2)**2 - Sf3028*MStop2(2))

	Sf3031 = 10*MTop2*MSbot2(2)**2 + Sf3030*MStop2(1)**2 + 
     &   MSbot2(2)*(-123*MSbot2(2)**2 - Sf3028*MStop2(1))

	Sf6910 = 10*MTop2*MSbot2(2)**2 + Sf6909*MStop2(2)**2 + 
     &   MSbot2(2)*(-227*MSbot2(2)**2 + Sf4714*MStop2(2))

	Sf4717 = 10*MTop2*MSbot2(2)**2 + Sf4716*MStop2(1)**2 + 
     &   MSbot2(2)*(-227*MSbot2(2)**2 + Sf4714*MStop2(1))

	Sf5435 = 10*MBot2*MSbot2(2)**2 + Sf5434*MStop2(2)**2 + 
     &   MSbot2(2)*(123*MSbot2(2)**2 - Sf3039*MStop2(2))

	Sf3042 = 10*MBot2*MSbot2(2)**2 + Sf3041*MStop2(1)**2 + 
     &   MSbot2(2)*(123*MSbot2(2)**2 - Sf3039*MStop2(1))

	Sf6644 = Sf4597*MSbot2(1)**2 + 3*Sf4599*MStop2(2) + 
     &   (Sf4601 - 260*MStop2(2))*MStop2(2)**2

	Sf4602 = Sf4597*MSbot2(1)**2 + 3*Sf4599*MStop2(1) + 
     &   (Sf4601 - 260*MStop2(1))*MStop2(1)**2

	Sf5322 = 10*MTop2*MSbot2(1)**2 + Sf5321*MStop2(2)**2 + 
     &   MSbot2(1)*(-123*MSbot2(1)**2 - Sf2901*MStop2(2))

	Sf2904 = 10*MTop2*MSbot2(1)**2 + Sf2903*MStop2(1)**2 + 
     &   MSbot2(1)*(-123*MSbot2(1)**2 - Sf2901*MStop2(1))

	Sf6659 = 10*MTop2*MSbot2(1)**2 + Sf6658*MStop2(2)**2 + 
     &   MSbot2(1)*(-227*MSbot2(1)**2 + Sf4619*MStop2(2))

	Sf4622 = 10*MTop2*MSbot2(1)**2 + Sf4621*MStop2(1)**2 + 
     &   MSbot2(1)*(-227*MSbot2(1)**2 + Sf4619*MStop2(1))

	Sf5329 = 10*MBot2*MSbot2(1)**2 + Sf5328*MStop2(2)**2 + 
     &   MSbot2(1)*(123*MSbot2(1)**2 - Sf2912*MStop2(2))

	Sf2915 = 10*MBot2*MSbot2(1)**2 + Sf2914*MStop2(1)**2 + 
     &   MSbot2(1)*(123*MSbot2(1)**2 - Sf2912*MStop2(1))

	Sf9286 = Sf9285*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf9284 + MHiggs2(1)*MSbot2(1)) - 
     &   Sf8885*Sf8886*Sf9149*MSbot2(2)

	Sf9107 = Sf9106*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf9105 + MHiggs2(2)*MSbot2(1)) - 
     &   Sf8885*Sf8886*Sf8959*MSbot2(2)

	Sf10258 = 2*Sf8681 + 1/2.D0*Sf8678 + (-1 - 2*Ss988)*MStop2(2)

	Sf10238 = 2*Sf7995 + 1/2.D0*Sf7992 + (-1 - 2*Ss969)*MStop2(1)

	Sf5996 = 2*Sf5574 - 2*Sf5570*MSbot2(2) + 
     &   (5*Sf5507 - MSbot2(2))*MSbot2(2)**2

	Sf5577 = 2*Sf5574 - 2*Sf5570*MSbot2(1) + 
     &   (5*Sf5507 - MSbot2(1))*MSbot2(1)**2

	Sf3591 = 2*Sf3252 - 2*Sf3247*MSbot2(2) + 
     &   (5*Sf3128 - MSbot2(2))*MSbot2(2)**2

	Sf3257 = 2*Sf3252 - 2*Sf3247*MSbot2(1) + 
     &   (5*Sf3128 - MSbot2(1))*MSbot2(1)**2

	Sf9370 = -(Sf9369*MHiggs2(2)) + 
     &   MZ2*(MZ2**2 + Ss143*MStau2(1)) + Sf9319*MStau2(2)

	Sf9371 = -(Sf9369*MHiggs2(2)) + Sf9350*MStau2(1) + 
     &   MZ2*(MZ2**2 + Ss143*MStau2(2))

	Sf10111 = MZ2**2*(MZ2 - MStau2(3)) - 
     &   MHiggs2(1)*(2*MZ2**2 - MZ2*(MStau2(1) + 3*MStau2(2)) + 
     &      MStau2(3)**2)

	Sf10112 = (MZ2 - MStau2(1))*
     &    (MZ2*(MZ2 - 2*MHiggs2(1)) + MHiggs2(1)*MStau2(1)) + 
     &   (MZ2*(MZ2 + MHiggs2(1)) + 2*MHiggs2(1)*MStau2(1))*
     &    MStau2(2) - MHiggs2(1)*MStau2(2)**2

	Sf8987 = -(Sf8986*MHiggs2(2)) + 
     &   MZ2*(MZ2**2 + Ss143*MSbot2(1)) + Sf8906*MSbot2(2)

	Sf8988 = -(Sf8986*MHiggs2(2)) + Sf8958*MSbot2(1) + 
     &   MZ2*(MZ2**2 + Ss143*MSbot2(2))

	Sf9174 = -(Sf8986*MHiggs2(1)) + 
     &   MZ2*(MZ2**2 + Ss145*MSbot2(1)) + Sf9110*MSbot2(2)

	Sf9175 = -(Sf8986*MHiggs2(1)) + Sf9148*MSbot2(1) + 
     &   MZ2*(MZ2**2 + Ss145*MSbot2(2))

	Sf10489 = -(Sf10383*MHiggs2(2)) + 
     &   MZ2*(MZ2**2 + Ss143*MStop2(1)) + Sf10438*MStop2(2)

	Sf10490 = -(Sf10383*MHiggs2(2)) + Sf10471*MStop2(1) + 
     &   MZ2*(MZ2**2 + Ss143*MStop2(2))

	Sf10384 = -(Sf10383*MHiggs2(1)) + 
     &   MZ2*(MZ2**2 + Ss145*MStop2(1)) + Sf10344*MStop2(2)

	Sf10385 = -(Sf10383*MHiggs2(1)) + Sf10371*MStop2(1) + 
     &   MZ2*(MZ2**2 + Ss145*MStop2(2))

	Sf189 = 4*ML2*MW2**2 - MW2**3*Sf168 + 
     &   MHiggs2(4)*(-3*MW2*Sf185 + 2*Sf176*MHiggs2(4))

	Sf1012 = 2*ML2*MW2**2 - 2*MW2**3*Sf243 + 
     &   MHiggs2(4)*(-3*MW2*Sf1010 + Sf1011*MHiggs2(4))

	Sf5091 = MW2*(7*MW2 - 9*Sf99) - 5*Sf96*MHiggs2(4) + 
     &   2*MSq2Diff(2,2)**2

	Sf4965 = MW2*(7*MW2 - 9*Sf91) - 5*Sf88*MHiggs2(4) + 
     &   2*MSq2Diff(1,2)**2

	Sf2651 = MW2*(7*MW2 - 9*Sf77) - 5*Sf74*MHiggs2(4) + 
     &   2*MSq2Diff(2,1)**2

	Sf2437 = MW2*(7*MW2 - 9*Sf69) - 5*Sf66*MHiggs2(4) + 
     &   2*MSq2Diff(1,1)**2

	Sf1272 = MW2*(7*MW2 - 9*Sf47) - 5*Sf44*MHiggs2(4) + 
     &   2*MSl2Diff(2,1)**2

	Sf312 = MW2*(7*MW2 - 9*Sf34) - 5*Sf31*MHiggs2(4) + 
     &   2*MSl2Diff(1,1)**2

	Sf6229 = MW2*(3*MW2 - 5*Sf99) - Sf96*MHiggs2(4) + 
     &   2*MSq2Diff(2,2)**2

	Sf5863 = MW2*(3*MW2 - 5*Sf91) - Sf88*MHiggs2(4) + 
     &   2*MSq2Diff(1,2)**2

	Sf4170 = MW2*(3*MW2 - 5*Sf77) - Sf74*MHiggs2(4) + 
     &   2*MSq2Diff(2,1)**2

	Sf4085 = MW2*(3*MW2 - 5*Sf69) - Sf66*MHiggs2(4) + 
     &   2*MSq2Diff(1,1)**2

	Sf1569 = MW2*(MW2 + Sf47) - 3*Sf44*MHiggs2(4) - 
     &   2*MSl2Diff(2,1)**2

	Sf701 = MW2*(MW2 + Sf34) - 3*Sf31*MHiggs2(4) - 
     &   2*MSl2Diff(1,1)**2

	tmp32 = MW2*(-33*MW2 + 12*MHiggs2(4))

	Sf1668 = 25*MW2**3 + tmp32*MHiggs2(4) - Sf1667*MStau2(2)

	Sf852 = 25*MW2**3 + tmp32*MHiggs2(4) - Sf851*MStau2(1)

	Sf3606 = 25*MW2**3 + tmp32*MHiggs2(4) - Sf3605*MSbot2(2)

	Sf3279 = 25*MW2**3 + tmp32*MHiggs2(4) - Sf3278*MSbot2(1)

	Sf1450 = -9*MW2**3 + 
     &   MW2*MHiggs2(4)*(MW2 + 4*MHiggs2(4)) + Sf1449*MStau2(2)

	Sf547 = -9*MW2**3 + 
     &   MW2*MHiggs2(4)*(MW2 + 4*MHiggs2(4)) + Sf546*MStau2(1)

	Sf3522 = -9*MW2**3 + 
     &   MW2*MHiggs2(4)*(MW2 + 4*MHiggs2(4)) + Sf3521*MSbot2(2)

	Sf3160 = -9*MW2**3 + 
     &   MW2*MHiggs2(4)*(MW2 + 4*MHiggs2(4)) + Sf3159*MSbot2(1)

	Sf2286 = 4*MW2**2 + 
     &   MHiggs2(4)*(-Sf2223 + 9*MHiggs2(4)) + 3*Sf2201*MSbot2(2)

	Sf2224 = 4*MW2**2 + 
     &   MHiggs2(4)*(-Sf2223 + 9*MHiggs2(4)) + 3*Sf2201*MSbot2(1)

	Sf1666 = 10*MW2**2 + 
     &   MHiggs2(4)*(-15*MW2 + 6*MHiggs2(4)) + Sf1445*MStau2(2)

	Sf849 = 10*MW2**2 + 
     &   MHiggs2(4)*(-15*MW2 + 6*MHiggs2(4)) + Sf541*MStau2(1)

	Sf3604 = 10*MW2**2 + 
     &   MHiggs2(4)*(-15*MW2 + 6*MHiggs2(4)) + Sf3519*MSbot2(2)

	Sf3277 = 10*MW2**2 + 
     &   MHiggs2(4)*(-15*MW2 + 6*MHiggs2(4)) + Sf3157*MSbot2(1)

	Sf1364 = 2*MW2**2 + MHiggs2(4)*(-2*Sf3 + 3*MHiggs2(4)) + 
     &   3*MStau2(2)**2

	Sf446 = 2*MW2**2 + MHiggs2(4)*(-2*Sf3 + 3*MHiggs2(4)) + 
     &   3*MStau2(1)**2

	Sf5060 = MHiggs2(4)*(2*MHiggs2(4) - MSbot2(2)) - 
     &   11*MSbot2(2)**2 + Sf5059*MStop2(2)

	Sf2607 = MHiggs2(4)*(2*MHiggs2(4) - MSbot2(2)) - 
     &   11*MSbot2(2)**2 + Sf2606*MStop2(1)

	Sf4934 = MHiggs2(4)*(2*MHiggs2(4) - MSbot2(1)) - 
     &   11*MSbot2(1)**2 + Sf4933*MStop2(2)

	Sf2369 = MHiggs2(4)*(2*MHiggs2(4) - MSbot2(1)) - 
     &   11*MSbot2(1)**2 + Sf2368*MStop2(1)

	Sf1448 = -13*MW2**2 + 
     &   MHiggs2(4)*(5*MW2 + 2*MHiggs2(4)) + Sf1447*MStau2(2)

	Sf544 = -13*MW2**2 + MHiggs2(4)*(5*MW2 + 2*MHiggs2(4)) + 
     &   Sf543*MStau2(1)

	Sf3524 = -13*MW2**2 + 
     &   MHiggs2(4)*(5*MW2 + 2*MHiggs2(4)) + Sf3523*MSbot2(2)

	Sf3162 = -13*MW2**2 + 
     &   MHiggs2(4)*(5*MW2 + 2*MHiggs2(4)) + Sf3161*MSbot2(1)

	Sf1446 = -2*MW2**2 + MHiggs2(4)*(-MW2 + 2*MHiggs2(4)) - 
     &   Sf1445*MStau2(2)

	Sf542 = -2*MW2**2 + MHiggs2(4)*(-MW2 + 2*MHiggs2(4)) - 
     &   Sf541*MStau2(1)

	Sf3520 = -2*MW2**2 + MHiggs2(4)*(-MW2 + 2*MHiggs2(4)) - 
     &   Sf3519*MSbot2(2)

	Sf3158 = -2*MW2**2 + MHiggs2(4)*(-MW2 + 2*MHiggs2(4)) - 
     &   Sf3157*MSbot2(1)

	Sf2239 = 2*MW2**2 + 3*MBot2*MHiggs2(4) - 
     &   3*Sf2238*MHiggs2(4) + MHiggs2(4)**2

	Sf8466 = (Sf8465*UStop21**2)/Sf4820**3 + 
     &   ((Sf2320*Sf8232 + Sf8239)*UStop22**2)/Sf2320**3

	Sf8464 = (Sf8463*UStop21**2)/Sf4906**3 + 
     &   ((Sf2334*Sf8203 + Sf8210)*UStop22**2)/Sf2334**3

	Sf8545 = MZ2**2*(2*MZ2*Sf8542 + 21*Sf8544) + 7*Sf8541*UStau11

	Sf8580 = MZ2**2*(-(MZ2*Sf8578) + 6*Sf8579) + 2*Sf8541*UStau11

	Sf8462 = (Sf8461*USbot21**2)/Sf4906**3 + 
     &   ((Sf4820*Sf8160 - Sf8173)*USbot22**2)/Sf4820**3

	Sf8438 = (Sf8437*USbot21**2)/Sf2334**3 + 
     &   ((Sf2320*Sf7951 - Sf7964)*USbot22**2)/Sf2320**3

	Sf10572 = MHiggs2(2)*MStop2(1)**2 - MStop2(1)**3 - 
     &   Sf10571*MStop2(2)**2

	Sf10422 = MHiggs2(1)*MStop2(1)**2 - MStop2(1)**3 - 
     &   Sf10421*MStop2(2)**2

	Sf10462 = 4*MZ2*MStop2(1) - MStop2(1)**2 + Sf10461*MStop2(2)

	Sf9847 = MHiggs2(2)*MStau2(1)**2 - MStau2(1)**3 - 
     &   Sf9846*MStau2(2)**2

	Sf10197 = MHiggs2(1)*MStau2(1)**2 - MStau2(1)**3 - 
     &   Sf10196*MStau2(2)**2

	Sf594 = Sf593*MSneu2(1)**2 + Sf118*MSneu2(1)**3 + 
     &   Sf592*MStau2(1)**2

	Sf298 = Sf297*MSneu2(1)**2 - Sf292*MSneu2(1)**3 - 
     &   Sf295*MStau2(1)**2

	Sf282 = Sf280*MSneu2(1)**2 - Sf276*MSneu2(1)**3 - 
     &   Sf281*MStau2(1)**2

	Sf223 = Sf222*MSneu2(1)**2 - Sf219*MSneu2(1)**3 - 
     &   Sf167**2*Sf220*MStau2(1)

	Sf1494 = Sf1493*MSneu2(1)**2 + Sf118*MSneu2(1)**3 + 
     &   Sf1492*MStau2(2)**2

	Sf1261 = Sf1260*MSneu2(1)**2 - Sf292*MSneu2(1)**3 - 
     &   Sf1259*MStau2(2)**2

	Sf1252 = Sf1250*MSneu2(1)**2 - Sf276*MSneu2(1)**3 - 
     &   Sf1251*MStau2(2)**2

	Sf1215 = Sf1214*MSneu2(1)**2 - Sf1212*MSneu2(1)**3 - 
     &   Sf1196**2*Sf1213*MStau2(2)

	Sf751 = Sf748*MSneu2(1) - Sf745*MSneu2(1)**2 + Sf750*MStau2(1)

	Sf7461 = -MSneu2(1) + Sf7457*MStau2(1) + Sf7460*MStau2(2)

	Sf728 = Sf726*MSneu2(1) - 2*C2B*MSneu2(1)**2 + Sf727*MStau2(1)

	Sf725 = Sf723*MSneu2(1) + Sf243*MSneu2(1)**2 + Sf724*MStau2(1)

	Sf718 = Sf716*MSneu2(1) + Sf299*MSneu2(1)**2 + Sf717*MStau2(1)

	Sf694 = Sf691*MSneu2(1) + Sf690*MSneu2(1)**2 + Sf693*MStau2(1)

	Sf687 = Sf685*MSneu2(1) - Sf684*MSneu2(1)**2 - 3*Sf686*MStau2(1)

	Sf678 = Sf676*MSneu2(1) - Sf415*MSneu2(1)**2 - Sf677*MStau2(1)

	Sf591 = Sf589*MSneu2(1) + Sf300*MSneu2(1)**2 + Sf590*MStau2(1)

	Sf530 = Sf529*MSneu2(1) + Sf526*MSneu2(1)**2 + Sf527*MStau2(1)

	Sf304 = Sf303*MSneu2(1) + 13*ML2*MStau2(1) - Sf300*MStau2(1)**2

	Sf291 = Sf289*MSneu2(1) - Sf288*MSneu2(1)**2 - Sf290*MStau2(1)

	Sf275 = Sf273*MSneu2(1) + Sf272*MSneu2(1)**2 - Sf274*MStau2(1)

	Sf2121 = Sf2119*MSneu2(1) + Sf299*MSneu2(1)**2 + 
     &   Sf2120*MStau2(2)

	Sf2118 = Sf2116*MSneu2(1) + Sf2033*MSneu2(1)**2 + 
     &   3*Sf2117*MStau2(2)

	Sf2114 = Sf2112*MSneu2(1) + 2*Sf2023*MSneu2(1)**2 + 
     &   Sf2113*MStau2(2)

	Sf2039 = Sf2037*MSneu2(1) + Sf299*MSneu2(1)**2 + 
     &   Sf2038*MStau2(1)

	Sf2036 = Sf2034*MSneu2(1) + Sf2033*MSneu2(1)**2 + 
     &   3*Sf2035*MStau2(1)

	Sf2030 = Sf2028*MSneu2(1) + 2*Sf2023*MSneu2(1)**2 + 
     &   Sf2029*MStau2(1)

	Sf1603 = Sf1601*MSneu2(1) - Sf745*MSneu2(1)**2 + 
     &   Sf1602*MStau2(2)

	Sf1589 = Sf1587*MSneu2(1) - 2*C2B*MSneu2(1)**2 + 
     &   Sf1588*MStau2(2)

	Sf1586 = Sf1584*MSneu2(1) + Sf243*MSneu2(1)**2 + 
     &   Sf1585*MStau2(2)

	Sf1581 = Sf1579*MSneu2(1) + Sf299*MSneu2(1)**2 + 
     &   Sf1580*MStau2(2)

	Sf1562 = Sf1560*MSneu2(1) + Sf690*MSneu2(1)**2 + 
     &   Sf1561*MStau2(2)

	Sf1557 = Sf1555*MSneu2(1) - Sf684*MSneu2(1)**2 - 
     &   3*Sf1556*MStau2(2)

	Sf1551 = Sf1549*MSneu2(1) - Sf415*MSneu2(1)**2 - 
     &   Sf1550*MStau2(2)

	Sf1491 = Sf1489*MSneu2(1) + Sf300*MSneu2(1)**2 + 
     &   Sf1490*MStau2(2)

	Sf1437 = Sf1436*MSneu2(1) + Sf526*MSneu2(1)**2 + 
     &   Sf1435*MStau2(2)

	Sf1264 = Sf1263*MSneu2(1) + 13*ML2*MStau2(2) - 
     &   Sf300*MStau2(2)**2

	Sf1258 = Sf1256*MSneu2(1) - Sf288*MSneu2(1)**2 - 
     &   Sf1257*MStau2(2)

	Sf1249 = Sf1247*MSneu2(1) + Sf272*MSneu2(1)**2 - 
     &   Sf1248*MStau2(2)

	Sf6151 = Sf3759*MSbot2(2)**3 + 6*MBot2*MTop2*MSq2Diff(2,2) + 
     &   Sf6150*MStop2(2)

	Sf3903 = Sf3759*MSbot2(2)**3 + 6*MBot2*MTop2*MSq2Diff(2,1) + 
     &   Sf3902*MStop2(1)

	Sf6170 = 3*Sf3778*MSbot2(2)**2 + Sf3923*MStop2(2) - 
     &   Sf3780*MStop2(2)**2

	Sf3924 = 3*Sf3778*MSbot2(2)**2 + Sf3923*MStop2(1) - 
     &   Sf3780*MStop2(1)**2

	Sf6171 = 12*MTop2*MSbot2(2)**2 + Sf3777*MSbot2(2)**3 - 
     &   Sf6170*MStop2(2)

	Sf3925 = 12*MTop2*MSbot2(2)**2 + Sf3777*MSbot2(2)**3 - 
     &   Sf3924*MStop2(1)

	Sf6801 = 343*MSbot2(2)**2 - 7*Sf4452*MStop2(2) + 
     &   Sf4328*MStop2(2)**2

	Sf4456 = 343*MSbot2(2)**2 - 7*Sf4452*MStop2(1) + 
     &   Sf4328*MStop2(1)**2

	Sf8712 = Sf4906*Sf8711 + Sf8402*MSbot2(2) + 
     &   Sf8710*MSq2Diff(2,2)**2

	Sf8675 = Sf2334*Sf8674 + Sf8378*MSbot2(2) - 
     &   Sf8673*MSq2Diff(2,1)**2

	Sf6800 = Sf4453*MSbot2(2) - 9*Sf4454*MStop2(2) + 
     &   Sf4332*MStop2(2)**2

	Sf4455 = Sf4453*MSbot2(2) - 9*Sf4454*MStop2(1) + 
     &   Sf4332*MStop2(1)**2

	Sf5936 = Sf3495*MSbot2(2) + 3*Sf3496*MStop2(2) + 
     &   Sf3119*MStop2(2)**2

	Sf3497 = Sf3495*MSbot2(2) + 3*Sf3496*MStop2(1) + 
     &   Sf3119*MStop2(1)**2

	Sf5935 = Sf3492*MSbot2(2) + 3*Sf3493*MStop2(2) + 
     &   Sf3114*MStop2(2)**2

	Sf3494 = Sf3492*MSbot2(2) + 3*Sf3493*MStop2(1) + 
     &   Sf3114*MStop2(1)**2

	Sf5077 = 3*Sf2634*MSbot2(2) + Sf2635*MStop2(2) + 
     &   Sf873*MStop2(2)**2

	Sf2636 = 3*Sf2634*MSbot2(2) + Sf2635*MStop2(1) + 
     &   Sf873*MStop2(1)**2

	Sf6908 = 30*MTop2*MSbot2(2) - Sf4614*MSbot2(2)**2 - 
     &   Sf6907*MStop2(2)

	Sf4713 = 30*MTop2*MSbot2(2) - Sf4614*MSbot2(2)**2 - 
     &   Sf4712*MStop2(1)

	Sf5426 = 10*MTop2*MSbot2(2) + Sf2875*MSbot2(2)**2 - 
     &   Sf5425*MStop2(2)

	Sf3027 = 10*MTop2*MSbot2(2) + Sf2875*MSbot2(2)**2 - 
     &   Sf3026*MStop2(1)

	Sf6169 = 6*MTop2*MSbot2(2) + Sf3795*MSbot2(2)**2 - 
     &   Sf6168*MStop2(2)

	Sf3922 = 6*MTop2*MSbot2(2) + Sf3795*MSbot2(2)**2 - 
     &   Sf3921*MStop2(1)

	Sf6916 = 30*MBot2*MSbot2(2) + Sf4614*MSbot2(2)**2 - 
     &   Sf6915*MStop2(2)

	Sf4724 = 30*MBot2*MSbot2(2) + Sf4614*MSbot2(2)**2 - 
     &   Sf4723*MStop2(1)

	Sf5067 = 21*MBot2*MSbot2(2) - Sf2387*MSbot2(2)**2 - 
     &   Sf5066*MStop2(2)

	Sf2617 = 21*MBot2*MSbot2(2) - Sf2387*MSbot2(2)**2 - 
     &   Sf2616*MStop2(1)

	Sf5433 = 10*MBot2*MSbot2(2) - Sf2875*MSbot2(2)**2 - 
     &   Sf5432*MStop2(2)

	Sf3038 = 10*MBot2*MSbot2(2) - Sf2875*MSbot2(2)**2 - 
     &   Sf3037*MStop2(1)

	Sf6176 = 6*MBot2*MSbot2(2) - 3*Sf3773*MSbot2(2)**2 - 
     &   Sf6175*MStop2(2)

	Sf3931 = 6*MBot2*MSbot2(2) - 3*Sf3773*MSbot2(2)**2 - 
     &   Sf3930*MStop2(1)

	Sf5768 = Sf3759*MSbot2(1)**3 + 6*MBot2*MTop2*MSq2Diff(1,2) + 
     &   Sf5767*MStop2(2)

	Sf3764 = Sf3759*MSbot2(1)**3 + 6*MBot2*MTop2*MSq2Diff(1,1) + 
     &   Sf3763*MStop2(1)

	Sf5792 = 3*Sf3778*MSbot2(1)**2 + Sf3799*MStop2(2) - 
     &   Sf3780*MStop2(2)**2

	Sf3800 = 3*Sf3778*MSbot2(1)**2 + Sf3799*MStop2(1) - 
     &   Sf3780*MStop2(1)**2

	Sf5793 = 12*MTop2*MSbot2(1)**2 + Sf3777*MSbot2(1)**3 - 
     &   Sf5792*MStop2(2)

	Sf3801 = 12*MTop2*MSbot2(1)**2 + Sf3777*MSbot2(1)**3 - 
     &   Sf3800*MStop2(1)

	Sf9085 = MHiggs2(2)*MSbot2(1)**2 - MSbot2(1)**3 - 
     &   Sf9084*MSbot2(2)**2

	Sf9265 = MHiggs2(1)*MSbot2(1)**2 - MSbot2(1)**3 - 
     &   Sf9264*MSbot2(2)**2

	Sf6546 = 343*MSbot2(1)**2 - 7*Sf4327*MStop2(2) + 
     &   Sf4328*MStop2(2)**2

	Sf4334 = 343*MSbot2(1)**2 - 7*Sf4327*MStop2(1) + 
     &   Sf4328*MStop2(1)**2

	Sf6545 = Sf4330*MSbot2(1) - 9*Sf4331*MStop2(2) + 
     &   Sf4332*MStop2(2)**2

	Sf4333 = Sf4330*MSbot2(1) - 9*Sf4331*MStop2(1) + 
     &   Sf4332*MStop2(1)**2

	Sf5503 = Sf3117*MSbot2(1) + 3*Sf3118*MStop2(2) + 
     &   Sf3119*MStop2(2)**2

	Sf3120 = Sf3117*MSbot2(1) + 3*Sf3118*MStop2(1) + 
     &   Sf3119*MStop2(1)**2

	Sf5502 = Sf3112*MSbot2(1) + 3*Sf3113*MStop2(2) + 
     &   Sf3114*MStop2(2)**2

	Sf3115 = Sf3112*MSbot2(1) + 3*Sf3113*MStop2(1) + 
     &   Sf3114*MStop2(1)**2

	Sf4951 = 3*Sf2416*MSbot2(1) + Sf2418*MStop2(2) + 
     &   Sf873*MStop2(2)**2

	Sf2419 = 3*Sf2416*MSbot2(1) + Sf2418*MStop2(1) + 
     &   Sf873*MStop2(1)**2

	Sf8899 = 4*MZ2*MSbot2(1) + MSbot2(1)**2 + Sf8898*MSbot2(2)

	Sf8949 = 4*MZ2*MSbot2(1) - MSbot2(1)**2 + Sf8948*MSbot2(2)

	Sf6657 = 30*MTop2*MSbot2(1) - Sf4614*MSbot2(1)**2 - 
     &   Sf6656*MStop2(2)

	Sf4618 = 30*MTop2*MSbot2(1) - Sf4614*MSbot2(1)**2 - 
     &   Sf4617*MStop2(1)

	Sf5320 = 10*MTop2*MSbot2(1) + Sf2875*MSbot2(1)**2 - 
     &   Sf5319*MStop2(2)

	Sf2900 = 10*MTop2*MSbot2(1) + Sf2875*MSbot2(1)**2 - 
     &   Sf2899*MStop2(1)

	Sf5791 = 6*MTop2*MSbot2(1) + Sf3795*MSbot2(1)**2 - 
     &   Sf5790*MStop2(2)

	Sf3798 = 6*MTop2*MSbot2(1) + Sf3795*MSbot2(1)**2 - 
     &   Sf3797*MStop2(1)

	Sf6665 = 30*MBot2*MSbot2(1) + Sf4614*MSbot2(1)**2 - 
     &   Sf6664*MStop2(2)

	Sf4629 = 30*MBot2*MSbot2(1) + Sf4614*MSbot2(1)**2 - 
     &   Sf4628*MStop2(1)

	Sf4941 = 21*MBot2*MSbot2(1) - Sf2387*MSbot2(1)**2 - 
     &   Sf4940*MStop2(2)

	Sf2391 = 21*MBot2*MSbot2(1) - Sf2387*MSbot2(1)**2 - 
     &   Sf2390*MStop2(1)

	Sf5327 = 10*MBot2*MSbot2(1) - Sf2875*MSbot2(1)**2 - 
     &   Sf5326*MStop2(2)

	Sf2911 = 10*MBot2*MSbot2(1) - Sf2875*MSbot2(1)**2 - 
     &   Sf2910*MStop2(1)

	Sf5798 = 6*MBot2*MSbot2(1) - 3*Sf3773*MSbot2(1)**2 - 
     &   Sf5797*MStop2(2)

	Sf3807 = 6*MBot2*MSbot2(1) - 3*Sf3773*MSbot2(1)**2 - 
     &   Sf3806*MStop2(1)

	Sf8493 = MZ2*Sf8492*SW2 + 
     &   MW2*USbot11*(4*MZ2*Sf8489*USbot12**2 + 
     &      Sf15**2*Sf4773*Sf7756*USbot21**2)

	Sf9951 = Sf9944 - Sf9948*MStau2(2) + 4*Sf9950*MStau2(2)**2

	Sf9858 = Sf8923*Sf9857 - Sf9855*MStau2(1) + 
     &   Sf9363*Sf9364*MHiggs2(2)**2*MStau2(2)

	Sf9848 = Sf8923*Sf9847 - 
     &   Sf9332*Sf9333*MHiggs2(2)**2*MStau2(1) + Sf9844*MStau2(2)

	Sf9387 = Sf9386 + Sf9383*MStau2(1) - Sf8911*MStau2(1)**2

	Sf9380 = Sf9378 + Sf9379*MStau2(2) - Sf8911*MStau2(2)**2

	Sf9277 = Sf9125*Sf9276 - Sf9274*MSbot2(1) + 
     &   Sf9161*Sf9162*MHiggs2(1)**2*MSbot2(2)

	Sf9266 = Sf9125*Sf9265 - 
     &   Sf9123*Sf9124*MHiggs2(1)**2*MSbot2(1) + Sf9262*MSbot2(2)

	Sf9273 = Sf9161*Sf9162*Sf9245 - 2*Sf9244*MSbot2(1) + 
     &   MSbot2(1)**2

	Sf9260 = Sf9123*Sf9124*Sf9223 - 2*Sf9224*MSbot2(2) + 
     &   MSbot2(2)**2

	Sf10433 = -2*MZ2*MStop2(3)**3 + 
     &   MHiggs2(1)**3*(-2*MStop2(1) + MStop2(3)) + 
     &   2*MHiggs2(1)**2*(MZ2*MStop2(1) - 
     &      (3*MStop2(1) + 2*MStop2(2))*MStop2(3)) + 
     &   MHiggs2(1)*MStop2(3)*
     &    (2*MZ2*(2*MStop2(1) + MStop2(2)) + 3*MStop2(3)**2)

	Sf10423 = Sf10422*Sf9125 - 
     &   Sf10335*Sf10336*MHiggs2(1)**2*MStop2(1) + 
     &   Sf10419*MStop2(2)

	Sf10209 = Sf10208*Sf9125 - Sf10206*MStau2(1) + 
     &   Sf10099*Sf10100*MHiggs2(1)**2*MStau2(2)

	Sf10198 = Sf10197*Sf9125 - 
     &   Sf10069*Sf10070*MHiggs2(1)**2*MStau2(1) + 
     &   Sf10194*MStau2(2)

	Sf9096 = Sf8923*Sf9095 - Sf9093*MSbot2(1) + 
     &   Sf8971*Sf8972*MHiggs2(2)**2*MSbot2(2)

	Sf9086 = Sf8923*Sf9085 - 
     &   Sf8921*Sf8922*MHiggs2(2)**2*MSbot2(1) + Sf9082*MSbot2(2)

	Sf9010 = Sf9009 + Sf9006*MSbot2(1) - Sf8911*MSbot2(1)**2

	Sf9003 = Sf9001 + Sf9002*MSbot2(2) - Sf8911*MSbot2(2)**2

	Sf10584 = Sf10583*Sf8923 - Sf10581*MStop2(1) + 
     &   Sf10479*Sf10480*MHiggs2(2)**2*MStop2(2)

	Sf10573 = Sf10572*Sf8923 - 
     &   Sf10445*Sf10446*MHiggs2(2)**2*MStop2(1) + 
     &   Sf10569*MStop2(2)

	Sf8857 = MZ2*(MZ2*Sf8856 + 
     &      672*Sf8820*Ss443*SW2*UStau11*UStau12) - 
     &   168*Sf8819*Ss443**2*UStau22**2

	Sf8718 = Sf7990*Sf8659 + Sf8532*Sf8659 + 4*Sf8717*Ss78

	Sf8255 = Sf8251 + Sf4820*Sf8254 + Sf8166*MStop2(2)

	Sf8226 = Sf8222 + Sf4906*Sf8225 + Sf8186*MStop2(2)

	Sf8197 = Sf8193 - Sf4906*Sf8196 - Sf8186*MSbot2(2)

	Sf7988 = Sf7984 - Sf2334*Sf7987 - Sf7977*MSbot2(2)

	Sf7228 = Sf1749*Sf7215 + Sf1743*Sf7218 + 
     &   2*SB*Sf7227*Ss25**2*UStop21**2

	Sf7209 = Sf1749*Sf7193 + Sf1743*Sf7197 + 
     &   2*SB*Sf7208*Ss25**2*UStop22**2

	Sf6802 = Sf6800*MHiggs2(4) - 3*Sf4451*MSbot2(2)**2 - 
     &   Sf6801*MStop2(2)

	Sf6547 = Sf6545*MHiggs2(4) - 3*Sf4326*MSbot2(1)**2 - 
     &   Sf6546*MStop2(2)

	Sf6211 = 2*Sf6209*MHiggs2(4)**2 + 
     &   2*Sf6208*MSbot2(2)**2 - Sf6210*MStop2(2)**2

	Sf5834 = Sf5832*MHiggs2(4)**2 + 2*Sf5830*MSbot2(1)**2 - 
     &   Sf5833*MStop2(2)**2

	Sf5422 = Sf5419 + 72*Sf2552*MSq2Diff(2,2)**5 + 
     &   76*MSq2Diff(2,2)**6

	Sf5316 = Sf5313 + 72*Sf2552*MSq2Diff(1,2)**5 + 
     &   76*MSq2Diff(1,2)**6

	Sf1255 = 9*ML2*Sf47 + Sf284*MSneu2(1)**2 - Sf1254*MStau2(2)

	Sf1502 = 3*ML2*Sf47 - Sf601*MSneu2(1)**2 + Sf1501*MStau2(2)

	Sf1488 = ML2*Sf47 + 3*Sf174*MSneu2(1)**2 - Sf1487*MStau2(2)

	Sf7031 = MW2**2*Sf4520 - Sf4521*MHiggs2(4) - Sf4523*MStop2(2)

	Sf4524 = MW2**2*Sf4520 - Sf4521*MHiggs2(4) - Sf4523*MStop2(1)

	Sf4457 = Sf4455*MHiggs2(4) - 3*Sf4451*MSbot2(2)**2 - 
     &   Sf4456*MStop2(1)

	Sf6784 = Sf4426 - Sf4427*MStop2(2) - 3*Sf4298*MStop2(2)**2

	Sf4428 = Sf4426 - Sf4427*MStop2(1) - 3*Sf4298*MStop2(1)**2

	Sf7009 = MW2**2*Sf4416 - Sf4417*MHiggs2(4) - Sf4419*MStop2(2)

	Sf4420 = MW2**2*Sf4416 - Sf4417*MHiggs2(4) - Sf4419*MStop2(1)

	Sf4335 = Sf4333*MHiggs2(4) - 3*Sf4326*MSbot2(1)**2 - 
     &   Sf4334*MStop2(1)

	Sf709 = ML2*Sf430 - Sf301*MSneu2(1)**2 - Sf708*MStau2(1)

	Sf6530 = Sf4294 - Sf4297*MStop2(2) - 3*Sf4298*MStop2(2)**2

	Sf4299 = Sf4294 - Sf4297*MStop2(1) - 3*Sf4298*MStop2(1)**2

	Sf4150 = 2*Sf4148*MHiggs2(4)**2 + 
     &   2*Sf4147*MSbot2(2)**2 - Sf4149*MStop2(1)**2

	Sf4039 = 2*Sf4037*MHiggs2(4)**2 + 
     &   2*Sf4036*MSbot2(1)**2 - Sf4038*MStop2(1)**2

	Sf287 = 9*ML2*Sf34 + Sf284*MSneu2(1)**2 - Sf286*MStau2(1)

	Sf605 = 3*ML2*Sf34 - Sf601*MSneu2(1)**2 + Sf604*MStau2(1)

	Sf587 = ML2*Sf34 + 3*Sf174*MSneu2(1)**2 - Sf586*MStau2(1)

	Sf3022 = Sf3019 + 72*Sf2552*MSq2Diff(2,1)**5 + 
     &   76*MSq2Diff(2,1)**6

	Sf2895 = Sf2892 + 72*Sf2552*MSq2Diff(1,1)**5 + 
     &   76*MSq2Diff(1,1)**6

	Sf6248 = Sf2722**3 - Sf6138*MStop2(2) + Sf6139*MStop2(2)**2

	Sf6140 = Sf2722**3 - Sf6138*MStop2(1) + Sf6139*MStop2(1)**2

	Sf5884 = Sf2543**3 - Sf5755*MStop2(2) + Sf5756*MStop2(2)**2

	Sf5757 = Sf2543**3 - Sf5755*MStop2(1) + Sf5756*MStop2(1)**2

	Sf225 = -((-MW2 + MHiggs2(4))*
     &     (MHiggs2(4)*(-MSl2Diff(1,1)**4 + 
     &          3*MW2**3*(MSneu2(1) + MStau2(1)) + 
     &          5*MW2*MSl2Diff(1,1)**2*(MSneu2(1) + MStau2(1)) - 
     &          7*MW2**2*(MSneu2(1)**2 + MStau2(1)**2)) + 
     &       2*MW2*(MSl2Diff(1,1)**4 - 
     &          2*MW2**3*(MSneu2(1) + MStau2(1)) - 
     &          4*MW2*MSl2Diff(1,1)**2*(MSneu2(1) + MStau2(1)) + 
     &          MW2**2*(5*MStau2(1)**2 + 
     &             MSneu2(1)*(5*MSneu2(1) + MStau2(1))))))

	Sf1217 = -((-MW2 + MHiggs2(4))*
     &     (MHiggs2(4)*(-MSl2Diff(2,1)**4 + 
     &          3*MW2**3*(MSneu2(1) + MStau2(2)) + 
     &          5*MW2*MSl2Diff(2,1)**2*(MSneu2(1) + MStau2(2)) - 
     &          7*MW2**2*(MSneu2(1)**2 + MStau2(2)**2)) + 
     &       2*MW2*(MSl2Diff(2,1)**4 - 
     &          2*MW2**3*(MSneu2(1) + MStau2(2)) - 
     &          4*MW2*MSl2Diff(2,1)**2*(MSneu2(1) + MStau2(2)) + 
     &          MW2**2*(5*MStau2(2)**2 + 
     &             MSneu2(1)*(5*MSneu2(1) + MStau2(2))))))

	Sf1883 = Sf1882*MHiggs2(4)**2 + 
     &   2*Sf1881*MHiggs2(4)*MStau2(2) + 
     &   Sf1137*MStau2(2)**2*UStau21

	Sf1878 = Sf1877*MHiggs2(4)**2 + 
     &   Sf1876*MHiggs2(4)*MStau2(2) - 
     &   Sf1125*MStau2(2)**2*UStau21

	Sf1203 = 2*MW2**2*Sf187 + Sf176*Sf203*MHiggs2(4) + 
     &   10*MStau2(2)**3

	Sf204 = 2*MW2**2*Sf187 + Sf176*Sf203*MHiggs2(4) + 
     &   10*MStau2(1)**3

	Sf1575 = ML2*Sf1353 - Sf301*MSneu2(1)**2 - Sf1574*MStau2(2)

	Sf1142 = Sf1141*MHiggs2(4)**2 + 
     &   2*Sf1139*MHiggs2(4)*MStau2(1) + 
     &   Sf1137*MStau2(1)**2*UStau11

	Sf1132 = Sf1131*MHiggs2(4)**2 + 
     &   Sf1128*MHiggs2(4)*MStau2(1) - 
     &   Sf1125*MStau2(1)**2*UStau11

	Sf7272 = Sf108*Sf110*Sf112**2 + Sf21*Sf23*Sf25**2 + 
     &   3*Sf56*Sf58*Sf60**2

	Sf10580 = Sf10479*Sf10480*Sf10540 - 
     &   2*Sf10539*MStop2(1) + MStop2(1)**2

	Sf10568 = Sf10445*Sf10446*Sf10519 - 
     &   2*Sf10520*MStop2(2) + MStop2(2)**2

	Sf10205 = Sf10099*Sf10100*Sf10179 - 
     &   2*Sf10178*MStau2(1) + MStau2(1)**2

	Sf10193 = Sf10069*Sf10070*Sf10158 - 
     &   2*Sf10159*MStau2(2) + MStau2(2)**2

	Sf1161 = 5*CA2*SB2 - CA*Sf1158 - SA*Sf1160

	Sf10 = MZ2**2 - MZ2*Sf9 + 2*MZ*Sqrt(Sf1)*MSneu2(1)

	Sf85 = MZ2**2 - MZ2*Sf84 + 2*MZ*Sqrt(Sf82)*MStop2(1)

	Sf64 = MZ2**2 - MZ2*Sf63 + 2*MZ*Sqrt(Sf61)*MStau2(2)

	Sf42 = MZ2**2 - MZ2*Sf41 + 2*MZ*Sqrt(Sf39)*MStau2(1)

	Sf29 = MZ2**2 - MZ2*Sf28 + 2*MZ*Sqrt(Sf26)*MSbot2(2)

	Sf14 = MZ2**2 - MZ2*Sf13 + 2*MZ*Sqrt(Sf11)*MSbot2(1)

	Sf116 = MZ2**2 - MZ2*Sf115 + 2*MZ*Sqrt(Sf113)*MStop2(2)

	Sf10028 = MZ2*(-MZ2 + 3*MHiggs2(2)) + 2*Sf10027*MStau2(2)

	Sf9930 = MZ2*(-MZ2 + 3*MHiggs2(2)) + 2*Sf9929*MStau2(1)

	Sf54 = Sf52 - Sqrt(Sf53)

	Sf19 = Sf17 - Sqrt(Sf18)

	Sf106 = Sf104 - Sqrt(Sf105)

	Sf1102 = 29*ML2*MW2 + MW2**2*Sf1099 - Sf1101*MHiggs2(4)

	Sf1113 = 13*ML2*MW2 - MW2**2*Sf1111 - Sf1112*MHiggs2(4)

	Sf1944 = 8*ML2*MW2 - 4*MW2**2*Sf720 - 2*Sf1943*MHiggs2(4)

	Sf1952 = 5*ML2*MW2 - MW2**2*Sf1950 + Sf1951*MHiggs2(4)

	Sf1947 = 4*ML2*MW2 - 4*C2B*MW2**2 - Sf1946*MHiggs2(4)

	Sf1117 = 3*ML2*MW2 + MW2**2*Sf1116 + Sf963*MHiggs2(4)

	Sf202 = 2*ML2*MW2 + MW2**2*Sf184 - Sf176*MHiggs2(4)

	Sf177 = ML2*MW2 - MW2**2*Sf174 - Sf176*MHiggs2(4)

	Sf1097 = ML2*MW2 + MW2**2*Sf1095 - 3*Sf1096*MHiggs2(4)

	Sf1032 = -Sf1014 + 8*MHiggs2(4)

	Sf192 = -Sf176 + MHiggs2(4)

	Sf5693 = 12*MBot2 - Sf5691

	Sf3371 = 6*MBot2 - Sf3369

	Sf3867 = 6*MBot2 - Sf3865

	Sf5683 = 3*MBot2 - Sf5681

	Sf5668 = 2*MBot2 - Sf5666

	Sf3854 = 2*MBot2 - Sf3852

	Sf3351 = MBot2 - Sf2238

	Sf8817 = 2 - Re((Sf116*Log(-(1/2.D0*Sf114/MStop2(2))))/
     &     (MZ2*Sf114)) - Re(Log(MStop2(2)))

	Sf8800 = 2 - Re((Sf85*Log(-(1/2.D0*Sf83/MStop2(1))))/
     &     (MZ2*Sf83)) - Re(Log(MStop2(1)))

	Sf8787 = 2 - Re((Sf64*Log(-(1/2.D0*Sf62/MStau2(2))))/
     &     (MZ2*Sf62)) - Re(Log(MStau2(2)))

	Sf8776 = 2 - Re((Sf42*Log(-(1/2.D0*Sf40/MStau2(1))))/
     &     (MZ2*Sf40)) - Re(Log(MStau2(1)))

	Sf8763 = 2 - Re((Sf10*Log(-(1/2.D0*Sf8/MSneu2(1))))/
     &     (MZ2*Sf8)) - Re(Log(MSneu2(1)))

	Sf8769 = 2 - Re((Sf29*Log(-(1/2.D0*Sf27/MSbot2(2))))/
     &     (MZ2*Sf27)) - Re(Log(MSbot2(2)))

	Sf8764 = 2 - Re((Sf14*Log(-(1/2.D0*Sf12/MSbot2(1))))/
     &     (MZ2*Sf12)) - Re(Log(MSbot2(1)))

	Sf491 = 3*ML2*Sf487 - 10*MSl2Diff(1,1)**8 + 
     &   C2B*MW2**3*(MW2**2*(-6*MW2**2*Sf470 - 9*Sf475) + 
     &      18*Sf472*MSl2Diff(1,1)**3) + 
     &   MW2*(Sf490 + MW2*(MW2**2*
     &          (MW2**2*(2*MW2**2 - Sf473) - 2*Sf489) - 
     &         Sf488*MSl2Diff(1,1)**4))

	Sf1408 = 3*ML2*Sf1404 - 10*MSl2Diff(2,1)**8 + 
     &   C2B*MW2**3*(MW2**2*(-6*MW2**2*Sf1387 - 9*Sf1392) + 
     &      18*Sf1389*MSl2Diff(2,1)**3) + 
     &   MW2*(Sf1407 + MW2*
     &       (MW2**2*(MW2**2*(2*MW2**2 - Sf1390) - 2*Sf1406) - 
     &         Sf1405*MSl2Diff(2,1)**4))

	tmp39 = -3*MBot2 + 3*MTop2 - 53*MW2 + 3*C2B*MW2

	tmp40 = -3*MBot2 + 3*MTop2 + 19*MW2 + 3*C2B*MW2

	tmp41 = -6*MBot2 + 6*MTop2 + 11*MW2 + 6*C2B*MW2

	Sf6185 = (12*MSq2Diff(2,2)*(-MW2 + MSq2Diff(2,2))**5 - 
     &      MW2*(MW2 - MSq2Diff(2,2))*
     &       (-112*MSq2Diff(2,2)**3 + 
     &         MW2*(MW2*tmp41 + 2*tmp39*MSq2Diff(2,2) + 
     &            219*MSq2Diff(2,2)**2))*MStop2(2) + 
     &      4*MW2**2*(MW2*tmp40 + 
     &         3*(-MBot2 + MTop2 - 22*MW2 + C2B*MW2)*MSq2Diff(2,2) + 
     &         75*MSq2Diff(2,2)**2)*MStop2(2)**2 - 
     &      128*MW2**3*MStop2(2)**3)*USbot21 - 
     &   6*MBot*MW2**2*(reMUE*TB - Re(Abot))*
     &    (MW2 + MSq2Diff(2,2))*(MW2 - MSbot2(2) - MStop2(2))*
     &    MStop2(2)*USbot22

	Sf5807 = (12*MSq2Diff(1,2)*(-MW2 + MSq2Diff(1,2))**5 - 
     &      MW2*(MW2 - MSq2Diff(1,2))*
     &       (-112*MSq2Diff(1,2)**3 + 
     &         MW2*(MW2*tmp41 + 2*tmp39*MSq2Diff(1,2) + 
     &            219*MSq2Diff(1,2)**2))*MStop2(2) + 
     &      4*MW2**2*(MW2*tmp40 + 
     &         3*(-MBot2 + MTop2 - 22*MW2 + C2B*MW2)*MSq2Diff(1,2) + 
     &         75*MSq2Diff(1,2)**2)*MStop2(2)**2 - 
     &      128*MW2**3*MStop2(2)**3)*USbot11 - 
     &   6*MBot*MW2**2*(reMUE*TB - Re(Abot))*
     &    (MW2 + MSq2Diff(1,2))*(MW2 - MSbot2(1) - MStop2(2))*
     &    MStop2(2)*USbot12

	tmp36 = -6*MBot2 + 6*MTop2 - 89*MStop2(1)

	Sf6103 = (12*MSq2Diff(2,1)**6 - 
     &      4*MW2*MSq2Diff(2,1)**4*(15*MSbot2(2) + 13*MStop2(1)) + 
     &      MW2**5*(-12*MSbot2(2) + (1 - 6*C2B)*MStop2(1)) + 
     &      MW2**3*(-120*MSbot2(2)**3 + 
     &         (35 + 6*C2B)*MSbot2(2)**2*MStop2(1) + 
     &         2*(-6*MBot2 + 6*MTop2 + 13*MSbot2(2))*MStop2(1)**2 + 
     &         (-69 - 6*C2B)*MStop2(1)**3) + 
     &      MW2**4*(MSbot2(2)*(60*MSbot2(2) - 3*MStop2(1)) + 
     &         6*MBot2*MStop2(1) + 
     &         MStop2(1)*(-6*MTop2 + (19 + 12*C2B)*MStop2(1))) + 
     &      MW2**2*MSq2Diff(2,1)*
     &       (tmp36*MStop2(1)**2 + 
     &         MSbot2(2)*(MSbot2(2)*
     &             (120*MSbot2(2) - 29*MStop2(1)) - 
     &            6*MBot2*MStop2(1) + (6*MTop2 - 2*MStop2(1))*MStop2(1)
     &            )))*USbot21 - 
     &   6*MBot*MW2**2*(reMUE*TB - Re(Abot))*
     &    (-MSbot2(2)**2 + (MW2 - MStop2(1))**2)*MStop2(1)*
     &    USbot22

	Sf5714 = (12*MSq2Diff(1,1)**6 - 
     &      4*MW2*MSq2Diff(1,1)**4*(15*MSbot2(1) + 13*MStop2(1)) + 
     &      MW2**5*(-12*MSbot2(1) + (1 - 6*C2B)*MStop2(1)) + 
     &      MW2**3*(-120*MSbot2(1)**3 + 
     &         (35 + 6*C2B)*MSbot2(1)**2*MStop2(1) + 
     &         2*(-6*MBot2 + 6*MTop2 + 13*MSbot2(1))*MStop2(1)**2 + 
     &         (-69 - 6*C2B)*MStop2(1)**3) + 
     &      MW2**4*(MSbot2(1)*(60*MSbot2(1) - 3*MStop2(1)) + 
     &         6*MBot2*MStop2(1) + 
     &         MStop2(1)*(-6*MTop2 + (19 + 12*C2B)*MStop2(1))) + 
     &      MW2**2*MSq2Diff(1,1)*
     &       (tmp36*MStop2(1)**2 + 
     &         MSbot2(1)*(MSbot2(1)*
     &             (120*MSbot2(1) - 29*MStop2(1)) - 
     &            6*MBot2*MStop2(1) + (6*MTop2 - 2*MStop2(1))*MStop2(1)
     &            )))*USbot11 - 
     &   6*MBot*MW2**2*(reMUE*TB - Re(Abot))*
     &    (-MSbot2(1)**2 + (MW2 - MStop2(1))**2)*MStop2(1)*
     &    USbot12

	Sf658 = Sf657*MSneu2(1) + 
     &   MSneu2(1)**2*(-Sf654 - Sf326*MSneu2(1)) - 
     &   Sf651*MStau2(1)**3 + 
     &   MHiggs2(4)*(-2*Sf650*MStau2(1)**2 + 
     &      MHiggs2(4)*(6*MHiggs2(4) + Sf649*MStau2(1)))

	Sf1538 = Sf1537*MSneu2(1) + 
     &   MSneu2(1)**2*(-Sf1536 - Sf326*MSneu2(1)) - 
     &   Sf651*MStau2(2)**3 + 
     &   MHiggs2(4)*(-2*Sf650*MStau2(2)**2 + 
     &      MHiggs2(4)*(6*MHiggs2(4) + Sf649*MStau2(2)))

	Sf5142 = TB*(12*Sf2731*MSbot2(2)*MStop2(2) - 
     &      (6*Sf2552*Sf5141 + 5*MSq2Diff(2,2)**3)*USbot21)*
     &    UStop21 + 6*MBot*Sf2225*TB*
     &    (MSbot2(2)**2 + MStop2(2)**2)*USbot22*UStop21 - 
     &   6*MTop*Sf2245*MSq2Diff(2,2)**2*USbot21*UStop22

	Sf5023 = TB*(12*Sf2553*MSbot2(1)*MStop2(2) - 
     &      (6*Sf2552*Sf5022 + 5*MSq2Diff(1,2)**3)*USbot11)*
     &    UStop21 + 6*MBot*Sf2225*TB*
     &    (MSbot2(1)**2 + MStop2(2)**2)*USbot12*UStop21 - 
     &   6*MTop*Sf2245*MSq2Diff(1,2)**2*USbot11*UStop22

	Sf2734 = TB*(12*Sf2731*MSbot2(2)*MStop2(1) - 
     &      (6*Sf2552*Sf2733 + 5*MSq2Diff(2,1)**3)*USbot21)*
     &    UStop11 + 6*MBot*Sf2225*TB*
     &    (MSbot2(2)**2 + MStop2(1)**2)*USbot22*UStop11 - 
     &   6*MTop*Sf2245*MSq2Diff(2,1)**2*USbot21*UStop12

	Sf2556 = TB*(12*Sf2553*MSbot2(1)*MStop2(1) - 
     &      (6*Sf2552*Sf2555 + 5*MSq2Diff(1,1)**3)*USbot11)*
     &    UStop11 + 6*MBot*Sf2225*TB*
     &    (MSbot2(1)**2 + MStop2(1)**2)*USbot12*UStop11 - 
     &   6*MTop*Sf2245*MSq2Diff(1,1)**2*USbot11*UStop12

	Sf8480 = Sf8464*Ss906*MSbot2(2)*USbot21**2 + 
     &   Sf8466*Ss885*MSbot2(1)*USbot22**2 + 
     &   Sf8462*Ss988*MStop2(2)*UStop21**2 + 
     &   Sf8438*Ss969*MStop2(1)*UStop22**2

	Sf715 = Sf714*MSneu2(1) - Sf712*MSneu2(1)**2 + 
     &   11*ML2*MStau2(1) - Sf713*MStau2(1)**2

	Sf1578 = Sf1577*MSneu2(1) - Sf712*MSneu2(1)**2 + 
     &   11*ML2*MStau2(2) - Sf713*MStau2(2)**2

	Sf944 = CA2*ML*Sf940 + SA*Sf943 + 
     &   ML*SB*Re(Atau)*UStau11 - 
     &   2*CA*CB*MZ2*Sf127*SW2*UStau12

	Sf8827 = 72*Sf8819*Ss443**2*UStau22**2 + 
     &   MZ2*(MZ2*(MZ2*Sf8822 + 6*Sf8826) + 
     &      288*Sf8820*Ss443*SW2*UStau21*UStau22)

	Sf6981 = 3*CA2*CB*MTop*Sf6978 + SA*Sf6980 - 
     &   CA*SB*Sf127*Sf86*UStop21 + 3*MTop*reMUE*SB*UStop22

	Sf688 = Sf687*MHiggs2(4) - 7*Sf681*MSneu2(1)**2 + 
     &   Sf680*MSneu2(1)**3 + Sf683*MStau2(1)**2

	Sf4264 = 3*CA2*CB*MTop*Sf4261 + SA*Sf4263 - 
     &   CA*SB*Sf127*Sf86*UStop11 + 3*MTop*reMUE*SB*UStop12

	Sf1732 = CA2*ML*Sf1729 + SA*Sf1731 + 
     &   ML*SB*Re(Atau)*UStau21 - 
     &   2*CA*CB*MZ2*Sf127*SW2*UStau22

	Sf1558 = Sf1557*MHiggs2(4) - 7*Sf1553*MSneu2(1)**2 + 
     &   Sf680*MSneu2(1)**3 + Sf1554*MStau2(2)**2

	Sf147 = 3*CA2*MBot*SB*Sf144 + SA*Sf146 - 
     &   CA*CB*Sf127*Sf15*USbot21 + 3*CB*MBot*reMUE*USbot22

	Sf143 = 3*CA2*MBot*Sf140 + SA*Sf142 + 
     &   3*MBot*SB*Re(Abot)*USbot21 - 
     &   2*CA*CB*MZ2*Sf127*SW2*USbot22

	Sf136 = 3*CA2*MBot*SB*Sf133 + SA*Sf135 - 
     &   CA*CB*Sf127*Sf15*USbot11 + 3*CB*MBot*reMUE*USbot12

	Sf132 = 3*CA2*MBot*Sf128 + SA*Sf131 + 
     &   3*MBot*SB*Re(Abot)*USbot11 - 
     &   2*CA*CB*MZ2*Sf127*SW2*USbot12

	Sf1201 = 4*ML2*MW2**2 + MW2**3*Sf190 + 
     &   Sf176*Sf191*MHiggs2(4) - 2*Sf192*Ss25*MStau2(2)

	Sf193 = 4*ML2*MW2**2 + MW2**3*Sf190 + 
     &   Sf176*Sf191*MHiggs2(4) - 2*Sf192*Ss25*MStau2(1)

	Sf98 = Sf96 - Sqrt(Sf97)

	Sf90 = Sf88 - Sqrt(Sf89)

	Sf76 = Sf74 - Sqrt(Sf75)

	Sf68 = Sf66 - Sqrt(Sf67)

	Sf46 = Sf44 - Sqrt(Sf45)

	Sf33 = Sf31 - Sqrt(Sf32)

	Sf5513 = Sf3371 - 3*MStop2(2)

	Sf3135 = Sf3371 - 3*MStop2(1)

	Sf8870 = 31076 - 2*CW2*Sf8869 - 
     &   6*USbot22**2*(-319 + Sf8866 + 241*USbot22**2) - 
     &   5778*UStau22**2 + 11565*UStau22**4

	Sf8594 = -3*MStau2(2)*UStau22**2 + 
     &   3*MStau2(2)*UStau22**4 - 
     &   2*Sf8534*SW2**2*(-1 + UStau22**2)**2

	Sf7492 = -25*MSbot2(2) + Sf7488*MStop2(1) + 
     &   Sf7491*MStop2(2) + MSbot2(3)*(-16 + 57*USbot22**2)

	tmp35 = Sf55**2*Sf57**4*Sf9113 - 2*Sf57**2*Sf9169*MStau2(4)

	Sf10118 = Sf55**2*tmp35 - Sf9167*MStau2(1) + 
     &   Sf9168*MStau2(1)**2 + MZ2**3*Sf10113*MStau2(2) + 
     &   MHiggs2(1)*(-5*MZ2**2*Sf10114*MStau2(2) + 
     &      MHiggs2(1)*(MZ2*Sf10117 - Sf10116*MHiggs2(1)*MStau2(2))
     &      )

	Sf10110 = Sf55**2*tmp35 + MZ2**3*Sf10105*MStau2(1) + 
     &   MHiggs2(1)*(-5*MZ2**2*Sf10106*MStau2(1) + 
     &      MHiggs2(1)*(MZ2*Sf10109 + Sf10108*MHiggs2(1)*MStau2(1))
     &      ) - Sf9167*MStau2(2) + Sf9168*MStau2(2)**2

	tmp34 = Sf20**2*Sf22**4*Sf9113 - 2*Sf22**2*Sf9169*MSbot2(4)

	Sf9179 = Sf20**2*tmp34 - Sf9167*MSbot2(1) + 
     &   Sf9168*MSbot2(1)**2 + MZ2**3*Sf9176*MSbot2(2) + 
     &   MHiggs2(1)*(-5*MZ2**2*Sf8990*MSbot2(2) + 
     &      MHiggs2(1)*(MZ2*Sf9178 - Sf8992*MHiggs2(1)*MSbot2(2)))

	Sf9173 = Sf20**2*tmp34 + MZ2**3*Sf9170*MSbot2(1) + 
     &   MHiggs2(1)*(-5*MZ2**2*Sf8981*MSbot2(1) + 
     &      MHiggs2(1)*(MZ2*Sf9172 + Sf8983*MHiggs2(1)*MSbot2(1)))-
     &     Sf9167*MSbot2(2) + Sf9168*MSbot2(2)**2

	tmp38 = Sf107**2*Sf109**4*Sf9113 - 2*Sf109**2*Sf9169*MStop2(4)

	Sf10302 = Sf107**2*tmp38 - Sf9167*MStop2(1) + 
     &   Sf9168*MStop2(1)**2 + MZ2**3*Sf10297*MStop2(2) + 
     &   MHiggs2(1)*(-5*MZ2**2*Sf10298*MStop2(2) + 
     &      MHiggs2(1)*(MZ2*Sf10301 - Sf10300*MHiggs2(1)*MStop2(2))
     &      )

	Sf10292 = Sf107**2*tmp38 + MZ2**3*Sf10287*MStop2(1) + 
     &   MHiggs2(1)*(-5*MZ2**2*Sf10288*MStop2(1) + 
     &      MHiggs2(1)*(MZ2*Sf10291 + Sf10290*MHiggs2(1)*MStop2(1))
     &      ) - Sf9167*MStop2(2) + Sf9168*MStop2(2)**2

	tmp33 = Sf20**2*Sf22**4*Sf8911 - 2*Sf22**2*Sf8979*MSbot2(4)

	Sf8985 = Sf20**2*tmp33 + MZ2**3*Sf8980*MSbot2(1) + 
     &   MHiggs2(2)*(-5*MZ2**2*Sf8981*MSbot2(1) + 
     &      MHiggs2(2)*(MZ2*Sf8984 + Sf8983*MHiggs2(2)*MSbot2(1)))-
     &     Sf8977*MSbot2(2) + Sf8978*MSbot2(2)**2

	Sf8994 = Sf20**2*tmp33 - Sf8977*MSbot2(1) + 
     &   Sf8978*MSbot2(1)**2 + MZ2**3*Sf8989*MSbot2(2) + 
     &   MHiggs2(2)*(-5*MZ2**2*Sf8990*MSbot2(2) + 
     &      MHiggs2(2)*(MZ2*Sf8993 - Sf8992*MHiggs2(2)*MSbot2(2)))

	tmp37 = Sf107**2*Sf109**4*Sf8911 - 2*Sf109**2*Sf8979*MStop2(4)

	Sf10494 = Sf107**2*tmp37 - Sf8977*MStop2(1) + 
     &   Sf8978*MStop2(1)**2 + MZ2**3*Sf10491*MStop2(2) + 
     &   MHiggs2(2)*(-5*MZ2**2*Sf10298*MStop2(2) + 
     &      MHiggs2(2)*(MZ2*Sf10493 - Sf10300*MHiggs2(2)*MStop2(2))
     &      )

	Sf10488 = Sf107**2*tmp37 + MZ2**3*Sf10485*MStop2(1) + 
     &   MHiggs2(2)*(-5*MZ2**2*Sf10288*MStop2(1) + 
     &      MHiggs2(2)*(MZ2*Sf10487 + Sf10290*MHiggs2(2)*MStop2(1))
     &      ) - Sf8977*MStop2(2) + Sf8978*MStop2(2)**2

	Sf6894 = 30*MBot2*MTop2*MSbot2(2) + 
     &   MSbot2(2)**2*(-3*Sf4589*Ss30 - 645*C2B*MSbot2(2)) - 
     &   Sf4691*MStop2(2) - Sf4692*MStop2(2)**2 - 
     &   915*C2B*MStop2(2)**3

	Sf4693 = 30*MBot2*MTop2*MSbot2(2) + 
     &   MSbot2(2)**2*(-3*Sf4589*Ss30 - 645*C2B*MSbot2(2)) - 
     &   Sf4691*MStop2(1) - Sf4692*MStop2(1)**2 - 
     &   915*C2B*MStop2(1)**3

	Sf5414 = 10*MBot2*MTop2*MSbot2(2) + 
     &   MSbot2(2)**2*(Sf2875*Ss30 - 47*C2B*MSbot2(2)) - 
     &   Sf3009*MStop2(2) - Sf3010*MStop2(2)**2 - 
     &   89*C2B*MStop2(2)**3

	Sf3011 = 10*MBot2*MTop2*MSbot2(2) + 
     &   MSbot2(2)**2*(Sf2875*Ss30 - 47*C2B*MSbot2(2)) - 
     &   Sf3009*MStop2(1) - Sf3010*MStop2(1)**2 - 
     &   89*C2B*MStop2(1)**3

	Sf6643 = 30*MBot2*MTop2*MSbot2(1) + 
     &   MSbot2(1)**2*(-3*Sf4589*Ss30 - 645*C2B*MSbot2(1)) - 
     &   Sf4590*MStop2(2) - Sf4592*MStop2(2)**2 - 
     &   915*C2B*MStop2(2)**3

	Sf4593 = 30*MBot2*MTop2*MSbot2(1) + 
     &   MSbot2(1)**2*(-3*Sf4589*Ss30 - 645*C2B*MSbot2(1)) - 
     &   Sf4590*MStop2(1) - Sf4592*MStop2(1)**2 - 
     &   915*C2B*MStop2(1)**3

	Sf5308 = 10*MBot2*MTop2*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf2875*Ss30 - 47*C2B*MSbot2(1)) - 
     &   Sf2877*MStop2(2) - Sf2879*MStop2(2)**2 - 
     &   89*C2B*MStop2(2)**3

	Sf2880 = 10*MBot2*MTop2*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf2875*Ss30 - 47*C2B*MSbot2(1)) - 
     &   Sf2877*MStop2(1) - Sf2879*MStop2(1)**2 - 
     &   89*C2B*MStop2(1)**3

	Sf648 = Sf647*MHiggs2(4)**2 + Sf644*MSneu2(1)**2 - 
     &   2*ML2*MStau2(1)**2 - Sf642*MStau2(1)**3 + 
     &   MSneu2(1)*(-5*Sf640*MSneu2(1)**2 - Sf641*MStau2(1)**2)

	Sf1535 = Sf1534*MHiggs2(4)**2 + Sf1533*MSneu2(1)**2 - 
     &   2*ML2*MStau2(2)**2 - Sf642*MStau2(2)**3 + 
     &   MSneu2(1)*(-5*Sf640*MSneu2(1)**2 - Sf641*MStau2(2)**2)

	Sf7643 = 32*MW2**4 - 88*MW2**3*MZ2 + 92*MW2**2*MZ2**2 - 
     &   33*MW2*MZ2**3 + MZ2*
     &    (Sf7639*UStau22**2 + Sf7642*(1 - UStau22**2))

	Sf8737 = MW2*(16*MW2*MZ2**2 - 24*MZ2**3) + 
     &   2*Sf7756*Sf8736*MStop2(2)*UStop22**2 + 
     &   MZ2*(11*MZ2**3 - 2*Sf7756*MStop2(2)*UStop22**4)

	Sf7766 = MW2*(16*MW2*MZ2**2 - 24*MZ2**3) + 
     &   2*Sf7756*Sf7765*MStau2(2)*UStau22**2 + 
     &   MZ2*(11*MZ2**3 - 2*Sf7756*MStau2(2)*UStau22**4)

	Sf8505 = MW2*(16*MW2*MZ2**2 - 24*MZ2**3) + 
     &   2*Sf7756*Sf8504*MSbot2(2)*USbot22**2 + 
     &   MZ2*(11*MZ2**3 - 2*Sf7756*MSbot2(2)*USbot22**4)

	Sf630 = 6*ML2*MW2 - 6*C2B*MW2**2 + 6*Sf629*MSneu2(1) - 
     &   3*Sf198*MStau2(1) - 4*MStau2(1)**2

	Sf1524 = 6*ML2*MW2 - 6*C2B*MW2**2 + 6*Sf1523*MSneu2(1) - 
     &   3*Sf198*MStau2(2) - 4*MStau2(2)**2

	Sf8758 = -81*MZ2**6*Sf8741*UStop11 + 
     &   MW2*(MW2*(MW2*(MW2*
     &             (48*MZ2**2*Sf8745 + 
     &               MW2*(-64*MZ2*Sf8747 + 256*MW2*Sf8749)) + 
     &            4*MZ2**3*Sf8752) - Sf8757) + 
     &      9*MZ2**5*Sf8743*UStop11)

	Sf10035 = MZ2*(-MZ2 + MHiggs2(2))*(MZ2 - MStau2(1))**2*
     &    MStau2(1)**2*(MZ2 + 4*MStau2(1)) + 
     &   MStau2(1)*(25*MZ2**4*(-MZ2 + MHiggs2(2)) - 
     &      36*MZ2**3*(-4*MZ2 + 5*MHiggs2(2))*MStau2(1) + 
     &      7*MZ2**2*(-13*MZ2 + 29*MHiggs2(2))*MStau2(1)**2 + 
     &      2*MZ2*(-111*MZ2 + 101*MHiggs2(2))*MStau2(1)**3 + 
     &      (152*MZ2 - 144*MHiggs2(2))*MStau2(1)**4)*MStau2(2) + 
     &   MStau2(1)*(-46*MZ2**4 - 54*MZ2**3*MHiggs2(2) + 
     &      3*MZ2**2*(-19*MZ2 + 119*MHiggs2(2))*MStau2(1) + 
     &      4*MZ2*(MZ2 - 107*MHiggs2(2))*MStau2(1)**2 + 
     &      (-604*MZ2 + 596*MHiggs2(2))*MStau2(1)**3)*MStau2(2)**2+
     &     (8*MZ2**4 - 24*MZ2**3*MHiggs2(2) + 
     &      MZ2**2*(193*MZ2 + 139*MHiggs2(2))*MStau2(1) + 
     &      2*MZ2*(129*MZ2 - 67*MHiggs2(2))*MStau2(1)**2 + 
     &      (936*MZ2 - 904*MHiggs2(2))*MStau2(1)**3)*MStau2(2)**3+
     &     (24*MZ2**3 - 8*MZ2**2*MHiggs2(2) + 
     &      4*MZ2*(-81*MZ2 + 5*MHiggs2(2))*MStau2(1) + 
     &      (-684*MZ2 + 596*MHiggs2(2))*MStau2(1)**2)*MStau2(2)**4+
     &     8*(-2*MZ2**2 + 4*MZ2*MHiggs2(2) + 
     &      (27*MZ2 - 18*MHiggs2(2))*MStau2(1))*MStau2(2)**5 - 
     &   16*MZ2*MStau2(2)**6

	Sf7283 = Sf7240*MSbot2(1) + Sf7241*MSbot2(2) + 
     &   3*Sf7244*MStau2(1) + 3*Sf7248*MStau2(2) + 
     &   Sf7253*MStop2(1) + Sf7271*MStop2(2)

	Sf8859 = Sf11*Sf16**2*Sf8764 + Sf26*Sf30**2*Sf8769 + 
     &   3*Sf39*Sf43**2*Sf8776 + 3*Sf61*Sf65**2*Sf8787 + 
     &   Sf82*Sf87**2*Sf8800 + Sf113*Sf117**2*Sf8817

	Sf1159 = 12 - 2*CB2*SA2 - 2*CAB*SA*SB + 5*CA2*SB2 - 
     &   3*SA*SB*Sf1156 - CA*Sf1158

	Sf1165 = 6 - CB2*SA2 - CAB*SA*SB - 2*CA2*SB2 + 
     &   3*SA*SB*Sf1156 - CA*Sf1164

	Sf9221 = Sf9067*Sf9191 - 2*Sf9217*MSbot2(1)**3 + 
     &   MZ2**2*(MZ2*Sf9215*MSbot2(1) - Sf9219*MSbot2(1)) + 
     &   MHiggs2(1)**2*(Sf9220 + 3*Sf8983*MHiggs2(1)*MSbot2(1)) - 
     &   Sf9182*MSbot2(2) + Sf9183*MSbot2(2)**2 - 
     &   Sf9184*MSbot2(2)**3

	Sf9193 = Sf9058*Sf9191 - Sf9182*MSbot2(1) + 
     &   Sf9183*MSbot2(1)**2 - Sf9184*MSbot2(1)**3 - 
     &   2*Sf9187*MSbot2(2)**3 + 
     &   MZ2**2*(MZ2*Sf9185*MSbot2(2) - Sf9190*MSbot2(2)) + 
     &   MHiggs2(1)**2*(Sf9192 - 3*Sf8992*MHiggs2(1)*MSbot2(2))

	Sf10157 = MHiggs2(1)**3*
     &    (3*MZ2**2 - 5*MZ2*MStau2(2) - 3*MStau2(1)*MStau2(2) + 
     &      3*MStau2(2)**2 + MStau2(1)*(-7*MZ2 - 3*MStau2(3))) + 
     &   4*MZ2**3*MStau2(1)*MStau2(3) - 
     &   4*MZ2*MStau2(3)*(MZ2**2*MStau2(2) + 
     &      2*MZ2*(MStau2(1)**2 - MStau2(2)**2) + MStau2(3)**3) + 
     &   MHiggs2(1)*(MZ2**3*(7*MStau2(1) + 5*MStau2(2)) - 
     &      MZ2**2*(7*MStau2(2)**2 + 
     &         MStau2(1)*(11*MStau2(1) + 30*MStau2(2))) + 
     &      3*MStau2(3)**4) - 
     &   MHiggs2(1)**2*(6*MStau2(2)**3 + 
     &      MZ2*(3*MZ2**2 + 
     &         MStau2(1)*(-11*MStau2(1) - 30*MStau2(2)) - 
     &         7*MStau2(2)**2) + 
     &      MStau2(1)*(6*MStau2(1)**2 - 6*MStau2(2)*MStau2(4)))

	Sf10128 = -5*MZ2*MHiggs2(1)*(-MZ2 + MHiggs2(1))*
     &    (MZ2 + MHiggs2(1))*MStau2(1) + 
     &   (3*MHiggs2(1)**3 + 
     &      MZ2*(MZ2*(-4*MZ2 - 7*MHiggs2(1)) + 7*MHiggs2(1)**2))*
     &    MStau2(1)**2 + (8*MZ2**2 - 6*MHiggs2(1)**2)*
     &    MStau2(1)**3 - 2*
     &    (3*MHiggs2(1)**2 + 6*MHiggs2(1)*MStau2(1) - 
     &      4*MZ2*(MZ2 + 2*MStau2(1)))*MStau2(2)**3 + 
     &   (-4*MZ2 + 3*MHiggs2(1))*
     &    (6*MStau2(1)**2*MStau2(2)**2 + MStau2(2)**4 - 
     &      MStau2(1)**3*(3*MStau2(2) + MStau2(3))) - 
     &   MZ2**2*MStau2(2)*(MHiggs2(1)*
     &       (-7*MZ2 + 30*MStau2(1) + 11*MStau2(2)) - 
     &      4*MZ2*(MStau2(1) - MStau2(3)) + 8*MStau2(1)*MStau2(4))+
     &     MHiggs2(1)**2*(-3*MZ2**3 + 3*MZ2**2*MHiggs2(1) + 
     &      11*MZ2*MStau2(2)**2 + 
     &      MHiggs2(1)*MStau2(2)*(-7*MZ2 + 3*MStau2(3)) + 
     &      MStau2(1)*MStau2(2)*
     &       (30*MZ2 - 3*MHiggs2(1) + 6*MStau2(4)))

	Sf7700 = 21*MW2 - 6*Sf7699*Ss398 - 
     &   (24*MW2**3*Sf7497*MSneu2(1))/(CW2**2*MZ2**3) + 
     &   ((Sf7698 + 6*Sf7661*Ss949*MStau2(2) + 
     &        (Sf7691*MStau2(2))/MSl2Diff(2,1)**2)*UStau21**2)/
     &    Sf1195**3 + ((Sf7684 + 6*Sf7649*Ss928*MStau2(1) + 
     &        (Sf7677*MStau2(1))/MSl2Diff(1,1)**2)*UStau22**2)/
     &    Sf166**3

	Sf668 = MW2*(2*MW2**2 - Sf663*MHiggs2(4)) + 
     &   Sf665*MSneu2(1)**2 + Sf659*MSneu2(1)**3 - 
     &   Sf667*MSneu2(1)*MStau2(1) + 5*ML2*MStau2(1)**2 - 
     &   Sf660*MStau2(1)**3

	Sf1542 = MW2*(2*MW2**2 - Sf1539*MHiggs2(4)) + 
     &   Sf1540*MSneu2(1)**2 + Sf659*MSneu2(1)**3 - 
     &   Sf1541*MSneu2(1)*MStau2(2) + 5*ML2*MStau2(2)**2 - 
     &   Sf660*MStau2(2)**3

	Sf5230 = reMUE*Sf5229 + 
     &   3*MBot*SA*SB*Sf5224*Re(Abot) + 6*MBot2*SA*SB*MSbot2(1) + 
     &   CB*(2*MW2 - 2*MZ2)*SAB*SB*MSbot2(1) + 
     &   6*MBot2*SA*SB*MSbot2(2) + 
     &   CB*(-2*MW2 - MZ2)*SAB*SB*MSbot2(2) + 
     &   CB*MZ2*SAB*SB*MSneu2(1) + 2*ML2*SA*SB*MStau2(2) + 
     &   CB*(-2*MW2 + MZ2)*SAB*SB*MStau2(2) - 
     &   6*CA*CB*MTop2*MStop2(1) + 
     &   CB*(-4*MW2 + 4*MZ2)*SAB*SB*MStop2(1) - 
     &   6*CA*CB*MTop2*MStop2(2) + 
     &   CB*(4*MW2 - MZ2)*SAB*SB*MStop2(2) - 
     &   SB*MSbot2(3)*(-3*MBot*(CA*reMUE + SA*Re(Abot))*
     &       USbot21*USbot22 - 
     &      CB*(4*MW2 - MZ2)*SAB*USbot22**2) + 
     &   SB*Sf1746*MStau2(1)*UStau11 + 
     &   SB*Sf1745*MStau2(1)*UStau12 + 
     &   ML*SB*(CA*reMUE + 2*SA*Re(Re(Atau)))*MStau2(2)*
     &    UStau21*UStau22 + 
     &   CB*(4*MW2 - 3*MZ2)*SAB*SB*MStau2(2)*UStau22**2 - 
     &   CB*MStop2(3)*(3*MTop*(reMUE*SA + 2*CA*Re(Re(Atop)))*
     &       UStop21*UStop22 + 
     &      (8*MW2 - 5*MZ2)*SAB*SB*UStop22**2)

	Sf5227 = -(reMUE*Sf5226) + 
     &   3*CA*MBot*SB*Sf5224*Re(Abot) + 6*CA*MBot2*SB*MSbot2(1) + 
     &   CAB*CB*(2*MW2 - 2*MZ2)*SB*MSbot2(1) + 
     &   6*CA*MBot2*SB*MSbot2(2) + 
     &   CAB*CB*(-2*MW2 - MZ2)*SB*MSbot2(2) + 
     &   CAB*CB*MZ2*SB*MSneu2(1) + 2*CA*ML2*SB*MStau2(2) + 
     &   CAB*CB*(-2*MW2 + MZ2)*SB*MStau2(2) + 
     &   6*CB*MTop2*SA*MStop2(1) + 
     &   CAB*CB*(-4*MW2 + 4*MZ2)*SB*MStop2(1) + 
     &   6*CB*MTop2*SA*MStop2(2) + 
     &   CAB*CB*(4*MW2 - MZ2)*SB*MStop2(2) - 
     &   SB*MSbot2(3)*(3*MBot*(reMUE*SA - CA*Re(Abot))*
     &       USbot21*USbot22 - 
     &      CAB*CB*(4*MW2 - MZ2)*USbot22**2) + 
     &   SB*Sf1752*MStau2(1)*UStau11 - 
     &   SB*Sf1751*MStau2(1)*UStau12 + 
     &   ML*SB*(-(reMUE*SA) + 2*CA*Re(Re(Atau)))*MStau2(2)*
     &    UStau21*UStau22 + 
     &   CAB*CB*(4*MW2 - 3*MZ2)*SB*MStau2(2)*UStau22**2 - 
     &   CB*MStop2(3)*(3*MTop*(CA*reMUE - 2*SA*Re(Re(Atop)))*
     &       UStop21*UStop22 + 
     &      CAB*(8*MW2 - 5*MZ2)*SB*UStop22**2)

	Sf7404 = Sf7403

	Sf8759 = UStop12*
     &   (4*MZ2*Sf7204*SW2 - Sf7199*Sf86*UStop22**2)

	Sf8548 = UStau12*
     &   (2*MZ2*Sf8547*SW2 + Sf8536*Ss443*UStau22**2)

	Sf7785 = UStau12*
     &   (2*MZ2*Sf1760*SW2 - Sf1755*Ss443*UStau22**2)

	Sf8487 = USbot12*
     &   (2*MZ2*Sf4778*SW2 - Sf15*Sf4773*USbot22**2)

	Sf6181 = MW2**2*Sf6180*MSq2Diff(2,2) + 
     &   (MW2**3*(MW2*(MW2*(4*MW2 - Sf6154) - Sf6176) + Sf6179) - 
     &      4*MW2*Sf6167*MSq2Diff(2,2)**4 + 24*MSq2Diff(2,2)**6)*
     &    USbot21

	Sf3936 = MW2**2*Sf3935*MSq2Diff(2,1) + 
     &   (MW2**3*(MW2*(MW2*(4*MW2 - Sf3906) - Sf3931) + Sf3934) - 
     &      4*MW2*Sf3919*MSq2Diff(2,1)**4 + 24*MSq2Diff(2,1)**6)*
     &    USbot21

	Sf5803 = MW2**2*Sf5802*MSq2Diff(1,2) + 
     &   (MW2**3*(MW2*(MW2*(4*MW2 - Sf5775) - Sf5798) + Sf5801) - 
     &      4*MW2*Sf5789*MSq2Diff(1,2)**4 + 24*MSq2Diff(1,2)**6)*
     &    USbot11

	Sf3812 = MW2**2*Sf3811*MSq2Diff(1,1) + 
     &   (MW2**3*(MW2*(MW2*(4*MW2 - Sf3772) - Sf3807) + Sf3810) - 
     &      4*MW2*Sf3794*MSq2Diff(1,1)**4 + 24*MSq2Diff(1,1)**6)*
     &    USbot11

	Sf572 = 5*MSl2Diff(1,1)**7 + 
     &   MW2*(-(Sf565*MSl2Diff(1,1)**5) + 
     &      MW2*(10*Sf568*MSl2Diff(1,1)**3 + 
     &         MW2*(Sf571 + 
     &            MW2*(-Sf570 + MW2*(Sf567 + 2*MW2*MStau2(1))))))

	Sf1473 = 5*MSl2Diff(2,1)**7 + 
     &   MW2*(-(Sf1468*MSl2Diff(2,1)**5) + 
     &      MW2*(10*Sf1470*MSl2Diff(2,1)**3 + 
     &         MW2*(Sf1472 + 
     &            MW2*(-Sf1471 + MW2*(Sf1469 + 2*MW2*MStau2(2))))))

	Sf7761 = MW2*(-(MZ2**4*Sf7759) + 
     &      MW2*(MW2*(-268*MZ2**3 + 
     &            MW2*(MW2*(192*MW2 - 576*MZ2) + 588*MZ2**2)) + 
     &         Sf7760)) + 9*MZ2**5*MSl2Diff(2,1)

	Sf6913 = 120*MSq2Diff(2,2)**6 + 
     &   MW2*(-12*Sf6888*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf6906) - Sf6908) + 
     &            3*Sf6912) - 3*Sf6910*MSq2Diff(2,2)))

	Sf6662 = 120*MSq2Diff(1,2)**6 + 
     &   MW2*(-12*Sf6637*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf6655) - Sf6657) + 
     &            3*Sf6661) - 3*Sf6659*MSq2Diff(1,2)))

	Sf4720 = 120*MSq2Diff(2,1)**6 + 
     &   MW2*(-12*Sf4683*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf4710) - Sf4713) + 
     &            3*Sf4719) - 3*Sf4717*MSq2Diff(2,1)))

	Sf4625 = 120*MSq2Diff(1,1)**6 + 
     &   MW2*(-12*Sf4578*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf4613) - Sf4618) + 
     &            3*Sf4624) - 3*Sf4622*MSq2Diff(1,1)))

	Sf6172 = 24*MSq2Diff(2,2)**6 + 
     &   MW2*(-4*Sf6167*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf6154) + Sf6169) - 
     &            Sf6171) + Sf6161*MSq2Diff(2,2)))

	Sf5794 = 24*MSq2Diff(1,2)**6 + 
     &   MW2*(-4*Sf5789*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf5775) + Sf5791) - 
     &            Sf5793) + Sf5782*MSq2Diff(1,2)))

	Sf5430 = -8*MSq2Diff(2,2)**6 + 
     &   MW2*(4*Sf5408*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf5424) + Sf5426) - 
     &            Sf5429) + Sf5428*MSq2Diff(2,2)))

	Sf5324 = -8*MSq2Diff(1,2)**6 + 
     &   MW2*(4*Sf5302*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf5318) + Sf5320) - 
     &            Sf5323) + Sf5322*MSq2Diff(1,2)))

	Sf3926 = 24*MSq2Diff(2,1)**6 + 
     &   MW2*(-4*Sf3919*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf3906) + Sf3922) - 
     &            Sf3925) + Sf3913*MSq2Diff(2,1)))

	Sf3802 = 24*MSq2Diff(1,1)**6 + 
     &   MW2*(-4*Sf3794*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf3772) + Sf3798) - 
     &            Sf3801) + Sf3785*MSq2Diff(1,1)))

	Sf3034 = -8*MSq2Diff(2,1)**6 + 
     &   MW2*(4*Sf3000*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf3024) + Sf3027) - 
     &            Sf3033) + Sf3031*MSq2Diff(2,1)))

	Sf2907 = -8*MSq2Diff(1,1)**6 + 
     &   MW2*(4*Sf2862*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf2897) + Sf2900) - 
     &            Sf2906) + Sf2904*MSq2Diff(1,1)))

	Sf1966 = -5*MSl2Diff(2,1)**6 + 
     &   MW2*(Sf1961*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1962 - Sf1964) + Sf1965) - 
     &         Sf1963*MSl2Diff(2,1)**2))

	Sf1906 = -5*MSl2Diff(1,1)**6 + 
     &   MW2*(Sf1899*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1900 - Sf1902) + Sf1905) - 
     &         Sf1901*MSl2Diff(1,1)**2))

	Sf5619 = 5*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2 - Sf5614) - 2*Sf5616) + Sf5618) - 
     &      Sf5615*MSq2Diff(1,2)**3)

	Sf599 = MSl2Diff(1,1)**5 + 
     &   MW2*(-(Sf595*MSl2Diff(1,1)**3) + 
     &      MW2*(-Sf598 + MW2*(Sf596 + MW2*MStau2(1))))

	Sf3638 = Sf3637*MSbot2(2)**4 + 
     &   MW2*(20*Sf3618*MSbot2(2)**3 + 
     &      MW2*(10*Sf3633*MSbot2(2)**2 + 
     &         MW2*(MW2*Sf3636 + 2*Sf3634*MSbot2(2))))

	Sf3337 = Sf3336*MSbot2(1)**4 + 
     &   MW2*(20*Sf3294*MSbot2(1)**3 + 
     &      MW2*(10*Sf3326*MSbot2(1)**2 + 
     &         MW2*(MW2*Sf3335 + 2*Sf3329*MSbot2(1))))

	Sf3563 = -14*MSbot2(2)**5 + 
     &   MW2*(76*MSbot2(2)**4 + 
     &      MW2*(-86*MSbot2(2)**3 + 
     &         MW2*(6*MW2*Sf3206 + 3*Sf3562*MSbot2(2))))

	Sf3209 = -14*MSbot2(1)**5 + 
     &   MW2*(76*MSbot2(1)**4 + 
     &      MW2*(-86*MSbot2(1)**3 + 
     &         MW2*(6*MW2*Sf3206 + 3*Sf3208*MSbot2(1))))

	Sf2122 = ML2*MSl2Diff(2,1)**3 + 
     &   MW2*(MW2*(-Sf2118 + MW2*(Sf2115 - MW2*Sf991)) - 
     &      Sf2121*MSl2Diff(2,1))

	Sf213 = MHiggs2(4)*
     &    (-(Sf176*Sf212) - 12*ML2*MW2*MStau2(1)) + 
     &   MW2**2*(MW2*(MW2*Sf206 + Sf210) + 2*Sf208*MStau2(1))

	Sf2040 = ML2*MSl2Diff(1,1)**3 + 
     &   MW2*(MW2*(-Sf2036 + MW2*(Sf2032 - MW2*Sf991)) - 
     &      Sf2039*MSl2Diff(1,1))

	Sf719 = ML2*MSl2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf197 - Sf711) + Sf715) - 
     &      Sf718*MSl2Diff(1,1))

	Sf729 = ML2*MSl2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf169 - Sf722) + Sf728) - 
     &      Sf725*MSl2Diff(1,1))

	Sf1590 = -(ML2*MSl2Diff(2,1)**3) + 
     &   MW2*(MW2*(-Sf1589 + MW2*(Sf1583 - MW2*Sf169)) + 
     &      Sf1586*MSl2Diff(2,1))

	Sf1582 = -(ML2*MSl2Diff(2,1)**3) + 
     &   MW2*(MW2*(-Sf1578 + MW2*(Sf1576 - MW2*Sf197)) + 
     &      Sf1581*MSl2Diff(2,1))

	Sf1499 = MSl2Diff(2,1)**5 + 
     &   MW2*(-(Sf1495*MSl2Diff(2,1)**3) + 
     &      MW2*(-Sf1498 + MW2*(Sf1496 + MW2*MStau2(2))))

	Sf1209 = MHiggs2(4)*
     &    (-(Sf1208*Sf176) - 12*ML2*MW2*MStau2(2)) + 
     &   MW2**2*(MW2*(Sf1206 + MW2*Sf206) + 2*Sf1205*MStau2(2))

	Sf8810 = MW2*(-(MZ2*Sf8808) + 
     &      MW2*(MW2*(8*MW2 - Sf8807) + Sf8809)) + 
     &   9*MZ2**2*MSq2Diff(2,2)**2

	Sf8804 = MW2*(-(MZ2*Sf8802) + 
     &      MW2*(MW2*(8*MW2 - Sf8801) + Sf8803)) + 
     &   9*MZ2**2*MSq2Diff(1,2)**2

	Sf8797 = MW2*(-(MZ2*Sf8795) + 
     &      MW2*(MW2*(8*MW2 - Sf8794) + Sf8796)) + 
     &   9*MZ2**2*MSq2Diff(2,1)**2

	Sf8791 = MW2*(-(MZ2*Sf8789) + 
     &      MW2*(MW2*(8*MW2 - Sf8788) + Sf8790)) + 
     &   9*MZ2**2*MSq2Diff(1,1)**2

	Sf8780 = MW2*(-(MZ2*Sf8778) + 
     &      MW2*(MW2*(8*MW2 - Sf8777) + Sf8779)) + 
     &   9*MZ2**2*MSl2Diff(2,1)**2

	Sf8773 = MW2*(-(MZ2*Sf8771) + 
     &      MW2*(MW2*(8*MW2 - Sf8770) + Sf8772)) + 
     &   9*MZ2**2*MSl2Diff(1,1)**2

	Sf5882 = 5*MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(MW2 - Sf5877) + Sf5881) - 
     &      Sf5879*MSq2Diff(1,2))

	Sf8581 = MW2*(MW2*(8*MW2*Sf8540 - 12*MZ2*Sf8577) + 
     &      Sf8580) + 3*MZ2**3*MSneu2(1)*UStau11

	Sf1977 = -(Sf1976*MStau2(2)**2) + 
     &   MSneu2(1)*(-8*MStau2(2)**2 + 
     &      MSneu2(1)*(50*MSneu2(1) - Sf1903*MStau2(2)))

	Sf1917 = -(Sf1916*MStau2(1)**2) + 
     &   MSneu2(1)*(-8*MStau2(1)**2 + 
     &      MSneu2(1)*(50*MSneu2(1) - Sf1903*MStau2(1)))

	Sf1483 = -45*MStau2(2)**3 + 
     &   MSneu2(1)*(-(Sf1482*MStau2(2)) + 
     &      MSneu2(1)*(25*MSneu2(1) - Sf569*MStau2(2)))

	Sf582 = -45*MStau2(1)**3 + 
     &   MSneu2(1)*(-(Sf581*MStau2(1)) + 
     &      MSneu2(1)*(25*MSneu2(1) - Sf569*MStau2(1)))

	Sf3567 = -70*MSbot2(2)**4 + 
     &   MW2*(20*MSbot2(2)**3 + 
     &      MW2*(-(MW2*Sf3566) + 108*MSbot2(2)**2))

	Sf3215 = -70*MSbot2(1)**4 + 
     &   MW2*(20*MSbot2(1)**3 + 
     &      MW2*(-(MW2*Sf3214) + 108*MSbot2(1)**2))

	Sf8546 = -(MW2*(MW2*(-3*MZ2*Sf8537 + 28*MW2*Sf8540) + 
     &        Sf8545)) + 21*MZ2**3*MSneu2(1)*UStau11

	Sf3628 = Sf3627*MSbot2(2)**3 + 
     &   MW2*(Sf3625*MSbot2(2)**2 + 
     &      MW2*(2*MW2*Sf3624 + 8*Sf3623*MSbot2(2)))

	Sf3315 = Sf3314*MSbot2(1)**3 + 
     &   MW2*(Sf3310*MSbot2(1)**2 + 
     &      MW2*(2*MW2*Sf3307 + 8*Sf3304*MSbot2(1)))

	Sf283 = -(Sf275*MSl2Diff(1,1)**2) + 
     &   MHiggs2(4)*(Sf282 + 
     &      MHiggs2(4)*(Sf271*MHiggs2(4) + Sf269*MSl2Diff(1,1)))

	Sf2126 = -(Sf2124*MSl2Diff(2,1)**2) + 
     &   MW2*(MW2*(MW2*Sf169 - Sf2123) - Sf2125*MSl2Diff(2,1))

	Sf2046 = -(Sf2043*MSl2Diff(1,1)**2) + 
     &   MW2*(MW2*(MW2*Sf169 - Sf2042) - Sf2045*MSl2Diff(1,1))

	Sf1253 = -(Sf1249*MSl2Diff(2,1)**2) + 
     &   MHiggs2(4)*(Sf1252 + 
     &      MHiggs2(4)*(Sf1246*MHiggs2(4) + Sf1245*MSl2Diff(2,1)))

	Sf9735 = Sf9734*MStau2(2)**3 + 
     &   MHiggs2(2)*(MZ2*Sf9434*MStau2(2)**2 + 
     &      MZ2**2*MHiggs2(2)*
     &       (336*MZ2*MHiggs2(2) - Sf9433*MStau2(2)))

	Sf9440 = Sf9439*MStau2(1)**3 + 
     &   MHiggs2(2)*(MZ2*Sf9434*MStau2(1)**2 + 
     &      MZ2**2*MHiggs2(2)*
     &       (336*MZ2*MHiggs2(2) - Sf9433*MStau2(1)))

	Sf9732 = -(Sf9731*MStau2(2)**3) + 
     &   MHiggs2(2)*(-10*MZ2*Sf9442*MStau2(2)**2 + 
     &      MZ2**2*MHiggs2(2)*
     &       (116*MZ2*MHiggs2(2) - Sf9441*MStau2(2)))

	Sf9448 = -(Sf9447*MStau2(1)**3) + 
     &   MHiggs2(2)*(-10*MZ2*Sf9442*MStau2(1)**2 + 
     &      MZ2**2*MHiggs2(2)*
     &       (116*MZ2*MHiggs2(2) - Sf9441*MStau2(1)))

	Sf9764 = Sf9763*MStau2(2)**3 + 
     &   MHiggs2(2)*(-2*MZ2*Sf9500*MStau2(2)**2 + 
     &      MZ2**2*MHiggs2(2)*
     &       (2*MZ2*MHiggs2(2) + 4*Sf8940*MStau2(2)))

	Sf9506 = Sf9505*MStau2(1)**3 + 
     &   MHiggs2(2)*(-2*MZ2*Sf9500*MStau2(1)**2 + 
     &      MZ2**2*MHiggs2(2)*
     &       (2*MZ2*MHiggs2(2) + 4*Sf8940*MStau2(1)))

	Sf9147 = -(Sf20**2*Sf22**2*Sf9113) + 
     &   MZ2*(MZ2*(3*MZ2 - Sf9145) + Sf9146)

	Sf9144 = -(Sf20**2*Sf22**2*Sf9113) + 
     &   MZ2*(MZ2*(3*MZ2 - Sf9141) + Sf9143)

	Sf9348 = -(Sf55**2*Sf57**2*Sf8911) + 
     &   MZ2*(MZ2*(3*MZ2 - Sf9346) + Sf9347)

	Sf9340 = -(Sf55**2*Sf57**2*Sf8911) + 
     &   MZ2*(MZ2*(3*MZ2 - Sf9338) + Sf9339)

	Sf8956 = -(Sf20**2*Sf22**2*Sf8911) + 
     &   MZ2*(MZ2*(3*MZ2 - Sf8954) + Sf8955)

	Sf8942 = -(Sf20**2*Sf22**2*Sf8911) + 
     &   MZ2*(MZ2*(3*MZ2 - Sf8939) + Sf8941)

	Sf10469 = MZ2*(MZ2*(3*MZ2 - Sf10467) + Sf10468) - 
     &   Sf107**2*Sf109**2*Sf8911

	Sf10459 = MZ2*(MZ2*(3*MZ2 - Sf10457) + Sf10458) - 
     &   Sf107**2*Sf109**2*Sf8911

	Sf10085 = MZ2*(MZ2*(3*MZ2 - Sf10083) + Sf10084) - 
     &   Sf55**2*Sf57**2*Sf9113

	Sf10082 = MZ2*(MZ2*(3*MZ2 - Sf10080) + Sf10081) - 
     &   Sf55**2*Sf57**2*Sf9113

	tmp31 = MZ2*(464*MZ2 + 466*MHiggs2(2))

	Sf9726 = MHiggs2(2)*(-288*MZ2**3 + tmp31*MHiggs2(2)) + 
     &   Sf9725*MStau2(2)

	Sf9462 = MHiggs2(2)*(-288*MZ2**3 + tmp31*MHiggs2(2)) + 
     &   Sf9461*MStau2(1)

	Sf5852 = -20*MStop2(2)**4 + 
     &   MW2*(26*MStop2(2)**3 + MW2*(8*MW2**2 + Sf5851*MStop2(2)))

	Sf4069 = -20*MStop2(1)**4 + 
     &   MW2*(26*MStop2(1)**3 + MW2*(8*MW2**2 + Sf4068*MStop2(1)))

	Sf679 = MW2*(MW2*(8*MW2 - 32*Sf34) + 16*Sf674) + 
     &   Sf678*MHiggs2(4)

	Sf1552 = MW2*(16*Sf1548 + MW2*(8*MW2 - 32*Sf47)) + 
     &   Sf1551*MHiggs2(4)

	tmp46 = 4*MW2**2 + MHiggs2(4)*(-3*Sf2228 + 5*MHiggs2(4))

	Sf2288 = tmp46*MHiggs2(4) - Sf2201*MSbot2(2)**2

	Sf2229 = tmp46*MHiggs2(4) - Sf2201*MSbot2(1)**2

	Sf165 = -2*MW2**2*Sf34 + 
     &   MHiggs2(4)*(Sf164 + MHiggs2(4)*(-4*Sf163 + 3*MHiggs2(4)))

	Sf1194 = -2*MW2**2*Sf47 + 
     &   MHiggs2(4)*(Sf1193 + 
     &      MHiggs2(4)*(-4*Sf1192 + 3*MHiggs2(4)))

	Sf6144 = MHiggs2(4)*
     &    (-Sf6143 + MHiggs2(4)*(MHiggs2(4) - 2*MSbot2(2))) + 
     &   Sf173*MSq2Diff(2,2)**2

	Sf3742 = MHiggs2(4)*
     &    (-Sf3741 + MHiggs2(4)*(MHiggs2(4) - 2*MSbot2(2))) + 
     &   Sf173*MSq2Diff(2,1)**2

	Sf5761 = MHiggs2(4)*
     &    (-Sf5760 + MHiggs2(4)*(MHiggs2(4) - 2*MSbot2(1))) + 
     &   Sf173*MSq2Diff(1,2)**2

	Sf3733 = MHiggs2(4)*
     &    (-Sf3732 + MHiggs2(4)*(MHiggs2(4) - 2*MSbot2(1))) + 
     &   Sf173*MSq2Diff(1,1)**2

	Sf6961 = 36*MTop2*Sf6956*TB*UStop21 + 
     &   Sf6960*TB*UStop21 - 36*MTop*Sf2245*Sf6956*UStop22

	Sf4897 = (Sf2309*TB + MTop2*USbot21)*UStop21 + 
     &   MTop*Sf2310*UStop22

	Sf4863 = (Sf2258*TB + MTop2*USbot11)*UStop21 + 
     &   MTop*Sf2259*UStop22

	Sf7205 = 4*SW2*UStop12*(Sf7204 - MZ2*UStop12**2) + 
     &   (Sf7199*Sf86*UStop11*UStop21*UStop22)/MZ2

	Sf8413 = Sf4906**3*Sf8411*MSq2Diff(2,2)**2*
     &    UStop22**2 + Sf2334*
     &    (Sf2334**2*Sf8412*MSq2Diff(2,1)**2*UStop21**2 + 
     &      Sf4906**3*Sf8410*MSq2Diff(2,2)**2*UStop22**2)

	Sf8317 = Sf4820**3*Sf8315*MSq2Diff(1,2)**2*
     &    UStop22**2 + Sf2320*
     &    (Sf2320**2*Sf8316*MSq2Diff(1,1)**2*UStop21**2 + 
     &      Sf4820**3*Sf8314*MSq2Diff(1,2)**2*UStop22**2)

	Sf8256 = Sf2320**3*Sf8255 + 
     &   (Sf4820**3*Sf8253 + 
     &      Sf2320*(Sf4820**3*Sf8252 - Sf2320**2*Sf8255))*
     &    UStop22**2

	Sf8227 = Sf2334**3*Sf8226 + 
     &   (Sf4906**3*Sf8224 + 
     &      Sf2334*(Sf4906**3*Sf8223 - Sf2334**2*Sf8226))*
     &    UStop22**2

	Sf2311 = (Sf2309*TB + MTop2*USbot21)*UStop11 + 
     &   MTop*Sf2310*UStop12

	Sf2260 = (Sf2258*TB + MTop2*USbot11)*UStop11 + 
     &   MTop*Sf2259*UStop12

	Sf1383 = 9*ML2*Sf1380*UStau21 - Sf1382*UStau21 - 
     &   9*ML*Sf1380*Sf215*UStau22

	Sf1376 = 27*ML2*Sf1373*UStau21 - Sf1375*UStau21 - 
     &   27*ML*Sf1373*Sf215*UStau22

	Sf1761 = 2*SW2*UStau12*(Sf1760 - MZ2*UStau12**2) + 
     &   (Sf1755*Ss443*UStau11*UStau21*UStau22)/MZ2

	Sf8349 = -(Sf166**3*Sf8348*USbot21**2*
     &      UStop21**2) + 
     &   Sf4906*(Sf4906**2*Sf8346*UStau22**2 - 
     &      Sf166**3*Sf8347*USbot21**2*UStop21**2)

	Sf8605 = MZ2**4*(-(MZ2*Sf8603) - 3*Sf8604) + Sf8602*UStau11

	Sf466 = 9*ML2*Sf463*UStau11 - Sf465*UStau11 - 
     &   9*ML*Sf215*Sf463*UStau12

	Sf459 = 27*ML2*Sf456*UStau11 - Sf458*UStau11 - 
     &   27*ML*Sf215*Sf456*UStau12

	Sf4779 = 2*SW2*USbot12*(Sf4778 - MZ2*USbot12**2) + 
     &   (Sf15*Sf4773*USbot11*USbot21*USbot22)/MZ2

	Sf8198 = -(Sf4820**3*Sf8197) + 
     &   (-(Sf4906**3*Sf8195) + 
     &      Sf4820*(Sf4906**3*Sf8194 + Sf4820**2*Sf8197))*
     &    USbot22**2

	Sf7989 = -(Sf2320**3*Sf7988) + 
     &   (-(Sf2334**3*Sf7986) + 
     &      Sf2320*(Sf2334**3*Sf7985 + Sf2320**2*Sf7988))*
     &    USbot22**2

	Sf7727 = MZ2*Sf7726*MStau2(1) + 
     &   CW2*(2*MZ2**2*Sf7579*Ss78*SW2*MStau2(3) - 
     &      32*Sf7579*Ss78*SW2*MStau2(1)**2*MStau2(3))

	Sf7749 = MZ2*Sf7748*MStau2(2) + 
     &   CW2*(-2*MZ2**2*Sf7713*Ss78*SW2*MStau2(3) + 
     &      32*Sf7713*Ss78*SW2*MStau2(2)**2*MStau2(3))

	Sf5050 = -3*Sf5048*MSq2Diff(2,2)**2 + 
     &   MW2*(Sf5049*MSq2Diff(2,2) - 4*MW2*MStop2(2))

	Sf2596 = -3*Sf2594*MSq2Diff(2,1)**2 + 
     &   MW2*(Sf2595*MSq2Diff(2,1) - 4*MW2*MStop2(1))

	Sf4924 = -3*Sf4922*MSq2Diff(1,2)**2 + 
     &   MW2*(Sf4923*MSq2Diff(1,2) - 4*MW2*MStop2(2))

	Sf2358 = -3*Sf2355*MSq2Diff(1,1)**2 + 
     &   MW2*(Sf2357*MSq2Diff(1,1) - 4*MW2*MStop2(1))

	Sf1979 = -2*Sf1978*MStau2(2) + 
     &   MSneu2(1)*(25*MSneu2(1) + 3*MStau2(2))

	Sf1919 = -2*Sf1918*MStau2(1) + 
     &   MSneu2(1)*(25*MSneu2(1) + 3*MStau2(1))

	Sf7316 = -(MSl2Diff(1,1)*MStau2(2)) - 
     &   MSneu2(1)*MStau2(3)*UStau22**2

	Sf6038 = -14*MStop2(2)**3 + 
     &   MSbot2(2)*(7*MSbot2(2)**2 + Sf6037*MStop2(2))

	Sf3698 = -14*MStop2(1)**3 + 
     &   MSbot2(2)*(7*MSbot2(2)**2 + Sf3697*MStop2(1))

	Sf6091 = 45*Sf6090*MSbot2(2)**4 + 
     &   MW2*(MW2*Sf2731*Sf6089*MSbot2(2) - 8*Sf6088*MSbot2(2)**3)

	Sf3420 = -14*MStop2(1)**3 + 
     &   MSbot2(1)*(7*MSbot2(1)**2 + Sf3419*MStop2(1))

	Sf5702 = 45*Sf5701*MSbot2(1)**4 + 
     &   MW2*(MW2*Sf2553*Sf5699*MSbot2(1) - 8*Sf5698*MSbot2(1)**3)

	Sf9758 = Sf9757*MStau2(2)**2 + 
     &   MZ2*MHiggs2(2)*(MZ2*Sf9507*MHiggs2(2) + Sf9508*MStau2(2))

	Sf9514 = Sf9513*MStau2(1)**2 + 
     &   MZ2*MHiggs2(2)*(MZ2*Sf9507*MHiggs2(2) + Sf9508*MStau2(1))

	Sf9729 = -(Sf9728*MStau2(2)**2) + 
     &   MZ2*MHiggs2(2)*(2*MZ2*Sf9449*MHiggs2(2) - 
     &      Sf9450*MStau2(2))

	Sf9456 = -(Sf9455*MStau2(1)**2) + 
     &   MZ2*MHiggs2(2)*(2*MZ2*Sf9449*MHiggs2(2) - 
     &      Sf9450*MStau2(1))

	Sf9073 = MZ2**2*(MZ2*Sf9068 - Sf9070) - Sf9072*MHiggs2(2)**2

	Sf9064 = MZ2**2*(MZ2*Sf9059 - Sf9061) - Sf9063*MHiggs2(2)**2

	Sf8481 = -(Sf8479*USbot21**2*UStop21**2) + 
     &   Sf4906*(6*Sf4906**2*Sf8480 - 
     &      Sf8478*USbot21**2*UStop21**2)

	Sf7433 = MW2*(16*MW2*Sf7426 - 8*MZ2*Sf7429) + MZ2**2*Sf7432

	Sf7417 = MW2*(16*MW2*Sf7409 - 8*MZ2*Sf7413) + MZ2**2*Sf7416

	Sf7385 = MZ2**2*Sf7381 + MW2*(4*MZ2*Sf7378 + 4*MW2*Sf7384)

	Sf7365 = MW2*(4*MZ2*Sf7358 + 4*MW2*Sf7361) + MZ2**2*Sf7364

	Sf7348 = MW2*(4*MW2*Sf7341 - 4*MZ2*Sf7344) + MZ2**2*Sf7347

	Sf7336 = MW2*(4*MW2*Sf7328 - 4*MZ2*Sf7332) + MZ2**2*Sf7335

	Sf6070 = 5*Sf6069*MSbot2(2)**2 + 
     &   MW2*(MW2*Sf6067 + Sf6068*MSbot2(2))

	Sf5857 = 60*MStop2(2)**4 + 
     &   MW2*(MW2**2*Sf5855 + 3*Sf5856*MStop2(2)**2)

	Sf5655 = 5*Sf5654*MSbot2(1)**2 + 
     &   MW2*(MW2*Sf5650 + Sf5653*MSbot2(1))

	Sf5514 = -MStop2(2)**3 + MW2*(MW2*Sf5513 + 9*MStop2(2)**2)

	Sf5208 = 4*MW2**2*Sf4870 + 
     &   MHiggs2(4)*(-(MW2*Sf5206) + 2*Sf5207*MHiggs2(4))

	Sf5196 = 4*MW2**2*Sf4835 + 
     &   MHiggs2(4)*(-(MW2*Sf5194) + 2*Sf5195*MHiggs2(4))

	Sf4079 = 60*MStop2(1)**4 + 
     &   MW2*(MW2**2*Sf4075 + 3*Sf4078*MStop2(1)**2)

	Sf6281 = Sf3942 + MStop2(2)*(-Sf3939 + 4*MStop2(2)*USbot21)

	Sf3970 = Sf3942 + MStop2(1)*(-Sf3939 + 4*MStop2(1)*USbot21)

	Sf6267 = Sf3819 + MStop2(2)*(-Sf3815 + 4*MStop2(2)*USbot11)

	Sf3873 = Sf3819 + MStop2(1)*(-Sf3815 + 4*MStop2(1)*USbot11)

	Sf3951 = 5*Sf3950*MSbot2(2)**2 + 
     &   MW2*(70*MW2*Sf3618 + 2*Sf3949*MSbot2(2))

	Sf3622 = Sf3621*MSbot2(2)**2 + 
     &   MW2*(70*MW2*Sf3618 + 10*Sf3619*MSbot2(2))

	Sf3836 = 5*Sf3835*MSbot2(1)**2 + 
     &   MW2*(70*MW2*Sf3294 + 2*Sf3834*MSbot2(1))

	Sf3302 = Sf3301*MSbot2(1)**2 + 
     &   MW2*(70*MW2*Sf3294 + 10*Sf3297*MSbot2(1))

	Sf3136 = -MStop2(1)**3 + MW2*(MW2*Sf3135 + 9*MStop2(1)**2)

	Sf4815 = MZ2*(Sf21*Sf23*Sf25**2 + MZ2*Sf26*Sf30**2) + 
     &   Sf21*Sf23*Sf25**2*MSbot2(3)

	Sf3748 = 4*MW2**2*Sf2273 + 
     &   MHiggs2(4)*(-(MW2*Sf3746) + 2*Sf3747*MHiggs2(4))

	Sf3739 = 4*MW2**2*Sf2210 + 
     &   MHiggs2(4)*(-(MW2*Sf3737) + 2*Sf3738*MHiggs2(4))

	Sf201 = 2*MW2**2*Sf196 + 
     &   MHiggs2(4)*(-(MW2*Sf200) + Sf198*MHiggs2(4))

	Sf1884 = Sf1883 + 
     &   MSneu2(1)*(Sf1880 - 5*MHiggs2(4)*MSneu2(1)*UStau21)

	Sf188 = -4*MW2**2*Sf185 + 
     &   MHiggs2(4)*(13*MW2*Sf176 - Sf187*MHiggs2(4))

	Sf4810 = MZ2*(MZ2*Sf11*Sf16**2 + Sf21*Sf23*Sf25**2) - 
     &   Sf21*Sf23*Sf25**2*MSbot2(3)

	Sf1204 = -6*Sf202*Ss241*Ss242*MStau2(2)**2 + 
     &   MW2*(MW2*Sf1203 + Sf201*MStau2(2))

	Sf1143 = Sf1142 + 
     &   MSneu2(1)*(Sf1136 - 5*MHiggs2(4)*MSneu2(1)*UStau11)

	Sf7230 = MZ2*(Sf108*Sf110*Sf112**2 + MZ2*Sf82*Sf87**2) - 
     &   Sf108*Sf110*Sf112**2*MStop2(3)

	Sf7235 = MZ2*(Sf108*Sf110*Sf112**2 + 
     &      MZ2*Sf113*Sf117**2) + Sf108*Sf110*Sf112**2*MStop2(3)

	Sf10566 = MZ2**2*(MZ2*Sf10561 - Sf10563) - Sf10565*MHiggs2(2)**2

	Sf10558 = MZ2**2*(MZ2*Sf10553 - Sf10555) - Sf10557*MHiggs2(2)**2

	Sf10319 = MZ2**2*(MZ2*Sf10314 - Sf10316) - Sf10318*MHiggs2(1)**2

	Sf10310 = MZ2**2*(MZ2*Sf10305 - Sf10307) - Sf10309*MHiggs2(1)**2

	Sf9761 = Sf9760*MStau2(2)**2 + 
     &   MZ2*MHiggs2(2)*(8*MZ2**2*MHiggs2(2) + Sf9493*MStau2(2))

	Sf9499 = Sf9498*MStau2(1)**2 + 
     &   MZ2*MHiggs2(2)*(8*MZ2**2*MHiggs2(2) + Sf9493*MStau2(1))

	tmp30 = MZ2*(-6*MZ2 + 10*MHiggs2(2))

	Sf9780 = tmp30*MHiggs2(2) + Sf9779*MStau2(2)

	Sf9547 = tmp30*MHiggs2(2) + Sf9546*MStau2(1)

	Sf6308 = 20*MSbot2(1)**3 + MW2*(14*MW2**2 + Sf6306*MSbot2(1))

	Sf6307 = -5*MSbot2(1)**3 + MW2*(3*MW2**2 - Sf6306*MSbot2(1))

	Sf1030 = 3*MW2**3 + MHiggs2(4)*(MW2*Sf1029 - Sf439*MHiggs2(4))

	Sf1162 = MW2*(16*MW2 - MZ2*Sf1159) + MZ2**2*Sf1161

	Sf1027 = MW2*Sf1026 + MHiggs2(4)*(-Sf1025 + 11*MHiggs2(4))

	Sf443 = 4*MW2**2 + MHiggs2(4)*(-Sf442 + 6*MHiggs2(4))

	Sf1019 = -(MW2*Sf1018) + MHiggs2(4)*(Sf1017 + 5*MHiggs2(4))

	Sf236 = 2*MW2**2 + MHiggs2(4)*(-Sf235 + 4*MHiggs2(4))

	Sf1021 = -(MW2*Sf170) + MHiggs2(4)*(Sf963 + 2*MHiggs2(4))

	Sf1015 = MW2*Sf3 + MHiggs2(4)*(-Sf1014 + 2*MHiggs2(4))

	Sf440 = 2*MW2**2 + MHiggs2(4)*(Sf439 + 2*MHiggs2(4))

	Sf447 = MW2**2 - Sf3*MHiggs2(4) + MHiggs2(4)**2

	Sf245 = MW2**2 + Sf244*MHiggs2(4) + MHiggs2(4)**2

	Sf10545 = MHiggs2(2)*(-4*MZ2 + 5*MHiggs2(2)) - Sf10544*MStop2(2)

	Sf10528 = MHiggs2(2)*(-4*MZ2 + 5*MHiggs2(2)) - Sf10527*MStop2(1)

	Sf9621 = MHiggs2(2)*(-4*MZ2 + 5*MHiggs2(2)) - Sf9620*MStau2(2)

	Sf9602 = MHiggs2(2)*(-4*MZ2 + 5*MHiggs2(2)) - Sf9601*MStau2(1)

	Sf9047 = MHiggs2(2)*(-4*MZ2 + 5*MHiggs2(2)) - Sf9046*MSbot2(2)

	Sf9028 = MHiggs2(2)*(-4*MZ2 + 5*MHiggs2(2)) - Sf9027*MSbot2(1)

	Sf10547 = MHiggs2(2)*(-2*MZ2 + 2*MHiggs2(2)) + Sf10546*MStop2(2)

	Sf10526 = MHiggs2(2)*(-2*MZ2 + 2*MHiggs2(2)) + Sf10525*MStop2(1)

	Sf9619 = MHiggs2(2)*(-2*MZ2 + 2*MHiggs2(2)) + Sf9618*MStau2(2)

	Sf9600 = MHiggs2(2)*(-2*MZ2 + 2*MHiggs2(2)) + Sf9599*MStau2(1)

	Sf9045 = MHiggs2(2)*(-2*MZ2 + 2*MHiggs2(2)) + Sf9044*MSbot2(2)

	Sf9026 = MHiggs2(2)*(-2*MZ2 + 2*MHiggs2(2)) + Sf9025*MSbot2(1)

	Sf10405 = MHiggs2(1)*(-4*MZ2 + 5*MHiggs2(1)) - Sf10404*MStop2(2)

	Sf10393 = MHiggs2(1)*(-4*MZ2 + 5*MHiggs2(1)) - Sf10392*MStop2(1)

	Sf10184 = MHiggs2(1)*(-4*MZ2 + 5*MHiggs2(1)) - Sf10183*MStau2(2)

	Sf10167 = MHiggs2(1)*(-4*MZ2 + 5*MHiggs2(1)) - Sf10166*MStau2(1)

	Sf9250 = MHiggs2(1)*(-4*MZ2 + 5*MHiggs2(1)) - Sf9249*MSbot2(2)

	Sf9233 = MHiggs2(1)*(-4*MZ2 + 5*MHiggs2(1)) - Sf9232*MSbot2(1)

	Sf10403 = MHiggs2(1)*(-2*MZ2 + 2*MHiggs2(1)) + Sf10402*MStop2(2)

	Sf10391 = MHiggs2(1)*(-2*MZ2 + 2*MHiggs2(1)) + Sf10390*MStop2(1)

	Sf10186 = MHiggs2(1)*(-2*MZ2 + 2*MHiggs2(1)) + Sf10185*MStau2(2)

	Sf10165 = MHiggs2(1)*(-2*MZ2 + 2*MHiggs2(1)) + Sf10164*MStau2(1)

	Sf9252 = MHiggs2(1)*(-2*MZ2 + 2*MHiggs2(1)) + Sf9251*MSbot2(2)

	Sf9231 = MHiggs2(1)*(-2*MZ2 + 2*MHiggs2(1)) + Sf9230*MSbot2(1)

	Sf6886 = Sf3578*TB2*UStop21**2 + 
     &   MTop2*Sf2245**2*USbot21**2*UStop22**2

	Sf6635 = Sf3226*TB2*UStop21**2 + 
     &   MTop2*Sf2245**2*USbot22**2*UStop22**2

	Sf4899 = Sf2317*TB*UStop21**2 + Sf2318*UStop22**2

	Sf7224 = Sf7219*Sf86*UStop11*UStop12*UStop21 + 
     &   4*MZ2*Sf7223*SW2*UStop22

	Sf6187 = MTop*Sf6185*TB*UStop21 + 
     &   3*MTop2*MW2**2*Sf2245*Sf6186*MStop2(2)*USbot21*UStop22

	Sf6182 = MTop*Sf6181*TB*UStop21 - 
     &   3*MTop2*MW2**2*Sf2245*Sf96**2*MSq2Diff(2,2)*USbot21*
     &    UStop22

	Sf5809 = MTop*Sf5807*TB*UStop21 + 
     &   3*MTop2*MW2**2*Sf2245*Sf5808*MStop2(2)*USbot11*UStop22

	Sf5804 = MTop*Sf5803*TB*UStop21 - 
     &   3*MTop2*MW2**2*Sf2245*Sf88**2*MSq2Diff(1,2)*USbot11*
     &    UStop22

	Sf5144 = Sf5143*USbot21 + 4*Sf2731*Sf5140*TB*UStop21

	Sf5025 = Sf5024*USbot11 + 4*Sf2553*Sf5021*TB*UStop21

	Sf5442 = 2*MTop*Sf3051*TB*UStop21 - 
     &   MTop2*Sf2245*USbot21*UStop22

	Sf5371 = 2*MBot*MTop*Sf2962*TB*UStop21 + 
     &   MTop2*Sf2245*USbot21*UStop22

	Sf5336 = 2*MTop*Sf2924*TB*UStop21 - 
     &   MTop2*Sf2245*USbot11*UStop22

	Sf5263 = 2*MBot*MTop*Sf2810*TB*UStop21 + 
     &   MTop2*Sf2245*USbot11*UStop22

	Sf5148 = Sf2731*TB*UStop21 + MTop*Sf2245*USbot21*UStop22

	Sf4827 = Sf2553*TB*UStop21 + MTop*Sf2245*USbot11*UStop22

	Sf4864 = Sf2261*UStop21 + MTop*Sf2262*UStop22

	Sf4828 = Sf2195*UStop21 + MTop*Sf2197*UStop22

	Sf6976 = 3*MTop*reMUE*SA*SB*UStop21 + 2*Sf4258*UStop22

	Sf4681 = MTop2*Sf2245**2*USbot21**2*UStop12**2 + 
     &   Sf3578*TB2*UStop22**2

	Sf4576 = MTop2*Sf2245**2*USbot22**2*UStop12**2 + 
     &   Sf3226*TB2*UStop22**2

	Sf2319 = Sf2318*UStop12**2 + Sf2317*TB*UStop22**2

	Sf6105 = MTop*Sf6103*TB*UStop11 + 
     &   3*MTop2*MW2**2*Sf2245*Sf6104*MStop2(1)*USbot21*UStop12

	Sf5716 = MTop*Sf5714*TB*UStop11 + 
     &   3*MTop2*MW2**2*Sf2245*Sf5715*MStop2(1)*USbot11*UStop12

	Sf3937 = MTop*Sf3936*TB*UStop11 - 
     &   3*MTop2*MW2**2*Sf2245*Sf74**2*MSq2Diff(2,1)*USbot21*
     &    UStop12

	Sf3813 = MTop*Sf3812*TB*UStop11 - 
     &   3*MTop2*MW2**2*Sf2245*Sf66**2*MSq2Diff(1,1)*USbot11*
     &    UStop12

	Sf3052 = 2*MTop*Sf3051*TB*UStop11 - 
     &   MTop2*Sf2245*USbot21*UStop12

	Sf2963 = 2*MBot*MTop*Sf2962*TB*UStop11 + 
     &   MTop2*Sf2245*USbot21*UStop12

	Sf2925 = 2*MTop*Sf2924*TB*UStop11 - 
     &   MTop2*Sf2245*USbot11*UStop12

	Sf2811 = 2*MBot*MTop*Sf2810*TB*UStop11 + 
     &   MTop2*Sf2245*USbot11*UStop12

	Sf2736 = Sf2735*USbot21 + 4*Sf2731*Sf2732*TB*UStop11

	Sf2740 = Sf2731*TB*UStop11 + MTop*Sf2245*USbot21*UStop12

	Sf2558 = Sf2557*USbot11 + 4*Sf2553*Sf2554*TB*UStop11

	Sf2562 = Sf2553*TB*UStop11 + MTop*Sf2245*USbot11*UStop12

	Sf7221 = (Sf112*Sf7219*UStop11 - MZ2*Sf7220*UStop21)/MZ2

	Sf2263 = Sf2261*UStop11 + MTop*Sf2262*UStop12

	Sf2198 = Sf2195*UStop11 + MTop*Sf2197*UStop12

	Sf4259 = 3*MTop*reMUE*SA*SB*UStop11 + 2*Sf4258*UStop12

	Sf1298 = Sf345*UStau21**2 - Sf1297*Sf215*UStau22

	Sf1781 = Sf1775*Ss443*UStau11*UStau12*UStau21 + 
     &   2*MZ2*Sf1780*SW2*UStau22

	Sf1365 = Sf447*UStau21 + ML*Sf215*MHiggs2(4)*UStau22

	Sf1228 = Sf245*UStau21 - ML*Sf215*MHiggs2(4)*UStau22

	Sf2003 = Sf1952*UStau21 - ML*Sf1953*Sf215*UStau22

	Sf2002 = Sf1947*UStau21 + ML*Sf1948*Sf215*UStau22

	Sf1737 = Sf1736*UStau21 + Sf1732*UStau22

	Sf1456 = Sf1408*UStau21 - 3*ML*Sf1412*Sf215*UStau22

	Sf1366 = Sf1364*UStau21 + 2*ML*Sf215*MHiggs2(4)*UStau22

	Sf1273 = ML*Sf1271*UStau21 - Sf1272*Sf215*UStau22

	Sf1869 = Sf1117*UStau21 - 3*ML*Sf215*Ss115*UStau22

	Sf1868 = Sf1113*UStau21 + ML*Sf1114*Sf215*UStau22

	Sf1865 = Sf1102*UStau21 + ML*Sf1103*Sf215*UStau22

	Sf1864 = Sf1097*UStau21 + ML*Sf215*Sf733*UStau22

	Sf1825 = Sf1030*UStau21 + ML*Sf215*Ss25*MHiggs2(4)*UStau22

	Sf1824 = Sf1027*UStau21 + 3*ML*Sf215*Ss241*Ss242*UStau22

	Sf1823 = Sf1021*UStau21 - 3*ML*Sf215*Ss241*Ss242*UStau22

	Sf1822 = Sf1019*UStau21 - ML*Sf215*Ss241*Ss242*UStau22

	Sf1821 = Sf1015*UStau21 + ML*Sf215*Ss241*Ss242*UStau22

	Sf1820 = Sf1012*UStau21 - ML*Sf215*Sf217*UStau22

	Sf8635 = -4*Sf8548*Sf8611*UStau12 + 
     &   Sf8534*Sf8597*Ss443**2*UStau22**4

	Sf7670 = Sf166**3*Sf7669*UStau21**2 + 
     &   Sf1195**3*Sf7668*UStau22**2

	Sf347 = -(Sf215*Sf346*UStau12) + Sf345*UStau22**2

	Sf949 = Sf948*UStau11 + Sf944*UStau12

	Sf8582 = -3*MZ2**3*Sf8551 + Sf8581*UStau11

	Sf8552 = -21*MZ2**3*Sf8551 + Sf8546*UStau11

	Sf7763 = MW2**2*Sf1775*Sf7762*UStau11 - MZ2*Sf7761*UStau21

	Sf1778 = (Sf1775*Sf60*UStau11 - MZ2*Sf1777*UStau21)/MZ2

	Sf553 = Sf491*UStau11 - 3*ML*Sf215*Sf495*UStau12

	Sf448 = Sf447*UStau11 + ML*Sf215*MHiggs2(4)*UStau12

	Sf449 = Sf446*UStau11 + 2*ML*Sf215*MHiggs2(4)*UStau12

	Sf313 = ML*Sf311*UStau11 - Sf215*Sf312*UStau12

	Sf246 = Sf245*UStau11 - ML*Sf215*MHiggs2(4)*UStau12

	Sf1954 = Sf1952*UStau11 - ML*Sf1953*Sf215*UStau12

	Sf1949 = Sf1947*UStau11 + ML*Sf1948*Sf215*UStau12

	end

