* gm2_2L_1.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:39
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_1
	implicit none

#include "gm2_2L.h"

	Ss55 = SA*(1 - 4*SW2)

	Ss133 = MZ2*(1 + 3*SA2)

	Ss53 = MM*(1 - 4*SW2)

	Ss824 = ML*(1 - 4*SW2)

	Ss51 = CA*(1 - 4*SW2)

	Ss595 = -4313*MW2**6 + 
     &   MTop2*(9217*MW2**5 + 
     &      MTop2*(-7902*MW2**4 + 
     &         MTop2*(-12502*MW2**3 + 
     &            MTop2*(MTop2*(7456*MTop2 - 39579*MW2) + 39639*MW2**2)))
     &      )

	Ss586 = -5741*MW2**6 + 
     &   MTop2*(-2921*MW2**5 + 
     &      MTop2*(-5192*MW2**4 + 
     &         MTop2*(8458*MW2**3 + 
     &            MTop2*(MTop2*(5488*MTop2 - 38853*MW2) + 40089*MW2**2)))
     &      )

	Ss604 = 19*MW2**6 + 
     &   MTop2*(-173*MW2**5 + 
     &      MTop2*(-1004*MW2**4 + 
     &         MTop2*(-2522*MW2**3 + 
     &            MTop2*(MTop2*(2128*MTop2 - 12573*MW2) + 11805*MW2**2)))
     &      )

	Ss772 = 2*MW2**6 + 
     &   MTop2*(33*MW2**5 + MTop2*
     &       (2*MW2**4 + MTop2*
     &          (-170*MW2**3 + 
     &            MTop2*(MTop2*(532*MTop2 - 1083*MW2) + 700*MW2**2))))

	Ss594 = -11565*MW2**5 + 
     &   MTop2*(6232*MW2**4 + 
     &      MTop2*(-8084*MW2**3 + 
     &         MTop2*(MTop2*(16904*MTop2 - 45583*MW2) + 852*MW2**2)))

	Ss585 = -7921*MW2**5 + 
     &   MTop2*(-19202*MW2**4 + 
     &      MTop2*(-11110*MW2**3 + 
     &         MTop2*(MTop2*(10976*MTop2 - 44755*MW2) - 13504*MW2**2)))

	Ss602 = 899*MW2**5 + 
     &   MTop2*(-74*MW2**4 + 
     &      MTop2*(-874*MW2**3 + 
     &         MTop2*(MTop2*(4256*MTop2 - 10555*MW2) - 640*MW2**2)))

	Ss593 = -652*MW2**5 + 
     &   MTop2*(1463*MW2**4 + 
     &      MTop2*(781*MW2**3 + 
     &         MTop2*(MTop2*(1656*MTop2 - 11293*MW2) + 14017*MW2**2)))

	Ss584 = -1486*MW2**5 + 
     &   MTop2*(-3334*MW2**4 + 
     &      MTop2*(10945*MW2**3 + 
     &         MTop2*(MTop2*(1568*MTop2 - 13815*MW2) + 19894*MW2**2)))

	Ss395 = 195*MW2**5 + 
     &   MTop2*(-1522*MW2**4 + 
     &      MTop2*(-1938*MW2**3 + 
     &         MTop2*(MTop2*(1452*MTop2 - 7263*MW2) + 924*MW2**2)))

	Ss406 = -1029*MW2**5 + 
     &   MTop2*(758*MW2**4 + 
     &      MTop2*(1728*MW2**3 + 
     &         MTop2*(MTop2*(1236*MTop2 - 5667*MW2) + 1038*MW2**2)))

	Ss771 = 29*MW2**5 + 
     &   MTop2*(70*MW2**4 + MTop2*
     &       (42*MW2**3 + MTop2*
     &          (MTop2*(1064*MTop2 - 765*MW2) + 200*MW2**2)))

	Ss603 = 8*MW2**5 + 
     &   MTop2*(-658*MW2**4 + 
     &      MTop2*(-161*MW2**3 + 
     &         MTop2*(MTop2*(608*MTop2 - 5055*MW2) + 4342*MW2**2)))

	Ss396 = 77*MW2**5 + 
     &   MTop2*(133*MW2**4 + 
     &      MTop2*(215*MW2**3 + 
     &         MTop2*(MTop2*(508*MTop2 - 4238*MW2) + 3121*MW2**2)))

	Ss752 = 29*MW2**5 + 
     &   MTop2*(-5*MW2**4 + MTop2*
     &       (-22*MW2**3 + 
     &         MTop2*(MTop2*(399*MTop2 - 435*MW2) + 122*MW2**2)))

	Ss762 = 42*MW2**5 + 
     &   MTop2*(-25*MW2**4 + 
     &      MTop2*(-46*MW2**3 + 
     &         MTop2*(MTop2*(399*MTop2 - 720*MW2) + 366*MW2**2)))

	Ss753 = -2*MW2**5 + 
     &   MTop2*(7*MW2**4 + MTop2*
     &       (-20*MW2**3 + 
     &         MTop2*(MTop2*(133*MTop2 - 230*MW2) + 104*MW2**2)))

	Ss592 = -16046*MW2**4 + 
     &   MTop2*(-5458*MW2**3 + 
     &      MTop2*(MTop2*(21960*MTop2 - 20927*MW2) - 10785*MW2**2))

	Ss582 = -3728*MW2**4 + 
     &   MTop2*(-21506*MW2**3 + 
     &      MTop2*(MTop2*(13720*MTop2 - 24755*MW2) - 27921*MW2**2))

	Ss333 = 3015*MW2**4 + 
     &   MTop2*(-557*MW2**3 + 
     &      MTop2*(MTop2*(5700*MTop2 - 3594*MW2) + 864*MW2**2))

	Ss601 = 3082*MW2**4 + 
     &   MTop2*(2548*MW2**3 + 
     &      MTop2*(MTop2*(5320*MTop2 - 155*MW2) + 795*MW2**2))

	Ss332 = -1125*MW2**4 + 
     &   MTop2*(1424*MW2**3 + 
     &      MTop2*(MTop2*(2280*MTop2 - 2961*MW2) + 474*MW2**2))

	Ss770 = 115*MW2**4 + 
     &   MTop2*(248*MW2**3 + 
     &      MTop2*(350*MW2**2 + MTop2*(1330*MTop2 + 275*MW2)))

	Ss394 = 641*MW2**4 + 
     &   MTop2*(-691*MW2**3 + 
     &      MTop2*(MTop2*(1150*MTop2 - 1966*MW2) - 1287*MW2**2))

	Ss751 = 115*MW2**4 + 
     &   MTop2*(136*MW2**3 + 
     &      MTop2*(162*MW2**2 + MTop2*(665*MTop2 + 50*MW2)))

	Ss761 = 75*MW2**4 + 
     &   MTop2*(64*MW2**3 + MTop2*
     &       (MTop2*(665*MTop2 - 330*MW2) + 78*MW2**2))

	Ss405 = -70*MW2**4 + 
     &   MTop2*(-406*MW2**3 + 
     &      MTop2*(MTop2*(545*MTop2 - 547*MW2) - 927*MW2**2))

	Ss741 = 95*MW2**4 + 
     &   MTop2*(-61*MW2**3 + 
     &      MTop2*(MTop2*(285*MTop2 - 259*MW2) - 8*MW2**2))

	Ss740 = -71*MW2**4 + 
     &   MTop2*(149*MW2**3 + 
     &      MTop2*(MTop2*(228*MTop2 - 409*MW2) + 95*MW2**2))

	Ss583 = -2*MW2**4 + 
     &   MTop2*(-1942*MW2**3 + 
     &      MTop2*(MTop2*(196*MTop2 - 2209*MW2) + 5841*MW2**2))

	Ss404 = -332*MW2**4 + 
     &   MTop2*(525*MW2**3 + 
     &      MTop2*(MTop2*(194*MTop2 - 1493*MW2) + 552*MW2**2))

	Ss769 = 10*MW2**4 + 
     &   MTop2*(-23*MW2**3 + 
     &      MTop2*(MTop2*(152*MTop2 - 265*MW2) + 118*MW2**2))

	Ss760 = 12*MW2**4 + 
     &   MTop2*(-7*MW2**3 + MTop2*
     &       (MTop2*(133*MTop2 - 211*MW2) + 65*MW2**2))

	Ss393 = 16*MW2**4 + 
     &   MTop2*(-38*MW2**3 + 
     &      MTop2*(MTop2*(76*MTop2 - 889*MW2) + 1101*MW2**2))

	Ss430 = 16*MW2**4 + 
     &   ML2*(56*MW2**3 + ML2*
     &       (ML2*(76*ML2 - 535*MW2) + 291*MW2**2))

	Ss609 = 255*MTop2**4 + 
     &   MBot2*(-174*MTop2**3 + 
     &      MBot2*(MBot2*(115*MBot2 - 70*MTop2) + 2*MTop2**2))

	Ss129 = -41*MZ2**3 + MW2*(MW2*(192*MW2 - 288*MZ2) + 137*MZ2**2)

	Ss692 = 6*ML2*MZ2**2 + MW2*(MW2*(96*MW2 - 48*MZ2) - 5*MZ2**2)

	Ss683 = 5*MZ2**3 + MW2*(MW2*(96*MW2 - 144*MZ2) + 43*MZ2**2)

	Ss591 = -11974*MW2**3 + 
     &   MTop2*(MTop2*(16904*MTop2 - 627*MW2) - 6087*MW2**2)

	Ss581 = 1376*MW2**3 + 
     &   MTop2*(MTop2*(10976*MTop2 - 2853*MW2) - 5316*MW2**2)

	Ss331 = 4220*MW2**3 + 
     &   MTop2*(4101*MW2**2 + MTop2*(7600*MTop2 + 3294*MW2))

	Ss600 = 4334*MW2**3 + 
     &   MTop2*(6252*MW2**2 + MTop2*(4256*MTop2 + 6147*MW2))

	Ss403 = 1410*MW2**3 + 
     &   MTop2*(921*MW2**2 + MTop2*(2300*MTop2 + 2673*MW2))

	Ss392 = 2198*MW2**3 + 
     &   MTop2*(1536*MW2**2 + MTop2*(2180*MTop2 + 1077*MW2))

	Ss768 = 210*MW2**3 + MTop2*(520*MW2**2 + MTop2*(1064*MTop2 + 789*MW2))

	Ss330 = -309*MW2**3 + 
     &   MTop2*(MTop2*(760*MTop2 - 1194*MW2) + 477*MW2**2)

	Ss750 = 210*MW2**3 + MTop2*(394*MW2**2 + MTop2*(665*MTop2 + 510*MW2))

	Ss759 = 112*MW2**3 + MTop2*(188*MW2**2 + MTop2*(665*MTop2 + 225*MW2))

	Ss739 = 161*MW2**3 + MTop2*(128*MW2**2 + MTop2*(380*MTop2 + 119*MW2))

	Ss590 = -194*MW2**3 + 
     &   MTop2*(MTop2*(124*MTop2 - 1153*MW2) + 3107*MW2**2)

	Ss599 = -16*MW2**3 + MTop2*(MTop2*(76*MTop2 - 813*MW2) + 312*MW2**2)

	Ss738 = -14*MW2**3 + MTop2*(MTop2*(38*MTop2 - 89*MW2) + 69*MW2**2)

	Ss749 = 2*MW2**3 + MTop2*(MTop2*(19*MTop2 - 35*MW2) + 18*MW2**2)

	Ss308 = 16*MW2**3 + ML2*(ML2*(52*ML2 - 295*MW2) + 131*MW2**2)

	Ss728 = -2*MW2**3 + ML2*(ML2*(3*ML2 - 27*MW2) - 7*MW2**2)

	Ss385 = 263*MTop2**3 + MBot2*(27*MTop2**2 + MBot2*(161*MBot2 + 61*MTop2))

	Ss459 = -83*MTop2**3 + MBot2*(27*MTop2**2 + MBot2*(115*MBot2 + 5*MTop2))

	Ss608 = 175*MTop2**3 + MBot2*(88*MTop2**2 + MBot2*(105*MBot2 + 86*MTop2))

	Ss458 = 88*MTop2**3 + MBot2*(26*MTop2**2 + MBot2*(105*MBot2 + 37*MTop2))

	Ss467 = -255*MTop2**3 + MBot2*(91*MTop2**2 + MBot2*(75*MBot2 + 25*MTop2))

	Ss466 = 175*MTop2**3 + MBot2*(MTop2**2 + MBot2*(56*MBot2 + 24*MTop2))

	Ss610 = 99*MTop2**3 + MBot2*(MBot2*(29*MBot2 - 33*MTop2) - 43*MTop2**2)

	Ss160 = -8*MW2**3 + MBot2*(MBot2*(9*MBot2 - 36*MW2) + 35*MW2**2)

	Ss127 = MW2*(64*MW2 - 80*MZ2) + 7*MZ2**2

	Ss497 = MW2*(16*MW2 - 8*MZ2) - 17*MZ2**2

	Ss502 = MW2*(16*MW2 - 24*MZ2) + 7*MZ2**2

	Ss696 = MW2*(8*MW2 - 45*MZ2) + 10*MZ2**2

	Ss702 = -MW2**3 + MTop2*(2*MTop2**2 - MW2**2)

	Ss589 = -4293*MW2**2 + MTop2*(7456*MTop2 + 1003*MW2)

	Ss580 = 2251*MW2**2 + MTop2*(5488*MTop2 + 3049*MW2)

	Ss329 = 2335*MW2**2 + MTop2*(3800*MTop2 + 3394*MW2)

	Ss598 = 2935*MW2**2 + MTop2*(2128*MTop2 + 4129*MW2)

	Ss737 = 326*MW2**2 + MTop2*(570*MTop2 + 497*MW2)

	Ss767 = 200*MW2**2 + MTop2*(532*MTop2 + 471*MW2)

	Ss391 = 557*MW2**2 + MTop2*(412*MTop2 + 642*MW2)

	Ss748 = 200*MW2**2 + MTop2*(399*MTop2 + 393*MW2)

	Ss402 = 314*MW2**2 + MTop2*(242*MTop2 + 485*MW2)

	Ss758 = 42*MW2**2 + MTop2*(133*MTop2 + 93*MW2)

	Ss42 = MTop2*(81*MTop2 - 183*MW2) + 56*MW2**2

	Ss401 = MTop2*(52*MTop2 - 597*MW2) + 176*MW2**2

	Ss39 = 28*MW2**2 + MTop2*(27*MTop2 + 102*MW2)

	Ss41 = MTop2*(27*MTop2 - 78*MW2) + 28*MW2**2

	Ss277 = 4*MW2**2 + MTop2*(9*MTop2 + 18*MW2)

	Ss351 = -8*MW2**2 + MTop2*(9*MTop2 + 9*MW2)

	Ss246 = MW2**2 + MTop2*(6*MTop2 + 4*MW2)

	Ss259 = -MW2**2 + MTop2*(3*MTop2 + 2*MW2)

	Ss258 = MTop2*(3*MTop2 - 2*MW2) + MW2**2

	Ss147 = -(MBot2**2*TB2) + MTop2*(MTop2 + 2*MBot2*TB2)

	Ss727 = ML2*(76*ML2 - 511*MW2) + 75*MW2**2

	Ss654 = 8*MW2**2 + ML2*(ML2 + 4*MW2)

	Ss232 = 12*MTop2**2 + MHiggs2(4)*(-17*MTop2 + 7*MHiggs2(4))

	Ss297 = 4*MW2**2 + MHiggs2(4)*(-8*MW2 + 5*MHiggs2(4))

	Ss376 = 5*MW2**2 + MHiggs2(4)*(-16*MW2 + 5*MHiggs2(4))

	Ss231 = -12*MTop2**2 + MHiggs2(4)*(-18*MTop2 + 5*MHiggs2(4))

	Ss377 = 4*MW2**2 + MHiggs2(4)*(-7*MW2 + 4*MHiggs2(4))

	Ss227 = 12*MTop2**2 + MHiggs2(4)*(MTop2 + 4*MHiggs2(4))

	Ss223 = 12*MTop2**2 + MHiggs2(4)*(-15*MTop2 + 4*MHiggs2(4))

	Ss245 = MW2**2 + MHiggs2(4)*(-2*MW2 + 2*MHiggs2(4))

	Ss296 = MHiggs2(4)**2 + MW2*(-5*MW2 + 10*MHiggs2(4))

	Ss255 = MW2*(MW2 - MHiggs2(4)) + MHiggs2(4)**2

	Ss409 = MW2*(-MW2 - MHiggs2(4)) + MHiggs2(4)**2

	Ss72 = MW2*(3*MW2 - 3*MHiggs2(4)) + MHiggs2(4)**2

	Ss435 = MW2*(-2*MW2 - 5*MHiggs2(4)) + MHiggs2(4)**2

	Ss222 = MTop2*(-12*MTop2 - 38*MHiggs2(4)) + MHiggs2(4)**2

	Ss170 = 8*MW2**3 + MBot2*(9*MBot2**2 - 17*MW2**2)

	Ss384 = 225*MTop2**2 + MBot2*(163*MBot2 + 198*MTop2)

	Ss457 = 90*MTop2**2 + MBot2*(100*MBot2 + 103*MTop2)

	Ss386 = 169*MTop2**2 + MBot2*(95*MBot2 + 110*MTop2)

	Ss465 = 135*MTop2**2 + MBot2*(63*MBot2 + 95*MTop2)

	Ss468 = -99*MTop2**2 + MBot2*(42*MBot2 + 31*MTop2)

	Ss460 = -12*MTop2**2 + MBot2*(29*MBot2 + 9*MTop2)

	Ss607 = 27*MTop2**2 + MBot2*(20*MBot2 + 28*MTop2)

	Ss153 = MBot2*(9*MBot2 - 36*MW2) - 26*MW2**2

	Ss166 = MBot2*(9*MBot2 - 45*MW2) + 16*MW2**2

	Ss611 = MBot2*(MBot2 - 5*MTop2) + 8*MTop2**2

	Ss711 = 4 + Log(MBot2)*(Log(MBot2) - Log(MHiggs2(1)))

	Ss56 = 3 + SW2*(-16 + 16*SW2)

	Ss63 = 3 + SW2*(-20 + 32*SW2)

	Ss859 = Log(MTop2) - Log(MZ2)

	Ss846 = Log(ML2) - Log(MZ2)

	Ss799 = Log(ML2) - 2*Log(MZ2)

	Ss532 = Log(ML2) + Log(MW2)

	Ss531 = Log(ML2) - Log(MW2)

	Ss527 = Log(ML2) + Log(MHiggs2(4))

	Ss526 = -Log(ML2) + Log(MHiggs2(4))

	Ss808 = Log(MBot2) - Log(MZ2)

	Ss513 = Log(MBot2) + Log(MTop2)

	Ss512 = Log(MBot2) - Log(MTop2)

	Ss612 = Ss610 - 2*MW2*Ss611

	Ss613 = Ss609 - MW2*Ss612

	Ss881 = reMUE*SA - CA*Re(Abot)

	Ss978 = reMUE*SA + CA*Re(Atop)

	Ss924 = reMUE*SA - CA*Re(Atau)

	Ss967 = 3*CW*MTop2*SA + 2*CAB*MW*MZ*SB*SW2

	Ss938 = CW*ML2*SA - CB*MW*MZ*SAB*SW2

	Ss895 = 3*CW*MBot2*SA - CB*MW*MZ*SAB*SW2

	Ss303 = 11*MTop2*MW2**2 + MBot2*Ss296

	Ss304 = 22*MTop2*MW2 + 3*Ss297

	Ss695 = 71*MW2 - 59*MZ2

	Ss663 = 70*MW2 - 40*MZ2

	Ss666 = 35*MW2 - 53*MZ2

	Ss693 = 32*MW2 - 427*MZ2

	Ss668 = 16*MW2 - 7*MZ2

	Ss664 = 13*MW2 - 7*MZ2

	Ss433 = 8*MW2 - 5*MZ2

	Ss670 = 7*MW2 + 11*MZ2

	Ss694 = 7*MW2 - MZ2

	Ss635 = 7*MW2 - 4*MZ2

	Ss661 = 5*MW2 + MZ2

	Ss154 = 4*MW2 + MBot2*TB2

	Ss675 = 4*MW2 + 5*MZ2

	Ss450 = 4*MW2 - 3*MZ2

	Ss161 = 2*MW2 + MBot2*TB2

	Ss171 = 2*MW2 - MBot2*TB2

	Ss677 = 2*MW2 + 7*MZ2

	Ss443 = 2*MW2 - MZ2

	Ss699 = MW2 - 4*MZ2

	Ss650 = MW + MZ

	Ss649 = MW - MZ

	Ss894 = CA*reMUE + SA*Re(Abot)

	Ss966 = CA*reMUE - SA*Re(Atop)

	Ss937 = CA*reMUE + SA*Re(Atau)

	Ss350 = 9*MTop2**2 - 4*MW2**2

	Ss701 = 6*MTop2**2 - 3*MW2**2

	Ss979 = 3*CA*CW*MTop2 - 2*MW*MZ*SAB*SB*SW2

	Ss588 = 1656*MTop2 - 337*MW2

	Ss579 = 1568*MTop2 + 1007*MW2

	Ss328 = 760*MTop2 + 569*MW2

	Ss681 = 639*MTop2 - 182*MZ2

	Ss597 = 608*MTop2 + 887*MW2

	Ss400 = 508*MTop2 + 811*MW2

	Ss720 = 392*MTop2 - 885*MW2

	Ss390 = 388*MTop2 + 551*MW2

	Ss736 = 228*MTop2 + 175*MW2

	Ss721 = 196*MTop2 - 1101*MW2

	Ss766 = 152*MTop2 + 97*MW2

	Ss724 = 152*MTop2 - 405*MW2

	Ss747 = 133*MTop2 + 97*MW2

	Ss757 = 133*MTop2 + 78*MW2

	Ss725 = 76*MTop2 + 99*MW2

	Ss281 = 63*MTop2 - 52*MW2

	Ss682 = 60*MTop2 - 59*MZ2

	Ss38 = 39*MTop2 + 28*MW2

	Ss130 = 32*MTop2 + 41*MZ2

	Ss687 = 27*MTop2 + 200*MZ2

	Ss76 = 27*MTop2 + MW2

	Ss506 = 19*MTop2 + MW2

	Ss507 = 19*MTop2 - 16*MW2

	Ss73 = 9*MTop2 + 26*MW2

	Ss276 = 9*MTop2 + 4*MW2

	Ss75 = 9*MTop2 - 4*MW2

	Ss378 = 6*MTop2 + MW2

	Ss253 = 6*MTop2 - 2*MW2

	Ss247 = 6*MTop2 - 4*MW2

	Ss9 = 4*MTop2 - MZ2

	Ss244 = 4*MTop2 + MW2

	Ss337 = 3*MTop2 - MZ2

	Ss256 = 3*MTop2 + 2*MW2

	Ss243 = 3*MTop2 + MW2

	Ss248 = 3*MTop2 - MW2

	Ss12 = 2*MTop2 + MZ2

	Ss45 = 2*MTop2 - MZ2

	Ss688 = 2*MTop2 - 3*MZ2

	Ss354 = MTop2 + MW2*TB2

	Ss353 = MTop2 - MW2*TB2

	Ss341 = MTop2 - 2*MW2*TB2

	Ss516 = MTop2 + MBot2*TB2

	Ss410 = MTop2 + 3*MW2

	Ss279 = MTop2 + 2*MW2

	Ss36 = MTop2 + MW2

	Ss35 = MTop2 - MW2

	Ss33 = MTop + MW

	Ss32 = MTop - MW

	Ss925 = CA*CW*ML2 - CAB*CB*MW*MZ*SW2

	Ss320 = 76*ML2 + 453*MW2

	Ss837 = 23*ML2 - 7*MZ2

	Ss325 = 19*ML2 - 6*MW2

	Ss8 = 4*ML2 - MZ2

	Ss730 = 4*ML2 + 21*MW2

	Ss309 = 3*ML2 - MZ2

	Ss11 = 2*ML2 + MZ2

	Ss52 = 2*ML2 - MZ2

	Ss321 = ML2*(3/2.D0 - Log(ML2))

	Ss651 = ML2 + MZ2

	Ss655 = ML2 + 2*MW2

	Ss31 = ML2 + MW2

	Ss28 = ML2 - MW2

	Ss638 = ML + MW

	Ss637 = ML - MW

	Ss340 = MW2**2 + MHiggs2(4)**2

	Ss230 = 12*MTop2 + 13*MHiggs2(4)

	Ss234 = 12*MTop2 + 5*MHiggs2(4)

	Ss226 = 6*MTop2 + 4*MHiggs2(4)

	Ss122 = 7*MW2 + 3*MHiggs2(4)

	Ss115 = MW2 + MHiggs2(4)

	Ss25 = -MW2 + MHiggs2(4)

	Ss300 = 6*MTop2 + MHiggs2(4)

	Ss221 = 3*MTop2 + MHiggs2(4)

	Ss23 = MTop2 + MHiggs2(4)

	Ss22 = -MTop2 + MHiggs2(4)

	Ss239 = -3*MTop2 + MHiggs2(4)

	Ss70 = ML2 + MHiggs2(4)

	Ss525 = -ML2 + MHiggs2(4)

	Ss242 = MW + MHiggs(4)

	Ss241 = -MW + MHiggs(4)

	Ss27 = MTop + MHiggs(4)

	Ss26 = -MTop + MHiggs(4)

	Ss991 = MHiggs2(2) + 2*MStop2(2)

	Ss986 = MHiggs2(2) - 4*MStop2(2)

	Ss972 = MHiggs2(2) + 2*MStop2(1)

	Ss964 = MHiggs2(2) - 4*MStop2(1)

	Ss952 = MHiggs2(2) + 2*MStau2(2)

	Ss946 = MHiggs2(2) - 4*MStau2(2)

	Ss931 = MHiggs2(2) + 2*MStau2(1)

	Ss921 = MHiggs2(2) - 4*MStau2(1)

	Ss909 = MHiggs2(2) + 2*MSbot2(2)

	Ss903 = MHiggs2(2) - 4*MSbot2(2)

	Ss888 = MHiggs2(2) + 2*MSbot2(1)

	Ss878 = MHiggs2(2) - 4*MSbot2(1)

	Ss143 = MZ2 + MHiggs2(2)

	Ss49 = -MZ2 + MHiggs2(2)

	Ss20 = 2*MTop2 + MHiggs2(2)

	Ss65 = -2*MTop2 + MHiggs2(2)

	Ss18 = -4*MTop2 + MHiggs2(2)

	Ss4 = 2*ML2 + MHiggs2(2)

	Ss50 = -2*ML2 + MHiggs2(2)

	Ss3 = -4*ML2 + MHiggs2(2)

	Ss854 = 4*MTop2**2 + MHiggs2(1)**2

	Ss840 = 4*ML2**2 + MHiggs2(1)**2

	Ss998 = MHiggs2(1) + 2*MStop2(2)

	Ss995 = MHiggs2(1) - 4*MStop2(2)

	Ss982 = MHiggs2(1) + 2*MStop2(1)

	Ss976 = MHiggs2(1) - 4*MStop2(1)

	Ss959 = MHiggs2(1) + 2*MStau2(2)

	Ss956 = MHiggs2(1) - 4*MStau2(2)

	Ss941 = MHiggs2(1) + 2*MStau2(1)

	Ss935 = MHiggs2(1) - 4*MStau2(1)

	Ss916 = MHiggs2(1) + 2*MSbot2(2)

	Ss913 = MHiggs2(1) - 4*MSbot2(2)

	Ss898 = MHiggs2(1) + 2*MSbot2(1)

	Ss892 = MHiggs2(1) - 4*MSbot2(1)

	Ss145 = MZ2 + MHiggs2(1)

	Ss44 = -MZ2 + MHiggs2(1)

	Ss19 = 2*MTop2 + MHiggs2(1)

	Ss64 = -2*MTop2 + MHiggs2(1)

	Ss17 = -4*MTop2 + MHiggs2(1)

	Ss16 = 2*ML2 + MHiggs2(1)

	Ss54 = -2*ML2 + MHiggs2(1)

	Ss15 = -4*ML2 + MHiggs2(1)

	Ss142 = 168*MBot2**2 - 36*MZ2**2

	Ss123 = 3*MBot2**2 + ML2**2

	Ss640 = MBot2**2 + MTop2**2

	Ss882 = 3*CA*CW*MBot2 - CAB*CB*MW*MZ*SW2

	Ss606 = 97*MBot2 + 111*MTop2

	Ss456 = 97*MBot2 + 92*MTop2

	Ss387 = 71*MBot2 + 111*MTop2

	Ss464 = 26*MBot2 + 37*MTop2

	Ss383 = 25*MBot2 + 29*MTop2

	Ss148 = 17*MBot2 + 16*MW2

	Ss480 = 13*MBot2 - 6*MZ2

	Ss150 = 13*MBot2 + 2*MW2

	Ss478 = 5*MBot2 + 3*MZ2

	Ss7 = 4*MBot2 - MZ2

	Ss1 = 4*MBot2 - MHiggs2(2)

	Ss13 = 4*MBot2 - MHiggs2(1)

	Ss469 = 3*MBot2 - 4*MTop2

	Ss10 = 2*MBot2 + MZ2

	Ss58 = 2*MBot2 - MZ2

	Ss2 = 2*MBot2 + MHiggs2(2)

	Ss57 = 2*MBot2 - MHiggs2(2)

	Ss14 = 2*MBot2 + MHiggs2(1)

	Ss60 = 2*MBot2 - MHiggs2(1)

	Ss643 = MBot2 + Ss32*Ss33

	Ss172 = MBot2 + 3*MW2

	Ss167 = MBot2 - 3*MW2

	Ss155 = MBot2 - 4*MW2

	Ss29 = MBot2 + MTop2

	Ss30 = MBot2 - MTop2

	Ss163 = MBot + MW

	Ss162 = MBot - MW

	Ss218 = MBot + MTop

	Ss215 = MBot - MTop

	Ss299 = MBot + MHiggs(4)

	Ss298 = MBot - MHiggs(4)

	Ss62 = -MZ2 + MHiggs2(3)

	Ss6 = 2*MTop2 + MHiggs2(3)

	Ss5 = -4*MTop2 + MHiggs2(3)

	Ss69 = 2*ML2 + MHiggs2(3)

	Ss68 = -4*ML2 + MHiggs2(3)

	Ss67 = 2*MBot2 + MHiggs2(3)

	Ss66 = -4*MBot2 + MHiggs2(3)

	Ss114 = 410 + 479*TB2

	Ss80 = 398 + 301*TB2

	Ss188 = 178 + 19*TB2

	Ss113 = 147 + 32*TB2

	Ss205 = 139 + 238*TB2

	Ss181 = 111 - 8*TB2

	Ss116 = 101 + 32*TB2

	Ss109 = 97 - 110*TB2

	Ss112 = 55 + 32*TB2

	Ss179 = 49 - 8*TB2

	Ss509 = 42 + Pi**2

	Ss118 = 35 + 173*TB2

	Ss86 = 29 + 35*TB2

	Ss117 = 28 + 97*TB2

	Ss110 = 28 + 5*TB2

	Ss79 = 25 - 44*TB2

	Ss853 = 23 - 12*Log(MTop2)

	Ss178 = 23 - 8*TB2

	Ss558 = 13 - 5*TB2

	Ss187 = 11 + 5*TB2

	Ss95 = 8 - 61*TB2

	Ss98 = 8 - 79*TB2

	Ss103 = 8 - 115*TB2

	Ss132 = 8 + 15*SA2

	Ss829 = 7 - 4*Log(ML2)

	Ss860 = 6 - Log(MTop2)

	Ss201 = 6 - 4*TB2

	Ss517 = 5 - 2*TB2

	Ss89 = 5 - 7*TB2

	Ss830 = 4 - 8*Log(MZ2)

	Ss790 = 4 - Log(MTop2)

	Ss781 = 4 - Log(ML2)

	Ss709 = 4 - Log(MBot2)

	Ss195 = 4 + 183*TB2

	Ss203 = 4 + 91*TB2

	Ss190 = 4 + 35*TB2

	Ss111 = 4 + 27*TB2

	Ss343 = 4 + 11*TB2

	Ss184 = 4 - TB2

	Ss193 = 4 - 6*TB2

	Ss284 = 4 - 19*TB2

	Ss81 = 4 - 73*TB2

	Ss866 = 4 - (12*MTop2)/MZ2

	Ss993 = 3 - Log(MStop2(2))

	Ss974 = 3 - Log(MStop2(1))

	Ss954 = 3 - Log(MStau2(2))

	Ss933 = 3 - Log(MStau2(1))

	Ss911 = 3 - Log(MSbot2(2))

	Ss890 = 3 - Log(MSbot2(1))

	Ss550 = 3 - Log(MZ2)

	Ss530 = 3 - Log(MW2)

	Ss511 = 3 - Log(MTop2)

	Ss326 = 3 - Log(ML2)

	Ss524 = 3 - Log(MHiggs2(4))

	Ss548 = 3 - Log(MHiggs2(2))

	Ss546 = 3 - Log(MHiggs2(1))

	Ss510 = 3 - Log(MBot2)

	Ss535 = 3 - Log(MHiggs2(3))

	Ss559 = 3 - TB2

	Ss519 = 3 - 6*TB2

	Ss555 = 3 - 9*TB2

	Ss879 = 3 - 2*SW2

	Ss211 = 3 - 4*SW2

	Ss48 = 3 - 8*SW2

	Ss989 = 2 - Log(MStop2(2))

	Ss970 = 2 - Log(MStop2(1))

	Ss950 = 2 - Log(MStau2(2))

	Ss929 = 2 - Log(MStau2(1))

	Ss907 = 2 - Log(MSbot2(2))

	Ss886 = 2 - Log(MSbot2(1))

	Ss833 = 2 - Log(ML2)*Log(MZ2)

	Ss809 = 2 + Log(MZ2)

	Ss483 = 2 - Log(MZ2)

	Ss454 = 2 - Log(MW2)

	Ss374 = 2 - Log(MTop2)

	Ss323 = 2 - Log(ML2)

	Ss294 = 2 - Log(MHiggs2(4))

	Ss273 = 2 - Log(MHiggs2(2))

	Ss267 = 2 - Log(MHiggs2(1))

	Ss213 = 2 - Log(MBot2)

	Ss139 = 2 - Log(MHiggs2(3))

	Ss84 = 2 + 35*TB2

	Ss185 = 2 - TB2

	Ss419 = 2 - 3*TB2

	Ss723 = 3/2.D0 - Log(ML2)

	Ss988 = 1 - Log(MStop2(2))

	Ss969 = 1 - Log(MStop2(1))

	Ss949 = 1 - Log(MStau2(2))

	Ss928 = 1 - Log(MStau2(1))

	Ss906 = 1 - Log(MSbot2(2))

	Ss885 = 1 - Log(MSbot2(1))

	Ss474 = 1 - Log(MZ2)

	Ss398 = 1 - Log(MW2)

	Ss335 = 1 - Log(MTop2)

	Ss306 = 1 - Log(ML2)

	Ss275 = 1 - Log(MHiggs2(4))

	Ss269 = 1 - Log(MHiggs2(2))

	Ss263 = 1 - Log(MHiggs2(1))

	Ss141 = 1 - Log(MBot2)

	Ss138 = 1 - Log(MHiggs2(3))

	Ss288 = 1 + 3*TB2

	Ss560 = 1 + 2*TB2

	Ss199 = 1 + TB2

	Ss96 = 1 - TB2

	Ss200 = 1 - 2*TB2

	Ss285 = 1 - 3*TB2

	Ss176 = 1 - 4*TB2

	Ss82 = 1 - 5*TB2

	Ss105 = 1 - 11*TB2

	Ss922 = 1 - 2*SW2

	Ss78 = 1 - 4*SW2

	Ss47 = 1 + 2*SW

	Ss46 = 1 - 2*SW

	Ss126 = 1 - (4*MTop2)/MZ2

	Ss501 = 1 - (4*ML2)/MZ2

	Ss496 = 1 - (4*MBot2)/MZ2

	Ss43 = 28*MW2**4 + MTop2*(81*MTop2*MW2**2 - 132*MW2**3) + 
     &   MHiggs2(4)*(-(MW2*Ss42) + Ss41*MHiggs2(4))

	Ss352 = 4*MW2**4 + MTop2*(9*MTop2*MW2**2 - 18*MW2**3) + 
     &   MHiggs2(4)*(MW2*Ss351 - Ss350*MHiggs2(4))

	Ss260 = MW2**4 + MTop2*(3*MTop2*MW2**2 - 2*MW2**3) + 
     &   MHiggs2(4)*(-(MW2*Ss259) + Ss258*MHiggs2(4))

	Ss235 = 6*MTop2**2 + MBot2*(6*MBot2 - Ss234) + 
     &   MHiggs2(4)*(-7*MTop2 + 3*MHiggs2(4))

	Ss379 = MTop2*(3*MTop2 - 5*MW2) + 2*MW2**2 + MBot2*(3*MBot2 - Ss378)

	Ss254 = MTop2*(3*MTop2 - 8*MW2) + 3*MW2**2 + MBot2*(3*MBot2 - Ss253)

	Ss240 = 3*MTop2**2 + MBot2*(3*MBot2 + 2*Ss239) + 
     &   MHiggs2(4)*(-8*MTop2 + 3*MHiggs2(4))

	Ss301 = 3*MTop2**2 + MBot2*(3*MBot2 - Ss300) + 
     &   MHiggs2(4)*(-5*MTop2 + 2*MHiggs2(4))

	Ss128 = 17*MZ2**3 + MW2*(32*MW2*MZ2 - 40*MZ2**2) + MTop2*Ss127

	Ss498 = 5*MZ2**3 + MW2*(8*MW2*MZ2 - 4*MZ2**2) + MBot2*Ss497

	Ss503 = 5*MZ2**3 + MW2*(8*MW2*MZ2 - 12*MZ2**2) + ML2*Ss502

	Ss77 = 27*MTop2*MW2**2 + 14*MW2**3 + 
     &   MHiggs2(4)*(-(MW2*Ss76) + Ss75*MHiggs2(4))

	Ss257 = 3*MTop2*MW2**2 + 2*MW2**3 + 
     &   MHiggs2(4)*(-(MW2*Ss243) + Ss256*MHiggs2(4))

	Ss74 = 26*MW2**2 + 9*MBot2*MHiggs2(4) + 
     &   MHiggs2(4)*(-2*Ss73 + 8*MHiggs2(4))

	Ss795 = 1 + (-2 - Log(MTop2))*Log(MTop2) - 2*Ss335*Log(MHiggs2(1))

	Ss786 = 1 + (-2 - Log(ML2))*Log(ML2) - 2*Ss306*Log(MHiggs2(1))

	Ss1000 = 3/2.D0*Ss998 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Log(MStop2(2))*MStop2(2)

	Ss992 = 3/2.D0*Ss991 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Log(MStop2(2))*MStop2(2)

	Ss984 = 3/2.D0*Ss982 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Log(MStop2(1))*MStop2(1)

	Ss973 = 3/2.D0*Ss972 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Log(MStop2(1))*MStop2(1)

	Ss961 = 3/2.D0*Ss959 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Log(MStau2(2))*MStau2(2)

	Ss953 = 3/2.D0*Ss952 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Log(MStau2(2))*MStau2(2)

	Ss943 = 3/2.D0*Ss941 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Log(MStau2(1))*MStau2(1)

	Ss932 = 3/2.D0*Ss931 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Log(MStau2(1))*MStau2(1)

	Ss918 = 3/2.D0*Ss916 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Log(MSbot2(2))*MSbot2(2)

	Ss910 = 3/2.D0*Ss909 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Log(MSbot2(2))*MSbot2(2)

	Ss900 = 3/2.D0*Ss898 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Log(MSbot2(1))*MSbot2(1)

	Ss889 = 3/2.D0*Ss888 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Log(MSbot2(1))*MSbot2(1)

	Ss322 = 3/2.D0*Ss70 - ML2*Log(ML2) - Log(MHiggs2(4))*MHiggs2(4)

	Ss317 = 3/2.D0*Ss69 - 2*ML2*Log(ML2) - Log(MHiggs2(3))*MHiggs2(3)

	Ss212 = 3/2.D0*Ss67 - 2*MBot2*Log(MBot2) - Log(MHiggs2(3))*MHiggs2(3)

	Ss368 = 3/2.D0*Ss6 - 2*MTop2*Log(MTop2) - Log(MHiggs2(3))*MHiggs2(3)

	Ss319 = 3/2.D0*Ss4 - 2*ML2*Log(ML2) - Log(MHiggs2(2))*MHiggs2(2)

	Ss431 = 3/2.D0*Ss31 - ML2*Log(ML2) - MW2*Log(MW2)

	Ss381 = 3*Ss255*Ss30**2 + MBot2*MW2*Ss376 + 3*MW2**2*Ss377

	Ss371 = 3/2.D0*Ss29 - MBot2*Log(MBot2) - MTop2*Log(MTop2)

	Ss370 = 3/2.D0*Ss20 - 2*MTop2*Log(MTop2) - Log(MHiggs2(2))*MHiggs2(2)

	Ss272 = 3/2.D0*Ss2 - 2*MBot2*Log(MBot2) - Log(MHiggs2(2))*MHiggs2(2)

	Ss369 = 3/2.D0*Ss19 - 2*MTop2*Log(MTop2) - Log(MHiggs2(1))*MHiggs2(1)

	Ss318 = 3/2.D0*Ss16 - 2*ML2*Log(ML2) - Log(MHiggs2(1))*MHiggs2(1)

	Ss266 = 3/2.D0*Ss14 - 2*MBot2*Log(MBot2) - Log(MHiggs2(1))*MHiggs2(1)

	Ss482 = 3/2.D0*Ss12 - 2*MTop2*Log(MTop2) - MZ2*Log(MZ2)

	Ss481 = 3/2.D0*Ss11 - 2*ML2*Log(ML2) - MZ2*Log(MZ2)

	Ss479 = 3/2.D0*Ss10 - 2*MBot2*Log(MBot2) - MZ2*Log(MZ2)

	Ss280 = 2*MW2**2 + MBot2*MHiggs2(4) - 2*Ss279*MHiggs2(4)

	Ss367 = -Ss9 + MHiggs2(2)

	Ss366 = -Ss9 + MHiggs2(1)

	Ss121 = 3*MTop2**2 + Ss123

	Ss59 = 4*MBot2 - Ss143

	Ss61 = 4*MBot2 - Ss145

	Ss633 = 3*MBot2 + ML2 + 3*MTop2

	Ss653 = 2*MBot2 - 9*ML2 + MW2

	Ss34 = MBot2 + Ss36

	Ss486 = MBot2 + Ss35

	Ss575 = MBot2 - Ss35

	Ss576 = MBot2 - Ss36

	Ss21 = MBot2 + Ss23

	Ss564 = MBot2 + Ss22

	Ss565 = MBot2 - Ss22

	Ss566 = MBot2 - Ss23

	Ss220 = MBot + Ss27

	Ss217 = MBot + Ss26

	Ss219 = MBot - Ss26

	Ss216 = MBot - Ss27

	Ss327 = 14 + Pi**2 - 4*Ss326*Log(ML2)

	Ss990 = 12 + Pi**2 - 6*Ss989*Log(MStop2(2))

	Ss971 = 12 + Pi**2 - 6*Ss970*Log(MStop2(1))

	Ss951 = 12 + Pi**2 - 6*Ss950*Log(MStau2(2))

	Ss930 = 12 + Pi**2 - 6*Ss929*Log(MStau2(1))

	Ss908 = 12 + Pi**2 - 6*Ss907*Log(MSbot2(2))

	Ss887 = 12 + Pi**2 - 6*Ss886*Log(MSbot2(1))

	Ss484 = 12 + Pi**2 - 6*Ss483*Log(MZ2)

	Ss455 = 12 + Pi**2 - 6*Ss454*Log(MW2)

	Ss375 = 12 + Pi**2 - 6*Ss374*Log(MTop2)

	Ss324 = 12 + Pi**2 - 6*Ss323*Log(ML2)

	Ss295 = 12 + Pi**2 - 6*Ss294*Log(MHiggs2(4))

	Ss274 = 12 + Pi**2 - 6*Ss273*Log(MHiggs2(2))

	Ss268 = 12 + Pi**2 - 6*Ss267*Log(MHiggs2(1))

	Ss214 = 12 + Pi**2 - 6*Ss213*Log(MBot2)

	Ss140 = 12 + Pi**2 - 6*Ss139*Log(MHiggs2(3))

	Ss801 = 5 - Log(MTop2)**2 - 2*Ss374*Log(MZ2)

	Ss831 = 5 - Ss830*Log(ML2) - 8*Log(MZ2)

	Ss716 = 5 - Log(MBot2)**2 - 2*Ss213*Log(MZ2)

	Ss861 = 3 + Ss860*Log(MTop2) - 6*Log(MZ2)

	Ss746 = Ss454 - Log(MTop2)

	Ss745 = Ss454 - Log(MBot2)

	Ss735 = 2 - Ss513

	Ss712 = -(Ss14*Ss509) + 
     &   MBot2*(2*Ss214 + 24*Ss510*Log(MBot2)) + 
     &   (Ss268 - 12*Ss711 + 6*Log(MBot2/MHiggs2(1))**2 + 
     &      12*Ss546*Log(MHiggs2(1)))*MHiggs2(1) - 
     &   12*TF(MBot2,MBot2,MHiggs2(1))

	Ss629 = -(Ss12*Ss509) + 
     &   MTop2*(2*Ss375 + 24*Ss511*Log(MTop2)) + 
     &   MZ2*(Ss484 + 6*Log(MTop2/MZ2)**2 + 12*Ss550*Log(MZ2)) - 
     &   12*TF(MTop2,MTop2,MZ2)

	Ss424 = 8*MW2**2 + MBot2*(27*MW2 + 18*MBot2*SA2) + 
     &   ML2*(9*MW2 + 6*ML2*SA2) + MTop2*(27*MW2 + 18*MTop2*SA2)

	Ss620 = -(Ss4*Ss509) + 24*ML2*Ss326*Log(ML2) + 
     &   (Ss274 + 12*Ss548*Log(MHiggs2(2)) + 
     &      6*Log(MHiggs2(2)/ML2)**2)*MHiggs2(2) - 
     &   12*TF(MHiggs2(2),ML2,ML2)

	Ss616 = -(Ss16*Ss509) + 24*ML2*Ss326*Log(ML2) + 
     &   (Ss268 + 12*Ss546*Log(MHiggs2(1)) + 
     &      6*Log(MHiggs2(1)/ML2)**2)*MHiggs2(1) - 
     &   12*TF(MHiggs2(1),ML2,ML2)

	Ss533 = -3*Ss28*Ss531*Ss532 - 12*ML2*Ss326*Log(ML2) + 
     &   Ss31*(Ss509 - 3*Log(ML2/MW2)**2) - 12*MW2*Ss530*Log(MW2)

	Ss514 = -3*Ss30*Ss512*Ss513 - 12*MBot2*Ss510*Log(MBot2) + 
     &   Ss29*(Ss509 - 3*Log(MBot2/MTop2)**2) - 12*MTop2*Ss511*Log(MTop2)

	Ss634 = 18*MBot2**2 + 6*ML2**2 + 18*MTop2**2 + MW2*(8*MW2 + 9*Ss633)

	Ss432 = 3/2.D0*Ss34 - MBot2*Log(MBot2) - MTop2*Log(MTop2) - MW2*Log(MW2)

	Ss372 = 3/2.D0*Ss21 - MBot2*Log(MBot2) - MTop2*Log(MTop2) - 
     &   Log(MHiggs2(4))*MHiggs2(4)

	Ss491 = 2*MBot*MTop + Ss486

	Ss490 = -2*MBot*MTop + Ss486

	Ss847 = 5 + (2*Ss809 - Log(ML2))*Log(ML2) - 8*Log(MZ2)

	Ss810 = 5 + (2*Ss809 - Log(MBot2))*Log(MBot2) - 8*Log(MZ2)

	Ss713 = 1 + (-2*Ss263 - Log(MBot2))*Log(MBot2) - 2*Log(MHiggs2(1))

	Ss733 = 1/12.D0*(MW2*Ss455 - Ss31*Ss509) + 
     &   1/4.D0*(Ss28*Ss531*Ss532 + Ss31*Log(ML2/MW2)**2) + 
     &   ML2*Ss326*Log(ML2) + MW2*(Ss306*Ss398 + Ss530*Log(MW2)) - 
     &   TF(0.D0,ML2,MW2)

	Ss1001 = 1/12.D0*(Ss509*Ss998) - 
     &   1/2.D0*(Log(MHiggs2(1)/MStop2(2))**2*MHiggs2(1)) - 
     &   Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Ss993*Log(MStop2(2))*MStop2(2) + 
     &   TF(MHiggs2(1),MStop2(2),MStop2(2))

	Ss994 = 1/12.D0*(Ss509*Ss991) - 
     &   1/2.D0*(Log(MHiggs2(2)/MStop2(2))**2*MHiggs2(2)) - 
     &   Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Ss993*Log(MStop2(2))*MStop2(2) + 
     &   TF(MHiggs2(2),MStop2(2),MStop2(2))

	Ss985 = 1/12.D0*(Ss509*Ss982) - 
     &   1/2.D0*(Log(MHiggs2(1)/MStop2(1))**2*MHiggs2(1)) - 
     &   Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Ss974*Log(MStop2(1))*MStop2(1) + 
     &   TF(MHiggs2(1),MStop2(1),MStop2(1))

	Ss975 = 1/12.D0*(Ss509*Ss972) - 
     &   1/2.D0*(Log(MHiggs2(2)/MStop2(1))**2*MHiggs2(2)) - 
     &   Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Ss974*Log(MStop2(1))*MStop2(1) + 
     &   TF(MHiggs2(2),MStop2(1),MStop2(1))

	Ss962 = 1/12.D0*(Ss509*Ss959) - 
     &   1/2.D0*(Log(MHiggs2(1)/MStau2(2))**2*MHiggs2(1)) - 
     &   Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Ss954*Log(MStau2(2))*MStau2(2) + 
     &   TF(MHiggs2(1),MStau2(2),MStau2(2))

	Ss955 = 1/12.D0*(Ss509*Ss952) - 
     &   1/2.D0*(Log(MHiggs2(2)/MStau2(2))**2*MHiggs2(2)) - 
     &   Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Ss954*Log(MStau2(2))*MStau2(2) + 
     &   TF(MHiggs2(2),MStau2(2),MStau2(2))

	Ss944 = 1/12.D0*(Ss509*Ss941) - 
     &   1/2.D0*(Log(MHiggs2(1)/MStau2(1))**2*MHiggs2(1)) - 
     &   Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Ss933*Log(MStau2(1))*MStau2(1) + 
     &   TF(MHiggs2(1),MStau2(1),MStau2(1))

	Ss934 = 1/12.D0*(Ss509*Ss931) - 
     &   1/2.D0*(Log(MHiggs2(2)/MStau2(1))**2*MHiggs2(2)) - 
     &   Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Ss933*Log(MStau2(1))*MStau2(1) + 
     &   TF(MHiggs2(2),MStau2(1),MStau2(1))

	Ss919 = 1/12.D0*(Ss509*Ss916) - 
     &   1/2.D0*(Log(MHiggs2(1)/MSbot2(2))**2*MHiggs2(1)) - 
     &   Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Ss911*Log(MSbot2(2))*MSbot2(2) + 
     &   TF(MHiggs2(1),MSbot2(2),MSbot2(2))

	Ss912 = 1/12.D0*(Ss509*Ss909) - 
     &   1/2.D0*(Log(MHiggs2(2)/MSbot2(2))**2*MHiggs2(2)) - 
     &   Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Ss911*Log(MSbot2(2))*MSbot2(2) + 
     &   TF(MHiggs2(2),MSbot2(2),MSbot2(2))

	Ss901 = 1/12.D0*(Ss509*Ss898) - 
     &   1/2.D0*(Log(MHiggs2(1)/MSbot2(1))**2*MHiggs2(1)) - 
     &   Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   2*Ss890*Log(MSbot2(1))*MSbot2(1) + 
     &   TF(MHiggs2(1),MSbot2(1),MSbot2(1))

	Ss891 = 1/12.D0*(Ss509*Ss888) - 
     &   1/2.D0*(Log(MHiggs2(2)/MSbot2(1))**2*MHiggs2(2)) - 
     &   Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   2*Ss890*Log(MSbot2(1))*MSbot2(1) + 
     &   TF(MHiggs2(2),MSbot2(1),MSbot2(1))

	Ss540 = 1/12.D0*(Ss509*Ss69) - 
     &   1/2.D0*(Log(MHiggs2(3)/ML2)**2*MHiggs2(3)) - 
     &   2*ML2*Ss326*Log(ML2) - Ss535*Log(MHiggs2(3))*MHiggs2(3) + 
     &   TF(MHiggs2(3),ML2,ML2)

	Ss536 = 1/12.D0*(Ss509*Ss67) - 
     &   1/2.D0*(Log(MHiggs2(3)/MBot2)**2*MHiggs2(3)) - 
     &   2*MBot2*Ss510*Log(MBot2) - Ss535*Log(MHiggs2(3))*MHiggs2(3) + 
     &   TF(MHiggs2(3),MBot2,MBot2)

	Ss543 = 1/12.D0*(Ss509*Ss6) - 
     &   1/2.D0*(Log(MHiggs2(3)/MTop2)**2*MHiggs2(3)) - 
     &   2*MTop2*Ss511*Log(MTop2) - Ss535*Log(MHiggs2(3))*MHiggs2(3) + 
     &   TF(MHiggs2(3),MTop2,MTop2)

	Ss528 = 3*Ss525*Ss526*Ss527 - Ss509*Ss70 + 
     &   12*ML2*Ss326*Log(ML2) + 3*Ss70*Log(MHiggs2(4)/ML2)**2 + 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4)

	Ss621 = 1/12.D0*(Ss4*Ss509) - 
     &   1/2.D0*(Log(MHiggs2(2)/ML2)**2*MHiggs2(2)) - 
     &   2*ML2*Ss326*Log(ML2) - Ss548*Log(MHiggs2(2))*MHiggs2(2) + 
     &   TF(MHiggs2(2),ML2,ML2)

	Ss623 = 1/12.D0*(Ss20*Ss509) - 
     &   1/2.D0*(Log(MHiggs2(2)/MTop2)**2*MHiggs2(2)) - 
     &   2*MTop2*Ss511*Log(MTop2) - Ss548*Log(MHiggs2(2))*MHiggs2(2) + 
     &   TF(MHiggs2(2),MTop2,MTop2)

	Ss549 = 1/12.D0*(Ss2*Ss509) - 
     &   1/2.D0*(Log(MBot2/MHiggs2(2))**2*MHiggs2(2)) - 
     &   2*MBot2*Ss510*Log(MBot2) - Ss548*Log(MHiggs2(2))*MHiggs2(2) + 
     &   TF(MBot2,MBot2,MHiggs2(2))

	Ss619 = 1/12.D0*(Ss19*Ss509) - 
     &   1/2.D0*(Log(MHiggs2(1)/MTop2)**2*MHiggs2(1)) - 
     &   2*MTop2*Ss511*Log(MTop2) - Ss546*Log(MHiggs2(1))*MHiggs2(1) + 
     &   TF(MHiggs2(1),MTop2,MTop2)

	Ss617 = 1/12.D0*(Ss16*Ss509) - 
     &   1/2.D0*(Log(MHiggs2(1)/ML2)**2*MHiggs2(1)) - 
     &   2*ML2*Ss326*Log(ML2) - Ss546*Log(MHiggs2(1))*MHiggs2(1) + 
     &   TF(MHiggs2(1),ML2,ML2)

	Ss547 = 1/12.D0*(Ss14*Ss509) - 
     &   1/2.D0*(Log(MBot2/MHiggs2(1))**2*MHiggs2(1)) - 
     &   2*MBot2*Ss510*Log(MBot2) - Ss546*Log(MHiggs2(1))*MHiggs2(1) + 
     &   TF(MBot2,MBot2,MHiggs2(1))

	Ss628 = 1/12.D0*(Ss12*Ss509) - 
     &   1/2.D0*(MZ2*Log(MTop2/MZ2)**2) - 2*MTop2*Ss511*Log(MTop2) - 
     &   MZ2*Ss550*Log(MZ2) + TF(MTop2,MTop2,MZ2)

	Ss625 = 1/12.D0*(Ss11*Ss509) - 
     &   1/2.D0*(MZ2*Log(ML2/MZ2)**2) - 2*ML2*Ss326*Log(ML2) - 
     &   MZ2*Ss550*Log(MZ2) + TF(ML2,ML2,MZ2)

	Ss551 = 1/12.D0*(Ss10*Ss509) - 
     &   1/2.D0*(MZ2*Log(MBot2/MZ2)**2) - 2*MBot2*Ss510*Log(MBot2) - 
     &   MZ2*Ss550*Log(MZ2) + TF(MBot2,MBot2,MZ2)

	Ss622 = 2*Ss621 - 1/3.D0*(ML2*Ss324) + 3/2.D0*Ss4 - 
     &   1/6.D0*(Ss274*MHiggs2(2)) + ML2*(-6*Ss306 - 2*Log(ML2)) + 
     &   (Ss269 - Log(MHiggs2(2)))*MHiggs2(2)

	Ss624 = 2*Ss623 + 3/2.D0*Ss20 - 1/3.D0*(MTop2*Ss375) - 
     &   1/6.D0*(Ss274*MHiggs2(2)) + MTop2*(-6*Ss335 - 2*Log(MTop2)) + 
     &   (Ss269 - Log(MHiggs2(2)))*MHiggs2(2)

	Ss618 = 2*Ss617 + 3/2.D0*Ss16 - 1/3.D0*(ML2*Ss324) - 
     &   1/6.D0*(Ss268*MHiggs2(1)) + ML2*(-6*Ss306 - 2*Log(ML2)) + 
     &   (Ss263 - Log(MHiggs2(1)))*MHiggs2(1)

	Ss963 = 2*Ss962 + 3/2.D0*Ss959 - 1/6.D0*(Ss268*MHiggs2(1)) - 
     &   1/3.D0*(Ss951*MStau2(2)) + 
     &   (Ss263 - Log(MHiggs2(1)))*MHiggs2(1) - 
     &   6*Ss949*MStau2(2) - 2*Log(MStau2(2))*MStau2(2)

	Ss945 = 2*Ss944 + 3/2.D0*Ss941 - 1/6.D0*(Ss268*MHiggs2(1)) - 
     &   1/3.D0*(Ss930*MStau2(1)) + 
     &   (Ss263 - Log(MHiggs2(1)))*MHiggs2(1) - 
     &   6*Ss928*MStau2(1) - 2*Log(MStau2(1))*MStau2(1)

	Ss920 = 2*Ss919 + 3/2.D0*Ss916 - 1/6.D0*(Ss268*MHiggs2(1)) - 
     &   1/3.D0*(Ss908*MSbot2(2)) + 
     &   (Ss263 - Log(MHiggs2(1)))*MHiggs2(1) - 
     &   6*Ss906*MSbot2(2) - 2*Log(MSbot2(2))*MSbot2(2)

	Ss902 = 2*Ss901 + 3/2.D0*Ss898 - 1/6.D0*(Ss268*MHiggs2(1)) - 
     &   1/3.D0*(Ss887*MSbot2(1)) + 
     &   (Ss263 - Log(MHiggs2(1)))*MHiggs2(1) - 
     &   6*Ss885*MSbot2(1) - 2*Log(MSbot2(1))*MSbot2(1)

	Ss577 = Ss34*Ss509 - 12*MBot2*Ss510*Log(MBot2) - 
     &   3*Ss486*Log(MBot2/MTop2)**2 - 12*MTop2*Ss511*Log(MTop2) - 
     &   3*Ss575*Log(MBot2/MW2)**2 + 3*Ss576*Log(MTop2/MW2)**2 - 
     &   12*MW2*Ss530*Log(MW2)

	Ss567 = Ss21*Ss509 - 12*MBot2*Ss510*Log(MBot2) - 
     &   3*Ss565*Log(MBot2/MTop2)**2 - 12*MTop2*Ss511*Log(MTop2) - 
     &   3*Ss564*Log(MBot2/MHiggs2(4))**2 + 
     &   3*Ss566*Log(MHiggs2(4)/MTop2)**2 - 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4)

	Ss587 = Ss32**4*Ss33**4*Ss583 + 
     &   MBot2*(-(Ss32**2*Ss33**2*Ss584) + 
     &      MBot2*(MBot2*(MBot2*(MBot2*
     &                (MBot2*(MBot2*(196*MBot2 - Ss579) + Ss580) - Ss581)
     &                 + Ss582) - Ss585) + Ss586))

	Ss596 = MTop2*Ss32**4*Ss33**4*Ss590 + 
     &   MBot2*(-(Ss32**2*Ss33**2*Ss593) + 
     &      MBot2*(MBot2*(MBot2*(MBot2*
     &                (MBot2*(MBot2*(124*MBot2 - Ss588) + Ss589) - Ss591)
     &                 + Ss592) - Ss594) + Ss595))

	Ss605 = Ss32**5*Ss33**5*Ss599 + 
     &   MBot2*(-(Ss32**2*Ss33**2*Ss603) + 
     &      MBot2*(MBot2*(MBot2*(MBot2*
     &                (MBot2*(MBot2*(76*MBot2 - Ss597) + Ss598) - Ss600)+
     &                 Ss601) - Ss602) + Ss604))

	Ss773 = MTop2**2*Ss32**5*Ss33**5*Ss507 + 
     &   MBot2*(-(MTop2*Ss32**2*Ss33**2*Ss769) + 
     &      MBot2*(MBot2*(MBot2*(MBot2*
     &                (MBot2*(MBot2*(19*MBot2 - Ss766) + Ss767) - Ss768)+
     &                 Ss770) - Ss771) + Ss772))

	Ss407 = -(Ss32**5*Ss33**5*Ss401) + 
     &   MBot2*(2*Ss32**2*Ss33**2*Ss404 + 
     &      MBot2*(MBot2*(MBot2*(MBot2*(MBot2*(76*MBot2 - Ss400) + 6*Ss402) - 
     &               Ss403) + 4*Ss405) - Ss406))

	Ss397 = -(Ss32**3*Ss33**3*Ss393) + 
     &   MBot2*(MBot2*(MBot2*(MBot2*
     &             (MBot2*(MBot2*(52*MBot2 - Ss390) + 3*Ss391) - Ss392)+
     &              2*Ss394) - Ss395) + Ss32*Ss33*Ss396)

	Ss754 = -(MTop2*Ss32**3*Ss33**3*Ss749) + 
     &   MBot2*(MBot2*(MBot2*(MBot2*
     &             (MBot2*(MBot2*(19*MBot2 - Ss747) + Ss748) - Ss750) + 
     &            Ss751) - Ss752) + Ss32*Ss33*Ss753)

	Ss763 = -(MTop2*Ss32**5*Ss33**5*Ss507) + 
     &   MBot2*(Ss32**2*Ss33**2*Ss760 + 
     &      MBot2*(MBot2*(MBot2*(MBot2*(MBot2*(19*MBot2 - Ss757) + 3*Ss758) - 
     &               Ss759) + Ss761) - Ss762))

	Ss388 = 38*Ss215**6*Ss218**6 + 
     &   MW2*(-7*Ss215**4*Ss218**4*Ss383 + 
     &      MW2*(2*Ss215**2*Ss218**2*Ss384 + 
     &         MW2*(-2*Ss385 + 
     &            MW2*(2*Ss386 + MW2*(14*MW2 - Ss387)))))

	Ss729 = ML2*(-218*MW2**3 + 
     &      ML2*(ML2*(1564*ML2 - 8077*MW2) - 7959*MW2**2 - 
     &         12*Ss727*Log(ML2))) + 
     &   MW2*(-286*MW2**3 - 96*Ss728*Log(MW2))

	Ss334 = Ss32**3*Ss33**3*Ss330 + 
     &   MBot2*(-2*Ss32*Ss33*Ss332 + 
     &      MBot2*(MBot2*(MBot2*(MBot2*(760*MBot2 - 6*Ss328) + 3*Ss329) - 
     &            2*Ss331) + 2*Ss333))

	Ss742 = Ss32**3*Ss33**3*Ss738 + 
     &   MBot2*(-(Ss32*Ss33*Ss740) + 
     &      MBot2*(MBot2*(MBot2*(MBot2*(38*MBot2 - Ss736) + Ss737) - 
     &            2*Ss739) + 2*Ss741))

	Ss671 = MW2**4*(2*MW2 - 11*MZ2) + 
     &   MTop2*(-34*MW2**3*Ss649*Ss650 + 
     &      MTop2*(MTop2*(5*MTop2*Ss635 - 4*MW2*Ss635) - MW2**2*Ss670))

	Ss614 = 19*Ss215**8*Ss218**8 + 
     &   MW2*(-(Ss215**6*Ss218**6*Ss606) + 
     &      MW2*(10*Ss215**4*Ss218**4*Ss607 + 
     &         MW2*(-2*Ss215**2*Ss218**2*Ss608 + MW2*Ss613)))

	Ss131 = 60*MTop2*MZ2**3 + 
     &   MW2*(MW2*(MW2*(64*MW2 - 32*MZ2) + 91*MZ2**2) - 
     &      3*MZ2**2*Ss130)

	Ss494 = MW2**3*(28*MW2 - 13*MZ2) + 
     &   MTop2*(-2*MW2**2*Ss433 + MTop2*(-4*MTop2*Ss433 + MW2*Ss433))

	Ss679 = MW2**3*(25*MW2 - 7*MZ2) + 
     &   MTop2*(MTop2*(-6*MW2*Ss649*Ss650 - MTop2*Ss663) + MW2**2*Ss670)

	Ss667 = MW2**3*(11*MW2 + 7*MZ2) + 
     &   MTop2*(MTop2*(-10*MTop2*Ss635 + 6*MW2*Ss649*Ss650) - 
     &      MW2**2*Ss666)

	Ss676 = MTop2*(-50*MW2**3*Ss649*Ss650 + 
     &      MTop2*(MTop2*(12*MW2**2 - 5*MTop2*Ss635) - MW2**2*Ss666)) + 
     &   MW2**4*Ss675

	Ss703 = MTop2*Ss279*Ss35**2 + 
     &   MBot2*(MBot2*(MBot2*(MBot2 - 4*MTop2) + Ss701) - 2*Ss702)

	Ss697 = MW2**3*Ss649*Ss650*Ss693 + 
     &   MTop2*(MTop2*(-36*MTop2*MZ2*Ss694 + 9*MW2*MZ2*Ss695) + 
     &      4*MW2**2*Ss696)

	Ss669 = 2*MW2**3*Ss649*Ss650 + 
     &   MTop2*(-(MW2**2*Ss668) + MTop2*(MTop2*Ss635 + MW2*MZ2*SW2))

	Ss678 = 6*MW2**3*Ss649*Ss650 + 
     &   MTop2*(MTop2*(MTop2*Ss635 + 9*MW2*Ss649*Ss650) - MW2**2*Ss677)

	Ss437 = -6*MW2**3*MZ2*SW2 + 
     &   MTop2*(-9*MW2**2*MZ2 + MTop2*(MTop2*Ss433 - 12*MW2*MZ2*SA2*SW2))

	Ss438 = 6*MW2**3*MZ2*SW2 + 
     &   MTop2*(9*MW2**2*MZ2 + 
     &      MTop2*(-(MTop2*Ss433) + 30*MW2*MZ2*SA2*SW2))

	Ss652 = ML2**2*(ML2*Ss635 + 4*MW2*Ss649*Ss650) + 
     &   MW2**3*(12*MW2 - 9*Ss651)

	Ss485 = MW2**3*(10*MW2 - 4*MZ2) + 
     &   ML2*(ML2**2*Ss433 - 9*MW2**2*Ss443)

	Ss451 = -3*MW2**3*Ss450 + 
     &   ML2*(9*MW2**3 + ML2*(-(ML2*Ss433) + 6*MW2*MZ2*SW2))

	Ss493 = MW2**2*(18*MW2 - 9*MZ2) + MTop2*(-6*MTop2*Ss433 - MW2*Ss433)

	Ss636 = MW2**2*(10*MW2 - 4*MZ2) + ML2*(-3*MW2*Ss443 - ML2*Ss635)

	Ss492 = MW2**2*(10*MW2 - 4*MZ2) + MTop2*(-(MTop2*Ss433) - MW2*Ss433)

	Ss674 = MW2**2*(7*MW2 - MZ2) + MTop2*(-10*MTop2*Ss635 - 4*MW2*Ss635)

	Ss87 = MTop2*(MW2**2*Ss81 + MTop2*(162*MTop2 + 3*MW2*Ss86)) + 
     &   32*MW2**3*TB2

	Ss85 = MTop2*(MW2**2*Ss81 + MTop2*(81*MTop2 + 3*MW2*Ss84)) + 
     &   32*MW2**3*TB2

	Ss83 = MTop2*(MW2**2*Ss81 + MTop2*(54*MTop2 - 21*MW2*Ss82)) + 
     &   32*MW2**3*TB2

	Ss418 = MTop2*(MW2**2*Ss284 + 
     &      MTop2*(18*MTop2 + 9*MW2*Ss288)) + 8*MW2**3*TB2

	Ss417 = MTop2*(-(MW2**2*Ss284) + 
     &      MTop2*(18*MTop2 + 27*MW2*Ss96)) - 8*MW2**3*TB2

	Ss420 = MTop2*(-(MW2**2*Ss284) + 
     &      MTop2*(9*MTop2 + 9*MW2*Ss419)) - 8*MW2**3*TB2

	Ss40 = -51*MTop2*MW2**2 + 
     &   MHiggs2(4)*(Ss39 + MHiggs2(4)*(-2*Ss38 + 28*MHiggs2(4)))

	Ss249 = 3*MTop2*MW2**2 + 
     &   MHiggs2(4)*(-2*MW2*Ss248 + 
     &      MHiggs2(4)*(Ss247 + 5*MHiggs2(4)))

	Ss278 = -9*MTop2*MW2**2 + 
     &   MHiggs2(4)*(Ss277 + MHiggs2(4)*(-2*Ss276 + 4*MHiggs2(4)))

	Ss233 = MBot2*(MBot2*(12*MBot2 - Ss230) + Ss231) + MTop2*Ss232

	Ss107 = MTop2*(27*MTop2*Ss96 - 3*MW2*Ss98) + 83*MW2**2*TB2

	Ss413 = 17*MW2**2*TB2 + MTop2*(9*MTop2*Ss96 + 63*MW2*TB2)

	Ss414 = -35*MW2**2*TB2 + MTop2*(9*MTop2*Ss96 - 27*MW2*TB2)

	Ss415 = -61*MW2**2*TB2 + MTop2*(9*MTop2*Ss96 - 72*MW2*TB2)

	Ss282 = 26*MW2**2*TB2 + MTop2*(9*MTop2*Ss96 - 90*MW2*TB2)

	Ss97 = MTop2*(6*MW2*Ss95 + 27*MTop2*Ss96) + 110*MW2**2*TB2

	Ss88 = MW2**2*Ss87 + MHiggs2(4)*(-2*MW2*Ss85 + Ss83*MHiggs2(4))

	Ss92 = MTop2*(2*MW2*Ss81 + 21*MTop2*Ss82) - 32*MW2**2*TB2

	Ss571 = MTop2*(MW2*Ss558 - 3*MTop2*Ss559) + MW2**2*Ss560

	Ss561 = -3*MTop2**2*Ss559 + 
     &   MHiggs2(4)*(MTop2*Ss558 + Ss560*MHiggs2(4))

	Ss444 = MTop2*(MTop2**2*Ss433 - 9*MW2**2*Ss443) - 2*MW2**3*MZ2*SW2

	Ss421 = -(MW2**2*Ss418) + 
     &   MHiggs2(4)*(-2*MW2*Ss420 + Ss417*MHiggs2(4))

	Ss416 = MW2**2*Ss413 + MHiggs2(4)*(MW2*Ss415 - Ss414*MHiggs2(4))

	Ss355 = MW2**2*Ss354 + MHiggs2(4)*(MW2*Ss341 - Ss353*MHiggs2(4))

	Ss832 = 4*ML2*Ss324 + MZ2*(-2*Ss484 - 3*Ss831)

	Ss461 = Ss459 + MW2*(2*MW2*Ss29 - Ss460)

	Ss289 = MTop2*(2*MW2*Ss284 - 9*MTop2*Ss288) - 8*MW2**2*TB2

	Ss261 = MBot2*(3*MBot2*Ss255 - 2*Ss257) + Ss260

	Ss202 = MW2**2*Ss200 + MTop2*(2*MW2*Ss199 + MTop2*Ss201)

	Ss206 = -(MW2**2*Ss203) + MTop2*(72*MTop2*Ss199 + MW2*Ss205)

	Ss194 = MTop2*(MTop2*Ss193 + 2*MW2*Ss96) - MW2**2*TB2

	Ss189 = 5*MW2**2*Ss187 + MTop2*(-(MW2*Ss188) - 72*MTop2*Ss96)

	Ss570 = MTop2*(MW2*Ss185 - 9*MTop2*Ss96) + MW2**2*TB2

	Ss191 = MTop2*(36*MTop2*Ss185 + 5*MW2*Ss187) - MW2**2*Ss190

	Ss557 = -9*MTop2**2*Ss96 + MHiggs2(4)*(MTop2*Ss185 + TB2*MHiggs2(4))

	Ss186 = MW2**2*Ss185 + MTop2*(MTop2*Ss184 - 2*MW2*Ss96)

	Ss119 = -(MBot2**2*Ss117) + MTop2*(MTop2*Ss116 + MBot2*Ss118)

	Ss106 = MTop2*(48*MW2*Ss105 - 81*MTop2*Ss96) - 139*MW2**2*TB2

	Ss104 = MTop2*(3*MW2*Ss103 - 81*MTop2*Ss96) - 29*MW2**2*TB2

	Ss439 = MW2**3*(2*MW2*SA2 - 2*MZ2*SA2) + CA2*MW2*MZ2*Ss435*SW2

	Ss434 = MW2**2*(6*MW2 - 3*MZ2) - MTop2**2*Ss433

	Ss665 = MTop2*(12*MW2**2 + MTop2*Ss663) - MW2**2*Ss664

	Ss673 = MW2*(15*MW2 - 9*MZ2) + 5*MTop2*Ss635

	Ss641 = MW2*(6*MW2 - 3*MZ2) - 2*MTop2*Ss635

	Ss700 = MW2*(MW2 - 13*MZ2) - 4*MTop2*Ss694

	Ss182 = MW2**2*Ss179 + MTop2*(72*MTop2 - MW2*Ss181)

	Ss91 = MW2**2*Ss81 + MTop2*(54*MTop2 + 30*MW2*Ss89)

	Ss180 = MW2**2*Ss178 + MTop2*(36*MTop2 - MW2*Ss179)

	Ss345 = MTop2*(18*MTop2 - 18*MW2*Ss288) - MW2**2*Ss343

	Ss356 = MTop2*(18*MTop2 + 9*MW2*Ss285) + MW2**2*Ss343

	Ss286 = MW2**2*Ss284 + MTop2*(18*MTop2 + 18*MW2*Ss285)

	Ss339 = MTop2*(12*MTop2 - 3*MZ2) - Ss337*MHiggs2(2)

	Ss338 = MTop2*(12*MTop2 - 3*MZ2) - Ss337*MHiggs2(1)

	Ss357 = MW2**2*Ss343 + MTop2*(9*MTop2 - 27*MW2*TB2)

	Ss869 = 4*MTop2**2 + MZ2*(MZ2 + 2*Ss45*Log(MZ2))

	Ss310 = ML2*(12*ML2 - 3*MZ2) - Ss309*MHiggs2(2)

	Ss311 = ML2*(12*ML2 - 3*MZ2) - Ss309*MHiggs2(1)

	Ss722 = MBot2*(196*MBot2 - Ss720) + MTop2*Ss721

	Ss726 = MBot2*(76*MBot2 - Ss724) + MTop2*Ss725

	Ss805 = MBot2*(20*MBot2 - 5*Ss145) + MZ2*MHiggs2(1)

	Ss508 = MBot2*(19*MBot2 - 2*Ss506) + MTop2*Ss507

	Ss144 = MBot2*(12*MBot2 - 3*Ss143) + MZ2*MHiggs2(2)

	Ss146 = MBot2*(12*MBot2 - 3*Ss145) + MZ2*MHiggs2(1)

	Ss817 = 4*MBot2**2 + MZ2*(MZ2 + 2*Ss58*Log(MZ2))

	Ss814 = 4*MBot2**2 + 
     &   MHiggs2(1)*(2*Ss60*Log(MHiggs2(1)) + MHiggs2(1))

	Ss37 = Ss35**2 + MBot2*(MBot2 - 2*Ss36)

	Ss24 = Ss22**2 + MBot2*(MBot2 - 2*Ss23)

	Ss996 = 3*CW*MTop*Ss978*UStop21 + Ss979*UStop22

	Ss980 = 3*CW*MTop*Ss978*UStop11 + Ss979*UStop12

	Ss957 = CW*ML*Ss937*UStau21 + Ss938*UStau22

	Ss947 = CW*ML*Ss924*UStau21 - Ss925*UStau22

	Ss939 = CW*ML*Ss937*UStau11 + Ss938*UStau12

	Ss926 = CW*ML*Ss924*UStau11 - Ss925*UStau12

	Ss914 = 3*CW*MBot*Ss894*USbot21 + Ss895*USbot22

	Ss904 = 3*CW*MBot*Ss881*USbot21 - Ss882*USbot22

	Ss896 = 3*CW*MBot*Ss894*USbot11 + Ss895*USbot12

	Ss883 = 3*CW*MBot*Ss881*USbot11 - Ss882*USbot12

	Ss99 = 3*MTop2*Ss98 + 220*MW2*TB2

	Ss168 = 8*MW2**2*Ss96 + 9*MBot2*Ss167*TB2

	Ss156 = 8*MW2**2*Ss96 - 18*MBot2*Ss155*TB2

	Ss173 = 4*MW2**2*Ss96 - 9*MBot2*Ss172*TB2

	Ss164 = 4*MW2**2*Ss96 + 9*MBot2*Ss162*Ss163*TB2

	Ss151 = 8*MW2*Ss96 + 27*MBot2*TB2

	Ss149 = 8*MTop2*Ss96 - Ss148*TB2

	Ss90 = MTop2*Ss81 - 64*MW2*TB2

	Ss806 = Ss805 - MBot2*Ss61*Log(MBot2)

	Ss704 = 2*MBot2*Ss699 - Ss700

	Ss698 = 6*ML2*Ss36*Ss649*Ss650*Ss692 + MZ2*Ss697

	Ss690 = ML2*Ss683 + 3*MTop2*MZ2**2*Ss688

	Ss685 = 9*MTop2*MZ2**2*Ss682 - 3*ML2*Ss683

	Ss662 = 5*MTop2*Ss635 + MW2*Ss661

	Ss578 = 1/12.D0*Ss577 + TF(MBot2,MTop2,MW2)

	Ss568 = 1/12.D0*Ss567 + TF(MBot2,MHiggs2(4),MTop2)

	Ss569 = MTop2*Ss555 - 4*MW2*TB2

	Ss556 = MTop2*Ss555 - 4*TB2*MHiggs2(4)

	Ss792 = Ss546 + Log(MTop2)

	Ss783 = Ss546 + Log(ML2)

	Ss534 = 1/12.D0*Ss533 + TF(0.D0,ML2,MW2)

	Ss529 = 1/12.D0*Ss528 - TF(0.D0,MHiggs2(4),ML2)

	Ss520 = MTop2*Ss519 - TB2*MHiggs2(4)

	Ss515 = 1/12.D0*Ss514 + TF(0.D0,MBot2,MTop2)

	Ss848 = Ss484 + 6*Ss847

	Ss834 = Ss484 - 12*Ss833

	Ss803 = Ss484 + 6*Ss801

	Ss718 = Ss484 + 6*Ss716

	Ss811 = Ss810 + 1/6.D0*Ss484

	Ss630 = Ss482 - 6*Ss628

	Ss626 = Ss481 - 6*Ss625

	Ss553 = Ss479 - 6*Ss551

	Ss822 = Ss335*Ss474 + 1/12.D0*Ss484

	Ss470 = Ss468 - 4*MW2*Ss469

	Ss471 = Ss467 - MW2*Ss470

	Ss452 = CB*Ss25**2*Ss451 + 18*ML2**2*MW2**2*MZ2*SB*Ss115*SW2*TB

	Ss445 = MTop2*Ss433 + 4*MW2*MZ2*SW2

	Ss764 = Ss335*Ss398 + 1/12.D0*Ss455

	Ss755 = Ss141*Ss398 + 1/12.D0*Ss455

	Ss823 = Ss375 + 12*Ss822

	Ss765 = Ss375 + 12*Ss764

	Ss544 = Ss368 - 6*Ss543

	Ss342 = Ss340*TB2 + Ss341*MHiggs2(4)

	Ss865 = 6*Ss335**2 + Ss375

	Ss743 = Ss141*Ss335 + 1/12.D0*Ss375

	Ss477 = MTop2*Ss335 + MZ2*Ss474

	Ss408 = Ss335 + Ss398

	Ss731 = 72*Ss325*Ss327 + 49*Ss730

	Ss734 = Ss733 + 1/12.D0*(MW2*Ss324)

	Ss541 = Ss317 - 6*Ss540

	Ss71 = -(MW2*Ss70) + Ss31*MHiggs2(4)

	Ss838 = -3*ML2*Ss8 + Ss309*MHiggs2(1)

	Ss825 = 6*Ss306**2 + Ss324

	Ss476 = ML2*Ss306 + MZ2*Ss474

	Ss302 = MBot2*Ss301*TB2 + 4*Ss298*Ss299*MHiggs2(4)**2

	Ss287 = MTop2*Ss284 - 16*MW2*TB2

	Ss411 = -(MW2*Ss410) + Ss279*MHiggs2(4)

	Ss365 = 2*MTop2*Ss335 + Ss269*MHiggs2(2)

	Ss363 = MTop2*Ss335 + Ss269*MHiggs2(2)

	Ss316 = 2*ML2*Ss306 + Ss269*MHiggs2(2)

	Ss314 = ML2*Ss306 + Ss269*MHiggs2(2)

	Ss857 = Ss268 + Ss375

	Ss843 = Ss268 + Ss324

	Ss796 = Ss795 + 1/6.D0*Ss268

	Ss787 = Ss786 + 1/6.D0*Ss268

	Ss999 = Ss263*MHiggs2(1) + Ss988*MStop2(2)

	Ss983 = Ss263*MHiggs2(1) + Ss969*MStop2(1)

	Ss960 = Ss263*MHiggs2(1) + Ss949*MStau2(2)

	Ss942 = Ss263*MHiggs2(1) + Ss928*MStau2(1)

	Ss917 = Ss263*MHiggs2(1) + Ss906*MSbot2(2)

	Ss899 = Ss263*MHiggs2(1) + Ss885*MSbot2(1)

	Ss364 = 2*MTop2*Ss335 + Ss263*MHiggs2(1)

	Ss362 = MTop2*Ss335 + Ss263*MHiggs2(1)

	Ss315 = 2*ML2*Ss306 + Ss263*MHiggs2(1)

	Ss313 = ML2*Ss306 + Ss263*MHiggs2(1)

	Ss262 = MTop2*Ss241**2*Ss242**2*Ss254 + MBot2*Ss261*TB2

	Ss307 = 3*ML2*Ss245*Ss306 + 4*Ss25**2*MHiggs2(4)

	Ss250 = 3*MBot2*Ss245 - 2*Ss249

	Ss236 = Ss233 - MW2*Ss235

	Ss228 = MBot2*Ss226 - Ss227

	Ss224 = 4*MBot2*Ss221 - Ss223

	Ss237 = 6*Ss215**2*Ss216*Ss217*Ss218**2*Ss219*Ss220 - MW2*Ss236

	Ss756 = Ss214 + 12*Ss755

	Ss744 = Ss214 + 12*Ss743

	Ss538 = Ss212 - 6*Ss536

	Ss204 = 36*MTop2*Ss200 - MW2*Ss203

	Ss344 = 27*MTop2*Ss199 + 2*MW2*Ss343

	Ss196 = MW2*Ss195 - 72*MTop2*Ss96

	Ss518 = 3*MTop2*Ss185 - Ss517*MHiggs2(4)

	Ss177 = MTop2*Ss176 - 2*Ss115*TB2

	Ss380 = 4*MW2**2*Ss162*Ss163 + MBot2*Ss379*TB2

	Ss271 = 2*MBot2*Ss141 + Ss269*MHiggs2(2)

	Ss265 = 2*MBot2*Ss141 + Ss263*MHiggs2(1)

	Ss475 = MBot2*Ss141 + MZ2*Ss474

	Ss270 = MBot2*Ss141 + Ss269*MHiggs2(2)

	Ss264 = MBot2*Ss141 + Ss263*MHiggs2(1)

	Ss399 = Ss141 + Ss398

	Ss336 = Ss141 + Ss335

	Ss361 = MTop2*Ss335 + Ss138*MHiggs2(3)

	Ss312 = ML2*Ss306 + Ss138*MHiggs2(3)

	Ss210 = MBot2*Ss141 + Ss138*MHiggs2(3)

	Ss134 = MW2*Ss132 - 5*Ss133

	Ss965 = 6*CW*MTop2*SA + CAB*MW*MZ*SB*Ss211

	Ss936 = 2*CW*ML2*SA - CB*MW*MZ*SAB*Ss922

	Ss893 = 6*CW*MBot2*SA - CB*MW*MZ*SAB*Ss879

	Ss436 = 9*MW2**4*MZ2 - 3*MBot2**2*Ss25**2*Ss433

	Ss442 = 3*MW2**2*MZ2 - MTop2**2*Ss433

	Ss684 = 12*MTop2**2*MZ2 + ML2**2*Ss649*Ss650

	Ss689 = 4*MTop2**2*MZ2 - ML2**2*Ss649*Ss650

	Ss572 = 4*MW2**3 - MBot2*Ss571

	Ss412 = MTop2*MW2**2 - Ss411*MHiggs2(4)

	Ss656 = 20*MW2**2 - Ss649*Ss650*Ss655*myB0(MW2,0.D0,ML2,1.D0)

	Ss639 = 2*MW2 + Ss637*Ss638*Log(Abs(Ss637*Ss638))

	Ss124 = 3*CB2*MTop2**2 - SB2*Ss123

	Ss644 = 4*MBot2*MTop2 - Ss643**2

	Ss487 = 4*MBot2*MTop2 - Ss486**2

	Ss977 = 6*CA*CW*MTop2 - MW*MZ*SAB*SB*Ss211

	Ss923 = 2*CA*CW*ML2 - CAB*CB*MW*MZ*Ss922

	Ss867 = 2*MZ2*Ss865 + 69*MHiggs2(1)

	Ss855 = -(Ss64*Log(MTop2)) + MHiggs2(1)

	Ss841 = -(Ss54*Log(ML2)) + MHiggs2(1)

	Ss880 = 6*CA*CW*MBot2 - CAB*CB*MW*MZ*Ss879

	Ss552 = 23 + 3*Ss484

	Ss537 = 23 + 3*Ss140

	Ss791 = 5 - Ss790*Log(MTop2)

	Ss782 = 5 - Ss781*Log(ML2)

	Ss710 = 5 - Ss709*Log(MBot2)

	Ss632 = 1 - atan(Re(1/Sqrt(-Ss126)))*Sqrt(-Ss126)

	Ss505 = 1 + Sqrt(Ss501)

	Ss504 = 1 - Sqrt(Ss501)

	Ss500 = 1 + Sqrt(Ss496)

	Ss499 = 1 - Sqrt(Ss496)

	Ss856 = MTop2*(28*MTop2 + 16*Ss619) - 4*Ss854*Log(MTop2) + 
     &   MHiggs2(1)*(32*MTop2 - 8*Ss619 - 8*Ss855*Log(MHiggs2(1)) + 
     &      5*MHiggs2(1))

	Ss842 = ML2*(28*ML2 + 16*Ss617) - 4*Ss840*Log(ML2) + 
     &   MHiggs2(1)*(32*ML2 - 8*Ss617 - 8*Ss841*Log(MHiggs2(1)) + 
     &      5*MHiggs2(1))

	Ss793 = Ss19*Ss509 + 
     &   MTop2*(-2*Ss375 - 24*Ss511*Log(MTop2)) + 
     &   (48 - Ss268 + 12*Log(MTop2)**2 - 6*Log(MHiggs2(1)/MTop2)**2)*
     &    MHiggs2(1)

	Ss784 = Ss16*Ss509 + 
     &   ML2*(-2*Ss324 - 24*Ss326*Log(ML2)) + 
     &   (48 - Ss268 + 12*Log(ML2)**2 - 6*Log(MHiggs2(1)/ML2)**2)*
     &    MHiggs2(1)

	Ss251 = 3*MTop2**2*MW2**2 + MBot2*Ss250 + 
     &   MHiggs2(4)*(-2*MTop2*MW2*Ss243 + 
     &      MHiggs2(4)*(Ss246 + 
     &         MHiggs2(4)*(-2*Ss244 + 4*MHiggs2(4))))

	Ss358 = MTop2*(18*MTop2*MW2**2 + 27*MW2**3*Ss199) - 
     &   MW2**4*Ss343 + MHiggs2(4)*(2*MW2*Ss357 - Ss356*MHiggs2(4))

	Ss382 = -(Ss241**2*Ss242**2*Ss380) + 
     &   MTop2*(11*MTop2*MW2**3 - Ss381) + 
     &   MTop2**2*MW2*MHiggs2(4)*(-16*MW2 + 11*MHiggs2(4))

	Ss305 = -3*MTop2*Ss245*Ss30**2 - Ss241**2*Ss242**2*Ss302 + 
     &   MHiggs2(4)*(MTop2*Ss303 + 
     &      MTop2*MHiggs2(4)*(-Ss304 + 17*MTop2*MHiggs2(4)))

	Ss225 = MTop2*(12*MTop2**3 + MBot2*Ss222) + MBot2**2*Ss224 + 
     &   MTop2**2*MHiggs2(4)*(-29*MTop2 + 15*MHiggs2(4))

	Ss446 = MTop2*MW2**2*Ss433 + 4*MW2**3*MZ2*SW2 + 
     &   MHiggs2(4)*(-2*MW2*Ss445 + Ss445*MHiggs2(4))

	Ss714 = 2*MBot2*Ss214 - 12*Ss547 + Ss268*MHiggs2(1) + 
     &   6*Ss713*MHiggs2(1)

	Ss657 = 34*MBot2*MW2**2 - 40*MW2**3 + ML2*(2*ML2*MW2 + Ss656)

	Ss229 = 6*MTop2**3 + MBot2*Ss228 + 
     &   MTop2*MHiggs2(4)*(-15*MTop2 + 7*MHiggs2(4))

	Ss839 = ML2*(92*ML2 - 23*MZ2) + 4*Ss838*Log(ML2) - 
     &   Ss837*MHiggs2(1)

	Ss968 = 6*CW*MTop*Ss966*UStop21*UStop22 + 
     &   Ss965*UStop22**2 + Ss967*(2 - 2*UStop22**2)

	Ss794 = -(1/12.D0*Ss793) + 
     &   Ss792*Log(MHiggs2(1))*MHiggs2(1) - TF(MHiggs2(1),MTop2,MTop2)

	Ss785 = -(1/12.D0*Ss784) + 
     &   Ss783*Log(MHiggs2(1))*MHiggs2(1) - TF(MHiggs2(1),ML2,ML2)

	Ss373 = 5*Ss124*TB + 2*ML2**2*SB2*TB*Log(ML2) - 
     &   6*CB*MTop2**2*SB*Log(MTop2)

	Ss715 = MBot2*(4*MBot2*Ss710 + 1/3.D0*Ss712) - 1/6.D0*(Ss714*MHiggs2(1))

	Ss797 = -2*Ss619 + 1/3.D0*(MTop2*Ss375) + Ss796*MHiggs2(1)

	Ss473 = Ss375 + 12*Ss335*Ss398 + Ss455

	Ss849 = -2*Ss625 + 1/3.D0*(ML2*Ss324) + 1/6.D0*(MZ2*Ss848)

	Ss788 = -2*Ss617 + 1/3.D0*(ML2*Ss324) + Ss787*MHiggs2(1)

	Ss812 = -2*Ss551 + MZ2*Ss811 + 1/3.D0*(MBot2*Ss214)

	Ss463 = Ss214 + 12*Ss141*Ss398 + Ss455

	Ss389 = Ss214 + 12*Ss141*Ss335 + Ss375

	Ss428 = 3*MBot2**2*SB2 + ML2**2*SB2 - Ss124

	Ss645 = Ss643 + Sqrt(-Ss644)

	Ss868 = 7 - 1/3.D0*(MTop2*Ss867)/(MZ2*MHiggs2(1)) - Ss866*Log(MTop2)

	Ss135 = MW2**2*(680*MZ2**3 + 
     &      MW2*(640*MW2*MZ2 - 1320*MZ2**2)) + 
     &   MTop2*(MW2*(-369*MZ2**3 + 
     &         MW2*(MW2*(768*MW2 - 960*MZ2) + 561*MZ2**2)) + 
     &      MTop2*(MZ2**3*(-90 - 270*SA2) + 
     &         MW2*MZ2**2*(144 + 270*SA2))) + 
     &   MBot2*(3*Ss131 + 18*MBot2*MZ2**2*Ss134) + 
     &   ML2*(3*MW2*Ss129 + 6*ML2*MZ2**2*Ss134)

	Ss440 = MW2**5*(6*MW2 - 6*MZ2) + 
     &   MBot2*(-(MBot2**2*Ss25**2*Ss433) + 3*Ss25**2*Ss434) + 
     &   MTop2*(-Ss436 + MTop2*(MTop2*MW2**2*Ss433 - 12*Ss439)) + 
     &   MHiggs2(4)*(2*MW2*Ss438 + Ss437*MHiggs2(4))

	Ss870 = MTop2*(28*MTop2 + 32*MZ2 + 16*Ss628) - 
     &   2/3.D0*(MZ2*Ss45*Ss484) - 4*Ss869*Log(MTop2) + 
     &   MZ2*(-8*Ss628 + MZ2*(5 - 8*Log(MZ2)))

	Ss818 = MBot2*(28*MBot2 + 32*MZ2 + 16*Ss551) - 
     &   2/3.D0*(MZ2*Ss484*Ss58) - 4*Ss817*Log(MBot2) + 
     &   MZ2*(-8*Ss551 + MZ2*(5 - 8*Log(MZ2)))

	Ss815 = -(2/3.D0*(Ss268*Ss60*MHiggs2(1))) - 
     &   4*Ss814*Log(MBot2) + 
     &   MBot2*(28*MBot2 + 16*Ss547 + 32*MHiggs2(1)) + 
     &   MHiggs2(1)*(-8*Ss547 + (5 - 8*Log(MHiggs2(1)))*MHiggs2(1))

	Ss120 = -(MTop2**3*Ss113) + Ss115*Ss119 - 
     &   MTop2*MW2*Ss112*MHiggs2(4) + 
     &   MBot2*(MBot2*(MTop2*Ss109 + 7*MBot2*Ss111) - MTop2**2*Ss114 + 
     &      MW2*Ss110*MHiggs2(4))

	Ss807 = 2/3.D0*(MBot2*Ss214*Ss61) - 4*Ss806*Log(MBot2) + 
     &   MBot2*(108*MBot2 - 27*MZ2 - 27*MHiggs2(1)) + 7*MZ2*MHiggs2(1)

	Ss642 = MW2**2*(10*MW2 - 4*MZ2) - 3*MTop2*MW2*Ss443 - 
     &   Ss635*Ss640 - MBot2*Ss641

	Ss615 = -(MW2*Ss596) + 6*Ss432*Ss605 + 144*Ss578*Ss614 + 
     &   1/2.D0*(Ss34*Ss587)

	Ss488 = Ss486 + Sqrt(-Ss487)

	Ss800 = -2*Ss625 + 2*ML2*Ss782 + 1/3.D0*(ML2*Ss324) + 
     &   1/6.D0*(MZ2*Ss484) + MZ2*(5 - Ss799*Log(ML2) - 4*Log(MZ2))

	Ss862 = 1/6.D0*(MTop2*Ss375) - 1/12.D0*(Ss12*Ss509) + 
     &   1/2.D0*(MZ2*Log(MTop2/MZ2)**2) + 
     &   MTop2*(Ss791 + 2*Ss511*Log(MTop2)) + MZ2*Ss550*Log(MZ2) - 
     &   TF(MTop2,MTop2,MZ2)

	Ss774 = MBot2*(MTop2*(Ss334*Ss735 - 2*Ss742*Ss744) - 
     &      Ss397*Ss745 - 2*Ss754*Ss756) + 
     &   MTop2*(Ss407*Ss746 + 2*Ss763*Ss765)

	Ss680 = -(Ss32**2*Ss33**2*Ss678) + 
     &   MBot2*(-Ss676 + MBot2*
     &       (MBot2*(MBot2*(MBot2*Ss635 - Ss673) - Ss674) + Ss679))

	Ss672 = -(Ss32**2*Ss33**2*Ss669) + 
     &   MBot2*(MBot2*(MBot2*(MBot2*(MBot2*Ss635 - Ss662) + Ss665) + 
     &         Ss667) + Ss671)

	Ss495 = -(Ss35**2*Ss492) + 
     &   MBot2*(MBot2*(MBot2*(MBot2*Ss433 - Ss244*Ss433) - Ss493) + Ss494)

	Ss472 = 19*Ss215**7*Ss218**7 + 
     &   MW2*(-3*Ss215**5*Ss218**5*Ss464 + 
     &      MW2*(2*Ss215**3*Ss218**3*Ss465 + 
     &         MW2*(-2*Ss215*Ss218*Ss466 + MW2*Ss471)))

	Ss462 = 19*Ss215**7*Ss218**7 + 
     &   MW2*(-(Ss215**5*Ss218**5*Ss456) + 
     &      MW2*(2*Ss215**3*Ss218**3*Ss457 + 
     &         MW2*(-2*Ss215*Ss218*Ss458 + MW2*Ss461)))

	Ss691 = -18*MZ2**2*Ss689 + 
     &   MW2*(MW2*(MW2*(16*MW2*MZ2 - 216*MZ2**2) + MZ2**2*Ss687) + 
     &      3*Ss690)

	Ss686 = 18*MZ2**2*Ss684 + 
     &   MW2*(MW2*(MW2*(16*MW2*MZ2 - 252*MZ2**2) - MZ2**2*Ss681) - 
     &      Ss685)

	Ss447 = -(Ss25**2*Ss444) + 
     &   MBot2*(-3*Ss25**2*Ss442 + MBot2*(MBot2*Ss25**2*Ss433 - 3*Ss446))

	Ss422 = Ss421 + MBot2*
     &    (Ss416 + MBot2*(9*MBot2*Ss409*TB2 + 18*Ss412*TB2))

	Ss283 = 45*MTop2*MW2**2*TB2 + 
     &   MHiggs2(4)*(Ss282 + 
     &      MHiggs2(4)*(Ss281*TB2 + 17*TB2*MHiggs2(4)))

	Ss987 = Ss965*UStop21**2 + 
     &   2*UStop22*(-3*CW*MTop*Ss966*UStop21 + 
     &      Ss967*UStop22)

	Ss93 = -Ss92 + MHiggs2(4)*(Ss90 + 32*TB2*MHiggs2(4))

	Ss94 = 15*MTop2**2*MW2**2*Ss89 + 
     &   MHiggs2(4)*(-(MTop2*Ss91) - Ss93*MHiggs2(4))

	Ss798 = MTop2*(4*MTop2*Ss791 + 4*Ss794) - Ss797*MHiggs2(1)

	Ss789 = ML2*(4*ML2*Ss782 + 4*Ss785) - Ss788*MHiggs2(1)

	Ss850 = -(MZ2**2*Ss846) + ML2*(2*ML2*Ss782 + Ss849)

	Ss813 = -(MZ2**2*Ss808) + MBot2*(2*MBot2*Ss710 + Ss812)

	Ss573 = Ss570 + MBot2*(Ss569 + 3*MBot2*TB2)

	Ss562 = Ss557 + MBot2*(Ss556 + 3*MBot2*TB2)

	Ss346 = -Ss345 + MHiggs2(4)*(-Ss344 + Ss343*MHiggs2(4))

	Ss290 = -Ss289 + MHiggs2(4)*(Ss287 + 8*TB2*MHiggs2(4))

	Ss291 = 9*MTop2**2*MW2**2*Ss285 + 
     &   MHiggs2(4)*(-(MTop2*Ss286) - Ss290*MHiggs2(4))

	Ss292 = Ss291 + MBot2*(Ss283 - 9*MBot2*Ss280*TB2)

	Ss192 = Ss191*MHiggs2(4)**2 + 
     &   MTop2*(36*MTop2*Ss186 + Ss189*MHiggs2(4))

	Ss183 = 36*MTop2**2*Ss32**2*Ss33**2 + 
     &   MHiggs2(4)*(-(MTop2*Ss182) + Ss180*MHiggs2(4))

	Ss174 = MTop2*(9*MTop2*Ss171 - 2*Ss173) - Ss170*TB2

	Ss165 = MTop2*(9*MTop2*Ss161 - 2*Ss164) + Ss160*TB2

	Ss169 = MTop2*(9*MTop2*Ss154 - 2*Ss168) + Ss162*Ss163*Ss166*TB2

	Ss157 = MTop2*(9*MTop2*Ss154 + Ss156) + MBot2*Ss153*TB2

	Ss108 = MW2**2*Ss104 + 
     &   MHiggs2(4)*(-(MW2*Ss106) - Ss107*MHiggs2(4))

	Ss152 = MTop2*(9*MTop2 + Ss151) - 2*MW2*Ss150*TB2

	Ss997 = Ss977*UStop21**2 + 2*Ss996*UStop22

	Ss981 = 2*Ss980*UStop12 + Ss977*UStop22**2

	Ss958 = Ss936*UStau21**2 + 2*Ss957*UStau22

	Ss948 = Ss923*UStau21**2 - 2*Ss947*UStau22

	Ss940 = 2*Ss939*UStau12 + Ss936*UStau22**2

	Ss927 = -2*Ss926*UStau12 + Ss923*UStau22**2

	Ss915 = Ss893*USbot21**2 + 2*Ss914*USbot22

	Ss905 = Ss880*USbot21**2 - 2*Ss904*USbot22

	Ss897 = 2*Ss896*USbot12 + Ss893*USbot22**2

	Ss884 = -2*Ss883*USbot12 + Ss880*USbot22**2

	Ss100 = Ss99 - 83*TB2*MHiggs2(4)

	Ss101 = Ss97 - Ss100*MHiggs2(4)

	Ss102 = 3*MTop2*MW2**2*Ss95 - Ss101*MHiggs2(4)

	Ss858 = 3*Ss856 + 2*Ss64*Ss857*MHiggs2(1)

	Ss844 = 3*Ss842 + 2*Ss54*Ss843*MHiggs2(1)

	Ss845 = ML2*((4*ML2*Ss839)/(MZ2*Ss8) + 1/3.D0*Ss844/Ss44)

	Ss835 = Ss832 + 4*ML2*Ss834

	Ss826 = MZ2*Ss825 + Ss15*Ss825

	Ss778 = Ss729 - 288*ML2**2*Ss325*Ss734

	Ss705 = MBot2**3*Ss704 + 2*Ss649*Ss650*Ss703*myB0(MW2,MBot2,MTop2,1.D0)

	Ss658 = MZ2*Ss654 + Ss657

	Ss659 = 8*MW2**3*Ss653 + MZ2*Ss658

	Ss827 = -(ML2**2*Ss44*Ss826) + 3*Ss15*Ss52*Ss625*MHiggs2(1)

	Ss776 = Ss432*Ss605 + 24*Ss578*Ss773

	Ss816 = -Ss815 + 2/3.D0*(Ss214*Ss60*MHiggs2(1))

	Ss775 = Ss596 - 6*Ss774

	Ss819 = -Ss818 + 2/3.D0*(MZ2*Ss214*Ss58)

	Ss631 = MZ2*Ss552 + 6*Ss630

	Ss627 = MZ2*Ss552 + 6*Ss626

	Ss554 = MZ2*Ss552 + 6*Ss553

	Ss545 = 6*Ss544 + Ss537*MHiggs2(3)

	Ss542 = 6*Ss541 + Ss537*MHiggs2(3)

	Ss539 = 6*Ss538 + Ss537*MHiggs2(3)

	Ss521 = Ss520 + 3*MBot2*TB2

	Ss732 = 864*Ss508*Ss515 + ML2**2*Ss731

	Ss646 = ((2*MBot2*MTop)/Ss488 - 1/2.D0*Ss645/MTop)/MBot

	Ss489 = ((2*MBot2*MTop)/Ss488 - 1/2.D0*Ss488/MTop)/MBot

	Ss453 = Ss25*Ss452 + 
     &   36*ML2**2*MW2**3*MZ2*SB*SW2*TB*Log(MW2)*MHiggs2(4)

	Ss871 = -Ss870 + 2/3.D0*(MZ2*Ss375*Ss45)

	Ss448 = CB*Ss447 - 36*MBot2**2*MW2**2*MZ2*SB*Ss115*SW2*TB

	Ss441 = Ss25*Ss440 + 
     &   72*MTop2**2*MW2**3*MZ2*SW2*Log(MW2)*MHiggs2(4)

	Ss802 = Ss375 + 6*Ss791

	Ss359 = 9*MBot2*Ss355 - Ss358

	Ss423 = MTop2*Ss32*Ss33*Ss352 + MBot2*Ss422

	Ss360 = Ss32*Ss33*Ss352 + MBot2*Ss359

	Ss347 = 9*MTop2*MW2**2*Ss288 - Ss346*MHiggs2(4)

	Ss293 = MTop2*Ss22*Ss278 - MBot2*Ss292

	Ss425 = CB*SB*Ss121*Ss25**3 - 3*MW2*Ss124*Ss409*TB

	Ss125 = 5*CB*SB*Ss121*Ss25**2 - 3*MW2*Ss122*Ss124*TB

	Ss252 = MTop2*Ss240*Ss241**2*Ss242**2 + MBot2*Ss251*TB2

	Ss717 = Ss214 + 6*Ss710

	Ss207 = Ss206 + Ss204*MHiggs2(4)

	Ss208 = 36*MTop2*Ss202 + Ss207*MHiggs2(4)

	Ss197 = Ss196 + 36*TB2*MHiggs2(4)

	Ss198 = 36*Ss194 - Ss197*MHiggs2(4)

	Ss136 = CA2*Ss125 - 3*MW2*SA2*Ss122*Ss124*TB

	Ss647 = 2*MW2 - MBot*MTop*Ss646*Log(1/2.D0*Abs(Ss645/(MBot*MTop)))

	Ss648 = 2 - Sqrt(Ss501)*Log(Abs(Ss505/Ss504))

	Ss563 = MTop2*(-3*MTop2**3 - MBot2*Ss561) + MBot2**2*Ss562 + 
     &   MHiggs2(4)*(11*MTop2**3 + 
     &      MTop2*MHiggs2(4)*(-12*MTop2 + 4*MHiggs2(4)))

	Ss706 = Ss32**2*Ss33**2*Ss691 + 
     &   MBot2*(MBot2*Ss686 - Ss698) + 9*MZ2**2*Ss705

	Ss238 = 6*MTop2*Ss215**2*Ss216*Ss217*Ss218**2*Ss219*
     &    Ss220 + MW2*(-Ss225 + MW2*Ss229) + MBot2*Ss237*TB2

	Ss522 = MTop2*(-3*MTop2**2 + MBot2*Ss518) - MBot2**2*Ss521 + 
     &   8*MTop2**2*MHiggs2(4)

	Ss804 = 2*Ss628 - 1/3.D0*(MTop2*Ss802) - 1/6.D0*(MZ2*Ss803)

	Ss777 = -(MW2*Ss775) + 6*Ss776 + 1/2.D0*(Ss34*Ss587)

	Ss863 = (-2*Ss862 + MZ2*(-Ss861 + 1/6.D0*Ss375))/MZ2

	Ss137 = -18*MZ2**3*Ss136*SW2 + 
     &   CB*SB*Ss25**2*(Ss135 - 
     &      12*MW2**2*Ss128*atan(Re(1/Sqrt(-Ss126)))*
     &       Sqrt(-Ss126))

	Ss660 = 32*MBot2*MW2**4 - 48*ML2*MW2**4 - MZ2*Ss659

	Ss874 = MTop2*(12*MTop2**2*Ss44*Ss853 + MZ2*Ss858)

	Ss209 = MTop2*Ss183 + 
     &   MBot2*(-Ss192 + MBot2*
     &       (Ss208 + MBot2*(-Ss198 + MBot2*(36*Ss177 + 36*MBot2*TB2))))

	Ss574 = MTop2*(MTop2*(MTop2*(3*MTop2 - 11*MW2) + 12*MW2**2) - 
     &      Ss572) - MBot2**2*Ss573

	Ss158 = Ss157 + MHiggs2(4)*
     &    (-2*Ss152 + MHiggs2(4)*(Ss149 + 8*TB2*MHiggs2(4)))

	Ss707 = MW2*Ss706 + 
     &   MZ2**2*(-9*Ss672*Log(MBot2) + 9*Ss680*Log(MTop2))

	Ss719 = 4*Ss551 - 2/3.D0*(MBot2*Ss717) - 1/3.D0*(MZ2*Ss718)

	Ss449 = Ss25*Ss448 + 
     &   MBot2**2*MW2**3*MZ2*SB*SW2*TB*
     &    (-72*Log(MW2) + 72*Log(MHiggs2(4)))*MHiggs2(4)

	Ss836 = 12*ML2**2*Ss829 + MZ2*(2*MZ2*Ss324 - Ss835)

	Ss175 = MW2**2*Ss174 + 
     &   MHiggs2(4)*(-(MW2*Ss169) + Ss165*MHiggs2(4))

	Ss872 = (4*MTop2*Ss868)/Ss17 - Ss871/(MZ2*Ss44)

	Ss873 = MTop2*(-Ss872 + 8/3.D0*(MTop2**2*Ss865)/(MZ2*MHiggs2(1)))

	Ss864 = MZ2*Ss859 + MTop2*Ss863

	Ss851 = Ss8*Ss836 + 24*Ss44*Ss850

	Ss828 = -((MTop2**2*MZ2*Ss63*Ss823)/Ss9**2) + 
     &   (3*Ss824**2*Ss827)/(Ss15*Ss44*Ss8*MHiggs2(1))

	Ss820 = (-(Ss819/Ss7) + (MZ2*Ss816)/(Ss13*MHiggs2(1)))/MZ2

	Ss523 = -12*MW2*Ss522 + 36*Ss215**2*Ss218**2*Ss516*MHiggs2(4)

	Ss426 = CA2*Ss425 - 3*MW2*SA2*Ss124*Ss409*TB

	Ss427 = CB*SB*Ss25**3*Ss424 + 6*Ss426

	Ss348 = 9*MBot2*Ss342 - Ss347

	Ss349 = Ss26*Ss27*Ss278 - MBot2*Ss348

	Ss159 = 18*MW2**2*Ss147 - Ss158*MHiggs2(4)

	Ss779 = 3*Ss29*Ss37**3*Ss722 + 
     &   Ss37**3*(12*ML2**2*Ss320*Ss723 + 36*Ss371*Ss726 + 
     &      Ss732) - 6*Ss777

	Ss821 = ((-8*Ss813)/Ss7**2 - (MZ2*Ss820)/Ss44)/MZ2 + 
     &   (4*MBot2*Ss807)/(MZ2*Ss13*Ss7*MHiggs2(1))

	Ss429 = Ss427 + MW2**2*(108*CB*MTop2**2*SB - 36*Ss428*TB)*
     &    MHiggs2(4)

	Ss708 = (6*Ss660)/MW2**2 + Ss707/(MW2**3*Ss37) - 
     &   (6*MZ2**2*Ss652*Log(ML2))/MW2**3 - 18*MZ2**3*Log(MZ2) - 
     &   6*Ss498*Sqrt(Ss496)*Log(Abs(Ss500/Ss499))

	Ss875 = 8*MTop2*Ss864 - Ss873*Ss9

	Ss852 = 1/3.D0*(ML2*Ss851)/(MZ2*Ss44*Ss8**2) + 
     &   Ss845/(Ss15*MHiggs2(1))

	Ss780 = Ss37**3*Ss778 - Ss28*Ss779

	Ss876 = ((MZ2*Ss875)/Ss9**2 + 
     &     1/3.D0*Ss874/(Ss17*Ss44*MHiggs2(1)))/MZ2

	Ss877 = 3*MBot2*MM2*MZ2*Ss56*Ss821 - 9*MZ2*Ss53**2*Ss852 + 
     &   MM2*(-8*Ss828 - 6*MZ2*Ss63*Ss876)

	end

