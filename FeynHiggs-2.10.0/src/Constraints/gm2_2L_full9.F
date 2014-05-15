* gm2_2L_full9.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full9
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp73, tmp74, tmp75, tmp76, tmp77, tmp78, tmp79
	ComplexType tmp80, tmp81, tmp82, tmp83, tmp84, tmp85, tmp86
	ComplexType tmp87, tmp88, tmp89, tmp90, tmp91, tmp92, tmp93
	ComplexType tmp94, tmp95, tmp96, tmp97

	Sf6005 = 3*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2*Sf3260 - Sf5999) + Sf6001) - 
     &         Sf6004) + Sf5943*MSq2Diff(2,2)**2)

	Sf5585 = (MW2 - MSq2Diff(1,2))**2*
     &    (MW2**2*(-45*MBot2 + 45*MTop2 + MW2 + 45*C2B*MW2) + 
     &      MW2*(9*MBot2 - 9*MTop2 + 70*MW2 - 9*C2B*MW2)*
     &       MSq2Diff(1,2) - 38*MW2*MSq2Diff(1,2)**2 + 
     &      3*MSq2Diff(1,2)**3) + 
     &   2*MW2*(MW2**2*(99*MBot2 - 99*MTop2 + 2*MW2 - 99*C2B*MW2) + 
     &      MW2*(-126*MBot2 + 126*MTop2 - 181*MW2 + 126*C2B*MW2)*
     &       MSq2Diff(1,2) + 
     &      (9*MBot2 - 9*MTop2 + 236*MW2 - 9*C2B*MW2)*
     &       MSq2Diff(1,2)**2 - 39*MSq2Diff(1,2)**3)*MStop2(2) + 
     &   4*MW2**2*(-63*MBot2 + 63*MTop2 - 8*MW2 + 63*C2B*MW2 + 
     &      111*MSq2Diff(1,2))*MStop2(2)**2

	Sf3602 = 3*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2*Sf3260 - Sf3594) + Sf3598) - 
     &         Sf3601) + Sf3509*MSq2Diff(2,1)**2)

	Sf3275 = 3*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2*Sf3260 - Sf3263) + Sf3269) - 
     &         Sf3274) + Sf3143*MSq2Diff(1,1)**2)

	Sf3155 = 3*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(-Sf3154 + 
     &         MW2*(Sf3149 + MW2*(Sf3139 + MW2*Sf601))) + 
     &      Sf3143*MSq2Diff(1,1)**2)

	Sf1980 = -5*MSl2Diff(2,1)**6 + 
     &   MW2*(Sf1961*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(Sf1977 + MW2*(MW2*Sf1962 - Sf1979)) + 
     &         Sf1975*MSl2Diff(2,1)))

	Sf1920 = -5*MSl2Diff(1,1)**6 + 
     &   MW2*(Sf1899*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(Sf1917 + MW2*(MW2*Sf1900 - Sf1919)) + 
     &         Sf1915*MSl2Diff(1,1)))

	Sf1842 = MW2*(MW2*
     &       (MW2*(MW2*(MW2*Sf1047 - Sf1837) + Sf1839) - Sf1841) - 
     &      Sf1840*MSl2Diff(2,1)) - 
     &   4*Sf1836*MSl2Diff(2,1)**3*MStau2(2)

	Sf1063 = MW2*(MW2*
     &       (MW2*(MW2*(MW2*Sf1047 - Sf1051) + Sf1056) - Sf1062) - 
     &      Sf1059*MSl2Diff(1,1)) - 
     &   4*Sf1048*MSl2Diff(1,1)**3*MStau2(1)

	Sf6039 = 5*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2 - Sf6034) - 2*Sf6036) + 
     &         2*Sf6038) - Sf6035*MSq2Diff(2,2)**3)

	Sf3699 = 5*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2 - Sf3693) - 2*Sf3695) + 
     &         2*Sf3698) - Sf3694*MSq2Diff(2,1)**3)

	Sf3421 = 5*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2 - Sf3415) - 2*Sf3417) + 
     &         2*Sf3420) - Sf3416*MSq2Diff(1,1)**3)

	Sf5051 = 4*MW2**2*MSq2Diff(2,2)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf5050 + MHiggs2(4)*
     &          (Sf5046 - Sf4918*MHiggs2(4) + MHiggs2(4)**2)) - 
     &      MW2*Sf5047*MSq2Diff(2,2))

	Sf4925 = 4*MW2**2*MSq2Diff(1,2)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf4924 + MHiggs2(4)*
     &          (Sf4920 - Sf4918*MHiggs2(4) + MHiggs2(4)**2)) - 
     &      MW2*Sf4921*MSq2Diff(1,2))

	Sf2597 = 4*MW2**2*MSq2Diff(2,1)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf2596 + MHiggs2(4)*
     &          (Sf2592 - Sf2347*MHiggs2(4) + MHiggs2(4)**2)) - 
     &      MW2*Sf2593*MSq2Diff(2,1))

	Sf2359 = 4*MW2**2*MSq2Diff(1,1)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf2358 + MHiggs2(4)*
     &          (Sf2353 - Sf2347*MHiggs2(4) + MHiggs2(4)**2)) - 
     &      MW2*Sf2354*MSq2Diff(1,1))

	Sf3668 = -(MW2**3*Sf3667*MSbot2(2)) + 
     &   (63*MSbot2(2)**5 + 
     &      MW2*(-70*MSbot2(2)**4 + 
     &         MW2*(8*MW2**3 + 20*MSbot2(2)**3)))*USbot21

	Sf3388 = -(MW2**3*Sf3387*MSbot2(1)) + 
     &   (63*MSbot2(1)**5 + 
     &      MW2*(-70*MSbot2(1)**4 + 
     &         MW2*(8*MW2**3 + 20*MSbot2(1)**3)))*USbot11

	Sf5162 = -(Sf5148*Sf5152*Ss25**2) + 
     &   MW2*(Sf5161 + MW2*
     &       (-(Sf5148*Sf5159*Ss25**2) + 
     &         MW2*(12*MW2*Sf5148*Ss25**2*MSbot2(2) - 
     &            Sf5148*Sf5157*Ss25**2*MSbot2(2))))

	Sf5155 = -(Sf5148*Sf5152*Ss25**2) + 
     &   MW2*(Sf5154 + MW2*
     &       (-(Sf5148*Sf5151*Ss25**2) + 
     &         MW2*(-(Sf5148*Sf5150*Ss25**2) + 
     &            4*MW2*Sf5148*Ss25**2*MSbot2(2))))

	Sf2755 = -(Sf2740*Sf2745*Ss25**2) + 
     &   MW2*(Sf2754 + MW2*
     &       (-(Sf2740*Sf2752*Ss25**2) + 
     &         MW2*(12*MW2*Sf2740*Ss25**2*MSbot2(2) - 
     &            Sf2740*Sf2750*Ss25**2*MSbot2(2))))

	Sf2748 = -(Sf2740*Sf2745*Ss25**2) + 
     &   MW2*(Sf2747 + MW2*
     &       (-(Sf2740*Sf2744*Ss25**2) + 
     &         MW2*(-(Sf2740*Sf2743*Ss25**2) + 
     &            4*MW2*Sf2740*Ss25**2*MSbot2(2))))

	Sf3671 = MW2**3*Sf3670*MSbot2(2) + 
     &   (42*MSbot2(2)**5 + 
     &      MW2*(-76*MSbot2(2)**4 + 
     &         MW2*(-2*MW2**3 + 34*MSbot2(2)**3)))*USbot21

	Sf5044 = -(Sf4827*Sf5033*Ss25**2) + 
     &   MW2*(Sf5043 + MW2*
     &       (-(Sf4827*Sf5041*Ss25**2) + 
     &         MW2*(12*MW2*Sf4827*Ss25**2*MSbot2(1) - 
     &            Sf4827*Sf5039*Ss25**2*MSbot2(1))))

	Sf5036 = -(Sf4827*Sf5033*Ss25**2) + 
     &   MW2*(Sf5035 + MW2*
     &       (-(Sf4827*Sf5032*Ss25**2) + 
     &         MW2*(-(Sf4827*Sf5031*Ss25**2) + 
     &            4*MW2*Sf4827*Ss25**2*MSbot2(1))))

	Sf2579 = -(Sf2562*Sf2568*Ss25**2) + 
     &   MW2*(Sf2578 + MW2*
     &       (-(Sf2562*Sf2576*Ss25**2) + 
     &         MW2*(12*MW2*Sf2562*Ss25**2*MSbot2(1) - 
     &            Sf2562*Sf2574*Ss25**2*MSbot2(1))))

	Sf2571 = -(Sf2562*Sf2568*Ss25**2) + 
     &   MW2*(Sf2570 + MW2*
     &       (-(Sf2562*Sf2567*Ss25**2) + 
     &         MW2*(-(Sf2562*Sf2566*Ss25**2) + 
     &            4*MW2*Sf2562*Ss25**2*MSbot2(1))))

	Sf3391 = MW2**3*Sf3390*MSbot2(1) + 
     &   (42*MSbot2(1)**5 + 
     &      MW2*(-76*MSbot2(1)**4 + 
     &         MW2*(-2*MW2**3 + 34*MSbot2(1)**3)))*USbot11

	Sf6842 = -3*Sf6836*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf6445 + Sf6837) - Sf6841) + 
     &      Sf6839*MSq2Diff(2,2))

	Sf6724 = -3*Sf4507*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf6445 + Sf6715) - Sf6723) + 
     &      Sf6719*MSq2Diff(2,1))

	Sf6587 = -3*Sf6581*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf6445 + Sf6582) - Sf6586) + 
     &      Sf6584*MSq2Diff(1,2))

	Sf6457 = -3*Sf4400*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf6445 + Sf6447) - Sf6456) + 
     &      Sf6451*MSq2Diff(1,1))

	Sf6835 = -(Sf6829*MSq2Diff(2,2)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf6434 + Sf6830) - Sf6834) + 
     &      Sf6832*MSq2Diff(2,2))

	Sf6714 = -(Sf4496*MSq2Diff(2,1)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf6434 + Sf6705) - Sf6713) + 
     &      Sf6709*MSq2Diff(2,1))

	Sf6580 = -(Sf6574*MSq2Diff(1,2)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf6434 + Sf6575) - Sf6579) + 
     &      Sf6577*MSq2Diff(1,2))

	Sf6444 = -(Sf4387*MSq2Diff(1,1)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf6434 + Sf6435) - Sf6443) + 
     &      Sf6439*MSq2Diff(1,1))

	Sf867 = MW2*(MW2*(MW2*(MW2*Sf601 + Sf860) - Sf863) + 
     &      Sf866) + Sf455*MSl2Diff(1,1)**3

	Sf7029 = -3*Sf6836*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf4399 - Sf7024) + Sf7028) + 
     &      Sf7026*MSq2Diff(2,2))

	Sf7007 = -3*Sf6581*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf4399 - Sf7002) + Sf7006) + 
     &      Sf7004*MSq2Diff(1,2))

	Sf4516 = -3*Sf4507*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf4399 - Sf4508) + Sf4515) + 
     &      Sf4511*MSq2Diff(2,1))

	Sf4412 = -3*Sf4400*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf4399 - Sf4403) + Sf4411) + 
     &      Sf4407*MSq2Diff(1,1))

	Sf536 = MW2*(MW2*(MW2*(MW2*Sf438 - Sf524) - Sf530) + 
     &      Sf535) - Sf525*MSl2Diff(1,1)**2

	Sf3947 = -(Sf3946*MSbot2(2)**4) + 
     &   MW2*(Sf3944*MSbot2(2)**3 + 
     &      MW2*(3*Sf3666*MSbot2(2)**2 + 
     &         MW2*(MW2*Sf3641 + Sf3943*MSbot2(2))))

	Sf3830 = -(Sf3829*MSbot2(1)**4) + 
     &   MW2*(Sf3825*MSbot2(1)**3 + 
     &      MW2*(3*Sf3386*MSbot2(1)**2 + 
     &         MW2*(MW2*Sf3345 + Sf3822*MSbot2(1))))

	Sf5133 = 3*Sf5099*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2525 - Sf5128) + Sf5130) - 
     &      Sf5132*MSq2Diff(2,2))

	Sf5007 = 3*Sf4973*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2525 - Sf5002) + Sf5004) - 
     &      Sf5006*MSq2Diff(1,2))

	Sf2720 = 3*Sf2663*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2525 - Sf2711) + Sf2715) - 
     &      Sf2719*MSq2Diff(2,1))

	Sf2541 = 3*Sf2458*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2525 - Sf2528) + Sf2534) - 
     &      Sf2540*MSq2Diff(1,1))

	Sf5127 = Sf5099*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2513 - Sf5122) + Sf5126) - 
     &      Sf5124*MSq2Diff(2,2))

	Sf5001 = Sf4973*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2513 - Sf4996) + Sf5000) - 
     &      Sf4998*MSq2Diff(1,2))

	Sf2709 = Sf2663*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2513 - Sf2701) + Sf2708) - 
     &      Sf2705*MSq2Diff(2,1))

	Sf2524 = Sf2458*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2513 - Sf2515) + Sf2523) - 
     &      Sf2519*MSq2Diff(1,1))

	Sf5111 = 3*Sf5099*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2472 + Sf5106) - Sf5108) - 
     &      Sf5110*MSq2Diff(2,2))

	Sf4985 = 3*Sf4973*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2472 + Sf4980) - Sf4982) - 
     &      Sf4984*MSq2Diff(1,2))

	Sf2682 = 3*Sf2663*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2472 + Sf2673) - Sf2677) - 
     &      Sf2681*MSq2Diff(2,1))

	Sf2489 = 3*Sf2458*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2472 + Sf2475) - Sf2482) - 
     &      Sf2488*MSq2Diff(1,1))

	Sf5105 = Sf5099*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2457 + Sf5100) + Sf5104) - 
     &      Sf5102*MSq2Diff(2,2))

	Sf4979 = Sf4973*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2457 + Sf4974) + Sf4978) - 
     &      Sf4976*MSq2Diff(1,2))

	Sf2672 = Sf2663*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2457 + Sf2664) + Sf2671) - 
     &      Sf2668*MSq2Diff(2,1))

	Sf2471 = Sf2458*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf2457 + Sf2460) + Sf2470) - 
     &      Sf2465*MSq2Diff(1,1))

	Sf6342 = -(Sf6245*MSq2Diff(2,2)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf2371 - Sf6337) - Sf6339) + 
     &      Sf6341*MSq2Diff(2,2))

	Sf6319 = -(Sf5870*MSq2Diff(1,2)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf2371 - Sf6314) - Sf6316) + 
     &      Sf6318*MSq2Diff(1,2))

	Sf4198 = -(Sf4197*MSq2Diff(2,1)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf2371 - Sf4188) - Sf4192) + 
     &      Sf4196*MSq2Diff(2,1))

	Sf4118 = -(Sf4117*MSq2Diff(1,1)**3) + 
     &   MW2*(MW2*(MW2*(MW2*Sf2371 - Sf4105) - Sf4110) + 
     &      Sf4116*MSq2Diff(1,1))

	Sf1993 = -(Sf1836*MSl2Diff(2,1)**3) + 
     &   MW2*(MW2*(-Sf1992 + MW2*(2*Sf1988 - MW2*Sf410)) + 
     &      2*Sf1990*MSl2Diff(2,1))

	Sf1935 = -(Sf1048*MSl2Diff(1,1)**3) + 
     &   MW2*(MW2*(-Sf1934 + MW2*(2*Sf1929 - MW2*Sf410)) + 
     &      2*Sf1932*MSl2Diff(1,1))

	Sf1680 = MW2*(Sf1679 + 
     &      MW2*(-Sf1676 + MW2*(Sf1674 + MW2*Sf601))) + 
     &   Sf1372*MSl2Diff(2,1)**3

	Sf1441 = MW2*(-Sf1440 + 
     &      MW2*(Sf1437 + MW2*(Sf1433 - MW2*Sf438))) + 
     &   Sf1434*MSl2Diff(2,1)**2

	Sf6336 = 5*MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(3*MW2 - Sf6331) - Sf6333) + 
     &      Sf6335*MSq2Diff(2,2))

	Sf4187 = 5*MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(3*MW2 - Sf4180) - Sf4183) + 
     &      Sf4186*MSq2Diff(2,1))

	Sf4103 = 5*MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(3*MW2 - Sf4096) - Sf4099) + 
     &      Sf4102*MSq2Diff(1,1))

	Sf3986 = MW2**2*Sf3985*MSbot2(2) + 
     &   4*(45*MSbot2(2)**4 + MW2*(25*MW2**3 - 8*MSbot2(2)**3))*
     &    USbot21

	Sf3983 = -(MW2**2*Sf3982*MSbot2(2)) + 
     &   4*(3*MSbot2(2)**4 + MW2*(MW2**3 - 6*MSbot2(2)**3))*
     &    USbot21

	Sf3889 = MW2**2*Sf3888*MSbot2(1) + 
     &   4*(45*MSbot2(1)**4 + MW2*(25*MW2**3 - 8*MSbot2(1)**3))*
     &    USbot11

	Sf3886 = -(MW2**2*Sf3885*MSbot2(1)) + 
     &   4*(3*MSbot2(1)**4 + MW2*(MW2**3 - 6*MSbot2(1)**3))*
     &    USbot11

	Sf8414 = 4*MW2**2*Sf2334**3*Sf4906**3*
     &    (MW2*Sf8409*MSbot2(2) - 
     &      2*CW2**2*MZ2**3*Sf7521*Ss78*SW2*MSbot2(3)**2)*
     &    MSq2Diff(2,1)**2*MSq2Diff(2,2)**2 - 
     &   3*CW2**2*MZ2**3*MSbot2(2)*MSbot2(3)**2*
     &    (Sf8413 - 63*MW2*Sf2334**3*Sf4906**3*MSq2Diff(2,1)**2*
     &       MSq2Diff(2,2)**2)*USbot21**2

	Sf8318 = 4*MW2**2*Sf2320**3*Sf4820**3*
     &    (MW2*Sf8313*MSbot2(1) - 
     &      2*CW2**2*MZ2**3*Sf7501*Ss78*SW2*MSbot2(3)**2)*
     &    MSq2Diff(1,1)**2*MSq2Diff(1,2)**2 - 
     &   3*CW2**2*MZ2**3*MSbot2(1)*MSbot2(3)**2*
     &    (Sf8317 - 63*MW2*Sf2320**3*Sf4820**3*MSq2Diff(1,1)**2*
     &       MSq2Diff(1,2)**2)*USbot22**2

	Sf1807 = -(Sf1806*MStau2(2)**4) + 
     &   MW2*(2*Sf990*MStau2(2)**3 + 
     &      MW2*(MW2*Sf988 - Sf1805*Sf241*MStau2(2)))

	Sf995 = -(Sf994*MStau2(1)**4) + 
     &   MW2*(2*Sf990*MStau2(1)**3 + 
     &      MW2*(MW2*Sf988 - Sf241*Sf993*MStau2(1)))

	Sf1804 = Sf1803*MStau2(2)**4 + 
     &   MW2*(-2*Sf980*MStau2(2)**3 + 
     &      MW2*(3*MW2*Sf978 + 3*Sf1802*MStau2(2)))

	Sf986 = Sf985*MStau2(1)**4 + 
     &   MW2*(-2*Sf980*MStau2(1)**3 + 
     &      MW2*(3*MW2*Sf978 + 3*Sf984*MStau2(1)))

	Sf8596 = 3*Sf8582 + 
     &   MW2*(MW2*(-4*MZ2*Sf8585 + 4*MW2*Sf8587) + Sf8595)

	Sf8575 = 12*MW2*MZ2*
     &    (7*MW2**2 - 12*MW2*MZ2 + 5*MZ2**2 - 
     &      7*MZ2*SW2**2*(MZ2 - 6*MStau2(1))) - 
     &   168*MW2**3*(MZ2 - 3*MStau2(2))*UStau11*UStau12*
     &    UStau21*UStau22 + 
     &   6*MZ2**2*(-MW2**2 + MW2*MZ2 + 42*MW2*MSneu2(1) - 
     &      21*MZ2*MSneu2(1) + 28*MW2*SW2**2*(MZ2 - 6*MStau2(1)) - 
     &      21*MW2*MStau2(2))*UStau22**2 + 
     &   21*MW2*(-MZ2**3 - (-2*MW2 + MZ2)**2*(MZ2 - 6*MStau2(1)) - 
     &      4*MZ2**2*SW2**2*(MZ2 - 6*MStau2(1)) + 
     &      6*MZ2**2*MStau2(2))*UStau22**4 + 
     &   UStau11**2*(-18*
     &       (-7*MZ2**3*MSneu2(1) + 56*MW2**3*MStau2(1) + 
     &         MW2**2*MZ2*(9*MZ2 - 70*MStau2(1) - 42*MStau2(2)) + 
     &         MW2*MZ2**2*(-2*MZ2 + 
     &            7*(2*MSneu2(1) + 3*MStau2(1) + 2*MStau2(2)))) - 
     &      84*MW2*MZ2*(-2*MW2 + MZ2)*SW2*(MZ2 - 6*MStau2(1))*
     &       UStau12**2 - 
     &      42*MW2*(-2*MW2 + MZ2)**2*(MZ2 - 3*MStau2(4))*
     &       UStau21**2 + 
     &      84*MW2*MZ2*(MW2*
     &          (3*MZ2 - 6*SW2*MStau2(1) - 9*MStau2(2)) + 
     &         MZ2*(-MZ2 + 3*SW2*MStau2(1) + 3*MStau2(2)))*
     &       UStau22**2) + 
     &   UStau12**2*(-504*MW2*MZ2**2*SW2**2*MStau2(1) + 
     &      84*MW2*MZ2*(-2*MW2 + MZ2*(1 - 2*SW2))*SW2*
     &       (MZ2 - 3*MStau2(4))*UStau22**2)

	Sf7494 = 27*MZ2**3*Sf7464 + 
     &   MW2*(MW2*(-8*MW2*Sf7473 + 8*MZ2*Sf7481) - MZ2**2*Sf7493)

	Sf7463 = 3*MZ2**3*Sf7436 + 
     &   MW2*(MW2*(-8*MW2*Sf7444 + 8*MZ2*Sf7451) - MZ2**2*Sf7462)

	Sf5072 = Sf5071*MSq2Diff(2,2)**2 + 
     &   MHiggs2(4)*(-Sf5070 + 
     &      MHiggs2(4)*(Sf5065*MHiggs2(4) + Sf5064*MSq2Diff(2,2)))

	Sf4946 = Sf4945*MSq2Diff(1,2)**2 + 
     &   MHiggs2(4)*(-Sf4944 + 
     &      MHiggs2(4)*(Sf4939*MHiggs2(4) + Sf4938*MSq2Diff(1,2)))

	Sf3955 = 20*Sf3954*MSbot2(2)**3 + 
     &   MW2*(48*Sf3952*MSbot2(2)**2 + 
     &      MW2*(10*MW2*Sf3624 + Sf3953*MSbot2(2)))

	Sf6013 = Sf3622 - 
     &   MStop2(2)*(Sf3632 + 
     &      MStop2(2)*(-Sf3615 + 4*MStop2(2)*USbot21))

	Sf3657 = Sf3622 - 
     &   MStop2(1)*(Sf3632 + 
     &      MStop2(1)*(-Sf3615 + 4*MStop2(1)*USbot21))

	Sf3842 = 20*Sf3841*MSbot2(1)**3 + 
     &   MW2*(48*Sf3837*MSbot2(1)**2 + 
     &      MW2*(10*MW2*Sf3307 + Sf3840*MSbot2(1)))

	Sf5593 = Sf3302 - 
     &   MStop2(2)*(Sf3323 + 
     &      MStop2(2)*(-Sf3290 + 4*MStop2(2)*USbot11))

	Sf3377 = Sf3302 - 
     &   MStop2(1)*(Sf3323 + 
     &      MStop2(1)*(-Sf3290 + 4*MStop2(1)*USbot11))

	Sf306 = -(Sf287*MSl2Diff(1,1)**3) + 
     &   MHiggs2(4)*(MHiggs2(4)*(-Sf298 + Sf291*MHiggs2(4)) + 
     &      Sf305*MSl2Diff(1,1))

	Sf2623 = Sf2622*MSq2Diff(2,1)**2 + 
     &   MHiggs2(4)*(-Sf2621 + 
     &      MHiggs2(4)*(Sf2614*MHiggs2(4) + Sf2613*MSq2Diff(2,1)))

	Sf2401 = Sf2400*MSq2Diff(1,1)**2 + 
     &   MHiggs2(4)*(-Sf2399 + 
     &      MHiggs2(4)*(Sf2386*MHiggs2(4) + Sf2384*MSq2Diff(1,1)))

	Sf194 = 2*Sf189*MStau2(1)**3 + 
     &   MW2*(Sf188*MStau2(1)**2 + 
     &      MW2*(MW2*Sf193 + 5*MStau2(1)**4))

	Sf1266 = -(Sf1255*MSl2Diff(2,1)**3) + 
     &   MHiggs2(4)*(MHiggs2(4)*(-Sf1261 + Sf1258*MHiggs2(4)) + 
     &      Sf1265*MSl2Diff(2,1))

	Sf1202 = 2*Sf189*MStau2(2)**3 + 
     &   MW2*(Sf188*MStau2(2)**2 + 
     &      MW2*(MW2*Sf1201 + 5*MStau2(2)**4))

	Sf6794 = MW2*(32*Sf6203 + MW2*(16*MW2 - 64*Sf99)) - 
     &   Sf6793*MHiggs2(4)

	Sf6539 = MW2*(32*Sf5825 + MW2*(16*MW2 - 64*Sf91)) - 
     &   Sf6538*MHiggs2(4)

	Sf4441 = MW2*(32*Sf4137 + MW2*(16*MW2 - 64*Sf77)) - 
     &   Sf4440*MHiggs2(4)

	Sf4314 = MW2*(32*Sf4018 + MW2*(16*MW2 - 64*Sf69)) - 
     &   Sf4313*MHiggs2(4)

	Sf6206 = MW2*(16*Sf6203 + MW2*(8*MW2 - 32*Sf99)) + 
     &   Sf6205*MHiggs2(4)

	Sf5828 = MW2*(16*Sf5825 + MW2*(8*MW2 - 32*Sf91)) + 
     &   Sf5827*MHiggs2(4)

	Sf4142 = MW2*(16*Sf4137 + MW2*(8*MW2 - 32*Sf77)) + 
     &   Sf4141*MHiggs2(4)

	Sf4025 = MW2*(16*Sf4018 + MW2*(8*MW2 - 32*Sf69)) + 
     &   Sf4024*MHiggs2(4)

	Sf8350 = Sf166**3*Sf4906**3*Sf8345*USbot22**2*
     &    UStop21**2 + Sf4820*
     &    (-(Sf4820**2*Sf8349) + 
     &      Sf166**3*Sf4906**3*Sf8344*USbot22**2*UStop21**2)

	Sf8351 = Sf166**3*Sf4820**3*Sf4906**3*Sf8343*
     &    USbot21**2*UStop22**2 + 
     &   Sf2334*(Sf2334**2*Sf8350 + 
     &      Sf166**3*Sf4820**3*Sf4906**3*Sf8342*USbot21**2*
     &       UStop22**2)

	Sf8352 = Sf166**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf8341*
     &    USbot22**2*UStop22**2 + 
     &   Sf2320*(Sf2320**2*Sf8351 + 
     &      Sf166**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf8340*
     &       USbot22**2*UStop22**2)

	Sf1361 = (Sf440 - MHiggs2(4)*MStau2(2))*UStau21 - 
     &   ML*Sf215*MHiggs2(4)*UStau22

	Sf1994 = ML2*Sf1987*UStau21 - Sf1993*UStau21 - 
     &   ML*Sf1987*Sf215*UStau22

	Sf1843 = ML2*Sf1835*UStau21 + Sf1842*UStau21 - 
     &   ML*Sf1835*Sf215*UStau22

	Sf441 = (Sf440 - MHiggs2(4)*MStau2(1))*UStau11 - 
     &   ML*Sf215*MHiggs2(4)*UStau12

	Sf1936 = ML2*Sf1928*UStau11 - Sf1935*UStau11 - 
     &   ML*Sf1928*Sf215*UStau12

	Sf1064 = ML2*Sf1046*UStau11 + Sf1063*UStau11 - 
     &   ML*Sf1046*Sf215*UStau12

	Sf8713 = Sf4906**3*Sf8709*MSq2Diff(2,2)**2*
     &    USbot22**2 + Sf4820*
     &    (Sf4820**2*Sf8712*MSq2Diff(1,2)**2*USbot21**2 + 
     &      Sf4906**3*Sf8707*MSq2Diff(2,2)**2*USbot22**2)

	Sf8676 = Sf2334**3*Sf8672*MSq2Diff(2,1)**2*
     &    USbot22**2 + Sf2320*
     &    (Sf2320**2*Sf8675*MSq2Diff(1,1)**2*USbot21**2 + 
     &      Sf2334**3*Sf8670*MSq2Diff(2,1)**2*USbot22**2)

	Sf7755 = 4*MW2**3*Sf1195**3*Sf7753*MSl2Diff(2,1)**2 - 
     &   CW2**2*MZ2**3*(Sf7754 - 
     &      21*MW2*Sf1195**3*MSl2Diff(2,1)**2)*MStau2(2)*
     &    MStau2(3)**2*UStau21**2

	Sf7733 = 4*MW2**3*Sf166**3*Sf7731*MSl2Diff(1,1)**2 - 
     &   CW2**2*MZ2**3*(Sf7732 - 21*MW2*Sf166**3*MSl2Diff(1,1)**2)*
     &    MStau2(1)*MStau2(3)**2*UStau22**2

	Sf5146 = MW2*(3*Sf5145 + 12*MW2*Sf4864*Sf5140*TB) - 
     &   2*Sf5142*MHiggs2(4)

	Sf5027 = MW2*(3*Sf5026 + 12*MW2*Sf4828*Sf5021*TB) - 
     &   2*Sf5023*MHiggs2(4)

	Sf2738 = MW2*(3*Sf2737 + 12*MW2*Sf2263*Sf2732*TB) - 
     &   2*Sf2734*MHiggs2(4)

	Sf2560 = MW2*(3*Sf2559 + 12*MW2*Sf2198*Sf2554*TB) - 
     &   2*Sf2556*MHiggs2(4)

	Sf2011 = -(Sf1936*Ss25**3) + 
     &   Sf166*Ss241*Ss242*
     &    (Sf1956*Ss241*Ss242 - 2*Sf166*Sf1923*UStau11)

	Sf898 = MW2*Sf897 + 
     &   MHiggs2(4)*(-(MW2*Sf882) + 2*Sf867*MHiggs2(4))

	Sf8482 = -(Sf4906**3*Sf8477*USbot22**2*
     &      UStop21**2) + 
     &   Sf4820*(Sf4820**2*Sf8481 - 
     &      Sf4906**3*Sf8476*USbot22**2*UStop21**2)

	Sf7275 = 2*Sf7272 + MZ2*(3*MZ2**2*Sf7273 + Sf7274)

	Sf205 = -6*Sf202*Ss241*Ss242*MStau2(1)**2 + 
     &   MW2*(MW2*Sf204 + Sf201*MStau2(1))

	Sf1697 = MW2*Sf1696 + 
     &   MHiggs2(4)*(-(MW2*Sf1688) + 2*Sf1680*MHiggs2(4))

	Sf1169 = MW2*(8*MW2 - MZ2*Sf1165) - MZ2**2*Sf1168

	Sf9373 = Sf8923*Sf9372 + 
     &   MHiggs2(2)**2*(-2*Sf9315 + 3*MHiggs2(2))

	Sf8996 = Sf8923*Sf8995 + 
     &   MHiggs2(2)**2*(-2*Sf8902 + 3*MHiggs2(2))

	Sf10496 = Sf10495*Sf8923 + 
     &   MHiggs2(2)**2*(-2*Sf10265 + 3*MHiggs2(2))

	Sf10120 = Sf10119*Sf9125 + 
     &   MHiggs2(1)**2*(-2*Sf9315 + 3*MHiggs2(1))

	Sf9181 = Sf9125*Sf9180 + 
     &   MHiggs2(1)**2*(-2*Sf8902 + 3*MHiggs2(1))

	Sf10387 = Sf10386*Sf9125 + 
     &   MHiggs2(1)**2*(-2*Sf10265 + 3*MHiggs2(1))

	Sf5989 = Sf3578*TB2*UStop21**2 - 
     &   Sf2245*Sf5442*USbot21*UStop22

	Sf5563 = Sf3226*TB2*UStop21**2 - 
     &   Sf2245*Sf5336*USbot11*UStop22

	Sf5443 = Sf2961*TB2*UStop21**2 - 
     &   Sf2245*Sf5442*USbot21*UStop22

	Sf5372 = Sf2961*TB2*UStop21**2 + 
     &   Sf2245*Sf5371*USbot21*UStop22

	Sf5337 = Sf2809*TB2*UStop21**2 - 
     &   Sf2245*Sf5336*USbot11*UStop22

	Sf5264 = Sf2809*TB2*UStop21**2 + 
     &   Sf2245*Sf5263*USbot11*UStop22

	Sf4898 = -2*Sf2314*UStop22**2 + Sf2313*(1 - UStop22**2)

	Sf8734 = -(MW2**2*Sf7219*Sf8733*UStop11) + 
     &   MZ2*Sf8732*UStop21

	Sf6982 = Sf6981*UStop21 + Sf6977*UStop22

	Sf3579 = -(Sf2245*Sf3052*USbot21*UStop12) + 
     &   Sf3578*TB2*UStop22**2

	Sf3227 = -(Sf2245*Sf2925*USbot11*UStop12) + 
     &   Sf3226*TB2*UStop22**2

	Sf3053 = -(Sf2245*Sf3052*USbot21*UStop12) + 
     &   Sf2961*TB2*UStop22**2

	Sf2964 = Sf2245*Sf2963*USbot21*UStop12 + 
     &   Sf2961*TB2*UStop22**2

	Sf2926 = -(Sf2245*Sf2925*USbot11*UStop12) + 
     &   Sf2809*TB2*UStop22**2

	Sf2812 = Sf2245*Sf2811*USbot11*UStop12 + 
     &   Sf2809*TB2*UStop22**2

	Sf2315 = Sf2313*UStop22**2 + Sf2314*(-2 + 2*UStop22**2)

	Sf7180 = Sf2263*MSq2Diff(2,1)*USbot21*UStop11 + 
     &   Sf7179*UStop21

	Sf6384 = Sf2263*MSq2Diff(2,1)*USbot21*UStop11 + 
     &   Sf6383*UStop21

	Sf6400 = Sf2263*Sf4867*Sf6259*UStop11 + 
     &   Sf2270*Sf4864*Sf6262*UStop21

	Sf6398 = Sf2198*Sf4831*Sf5895*UStop11 + 
     &   Sf2206*Sf4828*Sf5898*UStop21

	Sf4265 = Sf4264*UStop11 + Sf4260*UStop12

	Sf6251 = Sf2263*UStop11 + Sf4864*UStop21

	Sf5887 = Sf2198*UStop11 + Sf4828*UStop21

	Sf1981 = ML*Sf1980*UStau21 - 
     &   ML2*MW2**2*Sf1459*Sf215*Sf44*MStau2(2)*UStau22

	Sf1727 = Sf1697*UStau21 - ML*MW2*Sf1673*Sf215*UStau22

	Sf1516 = Sf1479*UStau21 - 3*ML*Sf1484*Sf215*UStau22

	Sf1863 = Sf1093*UStau21 - 3*ML*Sf215*MHiggs2(4)*UStau22

	Sf1826 = Sf1033*UStau21 - ML*Sf215*Ss115*Ss25*UStau22

	Sf938 = Sf898*UStau11 - ML*MW2*Sf215*Sf857*UStau12

	Sf8612 = -15*MZ2**5*Sf8551 + Sf8610*UStau11

	Sf621 = Sf578*UStau11 - 3*ML*Sf215*Sf583*UStau12

	Sf1921 = ML*Sf1920*UStau11 - 
     &   ML2*MW2**2*Sf215*Sf31*Sf556*MStau2(1)*UStau12

	Sf1094 = Sf1093*UStau11 - 3*ML*Sf215*MHiggs2(4)*UStau12

	Sf1034 = Sf1033*UStau11 - ML*Sf215*Ss115*Ss25*UStau12

	Sf6087 = Sf5692*USbot21**2 - Sf2225*Sf6086*USbot22

	Sf6082 = Sf5682*USbot21**2 - Sf2225*Sf6081*USbot22

	Sf6075 = Sf5667*USbot21**2 - 3*Sf2225*Sf6074*USbot22

	Sf3966 = Sf3866*USbot21**2 - Sf2225*Sf3965*USbot22

	Sf3960 = Sf3853*USbot21**2 - Sf2225*Sf3959*USbot22

	Sf3653 = Sf3370*USbot21**2 - Sf2225*Sf3652*USbot22

	Sf3644 = Sf3350*USbot21**2 - 2*Sf2225*Sf3643*USbot22

	Sf138 = -(Sf121*USbot22**2) + Sf120*(1 - USbot22**2)

	Sf6914 = Sf6905*USbot21 + MBot*Sf2225*Sf6913*USbot22

	Sf5431 = Sf5423*USbot21 + 9*MBot*Sf2225*Sf5430*USbot22

	Sf4721 = Sf4709*USbot21 + MBot*Sf2225*Sf4720*USbot22

	Sf3035 = Sf3023*USbot21 + 9*MBot*Sf2225*Sf3034*USbot22

	Sf5695 = -(Sf2225*Sf5694*USbot12) + Sf5692*USbot22**2

	Sf5685 = -(Sf2225*Sf5684*USbot12) + Sf5682*USbot22**2

	Sf5670 = -3*Sf2225*Sf5669*USbot12 + Sf5667*USbot22**2

	Sf3869 = -(Sf2225*Sf3868*USbot12) + Sf3866*USbot22**2

	Sf3856 = -(Sf2225*Sf3855*USbot12) + Sf3853*USbot22**2

	Sf3373 = -(Sf2225*Sf3372*USbot12) + Sf3370*USbot22**2

	Sf3353 = -2*Sf2225*Sf3352*USbot12 + Sf3350*USbot22**2

	Sf122 = Sf120*USbot22**2 + Sf121*(-1 + USbot22**2)

	Sf6663 = Sf6654*USbot11 + MBot*Sf2225*Sf6662*USbot12

	Sf6414 = Sf4828*Sf4867*Sf6411*USbot11 + 
     &   Sf4831*Sf4864*Sf6413*USbot21

	Sf5325 = Sf5317*USbot11 + 9*MBot*Sf2225*Sf5324*USbot12

	Sf5218 = Sf4828*Sf4867*Sf5215*USbot11 + 
     &   Sf4831*Sf4864*Sf5217*USbot21

	Sf6402 = Sf4828*Sf4867*Sf5214*USbot11 + 
     &   Sf4831*Sf4864*Sf5216*USbot21

	Sf4626 = Sf4612*USbot11 + MBot*Sf2225*Sf4625*USbot12

	Sf6396 = Sf2198*Sf2270*Sf4209*USbot11 + 
     &   Sf2206*Sf2263*Sf4212*USbot21

	Sf2908 = Sf2896*USbot11 + 9*MBot*Sf2225*Sf2907*USbot12

	Sf6289 = Sf3986*MStop2(2) - Sf6288*USbot21

	Sf6275 = Sf3889*MStop2(2) - Sf6274*USbot11

	Sf3987 = Sf3986*MStop2(1) - Sf3980*USbot21

	Sf3890 = Sf3889*MStop2(1) - Sf3883*USbot11

	Sf10590 = Sf10225*Sf10226*Sf10437*MStop2(1) - Sf10589*MStop2(2)

	Sf10326 = Sf10225*Sf10226*Sf10321*MStop2(1) - Sf10325*MStop2(2)

	Sf9589 = 2*MZ2*Sf9290*Sf9291*Sf9332*Sf9333*Sf9547*
     &    MHiggs2(2)*MStau2(1)**2 - Sf9588*MStau2(2)

	Sf9530 = Sf9499*MStau2(1)**2 - Sf9529*MStau2(2)

	Sf9531 = MZ2*Sf9290*Sf9291*Sf9332*Sf9333*Sf9492*
     &    MHiggs2(2)*MStau2(1)**2 - Sf9530*MStau2(2)

	Sf9668 = 16*Sf9329*Sf9332**3*Sf9333**3*MHiggs2(2)*
     &    MStau2(1)**2 + Sf9667*MStau2(2)

	Sf9638 = 4*Sf9329*Sf9332**3*Sf9333**3*MHiggs2(2)*
     &    MStau2(1)**2 + Sf9637*MStau2(2)

	Sf9953 = Sf9935*MStau2(1) + Sf9952*MStau2(2)

	Sf9475 = Sf9440*MStau2(1) + Sf9474*MStau2(2)

	Sf9864 = MZ2*MHiggs2(2)*(MStau2(1) + 3*MStau2(2))*
     &    MStau2(3) - MHiggs2(2)*MStau2(3)**3 - 
     &   MZ2**2*(2*MHiggs2(2)*MStau2(2) + MStau2(3)**2) + 
     &   MZ2**3*MStau2(4)

	Sf10215 = Sf10056*Sf9290*Sf9291*MStau2(1) - Sf10214*MStau2(2)

	Sf1146 = (MW2*Sf1145*Sf162)/(Sf154*MStau2(1)) + 
     &   8*Sf162*Ss72*UStau11

	Sf1715 = Sf919*MSneu2(1)**3 + Sf1714*MStau2(2)

	Sf925 = Sf919*MSneu2(1)**3 + Sf924*MStau2(1)

	Sf1619 = Sf783*MSneu2(1)**3 + Sf1618*MStau2(2)

	Sf789 = Sf783*MSneu2(1)**3 + Sf788*MStau2(1)

	Sf1323 = 8*Sf381*MSneu2(1)**3 + Sf1322*MStau2(2)

	Sf386 = 8*Sf381*MSneu2(1)**3 + Sf385*MStau2(1)

	Sf1291 = Sf333*MSneu2(1)**3 + Sf1290*MStau2(2)

	Sf339 = Sf333*MSneu2(1)**3 + Sf338*MStau2(1)

	Sf2159 = 4*Sf2084*MSneu2(1)**3 + Sf2158*MStau2(2)

	Sf2090 = 4*Sf2084*MSneu2(1)**3 + Sf2089*MStau2(1)

	Sf1858 = Sf1079*MSneu2(1)**3 - Sf1857*MStau2(2)

	Sf1085 = Sf1079*MSneu2(1)**3 - Sf1084*MStau2(1)

	Sf1998 = Sf1060*MSneu2(1)**3 + Sf1997*MStau2(2)

	Sf1940 = Sf1060*MSneu2(1)**3 + Sf1939*MStau2(1)

	Sf1150 = Sf1146*Sf154*Ss25*MSneu2(1) + 
     &   6*MW2*Sf1042*Sf162*MSl2Diff(1,1)*UStau11

	Sf6407 = Sf165*MSneu2(1) - Sf159*MStau2(1)

	Sf6406 = Sf1194*MSneu2(1) - Sf1189*MStau2(2)

	Sf6791 = Sf6788*MSbot2(2)**3 - Sf6790*MStop2(2)**2

	Sf6868 = 2*Sf6494*MSbot2(2)**3 + Sf6867*MStop2(2)

	Sf6754 = 2*Sf6494*MSbot2(2)**3 + Sf6753*MStop2(1)

	Sf7095 = Sf4554*MSbot2(2)**3 + Sf7094*MStop2(2)

	Sf4664 = Sf4554*MSbot2(2)**3 + Sf4663*MStop2(1)

	Sf4436 = Sf4433*MSbot2(2)**3 - Sf4435*MStop2(1)**2

	Sf6159 = Sf3777*MSbot2(2)**3 - Sf6158*MStop2(2)

	Sf3911 = Sf3777*MSbot2(2)**3 - Sf3910*MStop2(1)

	Sf5976 = Sf3190*MSbot2(2)**3 + Sf5975*MStop2(2)

	Sf3551 = Sf3190*MSbot2(2)**3 + Sf3550*MStop2(1)

	Sf5920 = Sf3089*MSbot2(2)**3 + Sf5919*MStop2(2)

	Sf3476 = Sf3089*MSbot2(2)**3 + Sf3475*MStop2(1)

	Sf5360 = 14*Sf2791*MSbot2(2)**3 + Sf5359*MStop2(2)

	Sf2949 = 14*Sf2791*MSbot2(2)**3 + Sf2948*MStop2(1)

	Sf5079 = Sf2632*MSbot2(2)**3 + Sf5078*MStop2(2)

	Sf2638 = Sf2632*MSbot2(2)**3 + Sf2637*MStop2(1)

	Sf5438 = Sf5436*MSbot2(2)**2 - Sf5437*MStop2(2)**2

	Sf6919 = Sf4596*MSbot2(2)**2 - Sf6918*MStop2(2)

	Sf4731 = Sf4596*MSbot2(2)**2 - Sf4730*MStop2(1)

	Sf6808 = Sf4463*MSbot2(2)**2 - Sf6807*MStop2(2)

	Sf4468 = Sf4463*MSbot2(2)**2 - Sf4467*MStop2(1)

	Sf5995 = 3*Sf3585*MSbot2(2)**2 + Sf5994*MStop2(2)

	Sf3590 = 3*Sf3585*MSbot2(2)**2 + Sf3589*MStop2(1)

	Sf5939 = 3*Sf3498*MSbot2(2)**2 + Sf5938*MStop2(2)

	Sf3503 = 3*Sf3498*MSbot2(2)**2 + Sf3502*MStop2(1)

	Sf3046 = Sf3044*MSbot2(2)**2 - Sf3045*MStop2(1)**2

	Sf6358 = 2*Sf6144*Ss241*Ss242*Ss25*MSbot2(2) + 
     &   MW2*Sf5216*MSq2Diff(2,2)

	Sf4213 = 2*Sf3742*Ss241*Ss242*Ss25*MSbot2(2) + 
     &   MW2*Sf4212*MSq2Diff(2,1)

	Sf6346 = Sf3560*Sf3561*Ss25*MSbot2(2) + Sf6345*MStop2(2)

	Sf4205 = Sf3560*Sf3561*Ss25*MSbot2(2) + Sf4204*MStop2(1)

	Sf6217 = Sf4155*MSbot2(2) + Sf6216*MStop2(2)

	Sf4158 = Sf4155*MSbot2(2) + Sf4157*MStop2(1)

	Sf6616 = Sf6613*MSbot2(1)**3 + Sf6615*MStop2(2)

	Sf6536 = Sf6533*MSbot2(1)**3 - Sf6535*MStop2(2)**2

	Sf6499 = 2*Sf6494*MSbot2(1)**3 + Sf6498*MStop2(1)

	Sf7055 = Sf4554*MSbot2(1)**3 + Sf7054*MStop2(2)

	Sf4559 = Sf4554*MSbot2(1)**3 + Sf4558*MStop2(1)

	Sf4308 = Sf4305*MSbot2(1)**3 - Sf4307*MStop2(1)**2

	Sf5780 = Sf3777*MSbot2(1)**3 - Sf5779*MStop2(2)

	Sf3783 = Sf3777*MSbot2(1)**3 - Sf3782*MStop2(1)

	Sf5773 = Sf3765*MSbot2(1)**3 + Sf5772*MStop2(2)

	Sf5550 = Sf3190*MSbot2(1)**3 + Sf5549*MStop2(2)

	Sf3195 = Sf3190*MSbot2(1)**3 + Sf3194*MStop2(1)

	Sf5487 = Sf3089*MSbot2(1)**3 + Sf5486*MStop2(2)

	Sf3093 = Sf3089*MSbot2(1)**3 + Sf3092*MStop2(1)

	Sf5252 = 14*Sf2791*MSbot2(1)**3 + Sf5251*MStop2(2)

	Sf2796 = 14*Sf2791*MSbot2(1)**3 + Sf2795*MStop2(1)

	Sf4953 = Sf2413*MSbot2(1)**3 + Sf4952*MStop2(2)

	Sf2421 = Sf2413*MSbot2(1)**3 + Sf2420*MStop2(1)

	Sf5332 = Sf5330*MSbot2(1)**2 - Sf5331*MStop2(2)**2

	Sf6668 = Sf4596*MSbot2(1)**2 - Sf6667*MStop2(2)

	Sf4636 = Sf4596*MSbot2(1)**2 - Sf4635*MStop2(1)

	Sf6553 = Sf4343*MSbot2(1)**2 - Sf6552*MStop2(2)

	Sf4351 = Sf4343*MSbot2(1)**2 - Sf4350*MStop2(1)

	Sf5569 = 3*Sf3239*MSbot2(1)**2 + Sf5568*MStop2(2)

	Sf3246 = 3*Sf3239*MSbot2(1)**2 + Sf3245*MStop2(1)

	Sf5506 = 3*Sf3121*MSbot2(1)**2 + Sf5505*MStop2(2)

	Sf3127 = 3*Sf3121*MSbot2(1)**2 + Sf3126*MStop2(1)

	Sf2919 = Sf2917*MSbot2(1)**2 - Sf2918*MStop2(1)**2

	Sf7125 = 2*Sf6519*Ss25**3*MSbot2(1) - 
     &   MW2*Sf6410*Ss241*Ss242*MHiggs2(4)*MSq2Diff(1,2)

	Sf4747 = 2*Sf4267*Ss25**3*MSbot2(1) - 
     &   MW2*Sf4746*Ss241*Ss242*MHiggs2(4)*MSq2Diff(1,1)

	Sf6357 = 2*Sf5761*Ss241*Ss242*Ss25*MSbot2(1) + 
     &   MW2*Sf5214*MSq2Diff(1,2)

	Sf4210 = 2*Sf3733*Ss241*Ss242*Ss25*MSbot2(1) + 
     &   MW2*Sf4209*MSq2Diff(1,1)

	Sf6323 = Sf3204*Sf3205*Ss25*MSbot2(1) + Sf6322*MStop2(2)

	Sf4125 = Sf3204*Sf3205*Ss25*MSbot2(1) + Sf4124*MStop2(1)

	Sf9283 = Sf8873*Sf8874*Sf9109*MSbot2(1) - Sf9282*MSbot2(2)

	Sf9102 = Sf8873*Sf8874*Sf8905*MSbot2(1) - Sf9101*MSbot2(2)

	Sf5840 = Sf4047*MSbot2(1) + Sf5839*MStop2(2)

	Sf4053 = Sf4047*MSbot2(1) + Sf4052*MStop2(1)

	Sf6965 = 9*Sf6964*TB - Sf6962*USbot11

	Sf7226 = 8*MZ2**3*SW2**2 - Sf7225*MStop2(2)

	Sf7207 = 8*MZ2**3*SW2**2 - Sf7206*MStop2(1)

	Sf1784 = 6*MZ2**3*SW2**2 - 3*Sf1782*MStau2(2)

	Sf1764 = 6*MZ2**3*SW2**2 - 3*Sf1762*MStau2(1)

	Sf4806 = 2*MZ2**3*SW2**2 + Sf4801*MSbot2(2)

	Sf4785 = 2*MZ2**3*SW2**2 + Sf4780*MSbot2(1)

	Sf8761 = 4*Sf8760*SW2 + Sf7756*Sf86**2*MStop2(1)*UStop22**4

	Sf8739 = 4*Sf8738*SW2 - Sf7756*Sf86**2*MStop2(2)*UStop21**4

	Sf8518 = 2*Sf8506*SW2 - Sf15**2*Sf7756*MSbot2(2)*USbot21**4

	Sf8502 = 2*Sf8488*SW2 + Sf15**2*Sf7756*MSbot2(1)*USbot22**4

	Sf7787 = 2*Sf7786*SW2 + Sf7756*Ss443**2*MStau2(1)*UStau22**4

	Sf7768 = 2*Sf7767*SW2 - Sf7756*Ss443**2*MStau2(2)*UStau21**4

	Sf10053 = (MZ2*Sf10052*Ss49 + 
     &     Sf10035/(Sf9344*MStau2(1)*MStau2(2)))/MZ2

	Sf2181 = Sf4*Ss25**3 + 3*MW2*Sf7*Ss115*Ss25*TB

	Sf2058 = Sf2053*Sf976*Sf977*Ss25 + Sf2057*Ss241*Ss242*MSneu2(1)

	Sf738 = -((-MW2 + MHiggs2(4))*(MW2 - MStau2(1))**2*
     &      (MW2*(-2*MW2 + 3*MHiggs2(4)) - MHiggs2(4)*MStau2(1)))-
     &     (-MW + MHiggs(4))*(MW + MHiggs(4))*MSneu2(1)*
     &    (MW2**2*(MW2 + MSneu2(1) + MStau2(1)) + 
     &      MHiggs2(4)*(MW2*(-2*MW2 - 2*MSneu2(1) - 3*MStau2(1)) + 
     &         MSneu2(1)*(MSneu2(1) - 3*MStau2(1)) + 3*MStau2(1)**2
     &         ))

	Sf6846 = (MW2 - MHiggs2(4))*
     &   (-((MW2 - MSq2Diff(2,2))**2*
     &        (MW2*(2*MW2 - 3*MHiggs2(4)) + 
     &          MHiggs2(4)*MSq2Diff(2,2))) + 
     &     MW2*(5*MW2**2 - 9*MW2*MHiggs2(4) + 
     &        (-3*MW2 + 7*MHiggs2(4))*MSq2Diff(2,2))*MStop2(2))

	Sf6732 = -((-MW2 + MHiggs2(4))*
     &     (MW2**2*(-2*(MW2 - MSbot2(2))**2 + 
     &          (MW2 + MSbot2(2))*MStop2(1) + MStop2(1)**2) + 
     &       MHiggs2(4)*(3*MW2**3 - MSq2Diff(2,1)**3 + 
     &          MW2*MSq2Diff(2,1)*(5*MSbot2(2) + 2*MStop2(1)) - 
     &          MW2**2*(7*MSbot2(2) + 2*MStop2(1)))))

	Sf6591 = Sf2543**2*Sf6460*Ss25 + Sf6590*Ss241*Ss242*MStop2(2)

	Sf6465 = -((-MW2 + MHiggs2(4))*
     &     (MW2**2*(-2*(MW2 - MSbot2(1))**2 + 
     &          (MW2 + MSbot2(1))*MStop2(1) + MStop2(1)**2) + 
     &       MHiggs2(4)*(3*MW2**3 - MSq2Diff(1,1)**3 + 
     &          MW2*MSq2Diff(1,1)*(5*MSbot2(1) + 2*MStop2(1)) - 
     &          MW2**2*(7*MSbot2(1) + 2*MStop2(1)))))

	Sf7033 = -((-MW2 + MHiggs2(4))*
     &     (-2*MW2**4 + 3*MW2**3*MHiggs2(4) + 
     &       MW2**2*(MW2 - 2*MHiggs2(4))*MSq2Diff(2,2) + 
     &       MW2*(MW2 - 2*MHiggs2(4))*MSq2Diff(2,2)**2 + 
     &       MHiggs2(4)*MSq2Diff(2,2)**3 + 
     &       MW2**2*(5*MW2 - 9*MHiggs2(4))*MStop2(2) + 
     &       MW2*(3*MW2 - 7*MHiggs2(4))*MSq2Diff(2,2)*MStop2(2)))

	Sf4526 = Sf3560*Sf3561*Sf4519*Ss25 - 
     &   Sf4525*Ss241*Ss242*MStop2(1)

	Sf7011 = Sf3204*Sf3205*Sf4415*Ss25 - 
     &   Sf7010*Ss241*Ss242*MStop2(2)

	Sf4422 = -((-MW2 + MHiggs2(4))*
     &     (MW2**2*(-2*MStop2(1)**2 + 
     &          MSbot2(1)*(MSbot2(1) + MStop2(1)) + 
     &          MW2*(-2*MW2 + MSbot2(1) + 4*MStop2(1))) + 
     &       MHiggs2(4)*(3*MW2**3 + MSq2Diff(1,1)**3 - 
     &          MW2*MSq2Diff(1,1)*(2*MSbot2(1) + 5*MStop2(1)) - 
     &          MW2**2*(2*MSbot2(1) + 7*MStop2(1)))))

	Sf2138 = Sf1800*Sf1801*Sf2133*Ss25 + 
     &   Sf2137*Ss241*Ss242*MSneu2(1)

	Sf1598 = Sf1196**2*Sf1593*Ss25 + Sf1597*Ss241*Ss242*MSneu2(1)

	Sf6347 = (-MW2 + MHiggs2(4))*
     &   (-MW2**3 + MW2*MSq2Diff(2,2)*
     &      (MHiggs2(4) + 2*MSq2Diff(2,2)) - 
     &     MSq2Diff(2,2)*(MHiggs2(4)*MSbot2(2) + 
     &        MSq2Diff(2,2)**2) + 3*MW2**2*MStop2(2) + 
     &     MW2*(MHiggs2(4) + 5*MSq2Diff(2,2))*MStop2(2))

	Sf6324 = Sf6323*Ss241*Ss242 - Sf2543*Sf3204*Sf3205*Sf4120*Ss25

	Sf4206 = Sf4205*Ss241*Ss242 - Sf2722*Sf3560*Sf3561*Sf4200*Ss25

	Sf4126 = Sf4125*Ss241*Ss242 - Sf2543*Sf3204*Sf3205*Sf4120*Ss25

	Sf10004 = Sf8929*Sf8930*Sf9979 - 
     &   (Sf10003*Sf55*Sf57)/
     &    (Sf9344**2*Sf9362**2*MStau2(1)*MStau2(2))

	Sf9901 = Sf9317*Sf9881 - 2*Sf9890*MStau2(1)*MStau2(2)

	Sf9837 = Sf55*Sf57*Sf9836 + Sf9831/(Sf9344**2*MStau2(1))

	Sf9816 = Sf55*Sf57*Sf9815 - Sf9810/(Sf9344**2*MStau2(1))

	Sf9794 = Sf9362*Sf9793 - MZ2*Sf9344*Sf9774*Ss49

	Sf1004 = 3*Sf1001*Sf976*Sf977 - Sf1003*MSneu2(1)

	Sf9707 = (Sf9706 - 
     &     (MZ2*Sf9692)/(Sf9362**2*MHiggs2(2)*MStau2(1)))/MZ2

	Sf9684 = Sf55*Sf57*Sf9683 + Sf9677/(Sf9317**2*MStau2(2))

	Sf9654 = Sf55*Sf57*Sf9653 - Sf9646/(Sf9317**2*MStau2(2))

	Sf9590 = 2*Sf9331*Sf9589 - 2*MZ2*Sf9317*Sf9542*Ss49

	Sf9416 = (Sf9415 - 
     &     (MZ2*Sf9397)/(Sf9331**2*MHiggs2(2)*MStau2(2)))/MZ2

	Sf10144 = Sf10141*Sf9302*Sf9303 - Sf10143*MStau2(1)

	Sf9255 = Sf9161*Sf9162*Sf9252 - Sf9254*MSbot2(1)

	Sf9209 = Sf8885*Sf8886*Sf9206 - Sf9208*MSbot2(1)

	Sf5774 = Sf5773*MSq2Diff(1,2) + 12*MBot2*MTop2*Sf91*MSq2Diff(1,2)

	Sf8353 = Sf1195**3*Sf8352 - 
     &   Sf166**3*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf8339*
     &    UStau21**2

	Sf842 = ML*Sf341*Sf834 + Sf841*UStau11

	Sf795 = 3*ML2*Sf789 + Sf794*MSl2Diff(1,1)

	Sf7495 = 3*Sf7463 + Sf7494

	Sf7276 = (6*MZ2**3*Sf7270 - Sf7275)/MZ2**3

	Sf6817 = Sf6814 - Sf6794*MSq2Diff(2,2)**2

	Sf6818 = Sf6799 - 2*Sf6786*MSq2Diff(2,2)**2

	Sf6815 = 2*Sf6791 + MW2*Sf6802

	Sf6562 = Sf6559 - Sf6539*MSq2Diff(1,2)**2

	Sf6563 = Sf6544 - 2*Sf6527*MSq2Diff(1,2)**2

	Sf6560 = 2*Sf6536 + MW2*Sf6547

	Sf6250 = 3*MW2*Sf5056*Sf6249 + Sf6247*TB*USbot21*UStop21

	Sf618 = Sf616 - Sf617*MSneu2(1)

	Sf6142 = 3*MW2*Sf2602*Sf6141 + Sf6136*TB*USbot21*UStop11

	Sf6093 = MW2**3*Sf6087 - Sf6092*USbot21

	Sf6084 = Sf6082 + 33*Sf2731*MSbot2(2)*USbot21

	Sf6077 = Sf6075 - Sf6072*MSbot2(2)*USbot21

	Sf6059 = 4*Sf5056*Sf6058 + Sf6057*TB*USbot21*UStop21

	Sf6060 = (-6*Sf6059)/TB + (6*MW2*Sf4864*Sf6054)/MHiggs2(4)**2

	Sf5886 = 3*MW2*Sf4930*Sf5885 + Sf5883*TB*USbot11*UStop21

	Sf5759 = 3*MW2*Sf2364*Sf5758 + Sf5753*TB*USbot11*UStop11

	Sf10054 = ((MZ2*Sf10053)/Ss49**2 + 
     &     (8*Sf55**2*Sf57**2)/MHiggs2(2))/MZ2

	Sf5704 = MW2**3*Sf5695 - Sf5703*USbot11

	Sf5687 = Sf5685 + 33*Sf2553*MSbot2(1)*USbot11

	Sf5672 = Sf5670 - Sf5661*MSbot2(1)*USbot11

	Sf550 = Sf547 - Sf549*MSneu2(1)

	Sf5467 = 4*Sf4930*Sf5466 + Sf5463*TB*USbot11*UStop21

	Sf5468 = (-6*Sf5467)/TB + (6*MW2*Sf4828*Sf5460)/MHiggs2(4)**2

	Sf551 = Sf167**2*Sf542 - Sf550*MSneu2(1)

	Sf521 = Sf516*Sf518 - Sf520*MSneu2(1)**2

	Sf5057 = 3*Sf5052*Sf5053*Sf5054*Sf5055*Sf5056*
     &    MHiggs2(4)**2 + Sf5051*TB*USbot21*UStop21

	tmp92 = (3 + 2*C2B)*MSneu2(1) + 15*MStau2(1)

	tmp93 = 60*MStau2(1)**2 + 
     &   MSneu2(1)*(29*MSneu2(1) + 39*MStau2(1))

	Sf1005 = 3*(2*MW2**6 - MW2**5*tmp92 + 
     &     5*MSl2Diff(1,1)**6 - 
     &     MW2*MSl2Diff(1,1)**4*(20*MSneu2(1) + 27*MStau2(1)) - 
     &     MW2**2*MSl2Diff(1,1)*
     &      (-(tmp93*MSl2Diff(1,1)) + 
     &        ML2*MSneu2(1)*(2*MSneu2(1) + 2*MStau2(1))) + 
     &     MW2**4*(2*ML2*MSneu2(1) + 45*MStau2(1)**2 + 
     &        MSneu2(1)*((4 + 4*C2B)*MSneu2(1) + 7*MStau2(1))) - 
     &     MW2**3*(4*ML2*MSneu2(1)**2 + 70*MStau2(1)**3 + 
     &        MSneu2(1)*(MSneu2(1)*
     &            ((17 + 2*C2B)*MSneu2(1) - 4*MStau2(1)) - 
     &           19*MStau2(1)**2 - 2*C2B*MStau2(1)**2)))

	Sf4931 = 3*Sf4926*Sf4927*Sf4928*Sf4929*Sf4930*
     &    MHiggs2(4)**2 + Sf4925*TB*USbot11*UStop21

	Sf6403 = Sf4831*Sf4867*Sf4899*MHiggs2(4)**2 + 
     &   2*CW2*Sf6402*UStop21

	Sf6401 = Sf139*Sf2270*Sf4867*MHiggs2(4)**2 - 
     &   2*CW2*Sf6400*USbot21

	Sf6399 = Sf126*Sf2206*Sf4831*MHiggs2(4)**2 - 
     &   2*CW2*Sf6398*USbot11

	Sf4478 = Sf4475 - Sf4441*MSq2Diff(2,1)**2

	Sf4479 = Sf4450 - 2*Sf4430*MSq2Diff(2,1)**2

	Sf451 = Sf445 - MHiggs2(4)*MSneu2(1)*UStau11

	Sf4476 = 2*Sf4436 + MW2*Sf4457

	Sf4367 = Sf4364 - Sf4314*MSq2Diff(1,1)**2

	Sf4368 = Sf4324 - 2*Sf4301*MSq2Diff(1,1)**2

	Sf4365 = 2*Sf4308 + MW2*Sf4335

	Sf436 = 9*ML*Sf341*Sf429 + Sf435*UStau11

	tmp90 = 2*MBot2 - 2*MTop2 + (-3 + 4*C2B)*MSbot2(2)

	tmp91 = 12*MBot2 - 12*MTop2 + (-43 + 6*C2B)*MSbot2(2)

	Sf6290 = (4*MW2**6 + 12*MSq2Diff(2,2)**6 + 
     &      MW2**5*((-5 + 6*C2B)*MSbot2(2) - 32*MStop2(2)) - 
     &      16*MW2*MSq2Diff(2,2)**4*(3*MSbot2(2) + 4*MStop2(2)) + 
     &      MW2**4*(-3*tmp90*MSbot2(2) + 7*MSbot2(2)*MStop2(2) + 
     &         100*MStop2(2)**2) + 
     &      MW2**2*MSq2Diff(2,2)*
     &       (MSq2Diff(2,2)**2*
     &          (-6*MBot2 + 6*MTop2 + 71*MSq2Diff(2,2)) + 
     &         3*MSq2Diff(2,2)*(-6*MBot2 + 6*MTop2 + 77*MSq2Diff(2,2))*
     &          MStop2(2) + 
     &         12*(-MBot2 + MTop2 + 25*MSq2Diff(2,2))*MStop2(2)**2) + 
     &      MW2**3*(tmp91*MSbot2(2)**2 - 160*MStop2(2)**3 + 
     &         MSbot2(2)*(22*MSbot2(2)*MStop2(2) + 
     &            (53 - 6*C2B)*MStop2(2)**2)))*USbot21 + 
     &   6*MBot*MW2**2*(reMUE*TB - Re(Abot))*MSbot2(2)*
     &    ((MW2 - MSbot2(2))**2 - MStop2(2)**2)*USbot22

	Sf3988 = (4*MW2**6 + 12*MSq2Diff(2,1)**6 + 
     &      MW2**5*((-5 + 6*C2B)*MSbot2(2) - 32*MStop2(1)) - 
     &      16*MW2*MSq2Diff(2,1)**4*(3*MSbot2(2) + 4*MStop2(1)) + 
     &      MW2**4*(-3*tmp90*MSbot2(2) + 7*MSbot2(2)*MStop2(1) + 
     &         100*MStop2(1)**2) + 
     &      MW2**2*MSq2Diff(2,1)*
     &       (MSbot2(2)**2*(-6*MBot2 + 6*MTop2 + 71*MSbot2(2)) + 
     &         6*MSbot2(2)*(-MBot2 + MTop2 + 3*MSbot2(2))*MStop2(1) + 
     &         51*MSbot2(2)*MStop2(1)**2 - 140*MStop2(1)**3) + 
     &      MW2**3*(tmp91*MSbot2(2)**2 - 160*MStop2(1)**3 + 
     &         MSbot2(2)*(22*MSbot2(2)*MStop2(1) + 
     &            (53 - 6*C2B)*MStop2(1)**2)))*USbot21 + 
     &   6*MBot*MW2**2*(reMUE*TB - Re(Abot))*MSbot2(2)*
     &    ((MW2 - MSbot2(2))**2 - MStop2(1)**2)*USbot22

	Sf3968 = Sf3966 + 33*Sf2731*MSbot2(2)*USbot21

	Sf3962 = 3*Sf3960 - Sf3957*MSbot2(2)*USbot21

	Sf6283 = Sf3955 - Sf6282*MStop2(2)

	Sf3972 = Sf3955 - Sf3971*MStop2(1)

	Sf6284 = Sf3560*Sf3561*Sf3947 + Sf6283*MStop2(2)**2

	Sf3973 = Sf3560*Sf3561*Sf3947 + Sf3972*MStop2(1)**2

	Sf6276 = Sf2543**2*Sf3886 + Sf6275*MStop2(2)

	Sf3891 = (4*MW2**6 + 12*MSq2Diff(1,1)**6 + 
     &      MW2**5*((-5 + 6*C2B)*MSbot2(1) - 32*MStop2(1)) - 
     &      16*MW2*MSq2Diff(1,1)**4*(3*MSbot2(1) + 4*MStop2(1)) + 
     &      MW2**4*(-3*MSbot2(1)*
     &          (2*MBot2 - 2*MTop2 + (-3 + 4*C2B)*MSbot2(1)) + 
     &         7*MSbot2(1)*MStop2(1) + 100*MStop2(1)**2) + 
     &      MW2**3*(MSbot2(1)**2*
     &          (12*MBot2 - 12*MTop2 + (-43 + 6*C2B)*MSbot2(1)) + 
     &         22*MSbot2(1)**2*MStop2(1) - 
     &         (-53 + 6*C2B)*MSbot2(1)*MStop2(1)**2 - 
     &         160*MStop2(1)**3) + 
     &      MW2**2*MSq2Diff(1,1)*
     &       (MSbot2(1)**2*(-6*MBot2 + 6*MTop2 + 71*MSbot2(1)) + 
     &         6*MSbot2(1)*(-MBot2 + MTop2 + 3*MSbot2(1))*MStop2(1) + 
     &         51*MSbot2(1)*MStop2(1)**2 - 140*MStop2(1)**3))*
     &    USbot11 + 6*MBot*MW2**2*(reMUE*TB - Re(Abot))*
     &    MSbot2(1)*((MW2 - MSbot2(1))**2 - MStop2(1)**2)*
     &    USbot12

	Sf3871 = Sf3869 + 33*Sf2553*MSbot2(1)*USbot11

	Sf3858 = 3*Sf3856 - Sf3848*MSbot2(1)*USbot11

	Sf6269 = Sf3842 - Sf6268*MStop2(2)

	Sf3875 = Sf3842 - Sf3874*MStop2(1)

	Sf6270 = Sf3204*Sf3205*Sf3830 + Sf6269*MStop2(2)**2

	Sf3876 = Sf3204*Sf3205*Sf3830 + Sf3875*MStop2(1)**2

	Sf3728 = 4*Sf2602*Sf3727 + Sf3726*TB*USbot21*UStop11

	Sf3729 = (-6*Sf3728)/TB + (6*MW2*Sf2263*Sf3722)/MHiggs2(4)**2

	Sf6050 = Sf3713 - Sf6049*MStop2(2)

	Sf3718 = Sf3713 - Sf3717*MStop2(1)

	Sf3655 = 2*Sf3653 - Sf3648*MSbot2(2)*USbot21

	Sf3646 = 3*Sf3644 - 4*Sf2731*MSbot2(2)*USbot21

	Sf6014 = 5*Sf3628 - Sf6013*MStop2(2)

	Sf3658 = 5*Sf3628 - Sf3657*MStop2(1)

	Sf5984 = Sf3560*Sf3561*Sf3563 - Sf5983*MStop2(2)

	Sf3572 = Sf3560*Sf3561*Sf3563 - Sf3571*MStop2(1)

	tmp75 = 41*MBot2 - 41*MTop2 + 154*MW2 - 41*C2B*MW2

	tmp77 = 2*MBot2 - 2*MTop2 + 5*MW2 - 2*C2B*MW2

	tmp78 = -2*MBot2 + 2*MTop2 - 5*MW2 + 2*C2B*MW2

	tmp80 = -22*MBot2 + 22*MTop2 + 15*MW2 + 22*C2B*MW2

	tmp81 = -26*MBot2 + 26*MTop2 - 93*MW2 + 26*C2B*MW2

	Sf5985 = (MW2 - MSq2Diff(2,2))**3*
     &    (-2*MSq2Diff(2,2)**4 + 
     &      MW2*(16*MSq2Diff(2,2)**3 + 
     &         MW2*(-54*MSq2Diff(2,2)**2 + 
     &            MW2*(tmp77 + 33*MSq2Diff(2,2))))) + 
     &   3*MW2*(MW2 - MSq2Diff(2,2))*
     &    (14*MSq2Diff(2,2)**4 + 
     &      MW2*(-112*MSq2Diff(2,2)**3 + 
     &         MW2*(2*MW2*tmp78 + tmp81*MSq2Diff(2,2) + 
     &            231*MSq2Diff(2,2)**2)))*MStop2(2) + 
     &   3*MW2**2*(120*MSq2Diff(2,2)**3 + 
     &      MW2*(MW2*tmp80 + 2*tmp75*MSq2Diff(2,2) - 
     &         547*MSq2Diff(2,2)**2))*MStop2(2)**2 - 
     &   4*MW2**3*(-41*MBot2 + 41*MTop2 + 5*MW2 + 41*C2B*MW2 + 
     &      265*MSq2Diff(2,2))*MStop2(2)**3

	Sf3573 = Sf2722**3*Sf3559 - Sf3572*MStop2(1)

	Sf3457 = 4*Sf2364*Sf3456 + Sf3455*TB*USbot11*UStop11

	Sf3458 = (-6*Sf3457)/TB + (6*MW2*Sf2198*Sf3451)/MHiggs2(4)**2

	Sf5631 = Sf3441 - Sf5630*MStop2(2)

	Sf3446 = Sf3441 - Sf3445*MStop2(1)

	Sf3375 = 2*Sf3373 - Sf3359*MSbot2(1)*USbot11

	Sf3355 = 3*Sf3353 - 4*Sf2553*MSbot2(1)*USbot11

	Sf5594 = 5*Sf3315 - Sf5593*MStop2(2)

	Sf3378 = 5*Sf3315 - Sf3377*MStop2(1)

	Sf5558 = Sf3204*Sf3205*Sf3209 - Sf5557*MStop2(2)

	Sf3220 = Sf3204*Sf3205*Sf3209 - Sf3219*MStop2(1)

	Sf5559 = (MW2 - MSq2Diff(1,2))**3*
     &    (-2*MSq2Diff(1,2)**4 + 
     &      MW2*(16*MSq2Diff(1,2)**3 + 
     &         MW2*(-54*MSq2Diff(1,2)**2 + 
     &            MW2*(tmp77 + 33*MSq2Diff(1,2))))) + 
     &   3*MW2*(MW2 - MSq2Diff(1,2))*
     &    (14*MSq2Diff(1,2)**4 + 
     &      MW2*(-112*MSq2Diff(1,2)**3 + 
     &         MW2*(2*MW2*tmp78 + tmp81*MSq2Diff(1,2) + 
     &            231*MSq2Diff(1,2)**2)))*MStop2(2) + 
     &   3*MW2**2*(120*MSq2Diff(1,2)**3 + 
     &      MW2*(MW2*tmp80 + 2*tmp75*MSq2Diff(1,2) - 
     &         547*MSq2Diff(1,2)**2))*MStop2(2)**2 - 
     &   4*MW2**3*(-41*MBot2 + 41*MTop2 + 5*MW2 + 41*C2B*MW2 + 
     &      265*MSq2Diff(1,2))*MStop2(2)**3

	Sf3221 = Sf2543**3*Sf3203 - Sf3220*MStop2(1)

	tmp89 = -2*MHiggs2(4)**2 + MW2*(6*MW2 + 12*MHiggs2(4))

	Sf6051 = (MW2 - MSq2Diff(2,2))**3*
     &    (MHiggs2(4) - MSq2Diff(2,2)) + 
     &   (-2*MW2**2*(MW2 + 3*MHiggs2(4)) + tmp89*MSq2Diff(2,2) - 
     &      2*(3*MW2 + MHiggs2(4))*MSq2Diff(2,2)**2)*MStop2(2) + 
     &   tmp89*MStop2(2)**2

	Sf3719 = Sf2264*Sf2265*Sf2722**3 - Sf3718*MStop2(1)

	Sf267 = Sf265 - Sf266*MSneu2(1)

	Sf268 = Sf259*Sf260*Sf263 + Sf267*MSneu2(1)

	Sf2603 = 3*Sf2598*Sf2599*Sf2600*Sf2601*Sf2602*
     &    MHiggs2(4)**2 + Sf2597*TB*USbot21*UStop11

	Sf619 = (MW2 - MStau2(1))**3*(MHiggs2(4) - MStau2(1)) - 
     &   MSneu2(1)*(MW2**3 - (3*MW2 + MHiggs2(4))*MStau2(1)**2 + 
     &      4*MStau2(1)**3 + 
     &      MHiggs2(4)*(MW2*(3*MW2 - 6*MStau2(1)) + 
     &         2*MHiggs2(4)*MStau2(1)) - 
     &      MSneu2(1)*(3*MW2*(MW2 + MHiggs2(4)) + MSneu2(1)**2 - 
     &         MSneu2(1)*(3*MW2 + MHiggs2(4) + 4*MStau2(1)) + 
     &         MStau2(1)*(3*MW2 + MHiggs2(4) + 6*MStau2(1))))

	Sf5632 = Sf2199*Sf2200*Sf2543**3 - Sf5631*MStop2(2)

	Sf3447 = Sf2199*Sf2200*Sf2543**3 - Sf3446*MStop2(1)

	Sf2365 = 3*Sf2360*Sf2361*Sf2362*Sf2363*Sf2364*
     &    MHiggs2(4)**2 + Sf2359*TB*USbot11*UStop11

	Sf6397 = Sf2206*Sf2270*Sf2319*MHiggs2(4)**2 + 
     &   2*CW2*Sf6396*UStop11

	Sf4884 = Sf2293 + Sf4883*MStop2(2)

	Sf2295 = Sf2293 + Sf2294*MStop2(1)

	Sf4885 = Sf2264*Sf2265*Sf2290 + Sf4884*MStop2(2)

	Sf2296 = Sf2264*Sf2265*Sf2290 + Sf2295*MStop2(1)

	Sf4850 = Sf2241 + Sf4849*MStop2(2)

	Sf2243 = Sf2241 + Sf2242*MStop2(1)

	Sf4851 = Sf2199*Sf2200*Sf2234 + Sf4850*MStop2(2)

	Sf2244 = Sf2199*Sf2200*Sf2234 + Sf2243*MStop2(1)

	Sf1886 = Sf1885 - 2*Sf1875*Ss241*Ss242*UStau21

	Sf1871 = MW2*Sf1869 - Sf1870*MStau2(2)

	Sf1815 = 3*Sf1800*Sf1801*Sf1812 - Sf1814*MSneu2(1)

	Sf1659 = ML*Sf1293*Sf1651 + Sf1658*UStau21

	Sf1621 = 3*ML2*Sf1619 + Sf1620*MSl2Diff(2,1)

	Sf1513 = Sf1511 - Sf1512*MSneu2(1)

	Sf1453 = Sf1450 - Sf1452*MSneu2(1)

	Sf1454 = Sf1196**2*Sf1446 - Sf1453*MSneu2(1)

	Sf1432 = Sf1427*Sf1429 + Sf1431*MSneu2(1)**2

	tmp94 = (3 + 2*C2B)*MSneu2(1) + 15*MStau2(2)

	tmp95 = 60*MStau2(2)**2 + 
     &   MSneu2(1)*(29*MSneu2(1) + 39*MStau2(2))

	Sf1816 = 3*(2*MW2**6 - MW2**5*tmp94 + 
     &     5*MSl2Diff(2,1)**6 - 
     &     MW2*MSl2Diff(2,1)**4*(20*MSneu2(1) + 27*MStau2(2)) - 
     &     MW2**2*MSl2Diff(2,1)*
     &      (-(tmp95*MSl2Diff(2,1)) + 
     &        ML2*MSneu2(1)*(2*MSneu2(1) + 2*MStau2(2))) + 
     &     MW2**4*(2*ML2*MSneu2(1) + 45*MStau2(2)**2 + 
     &        MSneu2(1)*((4 + 4*C2B)*MSneu2(1) + 7*MStau2(2))) - 
     &     MW2**3*(4*ML2*MSneu2(1)**2 + 70*MStau2(2)**3 + 
     &        MSneu2(1)*(MSneu2(1)*
     &            ((17 + 2*C2B)*MSneu2(1) - 4*MStau2(2)) - 
     &           19*MStau2(2)**2 - 2*C2B*MStau2(2)**2)))

	Sf1368 = Sf1363 - MHiggs2(4)*MSneu2(1)*UStau21

	Sf1359 = 9*ML*Sf1293*Sf1352 - Sf1358*UStau21

	Sf1514 = (MW2 - MStau2(2))**3*(MHiggs2(4) - MStau2(2)) - 
     &   MSneu2(1)*(MW2**3 - (3*MW2 + MHiggs2(4))*MStau2(2)**2 + 
     &      4*MStau2(2)**3 + 
     &      MHiggs2(4)*(MW2*(3*MW2 - 6*MStau2(2)) + 
     &         2*MHiggs2(4)*MStau2(2)) - 
     &      MSneu2(1)*(3*MW2*(MW2 + MHiggs2(4)) + MSneu2(1)**2 - 
     &         MSneu2(1)*(3*MW2 + MHiggs2(4) + 4*MStau2(2)) + 
     &         MStau2(2)*(3*MW2 + MHiggs2(4) + 6*MStau2(2))))

	Sf1120 = MW2*Sf1118 - Sf1119*MStau2(1)

	Sf10550 = Sf10479*Sf10480*Sf10547 - Sf10549*MStop2(1)

	Sf10513 = Sf10245*Sf10246*Sf10510 - Sf10512*MStop2(1)

	Sf10189 = Sf10099*Sf10100*Sf10186 - Sf10188*MStau2(1)

	Sf7211 = 3*Sf5231 + 
     &   SB*(Sf137*(2 - 2*Log(MSbot2(1)))*MSbot2(1) + 
     &      Sf148*(2 - 2*Log(MSbot2(2)))*MSbot2(2) + 
     &      Sf949*(2 - 2*Log(MStau2(1)))*MStau2(1) + 
     &      Sf1737*(2 - 2*Log(MStau2(2)))*MStau2(2)) + 
     &   CB*(MZ2*SB*Sf5*(2 - 2*Log(MSneu2(1)))*MSneu2(1) + 
     &      Sf4265*(-2 + 2*Log(MStop2(1)))*MStop2(1) + 
     &      Sf6982*(-2 + 2*Log(MStop2(2)))*MStop2(2))

	tmp96 = (MW2 - MSneu2(1))**2 - MStau2(1)**2

	Sf996 = 3*ML2**2*MW2**2*tmp96*MSneu2(1) + 
     &   3*C2B**2*MW2**4*tmp96*MSneu2(1) + 
     &   (MSl2Diff(1,1)**2 + 
     &      MW2*(MW2 - 2*MSneu2(1) - 2*MStau2(1)))*
     &    (4*MSl2Diff(1,1)**5 + 3*MW2**4*MStau2(1) - 
     &      MW2*MSl2Diff(1,1)**3*(11*MSneu2(1) + 15*MStau2(1)) + 
     &      MW2**3*(3*MSneu2(1)**2 - 13*MStau2(1)**2) + 
     &      MW2**2*MSl2Diff(1,1)*
     &       (21*MStau2(1)**2 + 
     &         MSneu2(1)*(10*MSneu2(1) + 7*MStau2(1)))) + 
     &   3*C2B*MW2*(-2*MW2**6 - MW2**2*tmp93*MSl2Diff(1,1)**2 - 
     &      5*MSl2Diff(1,1)**6 + 
     &      3*MW2**5*(MSneu2(1) + 5*MStau2(1)) + 
     &      MW2*MSl2Diff(1,1)**4*(20*MSneu2(1) + 27*MStau2(1)) - 
     &      MW2**4*(45*MStau2(1)**2 + 
     &         MSneu2(1)*(4*MSneu2(1) + 7*MStau2(1))) + 
     &      MW2**3*(70*MStau2(1)**3 + 
     &         MSneu2(1)*(MSneu2(1)*(17*MSneu2(1) - 4*MStau2(1)) - 
     &            19*MStau2(1)**2))) + 
     &   3*ML2*(2*MW2**6 - MW2**5*tmp92 + 
     &      MW2**2*tmp93*MSl2Diff(1,1)**2 + 5*MSl2Diff(1,1)**6 - 
     &      MW2*MSl2Diff(1,1)**4*(20*MSneu2(1) + 27*MStau2(1)) + 
     &      MW2**4*(45*MStau2(1)**2 + 
     &         MSneu2(1)*((4 + 4*C2B)*MSneu2(1) + 7*MStau2(1))) - 
     &      MW2**3*(70*MStau2(1)**3 + 
     &         MSneu2(1)*(MSneu2(1)*
     &             ((17 + 2*C2B)*MSneu2(1) - 4*MStau2(1)) - 
     &            19*MStau2(1)**2 - 2*C2B*MStau2(1)**2)))

	Sf845 = 3*MW2**2*Sf844 + Sf842*UStau11 + 
     &   (-54*C2B*MW2*(MW2**6 + 3*MSl2Diff(1,1)**6) + 
     &      MW2**3*(Sf840 + 
     &         MW2*(3*Sf839 + MW2*(3*Sf838 + 7*MW2*MSl2Diff(1,1))))
     &      )*UStau22**2

	Sf807 = ML2*Sf806 + 9*MSl2Diff(1,1)**6 + 
     &   MW2*(-3*Sf798*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(-3*Sf789 + 
     &            MW2*(MW2*(MW2*Sf741 - 6*Sf766) + 3*Sf801)) + 
     &         3*Sf777*MSl2Diff(1,1)**2))

	Sf6876 = 2*MBot2*Sf6875 + 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf6861*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(-Sf6868 + 
     &            MW2*(Sf6863 + MW2*(MW2*Sf6479 - Sf6870))) + 
     &         Sf6865*MSq2Diff(2,2)**2))

	Sf6762 = 2*MBot2*Sf6761 + 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf6747*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(-Sf6754 + 
     &            MW2*(Sf6749 + MW2*(MW2*Sf6479 - Sf6756))) + 
     &         Sf6751*MSq2Diff(2,1)**2))

	Sf6625 = 2*MBot2*Sf6624 + 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf6607*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(-Sf6616 + 
     &            MW2*(Sf6610 + MW2*(MW2*Sf6479 - Sf6619))) + 
     &         Sf6612*MSq2Diff(1,2)**2))

	Sf6507 = 2*MBot2*Sf6506 + 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf6483*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(-Sf6499 + 
     &            MW2*(Sf6488 + MW2*(MW2*Sf6479 - Sf6501))) + 
     &         Sf6493*MSq2Diff(1,1)**2))

	Sf5404 = Sf5395 + 2*MTop2*Sf5400

	Sf5401 = Sf5395 - 2*MBot2*Sf5400

	Sf5298 = Sf5289 + 2*MTop2*Sf5294

	Sf5295 = Sf5289 - 2*MBot2*Sf5294

	Sf2996 = Sf2987 + 2*MTop2*Sf2992

	Sf2993 = Sf2987 - 2*MBot2*Sf2992

	Sf2858 = Sf2849 + 2*MTop2*Sf2854

	Sf2855 = Sf2849 - 2*MBot2*Sf2854

	Sf314 = MW2*(MW2*
     &       (MW2*(MW2*(2*MW2*Sf250 - Sf254) + Sf268) - Sf283) - 
     &      Sf306) + 2*MSl2Diff(1,1)**6 + 
     &   MHiggs2(4)*(-(Sf251*Sf252*MHiggs2(4)*MSl2Diff(1,1)**2) + 
     &      Sf34*MSl2Diff(1,1)**4)

	tmp97 = (MW2 - MSneu2(1))**2 - MStau2(2)**2

	Sf1808 = 3*ML2**2*MW2**2*tmp97*MSneu2(1) + 
     &   3*C2B**2*MW2**4*tmp97*MSneu2(1) + 
     &   (MSl2Diff(2,1)**2 + 
     &      MW2*(MW2 - 2*MSneu2(1) - 2*MStau2(2)))*
     &    (4*MSl2Diff(2,1)**5 + 3*MW2**4*MStau2(2) - 
     &      MW2*MSl2Diff(2,1)**3*(11*MSneu2(1) + 15*MStau2(2)) + 
     &      MW2**3*(3*MSneu2(1)**2 - 13*MStau2(2)**2) + 
     &      MW2**2*MSl2Diff(2,1)*
     &       (21*MStau2(2)**2 + 
     &         MSneu2(1)*(10*MSneu2(1) + 7*MStau2(2)))) + 
     &   3*C2B*MW2*(-2*MW2**6 - MW2**2*tmp95*MSl2Diff(2,1)**2 - 
     &      5*MSl2Diff(2,1)**6 + 
     &      3*MW2**5*(MSneu2(1) + 5*MStau2(2)) + 
     &      MW2*MSl2Diff(2,1)**4*(20*MSneu2(1) + 27*MStau2(2)) - 
     &      MW2**4*(45*MStau2(2)**2 + 
     &         MSneu2(1)*(4*MSneu2(1) + 7*MStau2(2))) + 
     &      MW2**3*(70*MStau2(2)**3 + 
     &         MSneu2(1)*(MSneu2(1)*(17*MSneu2(1) - 4*MStau2(2)) - 
     &            19*MStau2(2)**2))) + 
     &   3*ML2*(2*MW2**6 - MW2**5*tmp94 + 
     &      MW2**2*tmp95*MSl2Diff(2,1)**2 + 5*MSl2Diff(2,1)**6 - 
     &      MW2*MSl2Diff(2,1)**4*(20*MSneu2(1) + 27*MStau2(2)) + 
     &      MW2**4*(45*MStau2(2)**2 + 
     &         MSneu2(1)*((4 + 4*C2B)*MSneu2(1) + 7*MStau2(2))) - 
     &      MW2**3*(70*MStau2(2)**3 + 
     &         MSneu2(1)*(MSneu2(1)*
     &             ((17 + 2*C2B)*MSneu2(1) - 4*MStau2(2)) - 
     &            19*MStau2(2)**2 - 2*C2B*MStau2(2)**2)))

	Sf1662 = 3*MW2**2*Sf1661 + Sf1659*UStau21 - 
     &   MW2*(54*C2B*(MW2**6 + 3*MSl2Diff(2,1)**6) - 
     &      MW2**2*(Sf1657 + 
     &         MW2*(3*Sf1656 + 
     &            MW2*(-3*Sf1655 + 7*MW2*MSl2Diff(2,1)))))*
     &    UStau21**2

	Sf1632 = -(ML2*Sf1631) - 9*MSl2Diff(2,1)**6 + 
     &   MW2*(3*Sf1624*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(3*Sf1619 + 
     &            MW2*(-3*Sf1626 + MW2*(6*Sf1609 - MW2*Sf741))) - 
     &         3*Sf1614*MSl2Diff(2,1)**2))

	Sf1274 = MW2*(MW2*
     &       (MW2*(MW2*(2*MW2*Sf1232 - Sf1236) + Sf1244) - Sf1253)-
     &        Sf1266) + 2*MSl2Diff(2,1)**6 + 
     &   MHiggs2(4)*(-(Sf1233*Sf1234*MHiggs2(4)*
     &         MSl2Diff(2,1)**2) + Sf47*MSl2Diff(2,1)**4)

	Sf6006 = 30*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf6005) - Sf5998*MHiggs2(4)) + 
     &   MW2**2*(-(Sf5933*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf5990) + Sf5992) - Sf5995) + 
     &         Sf5935*MSq2Diff(2,2)**2))

	Sf5586 = 30*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf5585) - Sf5579*MHiggs2(4)) + 
     &   MW2**2*(-(Sf5500*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf5564) + Sf5566) - Sf5569) + 
     &         Sf5502*MSq2Diff(1,2)**2))

	Sf3603 = 30*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf3602) - Sf3593*MHiggs2(4)) + 
     &   MW2**2*(-(Sf3490*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf3580) + Sf3584) - Sf3590) + 
     &         Sf3494*MSq2Diff(2,1)**2))

	Sf3276 = 30*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf3275) - Sf3259*MHiggs2(4)) + 
     &   MW2**2*(-(Sf3107*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf3230) + Sf3237) - Sf3246) + 
     &         Sf3115*MSq2Diff(1,1)**2))

	Sf5950 = 6*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf5949) + Sf5940*MHiggs2(4)) + 
     &   MW2**2*(-(Sf5933*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf5934) + Sf5936) - Sf5939) + 
     &         Sf5935*MSq2Diff(2,2)**2))

	Sf5524 = 6*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf5523) + Sf5515*MHiggs2(4)) + 
     &   MW2**2*(-(Sf5500*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf5501) + Sf5503) - Sf5506) + 
     &         Sf5502*MSq2Diff(1,2)**2))

	Sf3518 = 6*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf3517) + Sf3504*MHiggs2(4)) + 
     &   MW2**2*(-(Sf3490*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf3491) + Sf3497) - Sf3503) + 
     &         Sf3494*MSq2Diff(2,1)**2))

	Sf3156 = 6*C2B*MW2**7 + 
     &   MHiggs2(4)*(-(MW2*Sf3155) + Sf3137*MHiggs2(4)) + 
     &   MW2**2*(-(Sf3107*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(-(MW2*Sf3110) + Sf3120) - Sf3127) + 
     &         Sf3115*MSq2Diff(1,1)**2))

	Sf700 = MHiggs2(4)*
     &    (-(MW2*Sf699) + MW2*MHiggs2(4)*
     &       (MW2*Sf697 + Sf696*MHiggs2(4))) - 
     &   MW2**3*Sf679*MSl2Diff(1,1)**2 - 8*Sf698*MSl2Diff(1,1)**4

	Sf6223 = MW2**3*(MW2*(4*MW2*Sf6213 + 4*Sf6215) + 
     &      Sf6222) - Sf6217*MSq2Diff(2,2)**3 + 
     &   6*Sf6214*MSq2Diff(2,2)**5

	Sf5846 = MW2**3*(MW2*(4*MW2*Sf5836 + 4*Sf5838) + 
     &      Sf5845) - Sf5840*MSq2Diff(1,2)**3 + 
     &   6*Sf5837*MSq2Diff(1,2)**5

	Sf4164 = MW2**3*(MW2*(4*MW2*Sf4152 + 4*Sf4154) + 
     &      Sf4163) - Sf4158*MSq2Diff(2,1)**3 + 
     &   6*Sf4153*MSq2Diff(2,1)**5

	Sf4062 = MW2**3*(MW2*(4*MW2*Sf4041 + 4*Sf4043) + 
     &      Sf4061) - Sf4053*MSq2Diff(1,1)**3 + 
     &   6*Sf4042*MSq2Diff(1,1)**5

	Sf1568 = MHiggs2(4)*
     &    (-(MW2*Sf1567) + 
     &      MW2*MHiggs2(4)*(MW2*Sf1565 + Sf1564*MHiggs2(4))) - 
     &   MW2**3*Sf1552*MSl2Diff(2,1)**2 - 8*Sf1566*MSl2Diff(2,1)**4

	Sf6405 = Sf2206*Sf4831*Sf6401*Ss906*MSbot2(2) + 
     &   Sf4831*Sf4867*Sf6397*Ss969*MStop2(1) + 
     &   Sf2270*(Sf4867*Sf6399*Ss885*MSbot2(1) + 
     &      Sf2206*Sf6403*Ss988*MStop2(2))

	Sf7149 = MW2*Sf2263*Sf6142 + 
     &   (-(Sf2340*Sf6110*Ss25**3) + 
     &      MW2*Sf2263*Sf2334**2*Ss25*TB*MHiggs2(4))*USbot21*
     &    UStop11 + MW2*Sf2334**2*Sf4864*Ss25*TB*MHiggs2(4)*
     &    USbot21*UStop21

	Sf7144 = MW2*Sf2198*Sf5759 + 
     &   (-(Sf2327*Sf5721*Ss25**3) + 
     &      MW2*Sf2198*Sf2320**2*Ss25*TB*MHiggs2(4))*USbot11*
     &    UStop11 + MW2*Sf2320**2*Sf4828*Ss25*TB*MHiggs2(4)*
     &    USbot11*UStop21

	Sf6085 = MW2**3*Sf3653 + MW2**2*Sf6084*MSbot2(2) + 
     &   2*MSbot2(2)**3*(-(MW2*Sf6079) + 2*Sf6083*MSbot2(2))*
     &    USbot21

	Sf5688 = MW2**3*Sf3373 + MW2**2*Sf5687*MSbot2(1) + 
     &   2*MSbot2(1)**3*(-(MW2*Sf5676) + 2*Sf5686*MSbot2(1))*
     &    USbot11

	tmp76 = 6*MBot2 - 6*MTop2 + 13*MW2 - 6*C2B*MW2

	tmp79 = -3*MBot2 + 3*MTop2 + 14*MW2 + 3*C2B*MW2

	Sf6022 = (6*MSq2Diff(2,2)**2*(-MW2 + MSq2Diff(2,2))**5 - 
     &      MW2*(MW2 - MSq2Diff(2,2))*
     &       (2*MW2**4 + (11 + 6*C2B)*MW2**3*MSq2Diff(2,2) + 
     &         MW2**2*(-6*MBot2 + 6*MTop2 - 73*MSq2Diff(2,2))*
     &          MSq2Diff(2,2) + 128*MW2*MSq2Diff(2,2)**3 - 
     &         62*MSq2Diff(2,2)**4)*MStop2(2) + 
     &      MW2**2*(200*MSq2Diff(2,2)**3 + 
     &         MW2*(MW2*tmp76 + 6*tmp79*MSq2Diff(2,2) - 
     &            241*MSq2Diff(2,2)**2))*MStop2(2)**2 + 
     &      4*MW2**3*(-3*MBot2 + 3*MTop2 - 5*MW2 + 3*C2B*MW2 - 
     &         45*MSq2Diff(2,2))*MStop2(2)**3)*USbot21 + 
     &   6*MBot*MW2**3*(reMUE*TB - Re(Abot))*MSbot2(2)*
     &    MStop2(2)*(-MW2 + MSbot2(2) + MStop2(2))*USbot22

	Sf5602 = (6*MSq2Diff(1,2)**2*(-MW2 + MSq2Diff(1,2))**5 - 
     &      MW2*(MW2 - MSq2Diff(1,2))*
     &       (2*MW2**4 + (11 + 6*C2B)*MW2**3*MSq2Diff(1,2) + 
     &         MW2**2*(-6*MBot2 + 6*MTop2 - 73*MSq2Diff(1,2))*
     &          MSq2Diff(1,2) + 128*MW2*MSq2Diff(1,2)**3 - 
     &         62*MSq2Diff(1,2)**4)*MStop2(2) + 
     &      MW2**2*(200*MSq2Diff(1,2)**3 + 
     &         MW2*(MW2*tmp76 + 6*tmp79*MSq2Diff(1,2) - 
     &            241*MSq2Diff(1,2)**2))*MStop2(2)**2 + 
     &      4*MW2**3*(-3*MBot2 + 3*MTop2 - 5*MW2 + 3*C2B*MW2 - 
     &         45*MSq2Diff(1,2))*MStop2(2)**3)*USbot11 + 
     &   6*MBot*MW2**3*(reMUE*TB - Re(Abot))*MSbot2(1)*
     &    MStop2(2)*(-MW2 + MSbot2(1) + MStop2(2))*USbot12

	tmp73 = 6*MBot2 - 6*MTop2 + 61*MStop2(1)

	tmp74 = MBot2 - MTop2 + (1 + C2B)*MStop2(1)

	Sf3676 = (6*MSq2Diff(2,1)**7 - 2*MW2**6*MStop2(1) - 
     &      2*MW2*MSq2Diff(2,1)**5*(15*MSbot2(2) + 16*MStop2(1)) + 
     &      MW2**5*(16*MStop2(1)**2 + 
     &         MSbot2(2)*(-6*MSbot2(2) + (3 - 6*C2B)*MStop2(1))) + 
     &      MW2**4*(-50*MStop2(1)**3 + 
     &         MSbot2(2)*(6*tmp74*MStop2(1) + 
     &            MSbot2(2)*(30*MSbot2(2) + (-6 + 6*C2B)*MStop2(1))
     &            )) + MW2**3*
     &       (80*MStop2(1)**4 + 
     &         MSbot2(2)*(-(tmp73*MStop2(1)**2) + 
     &            MSbot2(2)*
     &             (-6*MBot2*MStop2(1) + 
     &               MStop2(1)*(6*MTop2 + 2*MStop2(1)) + 
     &               MSbot2(2)*(-60*MSbot2(2) + 39*MStop2(1))))) + 
     &      10*MW2**2*MSq2Diff(2,1)**3*
     &       (6*MSq2Diff(2,2)**2 + 7*MStop2(1)**2 + 
     &         7*MStop2(1)*MStop2(2) + 6*MStop2(2)**2 + 
     &         MSq2Diff(2,2)*(7*MStop2(1) + 12*MStop2(2))))*
     &    USbot21 + 6*MBot*MW2**3*(reMUE*TB - Re(Abot))*
     &    MSbot2(2)*MStop2(1)*(-MW2 + MSbot2(2) + MStop2(1))*
     &    USbot22

	Sf3396 = (6*MSq2Diff(1,1)**7 - 2*MW2**6*MStop2(1) - 
     &      2*MW2*MSq2Diff(1,1)**5*(15*MSbot2(1) + 16*MStop2(1)) + 
     &      MW2**5*(16*MStop2(1)**2 + 
     &         MSbot2(1)*(-6*MSbot2(1) + (3 - 6*C2B)*MStop2(1))) + 
     &      MW2**4*(-50*MStop2(1)**3 + 
     &         MSbot2(1)*(6*tmp74*MStop2(1) + 
     &            MSbot2(1)*(30*MSbot2(1) + (-6 + 6*C2B)*MStop2(1))
     &            )) + MW2**3*
     &       (80*MStop2(1)**4 + 
     &         MSbot2(1)*(-(tmp73*MStop2(1)**2) + 
     &            MSbot2(1)*
     &             (-6*MBot2*MStop2(1) + 
     &               MStop2(1)*(6*MTop2 + 2*MStop2(1)) + 
     &               MSbot2(1)*(-60*MSbot2(1) + 39*MStop2(1))))) + 
     &      10*MW2**2*MSq2Diff(1,1)**3*
     &       (6*MSq2Diff(1,2)**2 + 7*MStop2(1)**2 + 
     &         7*MStop2(1)*MStop2(2) + 6*MStop2(2)**2 + 
     &         MSq2Diff(1,2)*(7*MStop2(1) + 12*MStop2(2))))*
     &    USbot11 + 6*MBot*MW2**3*(reMUE*TB - Re(Abot))*
     &    MSbot2(1)*MStop2(1)*(-MW2 + MSbot2(1) + MStop2(1))*
     &    USbot12

	Sf9766 = -(Sf9344*Sf9362*Sf9765) + 
     &   MZ2*(Sf9344**2*Sf9747 - 
     &      Sf9362**2*Sf9709**2*Sf9740*MHiggs2(2)*MStau2(3)**3)

	Sf9737 = -(Sf9344*Sf9362*Sf9736) + 
     &   MZ2*(Sf9344**2*Sf9718 - 
     &      16*Sf9362**2*Sf9709*Sf9711*MHiggs2(2)*MStau2(3)**3)

	Sf9532 = Sf9317*Sf9331*Sf9531 + 
     &   MZ2*(Sf9317**2*Sf9489 - 
     &      Sf9331**2*Sf9418**2*Sf9480*MHiggs2(2)*MStau2(3)**3)

	Sf9256 = 1/2.D0*(Sf9240*Sf9255) + 
     &   Sf9242*Sf9243*MSbot2(1) - Sf9222*Sf9247*MSbot2(2)

	Sf8858 = Sf25**2*Sf8765*Sf8768 + 
     &   3*Sf60**2*Sf8783*Sf8786 + Sf112**2*Sf8813*Sf8816

	Sf537 = ML2*Sf512 - Sf536*MHiggs2(4) - MSl2Diff(1,1)**5

	Sf7119 = Sf4898*Ss25**3 + 
     &   MW2*(-3*MW2**2*Sf4899 + 3*Sf4899*MHiggs2(4)**2)

	Sf4745 = Sf2315*Ss25**3 + 
     &   MW2*(-3*MW2**2*Sf2319 + 3*Sf2319*MHiggs2(4)**2)

	Sf1827 = Sf1826 + Sf1823*MSneu2(1) + Sf1824*MStau2(2)

	Sf1442 = ML2*Sf1423 + Sf1441*MHiggs2(4) - MSl2Diff(2,1)**5

	Sf6926 = Sf138*Ss25**3 + 
     &   MW2*(-3*MW2**2*Sf139 + 3*Sf139*MHiggs2(4)**2)

	Sf6675 = Sf122*Ss25**3 + 
     &   MW2*(-3*MW2**2*Sf126 + 3*Sf126*MHiggs2(4)**2)

	Sf10551 = 1/2.D0*(Sf10535*Sf10550) + 
     &   Sf10537*Sf10538*MStop2(1) - Sf10542*Sf9017*MStop2(2)

	Sf1035 = Sf1034 + Sf1022*MSneu2(1) + Sf1028*MStau2(1)

	Sf10190 = 1/2.D0*(Sf10174*Sf10189) + 
     &   Sf10176*Sf10177*MStau2(1) - Sf10181*Sf9222*MStau2(2)

	Sf9591 = (Sf9590/
     &      (Sf9317*Sf9331**2*Ss49**2*MHiggs2(2)*MStau2(2)) + 
     &     (MZ2*(25*MStau2(1) + MStau2(2)))/MStau2(1))/MZ2

	Sf9795 = ((-2*Sf9794)/
     &      (Sf9344*Sf9362**2*Ss49**2*MHiggs2(2)*MStau2(1)) + 
     &     MZ2*(25 + MStau2(1)/MStau2(2)))/MZ2

	Sf5403 = -9*MBot2*Sf5395 + 9*MTop2*Sf5401 - 
     &   22*MSq2Diff(2,2)**7 + 
     &   MW2*(Sf5402 + MW2*
     &       (-3*Sf5377*MSq2Diff(2,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf5376) + 3*Sf5378) - 
     &               2*Sf5380) + 2*Sf5379*MSq2Diff(2,2))))

	Sf5297 = -9*MBot2*Sf5289 + 9*MTop2*Sf5295 - 
     &   22*MSq2Diff(1,2)**7 + 
     &   MW2*(Sf5296 + MW2*
     &       (-3*Sf5269*MSq2Diff(1,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf5268) + 3*Sf5270) - 
     &               2*Sf5272) + 2*Sf5271*MSq2Diff(1,2))))

	Sf2995 = -9*MBot2*Sf2987 + 9*MTop2*Sf2993 - 
     &   22*MSq2Diff(2,1)**7 + 
     &   MW2*(Sf2994 + MW2*
     &       (-3*Sf2969*MSq2Diff(2,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf2968) + 3*Sf2970) - 
     &               2*Sf2972) + 2*Sf2971*MSq2Diff(2,1))))

	Sf2857 = -9*MBot2*Sf2849 + 9*MTop2*Sf2855 - 
     &   22*MSq2Diff(1,1)**7 + 
     &   MW2*(Sf2856 + MW2*
     &       (-3*Sf2817*MSq2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf2816) + 3*Sf2818) - 
     &               2*Sf2820) + 2*Sf2819*MSq2Diff(1,1))))

	Sf6228 = MHiggs2(4)*
     &    (-2*MW2*Sf6223 + 
     &      MHiggs2(4)*(MW2**2*Sf6212 + 
     &         MHiggs2(4)*(2*Sf6227 + MW2*Sf6202*MHiggs2(4)))) - 
     &   2*MW2**3*Sf6206*MSq2Diff(2,2)**2 - 
     &   Sf6197*MSq2Diff(2,2)**4 - 4*Sf4007*MSq2Diff(2,2)**6

	Sf5862 = MHiggs2(4)*
     &    (-2*MW2*Sf5846 + 
     &      MHiggs2(4)*(MW2**2*Sf5835 + 
     &         MHiggs2(4)*(2*Sf5861 + MW2*Sf5824*MHiggs2(4)))) - 
     &   2*MW2**3*Sf5828*MSq2Diff(1,2)**2 - 
     &   Sf5819*MSq2Diff(1,2)**4 - 4*Sf4007*MSq2Diff(1,2)**6

	Sf4169 = MHiggs2(4)*
     &    (-2*MW2*Sf4164 + 
     &      MHiggs2(4)*(MW2**2*Sf4151 + 
     &         MHiggs2(4)*(2*Sf4168 + MW2*Sf4136*MHiggs2(4)))) - 
     &   2*MW2**3*Sf4142*MSq2Diff(2,1)**2 - 
     &   Sf4128*MSq2Diff(2,1)**4 - 4*Sf4007*MSq2Diff(2,1)**6

	Sf4084 = MHiggs2(4)*
     &    (-2*MW2*Sf4062 + 
     &      MHiggs2(4)*(MW2**2*Sf4040 + 
     &         MHiggs2(4)*(2*Sf4083 + MW2*Sf4017*MHiggs2(4)))) - 
     &   2*MW2**3*Sf4025*MSq2Diff(1,1)**2 - 
     &   Sf4008*MSq2Diff(1,1)**4 - 4*Sf4007*MSq2Diff(1,1)**6

	Sf1872 = MW2**2*Sf1863 - 
     &   MSneu2(1)*(Sf1871 + Sf1867*MSneu2(1)) + 
     &   MW2*(3*Sf1864*MStau2(2) - 5*MSneu2(1)**3*UStau21) - 
     &   MStau2(2)**2*(Sf1865 + Sf1090*MStau2(2)*UStau21)

	Sf1121 = MW2**2*Sf1094 - 
     &   MSneu2(1)*(Sf1120 + Sf1109*MSneu2(1)) + 
     &   MW2*(3*Sf1098*MStau2(1) - 5*MSneu2(1)**3*UStau11) - 
     &   MStau2(1)**2*(Sf1104 + Sf1090*MStau2(1)*UStau11)

	Sf7496 = -(MW2*MZ2*
     &      (-2*Sf7495 + 3*MW2*MZ2**2*
     &         (3*Sf2022*Sf7318*MSneu2(1) + 
     &           4*MZ2*SW2**2*
     &            (Sf7375*Ss885*MSbot2(1) + 
     &              Sf7391*Ss906*MSbot2(2) + 
     &              3*Sf7338*Ss928*MStau2(1) + 
     &              3*Sf7350*Ss949*MStau2(2) + 
     &              Sf7419*Ss969*MStop2(1) + Sf7435*Ss988*MStop2(2)
     &              )))) - 
     &   6*MW2**2*MZ2*Sf11*Sf7300**2*
     &    myB0(MZ2,MSbot2(1),MSbot2(1),1.D0) - 
     &   6*(MW2**2*(2*Sf21*Sf23*Sf7303**2*
     &          myB0(MZ2,MSbot2(1),MSbot2(2),1.D0) + 
     &         6*Sf56*Sf58*Sf7321**2*
     &          myB0(MZ2,MStau2(1),MStau2(2),1.D0) + 
     &         2*Sf108*Sf110*Sf7394**2*
     &          myB0(MZ2,MStop2(1),MStop2(2),1.D0) + 
     &         MZ2*(Sf26*Sf7306**2*
     &             myB0(MZ2,MSbot2(2),MSbot2(2),1.D0) + 
     &            3*MZ2**2*Sf1*myB0(MZ2,MSneu2(1),MSneu2(1),1.D0) + 
     &            3*Sf39*Sf7309**2*
     &             myB0(MZ2,MStau2(1),MStau2(1),1.D0) + 
     &            3*Sf61*Sf7324**2*
     &             myB0(MZ2,MStau2(2),MStau2(2),1.D0) + 
     &            Sf7353**2*Sf82*
     &             myB0(MZ2,MStop2(1),MStop2(1),1.D0) + 
     &            Sf113*Sf7397**2*myB0(MZ2,MStop2(2),MStop2(2),1.D0))
     &         ) + MZ2**3*(-6*Sf1195*Ss443*
     &          myB0(MW2,MSneu2(1),MStau2(2),1.D0)*UStau21**2 - 
     &         6*Sf166*Ss443*myB0(MW2,MSneu2(1),MStau2(1),1.D0)*
     &          UStau22**2)) + 
     &   108*MZ2**3*Ss443*((Sf4906*
     &          myB0(MW2,MStop2(2),MSbot2(2),1.D0)*USbot21**2 + 
     &         Sf4820*myB0(MW2,MSbot2(1),MStop2(2),1.D0)*
     &          USbot22**2)*UStop21**2 + 
     &      (Sf2334*myB0(MW2,MStop2(1),MSbot2(2),1.D0)*
     &          USbot21**2 + 
     &         Sf2320*myB0(MW2,MStop2(1),MSbot2(1),1.D0)*
     &          USbot22**2)*UStop22**2)

	Sf214 = Sf194*MSneu2(1) - Sf205*MSneu2(1)**2 + 
     &   Sf213*MSneu2(1)**3 - Sf183*MSneu2(1)**4 + 
     &   MW2**2*MSneu2(1)**5 + Sf167**2*Sf179*MStau2(1)

	Sf1210 = MW2**2*MSneu2(1)**5 - 
     &   MSneu2(1)**4*((ML2 + MW2 - C2B*MW2)*MHiggs2(4)*
     &       (-3*MW2 + MHiggs2(4)) + 
     &      MW2**2*(2*ML2 + 6*MW2 - 2*C2B*MW2 + 5*MStau2(2))) + 
     &   (MW2 - MStau2(2))**2*MStau2(2)*
     &    (4*ML2*MW2**3 - MW2*(-3*ML2 + MW2 + 3*C2B*MW2)*
     &       MHiggs2(4)**2 + 
     &      (-2*MW2 + MHiggs2(4))*
     &       (MW2*(-MW2 - MHiggs2(4)) + ML2*(MW2 - MHiggs2(4)) + 
     &         C2B*MW2*(-MW2 + MHiggs2(4)))*MStau2(2) - 
     &      MW2**2*(C2B*MW2*(4*MW2 - 7*MHiggs2(4)) + 
     &         MW2*(3*MW2 - 3*MHiggs2(4)) + 7*ML2*MHiggs2(4) + 
     &         MStau2(2)**2)) + 
     &   MSneu2(1)**3*(2*MW2**2*
     &       ((7 - 4*C2B)*MW2**2 + (7 - 4*C2B)*MW2*MStau2(2) + 
     &         4*ML2*(MW2 + MStau2(2)) + 
     &         MStau2(2)*((-5 + 6*C2B)*MHiggs2(4) + 5*MStau2(2)))+
     &        MHiggs2(4)*(-12*ML2*MW2*MStau2(2) + 
     &         (ML2 + MW2 - C2B*MW2)*
     &          (MW2*(-13*MW2 + 5*MHiggs2(4)) + 
     &            4*MHiggs2(4)*MStau2(2)))) - 
     &   MSneu2(1)**2*(-6*(-MW2 + MHiggs2(4))*
     &       (MW2*(MW2 - MHiggs2(4)) + ML2*(2*MW2 - MHiggs2(4)) + 
     &         C2B*MW2*(-2*MW2 + MHiggs2(4)))*MStau2(2)**2 + 
     &      MW2*((2*MW2**2*(4*ML2 + MW2 - 4*C2B*MW2) + 
     &            MW2*(-13*ML2 + MW2 + 13*C2B*MW2)*MHiggs2(4) + 
     &            (5*ML2 - MW2 - 5*C2B*MW2)*MHiggs2(4)**2)*
     &          MStau2(2) + 
     &         MW2*(2*MW2**2*(5*ML2 + 7*MW2 - 5*C2B*MW2) + 
     &            (ML2 + MW2 - C2B*MW2)*MHiggs2(4)*
     &             (-17*MW2 + 7*MHiggs2(4)) + 10*MStau2(2)**3))) + 
     &   MSneu2(1)*(2*(MW2**2*(4*ML2 - 3*MW2 - 4*C2B*MW2) - 
     &         3*MW2*(2*ML2 + MW2 - 2*C2B*MW2)*MHiggs2(4) + 
     &         2*(ML2 + MW2 - C2B*MW2)*MHiggs2(4)**2)*MStau2(2)**3+
     &        MW2*((4*MW2**2*(-2*ML2 - MW2 + 2*C2B*MW2) + 
     &            13*MW2*(ML2 + MW2 - C2B*MW2)*MHiggs2(4) + 
     &            (-5*ML2 - 7*MW2 + 5*C2B*MW2)*MHiggs2(4)**2)*
     &          MStau2(2)**2 + 
     &         MW2*(5*MStau2(2)**4 + 
     &            MW2*(4*ML2*MW2**2 + (5 - 4*C2B)*MW2**3 + 
     &               (ML2 + MW2 - C2B*MW2)*MHiggs2(4)*
     &                (-7*MW2 + 3*MHiggs2(4)) - 
     &               2*(-MW2 + MHiggs2(4))*
     &                (-ML2 - MW2 + C2B*MW2 + MHiggs2(4))*MStau2(2)
     &               ))))

	Sf7172 = Sf122*MSbot2(1) + Sf138*MSbot2(2) - 
     &   Sf4*MSneu2(1) - Sf151*MStau2(1) - Sf1182*MStau2(2) + 
     &   Sf2315*MStop2(1) + Sf4898*MStop2(2)

	Sf6819 = MW2*(MW2**3*
     &      (MW2*(4*MW2*Sf6782 - 4*Sf6783) + Sf6785) - 
     &     Sf6818*MSq2Diff(2,2)**3)

	Sf6564 = MW2*(MW2**3*
     &      (MW2*(4*MW2*Sf6528 - 4*Sf6529) + Sf6531) - 
     &     Sf6563*MSq2Diff(1,2)**3)

	Sf4480 = MW2*(MW2**3*
     &      (MW2*(4*MW2*Sf4424 - 4*Sf4425) + Sf4429) - 
     &     Sf4479*MSq2Diff(2,1)**3)

	Sf4369 = MW2*(MW2**3*
     &      (MW2*(4*MW2*Sf4287 - 4*Sf4288) + Sf4300) - 
     &     Sf4368*MSq2Diff(1,1)**3)

	Sf797 = -(Sf796*MSl2Diff(1,1)**6) + 
     &   MW2*(Sf751*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(-Sf763 + 
     &               MW2*(MW2*(MW2*Sf740 - Sf744) + Sf772)) + Sf795
     &            ) - 3*Sf782*MSl2Diff(1,1)**2))

	Sf927 = -3*MSl2Diff(1,1)**7 + 
     &   MW2*(3*Sf899*MSl2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(4*MW2*Sf326 - 2*Sf913) - 3*Sf918) + 
     &               Sf925) - 2*Sf926) - 54*Sf901*MSl2Diff(1,1)**3)
     &      )

	tmp88 = MW2*(MW2*Sf3184 + 6*Sf5971) - 3*Sf5973

	Sf5981 = 2*MSq2Diff(2,2)**7 + 
     &   MW2*(-2*Sf5965*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(-Sf5980 + MW2*(2*Sf5976 + MW2*tmp88)) + 
     &         18*Sf5966*MSq2Diff(2,2)**3))

	Sf5977 = 2*MSq2Diff(2,2)**7 + 
     &   MW2*(-2*Sf5965*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(-Sf5969 + MW2*(2*Sf5976 + MW2*tmp88)) + 
     &         18*Sf5966*MSq2Diff(2,2)**3))

	tmp86 = MW2*(MW2*Sf3184 + 6*Sf5545) - 3*Sf5547

	Sf5555 = 2*MSq2Diff(1,2)**7 + 
     &   MW2*(-2*Sf5539*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(-Sf5554 + MW2*(2*Sf5550 + MW2*tmp86)) + 
     &         18*Sf5541*MSq2Diff(1,2)**3))

	Sf5551 = 2*MSq2Diff(1,2)**7 + 
     &   MW2*(-2*Sf5539*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(-Sf5543 + MW2*(2*Sf5550 + MW2*tmp86)) + 
     &         18*Sf5541*MSq2Diff(1,2)**3))

	tmp84 = MW2*(MW2*Sf3184 + 6*Sf3546) - 3*Sf3548

	Sf3556 = 2*MSq2Diff(2,1)**7 + 
     &   MW2*(-2*Sf3540*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(-Sf3555 + MW2*(2*Sf3551 + MW2*tmp84)) + 
     &         18*Sf3541*MSq2Diff(2,1)**3))

	Sf3552 = 2*MSq2Diff(2,1)**7 + 
     &   MW2*(-2*Sf3540*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(-Sf3544 + MW2*(2*Sf3551 + MW2*tmp84)) + 
     &         18*Sf3541*MSq2Diff(2,1)**3))

	Sf3200 = 2*MSq2Diff(1,1)**7 + 
     &   MW2*(-2*Sf3178*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(MW2*Sf3184 + 6*Sf3185) - 3*Sf3189) + 
     &               2*Sf3195) - Sf3199) + 
     &         18*Sf3179*MSq2Diff(1,1)**3))

	Sf3196 = -2*MSq2Diff(1,1)**7 + 
     &   MW2*(2*Sf3178*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(Sf3182 + 
     &            MW2*(MW2*
     &                (MW2*(-(MW2*Sf3184) - 6*Sf3185) + 3*Sf3189)-
     &                 2*Sf3195)) - 18*Sf3179*MSq2Diff(1,1)**3))

	tmp87 = MW2*(MW2*Sf184 + 2*Sf5915) - Sf5917

	Sf5929 = 2*MSq2Diff(2,2)**7 + 
     &   MW2*(-2*Sf5909*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(-Sf5928 + MW2*(2*Sf5920 + MW2*tmp87)) + 
     &         2*Sf5911*MSq2Diff(2,2)**3))

	Sf5925 = 2*MSq2Diff(2,2)**7 + 
     &   MW2*(-2*Sf5909*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(-Sf5924 + MW2*(2*Sf5920 + MW2*tmp87)) + 
     &         2*Sf5911*MSq2Diff(2,2)**3))

	Sf5921 = 2*MSq2Diff(2,2)**7 + 
     &   MW2*(-2*Sf5909*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(-Sf5913 + MW2*(2*Sf5920 + MW2*tmp87)) + 
     &         2*Sf5911*MSq2Diff(2,2)**3))

	tmp85 = MW2*(MW2*Sf184 + 2*Sf5481) - Sf5484

	Sf5496 = 2*MSq2Diff(1,2)**7 + 
     &   MW2*(-2*Sf5475*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(-Sf5495 + MW2*(2*Sf5487 + MW2*tmp85)) + 
     &         2*Sf5477*MSq2Diff(1,2)**3))

	Sf5492 = 2*MSq2Diff(1,2)**7 + 
     &   MW2*(-2*Sf5475*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(-Sf5491 + MW2*(2*Sf5487 + MW2*tmp85)) + 
     &         2*Sf5477*MSq2Diff(1,2)**3))

	Sf5488 = 2*MSq2Diff(1,2)**7 + 
     &   MW2*(-2*Sf5475*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(-Sf5479 + MW2*(2*Sf5487 + MW2*tmp85)) + 
     &         2*Sf5477*MSq2Diff(1,2)**3))

	tmp83 = MW2*(MW2*Sf184 + 2*Sf3471) - Sf3473

	Sf3485 = 2*MSq2Diff(2,1)**7 + 
     &   MW2*(-2*Sf3465*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(-Sf3484 + MW2*(2*Sf3476 + MW2*tmp83)) + 
     &         2*Sf3467*MSq2Diff(2,1)**3))

	Sf3481 = 2*MSq2Diff(2,1)**7 + 
     &   MW2*(-2*Sf3465*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(-Sf3480 + MW2*(2*Sf3476 + MW2*tmp83)) + 
     &         2*Sf3467*MSq2Diff(2,1)**3))

	Sf3477 = 2*MSq2Diff(2,1)**7 + 
     &   MW2*(-2*Sf3465*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(-Sf3469 + MW2*(2*Sf3476 + MW2*tmp83)) + 
     &         2*Sf3467*MSq2Diff(2,1)**3))

	tmp82 = MW2*(MW2*Sf184 + 2*Sf3085) - Sf3088

	Sf3102 = 2*MSq2Diff(1,1)**7 + 
     &   MW2*(-2*Sf3078*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(-Sf3101 + MW2*(2*Sf3093 + MW2*tmp82)) + 
     &         2*Sf3080*MSq2Diff(1,1)**3))

	Sf3098 = 2*MSq2Diff(1,1)**7 + 
     &   MW2*(-2*Sf3078*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(-Sf3097 + MW2*(2*Sf3093 + MW2*tmp82)) + 
     &         2*Sf3080*MSq2Diff(1,1)**3))

	Sf3094 = 2*MSq2Diff(1,1)**7 + 
     &   MW2*(-2*Sf3078*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(-Sf3082 + MW2*(2*Sf3093 + MW2*tmp82)) + 
     &         2*Sf3080*MSq2Diff(1,1)**3))

	Sf1717 = 3*MSl2Diff(2,1)**7 + 
     &   MW2*(-3*Sf1698*MSl2Diff(2,1)**5 + 
     &      MW2*(MW2*(2*Sf1716 + 
     &            MW2*(-Sf1715 + 
     &               MW2*(3*Sf1712 + MW2*(2*Sf1710 - 4*MW2*Sf326)))
     &            ) + 54*Sf1700*MSl2Diff(2,1)**3))

	Sf1623 = Sf1622*MSl2Diff(2,1)**6 + 
     &   MW2*(-(Sf1603*MSl2Diff(2,1)**4) + 
     &      MW2*(MW2*(-Sf1621 + 
     &            MW2*(Sf1608 + 
     &               MW2*(-Sf1612 + MW2*(Sf1600 - MW2*Sf740)))) + 
     &         3*Sf1616*MSl2Diff(2,1)**2))

	Sf6869 = 12*MSq2Diff(2,2)**6 + 
     &   MW2*(-3*Sf6861*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(3*MW2*Sf6479 - Sf6860) + 
     &               3*Sf6863) - 3*Sf6868) + 
     &         3*Sf6865*MSq2Diff(2,2)**2))

	Sf6755 = 12*MSq2Diff(2,1)**6 + 
     &   MW2*(-3*Sf6747*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(3*MW2*Sf6479 - Sf6746) + 
     &               3*Sf6749) - 3*Sf6754) + 
     &         3*Sf6751*MSq2Diff(2,1)**2))

	Sf6617 = 12*MSq2Diff(1,2)**6 + 
     &   MW2*(-3*Sf6607*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(3*MW2*Sf6479 - Sf6606) + 
     &               3*Sf6610) - 3*Sf6616) + 
     &         3*Sf6612*MSq2Diff(1,2)**2))

	Sf6500 = 12*MSq2Diff(1,1)**6 + 
     &   MW2*(-3*Sf6483*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(3*MW2*Sf6479 - Sf6481) + 
     &               3*Sf6488) - 3*Sf6499) + 
     &         3*Sf6493*MSq2Diff(1,1)**2))

	Sf6048 = MHiggs2(4)*
     &    (-(MW2*Sf6047) + 2*Sf6039*MHiggs2(4)) + 
     &   MW2**2*(-6*Sf6025*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf6026 + 2*Sf6028) - 6*Sf6033) + 
     &         2*Sf6030*MSq2Diff(2,2)**2))

	Sf5629 = MHiggs2(4)*
     &    (-(MW2*Sf5628) + 2*Sf5619*MHiggs2(4)) + 
     &   MW2**2*(-6*Sf5605*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf5606 + 2*Sf5608) - 6*Sf5613) + 
     &         2*Sf5610*MSq2Diff(1,2)**2))

	Sf7096 = 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf7088*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf4540 - Sf7087) + Sf7092) - 
     &            Sf7095) + Sf7090*MSq2Diff(2,2)**2))

	Sf7056 = 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf7048*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf4540 - Sf7047) + Sf7052) - 
     &            Sf7055) + Sf7050*MSq2Diff(1,2)**2))

	Sf4665 = 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf4657*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf4540 - Sf4656) + Sf4661) - 
     &            Sf4664) + Sf4659*MSq2Diff(2,1)**2))

	Sf4560 = 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf4544*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf4540 - Sf4541) + Sf4553) - 
     &            Sf4559) + Sf4548*MSq2Diff(1,1)**2))

	Sf7030 = MHiggs2(4)*
     &    (-(MW2*Sf7029) + Sf7023*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf4376 - Sf7013) + Sf7017) + 
     &         Sf7015*MSq2Diff(2,2)))

	Sf7008 = MHiggs2(4)*
     &    (-(MW2*Sf7007) + Sf7001*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf4376 - Sf6991) + Sf6995) + 
     &         Sf6993*MSq2Diff(1,2)))

	Sf4517 = MHiggs2(4)*
     &    (-(MW2*Sf4516) + Sf4506*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf4376 - Sf4487) + Sf4495) + 
     &         Sf4491*MSq2Diff(2,1)))

	Sf4413 = MHiggs2(4)*
     &    (-(MW2*Sf4412) + Sf4398*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf4376 - Sf4378) + Sf4386) + 
     &         Sf4382*MSq2Diff(1,1)))

	Sf3712 = MHiggs2(4)*
     &    (-(MW2*Sf3711) + 2*Sf3699*MHiggs2(4)) + 
     &   MW2**2*(-6*Sf3679*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf3680 + 2*Sf3684) - 6*Sf3692) + 
     &         2*Sf3688*MSq2Diff(2,1)**2))

	Sf387 = 3*MSl2Diff(1,1)**6 + 
     &   MW2*(-(Sf367*MSl2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf364 - 2*Sf370) + Sf380) - 
     &            Sf386) + Sf375*MSl2Diff(1,1)**2))

	Sf3440 = MHiggs2(4)*
     &    (-(MW2*Sf3439) + 2*Sf3421*MHiggs2(4)) + 
     &   MW2**2*(-6*Sf3399*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf3401 + 2*Sf3405) - 6*Sf3414) + 
     &         2*Sf3410*MSq2Diff(1,1)**2))

	Sf5361 = 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf5352*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2775 - Sf5353) + Sf5357) - 
     &            Sf5360) + Sf5355*MSq2Diff(2,2)**2))

	Sf5253 = 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf5244*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2775 - Sf5245) + Sf5249) - 
     &            Sf5252) + Sf5247*MSq2Diff(1,2)**2))

	Sf2950 = 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf2941*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2775 - Sf2942) + Sf2946) - 
     &            Sf2949) + Sf2944*MSq2Diff(2,1)**2))

	Sf2797 = 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf2778*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2775 - Sf2781) + Sf2790) - 
     &            Sf2796) + Sf2785*MSq2Diff(1,1)**2))

	Sf5134 = MHiggs2(4)*
     &    (-(MW2*Sf5133) + Sf5127*MHiggs2(4)) + 
     &   MW2**2*(Sf5093*MSq2Diff(2,2)**3 + 
     &      MW2*(MW2*(MW2*(3*MW2*Sf2499 - Sf5117) + Sf5121) - 
     &         Sf5119*MSq2Diff(2,2)))

	Sf5008 = MHiggs2(4)*
     &    (-(MW2*Sf5007) + Sf5001*MHiggs2(4)) + 
     &   MW2**2*(Sf4967*MSq2Diff(1,2)**3 + 
     &      MW2*(MW2*(MW2*(3*MW2*Sf2499 - Sf4991) + Sf4995) - 
     &         Sf4993*MSq2Diff(1,2)))

	Sf2721 = MHiggs2(4)*
     &    (-(MW2*Sf2720) + Sf2709*MHiggs2(4)) + 
     &   MW2**2*(Sf2653*MSq2Diff(2,1)**3 + 
     &      MW2*(MW2*(MW2*(3*MW2*Sf2499 - Sf2692) + Sf2700) - 
     &         Sf2696*MSq2Diff(2,1)))

	Sf2542 = MHiggs2(4)*
     &    (-(MW2*Sf2541) + Sf2524*MHiggs2(4)) + 
     &   MW2**2*(Sf2440*MSq2Diff(1,1)**3 + 
     &      MW2*(MW2*(MW2*(3*MW2*Sf2499 - Sf2501) + Sf2512) - 
     &         Sf2506*MSq2Diff(1,1)))

	Sf5112 = MHiggs2(4)*
     &    (-(MW2*Sf5111) + Sf5105*MHiggs2(4)) + 
     &   MW2**2*(Sf5093*MSq2Diff(2,2)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2439 + Sf5094) - Sf5096) - 
     &         Sf5098*MSq2Diff(2,2)))

	Sf4986 = MHiggs2(4)*
     &    (-(MW2*Sf4985) + Sf4979*MHiggs2(4)) + 
     &   MW2**2*(Sf4967*MSq2Diff(1,2)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2439 + Sf4968) - Sf4970) - 
     &         Sf4972*MSq2Diff(1,2)))

	Sf2683 = MHiggs2(4)*
     &    (-(MW2*Sf2682) + Sf2672*MHiggs2(4)) + 
     &   MW2**2*(Sf2653*MSq2Diff(2,1)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2439 + Sf2654) - Sf2658) - 
     &         Sf2662*MSq2Diff(2,1)))

	Sf2490 = MHiggs2(4)*
     &    (-(MW2*Sf2489) + Sf2471*MHiggs2(4)) + 
     &   MW2**2*(Sf2440*MSq2Diff(1,1)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2439 + Sf2443) - Sf2449) - 
     &         Sf2456*MSq2Diff(1,1)))

	Sf6343 = MHiggs2(4)*
     &    (-(MW2*Sf6342) + Sf6336*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf6025*MSq2Diff(2,2)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2371 - Sf6326) + Sf6330) + 
     &         Sf6328*MSq2Diff(2,2)))

	Sf6320 = MHiggs2(4)*
     &    (-(MW2*Sf6319) + Sf6313*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf5605*MSq2Diff(1,2)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2371 - Sf6301) + Sf6305) + 
     &         Sf6303*MSq2Diff(1,2)))

	Sf4199 = MHiggs2(4)*
     &    (-(MW2*Sf4198) + Sf4187*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf3679*MSq2Diff(2,1)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2371 - Sf4172) + Sf4179) + 
     &         Sf4176*MSq2Diff(2,1)))

	Sf4119 = MHiggs2(4)*
     &    (-(MW2*Sf4118) + Sf4103*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf3399*MSq2Diff(1,1)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2371 - Sf4087) + Sf4095) + 
     &         Sf4091*MSq2Diff(1,1)))

	Sf6843 = MHiggs2(4)*
     &    (-(MW2*Sf6842) + Sf6835*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2350 + Sf6824) - Sf6828) + 
     &         Sf6826*MSq2Diff(2,2)))

	Sf6725 = MHiggs2(4)*
     &    (-(MW2*Sf6724) + Sf6714*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2350 + Sf6696) - Sf6704) + 
     &         Sf6700*MSq2Diff(2,1)))

	Sf6588 = MHiggs2(4)*
     &    (-(MW2*Sf6587) + Sf6580*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2350 + Sf6569) - Sf6573) + 
     &         Sf6571*MSq2Diff(1,2)))

	Sf6458 = MHiggs2(4)*
     &    (-(MW2*Sf6457) + Sf6444*MHiggs2(4)) + 
     &   MW2**2*(5*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf2350 + Sf6424) - Sf6433) + 
     &         Sf6428*MSq2Diff(1,1)))

	Sf2160 = -3*MSl2Diff(2,1)**6 + 
     &   MW2*(Sf2151*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(Sf2159 + 
     &            MW2*(-Sf2154 + MW2*(2*Sf2152 - MW2*Sf294))) - 
     &         Sf2156*MSl2Diff(2,1)**2))

	Sf2091 = -3*MSl2Diff(1,1)**6 + 
     &   MW2*(Sf2072*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(Sf2090 + 
     &            MW2*(-Sf2079 + MW2*(2*Sf2075 - MW2*Sf294))) - 
     &         Sf2083*MSl2Diff(1,1)**2))

	Sf538 = 2*Sf521*MHiggs2(4)**2 + 
     &   MW2*(Sf537 + MW2*(Sf500*MSl2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(2*MW2*Sf175 + Sf502) - Sf507) - 
     &            Sf510*MSl2Diff(1,1))))

	Sf1443 = -2*Sf1432*MHiggs2(4)**2 + 
     &   MW2*(Sf1442 + MW2*
     &       (Sf1416*MSl2Diff(2,1)**3 + 
     &         MW2*(MW2*(-Sf1419 + MW2*(Sf1417 + 2*MW2*Sf175)) - 
     &            Sf1421*MSl2Diff(2,1))))

	Sf1324 = -3*MSl2Diff(2,1)**6 + 
     &   MW2*(Sf1315*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(Sf1323 + 
     &            MW2*(-Sf1320 + MW2*(2*Sf1316 - MW2*Sf364))) - 
     &         Sf1318*MSl2Diff(2,1)**2))

	Sf340 = -6*MSl2Diff(1,1)**6 + 
     &   MW2*(3*Sf321*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(6*MW2 - 2*Sf328) + Sf332) - 
     &            Sf339) - Sf323*MSl2Diff(1,1)**2))

	Sf1292 = -6*MSl2Diff(2,1)**6 + 
     &   MW2*(3*Sf1281*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(6*MW2 - 2*Sf1286) + Sf1288) - 
     &            Sf1291) - Sf1283*MSl2Diff(2,1)**2))

	Sf5439 = -8*MSq2Diff(2,2)**6 + 
     &   MW2*(4*Sf5408*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf5424) - Sf5433) + 
     &            Sf5438) - Sf5435*MSq2Diff(2,2)))

	end

