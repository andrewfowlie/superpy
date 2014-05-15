* TLhr_abbr1_6.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_6
	implicit none

#include "TLhr.h"


	Cd(9665) = MTy2*Cd(1739) + (Cd(6005) - 5*MStop2(2))*MStop2(2)

	Cd(9666) = 6*Cd(874)*MStop2(2) + MTy2*(8*MTy2 - 6*MStop2(3))

	Cd(9673) = MTy2*Cd(3250) + (Cd(9672) - 5*MStop2(2))*MStop2(2)

	Cd(9677) = -3*MStop2(2)**2 + MTy2*(9*MTy2 + 2*MStop2(2))

	Cd(9691) = MTy2*(32*MTy2 - 4*MStop2(1)) - 11*MStop2(1)**2

	Cd(9692) = MTy2*(2*MTy2 - MStop2(1)) + MStop2(1)**2

        Cd(9718) = 2*(Cd(9717)*MStop2(2) + 
     &     UCStopC(1,3)*(Cd(9716)*MStop2(3) - 
     &        Cd(369)*MStop2(4)*UCStopC(1,3)))

        Cd(9720) = 2*(Cd(9717)*MStop2(1) + 
     &     UCStopC(1,3)*(Cd(9719)*MStop2(3) - 
     &        Cd(369)*MStop2(4)*UCStopC(1,3)))

        Cd(963) = 6*MStop2(2)**2 + 
     &   2*MStop2(3)*(-2*MStop2(3) + MStop2(4))

        Cd(9724) = 2*(Cd(9723)*MStop2(2) + 
     &     UCStop(1,3)*(Cd(9722)*MStop2(3) - 
     &        Cd(368)*MStop2(4)*UCStop(1,3)))

        Cd(9726) = 2*(Cd(9723)*MStop2(1) + 
     &     UCStop(1,3)*(Cd(9725)*MStop2(3) - 
     &        Cd(368)*MStop2(4)*UCStop(1,3)))

	Cd(9728) = SA2*Cd(425) - S2A*Cd(426) + CA2*Cd(427)

        Cd(964) = -(Cd(963)*UStop2(1,2)**2) + 
     &   MStop2(1)*(6*MStop2(3) - 14*MStop2(1)*UStop2(1,2))

        Cd(9732) = 8*MTy*SA2 + 3*Cd(9730)*UCStop(1,3) + 
     &   3*Cd(9731)*UCStopC(1,3)

        Cd(9735) = 8*MTy*SA2 - 
     &   3*(Cd(9730)*UCStop(1,3) + Cd(9731)*UCStopC(1,3))

	Cd(9755) = 2*CA2*MUE2 + 2*At*AtC*SA2 - S2A*Cd(424)

	Cd(9776) = -Cd(538) + Cr(85)*UStop2(1,1)**3

        Cd(9781) = 2*MStop2(1)*
     &    (MStop2(1) - MStop2(2)*UStop2(1,1)) - 
     &   2*MStop2(2)**2*UStop2(1,2)

	Cd(9791) = -2*SA2*Cd(471) - 2*CA2*Cd(474) - S2A*Cd(567)

	Cd(9796) = CA2*MUE**2 + AtC*(-(MUE*S2A) + AtC*SA2)

        Cd(9801) = SA2*MStop2(3)*UCStopC(3,1) + 
     &   Cd(9792)**2*UUStopC(3,1) - Cd(9793)*Cd(9795)*UUStopC(3,3)

        Cd(9816) = SA2*MStop2(3)*UCStop(3,1) - 
     &   Cd(9807)**2*UUStopC(3,1) + Cd(9808)*Cd(9810)*UUStopC(3,4)

	Cd(9834) = 2 - 10*UStop2(1,1) + 6*UStop2(1,2)

        Cd(9853) = 2*(Cd(9852)*MStop2(1)**2 + 
     &     MStop2(2)*(4*Cd(779)*MStop2(1) + Cd(837)*MStop2(2)))

        Cd(9857) = Cd(9856)*MStop2(1)**2 + 6*MStop2(2)**2 + 
     &   4*MStop2(1)*MStop2(2)*UStop2(1,1)

        Cd(9858) = 2*(Cd(782)*MStop2(1)**2 + 
     &     MStop2(2)*(Cd(781)*MStop2(2) - 2*MStop2(1)*UStop2(1,2)))

        Cd(9861) = -6*Cd(290)*MStop2(1)**2 - 
     &   2*Cd(787)*MStop2(2)**2 + 4*MStop2(1)*MStop2(2)*UStop2(1,2)

        Cd(9865) = Cd(9864)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(9861)*UStop2(1,1) - Cd(853)*UStop2(1,2))

        Cd(9876) = 2*At*AtC*Cd(577) + 2*Cd(9875)*MStop2(3) - 
     &   Cd(425)*MStop2(4)

        Cd(9891) = Cd(424)*Cr(28) - 
     &   2*(AtC*MUE*UCStop(1,3)**2 + At*MUEC*UCStopC(1,3)**2)

        Cd(9896) = Cd(426)*Cd(9895) + 
     &   Cr(75)**2*(-2*MTy*Cd(13) + 
     &      Cd(424)*UStop2(1,1)*UStop2(1,2))

        Cd(9897) = MUE**2*Cd(9895)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(9895)*UCStopC(1,3)**2 + 
     &   2*MUE2*Cr(75)**2*UStop2(1,1)*UStop2(1,2)

        Cd(9902) = 2*MUE2*Cd(577)*Cr(75)**2 - 
     &   Cd(9901)*(MUE**2*UCStop(1,3)**2 + MUEC**2*UCStopC(1,3)**2)

        Cd(9905) = -(AtC**2*Cd(9901)*UCStop(1,3)**2) + 
     &   At*(2*AtC*Cd(9904) - At*Cd(9901)*UCStopC(1,3)**2)

        Cd(9906) = 4*(-2*MTy*Cd(585) + Cd(2036))*Cr(75)**2 + 
     &   (2*Cd(9905))/MStop2(3)

        Cd(9907) = Cd(424)*Cd(9904) - 
     &   Cd(9901)*(AtC*MUE*UCStop(1,3)**2 + 
     &      At*MUEC*UCStopC(1,3)**2)

        Cd(9919) = SA2*MStop2(3)*UStop(1,2) + 
     &   Cd(9792)*Cd(9808)*UStopC(2,1) - 
     &   Cd(9793)*Cd(9808)*UStopC(2,2)

        Cd(9920) = SA2*MStop2(3)*UStop(1,1) - 
     &   Cd(9792)*Cd(9807)*UStopC(2,1) + 
     &   Cd(9793)*Cd(9810)*UStopC(2,2)

        Cd(9932) = 2*Cd(9930) - Cd(9793)*Cd(9931)*MStop2(1) + 
     &   8*SA2*MStop2(3)**2*UStop2(1,1)

	Cd(9966) = Cd(1313) + UStop2(1,1)*UStop2(1,2)

        Cd(9998) = 2*At*UCStopC(1,3)**2 + 
     &   2*UStop2(1,1)*(M_3C + 2*AtC*UStop2(1,1)) + 
     &   AtC*UStop2(1,2)**2

        Cd(10052) = MUE*(M_3C*UCStop(1,3)**2 + 
     &      M_3*UStop2(1,1)*UStop2(1,2)) + 
     &   MUEC*(M_3*UCStopC(1,3)**2 + 
     &      M_3C*UStop2(1,1)*UStop2(1,2))

	Cd(10177) = Cd(12016) + At*(-Cd(9948) + UStop2(1,2))

        Cd(10193) = Cd(1313) + 
     &   3*UStop2(1,1)*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

        Cd(10204) = (-(Cd(152)*Cd(2597)) + 
     &     8*MGlpTmSt2(1)*MStop2(1)**2 + 
     &     Cd(164)*(4*Cd(9948)*MStop2(1)**2 + 
     &        4*MStop2(1)**2*UStop2(1,1)))/MStop2(1)**2

        Cd(10207) = Cd(10353) + M_3C*At*Cd(9948)*UStop2(1,1) + 
     &   M_3*AtC*Cd(9948)*UStop2(1,1)

        Cd(10208) = -Cd(599) + 
     &   3*UStop2(1,1)*(UStop2(1,1) + UStop2(1,2)**2)

        Cd(10295) = (-Cd(471) + 2*Cd(9948)*MStop2(3) + 
     &     MStop2(3)*(UStop2(1,1) - UStop2(1,2)))/MStop2(3)

	Cd(10297) = Cd(471) + MStop2(3)

	Cd(10300) = Cd(10294) + Cd(10299)*MStop2(3)

        Cd(10342) = 2*Cd(6277)*UStop2(1,1) + 
     &   UStop2(1,1)**2*(-Cd(10340) + 4*UStop2(1,1)) + 
     &   Cd(10341)*UStop2(1,2)

        Cd(10344) = -((Cd(172) - UStop2(1,1))*
     &     (UStop2(1,1)**2 + UStop2(1,2)**2))

        Cd(10355) = -(Cd(10164)*UStop2(1,1)**2) + 
     &   4*UStop2(1,1)**3 - 
     &   (Cd(2361) - 12*UStop2(1,1))*UStop2(1,2)**2

        Cd(10358) = -(AtC*Cd(10356)) + Cd(10355)*MStop2(3) + 
     &   At*(Cd(10357) - At*Cr(85)*UCStopC(1,3)**2)

        Cd(10362) = AtC*Cd(10359) - 
     &   At*(Cd(10361) + At*UCStopC(1,3)**2) + 
     &   10*MStop2(3)*UStop2(1,1)

        Cd(10414) = M_3*AtC*UCStop(1,3) + 
     &   MStop2(1)*UCStop(1,3) + 
     &   At*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))

        Cd(10471) = -(Cd(2429)*UStop2(1,1)) - 
     &   Cd(2364)*UStop2(1,2) + UStop2(1,1)**2*(1 + 2*UStop2(1,2))

        Cd(10475) = Cd(425) + 4*Cd(152)*Cr(28) + 
     &   2*At*AtC*UStop2(1,1)*UStop2(1,2)**2

	Cd(10476) = 4*Cd(152) + Cd(10191)

        Cd(10502) = 3*MStop2(3) - 3*At**2*UCStopC(1,3)**2 + 
     &   4*MTy2*UStop2(1,1)*(-Cd(9948) + UStop2(1,2))

        Cd(10507) = -2*MTy2*Cd(10506) + 
     &   (2*MGlpT2*MTy2 + Cd(152)*MGlpTmSt2(1))*UCStop(1,3)

	Cd(10563) = Cd(11574) + (-2 - 2*UStop2(1,1))*UStop2(1,2)

        Cd(10577) = Cd(47) + Cd(10383)*MStop2(2) + 
     &   MGlmT2*(Cd(9948)*MStop2(1) - Cd(1272)*UStop2(1,2))

	Cd(10579) = Cd(10448) + 2*MUE2*Cd(9948)*UStop2(1,1)**2

        Cd(10584) = (M_3*MUE**2 + 2*M_3C*MUE2)*UCStop(1,3) + 
     &   (2*M_3*MUE2 + M_3C*MUEC**2)*UCStopC(1,3)

        Cd(10616) = 2*M_3C*UCStop(1,3)**2 + 
     &   Cd(10614)*UStop2(1,1) + 
     &   At*UStop2(1,2)*(-UStop2(1,1)**2 + Cd(10615)*UStop2(1,2))

	Cd(10623) = Cd(10608) + MStop2(3)

	Cd(10630) = Cd(12517) + 2*Cd(585)*Cd(9955)*Cr(75)*MStop2(1)

        Cd(10649) = -(Cd(10647)*UStop2(1,1)**2) + 
     &   UStop2(1,1)**3 + (-Cd(10648) + UStop2(1,1))*UStop2(1,2)**2

        Cd(10658) = Cd(1313) + 
     &   UStop2(1,1)*(UStop2(1,2) + UStop2(1,1)*UStop2(1,2))

        Cd(10692) = MGlpTmSt2(1) + 6*M_3*MTy*UCStopC(1,3) + 
     &   (MTy2 - MStop2(1))*UStop2(1,2)

        Cd(10786) = Cd(425) + 4*Cd(179)*Cr(28) + 
     &   2*At*AtC*UStop2(1,1)*UStop2(1,2)**2

	Cd(10788) = 4*Cd(179) + Cd(10191)

        Cd(10886) = MUE**2*(MTy2 - MStop2(1))*UCStop(1,3)**2 + 
     &   MUEC**2*MStop2(3)*UCStopC(1,3)**2 + 
     &   MTy2*MUE2*UStop2(1,1)**2

	Cd(10902) = 4*MTy2 - Cd(10191)

        Cd(11224) = Cd(11316) + 
     &   2*(UStop2(1,1) - 3*UStop2(1,2))*UStop2(1,2)

        Cd(11248) = Cd(11054)*(MStop2(2) + MStop2(3)) + 
     &   MStop2(1)*(At*Cd(11247) - 4*MTy*UStop(1,1)*UStop2(1,2))

        Cd(11387) = At*(UCStop(3,1)*UStop2(1,2)*UUStop(2,2) - 
     &      Cd(970)*UUStop(3,2)) + 
     &   MTy*(-3*Cd(2364)*UUStop(3,3) - Cd(970)*UUStop(3,4))

	Cd(11450) = Cd(11307) + 2*Cd(11449)

	Cd(11497) = Cd(11351) + 2*Cd(11496)

	Cd(11519) = Cd(11366) + 2*Cd(11518)

        Cd(11528) = At*(UCStop(3,1)*UStop2(1,2)*UUStop(2,2) - 
     &      Cd(2364)*UUStop(3,2)) - 
     &   MTy*(Cd(2429)*UUStop(3,3) + Cd(2364)*UUStop(3,4))

        Cd(1152) = 4*Cr(54)*MStop2(1)**2 - 
     &   4*Cd(1150)*MStop2(1)*MStop2(2) + 12*MStop2(2)**2 + 
     &   Cd(1151)*MStop2(3)

        Cd(11624) = UStop2(1,2)**2*(Cd(11623) + UStop2(1,2)) + 
     &   UStop2(1,1)*(2*UStop2(1,1)**2 + 2*Cd(11622)*UStop2(1,2))

        Cd(11632) = -8*MTy2*MStop2(3)*UStop2(1,2) + 
     &   Cd(11630)*UStop2(1,2)**2 + 
     &   MStop2(2)*(2*At*AtC*Cd(9948)*UStop2(1,1)**2 + 
     &      Cd(11631)*UStop2(1,2)**3)

        Cd(11715) = 3*M_3C*Cd(47)*Cd(118) + 
     &   4*M_3C*MGl2*MTy2*Cd(152) + 
     &   AtC*(Cd(47)*Cd(118) + MGl2*Cd(152)**2)

        Cd(11717) = 3*M_3*Cd(47)*Cd(118) + 
     &   4*M_3*MGl2*MTy2*Cd(152) + 
     &   At*(Cd(47)*Cd(118) + MGl2*Cd(152)**2)

        Cd(11728) = At*(2*Cd(11727) + At*UCStopC(1,3)**2) + 
     &   AtC*(Cd(9967)*UCStop(1,3)**2 + 
     &      2*M_3*UStop2(1,1)*UStop2(1,2))

	Cd(11730) = Cd(471) + 2*MStop2(3)

        Cd(11735) = Cd(10353) + 
     &   M_3C*At*UStop2(1,1)*UStop2(1,2) + 
     &   M_3*AtC*UStop2(1,1)*UStop2(1,2)

        Cd(11736) = -4*Cd(152)**2*Cd(11735) + 
     &   (MGlmT2*Cd(47)*(-Cd(10294) + Cr(28)*MStop2(3) + 
     &        2*Cd(9948)*MStop2(3)*UStop2(1,1)))/MStop2(3)

        Cd(11742) = (Cr(85) + UStop2(1,1))*
     &   (UStop2(1,1)**2 + UStop2(1,2)**2)

	Cd(11752) = Cd(1313) - Cd(9948)*Cd(11502) - Cd(2429)*UStop2(1,1)

        Cd(11755) = 4*Cd(152) - AtC*Cd(11753) + 
     &   At*(2*Cd(11754) + At*UCStopC(1,3)**2)

	Cd(11783) = Cd(766) + At*MUE*UStop2(1,1)*UStop2(1,2)

        Cd(11811) = -(MGlmT2*UStop2(1,1)) + MGlpT2*UStop2(1,1) + 
     &   M_3C*At*UStop2(1,2)

	Cd(11817) = Cd(10510) + M_3*AtC*UCStop(1,3)

        Cd(11859) = -MGlmT2 - 2*MStop2(2) + 
     &   Cd(138)*(5*Cd(9948) + Cd(11858)*UStop2(1,1))

        Cd(11874) = 4*Cd(179) - AtC*Cd(11753) + 
     &   At*(2*Cd(11754) + At*UCStopC(1,3)**2)

        Cd(11878) = Cd(11876)*UStop2(1,1)**2 + 
     &   Cd(11877)*UStop2(1,2)**2 + 
     &   UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Cd(11908) = Cd(11907)*UStop2(1,1)**2 + 
     &   Cr(85)*UStop2(1,2)**2 + 
     &   UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Cd(11946) = MStop2(2) - 
     &   2*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3)) - 
     &   MGl2*UStop2(1,2)

        Cd(11965) = AtC*(M_3*Cr(28) - 
     &      2*M_3C*UCStop(1,3)**2) + 
     &   At*(M_3C*Cr(28) - 2*M_3*UCStopC(1,3)**2)

        Cd(11990) = 4*Cd(152) + Cd(425) + 
     &   2*At*AtC*UStop2(1,1)*UStop2(1,2)**2

        Cd(11991) = 4*Cd(179) + Cd(425) + 
     &   2*At*AtC*UStop2(1,1)*UStop2(1,2)**2

        Cd(12014) = UStop2(1,1)*
     &    (Cd(12013) + 2*UStop2(1,1)*(-7 + 4*UStop2(1,1))) - 
     &   2*Cd(2368)*UStop2(1,2)**2

        Cd(12018) = 4*Cd(152) - At*Cd(12017) + 
     &   AtC*(-Cd(12015) + AtC*UCStop(1,3)**2)

        Cd(12020) = UStop2(1,1)*
     &    (Cd(12013) - 2*(Cd(12019) - 4*UStop2(1,1))*UStop2(1,1))-
     &     2*Cd(2368)*UStop2(1,2)**2

        Cd(12021) = 4*Cd(179) - At*Cd(12017) + 
     &   AtC*(-Cd(12015) + AtC*UCStop(1,3)**2)

        Cd(12025) = -(Cd(12024)*UStop2(1,1)**2) + 
     &   4*UStop2(1,1)**3 - 
     &   2*(Cd(290) - 2*UStop2(1,1))*UStop2(1,2)**2

	Cd(12030) = Cd(12035) - 2*UStop2(1,1)

        Cd(12168) = Cd(138)**2 - 
     &   MTy2*(MGl2 + MStop2(2) - MGlpTmSt2(2)*UStop2(1,1))

        Cd(12175) = 4*MTy2*(MTy2 + MStop2(2)) - 
     &   (5*MTy2 + 2*MStop2(1))*MStop2(3)

        Cd(12202) = 2*MGlmT2**2 - 5*MGlpT2*MStop2(1) + 
     &   3*MStop2(1)**2 + MGlpTmSt2(1)*MStop2(2)

	Cd(12210) = Cd(179) + MTy2*(4*UStop2(1,1) - 2*UStop2(1,2))

        Cd(12231) = MGlmT2**2 - 3*MGlpT2*MStop2(1) + 
     &   2*MStop2(1)**2 + MGlpTmSt2(1)*MStop2(2)

	Cd(12280) = 4*MGl2 + 4*MTy2 + Cd(12279)*MStop2(1) - 4*MStop2(3)

        Cd(12314) = 2*MTy2 + 2*MGl2*Cd(10112) + 
     &   MTy*(Cd(10398)*UCStop(1,3) + Cd(10140)*UCStopC(1,3))

        Cd(12350) = 2*MTy2 + 2*MGl2*Cd(10112) - 
     &   MTy*(Cd(10398)*UCStop(1,3) + Cd(10140)*UCStopC(1,3))

        Cd(1238) = -(Cr(85)*UStop2(1,2)**2) + 
     &   UStop2(1,1)*((4 - 6*UStop2(1,1))*UStop2(1,1) - 
     &      Cd(1237)*UStop2(1,2))

        Cd(12585) = 4*MTy2**2*Cd(10193) + 
     &   MGl2*(AtC**2*UCStop(1,3)**2 + 
     &      At*(At*UCStopC(1,3)**2 + 2*AtC*Cd(9981)*UStop2(1,1)))

        Cd(12669) = Cd(11574) - 
     &   2*(2*Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

	Cd(124) = MTy2**2 + Cd(12266) - 2*MGlpT2*MStop2(1)

        Cd(1267) = UStop2(1,1)*
     &    (2*Cd(1263)*UStop2(1,1) + Cd(1266)*UStop2(1,2)) + 
     &   4*MStop2(2)*(Cd(1213)*MStop2(2) + 
     &      3*MStop2(3)*UStop2(1,2)**2)

        Cd(12750) = 2*Cd(12700)*Cd(12701)*Cd(12749) + 
     &   2*MStop2(2)*MStop2(3) + 
     &   MStop2(4)*(Cd(12700)**2*UCStop(1,3)**2 + 
     &      Cd(12701)**2*UCStopC(1,3)**2)

        Cd(1275) = Cd(1274)*MStop2(3) + Cd(778)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(-4*MStop2(2)*UStop2(1,1) + 
     &      MStop2(4)*UStop2(1,2))

        Cd(12858) = 4*Cd(12739) + 
     &   2*Cd(1272)*Cd(12700)**2*UCStop(1,3)**2 - 
     &   2*Cd(12701)**2*MStop2(3)*UCStopC(1,3)**2 + 
     &   4*Cd(12736)*MStop2(2)*UStop2(1,1)

        Cd(13147) = Cd(1213) + 
     &   UStop2(1,1)*(-Cd(13146) + UStop2(1,1)*UStop2(1,2)**2)

        Cd(13260) = MStop2(3) + Cd(13256)*UUStopC(1,1) + 
     &   Cd(13258)*UUStopC(1,2) + 2*Cd(13259)*UUStopC(1,3)

	Cd(13261) = 2*MGl2 + Cd(1701)

        Cd(13378) = 1/(512.D0*Pi**4)*
     &   (-(Cd(12700)*Cd(13377)) + 
     &      MStop2(2)*(Cd(12701)*
     &          (Cd(12701)*Cd(13372) + 512*Pi**4*Cd(13375))*
     &          MStop2(2) + 4*MStop2(3)))/MStop2(2)**2

	Cd(13474) = Cd(13479) + Cd(13805)

        Cd(13480) = MTy*Cd(13473) + Cd(14305) + 
     &   AtC*Cd(13469)*UStopC(1,2)

        Cd(13513) = AtC*Cd(13506) + MUEC*TB*Cd(13507) + 
     &   At*Cd(13508) + MTy*MUE*TB*UCStop(1,3)

        Cd(13532) = MTy*(Cd(13526) - MUE*TB*UCStop(1,3)) + 
     &   Cd(13479)*UStop2(1,1) + AtC*Cd(13504)*UStopC(2,2)

        Cd(13542) = -(Cd(13541)*UCStopC(3,1)) + 
     &   MTy*(Cd(13540) + MUE*TB*UCStopC(3,4)) + 
     &   AtC*Cd(13504)*UStopC(1,2)

        Cd(13544) = -(Cd(13541)*UCStop(3,1)) + 
     &   MTy*(Cd(13543) + MUE*TB*UCStop(3,3)) + 
     &   AtC*Cd(13469)*UStopC(2,2)

        Cd(13551) = -(Cd(13550)*UCStopC(3,1)) + 
     &   MTy*(Cd(13540) + MUE*TB*UCStopC(3,4)) + 
     &   AtC*Cd(13504)*UStopC(1,2)

        Cd(13552) = -(Cd(13550)*UCStop(3,1)) + 
     &   MTy*(Cd(13543) + MUE*TB*UCStop(3,3)) + 
     &   AtC*Cd(13469)*UStopC(2,2)

        Cd(13591) = MTy*(Cd(13584) + 3*MUE*TB*UCStop(1,3)) + 
     &   2*Cd(13541)*UStop2(1,2) + AtC*Cd(13590)*UStopC(1,2)

        Cd(13593) = MTy*(Cd(13586) - 3*MUE*TB*UCStop(1,3)) + 
     &   2*Cd(13550)*UStop2(1,1) + AtC*Cd(13592)*UStopC(2,2)

        Cd(1356) = Cd(1355) + 4*Cr(28)*MStop2(3)**2 - 
     &   8*Cd(852)*UStop2(1,1)*UStop2(1,2) + 
     &   12*MStop2(1)*MStop2(2)*UStop2(1,2)**2

        Cd(13629) = At*Cd(13627) + MUEC*TB*Cd(13627) + 
     &   AtC*Cd(13628) + MUE*TB*Cd(13628)

        Cd(13697) = 2*MGl2*MGlmT2*MSbot2(1) - 
     &   MStop2(2)*(Cd(13695) + 
     &      MStop2(2)*(-2*Cd(13696) + MStop2(2)))

	Cd(13704) = Cd(13702) + 4*MSbot2(1)

	Cd(13706) = At*Cd(12700) + MUEC*TB*Cd(12700) - Cd(13642)

	Cd(13722) = At*AtC + AtC*MUEC*TB + Cd(13479)

	Cd(13743) = 1 + Cd(13737) + UStop2(1,1)**2

	Cd(13746) = Cd(1179) - UStop2(1,1)*(1 + 2*UStop2(1,1))

	Cd(13748) = Cd(2366) + (-2 - 2*UStop2(1,1))*UStop2(1,1)

	Cd(13782) = Cd(13736) + UStop2(1,1)**2

	Cd(13786) = Cd(6268) - 2*UStop2(1,1)*UStop2(1,2)

	Cd(13788) = Cd(9948) + UStop2(1,1) - UStop2(1,2)

	Cd(13827) = MGl2*(MGl2 - 3*MTy2) - MTy2**2 + Cd(153)*MSbot2(1)

        Cd(13829) = MGl2*MGlmT2*MSbot2(1) - 
     &   MStop2(2)*(Cd(13827) + MStop2(2)*(-Cd(13828) + MStop2(2)))

	Cd(13837) = MSq2Diff(1,2)*(-Cd(13728) + Cr(56)*MSq2Diff(1,2))

        Cd(1378) = 2*Cd(1377) - 12*MTy2*MStop2(1)*UStop2(1,2) + 
     &   MStop2(2)*(3*Cd(579) + 8*MTy2*UStop2(1,1)*UStop2(1,2)**2)

        Cd(13844) = 2*Cd(13744)*MSbot2(1)*MStop2(2) - 
     &   Cd(13843)*MStop2(2)**2 + 
     &   (-(Cr(56)*MSbot2(1)**2) + Cd(13728)*MSq2Diff(1,2))*
     &    UStop2(1,1)

        Cd(13949) = -(Cd(13948)*MSq2Diff(1,2)) - 
     &   MStop2(2)*(2*Cr(75)*MSbot2(1) + Cd(13947)*MStop2(2)) + 
     &   Cr(56)*MSbot2(1)**2*UStop2(1,1)

        Cd(13952) = Cd(13951)*MSbot2(1)*MStop2(2) - 
     &   Cd(13744)*MStop2(2)**2 + 
     &   (-(Cr(56)*MSbot2(1)**2) + Cd(13942)*MSq2Diff(1,2))*
     &    UStop2(1,1)

	Cd(13965) = Cd(13963) + 4*MSbot2(1)

        Cd(13977) = At*Cd(12700)*Cr(56)*UStop2(1,1) + 
     &   AtC*MUEC*TB*UStop2(1,1)**2 + MUE2*TB2*UStop2(1,1)**2 + 
     &   2*MSq2Diff(1,2)*UStop2(1,2)

	Cd(13979) = Cd(9948) + 2*UStop2(1,1) - 2*UStop2(1,2)

	Cd(13986) = Cd(2380) + 2*UStop2(1,1)*(-1 + 2*UStop2(1,1))

	Cd(13989) = Cd(1179) + UStop2(1,1)*(1 + 2*UStop2(1,1))

        Cd(13995) = MStop2(2)**2 + 
     &   MGl2*(MGl2 - 4*MSbot2(1) + 2*MStop2(2))

        Cd(14071) = MTy2**2 + 2*At*MTy2*MUE*TB + 
     &   MGl2*Cd(14069) - Cd(14070)*MStop2(1)

        Cd(14109) = MUEC*TB*Cd(14037)*UStop(2,2) + 
     &   MTy*MUE*TB*UUStop(3,1) + At*MTy*UUStop(3,2) + 
     &   Cd(13698)*UUStop(3,3)

	Cd(14127) = MTy2*TB + Cd(14126)

	Cd(14275) = MGl2*Cd(13638) - Cd(14280)

	Cd(14284) = Cd(99) + MSbot2(1)

        Cd(14310) = Cd(13526)*(Cd(14308) - Cd(14309)) + 
     &   MSq2Diff(1,2)*(Cd(12700)*Cd(13561) + MSq2Diff(1,2))

        Cd(14335) = Cd(13526)*Cd(14308) + 
     &   6*Cd(14328)*MSq2Diff(1,2) + 
     &   2*Cd(13504)*MSbot2(1)*
     &    (-(MTy*UStop2(1,2)*UStopC(2,1)) + 
     &      Cd(12700)*UStop2(1,1)*UStopC(2,2))

        Cd(1454) = Cd(1265) - Cd(1328) + Cd(1416)*MStop2(3) - 
     &   Cd(1453)*UStop2(1,1)**2

	Cd(1460) = Cd(1379) + 2*Cd(1459)

	Cd(144) = MGl2**2 + Cd(4591) - 2*MGlpT2*MStop2(2)

        Cd(1504) = -2*Cd(1502) - MGl2*Cd(1503) + 
     &   2*(MGl2**2 + Cd(45)*MStop2(2))

        Cd(1536) = MGl2*(-4*MGl2 + Cd(1530)) - MTy2*Cd(1535) - 
     &   3*Cd(637)*MStop2(2)

        Cd(1547) = MTy2*Cd(1544) + 
     &   MStop2(2)*(-(MTy2*Cd(1545)) + 
     &      MStop2(2)*(-Cd(1546) + 15*MStop2(2)))

	Cd(154) = Cd(13759) - MTy*Cd(152)*Cr(5)

        Cd(1586) = -13*MStop2(1)**2 + 
     &   MTy2*(8*MTy2 + 11*MStop2(1)) - 3*Cd(45)*MStop2(2)

        Cd(1587) = MGl2*(-4*MGl2 + Cd(1582)) - MTy2*Cd(1584) - 
     &   3*Cd(1583)*MStop2(1)

        Cd(1593) = 9*MStop2(1)**2 + 
     &   MTy2*(12*MTy2 + 7*MStop2(1)) + Cd(1592)*MStop2(2)

        Cd(1597) = MTy2*(MTy2*(24*MTy2 - 3*Cd(1583)) - 
     &      2*Cd(1596)) + 3*Cd(1594)*MStop2(1)**2

        Cd(1621) = -8*MStop2(1)**2 + 
     &   MTy2*(3*MTy2 + 7*MStop2(1)) - Cd(45)*MStop2(2)

        Cd(1645) = -3*MStop2(1)**2 + 
     &   MTy2*(2*MTy2 + 3*MStop2(1)) - Cd(45)*MStop2(2)

        Cd(1712) = 2*(MTy2*Cd(45)*Cd(1592) + 
     &     MStop2(2)*(-3*Cd(1688) + 
     &        2*(3*Cd(45) - 7*MStop2(2))*MStop2(2)))

        Cd(1716) = MTy2*Cd(1713) + 
     &   MStop2(2)*(-Cd(1714) + 
     &      MStop2(2)*(-Cd(1715) + 15*MStop2(2)))

        Cd(1733) = MTy2*(3*MTy2 - 4*MStop2(1)) - 
     &   (-6*MTy2 + Cd(1732))*MStop2(2)

        Cd(1736) = 3*Cd(1735)*MStop2(2)**2 + MStop2(2)**3 + 
     &   MTy2**2*(Cd(1735) + 3*MStop2(2))

        Cd(1740) = MTy2*Cd(1737) + 
     &   MStop2(2)*(-(MTy2*Cd(1738)) + 
     &      MStop2(2)*(-Cd(1739) + 5*MStop2(2)))

        Cd(1741) = MGl2*(MGl2*(MGl2 - Cd(1731)) + 2*Cd(1734)) + 
     &   2*Cd(1736)

        Cd(1766) = -(MGl2*
     &      (MGl2*(-3*MGl2 + Cd(1763)) + Cd(1765))) + 
     &   Cd(1760)*MStop2(2)

        Cd(1769) = -(Cd(1122)*MStop2(2)**2) + 
     &   MTy2*(MTy2*(-8*MTy2 + Cd(637)) + Cd(1768)*MStop2(2))

        Cd(1808) = Cd(1805)*Cd(1807) - 
     &   MGl2*(MGl2*(MGl2 + Cd(1806)) + 3*Cr(26)*MStop2(2))

        Cd(1812) = Cd(1799)*Cd(1807) + 
     &   MGl2*(MGl2*(3*MGl2 - Cd(1809)) - Cd(1811))

        Cd(1815) = -(Cd(1594)*MStop2(1)**2) + 
     &   MTy2*(MTy2*(-8*MTy2 + Cd(1583)) + Cd(1814)*MStop2(1))

	Cd(1820) = Cd(1502) - Cd(45)*MStop2(2)

	Cd(180) = Cd(6391) + MTy*Cd(179)*Cr(5)

        Cd(1846) = 4*MTy2**2 + 
     &   2*MStop2(2)*(Cd(67) + 9*MStop2(2)) + MStop2(3)*MStop2(4)

	Cd(1847) = Cd(949) + 2*MTy2*(MTy2 - 2*MStop2(2))

	Cd(1865) = 4*MTy2**2 + Cd(982) + 2*Cd(1864)*MStop2(2)

        Cd(1944) = -3*MStop2(1)**2 + MTy2*(6*MTy2 + MStop2(1)) + 
     &   3*Cd(1592)*MStop2(2)

        Cd(1949) = MStop2(1)**3 + 
     &   MTy2*(-15*MStop2(1)**2 + MTy2*(15*MTy2 + 15*MStop2(1)))

	Cd(1963) = Cd(1894) - 4*MTy2*MStop2(2)

        Cd(1964) = 3*Cd(44)**3 + 
     &   MGl2*(3*MGl2*(MGl2 - Cd(1592)) - Cd(1963))

        Cd(1976) = MTy2*(MTy2*(22*MTy2 + Cd(1975)) - 
     &      2*Cd(1974)*MStop2(1)) + 3*MStop2(1)**2*MStop2(3)

	Cd(1977) = 6*MGl2 + 2*MTy2 - 26*MStop2(2) + 23*MStop2(3)

        Cd(1984) = 10*MTy2*(MTy2 - MStop2(2)) + 
     &   (7*MTy2 + MStop2(1))*MStop2(3)

	Cd(1986) = 3*MGl2 + Cd(2812)

        Cd(1997) = 12*MTy2*(MTy2 - MStop2(2)) + 
     &   (5*MTy2 - MStop2(2))*MStop2(3)

	Cd(2003) = Cd(2957) + 4*(MGl2 - 5*MStop2(2))

        Cd(2007) = 10*MTy2*(MTy2 - MStop2(2)) + 
     &   (3*MTy2 - MStop2(2))*MStop2(3)

	Cd(2009) = 3*MGl2 + Cd(319)

	Cd(2019) = Cd(1922) + Cd(1496)*MStop2(1)

        Cd(2038) = -9*MStop2(2)**2 + 
     &   MTy2*(3*MTy2 - 4*MStop2(1) + 10*MStop2(2))

        Cd(2039) = 3*Cd(48)**3 + 
     &   MGl2*(3*MGl2*(MGl2 - Cd(1559)) - Cd(2038))

        Cd(2040) = MTy2**2 + MGl2*(7*MGl2 + 3*MTy2) - 
     &   3*Cd(1496)*MStop2(1)

	Cd(2048) = Cd(523) + 2*MTy2*(MTy2 - MStop2(4))

        Cd(2051) = 5*MStop2(1)**3 + 
     &   MTy2*(-7*MStop2(1)**2 + MTy2*(6*MTy2 + 2*MStop2(1)))

	Cd(2053) = Cd(1887) - 3*MGlpTmSt2(1)*MStop2(1)

        Cd(2082) = MTy2*Cd(2079) + 
     &   MStop2(2)*(Cd(2080) + MStop2(2)*(-Cd(2081) + 5*MStop2(2)))

        Cd(2083) = 2*MTy2*Cd(2079) + 
     &   2*MStop2(2)*(Cd(2080) + 
     &      MStop2(2)*(-Cd(2081) + 5*MStop2(2)))

        Cd(2115) = 4*MStop2(2)**2 - 2*Cd(2114)*MStop2(3) + 
     &   4*MTy2*(MTy2 - MStop2(4))

        Cd(2139) = -5*MStop2(1)**2 + 
     &   MTy2*(2*MTy2 + 5*MStop2(1)) - Cd(45)*MStop2(2)

        Cd(2152) = MTy2*Cd(2149) + 
     &   MStop2(2)*(Cd(2150) + 
     &      (3*Cd(2151) - 5*MStop2(2))*MStop2(2))

        Cd(2156) = MTy2*Cd(2153) - 
     &   MStop2(2)*(-Cd(2154) + Cd(2155)*MStop2(2) + MStop2(2)**2)

        Cd(2223) = -17*MStop2(1)**2 + 
     &   MTy2*(4*MTy2 + 9*MStop2(1)) - Cd(2221)*MStop2(2)

        Cd(2225) = 8*MTy2*(MTy2 - MStop2(2)) + 
     &   (5*MTy2 + MStop2(1))*MStop2(3)

        Cd(2226) = MTy2*(-2*Cd(1421)*MStop2(1) + 
     &      MTy2*(16*MTy2 - 5*MStop2(3))) + 
     &   3*MStop2(1)**2*MStop2(3)

        Cd(2249) = MTy2*Cd(2246) + 
     &   MStop2(2)*(Cd(2247) + (Cd(2248) - 5*MStop2(2))*MStop2(2))

        Cd(2253) = MTy2**2*Cd(2251) + 
     &   MStop2(2)*(-3*MTy2*Cd(155) + 
     &      MStop2(2)*(-3*Cd(2252) + 5*MStop2(2)))

        Cd(2270) = -(MGl2*
     &      (MGl2*(2*MGl2 - Cd(1540)) + Cd(2250))) + Cd(2253)

        Cd(2281) = -3*MStop2(2)**2*MStop2(3) + 
     &   MTy2*(-2*Cd(1122)*MStop2(2) + 
     &      MTy2*(16*MTy2 + 5*MStop2(3)))

        Cd(2305) = MGl2*(Cd(44)**2*Cd(508) + 
     &      MGl2*(-(Cd(44)*Cd(519)) - MGl2*Cr(85))) + 
     &   Cd(44)**3*UStop2(1,1)

        Cd(2343) = 20*MStop2(1)**2 + Cd(1583)*MStop2(3) + 
     &   MTy2*(4*MTy2 + 4*MStop2(4))

        Cd(2348) = Cd(1807)*Cd(2345) + 
     &   MTy2*(-Cd(2347) + MTy2*(6*MTy2 - 2*Cr(1)))

        Cd(2353) = MTy2*Cd(2352) + 
     &   (MTy2**2*(4*MTy2 + 4*Cr(26)) - Cd(2350)*MStop2(1))*
     &    UStop2(1,1)

        Cd(2456) = MGl2*(Cd(48)**2*Cd(508) + 
     &      MGl2*(-(Cd(48)*Cd(519)) - MGl2*Cr(85))) + 
     &   Cd(48)**3*UStop2(1,1)

        Cd(2522) = 2*MTy2*Cd(2149) + 
     &   2*MStop2(2)*(Cd(2150) + 
     &      (3*Cd(2151) - 5*MStop2(2))*MStop2(2))

        Cd(2545) = 20*MStop2(2)**2 - Cd(637)*MStop2(3) + 
     &   MTy2*(4*MTy2 + 4*MStop2(4))

        Cd(2549) = -(MTy2*
     &      (Cd(2547) + MTy2*(-6*MTy2 + Cr(27)))) + 
     &   Cd(2548)*MStop2(2)

        Cd(2552) = MTy2*Cd(2551) + 
     &   (MTy2**2*(4*MTy2 + 4*Cr(1)) - Cd(2550)*MStop2(2))*
     &    UStop2(1,1)

        Cd(2565) = 2*MGl2*Cd(2564) + 
     &   2*M_3*Cd(2377)*MStop2(2)*MStop2(3) + 
     &   M_3C*UCStop(1,3)**2*(Cd(2562) - 48*MGl2**2*UStop2(1,1))

        Cd(2822) = 6*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 8*MStop2(1)) - 9*MStop2(1)**2)

        Cd(2895) = -2*MStop2(1)**2*
     &    (-2*MTy2**2 + Cd(2894)*MStop2(1)) - 
     &   MTy2*(MTy2**2 - 3*MStop2(1)**2)*MStop2(4)

        Cd(2913) = 18*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 4*MStop2(1)) - 9*MStop2(1)**2)

        Cd(2923) = 17*MStop2(1)**2 + 
     &   MTy2*(6*MTy2 + 8*MStop2(1)) - Cd(2922)*MStop2(2)

        Cd(2925) = -14*MStop2(1)**2 + 
     &   MTy2*(4*MTy2 + 7*MStop2(1)) + Cd(155)*MStop2(2)

        Cd(2928) = 30*MStop2(1)**2 + 
     &   MTy2*(12*MTy2 + 23*MStop2(1)) + Cd(1735)*MStop2(2)

        Cd(289) = Cd(287)*UStop2(1,1) + 
     &   MStop2(1)*(MStop2(1) + 2*MStop2(2)*UStop2(1,1)) - 
     &   MStop2(2)**2*UStop2(1,2)

        Cd(2959) = 6*MStop2(1)**2 + 
     &   MTy2*(12*MTy2 + 7*MStop2(1)) + Cd(2958)*MStop2(2)

        Cd(2960) = MGl2*(Cd(2956) - MGl2*Cd(2957)) - 
     &   MTy2*Cd(2959) - 9*MStop2(1)**2*MStop2(3)

	Cd(2974) = Cd(3524) + Cd(2973)*MStop2(2)

	Cd(2981) = MGl2*(-MGl2 + Cd(2200)) - Cd(7559)

        Cd(2996) = 19*MStop2(1)**2 + 
     &   MTy2*(18*MTy2 + 13*MStop2(1)) + Cd(2995)*MStop2(2)

        Cd(2997) = MGl2*(-(MGl2*Cd(2875)) + Cd(2993)) - 
     &   MTy2*Cd(2996) + 6*MStop2(1)**2*MStop2(3)

        Cd(3017) = (-18*MGl2**3*Cd(607) + 
     &     MStop2(1)*(MGl2**2*Cd(3014) + 
     &        MStop2(1)*(8*MGl2*Cd(3015) - 9*Cd(3016)*MStop2(1))))*
     &   MStop2(3)

        Cd(3027) = MStop2(3)*
     &   (MTy2*(MTy2*(12*MTy2 + Cd(1583)) + 
     &        2*Cd(3026)*MStop2(1)) + 9*MStop2(1)**2*MStop2(3))

        Cd(3050) = -4*MStop2(1)**3 + 
     &   MTy2*(31*MStop2(1)**2 + MTy2*(3*MTy2 + 2*MStop2(1)))

        Cd(3063) = MGl2*(2*MGl2*(MGl2 - 3*Cd(45)) + Cd(2902)) - 
     &   Cd(44)*Cd(3062)

        Cd(3071) = 7*MStop2(1)**2 + 
     &   MTy2*(6*MTy2 + 4*MStop2(1)) + Cd(2922)*MStop2(2)

        Cd(3119) = -(Cd(48)*Cd(3117)) + 
     &   MGl2*(MGl2*(MGl2 - 3*Cd(49)) + Cd(3118))

        Cd(3131) = MTy2**2*Cd(2973) + 
     &   MStop2(2)*(-(MTy2*Cd(3129)) + 
     &      MStop2(2)*(-Cd(3130) + 13*MStop2(2)))

        Cd(3133) = MGl2*(15*MGl2 - 8*MTy2) - 11*MTy2**2 - 
     &   Cd(3132)*MStop2(1)

        Cd(3151) = -2*MStop2(2)**2*
     &    (-2*MTy2**2 + Cd(1468)*MStop2(2)) - 
     &   MTy2*(MTy2**2 - 3*MStop2(2)**2)*MStop2(4)

        Cd(3162) = MGl2*(Cd(48)*Cd(3124) + 
     &      MGl2*(-3*Cd(3123) + MGl2*Cd(3161))) + 
     &   Cd(48)**3*MStop2(2)

        Cd(3170) = -(Cd(48)**2*Cd(49)) + 
     &   MGl2*(MGl2*(MGl2 - 3*Cd(49)) + Cd(3113))

        Cd(3202) = 8*MTy2*Cd(2200)*MStop2(2) + 
     &   2*MStop2(2)**2*(Cd(3127) + 21*MStop2(2)) - 
     &   2*Cd(2682)*MStop2(3)

	Cd(3204) = Cd(2626) - 2*MGlmT2*MStop2(1)

	Cd(3212) = 3*MGl2 - 9*MTy2 - 16*MStop2(2) + 6*MStop2(3)

        Cd(3239) = MTy2*Cd(3235) + 
     &   MStop2(2)*(Cd(3237) + (Cd(3238) - 15*MStop2(2))*MStop2(2))

        Cd(3241) = MGl2**2*(Cd(3232) - MGl2*Cd(3233)) - 
     &   MTy2*Cd(3239) - Cd(3240)*MStop2(2)**2*MStop2(3)

        Cd(3244) = MGl2*(3*MGl2*(MGl2 - 3*Cd(49)) + Cd(3242)) - 
     &   Cd(48)*Cd(3243)

	Cd(320) = 3*MTy2**2 - Cd(875) + 6*Cd(44)*MStop2(2)

        Cd(3251) = MTy2**2*(Cd(529) - 5*MStop2(2)) - 
     &   Cd(3250)*MStop2(2)**2 + MStop2(2)**3

        Cd(3273) = MTy2**2*MStop2(1) - MTy2*Cd(2086)*MStop2(2) - 
     &   4*MStop2(2)**2*(Cd(44) + 2*MStop2(2))

        Cd(3289) = MTy2**2*(Cd(2092) - 3*MStop2(2)) - 
     &   Cd(3288)*MStop2(2)**2 + 15*MStop2(2)**3

        Cd(3295) = MTy2*(MTy2*(4*MTy2 + Cd(637)) - 
     &      2*Cd(3294)*MStop2(2)) - MStop2(2)**2*MStop2(3)

        Cd(3310) = 3*MStop2(2)**2 + 
     &   MTy2*(3*MTy2 - MStop2(1) + 3*MStop2(2))

        Cd(3322) = MGl2*(-(MGl2*Cd(2964)) + Cd(3320)) - 
     &   MTy2*Cd(3321) + 9*MStop2(2)**2*MStop2(3)

        Cd(3327) = MGl2*(17*MGl2 - 9*MTy2) - 5*MTy2**2 - 
     &   Cd(3326)*MStop2(1)

        Cd(3333) = MTy2**2*Cd(160) + 
     &   MStop2(2)*(-(MTy2*Cd(2081)) + 
     &      MStop2(2)*(-Cd(1532) + 5*MStop2(2)))

        Cd(3352) = MTy2**2*Cd(529) + 
     &   MStop2(2)*(-(MTy2*Cd(2151)) + 
     &      MStop2(2)*(-Cd(3351) + 3*MStop2(2)))

        Cd(3379) = (-18*MGl2**3*Cd(843) + 
     &     MStop2(2)*(MGl2**2*Cd(3377) + 
     &        MStop2(2)*(4*MGl2*Cd(3378) - 9*Cd(2201)*MStop2(2))))*
     &   MStop2(3)

        Cd(3390) = MStop2(3)*
     &   (MTy2*(MTy2*(12*MTy2 + Cd(637)) + 2*Cd(3389)*MStop2(2)) - 
     &     9*MStop2(2)**2*MStop2(3))

        Cd(3404) = MTy2*Cd(3401) + 
     &   MStop2(2)*(Cd(3402) + (Cd(3403) - 5*MStop2(2))*MStop2(2))

        Cd(3405) = MGl2**2*(-(MGl2*Cd(2957)) + Cd(3400)) - 
     &   MTy2*Cd(3404) - 3*Cd(3240)*MStop2(2)**2*MStop2(3)

        Cd(3411) = -4*MStop2(2)**3 + 
     &   MTy2*(31*MStop2(2)**2 + MTy2*(3*MTy2 + 2*MStop2(2)))

        Cd(3421) = MTy2*Cd(529) + MGl2*(-3*MGl2 + Cd(1559)) - 
     &   MStop2(2)**2

        Cd(3430) = MTy2*Cd(3426) + 
     &   MStop2(2)*(-Cd(3428) + 
     &      MStop2(2)*(-Cd(3429) + 5*MStop2(2)))

        Cd(3452) = MGl2**2*(Cd(3232) - MGl2*Cd(3233)) - 
     &   MTy2*Cd(3239) + Cd(3451)*MStop2(2)*MStop2(3)

        Cd(3496) = 7*MStop2(1)**3 + 
     &   MTy2*(MTy2*(MTy2 - 7*MStop2(1)) - 9*MStop2(1)**2)

        Cd(3497) = 6*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 12*MStop2(1)) - 13*MStop2(1)**2)

        Cd(3518) = -(Cd(44)**2*Cd(45)) + 
     &   MGl2*(MGl2*(MGl2 - 3*Cd(45)) + Cd(2314))

        Cd(3526) = -5*MStop2(1)**3 + 
     &   MTy2*(6*MStop2(1)**2 + MTy2*(4*MTy2 + 3*MStop2(1)))

        Cd(3527) = 5*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - 7*MStop2(1)) - 8*MStop2(1)**2)

        Cd(3553) = -(Cd(44)*Cd(3080)) + 
     &   MGl2*(3*MGl2*(MGl2 - 3*Cd(45)) + Cd(3081))

        Cd(3559) = 3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 3*MStop2(1)) + 2*MStop2(1)**2)

        Cd(3560) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(MTy2 - 3*MStop2(1)) - 7*MStop2(1)**2)

        Cd(3583) = -3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 7*MStop2(1)) - 4*MStop2(1)**2)

        Cd(3584) = 7*MStop2(1)**3 + 
     &   MTy2*(-MStop2(1)**2 + MTy2*(MTy2 + MStop2(1)))

        Cd(3595) = -3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 5*MStop2(1)) - 6*MStop2(1)**2)

        Cd(3596) = 7*MStop2(1)**3 + 
     &   MTy2*(MTy2*(MTy2 - MStop2(1)) + MStop2(1)**2)

        Cd(3622) = -6*MTy2**2 + MGl2*(2*MTy2 - MStop2(1)) + 
     &   3*MStop2(1)**2

	Cd(3641) = 4*MGl2 - Cd(2767)

	Cd(3663) = Cd(7559) - 6*MTy2*MStop2(2) + Cd(3662)*MStop2(4)

        Cd(3667) = -3*MStop2(1)**3 + 
     &   MTy2*(5*MStop2(1)**2 + MTy2*(3*MTy2 + 3*MStop2(1)))

        Cd(3677) = MStop2(3)*
     &   (MTy2*(MTy2*(12*MTy2 + Cd(1583)) + 
     &        6*Cd(1272)*MStop2(1)) + 9*MStop2(1)**2*MStop2(3))

        Cd(3701) = Cd(3700)*MStop2(1)**2 + 
     &   MGl2*(2*Cd(3694)*MStop2(1) + 
     &      MGl2*(Cd(3697) + 24*MGl2*Cd(652)*UUStop(1,1)))

        Cd(3722) = -4*MStop2(1)**3 + 
     &   MTy2*(15*MStop2(1)**2 + MTy2*(3*MTy2 + 18*MStop2(1)))

        Cd(3738) = 27*MStop2(1)**2 + 
     &   MTy2*(18*MTy2 + 13*MStop2(1)) + Cd(2171)*MStop2(2)

        Cd(3739) = MGl2*(-(MGl2*Cd(2875)) + Cd(2993)) - 
     &   MTy2*Cd(3738) + 6*MStop2(1)**2*MStop2(3)

        Cd(3756) = MTy2*(6*MStop2(1) - 8*MStop2(2)) - 
     &   2*MStop2(2)**2 - Cd(637)*MStop2(3)

        Cd(3759) = 8*MTy2*Cd(3757)*MStop2(1)*MStop2(2) - 
     &   4*MStop2(2)**2*(Cd(3758) + 
     &      2*MStop2(2)*(Cd(2084) + 6*MStop2(2)))

        Cd(3772) = 6*MTy2**2*MStop2(1) + 
     &   6*MTy2*MStop2(1)*MStop2(2) - 
     &   2*MStop2(2)**2*(Cd(3665) + 3*MStop2(2))

        Cd(3797) = 4*MTy2*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(3796) + MStop2(2)*(2*Cd(160) + 9*MStop2(2)))

        Cd(3805) = MGl2*(3*MGl2*(MGl2 - 3*Cd(49)) + Cd(3140)) - 
     &   Cd(48)*Cd(3804)

        Cd(3828) = Cd(3825)*MStop2(1) + 
     &   MStop2(2)*(-Cd(3826) + 
     &      MStop2(2)*(-Cd(3827) + 3*MStop2(2)))

        Cd(3834) = MTy2**2*Cd(160)*MStop2(1) + 
     &   MStop2(2)*(2*MTy2*Cd(3832) + 
     &      MStop2(2)*(-Cd(3833) + Cd(2802)*MStop2(2)))

        Cd(3843) = 2*MTy2*Cd(2086)*MStop2(1) + 
     &   MStop2(2)*(Cd(3841) - 
     &      MStop2(2)*(2*Cd(3842) + 3*MStop2(2)))

	Cd(3900) = 4*MGl2 - Cd(3126)

	Cd(3915) = 3*(MTy2 - MStop2(2))**2 - 2*MStop2(3)**2

        Cd(3918) = 4*MTy2*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(3917) + 
     &      2*MStop2(2)*(Cd(1592) + 3*MStop2(2)))

        Cd(3934) = MStop2(3)*
     &   (MTy2*(MTy2*(12*MTy2 + Cd(637)) + 6*Cd(843)*MStop2(2)) - 
     &     9*MStop2(2)**2*MStop2(3))

        Cd(3958) = Cd(3957)*MStop2(2)**2 + 
     &   MGl2*(-2*Cd(3951)*MStop2(2) + 
     &      MGl2*(Cd(3954) + 12*MGl2*Cd(1160)*UUStop(1,1)))

        Cd(3973) = -4*MStop2(2)**3 + 
     &   MTy2*(15*MStop2(2)**2 + MTy2*(3*MTy2 + 18*MStop2(2)))

        Cd(3994) = MGl2*(-(MGl2*Cd(3233)) + Cd(3348)) - 
     &   MTy2*Cd(3993) - 6*MStop2(2)**2*MStop2(3)

        Cd(3999) = MGl2*(17*MGl2 - 9*MTy2) - 5*MTy2**2 + 
     &   3*Cd(153)*MStop2(1)

	Cd(4074) = -3*MTy2**2 + MGl2*(5*MGl2 + 6*MTy2) - MStop2(1)**2

	Cd(4172) = -3*MTy2**2 + MGl2*(5*MGl2 + 6*MTy2) - MStop2(2)**2

        Cd(4194) = 3*MGl2*Cd(48)**2 + Cd(48)**3 + 
     &   MGl2**2*(MGl2 - Cd(4193))

        Cd(4231) = MTy2*(9*MTy2 - 5*MStop2(1)) - 
     &   2*MStop2(1)**2 - Cd(45)*MStop2(2)

        Cd(4258) = (-(MGl2**3*Cd(594)) + 
     &     MGl2*(10*MGl2 + 3*Cd(594))*MStop2(1)**2 - 
     &     10*Cd(1272)*MStop2(1)**3)*MStop2(3)

        Cd(4354) = MTy2*(6*MTy2 - 3*MStop2(1)) - MStop2(1)**2 - 
     &   Cd(45)*MStop2(2)

        Cd(433) = 8*CA2*MTy - 6*AtC*CA2*UCStop(1,3) - 
     &   3*MUE*S2A*UCStop(1,3) - 3*Cd(432)*UCStopC(1,3)

        Cd(4413) = MTy2*(9*MTy2 - 4*MStop2(1)) - MStop2(1)**2 - 
     &   2*Cd(45)*MStop2(2)

        Cd(4414) = 3*Cd(44)**3 - 
     &   MGl2*(MGl2*(3*MGl2 - Cd(4412)) + Cd(4413))

        Cd(4445) = 4*MTy2*(MTy2 - MStop2(2)) + 
     &   (7*MTy2 - 3*MStop2(1))*MStop2(3)

        Cd(4448) = 3*Cd(67)*Cd(1735) + 
     &   MGl2*(4*MGl2 - Cd(4447)) + Cd(4446)*MStop2(2)

        Cd(437) = 8*CA2*MTy + 6*AtC*CA2*UCStop(1,3) + 
     &   3*MUE*S2A*UCStop(1,3) + 3*Cd(432)*UCStopC(1,3)

	Cd(4466) = 3*MTy2**2 - MGl2*Cd(2200) + Cd(12266)

        Cd(4472) = MTy2*(-2*Cd(4471)*MStop2(1) + 
     &      MTy2*(18*MTy2 - 31*MStop2(3))) + 
     &   3*MStop2(1)**2*MStop2(3)

        Cd(4501) = (-6*MGl2**3*Cd(607) + 
     &     MStop2(1)*(MGl2**2*Cd(4499) + 
     &        MStop2(1)*(MGl2*Cd(4498) - 6*Cd(4500)*MStop2(1))))*
     &   MStop2(3)

        Cd(4514) = MStop2(3)*
     &   (MTy2*(MTy2*(12*MTy2 - Cd(4513)) - 
     &        2*Cd(4512)*MStop2(1)) + 13*MStop2(1)**2*MStop2(3))

        Cd(4521) = MGl2*(5*MGl2*(MGl2 - Cd(1592)) - Cd(4475)) + 
     &   Cd(44)**2*Cd(4520)

        Cd(4525) = Cd(44)**3 + MGl2**2*(2*MGl2 - Cd(3252)) + 
     &   MGl2*Cd(44)*Cd(3770)

        Cd(4526) = MTy2*Cd(2626) - 
     &   3*MStop2(1)*(2*MTy2**2 + (MGlpT2 - 2*MStop2(1))*MStop2(1))

        Cd(4606) = (-(MGl2**3*Cd(874)) + 
     &     MGl2*(8*MGl2 + 3*Cd(874))*MStop2(2)**2 - 
     &     2*Cd(1421)*MStop2(2)**3)*MStop2(3)

	Cd(4609) = -Cd(48)**3 + MGl2*(MGl2*(MGl2 - Cd(2254)) + Cd(4608))

        Cd(4646) = (Cd(4645) - 9*MStop2(2))*MStop2(2)**2 + 
     &   MTy2*Cd(2617)*(4*MStop2(2) - MStop2(3))

	Cd(457) = Cd(429) + Cd(428)*UStop2(1,1)*UStop2(1,2)

        Cd(4712) = MTy2*(3*MTy2 + 4*MStop2(1) - 2*MStop2(2)) + 
     &   3*MStop2(1)*MStop2(2)

	Cd(4781) = (MTy2 - MStop2(2))*(4*MTy2 - 3*MStop2(3))

        Cd(4788) = 4*MTy2*(MTy2 - MStop2(2)) + 
     &   (7*MTy2 + MStop2(2))*MStop2(3)

	Cd(4797) = Cd(2626) - Cd(4467)*MStop2(1)

        Cd(4799) = Cd(2626)*MStop2(1) + 
     &   MStop2(2)*(Cd(4797) + MStop2(2)*(-Cd(4798) + MStop2(2)))

        Cd(4831) = (6*MGl2**3*Cd(843) + 
     &     MStop2(2)*(-(MGl2**2*Cd(4830)) + 
     &        MStop2(2)*(-(MGl2*Cd(4829)) + 6*Cd(1357)*MStop2(2))))
     &    *MStop2(3)

        Cd(4846) = -(MStop2(3)*
     &     (MTy2*(MTy2*(-12*MTy2 + Cd(4845)) + 
     &          2*Cd(4844)*MStop2(2)) + 13*MStop2(2)**2*MStop2(3)))

        Cd(4854) = 3*Cd(48)**3 + 
     &   MGl2*(MGl2*(3*MGl2 + Cd(4852)) - Cd(4853))

        Cd(4864) = Cd(48)**3 + MGl2*Cd(48)*Cd(4862) + 
     &   MGl2**2*(2*MGl2 - Cd(4863))

        Cd(4865) = MTy2*Cd(2626) - 
     &   3*MStop2(2)*(2*MTy2**2 + (MGlpT2 - 2*MStop2(2))*MStop2(2))

        Cd(4903) = -6*MTy2**3*Cd(874) + 
     &   MStop2(2)*(MTy2**2*Cd(4507) + 
     &      MStop2(2)*(Cd(1148)*MStop2(2) + 4*MTy2*MStop2(3)))

        Cd(4908) = 2*MTy2*Cd(2617)*MStop2(1) + 
     &   MStop2(2)*(Cd(4906) + Cd(4907)*MStop2(2) + MStop2(2)**2)

        Cd(4922) = 80*MTy2**3 + 4*MTy2*Cd(4919) + 
     &   2*(12*MTy2**2 + Cd(4921))*MStop2(4)

        Cd(4923) = -3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(8*MTy2 - 18*MStop2(1)) + 9*MStop2(1)**2)

        Cd(4926) = -7*MStop2(1)**3 + 
     &   MTy2*(21*MStop2(1)**2 + MTy2*(12*MTy2 + 46*MStop2(1)))

        Cd(4931) = 6*MStop2(1)**3 + 
     &   MTy2*(-12*MStop2(1)**2 + MTy2*(14*MTy2 + 21*MStop2(1)))

        Cd(4946) = MGl2*(MGl2*(3*MGl2 - 8*Cd(529)) + Cd(4944)) - 
     &   2*Cd(4945)

        Cd(4960) = 12*MStop2(1)**3 + 
     &   MTy2*(MTy2*(9*MTy2 - 9*MStop2(1)) - 32*MStop2(1)**2)

        Cd(4966) = 7*MStop2(1)**3 + 
     &   MTy2*(MTy2*(13*MTy2 - 13*MStop2(1)) - 31*MStop2(1)**2)

        Cd(4970) = -17*MStop2(1)**3 + 
     &   MTy2*(18*MStop2(1)**2 + MTy2*(14*MTy2 + 81*MStop2(1)))

        Cd(4981) = 60*MStop2(2)**2 - Cd(4980)*MStop2(3) + 
     &   MTy2*(60*MTy2 + 36*MStop2(4))

	Cd(494) = Cd(10191) + 4*(Cd(493) + At*MTy*UCStopC(1,3))

        Cd(5036) = -17*MStop2(1)**3 + 
     &   MTy2*(39*MStop2(1)**2 + MTy2*(3*MTy2 + 8*MStop2(1)))

        Cd(5041) = -29*MStop2(1)**3 + 
     &   MTy2*(60*MStop2(1)**2 + MTy2*(3*MTy2 + 69*MStop2(1)))

        Cd(5049) = -4*MStop2(1)**3 + 
     &   MTy2*(MStop2(1)**2 + MTy2*(4*MTy2 + 19*MStop2(1)))

        Cd(5060) = -8*MStop2(1)**3 + 
     &   MTy2*(14*MStop2(1)**2 + MTy2*(11*MTy2 + 39*MStop2(1)))

        Cd(5066) = 4*MStop2(1)**3 + 
     &   MTy2*(-10*MStop2(1)**2 + MTy2*(12*MTy2 + 23*MStop2(1)))

        Cd(499) = MUE*Cd(498)*UCStop(1,3) + 
     &   MUEC*UCStopC(1,3)*(2*MTy + At*UCStopC(1,3)) + 
     &   Cd(424)*UStop2(1,1)*UStop2(1,2)

        Cd(5082) = 6*Cd(5081) + 54*MTy2*(MTy2 + 2*MStop2(2)) - 
     &   88*MTy2*MStop2(3)

        Cd(5085) = 19*MStop2(1)**3 + 
     &   MTy2*(40*MStop2(1)**2 + MTy2*(15*MTy2 + 50*MStop2(1)))

        Cd(5088) = 2*(Cd(5085) - 
     &     MStop2(2)*(Cd(5086) + 
     &        MStop2(2)*(Cd(5087) + 4*MStop2(2))))

	Cd(501) = Cd(10191) + 4*At*MTy*UCStopC(1,3)

        Cd(5089) = -12*MStop2(1)**3 + 
     &   MTy2*(27*MStop2(1)**2 + MTy2*(13*MTy2 + 6*MStop2(1)))

        Cd(5094) = 20*MStop2(1)**3 + 
     &   MTy2*(MTy2*(5*MTy2 - 23*MStop2(1)) - 45*MStop2(1)**2)

        Cd(5101) = -3*MStop2(1)**3 + 
     &   MTy2*(5*MStop2(1)**2 + MTy2*(6*MTy2 + 12*MStop2(1)))

        Cd(5107) = 10*MTy2**3 + 3*MTy2*Cd(523) + 
     &   (3*MTy2**2 + Cd(831))*MStop2(4)

	Cd(5119) = 3*(MGl2 - 5*MTy2 - 7*MStop2(2)) + 14*MStop2(3)

        Cd(5132) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 7*MStop2(1)) + 2*MStop2(1)**2)

        Cd(5134) = 2*(Cd(5132) - 
     &     MStop2(2)*(Cd(5133) - 6*MTy2*MStop2(2) + 2*MStop2(2)**2)
     &     )

        Cd(5159) = 2*(Cd(5158) + 
     &     MTy2*(MTy2*(-4*MTy2 + 3*Cd(874)) - 
     &        6*MStop2(2)*MStop2(4)))

	Cd(5163) = -16*MTy2 - Cd(1572)

        Cd(5188) = 3*Cd(5183) + 63*MGl2*(MGl2 + 2*MStop2(2)) - 
     &   76*MGl2*MStop2(3)

        Cd(5194) = MGl2*(MGl2*(35*MGl2 + Cd(5189)) + Cd(5193)) + 
     &   5*Cd(5191)*MStop2(1)

        Cd(5213) = (-(Cd(595)*MStop2(2)**2) + 
     &     MTy2*(MTy2*(MTy2 + Cr(2)) + Cd(5211)*MStop2(2)))*
     &   MStop2(3)

        Cd(5231) = 3*Cd(1163) + 
     &   2*(MTy2*(9*MTy2 + 5*MStop2(1)) + Cd(5230)*MStop2(2))

        Cd(5237) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 8*MStop2(1)) + 9*MStop2(1)**2)

        Cd(5241) = -10*MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 75*MStop2(1)) + 60*MStop2(1)**2)

        Cd(5248) = 4*MStop2(1)**3 + 
     &   MTy2*(-22*MStop2(1)**2 + MTy2*(9*MTy2 + 29*MStop2(1)))

        Cd(5251) = 60*MStop2(2)**2 - Cd(4936)*MStop2(3) + 
     &   MTy2*(60*MTy2 + 36*MStop2(4))

        Cd(5267) = -23*MStop2(1)**3 + 
     &   MTy2*(53*MStop2(1)**2 + MTy2*(36*MTy2 + 110*MStop2(1)))

        Cd(5272) = 14*MStop2(1)**3 + 
     &   MTy2*(-28*MStop2(1)**2 + MTy2*(38*MTy2 + 53*MStop2(1)))

	Cd(5286) = Cd(1227) + 6*MTy2*(MTy2 - MStop2(4))

        Cd(5288) = -90*MStop2(2)**2 + Cd(5287)*MStop2(3) + 
     &   6*MTy2*(MTy2 - MStop2(4))

        Cd(5292) = 5*MStop2(1)**3 + 
     &   MTy2*(-27*MStop2(1)**2 + MTy2*(3*MTy2 + 46*MStop2(1)))

        Cd(5298) = 4*MStop2(1)**3 + 
     &   MTy2*(-24*MStop2(1)**2 + MTy2*(5*MTy2 + 27*MStop2(1)))

        Cd(5322) = 20*MTy2**3 + 2*MTy2*Cd(5321) - 
     &   2*(11*MTy2**2 + MStop2(1)*MStop2(2))*MStop2(4)

        Cd(5323) = -4*MStop2(1)**3 + 
     &   MTy2*(34*MStop2(1)**2 + MTy2*(25*MTy2 + 20*MStop2(1)))

        Cd(5329) = -9*MStop2(1)**3 + 
     &   MTy2*(MTy2*(29*MTy2 - 58*MStop2(1)) + 54*MStop2(1)**2)

        Cd(5335) = 2*MStop2(1)**3 + 
     &   MTy2*(-11*MStop2(1)**2 + MTy2*(13*MTy2 + 28*MStop2(1)))

        Cd(5338) = 135*MTy2**3 + 
     &   MGl2*(-55*MTy2**2 + MGl2*(MGl2 + 7*MTy2))

        Cd(5341) = 2*Cd(5338) - 
     &   MStop2(1)*(Cd(5339) + MStop2(1)*(-Cd(5340) + 5*MStop2(1)))

	Cd(5351) = Cd(8360) + MTy2*(62*MTy2 + 34*MStop2(4))

        Cd(5354) = 22*MTy2**3 + 9*MTy2*Cd(523) + 
     &   (7*MTy2**2 + Cd(831))*MStop2(4)

        Cd(5372) = Cd(1272)*Cd(1280) + 
     &   MTy2*(-Cd(523) + MTy2*(2*MTy2 - MStop2(4)))

        Cd(5393) = 78*MStop2(2)**2 - Cd(5392)*MStop2(3) + 
     &   MTy2*(126*MTy2 + 74*MStop2(4))

        Cd(5423) = Cd(1272)*Cd(5420) + 
     &   MTy2*(Cd(5422) + MTy2*(210*MTy2 + 107*MStop2(4)))

        Cd(5431) = 10*MStop2(1)**3 + 
     &   MTy2*(MTy2*(16*MTy2 - 20*MStop2(1)) + 15*MStop2(1)**2)

        Cd(5436) = 34*MStop2(1)**3 + 
     &   MTy2*(23*MStop2(1)**2 + MTy2*(24*MTy2 + 20*MStop2(1)))

        Cd(5446) = 40*MStop2(1)**3 + 
     &   MTy2*(-120*MStop2(1)**2 + MTy2*(17*MTy2 + 26*MStop2(1)))

        Cd(52) = 3*MGlpT2 - 2*MGl2*Log(MGl2/Q) - 
     &   2*MTy2*Log(MTy2/Q) + Cd(51)*MStop2(1)

        Cd(5473) = -9*MStop2(1)**3 + 
     &   MTy2*(18*MStop2(1)**2 + MTy2*(67*MTy2 + 6*MStop2(1)))

        Cd(5476) = MTy2*Cd(5473) + 
     &   MStop2(2)*(Cd(5475)*MStop2(2) + 
     &      Cd(45)*(Cd(5474) - 9*MStop2(2)**2))

        Cd(5479) = 10*MStop2(1)**3 + 
     &   MTy2*(-47*MStop2(1)**2 + MTy2*(29*MTy2 + 92*MStop2(1)))

        Cd(5480) = -9*MStop2(1)**3 + 
     &   MTy2*(10*MStop2(1)**2 + MTy2*(36*MTy2 + 47*MStop2(1)))

        Cd(5482) = 2*MTy2**2*Cd(5477)*Cd(5478) - 
     &   2*MStop2(2)*(MTy2*Cd(5479) + 
     &      MStop2(2)*(-Cd(5480) + Cd(5481)*MStop2(2)))

        Cd(5491) = 5*MStop2(4)**3 + 
     &   MTy2*(Cd(5490) + MTy2*(40*MTy2 + 12*MStop2(4)))

        Cd(5493) = -18*MStop2(1)**3 + 
     &   MTy2*(MTy2*(16*MTy2 - 28*MStop2(1)) + 9*MStop2(1)**2)

        Cd(5498) = 17*MStop2(1)**3 + 
     &   MTy2*(MTy2*(60*MTy2 - 24*MStop2(1)) + 42*MStop2(1)**2)

        Cd(5499) = -66*MStop2(1)**3 + 
     &   MTy2*(41*MStop2(1)**2 + MTy2*(24*MTy2 + 108*MStop2(1)))

        Cd(5502) = MTy2*Cd(5498) + 
     &   MStop2(2)*(-Cd(5499) + 
     &      MStop2(2)*(Cd(5500) + Cd(5501)*MStop2(2)))

        Cd(5505) = -24*MStop2(1)**3 + 
     &   MTy2*(88*MStop2(1)**2 + MTy2*(31*MTy2 + 38*MStop2(1)))

        Cd(5527) = 3*Cd(5526) + MTy2*(27*MTy2 + 10*MStop2(1)) + 
     &   4*Cd(5525)*MStop2(2)

        Cd(5539) = -7*MStop2(1)**3 + 
     &   MTy2*(MTy2*(8*MTy2 - 53*MStop2(1)) + 45*MStop2(1)**2)

        Cd(5546) = 3*MStop2(1)**3 + 
     &   MTy2*(-16*MStop2(1)**2 + MTy2*(9*MTy2 + 24*MStop2(1)))

        Cd(5553) = 20*MStop2(2)**2 - Cd(5552)*MStop2(3) + 
     &   MTy2*(20*MTy2 - 4*MStop2(4))

	Cd(5562) = 3*(MGl2 - 5*MTy2 - 7*MStop2(2)) + 7*MStop2(3)

        Cd(5576) = 60*MStop2(2)**2 - Cd(4980)*MStop2(3) + 
     &   4*MTy2*(7*MTy2 + 5*MStop2(4))

        Cd(5577) = 2*MTy2*(8*MTy2**2 + 7*MStop2(3)**2) + 
     &   (4*MTy2**2 + Cd(5429))*MStop2(4)

        Cd(5579) = -7*MStop2(1)**3 + 
     &   MTy2*(MTy2*(4*MTy2 - 38*MStop2(1)) + 31*MStop2(1)**2)

        Cd(5597) = -42*MStop2(2)**2 + 3*Cd(5596)*MStop2(3) + 
     &   MTy2*(98*MTy2 + 64*MStop2(4))

        Cd(5598) = 32*MTy2**3 - 
     &   2*(17*MTy2**2 + MStop2(1)*MStop2(2))*MStop2(4) + 
     &   10*MTy2*MStop2(4)**2

        Cd(5602) = 336*MTy2**3 + 2*MTy2*Cd(5600) + 
     &   2*(117*MTy2**2 - 5*Cd(5601))*MStop2(4)

        Cd(5603) = -8*MStop2(1)**3 + 
     &   MTy2*(MTy2*(20*MTy2 - 44*MStop2(1)) + 95*MStop2(1)**2)

        Cd(5609) = 9*MStop2(1)**3 + 
     &   MTy2*(-81*MStop2(1)**2 + MTy2*(13*MTy2 + 127*MStop2(1)))

        Cd(5615) = 4*MStop2(1)**3 + 
     &   MTy2*(-31*MStop2(1)**2 + MTy2*(29*MTy2 + 62*MStop2(1)))

        Cd(5619) = 2*(MTy2*Cd(5473) + 
     &     MStop2(2)*(Cd(5475)*MStop2(2) + 
     &        Cd(45)*(Cd(5474) - 9*MStop2(2)**2)))

        Cd(5645) = 11*MStop2(1)**3 + 
     &   MTy2*(-35*MStop2(1)**2 + MTy2*(12*MTy2 + 22*MStop2(1)))

        Cd(5646) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(58*MTy2 - 24*MStop2(1)) - 5*MStop2(1)**2)

        Cd(5649) = MTy2*Cd(5645) - 
     &   MStop2(2)*(Cd(5646) + 
     &      MStop2(2)*(-(MTy2*Cd(5647)) + Cd(5648)*MStop2(2)))

        Cd(5650) = 25*MStop2(1)**3 + 
     &   MTy2*(-37*MStop2(1)**2 + MTy2*(390*MTy2 + 48*MStop2(1)))

        Cd(5655) = 10*MStop2(1)**3 + 
     &   MTy2*(MTy2*(50*MTy2 - 387*MStop2(1)) + 108*MStop2(1)**2)

        Cd(5675) = 52*MStop2(1)**3 + 
     &   MTy2*(MTy2*(89*MTy2 - 922*MStop2(1)) + 360*MStop2(1)**2)

        Cd(5680) = 60*MStop2(2)**2 - Cd(5464)*MStop2(3) + 
     &   MTy2*(60*MTy2 + 36*MStop2(4))

        Cd(5725) = MTy2*(MTy2*(-5*MTy2 + Cd(5211)) + Cd(5723)) - 
     &   Cd(5724)*MStop2(1)*MStop2(2)

        Cd(5759) = 7*MStop2(1)**3 + 
     &   MTy2*(-35*MStop2(1)**2 + MTy2*(6*MTy2 + 8*MStop2(1)))

        Cd(5763) = 3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(5*MTy2 - 3*MStop2(1)) - 9*MStop2(1)**2)

        Cd(5765) = MStop2(1)**4 + 
     &   MTy2*(-36*MStop2(1)**3 + 
     &      MTy2*(18*MTy2**2 + 70*MStop2(1)**2))

        Cd(5766) = -6*MStop2(1)**3 + 
     &   MTy2*(-2*MStop2(1)**2 + MTy2*(6*MTy2 + 43*MStop2(1)))

        Cd(5787) = -MStop2(1)**3 + 
     &   MTy2*(MTy2*(5*MTy2 - 41*MStop2(1)) + 3*MStop2(1)**2)

	Cd(5797) = Cd(6534) + MGl2*(MGl2 - 3*MStop2(4))

        Cd(5808) = Cd(3852)*MStop2(4) - 
     &   2*MTy2*(Cd(5807) + MTy2*(-4*MTy2 + 2*MStop2(4)))

	Cd(5818) = Cd(6534) + MTy2*(7*MTy2 - 3*MStop2(4))

        Cd(5827) = Cd(1272)*Cd(5826) - 
     &   2*MTy2*(Cd(835) + MTy2*(-14*MTy2 + MStop2(4)))

        Cd(5842) = 42*MStop2(2)**2 - Cd(5841)*MStop2(3) + 
     &   MTy2*(18*MTy2 + 18*MStop2(4))

        Cd(5845) = 20*MTy2**3 + MTy2*Cd(5843) + 
     &   (18*MTy2**2 + 7*Cd(5844))*MStop2(4)

        Cd(5846) = -3*MStop2(1)**3 + 
     &   MTy2*(5*MStop2(1)**2 + MTy2*(MTy2 + 6*MStop2(1)))

        Cd(5851) = 26*MStop2(1)**3 + 
     &   MTy2*(-92*MStop2(1)**2 + MTy2*(12*MTy2 + 21*MStop2(1)))

        Cd(5858) = -3*MStop2(1)**3 + 
     &   MTy2*(10*MStop2(1)**2 + MTy2*(MTy2 + 2*MStop2(1)))

        Cd(56) = 2*MGlmT2**2 - 3*MGlpT2*MStop2(1) + 
     &   MStop2(1)**2 - MGlpTmSt2(1)*MStop2(2)

        Cd(5866) = -2*MStop2(1)**3 + 
     &   MTy2*(6*MStop2(1)**2 + MTy2*(3*MTy2 + 16*MStop2(1)))

        Cd(5890) = 42*MStop2(2)**2 - Cd(5889)*MStop2(3) + 
     &   MTy2*(18*MTy2 + 18*MStop2(4))

	Cd(5919) = Cd(6534) + MTy2*(13*MTy2 + 2*MStop2(4))

        Cd(5922) = 64*MTy2**3 + 
     &   2*(2*MTy2**2 + Cd(5921))*MStop2(4) - 12*MTy2*MStop2(4)**2

        Cd(5924) = -5*MStop2(1)**3 + 
     &   MTy2*(MTy2*(8*MTy2 - 2*MStop2(1)) + 23*MStop2(1)**2)

        Cd(5936) = -MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 6*MStop2(1)) + MStop2(1)**2)

        Cd(5941) = -2*MStop2(1)**3 + 
     &   MTy2*(MTy2*(MTy2 - MStop2(1)) + 14*MStop2(1)**2)

	Cd(5943) = Cd(5559) + MTy2*(9*MTy2 + 3*MStop2(4))

	Cd(5944) = Cd(5112) - 3*MTy2*MStop2(4)

        Cd(5946) = 24*MTy2**3 - 24*MTy2*Cd(1231) + 
     &   2*(-6*MTy2**2 + Cd(5945))*MStop2(4)

	Cd(5956) = Cd(1227) + 2*MTy2*(11*MTy2 + 5*MStop2(4))

        Cd(5959) = MTy2*(Cd(5958) + 
     &      MTy2*(6*MTy2 - 5*MStop2(4))) + 
     &   MStop2(1)*MStop2(2)*MStop2(4)

	Cd(5961) = Cd(1280) + MTy2*(MTy2 - MStop2(4))

	Cd(5968) = Cd(4243) + Cd(256)*MStop2(3) - MTy2*MStop2(4)

	Cd(5995) = Cd(831) + 5*MTy2*(MTy2 - MStop2(4))

        Cd(5999) = MTy2*(-12*Cd(5996) + 
     &      MTy2*(4*MTy2 - 2*MStop2(4))) + Cd(5998)*MStop2(4)

        Cd(6001) = 17*MStop2(1)**3 + 
     &   MTy2*(-85*MStop2(1)**2 + MTy2*(2*MTy2 + 130*MStop2(1)))

        Cd(6006) = -2*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - MStop2(1)) + 3*MStop2(1)**2)

        Cd(6008) = 13*MStop2(1)**3 + 
     &   MTy2*(-74*MStop2(1)**2 + MTy2*(24*MTy2 + 85*MStop2(1)))

        Cd(6019) = 3*Cd(6018) + MTy2*(3*MTy2 - 2*MStop2(1)) + 
     &   4*Cd(529)*MStop2(2)

        Cd(6020) = MStop2(1)**3 + 
     &   MTy2*(-4*MStop2(1)**2 + MTy2*(MTy2 + 3*MStop2(1)))

        Cd(6022) = 2*(Cd(6020) - 
     &     MStop2(2)*(Cd(6021) + MStop2(2)*(-Cd(1852) + MStop2(2)))
     &     )

        Cd(6027) = -6*MStop2(1)**3 + 
     &   MTy2*(MTy2*(4*MTy2 - 2*MStop2(1)) + 21*MStop2(1)**2)

        Cd(6032) = 13*MStop2(1)**3 + 
     &   MTy2*(-62*MStop2(1)**2 + MTy2*(12*MTy2 + 10*MStop2(1)))

        Cd(6033) = 5*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - 77*MStop2(1)) + 15*MStop2(1)**2)

        Cd(6039) = MStop2(1)**3 + 
     &   MTy2*(-14*MStop2(1)**2 + MTy2*(7*MTy2 + 26*MStop2(1)))

	Cd(6041) = -5*MTy2 + Cd(10678) - 7*MStop2(2)

        Cd(6072) = 3*Cd(6058) + 
     &   MGl2*(9*MGl2 + 2*(9*MStop2(2) + MStop2(3)))

        Cd(6124) = -MStop2(1)**3 + 
     &   MTy2*(MTy2*(4*MTy2 - 5*MStop2(1)) + MStop2(1)**2)

        Cd(6129) = -MStop2(1)**3 + 
     &   MTy2*(12*MStop2(1)**2 + MTy2*(16*MTy2 + 10*MStop2(1)))

        Cd(6139) = -MStop2(1)**3 + 
     &   MTy2*(-17*MStop2(1)**2 + MTy2*(14*MTy2 + 32*MStop2(1)))

        Cd(6143) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(10*MTy2 - 12*MStop2(1)) - 3*MStop2(1)**2)

        Cd(6145) = -3*MStop2(1)**4 + 
     &   MTy2*(-4*MStop2(1)**3 + 
     &      MTy2*(42*MTy2**2 + 34*MStop2(1)**2))

        Cd(6146) = -2*MStop2(1)**3 + 
     &   MTy2*(-8*MStop2(1)**2 + MTy2*(14*MTy2 + 49*MStop2(1)))

        Cd(6170) = 6*MTy2*Cd(286) + 
     &   15*(MTy2**2 + MStop2(1)**2) + Cd(6169)*MStop2(3)

        Cd(6175) = -3*MTy2**3 + MTy2*Cd(6174) + 
     &   (3*MTy2**2 - Cd(535)*MStop2(1))*MStop2(2)

        Cd(6180) = 21*MStop2(2)**2 - Cd(6179)*MStop2(3) + 
     &   MTy2*(9*MTy2 + 9*MStop2(4))

        Cd(6183) = 40*MTy2*(MTy2**2 + 3*MStop2(1)*MStop2(2)) + 
     &   2*(18*MTy2**2 + 5*Cd(6182))*MStop2(4)

        Cd(6196) = -MStop2(1)**3 + 
     &   MTy2*(MTy2*(5*MTy2 - 42*MStop2(1)) + 18*MStop2(1)**2)

        Cd(6197) = -6*MStop2(1)**3 + 
     &   MTy2*(18*MStop2(1)**2 + MTy2*(5*MTy2 + 3*MStop2(1)))

        Cd(6212) = 2*Cd(5921)*MStop2(4) - 
     &   4*MTy2*(Cd(835) + MTy2*(-14*MTy2 + MStop2(4)))

        Cd(6216) = -5*MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 22*MStop2(1)) + 32*MStop2(1)**2)

        Cd(6222) = MStop2(1)**3 + 
     &   MTy2*(-16*MStop2(1)**2 + MTy2*(34*MTy2 + 59*MStop2(1)))

        Cd(612) = AtC*(MUEC*Cd(601) - 
     &      4*MUE*Cd(609)*UCStop(1,3)**2) + 
     &   At*(2*MUE*Cd(606) - 4*MUEC*Cd(611)*UCStopC(1,3)**2)

        Cd(6243) = -2*MTy2**3 + MTy2*Cd(6241) + 
     &   (MTy2**2 - Cd(6242))*MStop2(4)

	Cd(60) = Cd(12162) + Cd(56)*UStop2(1,1)**2

        Cd(6273) = (Cd(1272)*Cd(6259)*Cd(6268) + 
     &     2*MGl2*(MGl2*(Cd(843)*Cd(6269) + 2*MGl2*Cd(6270)) - 
     &        2*Cd(6272)))*MStop2(3)

        Cd(6288) = -(Cd(1272)*Cd(5062)) + 
     &   MTy2*(4*Cd(832) + MTy2*(4*MTy2 - 6*MStop2(4)))

        Cd(6301) = -4*MTy2*Cd(6300) + Cd(6298)*MStop2(4) + 
     &   2*MTy2**2*(10*MTy2 + 9*MStop2(4))*UStop2(1,1)

        Cd(621) = UStop2(1,1)**3*
     &    (Cd(597) - 2*Cd(607)*UStop2(1,1)) - 
     &   2*Cd(537)*MStop2(2)*UStop2(1,2)**2 + 
     &   4*Cd(286)*UStop2(1,1)*UStop2(1,2)**3

        Cd(6357) = -3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - 46*MStop2(1)) + 21*MStop2(1)**2)

	Cd(6382) = Cd(8052) - 4*MGlpT2*MStop2(2) + MStop2(2)*MStop2(4)

        Cd(6397) = -10*MStop2(1)**3 + 
     &   MTy2*(40*MStop2(1)**2 + MTy2*(35*MTy2 + 88*MStop2(1)))

        Cd(6398) = 3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(16*MTy2 - 19*MStop2(1)) + 6*MStop2(1)**2)

        Cd(6404) = -2*MStop2(1)**3 + 
     &   MTy2*(-2*MStop2(1)**2 + MTy2*(9*MTy2 + 34*MStop2(1)))

        Cd(6408) = -MStop2(1)**4 + 
     &   MTy2**2*(MTy2*(15*MTy2 - 21*MStop2(1)) + 2*MStop2(1)**2)

        Cd(6414) = MGl2*(MGl2 - 2*Cd(6413)) + 15*MStop2(2)**2 - 
     &   Cd(5441)*MStop2(3)

	Cd(6438) = Cd(832) + 3*MTy2*(MTy2 - MStop2(4))

        Cd(6443) = -7*MStop2(1)**3 + 
     &   MTy2*(11*MStop2(1)**2 + MTy2*(26*MTy2 + 26*MStop2(1)))

        Cd(6449) = 5*MStop2(1)**3 + 
     &   MTy2*(-18*MStop2(1)**2 + MTy2*(28*MTy2 + 33*MStop2(1)))

	Cd(6458) = -6*MGl2*MTy2 + Cd(58)

        Cd(6466) = 19*MStop2(1)**3 + 
     &   MTy2*(MTy2*(7*MTy2 - 7*MStop2(1)) - 37*MStop2(1)**2)

        Cd(6472) = 14*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 22*MStop2(1)) - 19*MStop2(1)**2)

        Cd(6477) = -39*MStop2(1)**3 + 
     &   MTy2*(102*MStop2(1)**2 + MTy2*(6*MTy2 + 70*MStop2(1)))

        Cd(6481) = 6*MTy2*Cd(1226) + 
     &   30*(MTy2**2 + MStop2(2)**2) - Cd(5584)*MStop2(3)

        Cd(6484) = 7*MStop2(1)**3 + 
     &   MTy2*(-37*MStop2(1)**2 + MTy2*(10*MTy2 + 16*MStop2(1)))

        Cd(6487) = 3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(8*MTy2 - 6*MStop2(1)) - 9*MStop2(1)**2)

        Cd(6489) = MStop2(1)**4 + 
     &   MTy2*(-36*MStop2(1)**3 + 
     &      MTy2*(30*MTy2**2 + 74*MStop2(1)**2))

        Cd(6490) = -6*MStop2(1)**3 + 
     &   MTy2*(-4*MStop2(1)**2 + MTy2*(10*MTy2 + 53*MStop2(1)))

	Cd(6492) = 6*MGl2 - Cd(6506)

        Cd(6515) = 19*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - 19*MStop2(1)) - 37*MStop2(1)**2)

        Cd(6520) = -37*MStop2(1)**3 + 
     &   MTy2*(94*MStop2(1)**2 + MTy2*(5*MTy2 + 55*MStop2(1)))

        Cd(6528) = 10*MStop2(1)**3 + 
     &   MTy2*(MTy2*(5*MTy2 - 12*MStop2(1)) - 51*MStop2(1)**2)

        Cd(6533) = 14*MStop2(1)**2 + 
     &   MTy2*(30*MTy2 + 37*MStop2(1)) - Cd(6532)*MStop2(2)

        Cd(6539) = -MStop2(1)**3 + 
     &   MTy2*(6*MStop2(1)**2 + MTy2*(10*MTy2 + 6*MStop2(1)))

        Cd(6541) = MStop2(1)**3 + 
     &   MTy2*(-5*MStop2(1)**2 + MTy2*(2*MTy2 + MStop2(1)))

        Cd(6542) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 2*MStop2(1)) - 3*MStop2(1)**2)

        Cd(6544) = MTy2*Cd(6541) + 
     &   MStop2(2)*(-3*Cd(6542) + 
     &      MStop2(2)*(MTy2*Cd(6543) - Cd(4362)*MStop2(2)))

        Cd(6546) = 15*MStop2(1)**4 + 
     &   MTy2*(-24*MStop2(1)**3 + 
     &      MTy2*(40*MTy2**2 + 50*MStop2(1)**2))

        Cd(6547) = 16*MStop2(1)**3 + 
     &   MTy2*(-8*MStop2(1)**2 + MTy2*(37*MTy2 + 82*MStop2(1)))

        Cd(6555) = 5*Cd(6554) + MTy2*(27*MTy2 + 34*MStop2(1)) + 
     &   2*Cd(6553)*MStop2(2)

        Cd(6565) = -9*MStop2(1)**4 + 
     &   MTy2*(60*MStop2(1)**3 + 
     &      MTy2*(5*MTy2**2 - 88*MStop2(1)**2))

        Cd(6566) = -8*MStop2(1)**3 + 
     &   MTy2*(13*MStop2(1)**2 + MTy2*(3*MTy2 + 16*MStop2(1)))

        Cd(6570) = 21*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 34*MStop2(1)) - 57*MStop2(1)**2)

        Cd(6600) = 3*(MTy2 - MStop2(1))*(MTy2 + 3*MStop2(1)) - 
     &   Cd(1583)*MStop2(3)

	Cd(6604) = 3*(MTy2 - MStop2(1))**2 - MStop2(3)**2

        Cd(6609) = -8*MStop2(1)**3 + 
     &   MTy2*(13*MStop2(1)**2 + MTy2*(9*MTy2 + MStop2(1)))

        Cd(6615) = -16*MStop2(1)**3 + 
     &   MTy2*(43*MStop2(1)**2 + MTy2*(11*MTy2 + 10*MStop2(1)))

        Cd(6633) = -MStop2(1)**3 + 
     &   MTy2*(MTy2*(8*MTy2 - 34*MStop2(1)) + 36*MStop2(1)**2)

        Cd(6634) = 6*MStop2(1)**3 + 
     &   MTy2*(-16*MStop2(1)**2 + MTy2*(22*MTy2 + 41*MStop2(1)))

        Cd(64) = 3*MGlpT2 - 2*MGl2*Log(MGl2/Q) - 
     &   2*MTy2*Log(MTy2/Q) + Cd(63)*MStop2(2)

	Cd(6680) = -6*MGl2*MTy2 + Cd(71)

        Cd(6682) = 12*MTy2*Cd(286) + 
     &   30*(MTy2**2 + MStop2(2)**2) - Cd(6681)*MStop2(3)

        Cd(6689) = -MStop2(1)**3 + 
     &   MTy2*(MTy2*(8*MTy2 - 50*MStop2(1)) + 36*MStop2(1)**2)

        Cd(6690) = 6*MStop2(1)**3 + 
     &   MTy2*(-8*MStop2(1)**2 + MTy2*(6*MTy2 + MStop2(1)))

        Cd(6697) = -2*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 16*MStop2(1)) + 22*MStop2(1)**2)

        Cd(6756) = -9*MStop2(1)**3 + 
     &   MTy2*(7*MStop2(1)**2 + MTy2*(45*MTy2 + 19*MStop2(1)))

        Cd(6760) = 20*MStop2(1)**3 + 
     &   MTy2*(-94*MStop2(1)**2 + MTy2*(26*MTy2 + 159*MStop2(1)))

        Cd(6766) = -3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 10*MStop2(1)) + 20*MStop2(1)**2)

        Cd(6776) = 2*MStop2(1)**3 + 
     &   MTy2*(MTy2*(9*MTy2 - 14*MStop2(1)) + 9*MStop2(1)**2)

        Cd(6780) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 3*MStop2(1)) - 2*MStop2(1)**2)

        Cd(6781) = MTy2**2*Cd(6779) + 
     &   MStop2(2)*(-3*Cd(6780) + 
     &      MStop2(2)*(5*MTy2*Cd(2084) - 3*Cd(1735)*MStop2(2)))

        Cd(6783) = 13*MStop2(1)**4 + 
     &   MTy2*(-16*MStop2(1)**3 + 
     &      MTy2*(38*MTy2**2 + 30*MStop2(1)**2))

        Cd(6784) = 2*MStop2(1)**3 + 
     &   MTy2*(-6*MStop2(1)**2 + MTy2*(15*MTy2 + 39*MStop2(1)))

        Cd(6809) = -14*MStop2(1)**3 + 
     &   MTy2*(MTy2*(15*MTy2 - 35*MStop2(1)) + 30*MStop2(1)**2)

        Cd(6810) = -14*MStop2(1)**3 + 
     &   MTy2*(MTy2*(11*MTy2 - 47*MStop2(1)) + 62*MStop2(1)**2)

        Cd(6812) = MTy2*Cd(6809) + 
     &   MStop2(2)*(Cd(6810) + 
     &      MStop2(2)*(-Cd(6811) + Cd(2200)*MStop2(2)))

        Cd(6813) = -14*MStop2(1)**3 + 
     &   MTy2*(MTy2*(15*MTy2 - 58*MStop2(1)) - 20*MStop2(1)**2)

        Cd(6814) = -7*MStop2(1)**3 + 
     &   MTy2*(MStop2(1)**2 + MTy2*(35*MTy2 + 7*MStop2(1)))

        Cd(6817) = -42*MStop2(1)**3 + 
     &   MTy2*(MTy2*(54*MTy2 - 109*MStop2(1)) + 40*MStop2(1)**2)

        Cd(6818) = -28*MStop2(1)**3 + 
     &   MTy2*(MTy2*(55*MTy2 - 28*MStop2(1)) + 166*MStop2(1)**2)

        Cd(6819) = 21*MStop2(1)**3 + 
     &   MTy2*(-102*MStop2(1)**2 + MTy2*(24*MTy2 + 134*MStop2(1)))

        Cd(6826) = -(Cd(1272)*Cd(1280)) + 
     &   MTy2*(3*Cd(523) + MTy2*(2*MTy2 - 3*MStop2(4)))

        Cd(6830) = -5*MStop2(1)**3 + 
     &   MTy2*(MStop2(1)**2 + MTy2*(5*MTy2 + 9*MStop2(1)))

        Cd(6833) = -5*MStop2(1)**3 + 
     &   MTy2*(MTy2*(15*MTy2 - 19*MStop2(1)) + 17*MStop2(1)**2)

        Cd(6855) = -16*MStop2(1)**3 + 
     &   MTy2*(MTy2*(15*MTy2 - 40*MStop2(1)) + 33*MStop2(1)**2)

        Cd(6856) = -8*MStop2(1)**3 + 
     &   MTy2*(MTy2*(8*MTy2 - 23*MStop2(1)) + 35*MStop2(1)**2)

        Cd(6857) = MTy2*Cd(6855) + 
     &   MStop2(2)*(2*Cd(6856) + 
     &      MStop2(2)*(-(Cd(45)*Cd(2065)) + 2*Cd(2200)*MStop2(2)))

        Cd(6858) = -16*MStop2(1)**3 + 
     &   MTy2*(MTy2*(15*MTy2 - 80*MStop2(1)) - 22*MStop2(1)**2)

        Cd(6859) = -4*MStop2(1)**3 + 
     &   MTy2*(MStop2(1)**2 + MTy2*(23*MTy2 + MStop2(1)))

        Cd(6863) = -24*MStop2(1)**3 + 
     &   MTy2*(MTy2*(27*MTy2 - 70*MStop2(1)) + 22*MStop2(1)**2)

        Cd(6864) = -16*MStop2(1)**3 + 
     &   MTy2*(MTy2*(43*MTy2 - 4*MStop2(1)) + 94*MStop2(1)**2)

        Cd(6865) = 12*MStop2(1)**3 + 
     &   MTy2*(-57*MStop2(1)**2 + MTy2*(18*MTy2 + 68*MStop2(1)))

        Cd(6870) = MStop2(1)**3 + 
     &   MTy2*(2*MStop2(1)**2 + MTy2*(6*MTy2 + 2*MStop2(1)))

        Cd(6872) = 5*Cd(6870) + 
     &   MStop2(2)*(Cd(6871) + Cd(6507)*MStop2(2) + MStop2(2)**2)

        Cd(6873) = -5*MStop2(1)**3 + 
     &   MTy2*(-3*MStop2(1)**2 + MTy2*(8*MTy2 + 24*MStop2(1)))

        Cd(6878) = -5*MStop2(1)**4 + 
     &   MTy2*(12*MStop2(1)**3 + 
     &      MTy2*(24*MTy2**2 + 6*MStop2(1)**2))

        Cd(6900) = MStop2(1)**3 + 
     &   MTy2*(-2*MStop2(1)**2 + MTy2*(6*MTy2 + 7*MStop2(1)))

        Cd(6903) = -2*MStop2(1)**3 + 
     &   MTy2*(4*MStop2(1)**2 + MTy2*(7*MTy2 + 10*MStop2(1)))

        Cd(6914) = 3*Cd(831)*MStop2(1)*MStop2(2) + 
     &   MTy2*(MTy2*(-25*MTy2**2 + 2*Cd(6911)) - 
     &      Cd(6913)*MStop2(4))

        Cd(6930) = -3*MStop2(1)**3 + 
     &   MTy2*(MTy2*(4*MTy2 - 11*MStop2(1)) + 11*MStop2(1)**2)

        Cd(6935) = MStop2(1)**3 + 
     &   MTy2*(MTy2*(4*MTy2 - 42*MStop2(1)) + 12*MStop2(1)**2)

        Cd(6938) = -19*MStop2(1)**4 + 
     &   MTy2*(56*MStop2(1)**3 + 
     &      MTy2*(20*MTy2**2 - 74*MStop2(1)**2))

        Cd(6939) = -8*MStop2(1)**3 + 
     &   MTy2*(MTy2*(3*MTy2 - 10*MStop2(1)) - 16*MStop2(1)**2)

        Cd(6943) = -5*MStop2(1)**3 + 
     &   MTy2*(MTy2*(6*MTy2 - 12*MStop2(1)) + 15*MStop2(1)**2)

        Cd(6945) = 13*MStop2(1)**3 + 
     &   MTy2*(MTy2*(30*MTy2 - 14*MStop2(1)) + 9*MStop2(1)**2)

        Cd(6946) = -15*MStop2(1)**3 + 
     &   MTy2*(49*MStop2(1)**2 + MTy2*(10*MTy2 + 32*MStop2(1)))

        Cd(6949) = -5*MStop2(1)**4 + 
     &   MTy2*(68*MStop2(1)**3 + 
     &      MTy2*(30*MTy2**2 - 98*MStop2(1)**2))

        Cd(6950) = 10*MStop2(1)**3 + 
     &   MTy2*(MTy2*(10*MTy2 - 17*MStop2(1)) - 8*MStop2(1)**2)

        Cd(6952) = -10*MGlmT2*MTy2 + 
     &   MGl2*(5*MGlmT2 - 16*MStop2(1)) + 23*MStop2(1)**2

        Cd(6963) = Cd(6960)*MStop2(1) + 
     &   MStop2(2)*(Cd(6961) + MStop2(2)*(-Cd(6962) + 5*MStop2(2)))

	Cd(686) = Cd(10191) - 4*At*MTy*UCStopC(1,3)

        Cd(6987) = 2*(Cd(6985) + 
     &     MStop2(2)*(-3*Cd(6986) + 
     &        4*(3*MTy2 - MStop2(2))*MStop2(2)))

        Cd(6997) = MStop2(1)**3 + 
     &   MTy2*(-2*MStop2(1)**2 + MTy2*(2*MTy2 + MStop2(1)))

        Cd(6999) = 2*(Cd(6997) - 
     &     MStop2(2)*(Cd(6998) - 6*MTy2*MStop2(2) + 2*MStop2(2)**2)
     &     )

        Cd(7030) = -2*MStop2(1)**3 + 
     &   MTy2*(5*MStop2(1)**2 + MTy2*(MTy2 + 4*MStop2(1)))

        Cd(7038) = 4*MTy2*(MTy2**2 + Cd(5456)) - 
     &   (10*MTy2**2 + Cd(1227))*MStop2(4)

        Cd(7043) = MGl2*(-2*Cd(5454) + 
     &      MGl2*(3*MGl2 + 2*Cd(5731))) + Cd(7038)

	Cd(7048) = Cd(5015) + 2*MGlpTmSt2(2)*MStop2(3)

        Cd(7060) = 9*Cd(256)*Cd(286) + 
     &   MTy2*(27*MTy2 + 52*MStop2(1)) + 2*Cd(7059)*MStop2(2)

        Cd(7065) = -6*MStop2(1)**3 + 
     &   MTy2*(7*MStop2(1)**2 + MTy2*(3*MTy2 + 2*MStop2(1)))

        Cd(7069) = 24*MStop2(1)**3 + 
     &   MTy2*(-45*MStop2(1)**2 + MTy2*(MTy2 + 35*MStop2(1)))

	Cd(7079) = -Cd(2065) + Cd(10678)

	Cd(7098) = 4*MGl2 - 16*MTy2 - 16*MStop2(2) + 5*MStop2(3)

	Cd(7107) = 7*MGl2*(3*MGl2 + 2*Cd(286)) + Cd(7106)

        Cd(7148) = -7*Cd(523) + 2*MTy2*(MTy2 - MStop2(1)) - 
     &   2*Cd(7147)*MStop2(2)

        Cd(7149) = 2*MTy2**3 + 3*MTy2*Cd(523) + 
     &   (MTy2**2 - Cd(831))*MStop2(4)

        Cd(7150) = MStop2(1)**3 + 
     &   MTy2*(-12*MStop2(1)**2 + MTy2*(4*MTy2 + 26*MStop2(1)))

        Cd(7151) = 2*MStop2(1)**3 + 
     &   MTy2*(MTy2*(MTy2 - 8*MStop2(1)) - 4*MStop2(1)**2)

        Cd(7155) = Cd(832)*MStop2(4) + 
     &   MTy2*(-5*Cd(523) + MTy2*(2*MTy2 + MStop2(4)))

        Cd(7157) = MTy2**2*Cd(2193) + 
     &   MStop2(2)*(2*MTy2*Cd(2092) + 
     &      (Cd(2994) - 2*MStop2(2))*MStop2(2))

        Cd(7161) = MGl2*(MGl2*(MGl2 - Cd(2088)) + 2*Cd(7156)) + 
     &   2*Cd(7157)

        Cd(7176) = -5*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - 8*MStop2(1)) + 3*MStop2(1)**2)

        Cd(7180) = -6*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - 17*MStop2(1)) + 60*MStop2(1)**2)

        Cd(7188) = 4*MStop2(1)**3 + 
     &   MTy2*(-18*MStop2(1)**2 + MTy2*(3*MTy2 + 31*MStop2(1)))

	Cd(7191) = Cd(5550) + 12*MTy2*(5*MTy2 + 3*MStop2(4))

	Cd(7192) = 2*(MGl2 - 5*MTy2 - 7*MStop2(2)) + 3*MStop2(3)

        Cd(7207) = 15*MStop2(4)**3 + 
     &   MTy2*(Cd(7206) + MTy2*(120*MTy2 + 36*MStop2(4)))

	Cd(7237) = Cd(5921) + 10*MTy2*(MTy2 - MStop2(4))

        Cd(7240) = -(Cd(594)*Cd(5152)*MStop2(4)) + 
     &   MTy2*(Cd(7239) + MTy2*(10*MTy2 + 5*MStop2(4)))

        Cd(7244) = MStop2(1)**4 + 
     &   MTy2**2*(MTy2*(10*MTy2 - 5*MStop2(1)) - 2*MStop2(1)**2)

        Cd(7246) = -6*MStop2(1)**3 + 
     &   MTy2*(MTy2*(MTy2 - 35*MStop2(1)) + 28*MStop2(1)**2)

        Cd(7256) = 10*MStop2(1)**3 + 
     &   MTy2*(-20*MStop2(1)**2 + MTy2*(34*MTy2 + 43*MStop2(1)))

        Cd(7275) = 26*MTy2**3 + MTy2*Cd(7274) - 
     &   (27*MTy2**2 + MStop2(1)*MStop2(2))*MStop2(4)

        Cd(7282) = -9*MStop2(1)**3 + 
     &   MTy2*(MTy2*(53*MTy2 - 90*MStop2(1)) + 102*MStop2(1)**2)

        Cd(7316) = 78*MStop2(2)**2 - Cd(5392)*MStop2(3) + 
     &   MTy2*(174*MTy2 + 82*MStop2(4))

        Cd(7330) = Cd(1272)*Cd(5420) + 
     &   MTy2*(Cd(7329) + 5*MTy2*(66*MTy2 + 23*MStop2(4)))

        Cd(7349) = 2*MStop2(1)**3 + 
     &   MTy2*(MTy2*(16*MTy2 - 68*MStop2(1)) + 39*MStop2(1)**2)

        Cd(7353) = 74*MStop2(1)**3 + 
     &   MTy2*(31*MStop2(1)**2 + MTy2*(24*MTy2 + 68*MStop2(1)))

        Cd(7372) = -5*MStop2(1)**2*MStop2(2)**2 + 
     &   MTy2**2*(3*Cd(5643) + MTy2*(27*MTy2 - 32*MStop2(4)))

        Cd(7375) = 76*MTy2**3 - 4*MTy2*Cd(7374) - 
     &   (22*MTy2**2 + 5*Cd(5442))*MStop2(4)

        Cd(7376) = -19*MStop2(1)**3 + 
     &   MTy2*(MTy2*(149*MTy2 - 46*MStop2(1)) + 14*MStop2(1)**2)

        Cd(7377) = 15*MStop2(1)**3 + 
     &   MTy2*(-5*MStop2(1)**2 + MTy2*(46*MTy2 + 46*MStop2(1)))

        Cd(7380) = MTy2*Cd(7376) - 
     &   MStop2(2)*(Cd(7377) + 
     &      MStop2(2)*(-Cd(7378) + Cd(7379)*MStop2(2)))

        Cd(7381) = -23*MStop2(1)**3 + 
     &   MTy2*(MTy2*(106*MTy2 - 107*MStop2(1)) + 76*MStop2(1)**2)

        Cd(7382) = 14*MStop2(1)**3 + 
     &   MTy2*(-49*MStop2(1)**2 + MTy2*(107*MTy2 + 84*MStop2(1)))

        Cd(7383) = -15*MStop2(1)**3 + 
     &   MTy2*(46*MStop2(1)**2 + MTy2*(76*MTy2 + 49*MStop2(1)))

        Cd(7385) = MTy2**2*Cd(7381) - MTy2*Cd(7382)*MStop2(2) + 
     &   MStop2(2)**2*(Cd(7383) - Cd(7384)*MStop2(2))

        Cd(7396) = -10*MStop2(1)**3 + 
     &   MTy2*(-15*MStop2(1)**2 + MTy2*(16*MTy2 + 20*MStop2(1)))

        Cd(7401) = -106*MStop2(1)**3 + 
     &   MTy2*(33*MStop2(1)**2 + MTy2*(24*MTy2 + 60*MStop2(1)))

        Cd(7407) = -40*MStop2(1)**3 + 
     &   MTy2*(120*MStop2(1)**2 + MTy2*(47*MTy2 + 46*MStop2(1)))

	Cd(7423) = Cd(7211) + 8*MTy2*(MTy2 - MStop2(4))

        Cd(7424) = 40*MStop2(2)**2 - Cd(6681)*MStop2(3) + 
     &   MTy2*(8*MTy2 + 8*MStop2(4))

        Cd(7428) = 3*Cd(5442)*MStop2(4)**2 + 
     &   MTy2*(MTy2*(24*MTy2**2 + 6*Cd(7211)) - 
     &      2*Cd(7359)*MStop2(4))

        Cd(7440) = 9*Cd(535)*Cr(1) + 
     &   MTy2*(27*MTy2 + 2*MStop2(1)) + 4*Cd(3992)*MStop2(2)

        Cd(7448) = -2*MStop2(1)**3 + 
     &   MTy2*(MTy2*(2*MTy2 - 2*MStop2(1)) + MStop2(1)**2)

        Cd(7454) = -3*MStop2(1)**3 + 
     &   MTy2*(45*MStop2(1)**2 + MTy2*(4*MTy2 + 5*MStop2(1)))

        Cd(7461) = 3*MStop2(1)**3 + 
     &   MTy2*(-12*MStop2(1)**2 + MTy2*(3*MTy2 + 26*MStop2(1)))

        Cd(7465) = 60*MStop2(2)**2 - Cd(7464)*MStop2(3) + 
     &   MTy2*(60*MTy2 + 36*MStop2(4))

        Cd(7466) = 40*MTy2**3 + 2*MTy2*Cd(4919) + 
     &   (12*MTy2**2 + Cd(4921))*MStop2(4)

	Cd(7468) = MGl2 - 5*MTy2 - Cd(1060)

        Cd(7481) = -42*MStop2(2)**2 + 3*Cd(5596)*MStop2(3) + 
     &   MTy2*(34*MTy2 + 32*MStop2(4))

        Cd(7484) = 32*MTy2**3 + MTy2*Cd(7483) - 
     &   (33*MTy2**2 + MStop2(1)*MStop2(2))*MStop2(4)

        Cd(7487) = 88*MTy2**3 + MTy2*Cd(7486) + 
     &   (69*MTy2**2 - 5*Cd(5601))*MStop2(4)

        Cd(7488) = -8*MStop2(1)**3 + 
     &   MTy2*(MTy2*(20*MTy2 - 108*MStop2(1)) + 111*MStop2(1)**2)

        Cd(7494) = -9*MStop2(1)**3 + 
     &   MTy2*(MTy2*(11*MTy2 - 159*MStop2(1)) + 129*MStop2(1)**2)

        Cd(7500) = 4*MStop2(1)**3 + 
     &   MTy2*(-47*MStop2(1)**2 + MTy2*(13*MTy2 + 94*MStop2(1)))

        Cd(7518) = -Cd(7375) - 
     &   MGl2*(Cd(7503) + MGl2*(5*MGl2 - Cd(7510)))

	end
