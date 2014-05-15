* TLhr_abbr1_4.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_4
	implicit none

#include "TLhr.h"


        Cd(1012) = -4*Cd(169)*MStop2(1)**2 + 
     &   2*Cd(1011)*MStop2(2) + 4*Cd(1010)*UStop2(1,2)**2

        Cd(10210) = -(Cd(10209)*UCStop(1,3)**2) + 
     &   8*M_3*UStop2(1,1)*
     &    (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

        Cd(10211) = M_3*UCStopC(1,3)**2 + 
     &   M_3C*UStop2(1,1)*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

	Cd(10222) = MTy2 - 2*M_3C*MTy*UCStop(1,3) + MGl2*UStop2(1,2)

        Cd(10227) = MTy2*
     &    (2*M_3C*Cd(795) + 2*M_3*Cd(10219)) - 
     &   MTy*Cd(13)*MGlpTmSt2(1)

	Cd(10239) = Cd(10015) + 2*MStop2(2)

	Cd(10242) = 2*MTy2 - MTy*(AtC*UCStop(1,3) + At*UCStopC(1,3))

        Cd(10248) = 2*MStop2(2) + 
     &   At*(M_3*UCStopC(1,3)**2 + 
     &      M_3C*UStop2(1,1)*UStop2(1,2))

        Cd(10250) = AtC*(AtC*Cd(10032) + 2*Cd(10249)) + 
     &   At*Cd(10248)*UCStopC(1,3)

        Cd(10252) = Cd(10251)*UCStopC(1,3) + 
     &   2*M_3C*UStop2(1,1)*
     &    (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

	Cd(10264) = MTy2 + 2*M_3C*MTy*UCStop(1,3) + MGl2*UStop2(1,2)

        Cd(10294) = AtC**2*Cd(10077)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(9969)*Cr(28) + At*Cd(10077)*UCStopC(1,3)**2)

        Cd(10299) = 2*Cd(10298)*UStop2(1,1) - UStop2(1,1)**2 + 
     &   Cd(172)*UStop2(1,2)

        Cd(10306) = (Cd(6306) - 2*UStop2(1,1))*UStop2(1,1) + 
     &   Cd(1052)*UStop2(1,2)

	Cd(10311) = -Cd(10141) + AtC*UCStop(1,3)

        Cd(10313) = -2*M_3C*UCStop(1,3) + AtC*UCStop(1,3) + 
     &   Cd(10312)*UCStopC(1,3)

        Cd(10322) = At*Cd(10321) + 
     &   UStop2(1,1)*(Cd(9948)*Cd(10320) - 2*MGlmT2*UStop2(1,2))

        Cd(10324) = At*Cr(5)*Cr(28) + 
     &   UStop2(1,1)*(Cd(9948)*Cd(10323) - 
     &      MGlmT2*UCStop(1,3)*UStop2(1,2))

        Cd(10325) = AtC*(-2*Cd(10324) + 
     &      AtC*Cd(10088)*UCStop(1,3)) - At*Cd(10322)*UCStopC(1,3)

	Cd(10341) = (-2 + UStop2(1,2))*(-5 + 4*UStop2(1,2))

        Cd(10357) = 8*M_3*Cr(85)*UCStopC(1,3)**2 - 
     &   4*M_3C*Cr(28)*UStop2(1,1) + AtC*Cr(28)*UStop2(1,1)

	Cd(10367) = Cd(13737) - UStop2(1,1)**2

	Cd(10384) = Cd(152) + MGlmT2*UStop2(1,2)

        Cd(10393) = 8*M_3C*Cd(795) + 8*M_3*Cd(10219) - 
     &   2*AtC*Cd(10382)

	Cd(10401) = 2 - Cd(9969)

	Cd(10413) = Cd(1052) - 2*UStop2(1,1)

	Cd(10418) = MGlpT2*UStop2(1,1) + MStop2(1)*UStop2(1,2)

	Cd(10420) = MGlpT2*UStop2(1,1) + Cd(10419)*UStop2(1,2)

        Cd(10421) = At*Cd(10329) + Cd(9948)*Cd(10420) - 
     &   Cd(10418)*UStop2(1,2)

        Cd(1034) = -(Cd(1033)*UStop2(1,2)**2) + 
     &   MStop2(1)*(5*MStop2(3) - 14*MStop2(1)*UStop2(1,2))

	Cd(10435) = Cd(427) + 2*MUE2*Cd(9981)*UStop2(1,1)

        Cd(10438) = Cr(28)*MStop2(1)**2 + 
     &   UStop2(1,1)*(-(Cd(10437)*MStop2(1)) + 
     &      2*MGlmT2*MStop2(3)*UStop2(1,2)**2)

        Cd(10440) = MUE2*Cd(1281)*UStop2(1,1)**2 + 
     &   MStop2(1)*(Cd(427) + MUE2*UStop2(1,2)**2)

	Cd(10448) = Cd(427) + 2*MUE2*UStop2(1,1)*UStop2(1,2)**2

	Cd(10450) = Cd(427) + MUE2*Cr(28)

        Cd(10466) = Cd(9948)*(M_3 + At*UStop2(1,1)) + 
     &   At*UStop2(1,2)**2

        Cd(10474) = 13*Cd(47) + 4*Cd(152)*Cd(164) + 
     &   At*AtC*Cd(164)*UStop2(1,1)

	Cd(10479) = Cd(10356) + Cd(9948)*Cd(10359)

        Cd(10480) = Cd(9948)*Cd(10360) - 
     &   2*M_3*Cr(85)*UCStopC(1,3)**2 + M_3C*Cr(28)*UStop2(1,1)

        Cd(10485) = Cd(47)*
     &    (-(M_3*Cd(9948)*Cd(10453)) + Cd(10456)) + 
     &   MTy2*Cd(164)*Cd(10086)*MStop2(3)

	Cd(10488) = Cd(10200) + 4*MStop2(1)*UCStop(1,3)

        Cd(10489) = 4*M_3*Cd(45) + 
     &   At*(-2*MGlpTmSt2(1) + M_3*At*UCStopC(1,3)**2)

	Cd(10493) = (MGl2 - MStop2(1))*UStop2(1,1) - MTy2*UStop2(1,2)

        Cd(10501) = 2*Cd(9948)*Cd(10500) - Cd(10306)*MStop2(3) - 
     &   2*Cd(10500)*UStop2(1,2)

        Cd(10508) = 2*Cd(9948)*Cd(10500) - Cd(10306)*MStop2(3) - 
     &   2*Cd(10500)*UStop2(1,2)

        Cd(10510) = -(MTy2*UCStop(1,3)) + 
     &   At*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))

        Cd(10511) = (MGl2 - MStop2(1))*UCStop(1,3)*UStop2(1,1) + 
     &   Cd(10510)*UStop2(1,2)

        Cd(10512) = Cd(9948)*Cd(10511) + 
     &   At*Cr(5)*Cr(75)*UStop2(1,1) - 
     &   Cd(10493)*UCStop(1,3)*UStop2(1,2)

	Cd(1044) = Cd(971) + 2*Cd(1043)

        Cd(10532) = -4*M_3*Cd(10219) + At*Cd(10219) + 
     &   Cd(795)*Cd(10398)

	Cd(10533) = MGlmT2*Cd(9948) - Cd(10383)

        Cd(10539) = MUE*Cd(118)*UCStop(1,3)**2 + 
     &   MUEC*UStop2(1,1)*(MGl2*Cd(9948)*UStop2(1,1) - 
     &      Cd(10538)*UStop2(1,2))

        Cd(10541) = MUEC*Cd(118)*UCStopC(1,3)**2 + 
     &   MUE*UStop2(1,1)*(MGl2*Cd(9948)*UStop2(1,1) - 
     &      Cd(10538)*UStop2(1,2))

        Cd(10548) = -(Cd(6320)*UStop2(1,1)) + 2*UStop2(1,1)**2 - 
     &   Cd(290)*UStop2(1,2)

        Cd(10552) = M_3*Cd(9948)*Cd(10303) - 
     &   MGlmT2*Cd(10454) + Cd(10455)*MStop2(1)

        Cd(10553) = Cd(9948)*Cd(10303) - MGlmT2*Cr(28) + 
     &   Cr(75)*MStop2(1)*UStop2(1,1)

	Cd(10555) = Cd(10033) + MStop2(1)*UCStop(1,3)

        Cd(10556) = MStop2(1) + 
     &   At*(M_3*UCStopC(1,3)**2 + 
     &      M_3C*UStop2(1,1)*UStop2(1,2))

        Cd(10585) = MUE**2*UCStop(1,3)**3 + 
     &   2*MUE2*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) + 
     &   MUEC**2*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Cd(10586) = Cd(12609) + 
     &   2*MUE2*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Cd(10588) = 2*MUE2*Cr(28)*UCStop(1,3) - 
     &   3*MUE**2*UCStop(1,3)**3 + MUEC**2*Cr(28)*UCStopC(1,3)

        Cd(10592) = -3*MUE**2*UCStop(1,3)**3 + 
     &   Cr(75)*(2*MUE2*UCStop(1,3) + MUEC**2*UCStopC(1,3))*
     &    UStop2(1,1)

	Cd(10605) = Cd(12666) - UStop2(1,2)

        Cd(10608) = -(AtC**2*UCStop(1,3)**2) + 
     &   At*(AtC*Cd(10070) + At*Cd(10077)*UCStopC(1,3)**2)

        Cd(1053) = -4*Cd(1052)*MStop2(1)**2 + 
     &   2*Cd(1011)*MStop2(2) - 4*Cd(1051)*UStop2(1,2)**2

	Cd(1055) = 5 + UStop2(1,2)*(6 + 2*UStop2(1,2)**2)

	Cd(10644) = Cd(10349) + M_3*Cd(9948)*UStop2(1,1)**2

	Cd(10647) = 2 - 3*Cd(9948) + 2*UStop2(1,2)

        Cd(10697) = 2*MUEC*UCStopC(1,3)**2 + 
     &   MUE*(3*UStop2(1,1)**2 + UStop2(1,2)**2)

        Cd(10709) = MGl2*MTy*Cd(13) + 
     &   2*MTy2*(M_3C*Cd(795) + M_3*Cd(10219))

        Cd(10710) = MTy2*
     &    (2*M_3C*Cd(795) + 2*M_3*Cd(10219)) + 
     &   MTy*Cd(13)*MGlpTmSt2(2)

        Cd(10712) = UStop2(1,1)**3 + 
     &   (Cd(172) - UStop2(1,1))*UStop2(1,2)**2

        Cd(10715) = -(Cd(788)*UStop2(1,1)) + 2*UStop2(1,1)**2 - 
     &   3*UStop2(1,2)

        Cd(10722) = -(Cd(10721)*UStop2(1,1)) + 
     &   2*UStop2(1,1)**2 - UStop2(1,2)

        Cd(10743) = Cd(9948)*Cd(10741) - Cd(10742)*MStop2(1) + 
     &   MStop2(2)*MStop2(3)*UStop2(1,2)

        Cd(10756) = -(MUE**2*UCStop(1,3)**2) - 
     &   MUEC**2*Cr(85)*UCStopC(1,3)**2 + 
     &   MUE2*Cr(75)**2*UStop2(1,1)

        Cd(10780) = 2*M_3C*UCStop(1,3) + AtC*UCStop(1,3) + 
     &   Cd(10001)*UCStopC(1,3)

        Cd(10808) = -(At*Cd(10806)) - 
     &   2*AtC*Cd(10807)*UCStop(1,3)**2 + 
     &   2*M_3*MStop2(2)*MStop2(3)*UStop2(1,1)

	Cd(10811) = AtC*MStop2(1) - MStop2(3)*(M_3C + AtC*UStop2(1,2))

        Cd(10815) = -(AtC*Cd(10806)) - 
     &   2*At*Cd(10741)*UCStopC(1,3)**2 + 
     &   2*M_3C*MStop2(2)*MStop2(3)*UStop2(1,1)

	Cd(10826) = MGlmT2*Cd(9948) - Cd(10742)

        Cd(10834) = M_3*MUE*UCStop(1,3) + 
     &   MUEC*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))

        Cd(10840) = M_3C*MUEC*UCStopC(1,3) + 
     &   MUE*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))

	Cd(10865) = 4*M_3 + At*(1 + UStop2(1,2))

        Cd(10880) = -3*MUEC**2*UCStopC(1,3)**3 + 
     &   MUE**2*UCStop(1,3)*UStop2(1,1)**2 + 
     &   2*MUE2*UCStopC(1,3)*UStop2(1,1)**2

        Cd(10910) = -(MUE2*UUStop(1,2)*UUStop(2,1)*
     &      UUStopC(2,2)) + 
     &   MUEC*(MUEC*UUStop(3,1)*UUStop(3,2)*UUStopC(2,1) + 
     &      Cd(10909)*UUStopC(2,3))

	Cd(10911) = Cd(11672) - 2*MUE2*UUStop(3,1)*UUStop(3,2)

        Cd(10912) = MUEC*Cd(10906)*UUStopC(2,1) + 
     &   MUE*Cd(10909)*UUStopC(2,2) - Cd(10911)*UUStopC(2,3)

        Cd(10913) = Cd(10910)*UUStopC(1,1) + 
     &   Cd(10908)*UUStopC(1,2) - Cd(10912)*UUStopC(1,3)

	Cd(10934) = -Cd(10938) + 5*MStop2(3)**2*UUStop(1,3)

	Cd(10943) = Cd(11263) - Cd(10929)*MStop2(3)

        Cd(10965) = Cd(10959) + 5*MStop2(3)**2 - 
     &   8*MStop2(1)**2*UStop2(1,1)**2

	Cd(10970) = Cd(11266) - Cd(10950)*MStop2(3)

        Cd(10995) = UUStop(1,2)*UUStop(2,1) - 
     &   UUStop(1,1)*UUStop(2,2) + 2*UUStop(3,1)*UUStop(3,2)

	Cd(11000) = Cd(11132) - AtC*UUStop(3,1)

        Cd(11001) = UUStop(1,2)*UUStop(2,1) - 
     &   UUStop(1,1)*UUStop(2,2) - 2*UUStop(3,1)*UUStop(3,2)

        Cd(11058) = 2*Cd(9794)**2 + AtC*Cd(11057)*UCStop(1,3) + 
     &   Cd(11056)*UStop2(1,2)**2

        Cd(11063) = At*UUStop(3,2) + 
     &   MTy*(-6*UUStop(3,3) - 4*UUStop(3,4))

        Cd(11068) = 8*MTy2*Cr(58)*UStop(1,1) + 
     &   At*(4*MTy*Cr(58)*UStop(1,2) + At*UStop(2,2)*UUStop(1,2))

        Cd(11113) = 2*MTy*Cd(9809)*UCStop(3,1) + 
     &   MStop2(3)*UStop(1,1) - Cd(9794)*Cd(9809)*UStopC(2,1)

        Cd(11114) = 8*MTy2*UUStop(3,1) + 
     &   At*(2*MTy*Cr(91) - At*UUStop(3,2))

        Cd(11121) = At*(MTy*Cd(10995) + 
     &      At*UUStop(1,2)*UUStop(2,3)) - 2*MTy2*Cr(58)*UUStop(3,1)

	Cd(11151) = Cd(11256) - Cd(11110)*MStop2(3)

        Cd(11163) = MUEC*(2*Cd(11162) + Cd(11160)*UCStop(3,1)) - 
     &   MUE*Cd(9794)*UCStop(1,3)*UStop(1,1)

        Cd(11186) = MTy*Cd(11106)*UCStop(3,1) + 
     &   At*Cd(10925)*UCStop(3,4) + 
     &   Cd(10980)*UStop(1,2)*UStop2(1,1)

        Cd(11200) = -(Cd(11198)*MStop2(2)) + 
     &   MStop2(1)*(-(At*Cd(11199)) + 8*MTy*UStop(1,1)*UStop2(1,2))

        Cd(11205) = At*(MTy*Cd(11001) - 
     &      At*UUStop(1,3)*UUStop(2,2)) - 2*MTy2*Cr(58)*UUStop(3,1)

        Cd(11206) = Cd(11204)*UCStop(3,1)*UStopC(2,2) + 
     &   Cd(698)**2*UUStopC(2,1) + 4*Cd(11205)*UUStopC(2,3)

        Cd(11221) = 2*MUE2*Cr(91)*UUStop(1,3) - 
     &   MUE**2*UUStop(1,1)*UUStop(3,1) - 
     &   MUEC**2*UUStop(1,2)*UUStop(3,2)

	Cd(109) = -MGl2**2 + Cd(44)**2 - MGlpTmSt4(1)*UStop2(1,2)**2

        Cd(11278) = -(Cd(11277)*UUStopC(2,1)) + 
     &   MUE*(2*Cd(11275)*UCStop(3,1)*UStopC(2,1) - 
     &      MUE*Cr(91)**2*UUStopC(2,2))

        Cd(11287) = 2*MUE2*Cd(11286)*UStop2(1,2) + 
     &   MUE**2*Cd(970)*UUStop(3,1) - 
     &   MUEC**2*UCStopC(1,3)**2*UUStop(3,2)

	Cd(11288) = Cd(11672) - 6*MUE2*UUStop(3,1)*UUStop(3,2)

        Cd(11301) = 2*MUE2*Cr(91)*UUStop(2,3) - 
     &   MUE**2*UUStop(2,1)*UUStop(3,1) - 
     &   MUEC**2*UUStop(2,2)*UUStop(3,2)

        Cd(11307) = MUEC*Cd(11306)*UCStopC(1,3) + 
     &   MUE*(Cd(10925)*UStop2(1,2)**2 - 3*UUStop(3,3))

        Cd(11322) = 2*MTy*Cd(11321)*MStop2(2)*UCStopC(1,3) + 
     &   Cd(11319)*UStop(2,2) - 
     &   Cd(11204)*MStop2(2)**2*UCStop(3,1)*UUStopC(1,1)

        Cd(11326) = 8*Cd(11325)*MStop2(2)*UCStopC(3,4) - 
     &   Cd(11323)*UCStop(1,3)*UStop(2,1) - 
     &   4*Cr(91)**2*MStop2(2)**2*UStopC(2,2)*UUStopC(1,2)

	Cd(110) = Cd(164) + MGlpTmSt2(1)*UStop2(1,2)

        Cd(11329) = Cd(11327)*UStop(2,1) + 
     &   At*(MStop2(2)*UCStopC(3,1)*UStop(1,2) + 
     &      Cd(11328)*UStop(2,2))

        Cd(1126) = 2*MStop2(1) - Cd(1125)*UStop2(1,1) + 
     &   Cd(1121)*UStop2(1,1)**3

        Cd(11342) = Cd(11316)*MStop2(1)**2 + 
     &   2*MStop2(2)*(-(Cd(11341)*MStop2(1)) + MStop2(2))

        Cd(11351) = 2*MTy*Cd(11343)*UUStop(3,1) + 
     &   At*(3*UUStop(3,3) - UStop2(1,1)**2*UUStop(3,4))

        Cd(1128) = 4*Cd(652) - Cd(1122)*Cr(28) + 
     &   2*Cd(1121)*UStop2(1,1)*UStop2(1,2)

        Cd(11366) = Cd(11017)*Cd(11106)*UStop2(1,1)*
     &    UStopC(2,2) + 6*MTy*UUStop(3,1) + 3*At*UUStop(3,4)

	Cd(11375) = Cd(970) + Cr(7)*UStop2(1,1)

        Cd(11377) = MUE*Cd(11376) - AtC*Cd(11175)*UUStop(1,1) - 
     &   MUEC*Cd(11187)*UUStop(3,2)

	Cd(11383) = Cd(11025) + MUE*Cr(91)**2*UUStopC(2,2)

        Cd(11395) = 6*MTy*Cd(11394) + 3*Cd(424)*UUStop(3,3) + 
     &   Cd(11393)*UUStopC(2,1)

	Cd(11427) = Cd(11283) + 2*Cd(11426)

        Cd(11429) = 2*MUE2*Cd(11286)*UStop2(1,2) + 
     &   MUE**2*Cd(2364)*UUStop(3,1) - 
     &   MUEC**2*UCStopC(1,3)**2*UUStop(3,2)

        Cd(11445) = Cd(11444) + MUE**2*Cd(970)*UUStop(2,1) + 
     &   Cd(11297)*UUStop(2,2)

	Cd(1137) = Cd(1042) + 2*Cd(1136)

        Cd(11449) = MUEC*Cd(11306)*UCStopC(1,3) + 
     &   MUE*(Cd(10925)*UStop2(1,2)**2 - UUStop(3,3))

	Cd(11460) = Cd(11317) + 2*Cd(11459)

	Cd(11464) = Cd(11320) + 2*Cd(11463)

        Cd(1139) = Cd(1134) + 
     &   MTy2*(2*Cd(1135)*UStop2(1,1) - 5*UStop2(1,1)**3)

	Cd(11471) = Cd(11324) + 2*Cd(11470)

	Cd(11484) = Cd(11483) - 3*UStop2(1,1)*UStop2(1,2)

        Cd(11485) = Cd(11483)*MStop2(1)**2 + 
     &   MStop2(2)*(-2*Cd(11484)*MStop2(1) + Cd(11482)*MStop2(2))

        Cd(11492) = -(Cd(11161)*UStop2(1,1)**2) + UUStop(3,3) + 
     &   UUStop(3,4)

	Cd(11493) = Cd(11348) + 2*Cd(11492)

        Cd(11496) = 2*MTy*Cr(53)*UUStop(3,1) + 
     &   At*(UUStop(3,3) - UStop2(1,1)**2*UUStop(3,4))

	Cd(11504) = Cd(11355) + 2*Cd(11503)

	Cd(11507) = Cd(11357) + 2*Cd(11506)

	Cd(1144) = 2 + UStop2(1,2)*(10 - UStop2(1,2)**2)

        Cd(11518) = Cd(11017)*Cd(11106)*UStop2(1,1)*
     &    UStopC(2,2) + 2*MTy*UUStop(3,1) + At*UUStop(3,4)

        Cd(11529) = 2*MTy*Cd(11394) + Cd(424)*UUStop(3,3) + 
     &   Cd(11393)*UUStopC(2,1)

	Cd(11541) = -Cd(11412) + 2*Cd(11540)

	Cd(11555) = 6 + 3*Cd(1136) - 2*UStop2(1,2)**2

        Cd(11576) = MStop2(2) + 
     &   MTy2*(2*Cd(11575)*UStop2(1,1) - 3*UStop2(1,1)**3)

	Cd(11579) = 6 + 3*Cd(11490) - 2*UStop2(1,1)**2

        Cd(1151) = 2*Cd(1149)*UStop2(1,1) - 
     &   Cd(1148)*UStop2(1,1)**2 - Cd(1141)*UStop2(1,2)**2

        Cd(11593) = 2*Cd(11592) + 16*MTy2*UStop(2,1) + 
     &   At*MTy*Cd(11458)*UStop(2,2)

	Cd(11598) = Cd(11365) + 3*Cd(11516)

        Cd(11627) = UStop2(1,1)*
     &    (3*UStop2(1,1) - 2*UStop2(1,2)) + 3*UStop2(1,2)**2

        Cd(11640) = UStop2(1,1)*
     &    (6*UStop2(1,1) - 4*UStop2(1,2)) + 6*UStop2(1,2)**2

        Cd(11642) = UStop2(1,1)*
     &    (2*UStop2(1,1) - 3*UStop2(1,2)) + 2*UStop2(1,2)**2

        Cd(11692) = Cd(11691)*UStop2(1,2)**2 + 
     &   MStop2(2)*(2*Cd(9948)*UStop2(1,1)**2 + UStop2(1,2)**3)

	Cd(11694) = Cd(11633) - 4*UStop2(1,1)*UStop2(1,2)

        Cd(11698) = UStop2(1,2)**2 + 
     &   UStop2(1,1)*(UStop2(1,1) + 4*UStop2(1,2))

        Cd(11699) = -2*Cd(11698)*MStop2(3) + 
     &   UStop2(1,1)*(Cd(9948)*MStop2(2) - Cd(674)*UStop2(1,2))

        Cd(11723) = -(Cd(2391)*UStop2(1,1)) + 
     &   4*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Cd(11726) = AtC*Cd(11725) - 5*Cd(492)*MStop2(3) - 
     &   At**2*UCStopC(1,3)**2

	Cd(11750) = AtC*Cd(11748) + At*(-Cd(11749) + At*UCStopC(1,3)**2)

        Cd(11758) = MGl2*Cd(152)*Cd(10086)*MStop2(3) + 
     &   Cd(47)*Cd(164)*(-Cd(10454) + M_3*Cd(9948)*UStop2(1,1))

        Cd(11761) = Cd(11760)*MStop2(1) + 
     &   Cd(9955)*(Cd(47) + MStop2(1)**2)*UStop2(1,1)

        Cd(11762) = 2*At*Cd(11761) + 
     &   M_3*MStop2(1)*(2*Cd(11759) + 3*At**2*UCStopC(1,3)**2)

        Cd(11772) = At*Cd(10091)*Cr(5) + 
     &   UCStop(1,3)*(MStop2(3) + MGlpT2*Cd(9955)*UStop2(1,2))

        Cd(11773) = Cd(11772)*MStop2(1) + 
     &   Cd(9955)*(Cd(47) + MStop2(1)**2)*UCStop(1,3)*UStop2(1,1)

        Cd(11779) = 4*Cd(118)*Cd(10169) + 
     &   At*MGlmT2*Cd(9969)*Cr(28) + 
     &   2*AtC*MGlmT2*Cd(10077)*UCStop(1,3)**2

        Cd(11782) = -(MGl2*Cd(118)*Cd(152)**2) + 
     &   Cd(47)*(Cd(47) + Cd(118)*MStop2(1))

        Cd(11787) = 4*Cd(118)*Cd(10176) + 
     &   AtC*MGlmT2*Cd(9969)*Cr(28) + 
     &   2*At*MGlmT2*Cd(10077)*UCStopC(1,3)**2

	Cd(11790) = MStop2(1) + MGlmT2*(Cd(9948) + UStop2(1,1))

	Cd(11792) = Cd(164) + MGlmT2*UStop2(1,2)

	Cd(11800) = 4*M_3*MUE + 4*M_3C*MUEC - AtC*MUEC

	Cd(11809) = M_3C*At*MGlmT2 - MGlmT2*MTy2 + Cd(47)

        Cd(11844) = At*Cd(10070)*MStop2(1) + 
     &   Cd(138)*(4*M_3*Cr(28) - 8*M_3C*UCStop(1,3)**2)

	Cd(1177) = Cd(1095) + Cd(1176)*MStop2(3)

        Cd(11852) = AtC + 2*M_3*UCStopC(1,3)**2 + 
     &   2*M_3C*UStop2(1,1)*UStop2(1,2)

	Cd(11858) = Cd(2391) - 4*UStop2(1,1)

	Cd(11866) = Cd(10147) + AtC*UCStop(1,3)

        Cd(11869) = -(AtC**2*Cd(10087)*UCStop(1,3)) + 
     &   At*(2*AtC*Cd(10092)*Cr(5) - At*Cd(10095)*UCStopC(1,3))

	Cd(11876) = Cd(290) + 3*Cd(9948)

	Cd(11881) = 2 - Cd(9955)

        Cd(1181) = Cd(1179)*UStop2(1,1) - UStop2(1,1)**3 + 
     &   Cd(1180)*UStop2(1,2)

	Cd(11886) = 5 - 24*Cd(9948) + 12*UStop2(1,2)

        Cd(11887) = 4*M_3*Cd(9948)*Cd(9983) + 
     &   Cd(10209)*UCStop(1,3)**2 + 4*M_3*Cr(75)**2*UStop2(1,1)

        Cd(11888) = Cd(9948)*Cd(10181) - 
     &   2*M_3*Cr(85)*UCStopC(1,3)**2 + 
     &   M_3C*Cr(75)**2*UStop2(1,1)

        Cd(11896) = (MStop2(2) - 5*MStop2(3))*UStop2(1,1)**2 - 
     &   MStop2(3)*UStop2(1,2)**2

        Cd(11897) = -(At*Cd(11896)) - 
     &   2*AtC*MStop2(1)*UCStop(1,3)**2 + 
     &   2*M_3*MStop2(3)*UStop2(1,1)

        Cd(11902) = -(AtC*Cd(11896)) + 
     &   2*MStop2(3)*(At*UCStopC(1,3)**2 + M_3C*UStop2(1,1))

	Cd(11907) = Cd(10298) + 2*UStop2(1,2)

	Cd(11909) = 1 + Cd(9969)

	Cd(1184) = 8*Cd(832) - 4*Cd(842) - Cd(1122)*MStop2(3)

	Cd(116) = 2*MTy2 + Cd(118)

	Cd(11926) = MGl2*Cd(251) - MGlmT2*MStop2(1) + Cd(1952)*MStop2(2)

        Cd(11952) = -3*MUE**2*UCStop(1,3)**3 + 
     &   2*MUE2*UCStop(1,3)*UStop2(1,1)**2 + 
     &   MUEC**2*UCStopC(1,3)*UStop2(1,1)**2

	Cd(11971) = Cd(11975) - 2*UStop2(1,1)

        Cd(11976) = (Cd(11975) - 2*UStop2(1,1))*UStop2(1,1) + 
     &   Cd(290)*UStop2(1,2)

        Cd(11989) = Cd(6320)*UStop2(1,1) - 
     &   2*(UStop2(1,1)**2 + UStop2(1,2)**2)

	Cd(1192) = Cd(1111) + 5*Cd(607)*MStop2(3)

	Cd(11995) = Cd(10162) - 4*UStop2(1,1)

	Cd(1194) = -4*Cd(982) + 4*Cd(1063) - Cd(1141)*MStop2(3)

	Cd(12013) = 5 + 4*UStop2(1,2)*(-7 + 6*UStop2(1,2))

        Cd(12015) = 8*M_3C*UCStop(1,3)**2 - 
     &   2*UStop2(1,1)*UStop2(1,2)*(-4*M_3 + At*UStop2(1,2))

	Cd(12024) = 2 + 4*Cd(9948) - 8*UStop2(1,2)

	Cd(12026) = 4 - 8*Cd(9948) - 8*UStop2(1,1)

	Cd(12029) = Cd(10027) - 2*UStop2(1,1)

	Cd(12035) = Cd(6320) + 2*Cd(9948)

        Cd(12055) = -(Cd(12054)*UStop2(1,1)) + 
     &   4*UStop2(1,1)**2 - Cd(6306)*UStop2(1,2)

	Cd(12056) = Cd(788) - 6*UStop2(1,1)

	Cd(12062) = 3 + 4*Cd(9948) - 12*UStop2(1,2)

        Cd(12063) = -(Cd(12062)*UStop2(1,1)) + 
     &   4*UStop2(1,1)**2 - Cd(12061)*UStop2(1,2)

	Cd(12070) = Cd(44)**2 + MGl2*(MGl2 - 2*Cd(2200))

        Cd(12072) = Cd(44)**3 + 
     &   MGl2*(-(Cd(44)*Cd(1592)) - MGl2*Cd(2129))

        Cd(12073) = MUE*Cr(91)*UCStop(3,3) + 
     &   MUEC*(Cd(10925)*UCStop(3,4) + UStop2(1,1)*UUStop(1,2))

        Cd(12081) = -(MUE*Cd(12079)*UUStopC(2,2)) + 
     &   MUEC*(MUEC*Cd(11106)*UCStopC(1,3)*UStopC(2,1) - 
     &      Cd(12080)*UUStopC(2,3))

	Cd(1201) = Cd(10233) + MUEC**2*UCStopC(1,3)**2

        Cd(12087) = MUE**2*UCStop(1,3)*UCStop(3,3) + 
     &   3*MUEC**2*UCStop(3,4)*UCStopC(1,3) - 
     &   2*MUE2*Cr(91)*UUStopC(2,3)

        Cd(1205) = Cd(1204)*UStop2(1,1) - UStop2(1,1)**3 + 
     &   Cr(85)*UStop2(1,2)**2

	Cd(12126) = -Cd(427) + 2*MUE2*Cd(12125)*UStop2(1,1)*UStop2(1,2)

	Cd(12130) = Cd(427) + 2*MUE2*Cd(12107)*UStop2(1,1)*UStop2(1,2)

	Cd(12142) = Cd(44)**2 + MGl2*(3*MGl2 - 4*Cd(45))

        Cd(12152) = MUEC*Cr(91)*UCStop(3,4) + 
     &   MUE*(2*Cd(11106)*UCStop(3,3) + 2*UStop2(1,2)*UUStop(1,1))

	Cd(12162) = MGlmT2**2 - 2*MGlpT2*MStop2(1) + MStop2(1)**2

	Cd(12180) = Cd(179) + 6*MTy2*UStop2(1,1)

	Cd(12205) = 4*MTy2 + Cd(1442)

        Cd(12228) = MGl2*Cd(1442) - MStop2(1)*MStop2(4) + 
     &   MGlpTmSt2(1)*MStop2(3)*UStop2(1,2)

	Cd(12235) = 4*MTy2 + Cd(1468)

        Cd(1219) = 2*(Cd(1218)*UStop2(1,1) + 
     &     UStop2(1,2)*(-(MStop2(1)*MStop2(3)) + 
     &        Cd(286)*MStop2(2)*UStop2(1,2)))

	Cd(12268) = MStop2(1)*MStop2(2) + MGl2*(MGl2 - MStop2(3))

        Cd(12272) = -UStop2(1,2)**2 + 
     &   UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))

	Cd(12273) = Cd(6320) - 4*UStop2(1,1)

	Cd(12274) = 4*MGl2 - MTy2*Cd(12273) - 4*MStop2(2)

	Cd(12279) = Cd(6320) + 8*UStop2(1,1)

        Cd(1221) = MStop2(2)*(36*MStop2(2) - 34*MStop2(3)) + 
     &   6*MStop2(3)**2

	Cd(1225) = Cd(523) + MStop2(3)**2*UStop2(1,2)

	Cd(12330) = Cd(10271) - 2*M_3C*MTy*UCStop(1,3)

	Cd(12349) = Cd(9872) - At*MTy*UCStopC(1,3)

        Cd(12352) = 4*M_3*MTy2 + 
     &   At*(-2*MGlpTmSt2(2) + Cd(10137)*UStop2(1,1)*UStop2(1,2))

	Cd(12361) = Cd(10271) + 2*M_3C*MTy*UCStop(1,3)

        Cd(12399) = -8*M_3*UCStopC(1,3)**2 - 
     &   8*M_3C*UStop2(1,1)*UStop2(1,2) + 
     &   2*AtC*UStop2(1,1)*UStop2(1,2)

	Cd(12400) = AtC*Cd(11753) - At*(Cd(12399) + At*UCStopC(1,3)**2)

	Cd(1233) = Cd(523) - Cd(315)*Cd(778)*UStop2(1,2)**2

        Cd(12412) = At*Cd(10091)*Cr(5) + 
     &   (Cd(9955)*Cd(10500) + MStop2(3))*UCStop(1,3)

        Cd(1236) = Cd(1217)*MStop2(4)*UStop2(1,1) + 
     &   Cd(842)*UStop2(1,1)**2 + Cd(1235)*UStop2(1,2)

        Cd(12435) = Cd(9948)*Cd(11795) + 
     &   Cd(10398)*Cr(85)*UCStop(1,3)**2 + 
     &   2*M_3*Cr(28)*UStop2(1,1)

        Cd(12441) = Cd(9948)*Cd(11803) + 
     &   Cd(10140)*Cr(85)*UCStopC(1,3)**2 + 
     &   2*M_3C*Cr(28)*UStop2(1,1)

        Cd(12455) = (MGl2 - MStop2(1))*UStop2(1,1) + 
     &   Cd(12454)*UStop2(1,2)

        Cd(12456) = At*Cd(10329) + Cd(9948)*Cd(12455) - 
     &   Cd(10493)*UStop2(1,2)

        Cd(12517) = AtC**2*Cd(10087)*UCStop(1,3) + 
     &   At*(-2*AtC*Cd(10092)*Cr(5) + At*Cd(10095)*UCStopC(1,3))

        Cd(12530) = 2*UStop2(1,1)*(M_3 + 2*At*UStop2(1,1)) + 
     &   At*UStop2(1,2)**2

        Cd(1247) = UStop2(1,1)*
     &    (Cd(1246) - Cd(256)*UStop2(1,1)) - Cd(1245)*UStop2(1,2)

        Cd(12571) = MUE2*Cd(12506) - 
     &   MUE**2*Cd(11942)*UCStop(1,3)**2 - 
     &   MUEC**2*Cd(9972)*MStop2(1)*UCStopC(1,3)**2

        Cd(1252) = Cd(1251)*UStop2(1,1) - 2*UStop2(1,1)**3 + 
     &   Cd(519)*UStop2(1,2)**2

        Cd(1253) = -(Cd(511)*UStop2(1,1)**3) + 
     &   2*MStop2(3)*(Cd(1204)*UStop2(1,1) + Cr(85)*UStop2(1,2)**2)

	Cd(12625) = -MGl2**2 + Cd(48)**2 - MGlpTmSt4(2)*UStop2(1,1)

	Cd(12629) = Cd(251) + MGlpTmSt2(2)*UStop2(1,1)

        Cd(12640) = Cd(44)*Cd(9794)*UStop(1,1) + 
     &   2*MGl2*MTy*UUStop(1,1) + At*MGl2*UUStop(1,3)

        Cd(12654) = 8*M_3*Cd(164) + 8*M_3*MGlpTmSt2(1) - 
     &   At*MGlpTmSt2(1)

        Cd(12656) = 8*M_3*Cd(251) + 8*M_3*MGlpTmSt2(2) - 
     &   At*MGlpTmSt2(2)

	Cd(12667) = Cd(6268) - 2*Cd(12666)*UStop2(1,1)

	Cd(12682) = Cd(45) + 12*MStop2(1)*UStop2(1,1)

        Cd(1263) = 4*MStop2(3)**2*UStop2(1,2) + 
     &   MStop2(2)*(6*MStop2(3) + 4*Cd(594)*UStop2(1,2)**2)

        Cd(12720) = -Cd(12796) + 
     &   2*Cd(12700)*Cd(12701)*UStop2(1,1)*UStop2(1,2)

        Cd(12721) = (2*A0delStop(1)*Cd(12706) + 
     &     2*Cd(12720)*MStop2(3)**2 - Cd(12717)*Cr(13)*MStop2(3)**2
     &     )/MStop2(3)**2

        Cd(12723) = Cd(12700)*
     &    (Cd(12701)*Cr(75)*UCStop(3,1) + 
     &      Cd(12700)*UCStop(1,3)*UCStop(3,3)) + 
     &   Cd(12701)**2*UCStop(3,4)*UCStopC(1,3)

        Cd(12728) = Cd(12701)**2*UCStopC(1,3)*UCStopC(3,3) + 
     &   Cd(12700)*(Cd(12701)*Cr(75)*UCStopC(3,1) + 
     &      Cd(12700)*UCStop(1,3)*UCStopC(3,4))

	Cd(12747) = Cd(1359)*Cd(12700)*Cd(12701) + Cd(12796)

	Cd(12757) = AtC*(AtC + 2*MUE*TB) + MUE**2*TB2

	Cd(12758) = At*(At + 2*MUEC*TB) + MUEC**2*TB2

        Cd(12780) = -4*Cd(12701)*Cd(12779)*UCStop(3,1)*
     &    UCStopC(1,3) + 4*Cd(640)*Cd(12700)*UStop2(1,1)*
     &    UStop2(1,2)**2 - 
     &   4*Cd(12700)*Cd(12768)*MStop2(2)*UStop(1,2)*UUStop(3,1)

        Cd(12787) = MUEC**2*TB2*Cd(843) + 
     &   At*(At*Cd(843) + 2*MUEC*TB*MStop2(4))

	Cd(12794) = Cd(12796) + Cd(12700)*Cd(12701)*Cr(28)

        Cd(12797) = 4*Cd(485)*Cd(12700)*Cd(12701) - 
     &   4*MStop2(1)*MStop2(3) + 2*Cd(12796)*MStop2(4)

        Cd(12825) = Cd(12701)*Cd(12822)*UCStop(3,1)*
     &    UCStopC(1,3) + Cd(12700)*
     &    (Cd(12810)*MStop2(2)*UStop2(1,1)*UStop2(1,2) - 
     &      Cd(12824)*UUStopC(3,3))

        Cd(12828) = -(Cd(12700)*Cd(12827)) + 
     &   Cd(12758)*MStop2(1)*UCStopC(1,3)**2*UStop2(1,2) + 
     &   2*Cd(12809)*MStop2(3)**2*UStopC(1,1)

        Cd(12859) = 4*Cd(12712) - 2*Cd(12796)*MStop2(3) - 
     &   4*Cd(12714)*MStop2(1)*UStopC(1,2)

	Cd(12886) = Cd(898) + 2*Cr(85)*MStop2(1)*MStop2(2)

        Cd(1283) = -(Cd(1282)*UStop2(1,2)**2) + 
     &   MStop2(3)*(3*MStop2(3) - 3*MStop2(2)*UStop2(1,2))

	Cd(12900) = Cd(11656) + 4*UUStop(3,1)*UUStop(3,2)

        Cd(12909) = 2*(Cd(11106)*MStop2(2)**2*UStop2(1,1) + 
     &     MStop2(1)*(-(Cd(10925)*MStop2(2)*UStop2(1,2)) + 
     &        MStop2(1)*UStop2(1,1)*UUStop(3,3)))

        Cd(12910) = 2*Cd(11106)*MStop2(1)*MStop2(2)*
     &    UStop2(1,1) - 2*UStop2(1,2)*
     &    (Cd(10925)*MStop2(1)**2 + MStop2(2)**2*UUStop(3,4))

        Cd(12916) = UStop2(1,1)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2) - 
     &   Cd(10925)*UStop2(1,2)*UUStopC(3,3)

        Cd(12924) = UStop2(1,1)*
     &    (Cd(814) - Cd(7032)*UStop2(1,1)) + 
     &   MStop2(2)*UStop2(1,2)**2

        Cd(12940) = MUEC**2*TB2*Cd(12933) + 
     &   At*(At*Cd(12933) + 2*MUEC*TB*Cd(12934))

        Cd(1290) = MStop2(1)*(MStop2(1) - Cd(1237)*MStop2(2)) + 
     &   MStop2(2)**2*UStop2(1,1)

        Cd(12970) = MTy*Cd(12701)**2*UCStopC(1,3)*
     &    UStop2(1,1)**2 - 
     &   Cd(12700)*Cd(12967)*UCStop(1,3)*UStopC(1,2) + 
     &   Cd(12701)*Cd(12969)*UUStopC(1,3)

        Cd(1291) = -2*Cd(809)*MStop2(2)**2 + 
     &   MStop2(1)*(6*MStop2(1) - 4*Cr(85)*MStop2(2)*UStop2(1,2))

        Cd(12983) = Cd(12982)*UStop(1,2) + 
     &   AtC*Cr(91)*UStop(2,1) + MUE*TB*UStop(2,2)*UUStop(3,1)

        Cd(12993) = -(Cd(12701)**2*UCStopC(1,3)*UStop2(1,1)*
     &      UStop2(1,2)) + 
     &   Cd(12700)*(-(Cd(12700)*UCStop(1,3)**3) + 
     &      2*Cd(12701)*UCStop(1,3)*UStop2(1,1)*UStop2(1,2))

        Cd(12998) = Cd(12701)**2*UCStopC(1,3)**2 + 
     &   Cd(12700)*(-(Cd(12997)*UStop2(1,1)) + 
     &      Cd(12700)*UCStop(1,3)**2*UStop2(1,2))

        Cd(13021) = Cd(13020)*UStop(1,1) + 
     &   MGlpTmSt2(2)*(MUEC*TB*Cd(11106)*UStop(2,2) + 
     &      2*At*UStop(2,1)*UUStop(3,2))

        Cd(13034) = MUEC*TB*Cd(11106)*UStop(2,2) + 
     &   At*(UStop(1,1)*UUStop(2,2) + 2*UStop(2,1)*UUStop(3,2))

	Cd(13047) = 2*MTy2 + Cd(152)

        Cd(13051) = MTy*Cd(12701)**2*UCStopC(1,3)*
     &    UStop2(1,1)**2 - 
     &   Cd(12700)*Cd(12967)*UCStop(1,3)*UStopC(1,2) - 
     &   Cd(12701)*Cd(13050)*UUStopC(1,3)

        Cd(13054) = MTy*Cd(12700)**2*UCStop(1,3)*
     &    UStop2(1,2)**2 - 
     &   Cd(12701)*Cd(12961)*UCStopC(1,3)*UStopC(1,1) + 
     &   Cd(12700)*Cd(12963)*UUStopC(1,3)

        Cd(13059) = -(Cd(12701)*Cd(13004)*UCStopC(1,3)) + 
     &   MTy*Cd(12700)**2*Cr(75)*UCStop(1,3)*UStop2(1,2) - 
     &   Cd(12700)*Cd(13003)*UStopC(1,1)

        Cd(13062) = -(Cd(12700)*Cd(13009)*UCStop(1,3)**2) + 
     &   Cd(12701)*(Cd(13008) - 
     &      MTy*Cd(12701)*Cr(75)*UCStopC(1,3)*UStop2(1,2))

	Cd(13064) = Cd(179) + 4*MTy2*UStop2(1,2)

        Cd(13065) = 2*MTy2*Cr(28) + Cd(13064)*UStop2(1,1) + 
     &   Cd(179)*UStop2(1,2)

        Cd(13079) = 1/Pi**4*
     &   (MGl*MStop2(3)*UCStop(1,3)*
     &      (MTy*Cd(164)*Cd(13078) + MGlpTmSt2(1)*MStop2(1) + 
     &        MTy2*Cd(164)*UStop2(1,1)))/(MTy*MStop2(1))

	Cd(1303) = Cd(1213) - 3*UStop2(1,2)**2

        Cd(13102) = MStop2(2) + 
     &   MTy2*(Cr(28) + 2*UStop2(1,1)*UStop2(1,2))

	Cd(13108) = Cd(13422) - 2*UStop2(1,1)**2

	Cd(13109) = 7 + 2*UStop2(1,1)**2 + 2*UStop2(1,2)**2

	Cd(13123) = Cd(2384) + UStop2(1,1)**2

	Cd(13133) = Cd(11316) - UStop2(1,2)**2

	Cd(13134) = Cd(10918) - UStop2(1,2)**2

        Cd(13135) = UStop2(1,1)*(3*UStop2(1,1) - UStop2(1,2)) + 
     &   2*UStop2(1,2)**2

        Cd(13140) = 3*MStop2(2)*MStop2(3) - 
     &   Cd(951)*UStop2(1,1)**2 - Cd(842)*UStop2(1,2)**2

        Cd(1311) = 12*MStop2(3) + 8*MStop2(2)*UStop2(1,2) - 
     &   4*Cd(552)*UStop2(1,2)**2

        Cd(13180) = Cd(13178)*MStop2(1)*MStop2(3) - 
     &   Cd(12700)*Cd(13179)*UCStop(1,3) - 
     &   2*Cd(12701)*Cd(13177)*UStop2(1,2)*UUStop(3,2)

        Cd(13204) = 1/2.D0*
     &    (MStop2(1)*(Cr(85)*MStop2(1) - Cd(519)*MStop2(2))) + 
     &   MStop2(2)**2*UStop2(1,2)

        Cd(13213) = 1/2.D0*(Cd(13211)*MStop2(1)*MStop2(2)) - 
     &   Cd(13212)*MStop2(1)**2 + MStop2(2)**2*UUStop(3,4)

	Cd(129) = MGl2**2 - Cd(48)**2 + MGlpTmSt4(2)*UStop2(1,2)**2

        Cd(13256) = Cd(13255)*MStop2(1) - 
     &   MStop2(2)*UUStop(1,1) - MGlpT2*UUStop(2,1)

        Cd(13258) = Cd(13257)*MStop2(1) - 
     &   MStop2(2)*UUStop(1,2) - MGlpT2*UUStop(2,2)

        Cd(13259) = Cd(9335)*MStop2(1) - MStop2(2)*UUStop(1,3) - 
     &   MGlpT2*UUStop(2,3)

        Cd(13299) = -(Cd(13298)*UStop(1,1)) + 
     &   MStop2(2)*(MUEC*TB*Cd(11106)*UStop(2,2) + 
     &      2*At*UStop(2,1)*UUStop(3,2))

	Cd(130) = Cd(251) + MGlpTmSt2(2)*UStop2(1,2)

        Cd(13316) = -Cd(12796) + 
     &   2*Cd(12700)*Cd(12701)*UStop2(1,1)*UStop2(1,2)**2

        Cd(1326) = 8*MStop2(3)**2*UStop2(1,2) + 
     &   MStop2(2)*(5*MStop2(3) + 4*Cd(1325)*UStop2(1,2)**2)

        Cd(13331) = 2*Cd(12701)*Cd(13330) + 
     &   MGl*(-(Cd(12700)*Cd(13009)*UCStop(1,3)**2) - 
     &      MTy*Cd(12701)**2*Cr(75)*UCStopC(1,3)*UStop2(1,2))

        Cd(13342) = 1/(2.D0*Pi**4)*
     &   (MGl*MStop2(3)*UCStop(1,3)*
     &      (MTy*Cd(164)*Cd(13078) + MGlpTmSt2(1)*MStop2(1) + 
     &        MTy2*Cd(164)*UStop2(1,1)))/(MTy*MStop2(1))

        Cd(13386) = UStop2(1,1)*
     &    (6*UStop2(1,1) - 2*UStop2(1,2)) + UStop2(1,2)**2

	Cd(13400) = Cd(1213) - UStop2(1,1)**3*UStop2(1,2)

        Cd(13401) = Cd(617)*UStop2(1,1)*UStop2(1,2) - 
     &   Cr(26)*UStop2(1,2)**2 - Cd(13398)*UUStopC(3,3)

	Cd(1336) = Cd(1271) - 8*Cd(1333) - 2*Cd(1335)*UStop2(1,1)

        Cd(13430) = -2*Cd(13429)*UStop2(1,1) + 
     &   Cd(6281)*UStop2(1,1)**2 + UStop2(1,2)**2

        Cd(13438) = Cd(12701)*
     &    (Cd(13242) - MTy*Cd(12701)*Cd(13234)*UStop(1,2)) + 
     &   2*Cd(12700)*Cd(13236)*UCStop(1,3)*UStopC(1,2)

        Cd(13439) = Cd(9335)*MStop2(2) - MStop2(1)*UUStop(1,3) - 
     &   MGlpT2*UUStop(2,3)

        Cd(13441) = Cd(13440)*MStop2(3) + 
     &   2*MGlpTmSt2(2)*UCStopC(3,1)**2 - 2*Cd(13439)*UUStopC(1,3)

        Cd(13483) = 1/6.D0*
     &   (-12*A0delq(MSbot2(1),Q) + (6 + Pi**2)*MSbot2(1))/
     &    MSbot2(1)

	Cd(13497) = Cd(13490) + Cd(13496)*MSq2Diff(1,2)

        Cd(1346) = 2*Cd(1345)*MStop2(1) - 
     &   4*Cd(169)*MStop2(2)**2 + Cd(1344)*UStop2(1,2)**2

	Cd(13541) = Cd(13479) - MSq2Diff(1,1)

	Cd(13550) = Cd(13479) - MSq2Diff(1,2)

	Cd(13567) = Cd(13562) + MSq2Diff(1,2)

        Cd(1351) = -(Cd(1350)*UStop2(1,2)**2) + 
     &   MStop2(3)*(22*MStop2(3) - 10*MStop2(2)*UStop2(1,2))

        Cd(13575) = A0delq(MSbot2(1),Q)*Cd(13574) - 
     &   A0delStop(2)*Cd(13561)*MSq2Diff(1,1) - 
     &   A0delStop(1)*Cd(13560)*MSq2Diff(1,2)

        Cd(13576) = 3/(4.D0*Pi)*
     &   (Alfa2L*(-(Cd(13574)*MSbot2(1)) + 
     &        Cd(13560)*MSq2Diff(1,2)*MStop2(1) + 
     &        Cd(13561)*MSq2Diff(1,1)*MStop2(2)))/
     &    (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2))

        Cd(13577) = 3/(8.D0*Pi)*
     &   (Alfa2L*(-(Cd(13466)*Cd(13574)*MSbot2(1)) + 
     &        Cd(13560)*Cr(13)*MSq2Diff(1,2)*MStop2(1) + 
     &        Cd(13561)*Cr(34)*MSq2Diff(1,1)*MStop2(2)))/
     &    (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2))

        Cd(13580) = A0delq(MSbot2(1),Q)*Cd(13579) - 
     &   A0delStop(2)*Cd(13578)*MSq2Diff(1,1) - 
     &   A0delStop(1)*Cd(13508)*MSq2Diff(1,2)

        Cd(1352) = Cd(1351) + Cd(1346)*UStop2(1,1) + 
     &   4*MStop2(2)**2*UStop2(1,1)**3

        Cd(13581) = 3/(4.D0*Pi)*
     &   (Alfa2L*(-(Cd(13579)*MSbot2(1)) + 
     &        Cd(13508)*MSq2Diff(1,2)*MStop2(1) + 
     &        Cd(13578)*MSq2Diff(1,1)*MStop2(2)))/
     &    (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2))

        Cd(13582) = 3/(8.D0*Pi)*
     &   (Alfa2L*(-(Cd(13466)*Cd(13579)*MSbot2(1)) + 
     &        Cd(13508)*Cr(13)*MSq2Diff(1,2)*MStop2(1) + 
     &        Cd(13578)*Cr(34)*MSq2Diff(1,1)*MStop2(2)))/
     &    (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2))

        Cd(13587) = Cd(12700)*Cd(13585) + MTy*Cd(13586) + 
     &   2*MSq2Diff(1,2)

        Cd(1354) = Cd(1352)*UStop2(1,1) + Cd(1353)*UStop2(1,2) + 
     &   4*MStop2(2)**2*UStop2(1,2)**4

        Cd(1355) = 12*(Cr(53)*MStop2(1) - MStop2(2))*MStop2(2) + 
     &   11*MStop2(3)*MStop2(4)

	Cd(13632) = Cd(11166) + MUEC*TB*Cr(91)

	Cd(13639) = Cd(1898) - 3*MSbot2(1)

        Cd(13644) = -(Cd(12902)*MSbot2(1)) + 
     &   Cd(11106)*MStop2(1) + 2*MStop2(2)*UUStop(3,3)

	Cd(1360) = Cd(1179) + UStop2(1,2)

        Cd(1361) = Cd(1359)*MStop2(2)**2 + 
     &   MStop2(1)*(11*MStop2(1) - 6*Cd(1360)*MStop2(2))

	Cd(13695) = MGl2**2 - MTy2**2 + 2*Cd(153)*MSbot2(1)

	Cd(13701) = MGl2 - Cd(14069)

	Cd(13702) = 4*MGl2 + Cd(13635)

	Cd(13703) = MGl2*Cd(13701) + Cd(13702)*MSbot2(1) + MSbot2(1)**2

	Cd(13714) = At*Cd(12700) + MUEC*TB*Cd(12700) - MSq2Diff(1,2)

	Cd(13716) = Cd(179)*(AtC + MUE*TB - Cd(12700)*UStop2(1,1))

	Cd(13725) = Cd(11574) - 2*UStop2(1,1)

	Cd(13734) = -2*MSbot2(1) + MStop2(2)

	Cd(13736) = Cd(2364) - 2*UStop2(1,1)

	Cd(13744) = UStop2(1,1)*(2 + UStop2(1,1)) + UStop2(1,2)**2

        Cd(13753) = AtC*MSq2Diff(1,2)*UCStop(1,3) + 
     &   MUE*TB*MSq2Diff(1,2)*UCStop(1,3) + Cd(13752)*UCStopC(1,3)

        Cd(13754) = 2*MSq2Diff(1,2)*
     &    (Cd(585) + MUE*TB*UCStop(1,3)) + 
     &   MUEC*TB*Cd(13719)*UCStopC(1,3)

        Cd(13755) = 2*MSq2Diff(1,2)*
     &    (Cd(585) + MUE*TB*UCStop(1,3)) - 
     &   MUEC*TB*Cd(13734)*UCStopC(1,3)

	Cd(13759) = MGl2*MGlmT2 + MStop2(1)*(-Cd(153) + MStop2(1))

	Cd(13771) = Cd(12701) + MTy*UCStop(1,3)

        Cd(13772) = AtC*Cd(13469) + MTy*MUE*TB*UStop(1,1) + 
     &   Cd(13541)*UStop(1,2)

        Cd(13776) = Cd(10976)*UStop(1,2) - MUE*TB*UUStop(1,1) + 
     &   MUEC*TB*UUStop(1,2)

	Cd(13784) = Cd(10918) + UStop2(1,2)**2

	Cd(1373) = Cd(2116) + 4*Cd(286)*UStop2(1,2)**2

	Cd(13793) = (6 - 2*UStop2(1,1))*UStop2(1,1) + 2*UStop2(1,2)**2

	Cd(13794) = (2 - UStop2(1,1))*UStop2(1,1) + UStop2(1,2)**2

	Cd(13797) = Cd(2364) + UStop2(1,1)**2

	Cd(13801) = Cd(585) + MUE*TB*UCStop(1,3)

        Cd(13807) = AtC*Cd(13806) + 
     &   MUE*TB*Cd(13805)*MSq2Diff(1,2) - 
     &   Cd(13077)*MSbot2(1)*MSq2Diff(1,1)*UStopC(1,1)

        Cd(13818) = -(Cd(13817)*MSbot2(1)) + 
     &   Cd(11877)*MSbot2(1)**2 - 2*MSq2Diff(1,2)*MStop2(1)

	Cd(13825) = Cd(153) - MStop2(2)

        Cd(13842) = Cr(87)*MSbot2(1) + Cd(13784)*MStop2(2) - 
     &   Cd(13841)*UStop2(1,1)**2

	Cd(13843) = Cd(12810) + 3*UStop2(1,1)

	Cd(13846) = Cd(13996) + AtC*UCStop(1,3)

        Cd(1379) = 3*UStop2(1,1)**4 + 
     &   UStop2(1,2)**2*(-8 + 3*UStop2(1,2)**2)

	Cd(13852) = Cd(153) - MStop2(1)

	Cd(13855) = AtC + Cd(13507)

	Cd(13858) = (At*(AtC + MUE*TB) - MSq2Diff(1,1))*MSq2Diff(1,1)

	Cd(13862) = (3 - UStop2(1,1))*UStop2(1,1) + UStop2(1,2)**2

        Cd(13864) = Cd(10367)*Cd(13863) + 
     &   MStop2(1)*(-2*Cd(13794)*MSbot2(1) + Cd(13862)*MStop2(1))

        Cd(13867) = MGl2*Cd(13866) + 
     &   MSq2Diff(1,1)*(MGl2**2 - Cd(44)*MStop2(1))

	Cd(13884) = Cd(14341) + Cd(13663)*UUStopC(1,3)

	Cd(13896) = Cd(13841) - MSq2Diff(1,1)

	Cd(13903) = MTy2 + Cd(13841)

	Cd(13907) = 2*M_3C*At + 2*M_3*AtC + MStop2(1)

	Cd(12) = 2*MTy - Cd(585)

	Cd(136) = 2*MTy2 + Cd(138)

	Cd(13915) = Cd(13642) - MSbot2(1)

	Cd(13916) = MGl2*MTy + MTy*Cd(13915) + 2*M_3*MTy2*UCStopC(1,3)

	Cd(13917) = 2*MTy*Cr(5) + MGlpT2*UStop2(1,1)

        Cd(13918) = Cd(13916)*UCStopC(1,3) + 
     &   UStop2(1,1)*(-(AtC*Cd(13917)) + 2*M_3C*MTy2*UStop2(1,2))

	Cd(13928) = MGlpT2 - At*MUE*TB + MSbot2(1)

	Cd(13941) = Cd(11515) + UStop2(1,1)

        Cd(1389) = -2*MStop2(2)*
     &    (Cd(1388)*MStop2(1) + Cd(1387)*MStop2(2)) + 
     &   8*MStop2(1)**2*UStop2(1,1)

	Cd(13947) = UStop2(1,2)*(1 + UStop2(1,2))

        Cd(13948) = -(MGl2*Cd(13737)) + 
     &   UStop2(1,1)*(At*Cd(12700)*Cr(56) + Cd(13726)*UStop2(1,1))

	Cd(13951) = Cd(537) + 3*UStop2(1,1)

	Cd(13963) = -(At*MUE*TB) + Cd(8515)

	Cd(13964) = MTy2*Cd(13962) - Cd(13963)*MSbot2(1) - MSbot2(1)**2

	Cd(13971) = Cd(11515) - UStop2(1,1)

	Cd(13973) = Cd(11341) + 2*UStop2(1,1)**2

	Cd(13983) = Cd(1387) + 4*UStop2(1,1)**2

	Cd(13984) = Cd(970) + 4*UStop2(1,1)**2

	Cd(13997) = MUEC*TB*MSbot2(1) + At*(2*MSbot2(1) - MStop2(2))

	Cd(13999) = MUEC*TB*(2*MSbot2(1) - MStop2(2)) + At*MStop2(2)

        Cd(14010) = MTy*MUE*TB*UCStop(1,3) + 
     &   Cd(13541)*UStop2(1,2) + AtC*Cd(13469)*UStopC(1,2)

	Cd(14020) = At*(AtC + MUE*TB) - MSbot2(1)

        Cd(14024) = Cd(10648)*MSbot2(1) + 4*MStop2(2) + 
     &   At*Cd(12700)*UStop2(1,2)

        Cd(14025) = MStop2(1)*(Cd(14024) - Cd(807)*MStop2(1)) - 
     &   2*MSbot2(1)*MStop2(2)

        Cd(14026) = 3*MStop2(2) + MSbot2(1)*UStop2(1,1) + 
     &   At*Cd(12700)*UStop2(1,2)

        Cd(14027) = MStop2(1)*(Cd(14026) - Cd(169)*MStop2(1)) - 
     &   MSbot2(1)*MStop2(2)

	Cd(14029) = Cd(14031) + AtC*UCStop(1,3)*UStop(1,1)

        Cd(14041) = 2*MTy2*Cd(13962) + 
     &   (Cd(13963) - 2*MStop2(1))*MStop2(1)

	Cd(14043) = Cd(13698) + 2*MStop2(1)

	Cd(1399) = Cd(1311) - Cd(1397) - 4*MStop2(2)*UStop2(1,1)**2

	Cd(14077) = MGl2**2 - MTy2**2 - MGlmT2*MStop2(1)

	Cd(14080) = MGl2**2 - Cd(4243)

	Cd(14084) = At*MUE*TB2 + MUE2*TB*TB2 - TB*MSq2Diff(1,1)

	Cd(14098) = At*AtC + Cd(13698)

	Cd(14100) = Cd(14316) - At*UStop(1,2)*UStop2(1,1)

        Cd(14105) = Cd(14104) - 
     &   Cd(12700)*MSq2Diff(1,2)*UCStop(1,3) + 
     &   MTy2*Cd(13504)*UStop2(1,2)*UStopC(2,1)

        Cd(14106) = MTy*Cd(14105) + MSq2Diff(1,2)**2 + 
     &   MTy2*Cd(12700)*Cd(13504)*UStop2(1,1)*UStopC(2,2)

	Cd(14126) = AtC*MUEC*TB2 + MUE2*TB*TB2 - TB*MSq2Diff(1,1)

	Cd(14130) = MTy2 - 2*MTy*Cr(5) + MGl2*UStop2(1,1)

        Cd(14135) = MGl2*Cd(14131) - MTy*Cd(14133) + 
     &   AtC*TB*Cd(14134)*UStop2(1,2)

        Cd(14143) = 2*M_3C*At + 2*M_3*AtC + 
     &   MSq2Diff(1,1)*UStop2(1,1)

        Cd(14144) = Cd(9948)*Cd(13896)*UStop2(1,1) + 
     &   UStop2(1,2)*(Cd(14143) + Cd(13841)*UStop2(1,2))

        Cd(14149) = 2*M_3*Cd(9950) + 
     &   At*(Cd(9950) + MGlpTmSt2(1)*MStop2(1))*UStop2(1,1)

        Cd(14153) = 2*M_3C*Cd(9950) + 
     &   AtC*(Cd(9950) + MGlpTmSt2(1)*MStop2(1))*UStop2(1,1)

        Cd(14166) = AtC*TB*UCStop(1,3) + MUE*TB2*UCStop(1,3) + 
     &   Cd(14147)*UCStopC(1,3)

        Cd(14172) = Cd(14171)*UCStopC(1,3) + 
     &   UCStop(1,3)*(AtC - 2*M_3C*UStop2(1,1))

	Cd(14189) = MTy2 + 2*MTy*Cr(5) + MGl2*UStop2(1,1)

	Cd(14191) = Cd(10669) + 2*M_3*MTy*UCStopC(1,3)

        Cd(14196) = Cd(9948)*MSbot2(1) - MSq2Diff(1,2) + 
     &   MSbot2(1)*UStop2(1,1)

	Cd(14217) = MGl2 - Cd(14235)

	Cd(14259) = -Cd(14255) + MSbot2(1)*UCStop(1,3)

        Cd(14260) = -2*At*Cr(5)*Cr(75)*MSq2Diff(1,2)**2 + 
     &   MStop2(3)*(Cd(14256) + 
     &      MSbot2(1)*MSq2Diff(1,1)*UCStop(1,3))

        Cd(14267) = 3*Cd(13696) - 2*MGl2*Log(MGl2/Q) - 
     &   2*Log(MSbot2(1)/Q)*MSbot2(1)

	Cd(14274) = Cd(13696) - 2*MStop2(1)

        Cd(14279) = MSbot2(1)*(-Cd(155) + 2*MSbot2(1)) + 
     &   2*MTy2*MStop2(1)

	Cd(14280) = MSbot2(1)*(-Cd(45) + MSbot2(1)) + 2*MTy2*MStop2(1)

	Cd(1423) = Cd(1271) - 8*Cd(1333) + Cd(1422)*MStop2(3)

	Cd(14283) = Cd(1496) - MSbot2(1)

	Cd(14285) = MGl2*Cd(13642) + Cd(14292)

	Cd(1424) = 4*Cd(1268) - 24*Cd(1280) - Cd(1357)*MStop2(3)

        Cd(14302) = MTy*Cd(14301)*UStop(1,1) + 
     &   5*Cd(12701)*MSq2Diff(1,2)*UStop(1,2) + 
     &   6*Cd(13469)*MSbot2(1)*UStop2(1,1)

        Cd(14304) = 2*MTy*Cd(14303)*UStop(1,1) + 
     &   6*Cd(12701)*MSq2Diff(1,2)*UStop(1,2) + 
     &   6*Cd(13469)*MSbot2(1)*UStop2(1,1)

        Cd(14317) = At*MTy*Cd(14315)*UStop(1,2) + 
     &   MSbot2(1)*(MTy*Cd(14316) + 
     &      At*MUE*TB*UStop(1,1)*UStop2(1,2))

        Cd(14318) = MUEC*TB*Cd(14312) + 2*Cd(14317) + 
     &   2*Cd(14314)*UStop(1,1)

	Cd(14329) = Cd(14308) - Cd(14309) + MTy*MSq2Diff(1,2)

	Cd(1428) = Cd(1136) + Cr(7)*UStop2(1,1)

        Cd(1441) = 4*Cd(1439) + 5*Cr(28)*MStop2(3) + 
     &   2*Cd(1440)*UStop2(1,2)

        Cd(1443) = 5*(Cd(290)*MStop2(2)**2 + 
     &     MStop2(1)*(MStop2(1) - 2*MStop2(2)*UStop2(1,1)))

	Cd(1453) = -8*Cd(852) + 4*Cd(1264) - Cd(1421)*MStop2(3)

	Cd(1456) = 2*MStop2(2) - 6*MStop2(3) + Cd(1442)*UStop2(1,2)**2

	Cd(1459) = Cd(1213) - 2*UStop2(1,2)**2

	Cd(1471) = -4*Cd(523) + 4*Cd(1390) + Cd(1468)*MStop2(3)

	Cd(1476) = Cd(1311) - Cd(1397) + Cd(1474)

	Cd(1494) = Cd(1496) - 2*MStop2(1)

	Cd(1502) = -2*MStop2(1)**2 + MTy2*(MTy2 + 3*MStop2(1))

	Cd(1503) = Cd(1647) + 6*MStop2(1)

	Cd(1509) = 2*MTy2 + Cd(99)

	Cd(1512) = Cd(6362) - 2*MStop2(1)

	Cd(1517) = 9*MGl2 + Cd(1739)

	Cd(1525) = MGl2 + Cd(1636)

	Cd(1526) = Cd(48)**3 + MGl2*(3*Cd(48)**2 + MGl2*Cd(1525))

	Cd(1529) = Cd(2157) + MStop2(3)

	Cd(1530) = 8*MTy2 + Cd(1749)

	Cd(1531) = 6*Cd(637)*MStop2(2) + MTy2*(16*MTy2 - 6*MStop2(3))

        Cd(1534) = 2*MTy2*Cd(1532) + 
     &   2*(Cd(1533) - 13*MStop2(2))*MStop2(2)

	Cd(1540) = 8*MTy2 + Cd(637)

	Cd(1543) = MTy2*Cd(1541) + MStop2(2)*(Cd(1542) + 9*MStop2(2))

	Cd(1544) = MTy2*(24*MTy2 - 3*MStop2(1)) + 2*MStop2(1)**2

	Cd(1548) = MGl2*(8*MGl2 - 2*Cd(1530)) + Cd(1531)

	Cd(1552) = MGl2*(2*MGl2 - Cd(1540)) + Cd(1543)

        Cd(1558) = MTy2*(3*MStop2(1) - 19*MStop2(2)) + 
     &   Cd(1557)*MStop2(2)

	Cd(1560) = 8*MTy2 - Cd(1122)

	Cd(1561) = Cd(2122) - MStop2(3)

	Cd(1562) = 4*MTy2 - Cd(1657)

	Cd(1563) = MGl2*(4*MGl2 - Cd(1556)) - Cd(1559)*Cd(1560)

	Cd(1567) = 6*MTy2 - Cd(1226)

	Cd(1568) = Cd(2875) + 18*MStop2(2)

	Cd(1571) = MTy2*(3*MTy2 + Cd(1569)) - Cd(1570)*MStop2(2)

	Cd(1582) = Cd(1718) + 16*MStop2(1)

	Cd(156) = MGl2*(-MGl2 + Cd(155)) + Cd(44)*MStop2(1)

	Cd(1591) = 8*MTy2 + Cd(1583)

	Cd(1598) = MGl2*(8*MGl2 - 2*Cd(1582)) + Cd(1585)

	Cd(1606) = MGl2 + Cd(1659)

	Cd(1607) = Cd(44)**3 + MGl2*(3*Cd(44)**2 + MGl2*Cd(1606))

	Cd(1613) = 8*MTy2 - Cd(1421)

	Cd(1614) = 19*MTy2 - 24*MStop2(2) + 15*MStop2(3)

	Cd(1619) = Cd(3233) + 18*MStop2(1)

	Cd(1620) = 6*MTy2 - 6*MStop2(2) + 5*MStop2(3)

	Cd(1644) = Cd(1636) + 2*MStop2(3)

	Cd(1646) = Cd(3233) + 12*MStop2(1)

	Cd(1660) = MGl2*(MGl2 - Cd(874)) - Cd(89)*Cd(1559)

	Cd(1666) = MTy2*(4*MTy2 + 2*Cd(315)) - 2*Cd(874)*MStop2(2)

	Cd(1669) = Cd(2875) + 12*MStop2(2)

	Cd(1670) = 3*Cd(1668) + MGl2*(3*MGl2 - Cd(1669))

        Cd(15) = 1/6.D0*(-12*A0delStop(2) + (6 + Pi**2)*MStop2(2))/
     &    MStop2(2)

	Cd(166) = MGl2*Cd(165) + Cd(44)*MStop2(1) - Cd(156)*UStop2(1,1)

	Cd(1686) = MTy2*(MTy2 - MStop2(1)) + 2*MStop2(1)**2

	Cd(1688) = -MStop2(1)**2 + MTy2*(MTy2 + 2*MStop2(1))

	Cd(1690) = Cd(1845) + MStop2(3)

	Cd(1691) = -(Cd(1594)*MStop2(1)) + MTy2*(8*MTy2 - MStop2(3))

	Cd(1701) = Cd(1647) - 2*MStop2(1)

        Cd(1703) = MTy2*(3*MStop2(1) - 11*MStop2(2)) + 
     &   Cd(1702)*MStop2(2)

	Cd(1704) = MGl2*(-2*MGl2 + Cd(1122)) + Cd(1559)*Cd(1562)

	Cd(1713) = -6*MStop2(1)**2 + MTy2*(3*MTy2 + 8*MStop2(1))

	Cd(1714) = MTy2*(9*MTy2 - 20*MStop2(1)) + 6*MStop2(1)**2

	Cd(1717) = Cd(1845) + 4*MStop2(3)

	Cd(170) = Cd(244) + Cd(45)*UStop2(1,1)

	Cd(1731) = Cd(1845) - 2*MStop2(3)

	Cd(1734) = MTy2*(MTy2 + 2*MStop2(1)) + Cd(602)*MStop2(2)

        Cd(171) = -(Cd(169)*MStop2(1)**2) + 
     &   2*MStop2(1)*MStop2(2) + Cd(44)**2*UStop2(1,1)

	Cd(1737) = MTy2*(8*MTy2 - MStop2(1)) + 2*MStop2(1)**2

	Cd(1748) = 3*MGl2 + Cd(1747)

	Cd(1767) = MTy2*(12*MTy2 + Cd(637)) + 3*Cd(874)*MStop2(2)

	Cd(1773) = 2*MTy2 + Cd(1226)

	Cd(1774) = MTy2*(MTy2 + Cd(315)) - Cd(535)*MStop2(2)

	Cd(1782) = MGlmT2**2 - 2*MGlpT2*MStop2(2) + MStop2(2)**2

	Cd(1786) = 5*MStop2(1)**2 + MTy2*(MTy2 + 2*MStop2(1))

	Cd(1787) = Cd(3006) + 6*MTy2*MStop2(1)

	Cd(1813) = MTy2*(12*MTy2 + Cd(1583)) + 3*Cd(607)*MStop2(1)

	Cd(1825) = 5*MTy2 - Cd(2802)

	Cd(1827) = 3*MTy2**2 - 3*MGlpT2*MStop2(1) + MStop2(1)**2

	Cd(1830) = 10*MTy2 - Cd(2803)

	Cd(1848) = MTy2*(2*MTy2 - 6*MStop2(1)) + 3*MStop2(1)**2

	Cd(1851) = MTy2*(4*MTy2 - 6*MStop2(1)) + 3*MStop2(1)**2

	Cd(1867) = -5*MStop2(1)**2 + MTy2*(4*MTy2 + 7*MStop2(1))

        Cd(1868) = -4*Cd(1867)*MStop2(2) + 
     &   6*Cd(45)*(Cd(1866) + MStop2(2)**2)

	Cd(1869) = MStop2(1)**2 + MTy2*(4*MTy2 + 8*MStop2(1))

	Cd(1870) = MTy2*(4*MTy2 - 3*MStop2(1)) + 9*MStop2(1)**2

        Cd(1871) = 2*Cd(45)*Cd(1869) - 
     &   2*MStop2(2)*(2*Cd(1870) + 3*Cd(45)*MStop2(2))

	Cd(1887) = MGl2*(MGl2 - 2*MTy2) + 2*MTy2**2

	Cd(1889) = 2*Cd(1887) + MStop2(1)*(-Cd(1888) + 9*MStop2(1))

	Cd(187) = MTy2 - Cd(256)

	Cd(1894) = -9*MStop2(1)**2 + MTy2*(3*MTy2 + 10*MStop2(1))

	Cd(1897) = Cd(1592) + MStop2(3)

	Cd(1899) = -(Cd(1898)*MStop2(1)) + MTy2*(MTy2 + 2*MStop2(1))

	Cd(1911) = 17*MGl2 + 14*MTy2 - 2*MStop2(1)

	Cd(1917) = -4*MTy2 + Cd(121)

	Cd(1919) = 11*MGl2 + 12*MTy2 - 4*MStop2(1)

	Cd(1922) = MGl2*(MGl2 - MTy2) - MTy2**2

	Cd(1925) = 27*MGl2 + Cd(1546)

	Cd(1928) = MTy2**2 + (2*Cd(44) - MStop2(2))*MStop2(2)

	Cd(1942) = Cd(48) - 2*MStop2(3)

	Cd(1943) = -9*MStop2(1)**2 + MTy2*(3*MTy2 + 14*MStop2(1))

	Cd(1947) = Cd(48) - 7*MStop2(3)

	Cd(1950) = -3*MStop2(1)**2 + MTy2*(13*MTy2 + 6*MStop2(1))

	Cd(1956) = Cd(2309) + 6*MStop2(1)

	Cd(194) = MGl2*(-MGl2 + Cd(181)) + Cd(48)*MStop2(2)

	Cd(1972) = 14*MTy2 - 14*MStop2(2) + 9*MStop2(3)

	Cd(1980) = Cd(2774) + 9*MStop2(1)

	Cd(1983) = Cd(164) + MGlpTmSt2(1)*UStop2(1,2)**2

	Cd(1987) = MGl2 - 5*MTy2 - 6*MStop2(1)

	Cd(1989) = MGl2*(MGl2 - 3*MStop2(1)) + 3*MStop2(1)**2

	Cd(2000) = MTy2*Cd(1998) + (Cd(1999) - 18*MStop2(2))*MStop2(2)

	Cd(2006) = Cd(251) + MGlpTmSt2(2)*UStop2(1,2)**2

	Cd(2010) = 13*MGl2 + 9*MTy2 - 6*MStop2(1)

	Cd(2012) = MGl2*(MGl2 - 3*MStop2(2)) + 3*MStop2(2)**2

	Cd(2023) = -13*MStop2(1)**2 + MTy2*(8*MTy2 + 28*MStop2(1))

	Cd(2024) = MTy2*(18*MTy2 - 13*MStop2(1)) + 3*MStop2(1)**2

	Cd(2025) = -MStop2(1)**2 + MTy2*(7*MTy2 + 8*MStop2(1))

	Cd(2028) = MGlmT2**2 - MGlpT2*MStop2(1) - MGlpTmSt2(1)*MStop2(2)

	Cd(2030) = 2*MGl2 + Cd(2973)

	Cd(2041) = 3*MGl2 + Cd(2926)

	Cd(2049) = Cd(1866) + 5*MTy2*MStop2(1)

	Cd(2063) = Cd(1667) - 2*MStop2(2)

	Cd(2069) = 11*MTy2 + 3*(-4*MStop2(2) + MStop2(3))

	Cd(2073) = 7*MTy2 - Cd(1148)

	Cd(2079) = MTy2*(6*MTy2 - 4*MStop2(1)) + 3*MStop2(1)**2

	Cd(2080) = MTy2*(2*MTy2 - 10*MStop2(1)) + 3*MStop2(1)**2

	Cd(206) = Cd(48) + 2*MStop2(3)

	Cd(2085) = 8*MTy2**2 + 6*MGlpT2*MStop2(1) - 3*MStop2(1)**2

	Cd(2087) = Cd(529)**2 + MStop2(2)*(4*MTy2 + 13*MStop2(2))

	Cd(2088) = Cd(1845) - MStop2(3)

        Cd(2099) = MTy2*(5*MStop2(1) - 13*MStop2(2)) + 
     &   Cd(2098)*MStop2(2)

	Cd(2101) = 3*MTy2*MStop2(1) - Cd(2100)*MStop2(2) - MStop2(2)**2

	Cd(2102) = 5*MGl2 + Cd(2100)

	Cd(2113) = Cd(1684) + MStop2(3)

	Cd(2117) = 2*Cd(2116)*MStop2(4) + MTy2*(8*MTy2 + 8*MStop2(4))

	Cd(2118) = 9*MStop2(1)**2 + MTy2*(2*MTy2 + 3*MStop2(1))

	Cd(2120) = MTy2*(2*MTy2 - MStop2(1)) - 7*MStop2(1)**2

	Cd(2123) = MTy2*(2*MTy2 + 2*Cd(814)) - 6*MStop2(1)**2

	Cd(2134) = 2*MTy2 - Cd(1432)

	Cd(2137) = Cd(2157) + 7*MStop2(3)

	Cd(2138) = 4*MTy2 + Cd(4225)

	Cd(2148) = 2*Cd(1022)*MStop2(4) + MTy2*(8*MTy2 + 8*MStop2(4))

	Cd(2149) = MTy2*(4*MTy2 - 4*MStop2(1)) + 3*MStop2(1)**2

	Cd(2150) = MTy2*(8*MTy2 - 6*MStop2(1)) + 3*MStop2(1)**2

	Cd(2153) = MTy2*(6*MTy2 - 4*MStop2(1)) - 3*MStop2(1)**2

	Cd(2154) = -3*MStop2(1)**2 + MTy2*(2*MTy2 + 2*MStop2(1))

	Cd(2158) = 2*MTy2*Cd(44) + 6*(MTy2 - MStop2(2))*MStop2(2)

	Cd(2168) = 2*MTy2 - Cd(1160)

	Cd(2170) = Cd(2157) - 3*MStop2(3)

        Cd(2172) = 3*MTy2*MStop2(1) + 
     &   (Cd(2171) - 17*MStop2(2))*MStop2(2)

	Cd(2173) = Cd(2930) + 12*MStop2(2)

	Cd(2174) = MTy2*(2*MTy2 + Cd(1141)) - Cd(1226)*MStop2(2)

	Cd(2186) = 3*Cd(44)**2 + MGl2*(MGl2 + 3*Cd(44))

	Cd(2192) = MTy2*(2*MTy2 - MStop2(1)) + 9*MStop2(1)**2

	Cd(2195) = 8*MTy2 + Cd(2098)

	Cd(2197) = 3*Cd(607)*MStop2(1) + MTy2*(20*MTy2 - 3*MStop2(3))

	Cd(2198) = -3*MStop2(1)**2 + MTy2*(5*MTy2 + 6*MStop2(1))

	Cd(218) = Cd(155) + M_3*MTy*UCStopC(1,3)

	Cd(2206) = 4*MTy2 - Cd(1442)

	Cd(2207) = Cd(2186) + MGl2*(MGl2 - 3*Cd(2200))

	Cd(2208) = Cd(44)**3 + MGl2*(-3*Cd(44)*Cd(2200) + Cd(2207))

	Cd(2224) = MTy2*(-8*MTy2 + Cd(1583)) - 3*MStop2(1)*MStop2(3)

	Cd(2227) = 2*Cd(2186) + MGl2*(MGl2 - 3*Cd(2200))

	Cd(2228) = Cd(44)**3 + MGl2*(-3*Cd(44)*Cd(2200) + Cd(2227))

	Cd(2240) = 3*Cd(48)**2 + MGl2*(MGl2 + 3*Cd(48))

	Cd(2246) = MTy2*(4*MTy2 - 4*MStop2(1)) - MStop2(1)**2

	Cd(2247) = 3*MStop2(1)**2 + MTy2*(8*MTy2 + 2*MStop2(1))

	Cd(2250) = 3*Cd(874)*MStop2(2) + MTy2*(20*MTy2 + 3*MStop2(3))

	Cd(2256) = 4*MTy2 - Cd(1141)

	Cd(2260) = Cd(2240) + MGl2*(MGl2 - 3*Cd(2254))

	Cd(2261) = Cd(48)**3 + MGl2*(-3*Cd(48)*Cd(2254) + Cd(2260))

	Cd(224) = Cd(181) - M_3*MTy*UCStopC(1,3)

	Cd(2276) = MTy2*Cd(2274) + (Cd(2275) - 17*MStop2(2))*MStop2(2)

        Cd(2278) = 5*MTy2*MStop2(1) + 
     &   (Cd(2277) - 15*MStop2(2))*MStop2(2)

	Cd(2280) = MTy2*(-8*MTy2 + Cd(637)) + 3*MStop2(2)*MStop2(3)

	Cd(2282) = 2*Cd(2240) + MGl2*(MGl2 - 3*Cd(2254))

        Cd(21) = 1/6.D0*(-12*A0delStop(1) + (6 + Pi**2)*MStop2(1))/
     &    MStop2(1)

	Cd(2283) = Cd(48)**3 + MGl2*(-3*Cd(48)*Cd(2254) + Cd(2282))

	Cd(2310) = MTy2*Cd(2309) + 3*(Cd(1496) - MStop2(1))*MStop2(1)

	Cd(2314) = 3*MStop2(1)**2 + MTy2*(3*MTy2 + 2*MStop2(1))

	Cd(2330) = Cd(3137) + 9*MStop2(1)

	Cd(2331) = 3*MTy2 - Cd(1595)

        Cd(2352) = (-Cd(519) + Cd(2351))*MStop2(1)**2 + 
     &   Cd(519)*MStop2(3)**2

        Cd(2362) = Cd(2360)*MStop2(1) + Cd(2361)*MStop2(2) + 
     &   16*MGl2*UStop2(1,1)

        Cd(2369) = Cd(2368)*MStop2(1)**2 + 
     &   3*MStop2(2)*(2*Cd(290)*MStop2(1) + MStop2(2))

        Cd(2372) = 2*MStop2(2)*
     &    (-(Cd(2361)*MStop2(1)) + MStop2(2)) - 
     &   5*MStop2(3)*MStop2(4)

        Cd(2374) = Cd(2373)*MStop2(1)**2 + 
     &   MStop2(2)*(-2*Cd(779)*MStop2(1) + MStop2(2))

        Cd(2379) = Cd(2378)*MStop2(1)**2 + 
     &   MStop2(2)*(6*Cd(2361)*MStop2(1) + MStop2(2))

        Cd(2389) = Cd(837)*MStop2(1)**2 - 
     &   3*MStop2(2)*(2*Cd(2361)*MStop2(1) + MStop2(2))

        Cd(2394) = -(Cd(2392)*MStop2(1)**2) + 
     &   MStop2(2)*(-2*Cd(2393)*MStop2(1) + Cd(2391)*MStop2(2))

	Cd(2398) = MGl2*(-2*MGl2 + Cd(2133)) + Cd(1592)*Cd(2137)

	Cd(2403) = 9*MTy2 - Cd(2965)

	Cd(2427) = M_3 - 4*M_3C*UCStop(1,3)**2 + M_3*UStop2(1,2)

        Cd(2447) = MGl2*(2*MGl2*Cd(2440) + Cd(2443)) + 
     &   Cd(2446)*MStop2(1)

	Cd(2462) = MGl2*(-2*MGl2 + Cd(2169)) + Cd(1559)*Cd(2170)

	Cd(2474) = Cd(2774) + 9*MStop2(2)

	Cd(2475) = 3*MTy2 - Cd(1074)

	Cd(23) = Cd(498) + At*UCStopC(1,3)

        Cd(2497) = M_3 - 4*M_3C*UCStop(1,3)**2 + 
     &   M_3*UStop2(1,2)**2

	Cd(246) = Cd(173) + Cd(45)*UStop2(1,1)

        Cd(2510) = MGl2*(2*MGl2*Cd(2505) + Cd(2509)) + 
     &   Cd(2506)*MStop2(2)

	Cd(2512) = 4*MTy2*Cd(1496) + 6*MGlpT2*MStop2(1) - 3*MStop2(1)**2

	Cd(2513) = 7*MGl2 + 21*MTy2 + 6*MStop2(1)

	Cd(2519) = MTy2*Cd(2518) + 3*Cd(45)*MStop2(2) + 3*MStop2(2)**2

	Cd(2520) = Cd(1022)*MStop2(4) + MTy2*(4*MTy2 + 4*MStop2(4))

	Cd(2523) = MGl2*(MGl2 - Cd(98)) + Cd(2519)

        Cd(2551) = Cd(2351)*MStop2(2)**2 + 
     &   Cd(519)*MStop2(1)*(-MStop2(2) - MStop2(3))

	Cd(2558) = MStop2(1) - Cd(2557)*MStop2(2) - 16*MGl2*UStop2(1,1)

	Cd(2573) = (3 + Cd(2572))*MStop2(2)**2 - 3*MStop2(3)**2

        Cd(2585) = Cd(2391)*MStop2(1)**2 + 
     &   MStop2(2)*(2*Cd(2584)*MStop2(1) + Cd(2583)*MStop2(2))

	Cd(2598) = 9*MStop2(1)**2 + MTy2*(3*MTy2 + 2*MStop2(1))

	Cd(2599) = MTy2*(MTy2 - 2*MStop2(1)) - 7*MStop2(1)**2

	Cd(2604) = 7*MStop2(1)**2 + MTy2*(3*MTy2 + 2*MStop2(1))

	Cd(2612) = 2*MTy2**2 + MGl2*(9*MGl2 + 5*MTy2)

        Cd(258) = M_3*MTy*MStop2(3)*UCStopC(1,3) + 
     &   Cd(49)*(Cd(188) - MStop2(1)*UStop2(1,1))

	Cd(2626) = MGl2*(MGl2 - 3*MTy2) + 3*MTy2**2

	Cd(2634) = 11*MGl2 + 5*MTy2 - 7*MStop2(1)

        Cd(2635) = Cd(2634)*MStop2(1)**2 + 
     &   MGlmT2*(MGlmT2**2 - Cd(2600)*MStop2(1))

	Cd(2650) = 10*MTy2**2 + MGl2*(31*MGl2 + 23*MTy2)

	Cd(2661) = MGl2*(MGl2 - Cd(1684)) + Cd(2604)

	Cd(2663) = MGl2*(MGl2 - Cd(1684)) + Cd(2598)

        Cd(2664) = MGl2*(Cd(44)*Cd(2599) - MGl2*Cd(2663)) - 
     &   2*Cd(44)**2*MStop2(1)**2

	Cd(263) = Cd(49) + MStop2(3)

	Cd(2682) = 9*MStop2(2)**2 + MTy2*(3*MTy2 + 2*MStop2(2))

	Cd(2683) = MTy2*(MTy2 - 2*MStop2(2)) - 7*MStop2(2)**2

	Cd(2687) = 7*MStop2(2)**2 + MTy2*(3*MTy2 + 2*MStop2(2))

	Cd(2710) = 11*MGl2 + 5*MTy2 - 7*MStop2(2)

        Cd(2711) = Cd(2710)*MStop2(2)**2 + 
     &   MGlmT2*(MGlmT2**2 - Cd(2600)*MStop2(2))

	Cd(2737) = MGl2*(MGl2 - Cd(1845)) + Cd(2687)

	Cd(2739) = MGl2*(MGl2 - Cd(1845)) + Cd(2682)

        Cd(2740) = MGl2*(Cd(48)*Cd(2683) - MGl2*Cd(2739)) - 
     &   2*Cd(48)**2*MStop2(2)**2

	Cd(2759) = 3*(MTy2**2 + MStop2(1)**2) + MTy2*MStop2(4)

	Cd(2764) = 2*MStop2(1)**2 + MTy2*(MTy2 + MStop2(1))

	Cd(2765) = MTy2*(MTy2 - MStop2(1)) - 4*MStop2(1)**2

	Cd(2767) = 12*MTy2 + Cd(5884)

	Cd(2768) = MTy2*(12*MTy2 + 4*Cd(286)) + 3*Cd(1583)*MStop2(1)

	Cd(2773) = MTy2*Cd(786) + 9*(MTy2**2 + MStop2(1)**2)

	Cd(2776) = MGl2*(4*MGl2 - Cd(2767)) + Cd(2768)

        Cd(26) = -(Conjugate(dAf133eps(0))*Cr(139)) + 
     &   1/2.D0*(Cd(25)*Conjugate(dAf133eps(-1))) - 
     &   2*Conjugate(dAf133eps(1))*MStop2(3)

	Cd(2785) = 2*MStop2(1)**2 + MTy2*(2*MTy2 + MStop2(1))

	Cd(2786) = 4*MStop2(1)**2 + MTy2*(MTy2 + 2*MStop2(1))

	Cd(2788) = Cd(4243) - 4*MTy2*MStop2(1)

	Cd(2789) = 10*MStop2(1)**3 + MTy2*(MTy2**2 - 3*MStop2(1)**2)

        Cd(2805) = (3*MGl2**2*Cd(594) + 
     &     (10*MGl2 + 3*Cd(594))*MStop2(1)**2)*MStop2(3)

	Cd(2807) = MGlmT2*(3*MGl2 + MTy2)

	Cd(2812) = Cd(1659) - MStop2(3)

        Cd(2821) = MTy2**2*Cd(607) + 
     &   2*MStop2(1)*(MTy2 + 2*MStop2(1))*MStop2(4)

	Cd(2826) = Cd(2617) - 6*MStop2(1)

	Cd(2830) = 5*Cd(2626) + MStop2(1)*(-Cd(2829) + 11*MStop2(1))

	Cd(2853) = MTy2**2 - 2*(MGlmT2 - 2*MStop2(1))*MStop2(1)

	Cd(2854) = MGl2*(MGl2 - 3*MTy2) + 3*Cd(2853)

	Cd(2858) = 3*MGl2 - Cd(2834)

	Cd(2861) = 5*MGl2 - Cd(6824)

        Cd(2870) = MGl2*MStop2(1)**2 + 
     &   MGlmT2*(MGlmT2**2 - Cd(153)*MStop2(1))

	Cd(2871) = 7*MGl2 + Cd(1659)

        Cd(2872) = Cd(2871)*MStop2(1)**2 + 
     &   MGlmT2*(MGlmT2**2 - Cd(2321)*MStop2(1))

        Cd(27) = -(Cr(139)*dAf133eps(0)) + 
     &   1/2.D0*(Cd(25)*dAf133eps(-1)) - 2*dAf133eps(1)*MStop2(3)

	Cd(2879) = MGlmT2*MTy2 + MStop2(1)*(-Cd(153) + MStop2(1))

	Cd(2881) = 8*MTy2*MGlpTmSt2(1) - 3*MGl2*MStop2(1) + MStop2(1)**2

	Cd(2889) = MGl2 - Cd(2922)

	Cd(2896) = Cd(2617) - 2*MStop2(1)

	Cd(285) = Cd(1559) - MStop2(3)

	Cd(2902) = MTy2*Cd(607) + 6*(MTy2**2 + MStop2(1)**2)

	Cd(2905) = 3*MGl2 - Cd(4928)

	Cd(2911) = 8*MStop2(1)**2 + MTy2*(3*MTy2 + 4*MStop2(1))

	Cd(287) = MTy2*(MTy2 + 2*Cr(1)) - Cd(286)*MStop2(2)

	Cd(2921) = 6*MTy2 + Cd(2220)

	Cd(288) = 2*MTy2*Cd(283) + Cd(523)

        Cd(2944) = MTy2**2*Cd(2943) + 
     &   MStop2(1)**2*(-4*MStop2(1) + 6*MStop2(3))

	Cd(2946) = Cd(4243) + MTy2*MStop2(1)

        Cd(2954) = MTy2*(4*MTy2**2 - 4*MStop2(1)*MStop2(2)) - 
     &   Cd(2953)*MStop2(3)

        Cd(291) = -(Cd(172)*MStop2(1)) - Cd(290)*MStop2(2) + 
     &   Cd(285)*UStop2(1,1)

	Cd(2968) = MGlpT2*MTy2 + (3*Cd(153) - MStop2(1))*MStop2(1)

        Cd(292) = Cd(172)*MStop2(1) + Cd(290)*MStop2(2) + 
     &   MGl2*UStop2(1,1)

	Cd(2978) = 2*MTy2 - Cd(2114)

	Cd(28) = 4*CA*MTy - Cd(5)*UCStop(1,3) - Cd(6)*UCStopC(1,3)

	Cd(2986) = MTy2*(6*MTy2 - Cd(2985)) + MStop2(1)*MStop2(3)

	Cd(3005) = 7*MStop2(1)**2 + MTy2*(3*MTy2 + 4*MStop2(1))

	Cd(3010) = 13*MGl2 + 5*MTy2 - 5*MStop2(1)

        Cd(3011) = Cd(3010)*MStop2(1)**2 + 
     &   MGlmT2*(3*MGlmT2**2 - Cd(3009)*MStop2(1))

        Cd(3020) = 6*MGl2**2*Cd(607) + 
     &   MStop2(1)*(-(MGl2*Cd(3018)) + 3*Cd(3019)*MStop2(1))

        Cd(3023) = 18*MGl2**2*Cd(607) + 
     &   MStop2(1)*(3*MGl2*Cd(3021) + Cd(3022)*MStop2(1))

        Cd(3025) = MStop2(3)*
     &   (MTy2*(12*MTy2 + Cd(1148)) + 9*MStop2(1)*MStop2(3))

	Cd(298) = 3*MTy2 + Cd(118)

	Cd(3038) = Cd(121) + Cd(164)

	Cd(299) = MGl2*MGlmT2 + 3*MGlmT2*MTy2 + Cr(86)*MStop2(1)

	Cd(3040) = 6*MGl2 + Cd(4081)

        Cd(3043) = 16*MTy2*MGlpTmSt2(1) - 9*MGl2*MStop2(1) + 
     &   3*MStop2(1)**2

        Cd(3045) = -(MGlpTmSt2(1)*MStop2(1)**2) + 
     &   MGlmT2*(MGlmT2**2 - MGlpT2*MStop2(1))

	Cd(3049) = -6*MStop2(1)**2 + MTy2*(5*MTy2 + 17*MStop2(1))

	Cd(3070) = 6*MTy2 + Cd(1583)

	Cd(29) = 4*CA*MTy + Cd(5)*UCStop(1,3) + Cd(6)*UCStopC(1,3)

	Cd(3080) = -3*MStop2(1)**2 + MTy2*(3*MTy2 + 2*MStop2(3))

	Cd(3081) = 2*MTy2*Cd(286) + 9*(MTy2**2 + MStop2(1)**2)

	Cd(3084) = Cd(2314) + MGl2*(MGl2 - Cd(2973))

	Cd(3113) = 3*MStop2(2)**2 + MTy2*(3*MTy2 + 2*MStop2(2))

	Cd(3117) = -MStop2(2)**2 + MTy2*(MTy2 - MStop2(3))

	Cd(3118) = 3*(MTy2**2 + MStop2(2)**2) + MTy2*MStop2(4)

	Cd(3123) = 2*MStop2(2)**2 + MTy2*(MTy2 + MStop2(2))

	Cd(3124) = MTy2*(MTy2 - MStop2(2)) - 4*MStop2(2)**2

	Cd(3126) = 12*MTy2 + Cd(1749)

	Cd(3128) = 4*MTy2*Cd(2200) + MStop2(2)*(Cd(3127) + 21*MStop2(2))

	Cd(3135) = 21*MGl2 + Cd(3130)

	Cd(3140) = MTy2*Cd(3139) + MStop2(2)*(5*MTy2 + 9*MStop2(2))

	Cd(3141) = MGl2*(4*MGl2 - Cd(3126)) + Cd(3128)

        Cd(3168) = (3*MGl2**2*Cd(874) + 
     &     (8*MGl2 + 3*Cd(874))*MStop2(2)**2)*MStop2(3)

	Cd(3173) = 15*MGl2 + Cd(1533)

	Cd(313) = 4*MGl2 + Cd(607)

	Cd(3184) = 4*MStop2(2)**2 + MTy2*(MTy2 + 2*MStop2(2))

	Cd(3185) = 10*MStop2(2)**3 + MTy2*(MTy2**2 - 3*MStop2(2)**2)

	Cd(3205) = 5*MGl2 - 3*(MTy2 + 2*MStop2(1))

	Cd(3215) = -3*MTy2 + Cd(118)

        Cd(3216) = Cd(2617)*MStop2(1) + 
     &   (5*Cd(3215) - 13*MStop2(2))*MStop2(2)

	Cd(3223) = MTy2*Cd(3221) + MStop2(2)*(Cd(3222) + 42*MStop2(2))

	Cd(3224) = 24*MTy2 + 32*MStop2(2) + MStop2(3)

	Cd(318) = 3*MTy2 + Cd(3378)

	Cd(3238) = 36*MGl2 + 23*MTy2 + MStop2(1)

	Cd(319) = Cd(1636) + MStop2(3)

	Cd(3242) = 2*MTy2*Cd(535) + 9*(MTy2**2 + MStop2(2)**2)

	Cd(3243) = -3*MStop2(2)**2 + MTy2*(3*MTy2 - 2*MStop2(3))

	Cd(3253) = Cd(2875) + 8*MStop2(2)

        Cd(3263) = 2*MTy2**2*Cd(874) + 
     &   MStop2(2)**2*(-8*MStop2(2) + 6*MStop2(3))

	Cd(3264) = MTy2*Cd(2518) + 3*MStop2(2)*(MTy2 + 2*MStop2(2))

	Cd(3267) = 3*MGl2 - Cd(2275)

        Cd(3268) = Cd(2617)*MStop2(1) + 
     &   (Cd(3267) - 7*MStop2(2))*MStop2(2)

        Cd(3272) = 3*MTy2*MStop2(1) - 
     &   MStop2(2)*(3*Cd(2086) + 14*MStop2(2))

	Cd(3285) = 6*MTy2 + Cd(2273)

	Cd(3287) = MTy2*(6*MTy2 + 2*Cd(1074)) + Cd(3286)*MStop2(2)

	Cd(3291) = MTy2*Cd(1541) + MStop2(2)*(Cd(3290) + 30*MStop2(2))

        Cd(3317) = MTy2*(4*MTy2**2 - 4*MStop2(1)*MStop2(2)) + 
     &   Cd(3316)*MStop2(3)

	Cd(3321) = MTy2*Cd(1541) + MStop2(2)*(Cd(1738) + 6*MStop2(2))

	Cd(3328) = 2*MGl2 - 9*MTy2 - 15*MStop2(1)

	Cd(3331) = MTy2*(6*MTy2 + 2*Cd(286)) + 3*Cd(874)*MStop2(2)

        Cd(3336) = MStop2(2)**2*(4*MStop2(2) - 6*MStop2(3)) - 
     &   2*MTy2**2*MStop2(4)

	Cd(3359) = MTy2*Cd(3357) + MStop2(2)*(Cd(3358) + 19*MStop2(2))

        Cd(3360) = MGl2**2*Cd(3233) + MTy2*Cd(3359) + 
     &   6*MStop2(2)**2*MStop2(3)

	Cd(3367) = 7*MStop2(2)**2 + MTy2*(3*MTy2 + 4*MStop2(2))

        Cd(3381) = Cd(3380)*MStop2(2)**2 + 
     &   MGl2*(3*Cd(2273)*MStop2(2) + 18*MGl2*MStop2(4))

        Cd(3384) = 3*Cd(3382)*MStop2(2)**2 + 
     &   MGl2*(-(Cd(3383)*MStop2(2)) + 6*MGl2*MStop2(4))

        Cd(3388) = MStop2(3)*
     &   (MTy2*(12*MTy2 + Cd(3318)) - 9*MStop2(2)*MStop2(3))

	Cd(3403) = 18*MGl2 + Cd(4530)

	Cd(3410) = -6*MStop2(2)**2 + MTy2*(5*MTy2 + 17*MStop2(2))

	Cd(3429) = 6*MGl2 + 13*MTy2 + MStop2(1)

	Cd(3437) = 42*MGl2 + 21*MTy2 + MStop2(1)

	Cd(3445) = 2*MTy2*Cd(2200) + MStop2(2)*(Cd(3252) + 7*MStop2(2))

	Cd(340) = Cd(45) - M_3*MTy*UCStopC(1,3)

	Cd(3453) = 9*MGl2 - 8*MTy2 + 5*MStop2(1)

	Cd(3458) = Cd(3113) + MGl2*(MGl2 - Cd(3355))

        Cd(3492) = Cd(3491)*MStop2(1) - Cd(2834)*MStop2(2) - 
     &   MStop2(2)**2

	Cd(3493) = 4*MStop2(1)**2 + MTy2*(2*MTy2 + MStop2(1))

	Cd(3494) = 12*MStop2(1)**2 + MTy2*(9*MTy2 + MStop2(1))

        Cd(3495) = 3*Cd(3493)*MStop2(1) - 
     &   MStop2(2)*(Cd(3494) + Cd(2973)*MStop2(2))

        Cd(3508) = -(Cd(3507)*MStop2(1)) + 
     &   MGl2*(MGl2*Cd(3506) - 4*MStop2(3)**2)

	Cd(3513) = MGl2*(3*MGl2*Cd(802) + Cd(3512)) + Cd(3510)*MStop2(1)

	Cd(3524) = 9*MStop2(1)**2 + MTy2*(6*MTy2 + 2*MStop2(1))

        Cd(3525) = 9*MStop2(1)**3 - 
     &   MStop2(2)*(2*Cd(3524) + Cd(2973)*MStop2(2))

	Cd(3529) = 2*MGl2 - Cd(1864)

	Cd(3530) = 2*MGl2 - Cd(3288)

        Cd(3531) = Cd(3529)*MStop2(1) + 
     &   MStop2(2)*(-2*Cd(3530) + MStop2(2))

        Cd(3537) = Cd(1864)*MStop2(1) - 
     &   MStop2(2)*(2*Cd(3288) + MStop2(2))

	Cd(3538) = MGl2*(MGl2 - 3*Cd(45)) + Cd(2759)

        Cd(3547) = -(Cd(44)*Cd(3546)*MStop2(1)) + 
     &   MGl2*(Cd(202) - 2*MGl2*MStop2(2))

	Cd(3563) = 2*MTy2*Cd(2617) + (Cd(3562) - 3*MStop2(1))*MStop2(1)

	Cd(3565) = MGl2 - Cd(160)

        Cd(3571) = MTy2*Cd(3570) + 
     &   3*MStop2(4)*(2*MTy2**2 + MStop2(1)*MStop2(4))

	Cd(3572) = MTy2**2*(6*MTy2 - MStop2(1)) + 3*MStop2(1)**3

	Cd(3573) = Cd(4243) - 6*MTy2*MStop2(1)

        Cd(3574) = Cd(3572)*MStop2(1) + 
     &   MStop2(2)*(2*Cd(45)*Cd(3573) + Cd(2314)*MStop2(2))

        Cd(3579) = Cd(2092)*MStop2(1) + 
     &   (2*Cd(1684) - MStop2(2))*MStop2(2)

	Cd(3580) = MTy2*(6*MTy2 - 2*MStop2(1)) - 3*MStop2(1)**2

	Cd(3581) = 9*MStop2(1)**2 + MTy2*(3*MTy2 + 7*MStop2(1))

        Cd(3582) = Cd(3580)*MStop2(1) + 
     &   MStop2(2)*(2*Cd(3581) - Cd(2973)*MStop2(2))

	Cd(3593) = 9*MStop2(1)**2 + MTy2*(3*MTy2 + 5*MStop2(1))

        Cd(3594) = 6*Cd(1866)*MStop2(1) + 4*Cd(3593)*MStop2(2) - 
     &   2*Cd(160)*MStop2(2)**2

        Cd(3609) = Cd(2092)*MStop2(1) + 
     &   (2*Cd(1789) - MStop2(2))*MStop2(2)

        Cd(3614) = 2*(Cd(529)*MStop2(1) + 
     &     MStop2(2)*(6*Cd(45) + MStop2(2)))

	Cd(3637) = -4*MTy2 + Cd(1954)

        Cd(3664) = Cd(1684)*MStop2(1) - Cd(2834)*MStop2(2) - 
     &   MStop2(2)**2

        Cd(3666) = Cd(3665)*MStop2(1)**2 + 
     &   MStop2(2)*(Cd(3494) + Cd(2973)*MStop2(2))

	Cd(3675) = -MStop2(1)**2 + MTy2*(MTy2 - MStop2(3))

        Cd(3687) = Cd(3686)*MStop2(1)**2 - 
     &   MStop2(2)*(2*Cd(3685)*MStop2(1) + 
     &      11*MStop2(2)*UUStop(1,1))

        Cd(3690) = Cd(3689)*MStop2(1)**2 + 
     &   MStop2(2)*(-2*Cd(3688)*MStop2(1) + 
     &      21*MStop2(2)*UUStop(1,1))

        Cd(3691) = Cd(3687)*MStop2(1) + 
     &   MGl2*(Cd(3690) - 36*MGl2*Cr(1)*UUStop(1,1))

        Cd(3694) = Cd(3692)*MStop2(1)**2 + 
     &   MStop2(2)*(-(Cd(3693)*MStop2(1)) + 
     &      15*MStop2(2)*UUStop(1,1))

        Cd(3697) = Cd(3695)*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(3696)*MStop2(1) + 29*MStop2(2)*UUStop(1,1))

        Cd(3700) = -(Cd(3698)*MStop2(1)**2) + 
     &   MStop2(2)*(Cd(3699)*MStop2(1) + 23*MStop2(2)*UUStop(1,1))

        Cd(3704) = -(Cd(3703)*MStop2(1)**2) + 
     &   MStop2(2)*(2*Cd(3702)*MStop2(1) + 7*MStop2(2)*UUStop(1,1))

        Cd(3707) = Cd(3706)*MStop2(1)**2 + 
     &   MStop2(2)*(2*Cd(3705)*MStop2(1) - 
     &      13*MStop2(2)*UUStop(1,1))

        Cd(3708) = Cd(3707)*MStop2(1) + 
     &   MGl2*(Cd(3704) - 6*MGl2*Cr(2)*UUStop(1,1))

	Cd(3720) = -6*MStop2(1)**2 + MTy2*(5*MTy2 + 9*MStop2(1))

	Cd(3736) = MTy2*(-6*MTy2 + Cd(1226)) - MStop2(1)*MStop2(3)

	Cd(3758) = -3*MStop2(1)**2 + MTy2*(9*MTy2 + 2*MStop2(1))

	Cd(3760) = MTy2*(9*MTy2 - 10*MStop2(1)) + 3*MStop2(1)**2

	Cd(3766) = MGl2*(MGl2 - 3*Cd(49)) + Cd(3113)

        Cd(3771) = 6*MTy2*MStop2(1) - 
     &   2*MStop2(2)*(Cd(3770) + 2*MStop2(2))

        Cd(3773) = 4*Cr(26)*MStop2(2)**2 - 
     &   2*MTy2*(MStop2(2)*MStop2(3) + MTy2*MStop2(4))

        Cd(36) = 8*CA*MTy*MTy2*Cd(35) - 
     &   Cd(29)*Cr(13)*MStop2(1) - Cd(28)*Cr(34)*MStop2(2)

        Cd(3786) = -(Cd(3785)*MStop2(2)) + 
     &   MGl2*(MGl2*Cd(3784) + 2*MStop2(3)**2)

	Cd(3791) = MGl2*(3*MGl2*Cd(802) + Cd(3788)) + Cd(3790)*MStop2(2)

	Cd(3796) = -3*MStop2(1)**2 + MTy2*(12*MTy2 + 8*MStop2(1))

	Cd(3798) = MTy2*(11*MTy2 - 4*MStop2(1)) + 3*MStop2(1)**2

	Cd(3800) = -4*MTy2 + Cd(118)

	Cd(3801) = -Cd(2031) + 2*Cd(3800)*MStop2(2)

	Cd(3804) = -3*MStop2(2)**2 + MTy2*(3*MTy2 - MStop2(3))

	Cd(373) = 4*MTy*SA - Cd(368)*UCStop(1,3) - Cd(369)*UCStopC(1,3)

	Cd(3806) = Cd(2031) + 2*Cd(67)*MStop2(2)

	Cd(3807) = MGl2*(MGl2 - 3*Cd(49)) + Cd(3118)

        Cd(3815) = MGl2*(Cd(202) - 2*MGl2*MStop2(1)) - 
     &   Cd(48)*Cd(3814)*MStop2(2)

	Cd(374) = 4*MTy*SA + Cd(368)*UCStop(1,3) + Cd(369)*UCStopC(1,3)

	Cd(3822) = MTy2*(3*MTy2 - 3*MStop2(1)) + 2*MStop2(1)**2

	Cd(3823) = -3*MStop2(1)**2 + MTy2*(3*MTy2 + 14*MStop2(1))

	Cd(3825) = 2*(MGl2**2 + MTy2*Cd(2200)) - MGl2*Cd(2518)

        Cd(3826) = 2*MTy2*Cd(2617) + 
     &   (2*Cd(3234) - 3*MStop2(1))*MStop2(1)

	Cd(3827) = Cd(3562) - 6*MStop2(1)

        Cd(3831) = MTy2*Cd(3830) + 
     &   3*MStop2(4)*(2*MTy2**2 + MStop2(2)*MStop2(4))

	Cd(3832) = MTy2*(3*MTy2 - 5*MStop2(1)) + MStop2(1)**2

	Cd(3833) = Cd(3006) + 10*MTy2*MStop2(1)

	Cd(3835) = 2*MGl2 - Cd(2518)

        Cd(375) = -(A0delStop(1)*Cd(373)) - 
     &   A0delStop(2)*Cd(374) - 8*MTy*SA*Cd(12862)

        Cd(3836) = Cd(3835)*MStop2(1) - 
     &   MStop2(2)*(2*Cd(1592) + MStop2(2))

        Cd(3840) = 4*MStop2(2)*
     &   (Cd(2252)*MStop2(1) + (2*Cd(2987) - MStop2(2))*MStop2(2))

	Cd(3841) = -3*MStop2(1)**2 + MTy2*(6*MTy2 + 14*MStop2(1))

	Cd(3844) = Cd(7559) + MTy2*MStop2(1)

	Cd(3845) = 3*MStop2(1)**2 + MTy2*(7*MTy2 + 2*MStop2(1))

        Cd(3851) = Cd(2252)*MStop2(1) + 
     &   (2*Cd(2987) - MStop2(2))*MStop2(2)

        Cd(3853) = 6*Cd(3852)*MStop2(2) + 
     &   4*MTy2*(Cd(1141)*MStop2(1) + 3*MTy2*MStop2(4))

	Cd(3854) = (MTy2 - MStop2(1))*(3*MTy2 + 2*MStop2(1))

	Cd(3855) = MTy2*(5*MTy2 - 2*MStop2(1)) + 3*MStop2(1)**2

	Cd(3859) = MTy2*Cd(67) + (Cd(3858) - 14*MStop2(2))*MStop2(2)

        Cd(378) = 8*MTy*MTy2*SA*Cd(35) - 
     &   Cd(373)*Cr(13)*MStop2(1) - Cd(374)*Cr(34)*MStop2(2)

        Cd(3876) = Cd(2518)*MStop2(1) + 
     &   MStop2(2)*(2*Cd(1592) + MStop2(2))

	Cd(3887) = 2*MTy2*Cd(2617) + MStop2(1)*(10*MGlmT2 + 3*MStop2(1))

        Cd(3916) = 2*MStop2(2)*(Cd(2200) + 2*MStop2(2)) + 
     &   MStop2(3)*MStop2(4)

	Cd(3917) = -3*MStop2(1)**2 + MTy2*(9*MTy2 + 8*MStop2(1))

	Cd(3919) = Cd(7559) + 3*MTy2*MStop2(1)

	Cd(3920) = MTy2*(9*MTy2 - 4*MStop2(1)) + 3*MStop2(1)**2

	Cd(3932) = MTy2*Cd(315) + 3*(MTy2**2 + MStop2(2)**2)

        Cd(3944) = 2*Cd(3942)*MStop2(1)*MStop2(2) - 
     &   Cd(3943)*MStop2(2)**2 + 11*MStop2(1)**2*UUStop(1,1)

        Cd(3947) = 2*Cd(3945)*MStop2(1)*MStop2(2) - 
     &   Cd(3946)*MStop2(2)**2 - 21*MStop2(1)**2*UUStop(1,1)

        Cd(3948) = Cd(3944)*MStop2(2) + 
     &   MGl2*(Cd(3947) - 36*MGl2*Cr(1)*UUStop(1,1))

        Cd(3951) = -(Cd(3950)*MStop2(1)*MStop2(2)) + 
     &   Cd(3949)*MStop2(2)**2 + 15*MStop2(1)**2*UUStop(1,1)

        Cd(3954) = -(MStop2(2)*
     &      (Cd(3953)*MStop2(1) + Cd(3952)*MStop2(2))) + 
     &   29*MStop2(1)**2*UUStop(1,1)

        Cd(3957) = Cd(3955)*MStop2(1)*MStop2(2) + 
     &   Cd(3956)*MStop2(2)**2 - 23*MStop2(1)**2*UUStop(1,1)

        Cd(3960) = 2*Cd(3959)*MStop2(1)*MStop2(2) - 
     &   Cd(3703)*MStop2(2)**2 + 7*MStop2(1)**2*UUStop(1,1)

        Cd(3962) = 2*Cd(3961)*MStop2(1)*MStop2(2) + 
     &   Cd(3706)*MStop2(2)**2 - 13*MStop2(1)**2*UUStop(1,1)

        Cd(3963) = Cd(3962)*MStop2(2) + 
     &   MGl2*(Cd(3960) + 6*MGl2*Cd(652)*UUStop(1,1))

	Cd(3971) = -6*MStop2(2)**2 + MTy2*(5*MTy2 + 9*MStop2(2))

	Cd(38) = Cd(11105) + UStop2(1,2)

	Cd(3989) = MTy2*(-6*MTy2 + Cd(786)) + MStop2(2)*MStop2(3)

	Cd(3993) = MTy2*Cd(3357) + MStop2(2)*(Cd(3992) + 27*MStop2(2))

	Cd(4018) = 2*MStop2(1)**2 + MTy2*(3*MTy2 + 3*MStop2(1))

	Cd(4028) = MTy2*Cd(2617) - MStop2(1)*(3*MTy2 + 2*MStop2(1))

	Cd(4046) = MTy2*(MTy2 - 4*MStop2(1)) - MStop2(1)**2

	Cd(4049) = 3*Cd(44)**2 + MGl2*(MGl2 - 3*Cd(45))

	Cd(4054) = Cd(44)**3 + MGl2*(-3*Cd(44)**2 + MGl2*Cd(1789))

	end
