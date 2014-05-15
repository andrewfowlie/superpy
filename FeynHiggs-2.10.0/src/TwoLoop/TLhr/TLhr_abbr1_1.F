* TLhr_abbr1_1.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_1
	implicit none

#include "TLhr.h"


	Opt(3) = MGlmT2*(-3*MTy2**2 + MGl2*(MGl2 + 6*MTy2))

	Opt(4) = MGl2*(2*MGl2 - 8*MTy2) + 12*MTy2**2 + 7*MStop2(2)**2

        Opt(6) = MTy*(MUE*UCStop(1,3)*
     &      (-2*M_3C*UCStop(1,3)**2 - 
     &        2*M_3*UStop2(1,1)*UStop2(1,2)) + 
     &     MUEC*(-2*M_3C*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) - 
     &        2*M_3*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)))

        Opt(7) = (-5*MTy2**2 + MGl2*(-3*MGl2 + 7*MTy2))*
     &    MStop2(1) - 11*MStop2(2)**3

        Opt(13) = MStop2(1)*
     &   (4*MTy2*MStop2(1) - 11*MStop2(1)**2 + 
     &     MGl2*(-MGl2 + 7*MStop2(1)))

        Opt(15) = MTy*(-2*M_3C*MUE**2*UCStop(1,3)**3 - 
     &     2*M_3*MUEC**2*UCStopC(1,3)**3 - 
     &     2*(M_3*MUE**2 + 2*M_3C*MUE2)*UCStop(1,3)*
     &      UStop2(1,1)*UStop2(1,2))

        Opt(16) = UStop2(1,1)*
     &   (-(At*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))*
     &        (UStop2(1,1) - UStop2(1,2))) + 
     &     MGl2*UCStop(1,3)*UStop2(1,2))

        Opt(17) = 2*(2*MTy2**3 + 
     &     MGl2*(MGl2*(MGl2 - 3*MTy2) + 
     &        MTy2*(2*MTy2 + 2*MStop2(1))))

        Opt(19) = At*(MUEC*UCStopC(1,3)*
     &       (-(M_3*UCStopC(1,3)**2) - 
     &         M_3C*UStop2(1,1)*UStop2(1,2)) + 
     &      MUE*(-(M_3C*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)) - 
     &         M_3*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2))) + 
     &   AtC*(MUE*UCStop(1,3)*
     &       (-(M_3C*UCStop(1,3)**2) - 
     &         M_3*UStop2(1,1)*UStop2(1,2)) + 
     &      MUEC*(-(M_3C*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)) - 
     &         M_3*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)))

	Opt(20) = Opt(19) + MGlpT2*(MUE*UCStop(1,3) + MUEC*UCStopC(1,3))

	Cd(1) = dMSfsq1233eps(0) - 2*dMSfsq1233eps(1)

	Cd(9948) = 1 - UStop2(1,1)

	Cd(9950) = 4*MGl2*MTy2 - MGlpTmSt4(1)

	Cd(987) = 6*MStop2(3) - 2*MStop2(4)*UStop2(1,2)**2

	Cd(9955) = Cd(9948) - UStop2(1,2)

        Cd(9956) = Cr(28)*T134(MGl2,MTy2,MStop2(1)) + 
     &   T134(MGl2,MTy2,MStop2(2))*UStop2(1,2)**2

	Cd(9957) = AtC*UCStop(1,3)**2 + At*UStop2(1,1)*UStop2(1,2)

	Cd(9959) = 4*UStop2(1,1)**2 + UStop2(1,2)**2

	Cd(9960) = At*Cd(9959) + 2*M_3*UStop2(1,1)

	Cd(9962) = M_3 + At*UStop2(1,1)

	Cd(9963) = M_3C + AtC

	Cd(9967) = 2*M_3C + AtC

	Cd(9969) = Cd(9948) + UStop2(1,1)

	Cd(9975) = M_3*Cr(28) - Cd(9967)*UCStop(1,3)**2

	Cd(9976) = 2*Cd(9975) + At*Cd(9969)*Cr(28)

        Cd(9979) = Cr(28)*T134(MGl2,MTy2,MStop2(2)) + 
     &   T134(MGl2,MTy2,MStop2(1))*UStop2(1,2)**2

	Cd(9981) = Cd(9948)*UStop2(1,1) + UStop2(1,2)**2

	Cd(9982) = AtC*UCStop(1,3)**2 + At*Cd(9981)*UStop2(1,1)

	Cd(9983) = 3*UStop2(1,1)**2 + UStop2(1,2)**2

	Cd(990) = 4*MStop2(3) + MStop2(4)

	Cd(9988) = At*UCStopC(1,3)**2 + AtC*Cd(9981)*UStop2(1,1)

	Cd(9993) = At*UCStopC(1,3)**2 + AtC*UStop2(1,1)*UStop2(1,2)

	Cd(9994) = M_3C + AtC*UStop2(1,1)

	Cd(991) = 6*MStop2(2)**2 - Cd(990)*MStop2(3)

	Cd(9995) = M_3 + At

	Cd(9997) = M_3 - At*UStop2(1,2)

	Cd(10001) = 2*M_3 + At

	Cd(10005) = M_3C*Cr(28) - Cd(10001)*UCStopC(1,3)**2

	Cd(10006) = 2*Cd(10005) + AtC*Cd(9969)*Cr(28)

        Cd(10014) = T134(MGl2,MTy2,MStop2(1)) - 
     &   T134(MGl2,MTy2,MStop2(2))

	Cd(10015) = M_3C*At + M_3*AtC

	Cd(10016) = UStop2(1,1)**2 - Cr(75)*UStop2(1,2)

	Cd(10027) = 1 - 2*Cd(9948)

        Cd(10032) = UCStop(1,3)*
     &   (M_3C*UCStop(1,3)**2 + M_3*UStop2(1,1)*UStop2(1,2))

	Cd(10033) = AtC*Cd(10032) + At*Cr(5)*UStop2(1,1)*UStop2(1,2)

        Cd(10039) = M_3*UCStopC(1,3)**2 + 
     &   M_3C*UStop2(1,1)*UStop2(1,2)

	Cd(10054) = MUEC*Cr(28) - 2*MUE*UCStop(1,3)**2

	Cd(997) = 8*MStop2(2)**2 - 4*MStop2(3)**2

	Cd(10055) = MUE*Cr(28) - 2*MUEC*UCStopC(1,3)**2

	Cd(10056) = M_3C*Cd(10054) + M_3*Cd(10055)

	Cd(10058) = MUEC*Cr(28) + 2*MUE*UCStop(1,3)**2

	Cd(10059) = MUE*Cr(28) + 2*MUEC*UCStopC(1,3)**2

	Cd(10060) = M_3C*Cd(10058) + M_3*Cd(10059)

	Cd(10063) = MUE*Cd(9983) + 2*MUEC*UCStopC(1,3)**2

	Cd(998) = 12*MStop2(3)**2 - Cd(997)*UStop2(1,2)**2

	Cd(10070) = Cd(9948)*Cd(9983) + Cr(75)**2*UStop2(1,1)

	Cd(1000) = 12*MStop2(3) - 4*MStop2(4)*UStop2(1,2)**2

	Cd(10085) = M_3*Cd(10016) - 3*M_3C*UCStop(1,3)**2

        Cd(10086) = M_3C*UCStop(1,3)**2 + 
     &   M_3*UStop2(1,1)*UStop2(1,2)

	Cd(10092) = -(Cd(9948)*Cr(75)) + UStop2(1,1)**2

	Cd(10094) = M_3C*Cd(10016) - 3*M_3*UCStopC(1,3)**2

	Cd(10103) = 4*MGl2*MTy2 - MGlpTmSt4(2)

	Cd(10110) = 4*M_3C + AtC

	Cd(1003) = 3*MStop2(3) - MStop2(4)*UStop2(1,2)**2

	Cd(10135) = 2*M_3C*At + M_3*AtC

	Cd(10136) = Cd(10016)*Cd(10135) - 3*M_3C*AtC*UCStop(1,3)**2

	Cd(10137) = M_3C*At + 2*M_3*AtC

	Cd(10140) = -4*M_3 + At

	Cd(10141) = 4*M_3C*UCStop(1,3) - Cd(10140)*UCStopC(1,3)

	Cd(10142) = 2*M_3*UStop2(1,1) + At*UStop2(1,2)

	Cd(10146) = 3*M_3 + At

	Cd(10147) = 3*M_3C*UCStop(1,3) + Cd(10146)*UCStopC(1,3)

	Cd(10162) = 5 - 8*Cd(9948)

	Cd(10164) = 5 - 12*UStop2(1,2)

	Cd(10169) = M_3*Cr(28) - 2*M_3C*UCStop(1,3)**2

	Cd(10170) = -4*Cd(10169) + At*Cd(9969)*Cr(28)

	Cd(1009) = MStop2(3) - 2*MStop2(4)

	Cd(10176) = M_3C*Cr(28) - 2*M_3*UCStopC(1,3)**2

	Cd(10180) = At - 8*M_3*Cr(85)

	Cd(10181) = M_3C*Cd(9983) + 2*M_3*UCStopC(1,3)**2

	Cd(1010) = MStop2(2)**2 - Cd(1009)*MStop2(3)

	Cd(10189) = AtC*UCStop(1,3)**2 + 2*At*UStop2(1,1)*UStop2(1,2)

	Cd(10190) = AtC*Cd(10189) + At**2*UCStopC(1,3)**2

	Cd(1011) = 6*MStop2(2) - 17*MStop2(3)

        Cd(10200) = 4*MTy2*UCStop(1,3) + 
     &   At**2*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

	Cd(10209) = -8*M_3C + AtC

	Cd(10219) = MUEC*UCStopC(1,3)**2 + MUE*UStop2(1,1)*UStop2(1,2)

	Cd(10224) = M_3*MUE + M_3C*MUEC

        Cd(10233) = MUE**2*UCStop(1,3)**2 + 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2)

	Cd(10234) = 2*M_3*MUE2 + M_3C*MUEC**2

	Cd(99) = MGl2 - 2*MStop2(2)

	Cd(10243) = MTy + 2*M_3C*UCStop(1,3)

        Cd(10249) = MStop2(2)*UCStop(1,3) + 
     &   At*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(10251) = MTy + 2*M_3*UCStopC(1,3)

	Cd(10271) = MGl2 + MTy2*UStop2(1,2)

	Cd(10272) = MTy2*Cd(9948)*UStop2(1,1) + Cd(10271)*UStop2(1,2)

	Cd(10273) = M_3*MUE**2 + 2*M_3C*MUE2

	Cd(10277) = 4*MGlpTmSt2(1) + MStop2(1)

	Cd(10278) = MGlmT2**2 - Cd(10277)*MStop2(1)

	Cd(10298) = 1 - Cd(9948)

	Cd(10303) = MGlmT2*UStop2(1,1) + MStop2(1)*UStop2(1,2)

	Cd(1022) = 10*MStop2(2) + MStop2(3)

	Cd(10312) = -2*M_3 + At

	Cd(10320) = -(M_3C*At) + 2*MGlmT2

	Cd(10321) = M_3C*Cr(28) - 3*M_3*UCStopC(1,3)**2

	Cd(10323) = -(At*Cr(5)) + MGlmT2*UCStop(1,3)

	Cd(10329) = -3*M_3*UCStopC(1,3)**2 + M_3C*Cr(75)*UStop2(1,1)

	Cd(1025) = MStop2(2) - 12*MStop2(3)

	Cd(10340) = 13 - 12*UStop2(1,2)

	Cd(1026) = 2 - 16*UStop2(1,2)

	Cd(10347) = 4*Cd(10176) + AtC*Cd(9969)*Cr(28)

        Cd(10349) = M_3C*UCStop(1,3)**2 + 
     &   M_3*UStop2(1,1)*UStop2(1,2)**2

	Cd(10353) = M_3C*AtC*UCStop(1,3)**2 + M_3*At*UCStopC(1,3)**2

        Cd(10356) = Cd(10209)*Cr(85)*UCStop(1,3)**2 + 
     &   4*M_3*Cr(28)*UStop2(1,1)

	Cd(10359) = 4*M_3*Cr(28) - Cd(10209)*UCStop(1,3)**2

	Cd(10360) = M_3C*Cr(28) + 2*M_3*UCStopC(1,3)**2

	Cd(10361) = -4*Cd(10360) + AtC*Cr(28)

	Cd(10369) = At*Cr(28) - 2*AtC*UCStop(1,3)**2

        Cd(10372) = Cd(9963)*UCStop(1,3)**2 + 
     &   M_3*UStop2(1,1)*UStop2(1,2)

	Cd(10376) = AtC*Cr(28) - 2*At*UCStopC(1,3)**2

        Cd(10378) = Cd(9995)*UCStopC(1,3)**2 + 
     &   M_3C*UStop2(1,1)*UStop2(1,2)

	Cd(10382) = MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)**2

	Cd(10383) = MStop2(1) - MGlmT2*UStop2(1,1)

        Cd(10388) = MUEC*UCStopC(1,3)**2 + 
     &   MUE*UStop2(1,1)*UStop2(1,2)**2

	Cd(10398) = -4*M_3C + AtC

	Cd(10399) = -4*M_3*Cd(10055) + Cd(10054)*Cd(10398)

	Cd(1032) = MStop2(3) - 3*MStop2(4)

	Cd(10407) = MStop2(3) - 2*Cd(9955)*MStop2(1)*UStop2(1,2)

	Cd(1033) = 6*MStop2(1)**2 - 4*Cd(1032)*MStop2(3)

	Cd(10419) = M_3C*At + MStop2(1)

	Cd(1035) = 5*MStop2(2) + 6*MStop2(3)

	Cd(10437) = MGlmT2*Cr(28) + 2*MStop2(3)*UStop2(1,2)

	Cd(10443) = MGlmT2 - Cd(9969)*MStop2(1)

	Cd(10453) = MStop2(1)*UStop2(1,1) + MGlmT2*UStop2(1,2)

	Cd(10454) = M_3*Cr(28) - 3*M_3C*UCStop(1,3)**2

	Cd(10455) = -3*M_3C*UCStop(1,3)**2 + M_3*Cr(75)*UStop2(1,1)

	Cd(10456) = -(MGlmT2*Cd(10455)) + Cd(10454)*MStop2(1)

	Cd(10458) = Cr(28)*MStop2(1) - MGlmT2*Cr(75)*UStop2(1,1)

	Cd(10460) = -(MGlmT2*Cd(10329)) + Cd(10321)*MStop2(1)

	Cd(1040) = 11 + 8*UStop2(1,2)**2

	Cd(10486) = 4*M_3*MTy2 - At*MGlpTmSt2(1)

	Cd(10487) = 4*M_3C*MTy2*UCStop(1,3) + Cd(10486)*UCStopC(1,3)

	Cd(10492) = 3*MStop2(3) - 4*MTy2*Cd(9955)*UStop2(1,2)

	Cd(10496) = 4*MTy2*UCStop(1,3) + At**2*UCStopC(1,3)

	Cd(10500) = MStop2(1)*UStop2(1,1) + MGlpT2*UStop2(1,2)

	Cd(1042) = 8 - 3*UStop2(1,2)**2

        Cd(10506) = MStop2(1)*UCStop(1,3) + 
     &   At*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(1043) = MStop2(1) + MTy2*Cd(1042)*UStop2(1,2)

	Cd(10520) = MStop2(2) - 3*MStop2(3)

        Cd(10522) = Cr(28)*MStop2(1) - 
     &   2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(1045) = Cd(1040) - 4*UStop2(1,1)**2

	Cd(10538) = MStop2(1) - MGl2*UStop2(1,2)

	Cd(10589) = 2*MUE2*Cr(28) - 3*MUEC**2*UCStopC(1,3)**2

	Cd(10590) = MUE**2*Cr(28)*UCStop(1,3) + Cd(10589)*UCStopC(1,3)

	Cd(1051) = 2*MStop2(2)**2 - MStop2(1)*MStop2(3)

        Cd(10593) = -3*MUEC**2*UCStopC(1,3)**2 + 
     &   2*MUE2*Cr(75)*UStop2(1,1)

        Cd(10594) = Cd(10593)*UCStopC(1,3) + 
     &   MUE**2*Cr(75)*UCStop(1,3)*UStop2(1,1)

	Cd(10601) = 4*MGlpTmSt2(2) + MStop2(2)

	Cd(10602) = MGlmT2**2 - Cd(10601)*MStop2(2)

	Cd(1052) = 3 - 2*UStop2(1,2)

	Cd(10606) = AtC**2*UCStop(1,3)**2 + 2*MGlmT2*UStop2(1,1)

	Cd(10613) = M_3 + 2*At*UStop2(1,1)

	Cd(10614) = 2*Cd(9948)*Cd(10613) + At*UStop2(1,1)**2

	Cd(10615) = Cd(9948) + 3*UStop2(1,1)

	Cd(10617) = 2*M_3 + At*Cd(9955)

	Cd(1054) = 11 - 4*UStop2(1,2)**2

	Cd(10648) = 2 - Cd(9948)

	Cd(10650) = 8*M_3C + AtC*Cd(9969)

	Cd(10666) = Cr(28)*MStop2(2) + Cd(9981)*MStop2(3)*UStop2(1,1)

        Cd(10668) = -(Cd(10382)*MStop2(2)) + 
     &   MTy2*MUEC*UStop2(1,1)**2*UStop2(1,2)

	Cd(10669) = MGlpT2 + 2*M_3C*MTy*UCStop(1,3)

        Cd(10670) = Cd(10669)*UCStop(1,3) + 
     &   2*M_3*MTy*UStop2(1,1)*UStop2(1,2)

	Cd(10674) = -2*MUE*UCStop(1,3)**2 + MUEC*Cr(75)**2*UStop2(1,1)

	Cd(10678) = MGl2 - MStop2(3)

	Cd(10679) = 2*M_3*MTy*Cr(75) + Cd(10678)*UCStop(1,3)

	Cd(1060) = 7*MStop2(2) - MStop2(3)

        Cd(10687) = -(Cd(10388)*MStop2(2)) + 
     &   MTy2*MUE*UStop2(1,1)**2*UStop2(1,2)

	Cd(10688) = MGlpT2 + 2*M_3*MTy*UCStopC(1,3)

        Cd(10689) = Cd(10688)*UCStopC(1,3) + 
     &   2*M_3C*MTy*UStop2(1,1)*UStop2(1,2)

	Cd(1061) = 8*MStop2(2)**2 - Cd(1060)*MStop2(3)

        Cd(10696) = -2*MUEC*Cr(85)*UCStopC(1,3)**2 + 
     &   MUE*Cr(75)**2*UStop2(1,1)

	Cd(10702) = MStop2(3) - 2*MTy2*UStop2(1,1)

	Cd(1062) = 9*MStop2(2) - 5*MStop2(3)

        Cd(10703) = 2*MGlmT2*MTy2*Cd(9983) - 
     &   Cd(10702)*MStop2(2)*UStop2(1,1)

	Cd(10704) = MUEC*Cd(10703) + 4*MTy2*MUE*MStop2(2)*UCStop(1,3)**2

	Cd(10705) = MGlpT2 - MStop2(3)

	Cd(10707) = MUE*Cd(10703) + 4*MGlmT2*MTy2*MUEC*UCStopC(1,3)**2

	Cd(1063) = 6*MStop2(2)**2 - Cd(1062)*MStop2(3)

	Cd(10721) = 1 - 6*UStop2(1,2)

	Cd(1064) = 2*MStop2(2) + MStop2(3)

	Cd(10724) = 2*MTy2*Cd(10696) + Cd(10219)*MStop2(3)

	Cd(10729) = MGlmT2 + 3*MStop2(2)

	Cd(10741) = MGlmT2*MStop2(1) - MStop2(2)*MStop2(3)

	Cd(10742) = MStop2(2) - MGlmT2*UStop2(1,1)

	Cd(10745) = MStop2(1) - Cd(9955)*MStop2(3)

	Cd(10749) = MUE**2*UCStop(1,3)**2 + MUE2*UStop2(1,1)**2

	Cd(1067) = 2*MStop2(2) + MStop2(3)*UStop2(1,2)

	Cd(10754) = MUE2*Cd(9983) + MUEC**2*UCStopC(1,3)**2

	Cd(10755) = MGlmT2*Cd(10754) + Cd(10749)*MStop2(2)

	Cd(10759) = -3*M_3C*UCStop(1,3)**2 + M_3*UStop2(1,1)**2

	Cd(10764) = -3*M_3*UCStopC(1,3)**2 + M_3C*UStop2(1,1)**2

	Cd(10792) = -(Cd(10209)*UCStop(1,3)**2) + 4*M_3*UStop2(1,1)**2

	Cd(10802) = MStop2(3) + Cr(75)**2*MStop2(1)*UStop2(1,1)

	Cd(10805) = -(Cd(9959)*MStop2(3)) + MStop2(1)*UStop2(1,1)**2

	Cd(10806) = MGlmT2*Cd(9983)*MStop2(1) + Cd(10805)*MStop2(2)

	Cd(10807) = MStop2(1)*MStop2(2) - MGlmT2*MStop2(3)

	Cd(10816) = MStop2(1) + MStop2(3)*UStop2(1,2)

	Cd(10819) = -(At*Cr(85)*MStop2(1)) + M_3*MStop2(3)

        Cd(10820) = Cd(10819)*UCStopC(1,3)**2 + 
     &   M_3C*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(1074) = 3*MStop2(2) + MStop2(3)

	Cd(10824) = Cd(9983)*MStop2(2) + MGlmT2*UStop2(1,1)**2

	Cd(10835) = -3*MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)**2

	Cd(10841) = -3*MUEC*UCStopC(1,3)**2 + MUE*UStop2(1,1)**2

	Cd(10851) = 2*M_3C*UCStop(1,3)**2 + M_3*UStop2(1,1)**2

	Cd(10852) = M_3*MGlmT2*Cd(9983) + Cd(10851)*MStop2(2)

        Cd(10853) = Cd(10398)*UCStop(1,3)**2 + 
     &   2*M_3*Cr(75)**2*UStop2(1,1)

	Cd(10860) = At*UCStopC(1,3)**2 + 2*M_3C*UStop2(1,1)**2

	Cd(1078) = 15*MStop2(2) - 13*MStop2(3)

	Cd(10861) = 2*M_3C*Cd(9983) - Cd(10140)*UCStopC(1,3)**2

	Cd(10862) = MGlmT2*Cd(10861) + Cd(10860)*MStop2(2)

        Cd(10864) = Cd(10140)*Cr(85)*UCStopC(1,3)**2 + 
     &   2*M_3C*Cr(75)**2*UStop2(1,1)

	Cd(1079) = 4*(6*MStop2(1)**2 + Cd(1078)*MStop2(3))

	Cd(1080) = 18*MStop2(1)*MStop2(3) + 4*Cd(1051)*UStop2(1,2)

	Cd(10887) = MUE**2*UCStop(1,3) + 2*MUE2*UCStopC(1,3)

	Cd(1081) = 12*MStop2(2)**2 - MStop2(3)*MStop2(4)

	Cd(10904) = Cd(9948) + UStop2(1,2)

        Cd(10906) = MUEC*UUStop(1,3)*UUStop(2,2) - 
     &   MUE*UUStop(1,1)*UUStop(2,3)

        Cd(10907) = MUE*UUStop(3,1)*UUStop(3,2)*UUStopC(2,2) + 
     &   Cd(10906)*UUStopC(2,3)

        Cd(10908) = MUE*Cd(10907) - 
     &   MUE2*UUStop(1,1)*UUStop(2,2)*UUStopC(2,1)

        Cd(10909) = MUE*UUStop(1,3)*UUStop(2,1) - 
     &   MUEC*UUStop(1,2)*UUStop(2,3)

	Cd(10918) = 1 - UStop2(1,1)**2

	Cd(10919) = 2*MStop2(1) - 3*MStop2(4)*UStop2(1,1)

	Cd(10920) = 5*MStop2(3) + Cd(10919)*UStop2(1,1)

	Cd(10921) = 3*UUStop(3,3) + 4*UUStop(3,4)

	Cd(10922) = 12 - 5*UStop2(1,1)**2

        Cd(10923) = Cd(10922)*MStop2(3)*UStop2(1,1) - 
     &   2*MStop2(1)*UStop2(1,2)

        Cd(10924) = Cd(10923)*UStop(1,2) + 
     &   2*Cd(10921)*MStop2(1)*UCStop(3,1)*UStop2(1,1)*UStopC(1,1)

	Cd(10925) = 2*UUStop(3,3) + UUStop(3,4)

	Cd(10926) = UStop(1,2)*UStop2(1,1) - Cd(10925)*UStopC(2,1)

	Cd(10927) = Cr(7)*MStop2(1)**2 - 2*MStop2(3)**2*UStop2(1,1)

        Cd(10928) = 6*Cd(10926)*MStop2(1)*MStop2(3) - 
     &   4*Cd(10927)*UStop(1,2)

        Cd(10929) = 2*Cd(10926)*MStop2(1) + 
     &   5*MStop2(3)*UStop(1,2)*UStop2(1,1)

	Cd(10930) = Cr(7)*MStop2(1)**2 - MStop2(3)**2*UStop2(1,1)

	Cd(10931) = Cd(10926)*MStop2(1)*MStop2(3) - Cd(10930)*UStop(1,2)

        Cd(10932) = 6*MStop2(3)**2*UStop(1,1) - 
     &   2*Cr(91)*MStop2(1)**2*UStop2(1,1)*UStopC(2,2)

	Cd(10933) = Cr(91)*UCStop(3,1) - 6*UStop2(1,1)*UUStop(1,3)

	Cd(10935) = 4*MStop2(3)*UStop2(1,1)**3 + 5*MStop2(1)*UStop2(1,2)

        Cd(10938) = -(Cd(10932)*UStop(1,2)) + 
     &   2*Cd(10933)*MStop2(1)**2*UStop2(1,2)

	Cd(10939) = MStop2(3)*UStop2(1,1)**3 + 3*MStop2(1)*UStop2(1,2)

        Cd(10940) = Cd(10938) - 
     &   2*Cd(10921)*MStop2(1)*MStop2(3)*UCStop(3,1)*UStop2(1,1)

        Cd(10941) = 4*Cd(10939)*MStop2(3)*UStop(1,2) + 
     &   2*Cd(10940)*UStopC(1,1)

	Cd(10946) = MUE*UUStop(1,1) - MUEC*UUStop(1,2)

	Cd(10947) = 4*MUE2*UUStop(1,1) - MUEC**2*UUStop(1,2)

	Cd(10948) = MUE**2*UUStop(1,1)**2 - Cd(10947)*UUStop(1,2)

	Cd(10949) = 2*MUE**2*UUStop(1,1) - Cd(10948)*UUStopC(1,1)

	Cd(10950) = -UStop(2,2) + UCStopC(3,1)*UStop(1,2)*UStop2(1,1)

	Cd(10951) = Cd(9948)*Cr(91) - UStop2(1,1)*UUStop(3,3)

        Cd(10952) = -5*UStop(2,2) + 
     &   3*UCStopC(3,1)*UStop(1,2)*UStop2(1,1)

        Cd(10953) = Cd(10952)*MStop2(3) - 
     &   2*Cd(10951)*MStop2(1)*UStopC(1,1)

        Cd(10954) = -6*UStop(2,2) + 
     &   2*UCStopC(3,1)*UStop(1,2)*UStop2(1,1)

        Cd(10955) = Cd(10954)*MStop2(3) - 
     &   2*Cd(10951)*MStop2(1)*UStopC(1,1)

	Cd(10956) = UStop(1,2)*UStop2(1,1) - Cr(91)*UStopC(2,1)

	Cd(10957) = 4*MStop2(1) + Cd(10922)*MStop2(3)

        Cd(10958) = 4*Cd(10956)*MStop2(1)*UStop(1,2)*
     &    UStop2(1,1) - Cd(10957)*UUStop(1,2)

	Cd(10959) = 6*MStop2(3)**2 - 2*MStop2(1)**2*UStop2(1,1)**2

        Cd(10960) = 4*UStop2(1,1)**2*UUStop(1,2) + 
     &   Cr(91)**2*UUStopC(2,1)

	Cd(10961) = -2*Cd(10960)*MStop2(1)**2 + Cd(10959)*UUStop(1,2)

	Cd(10962) = 6*MStop2(1) - MStop2(3)*UStop2(1,1)**2

        Cd(10963) = 2*Cd(10956)*MStop2(1)*UStop(1,2)*
     &    UStop2(1,1) - Cd(10962)*UUStop(1,2)

	Cd(10964) = 2*(Cd(10961) - 2*Cd(10963)*MStop2(3))

	Cd(10967) = 5*MStop2(1) - 2*MStop2(3)*UStop2(1,1)**2

        Cd(10968) = 3*Cd(10956)*MStop2(1)*UStop(1,2)*
     &    UStop2(1,1) - Cd(10967)*UUStop(1,2)

	Cd(10974) = Cd(10949) + 2*MUE*Cd(10946)*UStop2(1,1)*UStop2(1,2)

	Cd(10976) = MTy*UStop(1,1) + At*UStop(1,2)

	Cd(10977) = AtC*UStop(1,2)*UUStop(1,1) - 2*Cd(10976)*UUStop(1,2)

	Cd(10979) = MUE*UUStop(1,1) - 2*MUEC*UUStop(1,2)

	Cd(10980) = 2*MTy*UStop(1,1) - At*UStop(1,2)

	Cd(10981) = At*UUStop(1,2) + 2*MTy*UUStop(1,3)

	Cd(10982) = 2*Cd(10981) + Cd(10980)*UStop(1,2)*UStop2(1,2)**2

        Cd(10983) = AtC*Cd(10979)*UCStop(1,3) - 
     &   2*MUE*Cd(10976)*UStop(1,1)*UStop2(1,2)

	Cd(10984) = MUEC*Cd(10982) - Cd(10983)*UCStop(1,3)

	Cd(10985) = AtC*UStop(1,1) + 2*MTy*UStop(1,2)

	Cd(10987) = MTy*Cr(91) + At*UUStop(3,2)

	Cd(10989) = MUE*Cd(10980)*UStop(1,1) - 2*MTy*MUEC*UUStop(1,2)

	Cd(1091) = Cd(997) + 4*MStop2(3)*MStop2(4)

	Cd(10991) = 2*MTy*UUStop(1,2) + AtC*UUStop(1,3)

	Cd(10993) = 6*MStop2(3)**2 - 4*MStop2(1)**2*UStop2(1,1)**2

	Cd(10994) = 2*MTy*Cr(91) - At*UUStop(3,2)

	Cd(1092) = 22*MStop2(3)**2 - Cd(1091)*UStop2(1,2)**2

	Cd(11010) = 2*MUE*UCStop(3,3) + MUEC*UCStop(3,4)

	Cd(11011) = 2*MUE*UUStop(3,1) - MUEC*UUStop(3,2)

	Cd(11012) = Cd(11011)*UStopC(2,1) - 2*MUE*Cr(91)*UStopC(2,2)

        Cd(11013) = Cd(11010)*UStop(2,1) + 
     &   Cd(11012)*UCStop(1,3)*UStop2(1,2)

	Cd(11017) = 2*MTy*UStop(2,1) - At*UStop(2,2)

	Cd(11018) = MUE*Cd(11017)*UStop(1,1) - 2*MTy*MUEC*UUStop(3,2)

	Cd(11019) = 3*UUStop(3,3) + UUStop(3,4)

        Cd(11020) = Cd(11018)*UStopC(2,1) - 
     &   2*MTy*MUE*Cd(11019)*UStopC(2,2)

        Cd(11023) = MUEC*Cr(91)*UStop(2,2) - 
     &   4*MUE*UStop(1,1)*UUStop(2,1)

        Cd(11024) = 6*MStop2(3)**2*UUStop(1,1) - 
     &   2*Cr(91)**2*MStop2(1)**2*UUStopC(2,2)

        Cd(11025) = Cd(10946)*UStop2(1,2)**2 - 
     &   Cd(11023)*UCStop(3,1)*UStopC(2,1)

	Cd(11026) = Cd(11024) + 5*MStop2(3)**2*UUStop(1,1)

	Cd(11027) = -(MUE*Cd(11026)) + 2*Cd(11025)*MStop2(1)**2

        Cd(11028) = 3*Cd(11012)*MStop2(1) + 
     &   4*MUE*MStop2(3)*UStop(1,1)*UStop2(1,2)

	Cd(1095) = Cd(998) - Cd(1092)

        Cd(11029) = 5*Cd(11010)*MStop2(1)*UStop(2,1) + 
     &   Cd(11028)*UCStop(1,3)*UStop2(1,2)

	Cd(11030) = 2*(Cd(11029)*MStop2(3) + Cd(11027)*UStopC(1,2))

	Cd(11031) = MUE*Cd(11024) - 2*Cd(11025)*MStop2(1)**2

        Cd(11032) = Cd(11012)*MStop2(1) + 
     &   MUE*MStop2(3)*UStop(1,1)*UStop2(1,2)

        Cd(11033) = 3*Cd(11010)*MStop2(1)*UStop(2,1) + 
     &   Cd(11032)*UCStop(1,3)*UStop2(1,2)

	Cd(11034) = 4*Cd(11033)*MStop2(3) - 2*Cd(11031)*UStopC(1,2)

        Cd(11035) = 6*MStop2(3)**2*UStop(1,1) - 
     &   4*Cr(91)*MStop2(1)**2*UStop2(1,1)*UStopC(2,2)

	Cd(11036) = 2*MTy*Cr(58) + At*UUStop(3,2)

	Cd(11037) = 2*MTy*Cr(58)*UStop(1,1) - At*UStop(2,1)*UUStop(1,2)

        Cd(11038) = MUE*Cd(11037)*UStop(2,1) + 
     &   MUEC*Cd(11036)*UUStop(3,2)

	Cd(1097) = 10*MStop2(3) - 4*MStop2(2)*UStop2(1,2)**2

        Cd(11054) = AtC*UCStop(1,3)*UStop(1,1) - 
     &   At*UStop(1,2)*UStop2(1,2)

	Cd(11055) = At*UUStop(1,2) - 4*MTy*UUStop(1,3)

	Cd(11056) = -(At*Cd(11055)) + 8*MTy2*UUStop(1,1)

        Cd(11057) = 4*Cd(10976)*UStop(1,1)*UStop2(1,2) - 
     &   AtC*UCStop(1,3)*UUStop(1,1)

	Cd(1098) = 5*MStop2(3) + 2*MStop2(2)*UStop2(1,2)**2

	Cd(11059) = MStop2(3) + 2*At*MTy*UCStopC(1,3)

	Cd(11060) = At*UStopC(2,1) + 4*MTy*UStopC(2,2)

	Cd(11067) = 6*MStop2(3)**2 - 8*MStop2(1)**2*UStop2(1,1)**2

	Cd(1099) = Cd(1098) - 2*Cr(75)*MStop2(2)*UStop2(1,1)

	Cd(11069) = MTy*Cr(58) + At*UUStop(3,2)

        Cd(11072) = 3*Cd(11059)*MStop2(1) - 
     &   2*MTy2*MStop2(3)*UStop2(1,2)**3

        Cd(11073) = Cd(11072)*UStop(1,2) + 
     &   2*MTy*Cd(11060)*MStop2(1)*UStop2(1,2)**2*UUStop(3,2)

	Cd(11075) = Cd(11067) + 5*MStop2(3)**2

        Cd(11077) = 5*Cd(11059)*MStop2(1) - 
     &   8*MTy2*MStop2(3)*UStop2(1,2)**3

        Cd(11078) = Cd(11077)*UStop(1,2) + 
     &   6*MTy*Cd(11060)*MStop2(1)*UStop2(1,2)**2*UUStop(3,2)

	Cd(1100) = 2*Cd(1099)*UStop2(1,1) + Cd(1097)*UStop2(1,2)

	Cd(11080) = UStop(1,1)*UStop2(1,2) - Cr(91)*UStopC(2,2)

        Cd(11083) = 5*UStop2(1,2)**2*UUStop(1,1) + 
     &   Cr(91)**2*UUStopC(2,2)

        Cd(11084) = -2*Cd(11083)*MStop2(1)**2 + 
     &   6*MStop2(3)**2*UUStop(1,1)

	Cd(11085) = 6*MStop2(1) - MStop2(3)*UStop2(1,2)**2

        Cd(11086) = 2*Cd(9948)*Cd(11080)*MStop2(1)*UStop(1,1) - 
     &   Cd(11085)*UUStop(1,1)

	Cd(11087) = 2*(Cd(11084) - 2*Cd(11086)*MStop2(3))

        Cd(11088) = 2*Cd(11083)*MStop2(1)**2 - 
     &   11*MStop2(3)**2*UUStop(1,1)

	Cd(1101) = 2*MStop2(2) - 7*MStop2(3)

	Cd(11089) = 5*MStop2(1) - 2*MStop2(3)*UStop2(1,2)**2

        Cd(11090) = 3*Cd(9948)*Cd(11080)*MStop2(1)*UStop(1,1) - 
     &   Cd(11089)*UUStop(1,1)

	Cd(11091) = Cd(11088) + 2*Cd(11090)*MStop2(3)

	Cd(11105) = 1 - 3*UStop2(1,1)

	Cd(11106) = UUStop(3,3) + 2*UUStop(3,4)

	Cd(11107) = UStop(1,1)*UStop2(1,2) - Cd(11106)*UStopC(2,2)

	Cd(11108) = Cd(11105)*MStop2(1)**2 + 2*Cd(9948)*MStop2(3)**2

        Cd(11109) = 6*Cd(11107)*MStop2(1)*MStop2(3) + 
     &   4*Cd(11108)*UStop(1,1)

        Cd(11110) = 2*Cd(11107)*MStop2(1) + 
     &   5*Cd(9948)*MStop2(3)*UStop(1,1)

	Cd(11111) = Cd(11105)*MStop2(1)**2 + Cd(9948)*MStop2(3)**2

	Cd(11112) = Cd(11107)*MStop2(1)*MStop2(3) + Cd(11111)*UStop(1,1)

	Cd(11119) = -(MTy*Cr(58)) + At*UUStop(3,2)

	Cd(11122) = MTy2*Cr(91) + At*MTy*UUStop(3,2)

	Cd(11123) = 4*Cd(11122)*Cr(91) - At**2*UUStop(1,2)*UUStop(2,2)

	Cd(11132) = 2*MTy*Cr(91) + 3*At*UUStop(3,2)

	Cd(11137) = 2*UUStop(3,3) - UUStop(3,4)

	Cd(11138) = 4*UUStop(3,3) + 3*UUStop(3,4)

	Cd(11139) = 2*MTy*Cd(11137)*UStop(1,1) - At*Cd(11138)*UStop(1,2)

	Cd(11158) = 2*Cd(10976)*UStop(1,2) - AtC*UUStop(1,1)

	Cd(11160) = 2*Cd(10976)*UStop(2,2) - AtC*UUStop(3,1)

	Cd(11161) = UUStop(3,3) + 3*UUStop(3,4)

	Cd(11162) = MTy*Cd(11161)*UCStop(3,1) + At*Cr(91)*UCStop(3,4)

	Cd(11165) = MTy*Cr(91) - At*UUStop(3,2)

	Cd(11166) = At*Cr(91) + 2*MTy*UUStop(3,1)

	Cd(11167) = MUE*Cd(11166)*UStop(1,1) + MUEC*Cd(11165)*UStop(1,2)

        Cd(11168) = MUE*Cd(10925)*Cd(10980) - 
     &   6*MTy*MUEC*UStop(2,2)*UUStop(1,2)

        Cd(11169) = 2*MUE*Cd(11161)*UUStop(1,1) - 
     &   MUEC*Cd(11138)*UUStop(1,2)

	Cd(11170) = AtC*Cd(11169) + Cd(11168)*UStop(1,2)

	Cd(11171) = 2*Cd(11167)*UCStop(3,4) + Cd(11170)*UStopC(2,2)

	Cd(11175) = MUE*UUStop(2,1) - MUEC*UUStop(2,2)

	Cd(11176) = MUEC*Cd(9948)*Cr(91) + Cd(11175)*UCStop(3,3)

        Cd(11177) = MUEC*Cd(10921)*UStop(1,2) - 
     &   2*MUE*UStop(2,1)*UUStop(1,1)

        Cd(11178) = Cd(10979)*MStop2(3)*UCStop(1,3) + 
     &   Cd(11177)*MStop2(1)*UStopC(2,2)

        Cd(11179) = Cd(11178)*MStop2(3)*UCStop(1,3) + 
     &   Cd(11176)*MStop2(1)**2*UCStop(3,3)

	Cd(11180) = UUStop(3,3) + 4*UUStop(3,4)

	Cd(11181) = 2*MTy*Cd(11180)*UStop(1,1) + At*Cd(10921)*UStop(1,2)

	Cd(11182) = 2*MTy*Cr(58)*UStop(2,1) + At*UStop(1,1)*UUStop(2,2)

        Cd(11184) = 2*Cd(10976)*MStop2(3)*UStop(1,1)*
     &    UStop2(1,2) - Cd(11181)*MStop2(1)*UStopC(2,2)

	Cd(11187) = -2*MTy*Cr(58) + At*UUStop(3,2)

	Cd(1111) = 4*Cd(991) - Cd(1079)

        Cd(11198) = 2*AtC*UCStop(1,3)*UStop(1,1) - 
     &   2*At*UStop(1,2)*UStop2(1,2)

	Cd(11199) = -2*UCStop(3,1)*UStop(2,2) + 2*UStop(1,2)*UStop2(1,1)

        Cd(11201) = 2*MTy*Cd(11106)*UCStop(3,1) + 
     &   At*Cd(10925)*UCStop(3,4)

	Cd(11202) = At*Cd(11201) + Cd(11056)*UStop2(1,1)

	Cd(11203) = 2*MTy2*UStop(2,1) + At*MTy*UStop(2,2)

	Cd(11204) = 4*Cd(11203)*Cr(58) - At**2*UStop(1,2)*UUStop(2,2)

	Cd(11215) = MUE*UUStop(3,1) - MUEC*UUStop(3,2)

	Cd(11216) = Cd(11215)*UStopC(2,1) - MUE*Cr(91)*UStopC(2,2)

        Cd(11217) = 2*Cd(11216)*MStop2(1) - 
     &   2*Cd(10946)*MStop2(3)*UStopC(1,1)

	Cd(1114) = 6*MStop2(2) + MStop2(3)*UStop2(1,2)

        Cd(11219) = MUE**2*MStop2(1)*UCStop(1,3)**2*
     &    UStop2(1,1) - MUE2*Cd(10955)*UStopC(2,2)

        Cd(11220) = 4*Cd(11219)*MStop2(1)*UStop2(1,1) - 
     &   MUEC**2*Cd(10964)*UUStopC(1,1)

        Cd(11222) = MUE*Cd(11106)*UUStop(1,1) - 
     &   MUEC*Cd(10925)*UUStop(1,2)

	Cd(11233) = 5*Cd(10981) + 2*Cd(10980)*UStop(1,2)*UStop2(1,2)**2

        Cd(11235) = 4*Cd(10976)*MStop2(3)*UStop(1,1)*
     &    UStop2(1,2) - Cd(11181)*MStop2(1)*UStopC(2,2)

        Cd(11236) = 2*Cd(10979)*MStop2(3)*UCStop(1,3) + 
     &   Cd(11177)*MStop2(1)*UStopC(2,2)

	Cd(11237) = -(MUE*Cd(11235)) + AtC*Cd(11236)

	Cd(11247) = -(UCStop(3,1)*UStop(2,2)) + UStop(1,2)*UStop2(1,1)

	Cd(11250) = Cd(11235) - AtC*MStop2(2)*UCStop(1,3)*UUStop(1,1)

	Cd(11253) = Cd(11087) + 2*Cd(11091)

	Cd(11256) = Cd(11109) - 4*Cd(11112)

        Cd(11261) = 2*Cd(11216)*MStop2(1) - 
     &   4*Cd(10946)*MStop2(3)*UStopC(1,1)

	Cd(11263) = Cd(10928) - 4*Cd(10931)

	Cd(11266) = Cd(10953) - Cd(10955)

	Cd(11269) = 5*MUE**2*UUStop(1,1) - 2*Cd(10948)*UUStopC(1,1)

        Cd(11275) = MUEC*Cr(91)*UStop(2,2) - 
     &   2*MUE*UStop(1,1)*UUStop(2,1)

	Cd(11276) = 2*MUE2*UUStop(2,1) - MUEC**2*UUStop(2,2)

        Cd(11277) = -(Cd(11276)*UUStop(1,2)) + 
     &   MUE**2*UUStop(1,1)*UUStop(2,1)

        Cd(11279) = MUEC*Cr(91)*UStop(1,2) - 
     &   2*MUE*UStop(2,1)*UUStop(1,1)

        Cd(11280) = MUE*Cd(10925)*UUStop(2,1) - 
     &   MUEC*Cd(11106)*UUStop(2,2)

        Cd(11281) = 2*Cd(11279)*MStop2(2)*UCStopC(3,3) - 
     &   2*Cd(11280)*MStop2(3)*UStopC(2,1)

	Cd(11282) = 3*MUE + MUEC*UCStopC(1,3)**2

        Cd(11283) = Cd(11282)*UStop(1,2) - 
     &   MUE*UCStop(3,1)*UStop(2,2)*UStop2(1,2)

        Cd(11284) = MUE*Cr(91)*UStop(1,1) - 
     &   2*MUEC*UStop(2,2)*UUStop(1,2)

        Cd(11285) = 2*Cd(11283)*MStop2(3) - 
     &   2*Cd(11284)*MStop2(2)*UUStopC(3,1)

	Cd(11286) = UCStop(3,1)*UUStop(2,2) + UStop2(1,2)*UUStop(3,2)

	Cd(1121) = MStop2(3) + 4*MStop2(4)

	Cd(11290) = MUEC*Cd(9948)*Cd(11285) - MUE*Cd(11281)*UUStopC(2,2)

        Cd(11293) = Cd(11215)*MStop2(2)*UStopC(1,1) + 
     &   Cd(11175)*MStop2(3)*UStopC(2,1)

	Cd(11297) = -(MUEC**2*UCStopC(1,3)**2) + 4*MUE2*UStop2(1,2)**2

	Cd(1122) = 8*MStop2(2) - 3*MStop2(3)

	Cd(11299) = 2*MUE2*UUStop(1,1) - MUEC**2*UUStop(1,2)

        Cd(11300) = MUE**2*UUStop(1,1)*UUStop(2,1) - 
     &   Cd(11299)*UUStop(2,2)

        Cd(11305) = MUE*Cr(91)*UStop(2,1) - 
     &   2*MUEC*UStop(1,2)*UUStop(2,2)

	Cd(11306) = UCStop(3,1)*UUStop(2,2) + 2*UStop2(1,2)*UUStop(3,2)

	Cd(1123) = MStop2(1) - Cd(1122)*UStop2(1,2)**2

	Cd(11315) = 6 - 4*UStop2(1,1)**2

	Cd(11316) = 3 - UStop2(1,1)**2

	Cd(11317) = 3*MStop2(2) + 2*MTy2*Cd(11316)*UStop2(1,1)

        Cd(11318) = 2*MTy2*Cd(9948)*Cd(11315)*MStop2(2) + 
     &   Cd(11317)*MStop2(3)

	Cd(1124) = MStop2(2) + 5*MStop2(3)

        Cd(11319) = 2*(Cd(11318)*MStop2(3) - 
     &     8*MTy2*MStop2(2)**2*UStop2(1,1)*UStop2(1,2)**2)

	Cd(11320) = 6*MTy*UStop(2,1) + At*Cd(11316)*UStop(2,2)

        Cd(11321) = 2*Cd(11320)*MStop2(3) - 
     &   4*Cd(11119)*MStop2(2)*UStop2(1,1)*UStopC(1,2)

        Cd(11323) = 4*Cd(11316)*MStop2(3)**2 - 
     &   4*Cr(26)*MStop2(2)*UStop2(1,1)**2

        Cd(11324) = Cd(11316)*UStop(1,1) + 
     &   UCStop(3,1)*UStop(2,1)*UStop2(1,1)

        Cd(11325) = Cd(11324)*MStop2(3) + 
     &   2*MStop2(2)*UStop(1,1)*UStop2(1,1)**2

        Cd(11327) = AtC*MStop2(1)*UCStop(1,3) + 
     &   4*MTy*MStop2(2)*UStop2(1,1)

	Cd(11328) = MStop2(1)*UStop2(1,1) - MStop2(2)*UStop2(1,2)

	Cd(1125) = 2*Cd(1124) - Cd(1122)*UStop2(1,2)**2

	Cd(11330) = MTy*UStop(2,1) + At*UStop(2,2)

        Cd(11331) = MStop2(3)*UStop(2,1)*UStop2(1,1) + 
     &   Cr(91)*MStop2(2)*UStopC(1,2)

        Cd(11332) = 2*Cd(11330)*Cd(11331)*MStop2(3)*
     &    UCStop(1,3) + Cd(11037)*MStop2(2)**2*UCStopC(3,4)*
     &    UStopC(1,2)

	Cd(11333) = 2*MTy*Cd(10925)*UStop(2,1) + At*Cd(11106)*UStop(2,2)

	Cd(11334) = At*Cd(10925) + 6*MTy*UUStop(3,1)

	Cd(11340) = 3*MStop2(3) + 4*MStop2(2)*UStop2(1,1)*UStop2(1,2)

	Cd(11341) = 3 - UStop2(1,1)

	Cd(11343) = 3 + UStop2(1,1)**2

	Cd(11344) = Cd(11343)*MStop2(3)**2 + MStop2(2)**2*UStop2(1,2)**2

	Cd(11348) = Cd(11316)*UUStop(3,3) + 3*Cd(10918)*UUStop(3,4)

        Cd(11349) = Cd(11348)*MStop2(3) - 
     &   2*Cr(91)*MStop2(2)*UStop2(1,1)*UStop2(1,2)

        Cd(11350) = Cd(11316)*MStop2(3)**2*UStop(2,1)*
     &    UStop2(1,1) + Cd(11349)*MStop2(2)*UStopC(1,2)

	Cd(11355) = 6*MTy*UStop(1,1) + 3*At*Cd(10918)*UStop(1,2)

	Cd(11356) = MStop2(2) + 4*MTy2*UStop2(1,1)

        Cd(11357) = 3*Cd(11356)*UStop(2,1) + 
     &   2*At*MTy*Cd(11316)*UStop(2,2)*UStop2(1,1)

        Cd(11358) = Cd(11355)*MStop2(3)*UCStopC(3,1) + 
     &   2*Cd(11069)*MStop2(2)*UStop2(1,1)*UStop2(1,2)*UStopC(1,2)

        Cd(11359) = -4*MTy*Cd(11358)*MStop2(2) + 
     &   2*Cd(11357)*MStop2(3)**2

        Cd(11361) = Cd(11166)*MStop2(3)*UStopC(1,1) + 
     &   Cd(11132)*MStop2(2)*UStopC(1,2)

        Cd(11363) = 2*Cd(11361)*UStop(1,1) - 
     &   AtC*MStop2(1)*UCStop(1,3)*UUStop(3,1)

        Cd(11365) = 4*MTy*Cd(11343)*UStop(2,1) + 
     &   2*At*Cd(11316)*UStop(2,2)

        Cd(11376) = At*UUStop(1,1)*UUStop(2,2) + 
     &   2*MTy*Cr(58)*UUStop(3,1)

        Cd(11381) = At*UUStop(1,2)*UUStop(2,1) - 
     &   2*MTy*Cr(58)*UUStop(3,1)

	Cd(11382) = MUE*Cd(11381) - MUEC*Cd(11036)*UUStop(3,2)

	Cd(11385) = MUE*Cd(10980)*UStop(2,1) - 2*MTy*MUEC*UUStop(3,2)

	Cd(11386) = AtC*Cd(11011)*UStop(2,1) + Cd(11385)*UStop(2,2)

	Cd(11388) = MUEC*Cr(91)*UStop(2,2) - MUE*UStop(1,1)*UUStop(2,1)

	Cd(11389) = MUE*Cd(11166)*UStop(2,1) + MUEC*Cd(11165)*UStop(2,2)

	Cd(11390) = AtC*Cd(11388)*UStop(2,1) + Cd(11389)*UStop(2,2)

	Cd(11391) = MTy*Cd(10925)*UStop(2,1) + At*Cd(11106)*UStop(2,2)

	Cd(11392) = 2*Cd(11391)*UStop(2,2) - AtC*Cd(10925)*UUStop(2,1)

	Cd(11393) = MUEC*Cd(11392) - MUE*Cd(11334)*UUStop(2,1)

	Cd(11394) = MUE*UUStop(3,1) + MUEC*UUStop(3,2)

	Cd(11398) = 2*Cd(11330)*UStop(1,2) - AtC*UUStop(3,1)

	Cd(11400) = AtC*Cd(11388) + MUE*Cd(11330)*Cr(91)

        Cd(11401) = Cd(11400)*UStop(2,1) + 
     &   MUEC*Cd(11017)*UStop(1,2)*UUStop(2,2)

        Cd(11405) = MUEC*Cr(91)*UStop(1,2) - 
     &   4*MUE*UStop(2,1)*UUStop(1,1)

        Cd(11406) = 2*MUE*Cd(11324)*MStop2(3)*UStop(2,1) - 
     &   Cd(11405)*MStop2(2)*UStop2(1,1)*UStopC(1,1)

        Cd(11409) = 2*MStop2(2)*UStop(1,1)*UStop2(1,1) - 
     &   2*Cd(10925)*MStop2(3)*UStopC(2,2)

        Cd(11410) = -(AtC*Cd(11409)*UStop(2,1)) + 
     &   2*Cd(11132)*MStop2(2)*UStop2(1,1)

	Cd(11412) = Cd(11106)*UStop2(1,1)**2 - 3*UUStop(3,4)

        Cd(11413) = 2*(Cd(11412)*MStop2(3) - 
     &     6*MStop2(2)*UStop2(1,1)**2*UUStop(3,4))

        Cd(11414) = 6*MTy*Cd(10918)*MStop2(3) + 
     &   2*Cd(10994)*MStop2(2)*UUStopC(3,3)

	Cd(11415) = AtC*Cd(11413) - 2*Cd(11414)*UUStop(3,2)

	Cd(1134) = MStop2(1) + 24*MTy2*UStop2(1,2)

        Cd(11424) = 2*Cd(11279)*MStop2(2)*UCStopC(3,3) - 
     &   3*Cd(11280)*MStop2(3)*UStopC(2,1)

	Cd(11425) = MUE + MUEC*UCStopC(1,3)**2

        Cd(11426) = Cd(11425)*UStop(1,2) - 
     &   MUE*UCStop(3,1)*UStop(2,2)*UStop2(1,2)

        Cd(111) = MGlpTmSt4(1)*UStop(1,1) + 
     &   2*M_3*MTy*MGlpTmSt2(1)*UStop(1,2)

	Cd(1135) = 6 + 5*UStop2(1,2)**2

	Cd(11434) = Cd(11281) - Cd(11424)

	Cd(1136) = 2 - UStop2(1,2)**2

        Cd(11443) = Cd(11215)*MStop2(2)*UStopC(1,1) + 
     &   2*Cd(11175)*MStop2(3)*UStopC(2,1)

	Cd(11444) = MUE**2*Cd(1136)*UUStop(2,1) + Cd(11297)*UUStop(2,2)

	Cd(11447) = MUE*Cd(1136) + MUEC*UCStopC(1,3)**2

	Cd(11457) = 5 - 6*UStop2(1,1)**2

	Cd(11458) = 8 - 3*UStop2(1,1)**2

	Cd(11459) = MStop2(2) + MTy2*Cd(11458)*UStop2(1,1)

        Cd(11461) = Cd(9948)*Cd(11457)*MStop2(3) - 
     &   4*MStop2(2)*UStop2(1,1)*UStop2(1,2)**2

	Cd(11463) = 2*MTy*UStop(2,1) + At*Cd(10918)*UStop(2,2)

	Cd(11467) = 10 - 4*UStop2(1,1)**2

	Cd(11468) = 4*Cd(11316) + Cd(11467)

        Cd(11470) = Cd(10918)*UStop(1,1) + 
     &   UCStop(3,1)*UStop(2,1)*UStop2(1,1)

        Cd(11474) = AtC*UCStop(1,3)*UStop(2,1) + 
     &   At*UStop(2,2)*UStop2(1,1)

        Cd(11476) = 2*MStop2(3)*UStop(2,1)*UStop2(1,1) + 
     &   Cr(91)*MStop2(2)*UStopC(1,2)

        Cd(11477) = AtC*Cd(11334)*UCStopC(3,4) + 
     &   At*Cd(11333)*UStopC(1,1)

        Cd(11478) = Cd(11332) + 
     &   Cd(11330)*Cd(11476)*MStop2(3)*UCStop(1,3)

	Cd(11479) = MStop2(2)**2 - 2*MStop2(3)**2

	Cd(11480) = 11*MStop2(3) + 12*MStop2(2)*UStop2(1,1)*UStop2(1,2)

	Cd(11482) = 9 - 2*UStop2(1,1)

	Cd(11483) = 11 - 4*UStop2(1,1)**2

	Cd(11486) = 11 + 4*UStop2(1,1)**2

        Cd(11487) = MStop2(2)*UStop2(1,2)*UUStop(2,3) + 
     &   Cr(2)*UCStopC(3,1)*UUStop(3,3)

	Cd(1141) = 4*MStop2(2) + MStop2(3)

        Cd(11488) = 2*Cd(11487)*MStop2(2)*UStop2(1,2) + 
     &   Cd(11486)*MStop2(3)**2*UUStop(2,3)

	Cd(11490) = 2 - UStop2(1,1)**2

	Cd(11491) = 3*Cd(11316) + Cd(11490)

	Cd(1142) = 2*Cd(1124) + Cd(1141)*UStop2(1,2)**2

	Cd(11502) = 1 - 3*UStop2(1,1)**2

	Cd(11503) = 2*MTy*UStop(1,1) + At*Cd(11502)*UStop(1,2)

	Cd(11505) = MStop2(2) + 16*MTy2*UStop2(1,1)

        Cd(11506) = Cd(11505)*UStop(2,1) + 
     &   At*MTy*Cd(11458)*UStop(2,2)*UStop2(1,1)

	Cd(1143) = 12 - 5*UStop2(1,2)**2

        Cd(11510) = 3*Cd(11166)*MStop2(3)*UStopC(1,1) + 
     &   2*Cd(11132)*MStop2(2)*UStopC(1,2)

	Cd(11515) = 2 + UStop2(1,1)**2

	Cd(11516) = 2*MTy*Cd(11515)*UStop(2,1) + At*Cd(11490)*UStop(2,2)

	Cd(11517) = 3*Cd(11365) + 2*Cd(11516)

	Cd(11530) = Cd(11386)*UStop2(1,1) + 2*Cd(11390)*UStopC(2,1)

        Cd(11532) = Cd(11386)*UStop2(1,1)*UStopC(1,1) + 
     &   2*Cd(11390)*UUStopC(3,1)

        Cd(11537) = 2*MStop2(2)*UStop(1,1)*UStop2(1,1) - 
     &   3*Cd(10925)*MStop2(3)*UStopC(2,2)

	Cd(1146) = Cd(1142) - Cd(1141)*Cr(75)*UStop2(1,1)

        Cd(11538) = -(AtC*Cd(11537)*UStop(2,1)) + 
     &   2*Cd(11132)*MStop2(2)*UStop2(1,1)

	Cd(11540) = -(Cd(11106)*UStop2(1,1)**2) + UUStop(3,4)

	Cd(11543) = 5 - 9*UStop2(1,1)**2

        Cd(11544) = MTy*Cd(11543)*MStop2(3) + 
     &   2*Cd(10994)*MStop2(2)*UUStopC(3,3)

	Cd(1148) = 12*MStop2(2) - 7*MStop2(3)

	Cd(11564) = Cd(11444) + 2*MUE*Cd(11175)*UStop2(1,1)*UStop2(1,2)

	Cd(1149) = 2*MStop2(1) + Cd(1121)*UStop2(1,2)

        Cd(11568) = Cd(11434)*UStopC(1,2) - 
     &   Cd(11280)*MStop2(3)*UUStopC(3,4)

	Cd(11571) = Cd(11467) + Cd(11490)

	Cd(11572) = Cd(11571)*MStop2(3) + 4*MStop2(2)*UStop2(1,1)**2

        Cd(11573) = 4*Cd(11470)*MStop2(2)*UCStopC(3,4) - 
     &   Cd(11572)*UCStop(1,3)*UStop(2,1)

	Cd(11574) = 1 - 2*UStop2(1,1)**2

	Cd(11575) = 4 + Cd(11490)

	Cd(1150) = 5 + UStop2(1,2)

	Cd(11580) = Cd(10918)*UUStop(3,3) + Cd(11502)*UUStop(3,4)

	Cd(11584) = 2*Cd(11476) + MStop2(3)*UStop(2,1)*UStop2(1,1)

	Cd(11585) = 3*MStop2(3) + 2*MStop2(2)*UStop2(1,1)*UStop2(1,2)

        Cd(11586) = 5*MStop2(3)*UCStopC(3,1)**2*UUStop(1,1) + 
     &   2*Cd(11585)*UUStop(2,1)

        Cd(11587) = Cd(10922)*MStop2(3) - 
     &   4*MStop2(2)*UStop2(1,1)*UStop2(1,2)

	Cd(11588) = 12 + 5*UStop2(1,1)**2

        Cd(11589) = Cd(11588)*MStop2(3)*UUStop(2,3) + 
     &   4*MStop2(2)*UCStopC(3,1)*UStop2(1,2)*UUStop(3,3)

	Cd(11592) = 4*MTy2*UStop(2,1) + At*MTy*Cd(11490)*UStop(2,2)

        Cd(11596) = 2*Cd(11166)*UStop2(1,1) + 
     &   AtC*UCStop(1,3)*UUStop(3,1)

	Cd(1153) = 4 - 3*UStop2(1,2)

	Cd(11610) = MTy*Cd(11580) + AtC*Cd(11470)*UStop(2,1)

	Cd(11611) = -(AtC*Cd(11540)) - 2*MTy*Cd(11502)*UUStop(3,2)

	Cd(11614) = MUE*Cd(11610) + Cd(11401)*UUStopC(2,1)

	Cd(1154) = Cd(1153)*MStop2(1) - 5*MStop2(2)*UStop2(1,2)

	Cd(11622) = Cd(9948) - UStop2(1,1)

	Cd(11623) = Cd(9948) - 4*UStop2(1,1)

	Cd(1155) = 12*MStop2(3) - Cd(1154)*UStop2(1,2)

	Cd(11629) = At*AtC + 4*MTy2

        Cd(11630) = Cd(9948)*Cd(11629)*MStop2(2) - 
     &   4*At*AtC*MStop2(1)*UStop2(1,1)

	Cd(11631) = At*AtC - 4*MTy2

	Cd(11633) = 2*UStop2(1,1)**2 + Cd(10904)*UStop2(1,2)

        Cd(11634) = Cd(11633)*MStop2(2) - 
     &   4*MStop2(1)*UStop2(1,1)*UStop2(1,2)

	Cd(1156) = 10*MStop2(2) - 5*MStop2(3)

        Cd(11637) = Cr(75)**2*MStop2(2) + 
     &   2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(1157) = 2*MStop2(1) - Cd(1156)*UStop2(1,2)

	Cd(11652) = MUEC*Cr(91)*UStopC(2,1) + Cd(11215)*UStopC(2,2)

	Cd(11653) = MUE**2*UUStop(1,1) - 2*MUE2*UUStop(1,2)

        Cd(11654) = MUE2*MStop2(3)**2*UUStop(1,1) + 
     &   MUEC**2*MStop2(1)*MStop2(2)*UUStop(1,2)

        Cd(11655) = 4*(Cd(11653)*MStop2(1)*MStop2(2)*
     &      UUStop(2,1) + Cd(11654)*UUStop(2,2))

	Cd(11656) = UUStop(1,2)*UUStop(2,1) + UUStop(1,1)*UUStop(2,2)

	Cd(1158) = 4 + Cd(1136)

        Cd(11659) = MUE*MStop2(3)**2*UUStop(1,1) + 
     &   2*MUEC*MStop2(1)*MStop2(2)*UUStop(1,2)

        Cd(11666) = 2*MUE*MStop2(1)*MStop2(2)*UUStop(2,1) + 
     &   MUEC*MStop2(3)**2*UUStop(2,2)

        Cd(11669) = MUE*MStop2(3)**2*UUStop(2,1) + 
     &   2*MUEC*MStop2(1)*MStop2(2)*UUStop(2,2)

        Cd(11672) = MUE**2*UUStop(1,1)*UUStop(2,1) + 
     &   MUEC**2*UUStop(1,2)*UUStop(2,2)

	Cd(1160) = 2*MStop2(2) + 3*MStop2(3)

        Cd(11677) = MUE**2*MStop2(1)*MStop2(2)*UUStop(1,1) + 
     &   MUE2*MStop2(3)**2*UUStop(1,2)

        Cd(11678) = 4*Cd(11677)*UUStop(2,1) - 
     &   4*Cd(11299)*MStop2(1)*MStop2(2)*UUStop(2,2)

        Cd(11680) = 2*MUE*MStop2(1)*MStop2(2)*UUStop(1,1) + 
     &   MUEC*MStop2(3)**2*UUStop(1,2)

	Cd(1161) = Cd(1157) + Cd(1160)*UStop2(1,1)

	Cd(11691) = Cd(9948)*MStop2(2) - 4*Cr(26)*UStop2(1,1)

	Cd(11693) = 4*MStop2(3) - Cd(9955)*MStop2(2)*UStop2(1,2)

	Cd(11700) = 1 + 2*Cd(9948)*UStop2(1,1)

	Cd(1163) = 6*MStop2(2)**2 - MStop2(3)*MStop2(4)

	Cd(11725) = At*Cr(28) - AtC*UCStop(1,3)**2

	Cd(114) = 2*MGlpTmSt4(1) - Cd(111)*UStopC(1,1)

        Cd(11727) = M_3*UCStopC(1,3)**2 + 
     &   Cd(9963)*UStop2(1,1)*UStop2(1,2)

	Cd(11738) = 5 + 2*UStop2(1,1)

	Cd(11739) = 1 - 12*UStop2(1,2)

	Cd(11740) = 5 - 12*UStop2(1,2)**2

	Cd(11745) = -8*Cd(10176) + AtC*Cd(9969)*Cr(28)

	Cd(11747) = 5 - 12*Cd(9948)

	Cd(11748) = 4*M_3*Cr(28) + Cd(10209)*UCStop(1,3)**2

	Cd(11749) = -4*Cd(10176) + AtC*Cr(28)

        Cd(11753) = -(Cd(10209)*UCStop(1,3)**2) + 
     &   8*M_3*UStop2(1,1)*UStop2(1,2)

	Cd(11754) = -4*Cd(10039) + AtC*UStop2(1,1)*UStop2(1,2)

	Cd(11759) = 3*MStop2(3) + 2*MTy2*Cd(9955)*UStop2(1,1)

	Cd(11760) = MStop2(3) + MGlpT2*Cd(9955)*UStop2(1,2)

	Cd(11764) = 3*MStop2(3) + 4*MTy2*Cd(9955)*UStop2(1,1)

	Cd(11791) = Cr(56)*MStop2(1) + MGlmT2*UStop2(1,1)

	Cd(11793) = MGlmT2*Cr(85) + Cr(54)*MStop2(1)

	Cd(11795) = 2*M_3*Cr(28) - Cd(10398)*UCStop(1,3)**2

	Cd(11796) = M_3*Cr(28) + 2*M_3C*UCStop(1,3)**2

	Cd(11797) = MGlmT2*Cd(11795) + 2*Cd(11796)*MStop2(1)

	Cd(11803) = 2*M_3C*Cr(28) - Cd(10140)*UCStopC(1,3)**2

	Cd(11804) = MGlmT2*Cd(11803) + 2*Cd(10360)*MStop2(1)

	Cd(1174) = 12*MStop2(3) - Cd(1121)*UStop2(1,2)**2

	Cd(1176) = Cd(1174) - Cd(1121)*UStop2(1,1)**2

	Cd(11850) = 2*M_3C - AtC*Cd(9955)

        Cd(11851) = Cd(11850)*UCStop(1,3)**2 + 
     &   2*M_3*UStop2(1,1)*UStop2(1,2)

	Cd(1179) = 2 + UStop2(1,2)**2

	Cd(1180) = Cd(1136) + UStop2(1,1)**2

	Cd(11877) = 1 + Cd(9948)

	Cd(11879) = -16*M_3C + AtC*Cd(9969)

	Cd(11912) = -(Cd(10398)*UCStop(1,3)**2) + 2*M_3*UStop2(1,1)**2

	Cd(1186) = 12*MStop2(3) - Cd(1141)*UStop2(1,2)**2

	Cd(11942) = MStop2(1) + Cd(9955)*MStop2(3)

	Cd(11968) = 3*M_3C + AtC

	Cd(11969) = Cd(11968)*UCStop(1,3) + Cd(10146)*UCStopC(1,3)

	Cd(11975) = 1 - 4*Cd(9948)

	Cd(12005) = M_3*UCStopC(1,3)**2 + M_3C*Cd(9948)*UStop2(1,1)

	Cd(12009) = 5 + 24*UStop2(1,1)**2

	Cd(12010) = 5 + 48*UStop2(1,1)**2

	Cd(12016) = -8*M_3 + At

        Cd(12017) = -(Cd(12016)*UCStopC(1,3)**2) + 
     &   8*M_3C*UStop2(1,1)*UStop2(1,2)

	Cd(12019) = 7 + 12*UStop2(1,2)

	Cd(12046) = At*UStop2(1,1) + 2*M_3*UStop2(1,2)

        Cd(12047) = Cd(12046)*UCStopC(1,3) + 
     &   2*M_3C*UCStop(1,3)*UStop2(1,2)

	Cd(12054) = 3 - 8*UStop2(1,2)

	Cd(12061) = 3 - 4*Cd(9948)

	Cd(12071) = 2*MUE*UUStop(1,1) + MUEC*UUStop(1,2)

	Cd(12075) = MUE*UUStop(1,1) + 2*MUEC*UUStop(1,2)

        Cd(12076) = MUEC*Cd(11019)*UStop(1,2) + 
     &   2*MUE*UStop(2,1)*UUStop(1,1)

        Cd(12077) = MUE*Cd(11161)*UStop(1,1) + 
     &   2*MUEC*UStop(2,2)*UUStop(1,2)

        Cd(12078) = MUEC*Cd(12077)*UStopC(2,1) + 
     &   MUE*Cd(12076)*UStopC(2,2)

	Cd(12079) = MUEC*Cr(91)*UStop(2,2) + MUE*UStop(1,1)*UUStop(2,1)

        Cd(12080) = MUE*Cd(11019)*UStop(2,1) + 
     &   2*MUEC*UStop(1,2)*UUStop(2,2)

        Cd(12085) = 2*MUEC*UCStopC(1,3)*UStop(1,2) - 
     &   MUE*Cr(91)*UStopC(2,2)

	Cd(12086) = MUE*UCStop(3,3) + 2*MUEC*UCStop(3,4)

	Cd(1202) = 1 + 3*UStop2(1,2)**2

        Cd(12094) = 3*MUE*UCStop(1,3)*UStop(1,1) - 
     &   2*MUEC*Cd(9948)*UStop(1,2)

	Cd(12107) = 1 + 4*UStop2(1,1)

	Cd(12108) = 4*Cd(9948)*UStop2(1,1) - Cd(12107)*UStop2(1,2)

	Cd(12109) = Cd(9948) + 4*UStop2(1,1)

	Cd(12110) = Cd(9948)*UStop2(1,1)**2 + Cd(12109)*UStop2(1,2)**2

	Cd(1204) = 1 + 2*UStop2(1,2)**2

	Cd(12114) = UStop2(1,1) - 2*UStop2(1,2)

	Cd(12115) = Cd(12114)*MStop2(1) - MStop2(3)*UStop2(1,2)

	Cd(12117) = Cd(11877)*MStop2(3) + MStop2(1)*UStop2(1,2)

	Cd(12122) = -(MUE**2*UCStop(1,3)**2) + MUE2*UStop2(1,2)**2

	Cd(118) = MGl2 - MStop2(1)

        Cd(12124) = MUE**2*Cr(56)*UCStop(1,3)**2 + 
     &   MUE2*Cd(12114)*UStop2(1,1)*UStop2(1,2)

	Cd(12125) = 4*UStop2(1,1) - UStop2(1,2)

	Cd(12137) = -(At*UCStopC(1,3)) + 2*MTy*UStop2(1,2)

	Cd(12139) = 4*MTy*UStop(1,1) + At*UStop(1,2)

	Cd(12140) = Cd(12139)*UStop(1,2) + AtC*UUStop(1,1)

	Cd(12143) = 2*MTy*Cr(91) + At*UUStop(3,2)

	Cd(12147) = MUE*Cd(12139)*UStop(1,1) + 2*MTy*MUEC*UUStop(1,2)

	Cd(12148) = MUE*UUStop(1,1) + MUEC*UUStop(1,2)

	Cd(12149) = AtC*Cd(12148)*UStop(1,1) + Cd(12147)*UStop(1,2)

        Cd(12150) = MUE*Cd(12139)*Cr(91)*UStopC(2,1) - 
     &   4*MTy*MUEC*UCStopC(1,3)*UUStop(1,2)

	Cd(12156) = MTy*Cd(11019) + At*UUStop(3,2)

	Cd(12158) = At*Cd(11161) + 8*MTy*UUStop(3,1)

	Cd(1209) = MStop2(1)**2 + MStop2(3)**2

        Cd(12164) = Cd(12158)*UCStop(3,1) - 
     &   4*AtC*UCStop(1,3)*UUStop(1,1)

	Cd(12165) = 3*UStop2(1,1) - UStop2(1,2)

	Cd(1210) = MStop2(3) - 5*MStop2(2)*UStop2(1,2)

	Cd(12179) = Cd(9948)*Cd(10987) - AtC*Cr(91)*UCStop(1,3)

	Cd(1211) = 2*MStop2(1)*MStop2(3) - 2*MStop2(2)**2*UStop2(1,2)

	Cd(12186) = 3*MTy2 - MStop2(4)

	Cd(12190) = MTy*Cd(11019)*UStop(2,1) + At*Cr(91)*UStop(2,2)

	Cd(1212) = 2*Cd(1209)*UStop2(1,1)**2 - Cd(1211)*UStop2(1,2)

	Cd(12194) = 3*MUE*UCStop(1,3)*UUStop(1,1) - MUEC*UUStop(1,3)

	Cd(1213) = UStop2(1,1)**4 + UStop2(1,2)**4

	Cd(12203) = 8*MTy*UStop(1,1) - 2*At*UStop(1,2)

        Cd(1214) = Cd(1213)*MStop2(2) - 
     &   2*Cd(1210)*UStop2(1,1)**2*UStop2(1,2)

	Cd(12218) = 2*MTy*Cd(11106) + 3*At*UUStop(3,2)

	Cd(12219) = AtC*Cr(91) + 4*MTy*UUStop(3,2)

	Cd(12220) = MUE*Cd(12190) + MUEC*Cd(12219)*UStop(2,2)

	Cd(1215) = Cd(1214)*MStop2(1) - Cd(1212)*UStop2(1,1)*UStop2(1,2)

	Cd(12225) = 4*MTy*UStop(1,1) + 3*At*UStop(1,2)

	Cd(12226) = Cd(12225)*UStop(1,2) - AtC*UUStop(1,1)

	Cd(1216) = 3*MStop2(1)**2 + MStop2(3)*MStop2(4)

	Cd(12234) = At*MUE*UStop(1,1) + 2*MTy*MUEC*UStop(1,2)

	Cd(12241) = MTy*Cr(91) + AtC*UUStop(3,1)

	Cd(12242) = AtC*UStop(2,1) + 2*MTy*UStop(2,2)

	Cd(1217) = MStop2(3) - 2*MStop2(2)*UStop2(1,2)

        Cd(12243) = MUE*Cd(10987)*UStop(2,1) + 
     &   MUEC*Cd(12242)*UUStop(3,2)

	Cd(1218) = Cd(1217)*MStop2(1) + Cd(1216)*UStop2(1,1)

	Cd(12253) = AtC*Cr(91)*UStop(2,1) + 2*Cd(12143)*UStop(2,2)

        Cd(12259) = MUEC*Cd(11161)*UStop(2,2) + 
     &   4*MUE*UStop(1,1)*UUStop(2,1)

	Cd(12260) = 4*MTy*Cd(11019)*UStop(2,1) + At*Cd(11161)*UStop(2,2)

	Cd(12261) = MUE*Cd(12260) + 8*MTy*MUEC*UStop(1,2)*UUStop(2,2)

	Cd(12262) = AtC*Cd(12259) + Cd(12261)

	Cd(12266) = MGl2**2 + MStop2(1)**2

	Cd(12277) = 1 + 12*UStop2(1,1)

	Cd(12278) = MTy2 + Cd(12277)*MStop2(1)

	Cd(12286) = 3*Cd(9948)*MStop2(1) + MStop2(2)

	Cd(1222) = 7*MStop2(2)**2 - 3*MStop2(3)*MStop2(4)

	Cd(1223) = 42*MStop2(2) - 2*MStop2(3)*UStop2(1,2)

	Cd(1224) = 2*Cd(1222) + Cd(1223)*MStop2(1)

	Cd(12313) = At*MGlpT2 - 8*M_3*MTy2

	Cd(1226) = 6*MStop2(2) - MStop2(3)

	Cd(1227) = 6*MStop2(2)**2 - Cd(1226)*MStop2(3)

	Cd(1228) = 6*MStop2(1)*MStop2(2) + Cd(1227)*UStop2(1,2)**2

	Cd(1231) = MStop2(2)**2 - Cr(1)*MStop2(3)

	Cd(12384) = (4*At*AtC + 36*MStop2(3))/MStop2(3)

	Cd(1232) = 8*MStop2(1)*MStop2(2) - 8*Cd(1231)*UStop2(1,2)

	Cd(12405) = MStop2(3) + 2*Cd(9955)*MStop2(1)*UStop2(1,1)

	Cd(12406) = Cd(9955)*Cd(10500) + MStop2(3)

	Cd(121) = MGl2 + MStop2(1)

	Cd(1235) = MStop2(3)**2 + Cd(1227)*UStop2(1,2)

	Cd(1237) = 2 - 5*UStop2(1,2)

	Cd(12454) = M_3C*At - MTy2

	Cd(12470) = MGlmT2 + Cd(9969)*MStop2(1)

	Cd(12506) = Cd(10070)*MStop2(1) + MStop2(3)

	Cd(12515) = 8*M_3C*UCStop(1,3) - 2*Cd(10140)*UCStopC(1,3)

	Cd(122) = 2*Cd(121) + MGlpTmSt2(1)*UStop2(1,2)

	Cd(1245) = MStop2(1) - MStop2(2)*UStop2(1,2)

	Cd(1246) = MStop2(1) + 2*MStop2(3)*UStop2(1,2)

	Cd(1251) = 2 + 5*UStop2(1,2)**2

        Cd(12609) = MUEC**2*UCStopC(1,3)**3 + 
     &   MUE**2*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)

	Cd(12611) = MUE**2*UCStop(1,3)**2 + 2*MUE2*Cd(9948)*UStop2(1,1)

        Cd(12612) = Cd(12611)*UStop(2,1) - 
     &   MUEC**2*UCStopC(1,3)*UStop(2,2)*UStop2(1,2)

	Cd(12617) = MGlpT2 - 2*MStop2(1)

	Cd(12618) = 2*MStop2(2) - MGlpT2*UStop2(1,2)

	Cd(12619) = MGlpT2 - 2*MStop2(2)

        Cd(123) = 2*(M_3*Cd(122)*UCStopC(1,3) + 
     &     M_3C*MGlpTmSt2(1)*UCStop(1,3)*UStop2(1,1))

        Cd(12623) = MUEC*Cd(9948)*Cd(12071)*UUStopC(1,1) + 
     &   MUE*Cd(12075)*UStop2(1,1)*UUStopC(1,2)

	Cd(12626) = Cd(12137) - AtC*UCStop(1,3)

	Cd(12632) = 2*MTy*UUStop(1,1) + At*UStop2(1,2)*UUStop(1,3)

	Cd(12633) = MTy*Cd(12632) - 2*MTy2*UStop2(1,1)*UUStop(1,1)

	Cd(12648) = MUE*UCStop(1,3)*UStop2(1,2) + Cd(12148)*UUStopC(1,3)

	Cd(12666) = 2*Cd(9948) + UStop2(1,1)

	Cd(12679) = -4*Cd(9955)*UStop2(1,1) + UStop2(1,2)

	Cd(1262) = 2*MStop2(1)**2 + MStop2(3)**2

	Cd(12700) = AtC + MUE*TB

	Cd(12701) = At + MUEC*TB

        Cd(12702) = MStop2(3)**2 - 
     &   Cd(12700)*Cd(12701)*MStop2(2)*UStop2(1,2)

        Cd(12703) = -2*Cd(12700)**2*MStop2(2)*UCStop(1,3)**2 + 
     &   2*Cd(12702)*UStop2(1,2)

	Cd(12704) = Cd(12701)*UCStopC(1,3)**2 + Cd(12700)*UStop2(1,1)**2

        Cd(12705) = -2*Cd(12701)*Cd(12704)*MStop2(2) + 
     &   2*MStop2(3)**2*UStop2(1,1)

	Cd(12706) = Cd(12703) + Cd(12705)

	Cd(12707) = Cd(12700)*UCStop(3,3) - Cd(12701)*UCStop(3,4)

	Cd(1264) = 2*MStop2(1)**2 - MStop2(3)**2

	Cd(12708) = Cd(12701)*UCStopC(3,3) - Cd(12700)*UCStopC(3,4)

	Cd(12709) = Cd(12706) + 2*Cd(12707)*Cd(12708)*MStop2(2)

	Cd(12710) = MStop2(3) - Cd(12700)*Cd(12701)*UStop2(1,1)

        Cd(12711) = Cd(12710)*UStop(1,1) - 
     &   Cd(12701)**2*UCStopC(1,3)*UStop(1,2)

        Cd(12712) = Cd(12700)*Cd(12701)*MStop2(3)*UStop2(1,1)*
     &    UStop2(1,2) + Cd(12711)*MStop2(1)*UStopC(1,1)

	Cd(12713) = MStop2(3) - Cd(12700)*Cd(12701)*UStop2(1,2)

        Cd(12714) = Cd(12700)**2*UCStop(1,3)*UStop(1,1) - 
     &   Cd(12713)*UStop(1,2)

        Cd(12715) = Cd(12700)**2*MStop2(3)*UCStop(1,3)**2 + 
     &   2*Cd(12714)*MStop2(1)*UStopC(1,2)

	Cd(12716) = 4*Cd(12712) - 2*Cd(12715)

        Cd(12717) = (Cd(12716) - 
     &     2*Cd(12701)**2*MStop2(3)*UCStopC(1,3)**2)/MStop2(3)

	Cd(1265) = 12*MStop2(3)**2 - 4*Cd(1264)*UStop2(1,2)**2

        Cd(12718) = (2*MStop2(1)*
     &     (-2*Cd(12706)*Cr(13)*MStop2(1) + Cd(12717)*MStop2(3)**2)
     &     )/MStop2(3)**2

        Cd(12719) = 3/(64.D0*Pi)*
     &   (Alfa2L*(-8*Cd(12707)*Cd(12708)*Cr(34)*MStop2(1)**2*
     &         MStop2(2) + Cd(12718)*MStop2(3)**2))/
     &    (SW2*TB2*MStop2(1)**2*MStop2(3)**2)

	Cd(12724) = A0delStop(2)*MStop2(1) - A0delStop(1)*MStop2(2)

	Cd(12726) = Cr(13) - Cr(34)

	Cd(1266) = Cd(1265) - 4*Cd(1262)*UStop2(1,1)**2

        Cd(12731) = MStop2(3)**2 - 
     &   Cd(12700)*Cd(12701)*MStop2(1)*UStop2(1,1)

        Cd(12732) = -8*Cd(12700)**2*MStop2(1)*UCStop(1,3)**2 + 
     &   8*Cd(12731)*UStop2(1,1)

	Cd(12733) = Cd(12701)*UCStopC(1,3)**2 + Cd(12700)*UStop2(1,2)**2

        Cd(12734) = 4*Cd(12701)*Cd(12733)*MStop2(1) - 
     &   4*MStop2(3)**2*UStop2(1,2)

	Cd(12735) = Cd(12732) - 2*Cd(12734)

	Cd(12736) = MStop2(3) + Cd(12700)*Cd(12701)*UStop2(1,1)

	Cd(12738) = MStop2(2) + Cd(12700)*Cd(12701)*UStop2(1,1)

        Cd(12739) = Cd(12701)*Cd(12733)*MStop2(2) + 
     &   Cd(12738)*MStop2(3)*UStop2(1,2)

	Cd(12744) = Cd(12735) + 8*Cd(12707)*Cd(12708)*MStop2(1)

	Cd(1268) = 6*MStop2(2)**2 + Cd(1009)*MStop2(3)

	Cd(12749) = MStop2(2) - MStop2(4)*UStop2(1,1)*UStop2(1,2)

        Cd(12753) = -(Cd(12701)*UCStopC(1,3)**2) + 
     &   Cd(12700)*UStop2(1,1)*UStop2(1,2)

        Cd(12754) = UCStopC(1,3)*UStop2(1,1)*UUStop(2,1) - 
     &   UCStopC(3,1)*UStop2(1,2)*UUStop(3,3)

        Cd(12755) = Cd(12700)*UCStop(1,3)*UUStop(3,1) + 
     &   Cd(12701)*UStop2(1,2)*UUStop(3,3)

        Cd(12756) = Cd(12753)*UStop2(1,1)*UStop2(1,2) - 
     &   Cd(12700)*Cd(12754)*UUStopC(2,3)

	Cd(12762) = MStop2(2)*UStop2(1,1) - MStop2(1)*UStop2(1,2)

	Cd(12763) = 2*MStop2(1)*UStop2(1,1) - MStop2(2)*UStop2(1,2)

	Cd(12764) = 2*UStop(1,1) + Cr(58)*UStopC(2,2)

	Cd(12765) = UStop2(1,1) + Cd(12764)*UStopC(1,1)

        Cd(12767) = 2*MStop2(2)*UCStopC(1,3)*UCStopC(3,1) - 
     &   MStop2(1)*UCStopC(3,3)*UStop2(1,2)

	Cd(1270) = 4 - 8*UStop2(1,2)

        Cd(12768) = 2*UStop2(1,1)*UStopC(2,2)*UUStopC(1,1) - 
     &   UCStopC(3,1)*UStopC(1,2)*UUStopC(2,1)

        Cd(12770) = 4*MStop2(2)*
     &   (2*Cr(28)*MStop2(2)*UUStop(2,3) - 
     &     Cd(11656)*MStop2(1)*UUStopC(1,3))

	Cd(12771) = MStop2(2)*UUStop(3,3) + MStop2(1)*UUStop(3,4)

        Cd(12772) = Cd(12771)*UStop2(1,2)**2 + 
     &   Cd(12762)*UStop2(1,1)*UUStop(3,3)

        Cd(12773) = Cd(12701)*Cd(12772) + 
     &   2*Cd(12700)*MStop2(4)*UCStop(1,3)*UStop2(1,2)*UUStop(3,1)

	Cd(12775) = At + 2*MUEC*TB

	Cd(1271) = Cd(1270)*MStop2(3)**2 + 8*MStop2(1)**2*UStop2(1,2)

        Cd(12778) = MStop2(3)*UCStop(3,1)*UStop(2,2) + 
     &   4*MStop2(2)*UStop(1,2)*UStop2(1,1)

        Cd(12779) = MStop2(3)*UCStopC(1,3)*UCStopC(3,1) + 
     &   MStop2(2)*UCStopC(3,3)*UStop2(1,2)

	Cd(12781) = UStop2(1,2)*UUStop(3,3) - UStop2(1,1)*UUStop(3,4)

        Cd(12782) = Cd(12781)*MStop2(3)*UCStopC(3,1) - 
     &   Cd(11656)*MStop2(2)*UUStopC(1,3)

	Cd(12783) = MStop2(3)*UUStop(3,3) - MStop2(2)*UUStop(3,4)

        Cd(12784) = Cd(12783)*UStop2(1,2) + 
     &   MStop2(2)*UStop2(1,1)*UUStop(3,3)

        Cd(12785) = Cd(12701)*Cd(12784) - 
     &   Cd(12700)*MStop2(4)*UCStop(1,3)*UUStop(3,1)

	Cd(1272) = 2*MStop2(1) + MStop2(3)

        Cd(12796) = Cd(12700)**2*UCStop(1,3)**2 + 
     &   Cd(12701)**2*UCStopC(1,3)**2

        Cd(12802) = 2*UCStopC(1,3)*UCStopC(3,1) - 
     &   UCStopC(3,3)*UStop2(1,2)

        Cd(12803) = Cd(12701)*Cd(12802)*UStop(1,2) - 
     &   2*Cd(12700)*UStop(2,1)*UStop2(1,1)*UStop2(1,2)

        Cd(12804) = Cd(12700)*UCStop(1,3)*UStop(1,1) - 
     &   Cd(12701)*UStop(1,2)*UStop2(1,1)

	Cd(12807) = 2*MStop2(1)*UStop2(1,1) - MStop2(3)*UStop2(1,2)

	Cd(1274) = 3 + 2*UStop2(1,1)

        Cd(12808) = Cd(12700)*MStop2(4)*UCStop(1,3)**2 + 
     &   Cd(12701)*Cd(12807)*UStop2(1,1)

	Cd(12809) = UStop(1,1) + Cr(58)*UStopC(2,2)

	Cd(12810) = UStop2(1,1)**2 + 2*UStop2(1,2)**2

        Cd(12812) = -(Cd(12810)*MStop2(1)) + 
     &   MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(12813) = Cd(12810)*MStop2(1)*UStop2(1,1)*
     &    UStop2(1,2) + Cd(12812)*UUStop(3,4)*UUStopC(3,3)

        Cd(12819) = 2*MStop2(2)*UStop2(1,1)**3 - 
     &   Cr(28)*MStop2(1)*UStop2(1,2)

        Cd(12820) = Cd(12701)*Cd(12819) + 
     &   2*Cd(12700)*MStop2(4)*UCStop(1,3)**2*UStop2(1,1)

        Cd(12821) = -(MStop2(1)*UCStop(3,1)*UStop(2,1)) + 
     &   2*MStop2(2)*UStop(1,1)*UStop2(1,1)

        Cd(12822) = 3*MStop2(2)*UCStopC(1,3)*UCStopC(3,1) + 
     &   Cd(12114)*MStop2(1)*UCStopC(3,3)

	Cd(12823) = UStop2(1,1)*UUStop(3,3) - UStop2(1,2)*UUStop(3,4)

        Cd(12824) = Cd(12823)*MStop2(1)*UStop2(1,1) + 
     &   Cd(12810)*MStop2(2)*UUStop(3,4)

        Cd(12827) = Cd(12820) + 
     &   Cd(12700)*Cd(12821)*UCStop(1,3)*UStopC(1,2)

	Cd(12832) = Cd(12700)*UCStop(3,3) - 2*Cd(12701)*UCStop(3,4)

        Cd(12833) = -2*Cd(12701)*UCStopC(1,3)**2 + 
     &   2*Cd(12700)*UStop2(1,1)*UStop2(1,2)

	Cd(12834) = Cd(12832)*UCStopC(3,3) + Cd(12700)*UStop2(1,2)**2

	Cd(12835) = 2*Cd(12834)*MStop2(1) - Cd(12833)*MStop2(3)

	Cd(12836) = Cd(12700)*UCStop(1,3) - Cd(12701)*UCStopC(1,3)

        Cd(12837) = 2*(Cd(12834)*MStop2(2) + 
     &     Cd(12836)*MStop2(3)*UCStopC(1,3))

	Cd(12838) = Cd(12835) - Cd(12837)

        Cd(12839) = 1/2.D0*(A0delStop(1)*Cd(12835)) + 
     &   Cd(12753)*Cr(13)*MStop2(1)*MStop2(3)

        Cd(12840) = 1/2.D0*(A0delStop(2)*Cd(12837)) - 
     &   Cd(12836)*Cr(34)*MStop2(2)*MStop2(3)*UCStopC(1,3)

        Cd(12841) = 3/(16.D0*Pi)*
     &   (Alfa2L*(-(Cd(12840)*MStop2(1)) + Cd(12839)*MStop2(2)))/
     &    (SW2*TB2*MStop2(1)*MStop2(2)*MStop2(3))

	Cd(12842) = Cd(12835)*Cr(13) + 2*Cd(12753)*MStop2(3)

        Cd(12843) = Cd(12837)*Cr(34) - 
     &   2*Cd(12836)*MStop2(3)*UCStopC(1,3)

	Cd(12844) = Cd(12842) - Cd(12843)

	Cd(12845) = 2*Cd(12700)*UCStop(3,3) - Cd(12701)*UCStop(3,4)

        Cd(12846) = -(Cd(12700)*UCStop(1,3)**2) + 
     &   Cd(12701)*UStop2(1,1)*UStop2(1,2)

        Cd(12847) = -2*Cd(12845)*UCStopC(3,4) + 
     &   2*Cd(12701)*UStop2(1,1)**2

	Cd(12848) = -(Cd(12847)*MStop2(1)) + 2*Cd(12846)*MStop2(3)

        Cd(12849) = Cd(12700)*MStop2(3)*UCStop(1,3)**2 + 
     &   Cd(12845)*MStop2(2)*UCStopC(3,4)

        Cd(12852) = 1/2.D0*(A0delStop(1)*Cd(12848)) - 
     &   Cd(12846)*Cr(13)*MStop2(1)*MStop2(3)

	Cd(12855) = Cd(12848)*Cr(13) - 2*Cd(12846)*MStop2(3)

	Cd(12862) = -2*A0delT + MTy2*Cr(15)

	Cd(12867) = -2*MTy2 + 4*Cr(66)

	Cd(1280) = MStop2(1)**2 + MStop2(2)*MStop2(3)

	Cd(12872) = UStop2(1,1)**2 + Cd(11106)*UUStopC(3,3)

	Cd(12873) = Cd(12872) + Cr(91)*UUStopC(3,3)

	Cd(1281) = MStop2(2) - 3*MStop2(3)

	Cd(12878) = UStop2(1,1)*UStop2(1,2) - Cd(10925)*UUStopC(3,3)

	Cd(12879) = Cd(12878)*MStop2(1) - Cd(12872)*MStop2(2)

	Cd(12880) = Cd(12879) - Cr(91)*MStop2(2)*UUStopC(3,3)

        Cd(12882) = Cd(12701)**2*UCStopC(1,3)*UUStopC(3,1) - 
     &   Cd(12700)**2*UCStop(1,3)*UUStopC(3,2)

	Cd(1282) = Cd(1280) + Cd(1281)*MStop2(2)*UStop2(1,2)

	Cd(12888) = Cd(10904) - 2*UStop2(1,2)

	Cd(12889) = MStop2(4)*UStop2(1,1) - MStop2(1)*UStop2(1,2)

        Cd(12890) = 2*Cd(12889)*MStop2(2)*UStop(1,1) - 
     &   2*MStop2(1)**2*UCStop(3,1)*UStop(2,1)

        Cd(12892) = UStop(1,2)*UStop2(1,2)**2 - 
     &   Cd(11106)*UCStop(3,1)*UStopC(1,1)

	Cd(12895) = Cd(10925)*UCStop(1,3) + UStop2(1,2)*UUStop(3,1)

	Cd(12902) = 3*UUStop(3,3) + 2*UUStop(3,4)

        Cd(12903) = Cd(12902)*MStop2(1)*UStop2(1,2) + 
     &   2*MStop2(2)*UUStop(3,3)

	Cd(12905) = UStop2(1,1)*UStop2(1,2) - 2*Cr(91)*UUStopC(3,3)

	Cd(12908) = 2*MStop2(4)*(MStop2(1) - MStop2(4)*UStop2(1,1))

	Cd(11) = GS2L*(1 - 2*UStop2(1,1))

        Cd(12925) = MStop2(1)*UCStop(3,1)*UStop(2,1) + 
     &   MStop2(2)*UStop(1,1)*UStop2(1,1)

	Cd(12926) = MStop2(1)*UStop2(1,1) + MStop2(2)*UStop2(1,2)

	Cd(12927) = MStop2(3)*UUStop(3,3) - 2*MStop2(1)*UUStop(3,4)

	Cd(12928) = 2*MStop2(2)*UUStop(3,3) + MStop2(3)*UUStop(3,4)

        Cd(12929) = Cd(12700)*Cd(12928)*UCStop(1,3) - 
     &   Cd(12701)*Cd(12927)*UCStopC(1,3)

	Cd(12931) = Cd(12928)*UStop2(1,1) - 2*MStop2(2)*UUStop(3,3)

	Cd(12933) = MStop2(2)*UStop2(1,1) + MStop2(1)*UStop2(1,2)

	Cd(12934) = MStop2(1) + MStop2(3)*UStop2(1,1)

	Cd(12935) = MStop2(3)*UUStop(3,3) + 2*MStop2(2)*UUStop(3,4)

	Cd(12936) = 2*MStop2(1)*UUStop(3,3) - MStop2(3)*UUStop(3,4)

	Cd(12937) = Cd(12935)*UStop2(1,1) + Cd(12936)*UStop2(1,2)

	Cd(12938) = Cd(12936)*UUStopC(3,3) - Cd(12935)*UUStopC(3,4)

        Cd(12939) = Cd(12700)*Cd(12938)*UStop2(1,1)*
     &    UStop2(1,2) - Cd(12701)*Cd(12937)*UCStopC(1,3)*
     &    UUStopC(3,1)

        Cd(12941) = Cd(12939) + 
     &   Cd(12701)*Cd(12926)*Cr(91)*UCStopC(1,3)*UUStopC(3,1)

        Cd(12955) = Cd(12700)*MGlpTmSt2(1)*UUStop(1,1) + 
     &   Cd(12701)*MGlpTmSt2(2)*UUStop(1,2)

        Cd(12956) = 2*Cd(12700)*MStop2(3)*UStop2(1,1)*
     &    UStop2(1,2) + Cd(12955)*UUStopC(1,1)

        Cd(12957) = Cd(12700)*MGlpTmSt2(2)*UUStop(1,1) + 
     &   Cd(12701)*MGlpTmSt2(1)*UUStop(1,2)

        Cd(12958) = Cd(12701)*Cd(12956) + 
     &   Cd(12700)*Cd(12957)*UUStopC(1,2)

	Cd(12960) = 2*Cd(12700)*UUStop(1,1) + 3*Cd(12701)*UUStop(1,2)

	Cd(12961) = MStop2(3)*UStop(1,2) + MTy*Cd(12960)*UStopC(1,1)

	Cd(12962) = Cd(12700)*UUStop(1,1) + 2*Cd(12701)*UUStop(1,2)

	Cd(12963) = MTy*Cd(12962)*Cr(75) + MStop2(3)*UUStop(1,3)

	Cd(12966) = 3*Cd(12700)*UUStop(1,1) + 2*Cd(12701)*UUStop(1,2)

	Cd(12967) = MStop2(3)*UStop(1,1) + MTy*Cd(12966)*UStopC(1,2)

	Cd(12968) = 2*Cd(12700)*UUStop(1,1) + Cd(12701)*UUStop(1,2)

	Cd(12969) = MTy*Cd(12968)*Cr(9) + MStop2(3)*UUStop(1,3)

        Cd(12974) = Cd(12700)*MStop2(3)*UUStop(1,1) + 
     &   Cd(12701)*MStop2(2)*UUStop(1,2)

	Cd(12976) = 2*Cd(12700)*UUStop(3,1) + Cd(12701)*UUStop(3,2)

	Cd(12979) = Cd(12700)*UUStop(1,1) + Cd(12701)*UUStop(1,2)

	Cd(12982) = MUE*TB*UUStop(2,1) - Cd(12701)*UUStop(2,2)

	Cd(12984) = Cd(12700)*UUStop(2,1) - 2*Cd(12701)*UUStop(2,2)

	Cd(12994) = 2*Cd(12700)*UCStop(1,3) - Cd(12701)*UCStopC(1,3)

        Cd(12995) = -(Cd(12701)*Cd(12994)*UCStopC(1,3)**2) + 
     &   Cd(12700)**2*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)

        Cd(12996) = Cd(12700)*UCStop(1,3)*UStop(2,1) + 
     &   2*Cd(12701)*UStop(2,2)*UStop2(1,2)

	Cd(12997) = 2*Cd(12701)*UStop2(1,2)**2 + Cd(12996)*UStopC(2,2)

        Cd(13001) = -(MStop2(3)*UCStopC(1,3)) + 
     &   2*MTy*Cd(12700)*Cr(75)*UStop2(1,2)

	Cd(13002) = Cd(12700)*UUStop(3,1) + 2*Cd(12701)*UUStop(3,2)

        Cd(13003) = MStop2(3)*UStop(1,1)*UStop2(1,2) - 
     &   MTy*Cd(13002)*UStop2(1,1)*UStopC(2,2)

	Cd(13004) = Cd(13001) - 3*MTy*Cd(12701)*UCStopC(1,3)**2

        Cd(13007) = MStop2(3)*UStop(1,1)*UStop2(1,2) + 
     &   MTy*Cd(12976)*UStop2(1,1)*UStopC(2,2)

        Cd(13008) = 2*MTy*Cd(12700)*Cr(75)*UCStop(1,3)*
     &    UStop2(1,2) + Cd(13007)*UStopC(1,1)

	Cd(13009) = MStop2(3) + 3*MTy*Cd(12700)*UCStop(1,3)

	Cd(13014) = MStop2(3) + MGlpTmSt2(2)*UStop2(1,2)

	Cd(13015) = Cd(12700)*UUStop(1,1) - 3*Cd(12701)*UUStop(1,2)

	Cd(13016) = MGlpTmSt2(1) + MStop2(3)

        Cd(13017) = 2*Cd(12700)*MGlpTmSt2(2)*UUStop(3,1) + 
     &   Cd(12701)*Cd(13016)*UUStop(3,2)

        Cd(13018) = Cd(13015)*MGlpTmSt2(2)*UStopC(1,2) + 
     &   Cd(13017)*UStopC(2,2)

	Cd(13019) = 3*UStop(1,1)*UStop2(1,2) + Cr(91)*UStopC(2,2)

        Cd(13020) = 2*Cd(12700)*MStop2(3)*UUStop(2,1) + 
     &   At*MGlpTmSt2(2)*UUStop(2,2)

        Cd(13025) = Cd(12700)*Cd(13014)*UCStop(1,3)**2 - 
     &   Cd(13018)*UStop2(1,1)*UStopC(1,2)

	Cd(13028) = Cr(85)*UStop2(1,1) + UStop2(1,2)**2

	Cd(13030) = 4*Cd(12700)*UUStop(3,1) + 2*Cd(12701)*UUStop(3,2)

	Cd(1297) = Cd(1265) - 4*Cd(1264)*UStop2(1,1)**2

	Cd(13044) = -2*MTy*Cr(5)*MGlpTmSt2(1) + MGlpTmSt4(1)

	Cd(13050) = MTy*Cd(12968)*Cr(75) - MStop2(3)*UUStop(1,3)

	Cd(1299) = MStop2(4)*(6 + 4*UStop2(1,2)**2)

	Cd(1300) = 3 + 2*UStop2(1,1)*UStop2(1,2)**2

	Cd(1301) = Cd(1300)*MStop2(2) - 3*MStop2(1)*UStop2(1,2)

	Cd(13077) = MUE*TB*UStop(1,1) - MTy*UStop(1,2)

	Cd(13078) = MUE*TB*UCStop(1,3) - Cd(13077)*UStopC(1,2)

	Cd(1302) = 4*MTy2*Cd(1301) + 3*MStop2(2)*MStop2(3)

	Cd(13106) = 5*MStop2(2) - 2*Cr(28)*MStop2(4)

	Cd(13107) = ((Cd(13106) - 7*MStop2(1))*UCStop(1,3))/MStop2(1)**2

	Cd(13111) = 5*MStop2(3) + 4*MStop2(4)*UStop2(1,1)*UStop2(1,2)

        Cd(13112) = 1/Pi**4*
     &   (Cr(28)*MStop2(1) + Cd(13111)*UStop2(1,1)*UStop2(1,2))/
     &    MStop2(1)**2

	Cd(13113) = (AtC + MUE*TB)*MStop2(4)

        Cd(13114) = Cd(13113)*UUStop(1,1) + 
     &   Cd(12701)*MStop2(4)*UUStop(1,2)

        Cd(13115) = 1/(256.D0*Pi**4)*
     &   (Pi**4*(Cd(12700)*Cd(13112)*MStop2(1)**2) - 
     &      2*Cd(13114)*UStop2(1,1)*UStop2(1,2)*UUStopC(1,1))/
     &    MStop2(1)**2

	Cd(1305) = 4*(3*MStop2(3) + Cr(75)**2*MStop2(4))

	Cd(13120) = 2*MStop2(2) - MStop2(4)

	Cd(13121) = Cd(13120)*Cr(28) - 2*MStop2(2)

        Cd(13122) = 1/2.D0*
     &   ((Cd(13121) - 6*MStop2(1))*UCStop(1,3))/MStop2(1)**2

	Cd(13125) = MStop2(3) + Cd(13120)*UStop2(1,1)*UStop2(1,2)

	Cd(1306) = Cd(1302) - MTy2*Cd(1299)*UStop2(1,1)

        Cd(13126) = 1/Pi**4*
     &   (-2*Cr(28)*MStop2(1) + Cd(13125)*UStop2(1,1)*UStop2(1,2))/
     &    MStop2(1)**2

	Cd(13127) = AtC*Cd(13120) + MUE*TB*MStop2(3)

        Cd(13128) = Cd(13127)*UUStop(1,1) + 
     &   Cd(12701)*MStop2(3)*UUStop(1,2)

        Cd(13129) = 1/(512.D0*Pi**4)*
     &   (2*Pi**4*(Cd(12700)*Cd(13126)*MStop2(1)**2) - 
     &      Cd(13128)*UStop2(1,1)*UStop2(1,2)*UUStopC(1,1))/
     &    MStop2(1)**2

        Cd(13130) = Cd(12700)*MStop2(3)*UUStop(1,1) + 
     &   Cd(12701)*Cd(13120)*UUStop(1,2)

        Cd(13146) = UStop2(1,2)*
     &   (UStop2(1,2)**2 + Cd(11106)*UUStopC(3,3))

        Cd(13154) = Cd(12700)*UCStop(1,3)**2*UStop(1,1) + 
     &   Cd(12701)*UStop2(1,1)*UStopC(2,2)*UUStop(3,4)

        Cd(13155) = 1/4.D0*
     &   (MStop2(3)*(2*Cd(12701)*Cd(13120)*UStop(2,2)*
     &        UUStop(1,2) + 
     &       Cd(12700)*MStop2(3)*UStop(1,2)*UUStop(3,1)))

        Cd(13156) = 1/2.D0*
     &    (Cd(13154)*MStop2(3)**2*UStop2(1,2)*UStopC(1,1)) - 
     &   Cd(13155)*UStop2(1,1)*UStopC(2,2)*UUStopC(1,2)

	Cd(13157) = UUStop(1,2)*UUStopC(1,1) + UUStop(2,2)*UUStopC(2,1)

        Cd(13158) = Cd(12701)*Cd(13157)*MStop2(3) - 
     &   4*Cd(12700)*Cd(13120)*UStop2(1,1)*UStop2(1,2)

	Cd(13159) = UCStopC(1,3)*UStop2(1,1) - UUStop(3,3)*UUStopC(3,1)

        Cd(13160) = Cd(12701)*Cd(13158) + 
     &   Cd(12700)**2*Cd(13120)*UCStop(1,3)**2

        Cd(13161) = Cd(12701)**2*Cd(13159)*MStop2(3)*
     &    UCStopC(1,3) + Cd(13160)*UStop2(1,1)

        Cd(13162) = 1/(512.D0*Pi**4)*
     &   (4*Cd(12700)*Cd(13156) + Cd(13161)*MStop2(3)*UStop2(1,2))/
     &    MStop2(3)

	Cd(13170) = 2*MStop2(2)*UStop2(1,1) + MStop2(1)*UStop2(1,2)

        Cd(13173) = Cd(12700)*MStop2(3)*UUStop(1,1) + 
     &   Cd(12701)*Cr(1)*UUStop(1,2)

        Cd(13174) = Cd(12700)*Cr(91)*UStop(2,1) + 
     &   Cd(12701)*UStop(1,2)*UUStop(2,2)

	Cd(13177) = MStop2(2)*(MStop2(2) - 3*MStop2(3)) + MStop2(3)**2

        Cd(13178) = Cd(12700)*Cr(91)*UCStop(1,3) - 
     &   Cd(12701)*UStop2(1,1)*UUStop(3,2)

        Cd(13179) = Cd(11106)*MStop2(1)*MStop2(2) - 
     &   MStop2(3)**2*UUStop(3,4)

	Cd(13184) = 7*MStop2(2) - 5*MStop2(3)

	Cd(13185) = Cd(12700)*UUStop(1,1) + 3*Cd(12701)*UUStop(1,2)

	Cd(1312) = 3*MStop2(3) + Cr(85)*MStop2(2)*UStop2(1,2)

	Cd(13187) = 2*Cd(12700)*UUStop(2,1) + Cd(12701)*UUStop(2,2)

	Cd(13188) = Cd(12700)*UUStop(3,1) + Cd(12701)*UUStop(3,2)

        Cd(13189) = 2*Cd(12700)*MStop2(1)**2*UUStop(2,1) - 
     &   Cd(12701)*MStop2(3)*MStop2(4)*UUStop(2,2)

        Cd(13190) = Cd(13189)*UStop(1,1) + 
     &   4*Cd(12701)*MStop2(1)*MStop2(2)*UStop(2,1)*UUStop(3,2)

        Cd(13191) = Cd(12701)*UUStop(2,2)*UUStop(3,3) + 
     &   Cd(12700)*UUStop(2,1)*UUStop(3,4)

        Cd(13193) = Cd(13188)*UCStopC(3,1)*UStop2(1,2) - 
     &   Cd(13191)*UUStopC(3,3)

        Cd(13194) = Cd(13193)*MStop2(3) + 
     &   Cd(13187)*MStop2(2)*UStop2(1,2)**2

        Cd(13195) = 2*Cd(13185)*MStop2(2)*UStop2(1,1)**2*
     &    UUStopC(1,1) + Cd(13194)*UUStopC(2,1)

	Cd(1313) = UStop2(1,1)**3 + UStop2(1,2)**3

        Cd(13201) = -(MStop2(1)*UCStop(3,1)*UStop(2,1)) + 
     &   MStop2(3)*UStop(1,1)*UStop2(1,1)

	Cd(13202) = Cr(28)*MStop2(2)*UStop(1,1) + Cd(13201)*UStop2(1,2)

        Cd(13203) = MStop2(2)*UStop(1,2)*UStop2(1,2) - 
     &   MStop2(3)*UStopC(2,1)*UUStop(3,3)

        Cd(13208) = 1/2.D0*
     &   (MStop2(3)*(Cd(12700)*MStop2(1)*UUStop(1,1) + 
     &       2*Cd(12701)*MStop2(2)*UUStop(1,2)))

        Cd(13209) = 2*Cd(12701)*MStop2(2)*UStop(2,2)*
     &    UUStop(1,2) + Cd(12700)*MStop2(3)*UStop(1,2)*UUStop(3,1)

        Cd(13210) = 1/2.D0*
     &   (MStop2(3)*(Cd(13178)*MStop2(1) + 
     &       2*Cd(13209)*UStopC(1,2)))

	Cd(13211) = UUStop(3,3) - 4*UUStop(3,4)

	Cd(13212) = 1/2.D0*(UUStop(3,3) - 2*UUStop(3,4))

        Cd(13217) = 2*Cd(12700)*MStop2(2)*UUStop(1,1) + 
     &   Cd(12701)*MStop2(4)*UUStop(1,2)

        Cd(13218) = Cd(12700)*MStop2(1)*UUStop(2,1) + 
     &   Cd(12701)*MStop2(2)*UUStop(2,2)

	Cd(13219) = Cd(13217)*UUStopC(1,1) - Cd(13218)*UUStopC(2,1)

        Cd(13220) = -(1/2.D0*(Cd(13193)*MStop2(1)*UUStopC(2,1))) + 
     &   Cd(13219)*UStop2(1,1)*UStop2(1,2)

        Cd(13221) = 1/(128.D0*Pi**4)*
     &   (Cd(13220) - 2*Cd(12700)*MStop2(1)*UStop2(1,1)**2*
     &       UStop2(1,2)**2)/(MStop2(1)*MStop2(2))

	Cd(1316) = 4*MStop2(2)*MStop2(3) - 4*Cd(1264)*UStop2(1,2)

        Cd(13234) = UStop2(1,1)**2*UStopC(1,1) + 
     &   UCStopC(3,1)*UStop2(1,2)*UStopC(2,1)

	Cd(13235) = 3*Cd(12700)*UUStop(1,1) + Cd(12701)*UUStop(1,2)

	Cd(13236) = MStop2(3)*UStop(1,1) + MTy*Cd(13235)*UStopC(1,2)

	Cd(13237) = MStop2(3)*UStop(1,1) + 2*MTy*Cd(12979)*UStopC(1,2)

	Cd(13238) = 4*Cd(12700)*UUStop(1,1) + Cd(12701)*UUStop(1,2)

        Cd(13239) = Cd(13238)*UStop2(1,2)*UStopC(1,2) + 
     &   Cd(12976)*UStop2(1,1)*UStopC(2,2)

	Cd(13240) = MTy*Cd(13239) - Cd(13237)*UStop2(1,1)

	Cd(13241) = MStop2(3)*UStop2(1,2)**2 - Cd(13240)*UStopC(1,1)

        Cd(13242) = Cd(13241) + 
     &   2*MTy*Cd(12700)*UCStop(1,3)*UStop2(1,2)**2

	Cd(13243) = Cd(13242) - MTy*Cd(12701)*Cd(13234)*UStop(1,2)

        Cd(13244) = Cd(12701)*Cd(13243) + 
     &   2*Cd(12700)*Cd(13236)*UCStop(1,3)*UStopC(1,2)

	Cd(13245) = 2*MStop2(3)*UStop(1,2) + MTy*Cd(12960)*UStopC(1,1)

        Cd(13247) = Cd(13002)*UStop2(1,1)*UStopC(2,2) + 
     &   2*Cd(12979)*UStop(1,2)*UUStopC(1,2)

	Cd(13248) = Cr(28)*MStop2(3) + MTy*Cd(13247)*UStopC(1,1)

        Cd(13251) = Cd(12700)*MStop2(3)*UUStop(1,1) - 
     &   Cd(12701)*MGlpTmSt2(2)*UUStop(1,2)

	Cd(13252) = 2*MGlpTmSt2(1) + MStop2(3)

	Cd(13255) = UUStop(1,1) + UUStop(2,1)

	Cd(13257) = UUStop(1,2) + UUStop(2,2)

	Cd(1320) = Cd(1316) + 12*MStop2(1)*MStop2(2)*UStop2(1,1)

	Cd(13266) = UStop2(1,1)*UStop2(1,2) - Cr(91)*UUStopC(3,3)

        Cd(13276) = Cd(12700)*MStop2(4)*UUStop(1,1) + 
     &   2*Cd(12701)*MStop2(1)*UUStop(1,2)

        Cd(13277) = Cd(12700)*MStop2(4)*UUStop(1,1) + 
     &   2*Cd(12701)*Cr(26)*UUStop(1,2)

	Cd(13278) = Cd(13277)*UStop2(1,1) + Cd(13276)*UStop2(1,2)

        Cd(13279) = 2*Cd(12700)*MStop2(1)*UUStop(1,1) + 
     &   Cd(12701)*MStop2(4)*UUStop(1,2)

	Cd(13280) = (AtC + MUE*TB)*Cr(26)

        Cd(13281) = 2*Cd(13280)*UUStop(1,1) + 
     &   Cd(12701)*MStop2(4)*UUStop(1,2)

	Cd(13282) = Cd(13279)*UStop2(1,1) + Cd(13281)*UStop2(1,2)

        Cd(13283) = Cd(12701)*Cd(13282)*UUStopC(1,1) + 
     &   Cd(12700)*Cd(13278)*UUStopC(1,2)

        Cd(13284) = 1/(256.D0*Pi**4)*
     &   (Cd(13283) - 2*MStop2(1)*MStop2(3))/MStop2(1)**2

        Cd(13295) = 2*Cd(12700)*MStop2(2)*UUStop(3,1) + 
     &   Cd(12701)*MStop2(1)*UUStop(3,2)

        Cd(13296) = Cd(13015)*MStop2(2)*UStopC(1,2) + 
     &   Cd(13295)*UStopC(2,2)

        Cd(13297) = -(Cd(13296)*UStop2(1,1)) + 
     &   Cd(12700)*MStop2(2)*UCStop(1,3)*UStop(1,1)*UStop2(1,2)

        Cd(13298) = Cd(12700)*MStop2(3)*UUStop(2,1) - 
     &   At*MStop2(2)*UUStop(2,2)

        Cd(13303) = Cd(12796)*MStop2(3) - 
     &   2*Cd(12700)*Cd(13297)*UStopC(1,2)

	Cd(13308) = Cd(12700)*UUStop(1,1) - 2*Cd(12701)*UUStop(1,2)

        Cd(13309) = Cd(12700)*MStop2(3)*UCStop(1,3)*UStop(1,1)*
     &    UStop2(1,2) + Cd(13120)*Cd(13308)*UStop2(1,1)*UStopC(1,2)

        Cd(13310) = 2*AtC*Cd(13120)*UUStop(1,1) - 
     &   Cd(12701)*MStop2(3)*UUStop(1,2)

        Cd(13311) = 1/2.D0*
     &   (MStop2(3)*(Cd(13310) + 2*MUE*TB*MStop2(3)*UUStop(1,1)))

        Cd(13312) = 1/2.D0*
     &   ((2*Cd(13311)*UStop2(1,2) - 
     &       Cd(12701)*MStop2(3)**2*UStop2(1,1)*UUStop(1,2))*
     &     UUStopC(1,1))

        Cd(13313) = 1/(512.D0*Pi**4)*
     &   (2*Cd(12701)*Cd(13312) - 
     &      Cd(12700)*Cd(13309)*MStop2(3)*UStopC(1,2))/MStop2(3)

	Cd(1325) = 5*MStop2(2) - 6*MStop2(3)

        Cd(13330) = MGl*(1/2.D0*(Cd(13007)*UStopC(1,1)) + 
     &     MTy*Cd(12700)*Cr(75)*UCStop(1,3)*UStop2(1,2))

	Cd(1327) = Cd(1264) - MStop2(3)*MStop2(4)

	Cd(1328) = 22*MStop2(3)**2 - 4*Cd(1327)*UStop2(1,2)**2

	end
