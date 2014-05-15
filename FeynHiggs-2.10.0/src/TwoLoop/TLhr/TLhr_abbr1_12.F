* TLhr_abbr1_12.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_12
	implicit none

#include "TLhr.h"


	Cd(8739) = Cd(8732) - MGl2*MStop2(1)

	Cd(8740) = Cd(8734) - 8*MGl2*Cd(8739)

	Cd(864) = 2*Cd(831)*UStop2(1,1)**2 + Cd(863)*UStop2(1,2)

	Cd(8742) = -2*MGl2*Cd(1529) + Cd(8736)

	Cd(8743) = 2*Cd(8738) - MGl2*Cd(8742)

	Cd(8744) = Cd(44)**2*Cd(8737) - MGl2*Cd(8743)

	Cd(8748) = -2*MGl2*Cd(8722) + Cd(8747)

        Cd(84) = Cd(64)*MGlpTmSt4(2) - 
     &   2*Cd(71)*T134(MGl2,MTy2,MStop2(2))

	Cd(865) = -2*Cd(864) + Cd(862)*UStop2(1,1)

	Cd(8750) = -2*MGl2 + Cd(4366)

	Cd(8751) = Cd(4354) - MGl2*Cd(8750)

	Cd(8752) = 2*Cd(44)**3 - MGl2*Cd(8751)

	Cd(8768) = Cd(4329)*UStop2(1,1)**3 - MGl2**3*UStop2(1,2)**2

        Cd(867) = 2*MUEC*Cd(865)*UCStopC(1,3)**2 + 
     &   MUE*Cd(841)*UStop2(1,1)*UStop2(1,2)

	Cd(8784) = Cd(4354) - MGl2*Cd(8726)

	Cd(8785) = 2*Cd(44)**3 - MGl2*Cd(8784)

        Cd(8786) = Cd(4415)*MStop2(3)*UStop(1,1) + 
     &   M_3*MTy*Cd(8785)*UStop(1,2)

	Cd(8798) = -(MGl2*Cd(4561)) + Cd(8795)

	Cd(8799) = Cd(44)**3 - MGl2*Cd(8798)

        Cd(8808) = 4*M_3C*MTy2*Cd(4493)*MStop2(3)*UCStop(1,3)*
     &    UUStop(1,1) + M_3*Cd(8807)*UUStop(1,3)

	Cd(85) = Cd(50)*Cd(83) - Cd(47)*Cd(84)

	Cd(8853) = Cd(44)**3 + MGl2*Cd(8852)

        Cd(8856) = 2*Cd(8855)*MStop2(3)*UUStop(1,2) - 
     &   4*Cd(44)*Cd(45)*MStop2(1)*UUStop(2,2)

        Cd(8858) = 2*Cd(8857)*MStop2(3)*UUStop(1,2) - 
     &   2*Cd(2314)*MStop2(1)*UUStop(2,2)

	Cd(8861) = Cd(8850)*MStop2(1) + MTy2*Cd(8853)*UStop2(1,2)**3

	Cd(8866) = -2*MGl2*Cd(45) + Cd(4460)

	Cd(8875) = Cd(8842) - MGl2*Cd(8874)

	Cd(8876) = 3*Cd(44)**3 - MGl2*Cd(8875)

	Cd(8898) = 2*Cd(8896) + MGl2*Cd(8897)

	Cd(8902) = MGl2*Cd(8901) + 2*Cd(3113)*MStop2(1)

	Cd(8904) = -MGl2 + Cd(4597)

	Cd(8905) = Cd(4253) - MGl2*Cd(8904)

	Cd(8906) = Cd(48)**3 - MGl2*Cd(8905)

	Cd(8925) = Cd(8924) + 2*MGl2*Cd(3976)*MStop2(2)

        Cd(8931) = 2*(Cd(271)*UStop(1,1) + 
     &     M_3*MTy*Cd(254)*MStop2(3)*UStop(1,2))

        Cd(8932) = -(M_3C*Cd(4630)*UCStop(1,3)**2) + 
     &   2*M_3*Cd(4562)*UStop2(1,1)**2

	Cd(8936) = Cd(4578) + MGl2*Cd(4582)

	Cd(8937) = 3*Cd(48)**3 + MGl2*Cd(8936)

        Cd(885) = Cd(882)*UStop2(1,1) + 
     &   2*Cr(85)*MStop2(3)*MStop2(4)*UStop2(1,2)

	Cd(8954) = Cd(8952) + 2*MGl2*Cd(8953)

        Cd(8957) = Cd(4752)*MStop2(3)*UStop(1,1) + 
     &   M_3*MTy*Cd(4894)*UStop(1,2)

	Cd(8968) = Cd(8958) - 8*MGl2*MStop2(2)

	Cd(8969) = 4*Cd(8963) - MGl2*Cd(8968)

	Cd(8971) = 2*MGl2*Cd(1562) + Cd(8966)

	Cd(8972) = 2*Cd(8967) + MGl2*Cd(8971)

	Cd(8973) = 2*Cd(48)**2*Cd(8964) - MGl2*Cd(8972)

	Cd(8977) = 2*MGl2*Cd(8953) + Cd(8976)

	Cd(8992) = Cd(8989) + MGl2*Cd(8991)

	Cd(8993) = Cd(48)**3*Cd(4659) - MGl2*Cd(8992)

	Cd(9003) = 2*Cd(4565)*UStop2(1,1)**3 - MGl2**3*UStop2(1,2)**2

	Cd(892) = Cd(891) - 2*Cd(888)*UStop2(1,1)

	Cd(9018) = Cd(4691) + MGl2*Cd(4892)

	Cd(9019) = 2*Cd(48)**3 + MGl2*Cd(9018)

        Cd(9020) = Cd(4752)*MStop2(3)*UStop(1,1) + 
     &   M_3*MTy*Cd(9019)*UStop(1,2)

	Cd(893) = Cd(892)*UStop2(1,1) - Cd(889)*UStop2(1,2)

	Cd(9033) = Cd(9032)*MStop2(2) + 2*MTy2*MStop2(3)**2

	Cd(9035) = 2*Cd(9031) + 3*MGl2*Cd(3976)*MStop2(2)

	Cd(9037) = -(MGl2*Cd(4216)) + Cd(9034)

	Cd(9038) = Cd(48)**3 - MGl2*Cd(9037)

        Cd(9047) = 4*M_3C*MTy2*Cd(4827)*MStop2(3)*UCStop(1,3)*
     &    UUStop(1,1) + M_3*Cd(9046)*UUStop(1,3)

	Cd(9079) = Cd(48)**3 + MGl2*Cd(9078)

        Cd(9082) = 2*Cd(9081)*MStop2(3)*UUStop(1,2) + 
     &   4*Cd(48)*Cd(49)*MStop2(2)*UUStop(2,2)

        Cd(9084) = 2*(Cd(9083)*MStop2(3)*UUStop(1,2) + 
     &     Cd(3113)*MStop2(2)*UUStop(2,2))

	Cd(9088) = Cd(9077)*MStop2(2) + MTy2*Cd(9079)*UStop2(1,2)**3

	Cd(9093) = -2*MGl2*Cd(49) + Cd(9092)

	Cd(9099) = MTy2*Cd(9098) + Cd(2002)*MStop2(2)*MStop2(3)

	Cd(9100) = Cd(48)*Cd(9096) + 2*MGl2*Cd(9099)

	Cd(9101) = 2*Cd(48)**3*Cd(3989) + MGl2*Cd(9100)

	Cd(9103) = Cd(4853) - MGl2*Cd(9102)

	Cd(9104) = 3*Cd(48)**3 - MGl2*Cd(9103)

        Cd(9105) = 2*(MTy*Cd(9101)*UStop(1,1) + 
     &     M_3*Cd(9104)*MStop2(2)*MStop2(3)*UStop(1,2))

	Cd(9126) = 2*Cd(44)**3 - MGl2*Cd(2130)

        Cd(9127) = Cd(2128)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(9126)*MStop2(3)*UStop(1,2)

	Cd(9133) = MGl2*Cd(2179) + 2*Cd(9132)

	Cd(9134) = Cd(48)**2*Cd(1529) - MGl2*Cd(9133)

	Cd(9135) = Cd(2178)*UStop(1,1) - 2*M_3*MTy*Cd(9134)*UStop(1,2)

	Cd(9137) = 2*MGl2*Cd(2164) + Cd(9136)

	Cd(9138) = Cd(48)**2*Cd(2157) - 2*MGl2*Cd(9137)

        Cd(9139) = Cd(2163)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(9138)*MStop2(3)*UStop(1,2)

        Cd(9140) = Cd(64)*Cd(9135) + 
     &   2*Cd(96)**2*Cd(97)*Cd(1508)*MStop2(3)

        Cd(9141) = Cd(9140)*MStop2(3) - 
     &   Cd(9139)*T134(MGl2,MTy2,MStop2(2))

        Cd(9144) = -(Cd(86)*Cd(167)*Cd(2190)) + 
     &   2*MTy*Cd(9143)*MGlpTmSt4(1)*UStop2(1,2)

	Cd(9147) = Cd(9143)*MGlpTmSt4(1) + Cd(86)*Cd(9146)*UStop(1,1)

	Cd(9148) = Cd(9147)*MStop2(3) + M_3*MTy*Cd(2212)*UStop(1,2)

	Cd(9150) = 2*Cd(9143)*MGlpTmSt4(1) + Cd(86)*Cd(9146)*UStop(1,1)

	Cd(9151) = Cd(9150)*MStop2(3) + M_3*MTy*Cd(2212)*UStop(1,2)

	Cd(9155) = Cd(96)*UStop(1,1) + M_3*MTy*UStop(1,2)

        Cd(9156) = -(Cd(96)*Cd(254)*Cd(2244)) + 
     &   2*MTy*Cd(9155)*MGlpTmSt4(2)*UStop2(1,2)

        Cd(9172) = Cd(9171)*MGlpTmSt4(1) - 
     &   2*M_3C*MTy*Cd(1495)*UCStop(1,3)

        Cd(9174) = Cd(9172)*UUStop(1,1) - 
     &   4*M_3*MTy*Cd(9173)*UUStop(1,3)

	Cd(9176) = Cd(44)**3 - MGl2*Cd(9175)

        Cd(9177) = 2*Cd(86)**2*Cd(167) - 
     &   MGlpTmSt2(1)*MGlpTmSt4(1)*UStop2(1,2)

	Cd(9178) = Cd(9176) - Cd(9177)*UStop2(1,2)**2

	Cd(9179) = Cd(164)**2 + Cd(1495)*UStop2(1,2)**2

        Cd(9180) = Cd(9179)*UCStop(1,3) - 
     &   2*Cd(86)*Cd(167)*UUStop(1,3)*UUStopC(1,2)

        Cd(9183) = 2*Cd(167)*Cd(2647) - 
     &   Cd(86)**2*Cd(9181)*UStop2(1,2)**2

        Cd(9184) = M_3C*Cd(9182)*UUStop(1,1) + 
     &   M_3*Cd(86)*Cd(167)*UUStop(1,2)

	Cd(9185) = Cd(9183)*UStopC(1,1) - 4*MTy*Cd(9184)*UUStopC(1,3)

	Cd(9191) = -2*MGl2 + Cd(9188)

	Cd(9192) = 2*Cd(9189) + MGl2*Cd(9191)

	Cd(9193) = Cd(44)**2*Cd(1701) - MGl2*Cd(9192)

	Cd(9196) = Cd(44)**2*Cd(1701) + MGl2*Cd(9195)

        Cd(9213) = 4*M_3C*MTy2*Cd(9196)*MStop2(3)*UCStop(1,3)*
     &    UUStop(1,1) - 2*M_3*Cd(9212)*UUStop(1,3)

	Cd(9217) = Cd(44)**3*Cd(8310) - MGl2*Cd(9216)

	Cd(9221) = 2*Cd(9219) + MGl2*Cd(9220)

	Cd(9222) = 4*Cd(44)*Cd(9218) - MGl2*Cd(9221)

	Cd(9223) = Cd(44)**3*Cd(2064) - MGl2*Cd(9222)

	Cd(9224) = Cd(9223) - 6*MTy2*Cd(1985)*MGlpTmSt2(1)*UStop2(1,2)

	Cd(9225) = Cd(9217) - Cd(9224)*UStop2(1,2)**2

	Cd(9227) = MTy2*Cd(2170) + 3*MStop2(1)*MStop2(3)

	Cd(9228) = MTy2*Cd(2138) - 9*MStop2(1)*MStop2(3)

	Cd(9231) = MTy2*Cd(9230) + 9*MStop2(1)**2*MStop2(3)

	Cd(9232) = -(MGl2*Cd(2964)) + Cd(9228)

	Cd(9233) = Cd(9231) + MGl2*Cd(9232)

	Cd(9234) = Cd(44)**2*Cd(9227) - MGl2*Cd(9233)

        Cd(9235) = Cd(1985)*Cd(9226)*MStop2(3)**2 + 
     &   2*MTy2*Cd(1505)*UStop2(1,2)**2

        Cd(9236) = Cd(9235)*UCStop(1,3) - 
     &   Cd(9234)*MStop2(3)*UUStop(1,3)*UUStopC(1,2)

	Cd(9238) = Cd(44)*Cd(2068) + MGl2*Cd(9237)

	Cd(9239) = Cd(44)**3*Cd(2978) - MGl2*Cd(9238)

        Cd(9240) = 3*Cd(9239)*MStop2(3)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(1823)*UStop(1,2)

	Cd(9251) = Cd(9246) - MGl2*Cd(9250)

	Cd(916) = Cd(915) + 2*Cd(359)*Cd(360)*MStop2(2)

	Cd(9252) = Cd(44)**2*Cd(2994) - MGl2*Cd(9251)

	Cd(9256) = MTy2*Cd(9255) - 3*MStop2(1)*MStop2(3)

	Cd(9258) = MTy2*Cd(9257) - 9*MStop2(1)**2*MStop2(3)

	Cd(9260) = Cd(9258) + MGl2*Cd(9259)

	Cd(9273) = -3*MGl2*Cd(1873) + 2*Cd(9269)

	Cd(9274) = 3*Cd(9270) + MGl2*Cd(9273)

	Cd(9275) = 3*Cd(44)*Cd(9271) - MGl2*Cd(9274)

	Cd(9276) = Cd(44)**3*Cd(3663) - MGl2*Cd(9275)

	Cd(9282) = Cd(44)**3*Cd(89) + MGl2*Cd(9281)

        Cd(920) = Cd(919)*MStop2(1)*UStopC(1,1) + 
     &   Cd(693)*Cd(694)*MStop2(3)*UUStopC(1,3)

        Cd(9360) = M_3*Cd(118)*Cd(9359) - 
     &   4*M_3C*Cd(9282)*UCStop(1,3)*UUStop(1,1)

	Cd(9365) = Cd(44)**3*Cd(4081) - MGl2*Cd(9364)

	Cd(9367) = 6*Cd(44)*MStop2(3)*UUStop(1,2) + Cd(1786)*UUStop(2,2)

        Cd(9369) = 4*Cd(9368)*MStop2(3)*UUStop(1,2) + 
     &   Cd(44)**2*Cd(1684)*UUStop(2,2)

	Cd(9370) = -2*MGl2*Cd(9367) + Cd(9369)

        Cd(9373) = 2*Cd(9371)*MStop2(3)*UUStop(1,2) + 
     &   Cd(9372)*UUStop(2,2)

	Cd(9374) = Cd(44)*Cd(9370) - MGl2**2*Cd(9373)

	Cd(9375) = Cd(44)**3*Cd(3719) - MGl2*Cd(9374)

        Cd(9376) = Cd(9365) + 
     &   2*MTy2*Cd(9366)*MGlpTmSt2(1)*UStop2(1,2)**3

	Cd(9377) = 2*Cd(9376)*MStop2(3) - Cd(9375)*UUStopC(1,2)

	Cd(9379) = 3*MGlmT2**3 - Cd(9378)*MStop2(1)

	Cd(9380) = 4*MTy2*Cd(1591) - 3*MStop2(3)**2

	Cd(9381) = Cd(9380)*UUStop(1,3) + 4*Cd(45)*MStop2(3)*UUStop(2,3)

	Cd(9383) = MTy2*Cd(9382) - 3*MStop2(1)*MStop2(3)**2

        Cd(9384) = Cd(9383)*UUStop(1,3) + 
     &   Cd(44)**2*MStop2(3)*UUStop(2,3)

	Cd(9386) = MTy2*Cd(9385) - 9*MStop2(1)*MStop2(3)**2

        Cd(9387) = Cd(9386)*UUStop(1,3) + 
     &   2*Cd(2314)*MStop2(3)*UUStop(2,3)

	Cd(9390) = MTy2*Cd(9389) - 9*MStop2(1)**2*MStop2(3)**2

        Cd(9391) = Cd(9390)*UUStop(1,3) + 
     &   4*Cd(44)**2*Cd(45)*MStop2(3)*UUStop(2,3)

        Cd(9396) = Cd(9379)*MStop2(3)**2 - 
     &   4*MTy2*Cd(1499)*UStop2(1,2)**2

	Cd(931) = Cd(930) + 8*Cd(359)*Cd(360)*MStop2(1)

	Cd(9401) = -2*MGl2 + Cd(9398)

	Cd(9402) = Cd(9400) - 2*MGl2*Cd(9401)

	Cd(9405) = MTy2*Cd(1646) - 2*MStop2(1)*MStop2(3)

	Cd(9409) = MTy2*Cd(9408) + 2*MStop2(1)**2*MStop2(3)

	Cd(9413) = Cd(3758) - 3*MGl2*Cd(4216)

	Cd(9414) = Cd(44)**2*Cd(1659) - MGl2*Cd(9413)

	Cd(9418) = MTy2*Cd(9417) + 9*MStop2(1)*MStop2(3)**2

        Cd(9419) = Cd(9418)*UUStop(1,1) - 
     &   2*Cd(2314)*MStop2(3)*UUStop(2,1)

	Cd(9421) = MTy2*Cd(9420) + 3*MStop2(1)*MStop2(3)**2

        Cd(9422) = Cd(9421)*UUStop(1,1) - 
     &   Cd(44)**2*MStop2(3)*UUStop(2,1)

	Cd(9425) = MTy2*Cd(9424) + 9*MStop2(1)**2*MStop2(3)**2

        Cd(9426) = Cd(9425)*UUStop(1,1) - 
     &   4*Cd(44)**2*Cd(45)*MStop2(3)*UUStop(2,1)

	Cd(9428) = 4*MTy2*Cd(1591) + 3*MStop2(3)**2

        Cd(9429) = MGl2**2*Cd(9428) + 
     &   16*MTy2**2*MStop2(3)**2*UStop2(1,2)

        Cd(9432) = Cd(9383)*UUStop(1,2) + 
     &   Cd(44)**2*MStop2(3)*UUStop(2,2)

        Cd(9433) = Cd(9386)*UUStop(1,2) + 
     &   2*Cd(2314)*MStop2(3)*UUStop(2,2)

        Cd(9434) = Cd(9390)*UUStop(1,2) + 
     &   4*Cd(44)**2*Cd(45)*MStop2(3)*UUStop(2,2)

	Cd(9446) = Cd(96)*UStop(1,1) - M_3*MTy*UStop(1,2)

        Cd(9448) = Cd(9447)*MGlpTmSt4(2) + 
     &   2*M_3C*MTy*Cd(1510)*UCStop(1,3)

	Cd(9452) = Cd(48)**3 - MGl2*Cd(9451)

        Cd(9453) = 2*Cd(96)**2*Cd(254) - 
     &   MGlpTmSt2(2)*MGlpTmSt4(2)*UStop2(1,2)

	Cd(9454) = Cd(9452) - Cd(9453)*UStop2(1,2)**2

	Cd(9455) = Cd(251)**2 + Cd(1510)*UStop2(1,2)**2

        Cd(9456) = Cd(9455)*UCStop(1,3) - 
     &   2*Cd(96)*Cd(254)*UUStop(1,3)*UUStopC(1,2)

	Cd(9457) = 3*Cd(96)*UStop(1,1) + 2*M_3*MTy*UStop(1,2)

        Cd(9459) = 2*Cd(254)*Cd(2725) - 
     &   Cd(96)**2*Cd(9457)*UStop2(1,2)**2

        Cd(9460) = M_3C*Cd(9458)*UUStop(1,1) + 
     &   M_3*Cd(96)*Cd(254)*UUStop(1,2)

	Cd(937) = Cd(930)*Cr(34)*MStop2(2) - 2*Cd(936)*MStop2(3)

	Cd(9461) = Cd(9459)*UStopC(1,1) + 4*MTy*Cd(9460)*UUStopC(1,3)

	Cd(9465) = -2*MGl2 + Cd(9463)

	Cd(9466) = 2*Cd(9464) + MGl2*Cd(9465)

	Cd(9467) = Cd(48)**2*Cd(2063) - MGl2*Cd(9466)

        Cd(938) = 3/(64.D0*Pi)*
     &   (Alfa2L*(Cd(937) + 
     &        8*Cd(359)*Cd(360)*Cr(13)*MStop2(1)*MStop2(2)))/
     &    (SB2*SW2*MStop2(2)*MStop2(3)**2)

	Cd(9471) = Cd(48)**2*Cd(2063) - MGl2*Cd(9470)

        Cd(9483) = 4*M_3C*MTy2*Cd(9471)*MStop2(3)*UCStop(1,3)*
     &    UUStop(1,1) + 2*M_3*Cd(9482)*UUStop(1,3)

	Cd(9487) = Cd(48)**3*Cd(9485) - MGl2*Cd(9486)

	Cd(940) = Cd(936)*Cr(34) - 2*Cd(939)*MStop2(3)

	Cd(9491) = 2*Cd(9488) + MGl2*Cd(9490)

	Cd(9492) = 4*Cd(48)*Cd(9489) - MGl2*Cd(9491)

	Cd(9493) = Cd(48)**3*Cd(3409) - MGl2*Cd(9492)

	Cd(9494) = Cd(9493) - 6*MTy2*Cd(2008)*MGlpTmSt2(2)*UStop2(1,2)

	Cd(9495) = Cd(9487) - Cd(9494)*UStop2(1,2)**2

	Cd(9498) = MTy2*Cd(2173) + 9*MStop2(2)*MStop2(3)

	Cd(941) = 1/2.D0*(A0delStop(2)*Cd(930)) - Cd(940)*MStop2(3)

	Cd(9502) = MTy2*Cd(9501) - 9*MStop2(2)**2*MStop2(3)

	Cd(9503) = -(MGl2*Cd(2957)) + Cd(9498)

	Cd(9504) = Cd(9502) + MGl2*Cd(9503)

	Cd(9505) = Cd(48)**2*Cd(9497) - MGl2*Cd(9504)

        Cd(9506) = Cd(2008)*Cd(9496)*MStop2(3)**2 + 
     &   2*MTy2*Cd(1514)*UStop2(1,2)**2

        Cd(9507) = Cd(9506)*UCStop(1,3) + 
     &   Cd(9505)*MStop2(3)*UUStop(1,3)*UUStopC(1,2)

	Cd(9508) = MGl2*Cd(1705) + Cd(48)**3*Cd(3335)

        Cd(9509) = 3*Cd(9508)*MStop2(3)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(1779)*UStop(1,2)

        Cd(942) = 3/(64.D0*Pi)*
     &   (Alfa2L*(Cd(941) + 
     &        4*A0delStop(1)*Cd(359)*Cd(360)*MStop2(2)))/
     &    (SB2*SW2*MStop2(2))

	Cd(9513) = MTy2*Cd(9512) + MStop2(2)**2*MStop2(3)

	Cd(9518) = 2*Cd(9511) + MGl2*Cd(9517)

	Cd(9521) = Cd(9516) - MGl2*Cd(9520)

	Cd(9522) = Cd(48)**2*Cd(3356) - MGl2*Cd(9521)

	Cd(9525) = MTy2*Cd(9524) - 9*MStop2(2)*MStop2(3)

	Cd(9527) = MTy2*Cd(9526) + 3*MStop2(2)*MStop2(3)

	Cd(9529) = MTy2*Cd(9528) + 9*MStop2(2)**2*MStop2(3)

	Cd(944) = -2*A0delT + MTy2*Cd(943)

	Cd(92) = -2*MGl2 + Cd(88)

	Cd(9545) = -3*MGl2*Cd(1855) + 2*Cd(9540)

	Cd(9546) = 3*Cd(9541) + MGl2*Cd(9545)

	Cd(9547) = 3*Cd(48)*Cd(9542) - MGl2*Cd(9546)

	Cd(9548) = Cd(48)**3*Cd(3915) - MGl2*Cd(9547)

	Cd(9549) = -(MGl2*Cd(1552)) + Cd(9544)

	Cd(9550) = Cd(48)**3*Cd(89) - MGl2*Cd(9549)

	Cd(946) = 9*Cd(944) - 10*MTy2*Cd(945)

	Cd(9551) = Cd(9548)*UStop(1,1) - 2*M_3*MTy*Cd(9550)*UStop(1,2)

	Cd(9553) = -(MGl2*Cd(1741)) + Cd(48)**2*Cd(9552)

	Cd(9554) = Cd(48)**4*Cd(187) - MGl2*Cd(9553)

        Cd(9555) = 3*Cd(9554)*UStop(1,1) - 
     &   8*M_3*MGl2*MTy*MTy2*MStop2(3)**2*UStop(1,2)

	Cd(9559) = Cd(48)**3*Cd(89) + MGl2*Cd(9558)

        Cd(9589) = -(Cd(138)**2*Cd(9588)) + 
     &   4*Cd(9555)*UStop(1,1)*UStop2(1,2)

        Cd(9624) = M_3*Cd(138)*Cd(9623) - 
     &   4*M_3C*Cd(9559)*UCStop(1,3)*UUStop(1,1)

	Cd(9628) = Cd(48)**3*Cd(4178) - MGl2*Cd(9627)

	Cd(9631) = 6*Cd(48)*MStop2(3)*UUStop(1,2) - Cd(9630)*UUStop(2,2)

        Cd(9633) = 4*Cd(9632)*MStop2(3)*UUStop(1,2) - 
     &   Cd(48)**2*Cd(1845)*UUStop(2,2)

	Cd(9634) = -2*MGl2*Cd(9631) + Cd(9633)

        Cd(9637) = 2*Cd(9635)*MStop2(3)*UUStop(1,2) - 
     &   Cd(9636)*UUStop(2,2)

	Cd(9638) = Cd(48)*Cd(9634) - MGl2**2*Cd(9637)

	Cd(9639) = Cd(48)**3*Cd(3970) - MGl2*Cd(9638)

        Cd(9640) = Cd(9628) + 
     &   2*MTy2*Cd(9629)*MGlpTmSt2(2)*UStop2(1,2)**3

	Cd(93) = Cd(44)*Cd(89) - MGl2*Cd(92)

	Cd(955) = Cd(950) + 2*Cd(954)*UStop2(1,1)

	Cd(9641) = 2*Cd(9640)*MStop2(3) - Cd(9639)*UUStopC(1,2)

	Cd(9643) = 3*MGlmT2**3 - Cd(9642)*MStop2(2)

	Cd(9644) = 4*MTy2*Cd(1540) - 3*MStop2(3)**2

	Cd(9645) = Cd(9644)*UUStop(1,3) + 4*Cd(49)*MStop2(3)*UUStop(2,3)

	Cd(9647) = MTy2*Cd(9646) - 3*MStop2(2)*MStop2(3)**2

        Cd(9648) = Cd(9647)*UUStop(1,3) + 
     &   Cd(48)**2*MStop2(3)*UUStop(2,3)

	Cd(9651) = MTy2*Cd(9650) - 9*MStop2(2)*MStop2(3)**2

        Cd(9652) = Cd(9651)*UUStop(1,3) + 
     &   2*Cd(3113)*MStop2(3)*UUStop(2,3)

	Cd(9657) = MTy2*Cd(9656) - 9*MStop2(2)**2*MStop2(3)**2

        Cd(9658) = Cd(9657)*UUStop(1,3) + 
     &   4*Cd(48)**2*Cd(49)*MStop2(3)*UUStop(2,3)

        Cd(9662) = Cd(9643)*MStop2(3)**2 - 
     &   4*MTy2*Cd(1519)*UStop2(1,2)**2

	Cd(9667) = -2*MGl2 + Cd(9664)

	Cd(9668) = Cd(9666) - 2*MGl2*Cd(9667)

	Cd(9671) = MTy2*Cd(1669) + 2*MStop2(2)*MStop2(3)

	Cd(9674) = MTy2*Cd(9673) - 2*MStop2(2)**2*MStop2(3)

	Cd(9681) = -3*MGl2*Cd(4561) + Cd(9677)

	Cd(9682) = Cd(48)**2*Cd(1636) - MGl2*Cd(9681)

	Cd(9686) = MTy2*Cd(9685) + 9*MStop2(2)*MStop2(3)**2

        Cd(9687) = Cd(9686)*UUStop(1,1) - 
     &   2*Cd(3113)*MStop2(3)*UUStop(2,1)

	Cd(9689) = MTy2*Cd(9688) + 3*MStop2(2)*MStop2(3)**2

        Cd(9690) = Cd(9689)*UUStop(1,1) - 
     &   Cd(48)**2*MStop2(3)*UUStop(2,1)

	Cd(9695) = MTy2*Cd(9694) + 9*MStop2(2)**2*MStop2(3)**2

        Cd(9696) = Cd(9695)*UUStop(1,1) - 
     &   4*Cd(48)**2*Cd(49)*MStop2(3)*UUStop(2,1)

	Cd(9697) = 4*MTy2*Cd(1540) + 3*MStop2(3)**2

        Cd(9698) = MGl2**2*Cd(9697) + 
     &   16*MTy2**2*MStop2(3)**2*UStop2(1,2)

        Cd(9701) = Cd(9647)*UUStop(1,2) + 
     &   Cd(48)**2*MStop2(3)*UUStop(2,2)

        Cd(9702) = Cd(9651)*UUStop(1,2) + 
     &   2*Cd(3113)*MStop2(3)*UUStop(2,2)

        Cd(9703) = Cd(9657)*UUStop(1,2) + 
     &   4*Cd(48)**2*Cd(49)*MStop2(3)*UUStop(2,2)

	Cd(9721) = Cd(9718) - Cd(9720)

	Cd(9727) = Cd(9724) - Cd(9726)

        Cd(9733) = (Cd(9729) + 2*MTy*Cd(9732)*MStop2(3))/
     &   (MTy*MZ2*SW2*MStop2(3))

        Cd(9736) = -((Cd(9734) - 2*MTy*Cd(9735)*MStop2(3))/
     &     (MTy*MZ2*SW2*MStop2(3)))

	Cd(9737) = Cd(9733) + Cd(9736)

        Cd(9738) = 3/(8.D0*Pi)*
     &   (Alfa2L*(-48*MTy*SA2 + MZ2*SW2*Cd(9737)))/(SB2*SW2)

        Cd(9740) = 1/2.D0*(A0delStop(2)*Cd(9718)) - 
     &   Cd(9739)*Cr(34)*MStop2(2)*MStop2(3)*UCStopC(1,3)

        Cd(94) = 2*(Cd(91)*MStop2(3)*UStop(1,1) + 
     &     M_3*MTy*Cd(93)*UStop(1,2))

	Cd(965) = Cd(964) + 2*Cd(954)*UStop2(1,1)

        Cd(9742) = 1/2.D0*(A0delStop(1)*Cd(9720)) - 
     &   Cd(9741)*Cr(13)*MStop2(1)*MStop2(3)*UCStopC(1,3)

        Cd(9743) = 3/(16.D0*Pi)*
     &   (Alfa2L*SA*(-(Cd(9740)*MStop2(1)) + Cd(9742)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(2)*MStop2(3))

	Cd(9744) = Cd(9718)*Cr(34) - 2*Cd(9739)*MStop2(3)*UCStopC(1,3)

	Cd(9745) = Cd(9720)*Cr(13) - 2*Cd(9741)*MStop2(3)*UCStopC(1,3)

	Cd(9746) = Cd(9744) - Cd(9745)

        Cd(9748) = 1/2.D0*(A0delStop(2)*Cd(9724)) - 
     &   Cd(9747)*Cr(34)*MStop2(2)*MStop2(3)*UCStop(1,3)

        Cd(9750) = 1/2.D0*(A0delStop(1)*Cd(9726)) - 
     &   Cd(9749)*Cr(13)*MStop2(1)*MStop2(3)*UCStop(1,3)

        Cd(9751) = 3/(16.D0*Pi)*
     &   (Alfa2L*SA*(-(Cd(9748)*MStop2(1)) + Cd(9750)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(2)*MStop2(3))

	Cd(9752) = Cd(9724)*Cr(34) - 2*Cd(9747)*MStop2(3)*UCStop(1,3)

	Cd(9753) = Cd(9726)*Cr(13) - 2*Cd(9749)*MStop2(3)*UCStop(1,3)

	Cd(9754) = Cd(9752) - Cd(9753)

	Cd(9757) = (MTy*Cd(9732) + Cd(9756))/(MTy*MZ2*SW2)

        Cd(9758) = 1/2.D0*(A0delStop(2)*Cd(9733)) - 
     &   Cd(9757)*Cr(34)*MStop2(2)

	Cd(9759) = (MTy*Cd(9735) + Cd(9756))/(MTy*MZ2*SW2)

        Cd(9760) = 1/2.D0*(A0delStop(1)*Cd(9736)) - 
     &   Cd(9759)*Cr(13)*MStop2(1)

	Cd(9761) = Cd(9758)*MStop2(1) + Cd(9760)*MStop2(2)

        Cd(9762) = 3/(16.D0*Pi)*
     &   (Alfa2L*(MTy*MZ2*SW2*Cd(9761) + 
     &        4*SA2*Cd(456)*MStop2(1)*MStop2(2)))/
     &    (MTy*SB2*SW2*MStop2(1)*MStop2(2))

	Cd(9763) = 2*Cd(9757) - Cd(9733)*Cr(34)

	Cd(9764) = 2*Cd(9759) - Cd(9736)*Cr(13)

	Cd(9765) = Cd(9763) + Cd(9764)

        Cd(9766) = 3/(16.D0*Pi)*
     &   (Alfa2L*(8*SA2*Cd(463) + MTy*MZ2*SW2*Cd(9765)))/
     &    (MTy*SB2*SW2)

	Cd(9771) = SA2*Cd(484) + S2A*Cd(9770)

        Cd(9772) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*Cr(9)**2*(CA2*Cd(488) + Cd(9771)*MStop2(3)))/
     &    (SB2*SW2*MStop2(3))

	Cd(9774) = CA2*Cd(496) - S2A*Cd(499)

	Cd(9775) = SA2*Cd(502) + Cd(9774)*UStop2(1,2)

        Cd(9777) = -(Cd(9776)*MStop2(1)) + 
     &   2*MStop2(3)*UStop2(1,1)**2*UStop2(1,2)**2

	Cd(9782) = Cd(9781) + 2*MStop2(1)*MStop2(3)*UStop2(1,1)

	Cd(9797) = -(Cd(9796)*UCStop(1,3)) + 2*MTy*Cd(9730)*UStop2(1,2)

	Cd(9804) = Cd(9801) + Cd(9803)

	Cd(9805) = 2*(Cd(402)*Cd(9739)*MStop2(1) + Cd(9804)*MStop2(2))

	Cd(9811) = Cd(9796)*UCStop(1,3) + 2*MTy*Cd(9730)*UStop2(1,1)

	Cd(9817) = Cd(9815) - Cd(9816)

	Cd(9818) = 2*(Cd(9817)*MStop2(1) + Cd(403)*Cd(9741)*MStop2(2))

        Cd(9821) = 1/2.D0*(A0delStop(2)*Cd(9805)) + 
     &   Cd(402)*Cd(9739)*Cr(34)*MStop2(2)*MStop2(3)

	Cd(9824) = Cd(9805)*Cr(34) + 2*Cd(402)*Cd(9739)*MStop2(3)

        Cd(9827) = 1/2.D0*(A0delStop(1)*Cd(9818)) - 
     &   Cd(403)*Cd(9741)*Cr(13)*MStop2(1)*MStop2(3)

	Cd(9830) = Cd(9818)*Cr(13) - 2*Cd(403)*Cd(9741)*MStop2(3)

	Cd(9835) = 8*Cr(75)*MStop2(1) - Cd(9834)*MStop2(2)

        Cd(9836) = 2*MTy*Cd(764)*UCStop(1,3) + 
     &   At*Cd(9835)*UStop2(1,1)*UStop2(1,2)

        Cd(95) = Cd(52)*Cd(86)*Cd(87)*MStop2(3) - 
     &   Cd(94)*T134(MGl2,MTy2,MStop2(1))

	Cd(9843) = Cd(9842)*MStop2(3) - Cd(783)*UStop2(1,1)*UStop2(1,2)

        Cd(9844) = -(Cd(791)*MStop2(3)) - 
     &   Cd(789)*UStop2(1,1)*UStop2(1,2)

	Cd(9845) = MUE*Cd(9843)*UCStop(1,3) + MUEC*Cd(9844)*UCStopC(1,3)

        Cd(9846) = -2*MUEC*Cd(794)*UCStopC(1,3)**2 + 
     &   MUE*Cd(9840)*UStop2(1,1)*UStop2(1,2)

	Cd(9854) = Cd(9853) - 2*Cd(835)*UStop2(1,1)

        Cd(9855) = Cd(9854)*UStop2(1,1) - 
     &   2*Cd(836)*MStop2(4)*UStop2(1,2)

	Cd(9859) = Cd(9857) - 2*Cd(842)*UStop2(1,1)

        Cd(9877) = (2*(Cd(9876) - 
     &       4*At*MTy*MStop2(3)*UCStopC(1,3))*UStop2(1,1)*
     &     UStop2(1,2))/(MZ2*SW2*MStop2(3))

	Cd(9881) = -(MZ2*SA2*SW2*Cd(9877)) + S2A*Cd(9879)

        Cd(9882) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(9881)*MStop2(3) - 
     &        2*CA2*Cd(9880)*UStop2(1,1)*UStop2(1,2)))/
     &    (SB2*SW2*MStop2(3))

	Cd(9893) = -2*CA2*Cd(474) - S2A*Cd(9891)

        Cd(9894) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*(-4*Cd(9893)*Cr(75)**2*MStop2(1) + 
     &        SA2*Cd(9892)*MStop2(3)**2))/(SB2*SW2*MStop2(3)**2)

        Cd(9908) = (-2*Cd(9907) + 
     &     4*MTy*Cd(13)*Cr(75)**2*MStop2(3))/(MZ2*SW2*MStop2(3))

	Cd(9909) = SA2*Cd(9906) + MZ2*S2A*SW2*Cd(9908)

        Cd(9910) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*(2*CA2*Cd(9902) + Cd(9909)*MStop2(3)))/
     &    (SB2*SW2*MStop2(3))

        Cd(9921) = Cd(9920)*MStop2(1)*UStopC(1,1) - 
     &   Cd(9793)*Cd(9795)*MStop2(3)*UUStopC(1,3)

	Cd(9933) = Cd(9932) + 8*Cd(402)*Cd(403)*MStop2(1)

        Cd(9939) = 2*Cd(9936) + Cd(9793)*Cd(9938)*MStop2(2) + 
     &   2*Cd(9937)*MStop2(3)

        Cd(9951) = Cd(52)*Cd(1985) - 2*Cd(1493)*MGlpTmSt4(1) + 
     &   4*MGl2*MTy2*T134(MGl2,MTy2,MStop2(1))

        Cd(9952) = Cd(64)*Cd(2008) - 2*Cd(1508)*MGlpTmSt4(2) + 
     &   4*MGl2*MTy2*T134(MGl2,MTy2,MStop2(2))

        Cd(9953) = Cd(52)*Cd(1985) + 
     &   2*Cd(164)*Cd(1493)*MGlpTmSt2(1) - 
     &   2*MGl2*Cd(152)*T134(MGl2,MTy2,MStop2(1))

        Cd(9954) = Cd(64)*Cd(2008) + 
     &   2*Cd(251)*Cd(1508)*MGlpTmSt2(2) - 
     &   2*MGl2*Cd(179)*T134(MGl2,MTy2,MStop2(2))

        Cd(9968) = At*Cd(9966) + 2*M_3*Cr(28) - 
     &   2*Cd(9967)*UCStop(1,3)**2

        Cd(989) = UStop2(1,1)*
     &    (Cd(986) - 2*Cd(988)*UStop2(1,1)) + 2*Cd(987)*UStop2(1,2)

        Cd(9984) = -(Cd(50)*Cd(9947)*Cd(9966)*MGlpTmSt2(1)) + 
     &   Cd(47)*Cd(9949)*MGlpTmSt2(2)*
     &    (Cd(9948)*Cd(9983) + Cr(75)**2*UStop2(1,1))

        Cd(9999) = Cd(9948)*Cd(9998) + 
     &   2*Cd(9997)*UCStopC(1,3)**2 + AtC*Cd(9958)*UStop2(1,1)

        Cd(10002) = AtC*Cd(9966) + 2*M_3C*Cr(28) - 
     &   2*Cd(10001)*UCStopC(1,3)**2

        Cd(10101) = 2*Cd(164)**2*Cd(1493) - Cd(52)*Cd(1985) - 
     &   4*MGl2*MStop2(1)*T134(MGl2,MTy2,MStop2(1))

        Cd(10102) = 2*Cd(251)**2*Cd(1508) - Cd(64)*Cd(2008) - 
     &   4*MGl2*MStop2(2)*T134(MGl2,MTy2,MStop2(2))

        Cd(1004) = UStop2(1,1)*
     &    (Cd(986) - 2*Cd(988)*UStop2(1,1)) + 
     &   4*Cd(1003)*UStop2(1,2)

        Cd(10150) = Cd(164)*Cd(10147) + Cd(10149)*MGlpTmSt2(1) + 
     &   AtC*Cd(10026)*UCStop(1,3)

        Cd(10152) = Cd(251)*Cd(10147) + Cd(10149)*MGlpTmSt2(2) + 
     &   AtC*Cd(10151)*UCStop(1,3)

        Cd(10197) = 8*MTy2*MStop2(1)**2*
     &    (Cd(164)*Cr(5) + 2*M_3C*MGlpTmSt2(1)*UCStop(1,3)) + 
     &   Cd(10196)*UCStopC(1,3)

        Cd(10199) = Cd(152)*Cd(2597)*Cd(10198) + 
     &   (Cd(47)**2 - 4*MGlpTmSt4(1)*MStop2(1)**2)*UCStop(1,3)

        Cd(10206) = Cd(9948)*Cd(10205) + 2*Cr(28)*MGlpTmSt2(1) + 
     &   Cd(164)*UStop2(1,1)

        Cd(10220) = At*Cd(10218)*Cd(10219) + 
     &   AtC*(MUE*Cd(10218)*UCStop(1,3)**2 + 
     &      MUEC*Cd(10217)*UStop2(1,1))

        Cd(1014) = -8*Cd(1013)*UStop2(1,2)**2 + 
     &   MStop2(1)*(5*MStop2(3) - 14*MStop2(1)*UStop2(1,2))

        Cd(10232) = -(Cd(10218)*
     &      (MUE**2*UCStop(1,3)**2 + MUEC**2*UCStopC(1,3)**2)) + 
     &   2*MUE2*Cd(10231)*UStop2(1,1)

        Cd(10262) = Cd(50)**2 + Cd(179)*Cd(2681)*MGlpTmSt2(2) - 
     &   4*MGlpTmSt4(2)*MStop2(2)**2

        Cd(10281) = Cd(10193)*Cd(10279) + 
     &   4*MStop2(1)*(Cd(10278) - Cd(164)*MStop2(1))

        Cd(10296) = AtC**2*Cd(10279)*UCStop(1,3)**2 + 
     &   At*(At*Cd(10279)*UCStopC(1,3)**2 + 
     &      4*AtC*Cd(10280)*UStop2(1,1)*UStop2(1,2))

        Cd(10310) = Cd(47)*Cd(10309) + 
     &   AtC*Cd(10307)*UCStop(1,3) + At*Cd(10278)*UCStopC(1,3)

        Cd(10315) = -4*Cd(10280)*Cr(5) + 
     &   AtC*Cd(10279)*UCStop(1,3) + At*Cd(10279)*UCStopC(1,3)

        Cd(10319) = AtC*(Cd(10317) + 
     &      AtC*Cd(10086)*Cd(10280)*UCStop(1,3)) + 
     &   At*Cd(10318)*UCStopC(1,3)

        Cd(1024) = Cd(1022)*MStop2(3) + 4*Cd(1023)*UStop2(1,1) - 
     &   4*Cd(832)*UStop2(1,2)**2

        Cd(10327) = 2*MGlmT2*Cd(47)*Cd(10149) + 
     &   AtC*Cd(10326)*UCStop(1,3) + At*Cd(10279)*UCStopC(1,3)

        Cd(10330) = Cd(9948)*Cd(10328) + At*MGlmT2*Cd(10329) - 
     &   2*Cd(47)*UStop2(1,1)*UStop2(1,2)

        Cd(10334) = AtC*(2*Cd(10333) + 
     &      AtC*MGlmT2*Cd(10019)*UCStop(1,3)) + 
     &   At*Cd(10330)*UCStopC(1,3)

        Cd(10343) = -(Cd(47)*Cd(10342)) + 4*Cd(10278)*Cr(28) + 
     &   2*Cd(9948)*Cd(10339)*UStop2(1,1)

        Cd(10348) = AtC*Cd(10346) - 
     &   At*(Cd(10347) + At*Cd(10077)*UCStopC(1,3)**2)

        Cd(10350) = -(AtC*Cd(10279)*UCStop(1,3)**2) + 
     &   16*Cd(10280)*(Cd(10349) + M_3*Cd(9948)*UStop2(1,1)**2)

        Cd(10364) = MGlmT2*Cd(47)*Cd(10358) + 
     &   Cd(9948)*Cd(10363) + Cd(10354)*MStop2(3)

        Cd(10371) = 4*MGlmT2*Cd(47) + 2*Cd(10280)*UStop2(1,1) + 
     &   Cd(10279)*UStop2(1,2)

        Cd(10373) = 4*MGlmT2*Cd(47)*Cd(10372) + 
     &   AtC*Cd(10279)*UCStop(1,3)**2 + 
     &   At*Cd(10371)*UStop2(1,1)*UStop2(1,2)

        Cd(10379) = 4*MGlmT2*Cd(47)*Cd(10378) + 
     &   At*Cd(10279)*UCStopC(1,3)**2 + 
     &   AtC*Cd(10371)*UStop2(1,1)*UStop2(1,2)

        Cd(10387) = -(Cd(9948)*Cd(10386)) + 
     &   Cd(47)*Cd(10385)*MStop2(1) + Cd(10280)*Cd(10382)*MStop2(3)

        Cd(10391) = -(Cd(9948)*Cd(10390)) + 
     &   Cd(47)*Cd(10389)*MStop2(1) + Cd(10280)*Cd(10388)*MStop2(3)

        Cd(10394) = Cd(10280)*Cd(10393) + 
     &   8*Cd(10052)*Cd(10278)*MStop2(1) - 
     &   AtC*MUEC*Cd(10279)*UStop2(1,1)**2*UStop2(1,2)

        Cd(10404) = AtC*Cd(10400) + 4*M_3C*Cd(10402) + 
     &   4*M_3*Cd(10403)

	Cd(10411) = MGl2*Cd(13) - AtC*Cd(10409) - At*Cd(10410)

        Cd(10415) = MGlpT2*UCStop(1,3)*UStop2(1,1) + 
     &   Cd(10414)*UStop2(1,2)

	Cd(10426) = -Cd(11770) + AtC*Cd(10088)*UCStop(1,3)

        Cd(10441) = MUE2*Cd(10438) + Cd(427)*Cd(10439) - 
     &   MGlmT2*Cd(9948)*Cd(10440)

        Cd(10442) = 2*Cd(47)**2*Cd(10436) - 
     &   2*MStop2(1)*(2*Cd(47)*Cd(10441) + 
     &      Cd(10280)*Cd(10435)*MStop2(3))

        Cd(10447) = 4*MUE2*Cd(10444) + 
     &   MUE**2*Cd(10446)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(10446)*UCStopC(1,3)**2

        Cd(10452) = Cd(9948)*Cd(10451) - 
     &   4*Cd(47)*Cd(152)*Cd(10449)*MStop2(1) + 
     &   Cd(10279)*Cd(10448)*MStop2(3)

        Cd(10457) = -2*Cd(47)*
     &    (-(M_3*Cd(9948)*Cd(10453)) + Cd(10456))*MStop2(1) - 
     &   Cd(10086)*Cd(10280)*MStop2(3)

        Cd(10459) = 2*Cd(47)*(-(Cd(9948)*Cd(10453)) + Cd(10458))*
     &    MStop2(1) + Cd(10280)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(10461) = -2*Cd(47)*
     &    (-(M_3C*Cd(9948)*Cd(10453)) + Cd(10460))*MStop2(1) - 
     &   Cd(10039)*Cd(10280)*MStop2(3)

        Cd(10462) = -2*MUE2*Cd(10459)*Cr(5) + 
     &   MUE**2*Cd(10457)*UCStop(1,3) + 
     &   MUEC**2*Cd(10461)*UCStopC(1,3)

	Cd(10463) = MGl2*Cd(10447) + MTy2*Cd(10452) + 4*MTy*Cd(10462)

        Cd(1039) = 2*Cd(1037)*UStop2(1,1)**2 + 
     &   Cd(1038)*UStop2(1,2) - 4*MStop2(1)**2*UStop2(1,2)**4

        Cd(10482) = AtC*Cd(10479) - 
     &   At*(Cd(10481) + At*Cd(10077)*UCStopC(1,3)**2)

        Cd(10495) = M_3*Cd(10492) + 
     &   At*(-Cd(10494) + 3*M_3*At*UCStopC(1,3)**2)

        Cd(10515) = Cd(164)*Cd(10507)*MStop2(3) + 
     &   Cd(47)*(Cd(10514) + 
     &      2*Cd(47)*Cd(9955)*UCStop(1,3)*UStop2(1,1))

        Cd(10523) = -(MGlmT2*Cd(10522)) + 
     &   Cd(10519)*UStop2(1,1)**3 + 
     &   Cd(10521)*UStop2(1,1)*UStop2(1,2)**2

	Cd(10524) = Cd(47) - Cd(152)*MStop2(2) + Cd(10519)*UStop2(1,2)

        Cd(10543) = AtC*Cd(10054) + 4*M_3C*Cd(10402) + 
     &   4*M_3*Cd(10403)

        Cd(10551) = Cd(152)*Cd(164)*MStop2(1)*MStop2(3) + 
     &   Cd(47)*(Cd(10550) + 2*Cd(47)*Cd(9955)*UStop2(1,2))

        Cd(10554) = 2*At*Cd(47)*Cd(10553)*Cr(5) + 
     &   (Cd(10551) + 2*AtC*Cd(47)*Cd(10552))*UCStop(1,3)

	Cd(10558) = MGl2*Cd(13) - MUE*Cd(10555) - MUEC*Cd(10557)

        Cd(10567) = Cd(152)*Cd(164)*MStop2(1)*MStop2(3) + 
     &   Cd(47)*(Cd(10566) + 2*Cd(47)*Cd(9955)*UStop2(1,2))

        Cd(10580) = MUE2*Cd(10443)*Cr(28) - 
     &   Cd(10445)*(MUE**2*UCStop(1,3)**2 + 
     &      MUEC**2*UCStopC(1,3)**2)

        Cd(10583) = Cd(47)*Cd(152)*Cd(10449) + 
     &   Cd(9948)*Cd(10582) - Cd(164)*Cd(10581)*MStop2(3)

        Cd(10597) = Cd(47)*
     &    (Cd(9948)*Cd(10453)*Cd(10584) - Cd(10596)) - 
     &   MTy2*Cd(164)*Cd(10587)*MStop2(3)

        Cd(10618) = AtC*Cd(10616) + 
     &   At*(Cd(10617)*UCStopC(1,3)**2 + 
     &      2*M_3C*Cd(9948)*UStop2(1,1))

        Cd(10620) = Cd(10193)*Cd(10603) + 
     &   4*MStop2(2)*(Cd(10602) - Cd(251)*MStop2(2))

        Cd(10621) = AtC**2*Cd(10603)*UCStop(1,3)**2 + 
     &   At*(At*Cd(10603)*UCStopC(1,3)**2 + 
     &      4*AtC*Cd(10604)*UStop2(1,1)*UStop2(1,2))

        Cd(10631) = -4*Cd(10604)*Cr(5) + 
     &   AtC*Cd(10603)*UCStop(1,3) + At*Cd(10603)*UCStopC(1,3)

        Cd(10673) = Cd(10603)*Cd(10668) + 
     &   2*Cd(10604)*(Cd(10672) + 
     &      MGl2*MUEC*Cd(9948)*UStop2(1,1)**2)

        Cd(10680) = Cd(9948)*Cd(10679) + Cd(10677)*UCStop(1,3) - 
     &   2*M_3*MTy*UStop2(1,1)**2

        Cd(10691) = Cd(10603)*Cd(10687) - 
     &   2*Cd(10604)*Cd(10690) + 
     &   2*MGl2*MUE*Cd(9948)*Cd(10604)*UStop2(1,1)**2

        Cd(10694) = Cd(9948)*Cd(10693) + 
     &   Cd(10692)*UCStopC(1,3) - 2*M_3C*MTy*UStop2(1,1)**2

        Cd(10708) = M_3C*Cd(10704) + M_3*Cd(10707) - 
     &   2*MTy*Cd(13)*Cd(10706)*Cr(75)

        Cd(10711) = MTy*Cd(13)*Cd(48)*Cd(10603) + 
     &   2*Cd(10604)*Cd(10709) + 2*Cd(10602)*Cd(10710)*MStop2(2)

        Cd(10719) = M_3C*Cd(10714) + MTy*Cd(13)*Cd(10716) + 
     &   M_3*Cd(10718)

        Cd(10725) = M_3C*Cd(10720) - MTy*Cd(13)*Cd(10723) + 
     &   M_3*Cd(10724)

        Cd(10736) = Cd(50) + 
     &   Cd(9948)*(MStop2(1)*MStop2(2) - MGlmT2*MStop2(3))

        Cd(10744) = Cd(50)**2*Cd(1272) + 
     &   MStop2(2)*(-2*Cd(50)*Cd(10743) + Cd(10604)*MStop2(3))

        Cd(10747) = Cd(50)**2*Cd(843) + 
     &   MStop2(2)*(-2*Cd(50)*Cd(10746) + Cd(10604)*MStop2(3))

        Cd(10751) = MUE2*Cd(10658) - MUEC**2*UCStopC(1,3)**2 - 
     &   MUE**2*UCStop(1,3)**2*UStop2(1,1)

        Cd(10753) = MUE**2*Cd(10603)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(10603)*UCStopC(1,3)**2 + 
     &   4*MUE2*Cd(10604)*UStop2(1,1)*UStop2(1,2)

        Cd(10761) = 2*Cd(50)*MStop2(2)*
     &    (Cd(10760) + M_3*Cd(9948)*Cr(75)*MStop2(2)) + 
     &   Cd(10086)*Cd(10604)*MStop2(3)

        Cd(10763) = 2*Cd(50)*MStop2(2)*
     &    (Cd(10762) + Cd(9948)*Cr(75)*MStop2(2)) + 
     &   Cd(10604)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(10766) = 2*Cd(50)*MStop2(2)*
     &    (Cd(10765) + M_3C*Cd(9948)*Cr(75)*MStop2(2)) + 
     &   Cd(10039)*Cd(10604)*MStop2(3)

        Cd(10767) = 2*MUE2*Cd(10763)*Cr(5) + 
     &   MUE**2*Cd(10761)*UCStop(1,3) + 
     &   MUEC**2*Cd(10766)*UCStopC(1,3)

        Cd(10785) = 13*Cd(50) + 4*Cd(179)*Cd(251) + 
     &   At*AtC*Cd(251)*UStop2(1,1)

        Cd(10787) = -(Cd(50)*Cd(10342)) + Cd(251)*Cd(10786) + 
     &   2*Cd(9948)*Cd(10785)*UStop2(1,1)

        Cd(10793) = 4*M_3*Cd(10658) - Cd(9948)*Cd(10792) + 
     &   Cd(10209)*UCStop(1,3)**2*UStop2(1,1)

	Cd(10794) = Cd(10659) - M_3C*Cd(9948)*UStop2(1,1)**2

	Cd(10798) = AtC*Cd(10651) + At*(Cd(10653) + At*UCStopC(1,3)**2)

        Cd(10804) = Cd(50)*Cd(9966) - MGlmT2*Cd(10802) + 
     &   Cd(10803)*MStop2(2)

	Cd(1073) = Cd(1003) + Cd(1072)

        Cd(10825) = -(MGlmT2*Cd(10658)) + Cd(9948)*Cd(10824) + 
     &   Cr(75)**2*MStop2(2)*UStop2(1,1)

        Cd(10833) = MGlmT2*Cd(10723) + 
     &   2*Cd(9948)*Cd(10706)*Cr(75) + Cd(10832)*MStop2(2)

        Cd(10839) = MGlmT2*
     &    (-Cd(10836) + Cd(9948)*Cd(10834)*Cr(75)) + 
     &   Cd(10838)*MStop2(2)

        Cd(10845) = MGlmT2*
     &    (-Cd(10842) + Cd(9948)*Cd(10840)*Cr(75)) + 
     &   Cd(10844)*MStop2(2)

        Cd(10846) = -(Cd(13)*Cd(10833)) + 2*AtC*Cd(10839) + 
     &   2*At*Cd(10845)

        Cd(10847) = Cd(50)*Cd(10846) + 
     &   Cd(251)*(2*MTy2*Cd(10831) + Cd(13)*Cd(179)*MGlpTmSt2(2))*
     &    MStop2(3)

        Cd(10857) = Cd(50)*(2*Cd(9948)*Cd(10852) + Cd(10856)) + 
     &   Cd(251)*Cd(10850)*MStop2(3)

        Cd(10868) = Cd(50)*(Cd(9948)*Cd(10862) + Cd(10867)) + 
     &   Cd(251)*Cd(10863)*MStop2(3)

        Cd(10889) = MGl2*Cd(10755) - MGlmT2*Cd(10886) + 
     &   Cd(10888)*MStop2(2)

	Cd(10942) = Cd(11262) + Cd(10918)*MStop2(3)

	Cd(10944) = Cd(11264) + Cd(10924)*MStop2(3)

        Cd(10945) = 2*MUE**2*Cd(10942)*MStop2(1)*
     &    UCStop(1,3)**2 - 
     &   MUEC**2*Cd(10943)*UCStopC(1,3)*UStop2(1,2)*UStopC(1,1) + 
     &   MUE2*Cd(10944)*UStopC(1,2)

	Cd(10971) = Cd(11267) - Cd(10958)*MStop2(3)

	Cd(11002) = Cd(11133) - 4*AtC*UUStop(1,1)*UUStop(2,3)

	Cd(1094) = -Cd(964) + Cd(1034) + 2*Cd(1093)*UStop2(1,1)

        Cd(1096) = 2*Cd(1094)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(1095)*UStop2(1,1) - 
     &      2*MStop2(1)*MStop2(3)*UStop2(1,2))

	Cd(11042) = -Cd(11046) + 5*MTy*MUE*MStop2(3)**2*UUStop(1,3)

        Cd(11048) = 2*(Cd(11047)*MStop2(1)*MStop2(3) + 
     &     2*MTy*MUE*MStop2(3)**2*UStop(1,1)*UStop2(1,2)**3 + 
     &     Cd(11046)*UStopC(1,2))

	Cd(11093) = -Cd(11097) + 5*MTy*MStop2(3)**2*UUStop(1,3)

        Cd(11099) = 2*(Cd(11098)*MStop2(1)*MStop2(3) + 
     &     2*MTy*MStop2(3)**2*UStop(1,1)*UStop2(1,2)**3 + 
     &     Cd(11097)*UStopC(1,2))

	Cd(11103) = Cd(11253) - Cd(11082)*MStop2(3)

	Cd(1103) = -2*Cd(952) + 2*Cd(1021) - Cd(1102)*UStop2(1,1)

        Cd(11128) = Cd(11127)*MStop2(1)*MStop2(3) - 
     &   8*At*MTy*MStop2(3)**2*UStop(1,2)*UStop2(1,2)**3 + 
     &   2*Cd(11126)*UStopC(1,2)

        Cd(1105) = Cd(1101)*MStop2(3) - 2*Cd(1103)*UStop2(1,1) + 
     &   2*Cd(1104)*UStop2(1,2)**2

        Cd(11131) = 2*(Cd(11130)*MStop2(1)*MStop2(3) - 
     &     2*At*MTy*MStop2(3)**2*UStop(1,2)*UStop2(1,2)**3 + 
     &     Cd(11129)*UStopC(1,2))

        Cd(11141) = 2*(Cd(11140)*MStop2(1)*MStop2(3) - 
     &     4*Cd(10980)*MStop2(3)**2*UCStop(1,3)*UStop2(1,2)**2 + 
     &     Cd(11136)*UStopC(1,2))

        Cd(11147) = -4*Cd(11146)*MStop2(1)*MStop2(3) + 
     &   4*Cd(10980)*MStop2(3)**2*UCStop(1,3)*UStop2(1,2)**2 + 
     &   2*Cd(11145)*UStopC(1,2)

	Cd(1107) = 2*Cd(985) - Cd(1071) - 2*Cr(1)*UStop2(1,1)

	Cd(1108) = Cd(986) - Cd(1070) - Cd(1107)*UStop2(1,1)

        Cd(11189) = 2*Cd(11187)*UCStop(3,1)*UStop2(1,1) - 
     &   2*Cd(698)*Cr(91)*UUStopC(2,1) - 4*Cd(11188)*UUStopC(2,3)

	Cd(1112) = 2*Cd(992) - 2*Cd(1080) - Cd(1111)*UStop2(1,1)

        Cd(11218) = 2*MUE**2*Cd(10917)*MStop2(1)*
     &    UCStop(1,3)**2 + 
     &   4*MUEC**2*Cd(10931)*UCStopC(1,3)*UStop2(1,2)*
     &    UStopC(1,1) - MUE2*Cd(10941)*UStopC(1,2)

        Cd(11251) = 4*Cd(9948)*Cd(11202)*MStop2(1) - 
     &   2*Cd(11249)*MStop2(3) - 4*AtC*Cd(11250)*UCStop(1,3)

        Cd(11265) = 2*MUE**2*Cd(11262)*MStop2(1)*
     &    UCStop(1,3)**2 - 
     &   MUEC**2*Cd(11263)*UCStopC(1,3)*UStop2(1,2)*UStopC(1,1) + 
     &   MUE2*Cd(11264)*UStopC(1,2)

        Cd(11271) = 2*Cd(11265) - Cd(11268) + 
     &   2*Cd(11270)*MStop2(3)*UUStopC(1,2)

        Cd(11272) = (Cd(11271) + 
     &     2*MUE*MStop2(3)*UCStop(1,3)*UStop2(1,2)*
     &      (2*MUE*MStop2(3)*UCStop(1,3)*UStop2(1,2) + 
     &        Cd(11261)*UStopC(1,2)))/(MStop2(1)**2*MStop2(3)**2)

        Cd(11312) = (-Cd(11311) + 
     &     4*MUE*MStop2(3)*UCStop(1,3)*UStop2(1,1)*
     &      (MUE*MStop2(3)*UCStop(1,3)*UStop2(1,1) + 
     &        2*Cd(11293)*UStopC(2,2)))/MStop2(3)**2

        Cd(11364) = Cd(11362)*MStop2(3) - 
     &   AtC*Cd(11363)*UCStop(1,3) + 
     &   Cd(11123)*MStop2(2)*UUStopC(3,2)

        Cd(11367) = 4*Cd(11366)*MStop2(2)*MStop2(3) - 
     &   2*Cd(11365)*MStop2(3)**2*UStop(1,1) + 
     &   4*Cd(11039)*MStop2(2)**2*UUStopC(3,2)

        Cd(11380) = MUEC*Cd(698)*Cr(91)*UUStopC(2,1) - 
     &   Cd(11377)*UUStopC(2,2) - Cd(11379)*UUStopC(2,3)

	Cd(1131) = -4*Cd(958) + 2*Cd(1024) - Cd(1128)*MStop2(3)

        Cd(11397) = Cd(11386)*UStop2(1,1)*UStopC(1,1)*
     &    UUStopC(2,2) - 2*Cd(11396)*UUStopC(3,1) + 
     &   Cd(11395)*UUStopC(3,3)

        Cd(11408) = AtC*Cd(11406) + 
     &   2*MTy*MUE*Cd(11348)*MStop2(3) - 
     &   Cd(11407)*MStop2(2)*UUStopC(3,3)

        Cd(11418) = 2*Cd(11417)*UStopC(1,2)*UUStopC(2,1) + 
     &   2*Cd(11408)*UUStopC(3,2) + Cd(11416)*UUStopC(3,4)

        Cd(11431) = 2*Cd(11430) + MUE**2*Cd(970)*UUStop(3,1) - 
     &   MUEC**2*UCStopC(1,3)**2*UUStop(3,2)

	Cd(11439) = Cd(11294) + 2*Cd(11438)

        Cd(11466) = 2*MTy*Cd(11465)*MStop2(2)*UCStopC(1,3) + 
     &   Cd(11462)*UStop(2,2) - 
     &   Cd(11204)*MStop2(2)**2*UCStop(3,1)*UUStopC(1,1)

        Cd(11473) = 4*Cd(11472)*MStop2(2)*UCStopC(3,4) - 
     &   Cd(11469)*UCStop(1,3)*UStop(2,1) - 
     &   4*Cr(91)**2*MStop2(2)**2*UStopC(2,2)*UUStopC(1,2)

	Cd(11512) = Cd(11362) + 2*Cd(11511)

        Cd(11514) = Cd(11512)*MStop2(3) - 
     &   AtC*Cd(11513)*UCStop(1,3) + 
     &   2*Cd(11123)*MStop2(2)*UUStopC(3,2)

        Cd(11520) = 2*Cd(11519)*MStop2(2)*MStop2(3) - 
     &   Cd(11517)*MStop2(3)**2*UStop(1,1) + 
     &   4*Cd(11039)*MStop2(2)**2*UUStopC(3,2)

        Cd(11524) = 4*Cd(11329) - Cd(11475) - 
     &   2*Cd(9809)*MStop2(3)*UStop2(1,2)

        Cd(11531) = 2*MUEC*Cd(11528)*UStopC(2,1) - 
     &   Cd(11529)*UStopC(2,2) - Cd(11530)*UUStopC(2,2)

        Cd(11536) = AtC*Cd(11535) + 
     &   MTy*MUE*Cd(11493)*MStop2(3) - 
     &   Cd(11407)*MStop2(2)*UUStopC(3,3)

        Cd(11548) = Cd(11547)*UStopC(1,2)*UUStopC(2,1) + 
     &   2*Cd(11536)*UUStopC(3,2) + Cd(11546)*UUStopC(3,4)

	Cd(11556) = -4*Cd(970) + 2*Cd(11437) - Cd(11555)

	Cd(11559) = Cd(11304) + 3*Cd(11447)

	Cd(11560) = Cd(11451) + Cd(11559)

	Cd(11561) = Cd(11449) + Cd(11452)

	Cd(11566) = Cd(11433) + Cd(11426)*MStop2(3)

	Cd(11601) = 4*Cd(11578) + AtC*(AtC*Cd(11573) + 4*Cd(11582))

        Cd(11609) = Cd(11386)*UStop2(1,1)*UStopC(1,1)*
     &    UUStopC(2,2) - 2*Cd(11608)*UUStopC(3,1) + 
     &   Cd(11529)*UUStopC(3,3)

        Cd(11615) = Cd(11399)*UStop2(1,2)*UStopC(1,2)*
     &    UUStopC(2,1) + 2*Cd(11614)*UUStopC(3,2) + 
     &   Cd(11613)*UUStopC(3,4)

	Cd(11616) = 4*Cd(11375) - 2*Cd(11527) + Cd(11607)

        Cd(11628) = -(AtC**2*Cd(11627)*UCStop(1,3)**2) + 
     &   At*(AtC*Cd(11626) - At*Cd(11627)*UCStopC(1,3)**2)

        Cd(11635) = 4*(AtC**2*Cd(11634)*UCStop(1,3)**2 + 
     &     At**2*Cd(11634)*UCStopC(1,3)**2 + 
     &     2*Cd(11632)*UStop2(1,1))

        Cd(11644) = Cd(11641)*MStop2(1) + 
     &   MStop2(3)*(-Cd(11643) + Cd(9948)*MStop2(2)*UStop2(1,1))

        Cd(11649) = Cd(11641)*MStop2(1) + 
     &   MStop2(3)*(-Cd(11648) + Cd(9948)*Cr(26)*UStop2(1,1))

        Cd(11674) = MUEC*Cd(11668)*UUStopC(2,1) + 
     &   MUE*Cd(11671)*UUStopC(2,2) + Cd(11673)*UUStopC(2,3)

        Cd(1162) = UStop2(1,1)*
     &    (Cd(1155) - Cd(1161)*UStop2(1,1)) + Cd(1159)*UStop2(1,2)

        Cd(11703) = -2*AtC*MUE*Cd(11699)*UCStop(1,3)**2 + 
     &   2*At*MUEC*Cd(11702)*UCStopC(1,3)**2 + 
     &   4*Cd(424)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(1167) = Cd(1084) - Cd(1140)

        Cd(11746) = AtC*Cd(11744) + 
     &   At*(-Cd(11745) + At*Cd(9972)*UCStopC(1,3)**2)

	Cd(1169) = Cd(1085) + Cd(1147)*MStop2(3)

        Cd(11766) = M_3*Cd(11764) + 
     &   At*(-Cd(11765) + 3*M_3*At*UCStopC(1,3)**2)

        Cd(11775) = Cd(47)*Cd(11774) + 
     &   MGl2*Cd(152)*MStop2(3)*(2*Cd(10198) + Cd(152)*UCStop(1,3))

        Cd(11798) = Cd(9948)*Cd(11797) + 
     &   2*M_3*Cd(11791)*Cr(28) + Cd(11794)*UCStop(1,3)**2

        Cd(11805) = Cd(9948)*Cd(11804) + 
     &   2*M_3C*Cd(11791)*Cr(28) + Cd(11802)*UCStopC(1,3)**2

	Cd(115) = 2*MTy2*Cd(113) + MTy*(Cd(109) - Cd(114)*UStop2(1,1))

        Cd(11826) = 8*MUE2*Cd(11824) - 
     &   8*Cd(11825)*(MUE**2*UCStop(1,3)**2 + 
     &      MUEC**2*UCStopC(1,3)**2)

        Cd(11827) = Cd(47)*Cd(164)*
     &    (Cd(10449) + Cd(9948)*Cd(10450)) + 
     &   Cd(152)*Cd(1201)*MStop2(1)*MStop2(3)

	Cd(11829) = Cd(10449) + Cd(9948)*Cd(10450)

        Cd(11861) = At*Cd(11852) + 
     &   AtC*Cd(11850)*UCStop(1,3)**2 + 2*Cd(11860)*UStop2(1,2)

        Cd(11863) = AtC**2*Cd(9969)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(9970) + At*UCStopC(1,3)**2)

        Cd(11867) = (Cd(50)*Cd(9955)*Cr(75)*
     &      (-(MGlmT2*Cd(10313)) + 2*Cr(5)*MStop2(2)) + 
     &     MGl2*Cd(179)*Cd(10313)*MStop2(3))/(MStop2(2)*MStop2(3))

	Cd(11884) = AtC*Cd(11880) + At*(Cd(11883) + At*UCStopC(1,3)**2)

        Cd(11890) = AtC*Cd(11887) - 
     &   At*(Cd(11889) + At*Cd(10077)*UCStopC(1,3)**2)

	Cd(1182) = -4*Cd(947) + 8*Cd(1019) - Cd(1122)*MStop2(3)

        Cd(11894) = -(MGl2*Cd(138)*Cd(179)**2) + 
     &   Cd(50)*(Cd(50) + Cd(138)*MStop2(2))

	Cd(1183) = -8*Cd(952) + 8*Cd(1021) - Cd(1182)*UStop2(1,1)

        Cd(11913) = 2*M_3*Cd(10658) - Cd(9948)*Cd(11912) + 
     &   Cd(10398)*UCStop(1,3)**2*UStop2(1,1)

        Cd(11919) = -(AtC*Cd(9970)) - 2*Cd(11918) + 
     &   4*M_3C*Cd(9948)*UStop2(1,1)**2

	Cd(11920) = (MGlmT2*Cd(11919) + Cd(12551)*MStop2(2))/MStop2(2)

        Cd(11924) = -(Cd(50)*Cd(251)*Cd(11923)) + 
     &   2*Cd(179)*MStop2(3)*UStop2(1,1)*
     &    (2*Cd(179)*Cd(10224) + Cd(424)*MStop2(2)*UStop2(1,1))

        Cd(11928) = At*Cd(10043) + 
     &   Cd(787)*MStop2(3)*UStop2(1,1) - MGlpT2*Cr(75)*UStop2(1,2)

	Cd(11931) = Cd(12557) - Cd(118)*Cr(75)*UCStop(1,3)*UStop2(1,2)

        Cd(11941) = Cd(50)*Cr(28) + 
     &   Cd(251)*(Cd(10802) + Cd(9948)*Cd(9983)*MStop2(1))

        Cd(11950) = MGl2*Cd(251)*Cd(10749) + MGlmT2*Cd(11949) + 
     &   Cd(11948)*MStop2(2)

	Cd(1189) = Cd(1109) + Cd(1186)

        Cd(11979) = Cd(6096)*Cd(9948) - Cd(50)*Cr(85) + 
     &   Cd(47)*UStop2(1,1)

        Cd(11983) = -(Cd(50)*Cd(9966)) + 
     &   Cd(47)*(Cd(9948)*Cd(9983) + Cr(75)**2*UStop2(1,1))

	Cd(1193) = 2*Cd(992) - 2*Cd(1191) - Cd(1192)*UStop2(1,1)

	Cd(12004) = Cd(50) + Cd(47)*(Cd(9948) - UStop2(1,2))

	Cd(117) = Cd(14038) + Cd(58)*UStop2(1,2)**2

        Cd(1195) = 2*Cd(607)*MStop2(3) + Cd(1193)*UStop2(1,1) + 
     &   Cd(1194)*UStop2(1,2)**2

        Cd(12058) = (-3*M_3C*Cd(10025) + AtC*Cd(12055))*
     &    UCStop(1,3) + Cd(12057)*UCStopC(1,3)

        Cd(12065) = (-3*M_3C*Cd(10025) + AtC*Cd(12063))*
     &    UCStop(1,3) + Cd(12064)*UCStopC(1,3)

	Cd(12111) = Cd(47) - Cd(50)*Cd(9948) + Cd(6096)*UStop2(1,1)

        Cd(12123) = -(MUEC**2*Cd(50)*UCStopC(1,3)**2) + 
     &   MUE2*(10*MStop2(2)*MStop2(3)*UStop2(1,1)**2 - 
     &      Cd(50)*UStop2(1,2)**2)

	Cd(12136) = Cd(1775) - 2*MGl2*Cd(2254)

        Cd(12145) = 2*Cd(11201)*Cd(12072) - 2*Cd(12144) + 
     &   AtC*Cd(12072)*Cr(91)*UCStop(3,3)

        Cd(12146) = MUEC*Cd(12145) + 
     &   Cd(12070)*Cd(12141)*MStop2(3) + 
     &   MUE*Cd(698)*Cd(12072)*UCStop(3,3)

        Cd(12157) = MTy*Cd(12155)*UStop(1,2) - 
     &   Cd(12072)*(2*Cd(12156)*UStopC(2,1) + 
     &      AtC*Cr(58)*UStopC(2,2))

	Cd(12187) = MGl2*Cd(303) - MTy2*Cd(12186) + MStop2(1)*MStop2(2)

        Cd(12199) = MGl2*Cd(12195)*UStop(1,1) - 
     &   Cd(12198)*UUStop(1,1) + 2*At*MTy*Cd(45)*UUStop(1,3)

        Cd(12206) = MGl2*Cd(12203)*UStop(1,1) - 
     &   2*MTy*Cd(1591)*UUStop(1,1) + At*Cd(12205)*UUStop(1,3)

        Cd(12237) = MGl2*Cd(12172) - Cd(12236)*UStop(1,1) - 
     &   2*At*MTy*Cd(45)*UStop(1,2)

        Cd(12264) = Cd(47)*Cd(12254)*MGlpTmSt2(2)*UStopC(1,2)*
     &    UUStopC(2,1) - 2*MUE*Cd(12258)*UUStopC(3,2) - 
     &   Cd(12263)*UUStopC(3,4)

        Cd(12267) = Cd(50)*Cd(12266) - 
     &   Cd(47)*MStop2(2)*MStop2(3) + MGl2*Cd(47)*MStop2(4)

        Cd(12312) = 4*M_3*MTy2*Cd(164)*Cd(2597) + 
     &   At*(2*Cd(12311) + 
     &      Cd(164)*Cd(2597)*Cd(10137)*UStop2(1,1)*UStop2(1,2))

        Cd(12324) = 4*M_3C*Cd(9948)*Cd(12298)*MStop2(1) + 
     &   AtC*Cd(12323)*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

        Cd(12398) = AtC*Cd(11744) - 
     &   At*(Cd(11745) + At*Cd(10077)*UCStopC(1,3)**2)

        Cd(12404) = Cd(10086)*Cd(12371)*MStop2(3) - 
     &   4*Cd(47)*Cd(164)*MStop2(1)*
     &    (Cd(10454) - M_3*Cd(9948)*UStop2(1,1))

        Cd(12408) = At*Cd(12407) + 
     &   M_3*MStop2(1)**2*(2*Cd(11759) + 3*At**2*UCStopC(1,3)**2)

        Cd(1234) = 4*Cd(1233)*UStop2(1,1) - 
     &   Cd(1230)*UStop2(1,1)**2 + Cd(1232)*UStop2(1,2)**2

        Cd(12414) = -(Cd(12411)*MStop2(3)) + 
     &   Cd(47)*(-2*Cd(12413)*MStop2(1) + 
     &      2*Cd(47)*Cd(12405)*UCStop(1,3))

        Cd(12418) = 4*MGlmT2*Cd(47)*Cd(9962) + 
     &   8*M_3*Cd(12370) + At*Cd(12371)*UStop2(1,1)

        Cd(12419) = 4*MGlmT2*Cd(47)*Cd(9963) + 
     &   8*M_3C*Cd(12370) + AtC*Cd(12371)

        Cd(12420) = Cd(12419)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(9948)*Cd(12418) + 
     &      At*Cd(12391)*UStop2(1,2)**2)

	Cd(12423) = Cd(12374) + Cd(12375)*MStop2(1)

        Cd(12427) = 4*MGlmT2*Cd(47)*Cd(9994) + 
     &   8*M_3C*Cd(12370) + AtC*Cd(12371)*UStop2(1,1)

        Cd(12428) = 4*MGlmT2*Cd(47)*Cd(9995) + 
     &   8*M_3*Cd(12370) + At*Cd(12371)

        Cd(12429) = Cd(12428)*UCStopC(1,3)**2 + 
     &   UStop2(1,1)*(Cd(9948)*Cd(12427) + 
     &      AtC*Cd(12391)*UStop2(1,2)**2)

        Cd(12433) = -(At*Cd(9966)) + 4*M_3*Cd(11909)*Cr(28) + 
     &   2*Cd(12432)*UCStop(1,3)**2

        Cd(12437) = (Cd(12371)*Cd(12434) + 
     &     (2*Cd(47)*MStop2(1)*
     &        (-(MGlmT2*Cd(12436)) + Cd(12433)*MStop2(1)))/
     &      MStop2(3))/MStop2(1)**2

        Cd(12439) = -(AtC*Cd(9966)) + 
     &   4*M_3C*Cd(11909)*Cr(28) + 2*Cd(12438)*UCStopC(1,3)**2

        Cd(12443) = (Cd(12371)*Cd(12440) + 
     &     (2*Cd(47)*MStop2(1)*
     &        (-(MGlmT2*Cd(12442)) + Cd(12439)*MStop2(1)))/
     &      MStop2(3))/MStop2(1)**2

        Cd(12451) = (MGl2 - MStop2(1))*UCStop(1,3)*UStop2(1,1) + 
     &   Cd(11817)*UStop2(1,2)

        Cd(12464) = 2*Cd(47)**2*Cd(10436) + 
     &   MStop2(1)*(-4*Cd(47)*Cd(164)*Cd(12463) + 
     &      Cd(1201)*Cd(12371)*MStop2(3))

        Cd(12469) = 2*MUE2*Cd(12467) - 
     &   Cd(12468)*(MUE**2*UCStop(1,3)**2 + 
     &      MUEC**2*UCStopC(1,3)**2)

        Cd(12474) = 2*MUE2*Cd(12471) - 
     &   Cd(12473)*(MUE**2*UCStop(1,3)**2 + 
     &      MUEC**2*UCStopC(1,3)**2)

        Cd(12475) = -4*Cd(47)*Cd(164)*MStop2(1)*
     &    (Cr(28) - Cd(9948)*UStop2(1,1)) + 
     &   Cd(12371)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(12476) = Cd(10039)*Cd(12371)*MStop2(3) - 
     &   4*Cd(47)*Cd(164)*MStop2(1)*
     &    (Cd(10321) - M_3C*Cd(9948)*UStop2(1,1))

        Cd(12477) = 2*MUE2*Cd(12475)*Cr(5) + 
     &   MUE**2*Cd(12404)*UCStop(1,3) + 
     &   MUEC**2*Cd(12476)*UCStopC(1,3)

	Cd(12478) = MGl2*Cd(12469) + MTy2*Cd(12474) + 2*MTy*Cd(12477)

        Cd(12495) = 5*Cd(50)*Cd(138) + Cd(12481) + 
     &   AtC**2*Cd(50)*UCStop(1,3)**2

        Cd(12507) = 2*Cd(50)*
     &    (MGlmT2*Cd(12506) + Cd(50)*Cr(28)) + 
     &   Cd(12483)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(12510) = (Cd(10313)*Cd(12483) + 
     &     (4*Cd(50)*Cd(9955)*Cr(75)*MStop2(2)*
     &        (-(MGlmT2*Cd(10313)) + 2*Cr(5)*MStop2(2)))/MStop2(3))
     &    /MStop2(2)**2

        Cd(12523) = AtC*Cd(11887) + 
     &   At*(-Cd(11889) + At*Cd(9972)*UCStopC(1,3)**2)

        Cd(12539) = 2*Cd(50)*Cd(251)*
     &    (MUEC*Cd(9948)*Cd(9983) + Cd(10674))*MStop2(2) + 
     &   Cd(795)*Cd(12483)*MStop2(3)

        Cd(12540) = 2*Cd(50)*Cd(251)*
     &    (Cd(9948)*Cd(10063) + Cd(10696))*MStop2(2) + 
     &   Cd(10219)*Cd(12483)*MStop2(3)

        Cd(12552) = (Cd(12440)*Cd(12483) + 
     &     (2*Cd(50)*MStop2(2)*
     &        (MGlmT2*Cd(11919) + Cd(12551)*MStop2(2)))/MStop2(3))/
     &   MStop2(2)**2

	Cd(12560) = Cd(12557) - Cd(9955)*Cr(75)*MStop2(1)*UCStop(1,3)

        Cd(1249) = Cd(1248)*UStop2(1,1) + 
     &   UStop2(1,2)*(-MStop2(1) + MStop2(2)*UStop2(1,2))

	Cd(12564) = Cd(12563) - Cd(9955)*Cr(75)*MStop2(1)*UCStopC(1,3)

        Cd(1250) = MUE**2*Cd(1249)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(1249)*UCStopC(1,3)**2 + 
     &   2*MUE2*Cd(1247)*UStop2(1,1)*UStop2(1,2)

        Cd(12573) = 2*Cd(50)**2*Cd(12572) + 
     &   MStop2(2)*(4*Cd(50)*Cd(251)*Cd(12571) + 
     &      Cd(1201)*Cd(12483)*MStop2(3))

	Cd(12575) = MTy2*Cd(1201) + MGl2*Cd(10448) + 2*MTy*Cd(12574)

        Cd(12580) = (-(Cd(50)**2*Cd(118)*Cd(152)*Cd(4095)) + 
     &     Cd(47)**2*(2*Cd(50)**2 - Cd(138)*Cd(179)*Cd(4190)))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12595) = 4*M_3C*Cd(12594)*UCStop(1,3)**2 + 
     &   MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2)*
     &    (4*M_3 - At*UStop2(1,2))

        Cd(12597) = 4*M_3C*Cd(12596)*UCStop(1,3)**2 + 
     &   MGlpTmSt2(2)*UStop2(1,1)*UStop2(1,2)*
     &    (4*M_3 - At*UStop2(1,2))

        Cd(12621) = Cd(47)**2*
     &    (Cd(12620) + 2*MGlpT2*MUE2*Cd(9948)*UStop2(1,1)**2) + 
     &   MUE**2*Cd(50)**2*Cd(12617)*UCStop(1,3)**2*UStop2(1,2)

        Cd(12649) = 2*MTy*Cd(164)*Cd(12648)*MGlpTmSt2(1) + 
     &   Cd(12647)*MGlpTmSt4(1) + Cd(12628)*UUStopC(1,1)

        Cd(12652) = 8*M_3C*Cd(12594)*UCStop(1,3)**2 - 
     &   2*MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2)*
     &    (-4*M_3 + At*UStop2(1,2))

        Cd(1259) = 2*Cd(426)*Cd(1256) + 
     &   2*MTy*Cd(13)*Cd(1238)*MStop2(4) - 
     &   Cd(424)*Cd(1253)*UStop2(1,1)

        Cd(12685) = 2*MUE*MStop2(2)*UCStop(1,3)**2 + 
     &   MUEC*(Cd(12286)*UStop2(1,1) - Cd(12284)*UStop2(1,2))

        Cd(12686) = Cd(50)*Cd(12070)*Cd(12289)*MStop2(3) + 
     &   Cd(47)**2*(At*Cd(12288) + AtC*Cd(12685)*UStop2(1,1))

        Cd(12742) = 1/2.D0*
     &   (A0delStop(2)*Cd(12735) + 4*Cd(12720)*MStop2(3)**2 - 
     &      2*Cd(12741)*Cr(34)*MStop2(3)**2)/MStop2(3)**2

        Cd(12769) = -8*Cd(12701)*Cd(12767)*UCStop(3,1)*
     &    UCStopC(1,3) + Cd(12700)*
     &    (2*Cd(12766)*UStop2(1,2)**2 + 
     &      8*Cd(12768)*MStop2(1)*UStop(1,2)*UUStop(3,1))

        Cd(12790) = Cd(12789)*UStop2(1,1) + 
     &   2*UStop2(1,2)*(Cd(12765)*MStop2(2)*MStop2(3) + 
     &      Cd(12700)*Cd(12785)*UUStopC(3,3))

        Cd(12791) = 2*Cd(12790) + 
     &   Cd(12701)*(Cd(12786) + 
     &      2*Cd(12701)*Cd(12778)*UCStopC(1,3)*UStop2(1,2)*
     &       UStopC(1,1))

        Cd(12806) = -(Cd(12805)*UStopC(1,1)) + 
     &   Cd(12700)*UStop2(1,1)*UStop2(1,2)*
     &    (Cd(12700)*UCStop(1,3)**2 + 2*Cd(12804)*UStopC(1,2))

        Cd(12816) = Cd(12759)*MStop2(3) + 
     &   Cd(12700)*(2*Cd(12808) + Cd(607)*Cd(12700)*UCStop(1,3)**2)

	Cd(12870) = Cd(12868) + Cd(12869) - 4*Cd(12867)*MStop2(3)

        Cd(12876) = -(Cd(12701)*Cd(12875)) + 
     &   Cd(12700)*(Cd(12871) - 
     &      Cd(12700)*Cr(91)*UCStop(1,3)*UUStopC(3,2))

        Cd(12883) = 2*(Cd(12700)*Cd(12877) - 
     &     Cd(12701)*Cd(12881) + Cd(12882)*Cr(91)*MStop2(4))

        Cd(12887) = 2*UStop2(1,1)*
     &    (Cd(12886) - Cd(12885)*UStop2(1,1)) - 
     &   2*MStop2(2)*MStop2(4)*UStop2(1,2)**2

        Cd(1284) = Cd(1283) + Cd(1279)*UStop2(1,1) + 
     &   MStop2(2)**2*UStop2(1,1)**3

        Cd(1286) = Cd(1284)*UStop2(1,1) + Cd(1285)*UStop2(1,2) + 
     &   MStop2(2)**2*UStop2(1,2)**4

        Cd(12948) = -(Cd(12701)*Cd(12947)) + 
     &   Cd(12700)*(Cd(12945) - 
     &      Cd(12700)*Cr(91)*UCStop(1,3)*UUStopC(3,2))

        Cd(1289) = Cd(1287) + 2*Cr(28)*MStop2(3)**2 - 
     &   4*Cd(1288)*UStop2(1,2)

        Cd(12965) = Cd(12964)*UCStopC(1,3) + 
     &   Cd(152)*Cd(12700)*
     &    (MTy*Cd(12700)*UCStop(1,3)*UStop2(1,2)**2 + 
     &      Cd(12963)*UUStopC(1,3))

        Cd(12990) = 1/(256.D0*Pi**4)*
     &   (-Cd(12989) + Cd(12701)*
     &       (256*Pi**4*(Cd(12988)*MStop2(1)*MStop2(3)) - 
     &         4*Cd(11106)*Cd(12701)*MStop2(1)*UCStopC(1,3)*
     &          UUStopC(3,1)))/(MStop2(1)*MStop2(3))

        Cd(1293) = Cd(1291) - 2*Cr(28)*MStop2(3)**2 + 
     &   4*Cd(1292)*UStop2(1,1)

	Cd(12999) = -2*M_3C*MTy*Cd(12993) + Cd(13289)

        Cd(13006) = Cd(13005)*UCStopC(1,3) + 
     &   Cd(179)*Cd(12700)*
     &    (MTy*Cd(12700)*Cr(75)*UCStop(1,3)*UStop2(1,2) - 
     &      Cd(13003)*UStopC(1,1))

        Cd(1294) = 8*MUE2*Cd(1286) - 
     &   2*MUE**2*Cd(1289)*UCStop(1,3)**2 + 
     &   2*MUEC**2*Cd(1293)*UCStopC(1,3)**2

        Cd(13011) = -(Cd(13010)*UCStop(1,3)) + 
     &   Cd(179)*Cd(12701)*
     &    (Cd(13008) - MTy*Cd(12701)*Cr(75)*UCStopC(1,3)*
     &       UStop2(1,2))

        Cd(13013) = Cd(179)*Cr(28) + 2*Cd(13012)*UStop2(1,1) + 
     &   2*MStop2(2)*UStop2(1,2)

        Cd(13068) = Cd(12700)*Cd(13025) + 
     &   Cd(12701)*(-Cd(13067) + 
     &      Cd(12701)*MStop2(3)*UCStopC(1,3)**2)

        Cd(13118) = -(Cd(12700)**2*Cd(13107)*UCStop(1,3)) + 
     &   Cd(12701)*(512*Pi**4*Cd(13115) - 
     &      Cd(12701)*Cd(13110)*UStop(1,2))

        Cd(13119) = 1/(512.D0*Pi**4)*
     &   (-(Cd(13118)*MStop2(1)**2) + 2*MStop2(1)*MStop2(3) + 
     &      4*Cd(12700)*Cd(13117)*UStop2(1,1)*UStop2(1,2)*
     &       UUStopC(1,2))/MStop2(1)**2

        Cd(1307) = 2*Cd(1304)*MStop2(3) + 
     &   AtC*MTy*Cd(1305)*UCStop(1,3) - 2*Cd(1306)*UStop2(1,1)

        Cd(13138) = Cd(13137)/MStop2(3) - 
     &   (2*(3*Cd(13134)*MStop2(1) + Cd(13133)*MStop2(2))*
     &      UCStopC(1,3)*UStopC(1,1))/MStop2(1)**2

        Cd(13142) = (-6*UCStop(1,3))/MStop2(1) + 
     &   (-2*Cd(13140)*UCStop(1,3) + 
     &      2*Cd(13141)*MStop2(1)**2*UStop2(1,1))/
     &    (MStop2(1)**2*MStop2(3))

        Cd(13148) = 2/Pi**4*
     &   (Cd(13147)*MStop2(1)**2 - 3*Cr(28)*MStop2(1)*MStop2(3) + 
     &      Cd(13145)*UStop2(1,1)*UStop2(1,2))/
     &    (MStop2(1)**2*MStop2(3))

        Cd(1310) = Cd(1287) + 2*MStop2(3)**2*UStop2(1,1)**2 + 
     &   Cd(1309)*UStop2(1,2)

        Cd(13222) = -(Cd(12701)**2*Cd(13207)*UCStopC(1,3)) + 
     &   2*Cd(12700)*(Cd(13216) - 
     &      Cd(12700)*Cd(13202)*MStop2(3)*UCStop(1,3)*UStopC(1,2))

	Cd(13225) = Cd(12958) - Cd(12594)*MStop2(3)

        Cd(13230) = Cd(13229)*UCStopC(1,3) + 
     &   Cd(12700)*MGlpTmSt2(1)*
     &    (MTy*Cd(12700)*UCStop(1,3)*UStop2(1,2)**2 + 
     &      Cd(12963)*UUStopC(1,3))

        Cd(13250) = Cd(12700)*
     &    (Cd(13248) - MTy*Cd(12700)*Cd(12810)*UCStop(1,3)) - 
     &   Cd(12701)*Cd(13249)*UCStopC(1,3)

        Cd(1321) = Cd(1319) - 2*Cr(28)*MStop2(3)**2 - 
     &   Cd(1320)*UStop2(1,1)

        Cd(13285) = -(Cd(12954)*MGlpTmSt2(1)*MStop2(3)) + 
     &   Cd(152)*(Cd(12701)*Cd(12956) + 
     &      Cd(12700)*Cd(12957)*UUStopC(1,2))

	Cd(13286) = 2*M_3*Cd(12965) + 2*M_3C*Cd(12971) + Cd(13285)

        Cd(13287) = 1/(256.D0*Pi**4)*
     &   (-2*(Cd(13286) + 128*Pi**4*
     &          (Cd(13284)*MGlpTmStxGlT4(1)*MStop2(1)))*
     &       MGlpTmStxGlT4(1)*MStop2(1) + 
     &      Cd(12836)**2*Cd(13275)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(1)**2*MStop2(1)**2*MStop2(3))

        Cd(13291) = Cd(13005)*UCStopC(1,3) + 
     &   Cd(179)*Cd(12700)*
     &    (MTy*Cd(12700)*Cr(75)*UCStop(1,3)*UStop2(1,2) - 
     &      Cd(13003)*UStopC(1,1))

	Cd(13300) = Cd(12986) + 3*Cd(12700)*UStop2(1,1)**2*UStop2(1,2)

        Cd(13348) = 1/(512.D0*Pi**4)*
     &   (512*Pi**4*(Cd(13340)*MGlpTmStxGlT4(1)**2*MStop2(1)**2) + 
     &      512*Pi**4*(TB2*Cd(13347)*MGlpTmStxGlT4(1)*MStop2(1)**2*
     &         MStop2(3)) - Cd(12371)*Cd(12836)**2*Cr(45)*MStop2(3)
     &      )/(TB2*MGlpTmStxGlT4(1)**2*MStop2(1)**2*MStop2(3))

	Cd(1331) = Cd(1265) - Cd(1328) - 4*Cd(1330)*UStop2(1,1)**2

        Cd(1332) = 2*Cd(1329)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(1331)*UStop2(1,1) + 
     &      2*MStop2(2)*MStop2(3)*UStop2(1,2))

        Cd(13390) = (1/2.D0*
     &      (Cd(13385) + 8*Cd(13389)*MStop2(2))/MStop2(2)**2 + 
     &     (2*Cd(13387))/MStop2(3))*UCStopC(1,3)

        Cd(13396) = 2*UCStop(1,3)*
     &   ((Cd(13393)*UCStop(1,3) + 
     &        Cd(13395)*MStop2(2)**2*UStop2(1,1))/
     &      (MStop2(2)**2*MStop2(3)) - 
     &     (2*UCStop(1,3)*UStop2(1,2)**2)/MStop2(2))

        Cd(13403) = 4/Pi**4*
     &   (Cd(13400)*MStop2(2)**2 + Cd(13402)*MStop2(2)*MStop2(3) + 
     &      Cd(13399)*UStop2(1,1)*UStop2(1,2))/
     &    (MStop2(2)**2*MStop2(3))

        Cd(1338) = Cd(1334)*MStop2(3) - Cd(1336)*UStop2(1,1) + 
     &   2*Cd(1337)*UStop2(1,2)**2

        Cd(13447) = -(Cd(13446)*UCStopC(1,3)) + 
     &   Cd(12700)*MGlpTmSt2(2)*
     &    (MTy*Cd(12700)*Cr(75)*UCStop(1,3)*UStop2(1,2) - 
     &      Cd(13003)*UStopC(1,1))

        Cd(13455) = Cd(12701)*Cd(13453) - 
     &   Cd(12700)*Cd(13025)*MGlpTmSt2(2) - Cd(13454)*MStop2(3)

        Cd(13461) = -(Cd(13446)*UCStopC(1,3)) + 
     &   Cd(12700)*MGlpTmSt2(2)*
     &    (MTy*Cd(12700)*Cr(75)*UCStop(1,3)*UStop2(1,2) - 
     &      Cd(13003)*UStopC(1,1))

        Cd(1340) = MTy*Cd(13)*Cd(1339)*MStop2(3) + 
     &   AtC*MUE*Cd(1338)*UCStop(1,3)**2 + 
     &   At*MUEC*Cd(1338)*UCStopC(1,3)**2

        Cd(13476) = MTy*Cd(13475) + Cd(13474)*UStop2(1,2) + 
     &   AtC*Cd(13469)*UStopC(1,2)

        Cd(13595) = A0delq(MSbot2(1),Q)*Cd(13589) - 
     &   A0delStop(2)*Cd(13593)*MSq2Diff(1,1) - 
     &   Cd(13594)*MSq2Diff(1,2)

        Cd(13596) = Cd(13589)*MSbot2(1) - 
     &   Cd(13591)*MSq2Diff(1,2)*MStop2(1) - 
     &   Cd(13593)*MSq2Diff(1,1)*MStop2(2)

        Cd(13599) = Cd(13466)*Cd(13589)*MSbot2(1) - 
     &   Cd(13598)*MSq2Diff(1,2) - 
     &   Cd(13593)*Cr(34)*MSq2Diff(1,1)*MStop2(2)

        Cd(1363) = Cd(1361) - 4*Cr(28)*MStop2(3)**2 + 
     &   Cd(1362)*UStop2(1,1)

        Cd(13708) = -2*MSq2Diff(1,2)*MStop2(2) + 
     &   Cd(179)*(Cd(13707) - MTy*Cd(12700)*UCStop(1,3))

        Cd(13713) = -2*Cd(13712)*UStop(1,1) + 
     &   MStop2(2)*(Cd(48)*Cd(13710) - 
     &      MGl2*MTy*Cd(1909)*UStop(1,2))

        Cd(13733) = 2*Cd(13729)*MSbot2(1) - 
     &   Cd(13725)*MSbot2(1)**2 - Cd(13732)*MStop2(2)

        Cd(13742) = Cr(56)*MSbot2(1)**2 - 
     &   At*MUE*TB*MSq2Diff(1,2) - Cd(13722)*MSbot2(1)*UStop2(1,1)

        Cd(13745) = -2*Cd(13744)*MSbot2(1)*MStop2(2) + 
     &   Cd(13743)*MStop2(2)**2 + Cd(13742)*UStop2(1,1)

        Cd(13750) = -2*Cd(13747)*MSbot2(1) + 
     &   2*Cr(87)*MSbot2(1)**2 + Cd(13749)*MStop2(2)

        Cd(13751) = MGl2*(-(MGl2*Cd(13741)) + 2*Cd(13745)) + 
     &   Cd(13750)*MStop2(2)

        Cd(13764) = Cd(47)*MStop2(1) - Cd(13763)*UStop2(1,1) + 
     &   Cd(13760)*UStop2(1,2)

	Cd(1371) = Cd(1265) - Cd(1328) - 4*Cd(1370)*UStop2(1,1)**2

        Cd(1372) = 2*Cd(1369)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(1371)*UStop2(1,1) + 
     &      2*MStop2(2)*MStop2(3)*UStop2(1,2))

        Cd(13787) = At*(AtC*Cd(13785) + MUE*TB*Cd(13785)) - 
     &   Cd(13786)*MSq2Diff(1,1)**2

        Cd(13796) = -2*Cd(10367)*Cd(13792)*MSbot2(1) + 
     &   2*Cd(13795)*MStop2(1) + Cd(13793)*MStop2(1)**2

	Cd(135) = 2*MTy2*Cd(133) + MTy*(Cd(129) + Cd(134)*UStop2(1,1))

        Cd(13822) = Cd(13560)*Cd(13807) - Cd(13820) - 
     &   Cd(13469)*Cd(13821)*MSbot2(1)*MSq2Diff(1,2)

        Cd(13859) = At*AtC*Cd(10367) - Cd(13788)*MStop2(1) + 
     &   MSbot2(1)*UStop2(1,1)

	Cd(13861) = At*AtC*Cd(10367) - Cd(13798)

	Cd(1381) = -Cd(1304) + 2*Cd(1380)

        Cd(13885) = 2*Cd(13469)*Cd(13884)*MSq2Diff(1,2) + 
     &   MSq2Diff(1,1)*(Cd(13883)*UStopC(1,1) + 
     &      Cd(13882)*UStopC(1,2))

        Cd(13898) = Cd(13897)*UStop2(1,1) + 
     &   UStop2(1,2)*(MSq2Diff(1,1) + 2*Cd(13896)*UStop2(1,2))

        Cd(13906) = MTy2*Cd(13902) + MGl2*Cd(13905) + 
     &   2*TB2*Cd(424)*MSq2Diff(1,2)*UStop2(1,1)

	Cd(1385) = 2*Cd(1302) - Cd(1378) - 2*MTy2*Cd(1384)*UStop2(1,1)

        Cd(13910) = MTy2*Cd(13909) + MGl2*TB*MSq2Diff(1,1) - 
     &   2*TB2*Cd(424)*MStop2(1)*UStop2(1,1)

        Cd(13931) = -(MGlpTmSt2(2)*MSq2Diff(1,2)*MStop2(2)) + 
     &   MTy2*Cd(251)*(Cd(13707) - MTy*Cd(12700)*UCStop(1,3))

        Cd(13961) = MGlmT2*MGlpT2*Cd(13698) + 
     &   2*MSbot2(1)*(-Cd(13960) + Cd(1496)*MSbot2(1))

	Cd(1392) = Cd(1308) + Cd(1389) - 2*Cd(1391)*UStop2(1,2)

        Cd(13976) = -(Cd(13975)*MSbot2(1)) + 
     &   Cr(56)**2*MSbot2(1)**2 + Cd(11343)*MStop2(2)**2

        Cd(13982) = -(Cd(13981)*MStop2(2)) + 
     &   2*MSbot2(1)*(Cd(13980) - Cr(56)*MSbot2(1)*UStop2(1,1))

        Cd(13988) = -2*Cd(13985)*MSbot2(1) + 
     &   Cd(13983)*MSbot2(1)**2 + Cd(13987)*MStop2(2)

        Cd(13991) = -(Cd(13990)*MSbot2(1)) + 
     &   Cd(13971)*MStop2(2)**2 + Cr(56)*MSbot2(1)**2*UStop2(1,1)

        Cd(13992) = 2*MGl2*Cd(13991) + 
     &   (MGl2**2*Cd(13737) + Cd(13988))*MStop2(2)

	Cd(137) = Cd(344) + Cd(71)*UStop2(1,2)**2

        Cd(14028) = MGl2*(MGl2*Cd(14023) + 2*Cd(14027)) - 
     &   Cd(14025)*MStop2(1)

        Cd(14044) = 2*Cd(14038)*MSbot2(1)*
     &    (-Cd(14043) + MSbot2(1)) + Cd(14042)*MStop2(1)

	Cd(1401) = 4*Cd(1312) - 2*Cd(1398) + Cd(1400)*UStop2(1,1)**2

        Cd(1402) = 4*Cd(1313)*MStop2(3) + Cd(1399)*UStop2(1,1) + 
     &   Cd(1401)*UStop2(1,2)

        Cd(14110) = AtC*Cd(14107)*UCStop(1,3) + 
     &   MUE*TB*Cd(14108)*UStopC(1,2) + MTy*Cd(14109)*UUStopC(3,4)

        Cd(14124) = Cd(14123)*UStop2(1,1) + 
     &   UStop2(1,2)*(Cd(164)*Cd(9950)*MSq2Diff(1,1) + 
     &      Cd(47)*Cd(13896)*MGlpTmSt2(1)*UStop2(1,2))

        Cd(14128) = At*MGlpT2*Cd(14125) + MTy2*Cd(14126) + 
     &   MGl2*Cd(14127)

        Cd(1409) = Cd(1316) - Cd(1403) - 
     &   6*Cd(1408)*MStop2(2)*UStop2(1,1)

        Cd(14146) = (Cd(9950)*Cd(14145) + 
     &      Cd(47)*Cd(9981)*Cd(13896)*MGlpTmSt2(1))*MSq2Diff(1,1)+
     &     Cd(47)*Cd(9950)*Cd(13841)*UStop2(1,2)

        Cd(14156) = AtC*MUEC*TB2*Cd(14142) + 
     &   At*Cd(14125)*Cd(14142) + TB*Cd(14155)

        Cd(14167) = (2*At*Cr(5)*Cr(75)*MSq2Diff(1,1))/
     &    MStop2(3) - (Cd(14142)*
     &      (MStop2(1)*UCStop(1,3) + 2*At*Cr(5)*UStop2(1,2)))/
     &    Cd(47)

        Cd(14194) = Cd(9948)*Cd(14188) - 
     &   Cd(14193)*MGlpTmSt2(2) + 
     &   TB*Cd(251)*Cd(14189)*MSq2Diff(1,2)

        Cd(14206) = -(MUE2*TB2*MGlpTmSt2(2)) + 
     &   Cd(50)*(Cd(9948) - UStop2(1,2))

	Cd(1417) = Cd(1329) + Cd(1414)*MStop2(3)

        Cd(14224) = -(Cd(50)*Cd(14223)) + 
     &   MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)*
     &    (M_3C*UCStop(1,3)**2 + M_3*Cd(9948)*UStop2(1,1))

        Cd(14229) = -(Cd(50)*Cd(14228)) + 
     &   MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)*
     &    (M_3*UCStopC(1,3)**2 + M_3C*Cd(9948)*UStop2(1,1))

	Cd(1418) = 4*Cd(1262) - 8*Cd(1324) - Cd(1357)*MStop2(3)

        Cd(14282) = -(Cd(14281)*UCStopC(1,3)) + 
     &   MTy*Cd(14279)*UStop2(1,1) - MGl2*Cd(14278)*UStopC(1,1)

	Cd(14287) = MTy2*Cd(14284) + (MGl2 - MSbot2(1))*MSq2Diff(1,2)

	Cd(1426) = -8*Cd(832) + 4*Cd(1269) - Cd(1421)*MStop2(3)

        Cd(1427) = 2*Cd(874)*MStop2(3) - Cd(1425)*UStop2(1,1) + 
     &   Cd(1426)*UStop2(1,2)**2

	Cd(14326) = Cd(14306) - Cd(13469)*MSq2Diff(1,2)

	Cd(14332) = Cd(14323) + Cd(14331)*MSq2Diff(1,2)

        Cd(1445) = 4*Cd(1444) + 5*Cr(28)*MStop2(3) - 
     &   2*Cd(1442)*UStop2(1,1)**2

	Cd(1447) = 4*Cd(1286) - Cd(1354) + Cd(1438)*MStop2(3)

	Cd(1448) = 2*Cd(1289) - 2*Cd(1356) + Cd(1441)*MStop2(3)

	Cd(1452) = Cd(1369) + Cd(1451)*MStop2(3)

        Cd(1455) = 2*Cd(1452)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(1454)*UStop2(1,1) + 
     &      4*MStop2(2)*MStop2(3)*UStop2(1,2))

	Cd(1463) = Cd(1384) + Cd(1456)

	Cd(1470) = Cd(1308) + Cd(1389) - Cd(1467)*MStop2(3)

	Cd(1477) = Cd(1400) - Cd(1421)

	Cd(1484) = Cd(1408) + MStop2(3)

        Cd(1507) = Cd(52)*Cd(1506) + 
     &   2*Cd(1500)*T134(MGl2,MTy2,MStop2(1)) + 
     &   2*Cd(1493)*Cd(1495)*MStop2(3)**2*UStop2(1,2)

	Cd(1602) = MGl2*(2*MGl2 - Cd(1591)) + Cd(1593)

	Cd(1626) = MGl2*(-6*MGl2 + Cd(1619)) + 2*Cd(1621)

	Cd(1649) = MGl2*(3*MGl2 - Cd(1646)) + 3*Cd(1648)

        Cd(1679) = Cd(44)**3*Cd(48) + 
     &   MGl2*(-(Cd(44)*Cd(1635)) - MGl2*Cd(1678))

	Cd(1681) = Cd(44)**2 + MGl2*(MGl2 - Cd(1680))

        Cd(1697) = Cd(44)**3*Cd(89) + 
     &   MGl2*(-(Cd(44)*Cd(1691)) + MGl2*Cd(1696))

	Cd(1708) = -(MGl2*Cd(1707)) + Cd(1775)

        Cd(1725) = Cd(48)**3*Cd(89) + 
     &   MGl2*(-(Cd(48)*Cd(1719)) + MGl2*Cd(1724))

        Cd(1844) = 3*Cd(96)*Cd(1708)*MStop2(3)**2*UCStop(1,3) + 
     &   MTy*(4*M_3*Cd(50)*MGlpTmSt2(2)*MStop2(3) - 
     &      2*M_3*Cd(1779)*UStop2(1,2))

	Cd(1879) = MGl2*(-2*MGl2 + Cd(1503)) + 2*Cd(1820)

        Cd(1958) = 2*MGl2**2*Cd(1952) - MTy2*Cd(1957) + 
     &   Cd(1955)*MStop2(3)

	Cd(1967) = Cd(2320) - 2*MTy2*Cd(1964)*UStop2(1,2)**2

        Cd(1978) = Cd(1976) + 
     &   MGl2*(-(MTy2*Cd(1977)) + Cd(1954)*MStop2(3))

	Cd(2004) = Cd(2001) - MGl2*(MTy2*Cd(2003) + Cd(2002)*MStop2(3))

        Cd(2005) = Cd(48)**2*Cd(1997) - MGl2*Cd(2004) + 
     &   2*MTy2*Cd(50)*MGlpTmSt2(2)*UStop2(1,2)

        Cd(2033) = 2*MGl2**2*Cd(1952) - MTy2*Cd(2032) - 
     &   Cd(2029)*MStop2(3)

	Cd(2043) = Cd(2532) - 2*MTy2*Cd(2039)*UStop2(1,2)**2

        Cd(2052) = 2*(Cd(2051) + 
     &     MStop2(2)*(-2*Cd(2049) + 3*Cd(45)*MStop2(2)))

	Cd(2062) = Cd(2300) - Cd(1607)*UStop2(1,2)

        Cd(2076) = Cd(44)**3*Cd(2063) - MGl2*Cd(2075) - 
     &   Cd(2072)*UStop2(1,2)

        Cd(2096) = 2*Cd(48)**3*Cd(2048) - 2*MGl2*Cd(2095) - 
     &   Cd(2091)*UStop2(1,2)

	Cd(2097) = Cd(2451) - Cd(1526)*UStop2(1,2)

        Cd(2109) = Cd(48)**3*Cd(1701) + MGl2*Cd(2108) - 
     &   Cd(2105)*UStop2(1,2)

        Cd(212) = Cd(201)*MStop2(2) + Cd(205)*UStop2(1,1) + 
     &   Cd(210)*UStop2(1,2)

	Cd(2143) = MGl2*(-4*MGl2 + Cd(2138)) + 2*Cd(2139)

        Cd(2272) = -(M_3*Cd(2271)*MStop2(3)*UStop(1,2)) + 
     &   MTy*(Cd(2269)*UStop(1,1) + 
     &      2*Cd(2266)*MStop2(3)*UStop2(1,2))

        Cd(2291) = Cd(48)**2*Cd(2279) + 
     &   MGl2*(-Cd(2290) + MGl2**2*MStop2(3))

        Cd(2292) = -(M_3*Cd(2291)*UStop(1,2)) + 
     &   MTy*(Cd(2289)*UStop(1,1) + 2*Cd(2285)*UStop2(1,2))

        Cd(2303) = 2*M_3*MTy2*Cd(86)*Cd(167) + 
     &   MTy*(-(Cd(2300)*UCStop(1,3)) + 2*Cd(2302)*UUStopC(1,2))

        Cd(228) = Cd(97)*MStop2(2) - Cd(227)*UStop2(1,1) - 
     &   Cd(226)*UStop2(1,2)

        Cd(2317) = Cd(44)**4 + 
     &   MGl2*(-4*Cd(44)**2*Cd(45) + MGl2*Cd(2316))

	Cd(229) = Cd(211) + Cd(228)*UStopC(1,1)

        Cd(2402) = Cd(44)**2*Cd(1984) + 
     &   MGl2*(MTy2*Cd(1988) - Cd(1989)*MStop2(3))

	Cd(238) = Cd(239) + MTy*Cd(164)*Cr(5)

        Cd(2420) = -(Cd(2419)*MStop2(1)) + 
     &   MTy2*(Cd(2416) + 4*MTy2*UStop2(1,1))

        Cd(2434) = -(Cd(2433)*MStop2(1)) + 
     &   MGl2*(Cd(290)*Cd(1583) + 8*MGl2*UStop2(1,2))

        Cd(2439) = 2*MGl2*Cd(2438) + 
     &   6*M_3*Cd(2364)*MStop2(1) - 
     &   M_3C*Cd(2437)*UCStop(1,3)**2

        Cd(2454) = 2*M_3*MTy2*Cd(96)*Cd(254) + 
     &   MTy*(Cd(2451)*UCStop(1,3) - 2*Cd(2453)*UUStopC(1,2))

        Cd(2467) = Cd(48)**2*Cd(2007) + 
     &   MGl2*(MTy2*Cd(2011) + Cd(2012)*MStop2(3))

        Cd(2488) = -(Cd(2487)*MStop2(2)) + 
     &   MTy2*(Cd(2486) + 4*MTy2*UStop2(1,1))

        Cd(2495) = 2*MGl2*Cd(2494) + 
     &   6*M_3*Cd(2364)*MStop2(2) + 
     &   M_3C*Cd(1052)*MStop2(3)*UCStop(1,3)**2

        Cd(2496) = Cd(2364)*MStop2(2)**3*MStop2(3) + 
     &   Cd(1996)*UStop2(1,1)**3

        Cd(2544) = M_3*Cd(2532)*MStop2(3) - 
     &   MTy*(Cd(2537)*UCStop(1,3) + 2*Cd(2543)*UUStopC(1,2))

	Cd(252) = Cd(344) + MTy*Cd(251)*Cr(5)

        Cd(2566) = -(Cd(970)*MStop2(2)**3*MStop2(3)) + 
     &   Cd(2022)*UStop2(1,1)**3

        Cd(2606) = MGl2*(Cd(44)*Cd(2599) - MGl2*Cd(2605)) - 
     &   2*Cd(44)**2*MStop2(1)**2

        Cd(2665) = MTy*(M_3C*Cd(152)*Cd(2597)*UUStop(1,1) - 
     &      M_3*Cd(2662)*UUStop(1,2)) + 2*Cd(2664)*UUStop(1,3)

        Cd(2689) = MGl2*(Cd(48)*Cd(2683) - MGl2*Cd(2688)) - 
     &   2*Cd(48)**2*MStop2(2)**2

        Cd(2741) = MTy*(M_3C*Cd(179)*Cd(2681)*UUStop(1,1) - 
     &      M_3*Cd(2738)*UUStop(1,2)) - 2*Cd(2740)*UUStop(1,3)

        Cd(2840) = -2*Cd(2839)*UCStop(1,3) + 
     &   MTy*(4*M_3*Cd(2757)*MStop2(1) - 
     &      M_3*Cd(2801)*UStop2(1,2))

        Cd(2855) = 2*Cd(2852)*MStop2(1) + 
     &   MGl2*(Cd(2850) - 2*MGl2*Cd(2854)*MStop2(2))

        Cd(2862) = Cd(2856) + 
     &   MGl2*(2*Cd(2861)*MStop2(1) - 2*Cd(1789)*MStop2(2))

        Cd(2906) = Cd(2903) + 
     &   MGl2*(2*Cd(2905)*MStop2(1) - 2*Cd(1789)*MStop2(2))

        Cd(3086) = MTy*(2*M_3C*Cd(3063)*UUStop(1,1) - 
     &      M_3*Cd(3083)*UUStop(1,2)) + 
     &   2*Cd(3085)*MStop2(3)*UUStop(1,3)

        Cd(3155) = MGl2*Cd(3154) + 
     &   MStop2(2)*(MGl2**4 - Cd(48)**2*Cd(1701)*MStop2(2))

        Cd(3175) = MGlmT2*(3*MGlmT2**2 - Cd(2617)*MStop2(1)) - 
     &   Cd(3174)*MStop2(2)

        Cd(3188) = -3*Cd(48)**2*MStop2(2)**2 + 
     &   MGl2*(Cd(3187) + 6*MGl2**2*MStop2(2))

        Cd(3191) = MGlmT2*(3*MGlmT2**2 + MTy2*MStop2(1)) - 
     &   Cd(3190)*MStop2(2)

        Cd(3255) = MGl2*(Cd(3251) - MGl2*Cd(3254)) + 
     &   Cd(48)**2*Cd(1701)*MStop2(2)

        Cd(3260) = MGl2*Cd(3257) + 
     &   MStop2(2)*(MGl2**3*Cd(3258) - Cd(3259)*MStop2(2))

        Cd(3280) = MGl2*Cd(3279) + 
     &   MStop2(2)*(3*MGl2**4 - Cd(48)**2*Cd(2168)*MStop2(2))

        Cd(328) = Cd(47)*MGlpTmSt2(2) - 
     &   Cd(50)*MGlpTmSt2(1)*UStop2(1,2) - Cd(327)*UStopC(1,1)

        Cd(3340) = MGl2*Cd(3338) + 
     &   2*MStop2(2)*(-(MGl2**3*Cd(3258)) + Cd(3339)*MStop2(2))

        Cd(330) = Cd(47)*MGlpTmSt2(2) + 
     &   2*M_3C*MTy*Cd(295)*UCStop(1,3) - 
     &   Cd(50)*MGlpTmSt2(1)*UStop2(1,2)

        Cd(3422) = MGlmT2*(MGlmT2**2 + MTy2*MStop2(1)) + 
     &   Cd(3421)*MStop2(2)

        Cd(3455) = MGl2*(Cd(3289) - MGl2*Cd(3454)) + 
     &   Cd(48)**2*Cd(2168)*MStop2(2)

	Cd(341) = Cd(44)**2 + MTy*Cd(339) - MGl2*Cd(340)

        Cd(3460) = MTy*(2*M_3C*MTy2*Cd(251)*MStop2(3)*
     &       UUStop(1,1) + M_3*Cd(3457)*UUStop(1,2)) - 
     &   2*Cd(3459)*MStop2(3)*UUStop(1,3)

        Cd(3498) = Cd(3496)*MStop2(1) - 
     &   MStop2(2)*(Cd(3497) + Cd(2314)*MStop2(2))

        Cd(3528) = Cd(3526)*MStop2(1) + 
     &   MStop2(2)*(2*Cd(3527) + Cd(2314)*MStop2(2))

	Cd(351) = Cd(273) + Cd(350)*UStopC(1,1)

        Cd(3561) = 2*(Cd(3559)*MStop2(1) + 
     &     MStop2(2)*(2*Cd(3560) + Cd(45)*Cd(1592)*MStop2(2)))

        Cd(3566) = Cd(3563)*MStop2(1) + 
     &   MStop2(2)*(-2*Cd(3564) + Cd(3565)*MStop2(2))

        Cd(3575) = Cd(529)*MStop2(1) + 
     &   MStop2(2)*(-2*Cd(1890) + MStop2(2))

        Cd(3585) = Cd(3583)*MStop2(1) + 
     &   MStop2(2)*(2*Cd(3584) - Cd(2314)*MStop2(2))

        Cd(3597) = 2*Cd(3595)*MStop2(1) + 4*Cd(3596)*MStop2(2) - 
     &   2*Cd(45)*Cd(1592)*MStop2(2)**2

        Cd(3602) = MGl2*(Cd(3597) - MGl2*Cd(3601)) - 
     &   Cd(44)**2*Cd(2115)*MStop2(1)

        Cd(3615) = 4*MTy2*Cd(1019) + 
     &   6*MStop2(4)*(2*MTy2**2 + MStop2(1)*MStop2(4))

        Cd(3624) = Cd(3622)*MStop2(1) - 
     &   MStop2(2)*(2*Cd(3623) + Cd(3565)*MStop2(2))

        Cd(3625) = MGl2*(Cd(3597) + 2*MGl2*Cd(3624)) - 
     &   Cd(44)**2*Cd(2115)*MStop2(1)

        Cd(3626) = MTy*(2*M_3C*Cd(3063)*UUStop(1,1) - 
     &      M_3*Cd(3083)*UUStop(1,2)) + 
     &   2*Cd(3101)*MStop2(3)*UUStop(1,3)

        Cd(3668) = Cd(3667)*MStop2(1) + 
     &   MStop2(2)*(Cd(3497) + Cd(2314)*MStop2(2))

        Cd(3713) = MGl2*MGlmT2*Cd(2617) + 
     &   MStop2(1)*(-Cd(3712) + 3*MGl2*MStop2(1))

        Cd(3725) = Cd(3721) - 
     &   MGl2*(2*Cd(3724)*MStop2(3)*UUStop(1,2) + 
     &      Cd(2315)*MStop2(1)*UUStop(2,2))

        Cd(3740) = MGl2*(Cd(3737) + MGl2*Cd(3739)) - 
     &   Cd(44)**2*Cd(3736)*MStop2(1)

        Cd(3765) = MGl2*(2*Cd(3762) + 2*MGl2**2*Cd(3763)) + 
     &   Cd(3764)*MStop2(2)

        Cd(3777) = MGl2*(Cd(48)*Cd(3773) + MGl2*Cd(3775)) + 
     &   Cd(3776)*MStop2(2)

        Cd(3803) = MGl2*(Cd(3799) - MGl2*Cd(3802)) + 
     &   Cd(48)**2*Cd(1847)*MStop2(2)

	Cd(3812) = MGl2*(Cd(3799) - MGl2*Cd(3810)) + Cd(3811)*MStop2(2)

        Cd(3829) = MGl2*(Cd(3824) - 2*MGl2*Cd(3828)) - 
     &   2*Cd(48)**2*Cd(2048)*MStop2(2)

        Cd(3838) = MGl2*(Cd(3834) - MGl2*Cd(3837)) - 
     &   Cd(48)**2*Cd(2048)*MStop2(2)

        Cd(3864) = MGl2*(Cd(3857) - MGl2*Cd(3863)) - 
     &   Cd(48)**2*Cd(2115)*MStop2(2)

        Cd(37) = Cd(32)*dMf133eps(-1) - Cd(36)*dMf133eps(0) - 
     &   Cd(34)*dMf133eps(1)

	end
