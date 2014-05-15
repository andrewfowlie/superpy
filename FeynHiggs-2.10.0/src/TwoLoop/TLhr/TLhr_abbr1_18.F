* TLhr_abbr1_18.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_18
	implicit none

#include "TLhr.h"


        Cd(14231) = At*Cd(14230) + 
     &   Cd(251)*MSq2Diff(1,1)*MSq2Diff(1,2)**2 + 
     &   2*M_3C*AtC*(-(Cd(50)*MSq2Diff(1,1)) + 
     &      (2*Cd(50) + MGlpTmSt2(2)*MSq2Diff(1,1))*MSq2Diff(1,2))*
     &    UCStop(1,3)**2 + 
     &   Cd(9948)*Cd(14227)*MSq2Diff(1,1)*UStop2(1,1) + 
     &   MUE2*TB2*Cd(50)*MSq2Diff(1,1)*UStop2(1,1)**2 - 
     &   2*M_3*AtC*Cd(50)*MSq2Diff(1,2)*UStop2(1,1)**2 - 
     &   MUE2*TB2*MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)*
     &    UStop2(1,1)**2 + 
     &   MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)**2*
     &    UStop2(1,1)**2 + 
     &   2*Cd(50)*MSq2Diff(1,1)*MSq2Diff(1,2)*UStop2(1,2) - 
     &   2*M_3*AtC*Cd(50)*MSq2Diff(1,1)*UStop2(1,1)*
     &    UStop2(1,2) - MUE2*TB2*MGlpTmSt2(2)*MSq2Diff(1,1)*
     &    MSq2Diff(1,2)*UStop2(1,1)*UStop2(1,2) - 
     &   2*M_3*AtC*Cd(50)*MSq2Diff(1,2)*UStop2(1,2)**2 - 
     &   MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)**2*UStop2(1,2)**2

        Cd(5427) = 2*MGl2*Cd(5418) + 
     &   MGl2**2*(-2*Cd(5409) + 
     &      MGl2*(2*Cd(5398) + 
     &         MGl2*(MGl2*(MGl2*(12*MGl2 - Cd(5391)) + 
     &               2*Cd(5393)) - 2*Cd(5423)))) - 
     &   2*MTy2*Cd(5426) + 
     &   2*(-(MTy2**3*Cd(5390)) + 
     &      Cd(1272)*(MStop2(1)**3*MStop2(2)**3 + 
     &         MTy2**2*(MTy2**2*Cd(5387) + 
     &            Cd(5385)*MStop2(1)*MStop2(2))) - 
     &      MTy2**5*(Cd(5400) + MTy2*(6*MTy2 - 19*MStop2(4)))) + 
     &   Cd(5348)*UStop2(1,2)**2

        Cd(7344) = -14*MTy2**7 - 3*MTy2**5*Cd(7309) + 
     &   MTy2**4*Cd(1272)*Cd(7311) - MTy2**3*Cd(7314) + 
     &   MGl2*(MGl2*(-Cd(7327) + 
     &         MGl2*(Cd(7321) + 
     &            MGl2*(MGl2*
     &                (MGl2*(6*MGl2 - Cd(7315)) + Cd(7316)) - 
     &               Cd(7330)))) + Cd(7340)) - MTy2*Cd(7343) + 
     &   3*MTy2**2*Cd(1272)*Cd(7307)*MStop2(1)*MStop2(2) + 
     &   Cd(1272)*MStop2(1)**3*MStop2(2)**3 + 
     &   43*MTy2**6*MStop2(4) + Cd(7295)*UStop2(1,2)**2

        Cd(6091) = 8*MTy*Cd(6090) + 
     &   2*(MGl2**2*Cd(6057) - 
     &      MGl2**3*(MGl2*(Cd(6054) + 
     &            MGl2*(MGl2*(MGl2 - Cd(6050)) + 3*Cd(6058))) + 
     &         Cd(6063)) + 
     &      MTy2*(Cd(6082) - 
     &         MTy2*(MTy2*(-Cd(6068) + 
     &               MTy2*(MTy2*
     &                   (MTy2*(MTy2 - Cd(6051)) + Cd(6072)) - 
     &                  Cd(6075))) + Cd(6086))) - 
     &      MGl2*Cd(6071)*MStop2(1)*MStop2(2) + 
     &      Cd(6060)*MStop2(1)**2*MStop2(2)**2)*UStop(1,1) + 
     &   Cd(6049)*UStop2(1,2)

        Cd(7839) = 4*MTy*Cd(6090) + 
     &   (3*MGl2**2*Cd(7815) + 
     &      MTy2*(MTy2*(Cd(7824) + 
     &            MTy2*(MTy2*
     &                (Cd(7806) - 
     &                  MTy2*(MTy2*(MTy2 - Cd(7802)) + Cd(7812)))+
     &                 Cd(7831))) + Cd(7838)) - 
     &      MGl2**3*(Cd(7808) + 
     &         MGl2*(-Cd(7804) + 
     &            MGl2*(3*Cd(7813) + MGl2*(MGl2 - 7*Cr(1))))) - 
     &      MGl2*Cd(7810)*MStop2(1)*MStop2(2) + 
     &      Cd(7816)*MStop2(1)**2*MStop2(2)**2)*UStop(1,1) + 
     &   Cd(7801)*UStop2(1,2)

        Cd(8004) = -(M_3*MTy2**2*Cd(6287)) + 
     &   MTy2**3*(M_3*Cd(6273) + 
     &      MTy2*(-(M_3*Cd(6264)) + MTy2*Cd(7966))) + 
     &   M_3*MTy2*Cd(7967) + 
     &   M_3*(MGl2**6 - 3*MGl2**5*Cd(843) - 
     &      MGl2**3*Cd(1272)*Cd(6259) - MTy2**6*Cd(6260) - 
     &      3*MGl2*Cd(843)*MStop2(1)**2*MStop2(2)**2 + 
     &      MStop2(1)**3*MStop2(2)**3 + 
     &      3*MGl2**2*Cd(831)*(MGl2**2 + MStop2(1)*MStop2(2)))*
     &    MStop2(3) - 2*M_3C*Cd(7906)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(7965) + Cd(7936)*UStop2(1,1)) + 
     &   MTy*Cd(8003)*UCStop(1,3)*UStop2(1,2)

        Opt(2) = Cd(12580)/MTy2 - Cd(12583)/Cd(47)**2 - 
     &   Cd(12584)/Cd(50)**2 + 
     &   ((-2*Cd(12582) + MTy*Cd(12586))/(Cd(47)**2*Cd(50)**2) + 
     &      MTy2*Cd(12589))/MTy

        Cd(10571) = (-(Cd(47)*Cd(10531)) + 
     &     4*Cd(10570)*MStop2(3))/(Cd(47)**2*MStop2(1)*MStop2(3))

        Cd(10769) = (Cd(10748) + Cd(10768))/
     &   (Cd(50)**2*MStop2(2)**2*MStop2(3))

	Cd(10770) = SA2*Cd(10665) + CA2*Cd(10769)

        Cd(10771) = 1/(128.D0*Pi**3)*
     &   (Alfas2L*(-(S2A*Cd(10728)) + 
     &        Cd(50)**2*Cd(10770)*MStop2(2)**2*MStop2(3)))/
     &    (SB2*SW2*Cd(50)**2*MStop2(2)**2*MStop2(3))

	Cd(10894) = 2*Cd(10237) + Cd(10465)*Cr(13)

        Cd(10900) = 1/(128.D0*Pi**3)*
     &   (-(Alfas2L*Cd(10893)*Cr(34)) + 
     &      256*Pi**3*(SB2*SW2*Cd(10771)))/(SB2*SW2)

	Cd(11051) = Cd(11009) - Cd(11050)*UStopC(1,2)

        Cd(11197) = (Cd(11196) + 
     &     4*MUEC*Cd(9794)*MStop2(3)**2*UCStopC(1,3)*
     &      UStop2(1,1)**2*UStopC(1,1))/(MStop2(1)**2*MStop2(3)**2)

	Cd(11230) = -(S2A*Cd(11197)) + CA2*Cd(11229)

        Cd(11231) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*(SA2*Cd(11214) + Cd(11230)*MStop2(1)**2))/
     &    (SB2*SW2*MStop2(1)**2)

        Cd(11246) = (-Cd(11245) + 
     &     4*MUEC*Cd(9794)*MStop2(3)**2*UCStopC(1,3)*
     &      UStop2(1,1)**2*UStopC(1,1))/(MStop2(1)**2*MStop2(3)**2)

	Cd(11273) = -(S2A*Cd(11246)) + CA2*Cd(11272)

        Cd(11274) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*(-(SA2*Cd(11260)) + Cd(11273)*MStop2(1)**2))/
     &    (SB2*SW2*MStop2(1)**2)

	Cd(11422) = SA2*Cd(11374) + S2A*Cd(11421)

        Cd(11423) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*(Cd(11422) + 
     &        CA2*Cd(11314)*MStop2(2)**2*MStop2(3)**2))/
     &    (SB2*SW2*MStop2(2)**2*MStop2(3)**2)

        Cd(11570) = (-4*Cd(11569)*MStop2(2) + 
     &     Cd(11565)*MStop2(3)**2)/(MStop2(2)**2*MStop2(3)**2)

        Cd(11619) = 2*Cd(11618)*MStop2(2) - 
     &   Cd(12)*Cd(13)*Cd(11616)*MStop2(3)**2

	Cd(11620) = SA2*Cd(11606) + S2A*Cd(11619)

        Cd(11621) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*(Cd(11620) + 
     &        CA2*Cd(11570)*MStop2(2)**2*MStop2(3)**2))/
     &    (SB2*SW2*MStop2(2)**2*MStop2(3)**2)

        Cd(11940) = (4*(-2*Cd(13)*Cd(11894) + 
     &       MTy*Cd(11939)*MStop2(2)))/(MTy*Cd(50)**2*MStop2(2))

	Cd(12291) = -(M_3C*Cd(12217)) + M_3*Cd(12252)

        Cd(12292) = 2*Cd(47)**2*Cd(12291) + 
     &   Cd(12070)*Cd(12149)*UStopC(1,1)*UUStopC(1,2)

	Cd(13358) = 2*Cd(13357) - Cd(13104)*Cr(34)

	Cd(13925) = Cd(13924)*Cr(65) + 1/2.D0*(A0delGl*Cd(13922))

        Cd(14004) = 1/2.D0*(A0delStop(2)*Cd(14003)) - 
     &   Cd(13957)*Cr(34)*MStop2(2)

        Cd(14095) = 1/2.D0*(A0delStop(1)*Cd(14067)) + 
     &   Cd(14094)*Cr(13)*MStop2(1)

	Cd(14113) = 2*Cd(13925) + Cd(14095)

	Cd(14114) = -(A0delq(MSbot2(1),Q)*Cd(14112)) + TB2*Cd(14113)

        Cd(14115) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(2*Cd(14004) + 
     &        MTy2*Cd(50)*Cd(14114)*MSq2Diff(1,2)**2*MStop2(2)))/
     &    (MTy2*SW2*TB2*Cd(50)*MSq2Diff(1,2)**2*MStop2(2))

	Cd(14118) = 2*MGl2*Cd(13924) + Cd(13922)*Cr(65)

	Cd(14119) = 2*Cd(13957) - Cd(14003)*Cr(34)

	Cd(14120) = 2*Cd(14094) + Cd(14067)*Cr(13)

        Cd(14184) = Cd(52)*Cd(14141) + 
     &   Cd(14183)*T134(MGl2,MTy2,MStop2(1))

        Cd(1837) = (Cd(47)**2*Cd(1785) + Cd(50)**2*Cd(1836))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(1840) = Cd(1837) + Cd(1729)*UStop2(1,1)

	Cd(2296) = MTy*Cd(2112)*UStop2(1,2) + Cd(2295)*UStopC(1,1)

	Cd(2297) = M_3C*Cd(2047) + Cd(2296)*UCStopC(1,3)

        Cd(2593) = (Cd(50)**2*Cd(2449) + Cd(47)**2*Cd(2592))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(2594) = 2*Cd(2593) + AtC*Cd(2298)*UCStop(1,3)

	Cd(2820) = Cd(2819)*MStop2(3) + Cd(2811)*UStop2(1,1)*UStopC(1,1)

	Cd(2844) = 2*MUE2*Cd(2820) + Cd(2843)*MStop2(3)

	Cd(3059) = Cd(3058)*MStop2(3) + Cd(3032)*UUStopC(1,1)

	Cd(3111) = Cd(424)*Cd(2820) + Cd(3110)*MStop2(3)

	Cd(3200) = 2*MUE2*Cd(3183) - Cd(3199)*MStop2(3)

	Cd(324) = 2*Cd(323) + AtC*Cd(297)*UCStop(1,3)

	Cd(325) = Cd(324) + At*Cd(284)*UCStopC(1,3)

	Cd(337) = CA*Cd(325) + SA*Cd(336)

	Cd(3488) = Cd(424)*Cd(3183) + Cd(3487)*MStop2(3)

	Cd(354) = CA*Cd(325) + SA*Cd(336)

	Cd(3660) = 2*MTy*Cd(3578)*UStop2(1,2) + Cd(3659)*UStopC(1,1)

	Cd(3661) = 2*M_3C*Cd(3651)*MStop2(3) - Cd(3660)*UCStopC(1,3)

	Cd(3748) = Cd(3729) + 2*Cd(3747)

        Cd(3749) = Cd(3748)*MStop2(3) + 
     &   4*MTy2*Cd(3673)*UStop2(1,1)**2*UStopC(1,1)

	Cd(3750) = Cd(3749)*MStop2(1) - Cd(3710)*UUStopC(1,1)

	Cd(3753) = SA2*Cd(3551) + S2A*Cd(3635)

	Cd(3896) = Cd(3895) + 2*AtC*Cd(3818)*MStop2(2)*UCStop(1,3)

	Cd(393) = -(SA*Cd(325)) + CA*Cd(336)

	Cd(397) = -(SA*Cd(325)) + CA*Cd(336)

	Cd(4106) = ((4*Cd(4104) + MTy*Cd(4105))*MStop2(3)**2)/MTy

	Cd(4113) = SA2*Cd(4041) + S2A*Cd(4112)

        Cd(4114) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(CA2*Cd(4106) + Cd(4113)*MStop2(3)**2))/
     &    (SB2*SW2*Cd(47)**2*MStop2(3)**2)

	Cd(4214) = SA2*Cd(4136) + S2A*Cd(4213)

        Cd(4215) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(CA2*Cd(4201) + Cd(4214)*MStop2(3)**2))/
     &    (SB2*SW2*Cd(50)**2*MStop2(3)**2)

        Cd(5753) = 2*(Cd(5752) - 
     &     2*M_3C*MTy*MTy2*Cd(5638)*UCStop(1,3)*UStop(1,1))

	Cd(5755) = Cd(5753) + 4*MTy2*Cd(5575)*UStop2(1,1)

	Cd(5756) = Cd(5755)*UStop2(1,1) + 2*MTy2*Cd(5595)*UStop2(1,2)**2

	Cd(5757) = -(MTy*Cd(5756)) + Cd(5633)*MStop2(3)

        Cd(5991) = Cd(5953)*UStop2(1,1) + 
     &   2*Cd(5984)*MStop2(3)*UStop2(1,2)

        Cd(6257) = Cd(302)*Cd(6256)*MStop2(3) + 
     &   M_3C*MTy*Cd(6156)*UUStop(1,1)

	Cd(6984) = 4*Cd(6983) - 2*AtC*Cd(6728)*UCStop(1,3)

	Cd(7556) = Cd(7554) + 4*MTy2*Cd(7480)*UStop2(1,1)

	Cd(7753) = M_3C*Cd(7630) + Cd(7752)*UCStopC(1,3)

	Cd(8273) = 2*Cd(8261) - MTy*Cd(8272)*UStop2(1,1)

	Cd(8306) = Cd(8290) - MTy*Cd(8305)*UStop2(1,2)

	Cd(8471) = 2*MGl2**5*Cd(7527) + Cd(8470)*MStop2(3)

        Cd(8472) = 2*MTy2*
     &   (Cd(8471)*UStop(1,1) - 
     &     4*M_3*MTy*Cd(7547)*MStop2(3)**2*UStop(1,2))

	Cd(8475) = Cd(8472) + 4*MTy2*Cd(8454)*UStop2(1,1)

	Cd(8477) = Cd(8476)*MStop2(3) - MTy*Cd(8475)*UStop2(1,1)

	Cd(8648) = 2*Cd(8647) - AtC*Cd(8641)*MStop2(3)*UCStop(1,3)

	Cd(8779) = Cd(8778) + 2*AtC*Cd(8701)*UCStop(1,3)

	Cd(8885) = Cd(8778) + AtC*Cd(8701)*UCStop(1,3)

	Cd(9013) = Cd(9012) + AtC*Cd(8933)*UCStop(1,3)

	Cd(9115) = 2*Cd(9026) + AtC*Cd(8933)*UCStop(1,3)

        Cd(9121) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*Cd(4560)*Cr(13) + 
     &      128*Pi**3*(SB2*SW2*Cd(47)**2*Cd(4114)*MStop2(3)**2))/
     &    (SB2*SW2*Cd(47)**2*MStop2(3)**2)

	Cd(9168) = MTy*Cd(2112)*UStop2(1,2) + Cd(9167)*UStopC(1,1)

	Cd(9169) = M_3C*Cd(2047) + Cd(9168)*UCStopC(1,3)

	Cd(10276) = -(SA2*Cd(10257)) + S2A*Cd(10268) + CA2*Cd(10275)

	Cd(10600) = -(SA2*Cd(10518)) - S2A*Cd(10571) + CA2*Cd(10599)

        Cd(10895) = A0delStop(1)*Cd(10600) + 
     &   A0delStop(2)*Cd(10893) - 2*Cd(10894)*MStop2(1)

        Cd(11053) = (Cd(11051) + 
     &     MUEC*Cd(9794)*MStop2(3)**2*UCStopC(1,3)*UStop2(1,1)**2*
     &      UStopC(1,1) - Cd(11052)*MStop2(3)**2*UUStopC(1,1))/
     &   MStop2(3)**2

	Cd(11157) = CA2*Cd(10975) + S2A*Cd(11053) - SA2*Cd(11156)

	Cd(11959) = SA2*Cd(11893) + S2A*Cd(11940) + CA2*Cd(11958)

	Cd(12480) = -(SA2*Cd(12417)) + S2A*Cd(12462) + CA2*Cd(12479)

	Cd(12579) = SA2*Cd(12526) + S2A*Cd(12570) + CA2*Cd(12578)

	Cd(12687) = Cd(12291) + Cd(12650)

        Cd(13359) = Cd(13339) - Cd(13349)*MStop2(1) + 
     &   Cd(13358)*MStop2(2)

        Cd(14117) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(2*Cd(14003) + 
     &        MTy2*Cd(50)*(2*MGl2*TB2*Cd(13922) + Cd(14116))*
     &         MSq2Diff(1,2)**2))/
     &    (MTy2*SW2*TB2*Cd(50)*MSq2Diff(1,2)**2)

        Cd(14121) = -2*Cd(13466)*Cd(14112)*MSbot2(1) + 
     &   TB2*(2*Cd(14118) + Cd(14120)*MStop2(1))

	Cd(1489) = S2A*Cd(1429) + SA2*Cd(1450) + CA2*Cd(1488)

        Cd(1841) = (-(Cd(47)**2*Cd(1676)) + 
     &     Cd(50)**2*(Cd(1655) + Cd(47)**2*Cd(1840)*UStop2(1,1)))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(2679) = SA2*Cd(2625) + CA2*Cd(2672) + S2A*Cd(2678)

	Cd(2755) = SA2*Cd(2703) + CA2*Cd(2748) + S2A*Cd(2754)

	Cd(3095) = 2*At*AtC*Cd(2820) + 4*Cd(3059) + Cd(3094)*MStop2(3)

	Cd(3112) = SA2*Cd(2844) + CA2*Cd(3095) + S2A*Cd(3111)

	Cd(338) = -(Cd(50)*Cd(250)) + Cd(47)*Cd(279) + 2*MGl2*Cd(337)

	Cd(3489) = SA2*Cd(3200) + CA2*Cd(3469) + S2A*Cd(3488)

	Cd(3898) = Cd(424)*Cd(3795) + 2*MUE*Cd(3896) + 2*MUEC*Cd(3897)

	Cd(4009) = 2*SA2*Cd(3820) + S2A*Cd(3898) + CA2*Cd(4008)

	Cd(394) = Cd(50)*Cd(391) + Cd(47)*Cd(392) - 2*MGl2*Cd(393)

	Cd(4917) = SA2*Cd(4638) + CA2*Cd(4879) - S2A*Cd(4916)

        Cd(5283) = 2*MUE2*Cd(5223) - 
     &   MUE**2*Cd(5282)*UCStop(1,3) - 
     &   MUEC**2*Cd(5012)*UCStopC(1,3)

        Cd(5992) = 2*Cd(5990)*UStop2(1,1) + 
     &   Cd(5918)*UStop2(1,2) - Cd(302)*Cd(5991)*UStopC(1,1)

        Cd(6258) = 2*M_3*Cd(6230)*MStop2(3) + 
     &   MTy*(-(Cd(6237)*UCStop(1,3)) + 2*Cd(6257)*UUStopC(1,2))

        Cd(7231) = 2*MUE2*Cd(7168) - 
     &   MUE**2*Cd(7230)*UCStop(1,3) - 
     &   MUEC**2*Cd(7047)*UCStopC(1,3)

        Cd(7557) = Cd(7522)*MStop2(3) - 
     &   MTy*(Cd(7556)*UStop2(1,1) + 
     &      2*MTy2*Cd(7439)*UStop2(1,2)**2)

        Cd(8478) = 2*MTy*MTy2*Cd(8419)*MStop2(3) + 
     &   Cd(8474)*MStop2(3)**2 - Cd(8477)*UStopC(1,1)

	Cd(8654) = 2*MUEC*Cd(8638) + MUE*Cd(8648) - Cd(8653)*MStop2(3)

	Cd(8781) = Cd(424)*Cd(8687) + 2*MUE*Cd(8779) + 2*MUEC*Cd(8780)

        Cd(8886) = 2*AtC*Cd(8885) + 
     &   At*(Cd(8884) + At*Cd(8783)*MStop2(3)*UCStopC(1,3))

        Cd(9014) = -(MUEC*Cd(8983)) + 2*MUE*Cd(9013) + 
     &   At*MUEC*Cd(8939)*MStop2(3)*UCStopC(1,3)

        Cd(9116) = -2*AtC*Cd(9115) + 
     &   At*(Cd(9114) - At*Cd(9017)*MStop2(3)*UCStopC(1,3))

        Cd(9445) = Cd(52)*Cd(9268) - 
     &   4*MTy2*Cd(1493)*Cd(9186)*MStop2(3)**2 + 
     &   Cd(9444)*T134(MGl2,MTy2,MStop2(1))

        Cd(9713) = Cd(64)*Cd(9539) - 
     &   4*MTy2*Cd(1508)*Cd(9462)*MStop2(3)**2 - 
     &   Cd(9712)*T134(MGl2,MTy2,MStop2(2))

        Cd(10048) = (2*MGlpT2*Cd(13)*Cd(10011) + 
     &     2*MUE*Cd(47)**2*Cd(50)**2*Cd(10035)*MStop2(3) + 
     &     MUEC*Cd(47)**2*Cd(50)**2*Cd(10047)*MStop2(3))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(10187) = (MGl2*Cd(9950)*Cd(10127) - 
     &     2*MTy*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10156) + 
     &     MTy2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10186) + 
     &     Cd(47)**2*Cd(50)**2*Cd(10132)*MStop2(1))/
     &   (Cd(47)**2*Cd(50)**2*Cd(9950))

        Cd(12293) = Cd(12292)*MStop2(3) + 
     &   Cd(12146)*UUStopC(1,1) + Cd(12154)*UUStopC(1,2) + 
     &   Cd(12167)*UUStopC(1,3)

        Cd(356) = A0delGl*Cd(337) + Cd(50)*Cd(343)*Cr(13) - 
     &   Cd(47)*Cd(353)*Cr(34) - Cd(354)*Cr(65)

        Cd(399) = A0delGl*Cd(393) + Cd(50)*Cd(395)*Cr(13) + 
     &   Cd(47)*Cd(396)*Cr(34) - Cd(397)*Cr(65)

        Cd(5758) = MTy*Cd(5427)*MStop2(3) - 
     &   M_3C*Cd(5754)*UCStop(1,3) - 
     &   2*M_3*MTy2*Cd(5516)*UCStopC(1,3)*UStop2(1,2) + 
     &   Cd(5757)*UStopC(1,1)

        Cd(7558) = MTy*Cd(7344)*MStop2(3) - 
     &   M_3C*Cd(7555)*UCStop(1,3) + 
     &   2*M_3*MTy2*Cd(7418)*UCStopC(1,3)*UStop2(1,2) + 
     &   Cd(7557)*UStopC(1,1)

        Cd(8599) = (2*Cd(8478) - MTy*Cd(8598) + 
     &      At**2*MTy*Cd(8397)*MStop2(3)*UCStopC(1,3))/
     &    (MTy*Cd(47)**2*Cd(50)**2*MStop2(3)) - 
     &   AtC**2*Cd(8559)*UCStop(1,3)*UStopC(1,2)

        Cd(10069) = (-(MTy2*Cd(10068)) + 
     &     MUE*Cd(47)**2*Cd(50)**2*Cd(9978)*MStop2(3)**2 + 
     &     MGl2*Cd(47)**2*Cd(50)**2*Cd(9992)*MStop2(3)**2 + 
     &     MUEC*Cd(47)**2*Cd(50)**2*Cd(10008)*MStop2(3)**2 + 
     &     MTy*Cd(47)**2*Cd(50)**2*Cd(10048)*MStop2(3)**2)/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3)**2)

        Cd(10188) = (Cd(10187) + 
     &     (4*MTy2**2*Cd(10103)*Cd(10109) - 
     &        2*Cd(47)*Cd(50)*Cd(10103)*Cd(10120) + 
     &        Cd(47)**2*Cd(10115)*MStop2(2))/
     &      (Cd(47)**2*Cd(50)**2*Cd(10103)) + 
     &     (2*Cd(10014)*Cd(10111))/MStop2(3))/(MZ2*SW2)

        Cd(3751) = AtC*(-4*Cd(3632) + 
     &      AtC*Cd(3550)*UCStop(1,3)) + 
     &   At*(2*AtC*Cd(3523) - 4*Cd(3661) + 
     &      2*At*Cd(3648)*MStop2(3)*UCStopC(1,3))

        Cd(6094) = 2*MTy*Cd(6093) + 
     &   UStop2(1,1)*(Cd(6092)*MStop2(3) - 
     &      (Cd(6091) + Cd(6048)*UStop2(1,1))*UStopC(1,1)) - 
     &   2*Cd(47)*Cd(50)*Cd(5220)*UStop2(1,2)*UStopC(1,2)

        Cd(7015) = 2*MUEC*Cd(6598) - MUE*Cd(6984) + 
     &   AtC*MUEC*Cd(7014) + 
     &   At*(2*MUE*Cd(6678) - MUEC*Cd(6378)*UCStopC(1,3))

        Cd(7842) = MTy*Cd(7841) + 
     &   UStop2(1,1)*(Cd(7840)*MStop2(3) - 
     &      (Cd(7839) + Cd(7800)*UStop2(1,1))*UStopC(1,1)) - 
     &   Cd(47)*Cd(50)*Cd(5220)*UStop2(1,2)*UStopC(1,2)

        Cd(14233) = (TB*(Cd(14231) + 
     &        (MGl2*Cd(14215) + Cd(14219)*MSbot2(1))*MSq2Diff(1,1))
     &       *MStop2(3) - Cd(50)*MSq2Diff(1,2)**2*
     &      (2*TB*Cd(14177) + 
     &        TB2*(MUE*Cd(14226) + MUEC*Cd(14232))*MStop2(3)))/
     &   (TB*TB2*Cd(50)*MSq2Diff(1,1)*MSq2Diff(1,2)**2*MStop2(3))

        Cd(355) = 2*MGl2*Cd(354) + 
     &   Cd(50)*(-2*Cd(343) + Cd(250)*Cr(13)) + 
     &   Cd(47)*(2*Cd(353) - Cd(279)*Cr(34)) - 2*Cd(337)*Cr(65)

        Cd(398) = -2*MGl2*Cd(397) + 
     &   Cd(50)*(2*Cd(395) - Cd(391)*Cr(13)) + 
     &   Cd(47)*(2*Cd(396) - Cd(392)*Cr(34)) + 2*Cd(393)*Cr(65)

        Cd(8274) = MTy*Cd(47)**2*Cd(8260)*MStop2(3) + 
     &   2*M_3*(-(Cd(44)*Cd(48)*Cd(6782)) + MGl2*Cd(6798))*
     &    UCStopC(1,3) + MTy*
     &    (-(Cd(44)**2*Cd(48)**2*Cd(6755)) + MGl2*Cd(6793))*
     &    UStop2(1,1) - Cd(8273)*UStopC(1,1)

	Cd(12607) = MTy2*Cd(12606) + Opt(2)

	Cd(12661) = MTy2*Cd(12660) - Opt(2)

        Cd(5993) = MTy*Cd(5992) - 
     &   2*M_3*(Cd(44)*Cd(48)*Cd(5923) + 
     &      MGl2*(-2*Cd(5932) + 
     &         MGl2*(2*Cd(5926) + 
     &            MGl2*(-Cd(5922) + 
     &               2*MGl2*(Cd(5919) + MGl2*(MGl2 - Cd(5933)))))))
     &     *MStop2(3)*UCStopC(1,3)

        Cd(14263) = (-(MTy*MTy2*TB*Cd(14185)*MSq2Diff(1,2)**2*
     &        MStop2(2)) + 
     &     Cd(10103)*(MGl2*MTy*MTy2*Cd(14204) + 
     &        MTy*MTy2*TB*Cd(14216) + 
     &        MSq2Diff(1,2)*
     &         (MTy*MTy2*TB*Cd(14208) - 
     &           MTy2*Cd(50)*Cd(14148)*Cd(14217) + 
     &           TB*Cd(50)*
     &            (MTy*Cd(14200) + 
     &              MTy2**2*TB2*(MTy*Cd(14233) + Cd(14262)))*
     &            MSq2Diff(1,2)) - 
     &        MTy*MTy2*TB2*(MUE*Cd(14198) + MUEC*Cd(14199))*
     &         UStop2(1,1)))/
     &   (MTy*MTy2*TB*TB2*Cd(50)*Cd(10103)*MSq2Diff(1,2)**2)

        Cd(6330) = -(MTy2**3*
     &      (MTy2*(-2*M_3*Cd(6264) + MTy2*Cd(6267)) + 
     &        2*M_3*Cd(6273))) + 2*M_3*MTy2*Cd(6279) + 
     &   2*M_3*MTy2**2*Cd(6287) - 
     &   2*M_3*(MGl2**6 - MTy2**6*Cd(6260) + 
     &      MStop2(1)**3*MStop2(2)**3 + 
     &      3*MGl2**2*Cd(831)*(MGl2**2 + MStop2(1)*MStop2(2)) - 
     &      Cd(1272)*(MGl2**3*Cd(6259) + 
     &         3*MGl2*(MGl2**4 + MStop2(1)**2*MStop2(2)**2)))*
     &    MStop2(3) - UStop2(1,1)*
     &    (Cd(6258) + Cd(6211)*UStop2(1,1)) + 
     &   2*UCStop(1,3)*(M_3C*Cd(6161)*UCStop(1,3) + 
     &      MTy*Cd(6329)*UStop2(1,2))

        Cd(12664) = AtC*(Cd(50)**2*Cd(152)*
     &       (-4*M_3C*MGlpTmSt2(1)*UCStop(1,3)**2 - 
     &         4*M_3*Cd(9948)*MGlpTmSt2(1)*UStop2(1,1)) + 
     &      Cd(47)*(Cd(50)**2*
     &          (-4*M_3C*UCStop(1,3)**2 - 
     &            4*M_3*Cd(9948)*UStop2(1,1)) + 
     &         Cd(47)*(Cd(50)*
     &             (-4*M_3C*UCStop(1,3)**2 - 
     &               4*M_3*Cd(9948)*UStop2(1,1)) + 
     &            Cd(179)*(-4*M_3C*MGlpTmSt2(2)*UCStop(1,3)**2 - 
     &               4*M_3*Cd(9948)*MGlpTmSt2(2)*UStop2(1,1))))+
     &        AtC*Cd(47)*Cd(50)**2*UCStop(1,3)**2*
     &       (-2*UStop2(1,1) - 2*UStop2(1,2))) + 
     &   Cd(47)*(Cd(47)*Cd(50)*
     &       (Cd(138)*Cd(11996) + 2*Cd(50)*Cd(12661)) + 
     &      Cd(50)**2*Cd(118)*
     &       (-8*Cd(9948)*UStop2(1,1) + 
     &         (5 - 4*UStop2(1,1))*UStop2(1,1) + 
     &         (5 - 4*UStop2(1,2))*UStop2(1,2))) + 
     &   At*(-2*Cd(12663)*UCStopC(1,3)**2 - 
     &      4*M_3C*Cd(50)**2*Cd(152)*Cd(9948)*MGlpTmSt2(1)*
     &       UStop2(1,1) + 
     &      Cd(47)*(-4*M_3C*Cd(50)**2*Cd(9948)*UStop2(1,1) + 
     &         Cd(47)*(-4*M_3C*Cd(50)*Cd(9948)*UStop2(1,1) - 
     &            4*M_3C*Cd(179)*Cd(9948)*MGlpTmSt2(2)*
     &             UStop2(1,1))) + 
     &      AtC*Cd(47)*(Cd(50)**2*
     &          (-4*Cd(9948)*UStop2(1,1)**2 - 
     &            4*UStop2(1,1)*UStop2(1,2)**2) + 
     &         Cd(47)*Cd(50)*
     &          (Cd(9948)*(-8*UStop2(1,1)**2 - 2*UStop2(1,2)**2) + 
     &            UStop2(1,1)*
     &             (-6*UStop2(1,2)**2 + 
     &               UStop2(1,1)*(-2*UStop2(1,1) + 2*UStop2(1,2))))
     &         ))

        Cd(10896) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*Cd(10276) + 64*Pi**3*(SB2*SW2*Cd(10771)*Cr(34)))/
     &    (SB2*SW2)

        Cd(10897) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*Cd(10895) - 
     &      256*Pi**3*(SB2*SW2*Cd(10896)*MStop2(2)))/(SB2*SW2)

	Cd(10898) = Cd(10600)*MStop2(1) + Cd(10893)*MStop2(2)

	Cd(10899) = 2*Cd(10465) - Cd(10600)*Cr(13)

        Cd(10901) = 1/(128.D0*Pi**3)*
     &   (128*Pi**3*(SB2*SW2*Cd(10900)*MStop2(2)) + 
     &      Alfas2L*Cd(10899)*MStop2(1))/(SB2*SW2)

        Cd(12688) = 2*Cd(47)**2*Cd(12687) + 
     &   Cd(12070)*Cd(12149)*UStopC(1,1)*UUStopC(1,2)

	Cd(12693) = 2*Cd(12337) - Cd(12480)*Cr(13)

	Cd(12694) = 2*Cd(12368) - Cd(12579)*Cr(34)

	Cd(12697) = 2*Cd(12480) + Cd(11834)*Cr(13)

	Cd(12698) = 2*Cd(12579) + Cd(11959)*Cr(34)

        Cd(14122) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(-2*Cd(14119) + 
     &        MTy2*Cd(50)*Cd(14121)*MSq2Diff(1,2)**2))/
     &    (MTy2*SW2*TB2*Cd(50)*MSq2Diff(1,2)**2)

        Cd(14264) = Cd(64)*Cd(14195) + 
     &   Cd(14263)*T134(MGl2,MTy2,MStop2(2))

	Cd(14354) = Cd(14184) + Cd(14264)

        Cd(1842) = (Cd(47)*Cd(50)*Cd(1841)*UStopC(1,1) - 
     &     Cd(1839)*UStop2(1,2)*UStopC(1,2))/(Cd(47)*Cd(50))

	Cd(357) = A0delStop(1)*Cd(50)*Cd(250) - 2*Cd(356)*MStop2(1)

        Cd(358) = A0delStop(2)*Cd(47)*Cd(279)*MStop2(1) - 
     &   Cd(357)*MStop2(2)

	Cd(3752) = 4*Cd(3750) - Cd(3751)*MStop2(1)

        Cd(3754) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(CA2*Cd(3752) - Cd(3753)*MStop2(1)))/
     &    (SB2*SW2*Cd(47)**2*MStop2(1)**2*MStop2(3)**2)

	Cd(4010) = Cd(3112)*Cr(13) - 2*Cd(2679)*MStop2(3)**2

	Cd(4011) = Cd(3489)*Cr(34) + 2*Cd(2755)*MStop2(3)**2

        Cd(4013) = 1/(32.D0*Pi**3)*
     &   (-(Alfas2L*Cd(4009)) + 
     &      64*Pi**3*(SB2*SW2*Cd(50)**2*Cd(3754)*MStop2(1)*
     &         MStop2(2)*MStop2(3)**2))/
     &    (SB2*SW2*Cd(50)**2*MStop2(2)*MStop2(3)**2)

	Cd(4014) = 2*Cd(3489) - Cd(4009)*Cr(34)

	Cd(400) = A0delStop(1)*Cd(50)*Cd(391) - 2*Cd(399)*MStop2(1)

        Cd(401) = A0delStop(2)*Cd(47)*Cd(392)*MStop2(1) + 
     &   Cd(400)*MStop2(2)

	Cd(5994) = M_3C*Cd(5840) - Cd(5993)*UCStopC(1,3)

	Cd(7016) = -(SA2*Cd(5283)) + S2A*Cd(7015)

	Cd(8655) = SA2*Cd(8383) + S2A*Cd(8654)

        Cd(8656) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(-Cd(8655) + 
     &        CA2*Cd(47)**2*Cd(50)**2*Cd(8599)*MStop2(3)))/
     &    (SB2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))

	Cd(8887) = (4*Cd(8883) + 2*MTy*Cd(8886)*MStop2(1))/MTy

	Cd(8888) = 2*SA2*Cd(8702) + S2A*Cd(8781)

        Cd(8889) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(CA2*Cd(8887) + Cd(8888)*MStop2(1)))/
     &    (SB2*SW2*Cd(47)**2*MStop2(1)**2*MStop2(3)**2)

	Cd(9015) = Cd(424)*Cd(8988) + 2*Cd(9014)*MStop2(2)

	Cd(9117) = (2*(2*Cd(9113) + MTy*Cd(9116))*MStop2(2))/MTy

	Cd(9118) = -(S2A*Cd(9015)) + CA2*Cd(9117)

        Cd(9119) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(9118) + 2*SA2*Cd(8934)*MStop2(2)))/
     &    (SB2*SW2*Cd(50)**2*MStop2(2)**2*MStop2(3)**2)

        Cd(9122) = 1/(64.D0*Pi**3)*
     &   (-(Alfas2L*Cd(4917)*Cr(34)) + 
     &      128*Pi**3*(SB2*SW2*Cd(50)**2*Cd(4215)*MStop2(3)**2))/
     &    (SB2*SW2*Cd(50)**2*MStop2(3)**2)

        Cd(9714) = (Cd(50)**2*Cd(9445) + Cd(47)**2*Cd(9713))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12294) = -(Cd(50)**2*Cd(12293)) + 
     &   Cd(47)*(Cd(50)*Cd(12264) + 
     &      Cd(13)*Cd(47)*Cd(12138)*MStop2(3))

	Cd(12651) = SA2*Cd(12607) + CA2*Cd(12624) - S2A*Cd(12650)

        Cd(12695) = 8*MGl2*Cd(12651) - Cd(12693)*MStop2(1) - 
     &   Cd(12694)*MStop2(2)

        Cd(14356) = 4*Cd(14354) - 
     &   (4*Cd(14355) + MTy2*Cd(14344)*MSbot2(1)*MSq2Diff(1,1))/
     &    (MTy2*TB2*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

        Cd(1939) = 2*MUE2*Cd(1842) + 
     &   MUE**2*Cd(1938)*UCStop(1,3) + 
     &   MUEC**2*Cd(1632)*UCStopC(1,3)

        Cd(4015) = 1/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)**2*Cd(4014)*MStop2(1) - 
     &         2*Cd(50)**2*Cd(3112)*MStop2(2)))/
     &     (SB2*SW2*Cd(47)**2*Cd(50)**2*MStop2(1)*MStop2(2)*
     &       MStop2(3)**2) + Cd(3754)*Cr(13)*MStop2(1)

        Cd(8307) = -(Cd(8274)*UCStop(1,3)) + 
     &   M_3*Cd(50)*Cd(8259)*UStop2(1,1)**2 - 
     &   Cd(8306)*UUStopC(1,2)

        Cd(12295) = (Cd(47)*Cd(50)*Cd(12271) + 
     &     MTy*(Cd(47)*Cd(50)*Cd(12283) - Cd(50)*Cd(12290) + 
     &        Cd(12294)))/(MTy*Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(12689) = Cd(12688)*MStop2(3) + 
     &   Cd(12146)*UUStopC(1,1) + Cd(12154)*UUStopC(1,2) + 
     &   Cd(12167)*UUStopC(1,3)

        Cd(2595) = Cd(424)*Cd(1842) + 2*MUEC*Cd(2297) + 
     &   MUE*Cd(2594) + At*MUEC*Cd(1941)*UCStopC(1,3)

        Cd(4012) = 1/128.D0*
     &   (64*A0delStop(1)*Cd(3754) - 
     &     2/Pi**3*(Alfas2L*Cd(4010))/
     &       (SB2*SW2*Cd(47)**2*MStop2(1)*MStop2(3)**2) + 
     &     1/Pi**3*(Alfas2L*
     &         (-(A0delStop(2)*Cd(4009)) + 2*Cd(4011)*MStop2(2)))/
     &       (SB2*SW2*Cd(50)**2*MStop2(2)**2*MStop2(3)**2))

        Cd(6331) = -4*AtC*Cd(6330) + 
     &   At*(4*Cd(5994) - 2*AtC*Cd(6094) + 
     &      At*Cd(5319)*UCStopC(1,3))

        Cd(8005) = 4*AtC*Cd(8004) + 
     &   At*(4*Cd(7753) - 2*AtC*Cd(7842) - 
     &      At*Cd(7272)*UCStopC(1,3))

        Cd(9715) = 2*Cd(9714) + 
     &   AtC*(4*Cd(2593) + AtC*Cd(2298)*UCStop(1,3)) + 
     &   At*(2*AtC*Cd(1842) + 4*Cd(9169) + 
     &      At*Cd(1941)*UCStopC(1,3))

        Cd(12675) = MTy*Cd(12068) + 
     &   Cd(12662)/(MTy2*Cd(47)*Cd(50)) + 
     &   Cd(12664)/(Cd(47)**2*Cd(50)**2) + (2*Cd(12665))/MTy + 
     &   MGl2*Cd(12672) - MTy2*Cd(12674) + 
     &   (MStop2(1)*(Cd(12670) - 4*MGlpTmSt2(1)*MStop2(1)))/
     &    Cd(47)**2 + (4*MTy2**2*Cd(11963) + Cd(12668)*MStop2(2))/
     &    Cd(50)**2 - (4*Cd(11964)*Cd(11965))/
     &    (Cd(47)*Cd(50)*MStop2(3))

	Cd(8308) = 2*Cd(8307) - AtC*Cd(8254)*UCStop(1,3)

	Cd(12296) = SA2*Cd(12069) + CA2*Cd(12134) + S2A*Cd(12295)

        Cd(12297) = 8*MGl2*Cd(12296) + Cd(11834)*MStop2(1) + 
     &   Cd(11959)*MStop2(2)

        Cd(12690) = -(Cd(50)**2*Cd(12689)) + 
     &   Cd(47)*(Cd(50)*Cd(12264) + 
     &      Cd(13)*Cd(47)*Cd(12138)*MStop2(3))

        Cd(6332) = (-2*Cd(5758) + MTy*Cd(6331))/
     &    (MTy*Cd(47)**2*Cd(50)**2*MStop2(3)) - 
     &   AtC**2*Cd(6120)*UCStop(1,3)

        Cd(8006) = (2*Cd(7558) - MTy*Cd(8005))/
     &    (MTy*Cd(47)**2*Cd(50)**2*MStop2(3)) + 
     &   AtC**2*Cd(7863)*UCStop(1,3)

        Cd(12691) = (Cd(47)*Cd(50)*Cd(12681) - 
     &     MTy*(Cd(47)*Cd(50)*Cd(12684) - Cd(50)*Cd(12686) + 
     &        Cd(12690)))/(MTy*Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(8340) = 2*MUEC*Cd(8226) - MUE*Cd(8308) + 
     &   AtC*MUEC*Cd(8339) + 
     &   At*(MUE*Cd(8241) + MUEC*Cd(8051)*UCStopC(1,3))

        Cd(7017) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(Cd(7016) + 
     &        CA2*Cd(47)**2*Cd(50)**2*Cd(6332)*MStop2(3)))/
     &    (SB2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))

	Cd(8341) = -(SA2*Cd(7231)) + S2A*Cd(8340)

        Cd(8342) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(-Cd(8341) + 
     &        CA2*Cd(47)**2*Cd(50)**2*Cd(8006)*MStop2(3)))/
     &    (SB2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))

	Cd(12692) = -(SA2*Cd(12675)) + CA2*Cd(12680) - S2A*Cd(12691)

        Cd(9124) = MGl2*Cd(7017) + Cd(8889)*MStop2(1) - 
     &   Cd(9119)*MStop2(2)

        Cd(12699) = -8*MGl2*Cd(12692) + 8*Cd(12296)*Cr(65) + 
     &   Cd(12697)*MStop2(1) + Cd(12698)*MStop2(2)

        Cd(9120) = A0delGl*Cd(7017) + A0delStop(1)*Cd(8889) - 
     &   A0delStop(2)*Cd(9119) + 2*Cd(8342)*Cr(65)

        Cd(9123) = 1/2.D0*(2*Cd(9122)*MStop2(1) + 
     &      2*Cd(9121)*MStop2(2) + 
     &      4*MGl2*Cd(8656)*MStop2(1)*MStop2(2) + 
     &      Cd(9120)*MStop2(1)*MStop2(2))/(MStop2(1)*MStop2(2))

        Cd(12696) = -(A0delStop(1)*Cd(11834)) - 
     &   A0delStop(2)*Cd(11959) - 8*A0delGl*Cd(12296) - 
     &   2*Cd(12695) + 8*Cd(12692)*Cr(65)

        Cd(9125) = 2*MGl2*Cd(8342) + Cd(7017)*Cr(65) + 
     &   1/(32.D0*Pi**3)*(Alfas2L*
     &       (-(Cd(47)**2*Cd(4917)*MStop2(1)) + 
     &         Cd(50)**2*Cd(4560)*MStop2(2)))/
     &     (SB2*SW2*Cd(47)**2*Cd(50)**2*MStop2(1)*MStop2(2)*
     &       MStop2(3)**2) + Cd(8889)*Cr(13)*MStop2(1) - 
     &   Cd(9119)*Cr(34)*MStop2(2)

	end
