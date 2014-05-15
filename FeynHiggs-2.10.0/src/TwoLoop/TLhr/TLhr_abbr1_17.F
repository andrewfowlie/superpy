* TLhr_abbr1_17.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_17
	implicit none

#include "TLhr.h"


        Cd(14067) = (2*MTy*MTy2*Cd(14065) + 
     &     TB*(Cd(9950)*(-2*MTy*Cd(14049) + 
     &           MTy2*Cd(47)*
     &            (MTy*Cd(14066) + 2*TB2*Cd(14051)*MSq2Diff(1,1))*
     &            MStop2(1)*MStop2(3)) - 
     &        2*MTy*MTy2*MSq2Diff(1,1)*MStop2(1)*
     &         (Cd(14050)*MStop2(1)*MStop2(3) + 
     &           Cd(47)*Cd(9950)*Cd(14047)*UStop2(1,1))*UStop2(1,2)
     &          + 2*MTy*MTy2*Cd(47)*Cd(9950)*Cd(14047)*
     &         MSq2Diff(1,1)*MStop2(3)*UStop2(1,2)**2))/
     &   (MTy*MTy2*TB*TB2*Cd(47)*Cd(9950)*MSq2Diff(1,1)**2*
     &     MStop2(1)*MStop2(3))

        Cd(14262) = (MUEC*TB2*Cd(14239) + MUE*TB2*Cd(14261) + 
     &     (TB*Cd(14258) + 
     &        Cd(14234)*(-(MGl2*Cd(14166)) + TB*Cd(585)*MSbot2(1))*
     &         MSq2Diff(1,1))*MStop2(3) + 
     &     TB*Cd(50)*MSq2Diff(1,2)**2*
     &      (-2*Cd(13896)*Cr(5)*Cr(75) + 
     &        MTy2*TB2*Cd(14246)*MSq2Diff(1,1)*MStop2(3)))/
     &   (TB*TB2*Cd(50)*MSq2Diff(1,1)*MSq2Diff(1,2)**2*MStop2(3))

        Cd(147) = -(MTy2*Cd(143)) - MTy*Cd(146) + 
     &   M_3C*Cd(139)*UCStop(1,3) + 
     &   M_3*(3*MTy2**2 + (MGl2 - MStop2(2))**2)*UCStopC(1,3)

        Cd(8259) = Cd(44)*Cd(6732) + 2*MTy2*Cd(8258) + 
     &   MGl2*(MGl2**2*(-(MGl2*Cd(1647)) + Cd(6738)) + 
     &      MTy2*Cd(6741) - MGl2*Cd(6743) + 
     &      Cd(874)*MStop2(1)**2*MStop2(3))

        Cd(10159) = 2*Cd(9954) + Cd(9948)*Cd(10158) + 
     &   Cd(10102)*UStop2(1,1) - 
     &   Cd(9952)*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,1)**2 + 3*UStop2(1,2)**2))

        Cd(10784) = Cd(50)*Cd(10309) - MTy2*Cd(10783) + 
     &   (-2*MTy2**2*Cd(251)*Cd(10313) - MGl2*Cd(585)*Cd(10778) + 
     &      MGlmT2*Cd(50)*Cd(10779))/MStop2(2) + 
     &   (Cd(50)*Cd(10630))/MStop2(3) + 
     &   AtC*Cd(10777)*UCStop(1,3) + 
     &   At*Cd(179)*Cd(251)*UCStopC(1,3)

        Cd(10801) = (-8*Cd(50)*
     &      (MGl2*Cd(10624) + Cd(10608)*MStop2(1)*MStop2(2)) + 
     &     (-32*MTy2**3*Cd(251)*Cd(10193) - 8*Cd(50)*Cd(10774) + 
     &        (8*MTy2**2*Cd(251)*Cd(10776) - 16*MTy*Cd(10784) + 
     &           Cd(50)*(Cd(10773) + MTy2*Cd(50)*Cd(10800)))*
     &         MStop2(2))*MStop2(3))/
     &   (Cd(50)**2*MStop2(2)*MStop2(3))

        Cd(12359) = At*Cd(10388)*Cd(12357) + 
     &   4*Cd(4187)*MStop2(2)*
     &    (M_3C*MUE*UCStop(1,3)**2 + M_3*MUEC*UCStopC(1,3)**2)+
     &     Cd(9948)*Cd(12358)*UStop2(1,1) + 
     &   AtC*(Cd(50)**2 + Cd(251)*Cd(2681)*MStop2(2))*
     &    (MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)**2)

        Cd(13849) = -(MTy2*Cd(13846)*Cd(13848)) + 
     &   MGl2*Cd(13846)*MSq2Diff(1,2)*(MGl2 - 2*MStop2(2)) + 
     &   MSq2Diff(1,2)*(M_3C*Cd(13709) + M_3*Cd(13718) + 
     &      MStop2(2)**2*(Cd(585) + MUE*TB*UCStop(1,3) + 
     &         MUEC*TB*UCStopC(1,3)))

        Cd(13869) = MUEC*TB*Cd(13856) - 
     &   MGl2**2*Cd(10367)*Cd(13858) - MTy2*Cd(13865) - 
     &   MTy*Cd(13868) - MSq2Diff(1,1)*MStop2(1)*
     &    (-(MGl2*Cd(13860)) + 
     &      (At*MUE*TB*Cd(10367) + Cd(13859))*MStop2(1)) + 
     &   MTy2**2*Cd(13857)*UStop2(1,1)

        Cd(6279) = MStop2(3)*
     &    (MGl2*(MGl2*(-(Cd(843)*Cd(6274)) + 
     &            MGl2*(MGl2*
     &                (2*MGl2*Cd(6268) - Cd(843)*Cd(6271)) + 
     &               4*Cd(6276))) + 2*Cd(6278)*MStop2(1)*MStop2(2))
     &        + Cd(1204)*MStop2(1)**2*MStop2(2)**2*MStop2(4)) + 
     &   2*Cd(5776)*UStop2(1,1)**3

        Cd(7967) = MStop2(3)*
     &    (MGl2*(MGl2*(Cd(843)*Cd(6274) + 
     &            MGl2*(MGl2*
     &                (-2*MGl2*Cd(6268) + Cd(843)*Cd(6271)) - 
     &               4*Cd(6276))) - 2*Cd(6278)*MStop2(1)*MStop2(2))
     &        - Cd(1204)*MStop2(1)**2*MStop2(2)**2*MStop2(4)) - 
     &   2*Cd(7577)*UStop2(1,1)**3

        Cd(10366) = (MGl2*Cd(10302) + 
     &     MStop2(1)*(Cd(10293) + 
     &        4*MGlmT2*Cd(47)*Cd(10294)*MStop2(1)) + 
     &     (2*Cd(47)**2*Cd(10191) + 4*MTy2**2*Cd(10281) + 
     &        MStop2(1)**2*
     &         (MTy2*Cd(10365) + 
     &           Cd(47)*(-4*Cd(10282) + MTy*Cd(47)*Cd(10336) - 
     &              4*Cd(10295)*MStop2(1))))*MStop2(3))/
     &   (Cd(47)**2*MStop2(1)**2*MStop2(3))

        Cd(13756) = -2*MTy2*Cd(141)*Cd(13753) + 
     &   MGl2**2*Cd(13754) + 
     &   2*(M_3C*Cd(13709) + M_3*Cd(13718))*MSq2Diff(1,2) + 
     &   MStop2(2)*(-2*MGl2*Cd(13755) + 
     &      2*MSq2Diff(1,2)*MStop2(2)*
     &       (Cd(585) + MUE*TB*UCStop(1,3)) + 
     &      MUEC*TB*(MTy2**2 + Cd(13719)*MStop2(2))*UCStopC(1,3))

        Cd(2448) = MTy2*Cd(2435) + 
     &   MStop2(3)*(-Cd(2428) + 
     &      M_3C*MGl2*(MGl2*(MGl2 - 3*MStop2(1)) + 
     &         3*MStop2(1)**2)*MStop2(3)*UCStop(1,3)**2 + 
     &      MTy2**2*(Cd(2447) + 
     &         MTy2*(-Cd(2439) + 2*MTy2*Cd(2427)*UStop2(1,1)))) + 
     &   MTy*Cd(2426)*UCStop(1,3)*UStop2(1,2)

        Cd(2511) = MTy2*Cd(2504) + 
     &   MStop2(3)*(Cd(2493) + 
     &      MTy2**3*(-2*M_3*MTy2*Cd(2364) + Cd(2495)) - 
     &      MTy2**2*Cd(2510) + 
     &      M_3C*MGl2*(MGl2*(MGl2 - 3*MStop2(2)) + 
     &         3*MStop2(2)**2)*MStop2(3)*UCStop(1,3)**2) + 
     &   MTy*Cd(2492)*UCStop(1,3)*UStop2(1,2)

        Cd(2811) = -(MTy*Cd(2810)) + 
     &   2*(MTy2**2*(-(MTy2*Cd(2804)) + Cd(2805)) + 
     &      (-3*MTy2*Cd(118)*Cd(121)*Cd(2806) + 
     &         Cd(2802)*MStop2(1)**4 + 
     &         MGl2*(-10*Cd(843)*MStop2(1)**3 + 
     &            MGl2*(6*Cd(595)*MStop2(1)**2 + 
     &               MGl2*(MGl2*Cd(594) - 2*Cd(2803)*MStop2(1)))))*
     &       MStop2(3))*UStop(1,1)

        Cd(4266) = MTy*Cd(4265) + 
     &   2*(MTy2*Cd(4258) + 
     &      (MTy2**2*(3*Cd(4259) + 
     &            MTy2*(MTy2*Cd(594) - Cd(4261))) + 
     &         MGl2**2*(-5*MGl2 + 3*Cd(594))*MStop2(1)**2 - 
     &         3*MGl2*Cd(595)*MStop2(1)**3 + Cd(2802)*MStop2(1)**4)
     &        *MStop2(3))*UStop(1,1)

        Cd(4435) = M_3*MTy2*Cd(4420) + 
     &   MStop2(3)*(Cd(4434) + 
     &      2*MTy*Cd(4433)*UCStop(1,3)*UStop2(1,2) + 
     &      2*M_3*MStop2(1)*
     &       (-(MGl2*MStop2(1)*
     &            (MGl2*(MGl2 - 3*MStop2(1)) + 3*MStop2(1)**2)) + 
     &         MTy2**2*(Cd(4421) + 
     &            MTy2*(-Cd(4422) + MTy2*UStop2(1,2)**2))))

	Cd(4772) = -(Cd(4771)*MStop2(3)) + Opt(8)

	Cd(4913) = -(Cd(4912)*MStop2(3)) + Opt(8)

        Cd(10100) = (MGlpT2*Cd(10084))/(Cd(47)**2*Cd(50)**2) + 
     &   2*MTy*Cd(10099) + Cd(10072)*MStop2(1) + 
     &   MUE**2*Cd(10075)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(10079)*UCStopC(1,3)**2 + 
     &   2*MUE2*((Cd(9949)*MGlpTmSt2(2))/Cd(50) - 
     &      (Cd(10014)*Cr(28))/MStop2(3) + 
     &      (UStop2(1,1)*(-(Cd(47)*Cd(9952)*Cd(9981)*MStop2(2)) + 
     &           2*Cd(50)**2*Cd(9947)*MGlpTmSt2(1)*UStop2(1,2)))/
     &       (Cd(47)*Cd(50)**2))

        Cd(11893) = (32*MTy*MTy2**3*Cd(11836) - 
     &     8*MTy*MTy2*Cd(50)*Cd(11846)*MStop2(2) + 
     &     (16*MTy2*Cd(11843) - 8*MTy*MTy2*Cd(50)*Cd(11857) + 
     &        (-8*MTy*Cd(50)*Cd(11835) + 
     &           MTy2*(16*MTy2*Cd(11873) + 
     &              MTy*(8*MGl2*
     &                  (MGl2*Cd(179)*Cd(11849) - Cd(11865)) + 
     &                 Cd(50)*(Cd(11862) + MTy2*Cd(50)*Cd(11892))))
     &           )*MStop2(2))*MStop2(3))/
     &   (MTy*MTy2*Cd(50)**2*MStop2(2)*MStop2(3))

        Cd(13850) = MTy2*
     &    (MTy2*Cd(13737)*Cd(13840) + Cd(13845)) + MTy*Cd(13849) + 
     &   MSq2Diff(1,2)*MStop2(2)*
     &    (MGl2*Cd(13838) + Cd(13839)*MStop2(2)) + 
     &   (MUEC*TB*Cd(13836) + MGl2**2*Cd(13837) + 
     &      At*(AtC*Cd(13835) - 
     &         MUE*TB*Cr(56)*MSq2Diff(1,2)*MStop2(2)**2))*
     &    UStop2(1,1) + Cd(13831)*UStop2(1,1)**2

        Cd(14002) = MTy2**3*Cd(13734)*Cd(13846) + 
     &   MTy2**2*Cd(13994) + 
     &   2*(M_3C*Cd(13932) + M_3*Cd(13939))*MSq2Diff(1,2) - 
     &   MStop2(2)*(MTy2*Cd(14001) + 
     &      MSq2Diff(1,2)*(MGl2*Cd(13801)*(MGl2 - 2*MStop2(2)) + 
     &         MStop2(2)**2*(Cd(585) + MUE*TB*UCStop(1,3)))) - 
     &   MUEC*TB*Cd(13970)*UStopC(1,1)

        Cd(14204) = Cd(9948)*Cd(14203) + 
     &   UStop2(1,1)*(-(At*Cd(14125)*Cd(14201)) + 
     &      MSq2Diff(1,2)*(-(AtC*MUEC*TB2*MGlpTmSt2(2)) - 
     &         MUE2*TB*TB2*MGlpTmSt2(2) + TB*Cd(251)*MSq2Diff(1,2))
     &        + (AtC*MUEC*TB2*Cd(50) + MUE2*TB*TB2*Cd(50) + 
     &         TB*MGlpTmSt2(2)*MSq2Diff(1,2)**2)*UStop2(1,1)) + 
     &   2*TB*Cd(50)*MSq2Diff(1,2)*UStop2(1,2)

        Cd(3032) = MTy*Cd(3031) + 
     &   (-6*MTy2**4*Cd(3013) + MTy2**2*Cd(3017) + 
     &      MStop2(3)*(MTy2*Cd(118)**2*Cd(3020) + 
     &         MTy2**3*Cd(3023) - 
     &         (MGl2 - MStop2(1))**4*MStop2(1)*MStop2(3)))*
     &    UUStop(1,1)

        Cd(4519) = MTy*Cd(4518) + 
     &   (MTy2**2*Cd(4501) + MTy2**3*Cd(4506) + 
     &      MStop2(3)*(6*MTy2**5*Cd(607) - 3*MTy2**4*Cd(4503) - 
     &         MTy2*Cd(118)*Cd(4497)*MStop2(1) + 
     &         (MGl2 - MStop2(1))**3*MStop2(1)**2*MStop2(3)))*
     &    UUStop(1,1)

        Cd(4611) = 2*MTy2*Cd(4606)*UStop(1,1) + 
     &   MStop2(3)*(-4*MTy*Cd(4610)*MStop2(2) + 
     &      2*(MTy2**2*(MTy2*(MTy2*Cd(874) - Cd(4605)) + 
     &            3*Cd(4607)) + 
     &         MGl2**2*(-4*MGl2 + 3*Cd(874))*MStop2(2)**2 - 
     &         6*MGl2*Cd(843)*MStop2(2)**3 + Cd(607)*MStop2(2)**4)*
     &       UStop(1,1) + Cd(4604)*UStop2(1,2))

        Cd(4851) = MTy*Cd(4850) - MTy2*Cd(4840)*MStop2(3) + 
     &   (MTy2**2*Cd(4831) - 
     &      MTy2**3*(MTy2*(6*MTy2*Cd(1272) - 3*Cd(4833)) + 
     &         Cd(4836))*MStop2(3) + 
     &      (MGl2 - MStop2(2))**3*MStop2(2)**2*MStop2(3)**2)*
     &    UUStop(1,1)

        Cd(4900) = -2*MTy2*Cd(4606)*UStop(1,1) + 
     &   MStop2(3)*(4*MTy*Cd(4610)*MStop2(2) - 
     &      2*(MTy2**2*(MTy2*(MTy2*Cd(874) - Cd(4605)) + 
     &            3*Cd(4607)) + 
     &         MGl2**2*(-4*MGl2 + 3*Cd(874))*MStop2(2)**2 - 
     &         6*MGl2*Cd(1272)*MStop2(2)**3 + Cd(607)*MStop2(2)**4)
     &        *UStop(1,1) - Cd(4604)*UStop2(1,2))

        Cd(8986) = -2*(Cd(8985) + 
     &     MStop2(2)*(4*MTy*Cd(4610)*MStop2(3) + 
     &        (MTy2*(MTy2*(MTy2*(-(MTy2*Cd(3781)) + 2*Cd(8908)) - 
     &                 2*Cd(8912)) + 2*Cd(8914)) - 
     &           MGl2*MStop2(2)*
     &            (MGl2*(MGl2*(MGl2 - 4*Cd(652)) + 2*Cd(3790)) - 
     &              4*Cd(3785)*MStop2(2)))*UStop(1,1)))

        Cd(13758) = -(MTy2**2*Cd(13740)) - MTy2*Cd(13751) + 
     &   MTy*Cd(13756) + MUEC*TB*Cd(13757) - 
     &   2*MGl2*Cd(13733)*MStop2(2) + 
     &   MTy2**3*Cd(13736)*MStop2(2) + Cd(13735)*MStop2(2)**2 + 
     &   (-(MGl2**2*Cd(13724)) + 
     &      At*(AtC*Cd(13720) - MUE*TB*Cd(13721)*MStop2(2)**2))*
     &    UStop2(1,1) + Cd(13705)*UStop2(1,1)**2

        Cd(3520) = 4*MTy*Cd(3519)*MStop2(3) + 
     &   (-(MGl2**4*Cd(3504)) + 2*MGl2**3*Cd(3505) + 
     &      2*MTy2*(Cd(118)*Cd(3508) + 
     &         MTy2*(MTy2*Cd(121)*Cd(1160) - Cd(3513))) - 
     &      MStop2(1)*(MTy2**4 + Cd(1803)*MStop2(1)**2 + 
     &         MGl2*(6*MGl2*Cd(3515) - 2*Cd(3517)*MStop2(1))))*
     &    UStop(1,1) + Cd(3503)*UStop2(1,2)

        Cd(3710) = 2*MTy2*Cd(3709) + 
     &   MStop2(1)*(2*MTy*Cd(3679) + 
     &      Cd(3681)*(MGl2 - MStop2(1))**4*MStop2(3) + 
     &      MTy2**2*(2*Cd(3701) + 
     &         MTy2*(2*Cd(3691) + 
     &            MTy2*(Cd(3684) - 12*MTy2*MStop2(1)*UUStop(1,1))))
     &      )

        Cd(3792) = 4*MTy*Cd(3171)*MStop2(3) + 
     &   (MGl2**4*Cd(3781) - 2*MGl2**3*Cd(3782) + 
     &      2*MTy2*(-2*Cd(138)*Cd(3786) + 
     &         MTy2*(-2*MTy2*Cd(141)*Cd(652) + Cd(3791))) + 
     &      MStop2(2)*(MTy2**4 + Cd(1756)*MStop2(2)**2 + 
     &         MGl2*(6*MGl2*Cd(3783) - 2*Cd(3787)*MStop2(2))))*
     &    UStop(1,1) - Cd(3780)*UStop2(1,2)

        Cd(3965) = 2*MTy*Cd(3936) - 
     &   MTy2*(MTy2*(MTy2*(MTy2*Cd(3941) + 2*Cd(3948)) + 
     &         2*Cd(3958)) + 2*Cd(3964)) - 
     &   Cd(3937)*(MGl2 - MStop2(2))**4*MStop2(3) + 
     &   12*MTy2**5*MStop2(2)*UUStop(1,1)

        Cd(6965) = Cd(44)*Cd(48)*Cd(6933) - MGl2*Cd(6942) - 
     &   MTy2*Cd(6964) + MGl2**2*
     &    (Cd(6937) + MGl2*
     &       (MGl2**3*Cd(579) - 
     &         MTy2**2*(MTy2*(80*MTy2 - 4*Cd(6958)) + Cd(6959)) + 
     &         MStop2(3)*(MGl2*Cd(6534) - 
     &            (3*MGl2**2 + Cd(5921))*MStop2(4))))

        Cd(8684) = -4*MTy*Cd(8683)*MStop2(3) + 
     &   (-(MTy2**2*(MTy2*(MTy2*Cd(3504) + 2*Cd(8673)) + 
     &           2*Cd(8677))) + 2*MTy2*Cd(8680) - 
     &      MStop2(1)*(Cd(1803)*MStop2(1)**2 + 
     &         MGl2*(MGl2*
     &             (MGl2*(MGl2 - 2*Cd(1160)) + 2*Cd(8674)) - 
     &            2*Cd(8672)*MStop2(1))))*UStop(1,1) + 
     &   Cd(8671)*UStop2(1,2)

        Cd(8841) = 2*MTy2*Cd(8840) + 
     &   MStop2(1)*(-2*MTy*Cd(8808) + 
     &      MTy2**2*(2*Cd(8831) + 
     &         MTy2*(2*Cd(8821) + 
     &            MTy2*(Cd(8814) - 12*MTy2*Cr(2)*UUStop(1,1)))) + 
     &      MStop2(1)*MStop2(3)*
     &       (Cd(3681)*MStop2(1)**3 + 
     &         MGl2*(-2*Cd(8816)*MStop2(1)**2 + 
     &            MGl2*(6*Cd(8811)*MStop2(1) + 
     &               MGl2*(-2*Cd(8810) + MGl2*UUStop(2,1))))))

        Cd(8915) = 4*MTy*Cd(4610)*MStop2(3) + 
     &   (MTy2*(MTy2*(MTy2*(-(MTy2*Cd(3781)) + 2*Cd(8908)) - 
     &            2*Cd(8912)) + 2*Cd(8914)) - 
     &      MStop2(2)*(Cd(1756)*MStop2(2)**2 + 
     &         MGl2*(MGl2*(MGl2*(MGl2 - 4*Cd(652)) + 2*Cd(8909)) - 
     &            4*Cd(3785)*MStop2(2))))*UStop(1,1) + 
     &   Cd(8907)*UStop2(1,2)

        Cd(9074) = 2*MTy*Cd(9047) + 2*MTy2*Cd(9073) + 
     &   MTy2**2*(2*Cd(9066) - 
     &      MTy2*(-2*Cd(9057) + 
     &         MTy2*(Cd(9051) + 12*MTy2*Cd(652)*UUStop(1,1)))) + 
     &   MStop2(2)*MStop2(3)*
     &    (Cd(3937)*MStop2(2)**3 + 
     &      MGl2*(-2*Cd(9052)*MStop2(2)**2 + 
     &         MGl2*(6*Cd(9049)*MStop2(2) + 
     &            MGl2*(-2*Cd(9048) + MGl2*UUStop(2,1)))))

        Cd(9484) = -(MTy2*Cd(9478)) + MTy*Cd(9483) + 
     &   (2*MTy2**3*Cd(9475) + 
     &      MStop2(3)*(MTy2**4*(-12*MTy2 + Cd(9472)) + 
     &         16*MTy2**2*Cd(9476) - 
     &         (MGl2 - MStop2(2))**4*MStop2(3)))*UUStop(1,1)

        Cd(12524) = 4*Cd(12482) - 
     &   (MGlmT2*Cd(50)*Cd(12523))/MStop2(3) - 
     &   UStop2(1,1)**2*(-5*MGlmT2*Cd(50) + 3*Cd(9948)*Cd(12522) + 
     &      (4*MGlmT2*Cd(50) + Cd(12483))*UStop2(1,1)) + 
     &   2*MGlmT2*Cd(50)*UStop2(1,1)*(5 + 6*UStop2(1,1))*
     &    UStop2(1,2) + (5*MGlmT2*Cd(50) - 
     &      3*(4*MGlmT2*Cd(50) + Cd(12483))*UStop2(1,1))*
     &    UStop2(1,2)**2 - 
     &   (4*MGlmT2*Cd(50) + Cd(12483))*UStop2(1,2)**3

        Cd(1772) = -4*MTy*Cd(1771) - 
     &   (MTy2*(MTy2*(2*Cd(1754) + 
     &            MTy2*(MTy2*(MTy2 - Cd(1748)) + 2*Cd(1762))) - 
     &         Cd(138)*Cd(1766)) - Cd(1756)*MStop2(2)**3 + 
     &      MGl2*(Cd(1760)*MStop2(2)**2 + 
     &         MGl2*(MGl2*(MGl2*(MGl2 - Cd(1747)) + 2*Cd(1750)) - 
     &            2*Cd(1758)*MStop2(2))))*UStop(1,1) + 
     &   Cd(1746)*UStop2(1,2)

        Cd(1818) = -4*MTy*Cd(1817) + 
     &   (MTy2*(MTy2*(MTy2*
     &             (-2*Cd(1796) + MTy2*(-MTy2 + Cd(1797))) + 
     &            2*Cd(1808)) + Cd(118)*Cd(1812)) - 
     &      MGl2**3*(2*Cd(1801) + MGl2*(MGl2 - 5*Cr(1))) + 
     &      2*MGl2**2*Cd(1805)*MStop2(1) - 
     &      MGl2*Cd(1799)*MStop2(1)**2 + Cd(1803)*MStop2(1)**3)*
     &    UStop(1,1) + Cd(1794)*UStop2(1,2)

        Cd(7627) = MGl2**2*
     &    (Cd(7615) + MGl2*
     &       (MGl2*(MGl2*(MGl2 - 3*Cd(198)) + Cd(7612)) - Cd(7618))
     &      ) - MGl2*Cd(7625) - MTy2*Cd(7626) + 
     &   Cd(1280)*MStop2(1)**2*MStop2(2)**2 + 
     &   MTy2**2*(Cd(7611) + 
     &      MTy2*(-2*Cd(1272)*Cd(7608) + 
     &         MTy2*(15*Cd(286)*Cd(1898) + 
     &            MTy2*(13*MTy2 - 33*MStop2(4)))))

        Cd(9362) = MTy*Cd(9361) + 
     &   MTy2*(Cd(9323) + MTy2*
     &       (Cd(9307) + MTy2*
     &          (-Cd(9314) + 
     &            MTy2*(Cd(9298) + 
     &               MTy2*(Cd(9286) - 12*MTy2*UUStop(1,1)))))) + 
     &   MStop2(3)*(-(Cd(3681)*MStop2(1)**4) + 
     &      MGl2*(Cd(9288)*MStop2(1)**3 + 
     &         MGl2*(-2*Cd(9290)*MStop2(1)**2 + 
     &            MGl2*(2*Cd(9292)*MStop2(1) + 
     &               MGl2*(Cd(9284) + MGl2*UUStop(2,1))))))

        Cd(11925) = At*Cd(11922) + Cd(9948)*Cd(11924) + 
     &   2*MUE*(AtC*Cd(50)*(MGlmT2 + MStop2(2)) + 
     &      Cd(179)*(2*M_3C*Cd(179) + AtC*MStop2(2))*MStop2(3))*
     &    UCStop(1,3)**2 + 
     &   MUEC*(4*M_3*Cd(179)**2*MStop2(3)*UCStopC(1,3)**2 + 
     &      AtC*UStop2(1,1)*
     &       (-(Cd(50)*(MGlmT2 + MStop2(2))*
     &            (UStop2(1,1) - UStop2(1,2))**2) + 
     &         2*Cd(179)*MStop2(2)*MStop2(3)*UStop2(1,2)**2))

        Cd(8373) = -(MTy*Cd(8372)) + Opt(21) + 
     &   (-(MTy2**4*(2*MTy2*Cd(8359) + Cd(8362))) + 
     &      4*MTy2**3*Cd(8367) + 2*MTy2*Cd(8370) + 
     &      MGl2**4*(-6*MGl2*Cd(843) + Cd(8360))*MStop2(3))*
     &    UStop(1,1)

        Cd(8405) = -Opt(21) + 
     &   (2*MTy2*Cd(8403) + 
     &      MTy2**3*(4*Cd(8402) + 
     &         MTy2*(2*MTy2*Cd(8399) + Cd(8404))) + 
     &      MGl2**4*(6*MGl2*Cd(1272) - Cd(8401))*MStop2(3))*
     &    UStop(1,1) - 2*M_3*MTy*Cd(8238)*MStop2(3)**2*
     &    UStop(1,2)*UStop2(1,2)

        Cd(12516) = Cd(12482)*Cd(12515) + 
     &   (6*M_3C*MGlmT2*Cd(50) - AtC*Cd(12514))*UCStop(1,3) + 
     &   UCStopC(1,3)*(6*M_3*MGlmT2*Cd(50) + 
     &      At*(Cd(12483) + 
     &         2*MGlmT2*Cd(50)*
     &          (UStop2(1,1)*(-1 - 2*Cd(9948) + 2*UStop2(1,1)) + 
     &            (-1 + 2*Cd(9948) + 6*UStop2(1,1))*UStop2(1,2))))

        Cd(2397) = MStop2(3)*
     &    (MGl2**4*Cd(2355) - MGl2**3*Cd(2358) + 
     &      MStop2(1)*(MGl2*
     &          (3*MGl2*Cd(2357) - Cd(2356)*MStop2(1)) + 
     &         MStop2(1)**2*
     &          (M_3*MStop2(1) - M_3C*MStop2(4)*UCStop(1,3)**2)
     &         )) - UStop2(1,1)*
     &    (Cd(2342) + Cd(2319)*UStop2(1,1)) + 
     &   MTy2*(Cd(2376) + MTy2*
     &       (-Cd(2396) + MTy2*
     &          (-Cd(2383) + 
     &            MTy2*(Cd(2363) + 
     &               8*M_3C*MTy2*UCStop(1,3)**2*UStop2(1,1)))))-
     &     2*MTy*Cd(2354)*UCStop(1,3)*UStop2(1,2)

        Cd(2591) = MStop2(3)*
     &    (MGl2**4*Cd(2355) - MGl2**3*Cd(2556) + 
     &      MStop2(2)*(MGl2*
     &          (3*MGl2*Cd(2555) - Cd(2554)*MStop2(2)) + 
     &         MStop2(2)**2*
     &          (M_3*MStop2(2) - M_3C*MStop2(4)*UCStop(1,3)**2)
     &         )) - UStop2(1,1)*
     &    (Cd(2544) + Cd(2531)*UStop2(1,1)) + 
     &   MTy2*(Cd(2579) + MTy2*
     &       (Cd(2590) + MTy2*
     &          (-Cd(2565) + 
     &            MTy2*(Cd(2559) + 
     &               8*M_3C*MTy2*UCStop(1,3)**2*UStop2(1,1)))))-
     &     2*MTy*Cd(2553)*UCStop(1,3)*UStop2(1,2)

        Cd(8419) = MTy2**2*Cd(8361) + 
     &   MGl2*(-2*Cd(44)*Cd(48)*Cd(8416) + 
     &      MGl2*(Cd(8417) + 
     &         MGl2*(MGl2*(2*MGl2*(MGl2 - Cd(5954)) + Cd(8415)) - 
     &            4*Cd(8418)))) + 
     &   Cd(5062)*MStop2(1)**2*MStop2(2)**2 - 
     &   2*MTy2*Cd(1272)*(2*MTy2**2*Cd(5444) + 
     &      Cd(5397)*MStop2(1)*MStop2(2)) + 
     &   MTy2**4*(Cd(8360) + 2*MTy2*(MTy2 - 3*MStop2(4))) - 
     &   2*Cd(8414)*UStop2(1,2)**2 + Cd(7300)*UStop2(1,2)**3

        Cd(7696) = 4*MGl2**7 + 6*MGl2**5*Cd(5643) + 
     &   MGl2**3*Cd(7669) - 
     &   MTy2*(Cd(7685) + MTy2*
     &       (MTy2*(2*Cd(7674) + 
     &            MTy2*(2*MTy2*
     &                (MTy2*(-2*MTy2 + Cd(7670)) - 6*Cd(7671)) + 
     &               Cd(7688))) - Cd(7695))) + 
     &   2*MGl2*Cd(7667)*MStop2(1)*MStop2(2) - 
     &   Cd(1272)*(MGl2**2*(MGl2**2*Cd(7663) + Cd(7665)) + 
     &      Cd(7686)*MStop2(1)**2*MStop2(2)**2) - 
     &   14*MGl2**6*MStop2(4) - Cd(302)*Cd(7661)*UStop2(1,2)

        Cd(12588) = -(Cd(47)**2*Cd(251)*Cd(10141)) + 
     &   (AtC*Cd(12587) + 2*M_3C*Cd(50)**2*
     &       (2*Cd(164) + MGlpTmSt2(1)) - 
     &      2*M_3C*Cd(47)**2*MGlpTmSt2(2))*UCStop(1,3) + 
     &   (-(Cd(50)**2*(At*(Cd(164) + MGlpTmSt2(1)) - 
     &           2*M_3*(2*Cd(164) + MGlpTmSt2(1)))) + 
     &      (-2*M_3 + At)*Cd(47)**2*MGlpTmSt2(2))*UCStopC(1,3)

        Cd(5219) = MTy*Cd(5218) + 
     &   2*(-(MTy2*(MTy2*(MTy2*
     &               (-Cd(5187) + 
     &                 MTy2*
     &                  (MTy2*
     &                     (MTy2*(7*MTy2 - Cd(5169)) + Cd(5188)) - 
     &                    Cd(5194))) + Cd(5201)) + Cd(5209))) + 
     &      MGl2**2*(Cd(5181) + 
     &         MGl2*(-Cd(5175) + 
     &            MGl2*(MGl2*
     &                (MGl2*(-7*MGl2 + Cd(5168)) - 3*Cd(5183)) + 
     &               5*Cd(5173)*MStop2(1)))) - 
     &      MGl2*Cd(5178)*MStop2(1)*MStop2(2) + 
     &      Cd(5171)*MStop2(1)**2*MStop2(2)**2)*UStop(1,1) + 
     &   Cd(47)*Cd(5167)*UStop2(1,2)

        Cd(7145) = 2*MTy*Cd(7144) + 
     &   (3*MGl2**2*Cd(7113) - MGl2**3*Cd(7117) + 
     &      MGl2**4*(MGl2*(MGl2*(-7*MGl2 + Cd(7103)) - 
     &            3*Cd(7106)) + Cd(7120)) + 
     &      MTy2*(MTy2*(MTy2*
     &             (Cd(7125) + 
     &               MTy2*(MTy2*
     &                   (MTy2*(-7*MTy2 + Cd(7104)) - 3*Cd(7107))+
     &                    Cd(7129))) + Cd(7136)) + Cd(7143)) + 
     &      MGl2*Cd(7110)*MStop2(1)*MStop2(2) + 
     &      Cd(7114)*MStop2(1)**2*MStop2(2)**2)*UStop(1,1) + 
     &   Cd(47)*Cd(7102)*UStop2(1,2)

        Cd(10127) = 2*Cd(47)*Cd(50)*Cd(9948)*Cd(10123) + 
     &   2*Cd(47)**2*Cd(50)*Cd(9949)*Cd(10124) - 
     &   (2*Cd(47)*Cd(50)*Cd(10126))/MStop2(3) - 
     &   At**2*(Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952))*
     &    UCStopC(1,3)**2 + 
     &   2*Cd(47)*Cd(50)**2*Cd(164)*Cd(9947)*UStop2(1,1) - 
     &   AtC**2*Cd(50)**2*Cd(9951)*UCStop(1,3)**2*UStop2(1,1) - 
     &   AtC**2*Cd(47)**2*Cd(9952)*UCStop(1,3)**2*UStop2(1,1) + 
     &   2*Cd(47)*Cd(50)**2*Cd(9947)*MGlpTmSt2(1)*UStop2(1,1)**2 - 
     &   AtC**2*Cd(50)**2*Cd(9951)*UCStop(1,3)**2*UStop2(1,2) - 
     &   AtC**2*Cd(47)**2*Cd(9952)*UCStop(1,3)**2*UStop2(1,2) - 
     &   2*At*AtC*Cd(50)**2*Cd(9951)*UStop2(1,1)**2*UStop2(1,2) - 
     &   2*At*AtC*Cd(47)**2*Cd(9952)*UStop2(1,1)**2*UStop2(1,2) + 
     &   2*Cd(47)*Cd(50)**2*Cd(9947)*MGlpTmSt2(1)*UStop2(1,2)**2 - 
     &   2*At*AtC*Cd(50)**2*Cd(9951)*UStop2(1,1)*UStop2(1,2)**2 - 
     &   2*At*AtC*Cd(47)**2*Cd(9952)*UStop2(1,1)*UStop2(1,2)**2

        Cd(12496) = Cd(9948)*Cd(12495) - 
     &   2*M_3C*AtC*Cd(50)*UCStop(1,3)**2 - 
     &   2*M_3*At*Cd(50)*UCStopC(1,3)**2 + 
     &   Cd(12481)*UStop2(1,1) - MGlmT2*Cd(50)*UStop2(1,1)**2 - 
     &   4*Cd(12482)*UStop2(1,1)**2 - 
     &   At*AtC*Cd(50)*UStop2(1,1)**3 - 
     &   AtC**2*Cd(50)*UCStop(1,3)**2*UStop2(1,2) - 
     &   2*M_3C*At*Cd(50)*UStop2(1,1)*UStop2(1,2) - 
     &   2*M_3*AtC*Cd(50)*UStop2(1,1)*UStop2(1,2) - 
     &   2*MGlmT2*Cd(50)*UStop2(1,1)*UStop2(1,2) - 
     &   8*Cd(12482)*UStop2(1,1)*UStop2(1,2) - 
     &   3*At*AtC*Cd(50)*UStop2(1,1)**2*UStop2(1,2) - 
     &   MGlmT2*Cd(50)*UStop2(1,2)**2 - 
     &   4*Cd(12482)*UStop2(1,2)**2 - 
     &   3*At*AtC*Cd(50)*UStop2(1,1)*UStop2(1,2)**2 - 
     &   At*AtC*Cd(50)*UStop2(1,2)**3 - 
     &   Cd(50)*Cd(138)*(4*UStop2(1,1)**2 + 4*UStop2(1,2)**2 + 
     &      UStop2(1,1)*(-5 + 8*UStop2(1,2)))

        Cd(10145) = Cd(47)**2*Cd(9954)*Cd(10141) + 
     &   (-2*Cd(47)*Cd(50)*Cd(9955)*Cd(10144) + 
     &      MStop2(3)*((-2*M_3C*
     &             (Cd(50)**2*
     &                (-Cd(9951) + 2*Cd(9953) + Cd(10101)) + 
     &               Cd(47)**2*(Cd(9952) - Cd(10102))) + 
     &            AtC*Cd(10139))*UCStop(1,3) + 
     &         (At*(Cd(47)**2*Cd(9952) + 
     &               Cd(50)**2*(-Cd(9951) + Cd(9953))) - 
     &            2*M_3*(Cd(50)**2*
     &                (-Cd(9951) + 2*Cd(9953) + Cd(10101)) + 
     &               Cd(47)**2*(Cd(9952) - Cd(10102))))*
     &          UCStopC(1,3)))/MStop2(3)

        Cd(13921) = Cd(9948)*Cd(13911) - 2*At*Cd(13920) + 
     &   MGl2*(MSq2Diff(1,1)*UStop2(1,1)*
     &       (TB*(Cd(47) + Cd(50))*MSq2Diff(1,2) + 
     &         2*Cd(47)*(AtC*MUEC*TB2 + MUE2*TB*TB2 - 
     &            TB*MSq2Diff(1,2))*UStop2(1,1)) + 
     &      2*MSq2Diff(1,2)*
     &       (AtC*MUEC*TB2*Cd(50) + MUE2*TB*TB2*Cd(50) + 
     &         TB*(-(Cd(47)*MSq2Diff(1,1)*UStop2(1,1)) - 
     &            Cd(50)*MSq2Diff(1,1)*UStop2(1,1)))*UStop2(1,2) - 
     &      2*TB*Cd(50)*MSq2Diff(1,1)*MSq2Diff(1,2)*UStop2(1,2)**2)
     &     - 2*AtC*MUEC*TB2*
     &    (Cd(47)*MSq2Diff(1,1)*UStop2(1,1)*
     &       (MSbot2(1) - MSq2Diff(1,2)*UStop2(1,1)) + 
     &      Cd(50)*MSq2Diff(1,2)*MStop2(1)*UStop2(1,2)**2) + 
     &   2*MTy*(2*M_3*AtC*MUEC*TB2*UCStopC(1,3)*
     &       (Cd(47)*MSq2Diff(1,1)*UStop2(1,1) - 
     &         Cd(50)*MSq2Diff(1,2)*UStop2(1,2)) + 
     &      AtC*UCStop(1,3)*
     &       (Cd(47)*MSq2Diff(1,1)*
     &          (-(TB*(MGlpT2 - MSbot2(1) + MSq2Diff(1,2))) + 
     &            2*M_3C*MUEC*TB2*UStop2(1,1)) + 
     &         Cd(50)*MSq2Diff(1,2)*
     &          (TB*MGlpTmSt2(1) - 2*M_3C*MUEC*TB2*UStop2(1,2)))+
     &        Cd(50)*MSq2Diff(1,2)*
     &       (UCStopC(1,3)*
     &          (MGl2*MUEC*TB2 - MUEC*TB2*MStop2(1) + 
     &            MTy2*(MUEC*TB2 - 2*M_3*TB*UStop2(1,1)) - 
     &            2*M_3*MUE2*TB*TB2*UStop2(1,2) + 
     &            TB*MSq2Diff(1,1)*(-M_3 + 2*M_3*UStop2(1,2)))+
     &           UCStop(1,3)*
     &          (MGl2*MUE*TB2 - MUE*TB2*MStop2(1) + 
     &            MTy2*(MUE*TB2 - 2*M_3C*TB*UStop2(1,1)) - 
     &            2*M_3C*MUE2*TB*TB2*UStop2(1,2) + 
     &            TB*MSq2Diff(1,1)*(-M_3C + 2*M_3C*UStop2(1,2))
     &            )) + Cd(47)*MSq2Diff(1,1)*
     &       (UCStopC(1,3)*
     &          (-(MGl2*MUEC*TB2) + 
     &            MUEC*TB2*(MSbot2(1) - MSq2Diff(1,2)) + 
     &            M_3*TB*MSq2Diff(1,2) + 
     &            2*M_3*TB*(MUE2*TB2 - MSq2Diff(1,2))*
     &             UStop2(1,1) + 
     &            MTy2*(-(MUEC*TB2) + 2*M_3*TB*UStop2(1,2))) + 
     &         UCStop(1,3)*
     &          (-(MGl2*MUE*TB2) + 
     &            MUE*TB2*(MSbot2(1) - MSq2Diff(1,2)) + 
     &            M_3C*TB*MSq2Diff(1,2) + 
     &            2*M_3C*TB*(MUE2*TB2 - MSq2Diff(1,2))*
     &             UStop2(1,1) + 
     &            MTy2*(-(MUE*TB2) + 2*M_3C*TB*UStop2(1,2))))) + 
     &   MTy2*(Cd(50)*MSq2Diff(1,2)*
     &       (-4*TB2*(M_3C*MUE*UCStop(1,3)**2 + 
     &            M_3*MUEC*UCStopC(1,3)**2) + 
     &         2*MGl2*TB*UStop2(1,1) + 
     &         2*AtC*MUEC*TB2*UStop2(1,2) + 
     &         TB*(-4*M_3C*AtC*UCStop(1,3)**2 - 
     &            2*MStop2(1)*UStop2(1,1)**2 + 
     &            2*MUE2*TB2*UStop2(1,2) + 
     &            MSq2Diff(1,1)*
     &             (UStop2(1,1)*(1 - 2*UStop2(1,2)) - 
     &               2*UStop2(1,2)**2))) - 
     &      Cd(47)*MSq2Diff(1,1)*
     &       (TB*(4*M_3C*AtC*UCStop(1,3)**2 - 
     &            2*MUE2*TB2*UStop2(1,1)**2 + 
     &            2*UStop2(1,2)*
     &             (MSbot2(1) + 2*M_3*AtC*UStop2(1,1) - 
     &               MGl2*UStop2(1,2)) + 
     &            MSq2Diff(1,2)*
     &             (UStop2(1,1)*(-1 + 2*UStop2(1,1)) - 
     &               2*UStop2(1,2)**2)) + 
     &         2*TB2*(-(AtC*MUEC*UStop2(1,1)**2) + 
     &            2*(MUE*(M_3C*UCStop(1,3)**2 + 
     &                  M_3*UStop2(1,1)*UStop2(1,2)) + 
     &               MUEC*(M_3*UCStopC(1,3)**2 + 
     &                  M_3C*UStop2(1,1)*UStop2(1,2))))))

        Opt(18) = 2*M_3*MTy2*Cd(9004) - 
     &   M_3C*Cd(4735)*MStop2(3)*UCStop(1,3)**2 + 
     &   Cd(9002)*UStop2(1,1) + 
     &   2*MTy*Cd(9011)*UCStop(1,3)*UStop2(1,2) + 
     &   2*M_3*MStop2(3)*
     &    (-((MGl2 - MStop2(2))**3*MStop2(2)) + 
     &      MTy2**2*(Cd(4758) + 
     &         MTy2*(-Cd(4759) + MTy2*UStop2(1,2)**2)))

	Cd(9992) = MUE*Cd(9987) + MUEC*Cd(9991)

        Cd(10041) = Cd(10040)*UCStopC(1,3) + 
     &   AtC*Cd(10031)*Cr(5)*UStop2(1,1)*UStop2(1,2)

        Cd(10132) = (-2*Cd(47)**2*Cd(10129) + 
     &     Cd(50)*Cd(10131)*MStop2(3))/(Cd(47)**2*Cd(50)*MStop2(3))

	Cd(10230) = (Cd(10220) + Cd(10229))/(Cd(47)**2*MStop2(1)**2)

	Cd(10268) = (Cd(10260) + Cd(10267))/(Cd(50)**2*MStop2(2)**2)

        Cd(10275) = (Cd(10270) - Cd(179)*Cd(2681)*Cd(10274))/
     &   (Cd(50)**2*MStop2(2)**2)

	Cd(10531) = 4*(MUE*Cd(10527) + MUEC*Cd(10530))

        Cd(10599) = (Cd(47)*Cd(10578) + 8*Cd(10598)*MStop2(3))/
     &   (Cd(47)**2*MStop2(1)*MStop2(3))

	Cd(10823) = 4*(MUE*Cd(10814) + MUEC*Cd(10822))

        Cd(10872) = (Cd(50)*Cd(10823) + 4*Cd(10871)*MStop2(3))/
     &   (Cd(50)**2*MStop2(2)*MStop2(3))

	Cd(107) = 2*M_3C*MTy*Cd(82) + Cd(106)*UCStopC(1,3)

        Cd(11104) = 4*Cd(11102)*UStopC(1,2) - 
     &   AtC**2*Cd(11103)*UUStopC(1,2)

	Cd(108) = MUE*Cd(76) + MUEC*Cd(107)

        Cd(11152) = Cd(11150)*UStopC(1,1) + 
     &   AtC**2*Cd(11151)*UCStop(1,3)*UStop2(1,1)*UStopC(1,2)

	Cd(11153) = Cd(11104) + 2*Cd(11152)

        Cd(11156) = (Cd(11153) - 
     &     Cd(11155)*MStop2(3)**2*UUStopC(1,1))/MStop2(3)**2

	Cd(11195) = Cd(11193) - 4*Cd(11164)*MStop2(3)*UStop2(1,1)

	Cd(11241) = Cd(11006) - Cd(11007)

        Cd(11242) = 2*(Cd(11240)*MStop2(3)*UStopC(1,1) + 
     &     Cd(11241)*UUStopC(1,3))

        Cd(11254) = 4*Cd(11252)*UStopC(1,2) - 
     &   AtC**2*Cd(11253)*UUStopC(1,2)

        Cd(11257) = Cd(11255)*UStopC(1,1) + 
     &   AtC**2*Cd(11256)*UCStop(1,3)*UStop2(1,1)*UStopC(1,2)

	Cd(11258) = Cd(11254) + 2*Cd(11257)

        Cd(11455) = (Cd(11454) + 
     &     4*MUE**2*MStop2(3)*UCStop(1,3)**2*UStop2(1,1)**2)/
     &   MStop2(3)

        Cd(11456) = (4*Cd(11436)*MStop2(2) + 
     &     Cd(11455)*MStop2(3)**2)/(MStop2(2)**2*MStop2(3)**2)

        Cd(11552) = 2*Cd(11551)*MStop2(2) - 
     &   2*Cd(12)*Cd(13)*Cd(11549)*MStop2(3)**2

	Cd(11553) = SA2*Cd(11526) + S2A*Cd(11552)

        Cd(11554) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*(-Cd(11553) + 
     &        CA2*Cd(11456)*MStop2(2)**2*MStop2(3)**2))/
     &    (SB2*SW2*MStop2(2)**2*MStop2(3)**2)

        Cd(11558) = Cd(11556)*MStop2(3)*UCStopC(1,3) + 
     &   4*Cd(11557)*MStop2(2)*UStopC(1,1)

	Cd(11563) = 2*Cd(11562) - MUEC*Cd(11558)*UCStopC(1,3)

        Cd(11690) = 1/(128.D0*Pi**3)*
     &   (Alfas2L*(-4*S2A*Cd(11651) + 
     &        Cd(11689)*MStop2(1)*MStop2(2)*MStop2(3)**2))/
     &    (SB2*SW2*MStop2(1)*MStop2(2)*MStop2(3)**2)

        Cd(11833) = (Cd(47)*Cd(11826) + 8*Cd(11832)*MStop2(3))/
     &   (Cd(47)**2*MStop2(1)*MStop2(3))

        Cd(11958) = (-(Cd(50)*Cd(11945)) + 8*Cd(11957))/
     &   (Cd(50)**2*MStop2(2)*MStop2(3))

	Cd(1196) = MUEC*Cd(1178) + MUE*Cd(1185)*UCStop(1,3)**2

        Cd(12049) = (2*Cd(47)*Cd(50)*Cd(9955)*Cd(12048) + 
     &     Cd(10313)*Cd(12045)*MStop2(3))/MStop2(3)

	Cd(1198) = MUE*Cd(1178) + MUEC*Cd(1195)*UCStopC(1,3)**2

	Cd(12106) = M_3*Cd(12092) - M_3C*Cd(12105)

	Cd(12209) = 2*Cd(12201) - MTy*Cd(12208)*UStop2(1,1)

	Cd(12249) = Cd(12233)*UUStopC(1,1) + 2*Cd(12240)*UUStopC(1,3)

	Cd(12251) = Cd(50)*Cd(12249) + Cd(47)**2*Cd(12250)

        Cd(12410) = Cd(47)*Cd(12409) + 
     &   Cd(10141)*Cd(12370)*MStop2(1)*MStop2(3)

        Cd(12518) = (4*(-2*MGlmT2*Cd(50)*Cd(12517) + 
     &       Cd(12516)*MStop2(3)))/(Cd(50)**2*MStop2(3))

	Cd(12553) = MUE*Cd(12550) + MUEC*Cd(12552)

        Cd(12561) = (4*MGlmT2*Cd(50)*Cd(12560) - 
     &     Cd(12559)*MStop2(3)*UCStop(1,3))/MStop2(3)

        Cd(12565) = (4*MGlmT2*Cd(50)*Cd(12564) - 
     &     Cd(12559)*MStop2(3)*UCStopC(1,3))/MStop2(3)

        Cd(12578) = (Cd(12573) + Cd(12577))/
     &   (Cd(50)**2*MStop2(2)**2*MStop2(3))

        Cd(12637) = (Cd(50)**2*Cd(12636) - 
     &     Cd(13)*Cd(47)**2*Cd(12630)*UCStopC(1,3))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(12677) = 2*M_3*Cd(12092) - 2*M_3C*Cd(12105)

        Cd(12793) = Cd(20)*Cd(12777) - 
     &   2*Cd(12792)*MStop2(1)*MStop2(2)*MStop2(3)

	Cd(12830) = Cd(17)*Cd(12818) - Cd(12806)*MStop2(3)

	Cd(12831) = 8*Cd(20)*Cd(12829) - 4*Cd(12830)*MStop2(2)*MStop2(3)

        Cd(12944) = Cd(20)*Cd(12923) + 
     &   2*Cd(17)*Cd(12943)*MStop2(1)*MStop2(2)*MStop2(3)

	Cd(13043) = Cd(12992)*MStop2(1) + Cd(13042)*MStop2(2)

        Cd(13074) = 1/(64.D0*Pi**4)*
     &   (2 + Pi**4*(MTy2*TB2*Cd(13073)))/(MTy2*TB2)

        Cd(13099) = 1/(64.D0*Pi**4)*
     &   (MTy2*Cd(13098) + Cd(13092)*MStop2(3))/(MTy2*MStop2(2))

	Cd(13101) = (MGl*Cd(13099) - Cd(13100))/(MGl*TB2*MStop2(3))

	Cd(128) = Cd(52)*Cd(115) + Cd(127)*T134(MGl2,MTy2,MStop2(1))

	Cd(13153) = Cd(22)*Cd(13151) + 4*Cd(13152)*MStop2(1)

	Cd(13306) = 2*Cd(13305) + Cd(13042)*Cr(34)

	Cd(13307) = Cd(13288)*MStop2(1) + Cd(13306)*MStop2(2)

	Cd(13320) = Cd(13315) + Cd(13319)

	Cd(13339) = -2*MGl2*Cd(13338) + Cd(13074)*Cr(65)

	Cd(13349) = 2*Cd(13348) + Cd(13086)*Cr(13)

        Cd(13353) = 1/(128.D0*Pi**4)*
     &   (MTy2*Cd(13098) + Cd(13092)*MStop2(3))/(MTy2*MStop2(2))

	Cd(13355) = (MGl*Cd(13353) - Cd(13354))/(MGl*TB2*MStop2(3))

	Cd(13437) = Cd(16)*Cd(13412) + 4*Cd(13436)*MStop2(2)

        Cd(13851) = 1/2.D0*(A0delStop(2)*Cd(13758)) - 
     &   Cd(13850)*Cr(34)*MStop2(2)

        Cd(13870) = 1/2.D0*(A0delStop(1)*Cd(13804)) - 
     &   Cd(13869)*Cr(13)*MStop2(1)

	Cd(13878) = 2*Cd(13850) - Cd(13758)*Cr(34)

	Cd(13879) = 2*Cd(13869) - Cd(13804)*Cr(13)

	Cd(14116) = -2*Cd(14112)*MSbot2(1) + TB2*Cd(14067)*MStop2(1)

        Cd(14208) = -(Cd(14207)*UStop2(1,1)) + 
     &   2*Cd(50)*MSq2Diff(1,2)*UStop2(1,2)

	Cd(14225) = 2*Cd(14224) + At*Cd(14220)*MSq2Diff(1,1)*UStop2(1,1)

        Cd(14226) = (2*Cd(50)*Cd(10169)*MSq2Diff(1,2)**2 - 
     &     Cd(14225)*MStop2(3))/(Cd(50)*MSq2Diff(1,2)**2*MStop2(3))

        Cd(14230) = 2*Cd(14229) + 
     &   AtC*Cd(14220)*MSq2Diff(1,1)*UStop2(1,1)

        Cd(14232) = (2*Cd(50)*Cd(10176)*MSq2Diff(1,2)**2 - 
     &     Cd(14230)*MStop2(3))/(Cd(50)*MSq2Diff(1,2)**2*MStop2(3))

        Cd(14297) = Cd(14295)*MSq2Diff(1,1)**2 + 
     &   MTy*Cd(14296)*MSq2Diff(1,2)**2

        Cd(14298) = Cd(14267)*Cd(14273)*MSq2Diff(1,1)*
     &    MSq2Diff(1,2) - Cd(14297)*T134(0.D0,MGl2,MSbot2(1))

	Cd(14352) = Cd(14106)*MSq2Diff(1,1)**2 + Cd(14351)*MSq2Diff(1,2)

	Cd(1466) = -(Cd(1465)*MStop2(3)) + Cd(1464)*UStop2(1,1)

	Cd(1487) = AtC*Cd(1455) + 2*MTy*Cd(1479)*MStop2(3)*UCStopC(1,3)

	Cd(148) = Cd(64)*Cd(135) + Cd(147)*T134(MGl2,MTy2,MStop2(2))

	Cd(149) = -(Cd(50)*Cd(128)) + Cd(47)*Cd(148)

        Cd(1522) = (Cd(50)**2*Cd(1507) + Cd(47)**2*Cd(1521))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(1524) = -2*Cd(1492)*Cd(1523) + Cd(1522)*UStop2(1,1)

	Cd(150) = AtC*Cd(76) + At*Cd(107)

	Cd(151) = (Cd(150) + 2*Cd(149)*MStop2(3))/MStop2(3)

	Cd(1580) = Cd(64)*Cd(1579) + 2*Cd(1508)*Cd(1528)*MStop2(3)

        Cd(1581) = Cd(1580)*MStop2(3) - 
     &   Cd(1555)*T134(MGl2,MTy2,MStop2(2))

	Cd(1629) = Cd(52)*Cd(1628) - 2*Cd(1493)*Cd(1609)*MStop2(3)

        Cd(1630) = Cd(1629)*MStop2(3) - 
     &   Cd(1605)*T134(MGl2,MTy2,MStop2(1))

        Cd(1631) = (-(Cd(47)**2*Cd(1581)) + Cd(50)**2*Cd(1630))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(1632) = 2*M_3C*MTy*Cd(1524) + Cd(1631)*UCStopC(1,3)

	Cd(1784) = Cd(1783)*MStop2(3) + Cd(1780)*UStop2(1,2)

        Cd(1883) = (Cd(47)**2*Cd(1861) + Cd(50)**2*Cd(1882))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(1936) = Cd(1883) - 8*M_3*MTy*Cd(1492)*UStop2(1,1)

	Cd(1938) = Cd(1936)*UStop2(1,1) + Cd(1937)*UStopC(1,2)

        Cd(1940) = (-2*Cd(1523)*Cd(1838) + 
     &     Cd(47)*Cd(50)*Cd(1522)*UStop2(1,1))/(Cd(47)*Cd(50))

	Cd(1941) = 2*M_3C*MTy*Cd(1940) + Cd(1631)*UCStopC(1,3)

        Cd(1995) = -(Cd(1994)*MStop2(3)) + 
     &   Cd(1969)*T134(MGl2,MTy2,MStop2(1))

        Cd(2047) = (Cd(50)**2*Cd(1995) + Cd(47)**2*Cd(2046))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(2078) = Cd(2077)*MStop2(3) + 
     &   Cd(2061)*T134(MGl2,MTy2,MStop2(1))

        Cd(2112) = (Cd(50)**2*Cd(2078) + Cd(47)**2*Cd(2111))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(2146) = Cd(52)*Cd(2145) - 
     &   2*Cd(86)**2*Cd(87)*Cd(1493)*MStop2(3)

        Cd(2147) = Cd(2146)*MStop2(3) + 
     &   Cd(2132)*T134(MGl2,MTy2,MStop2(1))

        Cd(2185) = (Cd(50)**2*Cd(2147) + Cd(47)**2*Cd(2184))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(2239) = Cd(2238)*MStop2(3) - 
     &   Cd(2219)*T134(MGl2,MTy2,MStop2(1))

	Cd(250) = 2*SA*Cd(241) + CA*Cd(249)

        Cd(2610) = Cd(2609)*UCStopC(1,3) + 
     &   M_3C*MTy*Cd(152)*Cd(2597)*UStop2(1,1)*UStop2(1,2)

        Cd(2633) = Cd(2632)*UCStopC(1,3) + 
     &   M_3C*MTy*Cd(152)*Cd(2597)*UStop2(1,1)*UStop2(1,2)

	Cd(2642) = MTy*Cd(2639)*UStop2(1,2) + Cd(2641)*UStopC(1,1)

        Cd(2643) = -(Cd(2642)*UCStopC(1,3)) + 
     &   M_3C*MTy2*Cd(2638)*UStop2(1,1)

	Cd(2671) = 2*Cd(2643) + AtC*Cd(2616)*UStop2(1,1)

	Cd(2675) = MTy*Cd(2639)*UStop2(1,2) + Cd(2674)*UStopC(1,1)

        Cd(2676) = -(Cd(2675)*UCStopC(1,3)) + 
     &   M_3C*MTy2*Cd(2638)*UStop2(1,1)

	Cd(2677) = 2*Cd(2670) - AtC*Cd(2624)*UCStop(1,3)

        Cd(2693) = -(Cd(2692)*UCStopC(1,3)) + 
     &   M_3C*MTy*Cd(179)*Cd(2681)*UStop2(1,1)*UStop2(1,2)

        Cd(2709) = -(Cd(2708)*UCStopC(1,3)) + 
     &   M_3C*MTy*Cd(179)*Cd(2681)*UStop2(1,1)*UStop2(1,2)

	Cd(2718) = MTy*Cd(2715)*UStop2(1,2) + Cd(2717)*UStopC(1,1)

        Cd(2719) = Cd(2718)*UCStopC(1,3) + 
     &   M_3C*MTy2*Cd(2714)*UStop2(1,1)

	Cd(2747) = 2*Cd(2719) + AtC*Cd(2697)*UStop2(1,1)

	Cd(2751) = MTy*Cd(2715)*UStop2(1,2) + Cd(2750)*UStopC(1,1)

        Cd(2752) = Cd(2751)*UCStopC(1,3) + 
     &   M_3C*MTy2*Cd(2714)*UStop2(1,1)

	Cd(2753) = 2*Cd(2746) - AtC*Cd(2702)*UCStop(1,3)

	Cd(2783) = -(Cd(2782)*MStop2(1)) + 2*Cd(2610)*MStop2(3)

	Cd(2816) = Cd(2815)*MStop2(1) + Cd(2611)*MStop2(3)*UStop2(1,2)

        Cd(2817) = 2*Cd(2816)*UStop2(1,1) - 
     &   4*M_3*MTy*Cd(2757)*MStop2(1)*UCStopC(1,3)*UStop2(1,2)

        Cd(2843) = MUE**2*Cd(2842)*UCStop(1,3) - 
     &   MUEC**2*Cd(2783)*UCStopC(1,3)

	Cd(2865) = M_3C*MTy*Cd(2848) + Cd(2864)*UCStopC(1,3)

	Cd(2941) = 2*MTy*Cd(2899)*UStop2(1,2) - Cd(2940)*UStopC(1,1)

	Cd(2942) = M_3C*Cd(2886) + Cd(2941)*UCStopC(1,3)

	Cd(3047) = Cd(3044)*MStop2(1) + MTy2*Cd(3046)*UStop2(1,2)

        Cd(3048) = 2*Cd(3039)*MStop2(1)**2*MStop2(3) + 
     &   Cd(3047)*UStop2(1,2)

	Cd(3057) = Cd(3056)*MStop2(3) - M_3C*MTy*Cd(3048)*UCStop(1,3)

	Cd(3093) = 4*Cd(3092) + AtC*Cd(2842)*UCStop(1,3)

	Cd(3107) = 2*MTy*Cd(2899)*UStop2(1,2) - Cd(3106)*UStopC(1,1)

	Cd(3108) = M_3C*Cd(3098) + Cd(3107)*UCStopC(1,3)

	Cd(3109) = 2*Cd(3092) + AtC*Cd(2842)*UCStop(1,3)

	Cd(3147) = M_3C*MTy*Cd(3121) + 2*Cd(3146)*UCStopC(1,3)

	Cd(3148) = Cd(3147)*MStop2(2) - 2*Cd(2693)*MStop2(3)

	Cd(3179) = Cd(3178)*MStop2(2) - Cd(2694)*MStop2(3)*UStop2(1,2)

        Cd(3180) = 2*Cd(3179)*UStop2(1,1) + 
     &   4*M_3*MTy*Cd(3116)*MStop2(2)*UCStopC(1,3)*UStop2(1,2)

        Cd(3199) = MUE**2*Cd(3198)*UCStop(1,3) - 
     &   MUEC**2*Cd(3148)*UCStopC(1,3)

	Cd(3220) = M_3C*MTy*Cd(3201) + Cd(3219)*UCStopC(1,3)

	Cd(322) = Cd(301) + Cd(321)*UStop2(1,1)

	Cd(3306) = 2*MTy*Cd(3261)*UStop2(1,2) - Cd(3305)*UStopC(1,1)

	Cd(3307) = M_3C*Cd(3249) - Cd(3306)*UCStopC(1,3)

	Cd(336) = MUEC*Cd(329) + MUE*Cd(335)

	Cd(3420) = Cd(3419)*MStop2(3) - Cd(3393)*UUStopC(1,1)

	Cd(3467) = 4*Cd(3466) - AtC*Cd(3198)*UCStop(1,3)

	Cd(343) = SA*Cd(241) + CA*Cd(342)

	Cd(3484) = 2*MTy*Cd(3261)*UStop2(1,2) - Cd(3483)*UStopC(1,1)

	Cd(3485) = M_3C*Cd(3249) - Cd(3484)*UCStopC(1,3)

	Cd(3486) = 2*Cd(3466) - AtC*Cd(3198)*UCStop(1,3)

        Cd(3490) = Cd(2815)*UStop2(1,1) - 
     &   2*M_3*MTy*Cd(2757)*UCStopC(1,3)*UStop2(1,2)

	Cd(3521) = Cd(3520) + Cd(3500)*UStop2(1,1)

	Cd(3522) = Cd(3490) + Cd(47)*Cd(2818)*UStop2(1,2)*UStopC(1,2)

        Cd(3523) = -4*Cd(3522)*MStop2(3) + 
     &   2*Cd(3521)*UStop2(1,1)*UStopC(1,1)

	Cd(353) = -(SA*Cd(13)*Cd(252)) + CA*Cd(352)

	Cd(3578) = Cd(3568) + 2*Cd(3577)*UStop2(1,2)

	Cd(3607) = 2*MTy*Cd(3578)*UStop2(1,2) + Cd(3606)*UStopC(1,1)

	Cd(3608) = 2*M_3C*Cd(3558)*MStop2(3) - Cd(3607)*UCStopC(1,3)

        Cd(3619) = Cd(3618)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(2757)*MStop2(3)*UStop(1,2)

        Cd(3620) = 2*Cd(3577)*UCStop(1,3)*UStop2(1,2) - 
     &   Cd(3619)*UStopC(1,2)

        Cd(3621) = 4*M_3*MTy*Cd(47)*MStop2(1)*MStop2(3) - 
     &   Cd(3620)*UStop2(1,2)

	Cd(3632) = 2*MTy*Cd(3621) - 2*Cd(3631)

	Cd(3633) = 2*Cd(3632) - AtC*Cd(3550)*UCStop(1,3)

	Cd(3634) = Cd(3608) + At*Cd(2782)*MStop2(3)*UCStopC(1,3)

        Cd(3658) = Cd(3657)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(3103)*MStop2(3)*UStop(1,2)

        Cd(3745) = M_3*MTy*Cd(3744)*UUStop(1,2) + 
     &   Cd(47)**2*MStop2(1)*UUStop(2,3)

        Cd(3755) = 2*Cd(3178)*UStop2(1,1) + 
     &   4*M_3*MTy*Cd(3116)*UCStopC(1,3)*UStop2(1,2)

	Cd(3793) = Cd(3755) - 2*Cd(50)*Cd(3181)*UStop2(1,2)*UStopC(1,2)

        Cd(3794) = 2*(Cd(3793)*MStop2(3) + 
     &     Cd(3792)*UStop2(1,1)*UStopC(1,1))

        Cd(3795) = Cd(3794)*MStop2(2) + 
     &   Cd(3769)*UStop2(1,1)**2*UStopC(1,1)

        Cd(3819) = 2*MUE**2*Cd(3818)*UCStop(1,3) + 
     &   MUEC**2*Cd(3147)*MStop2(3)*UCStopC(1,3)

	Cd(3820) = MUE2*Cd(3795) + Cd(3819)*MStop2(2)

        Cd(3882) = Cd(3881)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(3116)*MStop2(3)*UStop(1,2)

        Cd(3883) = 2*Cd(3838)*UCStop(1,3)*UStop2(1,2) - 
     &   Cd(3882)*UStopC(1,2)

        Cd(3884) = 4*M_3*MTy*Cd(50)*MStop2(2)*MStop2(3) - 
     &   Cd(3883)*UStop2(1,2)

	Cd(3893) = Cd(3892)*MStop2(3) + MTy*Cd(3890)*UCStop(1,3)

        Cd(3894) = -2*M_3C*Cd(3872)*MStop2(3)*UCStop(1,3)**2 + 
     &   2*Cd(3893)*UStop2(1,1)

        Cd(3897) = Cd(3871) + 
     &   At*Cd(3147)*MStop2(2)*MStop2(3)*UCStopC(1,3)

	Cd(391) = 2*CA*Cd(241) - SA*Cd(249)

	Cd(4005) = Cd(3981) + Cd(4004)

        Cd(4027) = -(Cd(4026)*UCStopC(1,3)) + 
     &   M_3C*MTy*Cd(4017)*UStop2(1,1)*UStop2(1,2)

	Cd(395) = CA*Cd(241) - SA*Cd(342)

	Cd(396) = CA*Cd(13)*Cd(252) + SA*Cd(352)

        Cd(4045) = -(Cd(4044)*UCStopC(1,3)) + 
     &   M_3C*MTy*Cd(4017)*UStop2(1,1)*UStop2(1,2)

	Cd(4111) = 2*Cd(4077) - AtC*Cd(4040)*UCStop(1,3)

        Cd(4127) = Cd(4126)*UCStopC(1,3) + 
     &   M_3C*MTy*Cd(4116)*UStop2(1,1)*UStop2(1,2)

	Cd(4201) = ((-4*Cd(4199) + MTy*Cd(4200))*MStop2(3)**2)/MTy

	Cd(4212) = 2*Cd(4211) - AtC*Cd(4135)*UCStop(1,3)

	Cd(4238) = M_3C*MTy*Cd(4221) + 2*Cd(4237)*UCStopC(1,3)

	Cd(4239) = Cd(4238)*MStop2(1) + 2*Cd(4027)*MStop2(3)

        Cd(4278) = -(Cd(4277)*MStop2(3)) + 
     &   Cd(4266)*UStop2(1,1)*UStopC(1,1)

        Cd(4301) = MUE**2*Cd(4300)*UCStop(1,3) + 
     &   MUEC**2*Cd(4239)*UCStopC(1,3)

	Cd(4302) = 2*MUE2*Cd(4278) - Cd(4301)*MStop2(3)

	Cd(4322) = M_3C*MTy*Cd(4305) + Cd(4321)*UCStopC(1,3)

	Cd(4539) = Cd(4538)*MStop2(3) - MTy*Cd(4519)*UUStopC(1,1)

	Cd(4540) = At*Cd(4278) + 2*Cd(4435)

	Cd(4547) = M_3C*MTy*Cd(4305) + Cd(4546)*UCStopC(1,3)

	Cd(4557) = 2*Cd(4435) - AtC*Cd(4300)*MStop2(3)*UCStop(1,3)

	Cd(4558) = 2*Cd(4556) + At*Cd(4547)*UCStopC(1,3)

	Cd(4587) = Cd(4586)*MStop2(2) - 2*Cd(4127)*MStop2(3)

        Cd(4624) = -(Cd(4623)*MStop2(3)) + 
     &   Cd(4611)*UStop2(1,1)*UStopC(1,1)

        Cd(4637) = MUE**2*Cd(4636)*UCStop(1,3) + 
     &   MUEC**2*Cd(4587)*UCStopC(1,3)

	Cd(4638) = 2*MUE2*Cd(4624) + Cd(4637)*MStop2(3)

	Cd(4658) = M_3C*MTy*Cd(4641) - Cd(4657)*UCStopC(1,3)

	Cd(4877) = Cd(4876)*MStop2(3) - MTy*Cd(4851)*UUStopC(1,1)

	Cd(4886) = M_3C*MTy*Cd(4641) - Cd(4885)*UCStopC(1,3)

	Cd(4902) = Cd(4901)*MStop2(3) + Cd(4900)*UStop2(1,1)*UStopC(1,1)

	Cd(4914) = 2*Cd(4913) - AtC*Cd(4636)*MStop2(3)*UCStop(1,3)

	Cd(4915) = 2*Cd(4898) - At*Cd(4886)*UCStopC(1,3)

	Cd(5012) = 2*M_3C*MTy*Cd(4942) - 2*Cd(5011)*UCStopC(1,3)

	Cd(5221) = Cd(5219) + 2*Cd(5079)*UStop2(1,1)

	Cd(5222) = Cd(5221)*UStop2(1,1) - Cd(5131)*UStop2(1,2)**2

	Cd(5383) = 2*Cd(44)**2*Cd(48)**2*Cd(288) + MGl2*Cd(5382)

        Cd(5426) = Cd(5425)*MStop2(1)**2*MStop2(2)**2 + 
     &   Cd(5383)*UStop2(1,2)**3

        Cd(5740) = 3*Cd(5736)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(5739)*MStop2(3)**2*UStop(1,2)

	Cd(5754) = Cd(5524) + 2*MTy2*Cd(5469)*UStop2(1,1)

	Cd(5948) = 2*Cd(5942) - MGl2*Cd(5947)

	Cd(5949) = Cd(44)**2*Cd(48)**2*Cd(2521) - MGl2*Cd(5948)

        Cd(5953) = 2*Cd(5949)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(5952)*MStop2(3)*UStop(1,2)

	Cd(5979) = 2*MGl2**3*Cd(1511) - Cd(5978)*MStop2(3)

        Cd(6049) = 2*Cd(5736)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(5264)*UStop(1,2)

        Cd(6210) = -(Cd(6209)*UCStop(1,3)) + 
     &   4*M_3*MTy*Cd(6202)*UStop2(1,2)

	Cd(6211) = 4*M_3*Cd(6167) - MTy*Cd(6210)

	Cd(6378) = 2*M_3C*MTy*Cd(6347) + Cd(6377)*UCStopC(1,3)

	Cd(6597) = 2*Cd(6595) - MTy*Cd(6596)

	Cd(6598) = M_3C*Cd(6432) - Cd(6597)*UCStopC(1,3)

	Cd(6677) = Cd(6676)*UStop2(1,1) - Cd(6673)*UStop2(1,2)

        Cd(6678) = Cd(6677)*UStopC(1,1) - 
     &   Cd(47)*Cd(50)*Cd(6675)*UStop2(1,2)*UStopC(1,2)

	Cd(6967) = Cd(6966)*UStop(1,2) - 2*M_3C*Cd(6965)*UUStop(1,1)

	Cd(6982) = Cd(6967) + MTy*Cd(6981)*UStop2(1,2)

	Cd(7047) = 2*M_3C*MTy*Cd(7023) - Cd(7046)*UCStopC(1,3)

	Cd(7166) = Cd(7145) - Cd(7165)*UStop2(1,1)

	Cd(7167) = Cd(7166)*UStop2(1,1) - Cd(7086)*UStop2(1,2)**2

	Cd(716) = 2*Cd(711) + Cd(692)*UStopC(2,1)

        Cd(718) = -(Cd(717)*Conjugate(dY33eps(-1))) + 
     &   Cd(716)*dY33eps(-1)

        Cd(7343) = Cd(7342)*MStop2(1)**2*MStop2(2)**2 + 
     &   2*Cd(7306)*UStop2(1,2)**3

	Cd(732) = 2*Cd(731) + Cd(728)*Cr(34)

	Cd(733) = 1/2.D0*(A0delStop(2)*Cd(718)) + Cd(732)*MStop2(2)

	Cd(7542) = 2*Cd(44)*Cd(48)*Cd(7538) - MGl2*Cd(7541)

	Cd(7543) = 5*Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(7542)

	Cd(7544) = 2*MGl2**5*Cd(7527) + Cd(7543)*MStop2(3)

        Cd(7548) = 2*MTy2*
     &   (Cd(7544)*UStop(1,1) - 
     &     4*M_3*MTy*Cd(7547)*MStop2(3)**2*UStop(1,2))

	Cd(7555) = Cd(7422) + 2*MTy2*Cd(7371)*UStop2(1,1)

	Cd(748) = 2*Cd(747) + Cd(745)*Cr(13)

	Cd(749) = 1/2.D0*(A0delStop(1)*Cd(742)) + Cd(748)*MStop2(1)

	Cd(7628) = Cd(7601) - MTy2*Cd(7577)*UStop2(1,1)

	Cd(7629) = Cd(7627)*MStop2(3) + Cd(7628)*UStop2(1,1)

        Cd(7630) = Cd(47)*Cd(50)*Cd(5805)*MStop2(3) + 
     &   2*Cd(7629)*UStop2(1,1)

        Cd(7750) = Cd(7730)*UStop2(1,1) + 
     &   2*Cd(7743)*MStop2(3)*UStop2(1,2)

	Cd(7801) = Cd(5736)*UStop(1,1) + 2*M_3*MTy*Cd(7203)*UStop(1,2)

	Cd(7933) = 2*Cd(7931) - MGl2*Cd(7932)

	Cd(7934) = Cd(44)**2*Cd(48)**2*Cd(89)*Cd(7647) - MGl2*Cd(7933)

        Cd(7935) = -(Cd(7934)*UCStop(1,3)) + 
     &   2*M_3*MTy*Cd(6202)*UStop2(1,2)

	Cd(7936) = 2*M_3*Cd(7912) - MTy*Cd(7935)

        Cd(7964) = Cd(302)*Cd(7963)*MStop2(3) + 
     &   M_3C*MTy*Cd(7899)*UUStop(1,1)

	Cd(8044) = Cd(8038) - Cd(8043)*UStop2(1,1)

        Cd(8050) = Cd(8044)*MStop2(3) - 
     &   2*M_3*MTy*Cd(8049)*UCStopC(1,3)

	Cd(8051) = 2*M_3C*MTy*Cd(8019) - Cd(8050)*UCStopC(1,3)

	Cd(8226) = M_3C*Cd(8101) + Cd(8225)*UCStopC(1,3)

	Cd(8240) = Cd(8236) + Cd(8239)*UStopC(1,1)

        Cd(8271) = Cd(8268)*MStop2(3)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(8270)*UStop(1,2)

        Cd(8277) = Cd(8276)*MStop2(3)*UStop(1,1) + 
     &   M_3*MTy*Cd(8214)*UStop(1,2)

        Cd(8288) = 2*M_3C*MTy2*Cd(8287)*UUStop(1,1) - 
     &   MTy*Cd(8275)*UUStop(1,3)

        Cd(8289) = Cd(8288)*MStop2(3) + 
     &   2*MTy*Cd(8277)*UStop(1,2)*UStop2(1,1)

        Cd(8343) = Cd(7045)*UCStopC(1,3) + 
     &   2*M_3C*MTy*Cd(8238)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(8351) = Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(8350)

	Cd(8352) = MGl2**5*Cd(8348) - Cd(8351)*MStop2(3)

	Cd(8374) = Cd(8373) + Cd(7164)*UStop2(1,1)

	Cd(8377) = 2*MGl2**5*Cd(8376) - Cd(8351)*MStop2(3)

        Cd(8378) = Cd(8377)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(8238)*MStop2(3)**2*UStop(1,2)

	Cd(8380) = 2*MGl2**5*Cd(8376) - Cd(8379)*MStop2(3)

        Cd(8381) = 2*M_3C*MTy*Cd(8238)*MStop2(3)**2*
     &    UUStop(1,1) - Cd(8380)*UUStop(1,3)

        Cd(8382) = Cd(8378)*UStop2(1,1)*UStopC(1,2) - 
     &   Cd(8381)*UUStopC(1,2)

        Cd(8397) = Cd(8396)*UCStopC(1,3) + 
     &   2*M_3C*MTy*Cd(8238)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(8454) = Cd(8352)*UStop(1,1) + 
     &   M_3*MTy*Cd(8453)*MStop2(3)**2*UStop(1,2)

	Cd(8463) = Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(8462)

        Cd(8469) = 8*MTy*MTy2*Cd(8463)*UStop(1,1) - 
     &   M_3*Cd(8468)*MStop2(3)*UStop(1,2)

	Cd(8473) = Cd(8442) - 8*MTy2**2*Cd(8420)*UStop2(1,1)

        Cd(8474) = M_3C*Cd(8473)*UCStop(1,3) - 
     &   8*M_3*MTy2**2*Cd(8440)*UCStopC(1,3)*UStop2(1,2)

	Cd(8476) = Cd(8469) - 2*MTy*MTy2*Cd(8450)*UStop2(1,2)**2

	Cd(8500) = Cd(8499) - Cd(8489)*UStop2(1,1)**2

	Cd(8501) = Cd(8500)*MStop2(3) - 2*Cd(8487)*UStop2(1,1)

	Cd(8542) = -(Cd(8541)*MStop2(3)) + Cd(8522)*UStop2(1,2)

        Cd(8543) = MTy*Cd(8542) + 
     &   2*M_3*MTy2*Cd(8534)*MStop2(3)*UCStopC(1,3)

        Cd(8544) = Cd(8543)*UCStopC(1,3) + 
     &   2*M_3C*MTy2*Cd(8501)*UStop2(1,1)

	Cd(8559) = Cd(8550) + Cd(8558)*UStop2(1,1)

        Cd(8573) = 4*M_3*MTy2*Cd(8566) + 
     &   MTy*Cd(8572)*MStop2(3)**2*UCStop(1,3)

	Cd(8583) = Cd(8582)*MStop2(3) + Cd(8573)*UStop2(1,1)

        Cd(8596) = -2*M_3C*Cd(8562)*UCStop(1,3)**2 + 
     &   2*M_3*Cd(8489)*MStop2(3)*UStop2(1,1)**3

        Cd(8604) = Cd(8603)*UCStopC(1,3) + 
     &   2*M_3C*MTy*Cd(8600)*UStop2(1,1)

        Cd(8638) = -(Cd(8637)*UCStopC(1,3)) + 
     &   2*M_3C*MTy2*Cd(8618)*MStop2(3)*UStop2(1,1)

        Cd(8641) = Cd(8640)*UCStop(1,3) + 
     &   Cd(8639)*UStop2(1,1)*UStopC(1,2)

	Cd(8642) = Cd(8218)*MStop2(3) - Cd(8271)*UStopC(1,1)

        Cd(8643) = Cd(47)**2*Cd(251)*MGlpTmSt2(2)*MStop2(3) + 
     &   Cd(8642)*UStop2(1,1)

	Cd(8644) = MTy*Cd(8643) - 2*M_3*MTy2*Cd(8222)*UCStopC(1,3)

	Cd(8646) = Cd(8289) + MTy*Cd(8645)*MStop2(3)*UStop2(1,2)

	Cd(8685) = Cd(8684) - Cd(8667)*UStop2(1,1)

        Cd(8687) = 4*Cd(8686)*MStop2(3) + 
     &   2*Cd(8685)*UStop2(1,1)*UStopC(1,1)

	Cd(8708) = M_3C*MTy*Cd(4221) + Cd(8707)*UCStopC(1,3)

	Cd(8755) = MTy*Cd(8725)*UStop2(1,2) - Cd(8754)*UStopC(1,1)

	Cd(8756) = 2*M_3C*Cd(8709)*MStop2(3) + Cd(8755)*UCStopC(1,3)

	Cd(8766) = M_3*Cd(4334) - 2*MTy*Cd(8765)*UUStopC(1,2)

	Cd(8767) = 2*Cd(8766)*MStop2(3) + MTy*Cd(8762)*UCStop(1,3)

	Cd(8780) = Cd(8756) + At*Cd(8708)*MStop2(3)*UCStopC(1,3)

	Cd(8789) = MTy*Cd(8725)*UStop2(1,2) - Cd(8788)*UStopC(1,1)

	Cd(8790) = 2*M_3C*Cd(8709)*MStop2(3) + Cd(8789)*UCStopC(1,3)

	Cd(8881) = Cd(8863) + Cd(8880)

        Cd(8882) = Cd(8881)*MStop2(3) - 
     &   4*MTy*MTy2*Cd(8800)*UStop2(1,1)**2*UStopC(1,1)

	Cd(8883) = Cd(8882)*MStop2(1) + MTy*Cd(8841)*UUStopC(1,1)

	Cd(8884) = AtC*Cd(8687) + 2*Cd(8790)

	Cd(8916) = Cd(8915) - Cd(8900)*UStop2(1,1)

        Cd(8918) = 2*(Cd(8917)*MStop2(3) + 
     &     Cd(8916)*UStop2(1,1)*UStopC(1,1))

	Cd(8939) = M_3C*MTy*Cd(4567) - Cd(8938)*UCStopC(1,3)

	Cd(8956) = Cd(8948) - 4*Cd(8955)*UStop2(1,2)

	Cd(8982) = MTy*Cd(8956)*UStop2(1,2) - Cd(8981)*UStopC(1,1)

	Cd(8983) = 2*M_3C*Cd(8940)*MStop2(3) + Cd(8982)*UCStopC(1,3)

	Cd(8988) = Cd(8987)*MStop2(2) + Cd(8986)*UStop2(1,1)*UStopC(1,1)

	Cd(9023) = MTy*Cd(8956)*UStop2(1,2) - Cd(9022)*UStopC(1,1)

	Cd(9024) = 2*M_3C*Cd(8940)*MStop2(3) + Cd(9023)*UCStopC(1,3)

	Cd(9111) = Cd(9090) + Cd(9110)

        Cd(9112) = 4*MTy2*Cd(9039)*UStop2(1,1)**2*UStopC(1,1) + 
     &   Cd(9074)*UUStopC(1,1)

	Cd(9113) = MTy*Cd(9112) - Cd(9111)*MStop2(3)

	Cd(9114) = AtC*Cd(8918) + 2*Cd(9024)

        Cd(9130) = Cd(52)*Cd(9129) - 
     &   2*Cd(86)**2*Cd(87)*Cd(1493)*MStop2(3)

        Cd(9131) = Cd(9130)*MStop2(3) + 
     &   Cd(9127)*T134(MGl2,MTy2,MStop2(1))

        Cd(9142) = (Cd(50)**2*Cd(9131) + Cd(47)**2*Cd(9141))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(9154) = Cd(9153)*MStop2(3) - 
     &   Cd(9149)*T134(MGl2,MTy2,MStop2(1))

	Cd(9165) = Cd(64)*Cd(9164) - 2*Cd(1508)*Cd(9156)*MStop2(3)

        Cd(9166) = Cd(9165)*MStop2(3) - 
     &   Cd(9161)*T134(MGl2,MTy2,MStop2(2))

        Cd(9537) = M_3C*Cd(9532)*UUStop(1,1) - 
     &   M_3*Cd(9535)*UUStop(1,2)

        Cd(9538) = Cd(9536)*UStopC(1,1) + 
     &   MTy*Cd(9537)*MStop2(3)*UUStopC(1,3)

	Cd(1006) = -(SA2*Cd(961)) + 2*CA2*Cd(996) - S2A*Cd(1005)

	Cd(10237) = SA2*Cd(10216) + S2A*Cd(10230) + CA2*Cd(10236)

        Cd(10256) = Cd(50)**2*Cd(10238) - 
     &   Cd(179)*Cd(2681)*Cd(10255) + 16*Cd(10246)*MStop2(2)**2

	Cd(10465) = SA2*Cd(10366) - S2A*Cd(10434) + CA2*Cd(10464)

	Cd(10728) = AtC*Cd(10685) + At*Cd(10701) - 2*Cd(10727)

        Cd(10748) = 4*MUE2*Cd(10740) - 
     &   2*(MUE**2*Cd(10744)*UCStop(1,3)**2 + 
     &      MUEC**2*Cd(10747)*UCStopC(1,3)**2)

	Cd(10893) = -(SA2*Cd(10801)) - S2A*Cd(10872) + CA2*Cd(10892)

	Cd(11008) = Cd(11241) - 6*Cd(10992)*MStop2(3)**2

        Cd(11009) = MStop2(3)*UStopC(1,1)*
     &    (4*Cd(10988)*MStop2(1) - 
     &      Cd(10990)*MStop2(3)*UStop2(1,2)*UUStopC(1,2)) + 
     &   2*Cd(11008)*UUStopC(1,3)

        Cd(11196) = Cd(11174) + Cd(11194)*UStopC(1,2) + 
     &   Cd(11195)*UUStopC(1,1)

        Cd(11245) = Cd(11242) - Cd(11243)*UStopC(1,2) - 
     &   Cd(11244)*MStop2(3)*UUStopC(1,1)

        Cd(11260) = (Cd(11258) - 
     &     MStop2(3)*(Cd(11259) + 
     &        4*Cd(9794)**2*MStop2(3)*UStop2(1,1)**2)*UUStopC(1,1))
     &    /MStop2(3)**2

	Cd(1120) = SA2*Cd(1031) + CA2*Cd(1090) + S2A*Cd(1119)

        Cd(11374) = 2*Cd(11369) - Cd(11371)*UStopC(2,2) + 
     &   4*Cd(11373)*UUStopC(2,1)

        Cd(11565) = (-(MUEC*Cd(11563)) + 
     &     MUE**2*MStop2(3)*UCStop(1,3)**2*UStop2(1,1)**2 - 
     &     Cd(11564)*MStop2(3)*UUStopC(2,2))/MStop2(3)

	Cd(11618) = Cd(11418) - Cd(11548) + Cd(11617)*MStop2(3)

	Cd(11834) = SA2*Cd(11778) + S2A*Cd(11823) + CA2*Cd(11833)

        Cd(12007) = At*Cd(12006) + 
     &   AtC*(Cd(12003) + AtC*Cd(47)*Cd(50)**2*UCStop(1,3)**2)

        Cd(12008) = -(Cd(47)*Cd(50)*
     &      (Cd(50)*Cd(118) + Cd(47)*Cd(138))*Cd(11996)) + 
     &   2*Cd(12007)

	Cd(1199) = AtC*Cd(1196) + At*Cd(1198) + 2*MTy*Cd(1197)*MStop2(3)

	Cd(1200) = -(SA2*Cd(1133)) + CA2*Cd(1173) - S2A*Cd(1199)

        Cd(12252) = (Cd(12251) + MTy*Cd(50)*Cd(12248)*MStop2(3))/
     &    (Cd(47)**2*Cd(50)*MStop2(3)) + 
     &   (Cd(13)*Cd(12170)*UCStopC(1,3))/Cd(50)**2

	Cd(12337) = -(SA2*Cd(12316)) - S2A*Cd(12334) + CA2*Cd(12336)

        Cd(12365) = (MTy*(Cd(12359) - Cd(12364)) - 
     &     2*Cd(13)*Cd(12356)*MStop2(2))/
     &   (MTy*Cd(50)**2*MStop2(2)**2)

	Cd(12368) = -(SA2*Cd(12355)) - S2A*Cd(12365) + CA2*Cd(12367)

        Cd(12394) = Cd(12393)/(Cd(47)**2*MStop2(3)) + 
     &   (UStop2(1,1)*(-2*Cd(12388) + 
     &        Cd(47)**2*Cd(12384)*UStop2(1,1)))/Cd(47)**2

        Cd(12678) = 2*(Cd(50)*Cd(12676) + 
     &     MUE2*(Cd(9948)*UStop2(1,1)**2 + 
     &        Cd(12109)*UStop2(1,2)**2))

        Cd(13104) = 1/(64.D0*Pi**4)*
     &    (64*Pi**4*(Cd(13103)*MGlpTmStxGlT4(2)**2*MStop2(2)) + 
     &       MGl2*Cd(179)*Cd(12999)*MStop2(3))/
     &     (TB2*MGlpTmStxGlT4(2)**2*MStop2(2)*MStop2(3)) + 
     &   Cd(13101)/MGlpTmStxGlT4(2)

        Cd(13105) = MGl2*Cd(13074) - Cd(13086)*MStop2(1) - 
     &   Cd(13104)*MStop2(2)

        Cd(13274) = 1/(64.D0*Pi**4)*
     &   (TB2*(64*Pi**4*(Cd(52)*Cd(13233)*MGlpTmStxGlT4(1)**2) + 
     &         64*Pi**4*(Cd(13273)*MGlpTmStxGlT4(1)**2*
     &            T134(MGl2,MTy2,MStop2(1)))) - 
     &      Cd(1493)*Cd(12836)**2*Cr(45)*MGlpTmSt4(1))/
     &    (TB2*MGlpTmStxGlT4(1)**2)

        Cd(13357) = 1/(512.D0*Pi**4)*
     &   (512*Pi**4*(Cd(13356)*MGlpTmStxGlT4(2)**2*MStop2(2)**2) + 
     &      512*Pi**4*(TB2*Cd(13355)*MGlpTmStxGlT4(2)*MStop2(2)**2*
     &         MStop2(3)) + Cd(12483)*Cd(12999)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(2)**2*MStop2(2)**2*MStop2(3))

        Cd(13361) = 2*MGl2*Cd(13360) + Cd(13338)*Cr(65) - 
     &   1/2.D0*(A0delGl*Cd(13074))

        Cd(13365) = 1/2.D0*
     &   (4*Cd(13364) + A0delStop(2)*Cd(13104)*MStop2(2) - 
     &      2*Cd(13357)*Cr(34)*MStop2(2)**2)/MStop2(2)

	Cd(13366) = Cd(13361) + Cd(13363) + Cd(13365)

        Cd(13824) = 1/(8.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*(Cd(50)*Cd(13823) + 
     &           Cd(13758)*MSq2Diff(1,1)**2) + 
     &        Cd(50)*Cd(13804)*MSq2Diff(1,2)**2))/
     &    (SW2*TB2*Cd(47)*Cd(50)*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

        Cd(13877) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(Cd(50)*Cd(13870)*MSbot2(1)*MSq2Diff(1,2)**2*
     &         MStop2(2) + 
     &        Cd(47)*MStop2(1)*
     &         (Cd(13851)*MSbot2(1)*MSq2Diff(1,1)**2 + 
     &           Cd(50)*Cd(13876)*MStop2(2))))/
     &    (SW2*TB2*Cd(47)*Cd(50)*MSbot2(1)*MSq2Diff(1,1)**2*
     &      MSq2Diff(1,2)**2*MStop2(1)*MStop2(2))

        Cd(13881) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*(-(Cd(50)*Cd(13880)) + 
     &           Cd(13878)*MSq2Diff(1,1)**2) + 
     &        Cd(50)*Cd(13879)*MSq2Diff(1,2)**2))/
     &    (SW2*TB2*Cd(47)*Cd(50)*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

        Cd(14089) = Cd(9950)*Cd(14088) + 
     &   MTy*MUE*TB2*Cd(47)*Cd(9948)*Cd(14080)*UCStop(1,3) + 
     &   Cd(47)*Cd(14079)*UStop2(1,2)**2

	Cd(1412) = -(S2A*Cd(1341)) + 2*SA2*Cd(1366) + CA2*Cd(1411)

        Cd(14353) = MTy2*
     &    (2*M_3C*Cd(13574)*Cd(14346) - 
     &      2*M_3*Cd(14348)*Cd(14349)) + Cd(14350)*Cd(14352)

        Cd(1785) = Cd(64)*Cd(1784) - 
     &   2*Cd(1508)*Cd(1730)*MStop2(3)**2 - 
     &   Cd(1772)*T134(MGl2,MTy2,MStop2(2))

        Cd(1836) = Cd(52)*Cd(1835) - 
     &   2*Cd(1493)*Cd(1819)*MStop2(3)**2 + 
     &   Cd(1818)*T134(MGl2,MTy2,MStop2(1))

        Cd(2295) = (-(Cd(47)**2*Cd(2294)) + 
     &     Cd(50)**2*(Cd(2239) + 
     &        MTy*Cd(47)**2*Cd(2185)*UStop2(1,1)))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(2298) = (UStop2(1,1)*
     &      (Cd(47)*Cd(50)*Cd(1883) - 
     &        8*M_3*MTy*Cd(1838)*UStop2(1,1)))/(Cd(47)*Cd(50)) + 
     &   Cd(1937)*UStopC(1,2)

        Cd(2449) = Cd(52)*Cd(2448) - 
     &   2*Cd(1493)*Cd(2308)*MStop2(3)**2 - 
     &   Cd(2397)*T134(MGl2,MTy2,MStop2(1))

        Cd(2592) = Cd(64)*Cd(2511) - 
     &   2*Cd(1508)*Cd(2459)*MStop2(3)**2 + 
     &   Cd(2591)*T134(MGl2,MTy2,MStop2(2))

        Cd(2625) = -(MUE**2*Cd(2624)*UCStop(1,3)) - 
     &   MUEC**2*Cd(2610)*UCStopC(1,3) + 
     &   2*MUE2*Cd(2616)*UStop2(1,1)

        Cd(2703) = -(MUE**2*Cd(2702)*UCStop(1,3)) + 
     &   MUEC**2*Cd(2693)*UCStopC(1,3) + 
     &   2*MUE2*Cd(2697)*UStop2(1,1)

        Cd(2819) = Cd(2817) + 
     &   2*Cd(2795)*UStop2(1,1)**2*UStopC(1,1) + 
     &   2*Cd(47)*Cd(2818)*MStop2(1)*UStop2(1,2)*UStopC(1,2)

        Cd(3058) = Cd(3004) + Cd(3057) + 
     &   MTy2*Cd(2952)*UStop2(1,1)**2*UStopC(1,1)

        Cd(3094) = AtC*Cd(3093) + 
     &   At*(4*Cd(2942) + At*Cd(2865)*UCStopC(1,3))

        Cd(3110) = 2*MUEC*Cd(3108) + MUE*Cd(3109) - 
     &   At*MUEC*Cd(2783)*UCStopC(1,3)

        Cd(3183) = -(Cd(3182)*UStop2(1,1)*UStopC(1,1)) + 
     &   MStop2(3)*(Cd(3180) - 
     &      2*Cd(50)*Cd(3181)*MStop2(2)*UStop2(1,2)*UStopC(1,2))

        Cd(3468) = AtC*Cd(3467) + 
     &   At*(4*Cd(3307) + At*Cd(3220)*UCStopC(1,3))

	Cd(3469) = 2*At*AtC*Cd(3183) - 4*Cd(3420) + Cd(3468)*MStop2(3)

        Cd(3487) = 2*MUEC*Cd(3485) + MUE*Cd(3486) + 
     &   At*MUEC*Cd(3148)*UCStopC(1,3)

        Cd(3551) = 2*MUE2*Cd(3523) + 
     &   MUE**2*Cd(3550)*UCStop(1,3) - 
     &   2*MUEC**2*Cd(2782)*MStop2(3)*UCStopC(1,3)

	Cd(3635) = Cd(424)*Cd(3523) - MUE*Cd(3633) - 2*MUEC*Cd(3634)

        Cd(3659) = Cd(3605) + 
     &   MTy*(-2*Cd(3658)*UStop2(1,1) + 
     &      8*Cd(3654)*MStop2(3)*UStop2(1,2))

        Cd(3747) = M_3C*MTy*Cd(3730)*UCStop(1,3)*UStop2(1,1) + 
     &   Cd(3746)*UStopC(1,1) - Cd(3745)*UUStopC(1,3)

        Cd(4006) = Cd(4005)*MStop2(3) + 
     &   4*MTy2*Cd(3928)*UStop2(1,1)**2*UStopC(1,1) - 
     &   Cd(3965)*UUStopC(1,1)

        Cd(4007) = 4*Cd(4006) - 4*AtC**2*Cd(3818)*UCStop(1,3) - 
     &   2*At**2*Cd(3907)*MStop2(3)*UCStopC(1,3)

        Cd(4041) = -(MUE**2*Cd(4040)*UCStop(1,3)) + 
     &   MUEC**2*Cd(4027)*UCStopC(1,3) + 
     &   2*MUE2*Cd(4037)*UStop2(1,1)

        Cd(4136) = -(MUE**2*Cd(4135)*UCStop(1,3)) - 
     &   MUEC**2*Cd(4127)*UCStopC(1,3) + 
     &   2*MUE2*Cd(4132)*UStop2(1,1)

        Cd(4541) = AtC**2*Cd(4300)*UCStop(1,3) + 
     &   At*(4*Cd(4392) + At*Cd(4322)*UCStopC(1,3))

        Cd(4542) = (-4*Cd(4539) + 2*AtC*MTy*Cd(4540) - 
     &     MTy*Cd(4541)*MStop2(3))/MTy

        Cd(4559) = Cd(424)*Cd(4278) + MUE*Cd(4557) - 
     &   MUEC*Cd(4558)*MStop2(3)

	Cd(4560) = SA2*Cd(4302) + CA2*Cd(4542) + S2A*Cd(4559)

        Cd(4878) = -(AtC**2*Cd(4636)*UCStop(1,3)) + 
     &   At*(4*Cd(4728) - At*Cd(4658)*UCStopC(1,3))

        Cd(4916) = Cd(424)*Cd(4902) + MUE*Cd(4914) + 
     &   MUEC*Cd(4915)*MStop2(3)

        Cd(5752) = -(Cd(5746)*UStop(1,1)) + 
     &   2*M_3*MTy*Cd(5751)*UStop(1,2) + 
     &   4*MTy2*Cd(5740)*UStop2(1,2)

        Cd(5984) = 2*MGl*MTy*Cd(5974) + Cd(5979)*UStop(1,1) - 
     &   M_3*MTy*Cd(5983)*UStop(1,2)

        Cd(6256) = MGl*MTy*Cd(6247)*UStop(1,2) - 
     &   M_3*MTy*Cd(6255)*UUStop(1,2) + Cd(5979)*UUStop(1,3)

        Cd(6983) = Cd(6846)*UCStop(1,3) + 
     &   M_3*Cd(50)*Cd(6744)*UStop2(1,1)**2 + 
     &   Cd(6982)*UUStopC(1,2)

        Cd(7751) = Cd(7749)*UStop2(1,1) - Cd(7696)*UStop2(1,2) - 
     &   Cd(302)*Cd(7750)*UStopC(1,1)

        Cd(7965) = 2*M_3*Cd(7952)*MStop2(3) + 
     &   MTy*(-(Cd(7959)*UCStop(1,3)) + 2*Cd(7964)*UUStopC(1,2))

        Cd(8241) = Cd(6678) + 
     &   2*Cd(8240)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(8272) = Cd(6844) + 2*Cd(8271)

	Cd(8290) = -Cd(6967) + 2*Cd(8289)

        Cd(8375) = -2*M_3C*MTy*Cd(8347)*MStop2(3)*
     &    UCStop(1,3) + Cd(8352)*UStop2(1,2)**2 - 
     &   Cd(8374)*UStopC(1,1)

        Cd(8383) = MUE**2*Cd(8382)*UCStop(1,3) - 
     &   MUEC**2*Cd(8343)*MStop2(3)*UCStopC(1,3) + 
     &   2*MUE2*Cd(8375)*UStop2(1,1)

        Cd(8406) = Cd(8352)*UStop2(1,1)**2 + 
     &   Cd(8398)*MStop2(3)*UStop2(1,2) - Cd(8405)*UStopC(1,1)

        Cd(8470) = Cd(7543) + 
     &   2*M_3C*MTy*Cd(8238)*MStop2(3)*UCStop(1,3)

        Cd(8597) = -(MTy2*Cd(8596)) + Cd(8583)*UStop2(1,1) + 
     &   MTy*Cd(8595)*MStop2(3)**2*UCStop(1,3)*UStop2(1,2)

        Cd(8598) = 4*AtC*Cd(8597) + 
     &   At*(4*Cd(8544) - 2*AtC*Cd(8406)*UStop2(1,1))

        Cd(8647) = Cd(8644)*UCStop(1,3) + 
     &   2*M_3*MTy2*Cd(50)*Cd(8257)*MStop2(3)*UStop2(1,1)**2 - 
     &   Cd(8646)*UUStopC(1,2)

        Cd(8653) = AtC*MUEC*Cd(8652) + 
     &   At*(MUEC*Cd(8604)*UCStopC(1,3) + 
     &      MUE*Cd(8240)*UStop2(1,1)*UStop2(1,2))

        Cd(8702) = MUE2*Cd(8687) + 
     &   2*MUE**2*Cd(8701)*UCStop(1,3) + 
     &   MUEC**2*Cd(4238)*MStop2(3)*UCStopC(1,3)

        Cd(8934) = MUE2*Cd(8918) - 
     &   2*MUE**2*Cd(8933)*UCStop(1,3) - 
     &   MUEC**2*Cd(4586)*MStop2(3)*UCStopC(1,3)

        Cd(9167) = (-(Cd(47)**2*Cd(9166)) + 
     &     Cd(50)**2*(Cd(9154) + 
     &        MTy*Cd(47)**2*Cd(9142)*UStop2(1,1)))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(910) = S2A*Cd(868) + 2*CA2*Cd(909) + 
     &   8*SA2*Cd(13)**2*Cd(833)*Cr(75)

        Cd(9871) = 2*SA2*Cd(909) - S2A*Cd(9870) + 
     &   8*CA2*Cd(13)**2*Cd(833)*Cr(75)

        Cd(10257) = (Cd(50)**2*Cd(10190) - Cd(10256) + 
     &     Cd(179)*Cd(2681)*(MGl2*Cd(10191) - Cd(10114)*MStop2(2)))
     &    /(Cd(50)**2*MStop2(2)**2)

        Cd(10570) = 2*MTy2**2*Cd(164)*Cd(10532) + 
     &   2*MTy*Cd(10569) - 
     &   (MGl2*Cd(47)*Cd(10537) + MTy2*Cd(10545))/MStop2(3)

        Cd(10664) = 4*MTy*Cd(10639) + MTy2*Cd(10663) + 
     &   (MGl2*Cd(10625) + 4*MTy2**2*Cd(10620)*MStop2(3))/
     &    (MStop2(2)**2*MStop2(3))

        Cd(12217) = MUEC*Cd(12215) + 
     &   (MUE*Cd(47)**2*Cd(12216) - 
     &      MTy*MUE*Cd(10985)*MGlpTmSt2(1)*MStop2(3)*
     &       UCStop(1,3)**2*UStopC(1,2) - Cd(12209)*UUStopC(1,2))/
     &    (Cd(47)**2*MStop2(3))

        Cd(12415) = 4*Cd(12410) + 
     &   Cd(10202)*Cd(12371)*MStop2(3) + 
     &   AtC*(2*Cd(12414) + AtC*Cd(12404)*UCStop(1,3))

        Cd(12416) = -(MGl2*Cd(12397)) - MTy2*Cd(12403) + 
     &   (2*MTy*Cd(12415) + 4*Cd(47)*Cd(11730)*MStop2(1)**3)/
     &    (MStop2(1)**2*MStop2(3))

        Cd(12562) = (Cd(12558) + Cd(12561)*MStop2(2))/
     &    MStop2(2)**2 + (4*Cd(50)*
     &      (Cd(12557) + Cd(12556)*MStop2(3)*UCStop(1,3)))/
     &    MStop2(3)

        Cd(12568) = (Cd(12567) + Cd(12565)*MStop2(2))/
     &    MStop2(2)**2 + (4*Cd(50)*
     &      (Cd(12563) + Cd(12556)*MStop2(3)*UCStopC(1,3)))/
     &    MStop2(3)

        Cd(12569) = MUE*Cd(12562) + MUEC*Cd(12568) - 
     &   (Cd(13)*(MTy2*Cd(12554) + MGl2*Cd(12555)))/
     &    (MStop2(2)**2*MStop2(3))

        Cd(12650) = (Cd(13)*Cd(47)**2*Cd(12627) + 
     &     Cd(47)**2*Cd(50)**2*
     &      (2*M_3*Cd(12637) - 2*M_3C*Cd(12645)) + 
     &     Cd(50)**2*Cd(12649))/(Cd(47)**2*Cd(50)**2)

        Cd(13946) = MGl2*Cd(13943)*MSq2Diff(1,2)*MStop2(2) + 
     &   (MUEC*TB*Cd(13940) - Cd(13945) + 
     &      At*MUE*TB*Cr(56)*MSq2Diff(1,2)*MStop2(2)**2)*
     &    UStop2(1,1)

        Cd(14355) = Cd(14298) + 
     &   Cd(118)*Cd(14265)*Cd(14266)*MSq2Diff(1,2)**2*
     &    T134(0.D0,MGl2,MStop2(1)) + 
     &   Cd(138)*Cd(14299)*Cd(14300)*MSq2Diff(1,1)**2*
     &    T134(0.D0,MGl2,MStop2(2)) + 
     &   Cd(14353)*T134(MGl2,MTy2,MSbot2(1))

        Cd(1488) = 4*Cd(1466)*MStop2(3) - 
     &   AtC**2*Cd(1472)*UCStop(1,3)**2 + 
     &   At*(2*Cd(1487) + At*Cd(1486)*UCStopC(1,3)**2)

        Cd(2678) = 2*MUEC*Cd(2676) + MUE*Cd(2677) - 
     &   At*MUEC*Cd(2610)*UCStopC(1,3) + 
     &   Cd(424)*Cd(2616)*UStop2(1,1)

        Cd(2754) = 2*MUEC*Cd(2752) + MUE*Cd(2753) + 
     &   At*MUEC*Cd(2693)*UCStopC(1,3) + 
     &   Cd(424)*Cd(2697)*UStop2(1,1)

        Cd(3895) = Cd(3875)*UStop2(1,1)**2 + 
     &   MStop2(2)*(2*MTy*Cd(3884) - Cd(3894) + 
     &      8*M_3*MTy2*Cd(3119)*MStop2(3)*UStop2(1,1)**3)

        Cd(4008) = 4*AtC*Cd(3895) + 
     &   At*(2*AtC*Cd(3795) + 4*Cd(3914)) - Cd(4007)*MStop2(2)

        Cd(4112) = 2*MUEC*Cd(4110) + MUE*Cd(4111) + 
     &   At*MUEC*Cd(4045)*UCStopC(1,3) + 
     &   Cd(424)*Cd(4037)*UStop2(1,1)

        Cd(4213) = 2*MUEC*Cd(4205) + MUE*Cd(4212) - 
     &   At*MUEC*Cd(4143)*UCStopC(1,3) + 
     &   Cd(424)*Cd(4132)*UStop2(1,1)

        Cd(4879) = (2*At*AtC*MTy*Cd(4624) - 4*AtC*MTy*Cd(4772) + 
     &     4*Cd(4877) - MTy*Cd(4878)*MStop2(3))/MTy

        Cd(5223) = 2*(3*Cd(47)**2*Cd(5013) - 
     &      2*M_3C*MTy*Cd(5027)*UCStop(1,3))*UStop2(1,1) + 
     &   Cd(5222)*UStopC(1,1) - 
     &   2*Cd(47)*Cd(50)*Cd(5220)*UStop2(1,2)*UStopC(1,2)

        Cd(7168) = (3*Cd(47)**2*Cd(5013) - 
     &      2*M_3C*MTy*Cd(7058)*UCStop(1,3))*UStop2(1,1) + 
     &   Cd(7167)*UStopC(1,1) - 
     &   Cd(47)*Cd(50)*Cd(5220)*UStop2(1,2)*UStopC(1,2)

        Cd(7554) = 2*Cd(7548) - Cd(7553)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(5751)*UStop(1,2) + 
     &   4*MTy2*Cd(5740)*UStop2(1,2)

        Cd(10047) = (2*(4*Cd(10014)*Cd(10036) + 
     &       (MStop2(3)*(2*Cd(47)*Cd(50)*Cd(10046) - 
     &            MStop2(3)*
     &             (Cd(10041) + 
     &               (-(Cd(50)**2*Cd(10012)*MStop2(1)) + 
     &                  Cd(47)**2*Cd(10013)*MStop2(2))*UCStopC(1,3)
     &               )))/(Cd(47)**2*Cd(50)**2)))/MStop2(3)**2

        Cd(10186) = -4*Cd(10161)*MStop2(1) - 
     &   (2*Cd(47)*Cd(50)*Cd(10185) + 
     &      (-4*MGl2*Cd(10157) + Cd(10168) + 
     &         4*Cd(47)**2*Cd(10159)*MStop2(2))*MStop2(3))/
     &    (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(10665) = (2*Cd(10191))/MStop2(2)**2 + 
     &   (Cd(10664) + (2*Cd(10619))/MStop2(2) + 
     &      (4*Cd(50)*(Cd(10608)*MStop2(1) + Cd(10607)*MStop2(3)))/
     &       MStop2(3))/Cd(50)**2

        Cd(12525) = 4*Cd(50)*Cd(11741) + 
     &   (4*MGl2*Cd(10193)*Cd(12483) - Cd(12520) + 
     &      4*Cd(12524)*MStop2(2))/MStop2(2)**2 - 
     &   (4*Cd(50)*Cd(12521))/MStop2(3)

        Cd(12589) = (-2*Cd(10153)*Cd(12045) - 4*MTy2*Cd(12588) + 
     &     4*Cd(585)*(MGl2*Cd(12045) - 
     &        Cd(50)**2*MGlpTmSt2(1)*MStop2(1) + 
     &        Cd(47)**2*MGlpTmSt2(2)*MStop2(2)))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(13957) = MTy*Cd(13956) + 
     &   MTy2*(Cd(13946) - MTy2*(MTy2*Cd(13949) + Cd(13953)) + 
     &      MTy2**3*Cd(13737)*MSq2Diff(1,2))

        Cd(14094) = (MTy*MTy2*TB*Cd(14093) - 
     &     MSq2Diff(1,1)*MStop2(3)*
     &      (2*MTy2*TB*Cd(9950)*Cd(14012)*Cr(5) + 
     &        At*MTy*MTy2*MUE*TB2*Cd(47)*Cd(14076)*
     &         UStop2(1,2)**2 + 
     &        MTy*(MTy2*Cd(14089) + 
     &           TB*Cd(9950)*Cd(14048)*MSq2Diff(1,1)*MStop2(1)*
     &            UStop2(1,2))))/
     &   (MTy*MTy2*TB*TB2*Cd(47)*Cd(9950)*MSq2Diff(1,1)**2*
     &     MStop2(1)*MStop2(3))

        Cd(14182) = -((MUE*Cd(14152) + MUEC*Cd(14154))/
     &      (TB*Cd(9950))) + (MGl2*Cd(14162))/(TB*TB2*Cd(47)) + 
     &   MTy*Cd(14174) + MTy2*Cd(14181)

        Cd(2672) = 4*MTy2*Cd(2658) + 
     &   AtC*(4*Cd(2670) - AtC*Cd(2624)*UCStop(1,3)) + 
     &   At*(2*Cd(2671) - At*Cd(2633)*UCStopC(1,3))

        Cd(2748) = -4*MTy2*Cd(2734) + 
     &   AtC*(4*Cd(2746) - AtC*Cd(2702)*UCStop(1,3)) + 
     &   At*(2*Cd(2747) + At*Cd(2709)*UCStopC(1,3))

        Cd(4105) = AtC*(4*Cd(4077) - AtC*Cd(4040)*UCStop(1,3)) + 
     &   At*(4*Cd(4063) + At*Cd(4045)*UCStopC(1,3) + 
     &      2*AtC*Cd(4037)*UStop2(1,1))

        Cd(9444) = -2*M_3C*MTy*Cd(9397) + 2*Cd(9443) - 
     &   Cd(9377)*MStop2(3) + 
     &   4*MTy2*Cd(9279)*UStop2(1,1)**2*UStopC(1,1) + 
     &   Cd(9362)*UUStopC(1,1)

        Cd(9539) = -2*M_3C*MTy*Cd(9507) + 2*Cd(9538) + 
     &   Cd(9495)*MStop2(3)**2 - 
     &   4*MTy2*Cd(9468)*MStop2(3)*UStop2(1,1)**2*UStopC(1,1) - 
     &   Cd(9484)*UUStopC(1,1)

        Cd(11939) = -(MTy2*Cd(11921)) + 
     &   (MUE*Cd(50)*Cd(11901) + MUEC*Cd(50)*Cd(11906) + 
     &      2*MTy*Cd(11938) + 
     &      MGl2*(Cd(11925) - 2*MGl2*Cd(179)*Cd(11783)*MStop2(3)))/
     &    (MStop2(2)*MStop2(3))

        Cd(12068) = (4*Cd(47)*Cd(50)*Cd(12053) - 
     &     4*Cd(50)**2*Cd(585)*Cd(12044)*MStop2(1) + 
     &     2*(2*MGl2*Cd(585)*Cd(12043) + 2*MTy2*Cd(12049) + 
     &        Cd(12067) + 2*Cd(47)**2*Cd(585)*MGlpTmSt2(2)*
     &         MStop2(2))*MStop2(3))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(12519) = (-16*Cd(50)*Cd(11866) + 4*MTy2*Cd(12510) - 
     &     (-4*MGl2*Cd(585)*Cd(12483) + 2*Cd(12513) + 
     &        Cd(50)**2*Cd(12518)*MStop2(2))/MStop2(2)**2 + 
     &     (8*Cd(50)*Cd(10630))/MStop2(3))/Cd(50)**2

        Cd(12570) = -(MUE*Cd(12537)) - MUEC*Cd(12548) + 
     &   (2*MTy*Cd(12569))/Cd(50)**2 + 
     &   (MTy2*Cd(12553) + 
     &      (MGl2*MTy*Cd(12541) + 
     &         4*Cd(13)*Cd(12538)*MStop2(2)*MStop2(3))/
     &       (MTy*MStop2(2)**2*MStop2(3)))/Cd(50)**2

        Cd(14183) = (MTy2*Cd(47)*Cd(9950)*
     &      (MTy*TB*TB2*Cd(14182) + 
     &        (-(MTy*TB*Cd(14144)) + Cd(118)*Cd(14148))*
     &         MSq2Diff(1,1)) - MTy*MTy2*TB*Cd(14146)*MStop2(1) + 
     &     MTy*TB*Cd(9950)*
     &      (MTy2**3*Cd(9969)*Cd(14142) + 
     &        Cd(47)*Cd(118)*MSq2Diff(1,1)**2)*UStop2(1,1))/
     &   (MTy*MTy2*TB*TB2*Cd(47)*Cd(9950)*MSq2Diff(1,1)**2)

        Cd(10156) = (-4*Cd(10014)*Cd(10138) + 
     &     (MStop2(3)*(2*Cd(47)*Cd(50)*Cd(10155) + 
     &          (Cd(10154) + 
     &             2*(-(MTy2*Cd(10145)) + 
     &                Cd(50)**2*Cd(585)*
     &                 (MGl2*Cd(10133) - Cd(10134)*MStop2(1)) + 
     &                Cd(47)**2*Cd(585)*Cd(10013)*MStop2(2)))*
     &           MStop2(3)))/(Cd(47)**2*Cd(50)**2))/MStop2(3)**2

        Cd(323) = -(MTy2*Cd(317)) - MTy*Cd(322) - 
     &   Cd(314)*MStop2(3) + 
     &   M_3*(MTy2**2*(2*MTy2 - Cd(313)) + 
     &      MGl2*MStop2(3)*(-MGl2 + MStop2(4)))*UCStopC(1,3)

        Cd(12417) = -(Cd(12373)/MTy2) - 
     &   (2*Cd(10191))/MStop2(1)**2 + 
     &   (4*Cd(12383) - Cd(12416) + 
     &      (4*MTy*MTy2**2*Cd(12371) + 
     &         (-Cd(12381) + MTy*Cd(47)**2*Cd(12394))*MStop2(1))/
     &       (MTy*MStop2(1)**2) + (4*Cd(47)*Cd(11722))/MStop2(3))/
     &    Cd(47)**2

        Cd(13922) = (MTy*MTy2*Cd(13921) + 
     &     TB*(2*MTy*MTy2**3*Cd(13900) - 
     &        MSq2Diff(1,2)*
     &         (MTy*MTy2*Cd(50)*Cd(13898)*MStop2(1) + 
     &           MSq2Diff(1,1)*
     &            (2*MTy2*Cd(11964)*Cr(5) + 
     &              MTy*(-Cd(11961) + 
     &                 MTy2*
     &                  (Cd(13895) + Cd(47)*Cd(13893)*MStop2(2)))))
     &          - 2*MTy*MTy2*Cd(47)*Cd(13841)*MSbot2(1)*
     &         MSq2Diff(1,1)*UStop2(1,1)))/
     &   (MTy*MTy2*TB*TB2*Cd(47)*Cd(50)*MSq2Diff(1,1)*
     &     MSq2Diff(1,2))

        Cd(13924) = 1/2.D0*
     &   (-(MTy*MTy2*Cd(13921)) + 
     &      TB*(-2*MTy*MTy2**3*Cd(13900) + 
     &         MSq2Diff(1,2)*
     &          (MTy*MTy2*Cd(50)*Cd(13898)*MStop2(1) + 
     &            MSq2Diff(1,1)*
     &             (2*MTy2*Cd(11964)*Cr(5) + 
     &               MTy*(Cd(12662) + 
     &                  MTy2*
     &                   (Cd(13923) + Cd(47)*Cd(13893)*MStop2(2))))
     &            ) + 2*MTy*MTy2*Cd(47)*Cd(13841)*MSbot2(1)*
     &          MSq2Diff(1,1)*UStop2(1,1)))/
     &    (MTy*MTy2*TB*TB2*Cd(47)*Cd(50)*MSq2Diff(1,1)*
     &      MSq2Diff(1,2))

        Cd(14003) = MTy2*
     &    (Cd(13978) + MTy2*
     &       (MTy2*(-(MTy2*Cd(13734)*Cd(13737)) + Cd(13982)) + 
     &         Cd(13992))) + MTy*Cd(14002) + 
     &   MSq2Diff(1,2)**2*(MGl2 - MStop2(2))**2*MStop2(2)*
     &    UStop2(1,2)

        Cd(7752) = -(MTy*Cd(7751)) + 
     &   2*M_3*(Cd(44)*Cd(48)*Cd(7700) + 
     &      MGl2*(-Cd(7708) + 
     &         MGl2*(Cd(7703) + 
     &            MGl2*(-Cd(7699) + 
     &               MGl2*(Cd(7697) + MGl2*(MGl2 - Cd(7709)))))))*
     &    MStop2(3)*UCStopC(1,3)

        Cd(12069) = Cd(11961)/(MTy2*Cd(47)*Cd(50)) - 
     &   Cd(11970)/MTy + MGl2*Cd(11994) + 
     &   Cd(12008)/(Cd(47)**2*Cd(50)**2) + MTy2*Cd(12041) - 
     &   MTy*Cd(12068) + (MStop2(1)*
     &      (-Cd(11982) + 4*MGlpTmSt2(1)*MStop2(1)))/Cd(47)**2 - 
     &   (4*MTy2**2*Cd(11963) + Cd(11974)*MStop2(2))/Cd(50)**2 + 
     &   (4*Cd(11964)*Cd(11965))/(Cd(47)*Cd(50)*MStop2(3))

        Cd(12526) = Cd(12485)/MTy2 + MTy*Cd(12519) + 
     &   (2*Cd(10191))/MStop2(2)**2 + 
     &   (-4*MTy*MTy2**2*Cd(12486) + 
     &      MStop2(2)*(MTy*Cd(12509) + 
     &         4*MTy*Cd(50)*Cd(11846)*MStop2(2) - 
     &         (Cd(12494) + 
     &            MTy*MStop2(2)*
     &             (4*Cd(12496) + MGl2*Cd(12499) + 
     &               MTy2*Cd(12525) - 8*Cd(50)*MStop2(2)))*
     &          MStop2(3)))/(MTy*Cd(50)**2*MStop2(2)**2*MStop2(3))

        Cd(8778) = 2*M_3C*Cd(4401)*MStop2(3)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(-Cd(8767) + 2*Cd(8758)*UStop2(1,1)) + 
     &   4*MTy*Cd(8777)*UCStop(1,3)*UStop2(1,2) - 
     &   4*M_3*(MTy2*Cd(8769) + 
     &      MStop2(3)*(-((MGl2 - MStop2(1))**3*MStop2(1)) + 
     &         MTy2**2*(Cd(4421) + 
     &            MTy2*(-Cd(4422) + MTy2*UStop2(1,2)**2))))

	Cd(9012) = Opt(18) - Cd(8995)*UStop2(1,1)**2

	Cd(9026) = Opt(18) - Cd(9025)*UStop2(1,1)**2

        Cd(12134) = ((MGl2**2 + MTy2**2)*Cd(12121))/Cd(50)**2 + 
     &   (2*MTy*Cd(47)**2*Cd(12106) + 
     &      (2*MGl2*Cd(12128) + 2*MTy2*Cd(12132) + 
     &         Cd(47)**2*Cd(50)**2*Cd(12133) - 
     &         MUE**2*Cd(12116)*UCStop(1,3)**2 + 
     &         MUEC**2*Cd(12119)*UCStopC(1,3)**2 + 
     &         2*MUE2*Cd(47)*
     &          (-(Cd(50)*Cd(12112)*MStop2(1)) + 
     &            Cd(47)*MStop2(3)*
     &             (Cd(50)*Cd(12110) + 
     &               Cd(12108)*MStop2(2)**2*UStop2(1,1))) + 
     &         MUEC*Cd(50)**2*Cd(12074)*UUStopC(1,1))/
     &       (Cd(50)**2*MStop2(3)) + 
     &      MUE*Cd(6458)*Cd(12075)*UStop2(1,1)*UUStopC(1,2))/
     &    Cd(47)**2

        Cd(12680) = Cd(12678)/Cd(50) + 
     &   ((2*MGl2*Cd(12128) + 2*MTy2*Cd(12132) - 
     &         MUE**2*Cd(12116)*UCStop(1,3)**2 + 
     &         MUEC**2*Cd(12119)*UCStopC(1,3)**2 + 
     &         Cd(47)*(Cd(47)*Cd(50)**2*Cd(12133) - 
     &            2*MUE2*Cd(50)*Cd(12112)*MStop2(1) + 
     &            Cd(47)*MStop2(3)*
     &             (MGl2**2*Cd(12121) + MTy2**2*Cd(12121) + 
     &               MTy*Cd(50)**2*Cd(12677) - 
     &               2*MUE2*Cd(12679)*MStop2(2)**2*UStop2(1,1))) + 
     &         MUEC*Cd(50)**2*Cd(12074)*UUStopC(1,1))/
     &       (Cd(50)**2*MStop2(3)) + 
     &      MUE*Cd(6458)*Cd(12075)*UStop2(1,1)*UUStopC(1,2))/
     &    Cd(47)**2

	end
