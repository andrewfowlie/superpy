* TLhr_abbr3_4.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr3_4
	implicit none

#include "TLhr.h"


        Co(306) = (MUE*Co(298) + MUEC*Co(305))/TB + 
     &   MUE**2*Co(293)*UCStop(1,3) + MUEC**2*Co(301)*UCStopC(1,3)

        Co(356) = (MGl2*TB*Cd(47)*Co(330) + 4*MTy*Co(355) + 
     &     MTy2*(Co(335) + 2*MTy2*Cd(164)*Co(327)*MStop2(3)))/
     &   (TB*MStop2(3))

        Co(428) = (MUE*Co(426) + 
     &     MUE**2*TB*Cd(12473)*UCStop(1,3)**2 - 
     &     MUEC*(Co(427) + MUEC*TB*Ce(1104)*UCStopC(1,3)**2))/TB

        Co(518) = (MUE*Co(510) + 
     &     MUE**2*TB*Co(507)*UCStop(1,3)**2 - 
     &     MUEC*(Co(517) + MUEC*TB*Co(512)*UCStopC(1,3)**2))/TB

        Co(528) = (MUE*Co(524) + MUEC*Co(527))/TB + 
     &   MUE**2*Ce(1364)*UCStop(1,3) + MUEC**2*Co(525)*UCStopC(1,3)

        Co(545) = (-2*Cd(50)*Co(542) + 
     &     2*Cd(9948)*Co(544)*MStop2(3) - 
     &     Co(541)*MStop2(3)*UCStop(1,3)**2 + 
     &     2*At*Co(539)*MStop2(3)*UStop2(1,1))/MStop2(3)

        Co(546) = Co(536)/MStop2(2)**2 + 
     &   (-(Co(545)/MStop2(2)) + 
     &      (2*Cd(50)*(Co(538)*MStop2(1) + Co(537)*MStop2(3)))/
     &       MStop2(3))/Cd(50)**2

        Co(583) = (MUE*Co(579) + MUEC*Co(582))/TB - 
     &   MUE**2*Ce(1312)*UCStop(1,3) + 
     &   MUEC**2*Cd(10766)*UCStopC(1,3)

        Co(72) = AtC*(MUEC*Co(70) - 
     &      2*MUE*Cd(639)*MStop2(3)*UCStop(1,3)**2*UStop2(1,1)) + 
     &   At*(MUE*Cd(636) + 
     &      4*MUEC*Cd(642)*UCStopC(1,3)**2*UStop2(1,1))

        Co(815) = 2*(Co(807) - 2*AtC*Co(814) + 
     &     2*Co(810)*MStop2(1)*MStop2(3) - 
     &     2*Co(749)*MStop2(3)**2*UUStopC(1,1))

        Co(816) = Co(802) - Co(803) + 
     &   2*Co(791)*MStop2(3)*UCStop(1,3)*UStop2(1,2)*UStopC(1,2) - 
     &   Co(815)*UUStopC(1,2)

        Co(886) = -(MUEC*Cd(11296)*Cd(12701)) + Co(883) - 
     &   Co(885)*MStop2(2) - 
     &   4*Co(837)*MStop2(3)*UCStop(1,3)*UStop2(1,1)*UStopC(2,2)

        Co(1019) = AtC*Co(1006) + 4*Co(1018)*MStop2(3) - 
     &   Co(1015)*UUStopC(1,1) + Co(971)*UUStopC(1,2) - 
     &   Co(980)*UUStopC(1,3)

        Co(1022) = Co(1005) + 
     &   Co(1020)*UStopC(2,2)*UUStopC(1,1) + 
     &   2*Co(958)*MStop2(3)*UUStopC(1,2) - Co(1021)*UUStopC(1,3)

        Co(1157) = 1/(16.D0*Pi**3)*
     &    (16*Pi**3*(SB*SW2*Cd(47)**2*Cd(50)**2*Co(1156)*
     &          MStop2(3)**2) + 
     &       cI*(Alfas2L*Cd(585)*MStop2(2)*
     &          (Cd(47)*Cd(50)*
     &             (4*Cd(47)*Cd(50)*Cd(10014)*Cd(10016) + 
     &               Co(1153)*MStop2(3)) - 
     &            Co(1129)*MStop2(3)**2*UStop2(1,1)*UStop2(1,2))))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(2)*MStop2(3)**2) + 
     &   Co(1155)/MStop2(1)

        Co(1269) = -2*MUEC*TB*Cd(10998) + AtC*Co(1260) - 
     &   2*At*Co(1268) + 4*
     &    (MTy*MUE*TB*Co(1262) + 
     &      AtC**2*Cd(11105)*MStop2(1)**2*UCStop(1,3))*UUStop(1,1)

        Co(1671) = AtC*Co(1647) + 4*Co(1670)*MStop2(3) + 
     &   Co(1667)*UUStopC(1,1) + Co(1668)*UUStopC(1,2) + 
     &   Co(1607)*UUStopC(1,3)

        Co(2217) = 1/(16.D0*Pi**3)*
     &    (16*Pi**3*(SB*SW2*Cd(47)**2*Cd(50)**2*Co(2216)*
     &          MStop2(3)**2) + 
     &       cI*(Alfas2L*Cd(585)*MStop2(2)*
     &          (Cd(47)*Cd(50)*
     &             (4*Cd(47)*Cd(50)*Cd(10014)*Cd(10016) + 
     &               Co(2210)*MStop2(3)) - 
     &            Co(1129)*MStop2(3)**2*UStop2(1,1)*UStop2(1,2))))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(2)*MStop2(3)**2) + 
     &   Co(2214)/MStop2(1)

        Co(2271) = 1/(16.D0*Pi**3)*
     &    (16*Pi**3*(SB*SW2*Cd(47)**2*Cd(50)**2*Co(2270)*
     &          MStop2(3)**2) - 
     &       cI*(Alfas2L*Cd(585)*MStop2(2)*
     &          (Cd(47)*Cd(50)*
     &             (4*Cd(47)*Cd(50)*Cd(10014)*Cd(10016) - 
     &               Co(2196)*MStop2(3)) - 
     &            Co(1129)*MStop2(3)**2*UStop2(1,1)*UStop2(1,2))))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(2)*MStop2(3)**2) + 
     &   Co(2269)/MStop2(1)

        Co(325) = AtC*Ce(618) + MGlmT2*Cd(9948)*Co(322) + 
     &   2*(Co(324) - At*Cd(9948)*MStop2(1)**2)*UCStopC(1,3)**2 - 
     &   2*M_3C*MStop2(1)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Co(367) = (Co(360) - Co(366)*MStop2(1))/MStop2(1)**2 + 
     &   (2*MGlmT2*Co(363) - 2*Co(361)*MStop2(1) + 
     &      4*Co(362)*MStop2(3))/(Cd(47)*MStop2(3))

        Co(376) = (Co(368) - Co(375)*MStop2(1))/MStop2(1)**2 + 
     &   (2*MGlmT2*Co(369) - 2*Co(371)*MStop2(1) + 
     &      4*Co(370)*MStop2(3))/(Cd(47)*MStop2(3))

        Co(455) = (MGl2*Cd(152)*Co(454)*MStop2(3) + 
     &     Cd(47)*Cd(164)*(MUE*Co(452) + 
     &        MUE**2*TB*Cd(10088)*UCStop(1,3) + 
     &        MUEC*(Co(453) + MUEC*TB*Cd(10042)*UCStopC(1,3))))/
     &   (TB*MStop2(3))

        Co(483) = (MGl2*Cd(179)*Co(454)*MStop2(3) + 
     &     Cd(50)*Cd(251)*(-(MUE*Co(481)) + 
     &        MUE**2*TB*Cd(10021)*UCStop(1,3) + 
     &        MUEC*(Co(482) + MUEC*TB*Cd(10095)*UCStopC(1,3))))/
     &   (TB*MStop2(3))

        Co(1367) = 2*MUEC*TB*Cd(10998) - AtC*Co(1356) - 
     &   MTy*Co(1362) + At*MUE*TB*Co(1366) + 
     &   MStop2(1)**2*(-4*AtC**2*Cd(11105)*UCStop(1,3)*
     &       UUStop(1,1) + 
     &      4*At**2*Cd(12114)*UCStopC(1,3)*UUStop(1,2))

        Co(166) = At*MUE*Cd(423)*Co(138) - 
     &   2*AtC*MUEC*Cd(577)*Co(138) + Opt(32)

        Co(170) = At*MUE*Cd(435)*Co(138) - 
     &   AtC*MUEC*Cd(435)*Co(138) + Opt(32)

        Co(1831) = MUE*Co(1824) + MUEC*Co(1830) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(1825) + 
     &         Cd(9950)*Cd(10103)*Co(163)*Co(1821)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(1828)/MStop2(1) + Co(1826)/MStop2(2)

        Co(1839) = MUE*Co(1833) + MUEC*Co(1838) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(1825) + 
     &         Cd(9950)*Cd(10103)*Co(163)*Co(1821)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(1835)/MStop2(1) + Co(1834)/MStop2(2)

        Co(1956) = At*MUE*Cd(435)*Co(1929) + 
     &   AtC*MUEC*Co(159)*Co(1929) + Opt(33)

        Co(1958) = -(AtC*MUEC*Cd(423)*Co(1929)) + 
     &   2*At*MUE*Cd(577)*Co(1929) + Opt(33)

        Co(2416) = MUE*Co(2407) + MUEC*Co(2415) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(2408) + 
     &         Cd(9950)*Cd(10103)*Co(163)*Co(1821)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2412)/MStop2(1) + Co(2410)/MStop2(2)

        Co(2425) = MUE*Co(2419) + MUEC*Co(2424) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(2408) + 
     &         Cd(9950)*Cd(10103)*Co(163)*Co(1821)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2421)/MStop2(1) + Co(2420)/MStop2(2)

        Co(307) = -(MGlpT2*Co(280)) + 
     &   (MUE*Co(272) + MUEC*Co(292))/TB + 2*MTy*Co(306) + 
     &   MUE**2*Co(260)*UCStop(1,3)**2 - 
     &   MUEC**2*Cd(10079)*UCStopC(1,3)**2

        Co(396) = (-(MUE*Co(367)) + 
     &     (MTy2*Co(382) + MGl2*Co(391) + 4*MTy*Co(395) + 
     &        2*MUE**2*TB*Co(359)*UCStop(1,3)**2 + 
     &        MUEC*(Cd(47)**2*Co(376)*MStop2(1)**2*MStop2(3) + 
     &           MUEC*TB*Ce(1033)*UCStopC(1,3)**2))/
     &      (Cd(47)**2*MStop2(1)**2*MStop2(3)))/TB

        Co(584) = (MUE*Co(546) + 
     &     (-(MTy2*Co(555)) + MGl2*Co(576) + 4*MTy*TB*Co(583) + 
     &        MUE**2*TB*Ce(1304)*UCStop(1,3)**2 - 
     &        MUEC*(Cd(50)**2*Co(564)*MStop2(2)**2*MStop2(3) + 
     &           2*MUEC*TB*Co(548)*UCStopC(1,3)**2))/
     &      (Cd(50)**2*MStop2(2)**2*MStop2(3)))/TB

        Co(1024) = AtC*Co(1022) + 
     &   Cr(9)*MStop2(1)*(Co(982)*UStopC(1,1) + 
     &      Co(985)*UStopC(1,2)) - Co(1023)*UUStopC(1,1) + 
     &   MStop2(2)*UStopC(1,1)*
     &    (Co(962)*UCStopC(1,3)*UStop2(1,1) - Co(964)*UUStopC(1,2))
     &     + MUE*(Co(963)*MStop2(2)*UCStop(1,3)*UStop2(1,2)*
     &       UStopC(1,2) - 2*Co(946)*MStop2(3)*UUStopC(1,2)) - 
     &   2*Co(953)*MStop2(3)*UUStopC(1,3)

        Co(2050) = MGl2**2*Cd(10382) - 2*MGl2*Co(2049) + 
     &   MUE*(MTy2 - MStop2(2))**2*UCStop(1,3)**2 + 
     &   MUEC*UStop2(1,1)*(Cd(9948)*Co(653)*UStop2(1,1) - 
     &      2*MTy2*MStop2(2)*UStop2(1,2) + 
     &      (MTy2**2 + MStop2(2)**2)*UStop2(1,2)**2)

        Co(731) = Cd(9948)*Co(730) + 
     &   2*MUEC*UCStopC(1,3)**2*
     &    (MStop2(3) + MStop2(1)*UStop2(1,2)) + 
     &   MUE*UStop2(1,2)*(MStop2(1)*UStop2(1,2) + 
     &      UStop2(1,1)*(MGl2*Cd(12114) + MTy2*UStop2(1,1) + 
     &         (-2*MTy2 + MStop2(2) + 3*MStop2(3))*UStop2(1,2)))

        Co(2079) = MGl2**2*Co(2077) + 2*MGl2*Co(2078) + 
     &   MUE*(MTy2 - MStop2(2))**2*UCStop(1,3)**2 + 
     &   MUEC*UStop2(1,1)*(-2*Cd(9948)*Co(738)*UStop2(1,1) + 
     &      (MTy2**2 + MStop2(2)**2 + 
     &         4*(MTy2 + MStop2(2))**2*UStop2(1,1))*UStop2(1,2) - 
     &      2*MTy2*MStop2(2)*UStop2(1,2)**2)

        Co(2081) = -(MGl2*(MGl2*Co(740) + 2*Co(2080))) + 
     &   MUEC*(MTy2 - MStop2(2))**2*UCStopC(1,3)**2 + 
     &   MUE*UStop2(1,1)*(Cd(9948)*Co(733)*UStop2(1,1) + 
     &      UStop2(1,2)*(-2*MTy2*MStop2(2) - 
     &         4*(MTy2 + MStop2(2))**2*UStop2(1,1) + 
     &         (MTy2**2 + MStop2(2)**2)*UStop2(1,2)))

        Co(737) = -(MGl2*(MGl2*Co(734) + 2*Co(736))) + 
     &   MUE*(MTy2 - MStop2(2))**2*UCStop(1,3)**2 + 
     &   MUEC*UStop2(1,1)*(Cd(9948)*Co(733)*UStop2(1,1) + 
     &      UStop2(1,2)*(-2*MTy2*MStop2(2) - 
     &         4*(MTy2 + MStop2(2))**2*UStop2(1,1) + 
     &         (MTy2**2 + MStop2(2)**2)*UStop2(1,2)))

        Co(310) = Opt(38) + Opt(39) + 
     &   UStop2(1,1)*(AtC*(-(MGl2*MUEC*UStop2(1,2)) + 
     &         MTy2*MUEC*(-(Cd(9948)*UStop2(1,1)) - UStop2(1,2)**2)
     &         ) + At*(MTy2*MUE*UStop2(1,2) + 
     &         MGl2*MUE*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)))

	Co(400) = Opt(38) + Opt(40)

	Co(534) = -Opt(38) + Opt(40)

        Co(2051) = 2*MTy*TB*Cd(50)**2*Co(2048) - AtC*Co(2050) + 
     &   Opt(41) + At*MUEC*MGlpTmSt4(2)*UCStopC(1,3)**2 + 
     &   At*MUE*UStop2(1,1)*
     &    (-2*MGlpT2*Cd(9948)*MStop2(2)*UStop2(1,1) + 
     &      UStop2(1,2)*(MStop2(2)**2 + 
     &         MGlpT2*(MGlpT2 - 2*MStop2(2)*UStop2(1,2))))

        Co(657) = 2*MTy*TB*Cd(50)**2*Co(652) + At*Co(656) + 
     &   Opt(41) - AtC*MUE*MGlpTmSt4(2)*UCStop(1,3)**2 + 
     &   AtC*MUEC*UStop2(1,1)*
     &    (2*MGlpT2*Cd(9948)*MStop2(2)*UStop2(1,1) - 
     &      UStop2(1,2)*(MStop2(2)**2 + 
     &         MGlpT2*(MGlpT2 - 2*MStop2(2)*UStop2(1,2))))

        Co(1199) = 2*Cd(50)*Cd(9947)*Cd(9950)*Cd(10103)*
     &    Cd(12594) - Cd(9948)*Cd(9950)*
     &    (4*Cd(47)*Cd(50)*Ce(411) + 
     &      Cd(10103)*(Cd(50)*(4*Ce(411) - Co(1034)) + 
     &         Cd(47)*Co(1047))) - 
     &   2*Cd(47)*Cd(9949)*Cd(9950)*Cd(10103)*Co(1198) + Opt(42) + 
     &   Opt(43) + 2*(-2*Cd(9950)*
     &       (-(Cd(50)*Cd(10103)) + Cd(47)*(Cd(50) + 2*Cd(10103)))*
     &       Ce(411) + Cd(47)*Cd(50)*Cd(10103)*
     &       (Co(1034) + Co(1037)))*UStop2(1,1) + 
     &   Cd(10103)*(2*Cd(47)*Cd(50)*(Co(1034) + Co(1037)) + 
     &      Cd(9950)*(Cd(50)*Co(1037) + 
     &         Cd(47)*(-8*Ce(411) + 2*Co(1047) + Co(1049))))*
     &    UStop2(1,2)

        Co(2175) = Cd(50)*Cd(9948)*
     &    (Cd(10103)*(4*Cd(47)*Ce(411) + 
     &         Cd(9950)*(Co(1034) + Co(1037))) - 
     &      2*Cd(47)*Cd(9950)*(Co(1047) + Co(1049))) - 
     &   2*Cd(47)*Cd(9949)*Cd(9950)*Cd(10103)*Co(1115) + 
     &   2*Cd(50)*Cd(9947)*Cd(9950)*Cd(10103)*Co(1116) + Opt(42) - 
     &   Opt(43) - 2*(-2*Cd(47)*Cd(50)*Cd(10103)*Ce(411) + 
     &      Cd(9950)*(2*(Cd(47) - 2*Cd(50))*Cd(10103)*Ce(411) + 
     &         Cd(47)*Cd(50)*(Co(1047) + Co(1049))))*UStop2(1,1) - 
     &   Cd(9950)*Cd(10103)*
     &    (2*Cd(50)*(Co(1034) + Co(1037)) + 
     &      Cd(47)*(4*Ce(411) + Co(1047) + Co(1049)))*UStop2(1,2)

        Co(2254) = 2*Cd(50)*Cd(9947)*Cd(9950)*Cd(10103)*
     &    Cd(12594) - Cd(9948)*Cd(9950)*
     &    (4*Cd(47)*Cd(50)*Ce(411) + 
     &      Cd(10103)*(4*Cd(50)*Ce(411) + 
     &         Cd(47)*(Co(1047) + Co(1049)))) - 
     &   2*Cd(47)*Cd(9949)*Cd(9950)*Cd(10103)*Co(1198) + Opt(42) + 
     &   Opt(43) + 2*(-2*Cd(9950)*
     &       (-(Cd(50)*Cd(10103)) + Cd(47)*(Cd(50) + 2*Cd(10103)))*
     &       Ce(411) + Cd(47)*Cd(50)*Cd(10103)*
     &       (Co(1034) + Co(1037)))*UStop2(1,1) + 
     &   Cd(10103)*(2*Cd(47)*Cd(50)*(Co(1034) + Co(1037)) + 
     &      Cd(9950)*(Cd(50)*(Co(1034) + Co(1037)) + 
     &         2*Cd(47)*(-4*Ce(411) + Co(1047) + Co(1049))))*
     &    UStop2(1,2)

        Co(1111) = 2*Cd(47)*Cd(9949)*Co(1110) - 
     &   2*Cd(50)*Cd(9947)*Cd(12594)*Cr(28) + Opt(44) + 
     &   Cd(9948)*(Cd(50)*(-4*Ce(411) + Co(1034)) + 
     &      Cd(47)*(-4*Ce(411) - 3*Co(1047) + Co(1049)))*
     &    UStop2(1,1)**2 + 
     &   (Cd(50)*Co(1037) + 
     &      2*Cd(47)*(4*Ce(411) + Co(1047) - Co(1049)))*
     &    UStop2(1,1)**2*UStop2(1,2) + 
     &   Cd(9948)*(Cd(50)*(-4*Ce(411) + Co(1034)) - 
     &      Cd(47)*Co(1047))*UStop2(1,2)**2 + 
     &   (Cd(50)*Co(1037) + Cd(47)*(4*Ce(411) - Co(1049)))*
     &    UStop2(1,2)**3

        Co(1191) = 2*Cd(47)*Cd(9949)*Co(1190) - 
     &   2*Cd(50)*Cd(9947)*Co(1189)*Cr(28) + Opt(45) + 
     &   Cd(9948)*(Cd(50)*Co(1037) + 
     &      Cd(47)*(12*Ce(411) + Co(1047)))*UStop2(1,1)**2 - 
     &   2*(2*Cd(50)*Ce(411) + Cd(47)*(4*Ce(411) + Co(1047)))*
     &    UStop2(1,1)**2*UStop2(1,2) + 
     &   Cd(9948)*(4*Cd(47)*Ce(411) + Cd(50)*Co(1037))*
     &    UStop2(1,2)**2 - 
     &   (4*Cd(50)*Ce(411) + Cd(47)*Co(1047))*UStop2(1,2)**3

        Co(2167) = 2*Cd(47)*Cd(9949)*Co(1110) - 
     &   2*Cd(50)*Cd(9947)*Cd(12594)*Cr(28) + Opt(45) - 
     &   Cd(9948)*(4*Cd(50)*Ce(411) + 
     &      Cd(47)*(4*Ce(411) + 3*Co(1047)))*UStop2(1,1)**2 + 
     &   (Cd(50)*Co(1037) + Cd(47)*(8*Ce(411) + 2*Co(1047)))*
     &    UStop2(1,1)**2*UStop2(1,2) - 
     &   Cd(9948)*(4*Cd(50)*Ce(411) + Cd(47)*Co(1047))*
     &    UStop2(1,2)**2 + 
     &   (4*Cd(47)*Ce(411) + Cd(50)*Co(1037))*UStop2(1,2)**3

        Co(2247) = 2*Cd(47)*Cd(9949)*Co(1190) - 
     &   2*Cd(50)*Cd(9947)*Co(1189)*Cr(28) + Opt(44) + 
     &   Cd(9948)*(Cd(50)*(Co(1034) + Co(1037)) + 
     &      Cd(47)*(12*Ce(411) + Co(1047) + Co(1049)))*
     &    UStop2(1,1)**2 - 
     &   2*(2*Cd(50)*Ce(411) + 
     &      Cd(47)*(4*Ce(411) + Co(1047) + Co(1049)))*
     &    UStop2(1,1)**2*UStop2(1,2) + 
     &   Cd(9948)*(4*Cd(47)*Ce(411) + 
     &      Cd(50)*(Co(1034) + Co(1037)))*UStop2(1,2)**2 - 
     &   (4*Cd(50)*Ce(411) + Cd(47)*(Co(1047) + Co(1049)))*
     &    UStop2(1,2)**3

        Co(2389) = Opt(46) + 
     &   Cd(47)**2*Cd(50)*Co(1049)*(2 + 2*UStop2(1,1))*UStop2(1,2)

        Co(1809) = Opt(46) + 
     &   2*Cd(47)**2*Cd(50)*Co(1049)*
     &    (-(Cd(9948)*UStop2(1,1)) + 
     &      (Cd(9948) + 3*UStop2(1,1))*UStop2(1,2))

        Co(1795) = Opt(47) + 
     &   Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*
     &    (3*Co(1047) + Co(1049)*(-2 - 2*UStop2(1,1))*UStop2(1,2))

        Co(2378) = Opt(47) - 
     &   Cd(47)**2*Cd(50)*Cd(9950)*
     &    (Cd(50)*(Co(1047) + Co(1049)) + 
     &      Cd(10103)*(-2*Co(1047) + Co(1049) - 
     &         2*Cd(9948)*Co(1049)*UStop2(1,1) + 
     &         2*Co(1049)*(Cd(9948) + 3*UStop2(1,1))*UStop2(1,2)))

	Co(1025) = Cd(17)*Co(1024) - Co(960)*MStop2(3)

        Co(1026) = Cd(20)*Co(1019) - 
     &   4*Co(1025)*MStop2(1)*MStop2(2)*MStop2(3)

        Co(103) = cI/(128.D0*Pi**3)*
     &   (Alfas2L*Cr(9)*(4*SA*Cd(13)*Cd(12836)*Cr(75) + 
     &        CA*TB*Co(102)*MStop2(1)*MStop2(3)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(3))

        Co(106) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*Cr(9)*(2*SA*Cd(13)*Cd(12836)*Cr(75) + 
     &          CA*TB*Co(105)*MStop2(2)*MStop2(3)))/
     &      (SB*SW2*TB*MStop2(2)*MStop2(3)))

        Co(1091) = -(AtC*Co(1090)*UCStop(1,3)**2) + 
     &   At*Cd(9948)*Co(1086)*UStop2(1,1)

	Co(1118) = At*Cd(9950)*Co(1111) - 2*AtC*Co(1117)*UCStop(1,3)**2

	Co(1138) = At*Co(1136) + 3*AtC*Co(1137)*UCStop(1,3)**2

        Co(114) = cI/(128.D0*Pi**3)*
     &   (Alfas2L*Cr(9)*(CA*TB*Co(113)*MStop2(1)*MStop2(3)**2 - 
     &        8*SA*Cd(13)*Cd(12836)*Cr(75)*MStop2(4)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(3)**2)

        Co(1154) = 3*AtC**2*Co(1137)*UCStop(1,3)**3 + 
     &   At**2*Co(1153)*UCStopC(1,3)

        Co(1176) = -(At*Co(1090)*UCStopC(1,3)**2) + 
     &   AtC*Co(1086)*UStop2(1,1)*UStop2(1,2)

	Co(1179) = Cd(47)**2*Cd(50)**2*Co(1071) + Cd(9950)*Co(1178)

        Co(1192) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Cd(10070)*Co(1067) + 
     &        Co(1191)*MStop2(2)))/(SB*SW2*Cd(47)*Cd(50)*MStop2(2))

        Co(1200) = AtC*Cd(9950)*Cd(10103)*Co(1191) - 
     &   2*At*Co(1199)*UCStopC(1,3)**2

        Co(1201) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Cd(9950)*Cd(10103)*Co(1067)*
     &         Co(1197) + Co(1200)*MStop2(2)))/
     &    (SB*SW2*TB*Cd(47)*Cd(50)*Cd(9950)*Cd(10103)*MStop2(2))

        Co(1203) = cI/(8.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Cd(9950)*Cd(10103)*Co(1067) - 
     &        Co(1199)*MStop2(2))*UCStopC(1,3)**2)/
     &    (SB*SW2*Cd(47)*Cd(50)*Cd(9950)*Cd(10103)*MStop2(2))

        Co(1212) = AtC**2*Co(1136)*UCStop(1,3) - 
     &   3*At**2*Co(1137)*UCStopC(1,3)**3

        Co(1213) = -(cI/(16.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(47)*Co(1067)*Co(1211) - 
     &          Co(1212)*MStop2(2)))/
     &      (SB*SW2*TB*Cd(47)*Cd(50)*MStop2(2)))

        Co(1215) = -(cI/(16.D0*Pi**3)*
     &     (Alfas2L*Cd(585)*
     &        (MGlmT2*Cd(47)*Cd(10092)*Co(1067) - 
     &          Co(1136)*MStop2(2)))/
     &      (SB*SW2*Cd(47)*Cd(50)*MStop2(2)))

	Co(1225) = AtC*Co(1153) - 3*At*Co(1137)*UCStopC(1,3)**2

        Co(1226) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Co(1067)*Co(1224) + 
     &        Co(1225)*MStop2(2))*UCStopC(1,3))/
     &    (SB*SW2*Cd(47)*Cd(50)*MStop2(2))

	Co(1236) = Co(1235)*MStop2(1) + Co(1233)*MStop2(3)

        Co(1272) = 4*Co(1271)*MStop2(3)*UStopC(1,1) + 
     &   2*Co(1269)*UUStopC(1,3)

        Co(1308) = 4*(At*Co(1307) + 
     &     MUE*TB*Cd(12114)*MStop2(1)**2*UCStop(1,3)*UUStop(1,1))

	Co(1310) = Co(1284) + Co(1274)*UStop2(1,1)

	Co(1311) = 2*Co(1278)*MStop2(3) + Co(1310)*UStop2(1,1)

	Co(1339) = 2*Co(1326) + AtC*Co(1338)

	Co(1340) = Co(1339)*UStopC(1,2) - AtC**2*Cd(11087)*UUStopC(1,2)

	Co(1343) = Co(1235)*MStop2(1) + 2*Co(1233)*MStop2(3)

	Co(1370) = Co(1269) + Co(1367)

        Co(1371) = 2*(Co(1369)*MStop2(3)*UStopC(1,1) + 
     &     Co(1370)*UUStopC(1,3))

	Co(1380) = Co(1338) - Co(1379)

	Co(137) = CA*TB*Co(136) + 8*SA*Cd(13)*Cd(833)*Cd(12836)*Cr(75)

	Co(1384) = Co(1354) - 2*Co(1343)*UStop2(1,1)

	Co(1393) = 2*Co(1389) + Co(1392)*UUStopC(1,2)

	Co(1406) = Co(1393) + Co(1388)*UStop2(1,1)*UStopC(1,2)

	Co(1449) = Co(1427) - Co(1448)

        Co(1450) = 2*Co(1449)*MStop2(2) - 
     &   Cd(12)*Cd(11607)*Cd(12836)*MStop2(3)

	Co(1472) = 2*Co(1468) + Co(1414)*UStop2(1,1)*UUStopC(2,2)

        Co(1493) = Co(1492)*MStop2(3) + 
     &   2*Co(1256)*MStop2(2)*UStop2(1,1)*UStopC(2,2)

	Co(1494) = MUEC*TB*Cd(11542) - 2*Co(1493)*UStop(2,1)

	Co(1500) = Co(1494) - 2*AtC*Cd(11537)*UCStop(1,3)*UStop(2,1)

	Co(1505) = Co(1458)*UUStopC(1,1) - Co(1466)*UUStopC(1,2)

	Co(1528) = Co(1527) - 2*Co(1508)*MStop2(3)

	Co(1529) = 4*Co(1458)*UUStopC(1,1) - 4*Co(1466)*UUStopC(1,2)

        Co(1680) = -(Co(1679)*MStop2(3)) + 
     &   Co(1676)*UStopC(2,2)*UUStopC(1,1)

	Co(1682) = Co(1613)*UStopC(1,1) - Co(1592)*UStopC(1,2)

	Co(1683) = 4*Co(1538)*UUStopC(1,2) - Co(1675)*UUStopC(1,3)

        Co(1698) = 4*Cd(47)**2*Cd(50)**2*Ce(411)*Co(1690) + 
     &   Cd(10103)*Co(1697)

	Co(1701) = At*Co(1698) + 2*AtC*Cd(10103)*Co(1700)*UCStop(1,3)**2

	Co(1729) = At*Co(1724) + 2*AtC*Co(1728)*UCStop(1,3)**2

        Co(1758) = AtC**2*Co(1728)*UCStop(1,3)**2 + 
     &   At**2*Co(1757)*UCStopC(1,3)**2

        Co(1765) = -(Cd(9950)*Co(1764)) + 
     &   2*Cd(47)**2*Cd(50)**2*Cd(9981)*Co(1071)*UStop2(1,1)

	Co(1766) = AtC*Co(1765) + 2*At*Co(1743)*UCStopC(1,3)**2

        Co(1789) = AtC*Co(1783) + 
     &   2*At*Cd(10103)*Co(1788)*UCStopC(1,3)**2

        Co(1804) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*Co(1803) - 
     &        Cd(50)*Cd(9950)*Co(1796)*MStop2(3))*UCStop(1,3))/
     &    (SB*SW2*Cd(47)**2*Cd(50)*Cd(9950)*MStop2(3))

        Co(1810) = AtC*Co(1795)*UCStop(1,3) + 
     &   At*Cd(9950)*Cd(10103)*Co(1809)*UCStopC(1,3)

        Co(1857) = At*Cd(9950)*Co(1850) + 
     &   2*AtC*Cd(10103)*Co(1856)*UCStop(1,3)**2

        Co(1868) = AtC**2*Co(1856)*UCStop(1,3)**2 + 
     &   At**2*Cd(9950)*Co(1867)*UCStopC(1,3)**2

	Co(1885) = AtC*Co(1881) + 2*At*Cd(9950)*Co(1884)*UCStopC(1,3)**2

        Co(1903) = -(cI/(256.D0*Pi**3)*
     &     (Alfas2L*Cr(9)**2*(CA*Co(3) + SA*Co(49)*MStop2(3)**2))/
     &      (SB*SW2*MStop2(3)**2))

	Co(1907) = -(SA*Co(22)) + CA*Co(1906)

	Co(1911) = SA*Co(29) + CA*Co(1910)

	Co(1914) = CA*Co(42) + SA*Co(47)

        Co(1915) = -(cI/(256.D0*Pi**3)*
     &     (Alfas2L*Cr(9)**2*(CA*Co(3) + SA*Co(49)*MStop2(3)**2))/
     &      (SB*SW2*MStop2(3)**2))

        Co(1917) = cI/(256.D0*Pi**3)*
     &   (Alfas2L*Cr(9)**2*
     &      (-(CA*Co(54)) + SA*Co(58)*MStop2(3)**2))/
     &    (SB*SW2*MStop2(3)**2)

	Co(1918) = -4*SA*Co(61) + CA*Co(65)

	Co(1920) = (Co(1919) + 2*TB*Co(71)*UStop2(1,1))/TB

	Co(1922) = -4*SA*Co(90) + CA*Co(93)

        Co(1923) = cI/(128.D0*Pi**3)*
     &   (Alfas2L*Cr(9)*(-4*CA*Cd(13)*Cd(12836)*Cr(75) + 
     &        SA*TB*Co(102)*MStop2(1)*MStop2(3)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(3))

        Co(1924) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*Cr(9)*(-2*CA*Cd(13)*Cd(12836)*Cr(75) + 
     &          SA*TB*Co(105)*MStop2(2)*MStop2(3)))/
     &      (SB*SW2*TB*MStop2(2)*MStop2(3)))

        Co(1925) = cI/(128.D0*Pi**3)*
     &   (Alfas2L*Cr(9)*(SA*TB*Co(113)*MStop2(1)*MStop2(3)**2 + 
     &        8*CA*Cd(13)*Cd(12836)*Cr(75)*MStop2(4)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(3)**2)

        Co(1928) = -(SA*TB*Co(136)) + 
     &   8*CA*Cd(13)*Cd(833)*Cd(12836)*Cr(75)

	Co(2021) = 2*Co(2019)*Cr(13) - 2*Co(2020)*MStop2(3)

	Co(2023) = 2*(Co(2022)*MStop2(1) + Co(2020)*MStop2(2))

	Co(2027) = Co(2023)*Cr(34) + 2*Co(2026)*MStop2(3)

	Co(2028) = Co(2021) - Co(2027)

        Co(2029) = A0delStop(1)*Co(2019) - 
     &   Co(2020)*Cr(13)*MStop2(1)*MStop2(3)

        Co(2030) = 1/2.D0*(A0delStop(2)*Co(2023)) + 
     &   Co(2026)*Cr(34)*MStop2(2)*MStop2(3)

        Co(2031) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(Co(2030)*MStop2(1) - Co(2029)*MStop2(2)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2)*MStop2(3)**2)

	Co(2052) = -(Cd(50)**2*Co(651)) + Cd(47)**2*Co(2051)

        Co(2054) = (-(Cd(47)**2*Co(677)) + 
     &     Cd(6458)*Co(650)*MStop2(3)*UStop2(1,1)*UUStopC(1,2))/
     &   (Cd(47)**2*MStop2(3))

        Co(2057) = MTy*(Co(2056)*UStop(1,1) - 
     &     6*MUEC*Cd(47)*Cd(12701)*UCStopC(1,3)*UStop(1,2))

        Co(2058) = (Co(2057) - MUEC*Cd(47)*MStop2(3)*UStop(1,2))/
     &   (TB*MStop2(3))

        Co(2065) = (-(Co(2064)*UCStop(1,3)) + 
     &     MUEC*Cd(50)**2*Co(706)*UStop2(1,1))/Cd(50)**2

	Co(2066) = Co(2065) + MUEC*Co(702)*UStopC(1,2)

        Co(2110) = 2*Cd(50)*Cd(9947)*Co(362)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Co(2109)*MGlpTmSt2(2)

        Co(2176) = At*Cd(9950)*Cd(10103)*Co(2167) - 
     &   2*AtC*Co(2175)*UCStop(1,3)**2

	Co(2197) = At*Co(2196) + 3*AtC*Co(1137)*UCStop(1,3)**2

        Co(2211) = 3*AtC**2*Co(1137)*UCStop(1,3)**3 + 
     &   At**2*Co(2210)*UCStopC(1,3)

        Co(2255) = AtC*Cd(9950)*Cd(10103)*Co(2247) - 
     &   2*At*Co(2254)*UCStopC(1,3)**2

        Co(2268) = AtC**2*Co(2196)*UCStop(1,3) - 
     &   3*At**2*Co(1137)*UCStopC(1,3)**3

	Co(2274) = AtC*Co(2210) - 3*At*Co(1137)*UCStopC(1,3)**2

	Co(2293) = Cd(9948)*Co(2286) + Co(2292)*UStop2(1,1)

	Co(2296) = At*Co(2293) + 2*AtC*Cd(10103)*Co(2295)*UCStop(1,3)**2

	Co(2317) = At*Cd(10103)*Co(2310) + 2*AtC*Co(2316)*UCStop(1,3)**2

        Co(2334) = AtC**2*Cd(9950)*Cd(10103)*Co(2295)*
     &    UCStop(1,3)**2 + At**2*Co(2333)*UCStopC(1,3)**2

        Co(2346) = Cd(10103)*Co(2345) + 
     &   2*Cd(47)**2*Cd(50)**2*Co(1710)*UStop2(1,1)*UStop2(1,2)**2

        Co(2373) = AtC*Co(2372) + 
     &   2*At*Cd(10103)*Co(2338)*UCStopC(1,3)**2

        Co(2390) = AtC*Co(2378)*UCStop(1,3) + 
     &   At*Cd(9950)*Cd(10103)*Co(2389)*UCStopC(1,3)

	Co(2393) = AtC*Co(2377)*UCStop(1,3) - At*Co(2392)*UCStopC(1,3)

        Co(2434) = At*Cd(9950)*Co(2429) + 
     &   2*AtC*Cd(10103)*Co(2433)*UCStop(1,3)**2

        Co(2450) = AtC**2*Cd(10103)*Co(2433)*UCStop(1,3)**2 + 
     &   At**2*Cd(9950)*Co(2449)*UCStopC(1,3)**2

        Co(2463) = AtC*Cd(10103)*Co(2462) + 
     &   2*At*Cd(9950)*Co(2449)*UCStopC(1,3)**2

	Co(247) = Co(245)*Cr(13) + 2*Co(246)*MStop2(3)

	Co(250) = Co(248)*MStop2(1) + Co(249)*MStop2(2)

	Co(254) = Co(250)*Cr(34) + 2*Co(253)*MStop2(3)

	Co(255) = Co(247) + Co(254)

        Co(256) = 1/2.D0*(A0delStop(1)*Co(245)) + 
     &   Co(246)*Cr(13)*MStop2(1)*MStop2(3)

        Co(257) = 1/2.D0*(A0delStop(2)*Co(250)) + 
     &   Co(253)*Cr(34)*MStop2(2)*MStop2(3)

        Co(258) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(Co(257)*MStop2(1) + Co(256)*MStop2(2)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2)*MStop2(3)**2)

	Co(30) = SA*Co(27) + CA*Co(29)

	Co(326) = 4*MUE*Co(321) - 4*MUEC*Co(325)

	Co(35) = (Co(34) - Cd(9782)*Co(4)*UStop2(1,1))/TB

	Co(36) = SA*Co(32) + CA*Co(35)

	Co(423) = MUE*Co(415) - MUEC*Co(422)

	Co(476) = 4*MUE*Co(468) - 4*MUEC*Co(475)

	Co(506) = MUE*Co(499) - MUEC*Co(505)

	Co(603) = 4*MUE*Co(595) - 4*MUEC*Co(602)

	Co(658) = Cd(50)**2*Co(651) - Cd(47)**2*Co(657)

	Co(693) = Co(688) - Co(692)

	Co(73) = (Co(72) - 2*TB*Co(71)*UStop2(1,1))/TB

	Co(785) = Co(763) + Co(784)

        Co(817) = (Co(816) + 
     &     4*MUE*Cd(12700)*MStop2(3)**2*UCStop(1,3)**2*
     &      UStop2(1,2)**2)/MStop2(3)**2

        Co(831) = (Co(830) - 
     &     4*MUE*Cd(12700)*MStop2(3)**2*UCStop(1,3)**2*
     &      UStop2(1,2)**2)/MStop2(3)**2

	Co(82) = ((Co(78) + Co(81))*UStop2(1,1))/TB

	Co(83) = -(SA*Co(73)) + 4*CA*Co(82)

        Co(887) = (Co(886) - 
     &     4*MUE*Cd(12700)*MStop2(3)**2*UCStop(1,3)**2*
     &      UStop2(1,1)**2)/(TB*MStop2(2)**2*MStop2(3)**2)

	Co(896) = (TB*Co(893) - Co(895))/(TB*MStop2(2)**2*MStop2(3))

	Co(915) = Co(913)*UStopC(1,2) + 2*Co(914)*UUStopC(3,2)

        Co(932) = Cd(9948)*Co(930)*MStop2(3)*UStop2(1,1) + 
     &   2*Co(931)*MStop2(2)*UUStopC(3,3)

        Co(941) = (Co(932) - 
     &     Co(836)*MStop2(3)*UStop2(1,1)*UStop2(1,2)*UUStopC(2,2))/
     &   MStop2(3)

	Co(1381) = 2*Co(1326) - 2*Co(1376) + AtC*Co(1380)

	Co(1396) = Co(1370) + 6*Co(1394)*MStop2(3)**2

        Co(1501) = AtC*Co(1500) - Co(1496)*MStop2(3) - 
     &   4*Co(1499)*MStop2(2)*UUStopC(3,3)

        Co(1503) = Co(1502)*UStopC(1,2)*UUStopC(2,1) - 
     &   2*Co(1490)*UUStopC(3,2) + Co(1501)*UUStopC(3,4)

        Co(1530) = 2*Co(1528)*MStop2(2) - 
     &   Co(1529)*MStop2(2)**2 + 
     &   4*Cd(12)*Cd(11375)*Cd(12836)*MStop2(3)**2

        Co(167) = (Co(166) + 
     &     6*CA*MTy*(Co(163) + TB*Co(164))*MStop2(3))/
     &   (MTy*MZ2*SW2*TB*MStop2(3))

        Co(171) = (Co(170) + 
     &     MTy*(Co(168) + CA*TB*Co(169))*MStop2(3))/
     &   (MTy*MZ2*SW2*TB*MStop2(3))

        Co(1730) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)**2*Cd(9950)*Ce(708)*Co(1067) + 
     &        Cd(47)*Co(1729)*MStop2(2) - 
     &        2*Cd(50)*Cd(9950)*Cd(9957)*Co(1717)*MStop2(2)*
     &         MStop2(3)))/
     &    (SB*SW2*Cd(47)**2*Cd(50)*Cd(9950)*MStop2(2)*MStop2(3))

        Co(1759) = -(cI/(32.D0*Pi**3)*
     &     (Alfas2L*(-(Cd(47)*Co(1758)*MStop2(2)) + 
     &          Cd(9950)*Co(4)*
     &           (MGlmT2*Cd(47)**2*Co(1067) + 
     &             Cd(50)*Co(1717)*MStop2(2)*MStop2(3))))/
     &      (SB*SW2*TB*Cd(47)**2*Cd(50)*Cd(9950)*MStop2(2)*
     &        MStop2(3)))

        Co(1790) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*Co(1789)*MStop2(2) + 
     &        Cd(9950)*Cd(10103)*
     &         (MGlmT2*Cd(47)**2*Cd(12542)*Co(1067) + 
     &           2*Cd(50)*Cd(9988)*Co(1717)*MStop2(2)*MStop2(3))))/
     &    (SB*SW2*Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*MStop2(2)*
     &      MStop2(3))

        Co(1957) = (-Co(1956) + 
     &     MTy*(Co(1953) + TB*Co(1954))*MStop2(3))/
     &   (MTy*MZ2*SW2*TB*MStop2(3))

        Co(1959) = (-Co(1958) + 
     &     6*MTy*SA*(Co(163) + TB*Co(164))*MStop2(3))/
     &   (MTy*MZ2*SW2*TB*MStop2(3))

        Co(2053) = (Co(2052) - Cd(50)**2*Co(659)*MGlpTmSt4(1) + 
     &     MUEC**2*TB*Cd(47)**2*MStop2(2)**2*UCStopC(1,3)**2)/
     &   (TB*Cd(47)**2*Cd(50)**2)

        Co(2059) = (Cd(47)**2*Co(693) + 
     &     MTy*MUEC*Cd(12701)*MGlpTmSt2(1)*UCStopC(1,3)**3 - 
     &     TB*Co(2058)*UCStopC(1,3)*UStopC(1,1))/(TB*Cd(47)**2)

        Co(2068) = (Cd(47)**2*Co(2066) + 
     &     MTy*TB*Cd(47)**2*Co(2067) + 
     &     TB*Co(2062)*UCStop(1,3)*UStopC(1,2))/(TB*Cd(47)**2)

	Co(2082) = AtC*Co(2079) - At*Co(2081) + TB*Co(2)*MGlpTmSt4(2)

        Co(2318) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*Co(2317)*MStop2(2) + 
     &        Cd(9950)*Cd(10103)*
     &         (MGlmT2*Cd(47)**2*Co(2304) - 
     &           2*Cd(50)*Cd(9957)*Co(1717)*MStop2(2)*MStop2(3))))/
     &    (SB*SW2*Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*MStop2(2)*
     &      MStop2(3))

        Co(2352) = Cd(9950)*Co(2346) + Cd(9948)*Co(2351) + 
     &   2*Cd(47)**2*Cd(50)**2*Cd(10103)*Co(1071)*UStop2(1,1)*
     &    UStop2(1,2)**2

        Co(2374) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*Co(2373)*MStop2(2) + 
     &        Cd(9950)*Cd(10103)*
     &         (MGlmT2*Cd(47)**2*Co(2365) + 
     &           2*Cd(50)*Cd(9988)*Co(1717)*MStop2(2)*MStop2(3))))/
     &    (SB*SW2*Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*MStop2(2)*
     &      MStop2(3))

        Co(311) = -((-(TB*Cd(10218)*Co(2)) + Co(308) + 
     &       Cd(152)*Cd(2597)*Co(310))/(TB*Cd(47)**2*MStop2(1)**2))

        Co(357) = (4*Co(356) + 
     &     (Cd(47)*(TB*Co(314) + Co(326)))/(TB*MStop2(3)))/
     &   (Cd(47)**2*MStop2(1))

        Co(401) = (-(TB*Cd(12323)*Co(2)) + Cd(12323)*Co(399) - 
     &     Cd(164)*Cd(2597)*Co(400))/(TB*Cd(47)**2*MStop2(1)**2)

        Co(435) = (-(MTy2*TB*Co(428)) + MGl2*Co(430) + 
     &     2*MTy*TB*Co(434))/TB

        Co(436) = (Co(423) + TB*(Co(406) + Co(435)))/
     &   (TB*Cd(47)**2*MStop2(1)**2*MStop2(3))

        Co(529) = (-(MTy2*TB*Co(518)) + MGl2*Co(521) + 
     &     2*MTy*TB*Co(528))/TB

        Co(530) = (Co(506) + TB*(Co(497) + Co(529)))/
     &   (TB*Cd(50)**2*MStop2(2)**2*MStop2(3))

        Co(535) = (TB*Cd(10259)*Co(2) + Co(533) - 
     &     Cd(179)*Cd(2681)*Co(534))/(TB*Cd(50)**2*MStop2(2)**2)

        Co(632) = (TB*Cd(50)*Co(587) + Cd(50)*Co(603) + 
     &     4*Co(631))/(TB*Cd(50)**2*MStop2(2)*MStop2(3))

        Co(660) = (Co(658) + Cd(50)**2*Co(659)*MGlpTmSt4(1) - 
     &     MUEC**2*TB*Cd(47)**2*MStop2(2)**2*UCStopC(1,3)**2)/
     &   (TB*Cd(47)**2*Cd(50)**2)

        Co(694) = (Cd(47)**2*Co(693) + 
     &     MTy*MUEC*Cd(12701)*MGlpTmSt2(1)*UCStopC(1,3)**3 + 
     &     TB*Co(685)*UCStopC(1,3)*UStopC(1,1))/(TB*Cd(47)**2)

	Co(732) = AtC*Co(726) + 2*TB*Co(728) - At*Co(731)

	Co(743) = AtC*Co(737) + At*Co(742) + TB*Co(2)*MGlpTmSt4(2)

        Co(787) = (Co(785) + 
     &     MUE*Cd(12700)*MStop2(3)**2*UCStop(1,3)**2*
     &      UStop2(1,2)**2 + Co(786)*MStop2(3)**2*UUStopC(1,2))/
     &   MStop2(3)**2

        Co(1112) = -(cI/(16.D0*Pi**3)*
     &      (Alfas2L*(Cd(47)*Cd(50)*Co(1111) + 
     &           2*Co(1109)*MStop2(3)*UStop2(1,1)))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))) + 
     &   Co(1105)/MStop2(1) + Co(1102)/MStop2(2)

        Co(1181) = cI/(32.D0*Pi**3)*
     &   (Alfas2L*((Cd(50)**2*Co(1180) - Co(1179)*MStop2(1))*
     &         MStop2(2)*MStop2(3) + 
     &        Cd(47)*Cd(9950)*MStop2(1)*
     &         (4*Cd(50)*Co(1077)*MStop2(2) + 
     &           Cd(47)*Co(1068)*MStop2(3)))*UCStopC(1,3)**2)/
     &    (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(1)*MStop2(2)*
     &      MStop2(3))

        Co(1193) = -(cI/(16.D0*Pi**3)*
     &      (Alfas2L*UStop2(1,1)*
     &         (-2*Co(1185)*MStop2(2) + 
     &           Cd(47)**2*Co(1120)*UStop2(1,2)))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(2))) + 
     &   Co(1188)/MStop2(1) + Co(1192)/MStop2(3)

        Co(1204) = cI/(16.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)**2*Co(1120) + 2*Co(1113)*MStop2(2))*
     &       UCStopC(1,3)**2)/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(2)) + 
     &   Co(1202)/MStop2(1) + Co(1203)/MStop2(3)

        Co(1309) = 2*Co(1303) + MTy*Co(1308) - 
     &   TB*(MUEC*Co(1298) + 
     &      2*At*MUE*Co(1299)*MStop2(1)*UUStop(1,1))

        Co(1473) = Co(1472)*UStopC(1,1) + 
     &   2*(Cd(11387)*Cd(12701) - Co(1419)*UUStopC(2,2))*
     &    UUStopC(3,1) + Co(1471)*UUStopC(3,3)

        Co(1685) = AtC*Co(1680) - 
     &   MStop2(3)*(Co(1684) + 2*Co(1682)*Cr(9)*MStop2(1) + 
     &      Co(1683)*MStop2(3))

        Co(1863) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(1857) + 
     &         2*Cd(9950)*Cd(10103)*Co(1847)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(1862)/MStop2(1) + Co(1858)/MStop2(2)

        Co(1889) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(1885) + 
     &         2*Cd(9950)*Co(1876)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(3)) + 
     &   Co(1888)/MStop2(1) + Co(1886)/MStop2(2)

        Co(1893) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(1857) + 
     &         2*Cd(9950)*Cd(10103)*Co(1847)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(1891)/MStop2(1) + Co(1892)/MStop2(2)

        Co(1900) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(1885) + 
     &         2*Cd(9950)*Co(1876)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(3)) + 
     &   Co(1898)/MStop2(1) + Co(1899)/MStop2(2)

        Co(2111) = (Cd(10014)*Co(2103))/MStop2(3) + 
     &   (Cd(9952)*Co(274)*MStop2(2) + 
     &      (Cd(50)*(-(Co(2108)*MStop2(1)) - 
     &           Cd(47)*Cd(9950)*Co(2110)*MStop2(3)))/
     &       (Cd(47)**2*Cd(9950)*MStop2(3)))/Cd(50)**2

        Co(2174) = -(cI/(16.D0*Pi**3)*
     &      (Alfas2L*(Cd(47)*Cd(50)*Co(2167) + 
     &           2*Co(1109)*MStop2(3)*UStop2(1,1)))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))) + 
     &   Co(2166)/MStop2(1) + Co(2173)/MStop2(2)

        Co(2189) = cI/(8.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(2175) - 
     &         Cd(9950)*Cd(10103)*Co(1113)*MStop2(3))*
     &       UCStop(1,3)**2)/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2188)/MStop2(1) + Co(2187)/MStop2(2)

        Co(2253) = cI/(16.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(2247) + 
     &         2*Co(1185)*MStop2(3)*UStop2(1,1)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(3)) + 
     &   Co(2246)/MStop2(1) + Co(2252)/MStop2(2)

        Co(2266) = cI/(8.D0*Pi**3)*
     &    (Alfas2L*(-(Cd(47)*Cd(50)*Co(2254)) + 
     &         Cd(9950)*Cd(10103)*Co(1113)*MStop2(3))*
     &       UCStopC(1,3)**2)/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2264)/MStop2(1) + Co(2265)/MStop2(2)

        Co(2445) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(2434) + 
     &         2*Cd(9950)*Cd(10103)*Co(1847)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2444)/MStop2(1) + Co(2440)/MStop2(2)

        Co(2472) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(2463) + 
     &         2*Cd(9950)*Cd(10103)*Co(1876)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2467)/MStop2(1) + Co(2471)/MStop2(2)

        Co(2481) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(2434) + 
     &         2*Cd(9950)*Cd(10103)*Co(1847)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2480)/MStop2(1) + Co(2477)/MStop2(2)

        Co(2493) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)*Co(2463) + 
     &         2*Cd(9950)*Cd(10103)*Co(1876)*MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2489)/MStop2(1) + Co(2492)/MStop2(2)

        Co(456) = MTy2*Co(447) - 4*MTy*Co(455) + 
     &   (MGl2*(-Co(451) + 2*MGl2*Cd(152)*Co(327)*MStop2(3)))/
     &    (TB*MStop2(3))

        Co(494) = -4*MTy*Co(483) - 
     &   (MTy2*Co(493) + MGl2*
     &       (Co(480) - 2*MGl2*Cd(179)*Co(446)*MStop2(3)))/
     &    (TB*MStop2(3))

        Co(924) = (2*Co(915)*MStop2(2) + 
     &      MStop2(3)*(Co(923) - 
     &         4*MUE*Cd(12700)*MStop2(3)*UCStop(1,3)**2*
     &          UStop2(1,1)**2))/(TB*MStop2(2)**2*MStop2(3)**2) - 
     &   Co(896)*UUStopC(2,2)

        Co(1144) = -(cI/(16.D0*Pi**3)*
     &      (Alfas2L*(4*Cd(47)**2*Cd(50)**2*Cd(10014)*Co(1128) - 
     &           MStop2(3)*
     &            (Cd(47)*Cd(50)*Co(1138) + 
     &              Cd(9957)*Co(1129)*MStop2(3)))*UCStop(1,3))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(3)**2)) + 
     &   Co(1141)/MStop2(1) + Co(1143)/MStop2(2)

        Co(1216) = 1/(16.D0*Pi**3)*
     &    (cI*(Alfas2L*Cd(585)*
     &           (Cd(47)**2*Co(1099) + Co(1129)*MStop2(2))*
     &           UStop2(1,1)*UStop2(1,2))/
     &         (Cd(47)**2*Cd(50)**2*MStop2(2)) + 
     &       (4*(-(cI*(Alfas2L*Cd(585)*Cd(10014)*Cd(10016))) + 
     &            4*Pi**3*(SB*SW2*Co(1215)*MStop2(3))))/
     &        MStop2(3)**2)/(SB*SW2) + Co(1214)/MStop2(1)

        Co(1227) = cI/(16.D0*Pi**3)*
     &    (Alfas2L*(4*Cd(47)**2*Cd(50)**2*Cd(10014)*Co(1221)*
     &          MStop2(2) - 
     &         Cd(9993)*(Cd(47)**2*Co(1099) + Co(1129)*MStop2(2))*
     &          MStop2(3)**2)*UCStopC(1,3))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(2)*MStop2(3)**2) + 
     &   Co(1223)/MStop2(1) + Co(1226)/MStop2(3)

        Co(1403) = Co(1339) - 2*Co(1376) - AtC*Co(1379) + 
     &   Co(1402)*MStop2(3)

        Co(2205) = -(cI/(16.D0*Pi**3)*
     &      (Alfas2L*(4*Cd(47)**2*Cd(50)**2*Cd(10014)*Co(1128) - 
     &           MStop2(3)*
     &            (Cd(47)*Cd(50)*Co(2197) + 
     &              Cd(9957)*Co(1129)*MStop2(3)))*UCStop(1,3))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(3)**2)) + 
     &   Co(2201)/MStop2(1) + Co(2204)/MStop2(2)

        Co(2279) = cI/(16.D0*Pi**3)*
     &    (Alfas2L*(-(Cd(9993)*Co(1129)*MStop2(3)**2) + 
     &         Cd(47)*(4*Cd(47)*Cd(50)**2*Cd(10014)*Co(1221) + 
     &            Cd(50)*Co(2274)*MStop2(3)))*UCStopC(1,3))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(3)**2) + 
     &   Co(2276)/MStop2(1) + Co(2278)/MStop2(2)

        Co(942) = (2*Co(940)*MStop2(2) + 
     &      MStop2(3)*(-(MUEC*Cd(11558)*Cd(12701)*UCStopC(1,3)) + 
     &         MStop2(3)*(Co(941) + 
     &            MUE*Cd(12700)*UCStop(1,3)**2*UStop2(1,1)**2)))/
     &    (TB*MStop2(2)**2*MStop2(3)**2) - Co(928)*UUStopC(2,2)

        Co(1098) = MUEC*Co(1078) + MUE*Co(1096) + 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (-(Cd(50)**2*Co(1085)*MStop2(2)*MStop2(3)) + 
     &         Cd(9950)*MStop2(1)*
     &          (2*Cd(47)*Cd(50)*Co(1097)*MStop2(2) + 
     &            (Cd(47)**2*Co(1080) - 2*Co(1091)*MStop2(2))*
     &             MStop2(3))))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(1)*
     &       MStop2(2)*MStop2(3))

        Co(1127) = MUEC*Co(1112) + MUE*Co(1125) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(1118) + 
     &         2*Cd(9950)*Co(1114)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(3)) + 
     &   Co(1126)/MStop2(1) + Co(1122)/MStop2(2)

        Co(1183) = MUE*Co(1169) + MUEC*Co(1181) - 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (Cd(50)**2*Co(1175)*MStop2(2)*MStop2(3) + 
     &         Cd(9950)*MStop2(1)*
     &          (2*Cd(47)*Cd(50)*Co(1182)*MStop2(2) + 
     &            (Cd(47)**2*Co(1170) - 2*Co(1176)*MStop2(2))*
     &             MStop2(3))))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(1)*
     &       MStop2(2)*MStop2(3))

        Co(1205) = MUE*Co(1193) + MUEC*Co(1204) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)**2*Co(281)*Co(1120) - 2*Co(1194)*MStop2(2)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*MStop2(2)) + 
     &   Co(1196)/MStop2(1) + Co(1201)/MStop2(3)

        Co(1731) = -(cI/(64.D0*Pi**3)*
     &      (Alfas2L*(Cd(47)*Cd(50)**2*Cd(10103)*Co(1716) + 
     &           (Co(1701) + 
     &              2*Cd(47)**2*Cd(9957)*Cd(10103)*Co(1689)*
     &               MStop2(2))*MStop2(3)))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(10103)*MStop2(3))) + 
     &   Co(1709)/MStop2(1) + Co(1730)*MStop2(1) + 
     &   Co(1715)/MStop2(2)

        Co(1791) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(Cd(47)*Cd(50)**2*Cd(9950)*Cd(10376)*Co(1739) - 
     &         (Co(1766) - 
     &            2*Cd(47)**2*Cd(9950)*Cd(9988)*Co(1689)*MStop2(2))
     &           *MStop2(3)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(3)) + 
     &   Co(1767)/MStop2(1) + Co(1790)*MStop2(1) + 
     &   Co(1774)/MStop2(2)

        Co(1805) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*UCStop(1,3)*
     &       ((-2*Co(1795)*MStop2(1)*MStop2(2) + 
     &            Cd(9950)*Cd(10103)*
     &             (Cd(50)**2*Co(1793)*MStop2(2) - 
     &               Cd(47)**2*MStop2(1)*
     &                (Co(1794) - 4*Co(1792)*MStop2(2)**2)))*
     &          MStop2(3) + 
     &         4*MGlmT2*Cd(47)*Cd(50)**2*Cd(9950)*Cd(9955)*
     &          Cd(10103)*Co(1104)*MStop2(1)*MStop2(2)*UStop2(1,1))
     &       )/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(1)*MStop2(2)*MStop2(3)) + Co(1804)*MStop2(1)

        Co(1820) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*UCStopC(1,3)*
     &       ((-((Cd(50)**2*Co(1793) + 2*Co(1809)*MStop2(1))*
     &               MStop2(2)) + 
     &            Cd(47)**2*MStop2(1)*
     &             (Co(1807) - 4*Co(1792)*MStop2(2)**2))*MStop2(3)+
     &           4*MGlmT2*Cd(47)*Cd(50)**2*Cd(9955)*Co(1104)*
     &          MStop2(1)*MStop2(2)*UStop2(1,2)))/
     &     (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(1)*MStop2(2)*
     &       MStop2(3)) + Co(1819)*MStop2(1)

        Co(1890) = MUE*Co(1863) + MUEC*Co(1889) + 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(1868) + 
     &         Cd(9950)*Co(4)*Co(1846)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(3)) + 
     &   Co(1873)/MStop2(1) + Co(1869)/MStop2(2)

        Co(1901) = MUE*Co(1893) + MUEC*Co(1900) + 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(1868) + 
     &         Cd(9950)*Co(4)*Co(1846)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(3)) + 
     &   Co(1895)/MStop2(1) + Co(1897)/MStop2(2)

        Co(2161) = MUEC*Co(2152) + MUE*Co(2160) - 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (Cd(50)**2*Cd(10103)*Co(1085)*MStop2(2)*MStop2(3) + 
     &         Cd(9950)*MStop2(1)*
     &          (Cd(47)**2*Co(2157)*MStop2(3) + 
     &            2*Cd(10103)*MStop2(2)*
     &             (-(Cd(47)*Cd(50)*Co(1097)) + Co(1091)*MStop2(3))
     &            )))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(1)*MStop2(2)*MStop2(3))

        Co(2190) = MUEC*Co(2174) + MUE*Co(2189) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(2176) + 
     &         2*Cd(9950)*Cd(10103)*Co(1114)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2186)/MStop2(1) + Co(2181)/MStop2(2)

        Co(2242) = MUE*Co(2226) + MUEC*Co(2241) - 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (Cd(50)**2*Cd(10103)*Co(2229)*MStop2(2)*MStop2(3) + 
     &         Cd(9950)*MStop2(1)*
     &          (Cd(47)**2*Co(2235)*MStop2(3) + 
     &            2*Cd(10103)*MStop2(2)*
     &             (Cd(47)*Cd(50)*Co(1182) - Co(1176)*MStop2(3)))))
     &      /(SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(1)*MStop2(2)*MStop2(3))

        Co(2267) = MUE*Co(2253) + MUEC*Co(2266) + 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(2255) + 
     &         2*Cd(9950)*Cd(10103)*Co(1194)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2259)/MStop2(1) + Co(2263)/MStop2(2)

        Co(2324) = -(cI/(64.D0*Pi**3)*
     &      (Alfas2L*(Cd(47)*Cd(50)**2*Cd(10103)*Co(2303) + 
     &           (Co(2296) + 
     &              2*Cd(47)**2*Cd(9957)*Cd(10103)*Co(1689)*
     &               MStop2(2))*MStop2(3)))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(10103)*MStop2(3))) + 
     &   Co(2301)/MStop2(1) + Co(2318)*MStop2(1) + 
     &   Co(2323)/MStop2(2)

        Co(2387) = -(cI/(64.D0*Pi**3)*
     &      (Alfas2L*(4*MGlmT2*Cd(47)*Cd(50)**2*Cd(9950)*Cd(9955)*
     &            Cd(10103)*Co(2212)*MStop2(1)*MStop2(2) + 
     &           (2*Co(2378)*MStop2(1)*MStop2(2) + 
     &              Cd(9950)*Cd(10103)*
     &               (-(Cd(50)**2*Co(1793)*MStop2(2)) + 
     &                 Cd(47)**2*MStop2(1)*
     &                  (Co(2377) - 4*Co(1792)*MStop2(2)**2)))*
     &            MStop2(3))*UCStop(1,3))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &         MStop2(1)*MStop2(2)*MStop2(3))) + Co(2386)*MStop2(1)

        Co(2404) = -(cI/(64.D0*Pi**3)*
     &      (Alfas2L*(4*MGlmT2*Cd(47)*Cd(50)**2*Cd(9955)*Co(2198)*
     &            MStop2(1)*MStop2(2) + 
     &           (Cd(50)**2*Co(1793)*MStop2(2) + 
     &              MStop2(1)*
     &               (2*Co(2389)*MStop2(2) - 
     &                 Cd(47)**2*
     &                  (Co(2392) - 4*Co(1792)*MStop2(2)**2)))*
     &            MStop2(3))*UCStopC(1,3))/
     &       (SB*SW2*Cd(47)**2*Cd(50)**2*MStop2(1)*MStop2(2)*
     &         MStop2(3))) + Co(2403)*MStop2(1)

        Co(2473) = MUE*Co(2445) + MUEC*Co(2472) + 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(2450) + 
     &         Cd(9950)*Cd(10103)*Co(4)*Co(1846)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2454)/MStop2(1) + Co(2458)/MStop2(2)

        Co(2494) = MUE*Co(2481) + MUEC*Co(2493) + 
     &   cI/(32.D0*Pi**3)*(Alfas2L*
     &       (Cd(47)*Cd(50)*Co(2450) + 
     &         Cd(9950)*Cd(10103)*Co(4)*Co(1846)*MStop2(3)))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2483)/MStop2(1) + Co(2486)/MStop2(2)

        Co(1164) = MUE*Co(1144) + MUEC*Co(1157) + 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (-(Co(1145)*MStop2(3)**2) + 
     &         Cd(47)*(4*Cd(47)*Cd(50)**2*Cd(10014)*Co(1146) + 
     &            Cd(50)*Co(1154)*MStop2(3))))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*MStop2(3)**2) + 
     &   Co(1161)/MStop2(1) + Co(1163)/MStop2(2)

        Co(1228) = MUE*Co(1216) + MUEC*Co(1227) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (4*Cd(47)**2*Cd(50)**2*Cd(10014)*Co(1206)*MStop2(2) - 
     &         (Cd(47)**2*Co(1208) + Co(1210)*MStop2(2))*
     &          MStop2(3)**2))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*MStop2(2)*MStop2(3)**2)+
     &     Co(1220)/MStop2(1) + Co(1213)/MStop2(3)

        Co(2141) = (-(MUE*Co(2111)) + MUEC*Co(2131) + 
     &      2*MTy*TB*Co(2140))/TB - 
     &   (MGlpT2*Cd(9950)*Co(2119)*MStop2(3) - 
     &      MUE**2*Co(2102)*UCStop(1,3)**2 + 
     &      MUEC**2*Co(2112)*UCStopC(1,3)**2)/
     &    (Cd(47)**2*Cd(50)**2*Cd(9950)*MStop2(3))

        Co(2222) = MUE*Co(2205) + MUEC*Co(2217) + 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (-(Co(1145)*MStop2(3)**2) + 
     &         Cd(47)*(4*Cd(47)*Cd(50)**2*Cd(10014)*Co(1146) + 
     &            Cd(50)*Co(2211)*MStop2(3))))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*MStop2(3)**2) + 
     &   Co(2219)/MStop2(1) + Co(2221)/MStop2(2)

        Co(2282) = MUE*Co(2271) + MUEC*Co(2279) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*
     &       (4*Cd(47)**2*Cd(50)**2*Cd(10014)*Co(1206) - 
     &         MStop2(3)*(Cd(47)*Cd(50)*Co(2268) + 
     &            Co(1210)*MStop2(3))))/
     &     (SB*SW2*TB*Cd(47)**2*Cd(50)**2*MStop2(3)**2) + 
     &   Co(2273)/MStop2(1) + Co(2281)/MStop2(2)

        Co(1840) = 1/(64.D0*Pi**3)*
     &   (4*cI*(Alfas2L*MGlmT2*Cd(47)*Cd(50)**2*Cd(9950)*Cd(9955)*
     &         Cd(10103)*Co(1104)*Co(1806)*MStop2(1)*MStop2(2)) + 
     &      (64*Pi**3*(SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*
     &            Cd(10103)*MStop2(1)*
     &            (MUE*Co(1805) + MUEC*Co(1820) + MTy2*Co(1831) + 
     &              MGl2*Co(1839) + Co(1818)*MStop2(1))*MStop2(2))-
     &           cI*(Alfas2L*
     &            (2*Co(1810)*MStop2(1)*MStop2(2) + 
     &              Cd(9950)*Cd(10103)*
     &               (-(Cd(50)**2*Co(163)*Co(1793)*MStop2(2)) + 
     &                 Cd(47)**2*MStop2(1)*
     &                  (Co(1808) - 
     &                    4*Co(163)*Co(1792)*MStop2(2)**2)))))*
     &       MStop2(3))/
     &    (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &      MStop2(1)*MStop2(2)*MStop2(3))

        Co(2426) = 1/(64.D0*Pi**3)*
     &   (64*Pi**3*(SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*
     &         Cd(10103)*MStop2(1)*
     &         (MUE*Co(2387) + MUEC*Co(2404) + MTy2*Co(2416) + 
     &           MGl2*Co(2425) + Co(2402)*MStop2(1))*MStop2(2)*
     &         MStop2(3)) - 
     &      cI*(Alfas2L*(4*MGlmT2*Cd(47)*Cd(50)**2*Cd(9950)*
     &            Cd(9955)*Cd(10103)*Co(2388)*MStop2(1)*MStop2(2)+
     &             (2*Co(2390)*MStop2(1)*MStop2(2) + 
     &              Cd(9950)*Cd(10103)*
     &               (-(Cd(50)**2*Co(163)*Co(1793)*MStop2(2)) + 
     &                 Cd(47)**2*MStop2(1)*
     &                  (Co(2393) - 
     &                    4*Co(163)*Co(1792)*MStop2(2)**2)))*
     &            MStop2(3))))/
     &    (SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &      MStop2(1)*MStop2(2)*MStop2(3))

	Co(1382) = Co(1381)*UStopC(1,2) - AtC**2*Cd(11253)*UUStopC(1,2)

	Co(1383) = Co(1371) - Co(1382)

        Co(1397) = Co(1395)*MStop2(3)*UStopC(1,1) + 
     &   2*Co(1396)*UUStopC(1,3)

	Co(1404) = Co(1403)*UStopC(1,2) - AtC**2*Cd(11103)*UUStopC(1,2)

	Co(1405) = Co(1397) - Co(1404)

	Co(1504) = Co(1503) - Co(1473)*MStop2(3)

	Co(1686) = Cd(17)*Co(1685) - 2*Co(1565)*MStop2(3)**2

	Co(1687) = Cd(20)*Co(1671) + 2*Co(1686)*MStop2(1)*MStop2(2)

	Co(172) = Co(167) - Co(171)

	Co(187) = 1/2.D0*(A0delStop(2)*Co(167)) - Co(186)*Cr(34)*MStop2(2)

	Co(189) = 1/2.D0*(A0delStop(1)*Co(171)) - Co(188)*Cr(13)*MStop2(1)

        Co(190) = (3*cI)/(16.D0*Pi)*
     &   (Alfa2L*MZ2*(Co(187)*MStop2(1) - Co(189)*MStop2(2)))/
     &    (SB*MStop2(1)*MStop2(2))

	Co(1913) = SA*Co(35) + CA*Co(1912)

	Co(1921) = -4*SA*Co(82) + CA*Co(1920)

	Co(191) = 2*Co(186) - Co(167)*Cr(34)

	Co(192) = 2*Co(188) - Co(171)*Cr(13)

	Co(193) = Co(191) - Co(192)

	Co(1960) = Co(1957) - Co(1959)

        Co(1973) = 1/2.D0*(A0delStop(1)*Co(1957)) - 
     &   Co(1972)*Cr(13)*MStop2(1)

        Co(1977) = 1/2.D0*(A0delStop(2)*Co(1959)) - 
     &   Co(1976)*Cr(34)*MStop2(2)

        Co(1978) = (3*cI)/(16.D0*Pi)*
     &   (Alfa2L*MZ2*(Co(1977)*MStop2(1) - Co(1973)*MStop2(2)))/
     &    (SB*MStop2(1)*MStop2(2))

	Co(1979) = 2*Co(1972) - Co(1957)*Cr(13)

	Co(1980) = 2*Co(1976) - Co(1959)*Cr(34)

	Co(1981) = Co(1979) - Co(1980)

	Co(2076) = M_3*Co(2059) - M_3C*Co(2068)

	Co(2353) = AtC*Co(2352) + 2*At*Co(2333)*UCStopC(1,3)**2

	Co(397) = -2*MGl2*Co(311) + Ce(389)*Co(396)

	Co(398) = Ce(388)*Co(357) - 2*Co(397)*MStop2(1)

	Co(458) = 2*MTy2*Co(401) + Ce(392)*Co(436)

	Co(531) = 2*MTy2*Co(464) + Ce(401)*Co(530)

	Co(633) = 2*MGl2*Co(535) + Ce(399)*Co(584)

	Co(634) = Ce(398)*Co(632) - 2*Co(633)*MStop2(2)

	Co(722) = M_3*Co(694) + M_3C*Co(714)

	Co(832) = Co(787) + Co(831)*Cr(13)

	Co(833) = Cd(22)*Co(817) + 4*Co(832)*MStop2(1)

	Co(943) = Co(942) + Co(924)*Cr(34)

	Co(944) = Cd(16)*Co(887) - 4*Co(943)*MStop2(2)

        Co(1312) = 2*(Co(1309) - Co(1311)*MStop2(1) + 
     &     2*Co(1289)*UUStopC(1,2))

	Co(1341) = Co(1272) - Co(1340) - Co(1312)*UUStopC(1,1)

        Co(1342) = 4*Cd(9794)*Cd(12701)*UCStopC(1,3)*
     &    UStop2(1,1)**2*UStopC(1,1) + 
     &   (Co(1341) - 4*Co(1236)*MStop2(3)*UStop2(1,1)*
     &       UUStopC(1,1))/MStop2(3)**2

        Co(1385) = -4*Cd(9794)*Cd(12701)*UCStopC(1,3)*
     &    UStop2(1,1)**2*UStopC(1,1) + 
     &   (Co(1383) - Co(1384)*MStop2(3)*UUStopC(1,1))/MStop2(3)**2

        Co(1407) = (Co(1405) + 
     &     Cd(9794)*Cd(12701)*MStop2(3)**2*UCStopC(1,3)*
     &      UStop2(1,1)**2*UStopC(1,1) - 
     &     Co(1406)*MStop2(3)**2*UUStopC(1,1))/MStop2(3)**2

        Co(1506) = 2*MStop2(2)*
     &    (Co(1504) - 2*Co(1505)*MStop2(2)) + 
     &   2*Cd(12)*Cd(11527)*Cd(12836)*MStop2(3)**2

        Co(2094) = -(Cd(50)*Co(2093)) + 
     &   (2*TB*Cd(50)**2*Co(2076) + Co(2082))*MStop2(3)

        Co(2095) = Cd(47)*
     &    (Cd(50)**2*Co(2074) + Cd(47)*Co(2094)) + 
     &   Cd(50)**2*Co(2075)*MStop2(3)

        Co(457) = (-(Cd(47)*(TB*Co(438) + Co(443))) + 
     &     4*TB*Co(456)*MStop2(3))/
     &   (TB*Cd(47)**2*MStop2(1)*MStop2(3))

        Co(495) = (-(Cd(50)*(TB*Co(465) + Co(476))) + 
     &     4*TB*Co(494)*MStop2(3))/
     &   (TB*Cd(50)**2*MStop2(2)*MStop2(3))

        Co(744) = -(Cd(50)*Co(732)) + 
     &   (2*TB*Cd(50)**2*Co(722) + Co(743))*MStop2(3)

        Co(745) = Cd(47)*(Cd(50)**2*Co(720) - Cd(47)*Co(744)) - 
     &   Cd(50)**2*Co(721)*MStop2(3)

        Co(746) = (Co(678) + 
     &     (-Co(745) + Cd(50)**2*Co(664)*UUStopC(1,1))/
     &      (Cd(47)**2*Cd(50)**2*MStop2(3)))/TB

        Co(1165) = MGl*(Co(1098) + MTy2*Co(1127) + 
     &      MTy*Co(1164)) + 
     &   cI/(16.D0*Pi**3)*(Alfas2L*MGl*Co(1066)*UCStop(1,3))/
     &     (MTy*SB*SW2*TB*Cd(47)*Cd(50))

        Co(1229) = MGl*(Co(1183) + MTy2*Co(1205) + 
     &      MTy*Co(1228)) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*MGl*Co(1066)*UCStopC(1,3))/
     &     (MTy*SB*SW2*TB*Cd(47)*Cd(50))

        Co(2223) = MGl*(Co(2161) + MTy2*Co(2190) + 
     &      MTy*Co(2222)) + 
     &   cI/(16.D0*Pi**3)*(Alfas2L*MGl*Co(1066)*UCStop(1,3))/
     &     (MTy*SB*SW2*TB*Cd(47)*Cd(50))

        Co(2283) = MGl*(Co(2242) + MTy2*Co(2267) + 
     &      MTy*Co(2282)) - 
     &   cI/(16.D0*Pi**3)*(Alfas2L*MGl*Co(1066)*UCStopC(1,3))/
     &     (MTy*SB*SW2*TB*Cd(47)*Cd(50))

        Co(2375) = cI/(64.D0*Pi**3)*
     &    (Alfas2L*(-(Co(2353)*MStop2(3)) + 
     &         Cd(47)*Cd(9950)*Cd(10103)*
     &          (Cd(50)**2*Co(2359) + 
     &            2*Cd(47)*Cd(9988)*Co(1689)*MStop2(2)*MStop2(3))))
     &      /(SB*SW2*Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*
     &       MStop2(3)) + Co(2357)/MStop2(1) + 
     &   Co(2374)*MStop2(1) + Co(2364)/MStop2(2)

	Co(1029) = Co(1026) - Co(833)*MStop2(2)*MStop2(3)**2

	Co(1408) = Co(1407) + Co(1385)*Cr(13)

	Co(1409) = Cd(22)*Co(1342) + 4*Co(1408)*MStop2(1)

	Co(1531) = Cd(13679)*Co(1506) + Co(1450)*MStop2(3)

	Co(1532) = Cd(13687)*Co(1530) - 4*Co(1531)*MStop2(2)

	Co(1688) = M_3C*Co(1165) + M_3*Co(1229)

	Co(1749) = (Cd(50)*Co(1532) + 8*Co(1748)*MStop2(3))/Cd(50)

	Co(1842) = 8*Co(4)*Co(1841) + Co(1409)*MStop2(3)

	Co(1843) = Co(1687) + Co(1842)*MStop2(2)*MStop2(3)

        Co(1844) = Cd(47)*Co(1843) - 
     &   8*MGlmT2*Co(4)*Co(1104)*MStop2(2)*MStop2(3)**2

        Co(2096) = Cd(47)**2*Cd(50)**2*Co(2069)*Co(2070) - 
     &   Cd(10103)*Co(2095)

        Co(2097) = -(Cd(9950)*Co(2096)) + 
     &   Cd(47)**2*Cd(50)**2*Cd(10103)*Cd(11622)*Co(679)*
     &    UStop2(1,2)**2

	Co(2284) = M_3C*Co(2223) + M_3*Co(2283)

	Co(2341) = Co(1532) + 8*Co(4)*Co(1745)*MStop2(3)

        Co(2342) = Cd(50)*Co(2341) - 
     &   8*MGlmT2*Cd(425)*Cd(9955)*Co(1047)*MStop2(3)**2

	Co(459) = Ce(394)*Co(457) + 2*Co(458)*MStop2(1)

	Co(532) = Ce(403)*Co(495) + 2*Co(531)*MStop2(2)

	Co(747) = 4*Ce(409)*Co(660) + 2*Ce(411)*Co(746)

        Co(2098) = (Co(2054) - 
     &     (Co(2097) + Cd(50)**2*Cd(9950)*Cd(10103)*Co(664)*
     &         UUStopC(1,1))/
     &      (Cd(47)**2*Cd(50)**2*Cd(9950)*Cd(10103)*MStop2(3)))/TB

	Co(2100) = Co(459) - Co(532) + Co(634)

	Co(1027) = -4*Co(747) + Co(2100)

        Co(1028) = 8*Co(307) - 2*Co(398) + 2*Co(1027) - 
     &   Co(944)*MStop2(2)

        Co(1902) = 1/(256.D0*Pi**3)*
     &   (256*Pi**3*(SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*
     &         MStop2(1)*(MTy*Co(1688) + 
     &           MGl*(MTy2*Co(1840) + 
     &              MTy*(MUE*Co(1731) + MUEC*Co(1791) + 
     &                 MTy2*Co(1890) + MGl2*Co(1901) + 
     &                 Co(1759)*MStop2(1))))*MStop2(2)*MStop2(3)**2
     &         ) + cI*(Alfas2L*MGl*
     &         (MTy*Cd(47)*Cd(50)**2*Cd(9950)*
     &            (-Co(1844) + Cd(47)*Co(1749)*MStop2(1)) - 
     &           8*MTy*Cd(47)*Cd(50)**2*Cd(9950)*Co(1740)*
     &            MStop2(1)*MStop2(2)*MStop2(3) + 
     &           4*(-2*MTy*Co(1744)*MStop2(1)*MStop2(2) + 
     &              Cd(9950)*
     &               (Co(1736)*Co(1737) - 
     &                 2*MTy*Cd(47)**2*Co(4)*Co(1689)*MStop2(1)*
     &                  MStop2(2)**2))*MStop2(3)**2)))/
     &    (MGl*MTy*SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*
     &      MStop2(1)*MStop2(2)*MStop2(3)**2)

        Co(2495) = 1/(256.D0*Pi**3)*
     &   (256*Pi**3*(SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*
     &         Cd(10103)*MStop2(1)*
     &         (MTy*Co(2284) + 
     &           MGl*(MTy2*Co(2426) + 
     &              MTy*(MUE*Co(2324) + MUEC*Co(2375) + 
     &                 MTy2*Co(2473) + MGl2*Co(2494) + 
     &                 Co(2340)*MStop2(1))))*MStop2(2)*MStop2(3)**2
     &         ) + cI*(Alfas2L*MGl*
     &         (MTy*Cd(47)*Cd(50)*Cd(9950)*Cd(10103)*
     &            (-(Cd(50)*Co(1844)) + Cd(47)*Co(2342)*MStop2(1))+
     &             8*MTy*Cd(47)*Cd(50)**2*Cd(9950)*Cd(10103)*
     &            Co(2327)*MStop2(1)*MStop2(2)*MStop2(3) + 
     &           4*(-2*MTy*Co(2334)*MStop2(1)*MStop2(2) + 
     &              Cd(9950)*Cd(10103)*
     &               (Co(1736)*Co(1737) - 
     &                 2*MTy*Cd(47)**2*Co(4)*Co(1689)*MStop2(1)*
     &                  MStop2(2)**2))*MStop2(3)**2)))/
     &    (MGl*MTy*SB*SW2*TB*Cd(47)**2*Cd(50)**2*Cd(9950)*
     &      Cd(10103)*MStop2(1)*MStop2(2)*MStop2(3)**2)

        Co(1030) = cI/(256.D0*Pi**3)*
     &   (Alfas2L*(-Co(1029) + 
     &        TB*Co(1028)*MStop2(1)*MStop2(2)*MStop2(3)**2))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2)*MStop2(3)**2)

	Co(2099) = 4*Ce(409)*Co(2053) + 2*Ce(411)*Co(2098)

        Co(2101) = -2*Co(398) + 8*Co(2099) + 2*Co(2100) - 
     &   Co(944)*MStop2(2)

        Co(2142) = cI/(256.D0*Pi**3)*
     &   (Alfas2L*(Co(1029) - 
     &        TB*(Co(2101) - 8*Co(2141))*MStop2(1)*MStop2(2)*
     &         MStop2(3)**2))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2)*MStop2(3)**2)

	end
