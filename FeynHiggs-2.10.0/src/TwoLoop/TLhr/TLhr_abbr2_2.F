* TLhr_abbr2_2.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr2_2
	implicit none

#include "TLhr.h"


        Ce(656) = Cd(50)**2 - Cd(251)*Cd(2681)*MGlpTmSt2(2) - 
     &   Cd(4187)*MStop2(2)

        Ce(659) = Cd(50)*Cd(9966) + 
     &   Cd(251)*(Cd(10802) + Cd(9948)*Cd(9983)*MStop2(1))

        Ce(671) = Cd(50)*Cd(251)*(Cd(9948)*Ce(669) - Ce(670)) + 
     &   2*Cd(179)*Cd(795)*MStop2(2)*MStop2(3)

        Ce(673) = 2*Cd(179)*Cd(10219)*MStop2(2)*MStop2(3) + 
     &   Cd(50)*Cd(251)*(-Ce(672) + MUE*Cd(9948)*UStop2(1,1)**2)

        Ce(674) = AtC*Ce(671) + At*Ce(673) + 
     &   4*Cd(179)**2*Cd(10052)*MStop2(3)

        Ce(676) = Cd(9948)*Ce(675) + Cd(10658)*MStop2(2) + 
     &   MGlmT2*Cr(75)**2*UStop2(1,1)

        Ce(682) = AtC*UCStop(1,3)**2 + 
     &   2*M_3*(3*UStop2(1,1)**2 + UStop2(1,2)**2)

	Ce(686) = MStop2(2) + MGlmT2*(1 + UStop2(1,2))

        Ce(696) = 2*At*Cd(251)*Cd(10092)*Cr(5) + 
     &   (-2*AtC*Cd(251)*Cd(10087) + Ce(695))*UCStop(1,3)

        Ce(70) = 2*Ce(69)*UStop2(1,2) + 
     &   AtC*(Ce(62) - 2*AtC*UCStop(1,3)**2*UStop2(1,2))

        Ce(700) = AtC*Cd(10032) + MStop2(2)*UCStop(1,3) + 
     &   At*Cr(5)*UStop2(1,1)*UStop2(1,2)

        Ce(701) = MStop2(2) + 
     &   At*(M_3*UCStopC(1,3)**2 + 
     &      M_3C*UStop2(1,1)*UStop2(1,2))

        Ce(715) = -(Cd(9948)*Ce(714)) + 
     &   2*Cd(50)*Cd(251)*Ce(670)*MStop2(2) + 
     &   Cd(10382)*Cd(12483)*MStop2(3)

        Ce(717) = 2*Cd(50)*Cd(251)*Ce(672)*MStop2(2) + 
     &   Cd(10388)*Cd(12483)*MStop2(3) - 
     &   MUE*Cd(9948)*Ce(716)*UStop2(1,1)**2

	Ce(724) = -Ce(554) + At*UStop2(1,1)*UStop2(1,2)

        Ce(725) = 2*M_3*Cd(9948)*Cd(9983) + 
     &   Cd(10398)*UCStop(1,3)**2 + 2*M_3*Cr(75)**2*UStop2(1,1)

	Ce(731) = -Ce(560) + AtC*UStop2(1,1)*UStop2(1,2)

        Ce(734) = Cd(9948)*Cd(10861) + 
     &   Cd(10140)*Cr(85)*UCStopC(1,3)**2 + 
     &   2*M_3C*Cr(75)**2*UStop2(1,1)

        Ce(741) = 2*Cd(1952)*Cd(9948)*Cr(75) + 
     &   Cd(10722)*MStop2(3) - 2*MTy2*Cr(75)*UStop2(1,2)

        Ce(742) = 2*At*Cd(10092)*Cr(5) + 
     &   (-2*AtC*Cd(10087) + Ce(741))*UCStop(1,3)

	Ce(744) = -(At*Cd(10095)) + Ce(694)

        Ce(75) = Cd(435)*Ce(73) + 
     &   4*MTy*S2A*MStop2(3)*UCStop(1,3) - 2*Ce(74)*UCStop(1,3)**2

        Ce(754) = MTy2*(4*M_3*MUE*Ce(753) + 
     &      4*M_3C*MUEC*Ce(753)) + 
     &   Cd(48)*Cd(179)*Cd(424)*Cd(2681)*UStop2(1,1)

        Ce(760) = -2*MTy*
     &    (M_3C*UCStop(1,3) + M_3*UCStopC(1,3)) + 
     &   MStop2(2)*UStop2(1,2)

        Ce(761) = UCStopC(1,3)*
     &    (-MStop2(2) + 2*M_3*MTy*UCStopC(1,3)) + 
     &   2*M_3C*MTy*UStop2(1,1)*UStop2(1,2)

        Ce(762) = MUEC*Ce(761)*UCStopC(1,3) + 
     &   MUE*UStop2(1,1)*(-Ce(760) + MGl2*UStop2(1,1))*UStop2(1,2)

        Ce(763) = MGlpT2*Cd(10388)*Cd(10603) + 
     &   2*Cd(10604)*Ce(762) + MUE*Cd(9948)*Ce(759)*UStop2(1,1)**2

	Ce(764) = MGl2*MGlmT2 + MStop2(1)*MStop2(2) - MGlmT2*MStop2(3)

        Ce(766) = MGlpTmSt2(1) + 6*M_3*MTy*UCStopC(1,3) - 
     &   MStop2(1)*UStop2(1,2)

        Ce(779) = UCStop(1,3)*
     &    (-MStop2(2) + 2*M_3C*MTy*UCStop(1,3)) + 
     &   2*M_3*MTy*UStop2(1,1)*UStop2(1,2)

        Ce(780) = MUE*Ce(779)*UCStop(1,3) + 
     &   MUEC*UStop2(1,1)*(-Ce(760) + MGl2*UStop2(1,1))*UStop2(1,2)

        Ce(781) = MGlpT2*Cd(10382)*Cd(10603) + 
     &   2*Cd(10604)*Ce(780) + MUEC*Cd(9948)*Ce(759)*UStop2(1,1)**2

        Ce(791) = MGl2*MTy*Cd(13) + 
     &   2*MTy2*(M_3C*Cd(10382) + M_3*Cd(10388))

	Ce(792) = MGl2*MTy*Cd(13) + MTy*Cd(13)*Cd(48) + 2*MTy2*Ce(634)

        Ce(793) = MTy*Cd(13)*Cd(48)*Cd(10603) + 
     &   2*Cd(10604)*Ce(791) + 2*Cd(10602)*Ce(792)*MStop2(2)

        Ce(80) = UCStop(1,3)*(-2*MTy*S2A + Ce(79)*UCStop(1,3)) + 
     &   Ce(73)*UStop2(1,1)*UStop2(1,2)

        Ce(809) = M_3C*Ce(808) - 2*MTy*Cd(13)*Cd(118)*Cr(75) + 
     &   M_3*MUE*Cd(10702)*UStop2(1,1)

        Ce(810) = MTy2*(M_3C*MUEC*Cd(9983) + 
     &      M_3*Cd(10063)) + MTy*Cd(13)*Cr(75)*MStop2(1)

        Ce(82) = UCStop(1,3)*
     &    (S2A*Cd(757) + C2A*MUE*UCStop(1,3)) + 
     &   Ce(73)*UStop2(1,1)*UStop2(1,2)

        Ce(819) = MUE*Cd(179)*UCStop(1,3)**2 + 
     &   MTy2*(2*MUE*UCStop(1,3)**2 + 
     &      2*MUEC*UStop2(1,1)*UStop2(1,2)**2)

        Ce(821) = MTy2*(2*M_3C*Ce(819) + 2*M_3*Ce(820)) + 
     &   MTy*Cd(13)*Cd(13061)*MGlpTmSt2(2)

	Ce(833) = 2*Cd(11045) - 2*Cd(11048) + AtC*Ce(832)

        Ce(835) = Cd(11242) - Ce(833)*UStopC(1,2) - 
     &   Ce(834)*MStop2(3)*UUStopC(1,1)

        Ce(844) = -(MUE*Cd(12635)*UStop2(1,1)*UStop2(1,2)) + 
     &   MTy*MUEC*Cd(9794)*MGlpTmSt2(1)*UCStopC(1,3)**2*
     &    UStopC(1,1) + Cd(12634)*UUStopC(1,1)

        Ce(854) = (Cd(9948)*Cd(12170)*UStop(2,2) + 
     &     (Cd(50)*(Cd(47)*Cd(12191) - 
     &          2*MTy*Cd(50)*Cd(10987)*UStop2(1,1))*UStopC(1,1))/
     &      (Cd(47)*MStop2(3)))/Cd(50)**2

        Ce(865) = (MTy*Cd(50)**2*Ce(864)*MStop2(3) + 
     &     Cd(47)**2*(Cd(50)**2*Ce(862) + 
     &        Cd(13)*Cd(12170)*MStop2(3)*UCStopC(1,3)))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

	Ce(869) = Ce(868)*UStop2(1,1) - UStop2(1,2)

	Ce(873) = MTy2 + 4*Cd(800)*Cd(10705) + 5*MStop2(1)

        Ce(874) = Cd(9948)*
     &    (Ce(872) + 12*Cd(47)*Cd(10705)*UStop2(1,1)) - 
     &   Cd(47)*Ce(873)*UStop2(1,2)

        Ce(88) = -(AtC*S2A*Cd(423)) + 2*C2A*MUE*Cd(577) - 
     &   2*Ce(73)*MStop2(4)*UCStopC(1,3)**2

        Ce(879) = 2*MUE*MStop2(2)*UCStop(1,3)**2 + 
     &   MUEC*(Ce(877)*UStop2(1,1) - Ce(878)*UStop2(1,2))

        Ce(882) = Cd(50)*Cd(6458)*Cd(12289)*MStop2(3) + 
     &   Cd(47)**2*(At*Ce(881) + AtC*Ce(879)*UStop2(1,1))

        Ce(897) = (Cd(50)*MGlpTmSt2(1) + Cd(47)*MGlpTmSt2(2) + 
     &     2*MTy*Cd(302)*Cr(5)*MStop2(3))/(Cd(47)*Cd(50))

        Ce(899) = -(Cd(10112)*Ce(898)) + 
     &   Cd(9950)*Cd(10103)*
     &    (T134(MGl2,MTy2,MStop2(1)) + T134(MGl2,MTy2,MStop2(2)))

        Ce(913) = AtC**2*UCStop(1,3)**2 + 
     &   MUE2*(3*UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(917) = Cd(47)*Cd(50)*
     &    (Cd(50)*Cd(164)*Cd(9947) + Cd(47)*Cd(251)*Cd(9949)) + 
     &   MUE2*Cd(9980)*UStop2(1,1)**2

        Ce(923) = -(AtC**2*UCStop(1,3)**2) + 
     &   MUE2*(-Cd(10658) + Cd(9948)*UStop2(1,1)**2)

        Ce(925) = Cd(9948)*Ce(924) + 
     &   AtC**2*Cr(85)*UCStop(1,3)**2 + MUE2*Cr(28)*UStop2(1,1)

        Ce(929) = (2*(-4*Cd(47)**2*Cd(50)**2*Cd(10014)*
     &        Cd(10085) + MStop2(3)*
     &        (-2*Cd(47)*Cd(50)*Cd(10022) + 
     &          Cd(10031)*Cd(10086)*MStop2(3)))*UCStop(1,3))/
     &   (MZ2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3)**2)

	Ce(931) = At*(M_3C*At + 2*M_3*AtC) - 2*M_3*MUE2

        Ce(941) = Cd(164)*Cd(10147) + Ce(940)*MGlpTmSt2(1) + 
     &   AtC*Cd(10026)*UCStop(1,3)

        Ce(942) = Cd(251)*Cd(10147) + Ce(940)*MGlpTmSt2(2) + 
     &   AtC*Cd(10151)*UCStop(1,3)

	Ce(95) = -8*MTy*S2A + 6*C2A*Cd(13) - 6*S2A*Cd(585)

        Ce(951) = (2*(-4*Cd(47)**2*Cd(50)**2*Cd(10014)*
     &        Cd(10094) + MStop2(3)*
     &        (-2*Cd(47)*Cd(50)*Cd(10044) + 
     &          Cd(10031)*Cd(10039)*MStop2(3)))*UCStopC(1,3))/
     &   (MZ2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3)**2)

	Ce(966) = M_3C*At + M_3*AtC + 2*MStop2(1)

	Ce(98) = 8*MTy*S2A + 6*C2A*Cd(13) - 6*S2A*Cd(585)

        Ce(977) = 2*MGlmT2*MStop2(3)*UStop2(1,1)*UStop2(1,2) + 
     &   Cd(47)*Cd(9969)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(989) = 12*Cd(47) + 4*Cd(152)*Cd(164) + 
     &   At*AtC*Cd(164)*UStop2(1,1)

	Ce(990) = 5 + UStop2(1,2)*(-6 + 2*UStop2(1,2))

        Ce(1028) = MStop2(3)*
     &    (2*M_3C*MGlmT2*Cd(47)*UCStop(1,3)**2 + 
     &      Cd(9948)*Ce(1027)*UStop2(1,1)) + 
     &   At*(Cd(47)**2*UStop2(1,1)**2 + Cd(10287)*UStop2(1,2)**2)

        Ce(1029) = 4*Ce(1025) + 
     &   AtC*(4*Ce(1028) - 2*AtC*Ce(1026)*UCStop(1,3)**2) - 
     &   2*At*Cd(10292)*UCStopC(1,3)**2

        Ce(1044) = -(Cd(10280)*Ce(1043)) + 
     &   AtC*(Cd(10317) + AtC*Cd(10086)*Cd(10280)*UCStop(1,3)) - 
     &   2*At*Cd(47)**2*UCStopC(1,3)

        Ce(1048) = M_3C*Ce(944) + 
     &   AtC*(-2*Cd(10324) + AtC*Cd(10088)*UCStop(1,3)) + 
     &   Ce(1047)*UCStopC(1,3)

	Ce(1078) = MTy*MGlpTmSt2(1) - 2*M_3*MTy2*UCStopC(1,3)

        Ce(1099) = -(MUE2*Cd(9966)) + 
     &   AtC**2*Cd(10077)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(9969)*Cr(28) + At*Cd(10077)*UCStopC(1,3)**2)

        Ce(111) = 2*C2A*Cd(426) + S2A*(-Cd(425) + Cd(427)) + 
     &   2*Ce(110)*UStop2(1,1)*UStop2(1,2)

        Ce(1119) = -(Cd(12371)*Ce(1043)) + 
     &   AtC*(2*Ce(1118) + AtC*Cd(10086)*Cd(12371)*UCStop(1,3)) + 
     &   4*At*Cd(47)**2*UCStopC(1,3)

        Ce(1126) = MGlmT2*Ce(946) + 
     &   AtC*(-2*Ce(1125) + AtC*MGlmT2*Cd(10088)*UCStop(1,3)) + 
     &   2*At*Cd(47)*Cd(9955)*UCStopC(1,3)*UStop2(1,1)

        Ce(1162) = -(Cd(11740)*UStop2(1,1)) + 
     &   4*(UStop2(1,2)**3 + 
     &      UStop2(1,1)**2*(UStop2(1,1) + 3*UStop2(1,2)))

        Ce(1165) = 5*MStop2(3) - 
     &   AtC*(At*Cr(28) + AtC*UCStop(1,3)**2) - 
     &   At**2*UCStopC(1,3)**2

        Ce(1166) = MUE2*Cr(28) + MUE**2*UCStop(1,3)**2 + 
     &   MUEC**2*UCStopC(1,3)**2 - 12*MStop2(3)*UStop2(1,1)**2

        Ce(1197) = Cd(9948)*Ce(1196) + 
     &   At*(AtC*Cd(9970) + At*UCStopC(1,3)**2) + 
     &   Ce(1195)*UStop2(1,1)

        Ce(1199) = -(MUE2*Cd(9966)) + 
     &   AtC**2*Cd(10077)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(9969)*Cr(28) - At*Cd(9972)*UCStopC(1,3)**2)

        Ce(1203) = 2*Cd(47)*Cd(50)*Cd(6096)*Cd(9948) - 
     &   Cd(47)*(Cd(47)*Ce(1202) + 2*Cd(50)**2*Cr(85)) - 
     &   Cd(50)**2*MGlpTmSt2(1)*MStop2(3)

        Ce(1213) = UStop2(1,1)*
     &    (Ce(1212) + 4*UStop2(1,1)*(-3 + 2*UStop2(1,1))) - 
     &   4*Cd(1052)*UStop2(1,2)**2

        Ce(1215) = UStop2(1,1)*
     &    (Ce(1212) + 4*UStop2(1,1)*(-3*Cd(787) + 2*UStop2(1,1)))-
     &     4*Cd(1052)*UStop2(1,2)**2

        Ce(1262) = -(MUE2*Ce(1261)) - 
     &   AtC**2*Cd(10610)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(10609) - At*Cd(10611)*UCStopC(1,3)**2)

        Ce(1264) = 2*M_3C*UCStop(1,3)**2 + 
     &   Ce(1263)*UStop2(1,1) + 
     &   At*UStop2(1,2)*(-UStop2(1,1)**2 + Cd(10615)*UStop2(1,2))

        Ce(1266) = -MUE2 + AtC*Ce(1264) + 
     &   At*(Ce(1265)*UCStopC(1,3)**2 + 
     &      2*M_3C*Cd(9948)*UStop2(1,1))

        Ce(1271) = At*Cd(11852) + AtC*Cd(11850)*UCStop(1,3)**2 - 
     &   MUE2*Ce(1268)*UStop2(1,1)**2 + Ce(1270)*UStop2(1,2)

        Ce(1314) = -(Cd(10604)*Ce(1043)) + 
     &   AtC*(Cd(10636) + AtC*Cd(10086)*Cd(10604)*UCStop(1,3)) - 
     &   2*At*Cd(50)**2*UCStopC(1,3)

        Ce(1339) = MTy2 + 
     &   2*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3)) + 
     &   MGl2*UStop2(1,2)

	Ce(1341) = MTy*MGlpTmSt2(2) + 2*M_3*MTy2*UCStopC(1,3)

        Ce(1349) = Ce(1410) + AtC**2*UCStop(1,3)**2 + 
     &   MUE2*UStop2(1,2)**2

        Ce(1359) = MStop2(2)*UStop2(1,1)**3 + 
     &   MStop2(1)*UStop2(1,2)**3 + 
     &   UStop2(1,1)*UStop2(1,2)*
     &    (Cd(814)*UStop2(1,1) + Cd(652)*UStop2(1,2))

        Ce(1362) = -(MUE2*Ce(1359)) + Cd(9948)*Ce(1361) + 
     &   AtC*(-(Ce(1360)*UCStop(1,3)**2) + 
     &      2*M_3*MStop2(3)*UStop2(1,1)*UStop2(1,2))

        Ce(1366) = -(Cd(12483)*Ce(1043)) + 
     &   AtC*(2*Cd(12511) + AtC*Cd(10086)*Cd(12483)*UCStop(1,3)) + 
     &   4*At*Cd(50)**2*UCStopC(1,3)

        Ce(1368) = Ce(949) + 
     &   AtC*(-2*At*Cd(10092)*Cr(5) + AtC*Cd(10087)*UCStop(1,3))

        Ce(1389) = 8*M_3*UCStopC(1,3)**2 + 
     &   12*M_3C*UStop2(1,1)**2 + AtC*UStop2(1,1)**2 + 
     &   4*M_3C*UStop2(1,2)**2

        Ce(1405) = MStop2(2)*UStop2(1,1)**3 + 
     &   MStop2(1)*UStop2(1,2)**3 + 
     &   UStop2(1,1)*UStop2(1,2)*
     &    (Cd(778)*UStop2(1,1) + Cd(652)*UStop2(1,2))

        Ce(1406) = Cd(9948)*Ce(1361) - MUE2*Ce(1405) + 
     &   AtC*(-(Ce(1360)*UCStop(1,3)**2) + 
     &      2*M_3*MStop2(3)*UStop2(1,1)*UStop2(1,2))

        Ce(1409) = Cd(50)*Ce(1404) + MGlmT2*Ce(1406) + 
     &   At*(Ce(1408) - At*Cd(11855)*UCStopC(1,3)**2)

        Ce(1455) = 2*MTy*Cd(164)*Cr(5)*MStop2(3) + 
     &   Ce(1446)*UUStopC(1,1) + Ce(1450)*UUStopC(1,2) + 
     &   2*Ce(1454)*UUStopC(1,3)

        Ce(1467) = 2*MTy*Cd(152)*Cr(5)*MStop2(3) + 
     &   Ce(1462)*UUStopC(1,1) + Ce(1464)*UUStopC(1,2) + 
     &   2*Ce(1466)*UUStopC(1,3)

        Ce(1482) = 2*MTy*Cd(251)*Cr(5)*MStop2(3) - 
     &   Ce(1473)*UUStopC(2,1) - Ce(1477)*UUStopC(2,2) - 
     &   2*Ce(1481)*UUStopC(2,3)

        Ce(1497) = 2*MTy*Cd(179)*Cr(5)*MStop2(3) + 
     &   Ce(1490)*UUStopC(2,1) + Ce(1493)*UUStopC(2,2) + 
     &   2*Ce(1496)*UUStopC(2,3)

        Ce(1522) = 4*MTy2*Ce(1521) + 
     &   MUE**2*UUStop(1,1)*UUStop(2,1) - 
     &   Ce(186)*Ce(187)*UUStop(1,2)*UUStop(2,2) - 
     &   2*MUE2*UUStop(3,1)*UUStop(3,2)

        Ce(1527) = Ce(1526)*UCStop(3,1) - 
     &   MUEC**2*Cr(91)*UCStop(3,4) - 
     &   AtC*Cd(9794)*UCStop(1,3)*UStop(1,1) + 
     &   Cd(698)*Cd(9794)*UStopC(2,1)

        Ce(1529) = -(MUE**2*Cd(11106)*UUStop(1,1)) + 
     &   AtC*(Cd(10925)*Cd(10980)*UStop(1,2) + 
     &      AtC*Cd(11106)*UUStop(1,1)) + Ce(1528)*UUStop(1,2)

        Ce(1530) = 2*At*MTy*Cr(91) + 8*MTy2*UUStop(3,1) + 
     &   (-At**2 + MUEC**2)*UUStop(3,2)

        Ce(1536) = Ce(35)*Ce(36)*UStop(1,1)*UUStop(2,1) - 
     &   2*At*MTy*UStop(1,2)*UUStop(2,2) + 
     &   AtC*Cd(10980)*UUStop(2,3) + MUE2*UStop(2,1)*UUStop(3,2)

        Ce(1544) = Ce(35)*Ce(36)*UStop(2,1)*UUStop(1,1) - 
     &   2*At*MTy*UStop(2,2)*UUStop(1,2) + 
     &   AtC*Cd(11017)*UUStop(1,3) + MUE2*UStop(1,1)*UUStop(3,2)

        Ce(1604) = -(MUE**2*Cd(11657)) + 
     &   16*MTy2*Cd(3852)*UUStop(1,2)*UUStop(2,2) + 
     &   AtC*(AtC*Cd(11657) + 8*MTy*Cd(3852)*Cr(91)*UUStop(3,2))

        Ce(1634) = Ce(1633)*UStop2(1,1) + 
     &   (UStop2(1,1)**2 - 2*UStop2(1,2)**2)*UUStop(1,3) + 
     &   UCStop(3,1)*UUStop(3,4)

        Ce(1647) = UUStop(1,1)*
     &    (4*AtC*Cd(10976)*UStop(1,2) + 
     &      (-AtC**2 + MUE**2)*UUStop(1,1)) + Ce(1646)*UUStop(1,2)

        Ce(1687) = -(MUEC**2*Ce(1683)) + Cd(9794)*Ce(1685) + 
     &   2*Ce(1686)*MStop2(2) - 2*Ce(1643)*MStop2(3)

        Ce(1754) = 1 - UStop2(1,2)**2 + 
     &   UStop2(1,1)*(-UStop2(1,1) + UStop2(1,2))

        Ce(1773) = -(Ce(1772)*UStop2(1,1)) + 5*UStop2(1,1)**3 - 
     &   (Cd(1143) + 5*UStop2(1,1)**2)*UStop2(1,2)

        Ce(1775) = -(Ce(1774)*UStop2(1,1)) + 5*UStop2(1,1)**3 - 
     &   (Cd(1143) - 3*UStop2(1,1)**2)*UStop2(1,2)

	Ce(1777) = 12 - 5*(UStop2(1,1) - UStop2(1,2))**2

	Ce(1778) = 8 - (UStop2(1,1) - UStop2(1,2))*(1 + 2*UStop2(1,2))

        Ce(1779) = 10 - UStop2(1,2)**2 + 
     &   UStop2(1,1)*(-3*UStop2(1,1) + 4*UStop2(1,2))

        Ce(1782) = -3*Ce(1770)*UStop2(1,1) + 
     &   UStop2(1,1)**2*(1 + 4*UStop2(1,1)) - Ce(1781)*UStop2(1,2)

        Ce(1786) = Ce(1676) + 
     &   UStop2(1,1)*(-5*UStop2(1,1) + 6*UStop2(1,2))

	Ce(1792) = Ce(1836) - UStop2(1,2)**2

        Ce(180) = 2*Cd(286)*UStop2(1,1)**3 - 
     &   MStop2(2)*UStop2(1,2)**4 + 
     &   UStop2(1,1)*(-(Cd(617)*UStop2(1,1)**3) + 
     &      2*Cd(604)*UStop2(1,2)**2)

        Ce(1818) = 2*Ce(1817) - MUE2*Ce(1812)*UStop(2,2) + 
     &   AtC*(Cd(11352) + 2*AtC*MStop2(4)*UCStop(1,3)*UStop2(1,1)*
     &       UUStop(3,1))

        Ce(1831) = Cd(2366)*UStop2(1,1) - UStop2(1,1)**3 + 
     &   (Cd(970) + UStop2(1,1)**2)*UStop2(1,2)

        Ce(1833) = -2*Cd(6269)*UStop2(1,1) + 
     &   UStop2(1,1)**2*(1 + UStop2(1,1)) - Ce(1832)*UStop2(1,2)

        Ce(1834) = Cd(6269)*UStop2(1,1) + 
     &   UStop2(1,1)**2*(-1 - UStop2(1,2)) + Cd(2384)*UStop2(1,2)

	Ce(1864) = -2*Ce(1768) - Ce(1823)

	Ce(1867) = 11 - 4*(UStop2(1,1) - UStop2(1,2))**2

        Ce(1877) = -(Ce(1876)*UStop2(1,1)) + 4*UStop2(1,1)**3 - 
     &   (Cd(1054) + 4*UStop2(1,1)**2)*UStop2(1,2)

        Ce(1881) = -(Cd(1040)*UStop2(1,1)) + 4*UStop2(1,1)**3 - 
     &   (Ce(1880) - 5*UStop2(1,1)**2)*UStop2(1,2)

        Ce(1884) = 11 - UStop2(1,2)**2 + 
     &   UStop2(1,1)*(-4*UStop2(1,1) + 5*UStop2(1,2))

        Ce(1887) = Ce(1886)*UStop2(1,1) - 5*UStop2(1,1)**3 + 
     &   (5 + 4*UStop2(1,1)**2)*UStop2(1,2)

        Ce(1890) = -(Ce(1870)*UStop2(1,1)) - 
     &   Cd(6270)*UStop2(1,2) + UStop2(1,1)**2*(1 + UStop2(1,2))

        Ce(193) = AtC*(MUEC*Ce(192) - 
     &      4*MUE*Cd(609)*UCStop(1,3)**2) + 
     &   At*(2*MUE*Cd(606) - 4*MUEC*Cd(611)*UCStopC(1,3)**2)

        Ce(257) = CA*Ce(255) + C2A*MTy*MUE*Cr(91) - 
     &   1/2.D0*(AtC*(2*Ce(256) + AtC*S2A*UUStop(3,1)))

        Ce(260) = C2A*MUE*Cd(698) - 4*S2A*Cd(11122) + 
     &   2*MUEC*Ce(258) + AtC*Ce(259)

        Ce(31) = MUE**2*Ce(30)*UCStop(1,3)**2*UStop2(1,2) + 
     &   At*UCStopC(1,3)*UStop2(1,1)*
     &    (At*Cd(535)*UCStopC(1,3)*UStop2(1,1) - 
     &      2*Ce(26)*UStopC(1,1)) - Ce(29)*UStopC(1,2)*UUStopC(1,1)

        Ce(440) = (Cd(50)*Ce(438)*MStop2(3) + 
     &     Cd(47)*(-(Cd(9952)*Cd(9988)*MStop2(2)*MStop2(3)) + 
     &        Cd(50)**2*(Ce(439) + Ce(436)*MStop2(1)*MStop2(3))))/
     &   (Cd(47)*Cd(50)**2*MStop2(3))

	Ce(453) = -4*M_3C + AtC*(1 - Cd(9948) + UStop2(1,2))

	Ce(465) = -4*M_3 + At*(1 - Cd(9948) + UStop2(1,2))

        Ce(479) = 2*Cd(9948)*MStop2(2)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(-2*MStop2(1)*UStop2(1,1)**2 + 
     &      UStop2(1,2)*(Ce(478) - Cd(256)*UStop2(1,2)))

        Ce(482) = At*(MUEC*Ce(481)*UCStopC(1,3)**2 + 
     &      MUE*Ce(479)*UStop2(1,1)) + 
     &   AtC*(MUE*Ce(481)*UCStop(1,3)**2 + 
     &      MUEC*Ce(479)*UStop2(1,1))

        Ce(493) = 2*UStop2(1,1)**3 + 
     &   2*(Ce(492) - UStop2(1,1))*UStop2(1,1)*UStop2(1,2) + 
     &   Cd(10904)*UStop2(1,2)**2

        Ce(505) = MUE*(M_3*Cr(28) - 2*M_3C*UCStop(1,3)**2) + 
     &   MUEC*(M_3C*Cr(28) - 2*M_3*UCStopC(1,3)**2)

        Ce(524) = MStop2(1) + 
     &   2*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3)) - 
     &   MGl2*UStop2(1,2)

        Ce(597) = Ce(634) + 
     &   M_3*MUE*UStop2(1,1)*UStop2(1,2)**2 + 
     &   M_3C*MUEC*UStop2(1,1)*UStop2(1,2)**2

        Ce(641) = 2*MGlpT2*MTy2*Cd(13) - 2*MTy2*Ce(640) + 
     &   Cd(13)*Cd(152)*MGlpTmSt2(1)

        Ce(677) = Cd(290)*UStop2(1,1)**2 + UStop2(1,2)**2 + 
     &   UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

	Ce(7) = 4*MTy*Ce(6) + Ce(12) + 2*Ce(1)*UStop2(1,1)*UStop2(1,2)

        Ce(77) = -(At*S2A*Cd(423)) + S2A*Ce(76)*UCStop(1,3) + 
     &   C2A*(2*MUEC*Cd(577) - 2*MUE*MStop2(4)*UCStop(1,3)**2)

        Ce(784) = MGl2 + 6*M_3C*MTy*UCStop(1,3) + 
     &   (MTy2 - MStop2(1))*UStop2(1,1)

        Ce(795) = -2*Cr(85)*UStop2(1,1)**2 - 2*UStop2(1,2)**2 + 
     &   UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

	Ce(839) = Ce(833) - Cd(11049)*MStop2(3)

        Ce(855) = MGlmT2 - MStop2(2) + 
     &   MTy2*(-2*Cd(9948) + 4*UStop2(1,1))

        Ce(905) = -(Ce(904)*MStop2(1)) + 
     &   ((2*Cd(47)*Cd(50)**2*Cd(10014) + 
     &        (2*Cd(50)**2*Cd(9947)*MGlpTmSt2(1) - 
     &           Cd(47)*Cd(9952)*MStop2(2))*MStop2(3))*
     &      UCStop(1,3)**2)/(MZ2*SW2*Cd(47)*Cd(50)**2*MStop2(3))

        Ce(910) = -Ce(1838) + 
     &   AtC*(AtC*UCStop(1,3)**2 + 2*At*Cd(9981)*UStop2(1,1))

        Ce(922) = Ce(921)*MStop2(1) + 
     &   ((-2*Cd(47)*Cd(50)**2*Cd(10014) + 
     &        (-2*Cd(50)*Ce(920) + Cd(47)*Cd(9952)*MStop2(2))*
     &         MStop2(3))*UCStopC(1,3)**2)/
     &    (MZ2*SW2*Cd(47)*Cd(50)**2*MStop2(3))

        Ce(926) = -(Cd(50)*Cd(9947)*Ce(925)*MGlpTmSt2(1)) + 
     &   Cd(9949)*Cd(9950)*Ce(923)*MGlpTmSt2(2) + 
     &   At*(AtC*Cd(10071) + 
     &      At*Cd(10076)*Cd(10077)*UCStopC(1,3)**2)

        Ce(932) = -3*M_3C*AtC**2*UCStop(1,3)**3 - 
     &   3*M_3*At**2*UCStopC(1,3)**3 + 
     &   Cd(10016)*(Ce(930)*UCStop(1,3) + Ce(931)*UCStopC(1,3))

        Ce(950) = 2*At*AtC*Cd(10093)*Cr(5) + 
     &   Cd(50)*Cd(9947)*Ce(946)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Ce(949)*MGlpTmSt2(2) + 
     &   AtC**2*Cd(10089)*UCStop(1,3)

	Ce(956) = 2*Cd(9948) - Ce(996)

        Ce(965) = -Ce(1838) + 
     &   AtC*(AtC*UCStop(1,3)**2 + 2*At*UStop2(1,1)*UStop2(1,2))

        Ce(976) = -(MUE2*Cr(28)) - 
     &   AtC**2*Cd(9972)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(9969)*Cr(28) + At*Cd(10077)*UCStopC(1,3)**2)

        Ce(980) = 2*MGlmT2*MStop2(3)*
     &    (M_3C*UCStop(1,3)**2 + Cd(9948)*Cd(9962)*UStop2(1,1))+
     &     At*(Cd(47)*UStop2(1,1)**2 + Ce(979)*UStop2(1,2)**2)

        Ce(982) = -(Cd(9972)*(Cd(47) + MStop2(1)**2)) + 
     &   MGlmT2*MStop2(2) - MStop2(1)*MStop2(3)

        Ce(999) = MUE2*Cd(9966) + 
     &   AtC*(4*Ce(998) + AtC*Cd(9972)*UCStop(1,3)**2) + 
     &   At*(-Cd(10481) + At*Cd(9972)*UCStopC(1,3)**2)

        Ce(1017) = 4*MTy2 - Ce(1838) + 4*MStop2(1) - 
     &   MUEC**2*UCStopC(1,3)**2

        Ce(1030) = AtC*Cd(9967)*UCStop(1,3)**2 + 
     &   At*(2*Cd(11727) + At*UCStopC(1,3)**2) - 
     &   2*MUE2*Cd(9948)*UStop2(1,1)**2 + 
     &   2*Ce(974)*UStop2(1,1)*UStop2(1,2)

        Ce(1062) = -(MUE2*Cd(9966)) + 4*AtC*Cd(10169) + 
     &   AtC**2*Cd(10077)*UCStop(1,3)**2 + 
     &   At*(Cd(10347) + At*Cd(10077)*UCStopC(1,3)**2)

        Ce(1085) = At*Cd(11720)*Cr(28) + Ce(1084)*Cr(28) - 
     &   AtC*Cd(11721)*UCStop(1,3)**2 - 
     &   At*(At*MGlmT2 + 8*M_3*Cd(118))*UCStopC(1,3)**2

        Ce(1093) = AtC*Ce(1090)*UCStop(1,3)**2 + 
     &   At*(2*Ce(1092) + At*Cd(12371)*UCStopC(1,3)**2) - 
     &   2*MUE2*Cd(9948)*Cd(12371)*UStop2(1,1)**2 + 
     &   2*Ce(1089)*UStop2(1,1)*UStop2(1,2)

	Ce(1102) = 5*Cd(9948) - Ce(1162)

        Ce(1107) = MUE2*Cd(9966) + 
     &   AtC*(8*Cd(10169) + AtC*Cd(9972)*UCStop(1,3)**2) + 
     &   At*(-Cd(11745) + At*Cd(9972)*UCStopC(1,3)**2)

        Ce(1183) = AtC*(Cd(12653) - 
     &      AtC*Cd(11962)*UCStop(1,3)**2) + 
     &   At*(8*Ce(1182) - At*Cd(11962)*UCStopC(1,3)**2) + 
     &   2*MUE2*Cd(11962)*UStop2(1,1)*UStop2(1,2)

        Ce(1193) = 2*At*Cd(12002)*Cd(12005) - 
     &   MUE2*Cd(47)*Cd(50)*Ce(1192) + 
     &   AtC*(Cd(12003) + AtC*Cd(47)*Cd(50)**2*UCStop(1,3)**2) + 
     &   At**2*Cd(47)*Cd(50)*Cd(12004)*UCStopC(1,3)**2

        Ce(1222) = MUE2*Cd(11978) + 2*AtC*Ce(1221) - 
     &   AtC**2*Cd(11984)*UCStop(1,3)**2 + 
     &   At*(Cd(12039) + At*Cd(11985)*UCStopC(1,3)**2)

        Ce(1247) = 4*MTy2 - 4*AtC*MTy*UCStop(1,3) - 
     &   MUE**2*UCStop(1,3)**2 - 
     &   UCStopC(1,3)*(4*At*MTy + MUEC**2*UCStopC(1,3))

        Ce(1249) = 4*UCStopC(1,3)*(MTy + 2*M_3*UCStopC(1,3)) + 
     &   2*UStop2(1,1)*(-(AtC*Cd(9981)) + 
     &      4*M_3C*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Ce(1257) = Ce(923) + 
     &   At*(AtC*Cd(10070) + At*Cd(10077)*UCStopC(1,3)**2)

	Ce(1274) = -Ce(1368) - 2*Cd(585)*Cd(9955)*Cr(75)*MStop2(1)

        Ce(1294) = MUE2*Cd(10070) + 4*AtC*Ce(1293) + 
     &   AtC**2*Cd(9969)*UCStop(1,3)**2 + 
     &   At*(Cd(10653) + At*UCStopC(1,3)**2)

        Ce(1300) = At*Cd(11852) - 2*Cd(9948)*Ce(1299) + 
     &   AtC*Cd(11850)*UCStop(1,3)**2 - Ce(1298)*UStop2(1,1) + 
     &   Ce(1270)*UStop2(1,2)

        Ce(1363) = 2*At*Cd(12508) + 
     &   4*Cd(50)*(Cd(50)*Ce(1355) + MGlmT2*Ce(1362)) - 
     &   Ce(1358)*MStop2(3) + At**2*Cd(12505)*UCStopC(1,3)**2

        Ce(1377) = Ce(913) + MUEC**2*UCStopC(1,3)**2 + 
     &   At*AtC*UStop2(1,1)**2

        Ce(1417) = Ce(701) + 
     &   AtC*(M_3C*UCStop(1,3)**2 + 
     &      M_3*UStop2(1,1)*UStop2(1,2))

        Ce(1436) = MUE2*Cd(10070) + 4*AtC*Ce(1435) + 
     &   AtC**2*Cd(9969)*UCStop(1,3)**2 + 
     &   At*(Cd(11883) + At*UCStopC(1,3)**2)

        Ce(1531) = Ce(1530)*UUStop(1,2) + 
     &   2*(AtC*Cd(11166) - MUE2*Cr(91))*UUStop(1,3) + 
     &   (-AtC**2 + MUE**2)*UUStop(1,1)*UUStop(3,1)

        Ce(1553) = AtC*(Ce(1552) - AtC*Cd(11080)*UCStop(1,3)) + 
     &   Ce(1547)*UCStopC(1,3)*UStopC(2,1) + 
     &   Ce(1548)*UUStopC(2,2) + Ce(1550)*UUStopC(2,3)

        Ce(1582) = -(At*(At*Ce(1579) + 2*MTy*Ce(1581))) + 
     &   Ce(1577)*UUStop(1,3) + Ce(1578)*UUStop(2,3) - 
     &   4*MTy2*Cd(523)*Cr(58)*UUStop(3,1)

        Ce(1759) = Ce(1530)*UUStop(2,2) + 
     &   2*(AtC*Cd(11166) - MUE2*Cr(91))*UUStop(2,3) + 
     &   (-AtC**2 + MUE**2)*UUStop(2,1)*UUStop(3,1)

        Ce(1765) = 6 + 3*Cd(1428) + 
     &   (4 - 6*UStop2(1,1))*UStop2(1,1) - 2*UStop2(1,2)**2

        Ce(1771) = 3*Ce(1770)*UStop2(1,1) - 5*UStop2(1,1)**3 + 
     &   (12 + UStop2(1,1)**2)*UStop2(1,2) - 5*UStop2(1,2)**3

        Ce(1784) = Ce(1783)*UStop2(1,1) - 5*UStop2(1,1)**3 + 
     &   (13 + 2*UStop2(1,1)**2)*UStop2(1,2) - 4*UStop2(1,2)**3

        Ce(1826) = -(Cd(6269)*UStop2(1,1)**2) + UStop2(1,1)**4 + 
     &   2*UStop2(1,1)*(-Cd(2366) + UStop2(1,1)**2)*UStop2(1,2) - 
     &   Cd(970)*UStop2(1,2)**2

        Ce(1827) = -(Cd(2359)*UStop2(1,1)**2) + UStop2(1,1)**4 + 
     &   2*UStop2(1,1)*(-Cd(6269) + 2*UStop2(1,1)**2)*
     &    UStop2(1,2) - Cd(970)*UStop2(1,2)**2

        Ce(1871) = -(Cd(1040)*UStop2(1,1)**2) + 
     &   4*UStop2(1,1)**4 + 
     &   2*UStop2(1,1)*(-Ce(1870) + 3*UStop2(1,1)**2)*
     &    UStop2(1,2) - Cd(1054)*UStop2(1,2)**2

        Ce(1879) = 2*Ce(1878)*UStop2(1,1) - 8*UStop2(1,1)**3 - 
     &   (-27 + UStop2(1,1)**2)*UStop2(1,2) - 5*UStop2(1,2)**3

        Ce(230) = -4*Ce(228) - 2*Ce(229)*UCStop(1,3)**2 - 
     &   Cd(639)*Ce(186)*Ce(187)*UCStopC(1,3)**2 + 
     &   At*(2*AtC*Ce(214) - 4*MTy*Ce(227)*UCStopC(1,3))

        Ce(246) = 4*MTy2*Cd(492) - 4*MTy*Cd(585) + Ce(12) + 
     &   2*Ce(1)*UStop2(1,1)*UStop2(1,2)

        Ce(279) = MUE**2*S2A*UUStop(1,1) + 
     &   AtC*(2*Ce(248)*UStop(1,1) - AtC*S2A*UUStop(1,1)) - 
     &   4*MTy2*S2A*UUStop(1,2) + 4*C2A*MTy*MUE*UUStop(1,3)

        Ce(300) = MUE**2*S2A*UUStop(2,1) + 
     &   AtC*(2*Ce(284)*UStop(2,1) - AtC*S2A*UUStop(2,1)) - 
     &   4*MTy2*S2A*UUStop(2,2) + 4*C2A*MTy*MUE*UUStop(2,3)

        Ce(37) = At*AtC*Cr(28) - MUE2*Cr(28) - 
     &   Ce(35)*Ce(36)*UCStop(1,3)**2 + 
     &   (-At**2 + MUEC**2)*UCStopC(1,3)**2

        Ce(371) = Cd(9965)/(Cd(47)*Cd(50)) - 
     &   (Cd(9952)*Cd(9957)*MStop2(2))/Cd(50)**2 + 
     &   Ce(370)/(Cd(50)*Cd(9950)*MStop2(3)) + 
     &   (MStop2(1)*(-(At*Cd(47)**2*Cd(9955)*Cd(9956)) - 
     &        Cd(9951)*Cd(9957)*MStop2(3)**2))/
     &    (Cd(47)**2*MStop2(3)**2)

        Ce(380) = Cd(10000)/(Cd(47)*Cd(50)) - 
     &   (Cd(9952)*Cd(9993)*MStop2(2))/Cd(50)**2 + 
     &   Ce(379)/(Cd(50)*Cd(9950)*MStop2(3)) + 
     &   (MStop2(1)*(-(AtC*Cd(47)**2*Cd(9955)*Cd(9956)) - 
     &        Cd(9951)*Cd(9993)*MStop2(3)**2))/
     &    (Cd(47)**2*MStop2(3)**2)

        Ce(382) = (4*Cd(10017)*Ce(381) + 
     &     (MStop2(3)*(-2*Cd(47)*Cd(50)*Cd(10024) + 
     &          MStop2(3)*(Cd(10034) + 
     &             (-(Cd(50)**2*Cd(10012)*MStop2(1)) + 
     &                Cd(47)**2*Cd(10013)*MStop2(2))*UCStop(1,3))))
     &       /(Cd(47)**2*Cd(50)**2))/MStop2(3)**2

        Ce(383) = (4*Cd(10036)*Ce(381) + 
     &     (MStop2(3)*(-2*Cd(47)*Cd(50)*Cd(10046) + 
     &          MStop2(3)*(Cd(10041) + 
     &             (-(Cd(50)**2*Cd(10012)*MStop2(1)) + 
     &                Cd(47)**2*Cd(10013)*MStop2(2))*UCStopC(1,3)))
     &        )/(Cd(47)**2*Cd(50)**2))/MStop2(3)**2

        Ce(445) = (8*Cd(10017)*Ce(381))/MStop2(3)**2 + 
     &   (4*Ce(443))/(Cd(47)*Cd(50)*MStop2(3)) + 
     &   (2*Cd(10013)*MStop2(2)*UCStop(1,3))/Cd(50)**2 + 
     &   (2*Ce(444) - 2*Cd(50)**2*Cd(10134)*MStop2(1)*UCStop(1,3))/
     &    (Cd(47)**2*Cd(50)**2)

        Ce(771) = MGl2*Cd(9983) + 
     &   (-MTy2 + MStop2(2) - 5*MStop2(3))*UStop2(1,1)**2 - 
     &   MStop2(3)*UStop2(1,2)**2

        Ce(786) = MStop2(1) - 6*M_3C*MTy*UCStop(1,3) - 
     &   MGlpT2*UStop2(1,1) - MStop2(3)*UStop2(1,2)

        Ce(907) = Ce(906)*Cr(28) - 
     &   AtC*Cd(10110)*UCStop(1,3)**2 + 
     &   At*(Cd(9967)*Cr(28) - 4*M_3*UCStopC(1,3)**2 - 
     &      At*UCStopC(1,3)**2)

        Ce(94) = -(At*AtC*S2A*Cd(435)) + MUE2*S2A*Cd(435) + 
     &   C2A*Cd(424)*Cd(435) + 2*S2A*MStop2(1)*MStop2(3) - 
     &   Ce(93)*MStop2(4)

        Ce(938) = AtC**2*Cd(10032) - Ce(1043) + 
     &   2*At*AtC*Cr(5)*UStop2(1,1)*UStop2(1,2)

        Ce(962) = MUE2*Cd(9984) - 4*At*Ce(961) + 
     &   AtC*(Cd(10175) + AtC*Cd(10081)*UCStop(1,3)**2) + 
     &   At**2*Cd(10125)*UCStopC(1,3)**2

        Ce(97) = -(At*AtC*S2A*Cd(423)) + MUE2*S2A*Cd(423) + 
     &   2*C2A*Cd(424)*Cd(577) - 2*S2A*MStop2(2)*MStop2(3) - 
     &   Ce(93)*MStop2(4)

        Ce(969) = 2*UCStopC(1,3)*(MTy - 2*M_3*UCStopC(1,3)) + 
     &   UStop2(1,1)*(AtC*Cd(9981) - 
     &      4*M_3C*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Ce(981) = -(MUE2*Ce(977)) + AtC*Ce(980) - 
     &   AtC**2*Ce(978)*UCStop(1,3)**2 + 
     &   At*(2*MGlmT2*Cd(12005)*MStop2(3) - 
     &      At*Ce(978)*UCStopC(1,3)**2)

        Ce(994) = MUE2*Cd(9966) + 
     &   AtC*(-4*Cd(10169) + AtC*Cd(9972)*UCStop(1,3)**2) + 
     &   At*(-Cd(10347) + At*Cd(9972)*UCStopC(1,3)**2)

        Ce(1016) = (At - MUEC)*(At + MUEC)*UCStopC(1,3)*
     &    UStop2(1,1)*UStop2(1,2) + 
     &   UCStop(1,3)*(4*MTy2 + 4*MStop2(1) - 
     &      MUE**2*UCStop(1,3)**2 - 2*MUE2*UStop2(1,1)*UStop2(1,2))

        Ce(1141) = 4*Cd(152) + AtC**2*UCStop(1,3)**2 - 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2)**2 + 
     &   At*(At*UCStopC(1,3)**2 + 2*AtC*UStop2(1,1)*UStop2(1,2))

        Ce(1146) = M_3*Ce(1145) + 
     &   M_3C*UCStop(1,3)*
     &    (MUE**2*UCStop(1,3)**2 + 2*MUE2*UStop2(1,1)*UStop2(1,2))+
     &     UCStopC(1,3)*(-2*At*MStop2(1) + 
     &      M_3C*(-At**2 + MUEC**2)*UStop2(1,1)*UStop2(1,2))

        Ce(1176) = -Ce(910) + MUE**2*UCStop(1,3)**2 + 
     &   MUEC**2*UCStopC(1,3)**2

        Ce(1184) = (-4*MGl2*Cd(12592) + Ce(1183) + 
     &     4*Cd(50)**2*Cd(12590)*MStop2(1) + 
     &     4*Cd(47)**2*Cd(12591)*MStop2(2) + 
     &     MUE**2*Cd(11962)*UCStop(1,3)**2 + 
     &     MUEC**2*Cd(11962)*UCStopC(1,3)**2)/(Cd(47)**2*Cd(50)**2)

        Ce(1291) = MUE2*Cd(10070) + 4*AtC*Ce(1290) + 
     &   AtC**2*Cd(9969)*UCStop(1,3)**2 + 
     &   At*(Cd(10795) + At*UCStopC(1,3)**2) + 
     &   12*Cd(9955)*MStop2(1)*UStop2(1,1)**2

	Ce(1325) = Ce(1294) + 12*Cd(9955)*MStop2(1)*UStop2(1,1)**2

        Ce(1352) = AtC*Cd(11844) + 4*M_3C*At*Cd(138)*Cr(28) + 
     &   MUE2*Cd(9970)*MStop2(1) - 
     &   AtC**2*MStop2(1)*UCStop(1,3)**2 - 
     &   At*(8*M_3*Cd(138) + At*Cd(9972)*MStop2(1))*
     &    UCStopC(1,3)**2

        Ce(1371) = Ce(1838) + 
     &   (-AtC**2 + MUE**2)*UCStop(1,3)**2 + 
     &   MUEC**2*UCStopC(1,3)**2 - 
     &   2*At*AtC*UStop2(1,1)*UStop2(1,2)**2

	Ce(1378) = -Ce(1377) + MStop2(3)

        Ce(168) = 2*Cd(427)*Cr(75) + 4*MTy*Cd(585)*Cr(75) + 
     &   AtC**2*Cd(38)*UCStop(1,3)**2 + 
     &   4*MUE2*Cr(75)*UStop2(1,1)*UStop2(1,2) + 
     &   At*(-2*At*Cr(75)*UCStopC(1,3)**2 - 
     &      4*AtC*Cr(75)*UStop2(1,1)*UStop2(1,2))

	Ce(170) = Cd(427) - 2*MTy*Cd(585) - Ce(965)

        Ce(1762) = -(MUE2*Ce(1761)*UStop(2,2)) + 
     &   MUE**2*UCStop(1,3)*UCStopC(3,1)*UUStop(1,1) + 
     &   AtC*(Cd(11496) - 
     &      AtC*UCStop(1,3)*UStop2(1,1)*UUStop(3,1)) + 
     &   2*At*MTy*Cd(10918)*UUStop(3,2) + 4*MTy2*UUStop(3,4)

        Ce(1869) = UStop2(1,1)*
     &    (UStop2(1,1)*(11 + 
     &         UStop2(1,1)*(-4*UStop2(1,1) - 12*UStop2(1,2))) + 
     &      2*Ce(1868)*UStop2(1,2)) + 
     &   UStop2(1,2)**2*(11 - 4*UStop2(1,2)**2)

        Ce(188) = At*AtC*Ce(180) - MUE2*Ce(180) + 
     &   (MUE**2*Cd(609) - AtC**2*Cd(623))*UCStop(1,3)**2 - 
     &   Cd(618)*Ce(186)*Ce(187)*UCStopC(1,3)**2

        Ce(768) = MGlpTmSt2(1) + 6*M_3*MTy*UCStopC(1,3) + 
     &   (MGl2 - MStop2(3))*UStop2(1,2)

        Ce(975) = -(AtC*Cd(9967)*UCStop(1,3)**2) - 
     &   At*(2*M_3 + At)*UCStopC(1,3)**2 + 
     &   2*MUE2*Cd(9948)*UStop2(1,1)**2 - 
     &   2*(At*Cd(9963) + Ce(974))*UStop2(1,1)*UStop2(1,2)

        Ce(1049) = -(At**2*UCStopC(1,3)*
     &      (UStop2(1,1)**2 + (Cd(9948) - UStop2(1,1))*UStop2(1,2))
     &      ) + 2*UCStop(1,3)*
     &    (MStop2(3) + MUE2*
     &       (-(Cd(9948)*UStop2(1,1)) + UStop2(1,1)**2 + 
     &         UStop2(1,2)**2))

        Ce(1179) = (2*(MTy2*Cd(12588) + 
     &        Cd(50)**2*Cd(585)*MGlpTmSt2(1)*MStop2(1) - 
     &        Cd(47)**2*Cd(585)*MGlpTmSt2(2)*MStop2(2)) + 
     &     Cd(12045)*(-2*MGl2*Cd(585) + Ce(938) - 
     &        MUE**2*Cd(10086)*UCStop(1,3) - 
     &        MUEC**2*Cd(10039)*UCStopC(1,3)))/
     &   (Cd(47)**2*Cd(50)**2)

	Ce(1207) = 4*Cd(152) + Opt(22) - 2*MUE2*UStop2(1,1)*UStop2(1,2)

	Ce(1208) = 4*Cd(179) + Opt(22) - 2*MUE2*UStop2(1,1)*UStop2(1,2)

	Ce(1238) = Cd(152)*Cd(10141) + AtC*Cd(12059) + Opt(24)

	Ce(1241) = Cd(179)*Cd(10141) + AtC*Cd(10781) + Opt(24)

        Ce(1283) = 4*Cd(179)*Cr(28) + Opt(22) - 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2)

        Ce(1373) = Opt(25) - MUE2*Cr(75)**2*UStop2(1,1) + 
     &   MUEC**2*UCStopC(1,3)**2*(1 + UStop2(1,2))

        Ce(1419) = -Opt(25) + MUE2*Cr(75)**2*UStop2(1,1) + 
     &   MUEC**2*UCStopC(1,3)**2*(-1 - UStop2(1,2))

        Ce(1501) = 8*MTy2*UUStop(1,1) + 
     &   2*MUE2*UStop2(1,2)*UUStop(1,1) - 
     &   MUEC**2*UCStop(3,1)*UUStop(3,2) + 
     &   At**2*(-(UStop2(1,1)*UUStop(1,2)) + 
     &      UCStop(3,1)*UUStop(3,2)) + 
     &   2*At*MTy*(UStop2(1,2)*UUStop(1,3) + 
     &      UCStop(3,1)*UUStop(3,4))

        Ce(1594) = 4*MTy2*Cd(5519)*Cd(11656) + 
     &   8*At*MTy*Ce(1593)*MStop2(1)*MStop2(2) + 
     &   MUE**2*Cd(5519)*UUStop(1,1)*UUStop(2,1) + 
     &   (-At**2 + MUEC**2)*Cd(5519)*UUStop(1,2)*UUStop(2,2) + 
     &   2*(4*MTy2*Cd(523) - MUE2*Cd(5945))*UUStop(3,1)*UUStop(3,2)

        Ce(1800) = 2*AtC*Cd(11376) + At*MTy*Ce(1798) - 
     &   8*MTy2*Ce(1799) + 
     &   (-AtC**2 + MUE**2)*UUStop(1,1)*UUStop(2,1) + 
     &   (-Ce(1797) + MUEC**2*UUStop(1,2))*UUStop(2,2)

        Ce(541) = 8*Cd(10086)*Cd(12370) + 
     &   (4*M_3C*MGlmT2*Cd(47) + 
     &      AtC*(4*MGlmT2*Cd(47) + Cd(12371)))*UCStop(1,3)**2 + 
     &   (4*M_3*MGlmT2*Cd(47) + At*Cd(12391))*UStop2(1,1)*
     &    UStop2(1,2)

        Ce(546) = 8*Cd(10039)*Cd(12370) + 
     &   (4*M_3*MGlmT2*Cd(47) + 
     &      At*(4*MGlmT2*Cd(47) + Cd(12371)))*UCStopC(1,3)**2 + 
     &   (4*M_3C*MGlmT2*Cd(47) + AtC*Cd(12391))*UStop2(1,1)*
     &    UStop2(1,2)

        Ce(598) = At*Cd(10279)*Cd(10388) - 8*Cd(10280)*Ce(597) - 
     &   8*Cd(10278)*MStop2(1)*
     &    (M_3C*MUE*UCStop(1,3)**2 + M_3*MUEC*UCStopC(1,3)**2)+
     &     AtC*Cd(10279)*(MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*UStop2(1,2)**2)

        Ce(657) = 2*MTy*Ce(656)*
     &    (MUE*UCStop(1,3) + MUEC*UCStopC(1,3)) + 
     &   Cd(251)*Cd(2681)*(AtC*Cd(10263) + At*Cd(10266) - 
     &      4*MTy2*(M_3C*MUE*UCStop(1,3)**2 + 
     &         M_3*MUEC*UCStopC(1,3)**2 + 
     &         Cd(10224)*UStop2(1,1)*UStop2(1,2)))

        Ce(984) = -(At*AtC*Cd(9966)) + MUE2*Cd(9948)*Cr(28) + 
     &   AtC**2*UCStop(1,3)**2 + At**2*UCStopC(1,3)**2 + 
     &   MUE2*UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(992) = 4*Cd(152)*Cr(28) + Opt(22) - 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2)

        Ce(1012) = 4*MStop2(3) + 
     &   3*(-At**2 + MUEC**2)*UCStopC(1,3)**2 - 
     &   2*UStop2(1,1)*(2*MTy2*Cd(9948) + MUE2*UStop2(1,1)) + 
     &   2*(2*MTy2*UStop2(1,1) + MUE2*(-Cd(9948) + UStop2(1,1)))*
     &    UStop2(1,2)

        Ce(1034) = 2*MUE**2*Cd(10280)*UCStop(1,3)**2 + 
     &   (-(At**2*Cd(10279)) + 2*MUEC**2*Cd(10280))*
     &    UCStopC(1,3)**2 - 
     &   AtC**2*Cd(10279)*UCStop(1,3)**2*UStop2(1,1) + 
     &   4*MUE2*Cd(9948)*Cd(10280)*UStop2(1,1)**2 - 
     &   AtC**2*Cd(10279)*UCStop(1,3)**2*UStop2(1,2) - 
     &   4*At*AtC*Cd(10280)*UStop2(1,1)*UStop2(1,2) + 
     &   4*MUE2*Cd(10280)*UStop2(1,1)*UStop2(1,2)**2

        Ce(1073) = 2*MTy2 + (-4*M_3C + AtC)*MTy*UCStop(1,3) + 
     &   (-4*M_3 + At)*MTy*UCStopC(1,3) + 
     &   2*MGl2*(2*Cd(9948)*UStop2(1,1) + UStop2(1,1)**2 + 
     &      UStop2(1,2)**2)

        Ce(1307) = 2*MUE**2*Cd(10604)*UCStop(1,3)**2 + 
     &   (-(At**2*Cd(10603)) + 2*MUEC**2*Cd(10604))*
     &    UCStopC(1,3)**2 - 
     &   AtC**2*Cd(10603)*UCStop(1,3)**2*UStop2(1,1) + 
     &   4*MUE2*Cd(9948)*Cd(10604)*UStop2(1,1)**2 - 
     &   AtC**2*Cd(10603)*UCStop(1,3)**2*UStop2(1,2) - 
     &   4*At*AtC*Cd(10604)*UStop2(1,1)*UStop2(1,2) + 
     &   4*MUE2*Cd(10604)*UStop2(1,1)*UStop2(1,2)**2

        Ce(1326) = -Ce(913) + 8*M_3C*AtC*UCStop(1,3)**2 + 
     &   8*MStop2(3)*UStop2(1,1) - 
     &   (-4*M_3*AtC + At*Cd(10398) + 12*MStop2(1))*
     &    UStop2(1,1)**2

        Ce(1337) = 2*MTy2 - (-4*M_3C + AtC)*MTy*UCStop(1,3) - 
     &   (-4*M_3 + At)*MTy*UCStopC(1,3) + 
     &   2*MGl2*(2*Cd(9948)*UStop2(1,1) + UStop2(1,1)**2 + 
     &      UStop2(1,2)**2)

        Ce(745) = 2*At*Cd(10095) - 2*Cd(1952)*Cd(9948)*Cr(75) + 
     &   2*MTy2*(UStop2(1,1) - UStop2(1,2))*UStop2(1,2) + 
     &   MStop2(3)*(1 - 2*UStop2(1,1)*(1 + 2*UStop2(1,2)))

        Ce(773) = MGl2*UStop2(1,1)**2 + MTy2*UStop2(1,1)**2 + 
     &   4*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))*
     &    (UStop2(1,1) - UStop2(1,2)) + 
     &   MStop2(1)*(3*UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(954) = 2*Cd(9953) + Cd(9948)*Ce(953) + 
     &   Cd(10101)*UStop2(1,1) - 
     &   Cd(9951)*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,1)**2 + 3*UStop2(1,2)**2))

        Ce(1065) = -Cd(10279) + 
     &   (MGlmT2*Cd(47)*Ce(999))/MStop2(3) - 
     &   4*Cd(10278)*(M_3C*AtC*UCStop(1,3)**2 + 
     &      M_3*At*UCStopC(1,3)**2) + 
     &   4*(-(MGlmT2*Cd(47)*Ce(996)) + Cd(9948)*Ce(1064))*
     &    UStop2(1,1)

        Ce(1190) = 4*Cd(179)*Cr(28) + 4*MStop2(2) + 
     &   AtC**2*UCStop(1,3)**2 + 
     &   2*Cd(9948)*Cd(11847)*UStop2(1,1) - 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2) + 
     &   At*(At*UCStopC(1,3)**2 + 2*AtC*UStop2(1,1)*UStop2(1,2)**2)

        Ce(1198) = 4*Cd(152)*Cr(28) + Opt(22) + 
     &   2*UStop2(1,1)*(Cd(9948)*Cd(11731) - MUE2*UStop2(1,2))

        Ce(1285) = 4*Cd(179) + Opt(22) - 
     &   2*UStop2(1,1)*(-(At*AtC*UStop2(1,1)*UStop2(1,2)) + 
     &      MUE2*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Ce(1324) = 8*Cd(50)**2 + 2*At*Cd(10645) + 
     &   AtC*(AtC*Cd(10603) - 16*M_3C*Cd(10604))*
     &    UCStop(1,3)**2 + At**2*Cd(10603)*UCStopC(1,3)**2 - 
     &   2*UStop2(1,1)*(MUE2*Cd(10603)*UStop2(1,2) + 
     &      8*M_3*AtC*Cd(10604)*
     &       (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Ce(1390) = -Ce(913) - At*Ce(1389) + 5*MStop2(3) - 
     &   MUEC**2*UCStopC(1,3)**2 - 12*MStop2(1)*UStop2(1,1)**2 + 
     &   AtC*(-12*M_3*UStop2(1,1)**2 - 4*M_3*UStop2(1,2)**2)

        Ce(519) = MGl2*Ce(517) + Cd(9948)*Ce(518) + 
     &   12*M_3C*MTy*MUE*UCStop(1,3)**3 - 
     &   2*MTy2*MUE*UCStop(1,3)**2*(1 + UStop2(1,2)) - 
     &   4*M_3C*MTy*MUEC*UCStop(1,3)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2) + 
     &   MTy2*MUEC*UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(666) = MGlmT2*Ce(665) - At*Cd(50)*UCStopC(1,3)**2 + 
     &   MStop2(2)*(2*Cd(138)*Cd(10176) + 
     &      UCStopC(1,3)**2*
     &       (M_3*MStop2(3) + 
     &         At*MStop2(1)*(-1 + Cd(9948) - UStop2(1,2))) + 
     &      M_3C*MStop2(3)*UStop2(1,1)*UStop2(1,2))

        Ce(722) = (2*Cd(50)*Ce(721))/MStop2(3) + 
     &   (4*M_3*MGlmT2*Cd(50) + 8*M_3*Cd(12482) + 
     &      At*Cd(12483))*UCStopC(1,3)**2 + 
     &   4*M_3C*(MGlmT2*Cd(50) + 2*Cd(12482))*UStop2(1,1)*
     &    UStop2(1,2) + AtC*
     &    (2*MGlmT2*Cd(50) + Cd(12483)*UStop2(1,1)*UStop2(1,2))

        Ce(1004) = -(Cd(9948)*Ce(1003)) + Opt(26) + 
     &   UCStopC(1,3)*UStop2(1,2)*
     &    (At**2*UStop2(1,1) + MUEC**2*UStop2(1,2)) + 
     &   2*UCStop(1,3)*(MStop2(3) + 2*MTy2*UStop2(1,2)**2 + 
     &      MUE2*UStop2(1,2)**2)

        Ce(1011) = Cd(9948)*Ce(1010) + Opt(26) - 
     &   UCStopC(1,3)*UStop2(1,2)*
     &    (MUEC**2*UStop2(1,1) + At**2*UStop2(1,2)) - 
     &   2*UCStop(1,3)*(2*MStop2(3) + 
     &      2*MTy2*UStop2(1,1)*UStop2(1,2) + 
     &      MUE2*UStop2(1,1)*UStop2(1,2))

        Ce(1127) = Cd(12370)*Cd(12515) + 
     &   (2*Cd(47)*Ce(1126))/MStop2(3) + 
     &   (4*M_3C*MGlmT2*Cd(47) - AtC*Ce(1123))*UCStop(1,3) + 
     &   UCStopC(1,3)*(4*M_3*MGlmT2*Cd(47) + 
     &      At*(Cd(12371) + 
     &         2*MGlmT2*Cd(47)*(-1 + 2*UStop2(1,1) + 2*UStop2(1,2))
     &         ))

        Ce(1276) = -2*At*AtC*Cd(10020)*Cr(5) - 
     &   AtC**2*Cd(10021)*UCStop(1,3) + 
     &   2*(-(MUE2*Cd(9948)*Cr(5)*Cr(75)) + 
     &      Cd(585)*Cd(9955)*Cr(75)*MStop2(1) + 
     &      MUE2*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))*
     &       UStop2(1,1)**2) + 
     &   At**2*UCStopC(1,3)*
     &    (3*M_3*UCStopC(1,3)**2 - 
     &      M_3C*(UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &         UStop2(1,2)**2))

        Ce(1374) = Opt(25) + 
     &   UStop2(1,1)*(MStop2(3) - 
     &      MUE2*(UStop2(1,1) - UStop2(1,2))**2) + 
     &   MUEC**2*UCStopC(1,3)**2*(1 + UStop2(1,2))

        Ce(1418) = Opt(27) - 2*MUE2*Cd(9948)*UStop2(1,1)**2 + 
     &   4*MStop2(2)*(Cd(9948) + UStop2(1,1)) + 
     &   UStop2(1,2)*(AtC**2*UCStop(1,3)**2 + 
     &      2*UStop2(1,1)*(At*AtC - MUE2*UStop2(1,2)))

        Ce(1420) = Opt(25) + 
     &   UStop2(1,1)*(MStop2(3) - 
     &      MUE2*(UStop2(1,1) - UStop2(1,2))**2) + 
     &   MUEC**2*UCStopC(1,3)**2*(1 + UStop2(1,2))

        Ce(1809) = AtC**2*UUStop(1,1)*UUStop(2,1) - 
     &   MUE**2*UUStop(1,1)*UUStop(2,1) + 
     &   At**2*UUStop(1,2)*UUStop(2,2) - 
     &   MUEC**2*UUStop(1,2)*UUStop(2,2) - 
     &   4*MTy2*(UUStop(1,2)*UUStop(2,1) + 
     &      UUStop(1,1)*UUStop(2,2)) - 
     &   4*At*MTy*(UUStop(1,3)*UUStop(2,2) + 
     &      UUStop(1,2)*UUStop(2,3)) - 
     &   2*AtC*Cd(11132)*UUStop(3,1) - 
     &   2*Ce(1808)*UUStop(3,1)*UUStop(3,2)

        Ce(662) = MGlmT2*Ce(661) - AtC*Cd(50)*UCStop(1,3)**2 + 
     &   MStop2(2)*(2*Cd(138)*Cd(10169) - 
     &      AtC*MStop2(1)*UCStop(1,3)**2 + 
     &      MStop2(3)*(M_3*UStop2(1,1)*UStop2(1,2) + 
     &         UCStop(1,3)**2*
     &          (M_3C + AtC*(-Cd(9948) + UStop2(1,2)))))

        Ce(1035) = -(At*AtC*Cd(9966)) + Cd(9948)*Cd(10450) + 
     &   AtC**2*UCStop(1,3)**2 + At**2*UCStopC(1,3)**2 - 
     &   (MUE**2*UCStop(1,3)**2 + MUEC**2*UCStopC(1,3)**2)*
     &    (1 + UStop2(1,2)) + 
     &   MUE2*UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(1063) = 8*Cd(47)**2*Cd(10208) + 2*At*Cd(10351) + 
     &   AtC*(AtC*Cd(10279) - 16*M_3C*Cd(10280))*
     &    UCStop(1,3)**2 + At**2*Cd(10279)*UCStopC(1,3)**2 - 
     &   2*UStop2(1,1)*(MUE2*Cd(10279)*UStop2(1,2) + 
     &      8*M_3*AtC*Cd(10280)*
     &       (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Ce(1072) = -8*Cd(47)**2*Cd(492) - 
     &   AtC*(AtC*Cd(164)*Cd(2597) + 8*M_3C*Cd(12298))*
     &    UCStop(1,3)**2 - At*Ce(1071)*UCStopC(1,3)**2 + 
     &   2*MUE2*Cd(164)*Cd(2597)*Cd(9948)*UStop2(1,1)**2 - 
     &   2*At*(AtC*Cd(164)*Cd(2597) + 4*M_3C*Cd(12298))*
     &    UStop2(1,1)*UStop2(1,2) + 
     &   2*UStop2(1,1)*UStop2(1,2)*
     &    (-4*M_3*AtC*Cd(12298) + 
     &      MUE2*Cd(164)*Cd(2597)*UStop2(1,2))

        Ce(1088) = MUE2*Cd(9948)*Cr(28) + 8*Cd(492)*MStop2(3) + 
     &   AtC**2*UCStop(1,3)**2 + At**2*UCStopC(1,3)**2 - 
     &   At*AtC*(UStop2(1,1)**2 + UStop2(1,2)**2) + 
     &   MUE2*UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(1132) = -(Cd(118)*Cd(11723)) + Cd(9948)*Ce(1131) - 
     &   AtC*(2*M_3C + AtC)*UCStop(1,3)**2 - 
     &   At*Cd(10001)*UCStopC(1,3)**2 - 
     &   2*At*(M_3C + AtC)*UStop2(1,1)*UStop2(1,2) + 
     &   2*UStop2(1,1)*UStop2(1,2)*
     &    (-(M_3*AtC) + MUE2*UStop2(1,2))

        Ce(1133) = -4*Cd(492)*MStop2(3) - 
     &   AtC**2*UCStop(1,3)**2 - At**2*UCStopC(1,3)**2 + 
     &   At*AtC*(UStop2(1,1)**2 + UStop2(1,2)**2) - 
     &   MUE2*(Cd(9948)*Cr(28) + 
     &      UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2))

        Ce(1168) = 4*AtC*Cd(10169) - At*Cd(11749) + 
     &   AtC**2*UCStop(1,3)**2 + 
     &   UCStopC(1,3)**2*(At**2 + MUEC**2*(-1 - UStop2(1,2))) + 
     &   12*MStop2(3)*UStop2(1,1)**2*UStop2(1,2) - 
     &   MUE**2*UCStop(1,3)**2*(1 + UStop2(1,2)) + 
     &   MUE2*UStop2(1,1)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(1248) = Ce(1247) + AtC**2*UCStop(1,3)**2 + 
     &   At*(At*UCStopC(1,3)**2 + 2*AtC*UStop2(1,1)*UStop2(1,2)) + 
     &   MUE2*(-2*Cd(9948)*UStop2(1,1)**2 - 
     &      2*UStop2(1,1)*UStop2(1,2)**2)

        Ce(1336) = -8*Cd(50)**2*Cd(492) - 
     &   AtC*(AtC*Cd(251)*Cd(2681) + 8*M_3C*Cd(4187))*
     &    UCStop(1,3)**2 - At*Ce(1335)*UCStopC(1,3)**2 + 
     &   2*MUE2*Cd(251)*Cd(2681)*Cd(9948)*UStop2(1,1)**2 - 
     &   2*At*(AtC*Cd(251)*Cd(2681) + 4*M_3C*Cd(4187))*
     &    UStop2(1,1)*UStop2(1,2) + 
     &   2*UStop2(1,1)*UStop2(1,2)*
     &    (-4*M_3*AtC*Cd(4187) + 
     &      MUE2*Cd(251)*Cd(2681)*UStop2(1,2))

        Ce(1434) = 4*AtC*Cd(10172) - At*Cd(11889) + 
     &   AtC**2*UCStop(1,3)**2 + At**2*Cd(9972)*UCStopC(1,3)**2 + 
     &   UStop2(1,1)**2*(24*Cd(9948)*MStop2(3) + 
     &      MUE2*(-Cd(9948) + UStop2(1,1))) + 
     &   2*(MUE2 - 12*MStop2(3))*UStop2(1,1)**2*UStop2(1,2) + 
     &   MUE2*UStop2(1,1)*UStop2(1,2)**2 + MUE2*UStop2(1,2)**3

        Ce(508) = -(MUE*MStop2(1)*UCStop(1,3)**2) - 
     &   MUEC*MStop2(1)*UStop2(1,1)*UStop2(1,2) + 
     &   MTy2*(MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*UStop2(1,2)) + 
     &   MTy*(MUE*UCStop(1,3)*
     &       (-2*M_3C*UCStop(1,3)**2 - 
     &         2*M_3*UStop2(1,1)*UStop2(1,2)) + 
     &      MUEC*(-2*M_3C*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) - 
     &         2*M_3*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2))) + 
     &   MGl2*(MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*
     &       (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Ce(512) = Cd(9948)*Ce(511) + 
     &   12*M_3C*MTy*MUE*UCStop(1,3)**3 - 
     &   2*MUE*UCStop(1,3)**2*(MGlpT2 + MGl2*UStop2(1,2)) - 
     &   4*M_3C*MTy*MUEC*UCStop(1,3)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2) + 
     &   MUEC*(MTy2 + MGl2*UStop2(1,1))*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(577) = MUEC*MGlpTmSt2(1)*UCStopC(1,3)**2 - 
     &   2*M_3*MTy*MUEC*UCStopC(1,3)**3 - 
     &   2*MTy*(M_3*MUE + M_3C*MUEC)*UCStopC(1,3)*UStop2(1,1)*
     &    UStop2(1,2) + MUE*UStop2(1,1)*UStop2(1,2)*
     &    (MGl2 - 2*M_3C*MTy*UCStop(1,3) + 
     &      (MTy2 - MStop2(1))*UStop2(1,2))

        Ce(698) = -2*At*Cd(251)*Cd(10095) - 
     &   2*Cd(9948)*Ce(697)*Cr(75) + 
     &   2*MStop2(2)*(MStop2(3) + 
     &      MStop2(1)*(UStop2(1,1) - UStop2(1,2))*UStop2(1,2)) + 
     &   MGlmT2*MStop2(3)*(-UStop2(1,2) + 
     &      UStop2(1,1)*(-1 + 2*UStop2(1,1) + 6*UStop2(1,2)))

        Ce(711) = (2*Cd(50)*Ce(710))/MStop2(3) + 
     &   4*M_3*(MGlmT2*Cd(50) + 2*Cd(12482))*UStop2(1,1)*
     &    UStop2(1,2) + At*
     &    (2*MGlmT2*Cd(50) + Cd(12483)*UStop2(1,1)*UStop2(1,2)) + 
     &   UCStop(1,3)**2*(4*M_3C*MGlmT2*Cd(50) + 
     &      8*M_3C*Cd(12482) + 
     &      AtC*(Cd(12483) + 
     &         4*MGlmT2*Cd(50)*(-Cd(9948) + UStop2(1,2))))

        Ce(755) = MUEC*MGlpTmSt2(2)*UCStopC(1,3)**2 + 
     &   2*M_3*MTy*MUEC*UCStopC(1,3)**3 + 
     &   2*MTy*(M_3*MUE + M_3C*MUEC)*UCStopC(1,3)*UStop2(1,1)*
     &    UStop2(1,2) + MUE*UStop2(1,1)*UStop2(1,2)*
     &    (MGl2 + 2*M_3C*MTy*UCStop(1,3) + 
     &      (MTy2 - MStop2(2))*UStop2(1,2))

        Ce(968) = 4*MTy2 + Opt(27) + 
     &   4*MTy*(AtC*UCStop(1,3) + At*UCStopC(1,3)) - 
     &   2*MUE2*Cd(9948)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(AtC**2*UCStop(1,3)**2 + 
     &      2*UStop2(1,1)*(At*AtC - MUE2*UStop2(1,2)))

        Ce(987) = -(At*AtC*Cd(9966)) + Cd(9948)*Cd(10450) + 
     &   Opt(27) + (AtC**2*UCStop(1,3)**2 - 
     &      MUE**2*UCStop(1,3)**2 - MUEC**2*UCStopC(1,3)**2)*
     &    UStop2(1,2) + MUE2*UStop2(1,1)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(1105) = -(AtC*(-8*M_3C + AtC)*UCStop(1,3)**2) + 
     &   8*M_3*At*UCStopC(1,3)**2 - At**2*UCStopC(1,3)**2 + 
     &   2*MUE2*Cd(9948)*UStop2(1,1)**2 - 
     &   2*At*(-4*M_3C + AtC)*UStop2(1,1)*UStop2(1,2) + 
     &   2*UStop2(1,1)*UStop2(1,2)*(4*M_3*AtC + MUE2*UStop2(1,2))

        Ce(1258) = At*AtC*Cd(10070) + MUEC**2*UCStopC(1,3)**2 + 
     &   At**2*Cd(10077)*UCStopC(1,3)**2 + 
     &   UCStop(1,3)**2*(-AtC**2 + 
     &      MUE**2*(Cd(9948) + UStop2(1,1))) + 
     &   MUE2*((Cd(9948) - UStop2(1,1))*UStop2(1,1)**2 - 
     &      2*UStop2(1,1)**2*UStop2(1,2) - 
     &      UStop2(1,1)*UStop2(1,2)**2 - UStop2(1,2)**3)

        Ce(1428) = Cd(179)*Cd(10141) + AtC*Cd(11871) + Opt(24) - 
     &   M_3*MUEC**2*UCStopC(1,3)**3 - 
     &   MUE**2*UCStop(1,3)*
     &    (M_3C*UCStop(1,3)**2 + M_3*UStop2(1,1)*UStop2(1,2))+
     &     UCStopC(1,3)*(2*At*MStop2(2) - 
     &      M_3C*MUEC**2*UStop2(1,1)*UStop2(1,2))

        Ce(738) = -(MGlpT2*Cd(13)) + 
     &   M_3*At*MUEC*UCStopC(1,3)**3 + 
     &   UCStopC(1,3)*(MUEC*MStop2(2) + 
     &      M_3*AtC*MUEC*UStop2(1,1)*UStop2(1,2) + 
     &      At*(M_3*MUE*UStop2(1,1)*UStop2(1,2) + 
     &         M_3C*MUEC*UStop2(1,1)*UStop2(1,2))) + 
     &   UCStop(1,3)*(M_3C*AtC*MUEC*UStop2(1,1)*UStop2(1,2) + 
     &      MUE*(MStop2(2) + M_3C*At*UStop2(1,1)*UStop2(1,2) + 
     &         AtC*(M_3C*UCStop(1,3)**2 + 
     &            M_3*UStop2(1,1)*UStop2(1,2))))

        Ce(740) = -(AtC*Cd(10087)*UCStop(1,3)) + 
     &   MStop2(3)*UCStop(1,3) + 
     &   MStop2(1)*UCStop(1,3)*(Cd(9948) - UStop2(1,2))*
     &    (-UStop2(1,1) + UStop2(1,2)) + 
     &   At*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))*
     &    (UStop2(1,1)**2 + Cd(9948)*(-UStop2(1,1) + UStop2(1,2)))

	Ce(815) = Ce(755) + MUE*Cd(48)*Cd(9948)*UStop2(1,1)**2

        Ce(817) = -(MUE*MStop2(2)*UCStop(1,3)**2) + 
     &   MGl2*(MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*UStop2(1,2)) + 
     &   MUEC*UStop2(1,1)*(Cd(48)*Cd(9948)*UStop2(1,1) - 
     &      MStop2(2)*UStop2(1,2)**2) + 
     &   MTy2*(MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*UStop2(1,2)**2) + 
     &   2*MTy*(MUE*UCStop(1,3)*
     &       (M_3C*UCStop(1,3)**2 + 
     &         M_3*UStop2(1,1)*UStop2(1,2)) + 
     &      MUEC*(M_3C*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) + 
     &         M_3*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)))

        Ce(986) = At*AtC*Cd(9969)*Cr(28) + 2*MStop2(3) + 
     &   MUEC**2*UCStopC(1,3)**2 + 
     &   At**2*Cd(10077)*UCStopC(1,3)**2 + 
     &   UCStop(1,3)**2*(MUE**2 + 
     &      AtC**2*(-1 + Cd(9948) - UStop2(1,2))) - 
     &   MUE2*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(1108) = -(At*Cd(11749)) + MUE2*Cd(9948)*Cr(28) + 
     &   AtC*(-8*M_3C + AtC)*UCStop(1,3)**2 + 
     &   At**2*UCStopC(1,3)**2 + 
     &   (4*M_3*AtC + MUE2*UStop2(1,1))*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(1138) = 2*Cd(9948)*Ce(1137) + Opt(22) - 
     &   MUE**2*UCStop(1,3)**2 - MUEC**2*UCStopC(1,3)**2 + 
     &   4*UStop2(1,1)*(MStop2(1) + Cd(152)*UStop2(1,1)) - 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2) + 4*Cd(152)*UStop2(1,2)**2

        Ce(1206) = At*AtC*Cd(11983) + 
     &   AtC**2*(-Cd(47) + Cd(50))*UCStop(1,3)**2 + 
     &   At**2*Cd(11985)*UCStopC(1,3)**2 + 
     &   MUE2*(Cd(9948)*Cd(11977) + 
     &      (-Cd(47) + Cd(50))*UStop2(1,1)**3 - 
     &      2*Cd(47)*UStop2(1,1)**2*UStop2(1,2) + 
     &      (-Cd(47) + Cd(50))*UStop2(1,1)*UStop2(1,2)**2 - 
     &      Cd(47)*UStop2(1,2)**3)

        Ce(1250) = At*Ce(1249) + Ce(1838) + 4*MStop2(2) + 
     &   4*(-2*M_3C + AtC)*MTy*UCStop(1,3) - 
     &   AtC*(-8*M_3C + AtC)*UCStop(1,3)**2 + 
     &   MUE**2*UCStop(1,3)**2 - 8*M_3*MTy*UCStopC(1,3) + 
     &   MUEC**2*UCStopC(1,3)**2 + 
     &   8*M_3*AtC*UStop2(1,1)*
     &    (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

        Ce(1328) = At*(-Cd(10660) + At*UCStopC(1,3)**2) + 
     &   4*Ce(996)*MStop2(3)*UStop2(1,1) - 
     &   8*M_3C*AtC*UCStop(1,3)**2*UStop2(1,1) + 
     &   AtC**2*UCStop(1,3)**2*UStop2(1,1) + 
     &   UStop2(1,1)*(MUE2*(UStop2(1,1) - UStop2(1,2))**2 - 
     &      12*MStop2(1)*UStop2(1,1)*UStop2(1,2)) + 
     &   4*M_3*AtC*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))))

        Ce(1401) = 2*Cd(9948)*Ce(1400) + Opt(22) - 
     &   MUE**2*UCStop(1,3)**2 - MUEC**2*UCStopC(1,3)**2 + 
     &   4*UStop2(1,1)*(MStop2(2) + Cd(179)*UStop2(1,1)) - 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2) + 4*Cd(179)*UStop2(1,2)**2

        Ce(970) = -2*At*Ce(969) + Ce(1838) + 4*MStop2(1) - 
     &   4*(-2*M_3C + AtC)*MTy*UCStop(1,3) - 
     &   AtC*(-8*M_3C + AtC)*UCStop(1,3)**2 + 
     &   MUE**2*UCStop(1,3)**2 + 8*M_3*MTy*UCStopC(1,3) + 
     &   MUEC**2*UCStopC(1,3)**2 + 
     &   8*M_3*AtC*UStop2(1,1)*
     &    (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)

        Ce(993) = 4*Cd(152) + Opt(22) - MUE**2*UCStop(1,3)**2 - 
     &   MUEC**2*UCStopC(1,3)**2 + 
     &   4*MStop2(1)*(Cd(9948) + UStop2(1,1)) - 
     &   2*UStop2(1,1)*(MUE2*Cd(9948)*UStop2(1,1) + 
     &      UStop2(1,2)*(-(At*AtC*UStop2(1,1)) + MUE2*UStop2(1,2)))

        Ce(1164) = AtC*Ce(1163) - Ce(1162)*MStop2(3) + 
     &   (-MUE**2 + AtC**2*Cr(85))*UCStop(1,3)**2 + 
     &   (8*M_3C*At + MUE2)*(UStop2(1,1)**2 + UStop2(1,2)**2) + 
     &   UCStopC(1,3)**2*(-MUEC**2 + 
     &      At*(-16*M_3 + At*(1 + UStop2(1,2))))

        Ce(1259) = Opt(27) + Opt(28) - 
     &   MUEC**2*UCStopC(1,3)**2*UStop2(1,2)

        Ce(1308) = Opt(28) + 
     &   UCStop(1,3)**2*(-MUE**2 + AtC**2*UStop2(1,1)) + 
     &   UCStopC(1,3)**2*(At**2 + MUEC**2*(-1 - UStop2(1,2)))

        Ce(1388) = -(At*Ce(1387)) - 
     &   8*M_3C*AtC*UCStop(1,3)**2 + AtC**2*UCStop(1,3)**2 + 
     &   (-MUEC**2 + At**2*Cr(85))*UCStopC(1,3)**2 + 
     &   4*M_3*AtC*UStop2(1,1)*(UStop2(1,1) - UStop2(1,2))**2 - 
     &   UStop2(1,1)*(MUE**2*UCStop(1,3)**2 + 
     &      24*MStop2(3)*UStop2(1,1)*UStop2(1,2)) + 
     &   MUE2*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))))

        Ce(452) = ((-4*M_3C + AtC)*
     &       (Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952)) + 
     &      4*M_3C*Cd(50)**2*Cd(9953) + 
     &      4*M_3C*Cd(47)**2*Cd(9954))*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(4*M_3*Cd(9948)*Ce(451) + 
     &      UStop2(1,2)*(At*Cd(9980) - 
     &         4*M_3*(Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952))*
     &          UStop2(1,2)))

        Ce(464) = ((-4*M_3 + At)*
     &       (Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952)) + 
     &      4*M_3*Cd(50)**2*Cd(9953) + 4*M_3*Cd(47)**2*Cd(9954)
     &      )*UCStopC(1,3)**2 + 
     &   UStop2(1,1)*(4*M_3C*Cd(9948)*Ce(451) + 
     &      UStop2(1,2)*(AtC*Cd(9980) - 
     &         4*M_3C*(Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952))*
     &          UStop2(1,2)))

        Ce(566) = M_3*At*UCStopC(1,3)*
     &    ((Cd(9948) - UStop2(1,1))*UStop2(1,1) - UStop2(1,2)**2)+
     &     UCStop(1,3)*(AtC*Cd(10088) + MStop2(3) + 
     &      (Cd(9948) - UStop2(1,2))*
     &       (MStop2(1)*UStop2(1,1) + MGlpT2*UStop2(1,2)) + 
     &      M_3C*At*((Cd(9948) - UStop2(1,1))*UStop2(1,1) - 
     &         UStop2(1,2)**2))

        Ce(1367) = Cd(12482)*Cd(12515) + 
     &   (4*M_3C*MGlmT2*Cd(50) - AtC*Cd(12514))*UCStop(1,3) + 
     &   UCStopC(1,3)*(4*M_3*MGlmT2*Cd(50) + 
     &      At*(Cd(12483) + 
     &         2*MGlmT2*Cd(50)*
     &          (UStop2(1,1)*(-1 - 2*Cd(9948) + 2*UStop2(1,1)) + 
     &            (-1 + 2*Cd(9948) + 6*UStop2(1,1))*UStop2(1,2))))

        Ce(1383) = 4*MGl2*Ce(1382) - 
     &   8*M_3C*AtC*UCStop(1,3)**2 + AtC**2*UCStop(1,3)**2 - 
     &   MUE**2*UCStop(1,3)**2 + 
     &   (At*(-8*M_3 + At) - MUEC**2)*UCStopC(1,3)**2 + 
     &   2*Cd(9948)*Ce(1338)*UStop2(1,1)**2 - 
     &   8*M_3C*At*UStop2(1,1)*UStop2(1,2) - 
     &   8*M_3*AtC*UStop2(1,1)*UStop2(1,2) + 
     &   2*At*AtC*UStop2(1,1)*UStop2(1,2) - 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2)**2 - 
     &   4*MStop2(2)*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,1)**2 + 3*UStop2(1,2)**2))

        Ce(638) = At*Ce(637) + 
     &   4*Cd(152)*Cd(164)*Cd(10224)*MStop2(3)*UStop2(1,1) + 
     &   2*AtC*MUEC*Cd(164)*MStop2(1)*MStop2(3)*UStop2(1,1)**2 - 
     &   Cd(47)*(AtC*MStop2(1)*
     &       (2*MUE*UCStop(1,3)**2 + 
     &         MUEC*(UStop2(1,1)**2 + UStop2(1,2)**2)) + 
     &      MGlmT2*(8*M_3*MUEC*UCStopC(1,3)**2 + 
     &         4*M_3*MUE*(UStop2(1,1)**2 + UStop2(1,2)**2) - 
     &         (-4*M_3C + AtC)*
     &          (2*MUE*UCStop(1,3)**2 + 
     &            MUEC*(UStop2(1,1)**2 + UStop2(1,2)**2))))

        Ce(774) = Cd(9948)*Ce(773) + MStop2(3) - 
     &   4*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))*
     &    UStop2(1,1)**2 + 
     &   MStop2(1)*UStop2(1,1)*(UStop2(1,1) - UStop2(1,2))**2 - 
     &   MGl2*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2)))) - 
     &   MTy2*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))))

        Ce(1008) = M_3C*Ce(1004) - Cd(9948)*Ce(1007) + 
     &   3*M_3*(At - MUEC)*(At + MUEC)*UCStopC(1,3)**3 + 
     &   M_3*MUE**2*UCStop(1,3)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2) + 
     &   UCStopC(1,3)*(2*M_3*MStop2(3) + At*MStop2(3) + 
     &      2*M_3*MUE2*UStop2(1,1)**2 + 
     &      2*At*(MGl2 - MStop2(1))*UStop2(1,1)*UStop2(1,2) + 
     &      2*(2*M_3*MTy2 - At*MTy2 + M_3*MUE2)*UStop2(1,2)**2)

        Ce(600) = At*Ce(599) - 
     &   8*Cd(10224)*Cd(10278)*MStop2(1)*MStop2(3)*UStop2(1,1) + 
     &   (AtC*MUEC*Cd(10279) - 8*M_3*MUE*Cd(10280) - 
     &      8*M_3C*MUEC*Cd(10280))*MStop2(3)*UStop2(1,1)**2 + 
     &   2*Cd(47)*MStop2(1)*
     &    (AtC*MStop2(1)*(2*MUE*UCStop(1,3)**2 + 
     &         MUEC*(UStop2(1,1)**2 + UStop2(1,2)**2)) + 
     &      MGlmT2*(8*M_3*MUEC*UCStopC(1,3)**2 + 
     &         4*M_3*MUE*(UStop2(1,1)**2 + UStop2(1,2)**2) - 
     &         (-4*M_3C + AtC)*
     &          (2*MUE*UCStop(1,3)**2 + 
     &            MUEC*(UStop2(1,1)**2 + UStop2(1,2)**2))))

        Ce(918) = 2*Cd(9948)*Ce(917) + 
     &   2*Cd(47)*Cd(50)**2*Cd(164)*Cd(9947)*UStop2(1,1) + 
     &   2*Cd(47)**2*Cd(50)*Cd(251)*Cd(9949)*UStop2(1,1) - 
     &   (Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952))*
     &    (AtC**2*UCStop(1,3)**2 - 
     &      2*MUE2*UStop2(1,1)*UStop2(1,2)**2 + 
     &      At*(At*UCStopC(1,3)**2 + 2*AtC*UStop2(1,1)*UStop2(1,2))
     &      )

        Ce(1087) = Cd(9948)*Ce(1086) - 
     &   2*M_3C*AtC*Cd(47)*UCStop(1,3)**2 - 
     &   At*Cd(47)*Cd(10001)*UCStopC(1,3)**2 + 
     &   Cd(12369)*UStop2(1,1) - 
     &   AtC**2*Cd(47)*UCStop(1,3)**2*UStop2(1,1) - 
     &   4*Cd(12370)*UStop2(1,1)**2 - 
     &   AtC**2*Cd(47)*UCStop(1,3)**2*UStop2(1,2) - 
     &   2*M_3C*At*Cd(47)*UStop2(1,1)*UStop2(1,2) - 
     &   2*M_3*AtC*Cd(47)*UStop2(1,1)*UStop2(1,2) - 
     &   8*Cd(12370)*UStop2(1,1)*UStop2(1,2) - 
     &   2*At*AtC*Cd(47)*UStop2(1,1)**2*UStop2(1,2) - 
     &   4*Cd(12370)*UStop2(1,2)**2 - 
     &   2*At*AtC*Cd(47)*UStop2(1,1)*UStop2(1,2)**2 + 
     &   2*MUE2*Cd(47)*UStop2(1,1)*UStop2(1,2)**2 - 
     &   Cd(47)*Cd(118)*(4*UStop2(1,1)**2 + 4*UStop2(1,2)**2 + 
     &      UStop2(1,1)*(-5 + 8*UStop2(1,2)))

        Ce(1251) = 4*MTy2**2*Cd(10193) + MGl2*Ce(1248) - 
     &   MTy2*Ce(1250) + MUEC**2*MStop2(2)*UCStopC(1,3)**2 + 
     &   2*M_3*MTy*(At - MUEC)*(At + MUEC)*UCStopC(1,3)**3 + 
     &   UCStop(1,3)*(4*AtC*MTy*MStop2(2) + 
     &      MUE**2*MStop2(2)*UCStop(1,3) + 
     &      2*M_3C*MTy*(AtC - MUE)*(AtC + MUE)*UCStop(1,3)**2 + 
     &      2*MTy*(2*M_3C*At*AtC + M_3*AtC**2 - M_3*MUE**2 - 
     &         2*M_3C*MUE2)*UStop2(1,1)*UStop2(1,2)) + 
     &   2*MTy*UCStopC(1,3)*
     &    (2*At*MStop2(2) + 
     &      (At*(M_3C*At + 2*M_3*AtC) - 2*M_3*MUE2 - 
     &         M_3C*MUEC**2)*UStop2(1,1)*UStop2(1,2))

        Ce(772) = Cd(9948)*Ce(771) + 
     &   MGl2*Cr(75)**2*UStop2(1,1) + 
     &   4*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))*
     &    (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &      UStop2(1,2)**2) - 
     &   MStop2(3)*UStop2(1,1)*
     &    (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &      3*UStop2(1,2)**2) + 
     &   MTy2*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2)))) - 
     &   MStop2(1)*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))))

        Ce(971) = 4*MTy2**2*Cd(10193) + MGl2*Ce(968) - 
     &   MTy2*Ce(970) + MUEC**2*MStop2(1)*UCStopC(1,3)**2 + 
     &   2*M_3*MTy*(-At**2 + MUEC**2)*UCStopC(1,3)**3 + 
     &   UCStop(1,3)*(-4*AtC*MTy*MStop2(1) + 
     &      MUE**2*MStop2(1)*UCStop(1,3) + 
     &      2*M_3C*MTy*(-AtC**2 + MUE**2)*UCStop(1,3)**2 + 
     &      2*MTy*(-2*M_3C*At*AtC - M_3*AtC**2 + 
     &         M_3*MUE**2 + 2*M_3C*MUE2)*UStop2(1,1)*
     &       UStop2(1,2)) - 
     &   2*MTy*UCStopC(1,3)*
     &    (2*At*MStop2(1) + 
     &      (At*(M_3C*At + 2*M_3*AtC) - 2*M_3*MUE2 - 
     &         M_3C*MUEC**2)*UStop2(1,1)*UStop2(1,2))

        Ce(935) = Cd(47)**2*Cd(9954)*Cd(10141) + 
     &   (-2*M_3C*(Cd(50)**2*
     &          (-Cd(9951) + 2*Cd(9953) + Cd(10101)) + 
     &         Cd(47)**2*(Cd(9952) - Cd(10102))) + AtC*Cd(10139))*
     &    UCStop(1,3) + (At*
     &       (Cd(47)**2*Cd(9952) + 
     &         Cd(50)**2*(-Cd(9951) + Cd(9953))) - 
     &      2*M_3*(Cd(50)**2*
     &          (-Cd(9951) + 2*Cd(9953) + Cd(10101)) + 
     &         Cd(47)**2*(Cd(9952) - Cd(10102))))*UCStopC(1,3)

        Ce(1161) = 2*Cd(9948)*Ce(1159) - 4*Ce(1160)*MStop2(1) + 
     &   Opt(29) + 4*MGl2*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,1)**2 + 3*UStop2(1,2)**2))

	Ce(1433) = -2*Cd(9948)*Ce(1432) + Opt(29) + Opt(30)

        Ce(1256) = Opt(30) + Opt(31) + 
     &   4*MStop2(2)*(1 - 3*UStop2(1,1)**2)*UStop2(1,2)

	Ce(983) = 4*Cd(152) + Opt(31)

        Ce(909) = Cd(47)*Cd(50)*
     &    (T134(MGl2,MTy2,MStop2(1)) + T134(MGl2,MTy2,MStop2(2)))+
     &     AtC*(Cd(50)*Cd(9947)*MGlpTmSt2(1)*
     &       (-2*M_3C*UCStop(1,3)**2 - 
     &         2*M_3*UStop2(1,1)*UStop2(1,2)) + 
     &      Cd(47)*Cd(9949)*MGlpTmSt2(2)*
     &       (-2*M_3C*UCStop(1,3)**2 - 
     &         2*M_3*UStop2(1,1)*UStop2(1,2)) + 
     &      AtC*(Cd(50)*Cd(9947)*MGlpTmSt2(1)*UCStop(1,3)**2*
     &          (-UStop2(1,1) - UStop2(1,2)) + 
     &         Cd(47)*(Cd(9948)*Cd(9949)*MGlpTmSt2(2)*
     &             UCStop(1,3)**2 - 
     &            Cd(9949)*MGlpTmSt2(2)*UCStop(1,3)**2*UStop2(1,2))
     &         )) + MUE2*(Cd(50)*Cd(9947)*
     &       (2*Cd(9948)*MGlpTmSt2(1)*UStop2(1,1)**2 + 
     &         2*MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2)**2) + 
     &      Cd(47)*(Cd(9948)*Cd(9949)*MGlpTmSt2(2)*
     &          (4*UStop2(1,1)**2 + UStop2(1,2)**2) + 
     &         Cd(9949)*MGlpTmSt2(2)*UStop2(1,1)*
     &          (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &            3*UStop2(1,2)**2))) + 
     &   At*(-(Ce(908)*UCStopC(1,3)**2) - 
     &      2*M_3C*Cd(50)*Cd(9947)*MGlpTmSt2(1)*UStop2(1,1)*
     &       UStop2(1,2) - 
     &      2*M_3C*Cd(47)*Cd(9949)*MGlpTmSt2(2)*UStop2(1,1)*
     &       UStop2(1,2) + 
     &      AtC*(Cd(50)*Cd(9947)*MGlpTmSt2(1)*UStop2(1,1)*
     &          (-2*UStop2(1,1)*UStop2(1,2) - 2*UStop2(1,2)**2) + 
     &         Cd(47)*Cd(9949)*MGlpTmSt2(2)*
     &          (-UStop2(1,2)**3 + 
     &            UStop2(1,1)*
     &             (UStop2(1,1)*(-UStop2(1,1) - 3*UStop2(1,2)) - 
     &               3*UStop2(1,2)**2))))

        Opt(23) = -(AtC*Cd(12015)) - Ce(1838) + 
     &   AtC**2*UCStop(1,3)**2 - 
     &   8*At*(M_3*UCStopC(1,3)**2 + 
     &      M_3C*UStop2(1,1)*UStop2(1,2))

	Ce(10) = Ce(9)*UStop2(1,1) - MStop2(1)*UStop2(1,2)**2

        Ce(1000) = (MGlmT2*Cd(47)*Ce(999) + 4*Ce(997)*MStop2(3))/
     &   MStop2(3)

	Ce(1013) = MUE**2*Cd(10018)*UCStop(1,3) - Ce(1012)*UCStopC(1,3)

        Ce(1018) = -(Ce(1017)*UCStopC(1,3)) + 
     &   MUE**2*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)

	Ce(1036) = MGlmT2*Ce(1035) + Ce(986)*MStop2(1)

	Ce(1037) = 4*Cd(47)*Ce(1036)*MStop2(1) - Ce(1034)*MStop2(3)

        Ce(105) = 1/2.D0*(MZ2*SW2*A0delStop(1)*Ce(91) + 
     &      2*Ce(104)*Cr(13)*MStop2(1))/(MZ2*SW2)

	Ce(1052) = M_3C*Ce(1049) + Ce(1051)*UCStopC(1,3)

	Ce(108) = (2*Ce(104) + MZ2*SW2*Ce(91)*Cr(13))/(MZ2*SW2)

	Ce(1077) = MGl2*Cd(9948)*UStop2(1,1) + Ce(1076)*UStop2(1,2)

        Ce(1079) = -(Ce(1078)*UCStopC(1,3)) + 
     &   2*M_3C*MTy2*UStop2(1,1)*UStop2(1,2)

        Ce(1080) = Cd(164)*Cd(2597)*Ce(1079) + 
     &   MTy*Cd(47)**2*UCStopC(1,3)

        Ce(1081) = -2*Ce(1080) + 
     &   AtC*Cd(164)*Cd(2597)*Ce(1077)*UStop2(1,1)

	Ce(1097) = (Ce(984) + MStop2(3))/MStop2(3)

        Ce(1100) = (-(MGlmT2*Cd(47)*Ce(1099)) + 
     &     4*Cd(10112)*Cd(12370)*MStop2(3))/MStop2(3)

	Ce(1106) = 8*Cd(47)**2*Cd(10208) + Cd(12371)*Ce(1105)

	Ce(1114) = 2*At + M_3*Ce(1113)

	Ce(1121) = M_3*Ce(945) - 2*At*MGlmT2*Cd(9955)*UStop2(1,2)

	Ce(113) = 3*C2A*Cd(13) - S2A*Ce(112)

	Ce(1134) = Cd(47)*Ce(1133) + MGlmT2*Ce(975)*MStop2(3)

	Ce(114) = 1/2.D0*(Ce(111) + 2*MTy*Ce(113))/(MTy*MZ2*SW2)

        Ce(1142) = (MGlmT2*Cd(47)*Ce(1099) + 
     &     4*Cd(152)**2*Cd(11735)*MStop2(3))/MStop2(3)

	Ce(1152) = M_3C*Ce(1149) + Ce(1151)*UCStopC(1,3)

	Ce(117) = 3*C2A*Cd(13) + S2A*Ce(116)

	Ce(1167) = MGlmT2*Ce(1166) + Ce(1165)*MStop2(1)

	Ce(118) = 1/2.D0*(Ce(111) - 2*MTy*Ce(117))/(MTy*MZ2*SW2)

	Ce(1177) = Cd(152)*Cd(4095) + Ce(910)*MGlpTmSt2(1)*MStop2(1)

	Ce(1178) = Cd(179)*Cd(4190) + Ce(910)*MGlpTmSt2(2)*MStop2(2)

	Ce(1191) = Cd(50) - Ce(1190)*MGlpTmSt2(2)

	Ce(1200) = 2*Ce(1199) + MStop2(3)

	Ce(1209) = -(Ce(1208)*MGlpTmSt2(2)) + Cd(50)*UStop2(1,1)

        Ce(1225) = (2*Cd(47)*Cd(50)*Ce(1224) + 
     &     Cd(10086)*Cd(12045)*MStop2(3))/MStop2(3)

        Ce(1227) = (2*Cd(47)*Cd(50)*Ce(1226) + 
     &     Cd(10039)*Cd(12045)*MStop2(3))/MStop2(3)

	Ce(1233) = M_3*Ce(1230) + M_3C*Ce(1232)

	Ce(1234) = Cd(50)*Ce(1233) + 2*MUE2*Cd(47)*Cr(5)*Cr(75)

	Ce(1239) = At*Cd(12063) - M_3*Ce(1236)

	Ce(1260) = MGlmT2*Ce(1258) + Ce(1259)*MStop2(2)

	Ce(1267) = (-Ce(1262) + MGlmT2*Ce(1266)*MStop2(3))/MStop2(3)

	Ce(1272) = -8*(Ce(1271) - 2*Cd(9969)*MStop2(2))

        Ce(1286) = (Cd(251)*Ce(1285)*MStop2(3) + 
     &     12*MGlmT2*Cd(50)*Cd(9955)*UStop2(1,1)**2)/MStop2(3)

        Ce(1292) = (-(MGlmT2*Cd(50)*Ce(1291)) + 
     &     4*Ce(1289)*MStop2(3))/MStop2(3)

        Ce(1295) = (8*(Ce(1294) + 
     &       12*Cd(9955)*MStop2(1)*UStop2(1,1)**2))/Cd(50)

	Ce(1309) = MGlmT2*Ce(1258) + Ce(1308)*MStop2(2)

	Ce(1310) = 4*Cd(50)*Ce(1309)*MStop2(2) - Ce(1307)*MStop2(3)

        Ce(1327) = -(MGlmT2*Cd(50)*Ce(1326)) + 
     &   4*Cd(10015)*Cd(10602)*MStop2(3)*UStop2(1,1)

	Ce(1340) = MGl2*Cd(9948)*UStop2(1,1) + Ce(1339)*UStop2(1,2)

        Ce(1342) = Ce(1341)*UCStopC(1,3) + 
     &   2*M_3C*MTy2*UStop2(1,1)*UStop2(1,2)

	Ce(1343) = -2*Ce(1342) + AtC*Ce(1340)*UStop2(1,1)

	Ce(1350) = Cd(12481) + Cd(50)*Ce(1349)

        Ce(136) = MUE2*Ce(127) + 
     &   2*At*AtC*Cr(75)**2*MStop2(3)*UStop2(1,2)

	Ce(137) = AtC**2*Ce(128) + 2*MUE**2*Cr(75)**2

        Ce(1369) = (4*(-2*MGlmT2*Cd(50)*Ce(1368) + 
     &       Ce(1367)*MStop2(3)))/(Cd(50)**2*MStop2(3))

	Ce(1372) = -(Cd(12483)*Ce(1371)) + 16*Cd(12482)*Cr(28)*MStop2(2)

	Ce(1375) = MGlmT2*Ce(1373) + Ce(1374)*MStop2(2)

	Ce(1379) = MGlmT2*Ce(1377) - Ce(1378)*MStop2(2)

        Ce(1380) = 4*Cd(50)*Ce(1379)*MStop2(2) - 
     &   2*Ce(1376)*MStop2(3)*UStop2(1,1)

	Ce(139) = At**2*Ce(133) - MUEC**2*Ce(135)

	Ce(1384) = Cd(12483)*Ce(1383) + 16*Cd(12482)*MStop2(2)

	Ce(1392) = AtC*Cd(9983) - 4*Ce(1391)

	Ce(1416) = 8*(Ce(1412) - Ce(1415))

        Ce(1425) = 2*MUE2*Cd(10020)*UCStop(1,3) + 
     &   At**2*Ce(1424)*UCStopC(1,3)

        Ce(1437) = (8*(Ce(1436) + 
     &       12*Cd(9955)*MStop2(1)*UStop2(1,1)**2))/Cd(50)

        Ce(1456) = Ce(394)*Ce(1455) - 
     &   2*Ce(392)*Ce(1442)*MStop2(1)*MStop2(3)

        Ce(1468) = Ce(388)*Ce(1467) + 
     &   2*Ce(389)*Ce(1459)*MStop2(1)*MStop2(3)

        Ce(1483) = Ce(403)*Ce(1482) + 
     &   2*Ce(401)*Ce(1469)*MStop2(2)*MStop2(3)

        Ce(1498) = Ce(398)*Ce(1497) - 
     &   2*Ce(399)*Ce(1486)*MStop2(2)*MStop2(3)

	Ce(150) = Ce(148) - Cd(5643)*UStop2(1,1)

	Ce(1504) = AtC*Ce(328) - MUE2*UUStop(2,2)

	Ce(1505) = 2*Ce(1503)*UStop(2,2) + Ce(1504)*UUStop(1,1)

	Ce(1506) = 2*AtC*MTy*Cr(91) + Ce(35)*Ce(36)*UUStop(3,1)

	Ce(1507) = 4*MTy2*UUStop(1,2)*UUStop(2,2) + Ce(1506)*UUStop(3,2)

	Ce(1511) = Ce(35)*Ce(36)*UUStop(1,1) - 4*MTy2*UUStop(1,2)

	Ce(1512) = Ce(1510)*UUStop(2,2) + Ce(1511)*UUStop(2,3)

	Ce(151) = Ce(149) - Ce(150)*UStop2(1,1)

	Ce(1515) = 2*At*MTy*UUStop(1,1) + Ce(186)*Ce(187)*UUStop(1,3)

	Ce(1523) = Ce(1520) - AtC*Ce(1514)*UCStop(1,3)*UStop(1,1)

	Ce(1532) = Ce(1531)*UStopC(2,1) + Ce(1529)*UStopC(2,2)

	Ce(1546) = Ce(1544)*UStopC(2,1) - Ce(1545)*UStopC(2,2)

	Ce(155) = Ce(153) + Ce(152)*UStop2(1,1)

	Ce(1562) = Ce(1559)*UUStopC(3,2) - Ce(1561)*UUStopC(3,4)

	Ce(1564) = 2*Ce(1562) - AtC*Ce(1563)

        Ce(1567) = (Ce(1564)*MStop2(2) + 
     &     4*Ce(1566)*MStop2(3)*UStop2(1,2))/(MStop2(2)*MStop2(3))

        Ce(1571) = 4*Ce(1570)*UStop2(1,2) + 
     &   8*MUEC**2*Cd(9948)*UStop2(1,1)*UUStop(1,2)

	Ce(1572) = (Ce(1571) + Ce(1567)*MStop2(2))/(MStop2(1)*MStop2(2))

	Ce(157) = 2*Ce(155)*UStop2(1,1)**2 - Ce(156)*UStop2(1,2)**2

	Ce(159) = Ce(158) - 2*Cd(535)*MStop2(4)*UStop2(1,1)

        Ce(1602) = Ce(1601)*UStop(2,2) - 
     &   2*MUE**2*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,1)

	Ce(160) = Ce(159)*UStop2(1,1) - 2*Ce(143)*UStop2(1,2)

	Ce(1608) = 2*At*MTy*Cd(3852)*UUStop(1,2) + Ce(1607)*UUStop(1,3)

        Ce(1630) = Ce(1619) + 
     &   4*AtC*MStop2(2)**2*UCStop(1,3)**2*UUStop(1,1)

	Ce(1631) = Ce(1562)*MStop2(2) + Ce(1566)*MStop2(3)*UStop2(1,2)

	Ce(1632) = AtC*Ce(1630) - 4*Ce(1631)*MStop2(3)

	Ce(163) = AtC**2*Ce(151) + 2*MUE**2*Cd(833)*Cr(75)

        Ce(1635) = (2*Ce(1634)*MStop2(1)*MStop2(2) + 
     &     Ce(1569)*MStop2(3)**2*UStop2(1,2))/
     &   (MStop2(1)*MStop2(2)*MStop2(3)**2)

	Ce(164) = 2*AtC*MTy*Ce(156)*Cr(75) + Ce(163)*UCStop(1,3)

	Ce(1648) = 2*Ce(1644) + Ce(1647)*UUStopC(1,2)

        Ce(165) = -2*MTy*Ce(157)*UCStopC(1,3) + 
     &   4*AtC*Cd(833)*Cr(75)*UStop2(1,1)*UStop2(1,2)

        Ce(1653) = 2*(Ce(186)*Ce(187)*UCStopC(1,3)**2 + 
     &     MUE2*UStop2(1,1)**2)

        Ce(1657) = Ce(1656)*MStop2(3)**2*UStop2(1,1) + 
     &   AtC*Ce(1651)*UStopC(1,1)

        Ce(1674) = 2*Ce(1664)*UCStop(3,1)*UStopC(2,1) - 
     &   Ce(1665)*UUStopC(2,1)

	Ce(1679) = 2*Ce(1643) - Ce(1641)*UStop2(1,1)

	Ce(1680) = Ce(1648) + Ce(1679)*UStop2(1,1)

	Ce(1688) = Ce(1546)*MStop2(1) + Ce(1662)*MStop2(3)*UStop2(1,2)

	Ce(169) = 4*C2A*Cd(760) + S2A*Ce(168)

	Ce(1699) = Ce(1698)*MStop2(1) + Cd(13422)*MStop2(3)

        Ce(1700) = Ce(1697)*UStop2(1,1) + 
     &   2*At*MTy*Ce(1699)*UCStopC(1,3)*UStop2(1,2)

	Ce(1703) = Ce(1702)*MStop2(1) + Cd(13422)*MStop2(3)*UStop2(1,2)

	Ce(1709) = 5*Ce(1644) + 2*Ce(1647)*UUStopC(1,2)

	Ce(1711) = Ce(1687) + 2*Ce(1641)*MStop2(3)*UStop2(1,1)

	Ce(171) = C2A*Cd(755) + 2*S2A*Ce(170)*Cr(75)

	Ce(1717) = 2*Ce(1546)*MStop2(1) + Ce(1662)*MStop2(3)*UStop2(1,2)

        Ce(1728) = Cd(2365)*Cd(9948)*MStop2(1)*MStop2(3) - 
     &   Ce(1727)*UStop2(1,2)

	Ce(1733) = -4*Ce(1732) + AtC*Ce(1720)*UStopC(1,1)

	Ce(1741) = 4*AtC*MTy + Ce(35)*Ce(36)*UCStop(1,3)

	Ce(1742) = Ce(1741)*UCStop(1,3) + 8*MTy2*UStop2(1,2)

	Ce(1748) = Ce(1716) + Ce(1641)*MStop2(3)*UStop2(1,1)

        Ce(1755) = Cd(10925)*UCStop(3,3)*UStop2(1,2) + 
     &   Ce(1754)*UUStop(1,1)

        Ce(1756) = 4*Ce(1755)*MStop2(2) - 
     &   Cd(11607)*MStop2(3)*UUStop(1,1)

	end
