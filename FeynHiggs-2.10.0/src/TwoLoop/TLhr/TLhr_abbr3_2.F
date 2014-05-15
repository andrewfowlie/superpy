* TLhr_abbr3_2.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr3_2
	implicit none

#include "TLhr.h"


        Co(862) = At*UCStop(3,1)*UStop2(1,2)*UUStop(2,2) - 
     &   MUE*TB*Cd(970)*UUStop(3,1)

        Co(863) = -(MUEC*Cd(12701)*UCStopC(1,3)**2) + 
     &   At*MUE*UStop2(1,2)**2

	Co(864) = MUE*Co(862) + Co(863)*UUStop(3,2)

	Co(85) = Co(84) - 2*MStop2(4)*UStop2(1,1)

	Co(865) = MUEC*Cd(14311)*UCStopC(1,3) - MUE*Cd(970)*UStop(1,1)

	Co(866) = MUE*UUStop(3,1) - 3*MUEC*UUStop(3,2)

        Co(867) = Co(865)*MStop2(3) + 
     &   Co(866)*MStop2(2)*UCStop(3,1)*UStopC(1,1)

	Co(868) = MUE*TB*UUStop(3,1) + 3*At*UUStop(3,2)

        Co(869) = MUEC*Cd(12701)*UUStop(1,2)*UUStop(2,2) - 
     &   MUE*Co(868)*UUStop(3,1)

	Co(871) = Cd(9948)*Co(855) + Co(861)*UUStopC(2,2)

	Co(874) = 2*Cd(970)*Co(641) + 4*MUEC*Cd(12701)*UCStopC(1,3)**2

	Co(86) = Cd(599)*MStop2(2) + Co(85)*UStop2(1,1)*UStop2(1,2)

	Co(875) = AtC*Cr(91)*UStop(2,1) + 4*At*UStop(1,2)*UUStop(2,2)

	Co(876) = Co(875) + 4*MUEC*TB*UStop(1,2)*UUStop(2,2)

	Co(877) = -(MUEC*Co(876)) + At*MUE*Cr(91)*UStop(2,1)

	Co(878) = AtC*Cd(10925)*UUStop(2,1) + 2*At*Cd(11106)*UUStop(2,2)

	Co(879) = Co(878) + 2*MUEC*TB*Cd(11106)*UUStop(2,2)

	Co(880) = -(MUEC*Co(879)) + At*MUE*Cd(10925)*UUStop(2,1)

	Co(881) = 3*Co(641)*UUStop(3,3) - Co(880)*UUStopC(2,1)

        Co(882) = Cd(9948)*Co(874)*MStop2(3)*UStop2(1,1) + 
     &   2*Co(881)*MStop2(2)*UUStopC(3,3)

        Co(883) = -2*Co(882)*MStop2(3) + 
     &   4*Co(877)*MStop2(2)**2*UCStop(3,1)*UStopC(2,2)*
     &    UUStopC(1,1)

        Co(888) = 2*Co(789)*MStop2(2)*UStopC(1,1) + 
     &   4*Co(836)*MStop2(3)*UStopC(2,1)

	Co(889) = MUE*Cd(1136) + 2*MUEC*UCStopC(1,3)**2

        Co(892) = -(Cd(13422)*MStop2(3)) + 
     &   2*MStop2(2)*UStop2(1,1)*UStop2(1,2)

        Co(893) = 2*MUE**2*Co(892)*UUStop(2,1) + 
     &   4*MUEC**2*MStop2(1)*UCStopC(1,3)**2*UUStop(2,2)

	Co(7) = ((TB*Co(2) - Co(19))*MStop2(3)**2)/TB

        Co(88) = AtC**2*Cd(609)*UCStop(1,3)**2 - 
     &   At**2*Cd(611)*UCStopC(1,3)**2

	Co(897) = UStop(1,2) - UCStop(3,1)*UStop(2,2)*UStop2(1,2)

        Co(903) = 2*Co(857)*MStop2(2)*UCStopC(3,3) - 
     &   3*Co(860)*MStop2(3)*UStopC(2,1)

        Co(904) = -(MUEC*Cd(14311)*UCStopC(1,3)) + 
     &   MUE*Cd(2364)*UStop(1,1)

        Co(907) = At*UCStop(3,1)*UStop2(1,2)*UUStop(2,2) - 
     &   MUE*TB*Cd(2364)*UUStop(3,1)

	Co(908) = MUE*Co(907) + Co(863)*UUStop(3,2)

	Co(912) = Co(861) - Co(903)

	Co(916) = Cd(1136)*Co(641) + 2*MUEC*Cd(12701)*UCStopC(1,3)**2

	Co(917) = 3*Co(874) + 2*Co(916)

	Co(918) = Co(641)*UUStop(3,3) - Co(880)*UUStopC(2,1)

	Co(920) = 2*Co(874) - Co(917)

        Co(925) = MUE**2*Cd(1136)*UUStop(2,1) + 
     &   MUEC**2*UCStopC(1,3)**2*UUStop(2,2)

	Co(926) = -(MUEC*UCStopC(1,3)**2) + 2*MUE*UStop2(1,2)**2

	Co(927) = AtC*Co(889)*UUStop(2,1) - At*Co(926)*UUStop(2,2)

	Co(928) = (TB*Co(925) + Co(927))/(TB*MStop2(2)**2)

	Co(933) = AtC*Co(897) - 2*At*UCStopC(1,3)**2*UStop(1,2)

	Co(934) = Co(933) - 2*MUEC*TB*UCStopC(1,3)**2*UStop(1,2)

        Co(92) = MUE**2*Cd(618)*UCStop(1,3)**2 - 
     &   MUEC**2*Cd(619)*UCStopC(1,3)**2

	Co(935) = At*MUE*Co(897) - MUEC*Co(934)

	Co(936) = Co(908) - AtC*Co(904)*UStop(2,1)

	Co(939) = Co(912)*UStopC(1,2) - Co(860)*MStop2(3)*UUStopC(3,4)

        Co(945) = At*UUStop(1,3)*UUStop(2,2) + 
     &   2*MUE*TB*UUStop(1,1)*UUStop(2,3)

        Co(947) = At*MUE*UUStop(1,3)*UUStop(2,1) - 
     &   2*MUEC*Cd(12701)*UUStop(1,2)*UUStop(2,3)

        Co(949) = 2*MUEC*Cd(12701)*UUStop(1,3)*UUStop(2,2) - 
     &   At*MUE*UUStop(1,1)*UUStop(2,3)

        Co(950) = 2*MUE*TB*UUStop(1,3)*UUStop(2,1) + 
     &   At*UUStop(1,2)*UUStop(2,3)

	Co(951) = MUE*TB*UUStop(3,1) + At*UUStop(3,2)

        Co(952) = MUEC*Cd(12701)*UUStop(1,2)*UUStop(2,2) - 
     &   MUE*Co(951)*UUStop(3,1)

        Co(954) = UStopC(2,2)*UUStop(1,2)*UUStop(2,1) - 
     &   Cd(9948)*UStop(1,1)*UUStop(3,4)

        Co(955) = 2*MUE*UUStop(1,3)*UUStop(2,1) - 
     &   MUEC*UUStop(1,2)*UUStop(2,3)

        Co(957) = MUEC*UUStop(1,3)*UUStop(2,2) - 
     &   2*MUE*UUStop(1,1)*UUStop(2,3)

	Co(961) = At*MUE*UUStop(3,1) - 2*MUEC*Cd(12701)*UUStop(3,2)

        Co(962) = 2*MUEC*Cd(12701)*Cr(91)*UStopC(2,1) + 
     &   Co(961)*UStopC(2,2)

	Co(963) = Co(788)*UStopC(2,1) - 2*MUE*TB*Cr(91)*UStopC(2,2)

	Co(95) = Cd(169) - UStop2(1,1)

	Co(964) = 2*Co(808)*UStopC(2,1) - MUE*Co(809)*UStopC(2,2)

        Co(965) = 8*Cd(14040)*MStop2(1)*MStop2(2)*UUStop(2,1) - 
     &   4*At*MStop2(3)**2*UUStop(1,1)*UUStop(2,2)

        Co(966) = MUE*Co(965) - 
     &   8*MUEC*Cd(12701)*MStop2(1)*MStop2(2)*UUStop(1,2)*
     &    UUStop(2,2)

        Co(967) = 8*MUE*TB*MStop2(1)*MStop2(2)*UUStop(2,1) + 
     &   At*Cd(5519)*UUStop(2,2)

        Co(968) = MUE*TB*MStop2(3)**2*UUStop(1,1) - 
     &   At*MStop2(1)*MStop2(2)*UUStop(1,2)

	Co(969) = 4*Co(967)*UUStop(1,3) - 8*Co(968)*UUStop(2,3)

	Co(970) = 2*MUE*TB*Cd(11657)*UUStopC(2,2) - Co(969)*UUStopC(2,3)

	Co(971) = MUE*Co(970) + Co(966)*UUStopC(2,1)

        Co(972) = At*MUE*MStop2(1)*MStop2(2)*UUStop(2,1) + 
     &   MUEC*Cd(12701)*MStop2(3)**2*UUStop(2,2)

        Co(973) = At*MUE*Cd(5519)*UUStop(1,1) - 
     &   8*MUEC*Cd(12701)*MStop2(1)*MStop2(2)*UUStop(1,2)

	Co(96) = Co(95)*UStop2(1,1) - 2*UStop2(1,2)

	Co(974) = 8*Co(972)*UUStop(1,3) + 4*Co(973)*UUStop(2,3)

        Co(975) = MUE*TB*MStop2(3)**2*UUStop(2,1) - 
     &   At*MStop2(1)*MStop2(2)*UUStop(2,2)

        Co(976) = 8*MUE*TB*MStop2(1)*MStop2(2)*UUStop(1,1) + 
     &   At*Cd(5519)*UUStop(1,2)

	Co(977) = 8*Co(975)*UUStop(1,3) - 4*Co(976)*UUStop(2,3)

	Co(978) = MUE*TB*Cd(5519)*UUStop(3,1) + At*Cd(5945)*UUStop(3,2)

        Co(979) = 4*MUEC*Cd(5519)*Cd(12701)*UUStop(1,2)*
     &    UUStop(2,2) - 4*MUE*Co(978)*UUStop(3,1)

        Co(981) = At*MUE*Cd(10925)*UStop(2,1) - 
     &   2*MUEC*Cd(12701)*UStop(1,2)*UUStop(2,2)

        Co(983) = At*Cd(11106)*UStop(2,2) + 
     &   2*MUE*TB*UStop(1,1)*UUStop(2,1)

	Co(97) = MUE*TB*Co(96) + 2*AtC*Cr(75)

        Co(984) = 2*MUE*TB*Cr(91)*UCStop(1,3)*UStopC(2,2) + 
     &   Co(983)*UUStopC(2,3)

	Co(985) = MUE*Co(984) + Co(961)*UCStopC(1,3)*UStopC(2,1)

        Co(986) = 4*Cd(11680)*UUStop(2,1) - 
     &   8*MUEC*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,2)

        Co(987) = 8*MStop2(1)*MStop2(2)*UUStop(3,3) + 
     &   4*Cd(5519)*UUStop(3,4)

        Co(989) = 8*Cd(10946)*MStop2(1)*MStop2(2)*UUStop(2,1) + 
     &   4*MUEC*MStop2(3)**2*UUStop(1,1)*UUStop(2,2)

        Co(990) = 8*MUE*MStop2(1)*MStop2(2)*UUStop(2,1) - 
     &   MUEC*Cd(5519)*UUStop(2,2)

        Co(991) = MUE*MStop2(3)**2*UUStop(1,1) + 
     &   MUEC*MStop2(1)*MStop2(2)*UUStop(1,2)

	Co(992) = 4*Co(990)*UUStop(1,3) - 8*Co(991)*UUStop(2,3)

        Co(995) = MUEC*UCStop(3,1)*UStop(2,2)*UStop2(1,2) + 
     &   2*MUE*Cr(91)*UCStop(1,3)*UStopC(2,2) - 
     &   Co(994)*UUStopC(2,3)

	Co(1038) = -A0delGl + MGl2 + MGl2*Log(MGl2/Q)

	Co(1040) = A0delT*MGl2 - MTy2*Co(1038) + Co(1039)*Cr(66)

	Co(1052) = -A0delT - MTy2 + 2*Cr(66)

	Co(1072) = Cd(50)*(4*Ce(411) - Co(1049)) + 2*Co(1062)

        Co(1073) = 2*Cd(50)**2*(2*Cd(47)*Ce(411) + Co(1057)) + 
     &   Cd(47)**2*Co(1072)

        Co(1075) = -2*Cd(50)*Cd(9947)*MGlpTmSt2(1) + 
     &   Cd(47)*(Co(1047) - 2*Cd(9949)*MGlpTmSt2(2))

	Co(1083) = 2*MGlmT2*Cd(47)*Co(1037) + Co(1069)

        Co(108) = -(Co(107)*MStop2(2)) + 
     &   2*Cd(286)*UStop2(1,1)**2 - 4*MStop2(1)*UStop2(1,2)**2

        Co(1097) = Co(364)*Co(1066) + 
     &   4*AtC*Cd(47)*Cd(50)*
     &    (T134(MGl2,MTy2,MStop2(1)) - T134(MGl2,MTy2,MStop2(2)))*
     &    UCStop(1,3)**2

	Co(1103) = 2*Co(1042) - 2*Co(1043) - Co(1056)

	Co(1120) = 2*Co(1050) - 2*Co(1051) - Co(1061)

        Co(1132) = (4*Ce(411) - Co(1034))*UStop2(1,1) + 
     &   Co(1131)*UStop2(1,2)

        Co(1134) = -4*Ce(411)*Cr(28) + Co(1034)*Cr(28) + 
     &   Co(1131)*Cr(75)*UStop2(1,1)

        Co(1149) = (Co(1037) - 2*Cd(9947)*MGlpTmSt2(1))*
     &    UStop2(1,1) + Co(1148)*UStop2(1,2)

        Co(1151) = -(Co(1131)*Cr(28)) + 
     &   (4*Ce(411) - Co(1034))*Cr(75)*UStop2(1,1)

        Co(1167) = Cd(47)**2*Co(1072) + 
     &   Cd(50)**2*(2*Co(1057) + Cd(47)*Co(1148))

        Co(1171) = Cd(9948)*Cd(9950)*Co(1058) + 
     &   MGlmT2*Cd(47)*(Cd(47)*Cd(9948)*Co(1034) + 
     &      Cd(9950)*Co(1034)*UStop2(1,2))

        Co(1180) = -(MGlmT2*Cd(47)*
     &      (Cd(9950)*Co(1034) + Cd(47)*Co(1037))) + 
     &   Cd(9950)*Co(1058)

        Co(1182) = Co(372)*Co(1066) + 
     &   4*At*Cd(47)*Cd(50)*
     &    (T134(MGl2,MTy2,MStop2(1)) - T134(MGl2,MTy2,MStop2(2)))*
     &    UCStopC(1,3)**2

        Co(117) = (TB*Co(116) + 
     &     4*MTy*Cd(12836)*Cr(75)*MStop2(3) - 2*Co(99)*MStop2(4))/
     &   (TB*MStop2(2)*MStop2(3)**2)

        Co(1249) = At*Co(1247) - MTy*Co(1248) - 
     &   MUEC*TB*Cd(11165)*UStop(1,2)

        Co(1258) = Co(1256)*UCStop(3,3)*UStop2(1,2) - 
     &   2*MTy*UStop2(1,2)**2*UUStop(1,1) - Co(1257)*UUStopC(2,2)

	Co(1261) = 3 + UStop2(1,1)**2 - 5*UStop2(1,2)**2

        Co(126) = 2*MTy*TB*Co(123) + At*Co(125) + 
     &   4*AtC*Co(122)*Cr(75)**2*MStop2(3)

        Co(1274) = 4*MUEC*TB*Co(1273)*MStop2(3) + 
     &   At*(2*MUE*TB*Cd(802)*UUStop(1,1) + 
     &      4*MTy*MStop2(4)*UUStop(1,3))

        Co(1278) = Co(1277)*UCStop(1,3) + 
     &   2*At*UStop2(1,2)*(At*UUStop(1,2) + 2*MTy*UUStop(1,3))

	Co(1281) = Co(1273) + Co(1280)

        Co(1287) = MUE*TB*UUStop(1,1) + 2*MUEC*TB*UUStop(1,2) - 
     &   2*MTy*UUStop(1,3)

        Co(1298) = 4*Co(1296)*MStop2(1) - 
     &   2*Co(1297)*MStop2(1)**2 - 
     &   6*Cd(9794)*MStop2(3)**2*UStop(1,2)

	Co(1300) = 3 - 6*UStop2(1,1)**2 - UStop2(1,2)**2

        Co(129) = UStop2(1,1)*
     &    (Cd(839) - 2*Cd(835)*UStop2(1,1)) - 
     &   2*Cd(836)*MStop2(4)*UStop2(1,2)

	Co(1304) = 3 + UStop2(1,1)**2 - 2*UStop2(1,2)**2

        Co(1316) = MTy*(2*MUE*TB*UUStop(3,1) - 
     &      2*Cd(12701)*UUStop(3,2)) + At*MUE*TB*UUStop(3,4)

        Co(1320) = 2*MTy*Co(1319)*Cr(58) - 
     &   At*Co(1318)*UUStop(1,2) - MUEC*TB*Cd(11036)*UUStop(3,2)

        Co(1322) = MTy*(2*Co(1321) + 
     &      4*MUEC*TB*UUStop(1,2)*UUStop(2,2)) - 
     &   At*MUE*TB*Cr(91)*UUStop(3,2)

        Co(1332) = 2*MTy*Cd(10995) + 
     &   4*MUE*TB*UUStop(1,3)*UUStop(2,1) + 
     &   MUEC*TB*Cr(91)*UUStop(3,2)

        Co(1333) = Co(1330)*UStop2(1,2)*UStopC(2,1) + 
     &   Co(1331)*Cr(91)*UUStopC(2,2) - Co(1332)*UUStopC(2,3)

        Co(1335) = 2*MUE*TB*Cd(11106)*UStop(1,1) + 
     &   2*MTy*Cd(12902)*UStop(1,2) + 
     &   MUEC*TB*UStop(2,2)*UUStop(1,2)

        Co(1344) = 2*MUEC*TB*Co(1273) + 
     &   At*(5*MUE*TB*UUStop(1,1) + 2*MTy*UUStop(1,3))

	Co(1357) = 11 + 2*UStop2(1,1)**2 - 10*UStop2(1,2)**2

        Co(1363) = UStop2(1,1)*(UStop2(1,1) - 6*UStop2(1,2)) + 
     &   UStop2(1,2)**2

        Co(1414) = Co(1413)*UStop(2,2) + 
     &   AtC*(Co(1411)*UStop(2,1) + 2*AtC*UStop(1,1)*UUStop(2,1))

        Co(1418) = -(MTy*Co(671)) + At*Co(1417) - 
     &   MUEC*TB*Cd(11165)*UStop(2,2)

        Co(1419) = -(Co(1418)*UStop(2,2)) + 
     &   AtC*(Co(1416)*UStop(2,1) + AtC*UUStop(2,1)*UUStop(3,1))

        Co(1422) = At*Co(1421) + 2*MTy*Co(1420)*UStop(2,1) + 
     &   2*MUEC*TB*Cd(11391)*UStop(2,2)

	Co(1423) = MUE*TB + AtC*(1 - 3*UStop2(1,2)**2)

        Co(1431) = MUE*TB*Co(1430) + 
     &   MTy*(Cd(10918)*UUStop(3,3) + Cd(11502)*UUStop(3,4))

        Co(1438) = -2*MUE*TB*UCStop(1,3)*UStop2(1,1) + 
     &   MTy*(1 + UStop2(1,1)**2)

	Co(1444) = At + MUEC*TB*(1 - 3*UStop2(1,1)**2)

        Co(1463) = 2*Co(1461)*UUStop(2,1) - 
     &   2*MTy*UUStop(1,1)*UUStop(2,2) + Co(1462)*UUStop(3,2)

        Co(1464) = Co(1330)*UStop2(1,2)*UStopC(2,1) + 
     &   Co(1331)*Cr(91)*UUStopC(2,2) - Co(1463)*UUStopC(2,3)

	Co(1469) = MUE*TB + AtC*(1 - UStop2(1,2)**2)

        Co(1480) = 2*UUStop(1,2)*UUStop(2,1) - 
     &   2*UUStop(1,1)*UUStop(2,2) - 4*UUStop(3,1)*UUStop(3,2)

        Co(1481) = MTy*Co(1480) - 
     &   TB*(4*MUE*UUStop(1,1)*UUStop(2,3) + 
     &      MUEC*Cr(91)*UUStop(3,2))

	Co(1487) = 5 + UStop2(1,1)*(-5*UStop2(1,1) + 2*UStop2(1,2))

	Co(1495) = At - MUEC*TB*(-1 + UStop2(1,1)**2)

        Co(1511) = MUE*TB*Co(1478) + 
     &   MTy*(Cd(11316)*UUStop(3,3) + 3*Cd(10918)*UUStop(3,4))

        Co(1517) = 2*MUEC*MStop2(2)*UCStopC(1,3)*UStop2(1,1)*
     &    UUStop(3,2) + MUE*
     &    (Co(1514)*UUStop(3,3) - Co(1516)*UUStop(3,4))

	Co(1533) = Co(1673) + Ce(1514)*UStop(1,1)*UStop2(1,1)

        Co(1538) = 2*MTy*MUE*TB*Cr(91)*UCStop(3,1)*UStop2(1,1) + 
     &   Co(1535)*UUStopC(2,1) + Co(1537)*UUStop(2,2)*UUStopC(2,3)

        Co(1543) = At*Co(1542) + 
     &   2*MUEC*TB*Co(1541)*UUStop(1,2) + 
     &   2*MTy*Cd(14040)*UUStop(2,1)

        Co(1544) = 2*Cd(10987)*Cd(12701)*UUStop(3,1)*
     &    UUStopC(2,1) + Co(1540)*UUStopC(2,2) - 
     &   Co(1543)*UUStopC(2,3)

        Co(1547) = At*Co(1546) + 2*MTy*Co(1318)*UUStop(1,1) + 
     &   2*MUEC*TB*Co(1545)*UUStop(2,2)

        Co(1558) = 2*TB*Cd(9948)*Cd(11394)*UCStop(3,3) + 
     &   Co(1555)*UCStop(3,1)*UStop2(1,2) - Co(1557)*UUStopC(2,2)

        Co(1562) = Co(1559)*UUStop(1,1)*UUStopC(2,1) + 
     &   2*Co(1327)*UUStop(3,2)*UUStopC(2,2) - 
     &   Co(1561)*UUStopC(2,3)

        Co(1563) = Co(1554)*UStopC(2,2)*UUStopC(1,1) + 
     &   Co(1562)*UUStopC(1,2) + Co(1558)*UUStopC(1,3)

        Co(1568) = -(At*MUE*TB*Cd(10925)*UStop(1,2)) + 
     &   MTy*(2*Co(1250) + 6*MUEC*TB*UStop(2,2)*UUStop(1,2))

        Co(1574) = Cd(11107)*MStop2(2)*MStop2(3) + 
     &   Co(1573)*UStop(1,1) - MStop2(3)**2*UCStop(3,1)*UStop(2,1)

        Co(1583) = At*Co(1581) + 
     &   MTy*(-2*Co(1319)*Co(1580) + 8*MUEC*TB*Co(1582))

        Co(1586) = 2*MTy*Co(1584) - 4*TB*Co(1585)*UUStop(2,2) + 
     &   8*At*MUE*TB*MStop2(1)*MStop2(2)*UUStop(1,2)*UUStop(2,3)

        Co(1587) = 8*MTy*MUE*TB*Cd(3852)*Cr(91)*UCStop(3,1)*
     &    UStop2(1,1) + Co(1583)*UUStopC(2,1) - 
     &   Co(1586)*UUStopC(2,3)

        Co(1588) = 2*MTy*Co(951)*UStop(2,1) + 
     &   At*Co(1566)*UStop(2,2) + 2*MUEC*TB*Cd(11330)*UUStop(3,2)

        Co(1592) = 2*MTy*MUE*TB*Cd(11161)*UStop2(1,1)*
     &    UStopC(2,2) - Co(1588)*UUStopC(2,1) + 
     &   Co(1591)*UUStopC(2,3)

        Co(1602) = MTy*Co(1600) - 4*TB*Co(1601)*UUStop(2,2) + 
     &   2*At*MUE*TB*Cd(5519)*UUStop(1,2)*UUStop(2,3)

        Co(1607) = Co(1599)*UUStopC(2,1) + 
     &   2*Co(1602)*UUStopC(2,2) - 2*Co(1606)*UUStopC(2,3)

        Co(1612) = At*Co(1611) + 2*MTy*Co(1610)*UStop(1,1) - 
     &   2*MUEC*TB*Co(1609)*UUStop(2,2)

        Co(1615) = -(MUEC*TB*Ce(1617)) + 8*MTy*Co(1614) + 
     &   8*MUE*TB*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,1)

        Co(1619) = 4*Co(1618)*MStop2(2)*MStop2(3)*UCStop(1,3) + 
     &   Cd(9948)*Co(1617)*UStop(1,1) + Co(1615)*UStopC(2,2)

        Co(1622) = 4*MUEC*TB*Co(1620) - 8*MTy*Co(1621) + 
     &   8*MUE*TB*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,1)

        Co(1627) = Co(1622)*UUStopC(2,1) + 
     &   2*Co(1623)*UUStopC(2,2) + Co(1626)*UUStopC(2,3)

        Co(1630) = Co(1237)*UCStop(3,1)*UStop2(1,2) - 
     &   2*Co(1628)*UUStopC(2,2) + Co(1629)*UUStopC(2,3)

        Co(1653) = At*Co(1651) + 
     &   MTy*(-2*Co(1319)*Co(1650) + 8*MUEC*TB*Co(1652))

	Co(1660) = At*Co(1656) - 2*MTy*Co(1657) + 2*MUEC*TB*Co(1659)

        Co(1662) = MTy*Co(1661)*UStop(1,1) + 
     &   At*Co(1247)*UStop(1,2) - 
     &   MUEC*TB*Cd(10980)*UStop(2,2)*UUStop(1,2)

	Co(165) = MUEC**2*SA*TB + At*(At*CA + MUEC*Co(139))

        Co(1665) = At*Co(1664) + 2*MTy*Co(1663)*UStop(1,1) + 
     &   2*MUEC*TB*Co(1295)*UStop(1,2)

        Co(1672) = Cd(11107)*MStop2(2) - 2*Ce(1514)*MStop2(3) - 
     &   Cr(9)*MStop2(1)*UStop(1,1)

        Co(1677) = 4*Co(1559)*UUStop(1,1)*UUStopC(2,1) + 
     &   8*Co(1327)*UUStop(3,2)*UUStopC(2,2) - 
     &   4*Co(1561)*UUStopC(2,3)

        Co(1678) = 2*(Co(1636)*Cr(9)*MStop2(1) - 
     &     2*Co(1558)*MStop2(3) + Co(1642)*MStop2(2)*UUStopC(3,2))

	Co(1691) = Co(1148) - Co(1131)*UStop2(1,2)

        Co(1695) = Co(1049) + Cd(9948)*Cd(9959)*Co(1147) + 
     &   Co(1694)*UStop2(1,1)

	Co(1699) = Co(1099) - Cd(50)*Co(1130)*UStop2(1,2)

        Co(1716) = MGlmT2*Cd(10370)*Co(1104) + 
     &   2*Cd(47)*Cd(10369)*
     &    (T134(MGl2,MTy2,MStop2(1)) - T134(MGl2,MTy2,MStop2(2)))

	Co(1725) = -(Cd(9972)*Co(1037)) + Co(1148)

	Co(1737) = Co(109) + MUE*TB*UCStop(1,3)

        Co(1738) = MGlmT2*Cd(10077)*Co(1104) + 
     &   Cd(47)*(-2*T134(MGl2,MTy2,MStop2(1)) + 
     &      2*T134(MGl2,MTy2,MStop2(2)))

        Co(1739) = MGlmT2*Co(1104) + 
     &   Cd(47)*(2*T134(MGl2,MTy2,MStop2(1)) - 
     &      2*T134(MGl2,MTy2,MStop2(2)))

	Co(1741) = Co(1099) + Cd(50)*Co(1130)*UStop2(1,2)

        Co(1742) = Cd(50)**2*(Co(1140) - 2*Cd(47)*Co(1177)) + 
     &   Cd(47)**2*Co(1741)

        Co(1743) = Cd(47)**2*Cd(50)*
     &    (Cd(50)*Co(1071) + Cd(9948)*Cd(9950)*Co(1130)) - 
     &   Cd(9950)*Co(1742)

	Co(1745) = Co(1710) + Co(1712)*MStop2(3)

	Co(1750) = Co(1071) - 2*Cd(9947)*MGlpTmSt2(1)

	Co(1754) = Co(1071) - 2*Ce(411)*Cr(85)

        Co(1770) = MGlmT2*Cd(9958)*(Co(1047) - Co(1049)) + 
     &   2*Cd(50)*Co(1769)*UStop2(1,2)

        Co(1776) = Cd(10658)*Co(1049) + 
     &   Cr(75)**2*(Co(1047) - 2*Cd(9949)*MGlpTmSt2(2))*UStop2(1,1)

        Co(1779) = Cd(9983)*
     &    (Co(1047) - 2*Cd(9949)*MGlpTmSt2(2)) - 
     &   Co(1049)*UStop2(1,1)**2

	Co(1784) = Cd(9948)*Co(1750) + Co(1753)

	Co(1785) = -2*Cd(9972)*Ce(411) + Co(1071)

        Co(1793) = Co(1058) - 
     &   2*MGlmT2*Cd(47)*Cd(10025)*Co(1104) + 4*Cd(47)**2*Co(1706)

        Co(1794) = Co(1063) - 
     &   2*MGlmT2*Cd(50)*Cd(10037)*Co(1067) + 4*Cd(50)**2*Co(1712)

        Co(1797) = (4*Ce(411) + Co(1034))*UStop2(1,1) - 
     &   Co(1131)*UStop2(1,2)

        Co(1800) = (12*Ce(411) - Co(1034))*UStop2(1,1) + 
     &   Co(1131)*UStop2(1,2)

        Co(1807) = Co(1063) - 
     &   2*MGlmT2*Cd(50)*Cd(10025)*Co(1067) + 4*Cd(50)**2*Co(1712)

	Co(1841) = Co(1071) - Co(1706)*MStop2(3)

        Co(1846) = Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952) - 
     &   2*Cd(11962)*Co(1040)

	Co(1848) = Co(1131) - Co(1148)*UStop2(1,1)

	Co(1851) = 4*Cd(50)*Ce(411) + Cd(47)*(4*Ce(411) - Co(1049))

	Co(1866) = Co(1148) - Co(1131)*Cr(85)

        Co(1874) = Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952) - 
     &   2*Cd(11962)*Co(1040)*UStop2(1,2)

	Co(1882) = -(Cd(9972)*Co(1131)) + Co(1148)

        Co(17) = Co(13)/TB + 
     &   (2*(TB*Co(16) + Co(4)*MStop2(4)))/(TB*MStop2(3))

	Co(1955) = CA*MUEC**2*TB + At*(-(At*SA) + MUEC*Co(1930))

        Co(195) = CA*MTy*Cr(91) + AtC*CA*UUStop(3,1) + 
     &   MUE*SA*UUStop(3,1)

        Co(1985) = 2*MTy*SA*Co(197) + AtC*Co(1984)*UStop(1,1) + 
     &   At*MUE*Co(1929)*UUStop(3,3)

        Co(1987) = 2*MTy*SA*Co(197) + AtC*Co(1986)*UStop(2,1) + 
     &   At*MUE*Co(1929)*UUStop(3,4)

        Co(198) = 2*CA*MTy*Co(197) + AtC*Co(196)*UStop(1,1) - 
     &   At*MUE*Co(138)*UUStop(3,3)

        Co(1991) = -(AtC*UCStop(1,3)) + 2*MTy*UStop2(1,1) + 
     &   Cd(9809)*UStopC(2,1)

        Co(2001) = -4*MTy*SA*Co(197) + AtC*Co(2000) - 
     &   At*MUE*Co(1929)*Cr(91)

        Co(2002) = 2*Cd(12701)*Co(1982)*UUStopC(2,1) - 
     &   2*Cd(12700)*Co(1983)*UUStopC(2,2) + Co(2001)*UUStopC(2,3)

        Co(2003) = 2*MTy*SA*Co(213) + AtC*Co(1986)*UStop(1,1) + 
     &   At*MUE*Co(1929)*UUStop(1,3)

        Co(200) = 2*CA*MTy*Co(197) + AtC*Co(199)*UStop(2,1) - 
     &   At*MUE*Co(138)*UUStop(3,4)

        Co(2014) = 2*MTy*SA*Cd(12982) + 
     &   AtC*Co(1984)*UStop(2,1) + At*MUE*Co(1929)*UUStop(2,3)

        Co(2035) = Cd(9792)*Cd(12701)*UCStopC(1,3)*UStopC(1,1) + 
     &   Cd(9793)*Cd(12700)*UCStop(1,3)*UStopC(1,2) - 
     &   Co(2003)*UUStopC(1,3)

        Co(2038) = 2*MTy*SA*Co(213) + AtC*Co(2037) + 
     &   At*MUE*Co(1929)*UUStop(1,3)

        Co(2039) = Cd(9792)*Cd(12701)*UCStopC(1,3)*UStopC(1,1) + 
     &   Cd(9793)*Cd(12700)*UCStop(1,3)*UStopC(1,2) - 
     &   Co(2038)*UUStopC(1,3)

        Co(2087) = UStop2(1,1)**2 - 
     &   (Cd(9948) + 2*UStop2(1,1))*UStop2(1,2)

        Co(2089) = Cd(9948)*MStop2(1) + Cr(85)*MStop2(3) - 
     &   MGlpT2*UStop2(1,2)

        Co(2104) = Cd(47)**2*Cd(9949)*MGlpTmSt2(2)*
     &    (Cd(10658) - Cd(9948)*UStop2(1,1)**2) - 
     &   Cd(50)*Cd(9951)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Co(2116) = 2*Cd(47)*Cd(50)*
     &    (Cd(50)*Cd(9947)*MGlpTmSt2(1) - 
     &      Cd(47)*Cd(9949)*Cd(9972)*MGlpTmSt2(2)) - 
     &   Cd(9980)*MStop2(3)

	Co(2127) = M_3 + At*(Cd(9948) - UStop2(1,2))

        Co(2132) = M_3*Cr(28) - 3*M_3C*UCStop(1,3)**2 - 
     &   M_3*Cd(9948)*UStop2(1,1)

        Co(2137) = -((Cd(10031)*Cd(10039))/
     &      (Cd(47)**2*Cd(50)**2)) + 
     &   (2*(2*Cd(47)*Cd(50)*Cd(10014)*Cd(10094) + 
     &        Cd(10044)*MStop2(3)))/(Cd(47)*Cd(50)*MStop2(3)**2)

        Co(2148) = Cd(50)**2*(2*Cd(47)*Ce(411) + Co(1057)) + 
     &   Cd(47)**2*Co(2147)

	Co(2155) = -(MGlmT2*Cd(50)*(Co(1047) - 2*Co(1049))) + Co(1063)

        Co(214) = 2*CA*MTy*Co(213) + AtC*Co(199)*UStop(1,1) - 
     &   At*MUE*Co(138)*UUStop(1,3)

	Co(2162) = Cd(9948)*Co(1034) - Co(1037) + Co(1034)*UStop2(1,1)

        Co(2172) = -(Cd(10658)*Co(1047)) + Cd(9948)*Co(2171) + 
     &   Co(1049)*Cr(75)**2*UStop2(1,1)

        Co(216) = Cd(693)*Cd(12701)*UCStopC(1,3)*UStopC(1,1) - 
     &   Cd(694)*Cd(12700)*UCStop(1,3)*UStopC(1,2) - 
     &   Co(214)*UUStopC(1,3)

	Co(2194) = Co(1049) + Co(1147)

        Co(2199) = Cd(9948)*Co(2198) - Co(1034)*Cr(28) + 
     &   Co(1037)*Cr(75)*UStop2(1,1)

        Co(2202) = Cd(10020)*Co(1049) + 
     &   Cd(9948)*Co(1047)*Cr(75) - Co(1047)*UStop2(1,1)**2

        Co(219) = 2*CA*MTy*Co(213) + AtC*Co(218) - 
     &   At*MUE*Co(138)*UUStop(1,3)

        Co(2213) = Cd(9948)*Co(2212) - Co(1037)*Cr(28) + 
     &   Co(1034)*Cr(75)*UStop2(1,1)

        Co(2215) = Cd(10020)*Co(1047) + 
     &   Cd(9948)*Co(1049)*Cr(75) - Co(1049)*UStop2(1,1)**2

        Co(220) = Cd(693)*Cd(12701)*UCStopC(1,3)*UStopC(1,1) - 
     &   Cd(694)*Cd(12700)*UCStop(1,3)*UStopC(1,2) - 
     &   Co(219)*UUStopC(1,3)

	Co(2236) = Co(2147) + Cd(50)*Cd(9949)*MGlpTmSt2(2)

        Co(2237) = Cd(50)**2*(Co(1057) + Cd(47)*Co(1177)) + 
     &   Cd(47)**2*Co(2236)

        Co(2240) = -(MGlmT2*Cd(50)*
     &      (Cd(10103)*Co(1047) + Cd(50)*Co(1049))) + 
     &   Cd(10103)*Co(1063)

	Co(2244) = 2*Co(1042) - Co(1056) - 2*Co(1043)*UStop2(1,2)

	Co(2248) = 2*Co(1050) - Co(1061) - 2*Co(1051)*UStop2(1,2)

        Co(2251) = -(Cd(10658)*Co(1049)) + Cd(9948)*Co(2250) + 
     &   Co(1047)*Cr(75)**2*UStop2(1,1)

        Co(2285) = Cd(50)*
     &    (4*Cd(50)*Cd(537)*Ce(411) + Cd(9959)*Cd(10103)*Co(2194))-
     &     2*Cd(10103)*Co(1051)*UStop2(1,1)**2

        Co(2286) = Cd(47)**2*Co(2285) - 
     &   2*Cd(50)**2*Cd(10103)*(Co(1043) - Cd(47)*Co(1750))*
     &    UStop2(1,1)**2

	Co(2290) = Cd(9958)*Co(1049) + Co(1694)

        Co(2291) = Cd(50)*
     &    (4*Cd(50)*Cd(10020)*Ce(411) + Cd(10103)*Co(2290)) + 
     &   2*Cd(10103)*Co(2289)*UStop2(1,2)

	Co(2294) = Co(1099) - 4*Cd(50)*Ce(411)*UStop2(1,2)

	Co(228) = 4*CA*MTy*Co(197) + AtC*Co(227) - At*MUE*Co(138)*Cr(91)

        Co(2303) = MGlmT2*Co(2302) + 
     &   2*Cd(47)*Cd(10369)*
     &    (T134(MGl2,MTy2,MStop2(1)) - T134(MGl2,MTy2,MStop2(2)))

        Co(229) = 2*Cd(12701)*Co(194)*UUStopC(2,1) - 
     &   2*Cd(12700)*Co(226)*UUStopC(2,2) - Co(228)*UUStopC(2,3)

	Co(2311) = Co(1710) - Cd(9949)*MGlpTmSt2(2)

        Co(2320) = MGlmT2*
     &    (Co(1047) - Cd(9948)*Cd(9959)*Co(1049)) + 
     &   Co(2319)*UStop2(1,1)

        Co(2325) = MGlmT2*Co(2184) + 
     &   Cd(47)*(2*T134(MGl2,MTy2,MStop2(1)) - 
     &      2*T134(MGl2,MTy2,MStop2(2)))

        Co(2326) = MGlmT2*Co(2257) + 
     &   Cd(47)*(-2*T134(MGl2,MTy2,MStop2(1)) + 
     &      2*T134(MGl2,MTy2,MStop2(2)))

	Co(2329) = Co(1710) - 2*Ce(411)*UStop2(1,2)

	Co(2330) = Co(1099) - 2*Cd(50)*Co(2329)

        Co(2331) = Cd(50)**2*(Co(1140) - 2*Cd(47)*Co(1177)) + 
     &   Cd(47)**2*Co(2330)

        Co(2349) = Cd(47)**2*Co(2348) + 
     &   2*Cd(50)**2*(4*Cd(47)*Ce(411) - Co(1042))*UStop2(1,1)**2

        Co(2354) = Cd(47)*Co(1706) + 
     &   MGlmT2*(-Co(1037) + Co(1034)*UStop2(1,2))

        Co(2361) = MGlmT2*
     &    (Cd(9948)*Cd(9959)*Co(1047) - Co(1049)) + 
     &   Co(2360)*UStop2(1,1)

	Co(2366) = Co(1071) - 4*Ce(411)*UStop2(1,1)

        Co(2392) = Co(1063) + 
     &   2*Cd(50)*(2*Cd(50)*Co(1712) + MGlmT2*Co(2391))

        Co(240) = 2*CA*MTy*Cd(12982) + AtC*Co(196)*UStop(2,1) - 
     &   At*MUE*Co(138)*UUStop(2,3)

	Co(2475) = Cd(9948)*Co(1047) - Co(1049) + Co(1047)*UStop2(1,1)

        Co(246) = 2*Cd(12701)*Co(194)*UUStopC(1,1) - 
     &   2*Cd(12700)*Co(226)*UUStopC(1,2) - Co(228)*UUStopC(1,3)

        Co(2502) = -(Cr(139)*dMf133eps(0)) + 
     &   1/2.D0*(Cd(25)*dMf133eps(-1)) - 2*dMf133eps(1)*MStop2(3)

        Co(2506) = (Co(2505) + MUE*TB*Cd(154)*UCStop(1,3) + 
     &     MUEC*TB*Cd(158)*UCStopC(1,3))/TB

        Co(2515) = (2*Cd(10014)*Co(2512) + 
     &     (-(TB*Co(2513)) + Co(2514))/(Cd(47)*Cd(50)*MStop2(3)))/
     &   TB

        Co(2517) = (Co(2516) + MUE*TB*Cd(238)*UCStop(1,3) - 
     &     MUEC*TB*Cd(240)*UCStopC(1,3))/TB

	Co(2519) = MUE*Cd(335) - (MUEC*TB*Cd(329) + Co(2518))/TB

	Co(2522) = -A0delGl + MGl2 - MGl2*Log(MGl2/Q)

        Co(2526) = 2*Cd(47)*Cd(252)*Cd(12836)*Co(2508) + 
     &   TB*(2*Cd(50)*Co(2498)*Co(2517) + Co(2519)*Co(2525))

        Co(263) = (Cd(9951)*Co(261))/Cd(47)**2 - 
     &   (At*Cd(50)*Cd(9950)*Cd(9955)*Cd(9979) + 
     &      Co(262)*MStop2(3))/(Cd(50)*Cd(9950)*MStop2(3)**2)

        Co(25) = At*MUE*Cd(545)*Cr(28) + 
     &   AtC*(MUEC*Co(24)*Cr(28) + 2*MUE*Cd(545)*UCStop(1,3)**2)

        Co(268) = At*Cd(9966) + 2*M_3*Cr(28) + 
     &   2*Cd(9967)*UCStop(1,3)**2

        Co(273) = 2*Cd(47)*Cd(50)*
     &    (Cd(50)*Cd(9947)*MGlpTmSt2(1) - 
     &      Cd(47)*Cd(9949)*MGlpTmSt2(2)) - Cd(9980)*MStop2(3)

        Co(276) = (Cd(9980)*Co(274))/(Cd(47)**2*Cd(50)**2) + 
     &   (At*Cd(47)*Cd(50)*Cd(9955)*Cd(9956) + Co(275)*MStop2(3))/
     &    (Cd(47)*Cd(50)*MStop2(3)**2)

        Co(279) = (Cd(9980)*Co(277))/(Cd(47)**2*Cd(50)**2) + 
     &   (-(AtC*Cd(47)*Cd(50)*Cd(9955)*Cd(9979)) + 
     &      Co(278)*MStop2(3))/(Cd(47)*Cd(50)*MStop2(3)**2)

        Co(283) = Cd(9948)*Co(282) - 2*M_3*UCStopC(1,3)**2 + 
     &   AtC*Cd(9958)*UStop2(1,1)

        Co(284) = -(Cd(9995)*UCStopC(1,3)**2) + 
     &   UStop2(1,1)*(Cd(9948)*Cd(9994) + AtC*UStop2(1,2)**2)

        Co(287) = (Cd(9951)*Co(281))/Cd(47)**2 + 
     &   (AtC*Cd(50)*Cd(9950)*Cd(9955)*Cd(9956) - 
     &      Co(286)*MStop2(3))/(Cd(50)*Cd(9950)*MStop2(3)**2)

        Co(288) = AtC*Cd(9966) + 2*M_3C*Cr(28) + 
     &   2*Cd(10001)*UCStopC(1,3)**2

        Co(293) = (Cd(10031)*Cd(10086) - 
     &     (2*Cd(47)*Cd(50)*
     &        (2*Cd(47)*Cd(50)*Cd(10014)*Cd(10085) + 
     &          Cd(10022)*MStop2(3)))/MStop2(3)**2)/
     &   (Cd(47)**2*Cd(50)**2)

        Co(295) = Cd(10016)*Co(294)*UCStop(1,3) - 
     &   3*M_3C*AtC*UCStop(1,3)**3 - 
     &   M_3*At*Cd(10016)*UCStopC(1,3)

        Co(298) = (Cd(10031)*Co(296) - 
     &     (2*Cd(47)*Cd(50)*
     &        (2*Cd(47)*Cd(50)*Cd(10014)*Co(295) + 
     &          Co(297)*MStop2(3)))/MStop2(3)**2)/
     &   (Cd(47)**2*Cd(50)**2)

        Co(299) = M_3C*Cr(28) - 3*M_3*UCStopC(1,3)**2 - 
     &   M_3C*Cd(9948)*UStop2(1,1)

        Co(302) = M_3C*AtC*Cd(10016)*UCStop(1,3) + 
     &   Cd(10016)*Co(294)*UCStopC(1,3) + 
     &   3*M_3*At*UCStopC(1,3)**3

        Co(305) = (Cd(10031)*Co(303) - 
     &     (2*Cd(47)*Cd(50)*
     &        (2*Cd(47)*Cd(50)*Cd(10014)*Co(302) + 
     &          Co(304)*MStop2(3)))/MStop2(3)**2)/
     &   (Cd(47)**2*Cd(50)**2)

        Co(319) = Co(316)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(2*MGlmT2*MStop2(3)*UStop2(1,1) + 
     &      Co(318)*UStop2(1,2))

        Co(322) = 2*M_3C*MStop2(3)*UStop2(1,1) - 
     &   AtC*(Cd(1281)*UStop2(1,1)**2 + MStop2(1)*UStop2(1,2)**2)

	Co(323) = Cd(47) - Cd(152)*MStop2(2) + MStop2(1)**2*UStop2(1,2)

	Co(327) = AtC*Cd(10049) - At*Cd(10219) + TB*Co(2)

        Co(332) = Cd(47)*Cd(152)*Ce(517) + Cd(9948)*Co(331) - 
     &   2*Cd(164)*Ce(633)*MStop2(3)

	Co(335) = At*Cd(10542) - AtC*Co(332) + 2*TB*Co(334)

        Co(340) = Cd(47)*
     &    (Cd(9948)*Cd(10453)*Cd(10834) + Co(339)) - 
     &   MTy2*Cd(164)*Cd(10409)*MStop2(3)

        Co(344) = Cd(47)*
     &    (Cd(9948)*Cd(10303)*Cd(10840) - Co(343)) + 
     &   MTy2*Cd(164)*Cd(10410)*MStop2(3)

        Co(358) = -(MGlmT2*Cd(13817)) + 
     &   MStop2(1)*(MGlmT2*Cd(9948) + MStop2(2))

        Co(359) = Cd(47)**2*Cd(1272) + 
     &   MStop2(1)*(2*Cd(47)*Co(358) - Cd(10280)*MStop2(3))

        Co(362) = -(Cd(9963)*UCStop(1,3)**2) + 
     &   UStop2(1,1)*(Cd(9948)*Cd(9962) + At*UStop2(1,2)**2)

        Co(365) = 4*MGlmT2*Cd(47)*Co(265) - 
     &   2*AtC*Cd(10280)*UCStop(1,3)**2 + 
     &   At*Cd(10371)*UStop2(1,1)*UStop2(1,2)

	Co(370) = Co(2126) + AtC*UStop2(1,1)*UStop2(1,2)

        Co(374) = -(Co(373)*UCStopC(1,3)**2) + 
     &   2*UStop2(1,1)*(Cd(9948)*Ce(592) + 
     &      AtC*Cd(10284)*UStop2(1,2)**2)

        Co(378) = Cd(9948)*Co(377) - 
     &   2*Cd(47)*Cd(152)*Ce(517)*MStop2(1) + 
     &   Cd(10279)*Cd(10382)*MStop2(3)

	Co(382) = AtC*Co(378) + TB*Co(380) + At*Co(381)

        Co(383) = MUE*Cd(10279)*UCStop(1,3)**2 + 
     &   MUEC*UStop2(1,1)*UStop2(1,2)*
     &    (2*Cd(10280)*UStop2(1,1) + Cd(10279)*UStop2(1,2))

        Co(385) = 2*Cd(47)*MStop2(1)*
     &    (-Co(384) + Cd(9948)*Cd(10058)*MStop2(1)) + 
     &   Co(383)*MStop2(3)

        Co(390) = Cd(9948)*Cd(10390) - 
     &   Cd(47)*Co(389)*MStop2(1) - Cd(10280)*Cd(10388)*MStop2(3)

	Co(391) = AtC*Co(385) + TB*Co(388) + 2*At*Co(390)

        Co(392) = 2*Cd(47)*
     &    (Cd(9948)*Cd(10453)*Cd(10834) + Co(339))*MStop2(1) - 
     &   Cd(10280)*Cd(10409)*MStop2(3)

        Co(393) = 2*Cd(47)*
     &    (Cd(9948)*Cd(10303)*Cd(10840) - Co(343))*MStop2(1) + 
     &   Cd(10280)*Cd(10410)*MStop2(3)

        Co(399) = -(AtC*Cd(795)) + 
     &   At*(MUEC*UCStopC(1,3)**2 + MUE*Cd(9981)*UStop2(1,1))

        Co(405) = Cd(47)**2*Ce(1031) + 
     &   MStop2(1)*(-4*Cd(47)*Cd(164)*Co(404) + 
     &      Cd(12371)*MStop2(3))

        Co(408) = 2*Cd(47)*Cd(164)*Cd(1281)*MStop2(1) - 
     &   (2*Cd(47)**2 + Cd(12371)*MStop2(1))*MStop2(3)

        Co(411) = 2*M_3*MStop2(3)*UStop2(1,1) - 
     &   At*(Cd(1281)*UStop2(1,1)**2 + MStop2(1)*UStop2(1,2)**2)

        Co(412) = 2*Cd(47)*Cd(164)*Co(411)*MStop2(1) + 
     &   At*(2*Cd(47)**2 + Cd(12371)*MStop2(1))*MStop2(3)*
     &    UStop2(1,1)**2

        Co(414) = 2*AtC*Cd(47)**2*Cd(1272) + 
     &   MStop2(1)*(4*Cd(47)*Cd(164)*Co(413) + 
     &      AtC*Cd(12371)*MStop2(3))

        Co(419) = -(Co(418)*MStop2(3)*UStop2(1,1)*UStop2(1,2)) + 
     &   2*Cd(47)*Co(416)*MStop2(1)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2)

        Co(425) = Cd(9948)*Co(424) + 
     &   2*Cd(47)*Ce(601)*Cr(28)*MStop2(1) - 
     &   Cd(12371)*MStop2(3)*UStop2(1,1)*UStop2(1,2)**2

	Co(430) = -(At*Cd(12425)) + AtC*Ce(550) - TB*Co(429)

        Co(442) = AtC*Cd(11824) - 2*Cd(164)*Co(421) + 
     &   2*At*Cd(47)*UCStopC(1,3)**2

        Co(448) = Cd(47)*Cd(164)*
     &    (Cd(9948)*Cd(10058) + Ce(517)) + 
     &   2*Cd(152)*Cd(795)*MStop2(1)*MStop2(3)

        Co(449) = Cd(47)*Cd(164)*Cd(10057) + 
     &   2*Cd(152)*MStop2(1)*MStop2(3)*
     &    (Cd(10388) + MUE*Cd(9948)*UStop2(1,1)**2)

	Co(451) = AtC*Co(448) - At*Co(449) + 2*TB*Co(450)

        Co(460) = -(At*Cd(10219)) + 
     &   AtC*(MUE*UCStop(1,3)**2 + MUEC*Cd(9981)*UStop2(1,1))

        Co(468) = At*Cd(11941) + 2*Cd(251)*Co(467) + 
     &   2*AtC*Cd(50)*UCStop(1,3)**2

        Co(469) = Cd(778)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(652)*UStop2(1,1) + MStop2(1)*UStop2(1,2))

        Co(46) = AtC*(TB*Cd(12071) + 2*AtC*UUStop(1,1)) - 
     &   At*MUE*TB*UUStop(1,2)

        Co(478) = Cd(50)*Cd(251)*Cd(10674) + 
     &   MUEC*Cd(9948)*Co(477) - 
     &   2*Cd(179)*Cd(10382)*MStop2(2)*MStop2(3)

	Co(480) = At*Ce(673) + AtC*Co(478) - 2*TB*Co(2)*Co(479)

        Co(487) = Cd(50)*Co(484) + Cd(9948)*Co(486) + 
     &   2*MGl2*Cd(179)*Cd(10382)*MStop2(3)

        Co(489) = Cd(50)*(MGlmT2*Cd(427)*Cd(9948) + Co(488)) - 
     &   MGl2*Cd(179)*Co(2)*MStop2(3)

        Co(492) = Cd(50)*(Cd(9948)*Co(490) + Co(491)) + 
     &   2*MGl2*Cd(179)*Cd(10219)*MStop2(3)

	Co(493) = AtC*Co(487) - 2*TB*Co(489) - At*Co(492)

        Co(496) = -2*Cd(50)*
     &    (Cd(50)*Cd(843) + 2*Cd(251)*Cd(9972)*MStop2(1)*MStop2(2))
     &     + Cd(12483)*MStop2(2)*MStop2(3)

        Co(498) = 2*Cd(50)**2*Cd(577) + 
     &   MStop2(2)*(2*Cd(50)*Cd(251)*Cd(12506) + 
     &      Cd(12483)*MStop2(3)*UStop2(1,1)*UStop2(1,2))

        Co(503) = -2*At*Cd(50)**2*Cd(843) + 
     &   MStop2(2)*(4*Cd(50)*Cd(251)*Cd(10819) + 
     &      At*Cd(12483)*MStop2(3))

        Co(507) = 4*Cd(50)*(MGlmT2*Cd(9948) + Ce(797))*
     &    MStop2(2) - Cd(12483)*MStop2(3)

        Co(509) = 2*Cd(50)*(Cd(9948)*Ce(675) + Co(508))*
     &    MStop2(2) + Cd(12483)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Co(516) = Cd(9948)*Co(515) + 
     &   2*Cd(50)*Co(513)*MStop2(2) + 
     &   Cd(12483)*MStop2(3)*UStop2(1,1)*UStop2(1,2)**2

        Co(520) = -2*Cd(50)*Cd(251)*Ce(672)*MStop2(2) - 
     &   Cd(10388)*Cd(12483)*MStop2(3) + 
     &   MUE*Cd(9948)*Ce(716)*UStop2(1,1)**2

	Co(521) = AtC*Cd(12539) - TB*Co(2)*Co(519) + At*Co(520)

        Co(525) = -4*Cd(50)*Cd(251)*
     &    (Cd(10764) - M_3C*Cd(9948)*Cr(75))*MStop2(2) - 
     &   Cd(10039)*Cd(12483)*MStop2(3)

        Co(537) = At - 2*M_3C*UCStop(1,3)**2 + 
     &   2*M_3*UStop2(1,1)*UStop2(1,2)

        Co(548) = Cd(50)**2*Cd(1272) + 
     &   MStop2(2)*(-2*Cd(50)*Co(547) + Cd(10604)*MStop2(3))

        Co(550) = 2*Cd(50)*Cd(179)*
     &    (MUEC*Cd(9948)*Cd(9983) + Cd(10674))*MStop2(2) + 
     &   Co(549)*MStop2(3)

        Co(552) = 2*Cd(50)*Cd(179)*Ce(672)*MStop2(2) + 
     &   Cd(10388)*Cd(10603)*MStop2(3) - 
     &   MUE*Cd(9948)*Co(551)*UStop2(1,1)**2

	Co(555) = AtC*Co(550) - At*Co(552) - TB*Co(554)

        Co(568) = Cd(9948)*Co(567) - Cd(50)*Co(565)*MStop2(2) - 
     &   Cd(10382)*Cd(10604)*MStop2(3)

	Co(569) = -MStop2(2) + MGlmT2*(1 + UStop2(1,2))

        Co(572) = MUEC*Cd(10603)*UCStopC(1,3)**2 + 
     &   MUE*UStop2(1,1)*UStop2(1,2)*
     &    (2*Cd(10604)*UStop2(1,1) + Cd(10603)*UStop2(1,2))

        Co(575) = 2*Cd(50)*(Cd(9948)*Co(573) - Co(574))*
     &    MStop2(2) + Co(572)*MStop2(3)

        Co(578) = 2*Cd(50)*
     &    (Co(577) + M_3*MGlmT2*Cd(9948)*Cr(75))*MStop2(2) - 
     &   Cd(10086)*Cd(10604)*MStop2(3)

        Co(581) = 2*Cd(50)*(Co(580) + MGlmT2*Cd(9948)*Cr(75))*
     &    MStop2(2) - Cd(10604)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Co(588) = Cd(50)*Cd(9966) - MGlmT2*Cd(10803) + 
     &   Cd(10802)*MStop2(2)

        Co(58) = Co(13)/TB + 
     &   (2*(TB*Co(57) + Co(4)*MStop2(4)))/(TB*MStop2(3))

        Co(596) = -(AtC*Cd(10806)) - 
     &   2*At*Cd(10807)*UCStopC(1,3)**2 + 
     &   2*M_3C*MStop2(2)*MStop2(3)*UStop2(1,1)

	Co(602) = AtC*Cd(10804) + Cd(9948)*Co(596) - 2*Co(601)

	Co(613) = MGl2 - MStop2(2) + 2*M_3C*MTy*UCStop(1,3)

	Co(618) = MGl2*Co(565) + MGlmT2*Co(616) - Co(617)*MStop2(2)

        Co(620) = MGl2*Co(566) - MGlmT2*Co(619) + 
     &   MTy2*MUEC*Cd(9983)*MStop2(2)

	Co(626) = -(MGl2*Co(574)) + MGlmT2*Co(625) - Co(624)*MStop2(2)

        Co(628) = MGl2*Co(573) + Co(627)*MStop2(2) - 
     &   MGlmT2*MTy2*MUE*UStop2(1,1)**2

        Co(646) = ((Co(645)*MGlpTmSt2(1) - 
     &        TB*Co(642)*UCStop(1,3)*UStopC(1,2))/Cd(47)**2 - 
     &     (Co(644)*MGlpTmSt2(2)*UStopC(2,2))/Cd(50)**2)/TB

        Co(648) = -(At*MUE*UUStop(1,1)) + 
     &   MUEC*(Co(647) - MUEC*TB*UUStop(1,2))

        Co(662) = 2*Cd(10925)*Cd(12701)*UCStop(3,4) - 
     &   Co(661)*UStopC(2,2) + 2*MUEC*TB*UStop2(1,1)*UUStop(1,2)

        Co(669) = -(At*MUE*Cd(11161)*UStop(1,1)) + 
     &   MUEC*(Co(668) - 4*MUEC*TB*UStop(2,2)*UUStop(1,2))

        Co(675) = -(At*MUE*Cd(11019)*UStop(2,1)) + 
     &   MUEC*(Co(674) - 4*MUEC*TB*UStop(1,2)*UUStop(2,2))

        Co(676) = 2*MUEC*Cd(11106)*Cd(12701)*UCStopC(1,3)*
     &    UStopC(2,1) + Co(673)*UUStopC(2,2) + Co(675)*UUStopC(2,3)

        Co(683) = -2*Cd(47)*Co(641)*UStop2(1,1) + 
     &   MStop2(3)*(4*Cd(9948)*Co(679) + Co(682)*UStop2(1,2))

        Co(691) = MUE*Cd(12700)*UCStop(1,3)*UCStop(3,3) - 
     &   3*MUEC*Cd(12701)*UCStop(3,4)*UCStopC(1,3) + 
     &   Co(641)*Cr(91)*UUStopC(2,3)

        Co(692) = (MTy*MUE*Cd(50)*Cd(12700)*MGlpTmSt2(1)*
     &      MStop2(3)*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) - 
     &     Cd(47)*(2*MTy*Cd(47)*Co(691)*UCStopC(3,3) + 
     &        Cd(50)*Co(690)*UUStopC(1,3)))/
     &   (Cd(47)**2*Cd(50)*MStop2(3))

        Co(68) = Co(67)*MStop2(1)**2 + 
     &   2*MStop2(2)*(-2*Cr(85)*MStop2(1) + MStop2(2)*UStop2(1,1))

        Co(702) = -((At*MTy*Cd(50)*MGlpTmSt2(1)*MStop2(3)*
     &        UCStopC(1,3)*UStop(1,2)*UStop2(1,1) - 
     &       2*MTy*Cd(47)**2*Co(699)*UStop2(1,2) + 
     &       Cd(47)*Cd(50)*Co(701)*UStopC(1,1))/
     &     (Cd(47)**2*Cd(50)*MStop2(3)))

        Co(723) = Cd(11698)*MStop2(3) - 
     &   Cd(12114)*MStop2(1)*UStop2(1,1) - 
     &   MGlpT2*Cr(56)*UStop2(1,2)

	Co(729) = Co(2083) - MStop2(1)*UStop2(1,1)**2

        Co(750) = AtC*Co(748) - 2*MUE**2*TB*UUStop(1,1) + 
     &   Co(749)*UUStopC(1,1)

        Co(76) = -2*Cd(172)*MStop2(2)**2 + 
     &   2*MStop2(1)*(2*MStop2(2)*UStop2(1,1) + 
     &      MStop2(1)*UStop2(1,2))

	Co(782) = Co(822) - Co(772)*MStop2(3)

	Co(783) = Co(823) - At*Cd(10950)*MStop2(3)

        Co(784) = -2*MUE*Co(783)*MStop2(1)*UStop2(1,1)*
     &    UStopC(2,2) + MUEC*
     &    (Co(782)*UStopC(1,1) + MUEC*TB*Cd(10971)*UUStopC(1,1))

        Co(803) = -2*MUE*Co(768)*MStop2(1)*UStop2(1,1)*
     &    UStopC(2,2) + MUEC*
     &    (Co(781)*UStopC(1,1) - MUEC*TB*Cd(10964)*UUStopC(1,1))

        Co(811) = UStop2(1,1)*(UStop2(1,1) - 2*UStop2(1,2)) + 
     &   6*UStop2(1,2)**2

        Co(80) = TB*Co(79)*UStop2(1,2) + 
     &   UCStop(1,3)*(AtC*MUE*TB*Cd(651)*UCStop(1,3) + 
     &      4*MTy*Cd(12700)*MStop2(3)*UStop2(1,2))

        Co(824) = -2*MUE*Co(823)*MStop2(1)*UStop2(1,1)*
     &    UStopC(2,2) + MUEC*
     &    (Co(822)*UStopC(1,1) + MUEC*TB*Cd(11267)*UUStopC(1,1))

	Co(828) = AtC*Co(826) - Co(810)*MStop2(1) - Co(827)*MStop2(3)

        Co(840) = 8*Cd(9948)*Co(839)*MStop2(2)*MStop2(3) - 
     &   4*Cd(11175)*MStop2(2)**2*UStop2(1,1)**2 + 
     &   2*Co(838)*MStop2(3)**2*UUStop(2,1)

        Co(841) = -(Cd(970)*MStop2(3)**2) + 
     &   MStop2(2)*UStop2(1,1)*
     &    (MStop2(2)*UStop2(1,1) + 2*MStop2(3)*UStop2(1,2))

        Co(844) = -(MUEC*MStop2(1)**2*UCStopC(1,3)**2) + 
     &   MUE*(Co(843)*MStop2(2)*UStop2(1,1) + 
     &      2*MStop2(3)**2*UStop2(1,2)**2)

        Co(845) = AtC*Co(840) + 
     &   At*(8*MUE*MStop2(2)*MStop2(3)*UCStopC(3,1)**2*
     &       UUStop(1,2) - 4*Co(844)*UUStop(2,2))

        Co(850) = Co(847)*UCStop(3,1)*UStopC(2,1) - 
     &   Co(849)*UUStopC(2,1) - 
     &   MUE*Cd(12700)*Cr(91)**2*UUStopC(2,2)

        Co(870) = Co(864)*MStop2(3) + AtC*Co(867)*UStop(2,1) - 
     &   Co(869)*MStop2(2)*UUStopC(3,1)

        Co(87) = 4*MStop2(1)*UStop2(1,1)**3 - 
     &   Cd(594)*UStop2(1,1)**4 + Co(86)*UStop2(1,2)

	Co(890) = MUE*Cd(970) + Co(889) + 2*MUEC*UCStopC(1,3)**2

        Co(894) = MUEC*MStop2(1)*UCStopC(1,3)**2 + 
     &   MUE*(Cd(9948)*MStop2(2)*UStop2(1,1) + 
     &      2*MStop2(3)*UStop2(1,2)**2)

	Co(898) = Co(851) + 2*Co(897)

	Co(905) = -Co(865) + 2*Co(904)

	Co(909) = Co(864) + 2*Co(908)

	Co(919) = Co(881) + 2*Co(918)

	Co(929) = Co(874) + 3*Co(916)

	Co(930) = Co(920) + Co(929)

        Co(946) = At*UUStop(1,1)*UUStop(2,2)*UUStopC(2,1) + 
     &   2*MUE*TB*UUStop(3,1)*UUStop(3,2)*UUStopC(2,2) - 
     &   Co(945)*UUStopC(2,3)

        Co(948) = 2*MUEC*Cd(12701)*UUStop(3,1)*UUStop(3,2)*
     &    UUStopC(2,1) - At*MUE*UUStop(1,2)*UUStop(2,1)*
     &    UUStopC(2,2) + Co(947)*UUStopC(2,3)

        Co(953) = Co(949)*UUStopC(2,1) - 
     &   MUE*Co(950)*UUStopC(2,2) - 2*Co(952)*UUStopC(2,3)

        Co(956) = 2*Cd(9948)*Cd(11215)*UCStop(3,3) + 
     &   MUEC*UCStop(3,1)*UStop2(1,2)*UUStop(3,3) - 
     &   Co(955)*UUStopC(2,2)

        Co(958) = MUEC*UUStop(1,1)*UUStop(2,2)*UUStopC(2,1) - 
     &   2*MUE*UUStop(3,1)*UUStop(3,2)*UUStopC(2,2) - 
     &   Co(957)*UUStopC(2,3)

        Co(959) = MUEC*Co(954)*UStopC(2,2)*UUStopC(1,1) + 
     &   Co(958)*UUStopC(1,2) - Co(956)*UUStopC(1,3)

        Co(980) = Co(974)*UUStopC(2,1) - 
     &   MUE*Co(977)*UUStopC(2,2) + 2*Co(979)*UUStopC(2,3)

        Co(982) = 2*MUEC*Cd(12701)*Cr(91)*UCStopC(1,3)*
     &    UStopC(2,1) - MUE*Co(788)*UCStop(1,3)*UStopC(2,2) - 
     &   Co(981)*UUStopC(2,3)

        Co(988) = 4*Cd(11177)*MStop2(2)*MStop2(3)*UCStop(1,3) + 
     &   MUEC*Cd(9948)*Co(987)*UStop(1,1) + Co(986)*UStopC(2,2)

        Co(993) = Co(989)*UUStopC(2,1) + 
     &   2*MUE*Cd(11657)*UUStopC(2,2) - Co(992)*UUStopC(2,3)

        Co(1) = AtC*(-(MUEC*Cr(28)) + 2*MUE*UCStop(1,3)**2) + 
     &   At*(MUE*Cr(28) - 2*MUEC*UCStopC(1,3)**2)

        Co(1003) = -4*Co(1002)*MStop2(3) + 
     &   2*Cd(9948)*Co(996)*UCStop(3,3) + 
     &   MUEC*Co(997)*UCStop(3,1)*UStop2(1,2) + 
     &   Co(1001)*UUStopC(2,2)

        Co(104) = At*(-(MUEC*UCStopC(1,3)**2) + 
     &      MUE*UStop2(1,1)*UStop2(1,2)) + 
     &   AtC*(MUE*UCStop(1,3)**2 - MUEC*UStop2(1,1)*UStop2(1,2))

	Co(1087) = 4*Ce(411) - Co(1750)

	Co(1088) = Co(1130) - Co(1147)

	Co(9) = (MUE*TB*Cd(9957) - MUEC*TB*Cd(9993) + Co(4) + Co(8))/TB

        Co(1129) = -(Cd(50)**2*
     &      (Cd(9951) - 2*Co(1040)*MGlpTmSt2(1))) + 
     &   Cd(47)**2*(Cd(9952) - 2*Co(1040)*MGlpTmSt2(2))

        Co(1238) = AtC*Co(1237)*UStop(1,1) + 
     &   MTy*(2*MUE*TB*UUStop(3,1) - 2*Cd(12701)*UUStop(3,2)) + 
     &   At*MUE*TB*UUStop(3,3)

        Co(1263) = 3 + UStop2(1,2)**2 + 
     &   UStop2(1,1)*(-5*UStop2(1,1) + 2*UStop2(1,2))

	Co(1265) = 3 - Co(1363)

        Co(1289) = AtC*UUStop(1,1)*
     &    (Co(1287)*MStop2(3)**2 + AtC*MStop2(2)**2*UUStop(1,1)) - 
     &   Co(1288)*UUStop(1,2) + 2*MTy*Co(1286)*UUStop(1,3)

        Co(11) = AtC*(MUEC*Cd(485) - 
     &      MUE*MStop2(4)*UCStop(1,3)**2) + 
     &   At*(-(MUE*Cd(485)) + MUEC*MStop2(4)*UCStopC(1,3)**2)

        Co(1345) = 2*MUEC*TB*Co(1281) + 2*Co(1279)*UStopC(1,2) + 
     &   At*(7*MUE*TB*UUStop(1,1) + 2*MTy*UUStop(1,3))

        Co(1351) = AtC*UUStop(1,1)*
     &    (Co(1287)*MStop2(3) + AtC*MStop2(2)*UUStop(1,1)) - 
     &   Co(1350)*UUStop(1,2) - 2*MTy*Co(649)*MStop2(3)*UUStop(1,3)

        Co(1359) = 11 + 2*UStop2(1,2)**2 + 
     &   UStop2(1,1)*(-10*UStop2(1,1) + 4*UStop2(1,2))

        Co(1394) = AtC*Co(1355) + 
     &   MTy*(2*MUE*TB*UUStop(1,1) - 2*Cd(12701)*UUStop(1,2)) + 
     &   At*MUE*TB*UUStop(1,3)

        Co(1445) = Co(1443)*UStop2(1,1)*UStopC(2,2) + 
     &   MTy*(2*MUE*TB*UUStop(3,1) - 2*Co(1444)*UUStop(3,2)) + 
     &   At*MUE*TB*UUStop(3,4)

        Co(1456) = 2*MUE*TB*UUStop(1,1)*UUStop(2,1) + 
     &   At*(-(UUStop(1,2)*UUStop(2,1)) + 
     &      UUStop(1,1)*UUStop(2,2) + 2*UUStop(3,1)*UUStop(3,2))

        Co(1475) = 2*MUE*TB*UUStop(1,1)*UUStop(2,1) + 
     &   At*(UUStop(1,2)*UUStop(2,1) - UUStop(1,1)*UUStop(2,2) + 
     &      2*UUStop(3,1)*UUStop(3,2))

	Co(1479) = 2*Co(1431) + Co(1511)

        Co(1523) = Co(1443)*UStop2(1,1)*UStopC(2,2) + 
     &   MTy*(6*MUE*TB*UUStop(3,1) - 6*Co(1495)*UUStop(3,2)) + 
     &   3*At*MUE*TB*UUStop(3,4)

        Co(1567) = 2*Cd(11162)*Cd(12701) + 
     &   At*Co(1566)*UCStop(3,1) - 
     &   2*(MTy*Co(951)*UCStop(3,3) + 
     &      MUEC*TB*Cd(10976)*UStop(1,2)*UStop2(1,1))

        Co(1644) = 4*Co(1643)*MStop2(3) - 
     &   2*Cd(9948)*Co(1634)*UCStop(3,3) + 
     &   Co(1632)*UCStop(3,1)*UStop2(1,2) + 2*Co(1640)*UUStopC(2,2)

        Co(1674) = Co(1571)*Co(1575) + 2*Co(1673)*MStop2(3) + 
     &   UStop2(1,1)*(-(Co(1672)*UStop(1,1)) + 
     &      MStop2(2)*UStop2(1,1)*UUStop(1,1))

        Co(1700) = -(Cd(47)*Cd(50)**2*Co(1087)) + 
     &   Cd(50)**2*Co(1140) + 
     &   Cd(47)**2*(Cd(50)*Cd(9948)*Co(1130) + Co(1699))

	Co(1854) = Co(1754) - Cd(9947)*MGlpTmSt2(1)

        Co(1864) = Cd(47)*(4*Ce(411) - Co(1049)) + 
     &   Cd(50)*(Co(1037) - 2*Cd(9947)*MGlpTmSt2(1))

        Co(204) = 2*CA*MTy - AtC*CA*UCStop(1,3) - 
     &   MUE*SA*UCStop(1,3) - Cd(6)*UCStopC(1,3)

        Co(2177) = Cd(9948)*
     &    (Co(1061) - 2*Co(1050)*UStop2(1,1)) + 
     &   2*UStop2(1,2)*(Co(1051) - Co(1050)*UStop2(1,2))

        Co(2182) = Cd(9948)*
     &    (Co(1056) - 2*Co(1042)*UStop2(1,1)) + 
     &   2*UStop2(1,2)*(Co(1043) - Co(1042)*UStop2(1,2))

        Co(2206) = (Co(1034) + Co(1037) - 
     &      2*Cd(9947)*MGlpTmSt2(1))*UStop2(1,1) + 
     &   4*Ce(411)*UStop2(1,2)

        Co(2295) = Cd(50)*
     &    (4*Cd(47)**2*Cd(9948)*Ce(411) + 
     &      Cd(50)*(-(Cd(47)*Co(1087)) + Co(1140))) + 
     &   Cd(47)**2*Co(2294)

        Co(2376) = -(Cd(10037)*Co(1047)) + 
     &   Co(1049)*(1 - 2*UStop2(1,1) - 2*UStop2(1,2))

	Co(2427) = Co(1750) - 4*Ce(411)*UStop2(1,1)

	Co(2431) = Co(2194) - 4*Ce(411)*UStop2(1,1)

	Co(2446) = -2*Cd(9972)*Ce(411) + Co(2311)

        Co(260) = (2*Cd(50)**2*Cd(9947)*MGlpTmSt2(1)*MStop2(3) + 
     &     Cd(47)*(-(Cd(9952)*MStop2(2)*MStop2(3)) + 
     &        Cd(50)**2*(2*Cd(10014) - Co(259)*MStop2(1)*MStop2(3))
     &        ))/(Cd(47)*Cd(50)**2*MStop2(3))

        Co(280) = (MUE*Co(276) - MUEC*Co(279))/TB + 
     &   (MUE**2*Co(273)*UCStop(1,3)**2 + 
     &      MUEC**2*Cd(10083)*UCStopC(1,3)**2)/
     &    (Cd(47)**2*Cd(50)**2*MStop2(3))

        Co(28) = AtC*(MUEC*Cd(539) + 
     &      2*MUE*Cd(545)*UCStop(1,3)**2) + 
     &   At*(MUE*Cd(9787) + 2*MUEC*Cd(542)*UCStopC(1,3)**2)

        Co(308) = At*(MUEC*Cd(10218)*UCStopC(1,3)**2 + 
     &      MUE*Cd(10217)*UStop2(1,1)) + 
     &   AtC*(-(MUE*Cd(10218)*UCStop(1,3)**2) + 
     &      MUEC*Cd(10231)*UStop2(1,1))

        Co(309) = -(M_3*AtC*MUEC) + 
     &   At*(M_3*MUE + M_3C*MUEC) + M_3C*MUEC**2*TB

        Co(313) = -(Cd(47)*Cd(9972)) + Cd(9948)*MStop2(1)**2 + 
     &   Cd(152)*MStop2(2)*UStop2(1,1) + Co(312)*UStop2(1,2)

	Co(320) = Co(411) - 2*AtC*MGlmT2*UCStop(1,3)**2

        Co(321) = At*Co(319) - Cd(9948)*Co(320)*MStop2(1) + 
     &   2*Ce(622)*UCStop(1,3)**2 + 
     &   2*M_3*MGlmT2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Co(31) = AtC*(-(MUEC*Cd(485)) + 
     &      MUE*MStop2(4)*UCStop(1,3)**2) + 
     &   At*(MUE*Cd(485) - MUEC*MStop2(4)*UCStopC(1,3)**2)

        Co(330) = (-(MUEC*Co(329)) + 
     &     MUE*(Co(328) + 2*MUE*TB*Cd(10445)*UCStop(1,3)**2) + 
     &     2*MUEC**2*TB*Cd(10534)*UCStopC(1,3)**2)/TB

        Co(346) = MUE**2*UCStop(1,3)*
     &    (M_3C*UCStop(1,3)**2 + M_3*UStop2(1,1)*UStop2(1,2))-
     &     MUEC**2*UCStopC(1,3)*
     &    (M_3*UCStopC(1,3)**2 + M_3C*UStop2(1,1)*UStop2(1,2))

        Co(33) = At*MUE*TB*Cd(9777) + 
     &   2*MTy*Cd(9780)*Cd(12700)*UCStop(1,3)*UStop2(1,1) + 
     &   AtC*TB*(MUEC*Cd(507) + 
     &      2*MUE*MStop2(4)*UCStop(1,3)**2*UStop2(1,1)*UStop2(1,2))

        Co(422) = AtC*Co(419) - 
     &   4*Cd(47)*Cd(164)*Co(421)*MStop2(1) + 
     &   At*(2*Cd(47)**2*Cd(1272) + Cd(12371)*MStop2(1)*MStop2(3))*
     &    UCStopC(1,3)**2

        Co(434) = (MUE*Co(432) + MUEC*Co(433))/TB + 
     &   MUE**2*Cd(12404)*UCStop(1,3) + 
     &   MUEC**2*Ce(1117)*UCStopC(1,3)

        Co(441) = -(Cd(164)*Cd(9948)*Co(411)) + 
     &   2*Co(440)*UCStop(1,3)**2 + 
     &   At*(Co(407)*UStop2(1,1)**2 + Co(439)*UStop2(1,2)**2)

        Co(446) = -(AtC*Cd(795)) + At*Cd(10050) - 
     &   MUE**2*TB*UCStop(1,3)**2 + MUEC**2*TB*UCStopC(1,3)**2

        Co(462) = MUE**2*Cd(10670)*UCStop(1,3) - 
     &   MUEC**2*UCStopC(1,3)*
     &    (UCStopC(1,3)*(MGlpT2 + 2*M_3*MTy*UCStopC(1,3)) + 
     &      2*M_3C*MTy*UStop2(1,1)*UStop2(1,2))

        Co(474) = -2*At*MStop2(1)*UCStopC(1,3)**2 - 
     &   AtC*Cd(774)*UStop2(1,1)**2 + 
     &   MStop2(3)*(2*M_3C*UStop2(1,1) + AtC*UStop2(1,2)**2)

        Co(5) = AtC*(-(MUEC*Cr(28)) - 2*MUE*UCStop(1,3)**2) + 
     &   At*(MUE*Cr(28) + 2*MUEC*UCStopC(1,3)**2)

        Co(499) = At*Co(498) + 
     &   4*Cd(50)*Cd(251)*Co(467)*MStop2(2) + 
     &   AtC*(2*Cd(50)**2*Cd(1272) - 
     &      Cd(12483)*MStop2(2)*MStop2(3))*UCStop(1,3)**2

	Co(51) = MUE*Cd(9957) - MUEC*Cd(9993) + (Co(4) - Co(8))/TB

        Co(533) = AtC*(MUE*Cd(10259)*UCStop(1,3)**2 + 
     &      MUEC*Cd(10258)*UStop2(1,1)) + 
     &   At*(-(MUEC*Cd(10259)*UCStopC(1,3)**2) + 
     &      MUE*Cd(10269)*UStop2(1,1))

        Co(543) = 2*AtC*UCStop(1,3)**2 + 
     &   2*UStop2(1,1)*(M_3 + 2*At*UStop2(1,1)) + 
     &   At*UStop2(1,2)**2

        Co(53) = AtC*(-(MUEC*Cd(590)) + 
     &      MUE*MStop2(4)*UCStop(1,3)**2) + 
     &   At*(MUE*Cd(590) - MUEC*MStop2(4)*UCStopC(1,3)**2)

        Co(563) = (-2*Cd(50)*Co(560) + AtC*Co(559)*MStop2(3) + 
     &     4*MGlmT2*Cd(50)*Co(562)*MStop2(3) - 
     &     2*At*Cd(10604)*MStop2(3)*UCStopC(1,3)**2)/MStop2(3)

        Co(564) = Co(556)/MStop2(2)**2 - 
     &   (Co(563)*MStop2(3) + 
     &      Cd(50)*(Co(557)*MStop2(1)*MStop2(2) - 
     &         2*Co(283)*MStop2(2)*MStop2(3)))/
     &    (Cd(50)**2*MStop2(2)*MStop2(3))

        Co(585) = Cd(50) + Cd(9948)*Cd(10807) - 
     &   Cd(179)*MStop2(1)*UStop2(1,1) - MGlmT2*Cd(256)*UStop2(1,2)

        Co(586) = Cd(50) - Cd(9948)*Cd(10807) - 
     &   Cd(179)*MStop2(1)*UStop2(1,1) - 
     &   Ce(1461)*MStop2(2)*UStop2(1,2)

        Co(589) = Cd(50) + 
     &   MStop2(1)*(-MGlmT2 + MStop2(2)*UStop2(1,1)) + 
     &   MGlmT2*MStop2(3)*UStop2(1,2)

	Co(653) = MGlpT2**2 + MStop2(2)**2

        Co(686) = MGlmT2**2 + MTy*Cd(12836)*MGlpTmSt2(2) - 
     &   2*MGlpT2*MStop2(2) + MStop2(2)**2

        Co(705) = MGlmT2**2 - 2*MGlpT2*MStop2(2) + 
     &   MStop2(2)**2 + MTy*Cd(12701)*MGlpTmSt2(2)*UCStopC(1,3)

	Co(735) = Co(654) + 4*MTy2*UStop2(1,1) + 4*MStop2(2)*UStop2(1,1)

        Co(758) = MUE*Co(756) + 
     &   MUEC*(-Co(757) + 4*MUEC*TB*Cr(7)*MStop2(1)**2*
     &       UCStopC(1,3)*UUStop(1,2)) - 
     &   5*Co(641)*MStop2(3)**2*UUStop(1,3)

        Co(89) = AtC*(-(MUEC*Cd(606)) - 
     &      2*MUE*Cd(609)*UCStop(1,3)**2) + 
     &   At*(MUE*Co(87) + 2*MUEC*Cd(611)*UCStopC(1,3)**2)

        Co(91) = AtC*(-(MUEC*Cd(616)) + 
     &      2*MUE*Cd(618)*UCStop(1,3)**2) + 
     &   At*(MUE*Cd(616) - 2*MUEC*Cd(619)*UCStopC(1,3)**2)

        Co(960) = AtC*Co(959) + Co(948)*UUStopC(1,1) - 
     &   MUE*Co(946)*UUStopC(1,2) - Co(953)*UUStopC(1,3)

        Co(1235) = 2*Cd(11162)*Cd(12701) + 
     &   Co(1231)*UCStop(1,3) - 
     &   2*MUEC*TB*Cd(10976)*UStop(1,2)*UStop2(1,1) + 
     &   At*(Co(1234) - 2*MTy*UStop2(1,1)*UUStop(1,3))

        Co(1425) = AtC*MUEC*TB*Co(1424) + 
     &   MTy*(-2*Co(1423)*UUStop(3,1) + 2*Cd(12701)*UUStop(3,2)) - 
     &   At*MUE*TB*UUStop(3,3) + Co(1422)*UUStopC(2,1)

        Co(1429) = 2*MTy*MUE*TB*UStop(1,1)*UUStop(2,1) + 
     &   AtC*Co(1237)*UUStop(3,1) + 
     &   At*(2*At*UStop(1,2)*UUStop(2,2) + 
     &      MUE*TB*UStop(2,2)*UUStop(3,1)) + 2*Co(1428)*UUStop(3,2)

        Co(1457) = -2*MUEC*TB*Cd(11188) + 2*MTy*Co(1456) + 
     &   At*(Co(1455)*UUStop(1,3) + 
     &      MUE*TB*UUStop(1,1)*UUStop(2,3)) - 
     &   AtC*Co(1454)*UUStop(3,1)

        Co(1471) = -(AtC*MUEC*TB*Co(1470)) + 
     &   MTy*(6*Co(1469)*UUStop(3,1) - 6*Cd(12701)*UUStop(3,2)) + 
     &   3*At*MUE*TB*UUStop(3,3) - Co(1422)*UUStopC(2,1)

        Co(1476) = 2*MUEC*TB*Cd(10996) + 2*MTy*Co(1475) + 
     &   At*(MUE*TB*UUStop(1,3)*UUStop(2,1) + 
     &      Co(1474)*UUStop(2,3)) - AtC*Co(1454)*UUStop(3,1)

	Co(1496) = 2*Co(1445) + Co(1523)

        Co(1577) = 4*Co(1571)*Co(1575)*MStop2(3) + 
     &   4*Cd(9948)*Co(1570)*UCStop(3,3) - 
     &   4*Co(1574)*UStop(1,1)*UStop2(1,1) + 
     &   4*Co(1576)*MStop2(2)*UUStop(1,1) + Cd(11657)*UUStopC(2,2)

        Co(1721) = 4*Cd(9983)*Ce(411) + 
     &   (Co(1047) - 3*Co(1049) - 2*Cd(9949)*MGlpTmSt2(2))*
     &    UStop2(1,1)**2 - Co(1049)*UStop2(1,2)**2

        Co(1853) = Co(1047) - 2*Cd(9949)*MGlpTmSt2(2) + 
     &   (-4*Ce(411) + Co(1049))*UStop2(1,1)

        Co(2113) = (Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952))*
     &    MStop2(3) + 2*Cd(47)*Cd(50)*
     &    (Cd(50)*Cd(9947)*Cd(10077)*MGlpTmSt2(1) + 
     &      Cd(47)*(Cd(9948)*Cd(9949)*MGlpTmSt2(2) + 
     &         Cd(9949)*MGlpTmSt2(2)*UStop2(1,1)))

        Co(447) = (-2*MGl2*Cd(152)*Co(446)*MStop2(3) + 
     &     Cd(47)*(MUE*Co(444) + 
     &        2*MUE**2*TB*Cd(12472)*UCStop(1,3)**2 - 
     &        MUEC*(Co(445) + 2*MUEC*TB*Ce(1103)*UCStopC(1,3)**2)))
     &    /(TB*MStop2(3))

        Co(593) = -(At*MStop2(1)*UStop2(1,1)**2) + 
     &   MStop2(3)*(2*AtC*UCStop(1,3)**2 + 
     &      2*UStop2(1,1)*(M_3 + 2*At*UStop2(1,1)) + 
     &      At*UStop2(1,2)**2)

        Co(724) = Cr(56)*MStop2(1) - 2*MGlpT2*UStop2(1,1) + 
     &   MStop2(3)*(UStop2(1,1) - UStop2(1,2))

        Co(1077) = Cd(50)*Cd(118)*(-2*Ce(411) + Co(1037)) + 
     &   Cd(47)*Cd(138)*(2*Ce(411) - Co(1049)) + 
     &   Cd(47)*Cd(50)*(-T134(MGl2,MTy2,MStop2(1)) + 
     &      T134(MGl2,MTy2,MStop2(2)))

        Co(1113) = Cd(50)**2*
     &    (Cd(9951) - Cd(9953) - 2*Co(1040)*MGlpTmSt2(1)) + 
     &   Cd(47)**2*(Cd(9952) - Cd(9954) - 2*Co(1040)*MGlpTmSt2(2))

        Co(1135) = Cd(10020)*Co(1047) + 
     &   (-4*Ce(411) + Co(1049) - 2*Cd(9949)*MGlpTmSt2(2))*
     &    UStop2(1,1)**2 + 
     &   2*Cd(9949)*MGlpTmSt2(2)*UStop2(1,1)*UStop2(1,2) - 
     &   2*Cd(9949)*MGlpTmSt2(2)*UStop2(1,2)**2

        Co(1152) = 4*Cd(10020)*Ce(411) - 
     &   (Co(1047) + Co(1049) - 2*Cd(9949)*MGlpTmSt2(2))*
     &    UStop2(1,1)**2 + Co(1049)*UStop2(1,1)*UStop2(1,2) - 
     &   Co(1049)*UStop2(1,2)**2

        Co(1453) = 2*MTy*Co(1318)*UUStop(1,3) + 
     &   AtC*(Co(1452) + AtC*UUStop(1,1)*UUStop(2,1)) - 
     &   At*Cd(14040)*UUStop(2,2) - 2*MTy*Cd(14040)*UUStop(2,3) - 
     &   MUEC*TB*Cd(11187)*UUStop(3,2)

        Co(1726) = 4*Ce(411) - Co(1049) - Cd(9948)*Co(1147) - 
     &   (Co(1047) - 2*Cd(9949)*MGlpTmSt2(2))*UStop2(1,1)

        Co(1755) = 4*Ce(411) - Co(1049) - Co(1047)*Cr(85) + 
     &   Cd(9949)*MGlpTmSt2(2)*(2 + 2*UStop2(1,2))

	Co(1821) = Cd(50)**2*Cd(9953) - Cd(47)**2*Cd(9954) + Co(1129)

        Co(185) = -((MUE**2*SA*TB + AtC*(AtC*CA + MUE*Co(139)))*
     &      UCStop(1,3)**2) + Co(165)*UCStopC(1,3)**2 + 
     &   (At*MUE - AtC*MUEC)*Co(138)*UStop2(1,1)*UStop2(1,2)

        Co(1865) = Co(1047) - 4*Ce(411)*Cr(85) - 
     &   2*Cd(9949)*MGlpTmSt2(2) + Co(1049)*(1 + UStop2(1,2))

        Co(1971) = (-(CA*MUE**2*TB) + 
     &      AtC*(AtC*SA - MUE*Co(1930)))*UCStop(1,3)**2 + 
     &   Co(1955)*UCStopC(1,3)**2 + 
     &   (At*MUE - AtC*MUEC)*Co(1929)*UStop2(1,1)*UStop2(1,2)

        Co(715) = Cd(9948)*MStop2(1) + 
     &   MGlpTmSt2(1)*UStop2(1,1) - MGlpT2*UStop2(1,2)

        Co(725) = 2*MUE*(MGlpT2 - MStop2(3))*UCStop(1,3)**2 + 
     &   MGlpT2*MUEC*UStop2(1,1)**2 + MUEC*MStop2(1)*UStop2(1,2)**2

	Co(733) = 5*MStop2(2)**2 + MGlpT2*(5*MGlpT2 + 8*MStop2(2))

	Co(738) = (2*MGlpT2 + MStop2(2))*(MGlpT2 + 2*MStop2(2))

        Co(1437) = MTy*MUE*TB*Cr(91)*UUStop(2,1) - 
     &   MUEC*TB*Cd(11017)*UStop(1,2)*UUStop(2,2) + 
     &   AtC*(Co(1435)*UStop(2,1) - AtC*UUStop(2,1)*UUStop(3,1)) + 
     &   At*(MUE*TB*Co(1436) + At*UUStop(2,2)*UUStop(3,2) - 
     &      2*MTy*UUStop(2,2)*UUStop(3,4))

        Co(1786) = 4*Ce(411) - Cd(9972)*Co(1047) - Co(1049) + 
     &   2*Cd(9949)*MGlpTmSt2(2)*(1 - Cd(9948) + UStop2(1,2))

        Co(1883) = -4*Cd(9972)*Ce(411) + Co(1047) - 
     &   Cd(9948)*Co(1049) - 2*Cd(9949)*MGlpTmSt2(2) + 
     &   Co(1049)*(1 + UStop2(1,2))

	Co(2071) = -(MGl2*Cr(75)) + Cd(9948)*MStop2(1) + Opt(34)

	Co(2072) = -Opt(34) + MGl2*(UStop2(1,1) - UStop2(1,2))

        Co(2075) = Opt(35) + 
     &   MUEC*Cd(12070)*Cd(12701)*UCStopC(1,3)**2*UStop2(1,1)

        Co(21) = 2*MTy*(AtC + MUE*TB)*UCStop(1,3) - 
     &   2*MTy*Cd(12701)*UCStopC(1,3) + 
     &   MUE*TB*(AtC*UCStop(1,3)**2 + 
     &      At*UStop2(1,1)*UStop2(1,2)) - 
     &   MUEC*TB*(At*UCStopC(1,3)**2 + AtC*UStop2(1,1)*UStop2(1,2))

        Co(610) = -2*M_3*MTy*Co(609) + 
     &   MUE**2*UCStop(1,3)**2*
     &    (MGlpT2 + 6*M_3C*MTy*UCStop(1,3)) - 
     &   MGlpT2*MUEC**2*UCStopC(1,3)**2 + 
     &   2*M_3C*MTy*MUEC**2*UCStopC(1,3)*UStop2(1,1)**2

        Co(60) = -2*MTy*(AtC + MUE*TB)*UCStop(1,3) + 
     &   2*MTy*Cd(12701)*UCStopC(1,3) + 
     &   MUE*TB*(AtC*UCStop(1,3)**2 + 
     &      At*UStop2(1,1)*UStop2(1,2)) - 
     &   MUEC*TB*(At*UCStopC(1,3)**2 + AtC*UStop2(1,1)*UStop2(1,2))

        Co(651) = Opt(35) - 
     &   MUEC*Cd(12701)*MGlpTmSt4(1)*UCStopC(1,3)**2*UStop2(1,1)

        Co(721) = Opt(35) + 
     &   MUEC*Cd(6458)*Cd(12701)*UCStopC(1,3)**2*UStop2(1,1)

        Co(1184) = Cd(50)**2*
     &    (-Cd(9951) + Cd(9953) + 
     &      2*Co(1040)*MGlpTmSt2(1)*UStop2(1,2)) + 
     &   Cd(47)**2*(-Cd(9952) + Cd(9954) + 
     &      2*Co(1040)*MGlpTmSt2(2)*UStop2(1,2))

        Co(2092) = MGlpT2*Co(2091) - 
     &   MUE**2*UCStop(1,3)**2*
     &    (Cd(9948)*MStop2(1) + 
     &      (2*MStop2(1) + MStop2(3))*UStop2(1,1)) - 
     &   MUEC**2*UCStopC(1,3)**2*
     &    (Cd(9948)*MStop2(1) + MStop2(3)*(1 + UStop2(1,2)))

        Co(716) = MGlpTmSt2(1)*UStop2(1,1) + 
     &   MGlpT2*(Cd(9948) - 2*UStop2(1,2)) + MStop2(2)*UStop2(1,2)

        Co(1711) = MGlmT2*Co(1047) - MGlmT2*Co(1049) + 
     &   2*Cd(50)*UStop2(1,1)*
     &    (Cd(9948)*Co(1044)*UStop2(1,1) + 
     &      UStop2(1,2)*(Co(1045) + Co(1044)*UStop2(1,2)))

        Co(1719) = Cd(10658)*Co(1047) - 
     &   4*Ce(411)*Cr(75)**2*UStop2(1,1) + 
     &   Co(1049)*UStop2(1,1)*(UStop2(1,1) - UStop2(1,2))**2 - 
     &   2*Cd(9949)*MGlpTmSt2(2)*
     &    (UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))))

        Co(2086) = Cd(9948)*Co(2084) + MGl2*Co(2085) + 
     &   2*MUE*UCStop(1,3)**2*
     &    (MTy2 + (MTy2 - 2*MStop2(2) + MStop2(3))*UStop2(1,1)) + 
     &   MUEC*UStop2(1,1)*UStop2(1,2)*
     &    (MTy2*UStop2(1,1) - 2*(MTy2 - 2*MStop2(3))*UStop2(1,2))

        Co(461) = MGl2*(MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*UStop2(1,2)) + 
     &   2*MTy*(MUE*UCStop(1,3)*
     &       (M_3C*UCStop(1,3)**2 + 
     &         M_3*UStop2(1,1)*UStop2(1,2)) + 
     &      MUEC*(M_3C*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) + 
     &         M_3*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2))) + 
     &   MTy2*(MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*
     &       (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Co(656) = MGl2**2*Cd(10388) - 2*MGl2*Co(655) + 
     &   MUEC*(MTy2 - MStop2(2))**2*UCStopC(1,3)**2 + 
     &   MUE*UStop2(1,1)*(Cd(9948)*Co(653)*UStop2(1,1) - 
     &      2*MTy2*MStop2(2)*UStop2(1,2) + 
     &      (MTy2**2 + MStop2(2)**2)*UStop2(1,2)**2)

        Co(2343) = Cd(10658)*Co(1047) - 
     &   2*Cd(9949)*MGlpTmSt2(2) - 
     &   4*Cd(9958)*Ce(411)*UStop2(1,1) + 
     &   Co(1049)*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))))

        Co(604) = Opt(36) + MUEC**2*MStop2(2)*UCStopC(1,3)**2 - 
     &   2*M_3*MTy*MUEC**2*UCStopC(1,3)**3 + 
     &   MUE**2*UCStop(1,3)*
     &    (UCStop(1,3)*(MGl2 - MStop2(2) + 
     &         2*M_3C*MTy*UCStop(1,3)) + 
     &      2*M_3*MTy*UStop2(1,1)*UStop2(1,2))

        Co(608) = 2*M_3*MTy*Co(607) - Opt(36) - 
     &   MUE**2*UCStop(1,3)**2*
     &    (6*M_3C*MTy*UCStop(1,3) + MGl2*UStop2(1,1)) + 
     &   MGl2*MUEC**2*UCStopC(1,3)**2*UStop2(1,2) - 
     &   2*M_3C*MTy*MUEC**2*UCStopC(1,3)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2)

        Co(623) = MUEC*MGlpTmSt2(2)*UCStopC(1,3)**2 + 
     &   2*M_3*MTy*MUEC*UCStopC(1,3)**3 + 
     &   2*MTy*(M_3*MUE + M_3C*MUEC)*UCStopC(1,3)*UStop2(1,1)*
     &    UStop2(1,2) + MUE*UStop2(1,1)*UStop2(1,2)*
     &    (MGl2 + 2*M_3C*MTy*UCStop(1,3) + 
     &      (MTy2 - MStop2(2))*UStop2(1,2))

        Co(2102) = Opt(37) + 
     &   2*Cd(47)**2*Cd(50)*Cd(9949)*Cd(9950)*MGlpTmSt2(2)*
     &    MStop2(3)*(Cd(9948) - UStop2(1,2))

        Co(2112) = Opt(37) + 
     &   2*Cd(47)**2*Cd(50)*Cd(9949)*Cd(9950)*MGlpTmSt2(2)*
     &    MStop2(3)*(-Cd(9948) + UStop2(1,2))

        Co(1762) = Co(1047) - Cd(9948)*Cd(9959)*Co(1130) - 
     &   2*Cd(9949)*MGlpTmSt2(2) - 
     &   (4*Ce(411) - Co(1049))*UStop2(1,1)*
     &    (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &      3*UStop2(1,2)**2)

        Co(719) = MGlpT2*Cd(9948)*Co(718) - 
     &   MUEC*UStop2(1,2)**2*
     &    (MGlpTmSt2(1)*UStop2(1,1) + MStop2(1)*UStop2(1,2)) + 
     &   2*MUE*UCStop(1,3)**2*
     &    (MGlpTmSt2(1)*UStop2(1,1) - 2*MGlpT2*UStop2(1,2) + 
     &      MStop2(2)*UStop2(1,2))

        Co(742) = -(MGl2*(MGl2*Co(739) + 2*Co(741))) - 
     &   MUEC*(MTy2 - MStop2(2))**2*UCStopC(1,3)**2 + 
     &   MUE*UStop2(1,1)*(2*Cd(9948)*Co(738)*UStop2(1,1) - 
     &      (MTy2**2 + MStop2(2)**2 + 
     &         4*(MTy2 + MStop2(2))**2*UStop2(1,1))*UStop2(1,2) + 
     &      2*MTy2*MStop2(2)*UStop2(1,2)**2)

        Co(1117) = 4*Cd(47)*Cd(50)*Ce(411)*
     &    (Cd(9948) + UStop2(1,1)) + 
     &   Cd(9950)*(Cd(47)*(2*(Co(1049) - Cd(9949)*Co(1115)) + 
     &         Co(1047)*(-2 + Cd(9948) + UStop2(1,1)) + 
     &         Ce(411)*(2 - 4*UStop2(1,1) - 4*UStop2(1,2))) + 
     &      Cd(50)*(Co(1034) + 2*Cd(9947)*Co(1116) + 
     &         Co(1037)*(-2 + Cd(9948) - UStop2(1,2)) + 
     &         Ce(411)*(-2 + 8*UStop2(1,1) + 8*UStop2(1,2))))

        Co(1209) = At**2*
     &    (Cd(50)**2*(Cd(9951) - 2*Co(1040)*MGlpTmSt2(1)) - 
     &      Cd(47)**2*(Cd(9952) - 2*Co(1040)*MGlpTmSt2(2)))*
     &    UCStopC(1,3)**2 + 
     &   Cd(47)*Cd(50)*(-2*Cd(50)*Cd(9947)*Cd(12594) + 
     &      2*Cd(47)*Cd(9949)*Cd(12596) - 
     &      Cd(50)*(Co(1034) - 2*Co(1037) + 
     &         Ce(411)*(-2 + 4*UStop2(1,1) + 4*UStop2(1,2))) + 
     &      Cd(47)*(Co(1047) - 2*Co(1049) + 
     &         Ce(411)*(-2 + 4*UStop2(1,1) + 4*UStop2(1,2))))

        Co(1145) = -2*Cd(47)*Cd(50)**2*Cd(9947)*Cd(12594)*
     &    UCStop(1,3) + 2*Cd(47)**2*Cd(50)*Cd(9949)*Cd(12596)*
     &    UCStop(1,3) - 2*Cd(47)**2*Cd(50)*Ce(411)*UCStop(1,3) + 
     &   2*Cd(47)*Cd(50)**2*Ce(411)*UCStop(1,3) - 
     &   2*Cd(47)*Cd(50)**2*Co(1034)*UCStop(1,3) + 
     &   Cd(47)*Cd(50)**2*Co(1037)*UCStop(1,3) + 
     &   2*Cd(47)**2*Cd(50)*Co(1047)*UCStop(1,3) - 
     &   Cd(47)**2*Cd(50)*Co(1049)*UCStop(1,3) + 
     &   AtC**2*Cd(50)**2*Cd(9951)*UCStop(1,3)**3 - 
     &   AtC**2*Cd(47)**2*Cd(9952)*UCStop(1,3)**3 - 
     &   2*AtC**2*Cd(50)**2*Co(1040)*MGlpTmSt2(1)*UCStop(1,3)**3 + 
     &   2*AtC**2*Cd(47)**2*Co(1040)*MGlpTmSt2(2)*UCStop(1,3)**3 + 
     &   4*Cd(47)**2*Cd(50)*Ce(411)*UCStop(1,3)*UStop2(1,1) - 
     &   4*Cd(47)*Cd(50)**2*Ce(411)*UCStop(1,3)*UStop2(1,1) + 
     &   Cd(47)*Cd(50)**2*Co(1034)*UCStop(1,3)*UStop2(1,1) + 
     &   Cd(47)*Cd(50)**2*Co(1037)*UCStop(1,3)*UStop2(1,1) - 
     &   Cd(47)**2*Cd(50)*Co(1047)*UCStop(1,3)*UStop2(1,1) - 
     &   Cd(47)**2*Cd(50)*Co(1049)*UCStop(1,3)*UStop2(1,1) + 
     &   4*Cd(47)**2*Cd(50)*Ce(411)*UCStop(1,3)*UStop2(1,2) - 
     &   4*Cd(47)*Cd(50)**2*Ce(411)*UCStop(1,3)*UStop2(1,2) + 
     &   Cd(47)*Cd(50)**2*Co(1034)*UCStop(1,3)*UStop2(1,2) + 
     &   Cd(47)*Cd(50)**2*Co(1037)*UCStop(1,3)*UStop2(1,2) - 
     &   Cd(47)**2*Cd(50)*Co(1047)*UCStop(1,3)*UStop2(1,2) - 
     &   Cd(47)**2*Cd(50)*Co(1049)*UCStop(1,3)*UStop2(1,2) + 
     &   At**2*(-(Cd(50)**2*
     &         (Cd(9951) - 2*Co(1040)*MGlpTmSt2(1))) + 
     &      Cd(47)**2*(Cd(9952) - 2*Co(1040)*MGlpTmSt2(2)))*
     &    UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Opt(32) = 2*MStop2(4)*
     &   ((MUE**2*SA*TB + AtC*(AtC*CA + MUE*Co(139)))*
     &      UCStop(1,3)**2 - Co(165)*UCStopC(1,3)**2)

        Opt(33) = 2*MStop2(4)*
     &   ((CA*MUE**2*TB + AtC*(-(AtC*SA) + MUE*Co(1930)))*
     &      UCStop(1,3)**2 - Co(1955)*UCStopC(1,3)**2)

        Opt(38) = MTy*(2*M_3C*MUE*(AtC + MUE*TB)*
     &      UCStop(1,3)**3 - 
     &     2*M_3*MUEC*Cd(12701)*UCStopC(1,3)**3 + 
     &     2*(-(M_3C*At*MUE) + AtC*(M_3*MUE + M_3C*MUEC) + 
     &        M_3*MUE**2*TB)*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)-
     &       2*Co(309)*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2))

        Opt(42) = Cd(9950)*Cd(10103)*
     &   (Cd(47)*(2*Ce(411) - 2*Co(1047) + Co(1049)) + 
     &     Cd(50)*(2*Co(1034) - Co(1037) + 
     &        Ce(411)*(-2 + 8*UStop2(1,2))))

        Opt(43) = Cd(9950)*Cd(10103)*
     &   (Cd(50)*(Co(1034) + Co(1037)) + 
     &     Cd(47)*(Co(1047) + Co(1049)))*UStop2(1,1)

        Opt(44) = (Cd(50)*
     &      (2*Ce(411) - 2*Co(1034) + Co(1037) + 
     &        (-4*Ce(411) + Co(1034) + Co(1037))*UStop2(1,1)) + 
     &     Cd(47)*(-2*Ce(411) + 2*Co(1047) - Co(1049) + 
     &        (4*Ce(411) - Co(1047) - Co(1049))*UStop2(1,1)))*
     &   (UStop2(1,1)**2 + UStop2(1,2)**2)

        Opt(45) = Opt(44) + 
     &   (Cd(50)*Co(1034) + 2*Cd(47)*Co(1049))*UStop2(1,1)**2*
     &    UStop2(1,2) + Cd(50)*Co(1034)*UStop2(1,2)**3 - 
     &   Cd(47)*Cd(9948)*Co(1049)*
     &    (3*UStop2(1,1)**2 + UStop2(1,2)**2)

        Opt(46) = 2*Cd(47)*Cd(50)**2*Cd(9947)*Cd(10026) - 
     &   2*Cd(47)**2*Cd(50)*Cd(9949)*Cd(10151) + 
     &   6*Cd(47)**2*Cd(50)*Ce(411)*UStop2(1,1) - 
     &   6*Cd(47)*Cd(50)**2*Ce(411)*UStop2(1,1) + 
     &   8*Cd(47)**2*Cd(50)*Cd(9948)*Ce(411)*UStop2(1,1) + 
     &   3*Cd(47)*Cd(50)**2*Co(1034)*UStop2(1,1) - 
     &   3*Cd(47)**2*Cd(50)*Co(1047)*UStop2(1,1) + 
     &   Cd(50)**2*Co(1056)*UStop2(1,1) - 
     &   Cd(50)**2*Co(1057)*UStop2(1,1) - 
     &   Cd(47)**2*Co(1061)*UStop2(1,1) + 
     &   Cd(47)**2*Co(1062)*UStop2(1,1) - 
     &   8*Cd(47)**2*Cd(50)*Ce(411)*UStop2(1,1)**2 + 
     &   8*Cd(47)*Cd(50)**2*Ce(411)*UStop2(1,1)**2 - 
     &   2*Cd(47)*Cd(50)**2*Co(1034)*UStop2(1,1)**2 - 
     &   2*Cd(47)*Cd(50)**2*Co(1037)*UStop2(1,1)**2 - 
     &   2*Cd(50)**2*Co(1042)*UStop2(1,1)**2 + 
     &   2*Cd(50)**2*Co(1043)*UStop2(1,1)**2 + 
     &   2*Cd(47)**2*Cd(50)*Co(1047)*UStop2(1,1)**2 + 
     &   2*Cd(47)**2*Cd(50)*Co(1049)*UStop2(1,1)**2 + 
     &   2*Cd(47)**2*Co(1050)*UStop2(1,1)**2 - 
     &   2*Cd(47)**2*Co(1051)*UStop2(1,1)**2 - 
     &   (Cd(47)*Cd(50)**2*
     &       (6*Ce(411) - 3*Co(1034) + 
     &         4*(-4*Ce(411) + Co(1034) + Co(1037))*UStop2(1,1)) + 
     &      Cd(50)**2*(-Co(1056) + Co(1057) + 
     &         4*(Co(1042) - Co(1043))*UStop2(1,1)) + 
     &      Cd(47)**2*(-6*Cd(50)*Ce(411) + 
     &         8*Cd(50)*Cd(9948)*Ce(411) + 3*Cd(50)*Co(1047) + 
     &         Co(1061) - Co(1062) + 
     &         4*(Cd(50)*(6*Ce(411) - Co(1047)) - Co(1050) + 
     &            Co(1051))*UStop2(1,1)))*UStop2(1,2) + 
     &   2*(Cd(47)*Cd(50)**2*(4*Ce(411) - Co(1034) - Co(1037)) + 
     &      Cd(50)**2*(-Co(1042) + Co(1043)) + 
     &      Cd(47)**2*(Cd(50)*Co(1047) + Co(1050) - Co(1051)))*
     &    UStop2(1,2)**2

        Opt(47) = -2*Cd(47)*Cd(50)**2*Cd(9947)*Cd(9950)*
     &    Cd(10026)*Cd(10103) + 
     &   2*Cd(47)**2*Cd(50)*Cd(9949)*Cd(9950)*Cd(10030)*
     &    Cd(10103) - 16*Cd(47)**2*Cd(50)**2*Cd(9948)*Cd(9950)*
     &    Ce(411)*UStop2(1,1) - 
     &   6*Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*Ce(411)*
     &    UStop2(1,1) + 6*Cd(47)*Cd(50)**2*Cd(9950)*Cd(10103)*
     &    Ce(411)*UStop2(1,1) + 
     &   Cd(47)**2*Cd(50)**2*Cd(10103)*Co(1034)*UStop2(1,1) - 
     &   2*Cd(47)*Cd(50)**2*Cd(9950)*Cd(10103)*Co(1034)*
     &    UStop2(1,1) + Cd(47)**2*Cd(50)**2*Cd(10103)*Co(1037)*
     &    UStop2(1,1) + Cd(47)*Cd(50)**2*Cd(9950)*Cd(10103)*
     &    Co(1037)*UStop2(1,1) + 
     &   2*Cd(47)**2*Cd(50)*Cd(9948)*Cd(9950)*Cd(10103)*Co(1047)*
     &    UStop2(1,1) - Cd(50)**2*Cd(9950)*Cd(10103)*Co(1056)*
     &    UStop2(1,1) + Cd(50)**2*Cd(9950)*Cd(10103)*Co(1057)*
     &    UStop2(1,1) + Cd(47)**2*Cd(9950)*Cd(10103)*Co(1061)*
     &    UStop2(1,1) - Cd(47)**2*Cd(9950)*Cd(10103)*Co(1062)*
     &    UStop2(1,1) + 8*Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*
     &    Ce(411)*UStop2(1,1)**2 - 
     &   8*Cd(47)*Cd(50)**2*Cd(9950)*Cd(10103)*Ce(411)*
     &    UStop2(1,1)**2 + 
     &   2*Cd(47)*Cd(50)**2*Cd(9950)*Cd(10103)*Co(1034)*
     &    UStop2(1,1)**2 + 
     &   2*Cd(47)*Cd(50)**2*Cd(9950)*Cd(10103)*Co(1037)*
     &    UStop2(1,1)**2 + 
     &   2*Cd(50)**2*Cd(9950)*Cd(10103)*Co(1042)*UStop2(1,1)**2 - 
     &   2*Cd(50)**2*Cd(9950)*Cd(10103)*Co(1043)*UStop2(1,1)**2 - 
     &   2*Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*Co(1047)*
     &    UStop2(1,1)**2 - 
     &   2*Cd(47)**2*Cd(50)*Cd(9950)*Cd(10103)*Co(1049)*
     &    UStop2(1,1)**2 - 
     &   2*Cd(47)**2*Cd(9950)*Cd(10103)*Co(1050)*UStop2(1,1)**2 + 
     &   2*Cd(47)**2*Cd(9950)*Cd(10103)*Co(1051)*UStop2(1,1)**2 + 
     &   (Cd(9950)*Cd(10103)*
     &       (Cd(47)*Cd(50)**2*
     &          (6*Ce(411) - 2*Co(1034) + Co(1037)) + 
     &         Cd(50)**2*(-Co(1056) + Co(1057)) - 
     &         Cd(47)**2*(Cd(50)*
     &             (6*Ce(411) + 2*Cd(9948)*Co(1047)) - Co(1061) + 
     &            Co(1062))) - 
     &      2*Cd(9950)*(16*Cd(47)**2*Cd(50)**2*Ce(411) + 
     &         Cd(10103)*(-2*Cd(47)*Cd(50)**2*
     &             (Co(1034) + Co(1037)) + 
     &            2*Cd(50)**2*(-Co(1042) + Co(1043)) + 
     &            Cd(47)**2*
     &             (Cd(50)*(16*Ce(411) + 3*Co(1047)) + 
     &               2*Co(1050) - 2*Co(1051))))*UStop2(1,1) + 
     &      Cd(47)**2*Cd(50)**2*Cd(10103)*
     &       (Co(1034) + Co(1037) + 16*Ce(411)*UStop2(1,1)))*
     &    UStop2(1,2) + 2*Cd(9950)*Cd(10103)*
     &    (Cd(47)*Cd(50)**2*(-4*Ce(411) + Co(1034) + Co(1037)) + 
     &      Cd(50)**2*(Co(1042) - Co(1043)) + 
     &      Cd(47)**2*(4*Cd(50)*Ce(411) - Co(1050) + Co(1051)))*
     &    UStop2(1,2)**2

        Co(98) = 2*MTy*(Co(97)*UCStop(1,3) - 
     &     2*Cd(12701)*Cr(75)*UCStopC(1,3))

        Co(1005) = Co(1004)*MStop2(2)*UStop2(1,2) + 
     &   Co(995)*Cr(9)*MStop2(1)*UStopC(1,2)

	Co(1014) = 2*Co(808)*UUStopC(3,2) - Co(1013)*UUStopC(3,4)

        Co(1016) = Co(962)*UCStopC(1,3)*UStop2(1,1) - 
     &   Co(964)*UUStopC(1,2)

        Co(1017) = Co(1016)*UStopC(1,1) + 
     &   MUE*Co(963)*UCStop(1,3)*UStop2(1,2)*UStopC(1,2)

        Co(1020) = 2*MUEC*Co(954)*MStop2(3) + 
     &   Cd(11177)*MStop2(2)*UCStop(1,3)

	Co(1023) = Co(1014)*MStop2(2) - 2*Co(948)*MStop2(3)

	Co(1054) = Co(1052) - Co(1053)*Cr(13)

	Co(1055) = 4*MTy2*A0delStop(1) - 4*Co(1054)*MStop2(1)

	Co(1059) = Co(1052) - Co(1053)*Cr(34)

	Co(1060) = 4*MTy2*A0delStop(2) - 4*Co(1059)*MStop2(2)

	Co(1074) = Cd(47)**2*Cd(50)**2*Co(1071) + Cd(9950)*Co(1073)

        Co(1076) = Cd(9948)*Co(1074) - 
     &   Cd(47)*Cd(50)*Cd(9950)*Co(1075)*UStop2(1,2)

	Co(1084) = MGlmT2*Cd(47)**2*Co(1034) + Cd(9950)*Co(1083)

	Co(1085) = AtC*Co(1084)*UCStop(1,3)**2 + At*Co(1082)*UStop2(1,1)

	Co(1089) = Co(1062) + Cd(50)*Co(1088)

	Co(1092) = 2*Co(1062) + Cd(50)*Co(1088)

        Co(1105) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Cd(9966)*Co(1104) + 
     &        Co(1103)*MStop2(3)*UStop2(1,1)*UStop2(1,2)))/
     &    (SB*SW2*Cd(47)**2*MStop2(3))

	Co(1106) = 2*Co(1040)*MGlpTmSt2(1) - Cd(9951)*UStop2(1,2)

	Co(1107) = 2*Co(1040)*MGlpTmSt2(2) - Cd(9952)*UStop2(1,2)

	Co(1108) = Cd(50)**2*Co(1106) + Cd(47)**2*Co(1107)

	Co(1109) = Cd(9948)*Ce(451) + Co(1108)*UStop2(1,2)

        Co(110) = -4*Co(109)*Cr(75)*MStop2(3) + 
     &   MUE*TB*Co(108)*UCStop(1,3)

	Co(1114) = AtC*Co(1113)*UCStop(1,3)**2 + At*Co(1109)*UStop2(1,1)

	Co(1121) = AtC*Co(1120)*UCStop(1,3)**2 + At*Co(1101)*UStop2(1,1)

        Co(1122) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(50)*Co(1067)*Co(1119) - 
     &        Co(1121)*MStop2(3)))/(SB*SW2*TB*Cd(50)**2*MStop2(3))

        Co(1123) = -(cI/(16.D0*Pi**3)*
     &     (Alfas2L*(2*MGlmT2*Cd(50)*Co(1067) + 
     &          Co(1120)*MStop2(3))*UCStop(1,3)**2)/
     &      (SB*SW2*Cd(50)**2*MStop2(3)))

        Co(1124) = -(cI/(16.D0*Pi**3)*
     &     (Alfas2L*(2*MGlmT2*Cd(47)*Cd(10077)*Co(1104) + 
     &          Co(1103)*MStop2(3))*UCStop(1,3)**2)/
     &      (SB*SW2*Cd(47)**2*MStop2(3)))

	end
