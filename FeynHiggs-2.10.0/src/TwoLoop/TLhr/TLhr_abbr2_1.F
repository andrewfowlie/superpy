* TLhr_abbr2_1.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr2_1
	implicit none

#include "TLhr.h"


        Opt(22) = AtC**2*UCStop(1,3)**2 + 
     &   At*(At*UCStopC(1,3)**2 + 2*AtC*UStop2(1,1)*UStop2(1,2)**2)

        Opt(24) = AtC**2*Cd(10032) + 
     &   M_3*At**2*UCStopC(1,3)**3 - 
     &   2*M_3C*MUE2*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) + 
     &   (M_3C*At**2 - 2*M_3*MUE2)*UCStopC(1,3)*UStop2(1,1)*
     &    UStop2(1,2)

        Opt(25) = At*AtC*Cd(10658) - At**2*UCStopC(1,3)**2 + 
     &   UCStop(1,3)**2*(MUE**2 - AtC**2*UStop2(1,1))

        Opt(26) = -3*MUE**2*UCStop(1,3)**3 + 
     &   2*MUE2*UCStop(1,3)*UStop2(1,1)**2 + 
     &   (-At**2 + MUEC**2)*UCStopC(1,3)*UStop2(1,1)**2

        Opt(27) = (At - MUEC)*(At + MUEC)*UCStopC(1,3)**2 + 
     &   UCStop(1,3)**2*(-MUE**2 + AtC**2*UStop2(1,1))

        Opt(28) = At*AtC*Cd(9970) + 2*MStop2(3) + 
     &   MUE2*UStop2(1,1)*(UStop2(1,1) - UStop2(1,2))**2 + 
     &   Cd(9948)*(AtC**2*UCStop(1,3)**2 + 
     &      MUEC**2*UCStopC(1,3)**2 + 
     &      MUE2*(3*UStop2(1,1)**2 + UStop2(1,2)**2))

        Opt(29) = Opt(22) - 8*M_3C*AtC*UCStop(1,3)**2 - 
     &   MUE**2*UCStop(1,3)**2 - 
     &   (8*M_3*At + MUEC**2)*UCStopC(1,3)**2 + 
     &   2*AtC*UStop2(1,1)*(-4*M_3 + At*UStop2(1,1))*
     &    UStop2(1,2) - 2*UStop2(1,1)*UStop2(1,2)*
     &    (4*M_3C*At + MUE2*UStop2(1,2))

        Opt(30) = 4*Cd(179) - 
     &   4*MStop2(2)*(UStop2(1,1)**3 + UStop2(1,2)**3 + 
     &      UStop2(1,1)*(-1 + 3*UStop2(1,2)**2))

        Opt(31) = 4*MGl2 + Opt(27) - 
     &   8*M_3*At*UCStopC(1,3)**2 + 
     &   AtC*UCStop(1,3)**2*(-8*M_3C + AtC*UStop2(1,2)) + 
     &   2*At*(-4*M_3C + AtC)*UStop2(1,1)*
     &    (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2) - 
     &   2*UStop2(1,1)*(MUE2*UStop2(1,2) + 
     &      4*M_3*AtC*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

	Ce(1) = At*AtC - MUE2

	Ce(1002) = MUEC**2*UStop2(1,1) + At**2*UStop2(1,2)

	Ce(1005) = MUE2*UStop2(1,1) + 2*MTy2*UStop2(1,2)

        Ce(1006) = 2*Ce(1005)*UCStopC(1,3) + 
     &   MUE**2*UCStop(1,3)*UStop2(1,1)

	Ce(1007) = M_3*Ce(1006) + 2*At*Cd(10493)*UCStopC(1,3)

	Ce(1009) = At**2*UStop2(1,1) + MUEC**2*UStop2(1,2)

        Ce(1024) = 2*M_3C*At*MGlmT2*Cd(9948)*MStop2(3) - 
     &   MUE2*Cd(47)*Cd(9969)*UStop2(1,1)

	Ce(1026) = 2*Cd(47)**2 - Cd(10284)*MStop2(3)

	Ce(1027) = 2*M_3*MGlmT2*Cd(47) + At*Cd(10284)*UStop2(1,1)

	Ce(1031) = 2*(2*Cd(9972)*MStop2(1) + MStop2(3))

	Ce(1038) = MGlmT2*Cd(10454) - Cd(10455)*MStop2(1)

	Ce(1040) = 4*M_3 + At*Cd(10306)

	Ce(1041) = 4*M_3C*UCStop(1,3) + Ce(1040)*UCStopC(1,3)

	Ce(1045) = MGlmT2*Cd(10321) - Cd(10329)*MStop2(1)

	Ce(1050) = -(MUE2*Cd(10091)) + MStop2(3)

        Ce(1053) = 2*At*Cd(10018)*Cr(5) + 
     &   Cd(10025)*MStop2(3)*UCStop(1,3)

        Ce(1054) = MGlmT2*Ce(1053) + 
     &   2*Cd(47)*Cd(9955)*UCStop(1,3)*UStop2(1,1)

	Ce(1055) = -2*Cd(47)*Ce(1054) + Cd(10279)*MStop2(3)*UCStop(1,3)

	Ce(1058) = 3*Cd(47) + Cd(10278)

	Ce(1061) = 2*Cd(47)*Cd(152)*MStop2(1) + Cd(10280)*MStop2(3)

	Ce(1064) = 2*MGlmT2*Cd(47) - Cd(10015)*Cd(10278)

	Ce(1071) = At*Cd(164)*Cd(2597) + 8*M_3*Cd(12298)

	Ce(1074) = MGlpT2 - 2*M_3*MTy*UCStopC(1,3)

        Ce(1075) = -(Ce(1074)*UCStopC(1,3)) + 
     &   2*M_3C*MTy*UStop2(1,1)*UStop2(1,2)

	Ce(1084) = 4*M_3*AtC*Cd(118) - MGlmT2*MUE2*Cd(9969)

	Ce(1089) = 8*M_3*AtC*Cd(12370) - MUE2*Cd(12371)*UStop2(1,2)

	Ce(1090) = 16*M_3C*Cd(12370) + AtC*Cd(12371)

	Ce(1091) = 8*M_3C*Cd(12370) + AtC*Cd(12371)

        Ce(1092) = 8*M_3*Cd(12370)*UCStopC(1,3)**2 + 
     &   Ce(1091)*UStop2(1,1)*UStop2(1,2)

	Ce(11) = 2*(MTy*Cd(492) + AtC*UCStop(1,3))

	Ce(110) = C2A*Cd(424) - S2A*Ce(1)

	Ce(1098) = 4*Cd(47)*Cd(164)*MStop2(1) + Cd(12371)*MStop2(3)

	Ce(1103) = MGlmT2*Cd(9972) + MStop2(1)

	Ce(1104) = 4*Cd(47)*Ce(1103)*MStop2(1) + Cd(12371)*MStop2(3)

        Ce(1111) = 4*Cd(47)*Cd(164)*Cd(10019)*MStop2(1) - 
     &   Cd(10086)*Cd(12371)*MStop2(3)

        Ce(1112) = Cd(12371)*MStop2(3) - 
     &   4*Cd(47)*Cd(164)*Cd(9955)*MStop2(1)*UStop2(1,2)

        Ce(1117) = 4*Cd(47)*Cd(164)*Cd(10042)*MStop2(1) - 
     &   Cd(10039)*Cd(12371)*MStop2(3)

        Ce(1118) = 2*Cd(47)**2*UCStop(1,3) + 
     &   At*Cd(12371)*Cr(5)*UStop2(1,1)*UStop2(1,2)

        Ce(1120) = At*Cd(10091)*Cr(5) - 
     &   MGlmT2*Cd(9955)*UCStop(1,3)*UStop2(1,2)

	Ce(1124) = At*MGlmT2*Cr(5) + Cd(47)*UCStop(1,3)

	Ce(1131) = 5*Cd(118) + 2*MUE2*UStop2(1,1)**2

	Ce(1135) = Cd(9972)*MStop2(1) + MStop2(3)

	Ce(1136) = Cd(47)*Cd(9972) + Cd(164)*Ce(1135)

	Ce(1137) = 2*MStop2(1) + Cd(11731)*UStop2(1,1)

	Ce(1139) = Cd(47) + 2*MUE2*Cd(152)*UStop2(1,1)**2

	Ce(1140) = Cd(47)*Cd(164) - Cd(152)*MStop2(1)*MStop2(3)

	Ce(1150) = MUE2*Cd(10018) + MStop2(3)

        Ce(1153) = 2*At*Cd(10091)*Cr(5) - 
     &   Cd(10025)*MStop2(3)*UCStop(1,3)

        Ce(1154) = MGlmT2*Ce(1153) + 
     &   2*Cd(47)*Cd(9955)*UCStop(1,3)*UStop2(1,1)

	Ce(1157) = 6*MGl2 - MUE2

	Ce(1158) = 2 - 6*UStop2(1,1)**2

	Ce(1159) = Ce(1158)*MStop2(1) + Ce(1157)*UStop2(1,1)**2

        Ce(1175) = Cd(50)**2*MGlpTmSt2(1)*MStop2(1) + 
     &   Cd(47)**2*MGlpTmSt2(2)*MStop2(2)

        Ce(1180) = M_3*Cd(12594)*UCStopC(1,3)**2 + 
     &   M_3C*MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2)

        Ce(1181) = M_3*Cd(12596)*UCStopC(1,3)**2 + 
     &   M_3C*MGlpTmSt2(2)*UStop2(1,1)*UStop2(1,2)

        Ce(1186) = -2*Cd(50)*Cd(11942) + 
     &   MGlpTmSt2(2)*MStop2(2)*MStop2(3)

        Ce(1188) = -2*Cd(50)*Cd(9972)*MStop2(1) + 
     &   MGlpTmSt2(2)*MStop2(2)*MStop2(3)

	Ce(1192) = Cd(47) + 2*Cd(50)*UStop2(1,1)*UStop2(1,2)

	Ce(1195) = MUE2*Cr(75)**2 + AtC**2*UCStop(1,3)**2

	Ce(12) = Cd(425) - Cd(427)

	Ce(1196) = MUE2*Cd(9983) + AtC**2*UCStop(1,3)**2

	Ce(1202) = MGlpTmSt2(2)*MStop2(3) - 2*Cd(50)*UStop2(1,1)

	Ce(1204) = 2*Cd(50) - MGlpTmSt2(2)*MStop2(3)

	Ce(1219) = M_3*Cd(12036) - 2*M_3C*Cd(12029)*UCStop(1,3)**2

        Ce(1220) = M_3*Cd(12029)*Cr(28) - 
     &   2*M_3C*Cd(12030)*UCStop(1,3)**2

	Ce(1221) = Cd(47)*Ce(1219) - Cd(50)*Ce(1220)

	Ce(1228) = Cd(47)*Cd(10021) + Cd(50)*Cd(10455)

	Ce(1229) = Cd(47)*Cd(10020) + Cd(50)*Cr(75)*UStop2(1,1)

	Ce(1231) = 2*AtC*UCStop(1,3) + At*UCStopC(1,3)

        Ce(1232) = 2*MUE2*UCStop(1,3)*UStop2(1,1) + 
     &   At*Ce(1231)*UStop2(1,2)

	Ce(126) = MStop2(1) - Cr(1)*UStop2(1,2)

	Ce(1261) = Cd(50)*Cr(28) - MGlmT2*Cd(10070)*MStop2(1)

	Ce(1268) = 4*Cd(9948) + UStop2(1,1)

	Ce(1269) = UStop2(1,1)**2 - Cd(10615)*UStop2(1,2)

	Ce(1270) = MUE2*Ce(1269) + 2*M_3*AtC*UStop2(1,1)

	Ce(1280) = Cd(50)*Cd(179) + Cd(251)*MStop2(2)*MStop2(3)

	Ce(1281) = Cd(50)*Cd(179)*Cd(9972) + Cd(251)*MStop2(2)*MStop2(3)

	Ce(129) = 2*Cd(640)*MStop2(2) - Cd(523)*UStop2(1,2)

	Ce(1288) = 2*MGlmT2*Cd(50) - Cd(179)*Cd(251)*Cd(10015)

	Ce(1293) = M_3*Cd(10649) + 2*M_3C*UCStop(1,3)**2

        Ce(130) = Cd(874)*UStop2(1,1)**3 - 
     &   Cr(75)*MStop2(4)*UStop2(1,2)**2

	Ce(1298) = 2*MStop2(2) + MUE2*UStop2(1,1)**2

	Ce(1299) = MStop2(2) + 2*MUE2*UStop2(1,1)**2

	Ce(1302) = 4*MStop2(2) - 2*MStop2(3)

	Ce(1303) = -(MGlmT2*Cd(11942)) + Cd(9969)*MStop2(1)*MStop2(2)

	Ce(131) = -(Ce(130)*MStop2(3)) + Ce(129)*UStop2(1,1)**2

	Ce(1305) = MGlmT2*Cd(9972)*MStop2(1) - Cd(10745)*MStop2(2)

	Ce(1311) = MGlmT2*Cd(10759) - Cd(10021)*MStop2(2)

	Ce(132) = Cr(85)*MStop2(2) - 2*MStop2(1)*UStop2(1,2)

	Ce(1316) = MGlmT2*Cd(10764) - Cd(10043)*MStop2(2)

	Ce(1319) = 3*Cd(50) + Cd(10602)

	Ce(1321) = 2*Cd(50)*Cd(179)*MStop2(2) - Cd(10604)*MStop2(3)

	Ce(1322) = MGlmT2*Cd(9972) + Cd(10077)*MStop2(2)

	Ce(1323) = -2*Cd(50)*Ce(1322)*MStop2(2) + Cd(10604)*MStop2(3)

	Ce(134) = Ce(132) - Cr(1)*UStop2(1,1)

	Ce(1335) = At*Cd(251)*Cd(2681) + 8*M_3*Cd(4187)

	Ce(135) = Ce(134)*UStop2(1,1) - MStop2(3)*UStop2(1,2)**2

	Ce(1353) = 8*Cd(491)*MStop2(3) - MUE2*UStop2(1,1)

	Ce(1354) = MUE2 - 8*MStop2(3)

        Ce(1356) = MUE2*Cd(9981)*Cd(12483) - 
     &   8*M_3*AtC*Cd(12482)*UStop2(1,2)

	Ce(1357) = 16*M_3C*Cd(12482) + AtC*Cd(12483)

        Ce(1358) = -(AtC*Ce(1357)*UCStop(1,3)**2) + 
     &   2*Ce(1356)*UStop2(1,1)

        Ce(1364) = 4*Cd(50)*Cd(251)*Cd(10021)*MStop2(2) + 
     &   Cd(10086)*Cd(12483)*MStop2(3)

        Ce(1365) = 4*Cd(50)*Cd(251)*Cd(10043)*MStop2(2) + 
     &   Cd(10039)*Cd(12483)*MStop2(3)

	Ce(1376) = 16*Cd(12482)*MStop2(2) + At*AtC*Cd(12483)*UStop2(1,1)

	Ce(14) = Ce(11) + 2*At*UCStopC(1,3)

	Ce(1398) = Cd(50) + Cd(251)*Cd(9969)*MStop2(1)

	Ce(1399) = Cd(50) + Cd(251)*Cd(10745)

	Ce(1400) = 2*MStop2(2) + Cd(11847)*UStop2(1,1)

	Ce(1402) = 4*Cd(491)*MStop2(3) - MUE2*UStop2(1,1)

	Ce(1403) = MUE2 - 4*MStop2(3)

	Ce(1407) = MGlmT2*Cd(12506) + Cd(50)*Cr(28)

	Ce(1408) = AtC*Ce(1407) + 2*MGlmT2*Cd(10039)*MStop2(3)

	Ce(1410) = 5*Cd(138) + 4*MUE2*UStop2(1,1)**2

	Ce(1411) = Cd(138)*Cd(11858) + MUE2*UStop2(1,1)**2

	Ce(142) = MStop2(1)*MStop2(2) + 2*MStop2(4)**2*UStop2(1,2)

	Ce(1413) = 4*Cd(138) - MUE2*Cd(10615)

	Ce(143) = MStop2(1)*MStop2(2) + Cd(832)*UStop2(1,2)

	Ce(1423) = 2*M_3 + At*Cd(10029)

	Ce(1431) = MGlmT2 + Cd(9972)*MStop2(2)

	Ce(1432) = -(Ce(1158)*MStop2(2)) + MUE2*UStop2(1,1)**2

	Ce(144) = 9*MStop2(2)**2 - Cd(637)*MStop2(3)

	Ce(1435) = M_3*Cd(11878) - 4*M_3C*UCStop(1,3)**2

	Ce(1442) = Cd(239) + MTy*Cd(164)*Cr(5)

	Ce(145) = MStop2(1)**2 + Ce(144)*UStop2(1,2)

	Ce(1443) = 2*MStop2(3)*UUStop(1,1) - 2*MStop2(1)*UUStop(2,1)

	Ce(1444) = 2*MStop2(3)*UUStop(1,1) - MStop2(1)*UUStop(2,1)

	Ce(1445) = Cd(45)*Ce(1443) + MGl2*MStop2(1)*UUStop(2,1)

	Ce(1446) = Cd(44)**2*Ce(1444) - MGl2*Ce(1445)

	Ce(1447) = 2*MStop2(3)*UUStop(1,2) - 2*MStop2(1)*UUStop(2,2)

	Ce(1448) = 2*MStop2(3)*UUStop(1,2) - MStop2(1)*UUStop(2,2)

	Ce(1449) = Cd(45)*Ce(1447) + MGl2*MStop2(1)*UUStop(2,2)

	Ce(1450) = Cd(44)**2*Ce(1448) - MGl2*Ce(1449)

	Ce(1451) = 2*MStop2(3)*UUStop(1,3) - 2*MStop2(1)*UUStop(2,3)

	Ce(1452) = 2*MStop2(3)*UUStop(1,3) - MStop2(1)*UUStop(2,3)

	Ce(146) = Ce(145) - MStop2(1)**2*UStop2(1,1)

	Ce(1453) = Cd(45)*Ce(1451) + MGl2*MStop2(1)*UUStop(2,3)

	Ce(1454) = Cd(44)**2*Ce(1452) - MGl2*Ce(1453)

	Ce(1457) = 2*MGl2 + Cd(44)

	Ce(1458) = MGl2*MGlmT2 - Ce(1457)*MStop2(1)

	Ce(1459) = Ce(1458) - MTy*Cd(152)*Cr(5)

	Ce(1460) = 2*MStop2(3)*UUStop(1,1) + Cd(44)*UUStop(2,1)

	Ce(1461) = MGlmT2 - 2*MStop2(1)

	Ce(1463) = 2*MStop2(3)*UUStop(1,2) + Cd(44)*UUStop(2,2)

	Ce(1465) = 2*MStop2(3)*UUStop(1,3) + Cd(44)*UUStop(2,3)

	Ce(1469) = Cd(253) - MTy*Cd(251)*Cr(5)

	Ce(1470) = MStop2(2)*UUStop(1,1) + MStop2(3)*UUStop(2,1)

	Ce(1471) = MStop2(2)*UUStop(1,1) + 2*MStop2(3)*UUStop(2,1)

	Ce(1472) = 2*Cd(49)*Ce(1470) - MGl2*MStop2(2)*UUStop(1,1)

	Ce(1473) = Cd(48)**2*Ce(1471) - MGl2*Ce(1472)

	Ce(1474) = MStop2(2)*UUStop(1,2) + MStop2(3)*UUStop(2,2)

	Ce(1475) = MStop2(2)*UUStop(1,2) + 2*MStop2(3)*UUStop(2,2)

	Ce(1476) = 2*Cd(49)*Ce(1474) - MGl2*MStop2(2)*UUStop(1,2)

	Ce(1477) = Cd(48)**2*Ce(1475) - MGl2*Ce(1476)

	Ce(1478) = MStop2(2)*UUStop(1,3) + MStop2(3)*UUStop(2,3)

	Ce(1479) = MStop2(2)*UUStop(1,3) + 2*MStop2(3)*UUStop(2,3)

	Ce(1480) = 2*Cd(49)*Ce(1478) - MGl2*MStop2(2)*UUStop(1,3)

	Ce(1481) = Cd(48)**2*Ce(1479) - MGl2*Ce(1480)

	Ce(1484) = 2*MGl2 + Cd(48)

	Ce(1485) = MGl2*MGlmT2 - Ce(1484)*MStop2(2)

	Ce(1486) = Ce(1485) + MTy*Cd(179)*Cr(5)

	Ce(1487) = Cd(48)*UUStop(1,1) - 2*MStop2(3)*UUStop(2,1)

	Ce(1488) = Cd(48)*Ce(1487) + MGl2*Cd(1909)*UUStop(1,1)

	Ce(1489) = MGlmT2 - 2*MStop2(2)

        Ce(1490) = Ce(1488)*MStop2(2) + 
     &   2*MGl2*Ce(1489)*MStop2(3)*UUStop(2,1)

	Ce(1491) = Cd(48)*UUStop(1,2) - 2*MStop2(3)*UUStop(2,2)

	Ce(1492) = Cd(48)*Ce(1491) + MGl2*Cd(1909)*UUStop(1,2)

	Ce(149) = Cr(85)*MStop2(1)*MStop2(2) + 2*Cd(832)*UStop2(1,2)**2

        Ce(1493) = Ce(1492)*MStop2(2) + 
     &   2*MGl2*Ce(1489)*MStop2(3)*UUStop(2,2)

	Ce(1494) = Cd(48)*UUStop(1,3) - 2*MStop2(3)*UUStop(2,3)

	Ce(1495) = Cd(48)*Ce(1494) + MGl2*Cd(1909)*UUStop(1,3)

        Ce(1496) = Ce(1495)*MStop2(2) + 
     &   2*MGl2*Ce(1489)*MStop2(3)*UUStop(2,3)

	Ce(1499) = 2*MTy*UUStop(3,1) + At*UUStop(3,4)

        Ce(1500) = Cd(9948)*Ce(1499)*UCStop(3,3) - 
     &   Cd(10981)*UUStop(2,1)*UUStopC(2,2)

        Ce(1503) = 2*MTy2*UStop(1,2)*UUStop(3,1) + 
     &   At*MTy*UStop(1,1)*UUStop(3,2)

	Ce(1508) = 2*MTy*UUStop(1,1) - At*UUStop(1,3)

	Ce(1509) = 4*MTy2 + MUE2

	Ce(1514) = Cd(9948)*UStop(1,1) + UCStop(3,1)*UStop(2,1)

	Ce(1518) = 2*MTy*UUStop(3,1) + At*UUStop(3,3)

	Ce(1519) = 2*MTy*UUStop(2,1) - At*UUStop(2,3)

	Ce(152) = MStop2(2)**2 + Cd(535)*MStop2(3)

	Ce(1528) = MUE2*Cd(10925) - 6*At*MTy*UUStop(3,2)

	Ce(1535) = Cd(698)*Cd(9809) - MUEC**2*Cr(91)*UStop(2,2)

	Ce(1537) = At*UUStop(3,2) - 2*MTy*UUStop(3,3)

	Ce(1538) = -(At*Ce(1537)) + 8*MTy2*UUStop(3,1)

	Ce(1539) = -(MUE2*Cd(10925)) + AtC*Cd(11334)

	Ce(154) = 2*MStop2(1)*MStop2(2) - 2*Cd(1013)*UStop2(1,2)

        Ce(1548) = MUE**2*Cr(91)*UStop(2,1) - 
     &   8*MTy2*UStop(1,2)*UUStop(2,2)

	Ce(1549) = 4*MTy*UStop(2,1) + 3*At*UStop(2,2)

        Ce(1551) = At*Cd(11106)*UStop(2,2) + 
     &   2*MTy*UStop(1,2)*UUStop(2,1)

        Ce(1556) = 4*Cd(9948)*Ce(1499)*UCStop(3,3) - 
     &   4*Cd(10981)*UUStop(2,1)*UUStopC(2,2)

        Ce(1557) = -(Ce(1556)*MStop2(3)) + 
     &   2*Cd(11181)*MStop2(2)*UCStop(1,3)*UStopC(2,2)

	Ce(1558) = Cd(11056) + MUEC**2*UUStop(1,2)

	Ce(1560) = 2*MTy*Cd(11106)*UStop(1,1) + At*Cd(10925)*UStop(1,2)

	Ce(156) = Ce(154) - 2*Cd(852)*UStop2(1,1)

        Ce(1563) = (Ce(1557) + 
     &     2*AtC*MStop2(2)*UCStop(1,3)**2*UUStop(1,1))/MStop2(2)

	Ce(1565) = UStop2(1,1)*UUStop(1,2) - UCStop(3,1)*UUStop(3,2)

	Ce(1566) = -(At**2*Ce(1565)) + 2*MUE2*UStop2(1,2)*UUStop(1,1)

	Ce(1568) = UStop2(1,1)*UUStop(1,2) + UCStop(3,1)*UUStop(3,2)

	Ce(1569) = UStop2(1,2)*UUStop(1,3) + UCStop(3,1)*UUStop(3,4)

	Ce(1575) = Cd(5519)*Cd(9948) - 4*MStop2(1)*MStop2(2)*UStop2(1,1)

        Ce(1576) = 4*Ce(1575)*UStop(1,1) - 
     &   4*MStop2(3)**2*UCStop(3,1)*UStop(2,1)

        Ce(1577) = 2*MUE2*MStop2(1)*MStop2(2)*UUStop(2,1) + 
     &   MUEC**2*MStop2(3)**2*UUStop(2,2)

        Ce(1578) = MUE2*Cd(5519)*UUStop(1,1) - 
     &   4*MUEC**2*MStop2(1)*MStop2(2)*UUStop(1,2)

        Ce(1579) = MStop2(3)**2*UUStop(1,3)*UUStop(2,2) - 
     &   4*MStop2(1)*MStop2(2)*UUStop(1,2)*UUStop(2,3)

	Ce(1580) = UUStop(1,2)*UUStop(2,1) + 2*UUStop(3,1)*UUStop(3,2)

	Ce(1583) = 8*MTy*Cd(523)*Cr(58) - 4*At*Cd(3852)*UUStop(3,2)

        Ce(1584) = MUE**2*MStop2(3)**2*UUStop(2,1) + 
     &   2*MUE2*MStop2(1)*MStop2(2)*UUStop(2,2)

        Ce(1585) = 4*MUE**2*MStop2(1)*MStop2(2)*UUStop(1,1) - 
     &   MUE2*Cd(5519)*UUStop(1,2)

        Ce(1588) = 4*MTy*Cr(91)*MStop2(1)*MStop2(2) + 
     &   At*Cd(5945)*UUStop(3,2)

	Ce(1589) = UUStop(1,1)*UUStop(2,2) + 2*UUStop(3,1)*UUStop(3,2)

	Ce(1593) = UUStop(1,3)*UUStop(2,2) + UUStop(1,2)*UUStop(2,3)

        Ce(1597) = 8*MStop2(1)*MStop2(2)*UUStop(1,2)*
     &    UUStop(2,1) - 4*MStop2(3)**2*UUStop(1,1)*UUStop(2,2)

	Ce(1599) = 4*MTy*UStop(2,1) - At*UStop(2,2)

        Ce(1600) = 2*MStop2(1)*MStop2(2)*UStop(2,2)*
     &    UUStop(1,1) - Cd(5519)*UStop(1,2)*UUStop(3,1)

        Ce(1606) = 4*MStop2(1)*MStop2(2)*UUStop(1,3)*
     &    UUStop(2,1) - MStop2(3)**2*UUStop(1,1)*UUStop(2,3)

	Ce(1609) = 2*MTy2*Cd(523) + MUE2*MStop2(1)*MStop2(2)

        Ce(1610) = MUE**2*MStop2(3)**2*UUStop(1,1) + 
     &   2*Ce(1609)*UUStop(1,2)

	Ce(161) = 2*(MStop2(2)*MStop2(4) + Cd(869)*UStop2(1,2))

        Ce(1617) = 4*MStop2(3)**2*UUStop(1,2)*UUStop(2,1) - 
     &   8*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,2)

	Ce(1620) = MStop2(1)*MStop2(2)

        Ce(1621) = (MStop2(3)**2 - 2*Ce(1620)*UStop2(1,2))/
     &   (Ce(1620)*MStop2(3)**2)

	Ce(1622) = 4*Cd(9948) - 3*UStop2(1,2)

        Ce(1623) = 4*Ce(1622)*UStop2(1,1)*UUStop(1,2) - 
     &   4*UCStop(3,1)*UStop2(1,2)*UUStop(3,2)

	Ce(1624) = 6*UStop2(1,1)**2 + UStop2(1,2)**2

        Ce(1625) = 4*Ce(1624)*UUStop(1,2) + 
     &   8*UCStop(3,1)*UStop2(1,2)*UUStop(3,2)

        Ce(1626) = (Ce(1625)*MStop2(1)*MStop2(2) + 
     &     Ce(1623)*MStop2(3)**2)/
     &   (MStop2(1)*MStop2(2)*MStop2(3)**2)

        Ce(1627) = Ce(1624)*UUStop(1,2) + 
     &   2*UCStop(3,1)*UStop2(1,2)*UUStop(3,2)

	Ce(1628) = At**2*Ce(1627) + 2*MUE2*Cd(10020)*UUStop(1,1)

	Ce(1629) = -4*Ce(1628) + 4*MUE**2*UCStop(1,3)**2*UUStop(1,1)

	Ce(1633) = Cd(11106)*UCStop(3,1) + 2*UStop2(1,2)*UUStop(1,3)

	Ce(1641) = Cd(9794)**2 - MUEC**2*UUStop(1,2)

	Ce(1642) = At*Cd(9794)*UStop(1,2) - MUEC**2*UUStop(1,2)

	Ce(1645) = 2*MTy2 - MUE2

	Ce(1649) = 12*Cd(9794) + Cd(10980)*UStop2(1,2)**2

	Ce(1650) = Cd(11136) + Cd(11145)

	Ce(1652) = 4*MStop2(1) - Cd(1136)*MStop2(3)

	Ce(1654) = 8*MTy2 - MUE2*Cd(1136)

	Ce(1655) = AtC**2*Cd(9948) - MUE**2*UStop2(1,2)

	Ce(1656) = Ce(1655)*UCStop(1,3)**2 + Ce(1654)*UStop2(1,2)

	Ce(1659) = 2*MTy*Cd(9794)*UStop(2,2) - MUE2*UUStop(3,4)

	Ce(1660) = Cd(11064) - AtC*UCStop(1,3)*UStop(1,1)

        Ce(1663) = 2*MTy*Cd(11069)*UStop(2,2) + 
     &   MUE**2*UStop(1,1)*UUStop(2,1)

	Ce(1664) = 2*Ce(1663) - MUE2*Cr(91)*UStop(2,2)

	Ce(1666) = 2*Cd(11083)*MStop2(1)**2 - 6*MStop2(3)**2*UUStop(1,1)

	Ce(1667) = Ce(1666) - 5*MStop2(3)**2*UUStop(1,1)

	Ce(1668) = Cd(11097) - 5*MTy*MStop2(3)**2*UUStop(1,3)

	Ce(1671) = 3 - Cr(75)**2

	Ce(1672) = -(Ce(1671)*MStop2(1)**2) + 3*Cr(26)*MStop2(2)

	Ce(1673) = 3*MStop2(3)**2 - 4*MStop2(1)**2*UStop2(1,1)**2

	Ce(1676) = 12 - UStop2(1,2)**2

        Ce(1683) = 2*Cd(10956)*MStop2(1)*UStop(1,2) - 
     &   2*MStop2(3)*UStop2(1,2)*UUStop(1,2)

	Ce(1684) = At*UStop(1,2)*UStop2(1,1) - Cd(698)*UStopC(2,1)

        Ce(1685) = 2*Ce(1684)*MStop2(1) - 
     &   2*At*MStop2(3)*UStop(1,2)*UStop2(1,2)

        Ce(1686) = AtC*Cd(9794)*UCStop(1,3)*UStop(1,1) - 
     &   MUE2*Cd(9948)*UUStop(1,1)

	Ce(1691) = Cd(11107)*MStop2(1) + 2*Cd(9948)*MStop2(3)*UStop(1,1)

        Ce(1692) = Cd(9809)*UCStop(3,3) - 
     &   Cd(11139)*UStop2(1,2)*UUStopC(3,2)

	Ce(1694) = 5*MStop2(3)*UStop2(1,2) - 2*MStop2(1)*UStop2(1,2)**3

	Ce(1695) = UStop2(1,1) - 3*UStop2(1,2)

	Ce(1696) = Ce(1695)*MStop2(1) - 2*MStop2(3)*UStop2(1,2)

	Ce(1697) = 4*MTy2*Ce(1694) + MUE**2*Ce(1696)*UCStop(1,3)**2

	Ce(1701) = Cd(9948) - 2*UStop2(1,1)

	Ce(1707) = Cd(11181)*UStopC(2,2) + AtC*UCStop(1,3)*UUStop(1,1)

        Ce(1715) = 2*AtC*Cd(9794)*UCStop(1,3)*UStop(1,1) - 
     &   2*MUE2*Cd(9948)*UUStop(1,1)

	Ce(1721) = Cd(970)*MStop2(3) + Cd(12114)*MStop2(1)*UStop2(1,2)

	Ce(172) = -2*Cr(85) + 6*UStop2(1,1)

	Ce(1723) = Cd(9948)*MStop2(3) + Cd(1272)*UStop2(1,1)

	Ce(1724) = Ce(1723)*MStop2(1) - Cd(951)*UStop2(1,2)

	Ce(1726) = 5*UStop2(1,1)**2 - UStop2(1,2)**2

	Ce(173) = Ce(172)*MStop2(1) + 4*Cr(75)*MStop2(2)

	Ce(1736) = 2*Cd(10960)*MStop2(1)**2 - Cd(10959)*UUStop(1,2)

	Ce(1737) = 3*MStop2(3)**2 - MStop2(1)**2*UStop2(1,1)**2

	Ce(1738) = AtC*UUStop(1,1) + 2*MTy*UUStop(1,3)

	Ce(1739) = Ce(1738)*UStop2(1,1) + 2*MTy*UCStop(3,1)*UUStop(3,4)

        Ce(1740) = Ce(1739)*MStop2(1)**2*UStop2(1,1) + 
     &   6*MTy*MStop2(3)**2*UUStop(1,3)

	Ce(1761) = UStop(1,1) + UCStop(3,1)*UStop(2,1)*UStop2(1,1)

	Ce(1766) = 2*MTy*UStop2(1,1) + Cd(9809)*UStopC(2,1)

	Ce(1770) = 4 + 3*UStop2(1,2)**2

	Ce(1772) = 12 + 5*UStop2(1,2)**2

	Ce(1774) = 12 + 13*UStop2(1,2)**2

	Ce(1781) = 14 - 3*UStop2(1,2)**2

	Ce(1783) = 12 + 7*UStop2(1,2)**2

	Ce(1787) = 4 - UStop2(1,2)**2

	Ce(18) = At*Cd(548) + 2*MTy*Cd(535)*UCStop(1,3)

	Ce(1795) = Cd(698) - MUEC*Cr(91)

	Ce(1796) = Cd(698) + MUEC*Cr(91)

	Ce(1797) = 2*MUE2*UUStop(1,1) + At**2*UUStop(1,2)

	Ce(1798) = 4*UUStop(1,3)*UUStop(2,2) - 4*UUStop(1,2)*UUStop(2,3)

	Ce(1799) = UUStop(1,2)*UUStop(2,1) - UUStop(3,1)*UUStop(3,2)

	Ce(1801) = AtC*Cd(698) - MUE2*Cr(91)

	Ce(1804) = 4*MTy2*Cr(58) - MUE2*Cr(91)

	Ce(1808) = 4*MTy2 - 3*MUE2

	Ce(1811) = 3*UStop(1,1) + UCStop(3,1)*UStop(2,1)*UStop2(1,1)

        Ce(1812) = 2*Ce(1811)*MStop2(3) - 
     &   2*Cr(91)*MStop2(2)*UStop2(1,2)*UStopC(2,2)

        Ce(1813) = MStop2(3)*UCStopC(3,1)*UUStop(1,1) - 
     &   2*MStop2(2)*UStop2(1,1)*UUStop(3,1)

        Ce(1814) = Cd(11316)*MStop2(3) - 
     &   2*MStop2(2)*UStop2(1,1)*UStop2(1,2)

        Ce(1815) = Cr(58)*MStop2(2)*UStop2(1,1)*UStop2(1,2) + 
     &   3*MStop2(3)*UUStop(3,4)

	Ce(1816) = 2*MTy2*Ce(1815) + At*MTy*Ce(1814)*UUStop(3,2)

	Ce(1817) = 2*Ce(1816) + MUE**2*Ce(1813)*UCStop(1,3)

	Ce(182) = Cd(1270)*MStop2(1) + 8*MStop2(2)*UStop2(1,1)

	Ce(1832) = 9 - UStop2(1,2)**2

        Ce(1838) = -(At**2*UCStopC(1,3)**2) + 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2)

	Ce(1850) = Cd(11335) - Cd(11334)*MStop2(3)

	Ce(185) = MStop2(2)*UStop2(1,1) - 2*MStop2(1)*UStop2(1,2)

        Ce(1856) = 3*MStop2(3)*UCStopC(3,1)*UUStop(1,1) - 
     &   4*MStop2(2)*UStop2(1,1)*UUStop(3,1)

	Ce(1857) = 5 - 3*UStop2(1,1)**2

        Ce(1858) = Ce(1857)*MStop2(3) - 
     &   4*MStop2(2)*UStop2(1,1)*UStop2(1,2)

	Ce(1859) = 2*Cr(58)*UStop2(1,1)*UStop2(1,2) + 5*UUStop(3,4)

	Ce(1860) = -4*Ce(1859)*MStop2(2) + 20*MStop2(1)*UUStop(3,4)

	Ce(1868) = 11 + 6*UStop2(1,2)**2

	Ce(1870) = 11 + 3*UStop2(1,2)**2

	Ce(186) = At - MUEC

	Ce(1876) = 11 + 4*UStop2(1,2)**2

	Ce(1878) = 11 + 7*UStop2(1,2)**2

	Ce(1880) = 16 + UStop2(1,2)**2

	Ce(187) = At + MUEC

	Ce(1886) = 22 + UStop2(1,2)**2

	Ce(19) = 2*MTy*Cr(1)*UStop(1,1) - At*Cd(595)*UStop(1,2)

	Ce(190) = 4*UStop2(1,1)**3 - 2*UStop2(1,2)**3

	Ce(191) = Ce(190)*MStop2(2) + Cd(598)*UStop2(1,1)*UStop2(1,2)

	Ce(197) = Cd(787) - UStop2(1,1)

	Ce(201) = At*AtC*Cd(520) - 4*MUE2*UStop2(1,2)

	Ce(205) = Cd(596)*MStop2(2) - Cd(1595)*UStop2(1,1)

	Ce(206) = 4*MStop2(2)*UStop2(1,1)**2 + Ce(205)*UStop2(1,2)

        Ce(207) = -(MUEC*Cd(655)*UCStopC(1,3)**2) + 
     &   MUE*Cd(649)*MStop2(3)*UStop2(1,2)

	Ce(21) = 4*MTy*MStop2(1)*UStop(1,1) - At*Cd(535)*UStop(1,2)

        Ce(208) = -2*MUE*Cd(651)*UCStop(1,3)**2 + 
     &   MUEC*Ce(206)*UStop2(1,2)

	Ce(209) = AtC*Ce(208) - 4*MTy*Cd(658)*UStop2(1,2)

	Ce(210) = 2*At*Ce(207) + Ce(209)*MStop2(3)

	Ce(211) = MStop2(3) - Cd(814)*UStop2(1,2)

	Ce(212) = 4*MStop2(2) - MStop2(3)*UStop2(1,2)

	Ce(213) = 4*UStop2(1,1)**3 - 2*Cr(85)*UStop2(1,2)**2

	Ce(215) = Cd(519)*MStop2(2) - 4*MStop2(1)*UStop2(1,2)

	Ce(216) = Ce(211)*UStop2(1,1) + 2*Ce(215)*UStop2(1,2)

	Ce(217) = UStop2(1,1)**3 - 4*UStop2(1,2)**3

	Ce(218) = Ce(217)*MStop2(3) - Ce(216)*UStop2(1,1)

	Ce(219) = Cd(787)*MStop2(1) - 2*MStop2(2)

	Ce(220) = Ce(219) + Cr(1)*UStop2(1,1)

	Ce(221) = Cd(519)*MStop2(1) - 2*Cr(85)*MStop2(2)

	Ce(222) = Ce(220)*UStop2(1,1) + Ce(221)*UStop2(1,2)

	Ce(223) = MTy*Ce(218) - AtC*Ce(222)*UCStop(1,3)

	Ce(224) = MStop2(3) + Cd(1160)*UStop2(1,2)

	Ce(225) = Ce(224) - MStop2(3)*UStop2(1,1)

	Ce(226) = Cd(508)*MStop2(2) + 2*MStop2(1)*UStop2(1,2)

	Ce(227) = Ce(225)*UStop2(1,1) - 2*Ce(226)*UStop2(1,2)

	Ce(228) = MTy*Ce(223) + MUE2*Cd(649)*UStop2(1,2)

	Ce(229) = -(MUE**2*Cd(651)) + AtC**2*Cd(665)

	Ce(233) = 4 - 8*UStop2(1,1)*UStop2(1,2)

	Ce(234) = Cd(427) + At*AtC*Cr(28)

	Ce(238) = (Cd(587) - 4*MTy*Cd(13)*MStop2(3))/MStop2(3)

	Ce(239) = Cr(115)*MStop2(2) - MStop2(1)*UStop2(1,1)*UStop2(1,2)

	Ce(240) = Cr(1) - MStop2(2)

        Ce(242) = (AtC*Ce(240)*UCStop(1,3) + 
     &     At*MStop2(3)*UCStopC(1,3))/MStop2(3)

	Ce(25) = At*MStop2(4)*UUStop(1,2) + 2*MTy*MStop2(2)*UUStop(1,3)

	Ce(248) = C2A*MUE*UStop(1,1) - 2*MTy*S2A*UStop(1,2)

	Ce(249) = 2*Ce(248)*UStop(1,1) - AtC*S2A*UUStop(1,1)

	Ce(253) = Cd(9794)*Cd(9809) - MUEC**2*UUStop(3,2)

	Ce(254) = C2A*MUEC*Cd(10987) - 1/2.D0*(S2A*Ce(253))

	Ce(255) = SA*(MUE**2*UUStop(3,1) - 4*MTy2*UUStop(3,2))

	Ce(256) = MTy*S2A*Cr(91) - C2A*MUE*UUStop(3,1)

	Ce(26) = At*MStop2(2)*UStop(1,2) + Ce(25)*UStopC(1,2)

	Ce(258) = 1/2.D0*(MUE*S2A*Cr(91)) + 2*C2A*MTy*UUStop(3,2)

	Ce(259) = -(S2A*Cd(698)) + C2A*MUEC*Cr(91)

	Ce(262) = S2A*(Cd(9794)**2 - MUEC**2*UUStop(1,2))

	Ce(263) = -(1/2.D0*Ce(262)) + C2A*MUEC*Cd(9794)*UStop(1,2)

	Ce(264) = 1/2.D0*(MUE*S2A*UStop(1,1)) + C2A*MTy*UStop(1,2)

	Ce(265) = S2A*Cd(9794) - C2A*MUEC*UStop(1,2)

	Ce(27) = At*MStop2(1)*UUStop(1,2) - 4*MTy*MStop2(2)*UUStop(1,3)

	Ce(270) = S2A*(Cd(9794)*Cd(9809) - MUEC**2*UUStop(3,2))

	Ce(271) = C2A*MUEC*Cd(10987) - 1/2.D0*Ce(270)

	Ce(272) = MUE**2*UUStop(3,1) - 4*MTy2*UUStop(3,2)

	Ce(273) = 2*Ce(256) + AtC*S2A*UUStop(3,1)

	Ce(274) = S2A*Ce(272) - AtC*Ce(273)

	Ce(275) = C2A*MTy*MUE*Cr(91) + 1/2.D0*Ce(274)

	Ce(28) = MUE**2*MStop2(1)*UUStop(1,1)**2 - At*Ce(27)*UUStop(1,2)

	Ce(284) = C2A*MUE*UStop(2,1) - 2*MTy*S2A*UStop(2,2)

	Ce(285) = MUE**2*UUStop(2,1) - 4*MTy2*UUStop(2,2)

	Ce(286) = 2*Ce(284)*UStop(2,1) - AtC*S2A*UUStop(2,1)

	Ce(287) = S2A*Ce(285) + AtC*Ce(286)

        Ce(29) = Ce(28)*UStopC(1,2) - 
     &   4*At*MTy*MStop2(2)*UStop(1,1)*UUStop(1,2)

	Ce(288) = 1/2.D0*(Ce(287) + 4*C2A*MTy*MUE*UUStop(2,3))

        Ce(289) = S2A*MStop2(2)*MStop2(3)*UStop(1,2) + 
     &   2*Ce(288)*MStop2(1)*UStopC(1,2)

	Ce(291) = S2A*(Cd(9809)**2 - MUEC**2*UUStop(2,2))

	Ce(292) = -(1/2.D0*Ce(291)) + C2A*MUEC*Cd(9809)*UStop(2,2)

	Ce(293) = 1/2.D0*(MUE*S2A*UStop(2,1)) + C2A*MTy*UStop(2,2)

	Ce(294) = S2A*Cd(9809) - C2A*MUEC*UStop(2,2)

	Ce(30) = 2*Cd(535)*UStop2(1,1) - MStop2(1)*UStop2(1,2)

	Ce(327) = 2*S2A*MStop2(2)*UStop(2,2) - 2*Ce(288)*UStopC(2,2)

	Ce(328) = At*UUStop(2,2) + 2*MTy*UUStop(2,3)

	Ce(329) = C2A*MUEC*Ce(328) - 1/2.D0*Ce(291)

	Ce(334) = C2A*MUE*UUStop(2,1) - 2*MTy*S2A*UUStop(2,3)

	Ce(337) = S2A*Cd(9809)*UStop(2,1) - C2A*MUEC*UUStop(2,3)

	Ce(338) = 1/2.D0*(MUE*S2A*UUStop(2,3)) + C2A*MTy*UUStop(2,2)

	Ce(35) = AtC - MUE

	Ce(352) = C2A*MUEC*Cd(10981) - 1/2.D0*Ce(262)

	Ce(356) = C2A*MUE*UUStop(1,1) - 2*MTy*S2A*UUStop(1,3)

	Ce(36) = AtC + MUE

	Ce(359) = S2A*Cd(9794)*UStop(1,1) - C2A*MUEC*UUStop(1,3)

	Ce(360) = 1/2.D0*(MUE*S2A*UUStop(1,3)) + C2A*MTy*UUStop(1,2)

        Ce(368) = Cd(50)*Cd(9947)*Cd(10077)*MGlpTmSt2(1) - 
     &   Cd(9949)*Cd(9950)*Cd(9969)*MGlpTmSt2(2)

	Ce(369) = At*Cd(9971) + 2*AtC*Ce(368)*UCStop(1,3)**2

	Ce(378) = AtC*Cd(9971) + 2*At*Cd(10128)*UCStopC(1,3)**2

	Ce(381) = -T134(MGl2,MTy2,Cd(1807)) + T134(MGl2,MTy2,MStop2(2))

	Ce(387) = A0delGl + Cr(13)*Cr(65)

	Ce(388) = MGl2*A0delStop(1) + Ce(387)*MStop2(1)

	Ce(389) = MGl2*Cr(13) + Cr(65)

	Ce(390) = -2*MGl2*Cd(10230) + Cd(10434)*Ce(389)

	Ce(391) = Cd(10571)*Ce(388) - 2*Ce(390)*MStop2(1)

	Ce(392) = MTy2*Cr(13) + Cr(66)

	Ce(393) = A0delT + Cr(13)*Cr(66)

	Ce(394) = MTy2*A0delStop(1) + Ce(393)*MStop2(1)

	Ce(395) = 2*MTy2*Cd(12334) + Cd(12462)*Ce(392)

	Ce(396) = Cd(11823)*Ce(394) + 2*Ce(395)*MStop2(1)

	Ce(397) = A0delGl + Cr(34)*Cr(65)

	Ce(40) = 2*MTy*MStop2(2)*UStop(1,1) - At*MStop2(1)*UStop(1,2)

	Ce(398) = MGl2*A0delStop(2) + Ce(397)*MStop2(2)

	Ce(399) = MGl2*Cr(34) + Cr(65)

	Ce(401) = MTy2*Cr(34) + Cr(66)

	Ce(402) = A0delT + Cr(34)*Cr(66)

	Ce(403) = MTy2*A0delStop(2) + Ce(402)*MStop2(2)

	Ce(404) = 2*MTy2*Cd(12365) + Cd(12570)*Ce(401)

	Ce(405) = Cd(11940)*Ce(403) + 2*Ce(404)*MStop2(2)

	Ce(406) = Cd(22)*Cd(11197) - 4*Cd(11246)*Cr(13)*MStop2(1)

	Ce(407) = 4*Cd(11053) + Ce(406)*MStop2(1)

	Ce(408) = -MGl2 + Cr(65)

	Ce(409) = MTy2*Ce(408) + MGl2*Cr(66)

	Ce(410) = -A0delGl + MGl2*Cr(14)

	Ce(413) = Cd(11609) + Cd(11615)

        Ce(414) = 2*Ce(413)*MStop2(2) - 
     &   Cd(12)*Cd(13)*Cd(11607)*MStop2(3)

	Ce(415) = Cd(11548) - Cd(11533)*MStop2(3)

	Ce(42) = 2*MTy*MStop2(3)*UStop(1,1) + At*Cd(315)*UStop(1,2)

	Ce(422) = At*Cd(10071) + 2*AtC*Cd(10128)*UCStop(1,3)**2

        Ce(424) = -(Cd(9976)*T134(MGl2,MTy2,MStop2(1))) + 
     &   Cd(9968)*T134(MGl2,MTy2,MStop2(2))

	Ce(430) = At*Cd(10080) + 2*AtC*Cd(10081)*UCStop(1,3)**2

	Ce(432) = AtC*Cd(10080) + 2*At*Cd(10125)*UCStopC(1,3)**2

        Ce(435) = AtC*Cd(10071) + 
     &   2*At*Cd(10076)*Cd(10077)*UCStopC(1,3)**2

        Ce(437) = Cd(50)*Cd(9947)*Cd(10378)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(10039)*MGlpTmSt2(2)

	Ce(438) = AtC*Cd(10118) + 2*Ce(437)

        Ce(439) = -(Cd(10006)*T134(MGl2,MTy2,MStop2(1))) + 
     &   Cd(10002)*T134(MGl2,MTy2,MStop2(2))

	Ce(443) = At*Cd(10093)*Cr(5) + AtC*Cd(10089)*UCStop(1,3)

	Ce(448) = AtC*Cd(10093)*Cr(5) + At*Cd(10097)*UCStopC(1,3)

	Ce(457) = 4*M_3C - AtC*UStop2(1,1)

	Ce(458) = 2*M_3*Cd(9983) + AtC*UCStop(1,3)**2

	Ce(47) = At*MStop2(2)*UUStop(1,2) + 2*MTy*Cr(1)*UUStop(1,3)

	Ce(469) = At - 4*M_3*Cr(85)

	Ce(477) = 2*Cr(75)**2*MStop2(1) + Cd(11624)*MStop2(2)

	Ce(48) = At*Ce(47) + 4*MTy2*MStop2(3)*UUStop(1,1)

	Ce(478) = Cd(9948)*MStop2(2) + 4*MStop2(3)*UStop2(1,1)

	Ce(480) = Cd(9948)*MStop2(2) - Cd(256)*UStop2(1,2)

	Ce(481) = 2*Cr(85)*MStop2(3)*UStop2(1,1) + Ce(480)*UStop2(1,2)

	Ce(484) = -(Cd(12125)*MStop2(2)) + 4*MStop2(3)*UStop2(1,1)

        Ce(485) = Cd(11638)*MStop2(1) - 
     &   Ce(484)*MStop2(3)*UStop2(1,1)*UStop2(1,2)**2

	Ce(49) = At**2*MStop2(1)*UStop(1,2) + Ce(48)*UStopC(1,2)

	Ce(488) = Cd(11640)*MStop2(2) - 2*Cr(75)**2*MStop2(3)

        Ce(490) = MUEC*Cd(537)*UStop2(1,1) + 
     &   MUE*UCStop(1,3)**2*UStop2(1,2)

	Ce(492) = Cd(9948) - 2*UStop2(1,2)

        Ce(495) = MUE*Cd(537)*UStop2(1,1) + 
     &   MUEC*UCStopC(1,3)**2*UStop2(1,2)

	Ce(5) = MTy*Cd(492) + AtC*UCStop(1,3)

	Ce(50) = 4*MTy2*MStop2(3)*UUStop(1,1) + At*MTy*Cr(1)*UUStop(1,3)

	Ce(501) = Cd(47)**2 - Cd(12298)*MStop2(1)

	Ce(502) = Cd(164)*Cd(2597)*Cd(10227) + MTy*Cd(13)*Ce(501)

        Ce(51) = Ce(50)*UStop(1,2)*UStop2(1,2) - 
     &   At*MTy*MStop2(1)*UStop(1,1)*UUStop(1,2)

        Ce(511) = MGl2*Cd(10058) + 
     &   4*M_3C*MTy*MUEC*UCStop(1,3)*UStop2(1,1)

	Ce(513) = 2*Cd(9955)*Cd(10425) + Cd(10563)*MStop2(3)

	Ce(514) = 2*M_3C*MTy2*Cd(10054) + MTy*Cd(13)*Ce(513)

	Ce(515) = -(MTy2*Cd(10055)) + AtC*MTy*Cd(13)*Cd(10091)

	Ce(517) = -2*MUE*Cr(85)*UCStop(1,3)**2 + MUEC*Cr(28)*UStop2(1,1)

        Ce(52) = MUE**2*MStop2(4)*UUStop(1,1) - 
     &   4*MTy2*MStop2(3)*UUStop(1,2)

	Ce(520) = 2*M_3C*MTy2*Cd(10054) - MTy*Cd(13)*Cd(12406)

	Ce(521) = -2*MTy2*Cd(10055) + AtC*MTy*Cd(13)*Cd(10091)

        Ce(530) = -2*MUEC*Cr(85)*UCStopC(1,3)**2 + 
     &   MUE*Cr(28)*UStop2(1,1)

	Ce(531) = MUE*Cr(28) - 3*MUEC*UCStopC(1,3)**2

	Ce(532) = M_3C*Cd(13)*Cr(28) + M_3*Ce(531)*UCStopC(1,3)

        Ce(549) = 2*Cd(47)*Cd(164)*Cd(10058)*MStop2(1) - 
     &   MUEC*Cd(12371)*MStop2(3)*UStop2(1,1)**2

        Ce(551) = 2*Cd(47)*Cd(164)*Cd(10059)*MStop2(1) - 
     &   MUE*Cd(12371)*MStop2(3)*UStop2(1,1)**2

        Ce(554) = -(Cd(10398)*UCStop(1,3)**2) + 
     &   4*M_3*UStop2(1,1)*UStop2(1,2)

	Ce(555) = 2*M_3*Cr(28) + Cd(10398)*UCStop(1,3)**2

	Ce(556) = -8*M_3C + AtC*Cd(9972)

	Ce(557) = 2*M_3*Cd(11742) + Ce(556)*UCStop(1,3)**2

	Ce(558) = MGlmT2*Ce(555) + Ce(557)*MStop2(1)

        Ce(560) = -(Cd(10140)*UCStopC(1,3)**2) + 
     &   4*M_3C*UStop2(1,1)*UStop2(1,2)

	Ce(561) = 2*M_3C*Cr(28) + Cd(10140)*UCStopC(1,3)**2

	Ce(562) = 2*M_3C*Cd(11742) + Cd(11916)*UCStopC(1,3)**2

	Ce(563) = MGlmT2*Ce(561) + Ce(562)*MStop2(1)

        Ce(575) = 4*MGlpTmSt2(1)*MStop2(1)**2 - 
     &   Cd(152)*Cd(2597)*UStop2(1,1)

	Ce(580) = 2*MGlmT2*Cd(47)*Cd(9962) + At*Cd(10280)*UStop2(1,1)

	Ce(581) = 2*MGlmT2*Cd(47)*Cd(9963) + AtC*Cd(10280)

	Ce(585) = 2*Cd(10280)*UStop2(1,1) + Cd(10279)*UStop2(1,2)

        Ce(586) = -2*Cd(47)*Cd(10443)*Cr(28)*MStop2(1) + 
     &   Ce(585)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Ce(587) = MGlmT2 + Cd(10077)*MStop2(1)

	Ce(59) = 2*MTy*UCStop(1,3) + At*UStop2(1,2)

	Ce(588) = 4*Cd(47)*Ce(587)*MStop2(1) + Cd(10279)*MStop2(3)

	Ce(589) = At*Ce(586) + AtC*Ce(588)*UCStop(1,3)**2

	Ce(590) = AtC*Ce(586) + At*Ce(588)*UCStopC(1,3)**2

	Ce(591) = MUE*Ce(589) + MUEC*Ce(590)

	Ce(592) = 2*MGlmT2*Cd(47)*Cd(9994) + AtC*Cd(10280)*UStop2(1,1)

	Ce(593) = 2*MGlmT2*Cd(47)*Cd(9995) + At*Cd(10280)

	Ce(6) = Ce(5) + At*UCStopC(1,3)

	Ce(60) = 2*MTy*UStop(1,1) + 3*At*UStop(1,2)

        Ce(599) = -2*Cd(47)*Cd(152)*Cd(10059)*MStop2(1) + 
     &   MUE*Cd(10279)*MStop2(3)*UStop2(1,1)**2

	Ce(601) = MStop2(1) + MGlmT2*UStop2(1,1)

	Ce(602) = MUEC*Ce(601)*Cr(28) - 2*MUE*Cd(11792)*UCStop(1,3)**2

	Ce(603) = MUE*Ce(601)*Cr(28) - 2*MUEC*Cd(11792)*UCStopC(1,3)**2

	Ce(604) = AtC*Ce(517) + At*Ce(530)

	Ce(609) = At*Cd(10096) + Cd(9955)*Cd(10500)

	Ce(610) = 2*Ce(609) - Cd(10413)*MStop2(3)

	Ce(611) = 2*AtC*Cd(10091)*Cr(5) + Ce(610)*UCStopC(1,3)

	Ce(616) = Cd(47) + Cd(10383)*MStop2(1)

	Ce(619) = Cd(1281)*UStop2(1,1)**2 + MStop2(1)*UStop2(1,2)**2

	Ce(63) = 2*MTy2*UUStop(1,1) + At*MTy*UUStop(1,3)

	Ce(630) = At*Cd(10443)*Cr(28) - 2*AtC*Cd(10445)*UCStop(1,3)**2

	Ce(631) = AtC*Cd(10443)*Cr(28) - 2*At*Cd(10445)*UCStopC(1,3)**2

	Ce(632) = MUE*Ce(630) + MUEC*Ce(631)

	Ce(633) = MGl2*Cd(795) - Cd(10382)*MStop2(1)

	Ce(634) = M_3C*MUE*UCStop(1,3)**2 + M_3*MUEC*UCStopC(1,3)**2

	Ce(635) = MGl2*Cd(10219) - Cd(10388)*MStop2(1)

        Ce(637) = Cd(47)*Cd(152)*Cd(10059) + 
     &   2*MUE*Cd(164)*MStop2(1)*MStop2(3)*UStop2(1,1)**2

	Ce(64) = At*UUStop(1,2) + 4*MTy*UUStop(1,3)

	Ce(640) = MUE*Cd(10555) + MUEC*Cd(10557)

	Ce(642) = MUEC*Cr(28) - 3*MUE*UCStop(1,3)**2

	Ce(646) = 2*Cd(12452) + Cd(10563)*MStop2(3)*UCStop(1,3)

	Ce(647) = 2*Cd(12456) + Cd(10563)*MStop2(3)

	Ce(65) = At*Ce(64) + 4*MTy2*UUStop(1,1)

	Ce(648) = 2*AtC*Cd(10018)*Cr(5) + Ce(647)*UCStopC(1,3)

	Ce(649) = MUE*Ce(646) + MUEC*Ce(648)

	Ce(66) = MUE**2*UUStop(1,1) - 4*MTy2*UUStop(1,2)

	Ce(660) = -(AtC*Cd(11942)) + M_3C*MStop2(3)

        Ce(661) = Ce(660)*UCStop(1,3)**2 + 
     &   M_3*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Ce(664) = At*Cd(10077)*MStop2(1) + M_3*MStop2(3)

        Ce(665) = Ce(664)*UCStopC(1,3)**2 + 
     &   M_3C*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Ce(67) = Ce(65)*UStop2(1,1) + 4*Ce(63)*UStop2(1,2)

	Ce(669) = 2*MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)**2

	Ce(670) = MUEC*Cd(10658) - 2*MUE*UCStop(1,3)**2*UStop2(1,1)

	Ce(672) = MUE*Cd(10658) - 2*MUEC*UCStopC(1,3)**2

	Ce(675) = MGlmT2*Cd(9983) - MStop2(2)*UStop2(1,1)**2

	Ce(68) = Ce(67)*UUStopC(1,1) - Ce(66)*UStop2(1,1)*UUStopC(1,2)

	Ce(679) = MGlmT2 + 2*MStop2(2)

	Ce(680) = MGlmT2 + MStop2(2)*UStop2(1,1)

	Ce(681) = -4*M_3C*Ce(679) + AtC*Ce(680)

	Ce(687) = MGlmT2*Cr(85) + Cr(54)*MStop2(2)

        Ce(69) = MUEC**2*UCStopC(1,3)**2 + 
     &   2*MUE2*UStop2(1,1)*UStop2(1,2)

	Ce(689) = MGlmT2*Cd(10861) + 2*Cd(10181)*MStop2(2)

        Ce(693) = -(MGlmT2*Cd(50)*Cd(9955)*Cr(75)) + 
     &   MGl2*Cd(179)*MStop2(3)

	Ce(694) = -(Cd(9955)*Cr(75)*MStop2(1)) + MStop2(3)

	Ce(695) = 2*Ce(694)*MStop2(2) + MGlmT2*Cd(10029)*MStop2(3)

	Ce(697) = MStop2(1)*MStop2(2) + MGlmT2*MStop2(3)

	Ce(708) = At*Cd(10070) - 2*AtC*UCStop(1,3)**2

	Ce(709) = 4*Cd(138)*Cd(10169) + Ce(708)*MStop2(1)

	Ce(710) = At*Cd(11853) - 2*AtC*Cd(11854)*UCStop(1,3)**2

        Ce(714) = 2*Cd(50)*Cd(251)*Ce(669)*MStop2(2) - 
     &   MUEC*Cd(12483)*MStop2(3)*UStop2(1,1)**2

	Ce(716) = 2*Cd(50)*Cd(251)*MStop2(2) - Cd(12483)*MStop2(3)

	Ce(719) = AtC*Cd(10070) + 2*At*Cd(10077)*UCStopC(1,3)**2

	Ce(720) = 4*Cd(138)*Cd(10176) + Ce(719)*MStop2(1)

	Ce(721) = AtC*Cd(11853) - 2*At*Cd(11855)*UCStopC(1,3)**2

	Ce(727) = -8*M_3C + AtC*Cd(9969)

	Ce(73) = C2A*MUEC - At*S2A

	Ce(728) = 2*M_3*Cd(11878) + Ce(727)*UCStop(1,3)**2

	Ce(729) = At*Cd(9970) + 2*Ce(728)

	Ce(732) = 2*M_3C*Cd(11878) + Cd(12438)*UCStopC(1,3)**2

	Ce(733) = AtC*Cd(9970) + 2*Ce(732)

	Ce(74) = (C2A*MUE - AtC*S2A)*MStop2(4)

        Ce(753) = 4*MGlpTmSt2(2)*MStop2(2)**2 - 
     &   Cd(179)*Cd(2681)*UStop2(1,1)

	Ce(76) = 4*MTy*MStop2(3) + 2*AtC*MStop2(4)*UCStop(1,3)

	Ce(758) = MUE*Cd(577) - MUEC*Cd(843)*UCStopC(1,3)**2

	Ce(759) = MTy2*Cd(10603) - 2*Cd(10604)*MStop2(2)

	Ce(778) = MUEC*Cd(577) - MUE*Cd(843)*UCStop(1,3)**2

	Ce(782) = MGlmT2*Cd(44) - Cd(10705)*MStop2(2)

	Ce(783) = 2*M_3*MTy*Cr(75)*MStop2(2) - Ce(782)*UCStop(1,3)

	Ce(79) = C2A*MUE - AtC*S2A

        Ce(794) = -(MUE*Cd(179)*UCStop(1,3)**2) + 
     &   MUEC*MStop2(2)*UStop2(1,1)*UStop2(1,2)

	Ce(797) = MStop2(2) + MGlmT2*UStop2(1,1)

	Ce(800) = Cd(10413)*MStop2(3) - 2*Cr(75)*MStop2(1)*UStop2(1,2)

	Ce(801) = Cd(10025)*MStop2(3) + 2*Cd(118)*Cr(75)*UStop2(1,2)

	Ce(802) = -(MGlmT2*Ce(801)) + Ce(800)*MStop2(2)

        Ce(803) = -(MUEC*Cd(179)*UCStopC(1,3)**2) + 
     &   MUE*MStop2(2)*UStop2(1,1)*UStop2(1,2)

	Ce(807) = Cd(10602)*MStop2(2) + Cd(10604)*UStop2(1,1)

        Ce(808) = -4*MTy2*MUE*UCStop(1,3)**2 + 
     &   MUEC*Cd(10702)*UStop2(1,1)

	Ce(820) = 2*MTy2*Cd(10388) + MUEC*Cd(179)*UCStopC(1,3)**2

	Ce(822) = Cd(10563)*MStop2(3) + 2*Cd(118)*Cr(75)*UStop2(1,2)

	Ce(823) = MGlmT2*Ce(822) - Ce(800)*MStop2(2)

	Ce(825) = Cd(179) + 2*MTy2*UStop2(1,1)

	Ce(831) = Cd(11163)*MStop2(1) + 2*Cd(11159)*MStop2(3)

	Ce(832) = Cd(11030) - Cd(11034)

	Ce(834) = Cd(11238) - 2*Ce(831)*UStop2(1,1)

        Ce(837) = 4*Cd(10988)*MStop2(1) - 
     &   Cd(10990)*MStop2(3)*UStop2(1,2)*UUStopC(1,2)

        Ce(838) = Ce(837)*MStop2(3)*UStopC(1,1) + 
     &   2*Cd(11008)*UUStopC(1,3)

	Ce(848) = -2*AtC*Cr(91)*UCStop(1,3) + 2*Cd(10987)*UStop2(1,2)

	Ce(849) = Cd(50)*Cd(152) - Cd(47)*Cd(12180)

        Ce(850) = 2*MTy2*Cd(47)*UStop2(1,1)*UUStop(3,3) - 
     &   Ce(849)*UUStop(3,4)

	Ce(851) = (MTy*Cd(47)*Ce(848) + Ce(850))/(Cd(47)*Cd(50))

        Ce(852) = (2*Cd(12178)*UStopC(1,1) - 
     &     Cd(47)**2*Ce(851)*UStopC(2,1))/Cd(47)**2

        Ce(853) = (-(Cd(47)**2*Ce(852)) + 
     &     MTy*Cd(9794)*MGlpTmSt2(1)*MStop2(3)*UCStopC(1,3)*
     &      UStop2(1,1))/(Cd(47)**2*MStop2(3))

        Ce(862) = (Cd(50)*Cd(12249) + Cd(47)**2*Cd(12250))/
     &   (Cd(47)**2*Cd(50))

        Ce(863) = Cd(10985)*MGlpTmSt2(1)*MStop2(3)*UStop2(1,1)*
     &    UStopC(1,2) - 2*Cd(47)*Cd(12241)*UUStopC(3,4)

        Ce(864) = (MUE*Ce(863)*UStop2(1,2) + 
     &     MUEC*Cd(9794)*MGlpTmSt2(1)*MStop2(3)*UCStopC(1,3)**2*
     &      UStopC(1,1))/MStop2(3)

	Ce(866) = Cd(9948)*Cd(12269) - Cd(12267)*UStop2(1,2)

        Ce(867) = -(MUEC*Ce(866)*UCStopC(1,3)) - 
     &   MUE*Cd(12265)*UCStop(1,3)*UStop2(1,1)

	Ce(868) = 3*Cd(9948) + UStop2(1,1)

	Ce(872) = MTy2*Cd(50) + Cd(47)*MStop2(1)

	Ce(877) = MStop2(2) + 3*Cd(9948)*MStop2(3)

	Ce(878) = MStop2(2) - Cd(814)*UStop2(1,1)

	Ce(880) = Ce(877)*UStop2(1,1) - Ce(878)*UStop2(1,2)

        Ce(881) = -2*MUEC*Cr(85)*MStop2(2)*UCStopC(1,3)**2 + 
     &   MUE*Ce(880)*UStop2(1,1)

	Ce(884) = Cd(12154) + Cd(6458)*Cd(12149)*MStop2(3)*UStopC(1,1)

	Ce(891) = 8 - 3*Log(MStop2(1)/Q)

	Ce(892) = 9 - Ce(891)*Log(MStop2(1)/Q)

	Ce(893) = 6*A0delStop(1) + Ce(892)*MStop2(1)

	Ce(894) = 8 - 3*Log(MStop2(2)/Q)

	Ce(895) = 9 - Ce(894)*Log(MStop2(2)/Q)

	Ce(896) = 6*A0delStop(2) + Ce(895)*MStop2(2)

	Ce(9) = MStop2(3) - MStop2(2)*UStop2(1,1)

	Ce(90) = Cd(435)*Ce(79) - 2*Ce(73)*MStop2(4)*UCStopC(1,3)**2

        Ce(898) = Cd(9947)*Cd(10103)*MGlpTmSt2(1)*MStop2(1) + 
     &   Cd(9949)*Cd(9950)*MGlpTmSt2(2)*MStop2(2)

	Ce(900) = MGlpT2*Cd(10112) - 2*MTy*Cr(5)

	Ce(901) = MGlpT2*Cd(10112) + 2*MTy*Cr(5)

        Ce(902) = Cd(50)*Cd(9947)*Ce(900)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Ce(901)*MGlpTmSt2(2)

        Ce(904) = ((2*Cd(47)**2*Cd(9973) + 
     &       Cd(50)*Cd(9950)*Cd(9951)*MStop2(3))*UCStop(1,3)**2)/
     &   (MZ2*SW2*Cd(47)**2*Cd(50)*Cd(9950)*MStop2(3))

	Ce(906) = 2*M_3*AtC - MUE2

        Ce(91) = (Ce(90) + 4*MTy*S2A*MStop2(3)*UCStopC(1,3))/
     &   (MZ2*SW2*MStop2(3))

        Ce(908) = Cd(50)*Cd(9947)*Cd(10001)*MGlpTmSt2(1) + 
     &   2*M_3*Cd(47)*Cd(9949)*MGlpTmSt2(2)

        Ce(912) = ((-2*Cd(47)*Cd(50)*Cd(9985) + 
     &       Cd(9980)*MStop2(3))*UCStop(1,3)**2)/
     &   (MZ2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))

        Ce(916) = ((-2*Cd(47)*Cd(50)*Cd(9989) + 
     &       Cd(9980)*MStop2(3))*UCStopC(1,3)**2)/
     &   (MZ2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))

        Ce(920) = Cd(50)*Cd(9947)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(9955)*MGlpTmSt2(2)

        Ce(921) = ((2*Cd(47)**2*Cd(10003) + 
     &       Cd(50)*Cd(9950)*Cd(9951)*MStop2(3))*UCStopC(1,3)**2)/
     &   (MZ2*SW2*Cd(47)**2*Cd(50)*Cd(9950)*MStop2(3))

	Ce(924) = MUE2*Cr(28) - AtC**2*UCStop(1,3)**2

	Ce(93) = 2*C2A*Cd(426) - S2A*Ce(12)

	Ce(930) = M_3*AtC**2 + 2*M_3C*Ce(1)

        Ce(933) = (4*Cd(585)*
     &     (Cd(10139)*MStop2(3) - 
     &       2*Cd(47)*Cd(50)**2*Cd(9947)*Cd(9955)*MGlpTmSt2(1)*
     &        UStop2(1,2)))/(MZ2*SW2*Cd(47)**2*Cd(50)**2*MStop2(3))

        Ce(934) = (4*(-2*Cd(47)*Cd(585)*Cd(9955)*Cd(10009) + 
     &       Cd(50)*Cd(585)*Cd(10012)*MStop2(3)))/
     &   (MZ2*SW2*Cd(47)**2*Cd(50)*MStop2(3))

	Ce(937) = 2*MUE2*UCStop(1,3) - At**2*UCStopC(1,3)

	Ce(939) = 2*M_3 - At*Cd(10025)

	Ce(940) = 2*M_3C*UCStop(1,3) + Ce(939)*UCStopC(1,3)

        Ce(944) = 2*MUE2*Cd(10018)*UCStop(1,3) + 
     &   At**2*Cd(10091)*UCStopC(1,3)

	Ce(945) = 2*MUE2*Cd(10018) + 3*At**2*UCStopC(1,3)**2

	Ce(946) = M_3C*Ce(944) + M_3*Ce(945)*UCStopC(1,3)

        Ce(947) = 2*MUE2*Cd(10020)*UCStop(1,3) - 
     &   At**2*Cd(10092)*UCStopC(1,3)

	Ce(948) = 2*MUE2*Cd(10020) + 3*At**2*UCStopC(1,3)**2

	Ce(949) = M_3C*Ce(947) + M_3*Ce(948)*UCStopC(1,3)

	Ce(953) = Cd(10101) - 3*Cd(9951)*UStop2(1,1)**2

	Ce(960) = Cd(251)*Cd(10176) + Cd(11888)*MGlpTmSt2(2)

        Ce(961) = Cd(50)*Cd(9947)*Cd(10176)*Cd(12594) - 
     &   Cd(47)*Cd(9949)*Ce(960)

	Ce(974) = M_3*AtC - MUE2*UStop2(1,2)

	Ce(978) = Cd(47) - MGlmT2*MStop2(3)

	Ce(979) = Cd(47) + 2*MGlmT2*MStop2(3)*UStop2(1,1)

	Ce(991) = Cd(47)*Cd(152) - Cd(164)*MStop2(1)*MStop2(3)

	Ce(995) = 2*MGlmT2*Cd(47) - Cd(152)*Cd(164)*Cd(10015)

        Ce(996) = -(Cd(508)*UStop2(1,1)) + UStop2(1,1)**2 + 
     &   2*UStop2(1,2)**2

        Ce(997) = -(Cd(152)*Cd(164)*Cd(10353)) + 
     &   (Cd(9948)*Ce(995) - MGlmT2*Cd(47)*Ce(996))*UStop2(1,1)

        Ce(998) = Cd(9948)*Cd(11796) - 
     &   2*M_3C*Cr(85)*UCStop(1,3)**2 + M_3*Cr(28)*UStop2(1,1)

        Ce(1003) = Ce(1002)*UCStopC(1,3) + 
     &   2*UCStop(1,3)*(MUE2*UStop2(1,1) + 2*MTy2*UStop2(1,2))

        Ce(1010) = Ce(1009)*UCStopC(1,3) + 
     &   2*UCStop(1,3)*(2*MTy2*UStop2(1,1) + MUE2*UStop2(1,2))

        Ce(102) = UCStopC(1,3)*
     &    (2*MTy*S2A + Ce(73)*UCStopC(1,3)) + 
     &   Ce(79)*UStop2(1,1)*UStop2(1,2)

        Ce(1025) = -(MUE2*Cd(10284)*MStop2(3)*UStop2(1,1)*
     &      UStop2(1,2)) + 
     &   Cd(47)*(Ce(1024)*UStop2(1,1) - 
     &      MUE2*Cd(47)*Cd(9969)*UStop2(1,2)**2)

        Ce(1032) = -(Cd(13817)*MStop2(1)) + 
     &   Cd(9948)*MStop2(1)**2 + MGlmT2*MStop2(2)

        Ce(1033) = -(Cd(47)**2*Ce(1031)) + 
     &   MStop2(1)*(4*Cd(47)*Ce(1032) + Cd(10279)*MStop2(3))

        Ce(104) = UCStopC(1,3)*
     &    (-2*MTy*S2A + Ce(73)*UCStopC(1,3)) + 
     &   Ce(79)*UStop2(1,1)*UStop2(1,2)

        Ce(1039) = -2*Cd(47)*
     &    (-(M_3*Cd(9948)*Cd(10303)) + Ce(1038))*MStop2(1) + 
     &   Cd(10086)*Cd(10280)*MStop2(3)

        Ce(1042) = Cd(47)*Ce(1041) + AtC*Cd(10307)*UCStop(1,3) + 
     &   At*Cd(10278)*UCStopC(1,3)

        Ce(1043) = -(M_3*At**2*UCStopC(1,3)**3) + 
     &   M_3C*Ce(937)*UStop2(1,1)*UStop2(1,2) + 
     &   2*M_3*MUE2*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Ce(1046) = -2*Cd(47)*
     &    (-(M_3C*Cd(9948)*Cd(10303)) + Ce(1045))*MStop2(1) + 
     &   Cd(10039)*Cd(10280)*MStop2(3)

        Ce(1047) = 2*M_3*MUE2*Cd(10018) + 
     &   At*(3*M_3*At*UCStopC(1,3)**2 - 
     &      2*MGlmT2*Cd(9955)*UStop2(1,1))

        Ce(1051) = 2*M_3*Ce(1050) + 
     &   At*(-(Cd(10025)*MStop2(3)) + 3*M_3*At*UCStopC(1,3)**2)

	Ce(1076) = MTy2 - 2*MTy*Cr(5) + MGl2*UStop2(1,2)

        Ce(1086) = Cd(12369) + 
     &   Cd(47)*(5*Cd(118) + 2*MUE2*UStop2(1,1)**2)

	Ce(1095) = -MStop2(2) + MStop2(1)*(Cd(9948) - UStop2(1,2))

        Ce(1096) = Cd(47)**2*Ce(1031) + 
     &   MStop2(1)*(-4*Cd(47)*Cd(164)*Ce(1095) + 
     &      Cd(12371)*MStop2(3))

	Ce(1113) = 3 + 2*UStop2(1,1) + 2*UStop2(1,2)

	Ce(112) = 4*MTy + 3*AtC*UCStop(1,3) + 3*At*UCStopC(1,3)

        Ce(1116) = (Cd(10313)*Cd(12371)*MStop2(3) + 
     &     4*Cd(47)*Cd(9955)*MStop2(1)*
     &      (-(Cd(10143)*MStop2(1)) + MGlmT2*Cd(10313)*UStop2(1,1))
     &     )/(MStop2(1)**2*MStop2(3))

	Ce(1123) = 2*MGlmT2*Cd(47)*Cd(10025) + 2*Cd(12370) - Cd(12371)

        Ce(1125) = At*MGlmT2*Cr(5)*Cr(28) + 
     &   UStop2(1,1)*(-(Cd(9948)*Ce(1124)) + 
     &      Cd(47)*UCStop(1,3)*UStop2(1,2))

        Ce(1144) = (MGl2*Cd(152)*Cd(10313)*MStop2(3) + 
     &     Cd(47)*Cd(9955)*
     &      (-(Cd(10143)*MStop2(1)) + MGlmT2*Cd(10313)*UStop2(1,1))
     &     )/(MStop2(1)*MStop2(3))

        Ce(1145) = -(Ce(186)*Ce(187)*UCStopC(1,3)**3) + 
     &   MUE**2*UCStop(1,3)*UStop2(1,1)*UStop2(1,2) + 
     &   2*MUE2*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Ce(1147) = Cd(152)*MStop2(3)*
     &    (-2*Cd(10506) + Cd(152)*UCStop(1,3)) - 
     &   2*Cd(47)*Cd(164)*Cd(9955)*UCStop(1,3)*UStop2(1,2)

	Ce(1149) = Ce(944) + 2*MStop2(3)*UCStop(1,3)

        Ce(1151) = 2*M_3*Ce(1150) + 
     &   At*(-(Cd(10025)*MStop2(3)) + 3*M_3*At*UCStopC(1,3)**2)

	Ce(116) = 4*MTy - 3*(AtC*UCStop(1,3) + At*UCStopC(1,3))

        Ce(1160) = -(Cd(2429)*UStop2(1,1)) + UStop2(1,1)**3 + 
     &   UStop2(1,2)**3

        Ce(1163) = 8*M_3*Cr(28) - 16*M_3C*UCStop(1,3)**2 - 
     &   At*Cr(28)*UStop2(1,1)

        Ce(1182) = Cd(50)**2*Ce(1180) + Cd(47)**2*Ce(1181) + 
     &   M_3C*Cd(9948)*Cd(12599)*UStop2(1,1)

        Ce(1187) = Cd(47)*
     &    (Cd(47)*Ce(1186) + 2*Cd(50)**2*MStop2(2)) + 
     &   Cd(50)**2*MGlpTmSt2(1)*MStop2(1)*MStop2(3)

        Ce(1189) = Cd(47)*
     &    (Cd(47)*Ce(1188) + 2*Cd(50)**2*MStop2(2)) + 
     &   Cd(50)**2*MGlpTmSt2(1)*MStop2(1)*MStop2(3)

        Ce(1205) = Cd(47)*
     &    (2*Cd(50)**2*Cd(9972) - Cd(47)*Ce(1204)) + 
     &   Cd(50)**2*MGlpTmSt2(1)*MStop2(3)

	Ce(1212) = 5 - 24*UStop2(1,1)*UStop2(1,2)

        Ce(1224) = -(Cd(50)*Cd(10454)) - Cd(47)*Cd(10759) + 
     &   M_3*Cd(9948)*Cd(12042)

        Ce(1226) = -(Cd(50)*Cd(10321)) - Cd(47)*Cd(10764) + 
     &   M_3C*Cd(9948)*Cd(12042)

        Ce(1230) = 2*MUE2*UCStopC(1,3)*UStop2(1,1) + 
     &   AtC*(AtC*UCStop(1,3) + 2*At*UCStopC(1,3))*UStop2(1,2)

	Ce(1236) = 3 - 4*UStop2(1,1) - 4*UStop2(1,2)

        Ce(1237) = M_3C*UCStop(1,3) + AtC*UCStop(1,3) + 
     &   Cd(9995)*UCStopC(1,3)

	Ce(1254) = Cd(50) + Ce(1303)

	Ce(1255) = Cd(50) - Ce(1305)

        Ce(1263) = At*UStop2(1,1)**2 + 
     &   2*Cd(9948)*(M_3 + 2*At*UStop2(1,1))

	Ce(1265) = 2*M_3 + At*(Cd(9948) - UStop2(1,2))

        Ce(127) = UStop2(1,1)**2*
     &    (Ce(126) - MStop2(1)*UStop2(1,1)) + 
     &   2*Cd(800)*MStop2(3)*UStop2(1,2)**2

	Ce(1273) = Ce(1311) - M_3*MGlmT2*Cd(9948)*Cr(75)

	Ce(1275) = Ce(1316) - M_3C*MGlmT2*Cd(9948)*Cr(75)

        Ce(128) = (Cd(519) - 3*UStop2(1,1))*UStop2(1,1) - 
     &   2*UStop2(1,2)**2

        Ce(1282) = 12*Cd(50) + 4*Cd(179)*Cd(251) + 
     &   At*AtC*Cd(251)*UStop2(1,1)

        Ce(1289) = -(Cd(179)*Cd(251)*Cd(10353)) + 
     &   (-(MGlmT2*Cd(50)*Ce(996)) + Cd(9948)*Ce(1288))*UStop2(1,1)

        Ce(1290) = M_3*Cd(10658) - Cd(9948)*Cd(10851) - 
     &   2*M_3C*UCStop(1,3)**2*UStop2(1,1)

        Ce(1304) = Cd(50)**2*Ce(1302) + 
     &   MStop2(2)*(4*Cd(50)*Ce(1303) + Cd(10603)*MStop2(3))

        Ce(1306) = Cd(50)**2*Ce(1302) + 
     &   MStop2(2)*(-4*Cd(50)*Ce(1305) + Cd(10603)*MStop2(3))

        Ce(1312) = 2*Cd(50)*
     &    (Ce(1311) - M_3*MGlmT2*Cd(9948)*Cr(75))*MStop2(2) + 
     &   Cd(10086)*Cd(10604)*MStop2(3)

        Ce(1313) = Cd(50)*Ce(1041) + AtC*Cd(10628)*UCStop(1,3) + 
     &   At*Cd(10602)*UCStopC(1,3)

        Ce(1317) = 2*Cd(50)*
     &    (Ce(1316) - M_3C*MGlmT2*Cd(9948)*Cr(75))*MStop2(2) + 
     &   Cd(10039)*Cd(10604)*MStop2(3)

        Ce(133) = UStop2(1,1)*(Ce(132) - Cr(1)*UStop2(1,1)) - 
     &   MStop2(3)*UStop2(1,2)**2

	Ce(1338) = Ce(1157) - 6*MStop2(2)

        Ce(1347) = Cd(50)*
     &    (Cd(50)*Ce(1302) + 
     &      4*Cd(251)*Cd(9969)*MStop2(1)*MStop2(2)) - 
     &   Cd(12483)*MStop2(2)*MStop2(3)

        Ce(1348) = Cd(50)**2*Ce(1302) + 
     &   MStop2(2)*(4*Cd(50)*Cd(251)*Cd(10745) - 
     &      Cd(12483)*MStop2(3))

        Ce(1355) = -(AtC**2*UCStop(1,3)**2) + 
     &   Ce(1353)*UStop2(1,1) - Ce(1354)*UStop2(1,2)**2

        Ce(1360) = AtC*MStop2(1) - 
     &   MStop2(3)*(2*M_3C + AtC*UStop2(1,2))

        Ce(1361) = MUE2*Cd(774)*UStop2(1,1)**2 - 
     &   MStop2(3)*(AtC**2*UCStop(1,3)**2 + MUE2*UStop2(1,2)**2)

        Ce(1382) = UStop2(1,2)**3 + 
     &   UStop2(1,1)*(UStop2(1,1)**2 + 3*UStop2(1,2)**2)

        Ce(1387) = 8*M_3*Cr(85)*UCStopC(1,3)**2 - 
     &   4*M_3C*Cr(75)**2*UStop2(1,1) + AtC*Cr(75)**2*UStop2(1,1)

        Ce(1391) = 2*M_3*UCStopC(1,3)**2 + 
     &   M_3C*(3*UStop2(1,1)**2 + UStop2(1,2)**2)

        Ce(1404) = -(AtC**2*UCStop(1,3)**2) + 
     &   Ce(1402)*UStop2(1,1) - Ce(1403)*UStop2(1,2)**2

	Ce(1412) = Cd(9948)*Ce(1410) - MStop2(2) + Ce(1411)*UStop2(1,1)

        Ce(1414) = 2*M_3*AtC*UStop2(1,1) + 
     &   MUE2*UStop2(1,1)**2 + Ce(1413)*UStop2(1,2)

        Ce(1415) = At*Cd(11852) + AtC*Cd(11850)*UCStop(1,3)**2 + 
     &   Ce(1414)*UStop2(1,2)

        Ce(1424) = (Cd(9948) - UStop2(1,1))*UStop2(1,1) - 
     &   Cd(9948)*UStop2(1,2)

	Ce(1430) = MGlmT2*Cd(9948) + Ce(797)

        Ce(1462) = -2*MGl2*Ce(1461)*MStop2(3)*UUStop(1,1) + 
     &   MStop2(1)*(Cd(44)*Ce(1460) - MGl2*Cd(46)*UUStop(2,1))

        Ce(147) = 2*Ce(146)*UStop2(1,1)**2 + 
     &   4*UStop2(1,2)*(-(Ce(142)*UStop2(1,1)) + 
     &      Ce(143)*UStop2(1,2))

        Ce(1464) = -2*MGl2*Ce(1461)*MStop2(3)*UUStop(1,2) + 
     &   MStop2(1)*(Cd(44)*Ce(1463) - MGl2*Cd(46)*UUStop(2,2))

        Ce(1466) = -2*MGl2*Ce(1461)*MStop2(3)*UUStop(1,3) + 
     &   MStop2(1)*(Cd(44)*Ce(1465) - MGl2*Cd(46)*UUStop(2,3))

        Ce(148) = 7*MStop2(2)**2 - Cd(1060)*MStop2(3) + 
     &   3*MStop2(4)**2*UStop2(1,2)

        Ce(1510) = AtC*Ce(1508) - 2*At*MTy*UUStop(1,2) + 
     &   Ce(1509)*UUStop(1,3)

        Ce(1517) = MUE**2*UUStop(1,3)*UUStop(2,1) - 
     &   MUE2*UUStop(1,2)*UUStop(2,3) + 2*MTy*Cd(11036)*UUStop(3,2)

        Ce(1520) = Ce(1518)*UCStop(3,1)*UStop2(1,2) - 
     &   2*At*Cd(9948)*UStop2(1,1)*UUStop(1,3) - 
     &   Ce(1519)*UUStop(1,2)*UUStopC(2,2)

        Ce(1521) = UUStop(1,2)*UUStop(2,1) + 
     &   UUStop(1,1)*UUStop(2,2) - 2*UUStop(3,1)*UUStop(3,2)

        Ce(1526) = At*Cd(9794)*UStop(2,2) + MUE2*UUStop(3,1) - 
     &   MUEC**2*UUStop(3,2)

        Ce(153) = -(Cd(172)*MStop2(2)**2) + 
     &   MStop2(1)*(MStop2(1) + MStop2(3)*UStop2(1,1))

        Ce(1533) = 4*Ce(1518)*UCStop(3,1)*UStop2(1,2) - 
     &   8*At*Cd(9948)*UStop2(1,1)*UUStop(1,3) - 
     &   4*Ce(1519)*UUStop(1,2)*UUStopC(2,2)

        Ce(1540) = Ce(1539)*UStop(2,1) + Ce(1538)*UStop(2,2) + 
     &   MUEC**2*UStop(1,2)*UUStop(2,2)

        Ce(1545) = Ce(35)*Ce(36)*Cr(91)*UStop(1,1) + 
     &   2*AtC*MTy*Cd(11019)*UStop(1,2) + 
     &   8*MTy2*UStop(2,2)*UUStop(1,2)

        Ce(1547) = At*Cd(9809)*UStop(1,2) + MUE2*UUStop(3,1) - 
     &   MUEC**2*UUStop(3,2)

        Ce(1550) = MUE2*Cd(11106)*UStop(2,2) - 
     &   MUE**2*UStop(1,1)*UUStop(2,1) - 2*MTy*Ce(1549)*UUStop(3,2)

        Ce(1552) = Cd(9809)*UCStop(3,1)*UStop2(1,2) - 
     &   2*MTy*Cd(11161)*UStop2(1,1)*UStopC(2,2) - 
     &   Ce(1551)*UUStopC(2,3)

        Ce(1559) = 2*MUE2*Cr(91)*UUStop(1,3) - 
     &   MUE**2*UUStop(1,1)*UUStop(3,1) - Ce(1558)*UUStop(3,2)

        Ce(1561) = At*Ce(1560)*UStop(1,2) + 
     &   MUE2*Cd(11106)*UUStop(1,1) - MUEC**2*Cd(10925)*UUStop(1,2)

        Ce(1570) = -(MUEC**2*Ce(1568)) + 2*At*MTy*Ce(1569) + 
     &   8*MTy2*UUStop(1,1)

        Ce(1581) = -2*Ce(1580)*MStop2(1)*MStop2(2) + 
     &   (2*MStop2(2)**2 - MStop2(3)*MStop2(4))*UUStop(1,1)*
     &    UUStop(2,2)

        Ce(158) = Cd(9856)*MStop2(1)*MStop2(2) + 
     &   2*Cr(85)*MStop2(2)**2 + 4*MStop2(1)**2*UStop2(1,2)

        Ce(1586) = 4*Ce(1584)*UUStop(1,3) - 
     &   4*Ce(1585)*UUStop(2,3) + 2*MTy*Ce(1583)*UUStop(3,2)

        Ce(1587) = 8*MTy*Cd(3852)*UUStop(3,1) + 
     &   At*(4*Cd(5519)*UUStop(3,3) + 
     &      8*MStop2(1)*MStop2(2)*UUStop(3,4))

        Ce(1590) = -2*Ce(1589)*MStop2(1)*MStop2(2) + 
     &   (2*MStop2(2)**2 - MStop2(3)*MStop2(4))*UUStop(1,2)*
     &    UUStop(2,1)

        Ce(1591) = 8*MTy*Ce(1590) + 
     &   At*(8*MStop2(1)*MStop2(2)*UUStop(1,3)*UUStop(2,2) + 
     &      4*Cd(5519)*UUStop(1,2)*UUStop(2,3))

        Ce(1592) = -8*Cd(9948)*Ce(1588)*UCStop(3,3) + 
     &   Ce(1587)*UCStop(3,1)*UStop2(1,2) + Ce(1591)*UUStopC(2,2)

        Ce(1598) = At*Ce(1597) + 
     &   MTy*(16*MStop2(1)*MStop2(2)*UUStop(1,3)*UUStop(2,1) - 
     &      8*Cd(523)*UUStop(1,1)*UUStop(2,3))

        Ce(1601) = 8*MTy2*Ce(1600) + 
     &   At*(-2*Ce(1599)*MStop2(1)*MStop2(2)*UUStop(1,2) + 
     &      4*MTy*Cd(523)*UStop(1,1)*UUStop(3,2))

        Ce(1605) = 8*MTy*Ce(1581) + 
     &   At*(4*Cd(5519)*UUStop(1,3)*UUStop(2,2) + 
     &      8*MStop2(1)*MStop2(2)*UUStop(1,2)*UUStop(2,3))

        Ce(1607) = -4*MUE2*MStop2(1)*MStop2(2) + 
     &   2*Ce(1509)*MStop2(2)**2 - Ce(1509)*MStop2(3)*MStop2(4)

        Ce(1616) = 8*MTy*Cd(3852)*UUStop(3,1) + 
     &   At*(8*MStop2(1)*MStop2(2)*UUStop(3,3) + 
     &      4*Cd(5519)*UUStop(3,4))

        Ce(1618) = At*Ce(1617) + 
     &   8*MTy*(Cd(523)*UUStop(1,3)*UUStop(2,1) - 
     &      2*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,3))

        Ce(1619) = Cd(9948)*Ce(1616)*UCStop(3,3) + 
     &   4*Cd(11181)*MStop2(2)*MStop2(3)*UCStop(1,3)*UStopC(2,2) + 
     &   Ce(1618)*UUStopC(2,2)

        Ce(162) = UStop2(1,1)*
     &    (Ce(161) - 2*Cd(535)*MStop2(4)*UStop2(1,1)) - 
     &   2*Ce(143)*UStop2(1,2)

	Ce(1643) = -Ce(1686) + Ce(1642)*UStop2(1,2)

        Ce(1644) = 4*MTy2*UUStop(1,1) + 
     &   Ce(186)*Ce(187)*UUStop(1,2) + 4*At*MTy*UUStop(1,3)

        Ce(1646) = 4*Ce(1645)*UUStop(1,1) - 
     &   Ce(186)*Ce(187)*UUStop(1,2) + 4*At*MTy*UUStop(1,3)

        Ce(1651) = Ce(1649)*MStop2(3)**2*UCStop(1,3) - 
     &   4*Cd(9809)*MStop2(1)*MStop2(3)*UCStop(3,3) + 
     &   2*Ce(1650)*UStopC(1,2)

        Ce(1661) = AtC*Ce(1660) + Ce(1659)*UStopC(2,1) - 
     &   Ce(272)*UStopC(2,2)

        Ce(1662) = Ce(35)*Ce(36)*UUStop(1,1) + 
     &   4*MTy2*UUStop(1,2) + 4*AtC*MTy*UUStop(1,3)

        Ce(1665) = Cd(11068)*UStop(2,2) + 
     &   Cd(11276)*UUStop(1,2) - MUE**2*UUStop(1,1)*UUStop(2,1)

        Ce(1693) = 2*Ce(1692)*MStop2(1)*MStop2(3) + 
     &   4*Cd(10980)*MStop2(3)**2*UCStop(1,3)*UStop2(1,2)**2 - 
     &   2*Ce(1650)*UStopC(1,2)

	Ce(1698) = 1 + 4*UStop2(1,1)*UStop2(1,2) - UStop2(1,2)**2

        Ce(1702) = UStop2(1,1) - 
     &   UStop2(1,2)**2*(2*Ce(1701) + UStop2(1,2))

        Ce(1716) = -(MUEC**2*Ce(1683)) + Cd(9794)*Ce(1685) + 
     &   Ce(1715)*MStop2(2)

        Ce(1720) = 4*Cd(11146)*MStop2(1)*MStop2(3) - 
     &   4*Cd(10980)*MStop2(3)**2*UCStop(1,3)*UStop2(1,2)**2 - 
     &   2*Cd(11145)*UStopC(1,2)

        Ce(1722) = 3*MStop2(3)**2 + 
     &   2*MStop2(1)*(Ce(492)*MStop2(3) + 2*MStop2(1)*UStop2(1,1))*
     &    UStop2(1,2)

	Ce(1725) = 2 - 2*UStop2(1,1)*UStop2(1,2) + UStop2(1,2)**2

        Ce(1727) = Ce(1726)*MStop2(1)**2 + 
     &   MStop2(3)*(2*Ce(1725)*MStop2(1) - Cd(970)*MStop2(3))

        Ce(1729) = UStop2(1,1)*(UStop2(1,1) - 8*UStop2(1,2)) + 
     &   UStop2(1,2)**2

        Ce(1730) = 3*UStop2(1,1) + 
     &   (4*Cd(11622) - UStop2(1,2))*UStop2(1,2)**2

        Ce(1731) = Ce(1730)*MStop2(1)*MStop2(3) + 
     &   (Ce(1729)*MStop2(1)**2 - Cd(970)*MStop2(3)**2)*UStop2(1,2)

        Ce(1732) = At**2*Ce(1721)*MStop2(1)*UCStopC(1,3)**2 + 
     &   UStop2(1,1)*(MUE2*Ce(1731) + 
     &      AtC**2*Cd(11112)*UCStop(1,3)*UStopC(1,2))

        Ce(174) = Cd(427) - At**2*UCStopC(1,3)**2 - 
     &   2*Ce(1)*UStop2(1,1)*UStop2(1,2)

        Ce(175) = 8*MTy*Cd(585)*Cr(75)*MStop2(3) - 
     &   4*Ce(174)*Cr(75)*MStop2(4) + AtC**2*Ce(173)*UCStop(1,3)**2

        Ce(1757) = At*Cd(11333)*UStop(2,2) + 
     &   MUE2*Cd(10925)*UUStop(2,1) - MUEC**2*Cd(11106)*UUStop(2,2)

        Ce(1767) = 2*MTy*Cd(11502)*UStop(2,1) + 
     &   Cd(9809)*UStop2(1,2)**2 + 2*Cd(11330)*Cr(91)*UUStopC(3,4)

        Ce(1768) = Ce(1767)*UCStop(3,1) - 
     &   2*MTy*Cd(10918)*UStop(1,1)*UStop2(1,1) - 
     &   Cd(11518)*UStopC(2,1)

	Ce(177) = MStop2(1) - Cd(595)*UStop2(1,1) + Cd(535)*UStop2(1,2)

        Ce(178) = 4*MTy*Cd(585)*Cr(75)*MStop2(3) - 
     &   2*Ce(174)*Cr(75)*MStop2(4) - AtC**2*Ce(177)*UCStop(1,3)**2

        Ce(1788) = -3*Ce(1787) + 
     &   UStop2(1,1)*(5*UStop2(1,1) - 8*UStop2(1,2))

	Ce(179) = -(SA2*Cd(768)) + CA2*Cd(9837) - 2*S2A*Ce(178)

        Ce(1802) = 2*Cd(11205) + 
     &   2*MUEC**2*UUStop(1,3)*UUStop(2,2) + Ce(1801)*UUStop(3,1)

        Ce(1805) = 2*At*MTy*Cd(10995) + 
     &   2*Ce(186)*Ce(187)*UUStop(1,2)*UUStop(2,3) - 
     &   Ce(1804)*UUStop(3,1)

        Ce(1806) = AtC*Cd(11335)*UCStop(1,3)*UStop(2,1) - 
     &   2*(Ce(1805)*MStop2(2)*UStopC(1,2) + 
     &      Ce(1757)*MStop2(3)*UStopC(2,2))

        Ce(181) = Cd(12933) - MStop2(2)*UStop2(1,1) - 
     &   MStop2(1)*UStop2(1,2)

        Ce(1821) = 4*MTy*Cr(91)*UCStop(3,1)*UStop2(1,1) - 
     &   Cd(11037)*UStop2(1,2)*UStopC(2,1) + Cd(11039)*UUStopC(2,3)

        Ce(1822) = 6*MTy*Cd(10918)*UStop(2,1) + 
     &   Cd(9809)*UStop2(1,2)**2 + 2*Cd(11330)*Cr(91)*UUStopC(3,4)

        Ce(1823) = Ce(1822)*UCStop(3,1) - 
     &   2*MTy*Cd(11316)*UStop(1,1)*UStop2(1,1) - 
     &   Cd(11366)*UStopC(2,1)

	Ce(1829) = 3 + 2*UStop2(1,1)*UStop2(1,2) + 3*UStop2(1,2)**2

        Ce(1830) = Ce(1671)*MStop2(1)**2 + 
     &   MStop2(2)*(-(Ce(1829)*MStop2(2)) + 
     &      2*Cr(75)*MStop2(1)*UStop2(1,2))

        Ce(183) = Ce(182) - 4*Cd(535)*UStop2(1,1) + 
     &   4*MStop2(1)*UStop2(1,2)

	Ce(1836) = 3 + UStop2(1,1)*(-UStop2(1,1) + UStop2(1,2))

        Ce(1837) = Ce(1671)*MStop2(1)**2 + 
     &   MStop2(2)*(-2*Ce(1836)*MStop2(1) + Cd(11316)*MStop2(2))

        Ce(1839) = -(Cd(2384)*UStop2(1,1)) + UStop2(1,1)**3 - 
     &   3*UStop2(1,2)

	Ce(1842) = Cd(2366) - UStop2(1,1)**2*UStop2(1,2)

        Ce(1846) = Cd(13134)*UStop(1,1) + 
     &   UCStop(3,1)*UStop(2,1)*UStop2(1,1) + 
     &   Cd(10925)*UStop2(1,2)*UStopC(2,2)

        Ce(1851) = AtC*Ce(1850)*UCStop(1,3)*UStop(2,1) - 
     &   2*Ce(1805)*MStop2(2)*UStopC(1,2) - 
     &   3*Ce(1757)*MStop2(3)*UStopC(2,2)

	Ce(1854) = 2*Ce(1761) + Ce(1811)

	Ce(1873) = 1 - 6*UStop2(1,1)*UStop2(1,2) - 4*UStop2(1,2)**2

	Ce(1874) = -2*Cd(2366) + UStop2(1,1)

	Ce(1883) = 11 + 2*UStop2(1,1)*(-2*UStop2(1,1) + UStop2(1,2))

        Ce(192) = 8*MStop2(1)*UStop2(1,1)**3 - 
     &   2*Cd(594)*UStop2(1,1)**4 + Ce(191)*UStop2(1,2)

        Ce(196) = Cd(520)*UStop2(1,1)**2 + 4*UStop2(1,2)**2 + 
     &   4*UStop2(1,1)*UStop2(1,2)**2

	Ce(2) = Ce(12) - Ce(1)*Cr(28)

        Ce(20) = -2*(MTy*Cr(1)*UStop(1,1) + 
     &      At*MStop2(2)*UStop(1,2)) + Ce(19)*UStop2(1,2)

	Ce(198) = Ce(197)*UStop2(1,1) + UStop2(1,2)*(1 + 3*UStop2(1,2))

	Ce(200) = (Cd(519) - UStop2(1,1))*UStop2(1,1) + 4*UStop2(1,2)**2

        Ce(214) = -(Ce(213)*MStop2(2)) + 
     &   UStop2(1,1)*(Ce(212)*UStop2(1,1) + Ce(211)*UStop2(1,2))

        Ce(22) = -2*MTy*Cd(535)*UStop(1,1) - 
     &   At*MStop2(1)*UStop(1,2) + Ce(21)*UStop2(1,2)

	Ce(235) = 8*Cd(425) + 2*MUE2*Ce(233) - 8*Ce(234)

        Ce(241) = -4*MUE2*Cd(12749) + 4*At*AtC*Ce(239) - 
     &   2*Ce(12)*MStop2(4)

        Ce(243) = (-Ce(241) - 8*MTy2*Cd(492)*MStop2(3) + 
     &     8*MTy*Ce(242)*MStop2(3))/MStop2(3)

	Ce(244) = 2*SA2*Cd(588) + 2*CA2*Ce(238) + S2A*Ce(243)

        Ce(266) = Cd(9794)*Ce(248) - AtC*Ce(265)*UStop(1,1) + 
     &   2*MUEC*Ce(264)*UStop(1,2)

	Ce(277) = Ce(262) - 2*C2A*MUEC*Cd(9794)*UStop(1,2)

        Ce(278) = S2A*MStop2(3)*UCStopC(3,1) + 
     &   Ce(277)*UUStopC(3,1) - Ce(266)*UUStopC(3,3)

        Ce(295) = Cd(9809)*Ce(284) - AtC*Ce(294)*UStop(2,1) + 
     &   2*MUEC*Ce(293)*UStop(2,2)

	Ce(302) = Ce(291) - 2*C2A*MUEC*Cd(9809)*UStop(2,2)

        Ce(303) = S2A*MStop2(3)*UCStop(3,1) - 
     &   Ce(302)*UUStopC(3,1) + Ce(295)*UUStopC(3,4)

        Ce(321) = Cd(9809)*Ce(248) - AtC*Ce(294)*UStop(1,1) + 
     &   2*MUEC*Ce(264)*UStop(2,2)

        Ce(322) = Cd(9794)*Ce(284) + 2*MUEC*Ce(293)*UStop(1,2) - 
     &   AtC*Ce(265)*UStop(2,1)

	Ce(335) = S2A*Ce(285) + AtC*(2*Ce(334) - AtC*S2A*UUStop(2,1))

	Ce(339) = Cd(9809)*Ce(284) - AtC*Ce(337) + 2*MUEC*Ce(338)

        Ce(350) = S2A*MStop2(3)*UStop(1,1) + 
     &   2*Ce(271)*UStopC(2,1) + Ce(321)*UStopC(2,2)

        Ce(351) = S2A*MStop2(3)*UStop(1,2) + 
     &   Ce(322)*UStopC(2,1) + 2*Ce(275)*UStopC(2,2)

	Ce(357) = S2A*Ce(66) + AtC*(2*Ce(356) - AtC*S2A*UUStop(1,1))

	Ce(361) = Cd(9794)*Ce(248) - AtC*Ce(359) + 2*MUEC*Ce(360)

        Ce(370) = Ce(369)*MStop2(1) + 
     &   Cd(50)*Cd(9950)*(-(Cd(9968)*T134(MGl2,MTy2,MStop2(1))) + 
     &      Cd(9976)*T134(MGl2,MTy2,MStop2(2)))

        Ce(372) = Cd(50)*Cd(9947)*Cd(9966)*MGlpTmSt2(1) - 
     &   Cd(47)*Cd(9949)*MGlpTmSt2(2)*
     &    (Cd(9948)*Cd(9983) + Cr(75)**2*UStop2(1,1))

        Ce(379) = Ce(378)*MStop2(1) + 
     &   Cd(50)*Cd(9950)*(-(Cd(10002)*T134(MGl2,MTy2,MStop2(1))) + 
     &      Cd(10006)*T134(MGl2,MTy2,MStop2(2)))

        Ce(386) = (4*Cd(20)*Cd(11651) + 
     &     2*Cd(17)*Cd(11704)*MStop2(1)**2*MStop2(2)**2*
     &      MStop2(3)**2 - 
     &     8*Cd(10905)*MStop2(1)*MStop2(2)*MStop2(3)**2*UStop2(1,1)
     &     )/(MStop2(1)*MStop2(2)*MStop2(3)**2)

        Ce(39) = 2*MTy*MStop2(4)*UCStop(1,3) + 
     &   At*(2*MStop2(1) + MStop2(4)*UStop2(1,2))

        Ce(400) = (-2*Cd(10728)*Ce(399) + 
     &     Cd(50)**2*Cd(10872)*Ce(398)*MStop2(2)*MStop2(3) + 
     &     4*MGl2*Cd(50)**2*Cd(10268)*MStop2(2)**2*MStop2(3))/
     &   (Cd(50)**2*MStop2(2)*MStop2(3))

        Ce(41) = -2*MTy*Cr(26)*UStop(1,1) - 
     &   At*MStop2(1)*UStop(1,2) + Ce(40)*UStop2(1,2)

	Ce(411) = A0delT*MGl2 - MTy2*Ce(410) + Ce(408)*Cr(66)

        Ce(416) = 2*MStop2(2)*
     &    (Ce(415) - 2*Cd(11420)*MStop2(2)) - 
     &   2*Cd(12)*Cd(13)*Cd(11527)*MStop2(3)**2

        Ce(417) = -2*Cd(11419)*MStop2(2) + 
     &   4*Cd(11420)*MStop2(2)**2 + 
     &   4*Cd(12)*Cd(13)*Cd(11375)*MStop2(3)**2

        Ce(423) = -((Cd(9951)*Cd(9982))/Cd(47)**2) + 
     &   (At*Cd(50)*Cd(9950)*Cd(9955)*Cd(9979) + 
     &      Ce(422)*MStop2(3))/(Cd(50)*Cd(9950)*MStop2(3)**2)

	Ce(425) = -M_3C + AtC*(Cd(9948) - UStop2(1,2))

        Ce(43) = 2*MTy*MStop2(4)*UStop(1,1) + 
     &   At*MStop2(3)*UStop(1,2) + Ce(42)*UStop2(1,2)

        Ce(431) = (Cd(9957)*Cd(9980))/(Cd(47)**2*Cd(50)**2) + 
     &   (At*Cd(47)*Cd(50)*Cd(9955)*Cd(9956) + Ce(430)*MStop2(3))/
     &    (Cd(47)*Cd(50)*MStop2(3)**2)

        Ce(433) = (Cd(9980)*Cd(9993))/(Cd(47)**2*Cd(50)**2) + 
     &   (AtC*Cd(47)*Cd(50)*Cd(9955)*Cd(9956) + Ce(432)*MStop2(3))/
     &    (Cd(47)*Cd(50)*MStop2(3)**2)

        Ce(436) = -((Cd(9951)*Cd(9988))/Cd(47)**2) + 
     &   (AtC*Cd(50)*Cd(9950)*Cd(9955)*Cd(9979) + 
     &      Ce(435)*MStop2(3))/(Cd(50)*Cd(9950)*MStop2(3)**2)

        Ce(441) = (Cd(9951) - Cd(9953))*MStop2(3) + 
     &   2*Cd(47)*Cd(9947)*Cd(9955)*MGlpTmSt2(1)*UStop2(1,2)

        Ce(442) = -(Cd(50)**2*Ce(441)) + 
     &   Cd(47)**2*(Cd(9952) - Cd(9954))*MStop2(3)

        Ce(444) = Cd(10031)*Cd(10033) + 
     &   Cd(47)*Cd(50)*(Cd(50)*Cd(9947)*Cd(10026) - 
     &      Cd(47)*Cd(9949)*Cd(10151))*UCStop(1,3)

        Ce(446) = At*Cd(10031)*Cd(10039) + 
     &   Cd(47)*Cd(50)*(Cd(50)*Cd(9947)*Cd(10026) - 
     &      Cd(47)*Cd(9949)*Cd(10151))

        Ce(451) = Cd(50)**2*Cd(9953) + Cd(47)**2*Cd(9954) - 
     &   Cd(9980)*UStop2(1,1)

        Ce(459) = Cd(9948)*Ce(458) - Ce(457)*UCStop(1,3)**2 + 
     &   2*M_3*Cr(75)**2*UStop2(1,1)

        Ce(470) = 2*Cd(9948)*Cd(10181) + 
     &   Ce(469)*UCStopC(1,3)**2 + 2*M_3C*Cr(75)**2*UStop2(1,1)

        Ce(486) = UStop2(1,1)*(2*UStop2(1,1) - 4*UStop2(1,2)) + 
     &   UStop2(1,2)**2

        Ce(500) = AtC*Cd(795)*Cd(12323) + 
     &   At*Cd(10219)*Cd(12323) + 4*Cd(10052)*Cd(12298)*MStop2(1)

        Ce(503) = AtC*Cd(164)*Cd(2597)*Cd(10221) + 
     &   At*Cd(164)*Cd(2597)*Cd(10225) - 2*Ce(502)

        Ce(504) = (MTy*(Ce(500) - Ce(503)) + 
     &     2*Cd(13)*Cd(12317)*MStop2(1))/
     &   (MTy*Cd(47)**2*MStop2(1)**2)

	Ce(509) = -(MTy*Cd(13)) + 2*M_3C*Cd(795) + 2*M_3*Cd(10219)

	Ce(518) = MTy2*Cd(10058) + Ce(511)

	Ce(525) = MGl2 - MStop2(1) - 2*M_3*MTy*UCStopC(1,3)

        Ce(53) = At**2*MStop2(1)*UCStopC(1,3)**2*UStop2(1,1) - 
     &   Ce(49)*UCStopC(1,3)*UStopC(1,1) + 
     &   Ce(52)*UStop2(1,2)*UUStopC(1,2)

        Ce(528) = MGl2*Cd(164)*Cd(10059) + 
     &   MTy2*Cd(10059)*MStop2(1) + 
     &   4*MTy*Cd(164)*Cd(10840)*UStop2(1,1)

	Ce(533) = MTy2*Cd(10055) + MGl2*Ce(530) - 4*MTy*Ce(532)

	Ce(534) = MGlpT2*Ce(530) - 4*MTy*Ce(532)

        Ce(539) = Cd(9963)*UCStop(1,3)**2 + 
     &   M_3*UStop2(1,1)*UStop2(1,2) + At*UStop2(1,1)*UStop2(1,2)

        Ce(540) = 4*Cd(118)*Cd(10169) + At*MGlmT2*Cr(28) - 
     &   2*AtC*MGlmT2*UCStop(1,3)**2

        Ce(544) = Cd(9995)*UCStopC(1,3)**2 + 
     &   M_3C*UStop2(1,1)*UStop2(1,2) + 
     &   AtC*UStop2(1,1)*UStop2(1,2)

        Ce(545) = 4*Cd(118)*Cd(10176) + AtC*MGlmT2*Cr(28) - 
     &   2*At*MGlmT2*UCStopC(1,3)**2

        Ce(550) = Cd(9948)*Ce(549) + 
     &   2*Cd(47)*Cd(164)*Ce(517)*MStop2(1) - 
     &   Cd(10382)*Cd(12371)*MStop2(3)

        Ce(552) = Cd(9948)*Ce(551) + 
     &   2*Cd(47)*Cd(164)*Ce(530)*MStop2(1) - 
     &   Cd(10388)*Cd(12371)*MStop2(3)

        Ce(559) = At*Cd(12471) - 4*Cd(47)*Ce(558)*MStop2(1) + 
     &   Cd(12371)*Ce(554)*MStop2(3)

        Ce(564) = AtC*Cd(12471) - 4*Cd(47)*Ce(563)*MStop2(1) + 
     &   Cd(12371)*Ce(560)*MStop2(3)

	Ce(567) = Ce(609) + MStop2(3)

        Ce(576) = MTy2*(4*M_3*MUE*Ce(575) + 
     &      4*M_3C*MUEC*Ce(575)) + 
     &   Cd(44)*Cd(152)*Cd(424)*Cd(2597)*UStop2(1,1)

        Ce(582) = Ce(581)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(9948)*Ce(580) + 
     &      At*Cd(10284)*UStop2(1,2)**2)

        Ce(594) = Ce(593)*UCStopC(1,3)**2 + 
     &   UStop2(1,1)*(Cd(9948)*Ce(592) + 
     &      AtC*Cd(10284)*UStop2(1,2)**2)

	Ce(605) = 4*M_3C*Ce(602) + 4*M_3*Ce(603) - Cd(152)*Ce(604)

        Ce(607) = At*Cd(10091)*Cr(5) + 
     &   (AtC*Cd(10088) + Cd(9955)*Cd(10500))*UCStop(1,3)

	Ce(61) = 2*MTy*UStop(1,1) + At*UStop(1,2) + Ce(60)*UStop2(1,2)

	Ce(617) = Cd(47) + MStop2(1)**2 - MGlmT2*Cd(1281)*UStop2(1,1)

        Ce(62) = Ce(59)*UStop2(1,1)**2 - 
     &   4*MTy*Cr(85)*UCStop(1,3)*UStop2(1,2) - 
     &   Ce(61)*UStop2(1,1)*UStopC(1,2)

        Ce(618) = Ce(616)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(-2*MStop2(1)*MStop2(3)*UStop2(1,1) + 
     &      Ce(617)*UStop2(1,2))

        Ce(620) = -(At*Ce(619)) - 
     &   2*AtC*MStop2(1)*UCStop(1,3)**2 + 
     &   2*M_3*MStop2(3)*UStop2(1,1)

        Ce(621) = Cd(47) - Cd(152)*MStop2(2) - 
     &   MGlmT2*MStop2(1)*UStop2(1,2)

        Ce(625) = -(AtC*Ce(619)) - 
     &   2*At*MStop2(1)*UCStopC(1,3)**2 + 
     &   2*M_3C*MStop2(3)*UStop2(1,1)

	Ce(636) = AtC*Ce(633) - 2*Cd(152)*Ce(634) + At*Ce(635)

        Ce(643) = M_3*Cd(13)*Cr(28) + 
     &   M_3C*Ce(642)*UCStop(1,3) - 
     &   Cd(9948)*Cd(10834)*UStop2(1,1)

	Ce(644) = Ce(532) - Cd(9948)*Cd(10840)*UStop2(1,1)

	Ce(645) = Cd(13)*Cd(10508) - 2*(AtC*Ce(643) + At*Ce(644))

        Ce(655) = AtC*Cd(795)*Cd(12357) + 
     &   At*Cd(10219)*Cd(12357) + 4*Cd(4187)*Cd(10052)*MStop2(2)

	end
