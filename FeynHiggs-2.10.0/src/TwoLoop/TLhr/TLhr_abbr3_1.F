* TLhr_abbr3_1.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr3_1
	implicit none

#include "TLhr.h"


        Opt(34) = (-MTy2 + MStop2(1))*UStop2(1,1) + 
     &   (MTy2 - 2*MStop2(2) + MStop2(3))*UStop2(1,2)

        Opt(35) = MUE*(AtC*MGlpTmSt4(1) + 
     &     MUE*TB*(MGl2**2 + 2*MGl2*Cd(44) + (MTy2 - MStop2(1))**2)
     &     )*UCStop(1,3)**2*UStop2(1,2)

        Opt(36) = MTy2*MUE**2*UCStop(1,3)**2 - 
     &   MUEC**2*UCStopC(1,3)*
     &    (MGlpT2*UCStopC(1,3) + 
     &      2*M_3C*MTy*UStop2(1,1)*UStop2(1,2))

        Opt(37) = Cd(47)*Cd(9950)*
     &    (-2*Cd(50)**2*Cd(9947)*MGlpTmSt2(1) + 
     &      Cd(47)*Cd(9952)*MStop2(2))*MStop2(3) + 
     &   2*Cd(47)**2*Cd(50)**2*Cd(9950)*
     &    (-T134(MGl2,MTy2,MStop2(1)) + T134(MGl2,MTy2,MStop2(2)))+
     &     Cd(50)*MStop2(1)*
     &    (Cd(50)*Cd(9950)*Cd(9951)*MStop2(3) + 
     &      2*Cd(47)**2*(Cd(9949)*Cd(9950)*MGlpTmSt2(2) + 
     &         Cd(50)*Cd(9947)*MGlpTmSt2(1)*
     &          (1 - Cd(9948) + UStop2(1,2))))

        Opt(39) = -(AtC*MGlpT2*MUE*UCStop(1,3)**2) + 
     &   MGlpT2*(-(MUE**2*TB*UCStop(1,3)**2) + 
     &      MUEC*Cd(12701)*UCStopC(1,3)**2)

        Opt(40) = Opt(39) + 
     &   UStop2(1,1)*(AtC*(-(MTy2*MUEC*UStop2(1,2)) + 
     &         MGl2*MUEC*(-(Cd(9948)*UStop2(1,1)) - UStop2(1,2)**2)
     &         ) + At*(MGl2*MUE*UStop2(1,2) + 
     &         MTy2*MUE*(Cd(9948)*UStop2(1,1) + UStop2(1,2)**2)))

        Opt(41) = -(MUE**2*TB*MGlpTmSt4(2)*UCStop(1,3)**2) + 
     &   MGlpT2*MUEC**2*TB*(MGlpT2 - 2*MStop2(2))*UCStopC(1,3)**2

	Co(8) = 2*MTy*(Cd(12700)*UCStop(1,3) - Cd(12701)*UCStopC(1,3))

        Co(994) = MUEC*Cd(11106)*UStop(2,2) - 
     &   2*MUE*UStop(1,1)*UUStop(2,1)

        Co(996) = 4*MUE*Cd(5519)*UUStop(3,1) - 
     &   4*MUEC*Cd(5945)*UUStop(3,2)

        Co(997) = 4*Cd(5519)*UUStop(3,3) + 
     &   8*MStop2(1)*MStop2(2)*UUStop(3,4)

	Co(998) = MUEC*Cd(10925)*UCStop(3,1) + Cd(11011)*UCStop(3,3)

        Co(999) = MUE*MStop2(3)**2*UUStop(2,1) + 
     &   MUEC*MStop2(1)*MStop2(2)*UUStop(2,2)

        Co(1000) = 8*MUE*MStop2(1)*MStop2(2)*UUStop(1,1) - 
     &   MUEC*Cd(5519)*UUStop(1,2)

	Co(1001) = 8*Co(999)*UUStop(1,3) - 4*Co(1000)*UUStop(2,3)

        Co(1002) = Co(998)*Cr(9)*MStop2(1) - 
     &   Cd(11169)*MStop2(2)*UUStopC(3,2)

        Co(99) = AtC**2*Cd(38)*UCStop(1,3)**2 + 
     &   2*At**2*Cr(75)*UCStopC(1,3)**2

        Co(1004) = MUEC*UStop2(1,1)**3 + 
     &   Cd(11012)*UCStop(1,3)*UStopC(1,2)

        Co(1007) = 2*MUE*TB*MStop2(1)*MStop2(2)*UUStop(1,1) - 
     &   At*MStop2(3)**2*UUStop(1,2)

        Co(1008) = 4*Co(1007)*UUStop(2,1) + 
     &   8*At*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,2)

        Co(1009) = MUE*Co(1008) - 
     &   8*MUEC*Cd(12701)*MStop2(1)*MStop2(2)*UUStop(1,2)*
     &    UUStop(2,2)

        Co(1010) = At*MUE*Cd(5519)*UUStop(2,1) - 
     &   8*MUEC*Cd(12701)*MStop2(1)*MStop2(2)*UUStop(2,2)

        Co(1011) = At*MUE*MStop2(1)*MStop2(2)*UUStop(1,1) + 
     &   MUEC*Cd(12701)*MStop2(3)**2*UUStop(1,2)

	Co(1012) = 4*Co(1010)*UUStop(1,3) + 8*Co(1011)*UUStop(2,3)

        Co(1013) = At*MUE*Cd(11106)*UUStop(1,1) - 
     &   2*MUEC*Cd(10925)*Cd(12701)*UUStop(1,2)

        Co(100) = -4*AtC*Cr(75)*UCStop(1,3)**2 + 
     &   At*Cd(750)*UStop2(1,1)*UStop2(1,2)

	Co(101) = MUE*Co(100) + 4*MUEC*Cd(9993)*Cr(75)

	Co(1031) = -MTy2 + Cr(66)

	Co(1032) = Co(1031) + MTy2*Cr(13)

	Co(1033) = Ce(410) - Ce(408)*Cr(13)

	Co(1034) = 2*MGl2*A0delStop(1) - 2*Co(1033)*MStop2(1)

	Co(1035) = -A0delT + MTy2*Cr(15)

	Co(1036) = Co(1035) - Co(1031)*Cr(13)

	Co(1037) = 2*MTy2*A0delStop(1) - 2*Co(1036)*MStop2(1)

	Co(1039) = MGl2 + Cr(65)

	Co(1041) = Ce(408) + MGl2*Cr(13)

	Co(1042) = -(MTy2*Cd(164)*Co(1034)) + Cd(152)*Cd(2597)*Co(1041)

	Co(1043) = Cd(164)*Cd(2597)*Co(1032) - MGl2*Cd(152)*Co(1037)

	Co(1044) = Ce(408) + MGl2*Cr(34)

	Co(1045) = Co(1031) + MTy2*Cr(34)

	Co(1046) = Ce(410) - Ce(408)*Cr(34)

	Co(1047) = 2*MGl2*A0delStop(2) - 2*Co(1046)*MStop2(2)

	Co(1048) = Co(1035) - Co(1031)*Cr(34)

	Co(1049) = 2*MTy2*A0delStop(2) - 2*Co(1048)*MStop2(2)

	Co(1050) = Cd(179)*Cd(2681)*Co(1044) - MTy2*Cd(251)*Co(1047)

	Co(1051) = Cd(251)*Cd(2681)*Co(1045) - MGl2*Cd(179)*Co(1049)

	Co(1053) = -2*MTy2 + Cr(66)

        Co(1056) = Cd(152)*Cd(164)*Co(1034) + 
     &   8*Co(1041)*MGlpTmSt2(1)*MStop2(1)**2

        Co(1057) = 8*MTy2*Cd(164)*Ce(409) + 
     &   2*Cd(152)*Ce(411)*MGlpTmSt2(1)

        Co(1058) = 2*MGl2*Cd(152)**2*Co(1037) - 
     &   4*Cd(12298)*Co(1032)*MStop2(1)

        Co(1061) = Cd(179)*Cd(251)*Co(1047) + 
     &   8*Co(1044)*MGlpTmSt2(2)*MStop2(2)**2

        Co(1062) = 8*MTy2*Cd(251)*Ce(409) + 
     &   2*Cd(179)*Ce(411)*MGlpTmSt2(2)

        Co(1063) = 2*MGl2*Cd(179)**2*Co(1049) - 
     &   4*Cd(4187)*Co(1045)*MStop2(2)

	Co(1064) = 2*Ce(411) - Co(1037)

	Co(1065) = 2*Ce(411) - Co(1049)

	Co(1066) = Cd(50)*Cd(118)*Co(1064) - Cd(47)*Cd(138)*Co(1065)

	Co(1067) = Co(1047) - Co(1049)

	Co(1068) = Co(1063) - MGlmT2*Cd(50)*Co(1067)

	Co(1069) = -(MGlmT2*Cd(47)*Co(1034)) + Co(1058)

        Co(1070) = MGlmT2*Cd(47)**2*Cd(9948)*Co(1037) - 
     &   Cd(9950)*Co(1069)*UStop2(1,2)

	Co(1071) = Co(1034) + Co(1037)

	Co(1079) = Co(1063) - 2*MGlmT2*Cd(50)*Co(1067)

        Co(1080) = -(AtC*Co(1079)*UCStop(1,3)**2) + 
     &   At*Cd(9948)*Co(1063)*UStop2(1,1)

	Co(1081) = Cd(47) + Cd(9950)

        Co(1082) = MGlmT2*Cd(47)*Cd(9948)*Co(1037)*Co(1081) - 
     &   Cd(9950)*Co(1058)*UStop2(1,2)

	Co(107) = 2*Cd(169)*UStop2(1,1) - 4*UStop2(1,2)

	Co(1086) = Cd(50)**2*Co(1057) + Cd(47)**2*Co(1062)

	Co(1094) = MGlmT2*Cd(47)*Co(1037) + Co(1058)

	Co(1095) = MGlmT2*Cd(47)**2*Co(1034) + Cd(9950)*Co(1094)

	Co(1099) = Co(1050) - Co(1051)

	Co(1100) = Co(1061) - 2*Co(1099)*UStop2(1,1)

	Co(1101) = Cd(9948)*Co(1100) - 2*Co(1099)*UStop2(1,2)**2

        Co(1102) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(50)*Cd(9970)*Co(1067) - 
     &        Co(1101)*MStop2(3)*UStop2(1,1)))/
     &    (SB*SW2*Cd(50)**2*MStop2(3))

	Co(109) = AtC*UCStop(1,3) - Cd(12701)*UCStopC(1,3)

	Co(1104) = Co(1034) - Co(1037)

	Co(1110) = Cd(251)*Cr(28) + Cd(10070)*MGlpTmSt2(2)

	Co(1115) = Cd(251) + Cd(9969)*MGlpTmSt2(2)

	Co(1116) = Cd(164) + Cd(9972)*MGlpTmSt2(1)

	Co(1119) = At*Cd(9970) - 2*AtC*UCStop(1,3)**2

        Co(111) = 8*AtC*Cr(75)*MStop2(4)*UCStop(1,3)**2 + 
     &   At*Cd(9835)*UStop2(1,1)*UStop2(1,2)

	Co(1128) = At*Cd(10016) - 3*AtC*UCStop(1,3)**2

	Co(1130) = 4*Ce(411) - Co(1049)

	Co(1131) = Co(1037) - 2*Cd(9947)*MGlpTmSt2(1)

	Co(112) = MUE*Co(111) - 8*MUEC*Cd(9993)*Cr(75)*MStop2(4)

	Co(1139) = At*Cd(10018) - 3*AtC*UCStop(1,3)**2

	Co(1140) = Co(1042) - Co(1043)

        Co(1141) = -(cI/(16.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(47)*Co(1104)*Co(1139) + 
     &          Cd(9957)*Co(1140)*MStop2(3))*UCStop(1,3))/
     &      (SB*SW2*Cd(47)**2*MStop2(3)))

	Co(1142) = At*Cd(10092) - 3*AtC*UCStop(1,3)**2

        Co(1143) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(-(MGlmT2*Cd(50)*Co(1067)*Co(1142)) + 
     &        Cd(9957)*Co(1099)*MStop2(3))*UCStop(1,3))/
     &    (SB*SW2*Cd(50)**2*MStop2(3))

        Co(1146) = 3*AtC**2*UCStop(1,3)**3 + 
     &   At**2*Cd(10020)*UCStopC(1,3)

	Co(1147) = Co(1047) - 2*Cd(9949)*MGlpTmSt2(2)

	Co(1148) = 4*Ce(411) - Co(1034)

        Co(1155) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*Cd(585)*(-(MGlmT2*Cd(47)*Cd(10091)*Co(1104)) + 
     &        Co(1140)*MStop2(3)*UStop2(1,1)*UStop2(1,2)))/
     &    (SB*SW2*Cd(47)**2*MStop2(3))

        Co(1156) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*Cd(585)*(MGlmT2*Cd(50)*Cd(10020)*Co(1067) - 
     &        Co(1099)*MStop2(3)*UStop2(1,1)*UStop2(1,2)))/
     &    (SB*SW2*Cd(50)**2*MStop2(3))

        Co(1158) = 3*AtC**2*UCStop(1,3)**3 - 
     &   At**2*Cd(10091)*UCStopC(1,3)

        Co(1159) = AtC**2*UCStop(1,3)**3 - 
     &   At**2*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Co(1160) = -(Co(1140)*Co(1159)) + 
     &   MGlmT2*Cd(47)*Co(1104)*UCStop(1,3)

        Co(1161) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Co(1104)*Co(1158) + 
     &        Co(1160)*MStop2(3)))/(SB*SW2*TB*Cd(47)**2*MStop2(3))

        Co(1162) = -(Co(1099)*Co(1159)) + 
     &   MGlmT2*Cd(50)*Co(1067)*UCStop(1,3)

        Co(1163) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(50)*Co(1067)*Co(1146) - 
     &        Co(1162)*MStop2(3)))/(SB*SW2*TB*Cd(50)**2*MStop2(3))

        Co(115) = 4*AtC*Cr(75)*MStop2(4)*UCStop(1,3)**2 + 
     &   At*Cd(822)*UStop2(1,1)*UStop2(1,2)

	Co(1166) = Cd(50)*Co(1131) + Cd(47)*Co(1147)

        Co(1170) = -(At*Co(1079)*UCStopC(1,3)**2) + 
     &   AtC*Co(1063)*UStop2(1,1)*UStop2(1,2)

	Co(1172) = 2*Co(1034) - Co(1037)

	Co(1173) = Co(1058) - MGlmT2*Cd(47)*Co(1172)

	Co(116) = MUE*Co(115) - 4*MUEC*Cd(9993)*Cr(75)*MStop2(4)

	Co(1174) = MGlmT2*Cd(47)**2*Co(1037) - Cd(9950)*Co(1173)

	Co(1177) = 2*Ce(411) + Cd(9947)*MGlpTmSt2(1)

	Co(1186) = Co(1056) - 2*Co(1140)*UStop2(1,1)

	Co(1187) = Cd(9948)*Co(1186) - 2*Co(1140)*UStop2(1,2)**2

        Co(1188) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(-(MGlmT2*Cd(47)*Cd(9969)*Co(1104)*Cr(28)) + 
     &        Co(1187)*MStop2(3)*UStop2(1,1)))/
     &    (SB*SW2*Cd(47)**2*MStop2(3))

	Co(1189) = Cd(164) + Cd(9969)*MGlpTmSt2(1)

	Co(1190) = Cd(251)*Cr(28) - Cd(9970)*MGlpTmSt2(2)

	Co(1197) = AtC*Cd(10070) + 2*At*UCStopC(1,3)**2

	Co(1198) = Cd(251) + Cd(9972)*MGlpTmSt2(2)

	Co(119) = Cr(75)*MStop2(1) + Cr(9)*MStop2(2)

        Co(1206) = AtC**2*Cd(10020)*UCStop(1,3) + 
     &   3*At**2*UCStopC(1,3)**3

        Co(1207) = MGlmT2*Cd(50)*Co(1067) - 
     &   At**2*Co(1099)*UCStopC(1,3)**2

        Co(1208) = Co(1207)*UCStopC(1,3) + 
     &   AtC**2*Co(1099)*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)

        Co(1211) = AtC**2*Cd(10092)*UCStop(1,3) + 
     &   3*At**2*UCStopC(1,3)**3

        Co(120) = AtC**2*Cd(38)*Co(119)*UCStop(1,3)**2 - 
     &   At**2*Cd(794)*UCStopC(1,3)**2

        Co(1214) = -(cI/(16.D0*Pi**3)*
     &     (Alfas2L*Cd(585)*
     &        (MGlmT2*Cd(47)*Cd(10018)*Co(1104) + 
     &          Co(1140)*MStop2(3)*UStop2(1,1)*UStop2(1,2)))/
     &      (SB*SW2*Cd(47)**2*MStop2(3)))

        Co(1217) = AtC**2*Cd(10018)*UCStop(1,3) + 
     &   3*At**2*UCStopC(1,3)**3

        Co(1218) = MGlmT2*Cd(47)*Co(1104) - 
     &   At**2*Co(1140)*UCStopC(1,3)**2

        Co(1219) = Co(1218)*UCStopC(1,3) + 
     &   AtC**2*Co(1140)*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)

        Co(1220) = -(cI/(16.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(47)*Co(1104)*Co(1217) + 
     &          Co(1219)*MStop2(3)))/
     &      (SB*SW2*TB*Cd(47)**2*MStop2(3)))

	Co(1221) = AtC*Cd(10016) - 3*At*UCStopC(1,3)**2

	Co(1222) = AtC*Cd(10091) + 3*At*UCStopC(1,3)**2

        Co(1223) = cI/(16.D0*Pi**3)*
     &   (Alfas2L*(-(MGlmT2*Cd(47)*Co(1104)*Co(1222)) + 
     &        Cd(9993)*Co(1140)*MStop2(3))*UCStopC(1,3))/
     &    (SB*SW2*Cd(47)**2*MStop2(3))

	Co(121) = -(MUE*UCStop(1,3)**2) + MUEC*UStop2(1,1)*UStop2(1,2)

	Co(1224) = AtC*Cd(10020) - 3*At*UCStopC(1,3)**2

	Co(1230) = 2*MTy*UStop(1,1) - MUEC*TB*UStop(1,2)

	Co(1231) = AtC*Co(1230)*UStop(1,1) + 2*MTy*MUE*TB*UUStop(1,1)

	Co(122) = TB*Co(121)*MStop2(3) - MTy*MStop2(4)*UCStop(1,3)

        Co(1234) = MUE*TB*UStop2(1,2)*UUStop(1,1) - 
     &   2*At*UStop2(1,1)*UUStop(1,2)

	Co(1237) = 2*MTy*UStop(2,1) - MUEC*TB*UStop(2,2)

	Co(1240) = 2*Cd(14037)*UStop(1,1) + MUEC*TB*UUStop(1,2)

	Co(1241) = MUE*TB*UUStop(1,1) + Cd(12701)*UUStop(1,2)

	Co(1242) = 2*MTy*Co(1241) - At*MUE*TB*UUStop(1,3)

	Co(1243) = Co(1240) + 2*AtC*UUStop(1,1)

	Co(123) = MUE*Cd(9843)*UCStop(1,3) - MUEC*Cd(9844)*UCStopC(1,3)

	Co(1244) = AtC*Co(1243)*UStop(1,1) + Co(1242)*UStop(1,2)

        Co(1245) = MUE*TB*Cd(11161)*UStop(1,1) + 
     &   MTy*Cd(11137)*UStop(1,2)

	Co(1246) = 2*Co(1245)*UStop(1,1) + MUEC*TB*Cd(11138)*UUStop(1,2)

	Co(1247) = MUE*TB*Cr(91)*UStop(1,1) + At*UStop(2,2)*UUStop(1,2)

        Co(1248) = At*Cr(91)*UStop(1,2) - 
     &   2*MUE*TB*UStop(2,1)*UUStop(1,1)

        Co(1250) = MUE*TB*Cd(10925)*UStop(1,1) + 
     &   3*At*UStop(2,2)*UUStop(1,2)

	Co(1251) = Co(1250) + 3*MUEC*TB*UStop(2,2)*UUStop(1,2)

	Co(1252) = -2*MTy*Co(1251) + At*MUE*TB*Cd(10925)*UStop(1,2)

        Co(124) = 2*MUEC*Cd(794)*UCStopC(1,3)**2 + 
     &   MUE*Cd(9840)*UStop2(1,1)*UStop2(1,2)

	Co(1254) = Co(1246)*UStopC(2,2) - 2*AtC*Cd(11107)*UUStop(1,1)

	Co(1256) = 2*MTy*Cr(91) - MUE*TB*UUStop(3,1)

	Co(1257) = MTy*Cd(11001) - 2*MUE*TB*UUStop(1,1)*UUStop(2,3)

	Co(125) = TB*Co(124)*MStop2(3) - 2*MTy*Cd(9844)*UCStopC(1,3)

	Co(1273) = At*UUStop(1,2) - 2*MTy*UUStop(1,3)

	Co(1275) = -(MTy*Cd(11180)) + MUE*TB*UUStop(3,1)

	Co(1276) = 2*Co(1275)*UStop(1,1) + MUEC*TB*Cd(10921)*UStop(1,2)

	Co(1277) = AtC*Co(1276)*UStopC(2,2) + 8*MTy*MUE*TB*UUStop(1,1)

        Co(1279) = 5*MTy*MUE*TB*UStop(1,1)*UUStop(1,1) + 
     &   4*At**2*UStop(1,2)*UUStop(1,2)

	Co(1280) = MTy*Cd(11106)*UCStop(3,1) + At*Cd(10925)*UCStop(3,4)

	Co(1283) = At*MUE*TB*Cd(3784) + 4*AtC*MTy*MStop2(1)*UCStop(1,3)

        Co(1285) = MUE*TB*Cd(1752)*UUStop(1,1) + 
     &   At*MStop2(3)**2*UUStop(1,2)

        Co(1286) = MUE*TB*Cd(256)*MStop2(2)*UUStop(1,1) + 
     &   At*MStop2(3)**2*UUStop(1,2)

        Co(1288) = At*Co(1285) - 
     &   MUEC*TB*Cd(10980)*MStop2(3)**2*UStop(1,2)

	Co(1290) = MUE*TB*UStop(1,1) + 2*MTy*UStop(1,2)

	Co(1291) = UStop2(1,2)*UUStop(3,3) + UUStop(3,4)

	Co(1292) = MUEC*TB*Co(1291) + 2*MTy*UStop2(1,1)*UUStop(3,1)

	Co(1293) = Co(1290)*UCStop(1,3)*UStop(1,1) - 4*MTy*UUStop(1,1)

	Co(1294) = Co(1293)*UCStop(1,3) - Co(1292)*UCStop(3,3)

	Co(1295) = MTy*Cd(11106)*UStop(1,1) + At*Cd(10925)*UStop(1,2)

        Co(1296) = Co(1295)*MStop2(3)*UStopC(2,1) + 
     &   Cd(11188)*MStop2(1)*UUStopC(2,3)

        Co(1297) = Cd(11187)*UCStop(3,1)*UStop2(1,1) - 
     &   Cd(698)*Cr(91)*UUStopC(2,1)

	Co(1299) = 2*MStop2(3) + MStop2(1)*UStop2(1,1)**2

	Co(1305) = Cd(11106)*UStop2(1,1) + UUStop(3,4)

	Co(130) = Cr(75)*MStop2(1) - Cr(9)*MStop2(2)

	Co(1318) = MUE*TB*UUStop(2,1) + At*UUStop(2,2)

	Co(1319) = MUE*TB*UUStop(3,1) - At*UUStop(3,2)

	Co(1321) = MUE*TB*Cd(10995) + 2*At*UUStop(1,2)*UUStop(2,2)

	Co(1327) = MTy*Cr(91) + MUE*TB*UUStop(3,1)

	Co(1328) = Cd(9948)*Co(1230) + 2*Co(1327)*UStopC(2,2)

	Co(1329) = 2*MTy*Cr(58) + MUE*TB*UUStop(3,1)

	Co(1330) = Co(1329)*UStop(1,1) + MUEC*TB*UStop(2,1)*UUStop(1,2)

	Co(1331) = MUE*TB*Cr(91) + 4*MTy*UUStop(3,2)

	Co(132) = MUE*TB*Cd(9860) + 2*AtC*Ce(156)*Cr(75)

        Co(133) = 2*MTy*(Co(132)*UCStop(1,3) - 
     &     Cd(9865)*Cd(12701)*UCStopC(1,3))

	Co(1346) = 2*Co(1280)*MStop2(1) - 5*Cd(10981)*MStop2(3)

        Co(1347) = -5*Cd(10981)*MStop2(3) + 
     &   2*MUE*TB*MStop2(1)*UUStop(1,1)

	Co(1348) = MUEC*TB*Co(1346) + At*Co(1347)

	Co(1349) = MUE*TB*Cr(1)*UUStop(1,1) + At*MStop2(3)*UUStop(1,2)

	Co(1350) = At*Co(1349) - MUEC*TB*Cd(10980)*MStop2(3)*UStop(1,2)

        Co(134) = At*Cd(865)*UCStopC(1,3)**2 + 
     &   4*AtC*Cd(833)*Cr(75)*UStop2(1,1)*UStop2(1,2)

	Co(1355) = 2*MTy*UUStop(1,1) - MUEC*TB*UUStop(1,3)

        Co(135) = 4*AtC**2*Cd(38)*Co(130)*MStop2(4)*
     &    UCStop(1,3)**2 - 2*At**2*Cd(865)*UCStopC(1,3)**2

	Co(1372) = Cd(11035)*UStop(1,2) + 5*MStop2(3)**2*UUStop(1,3)

	Co(1386) = MUE*TB*Cd(9794) + AtC*Co(1230)

	Co(1387) = (At + MUEC*TB)*Cd(10976)

	Co(1388) = Co(1386)*UUStop(1,1) + 2*Co(1387)*UUStop(1,2)

	Co(1389) = (At + MUEC*TB)*Cd(10981)

	Co(138) = SA - CA*TB

	Co(139) = SA + CA*TB

	Co(1410) = MUE*TB*UStop(2,1) + MTy*UStop(2,2)

	Co(1411) = 2*Co(1410)*UStop(1,1) + MUEC*TB*UUStop(3,2)

	Co(1412) = MUE*TB*UUStop(3,1) + Cd(12701)*UUStop(3,2)

	Co(1413) = 2*MTy*Co(1412) - At*MUE*TB*UUStop(3,4)

	Co(140) = (2*AtC*CA + MUE*Co(139))*MStop2(4)

	Co(1415) = MUE*TB*UStop(2,1) - 2*MTy*UStop(2,2)

	Co(1416) = MUEC*TB*Cr(91)*UStop(2,2) + Co(1415)*UUStop(3,1)

	Co(1417) = MUE*TB*Cr(91)*UStop(2,1) + At*UStop(1,2)*UUStop(2,2)

        Co(1420) = At*Cd(10925)*UStop(2,2) + 
     &   3*MUE*TB*UStop(1,1)*UUStop(2,1)

        Co(1421) = MUE*TB*Cd(10925)*UUStop(2,1) + 
     &   2*At*Cd(11106)*UUStop(2,2)

	Co(141) = 2*CA*MTy*MStop2(3) + Co(140)*UCStop(1,3)

	Co(1424) = Cd(6268)*UUStop(3,3) - UStop2(1,2)**2*UUStop(3,4)

	Co(1428) = MUEC*TB*Cd(11330) + At*MTy*UStop(2,1)

        Co(1430) = UCStop(3,1)*UStop2(1,1)*UUStop(2,1) + 
     &   Cd(10918)*UUStop(3,1)

	Co(142) = MUEC*Cd(435)*Co(138) - 2*Co(141)*UCStop(1,3)

	Co(1434) = MTy*Cr(91) - MUE*TB*UUStop(3,1)

	Co(1435) = Co(1434)*UStop(2,1) - MUEC*TB*Cr(91)*UStop(2,2)

	Co(1436) = UUStop(2,2)*UUStop(3,1) + UUStop(2,1)*UUStop(3,2)

	Co(1439) = -(MUE*TB*UCStop(1,3)) + 2*MTy*UStop2(1,1)

        Co(1442) = MUE*TB*Cd(11106)*UStop(2,1) + 
     &   3*At*UStop(1,2)*UUStop(2,2)

	Co(1443) = 2*MTy*Co(1442) - At*MUE*TB*Cd(11106)*UStop(2,2)

	Co(143) = MUEC*Cd(423)*Co(138) - 2*Co(141)*UCStop(1,3)

	Co(1451) = -2*MTy*Cr(58) + MUE*TB*UUStop(3,1)

        Co(1452) = MUEC*TB*UUStop(1,1)*UUStop(2,2) + 
     &   Co(1451)*UUStop(3,1)

	Co(144) = Co(142) - Co(143)

	Co(1454) = MUEC*TB*Cr(91) - 4*MTy*UUStop(3,1)

	Co(1455) = MUE*TB*UUStop(2,1) + 4*At*UUStop(2,2)

	Co(1459) = Co(1321) + 2*MUEC*TB*UUStop(1,2)*UUStop(2,2)

	Co(1460) = 2*MTy*Co(1459) - At*MUE*TB*Cr(91)*UUStop(3,2)

	Co(1461) = MTy*UUStop(1,2) + 2*MUE*TB*UUStop(1,3)

	Co(1462) = MUEC*TB*Cr(91) + 4*MTy*UUStop(3,1)

	Co(145) = 2*AtC*CA + MUE*Co(139)

	Co(1470) = Cd(2359)*UUStop(3,3) - UStop2(1,2)**2*UUStop(3,4)

	Co(146) = 2*CA*MTy + Co(145)*UCStop(1,3)

	Co(1474) = MUE*TB*UUStop(1,1) + 4*At*UUStop(1,2)

        Co(1478) = UCStop(3,1)*UStop2(1,1)*UUStop(2,1) + 
     &   Cd(11316)*UUStop(3,1)

        Co(147) = Co(146)*UCStop(1,3) + 
     &   MUEC*Co(138)*UStop2(1,1)*UStop2(1,2)

	Co(1484) = -4*MTy*UCStopC(1,3) + MUE*TB*UStop2(1,1)

        Co(1485) = MUE*TB*UCStopC(3,1)*UUStop(1,2) + 
     &   Co(1484)*UUStop(3,2)

	Co(1486) = 5*MStop2(3) + Cd(814)*UStop2(1,1)**2

	Co(1491) = MTy*Cd(11343) - 2*MUE*TB*UCStop(1,3)*UStop2(1,1)

	Co(13) = 4*MTy*(Cd(12700)*UCStop(1,3) - Cd(12701)*UCStopC(1,3))

        Co(148) = 1/2.D0*(A0delStop(1)*Co(142)) + 
     &   Co(147)*Cr(13)*MStop2(1)*MStop2(3)

	Co(1497) = 3*MUE*TB*UUStop(3,1) + At*UUStop(3,2)

	Co(1498) = -(MUEC*TB*Cd(10994)) + At*Co(1497)

	Co(1499) = 2*MTy*Co(1319)*Cr(91) + Co(1498)*UUStop(3,2)

	Co(149) = 2*CA*MTy - Co(145)*UCStop(1,3)

        Co(1510) = 2*Cd(11324)*MStop2(3) + 
     &   4*MStop2(2)*UStop(1,1)*UStop2(1,1)**2

        Co(150) = -(Co(149)*UCStop(1,3)) + 
     &   MUEC*Co(138)*UStop2(1,1)*UStop2(1,2)

	Co(1514) = 3*MStop2(3) + MStop2(1)*UStop2(1,1)**2

	Co(1515) = 3 - Cr(75)*UStop2(1,1)

	Co(1516) = 3*Cd(10918)*MStop2(1) - Co(1515)*MStop2(2)

        Co(1519) = Cd(11412)*MStop2(3) - 
     &   6*MStop2(2)*UStop2(1,1)**2*UUStop(3,4)

	Co(1520) = Co(1491)*UStop(1,1) - Co(1439)*UCStop(3,1)*UStop(2,1)

        Co(1521) = Co(1520)*MStop2(3) + 
     &   Co(1256)*MStop2(2)*UStop2(1,1)*UStopC(2,2)

	Co(1522) = 2*MUEC*TB*Co(1519) + 4*Co(1521)*UStop(2,1)

        Co(151) = 1/2.D0*(A0delStop(2)*Co(143)) - 
     &   Co(150)*Cr(34)*MStop2(2)*MStop2(3)

	Co(1524) = Co(1522) + 2*AtC*Cd(11409)*UCStop(1,3)*UStop(2,1)

        Co(152) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(-(Co(151)*MStop2(1)) + Co(148)*MStop2(2)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2)*MStop2(3))

	Co(1534) = At*MUE*UUStop(1,1) - 2*MTy*MUEC*UUStop(1,3)

	Co(1535) = TB*Co(1534)*UUStop(2,2) + 2*MTy*Co(1319)*UUStop(3,3)

	Co(1536) = Cd(14040) + MUEC*TB*UUStop(1,2)

	Co(1537) = 2*MTy*Co(1536) - At*MUE*TB*UUStop(1,3)

	Co(1539) = At*MUE*UUStop(2,1) - 2*MTy*MUEC*UUStop(2,3)

	Co(1540) = TB*Co(1539)*UUStop(1,2) + 2*MTy*Co(1319)*UUStop(3,4)

	Co(1541) = MTy*UUStop(2,1) + At*UUStop(2,3)

        Co(1542) = MUE*TB*UUStop(1,3)*UUStop(2,1) + 
     &   2*At*UUStop(1,2)*UUStop(2,3)

	Co(153) = Co(142)*Cr(13) + 2*Co(147)*MStop2(3)

	Co(1545) = MTy*UUStop(1,1) + At*UUStop(1,3)

        Co(1546) = 2*At*UUStop(1,3)*UUStop(2,2) + 
     &   MUE*TB*UUStop(1,1)*UUStop(2,3)

	Co(1548) = Co(1318) + MUEC*TB*UUStop(2,2)

	Co(1549) = 2*MTy*Co(1548) - At*MUE*TB*UUStop(2,3)

	Co(1552) = MUEC*TB*UUStop(1,2) - 2*MTy*UUStop(1,3)

	Co(1553) = 2*MTy*UUStop(3,1) - MUEC*TB*UUStop(3,4)

	Co(154) = Co(143)*Cr(34) - 2*Co(150)*MStop2(3)

        Co(1554) = Cd(9948)*Co(1553)*UStop(1,1) + 
     &   Co(1552)*UStopC(2,2)*UUStop(2,1)

	Co(1555) = 2*MTy*UUStop(3,1) - MUEC*TB*UUStop(3,3)

	Co(1556) = MTy*UUStop(1,2) + MUE*TB*UUStop(1,3)

        Co(1557) = 2*Co(1556)*UUStop(2,1) + 
     &   MUEC*TB*UUStop(1,2)*UUStop(2,3)

	Co(1559) = MUEC*TB*UUStop(2,2) - 2*MTy*UUStop(2,3)

	Co(1560) = MTy*UUStop(2,2) + MUE*TB*UUStop(2,3)

        Co(1561) = 2*Co(1560)*UUStop(1,1) + 
     &   MUEC*TB*UUStop(1,3)*UUStop(2,2)

	Co(155) = Co(153) - Co(154)

	Co(1566) = MUE*TB*UUStop(3,1) + 2*At*UUStop(3,2)

        Co(1570) = MStop2(3)**2*UUStop(3,3) - 
     &   4*MStop2(1)*MStop2(2)*UUStop(3,4)

	Co(1571) = Cr(91)*UCStop(3,3) - UStop2(1,2)*UUStop(1,1)

	Co(1572) = Cr(9)*MStop2(3) + 4*MStop2(2)*UStop2(1,1)

	Co(1573) = Cd(5519)*Cd(9948) - Co(1572)*MStop2(1)

	Co(156) = 2*At*CA + MUEC*Co(139)

	Co(1575) = Cr(9)*MStop2(1) - Cd(9948)*MStop2(2)

	Co(1576) = MStop2(2)*UStop2(1,1)**2 + MStop2(1)*UStop2(1,2)**2

        Co(1578) = MUE*TB*MStop2(3)**2*UUStop(1,1) - 
     &   2*At*MStop2(1)*MStop2(2)*UUStop(1,2)

        Co(1579) = 8*MUE*TB*MStop2(1)*MStop2(2)*UUStop(1,2)*
     &    UUStop(2,1) - 4*Co(1578)*UUStop(2,2)

        Co(1580) = 4*Cd(523)*UUStop(3,3) - 
     &   8*MStop2(1)*MStop2(2)*UUStop(3,4)

        Co(1581) = Co(1579) + 
     &   8*MUEC*TB*MStop2(1)*MStop2(2)*UUStop(1,2)*UUStop(2,2)

        Co(1582) = Cd(523)*UUStop(1,3)*UUStop(2,2) - 
     &   2*MStop2(1)*MStop2(2)*UUStop(1,2)*UUStop(2,3)

        Co(157) = MUE*Cd(423)*Co(138) - 
     &   2*Co(156)*MStop2(4)*UCStopC(1,3)**2

        Co(1584) = 4*MUE*TB*Ce(1581) - 
     &   4*At*Cd(3852)*UUStop(1,2)*UUStop(2,2)

        Co(1585) = 2*MTy*MUEC*Cd(3852)*UUStop(1,2) - 
     &   At*MUE*Cd(5519)*UUStop(1,3)

	Co(1589) = MUE*TB*UUStop(2,1) - 3*At*UUStop(2,2)

	Co(1590) = Co(1589) - 3*MUEC*TB*UUStop(2,2)

        Co(1591) = 2*MTy*Co(1590)*UStop(1,2) + 
     &   At*MUE*TB*Cd(11106)*UStop(2,2)

        Co(1593) = MUE*TB*MStop2(1)*MStop2(2)*UUStop(2,1) - 
     &   At*MStop2(3)**2*UUStop(2,2)

	Co(14) = MUEC*Cd(477) + MUE*MStop2(4)*UCStop(1,3)**2

	Co(158) = (-Co(157) + 4*CA*MTy*MStop2(3)*UCStopC(1,3))/MStop2(3)

        Co(1594) = MUE*TB*Cd(5519)*UUStop(1,1) + 
     &   8*At*MStop2(1)*MStop2(2)*UUStop(1,2)

	Co(1595) = 8*Co(1593)*UUStop(1,3) + 4*Co(1594)*UUStop(2,3)

        Co(1596) = 4*At*Ce(1581) - 
     &   4*MUE*TB*Cd(3852)*UUStop(1,1)*UUStop(2,1)

	Co(1597) = 4*(At*Ce(1579) + MTy*Ce(1581))

	Co(1598) = MTy*Co(1596) + MUEC*TB*Co(1597)

	Co(1599) = At*Co(1595) - 2*Co(1598)

        Co(1600) = 4*MUE*TB*Ce(1590) - 
     &   4*At*Cd(3852)*UUStop(1,2)*UUStop(2,2)

        Co(1601) = MTy*MUEC*Cd(3852)*UUStop(1,2) - 
     &   At*MUE*MStop2(1)*MStop2(2)*UUStop(1,3)

	Co(1603) = MUE*TB*Cd(5945)*UUStop(3,1) + At*Cd(5519)*UUStop(3,2)

        Co(159) = -2*Cr(28)*MStop2(1) + 
     &   2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Co(1604) = 16*MTy*Cr(91)*MStop2(1)*MStop2(2) - 
     &   4*At*Cd(5519)*UUStop(3,2)

	Co(1605) = 4*At*Co(1603) - MUEC*TB*Co(1604)

        Co(1606) = 16*MTy*Co(1319)*Cr(91)*MStop2(1)*MStop2(2) + 
     &   Co(1605)*UUStop(3,2)

	Co(1609) = MTy*UStop(1,1) - At*UStop(1,2)

	Co(1610) = 3*MUE*TB*UUStop(2,1) - At*UUStop(2,2)

        Co(1611) = MUE*TB*Cd(10925)*UStop(2,1) + 
     &   2*At*UStop(1,2)*UUStop(2,2)

        Co(160) = MUE*Co(138)*Co(159) + 
     &   2*Co(156)*MStop2(4)*UCStopC(1,3)**2

        Co(1614) = Cd(523)*UUStop(1,3)*UUStop(2,1) - 
     &   2*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,3)

        Co(1616) = 2*MStop2(1)*MStop2(2)*UUStop(3,3) + 
     &   Cd(5519)*UUStop(3,4)

	Co(1617) = -4*MUEC*TB*Co(1616) + 8*MTy*Cd(3852)*UUStop(3,1)

	Co(1618) = -2*Co(1275)*UStop(1,1) - MUEC*TB*Cd(10921)*UStop(1,2)

        Co(1620) = 2*MStop2(1)*MStop2(2)*UUStop(1,2)*
     &    UUStop(2,1) - MStop2(3)**2*UUStop(1,1)*UUStop(2,2)

        Co(1621) = 2*MStop2(1)*MStop2(2)*UUStop(1,3)*
     &    UUStop(2,1) - Cd(523)*UUStop(1,1)*UUStop(2,3)

        Co(161) = (Co(160) + 4*CA*MTy*MStop2(3)*UCStopC(1,3))/
     &   (MZ2*SW2*MStop2(3))

	Co(1623) = MUE*TB*Cd(11657) + 4*MTy*Cd(3852)*Cr(91)*UUStop(3,2)

	Co(1624) = 4*MTy*Ce(1581) - 4*MUE*TB*Ce(1606)

        Co(1625) = Cd(5519)*UUStop(1,3)*UUStop(2,2) + 
     &   2*MStop2(1)*MStop2(2)*UUStop(1,2)*UUStop(2,3)

	Co(1626) = 2*(Co(1624) - 2*MUEC*TB*Co(1625))

	Co(1628) = MUE*TB*Cr(91)*UStop(2,1) + MTy*Cd(11161)*UStop(2,2)

        Co(1629) = MUEC*TB*Cd(11106)*UStop(2,2) + 
     &   2*Cd(13077)*UUStop(2,1)

        Co(1631) = Cd(5519)*UUStop(3,3) + 
     &   2*MStop2(1)*MStop2(2)*UUStop(3,4)

	Co(1632) = -4*MUEC*TB*Co(1631) + 8*MTy*Cd(3852)*UUStop(3,1)

        Co(162) = (3*cI)/(16.D0*Pi)*
     &   (Alfa2L*(-Co(158) + MZ2*SW2*Co(161)))/(SB*SW2*TB)

	Co(1633) = MUE*Cd(5519)*UUStop(3,1) + MUEC*Cd(5945)*UUStop(3,2)

	Co(1634) = -4*TB*Co(1633) + 16*MTy*Cr(91)*MStop2(1)*MStop2(2)

	Co(1635) = MUEC*TB*Cd(10925) - 6*MTy*UUStop(3,1)

	Co(1636) = Co(1635)*UCStop(3,1) - Co(1411)*UCStop(3,3)

        Co(1637) = MStop2(3)**2*UUStop(1,3)*UUStop(2,1) - 
     &   4*MStop2(1)*MStop2(2)*UUStop(1,1)*UUStop(2,3)

	Co(1638) = 4*(MTy*Ce(1590) + MUE*TB*Co(1637))

        Co(1639) = 2*MStop2(1)*MStop2(2)*UUStop(1,3)*
     &    UUStop(2,2) + Cd(5519)*UUStop(1,2)*UUStop(2,3)

	Co(1640) = Co(1638) - 2*MUEC*TB*Co(1639)

        Co(1641) = MUE*TB*Cd(11161)*UUStop(1,1) + 
     &   MTy*Cd(11137)*UUStop(1,3)

	Co(1642) = 2*Co(1641) + MUEC*TB*Cd(11138)*UUStop(1,2)

	Co(163) = AtC*UCStop(1,3) - At*UCStopC(1,3)

        Co(1643) = Co(1636)*Cr(9)*MStop2(1) + 
     &   Co(1642)*MStop2(2)*UUStopC(3,2)

	Co(1648) = At*Cd(11657) + 4*MTy*Cd(3852)*Cr(91)*UUStop(3,1)

        Co(1649) = 4*MUE*TB*MStop2(3)**2*UUStop(1,2)*
     &    UUStop(2,1) - 8*Cd(14040)*MStop2(1)*MStop2(2)*UUStop(2,2)

        Co(1650) = 8*MStop2(1)*MStop2(2)*UUStop(3,3) - 
     &   4*Cd(523)*UUStop(3,4)

        Co(1651) = Co(1649) - 
     &   8*MUEC*TB*MStop2(1)*MStop2(2)*UUStop(1,2)*UUStop(2,2)

        Co(1652) = 2*MStop2(1)*MStop2(2)*UUStop(1,3)*
     &    UUStop(2,2) - Cd(523)*UUStop(1,2)*UUStop(2,3)

	Co(164) = MUE*UCStop(1,3) - MUEC*UCStopC(1,3)

        Co(1654) = MUE*TB*Cd(5519)*UUStop(2,1) + 
     &   8*At*MStop2(1)*MStop2(2)*UUStop(2,2)

        Co(1655) = MUE*TB*MStop2(1)*MStop2(2)*UUStop(1,1) - 
     &   At*MStop2(3)**2*UUStop(1,2)

	Co(1656) = 4*Co(1654)*UUStop(1,3) + 8*Co(1655)*UUStop(2,3)

        Co(1657) = 4*At*Ce(1590) - 
     &   4*MUE*TB*Cd(3852)*UUStop(1,1)*UUStop(2,1)

        Co(1658) = 4*MStop2(1)*MStop2(2)*UUStop(1,3)*
     &    UUStop(2,2) - MStop2(3)**2*UUStop(1,2)*UUStop(2,3)

	Co(1659) = -4*MTy*Ce(1590) + 4*At*Co(1658)

        Co(1661) = MUE*TB*Cr(91)*UStop(1,1) - 
     &   2*At*UStop(2,2)*UUStop(1,2)

        Co(1663) = At*Cd(11106)*UStop(1,2) + 
     &   3*MUE*TB*UStop(2,1)*UUStop(1,1)

        Co(1664) = MUE*TB*Cd(11106)*UUStop(1,1) + 
     &   2*At*Cd(10925)*UUStop(1,2)

        Co(1673) = Cd(9948)*UCStop(3,3)*UUStop(3,3) - 
     &   UUStop(3,1)*UUStop(3,2)*UUStopC(2,2)

        Co(1676) = 2*MStop2(3)*
     &   (2*Co(1554)*MStop2(3) - Co(1618)*MStop2(2)*UCStop(1,3))

	Co(1690) = Cd(537)*Cd(9948) + Cd(10020)*UStop2(1,1)

	Co(15) = MUE*Cd(485) + MUEC*MStop2(4)*UCStopC(1,3)**2

	Co(168) = 6*CA*(AtC*UCStop(1,3) - At*UCStopC(1,3))

        Co(1693) = Cd(9949)*Cd(9958)*MGlpTmSt2(2) + 
     &   4*Ce(411)*UStop2(1,2)

	Co(1694) = Cd(9958)*Co(1047) - 2*Co(1693)

	Co(1702) = Cd(47)*Co(1032) + MGlmT2*Co(1104)

	Co(169) = 6*MUE*UCStop(1,3) - 6*MUEC*UCStopC(1,3)

	Co(1703) = Co(1041) + Co(1032)*UStop2(1,2)

	Co(1704) = Cd(47)*Co(1703) + MGlmT2*Co(1104)*UStop2(1,2)

	Co(1705) = Cd(9948)*Co(1702)*UStop2(1,1) + Co(1704)*UStop2(1,2)

	Co(1706) = Co(1032) + Co(1041)

	Co(1707) = MGlmT2*Co(1104) + Cd(47)*Co(1706)

	Co(1708) = AtC*Co(1707)*UCStop(1,3)**2 + At*Co(1705)*UStop2(1,1)

        Co(1709) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*Cd(10369)*Co(1071) + 
     &        2*Co(1708)*MStop2(3)))/(SB*SW2*Cd(47)*MStop2(3))

	Co(1710) = Co(1047) + Co(1049)

	Co(1712) = Co(1044) + Co(1045)

	Co(1713) = -(MGlmT2*Cd(9955)*Co(1067)) + Cd(50)*Co(1712)

	Co(1718) = Co(1034) + Co(1037)*UStop2(1,1)

	Co(1732) = Co(1058) - 2*Co(1057)*MStop2(1)

	Co(1734) = Co(1063) - 2*Co(1062)*MStop2(2)

	Co(173) = 2*CA*MTy - Co(156)*UCStopC(1,3)

	Co(1751) = 4*Cd(50)*Ce(411) + Cd(47)*Co(1147)

        Co(174) = -(Co(173)*UCStopC(1,3)) + 
     &   MUE*Co(138)*UStop2(1,1)*UStop2(1,2)

	Co(1753) = 2*Cd(9947)*Cr(85)*MGlpTmSt2(1) + Co(1071)*UStop2(1,1)

	Co(1760) = Cd(9947)*MGlpTmSt2(1) + 2*Ce(411)*UStop2(1,2)

	Co(1761) = 2*Cd(9948)*Ce(411)*UStop2(1,1) + Co(1760)*UStop2(1,2)

	Co(175) = 1/2.D0*(A0delStop(2)*Co(158)) + Co(174)*Cr(34)*MStop2(2)

        Co(1767) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*(Cd(47)*Cd(10376)*Co(1071) + 
     &          2*Cd(9993)*Co(1707)*MStop2(3)))/
     &      (SB*SW2*Cd(47)*MStop2(3)))

        Co(1768) = MGlmT2*Cd(9959)*Co(1067) + 
     &   2*Cd(50)*Co(1045)*UStop2(1,1)**2

	Co(1769) = Co(1044) + Co(1045)*UStop2(1,2)

	Co(1772) = MGlmT2*Cd(9955)*Co(1067) + Cd(50)*Co(1712)

	Co(176) = 2*CA*MTy + Co(156)*UCStopC(1,3)

	Co(1775) = Co(1037) - Co(1148)*UStop2(1,1)

        Co(177) = Co(176)*UCStopC(1,3) + 
     &   MUE*Co(138)*UStop2(1,1)*UStop2(1,2)

	Co(16) = AtC*Co(14) - At*Co(15)

        Co(178) = 1/2.D0*(MZ2*SW2*A0delStop(1)*Co(161) - 
     &      2*Co(177)*Cr(13)*MStop2(1))/(MZ2*SW2)

        Co(179) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(-(Co(175)*MStop2(1)) + 
     &        MZ2*SW2*Co(178)*MStop2(2)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2))

        Co(1806) = AtC*UCStop(1,3)*UStop2(1,1) + 
     &   At*UCStopC(1,3)*UStop2(1,2)

        Co(1811) = -(Cd(10103)*Co(1147)*Cr(75)) + 
     &   8*Cd(50)*Ce(411)*UStop2(1,1)

	Co(180) = 2*Co(174) + Co(158)*Cr(34)

	Co(1813) = Co(1147)*Cr(75) + 16*Ce(411)*UStop2(1,1)

	Co(181) = (-2*Co(177) + MZ2*SW2*Co(161)*Cr(13))/(MZ2*SW2)

        Co(1827) = At*UCStopC(1,3)*UStop2(1,1) + 
     &   AtC*UCStop(1,3)*UStop2(1,2)

        Co(182) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(-Co(180) + MZ2*SW2*Co(181)))/(SB*SW2*TB)

	Co(183) = 3*CA*(AtC*UCStop(1,3) - At*UCStopC(1,3))

	Co(184) = 3*MUE*UCStop(1,3) - 3*MUEC*UCStopC(1,3)

        Co(1858) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(50)*Ce(708)*Co(1067) - 
     &          2*Cd(9957)*Co(1099)*MStop2(3)))/
     &      (SB*SW2*Cd(50)**2*MStop2(3)))

	Co(1859) = 2*Co(1034) - Cd(9972)*Co(1037)

	Co(1860) = Cd(47)*Cd(9969)*Co(1034) + Cd(9950)*Co(1859)

        Co(1861) = At*Cd(9950)*Cd(9969)*Co(1104)*Cr(28) - 
     &   2*AtC*Co(1860)*UCStop(1,3)**2

        Co(1862) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Co(1861) + 
     &        2*Cd(9950)*Cd(9982)*Co(1140)*MStop2(3)))/
     &    (SB*SW2*Cd(47)**2*Cd(9950)*MStop2(3))

	Co(1870) = Co(1034) - 2*Co(1037)

	Co(1871) = Cd(47)*Cd(9969)*Co(1037) - Cd(9950)*Co(1870)

        Co(1872) = AtC**2*Co(1860)*UCStop(1,3)**2 + 
     &   At**2*Co(1871)*UCStopC(1,3)**2

	Co(1877) = Co(1034) + Co(1131)*UStop2(1,1)

        Co(1886) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(50)*Cd(12542)*Co(1067) + 
     &          2*Cd(9988)*Co(1099)*MStop2(3)))/
     &      (SB*SW2*Cd(50)**2*MStop2(3)))

        Co(1887) = AtC*Cd(9950)*Cd(9966)*Co(1104) + 
     &   2*At*Co(1871)*UCStopC(1,3)**2

        Co(1888) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(47)*Co(1887) + 
     &          2*Cd(9950)*Cd(9993)*Co(1140)*MStop2(3)))/
     &      (SB*SW2*Cd(47)**2*Cd(9950)*MStop2(3)))

        Co(1891) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Cd(10369)*Co(1104) + 
     &        2*Cd(9957)*Co(1140)*MStop2(3)))/
     &    (SB*SW2*Cd(47)**2*MStop2(3))

        Co(1892) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(50)*Cd(12527)*Co(1067) + 
     &        2*Cd(9982)*Co(1099)*MStop2(3)))/
     &    (SB*SW2*Cd(50)**2*MStop2(3))

        Co(1894) = AtC**2*UCStop(1,3)**2 + 
     &   At**2*Cd(10077)*UCStopC(1,3)**2

        Co(1896) = AtC**2*Cd(9969)*UCStop(1,3)**2 + 
     &   At**2*Cd(10077)*UCStopC(1,3)**2

        Co(1898) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(47)*Cd(10377)*Co(1104) + 
     &          2*Cd(9988)*Co(1140)*MStop2(3)))/
     &      (SB*SW2*Cd(47)**2*MStop2(3)))

        Co(1899) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(50)*Ce(719)*Co(1067) - 
     &        2*Cd(9993)*Co(1099)*MStop2(3)))/
     &    (SB*SW2*Cd(50)**2*MStop2(3))

	Co(1929) = CA + SA*TB

	Co(1930) = CA - SA*TB

	Co(1931) = (-2*AtC*SA + MUE*Co(1930))*MStop2(4)

	Co(1932) = 2*MTy*SA*MStop2(3) - Co(1931)*UCStop(1,3)

	Co(1933) = MUEC*Cd(435)*Co(1929) + 2*Co(1932)*UCStop(1,3)

	Co(1934) = 2*(MUEC*Cd(577)*Co(1929) + Co(1932)*UCStop(1,3))

	Co(1935) = Co(1933) - Co(1934)

	Co(1936) = -2*AtC*SA + MUE*Co(1930)

	Co(1937) = 2*MTy*SA - Co(1936)*UCStop(1,3)

        Co(1938) = -(Co(1937)*UCStop(1,3)) + 
     &   MUEC*Co(1929)*UStop2(1,1)*UStop2(1,2)

        Co(1939) = 1/2.D0*(A0delStop(1)*Co(1933)) + 
     &   Co(1938)*Cr(13)*MStop2(1)*MStop2(3)

	Co(1940) = 2*MTy*SA + Co(1936)*UCStop(1,3)

        Co(1941) = Co(1940)*UCStop(1,3) + 
     &   MUEC*Co(1929)*UStop2(1,1)*UStop2(1,2)

        Co(1942) = 1/2.D0*(A0delStop(2)*Co(1934)) - 
     &   Co(1941)*Cr(34)*MStop2(2)*MStop2(3)

        Co(1943) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(Co(1942)*MStop2(1) - Co(1939)*MStop2(2)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2)*MStop2(3))

	Co(1944) = Co(1933)*Cr(13) + 2*Co(1938)*MStop2(3)

	Co(1945) = Co(1934)*Cr(34) - 2*Co(1941)*MStop2(3)

	Co(1946) = Co(1944) - Co(1945)

	Co(1947) = -2*At*SA + MUEC*Co(1930)

        Co(1948) = 2*MUE*Cd(577)*Co(1929) - 
     &   2*Co(1947)*MStop2(4)*UCStopC(1,3)**2

        Co(1949) = (Co(1948) + 4*MTy*SA*MStop2(3)*UCStopC(1,3))/
     &   MStop2(3)

        Co(1950) = MUE*Cd(435)*Co(1929) - 
     &   2*Co(1947)*MStop2(4)*UCStopC(1,3)**2

        Co(1951) = (Co(1950) + 4*MTy*SA*MStop2(3)*UCStopC(1,3))/
     &   (MZ2*SW2*MStop2(3))

	Co(194) = CA*Cd(10987) + MUEC*SA*UUStop(3,2)

        Co(1952) = (3*cI)/(16.D0*Pi)*
     &   (Alfa2L*(-Co(1949) + MZ2*SW2*Co(1951)))/(SB*SW2*TB)

	Co(1953) = 6*SA*(AtC*UCStop(1,3) - At*UCStopC(1,3))

	Co(1954) = 6*SA*(MUE*UCStop(1,3) - MUEC*UCStopC(1,3))

	Co(1961) = 2*MTy*SA + Co(1947)*UCStopC(1,3)

        Co(1962) = Co(1961)*UCStopC(1,3) + 
     &   MUE*Co(1929)*UStop2(1,1)*UStop2(1,2)

        Co(1963) = 1/2.D0*(A0delStop(2)*Co(1949)) - 
     &   Co(1962)*Cr(34)*MStop2(2)

	Co(1964) = 2*MTy*SA - Co(1947)*UCStopC(1,3)

        Co(1965) = -(Co(1964)*UCStopC(1,3)) + 
     &   MUE*Co(1929)*UStop2(1,1)*UStop2(1,2)

        Co(1966) = 1/2.D0*(MZ2*SW2*A0delStop(1)*Co(1951) + 
     &      2*Co(1965)*Cr(13)*MStop2(1))/(MZ2*SW2)

        Co(1967) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(-(Co(1963)*MStop2(1)) + 
     &        MZ2*SW2*Co(1966)*MStop2(2)))/
     &    (SB*SW2*TB*MStop2(1)*MStop2(2))

	Co(1968) = 2*Co(1962) - Co(1949)*Cr(34)

	Co(1969) = (2*Co(1965) + MZ2*SW2*Co(1951)*Cr(13))/(MZ2*SW2)

        Co(1970) = (3*cI)/(32.D0*Pi)*
     &   (Alfa2L*(Co(1968) + MZ2*SW2*Co(1969)))/(SB*SW2*TB)

	Co(196) = 2*CA*MTy*UStop(2,1) + MUEC*Co(138)*UStop(2,2)

	Co(1974) = 3*SA*(AtC*UCStop(1,3) - At*UCStopC(1,3))

	Co(1975) = 3*SA*(MUE*UCStop(1,3) - MUEC*UCStopC(1,3))

	Co(197) = MUE*TB*UUStop(3,1) - Cd(12701)*UUStop(3,2)

	Co(1982) = SA*Cd(10987) - CA*MUEC*UUStop(3,2)

	Co(1983) = MTy*SA*Cr(91) - Cd(368)*UUStop(3,1)

	Co(1984) = 2*MTy*SA*UStop(2,1) - MUEC*Co(1929)*UStop(2,2)

	Co(1986) = 2*MTy*SA*UStop(1,1) - MUEC*Co(1929)*UStop(1,2)

	Co(2) = MUE**2*UCStop(1,3)**2 - MUEC**2*UCStopC(1,3)**2

	Co(2000) = MUEC*Co(1929)*Cr(91) - 4*MTy*SA*UUStop(3,1)

	Co(199) = 2*CA*MTy*UStop(1,1) + MUEC*Co(138)*UStop(1,2)

	Co(2037) = 2*MTy*SA*UUStop(1,1) - MUEC*Co(1929)*UUStop(1,3)

        Co(2063) = (Cd(47)**2*Cd(12096) - 
     &     MGlpTmSt2(1)*MStop2(3)*UCStopC(1,3)*UStop2(1,1))/
     &   (Cd(47)**2*MStop2(3))

	Co(2069) = At*MUE*Cd(10705) + AtC*MUEC*MStop2(1)

	Co(2070) = Cd(9948)*Cr(75) - 2*UStop2(1,1)*UStop2(1,2)**2

        Co(2077) = MUE*UCStop(1,3)**2 + 
     &   MUEC*Cd(12107)*UStop2(1,1)*UStop2(1,2)

	Co(2083) = Cr(28)*MStop2(3) - MGlpT2*UStop2(1,2)**2

	Co(2084) = MUEC*Co(2083) - 2*MUE*MStop2(1)*UCStop(1,3)**2

        Co(2085) = 2*MUE*Cr(56)*UCStop(1,3)**2 + 
     &   MUEC*Cd(12114)*UStop2(1,1)*UStop2(1,2)

	Co(19) = -(AtC*Cd(795)) + At*Cd(10219)

        Co(2091) = MUE**2*Cr(56)*UCStop(1,3)**2 + 
     &   MUEC**2*UCStopC(1,3)**2*UStop2(1,2)

	Co(2103) = Cd(10001)*Cr(28) + 2*Cd(9967)*UCStop(1,3)**2

        Co(2106) = 2*Cd(47)**2*Cd(9949)*MGlpTmSt2(2) + 
     &   Cd(50)*Cd(9951)*MStop2(3)

        Co(2107) = Cd(9950)*Co(2106) + 
     &   2*Cd(47)**2*Cd(50)*Cd(9947)*Cd(9972)*MGlpTmSt2(1)

	Co(2114) = At*Cd(9984) + 2*AtC*Cd(10081)*UCStop(1,3)**2

	Co(2117) = AtC*Cd(10080) + 2*At*Cd(9989)*UCStopC(1,3)**2

	Co(2120) = Cd(9967)*Cr(28) + 2*Cd(10001)*UCStopC(1,3)**2

        Co(2121) = Cd(47)**2*Cd(9949)*Cr(75)**2*MGlpTmSt2(2) - 
     &   Cd(50)*Cd(9951)*MStop2(3)*UStop2(1,2)**2

        Co(2122) = Cd(47)**2*Cd(9949)*Cd(9983)*MGlpTmSt2(2) - 
     &   Cd(50)*Cd(9951)*MStop2(3)*UStop2(1,1)**2

	Co(2123) = Cd(9948)*Co(2122) + Co(2121)*UStop2(1,1)

        Co(2124) = Cd(9950)*Co(2123) + 
     &   Cd(47)**2*Cd(50)*Cd(9947)*Cd(9969)*Cr(28)*MGlpTmSt2(1)

	Co(2125) = AtC*Co(2124) + At*Co(2107)*UCStopC(1,3)**2

        Co(2126) = -(Cd(9995)*UCStopC(1,3)**2) + 
     &   M_3C*UStop2(1,1)*UStop2(1,2)

	Co(2135) = At*Cd(10023)*Cr(5) + AtC*Cd(10089)*UCStop(1,3)

	Co(2138) = AtC*Cd(10093)*Cr(5) + At*Cd(10044)*UCStopC(1,3)

	Co(2143) = Cd(9950)*Co(1034) + Cd(47)*Co(1037)

        Co(2144) = MGlmT2*Cd(47)*Cd(9948)*Co(2143) - 
     &   Cd(9950)*Co(1058)*UStop2(1,2)

	Co(2145) = Cd(10103)*Co(1047) + Cd(50)*Co(1049)

        Co(2146) = MGlmT2*Cd(50)*Cd(9948)*Co(2145) - 
     &   Cd(10103)*Co(1063)*UStop2(1,2)

	Co(2147) = 2*Cd(50)*Ce(411) + Co(1062)

	Co(213) = MUE*TB*UUStop(1,1) - Cd(12701)*UUStop(1,2)

	Co(2153) = Cd(50) + Cd(10103)

        Co(2154) = MGlmT2*Cd(50)*Cd(9948)*Co(1049)*Co(2153) - 
     &   Cd(10103)*Co(1063)*UStop2(1,2)

	Co(2158) = MGlmT2*Cd(50)*Co(1049) + Co(1063)

	Co(2159) = MGlmT2*Cd(50)**2*Co(1047) + Cd(10103)*Co(2158)

	Co(2163) = Co(1056) - 2*Co(1042)*UStop2(1,1)

	Co(2164) = Co(1043) - Co(1042)*UStop2(1,2)

	Co(2165) = Cd(9948)*Co(2163) + 2*Co(2164)*UStop2(1,2)

	Co(2168) = Co(1061) - 2*Co(1050)*UStop2(1,1)

	Co(2169) = Co(1051) - Co(1050)*UStop2(1,2)

	Co(2170) = Cd(9948)*Co(2168) + 2*Co(2169)*UStop2(1,2)

	Co(2171) = Cd(9983)*Co(1049) + Co(1047)*UStop2(1,1)**2

	Co(2179) = Co(1047) - Cd(9969)*Co(1049)

	Co(2184) = Co(1034) - Cd(9972)*Co(1037)

        Co(20) = ((TB*Co(2) - Co(19))*MStop2(3)**2*UStop2(1,1)*
     &     UStop2(1,2))/TB

	Co(2198) = Co(1034)*UStop2(1,1) + Co(1037)*UStop2(1,2)

	Co(218) = 2*CA*MTy*UUStop(1,1) + MUEC*Co(138)*UUStop(1,3)

	Co(2212) = Co(1037)*UStop2(1,1) + Co(1034)*UStop2(1,2)

	Co(2227) = MGlmT2*Cd(47)*Co(1034) + Co(1058)

	Co(2228) = MGlmT2*Cd(47)**2*Co(1034) + Cd(9950)*Co(2227)

        Co(2229) = At*Co(1174)*UCStopC(1,3)**2 + 
     &   AtC*Cd(9948)*Co(2228)*UStop2(1,1)

	Co(2230) = MGlmT2*Cd(50)*Co(1047) + Co(1063)

	Co(2231) = MGlmT2*Cd(50)**2*Co(1047) + Cd(10103)*Co(2230)

	Co(2232) = 2*Co(1047) - Co(1049)

	Co(2233) = Co(1063) - MGlmT2*Cd(50)*Co(2232)

	Co(2234) = MGlmT2*Cd(50)**2*Co(1049) - Cd(10103)*Co(2233)

        Co(2235) = At*Co(2234)*UCStopC(1,3)**2 + 
     &   AtC*Cd(9948)*Co(2231)*UStop2(1,1)

	Co(2243) = Co(1034) - Cd(9969)*Co(1037)

	Co(2250) = Cd(9983)*Co(1047) + Co(1049)*UStop2(1,1)**2

	Co(2257) = -(Cd(9972)*Co(1034)) + Co(1037)

	Co(2258) = AtC*Co(2243)*Cr(28) - 2*At*Co(2257)*UCStopC(1,3)**2

	Co(2261) = Co(1047) - Cd(9972)*Co(1049)

	Co(226) = CA*MTy*Cr(91) + Cd(5)*UUStop(3,1)

	Co(2287) = Co(1042) - Co(1043)*UStop2(1,2)

	Co(2289) = Co(1050) - Co(1051)*UStop2(1,2)

	Co(227) = MUEC*Co(138)*Cr(91) + 4*CA*MTy*UUStop(3,1)

	Co(2297) = Co(1034) - Co(1037)*UStop2(1,2)

	Co(2298) = Cd(47)*Co(1706) + MGlmT2*Co(2297)

        Co(2299) = MGlmT2*Cd(9948)*Co(1037)*UStop2(1,1) - 
     &   Co(2298)*UStop2(1,2)

        Co(2300) = -(AtC*Co(1707)*UCStop(1,3)**2) + 
     &   At*Co(2299)*UStop2(1,1)

        Co(2301) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(47)*Cd(10369)*Co(1071) - 
     &        2*Co(2300)*MStop2(3)))/(SB*SW2*Cd(47)*MStop2(3))

	Co(2302) = At*Co(2243)*Cr(28) - 2*AtC*Co(2184)*UCStop(1,3)**2

	Co(2313) = 4*Ce(411) - Cd(9972)*Co(1071)

        Co(2319) = -(MGlmT2*Cd(9958)*Co(1049)) + 
     &   2*Cd(50)*Co(1712)*UStop2(1,2)

	Co(2321) = -(MGlmT2*Cd(9955)*Co(1047)) + Cd(50)*Co(1712)

	Co(2328) = 4*Cd(10103)*Ce(411) + Cd(50)*Co(1710)

        Co(2335) = AtC**2*Co(2179)*UCStop(1,3)**2 - 
     &   At**2*Co(2261)*UCStopC(1,3)**2

	Co(2347) = 4*Cd(9959)*Ce(411) + Co(1710)*UStop2(1,1)**2

	Co(2348) = Cd(50)*Co(2347) - 2*Co(1050)*UStop2(1,1)**2

        Co(2360) = MGlmT2*Cd(9958)*Co(1047) + 
     &   2*Cd(50)*Co(1712)*UStop2(1,2)

	Co(2362) = MGlmT2*Cd(9955)*Co(1047) + Cd(50)*Co(1712)

	Co(2388) = AtC*Co(2212)*UCStop(1,3) + At*Co(2198)*UCStopC(1,3)

	Co(2391) = -(Cd(10563)*Co(1047)) + Cd(10037)*Co(1049)

	Co(2394) = AtC*Co(1049)*UCStop(1,3) + At*Co(1047)*UCStopC(1,3)

	Co(2409) = AtC*Co(1047)*UCStop(1,3) + At*Co(1049)*UCStopC(1,3)

	Co(2411) = AtC*Co(2198)*UCStop(1,3) + At*Co(2212)*UCStopC(1,3)

	Co(2430) = 4*Cd(6096)*Cd(9950)*Ce(411) + Cd(47)*Cd(50)*Co(1071)

	Co(2435) = Cd(9948)*Co(1051)*UStop2(1,1) - Co(2289)*UStop2(1,2)

        Co(2436) = -(AtC*Co(1099)*UCStop(1,3)**2) + 
     &   At*Co(2435)*UStop2(1,1)

	Co(2437) = Co(1047) - 2*Cd(9969)*Co(1049)

	Co(2438) = Cd(50)*Cd(9969)*Co(1049) - Cd(10103)*Co(2437)

	Co(2441) = Cd(9948)*Co(1043)*UStop2(1,1) - Co(2287)*UStop2(1,2)

        Co(2442) = -(AtC*Co(1140)*UCStop(1,3)**2) + 
     &   At*Co(2441)*UStop2(1,1)

        Co(2443) = At*Cd(9950)*Co(2243)*Cr(28) - 
     &   2*AtC*Co(1860)*UCStop(1,3)**2

        Co(2444) = cI/(64.D0*Pi**3)*
     &   (Alfas2L*(MGlmT2*Cd(47)*Co(2443) - 
     &        2*Cd(9950)*Co(2442)*MStop2(3)))/
     &    (SB*SW2*Cd(47)**2*Cd(9950)*MStop2(3))

	Co(2451) = -(Cd(9972)*Co(1034)) + 2*Co(1037)

	Co(2452) = Cd(47)*Cd(9969)*Co(1037) + Cd(9950)*Co(2451)

        Co(2453) = AtC**2*Co(1860)*UCStop(1,3)**2 + 
     &   At**2*Co(2452)*UCStopC(1,3)**2

	Co(2455) = 2*Co(1047) - Cd(9972)*Co(1049)

	Co(2456) = Cd(50)*Cd(9969)*Co(1047) + Cd(10103)*Co(2455)

        Co(2457) = AtC**2*Co(2438)*UCStop(1,3)**2 + 
     &   At**2*Co(2456)*UCStopC(1,3)**2

	Co(2464) = Cd(9948)*Co(1042)*UStop2(1,1) - Co(2164)*UStop2(1,2)

        Co(2465) = At*Co(1140)*UCStopC(1,3)**2 + 
     &   AtC*Co(2464)*UStop2(1,1)

	Co(2468) = Cd(9948)*Co(1050)*UStop2(1,1) - Co(2169)*UStop2(1,2)

        Co(2469) = At*Co(1099)*UCStopC(1,3)**2 + 
     &   AtC*Co(2468)*UStop2(1,1)

	Co(2474) = AtC*Co(1099)*UCStop(1,3)**2 + At*Co(2468)*UStop2(1,1)

	Co(2478) = AtC*Co(1140)*UCStop(1,3)**2 + At*Co(2464)*UStop2(1,1)

        Co(2482) = AtC**2*Co(2257)*UCStop(1,3)**2 + 
     &   At**2*Co(2184)*UCStopC(1,3)**2

	Co(2484) = -(Cd(9972)*Co(1047)) + Co(1049)

        Co(2487) = -(At*Co(1140)*UCStopC(1,3)**2) + 
     &   AtC*Co(2441)*UStop2(1,1)

	Co(2488) = AtC*Co(2243)*Cr(28) - 2*At*Co(2184)*UCStopC(1,3)**2

        Co(2489) = -(cI/(64.D0*Pi**3)*
     &     (Alfas2L*(MGlmT2*Cd(47)*Co(2488) - 
     &          2*Co(2487)*MStop2(3)))/(SB*SW2*Cd(47)**2*MStop2(3))
     &     )

        Co(2490) = -(At*Co(1099)*UCStopC(1,3)**2) + 
     &   AtC*Co(2435)*UStop2(1,1)

	Co(2496) = 1 - Log(MStop2(2)/Q)**2

	Co(2497) = 2*A0delStop(2) - Co(2496)*MStop2(2)

	Co(2498) = 1 - Cr(13)

	Co(2499) = A0delStop(1) - Cr(13)*MStop2(1)

	Co(2500) = Co(2499) - Co(2498)*Cr(34)*MStop2(1)

	Co(2501) = 2*(A0delStop(2)*MStop2(1) + Co(2500)*MStop2(2))

	Co(2503) = 1 - Log(MStop2(1)/Q)**2

	Co(2504) = 2*A0delStop(1) - Co(2503)*MStop2(1)

	Co(2505) = AtC*Cd(154)*UCStop(1,3) + At*Cd(158)*UCStopC(1,3)

	Co(2508) = 1 - Cr(34)

	Co(2510) = A0delStop(2) - Cr(34)*MStop2(2)

	Co(2512) = M_3C*UCStop(1,3) - M_3*UCStopC(1,3)

	Co(2513) = MUE*Cd(76) - MUEC*Cd(107)

	Co(2514) = -(AtC*Cd(76)) + At*Cd(107)

	Co(2516) = AtC*Cd(238)*UCStop(1,3) - At*Cd(240)*UCStopC(1,3)

	Co(2518) = At*Cd(329) - AtC*Cd(335)

	Co(2525) = 2*(-MGl2 + Cr(65))

	Co(24) = -(Cr(28)*MStop2(1)) + 2*Cd(535)*UStop2(1,1)*UStop2(1,2)

        Co(259) = (2*Cd(47)**2*Cd(9973) + 
     &     Cd(50)*Cd(9950)*Cd(9951)*MStop2(3))/
     &   (Cd(47)**2*Cd(50)*Cd(9950)*MStop2(3))

	Co(261) = -(AtC*UCStop(1,3)**2) + At*Cd(9981)*UStop2(1,1)

	Co(262) = At*Cd(10071) - 2*AtC*Ce(368)*UCStop(1,3)**2

        Co(264) = -(M_3C*UCStop(1,3)**2) + 
     &   M_3*UStop2(1,1)*UStop2(1,2)

        Co(265) = -(Cd(9963)*UCStop(1,3)**2) + 
     &   M_3*UStop2(1,1)*UStop2(1,2)

        Co(266) = Cd(50)*Cd(9947)*Co(265)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Co(264)*MGlpTmSt2(2)

	Co(267) = At*Cd(10118) + 2*Co(266)

	Co(269) = M_3*Cr(28) + Cd(9967)*UCStop(1,3)**2

	Co(270) = 2*Co(269) + At*Cd(9969)*Cr(28)

	Co(274) = -(AtC*UCStop(1,3)**2) + At*UStop2(1,1)*UStop2(1,2)

	Co(275) = At*Cd(10080) + 2*AtC*Cd(9985)*UCStop(1,3)**2

	Co(277) = -(At*UCStopC(1,3)**2) + AtC*Cd(9981)*UStop2(1,1)

	Co(278) = AtC*Ce(372) - 2*At*Cd(10125)*UCStopC(1,3)**2

	Co(281) = -(At*UCStopC(1,3)**2) + AtC*UStop2(1,1)*UStop2(1,2)

	Co(282) = AtC*Cd(9959) + 2*M_3C*UStop2(1,1)

        Co(286) = AtC*Cd(9971) - 
     &   2*At*Cd(10076)*Cd(10077)*UCStopC(1,3)**2

	Co(289) = M_3C*Cr(28) + Cd(10001)*UCStopC(1,3)**2

	Co(290) = 2*Co(289) + AtC*Cd(9969)*Cr(28)

	Co(294) = -(M_3C*At) + M_3*AtC

	Co(296) = AtC*Cd(10032) - At*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Co(297) = At*Cd(10093)*Cr(5) + AtC*Cd(10022)*UCStop(1,3)

        Co(303) = -(At*Cd(10039)*UCStopC(1,3)) + 
     &   AtC*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Co(304) = AtC*Cd(10023)*Cr(5) + At*Cd(10097)*UCStopC(1,3)

	Co(312) = -(Cd(1272)*MStop2(1)) + MGlmT2*MStop2(2)

	Co(315) = MGlmT2 - MStop2(1)*UStop2(1,1)

	Co(316) = Cd(47) - Co(315)*MStop2(1)

	Co(317) = MGlmT2 - Cd(1281)*UStop2(1,1)

	Co(318) = Cd(47) - Co(317)*MStop2(1)

	Co(328) = At*Cd(10533)*Cr(28) + 2*AtC*Cd(10445)*UCStop(1,3)**2

	Co(329) = AtC*Cd(10443)*Cr(28) - 2*At*Cd(10534)*UCStopC(1,3)**2

        Co(331) = Cd(47)*Cd(152)*Cd(10058) + 
     &   2*MUEC*Cd(164)*MStop2(1)*MStop2(3)*UStop2(1,1)**2

        Co(333) = MUE**2*Cd(9972)*UCStop(1,3)**2 - 
     &   MUEC**2*UCStopC(1,3)**2

        Co(334) = Cd(47)*Cd(152)*Co(333) + 
     &   Cd(118)*Cd(164)*Co(2)*MStop2(3)

	Co(336) = -3*MUE*UCStop(1,3)**2 + MUEC*Cr(75)*UStop2(1,1)

        Co(337) = M_3C*Co(336)*UCStop(1,3) + 
     &   M_3*Cd(13)*Cr(75)*UStop2(1,1)

	Co(338) = M_3*Cd(13)*Cr(28) + M_3C*Ce(642)*UCStop(1,3)

	Co(339) = MGlmT2*Co(337) - Co(338)*MStop2(1)

	Co(341) = -3*MUEC*UCStopC(1,3)**2 + MUE*Cr(75)*UStop2(1,1)

        Co(342) = M_3*Co(341)*UCStopC(1,3) + 
     &   M_3C*Cd(13)*Cr(75)*UStop2(1,1)

	Co(343) = MGlmT2*Ce(532) - Co(342)*MStop2(1)

        Co(345) = M_3*MUE**2*Cd(10453)*UCStop(1,3) + 
     &   M_3C*MUEC**2*Cd(10303)*UCStopC(1,3)

	Co(347) = 3*MUE**2*UCStop(1,3)**3 + MUEC**2*Cr(28)*UCStopC(1,3)

        Co(348) = 3*MUEC**2*UCStopC(1,3)**3 + 
     &   MUE**2*Cr(75)*UCStop(1,3)*UStop2(1,1)

	Co(349) = M_3C*Co(347) - M_3*Co(348)

        Co(350) = 3*MUE**2*UCStop(1,3)**3 + 
     &   MUEC**2*Cr(75)*UCStopC(1,3)*UStop2(1,1)

	Co(351) = MUE**2*Cr(28)*UCStop(1,3) + 3*MUEC**2*UCStopC(1,3)**3

	Co(352) = -(M_3C*Co(350)) + M_3*Co(351)

	Co(353) = MGlmT2*Co(349) + Co(352)*MStop2(1)

	Co(360) = -2*AtC*UCStop(1,3)**2 + 2*At*Cd(9981)*UStop2(1,1)

	Co(361) = At*Cd(9969)*Cr(28) + 2*AtC*UCStop(1,3)**2

	Co(363) = At*Cd(9966) + 2*AtC*Cd(9972)*UCStop(1,3)**2

	Co(364) = At*Cr(28) + 2*AtC*UCStop(1,3)**2

	Co(368) = -2*At*UCStopC(1,3)**2 + 2*AtC*UStop2(1,1)*UStop2(1,2)

	Co(369) = AtC*Cd(9969)*Cr(28) + 2*At*UCStopC(1,3)**2

	Co(371) = AtC*Cd(9966) + 2*At*Cd(9972)*UCStopC(1,3)**2

	Co(372) = AtC*Cr(28) + 2*At*UCStopC(1,3)**2

	Co(373) = 4*MGlmT2*Cd(47)*Cd(9995) + At*Cd(10279)

        Co(377) = -2*Cd(47)*Cd(152)*Cd(10058)*MStop2(1) + 
     &   MUEC*Cd(10279)*MStop2(3)*UStop2(1,1)**2

        Co(379) = MUE**2*Cd(10279)*UCStop(1,3)**2 - 
     &   2*MUEC**2*Cd(10280)*UCStopC(1,3)**2

	Co(380) = 4*Cd(47)*Cd(152)*Co(333)*MStop2(1) + Co(379)*MStop2(3)

        Co(381) = 2*Cd(47)*Cd(152)*Cd(10057)*MStop2(1) - 
     &   Cd(10396)*MStop2(3)

	Co(384) = MGlmT2*Cd(10400) - Ce(517)*MStop2(1)

        Co(386) = MUE**2*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(10077)*UCStopC(1,3)**2

	Co(37) = AtC*UUStop(1,1) + 2*At*UUStop(1,2)

	Co(387) = MGlmT2*Co(386) - Co(333)*MStop2(1)

	Co(388) = 4*Cd(47)*Co(387)*MStop2(1) + Co(379)*MStop2(3)

	Co(389) = -(MGlmT2*Ce(530)) + Cd(10057)*MStop2(1)

	Co(4) = AtC**2*UCStop(1,3)**2 - At**2*UCStopC(1,3)**2

	Co(38) = Co(37) + 2*MUEC*TB*UUStop(1,2)

	Co(402) = 2*Cd(47)**2 + Cd(12371)*MStop2(1)

	Co(403) = 4*Cd(47)*Cd(11825)*MStop2(1) + Co(402)*MStop2(3)

	Co(404) = Cd(9955)*MStop2(1) - MStop2(2)

	Co(39) = -(MUEC*Co(38)) + At*MUE*UUStop(1,1)

	Co(407) = Cd(47) + Cd(164)*MStop2(1)*UStop2(1,1)

	Co(413) = AtC*MStop2(2) + M_3C*MStop2(3)

	Co(416) = Cd(47) + Cd(164)*MStop2(1)

	Co(40) = 2*MUE*TB*UUStop(1,1) + At*UUStop(1,2)

	Co(417) = Cd(47) + 2*Cd(164)*MStop2(1)

	Co(418) = 2*Cd(47)*Co(417) + Cd(12371)*MStop2(1)

	Co(420) = -(At*Co(404)) + M_3*MStop2(3)

        Co(421) = -(Co(420)*UCStopC(1,3)**2) + 
     &   M_3C*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Co(424) = 2*MGlmT2*Cd(47)*Cr(28)*MStop2(1) - 
     &   Cd(12371)*MStop2(3)*UStop2(1,1)**2

	Co(41) = AtC*Cd(565) + MUE*Co(40)

	Co(427) = AtC*Cd(12471) + At*Ce(1104)*UCStopC(1,3)**2

        Co(429) = 4*Cd(47)*Cd(164)*Co(333)*MStop2(1) + 
     &   Cd(12371)*Co(2)*MStop2(3)

        Co(431) = 4*Cd(47)*Cd(164)*Cd(10018)*MStop2(1) - 
     &   Cd(12371)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Co(432) = At*Co(431)*Cr(5) + AtC*Cd(12404)*UCStop(1,3)

	Co(433) = AtC*Cd(12475)*Cr(5) + At*Ce(1117)*UCStopC(1,3)

	Co(42) = Co(39)*UUStopC(1,1) + Co(41)*UUStopC(1,2)

	Co(437) = Cd(47)*Cd(9972) - Cd(164)*Co(404)

        Co(438) = 8*MUE**2*Cd(11825)*UCStop(1,3)**2 - 
     &   8*MUEC**2*Co(437)*UCStopC(1,3)**2

	Co(439) = Cd(47) + Cd(164)*Cd(1281)*UStop2(1,1)

	Co(440) = AtC*Cd(11825) + M_3C*Cd(164)*MStop2(3)

	Co(444) = At*Cd(11790)*Cr(28) + 2*AtC*Cd(12472)*UCStop(1,3)**2

	Co(445) = AtC*Cd(12470)*Cr(28) + 2*At*Ce(1103)*UCStopC(1,3)**2

	Co(43) = MUE*TB*UUStop(1,1) + 2*At*UUStop(1,2)

        Co(450) = -(Cd(47)*Cd(164)*Co(333)) + 
     &   Cd(152)*Co(2)*MStop2(1)*MStop2(3)

	Co(452) = At*Cd(10018)*Cr(5) + AtC*Cd(10088)*UCStop(1,3)

	Co(453) = AtC*Cd(10091)*Cr(5) + At*Cd(10042)*UCStopC(1,3)

	Co(44) = TB*(AtC*UUStop(1,1) - 2*At*UUStop(1,2))

        Co(465) = 8*MUE**2*Ce(1398)*UCStop(1,3)**2 - 
     &   8*MUEC**2*Cd(11944)*UCStopC(1,3)**2

	Co(45) = At*Co(43) - MUEC*Co(44)

	Co(466) = AtC*Cd(9969)*MStop2(1) - M_3C*MStop2(3)

        Co(467) = Co(466)*UCStop(1,3)**2 + 
     &   M_3*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Co(472) = Cd(50) + Cd(251)*Cr(85)*MStop2(1)

	Co(473) = At*Co(472) - M_3*Cd(251)*MStop2(3)

        Co(477) = Cd(50)*Cd(251)*Cd(9983) - 
     &   2*Cd(179)*MStop2(2)*MStop2(3)*UStop2(1,1)**2

	Co(479) = Cd(50)*Cd(251) + Cd(179)*MStop2(2)*MStop2(3)

	Co(481) = At*Cd(10092)*Cr(5) - AtC*Cd(10021)*UCStop(1,3)

	Co(482) = AtC*Cd(10020)*Cr(5) + At*Cd(10095)*UCStopC(1,3)

	Co(484) = MGlmT2*Ce(670) + Cd(10674)*MStop2(2)

	Co(485) = -(MGlmT2*Ce(669)) + MUEC*Cd(9983)*MStop2(2)

        Co(486) = Cd(50)*Co(485) + 
     &   2*MGl2*MUEC*Cd(179)*MStop2(3)*UStop2(1,1)**2

        Co(488) = MUE**2*Ce(797)*UCStop(1,3)**2 - 
     &   MUEC**2*Ce(686)*UCStopC(1,3)**2

	Co(490) = MGlmT2*Cd(10697) - MUE*MStop2(2)*UStop2(1,1)**2

	Co(491) = MGlmT2*Cd(10696) + Ce(672)*MStop2(2)

	Co(500) = MStop2(2) + MStop2(3)*UStop2(1,1)

	Co(501) = MStop2(2)*UStop2(1,1)**2 + Co(500)*UStop2(1,2)**2

	Co(508) = Cd(10658)*MStop2(2) + MGlmT2*Cr(75)**2*UStop2(1,1)

	Co(511) = MGlmT2*Cd(9972) + MStop2(2)

	Co(512) = 4*Cd(50)*Co(511)*MStop2(2) - Cd(12483)*MStop2(3)

	Co(513) = MGlmT2*Cd(10658) + Cr(75)**2*MStop2(2)*UStop2(1,1)

	Co(514) = Cd(9983)*MStop2(2) - MGlmT2*UStop2(1,1)**2

        Co(515) = 2*Cd(50)*Co(514)*MStop2(2) + 
     &   Cd(12483)*MStop2(3)*UStop2(1,1)**2

	Co(519) = 4*Cd(50)*Cd(251)*MStop2(2) - Cd(12483)*MStop2(3)

        Co(522) = 4*Cd(50)*Cd(251)*MStop2(2)*UStop2(1,1) + 
     &   Cd(12483)*MStop2(3)*UStop2(1,2)

        Co(523) = 4*Cd(50)*Cd(251)*Cd(9948)*Cr(75)*MStop2(2) - 
     &   Co(522)*UStop2(1,1)

	Co(524) = At*Co(523)*Cr(5) + AtC*Ce(1364)*UCStop(1,3)

        Co(526) = 4*Cd(50)*Cd(251)*Cd(10020)*MStop2(2) + 
     &   Cd(12483)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Co(536) = -2*AtC*UCStop(1,3)**2 + 2*At*UStop2(1,1)*UStop2(1,2)

	Co(538) = At*Cd(10070) + 2*AtC*Cd(9969)*UCStop(1,3)**2

	Co(539) = MGlmT2*Cd(50)*Cd(9958) + Cd(10604)*UStop2(1,2)**2

	Co(540) = M_3C + AtC*UStop2(1,2)

	Co(541) = AtC*Cd(10603) + 4*MGlmT2*Cd(50)*Co(540)

	Co(542) = At*Cd(10609) + 2*AtC*Cd(10611)*UCStop(1,3)**2

	Co(547) = MGlmT2*Cd(10745) + Cd(10077)*MStop2(1)*MStop2(2)

        Co(549) = 2*Cd(10382)*Cd(10604) + 
     &   MUEC*Cd(10603)*UStop2(1,1)**2*UStop2(1,2)

	Co(551) = 2*Cd(50)*Cd(179)*MStop2(2) - Cd(10603)*MStop2(3)

        Co(553) = 2*MUE**2*Cd(10604)*UCStop(1,3)**2 - 
     &   MUEC**2*Cd(10603)*UCStopC(1,3)**2

	Co(554) = 4*Cd(50)*Cd(179)*Co(2)*MStop2(2) - Co(553)*MStop2(3)

	Co(556) = -2*At*UCStopC(1,3)**2 + 2*AtC*Cd(9981)*UStop2(1,1)

	Co(557) = 2*AtC*Cd(9970) + 4*At*Cd(10077)*UCStopC(1,3)**2

	Co(558) = 2*Cd(10604)*UStop2(1,1) + Cd(10603)*UStop2(1,2)

	Co(559) = 2*MGlmT2*Cd(50) + Co(558)*UStop2(1,1)*UStop2(1,2)

	Co(560) = AtC*Ce(1261) + 2*At*Cd(10611)*UCStopC(1,3)**2

	Co(561) = M_3 + At*Cd(9955)

        Co(562) = -(Co(561)*UCStopC(1,3)**2) + 
     &   M_3C*UStop2(1,1)*UStop2(1,2)

	Co(565) = MGlmT2*Ce(670) - Cd(10674)*MStop2(2)

	Co(55) = MUEC*Cd(577) + MUE*MStop2(4)*UCStop(1,3)**2

	Co(566) = MGlmT2*Ce(669) + MUEC*Cd(9983)*MStop2(2)

        Co(567) = Cd(50)*Co(566)*MStop2(2) - 
     &   MUEC*Cd(10604)*MStop2(3)*UStop2(1,1)**2

	Co(573) = MGlmT2*Cd(10697) + MUE*MStop2(2)*UStop2(1,1)**2

	Co(574) = -(MGlmT2*Cd(10696)) + Ce(672)*MStop2(2)

	Co(56) = MUE*Cd(577) + MUEC*MStop2(4)*UCStopC(1,3)**2

	Co(577) = -(MGlmT2*Cd(10759)) + Cd(10021)*MStop2(2)

	Co(580) = Cd(10020)*MStop2(2) - MGlmT2*UStop2(1,1)**2

	Co(57) = AtC*Co(55) - At*Co(56)

	Co(592) = At*Cd(9983) + 2*AtC*UCStop(1,3)**2

	Co(597) = Cd(50) + Cr(85)*MStop2(1)*MStop2(2)

	Co(598) = -(At*Co(597)) + M_3*MStop2(2)*MStop2(3)

	Co(599) = At*Cd(10816) - M_3*MStop2(3)

        Co(600) = Co(599)*UCStopC(1,3)**2 + 
     &   M_3C*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Co(601) = MGlmT2*Co(600) + Co(598)*UCStopC(1,3)**2

	Co(605) = MGl2*Cd(427) + 2*M_3*MTy*MUE**2*Cr(75)*UCStop(1,3)

        Co(606) = MGlmT2*Co(605) + 
     &   2*M_3C*MTy*MUEC**2*Cr(75)*MStop2(2)*UCStopC(1,3)

        Co(607) = 3*MUEC**2*UCStopC(1,3)**3 + 
     &   MUE**2*UCStop(1,3)*UStop2(1,1)**2

        Co(609) = MUE**2*Cd(10020)*UCStop(1,3) + 
     &   3*MUEC**2*UCStopC(1,3)**3

	Co(616) = MTy2*Cd(10674) + 4*MTy*Cd(10836)

	Co(617) = MTy2*Cd(10674) + 4*MTy*Cd(10838)

	Co(619) = MTy2*MUEC*Cd(9983) - 4*MTy*Cd(10834)*Cr(75)

	Co(624) = 4*MTy*Cd(10842) + MTy2*Ce(672)

	Co(625) = 4*MTy*Cd(10844) + MTy2*Ce(672)

	Co(627) = 4*MTy*Cd(10840)*Cr(75) + MTy2*MUE*UStop2(1,1)**2

        Co(62) = MUE**2*Cd(520)*UCStop(1,3)**2 - 
     &   4*MUEC**2*UCStopC(1,3)**2*UStop2(1,2)

	Co(635) = AtC*MUEC*Cd(9948) - At*MUE*UStop2(1,2)

	Co(636) = At*MUE*Cd(9948) - AtC*MUEC*UStop2(1,2)

	Co(637) = ((MGlpT2*Co(635) + Co(636)*MStop2(1))*UStop(1,1))/TB

        Co(638) = -(MUEC*Cd(12701)*UCStopC(1,3)**3) + 
     &   MUE*Cd(12700)*UCStop(1,3)*UStop2(1,1)*UStop2(1,2)

        Co(639) = (Cd(13)*Cd(47)**2*Cd(12836)*MGlpTmSt2(2) - 
     &     TB*Cd(50)**2*Co(637)*UStopC(1,1))/
     &   (TB*Cd(47)**2*Cd(50)**2)

        Co(640) = (Co(638)*MGlpTmSt2(1) - 
     &     TB*Cd(47)**2*Co(639)*UCStopC(1,3))/(TB*Cd(47)**2)

	Co(641) = At*MUE - AtC*MUEC

	Co(642) = (Co(641)*MGlpTmSt2(1)*UStop(1,2)*UStop2(1,1))/TB

        Co(643) = -(MUE*Cd(12700)*UCStop(1,3)**2) + 
     &   Co(641)*UStop2(1,1)*UStop2(1,2)

        Co(644) = Co(643)*UStop(2,1) - 
     &   MUEC*Cd(9948)*Cd(12701)*UCStopC(1,3)*UStop(2,2)

	Co(63) = 4*MUEC*UCStopC(1,3)**2 + MUE*Cd(520)*UStop2(1,1)

        Co(645) = MUE*Cd(12700)*UCStop(1,3)**3 - 
     &   MUEC*Cd(12701)*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

	Co(647) = AtC*UUStop(1,1) - At*UUStop(1,2)

	Co(649) = MUE*TB*UUStop(1,1) - At*UUStop(1,2)

	Co(650) = AtC*Cd(12148) + MUE*Co(649)

	Co(654) = MStop2(2) - MTy2*UStop2(1,2)

	Co(64) = AtC*Cd(520)*Cd(795) - At*Co(63)*UStop2(1,2)

        Co(655) = -(MUEC*Cd(48)*UCStopC(1,3)**2) + 
     &   MUE*Co(654)*UStop2(1,1)*UStop2(1,2)

	Co(661) = AtC*Cr(91)*UStop(1,1) - 2*At*UStop(2,2)*UUStop(1,2)

	Co(65) = ((TB*Co(62) + Co(64))*MStop2(3)**2*UStop2(1,1))/TB

        Co(665) = At*Cd(11019)*UStop(1,2) - 
     &   4*MUE*TB*UStop(2,1)*UUStop(1,1)

        Co(666) = MUEC*Cd(11019)*UStop(1,2) + 
     &   4*MUE*UStop(2,1)*UUStop(1,1)

	Co(667) = -(MUE*Co(665)) + AtC*Co(666)

	Co(668) = AtC*Cd(11161)*UStop(1,1) - 4*At*UStop(2,2)*UUStop(1,2)

	Co(671) = At*Cr(91)*UStop(2,2) - 2*MUE*TB*UStop(1,1)*UUStop(2,1)

	Co(672) = MUEC*Cr(91)*UStop(2,2) + 2*MUE*UStop(1,1)*UUStop(2,1)

	Co(673) = -(MUE*Co(671)) + AtC*Co(672)

	Co(674) = AtC*Cd(11019)*UStop(2,1) - 4*At*UStop(1,2)*UUStop(2,2)

	Co(679) = AtC*MGlpT2*MUEC + At*MUE*MStop2(1)

	Co(680) = MGlpT2 - 5*MStop2(1)

	Co(681) = 5*MGlpT2 - MStop2(1)

	Co(682) = At*MUE*Co(680) - AtC*MUEC*Co(681)

	Co(67) = 2 + 6*UStop2(1,2)

        Co(687) = Co(641)*UCStopC(1,3)*UStop(1,2) + 
     &   MUE*Cd(12700)*Cr(91)*UStopC(2,2)

        Co(689) = MUE*MStop2(3)*UStop(1,1) - 
     &   2*MTy*Co(641)*Cr(91)*UStopC(2,1)

        Co(690) = 2*MTy*MUE*Cd(12700)*Cr(91)*UCStop(3,3) + 
     &   Co(689)*UStop(1,2)

	Co(695) = 2*Cd(47)*Cd(9948) - MGlpTmSt2(1)*MStop2(3)*UStop2(1,1)

	Co(696) = AtC*MUEC*Cd(12088) + At*MUE*Co(695)

        Co(697) = MTy*(6*MUE*Cd(47)*Cd(12700)*UCStop(1,3)*
     &      UStop(1,1) + Co(696)*UStop(1,2))

        Co(698) = (Co(697) + MUE*Cd(47)*MStop2(3)*UStop(1,1))/
     &   (TB*MStop2(3))

	Co(699) = At*UCStopC(1,3)*UStop(1,2) + AtC*Cr(91)*UStopC(2,2)

	Co(700) = MStop2(3)*UStop(1,1) + 2*At*MTy*Cr(91)*UStopC(2,1)

	Co(701) = 2*AtC*MTy*Cr(91)*UCStop(3,3) - Co(700)*UStop(1,2)

	Co(703) = MUE*Cd(50) + MTy*Co(641)*MGlpTmSt2(2)*UCStopC(1,3)

        Co(704) = (Co(703)*MStop2(3)*UStop(2,1) - 
     &     2*MTy*Cd(50)*Co(641)*UStop(2,2)*UStop2(1,1))/MStop2(3)

        Co(69) = -4*Cd(634)*MStop2(2)*UStop2(1,1) + 
     &   2*Cd(632)*UStop2(1,1)**3

        Co(707) = MTy*MUE*Cd(12700)*MGlpTmSt2(2)*
     &    UCStop(1,3)**3 + Co(704)*UCStop(1,3)*UStopC(2,2)

        Co(709) = 3*Cd(12700)*UCStop(1,3)*UStop(1,1) + 
     &   At*Cr(91)*UStopC(2,1)

        Co(710) = (-(Cd(47)*Co(709)*UStop(2,1)) + 
     &     At*Cd(50)*Cr(91)*UStop2(1,1))/(Cd(47)*Cd(50))

        Co(712) = (Cd(12700)*MGlpTmSt2(1)*MStop2(3)*
     &      UCStop(1,3)**3 + 2*Cd(47)**2*Co(710)*UUStopC(3,2))/
     &   (Cd(47)**2*MStop2(3))

	Co(713) = (-(MUE*Co(712)) + MUEC**2*TB*Cd(12097)*UStop2(1,2))/TB

	Co(718) = 2*MUE*UCStop(1,3)**2 + MUEC*UStop2(1,2)**2

        Co(71) = MUE**2*Cd(639)*MStop2(3)*UCStop(1,3)**2 - 
     &   2*MUEC**2*Cd(642)*UCStopC(1,3)**2

	Co(727) = MStop2(3) + MStop2(1)*UStop2(1,2)

        Co(734) = -(MUE*UCStop(1,3)**2) + 
     &   MUEC*Cd(12125)*UStop2(1,1)*UStop2(1,2)

        Co(739) = MUEC*UCStopC(1,3)**2 + 
     &   MUE*Cd(12107)*UStop2(1,1)*UStop2(1,2)

        Co(740) = -(MUEC*UCStopC(1,3)**2) + 
     &   MUE*Cd(12125)*UStop2(1,1)*UStop2(1,2)

	Co(741) = MTy2*Co(739) + Co(740)*MStop2(2)

	Co(748) = Cd(10979)*UStop2(1,1)**2 - 2*MUE*UUStop(1,1)

	Co(749) = MUE*Co(43)*UUStop(1,1) - MUEC*Cd(12701)*UUStop(1,2)**2

	Co(751) = 2*MUE*Cd(12700)*UCStop(3,3) - Co(641)*UCStop(3,4)

	Co(752) = Cd(12700)*UUStop(3,1) + 3*At*UUStop(3,2)

        Co(753) = At*Cr(91)*UStop(1,2) + 
     &   4*Cd(12700)*UStop(2,1)*UUStop(1,1)

	Co(754) = At*Cr(91)*UCStop(3,1) - 2*Co(752)*UCStop(3,3)

	Co(74) = 2*UStop2(1,1)**2 + Cr(85)*UStop2(1,2)

	Co(755) = Co(754)*UStop2(1,2) + Co(753)*UStop2(1,1)*UStopC(2,2)

	Co(756) = 2*Co(755)*MStop2(1)**2 - 6*At*MStop2(3)**2*UUStop(1,3)

        Co(757) = AtC*Cd(10938) - 
     &   4*At*Cr(7)*MStop2(1)**2*UCStopC(1,3)*UUStop(1,2)

        Co(759) = Co(757) - 
     &   4*MUEC*TB*Cr(7)*MStop2(1)**2*UCStopC(1,3)*UUStop(1,2)

	Co(760) = MUE*Co(756) - MUEC*Co(759)

	Co(761) = 12*Co(641) - Co(39)*UUStopC(1,1)

	Co(764) = Cd(12700)*UUStop(2,1) + At*UUStop(2,2)

	Co(75) = -(Co(74)*MStop2(2)) + Cr(26)*UStop2(1,1)*UStop2(1,2)

	Co(765) = -(Co(764)*UCStop(3,3)) + At*Cr(91)*UStop2(1,2)

        Co(766) = At*Cd(10952)*MStop2(3) - 
     &   2*Co(765)*MStop2(1)*UStopC(1,1)

	Co(767) = -3*UStop(2,2) + UCStopC(3,1)*UStop(1,2)*UStop2(1,1)

        Co(768) = 2*At*Co(767)*MStop2(3) - 
     &   2*Co(765)*MStop2(1)*UStopC(1,1)

	Co(769) = AtC*UCStop(3,3) - 2*At*UCStop(3,4)

	Co(770) = 2*At*Cd(10956) - AtC*UCStop(1,3)*UStop(1,1)

        Co(771) = -(At*Cd(10922)*MStop2(3)*UStop(1,2)) + 
     &   2*Co(770)*MStop2(1)*UStop2(1,1)

	Co(772) = Co(771)*UCStopC(1,3) - 2*Co(769)*MStop2(1)*UStop(2,2)

	Co(773) = Cd(10961) + 5*MStop2(3)**2*UUStop(1,2)

	Co(774) = At*Co(773) + 2*AtC*Cd(10951)*MStop2(1)**2*UCStop(3,3)

        Co(775) = 3*Co(770)*MStop2(1) + 
     &   4*At*MStop2(3)*UStop(1,2)*UStop2(1,1)

        Co(776) = -5*Co(769)*MStop2(1)*UStop(2,2) + 
     &   Co(775)*UCStopC(1,3)*UStop2(1,1)

	Co(777) = 2*Co(776)*MStop2(3) - 2*Co(774)*UStopC(1,1)

        Co(778) = At*Cd(10961) + 
     &   2*AtC*Cd(10951)*MStop2(1)**2*UCStop(3,3)

        Co(779) = Co(770)*MStop2(1) + 
     &   At*MStop2(3)*UStop(1,2)*UStop2(1,1)

        Co(780) = -3*Co(769)*MStop2(1)*UStop(2,2) + 
     &   Co(779)*UCStopC(1,3)*UStop2(1,1)

	Co(781) = 4*Co(780)*MStop2(3) - 2*Co(778)*UStopC(1,1)

	Co(788) = 2*MUE*TB*UUStop(3,1) + At*UUStop(3,2)

	Co(789) = AtC*Cd(11011) + MUE*Co(788)

        Co(790) = Co(789)*UStopC(2,1) - 
     &   2*MUE*Cd(12700)*Cr(91)*UStopC(2,2)

	Co(791) = 2*Co(790)*MStop2(1) - 2*Co(41)*MStop2(3)*UStopC(1,1)

	Co(792) = At*Cr(91)*UStop(1,2) + MUE*TB*UStop(2,1)*UUStop(1,1)

	Co(793) = MUEC*Cr(91)*UStop(1,2) - MUE*UStop(2,1)*UUStop(1,1)

	Co(794) = -(MUE*Co(792)) + AtC*Co(793)

        Co(78) = AtC**2*Cd(651)*MStop2(3)*UCStop(1,3)**2 - 
     &   At**2*Cd(655)*UCStopC(1,3)**2

        Co(795) = Co(794)*UStop(1,1) + 
     &   MUEC*Cd(12701)*UUStop(1,2)*UUStop(3,2)

	Co(796) = AtC*Cd(11106)*UUStop(1,1) + 2*At*Cd(10925)*UUStop(1,2)

	Co(797) = Co(796) + 2*MUEC*TB*Cd(10925)*UUStop(1,2)

	Co(798) = MUEC*Co(797) - At*MUE*Cd(11106)*UUStop(1,1)

	Co(799) = Co(798)*UStopC(2,1) - 2*Co(795)*UStopC(2,2)

        Co(800) = Co(39)*MStop2(3)*UStop2(1,1)*UStop2(1,2) + 
     &   Co(799)*MStop2(1)*UStopC(1,2)

        Co(801) = 3*Co(751)*MStop2(1)*UCStopC(3,4) - 
     &   Co(800)*UUStopC(1,1)

	Co(802) = 4*Co(801)*MStop2(3) - 2*Co(760)*UUStopC(1,3)

	Co(804) = At*Cr(91)*UStop(2,2) + 4*MUE*TB*UStop(1,1)*UUStop(2,1)

	Co(79) = AtC*MUEC*Cd(649) + At*MUE*Co(75)

        Co(805) = MUE*Cd(14040)*UUStop(2,1) - 
     &   MUEC*Cd(12701)*UUStop(1,2)*UUStop(2,2)

        Co(806) = MUE*Co(804)*UCStop(3,1)*UStopC(2,1) + 
     &   Co(805)*UUStopC(2,1)

	Co(807) = MUE**2*TB*Cd(11024) - 2*Co(806)*MStop2(1)**2

        Co(808) = MUE*Co(792)*UStop(1,1) - 
     &   MUEC*Cd(12701)*UUStop(1,2)*UUStop(3,2)

        Co(809) = 2*MUE*TB*Cd(11106)*UUStop(1,1) + 
     &   At*Cd(10925)*UUStop(1,2)

	Co(810) = 2*Co(808)*UUStopC(3,1) - MUE*Co(809)*UUStopC(3,3)

        Co(813) = MUEC*Cd(10020)*MStop2(1)*UUStop(1,2) - 
     &   Cd(11169)*MStop2(3)*UUStopC(3,3)

	Co(818) = 2*Co(790)*MStop2(1) - 4*Co(41)*MStop2(3)*UStopC(1,1)

        Co(819) = 2*Co(39)*MStop2(3)*UStop2(1,1)*UStop2(1,2) + 
     &   Co(799)*MStop2(1)*UStopC(1,2)

	Co(820) = Co(751)*MStop2(1)*UCStopC(3,4) + Co(819)*UUStopC(1,1)

	Co(822) = Co(777) - Co(781)

	Co(823) = Co(766) - Co(768)

	Co(825) = 2*Cd(10979)*UStop2(1,1)**2 - 5*MUE*UUStop(1,1)

	Co(826) = Co(825)*MStop2(3) + Cd(11169)*MStop2(1)*UUStopC(3,3)

	Co(827) = 5*MUE**2*TB*UUStop(1,1) - 2*Co(749)*UUStopC(1,1)

	Co(834) = 2*MUE*TB*UUStop(2,1) + At*UUStop(2,2)

	Co(835) = 2*MUE*UUStop(2,1) - MUEC*UUStop(2,2)

	Co(836) = MUE*Co(834) + AtC*Co(835)

        Co(837) = Co(789)*MStop2(2)*UStopC(1,1) + 
     &   Co(836)*MStop2(3)*UStopC(2,1)

	Co(838) = 2*MUE*Cd(970) + 4*MUEC*UCStopC(1,3)**2

	Co(839) = MUEC*Cr(91)*UCStopC(3,3) - MUE*UStop2(1,1)*UUStop(2,1)

	Co(843) = 2*Cd(9948)*MStop2(3) + MStop2(2)*UStop2(1,1)

	Co(847) = AtC*Cd(11023) - MUE*Co(804)

	Co(848) = AtC*Cd(10946) + MUE*Cd(14040)

        Co(849) = Co(848)*UUStop(2,1) - 
     &   MUEC*Cd(12701)*UUStop(1,2)*UUStop(2,2)

	Co(851) = 3*UStop(1,2) - UCStop(3,1)*UStop(2,2)*UStop2(1,2)

        Co(852) = 2*Co(851)*MStop2(3) - 
     &   2*Cr(91)*MStop2(2)*UCStop(3,1)*UStopC(1,1)

        Co(853) = AtC*Co(852) + 
     &   4*At*MStop2(4)*UCStopC(1,3)**2*UStop(1,2)

        Co(854) = Co(853) + 
     &   4*MUEC*TB*MStop2(4)*UCStopC(1,3)**2*UStop(1,2)

	Co(84) = 2*(MStop2(1) + MStop2(4)*UStop2(1,2))

	Co(855) = At*MUE*Co(852) - MUEC*Co(854)

	Co(856) = At*Cr(91)*UStop(1,2) + 4*MUE*TB*UStop(2,1)*UUStop(1,1)

	Co(857) = AtC*Cd(11405) - MUE*Co(856)

        Co(858) = 2*MUE*TB*Cd(10925)*UUStop(2,1) + 
     &   At*Cd(11106)*UUStop(2,2)

        Co(859) = 2*MUE*Cd(10925)*UUStop(2,1) - 
     &   MUEC*Cd(11106)*UUStop(2,2)

	Co(860) = MUE*Co(858) + AtC*Co(859)

        Co(861) = 2*Co(857)*MStop2(2)*UCStopC(3,3) - 
     &   2*Co(860)*MStop2(3)*UStopC(2,1)

	end
