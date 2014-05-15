#include "model.h"
#include "util.h"
#include "looptools.h"
#include "renconst.h"

	double precision S, T, U
	common /gzkinvars/ S, T, U

	integer Hel(4)
	common /gzkinvars/ Hel

	double complex F9, F11, F1, F3, F10, F5, F7, F12, F4, F15, F13
	double complex F8, F16, F14, F2, F6, Pair4, Pair5, Pair3
	double complex Pair1, Pair2, Eps1, Eps2, Abb29, Abb30, Abb21
	double complex Abb22, Abb35, Abb36, Abb23, Abb24, Abb1, Abb41
	double complex Abb42, Abb4, Abb2, Abb43, Abb44, Abb5, Abb13
	double complex Abb16, Abb3, Abb6, Abb31, Abb17, Abb37, Abb14
	double complex Abb32, Abb18, Abb38, Abb15, Abb8, Abb10, Abb11
	double complex Abb7, Abb19, Abb20, Abb9, Abb12, Abb33, Abb25
	double complex Abb39, Abb26, Abb34, Abb27, Abb40, Abb28
	double complex AbbSum777, AbbSum393, AbbSum159, AbbSum144
	double complex AbbSum533, AbbSum114, AbbSum147, AbbSum364
	double complex AbbSum177, AbbSum7, AbbSum327, AbbSum313
	double complex AbbSum64, AbbSum406, AbbSum160, AbbSum158
	double complex AbbSum146, AbbSum168, AbbSum545, AbbSum116
	double complex AbbSum112, AbbSum148, AbbSum145, AbbSum69
	double complex AbbSum386, AbbSum310, AbbSum179, AbbSum175
	double complex AbbSum8, AbbSum316, AbbSum348, AbbSum317
	double complex AbbSum453, AbbSum585, AbbSum189, AbbSum336
	double complex AbbSum431, AbbSum558, AbbSum591, AbbSum586
	double complex AbbSum61, AbbSum60, AbbSum100, AbbSum247
	double complex AbbSum208, AbbSum120, AbbSum880, AbbSum134
	double complex AbbSum118, AbbSum874, AbbSum132, AbbSum207
	double complex AbbSum99, AbbSum246, AbbSum846, AbbSum915
	double complex AbbSum845, AbbSum66, AbbSum150, AbbSum206
	double complex AbbSum97, AbbSum151, AbbSum262, AbbSum149
	double complex AbbSum104, AbbSum245, AbbSum595, AbbSum896
	double complex AbbSum751, AbbSum244, AbbSum680, AbbSum806
	double complex AbbSum230, AbbSum186, AbbSum111, AbbSum212
	double complex AbbSum137, AbbSum174, AbbSum893, AbbSum606
	double complex AbbSum254, AbbSum844, AbbSum259, AbbSum129
	double complex AbbSum960, AbbSum713, AbbSum471, AbbSum203
	double complex AbbSum82, AbbSum383, AbbSum200, AbbSum659
	double complex AbbSum385, AbbSum446, AbbSum430, AbbSum1036
	double complex AbbSum297, AbbSum1035, AbbSum172, AbbSum658
	double complex AbbSum96, AbbSum608, AbbSum226, AbbSum195
	double complex AbbSum63, AbbSum759, AbbSum588, AbbSum1095
	double complex AbbSum565, AbbSum243, AbbSum665, AbbSum895
	double complex AbbSum242, AbbSum218, AbbSum135, AbbSum228
	double complex AbbSum106, AbbSum201, AbbSum610, AbbSum136
	double complex AbbSum225, AbbSum173, AbbSum889, AbbSum593
	double complex AbbSum1094, AbbSum220, AbbSum204, AbbSum261
	double complex AbbSum253, AbbSum216, AbbSum224, AbbSum167
	double complex AbbSum1060, AbbSum1061, AbbSum258, AbbSum154
	double complex AbbSum843, AbbSum1058, AbbSum124, AbbSum941
	double complex AbbSum913, AbbSum1059, AbbSum103, AbbSum618
	double complex AbbSum914, AbbSum704, AbbSum470, AbbSum202
	double complex AbbSum81, AbbSum360, AbbSum197, AbbSum643
	double complex AbbSum362, AbbSum434, AbbSum413, AbbSum185
	double complex AbbSum110, AbbSum1032, AbbSum296, AbbSum1030
	double complex AbbSum171, AbbSum642, AbbSum93, AbbSum596
	double complex AbbSum805, AbbSum750, AbbSum255, AbbSum193
	double complex AbbSum62, AbbSum211, AbbSum758, AbbSum577
	double complex AbbSum1004, AbbSum1033, AbbSum1055, AbbSum1108
	double complex AbbSum918, AbbSum39, AbbSum40, AbbSum1003
	double complex AbbSum1031, AbbSum1054, AbbSum141, AbbSum221
	double complex AbbSum67, AbbSum156, AbbSum363, AbbSum108
	double complex AbbSum157, AbbSum80, AbbSum184, AbbSum155
	double complex AbbSum138, AbbSum182, AbbSum414, AbbSum410
	double complex AbbSum92, AbbSum130, AbbSum433, AbbSum79
	double complex AbbSum205, AbbSum181, AbbSum209, AbbSum458
	double complex AbbSum152, AbbSum566, AbbSum140, AbbSum567
	double complex AbbSum290, AbbSum819, AbbSum900, AbbSum1214
	double complex AbbSum654, AbbSum526, AbbSum271, AbbSum334
	double complex AbbSum771, AbbSum115, AbbSum907, AbbSum1007
	double complex AbbSum178, AbbSum249, AbbSum861, AbbSum101
	double complex AbbSum916, AbbSum1080, AbbSum1096, AbbSum229
	double complex AbbSum44, AbbSum1062, AbbSum213, AbbSum1056
	double complex AbbSum50, AbbSum222, AbbSum142, AbbSum902
	double complex AbbSum1005, AbbSum1109, AbbSum256, AbbSum1064
	double complex AbbSum1088, AbbSum65, AbbSum117, AbbSum113
	double complex AbbSum908, AbbSum449, AbbSum1008, AbbSum180
	double complex AbbSum176, AbbSum109, AbbSum250, AbbSum248
	double complex AbbSum429, AbbSum105, AbbSum867, AbbSum102
	double complex AbbSum98, AbbSum917, AbbSum153, AbbSum575
	double complex AbbSum616, AbbSum1083, AbbSum235, AbbSum605
	double complex AbbSum745, AbbSum492, AbbSum71, AbbSum927
	double complex AbbSum864, AbbSum882, AbbSum1052, AbbSum274
	double complex AbbSum928, AbbSum1097, AbbSum231, AbbSum107
	double complex AbbSum48, AbbSum1063, AbbSum217, AbbSum214
	double complex AbbSum210, AbbSum1057, AbbSum215, AbbSum68
	double complex AbbSum52, AbbSum128, AbbSum223, AbbSum219
	double complex AbbSum564, AbbSum656, AbbSum143, AbbSum139
	double complex AbbSum619, AbbSum572, AbbSum906, AbbSum1006
	double complex AbbSum183, AbbSum462, AbbSum496, AbbSum1110
	double complex AbbSum252, AbbSum428, AbbSum556, AbbSum574
	double complex AbbSum573, AbbSum257, AbbSum251, AbbSum1065
	double complex AbbSum227, AbbSum1093, AbbSum240, AbbSum587
	double complex AbbSum849, AbbSum740, AbbSum162, AbbSum84
	double complex AbbSum967, AbbSum1114, AbbSum88, AbbSum1130
	double complex AbbSum922, AbbSum1116, AbbSum1066, AbbSum761
	double complex AbbSum754, AbbSum741, AbbSum1140, AbbSum1151
	double complex AbbSum972, AbbSum1143, AbbSum72, AbbSum1184
	double complex AbbSum978, AbbSum979, AbbSum294, AbbSum1153
	double complex AbbSum1199, AbbSum980, AbbSum292, AbbSum919
	double complex AbbSum956, AbbSum945, AbbSum119, AbbSum755
	double complex AbbSum742, AbbSum746, AbbSum163, AbbSum87
	double complex AbbSum970, AbbSum1115, AbbSum90, AbbSum1131
	double complex AbbSum923, AbbSum1117, AbbSum1074, AbbSum786
	double complex AbbSum942, AbbSum121, AbbSum161, AbbSum75
	double complex AbbSum78, AbbSum1044, AbbSum852, AbbSum981
	double complex AbbSum1147, AbbSum1159, AbbSum973, AbbSum1148
	double complex AbbSum467, AbbSum73, AbbSum1193, AbbSum920
	double complex AbbSum988, AbbSum295, AbbSum1160, AbbSum1207
	double complex AbbSum989, AbbSum885, AbbSum887, AbbSum70
	double complex AbbSum122, AbbSum853, AbbSum1045, AbbSum89
	double complex AbbSum293, AbbSum847, AbbSum961, AbbSum83
	double complex AbbSum86, AbbSum329, AbbSum764, AbbSum996
	double complex AbbSum436, AbbSum692, AbbSum938, AbbSum507
	double complex AbbSum1132, AbbSum812, AbbSum752, AbbSum1009
	double complex AbbSum597, AbbSum626, AbbSum169, AbbSum578
	double complex AbbSum1119, AbbSum716, AbbSum765, AbbSum1103
	double complex AbbSum1010, AbbSum832, AbbSum396, AbbSum476
	double complex AbbSum1038, AbbSum1024, AbbSum872, AbbSum934
	double complex AbbSum909, AbbSum1142, AbbSum535, AbbSum328
	double complex AbbSum666, AbbSum474, AbbSum504, AbbSum691
	double complex AbbSum506, AbbSum718, AbbSum532, AbbSum394
	double complex AbbSum549, AbbSum831, AbbSum435, AbbSum1186
	double complex AbbSum813, AbbSum361, AbbSum1098, AbbSum762
	double complex AbbSum579, AbbSum1011, AbbSum662, AbbSum594
	double complex AbbSum623, AbbSum641, AbbSum703, AbbSum1015
	double complex AbbSum944, AbbSum1201, AbbSum1187, AbbSum1154
	double complex AbbSum1072, AbbSum624, AbbSum969, AbbSum190
	double complex AbbSum1126, AbbSum1022, AbbSum1073, AbbSum1127
	double complex AbbSum1016, AbbSum315, AbbSum547, AbbSum1111
	double complex AbbSum966, AbbSum53, AbbSum350, AbbSum788
	double complex AbbSum1001, AbbSum448, AbbSum700, AbbSum530
	double complex AbbSum1135, AbbSum822, AbbSum830, AbbSum935
	double complex AbbSum1185, AbbSum1200, AbbSum1018, AbbSum943
	double complex AbbSum609, AbbSum639, AbbSum170, AbbSum589
	double complex AbbSum1128, AbbSum1068, AbbSum731, AbbSum305
	double complex AbbSum326, AbbSum717, AbbSum534, AbbSum508
	double complex AbbSum627, AbbSum123, AbbSum661, AbbSum811
	double complex AbbSum995, AbbSum997, AbbSum814, AbbSum839
	double complex AbbSum408, AbbSum500, AbbSum1141, AbbSum1152
	double complex AbbSum395, AbbSum531, AbbSum897, AbbSum412
	double complex AbbSum576, AbbSum715, AbbSum452, AbbSum898
	double complex AbbSum660, AbbSum546, AbbSum349, AbbSum910
	double complex AbbSum681, AbbSum499, AbbSum528, AbbSum699
	double complex AbbSum529, AbbSum505, AbbSum690, AbbSum392
	double complex AbbSum732, AbbSum544, AbbSum407, AbbSum555
	double complex AbbSum939, AbbSum753, AbbSum447, AbbSum1194
	double complex AbbSum823, AbbSum384, AbbSum1039, AbbSum1102
	double complex AbbSum1099, AbbSum787, AbbSum590, AbbSum663
	double complex AbbSum607, AbbSum637, AbbSum657, AbbSum712
	double complex AbbSum1120, AbbSum763, AbbSum1021, AbbSum1025
	double complex AbbSum1067, AbbSum873, AbbSum930, AbbSum304
	double complex AbbSum625, AbbSum702, AbbSum475, AbbSum411
	double complex AbbSum638, AbbSum325, AbbSum971, AbbSum192
	double complex AbbSum1129, AbbSum1017, AbbSum1023, AbbSum308
	double complex AbbSum309, AbbSum188, AbbSum548, AbbSum640
	double complex AbbSum557, AbbSum359, AbbSum454, AbbSum489
	double complex AbbSum343, AbbSum1012, AbbSum697, AbbSum644
	double complex AbbSum485, AbbSum388, AbbSum390, AbbSum288
	double complex AbbSum1213, AbbSum278, AbbSum1137, AbbSum298
	double complex AbbSum523, AbbSum419, AbbSum191, AbbSum409
	double complex AbbSum32, AbbSum1076, AbbSum187, AbbSum553
	double complex AbbSum352, AbbSum424, AbbSum486, AbbSum420
	double complex AbbSum456, AbbSum306, AbbSum1107, AbbSum311
	double complex AbbSum377, AbbSum387, AbbSum332, AbbSum866
	double complex AbbSum696, AbbSum710, AbbSum381, AbbSum291
	double complex AbbSum792, AbbSum267, AbbSum289, AbbSum714
	double complex AbbSum826, AbbSum1029, AbbSum127, AbbSum51
	double complex AbbSum333, AbbSum990, AbbSum840, AbbSum904
	double complex AbbSum1113, AbbSum318, AbbSum701, AbbSum353
	double complex AbbSum469, AbbSum1037, AbbSum416, AbbSum307
	double complex AbbSum354, AbbSum952, AbbSum391, AbbSum37
	double complex AbbSum241, AbbSum1215, AbbSum793, AbbSum865
	double complex AbbSum421, AbbSum342, AbbSum59, AbbSum477
	double complex AbbSum478, AbbSum272, AbbSum389, AbbSum270
	double complex AbbSum351, AbbSum35, AbbSum335, AbbSum427
	double complex AbbSum422, AbbSum460, AbbSum312, AbbSum791
	double complex AbbSum415, AbbSum423, AbbSum459, AbbSum881
	double complex AbbSum1043, AbbSum319, AbbSum863, AbbSum688
	double complex AbbSum33, AbbSum689, AbbSum323, AbbSum264
	double complex AbbSum324, AbbSum695, AbbSum491, AbbSum368
	double complex AbbSum404, AbbSum631, AbbSum773, AbbSum344
	double complex AbbSum648, AbbSum512, AbbSum522, AbbSum671
	double complex AbbSum398, AbbSum331, AbbSum542, AbbSum345
	double complex AbbSum543, AbbSum403, AbbSum513, AbbSum537
	double complex AbbSum632, AbbSum466, AbbSum378, AbbSum31
	double complex AbbSum650, AbbSum633, AbbSum266, AbbSum774
	double complex AbbSum775, AbbSum958, AbbSum959, AbbSum635
	double complex AbbSum698, AbbSum636, AbbSum379, AbbSum655
	double complex AbbSum279, AbbSum495, AbbSum265, AbbSum653
	double complex AbbSum924, AbbSum1040, AbbSum417, AbbSum855
	double complex AbbSum559, AbbSum560, AbbSum570, AbbSum856
	double complex AbbSum372, AbbSum598, AbbSum337, AbbSum479
	double complex AbbSum517, AbbSum738, AbbSum330, AbbSum488
	double complex AbbSum490, AbbSum480, AbbSum487, AbbSum501
	double complex AbbSum705, AbbSum707, AbbSum369, AbbSum514
	double complex AbbSum876, AbbSum287, AbbSum1070, AbbSum49
	double complex AbbSum984, AbbSum857, AbbSum58, AbbSum126
	double complex AbbSum286, AbbSum442, AbbSum418, AbbSum373
	double complex AbbSum461, AbbSum951, AbbSum164, AbbSum1212
	double complex AbbSum706, AbbSum694, AbbSum339, AbbSum365
	double complex AbbSum735, AbbSum370, AbbSum371, AbbSum515
	double complex AbbSum338, AbbSum518, AbbSum645, AbbSum509
	double complex AbbSum516, AbbSum482, AbbSum340, AbbSum481
	double complex AbbSum708, AbbSum277, AbbSum269, AbbSum341
	double complex AbbSum709, AbbSum693, AbbSum1048, AbbSum833
	double complex AbbSum366, AbbSum367, AbbSum646, AbbSum510
	double complex AbbSum483, AbbSum649, AbbSum571, AbbSum613
	double complex AbbSum647, AbbSum511, AbbSum374, AbbSum561
	double complex AbbSum1112, AbbSum562, AbbSum260, AbbSum273
	double complex AbbSum749, AbbSum94, AbbSum76, AbbSum95
	double complex AbbSum1134, AbbSum859, AbbSum1026, AbbSum1034
	double complex AbbSum877, AbbSum380, AbbSum375, AbbSum34
	double complex AbbSum875, AbbSum784, AbbSum268, AbbSum711
	double complex AbbSum563, AbbSum652, AbbSum198, AbbSum878
	double complex AbbSum785, AbbSum196, AbbSum199, AbbSum879
	double complex AbbSum860, AbbSum165, AbbSum425, AbbSum614
	double complex AbbSum376, AbbSum275, AbbSum444, AbbSum445
	double complex AbbSum426, AbbSum125, AbbSum493, AbbSum36
	double complex AbbSum494, AbbSum525, AbbSum276, AbbSum670
	double complex AbbSum622, AbbSum651, AbbSum925, AbbSum634
	double complex AbbSum615, AbbSum1106, AbbSum484, AbbSum985
	double complex AbbSum54, AbbSum1166, AbbSum1169, AbbSum1014
	double complex AbbSum1164, AbbSum1189, AbbSum899, AbbSum56
	double complex AbbSum744, AbbSum1220, AbbSum1218, AbbSum1078
	double complex AbbSum1202, AbbSum1174, AbbSum1047, AbbSum836
	double complex AbbSum281, AbbSum9, AbbSum19, AbbSum779
	double complex AbbSum1, AbbSum1216, AbbSum15, AbbSum10
	double complex AbbSum17, AbbSum767, AbbSum1013, AbbSum43
	double complex AbbSum780, AbbSum768, AbbSum1177, AbbSum911
	double complex AbbSum1210, AbbSum835, AbbSum756, AbbSum770
	double complex AbbSum800, AbbSum1133, AbbSum1190, AbbSum801
	double complex AbbSum1191, AbbSum778, AbbSum766, AbbSum931
	double complex AbbSum41, AbbSum1176, AbbSum975, AbbSum883
	double complex AbbSum42, AbbSum1049, AbbSum1046, AbbSum950
	double complex AbbSum743, AbbSum999, AbbSum948, AbbSum949
	double complex AbbSum946, AbbSum947, AbbSum1155, AbbSum1086
	double complex AbbSum1087, AbbSum854, AbbSum1069, AbbSum817
	double complex AbbSum890, AbbSum1071, AbbSum848, AbbSum1085
	double complex AbbSum1192, AbbSum818, AbbSum1203, AbbSum13
	double complex AbbSum1206, AbbSum729, AbbSum1100, AbbSum21
	double complex AbbSum1084, AbbSum1146, AbbSum987, AbbSum901
	double complex AbbSum4, AbbSum1079, AbbSum936, AbbSum1162
	double complex AbbSum22, AbbSum29, AbbSum23, AbbSum797
	double complex AbbSum27, AbbSum827, AbbSum664, AbbSum472
	double complex AbbSum3, AbbSum828, AbbSum356, AbbSum473
	double complex AbbSum921, AbbSum810, AbbSum321, AbbSum503
	double complex AbbSum781, AbbSum986, AbbSum322, AbbSum991
	double complex AbbSum674, AbbSum687, AbbSum55, AbbSum727
	double complex AbbSum1165, AbbSum1020, AbbSum903, AbbSum1167
	double complex AbbSum1171, AbbSum580, AbbSum673, AbbSum284
	double complex AbbSum722, AbbSum283, AbbSum684, AbbSum1221
	double complex AbbSum57, AbbSum748, AbbSum1219, AbbSum1081
	double complex AbbSum686, AbbSum685, AbbSum302, AbbSum232
	double complex AbbSum437, AbbSum438, AbbSum955, AbbSum782
	double complex AbbSum783, AbbSum798, AbbSum954, AbbSum1175
	double complex AbbSum1051, AbbSum1188, AbbSum629, AbbSum723
	double complex AbbSum724, AbbSum539, AbbSum683, AbbSum842
	double complex AbbSum1136, AbbSum282, AbbSum796, AbbSum794
	double complex AbbSum976, AbbSum46, AbbSum11, AbbSum790
	double complex AbbSum1195, AbbSum912, AbbSum20, AbbSum795
	double complex AbbSum47, AbbSum1211, AbbSum888, AbbSum1053
	double complex AbbSum2, AbbSum841, AbbSum804, AbbSum789
	double complex AbbSum1050, AbbSum965, AbbSum1217, AbbSum1196
	double complex AbbSum933, AbbSum45, AbbSum16, AbbSum757
	double complex AbbSum1179, AbbSum1019, AbbSum12, AbbSum18
	double complex AbbSum1122, AbbSum675, AbbSum628, AbbSum280
	double complex AbbSum402, AbbSum397, AbbSum401, AbbSum672
	double complex AbbSum405, AbbSum725, AbbSum728, AbbSum540
	double complex AbbSum301, AbbSum554, AbbSum399, AbbSum667
	double complex AbbSum721, AbbSum726, AbbSum551, AbbSum581
	double complex AbbSum737, AbbSum734, AbbSum550, AbbSum1121
	double complex AbbSum285, AbbSum233, AbbSum300, AbbSum630
	double complex AbbSum974, AbbSum263, AbbSum91, AbbSum599
	double complex AbbSum668, AbbSum669, AbbSum166, AbbSum739
	double complex AbbSum131, AbbSum464, AbbSum463, AbbSum519
	double complex AbbSum620, AbbSum521, AbbSum520, AbbSum465
	double complex AbbSum621, AbbSum747, AbbSum676, AbbSum1002
	double complex AbbSum584, AbbSum541, AbbSum538, AbbSum677
	double complex AbbSum678, AbbSum600, AbbSum964, AbbSum962
	double complex AbbSum963, AbbSum1161, AbbSum1091, AbbSum582
	double complex AbbSum455, AbbSum604, AbbSum238, AbbSum439
	double complex AbbSum440, AbbSum601, AbbSum602, AbbSum603
	double complex AbbSum552, AbbSum400, AbbSum1092, AbbSum457
	double complex AbbSum239, AbbSum432, AbbSum568, AbbSum592
	double complex AbbSum569, AbbSum441, AbbSum862, AbbSum720
	double complex AbbSum1075, AbbSum583, AbbSum824, AbbSum1077
	double complex AbbSum1197, AbbSum850, AbbSum825, AbbSum1208
	double complex AbbSum1090, AbbSum776, AbbSum894, AbbSum799
	double complex AbbSum957, AbbSum1178, AbbSum772, AbbSum803
	double complex AbbSum837, AbbSum769, AbbSum802, AbbSum1156
	double complex AbbSum617, AbbSum133, AbbSum468, AbbSum237
	double complex AbbSum1000, AbbSum998, AbbSum1157, AbbSum14
	double complex AbbSum1101, AbbSum24, AbbSum992, AbbSum6
	double complex AbbSum937, AbbSum25, AbbSum30, AbbSum26
	double complex AbbSum28, AbbSum1209, AbbSum1082, AbbSum1149
	double complex AbbSum1163, AbbSum733, AbbSum905, AbbSum1089
	double complex AbbSum299, AbbSum1125, AbbSum682, AbbSum85
	double complex AbbSum1124, AbbSum730, AbbSum736, AbbSum443
	double complex AbbSum236, AbbSum612, AbbSum524, AbbSum194
	double complex AbbSum234, AbbSum536, AbbSum968, AbbSum1204
	double complex AbbSum982, AbbSum1144, AbbSum816, AbbSum815
	double complex AbbSum1145, AbbSum834, AbbSum719, AbbSum611
	double complex AbbSum983, AbbSum1158, AbbSum891, AbbSum1205
	double complex AbbSum884, AbbSum1027, AbbSum1170, AbbSum838
	double complex AbbSum320, AbbSum1139, AbbSum993, AbbSum1105
	double complex AbbSum953, AbbSum679, AbbSum1138, AbbSum303
	double complex AbbSum932, AbbSum820, AbbSum497, AbbSum355
	double complex AbbSum5, AbbSum829, AbbSum994, AbbSum1118
	double complex AbbSum74, AbbSum1028, AbbSum382, AbbSum450
	double complex AbbSum886, AbbSum858, AbbSum498, AbbSum357
	double complex AbbSum892, AbbSum926, AbbSum1104, AbbSum1150
	double complex AbbSum807, AbbSum1180, AbbSum871, AbbSum358
	double complex AbbSum821, AbbSum1198, AbbSum870, AbbSum502
	double complex AbbSum940, AbbSum1041, AbbSum346, AbbSum1182
	double complex AbbSum527, AbbSum977, AbbSum1168, AbbSum1172
	double complex AbbSum809, AbbSum929, AbbSum851, AbbSum314
	double complex AbbSum347, AbbSum38, AbbSum869, AbbSum1173
	double complex AbbSum77, AbbSum1123, AbbSum451, AbbSum760
	double complex AbbSum1181, AbbSum868, AbbSum1183, AbbSum1042
	double complex AbbSum808
	common /gzabbrev/ F9, F11, F1, F3, F10, F5, F7, F12, F4, F15
	common /gzabbrev/ F13, F8, F16, F14, F2, F6, Pair4, Pair5
	common /gzabbrev/ Pair3, Pair1, Pair2, Eps1, Eps2, Abb29
	common /gzabbrev/ Abb30, Abb21, Abb22, Abb35, Abb36, Abb23
	common /gzabbrev/ Abb24, Abb1, Abb41, Abb42, Abb4, Abb2, Abb43
	common /gzabbrev/ Abb44, Abb5, Abb13, Abb16, Abb3, Abb6, Abb31
	common /gzabbrev/ Abb17, Abb37, Abb14, Abb32, Abb18, Abb38
	common /gzabbrev/ Abb15, Abb8, Abb10, Abb11, Abb7, Abb19
	common /gzabbrev/ Abb20, Abb9, Abb12, Abb33, Abb25, Abb39
	common /gzabbrev/ Abb26, Abb34, Abb27, Abb40, Abb28, AbbSum777
	common /gzabbrev/ AbbSum393, AbbSum159, AbbSum144, AbbSum533
	common /gzabbrev/ AbbSum114, AbbSum147, AbbSum364, AbbSum177
	common /gzabbrev/ AbbSum7, AbbSum327, AbbSum313, AbbSum64
	common /gzabbrev/ AbbSum406, AbbSum160, AbbSum158, AbbSum146
	common /gzabbrev/ AbbSum168, AbbSum545, AbbSum116, AbbSum112
	common /gzabbrev/ AbbSum148, AbbSum145, AbbSum69, AbbSum386
	common /gzabbrev/ AbbSum310, AbbSum179, AbbSum175, AbbSum8
	common /gzabbrev/ AbbSum316, AbbSum348, AbbSum317, AbbSum453
	common /gzabbrev/ AbbSum585, AbbSum189, AbbSum336, AbbSum431
	common /gzabbrev/ AbbSum558, AbbSum591, AbbSum586, AbbSum61
	common /gzabbrev/ AbbSum60, AbbSum100, AbbSum247, AbbSum208
	common /gzabbrev/ AbbSum120, AbbSum880, AbbSum134, AbbSum118
	common /gzabbrev/ AbbSum874, AbbSum132, AbbSum207, AbbSum99
	common /gzabbrev/ AbbSum246, AbbSum846, AbbSum915, AbbSum845
	common /gzabbrev/ AbbSum66, AbbSum150, AbbSum206, AbbSum97
	common /gzabbrev/ AbbSum151, AbbSum262, AbbSum149, AbbSum104
	common /gzabbrev/ AbbSum245, AbbSum595, AbbSum896, AbbSum751
	common /gzabbrev/ AbbSum244, AbbSum680, AbbSum806, AbbSum230
	common /gzabbrev/ AbbSum186, AbbSum111, AbbSum212, AbbSum137
	common /gzabbrev/ AbbSum174, AbbSum893, AbbSum606, AbbSum254
	common /gzabbrev/ AbbSum844, AbbSum259, AbbSum129, AbbSum960
	common /gzabbrev/ AbbSum713, AbbSum471, AbbSum203, AbbSum82
	common /gzabbrev/ AbbSum383, AbbSum200, AbbSum659, AbbSum385
	common /gzabbrev/ AbbSum446, AbbSum430, AbbSum1036, AbbSum297
	common /gzabbrev/ AbbSum1035, AbbSum172, AbbSum658, AbbSum96
	common /gzabbrev/ AbbSum608, AbbSum226, AbbSum195, AbbSum63
	common /gzabbrev/ AbbSum759, AbbSum588, AbbSum1095, AbbSum565
	common /gzabbrev/ AbbSum243, AbbSum665, AbbSum895, AbbSum242
	common /gzabbrev/ AbbSum218, AbbSum135, AbbSum228, AbbSum106
	common /gzabbrev/ AbbSum201, AbbSum610, AbbSum136, AbbSum225
	common /gzabbrev/ AbbSum173, AbbSum889, AbbSum593, AbbSum1094
	common /gzabbrev/ AbbSum220, AbbSum204, AbbSum261, AbbSum253
	common /gzabbrev/ AbbSum216, AbbSum224, AbbSum167, AbbSum1060
	common /gzabbrev/ AbbSum1061, AbbSum258, AbbSum154, AbbSum843
	common /gzabbrev/ AbbSum1058, AbbSum124, AbbSum941, AbbSum913
	common /gzabbrev/ AbbSum1059, AbbSum103, AbbSum618, AbbSum914
	common /gzabbrev/ AbbSum704, AbbSum470, AbbSum202, AbbSum81
	common /gzabbrev/ AbbSum360, AbbSum197, AbbSum643, AbbSum362
	common /gzabbrev/ AbbSum434, AbbSum413, AbbSum185, AbbSum110
	common /gzabbrev/ AbbSum1032, AbbSum296, AbbSum1030, AbbSum171
	common /gzabbrev/ AbbSum642, AbbSum93, AbbSum596, AbbSum805
	common /gzabbrev/ AbbSum750, AbbSum255, AbbSum193, AbbSum62
	common /gzabbrev/ AbbSum211, AbbSum758, AbbSum577, AbbSum1004
	common /gzabbrev/ AbbSum1033, AbbSum1055, AbbSum1108
	common /gzabbrev/ AbbSum918, AbbSum39, AbbSum40, AbbSum1003
	common /gzabbrev/ AbbSum1031, AbbSum1054, AbbSum141, AbbSum221
	common /gzabbrev/ AbbSum67, AbbSum156, AbbSum363, AbbSum108
	common /gzabbrev/ AbbSum157, AbbSum80, AbbSum184, AbbSum155
	common /gzabbrev/ AbbSum138, AbbSum182, AbbSum414, AbbSum410
	common /gzabbrev/ AbbSum92, AbbSum130, AbbSum433, AbbSum79
	common /gzabbrev/ AbbSum205, AbbSum181, AbbSum209, AbbSum458
	common /gzabbrev/ AbbSum152, AbbSum566, AbbSum140, AbbSum567
	common /gzabbrev/ AbbSum290, AbbSum819, AbbSum900, AbbSum1214
	common /gzabbrev/ AbbSum654, AbbSum526, AbbSum271, AbbSum334
	common /gzabbrev/ AbbSum771, AbbSum115, AbbSum907, AbbSum1007
	common /gzabbrev/ AbbSum178, AbbSum249, AbbSum861, AbbSum101
	common /gzabbrev/ AbbSum916, AbbSum1080, AbbSum1096, AbbSum229
	common /gzabbrev/ AbbSum44, AbbSum1062, AbbSum213, AbbSum1056
	common /gzabbrev/ AbbSum50, AbbSum222, AbbSum142, AbbSum902
	common /gzabbrev/ AbbSum1005, AbbSum1109, AbbSum256
	common /gzabbrev/ AbbSum1064, AbbSum1088, AbbSum65, AbbSum117
	common /gzabbrev/ AbbSum113, AbbSum908, AbbSum449, AbbSum1008
	common /gzabbrev/ AbbSum180, AbbSum176, AbbSum109, AbbSum250
	common /gzabbrev/ AbbSum248, AbbSum429, AbbSum105, AbbSum867
	common /gzabbrev/ AbbSum102, AbbSum98, AbbSum917, AbbSum153
	common /gzabbrev/ AbbSum575, AbbSum616, AbbSum1083, AbbSum235
	common /gzabbrev/ AbbSum605, AbbSum745, AbbSum492, AbbSum71
	common /gzabbrev/ AbbSum927, AbbSum864, AbbSum882, AbbSum1052
	common /gzabbrev/ AbbSum274, AbbSum928, AbbSum1097, AbbSum231
	common /gzabbrev/ AbbSum107, AbbSum48, AbbSum1063, AbbSum217
	common /gzabbrev/ AbbSum214, AbbSum210, AbbSum1057, AbbSum215
	common /gzabbrev/ AbbSum68, AbbSum52, AbbSum128, AbbSum223
	common /gzabbrev/ AbbSum219, AbbSum564, AbbSum656, AbbSum143
	common /gzabbrev/ AbbSum139, AbbSum619, AbbSum572, AbbSum906
	common /gzabbrev/ AbbSum1006, AbbSum183, AbbSum462, AbbSum496
	common /gzabbrev/ AbbSum1110, AbbSum252, AbbSum428, AbbSum556
	common /gzabbrev/ AbbSum574, AbbSum573, AbbSum257, AbbSum251
	common /gzabbrev/ AbbSum1065, AbbSum227, AbbSum1093, AbbSum240
	common /gzabbrev/ AbbSum587, AbbSum849, AbbSum740, AbbSum162
	common /gzabbrev/ AbbSum84, AbbSum967, AbbSum1114, AbbSum88
	common /gzabbrev/ AbbSum1130, AbbSum922, AbbSum1116
	common /gzabbrev/ AbbSum1066, AbbSum761, AbbSum754, AbbSum741
	common /gzabbrev/ AbbSum1140, AbbSum1151, AbbSum972
	common /gzabbrev/ AbbSum1143, AbbSum72, AbbSum1184, AbbSum978
	common /gzabbrev/ AbbSum979, AbbSum294, AbbSum1153, AbbSum1199
	common /gzabbrev/ AbbSum980, AbbSum292, AbbSum919, AbbSum956
	common /gzabbrev/ AbbSum945, AbbSum119, AbbSum755, AbbSum742
	common /gzabbrev/ AbbSum746, AbbSum163, AbbSum87, AbbSum970
	common /gzabbrev/ AbbSum1115, AbbSum90, AbbSum1131, AbbSum923
	common /gzabbrev/ AbbSum1117, AbbSum1074, AbbSum786, AbbSum942
	common /gzabbrev/ AbbSum121, AbbSum161, AbbSum75, AbbSum78
	common /gzabbrev/ AbbSum1044, AbbSum852, AbbSum981, AbbSum1147
	common /gzabbrev/ AbbSum1159, AbbSum973, AbbSum1148, AbbSum467
	common /gzabbrev/ AbbSum73, AbbSum1193, AbbSum920, AbbSum988
	common /gzabbrev/ AbbSum295, AbbSum1160, AbbSum1207, AbbSum989
	common /gzabbrev/ AbbSum885, AbbSum887, AbbSum70, AbbSum122
	common /gzabbrev/ AbbSum853, AbbSum1045, AbbSum89, AbbSum293
	common /gzabbrev/ AbbSum847, AbbSum961, AbbSum83, AbbSum86
	common /gzabbrev/ AbbSum329, AbbSum764, AbbSum996, AbbSum436
	common /gzabbrev/ AbbSum692, AbbSum938, AbbSum507, AbbSum1132
	common /gzabbrev/ AbbSum812, AbbSum752, AbbSum1009, AbbSum597
	common /gzabbrev/ AbbSum626, AbbSum169, AbbSum578, AbbSum1119
	common /gzabbrev/ AbbSum716, AbbSum765, AbbSum1103, AbbSum1010
	common /gzabbrev/ AbbSum832, AbbSum396, AbbSum476, AbbSum1038
	common /gzabbrev/ AbbSum1024, AbbSum872, AbbSum934, AbbSum909
	common /gzabbrev/ AbbSum1142, AbbSum535, AbbSum328, AbbSum666
	common /gzabbrev/ AbbSum474, AbbSum504, AbbSum691, AbbSum506
	common /gzabbrev/ AbbSum718, AbbSum532, AbbSum394, AbbSum549
	common /gzabbrev/ AbbSum831, AbbSum435, AbbSum1186, AbbSum813
	common /gzabbrev/ AbbSum361, AbbSum1098, AbbSum762, AbbSum579
	common /gzabbrev/ AbbSum1011, AbbSum662, AbbSum594, AbbSum623
	common /gzabbrev/ AbbSum641, AbbSum703, AbbSum1015, AbbSum944
	common /gzabbrev/ AbbSum1201, AbbSum1187, AbbSum1154
	common /gzabbrev/ AbbSum1072, AbbSum624, AbbSum969, AbbSum190
	common /gzabbrev/ AbbSum1126, AbbSum1022, AbbSum1073
	common /gzabbrev/ AbbSum1127, AbbSum1016, AbbSum315, AbbSum547
	common /gzabbrev/ AbbSum1111, AbbSum966, AbbSum53, AbbSum350
	common /gzabbrev/ AbbSum788, AbbSum1001, AbbSum448, AbbSum700
	common /gzabbrev/ AbbSum530, AbbSum1135, AbbSum822, AbbSum830
	common /gzabbrev/ AbbSum935, AbbSum1185, AbbSum1200
	common /gzabbrev/ AbbSum1018, AbbSum943, AbbSum609, AbbSum639
	common /gzabbrev/ AbbSum170, AbbSum589, AbbSum1128, AbbSum1068
	common /gzabbrev/ AbbSum731, AbbSum305, AbbSum326, AbbSum717
	common /gzabbrev/ AbbSum534, AbbSum508, AbbSum627, AbbSum123
	common /gzabbrev/ AbbSum661, AbbSum811, AbbSum995, AbbSum997
	common /gzabbrev/ AbbSum814, AbbSum839, AbbSum408, AbbSum500
	common /gzabbrev/ AbbSum1141, AbbSum1152, AbbSum395, AbbSum531
	common /gzabbrev/ AbbSum897, AbbSum412, AbbSum576, AbbSum715
	common /gzabbrev/ AbbSum452, AbbSum898, AbbSum660, AbbSum546
	common /gzabbrev/ AbbSum349, AbbSum910, AbbSum681, AbbSum499
	common /gzabbrev/ AbbSum528, AbbSum699, AbbSum529, AbbSum505
	common /gzabbrev/ AbbSum690, AbbSum392, AbbSum732, AbbSum544
	common /gzabbrev/ AbbSum407, AbbSum555, AbbSum939, AbbSum753
	common /gzabbrev/ AbbSum447, AbbSum1194, AbbSum823, AbbSum384
	common /gzabbrev/ AbbSum1039, AbbSum1102, AbbSum1099
	common /gzabbrev/ AbbSum787, AbbSum590, AbbSum663, AbbSum607
	common /gzabbrev/ AbbSum637, AbbSum657, AbbSum712, AbbSum1120
	common /gzabbrev/ AbbSum763, AbbSum1021, AbbSum1025
	common /gzabbrev/ AbbSum1067, AbbSum873, AbbSum930, AbbSum304
	common /gzabbrev/ AbbSum625, AbbSum702, AbbSum475, AbbSum411
	common /gzabbrev/ AbbSum638, AbbSum325, AbbSum971, AbbSum192
	common /gzabbrev/ AbbSum1129, AbbSum1017, AbbSum1023
	common /gzabbrev/ AbbSum308, AbbSum309, AbbSum188, AbbSum548
	common /gzabbrev/ AbbSum640, AbbSum557, AbbSum359, AbbSum454
	common /gzabbrev/ AbbSum489, AbbSum343, AbbSum1012, AbbSum697
	common /gzabbrev/ AbbSum644, AbbSum485, AbbSum388, AbbSum390
	common /gzabbrev/ AbbSum288, AbbSum1213, AbbSum278, AbbSum1137
	common /gzabbrev/ AbbSum298, AbbSum523, AbbSum419, AbbSum191
	common /gzabbrev/ AbbSum409, AbbSum32, AbbSum1076, AbbSum187
	common /gzabbrev/ AbbSum553, AbbSum352, AbbSum424, AbbSum486
	common /gzabbrev/ AbbSum420, AbbSum456, AbbSum306, AbbSum1107
	common /gzabbrev/ AbbSum311, AbbSum377, AbbSum387, AbbSum332
	common /gzabbrev/ AbbSum866, AbbSum696, AbbSum710, AbbSum381
	common /gzabbrev/ AbbSum291, AbbSum792, AbbSum267, AbbSum289
	common /gzabbrev/ AbbSum714, AbbSum826, AbbSum1029, AbbSum127
	common /gzabbrev/ AbbSum51, AbbSum333, AbbSum990, AbbSum840
	common /gzabbrev/ AbbSum904, AbbSum1113, AbbSum318, AbbSum701
	common /gzabbrev/ AbbSum353, AbbSum469, AbbSum1037, AbbSum416
	common /gzabbrev/ AbbSum307, AbbSum354, AbbSum952, AbbSum391
	common /gzabbrev/ AbbSum37, AbbSum241, AbbSum1215, AbbSum793
	common /gzabbrev/ AbbSum865, AbbSum421, AbbSum342, AbbSum59
	common /gzabbrev/ AbbSum477, AbbSum478, AbbSum272, AbbSum389
	common /gzabbrev/ AbbSum270, AbbSum351, AbbSum35, AbbSum335
	common /gzabbrev/ AbbSum427, AbbSum422, AbbSum460, AbbSum312
	common /gzabbrev/ AbbSum791, AbbSum415, AbbSum423, AbbSum459
	common /gzabbrev/ AbbSum881, AbbSum1043, AbbSum319, AbbSum863
	common /gzabbrev/ AbbSum688, AbbSum33, AbbSum689, AbbSum323
	common /gzabbrev/ AbbSum264, AbbSum324, AbbSum695, AbbSum491
	common /gzabbrev/ AbbSum368, AbbSum404, AbbSum631, AbbSum773
	common /gzabbrev/ AbbSum344, AbbSum648, AbbSum512, AbbSum522
	common /gzabbrev/ AbbSum671, AbbSum398, AbbSum331, AbbSum542
	common /gzabbrev/ AbbSum345, AbbSum543, AbbSum403, AbbSum513
	common /gzabbrev/ AbbSum537, AbbSum632, AbbSum466, AbbSum378
	common /gzabbrev/ AbbSum31, AbbSum650, AbbSum633, AbbSum266
	common /gzabbrev/ AbbSum774, AbbSum775, AbbSum958, AbbSum959
	common /gzabbrev/ AbbSum635, AbbSum698, AbbSum636, AbbSum379
	common /gzabbrev/ AbbSum655, AbbSum279, AbbSum495, AbbSum265
	common /gzabbrev/ AbbSum653, AbbSum924, AbbSum1040, AbbSum417
	common /gzabbrev/ AbbSum855, AbbSum559, AbbSum560, AbbSum570
	common /gzabbrev/ AbbSum856, AbbSum372, AbbSum598, AbbSum337
	common /gzabbrev/ AbbSum479, AbbSum517, AbbSum738, AbbSum330
	common /gzabbrev/ AbbSum488, AbbSum490, AbbSum480, AbbSum487
	common /gzabbrev/ AbbSum501, AbbSum705, AbbSum707, AbbSum369
	common /gzabbrev/ AbbSum514, AbbSum876, AbbSum287, AbbSum1070
	common /gzabbrev/ AbbSum49, AbbSum984, AbbSum857, AbbSum58
	common /gzabbrev/ AbbSum126, AbbSum286, AbbSum442, AbbSum418
	common /gzabbrev/ AbbSum373, AbbSum461, AbbSum951, AbbSum164
	common /gzabbrev/ AbbSum1212, AbbSum706, AbbSum694, AbbSum339
	common /gzabbrev/ AbbSum365, AbbSum735, AbbSum370, AbbSum371
	common /gzabbrev/ AbbSum515, AbbSum338, AbbSum518, AbbSum645
	common /gzabbrev/ AbbSum509, AbbSum516, AbbSum482, AbbSum340
	common /gzabbrev/ AbbSum481, AbbSum708, AbbSum277, AbbSum269
	common /gzabbrev/ AbbSum341, AbbSum709, AbbSum693, AbbSum1048
	common /gzabbrev/ AbbSum833, AbbSum366, AbbSum367, AbbSum646
	common /gzabbrev/ AbbSum510, AbbSum483, AbbSum649, AbbSum571
	common /gzabbrev/ AbbSum613, AbbSum647, AbbSum511, AbbSum374
	common /gzabbrev/ AbbSum561, AbbSum1112, AbbSum562, AbbSum260
	common /gzabbrev/ AbbSum273, AbbSum749, AbbSum94, AbbSum76
	common /gzabbrev/ AbbSum95, AbbSum1134, AbbSum859, AbbSum1026
	common /gzabbrev/ AbbSum1034, AbbSum877, AbbSum380, AbbSum375
	common /gzabbrev/ AbbSum34, AbbSum875, AbbSum784, AbbSum268
	common /gzabbrev/ AbbSum711, AbbSum563, AbbSum652, AbbSum198
	common /gzabbrev/ AbbSum878, AbbSum785, AbbSum196, AbbSum199
	common /gzabbrev/ AbbSum879, AbbSum860, AbbSum165, AbbSum425
	common /gzabbrev/ AbbSum614, AbbSum376, AbbSum275, AbbSum444
	common /gzabbrev/ AbbSum445, AbbSum426, AbbSum125, AbbSum493
	common /gzabbrev/ AbbSum36, AbbSum494, AbbSum525, AbbSum276
	common /gzabbrev/ AbbSum670, AbbSum622, AbbSum651, AbbSum925
	common /gzabbrev/ AbbSum634, AbbSum615, AbbSum1106, AbbSum484
	common /gzabbrev/ AbbSum985, AbbSum54, AbbSum1166, AbbSum1169
	common /gzabbrev/ AbbSum1014, AbbSum1164, AbbSum1189
	common /gzabbrev/ AbbSum899, AbbSum56, AbbSum744, AbbSum1220
	common /gzabbrev/ AbbSum1218, AbbSum1078, AbbSum1202
	common /gzabbrev/ AbbSum1174, AbbSum1047, AbbSum836, AbbSum281
	common /gzabbrev/ AbbSum9, AbbSum19, AbbSum779, AbbSum1
	common /gzabbrev/ AbbSum1216, AbbSum15, AbbSum10, AbbSum17
	common /gzabbrev/ AbbSum767, AbbSum1013, AbbSum43, AbbSum780
	common /gzabbrev/ AbbSum768, AbbSum1177, AbbSum911, AbbSum1210
	common /gzabbrev/ AbbSum835, AbbSum756, AbbSum770, AbbSum800
	common /gzabbrev/ AbbSum1133, AbbSum1190, AbbSum801
	common /gzabbrev/ AbbSum1191, AbbSum778, AbbSum766, AbbSum931
	common /gzabbrev/ AbbSum41, AbbSum1176, AbbSum975, AbbSum883
	common /gzabbrev/ AbbSum42, AbbSum1049, AbbSum1046, AbbSum950
	common /gzabbrev/ AbbSum743, AbbSum999, AbbSum948, AbbSum949
	common /gzabbrev/ AbbSum946, AbbSum947, AbbSum1155, AbbSum1086
	common /gzabbrev/ AbbSum1087, AbbSum854, AbbSum1069, AbbSum817
	common /gzabbrev/ AbbSum890, AbbSum1071, AbbSum848, AbbSum1085
	common /gzabbrev/ AbbSum1192, AbbSum818, AbbSum1203, AbbSum13
	common /gzabbrev/ AbbSum1206, AbbSum729, AbbSum1100, AbbSum21
	common /gzabbrev/ AbbSum1084, AbbSum1146, AbbSum987, AbbSum901
	common /gzabbrev/ AbbSum4, AbbSum1079, AbbSum936, AbbSum1162
	common /gzabbrev/ AbbSum22, AbbSum29, AbbSum23, AbbSum797
	common /gzabbrev/ AbbSum27, AbbSum827, AbbSum664, AbbSum472
	common /gzabbrev/ AbbSum3, AbbSum828, AbbSum356, AbbSum473
	common /gzabbrev/ AbbSum921, AbbSum810, AbbSum321, AbbSum503
	common /gzabbrev/ AbbSum781, AbbSum986, AbbSum322, AbbSum991
	common /gzabbrev/ AbbSum674, AbbSum687, AbbSum55, AbbSum727
	common /gzabbrev/ AbbSum1165, AbbSum1020, AbbSum903
	common /gzabbrev/ AbbSum1167, AbbSum1171, AbbSum580, AbbSum673
	common /gzabbrev/ AbbSum284, AbbSum722, AbbSum283, AbbSum684
	common /gzabbrev/ AbbSum1221, AbbSum57, AbbSum748, AbbSum1219
	common /gzabbrev/ AbbSum1081, AbbSum686, AbbSum685, AbbSum302
	common /gzabbrev/ AbbSum232, AbbSum437, AbbSum438, AbbSum955
	common /gzabbrev/ AbbSum782, AbbSum783, AbbSum798, AbbSum954
	common /gzabbrev/ AbbSum1175, AbbSum1051, AbbSum1188
	common /gzabbrev/ AbbSum629, AbbSum723, AbbSum724, AbbSum539
	common /gzabbrev/ AbbSum683, AbbSum842, AbbSum1136, AbbSum282
	common /gzabbrev/ AbbSum796, AbbSum794, AbbSum976, AbbSum46
	common /gzabbrev/ AbbSum11, AbbSum790, AbbSum1195, AbbSum912
	common /gzabbrev/ AbbSum20, AbbSum795, AbbSum47, AbbSum1211
	common /gzabbrev/ AbbSum888, AbbSum1053, AbbSum2, AbbSum841
	common /gzabbrev/ AbbSum804, AbbSum789, AbbSum1050, AbbSum965
	common /gzabbrev/ AbbSum1217, AbbSum1196, AbbSum933, AbbSum45
	common /gzabbrev/ AbbSum16, AbbSum757, AbbSum1179, AbbSum1019
	common /gzabbrev/ AbbSum12, AbbSum18, AbbSum1122, AbbSum675
	common /gzabbrev/ AbbSum628, AbbSum280, AbbSum402, AbbSum397
	common /gzabbrev/ AbbSum401, AbbSum672, AbbSum405, AbbSum725
	common /gzabbrev/ AbbSum728, AbbSum540, AbbSum301, AbbSum554
	common /gzabbrev/ AbbSum399, AbbSum667, AbbSum721, AbbSum726
	common /gzabbrev/ AbbSum551, AbbSum581, AbbSum737, AbbSum734
	common /gzabbrev/ AbbSum550, AbbSum1121, AbbSum285, AbbSum233
	common /gzabbrev/ AbbSum300, AbbSum630, AbbSum974, AbbSum263
	common /gzabbrev/ AbbSum91, AbbSum599, AbbSum668, AbbSum669
	common /gzabbrev/ AbbSum166, AbbSum739, AbbSum131, AbbSum464
	common /gzabbrev/ AbbSum463, AbbSum519, AbbSum620, AbbSum521
	common /gzabbrev/ AbbSum520, AbbSum465, AbbSum621, AbbSum747
	common /gzabbrev/ AbbSum676, AbbSum1002, AbbSum584, AbbSum541
	common /gzabbrev/ AbbSum538, AbbSum677, AbbSum678, AbbSum600
	common /gzabbrev/ AbbSum964, AbbSum962, AbbSum963, AbbSum1161
	common /gzabbrev/ AbbSum1091, AbbSum582, AbbSum455, AbbSum604
	common /gzabbrev/ AbbSum238, AbbSum439, AbbSum440, AbbSum601
	common /gzabbrev/ AbbSum602, AbbSum603, AbbSum552, AbbSum400
	common /gzabbrev/ AbbSum1092, AbbSum457, AbbSum239, AbbSum432
	common /gzabbrev/ AbbSum568, AbbSum592, AbbSum569, AbbSum441
	common /gzabbrev/ AbbSum862, AbbSum720, AbbSum1075, AbbSum583
	common /gzabbrev/ AbbSum824, AbbSum1077, AbbSum1197, AbbSum850
	common /gzabbrev/ AbbSum825, AbbSum1208, AbbSum1090, AbbSum776
	common /gzabbrev/ AbbSum894, AbbSum799, AbbSum957, AbbSum1178
	common /gzabbrev/ AbbSum772, AbbSum803, AbbSum837, AbbSum769
	common /gzabbrev/ AbbSum802, AbbSum1156, AbbSum617, AbbSum133
	common /gzabbrev/ AbbSum468, AbbSum237, AbbSum1000, AbbSum998
	common /gzabbrev/ AbbSum1157, AbbSum14, AbbSum1101, AbbSum24
	common /gzabbrev/ AbbSum992, AbbSum6, AbbSum937, AbbSum25
	common /gzabbrev/ AbbSum30, AbbSum26, AbbSum28, AbbSum1209
	common /gzabbrev/ AbbSum1082, AbbSum1149, AbbSum1163
	common /gzabbrev/ AbbSum733, AbbSum905, AbbSum1089, AbbSum299
	common /gzabbrev/ AbbSum1125, AbbSum682, AbbSum85, AbbSum1124
	common /gzabbrev/ AbbSum730, AbbSum736, AbbSum443, AbbSum236
	common /gzabbrev/ AbbSum612, AbbSum524, AbbSum194, AbbSum234
	common /gzabbrev/ AbbSum536, AbbSum968, AbbSum1204, AbbSum982
	common /gzabbrev/ AbbSum1144, AbbSum816, AbbSum815, AbbSum1145
	common /gzabbrev/ AbbSum834, AbbSum719, AbbSum611, AbbSum983
	common /gzabbrev/ AbbSum1158, AbbSum891, AbbSum1205, AbbSum884
	common /gzabbrev/ AbbSum1027, AbbSum1170, AbbSum838, AbbSum320
	common /gzabbrev/ AbbSum1139, AbbSum993, AbbSum1105, AbbSum953
	common /gzabbrev/ AbbSum679, AbbSum1138, AbbSum303, AbbSum932
	common /gzabbrev/ AbbSum820, AbbSum497, AbbSum355, AbbSum5
	common /gzabbrev/ AbbSum829, AbbSum994, AbbSum1118, AbbSum74
	common /gzabbrev/ AbbSum1028, AbbSum382, AbbSum450, AbbSum886
	common /gzabbrev/ AbbSum858, AbbSum498, AbbSum357, AbbSum892
	common /gzabbrev/ AbbSum926, AbbSum1104, AbbSum1150, AbbSum807
	common /gzabbrev/ AbbSum1180, AbbSum871, AbbSum358, AbbSum821
	common /gzabbrev/ AbbSum1198, AbbSum870, AbbSum502, AbbSum940
	common /gzabbrev/ AbbSum1041, AbbSum346, AbbSum1182, AbbSum527
	common /gzabbrev/ AbbSum977, AbbSum1168, AbbSum1172, AbbSum809
	common /gzabbrev/ AbbSum929, AbbSum851, AbbSum314, AbbSum347
	common /gzabbrev/ AbbSum38, AbbSum869, AbbSum1173, AbbSum77
	common /gzabbrev/ AbbSum1123, AbbSum451, AbbSum760, AbbSum1181
	common /gzabbrev/ AbbSum868, AbbSum1183, AbbSum1042, AbbSum808

	double complex cint1(3), cint2(3), cint3(3), cint4(3)
	double complex cint5(3), cint6(3), cint7(3), cint8(3)
	double complex cint9(3), cint10(3), cint11(3), cint12(3)
	double complex cint13(2), cint14(2), cint15(2), cint16(2)
	double complex cint17(2), cint18, cint19, cint20, cint21
	double complex cint22, cint23, cint24, cint25, cint26(3)
	double complex cint27(3), cint28(3), cint29(3), cint30(3)
	double complex cint31(3), cint32(3), cint33(3), cint34(3)
	double complex cint35(3), cint36(3), cint37(3), cint38(3)
	double complex cint39(3), cint40(3), cint41(3), cint42(3)
	double complex cint43(3), cint44(2,2), cint45(2,2), cint46(3)
	double complex cint47(3), cint48(3)
	common /gzloopint/ cint1, cint2, cint3, cint4, cint5, cint6
	common /gzloopint/ cint7, cint8, cint9, cint10, cint11, cint12
	common /gzloopint/ cint13, cint14, cint15, cint16, cint17
	common /gzloopint/ cint18, cint19, cint20, cint21, cint22
	common /gzloopint/ cint23, cint24, cint25, cint26, cint27
	common /gzloopint/ cint28, cint29, cint30, cint31, cint32
	common /gzloopint/ cint33, cint34, cint35, cint36, cint37
	common /gzloopint/ cint38, cint39, cint40, cint41, cint42
	common /gzloopint/ cint43, cint44, cint45, cint46, cint47
	common /gzloopint/ cint48

	integer*8 iint1(3), iint2(3), iint3(3), iint4(3), iint5(3)
	integer*8 iint6(3), iint7(3), iint8(3), iint9(3), iint10(3)
	integer*8 iint11(3), iint12(3), iint13(3), iint14(3), iint15(3)
	integer*8 iint16(3), iint17(3), iint18(3), iint19(3), iint20(3)
	integer*8 iint21(3), iint22(3), iint23(3), iint24(3), iint25(2)
	integer*8 iint26(2), iint27(2), iint28(2), iint29(2), iint30(2)
	integer*8 iint31(2), iint32(2), iint33(2), iint34(2), iint35(2)
	integer*8 iint36(2), iint37(3), iint38(3), iint39(3)
	integer*8 iint40(2,2), iint41, iint42, iint43(3), iint44(3)
	integer*8 iint45(3), iint46(3), iint47(3), iint48(3), iint49(3)
	integer*8 iint50(3), iint51(3), iint52(3), iint53(3), iint54(3)
	integer*8 iint55(3), iint56(3), iint57(3), iint58(3), iint59(3)
	integer*8 iint60(3), iint61(3), iint62(3), iint63(3), iint64(3)
	integer*8 iint65(3), iint66(3), iint67(3), iint68(3), iint69(3)
	integer*8 iint70(3), iint71(3), iint72(3), iint73(3), iint74(3)
	integer*8 iint75(3), iint76(3), iint77(3), iint78(3), iint79(3)
	integer*8 iint80(3), iint81(3), iint82(3), iint83(3), iint84(3)
	integer*8 iint85(3), iint86(3), iint87(3), iint88(3), iint89(3)
	integer*8 iint90(3), iint91(3), iint92(3), iint93(3), iint94(3)
	integer*8 iint95(3), iint96(3), iint97(3), iint98(3), iint99(3)
	integer*8 iint100(3), iint101(3), iint102(3), iint103(3)
	integer*8 iint104(3), iint105(3), iint106(3), iint107(3)
	integer*8 iint108(3), iint109(3), iint110(3), iint111(3)
	integer*8 iint112(3), iint113(3), iint114(3), iint115(3)
	integer*8 iint116(3), iint117(3), iint118(3), iint119(3)
	integer*8 iint120(3), iint121(3), iint122(3), iint123(3)
	integer*8 iint124(3), iint125(3), iint126(3), iint127(3)
	integer*8 iint128(3), iint129(3), iint130(3), iint131(3)
	integer*8 iint132(3), iint133(3), iint134(3), iint135(3)
	integer*8 iint136(3), iint137(3), iint138(3), iint139(3)
	integer*8 iint140(3), iint141(3), iint142(3), iint143(3)
	integer*8 iint144(3), iint145(2), iint146(2), iint147(2)
	integer*8 iint148(2), iint149(2), iint150(2), iint151(2)
	integer*8 iint152(2), iint153(2), iint154(2), iint155(2)
	integer*8 iint156(2), iint157(2), iint158(2), iint159(2)
	integer*8 iint160(2,2), iint161(2,2), iint162(2,2), iint163(2,2)
	integer*8 iint164(2), iint165(2), iint166(2,2), iint167(2,2)
	integer*8 iint168(2), iint169(2), iint170(2,2), iint171(2,2)
	integer*8 iint172(2,2), iint173(2,2), iint174(2,2), iint175(2,2)
	integer*8 iint176(2,2), iint177(2,2)
	common /gzloopint/ iint1, iint2, iint3, iint4, iint5, iint6
	common /gzloopint/ iint7, iint8, iint9, iint10, iint11, iint12
	common /gzloopint/ iint13, iint14, iint15, iint16, iint17
	common /gzloopint/ iint18, iint19, iint20, iint21, iint22
	common /gzloopint/ iint23, iint24, iint25, iint26, iint27
	common /gzloopint/ iint28, iint29, iint30, iint31, iint32
	common /gzloopint/ iint33, iint34, iint35, iint36, iint37
	common /gzloopint/ iint38, iint39, iint40, iint41, iint42
	common /gzloopint/ iint43, iint44, iint45, iint46, iint47
	common /gzloopint/ iint48, iint49, iint50, iint51, iint52
	common /gzloopint/ iint53, iint54, iint55, iint56, iint57
	common /gzloopint/ iint58, iint59, iint60, iint61, iint62
	common /gzloopint/ iint63, iint64, iint65, iint66, iint67
	common /gzloopint/ iint68, iint69, iint70, iint71, iint72
	common /gzloopint/ iint73, iint74, iint75, iint76, iint77
	common /gzloopint/ iint78, iint79, iint80, iint81, iint82
	common /gzloopint/ iint83, iint84, iint85, iint86, iint87
	common /gzloopint/ iint88, iint89, iint90, iint91, iint92
	common /gzloopint/ iint93, iint94, iint95, iint96, iint97
	common /gzloopint/ iint98, iint99, iint100, iint101, iint102
	common /gzloopint/ iint103, iint104, iint105, iint106, iint107
	common /gzloopint/ iint108, iint109, iint110, iint111, iint112
	common /gzloopint/ iint113, iint114, iint115, iint116, iint117
	common /gzloopint/ iint118, iint119, iint120, iint121, iint122
	common /gzloopint/ iint123, iint124, iint125, iint126, iint127
	common /gzloopint/ iint128, iint129, iint130, iint131, iint132
	common /gzloopint/ iint133, iint134, iint135, iint136, iint137
	common /gzloopint/ iint138, iint139, iint140, iint141, iint142
	common /gzloopint/ iint143, iint144, iint145, iint146, iint147
	common /gzloopint/ iint148, iint149, iint150, iint151, iint152
	common /gzloopint/ iint153, iint154, iint155, iint156, iint157
	common /gzloopint/ iint158, iint159, iint160, iint161, iint162
	common /gzloopint/ iint163, iint164, iint165, iint166, iint167
	common /gzloopint/ iint168, iint169, iint170, iint171, iint172
	common /gzloopint/ iint173, iint174, iint175, iint176, iint177

	integer cha1, cha2, his1, lpd1, neu1, qud1, quu1, sld1, sle1
	integer sqd1, sqe1, squ1, sqv1
	common /gzindices/ cha1, cha2, his1, lpd1, neu1, qud1, quu1
	common /gzindices/ sld1, sle1, sqd1, sqe1, squ1, sqv1

	double complex Cloop(1)
	common /gzcoeff/ Cloop
