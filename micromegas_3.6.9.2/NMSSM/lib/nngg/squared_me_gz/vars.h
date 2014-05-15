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
	double complex AbbSum165, AbbSum150, AbbSum120, AbbSum153
	double complex AbbSum183, AbbSum7, AbbSum328, AbbSum70
	double complex AbbSum166, AbbSum164, AbbSum152, AbbSum174
	double complex AbbSum122, AbbSum118, AbbSum154, AbbSum151
	double complex AbbSum75, AbbSum185, AbbSum181, AbbSum8
	double complex AbbSum339, AbbSum471, AbbSum195, AbbSum331
	double complex AbbSum473, AbbSum64, AbbSum62, AbbSum106
	double complex AbbSum253, AbbSum214, AbbSum126, AbbSum804
	double complex AbbSum140, AbbSum124, AbbSum798, AbbSum138
	double complex AbbSum213, AbbSum105, AbbSum252, AbbSum774
	double complex AbbSum836, AbbSum773, AbbSum72, AbbSum156
	double complex AbbSum212, AbbSum103, AbbSum157, AbbSum268
	double complex AbbSum155, AbbSum110, AbbSum251, AbbSum425
	double complex AbbSum819, AbbSum692, AbbSum250, AbbSum742
	double complex AbbSum236, AbbSum192, AbbSum117, AbbSum218
	double complex AbbSum143, AbbSum180, AbbSum817, AbbSum1017
	double complex AbbSum260, AbbSum772, AbbSum265, AbbSum135
	double complex AbbSum877, AbbSum431, AbbSum443, AbbSum397
	double complex AbbSum209, AbbSum528, AbbSum88, AbbSum582
	double complex AbbSum206, AbbSum398, AbbSum529, AbbSum347
	double complex AbbSum603, AbbSum399, AbbSum530, AbbSum960
	double complex AbbSum303, AbbSum959, AbbSum178, AbbSum102
	double complex AbbSum1018, AbbSum232, AbbSum201, AbbSum432
	double complex AbbSum69, AbbSum580, AbbSum945, AbbSum1006
	double complex AbbSum1025, AbbSum581, AbbSum433, AbbSum566
	double complex AbbSum249, AbbSum818, AbbSum248, AbbSum224
	double complex AbbSum141, AbbSum234, AbbSum112, AbbSum207
	double complex AbbSum142, AbbSum231, AbbSum179, AbbSum812
	double complex AbbSum1010, AbbSum1024, AbbSum226, AbbSum210
	double complex AbbSum267, AbbSum259, AbbSum222, AbbSum230
	double complex AbbSum173, AbbSum984, AbbSum985, AbbSum264
	double complex AbbSum160, AbbSum771, AbbSum982, AbbSum130
	double complex AbbSum866, AbbSum834, AbbSum983, AbbSum109
	double complex AbbSum835, AbbSum422, AbbSum442, AbbSum384
	double complex AbbSum208, AbbSum522, AbbSum87, AbbSum576
	double complex AbbSum203, AbbSum385, AbbSum523, AbbSum345
	double complex AbbSum602, AbbSum386, AbbSum524, AbbSum191
	double complex AbbSum116, AbbSum954, AbbSum302, AbbSum952
	double complex AbbSum177, AbbSum99, AbbSum1011, AbbSum741
	double complex AbbSum691, AbbSum261, AbbSum199, AbbSum423
	double complex AbbSum68, AbbSum573, AbbSum217, AbbSum942
	double complex AbbSum1002, AbbSum921, AbbSum955, AbbSum979
	double complex AbbSum1040, AbbSum839, AbbSum45, AbbSum46
	double complex AbbSum920, AbbSum953, AbbSum978, AbbSum147
	double complex AbbSum227, AbbSum73, AbbSum162, AbbSum114
	double complex AbbSum163, AbbSum86, AbbSum190, AbbSum161
	double complex AbbSum144, AbbSum188, AbbSum381, AbbSum98
	double complex AbbSum136, AbbSum531, AbbSum85, AbbSum211
	double complex AbbSum187, AbbSum215, AbbSum158, AbbSum146
	double complex AbbSum574, AbbSum424, AbbSum575, AbbSum337
	double complex AbbSum296, AbbSum312, AbbSum313, AbbSum1146
	double complex AbbSum277, AbbSum121, AbbSum828, AbbSum934
	double complex AbbSum184, AbbSum255, AbbSum786, AbbSum107
	double complex AbbSum837, AbbSum1005, AbbSum1028, AbbSum235
	double complex AbbSum50, AbbSum986, AbbSum219, AbbSum980
	double complex AbbSum56, AbbSum228, AbbSum148, AbbSum825
	double complex AbbSum922, AbbSum1041, AbbSum262, AbbSum988
	double complex AbbSum1016, AbbSum71, AbbSum123, AbbSum119
	double complex AbbSum829, AbbSum418, AbbSum941, AbbSum186
	double complex AbbSum182, AbbSum115, AbbSum256, AbbSum254
	double complex AbbSum111, AbbSum791, AbbSum108, AbbSum104
	double complex AbbSum838, AbbSum159, AbbSum537, AbbSum571
	double complex AbbSum1009, AbbSum241, AbbSum784, AbbSum77
	double complex AbbSum805, AbbSum845, AbbSum280, AbbSum1029
	double complex AbbSum237, AbbSum113, AbbSum54, AbbSum561
	double complex AbbSum987, AbbSum223, AbbSum220, AbbSum216
	double complex AbbSum981, AbbSum221, AbbSum74, AbbSum58
	double complex AbbSum134, AbbSum229, AbbSum225, AbbSum149
	double complex AbbSum145, AbbSum827, AbbSum923, AbbSum189
	double complex AbbSum1042, AbbSum258, AbbSum536, AbbSum535
	double complex AbbSum263, AbbSum257, AbbSum989, AbbSum233
	double complex AbbSum1023, AbbSum246, AbbSum547, AbbSum840
	double complex AbbSum683, AbbSum168, AbbSum90, AbbSum882
	double complex AbbSum1046, AbbSum94, AbbSum1070, AbbSum1048
	double complex AbbSum990, AbbSum710, AbbSum684, AbbSum779
	double complex AbbSum991, AbbSum697, AbbSum1080, AbbSum1089
	double complex AbbSum887, AbbSum1083, AbbSum78, AbbSum778
	double complex AbbSum1118, AbbSum896, AbbSum897, AbbSum300
	double complex AbbSum1091, AbbSum992, AbbSum1130, AbbSum898
	double complex AbbSum298, AbbSum867, AbbSum841, AbbSum682
	double complex AbbSum125, AbbSum775, AbbSum698, AbbSum688
	double complex AbbSum169, AbbSum93, AbbSum885, AbbSum1047
	double complex AbbSum96, AbbSum1071, AbbSum1049, AbbSum998
	double complex AbbSum727, AbbSum127, AbbSum167, AbbSum81
	double complex AbbSum84, AbbSum689, AbbSum788, AbbSum999
	double complex AbbSum970, AbbSum899, AbbSum681, AbbSum1086
	double complex AbbSum1096, AbbSum888, AbbSum1087, AbbSum439
	double complex AbbSum79, AbbSum787, AbbSum1126, AbbSum906
	double complex AbbSum301, AbbSum1097, AbbSum1000, AbbSum1138
	double complex AbbSum907, AbbSum808, AbbSum810, AbbSum76
	double complex AbbSum128, AbbSum971, AbbSum95, AbbSum776
	double complex AbbSum299, AbbSum878, AbbSum89, AbbSum92
	double complex AbbSum538, AbbSum703, AbbSum706, AbbSum912
	double complex AbbSum664, AbbSum510, AbbSum857, AbbSum1072
	double complex AbbSum747, AbbSum404, AbbSum693, AbbSum924
	double complex AbbSum590, AbbSum175, AbbSum552, AbbSum1053
	double complex AbbSum1035, AbbSum762, AbbSum370, AbbSum926
	double complex AbbSum539, AbbSum662, AbbSum962, AbbSum946
	double complex AbbSum796, AbbSum502, AbbSum851, AbbSum830
	double complex AbbSum371, AbbSum511, AbbSum403, AbbSum540
	double complex AbbSum554, AbbSum1082, AbbSum763, AbbSum1132
	double complex AbbSum503, AbbSum913, AbbSum512, AbbSum761
	double complex AbbSum405, AbbSum1120, AbbSum748, AbbSum1030
	double complex AbbSum542, AbbSum470, AbbSum556, AbbSum553
	double complex AbbSum1054, AbbSum589, AbbSum640, AbbSum444
	double complex AbbSum604, AbbSum651, AbbSum620, AbbSum482
	double complex AbbSum933, AbbSum928, AbbSum859, AbbSum1133
	double complex AbbSum1121, AbbSum1092, AbbSum329, AbbSum355
	double complex AbbSum997, AbbSum884, AbbSum513, AbbSum406
	double complex AbbSum196, AbbSum1063, AbbSum1064, AbbSum943
	double complex AbbSum1114, AbbSum711, AbbSum587, AbbSum993
	double complex AbbSum925, AbbSum705, AbbSum548, AbbSum1108
	double complex AbbSum704, AbbSum1043, AbbSum881, AbbSum59
	double complex AbbSum918, AbbSum671, AbbSum517, AbbSum1075
	double complex AbbSum754, AbbSum415, AbbSum760, AbbSum852
	double complex AbbSum1119, AbbSum1131, AbbSum860, AbbSum935
	double complex AbbSum600, AbbSum176, AbbSum562, AbbSum1066
	double complex AbbSum129, AbbSum746, AbbSum911, AbbSum914
	double complex AbbSum749, AbbSum419, AbbSum768, AbbSum378
	double complex AbbSum936, AbbSum549, AbbSum1081, AbbSum1090
	double complex AbbSum670, AbbSum501, AbbSum820, AbbSum663
	double complex AbbSum507, AbbSum821, AbbSum402, AbbSum379
	double complex AbbSum518, AbbSum414, AbbSum550, AbbSum564
	double complex AbbSum320, AbbSum769, AbbSum831, AbbSum369
	double complex AbbSum927, AbbSum1065, AbbSum1139, AbbSum508
	double complex AbbSum919, AbbSum519, AbbSum858, AbbSum694
	double complex AbbSum416, AbbSum1127, AbbSum755, AbbSum963
	double complex AbbSum1034, AbbSum1031, AbbSum551, AbbSum472
	double complex AbbSum565, AbbSum563, AbbSum1067, AbbSum599
	double complex AbbSum648, AbbSum446, AbbSum606, AbbSum652
	double complex AbbSum628, AbbSum500, AbbSum940, AbbSum937
	double complex AbbSum947, AbbSum797, AbbSum847, AbbSum641
	double complex AbbSum541, AbbSum654, AbbSum555, AbbSum621
	double complex AbbSum483, AbbSum340, AbbSum366, AbbSum886
	double complex AbbSum520, AbbSum417, AbbSum198, AbbSum1068
	double complex AbbSum1069, AbbSum944, AbbSum728, AbbSum309
	double complex AbbSum317, AbbSum194, AbbSum445, AbbSum605
	double complex AbbSum382, AbbSum521, AbbSum1055, AbbSum322
	double complex AbbSum346, AbbSum598, AbbSum588, AbbSum509
	double complex AbbSum459, AbbSum294, AbbSum284, AbbSum368
	double complex AbbSum1077, AbbSum1145, AbbSum304, AbbSum392
	double complex AbbSum197, AbbSum650, AbbSum380, AbbSum601
	double complex AbbSum38, AbbSum193, AbbSum367, AbbSum341
	double complex AbbSum396, AbbSum460, AbbSum336, AbbSum362
	double complex AbbSum363, AbbSum344, AbbSum318, AbbSum412
	double complex AbbSum467, AbbSum394, AbbSum319, AbbSum393
	double complex AbbSum342, AbbSum311, AbbSum969, AbbSum297
	double complex AbbSum273, AbbSum295, AbbSum315, AbbSum951
	double complex AbbSum133, AbbSum413, AbbSum57, AbbSum649
	double complex AbbSum908, AbbSum395, AbbSum527, AbbSum770
	double complex AbbSum316, AbbSum1045, AbbSum314, AbbSum468
	double complex AbbSum700, AbbSum338, AbbSum961, AbbSum695
	double complex AbbSum469, AbbSum43, AbbSum661, AbbSum247
	double complex AbbSum799, AbbSum1147, AbbSum731, AbbSum462
	double complex AbbSum67, AbbSum365, AbbSum377, AbbSum278
	double complex AbbSum276, AbbSum41, AbbSum696, AbbSum343
	double complex AbbSum701, AbbSum364, AbbSum321, AbbSum478
	double complex AbbSum448, AbbSum350, AbbSum351, AbbSum39
	double complex AbbSum353, AbbSum354, AbbSum449, AbbSum479
	double complex AbbSum480, AbbSum270, AbbSum481, AbbSum646
	double complex AbbSum335, AbbSum647, AbbSum458, AbbSum463
	double complex AbbSum374, AbbSum594, AbbSum358, AbbSum719
	double complex AbbSum623, AbbSum489, AbbSum613, AbbSum373
	double complex AbbSum614, AbbSum595, AbbSum490, AbbSum496
	double complex AbbSum37, AbbSum272, AbbSum596, AbbSum285
	double complex AbbSum271, AbbSum615, AbbSum495, AbbSum506
	double complex AbbSum617, AbbSum964, AbbSum577, AbbSum426
	double complex AbbSum622, AbbSum488, AbbSum658, AbbSum452
	double complex AbbSum842, AbbSum642, AbbSum533, AbbSum430
	double complex AbbSum611, AbbSum453, AbbSum387, AbbSum608
	double complex AbbSum487, AbbSum843, AbbSum525, AbbSum526
	double complex AbbSum612, AbbSum578, AbbSum486, AbbSum609
	double complex AbbSum655, AbbSum388, AbbSum656, AbbSum800
	double complex AbbSum293, AbbSum55, AbbSum901, AbbSum66
	double complex AbbSum132, AbbSum292, AbbSum409, AbbSum390
	double complex AbbSum428, AbbSum170, AbbSum332, AbbSum643
	double complex AbbSum356, AbbSum357, AbbSum454, AbbSum484
	double complex AbbSum657, AbbSum333, AbbSum644, AbbSum610
	double complex AbbSum450, AbbSum283, AbbSum275, AbbSum359
	double complex AbbSum330, AbbSum678, AbbSum451, AbbSum376
	double complex AbbSum455, AbbSum389, AbbSum427, AbbSum844
	double complex AbbSum764, AbbSum456, AbbSum579, AbbSum410
	double complex AbbSum391, AbbSum429, AbbSum645, AbbSum673
	double complex AbbSum607, AbbSum1050, AbbSum1044, AbbSum360
	double complex AbbSum266, AbbSum361, AbbSum279, AbbSum854
	double complex AbbSum100, AbbSum82, AbbSum101, AbbSum855
	double complex AbbSum1074, AbbSum967, AbbSum948, AbbSum956
	double complex AbbSum801, AbbSum464, AbbSum494, AbbSum957
	double complex AbbSum40, AbbSum274, AbbSum1051, AbbSum616
	double complex AbbSum1039, AbbSum968, AbbSum1026, AbbSum204
	double complex AbbSum202, AbbSum958, AbbSum205, AbbSum803
	double complex AbbSum1027, AbbSum171, AbbSum493, AbbSum619
	double complex AbbSum281, AbbSum802, AbbSum131, AbbSum1144
	double complex AbbSum42, AbbSum465, AbbSum282, AbbSum437
	double complex AbbSum492, AbbSum675, AbbSum585, AbbSum659
	double complex AbbSum618, AbbSum411, AbbSum334, AbbSum457
	double complex AbbSum1038, AbbSum558, AbbSum569, AbbSum461
	double complex AbbSum60, AbbSum723, AbbSum1109, AbbSum1102
	double complex AbbSum1104, AbbSum872, AbbSum713, AbbSum1100
	double complex AbbSum822, AbbSum1150, AbbSum63, AbbSum1152
	double complex AbbSum1003, AbbSum1111, AbbSum932, AbbSum1106
	double complex AbbSum973, AbbSum287, AbbSum9, AbbSum15
	double complex AbbSum23, AbbSum724, AbbSum1, AbbSum1148
	double complex AbbSum16, AbbSum715, AbbSum929, AbbSum10
	double complex AbbSum21, AbbSum17, AbbSum714, AbbSum931
	double complex AbbSum49, AbbSum725, AbbSum716, AbbSum1113
	double complex AbbSum1142, AbbSum686, AbbSum699, AbbSum717
	double complex AbbSum832, AbbSum737, AbbSum1073, AbbSum1123
	double complex AbbSum738, AbbSum1124, AbbSum848, AbbSum47
	double complex AbbSum890, AbbSum806, AbbSum48, AbbSum974
	double complex AbbSum972, AbbSum685, AbbSum1058, AbbSum869
	double complex AbbSum870, AbbSum868, AbbSum1014, AbbSum1015
	double complex AbbSum782, AbbSum995, AbbSum994, AbbSum781
	double complex AbbSum1013, AbbSum780, AbbSum813, AbbSum752
	double complex AbbSum1136, AbbSum13, AbbSum1137, AbbSum669
	double complex AbbSum1032, AbbSum25, AbbSum33, AbbSum1012
	double complex AbbSum904, AbbSum823, AbbSum4, AbbSum1004
	double complex AbbSum853, AbbSum26, AbbSum27, AbbSum35
	double complex AbbSum735, AbbSum31, AbbSum1098, AbbSum447
	double complex AbbSum707, AbbSum634, AbbSum586, AbbSum323
	double complex AbbSum862, AbbSum629, AbbSum3, AbbSum758
	double complex AbbSum653, AbbSum894, AbbSum635, AbbSum630
	double complex AbbSum474, AbbSum475, AbbSum903, AbbSum477
	double complex AbbSum349, AbbSum505, AbbSum930, AbbSum667
	double complex AbbSum680, AbbSum636, AbbSum61, AbbSum544
	double complex AbbSum546, AbbSum560, AbbSum639, AbbSum732
	double complex AbbSum1101, AbbSum1116, AbbSum1103, AbbSum1105
	double complex AbbSum687, AbbSum720, AbbSum873, AbbSum975
	double complex AbbSum824, AbbSum290, AbbSum289, AbbSum1151
	double complex AbbSum631, AbbSum1153, AbbSum65, AbbSum1007
	double complex AbbSum1059, AbbSum1122, AbbSum900, AbbSum626
	double complex AbbSum308, AbbSum238, AbbSum1117, AbbSum627
	double complex AbbSum1107, AbbSum939, AbbSum666, AbbSum625
	double complex AbbSum1076, AbbSum288, AbbSum734, AbbSum891
	double complex AbbSum52, AbbSum11, AbbSum730, AbbSum1128
	double complex AbbSum18, AbbSum24, AbbSum733, AbbSum53
	double complex AbbSum1143, AbbSum811, AbbSum977, AbbSum2
	double complex AbbSum690, AbbSum740, AbbSum976, AbbSum1149
	double complex AbbSum1129, AbbSum850, AbbSum51, AbbSum19
	double complex AbbSum702, AbbSum729, AbbSum938, AbbSum12
	double complex AbbSum22, AbbSum20, AbbSum750, AbbSum1084
	double complex AbbSum766, AbbSum1135, AbbSum1057, AbbSum286
	double complex AbbSum514, AbbSum372, AbbSum375, AbbSum591
	double complex AbbSum624, AbbSum307, AbbSum441, AbbSum421
	double complex AbbSum1056, AbbSum291, AbbSum239, AbbSum306
	double complex AbbSum889, AbbSum269, AbbSum97, AbbSum491
	double complex AbbSum172, AbbSum668, AbbSum137, AbbSum434
	double complex AbbSum677, AbbSum583, AbbSum570, AbbSum584
	double complex AbbSum435, AbbSum436, AbbSum568, AbbSum765
	double complex AbbSum1134, AbbSum665, AbbSum504, AbbSum915
	double complex AbbSum676, AbbSum567, AbbSum407, AbbSum1093
	double complex AbbSum880, AbbSum879, AbbSum1021, AbbSum1033
	double complex AbbSum916, AbbSum1094, AbbSum559, AbbSum420
	double complex AbbSum679, AbbSum244, AbbSum592, AbbSum1022
	double complex AbbSum516, AbbSum440, AbbSum245, AbbSum593
	double complex AbbSum545, AbbSum408, AbbSum515, AbbSum1001
	double complex AbbSum790, AbbSum485, AbbSum871, AbbSum833
	double complex AbbSum783, AbbSum712, AbbSum543, AbbSum756
	double complex AbbSum996, AbbSum1140, AbbSum1020, AbbSum721
	double complex AbbSum789, AbbSum736, AbbSum876, AbbSum1115
	double complex AbbSum718, AbbSum739, AbbSum816, AbbSum139
	double complex AbbSum532, AbbSum243, AbbSum1112, AbbSum14
	double complex AbbSum28, AbbSum34, AbbSum909, AbbSum6
	double complex AbbSum856, AbbSum29, AbbSum30, AbbSum36
	double complex AbbSum32, AbbSum1141, AbbSum1008, AbbSum672
	double complex AbbSum826, AbbSum1019, AbbSum1110, AbbSum305
	double complex AbbSum1062, AbbSum91, AbbSum1061, AbbSum572
	double complex AbbSum674, AbbSum534, AbbSum557, AbbSum242
	double complex AbbSum200, AbbSum240, AbbSum883, AbbSum1085
	double complex AbbSum751, AbbSum902, AbbSum1125, AbbSum722
	double complex AbbSum1099, AbbSum874, AbbSum814, AbbSum1095
	double complex AbbSum807, AbbSum917, AbbSum466, AbbSum709
	double complex AbbSum383, AbbSum949, AbbSum767, AbbSum637
	double complex AbbSum1037, AbbSum597, AbbSum352, AbbSum327
	double complex AbbSum324, AbbSum849, AbbSum753, AbbSum632
	double complex AbbSum5, AbbSum759, AbbSum1052, AbbSum80
	double complex AbbSum950, AbbSum875, AbbSum400, AbbSum660
	double complex AbbSum401, AbbSum809, AbbSum785, AbbSum905
	double complex AbbSum438, AbbSum815, AbbSum726, AbbSum638
	double complex AbbSum1036, AbbSum1088, AbbSum795, AbbSum633
	double complex AbbSum910, AbbSum794, AbbSum348, AbbSum965
	double complex AbbSum497, AbbSum895, AbbSum743, AbbSum1078
	double complex AbbSum757, AbbSum892, AbbSum745, AbbSum846
	double complex AbbSum863, AbbSum777, AbbSum310, AbbSum498
	double complex AbbSum476, AbbSum499, AbbSum44, AbbSum864
	double complex AbbSum793, AbbSum893, AbbSum83, AbbSum1060
	double complex AbbSum325, AbbSum792, AbbSum1079, AbbSum326
	double complex AbbSum861, AbbSum865, AbbSum708, AbbSum744
	double complex AbbSum966
	common /gzabbrev/ F9, F11, F1, F3, F10, F5, F7, F12, F4, F15
	common /gzabbrev/ F13, F8, F16, F14, F2, F6, Pair4, Pair5
	common /gzabbrev/ Pair3, Pair1, Pair2, Eps1, Eps2, Abb29
	common /gzabbrev/ Abb30, Abb21, Abb22, Abb35, Abb36, Abb23
	common /gzabbrev/ Abb24, Abb1, Abb41, Abb42, Abb4, Abb2, Abb43
	common /gzabbrev/ Abb44, Abb5, Abb13, Abb16, Abb3, Abb6, Abb31
	common /gzabbrev/ Abb17, Abb37, Abb14, Abb32, Abb18, Abb38
	common /gzabbrev/ Abb15, Abb8, Abb10, Abb11, Abb7, Abb19
	common /gzabbrev/ Abb20, Abb9, Abb12, Abb33, Abb25, Abb39
	common /gzabbrev/ Abb26, Abb34, Abb27, Abb40, Abb28, AbbSum165
	common /gzabbrev/ AbbSum150, AbbSum120, AbbSum153, AbbSum183
	common /gzabbrev/ AbbSum7, AbbSum328, AbbSum70, AbbSum166
	common /gzabbrev/ AbbSum164, AbbSum152, AbbSum174, AbbSum122
	common /gzabbrev/ AbbSum118, AbbSum154, AbbSum151, AbbSum75
	common /gzabbrev/ AbbSum185, AbbSum181, AbbSum8, AbbSum339
	common /gzabbrev/ AbbSum471, AbbSum195, AbbSum331, AbbSum473
	common /gzabbrev/ AbbSum64, AbbSum62, AbbSum106, AbbSum253
	common /gzabbrev/ AbbSum214, AbbSum126, AbbSum804, AbbSum140
	common /gzabbrev/ AbbSum124, AbbSum798, AbbSum138, AbbSum213
	common /gzabbrev/ AbbSum105, AbbSum252, AbbSum774, AbbSum836
	common /gzabbrev/ AbbSum773, AbbSum72, AbbSum156, AbbSum212
	common /gzabbrev/ AbbSum103, AbbSum157, AbbSum268, AbbSum155
	common /gzabbrev/ AbbSum110, AbbSum251, AbbSum425, AbbSum819
	common /gzabbrev/ AbbSum692, AbbSum250, AbbSum742, AbbSum236
	common /gzabbrev/ AbbSum192, AbbSum117, AbbSum218, AbbSum143
	common /gzabbrev/ AbbSum180, AbbSum817, AbbSum1017, AbbSum260
	common /gzabbrev/ AbbSum772, AbbSum265, AbbSum135, AbbSum877
	common /gzabbrev/ AbbSum431, AbbSum443, AbbSum397, AbbSum209
	common /gzabbrev/ AbbSum528, AbbSum88, AbbSum582, AbbSum206
	common /gzabbrev/ AbbSum398, AbbSum529, AbbSum347, AbbSum603
	common /gzabbrev/ AbbSum399, AbbSum530, AbbSum960, AbbSum303
	common /gzabbrev/ AbbSum959, AbbSum178, AbbSum102, AbbSum1018
	common /gzabbrev/ AbbSum232, AbbSum201, AbbSum432, AbbSum69
	common /gzabbrev/ AbbSum580, AbbSum945, AbbSum1006, AbbSum1025
	common /gzabbrev/ AbbSum581, AbbSum433, AbbSum566, AbbSum249
	common /gzabbrev/ AbbSum818, AbbSum248, AbbSum224, AbbSum141
	common /gzabbrev/ AbbSum234, AbbSum112, AbbSum207, AbbSum142
	common /gzabbrev/ AbbSum231, AbbSum179, AbbSum812, AbbSum1010
	common /gzabbrev/ AbbSum1024, AbbSum226, AbbSum210, AbbSum267
	common /gzabbrev/ AbbSum259, AbbSum222, AbbSum230, AbbSum173
	common /gzabbrev/ AbbSum984, AbbSum985, AbbSum264, AbbSum160
	common /gzabbrev/ AbbSum771, AbbSum982, AbbSum130, AbbSum866
	common /gzabbrev/ AbbSum834, AbbSum983, AbbSum109, AbbSum835
	common /gzabbrev/ AbbSum422, AbbSum442, AbbSum384, AbbSum208
	common /gzabbrev/ AbbSum522, AbbSum87, AbbSum576, AbbSum203
	common /gzabbrev/ AbbSum385, AbbSum523, AbbSum345, AbbSum602
	common /gzabbrev/ AbbSum386, AbbSum524, AbbSum191, AbbSum116
	common /gzabbrev/ AbbSum954, AbbSum302, AbbSum952, AbbSum177
	common /gzabbrev/ AbbSum99, AbbSum1011, AbbSum741, AbbSum691
	common /gzabbrev/ AbbSum261, AbbSum199, AbbSum423, AbbSum68
	common /gzabbrev/ AbbSum573, AbbSum217, AbbSum942, AbbSum1002
	common /gzabbrev/ AbbSum921, AbbSum955, AbbSum979, AbbSum1040
	common /gzabbrev/ AbbSum839, AbbSum45, AbbSum46, AbbSum920
	common /gzabbrev/ AbbSum953, AbbSum978, AbbSum147, AbbSum227
	common /gzabbrev/ AbbSum73, AbbSum162, AbbSum114, AbbSum163
	common /gzabbrev/ AbbSum86, AbbSum190, AbbSum161, AbbSum144
	common /gzabbrev/ AbbSum188, AbbSum381, AbbSum98, AbbSum136
	common /gzabbrev/ AbbSum531, AbbSum85, AbbSum211, AbbSum187
	common /gzabbrev/ AbbSum215, AbbSum158, AbbSum146, AbbSum574
	common /gzabbrev/ AbbSum424, AbbSum575, AbbSum337, AbbSum296
	common /gzabbrev/ AbbSum312, AbbSum313, AbbSum1146, AbbSum277
	common /gzabbrev/ AbbSum121, AbbSum828, AbbSum934, AbbSum184
	common /gzabbrev/ AbbSum255, AbbSum786, AbbSum107, AbbSum837
	common /gzabbrev/ AbbSum1005, AbbSum1028, AbbSum235, AbbSum50
	common /gzabbrev/ AbbSum986, AbbSum219, AbbSum980, AbbSum56
	common /gzabbrev/ AbbSum228, AbbSum148, AbbSum825, AbbSum922
	common /gzabbrev/ AbbSum1041, AbbSum262, AbbSum988, AbbSum1016
	common /gzabbrev/ AbbSum71, AbbSum123, AbbSum119, AbbSum829
	common /gzabbrev/ AbbSum418, AbbSum941, AbbSum186, AbbSum182
	common /gzabbrev/ AbbSum115, AbbSum256, AbbSum254, AbbSum111
	common /gzabbrev/ AbbSum791, AbbSum108, AbbSum104, AbbSum838
	common /gzabbrev/ AbbSum159, AbbSum537, AbbSum571, AbbSum1009
	common /gzabbrev/ AbbSum241, AbbSum784, AbbSum77, AbbSum805
	common /gzabbrev/ AbbSum845, AbbSum280, AbbSum1029, AbbSum237
	common /gzabbrev/ AbbSum113, AbbSum54, AbbSum561, AbbSum987
	common /gzabbrev/ AbbSum223, AbbSum220, AbbSum216, AbbSum981
	common /gzabbrev/ AbbSum221, AbbSum74, AbbSum58, AbbSum134
	common /gzabbrev/ AbbSum229, AbbSum225, AbbSum149, AbbSum145
	common /gzabbrev/ AbbSum827, AbbSum923, AbbSum189, AbbSum1042
	common /gzabbrev/ AbbSum258, AbbSum536, AbbSum535, AbbSum263
	common /gzabbrev/ AbbSum257, AbbSum989, AbbSum233, AbbSum1023
	common /gzabbrev/ AbbSum246, AbbSum547, AbbSum840, AbbSum683
	common /gzabbrev/ AbbSum168, AbbSum90, AbbSum882, AbbSum1046
	common /gzabbrev/ AbbSum94, AbbSum1070, AbbSum1048, AbbSum990
	common /gzabbrev/ AbbSum710, AbbSum684, AbbSum779, AbbSum991
	common /gzabbrev/ AbbSum697, AbbSum1080, AbbSum1089, AbbSum887
	common /gzabbrev/ AbbSum1083, AbbSum78, AbbSum778, AbbSum1118
	common /gzabbrev/ AbbSum896, AbbSum897, AbbSum300, AbbSum1091
	common /gzabbrev/ AbbSum992, AbbSum1130, AbbSum898, AbbSum298
	common /gzabbrev/ AbbSum867, AbbSum841, AbbSum682, AbbSum125
	common /gzabbrev/ AbbSum775, AbbSum698, AbbSum688, AbbSum169
	common /gzabbrev/ AbbSum93, AbbSum885, AbbSum1047, AbbSum96
	common /gzabbrev/ AbbSum1071, AbbSum1049, AbbSum998, AbbSum727
	common /gzabbrev/ AbbSum127, AbbSum167, AbbSum81, AbbSum84
	common /gzabbrev/ AbbSum689, AbbSum788, AbbSum999, AbbSum970
	common /gzabbrev/ AbbSum899, AbbSum681, AbbSum1086, AbbSum1096
	common /gzabbrev/ AbbSum888, AbbSum1087, AbbSum439, AbbSum79
	common /gzabbrev/ AbbSum787, AbbSum1126, AbbSum906, AbbSum301
	common /gzabbrev/ AbbSum1097, AbbSum1000, AbbSum1138
	common /gzabbrev/ AbbSum907, AbbSum808, AbbSum810, AbbSum76
	common /gzabbrev/ AbbSum128, AbbSum971, AbbSum95, AbbSum776
	common /gzabbrev/ AbbSum299, AbbSum878, AbbSum89, AbbSum92
	common /gzabbrev/ AbbSum538, AbbSum703, AbbSum706, AbbSum912
	common /gzabbrev/ AbbSum664, AbbSum510, AbbSum857, AbbSum1072
	common /gzabbrev/ AbbSum747, AbbSum404, AbbSum693, AbbSum924
	common /gzabbrev/ AbbSum590, AbbSum175, AbbSum552, AbbSum1053
	common /gzabbrev/ AbbSum1035, AbbSum762, AbbSum370, AbbSum926
	common /gzabbrev/ AbbSum539, AbbSum662, AbbSum962, AbbSum946
	common /gzabbrev/ AbbSum796, AbbSum502, AbbSum851, AbbSum830
	common /gzabbrev/ AbbSum371, AbbSum511, AbbSum403, AbbSum540
	common /gzabbrev/ AbbSum554, AbbSum1082, AbbSum763, AbbSum1132
	common /gzabbrev/ AbbSum503, AbbSum913, AbbSum512, AbbSum761
	common /gzabbrev/ AbbSum405, AbbSum1120, AbbSum748, AbbSum1030
	common /gzabbrev/ AbbSum542, AbbSum470, AbbSum556, AbbSum553
	common /gzabbrev/ AbbSum1054, AbbSum589, AbbSum640, AbbSum444
	common /gzabbrev/ AbbSum604, AbbSum651, AbbSum620, AbbSum482
	common /gzabbrev/ AbbSum933, AbbSum928, AbbSum859, AbbSum1133
	common /gzabbrev/ AbbSum1121, AbbSum1092, AbbSum329, AbbSum355
	common /gzabbrev/ AbbSum997, AbbSum884, AbbSum513, AbbSum406
	common /gzabbrev/ AbbSum196, AbbSum1063, AbbSum1064, AbbSum943
	common /gzabbrev/ AbbSum1114, AbbSum711, AbbSum587, AbbSum993
	common /gzabbrev/ AbbSum925, AbbSum705, AbbSum548, AbbSum1108
	common /gzabbrev/ AbbSum704, AbbSum1043, AbbSum881, AbbSum59
	common /gzabbrev/ AbbSum918, AbbSum671, AbbSum517, AbbSum1075
	common /gzabbrev/ AbbSum754, AbbSum415, AbbSum760, AbbSum852
	common /gzabbrev/ AbbSum1119, AbbSum1131, AbbSum860, AbbSum935
	common /gzabbrev/ AbbSum600, AbbSum176, AbbSum562, AbbSum1066
	common /gzabbrev/ AbbSum129, AbbSum746, AbbSum911, AbbSum914
	common /gzabbrev/ AbbSum749, AbbSum419, AbbSum768, AbbSum378
	common /gzabbrev/ AbbSum936, AbbSum549, AbbSum1081, AbbSum1090
	common /gzabbrev/ AbbSum670, AbbSum501, AbbSum820, AbbSum663
	common /gzabbrev/ AbbSum507, AbbSum821, AbbSum402, AbbSum379
	common /gzabbrev/ AbbSum518, AbbSum414, AbbSum550, AbbSum564
	common /gzabbrev/ AbbSum320, AbbSum769, AbbSum831, AbbSum369
	common /gzabbrev/ AbbSum927, AbbSum1065, AbbSum1139, AbbSum508
	common /gzabbrev/ AbbSum919, AbbSum519, AbbSum858, AbbSum694
	common /gzabbrev/ AbbSum416, AbbSum1127, AbbSum755, AbbSum963
	common /gzabbrev/ AbbSum1034, AbbSum1031, AbbSum551, AbbSum472
	common /gzabbrev/ AbbSum565, AbbSum563, AbbSum1067, AbbSum599
	common /gzabbrev/ AbbSum648, AbbSum446, AbbSum606, AbbSum652
	common /gzabbrev/ AbbSum628, AbbSum500, AbbSum940, AbbSum937
	common /gzabbrev/ AbbSum947, AbbSum797, AbbSum847, AbbSum641
	common /gzabbrev/ AbbSum541, AbbSum654, AbbSum555, AbbSum621
	common /gzabbrev/ AbbSum483, AbbSum340, AbbSum366, AbbSum886
	common /gzabbrev/ AbbSum520, AbbSum417, AbbSum198, AbbSum1068
	common /gzabbrev/ AbbSum1069, AbbSum944, AbbSum728, AbbSum309
	common /gzabbrev/ AbbSum317, AbbSum194, AbbSum445, AbbSum605
	common /gzabbrev/ AbbSum382, AbbSum521, AbbSum1055, AbbSum322
	common /gzabbrev/ AbbSum346, AbbSum598, AbbSum588, AbbSum509
	common /gzabbrev/ AbbSum459, AbbSum294, AbbSum284, AbbSum368
	common /gzabbrev/ AbbSum1077, AbbSum1145, AbbSum304, AbbSum392
	common /gzabbrev/ AbbSum197, AbbSum650, AbbSum380, AbbSum601
	common /gzabbrev/ AbbSum38, AbbSum193, AbbSum367, AbbSum341
	common /gzabbrev/ AbbSum396, AbbSum460, AbbSum336, AbbSum362
	common /gzabbrev/ AbbSum363, AbbSum344, AbbSum318, AbbSum412
	common /gzabbrev/ AbbSum467, AbbSum394, AbbSum319, AbbSum393
	common /gzabbrev/ AbbSum342, AbbSum311, AbbSum969, AbbSum297
	common /gzabbrev/ AbbSum273, AbbSum295, AbbSum315, AbbSum951
	common /gzabbrev/ AbbSum133, AbbSum413, AbbSum57, AbbSum649
	common /gzabbrev/ AbbSum908, AbbSum395, AbbSum527, AbbSum770
	common /gzabbrev/ AbbSum316, AbbSum1045, AbbSum314, AbbSum468
	common /gzabbrev/ AbbSum700, AbbSum338, AbbSum961, AbbSum695
	common /gzabbrev/ AbbSum469, AbbSum43, AbbSum661, AbbSum247
	common /gzabbrev/ AbbSum799, AbbSum1147, AbbSum731, AbbSum462
	common /gzabbrev/ AbbSum67, AbbSum365, AbbSum377, AbbSum278
	common /gzabbrev/ AbbSum276, AbbSum41, AbbSum696, AbbSum343
	common /gzabbrev/ AbbSum701, AbbSum364, AbbSum321, AbbSum478
	common /gzabbrev/ AbbSum448, AbbSum350, AbbSum351, AbbSum39
	common /gzabbrev/ AbbSum353, AbbSum354, AbbSum449, AbbSum479
	common /gzabbrev/ AbbSum480, AbbSum270, AbbSum481, AbbSum646
	common /gzabbrev/ AbbSum335, AbbSum647, AbbSum458, AbbSum463
	common /gzabbrev/ AbbSum374, AbbSum594, AbbSum358, AbbSum719
	common /gzabbrev/ AbbSum623, AbbSum489, AbbSum613, AbbSum373
	common /gzabbrev/ AbbSum614, AbbSum595, AbbSum490, AbbSum496
	common /gzabbrev/ AbbSum37, AbbSum272, AbbSum596, AbbSum285
	common /gzabbrev/ AbbSum271, AbbSum615, AbbSum495, AbbSum506
	common /gzabbrev/ AbbSum617, AbbSum964, AbbSum577, AbbSum426
	common /gzabbrev/ AbbSum622, AbbSum488, AbbSum658, AbbSum452
	common /gzabbrev/ AbbSum842, AbbSum642, AbbSum533, AbbSum430
	common /gzabbrev/ AbbSum611, AbbSum453, AbbSum387, AbbSum608
	common /gzabbrev/ AbbSum487, AbbSum843, AbbSum525, AbbSum526
	common /gzabbrev/ AbbSum612, AbbSum578, AbbSum486, AbbSum609
	common /gzabbrev/ AbbSum655, AbbSum388, AbbSum656, AbbSum800
	common /gzabbrev/ AbbSum293, AbbSum55, AbbSum901, AbbSum66
	common /gzabbrev/ AbbSum132, AbbSum292, AbbSum409, AbbSum390
	common /gzabbrev/ AbbSum428, AbbSum170, AbbSum332, AbbSum643
	common /gzabbrev/ AbbSum356, AbbSum357, AbbSum454, AbbSum484
	common /gzabbrev/ AbbSum657, AbbSum333, AbbSum644, AbbSum610
	common /gzabbrev/ AbbSum450, AbbSum283, AbbSum275, AbbSum359
	common /gzabbrev/ AbbSum330, AbbSum678, AbbSum451, AbbSum376
	common /gzabbrev/ AbbSum455, AbbSum389, AbbSum427, AbbSum844
	common /gzabbrev/ AbbSum764, AbbSum456, AbbSum579, AbbSum410
	common /gzabbrev/ AbbSum391, AbbSum429, AbbSum645, AbbSum673
	common /gzabbrev/ AbbSum607, AbbSum1050, AbbSum1044, AbbSum360
	common /gzabbrev/ AbbSum266, AbbSum361, AbbSum279, AbbSum854
	common /gzabbrev/ AbbSum100, AbbSum82, AbbSum101, AbbSum855
	common /gzabbrev/ AbbSum1074, AbbSum967, AbbSum948, AbbSum956
	common /gzabbrev/ AbbSum801, AbbSum464, AbbSum494, AbbSum957
	common /gzabbrev/ AbbSum40, AbbSum274, AbbSum1051, AbbSum616
	common /gzabbrev/ AbbSum1039, AbbSum968, AbbSum1026, AbbSum204
	common /gzabbrev/ AbbSum202, AbbSum958, AbbSum205, AbbSum803
	common /gzabbrev/ AbbSum1027, AbbSum171, AbbSum493, AbbSum619
	common /gzabbrev/ AbbSum281, AbbSum802, AbbSum131, AbbSum1144
	common /gzabbrev/ AbbSum42, AbbSum465, AbbSum282, AbbSum437
	common /gzabbrev/ AbbSum492, AbbSum675, AbbSum585, AbbSum659
	common /gzabbrev/ AbbSum618, AbbSum411, AbbSum334, AbbSum457
	common /gzabbrev/ AbbSum1038, AbbSum558, AbbSum569, AbbSum461
	common /gzabbrev/ AbbSum60, AbbSum723, AbbSum1109, AbbSum1102
	common /gzabbrev/ AbbSum1104, AbbSum872, AbbSum713, AbbSum1100
	common /gzabbrev/ AbbSum822, AbbSum1150, AbbSum63, AbbSum1152
	common /gzabbrev/ AbbSum1003, AbbSum1111, AbbSum932
	common /gzabbrev/ AbbSum1106, AbbSum973, AbbSum287, AbbSum9
	common /gzabbrev/ AbbSum15, AbbSum23, AbbSum724, AbbSum1
	common /gzabbrev/ AbbSum1148, AbbSum16, AbbSum715, AbbSum929
	common /gzabbrev/ AbbSum10, AbbSum21, AbbSum17, AbbSum714
	common /gzabbrev/ AbbSum931, AbbSum49, AbbSum725, AbbSum716
	common /gzabbrev/ AbbSum1113, AbbSum1142, AbbSum686, AbbSum699
	common /gzabbrev/ AbbSum717, AbbSum832, AbbSum737, AbbSum1073
	common /gzabbrev/ AbbSum1123, AbbSum738, AbbSum1124, AbbSum848
	common /gzabbrev/ AbbSum47, AbbSum890, AbbSum806, AbbSum48
	common /gzabbrev/ AbbSum974, AbbSum972, AbbSum685, AbbSum1058
	common /gzabbrev/ AbbSum869, AbbSum870, AbbSum868, AbbSum1014
	common /gzabbrev/ AbbSum1015, AbbSum782, AbbSum995, AbbSum994
	common /gzabbrev/ AbbSum781, AbbSum1013, AbbSum780, AbbSum813
	common /gzabbrev/ AbbSum752, AbbSum1136, AbbSum13, AbbSum1137
	common /gzabbrev/ AbbSum669, AbbSum1032, AbbSum25, AbbSum33
	common /gzabbrev/ AbbSum1012, AbbSum904, AbbSum823, AbbSum4
	common /gzabbrev/ AbbSum1004, AbbSum853, AbbSum26, AbbSum27
	common /gzabbrev/ AbbSum35, AbbSum735, AbbSum31, AbbSum1098
	common /gzabbrev/ AbbSum447, AbbSum707, AbbSum634, AbbSum586
	common /gzabbrev/ AbbSum323, AbbSum862, AbbSum629, AbbSum3
	common /gzabbrev/ AbbSum758, AbbSum653, AbbSum894, AbbSum635
	common /gzabbrev/ AbbSum630, AbbSum474, AbbSum475, AbbSum903
	common /gzabbrev/ AbbSum477, AbbSum349, AbbSum505, AbbSum930
	common /gzabbrev/ AbbSum667, AbbSum680, AbbSum636, AbbSum61
	common /gzabbrev/ AbbSum544, AbbSum546, AbbSum560, AbbSum639
	common /gzabbrev/ AbbSum732, AbbSum1101, AbbSum1116
	common /gzabbrev/ AbbSum1103, AbbSum1105, AbbSum687, AbbSum720
	common /gzabbrev/ AbbSum873, AbbSum975, AbbSum824, AbbSum290
	common /gzabbrev/ AbbSum289, AbbSum1151, AbbSum631, AbbSum1153
	common /gzabbrev/ AbbSum65, AbbSum1007, AbbSum1059, AbbSum1122
	common /gzabbrev/ AbbSum900, AbbSum626, AbbSum308, AbbSum238
	common /gzabbrev/ AbbSum1117, AbbSum627, AbbSum1107, AbbSum939
	common /gzabbrev/ AbbSum666, AbbSum625, AbbSum1076, AbbSum288
	common /gzabbrev/ AbbSum734, AbbSum891, AbbSum52, AbbSum11
	common /gzabbrev/ AbbSum730, AbbSum1128, AbbSum18, AbbSum24
	common /gzabbrev/ AbbSum733, AbbSum53, AbbSum1143, AbbSum811
	common /gzabbrev/ AbbSum977, AbbSum2, AbbSum690, AbbSum740
	common /gzabbrev/ AbbSum976, AbbSum1149, AbbSum1129, AbbSum850
	common /gzabbrev/ AbbSum51, AbbSum19, AbbSum702, AbbSum729
	common /gzabbrev/ AbbSum938, AbbSum12, AbbSum22, AbbSum20
	common /gzabbrev/ AbbSum750, AbbSum1084, AbbSum766, AbbSum1135
	common /gzabbrev/ AbbSum1057, AbbSum286, AbbSum514, AbbSum372
	common /gzabbrev/ AbbSum375, AbbSum591, AbbSum624, AbbSum307
	common /gzabbrev/ AbbSum441, AbbSum421, AbbSum1056, AbbSum291
	common /gzabbrev/ AbbSum239, AbbSum306, AbbSum889, AbbSum269
	common /gzabbrev/ AbbSum97, AbbSum491, AbbSum172, AbbSum668
	common /gzabbrev/ AbbSum137, AbbSum434, AbbSum677, AbbSum583
	common /gzabbrev/ AbbSum570, AbbSum584, AbbSum435, AbbSum436
	common /gzabbrev/ AbbSum568, AbbSum765, AbbSum1134, AbbSum665
	common /gzabbrev/ AbbSum504, AbbSum915, AbbSum676, AbbSum567
	common /gzabbrev/ AbbSum407, AbbSum1093, AbbSum880, AbbSum879
	common /gzabbrev/ AbbSum1021, AbbSum1033, AbbSum916
	common /gzabbrev/ AbbSum1094, AbbSum559, AbbSum420, AbbSum679
	common /gzabbrev/ AbbSum244, AbbSum592, AbbSum1022, AbbSum516
	common /gzabbrev/ AbbSum440, AbbSum245, AbbSum593, AbbSum545
	common /gzabbrev/ AbbSum408, AbbSum515, AbbSum1001, AbbSum790
	common /gzabbrev/ AbbSum485, AbbSum871, AbbSum833, AbbSum783
	common /gzabbrev/ AbbSum712, AbbSum543, AbbSum756, AbbSum996
	common /gzabbrev/ AbbSum1140, AbbSum1020, AbbSum721, AbbSum789
	common /gzabbrev/ AbbSum736, AbbSum876, AbbSum1115, AbbSum718
	common /gzabbrev/ AbbSum739, AbbSum816, AbbSum139, AbbSum532
	common /gzabbrev/ AbbSum243, AbbSum1112, AbbSum14, AbbSum28
	common /gzabbrev/ AbbSum34, AbbSum909, AbbSum6, AbbSum856
	common /gzabbrev/ AbbSum29, AbbSum30, AbbSum36, AbbSum32
	common /gzabbrev/ AbbSum1141, AbbSum1008, AbbSum672, AbbSum826
	common /gzabbrev/ AbbSum1019, AbbSum1110, AbbSum305
	common /gzabbrev/ AbbSum1062, AbbSum91, AbbSum1061, AbbSum572
	common /gzabbrev/ AbbSum674, AbbSum534, AbbSum557, AbbSum242
	common /gzabbrev/ AbbSum200, AbbSum240, AbbSum883, AbbSum1085
	common /gzabbrev/ AbbSum751, AbbSum902, AbbSum1125, AbbSum722
	common /gzabbrev/ AbbSum1099, AbbSum874, AbbSum814, AbbSum1095
	common /gzabbrev/ AbbSum807, AbbSum917, AbbSum466, AbbSum709
	common /gzabbrev/ AbbSum383, AbbSum949, AbbSum767, AbbSum637
	common /gzabbrev/ AbbSum1037, AbbSum597, AbbSum352, AbbSum327
	common /gzabbrev/ AbbSum324, AbbSum849, AbbSum753, AbbSum632
	common /gzabbrev/ AbbSum5, AbbSum759, AbbSum1052, AbbSum80
	common /gzabbrev/ AbbSum950, AbbSum875, AbbSum400, AbbSum660
	common /gzabbrev/ AbbSum401, AbbSum809, AbbSum785, AbbSum905
	common /gzabbrev/ AbbSum438, AbbSum815, AbbSum726, AbbSum638
	common /gzabbrev/ AbbSum1036, AbbSum1088, AbbSum795, AbbSum633
	common /gzabbrev/ AbbSum910, AbbSum794, AbbSum348, AbbSum965
	common /gzabbrev/ AbbSum497, AbbSum895, AbbSum743, AbbSum1078
	common /gzabbrev/ AbbSum757, AbbSum892, AbbSum745, AbbSum846
	common /gzabbrev/ AbbSum863, AbbSum777, AbbSum310, AbbSum498
	common /gzabbrev/ AbbSum476, AbbSum499, AbbSum44, AbbSum864
	common /gzabbrev/ AbbSum793, AbbSum893, AbbSum83, AbbSum1060
	common /gzabbrev/ AbbSum325, AbbSum792, AbbSum1079, AbbSum326
	common /gzabbrev/ AbbSum861, AbbSum865, AbbSum708, AbbSum744
	common /gzabbrev/ AbbSum966

	double complex cint1(3), cint2(3), cint3(3), cint4(3)
	double complex cint5(3), cint6(3), cint7(3), cint8(3)
	double complex cint9(3), cint10(3), cint11(3), cint12(3)
	double complex cint13(2), cint14(2), cint15(2), cint16(2)
	double complex cint17(2), cint18, cint19, cint20, cint21
	double complex cint22, cint23, cint24, cint25, cint26(3)
	double complex cint27(3), cint28(3), cint29(3), cint30(3)
	double complex cint31(3), cint32(3), cint33(3), cint34(3)
	double complex cint35(3), cint36(3), cint37(3), cint38(3)
	double complex cint39(3), cint40(3), cint41(2,2), cint42(3)
	double complex cint43(3), cint44(3), cint45(3), cint46(3)
	double complex cint47(3), cint48(2,2)
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
	integer*8 iint36(2), iint37(3), iint38(3), iint39(3), iint40
	integer*8 iint41, iint42(3), iint43(3), iint44(3), iint45(3)
	integer*8 iint46(3), iint47(3), iint48(3), iint49(3), iint50(3)
	integer*8 iint51(3), iint52(3), iint53(3), iint54(2,2)
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

	integer cha1, cha2, hia1, his1, lpd1, neu1, qud1, quu1, sld1
	integer sle1, sqd1, sqe1, squ1, sqv1
	common /gzindices/ cha1, cha2, hia1, his1, lpd1, neu1, qud1
	common /gzindices/ quu1, sld1, sle1, sqd1, sqe1, squ1, sqv1

	double complex Cloop(1)
	common /gzcoeff/ Cloop
