#include "model.h"
#include "util.h"
#include "looptools.h"
#include "renconst.h"

	double precision S, T, U
	common /kinvars/ S, T, U

	integer Hel(4)
	common /kinvars/ Hel

	double complex F1, F2, F3, F4, F12, F5, F6, F14, F7, F15, F9
	double complex F8, F16, F10, F11, F13, Pair1, Pair4, Pair5
	double complex Pair2, Pair3, Eps1, Eps2, Abb11, Abb12, Abb28
	double complex Abb29, Abb35, Abb36, Abb13, Abb14, Abb1, Abb27
	double complex Abb34, Abb2, Abb41, Abb42, Abb3, Abb7, Abb43
	double complex Abb44, Abb4, Abb8, Abb15, Abb5, Abb9, Abb16
	double complex Abb17, Abb30, Abb37, Abb18, Abb19, Abb31, Abb38
	double complex Abb20, Abb21, Abb6, Abb10, Abb22, Abb23, Abb32
	double complex Abb39, Abb24, Abb25, Abb33, Abb40, Abb26
	double complex AbbSum24, AbbSum25, AbbSum30, AbbSum88
	double complex AbbSum621, AbbSum723, AbbSum703, AbbSum1
	double complex AbbSum704, AbbSum224, AbbSum168, AbbSum225
	double complex AbbSum87, AbbSum65, AbbSum89, AbbSum101
	double complex AbbSum94, AbbSum625, AbbSum748, AbbSum31
	double complex AbbSum765, AbbSum754, AbbSum294, AbbSum214
	double complex AbbSum574, AbbSum55, AbbSum397, AbbSum123
	double complex AbbSum85, AbbSum93, AbbSum57, AbbSum41
	double complex AbbSum110, AbbSum129, AbbSum63, AbbSum402
	double complex AbbSum124, AbbSum56, AbbSum593, AbbSum667
	double complex AbbSum125, AbbSum152, AbbSum628, AbbSum36
	double complex AbbSum82, AbbSum117, AbbSum295, AbbSum106
	double complex AbbSum594, AbbSum292, AbbSum74, AbbSum113
	double complex AbbSum26, AbbSum262, AbbSum564, AbbSum143
	double complex AbbSum419, AbbSum282, AbbSum368, AbbSum227
	double complex AbbSum68, AbbSum665, AbbSum668, AbbSum158
	double complex AbbSum159, AbbSum127, AbbSum154, AbbSum162
	double complex AbbSum161, AbbSum157, AbbSum6, AbbSum5
	double complex AbbSum644, AbbSum111, AbbSum137, AbbSum126
	double complex AbbSum153, AbbSum396, AbbSum620, AbbSum546
	double complex AbbSum38, AbbSum84, AbbSum120, AbbSum317
	double complex AbbSum107, AbbSum609, AbbSum315, AbbSum79
	double complex AbbSum116, AbbSum27, AbbSum281, AbbSum576
	double complex AbbSum146, AbbSum28, AbbSum29, AbbSum155
	double complex AbbSum61, AbbSum109, AbbSum365, AbbSum211
	double complex AbbSum80, AbbSum367, AbbSum595, AbbSum91
	double complex AbbSum131, AbbSum122, AbbSum59, AbbSum133
	double complex AbbSum108, AbbSum90, AbbSum449, AbbSum135
	double complex AbbSum420, AbbSum283, AbbSum385, AbbSum248
	double complex AbbSum73, AbbSum666, AbbSum335, AbbSum369
	double complex AbbSum336, AbbSum58, AbbSum204, AbbSum533
	double complex AbbSum356, AbbSum16, AbbSum10, AbbSum767
	double complex AbbSum472, AbbSum64, AbbSum393, AbbSum66
	double complex AbbSum32, AbbSum545, AbbSum392, AbbSum86
	double complex AbbSum18, AbbSum72, AbbSum544, AbbSum95
	double complex AbbSum132, AbbSum713, AbbSum231, AbbSum712
	double complex AbbSum739, AbbSum639, AbbSum711, AbbSum637
	double complex AbbSum640, AbbSum300, AbbSum555, AbbSum686
	double complex AbbSum112, AbbSum624, AbbSum138, AbbSum14
	double complex AbbSum562, AbbSum610, AbbSum363, AbbSum401
	double complex AbbSum543, AbbSum445, AbbSum62, AbbSum60
	double complex AbbSum92, AbbSum142, AbbSum78, AbbSum130
	double complex AbbSum350, AbbSum156, AbbSum134, AbbSum151
	double complex AbbSum136, AbbSum588, AbbSum547, AbbSum575
	double complex AbbSum351, AbbSum128, AbbSum605, AbbSum457
	double complex AbbSum460, AbbSum184, AbbSum46, AbbSum405
	double complex AbbSum261, AbbSum458, AbbSum104, AbbSum75
	double complex AbbSum475, AbbSum477, AbbSum185, AbbSum48
	double complex AbbSum259, AbbSum388, AbbSum47, AbbSum37
	double complex AbbSum410, AbbSum280, AbbSum418, AbbSum215
	double complex AbbSum476, AbbSum416, AbbSum105, AbbSum404
	double complex AbbSum98, AbbSum327, AbbSum254, AbbSum202
	double complex AbbSum692, AbbSum689, AbbSum208, AbbSum352
	double complex AbbSum722, AbbSum726, AbbSum488, AbbSum707
	double complex AbbSum228, AbbSum566, AbbSum590, AbbSum165
	double complex AbbSum631, AbbSum774, AbbSum448, AbbSum768
	double complex AbbSum258, AbbSum461, AbbSum255, AbbSum421
	double complex AbbSum366, AbbSum523, AbbSum229, AbbSum424
	double complex AbbSum479, AbbSum487, AbbSum706, AbbSum370
	double complex AbbSum725, AbbSum186, AbbSum524, AbbSum182
	double complex AbbSum293, AbbSum50, AbbSum43, AbbSum591
	double complex AbbSum753, AbbSum567, AbbSum670, AbbSum592
	double complex AbbSum102, AbbSum565, AbbSum629, AbbSum381
	double complex AbbSum163, AbbSum650, AbbSum747, AbbSum534
	double complex AbbSum521, AbbSum19, AbbSum721, AbbSum456
	double complex AbbSum750, AbbSum511, AbbSum718, AbbSum249
	double complex AbbSum578, AbbSum606, AbbSum166, AbbSum646
	double complex AbbSum776, AbbSum450, AbbSum770, AbbSum627
	double complex AbbSum626, AbbSum193, AbbSum489, AbbSum513
	double complex AbbSum218, AbbSum724, AbbSum459, AbbSum118
	double complex AbbSum596, AbbSum632, AbbSum535, AbbSum649
	double complex AbbSum581, AbbSum775, AbbSum769, AbbSum279
	double complex AbbSum478, AbbSum256, AbbSum519, AbbSum323
	double complex AbbSum263, AbbSum563, AbbSum333, AbbSum558
	double complex AbbSum648, AbbSum422, AbbSum216, AbbSum260
	double complex AbbSum611, AbbSum384, AbbSum531, AbbSum250
	double complex AbbSum447, AbbSum480, AbbSum510, AbbSum717
	double complex AbbSum355, AbbSum394, AbbSum203, AbbSum486
	double complex AbbSum705, AbbSum386, AbbSum749, AbbSum187
	double complex AbbSum532, AbbSum183, AbbSum316, AbbSum53
	double complex AbbSum45, AbbSum607, AbbSum764, AbbSum579
	double complex AbbSum687, AbbSum608, AbbSum103, AbbSum577
	double complex AbbSum51, AbbSum525, AbbSum42, AbbSum40
	double complex AbbSum33, AbbSum195, AbbSum556, AbbSum613
	double complex AbbSum630, AbbSum580, AbbSum752, AbbSum568
	double complex AbbSum425, AbbSum67, AbbSum334, AbbSum645
	double complex AbbSum226, AbbSum387, AbbSum164, AbbSum664
	double complex AbbSum200, AbbSum210, AbbSum115, AbbSum669
	double complex AbbSum522, AbbSum548, AbbSum291, AbbSum683
	double complex AbbSum507, AbbSum272, AbbSum432, AbbSum301
	double complex AbbSum197, AbbSum756, AbbSum15, AbbSum710
	double complex AbbSum310, AbbSum35, AbbSum762, AbbSum408
	double complex AbbSum763, AbbSum244, AbbSum684, AbbSum235
	double complex AbbSum22, AbbSum309, AbbSum506, AbbSum305
	double complex AbbSum2, AbbSum243, AbbSum347, AbbSum440
	double complex AbbSum277, AbbSum411, AbbSum188, AbbSum76
	double complex AbbSum469, AbbSum174, AbbSum407, AbbSum360
	double complex AbbSum99, AbbSum743, AbbSum344, AbbSum77
	double complex AbbSum642, AbbSum437, AbbSum242, AbbSum237
	double complex AbbSum761, AbbSum716, AbbSum715, AbbSum671
	double complex AbbSum337, AbbSum399, AbbSum675, AbbSum463
	double complex AbbSum426, AbbSum643, AbbSum494, AbbSum493
	double complex AbbSum622, AbbSum428, AbbSum201, AbbSum720
	double complex AbbSum199, AbbSum504, AbbSum439, AbbSum505
	double complex AbbSum308, AbbSum714, AbbSum342, AbbSum196
	double complex AbbSum212, AbbSum307, AbbSum378, AbbSum341
	double complex AbbSum695, AbbSum71, AbbSum678, AbbSum207
	double complex AbbSum213, AbbSum169, AbbSum198, AbbSum251
	double complex AbbSum252, AbbSum751, AbbSum647, AbbSum400
	double complex AbbSum659, AbbSum638, AbbSum176, AbbSum677
	double complex AbbSum232, AbbSum233, AbbSum234, AbbSum236
	double complex AbbSum273, AbbSum274, AbbSum318, AbbSum529
	double complex AbbSum673, AbbSum346, AbbSum343, AbbSum614
	double complex AbbSum206, AbbSum339, AbbSum438, AbbSum69
	double complex AbbSum550, AbbSum326, AbbSum498, AbbSum786
	double complex AbbSum398, AbbSum230, AbbSum302, AbbSum436
	double complex AbbSum376, AbbSum429, AbbSum597, AbbSum430
	double complex AbbSum431, AbbSum435, AbbSum755, AbbSum481
	double complex AbbSum758, AbbSum298, AbbSum495, AbbSum760
	double complex AbbSum685, AbbSum508, AbbSum70, AbbSum499
	double complex AbbSum121, AbbSum239, AbbSum359, AbbSum276
	double complex AbbSum278, AbbSum23, AbbSum96, AbbSum296
	double complex AbbSum766, AbbSum427, AbbSum634, AbbSum491
	double complex AbbSum181, AbbSum782, AbbSum172, AbbSum240
	double complex AbbSum377, AbbSum379, AbbSum303, AbbSum17
	double complex AbbSum719, AbbSum538, AbbSum144, AbbSum299
	double complex AbbSum496, AbbSum757, AbbSum271, AbbSum657
	double complex AbbSum709, AbbSum433, AbbSum177, AbbSum497
	double complex AbbSum672, AbbSum681, AbbSum340, AbbSum238
	double complex AbbSum345, AbbSum328, AbbSum542, AbbSum759
	double complex AbbSum708, AbbSum391, AbbSum679, AbbSum160
	double complex AbbSum297, AbbSum490, AbbSum514, AbbSum676
	double complex AbbSum319, AbbSum539, AbbSum583, AbbSum674
	double complex AbbSum492, AbbSum785, AbbSum304, AbbSum551
	double complex AbbSum442, AbbSum552, AbbSum623, AbbSum338
	double complex AbbSum549, AbbSum321, AbbSum253, AbbSum512
	double complex AbbSum119, AbbSum361, AbbSum362, AbbSum209
	double complex AbbSum471, AbbSum554, AbbSum349, AbbSum97
	double complex AbbSum443, AbbSum175, AbbSum619, AbbSum682
	double complex AbbSum382, AbbSum383, AbbSum348, AbbSum171
	double complex AbbSum584, AbbSum357, AbbSum444, AbbSum170
	double complex AbbSum320, AbbSum553, AbbSum680, AbbSum306
	double complex AbbSum3, AbbSum145, AbbSum173, AbbSum658
	double complex AbbSum205, AbbSum441, AbbSum241, AbbSum265
	double complex AbbSum641, AbbSum434, AbbSum598, AbbSum585
	double complex AbbSum557, AbbSum536, AbbSum20, AbbSum8
	double complex AbbSum9, AbbSum7, AbbSum329, AbbSum21
	double complex AbbSum699, AbbSum617, AbbSum740, AbbSum569
	double complex AbbSum732, AbbSum517, AbbSum660, AbbSum465
	double complex AbbSum772, AbbSum700, AbbSum559, AbbSum771
	double complex AbbSum192, AbbSum139, AbbSum582, AbbSum616
	double complex AbbSum371, AbbSum372, AbbSum618, AbbSum733
	double complex AbbSum635, AbbSum734, AbbSum178, AbbSum702
	double complex AbbSum697, AbbSum180, AbbSum390, AbbSum516
	double complex AbbSum787, AbbSum586, AbbSum13, AbbSum11
	double complex AbbSum12, AbbSum735, AbbSum736, AbbSum742
	double complex AbbSum167, AbbSum264, AbbSum652, AbbSum633
	double complex AbbSum654, AbbSum730, AbbSum179, AbbSum503
	double complex AbbSum784, AbbSum415, AbbSum409, AbbSum587
	double complex AbbSum44, AbbSum570, AbbSum189, AbbSum612
	double complex AbbSum466, AbbSum270, AbbSum114, AbbSum781
	double complex AbbSum515, AbbSum417, AbbSum737, AbbSum527
	double complex AbbSum526, AbbSum467, AbbSum636, AbbSum741
	double complex AbbSum530, AbbSum275, AbbSum358, AbbSum738
	double complex AbbSum470, AbbSum602, AbbSum140, AbbSum777
	double complex AbbSum655, AbbSum778, AbbSum731, AbbSum656
	double complex AbbSum100, AbbSum268, AbbSum190, AbbSum413
	double complex AbbSum412, AbbSum560, AbbSum81, AbbSum500
	double complex AbbSum269, AbbSum389, AbbSum600, AbbSum502
	double complex AbbSum501, AbbSum414, AbbSum464, AbbSum406
	double complex AbbSum661, AbbSum701, AbbSum690, AbbSum468
	double complex AbbSum653, AbbSum662, AbbSum573, AbbSum325
	double complex AbbSum267, AbbSum324, AbbSum603, AbbSum604
	double complex AbbSum571, AbbSum83, AbbSum599, AbbSum528
	double complex AbbSum364, AbbSum540, AbbSum266, AbbSum589
	double complex AbbSum373, AbbSum374, AbbSum572, AbbSum601
	double complex AbbSum375, AbbSum150, AbbSum148, AbbSum693
	double complex AbbSum353, AbbSum541, AbbSum615, AbbSum537
	double complex AbbSum779, AbbSum729, AbbSum191, AbbSum773
	double complex AbbSum147, AbbSum149, AbbSum694, AbbSum141
	double complex AbbSum663, AbbSum783, AbbSum696, AbbSum746
	double complex AbbSum744, AbbSum395, AbbSum330, AbbSum520
	double complex AbbSum331, AbbSum780, AbbSum54, AbbSum561
	double complex AbbSum380, AbbSum745, AbbSum217, AbbSum332
	double complex AbbSum727, AbbSum462, AbbSum651, AbbSum728
	double complex AbbSum698, AbbSum451, AbbSum423, AbbSum49
	double complex AbbSum285, AbbSum220, AbbSum287, AbbSum288
	double complex AbbSum484, AbbSum290, AbbSum222, AbbSum223
	double complex AbbSum455, AbbSum691, AbbSum219, AbbSum518
	double complex AbbSum473, AbbSum194, AbbSum688, AbbSum446
	double complex AbbSum284, AbbSum322, AbbSum52, AbbSum39
	double complex AbbSum354, AbbSum403, AbbSum311, AbbSum452
	double complex AbbSum245, AbbSum286, AbbSum312, AbbSum482
	double complex AbbSum313, AbbSum483, AbbSum509, AbbSum257
	double complex AbbSum289, AbbSum485, AbbSum314, AbbSum221
	double complex AbbSum246, AbbSum247, AbbSum4, AbbSum34
	double complex AbbSum453, AbbSum454, AbbSum474
	common /abbrev/ F1, F2, F3, F4, F12, F5, F6, F14, F7, F15, F9
	common /abbrev/ F8, F16, F10, F11, F13, Pair1, Pair4, Pair5
	common /abbrev/ Pair2, Pair3, Eps1, Eps2, Abb11, Abb12, Abb28
	common /abbrev/ Abb29, Abb35, Abb36, Abb13, Abb14, Abb1, Abb27
	common /abbrev/ Abb34, Abb2, Abb41, Abb42, Abb3, Abb7, Abb43
	common /abbrev/ Abb44, Abb4, Abb8, Abb15, Abb5, Abb9, Abb16
	common /abbrev/ Abb17, Abb30, Abb37, Abb18, Abb19, Abb31
	common /abbrev/ Abb38, Abb20, Abb21, Abb6, Abb10, Abb22, Abb23
	common /abbrev/ Abb32, Abb39, Abb24, Abb25, Abb33, Abb40
	common /abbrev/ Abb26, AbbSum24, AbbSum25, AbbSum30, AbbSum88
	common /abbrev/ AbbSum621, AbbSum723, AbbSum703, AbbSum1
	common /abbrev/ AbbSum704, AbbSum224, AbbSum168, AbbSum225
	common /abbrev/ AbbSum87, AbbSum65, AbbSum89, AbbSum101
	common /abbrev/ AbbSum94, AbbSum625, AbbSum748, AbbSum31
	common /abbrev/ AbbSum765, AbbSum754, AbbSum294, AbbSum214
	common /abbrev/ AbbSum574, AbbSum55, AbbSum397, AbbSum123
	common /abbrev/ AbbSum85, AbbSum93, AbbSum57, AbbSum41
	common /abbrev/ AbbSum110, AbbSum129, AbbSum63, AbbSum402
	common /abbrev/ AbbSum124, AbbSum56, AbbSum593, AbbSum667
	common /abbrev/ AbbSum125, AbbSum152, AbbSum628, AbbSum36
	common /abbrev/ AbbSum82, AbbSum117, AbbSum295, AbbSum106
	common /abbrev/ AbbSum594, AbbSum292, AbbSum74, AbbSum113
	common /abbrev/ AbbSum26, AbbSum262, AbbSum564, AbbSum143
	common /abbrev/ AbbSum419, AbbSum282, AbbSum368, AbbSum227
	common /abbrev/ AbbSum68, AbbSum665, AbbSum668, AbbSum158
	common /abbrev/ AbbSum159, AbbSum127, AbbSum154, AbbSum162
	common /abbrev/ AbbSum161, AbbSum157, AbbSum6, AbbSum5
	common /abbrev/ AbbSum644, AbbSum111, AbbSum137, AbbSum126
	common /abbrev/ AbbSum153, AbbSum396, AbbSum620, AbbSum546
	common /abbrev/ AbbSum38, AbbSum84, AbbSum120, AbbSum317
	common /abbrev/ AbbSum107, AbbSum609, AbbSum315, AbbSum79
	common /abbrev/ AbbSum116, AbbSum27, AbbSum281, AbbSum576
	common /abbrev/ AbbSum146, AbbSum28, AbbSum29, AbbSum155
	common /abbrev/ AbbSum61, AbbSum109, AbbSum365, AbbSum211
	common /abbrev/ AbbSum80, AbbSum367, AbbSum595, AbbSum91
	common /abbrev/ AbbSum131, AbbSum122, AbbSum59, AbbSum133
	common /abbrev/ AbbSum108, AbbSum90, AbbSum449, AbbSum135
	common /abbrev/ AbbSum420, AbbSum283, AbbSum385, AbbSum248
	common /abbrev/ AbbSum73, AbbSum666, AbbSum335, AbbSum369
	common /abbrev/ AbbSum336, AbbSum58, AbbSum204, AbbSum533
	common /abbrev/ AbbSum356, AbbSum16, AbbSum10, AbbSum767
	common /abbrev/ AbbSum472, AbbSum64, AbbSum393, AbbSum66
	common /abbrev/ AbbSum32, AbbSum545, AbbSum392, AbbSum86
	common /abbrev/ AbbSum18, AbbSum72, AbbSum544, AbbSum95
	common /abbrev/ AbbSum132, AbbSum713, AbbSum231, AbbSum712
	common /abbrev/ AbbSum739, AbbSum639, AbbSum711, AbbSum637
	common /abbrev/ AbbSum640, AbbSum300, AbbSum555, AbbSum686
	common /abbrev/ AbbSum112, AbbSum624, AbbSum138, AbbSum14
	common /abbrev/ AbbSum562, AbbSum610, AbbSum363, AbbSum401
	common /abbrev/ AbbSum543, AbbSum445, AbbSum62, AbbSum60
	common /abbrev/ AbbSum92, AbbSum142, AbbSum78, AbbSum130
	common /abbrev/ AbbSum350, AbbSum156, AbbSum134, AbbSum151
	common /abbrev/ AbbSum136, AbbSum588, AbbSum547, AbbSum575
	common /abbrev/ AbbSum351, AbbSum128, AbbSum605, AbbSum457
	common /abbrev/ AbbSum460, AbbSum184, AbbSum46, AbbSum405
	common /abbrev/ AbbSum261, AbbSum458, AbbSum104, AbbSum75
	common /abbrev/ AbbSum475, AbbSum477, AbbSum185, AbbSum48
	common /abbrev/ AbbSum259, AbbSum388, AbbSum47, AbbSum37
	common /abbrev/ AbbSum410, AbbSum280, AbbSum418, AbbSum215
	common /abbrev/ AbbSum476, AbbSum416, AbbSum105, AbbSum404
	common /abbrev/ AbbSum98, AbbSum327, AbbSum254, AbbSum202
	common /abbrev/ AbbSum692, AbbSum689, AbbSum208, AbbSum352
	common /abbrev/ AbbSum722, AbbSum726, AbbSum488, AbbSum707
	common /abbrev/ AbbSum228, AbbSum566, AbbSum590, AbbSum165
	common /abbrev/ AbbSum631, AbbSum774, AbbSum448, AbbSum768
	common /abbrev/ AbbSum258, AbbSum461, AbbSum255, AbbSum421
	common /abbrev/ AbbSum366, AbbSum523, AbbSum229, AbbSum424
	common /abbrev/ AbbSum479, AbbSum487, AbbSum706, AbbSum370
	common /abbrev/ AbbSum725, AbbSum186, AbbSum524, AbbSum182
	common /abbrev/ AbbSum293, AbbSum50, AbbSum43, AbbSum591
	common /abbrev/ AbbSum753, AbbSum567, AbbSum670, AbbSum592
	common /abbrev/ AbbSum102, AbbSum565, AbbSum629, AbbSum381
	common /abbrev/ AbbSum163, AbbSum650, AbbSum747, AbbSum534
	common /abbrev/ AbbSum521, AbbSum19, AbbSum721, AbbSum456
	common /abbrev/ AbbSum750, AbbSum511, AbbSum718, AbbSum249
	common /abbrev/ AbbSum578, AbbSum606, AbbSum166, AbbSum646
	common /abbrev/ AbbSum776, AbbSum450, AbbSum770, AbbSum627
	common /abbrev/ AbbSum626, AbbSum193, AbbSum489, AbbSum513
	common /abbrev/ AbbSum218, AbbSum724, AbbSum459, AbbSum118
	common /abbrev/ AbbSum596, AbbSum632, AbbSum535, AbbSum649
	common /abbrev/ AbbSum581, AbbSum775, AbbSum769, AbbSum279
	common /abbrev/ AbbSum478, AbbSum256, AbbSum519, AbbSum323
	common /abbrev/ AbbSum263, AbbSum563, AbbSum333, AbbSum558
	common /abbrev/ AbbSum648, AbbSum422, AbbSum216, AbbSum260
	common /abbrev/ AbbSum611, AbbSum384, AbbSum531, AbbSum250
	common /abbrev/ AbbSum447, AbbSum480, AbbSum510, AbbSum717
	common /abbrev/ AbbSum355, AbbSum394, AbbSum203, AbbSum486
	common /abbrev/ AbbSum705, AbbSum386, AbbSum749, AbbSum187
	common /abbrev/ AbbSum532, AbbSum183, AbbSum316, AbbSum53
	common /abbrev/ AbbSum45, AbbSum607, AbbSum764, AbbSum579
	common /abbrev/ AbbSum687, AbbSum608, AbbSum103, AbbSum577
	common /abbrev/ AbbSum51, AbbSum525, AbbSum42, AbbSum40
	common /abbrev/ AbbSum33, AbbSum195, AbbSum556, AbbSum613
	common /abbrev/ AbbSum630, AbbSum580, AbbSum752, AbbSum568
	common /abbrev/ AbbSum425, AbbSum67, AbbSum334, AbbSum645
	common /abbrev/ AbbSum226, AbbSum387, AbbSum164, AbbSum664
	common /abbrev/ AbbSum200, AbbSum210, AbbSum115, AbbSum669
	common /abbrev/ AbbSum522, AbbSum548, AbbSum291, AbbSum683
	common /abbrev/ AbbSum507, AbbSum272, AbbSum432, AbbSum301
	common /abbrev/ AbbSum197, AbbSum756, AbbSum15, AbbSum710
	common /abbrev/ AbbSum310, AbbSum35, AbbSum762, AbbSum408
	common /abbrev/ AbbSum763, AbbSum244, AbbSum684, AbbSum235
	common /abbrev/ AbbSum22, AbbSum309, AbbSum506, AbbSum305
	common /abbrev/ AbbSum2, AbbSum243, AbbSum347, AbbSum440
	common /abbrev/ AbbSum277, AbbSum411, AbbSum188, AbbSum76
	common /abbrev/ AbbSum469, AbbSum174, AbbSum407, AbbSum360
	common /abbrev/ AbbSum99, AbbSum743, AbbSum344, AbbSum77
	common /abbrev/ AbbSum642, AbbSum437, AbbSum242, AbbSum237
	common /abbrev/ AbbSum761, AbbSum716, AbbSum715, AbbSum671
	common /abbrev/ AbbSum337, AbbSum399, AbbSum675, AbbSum463
	common /abbrev/ AbbSum426, AbbSum643, AbbSum494, AbbSum493
	common /abbrev/ AbbSum622, AbbSum428, AbbSum201, AbbSum720
	common /abbrev/ AbbSum199, AbbSum504, AbbSum439, AbbSum505
	common /abbrev/ AbbSum308, AbbSum714, AbbSum342, AbbSum196
	common /abbrev/ AbbSum212, AbbSum307, AbbSum378, AbbSum341
	common /abbrev/ AbbSum695, AbbSum71, AbbSum678, AbbSum207
	common /abbrev/ AbbSum213, AbbSum169, AbbSum198, AbbSum251
	common /abbrev/ AbbSum252, AbbSum751, AbbSum647, AbbSum400
	common /abbrev/ AbbSum659, AbbSum638, AbbSum176, AbbSum677
	common /abbrev/ AbbSum232, AbbSum233, AbbSum234, AbbSum236
	common /abbrev/ AbbSum273, AbbSum274, AbbSum318, AbbSum529
	common /abbrev/ AbbSum673, AbbSum346, AbbSum343, AbbSum614
	common /abbrev/ AbbSum206, AbbSum339, AbbSum438, AbbSum69
	common /abbrev/ AbbSum550, AbbSum326, AbbSum498, AbbSum786
	common /abbrev/ AbbSum398, AbbSum230, AbbSum302, AbbSum436
	common /abbrev/ AbbSum376, AbbSum429, AbbSum597, AbbSum430
	common /abbrev/ AbbSum431, AbbSum435, AbbSum755, AbbSum481
	common /abbrev/ AbbSum758, AbbSum298, AbbSum495, AbbSum760
	common /abbrev/ AbbSum685, AbbSum508, AbbSum70, AbbSum499
	common /abbrev/ AbbSum121, AbbSum239, AbbSum359, AbbSum276
	common /abbrev/ AbbSum278, AbbSum23, AbbSum96, AbbSum296
	common /abbrev/ AbbSum766, AbbSum427, AbbSum634, AbbSum491
	common /abbrev/ AbbSum181, AbbSum782, AbbSum172, AbbSum240
	common /abbrev/ AbbSum377, AbbSum379, AbbSum303, AbbSum17
	common /abbrev/ AbbSum719, AbbSum538, AbbSum144, AbbSum299
	common /abbrev/ AbbSum496, AbbSum757, AbbSum271, AbbSum657
	common /abbrev/ AbbSum709, AbbSum433, AbbSum177, AbbSum497
	common /abbrev/ AbbSum672, AbbSum681, AbbSum340, AbbSum238
	common /abbrev/ AbbSum345, AbbSum328, AbbSum542, AbbSum759
	common /abbrev/ AbbSum708, AbbSum391, AbbSum679, AbbSum160
	common /abbrev/ AbbSum297, AbbSum490, AbbSum514, AbbSum676
	common /abbrev/ AbbSum319, AbbSum539, AbbSum583, AbbSum674
	common /abbrev/ AbbSum492, AbbSum785, AbbSum304, AbbSum551
	common /abbrev/ AbbSum442, AbbSum552, AbbSum623, AbbSum338
	common /abbrev/ AbbSum549, AbbSum321, AbbSum253, AbbSum512
	common /abbrev/ AbbSum119, AbbSum361, AbbSum362, AbbSum209
	common /abbrev/ AbbSum471, AbbSum554, AbbSum349, AbbSum97
	common /abbrev/ AbbSum443, AbbSum175, AbbSum619, AbbSum682
	common /abbrev/ AbbSum382, AbbSum383, AbbSum348, AbbSum171
	common /abbrev/ AbbSum584, AbbSum357, AbbSum444, AbbSum170
	common /abbrev/ AbbSum320, AbbSum553, AbbSum680, AbbSum306
	common /abbrev/ AbbSum3, AbbSum145, AbbSum173, AbbSum658
	common /abbrev/ AbbSum205, AbbSum441, AbbSum241, AbbSum265
	common /abbrev/ AbbSum641, AbbSum434, AbbSum598, AbbSum585
	common /abbrev/ AbbSum557, AbbSum536, AbbSum20, AbbSum8
	common /abbrev/ AbbSum9, AbbSum7, AbbSum329, AbbSum21
	common /abbrev/ AbbSum699, AbbSum617, AbbSum740, AbbSum569
	common /abbrev/ AbbSum732, AbbSum517, AbbSum660, AbbSum465
	common /abbrev/ AbbSum772, AbbSum700, AbbSum559, AbbSum771
	common /abbrev/ AbbSum192, AbbSum139, AbbSum582, AbbSum616
	common /abbrev/ AbbSum371, AbbSum372, AbbSum618, AbbSum733
	common /abbrev/ AbbSum635, AbbSum734, AbbSum178, AbbSum702
	common /abbrev/ AbbSum697, AbbSum180, AbbSum390, AbbSum516
	common /abbrev/ AbbSum787, AbbSum586, AbbSum13, AbbSum11
	common /abbrev/ AbbSum12, AbbSum735, AbbSum736, AbbSum742
	common /abbrev/ AbbSum167, AbbSum264, AbbSum652, AbbSum633
	common /abbrev/ AbbSum654, AbbSum730, AbbSum179, AbbSum503
	common /abbrev/ AbbSum784, AbbSum415, AbbSum409, AbbSum587
	common /abbrev/ AbbSum44, AbbSum570, AbbSum189, AbbSum612
	common /abbrev/ AbbSum466, AbbSum270, AbbSum114, AbbSum781
	common /abbrev/ AbbSum515, AbbSum417, AbbSum737, AbbSum527
	common /abbrev/ AbbSum526, AbbSum467, AbbSum636, AbbSum741
	common /abbrev/ AbbSum530, AbbSum275, AbbSum358, AbbSum738
	common /abbrev/ AbbSum470, AbbSum602, AbbSum140, AbbSum777
	common /abbrev/ AbbSum655, AbbSum778, AbbSum731, AbbSum656
	common /abbrev/ AbbSum100, AbbSum268, AbbSum190, AbbSum413
	common /abbrev/ AbbSum412, AbbSum560, AbbSum81, AbbSum500
	common /abbrev/ AbbSum269, AbbSum389, AbbSum600, AbbSum502
	common /abbrev/ AbbSum501, AbbSum414, AbbSum464, AbbSum406
	common /abbrev/ AbbSum661, AbbSum701, AbbSum690, AbbSum468
	common /abbrev/ AbbSum653, AbbSum662, AbbSum573, AbbSum325
	common /abbrev/ AbbSum267, AbbSum324, AbbSum603, AbbSum604
	common /abbrev/ AbbSum571, AbbSum83, AbbSum599, AbbSum528
	common /abbrev/ AbbSum364, AbbSum540, AbbSum266, AbbSum589
	common /abbrev/ AbbSum373, AbbSum374, AbbSum572, AbbSum601
	common /abbrev/ AbbSum375, AbbSum150, AbbSum148, AbbSum693
	common /abbrev/ AbbSum353, AbbSum541, AbbSum615, AbbSum537
	common /abbrev/ AbbSum779, AbbSum729, AbbSum191, AbbSum773
	common /abbrev/ AbbSum147, AbbSum149, AbbSum694, AbbSum141
	common /abbrev/ AbbSum663, AbbSum783, AbbSum696, AbbSum746
	common /abbrev/ AbbSum744, AbbSum395, AbbSum330, AbbSum520
	common /abbrev/ AbbSum331, AbbSum780, AbbSum54, AbbSum561
	common /abbrev/ AbbSum380, AbbSum745, AbbSum217, AbbSum332
	common /abbrev/ AbbSum727, AbbSum462, AbbSum651, AbbSum728
	common /abbrev/ AbbSum698, AbbSum451, AbbSum423, AbbSum49
	common /abbrev/ AbbSum285, AbbSum220, AbbSum287, AbbSum288
	common /abbrev/ AbbSum484, AbbSum290, AbbSum222, AbbSum223
	common /abbrev/ AbbSum455, AbbSum691, AbbSum219, AbbSum518
	common /abbrev/ AbbSum473, AbbSum194, AbbSum688, AbbSum446
	common /abbrev/ AbbSum284, AbbSum322, AbbSum52, AbbSum39
	common /abbrev/ AbbSum354, AbbSum403, AbbSum311, AbbSum452
	common /abbrev/ AbbSum245, AbbSum286, AbbSum312, AbbSum482
	common /abbrev/ AbbSum313, AbbSum483, AbbSum509, AbbSum257
	common /abbrev/ AbbSum289, AbbSum485, AbbSum314, AbbSum221
	common /abbrev/ AbbSum246, AbbSum247, AbbSum4, AbbSum34
	common /abbrev/ AbbSum453, AbbSum454, AbbSum474

	double complex cint1, cint2, cint3, cint4, cint5, cint6
	double complex cint7(3), cint8(3), cint9(3), cint10(3)
	double complex cint11(3), cint12(3), cint13(3), cint14(3)
	double complex cint15(3), cint16(3), cint17(3), cint18(3)
	double complex cint19(3), cint20(3), cint21(3), cint22(3)
	double complex cint23(3), cint24(3), cint25(2), cint26(2)
	common /loopint/ cint1, cint2, cint3, cint4, cint5, cint6
	common /loopint/ cint7, cint8, cint9, cint10, cint11, cint12
	common /loopint/ cint13, cint14, cint15, cint16, cint17
	common /loopint/ cint18, cint19, cint20, cint21, cint22
	common /loopint/ cint23, cint24, cint25, cint26

	integer*8 iint1(3), iint2(3), iint3(3), iint4(2), iint5, iint6
	integer*8 iint7(3), iint8(3), iint9(3), iint10(3), iint11(3)
	integer*8 iint12(3), iint13(3), iint14(3), iint15(3), iint16(3)
	integer*8 iint17(3), iint18(3), iint19(3), iint20(3), iint21(3)
	integer*8 iint22(3), iint23(3), iint24(3), iint25(3), iint26(3)
	integer*8 iint27(3), iint28(3), iint29(3), iint30(3), iint31(3)
	integer*8 iint32(3), iint33(3), iint34(3), iint35(3), iint36(3)
	integer*8 iint37(3), iint38(3), iint39(3), iint40(3), iint41(3)
	integer*8 iint42(3), iint43(3), iint44(3), iint45(3), iint46(3)
	integer*8 iint47(3), iint48(3), iint49(3), iint50(3), iint51(3)
	integer*8 iint52(3), iint53(3), iint54(3), iint55(3), iint56(3)
	integer*8 iint57(3), iint58(3), iint59(3), iint60(3), iint61(3)
	integer*8 iint62(3), iint63(3), iint64(3), iint65(3), iint66(3)
	integer*8 iint67(2), iint68(2), iint69(2), iint70(2), iint71(2)
	integer*8 iint72(2), iint73(2), iint74(2), iint75(2), iint76(2)
	integer*8 iint77(2), iint78(2), iint79(2), iint80(2), iint81(2)
	integer*8 iint82(2), iint83(2), iint84(2), iint85(2), iint86(2)
	integer*8 iint87(2), iint88(2), iint89(2), iint90(2), iint91(2)
	integer*8 iint92(2), iint93(2), iint94(2), iint95(2)
	common /loopint/ iint1, iint2, iint3, iint4, iint5, iint6
	common /loopint/ iint7, iint8, iint9, iint10, iint11, iint12
	common /loopint/ iint13, iint14, iint15, iint16, iint17
	common /loopint/ iint18, iint19, iint20, iint21, iint22
	common /loopint/ iint23, iint24, iint25, iint26, iint27
	common /loopint/ iint28, iint29, iint30, iint31, iint32
	common /loopint/ iint33, iint34, iint35, iint36, iint37
	common /loopint/ iint38, iint39, iint40, iint41, iint42
	common /loopint/ iint43, iint44, iint45, iint46, iint47
	common /loopint/ iint48, iint49, iint50, iint51, iint52
	common /loopint/ iint53, iint54, iint55, iint56, iint57
	common /loopint/ iint58, iint59, iint60, iint61, iint62
	common /loopint/ iint63, iint64, iint65, iint66, iint67
	common /loopint/ iint68, iint69, iint70, iint71, iint72
	common /loopint/ iint73, iint74, iint75, iint76, iint77
	common /loopint/ iint78, iint79, iint80, iint81, iint82
	common /loopint/ iint83, iint84, iint85, iint86, iint87
	common /loopint/ iint88, iint89, iint90, iint91, iint92
	common /loopint/ iint93, iint94, iint95

	integer cha1, his1, lpd1, qud1, quu1, sld1, sle1, sqd1, sqe1
	integer squ1, sqv1
	common /indices/ cha1, his1, lpd1, qud1, quu1, sld1, sle1
	common /indices/ sqd1, sqe1, squ1, sqv1

	double complex Cloop(1)
	common /coeff/ Cloop
