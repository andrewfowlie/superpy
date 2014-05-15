#include "model.h"
#include "util.h"
#include "looptools.h"
#include "renconst.h"

	double precision S, T, U
	common /kinvars/ S, T, U

	integer Hel(4)
	common /kinvars/ Hel

	double complex F1, F2, F3, F4, F12, F5, F6, F14, F7, F15, F9
	double complex F8, F16, F10, F11, F13, Pair1, Pair2, Pair5
	double complex Pair3, Pair4, Eps1, Eps2, Abb11, Abb12, Abb28
	double complex Abb29, Abb35, Abb36, Abb13, Abb14, Abb1, Abb27
	double complex Abb34, Abb2, Abb23, Abb15, Abb19, Abb24, Abb41
	double complex Abb42, Abb16, Abb20, Abb3, Abb30, Abb37, Abb7
	double complex Abb43, Abb44, Abb17, Abb21, Abb4, Abb31, Abb38
	double complex Abb8, Abb25, Abb18, Abb22, Abb26, Abb5, Abb32
	double complex Abb39, Abb9, Abb6, Abb33, Abb40, Abb10, AbbSum3
	double complex AbbSum4, AbbSum394, AbbSum252, AbbSum196
	double complex AbbSum19, AbbSum343, AbbSum26, AbbSum6
	double complex AbbSum220, AbbSum119, AbbSum332, AbbSum539
	double complex AbbSum5, AbbSum540, AbbSum333, AbbSum168
	double complex AbbSum188, AbbSum70, AbbSum201, AbbSum562
	double complex AbbSum171, AbbSum298, AbbSum181, AbbSum497
	double complex AbbSum399, AbbSum427, AbbSum517, AbbSum259
	double complex AbbSum495, AbbSum219, AbbSum415, AbbSum221
	double complex AbbSum251, AbbSum313, AbbSum253, AbbSum203
	double complex AbbSum202, AbbSum296, AbbSum336, AbbSum410
	double complex AbbSum144, AbbSum338, AbbSum277, AbbSum205
	double complex AbbSum409, AbbSum276, AbbSum493, AbbSum195
	double complex AbbSum416, AbbSum76, AbbSum254, AbbSum492
	double complex AbbSum69, AbbSum170, AbbSum461, AbbSum198
	double complex AbbSum293, AbbSum275, AbbSum71, AbbSum21
	double complex AbbSum235, AbbSum145, AbbSum462, AbbSum224
	double complex AbbSum264, AbbSum339, AbbSum411, AbbSum297
	double complex AbbSum274, AbbSum272, AbbSum385, AbbSum197
	double complex AbbSum273, AbbSum312, AbbSum294, AbbSum496
	double complex AbbSum223, AbbSum344, AbbSum160, AbbSum35
	double complex AbbSum346, AbbSum97, AbbSum99, AbbSum121
	double complex AbbSum310, AbbSum27, AbbSum96, AbbSum255
	double complex AbbSum494, AbbSum488, AbbSum222, AbbSum337
	double complex AbbSum169, AbbSum7, AbbSum173, AbbSum341
	double complex AbbSum200, AbbSum295, AbbSum77, AbbSum311
	double complex AbbSum256, AbbSum453, AbbSum340, AbbSum227
	double complex AbbSum518, AbbSum521, AbbSum546, AbbSum20
	double complex AbbSum14, AbbSum479, AbbSum545, AbbSum478
	double complex AbbSum342, AbbSum489, AbbSum371, AbbSum523
	double complex AbbSum473, AbbSum472, AbbSum240, AbbSum248
	double complex AbbSum249, AbbSum372, AbbSum236, AbbSum369
	double complex AbbSum544, AbbSum218, AbbSum519, AbbSum175
	double complex AbbSum520, AbbSum566, AbbSum401, AbbSum524
	double complex AbbSum543, AbbSum217, AbbSum474, AbbSum374
	double complex AbbSum373, AbbSum490, AbbSum542, AbbSum491
	double complex AbbSum334, AbbSum375, AbbSum285, AbbSum250
	double complex AbbSum335, AbbSum324, AbbSum443, AbbSum377
	double complex AbbSum549, AbbSum484, AbbSum358, AbbSum475
	double complex AbbSum241, AbbSum172, AbbSum176, AbbSum174
	double complex AbbSum458, AbbSum378, AbbSum239, AbbSum498
	double complex AbbSum177, AbbSum64, AbbSum305, AbbSum66
	double complex AbbSum32, AbbSum304, AbbSum86, AbbSum22
	double complex AbbSum72, AbbSum95, AbbSum132, AbbSum88
	double complex AbbSum112, AbbSum138, AbbSum18, AbbSum421
	double complex AbbSum449, AbbSum62, AbbSum60, AbbSum405
	double complex AbbSum92, AbbSum142, AbbSum286, AbbSum78
	double complex AbbSum130, AbbSum430, AbbSum156, AbbSum134
	double complex AbbSum151, AbbSum446, AbbSum136, AbbSum438
	double complex AbbSum128, AbbSum185, AbbSum556, AbbSum163
	double complex AbbSum164, AbbSum115, AbbSum550, AbbSum55
	double complex AbbSum123, AbbSum347, AbbSum125, AbbSum361
	double complex AbbSum152, AbbSum104, AbbSum467, AbbSum468
	double complex AbbSum422, AbbSum439, AbbSum28, AbbSum36
	double complex AbbSum82, AbbSum106, AbbSum117, AbbSum184
	double complex AbbSum47, AbbSum525, AbbSum165, AbbSum232
	double complex AbbSum288, AbbSum290, AbbSum186, AbbSum74
	double complex AbbSum353, AbbSum309, AbbSum460, AbbSum113
	double complex AbbSum1, AbbSum392, AbbSum393, AbbSum182
	double complex AbbSum308, AbbSum459, AbbSum143, AbbSum440
	double complex AbbSum423, AbbSum541, AbbSum485, AbbSum328
	double complex AbbSum515, AbbSum270, AbbSum407, AbbSum246
	double complex AbbSum46, AbbSum50, AbbSum247, AbbSum215
	double complex AbbSum68, AbbSum43, AbbSum102, AbbSum216
	double complex AbbSum271, AbbSum228, AbbSum367, AbbSum408
	double complex AbbSum38, AbbSum30, AbbSum48, AbbSum362
	double complex AbbSum29, AbbSum37, AbbSum263, AbbSum87
	double complex AbbSum306, AbbSum65, AbbSum155, AbbSum501
	double complex AbbSum389, AbbSum89, AbbSum101, AbbSum327
	double complex AbbSum558, AbbSum557, AbbSum199, AbbSum94
	double complex AbbSum229, AbbSum231, AbbSum242, AbbSum127
	double complex AbbSum84, AbbSum257, AbbSum107, AbbSum31
	double complex AbbSum166, AbbSum299, AbbSum207, AbbSum300
	double complex AbbSum187, AbbSum53, AbbSum258, AbbSum466
	double complex AbbSum158, AbbSum162, AbbSum417, AbbSum111
	double complex AbbSum126, AbbSum61, AbbSum80, AbbSum118
	double complex AbbSum289, AbbSum51, AbbSum40, AbbSum204
	double complex AbbSum536, AbbSum526, AbbSum120, AbbSum109
	double complex AbbSum548, AbbSum537, AbbSum486, AbbSum441
	double complex AbbSum243, AbbSum382, AbbSum368, AbbSum400
	double complex AbbSum91, AbbSum384, AbbSum85, AbbSum414
	double complex AbbSum131, AbbSum499, AbbSum79, AbbSum93
	double complex AbbSum122, AbbSum75, AbbSum280, AbbSum301
	double complex AbbSum360, AbbSum552, AbbSum551, AbbSum260
	double complex AbbSum57, AbbSum41, AbbSum59, AbbSum452
	double complex AbbSum133, AbbSum10, AbbSum9, AbbSum23
	double complex AbbSum110, AbbSum129, AbbSum349, AbbSum63
	double complex AbbSum206, AbbSum278, AbbSum476, AbbSum477
	double complex AbbSum469, AbbSum487, AbbSum329, AbbSum412
	double complex AbbSum431, AbbSum315, AbbSum330, AbbSum464
	double complex AbbSum448, AbbSum154, AbbSum116, AbbSum225
	double complex AbbSum105, AbbSum397, AbbSum398, AbbSum183
	double complex AbbSum45, AbbSum226, AbbSum56, AbbSum465
	double complex AbbSum159, AbbSum161, AbbSum317, AbbSum450
	double complex AbbSum137, AbbSum98, AbbSum153, AbbSum325
	double complex AbbSum108, AbbSum42, AbbSum391, AbbSum33
	double complex AbbSum193, AbbSum124, AbbSum73, AbbSum157
	double complex AbbSum58, AbbSum2, AbbSum90, AbbSum103
	double complex AbbSum67, AbbSum406, AbbSum314, AbbSum432
	double complex AbbSum522, AbbSum516, AbbSum447, AbbSum463
	double complex AbbSum279, AbbSum268, AbbSum146, AbbSum135
	double complex AbbSum433, AbbSum413, AbbSum424, AbbSum209
	double complex AbbSum348, AbbSum24, AbbSum12, AbbSum13
	double complex AbbSum11, AbbSum509, AbbSum357, AbbSum533
	double complex AbbSum560, AbbSum265, AbbSum511, AbbSum25
	double complex AbbSum429, AbbSum426, AbbSum456, AbbSum445
	double complex AbbSum568, AbbSum418, AbbSum514, AbbSum532
	double complex AbbSum480, AbbSum506, AbbSum482, AbbSum139
	double complex AbbSum192, AbbSum528, AbbSum554, AbbSum455
	double complex AbbSum237, AbbSum457, AbbSum481, AbbSum178
	double complex AbbSum437, AbbSum513, AbbSum180, AbbSum303
	double complex AbbSum527, AbbSum16, AbbSum17, AbbSum15
	double complex AbbSum179, AbbSum190, AbbSum44, AbbSum354
	double complex AbbSum387, AbbSum234, AbbSum238, AbbSum419
	double complex AbbSum114, AbbSum140, AbbSum531, AbbSum559
	double complex AbbSum355, AbbSum483, AbbSum316, AbbSum386
	double complex AbbSum553, AbbSum376, AbbSum565, AbbSum321
	double complex AbbSum189, AbbSum451, AbbSum403, AbbSum529
	double complex AbbSum530, AbbSum81, AbbSum388, AbbSum302
	double complex AbbSum100, AbbSum326, AbbSum167, AbbSum322
	double complex AbbSum435, AbbSum323, AbbSum436, AbbSum500
	double complex AbbSum320, AbbSum564, AbbSum291, AbbSum319
	double complex AbbSum262, AbbSum261, AbbSum233, AbbSum567
	double complex AbbSum83, AbbSum444, AbbSum434, AbbSum470
	double complex AbbSum396, AbbSum428, AbbSum292, AbbSum471
	double complex AbbSum284, AbbSum395, AbbSum370, AbbSum402
	double complex AbbSum283, AbbSum425, AbbSum150, AbbSum148
	double complex AbbSum191, AbbSum502, AbbSum404, AbbSum555
	double complex AbbSum149, AbbSum147, AbbSum282, AbbSum503
	double complex AbbSum141, AbbSum307, AbbSum267, AbbSum208
	double complex AbbSum420, AbbSum266, AbbSum561, AbbSum534
	double complex AbbSum535, AbbSum442, AbbSum54, AbbSum390
	double complex AbbSum356, AbbSum454, AbbSum563, AbbSum331
	double complex AbbSum350, AbbSum351, AbbSum210, AbbSum49
	double complex AbbSum538, AbbSum504, AbbSum510, AbbSum505
	double complex AbbSum363, AbbSum364, AbbSum366, AbbSum244
	double complex AbbSum345, AbbSum352, AbbSum269, AbbSum245
	double complex AbbSum359, AbbSum214, AbbSum211, AbbSum52
	double complex AbbSum39, AbbSum281, AbbSum318, AbbSum547
	double complex AbbSum287, AbbSum507, AbbSum383, AbbSum512
	double complex AbbSum508, AbbSum230, AbbSum379, AbbSum194
	double complex AbbSum380, AbbSum365, AbbSum381, AbbSum8
	double complex AbbSum34, AbbSum212, AbbSum213
	common /abbrev/ F1, F2, F3, F4, F12, F5, F6, F14, F7, F15, F9
	common /abbrev/ F8, F16, F10, F11, F13, Pair1, Pair2, Pair5
	common /abbrev/ Pair3, Pair4, Eps1, Eps2, Abb11, Abb12, Abb28
	common /abbrev/ Abb29, Abb35, Abb36, Abb13, Abb14, Abb1, Abb27
	common /abbrev/ Abb34, Abb2, Abb23, Abb15, Abb19, Abb24, Abb41
	common /abbrev/ Abb42, Abb16, Abb20, Abb3, Abb30, Abb37, Abb7
	common /abbrev/ Abb43, Abb44, Abb17, Abb21, Abb4, Abb31, Abb38
	common /abbrev/ Abb8, Abb25, Abb18, Abb22, Abb26, Abb5, Abb32
	common /abbrev/ Abb39, Abb9, Abb6, Abb33, Abb40, Abb10
	common /abbrev/ AbbSum3, AbbSum4, AbbSum394, AbbSum252
	common /abbrev/ AbbSum196, AbbSum19, AbbSum343, AbbSum26
	common /abbrev/ AbbSum6, AbbSum220, AbbSum119, AbbSum332
	common /abbrev/ AbbSum539, AbbSum5, AbbSum540, AbbSum333
	common /abbrev/ AbbSum168, AbbSum188, AbbSum70, AbbSum201
	common /abbrev/ AbbSum562, AbbSum171, AbbSum298, AbbSum181
	common /abbrev/ AbbSum497, AbbSum399, AbbSum427, AbbSum517
	common /abbrev/ AbbSum259, AbbSum495, AbbSum219, AbbSum415
	common /abbrev/ AbbSum221, AbbSum251, AbbSum313, AbbSum253
	common /abbrev/ AbbSum203, AbbSum202, AbbSum296, AbbSum336
	common /abbrev/ AbbSum410, AbbSum144, AbbSum338, AbbSum277
	common /abbrev/ AbbSum205, AbbSum409, AbbSum276, AbbSum493
	common /abbrev/ AbbSum195, AbbSum416, AbbSum76, AbbSum254
	common /abbrev/ AbbSum492, AbbSum69, AbbSum170, AbbSum461
	common /abbrev/ AbbSum198, AbbSum293, AbbSum275, AbbSum71
	common /abbrev/ AbbSum21, AbbSum235, AbbSum145, AbbSum462
	common /abbrev/ AbbSum224, AbbSum264, AbbSum339, AbbSum411
	common /abbrev/ AbbSum297, AbbSum274, AbbSum272, AbbSum385
	common /abbrev/ AbbSum197, AbbSum273, AbbSum312, AbbSum294
	common /abbrev/ AbbSum496, AbbSum223, AbbSum344, AbbSum160
	common /abbrev/ AbbSum35, AbbSum346, AbbSum97, AbbSum99
	common /abbrev/ AbbSum121, AbbSum310, AbbSum27, AbbSum96
	common /abbrev/ AbbSum255, AbbSum494, AbbSum488, AbbSum222
	common /abbrev/ AbbSum337, AbbSum169, AbbSum7, AbbSum173
	common /abbrev/ AbbSum341, AbbSum200, AbbSum295, AbbSum77
	common /abbrev/ AbbSum311, AbbSum256, AbbSum453, AbbSum340
	common /abbrev/ AbbSum227, AbbSum518, AbbSum521, AbbSum546
	common /abbrev/ AbbSum20, AbbSum14, AbbSum479, AbbSum545
	common /abbrev/ AbbSum478, AbbSum342, AbbSum489, AbbSum371
	common /abbrev/ AbbSum523, AbbSum473, AbbSum472, AbbSum240
	common /abbrev/ AbbSum248, AbbSum249, AbbSum372, AbbSum236
	common /abbrev/ AbbSum369, AbbSum544, AbbSum218, AbbSum519
	common /abbrev/ AbbSum175, AbbSum520, AbbSum566, AbbSum401
	common /abbrev/ AbbSum524, AbbSum543, AbbSum217, AbbSum474
	common /abbrev/ AbbSum374, AbbSum373, AbbSum490, AbbSum542
	common /abbrev/ AbbSum491, AbbSum334, AbbSum375, AbbSum285
	common /abbrev/ AbbSum250, AbbSum335, AbbSum324, AbbSum443
	common /abbrev/ AbbSum377, AbbSum549, AbbSum484, AbbSum358
	common /abbrev/ AbbSum475, AbbSum241, AbbSum172, AbbSum176
	common /abbrev/ AbbSum174, AbbSum458, AbbSum378, AbbSum239
	common /abbrev/ AbbSum498, AbbSum177, AbbSum64, AbbSum305
	common /abbrev/ AbbSum66, AbbSum32, AbbSum304, AbbSum86
	common /abbrev/ AbbSum22, AbbSum72, AbbSum95, AbbSum132
	common /abbrev/ AbbSum88, AbbSum112, AbbSum138, AbbSum18
	common /abbrev/ AbbSum421, AbbSum449, AbbSum62, AbbSum60
	common /abbrev/ AbbSum405, AbbSum92, AbbSum142, AbbSum286
	common /abbrev/ AbbSum78, AbbSum130, AbbSum430, AbbSum156
	common /abbrev/ AbbSum134, AbbSum151, AbbSum446, AbbSum136
	common /abbrev/ AbbSum438, AbbSum128, AbbSum185, AbbSum556
	common /abbrev/ AbbSum163, AbbSum164, AbbSum115, AbbSum550
	common /abbrev/ AbbSum55, AbbSum123, AbbSum347, AbbSum125
	common /abbrev/ AbbSum361, AbbSum152, AbbSum104, AbbSum467
	common /abbrev/ AbbSum468, AbbSum422, AbbSum439, AbbSum28
	common /abbrev/ AbbSum36, AbbSum82, AbbSum106, AbbSum117
	common /abbrev/ AbbSum184, AbbSum47, AbbSum525, AbbSum165
	common /abbrev/ AbbSum232, AbbSum288, AbbSum290, AbbSum186
	common /abbrev/ AbbSum74, AbbSum353, AbbSum309, AbbSum460
	common /abbrev/ AbbSum113, AbbSum1, AbbSum392, AbbSum393
	common /abbrev/ AbbSum182, AbbSum308, AbbSum459, AbbSum143
	common /abbrev/ AbbSum440, AbbSum423, AbbSum541, AbbSum485
	common /abbrev/ AbbSum328, AbbSum515, AbbSum270, AbbSum407
	common /abbrev/ AbbSum246, AbbSum46, AbbSum50, AbbSum247
	common /abbrev/ AbbSum215, AbbSum68, AbbSum43, AbbSum102
	common /abbrev/ AbbSum216, AbbSum271, AbbSum228, AbbSum367
	common /abbrev/ AbbSum408, AbbSum38, AbbSum30, AbbSum48
	common /abbrev/ AbbSum362, AbbSum29, AbbSum37, AbbSum263
	common /abbrev/ AbbSum87, AbbSum306, AbbSum65, AbbSum155
	common /abbrev/ AbbSum501, AbbSum389, AbbSum89, AbbSum101
	common /abbrev/ AbbSum327, AbbSum558, AbbSum557, AbbSum199
	common /abbrev/ AbbSum94, AbbSum229, AbbSum231, AbbSum242
	common /abbrev/ AbbSum127, AbbSum84, AbbSum257, AbbSum107
	common /abbrev/ AbbSum31, AbbSum166, AbbSum299, AbbSum207
	common /abbrev/ AbbSum300, AbbSum187, AbbSum53, AbbSum258
	common /abbrev/ AbbSum466, AbbSum158, AbbSum162, AbbSum417
	common /abbrev/ AbbSum111, AbbSum126, AbbSum61, AbbSum80
	common /abbrev/ AbbSum118, AbbSum289, AbbSum51, AbbSum40
	common /abbrev/ AbbSum204, AbbSum536, AbbSum526, AbbSum120
	common /abbrev/ AbbSum109, AbbSum548, AbbSum537, AbbSum486
	common /abbrev/ AbbSum441, AbbSum243, AbbSum382, AbbSum368
	common /abbrev/ AbbSum400, AbbSum91, AbbSum384, AbbSum85
	common /abbrev/ AbbSum414, AbbSum131, AbbSum499, AbbSum79
	common /abbrev/ AbbSum93, AbbSum122, AbbSum75, AbbSum280
	common /abbrev/ AbbSum301, AbbSum360, AbbSum552, AbbSum551
	common /abbrev/ AbbSum260, AbbSum57, AbbSum41, AbbSum59
	common /abbrev/ AbbSum452, AbbSum133, AbbSum10, AbbSum9
	common /abbrev/ AbbSum23, AbbSum110, AbbSum129, AbbSum349
	common /abbrev/ AbbSum63, AbbSum206, AbbSum278, AbbSum476
	common /abbrev/ AbbSum477, AbbSum469, AbbSum487, AbbSum329
	common /abbrev/ AbbSum412, AbbSum431, AbbSum315, AbbSum330
	common /abbrev/ AbbSum464, AbbSum448, AbbSum154, AbbSum116
	common /abbrev/ AbbSum225, AbbSum105, AbbSum397, AbbSum398
	common /abbrev/ AbbSum183, AbbSum45, AbbSum226, AbbSum56
	common /abbrev/ AbbSum465, AbbSum159, AbbSum161, AbbSum317
	common /abbrev/ AbbSum450, AbbSum137, AbbSum98, AbbSum153
	common /abbrev/ AbbSum325, AbbSum108, AbbSum42, AbbSum391
	common /abbrev/ AbbSum33, AbbSum193, AbbSum124, AbbSum73
	common /abbrev/ AbbSum157, AbbSum58, AbbSum2, AbbSum90
	common /abbrev/ AbbSum103, AbbSum67, AbbSum406, AbbSum314
	common /abbrev/ AbbSum432, AbbSum522, AbbSum516, AbbSum447
	common /abbrev/ AbbSum463, AbbSum279, AbbSum268, AbbSum146
	common /abbrev/ AbbSum135, AbbSum433, AbbSum413, AbbSum424
	common /abbrev/ AbbSum209, AbbSum348, AbbSum24, AbbSum12
	common /abbrev/ AbbSum13, AbbSum11, AbbSum509, AbbSum357
	common /abbrev/ AbbSum533, AbbSum560, AbbSum265, AbbSum511
	common /abbrev/ AbbSum25, AbbSum429, AbbSum426, AbbSum456
	common /abbrev/ AbbSum445, AbbSum568, AbbSum418, AbbSum514
	common /abbrev/ AbbSum532, AbbSum480, AbbSum506, AbbSum482
	common /abbrev/ AbbSum139, AbbSum192, AbbSum528, AbbSum554
	common /abbrev/ AbbSum455, AbbSum237, AbbSum457, AbbSum481
	common /abbrev/ AbbSum178, AbbSum437, AbbSum513, AbbSum180
	common /abbrev/ AbbSum303, AbbSum527, AbbSum16, AbbSum17
	common /abbrev/ AbbSum15, AbbSum179, AbbSum190, AbbSum44
	common /abbrev/ AbbSum354, AbbSum387, AbbSum234, AbbSum238
	common /abbrev/ AbbSum419, AbbSum114, AbbSum140, AbbSum531
	common /abbrev/ AbbSum559, AbbSum355, AbbSum483, AbbSum316
	common /abbrev/ AbbSum386, AbbSum553, AbbSum376, AbbSum565
	common /abbrev/ AbbSum321, AbbSum189, AbbSum451, AbbSum403
	common /abbrev/ AbbSum529, AbbSum530, AbbSum81, AbbSum388
	common /abbrev/ AbbSum302, AbbSum100, AbbSum326, AbbSum167
	common /abbrev/ AbbSum322, AbbSum435, AbbSum323, AbbSum436
	common /abbrev/ AbbSum500, AbbSum320, AbbSum564, AbbSum291
	common /abbrev/ AbbSum319, AbbSum262, AbbSum261, AbbSum233
	common /abbrev/ AbbSum567, AbbSum83, AbbSum444, AbbSum434
	common /abbrev/ AbbSum470, AbbSum396, AbbSum428, AbbSum292
	common /abbrev/ AbbSum471, AbbSum284, AbbSum395, AbbSum370
	common /abbrev/ AbbSum402, AbbSum283, AbbSum425, AbbSum150
	common /abbrev/ AbbSum148, AbbSum191, AbbSum502, AbbSum404
	common /abbrev/ AbbSum555, AbbSum149, AbbSum147, AbbSum282
	common /abbrev/ AbbSum503, AbbSum141, AbbSum307, AbbSum267
	common /abbrev/ AbbSum208, AbbSum420, AbbSum266, AbbSum561
	common /abbrev/ AbbSum534, AbbSum535, AbbSum442, AbbSum54
	common /abbrev/ AbbSum390, AbbSum356, AbbSum454, AbbSum563
	common /abbrev/ AbbSum331, AbbSum350, AbbSum351, AbbSum210
	common /abbrev/ AbbSum49, AbbSum538, AbbSum504, AbbSum510
	common /abbrev/ AbbSum505, AbbSum363, AbbSum364, AbbSum366
	common /abbrev/ AbbSum244, AbbSum345, AbbSum352, AbbSum269
	common /abbrev/ AbbSum245, AbbSum359, AbbSum214, AbbSum211
	common /abbrev/ AbbSum52, AbbSum39, AbbSum281, AbbSum318
	common /abbrev/ AbbSum547, AbbSum287, AbbSum507, AbbSum383
	common /abbrev/ AbbSum512, AbbSum508, AbbSum230, AbbSum379
	common /abbrev/ AbbSum194, AbbSum380, AbbSum365, AbbSum381
	common /abbrev/ AbbSum8, AbbSum34, AbbSum212, AbbSum213

	double complex cint1, cint2, cint3, cint4(3), cint5(3)
	double complex cint6(3), cint7(3), cint8(3), cint9(3)
	double complex cint10(3), cint11(3), cint12(3), cint13(2)
	double complex cint14, cint15, cint16, cint17(3), cint18(3)
	double complex cint19(3), cint20(3), cint21(3), cint22(3)
	double complex cint23(3), cint24(3), cint25(3), cint26(2)
	common /loopint/ cint1, cint2, cint3, cint4, cint5, cint6
	common /loopint/ cint7, cint8, cint9, cint10, cint11, cint12
	common /loopint/ cint13, cint14, cint15, cint16, cint17
	common /loopint/ cint18, cint19, cint20, cint21, cint22
	common /loopint/ cint23, cint24, cint25, cint26

	integer*8 iint1, iint2, iint3(3), iint4(3), iint5(3), iint6(3)
	integer*8 iint7(3), iint8(3), iint9(3), iint10(3), iint11(3)
	integer*8 iint12(2), iint13(3), iint14(3), iint15(3), iint16(3)
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

	integer cha1, hia1, his1, lpd1, qud1, quu1, sld1, sle1, sqd1
	integer sqe1, squ1, sqv1
	common /indices/ cha1, hia1, his1, lpd1, qud1, quu1, sld1
	common /indices/ sle1, sqd1, sqe1, squ1, sqv1

	double complex Cloop(1)
	common /coeff/ Cloop
