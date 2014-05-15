*     LanHEP output produced at Mon Apr 13 18:09:41 2009
*     Model named 'MSSMrc cha1 nans'

      double precision Sqrt2, pi, degree, hbar_c2,bogus
      parameter (Sqrt2=1.41421356237309504880168872421D0)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      parameter (degree = pi/180D0)
      parameter (hbar_c2 = 3.8937966D8)
      parameter (bogus = -1D123)
      double complex cI
      parameter (cI = (0D0, 1D0))

      double precision Divergence
      common /renorm/ Divergence

      double precision EE, SW, s12, s23, s13, GG, CW, S2W, C2W
      double precision MZ, MW, c12, c23, c13, Me, Mm, Mu, Md, Ml
      double precision Mc, Ms, Mt, Mb, wtop, wZ, wW, tb, sb, cb
      double precision t2b, c2b, s2b, mue, MG2, MG3, MG1, Zcsx
      double precision Zctx, MC01, MC02, Zcc2u, Zcc2v, Zcsigu
      double precision Zcsigv, Zccu, Zcsu, Zccv, Zcsv, ZcdetM
      double precision Zcsig1, Zcsig2, Zcsig3, Zcsig4, MC1, MC2
      double precision Zm11, Zm12, Zm21, Zm22, Zp11, Zp12, Zp21
      double precision Zp22, neutk, MNE1, MNE2, MNE3, MNE4, wC1
      double precision wC2, wNE2, wNE3, wNE4, MSG, wSG, Ml1, Ml2
      double precision Ml3, Mr1, Mr2, Mr3, Al, ls3, Am, ls2, Ae
      double precision ls1, MSne, MSnmu, MSntau, MSeLL, MSeRR
      double precision MSeLR, MSe1, MSe2, MSeth, ceth, seth, MSmuLL
      double precision MSmuLR, MSmuRR, MSmu1, MSmu2, MSmuth, smuth
      double precision cmuth, MSlLL, MSlLR, MSlRR, MStau1, MStau2
      double precision MSlth, wSne, wSnmu, wSntau, wSe1, wSe2
      double precision wSmu1, wSmu2, wStau1, wStau2, clth, slth
      double precision Mq1, Mq2, Mq3, Mu1, Mu2, Mu3, Md1, Md2
      double precision Md3, At, Ab, ds3, us3, As, ds2, Ac, us2
      double precision Au, us1, Ad, ds1, MSuLL, MSuLR, MSuRR, MSu1
      double precision MSu2, MSuth, cuth, suth, MScLL, MScLR, MScRR
      double precision MSc1, MSc2, MScth, ccth, scth, MStLL, MStLR
      double precision MStRR, MStop1, MStop2, MStth, ctth, stth
      double precision MSdLL, MSdLR, MSdRR, MSd1, MSd2, MSdth
      double precision cdth, sdth, MSsLL, MSsLR, MSsRR, MSs1, MSs2
      double precision MSsth, csth, ssth, MSbLL, MSbLR, MSbRR
      double precision MSbot1, MSbot2, MSbth, cbth, sbth, MH3
      double precision vevv, t2a1, t2a2, hha2, ca, sa, sapb, samb
      double precision capb, camb, s2a, c2a, MHH, Mh, MHc, wh
      double precision wHh, wH3, wHc, wSu1, wSu2, wSc1, wSc2, wStop1
      double precision wStop2, wSd1, wSd2, wSs1, wSs2, wSbot1
      double precision wSbot2, MSnl, MSl1, MSl2, MSnm, MSm1, MSm2
      double precision smth, cmth, MSt1, MSt2, MSb1, MSb2, nla
      double precision nlb, nld1, nld2, nlk, nlr, nle1, nle2, SW2
      double precision CW2, MW2, MZ2, EE2, GG2, Mesq, Mmsq, Mlsq
      double precision Musq, Mdsq, Mcsq, Mssq, Mbsq, Mtsq, MC1sq
      double precision MC2sq, MNE1sq, MNE2sq, MNE3sq, MNE4sq, MSe1sq
      double precision MSe2sq, MSmu1sq, MSmu2sq, MStau1sq, MStau2sq
      double precision MSu1sq, MSu2sq, MSd1sq, MSd2sq, MSc1sq
      double precision MSc2sq, MSs1sq, MSs2sq, MStop1sq, MStop2sq
      double precision MSbot1sq, MSbot2sq, MSnesq, MSnmusq, MSntausq
      double precision MHHsq, MHH2, Mhsq, MH3sq, MHcsq, Alfa, Alfa2

      double complex Zn11, Zn12, Zn13, Zn14, Zn21, Zn22, Zn23
      double complex Zn24, Zn31, Zn32, Zn33, Zn34, Zn41, Zn42
      double complex Zn43, Zn44

      double precision AAABR(8514)
      double complex AAABC(1076)

      double precision quuMass, qudMass, lpdMass, neuMass, chaMass
      double precision sluMass, sldMass, sleMass, squMass, sqvMass
      double precision sqdMass, sqeMass, hisMass, MTR001, MTR002, hisWidt
      double precision MTR003, MTR004, MTR005, MTR006, MTR007
      double precision MTR008, MTR009, MTR010, MTR011, MTR012
      double precision MTR013, MTR014, MTR015, MTR016, MTR017
      double precision MTR018, MTR019, MTR020, MTR021, MTR022
      double precision MTR023, MTR024, MTR025, MTR026, MTR027
      double precision MTR028, MTR029, MTR030, MTR031, MTR032
      double precision MTR033, MTR034, MTR035, MTR036, MTR037
      double precision MTR038, MTR039, MTR040, MTR041, MTR042
      double precision MTR043, MTR044, MTR045, MTR046, MTR047
      double precision MTR048, MTR049, MTR050, MTR051, MTR052
      double precision MTR053, MTR054, MTR055, MTR056, MTR057
      double precision MTR058, MTR059, MTR060, MTR061, MTR062
      double precision MTR063, MTR064, MTR065, MTR066, MTR067
      double precision MTR068, MTR073, MTR074, MTR075, MTR076
      double precision MTR077, MTR078, MTR079, MTR080, MTR081
      double precision MTR086, MTR087, MTR093, MTR094, MTR095
      double precision MTR096, MTR097, MTR098, MTR099, MTR100
      double precision MTR101, MTR102, MTR103, MTR104, MTR105
      double precision MTR110, MTR111, MTR112, MTR113, MTR114
      double precision MTR138, MTR139, MTR146, MTR147, MTR148
      double precision MTR149, MTR150, MTR151, MTR152, MTR153
      double precision MTR154, MTR155, MTR156, MTR157, MTR158
      double precision MTR159, MTR160, MTR161, MTR162, MTR163
      double precision MTR164, MTR165, MTR166, MTR167, MTR168
      double precision MTR169, MTR170, MTR171, MTR172, MTR173
      double precision MTR174, MTR175, MTR176, MTR177, MTR178
      double precision MTR179, MTR180, MTR181, MTR182, MTR183
      double precision MTR184, MTR185, MTR186, MTR187, MTR188
      double precision MTR189, MTR190, MTR191, MTR192, MTR193
      double precision MTR194, MTR195, MTR196, MTR197, MTR198
      double precision MTR199, MTR200, MTR201, MTR202, MTR203
      double precision MTR204, MTR205, MTR206, MTR207, MTR208
      double precision MTR209, MTR210, MTR211, MTR212, MTR213
      double precision MTR214, MTR215, MTR216, MTR217, MTR218
      double precision MTR219, MTR220, MTR221, MTR222, MTR223
      double precision MTR224, MTR225, MTR226, MTR227, MTR228
      double precision MTR229, MTR230, MTR231, MTR232, MTR233
      double precision MTR234, MTR235, MTR236, MTR237, MTR238
      double precision MTR239, MTR240, MTR241, MTR242, MTR243
      double precision MTR244, MTR245, MTR246, MTR247, MTR248
      double precision MTR249, MTR250, MTR251, MTR252, MTR253
      double precision MTR254, MTR255, MTR256, MTR257, MTR258
      double precision MTR259, MTR260, MTR261, MTR262, MTR263
      double precision MTR264, MTR265, MTR266, MTR267, MTR268
      double precision MTR269, MTR270, MTR271, MTR272, MTR273
      double precision MTR274, MTR275, MTR276, MTR277, MTR278
      double precision MTR279, MTR280, MTR281, MTR282, MTR283
      double precision MTR284, MTR285, MTR286, MTR287, MTR288
      double precision MTR289, MTR290, MTR291, MTR292, MTR293
      double precision MTR294, MTR295, MTR296, MTR297, MTR298
      double precision MTR299, MTR300, MTR301, MTR302, MTR303
      double precision MTR304, MTR305, MTR306, MTR307, MTR308
      double precision MTR309, MTR310, MTR311, MTR312, MTR313
      double precision MTR314, MTR315, MTR316, MTR317, MTR318
      double precision MTR319, MTR320, MTR321, MTR322, MTR323
      double precision MTR324, MTR325, MTR326, MTR327, MTR328
      double precision MTR329, MTR330, MTR331, MTR332, MTR333
      double precision MTR334, MTR335, MTR336, MTR337, MTR338
      double precision MTR339, MTR340, MTR341, MTR342, MTR343
      double precision MTR344, MTR345, MTR346, MTR347, MTR348
      double precision MTR349, MTR350, MTR351, MTR352, MTR353
      double precision MTR354, MTR355, MTR356, MTR357, MTR358
      double precision MTR359, MTR360, MTR361, MTR362, MTR363
      double precision MTR364, MTR365, MTR366, MTR367, MTR368
      double precision MTR369, MTR370, MTR371, MTR372, MTR373
      double precision MTR374, MTR375, MTR376, MTR377, MTR378
      double precision MTR379, MTR380, MTR381, MTR382, MTR383
      double precision MTR384, MTR385, MTR386, MTR387, MTR388
      double precision MTR389, MTR390, MTR391, MTR392, MTR393
      double precision MTR394, MTR395, MTR396, MTR397, MTR398
      double precision MTR399, MTR400, MTR401, MTR402, MTR403
      double precision MTR404, MTR405, MTR406, MTR407, MTR408
      double precision MTR409, MTR410, MTR411, MTR412, MTR413
      double precision MTR414, MTR415, MTR416, MTR417, MTR418
      double precision MTR419, MTR420, MTR421, MTR422, MTR423
      double precision MTR424, MTR425, MTR426, MTR427
      double complex MTR069, MTR070, MTR071, MTR072, MTR082, MTR083
      double complex MTR084, MTR085, MTR088, MTR089, MTR090, MTR091
      double complex MTR092, MTR106, MTR107, MTR108, MTR109, MTR115
      double complex MTR116, MTR117, MTR118, MTR119, MTR120, MTR121
      double complex MTR122, MTR123, MTR124, MTR125, MTR126, MTR127
      double complex MTR128, MTR129, MTR130, MTR131, MTR132, MTR133
      double complex MTR134, MTR135, MTR136, MTR137, MTR140, MTR141
      double complex MTR142, MTR143, MTR144, MTR145

      dimension
     & quuMass(3), qudMass(3), lpdMass(3), neuMass(4), chaMass(2), 
     & sluMass(3), sldMass(3), sleMass(3), squMass(3), sqvMass(3), 
     & sqdMass(3), sqeMass(3), hisMass(2), MTR001(2), MTR002(2), hisWidt(2),
     & MTR003(2), MTR004(3), MTR005(3), MTR006(3), MTR007(3), 
     & MTR008(3), MTR009(3), MTR010(2), MTR011(2), MTR012(3), 
     & MTR013(3), MTR014(3), MTR015(2), MTR016(3), MTR017(3), 
     & MTR018(2,3), MTR019(2,3), MTR020(3), MTR021(2,3), 
     & MTR022(2,3), MTR023(3), MTR024(3), MTR025(3), MTR026(2,3), 
     & MTR027(2,3), MTR028(3), MTR029(3), MTR030(2,3), MTR031(3), 
     & MTR032(2,3), MTR033(2,3), MTR034(2,3), MTR035(2), 
     & MTR036(2), MTR037(2,2,2), MTR038(2), MTR039(2), MTR040(3), 
     & MTR041(3), MTR042(3), MTR043(3), MTR044(3), MTR045(3), 
     & MTR046(3), MTR047(3), MTR048(3), MTR049(3), MTR050(3), 
     & MTR051(3), MTR052(3), MTR053(3), MTR054(3), MTR055(2), 
     & MTR056(2), MTR057(2), MTR058(2), MTR059(2), MTR060(2), 
     & MTR061(2,2), MTR062(2,2), MTR063(2,2), MTR064(2,2), 
     & MTR065(2,2,2), MTR066(2,2,2), MTR067(2,3), MTR068(2,3), 
     & MTR069(2,4), MTR070(2,4), MTR071(2,4), MTR072(2,4), 
     & MTR073(2,3), MTR074(2,3), MTR075(2,3), MTR076(2,3), 
     & MTR077(3), MTR078(3), MTR079(3,2), MTR080(3), MTR081(3), 
     & MTR082(4,3), MTR083(4,3), MTR084(4,3), MTR085(4,3), 
     & MTR086(2,3), MTR087(2,3), MTR088(4,3), MTR089(4,3), 
     & MTR090(4,3), MTR091(4,3), MTR092(4,3), MTR093(3), 
     & MTR094(3), MTR095(3,2), MTR096(3), MTR097(3), MTR098(3), 
     & MTR099(3), MTR100(2,3), MTR101(2,3), MTR102(2,3), 
     & MTR103(2,3), MTR104(3), MTR105(3), MTR106(4,3), MTR107(4,3), 
     & MTR108(4,3), MTR109(4,3), MTR110(3), MTR111(3), MTR112(3,2), 
     & MTR113(3), MTR114(3), MTR115(4,3), MTR116(4,3), MTR117(4,3), 
     & MTR118(4,3), MTR119(4,3), MTR120(4,3), MTR121(4,3), 
     & MTR122(4,3), MTR123(4,3), MTR124(4,3), MTR125(4,3), 
     & MTR126(4,3), MTR127(4,3), MTR128(2,4), MTR129(2,4), 
     & MTR130(2,4), MTR131(2,4), MTR132(4,4), MTR133(4,4), 
     & MTR134(4,4), MTR135(4,4), MTR136(4,4,2), MTR137(4,4,2), 
     & MTR138(2,2), MTR139(2,2), MTR140(2,4), MTR141(2,4), 
     & MTR142(2,4), MTR143(2,4), MTR144(4,4), MTR145(4,4), 
     & MTR146(2), MTR147(2), MTR148(2), MTR149(2,2), MTR150(2), 
     & MTR151(2), MTR152(2,2), MTR153(3), MTR154(3), MTR155(3), 
     & MTR156(3), MTR157(3), MTR158(3), MTR159(3), MTR160(3), 
     & MTR161(3), MTR162(3), MTR163(2,2), MTR164(3), MTR165(3), 
     & MTR166(3), MTR167(3), MTR168(3), MTR169(3), MTR170(2), 
     & MTR171(3), MTR172(3), MTR173(3), MTR174(3), MTR175(3), 
     & MTR176(3), MTR177(2,3), MTR178(2,3), MTR179(3), MTR180(3), 
     & MTR181(3), MTR182(3), MTR183(3), MTR184(3), MTR185(3), 
     & MTR186(2,3), MTR187(2,3), MTR188(3), MTR189(3), MTR190(3), 
     & MTR191(2,3), MTR192(2,3), MTR193(2), MTR194(2,2), 
     & MTR195(3), MTR196(3), MTR197(3), MTR198(3), MTR199(3), 
     & MTR200(3), MTR201(3), MTR202(3), MTR203(3), MTR204(2,2), 
     & MTR205(3), MTR206(3), MTR207(3), MTR208(3), MTR209(3), 
     & MTR210(3), MTR211(3), MTR212(3), MTR213(3), MTR214(3), 
     & MTR215(3), MTR216(3), MTR217(3), MTR218(2,2), MTR219(3,3), 
     & MTR220(3,3), MTR221(3,3), MTR222(3,3), MTR223(3,3), 
     & MTR224(3,3), MTR225(3,3), MTR226(3,3), MTR227(3,3), 
     & MTR228(3,3), MTR229(3,3), MTR230(3,3), MTR231(3,3), 
     & MTR232(3,3), MTR233(3,3), MTR234(3,3), MTR235(3,3), 
     & MTR236(3,3), MTR237(3,3), MTR238(3,3), MTR239(3,3), 
     & MTR240(3,3), MTR241(3,3), MTR242(3,3), MTR243(3,3), 
     & MTR244(3,3), MTR245(3,3), MTR246(3,3), MTR247(3,3), 
     & MTR248(3,3), MTR249(3,3), MTR250(3,3), MTR251(3), 
     & MTR252(3), MTR253(3), MTR254(2,3), MTR255(3), MTR256(3), 
     & MTR257(2,2,3), MTR258(2,2,3), MTR259(3,3), MTR260(3,3), 
     & MTR261(3,3), MTR262(3,3), MTR263(3,3), MTR264(3,3), 
     & MTR265(3,3), MTR266(3,3), MTR267(3,3), MTR268(3,3), 
     & MTR269(3,3), MTR270(3,3), MTR271(3,3), MTR272(3,3), 
     & MTR273(3,3), MTR274(3,3), MTR275(3), MTR276(3), MTR277(2,3), 
     & MTR278(3), MTR279(2,2,3), MTR280(3,3), MTR281(3,3), 
     & MTR282(3,3), MTR283(3,3), MTR284(3,3), MTR285(3,3), 
     & MTR286(3,3), MTR287(3,3), MTR288(3), MTR289(2,2,3), 
     & MTR290(3,3), MTR291(3,3), MTR292(3,3), MTR293(3,3), 
     & MTR294(3,3), MTR295(3,3), MTR296(3,3), MTR297(3,3), 
     & MTR298(3,3), MTR299(3,3), MTR300(3,3), MTR301(3,3), 
     & MTR302(3,3), MTR303(3,3), MTR304(3,3), MTR305(3,3), 
     & MTR306(3,3), MTR307(3,3), MTR308(3,3), MTR309(3,3), 
     & MTR310(3,3), MTR311(3,3), MTR312(3,3), MTR313(3,3), 
     & MTR314(3,3), MTR315(3,3), MTR316(3,3), MTR317(3,3), 
     & MTR318(3,3), MTR319(3,3), MTR320(3,3), MTR321(3), 
     & MTR322(3), MTR323(3), MTR324(3), MTR325(2,3), MTR326(2,3), 
     & MTR327(3), MTR328(3), MTR329(2,2,3), MTR330(2,2,3), 
     & MTR331(3,3), MTR332(3,3), MTR333(3,3), MTR334(3,3), 
     & MTR335(3,3), MTR336(3,3), MTR337(3,3), MTR338(3,3), 
     & MTR339(3,3), MTR340(3,3), MTR341(3,3), MTR342(3,3), 
     & MTR343(3), MTR344(3), MTR345(3), MTR346(2,3), MTR347(2,3), 
     & MTR348(3), MTR349(2,2,3), MTR350(3,3), MTR351(3,3), 
     & MTR352(3,3), MTR353(3,3), MTR354(3,3), MTR355(3,3), 
     & MTR356(3,3), MTR357(3,3), MTR358(3,3), MTR359(3,3), 
     & MTR360(3,3), MTR361(3,3), MTR362(3,3), MTR363(3), 
     & MTR364(3), MTR365(3), MTR366(2,2,3), MTR367(2,2,3), 
     & MTR368(3,3), MTR369(3), MTR370(3), MTR371(2,2,3), 
     & MTR372(2,2), MTR373(2,2), MTR374(2,2,2,2), MTR375(2), 
     & MTR376(2), MTR377(3), MTR378(3), MTR379(3), MTR380(3), 
     & MTR381(3), MTR382(3), MTR383(3), MTR384(3), MTR385(3), 
     & MTR386(3), MTR387(3), MTR388(3), MTR389(3), MTR390(3), 
     & MTR391(3), MTR392(3), MTR393(3), MTR394(3), MTR395(3), 
     & MTR396(3), MTR397(3), MTR398(3), MTR399(3), MTR400(3), 
     & MTR401(3), MTR402(3), MTR403(3), MTR404(3), MTR405(3), 
     & MTR406(3), MTR407(3), MTR408(3), MTR409(3), MTR410(3), 
     & MTR411(3), MTR412(3), MTR413(3), MTR414(3), MTR415(3), 
     & MTR416(3), MTR417(3), MTR418(3), MTR419(3), MTR420(3), 
     & MTR421(3), MTR422(3), MTR423(3), MTR424(3), MTR425(3), 
     & MTR426(2), MTR427(2)

      common /mdl_mtrces/
     & quuMass, qudMass, lpdMass, neuMass, chaMass, sluMass, 
     & sldMass, sleMass, squMass, sqvMass, sqdMass, sqeMass, 
     & hisMass, MTR001, MTR002, MTR003, MTR004, MTR005, hisWidt, 
     & MTR006, MTR007, MTR008, MTR009, MTR010, MTR011, MTR012, 
     & MTR013, MTR014, MTR015, MTR016, MTR017, MTR018, MTR019, 
     & MTR020, MTR021, MTR022, MTR023, MTR024, MTR025, MTR026, 
     & MTR027, MTR028, MTR029, MTR030, MTR031, MTR032, MTR033, 
     & MTR034, MTR035, MTR036, MTR037, MTR038, MTR039, MTR040, 
     & MTR041, MTR042, MTR043, MTR044, MTR045, MTR046, MTR047, 
     & MTR048, MTR049, MTR050, MTR051, MTR052, MTR053, MTR054, 
     & MTR055, MTR056, MTR057, MTR058, MTR059, MTR060, MTR061, 
     & MTR062, MTR063, MTR064, MTR065, MTR066, MTR067, MTR068, 
     & MTR069, MTR070, MTR071, MTR072, MTR073, MTR074, MTR075, 
     & MTR076, MTR077, MTR078, MTR079, MTR080, MTR081, MTR082, 
     & MTR083, MTR084, MTR085, MTR086, MTR087, MTR088, MTR089, 
     & MTR090, MTR091, MTR092, MTR093, MTR094, MTR095, MTR096, 
     & MTR097, MTR098, MTR099, MTR100, MTR101, MTR102, MTR103, 
     & MTR104, MTR105, MTR106, MTR107, MTR108, MTR109, MTR110, 
     & MTR111, MTR112, MTR113, MTR114, MTR115, MTR116, MTR117, 
     & MTR118, MTR119, MTR120, MTR121, MTR122, MTR123, MTR124, 
     & MTR125, MTR126, MTR127, MTR128, MTR129, MTR130, MTR131, 
     & MTR132, MTR133, MTR134, MTR135, MTR136, MTR137, MTR138, 
     & MTR139, MTR140, MTR141, MTR142, MTR143, MTR144, MTR145, 
     & MTR146, MTR147, MTR148, MTR149, MTR150, MTR151, MTR152, 
     & MTR153, MTR154, MTR155, MTR156, MTR157, MTR158, MTR159, 
     & MTR160, MTR161, MTR162, MTR163, MTR164, MTR165, MTR166, 
     & MTR167, MTR168, MTR169, MTR170, MTR171, MTR172, MTR173, 
     & MTR174, MTR175, MTR176, MTR177, MTR178, MTR179, MTR180, 
     & MTR181, MTR182, MTR183, MTR184, MTR185, MTR186, MTR187, 
     & MTR188, MTR189, MTR190, MTR191, MTR192, MTR193, MTR194, 
     & MTR195, MTR196, MTR197, MTR198, MTR199, MTR200, MTR201, 
     & MTR202, MTR203, MTR204, MTR205, MTR206, MTR207, MTR208, 
     & MTR209, MTR210, MTR211, MTR212, MTR213, MTR214, MTR215, 
     & MTR216, MTR217, MTR218, MTR219, MTR220, MTR221, MTR222, 
     & MTR223, MTR224, MTR225, MTR226, MTR227, MTR228, MTR229, 
     & MTR230, MTR231, MTR232, MTR233, MTR234, MTR235, MTR236, 
     & MTR237, MTR238, MTR239, MTR240, MTR241, MTR242, MTR243, 
     & MTR244, MTR245, MTR246, MTR247, MTR248, MTR249, MTR250, 
     & MTR251, MTR252, MTR253, MTR254, MTR255, MTR256, MTR257, 
     & MTR258, MTR259, MTR260, MTR261, MTR262, MTR263, MTR264, 
     & MTR265, MTR266, MTR267, MTR268, MTR269, MTR270, MTR271, 
     & MTR272, MTR273, MTR274, MTR275, MTR276, MTR277, MTR278, 
     & MTR279, MTR280, MTR281, MTR282, MTR283, MTR284, MTR285, 
     & MTR286, MTR287, MTR288, MTR289, MTR290, MTR291, MTR292, 
     & MTR293, MTR294, MTR295, MTR296, MTR297, MTR298, MTR299, 
     & MTR300, MTR301, MTR302, MTR303, MTR304, MTR305, MTR306, 
     & MTR307, MTR308, MTR309, MTR310, MTR311, MTR312, MTR313, 
     & MTR314, MTR315, MTR316, MTR317, MTR318, MTR319, MTR320, 
     & MTR321, MTR322, MTR323, MTR324, MTR325, MTR326, MTR327, 
     & MTR328, MTR329, MTR330, MTR331, MTR332, MTR333, MTR334, 
     & MTR335, MTR336, MTR337, MTR338, MTR339, MTR340, MTR341, 
     & MTR342, MTR343, MTR344, MTR345, MTR346, MTR347, MTR348, 
     & MTR349, MTR350, MTR351, MTR352, MTR353, MTR354, MTR355, 
     & MTR356, MTR357, MTR358, MTR359, MTR360, MTR361, MTR362, 
     & MTR363, MTR364, MTR365, MTR366, MTR367, MTR368, MTR369, 
     & MTR370, MTR371, MTR372, MTR373, MTR374, MTR375, MTR376, 
     & MTR377, MTR378, MTR379, MTR380, MTR381, MTR382, MTR383, 
     & MTR384, MTR385, MTR386, MTR387, MTR388, MTR389, MTR390, 
     & MTR391, MTR392, MTR393, MTR394, MTR395, MTR396, MTR397, 
     & MTR398, MTR399, MTR400, MTR401, MTR402, MTR403, MTR404, 
     & MTR405, MTR406, MTR407, MTR408, MTR409, MTR410, MTR411, 
     & MTR412, MTR413, MTR414, MTR415, MTR416, MTR417, MTR418, 
     & MTR419, MTR420, MTR421, MTR422, MTR423, MTR424, MTR425, 
     & MTR426, MTR427

      common /mdl_para/
     &    EE, SW, s12, s23, s13, GG, CW, S2W, C2W, MZ, MW, c12,
     &    c23, c13, Me, Mm, Mu, Md, Ml, Mc, Ms, Mt, Mb, wtop, wZ,
     &    wW, tb, sb, cb, t2b, c2b, s2b, mue, MG2, MG3, MG1, Zcsx,
     &    Zctx, MC01, MC02, Zcc2u, Zcc2v, Zcsigu, Zcsigv, Zccu,
     &    Zcsu, Zccv, Zcsv, ZcdetM, Zcsig1, Zcsig2, Zcsig3, Zcsig4,
     &    MC1, MC2, Zm11, Zm12, Zm21, Zm22, Zp11, Zp12, Zp21, Zp22,
     &    neutk, Zn11, Zn12, Zn13, Zn14, Zn21, Zn22, Zn23, Zn24,
     &    Zn31, Zn32, Zn33, Zn34, Zn41, Zn42, Zn43, Zn44, MNE1,
     &    MNE2, MNE3, MNE4, wC1, wC2, wNE2, wNE3, wNE4, MSG, wSG,
     &    Ml1, Ml2, Ml3, Mr1, Mr2, Mr3, Al, ls3, Am, ls2, Ae, ls1,
     &    MSne, MSnmu, MSntau, MSeLL, MSeRR, MSeLR, MSe1, MSe2,
     &    MSeth, ceth, seth, MSmuLL, MSmuLR, MSmuRR, MSmu1, MSmu2,
     &    MSmuth, smuth, cmuth, MSlLL, MSlLR, MSlRR, MStau1, MStau2,
     &    MSlth, wSne, wSnmu, wSntau, wSe1, wSe2, wSmu1, wSmu2,
     &    wStau1, wStau2, clth, slth, Mq1, Mq2, Mq3, Mu1, Mu2,
     &    Mu3, Md1, Md2, Md3, At, Ab, ds3, us3, As, ds2, Ac, us2,
     &    Au, us1, Ad, ds1, MSuLL, MSuLR, MSuRR, MSu1, MSu2, MSuth,
     &    cuth, suth, MScLL, MScLR, MScRR, MSc1, MSc2, MScth, ccth,
     &    scth, MStLL, MStLR, MStRR, MStop1, MStop2, MStth, ctth,
     &    stth, MSdLL, MSdLR, MSdRR, MSd1, MSd2, MSdth, cdth, sdth,
     &    MSsLL, MSsLR, MSsRR, MSs1, MSs2, MSsth, csth, ssth, MSbLL,
     &    MSbLR, MSbRR, MSbot1, MSbot2, MSbth, cbth, sbth, MH3,
     &    vevv, t2a1, t2a2, hha2, ca, sa, sapb, samb, capb, camb,
     &    s2a, c2a, MHH, Mh, MHc, wh, wHh, wH3, wHc, wSu1, wSu2,
     &    wSc1, wSc2, wStop1, wStop2, wSd1, wSd2, wSs1, wSs2, wSbot1,
     &    wSbot2, MSnl, MSl1, MSl2, MSnm, MSm1, MSm2, smth, cmth,
     &    MSt1, MSt2, MSb1, MSb2, nla, nlb, nld1, nld2, nlk, nlr,
     &    nle1, nle2, SW2, CW2, MW2, MZ2, EE2, GG2, Mesq, Mmsq,
     &    Mlsq, Musq, Mdsq, Mcsq, Mssq, Mbsq, Mtsq, MC1sq, MC2sq,
     &    MNE1sq, MNE2sq, MNE3sq, MNE4sq, MSe1sq, MSe2sq, MSmu1sq,
     &    MSmu2sq, MStau1sq, MStau2sq, MSu1sq, MSu2sq, MSd1sq,
     &    MSd2sq, MSc1sq, MSc2sq, MSs1sq, MSs2sq, MStop1sq, MStop2sq,
     &    MSbot1sq, MSbot2sq, MSnesq, MSnmusq, MSntausq, MHHsq,
     &    MHH2, Mhsq, MH3sq, MHcsq, Alfa, Alfa2, AAABR, AAABC

