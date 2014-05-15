*     LanHEP output produced at Wed Apr 10 22:10:25 2013
*     Model named 'MSSM&CPV'

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
      double precision MZ, MW, c12, c23, c13, Me, Mm, Mu, Md, Mc
      double precision Ms, wt, wZ, wW, tb, Ml, Mt, Mb, sb, cb
      double precision t2b, c2b, s2b, MG3, Ml1, Ml2, Ml3, Mr1
      double precision Mr2, Mr3, Ae, ls1, Am, ls2, MSne, MSnm
      double precision MSnl, MSeLL, MSeRR, MSeLR, MSe1, MSe2, MSeth
      double precision Zl11, Zl14, Zl41, Zl44, MSmuLL, MSmuLR
      double precision MSmuRR, MSm1, MSm2, MSmuth, Zl22, Zl25
      double precision Zl52, Zl55, Mq1, Mq2, Mq3, Mu1, Mu2, Mu3
      double precision Md1, Md2, Md3, Zual4, Zuas4, Zur0, ZuQMt
      double precision Zur1, Zulam5, Zur2, Zur3, Zulam6, Zur4
      double precision MtopS, Au, Ad, us1, ds1, As, ds2, Ac, us2
      double precision MSuLL, MSuLR, MSuRR, MSu1, MSu2, MSuth
      double precision Zu11, Zu14, Zu41, Zu44, MSdLL, MSdLR, MSdRR
      double precision MSd1, MSd2, MSdth, Zd11, Zd14, Zd41, Zd44
      double precision MScLL, MScLR, MScRR, MSc1, MSc2, MScth
      double precision Zu22, Zu25, Zu52, Zu55, MSsLL, MSsLR, MSsRR
      double precision MSs1, MSs2, MSsth, Zd22, Zd25, Zd52, Zd55
      double precision aMG1, fiMG1, aMG2, fiMG2, MG1r, MG1i, MG2r
      double precision MG2i, MG3i, amu, fimu, mur, mui, Alr, Ali
      double precision Atr, Ati, Abr, Abi, MHc, wHc, fMG1, fMG2
      double precision fmu, neutct, Zn11r, Zn11i, Zn12r, Zn12i
      double precision Zn13r, Zn13i, Zn14r, Zn14i, Zn21r, Zn21i
      double precision Zn22r, Zn22i, Zn23r, Zn23i, Zn24r, Zn24i
      double precision Zn31r, Zn31i, Zn32r, Zn32i, Zn33r, Zn33i
      double precision Zn34r, Zn34i, Zn41r, Zn41i, Zn42r, Zn42i
      double precision Zn43r, Zn43i, Zn44r, Zn44i, MNE1, wNE1
      double precision MNE2, wNE2, MNE3, wNE3, MNE4, wNE4, chatkn
      double precision Zm11r, Zp11r, Zm11i, Zp11i, Zm12r, Zp12r
      double precision Zm12i, Zp12i, Zm21r, Zp21r, Zm21i, Zp21i
      double precision Zm22r, Zp22r, Zm22i, Zp22i, MC1, wC1, MC2
      double precision wC2, MSG, wSG, hM11, hM12, hM22, hM33, htkn
      double precision Mh1, wh1, Mh2, wh2, Mh3, wh3, Zh11, Zh12
      double precision Zh13, Zh21, Zh22, Zh23, Zh31, Zh32, Zh33
      double precision MSbLL, MSbLRr, MSbLRi, MSbRR, MSb1, MSb2
      double precision MSbth, MSbfi, Zd33r, Zd36r, Zd63r, Zd66r
      double precision Zd33i, Zd36i, Zd63i, Zd66i, MStLL, MStLRr
      double precision MStLRi, MStRR, MSt1, MSt2, MStth, MStfi
      double precision Zu33r, Zu36r, Zu63r, Zu66r, Zu33i, Zu36i
      double precision Zu63i, Zu66i, MSlLL, MSlLRr, MSlLRi, MSlRR
      double precision MSl1, MSl2, MSlth, MSlfi, Zl33r, Zl36r
      double precision Zl63r, Zl66r, Zl33i, Zl36i, Zl63i, Zl66i
      double precision wSe1, wSe2, wSm1, wSm2, wSl1, wSl2, wSne
      double precision wSnm, wSnl, wSu1, wSu2, wSd1, wSd2, wSc1
      double precision wSc2, wSs1, wSs2, wSt1, wSt2, wSb1, wSb2
      double precision Mh1sq, Mh2sq, Mh3sq, Mesq, Mmsq, Mlsq, Musq
      double precision Mdsq, Mcsq, Mssq, Mbsq, Mtsq, MSe1sq, MSm1sq
      double precision MSl1sq, MSu1sq, MSd1sq, MSc1sq, MSs1sq
      double precision MSb1sq, MSt1sq, MSe2sq, MSm2sq, MSl2sq
      double precision MSu2sq, MSd2sq, MSc2sq, MSs2sq, MSb2sq
      double precision MSt2sq, MSnesq, MSnmsq, MSnlsq

      double precision AAABR(1274)

      double precision quuMass, qudMass, lpdMass, neuMass, chaMass
      double precision sluMass, sldMass, sleMass, squMass, sqvMass
      double precision sqdMass, sqeMass, hisMass, MTR001, MTR002
      double precision MTR003, MTR019, MTR021, MTR025, MTR027
      double precision MTR030, MTR033, MTR035, MTR040, MTR042
      double precision MTR045, MTR047, MTR052, MTR054, MTR055
      double precision MTR056, MTR057, MTR060, MTR064, MTR068
      double precision MTR073, MTR078, MTR083, MTR084, MTR085
      double precision MTR086, MTR087, MTR102, MTR105, MTR106
      double precision MTR120, MTR123, MTR124, MTR125, MTR126
      double precision MTR143, MTR197, MTR198, MTR199, MTR208
      double precision MTR209, MTR210, MTR217, MTR218, MTR219
      double precision MTR220, MTR235, MTR236, MTR237, MTR238
      double precision MTR251, MTR252, MTR253, MTR254, MTR269
      double precision MTR270, MTR271, MTR272, MTR273, MTR274, hisW
      double complex MTR004, MTR005, MTR006, MTR007, MTR008, MTR009
      double complex MTR010, MTR011, MTR012, MTR013, MTR014, MTR015
      double complex MTR016, MTR017, MTR018, MTR020, MTR022, MTR023
      double complex MTR024, MTR026, MTR028, MTR029, MTR031, MTR032
      double complex MTR034, MTR036, MTR037, MTR038, MTR039, MTR041
      double complex MTR043, MTR044, MTR046, MTR048, MTR049, MTR050
      double complex MTR051, MTR053, MTR058, MTR059, MTR061, MTR062
      double complex MTR063, MTR065, MTR066, MTR067, MTR069, MTR070
      double complex MTR071, MTR072, MTR074, MTR075, MTR076, MTR077
      double complex MTR079, MTR080, MTR081, MTR082, MTR088, MTR089
      double complex MTR090, MTR091, MTR092, MTR093, MTR094, MTR095
      double complex MTR096, MTR097, MTR098, MTR099, MTR100, MTR101
      double complex MTR103, MTR104, MTR107, MTR108, MTR109, MTR110
      double complex MTR111, MTR112, MTR113, MTR114, MTR115, MTR116
      double complex MTR117, MTR118, MTR119, MTR121, MTR122, MTR127
      double complex MTR128, MTR129, MTR130, MTR131, MTR132, MTR133
      double complex MTR134, MTR135, MTR136, MTR137, MTR138, MTR139
      double complex MTR140, MTR141, MTR142, MTR144, MTR145, MTR146
      double complex MTR147, MTR148, MTR149, MTR150, MTR151, MTR152
      double complex MTR153, MTR154, MTR155, MTR156, MTR157, MTR158
      double complex MTR159, MTR160, MTR161, MTR162, MTR163, MTR164
      double complex MTR165, MTR166, MTR167, MTR168, MTR169, MTR170
      double complex MTR171, MTR172, MTR173, MTR174, MTR175, MTR176
      double complex MTR177, MTR178, MTR179, MTR180, MTR181, MTR182
      double complex MTR183, MTR184, MTR185, MTR186, MTR187, MTR188
      double complex MTR189, MTR190, MTR191, MTR192, MTR193, MTR194
      double complex MTR195, MTR196, MTR200, MTR201, MTR202, MTR203
      double complex MTR204, MTR205, MTR206, MTR207, MTR211, MTR212
      double complex MTR213, MTR214, MTR215, MTR216, MTR221, MTR222
      double complex MTR223, MTR224, MTR225, MTR226, MTR227, MTR228
      double complex MTR229, MTR230, MTR231, MTR232, MTR233, MTR234
      double complex MTR239, MTR240, MTR241, MTR242, MTR243, MTR244
      double complex MTR245, MTR246, MTR247, MTR248, MTR249, MTR250
      double complex MTR255, MTR256, MTR257, MTR258, MTR259, MTR260
      double complex MTR261, MTR262, MTR263, MTR264, MTR265, MTR266
      double complex MTR267, MTR268

      dimension
     & quuMass(3), qudMass(3), lpdMass(3), neuMass(4), chaMass(2), 
     & sluMass(3), sldMass(3), sleMass(3), squMass(3), sqvMass(3), 
     & sqdMass(3), sqeMass(3), hisMass(3), MTR001(3), MTR002(3), 
     & MTR003(3), MTR004(3), MTR005(3), MTR006(3), MTR007(3), 
     & MTR008(3), MTR009(3), MTR010(3), MTR011(3), MTR012(3), 
     & MTR013(3), MTR014(3), MTR015(3), MTR016(3), MTR017(3), 
     & MTR018(3), MTR019(3), MTR020(3), MTR021(3,3), MTR022(3,3), 
     & MTR023(3), MTR024(3), MTR025(3), MTR026(3,3), MTR027(3,3), 
     & MTR028(3), MTR029(3), MTR030(3,3), MTR031(3), MTR032(3), 
     & MTR033(3), MTR034(3), MTR035(3,3), MTR036(3,3), MTR037(3), 
     & MTR038(3), MTR039(3), MTR040(3), MTR041(3,3), MTR042(3,3), 
     & MTR043(3), MTR044(3), MTR045(3), MTR046(3), MTR047(3,3), 
     & MTR048(3,3), MTR049(3), MTR050(3), MTR051(3), MTR052(3), 
     & MTR053(3,3), MTR054(3,3), MTR055(3), MTR056(3,3), 
     & MTR057(3,3,3), MTR058(3), MTR059(3), MTR060(3), MTR061(3), 
     & MTR062(3), MTR063(3), MTR064(3), MTR065(3), MTR066(3), 
     & MTR067(3), MTR068(3), MTR069(3), MTR070(3), MTR071(3), 
     & MTR072(3), MTR073(3), MTR074(3), MTR075(3), MTR076(3), 
     & MTR077(3), MTR078(3), MTR079(3), MTR080(3), MTR081(3), 
     & MTR082(3), MTR083(3), MTR084(3), MTR085(3), MTR086(3,3), 
     & MTR087(3), MTR088(2,2), MTR089(2,2), MTR090(2,2,3), 
     & MTR091(2,2,3), MTR092(2,3), MTR093(2,3), MTR094(2,4), 
     & MTR095(2,4), MTR096(2,4), MTR097(2,4), MTR098(2,3), 
     & MTR099(2,3), MTR100(2,3), MTR101(2,3), MTR102(3), 
     & MTR103(3,3), MTR104(3,3), MTR105(3), MTR106(3), MTR107(4,3), 
     & MTR108(4,3), MTR109(4,3), MTR110(4,3), MTR111(2,3), 
     & MTR112(2,3), MTR113(2,3), MTR114(2,3), MTR115(4,3), 
     & MTR116(4,3), MTR117(4,3), MTR118(4,3), MTR119(4,3), 
     & MTR120(3), MTR121(3,3), MTR122(3,3), MTR123(3), MTR124(3), 
     & MTR125(3), MTR126(3), MTR127(2,3), MTR128(2,3), MTR129(2,3), 
     & MTR130(2,3), MTR131(3), MTR132(3), MTR133(3), MTR134(3), 
     & MTR135(2,3), MTR136(2,3), MTR137(2,3), MTR138(2,3), 
     & MTR139(4,3), MTR140(4,3), MTR141(4,3), MTR142(4,3), 
     & MTR143(3), MTR144(3,3), MTR145(3,3), MTR146(3), MTR147(3), 
     & MTR148(3), MTR149(3), MTR150(2,3), MTR151(2,3), MTR152(2,3), 
     & MTR153(2,3), MTR154(4,3), MTR155(4,3), MTR156(4,3), 
     & MTR157(4,3), MTR158(3), MTR159(3), MTR160(3), MTR161(3), 
     & MTR162(4,3), MTR163(4,3), MTR164(4,3), MTR165(4,3), 
     & MTR166(3), MTR167(3), MTR168(3), MTR169(3), MTR170(2,3), 
     & MTR171(2,3), MTR172(4,3), MTR173(4,3), MTR174(4,3), 
     & MTR175(4,3), MTR176(4,3), MTR177(2,4), MTR178(2,4), 
     & MTR179(2,4), MTR180(2,4), MTR181(4,4), MTR182(4,4), 
     & MTR183(4,4,3), MTR184(4,4,3), MTR185(2,2), MTR186(2,2), 
     & MTR187(2,4), MTR188(2,4), MTR189(2,4), MTR190(2,4), 
     & MTR191(4,4), MTR192(4,4), MTR193(3), MTR194(3), MTR195(3), 
     & MTR196(3), MTR197(3), MTR198(3), MTR199(3), MTR200(3), 
     & MTR201(3), MTR202(3), MTR203(3), MTR204(3), MTR205(3), 
     & MTR206(3), MTR207(3), MTR208(3), MTR209(3), MTR210(3), 
     & MTR211(3), MTR212(3), MTR213(3), MTR214(3), MTR215(3), 
     & MTR216(3), MTR217(3), MTR218(3), MTR219(3), MTR220(3), 
     & MTR221(3), MTR222(3), MTR223(3), MTR224(3), MTR225(3), 
     & MTR226(3), MTR227(3), MTR228(3), MTR229(3), MTR230(3), 
     & MTR231(3), MTR232(3), MTR233(3), MTR234(3), MTR235(3), 
     & MTR236(3), MTR237(3), MTR238(3), MTR239(3), MTR240(3), 
     & MTR241(3), MTR242(3), MTR243(3), MTR244(3), MTR245(3), 
     & MTR246(3), MTR247(3), MTR248(3), MTR249(3), MTR250(3), 
     & MTR251(3), MTR252(3), MTR253(3), MTR254(3), MTR255(3), 
     & MTR256(3), MTR257(3), MTR258(3), MTR259(3), MTR260(3), 
     & MTR261(3), MTR262(3), MTR263(3), MTR264(3), MTR265(3), 
     & MTR266(3), MTR267(3), MTR268(3), MTR269(3), MTR270(3), 
     & MTR271(3), MTR272(3), MTR273(3), MTR274(3), hisW(3)

      common /mdl_mtrces/
     & quuMass, qudMass, lpdMass, neuMass, chaMass, sluMass, 
     & sldMass, sleMass, squMass, sqvMass, sqdMass, sqeMass, 
     & hisMass, MTR001, MTR002, MTR003, MTR004, MTR005, 
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
     & MTR272, MTR273, MTR274

      common /mdl_para/
     &    EE, SW, s12, s23, s13, GG, CW, S2W, C2W, MZ, MW, c12,
     &    c23, c13, Me, Mm, Mu, Md, Mc, Ms, wt, wZ, wW, tb, Ml,
     &    Mt, Mb, sb, cb, t2b, c2b, s2b, MG3, Ml1, Ml2, Ml3, Mr1,
     &    Mr2, Mr3, Ae, ls1, Am, ls2, MSne, MSnm, MSnl, MSeLL,
     &    MSeRR, MSeLR, MSe1, MSe2, MSeth, Zl11, Zl14, Zl41, Zl44,
     &    MSmuLL, MSmuLR, MSmuRR, MSm1, MSm2, MSmuth, Zl22, Zl25,
     &    Zl52, Zl55, Mq1, Mq2, Mq3, Mu1, Mu2, Mu3, Md1, Md2, Md3,
     &    Zual4, Zuas4, Zur0, ZuQMt, Zur1, Zulam5, Zur2, Zur3,
     &    Zulam6, Zur4, MtopS, Au, Ad, us1, ds1, As, ds2, Ac, us2,
     &    MSuLL, MSuLR, MSuRR, MSu1, MSu2, MSuth, Zu11, Zu14, Zu41,
     &    Zu44, MSdLL, MSdLR, MSdRR, MSd1, MSd2, MSdth, Zd11, Zd14,
     &    Zd41, Zd44, MScLL, MScLR, MScRR, MSc1, MSc2, MScth, Zu22,
     &    Zu25, Zu52, Zu55, MSsLL, MSsLR, MSsRR, MSs1, MSs2, MSsth,
     &    Zd22, Zd25, Zd52, Zd55, aMG1, fiMG1, aMG2, fiMG2, MG1r,
     &    MG1i, MG2r, MG2i, MG3i, amu, fimu, mur, mui, Alr, Ali,
     &    Atr, Ati, Abr, Abi, MHc, wHc, fMG1, fMG2, fmu, neutct,
     &    Zn11r, Zn11i, Zn12r, Zn12i, Zn13r, Zn13i, Zn14r, Zn14i,
     &    Zn21r, Zn21i, Zn22r, Zn22i, Zn23r, Zn23i, Zn24r, Zn24i,
     &    Zn31r, Zn31i, Zn32r, Zn32i, Zn33r, Zn33i, Zn34r, Zn34i,
     &    Zn41r, Zn41i, Zn42r, Zn42i, Zn43r, Zn43i, Zn44r, Zn44i,
     &    MNE1, wNE1, MNE2, wNE2, MNE3, wNE3, MNE4, wNE4, chatkn,
     &    Zm11r, Zp11r, Zm11i, Zp11i, Zm12r, Zp12r, Zm12i, Zp12i,
     &    Zm21r, Zp21r, Zm21i, Zp21i, Zm22r, Zp22r, Zm22i, Zp22i,
     &    MC1, wC1, MC2, wC2, MSG, wSG, hM11, hM12, hM22, hM33,
     &    htkn, Mh1, wh1, Mh2, wh2, Mh3, wh3, Zh11, Zh12, Zh13,
     &    Zh21, Zh22, Zh23, Zh31, Zh32, Zh33, MSbLL, MSbLRr, MSbLRi,
     &    MSbRR, MSb1, MSb2, MSbth, MSbfi, Zd33r, Zd36r, Zd63r,
     &    Zd66r, Zd33i, Zd36i, Zd63i, Zd66i, MStLL, MStLRr, MStLRi,
     &    MStRR, MSt1, MSt2, MStth, MStfi, Zu33r, Zu36r, Zu63r,
     &    Zu66r, Zu33i, Zu36i, Zu63i, Zu66i, MSlLL, MSlLRr, MSlLRi,
     &    MSlRR, MSl1, MSl2, MSlth, MSlfi, Zl33r, Zl36r, Zl63r,
     &    Zl66r, Zl33i, Zl36i, Zl63i, Zl66i, wSe1, wSe2, wSm1,
     &    wSm2, wSl1, wSl2, wSne, wSnm, wSnl, wSu1, wSu2, wSd1,
     &    wSd2, wSc1, wSc2, wSs1, wSs2, wSt1, wSt2, wSb1, wSb2,
     &    Mh1sq, Mh2sq, Mh3sq, Mesq, Mmsq, Mlsq, Musq, Mdsq, Mcsq,
     &    Mssq, Mbsq, Mtsq, MSe1sq, MSm1sq, MSl1sq, MSu1sq, MSd1sq,
     &    MSc1sq, MSs1sq, MSb1sq, MSt1sq, MSe2sq, MSm2sq, MSl2sq,
     &    MSu2sq, MSd2sq, MSc2sq, MSs2sq, MSb2sq, MSt2sq, MSnesq,
     &    MSnmsq, MSnlsq, AAABR, hisW

