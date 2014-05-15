*     LanHEP output produced at Tue May 18 00:00:25 2010
*     Model named 'NMSSM-Ug'

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

      double precision wtop, EE, SW, s12, s23, s13, GG, CW, S2W
      double precision C2W, MZ, MW, c12, c23, c13, Vud, Vus, Vcs
      double precision Vcd, Me, Mm, Mu, Md, Mc, Ms, Ml, Mt, Mb
      double precision wt, wZ, wW, tb, sb, cb, t2b, c2b, s2b, hL
      double precision hK, hLs, hKs, mue, xvev, MG2, MG3, MG1
      double precision Ml2, Ml3, Mr2, Mr3, Al, Ae, ls1, Am, ls2
      double precision ls3, MSne, MSnm, MSnl, MSeLL, MSeRR, MSeLR
      double precision MSe1, MSe2, MSeth, Zl11, Zl14, Zl41, Zl44
      double precision MSmuLL, MSmuLR, MSmuRR, MSm1, MSm2, MSmuth
      double precision Zl22, Zl25, Zl52, Zl55, MSlLL, MSlLR, MSlRR
      double precision MSl1, MSl2, MSlth, Zl33, Zl36, Zl63, Zl66
      double precision Mq2, Mq3, Mu2, Mu3, Md2, Md3, At, Ab, Au
      double precision Ad, us1, ds1, As, ds2, Ac, us2, ds3, us3
      double precision MSuLL, MSuLR, MSuRR, MSu1, MSu2, MSuth
      double precision Zu11, Zu14, Zu41, Zu44, MSdLL, MSdLR, MSdRR
      double precision MSd1, MSd2, MSdth, Zd11, Zd14, Zd41, Zd44
      double precision MScLL, MScLR, MScRR, MSc1, MSc2, MScth
      double precision Zu22, Zu25, Zu52, Zu55, MSsLL, MSsLR, MSsRR
      double precision MSs1, MSs2, MSsth, Zd22, Zd25, Zd52, Zd55
      double precision MStLL, MStLR, MStRR, MSt1, MSt2, MStth
      double precision Zu33, Zu36, Zu63, Zu66, MSbLL, MSbLR, MSbRR
      double precision MSb1, MSb2, MSbth, Zd33, Zd36, Zd63, Zd66
      double precision dlh2, wh1, wh2, wh3, wha, whb, wHc, ntkhs
      double precision Mh1, Mh2, Mh3, Zh11, Zh12, Zh13, Zh21, Zh22
      double precision Zh23, Zh31, Zh32, Zh33, ntkhp, Mha, Mhb
      double precision Za11, Za12, Za13, Za21, Za22, Za23, Za31
      double precision Za32, Za33, MHc, Zcsx, Zctx, MC01, MC02
      double precision Zcc2u, Zcc2v, Zcsigu, Zcsigv, Zccu, Zcsu
      double precision Zccv, Zcsv, Zcsig1, Zcsig2, Zcsig3, Zcsig4
      double precision MC1, MC2, Zm11, Zm12, Zm21, Zm22, Zp11
      double precision Zp12, Zp21, Zp22, neutk, Zn11r, Zn11i, Zn12r
      double precision Zn12i, Zn13r, Zn13i, Zn14r, Zn14i, Zn15r
      double precision Zn15i, Zn21r, Zn21i, Zn22r, Zn22i, Zn23r
      double precision Zn23i, Zn24r, Zn24i, Zn25r, Zn25i, Zn31r
      double precision Zn31i, Zn32r, Zn32i, Zn33r, Zn33i, Zn34r
      double precision Zn34i, Zn35r, Zn35i, Zn41r, Zn41i, Zn42r
      double precision Zn42i, Zn43r, Zn43i, Zn44r, Zn44i, Zn45r
      double precision Zn45i, Zn51r, Zn51i, Zn52r, Zn52i, Zn53r
      double precision Zn53i, Zn54r, Zn54i, Zn55r, Zn55i, MNE1
      double precision MNE2, MNE3, MNE4, MNE5, wC1, wC2, wNE2
      double precision wNE3, wNE4, MSG, wSG, wNE5, wSe1, wSe2
      double precision wSm1, wSm2, wSl1, wSl2, wSne, wSnm, wSnl
      double precision wSu1, wSu2, wSd1, wSd2, wSc1, wSc2, wSs1
      double precision wSs2, wSt1, wSt2, wSb1, wSb2, Mesq, Mmsq
      double precision Mlsq, Musq, Mdsq, Mcsq, Mssq, Mbsq, Mtsq
      double precision MSe1sq, MSm1sq, MSl1sq, MSu1sq, MSd1sq
      double precision MSc1sq, MSs1sq, MSb1sq, MSt1sq, MSe2sq
      double precision MSm2sq, MSl2sq, MSu2sq, MSd2sq, MSc2sq
      double precision MSs2sq, MSb2sq, MSt2sq, MSnesq, MSnmsq
      double precision MSnlsq

      double complex Zn11, Zn12, Zn13, Zn14, Zn15, Zn21, Zn22
      double complex Zn23, Zn24, Zn25, Zn31, Zn32, Zn33, Zn34
      double complex Zn35, Zn41, Zn42, Zn43, Zn44, Zn45, Zn51
      double complex Zn52, Zn53, Zn54, Zn55

      double precision AAABR(734)
      double complex AAABC(640)

      double precision quuMass, qudMass, lpdMass, neuMass, chaMass
      double precision sluMass, sldMass, sleMass, squMass, sqvMass
      double precision sqdMass, sqeMass, hisMass, hiaMass 
      double precision hisWidth,hiaWidth,MTR001
      double precision MTR002, MTR003, MTR004, MTR005, MTR006
      double precision MTR007, MTR008, MTR009, MTR010, MTR011
      double precision MTR012, MTR013, MTR014, MTR015, MTR016
      double precision MTR017, MTR018, MTR019, MTR020, MTR021
      double precision MTR022, MTR023, MTR024, MTR025, MTR026
      double precision MTR027, MTR028, MTR029, MTR030, MTR031
      double precision MTR032, MTR033, MTR034, MTR035, MTR036
      double precision MTR037, MTR038, MTR039, MTR040, MTR041
      double precision MTR042, MTR043, MTR044, MTR045, MTR046
      double precision MTR047, MTR048, MTR049, MTR050, MTR051
      double precision MTR052, MTR053, MTR054, MTR055, MTR056
      double precision MTR057, MTR058, MTR059, MTR060, MTR061
      double precision MTR062, MTR063, MTR064, MTR065, MTR066
      double precision MTR067, MTR068, MTR069, MTR070, MTR071
      double precision MTR072, MTR073, MTR074, MTR075, MTR076
      double precision MTR077, MTR078, MTR079, MTR080, MTR085
      double precision MTR086, MTR087, MTR088, MTR093, MTR094
      double precision MTR095, MTR096, MTR097, MTR102, MTR103
      double precision MTR113, MTR114, MTR115, MTR116, MTR117
      double precision MTR118, MTR119, MTR120, MTR121, MTR122
      double precision MTR123, MTR124, MTR125, MTR126, MTR127
      double precision MTR132, MTR133, MTR134, MTR135, MTR140
      double precision MTR141, MTR142, MTR143, MTR144, MTR145
      double precision MTR146, MTR147, MTR148, MTR149, MTR150
      double precision MTR203, MTR204, MTR209, MTR210, MTR219
      double precision MTR220, MTR221, MTR222, MTR223, MTR224
      double precision MTR225, MTR226, MTR227, MTR228, MTR229
      double precision MTR230, MTR231, MTR232, MTR233, MTR234
      double precision MTR235, MTR236, MTR237, MTR238, MTR239
      double precision MTR240, MTR241, MTR242, MTR243, MTR244
      double precision MTR245, MTR246, MTR247, MTR248, MTR249
      double precision MTR250, MTR251, MTR252, MTR253, MTR254
      double precision MTR255, MTR256, MTR257, MTR258, MTR259
      double precision MTR260, MTR261, MTR262, MTR263, MTR264
      double precision MTR265, MTR266, MTR267, MTR268, MTR269
      double precision MTR270, MTR271, MTR272, MTR273, MTR274
      double precision MTR275, MTR276, MTR277, MTR278, MTR279
      double precision MTR280, MTR281, MTR282, MTR283, MTR284
      double precision MTR285, MTR286, MTR287, MTR288, MTR289
      double precision MTR290, MTR291, MTR292, MTR293
      double complex MTR081, MTR082, MTR083, MTR084, MTR089, MTR090
      double complex MTR091, MTR092, MTR098, MTR099, MTR100, MTR101
      double complex MTR104, MTR105, MTR106, MTR107, MTR108, MTR109
      double complex MTR110, MTR111, MTR112, MTR128, MTR129, MTR130
      double complex MTR131, MTR136, MTR137, MTR138, MTR139, MTR151
      double complex MTR152, MTR153, MTR154, MTR155, MTR156, MTR157
      double complex MTR158, MTR159, MTR160, MTR161, MTR162, MTR163
      double complex MTR164, MTR165, MTR166, MTR167, MTR168, MTR169
      double complex MTR170, MTR171, MTR172, MTR173, MTR174, MTR175
      double complex MTR176, MTR177, MTR178, MTR179, MTR180, MTR181
      double complex MTR182, MTR183, MTR184, MTR185, MTR186, MTR187
      double complex MTR188, MTR189, MTR190, MTR191, MTR192, MTR193
      double complex MTR194, MTR195, MTR196, MTR197, MTR198, MTR199
      double complex MTR200, MTR201, MTR202, MTR205, MTR206, MTR207
      double complex MTR208, MTR211, MTR212, MTR213, MTR214, MTR215
      double complex MTR216, MTR217, MTR218

      dimension
     & quuMass(3), qudMass(3), lpdMass(3), neuMass(4), chaMass(2), 
     & sluMass(3), sldMass(3), sleMass(3), squMass(3), sqvMass(3), 
     & sqdMass(3), sqeMass(3), hisMass(3), hiaMass(2), 
     & hisWidth(3), hiaWidth(2), MTR001(3), 
     & MTR002(3), MTR003(3), MTR004(3), MTR005(3), MTR006(3,3), 
     & MTR007(3,3), MTR008(3,3), MTR009(3,3), MTR010(2), 
     & MTR011(3), MTR012(3,3), MTR013(3,3), MTR014(3,3), 
     & MTR015(3,3), MTR016(3), MTR017(3), MTR018(2,3), MTR019(3,3), 
     & MTR020(3,3), MTR021(3), MTR022(3,3), MTR023(3,3), 
     & MTR024(3,3), MTR025(3,3), MTR026(3), MTR027(2,3), 
     & MTR028(3,3), MTR029(3,3), MTR030(3,3), MTR031(3,3), 
     & MTR032(3,3), MTR033(3,3), MTR034(3,3), MTR035(3), 
     & MTR036(2,3), MTR037(3,3), MTR038(3,3), MTR039(3,3), 
     & MTR040(3,3), MTR041(3,3), MTR042(3), MTR043(3), MTR044(2,3), 
     & MTR045(2,2,3), MTR046(3,3,3), MTR047(2), MTR048(3), 
     & MTR049(3), MTR050(3), MTR051(3), MTR052(3), MTR053(3), 
     & MTR054(3), MTR055(3), MTR056(3,3), MTR057(3,3), MTR058(3), 
     & MTR059(3,3), MTR060(3,3), MTR061(3,3), MTR062(3,3), 
     & MTR063(3), MTR064(3), MTR065(3,3), MTR066(3,3), MTR067(3), 
     & MTR068(2), MTR069(3), MTR070(3), MTR071(2,3), MTR072(3), 
     & MTR073(2,2), MTR074(2,2), MTR075(2,2,2), MTR076(2,2,2), 
     & MTR077(2,2,3), MTR078(2,2,3), MTR079(2,3), MTR080(2,3), 
     & MTR081(2,4), MTR082(2,4), MTR083(2,4), MTR084(2,4), 
     & MTR085(2,3,3), MTR086(2,3,3), MTR087(2,3,3), MTR088(2,3,3), 
     & MTR089(2), MTR090(2), MTR091(2), MTR092(2), MTR093(3), 
     & MTR094(3,2), MTR095(3,3), MTR096(3), MTR097(3), MTR098(4,3), 
     & MTR099(4,3), MTR100(4,3), MTR101(4,3), MTR102(2,3), 
     & MTR103(2,3), MTR104(3), MTR105(3), MTR106(3), MTR107(3), 
     & MTR108(4,3), MTR109(4,3), MTR110(4,3), MTR111(4,3), 
     & MTR112(4,3), MTR113(3), MTR114(3,2), MTR115(3,3), 
     & MTR116(3,3), MTR117(3,3), MTR118(3,3), MTR119(3,3), 
     & MTR120(3,2,3), MTR121(3,2,3), MTR122(3,2,3), MTR123(3,2,3), 
     & MTR124(3), MTR125(3), MTR126(3), MTR127(3), MTR128(3), 
     & MTR129(3), MTR130(3), MTR131(3), MTR132(3,2,3), MTR133(3,2,3), 
     & MTR134(3,2,3), MTR135(3,2,3), MTR136(4,3), MTR137(4,3), 
     & MTR138(4,3), MTR139(4,3), MTR140(3,3), MTR141(3,3), 
     & MTR142(3,3), MTR143(3,3), MTR144(3), MTR145(3,2), 
     & MTR146(3,3), MTR147(3), MTR148(3), MTR149(3), MTR150(3), 
     & MTR151(3), MTR152(3), MTR153(3), MTR154(3), MTR155(4,3), 
     & MTR156(4,3), MTR157(4,3), MTR158(4,3), MTR159(3), 
     & MTR160(3), MTR161(3), MTR162(3), MTR163(4,3), MTR164(4,3), 
     & MTR165(4,3), MTR166(4,3), MTR167(3), MTR168(3), MTR169(3), 
     & MTR170(3), MTR171(4,3), MTR172(4,3), MTR173(4,3), 
     & MTR174(4,3), MTR175(3), MTR176(3), MTR177(3), MTR178(3), 
     & MTR179(2,4), MTR180(2,4), MTR181(2,4), MTR182(2,4), 
     & MTR183(2), MTR184(2), MTR185(2), MTR186(2), MTR187(4,4), 
     & MTR188(4,4), MTR189(4,4,2), MTR190(4,4,2), MTR191(4,4,3), 
     & MTR192(4,4,3), MTR193(4), MTR194(4), MTR195(4,2), 
     & MTR196(4,2), MTR197(4,3), MTR198(4,3), MTR199(2), 
     & MTR200(2), MTR201(3), MTR202(3), MTR203(2,2), MTR204(2,2), 
     & MTR205(2,4), MTR206(2,4), MTR207(2), MTR208(2), MTR209(3,3), 
     & MTR210(3,3), MTR211(2,4), MTR212(2,4), MTR213(2), 
     & MTR214(2), MTR215(4,4), MTR216(4,4), MTR217(4), MTR218(4), 
     & MTR219(2), MTR220(2), MTR221(3), MTR222(3), MTR223(3), 
     & MTR224(3), MTR225(3), MTR226(3), MTR227(3), MTR228(3), 
     & MTR229(3), MTR230(3), MTR231(3), MTR232(3), MTR233(3), 
     & MTR234(3), MTR235(3), MTR236(3), MTR237(3), MTR238(3), 
     & MTR239(3), MTR240(3), MTR241(3), MTR242(3), MTR243(3), 
     & MTR244(3,3), MTR245(3,3), MTR246(3,3), MTR247(3,3), 
     & MTR248(3,3), MTR249(3,3), MTR250(3), MTR251(3), MTR252(3), 
     & MTR253(3), MTR254(3,3), MTR255(3,3), MTR256(3,3), 
     & MTR257(3,3), MTR258(3,3), MTR259(3,3), MTR260(3,3), 
     & MTR261(3,3), MTR262(3,3), MTR263(3,3), MTR264(3,3), 
     & MTR265(3,3), MTR266(3), MTR267(3), MTR268(3), MTR269(3), 
     & MTR270(3), MTR271(3), MTR272(3), MTR273(3), MTR274(3,3), 
     & MTR275(3,3), MTR276(3,3), MTR277(3,3), MTR278(3,3), 
     & MTR279(3,3), MTR280(3), MTR281(3), MTR282(3), MTR283(3), 
     & MTR284(2), MTR285(2), MTR286(3), MTR287(3), MTR288(2), 
     & MTR289(2), MTR290(2,2), MTR291(2,2), MTR292(3,3), 
     & MTR293(3,3)

      common /mdl_mtrces/
     & quuMass, qudMass, lpdMass, neuMass, chaMass, sluMass, 
     & sldMass, sleMass, squMass, sqvMass, sqdMass, sqeMass, 
     & hisMass, hiaMass, hisWidth, hiaWidth, 
     & MTR001, MTR002, MTR003, MTR004, 
     & MTR005, MTR006, MTR007, MTR008, MTR009, MTR010, MTR011, 
     & MTR012, MTR013, MTR014, MTR015, MTR016, MTR017, MTR018, 
     & MTR019, MTR020, MTR021, MTR022, MTR023, MTR024, MTR025, 
     & MTR026, MTR027, MTR028, MTR029, MTR030, MTR031, MTR032, 
     & MTR033, MTR034, MTR035, MTR036, MTR037, MTR038, MTR039, 
     & MTR040, MTR041, MTR042, MTR043, MTR044, MTR045, MTR046, 
     & MTR047, MTR048, MTR049, MTR050, MTR051, MTR052, MTR053, 
     & MTR054, MTR055, MTR056, MTR057, MTR058, MTR059, MTR060, 
     & MTR061, MTR062, MTR063, MTR064, MTR065, MTR066, MTR067, 
     & MTR068, MTR069, MTR070, MTR071, MTR072, MTR073, MTR074, 
     & MTR075, MTR076, MTR077, MTR078, MTR079, MTR080, MTR081, 
     & MTR082, MTR083, MTR084, MTR085, MTR086, MTR087, MTR088, 
     & MTR089, MTR090, MTR091, MTR092, MTR093, MTR094, MTR095, 
     & MTR096, MTR097, MTR098, MTR099, MTR100, MTR101, MTR102, 
     & MTR103, MTR104, MTR105, MTR106, MTR107, MTR108, MTR109, 
     & MTR110, MTR111, MTR112, MTR113, MTR114, MTR115, MTR116, 
     & MTR117, MTR118, MTR119, MTR120, MTR121, MTR122, MTR123, 
     & MTR124, MTR125, MTR126, MTR127, MTR128, MTR129, MTR130, 
     & MTR131, MTR132, MTR133, MTR134, MTR135, MTR136, MTR137, 
     & MTR138, MTR139, MTR140, MTR141, MTR142, MTR143, MTR144, 
     & MTR145, MTR146, MTR147, MTR148, MTR149, MTR150, MTR151, 
     & MTR152, MTR153, MTR154, MTR155, MTR156, MTR157, MTR158, 
     & MTR159, MTR160, MTR161, MTR162, MTR163, MTR164, MTR165, 
     & MTR166, MTR167, MTR168, MTR169, MTR170, MTR171, MTR172, 
     & MTR173, MTR174, MTR175, MTR176, MTR177, MTR178, MTR179, 
     & MTR180, MTR181, MTR182, MTR183, MTR184, MTR185, MTR186, 
     & MTR187, MTR188, MTR189, MTR190, MTR191, MTR192, MTR193, 
     & MTR194, MTR195, MTR196, MTR197, MTR198, MTR199, MTR200, 
     & MTR201, MTR202, MTR203, MTR204, MTR205, MTR206, MTR207, 
     & MTR208, MTR209, MTR210, MTR211, MTR212, MTR213, MTR214, 
     & MTR215, MTR216, MTR217, MTR218, MTR219, MTR220, MTR221, 
     & MTR222, MTR223, MTR224, MTR225, MTR226, MTR227, MTR228, 
     & MTR229, MTR230, MTR231, MTR232, MTR233, MTR234, MTR235, 
     & MTR236, MTR237, MTR238, MTR239, MTR240, MTR241, MTR242, 
     & MTR243, MTR244, MTR245, MTR246, MTR247, MTR248, MTR249, 
     & MTR250, MTR251, MTR252, MTR253, MTR254, MTR255, MTR256, 
     & MTR257, MTR258, MTR259, MTR260, MTR261, MTR262, MTR263, 
     & MTR264, MTR265, MTR266, MTR267, MTR268, MTR269, MTR270, 
     & MTR271, MTR272, MTR273, MTR274, MTR275, MTR276, MTR277, 
     & MTR278, MTR279, MTR280, MTR281, MTR282, MTR283, MTR284, 
     & MTR285, MTR286, MTR287, MTR288, MTR289, MTR290, MTR291, 
     & MTR292, MTR293

      common /mdl_para/
     &    wtop, EE, SW, s12, s23, s13, GG, CW, S2W, C2W, MZ, MW,
     &    c12, c23, c13, Vud, Vus, Vcs, Vcd, Me, Mm, Mu, Md, Mc,
     &    Ms, Ml, Mt, Mb, wt, wZ, wW, tb, sb, cb, t2b, c2b, s2b,
     &    hL, hK, hLs, hKs, mue, xvev, MG2, MG3, MG1, Ml2, Ml3,
     &    Mr2, Mr3, Al, Ae, ls1, Am, ls2, ls3, MSne, MSnm, MSnl,
     &    MSeLL, MSeRR, MSeLR, MSe1, MSe2, MSeth, Zl11, Zl14, Zl41,
     &    Zl44, MSmuLL, MSmuLR, MSmuRR, MSm1, MSm2, MSmuth, Zl22,
     &    Zl25, Zl52, Zl55, MSlLL, MSlLR, MSlRR, MSl1, MSl2, MSlth,
     &    Zl33, Zl36, Zl63, Zl66, Mq2, Mq3, Mu2, Mu3, Md2, Md3,
     &    At, Ab, Au, Ad, us1, ds1, As, ds2, Ac, us2, ds3, us3,
     &    MSuLL, MSuLR, MSuRR, MSu1, MSu2, MSuth, Zu11, Zu14, Zu41,
     &    Zu44, MSdLL, MSdLR, MSdRR, MSd1, MSd2, MSdth, Zd11, Zd14,
     &    Zd41, Zd44, MScLL, MScLR, MScRR, MSc1, MSc2, MScth, Zu22,
     &    Zu25, Zu52, Zu55, MSsLL, MSsLR, MSsRR, MSs1, MSs2, MSsth,
     &    Zd22, Zd25, Zd52, Zd55, MStLL, MStLR, MStRR, MSt1, MSt2,
     &    MStth, Zu33, Zu36, Zu63, Zu66, MSbLL, MSbLR, MSbRR, MSb1,
     &    MSb2, MSbth, Zd33, Zd36, Zd63, Zd66, dlh2, wh1, wh2,
     &    wh3, wha, whb, wHc, ntkhs, Mh1, Mh2, Mh3, Zh11, Zh12,
     &    Zh13, Zh21, Zh22, Zh23, Zh31, Zh32, Zh33, ntkhp, Mha,
     &    Mhb, Za11, Za12, Za13, Za21, Za22, Za23, Za31, Za32,
     &    Za33, MHc, Zcsx, Zctx, MC01, MC02, Zcc2u, Zcc2v, Zcsigu,
     &    Zcsigv, Zccu, Zcsu, Zccv, Zcsv, Zcsig1, Zcsig2, Zcsig3,
     &    Zcsig4, MC1, MC2, Zm11, Zm12, Zm21, Zm22, Zp11, Zp12,
     &    Zp21, Zp22, neutk, Zn11r, Zn11i, Zn12r, Zn12i, Zn13r,
     &    Zn13i, Zn14r, Zn14i, Zn15r, Zn15i, Zn21r, Zn21i, Zn22r,
     &    Zn22i, Zn23r, Zn23i, Zn24r, Zn24i, Zn25r, Zn25i, Zn31r,
     &    Zn31i, Zn32r, Zn32i, Zn33r, Zn33i, Zn34r, Zn34i, Zn35r,
     &    Zn35i, Zn41r, Zn41i, Zn42r, Zn42i, Zn43r, Zn43i, Zn44r,
     &    Zn44i, Zn45r, Zn45i, Zn51r, Zn51i, Zn52r, Zn52i, Zn53r,
     &    Zn53i, Zn54r, Zn54i, Zn55r, Zn55i, MNE1, MNE2, MNE3,
     &    MNE4, MNE5, Zn11, Zn12, Zn13, Zn14, Zn15, Zn21, Zn22,
     &    Zn23, Zn24, Zn25, Zn31, Zn32, Zn33, Zn34, Zn35, Zn41,
     &    Zn42, Zn43, Zn44, Zn45, Zn51, Zn52, Zn53, Zn54, Zn55,
     &    wC1, wC2, wNE2, wNE3, wNE4, MSG, wSG, wNE5, wSe1, wSe2,
     &    wSm1, wSm2, wSl1, wSl2, wSne, wSnm, wSnl, wSu1, wSu2,
     &    wSd1, wSd2, wSc1, wSc2, wSs1, wSs2, wSt1, wSt2, wSb1,
     &    wSb2, Mesq, Mmsq, Mlsq, Musq, Mdsq, Mcsq, Mssq, Mbsq,
     &    Mtsq, MSe1sq, MSm1sq, MSl1sq, MSu1sq, MSd1sq, MSc1sq,
     &    MSs1sq, MSb1sq, MSt1sq, MSe2sq, MSm2sq, MSl2sq, MSu2sq,
     &    MSd2sq, MSc2sq, MSs2sq, MSb2sq, MSt2sq, MSnesq, MSnmsq,
     &    MSnlsq, AAABR, AAABC

