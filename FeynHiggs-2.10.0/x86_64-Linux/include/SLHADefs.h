#if 0
	SLHADefs.h
		declarations for SLHALib data
		generated 1 Jun 2012 9:36
#endif

#ifndef SLHADEFS_H
#define SLHADEFS_H

#define invalid (-999)

#define OffsetModSel 0
#define LengthModSel 6
#define BlockModSel(i) SlhaData(i)
#define ModSel_Model Slhadata(1)
#define ModSel_GridPts Slhadata(2)
#define ModSel_Content Slhadata(3)
#define ModSel_RPV Slhadata(4)
#define ModSel_CPV Slhadata(5)
#define ModSel_FV Slhadata(6)

#define OffsetSMInputs 6
#define LengthSMInputs 16
#define BlockSMInputs(i) SlhaData(6+i)
#define SMInputs_invAlfaMZ Slhadata(7)
#define SMInputs_GF Slhadata(8)
#define SMInputs_AlfasMZ Slhadata(9)
#define SMInputs_MZ Slhadata(10)
#define SMInputs_Mf(t,g) Slhadata(6+t+4*(g))
#define SMInputs_MfFlat(i) Slhadata(10+i)
#define   SMInputs_Mnu1 SMInputs_Mf(1,1)
#define   SMInputs_Me SMInputs_Mf(2,1)
#define   SMInputs_Mu SMInputs_Mf(3,1)
#define   SMInputs_Md SMInputs_Mf(4,1)
#define   SMInputs_Mnu2 SMInputs_Mf(1,2)
#define   SMInputs_Mmu SMInputs_Mf(2,2)
#define   SMInputs_Mc SMInputs_Mf(3,2)
#define   SMInputs_Ms SMInputs_Mf(4,2)
#define   SMInputs_Mnu3 SMInputs_Mf(1,3)
#define   SMInputs_Mtau SMInputs_Mf(2,3)
#define   SMInputs_Mt SMInputs_Mf(3,3)
#define   SMInputs_Mb SMInputs_Mf(4,3)

#define OffsetMinPar 22
#define LengthMinPar 6
#define BlockMinPar(i) SlhaData(22+i)
#define MinPar_M0 Slhadata(23)
#define   MinPar_Lambda MinPar_M0
#define MinPar_M12 Slhadata(24)
#define   MinPar_Mmess MinPar_M12
#define   MinPar_M32 MinPar_M12
#define MinPar_TB Slhadata(25)
#define MinPar_signMUE Slhadata(26)
#define MinPar_A Slhadata(27)
#define   MinPar_N5 MinPar_A
#define MinPar_cgrav Slhadata(28)

#define OffsetExtPar 28
#define LengthExtPar 42
#define BlockExtPar(i) SlhaData(28+i)
#define ExtPar_Q SlhaData(29)
#define ExtPar_M1 Slhadata(30)
#define ExtPar_M2 Slhadata(31)
#define ExtPar_M3 Slhadata(32)
#define ExtPar_Af(t) Slhadata(31+t)
#define   ExtPar_Atau ExtPar_Af(2)
#define   ExtPar_At ExtPar_Af(3)
#define   ExtPar_Ab ExtPar_Af(4)
#define ExtPar_MHu2 Slhadata(36)
#define ExtPar_MHd2 Slhadata(37)
#define ExtPar_MUE Slhadata(38)
#define ExtPar_MA02 Slhadata(39)
#define ExtPar_TB Slhadata(40)
#define ExtPar_MA0 Slhadata(41)
#define ExtPar_MHp Slhadata(42)
#define ExtPar_MSS(g,q) Slhadata(39+g+3*(q))
#define   ExtPar_MSL(g) ExtPar_MSS(g,1)
#define   ExtPar_MSE(g) ExtPar_MSS(g,2)
#define   ExtPar_MSQ(g) ExtPar_MSS(g,3)
#define   ExtPar_MSU(g) ExtPar_MSS(g,4)
#define   ExtPar_MSD(g) ExtPar_MSS(g,5)
#define ExtPar_N5(g) Slhadata(57+g)
#define ExtPar_lambda Slhadata(61)
#define ExtPar_kappa Slhadata(62)
#define ExtPar_Alambda Slhadata(63)
#define ExtPar_Akappa Slhadata(64)
#define ExtPar_lambdaS Slhadata(65)
#define ExtPar_xiF Slhadata(66)
#define ExtPar_xiS Slhadata(67)
#define ExtPar_MUEprime Slhadata(68)
#define ExtPar_mS2prime Slhadata(69)
#define ExtPar_mS2 Slhadata(70)

#define OffsetQExtPar 70
#define LengthQExtPar 16
#define BlockQExtPar(i) SlhaData(70+i)
#define QExtPar_QM1 SlhaData(71)
#define QExtPar_QM2 SlhaData(72)
#define QExtPar_QM3 SlhaData(73)
#define QExtPar_QAf(t) SlhaData(72+t)
#define   QExtPar_QAtau QExtPar_QAf(2)
#define   QExtPar_QAt QExtPar_QAf(3)
#define   QExtPar_QAb QExtPar_QAf(4)
#define QExtPar_QMHu2 SlhaData(77)
#define QExtPar_QMHd2 SlhaData(78)
#define QExtPar_QMUE SlhaData(79)
#define QExtPar_QMA02 SlhaData(80)
#define QExtPar_QTB SlhaData(81)
#define QExtPar_QMSS(q) SlhaData(81+q)
#define   QExtPar_QMSL QExtPar_QMSS(1)
#define   QExtPar_QMSE QExtPar_QMSS(2)
#define   QExtPar_QMSQ QExtPar_QMSS(3)
#define   QExtPar_QMSU QExtPar_QMSS(4)
#define   QExtPar_QMSD QExtPar_QMSS(5)

#define OffsetNMSSMRun 86
#define LengthNMSSMRun 11
#define BlockNMSSMRun(i) SlhaData(86+i)
#define NMSSMRun_Q SlhaData(87)
#define NMSSMRun_lambda Slhadata(88)
#define NMSSMRun_kappa Slhadata(89)
#define NMSSMRun_Alambda Slhadata(90)
#define NMSSMRun_Akappa Slhadata(91)
#define NMSSMRun_lambdaS Slhadata(92)
#define NMSSMRun_xiF Slhadata(93)
#define NMSSMRun_xiS Slhadata(94)
#define NMSSMRun_MUEprime Slhadata(95)
#define NMSSMRun_mS2prime Slhadata(96)
#define NMSSMRun_mS2 Slhadata(97)

#define OffsetMass 97
#define LengthMass 53
#define BlockMass(i) SlhaData(97+i)
#define Mass_Mf(t,g) Slhadata(93+t+4*(g))
#define Mass_MfFlat(i) Slhadata(97+i)
#define Mass_MSf(s,t,g) Slhadata(99+s+8*(g)+2*(t))
#define Mass_MSfFlat(i) Slhadata(109+i)
#define Mass_MZ Slhadata(134)
#define Mass_MW Slhadata(135)
#define Mass_Mh0 Slhadata(136)
#define Mass_MHH Slhadata(137)
#define Mass_MA0 Slhadata(138)
#define Mass_MHp Slhadata(139)
#define   Mass_MH1 Mass_Mh0
#define   Mass_MH2 Mass_MHH
#define Mass_MH3 Slhadata(140)
#define   Mass_MA1 Mass_MA0
#define Mass_MA2 Slhadata(141)
#define Mass_MNeu(n) Slhadata(141+n)
#define Mass_MCha(c) Slhadata(146+c)
#define Mass_MGl Slhadata(149)
#define Mass_MGrav Slhadata(150)

#define OffsetDMass 150
#define LengthDMass 5
#define BlockDMass(i) SlhaData(150+i)
#define DMass_Q SlhaData(151)
#define DMass_DeltaMh0 Slhadata(152)
#define DMass_DeltaMHH Slhadata(153)
#define DMass_DeltaMA0 Slhadata(154)
#define DMass_DeltaMHp Slhadata(155)

#define OffsetNMix 155
#define LengthNMix 16
#define BlockNMix(i) SlhaData(155+i)
#define NMix_ZNeu(n1,n2) Slhadata(151+n1+4*(n2))
#define NMix_ZNeuFlat(i) Slhadata(155+i)

#define OffsetUMix 171
#define LengthUMix 4
#define BlockUMix(i) SlhaData(171+i)
#define UMix_UCha(c1,c2) Slhadata(169+c1+2*(c2))
#define UMix_UChaFlat(i) Slhadata(171+i)

#define OffsetVMix 175
#define LengthVMix 4
#define BlockVMix(i) SlhaData(175+i)
#define VMix_VCha(c1,c2) Slhadata(173+c1+2*(c2))
#define VMix_VChaFlat(i) Slhadata(175+i)

#define OffsetSfMix 179
#define LengthSfMix 12
#define BlockSfMix(i) SlhaData(179+i)
#define SfMix_USf(s1,s2,t) Slhadata(169+s1+2*(s2)+4*(t))
#define SfMix_USfFlat(i,t) Slhadata(171+i+4*(t))

#define OffsetStauMix 179
#define LengthStauMix 4
#define BlockStauMix(i) SlhaData(179+i)
#define   StauMix_USf(s1,s2) SfMix_USf(s1,s2,2)
#define   StauMix_USfFlat(i) SfMix_USfFlat(i,2)

#define OffsetStopMix 183
#define LengthStopMix 4
#define BlockStopMix(i) SlhaData(183+i)
#define   StopMix_USf(s1,s2) SfMix_USf(s1,s2,3)
#define   StopMix_USfFlat(i) SfMix_USfFlat(i,3)

#define OffsetSbotMix 187
#define LengthSbotMix 4
#define BlockSbotMix(i) SlhaData(187+i)
#define   SbotMix_USf(s1,s2) SfMix_USf(s1,s2,4)
#define   SbotMix_USfFlat(i) SfMix_USfFlat(i,4)

#define OffsetAlpha 191
#define LengthAlpha 1
#define BlockAlpha(i) SlhaData(191+i)
#define Alpha_Alpha Slhadata(192)

#define OffsetDAlpha 192
#define LengthDAlpha 1
#define BlockDAlpha(i) SlhaData(192+i)
#define DAlpha_DeltaAlpha Slhadata(193)

#define OffsetHMix 193
#define LengthHMix 5
#define BlockHMix(i) SlhaData(193+i)
#define HMix_Q SlhaData(194)
#define HMix_MUE Slhadata(195)
#define HMix_TB Slhadata(196)
#define HMix_VEV Slhadata(197)
#define HMix_MA02 Slhadata(198)

#define OffsetGauge 198
#define LengthGauge 4
#define BlockGauge(i) SlhaData(198+i)
#define Gauge_Q SlhaData(199)
#define Gauge_g1 Slhadata(200)
#define Gauge_g2 Slhadata(201)
#define Gauge_g3 Slhadata(202)

#define OffsetMSoft 202
#define LengthMSoft 21
#define BlockMSoft(i) SlhaData(202+i)
#define MSoft_Q SlhaData(203)
#define MSoft_M1 Slhadata(204)
#define MSoft_M2 Slhadata(205)
#define MSoft_M3 Slhadata(206)
#define MSoft_MHu2 Slhadata(207)
#define MSoft_MHd2 Slhadata(208)
#define MSoft_MSS(g,q) Slhadata(205+g+3*(q))
#define   MSoft_MSL(g) MSoft_MSS(g,1)
#define   MSoft_MSE(g) MSoft_MSS(g,2)
#define   MSoft_MSQ(g) MSoft_MSS(g,3)
#define   MSoft_MSU(g) MSoft_MSS(g,4)
#define   MSoft_MSD(g) MSoft_MSS(g,5)

#define OffsetAf 223
#define LengthAf 30
#define BlockAf(i) SlhaData(223+i)
#define Af_Q(t) SlhaData(204+10*(t))
#define Af_Af(g1,g2,t) Slhadata(201+g1+3*(g2)+10*(t))
#define Af_AfFlat(i,t) Slhadata(204+i+10*(t))

#define OffsetAe 223
#define LengthAe 11
#define BlockAe(i) SlhaData(223+i)
#define   Ae_Q Af_Q(2)
#define   Ae_Af(g1,g2) Af_Af(g1,g2,2)
#define   Ae_AfFlat(i) Af_AfFlat(i,2)
#define   Ae_Atau Ae_Af(3,3)

#define OffsetAu 234
#define LengthAu 11
#define BlockAu(i) SlhaData(234+i)
#define   Au_Q Af_Q(3)
#define   Au_Af(g1,g2) Af_Af(g1,g2,3)
#define   Au_AfFlat(i) Af_AfFlat(i,3)
#define   Au_At Au_Af(3,3)

#define OffsetAd 245
#define LengthAd 11
#define BlockAd(i) SlhaData(245+i)
#define   Ad_Q Af_Q(4)
#define   Ad_Af(g1,g2) Af_Af(g1,g2,4)
#define   Ad_AfFlat(i) Af_AfFlat(i,4)
#define   Ad_Ab Ad_Af(3,3)

#define OffsetYf 256
#define LengthYf 30
#define BlockYf(i) SlhaData(256+i)
#define Yf_Q(t) SlhaData(237+10*(t))
#define Yf_Yf(g1,g2,t) Slhadata(234+g1+3*(g2)+10*(t))
#define Yf_YfFlat(i,t) Slhadata(237+i+10*(t))

#define OffsetYe 256
#define LengthYe 11
#define BlockYe(i) SlhaData(256+i)
#define   Ye_Q Yf_Q(2)
#define   Ye_Yf(g1,g2) Yf_Yf(g1,g2,2)
#define   Ye_YfFlat(i) Yf_YfFlat(i,2)
#define   Ye_Ytau Ye_Yf(3,3)

#define OffsetYu 267
#define LengthYu 11
#define BlockYu(i) SlhaData(267+i)
#define   Yu_Q Yf_Q(3)
#define   Yu_Yf(g1,g2) Yf_Yf(g1,g2,3)
#define   Yu_YfFlat(i) Yf_YfFlat(i,3)
#define   Yu_Yt Yu_Yf(3,3)

#define OffsetYd 278
#define LengthYd 11
#define BlockYd(i) SlhaData(278+i)
#define   Yd_Q Yf_Q(4)
#define   Yd_Yf(g1,g2) Yf_Yf(g1,g2,4)
#define   Yd_YfFlat(i) Yf_YfFlat(i,4)
#define   Yd_Yb Yd_Yf(3,3)

#define OffsetRVLamLLEIn 289
#define LengthRVLamLLEIn 27
#define BlockRVLamLLEIn(i) SlhaData(289+i)
#define RVLamLLEIn_lamLLE(i,j,k) Slhadata(277+i+3*(j)+9*(k))
#define RVLamLLEIn_lamLLEFlat(i) Slhadata(289+i)

#define OffsetRVLamLQDIn 316
#define LengthRVLamLQDIn 27
#define BlockRVLamLQDIn(i) SlhaData(316+i)
#define RVLamLQDIn_lamLQD(i,j,k) Slhadata(304+i+3*(j)+9*(k))
#define RVLamLQDIn_lamLQDFlat(i) Slhadata(316+i)

#define OffsetRVLamUDDIn 343
#define LengthRVLamUDDIn 27
#define BlockRVLamUDDIn(i) SlhaData(343+i)
#define RVLamUDDIn_lamUDD(i,j,k) Slhadata(331+i+3*(j)+9*(k))
#define RVLamUDDIn_lamUDDFlat(i) Slhadata(343+i)

#define OffsetRVLamLLE 370
#define LengthRVLamLLE 28
#define BlockRVLamLLE(i) SlhaData(370+i)
#define RVLamLLE_Q SlhaData(371)
#define RVLamLLE_lamLLE(i,j,k) Slhadata(359+i+3*(j)+9*(k))
#define RVLamLLE_lamLLEFlat(i) Slhadata(371+i)

#define OffsetRVLamLQD 398
#define LengthRVLamLQD 28
#define BlockRVLamLQD(i) SlhaData(398+i)
#define RVLamLQD_Q SlhaData(399)
#define RVLamLQD_lamLQD(i,j,k) Slhadata(387+i+3*(j)+9*(k))
#define RVLamLQD_lamLQDFlat(i) Slhadata(399+i)

#define OffsetRVLamUDD 426
#define LengthRVLamUDD 28
#define BlockRVLamUDD(i) SlhaData(426+i)
#define RVLamUDD_Q SlhaData(427)
#define RVLamUDD_lamUDD(i,j,k) Slhadata(415+i+3*(j)+9*(k))
#define RVLamUDD_lamUDDFlat(i) Slhadata(427+i)

#define OffsetRVTLLEIn 454
#define LengthRVTLLEIn 27
#define BlockRVTLLEIn(i) SlhaData(454+i)
#define RVTLLEIn_TLLE(i,j,k) Slhadata(442+i+3*(j)+9*(k))
#define RVTLLEIn_TLLEFlat(i) Slhadata(454+i)

#define OffsetRVTLQDIn 481
#define LengthRVTLQDIn 27
#define BlockRVTLQDIn(i) SlhaData(481+i)
#define RVTLQDIn_TLQD(i,j,k) Slhadata(469+i+3*(j)+9*(k))
#define RVTLQDIn_TLQDFlat(i) Slhadata(481+i)

#define OffsetRVTUDDIn 508
#define LengthRVTUDDIn 27
#define BlockRVTUDDIn(i) SlhaData(508+i)
#define RVTUDDIn_TUDD(i,j,k) Slhadata(496+i+3*(j)+9*(k))
#define RVTUDDIn_TUDDFlat(i) Slhadata(508+i)

#define OffsetRVTLLE 535
#define LengthRVTLLE 28
#define BlockRVTLLE(i) SlhaData(535+i)
#define RVTLLE_Q SlhaData(536)
#define RVTLLE_TLLE(i,j,k) Slhadata(524+i+3*(j)+9*(k))
#define RVTLLE_TLLEFlat(i) Slhadata(536+i)

#define OffsetRVTLQD 563
#define LengthRVTLQD 28
#define BlockRVTLQD(i) SlhaData(563+i)
#define RVTLQD_Q SlhaData(564)
#define RVTLQD_TLQD(i,j,k) Slhadata(552+i+3*(j)+9*(k))
#define RVTLQD_TLQDFlat(i) Slhadata(564+i)

#define OffsetRVTUDD 591
#define LengthRVTUDD 28
#define BlockRVTUDD(i) SlhaData(591+i)
#define RVTUDD_Q SlhaData(592)
#define RVTUDD_TUDD(i,j,k) Slhadata(580+i+3*(j)+9*(k))
#define RVTUDD_TUDDFlat(i) Slhadata(592+i)

#define OffsetRVKappaIn 619
#define LengthRVKappaIn 3
#define BlockRVKappaIn(i) SlhaData(619+i)
#define RVKappaIn_kappa(i) Slhadata(619+i)

#define OffsetRVKappa 622
#define LengthRVKappa 4
#define BlockRVKappa(i) SlhaData(622+i)
#define RVKappa_Q SlhaData(623)
#define RVKappa_kappa(i) Slhadata(623+i)

#define OffsetRVDIn 626
#define LengthRVDIn 3
#define BlockRVDIn(i) SlhaData(626+i)
#define RVDIn_D(i) Slhadata(626+i)

#define OffsetRVD 629
#define LengthRVD 4
#define BlockRVD(i) SlhaData(629+i)
#define RVD_Q SlhaData(630)
#define RVD_D(i) Slhadata(630+i)

#define OffsetRVSnVEVIn 633
#define LengthRVSnVEVIn 3
#define BlockRVSnVEVIn(i) SlhaData(633+i)
#define RVSnVEVIn_VEV(i) Slhadata(633+i)

#define OffsetRVSnVEV 636
#define LengthRVSnVEV 4
#define BlockRVSnVEV(i) SlhaData(636+i)
#define RVSnVEV_Q SlhaData(637)
#define RVSnVEV_VEV(i) Slhadata(637+i)

#define OffsetRVM2LH1In 640
#define LengthRVM2LH1In 3
#define BlockRVM2LH1In(i) SlhaData(640+i)
#define RVM2LH1In_M2LH1(i) Slhadata(640+i)

#define OffsetRVM2LH1 643
#define LengthRVM2LH1 4
#define BlockRVM2LH1(i) SlhaData(643+i)
#define RVM2LH1_Q SlhaData(644)
#define RVM2LH1_M2LH1(i) Slhadata(644+i)

#define OffsetRVNMix 647
#define LengthRVNMix 49
#define BlockRVNMix(i) SlhaData(647+i)
#define RVNMix_ZNeu(n1,n2) Slhadata(640+n1+7*(n2))
#define RVNMix_ZNeuFlat(i) Slhadata(647+i)

#define OffsetRVUMix 696
#define LengthRVUMix 25
#define BlockRVUMix(i) SlhaData(696+i)
#define RVUMix_UCha(c1,c2) Slhadata(691+c1+5*(c2))
#define RVUMix_UChaFlat(i) Slhadata(696+i)

#define OffsetRVVMix 721
#define LengthRVVMix 25
#define BlockRVVMix(i) SlhaData(721+i)
#define RVVMix_VCha(c1,c2) Slhadata(716+c1+5*(c2))
#define RVVMix_VChaFlat(i) Slhadata(721+i)

#define OffsetRVHMix 746
#define LengthRVHMix 25
#define BlockRVHMix(i) SlhaData(746+i)
#define RVHMix_UH(h1,h2) Slhadata(741+h1+5*(h2))
#define RVHMix_UHFlat(i) Slhadata(746+i)

#define OffsetRVAMix 771
#define LengthRVAMix 25
#define BlockRVAMix(i) SlhaData(771+i)
#define RVAMix_UA(h1,h2) Slhadata(766+h1+5*(h2))
#define RVAMix_UAFlat(i) Slhadata(771+i)

#define OffsetRVLMix 796
#define LengthRVLMix 64
#define BlockRVLMix(i) SlhaData(796+i)
#define RVLMix_CLep(l1,l2) Slhadata(788+l1+8*(l2))
#define RVLMix_CLepFlat(i) Slhadata(796+i)

#define OffsetVCKMIn 860
#define LengthVCKMIn 4
#define BlockVCKMIn(i) SlhaData(860+i)
#define VCKMIn_lambda Slhadata(861)
#define VCKMIn_A Slhadata(862)
#define VCKMIn_rhobar Slhadata(863)
#define VCKMIn_etabar Slhadata(864)

#define OffsetVCKM 864
#define LengthVCKM 10
#define BlockVCKM(i) SlhaData(864+i)
#define VCKM_Q SlhaData(865)
#define VCKM_VCKM(g1,g2) Slhadata(862+g1+3*(g2))
#define VCKM_VCKMFlat(i) Slhadata(865+i)

#define OffsetUPMNSIn 874
#define LengthUPMNSIn 6
#define BlockUPMNSIn(i) SlhaData(874+i)
#define UPMNSIn_theta12 Slhadata(875)
#define UPMNSIn_theta23 Slhadata(876)
#define UPMNSIn_theta13 Slhadata(877)
#define UPMNSIn_delta13 Slhadata(878)
#define UPMNSIn_alpha1 Slhadata(879)
#define UPMNSIn_alpha2 Slhadata(880)

#define OffsetUPMNS 880
#define LengthUPMNS 10
#define BlockUPMNS(i) SlhaData(880+i)
#define UPMNS_Q SlhaData(881)
#define UPMNS_UPMNS(g1,g2) Slhadata(878+g1+3*(g2))
#define UPMNS_UPMNSFlat(i) Slhadata(881+i)

#define OffsetMSS2In 890
#define LengthMSS2In 45
#define BlockMSS2In(i) SlhaData(890+i)
#define MSS2In_MSS2(g1,g2,q) Slhadata(878+g1+3*(g2)+9*(q))
#define MSS2In_MSS2Flat(i,q) Slhadata(881+i+9*(q))

#define OffsetMSL2In 890
#define LengthMSL2In 9
#define BlockMSL2In(i) SlhaData(890+i)
#define   MSL2In_MSL2(g1,g2) MSS2In_MSS2(g1,g2,1)
#define   MSL2In_MSL2Flat(i) MSS2In_MSS2Flat(i,1)

#define OffsetMSE2In 899
#define LengthMSE2In 9
#define BlockMSE2In(i) SlhaData(899+i)
#define   MSE2In_MSE2(g1,g2) MSS2In_MSS2(g1,g2,2)
#define   MSE2In_MSE2Flat(i) MSS2In_MSS2Flat(i,2)

#define OffsetMSQ2In 908
#define LengthMSQ2In 9
#define BlockMSQ2In(i) SlhaData(908+i)
#define   MSQ2In_MSQ2(g1,g2) MSS2In_MSS2(g1,g2,3)
#define   MSQ2In_MSQ2Flat(i) MSS2In_MSS2Flat(i,3)

#define OffsetMSU2In 917
#define LengthMSU2In 9
#define BlockMSU2In(i) SlhaData(917+i)
#define   MSU2In_MSU2(g1,g2) MSS2In_MSS2(g1,g2,4)
#define   MSU2In_MSU2Flat(i) MSS2In_MSS2Flat(i,4)

#define OffsetMSD2In 926
#define LengthMSD2In 9
#define BlockMSD2In(i) SlhaData(926+i)
#define   MSD2In_MSD2(g1,g2) MSS2In_MSS2(g1,g2,5)
#define   MSD2In_MSD2Flat(i) MSS2In_MSS2Flat(i,5)

#define OffsetMSS2 935
#define LengthMSS2 50
#define BlockMSS2(i) SlhaData(935+i)
#define MSS2_Q(q) SlhaData(926+10*(q))
#define MSS2_MSS2(g1,g2,q) Slhadata(923+g1+3*(g2)+10*(q))
#define MSS2_MSS2Flat(i,q) Slhadata(926+i+10*(q))

#define OffsetMSL2 935
#define LengthMSL2 10
#define BlockMSL2(i) SlhaData(935+i)
#define   MSL2_Q MSS2_Q(1)
#define   MSL2_MSL2(g1,g2) MSS2_MSS2(g1,g2,1)
#define   MSL2_MSL2Flat(i) MSS2_MSS2Flat(i,1)

#define OffsetMSE2 945
#define LengthMSE2 10
#define BlockMSE2(i) SlhaData(945+i)
#define   MSE2_Q MSS2_Q(2)
#define   MSE2_MSE2(g1,g2) MSS2_MSS2(g1,g2,2)
#define   MSE2_MSE2Flat(i) MSS2_MSS2Flat(i,2)

#define OffsetMSQ2 955
#define LengthMSQ2 10
#define BlockMSQ2(i) SlhaData(955+i)
#define   MSQ2_Q MSS2_Q(3)
#define   MSQ2_MSQ2(g1,g2) MSS2_MSS2(g1,g2,3)
#define   MSQ2_MSQ2Flat(i) MSS2_MSS2Flat(i,3)

#define OffsetMSU2 965
#define LengthMSU2 10
#define BlockMSU2(i) SlhaData(965+i)
#define   MSU2_Q MSS2_Q(4)
#define   MSU2_MSU2(g1,g2) MSS2_MSS2(g1,g2,4)
#define   MSU2_MSU2Flat(i) MSS2_MSS2Flat(i,4)

#define OffsetMSD2 975
#define LengthMSD2 10
#define BlockMSD2(i) SlhaData(975+i)
#define   MSD2_Q MSS2_Q(5)
#define   MSD2_MSD2(g1,g2) MSS2_MSS2(g1,g2,5)
#define   MSD2_MSD2Flat(i) MSS2_MSS2Flat(i,5)

#define OffsetTfIn 985
#define LengthTfIn 27
#define BlockTfIn(i) SlhaData(985+i)
#define TfIn_Tf(g1,g2,t) Slhadata(964+g1+3*(g2)+9*(t))
#define TfIn_TfFlat(i,t) Slhadata(967+i+9*(t))

#define OffsetTeIn 985
#define LengthTeIn 9
#define BlockTeIn(i) SlhaData(985+i)
#define   TeIn_Tf(g1,g2) TfIn_Tf(g1,g2,2)
#define   TeIn_TfFlat(i) TfIn_TfFlat(i,2)

#define OffsetTuIn 994
#define LengthTuIn 9
#define BlockTuIn(i) SlhaData(994+i)
#define   TuIn_Tf(g1,g2) TfIn_Tf(g1,g2,3)
#define   TuIn_TfFlat(i) TfIn_TfFlat(i,3)

#define OffsetTdIn 1003
#define LengthTdIn 9
#define BlockTdIn(i) SlhaData(1003+i)
#define   TdIn_Tf(g1,g2) TfIn_Tf(g1,g2,4)
#define   TdIn_TfFlat(i) TfIn_TfFlat(i,4)

#define OffsetTf 1012
#define LengthTf 30
#define BlockTf(i) SlhaData(1012+i)
#define Tf_Q(t) SlhaData(993+10*(t))
#define Tf_Tf(g1,g2,t) Slhadata(990+g1+3*(g2)+10*(t))
#define Tf_TfFlat(i,t) Slhadata(993+i+10*(t))

#define OffsetTe 1012
#define LengthTe 10
#define BlockTe(i) SlhaData(1012+i)
#define   Te_Q Tf_Q(2)
#define   Te_Tf(g1,g2) Tf_Tf(g1,g2,2)
#define   Te_TfFlat(i) Tf_TfFlat(i,2)

#define OffsetTu 1022
#define LengthTu 10
#define BlockTu(i) SlhaData(1022+i)
#define   Tu_Q Tf_Q(3)
#define   Tu_Tf(g1,g2) Tf_Tf(g1,g2,3)
#define   Tu_TfFlat(i) Tf_TfFlat(i,3)

#define OffsetTd 1032
#define LengthTd 10
#define BlockTd(i) SlhaData(1032+i)
#define   Td_Q Tf_Q(4)
#define   Td_Tf(g1,g2) Tf_Tf(g1,g2,4)
#define   Td_TfFlat(i) Tf_TfFlat(i,4)

#define OffsetASfMix 1042
#define LengthASfMix 144
#define BlockASfMix(i) SlhaData(1042+i)
#define ASfMix_UASf(s1,s2,t) Slhadata(1000+s1+6*(s2)+36*(t))
#define ASfMix_UASfFlat(i,t) Slhadata(1006+i+36*(t))

#define OffsetSnuMix 1042
#define LengthSnuMix 36
#define BlockSnuMix(i) SlhaData(1042+i)
#define   SnuMix_UASf(s1,s2) ASfMix_UASf(s1,s2,1)
#define   SnuMix_UASfFlat(i) ASfMix_UASfFlat(i,1)

#define OffsetSelMix 1078
#define LengthSelMix 36
#define BlockSelMix(i) SlhaData(1078+i)
#define   SelMix_UASf(s1,s2) ASfMix_UASf(s1,s2,2)
#define   SelMix_UASfFlat(i) ASfMix_UASfFlat(i,2)

#define OffsetUSqMix 1114
#define LengthUSqMix 36
#define BlockUSqMix(i) SlhaData(1114+i)
#define   USqMix_UASf(s1,s2) ASfMix_UASf(s1,s2,3)
#define   USqMix_UASfFlat(i) ASfMix_UASfFlat(i,3)

#define OffsetDSqMix 1150
#define LengthDSqMix 36
#define BlockDSqMix(i) SlhaData(1150+i)
#define   DSqMix_UASf(s1,s2) ASfMix_UASf(s1,s2,4)
#define   DSqMix_UASfFlat(i) ASfMix_UASfFlat(i,4)

#define OffsetSnSMix 1186
#define LengthSnSMix 9
#define BlockSnSMix(i) SlhaData(1186+i)
#define SnSMix_US(g1,g2) Slhadata(1183+g1+3*(g2))
#define SnSMix_USFlat(i) Slhadata(1186+i)

#define OffsetSnAMix 1195
#define LengthSnAMix 9
#define BlockSnAMix(i) SlhaData(1195+i)
#define SnAMix_UA(g1,g2) Slhadata(1192+g1+3*(g2))
#define SnAMix_UAFlat(i) Slhadata(1195+i)

#define OffsetCVHMix 1204
#define LengthCVHMix 16
#define BlockCVHMix(i) SlhaData(1204+i)
#define CVHMix_UH(h1,h2) Slhadata(1200+h1+4*(h2))
#define CVHMix_UHFlat(i) Slhadata(1204+i)

#define OffsetNMNMix 1220
#define LengthNMNMix 25
#define BlockNMNMix(i) SlhaData(1220+i)
#define NMNMix_ZNeu(n1,n2) Slhadata(1215+n1+5*(n2))
#define NMNMix_ZNeuFlat(i) Slhadata(1220+i)

#define OffsetNMHMix 1245
#define LengthNMHMix 9
#define BlockNMHMix(i) SlhaData(1245+i)
#define NMHMix_UH(h1,h2) Slhadata(1242+h1+3*(h2))
#define NMHMix_UHFlat(i) Slhadata(1245+i)

#define OffsetNMAMix 1254
#define LengthNMAMix 9
#define BlockNMAMix(i) SlhaData(1254+i)
#define NMAMix_UA(h1,h2) Slhadata(1251+h1+3*(h2))
#define NMAMix_UAFlat(i) Slhadata(1254+i)

#define OffsetPrecObs 1263
#define LengthPrecObs 15
#define BlockPrecObs(i) SlhaData(1263+i)
#define PrecObs_DeltaRho Slhadata(1264)
#define PrecObs_MWMSSM Slhadata(1265)
#define PrecObs_MWSM Slhadata(1266)
#define PrecObs_SW2effMSSM Slhadata(1267)
#define PrecObs_SW2effSM Slhadata(1268)
#define PrecObs_gminus2mu Slhadata(1269)
#define PrecObs_EDMeTh Slhadata(1270)
#define PrecObs_EDMn Slhadata(1271)
#define PrecObs_EDMHg Slhadata(1272)
#define PrecObs_bsgammaMSSM Slhadata(1273)
#define PrecObs_bsgammaSM Slhadata(1274)
#define PrecObs_DeltaMsMSSM Slhadata(1275)
#define PrecObs_DeltaMsSM Slhadata(1276)
#define PrecObs_BsmumuMSSM Slhadata(1277)
#define PrecObs_BsmumuSM Slhadata(1278)

#define OffsetSPInfo 1278
#define LengthSPInfo 92
#define BlockSPInfo(i) SlhaData(1278+i)
#define SPInfo_NLines SlhaData(1279)
#define SPInfo_Severity SlhaData(1280)
#define SPInfo_Code(n) SlhaData(1280+n)
#define SPInfo_Text(i,n) SlhaData(1290+i+5*(n))
#define SPInfo_TextFlat(i) SlhaData(1295+i)
#define   SPInfo_Len 80

#define OffsetDCInfo 1370
#define LengthDCInfo 92
#define BlockDCInfo(i) SlhaData(1370+i)
#define DCInfo_NLines SlhaData(1371)
#define DCInfo_Severity SlhaData(1372)
#define DCInfo_Code(n) SlhaData(1372+n)
#define DCInfo_Text(i,n) SlhaData(1382+i+5*(n))
#define DCInfo_TextFlat(i) SlhaData(1387+i)
#define   DCInfo_Len 80

#define OffsetDecays 1462
#define LengthDecays 4096
#define BlockDecays(i) SlhaData(1462+i)
#define Decays_Data(n) Slhadata(1462+n)

#define nslhadata 5558

#endif
