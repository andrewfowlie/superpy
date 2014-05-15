* gm2_2L_full8.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full8
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp52, tmp53, tmp54, tmp55, tmp56, tmp57, tmp58
	ComplexType tmp59, tmp60, tmp61, tmp62, tmp63, tmp64, tmp65
	ComplexType tmp66, tmp67, tmp68, tmp69, tmp70, tmp71, tmp72

	tmp56 = 2*MSl2Diff(1,1)**2 + MW2*(MSneu2(1) + MStau2(1))

	tmp57 = 3*MSl2Diff(1,1)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 4*MSneu2(1) - 4*MStau2(1)) + 
     &   MW2*(4*MSneu2(1) + 4*MStau2(1))

	tmp58 = MSl2Diff(1,1)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 2*MSneu2(1) - 2*MStau2(1))

	tmp59 = 2*MSl2Diff(2,1)**2 + MW2*(MSneu2(1) + MStau2(2))

	tmp60 = 3*MSl2Diff(2,1)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 4*MSneu2(1) - 4*MStau2(2)) + 
     &   MW2*(4*MSneu2(1) + 4*MStau2(2))

	tmp61 = MSl2Diff(2,1)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 2*MSneu2(1) - 2*MStau2(2))

	tmp69 = MHiggs2(4)*(-2*MW2*tmp56 + tmp57*MHiggs2(4)) + 
     &   2*MW2**2*MSl2Diff(1,1)**2

	tmp70 = MHiggs2(4)*(-2*MW2*tmp59 + tmp60*MHiggs2(4)) + 
     &   2*MW2**2*MSl2Diff(2,1)**2

	Sf2175 = (-(MW2*S2B) + ML2*TB)*tmp61*tmp69*
     &    UStau11**2 + ML*tmp61*tmp69*
     &    (reMUE + TB*Re(Atau))*UStau11*UStau12 + 
     &   (-(MW2*S2B) + ML2*TB)*tmp58*tmp70*UStau21**2 + 
     &   ML*tmp58*tmp70*(reMUE + TB*Re(Atau))*UStau21*
     &    UStau22

	Sf2177 = Sf162*UStau11 + Sf1191*UStau21

	Sf1118 = Sf1117*UStau11 - 3*ML*Sf215*Ss115*UStau12

	Sf1115 = Sf1113*UStau11 + ML*Sf1114*Sf215*UStau12

	Sf1104 = Sf1102*UStau11 + ML*Sf1103*Sf215*UStau12

	Sf1098 = Sf1097*UStau11 + ML*Sf215*Sf733*UStau12

	Sf1031 = Sf1030*UStau11 + ML*Sf215*Ss25*MHiggs2(4)*UStau12

	Sf1028 = Sf1027*UStau11 + 3*ML*Sf215*Ss241*Ss242*UStau12

	Sf1022 = Sf1021*UStau11 - 3*ML*Sf215*Ss241*Ss242*UStau12

	Sf1020 = Sf1019*UStau11 - ML*Sf215*Ss241*Ss242*UStau12

	Sf1016 = Sf1015*UStau11 + ML*Sf215*Ss241*Ss242*UStau12

	Sf1013 = Sf1012*UStau11 - ML*Sf215*Sf217*UStau12

	Sf6973 = MBot*Sf6972*USbot22 + MBot2*Sf6968*USbot21*UStop21

	Sf6952 = MBot*Sf6947*USbot22 + MBot2*Sf4224*USbot21*UStop21

	Sf4254 = MBot*Sf4253*USbot22 + MBot2*Sf4249*USbot21*UStop11

	Sf4234 = MBot*Sf4229*USbot22 + MBot2*Sf4224*USbot21*UStop11

	Sf2961 = Sf2807*USbot21**2 - Sf2225*Sf2960*USbot22

	Sf139 = Sf124*USbot21**2 + 2*Sf125*TB*USbot22**2

	Sf4800 = Sf15*Sf4795*USbot11*USbot12*USbot21 + 
     &   2*MZ2*Sf4799*SW2*USbot22

	Sf6953 = -(Sf6952*TB) + Sf6944*USbot21

	Sf6086 = MBot*Sf5693*USbot21 - 6*MBot2*Sf2225*USbot22

	Sf6081 = 2*MBot*Sf5683*USbot21 - 3*MBot2*Sf2225*USbot22

	Sf6074 = MBot*Sf5668*USbot21 - MBot2*Sf2225*USbot22

	Sf4911 = Sf4904*USbot21 - Sf2339*TB*UStop21

	Sf4235 = -(Sf4234*TB) + Sf4226*USbot21

	Sf3965 = MBot*Sf3867*USbot21 - 3*MBot2*Sf2225*USbot22

	Sf3984 = Sf3867*USbot21 - 6*MBot*Sf2225*USbot22

	Sf3959 = MBot*Sf3854*USbot21 - MBot2*Sf2225*USbot22

	Sf3981 = Sf3854*USbot21 - 2*MBot*Sf2225*USbot22

	Sf3652 = MBot*Sf3371*USbot21 - 3*MBot2*Sf2225*USbot22

	Sf3669 = Sf3371*USbot21 - 6*MBot*Sf2225*USbot22

	Sf3643 = 2*MBot*Sf3351*USbot21 - MBot2*Sf2225*USbot22

	Sf3666 = Sf3351*USbot21 - MBot*Sf2225*USbot22

	Sf2340 = Sf2325*USbot21 - Sf2339*TB*UStop11

	Sf2287 = Sf2286*USbot21 + 3*MBot*Sf2225*MHiggs2(4)*USbot22

	Sf2292 = Sf2239*USbot21 - 3*MBot*Sf2225*MHiggs2(4)*USbot22

	Sf148 = Sf147*USbot21 + Sf143*USbot22

	Sf6964 = MBot*Sf6963*USbot12 + MBot2*Sf6958*USbot11*UStop21

	Sf6948 = MBot*Sf6947*USbot12 + MBot2*Sf4224*USbot11*UStop21

	Sf4245 = MBot*Sf4244*USbot12 + MBot2*Sf4240*USbot11*UStop11

	Sf4230 = MBot*Sf4229*USbot12 + MBot2*Sf4224*USbot11*UStop11

	Sf2809 = -(Sf2225*Sf2808*USbot12) + Sf2807*USbot22**2

	Sf126 = 2*Sf125*TB*USbot12**2 + Sf124*USbot22**2

	Sf7179 = Sf4828*MSq2Diff(1,2)*USbot11 + 
     &   Sf4864*MSq2Diff(2,2)*USbot21

	Sf6383 = Sf7179

	Sf6949 = -(Sf6948*TB) + Sf6944*USbot11

	Sf5694 = MBot*Sf5693*USbot11 - 6*MBot2*Sf2225*USbot12

	Sf5684 = 2*MBot*Sf5683*USbot11 - 3*MBot2*Sf2225*USbot12

	Sf5669 = MBot*Sf5668*USbot11 - MBot2*Sf2225*USbot12

	Sf4905 = Sf4904*USbot11 - Sf2326*TB*UStop21

	Sf6349 = Sf4828*USbot11 + Sf4864*USbot21

	Sf4797 = (Sf25*Sf4795*USbot11 - MZ2*Sf4796*USbot21)/MZ2

	Sf4231 = -(Sf4230*TB) + Sf4226*USbot11

	Sf3868 = MBot*Sf3867*USbot11 - 3*MBot2*Sf2225*USbot12

	Sf3887 = Sf3867*USbot11 - 6*MBot*Sf2225*USbot12

	Sf3855 = MBot*Sf3854*USbot11 - MBot2*Sf2225*USbot12

	Sf3884 = Sf3854*USbot11 - 2*MBot*Sf2225*USbot12

	Sf3372 = MBot*Sf3371*USbot11 - 3*MBot2*Sf2225*USbot12

	Sf3389 = Sf3371*USbot11 - 6*MBot*Sf2225*USbot12

	Sf3352 = 2*MBot*Sf3351*USbot11 - MBot2*Sf2225*USbot12

	Sf3386 = Sf3351*USbot11 - MBot*Sf2225*USbot12

	Sf2327 = Sf2325*USbot11 - Sf2326*TB*UStop11

	Sf2240 = Sf2239*USbot11 - 3*MBot*Sf2225*MHiggs2(4)*USbot12

	Sf2226 = Sf2224*USbot11 + 3*MBot*Sf2225*MHiggs2(4)*USbot12

	Sf4208 = Sf2198*USbot11 + Sf2263*USbot21

	Sf137 = Sf136*USbot11 + Sf132*USbot12

	Sf8705 = Sf8520**2*Sf8704 + Sf8532**2*Sf8699*MStop2(3)**2

	Sf8668 = Sf8520**2*Sf8667 + Sf7990**2*Sf8660*MStop2(3)**2

	Sf8703 = MZ2*Sf8702 + 64*CW2*Sf8532*Ss78*SW2*MStop2(2)*MStop2(3)

	Sf8666 = -(MZ2*Sf8665) + 
     &   64*CW2*Sf7990*Ss78*SW2*MStop2(1)*MStop2(3)

	Sf8706 = 3*Sf7497*Sf8705 + Sf8703*MStop2(3)

	Sf8669 = 3*Sf7497*Sf8668 + Sf8666*MStop2(3)

	Sf10280 = MZ2**2*Sf10279 + Sf9113*MStop2(2)**2

	Sf10510 = MZ2**2*Sf10509 + Sf8911*MStop2(2)**2

	Sf10263 = MZ2**2*Sf10262 + Sf9113*MStop2(1)**2

	Sf10500 = MZ2**2*Sf10499 + Sf8911*MStop2(1)**2

	Sf10495 = MStop2(1)**2 - Sf10443*MStop2(2)

	Sf10386 = MStop2(1)**2 - Sf10295*MStop2(2)

	Sf111 = (1/2.D0*Sf106**2/(MStop(1)*MStop(2)) - 
     &     2*MStop(1)*MStop(2))/Sf106

	Sf103 = (1/2.D0*Sf98**2/(MSbot(2)*MStop(2)) - 
     &     2*MSbot(2)*MStop(2))/Sf98

	Sf95 = (1/2.D0*Sf90**2/(MSbot(1)*MStop(2)) - 
     &     2*MSbot(1)*MStop(2))/Sf90

	Sf81 = (1/2.D0*Sf76**2/(MSbot(2)*MStop(1)) - 
     &     2*MSbot(2)*MStop(1))/Sf76

	Sf73 = (1/2.D0*Sf68**2/(MSbot(1)*MStop(1)) - 
     &     2*MSbot(1)*MStop(1))/Sf68

	Sf7752 = Sf7701**2*Sf7751 + Sf7713**2*Sf7750*MStau2(3)**2

	Sf7730 = Sf7701**2*Sf7729 + Sf7579**2*Sf7728*MStau2(3)**2

	Sf7753 = 3*Sf7497*Sf7752*MStau2(2) + Sf7749*MStau2(3)

	Sf10141 = MZ2**2*Sf10140 + Sf9113*MStau2(2)**2

	Sf9385 = MZ2**2*Sf9384 + Sf8911*MStau2(2)**2

	Sf1367 = Sf1366*MHiggs2(4) + 2*Sf1365*MStau2(2)

	Sf1229 = 2*Sf1228*MStau2(2) + Sf236*MSneu2(1)*UStau21

	Sf9974 = 4*Sf9332**3*Sf9333**3*Sf9964*MHiggs2(2)*
     &    MStau2(1)**3 + Sf9973*MStau2(2)

	Sf9916 = 16*Sf9332**3*Sf9333**3*Sf9906*MHiggs2(2)*
     &    MStau2(1)**3 - Sf9915*MStau2(2)

	Sf9900 = 2*Sf9332**3*Sf9333**3*Sf9394*MHiggs2(2)*
     &    MStau2(1)**3 - Sf9899*MStau2(2)

	Sf9429 = 32*Sf9332**3*Sf9333**3*MHiggs2(2)*
     &    MStau2(1)**3 - Sf9428*MStau2(2)

	Sf9542 = 8*Sf9332**3*Sf9333**3*MHiggs2(2)*MStau2(1)**3 - 
     &   Sf9541*MStau2(2)

	Sf9926 = 16*MZ2*Sf7728*Sf9290**3*Sf9291**3*
     &    MStau2(1)**3 - Sf9925*MStau2(2)

	Sf9963 = MZ2**4*(8*MStau2(2)**2 + 
     &      MStau2(1)*(29*MStau2(1) + MStau2(2))) - 
     &   2*MZ2**3*(34*MStau2(2)**3 + 
     &      MStau2(1)*(-43*MStau2(2)**2 + 
     &         MStau2(1)*(33*MStau2(1) + 128*MStau2(2)))) + 
     &   MZ2**2*(16*MStau2(2)**4 + 
     &      MStau2(1)*(123*MStau2(2)**3 + 
     &         MStau2(1)*(51*MStau2(2)**2 + 
     &            MStau2(1)*(73*MStau2(1) + 345*MStau2(2))))) - 
     &   2*MZ2*(-22*MStau2(2)**3 + 
     &      MStau2(1)*(17*MStau2(2)**2 + 
     &         MStau2(1)*(32*MStau2(1) + 13*MStau2(2))))*
     &    MStau2(3)**2 + 4*MStau2(1)*(5*MStau2(1) - 2*MStau2(3))*
     &    MStau2(3)**4

	Sf9890 = MZ2**4*(MStau2(1)*(5*MStau2(1) - 5*MStau2(2)) + 
     &      2*MStau2(2)**2) - 
     &   MZ2**3*(17*MStau2(2)**3 + 
     &      MStau2(1)*(-41*MStau2(2)**2 + 
     &         MStau2(1)*(9*MStau2(1) + 31*MStau2(2)))) + 
     &   MZ2**2*(22*MStau2(2)**4 + 
     &      MStau2(1)*(-45*MStau2(2)**3 + 
     &         MStau2(1)*(15*MStau2(2)**2 + 
     &            MStau2(1)*(7*MStau2(1) + 33*MStau2(2))))) - 
     &   MZ2*(7*MStau2(2)**3 + 
     &      MStau2(1)*(-41*MStau2(2)**2 + 
     &         MStau2(1)*(7*MStau2(1) + 11*MStau2(2))))*
     &    MStau2(3)**2 + 4*MStau2(1)*(-MStau2(2) - MStau2(3))*
     &    MStau2(3)**4

	Sf9397 = Sf9329*Sf9332**3*Sf9333**3*MHiggs2(2)*
     &    MStau2(1)**2 + Sf9396*MStau2(2)

	Sf9677 = 16*MZ2**2*Sf9290**3*Sf9291**3*MStau2(1)**2 + 
     &   Sf9676*MStau2(2)

	Sf9413 = MZ2**2*Sf9290**3*Sf9291**3*MStau2(1)**2 + 
     &   Sf9412*MStau2(2)

	Sf10132 = MZ2**2*Sf10131 + Sf9113*MStau2(1)**2

	Sf9377 = MZ2**2*Sf9376 + Sf8911*MStau2(1)**2

	Sf9372 = MStau2(1)**2 - Sf9330*MStau2(2)

	Sf10119 = MStau2(1)**2 - Sf10067*MStau2(2)

	Sf9667 = Sf9656*MStau2(1) + Sf9666*MStau2(2)

	Sf9637 = Sf9626*MStau2(1) + Sf9636*MStau2(2)

	Sf9588 = Sf9555*MStau2(1) + Sf9587*MStau2(2)

	Sf9404 = Sf9399*MStau2(1) + Sf9403*MStau2(2)

	Sf7731 = 3*Sf7497*Sf7730*MStau2(1) - Sf7727*MStau2(3)

	Sf450 = Sf449*MHiggs2(4) + 2*Sf448*MStau2(1)

	Sf247 = 2*Sf246*MStau2(1) + Sf236*MSneu2(1)*UStau11

	Sf59 = (1/2.D0*Sf54**2/(MStau(1)*MStau(2)) - 
     &     2*MStau(1)*MStau(2))/Sf54

	Sf51 = (1/2.D0*Sf46**2/(MSneu(1)*MStau(2)) - 
     &     2*MSneu(1)*MStau(2))/Sf46

	Sf38 = (1/2.D0*Sf33**2/(MSneu(1)*MStau(1)) - 
     &     2*MSneu(1)*MStau(1))/Sf33

	Sf5145 = Sf5144*MHiggs2(4) + Sf4864*TB*MSq2Diff(2,2)**2

	Sf6263 = MW2*Sf6262*MSq2Diff(2,2) - 
     &   2*Sf6144*Ss241*Ss242*Ss25*MStop2(2)

	Sf2737 = Sf2736*MHiggs2(4) + Sf2263*TB*MSq2Diff(2,1)**2

	Sf6260 = MW2*Sf6259*MSq2Diff(2,1) - 
     &   2*Sf3742*Ss241*Ss242*Ss25*MStop2(1)

	Sf5026 = Sf5025*MHiggs2(4) + Sf4828*TB*MSq2Diff(1,2)**2

	Sf5899 = MW2*Sf5898*MSq2Diff(1,2) - 
     &   2*Sf5761*Ss241*Ss242*Ss25*MStop2(2)

	Sf2559 = Sf2558*MHiggs2(4) + Sf2198*TB*MSq2Diff(1,1)**2

	Sf5896 = MW2*Sf5895*MSq2Diff(1,1) - 
     &   2*Sf3733*Ss241*Ss242*Ss25*MStop2(1)

	Sf1714 = 3*Sf920*MSneu2(1)**2 - Sf1713*MStau2(2)

	Sf924 = 3*Sf920*MSneu2(1)**2 - Sf923*MStau2(1)

	Sf1712 = Sf914*MSneu2(1)**2 - Sf1711*MStau2(2)

	Sf918 = Sf914*MSneu2(1)**2 - Sf917*MStau2(1)

	Sf1693 = Sf890*MSneu2(1)**2 - Sf1692*MStau2(2)

	Sf894 = Sf890*MSneu2(1)**2 - Sf893*MStau2(1)

	Sf1691 = Sf886*MSneu2(1)**2 + Sf1690*MStau2(2)

	Sf889 = Sf886*MSneu2(1)**2 + Sf888*MStau2(1)

	Sf1686 = 3*Sf876*MSneu2(1)**2 - Sf1685*MStau2(2)

	Sf880 = 3*Sf876*MSneu2(1)**2 - Sf879*MStau2(1)

	Sf1684 = Sf872*MSneu2(1)**2 + Sf1683*MStau2(2)

	Sf875 = Sf872*MSneu2(1)**2 + Sf874*MStau2(1)

	Sf1649 = Sf825*MSneu2(1)**2 + 3*Sf1648*MStau2(2)

	Sf828 = Sf825*MSneu2(1)**2 + 3*Sf827*MStau2(1)

	Sf1618 = 6*Sf784*MSneu2(1)**2 + Sf1617*MStau2(2)

	Sf788 = 6*Sf784*MSneu2(1)**2 + Sf787*MStau2(1)

	Sf1614 = 2*Sf773*MSneu2(1)**2 + Sf1613*MStau2(2)

	Sf777 = 2*Sf773*MSneu2(1)**2 + Sf776*MStau2(1)

	Sf1626 = 28*Sf759*MSneu2(1)**2 + Sf1625*MStau2(2)

	Sf801 = 28*Sf759*MSneu2(1)**2 + Sf800*MStau2(1)

	Sf1997 = Sf650*MSneu2(1)**2 + Sf1996*MStau2(2)

	Sf1939 = Sf650*MSneu2(1)**2 + Sf1938*MStau2(1)

	Sf1439 = Sf531*MSneu2(1)**2 + Sf1438*MStau2(2)

	Sf534 = Sf531*MSneu2(1)**2 + Sf533*MStau2(1)

	Sf1421 = Sf508*MSneu2(1)**2 + Sf1420*MStau2(2)

	Sf510 = Sf508*MSneu2(1)**2 + Sf509*MStau2(1)

	Sf1419 = Sf503*MSneu2(1)**2 - Sf1418*MStau2(2)

	Sf507 = Sf503*MSneu2(1)**2 - Sf506*MStau2(1)

	Sf1350 = 4*Sf420*MSneu2(1)**2 + Sf1349*MStau2(2)

	Sf424 = 4*Sf420*MSneu2(1)**2 + Sf423*MStau2(1)

	Sf1348 = 4*Sf415*MSneu2(1)**2 + Sf1347*MStau2(2)

	Sf419 = 4*Sf415*MSneu2(1)**2 + Sf418*MStau2(1)

	Sf1322 = 3*Sf382*MSneu2(1)**2 + Sf1321*MStau2(2)

	Sf385 = 3*Sf382*MSneu2(1)**2 + Sf384*MStau2(1)

	Sf1320 = Sf376*MSneu2(1)**2 + Sf1319*MStau2(2)

	Sf380 = Sf376*MSneu2(1)**2 + Sf379*MStau2(1)

	Sf1318 = Sf371*MSneu2(1)**2 + Sf1317*MStau2(2)

	Sf375 = Sf371*MSneu2(1)**2 + Sf374*MStau2(1)

	Sf1290 = Sf334*MSneu2(1)**2 + Sf1289*MStau2(2)

	Sf338 = Sf334*MSneu2(1)**2 + Sf337*MStau2(1)

	Sf1288 = Sf329*MSneu2(1)**2 + Sf1287*MStau2(2)

	Sf332 = Sf329*MSneu2(1)**2 + Sf331*MStau2(1)

	Sf2158 = 3*Sf2085*MSneu2(1)**2 + Sf2157*MStau2(2)

	Sf2089 = 3*Sf2085*MSneu2(1)**2 + Sf2088*MStau2(1)

	Sf2156 = Sf2080*MSneu2(1)**2 + Sf2155*MStau2(2)

	Sf2083 = Sf2080*MSneu2(1)**2 + Sf2082*MStau2(1)

	Sf2154 = 11*Sf2076*MSneu2(1)**2 + Sf2153*MStau2(2)

	Sf2079 = 11*Sf2076*MSneu2(1)**2 + Sf2078*MStau2(1)

	Sf1990 = Sf2*MSneu2(1)**2 + Sf1989*MStau2(2)

	Sf1932 = Sf2*MSneu2(1)**2 + Sf1931*MStau2(1)

	Sf1992 = 3*Sf118*MSneu2(1)**2 + Sf1991*MStau2(2)

	Sf1676 = 3*Sf118*MSneu2(1)**2 - Sf1675*MStau2(2)

	Sf863 = 3*Sf118*MSneu2(1)**2 - Sf862*MStau2(1)

	Sf1934 = 3*Sf118*MSneu2(1)**2 + Sf1933*MStau2(1)

	Sf1857 = Sf1080*MSneu2(1)**2 + Sf1856*MStau2(2)

	Sf1084 = Sf1080*MSneu2(1)**2 + Sf1083*MStau2(1)

	Sf1855 = Sf1075*MSneu2(1)**2 + Sf1854*MStau2(2)

	Sf1078 = Sf1075*MSneu2(1)**2 + Sf1077*MStau2(1)

	Sf1839 = Sf1052*MSneu2(1)**2 + Sf1838*MStau2(2)

	Sf1056 = Sf1052*MSneu2(1)**2 + Sf1055*MStau2(1)

	Sf1995 = 3*Sf1050*MSneu2(1)**2 + Sf1838*MStau2(2)

	Sf1937 = 3*Sf1050*MSneu2(1)**2 + Sf1055*MStau2(1)

	Sf8468 = Sf7698*MSneu2(1) - Sf8467*MStau2(2)

	Sf8471 = Sf7684*MSneu2(1) - Sf8470*MStau2(1)

	Sf762 = 3*Sf755*MSneu2(1) + Sf758*MStau2(1)

	Sf1607 = 3*Sf1604*MSneu2(1) + Sf1605*MStau2(2)

	Sf1661 = Sf1298*Sf1654*MSl2Diff(2,1) + 
     &   Sf1660*MSl2Diff(2,1)**3*UStau21**2

	Sf1301 = 3*Sf1296*Sf1298*MSl2Diff(2,1) + 
     &   Sf1278*MSl2Diff(2,1)**3*UStau21**2

	Sf844 = Sf347*Sf837*MSl2Diff(1,1) + 
     &   Sf843*MSl2Diff(1,1)**3*UStau22**2

	Sf350 = 3*Sf344*Sf347*MSl2Diff(1,1) + 
     &   Sf318*MSl2Diff(1,1)**3*UStau22**2

	Sf8408 = Sf7509**2*Sf8407 + Sf7521**2*Sf8354*MSbot2(3)**2

	Sf8312 = Sf7509**2*Sf8311 + Sf7501**2*Sf8257*MSbot2(3)**2

	Sf8406 = MZ2*Sf8405 + 32*CW2*Sf7521*Ss78*SW2*MSbot2(2)*MSbot2(3)

	Sf8310 = -(MZ2*Sf8309) + 
     &   32*CW2*Sf7501*Ss78*SW2*MSbot2(1)*MSbot2(3)

	Sf8409 = 3*Sf7497*Sf8408 + Sf8406*MSbot2(3)

	Sf8313 = 3*Sf7497*Sf8312 + Sf8310*MSbot2(3)

	Sf9206 = MZ2**2*Sf9205 + Sf9113*MSbot2(2)**2

	Sf9008 = MZ2**2*Sf9007 + Sf8911*MSbot2(2)**2

	Sf6147 = 3*Sf784*MSbot2(2)**2 + Sf6146*MStop2(2)

	Sf3899 = 3*Sf784*MSbot2(2)**2 + Sf3898*MStop2(1)

	Sf6813 = Sf6811*MSbot2(2)**2 - Sf6812*MStop2(2)**2

	Sf6867 = Sf6495*MSbot2(2)**2 + Sf6866*MStop2(2)

	Sf6753 = Sf6495*MSbot2(2)**2 + Sf6752*MStop2(1)

	Sf6865 = 3*Sf6489*MSbot2(2)**2 + Sf6864*MStop2(2)

	Sf6751 = 3*Sf6489*MSbot2(2)**2 + Sf6750*MStop2(1)

	Sf6863 = 2*Sf6484*MSbot2(2)**2 + Sf6862*MStop2(2)

	Sf6749 = 2*Sf6484*MSbot2(2)**2 + Sf6748*MStop2(1)

	Sf6033 = Sf6031*MSbot2(2)**2 + Sf6032*MStop2(2)**2

	Sf6056 = 3*Sf569*MSbot2(2)**2 + Sf6055*MStop2(2)

	Sf3725 = 3*Sf569*MSbot2(2)**2 + Sf3724*MStop2(1)

	Sf5082 = Sf5080*MSbot2(2)**2 - Sf5081*MStop2(2)**2

	Sf5070 = Sf5068*MSbot2(2)**2 - Sf5069*MStop2(2)**2

	Sf6896 = 5*Sf4695*MSbot2(2)**2 - Sf6895*MStop2(2)

	Sf4700 = 5*Sf4695*MSbot2(2)**2 - Sf4699*MStop2(1)

	Sf7092 = Sf4549*MSbot2(2)**2 + Sf7091*MStop2(2)

	Sf4661 = Sf4549*MSbot2(2)**2 + Sf4660*MStop2(1)

	Sf4474 = Sf4472*MSbot2(2)**2 - Sf4473*MStop2(1)**2

	Sf6158 = 3*Sf3778*MSbot2(2)**2 + Sf6157*MStop2(2)

	Sf3910 = 3*Sf3778*MSbot2(2)**2 + Sf3909*MStop2(1)

	Sf6156 = 3*Sf3773*MSbot2(2)**2 + Sf6155*MStop2(2)

	Sf3908 = 3*Sf3773*MSbot2(2)**2 + Sf3907*MStop2(1)

	Sf6149 = Sf3755*MSbot2(2)**2 + Sf6148*MStop2(2)

	Sf3901 = Sf3755*MSbot2(2)**2 + Sf3900*MStop2(1)

	Sf3692 = Sf3690*MSbot2(2)**2 + Sf3691*MStop2(1)**2

	Sf6219 = Sf3423*MSbot2(2)**2 + Sf6218*MStop2(2)

	Sf4160 = Sf3423*MSbot2(2)**2 + Sf4159*MStop2(1)

	Sf6003 = Sf3270*MSbot2(2)**2 + Sf6002*MStop2(2)

	Sf3600 = Sf3270*MSbot2(2)**2 + Sf3599*MStop2(1)

	Sf5975 = 6*Sf3191*MSbot2(2)**2 - Sf5974*MStop2(2)

	Sf3550 = 6*Sf3191*MSbot2(2)**2 - Sf3549*MStop2(1)

	Sf5973 = 2*Sf3186*MSbot2(2)**2 - Sf5972*MStop2(2)

	Sf3548 = 2*Sf3186*MSbot2(2)**2 - Sf3547*MStop2(1)

	Sf5947 = Sf3150*MSbot2(2)**2 + Sf5946*MStop2(2)

	Sf3515 = Sf3150*MSbot2(2)**2 + Sf3514*MStop2(1)

	Sf5917 = 6*Sf299*MSbot2(2)**2 - Sf5916*MStop2(2)

	Sf3473 = 6*Sf299*MSbot2(2)**2 - Sf3472*MStop2(1)

	Sf5392 = Sf2839*MSbot2(2)**2 + Sf5391*MStop2(2)

	Sf2984 = Sf2839*MSbot2(2)**2 + Sf2983*MStop2(1)

	Sf5390 = 2*Sf2834*MSbot2(2)**2 + Sf5389*MStop2(2)

	Sf2982 = 2*Sf2834*MSbot2(2)**2 + Sf2981*MStop2(1)

	Sf5359 = 3*Sf2792*MSbot2(2)**2 + Sf5358*MStop2(2)

	Sf2948 = 3*Sf2792*MSbot2(2)**2 + Sf2947*MStop2(1)

	Sf5357 = 2*Sf2786*MSbot2(2)**2 + Sf5356*MStop2(2)

	Sf2946 = 2*Sf2786*MSbot2(2)**2 + Sf2945*MStop2(1)

	Sf5355 = Sf2782*MSbot2(2)**2 + Sf5354*MStop2(2)

	Sf2944 = Sf2782*MSbot2(2)**2 + Sf2943*MStop2(1)

	Sf6021 = 6*Sf2722**5*MSbot2(2)**2 + Sf6020*MStop2(2)**3

	Sf3675 = 6*Sf2722**5*MSbot2(2)**2 + Sf3674*MStop2(1)**3

	Sf2642 = Sf2640*MSbot2(2)**2 - Sf2641*MStop2(1)**2

	Sf5078 = Sf2633*MSbot2(2)**2 - Sf5077*MStop2(2)

	Sf2637 = Sf2633*MSbot2(2)**2 - Sf2636*MStop2(1)

	Sf2621 = Sf2619*MSbot2(2)**2 - Sf2620*MStop2(1)**2

	Sf5104 = 2*Sf2466*MSbot2(2)**2 - 3*Sf5103*MStop2(2)

	Sf2671 = 2*Sf2466*MSbot2(2)**2 - 3*Sf2670*MStop2(1)

	Sf7094 = 4*Sf207*MSbot2(2)**2 + Sf7093*MStop2(2)

	Sf4663 = 4*Sf207*MSbot2(2)**2 + Sf4662*MStop2(1)

	Sf5919 = 6*Sf174*MSbot2(2)**2 - Sf5918*MStop2(2)

	Sf3475 = 6*Sf174*MSbot2(2)**2 - Sf3474*MStop2(1)

	Sf7090 = 16*Sf1095*MSbot2(2)**2 + Sf7089*MStop2(2)

	Sf4659 = 16*Sf1095*MSbot2(2)**2 + Sf4658*MStop2(1)

	Sf6330 = 2*Sf1060*MSbot2(2)**2 + Sf6329*MStop2(2)

	Sf4179 = 2*Sf1060*MSbot2(2)**2 + Sf4178*MStop2(1)

	Sf6790 = 210*MSbot2(2)**2 + Sf6789*MStop2(2)

	Sf4435 = 210*MSbot2(2)**2 + Sf4434*MStop2(1)

	Sf7405 = MSbot2(1)*MSbot2(2) - Sf7404*MStop2(1)

	Sf6841 = 3*Sf6720*MSbot2(2) + Sf6840*MStop2(2)

	Sf6723 = 3*Sf6720*MSbot2(2) + Sf6722*MStop2(1)

	Sf6839 = Sf6716*MSbot2(2) + Sf6838*MStop2(2)

	Sf6719 = Sf6716*MSbot2(2) + Sf6718*MStop2(1)

	Sf6834 = 3*Sf6710*MSbot2(2) + Sf6833*MStop2(2)

	Sf6713 = 3*Sf6710*MSbot2(2) + Sf6712*MStop2(1)

	Sf6832 = Sf6706*MSbot2(2) + Sf6831*MStop2(2)

	Sf6709 = Sf6706*MSbot2(2) + Sf6708*MStop2(1)

	Sf6828 = 6*Sf6701*MSbot2(2) + Sf6827*MStop2(2)

	Sf6704 = 6*Sf6701*MSbot2(2) + Sf6703*MStop2(1)

	Sf6826 = 6*Sf6697*MSbot2(2) + Sf6825*MStop2(2)

	Sf6700 = 6*Sf6697*MSbot2(2) + Sf6699*MStop2(1)

	Sf6242 = 3*Sf6128*MSbot2(2) + Sf6241*MStop2(2)

	Sf6131 = 3*Sf6128*MSbot2(2) + Sf6130*MStop2(1)

	Sf6893 = 3*Sf4687*MSbot2(2) + Sf6892*MStop2(2)

	Sf4690 = 3*Sf4687*MSbot2(2) + Sf4689*MStop2(1)

	Sf7028 = Sf4512*MSbot2(2) + Sf7027*MStop2(2)

	Sf4515 = Sf4512*MSbot2(2) + Sf4514*MStop2(1)

	Sf7022 = Sf4502*MSbot2(2) + Sf7021*MStop2(2)

	Sf4505 = Sf4502*MSbot2(2) + Sf4504*MStop2(1)

	Sf7020 = 2*Sf4498*MSbot2(2) + 3*Sf7019*MStop2(2)

	Sf4501 = 2*Sf4498*MSbot2(2) + 3*Sf4500*MStop2(1)

	Sf7017 = Sf4492*MSbot2(2) + Sf7016*MStop2(2)

	Sf4495 = Sf4492*MSbot2(2) + Sf4494*MStop2(1)

	Sf7015 = 3*Sf4488*MSbot2(2) + Sf7014*MStop2(2)

	Sf4491 = 3*Sf4488*MSbot2(2) + Sf4490*MStop2(1)

	Sf6807 = Sf4464*MSbot2(2) - Sf6806*MStop2(2)

	Sf4467 = Sf4464*MSbot2(2) - Sf4466*MStop2(1)

	Sf6805 = Sf4459*MSbot2(2) + Sf6804*MStop2(2)

	Sf4462 = Sf4459*MSbot2(2) + Sf4461*MStop2(1)

	Sf6798 = 3*Sf4446*MSbot2(2) + Sf6797*MStop2(2)

	Sf4449 = 3*Sf4446*MSbot2(2) + Sf4448*MStop2(1)

	Sf6796 = Sf4442*MSbot2(2) + Sf6795*MStop2(2)

	Sf4445 = Sf4442*MSbot2(2) + Sf4444*MStop2(1)

	Sf6793 = 3*Sf4437*MSbot2(2) - Sf6792*MStop2(2)

	Sf4440 = 3*Sf4437*MSbot2(2) - Sf4439*MStop2(1)

	Sf6341 = 3*Sf4193*MSbot2(2) + Sf6340*MStop2(2)

	Sf4196 = 3*Sf4193*MSbot2(2) + Sf4195*MStop2(1)

	Sf6339 = Sf4189*MSbot2(2) - Sf6338*MStop2(2)

	Sf4192 = Sf4189*MSbot2(2) - Sf4191*MStop2(1)

	Sf6244 = 3*Sf4181*MSbot2(2) - Sf6243*MStop2(2)

	Sf6134 = 3*Sf4181*MSbot2(2) - Sf6133*MStop2(1)

	Sf6333 = Sf4181*MSbot2(2) + 3*Sf6332*MStop2(2)

	Sf4183 = Sf4181*MSbot2(2) + 3*Sf4182*MStop2(1)

	Sf6328 = Sf4173*MSbot2(2) + Sf6327*MStop2(2)

	Sf4176 = Sf4173*MSbot2(2) + Sf4175*MStop2(1)

	Sf6205 = Sf4138*MSbot2(2) - Sf6204*MStop2(2)

	Sf4141 = Sf4138*MSbot2(2) - Sf4140*MStop2(1)

	Sf6200 = 2*Sf4131*MSbot2(2) - Sf6199*MStop2(2)

	Sf4134 = 2*Sf4131*MSbot2(2) - Sf4133*MStop2(1)

	Sf6044 = Sf3705*MSbot2(2) + Sf6043*MStop2(2)

	Sf3708 = Sf3705*MSbot2(2) + Sf3707*MStop2(1)

	Sf6042 = Sf3701*MSbot2(2) + Sf6041*MStop2(2)

	Sf3704 = Sf3701*MSbot2(2) + Sf3703*MStop2(1)

	Sf7026 = Sf3685*MSbot2(2) + Sf7025*MStop2(2)

	Sf6233 = Sf3685*MSbot2(2) + Sf6232*MStop2(2)

	Sf6030 = Sf3685*MSbot2(2) + Sf6029*MStop2(2)

	Sf6114 = Sf3685*MSbot2(2) + Sf6113*MStop2(1)

	Sf4511 = Sf3685*MSbot2(2) + Sf4510*MStop2(1)

	Sf3688 = Sf3685*MSbot2(2) + Sf3687*MStop2(1)

	Sf6028 = 3*Sf3681*MSbot2(2) + Sf6027*MStop2(2)

	Sf3684 = 3*Sf3681*MSbot2(2) + Sf3683*MStop2(1)

	Sf6001 = 3*Sf3595*MSbot2(2) + Sf6000*MStop2(2)

	Sf3598 = 3*Sf3595*MSbot2(2) + Sf3597*MStop2(1)

	Sf5994 = Sf3586*MSbot2(2) + Sf5993*MStop2(2)

	Sf3589 = Sf3586*MSbot2(2) + Sf3588*MStop2(1)

	Sf5992 = Sf3581*MSbot2(2) + Sf5991*MStop2(2)

	Sf3584 = Sf3581*MSbot2(2) + Sf3583*MStop2(1)

	Sf5945 = 3*Sf3510*MSbot2(2) + Sf5944*MStop2(2)

	Sf3513 = 3*Sf3510*MSbot2(2) + Sf3512*MStop2(1)

	Sf5943 = Sf3506*MSbot2(2) + Sf5942*MStop2(2)

	Sf3509 = Sf3506*MSbot2(2) + Sf3508*MStop2(1)

	Sf5938 = Sf3499*MSbot2(2) + Sf5937*MStop2(2)

	Sf3502 = Sf3499*MSbot2(2) + Sf3501*MStop2(1)

	Sf5413 = 3*Sf3004*MSbot2(2) + Sf5412*MStop2(2)

	Sf3007 = 3*Sf3004*MSbot2(2) + Sf3006*MStop2(1)

	Sf5132 = Sf2716*MSbot2(2) + Sf5131*MStop2(2)

	Sf2719 = Sf2716*MSbot2(2) + Sf2718*MStop2(1)

	Sf5130 = 3*Sf2712*MSbot2(2) + Sf5129*MStop2(2)

	Sf2715 = 3*Sf2712*MSbot2(2) + Sf2714*MStop2(1)

	Sf5126 = 3*Sf2706*MSbot2(2) + Sf5125*MStop2(2)

	Sf2708 = 3*Sf2706*MSbot2(2) + Sf2707*MStop2(1)

	Sf5124 = Sf2702*MSbot2(2) + Sf5123*MStop2(2)

	Sf2705 = Sf2702*MSbot2(2) + Sf2704*MStop2(1)

	Sf5121 = 6*Sf2697*MSbot2(2) + Sf5120*MStop2(2)

	Sf2700 = 6*Sf2697*MSbot2(2) + Sf2699*MStop2(1)

	Sf5119 = 6*Sf2693*MSbot2(2) + Sf5118*MStop2(2)

	Sf2696 = 6*Sf2693*MSbot2(2) + Sf2695*MStop2(1)

	Sf5110 = Sf2678*MSbot2(2) + Sf5109*MStop2(2)

	Sf2681 = Sf2678*MSbot2(2) + Sf2680*MStop2(1)

	Sf5108 = Sf2674*MSbot2(2) + Sf5107*MStop2(2)

	Sf2677 = Sf2674*MSbot2(2) + Sf2676*MStop2(1)

	Sf5102 = Sf2665*MSbot2(2) + Sf5101*MStop2(2)

	Sf2668 = Sf2665*MSbot2(2) + Sf2667*MStop2(1)

	Sf5098 = 3*Sf2659*MSbot2(2) + Sf5097*MStop2(2)

	Sf2662 = 3*Sf2659*MSbot2(2) + Sf2661*MStop2(1)

	Sf5096 = Sf2655*MSbot2(2) + Sf5095*MStop2(2)

	Sf2658 = Sf2655*MSbot2(2) + Sf2657*MStop2(1)

	Sf5076 = Sf2628*MSbot2(2) - Sf5075*MStop2(2)

	Sf2631 = Sf2628*MSbot2(2) - Sf2630*MStop2(1)

	Sf5074 = 3*Sf2624*MSbot2(2) + Sf5073*MStop2(2)

	Sf2627 = 3*Sf2624*MSbot2(2) + Sf2626*MStop2(1)

	Sf5046 = Sf2589*MSbot2(2) + Sf5045*MStop2(2)

	Sf2592 = Sf2589*MSbot2(2) + Sf2591*MStop2(1)

	Sf9197 = MZ2**2*Sf9196 + Sf9113*MSbot2(1)**2

	Sf9000 = MZ2**2*Sf8999 + Sf8911*MSbot2(1)**2

	Sf5764 = 3*Sf784*MSbot2(1)**2 + Sf5763*MStop2(2)

	Sf3754 = 3*Sf784*MSbot2(1)**2 + Sf3753*MStop2(1)

	Sf6610 = Sf6608*MSbot2(1)**2 + Sf6609*MStop2(2)

	Sf6558 = Sf6556*MSbot2(1)**2 - Sf6557*MStop2(2)**2

	Sf6615 = Sf6495*MSbot2(1)**2 + Sf6614*MStop2(2)

	Sf6498 = Sf6495*MSbot2(1)**2 + Sf6497*MStop2(1)

	Sf6612 = 3*Sf6489*MSbot2(1)**2 + Sf6611*MStop2(2)

	Sf6493 = 3*Sf6489*MSbot2(1)**2 + Sf6492*MStop2(1)

	Sf6488 = 2*Sf6484*MSbot2(1)**2 + Sf6487*MStop2(1)

	Sf5772 = Sf5769*MSbot2(1)**2 + Sf5771*MStop2(2)

	Sf5462 = 3*Sf569*MSbot2(1)**2 + Sf5461*MStop2(2)

	Sf3454 = 3*Sf569*MSbot2(1)**2 + Sf3453*MStop2(1)

	Sf5613 = Sf5611*MSbot2(1)**2 + Sf5612*MStop2(2)**2

	Sf5284 = Sf5282*MSbot2(1)**2 + Sf5283*MStop2(2)

	Sf4956 = Sf4954*MSbot2(1)**2 - Sf4955*MStop2(2)**2

	Sf4944 = Sf4942*MSbot2(1)**2 - Sf4943*MStop2(2)**2

	Sf6645 = 5*Sf4595*MSbot2(1)**2 - Sf6644*MStop2(2)

	Sf4603 = 5*Sf4595*MSbot2(1)**2 - Sf4602*MStop2(1)

	Sf7052 = Sf4549*MSbot2(1)**2 + Sf7051*MStop2(2)

	Sf4553 = Sf4549*MSbot2(1)**2 + Sf4552*MStop2(1)

	Sf4363 = Sf4361*MSbot2(1)**2 - Sf4362*MStop2(1)**2

	Sf5785 = Sf3787*MSbot2(1)**2 + Sf5784*MStop2(2)

	Sf5779 = 3*Sf3778*MSbot2(1)**2 + Sf5778*MStop2(2)

	Sf3782 = 3*Sf3778*MSbot2(1)**2 + Sf3781*MStop2(1)

	Sf5777 = 3*Sf3773*MSbot2(1)**2 + Sf5776*MStop2(2)

	Sf3776 = 3*Sf3773*MSbot2(1)**2 + Sf3775*MStop2(1)

	Sf5766 = Sf3755*MSbot2(1)**2 + Sf5765*MStop2(2)

	Sf3758 = Sf3755*MSbot2(1)**2 + Sf3757*MStop2(1)

	Sf5626 = Sf3435*MSbot2(1)**2 + Sf5625*MStop2(2)

	Sf5842 = Sf3423*MSbot2(1)**2 + Sf5841*MStop2(2)

	Sf4058 = Sf3423*MSbot2(1)**2 + Sf4057*MStop2(1)

	Sf3414 = Sf3412*MSbot2(1)**2 + Sf3413*MStop2(1)**2

	Sf3273 = Sf3270*MSbot2(1)**2 + Sf3272*MStop2(1)

	Sf5549 = 6*Sf3191*MSbot2(1)**2 - Sf5548*MStop2(2)

	Sf3194 = 6*Sf3191*MSbot2(1)**2 - Sf3193*MStop2(1)

	Sf5547 = 2*Sf3186*MSbot2(1)**2 - Sf5546*MStop2(2)

	Sf3189 = 2*Sf3186*MSbot2(1)**2 - Sf3188*MStop2(1)

	Sf3153 = Sf3150*MSbot2(1)**2 + Sf3152*MStop2(1)

	Sf5484 = 6*Sf299*MSbot2(1)**2 - Sf5483*MStop2(2)

	Sf3088 = 6*Sf299*MSbot2(1)**2 - Sf3087*MStop2(1)

	Sf5286 = Sf2839*MSbot2(1)**2 + Sf5285*MStop2(2)

	Sf2843 = Sf2839*MSbot2(1)**2 + Sf2842*MStop2(1)

	Sf2838 = 2*Sf2834*MSbot2(1)**2 + Sf2837*MStop2(1)

	Sf5251 = 3*Sf2792*MSbot2(1)**2 + Sf5250*MStop2(2)

	Sf2795 = 3*Sf2792*MSbot2(1)**2 + Sf2794*MStop2(1)

	Sf2790 = 2*Sf2786*MSbot2(1)**2 + Sf2789*MStop2(1)

	Sf5247 = Sf2782*MSbot2(1)**2 + Sf5246*MStop2(2)

	Sf2785 = Sf2782*MSbot2(1)**2 + Sf2784*MStop2(1)

	Sf5601 = 6*Sf2543**5*MSbot2(1)**2 + Sf5600*MStop2(2)**3

	Sf3395 = 6*Sf2543**5*MSbot2(1)**2 + Sf3394*MStop2(1)**3

	Sf4978 = 2*Sf2466*MSbot2(1)**2 - 3*Sf4977*MStop2(2)

	Sf2470 = 2*Sf2466*MSbot2(1)**2 - 3*Sf2469*MStop2(1)

	Sf2428 = Sf2426*MSbot2(1)**2 - Sf2427*MStop2(1)**2

	Sf4952 = Sf2415*MSbot2(1)**2 - Sf4951*MStop2(2)

	Sf2420 = Sf2415*MSbot2(1)**2 - Sf2419*MStop2(1)

	Sf2399 = Sf2397*MSbot2(1)**2 - Sf2398*MStop2(1)**2

	Sf7054 = 4*Sf207*MSbot2(1)**2 + Sf7053*MStop2(2)

	Sf4558 = 4*Sf207*MSbot2(1)**2 + Sf4557*MStop2(1)

	Sf5486 = 6*Sf174*MSbot2(1)**2 - Sf5485*MStop2(2)

	Sf3092 = 6*Sf174*MSbot2(1)**2 - Sf3091*MStop2(1)

	Sf7050 = 16*Sf1095*MSbot2(1)**2 + Sf7049*MStop2(2)

	Sf4548 = 16*Sf1095*MSbot2(1)**2 + Sf4547*MStop2(1)

	Sf6305 = 2*Sf1060*MSbot2(1)**2 + Sf6304*MStop2(2)

	Sf4095 = 2*Sf1060*MSbot2(1)**2 + Sf4094*MStop2(1)

	Sf6535 = 210*MSbot2(1)**2 + Sf6534*MStop2(2)

	Sf4307 = 210*MSbot2(1)**2 + Sf4306*MStop2(1)

	Sf5466 = 3*MSbot2(1)**2 - Sf5465*MStop2(2)

	Sf9180 = MSbot2(1)**2 - Sf9121*MSbot2(2)

	Sf8995 = MSbot2(1)**2 - Sf8919*MSbot2(2)

	Sf8709 = Sf8306*MSbot2(1) + Sf8708*MSq2Diff(1,2)**2

	Sf8672 = Sf8282*MSbot2(1) - Sf8671*MSq2Diff(1,1)**2

	Sf6586 = 3*Sf6452*MSbot2(1) + Sf6585*MStop2(2)

	Sf6456 = 3*Sf6452*MSbot2(1) + Sf6455*MStop2(1)

	Sf6584 = Sf6448*MSbot2(1) + Sf6583*MStop2(2)

	Sf6451 = Sf6448*MSbot2(1) + Sf6450*MStop2(1)

	Sf6579 = 3*Sf6440*MSbot2(1) + Sf6578*MStop2(2)

	Sf6443 = 3*Sf6440*MSbot2(1) + Sf6442*MStop2(1)

	Sf6577 = Sf6436*MSbot2(1) + Sf6576*MStop2(2)

	Sf6439 = Sf6436*MSbot2(1) + Sf6438*MStop2(1)

	Sf6573 = 6*Sf6429*MSbot2(1) + Sf6572*MStop2(2)

	Sf6433 = 6*Sf6429*MSbot2(1) + Sf6432*MStop2(1)

	Sf6571 = 6*Sf6425*MSbot2(1) + Sf6570*MStop2(2)

	Sf6428 = 6*Sf6425*MSbot2(1) + Sf6427*MStop2(1)

	Sf5873 = 3*Sf5743*MSbot2(1) + Sf5872*MStop2(2)

	Sf5747 = 3*Sf5743*MSbot2(1) + Sf5746*MStop2(1)

	Sf6642 = 3*Sf4584*MSbot2(1) + Sf6641*MStop2(2)

	Sf4588 = 3*Sf4584*MSbot2(1) + Sf4587*MStop2(1)

	Sf7006 = Sf4408*MSbot2(1) + Sf7005*MStop2(2)

	Sf4411 = Sf4408*MSbot2(1) + Sf4410*MStop2(1)

	Sf7000 = Sf4393*MSbot2(1) + Sf6999*MStop2(2)

	Sf4397 = Sf4393*MSbot2(1) + Sf4396*MStop2(1)

	Sf6998 = 2*Sf4389*MSbot2(1) + 3*Sf6997*MStop2(2)

	Sf4392 = 2*Sf4389*MSbot2(1) + 3*Sf4391*MStop2(1)

	Sf6995 = Sf4383*MSbot2(1) + Sf6994*MStop2(2)

	Sf4386 = Sf4383*MSbot2(1) + Sf4385*MStop2(1)

	Sf6993 = 3*Sf4379*MSbot2(1) + Sf6992*MStop2(2)

	Sf4382 = 3*Sf4379*MSbot2(1) + Sf4381*MStop2(1)

	Sf6552 = Sf4345*MSbot2(1) - Sf6551*MStop2(2)

	Sf4350 = Sf4345*MSbot2(1) - Sf4349*MStop2(1)

	Sf6550 = Sf4338*MSbot2(1) + Sf6549*MStop2(2)

	Sf4341 = Sf4338*MSbot2(1) + Sf4340*MStop2(1)

	Sf6543 = 3*Sf4319*MSbot2(1) + Sf6542*MStop2(2)

	Sf4323 = 3*Sf4319*MSbot2(1) + Sf4322*MStop2(1)

	Sf6541 = Sf4315*MSbot2(1) + Sf6540*MStop2(2)

	Sf4318 = Sf4315*MSbot2(1) + Sf4317*MStop2(1)

	Sf6538 = 3*Sf4310*MSbot2(1) - Sf6537*MStop2(2)

	Sf4313 = 3*Sf4310*MSbot2(1) - Sf4312*MStop2(1)

	Sf6318 = 3*Sf4111*MSbot2(1) + Sf6317*MStop2(2)

	Sf4116 = 3*Sf4111*MSbot2(1) + Sf4115*MStop2(1)

	Sf6316 = Sf4107*MSbot2(1) - Sf6315*MStop2(2)

	Sf4110 = Sf4107*MSbot2(1) - Sf4109*MStop2(1)

	Sf5875 = 3*Sf4097*MSbot2(1) - Sf5874*MStop2(2)

	Sf5751 = 3*Sf4097*MSbot2(1) - Sf5750*MStop2(1)

	Sf4099 = Sf4097*MSbot2(1) + 3*Sf4098*MStop2(1)

	Sf6303 = Sf4088*MSbot2(1) + Sf6302*MStop2(2)

	Sf4091 = Sf4088*MSbot2(1) + Sf4090*MStop2(1)

	Sf5827 = Sf4020*MSbot2(1) - Sf5826*MStop2(2)

	Sf4024 = Sf4020*MSbot2(1) - Sf4023*MStop2(1)

	Sf5822 = 2*Sf4011*MSbot2(1) - Sf5821*MStop2(2)

	Sf4015 = 2*Sf4011*MSbot2(1) - Sf4014*MStop2(1)

	Sf5624 = Sf3430*MSbot2(1) + Sf5623*MStop2(2)

	Sf3434 = Sf3430*MSbot2(1) + Sf3433*MStop2(1)

	Sf5622 = Sf3425*MSbot2(1) + Sf5621*MStop2(2)

	Sf3428 = Sf3425*MSbot2(1) + Sf3427*MStop2(1)

	Sf7004 = Sf3406*MSbot2(1) + Sf7003*MStop2(2)

	Sf5867 = Sf3406*MSbot2(1) + Sf5866*MStop2(2)

	Sf5610 = Sf3406*MSbot2(1) + Sf5609*MStop2(2)

	Sf5726 = Sf3406*MSbot2(1) + Sf5725*MStop2(1)

	Sf4407 = Sf3406*MSbot2(1) + Sf4406*MStop2(1)

	Sf3410 = Sf3406*MSbot2(1) + Sf3409*MStop2(1)

	Sf5608 = 3*Sf3402*MSbot2(1) + Sf5607*MStop2(2)

	Sf3405 = 3*Sf3402*MSbot2(1) + Sf3404*MStop2(1)

	Sf5582 = 3*Sf3265*MSbot2(1) + Sf5581*MStop2(2)

	Sf3269 = 3*Sf3265*MSbot2(1) + Sf3268*MStop2(1)

	Sf5568 = Sf3241*MSbot2(1) + Sf5567*MStop2(2)

	Sf3245 = Sf3241*MSbot2(1) + Sf3244*MStop2(1)

	Sf5566 = Sf3232*MSbot2(1) + Sf5565*MStop2(2)

	Sf3237 = Sf3232*MSbot2(1) + Sf3236*MStop2(1)

	Sf5520 = 3*Sf3145*MSbot2(1) + Sf5519*MStop2(2)

	Sf3149 = 3*Sf3145*MSbot2(1) + Sf3148*MStop2(1)

	Sf5518 = Sf3140*MSbot2(1) + Sf5517*MStop2(2)

	Sf3143 = Sf3140*MSbot2(1) + Sf3142*MStop2(1)

	Sf5505 = Sf3123*MSbot2(1) + Sf5504*MStop2(2)

	Sf3126 = Sf3123*MSbot2(1) + Sf3125*MStop2(1)

	Sf5307 = 3*Sf2869*MSbot2(1) + Sf5306*MStop2(2)

	Sf2874 = 3*Sf2869*MSbot2(1) + Sf2873*MStop2(1)

	Sf5006 = Sf2535*MSbot2(1) + Sf5005*MStop2(2)

	Sf2540 = Sf2535*MSbot2(1) + Sf2539*MStop2(1)

	Sf5004 = 3*Sf2529*MSbot2(1) + Sf5003*MStop2(2)

	Sf2534 = 3*Sf2529*MSbot2(1) + Sf2533*MStop2(1)

	Sf4998 = Sf2516*MSbot2(1) + Sf4997*MStop2(2)

	Sf2519 = Sf2516*MSbot2(1) + Sf2518*MStop2(1)

	Sf4995 = 6*Sf2507*MSbot2(1) + Sf4994*MStop2(2)

	Sf2512 = 6*Sf2507*MSbot2(1) + Sf2511*MStop2(1)

	Sf4993 = 6*Sf2502*MSbot2(1) + Sf4992*MStop2(2)

	Sf2506 = 6*Sf2502*MSbot2(1) + Sf2505*MStop2(1)

	Sf4984 = Sf2484*MSbot2(1) + Sf4983*MStop2(2)

	Sf2488 = Sf2484*MSbot2(1) + Sf2487*MStop2(1)

	Sf4982 = Sf2477*MSbot2(1) + Sf4981*MStop2(2)

	Sf2482 = Sf2477*MSbot2(1) + Sf2481*MStop2(1)

	Sf4976 = Sf2461*MSbot2(1) + Sf4975*MStop2(2)

	Sf2465 = Sf2461*MSbot2(1) + Sf2464*MStop2(1)

	Sf4972 = 3*Sf2451*MSbot2(1) + Sf4971*MStop2(2)

	Sf2456 = 3*Sf2451*MSbot2(1) + Sf2455*MStop2(1)

	Sf4970 = Sf2445*MSbot2(1) + Sf4969*MStop2(2)

	Sf2449 = Sf2445*MSbot2(1) + Sf2448*MStop2(1)

	Sf4950 = Sf2408*MSbot2(1) - Sf4949*MStop2(2)

	Sf2412 = Sf2408*MSbot2(1) - Sf2411*MStop2(1)

	Sf4948 = 3*Sf2402*MSbot2(1) + Sf4947*MStop2(2)

	Sf2406 = 3*Sf2402*MSbot2(1) + Sf2405*MStop2(1)

	Sf4920 = Sf2349*MSbot2(1) + Sf4919*MStop2(2)

	Sf2353 = Sf2349*MSbot2(1) + Sf2352*MStop2(1)

	Sf24 = (1/2.D0*Sf19**2/(MSbot(1)*MSbot(2)) - 
     &     2*MSbot(1)*MSbot(2))/Sf19

	Sf6974 = 9*Sf6973*TB - Sf6971*USbot21

	Sf5161 = -(Sf5148*Sf5160*Ss25**2) + 
     &   Sf173*Sf4802*Sf4864*Sf4906**2*TB

	Sf5154 = Sf5148*Sf5153*Ss25**2 + 
     &   Sf173*Sf4802*Sf4864*Sf4906**2*TB

	Sf5043 = -(Sf4827*Sf5042*Ss25**2) + 
     &   Sf173*Sf4781*Sf4820**2*Sf4828*TB

	Sf5035 = Sf4827*Sf5034*Ss25**2 + 
     &   Sf173*Sf4781*Sf4820**2*Sf4828*TB

	Sf4255 = 9*Sf4254*TB - Sf4252*USbot21

	Sf4246 = 9*Sf4245*TB - Sf4243*USbot11

	Sf2754 = -(Sf2740*Sf2753*Ss25**2) + 
     &   Sf173*Sf2263*Sf2334**2*Sf2741*TB

	Sf2747 = Sf2740*Sf2746*Ss25**2 + 
     &   Sf173*Sf2263*Sf2334**2*Sf2741*TB

	Sf2578 = -(Sf2562*Sf2577*Ss25**2) + 
     &   Sf173*Sf2198*Sf2320**2*Sf2564*TB

	Sf2570 = Sf2562*Sf2569*Ss25**2 + 
     &   Sf173*Sf2198*Sf2320**2*Sf2564*TB

	Sf8738 = 2*MZ2*Sf8737*SW2 + 
     &   Sf7219*Sf8735*MStop2(2)*UStop11*UStop12*UStop21*
     &    UStop22

	Sf8506 = MZ2*Sf8505*SW2 + 
     &   Sf4795*Sf8491*MSbot2(2)*USbot11*USbot12*USbot21*
     &    USbot22

	Sf7767 = MZ2*Sf7766*SW2 + 
     &   Sf1775*Sf7764*MStau2(2)*UStau11*UStau12*UStau21*
     &    UStau22

	Sf5459 = -((Sf4911*Sf5458)/(Sf4906**2*TB)) + 
     &   (MW2*Sf4802*Sf4864*MHiggs2(4))/Ss25**2

	Sf5452 = -((Sf4905*Sf5451)/(Sf4820**2*TB)) + 
     &   (MW2*Sf4781*Sf4828*MHiggs2(4))/Ss25**2

	Sf3071 = -((Sf2340*Sf3070)/(Sf2334**2*TB)) + 
     &   (MW2*Sf2263*Sf2741*MHiggs2(4))/Ss25**2

	Sf3063 = -((Sf2327*Sf3062)/(Sf2320**2*TB)) + 
     &   (MW2*Sf2198*Sf2564*MHiggs2(4))/Ss25**2

	Sf6249 = Sf6248*Ss241*Ss242 + Ss25*MStop2(2)**3

	Sf6141 = Sf6140*Ss241*Ss242 + Ss25*MStop2(1)**3

	Sf5885 = Sf5884*Ss241*Ss242 + Ss25*MStop2(2)**3

	Sf5758 = Sf5757*Ss241*Ss242 + Ss25*MStop2(1)**3

	Sf5137 = Sf5136*Ss241*Ss242 + Sf217*MStop2(2)**3

	Sf5115 = 3*Sf5114*Ss241*Ss242 - 3*Sf217*MStop2(2)**3

	Sf5011 = Sf5010*Ss241*Ss242 + Sf217*MStop2(2)**3

	Sf4989 = 3*Sf4988*Ss241*Ss242 - 3*Sf217*MStop2(2)**3

	Sf2728 = (-MW2 + MHiggs2(4))*
     &   (MHiggs2(4)*(3*MW2**3 - 7*MW2**2*MSbot2(2) - 
     &        MSq2Diff(2,1)**3 + 
     &        MW2*MSq2Diff(2,1)*(5*MSbot2(2) + 4*MStop2(1))) + 
     &     MW2*(-4*MW2**3 + 2*MSq2Diff(2,1)**3 + 
     &        MW2**2*(10*MSbot2(2) + MStop2(1)) - 
     &        MW2*MSq2Diff(2,1)*(8*MSbot2(2) + 5*MStop2(1))))

	Sf2689 = 3*(-MW2 + MHiggs2(4))*
     &   (MHiggs2(4)*(3*MW2**3 + MSq2Diff(2,1)**3 - 
     &        7*MW2**2*MStop2(1) - 
     &        MW2*MSq2Diff(2,1)*(4*MSbot2(2) + 5*MStop2(1))) + 
     &     MW2*(-4*MW2**3 - 2*MSq2Diff(2,1)**3 + 
     &        MW2*MSq2Diff(2,1)*(5*MSbot2(2) + 8*MStop2(1)) + 
     &        MW2**2*(MSbot2(2) + 10*MStop2(1))))

	Sf2549 = Sf2548*Ss241*Ss242 + Sf217*MStop2(1)**3

	Sf2497 = 3*Sf2496*Ss241*Ss242 - 3*Sf217*MStop2(1)**3

	Sf4156 = 3*Sf4044*Ss215*Ss218 - 4*Sf4048*MSbot2(2)

	Sf4155 = 3*Sf4044*Ss215*Ss218 - Sf4046*MSbot2(2)

	Sf4049 = 3*Sf4044*Ss215*Ss218 - 4*Sf4048*MSbot2(1)

	Sf4047 = 3*Sf4044*Ss215*Ss218 - Sf4046*MSbot2(1)

	Sf1510 = 3*MW2*Ss115 + Sf1509*MStau2(2)

	Sf615 = 3*MW2*Ss115 + Sf614*MStau2(1)

	Sf3715 = 3*MW2*Ss115 + Sf3714*MSbot2(2)

	Sf3443 = 3*MW2*Ss115 + Sf3442*MSbot2(1)

	Sf10002 = Sf9344*Sf9993 + 2*Sf9986*MStau2(1)*MStau2(2)

	Sf9975 = ((2*Sf9963)/Sf9317**2 - 
     &     (MZ2*Sf9974)/(Sf9331**2*MHiggs2(2)*MStau2(1)*MStau2(2)))
     &    /MZ2

	Sf9952 = Sf9940 + Sf9951*MStau2(2)

	Sf1806 = 15*Sf992 + 28*MStau2(2)

	Sf994 = 15*Sf992 + 28*MStau2(1)

	Sf9927 = (-(Sf9926/Sf9317**2) + 
     &     (MZ2*Sf9916)/(Sf9331**2*MHiggs2(2)))/MZ2

	Sf9928 = ((Sf8929*Sf8930*Sf9905)/MHiggs2(2) + 
     &     (2*MZ2*Sf55*Sf57*Sf9927)/(MStau2(1)*MStau2(2)))/MZ2

	Sf9881 = MZ2*Sf9878 - Sf9880*MStau2(2)**2

	Sf9862 = Sf9861 + MHiggs2(2)*MStau2(2)

	Sf9863 = Sf9860 + Sf9862*MStau2(2)

	Sf9859 = (2*Sf9853*Sf9858)/MHiggs2(2)**2 + 
     &   Sf9611*Ss928*MStau2(1)

	Sf9854 = (2*Sf9848*Sf9853)/MHiggs2(2)**2 - 
     &   Sf9598*Ss949*MStau2(2)

	Sf9836 = (Sf9835 + (2*Sf9344*Sf9832)/MHiggs2(2))/Sf9344

	Sf1803 = Sf982 + 28*MStau2(2)

	Sf985 = Sf982 + 28*MStau2(1)

	Sf9815 = (-Sf9814 + (2*Sf9344*Sf9811)/MHiggs2(2))/Sf9344

	Sf1798 = Sf973 + 7*MStau2(2)

	Sf974 = Sf973 + 7*MStau2(1)

	Sf9705 = Sf55*Sf57*Sf9344*Sf9697 - Sf9704/MStau2(1)

	Sf9706 = (Sf9705 + 
     &     (2*Sf55*Sf57*Sf9344**2*Sf9693)/MHiggs2(2))/Sf9344**2

	Sf9683 = (-Sf9682 + (2*Sf9317*Sf9678)/MHiggs2(2))/Sf9317

	Sf9653 = (Sf9652 + (2*Sf9317*Sf9647)/MHiggs2(2))/Sf9317

	Sf1796 = 15*Sf963 + 28*MStau2(2)

	Sf964 = 15*Sf963 + 28*MStau2(1)

	Sf9623 = Sf9621 - Sf9622*MStau2(1)

	Sf9624 = Sf9363*Sf9364*Sf9619 - Sf9623*MStau2(1)

	Sf9616 = Sf9363*Sf9364*Sf9614 - Sf9615*MStau2(1)

	Sf9605 = Sf9602 - Sf9604*MStau2(2)

	Sf9606 = Sf9332*Sf9333*Sf9600 - Sf9605*MStau2(2)

	Sf1794 = 3*Sf957 + 14*MStau2(2)

	Sf958 = 3*Sf957 + 14*MStau2(1)

	Sf1795 = 2*MW2*Sf955 + Sf1794*MStau2(2)

	Sf959 = 2*MW2*Sf955 + Sf958*MStau2(1)

	Sf9471 = 2*Sf9466 - Sf9470*MStau2(2)

	Sf9472 = Sf9462 - Sf9471*MStau2(2)

	Sf9473 = Sf9456 - Sf9472*MStau2(2)

	Sf9474 = Sf9448 - Sf9473*MStau2(2)

	Sf9414 = Sf9413 - Sf55*Sf57*Sf9317*Sf9404*MStau2(2)

	Sf9415 = (2*Sf55*Sf57*Sf9398)/MHiggs2(2) - 
     &   Sf9414/(Sf9317**2*MStau2(2))

	Sf9388 = MHiggs2(2)*
     &    (2*MZ2**3 - 2*MZ2**2*MStau2(2) + 
     &      MZ2*(3*MStau2(1) + MStau2(2))*MStau2(3) - MStau2(3)**3)
     &     + MZ2*(-2*MZ2**3 - MZ2**2*(-2*MStau2(2) - MStau2(3)) - 
     &      MZ2*(5*MStau2(1) + 3*MStau2(2))*MStau2(3) + 
     &      2*MStau2(3)**3)

	Sf9381 = Sf9290*Sf9291*Sf9377 - Sf9380*MStau2(2)

	Sf9368 = Sf9363*Sf9364*Sf9365 - Sf9367*MStau2(1)

	Sf9357 = MZ2**3 - MZ2**2*(5*MStau2(1) + 3*MStau2(2)) - 
     &   MHiggs2(2)*MStau2(3)**2 + 
     &   MZ2*(MHiggs2(2)*(3*MStau2(1) + MStau2(2)) + 
     &      2*MStau2(3)**2)

	Sf10003 = 2*MZ2*Sf10001*Sf9344**2 - 
     &   2*Sf10002*Sf9362**2*MHiggs2(2)

	Sf10052 = ((2*Sf10042)/Sf9344**2 - 
     &     (MZ2*Sf10051)/(Sf9362**2*MHiggs2(2)*MStau2(1)*MStau2(2))
     &     )/MZ2

	Sf10024 = (Sf10023/Sf9344**2 - 
     &     (MZ2*Sf10015)/(Sf9362**2*MHiggs2(2)))/MZ2

	Sf9337 = Sf9332*Sf9333*Sf9334 - Sf9336*MStau2(2)

	Sf9326 = Sf9290*Sf9291*Sf9323 + Sf9325*MStau2(2)

	Sf10093 = Sf10091*Sf9302*Sf9303 + Sf10092*MStau2(1)

	Sf936 = Sf31*Sf347*Sf930 - 
     &   27*C2B*Sf901*MSl2Diff(1,1)**3*UStau22**2

	Sf10137 = MHiggs2(1)*
     &    (2*MZ2**3 - 2*MZ2**2*MStau2(1) - 
     &      MZ2*(MStau2(1) + 3*MStau2(2))*MStau2(3) + MStau2(3)**3)
     &     - MZ2*(2*MZ2**3 - 
     &      MZ2*(3*MStau2(1) + 5*MStau2(2))*MStau2(3) + 
     &      2*MStau2(3)**3 + MZ2**2*(-2*MStau2(1) + MStau2(3)))

	Sf10063 = MZ2**3 - MZ2**2*(3*MStau2(1) + 5*MStau2(2)) - 
     &   MHiggs2(1)*MStau2(3)**2 + 
     &   MZ2*(MHiggs2(1)*(MStau2(1) + 3*MStau2(2)) + 
     &      2*MStau2(3)**2)

	Sf9281 = Sf9280 + MHiggs2(1)*MSbot2(2)

	Sf9282 = Sf9279 + Sf9281*MSbot2(2)

	Sf9278 = (2*Sf9271*Sf9277)/MHiggs2(1)**2 - 
     &   Sf9273*Ss906*MSbot2(2)

	Sf9272 = (2*Sf9266*Sf9271)/MHiggs2(1)**2 + 
     &   Sf9260*Ss885*MSbot2(1)

	Sf9254 = Sf9250 - Sf9253*MSbot2(1)

	Sf9247 = Sf9161*Sf9162*Sf9245 - Sf9246*MSbot2(1)

	Sf9236 = Sf9233 - Sf9235*MSbot2(2)

	Sf9237 = Sf9123*Sf9124*Sf9231 - Sf9236*MSbot2(2)

	Sf9208 = Sf9204 + Sf9207*MSbot2(1)

	Sf9202 = Sf8873*Sf8874*Sf9197 - Sf9201*MSbot2(2)

	Sf9166 = Sf9161*Sf9162*Sf9163 - Sf9165*MSbot2(1)

	Sf9155 = Sf8885*Sf8886*Sf9153 + Sf9154*MSbot2(1)

	Sf9130 = Sf9123*Sf9124*Sf9126 - Sf9129*MSbot2(2)

	Sf9117 = Sf8873*Sf8874*Sf9114 + Sf9116*MSbot2(2)

	Sf9100 = Sf9099 + MHiggs2(2)*MSbot2(2)

	Sf9101 = Sf9098 + Sf9100*MSbot2(2)

	Sf9097 = (2*Sf9091*Sf9096)/MHiggs2(2)**2 + 
     &   Sf9037*Ss885*MSbot2(1)

	Sf9092 = (2*Sf9086*Sf9091)/MHiggs2(2)**2 - 
     &   Sf9023*Ss906*MSbot2(2)

	Sf9049 = Sf9047 - Sf9048*MSbot2(1)

	Sf9050 = Sf8971*Sf8972*Sf9045 - Sf9049*MSbot2(1)

	Sf9042 = Sf8971*Sf8972*Sf9040 - Sf9041*MSbot2(1)

	Sf9031 = Sf9028 - Sf9030*MSbot2(2)

	Sf9032 = Sf8921*Sf8922*Sf9026 - Sf9031*MSbot2(2)

	Sf9011 = Sf8885*Sf8886*Sf9008 - Sf9010*MSbot2(1)

	Sf9004 = Sf8873*Sf8874*Sf9000 - Sf9003*MSbot2(2)

	Sf8976 = Sf8971*Sf8972*Sf8973 - Sf8975*MSbot2(1)

	Sf8965 = Sf8885*Sf8886*Sf8963 + Sf8964*MSbot2(1)

	Sf10025 = ((Sf10006*Sf8929*Sf8930)/MHiggs2(2) + 
     &     (2*MZ2*Sf10024*Sf55*Sf57)/(MStau2(1)*MStau2(2)))/MZ2

	Sf8928 = Sf8921*Sf8922*Sf8924 - Sf8927*MSbot2(2)

	Sf8915 = Sf8873*Sf8874*Sf8912 + Sf8914*MSbot2(2)

	Sf8860 = 2*Sf1*Sf7756*Sf8763 - 161*Sf8819*Ss443**2*UStau22**2

	Sf8732 = 3*Sf8731 - 8*MW2**2*Sf8722*UStop22**2

	Sf8720 = CW2*MZ2*Sf8718 - 12*Sf8719*Ss922

	Sf8628 = MZ2*Sf8622

	Sf8620 = MZ2*Sf8614 + 6*Sf8619

	Sf8587 = 29*MZ2

	Sf8574 = 3*MZ2

	Sf855 = Sf854 + Sf548*MSneu2(1)

	Sf856 = Sf852 - Sf855*MSneu2(1)

	Sf8501 = MW2*Sf8493 - 9*MZ2**4*Sf8500*USbot11

	Sf857 = Sf167**2*Sf849 - Sf856*MSneu2(1)

	Sf782 = ML2*Sf777 + Sf781*MSl2Diff(1,1)

	Sf7493 = 200*MZ2 - 3*Sf7492

	Sf7462 = 40*MZ2 - 3*Sf7461

	Sf7418 = 3*MW2*Sf7405 - Sf7399*MSq2Diff(1,1)*MSq2Diff(2,1)

	Sf737 = Sf736 - MHiggs2(4)*MSneu2(1)**2

	Sf7317 = Sf7311 - 
     &   (3*MZ2*Sf7316*SW2)/(MSl2Diff(1,1)*MSl2Diff(2,1))

	Sf7282 = 3*Sf7270 - 2*Sf7281

	Sf7032 = Sf7031 - MHiggs2(4)*MStop2(2)**2

	Sf7010 = Sf7009 - MHiggs2(4)*MStop2(2)**2

	Sf2178 = 2*CW2*Sf1923*Sf2177 + 3*MW2*Sf7*MHiggs2(4)

	Sf6962 = -9*MW2*Sf6943*Sf6957*TB + Sf6961*MHiggs2(4)

	Sf6904 = Sf6900 + 4*Sf6903*MSq2Diff(2,2)**5

	Sf695 = MW2*Sf694 + 7*Sf34*Sf689*MHiggs2(4)

	Sf6845 = Sf6844 - MHiggs2(4)*MStop2(2)**2

	Sf6799 = MW2*Sf6798 + Sf6796*MHiggs2(4)

	Sf6816 = Sf6787 + 2*Sf4286*MSq2Diff(2,2)**2

	Sf6731 = Sf6730 - MHiggs2(4)*MStop2(1)**2

	Sf696 = Sf673 + MW2*Sf688

	Sf6653 = Sf6649 + 4*Sf6652*MSq2Diff(1,2)**5

	Sf6590 = Sf6589 - MHiggs2(4)*MStop2(2)**2

	Sf6544 = MW2*Sf6543 + Sf6541*MHiggs2(4)

	Sf6561 = Sf6532 + 2*Sf4286*MSq2Diff(1,2)**2

	Sf697 = MW2*Sf648 - Sf668*MSl2Diff(1,1)

	Sf6464 = Sf6463 - MHiggs2(4)*MStop2(1)**2

	Sf698 = Sf639 + Ss340*MSl2Diff(1,1)**2

	Sf6312 = Sf6308 - Sf6311*MStop2(2)

	Sf6313 = Sf3204*Sf3205*Sf6307 - Sf6312*MStop2(2)

	Sf6212 = MW2*Sf6211 + 2*Sf6207*MSq2Diff(2,2)

	Sf6235 = Sf6119 + 5*Sf6234*MStop2(2)

	Sf6122 = Sf6119 + 5*Sf6121*MStop2(1)

	Sf6236 = Sf6117 - Sf6235*MStop2(2)

	Sf6123 = Sf6117 - Sf6122*MStop2(1)

	Sf6237 = (MW2 - 5*MSq2Diff(2,2))*
     &    (MW2 - MSq2Diff(2,2))**3 - 
     &   MW2*(3*MW2*(-MBot2 + MTop2 + 4*MW2 + C2B*MW2) + 
     &      (-3*MBot2 + 3*MTop2 - 37*MW2 + 3*C2B*MW2)*MSq2Diff(2,2) + 
     &      31*MSq2Diff(2,2)**2)*MStop2(2) + 32*MW2**2*MStop2(2)**2

	Sf6124 = MW2**4 + 5*MSq2Diff(2,1)**4 - 
     &   MW2**3*(8*MSbot2(2) + (4 + 3*C2B)*MStop2(1)) + 
     &   MW2**2*(18*MSbot2(2)**2 + 
     &      (3*MBot2 - 3*MTop2 + (1 - 3*C2B)*MSbot2(2))*MStop2(1) + 
     &      (13 + 3*C2B)*MStop2(1)**2) - 
     &   MW2*MSq2Diff(2,1)*
     &    (16*MSbot2(2)**2 - MSbot2(2)*MStop2(1) - 
     &      3*MStop2(1)*(MBot2 - MTop2 + 5*MStop2(1)))

	Sf6092 = Sf6091 + 84*MSbot2(2)**5*USbot21

	Sf6065 = MW2*Sf6063 + 3*Sf6064*MSbot2(2)

	Sf6226 = Sf5857 - Sf6225*MSbot2(2)

	Sf5860 = Sf5857 - Sf5859*MSbot2(1)

	Sf6227 = 8*MW2**5 + 
     &   MW2**4*((-22 + 9*C2B)*MSq2Diff(2,2) - 48*MStop2(2)) + 
     &   20*MSq2Diff(2,2)**4*(MSbot2(2) + MStop2(2)) - 
     &   MW2**3*MSq2Diff(2,2)*
     &    (9*MBot2 - 9*MTop2 + (-31 + 9*C2B)*MSq2Diff(2,2) + 
     &      18*C2B*MStop2(2)) - 
     &   2*MW2*MSq2Diff(2,2)*
     &    (19*MSq2Diff(2,2)**3 + 
     &      12*MSq2Diff(2,2)*(-MBot2 + MTop2 + 10*MSq2Diff(2,2))*
     &       MStop2(2) + 4*(-3*MBot2 + 3*MTop2 + 31*MSq2Diff(2,2))*
     &       MStop2(2)**2) + 
     &   MW2**2*(MSq2Diff(2,2)**2*
     &       (9*MBot2 - 9*MTop2 + MSq2Diff(2,2)) + 
     &      MSq2Diff(2,2)*(18*MBot2 - 18*MTop2 + 
     &         (110 - 24*C2B)*MSq2Diff(2,2))*MStop2(2) + 
     &      (352 - 24*C2B)*MSq2Diff(2,2)*MStop2(2)**2 + 
     &      256*MStop2(2)**3)

	Sf5861 = Sf5849*Sf5850*Sf5852 - Sf5860*MSbot2(1)

	Sf5835 = MW2*Sf5834 + 2*Sf5829*MSq2Diff(1,2)

	Sf5786 = MTop2*Sf5783 - Sf5785*MSq2Diff(1,2)

	Sf5787 = MBot2*Sf5782 - Sf5786*MSq2Diff(1,2)

	Sf5735 = Sf5732 + 5*Sf5734*MStop2(1)

	Sf5736 = Sf5730 - Sf5735*MStop2(1)

	Sf5737 = MW2**4 + 5*MSq2Diff(1,1)**4 - 
     &   MW2**3*(8*MSbot2(1) + (4 + 3*C2B)*MStop2(1)) + 
     &   MW2**2*(18*MSbot2(1)**2 + 
     &      (3*MBot2 - 3*MTop2 + (1 - 3*C2B)*MSbot2(1))*MStop2(1) + 
     &      (13 + 3*C2B)*MStop2(1)**2) - 
     &   MW2*MSq2Diff(1,1)*
     &    (16*MSbot2(1)**2 - MSbot2(1)*MStop2(1) - 
     &      3*MStop2(1)*(MBot2 - MTop2 + 5*MStop2(1)))

	Sf5703 = Sf5702 + 84*MSbot2(1)**5*USbot11

	Sf5644 = MW2*Sf5642 + 3*Sf5643*MSbot2(1)

	Sf5997 = Sf5576 - Sf5996*MSbot2(2)

	Sf5578 = Sf5576 - Sf5577*MSbot2(1)

	Sf5998 = -((MW2 - MSq2Diff(2,2))**2*
     &      (-18*MBot2*MW2**2 + 18*MTop2*MW2**2 + 18*C2B*MW2**3 + 
     &        MSq2Diff(2,2)**3 + 
     &        MW2*(-13*MSq2Diff(2,2)**2 + 
     &           MW2*(MW2 + 29*MSq2Diff(2,2))))) + 
     &   2*MW2*(MW2**2*(-36*MBot2 + 36*MTop2 + MW2 + 36*C2B*MW2) + 
     &      MW2*(45*MBot2 - 45*MTop2 + 71*MW2 - 45*C2B*MW2)*
     &       MSq2Diff(2,2) - 94*MW2*MSq2Diff(2,2)**2 + 
     &      13*MSq2Diff(2,2)**3)*MStop2(2) - 
     &   2*MW2**2*(-45*MBot2 + 45*MTop2 - 4*MW2 + 45*C2B*MW2 + 
     &      89*MSq2Diff(2,2))*MStop2(2)**2

	Sf5579 = Sf5029**2*Sf5572 - Sf5578*MSbot2(1)

	Sf549 = Sf544 - Sf548*MSneu2(1)

	Sf5231 = SBA*Sf5227 + CBA*Sf5230

	Sf5216 = Sf5208 + Sf5200*Ss25

	Sf5214 = Sf5196 + Sf5188*Ss25

	Sf5170 = (Sf5148*Sf5168)/TB + 
     &   (2*MW2*Sf4864*Sf5169)/MHiggs2(4)**2

	Sf520 = Sf514 - Sf519*MSneu2(1)

	Sf5071 = Sf5067 - 21*MTop2*MSq2Diff(2,2)

	Sf5019 = (Sf4827*Sf5017)/TB + 
     &   (2*MW2*Sf4828*Sf5018)/MHiggs2(4)**2

	Sf4945 = Sf4941 - 21*MTop2*MSq2Diff(1,2)

	Sf6412 = Sf4871 + Sf4873*Ss25

	Sf6930 = Sf4866 + Sf4869*Ss25

	Sf6410 = Sf4836 + Sf4838*Ss25

	Sf6683 = Sf4830 + Sf4834*Ss25

	Sf4805 = Sf1743*Sf4791 - Sf1749*Sf4794

	Sf4784 = Sf1743*Sf4767 - Sf1749*Sf4771

	Sf6918 = Sf4729 - 5*Sf4600*MStop2(2)

	Sf4730 = Sf4729 - 5*Sf4600*MStop2(1)

	Sf4708 = Sf4704 + 4*Sf4707*MSq2Diff(2,1)**5

	Sf1440 = 5*ML2*Sf47**2 - Sf1439*MSl2Diff(2,1)

	Sf468 = 18*MW2**2*Sf162*Sf467*TB + Sf466*MHiggs2(4)

	Sf6667 = Sf4634 - 5*Sf4600*MStop2(2)

	Sf4635 = Sf4634 - 5*Sf4600*MStop2(1)

	Sf4611 = Sf4607 + 4*Sf4610*MSq2Diff(1,1)**5

	Sf460 = 54*MW2**2*Sf162*Sf455*TB + Sf459*MHiggs2(4)

	Sf4525 = Sf4524 - MHiggs2(4)*MStop2(1)**2

	Sf4450 = MW2*Sf4449 + Sf4445*MHiggs2(4)

	Sf445 = Sf444*MHiggs2(4) + Sf443*UStau11

	Sf4477 = Sf4431 + 2*Sf4286*MSq2Diff(2,1)**2

	Sf4421 = Sf4420 - MHiggs2(4)*MStop2(1)**2

	Sf1725 = Sf1298*Sf1720*Sf44 - 
     &   27*C2B*Sf1700*MSl2Diff(2,1)**3*UStau21**2

	Sf4324 = MW2*Sf4323 + Sf4318*MHiggs2(4)

	Sf4366 = Sf4303 + 2*Sf4286*MSq2Diff(1,1)**2

	Sf6345 = Sf4202 - Sf6344*MStop2(2)

	Sf4204 = Sf4202 - Sf4203*MStop2(1)

	Sf6216 = Sf4156 + Sf4051*MStop2(2)

	Sf4157 = Sf4156 + Sf4051*MStop2(1)

	Sf4151 = MW2*Sf4150 + 2*Sf4145*MSq2Diff(2,1)

	Sf6322 = Sf4122 - Sf6321*MStop2(2)

	Sf4124 = Sf4122 - Sf4123*MStop2(1)

	Sf4167 = Sf4079 - Sf4166*MSbot2(2)

	Sf4082 = Sf4079 - Sf4081*MSbot2(1)

	Sf4168 = Sf4065*Sf4066*Sf4069 - Sf4167*MSbot2(2)

	Sf4083 = 8*MW2**5 + 
     &   20*MSq2Diff(1,1)**4*(MSbot2(1) + MStop2(1)) + 
     &   MW2**4*((-22 + 9*C2B)*MSbot2(1) - 26*MStop2(1) - 
     &      9*C2B*MStop2(1)) - 
     &   2*MW2*MSq2Diff(1,1)*
     &    (19*MSbot2(1)**3 + 
     &      3*MSbot2(1)*(-4*MBot2 + 4*MTop2 + 21*MSbot2(1))*
     &       MStop2(1) - 59*MSbot2(1)*MStop2(1)**2 - 
     &      23*MStop2(1)**3) + 
     &   MW2**2*(MSbot2(1)**2*(9*MBot2 - 9*MTop2 + MSbot2(1)) - 
     &      (-107 + 24*C2B)*MSbot2(1)**2*MStop2(1) + 
     &      3*(-3*MBot2 + 3*MTop2 + (45 + 8*C2B)*MSbot2(1))*
     &       MStop2(1)**2 + 13*MStop2(1)**3) - 
     &   MW2**3*MSq2Diff(1,1)*
     &    (9*MBot2 - 9*MTop2 + (-31 + 9*C2B)*MSq2Diff(1,2) + 
     &      31*MStop2(1) - 31*MStop2(2) + 
     &      C2B*(9*MStop2(1) + 9*MStop2(2)))

	Sf5839 = Sf4049 + Sf4051*MStop2(2)

	Sf4052 = Sf4049 + Sf4051*MStop2(1)

	Sf4040 = MW2*Sf4039 + 2*Sf4031*MSq2Diff(1,1)

	Sf3985 = Sf3984 + 33*MSbot2(2)*USbot21

	Sf3982 = 3*Sf3981 - 11*MSbot2(2)*USbot21

	Sf6282 = Sf3951 - Sf6281*MStop2(2)

	Sf3971 = Sf3951 - Sf3970*MStop2(1)

	Sf3888 = Sf3887 + 33*MSbot2(1)*USbot11

	Sf3885 = 3*Sf3884 - 11*MSbot2(1)*USbot11

	Sf6268 = Sf3836 - Sf6267*MStop2(2)

	Sf3874 = Sf3836 - Sf3873*MStop2(1)

	Sf4212 = Sf3748 + Sf3744*Ss25

	Sf4209 = Sf3739 + Sf3735*Ss25

	Sf3954 = 2*Sf3669 + 7*MSbot2(2)*USbot21

	Sf3670 = Sf3669 - 5*MSbot2(2)*USbot21

	Sf3667 = 3*Sf3666 - MSbot2(2)*USbot21

	Sf6007 = Sf3608 + Sf5951*MStop2(2)

	Sf3609 = Sf3608 + Sf3526*MStop2(1)

	Sf6008 = Sf3606 - Sf6007*MStop2(2)

	Sf3610 = Sf3606 - Sf3609*MStop2(1)

	Sf6009 = Sf2722**2*Sf3604 - Sf6008*MStop2(2)

	Sf3611 = Sf2722**2*Sf3604 - Sf3610*MStop2(1)

	Sf5952 = Sf3524 - Sf5951*MStop2(2)

	Sf3527 = Sf3524 - Sf3526*MStop2(1)

	Sf5953 = Sf3522 - Sf5952*MStop2(2)

	Sf3528 = Sf3522 - Sf3527*MStop2(1)

	Sf5954 = Sf2722**2*Sf3520 - Sf5953*MStop2(2)

	Sf3529 = Sf2722**2*Sf3520 - Sf3528*MStop2(1)

	Sf535 = 5*ML2*Sf34**2 - Sf534*MSl2Diff(1,1)

	Sf3841 = 2*Sf3389 + 7*MSbot2(1)*USbot11

	Sf3390 = Sf3389 - 5*MSbot2(1)*USbot11

	Sf3387 = 3*Sf3386 - MSbot2(1)*USbot11

	Sf5587 = Sf3281 + Sf5525*MStop2(2)

	Sf3282 = Sf3281 + Sf3164*MStop2(1)

	Sf5588 = Sf3279 - Sf5587*MStop2(2)

	Sf3283 = Sf3279 - Sf3282*MStop2(1)

	Sf5589 = Sf2543**2*Sf3277 - Sf5588*MStop2(2)

	Sf3284 = Sf2543**2*Sf3277 - Sf3283*MStop2(1)

	Sf3592 = Sf3256 - Sf3591*MSbot2(2)

	Sf3258 = Sf3256 - Sf3257*MSbot2(1)

	Sf3593 = Sf2563**2*Sf3249 - Sf3592*MSbot2(2)

	Sf3259 = Sf2563**2*Sf3249 - Sf3258*MSbot2(1)

	Sf5526 = Sf3162 - Sf5525*MStop2(2)

	Sf3165 = Sf3162 - Sf3164*MStop2(1)

	Sf5527 = Sf3160 - Sf5526*MStop2(2)

	Sf3166 = Sf3160 - Sf3165*MStop2(1)

	Sf5528 = Sf2543**2*Sf3158 - Sf5527*MStop2(2)

	Sf3167 = Sf2543**2*Sf3158 - Sf3166*MStop2(1)

	Sf5437 = Sf3043 - Sf2886*MStop2(2)

	Sf3045 = Sf3043 - Sf2886*MStop2(1)

	Sf1811 = 2*Sf3 + 23*MStau2(2)

	Sf1000 = 2*Sf3 + 23*MStau2(1)

	Sf5331 = Sf2916 - Sf2886*MStop2(2)

	Sf2918 = Sf2916 - Sf2886*MStop2(1)

	Sf2762 = (Sf2740*Sf2760)/TB + 
     &   (2*MW2*Sf2263*Sf2761)/MHiggs2(4)**2

	Sf2622 = Sf2617 - 21*MTop2*MSq2Diff(2,1)

	Sf5061 = Sf2609 + Sf2379*MStop2(2)

	Sf2610 = Sf2609 + Sf2379*MStop2(1)

	Sf5062 = Sf2608 - Sf5061*MStop2(2)

	Sf2611 = Sf2608 - Sf2610*MStop2(1)

	Sf1144 = Sf1132*Sf259*Sf260 - Sf1143*MSneu2(1)

	Sf2588 = (Sf2562*Sf2585)/TB + 
     &   (2*MW2*Sf2198*Sf2587)/MHiggs2(4)**2

	Sf266 = Sf258 + Sf255*MSneu2(1)

	Sf2400 = Sf2391 - 21*MTop2*MSq2Diff(1,1)

	Sf4935 = Sf2378 + Sf2379*MStop2(2)

	Sf2380 = Sf2378 + Sf2379*MStop2(1)

	Sf4936 = Sf2375 - Sf4935*MStop2(2)

	Sf2381 = Sf2375 - Sf2380*MStop2(1)

	Sf4883 = Sf2287 - Sf2201*MStop2(2)*USbot21

	Sf2294 = Sf2287 - Sf2201*MStop2(1)*USbot21

	Sf4754 = Sf2274 + Sf2276*Ss25

	Sf6929 = Sf2269 + Sf2272*Ss25

	Sf4849 = Sf2226 - Sf2201*MStop2(2)*USbot11

	Sf2242 = Sf2226 - Sf2201*MStop2(1)*USbot11

	Sf4746 = Sf2211 + Sf2213*Ss25

	Sf6679 = Sf2205 + Sf2208*Ss25

	Sf2137 = Sf2136 + MHiggs2(4)*MSneu2(1)**2

	Sf2130 = MW2*Sf2126 - ML2*Sf2129

	Sf2057 = Sf2056 + MHiggs2(4)*MSneu2(1)**2

	Sf2050 = MW2*Sf2046 - ML2*Sf2049

	Sf2004 = Sf2003 + 8*ML*Sf215*MSneu2(1)*UStau22

	Sf1792 = 3*Sf198 + 4*MStau2(2)

	Sf952 = 3*Sf198 + 4*MStau2(1)

	Sf1955 = Sf1954 + 8*ML*Sf215*MSneu2(1)*UStau12

	Sf1885 = Sf1238*Sf1239*Sf1878 - Sf1884*MSneu2(1)

	Sf1870 = 2*Sf1868 + Sf1110*MStau2(2)*UStau21

	Sf183 = MW2**2*Sf182 + Sf176*Sf180*MHiggs2(4)

	Sf1783 = Sf1743*Sf1771 - Sf1749*Sf1774

	Sf1763 = Sf1743*Sf1747 - Sf1749*Sf1753

	Sf1197 = Sf172*MHiggs2(4) + MStau2(2)**2

	Sf178 = Sf172*MHiggs2(4) + MStau2(1)**2

	Sf1671 = Sf1670 + Sf1451*MSneu2(1)

	Sf1672 = Sf1668 - Sf1671*MSneu2(1)

	Sf1673 = Sf1196**2*Sf1666 - Sf1672*MSneu2(1)

	Sf1616 = ML2*Sf1614 + Sf1615*MSl2Diff(2,1)

	Sf1597 = Sf1596 - MHiggs2(4)*MSneu2(1)**2

	Sf159 = Sf158*MHiggs2(4) - 2*Sf157*MSneu2(1)

	Sf1563 = MW2*Sf1562 + 7*Sf1559*Sf47*MHiggs2(4)

	Sf1564 = Sf1547 + MW2*Sf1558

	Sf1565 = MW2*Sf1535 - Sf1542*MSl2Diff(2,1)

	Sf1566 = Sf1532 + Ss340*MSl2Diff(2,1)**2

	Sf1452 = Sf1448 - Sf1451*MSneu2(1)

	Sf1431 = Sf1425 - Sf1430*MSneu2(1)

	Sf1385 = 18*MW2**2*Sf1191*Sf1384*TB + Sf1383*MHiggs2(4)

	Sf1377 = 54*MW2**2*Sf1191*Sf1372*TB + Sf1376*MHiggs2(4)

	Sf1363 = Sf1362*MHiggs2(4) + Sf443*UStau21

	Sf1243 = Sf1241 - Sf1242*MSneu2(1)

	Sf1244 = Sf1238*Sf1239*Sf1240 + Sf1243*MSneu2(1)

	Sf1200 = MW2**2*Sf1199 + Sf176*Sf180*MHiggs2(4)

	Sf1189 = Sf1188*MHiggs2(4) - 2*Sf1187*MSneu2(1)

	Sf1145 = Sf1144 - 2*Sf1124*Ss241*Ss242*UStau11

	Sf1119 = 2*Sf1115 + Sf1110*MStau2(1)*UStau11

	Sf10601 = Sf10245*Sf10246*Sf10599 + Sf10600*MStop2(1)

	Sf10597 = Sf10225*Sf10226*Sf10594 + Sf10596*MStop2(2)

	Sf10588 = Sf10587 + MHiggs2(2)*MStop2(2)

	Sf10589 = Sf10586 + Sf10588*MStop2(2)

	Sf10585 = (2*Sf10578*Sf10584)/MHiggs2(2)**2 - 
     &   Sf10580*Ss988*MStop2(2)

	Sf10579 = (2*Sf10573*Sf10578)/MHiggs2(2)**2 + 
     &   Sf10568*Ss969*MStop2(1)

	Sf10549 = Sf10545 - Sf10548*MStop2(1)

	Sf10542 = Sf10479*Sf10480*Sf10540 - Sf10541*MStop2(1)

	Sf10531 = Sf10528 - Sf10530*MStop2(2)

	Sf10532 = Sf10445*Sf10446*Sf10526 - Sf10531*MStop2(2)

	Sf10512 = Sf10508 + Sf10511*MStop2(1)

	Sf10505 = Sf10225*Sf10226*Sf10500 - Sf10504*MStop2(2)

	Sf10484 = Sf10479*Sf10480*Sf10481 - Sf10483*MStop2(1)

	Sf10450 = Sf10445*Sf10446*Sf10447 - Sf10449*MStop2(2)

	Sf10434 = (2*Sf10428*Sf10433)/MHiggs2(1)**2 + 
     &   Sf10361*Ss969*MStop2(1)

	Sf10429 = (2*Sf10423*Sf10428)/MHiggs2(1)**2 - 
     &   Sf10341*Ss988*MStop2(2)

	Sf10407 = Sf10405 - Sf10406*MStop2(1)

	Sf10408 = Sf10366*Sf10367*Sf10403 - Sf10407*MStop2(1)

	Sf10396 = Sf10393 - Sf10395*MStop2(2)

	Sf10397 = Sf10335*Sf10336*Sf10391 - Sf10396*MStop2(2)

	Sf10382 = Sf10366*Sf10367*Sf10379 - Sf10381*MStop2(1)

	Sf10370 = Sf10366*Sf10367*Sf10368 - Sf10369*MStop2(1)

	Sf10376 = Sf10245*Sf10246*Sf10362 + Sf10375*MStop2(1)

	Sf10355 = Sf10335*Sf10336*Sf10352 - Sf10354*MStop2(2)

	Sf10350 = Sf10225*Sf10226*Sf10342 + Sf10349*MStop2(2)

	Sf10324 = Sf10323 + MHiggs2(1)*MStop2(2)

	Sf10325 = Sf10322 + Sf10324*MStop2(2)

	Sf10261 = MZ2*Sf10258 - 1/2.D0*(Sf10260*MStop2(2))

	Sf10241 = MZ2*Sf10238 - 1/2.D0*(Sf10240*MStop2(1))

	Sf10213 = Sf10212 + MHiggs2(1)*MStau2(2)

	Sf10214 = Sf10211 + Sf10213*MStau2(2)

	Sf10210 = (2*Sf10203*Sf10209)/MHiggs2(1)**2 - 
     &   Sf10205*Ss949*MStau2(2)

	Sf10204 = (2*Sf10198*Sf10203)/MHiggs2(1)**2 + 
     &   Sf10193*Ss928*MStau2(1)

	Sf10188 = Sf10184 - Sf10187*MStau2(1)

	Sf10181 = Sf10099*Sf10100*Sf10179 - Sf10180*MStau2(1)

	Sf10170 = Sf10167 - Sf10169*MStau2(2)

	Sf10171 = Sf10069*Sf10070*Sf10165 - Sf10170*MStau2(2)

	Sf10143 = Sf10139 + Sf10142*MStau2(1)

	Sf10104 = Sf10099*Sf10100*Sf10101 - Sf10103*MStau2(1)

	Sf1033 = MW2**2*Sf1032 + Sf1010*MHiggs2(4)**2

	Sf10074 = Sf10069*Sf10070*Sf10071 - Sf10073*MStau2(2)

	Sf1167 = CB2*SA + SB*Sf1166

	Sf1093 = 6*MW2**2 + Sf1092*MHiggs2(4)

	Sf4 = C2B*MW2 + 2*CW2*Sf3

	Sf2313 = C2B*MW2 - 2*CW2*Sf2312

	Sf120 = C2B*MW2 - 2*CW2*Sf119

	Sf1973 = -3*ML2*Sf1966 + 4*MSl2Diff(2,1)**7 + 
     &   MW2*(Sf1971 + MW2*
     &       (Sf1972 + MW2*
     &          (2*Sf1960 + 
     &            MW2*(-2*Sf1959 + MW2*(Sf1957 - 3*MW2*MSneu2(1))))
     &         ))

	Sf1913 = -3*ML2*Sf1906 + 4*MSl2Diff(1,1)**7 + 
     &   MW2*(-Sf1911 + MW2*
     &       (Sf1912 + MW2*
     &          (2*Sf1898 + 
     &            MW2*(-2*Sf1897 + MW2*(Sf1895 - 3*MW2*MSneu2(1))))
     &         ))

	Sf6905 = 36*C2B*MW2**7 - Sf6904*MSq2Diff(2,2) + 
     &   MW2*(12*Sf6902*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-9*MW2*Sf6889 + Sf6893) + Sf6894) - Sf6896
     &            ) - 3*Sf6901*MSq2Diff(2,2)**2))

	Sf6654 = 36*C2B*MW2**7 - Sf6653*MSq2Diff(1,2) + 
     &   MW2*(12*Sf6651*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-9*MW2*Sf6638 + Sf6642) + Sf6643) - Sf6645
     &            ) - 3*Sf6650*MSq2Diff(1,2)**2))

	Sf5423 = 36*C2B*MW2**7 + Sf5422*MSq2Diff(2,2) + 
     &   MW2*(-12*Sf5421*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(-9*MW2*Sf5409 + Sf5413) - 
     &               9*Sf5414) + Sf5415) + 
     &         3*Sf5420*MSq2Diff(2,2)**2))

	Sf5317 = 36*C2B*MW2**7 + Sf5316*MSq2Diff(1,2) + 
     &   MW2*(-12*Sf5315*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(-9*MW2*Sf5303 + Sf5307) - 
     &               9*Sf5308) + Sf5309) + 
     &         3*Sf5314*MSq2Diff(1,2)**2))

	Sf4709 = 36*C2B*MW2**7 - Sf4708*MSq2Diff(2,1) + 
     &   MW2*(12*Sf4706*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-9*MW2*Sf4684 + Sf4690) + Sf4693) - Sf4700
     &            ) - 3*Sf4705*MSq2Diff(2,1)**2))

	Sf4612 = 36*C2B*MW2**7 - Sf4611*MSq2Diff(1,1) + 
     &   MW2*(12*Sf4609*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-9*MW2*Sf4579 + Sf4588) + Sf4593) - Sf4603
     &            ) - 3*Sf4608*MSq2Diff(1,1)**2))

	Sf3023 = 36*C2B*MW2**7 + Sf3022*MSq2Diff(2,1) + 
     &   MW2*(-12*Sf3021*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(-9*MW2*Sf3001 + Sf3007) - 
     &               9*Sf3011) + Sf3015) + 
     &         3*Sf3020*MSq2Diff(2,1)**2))

	Sf2896 = 36*C2B*MW2**7 + Sf2895*MSq2Diff(1,1) + 
     &   MW2*(-12*Sf2894*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(-9*MW2*Sf2863 + Sf2874) - 
     &               9*Sf2880) + Sf2888) + 
     &         3*Sf2893*MSq2Diff(1,1)**2))

	Sf897 = ML2*Sf896 - MSl2Diff(1,1)**5 + 
     &   MW2*(Sf500*MSl2Diff(1,1)**3 + 
     &      MW2*(MW2*(MW2*(2*MW2*Sf627 + Sf885) - Sf894) - 
     &         Sf889*MSl2Diff(1,1)))

	Sf578 = 3*ML2*Sf572 + 4*MSl2Diff(1,1)**8 + 
     &   MW2*(-Sf577 + MW2*
     &       (5*Sf560*MSl2Diff(1,1)**4 + 
     &         MW2*(MW2*(3*MW2**2*Sf559 + 10*Sf561) + 
     &            6*Sf31*Sf345*MSneu2(1)*MStau2(1))))

	Sf1696 = ML2*Sf1695 - MSl2Diff(2,1)**5 + 
     &   MW2*(Sf1416*MSl2Diff(2,1)**3 + 
     &      MW2*(MW2*(-Sf1693 + MW2*(Sf1689 + 2*MW2*Sf627)) - 
     &         Sf1691*MSl2Diff(2,1)))

	Sf1479 = 3*ML2*Sf1473 + 4*MSl2Diff(2,1)**8 + 
     &   MW2*(-Sf1478 + MW2*
     &       (5*Sf1463*MSl2Diff(2,1)**4 + 
     &         MW2*(MW2*(3*MW2**2*Sf1462 + 10*Sf1464) + 
     &            6*Sf345*Sf44*MSneu2(1)*MStau2(2))))

	Sf9529 = Sf9506 - Sf9514*MStau2(2) + 
     &   MStau2(2)**2*(Sf9521 + 
     &      MStau2(2)*(-Sf9525 + 
     &         MStau2(2)*(Sf9527 - 2*Sf9528*MStau2(2))))

	Sf882 = 3*ML2*Sf881 + Sf869*MSl2Diff(1,1)**3 + 
     &   MW2*(MW2*(MW2*(MW2*Sf868 + Sf871) - Sf880) - 
     &      Sf875*MSl2Diff(1,1))

	tmp52 = 6*MBot2 - 6*MTop2 + MW2 - 6*C2B*MW2

	tmp54 = -3*MBot2 + 3*MTop2 + 23*MW2 + 3*C2B*MW2

	tmp55 = -12*MBot2 + 12*MTop2 - MW2 + 12*C2B*MW2

	Sf5940 = (MW2 - MSq2Diff(2,2))**2*
     &    (MSq2Diff(2,2)**3 + 
     &      MW2*(-13*MSq2Diff(2,2)**2 + 
     &         MW2*(tmp52 + 5*MSq2Diff(2,2)))) + 
     &   2*MW2*(-13*MSq2Diff(2,2)**3 + 
     &      MW2*(MW2*tmp55 - tmp54*MSq2Diff(2,2) + 
     &         46*MSq2Diff(2,2)**2))*MStop2(2) - 
     &   2*MW2**2*(-3*MBot2 + 3*MTop2 + 4*MW2 + 3*C2B*MW2 - 
     &      41*MSq2Diff(2,2))*MStop2(2)**2

	Sf5515 = (MW2 - MSq2Diff(1,2))**2*
     &    (MSq2Diff(1,2)**3 + 
     &      MW2*(-13*MSq2Diff(1,2)**2 + 
     &         MW2*(tmp52 + 5*MSq2Diff(1,2)))) + 
     &   2*MW2*(-13*MSq2Diff(1,2)**3 + 
     &      MW2*(MW2*tmp55 - tmp54*MSq2Diff(1,2) + 
     &         46*MSq2Diff(1,2)**2))*MStop2(2) - 
     &   2*MW2**2*(-3*MBot2 + 3*MTop2 + 4*MW2 + 3*C2B*MW2 - 
     &      41*MSq2Diff(1,2))*MStop2(2)**2

	Sf5983 = 6*Sf3565 + Sf3567*MStop2(2) + 
     &   MStop2(2)**2*(Sf3568 + 
     &      MStop2(2)*(-6*Sf3569 + 
     &         (2*Sf3570 - 2*MStop2(2))*MStop2(2)))

	Sf3571 = 6*Sf3565 + Sf3567*MStop2(1) + 
     &   MStop2(1)**2*(Sf3568 + 
     &      MStop2(1)*(-6*Sf3569 + 
     &         (2*Sf3570 - 2*MStop2(1))*MStop2(1)))

	Sf5557 = 6*Sf3212 + Sf3215*MStop2(2) + 
     &   MStop2(2)**2*(Sf3216 + 
     &      MStop2(2)*(-6*Sf3217 + 
     &         (2*Sf3218 - 2*MStop2(2))*MStop2(2)))

	Sf3219 = 6*Sf3212 + Sf3215*MStop2(1) + 
     &   MStop2(1)**2*(Sf3216 + 
     &      MStop2(1)*(-6*Sf3217 + 
     &         (2*Sf3218 - 2*MStop2(1))*MStop2(1)))

	Sf3504 = Sf2563**2*Sf3136 - Sf3134*MSbot2(2) + 
     &   MSbot2(2)**2*(2*Sf3131 + 
     &      MSbot2(2)*(2*Sf3129 - 5*Sf3128*MSbot2(2) + 
     &         MSbot2(2)**2))

	Sf3137 = Sf2563**2*Sf3136 - Sf3134*MSbot2(1) + 
     &   MSbot2(1)**2*(2*Sf3131 + 
     &      MSbot2(1)*(2*Sf3129 - 5*Sf3128*MSbot2(1) + 
     &         MSbot2(1)**2))

	Sf1688 = 3*ML2*Sf1687 + Sf1681*MSl2Diff(2,1)**3 + 
     &   MW2*(MW2*(-Sf1686 + MW2*(Sf1682 + MW2*Sf868)) - 
     &      Sf1684*MSl2Diff(2,1))

	Sf8760 = -2*MZ2**3*(MW2*(16*MW2 - 24*MZ2) + 11*MZ2**2)*
     &    SW2 - Sf7756*Sf8759*MStop2(1)*UStop12 - 
     &   4*MZ2*Sf7784*MStop2(1)*UStop12**4

	Sf610 = C2B*MW2**5 - 4*MSl2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(-(MW2*Sf600) + Sf605) - Sf609) + 
     &      Sf606*MSl2Diff(1,1)**2)

	Sf1506 = C2B*MW2**5 - 4*MSl2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(-(MW2*Sf1500) + Sf1502) - Sf1505) + 
     &      Sf1503*MSl2Diff(2,1)**2)

	tmp53 = 3*ML2*MW2**2 - 3*MW2**3*Sf966

	Sf1799 = ML2*tmp53 - 20*Sf1798*MStau2(2)**3 + 
     &   MW2*(-2*Sf971*MStau2(2)**2 + 
     &      MW2*(-(MW2**2*Sf967) - Sf969*MStau2(2)))

	Sf975 = ML2*tmp53 - 20*Sf974*MStau2(1)**3 + 
     &   MW2*(-2*Sf971*MStau2(1)**2 + 
     &      MW2*(-(MW2**2*Sf967) - Sf969*MStau2(1)))

	Sf699 = MW2**3*(MW2*(2*MW2*Sf635 - 2*Sf637) + Sf658) - 
     &   Sf695*MSl2Diff(1,1)**3 + 2*Sf636*MSl2Diff(1,1)**5

	Sf1567 = MW2**3*(MW2*(2*MW2*Sf1529 - 2*Sf1531) + 
     &      Sf1538) - Sf1563*MSl2Diff(2,1)**3 + 
     &   2*Sf1530*MSl2Diff(2,1)**5

	Sf7023 = 9*C2B*MW2**4 - Sf6829*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(-(MW2*Sf7018) + Sf7020) + Sf7022*MSq2Diff(2,2))

	Sf7001 = 9*C2B*MW2**4 - Sf6574*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(-(MW2*Sf6996) + Sf6998) + Sf7000*MSq2Diff(1,2))

	Sf4506 = 9*C2B*MW2**4 - Sf4496*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(-(MW2*Sf4497) + Sf4501) + Sf4505*MSq2Diff(2,1))

	Sf4398 = 9*C2B*MW2**4 - Sf4387*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(-(MW2*Sf4388) + Sf4392) + Sf4397*MSq2Diff(1,1))

	Sf6246 = 3*C2B*MW2**4 + Sf6245*MSq2Diff(2,2)**3 + 
     &   MW2*(MW2*(-(MW2*Sf6240) + Sf6244) - Sf6242*MSq2Diff(2,2))

	Sf6135 = 3*C2B*MW2**4 + Sf4197*MSq2Diff(2,1)**3 + 
     &   MW2*(MW2*(-(MW2*Sf6127) + Sf6134) - Sf6131*MSq2Diff(2,1))

	Sf5876 = 3*C2B*MW2**4 + Sf5870*MSq2Diff(1,2)**3 + 
     &   MW2*(MW2*(-(MW2*Sf5871) + Sf5875) - Sf5873*MSq2Diff(1,2))

	Sf5752 = 3*C2B*MW2**4 + Sf4117*MSq2Diff(1,1)**3 + 
     &   MW2*(MW2*(-(MW2*Sf5742) + Sf5751) - Sf5747*MSq2Diff(1,1))

	Sf1812 = 2*ML2*MW2**3 + 30*MStau2(2)**4 + 
     &   MW2*(-58*MStau2(2)**3 + 
     &      MW2*(-(MW2**2*Sf720) + Sf1811*MStau2(2)))

	Sf1001 = 2*ML2*MW2**3 + 30*MStau2(1)**4 + 
     &   MW2*(-58*MStau2(1)**3 + 
     &      MW2*(-(MW2**2*Sf720) + Sf1000*MStau2(1)))

	tmp68 = -2*Sf4290*MSbot2(2)**2 + 
     &   MHiggs2(4)*(12*MHiggs2(4) + Sf4289*MSbot2(2))

	Sf6785 = tmp68*MHiggs2(4) - Sf4291*MSbot2(2)**3 + 
     &   Sf6784*MStop2(2)

	Sf4429 = tmp68*MHiggs2(4) - Sf4291*MSbot2(2)**3 + 
     &   Sf4428*MStop2(1)

	tmp66 = -2*Sf4290*MSbot2(1)**2 + 
     &   MHiggs2(4)*(12*MHiggs2(4) + Sf4289*MSbot2(1))

	Sf6531 = tmp66*MHiggs2(4) - Sf4291*MSbot2(1)**3 + 
     &   Sf6530*MStop2(2)

	Sf4300 = tmp66*MHiggs2(4) - Sf4291*MSbot2(1)**3 + 
     &   Sf4299*MStop2(1)

	tmp67 = -14*Sf2371*MSbot2(2)**2 + 
     &   MHiggs2(4)*(4*MHiggs2(4) + Sf2370*MSbot2(2))

	Sf5063 = tmp67*MHiggs2(4) - 3*Sf2372*MSbot2(2)**3 + 
     &   Sf5062*MStop2(2)

	Sf2612 = tmp67*MHiggs2(4) - 3*Sf2372*MSbot2(2)**3 + 
     &   Sf2611*MStop2(1)

	tmp65 = -14*Sf2371*MSbot2(1)**2 + 
     &   MHiggs2(4)*(4*MHiggs2(4) + Sf2370*MSbot2(1))

	Sf4937 = tmp65*MHiggs2(4) - 3*Sf2372*MSbot2(1)**3 + 
     &   Sf4936*MStop2(2)

	Sf2382 = tmp65*MHiggs2(4) - 3*Sf2372*MSbot2(1)**3 + 
     &   Sf2381*MStop2(1)

	Sf1797 = 87*ML2*MW2**2 + 5*Sf1796*MStau2(2)**2 + 
     &   MW2*(MW2**2*Sf960 + 3*Sf962*MStau2(2))

	Sf965 = 87*ML2*MW2**2 + 5*Sf964*MStau2(1)**2 + 
     &   MW2*(MW2**2*Sf960 + 3*Sf962*MStau2(1))

	Sf6222 = Sf4054*Sf6221 + 
     &   MHiggs2(4)*(-2*Sf6219 + 4*MHiggs2(4)**2) - 
     &   Sf6220*MSbot2(2)*MStop2(2)

	Sf5845 = Sf4054*Sf5844 + 
     &   MHiggs2(4)*(-2*Sf5842 + 4*MHiggs2(4)**2) - 
     &   Sf5843*MSbot2(1)*MStop2(2)

	Sf4163 = Sf4054*Sf4162 + 
     &   MHiggs2(4)*(-2*Sf4160 + 4*MHiggs2(4)**2) - 
     &   Sf4161*MSbot2(2)*MStop2(1)

	Sf4061 = Sf4054*Sf4060 + 
     &   MHiggs2(4)*(-2*Sf4058 + 4*MHiggs2(4)**2) - 
     &   Sf4059*MSbot2(1)*MStop2(1)

	Sf7225 = MZ2*Sf7198*Sf7221*UStop21 - 
     &   Sf86**2*UStop21**4 + 4*Sf7224*SW2*UStop22

	Sf7206 = Sf7198*Sf7203*UStop11 + 
     &   4*MZ2*Sf7205*SW2*UStop12 - Sf86**2*UStop22**4

	Sf1782 = MZ2*Sf1754*Sf1778*UStau21 - 
     &   Ss443**2*UStau21**4 + 2*Sf1781*SW2*UStau22

	Sf1762 = Sf1754*Sf1759*UStau11 + 
     &   2*MZ2*Sf1761*SW2*UStau12 - Ss443**2*UStau22**4

	Sf4801 = MZ2*Sf4772*Sf4797*USbot21 + 
     &   Sf15**2*USbot21**4 - 2*Sf4800*SW2*USbot22

	Sf4780 = Sf4772*Sf4777*USbot11 - 
     &   2*MZ2*Sf4779*SW2*USbot12 + Sf15**2*USbot22**4

	Sf763 = 6*Sf761*MSneu2(1)**2 - Sf752*MSneu2(1)**3 + 
     &   Sf762*MStau2(1)

	Sf631 = 8*MSl2Diff(1,1)**3 + 12*C2B*MW2**2*MStau2(1) + 
     &   6*MW2*MSl2Diff(1,1)*
     &    ((-1 + 5*C2B)*MSneu2(1) - MStau2(1) - 5*C2B*MStau2(1)) + 
     &   6*ML2*(5*MSl2Diff(1,1)**2 - 2*MW2*MStau2(1))

	Sf305 = Sf302*MSneu2(1)**2 - Sf299*MSneu2(1)**3 + 
     &   Sf304*MStau2(1)

	Sf1608 = 6*Sf1606*MSneu2(1)**2 - Sf752*MSneu2(1)**3 + 
     &   Sf1607*MStau2(2)

	Sf1525 = 8*MSl2Diff(2,1)**3 + 12*C2B*MW2**2*MStau2(2) + 
     &   6*MW2*MSl2Diff(2,1)*
     &    ((-1 + 5*C2B)*MSneu2(1) - MStau2(2) - 5*C2B*MStau2(2)) + 
     &   6*ML2*(5*MSl2Diff(2,1)**2 - 2*MW2*MStau2(2))

	Sf1265 = Sf1262*MSneu2(1)**2 - Sf299*MSneu2(1)**3 + 
     &   Sf1264*MStau2(2)

	Sf5217 = Sf5204*MSbot2(2) + Sf5198*Ss25*MSbot2(2) - 
     &   Sf5216*MStop2(2)

	Sf5213 = Sf5184*MSbot2(2) + Sf5180*Ss25*MSbot2(2) - 
     &   Sf4212*MStop2(1)

	Sf6413 = Sf4866*MSbot2(2) + Sf4869*Ss25*MSbot2(2) - 
     &   Sf6412*MStop2(2)

	Sf2293 = 2*Sf2292*MSbot2(2) - Sf2291*USbot21 - 
     &   6*MBot*Sf2225*MHiggs2(4)**2*USbot22

	Sf2290 = Sf2289*MSbot2(2) + Sf2288*USbot21 + 
     &   3*MBot*Sf2225*MHiggs2(4)**2*USbot22

	Sf6409 = Sf2269*MSbot2(2) + Sf2272*Ss25*MSbot2(2) - 
     &   Sf4754*MStop2(1)

	Sf5215 = Sf5192*MSbot2(1) + Sf5186*Ss25*MSbot2(1) - 
     &   Sf5214*MStop2(2)

	Sf5212 = Sf5178*MSbot2(1) + Sf5174*Ss25*MSbot2(1) - 
     &   Sf4209*MStop2(1)

	Sf6411 = Sf4830*MSbot2(1) + Sf4834*Ss25*MSbot2(1) - 
     &   Sf6410*MStop2(2)

	Sf2241 = 2*Sf2240*MSbot2(1) - Sf2237*USbot11 - 
     &   6*MBot*Sf2225*MHiggs2(4)**2*USbot12

	Sf2234 = Sf2233*MSbot2(1) + Sf2229*USbot11 + 
     &   3*MBot*Sf2225*MHiggs2(4)**2*USbot12

	Sf6408 = Sf2205*MSbot2(1) + Sf2208*Ss25*MSbot2(1) - 
     &   Sf4746*MStop2(1)

	Sf634 = (12*MW2*Sf162*Sf633*TB)/MHiggs2(4)**2 + 
     &   Sf631*UStau11 - 6*ML*Sf215*Sf632*UStau12

	Sf1528 = (12*MW2*Sf1191*Sf1527*TB)/MHiggs2(4)**2 + 
     &   Sf1525*UStau21 - 6*ML*Sf1526*Sf215*UStau22

	Sf7294 = Sf1749*Sf5227*Ss25 - Sf1743*Sf5230*Ss25 - 
     &   3*MW2**2*Sf5231*TB*MHiggs2(4)

	Sf6004 = 63*MBot2*Sf99**2 - 63*MTop2*Sf99**2 - Sf6003*MSq2Diff(2,2)

	Sf5948 = 15*MBot2*Sf99**2 - 15*MTop2*Sf99**2 - Sf5947*MSq2Diff(2,2)

	Sf9625 = 1/2.D0*(Sf9609*Sf9624) + 
     &   Sf9017*Sf9611*MStau2(1) - Sf9612*Sf9616*MStau2(2)

	Sf9610 = 1/2.D0*(Sf9606*Sf9609) - 
     &   Sf9593*Sf9597*MStau2(1) + Sf9017*Sf9598*MStau2(2)

	Sf9241 = 1/2.D0*(Sf9237*Sf9240) - 
     &   Sf9222*Sf9226*MSbot2(1) + Sf9227*Sf9228*MSbot2(2)

	Sf5627 = 9*MBot2*Sf91**2 - 9*MTop2*Sf91**2 + Sf5626*MSq2Diff(1,2)

	Sf9051 = 1/2.D0*(Sf9035*Sf9050) + 
     &   Sf9017*Sf9037*MSbot2(1) - Sf9038*Sf9042*MSbot2(2)

	Sf9036 = 1/2.D0*(Sf9032*Sf9035) - 
     &   Sf9018*Sf9022*MSbot2(1) + Sf9017*Sf9023*MSbot2(2)

	Sf7786 = -(Sf7756*Sf7785*MStau2(1)*UStau12) + 
     &   MZ2*(MZ2*Sf7783 - 2*Sf7784*MStau2(1)*UStau12**4)

	Sf8488 = -(Sf7756*Sf8487*MSbot2(1)*USbot12) + 
     &   MZ2*(MZ2*Sf7783 - 2*Sf7784*MSbot2(1)*USbot12**4)

	Sf3601 = 63*MBot2*Sf77**2 - 63*MTop2*Sf77**2 - Sf3600*MSq2Diff(2,1)

	Sf3516 = 15*MBot2*Sf77**2 - 15*MTop2*Sf77**2 - Sf3515*MSq2Diff(2,1)

	Sf3274 = 63*MBot2*Sf69**2 - 63*MTop2*Sf69**2 - Sf3273*MSq2Diff(1,1)

	Sf3154 = 15*MBot2*Sf69**2 - 15*MTop2*Sf69**2 - Sf3153*MSq2Diff(1,1)

	Sf6202 = MW2*Sf6200 - 2*Sf6198*MSbot2(2)**2 + 
     &   Sf6201*MStop2(2)**2

	Sf617 = Sf615 - Sf612*MSneu2(1) + MSneu2(1)**2

	Sf6057 = MW2*Sf6056 - 4*Sf5099*MSq2Diff(2,2)**2 - 
     &   4*C2B*MW2**2*MStop2(2)

	Sf5824 = MW2*Sf5822 - 2*Sf5820*MSbot2(1)**2 + 
     &   Sf5823*MStop2(2)**2

	Sf5463 = MW2*Sf5462 - 4*Sf4973*MSq2Diff(1,2)**2 - 
     &   4*C2B*MW2**2*MStop2(2)

	Sf7277 = Sf1749*Sf5227 - Sf1743*Sf5230 - 6*MW2*Sf5231*Ss25*TB

	Sf4136 = MW2*Sf4134 - 2*Sf4129*MSbot2(2)**2 + 
     &   Sf4135*MStop2(1)**2

	Sf4017 = MW2*Sf4015 - 2*Sf4009*MSbot2(1)**2 + 
     &   Sf4016*MStop2(1)**2

	Sf3726 = MW2*Sf3725 - 4*Sf2663*MSq2Diff(2,1)**2 - 
     &   4*C2B*MW2**2*MStop2(1)

	Sf6049 = Sf3715 - Sf3716*MStop2(2) + MStop2(2)**2

	Sf3717 = Sf3715 - Sf3716*MStop2(1) + MStop2(1)**2

	Sf3455 = MW2*Sf3454 - 4*Sf2458*MSq2Diff(1,1)**2 - 
     &   4*C2B*MW2**2*MStop2(1)

	Sf5630 = Sf3443 - Sf3444*MStop2(2) + MStop2(2)**2

	Sf3445 = Sf3443 - Sf3444*MStop2(1) + MStop2(1)**2

	Sf2176 = MW2*S2B*Sf2174 - 2*ML2*TB + 
     &   (2*Sf2175)/(Sf1184*Sf154*MHiggs2(4)**2)

	Sf1512 = Sf1510 - Sf1508*MSneu2(1) + MSneu2(1)**2

	Sf10536 = 1/2.D0*(Sf10532*Sf10535) - 
     &   Sf10522*Sf9017*MStop2(1) + Sf10523*Sf10524*MStop2(2)

	Sf10175 = 1/2.D0*(Sf10171*Sf10174) - 
     &   Sf10161*Sf9222*MStau2(1) + Sf10162*Sf10163*MStau2(2)

	Sf1168 = 2*CA2*SB2 + CA*Sf1164 + SA*Sf1167

	Sf1793 = 6*ML2*MW2 - 6*C2B*MW2**2 - Sf1792*MStau2(2)

	Sf953 = 6*ML2*MW2 - 6*C2B*MW2**2 - Sf952*MStau2(1)

	Sf7318 = (MW2*(8*MW2 - 4*Sf7317) + 4*MZ2*Sf7311*SW2)/MW2**2

	Sf8816 = 4 - Re(Sf8815) - 
     &   2*Re((Sf111*Log(-(1/2.D0*Sf106/(MStop(1)*MStop(2))))*
     &        MStop(1)*MStop(2))/MZ2)

	Sf8812 = 4 - Re(Sf8811) - 
     &   2*Re((Sf103*Log(-(1/2.D0*Sf98/(MSbot(2)*MStop(2))))*
     &        MSbot(2)*MStop(2))/MW2)

	Sf8806 = 4 - Re(Sf8805) - 
     &   2*Re((Sf95*Log(-(1/2.D0*Sf90/(MSbot(1)*MStop(2))))*
     &        MSbot(1)*MStop(2))/MW2)

	Sf8799 = 4 - Re(Sf8798) - 
     &   2*Re((Sf81*Log(-(1/2.D0*Sf76/(MSbot(2)*MStop(1))))*
     &        MSbot(2)*MStop(1))/MW2)

	Sf8793 = 4 - Re(Sf8792) - 
     &   2*Re((Sf73*Log(-(1/2.D0*Sf68/(MSbot(1)*MStop(1))))*
     &        MSbot(1)*MStop(1))/MW2)

	Sf8786 = 4 - Re(Sf8785) - 
     &   2*Re((Sf59*Log(-(1/2.D0*Sf54/(MStau(1)*MStau(2))))*
     &        MStau(1)*MStau(2))/MZ2)

	Sf8782 = 4 - Re(Sf8781) - 
     &   2*Re((Sf51*Log(-(1/2.D0*Sf46/(MSneu(1)*MStau(2))))*
     &        MSneu(1)*MStau(2))/MW2)

	Sf8775 = 4 - Re(Sf8774) - 
     &   2*Re((Sf38*Log(-(1/2.D0*Sf33/(MSneu(1)*MStau(1))))*
     &        MSneu(1)*MStau(1))/MW2)

	Sf8768 = 4 - Re(Sf8767) - 
     &   2*Re((Sf24*Log(-(1/2.D0*Sf19/(MSbot(1)*MSbot(2))))*
     &        MSbot(1)*MSbot(2))/MZ2)

	Sf2005 = MW2*Sf2004 + 2*Sf2002*MStau2(2) + 
     &   (Sf2001*Ss241*Ss242 - Sf1944*MSneu2(1))*UStau21 - 
     &   2*ML*Sf215*MHiggs2(4)*MSneu2(1)*UStau22

	Sf1956 = MW2*Sf1955 + 2*Sf1949*MStau2(1) + 
     &   (Sf1945*Ss241*Ss242 - Sf1944*MSneu2(1))*UStau11 - 
     &   2*ML*Sf215*MHiggs2(4)*MSneu2(1)*UStau12

	Sf179 = 4*ML2*MW2**3 + 
     &   MW2**2*(-(MW2**2*Sf168) - Sf178) + 
     &   MW2*Sf170*MHiggs2(4)**2 + Sf173*Sf177*MStau2(1)

	Sf1198 = 4*ML2*MW2**3 + 
     &   MW2**2*(-Sf1197 - MW2**2*Sf168) + 
     &   MW2*Sf170*MHiggs2(4)**2 + Sf173*Sf177*MStau2(2)

	Sf8641 = 3*Sf8635 + 
     &   MZ2**2*(MZ2**2*(MZ2*Sf8637 - 6*Sf8640) + 
     &      12*Sf8534*Sf8597*SW2**2*UStau12**4)

	Sf8559 = MZ2*(MZ2*(-(MZ2*Sf8554) - 42*Sf8558) + 
     &      84*Sf8548*SW2*UStau12) + 
     &   21*Sf8534*Ss443**2*UStau22**4

	Sf8595 = MZ2*(MZ2*(MZ2*Sf8589 - 6*Sf8594) + 
     &      12*Sf8548*SW2*UStau12) + 
     &   3*Sf8534*Ss443**2*UStau22**4

	Sf249 = Sf248*MHiggs2(4)**2 + 
     &   MHiggs2(4)*(-(Sf247*MSneu2(1)) + 
     &      Sf240*MHiggs2(4)**2*MSneu2(1)) + Sf237*Sf238*UStau11

	Sf1231 = Sf1230*MHiggs2(4)**2 + 
     &   MHiggs2(4)*(-(Sf1229*MSneu2(1)) + 
     &      Sf1226*MHiggs2(4)**2*MSneu2(1)) + 
     &   Sf1225*Sf237*UStau21

	Sf7435 = Sf7433/(MZ2**4*SW2**2) - 8/MStop2(2) - 
     &   (9*((CW2*Sf7420)/SW2**2 + 
     &        Sf7434/(SW2*MSq2Diff(1,2)*MSq2Diff(2,2)))*
     &      UStop21**2)/MW2**2

	Sf7419 = Sf7417/(MZ2**4*SW2**2) - 8/MStop2(1) - 
     &   (9*((CW2*Sf7399)/SW2**2 + 
     &        Sf7418/(SW2*MSq2Diff(1,1)*MSq2Diff(2,1)))*
     &      UStop22**2)/MW2**2

	Sf7391 = Sf7385/(MZ2**4*SW2**2) - 2/MSbot2(2) - 
     &   (9*((CW2*Sf7386)/SW2**2 + 
     &        Sf7390/(SW2*MSq2Diff(2,1)*MSq2Diff(2,2)))*
     &      USbot21**2)/MW2**2

	Sf7375 = Sf7365/(MZ2**4*SW2**2) - 2/MSbot2(1) - 
     &   (9*((CW2*Sf7367)/SW2**2 + 
     &        Sf7374/(SW2*MSq2Diff(1,1)*MSq2Diff(1,2)))*
     &      USbot22**2)/MW2**2

	Sf7350 = Sf7348/(MZ2**4*SW2**2) - 2/MStau2(2) + 
     &   ((-((CW2*Sf1768)/SW2**2) + Sf7349/(SW2*MSl2Diff(2,1)))*
     &      UStau21**2)/MW2**2

	Sf7338 = Sf7336/(MZ2**4*SW2**2) - 2/MStau2(1) + 
     &   ((-((CW2*Sf1738)/SW2**2) + Sf7337/(SW2*MSl2Diff(1,1)))*
     &      UStau22**2)/MW2**2

	Sf6977 = 3*CA2*MTop*Sf6975 + SA*Sf6976 + 
     &   3*CB*MTop*Re(Atop)*UStop21 - 
     &   4*CA*MZ2*SB*Sf127*SW2*UStop22

	Sf4260 = 3*CA2*MTop*Sf4256 + SA*Sf4259 + 
     &   3*CB*MTop*Re(Atop)*UStop11 - 
     &   4*CA*MZ2*SB*Sf127*SW2*UStop12

	tmp62 = MStop2(1)**2 - MStop2(2)**2

	tmp63 = 7*MStop2(2)**2 + MStop2(1)*(11*MStop2(1) + 30*MStop2(2))

	tmp64 = 6*MStop2(2)**3 + 
     &   MZ2*(3*MZ2**2 + MStop2(1)*
     &       (-11*MStop2(1) - 30*MStop2(2)) - 7*MStop2(2)**2) + 
     &   MStop2(1)*(6*MStop2(1)**2 - 6*MStop2(2)*MStop2(4))

	tmp71 = 3*MZ2**2 - 5*MZ2*MStop2(2) - 
     &   3*MStop2(1)*MStop2(2) + 3*MStop2(2)**2 + 
     &   MStop2(1)*(-7*MZ2 - 3*MStop2(3))

	tmp72 = -(MZ2**2*tmp63) + 
     &   MZ2**3*(7*MStop2(1) + 5*MStop2(2)) + 3*MStop2(3)**4

	Sf10320 = tmp72*MHiggs2(1) - tmp64*MHiggs2(1)**2 + 
     &   tmp71*MHiggs2(1)**3 + 4*MZ2**3*MStop2(1)*MStop2(3) - 
     &   4*MZ2*MStop2(3)*(MZ2*(2*tmp62 + MZ2*MStop2(2)) + 
     &      MStop2(3)**3)

	Sf10311 = Sf10304*Sf9191 - Sf9182*MStop2(1) + 
     &   Sf9183*MStop2(1)**2 - Sf9184*MStop2(1)**3 + 
     &   Sf10310*MStop2(2) - 2*Sf10303*MStop2(2)**3 + 
     &   MHiggs2(1)**2*(3*MZ2**2*Sf9131*Sf9132 - 
     &      3*Sf10300*MHiggs2(1)*MStop2(2))

	Sf9074 = Sf9057*Sf9067 + Sf9073*MSbot2(1) - 
     &   2*Sf9066*MSbot2(1)**3 + 
     &   MHiggs2(2)**2*(3*MZ2**2*Sf8929*Sf8930 + 
     &      3*Sf8983*MHiggs2(2)*MSbot2(1)) - Sf9052*MSbot2(2) + 
     &   Sf9053*MSbot2(2)**2 - Sf9054*MSbot2(2)**3

	Sf9065 = -5*MZ2*MHiggs2(2)*(-MZ2 + MHiggs2(2))*
     &    (MZ2 + MHiggs2(2))*MSbot2(1) + 
     &   (3*MHiggs2(2)**3 + 
     &      MZ2*(MZ2*(-4*MZ2 - 7*MHiggs2(2)) + 7*MHiggs2(2)**2))*
     &    MSbot2(1)**2 - 2*(-4*MZ2**2 + 3*MHiggs2(2)**2)*
     &    MSbot2(1)**3 - 2*
     &    (3*MHiggs2(2)**2 + 6*MHiggs2(2)*MSbot2(1) - 
     &      4*MZ2*(MZ2 + 2*MSbot2(1)))*MSbot2(2)**3 + 
     &   (-4*MZ2 + 3*MHiggs2(2))*
     &    (6*MSbot2(1)**2*MSbot2(2)**2 + MSbot2(2)**4 - 
     &      MSbot2(1)**3*(3*MSbot2(2) + MSbot2(3))) + 
     &   3*MHiggs2(2)**2*(-MZ2**3 + 
     &      MHiggs2(2)*(MZ2**2 + 
     &         MSbot2(2)*(-MSbot2(1) + MSbot2(3)))) + 
     &   MSbot2(2)*(MZ2**2*
     &       (MHiggs2(2)*(7*MZ2 - 30*MSbot2(1) - 11*MSbot2(2)) - 
     &         4*MZ2*MSbot2(3) + 4*MSbot2(1)*(MZ2 - 2*MSbot2(4)))-
     &        MHiggs2(2)**2*
     &       (MZ2*(7*MHiggs2(2) - 30*MSbot2(1) - 11*MSbot2(2)) - 
     &         6*MSbot2(1)*MSbot2(4)))

	Sf10567 = tmp72*MHiggs2(2) - tmp64*MHiggs2(2)**2 + 
     &   tmp71*MHiggs2(2)**3 + 4*MZ2**3*MStop2(1)*MStop2(3) - 
     &   4*MZ2*MStop2(3)*(MZ2*(2*tmp62 + MZ2*MStop2(2)) + 
     &      MStop2(3)**3)

	Sf10559 = Sf10304*Sf9057 - Sf9052*MStop2(1) + 
     &   Sf9053*MStop2(1)**2 - Sf9054*MStop2(1)**3 + 
     &   Sf10558*MStop2(2) - 2*Sf10552*MStop2(2)**3 + 
     &   MHiggs2(2)**2*(3*MZ2**2*Sf8929*Sf8930 - 
     &      3*Sf10300*MHiggs2(2)*MStop2(2))

	Sf6363 = Sf126*TB*MSbot2(1) + Sf139*TB*MSbot2(2) - 
     &   Sf7*TB2*MSneu2(1) - Sf153*TB2*MStau2(1) - 
     &   Sf1183*TB2*MStau2(2) + Sf2319*TB*MStop2(1) + 
     &   Sf4899*TB*MStop2(2)

	Sf5209 = Sf126*MSbot2(1) + Sf139*MSbot2(2) - 
     &   Sf7*TB*MSneu2(1) - Sf153*TB*MStau2(1) - 
     &   Sf1183*TB*MStau2(2) + Sf2319*MStop2(1) + Sf4899*MStop2(2)

	Sf6404 = Sf5209 + 
     &   CW2*(-2*Sf6394*Ss928*TB*MStau2(1) - 
     &      2*Sf6395*Ss949*TB*MStau2(2))

	Sf9793 = Sf9792*MStau2(1)**2 - Sf9786*MStau2(1)**3 + 
     &   Sf9783*MStau2(1)**4 + 2*Sf9777*MStau2(1)**5 - 
     &   2*Sf9775*MStau2(1)**6 + 2*Sf9583*MStau2(1)**7 + 
     &   Sf9789*MStau2(1)*MStau2(2) - 
     &   2*MZ2*Sf9302*Sf9303*Sf9363*Sf9364*Sf9780*MHiggs2(2)*
     &    MStau2(2)**2

	Sf9765 = Sf9764*MStau2(1)**2 - Sf9758*MStau2(1)**3 + 
     &   Sf9755*MStau2(1)**4 - Sf9752*MStau2(1)**5 + 
     &   Sf9748*MStau2(1)**6 - 2*Sf9528*MStau2(1)**7 + 
     &   MZ2*Sf9302*Sf9303*Sf9363*Sf9364*Sf9750*MHiggs2(2)*
     &    MStau2(2)**2 - Sf9761*MStau2(1)*MStau2(2)**2

	Sf9736 = Sf9732*MStau2(1)**2 - Sf9729*MStau2(1)**3 + 
     &   Sf9726*MStau2(1)**4 - 2*Sf9723*MStau2(1)**5 + 
     &   Sf9719*MStau2(1)**6 - 2*Sf9469*MStau2(1)**7 + 
     &   Sf9735*MStau2(1)*MStau2(2) + 
     &   8*MZ2*Sf9302*Sf9303*Sf9363*Sf9364*Sf9721*MHiggs2(2)*
     &    MStau2(2)**2

	Sf7274 = Sf113*Sf117**2 + Sf11*Sf16**2 + Sf26*Sf30**2 + 
     &   3*Sf39*Sf43**2 + 3*Sf61*Sf65**2 + Sf82*Sf87**2 - 
     &   2*Sf4780*MSbot2(1) - 2*Sf4801*MSbot2(2) + 
     &   6*Sf1762*MStau2(1) + 6*Sf1782*MStau2(2) + 
     &   2*Sf7206*MStop2(1) + 2*Sf7225*MStop2(2)

	Sf6814 = MHiggs2(4)**2*(-3*Sf6813 + Sf6810*MHiggs2(4)**2)

	Sf6559 = MHiggs2(4)**2*(-3*Sf6558 + Sf6555*MHiggs2(4)**2)

	Sf4475 = MHiggs2(4)**2*(-3*Sf4474 + Sf4470*MHiggs2(4)**2)

	Sf4364 = MHiggs2(4)**2*(-3*Sf4363 + Sf4355*MHiggs2(4)**2)

	Sf611 = MHiggs2(4)*(-(MW2*Sf610) - 2*Sf599*MHiggs2(4)) + 
     &   MW2**2*(Sf584*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf585 - Sf587) + 3*Sf594) - 
     &         Sf591*MSl2Diff(1,1)**2))

	Sf583 = 5*MSl2Diff(1,1)**7 + 
     &   MW2*(-(Sf565*MSl2Diff(1,1)**5) + 
     &      MW2*(10*Sf568*MSl2Diff(1,1)**3 + 
     &         MW2*(2*Sf580 + 
     &            MW2*(-Sf582 + MW2*(Sf567 + 2*MW2*MStau2(1))))))

	Sf429 = 3*MSl2Diff(1,1)**6 + 
     &   MW2*(-(Sf411*MSl2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf408 - 2*Sf414) + Sf424) - 
     &            Sf428) + Sf419*MSl2Diff(1,1)**2))

	Sf5395 = 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf5387*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2827 - Sf5388) + Sf5392) - 
     &            Sf5394) + Sf5390*MSq2Diff(2,2)**2))

	Sf5289 = 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf5279*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2827 - Sf5281) + Sf5286) - 
     &            Sf5288) + Sf5284*MSq2Diff(1,2)**2))

	Sf2987 = 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf2979*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2827 - Sf2980) + Sf2984) - 
     &            Sf2986) + Sf2982*MSq2Diff(2,1)**2))

	Sf2849 = 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf2830*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf2827 - Sf2833) + Sf2843) - 
     &            Sf2848) + Sf2838*MSq2Diff(1,1)**2))

	Sf2131 = MHiggs2(4)*
     &    (-(MW2*Sf2122) - Sf2130*MHiggs2(4)) + 
     &   MW2**2*(-MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(-Sf2114 + MW2*(Sf2108 - MW2*Sf720)) - 
     &         Sf2111*MSl2Diff(2,1)))

	Sf2051 = MHiggs2(4)*
     &    (-(MW2*Sf2040) - Sf2050*MHiggs2(4)) + 
     &   MW2**2*(-MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(-Sf2030 + MW2*(Sf2024 - MW2*Sf720)) - 
     &         Sf2027*MSl2Diff(1,1)))

	Sf730 = MHiggs2(4)*(MW2*Sf719 - Sf729*MHiggs2(4)) + 
     &   MW2**2*(MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf184 + Sf705) - Sf709) + 
     &         Sf707*MSl2Diff(1,1)))

	Sf1591 = MHiggs2(4)*
     &    (-(MW2*Sf1582) + Sf1590*MHiggs2(4)) + 
     &   MW2**2*(MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(-Sf1575 + MW2*(Sf1571 + MW2*Sf184)) + 
     &         Sf1573*MSl2Diff(2,1)))

	Sf1507 = MHiggs2(4)*
     &    (-(MW2*Sf1506) - 2*Sf1499*MHiggs2(4)) + 
     &   MW2**2*(Sf1485*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1486 - Sf1488) + 3*Sf1494) - 
     &         Sf1491*MSl2Diff(2,1)**2))

	Sf1484 = 5*MSl2Diff(2,1)**7 + 
     &   MW2*(-(Sf1468*MSl2Diff(2,1)**5) + 
     &      MW2*(10*Sf1470*MSl2Diff(2,1)**3 + 
     &         MW2*(2*Sf1481 + 
     &            MW2*(-Sf1483 + MW2*(Sf1469 + 2*MW2*MStau2(2))))))

	Sf1352 = -3*MSl2Diff(2,1)**6 + 
     &   MW2*(Sf1345*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(Sf1351 + 
     &            MW2*(-Sf1350 + MW2*(2*Sf1346 - MW2*Sf408))) - 
     &         Sf1348*MSl2Diff(2,1)**2))

	Sf6247 = MHiggs2(4)*(MW2*Sf6246 + Sf6237*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf6025*MSq2Diff(2,2)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1091 + Sf6231) - Sf6239) + 
     &         Sf6233*MSq2Diff(2,2)))

	Sf6136 = MHiggs2(4)*(MW2*Sf6135 + Sf6124*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf3679*MSq2Diff(2,1)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1091 + Sf6111) - Sf6126) + 
     &         Sf6114*MSq2Diff(2,1)))

	Sf5883 = MHiggs2(4)*(MW2*Sf5876 + Sf5882*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf5605*MSq2Diff(1,2)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1091 + Sf5865) - Sf5869) + 
     &         Sf5867*MSq2Diff(1,2)))

	Sf5753 = MHiggs2(4)*(MW2*Sf5752 + Sf5737*MHiggs2(4)) + 
     &   MW2**2*(-3*Sf3399*MSq2Diff(1,1)**3 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1091 + Sf5722) - Sf5740) + 
     &         Sf5726*MSq2Diff(1,1)))

	Sf834 = 162*MSl2Diff(1,1)**6 + 
     &   MW2*(-27*Sf815*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(54*MW2 - 6*Sf824) + Sf828) - 
     &            3*Sf833) + 3*Sf819*MSl2Diff(1,1)**2))

	Sf1651 = 162*MSl2Diff(2,1)**6 + 
     &   MW2*(-27*Sf1640*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(54*MW2 - 6*Sf1647) + Sf1649) - 
     &            3*Sf1650) + 3*Sf1644*MSl2Diff(2,1)**2))

	Sf8610 = MW2*(Sf8605 + 
     &      2*MW2*(MZ2**3*Sf8599 + 
     &         MW2*(MW2*(16*MW2*Sf8540 - 8*MZ2*Sf8601) - 
     &            2*MZ2**2*Sf8609))) + 
     &   15*MZ2**5*MSneu2(1)*UStau11

	Sf5949 = 3*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(-Sf5948 + 
     &         MW2*(Sf5945 + MW2*(Sf5941 + MW2*Sf601))) + 
     &      Sf5943*MSq2Diff(2,2)**2)

	Sf5523 = 3*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(-Sf5522 + 
     &         MW2*(Sf5520 + MW2*(Sf5516 + MW2*Sf601))) + 
     &      Sf5518*MSq2Diff(1,2)**2)

	Sf3517 = 3*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(-Sf3516 + 
     &         MW2*(Sf3513 + MW2*(Sf3505 + MW2*Sf601))) + 
     &      Sf3509*MSq2Diff(2,1)**2)

	Sf6047 = 22*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2*Sf3422 + 2*Sf6040) - 2*Sf6044) + 
     &         2*Sf6046) - 2*Sf6042*MSq2Diff(2,2)**2)

	Sf5628 = 22*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2*Sf3422 + 2*Sf5620) - 2*Sf5624) + 
     &         2*Sf5627) - 2*Sf5622*MSq2Diff(1,2)**2)

	Sf3711 = 22*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2*Sf3422 + 2*Sf3700) - 2*Sf3708) + 
     &         2*Sf3710) - 2*Sf3704*MSq2Diff(2,1)**2)

	Sf3439 = 22*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(MW2*Sf3422 + 2*Sf3424) - 2*Sf3434) + 
     &         2*Sf3438) - 2*Sf3428*MSq2Diff(1,1)**2)

	end

