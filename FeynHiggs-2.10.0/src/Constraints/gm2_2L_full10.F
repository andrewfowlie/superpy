* gm2_2L_full10.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full10
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp100, tmp101, tmp102, tmp103, tmp104, tmp105
	ComplexType tmp106, tmp107, tmp108, tmp109, tmp110, tmp111
	ComplexType tmp112, tmp113, tmp114, tmp115, tmp116, tmp117
	ComplexType tmp118, tmp119, tmp120, tmp121, tmp98, tmp99

	Sf5333 = -8*MSq2Diff(1,2)**6 + 
     &   MW2*(4*Sf5302*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf5318) - Sf5327) + 
     &            Sf5332) - Sf5329*MSq2Diff(1,2)))

	Sf3047 = -8*MSq2Diff(2,1)**6 + 
     &   MW2*(4*Sf3000*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf3024) - Sf3038) + 
     &            Sf3046) - Sf3042*MSq2Diff(2,1)))

	Sf2920 = -8*MSq2Diff(1,1)**6 + 
     &   MW2*(4*Sf2862*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf2897) - Sf2911) + 
     &            Sf2919) - Sf2915*MSq2Diff(1,1)))

	Sf1859 = 10*MSl2Diff(2,1)**6 + 
     &   MW2*(-(Sf1848*MSl2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(2*MW2 - 2*Sf1853) + Sf1855) - 
     &            Sf1858) + Sf1850*MSl2Diff(2,1)**2))

	Sf1086 = 10*MSl2Diff(1,1)**6 + 
     &   MW2*(-(Sf1069*MSl2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(2*MW2 - 2*Sf1074) + Sf1078) - 
     &            Sf1085) + Sf1071*MSl2Diff(1,1)**2))

	Sf8871 = 6*Sf7756*(Sf8858 + MZ2*Sf8859) + 
     &   MZ2**3*(9*Sf8860 + 
     &      MZ2*(MZ2*(6*Sf8865 - MZ2*Sf8870) + 
     &         5796*Sf8820*Ss443*SW2*UStau11*UStau12))

	Sf8642 = 6*Sf8612 + 
     &   MW2*(MW2*(MW2*(4*MZ2**2*Sf8620 + 
     &            MW2*(32*MW2*Sf8587 - 32*MZ2*Sf8628)) - 
     &         4*MZ2**3*Sf8634) + Sf8641)

	Sf1999 = MW2*(MW2*
     &       (MW2*(MW2*(MW2*Sf1047 - Sf1837) + Sf1995) - Sf1998) - 
     &      Sf1840*MSl2Diff(2,1)) - 
     &   4*Sf1836*MSl2Diff(2,1)**3*MStau2(2)

	Sf1941 = MW2*(MW2*
     &       (MW2*(MW2*(MW2*Sf1047 - Sf1051) + Sf1937) - Sf1940) - 
     &      Sf1059*MSl2Diff(1,1)) - 
     &   4*Sf1048*MSl2Diff(1,1)**3*MStau2(1)

	Sf6160 = 24*MSq2Diff(2,2)**6 + 
     &   MW2*(MW2**2*(MW2*(MW2*(4*MW2 - Sf6154) + Sf6156) - 
     &         Sf6159) - 108*MSbot2(2)*MSq2Diff(2,2)**4)

	Sf5781 = 24*MSq2Diff(1,2)**6 + 
     &   MW2*(MW2**2*(MW2*(MW2*(4*MW2 - Sf5775) + Sf5777) - 
     &         Sf5780) - 108*MSbot2(1)*MSq2Diff(1,2)**4)

	Sf3912 = 24*MSq2Diff(2,1)**6 + 
     &   MW2*(MW2**2*(MW2*(MW2*(4*MW2 - Sf3906) + Sf3908) - 
     &         Sf3911) - 108*MSbot2(2)*MSq2Diff(2,1)**4)

	Sf3784 = 24*MSq2Diff(1,1)**6 + 
     &   MW2*(MW2**2*(MW2*(MW2*(4*MW2 - Sf3772) + Sf3776) - 
     &         Sf3783) - 108*MSbot2(1)*MSq2Diff(1,1)**4)

	Sf3656 = MW2**3*Sf3655*MSbot2(2) + 
     &   (Sf3654*MSbot2(2)**5 + 
     &      MW2*(-(Sf3649*MSbot2(2)**4) + 
     &         MW2*(-4*MW2**3*Sf2731 + Sf3650*MSbot2(2)**3)))*
     &    USbot21

	Sf3647 = MW2**3*Sf3646*MSbot2(2) - 
     &   (Sf3645*MSbot2(2)**5 + 
     &      MW2*(-35*Sf3640*MSbot2(2)**4 + 
     &         MW2*(MW2**3*Sf3641 + 10*Sf3639*MSbot2(2)**3)))*
     &    USbot21

	Sf3376 = MW2**3*Sf3375*MSbot2(1) + 
     &   (Sf3374*MSbot2(1)**5 + 
     &      MW2*(-(Sf3362*MSbot2(1)**4) + 
     &         MW2*(-4*MW2**3*Sf2553 + Sf3365*MSbot2(1)**3)))*
     &    USbot11

	Sf3356 = MW2**3*Sf3355*MSbot2(1) - 
     &   (Sf3354*MSbot2(1)**5 + 
     &      MW2*(-35*Sf3342*MSbot2(1)**4 + 
     &         MW2*(MW2**3*Sf3345 + 10*Sf3339*MSbot2(1)**3)))*
     &    USbot11

	Sf3969 = MW2**2*Sf3968*MSbot2(2) + 
     &   (3*Sf3967*MSbot2(2)**4 + 
     &      MW2*(MW2**3*Sf3636 - 2*Sf3964*MSbot2(2)**3))*USbot21

	Sf3963 = -(MW2**2*Sf3962*MSbot2(2)) + 
     &   (Sf3961*MSbot2(2)**4 + 
     &      MW2*(4*MW2**3*Sf2731 - 2*Sf3956*MSbot2(2)**3))*
     &    USbot21

	Sf3872 = MW2**2*Sf3871*MSbot2(1) + 
     &   (3*Sf3870*MSbot2(1)**4 + 
     &      MW2*(MW2**3*Sf3335 - 2*Sf3862*MSbot2(1)**3))*USbot11

	Sf3859 = -(MW2**2*Sf3858*MSbot2(1)) + 
     &   (Sf3857*MSbot2(1)**4 + 
     &      MW2*(4*MW2**3*Sf2553 - 2*Sf3845*MSbot2(1)**3))*
     &    USbot11

	Sf8714 = 2*MW2**2*Sf4820**3*Sf4906**3*MSq2Diff(1,2)**2*
     &    MSq2Diff(2,2)**2*
     &    (2*MW2*Sf8706*MStop2(2) - 
     &      CW2**2*MZ2**3*Sf8532*Sf8662*SW2*MStop2(3)**2) - 
     &   3*CW2**2*MZ2**3*(Sf8713 - 
     &      63*MW2*Sf4820**3*Sf4906**3*MSq2Diff(1,2)**2*
     &       MSq2Diff(2,2)**2)*MStop2(2)*MStop2(3)**2*UStop21**2

	Sf8677 = 2*MW2**2*Sf2320**3*Sf2334**3*MSq2Diff(1,1)**2*
     &    MSq2Diff(2,1)**2*
     &    (2*MW2*Sf8669*MStop2(1) - 
     &      CW2**2*MZ2**3*Sf7990*Sf8662*SW2*MStop2(3)**2) + 
     &   3*CW2**2*MZ2**3*(Sf8676 + 
     &      63*MW2*Sf2320**3*Sf2334**3*MSq2Diff(1,1)**2*
     &       MSq2Diff(2,1)**2)*MStop2(1)*MStop2(3)**2*UStop22**2

	Sf5083 = MHiggs2(4)*
     &    (-Sf5079 + MHiggs2(4)*(-Sf5082 + Sf5076*MHiggs2(4))) + 
     &   Sf5074*MSq2Diff(2,2)**3

	Sf4957 = MHiggs2(4)*
     &    (-Sf4953 + MHiggs2(4)*(-Sf4956 + Sf4950*MHiggs2(4))) + 
     &   Sf4948*MSq2Diff(1,2)**3

	Sf2643 = MHiggs2(4)*
     &    (-Sf2638 + MHiggs2(4)*(-Sf2642 + Sf2631*MHiggs2(4))) + 
     &   Sf2627*MSq2Diff(2,1)**3

	Sf2429 = MHiggs2(4)*
     &    (-Sf2421 + MHiggs2(4)*(-Sf2428 + Sf2412*MHiggs2(4))) + 
     &   Sf2406*MSq2Diff(1,1)**3

	Sf6809 = MW2*(4*MW2**3 - Sf6808) + 
     &   MHiggs2(4)*(3*MW2**2*Sf6803 + Sf6805*MHiggs2(4))

	Sf6554 = MW2*(4*MW2**3 - Sf6553) + 
     &   MHiggs2(4)*(3*MW2**2*Sf6548 + Sf6550*MHiggs2(4))

	Sf4469 = MW2*(4*MW2**3 - Sf4468) + 
     &   MHiggs2(4)*(3*MW2**2*Sf4458 + Sf4462*MHiggs2(4))

	Sf4352 = MW2*(4*MW2**3 - Sf4351) + 
     &   MHiggs2(4)*(3*MW2**2*Sf4337 + Sf4341*MHiggs2(4))

	Sf2000 = ML2*Sf1835*UStau21 + Sf1999*UStau21 - 
     &   ML*Sf1835*Sf215*UStau22

	Sf6385 = 3*Sf6384*Ss25*TB*MHiggs2(4) + 
     &   Sf1923*TB2*MSneu2(1)*
     &    (Sf162*UStau11 + Sf1191*UStau21)

	Sf1942 = ML2*Sf1046*UStau11 + Sf1941*UStau11 - 
     &   ML*Sf1046*Sf215*UStau12

	Sf7229 = -3*MZ2**3*Sf7228*SW2*MStop2(2) + 
     &   2*MW2*SB*Ss25**2*(MW2*Sf7226 + 
     &      9*MZ2**2*Sf7212*MStop2(2)*UStop21**2)

	Sf7210 = -3*MZ2**3*Sf7209*SW2*MStop2(1) + 
     &   2*MW2*SB*Ss25**2*(MW2*Sf7207 + 
     &      9*MZ2**2*Sf7189*MStop2(1)*UStop22**2)

	Sf8483 = -(Sf4820**3*Sf4906**3*Sf8475*USbot21**2*
     &      UStop22**2) + 
     &   Sf2334*(Sf2334**2*Sf8482 - 
     &      Sf4820**3*Sf4906**3*Sf8474*USbot21**2*UStop22**2)

	Sf6078 = MW2**2*Sf6077 - 
     &   2*MSbot2(2)**2*(MW2*Sf6071 + 10*Sf6076*MSbot2(2))*
     &    USbot21

	Sf5673 = MW2**2*Sf5672 - 
     &   2*MSbot2(1)**2*(MW2*Sf5658 + 10*Sf5671*MSbot2(1))*
     &    USbot11

	Sf4808 = -3*MZ2**3*Sf4807*SW2*MSbot2(2)*USbot21**2 + 
     &   MW2*(MW2*Sf4806 + 9*MZ2**2*Sf4804*MSbot2(2)*USbot21**2)

	Sf4787 = -3*MZ2**3*Sf4786*SW2*MSbot2(1)*USbot22**2 + 
     &   MW2*(MW2*Sf4785 + 9*MZ2**2*Sf4783*MSbot2(1)*USbot22**2)

	Sf1786 = -(MZ2**3*Sf1785*SW2*MStau2(2)*UStau21**2) + 
     &   MW2*(MW2*Sf1784 + 3*MZ2**2*Sf1768*MStau2(2)*UStau21**2)

	Sf1766 = -(MZ2**3*Sf1765*SW2*MStau2(1)*UStau22**2) + 
     &   MW2*(MW2*Sf1764 + 3*MZ2**2*Sf1738*MStau2(1)*UStau22**2)

	Sf1170 = MW2**2*Sf1169 + 
     &   MHiggs2(4)*(-(MW2*Sf1162) + Sf1155*MHiggs2(4))

	Sf8740 = MW2**2*MZ2*Sf8739 + Sf8734*MStop2(2)*UStop21

	Sf5444 = Sf5431*TB2*UStop21 - 
     &   9*MTop*Sf2245*Sf5439*TB*USbot21*UStop22

	Sf5338 = Sf5325*TB2*UStop21 - 
     &   9*MTop*Sf2245*Sf5333*TB*USbot11*UStop22

	Sf6293 = MTop*Sf6290*TB*UStop21 - 
     &   3*MTop2*MW2**2*Sf2245*Sf6291*Sf6292*MSbot2(2)*USbot21*
     &    UStop22

	Sf6279 = MTop*Sf6276*TB*UStop21 - 
     &   3*MTop2*MW2**2*Sf2245*Sf6277*Sf6278*MSbot2(1)*USbot11*
     &    UStop22

	Sf6023 = MTop*Sf6022*TB*UStop21 + 
     &   3*MTop2*MW2**3*Sf2245*Sf96*MSbot2(2)*MStop2(2)*USbot21*
     &    UStop22

	Sf5603 = MTop*Sf5602*TB*UStop21 + 
     &   3*MTop2*MW2**3*Sf2245*Sf88*MSbot2(1)*MStop2(2)*USbot11*
     &    UStop22

	Sf4886 = (Sf4885*UStop21)/(Sf4867*MHiggs2(4)) - 
     &   (3*MTop*Sf2245*USbot21*UStop22)/TB

	Sf4852 = (Sf4851*UStop21)/(Sf4831*MHiggs2(4)) - 
     &   (3*MTop*Sf2245*USbot11*UStop22)/TB

	Sf6415 = Sf2263*Sf4831*Sf4867*Sf6409*USbot21*
     &    UStop11 + Sf2270*Sf6414*UStop21

	Sf5219 = Sf2263*Sf4831*Sf4867*Sf5213*USbot21*
     &    UStop11 + Sf2270*Sf5218*UStop21

	Sf5220 = Sf2206*Sf5219 + 
     &   Sf2198*Sf2270*Sf4831*Sf4867*Sf5212*USbot11*UStop11

	Sf8762 = -(MW2**2*MZ2*Sf8761) + Sf8758*MStop2(1)*UStop11

	Sf3054 = Sf3035*TB2*UStop11 - 
     &   9*MTop*Sf2245*Sf3047*TB*USbot21*UStop12

	Sf2927 = Sf2908*TB2*UStop11 - 
     &   9*MTop*Sf2245*Sf2920*TB*USbot11*UStop12

	Sf3991 = MTop*Sf3988*TB*UStop11 - 
     &   3*MTop2*MW2**2*Sf2245*Sf3989*Sf3990*MSbot2(2)*USbot21*
     &    UStop12

	Sf3894 = MTop*Sf3891*TB*UStop11 - 
     &   3*MTop2*MW2**2*Sf2245*Sf3892*Sf3893*MSbot2(1)*USbot11*
     &    UStop12

	Sf3677 = MTop*Sf3676*TB*UStop11 + 
     &   3*MTop2*MW2**3*Sf2245*Sf74*MSbot2(2)*MStop2(1)*USbot21*
     &    UStop12

	Sf3397 = MTop*Sf3396*TB*UStop11 + 
     &   3*MTop2*MW2**3*Sf2245*Sf66*MSbot2(1)*MStop2(1)*USbot11*
     &    UStop12

	Sf2297 = (Sf2296*UStop11)/(Sf2270*MHiggs2(4)) - 
     &   (3*MTop*Sf2245*USbot21*UStop12)/TB

	Sf2246 = (Sf2244*UStop11)/(Sf2206*MHiggs2(4)) - 
     &   (3*MTop*Sf2245*USbot11*UStop12)/TB

	Sf1982 = Sf1973*UStau21**2 + 3*Sf1981*Sf215*UStau22

	Sf7769 = -(MW2**2*MZ2*Sf7768) + Sf7763*MStau2(2)*UStau21

	Sf2139 = Sf2131*UStau21 + ML*Sf2138*Sf215*UStau22

	Sf1817 = ML*Sf1816*UStau21 - 
     &   3*ML2*MW2**2*Sf1458*Sf215*Sf44*MSneu2(1)*UStau22

	Sf1637 = Sf1623*UStau21 + ML*Sf1632*Sf215*UStau22

	Sf1599 = Sf1591*UStau21 + ML*Sf1598*Sf215*UStau22

	Sf1570 = Sf1568*UStau21 - 
     &   ML*MW2*Sf1267*Sf1268*Sf1269*Sf1270*Sf1569*Sf215*
     &    MHiggs2(4)*MSl2Diff(2,1)*UStau22

	Sf1515 = Sf1507*UStau21 - ML*MW2**2*Sf1514*Sf215*UStau22

	Sf1455 = Sf1443*UStau21 + ML*MW2*Sf1454*Sf215*UStau22

	Sf1218 = Sf1210*UStau21 + ML*Sf1217*Sf215*UStau22

	Sf1922 = 3*Sf1921*Sf215*UStau12 + Sf1913*UStau22**2

	Sf7788 = MW2**2*MZ2*Sf7787 + Sf7782*MStau2(1)*UStau11

	Sf812 = Sf797*UStau11 + ML*Sf215*Sf807*UStau12

	Sf739 = Sf730*UStau11 + ML*Sf215*Sf738*UStau12

	Sf702 = Sf700*UStau11 - 
     &   ML*MW2*Sf215*Sf307*Sf308*Sf309*Sf310*Sf701*MHiggs2(4)*
     &    MSl2Diff(1,1)*UStau12

	Sf620 = Sf611*UStau11 - ML*MW2**2*Sf215*Sf619*UStau12

	Sf552 = Sf538*UStau11 + ML*MW2*Sf215*Sf551*UStau12

	Sf226 = Sf214*UStau11 + ML*Sf215*Sf225*UStau12

	Sf2059 = Sf2051*UStau11 + ML*Sf2058*Sf215*UStau12

	Sf1006 = ML*Sf1005*UStau11 - 
     &   3*ML2*MW2**2*Sf215*Sf31*Sf555*MSneu2(1)*UStau12

	Sf8519 = -(MW2**2*MZ2*Sf8518) + Sf8517*MSbot2(2)*USbot21

	Sf5405 = Sf5403*USbot21 + 9*MBot*Sf2225*Sf5404*USbot22

	Sf2997 = Sf2995*USbot21 + 9*MBot*Sf2225*Sf2996*USbot22

	Sf8503 = MW2**2*MZ2*Sf8502 + Sf8501*MSbot2(1)*USbot11

	Sf5299 = Sf5297*USbot11 + 9*MBot*Sf2225*Sf5298*USbot12

	Sf2859 = Sf2857*USbot11 + 9*MBot*Sf2225*Sf2858*USbot12

	Sf6285 = Sf3969*MStop2(2) - Sf6284*USbot21

	Sf6271 = Sf3872*MStop2(2) - Sf6270*USbot11

	Sf3974 = Sf3969*MStop2(1) - Sf3973*USbot21

	Sf3877 = Sf3872*MStop2(1) - Sf3876*USbot11

	Sf1787 = 2*CB*Sf1786*Ss25**2 + 3*MZ2**3*Sf1783*SW2*MStau2(2)

	Sf9592 = ((8*Sf9418*Sf9420)/(Ss49*MStau2(2)) + 
     &     (Sf9317**2*Sf9591)/MStau2(3)**3)/Sf9317**2

	Sf9767 = ((MZ2*Sf9739)/MStau2(2) - 
     &     (4*Sf9766*MStau2(3))/
     &      (Sf9344**2*Sf9362**2*Ss49*MHiggs2(2)*MStau2(1)))/MZ2

	Sf9738 = ((MZ2*Sf9708)/MStau2(2) + 
     &     (4*Sf9737*MStau2(3))/
     &      (Sf9344**2*Sf9362**2*Ss49*MHiggs2(2)*MStau2(1)))/MZ2

	Sf1887 = (MW2*Sf1191*Sf1886)/(Sf1184*MStau2(2)) + 
     &   8*Sf1191*Ss72*UStau21

	Sf9954 = 8*MZ2*Sf9290*Sf9291*Sf9930*MStau2(1)**3 - 
     &   Sf9953*MStau2(2)

	Sf9476 = 8*MZ2*Sf9290*Sf9291*Sf9332*Sf9333*Sf9432*
     &    MHiggs2(2)*MStau2(1)**2 + Sf9475*MStau2(2)

	Sf2182 = Sf1184*Sf1195**2*Sf162*Sf702*MSl2Diff(2,1)*
     &    MStau2(1) + Sf1191*Sf154*Sf1570*Sf166**2*MSl2Diff(1,1)*
     &    MStau2(2)

	Sf1767 = 2*CB*Sf1766*Ss25**2 + 3*MZ2**3*Sf1763*SW2*MStau2(1)

	Sf9796 = ((8*Sf9709*Sf9711)/(Ss49*MStau2(1)) - 
     &     (Sf9344**2*Sf9795)/MStau2(3)**3)/Sf9344**2

	Sf9533 = ((MZ2*Sf9479)/MStau2(1) + 
     &     (4*Sf9532*MStau2(3))/
     &      (Sf9317**2*Sf9331**2*Ss49*MHiggs2(2)*MStau2(2)))/MZ2

	Sf6257 = MW2*Sf4864*Sf6250*MSq2Diff(2,2) - 
     &   Sf4911*Sf6256*Ss25**3*USbot21*UStop21

	Sf6230 = 6*MW2*Sf5052*Sf5053*Sf5054*Sf5055*Sf5056*Sf6229*
     &    MHiggs2(4)*MSq2Diff(2,2) + 
     &   Sf6228*TB*USbot21*UStop21

	Sf5445 = 9*MW2**2*Sf5441*Sf5443*MSq2Diff(2,2) + 
     &   Sf5444*USbot21*UStop21

	Sf4171 = 6*MW2*Sf2598*Sf2599*Sf2600*Sf2601*Sf2602*Sf4170*
     &    MHiggs2(4)*MSq2Diff(2,1) + 
     &   Sf4169*TB*USbot21*UStop11

	Sf3055 = 9*MW2**2*Sf3050*Sf3053*MSq2Diff(2,1) + 
     &   Sf3054*USbot21*UStop11

	Sf5893 = MW2*Sf4828*Sf5886*MSq2Diff(1,2) - 
     &   Sf4905*Sf5892*Ss25**3*USbot11*UStop21

	Sf5864 = 6*MW2*Sf4926*Sf4927*Sf4928*Sf4929*Sf4930*Sf5863*
     &    MHiggs2(4)*MSq2Diff(1,2) + 
     &   Sf5862*TB*USbot11*UStop21

	Sf5339 = 9*MW2**2*Sf5335*Sf5337*MSq2Diff(1,2) + 
     &   Sf5338*USbot11*UStop21

	Sf4086 = 6*MW2*Sf2360*Sf2361*Sf2362*Sf2363*Sf2364*Sf4085*
     &    MHiggs2(4)*MSq2Diff(1,1) + 
     &   Sf4084*TB*USbot11*UStop11

	Sf2928 = 9*MW2**2*Sf2923*Sf2926*MSq2Diff(1,1) + 
     &   Sf2927*USbot11*UStop11

	Sf1891 = Sf1184*Sf1887*Ss25*MSneu2(1) + 
     &   6*MW2*Sf1191*Sf1831*MSl2Diff(2,1)*UStau21

	Sf7158 = Sf2139*MSneu2(1) + Sf1599*MStau2(2)

	Sf7160 = Sf2059*MSneu2(1) + Sf739*MStau2(1)

	Sf1275 = ML*Sf1267*Sf1268*Sf1269*Sf1270*Sf1273*
     &    MSl2Diff(2,1) - Sf1274*UStau21

	Sf315 = ML*Sf307*Sf308*Sf309*Sf310*Sf313*MSl2Diff(1,1) - 
     &   Sf314*UStau11

	Sf6920 = 5*Sf4728*MSbot2(2)**2 + Sf6919*MStop2(2)

	Sf4732 = 5*Sf4728*MSbot2(2)**2 + Sf4731*MStop2(1)

	Sf4809 = 2*CB*Sf4808*Ss25**2 + 3*MZ2**3*Sf4805*SW2*MSbot2(2)

	Sf6359 = Sf6230*Ss25**3*MSbot2(2) - 
     &   2*Sf4906**2*Sf6358*Ss241**2*Ss242**2*TB*USbot21*
     &    UStop21

	Sf4214 = Sf4171*Ss25**3*MSbot2(2) - 
     &   2*Sf2334**2*Sf4213*Ss241**2*Ss242**2*TB*USbot21*
     &    UStop11

	Sf6669 = 5*Sf4633*MSbot2(1)**2 + Sf6668*MStop2(2)

	Sf4637 = 5*Sf4633*MSbot2(1)**2 + Sf4636*MStop2(1)

	Sf4788 = 2*CB*Sf4787*Ss25**2 + 3*MZ2**3*Sf4784*SW2*MSbot2(1)

	Sf4211 = Sf4086*Ss25**3*MSbot2(1) - 
     &   2*Sf2320**2*Sf4210*Ss241**2*Ss242**2*TB*USbot11*
     &    UStop11

	Sf9976 = (MZ2*Sf9975*Ss49 - 
     &     Sf9954/(Sf9317*MStau2(1)*MStau2(2)))/MZ2

	Sf7116 = Sf6230*Ss25**2 - 
     &   4*Sf4906**2*Sf6144*Ss241**3*Ss242**3*TB*USbot21*
     &    UStop21

	Sf7113 = Sf5864*Ss25**2 - 
     &   4*Sf4820**2*Sf5761*Ss241**3*Ss242**3*TB*USbot11*
     &    UStop21

	Sf4742 = Sf4171*Ss25**2 - 
     &   4*Sf2334**2*Sf3742*Ss241**3*Ss242**3*TB*USbot21*
     &    UStop11

	Sf4739 = Sf4086*Ss25**2 - 
     &   4*Sf2320**2*Sf3733*Ss241**3*Ss242**3*TB*USbot11*
     &    UStop11

	Sf2014 = Sf1872*Ss25 + 2*Sf2005*Ss241*Ss242*MSl2Diff(2,1)

	Sf2105 = Sf1121*Sf166*Ss241*Ss242 - Sf1942*Ss25**2

	Sf9902 = 2*MZ2*Sf9317**2*Sf9900 - 2*Sf9331**2*Sf9901*MHiggs2(2)

	Sf9903 = Sf8929*Sf8930*Sf9871 - 
     &   (Sf55*Sf57*Sf9902)/
     &    (Sf9317**2*Sf9331**2*MStau2(1)*MStau2(2))

	Sf9838 = (Sf9837 + 
     &     (MZ2*Sf9824)/(Sf9362**2*MHiggs2(2)*MStau2(1)))/MZ2

	Sf9817 = (Sf9816 - 
     &     (MZ2*Sf9803)/(Sf9362**2*MHiggs2(2)*MStau2(1)))/MZ2

	Sf9685 = (Sf9684 + 
     &     (MZ2*Sf9668)/(Sf9331**2*MHiggs2(2)*MStau2(2)))/MZ2

	Sf9655 = (Sf9654 - 
     &     (MZ2*Sf9638)/(Sf9331**2*MHiggs2(2)*MStau2(2)))/MZ2

	Sf935 = 9*ML*Sf341*Sf927 - 22*MSl2Diff(1,1)**8*UStau11

	Sf813 = ML2*Sf810*Sf811 - Sf812*UStau11

	Sf7292 = Sf7276 + 2*MZ2*Sf7291*SW2

	Sf7285 = Sf7276 - 72*Sf7284

	Sf7238 = Sf7229 + 3*MZ2**2*SB*Sf7237*Ss25**2*UStop21**2

	Sf7233 = Sf7210 + 3*MZ2**2*SB*Sf7232*Ss25**2*UStop22**2

	Sf7034 = 3*Sf5056*Sf7033 - Sf7030*TB*USbot21*UStop21

	Sf7012 = 3*Sf4930*Sf7011 - Sf7008*TB*USbot11*UStop21

	Sf6847 = 3*Sf5056*Sf6846 + Sf6843*TB*USbot21*UStop21

	Sf6733 = 3*Sf2602*Sf6732 + Sf6725*TB*USbot21*UStop11

	Sf6592 = 3*Sf4930*Sf6591 + Sf6588*TB*USbot11*UStop21

	Sf6466 = 3*Sf2364*Sf6465 + Sf6458*TB*USbot11*UStop11

	Sf6416 = Sf2206*Sf6415 + 
     &   Sf2198*Sf2270*Sf4831*Sf4867*Sf6408*USbot11*UStop11

	Sf6348 = 3*MW2*Sf5056*Sf6347 + Sf6343*TB*USbot21*UStop21

	Sf6325 = 3*MW2*Sf4930*Sf6324 + Sf6320*TB*USbot11*UStop21

	Sf7148 = MW2*Sf4864*Sf6250 - 
     &   Sf4911*Sf6192*Ss25**3*USbot21*UStop21

	Sf6052 = 6*MW2**2*Sf5056*Sf6051 + 
     &   Sf6048*TB*USbot21*UStop21

	Sf6010 = 3*MW2**2*Sf5056*Sf6009 + 
     &   Sf6006*TB*USbot21*UStop21

	Sf5955 = 3*MW2**2*Sf5056*Sf5954 - 
     &   Sf5950*TB*USbot21*UStop21

	Sf7143 = MW2*Sf4828*Sf5886 - 
     &   Sf4905*Sf5814*Ss25**3*USbot11*UStop21

	Sf9977 = ((MZ2*Sf9976)/Ss49**2 + 
     &     (8*Sf55**2*Sf57**2)/MHiggs2(2))/MZ2

	Sf5633 = 6*MW2**2*Sf4930*Sf5632 + 
     &   Sf5629*TB*USbot11*UStop21

	Sf5590 = 3*MW2**2*Sf4930*Sf5589 + 
     &   Sf5586*TB*USbot11*UStop21

	Sf5529 = 3*MW2**2*Sf4930*Sf5528 - 
     &   Sf5524*TB*USbot11*UStop21

	Sf5221 = 2*Sf154*Sf5220 - 
     &   Sf162*Sf2206*Sf2270*Sf4831*Sf4867*Sf5211*TB*UStau11

	Sf5138 = 3*Sf5056*Sf5137 + Sf5134*TB*USbot21*UStop21

	Sf5116 = Sf5056*Sf5115 - Sf5112*TB*USbot21*UStop21

	Sf5012 = 3*Sf4930*Sf5011 + Sf5008*TB*USbot11*UStop21

	Sf4990 = Sf4930*Sf4989 - Sf4986*TB*USbot11*UStop21

	Sf4818 = Sf4809 + 3*CB*MZ2**2*Sf4817*Ss25**2*USbot21**2

	Sf4813 = Sf4788 + 3*CB*MZ2**2*Sf4812*Ss25**2*USbot22**2

	Sf4527 = 3*Sf2602*Sf4526 - Sf4517*TB*USbot21*UStop11

	Sf452 = Sf450 - Sf451*MSneu2(1)

	Sf4423 = 3*Sf2364*Sf4422 - Sf4413*TB*USbot11*UStop11

	Sf453 = Sf155**2*Sf441 + Sf452*MSneu2(1)

	Sf437 = 9*Sf347*Sf434 + Sf436*UStau11

	Sf4207 = 3*MW2*Sf2602*Sf4206 + Sf4199*TB*USbot21*UStop11

	Sf4127 = 3*MW2*Sf2364*Sf4126 + Sf4119*TB*USbot11*UStop11

	Sf6286 = (4*MBot2*MW2**6 - 4*MTop2*MW2**6 - 4*C2B*MW2**7 - 
     &     MW2**4*(-MBot2 + MTop2 + C2B*MW2)*
     &      (-3*MBot2 + 3*MTop2 - 5*MW2 + 3*C2B*MW2)*MSq2Diff(2,2) + 
     &     3*MW2**3*(2*(MBot2 - MTop2)**2 - 
     &        (-3 + 4*C2B)*(MBot2 - MTop2)*MW2 + 
     &        (-1 + C2B*(-3 + 2*C2B))*MW2**2)*MSq2Diff(2,2)**2 + 
     &     MW2**2*(-3*(MBot2 - MTop2)**2 + 
     &        (-43 + 6*C2B)*(MBot2 - MTop2)*MW2 - 
     &        (-16 + C2B*(-43 + 3*C2B))*MW2**2)*MSq2Diff(2,2)**3 - 
     &     MW2**2*(-71*MBot2 + 71*MTop2 + 34*MW2 + 71*C2B*MW2)*
     &      MSq2Diff(2,2)**4 + 
     &     12*MW2*(-4*MBot2 + 4*MTop2 + 3*MW2 + 4*C2B*MW2)*
     &      MSq2Diff(2,2)**5 + 
     &     (12*MBot2 - 12*MTop2 - 19*MW2 - 12*C2B*MW2)*
     &      MSq2Diff(2,2)**6 + 4*MSq2Diff(2,2)**7 - 
     &     MW2**4*(3*(MBot2 - MTop2)**2 - 
     &        (-37 + 6*C2B)*(MBot2 - MTop2)*MW2 + 
     &        (3 + C2B*(-37 + 3*C2B))*MW2**2)*MStop2(2) + 
     &     MW2**3*(-MBot2 + MTop2 + C2B*MW2)*
     &      (-12*MBot2 + 12*MTop2 - 25*MW2 + 12*C2B*MW2)*MSq2Diff(2,2)*
     &      MStop2(2) + MW2**2*
     &      (-9*(MBot2 - MTop2)**2 + (-107 + 18*C2B)*(MBot2 - MTop2)*MW2 - 
     &        (-48 + C2B*(-107 + 9*C2B))*MW2**2)*MSq2Diff(2,2)**2*
     &      MStop2(2) - 21*MW2**2*
     &      (-11*MBot2 + 11*MTop2 + 6*MW2 + 11*C2B*MW2)*
     &      MSq2Diff(2,2)**3*MStop2(2) + 
     &     MW2*(-112*MBot2 + 112*MTop2 + 123*MW2 + 112*C2B*MW2)*
     &      MSq2Diff(2,2)**4*MStop2(2) - 
     &     42*MW2*MSq2Diff(2,2)**5*MStop2(2) + 
     &     2*MW2**3*(3*(MBot2 - MTop2)**2 - 
     &        2*(-29 + 3*C2B)*(MBot2 - MTop2)*MW2 + 
     &        (11 + C2B*(-58 + 3*C2B))*MW2**2)*MStop2(2)**2 - 
     &     2*MW2**2*(3*(MBot2 - MTop2)**2 - 
     &        2*(-8 + 3*C2B)*(MBot2 - MTop2)*MW2 + 
     &        (-21 + C2B*(-16 + 3*C2B))*MW2**2)*MSq2Diff(2,2)*
     &      MStop2(2)**2 - 
     &     6*MW2**2*(-50*MBot2 + 50*MTop2 + 29*MW2 + 50*C2B*MW2)*
     &      MSq2Diff(2,2)**2*MStop2(2)**2 + 
     &     142*MW2**2*MSq2Diff(2,2)**3*MStop2(2)**2 + 
     &     8*MW2**3*(-16*MBot2 + 16*MTop2 - 5*MW2 + 16*C2B*MW2)*
     &      MStop2(2)**3 - 152*MW2**3*MSq2Diff(2,2)*MStop2(2)**3)*
     &   USbot21**2

	Sf6286 = Sf6286 + 
     &   MBot*(reMUE*TB - Re(Abot))*
     &    (-((MW2 - MSq2Diff(2,2))**2*
     &         (-6*MBot2*MW2**2*MSq2Diff(2,2) + 
     &           6*MTop2*MW2**2*MSq2Diff(2,2) + 
     &           6*C2B*MW2**3*MSq2Diff(2,2) + 
     &           12*MSq2Diff(2,2)**4 + 
     &           MW2*(-24*MSq2Diff(2,2)**3 + 
     &              MW2*(11*MSq2Diff(2,2)**2 + 
     &                 MW2*(4*MW2 + 3*MSq2Diff(2,2)))))) + 
     &      MW2*(MW2 - MSq2Diff(2,2))*
     &       (MW2**2*(6*MBot2 - 6*MTop2 + 37*MW2 - 6*C2B*MW2) + 
     &         6*MW2*(-3*MBot2 + 3*MTop2 + 2*MW2 + 3*C2B*MW2)*
     &          MSq2Diff(2,2) + 119*MW2*MSq2Diff(2,2)**2 - 
     &         112*MSq2Diff(2,2)**3)*MStop2(2) - 
     &      4*MW2**2*(MW2*(3*MBot2 - 3*MTop2 + 29*MW2 - 3*C2B*MW2) + 
     &         (-3*MBot2 + 3*MTop2 - 8*MW2 + 3*C2B*MW2)*
     &          MSq2Diff(2,2) + 75*MSq2Diff(2,2)**2)*MStop2(2)**2+
     &        128*MW2**3*MStop2(2)**3)*USbot21*USbot22 - 
     &   3*MBot2*MW2**2*(-(reMUE*TB) + Re(Abot))**2*MSbot2(2)*
     &    ((MW2 - MSbot2(2))**2 - MStop2(2)**2)*USbot22**2

	Sf3975 = Sf2722**2*Sf3963 + Sf3974*MStop2(1)

	Sf394 = 9*ML*Sf341*Sf387 - Sf393*UStau11

	Sf6272 = Sf2543**2*Sf3859 + Sf6271*MStop2(2)

	Sf3878 = Sf2543**2*Sf3859 + Sf3877*MStop2(1)

	Sf3720 = 6*MW2**2*Sf2602*Sf3719 + 
     &   Sf3712*TB*USbot21*UStop11

	Sf6015 = Sf3638 - Sf6014*MStop2(2)

	Sf3659 = Sf3638 - Sf3658*MStop2(1)

	Sf3612 = 3*MW2**2*Sf2602*Sf3611 + 
     &   Sf3603*TB*USbot21*UStop11

	Sf3530 = 3*MW2**2*Sf2602*Sf3529 - 
     &   Sf3518*TB*USbot21*UStop11

	Sf3448 = 6*MW2**2*Sf2364*Sf3447 + 
     &   Sf3440*TB*USbot11*UStop11

	Sf349 = 9*ML*Sf340*Sf341 - 76*MSl2Diff(1,1)**7*UStau11

	Sf5595 = Sf3337 - Sf5594*MStop2(2)

	Sf3379 = Sf3337 - Sf3378*MStop2(1)

	Sf3285 = 3*MW2**2*Sf2364*Sf3284 + 
     &   Sf3276*TB*USbot11*UStop11

	Sf3168 = 3*MW2**2*Sf2364*Sf3167 - 
     &   Sf3156*TB*USbot11*UStop11

	Sf2729 = 3*Sf2602*Sf2728 + Sf2721*TB*USbot21*UStop11

	Sf2690 = Sf2602*Sf2689 - Sf2683*TB*USbot21*UStop11

	Sf2550 = 3*Sf2364*Sf2549 + Sf2542*TB*USbot11*UStop11

	Sf2498 = Sf2364*Sf2497 - Sf2490*TB*USbot11*UStop11

	Sf2106 = Sf162*Sf2105*MHiggs2(4) + 
     &   Sf1146*Sf166**2*Ss241*Ss242*MStau2(1)

	Sf2013 = Sf2000 + 2*Sf1994*MSl2Diff(2,1)

	Sf2020 = Sf1787 + CB*MZ2**2*Sf50*Ss25**2*UStau21**2

	Sf1789 = Sf1767 + CB*MZ2**2*Sf37*Ss25**2*UStau22**2

	Sf1724 = 9*ML*Sf1293*Sf1717 + 22*MSl2Diff(2,1)**8*UStau21

	Sf1638 = ML2*Sf1635*Sf1636 + Sf1637*UStau21

	Sf1369 = Sf1367 - Sf1368*MSneu2(1)

	Sf1370 = Sf1185**2*Sf1361 + Sf1369*MSneu2(1)

	Sf1360 = 9*Sf1298*Sf1357 - Sf1359*UStau21

	Sf1331 = 9*ML*Sf1293*Sf1324 + Sf1330*UStau21

	Sf1300 = 9*ML*Sf1292*Sf1293 - 76*MSl2Diff(2,1)**7*UStau21

	Sf1176 = 3*Sf1170 - 2*Sf1175*Ss25**2

	Sf1726 = 18*MW2**3*Sf1725 + Sf1724*UStau21 + 
     &   MW2*(C2B*MW2**4*(-18*MW2**2*Sf1699 + 27*Sf1703) - 
     &      Sf1723 + MW2*(MW2**2*
     &          (Sf1722 + MW2**2*(Sf1702 + 2*MW2**2*Sf931)) + 
     &         Sf1721*MSl2Diff(2,1)**4))*UStau21**2

	Sf937 = 18*MW2**3*Sf936 - Sf935*UStau11 + 
     &   MW2*(C2B*MW2**4*(-18*MW2**2*Sf900 + 27*Sf904) - Sf934 + 
     &      MW2*(MW2**2*(MW2**2*(Sf903 + 2*MW2**2*Sf931) + 
     &            Sf933) + Sf932*MSl2Diff(1,1)**4))*UStau22**2

	Sf2162 = -3*ML2*Sf2160 + 10*MSl2Diff(2,1)**7 + 
     &   MW2*(-Sf2161 + MW2*
     &       (3*Sf2143*MSl2Diff(2,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 - Sf2140) + 3*Sf2142) - 
     &               Sf2145) - 2*Sf2144*MSl2Diff(2,1))))

	Sf2093 = -3*ML2*Sf2091 + 10*MSl2Diff(1,1)**7 + 
     &   MW2*(-Sf2092 + MW2*
     &       (3*Sf2063*MSl2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 - Sf2060) + 3*Sf2062) - 
     &               Sf2065) - 2*Sf2064*MSl2Diff(1,1))))

	Sf6880 = Sf6869 - 6*MTop2*Sf6875

	Sf6766 = Sf6755 - 6*MTop2*Sf6761

	Sf6629 = Sf6617 - 6*MTop2*Sf6624

	Sf6511 = Sf6500 - 6*MTop2*Sf6506

	Sf7106 = Sf7096 + 2*MTop2*Sf7100

	Sf7101 = Sf7096 - 2*MBot2*Sf7100

	Sf7068 = Sf7056 + 2*MTop2*Sf7060

	Sf7061 = Sf7056 - 2*MBot2*Sf7060

	Sf4675 = Sf4665 + 2*MTop2*Sf4669

	Sf4670 = Sf4665 - 2*MBot2*Sf4669

	Sf4570 = Sf4560 + 2*MTop2*Sf4564

	Sf4565 = Sf4560 - 2*MBot2*Sf4564

	Sf5369 = Sf5361 - 2*MTop2*Sf5366

	Sf5367 = Sf5361 + 2*MBot2*Sf5366

	Sf5261 = Sf5253 - 2*MTop2*Sf5258

	Sf5259 = Sf5253 + 2*MBot2*Sf5258

	Sf2958 = Sf2950 - 2*MTop2*Sf2955

	Sf2956 = Sf2950 + 2*MBot2*Sf2955

	Sf2805 = Sf2797 - 2*MTop2*Sf2802

	Sf2803 = Sf2797 + 2*MBot2*Sf2802

	Sf2167 = Sf2160 + 2*ML2*Sf2166

	Sf2098 = Sf2091 + 2*ML2*Sf2097

	Sf1302 = 3*MW2**2*Sf1301 - Sf1300*UStau21 + 
     &   MW2*(4*MW2**2*Sf1279 + 3*Sf1299 + 
     &      MW2**3*(49*MW2**2 + 6*Sf1277)*MSl2Diff(2,1) - 
     &      348*MW2**4*Sf47*MSl2Diff(2,1))*UStau21**2

	Sf6820 = -(MW2**3*Sf6817) + Sf6816*MSq2Diff(2,2)**4 + 
     &   MHiggs2(4)*(Sf6819 + 
     &      MW2*MHiggs2(4)*
     &       (Sf6815*MHiggs2(4) - Sf6809*MSq2Diff(2,2)))

	Sf6565 = -(MW2**3*Sf6562) + Sf6561*MSq2Diff(1,2)**4 + 
     &   MHiggs2(4)*(Sf6564 + 
     &      MW2*MHiggs2(4)*
     &       (Sf6560*MHiggs2(4) - Sf6554*MSq2Diff(1,2)))

	Sf4481 = -(MW2**3*Sf4478) + Sf4477*MSq2Diff(2,1)**4 + 
     &   MHiggs2(4)*(Sf4480 + 
     &      MW2*MHiggs2(4)*
     &       (Sf4476*MHiggs2(4) - Sf4469*MSq2Diff(2,1)))

	Sf4370 = -(MW2**3*Sf4367) + Sf4366*MSq2Diff(1,1)**4 + 
     &   MHiggs2(4)*(Sf4369 + 
     &      MW2*MHiggs2(4)*
     &       (Sf4365*MHiggs2(4) - Sf4352*MSq2Diff(1,1)))

	Sf1663 = Sf1184*(Sf1638*Ss25**3 + 
     &      9*MW2*Sf1191*Sf1195*Sf1599*TB)*MSl2Diff(2,1) - 
     &   Sf1184*Sf1662*Ss25**3*MSneu2(1) - 
     &   9*MW2**2*Sf1189*Sf1191*Sf1195**3*TB*MSl2Diff(2,1)*
     &    UStau21

	Sf846 = Sf154*(Sf813*Ss25**3 + 
     &      9*MW2*Sf162*Sf166*Sf739*TB)*MSl2Diff(1,1) - 
     &   Sf154*Sf845*Ss25**3*MSneu2(1) - 
     &   9*MW2**2*Sf159*Sf162*Sf166**3*TB*MSl2Diff(1,1)*UStau11

	Sf7173 = -3*MW2**3*Sf6363 + Sf7172*Ss25**3*TB + 
     &   MW2*MHiggs2(4)*(-3*MW2*Sf6363 + 3*Sf6363*MHiggs2(4))

	Sf2006 = Sf1184*Sf154*Sf1887*Ss25*MStau2(2) - 
     &   12*MW2*MSl2Diff(2,1)*
     &    (Sf1184*Sf162*Sf1926*UStau11 + 
     &      Sf1191*Sf154*Sf1985*UStau21)

	Sf1828 = Sf1820*MSneu2(1) + Sf1822*MSneu2(1)**2 - 
     &   Sf1827*MStau2(2)

	Sf1036 = Sf1013*MSneu2(1) + Sf1020*MSneu2(1)**2 - 
     &   Sf1035*MStau2(1)

	Sf1664 = -2*CW2*Sf1663 + 
     &   Sf1184*Sf1195**3*(3*Sf1182*Ss25**3 - 9*MW2**3*Sf1183*TB)*
     &    MSl2Diff(2,1)

	Sf847 = -2*CW2*Sf846 + 
     &   Sf154*Sf166**3*(3*Sf151*Ss25**3 - 9*MW2**3*Sf153*TB)*
     &    MSl2Diff(1,1)

	Sf7178 = Sf6466*MSbot2(1) + Sf4423*MStop2(1) - 
     &   Sf173*Sf2320**2*Ss25*TB*MSq2Diff(1,1)*USbot11*
     &    UStop11

	Sf6382 = Sf5759*MSbot2(1) - Sf4127*MStop2(1) + 
     &   Sf2320**2*Ss25*TB*MHiggs2(4)*MSq2Diff(1,1)*USbot11*
     &    UStop11

	Sf622 = (3*MW2*Sf162*Sf166*Sf620*TB)/Ss25**3 - 
     &   Sf621*UStau11 - 
     &   6*ML2*MW2**3*Sf215**2*Sf31*MSneu2(1)*MStau2(1)*
     &    UStau12**2

	Sf1517 = (3*MW2*Sf1191*Sf1195*Sf1515*TB)/Ss25**3 - 
     &   Sf1516*UStau21 - 
     &   6*ML2*MW2**3*Sf215**2*Sf44*MSneu2(1)*MStau2(2)*
     &    UStau22**2

	Sf7239 = -(MW2*MZ2**2*Sf7238) + 
     &   SB*Ss25**2*(MW2**3*Sf7235 + 3*MZ2**4*Sf7236*UStop21**2)

	Sf7234 = -(MW2*MZ2**2*Sf7233) + 
     &   SB*Ss25**2*(MW2**3*Sf7230 + 3*MZ2**4*Sf7231*UStop22**2)

	Sf4819 = -(MW2*MZ2**2*Sf4818) + 
     &   CB*Ss25**2*(MW2**3*Sf4815 - 3*MZ2**4*Sf4816*USbot21**2)

	Sf4814 = -(MW2*MZ2**2*Sf4813) + 
     &   CB*Ss25**2*(MW2**3*Sf4810 - 3*MZ2**4*Sf4811*USbot22**2)

	Sf2012 = Sf1121*Sf166*Ss241*Ss242*Ss25 - 
     &   Sf1942*Ss25**3 + 2*Sf2011*MSl2Diff(1,1)

	Sf9477 = Sf9317*Sf9331*Sf9476 + 
     &   MZ2*(Sf9317**2*Sf9429 - 
     &      16*Sf9331**2*Sf9418*Sf9420*MHiggs2(2)*MStau2(3)**3)

	Sf554 = MW2**2*(6*MW2*Sf31*Sf347*Sf498 - 
     &      (9*Sf162*Sf166*Sf552*TB)/Ss25**3) - Sf553*UStau11

	Sf1457 = MW2**2*(6*MW2*Sf1298*Sf1415*Sf44 - 
     &      (9*Sf1191*Sf1195*Sf1455*TB)/Ss25**3) - 
     &   Sf1456*UStau21

	Sf2021 = MW2*MZ2**2*Sf2020 + 
     &   CB*Ss25**2*(-3*MW2**3*Sf2019 + 
     &      MZ2**4*Sf50*MSl2Diff(2,1)*UStau21**2)

	Sf1790 = MW2*MZ2**2*Sf1789 + 
     &   CB*Ss25**2*(-3*MW2**3*Sf1788 + 
     &      MZ2**4*Sf37*MSl2Diff(1,1)*UStau22**2)

	Sf7103 = 3*MBot2*Sf7096 - 3*MTop2*Sf7101 + 
     &   10*MSq2Diff(2,2)**7 + 
     &   MW2*(-Sf7102 + MW2*
     &       (3*Sf7078*MSq2Diff(2,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 + Sf7075) - 3*Sf7077) + 
     &               Sf7079) - 2*Sf7080*MSq2Diff(2,2))))

	Sf7063 = 3*MBot2*Sf7056 - 3*MTop2*Sf7061 + 
     &   10*MSq2Diff(1,2)**7 + 
     &   MW2*(-Sf7062 + MW2*
     &       (3*Sf7038*MSq2Diff(1,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 + Sf7035) - 3*Sf7037) + 
     &               Sf7039) - 2*Sf7040*MSq2Diff(1,2))))

	Sf4672 = 3*MBot2*Sf4665 - 3*MTop2*Sf4670 + 
     &   10*MSq2Diff(2,1)**7 + 
     &   MW2*(-Sf4671 + MW2*
     &       (3*Sf4647*MSq2Diff(2,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 + Sf4644) - 3*Sf4646) + 
     &               Sf4648) - 2*Sf4649*MSq2Diff(2,1))))

	Sf4567 = 3*MBot2*Sf4560 - 3*MTop2*Sf4565 + 
     &   10*MSq2Diff(1,1)**7 + 
     &   MW2*(-Sf4566 + MW2*
     &       (3*Sf4531*MSq2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 + Sf4528) - 3*Sf4530) + 
     &               Sf4532) - 2*Sf4533*MSq2Diff(1,1))))

	Sf6878 = 6*MBot2*Sf6875 + Sf6880

	Sf6764 = 6*MBot2*Sf6761 + Sf6766

	Sf6627 = 6*MBot2*Sf6624 + Sf6629

	Sf6509 = 6*MBot2*Sf6506 + Sf6511

	Sf7104 = -2*MBot2*Sf7100 + Sf7106

	Sf7066 = -2*MBot2*Sf7065 + Sf7068

	Sf4673 = -2*MBot2*Sf4669 + Sf4675

	Sf4568 = -2*MBot2*Sf4564 + Sf4570

	Sf1860 = C2B*MW2*Sf1852 - 3*ML2*Sf1859 - 
     &   8*MSl2Diff(2,1)**7 + 
     &   MW2*(42*Sf47*MSl2Diff(2,1)**5 + 
     &      MW2*(-(Sf1844*MSl2Diff(2,1)**3) + 
     &         MW2*(-8*Sf1846 + 
     &            MW2*(6*Sf1845 + 
     &               MW2*(-3*MW2*MSl2Diff(2,1) + 
     &                  22*Sf47*MSl2Diff(2,1))))))

	Sf1087 = C2B*MW2*Sf1073 - 3*ML2*Sf1086 - 
     &   8*MSl2Diff(1,1)**7 + 
     &   MW2*(42*Sf34*MSl2Diff(1,1)**5 + 
     &      MW2*(-(Sf1065*MSl2Diff(1,1)**3) + 
     &         MW2*(-8*Sf1067 + 
     &            MW2*(6*Sf1066 + 
     &               MW2*(-3*MW2*MSl2Diff(1,1) + 
     &                  22*Sf34*MSl2Diff(1,1))))))

	Sf6166 = 4*C2B*MW2**7 - Sf6160*Ss30 - 
     &   8*MSq2Diff(2,2)**7 + 
     &   MW2*(2*Sf6147*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-(MW2*Sf6145) + Sf6149) - Sf6151) + Sf6153
     &            ) - Sf6165*MSq2Diff(2,2)))

	Sf5788 = 4*C2B*MW2**7 - Sf5781*Ss30 - 
     &   8*MSq2Diff(1,2)**7 + 
     &   MW2*(2*Sf5764*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-(MW2*Sf5762) + Sf5766) - Sf5768) + Sf5774
     &            ) - Sf5787*MSq2Diff(1,2)))

	Sf3918 = 4*C2B*MW2**7 - Sf3912*Ss30 - 
     &   8*MSq2Diff(2,1)**7 + 
     &   MW2*(2*Sf3899*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-(MW2*Sf3896) + Sf3901) - Sf3903) + Sf3905
     &            ) - Sf3917*MSq2Diff(2,1)))

	Sf3793 = 4*C2B*MW2**7 - Sf3784*Ss30 - 
     &   8*MSq2Diff(1,1)**7 + 
     &   MW2*(2*Sf3754*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*
     &             (MW2*(-(MW2*Sf3751) + Sf3758) - Sf3764) + Sf3770
     &            ) - Sf3792*MSq2Diff(1,1)))

	Sf5927 = -6*MBot2*Sf5921 + 6*MTop2*Sf5925 - 
     &   10*MSq2Diff(2,2)**8 + 
     &   MW2*(Sf5926 + MW2*
     &       (MW2**2*(MW2**2*(2*MW2**2 - Sf5903) - 2*Sf5905) - 
     &         Sf5904*MSq2Diff(2,2)**4))

	Sf5494 = -6*MBot2*Sf5488 + 6*MTop2*Sf5492 - 
     &   10*MSq2Diff(1,2)**8 + 
     &   MW2*(Sf5493 + MW2*
     &       (MW2**2*(MW2**2*(2*MW2**2 - Sf5469) - 2*Sf5471) - 
     &         Sf5470*MSq2Diff(1,2)**4))

	Sf3483 = -6*MBot2*Sf3477 + 6*MTop2*Sf3481 - 
     &   10*MSq2Diff(2,1)**8 + 
     &   MW2*(Sf3482 + MW2*
     &       (MW2**2*(MW2**2*(2*MW2**2 - Sf3459) - 2*Sf3461) - 
     &         Sf3460*MSq2Diff(2,1)**4))

	Sf3100 = -6*MBot2*Sf3094 + 6*MTop2*Sf3098 - 
     &   10*MSq2Diff(1,1)**8 + 
     &   MW2*(Sf3099 + MW2*
     &       (MW2**2*(MW2**2*(2*MW2**2 - Sf3072) - 2*Sf3074) - 
     &         Sf3073*MSq2Diff(1,1)**4))

	Sf5368 = -9*C2B*MW2*Sf5351 + 9*MBot2*Sf5361 - 
     &   9*MTop2*Sf5367 + 22*MSq2Diff(2,2)**7 + 
     &   MW2*(-(Sf5340*MSq2Diff(2,2)**5) + 
     &      MW2*(3*Sf5342*MSq2Diff(2,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf5341) + 3*Sf5343) - 
     &               2*Sf5345) - 2*Sf5344*MSq2Diff(2,2))))

	Sf5260 = -9*C2B*MW2*Sf5243 + 9*MBot2*Sf5253 - 
     &   9*MTop2*Sf5259 + 22*MSq2Diff(1,2)**7 + 
     &   MW2*(-(Sf5232*MSq2Diff(1,2)**5) + 
     &      MW2*(3*Sf5234*MSq2Diff(1,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf5233) + 3*Sf5235) - 
     &               2*Sf5237) - 2*Sf5236*MSq2Diff(1,2))))

	Sf2957 = -9*C2B*MW2*Sf2940 + 9*MBot2*Sf2950 - 
     &   9*MTop2*Sf2956 + 22*MSq2Diff(2,1)**7 + 
     &   MW2*(-(Sf2929*MSq2Diff(2,1)**5) + 
     &      MW2*(3*Sf2931*MSq2Diff(2,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf2930) + 3*Sf2932) - 
     &               2*Sf2934) - 2*Sf2933*MSq2Diff(2,1))))

	Sf2804 = -9*C2B*MW2*Sf2774 + 9*MBot2*Sf2797 - 
     &   9*MTop2*Sf2803 + 22*MSq2Diff(1,1)**7 + 
     &   MW2*(-(Sf2763*MSq2Diff(1,1)**5) + 
     &      MW2*(3*Sf2765*MSq2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(22*MW2 - Sf2764) + 3*Sf2766) - 
     &               2*Sf2768) - 2*Sf2767*MSq2Diff(1,1))))

	Sf6877 = 3*C2B*MW2*Sf6859 - MBot2*Sf6869 + 3*MTop2*Sf6876 - 
     &   10*MSq2Diff(2,2)**7 + 
     &   MW2*(Sf6849*MSq2Diff(2,2)**5 + 
     &      MW2*(-3*Sf6851*MSq2Diff(2,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 - Sf6848) + Sf6850) - 
     &               Sf6853) + 2*Sf6852*MSq2Diff(2,2))))

	Sf6763 = 3*C2B*MW2*Sf6745 - MBot2*Sf6755 + 3*MTop2*Sf6762 - 
     &   10*MSq2Diff(2,1)**7 + 
     &   MW2*(Sf6735*MSq2Diff(2,1)**5 + 
     &      MW2*(-3*Sf6737*MSq2Diff(2,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 - Sf6734) + Sf6736) - 
     &               Sf6739) + 2*Sf6738*MSq2Diff(2,1))))

	Sf6626 = 3*C2B*MW2*Sf6604 - MBot2*Sf6617 + 3*MTop2*Sf6625 - 
     &   10*MSq2Diff(1,2)**7 + 
     &   MW2*(Sf6594*MSq2Diff(1,2)**5 + 
     &      MW2*(-3*Sf6596*MSq2Diff(1,2)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 - Sf6593) + 3*Sf6595) - 
     &               Sf6598) + 2*Sf6597*MSq2Diff(1,2))))

	Sf6508 = 3*C2B*MW2*Sf6478 - MBot2*Sf6500 + 3*MTop2*Sf6507 - 
     &   10*MSq2Diff(1,1)**7 + 
     &   MW2*(Sf6468*MSq2Diff(1,1)**5 + 
     &      MW2*(-3*Sf6470*MSq2Diff(1,1)**3 + 
     &         MW2*(MW2*(MW2*(MW2*(2*MW2 - Sf6467) + Sf6469) - 
     &               Sf6472) + 2*Sf6471*MSq2Diff(1,1))))

	Sf5982 = -18*C2B*MW2*Sf5970 + 18*MBot2*Sf5977 - 
     &   18*MTop2*Sf5981 + 22*MSq2Diff(2,2)**8 + 
     &   MW2*(-(Sf5963*Sf99) + 
     &      MW2*(MW2**2*(MW2**2*(22*MW2**2 + Sf5959) + Sf5960) + 
     &         115*Sf5958*MSq2Diff(2,2)**4))

	Sf5556 = -18*C2B*MW2*Sf5544 + 18*MBot2*Sf5551 - 
     &   18*MTop2*Sf5555 + 22*MSq2Diff(1,2)**8 + 
     &   MW2*(-(Sf5537*Sf91) + 
     &      MW2*(MW2**2*(MW2**2*(22*MW2**2 + Sf5533) + Sf5534) + 
     &         115*Sf5532*MSq2Diff(1,2)**4))

	Sf3557 = -18*C2B*MW2*Sf3545 + 18*MBot2*Sf3552 - 
     &   18*MTop2*Sf3556 + 22*MSq2Diff(2,1)**8 + 
     &   MW2*(-(Sf3538*Sf77) + 
     &      MW2*(MW2**2*(MW2**2*(22*MW2**2 + Sf3534) + Sf3535) + 
     &         115*Sf3533*MSq2Diff(2,1)**4))

	tmp99 = (3 + C2B)*MSbot2(1) - 8*MStop2(1) + C2B*MStop2(1)

	tmp114 = 5*MStop2(1)**2 + MSbot2(1)*(6*MSbot2(1) + 9*MStop2(1))

	tmp115 = -122*MStop2(1)**2 + C2B*MStop2(1)**2 + 
     &   MSbot2(1)*((78 - 42*C2B)*MSbot2(1) - 102*MStop2(1) - 
     &      42*C2B*MStop2(1))

	tmp116 = -115*MStop2(1)**3 + 
     &   MSbot2(1)*(24*MStop2(1)**2 + 
     &      MSbot2(1)*(245*MSbot2(1) + 49*MStop2(1)))

	Sf3201 = 22*MSq2Diff(1,1)**8 + 
     &   MW2*(-((MSbot2(1) + MStop2(1))*
     &         (230*MW2**6 + 203*MSq2Diff(1,1)**6 + 
     &           MW2**4*(MSbot2(1)*
     &               (1907*MSbot2(1) - 80*MStop2(1)) + 
     &              1907*MStop2(1)**2) + 
     &           2*MW2**2*MSq2Diff(1,1)**2*
     &            (886*MStop2(1)**2 + 
     &              MSbot2(1)*(886*MSbot2(1) + 149*MStop2(1))))) + 
     &      MW2*(22*MW2**6 + 
     &         115*MSq2Diff(1,1)**4*
     &          (7*MStop2(1)**2 + 
     &            MSbot2(1)*(7*MSbot2(1) + 10*MStop2(1))) + 
     &         MW2**4*(913*MStop2(1)**2 + 
     &            MSbot2(1)*(913*MSbot2(1) + 1096*MStop2(1))) + 
     &         MW2**2*(2350*MStop2(1)**4 + 
     &            MSbot2(1)*
     &             (718*MStop2(1)**3 + 
     &               MSbot2(1)*
     &                (792*MStop2(1)**2 + 
     &                  MSbot2(1)*(2350*MSbot2(1) + 718*MStop2(1)))
     &               )))) - 
     &   18*MBot2*(-2*MSq2Diff(1,1)**7 + 
     &      MW2*(2*MSq2Diff(1,1)**5*
     &          (11*MSbot2(1) + 10*MStop2(1)) + 
     &         MW2*(-18*tmp114*MSq2Diff(1,1)**3 + 
     &            MW2*(tmp116*MSbot2(1) - 203*MStop2(1)**4 + 
     &               MW2*(-2*
     &                   ((136 + C2B)*MSbot2(1)**3 + 
     &                     tmp115*MStop2(1)) + 
     &                  MW2*
     &                   (MW2*(-5*MW2 - 6*tmp99) + 
     &                     MSbot2(1)*
     &                      (138*MSbot2(1) - 27*MStop2(1)) - 
     &                     156*MStop2(1)**2 + 
     &                     C2B*
     &                      (2*MW2**2 + 
     &                       MSbot2(1)*
     &                       (6*MSbot2(1) - 78*MStop2(1)) + 
     &                       6*MStop2(1)**2))))))) - 
     &   18*C2B*MW2*(2*MSq2Diff(1,1)**7 + 
     &      MW2*(-2*MSq2Diff(1,1)**5*
     &          (11*MSbot2(1) + 10*MStop2(1)) + 
     &         MW2*(18*tmp114*MSq2Diff(1,1)**3 + 
     &            MW2*(-(tmp116*MSbot2(1)) + 203*MStop2(1)**4 + 
     &               MW2*(MW2*
     &                   (5*MW2**2 + 
     &                     6*MW2*(3*MSbot2(1) - 8*MStop2(1)) + 
     &                     3*
     &                      (52*MStop2(1)**2 + 
     &                       MSbot2(1)*
     &                       (-46*MSbot2(1) + 9*MStop2(1)))) + 
     &                  4*(-61*MStop2(1)**3 + 
     &                     MSbot2(1)*
     &                      (-51*MStop2(1)**2 + 
     &                       MSbot2(1)*
     &                       (68*MSbot2(1) + 39*MStop2(1)))))))))-
     &     18*MTop2*(2*MSq2Diff(1,1)**7 + 
     &      MW2*(-2*MSq2Diff(1,1)**5*
     &          (11*MSbot2(1) + 10*MStop2(1)) + 
     &         MW2*(18*tmp114*MSq2Diff(1,1)**3 + 
     &            MW2*(-(tmp116*MSbot2(1)) + 203*MStop2(1)**4 + 
     &               2*MBot2*(MW2 - MSbot2(1) - MStop2(1))*
     &                ((MW2 - MSbot2(1))**2 - 
     &                  (2*MW2 + 43*MSbot2(1))*MStop2(1) + 
     &                  MStop2(1)**2) + 
     &               MW2*(2*
     &                   ((136 + C2B)*MSbot2(1)**3 + 
     &                     tmp115*MStop2(1)) + 
     &                  MW2*
     &                   (MW2*(5*MW2 + 6*tmp99) + 
     &                     156*MStop2(1)**2 + 
     &                     MSbot2(1)*
     &                      (-138*MSbot2(1) + 27*MStop2(1)) + 
     &                     C2B*
     &                      (-2*MW2**2 - 6*MStop2(1)**2 + 
     &                       MSbot2(1)*
     &                       (-6*MSbot2(1) + 78*MStop2(1)))))))))

	tmp102 = 6*MBot2 - 6*MTop2 + 325*MW2 - 6*C2B*MW2

	tmp103 = -3*MBot2 + 3*MTop2 - MW2 + 3*C2B*MW2

	tmp105 = -3*MBot2 + 3*MTop2 + 19*MW2 + 3*C2B*MW2

	tmp108 = -6*MBot2 + 6*MTop2 + 11*MW2 + 6*C2B*MW2

	tmp109 = -50*MBot2 + 50*MTop2 + 47*MW2 + 50*C2B*MW2

	tmp111 = Re(Abot)**2 + 
     &   reMUE*(reMUE*TB2 - 2*TB*Re(Abot))

	tmp112 = -12*MBot2 + 12*MTop2 + 3*MW2 + 12*C2B*MW2 - 4*MSq2Diff(1,2)

	tmp117 = 3*(MBot2 - MTop2)**2 - 
     &   3*(23 + 2*C2B)*(MBot2 - MTop2)*MW2 + 
     &   (-16 + C2B*(69 + 3*C2B))*MW2**2

	tmp118 = 3*(MBot2 - MTop2)**2 + 
     &   2*(1 + 3*C2B)*(-MBot2 + MTop2)*MW2 + 
     &   (-6 + C2B*(2 + 3*C2B))*MW2**2

	tmp120 = 6*(MBot2 - MTop2)**2 - 
     &   (19 + 12*C2B)*(MBot2 - MTop2)*MW2 + 
     &   (-3 + C2B*(19 + 6*C2B))*MW2**2

	Sf5806 = (tmp112*(MW2 - MSq2Diff(1,2))**5*
     &       MSq2Diff(1,2) + 
     &      MW2*(MW2 - MSq2Diff(1,2))*
     &       (3*MW2**3*(-4*MBot2 + 4*MTop2 + MW2 + 4*C2B*MW2) - 
     &         MW2**2*(-108*MBot2 + 108*MTop2 + 35*MW2 + 108*C2B*MW2)*
     &          MSq2Diff(1,2) + 
     &         MW2*(-219*MBot2 + 219*MTop2 + 115*MW2 + 219*C2B*MW2)*
     &          MSq2Diff(1,2)**2 - 
     &         7*(-16*MBot2 + 16*MTop2 + 17*MW2 + 16*C2B*MW2)*
     &          MSq2Diff(1,2)**3 + 42*MSq2Diff(1,2)**4)*MStop2(2)+
     &        MW2**2*(-(MW2**2*
     &            (-60*MBot2 + 60*MTop2 + 19*MW2 + 60*C2B*MW2)) + 
     &         2*MW2*(-134*MBot2 + 134*MTop2 + 75*MW2 + 134*C2B*MW2)*
     &          MSq2Diff(1,2) - 6*tmp109*MSq2Diff(1,2)**2 + 
     &         142*MSq2Diff(1,2)**3)*MStop2(2)**2 + 
     &      MW2**3*(-61*MBot2 + 61*MTop2 + 50*MW2 + 61*C2B*MW2 - 
     &         152*MSq2Diff(1,2))*MStop2(2)**3 + 4*MStop2(2)**7)*
     &    USbot11**2 + MBot*(reMUE*TB - Re(Abot))*
     &    (12*(MW2 - MSq2Diff(1,2))**5*MSq2Diff(1,2) + 
     &      MW2*(112*MSq2Diff(1,2)**4 + 
     &         MW2*(MW2**2*(tmp108 - 117*MSq2Diff(1,2)) + 
     &            tmp102*MSq2Diff(1,2)**2 - 331*MSq2Diff(1,2)**3))*
     &       MStop2(2) - 4*MW2**2*
     &       (MW2*tmp105 + 
     &         3*(-MBot2 + MTop2 - 22*MW2 + C2B*MW2)*MSq2Diff(1,2) + 
     &         75*MSq2Diff(1,2)**2)*MStop2(2)**2 + 
     &      128*MW2**3*MStop2(2)**3)*USbot11*USbot12 + 
     &   3*MBot2*MW2**2*tmp111*(MW2 + MSq2Diff(1,2))*
     &    (MW2 - MSbot2(1) - MStop2(2))*MStop2(2)*USbot12**2 + 
     &   MStop2(2)*(MW2**3*(-MBot2 + MTop2 + C2B*MW2)*tmp103*
     &       (MW - MSbot(1))*(MW + MSbot(1)) + 
     &      MW2**2*tmp118*(MW - MSbot(1))*(MW + MSbot(1))*
     &       MSbot2(1) - MW2**3*tmp120*MStop2(2) + 
     &      MW2**2*tmp117*MStop2(2)**2 - 4*MStop2(2)**6)*
     &    USbot22**2

	tmp98 = -6*MBot2 + 6*MTop2 - 89*MStop2(1)

	Sf5705 = ((12*MTop2 + 3*MW2 + 12*C2B*MW2 - 4*MSbot2(1))*
     &       (MW2 - MSbot2(1))**5*MSbot2(1) - 
     &      (MW2 - MSbot2(1))*MSbot2(1)**2*
     &       (33*MW2**2*(MTop2 + C2B*MW2) - 
     &         2*MW2*(58*MTop2 + 23*MW2 + 58*C2B*MW2)*MSbot2(1) + 
     &         4*(18*MTop2 + 17*MW2 + 18*C2B*MW2)*MSbot2(1)**2 - 
     &         28*MSbot2(1)**3)*MStop2(1) - 
     &      MSbot2(1)*(26*MW2**3*(MTop2 + C2B*MW2) + 
     &         27*MW2**2*(MTop2 + C2B*MW2)*MSbot2(1) - 
     &         8*MW2*(19*MTop2 + 6*MW2 + 19*C2B*MW2)*MSbot2(1)**2 + 
     &         45*(4*MTop2 + 3*MW2 + 4*C2B*MW2)*MSbot2(1)**3 - 
     &         84*MSbot2(1)**4)*MStop2(1)**2 + 
     &      MSbot2(1)*(MW2**2*(87*MTop2 - 10*MW2 + 87*C2B*MW2) + 
     &         2*MW2*(36*MTop2 - 5*MW2 + 36*C2B*MW2)*MSbot2(1) + 
     &         40*(6*MTop2 + MW2 + 6*C2B*MW2)*MSbot2(1)**2 - 
     &         140*MSbot2(1)**3)*MStop2(1)**3 - 
     &      (MW2**2*(89*MTop2 - 34*MW2 + 89*C2B*MW2) + 
     &         MW2*(148*MTop2 - 57*MW2 + 148*C2B*MW2)*MSbot2(1) + 
     &         15*(12*MTop2 - 5*MW2 + 12*C2B*MW2)*MSbot2(1)**2 - 
     &         140*MSbot2(1)**3)*MStop2(1)**4 + 
     &      4*(MW2*(13*MTop2 - 9*MW2 + 13*C2B*MW2) + 
     &         18*(MTop2 - MW2 + C2B*MW2)*MSbot2(1) - 21*MSbot2(1)**2
     &         )*MStop2(1)**5 - 
     &      (12*MTop2 - 19*MW2 + 12*C2B*MW2 - 28*MSbot2(1))*
     &       MStop2(1)**6 + 
     &      MBot2*(-12*MW2**5*MSbot2(1) + 60*MW2**4*MSbot2(1)**2 + 
     &         12*MSq2Diff(1,1)**6 - 
     &         4*MW2*MSq2Diff(1,1)**4*
     &          (15*MSbot2(1) + 13*MStop2(1)) + 
     &         MW2**3*MSbot2(1)*
     &          (26*MStop2(1)**2 + 
     &            MSbot2(1)*(-120*MSbot2(1) + 33*MStop2(1))) + 
     &         MW2**2*MSq2Diff(1,1)**2*
     &          (89*MStop2(1)**2 + 
     &            MSbot2(1)*(120*MSbot2(1) + 91*MStop2(1)))))*
     &    USbot11**2 + MBot*(reMUE*TB - Re(Abot))*
     &    (-12*MSq2Diff(1,1)**6 + 
     &      4*MW2*MSq2Diff(1,1)**4*(15*MSbot2(1) + 13*MStop2(1)) + 
     &      MW2**5*(12*MSbot2(1) + (-1 + 6*C2B)*MStop2(1)) + 
     &      MW2**3*(120*MSbot2(1)**3 - 
     &         (35 + 6*C2B)*MSbot2(1)**2*MStop2(1) + 
     &         2*(6*MBot2 - 6*MTop2 - 13*MSbot2(1))*MStop2(1)**2 + 
     &         (69 + 6*C2B)*MStop2(1)**3) - 
     &      MW2**4*(MSbot2(1)*(60*MSbot2(1) - 3*MStop2(1)) + 
     &         6*MBot2*MStop2(1) + 
     &         MStop2(1)*(-6*MTop2 + (19 + 12*C2B)*MStop2(1))) - 
     &      MW2**2*MSq2Diff(1,1)*
     &       (tmp98*MStop2(1)**2 + 
     &         MSbot2(1)*(MSbot2(1)*
     &             (120*MSbot2(1) - 29*MStop2(1)) - 
     &            6*MBot2*MStop2(1) + (6*MTop2 - 2*MStop2(1))*MStop2(1)
     &            )))*USbot11*USbot12 + 
     &   3*MBot2*MW2**2*tmp111*(MW2 + MSq2Diff(1,1))*
     &    (MW2 - MSbot2(1) - MStop2(1))*MStop2(1)*USbot12**2

	Sf5705 = Sf5705 + 
     &   MStop2(1)*(MW2**3*(-MBot2 + MTop2 + C2B*MW2)*tmp103*
     &       (MW - MSbot(1))*(MW + MSbot(1)) + 
     &      MW2**2*tmp118*(MW - MSbot(1))*(MW + MSbot(1))*
     &       MSbot2(1) - MW2**3*tmp120*MStop2(1) + 
     &      MW2**2*tmp117*MStop2(1)**2 - 4*MStop2(1)**6)*
     &    USbot22**2

	tmp113 = -12*MBot2 + 12*MTop2 + 3*MW2 + 12*C2B*MW2 - 4*MSq2Diff(2,2)

	tmp119 = 3*(MBot2 - MTop2)**2 - 
     &   (11 + 6*C2B)*(MBot2 - MTop2)*MW2 + 
     &   (3 + C2B*(11 + 3*C2B))*MW2**2

	Sf6184 = -4*MStop2(2)**7 + 
     &   (tmp113*(MW2 - MSq2Diff(2,2))**5*MSq2Diff(2,2) + 
     &      MW2*(MW2**3*tmp119 - 
     &         MW2**3*(-117*MBot2 + 117*MTop2 + 44*MW2 + 117*C2B*MW2)*
     &          MSq2Diff(2,2) + 
     &         MW2*(-3*(MBot2 - MTop2)**2 + 
     &            (-325 + 6*C2B)*(MBot2 - MTop2)*MW2 - 
     &            (-156 + C2B*(-325 + 3*C2B))*MW2**2)*
     &          MSq2Diff(2,2)**2 - 
     &         MW2*(-331*MBot2 + 331*MTop2 + 234*MW2 + 331*C2B*MW2)*
     &          MSq2Diff(2,2)**3 + 
     &         7*(-16*MBot2 + 16*MTop2 + 23*MW2 + 16*C2B*MW2)*
     &          MSq2Diff(2,2)**4 - 42*MSq2Diff(2,2)**5)*MStop2(2)-
     &        2*MW2**2*(MW2*
     &          (3*(MBot2 - MTop2)**2 - 
     &            2*(19 + 3*C2B)*(MBot2 - MTop2)*MW2 + 
     &            (11 + C2B*(38 + 3*C2B))*MW2**2) + 
     &         3*((MBot2 - MTop2)**2 - 2*(-22 + C2B)*(MBot2 - MTop2)*MW2 + 
     &            (-27 + (-44 + C2B)*C2B)*MW2**2)*MSq2Diff(2,2) + 
     &         3*tmp109*MSq2Diff(2,2)**2 - 71*MSq2Diff(2,2)**3)*
     &       MStop2(2)**2 + 
     &      8*MW2**3*(-16*MBot2 + 16*MTop2 + 5*MW2 + 16*C2B*MW2 - 
     &         19*MSq2Diff(2,2))*MStop2(2)**3 + 4*MStop2(2)**7)*
     &    USbot21**2 + MBot*(reMUE*TB - Re(Abot))*
     &    (12*(MW2 - MSq2Diff(2,2))**5*MSq2Diff(2,2) + 
     &      MW2*(112*MSq2Diff(2,2)**4 + 
     &         MW2*(MW2**2*(tmp108 - 117*MSq2Diff(2,2)) + 
     &            tmp102*MSq2Diff(2,2)**2 - 331*MSq2Diff(2,2)**3))*
     &       MStop2(2) - 4*MW2**2*
     &       (MW2*tmp105 + 
     &         3*(-MBot2 + MTop2 - 22*MW2 + C2B*MW2)*MSq2Diff(2,2) + 
     &         75*MSq2Diff(2,2)**2)*MStop2(2)**2 + 
     &      128*MW2**3*MStop2(2)**3)*USbot21*USbot22 + 
     &   MStop2(2)*(3*MBot2*MW2**2*tmp111*(MW2 + MSq2Diff(2,2))*
     &       (MW2 - MSbot2(2) - MStop2(2)) + 4*MStop2(2)**6)*
     &    USbot22**2

	tmp106 = 4*MBot2 - 4*MTop2 - MW2 - 4*C2B*MW2

	Sf6094 = -4*MStop2(1)**7 + 
     &   (-3*tmp106*(MW2 - MSbot2(2))**5*MSbot2(2) - 
     &      4*(MW2 - MSbot2(2))**5*MSbot2(2)**2 + 
     &      MW2**3*(-MBot2 + MTop2 + C2B*MW2)*tmp103*(MW - MSbot(2))*
     &       (MW + MSbot(2))*MStop2(1) + 
     &      MW2**2*tmp118*(MW - MSbot(2))*(MW + MSbot(2))*
     &       MSbot2(2)*MStop2(1) - 
     &      33*MW2**2*(-MBot2 + MTop2 + C2B*MW2)*(MW - MSbot(2))*
     &       (MW + MSbot(2))*MSbot2(2)**2*MStop2(1) + 
     &      2*MW2*(-58*MBot2 + 58*MTop2 + 23*MW2 + 58*C2B*MW2)*
     &       (MW - MSbot(2))*(MW + MSbot(2))*MSbot2(2)**3*MStop2(1)
     &        + 4*(18*MBot2 - 18*MTop2 - 17*MW2 - 18*C2B*MW2)*
     &       (MW - MSbot(2))*(MW + MSbot(2))*MSbot2(2)**4*MStop2(1)
     &        + 28*(MW - MSbot(2))*(MW + MSbot(2))*MSbot2(2)**5*
     &       MStop2(1) - MW2**3*tmp120*MStop2(1)**2 + 
     &      8*MW2*(-19*MBot2 + 19*MTop2 + 6*MW2 + 19*C2B*MW2)*
     &       MSbot2(2)**3*MStop2(1)**2 + 
     &      45*(4*MBot2 - 4*MTop2 - 3*MW2 - 4*C2B*MW2)*MSbot2(2)**4*
     &       MStop2(1)**2 + 84*MSbot2(2)**5*MStop2(1)**2 - 
     &      MW2**2*(-MBot2 + MTop2 + C2B*MW2)*MSbot2(2)*
     &       (26*MW2 + 27*MSbot2(2))*MStop2(1)**2 + 
     &      MW2**2*tmp117*MStop2(1)**3 + 
     &      MW2**2*(-87*MBot2 + 87*MTop2 - 10*MW2 + 87*C2B*MW2)*
     &       MSbot2(2)*MStop2(1)**3 + 
     &      2*MW2*(-36*MBot2 + 36*MTop2 - 5*MW2 + 36*C2B*MW2)*
     &       MSbot2(2)**2*MStop2(1)**3 - 
     &      40*(6*MBot2 - 6*MTop2 - MW2 - 6*C2B*MW2)*MSbot2(2)**3*
     &       MStop2(1)**3 - 140*MSbot2(2)**4*MStop2(1)**3 + 
     &      MW2**2*(89*MBot2 - 89*MTop2 + 34*MW2 - 89*C2B*MW2)*
     &       MStop2(1)**4 + 
     &      MW2*(148*MBot2 - 148*MTop2 + 57*MW2 - 148*C2B*MW2)*
     &       MSbot2(2)*MStop2(1)**4 + 
     &      15*(12*MBot2 - 12*MTop2 + 5*MW2 - 12*C2B*MW2)*MSbot2(2)**2*
     &       MStop2(1)**4 + 140*MSbot2(2)**3*MStop2(1)**4 + 
     &      4*MW2*(-13*MBot2 + 13*MTop2 - 9*MW2 + 13*C2B*MW2)*
     &       MStop2(1)**5 - 
     &      72*(MBot2 - MTop2 + MW2 - C2B*MW2)*MSbot2(2)*
     &       MStop2(1)**5 - 84*MSbot2(2)**2*MStop2(1)**5 + 
     &      (12*MBot2 - 12*MTop2 + 19*MW2 - 12*C2B*MW2 + 28*MSbot2(2))*
     &       MStop2(1)**6)*USbot21**2

	Sf6094 = Sf6094 + 
     &   MBot*(reMUE*TB - Re(Abot))*
     &    (-12*MSq2Diff(2,1)**6 + 
     &      4*MW2*MSq2Diff(2,1)**4*(15*MSbot2(2) + 13*MStop2(1)) + 
     &      MW2**5*(12*MSbot2(2) + (-1 + 6*C2B)*MStop2(1)) + 
     &      MW2**3*(120*MSbot2(2)**3 - 
     &         (35 + 6*C2B)*MSbot2(2)**2*MStop2(1) + 
     &         2*(6*MBot2 - 6*MTop2 - 13*MSbot2(2))*MStop2(1)**2 + 
     &         (69 + 6*C2B)*MStop2(1)**3) - 
     &      MW2**4*(MSbot2(2)*(60*MSbot2(2) - 3*MStop2(1)) + 
     &         6*MBot2*MStop2(1) + 
     &         MStop2(1)*(-6*MTop2 + (19 + 12*C2B)*MStop2(1))) - 
     &      MW2**2*MSq2Diff(2,1)*
     &       (tmp98*MStop2(1)**2 + 
     &         MSbot2(2)*(MSbot2(2)*
     &             (120*MSbot2(2) - 29*MStop2(1)) - 
     &            6*MBot2*MStop2(1) + (6*MTop2 - 2*MStop2(1))*MStop2(1)
     &            )))*USbot21*USbot22 + 
     &   MStop2(1)*(3*MBot2*MW2**2*tmp111*(MW2 + MSq2Diff(2,1))*
     &       (MW2 - MSbot2(2) - MStop2(1)) + 4*MStop2(1)**6)*
     &    USbot22**2

	Sf8851 = MZ2*(-4*Sf8642 + 
     &     MZ2**2*(Sf8575 + 
     &        MZ2*(60*Sf8819*Ss443**2*UStau22**2 + 
     &           MZ2*(MZ2*(20*MZ2*Sf8846 - 3*Sf8850) + 
     &              240*Sf8820*Ss443*SW2*UStau21*UStau22))))

	Sf8818 = 6*Sf7769 + 6*Sf7788 + 2*Sf8503 + 2*Sf8519 - 
     &   2*Sf8740 - 2*Sf8762 + 
     &   3*MZ2**3*(-2*Sf7643*MSneu2(1) + 
     &      MZ2*(Sf8780*Sf8782*UStau21**2 + 
     &         Sf8773*Sf8775*UStau22**2 + 
     &         3*(Sf8810*Sf8812*USbot21**2 + 
     &            Sf8804*Sf8806*USbot22**2)*UStop21**2 + 
     &         3*(Sf8797*Sf8799*USbot21**2 + 
     &            Sf8791*Sf8793*USbot22**2)*UStop22**2))

	Sf2179 = (2*Sf1922)/Sf166**3 + (2*Sf1982)/Sf1195**3 - 
     &   (MW2*((Sf162*Sf1936)/Sf166**2 + 
     &        (Sf1191*Sf1994)/Sf1195**2)*TB)/(Ss241*Ss242) + 
     &   (MW2*Sf1191*Sf2005*Ss241*Ss242*TB)/(Sf1195*Ss25**3) + 
     &   (MW2*(-Sf2178 + (CW2*Sf162*Sf1956*Ss241*Ss242)/Sf166)*TB)/
     &    (CW2*Ss25**3)

	Sf351 = -(Sf349*UStau11) + 
     &   MW2*(3*MW2*Sf350 + 
     &      (4*MW2**2*Sf319 + 3*Sf348 + 
     &         MW2**3*(49*MW2**2 + 6*Sf317)*MSl2Diff(1,1) - 
     &         348*MW2**4*Sf34*MSl2Diff(1,1))*UStau22**2)

	Sf8872 = 6*Sf8818 + 
     &   MW2*(-6*Sf8851 + MW2*
     &       (MW2*(MW2*(12*MZ2*Sf8827 + 
     &               MW2*(128*MW2*MZ2*Sf8835 - 128*MZ2**2*Sf8844))+
     &              12*MZ2**2*Sf8857) - Sf8871))

	Sf5090 = -2*Sf4867*MHiggs2(4)*MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(MW2*
     &             (MW2*(4*MW2*Sf5058 - 4*Sf5060) + Sf5063) - 
     &            Sf5072) + Sf5083) - Sf5089*MSq2Diff(2,2))

	Sf4964 = -2*Sf4831*MHiggs2(4)*MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(MW2*
     &             (MW2*(4*MW2*Sf4932 - 4*Sf4934) + Sf4937) - 
     &            Sf4946) + Sf4957) - Sf4963*MSq2Diff(1,2))

	Sf2650 = -2*Sf2270*MHiggs2(4)*MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(MW2*
     &             (MW2*(4*MW2*Sf2604 - 4*Sf2607) + Sf2612) - 
     &            Sf2623) + Sf2643) - Sf2649*MSq2Diff(2,1))

	Sf2436 = -2*Sf2206*MHiggs2(4)*MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(MW2*
     &             (MW2*(4*MW2*Sf2366 - 4*Sf2369) + Sf2382) - 
     &            Sf2401) + Sf2429) - Sf2435*MSq2Diff(1,1))

	Sf6921 = 120*MSq2Diff(2,2)**6 + 
     &   MW2*(-12*Sf6888*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf6906) + Sf6916) - 
     &            3*Sf6920) + 3*Sf6917*MSq2Diff(2,2)))

	Sf6670 = 120*MSq2Diff(1,2)**6 + 
     &   MW2*(-12*Sf6637*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf6655) + Sf6665) - 
     &            3*Sf6669) + 3*Sf6666*MSq2Diff(1,2)))

	Sf4733 = 120*MSq2Diff(2,1)**6 + 
     &   MW2*(-12*Sf4683*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf4710) + Sf4724) - 
     &            3*Sf4732) + 3*Sf4727*MSq2Diff(2,1)))

	Sf4638 = 120*MSq2Diff(1,1)**6 + 
     &   MW2*(-12*Sf4578*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(36*MW2 - 3*Sf4613) + Sf4629) - 
     &            3*Sf4637) + 3*Sf4632*MSq2Diff(1,1)))

	tmp101 = Sf4828*USbot11 + Sf4864*USbot21

	Sf7154 = MW2*Sf4828*Sf6325 + 
     &   (MW2*Sf4820**2*Ss25*TB*tmp101*MHiggs2(4) + 
     &      Sf4905*Sf6297*Ss25**3*USbot11)*UStop21

	tmp100 = Sf2198*USbot11 + Sf2263*USbot21

	Sf7139 = MW2*Sf2198*Sf4127 + 
     &   (MW2*Sf2320**2*Ss25*TB*tmp100*MHiggs2(4) + 
     &      Sf2327*Sf3999*Ss25**3*USbot11)*UStop11

	Sf1665 = MHiggs2(4)*
     &    (Sf1664 + 9*MW2*Sf1183*Sf1184*Sf1195**3*TB*MHiggs2(4)**2*
     &       MSl2Diff(2,1)) - 
     &   18*CW2*MW2*Sf1191*Sf1195*Sf1570*Ss25*TB*MSneu2(1)

	Sf848 = MHiggs2(4)*
     &    (Sf847 + 9*MW2*Sf153*Sf154*Sf166**3*TB*MHiggs2(4)**2*
     &       MSl2Diff(1,1)) - 
     &   18*CW2*MW2*Sf162*Sf166*Sf702*Ss25*TB*MSneu2(1)

	Sf2015 = Sf2013*Ss25**3 + 
     &   Sf1195*Ss241*Ss242*
     &    (-Sf2014 + 4*Sf1195*Sf1923*MSl2Diff(2,1)*UStau21)

	Sf8484 = -(Sf2334**3*Sf4820**3*Sf4906**3*Sf8473*
     &      USbot22**2*UStop22**2) + 
     &   Sf2320*(Sf2320**2*Sf8483 - 
     &      Sf2334**3*Sf4820**3*Sf4906**3*Sf8472*USbot22**2*
     &       UStop22**2)

	Sf7288 = 24*Sf7283 + MZ2**2*(Sf7285 + 2*MZ2*Sf7287*SW2)

	Sf7127 = Sf4864*Sf7034 + 
     &   Sf173*Sf4906**2*Ss25*TB*tmp101*UStop21

	Sf6936 = Sf2263*(Sf6733 - 
     &      Sf173*Sf2334**2*Ss25*TB*USbot21*UStop11) - 
     &   Sf173*Sf2334**2*Sf4864*Ss25*TB*USbot21*UStop21

	Sf6686 = Sf2198*(Sf6466 - 
     &      Sf173*Sf2320**2*Ss25*TB*USbot11*UStop11) - 
     &   Sf173*Sf2320**2*Sf4828*Ss25*TB*USbot11*UStop21

	Sf4756 = Sf2263*Sf4527 + 
     &   Sf173*Sf2334**2*Ss25*TB*tmp100*UStop11

	Sf6294 = Sf6286*TB2*UStop21**2 + 
     &   Sf2245*Sf6293*USbot21*UStop22

	Sf6280 = Sf6272*TB2*UStop21**2 + 
     &   Sf2245*Sf6279*USbot11*UStop22

	Sf6188 = Sf6184*TB2*UStop21**2 + 
     &   Sf2245*Sf6187*USbot21*UStop22

	Sf5810 = Sf5806*TB2*UStop21**2 + 
     &   Sf2245*Sf5809*USbot11*UStop22

	Sf7128 = MW2*Sf7127*MSq2Diff(2,2) + 
     &   Sf4911*Sf7126*Ss25**3*USbot21*UStop21

	Sf6360 = -(Sf4831*Sf4864*Sf6359*MSq2Diff(1,2)) + 
     &   2*Sf4828*Sf4867*Sf4906**2*Sf6357*Ss241**2*Ss242**2*TB*
     &    MSq2Diff(2,2)*USbot11*UStop21

	Sf6924 = Sf6914*TB2*UStop21 - 
     &   MTop*Sf2245*Sf6921*TB*USbot21*UStop22

	Sf6673 = Sf6663*TB2*UStop21 - 
     &   MTop*Sf2245*Sf6670*TB*USbot11*UStop22

	Sf5406 = Sf5405*TB2*UStop21 - 
     &   9*MTop*Sf2245*Sf5404*TB*USbot21*UStop22

	Sf5300 = Sf5299*TB2*UStop21 - 
     &   9*MTop*Sf2245*Sf5298*TB*USbot11*UStop22

	Sf6106 = Sf2245*Sf6105*USbot21*UStop12 + 
     &   Sf6094*TB2*UStop22**2

	Sf5717 = Sf2245*Sf5716*USbot11*UStop12 + 
     &   Sf5705*TB2*UStop22**2

	Sf3992 = Sf2245*Sf3991*USbot21*UStop12 + 
     &   Sf3975*TB2*UStop22**2

	Sf3895 = Sf2245*Sf3894*USbot11*UStop12 + 
     &   Sf3878*TB2*UStop22**2

	Sf4757 = MW2*Sf4756*MSq2Diff(2,1) + 
     &   Sf2340*Sf4755*Ss25**3*USbot21*UStop11

	Sf4737 = Sf4721*TB2*UStop11 - 
     &   MTop*Sf2245*Sf4733*TB*USbot21*UStop12

	Sf4642 = Sf4626*TB2*UStop11 - 
     &   MTop*Sf2245*Sf4638*TB*USbot11*UStop12

	Sf2998 = Sf2997*TB2*UStop11 - 
     &   9*MTop*Sf2245*Sf2996*TB*USbot21*UStop12

	Sf2860 = Sf2859*TB2*UStop11 - 
     &   9*MTop*Sf2245*Sf2858*TB*USbot11*UStop12

	Sf1818 = Sf1808*UStau21**2 - Sf1817*Sf215*UStau22

	Sf2172 = Sf2162*UStau21 + 3*ML*Sf215*Sf2167*UStau22

	Sf1861 = Sf1860*UStau21 + 3*ML*Sf1859*Sf215*UStau22

	Sf1007 = -(Sf1006*Sf215*UStau12) + Sf996*UStau22**2

	Sf8485 = Sf166**3*Sf8484 + 
     &   Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf8471*
     &    UStau22**2

	Sf6417 = -(Sf154*Sf6416) + 
     &   Sf162*Sf2206*Sf2270*Sf4831*Sf4867*Sf6407*TB*UStau11

	Sf2103 = Sf2093*UStau11 + 3*ML*Sf2098*Sf215*UStau12

	Sf1088 = Sf1087*UStau11 + 3*ML*Sf1086*Sf215*UStau12

	Sf7105 = Sf7103*USbot21 - 3*MBot*Sf2225*Sf7104*USbot22

	Sf6879 = Sf6877*USbot21 + MBot*Sf2225*Sf6878*USbot22

	Sf6765 = Sf6763*USbot21 + MBot*Sf2225*Sf6764*USbot22

	Sf6173 = Sf6166*USbot21 + MBot*Sf2225*Sf6172*USbot22

	Sf5986 = Sf5982*USbot21 - 18*MBot*Sf2225*Sf5985*USbot22

	Sf5930 = Sf5927*USbot21 + 6*MBot*Sf2225*Sf5929*USbot22

	Sf5370 = Sf5368*USbot21 - 9*MBot*Sf2225*Sf5369*USbot22

	Sf4674 = Sf4672*USbot21 - 3*MBot*Sf2225*Sf4673*USbot22

	Sf3927 = Sf3918*USbot21 + MBot*Sf2225*Sf3926*USbot22

	Sf3574 = Sf3557*USbot21 - 18*MBot*Sf2225*Sf3573*USbot22

	Sf3486 = Sf3483*USbot21 + 6*MBot*Sf2225*Sf3485*USbot22

	Sf2959 = Sf2957*USbot21 - 9*MBot*Sf2225*Sf2958*USbot22

	Sf7067 = Sf7063*USbot11 - 3*MBot*Sf2225*Sf7066*USbot12

	Sf6628 = Sf6626*USbot11 + MBot*Sf2225*Sf6627*USbot12

	Sf6510 = Sf6508*USbot11 + MBot*Sf2225*Sf6509*USbot12

	Sf5795 = Sf5788*USbot11 + MBot*Sf2225*Sf5794*USbot12

	Sf5560 = Sf5556*USbot11 - 18*MBot*Sf2225*Sf5559*USbot12

	Sf5497 = Sf5494*USbot11 + 6*MBot*Sf2225*Sf5496*USbot12

	Sf5262 = Sf5260*USbot11 - 9*MBot*Sf2225*Sf5261*USbot12

	Sf4569 = Sf4567*USbot11 - 3*MBot*Sf2225*Sf4568*USbot12

	Sf3803 = Sf3793*USbot11 + MBot*Sf2225*Sf3802*USbot12

	Sf3222 = Sf3201*USbot11 - 18*MBot*Sf2225*Sf3221*USbot12

	Sf3103 = Sf3100*USbot11 + 6*MBot*Sf2225*Sf3102*USbot12

	Sf2806 = Sf2804*USbot11 - 9*MBot*Sf2225*Sf2805*USbot12

	Sf6264 = Sf6230*Ss25**3*MStop2(2) + 
     &   2*Sf4906**2*Sf6263*Ss241**2*Ss242**2*TB*USbot21*
     &    UStop21

	Sf5900 = Sf5864*Ss25**3*MStop2(2) + 
     &   2*Sf4820**2*Sf5899*Ss241**2*Ss242**2*TB*USbot11*
     &    UStop21

	Sf6381 = Sf3895*TB*MStop2(1) + 
     &   Sf2320*Sf2327*Sf6380*TB2*USbot11*UStop11

	Sf6261 = Sf4171*Ss25**3*MStop2(1) + 
     &   2*Sf2334**2*Sf6260*Ss241**2*Ss242**2*TB*USbot21*
     &    UStop11

	Sf5897 = Sf4086*Ss25**3*MStop2(1) + 
     &   2*Sf2320**2*Sf5896*Ss241**2*Ss242**2*TB*USbot11*
     &    UStop11

	Sf2007 = Sf2006*MSl2Diff(1,1) + 
     &   Sf1146*Sf1184*Sf154*Ss25*MSl2Diff(2,1)*MStau2(1)

	Sf9478 = ((MZ2*Sf9417)/MStau2(1) - 
     &     (4*Sf9477*MStau2(3))/
     &      (Sf9317**2*Sf9331**2*Ss49*MHiggs2(2)*MStau2(2)))/MZ2

	Sf6174 = 3*MW2**2*Sf2961*Sf96**2*MSq2Diff(2,2) + 
     &   Sf6173*USbot21

	Sf6925 = 3*MW2**2*Sf5443*Sf6923*MSq2Diff(2,2) - 
     &   Sf6924*USbot21*UStop21

	Sf6933 = MW2*Sf4864*Sf6847*MSq2Diff(2,2) - 
     &   Sf4911*Sf6932*Ss25**3*USbot21*UStop21

	Sf6823 = 3*MW2*Sf5052*Sf5053*Sf5054*Sf5055*Sf5056*Sf6822*
     &    MHiggs2(4)*MSq2Diff(2,2) - 
     &   Sf6820*TB*USbot21*UStop21

	Sf6354 = MW2*Sf4864*Sf6348*MSq2Diff(2,2) + 
     &   Sf4911*Sf6353*Ss25**3*USbot21*UStop21

	Sf6265 = Sf2270*Sf2334**2*Sf4864*Sf6264*MSq2Diff(2,1) + 
     &   Sf2263*Sf4867*Sf4906**2*Sf6261*MSq2Diff(2,2)

	Sf5092 = 3*MW2*Sf5052*Sf5053*Sf5054*Sf5055*Sf5056*Sf5091*
     &    MSq2Diff(2,2) + Sf5090*TB*USbot21*UStop21

	Sf3928 = 3*MW2**2*Sf2961*Sf74**2*MSq2Diff(2,1) + 
     &   Sf3927*USbot21

	Sf4738 = 3*MW2**2*Sf3053*Sf4736*MSq2Diff(2,1) - 
     &   Sf4737*USbot21*UStop11

	Sf4486 = 3*MW2*Sf2598*Sf2599*Sf2600*Sf2601*Sf2602*Sf4485*
     &    MHiggs2(4)*MSq2Diff(2,1) - 
     &   Sf4481*TB*USbot21*UStop11

	Sf4220 = MW2*Sf2263*Sf4207*MSq2Diff(2,1) + 
     &   Sf2340*Sf4219*Ss25**3*USbot21*UStop11

	Sf2652 = 3*MW2*Sf2598*Sf2599*Sf2600*Sf2601*Sf2602*Sf2651*
     &    MSq2Diff(2,1) + Sf2650*TB*USbot21*UStop11

	Sf5796 = 3*MW2**2*Sf2809*Sf88**2*MSq2Diff(1,2) + 
     &   Sf5795*USbot11

	Sf6674 = 3*MW2**2*Sf5337*Sf6672*MSq2Diff(1,2) - 
     &   Sf6673*USbot11*UStop21

	Sf6682 = MW2*Sf4828*Sf6592*MSq2Diff(1,2) - 
     &   Sf4905*Sf6681*Ss25**3*USbot11*UStop21

	Sf6568 = 3*MW2*Sf4926*Sf4927*Sf4928*Sf4929*Sf4930*Sf6567*
     &    MHiggs2(4)*MSq2Diff(1,2) - 
     &   Sf6565*TB*USbot11*UStop21

	Sf6352 = MW2*Sf4828*Sf6325*MSq2Diff(1,2) + 
     &   Sf4905*Sf6351*Ss25**3*USbot11*UStop21

	Sf5901 = Sf2206*Sf2320**2*Sf4828*Sf5900*MSq2Diff(1,1) + 
     &   Sf2198*Sf4820**2*Sf4831*Sf5897*MSq2Diff(1,2)

	Sf4966 = 3*MW2*Sf4926*Sf4927*Sf4928*Sf4929*Sf4930*Sf4965*
     &    MSq2Diff(1,2) + Sf4964*TB*USbot11*UStop21

	Sf3804 = 3*MW2**2*Sf2809*Sf66**2*MSq2Diff(1,1) + 
     &   Sf3803*USbot11

	Sf4643 = 3*MW2**2*Sf2926*Sf4641*MSq2Diff(1,1) - 
     &   Sf4642*USbot11*UStop11

	Sf4752 = MW2*Sf2198*Sf4423*MSq2Diff(1,1) + 
     &   Sf2327*Sf4751*Ss25**3*USbot11*UStop11

	Sf4375 = 3*MW2*Sf2360*Sf2361*Sf2362*Sf2363*Sf2364*Sf4374*
     &    MHiggs2(4)*MSq2Diff(1,1) - 
     &   Sf4370*TB*USbot11*UStop11

	Sf4215 = Sf2206*Sf2263*Sf2320**2*Sf4214*MSq2Diff(1,1) + 
     &   Sf2198*Sf2270*Sf2334**2*Sf4211*MSq2Diff(2,1)

	Sf4218 = MW2*Sf2198*Sf4127*MSq2Diff(1,1) + 
     &   Sf2327*Sf4217*Ss25**3*USbot11*UStop11

	Sf2438 = 3*MW2*Sf2360*Sf2361*Sf2362*Sf2363*Sf2364*Sf2437*
     &    MSq2Diff(1,1) + Sf2436*TB*USbot11*UStop11

	Sf6364 = Sf1982*MSneu2(1) - Sf1818*MStau2(2)

	Sf6366 = Sf1922*MSneu2(1) - Sf1007*MStau2(1)

	Sf1179 = Sf1178*Ss25**2 + Sf1176*MSneu2(1)

	Sf1862 = 3*MW2**2*Sf1298*Sf44**2*MSl2Diff(2,1) + 
     &   Sf1861*UStau21

	Sf2009 = Sf1982*MSl2Diff(2,1) + Sf1862*MStau2(2)

	Sf1888 = Sf1818*MSl2Diff(2,1) - Sf1862*MSneu2(1)

	Sf1089 = 3*MW2**2*Sf31**2*Sf347*MSl2Diff(1,1) + 
     &   Sf1088*UStau11

	Sf2016 = Sf1191*Sf166**2*Sf2015*MSl2Diff(1,1) - 
     &   Sf1195**2*Sf162*Sf2012*MSl2Diff(2,1)

	Sf2008 = Sf1922*MSl2Diff(1,1) + Sf1089*MStau2(1)

	Sf1147 = Sf1007*MSl2Diff(1,1) - Sf1089*MSneu2(1)

	Sf6016 = Sf2722**5*Sf3617*MSbot2(2)**2 + Sf6015*MStop2(2)**3

	Sf3660 = Sf2722**5*Sf3617*MSbot2(2)**2 + Sf3659*MStop2(1)**3

	Sf7169 = Sf6847*MSbot2(2) + Sf7034*MStop2(2)

	Sf7165 = Sf6733*MSbot2(2) + Sf4527*MStop2(1)

	Sf6377 = Sf6250*MSbot2(2) - Sf6348*MStop2(2)

	Sf6376 = Sf6188*MSbot2(2) - Sf6294*MStop2(2)

	Sf6369 = Sf6142*MSbot2(2) - Sf4207*MStop2(1)

	Sf6368 = Sf6106*MSbot2(2) - Sf3992*MStop2(1)

	Sf5596 = Sf2543**5*Sf3292*MSbot2(1)**2 + Sf5595*MStop2(2)**3

	Sf3380 = Sf2543**5*Sf3292*MSbot2(1)**2 + Sf3379*MStop2(1)**3

	Sf7122 = Sf6568*Ss25**3*MSbot2(1) + 
     &   Sf4831*Sf7012*Ss241**2*Ss242**2*MHiggs2(4)*MSq2Diff(1,2)

	Sf4748 = Sf4375*Ss25**3*MSbot2(1) - 
     &   Sf2320**2*Sf4747*Ss241*Ss242*TB*USbot11*UStop11

	Sf7175 = Sf6592*MSbot2(1) + Sf7012*MStop2(2)

	Sf6373 = Sf5886*MSbot2(1) - Sf6325*MStop2(2)

	Sf6372 = Sf5810*MSbot2(1) - Sf6280*MStop2(2)

	Sf7155 = Sf6280*Ss25**3*TB + Sf4820*Sf7154*TB2

	Sf7150 = Sf6106*Ss25**3*TB + Sf2334*Sf7149*TB2

	Sf7145 = Sf5717*Ss25**3*TB + Sf2320*Sf7144*TB2

	Sf7140 = Sf3895*Ss25**3*TB + Sf2320*Sf7139*TB2

	Sf939 = Sf937*Ss25**3 - 27*MW2**2*Sf162*Sf166*Sf938*TB

	Sf1728 = Sf1726*Ss25**3 - 27*MW2**2*Sf1191*Sf1195*Sf1727*TB

	Sf1180 = -(MW2*Sf1179) + Sf1177*Ss25**2

	Sf7293 = MW2**2*Sf7292 - 6*CW2*MZ2**2*Sf7270*SW2

	Sf7181 = Sf2198*Sf7178 - Sf173*Sf2320**2*Sf7180*Ss25*TB

	Sf7177 = MW2*Sf4828*Sf7175 - 
     &   Sf4905*Sf7176*Ss25**3*USbot11*UStop21

	Sf7171 = MW2*Sf4864*Sf7169 - 
     &   Sf4911*Sf7170*Ss25**3*USbot21*UStop21

	Sf7167 = MW2*Sf2263*Sf7165 - 
     &   Sf2340*Sf7166*Ss25**3*USbot21*UStop11

	Sf7124 = MW2*Sf4828*Sf7122 + 
     &   Sf4831*Sf4905*Sf7123*Ss241**2*Ss242**2*Ss25**3*MHiggs2(4)*
     &    USbot11*UStop21

	Sf7107 = 2*MBot*Sf2962*Sf7100 - Sf7106*USbot21

	Sf7069 = 2*MBot*Sf2810*Sf7065 - Sf7068*USbot11

	Sf6881 = 6*MBot*Sf2962*Sf6875 + Sf6880*USbot21

	Sf6927 = Sf6823 - 
     &   2*Sf4906**2*Sf6774*Ss241*Ss242*TB*USbot21*UStop21

	Sf6767 = 6*MBot*Sf2962*Sf6761 + Sf6766*USbot21

	Sf6630 = 6*MBot*Sf2810*Sf6624 + Sf6629*USbot11

	Sf6677 = Sf6568 - 
     &   2*Sf4820**2*Sf6519*Ss241*Ss242*TB*USbot11*UStop21

	Sf6512 = 6*MBot*Sf2810*Sf6506 + Sf6511*USbot11

	Sf6386 = 3*Sf2198*Sf6382 + Sf2320**2*Sf6385

	Sf6379 = MW2*Sf4864*Sf6377 - 
     &   Sf4911*Sf6378*Ss25**3*USbot21*UStop21

	Sf6375 = MW2*Sf4828*Sf6373 - 
     &   Sf4905*Sf6374*Ss25**3*USbot11*UStop21

	Sf6371 = MW2*Sf2263*Sf6369 - 
     &   Sf2340*Sf6370*Ss25**3*USbot21*UStop11

	Sf6361 = Sf4820**2*Sf6360 - 
     &   Sf4828*Sf4867*Sf4906**2*Sf5864*Ss25**3*MSbot2(1)*
     &    MSq2Diff(2,2)

	Sf7153 = MW2*Sf4864*Sf6348 + 
     &   Sf4911*Sf6300*Ss25**3*USbot21*UStop21

	Sf5407 = 9*Sf5372*Sf5400 + Sf5406*USbot21*UStop21

	Sf5301 = 9*Sf5264*Sf5294 + Sf5300*USbot11*UStop21

	Sf5222 = Sf1184*Sf5221 - 
     &   Sf1191*Sf154*Sf2206*Sf2270*Sf4831*Sf4867*Sf5210*TB*
     &    UStau21

	Sf5223 = -(CW2*Sf5222) + 
     &   Sf1184*Sf154*Sf2206*Sf2270*Sf4831*Sf4867*Sf5209*
     &    MHiggs2(4)**2

	Sf4676 = 2*MBot*Sf2962*Sf4669 - Sf4675*USbot21

	Sf4571 = 2*MBot*Sf2810*Sf4564 - Sf4570*USbot11

	Sf4753 = Sf4486 - 
     &   2*Sf2334**2*Sf4269*Ss241*Ss242*TB*USbot21*UStop11

	Sf6676 = Sf4375 - 
     &   2*Sf2320**2*Sf4267*Ss241*Ss242*TB*USbot11*UStop11

	Sf7138 = MW2*Sf2263*Sf4207 + 
     &   Sf2340*Sf4006*Ss25**3*USbot21*UStop11

	Sf395 = 9*Sf347*Sf392 - Sf394*UStau11

	Sf2999 = 9*Sf2964*Sf2992 + Sf2998*USbot21*UStop11

	Sf2861 = 9*Sf2812*Sf2854 + Sf2860*USbot11*UStop11

	Sf2173 = Sf2170*Sf2171 - Sf2172*UStau21

	Sf2104 = Sf2101*Sf2102 - Sf2103*UStau11

	Sf1829 = Sf1196**3*Sf1821 - Sf1828*MSneu2(1)

	Sf1037 = Sf1016*Sf167**3 - Sf1036*MSneu2(1)

	Sf1332 = 9*Sf1298*Sf1329 + Sf1331*UStau21

	Sf2107 = 4*Sf1089 + 
     &   (MW2*Sf166*Sf2106*TB)/
     &    (Ss241*Ss242*Ss25**2*MHiggs2(4)*MStau2(1))

	Sf1181 = (3*MW2**3*Sf1 + Sf1180/Ss25**2)/MW2**3

	tmp104 = -3*MBot2 + 3*MTop2 + 14*MW2 + 3*C2B*MW2

	tmp107 = -6*MBot2 + 6*MTop2 - 13*MW2 + 6*C2B*MW2

	tmp110 = -241*MBot2 + 241*MTop2 + 135*MW2 + 241*C2B*MW2

	Sf6017 = (tmp113*(MW2 - MSq2Diff(2,2))**5*
     &       MSq2Diff(2,2)**2 + 
     &      2*MW2*(MW2 - MSq2Diff(2,2))*
     &       (2*MW2**4*(-MBot2 + MTop2 + C2B*MW2) + 
     &         MW2**2*tmp119*MSq2Diff(2,2) - 
     &         MW2**2*(-73*MBot2 + 73*MTop2 + 27*MW2 + 73*C2B*MW2)*
     &          MSq2Diff(2,2)**2 + 
     &         4*MW2*(-32*MBot2 + 32*MTop2 + 17*MW2 + 32*C2B*MW2)*
     &          MSq2Diff(2,2)**3 + 
     &         (62*MBot2 - 62*MTop2 - 67*MW2 - 62*C2B*MW2)*
     &          MSq2Diff(2,2)**4 + 23*MSq2Diff(2,2)**5)*MStop2(2)+
     &        2*MW2**2*(MW2**2*
     &          (3*(MBot2 - MTop2)**2 - 
     &            (-13 + 6*C2B)*(MBot2 - MTop2)*MW2 + 
     &            (3 + C2B*(-13 + 3*C2B))*MW2**2) - 
     &         3*MW2*(3*(MBot2 - MTop2)**2 - 
     &            2*(14 + 3*C2B)*(MBot2 - MTop2)*MW2 + 
     &            (11 + C2B*(28 + 3*C2B))*MW2**2)*MSq2Diff(2,2) + 
     &         MW2*tmp110*MSq2Diff(2,2)**2 - 
     &         5*(-40*MBot2 + 40*MTop2 + 39*MW2 + 40*C2B*MW2)*
     &          MSq2Diff(2,2)**3 + 90*MSq2Diff(2,2)**4)*
     &       MStop2(2)**2 - 
     &      4*MW2**3*(3*MBot2**2 + 3*MTop2**2 + 
     &         (11 + C2B*(-10 + 3*C2B))*MW2**2 + 
     &         2*MTop2*(-5*MW2 + 3*C2B*MW2 - 45*MSq2Diff(2,2)) - 
     &         10*(4 + 9*C2B)*MW2*MSq2Diff(2,2) + 
     &         65*MSq2Diff(2,2)**2 + 
     &         MBot2*(-6*MTop2 + 10*MW2 - 6*C2B*MW2 + 90*MSq2Diff(2,2))
     &         )*MStop2(2)**3 + 80*MW2**4*MStop2(2)**4)*
     &    USbot21**2 + 2*MBot*(reMUE*TB - Re(Abot))*
     &    (6*(MW2 - MSq2Diff(2,2))**5*MSq2Diff(2,2)**2 + 
     &      MW2*(MW2 - MSq2Diff(2,2))*
     &       (2*MW2**4 + (11 + 6*C2B)*MW2**3*MSq2Diff(2,2) + 
     &         MW2**2*(-6*MBot2 + 6*MTop2 - 73*MSq2Diff(2,2))*
     &          MSq2Diff(2,2) + 128*MW2*MSq2Diff(2,2)**3 - 
     &         62*MSq2Diff(2,2)**4)*MStop2(2) + 
     &      MW2**2*(-200*MSq2Diff(2,2)**3 + 
     &         MW2*(MW2*tmp107 - 6*tmp104*MSq2Diff(2,2) + 
     &            241*MSq2Diff(2,2)**2))*MStop2(2)**2 + 
     &      4*MW2**3*(3*MBot2 - 3*MTop2 + 5*MW2 - 3*C2B*MW2 + 
     &         45*MSq2Diff(2,2))*MStop2(2)**3)*USbot21*
     &    USbot22 + 6*MBot2*MW2**3*tmp111*MSbot2(2)*
     &    (MW2 - MSbot2(2) - MStop2(2))*MStop2(2)*USbot22**2

	Sf5597 = (tmp112*(MW2 - MSq2Diff(1,2))**5*
     &       MSq2Diff(1,2)**2 + 
     &      2*MW2*(2*MW2**5*(-MBot2 + MTop2 + C2B*MW2) + 
     &         MW2**4*(-10*MBot2 + 10*MTop2 + 3*MW2 + 10*C2B*MW2)*
     &          MSq2Diff(1,2) - 
     &         5*MW2**3*(-17*MBot2 + 17*MTop2 + 6*MW2 + 17*C2B*MW2)*
     &          MSq2Diff(1,2)**2 + 
     &         MW2**2*(-201*MBot2 + 201*MTop2 + 95*MW2 + 201*C2B*MW2)*
     &          MSq2Diff(1,2)**3 - 
     &         5*MW2*(-38*MBot2 + 38*MTop2 + 27*MW2 + 38*C2B*MW2)*
     &          MSq2Diff(1,2)**4 + 
     &         (-62*MBot2 + 62*MTop2 + 90*MW2 + 62*C2B*MW2)*
     &          MSq2Diff(1,2)**5 - 23*MSq2Diff(1,2)**6)*MStop2(2)+
     &        MW2**2*(6*MW2**3*
     &          (4*MBot2 - 4*MTop2 + MW2 - 4*C2B*MW2) - 
     &         MW2**2*(-160*MBot2 + 160*MTop2 + 63*MW2 + 160*C2B*MW2)*
     &          MSq2Diff(1,2) + 2*MW2*tmp110*MSq2Diff(1,2)**2 + 
     &         10*(40*MBot2 - 40*MTop2 - 39*MW2 - 40*C2B*MW2)*
     &          MSq2Diff(1,2)**3 + 180*MSq2Diff(1,2)**4)*
     &       MStop2(2)**2 + 
     &      MW2**3*(MW2*(-50*MBot2 + 50*MTop2 - 41*MW2 + 50*C2B*MW2) + 
     &         40*(-9*MBot2 + 9*MTop2 + 4*MW2 + 9*C2B*MW2)*
     &          MSq2Diff(1,2) - 260*MSq2Diff(1,2)**2)*MStop2(2)**3+
     &        80*MW2**4*MStop2(2)**4)*USbot11**2 + 
     &   2*MBot*(reMUE*TB - Re(Abot))*
     &    (6*(MW2 - MSq2Diff(1,2))**5*MSq2Diff(1,2)**2 + 
     &      MW2*(MW2 - MSq2Diff(1,2))*
     &       (2*MW2**4 + (11 + 6*C2B)*MW2**3*MSq2Diff(1,2) + 
     &         MW2**2*(-6*MBot2 + 6*MTop2 - 73*MSq2Diff(1,2))*
     &          MSq2Diff(1,2) + 128*MW2*MSq2Diff(1,2)**3 - 
     &         62*MSq2Diff(1,2)**4)*MStop2(2) + 
     &      MW2**2*(-200*MSq2Diff(1,2)**3 + 
     &         MW2*(MW2*tmp107 - 6*tmp104*MSq2Diff(1,2) + 
     &            241*MSq2Diff(1,2)**2))*MStop2(2)**2 + 
     &      4*MW2**3*(3*MBot2 - 3*MTop2 + 5*MW2 - 3*C2B*MW2 + 
     &         45*MSq2Diff(1,2))*MStop2(2)**3)*USbot11*
     &    USbot12 + 6*MBot2*MW2**3*tmp111*MSbot2(1)*
     &    (MW2 - MSbot2(1) - MStop2(2))*MStop2(2)*USbot12**2 + 
     &   MW2**3*MSbot2(1)*MStop2(2)*
     &    (2*(3*MBot2 - 3*MTop2 + MW2 - 3*C2B*MW2)*
     &       (MBot2 - MTop2 - C2B*MW2)*(MW2 - MSbot2(1)) - 
     &      3*(2*(MBot2 - MTop2)**2 - 4*(1 + C2B)*(MBot2 - MTop2)*MW2 + 
     &         (1 + C2B*(4 + 2*C2B))*MW2**2)*MStop2(2))*
     &    USbot22**2

	tmp121 = -3*tmp106*(MW2 - MSbot2(2))**5*MSbot2(2)**2 - 
     &   4*(MW2 - MSbot2(2))**5*MSbot2(2)**3 + 
     &   4*MW2**5*(-MBot2 + MTop2 + C2B*MW2)*(MW - MSbot(2))*
     &    (MW + MSbot(2))*MStop2(1) + 
     &   2*MW2**3*(-MBot2 + MTop2 + C2B*MW2)*tmp103*(MW - MSbot(2))*
     &    (MW + MSbot(2))*MSbot2(2)*MStop2(1) + 
     &   MW2**3*(-10*MBot2 + 10*MTop2 - 3*MW2 + 10*C2B*MW2)*
     &    (MW - MSbot(2))*(MW + MSbot(2))*MSbot2(2)**2*MStop2(1) - 
     &   MW2**2*(-68*MBot2 + 68*MTop2 + 13*MW2 + 68*C2B*MW2)*
     &    (MW - MSbot(2))*(MW + MSbot(2))*MSbot2(2)**3*MStop2(1) + 
     &   MW2*(-152*MBot2 + 152*MTop2 + 67*MW2 + 152*C2B*MW2)*
     &    (MW - MSbot(2))*(MW + MSbot(2))*MSbot2(2)**4*MStop2(1) + 
     &   (84*MBot2 - 84*MTop2 - 83*MW2 - 84*C2B*MW2)*(MW - MSbot(2))*
     &    (MW + MSbot(2))*MSbot2(2)**5*MStop2(1) + 
     &   32*(MW - MSbot(2))*(MW + MSbot(2))*MSbot2(2)**6*
     &    MStop2(1) + MW2**5*
     &    (32*MBot2 - 32*MTop2 + 3*MW2 - 32*C2B*MW2)*MStop2(1)**2 - 
     &   3*MW2**3*(2*(MBot2 - MTop2)**2 + 
     &      4*(1 + C2B)*(-MBot2 + MTop2)*MW2 + 
     &      (1 + C2B*(4 + 2*C2B))*MW2**2)*MSbot2(2)*MStop2(1)**2 - 
     &   4*MW2**3*(-MBot2 + MTop2 + C2B*MW2)*MSbot2(2)**2*
     &    MStop2(1)**2 - 10*MW2**2*
     &    (-8*MBot2 + 8*MTop2 + MW2 + 8*C2B*MW2)*MSbot2(2)**3*
     &    MStop2(1)**2 + 35*MW2*
     &    (-8*MBot2 + 8*MTop2 + 3*MW2 + 8*C2B*MW2)*MSbot2(2)**4*
     &    MStop2(1)**2 + 9*(28*MBot2 - 28*MTop2 - 23*MW2 - 28*C2B*MW2)*
     &    MSbot2(2)**5*MStop2(1)**2 + 
     &   112*MSbot2(2)**6*MStop2(1)**2 + 
     &   MW2**4*(-100*MBot2 + 100*MTop2 - 19*MW2 + 100*C2B*MW2)*
     &    MStop2(1)**3 + 2*MW2**3*
     &    (-61*MBot2 + 61*MTop2 - 5*MW2 + 61*C2B*MW2)*MSbot2(2)*
     &    MStop2(1)**3 + 10*MW2**2*
     &    (-12*MBot2 + 12*MTop2 - MW2 + 12*C2B*MW2)*MSbot2(2)**2*
     &    MStop2(1)**3 + 20*MW2*(-2*MBot2 + 2*MTop2 - MW2 + 2*C2B*MW2)*
     &    MSbot2(2)**3*MStop2(1)**3 - 
     &   5*(84*MBot2 - 84*MTop2 - 23*MW2 - 84*C2B*MW2)*MSbot2(2)**4*
     &    MStop2(1)**3 - 224*MSbot2(2)**5*MStop2(1)**3 + 
     &   10*MW2**3*(16*MBot2 - 16*MTop2 + 5*MW2 - 16*C2B*MW2)*
     &    MStop2(1)**4 + 40*MW2**2*
     &    (7*MBot2 - 7*MTop2 + 2*MW2 - 7*C2B*MW2)*MSbot2(2)*
     &    MStop2(1)**4 + 5*MW2*
     &    (68*MBot2 - 68*MTop2 + 21*MW2 - 68*C2B*MW2)*MSbot2(2)**2*
     &    MStop2(1)**4 + 5*(84*MBot2 - 84*MTop2 + 23*MW2 - 84*C2B*MW2)*
     &    MSbot2(2)**3*MStop2(1)**4 + 
     &   280*MSbot2(2)**4*MStop2(1)**4 - 
     &   70*MW2**2*(2*MBot2 - 2*MTop2 + MW2 - 2*C2B*MW2)*
     &    MStop2(1)**5 + 10*MW2*
     &    (-26*MBot2 + 26*MTop2 - 15*MW2 + 26*C2B*MW2)*MSbot2(2)*
     &    MStop2(1)**5

	tmp121 = tmp121 - 
     &   9*(28*MBot2 - 28*MTop2 + 23*MW2 - 28*C2B*MW2)*MSbot2(2)**2*
     &    MStop2(1)**5 - 224*MSbot2(2)**3*MStop2(1)**5 + 
     &   MW2*(64*MBot2 - 64*MTop2 + 55*MW2 - 64*C2B*MW2)*
     &    MStop2(1)**6 + (84*MBot2 - 84*MTop2 + 115*MW2 - 84*C2B*MW2)*
     &    MSbot2(2)*MStop2(1)**6 + 112*MSbot2(2)**2*MStop2(1)**6 - 
     &   (12*MBot2 - 12*MTop2 + 23*MW2 - 12*C2B*MW2 + 32*MSbot2(2))*
     &    MStop2(1)**7 + 4*MStop2(1)**8

	Sf3661 = tmp121*USbot21**2 + 
     &   2*MBot*(reMUE*TB - Re(Abot))*
     &    (-6*MSq2Diff(2,1)**7 + 2*MW2**6*MStop2(1) + 
     &      2*MW2*MSq2Diff(2,1)**5*(15*MSbot2(2) + 16*MStop2(1)) + 
     &      MW2**3*(60*MSbot2(2)**4 - 39*MSbot2(2)**3*MStop2(1) - 
     &         80*MStop2(1)**4 - 
     &         2*MSbot2(2)**2*MStop2(1)*
     &          (-3*MBot2 + 3*MTop2 + MStop2(1)) + 
     &         MSbot2(2)*MStop2(1)**2*
     &          (6*MBot2 - 6*MTop2 + 61*MStop2(1))) - 
     &      2*MW2**4*(15*MSbot2(2)**3 + 
     &         3*(-1 + C2B)*MSbot2(2)**2*MStop2(1) - 
     &         25*MStop2(1)**3 + 
     &         3*MSbot2(2)*MStop2(1)*
     &          (MBot2 - MTop2 + (1 + C2B)*MStop2(1))) + 
     &      MW2**5*(-16*MStop2(1)**2 + 
     &         MSbot2(2)*(6*MSbot2(2) + (-3 + 6*C2B)*MStop2(1))) - 
     &      10*MW2**2*MSq2Diff(2,1)**3*
     &       (6*MSq2Diff(2,2)**2 + 7*MStop2(1)**2 + 
     &         7*MStop2(1)*MStop2(2) + 6*MStop2(2)**2 + 
     &         MSq2Diff(2,2)*(7*MStop2(1) + 12*MStop2(2))))*
     &    USbot21*USbot22 + 
     &   6*MBot2*MW2**3*tmp111*MSbot2(2)*
     &    (MW2 - MSbot2(2) - MStop2(1))*MStop2(1)*USbot22**2

	Sf3381 = Sf3204*Sf3205*Sf3376*MStop2(1) - 
     &   Sf3356*MStop2(1)**2 - Sf3380*USbot11

	Sf6387 = -3*Sf6381*Ss25**3 + MW2*Sf2320*Sf6386*TB2 + 
     &   3*Sf5717*Ss25**3*TB*MSbot2(1)

	Sf2180 = (4*Sf1862)/Sf1195**3 + 
     &   (MW2*TB*(Sf1887/MHiggs2(4) + 
     &        (Sf1191*Sf1872)/(Sf1195*MStau2(2))))/Ss25**2 - 
     &   (MW2*Sf1191*Sf2000*TB)/(Sf1195**2*Ss241*Ss242*MStau2(2))

	Sf8486 = Sf1195**3*
     &    (21*MW2*Sf166**3*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*
     &       Sf8469 + Sf8485) + 
     &   Sf166**3*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf8468*
     &    UStau21**2

	Sf6388 = Sf4906**3*Sf6387 + 
     &   Sf2320**3*(3*Sf6376*Ss25**3*TB + 3*Sf4906*Sf6379*TB2)

	Sf6389 = Sf4820**3*Sf6388 + 
     &   Sf2320**3*(3*Sf4906**3*Sf6372*Ss25**3*TB + 
     &      3*Sf4820*Sf4906**3*Sf6375*TB2)

	Sf7156 = Sf4906**3*Sf7155 + 
     &   Sf4820**3*(Sf6294*Ss25**3*TB + Sf4906*Sf7153*TB2)

	Sf7151 = Sf4906**3*Sf7150 + 
     &   Sf2334**3*(Sf6188*Ss25**3*TB + Sf4906*Sf7148*TB2)

	Sf7146 = Sf4820**3*Sf7145 + 
     &   Sf2320**3*(Sf5810*Ss25**3*TB + Sf4820*Sf7143*TB2)

	Sf7141 = Sf2334**3*Sf7140 + 
     &   Sf2320**3*(Sf3992*Ss25**3*TB + Sf2334*Sf7138*TB2)

	Sf7289 = 24*MZ2**3*Sf7278*SW2 + 
     &   MW2*(-(MW2*Sf7288) + 6*MZ2**3*Sf7282*SW2)

	Sf6390 = Sf2334**3*Sf6389 + 
     &   Sf2320**3*(3*Sf4820**3*Sf4906**3*Sf6368*Ss25**3*TB + 
     &      3*Sf2334*Sf4820**3*Sf4906**3*Sf6371*TB2)

	Sf1830 = MW2*(Sf1829 - MW2*Sf1825*MSneu2(1)) + 
     &   Sf1008*Sf1819*UStau21

	Sf1038 = MW2*(Sf1037 - MW2*Sf1031*MSneu2(1)) + 
     &   Sf1008*Sf1009*UStau11

	Sf6183 = Sf6174*TB2*UStop21**2 - 
     &   Sf2245*Sf6182*USbot21*UStop22

	Sf6024 = Sf6017*TB2*UStop21**2 + 
     &   2*Sf2245*Sf6023*USbot21*UStop22

	Sf5805 = Sf5796*TB2*UStop21**2 - 
     &   Sf2245*Sf5804*USbot11*UStop22

	Sf5604 = Sf5597*TB2*UStop21**2 + 
     &   2*Sf2245*Sf5603*USbot11*UStop22

	Sf6011 = Sf5986*TB2*UStop21 + 
     &   18*MTop*Sf2245*Sf5981*TB*USbot21*UStop22

	Sf5956 = Sf5930*TB2*UStop21 - 
     &   6*MTop*Sf2245*Sf5925*TB*USbot21*UStop22

	Sf5591 = Sf5560*TB2*UStop21 + 
     &   18*MTop*Sf2245*Sf5555*TB*USbot11*UStop22

	Sf5530 = Sf5497*TB2*UStop21 - 
     &   6*MTop*Sf2245*Sf5492*TB*USbot11*UStop22

	Sf5373 = Sf5370*TB2*UStop21 + 
     &   9*MTop*Sf2245*Sf5369*TB*USbot21*UStop22

	Sf5265 = Sf5262*TB2*UStop21 + 
     &   9*MTop*Sf2245*Sf5261*TB*USbot11*UStop22

	Sf7108 = Sf7105*TB*UStop21 - 3*MTop*Sf2245*Sf7107*UStop22

	Sf7070 = Sf7067*TB*UStop21 - 3*MTop*Sf2245*Sf7069*UStop22

	Sf6882 = Sf6879*TB*UStop21 - MTop*Sf2245*Sf6881*UStop22

	Sf6631 = Sf6628*TB*UStop21 - MTop*Sf2245*Sf6630*UStop22

	Sf3938 = -(Sf2245*Sf3937*USbot21*UStop12) + 
     &   Sf3928*TB2*UStop22**2

	Sf3814 = -(Sf2245*Sf3813*USbot11*UStop12) + 
     &   Sf3804*TB2*UStop22**2

	Sf3678 = 2*Sf2245*Sf3677*USbot21*UStop12 + 
     &   Sf3661*TB2*UStop22**2

	Sf3398 = 2*Sf2245*Sf3397*USbot11*UStop12 + 
     &   Sf3381*TB2*UStop22**2

	Sf3613 = Sf3574*TB2*UStop11 + 
     &   18*MTop*Sf2245*Sf3556*TB*USbot21*UStop12

	Sf3531 = Sf3486*TB2*UStop11 - 
     &   6*MTop*Sf2245*Sf3481*TB*USbot21*UStop12

	Sf3286 = Sf3222*TB2*UStop11 + 
     &   18*MTop*Sf2245*Sf3200*TB*USbot11*UStop12

	Sf3169 = Sf3103*TB2*UStop11 - 
     &   6*MTop*Sf2245*Sf3098*TB*USbot11*UStop12

	Sf2965 = Sf2959*TB2*UStop11 + 
     &   9*MTop*Sf2245*Sf2958*TB*USbot21*UStop12

	Sf2813 = Sf2806*TB2*UStop11 + 
     &   9*MTop*Sf2245*Sf2805*TB*USbot11*UStop12

	Sf6768 = Sf6765*TB*UStop11 - MTop*Sf2245*Sf6767*UStop12

	Sf6513 = Sf6510*TB*UStop11 - MTop*Sf2245*Sf6512*UStop12

	Sf4677 = Sf4674*TB*UStop11 - 3*MTop*Sf2245*Sf4676*UStop12

	Sf4572 = Sf4569*TB*UStop11 - 3*MTop*Sf2245*Sf4571*UStop12

	Sf6418 = Sf1184*Sf6417 + 
     &   Sf1191*Sf154*Sf2206*Sf2270*Sf4831*Sf4867*Sf6406*TB*
     &    UStau21

	Sf6928 = Sf2263*Sf4753*Sf4867*Sf4906**2*MSq2Diff(2,2)*
     &    MStop2(1) + Sf2270*Sf2334**2*Sf4864*Sf6927*MSq2Diff(2,1)*
     &    MStop2(2)

	Sf6678 = Sf2198*Sf4820**2*Sf4831*Sf6676*MSq2Diff(1,2)*
     &    MStop2(1) + Sf2206*Sf2320**2*Sf4828*Sf6677*MSq2Diff(1,1)*
     &    MStop2(2)

	Sf6935 = Sf4738*TB*MStop2(1) + 
     &   3*Sf2334*Sf2340*Sf6934*TB2*USbot21*UStop11

	Sf6685 = Sf4643*TB*MStop2(1) + 
     &   3*Sf2320*Sf2327*Sf6684*TB2*USbot11*UStop11

	Sf6254 = Sf3938*TB*MStop2(1) + 
     &   Sf2334*Sf2340*Sf6253*TB2*USbot21*UStop11

	Sf5890 = Sf3814*TB*MStop2(1) + 
     &   Sf2320*Sf2327*Sf5889*TB2*USbot11*UStop11

	Sf6252 = Sf6188*MSq2Diff(2,2) - Sf6183*MStop2(2)

	Sf5888 = Sf5810*MSq2Diff(1,2) - Sf5805*MStop2(2)

	Sf6365 = Sf1994*Ss25**3*MSneu2(1) + 
     &   6*Sf1830*Ss241*Ss242*MStau2(2)

	Sf6367 = Sf1936*Ss25**3*MSneu2(1) + 
     &   6*Sf1038*Ss241*Ss242*MStau2(1)

	Sf7159 = Sf2173*MSneu2(1) + Sf1638*MStau2(2)

	Sf2185 = 2*Sf2179*Ss398*MSl2Diff(2,1) + Sf2180*Ss949*MStau2(2)

	Sf2183 = Sf2173*MSl2Diff(2,1) - Sf1662*MStau2(2)

	Sf2010 = Sf166**3*Sf2009*MSl2Diff(1,1) + 
     &   Sf1195**3*Sf2008*MSl2Diff(2,1)

	Sf6350 = Sf5805*MSbot2(1) + Sf6280*MSq2Diff(1,2)

	Sf4216 = Sf3814*MSbot2(1) + Sf3895*MSq2Diff(1,1)

	Sf7152 = -2*CW2*Sf7151 + 
     &   MW2**2*Sf139*Sf2334**3*Sf4906**3*TB*TB2*MHiggs2(4)

	Sf7147 = -2*CW2*Sf7146 + 
     &   MW2**2*Sf126*Sf2320**3*Sf4820**3*TB*TB2*MHiggs2(4)

	Sf7117 = Sf6183*TB + 
     &   Sf4906*Sf4911*Sf6196*TB2*USbot21*UStop21

	Sf7114 = Sf5805*TB + 
     &   Sf4820*Sf4905*Sf5818*TB2*USbot11*UStop21

	Sf4743 = Sf3938*TB + 
     &   Sf2334*Sf2340*Sf4003*TB2*USbot21*UStop11

	Sf4740 = Sf3814*TB + 
     &   Sf2320*Sf2327*Sf3996*TB2*USbot11*UStop11

	Sf7290 = CB*SB*Sf7289*Ss25**2 + 9*MZ2**3*Sf7277*SW2

	Sf2017 = 4*Sf2010*Ss241*Ss242*Ss25**3 - 
     &   MW2*Sf1195*Sf166*Sf2016*TB

	Sf7295 = CB*SB*Sf7293*Ss25**3 - 3*MZ2*Sf7294*SW2

	Sf2184 = 2*Sf845 + 3*Sf2107*Ss928

	Sf7112 = Sf6886*Sf7111 - Sf7108*TB*USbot21*UStop21

	Sf7074 = Sf6635*Sf7073 - Sf7070*TB*USbot11*UStop21

	Sf6887 = Sf6885*Sf6886 - Sf6882*TB*USbot21*UStop21

	Sf6772 = Sf4681*Sf6771 - Sf6768*TB*USbot21*UStop11

	Sf6636 = Sf6634*Sf6635 - Sf6631*TB*USbot11*UStop21

	Sf6517 = Sf4576*Sf6516 - Sf6513*TB*USbot11*UStop11

	Sf6053 = 2*Sf6024 + (MW2*Sf4864*Sf4906*Sf6052*TB)/Ss25**3

	Sf5634 = 2*Sf5604 + (MW2*Sf4820*Sf4828*Sf5633*TB)/Ss25**3

	Sf5374 = 9*Sf5366*Sf5372 + Sf5373*USbot21*UStop21

	Sf5266 = 9*Sf5258*Sf5264 + Sf5265*USbot11*UStop21

	Sf7157 = ((2*CW2*Sf7156)/(Sf4820**3*Sf4906**3*TB*TB2) + 
     &     MW2**2*Sf4899*MHiggs2(4))/CW2

	Sf4682 = Sf4680*Sf4681 - Sf4677*TB*USbot21*UStop11

	Sf4577 = Sf4575*Sf4576 - Sf4572*TB*USbot11*UStop11

	Sf3721 = 2*Sf3678 + (MW2*Sf2263*Sf2334*Sf3720*TB)/Ss25**3

	Sf3449 = 2*Sf3398 + (MW2*Sf2198*Sf2320*Sf3448*TB)/Ss25**3

	Sf2966 = 9*Sf2955*Sf2964 + Sf2965*USbot21*UStop11

	Sf2814 = 9*Sf2802*Sf2812 + Sf2813*USbot11*UStop11

	Sf7142 = ((2*CW2*Sf7141)/(Sf2320**3*Sf2334**3*TB*TB2) + 
     &     MW2**2*Sf2319*MHiggs2(4))/CW2

	Sf7135 = 2*CW2*Sf1038*Sf162 + MW2*Sf153*Sf166**2*MHiggs2(4)

	Sf7136 = (2*CW2*Sf1007 + (3*MW2*Sf166*Sf7135*TB)/Ss25**3)/CW2

	Sf7296 = Sf7295 + 
     &   MW2**2*MZ2*MHiggs2(4)*
     &    (18*Sf5231*SW2*TB + 
     &      SB*(Sf137*SW2*TB*(18 - 18*Log(MSbot2(1)))*MSbot2(1) + 
     &         Sf148*SW2*TB*(18 - 18*Log(MSbot2(2)))*MSbot2(2) + 
     &         Sf949*SW2*TB*(18 - 18*Log(MStau2(1)))*MStau2(1) + 
     &         Sf1737*SW2*TB*(18 - 18*Log(MStau2(2)))*MStau2(2)) + 
     &      CB*(MZ2*SB*Sf5*SW2*TB*(18 - 18*Log(MSneu2(1)))*
     &          MSneu2(1) + 
     &         Sf4265*SW2*TB*(-18 + 18*Log(MStop2(1)))*MStop2(1) + 
     &         Sf6982*SW2*TB*(-18 + 18*Log(MStop2(2)))*MStop2(2)))

	Sf6391 = -2*Sf166**3*Sf6390*Ss241*Ss242 + 
     &   2*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf6366*Ss241*
     &    Ss242*Ss25**3*TB*TB2 + 
     &   MW2*Sf162*Sf166*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*
     &    (-(Sf6367*TB2**2) + 
     &      Sf166*Sf1956*Ss241**2*Ss242**2*TB2**2*MSneu2(1))

	Sf6255 = -(Sf6254*Ss25**3) + 
     &   Sf6106*Ss25**3*TB*MSq2Diff(2,1) + 
     &   MW2*Sf2263*Sf2334*Sf6142*TB2*MSq2Diff(2,1)

	Sf5891 = -(Sf5890*Ss25**3) + 
     &   Sf5717*Ss25**3*TB*MSq2Diff(1,1) + 
     &   MW2*Sf2198*Sf2320*Sf5759*TB2*MSq2Diff(1,1)

	Sf1889 = Sf1195*Sf1872*Ss241*Ss242*Ss25*MSneu2(1) - 
     &   Sf1843*Ss25**3*MSneu2(1) - 
     &   12*Sf1830*Ss241*Ss242*MSl2Diff(2,1)*MStau2(2)

	Sf1148 = Sf1121*Sf166*Ss241*Ss242*Ss25*MSneu2(1) - 
     &   Sf1064*Ss25**3*MSneu2(1) - 
     &   12*Sf1038*Ss241*Ss242*MSl2Diff(1,1)*MStau2(1)

	Sf7129 = -3*Sf4906*Sf7128*TB2 + 
     &   Ss25**3*TB*(Sf6925*MSbot2(2) - Sf7112*MSq2Diff(2,2))

	Sf6355 = Sf4906*Sf6354*TB2 + 
     &   Ss25**3*TB*(Sf6183*MSbot2(2) + Sf6294*MSq2Diff(2,2))

	Sf4758 = -3*Sf2334*Sf4757*TB2 + 
     &   Ss25**3*TB*(Sf4738*MSbot2(2) - Sf4682*MSq2Diff(2,1))

	Sf4221 = Sf2334*Sf4220*TB2 + 
     &   Ss25**3*TB*(Sf3938*MSbot2(2) + Sf3992*MSq2Diff(2,1))

	Sf6937 = Sf6935*Ss25**3 - 
     &   Sf6772*Ss25**3*TB*MSq2Diff(2,1) - 
     &   3*MW2*Sf2334*Sf6936*TB2*MSq2Diff(2,1)

	Sf6687 = Sf6685*Ss25**3 - 
     &   Sf6517*Ss25**3*TB*MSq2Diff(1,1) - 
     &   3*MW2*Sf2320*Sf6686*TB2*MSq2Diff(1,1)

	end

