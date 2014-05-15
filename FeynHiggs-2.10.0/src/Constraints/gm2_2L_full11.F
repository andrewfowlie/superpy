* gm2_2L_full11.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full11
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp122

	Sf6012 = MW2*(18*MW2**2*Sf5988*Sf5989*Sf96 + 
     &      (9*Sf4864*Sf4906*Sf6010*TB)/Ss25**3) + 
     &   Sf6011*USbot21*UStop21

	Sf5957 = MW2*(6*MW2**2*Sf5443*Sf5932*Sf96 - 
     &      (3*Sf4864*Sf4906*Sf5955*TB)/Ss25**3) - 
     &   Sf5956*USbot21*UStop21

	Sf5592 = MW2*(18*MW2**2*Sf5562*Sf5563*Sf88 + 
     &      (9*Sf4820*Sf4828*Sf5590*TB)/Ss25**3) + 
     &   Sf5591*USbot11*UStop21

	Sf5531 = MW2*(6*MW2**2*Sf5337*Sf5499*Sf88 - 
     &      (3*Sf4820*Sf4828*Sf5529*TB)/Ss25**3) - 
     &   Sf5530*USbot11*UStop21

	Sf3614 = MW2*(18*MW2**2*Sf3577*Sf3579*Sf74 + 
     &      (9*Sf2263*Sf2334*Sf3612*TB)/Ss25**3) + 
     &   Sf3613*USbot21*UStop11

	Sf3532 = MW2*(6*MW2**2*Sf3053*Sf3489*Sf74 - 
     &      (3*Sf2263*Sf2334*Sf3530*TB)/Ss25**3) - 
     &   Sf3531*USbot21*UStop11

	Sf3287 = MW2*(18*MW2**2*Sf3225*Sf3227*Sf66 + 
     &      (9*Sf2198*Sf2320*Sf3285*TB)/Ss25**3) + 
     &   Sf3286*USbot11*UStop11

	Sf3170 = MW2*(6*MW2**2*Sf2926*Sf3106*Sf66 - 
     &      (3*Sf2198*Sf2320*Sf3168*TB)/Ss25**3) - 
     &   Sf3169*USbot11*UStop11

	Sf6419 = 2*CW2*Sf6418*MHiggs2(4) + 
     &   Sf1184*Sf154*(2*Sf6405 + 
     &      Sf2206*Sf2270*Sf4831*Sf4867*Sf6404*MHiggs2(4)**2)

	Sf7137 = (2*Sf1818)/Sf1195**3 + 
     &   (TB*((6*CW2*MW2*Sf1191*Sf1830)/Sf1195**2 + 
     &        3*MW2**2*Sf1183*MHiggs2(4)))/(CW2*Ss25**3)

	Sf2186 = 18*MW2*Sf162*Sf166*Sf2059*TB*MSl2Diff(1,1)*
     &    MSl2Diff(2,1) + Ss25**3*MSl2Diff(1,1)*
     &    (-3*Sf166**3*Sf2185 + 2*Sf2104*MSl2Diff(2,1)) - 
     &   Sf2184*Ss25**3*MSl2Diff(2,1)*MStau2(1)

	Sf6356 = Sf4906**3*Sf6350*Ss25**3*TB*MSq2Diff(2,2) + 
     &   Sf4820*(Sf4820**2*Sf6355*MSq2Diff(1,2) + 
     &      Sf4906**3*Sf6352*TB2*MSq2Diff(2,2))

	Sf6258 = Sf2334**3*
     &    (Sf6252*Ss25**3*TB*MSq2Diff(2,1) + 
     &      Sf4906*Sf6257*TB2*MSq2Diff(2,1)) + 
     &   Sf4906**3*Sf6255*MSq2Diff(2,2)

	Sf4222 = Sf2334**3*Sf4216*Ss25**3*TB*MSq2Diff(2,1) + 
     &   Sf2320*(Sf2320**2*Sf4221*MSq2Diff(1,1) + 
     &      Sf2334**3*Sf4218*TB2*MSq2Diff(2,1))

	Sf5894 = Sf2320**3*
     &    (Sf5888*Ss25**3*TB*MSq2Diff(1,1) + 
     &      Sf4820*Sf5893*TB2*MSq2Diff(1,1)) + 
     &   Sf4820**3*Sf5891*MSq2Diff(1,2)

	Sf6266 = -(MW2*Sf2334*Sf4906*Sf6265*TB2) + 
     &   2*Sf2270*Sf4867*Ss241**2*Ss242**2*MHiggs2(4)*
     &    (Sf6258 + MW2*Sf2334**3*Sf4906**3*Sf6251*Ss25*TB*TB2*
     &       MHiggs2(4)*MSq2Diff(2,1)*MSq2Diff(2,2)*USbot21)

	Sf5902 = -(MW2*Sf2320*Sf4820*Sf5901*TB2) + 
     &   2*Sf2206*Sf4831*Ss241**2*Ss242**2*MHiggs2(4)*
     &    (Sf5894 + MW2*Sf2320**3*Sf4820**3*Sf5887*Ss25*TB*TB2*
     &       MHiggs2(4)*MSq2Diff(1,1)*MSq2Diff(1,2)*USbot11)

	Sf6688 = Sf4831*Sf6687 + 
     &   3*MW2**2*Sf2320**3*Sf4828*Sf6683*TB*TB2*MSq2Diff(1,1)*
     &    USbot11*UStop21

	Sf6931 = Sf6887*MSq2Diff(2,2) - Sf6925*MStop2(2)

	Sf6680 = Sf6636*MSq2Diff(1,2) - Sf6674*MStop2(2)

	Sf7168 = Sf6887*MSbot2(2) + Sf7112*MStop2(2)

	Sf7163 = Sf6772*MSbot2(2) + Sf4682*MStop2(1)

	Sf7120 = Sf6674*MSbot2(1) - Sf7074*MSq2Diff(1,2)

	Sf7174 = Sf6636*MSbot2(1) + Sf7074*MStop2(2)

	Sf4750 = Sf4643*MSbot2(1) - Sf4577*MSq2Diff(1,1)

	Sf7118 = MW2*Sf4864*Sf4906*Sf7116*TB2 + 
     &   2*Sf4867*Sf7117*Ss241**2*Ss242**2*Ss25**2*MHiggs2(4)

	Sf7115 = MW2*Sf4820*Sf4828*Sf7113*TB2 + 
     &   2*Sf4831*Sf7114*Ss241**2*Ss242**2*Ss25**2*MHiggs2(4)

	Sf4744 = MW2*Sf2263*Sf2334*Sf4742*TB2 + 
     &   2*Sf2270*Sf4743*Ss241**2*Ss242**2*Ss25**2*MHiggs2(4)

	Sf4741 = MW2*Sf2198*Sf2320*Sf4739*TB2 + 
     &   2*Sf2206*Sf4740*Ss241**2*Ss242**2*Ss25**2*MHiggs2(4)

	Sf2018 = MW2*Sf1195**3*Sf166**3*Sf2007*Ss241*Ss242*TB + 
     &   Sf1184*Sf154*Sf2017*MHiggs2(4)

	Sf1890 = MW2*Sf1191*Sf1195*Sf1889*TB - 
     &   4*Sf1888*Ss241*Ss242*Ss25**3*MStau2(2)

	Sf1149 = MW2*Sf1148*Sf162*Sf166*TB - 
     &   4*Sf1147*Ss241*Ss242*Ss25**3*MStau2(1)

	Sf7161 = Sf813 + 3*Sf7136*Ss928

	Sf7130 = Sf4867*Sf7129 - 
     &   3*MW2**2*Sf4864*Sf4906**3*Sf6412*TB*TB2*MSq2Diff(2,2)*
     &    USbot21*UStop21

	Sf4759 = Sf2270*Sf4758 - 
     &   3*MW2**2*Sf2263*Sf2334**3*Sf4754*TB*TB2*MSq2Diff(2,1)*
     &    USbot21*UStop11

	Sf2187 = Sf154*Sf2186 + 
     &   18*MW2**2*Sf162*Sf165*Sf166**3*TB*MSl2Diff(1,1)*
     &    MSl2Diff(2,1)*UStau11

	Sf1892 = Sf1184*Sf1890*MHiggs2(4) + 
     &   MW2*Sf1195**3*Sf1891*Ss241*Ss242*TB*MStau2(2)

	Sf1151 = Sf1149*Sf154*MHiggs2(4) + 
     &   MW2*Sf1150*Sf166**3*Ss241*Ss242*TB*MStau2(1)

	Sf7182 = 3*MW2*Sf2320*Sf7181*TB2 + 
     &   TB*(Sf6517*Ss25**3*MSbot2(1) + Sf4577*Ss25**3*MStop2(1) + 
     &      Sf2320**3*TB2*(Sf7137*Ss25**3*Ss949*MStau2(2) - 
     &         3*Sf7142*Ss969*MStop2(1) - 3*Sf7157*Ss988*MStop2(2))
     &      ) + 3*Sf2320*Sf2327*Ss25**3*TB2*
     &    (-(Sf6422*MSbot2(1)) + Sf4273*MStop2(1))*USbot11*
     &    UStop11

	Sf6392 = Sf1195**3*Sf6391 + 
     &   2*Sf166**3*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf6364*
     &    Ss241*Ss242*Ss25**3*TB*TB2 + 
     &   MW2*Sf1191*Sf1195*Sf166**3*Sf2320**3*Sf2334**3*Sf4820**3*
     &    Sf4906**3*(-(Sf6365*TB2**2) + 
     &      Sf1195*Sf2005*Ss241**2*Ss242**2*TB2**2*MSneu2(1))

	Sf2188 = Sf1195**3*Sf2187 + 
     &   Sf154*Sf166**3*MSl2Diff(1,1)*
     &    (2*Sf2183*Ss25**3 + 
     &      18*MW2*Sf1191*Sf1195*Sf2139*TB*MSl2Diff(2,1))

	Sf7131 = Sf7130*Ss241**2*Ss242**2*MHiggs2(4) + 
     &   3*Sf7118*Ss25*Ss906*MSbot2(2) - 
     &   3*MW2*Sf4864*Sf4906*Sf6927*Ss25**3*TB2*MSbot2(2)

	Sf4760 = Sf4759*Ss241**2*Ss242**2*MHiggs2(4) + 
     &   3*Sf4744*Ss25*Ss906*MSbot2(2) - 
     &   3*MW2*Sf2263*Sf2334*Sf4753*Ss25**3*TB2*MSbot2(2)

	Sf6362 = -(MW2*Sf4820*Sf4906*Sf6361*TB2) + 
     &   2*Sf4831*Sf4867*Ss241**2*Ss242**2*MHiggs2(4)*
     &    (Sf6356 + MW2*Sf4820**3*Sf4906**3*Sf6349*Ss25*TB*TB2*
     &       MHiggs2(4)*MSq2Diff(1,2)*MSq2Diff(2,2)*UStop21)

	Sf4223 = MW2*Sf2320*Sf2334*Sf4215*TB2 + 
     &   2*Sf2206*Sf2270*Ss241**2*Ss242**2*MHiggs2(4)*
     &    (Sf4222 + MW2*Sf2320**3*Sf2334**3*Sf4208*Ss25*TB*TB2*
     &       MHiggs2(4)*MSq2Diff(1,1)*MSq2Diff(2,1)*UStop11)

	Sf4761 = Sf2320**3*Sf4760*MSq2Diff(1,1) + 
     &   Sf2270*MHiggs2(4)*
     &    (Sf2334**3*Sf4750*Ss241**2*Ss242**2*Ss25**3*TB*
     &       MSq2Diff(2,1) - 
     &      3*Sf2320*Sf2334**3*Sf4752*Ss241**2*Ss242**2*TB2*
     &       MSq2Diff(2,1))

	Sf6938 = Sf2334**3*
     &    (Sf6931*Ss25**3*TB*MSq2Diff(2,1) + 
     &      3*Sf4906*Sf6933*TB2*MSq2Diff(2,1)) - 
     &   Sf4906**3*Sf6937*MSq2Diff(2,2)

	Sf6689 = Sf2320**3*
     &    (Sf4831*Sf6680*Ss25**3*TB*MSq2Diff(1,1) + 
     &      3*Sf4820*Sf4831*Sf6682*TB2*MSq2Diff(1,1)) - 
     &   Sf4820**3*Sf6688*MSq2Diff(1,2)

	Sf7183 = Sf4820**3*Sf7182 + 
     &   Sf2320**3*(Sf7174*Ss25**3*TB + 3*Sf4820*Sf7177*TB2)

	Sf6939 = -(Sf4867*Sf6938) + 
     &   3*MW2**2*Sf2334**3*Sf4864*Sf4906**3*Sf6930*TB*TB2*
     &    MSq2Diff(2,1)*MSq2Diff(2,2)*USbot21*UStop21

	Sf7132 = Sf4831*Sf7131*MSq2Diff(1,2) + 
     &   3*MW2*Sf4828*Sf4867*Sf4906**3*Sf7125*Ss241*Ss242*TB*TB2*
     &    MSq2Diff(2,2)*USbot11*UStop21

	Sf6940 = Sf2270*Sf6939 + 
     &   3*MW2**2*Sf2263*Sf2334**3*Sf4867*Sf4906**3*Sf6929*TB*TB2*
     &    MSq2Diff(2,1)*MSq2Diff(2,2)*USbot21*UStop11

	Sf6690 = -(Sf2206*Sf6689) + 
     &   3*MW2**2*Sf2198*Sf2320**3*Sf4820**3*Sf4831*Sf6679*TB*TB2*
     &    MSq2Diff(1,1)*MSq2Diff(1,2)*USbot11*UStop11

	Sf7162 = Sf2104*MSneu2(1) + Sf7161*MStau2(1)

	Sf7164 = -(CW2*Sf4906**3*Sf7163*Ss25**3*TB) + 
     &   3*Sf7152*Ss906*MSbot2(2)

	Sf4749 = -(MW2*Sf2198*Sf2320*Sf4748*TB2) + 
     &   Sf4741*Ss25*Ss885*MSbot2(1)

	Sf7121 = Sf4831*Sf7120*Ss241**2*Ss242**2*Ss25**2*TB*
     &    MHiggs2(4) + 3*Sf7115*Ss885*MSbot2(1)

	Sf6393 = CW2*Sf6392 + 
     &   3*MW2**2*Sf1195**3*Sf166**3*Sf2320**3*Sf2334**3*Sf4820**3*
     &    Sf4906**3*Sf6363*Ss241*Ss242*TB2*MHiggs2(4)

	Sf6941 = -3*MW2*Sf2334*Sf4906*Sf6928*Ss25**3*TB2 + 
     &   Sf6940*Ss241**2*Ss242**2*MHiggs2(4)

	Sf6691 = -3*MW2*Sf2320*Sf4820*Sf6678*Ss25**3*TB2 + 
     &   Sf6690*Ss241**2*Ss242**2*MHiggs2(4)

	Sf6942 = 2*CW2*Sf6941 - 
     &   Sf2270*Sf2334**3*Sf4867*Sf4906**3*Sf6926*Ss241**2*
     &    Ss242**2*TB*TB2*MHiggs2(4)*MSq2Diff(2,1)*MSq2Diff(2,2)

	Sf6692 = 2*CW2*Sf6691 - 
     &   Sf2206*Sf2320**3*Sf4820**3*Sf4831*Sf6675*Ss241**2*
     &    Ss242**2*TB*TB2*MHiggs2(4)*MSq2Diff(1,1)*MSq2Diff(1,2)

	Sf2189 = Sf1184*Sf2188 + 
     &   18*MW2**2*Sf1191*Sf1194*Sf1195**3*Sf154*Sf166**3*TB*
     &    MSl2Diff(1,1)*MSl2Diff(2,1)*UStau21

	Sf7184 = -2*CW2*Sf7183 - 
     &   Sf2320**3*Sf4820**3*Sf7173*TB2 + 6*Sf7147*Ss885*MSbot2(1)

	Sf7133 = Sf4867*Sf4906**3*Sf7121*Ss25*MSq2Diff(2,2) + 
     &   Sf4820*(Sf4820**2*Sf7132 - 
     &      3*Sf4867*Sf4906**3*Sf7124*TB2*MSq2Diff(2,2))

	Sf2190 = 18*MW2*Sf1195*Sf166*Sf2182*Ss25*TB + 
     &   MHiggs2(4)*(Sf2189 + 
     &      18*MW2**2*Sf1184*Sf1195**3*Sf154*Sf166**3*Sf2176*Ss275*
     &       TB*MHiggs2(4)*MSl2Diff(1,1)*MSl2Diff(2,1))

	Sf7185 = Sf4906**3*Sf7184 + 
     &   CW2*Sf2320**3*Sf4820**3*
     &    (-2*Sf7168*Ss25**3*TB - 6*Sf4906*Sf7171*TB2)

	Sf4762 = Sf2206*Sf4761 + 3*Sf2270*Sf2334**3*Sf4749*MSq2Diff(2,1)

	Sf2191 = -(CW2*Sf2190) + 
     &   3*Sf1184*Sf1195**3*Sf154*Sf166**3*Sf2181*MHiggs2(4)*
     &    MSl2Diff(1,1)*MSl2Diff(2,1)

	Sf7134 = 2*CW2*Sf7133 - 
     &   Sf4820**3*Sf4831*Sf4867*Sf4906**3*Sf7119*Ss241**2*
     &    Ss242**2*TB*TB2*MHiggs2(4)*MSq2Diff(1,2)*MSq2Diff(2,2)

	Sf4763 = 2*CW2*Sf4762 - 
     &   Sf2206*Sf2270*Sf2320**3*Sf2334**3*Sf4745*Ss241**2*
     &    Ss242**2*TB*TB2*MHiggs2(4)*MSq2Diff(1,1)*MSq2Diff(2,1)

	Sf7186 = Sf2334**3*Sf7185 + 
     &   Sf2320**3*Sf4820**3*
     &    (2*Sf7164 - 6*CW2*Sf2334*Sf4906**3*Sf7167*TB2)

	Sf7187 = 3*Sf166**3*Sf7186 + 
     &   CW2*(-2*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*Sf7162*
     &       Ss25**3*TB*TB2 - 
     &      18*MW2*Sf162*Sf166*Sf2320**3*Sf2334**3*Sf4820**3*
     &       Sf4906**3*Sf7160*TB2**2)

	Sf7188 = Sf1195**3*Sf7187 + 
     &   CW2*(-2*Sf166**3*Sf2320**3*Sf2334**3*Sf4820**3*Sf4906**3*
     &       Sf7159*Ss25**3*TB*TB2 - 
     &      18*MW2*Sf1191*Sf1195*Sf166**3*Sf2320**3*Sf2334**3*
     &       Sf4820**3*Sf4906**3*Sf7158*TB2**2)

	gm2_2L = gm2_2L + 1/480.D0*(12*Sf2572*Sf3170 + Sf2561*Sf3287 + 
     &      36*Sf3058*Sf3449 + 
     &      (Sf2320**3*(12*Sf2749*Sf3532 + Sf2739*Sf3614 + 
     &           36*Sf3066*Sf3721))/Sf2334**3)/
     &    (MW2**3*Sf2320**3*SW2*TB2)

	gm2_2L = gm2_2L + ((1/40.D0*(Sf2583*(3*Sf2588 - Sf3458)) + 
     &        1/480.D0*(Sf2551*(-12*Sf2560 + 36*Sf4231*Sf4233) - 
     &            9*Sf2220*Sf4231*Sf69 + Sf4246*Sf69)/
     &          (TB*MHiggs2(4)))*USbot11*UStop11)/
     &    (MW2**3*SW2*MSq2Diff(1,1))

	gm2_2L = gm2_2L + 1/20.D0*
     &    (-(Sf626*Sf7571*UStau11**2) - 
     &       Sf1521*Sf7629*UStau21**2 + 
     &       UStop11*(-3*
     &           (Sf3058*Sf3063*USbot11 + 
     &             Sf3066*Sf3071*USbot21) - 
     &          (Sf2583*(Sf7852*USbot11**2 + 
     &                Sf7854*USbot11**2) + 
     &             Sf2759*(Sf7930*USbot21**2 + 
     &                Sf7932*USbot21**2))*UStop11) - 
     &       3*(Sf5447*Sf5452*USbot11 + 
     &          Sf5454*Sf5459*USbot21)*UStop21 - 
     &       (Sf5016*(Sf8061*USbot11**2 + 
     &             Sf8063*USbot11**2) + 
     &          Sf5166*(Sf8139*USbot21**2 + 
     &             Sf8141*USbot21**2))*UStop21**2)/
     &     (MW2**3*SW2)

	gm2_2L = gm2_2L + 1/1440.D0*
     &    (((12*Sf5037*Sf8016 - Sf5028*Sf8028 - 72*Sf5447*Sf8068)/
     &          Sf4820**2 + 
     &         (Sf5028*Sf8036 - 12*Sf5037*Sf8051 - 
     &            72*Sf5447*Sf8075)/Sf4820**3)*USbot11**2*
     &       UStop21**2)/(MW2**3*SW2)

	gm2_2L = gm2_2L + ((1/40.D0*(Sf5166*(3*Sf5170 - Sf6060)) + 
     &        1/480.D0*(Sf5139*(-12*Sf5146 + 36*Sf6953*Sf6955) - 
     &            9*Sf4880*Sf6953*Sf99 + Sf6974*Sf99)/
     &          (TB*MHiggs2(4)))*USbot21*UStop21)/
     &    (MW2**3*SW2*MSq2Diff(2,2))

	gm2_2L = gm2_2L + ((1/80.D0*(Sf4763*Ss969) + 
     &        1/160.D0*(CW2*Sf4223*Ss971))*MStop2(1))/
     &    (CW2*MW2**2*Sf2206*Sf2270*Sf2320**3*Sf2334**3*Ss241**2*
     &      Ss242**2*Ss25**3*SW2*TB*TB2*MHiggs2(4)*MSq2Diff(1,1)*
     &      MSq2Diff(2,1))

	gm2_2L = gm2_2L + 1/1440.D0*
     &    (((Sf2739*Sf7905 - 12*Sf2749*Sf7920 - 72*Sf3066*Sf7944)/
     &          Sf2334**3 + 
     &         (Sf2730*(-12*Sf7868 - 12*Sf7870) + 
     &            Sf77*(Sf7873 + Sf7875))/MSq2Diff(2,1)**2)*
     &       USbot21**2*UStop11**2)/(MW2**3*SW2)

	gm2_2L = gm2_2L + ((3/160.D0*(Sf2739*Sf2748)/(Ss25**2*TB) + 
     &        3/40.D0*(Sf2749*Sf2755)/(Ss25**2*TB) + 
     &        Sf2334**2*Sf2759*
     &         (3/40.D0*Sf2762/MSq2Diff(2,1) - 
     &           1/40.D0*Sf3729/MSq2Diff(2,1)))*USbot21*UStop11
     &      )/(MW2**3*Sf2334**2*SW2)

	gm2_2L = gm2_2L + ((-(1/720.D0*Sf7850/MW2**2) + 
     &        1/1440.D0*(Sf2561*Sf7827 - 12*Sf2572*Sf7842 - 
     &            72*Sf3058*Sf7866)/MW2**3)*USbot11**2*
     &      UStop11**2)/(Sf2320**3*SW2)

	gm2_2L = gm2_2L + ((1/80.D0*(Sf7134*Ss988) + 
     &        1/160.D0*(CW2*Sf6362*Ss990))*MStop2(2))/
     &    (CW2*MW2**2*Sf4820**3*Sf4831*Sf4867*Sf4906**3*Ss241**2*
     &      Ss242**2*Ss25**3*SW2*TB*TB2*MHiggs2(4)*MSq2Diff(1,2)*
     &      MSq2Diff(2,2))

	gm2_2L = gm2_2L + 1/360.D0*
     &    ((MW2*(-2*Sf7502*
     &             (Sf7503*MSbot2(1) + Sf7522*MSbot2(2) + 
     &               3*Sf7580*MStau2(1) + 3*Sf7714*MStau2(2) + 
     &               Sf7991*MStop2(1) + Sf8533*MStop2(2)) + 
     &            Sf7497*(7*Sf7501**2*MSbot2(1) + 
     &               7*Sf7521**2*MSbot2(2) - 45*MSneu2(1) + 
     &               7*(3*Sf7579**2*MStau2(1) + 
     &                  3*Sf7713**2*MStau2(2) + 
     &                  Sf7990**2*MStop2(1) + Sf8532**2*MStop2(2)))
     &            ))/MZ2**2 + 
     &       (Sf8677*Ss969)/
     &        (MW2**2*MZ2**3*Sf2320**3*Sf2334**3*MSq2Diff(1,1)**2*
     &          MSq2Diff(2,1)**2*MStop2(3)**2))/(CW2**2*SW2)

	gm2_2L = gm2_2L + 1/1440.D0*
     &    ((Sf1192*Sf7598*UStau21**2 + 
     &          Sf1386*(12*Sf1457 + 12*Sf7606*UStau21**2) + 
     &          Sf1461*(-72*Sf1517 - 72*Sf7636*UStau21**2))/
     &        Sf1195**3 + ((12*Sf2572*Sf7813 - Sf2561*Sf7819 - 
     &            72*Sf3058*Sf7859)*USbot11**2*UStop11**2)/
     &        Sf2320**2)/(MW2**3*SW2)

	gm2_2L = gm2_2L + ((-(1/720.D0*(MW2*Sf7928)/Sf2334**3) + 
     &        (1/120.D0*(Sf2749*Sf7891) - 1/1440.D0*(Sf2739*Sf7897) - 
     &           1/20.D0*(Sf3066*Sf7937))/Sf2334**2)*USbot21**2*
     &      UStop11**2)/(MW2**3*SW2)

	gm2_2L = gm2_2L + (-(1/160.D0*
     &         (CW2*Sf5902*Ss887*MSbot2(1))/
     &          (Sf2206*Sf2320**3*Sf4820**3*Sf4831*MSq2Diff(1,1)*
     &            MSq2Diff(1,2))) + 
     &      1/80.D0*(Sf6942*Ss906*MSbot2(2))/
     &        (Sf2270*Sf2334**3*Sf4867*Sf4906**3*MSq2Diff(2,1)*
     &          MSq2Diff(2,2)))/
     &    (CW2*MW2**2*Ss241**2*Ss242**2*Ss25**3*SW2*TB*TB2*
     &      MHiggs2(4))

	gm2_2L = gm2_2L + 1/1440.D0*
     &    (((Sf5147*Sf8114 - 12*Sf5156*Sf8129 - 72*Sf5454*Sf8153)/
     &          Sf4906**3 + 
     &         (Sf5139*(-12*Sf8077 - 12*Sf8079) + Sf8082*Sf99 + 
     &            Sf8084*Sf99)/MSq2Diff(2,2)**2)*USbot21**2*
     &       UStop21**2)/(MW2**3*SW2)

	gm2_2L = gm2_2L + 1/480.D0*
     &    (12*Sf5037*Sf5531 + Sf5028*Sf5592 + 36*Sf5447*Sf5634 + 
     &       (Sf4820**3*(12*Sf5156*Sf5957 + Sf5147*Sf6012 + 
     &            36*Sf5454*Sf6053))/Sf4906**3)/
     &     (MW2**3*Sf4820**3*SW2*TB2)

	gm2_2L = gm2_2L + (1/360.D0*
     &       (Sf8714*Ss988)/
     &        (MZ2**3*MSq2Diff(1,2)**2*MSq2Diff(2,2)**2*
     &          MStop2(3)**2) + 
     &      1/240.D0*(CW2**2*Sf8198*Ss990*MStop2(2)*UStop21**2))/
     &    (CW2**2*MW2**2*Sf4820**3*Sf4906**3*SW2)

	gm2_2L = gm2_2L + ((1/240.D0*(Sf2022*Sf2191) + 
     &        1/960.D0*(CW2*Sf1894*Sf2018)/(Ss241*Ss242))*MSneu2(1))/
     &    (CW2*MW2**2*Sf1184*Sf1195**3*Sf154*Sf166**3*Ss25**3*SW2*
     &      MHiggs2(4)*MSl2Diff(1,1)*MSl2Diff(2,1))

	gm2_2L = gm2_2L + 1/1440.D0*
     &    ((Sf163*Sf7540*UStau11**2 + 
     &          Sf469*(12*Sf554 + 12*Sf7548*UStau11**2) + 
     &          Sf558*(-72*Sf622 - 72*Sf7578*UStau11**2))/
     &        Sf166**3 + ((-12*Sf1378*Sf1385 + Sf1377*Sf47)*
     &          UStau21)/(MHiggs2(4)*MSl2Diff(2,1)))/
     &     (MW2**3*SW2)

	gm2_2L = gm2_2L + ((1/240.D0*
     &         (Sf7989*Ss971*MStop2(1))/Sf2334**3 - 
     &        1/720.D0*(Sf2320*Sf7802*USbot11**2))*UStop11**2)/
     &    (MW2**2*Sf2320**3*SW2)

	gm2_2L = gm2_2L + ((1/40.D0*(Sf5016*(3*Sf5019 - Sf5468)) + 
     &        1/480.D0*(Sf5020*(-12*Sf5027 + 36*Sf6949*Sf6951) - 
     &            9*Sf4844*Sf6949*Sf91 + Sf6965*Sf91)/
     &          (TB*MHiggs2(4)))*USbot11*UStop21)/
     &    (MW2**3*SW2*MSq2Diff(1,2))

	gm2_2L = gm2_2L + 1/240.D0*
     &    ((CA*Sf10435*Sf8520*Ss78*MStop2(4))/
     &        (CB*CW**3*MZ2*SB*MHiggs2(2)) + 
     &       (Sf8256*Ss887*MSbot2(1)*USbot11**2)/
     &        (MW2**2*Sf2320**3*Sf4820**3))/SW2

	gm2_2L = gm2_2L + ((1/1440.D0*
     &         (Sf5020*(-12*Sf7998 - 12*Sf8000) + Sf8003*Sf91 + 
     &            Sf8005*Sf91)/MW2**3 + 
     &        1/720.D0*(Sf8010*MSbot2(1)*MStop2(2))/
     &          (MW2**2*Sf4820**2))*USbot11**2*UStop21**2)/
     &    (SW2*MSq2Diff(1,2)**2)

	gm2_2L = gm2_2L + (1/360.D0*
     &       (Sf8414*Ss906)/
     &        (MZ2**3*MSbot2(3)**2*MSq2Diff(2,1)**2*
     &          MSq2Diff(2,2)**2) - 
     &      1/160.D0*(CW2**2*Sf6266*Ss908*MSbot2(2))/
     &        (Sf2270*Sf4867*Ss241**2*Ss242**2*Ss25**3*TB*TB2*
     &          MHiggs2(4)*MSq2Diff(2,1)*MSq2Diff(2,2)))/
     &    (CW2**2*MW2**2*Sf2334**3*Sf4906**3*SW2)

	gm2_2L = gm2_2L + 1/960.D0*
     &    ((Sf1151*Ss930)/(Sf154*Sf166**3*MSl2Diff(1,1)) + 
     &       (Sf1892*Ss951)/(Sf1184*Sf1195**3*MSl2Diff(2,1)))/
     &     (MW2**2*Ss241*Ss242*Ss25**3*SW2*MHiggs2(4))

	gm2_2L = gm2_2L + ((1/240.D0*
     &         (Sf8227*Ss908*USbot21**2)/Sf4906**3 - 
     &        3/80.D0*(Sf2334*Sf2338*Sf2340*USbot21*UStop11)/TB
     &        )*MSbot2(2))/(MW2**2*Sf2334**3*SW2)

	gm2_2L = gm2_2L + (Ss885*
     &      (1/360.D0*Sf8318/
     &          (MZ2**3*MSbot2(3)**2*MSq2Diff(1,1)**2*
     &            MSq2Diff(1,2)**2) + 
     &        1/80.D0*(CW2*Sf6692*MSbot2(1))/
     &          (Sf2206*Sf4831*Ss241**2*Ss242**2*Ss25**3*TB*TB2*
     &            MHiggs2(4)*MSq2Diff(1,1)*MSq2Diff(1,2))))/
     &    (CW2**2*MW2**2*Sf2320**3*Sf4820**3*SW2)

	gm2_2L = gm2_2L + ((-(1/40.D0*
     &           (Sf1521*Sf1528)/(MW2**3*MSl2Diff(2,1))) - 
     &        3/80.D0*(Sf1191*TB*
     &            (Sf1194*MSneu2(1) - Sf1189*MStau2(2)))/
     &          (Sf1184*Ss25**3))*UStau21)/SW2

	gm2_2L = gm2_2L + ((1/1440.D0*
     &         (Sf2551*(-12*Sf7790 - 12*Sf7792) + 
     &            Sf69*(Sf7795 + Sf7797))/MW2**3 - 
     &        1/720.D0*(Sf7807*MSbot2(1)*MStop2(1))/
     &          (MW2**2*Sf2320**2))*USbot11**2*UStop11**2)/
     &    (SW2*MSq2Diff(1,1)**2)

	gm2_2L = gm2_2L + 1/240.D0*
     &    (CW2*(2*Sf8486*Ss398 + Sf8353*Ss455) + 
     &       (Sf7188*Ss398)/(Ss25**3*TB*TB2))/
     &     (CW2*MW2**2*Sf1195**3*Sf166**3*Sf2320**3*Sf2334**3*
     &       Sf4820**3*Sf4906**3*SW2)

	gm2_2L = gm2_2L + ((3/20.D0*
     &         (Sf2251*Sf2256*Sf2260)/MHiggs2(4)**2 + 
     &        (Sf2198*(3/160.D0*(Sf2215*Sf2221) + 
     &             3/40.D0*(Sf2221*Sf2222)))/MHiggs2(4))*USbot11*
     &      UStop11)/(Sf2206*Ss25**2*SW2)

	gm2_2L = gm2_2L - 3/20.D0*
     &    ((Sf2198*Sf2256*Sf2365)/Sf2206 + 
     &       (Sf2263*Sf2307*Sf2603)/Sf2270)/
     &     (Ss25**3*SW2*TB*MHiggs2(4)**2)

	gm2_2L = gm2_2L - 3/20.D0*
     &    ((Sf4828*Sf4862*Sf4931)/Sf4831 + 
     &       (Sf4864*Sf4896*Sf5057)/Sf4867)/
     &     (Ss25**3*SW2*TB*MHiggs2(4)**2)

	gm2_2L = gm2_2L + ((3/160.D0*(Sf5147*Sf5155)/Ss25**2 + 
     &        3/40.D0*(Sf5156*Sf5162)/Ss25**2 - 
     &        3/80.D0*(MW2*Sf4910*Sf4911*MSbot2(2)))*USbot21*
     &      UStop21)/(MW2**3*Sf4906**2*SW2*TB)

	gm2_2L = gm2_2L + ((1/480.D0*
     &         (Sf2730*(-12*Sf2738 + 36*Sf4235*Sf4237) - 
     &            9*Sf2283*Sf4235*Sf77 + Sf4255*Sf77)/
     &          (MW2**3*MHiggs2(4)) - 
     &        3/80.D0*(Sf2340*Sf2346*MStop2(1))/(MW2**2*Sf2334**2))*
     &      USbot21*UStop11)/(SW2*TB*MSq2Diff(2,1))

	gm2_2L = gm2_2L + ((3/20.D0*
     &         (Sf4891*Sf4896*Sf4897)/MHiggs2(4)**2 + 
     &        (Sf4864*(3/160.D0*(Sf4875*Sf4881) + 
     &             3/40.D0*(Sf4881*Sf4882)))/MHiggs2(4))*USbot21*
     &      UStop21)/(Sf4867*Ss25**2*SW2)

	gm2_2L = gm2_2L + ((3/20.D0*
     &         (Sf2302*Sf2307*Sf2311)/MHiggs2(4)**2 + 
     &        (Sf2263*(3/160.D0*(Sf2278*Sf2284) + 
     &             3/40.D0*(Sf2284*Sf2285)))/MHiggs2(4))*USbot21*
     &      UStop11)/(Sf2270*Ss25**2*SW2)

	gm2_2L = gm2_2L + ((3/160.D0*(Sf2561*Sf2571)/Ss25**2 + 
     &        3/40.D0*(Sf2572*Sf2579)/Ss25**2 - 
     &        3/80.D0*(MW2*Sf2324*Sf2327*MSbot2(1)))*USbot11*
     &      UStop11)/(MW2**3*Sf2320**2*SW2*TB)

	gm2_2L = gm2_2L + ((-(1/40.D0*(Sf626*Sf634)) + 
     &        1/1440.D0*(Sf34*Sf460)/MHiggs2(4) - 
     &        1/120.D0*(Sf461*Sf468)/MHiggs2(4))*UStau11)/
     &    (MW2**3*SW2*MSl2Diff(1,1))

	gm2_2L = gm2_2L + ((-(1/720.D0*(MW2*Sf8137)/Sf4906**3) + 
     &        (1/120.D0*(Sf5156*Sf8100) - 1/1440.D0*(Sf5147*Sf8106) - 
     &           1/20.D0*(Sf5454*Sf8146))/Sf4906**2)*USbot21**2*
     &      UStop21**2)/(MW2**3*SW2)

	gm2_2L = gm2_2L + ((1/240.D0*
     &         (Sf1665*Ss949)/(Sf1184*Ss25**3*MHiggs2(4)) + 
     &        1/720.D0*(CW2*Sf1302*MSneu2(1)))*MStau2(2))/
     &    (CW2*MW2**2*Sf1195**3*SW2*MSl2Diff(2,1))

	gm2_2L = gm2_2L + ((1/240.D0*
     &         (Sf848*Ss928)/(Sf154*Ss25**3*MHiggs2(4)) + 
     &        1/720.D0*(CW2*Sf351*MSneu2(1)))*MStau2(1))/
     &    (CW2*MW2**2*Sf166**3*SW2*MSl2Diff(1,1))

	gm2_2L = gm2_2L + ((-(1/120.D0*(Sf1378*Sf7587)) + 
     &        1/1440.D0*(Sf47*Sf7590) - 
     &        1/720.D0*(MW2*Sf7613*MSneu2(1)*MStau2(2))/Sf1195**3)*
     &      UStau21**2)/(MW2**3*SW2*MSl2Diff(2,1)**2)

	gm2_2L = gm2_2L + ((3/20.D0*
     &         (Sf4857*Sf4862*Sf4863)/MHiggs2(4)**2 + 
     &        (Sf4828*(3/160.D0*(Sf4841*Sf4847) + 
     &             3/40.D0*(Sf4847*Sf4848)))/MHiggs2(4))*USbot11*
     &      UStop21)/(Sf4831*Ss25**2*SW2)

	gm2_2L = gm2_2L + (1/1440.D0*
     &       (CW2*Sf163*Sf939)/(MW2**3*Sf166**3) + 
     &      (-(3/80.D0*(Sf6419*Ss275)/MHiggs2(4)) - 
     &         1/160.D0*(Sf5223*Ss295)/MHiggs2(4))/
     &       (Sf1184*Sf154*Sf2206*Sf2270*Sf4831*Sf4867))/
     &    (CW2*Ss25**3*SW2)

	gm2_2L = gm2_2L + ((-(1/120.D0*(Sf461*Sf7529)) + 
     &        1/1440.D0*(Sf34*Sf7532) - 
     &        1/720.D0*(MW2*Sf7555*MSneu2(1)*MStau2(1))/Sf166**3)*
     &      UStau11**2)/(MW2**3*SW2*MSl2Diff(1,1)**2)

	gm2_2L = gm2_2L + (1/480.D0*
     &       (CB*Sf6393*Ss455)/
     &        (CW2*MW2**2*Sf1195**3*Sf166**3*Sf2320**3*Sf2334**3*
     &          Sf4820**3*Sf4906**3*Ss241*Ss242*TB*TB2) - 
     &      3/80.D0*(Sf5231*TB*MHiggs2(4))/SB)/(CB*Ss25**3*SW2)

	gm2_2L = gm2_2L + (Sf162*TB*
     &      (3/80.D0*Sf226/Sf166**2 + 
     &        3/20.D0*(MW2*Sf233*Sf249)/(Sf154*MHiggs2(4)**2)))/
     &    (MW2*Ss25**3*SW2)

	gm2_2L = gm2_2L + (1/1440.D0*
     &       (Sf1192*Sf1728)/(MW2**3*Sf1195**3) - 
     &      3/160.D0*(Sf2198*(Sf2215*Sf2246 + 4*Sf2222*Sf2246) + 
     &         Sf2263*(Sf2278*Sf2297 + 4*Sf2285*Sf2297) + 
     &         Sf4828*(Sf4841*Sf4852 + 4*Sf4848*Sf4852) + 
     &         Sf4864*(Sf4875*Sf4886 + 4*Sf4882*Sf4886)))/
     &    (Ss25**3*SW2)

	gm2_2L = gm2_2L + (CA*Sf7701*Sf9312*Ss46*Ss47*
     &      (-(1/5760.D0*(Sf10025 + Sf9928)/Ss49) + 
     &        1/5760.D0*(Sf9478 - Sf9738 + 3*Sf9767*Ss930 + 
     &           Sf9533*(-3*Ss930 - 36*Ss928*Ss949 - 3*Ss951) + 
     &           Sf9767*(36*Ss928*Ss949 + 3*Ss951)))*MStau2(1)*
     &      MStau2(2))/(CB2*CW**3*SW2*MStau2(3)**4)

	gm2_2L = gm2_2L + (Sf1191*TB*
     &      (3/80.D0*Sf1218/Sf1195**2 + 
     &        3/20.D0*(MW2*Sf1224*Sf1231)/(Sf1184*MHiggs2(4)**2)))/
     &    (MW2*Ss25**3*SW2)

	gm2_2L = gm2_2L + ((3/160.D0*(Sf5028*Sf5036)/Ss25**2 + 
     &        3/40.D0*(Sf5037*Sf5044)/Ss25**2 - 
     &        3/80.D0*(MW2*Sf4903*Sf4905*MSbot2(1)))*USbot11*
     &      UStop21)/(MW2**3*Sf4820**2*SW2*TB)

	gm2_2L = gm2_2L + (CA*Sf7701*Sf9312*
     &      (-(1/160.D0*(Ss78*
     &              (3*Sf9610 - Sf9327*Sf9337*MStau2(1) + 
     &                Sf9358*Sf9373*MStau2(2)))/(Sf9331*Ss49**2))+
     &          1/1920.D0*((Sf10004 + Sf9903)*Ss46*Ss47*
     &            (Ss930 + 12*Ss928*Ss949 + Ss951)*MStau2(1)*
     &            MStau2(2))/(MZ2*Ss49*MStau2(3)**4)))/
     &    (CB2*CW**3*SW2*MHiggs2(2))

	gm2_2L = gm2_2L + (CA*Sf7701*Sf9312*
     &      (-(1/480.D0*((Sf10054 - Sf9977)*Ss46*Ss47*
     &             (Ss928 + Ss949))) + 
     &        1/1440.D0*(Ss78*
     &            (Sf9685 + Sf9817 - 
     &              3*(Sf9707*Ss930 + 2*Sf9838*Ss949 + 
     &                 Sf9655*(2*Ss928 + 2*Ss949) + 
     &                 Ss928*
     &                  (2*Sf9838 + 12*Sf9416*Ss949 + 
     &                    12*Sf9707*Ss949) + Sf9707*Ss951 + 
     &                 Sf9416*(Ss930 + Ss951))))/Ss49)*MStau2(1)*
     &      MStau2(2))/(CB2*CW**3*SW2*MStau2(3)**3)

	gm2_2L = gm2_2L + ((-(1/240.D0*(Sf1894*Sf7670)/Sf166**3) + 
     &        1/120.D0*(Sf1195**3*Sf2022*Sf7700))*MSneu2(1))/
     &    (MW2**2*Sf1195**3*SW2)

	gm2_2L = gm2_2L + ((-(3/80.D0*
     &           (MW2*Sf4864*Sf5092)/
     &            (Sf4867*Sf4906**2*Ss241**2*Ss242**2*TB)) + 
     &        1/240.D0*Sf5445/(Sf4906**3*TB2))*MSbot2(2)*MStop2(2))/
     &    (MW2**2*SW2*MSq2Diff(2,2))

	gm2_2L = gm2_2L + 1/720.D0*
     &    ((Sf4906*Sf8094 + Sf8121)*MSbot2(2)*MStop2(2)*
     &       USbot21**2*UStop21**2)/
     &     (MW2**2*Sf4906**3*SW2*MSq2Diff(2,2)**2)

	gm2_2L = gm2_2L + 1/720.D0*
     &    ((-(Sf395*MSneu2(1)) - Sf437*MStau2(1) - 
     &          (Sf7570 - 3*Sf7562*Ss930*MStau2(1))*UStau11**2)/
     &        Sf166**3 + (-(Sf1332*MSneu2(1)) - Sf1360*MStau2(2) - 
     &          (Sf7628 - 3*Sf7620*Ss951*MStau2(2))*UStau21**2)/
     &        Sf1195**3)/(MW2**2*SW2)

	gm2_2L = gm2_2L + (1/180.D0*(Sf8642*Log(MZ2))/MW2 + 
     &      1/240.D0*(Sf4814*Log(MSbot2(1)) + 
     &          Sf4819*Log(MSbot2(2)) - Sf1790*Log(MStau2(1)) - 
     &          Sf2021*Log(MStau2(2)))/(CB*MW2**3*Ss25**2))/
     &    (MZ2**5*SW2**2)

	gm2_2L = gm2_2L + (CA*Sf7509*Sf8897*Ss78*
     &      (1/960.D0*(Sf9005*Sf9011 - 8*Sf9104*Sf9107 + 
     &            2*Sf8910*Sf8942*MSbot2(1) - 
     &            2*Sf8953*Sf8965*MSbot2(2))/(MZ2**2*Sf8946) - 
     &        1/1440.D0*(3*Sf9051 + Sf8916*Sf8998*MSbot2(1) - 
     &            Sf8966*Sf8976*MSbot2(2))/(Sf8970*MHiggs2(2))))/
     &    (CB2*CW**3*Ss49**2*SW2)

	gm2_2L = gm2_2L + (CA*Sf7509*Sf8897*Ss78*
     &      (1/960.D0*(Sf9004*Sf9005 + 8*Sf9102*Sf9104 - 
     &            2*Sf8910*Sf8915*MSbot2(1) + 
     &            2*Sf8953*Sf8956*MSbot2(2))/(MZ2**2*Sf8904) - 
     &        1/1440.D0*(Sf8900*Sf8988*MSbot2(1) + 
     &            Sf8957*Sf8961*MSbot2(2))/(MZ2*Sf8946)))/
     &    (CB2*CW**3*Ss49**2*SW2)

	gm2_2L = gm2_2L + ((-(3/80.D0*
     &           (MW2*Sf2198*Sf2438)/
     &            (Sf2206*Sf2320**2*Ss241**2*Ss242**2*TB)) + 
     &        1/240.D0*Sf2928/(Sf2320**3*TB2))*MSbot2(1)*MStop2(1))/
     &    (MW2**2*SW2*MSq2Diff(1,1))

	gm2_2L = gm2_2L - 1/360.D0*
     &    (MW2*Sf7497*(Sf7509**2*
     &          (Sf7510*Sf7511 + 12*Sf7514*Sf7515) + 
     &         3*Sf7701**2*(Sf7702*Sf7703 + 12*Sf7706*Sf7707) + 
     &         Sf8520**2*(Sf8521*Sf8522 + 12*Sf8525*Sf8526) + 
     &         (Sf7509**2*(12*Sf7516*Sf7518 + Sf7520*MSbot2(4)))/
     &          MSbot2(3)**2))/(CW2**2*MZ2**4*SW2)

	gm2_2L = gm2_2L + (CA*Ss78*
     &      (1/480.D0*(Sf10435*Sf8520*
     &            (Sf10506*Sf10513 - 8*Sf10328*Sf10593 + 
     &              2*Sf10456*Sf10459*MStop2(1)))/
     &          (MZ2**2*Sf10276*Ss49**2) + 
     &        1/48.D0*(Sf10243*Sf8532*MStop2(2))/(Sf113*Ss986)))/
     &    (CB*CW**3*SB*SW2)

	gm2_2L = gm2_2L + (Sf8520*Ss78*
     &      (1/360.D0*(CA*Sf10435*
     &            (Sf10559*(-6*Ss969 - 6*Ss988) + 
     &              12*Sf10494*Ss969*Ss988 + 
     &              Sf10494*(Ss971 + Ss990))*MStop2(1)*MStop2(2))/
     &          (Sf10267*Sf10444*Sf8929*Sf8930*MHiggs2(2)) - 
     &        1/240.D0*(SA*Sf10333*MStop2(3)*(4*Sf8527 - MStop2(4)))/
     &          MHiggs2(1)))/(CB*CW**3*MZ2*SB*SW2*MStop2(3))

	gm2_2L = gm2_2L + (CA*Sf10435*Sf8520*Ss78*
     &      (1/480.D0*(Sf10505*Sf10506 + 8*Sf10328*Sf10590 + 
     &            2*Sf10466*Sf10469*MStop2(2))/(MZ2**2*Sf10267) - 
     &        1/720.D0*(3*Sf10551 + Sf10441*Sf10498*MStop2(1) - 
     &            Sf10475*Sf10484*MStop2(2))/(Sf10478*MHiggs2(2))))
     &     /(CB*CW**3*SB*Ss49**2*SW2)

	gm2_2L = gm2_2L + 1/360.D0*
     &    (Sf8520*Ss78*((SA*Sf10270*Sf9131*Sf9132*
     &            (Sf10311*(-6*Ss969 - 6*Ss988) + 
     &              12*Sf10302*Ss969*Ss988 + 
     &              Sf10302*(Ss971 + Ss990)))/
     &          (Sf10267*Sf10296*Ss44**2*MHiggs2(1)) + 
     &         (CA*Sf10435*
     &            (-12*Sf10488*Ss969*Ss988 + 
     &              Sf10567*(6*Ss969 + 6*Ss988) + 
     &              Sf10488*(-Ss971 - Ss990)))/
     &          (Sf10276*Sf10478*Sf8929*Sf8930*MHiggs2(2)))*
     &       MStop2(1)*MStop2(2))/(CB*CW**3*MZ2*SB*SW2*MStop2(3))

	gm2_2L = gm2_2L + 3/80.D0*
     &    (TB*MSneu2(1)*((Sf162*Sf315*MStau2(1))/
     &          (Sf166**2*MSl2Diff(1,1)) + 
     &         (Sf1191*Sf1275*Sf154*MStau2(2))/
     &          (Sf1184*Sf1195**2*MSl2Diff(2,1))))/
     &     (Sf154*Ss25**2*SW2)

	gm2_2L = gm2_2L + (Sf8532*
     &      (-(1/240.D0*(CA*Sf10243*
     &              (Sf10245*Sf10246*
     &                 (Sf8678*Ss78 + 4*Sf8681*Ss78) - 
     &                8*Sf10250*Ss46*Ss47*MStop2(2)))/Ss49) + 
     &        1/120.D0*(SA*Ss997*
     &            (Sf10261*Ss78 - 4*Sf10250*Ss46*Ss47*MStop2(2)))/
     &          Ss44))/(CB*CW**3*MZ2*SB*Sf113*SW2)

	gm2_2L = gm2_2L + 1/120.D0*
     &    ((Sf7733*Ss928)/MSl2Diff(1,1)**2 + 
     &       (Sf166**3*Sf7755*Ss949)/(Sf1195**3*MSl2Diff(2,1)**2))/
     &     (CW2**2*MW2**2*MZ2**3*Sf166**3*SW2*MStau2(3)**2)

	gm2_2L = gm2_2L + (MW2*
     &      (1/360.D0*(-(Sf7501*Sf7504*Sf7506) - 
     &           Sf7521*Sf7523*Sf7525 - 3*Sf7579*Sf7581*Sf7583 - 
     &           3*Sf7713*Sf7734*Sf7736 - Sf7990*Sf7992*Sf7994 - 
     &           Sf8532*Sf8678*Sf8680 + 4*Sf8721*Ss474 - 
     &           Sf7497*(3*Sf7498*Sf7499 + 
     &              2*(-(Sf7509**2*Sf7513) - 3*Sf7701**2*Sf7705 - 
     &                 Sf8520**2*Sf8524 + 
     &                 19*Sf7501**2*MSbot2(1)**2 + 
     &                 19*Sf7521**2*MSbot2(2)**2 - 
     &                 18*Sf7500*MSneu2(1) + 
     &                 19*(3*MSneu2(1)**2 + 
     &                    3*Sf7579**2*MStau2(1)**2 + 
     &                    3*Sf7713**2*MStau2(2)**2 + 
     &                    Sf7990**2*MStop2(1)**2 + 
     &                    Sf8532**2*MStop2(2)**2)))) + 
     &        1/180.D0*(Sf7497*Sf8520**2*Sf8523*MStop2(1)*MStop2(2))/
     &          MStop2(3)**2))/(CW2**2*MZ2**3*SW2)

	gm2_2L = gm2_2L + (MW2*
     &      (-(1/360.D0*(Sf7497*Sf8520**2*
     &              (12*Sf8527*Sf8529 + Sf8531*MStop2(4)))/
     &            (MZ2**4*SW2*MStop2(3)**2)) - 
     &        1/60.D0*(Sf8715*Sf8720 + 
     &            6*Ss922*(2*Sf7501**2*MSbot2(1) + 
     &               2*Sf7521**2*MSbot2(2) + Sf7509**2*MSbot2(4) + 
     &               6*Sf7579**2*MStau2(1) + 
     &               6*Sf7713**2*MStau2(2) + 
     &               3*Sf7701**2*MStau2(4) + 
     &               2*Sf7990**2*MStop2(1) + 
     &               2*Sf8532**2*MStop2(2) + Sf8520**2*MStop2(4)))/
     &          MZ2**2))/CW2**2

	gm2_2L = gm2_2L + (Sf8532*
     &      (-(1/360.D0*(SA*Ss46*Ss47*
     &              (Ss268 + 12*Ss263*Ss988 + Ss990)*Ss997)/
     &            (Ss44*Ss995)) + 
     &        1/360.D0*(CA*Sf10243*Ss78*
     &            (6*Sf10249*Ss988 + 6*Sf10248*Ss988**2 + 
     &              Sf10248*Ss990 - 33*Sf10247*MStop2(2)))/
     &          (MZ2*Sf113*Ss986*MHiggs2(2)))*MStop2(2))/
     &    (CB*CW**3*SB*SW2)

	gm2_2L = gm2_2L + 3/160.D0*
     &    ((Sf1191*(Sf1222*Sf1370 + 4*Sf1370*Sf1371) + 
     &         (Sf1184*Sf162*(Sf231*Sf453 + 4*Sf453*Sf454))/Sf154)*
     &       TB)/(Sf1184*Ss25**3*SW2*MHiggs2(4))

	gm2_2L = gm2_2L + 1/240.D0*
     &    ((-(Sf7296*Log(MW2)) + 
     &          9*MW2**2*MZ2*SW2*TB*
     &           (Sf5231*Log(MW2)**2 + 
     &             Log(MHiggs2(4))*
     &              (Sf7211 - Sf5231*Log(MHiggs2(4))))*MHiggs2(4))/
     &        (CB*SB*Ss25**3) + (Sf1176*MSneu2(1))/Ss25**2)/
     &     (MW2**2*MZ2*SW2**2)

	gm2_2L = gm2_2L + (Sf7521*Ss78*
     &      ((SA*(1/480.D0*(Sf7523*Ss915)/Ss44 + 
     &             1/120.D0*(Sf7526*Ss915)/Ss44 + 
     &             1/60.D0*(Sf8890*Ss915)/Ss44))/(MZ2*Sf26) + 
     &        (CA*Ss905*(-(1/480.D0*Ss909/MHiggs2(2)) - 
     &             1/120.D0*Ss910/MHiggs2(2) - 1/60.D0*Ss912/MHiggs2(2)
     &             ))/(Ss49*Ss903))*MSbot2(2))/(CB2*CW**3*SW2)

	gm2_2L = gm2_2L + (MW2*
     &      (1/90.D0*(Sf7501*Sf7507*Sf7508) + 
     &        1/90.D0*(Sf7521*Sf7526*Sf7527) + 
     &        1/30.D0*(Sf7579*Sf7584*Sf7585) + 
     &        1/30.D0*(Sf7713*Sf7737*Sf7738) + 
     &        1/90.D0*(Sf7990*Sf7995*Sf7996) + 
     &        1/90.D0*(Sf8532*Sf8681*Sf8682) + 
     &        1/180.D0*(Sf7497*Sf7509**2*Sf7512*MSbot2(1)*MSbot2(2))/
     &          (SW2*MSbot2(3)**2)))/(CW2**2*MZ2**3)

	gm2_2L = gm2_2L + ((-(3/80.D0*
     &           (MW2*Sf2263*Sf2652)/
     &            (Sf2270*Sf2334**2*Ss241**2*Ss242**2*TB)) + 
     &        1/240.D0*Sf3055/(Sf2334**3*TB2))*MSbot2(2)*MStop2(1))/
     &    (MW2**2*SW2*MSq2Diff(2,1))

	gm2_2L = gm2_2L + (SA*Sf7509*Sf9108*Ss78*
     &      (-(1/960.D0*(Sf9005*Sf9209 - 8*Sf9104*Sf9286 + 
     &              2*Sf8910*Sf9144*MSbot2(1) - 
     &              2*Sf8953*Sf9155*MSbot2(2))/(MZ2**2*Sf8946)) + 
     &        1/1440.D0*(3*Sf9256 + Sf9118*Sf9195*MSbot2(1) - 
     &            Sf9156*Sf9166*MSbot2(2))/(Sf9160*MHiggs2(1))))/
     &    (CB2*CW**3*Ss44**2*SW2)

	gm2_2L = gm2_2L + (SA*Sf7509*Sf9108*Ss78*
     &      (-(1/960.D0*(Sf9005*Sf9202 + 8*Sf9104*Sf9283 - 
     &              2*Sf8910*Sf9117*MSbot2(1) + 
     &              2*Sf8953*Sf9147*MSbot2(2))/(MZ2**2*Sf8904)) + 
     &        1/1440.D0*(Sf8900*Sf9175*MSbot2(1) + 
     &            Sf8957*Sf9151*MSbot2(2))/(MZ2*Sf8946)))/
     &    (CB2*CW**3*Ss44**2*SW2)

	gm2_2L = gm2_2L + (CA*Sf7521*Ss78*Ss905*
     &      (1/720.D0*(MSbot2(2)*
     &            (6*Sf8889*Ss906 + 6*Sf8888*Ss906**2 + 
     &              Sf8888*Ss908 - 33*Sf8887*MSbot2(2)))/
     &          (Ss903*MHiggs2(2)) - 
     &        1/480.D0*(Sf7523*Sf8885*Sf8886 + 
     &            4*Sf7526*Sf8885*Sf8886 - 8*Sf8890*MSbot2(2))/Ss49
     &        ))/(CB2*CW**3*MZ2*Sf26*SW2)

	gm2_2L = gm2_2L + (Sf7501*Ss78*
     &      ((SA*(1/480.D0*(Sf7504*Ss897)/Ss44 + 
     &             1/120.D0*(Sf7507*Ss897)/Ss44 + 
     &             1/60.D0*(Sf8878*Ss897)/Ss44))/(MZ2*Sf11) + 
     &        (CA*Ss884*(-(1/480.D0*Ss888/MHiggs2(2)) - 
     &             1/120.D0*Ss889/MHiggs2(2) - 1/60.D0*Ss891/MHiggs2(2)
     &             ))/(Ss49*Ss878))*MSbot2(1))/(CB2*CW**3*SW2)

	gm2_2L = gm2_2L + ((-(3/80.D0*(Sf4911*Sf4917)/(Sf4906**2*TB)) - 
     &        3/40.D0*(MW2*Sf4864*Sf4874*MSbot2(2))/
     &          (Sf4867*Ss241*Ss242))*MStop2(2)*USbot21*
     &      UStop21)/(MW2**2*SW2*MSq2Diff(2,2))

	gm2_2L = gm2_2L - 1/720.D0*
     &    ((Sf2334*Sf7885 + Sf7912)*MSbot2(2)*MStop2(1)*
     &       USbot21**2*UStop11**2)/
     &     (MW2**2*Sf2334**3*SW2*MSq2Diff(2,1)**2)

	gm2_2L = gm2_2L + 1/240.D0*
     &    ((-6*Sf1 + Sf1181*Log(MSneu2(1)))/MZ2 + 
     &       (Sf7234*Log(MStop2(1)) + Sf7239*Log(MStop2(2)))/
     &        (MW2**3*MZ2**5*SB*Ss25**2))/SW2**2

	gm2_2L = gm2_2L + (CA*Sf7501*Ss78*Ss884*
     &      (1/720.D0*(MSbot2(1)*
     &            (6*Sf8877*Ss885 + 6*Sf8876*Ss885**2 + 
     &              Sf8876*Ss887 - 33*Sf8875*MSbot2(1)))/
     &          (Ss878*MHiggs2(2)) - 
     &        1/480.D0*(Sf7504*Sf8873*Sf8874 + 
     &            4*Sf7507*Sf8873*Sf8874 - 8*Sf8878*MSbot2(1))/Ss49
     &        ))/(CB2*CW**3*MZ2*Sf11*SW2)

	gm2_2L = gm2_2L + (1/120.D0*
     &       (-(Sf113*Sf117**2) - Sf11*Sf16**2 - Sf26*Sf30**2 - 
     &         3*Sf39*Sf43**2 - 3*Sf61*Sf65**2 - Sf82*Sf87**2) + 
     &      1/240.D0*(Sf1767 + Sf1787 + Sf4788 + Sf4809)/
     &        (CB*MW2**2*Ss25**2))/(MZ2**3*SW2**2)

	gm2_2L = gm2_2L + (CA*(1/240.D0*
     &         (Sf10243*Sf8532*
     &            (Sf10244*Ss78*(Ss991 + 4*Ss992) - 
     &              8*Ss46*Ss47*Ss994*MStop2(2)))/(Ss49*Ss986) + 
     &        1/120.D0*(Sf10435*(-Sf10462 + Sf10465)*Sf8520*Sf8527*
     &            Ss78)/(MZ2**2*MStop2(3))))/
     &    (CB*CW**3*SB*SW2*MHiggs2(2))

	gm2_2L = gm2_2L + (CA*(1/240.D0*
     &         (Sf10223*Sf7990*
     &            (Sf10224*Ss78*(Ss972 + 4*Ss973) - 
     &              8*Ss46*Ss47*Ss975*MStop2(1)))/
     &          (Ss964*MHiggs2(2)) + 
     &        1/360.D0*(Sf8532*
     &            (Ss78*Ss987*(-3 - 6*Ss269 - 6*Ss988) + 
     &              Sf10243*
     &               (Ss274*Ss46*Ss47 + 
     &                 Ss46*Ss47*(12*Ss269*Ss988 + Ss990)))*
     &            MStop2(2))/Ss986))/(CB*CW**3*SB*Ss49*SW2)

	gm2_2L = gm2_2L + ((1/720.D0*
     &         (Sf8043*MSbot2(1)*MStop2(2)*USbot11**2)/
     &          MSq2Diff(1,2)**2 - 
     &        1/720.D0*(Sf4820**3*Sf8089*USbot21**2)/Sf4906**2)*
     &      UStop21**2)/(MW2**2*Sf4820**3*SW2)

	gm2_2L = gm2_2L + (CA*Sf10223*Sf7990*
     &      (1/360.D0*(Ss78*MStop2(1)*
     &            (6*Sf10229*Ss969 + 6*Sf10228*Ss969**2 + 
     &              Sf10228*Ss971 - 33*Sf10227*MStop2(1)))/
     &          (Ss964*MHiggs2(2)) - 
     &        1/240.D0*(Sf10225*Sf10226*
     &             (Sf7992*Ss78 + 4*Sf7995*Ss78) - 
     &            8*Sf10230*Ss46*Ss47*MStop2(1))/Ss49))/
     &    (CB*CW**3*MZ2*SB*Sf82*SW2)

	gm2_2L = gm2_2L + ((-(3/80.D0*(Sf2327*Sf2333)/(Sf2320**2*TB)) - 
     &        3/40.D0*(MW2*Sf2198*Sf2214*MSbot2(1))/
     &          (Sf2206*Ss241*Ss242))*MStop2(1)*USbot11*
     &      UStop11)/(MW2**2*SW2*MSq2Diff(1,1))

	gm2_2L = gm2_2L + (SA*Sf7521*Ss78*Ss915*
     &      (1/96.D0 + (Sf26*
     &           (-(1/480.D0*Ss916/MHiggs2(1)) + 
     &             (-(1/120.D0*Ss918) - 1/60.D0*Ss919)/MHiggs2(1)))/
     &         Ss44)*MSbot2(2))/(CB2*CW**3*Sf26*Ss913*SW2)

	gm2_2L = gm2_2L + (CA*Sf7701*Sf9312*Ss78*
     &      (3/320.D0*(Sf9381*Sf9382 + 8*Sf9864*Sf9866 - 
     &            2*Sf9322*Sf9326*MStau2(1) + 
     &            2*Sf9345*Sf9348*MStau2(2))/(MZ2**2*Sf9317) - 
     &        1/160.D0*(Sf9313*Sf9371*MStau2(1) + 
     &            Sf9349*Sf9353*MStau2(2))/(MZ2*Sf9344)))/
     &    (CB2*CW**3*Ss49**2*SW2)

	gm2_2L = gm2_2L + (-(1/4320.D0*(CB*Sf7297*Sf7496)/MZ2**4) + 
     &      1/240.D0*Sf7295/(MZ2*SB*Ss25**3))/(CB*MW2**2*SW2**2)

	gm2_2L = gm2_2L + (Sf7701*
     &      (1/480.D0*(CA*Sf9312*(Sf9592 + Sf9796)*Ss46*Ss47*
     &           (Ss928 + Ss949)*MStau2(1)*MStau2(2)) + 
     &        1/160.D0*(SA*Sf10055*Ss78*
     &            (3*Sf10190 + Sf10064*Sf10130*MStau2(1) - 
     &              Sf10094*Sf10104*MStau2(2)))/
     &          (Sf10098*Ss44**2*MHiggs2(1))))/(CB2*CW**3*SW2)

	gm2_2L = gm2_2L + (SA*Sf8520*Ss78*
     &      (1/720.D0*(-(Sf10333*MStop2(1)*
     &              (6*Sf10363*Ss474 + Sf10377*(3 + 6*Ss969) + 
     &                Sf10385*(Ss484 + Ss971) + 
     &                Ss474*
     &                 (12*Sf10385*Ss969 - 6*Sf10277*MStop2(1) + 
     &                   6*Sf9113*MStop2(1)**2))) - 
     &           Sf10333*(6*Sf10245*Sf10246*Sf10362*Ss474 + 
     &              12*Sf10373*Ss474*Ss988 + 
     &              Sf10376*(3 + 6*Ss988) + Sf10373*(Ss484 + Ss990)
     &              )*MStop2(2) + 
     &           6*Sf10270*Sf8526*
     &            (-3*MHiggs2(1)*MStop2(1)**2 - 2*MStop2(1)**3 + 
     &              6*MStop2(1)**2*MStop2(2) - 
     &              2*(Sf10245*Sf10246 + 3*MStop2(1))*MStop2(2)**2)
     &           ) + 1/360.D0*
     &         (Sf10270*Sf9131*Sf9132*
     &            (-12*Sf10292*Ss969*Ss988 + 
     &              Sf10320*(6*Ss969 + 6*Ss988) + 
     &              Sf10292*(-Ss971 - Ss990))*MStop2(1)*MStop2(2))/
     &          (Sf10286*MHiggs2(1)*MStop2(3))))/
     &    (CB*CW**3*MZ2*SB*Sf10276*Ss44**2*SW2)

	gm2_2L = gm2_2L + (SA*Sf8532*Ss997*
     &      (-(1/120.D0*(Sf10257*Ss78 - 
     &              4*Ss1001*Ss46*Ss47*MStop2(2))/Ss44) + 
     &        1/360.D0*(MStop2(2)*
     &            (-6*Sf10254*Ss78*Ss988 - 
     &              Sf10253*Ss46*Ss47*(6*Ss988**2 + Ss990) + 
     &              33*Sf10251*Ss78*MStop2(2)))/(MZ2*Sf113)))/
     &    (CB*CW**3*SB*Ss995*SW2*MHiggs2(1))

	gm2_2L = gm2_2L + (Sf2263*
     &      (-(3/40.D0*(Sf2277*MSbot2(2)*MStop2(1))/
     &            (MW2*Ss241*Ss242*MSq2Diff(2,1))) + 
     &        3/80.D0*(Sf2269*MSbot2(2) - Sf2274*MStop2(1))/Ss25**3)*
     &      USbot21*UStop11)/(Sf2270*SW2)

	gm2_2L = gm2_2L + 1/60.D0*
     &    (Sf10417*Sf8520*Ss78*
     &       (SA*Sf10333*(Sf10414/MHiggs2(1)**2 - 
     &            Sf10418/MHiggs2(1)**2) + 
     &         CA*Sf10435*(Sf10602/MHiggs2(2)**2 - 
     &            Sf10603/MHiggs2(2)**2)))/
     &     (CB*CW**3*MZ2**2*SB*SW2*MStop2(3))

	gm2_2L = gm2_2L + 1/120.D0*
     &    (SA*Sf10270*Sf8520*Sf8526*Ss78*
     &       ((Sf10225*Sf10226*(2*Sf9131*Sf9132 + MStop2(1)) + 
     &            5*MStop2(2)**2)/Sf10267 + 
     &         (5*MStop2(1)**2 + 
     &            Sf10245*Sf10246*(2*Sf9131*Sf9132 + MStop2(2)))/
     &          Sf10276))/(CB*CW**3*SB*Ss44**2*SW2)

	gm2_2L = gm2_2L + (CA*Sf7701*Sf9312*Ss78*
     &      (3/320.D0*(Sf9382*Sf9388 - 8*Sf9866*Sf9869 + 
     &            2*Sf9322*Sf9340*MStau2(1) - 
     &            2*Sf9345*Sf9357*MStau2(2))/(MZ2**2*Sf9344) - 
     &        1/160.D0*(3*Sf9625 + Sf9327*Sf9375*MStau2(1) - 
     &            Sf9358*Sf9368*MStau2(2))/(Sf9362*MHiggs2(2))))/
     &    (CB2*CW**3*Ss49**2*SW2)

	gm2_2L = gm2_2L + (Sf8520*Ss78*
     &      (-(1/720.D0*(CA*Sf10435*
     &              (Sf10436*Sf10490*MStop2(1) + 
     &                (Sf10470*Sf10474 + 
     &                   Sf10601*(3 + 6*Ss474 + 6*Ss988))*MStop2(2)
     &                ))/(Sf10276*Ss49**2)) + 
     &        1/720.D0*(SA*(MStop2(1)*
     &               (-(Sf10333*
     &                    (6*Sf10225*Sf10226*Sf10342*Ss474 + 
     &                      12*Sf10346*Ss474*Ss969 + 
     &                      Sf10350*(3 + 6*Ss969) + 
     &                      Sf10346*(Ss484 + Ss971))) - 
     &                 12*Sf10225*Sf10226*Sf10270*Sf8526*MStop2(1))
     &                - (Sf10333*
     &                  (6*Sf10343*Ss474 + 
     &                    12*Sf10384*Ss474*Ss988 + 
     &                    Sf10348*(3 + 6*Ss988) + 
     &                    Sf10384*(Ss484 + Ss990)) + 
     &                 36*Sf10270*Sf8526*MStop2(1)**2)*MStop2(2) + 
     &              6*(Sf10272*Sf10333*Ss474 + 
     &                 Sf10270*Sf8526*(-3*MHiggs2(1) + 6*MStop2(1))
     &                 )*MStop2(2)**2 - 
     &              6*(2*Sf10270*Sf8526 + Sf10333*Sf9113*Ss474)*
     &               MStop2(2)**3))/(Sf10267*Ss44**2)))/
     &    (CB*CW**3*MZ2*SB*SW2)

	gm2_2L = gm2_2L + (SA*Sf7501*Ss78*Ss897*
     &      (1/96.D0 + (Sf11*
     &           (-(1/480.D0*Ss898/MHiggs2(1)) + 
     &             (-(1/120.D0*Ss900) - 1/60.D0*Ss901)/MHiggs2(1)))/
     &         Ss44)*MSbot2(1))/(CB2*CW**3*Sf11*Ss892*SW2)

	gm2_2L = gm2_2L + (CB*(1/240.D0*Sf7210 + 1/240.D0*Sf7229) - 
     &      1/1440.D0*Sf7290)/(CB*MW2**2*MZ2**3*SB*Ss25**2*SW2**2)

	gm2_2L = gm2_2L + 1/120.D0*
     &    ((Re(Sf21*Sf23*Sf24*Sf25**2*
     &             Log(-(1/2.D0*Sf19/(MSbot(1)*MSbot(2)))))*MSbot(1)*
     &           MSbot(2) + 
     &          3*Re(Sf56*Sf58*Sf59*Sf60**2*
     &             Log(-(1/2.D0*Sf54/(MStau(1)*MStau(2)))))*MStau(1)*
     &           MStau(2) + 
     &          Re(Sf108*Sf110*Sf111*Sf112**2*
     &             Log(-(1/2.D0*Sf106/(MStop(1)*MStop(2)))))*
     &           MStop(1)*MStop(2))/MZ2**5 + 
     &       (-(Re(Sf37*Sf38*
     &               Log(-(1/2.D0*Sf33/(MSneu(1)*MStau(1))))*
     &               UStau11**2)*MSneu(1)*MStau(1)) - 
     &          Re(Sf50*Sf51*
     &             Log(-(1/2.D0*Sf46/(MSneu(1)*MStau(2))))*
     &             UStau21**2)*MSneu(1)*MStau(2) - 
     &          3*(Re(Sf72*Sf73*
     &                Log(-(1/2.D0*Sf68/(MSbot(1)*MStop(1))))*
     &                USbot11**2*UStop11**2)*MSbot(1)*
     &              MStop(1) + 
     &             Re(Sf80*Sf81*
     &                Log(-(1/2.D0*Sf76/(MSbot(2)*MStop(1))))*
     &                USbot21**2*UStop11**2)*MSbot(2)*
     &              MStop(1) + 
     &             Re(Sf94*Sf95*
     &                Log(-(1/2.D0*Sf90/(MSbot(1)*MStop(2))))*
     &                USbot11**2*UStop21**2)*MSbot(1)*
     &              MStop(2) + 
     &             Re(Sf102*Sf103*
     &                Log(-(1/2.D0*Sf98/(MSbot(2)*MStop(2))))*
     &                USbot21**2*UStop21**2)*MSbot(2)*
     &              MStop(2)))/(MW2**3*MZ2))/SW2**2

	gm2_2L = gm2_2L + 1/480.D0*
     &    (SA*Sf10270*Sf8520*Ss78*
     &       ((-8*Sf10328*Sf10332 + 
     &            Sf10245*Sf10246*
     &             (Sf10280*Sf8521 + 
     &               4*Sf8526*MHiggs2(1)*MStop2(2)**2) + 
     &            MStop2(1)*
     &             (-(Sf10278*Sf8521) - 4*Sf10282*Sf8526 + 
     &               MStop2(1)*
     &                (-(Sf10277*Sf8521) + 
     &                  Sf8521*Sf9113*MStop2(1) + 
     &                  4*Sf8526*MHiggs2(1)*
     &                   (-2*MStop2(2) - MStop2(3)))))/Sf10276 + 
     &         (8*Sf10326*Sf10328 + 
     &            Sf10225*Sf10226*
     &             (Sf10263*Sf8521 + 
     &               4*Sf8526*MHiggs2(1)*MStop2(1)**2) + 
     &            MStop2(2)*
     &             (-(Sf10271*Sf8521) - 4*Sf10281*Sf8526 + 
     &               MStop2(2)*
     &                (-(Sf10272*Sf8521) + 
     &                  Sf8521*Sf9113*MStop2(2) + 
     &                  4*Sf8526*MHiggs2(1)*
     &                   (-2*MStop2(1) + MStop2(3)))))/Sf10267))/
     &     (CB*CW**3*MZ2**2*SB*Ss44**2*SW2)

	gm2_2L = gm2_2L - 1/720.D0*
     &    (((Sf7834*MSbot2(1)*MStop2(1)*USbot11**2)/
     &          MSq2Diff(1,1)**2 + 
     &         (Sf2320**3*Sf7880*USbot21**2)/Sf2334**2)*
     &       UStop11**2)/(MW2**2*Sf2320**3*SW2)

	gm2_2L = gm2_2L + 1/120.D0*
     &    (Sf7509*Sf9079*Ss78*
     &       (SA*Sf9108*(-(Sf9258/MHiggs2(1)**2) + 
     &            Sf9259/MHiggs2(1)**2) + 
     &         CA*Sf8897*(Sf9076/MHiggs2(2)**2 - 
     &            Sf9080/MHiggs2(2)**2)))/
     &     (CB2*CW**3*MZ2**2*SW2*MSbot2(3))

	gm2_2L = gm2_2L + ((-(3/80.D0*
     &           (MW2*Sf4828*Sf4966)/
     &            (Sf4820**2*Sf4831*Ss241**2*Ss242**2*TB)) + 
     &        1/240.D0*Sf5339/(Sf4820**3*TB2))*MSbot2(1)*MStop2(2))/
     &    (MW2**2*SW2*MSq2Diff(1,2))

	gm2_2L = gm2_2L + 1/720.D0*
     &    (SA*((Sf7501*Ss897*MSbot2(1)*
     &            (6*Sf8881*Ss78*Ss885 + 
     &              Sf8880*Ss46*Ss47*(6*Ss885**2 + Ss887) - 
     &              33*Sf8879*Ss78*MSbot2(1)))/(Sf11*Ss892) + 
     &         (Sf7521*Ss915*MSbot2(2)*
     &            (6*Sf8893*Ss78*Ss906 + 
     &              Sf8892*Ss46*Ss47*(6*Ss906**2 + Ss908) - 
     &              33*Sf8891*Ss78*MSbot2(2)))/(Sf26*Ss913)))/
     &     (CB2*CW**3*MZ2*SW2*MHiggs2(1))

	gm2_2L = gm2_2L - 3/80.D0*
     &    ((Sf4864*(Sf5138*MSbot2(2) + Sf5116*MStop2(2)))/
     &        (MW2*Sf4906**2*TB) + 
     &       (Sf162*TB*(Sf165*MSneu2(1) - Sf159*MStau2(1))*
     &          UStau11)/Sf154)/(Ss25**3*SW2)

	gm2_2L = gm2_2L - 1/720.D0*
     &    ((Sf4820*Sf8022 + Sf8059)*USbot11**2*UStop21**2)/
     &     (MW2**2*Sf4820**3*SW2)

	gm2_2L = gm2_2L + (Sf7713*Ss78*
     &      ((SA*(3/160.D0*(Sf7734*Ss958)/Ss44 + 
     &             3/40.D0*(Sf7737*Ss958)/Ss44 + 
     &             3/20.D0*(Sf9307*Ss958)/Ss44))/(MZ2*Sf61) + 
     &        (CA*Ss948*(-(3/160.D0*Ss952/MHiggs2(2)) - 
     &             3/40.D0*Ss953/MHiggs2(2) - 3/20.D0*Ss955/MHiggs2(2))
     &           )/(Ss49*Ss946))*MStau2(2))/(CB2*CW**3*SW2)

	gm2_2L = gm2_2L + (SA*Sf10055*Sf7701*Ss78*
     &      (1/160.D0*(3*Sf10175 - Sf10064*Sf10074*MStau2(1) + 
     &            Sf10094*Sf10120*MStau2(2))/(Sf10068*MHiggs2(1))-
     &          3/320.D0*(Sf10137*Sf9382 + 8*Sf10215*Sf9866 - 
     &            2*Sf10063*Sf9322*MStau2(1) + 
     &            2*Sf10085*Sf9345*MStau2(2))/(MZ2**2*Sf9317)))/
     &    (CB2*CW**3*Ss44**2*SW2)

	gm2_2L = gm2_2L + (MW2*Ss78*
     &      (1/9.D0*(Sf7501*MSbot2(1) + Sf7521*MSbot2(2) + 
     &            3*Sf7579*MStau2(1) + 3*Sf7713*MStau2(2) + 
     &            2*Sf7990*MStop2(1) + 2*Sf8532*MStop2(2))/MZ2**2+
     &          4/9.D0*(Sf7501*MSbot2(1)**2 + Sf7521*MSbot2(2)**2 + 
     &            3*Sf7579*MStau2(1)**2 + 3*Sf7713*MStau2(2)**2 + 
     &            2*Sf7990*MStop2(1)**2 + 2*Sf8532*MStop2(2)**2)/
     &          MZ2**3))/CW2

	gm2_2L = gm2_2L + (CA*Sf7509*Sf8897*Ss78*
     &      (1/240.D0*Sf9092 - 1/240.D0*(Sf8920*Sf9097)/Sf8970))/
     &    (CB2*CW**3*Sf8920*Ss49**2*SW2)

	gm2_2L = gm2_2L + (SA*Sf10055*Sf7701*Ss78*
     &      (-(3/320.D0*(Sf10144*Sf9382 - 8*Sf10218*Sf9866 + 
     &              2*Sf10082*Sf9322*MStau2(1) - 
     &              2*Sf10093*Sf9345*MStau2(2))/MZ2**2) + 
     &        1/160.D0*(Sf10112*Sf9313*MStau2(1) + 
     &            Sf10089*Sf9349*MStau2(2))/MZ2))/
     &    (CB2*CW**3*Sf9344*Ss44**2*SW2)

	gm2_2L = gm2_2L + (Sf7579*Ss78*
     &      ((SA*(3/160.D0*(Sf7581*Ss940)/Ss44 + 
     &             3/40.D0*(Sf7584*Ss940)/Ss44 + 
     &             3/20.D0*(Sf9295*Ss940)/Ss44))/(MZ2*Sf39) + 
     &        (CA*Ss927*(-(3/160.D0*Ss931/MHiggs2(2)) - 
     &             3/40.D0*Ss932/MHiggs2(2) - 3/20.D0*Ss934/MHiggs2(2))
     &           )/(Ss49*Ss921))*MStau2(1))/(CB2*CW**3*SW2)

	gm2_2L = gm2_2L - 3/80.D0*
     &    ((Sf2198*(Sf2550*MSbot2(1) + Sf2498*MStop2(1)))/
     &        Sf2320**2 + (Sf2263*
     &          (Sf2729*MSbot2(2) + Sf2690*MStop2(1)))/Sf2334**2)/
     &     (MW2*Ss25**3*SW2*TB)

	gm2_2L = gm2_2L + (CA*Sf7713*Ss78*Ss948*
     &      (1/80.D0*(MStau2(2)*
     &            (6*Sf9306*Ss949 + 6*Sf9305*Ss949**2 + 
     &              Sf9305*Ss951 - 33*Sf9304*MStau2(2)))/
     &          (Ss946*MHiggs2(2)) - 
     &        3/160.D0*(Sf7734*Sf9302*Sf9303 + 
     &            4*Sf7737*Sf9302*Sf9303 - 8*Sf9307*MStau2(2))/Ss49
     &        ))/(CB2*CW**3*MZ2*Sf61*SW2)

	gm2_2L = gm2_2L + 3/80.D0*
     &    ((Sf4864*(Sf4866*MSbot2(2) - Sf4871*MStop2(2))*
     &          USbot21*UStop21)/(Sf4867*Ss25**3) + 
     &       (Sf173*((Sf2198*MSq2Diff(1,1)*USbot11 + 
     &               Sf2263*MSq2Diff(2,1)*USbot21)*UStop11 + 
     &            (Sf4828*MSq2Diff(1,2)*USbot11 + 
     &               Sf4864*MSq2Diff(2,2)*USbot21)*UStop21))/
     &        (MW2*Ss25**2))/SW2

	gm2_2L = gm2_2L + 1/240.D0*
     &    ((-(Sf5266*MSbot2(1)) - Sf5301*MStop2(2))/Sf4820**3 + 
     &       (-(Sf5374*MSbot2(2)) - Sf5407*MStop2(2))/Sf4906**3)/
     &     (MW2**2*SW2*TB2)

	gm2_2L = gm2_2L + (CA*Sf7579*Ss78*Ss927*
     &      (1/80.D0*(MStau2(1)*
     &            (6*Sf9294*Ss928 + 6*Sf9293*Ss928**2 + 
     &              Sf9293*Ss930 - 33*Sf9292*MStau2(1)))/
     &          (Ss921*MHiggs2(2)) - 
     &        3/160.D0*(Sf7581*Sf9290*Sf9291 + 
     &            4*Sf7584*Sf9290*Sf9291 - 8*Sf9295*MStau2(1))/Ss49
     &        ))/(CB2*CW**3*MZ2*Sf39*SW2)

	gm2_2L = gm2_2L + 1/240.D0*
     &    (-4*(Sf108*Sf110*Sf112**2 + Sf21*Sf23*Sf25**2 + 
     &          3*Sf56*Sf58*Sf60**2) + 
     &       Re((Sf11*Sf14*Sf16**2*Log(-(1/2.D0*Sf12/MSbot2(1))))/
     &         Sf12) + Re((Sf26*Sf29*Sf30**2*
     &           Log(-(1/2.D0*Sf27/MSbot2(2))))/Sf27) + 
     &       3*Re((Sf39*Sf42*Sf43**2*
     &            Log(-(1/2.D0*Sf40/MStau2(1))))/Sf40) + 
     &       3*Re((Sf61*Sf64*Sf65**2*
     &            Log(-(1/2.D0*Sf62/MStau2(2))))/Sf62) + 
     &       Re((Sf82*Sf85*Sf87**2*Log(-(1/2.D0*Sf83/MStop2(1))))/
     &         Sf83) + Re((Sf113*Sf116*Sf117**2*
     &           Log(-(1/2.D0*Sf114/MStop2(2))))/Sf114))/
     &     (MZ2**4*SW2**2)

	gm2_2L = gm2_2L + (CA*Sf10435*Sf8520*Ss78*
     &      (1/120.D0*Sf10579 - 1/120.D0*(Sf10444*Sf10585)/Sf10478))/
     &    (CB*CW**3*SB*Sf10444*Ss49**2*SW2)

	gm2_2L = gm2_2L + ((3/80.D0*(Sf4826*Sf4827)/(Sf4820**2*TB) - 
     &        3/40.D0*(MW2*Sf4828*Sf4840*MSbot2(1))/
     &          (Sf4831*Ss241*Ss242))*MStop2(2)*USbot11*
     &      UStop21)/(MW2**2*SW2*MSq2Diff(1,2))

	gm2_2L = gm2_2L + 1/720.D0*
     &    (CA*Sf7509*Sf8897*Ss78*
     &       ((-12*Sf8985*Ss885*Ss906 + 
     &            Sf9074*(6*Ss885 + 6*Ss906) + 
     &            Sf8985*(-Ss887 - Ss908))/(Sf8946*Sf8970) + 
     &         (Sf9065*(-6*Ss885 - 6*Ss906) + 
     &            12*Sf8994*Ss885*Ss906 + Sf8994*(Ss887 + Ss908))/
     &          (Sf8904*Sf8920))*MSbot2(1)*MSbot2(2))/
     &     (CB2*CW**3*MZ2*Sf8929*Sf8930*SW2*MHiggs2(2)*MSbot2(3))

	gm2_2L = gm2_2L + (MW2*Sf7497*Sf7701**2*
     &      (1/60.D0*(Sf7704*MStau2(1)*MStau2(2))/MZ2**3 - 
     &        1/120.D0*(12*Sf7708*Sf7710 + Sf7712*MStau2(4))/MZ2**4))
     &     /(CW2**2*SW2*MStau2(3)**2)

	gm2_2L = gm2_2L + 1/240.D0*
     &    ((-(Sf2814*MSbot2(1)) - Sf2861*MStop2(1))/Sf2320**3 + 
     &       (-(Sf2966*MSbot2(2)) - Sf2999*MStop2(1))/Sf2334**3)/
     &     (MW2**2*SW2*TB2)

	gm2_2L = gm2_2L + (CA*Sf7521*Ss78*Ss905*
     &      (1/96.D0/Ss903 - 
     &        1/720.D0*(-3 + Ss484 - 6*Ss906 + 
     &            Ss474*(-6 + 12*Ss906) + Ss908)/Ss49)*MSbot2(2))/
     &    (CB2*CW**3*Sf26*SW2)

	gm2_2L = gm2_2L + (CA*Sf7501*Ss78*Ss884*
     &      (1/96.D0/Ss878 - 
     &        1/720.D0*(-3 + Ss484 - 6*Ss885 + 
     &            Ss474*(-6 + 12*Ss885) + Ss887)/Ss49)*MSbot2(1))/
     &    (CB2*CW**3*Sf11*SW2)

	gm2_2L = gm2_2L - 1/1440.D0*
     &    (CA*Sf7509*Sf8897*Ss78*
     &       ((Sf8900*Sf8908*MSbot2(1) + Sf8957*Sf8987*MSbot2(2))/
     &          (MZ2*Sf8904) + 
     &         (3*Sf9036 - Sf8916*Sf8928*MSbot2(1) + 
     &            Sf8966*Sf8996*MSbot2(2))/(Sf8920*MHiggs2(2))))/
     &     (CB2*CW**3*Ss49**2*SW2)

	tmp122 = Sf7*MSneu2(1) + Sf153*MStau2(1) + Sf1183*MStau2(2)

	gm2_2L = gm2_2L + (3/160.D0*
     &       (MHiggs2(4)**2*
     &          (-(TB*tmp122) + Sf126*MSbot2(1) + 
     &            Sf139*MSbot2(2) + Sf2319*MStop2(1) + 
     &            Sf4899*MStop2(2)))/CW2 - 
     &      3/80.D0*(Sf4828*(Sf5012*MSbot2(1) + Sf4990*MStop2(2)))/
     &        (Sf4820**2*TB))/(MW2*Ss25**3*SW2)

	gm2_2L = gm2_2L + (SA*Sf7713*Ss78*Ss958*
     &      (3/32.D0 + (Sf61*
     &           (-(3/160.D0*Ss959/MHiggs2(1)) + 
     &             (-(3/40.D0*Ss961) - 3/20.D0*Ss962)/MHiggs2(1)))/Ss44
     &        )*MStau2(2))/(CB2*CW**3*Sf61*Ss956*SW2)

	gm2_2L = gm2_2L + (Sf7701*Sf9842*Ss78*
     &      (SA*Sf10055*(-(3/40.D0*Sf10191/MHiggs2(1)**2) + 
     &           3/40.D0*Sf10192/MHiggs2(1)**2) + 
     &        CA*Sf9312*(3/40.D0*Sf9839/MHiggs2(2)**2 - 
     &           3/40.D0*Sf9843/MHiggs2(2)**2)))/
     &    (CB2*CW**3*MZ2**2*SW2*MStau2(3))

	gm2_2L = gm2_2L + (CA*Sf7990*
     &      (1/48.D0*(Sf10223*Ss78)/Sf82 + 
     &        1/360.D0*(Ss78*Ss968*(-3 - 6*Ss269 - 6*Ss969) + 
     &            Sf10223*(Ss274*Ss46*Ss47 + 
     &               Ss46*Ss47*(12*Ss269*Ss969 + Ss971)))/Ss49)*
     &      MStop2(1))/(CB*CW**3*SB*Ss964*SW2)

	gm2_2L = gm2_2L + 3/80.D0*
     &    (Sf2198*((Sf2205*MSbot2(1) - Sf2211*MStop2(1))/Ss25**3 + 
     &         (Sf2208*MSbot2(1) - Sf2213*MStop2(1))/Ss25**2)*
     &       USbot11*UStop11)/(Sf2206*SW2)

	gm2_2L = gm2_2L + 1/1440.D0*
     &    (SA*Sf7509*Sf9108*Ss78*
     &       ((Sf8900*Sf9112*MSbot2(1) + Sf8957*Sf9174*MSbot2(2))/
     &          (MZ2*Sf8904) + 
     &         (3*Sf9241 - Sf9118*Sf9130*MSbot2(1) + 
     &            Sf9156*Sf9181*MSbot2(2))/(Sf9122*MHiggs2(1))))/
     &     (CB2*CW**3*Ss44**2*SW2)

	gm2_2L = gm2_2L + (Sf7990*
     &      (-(1/360.D0*(CA*(Ss78*Ss968*(-3 - 6*Ss474 - 6*Ss969) + 
     &                Sf10223*Ss46*Ss47*
     &                 (Ss484 + 12*Ss474*Ss969 + Ss971)))/Ss49) + 
     &        1/360.D0*(SA*Ss46*Ss47*
     &            (Ss484 + 12*Ss474*Ss969 + Ss971)*Ss981)/Ss44)*
     &      MStop2(1))/(CB*CW**3*SB*Sf82*SW2)

	gm2_2L = gm2_2L + 29/240.D0*
     &    (MSneu2(1) + 3*MSbot2(1)*USbot11**2 + 
     &       3*MSbot2(2)*USbot21**2 + MStau2(1)*UStau11**2 + 
     &       MStau2(2)*UStau21**2 + 3*MStop2(1)*UStop11**2 + 
     &       3*MStop2(2)*UStop21**2)/(MW2*SW2)

	gm2_2L = gm2_2L + (Sf8532*
     &      (-(1/360.D0*(CA*(Ss78*Ss987*(-3 - 6*Ss474 - 6*Ss988) + 
     &                Sf10243*Ss46*Ss47*
     &                 (Ss484 + 12*Ss474*Ss988 + Ss990)))/Ss49) - 
     &        1/48.D0*(SA*Ss78*Ss997)/Ss995)*MStop2(2))/
     &    (CB*CW**3*SB*Sf113*SW2)

	gm2_2L = gm2_2L + (SA*Sf7579*Ss78*Ss940*
     &      (3/32.D0 + (Sf39*
     &           (-(3/160.D0*Ss941/MHiggs2(1)) + 
     &             (-(3/40.D0*Ss943) - 3/20.D0*Ss944)/MHiggs2(1)))/Ss44
     &        )*MStau2(1))/(CB2*CW**3*Sf39*Ss935*SW2)

	gm2_2L = gm2_2L + 1/80.D0*
     &    (SA*Sf10055*Sf7701*Ss78*
     &       ((-12*Sf10118*Ss928*Ss949 + 
     &            Sf10128*(6*Ss928 + 6*Ss949) + 
     &            Sf10118*(-Ss930 - Ss951))/(Sf10068*Sf9317) + 
     &         (Sf10157*(-6*Ss928 - 6*Ss949) + 
     &            12*Sf10110*Ss928*Ss949 + Sf10110*(Ss930 + Ss951))
     &           /(Sf10098*Sf9344))*MStau2(1)*MStau2(2))/
     &     (CB2*CW**3*MZ2*Sf9131*Sf9132*SW2*MHiggs2(1)*MStau2(3))

	gm2_2L = gm2_2L + (SA*Sf7509*Sf9108*Ss78*
     &      (-(1/240.D0*Sf9272) + 1/240.D0*(Sf9122*Sf9278)/Sf9160))/
     &    (CB2*CW**3*Sf9122*Ss44**2*SW2)

	gm2_2L = gm2_2L + (1/160.D0*
     &       (-(Sf122*MSbot2(1)) - Sf138*MSbot2(2) + 
     &          Sf4*MSneu2(1) + Sf151*MStau2(1) + 
     &          Sf1182*MStau2(2) - Sf2315*MStop2(1) - 
     &          Sf4898*MStop2(2))/MW2**2 - 
     &      3/160.D0*((-3*MW2 + 4*MHiggs2(4))*
     &          (-(TB*tmp122) + Sf126*MSbot2(1) + 
     &            Sf139*MSbot2(2) + Sf2319*MStop2(1) + 
     &            Sf4899*MStop2(2)))/Ss25**3)/(CW2*SW2)

	gm2_2L = gm2_2L + 1/60.D0*
     &    (Sf37*UStau11**2 + Sf50*UStau21**2 + 
     &       3*(Sf72*USbot11**2 + Sf80*USbot21**2)*
     &        UStop11**2 + 
     &       3*(Sf94*USbot11**2 + Sf102*USbot21**2)*
     &        UStop21**2)/(MW2**2*MZ2*SW2**2)

	gm2_2L = gm2_2L + 3/80.D0*
     &    (Sf4828*((Sf4830*MSbot2(1) - Sf4836*MStop2(2))/Ss25**3 + 
     &         (Sf4834*MSbot2(1) - Sf4838*MStop2(2))/Ss25**2)*
     &       USbot11*UStop21)/(Sf4831*SW2)

	gm2_2L = gm2_2L + 3/80.D0*
     &    (USbot21*((Sf2263*
     &            (Sf2272*MSbot2(2) - Sf2276*MStop2(1))*UStop11)
     &           /Sf2270 + 
     &         (Sf4864*(Sf4869*MSbot2(2) - Sf4873*MStop2(2))*
     &            UStop21)/Sf4867))/(Ss25**2*SW2)

	gm2_2L = gm2_2L + 1/80.D0*
     &    Re((Sf1*Sf10*Log(-(1/2.D0*Sf8/MSneu2(1))))/Sf8)/
     &     (MZ2**2*SW2**2)

	gm2_2L = gm2_2L + (Sf7509*Ss78*
     &      (1/480.D0*(SA*Sf9108*MSbot2(3)*(4*Sf7516 - MSbot2(4)))/
     &          (MZ2*MHiggs2(1)) + 
     &        1/960.D0*(CA*Sf8897*
     &            (Sf7516*(-4*Sf8949 + 4*Sf8952) + 
     &              Sf8895*MSbot2(4) - Sf8899*MSbot2(4)))/
     &          (MZ2**2*MHiggs2(2))))/(CB2*CW**3*SW2*MSbot2(3))

	gm2_2L = gm2_2L + 1/80.D0*
     &    (SA*((Sf7579*Ss940*MStau2(1)*
     &            (6*Sf9298*Ss78*Ss928 + 
     &              Sf9297*Ss46*Ss47*(6*Ss928**2 + Ss930) - 
     &              33*Sf9296*Ss78*MStau2(1)))/(Sf39*Ss935) + 
     &         (Sf7713*Ss958*MStau2(2)*
     &            (6*Sf9310*Ss78*Ss949 + 
     &              Sf9309*Ss46*Ss47*(6*Ss949**2 + Ss951) - 
     &              33*Sf9308*Ss78*MStau2(2)))/(Sf61*Ss956)))/
     &     (CB2*CW**3*MZ2*SW2*MHiggs2(1))

	gm2_2L = gm2_2L + (SA*Sf10333*Sf8520*Ss78*
     &      (1/120.D0*Sf10429 - 1/120.D0*(Sf10296*Sf10434)/Sf10286))/
     &    (CB*CW**3*SB*Sf10296*Ss44**2*SW2)

	gm2_2L = gm2_2L - 1/720.D0*
     &    (CA*Sf10435*Sf8520*Ss78*
     &       ((Sf10436*Sf10440*MStop2(1) + 
     &            Sf10597*(3 + 6*Ss474 + 6*Ss969)*MStop2(1) + 
     &            Sf10470*Sf10489*MStop2(2))/(MZ2*Sf10267) + 
     &         (3*Sf10536 - Sf10441*Sf10450*MStop2(1) + 
     &            Sf10475*Sf10496*MStop2(2))/(Sf10444*MHiggs2(2))))
     &      /(CB*CW**3*SB*Ss49**2*SW2)

	gm2_2L = gm2_2L + (Sf7701*Ss78*
     &      (1/160.D0*(SA*Sf10055*
     &            (Sf10059*Sf9313*MStau2(1) + 
     &              Sf10111*Sf9349*MStau2(2)))/Ss44**2 - 
     &        1/160.D0*(CA*Sf9312*
     &            (Sf9313*Sf9321*MStau2(1) + 
     &              Sf9349*Sf9370*MStau2(2)))/Ss49**2))/
     &    (CB2*CW**3*MZ2*Sf9317*SW2)

	gm2_2L = gm2_2L + (SA*(-(1/120.D0*
     &           (Sf7990*Ss981*
     &              (Sf10237*Ss78 - 4*Ss46*Ss47*Ss985*MStop2(1)))/
     &            (Ss976*MHiggs2(1))) + 
     &        1/360.D0*(Sf8532*Ss46*Ss47*
     &            (Ss484 + 12*Ss474*Ss988 + Ss990)*Ss997*MStop2(2))
     &           /Sf113))/(CB*CW**3*SB*Ss44*SW2)

	gm2_2L = gm2_2L + (CA*Sf7713*Ss78*Ss948*
     &      (3/32.D0/Ss946 - 
     &        1/80.D0*(-3 + Ss484 - 6*Ss949 + 
     &            Ss474*(-6 + 12*Ss949) + Ss951)/Ss49)*MStau2(2))/
     &    (CB2*CW**3*Sf61*SW2)

	gm2_2L = gm2_2L + (CA*Sf7701*Sf9312*Ss78*
     &      (3/80.D0*Sf9854 - 3/80.D0*(Sf9331*Sf9859)/Sf9362))/
     &    (CB2*CW**3*Sf9331*Ss49**2*SW2)

	gm2_2L = gm2_2L + (CA*Sf7579*Ss78*Ss927*
     &      (3/32.D0/Ss921 - 
     &        1/80.D0*(-3 + Ss484 - 6*Ss928 + 
     &            Ss474*(-6 + 12*Ss928) + Ss930)/Ss49)*MStau2(1))/
     &    (CB2*CW**3*Sf39*SW2)

	gm2_2L = gm2_2L + 1/1440.D0*
     &    (CA*Ss78*((Sf7501*Ss884*
     &            (3*Ss888 + 12*Ss889 + 
     &              (-6 + 2*Ss274 - 12*Ss885 + 
     &                 Ss269*(-12 + 24*Ss885) + 2*Ss887)*MSbot2(1))
     &            )/Ss878 + 
     &         (Sf7521*Ss905*
     &            (3*Ss909 + 12*Ss910 + 
     &              (-6 + 2*Ss274 - 12*Ss906 + 
     &                 Ss269*(-12 + 24*Ss906) + 2*Ss908)*MSbot2(2))
     &            )/Ss903))/(CB2*CW**3*Ss49*SW2)

	gm2_2L = gm2_2L - (Sf7701*Ss78*
     &      (-(3/160.D0*(SA*Sf10055)/MHiggs2(1)) + 
     &        3/160.D0*(CA*Sf9312)/MHiggs2(2))*(4*Sf7708 - MStau2(4))
     &      )/(CB2*CW**3*MZ2*SW2)

	gm2_2L = gm2_2L + (SA*Sf7990*Ss981*
     &      (1/120.D0*(Sf10241*Ss78 - 4*Sf10230*Ss46*Ss47*MStop2(1))/
     &          Ss44 + 1/360.D0*
     &         (MStop2(1)*(-6*Sf10234*Ss78*Ss969 - 
     &              Sf10233*Ss46*Ss47*(6*Ss969**2 + Ss971) + 
     &              33*Sf10231*Ss78*MStop2(1)))/(Ss976*MHiggs2(1)))
     &      )/(CB*CW**3*MZ2*SB*Sf82*SW2)

	gm2_2L = gm2_2L - 1/1440.D0*
     &    (SA*Sf10333*Sf8520*Ss78*
     &       ((3*Sf10400*Sf10408 + 2*Sf10351*Sf10389*MStop2(1) + 
     &            6*Sf10361*Sf9222*MStop2(1) - 
     &            6*Sf10364*Sf10370*MStop2(2) - 
     &            2*Sf10378*Sf10382*MStop2(2))/Sf10286 + 
     &         (3*Sf10397*Sf10400 - 6*Sf10334*Sf10340*MStop2(1) - 
     &            2*Sf10351*Sf10355*MStop2(1) + 
     &            2*Sf10378*Sf10387*MStop2(2) + 
     &            6*Sf10341*Sf9222*MStop2(2))/Sf10296))/
     &     (CB*CW**3*SB*Ss44**2*SW2*MHiggs2(1))

	gm2_2L = gm2_2L + 1/720.D0*
     &    (SA*Sf7509*Sf9108*Ss78*
     &       ((-12*Sf9179*Ss885*Ss906 + 
     &            Sf9193*(6*Ss885 + 6*Ss906) + 
     &            Sf9179*(-Ss887 - Ss908))/(Sf8904*Sf9122) + 
     &         (Sf9221*(-6*Ss885 - 6*Ss906) + 
     &            12*Sf9173*Ss885*Ss906 + Sf9173*(Ss887 + Ss908))/
     &          (Sf8946*Sf9160))*MSbot2(1)*MSbot2(2))/
     &     (CB2*CW**3*MZ2*Sf9131*Sf9132*SW2*MHiggs2(1)*MSbot2(3))

	gm2_2L = gm2_2L + (CA*Sf10435*Sf8520*Ss78*
     &      (1/120.D0*Sf10455/Sf10444 - 
     &        1/120.D0*(Sf10267*Sf10518)/(Sf10276*Sf10478))*
     &      MStop2(1)*MStop2(2))/
     &    (CB*CW**3*MZ2*SB*Sf10267*Sf8929*Sf8930*SW2*MStop2(3))

	gm2_2L = gm2_2L + 1/1440.D0*
     &    (SA*Sf7521*((Ss78*Ss915*
     &            (3*Ss916 + 12*Ss918 + 
     &              (-6 + 2*Ss268 - 12*Ss906 + 
     &                 Ss263*(-12 + 24*Ss906) + 2*Ss908)*MSbot2(2))
     &            )/Ss913 + 
     &         (-3*(Sf7523 + 4*Sf7526)*Ss78*Ss915 + 
     &            2*(Sf9289*Ss46*Ss47*
     &                (Ss484 + 12*Ss474*Ss906 + Ss908) + 
     &               3*Ss78*Ss915 + 6*Ss474*Ss78*Ss915 + 
     &               6*Ss78*Ss906*Ss915)*MSbot2(2))/Sf26))/
     &     (CB2*CW**3*Ss44*SW2)

	gm2_2L = gm2_2L + 1/1440.D0*
     &    (SA*Sf7501*((Ss78*Ss897*
     &            (3*Ss898 + 12*Ss900 + 
     &              (-6 + 2*Ss268 - 12*Ss885 + 
     &                 Ss263*(-12 + 24*Ss885) + 2*Ss887)*MSbot2(1))
     &            )/Ss892 + 
     &         (-3*(Sf7504 + 4*Sf7507)*Ss78*Ss897 + 
     &            2*(Sf8884*Ss46*Ss47*
     &                (Ss484 + 12*Ss474*Ss885 + Ss887) + 
     &               3*Ss78*Ss897 + 6*Ss474*Ss78*Ss897 + 
     &               6*Ss78*Ss885*Ss897)*MSbot2(1))/Sf11))/
     &     (CB2*CW**3*Ss44*SW2)

	gm2_2L = gm2_2L + 1/160.D0*
     &    (CA*Ss78*((Sf7579*Ss927*
     &            (3*Ss931 + 12*Ss932 + 
     &              (-6 + 2*Ss274 - 12*Ss928 + 
     &                 Ss269*(-12 + 24*Ss928) + 2*Ss930)*MStau2(1))
     &            )/Ss921 + 
     &         (Sf7713*Ss948*
     &            (3*Ss952 + 12*Ss953 + 
     &              (-6 + 2*Ss274 - 12*Ss949 + 
     &                 Ss269*(-12 + 24*Ss949) + 2*Ss951)*MStau2(2))
     &            )/Ss946))/(CB2*CW**3*Ss49*SW2)

	gm2_2L = gm2_2L + (SA*Sf10055*Sf7701*Ss78*
     &      (-(3/80.D0*Sf10204) + 3/80.D0*(Sf10068*Sf10210)/Sf10098))/
     &    (CB2*CW**3*Sf10068*Ss44**2*SW2)

	gm2_2L = gm2_2L + (1/4320.D0*
     &       (-Sf8872 + 12*Sf8503*Log(MSbot2(1)) + 
     &          12*Sf8519*Log(MSbot2(2)) + 
     &          36*Sf7788*Log(MStau2(1)) + 
     &          36*Sf7769*Log(MStau2(2)) - 
     &          12*Sf8762*Log(MStop2(1)) - 12*Sf8740*Log(MStop2(2))
     &          )/MZ2**6 - 
     &      1/120.D0*(Sf7643*Log(MSneu2(1))*MSneu2(1))/MZ2**3)/
     &    (MW2**2*SW2**2)

	gm2_2L = gm2_2L + (-(1/15.D0*(Sf8596*Log(MM2))) - 
     &      1/720.D0*(Sf8575*Log(MW2))/SW2**2)/(MW2*MZ2**3)

	gm2_2L = gm2_2L + (CA*Sf7509*Sf8897*Ss78*
     &      (1/240.D0*Sf8938/Sf8920 - 
     &        1/240.D0*(Sf8904*Sf9016)/(Sf8946*Sf8970))*MSbot2(1)*
     &      MSbot2(2))/
     &    (CB2*CW**3*MZ2*Sf8904*Sf8929*Sf8930*SW2*MSbot2(3))

	gm2_2L = gm2_2L + (SA*Sf10055*Sf7701*Ss78*
     &      (-(3/80.D0*Sf10079/Sf9317) + 
     &        3/80.D0*(Sf10068*Sf10149)/(Sf10098*Sf9344))*MStau2(1)*
     &      MStau2(2))/
     &    (CB2*CW**3*MZ2*Sf10068*Sf9131*Sf9132*SW2*MStau2(3))

	gm2_2L = gm2_2L + 1/160.D0*
     &    (SA*Sf7713*((Ss78*Ss958*
     &            (3*Ss959 + 12*Ss961 + 
     &              (-6 + 2*Ss268 - 12*Ss949 + 
     &                 Ss263*(-12 + 24*Ss949) + 2*Ss951)*MStau2(2))
     &            )/Ss956 + 
     &         (-3*(Sf7734 + 4*Sf7737)*Ss78*Ss958 + 
     &            2*(Sf10221*Ss46*Ss47*
     &                (Ss484 + 12*Ss474*Ss949 + Ss951) + 
     &               3*Ss78*Ss958 + 6*Ss474*Ss78*Ss958 + 
     &               6*Ss78*Ss949*Ss958)*MStau2(2))/Sf61))/
     &     (CB2*CW**3*Ss44*SW2)

	gm2_2L = gm2_2L + 1/160.D0*
     &    (SA*Sf7579*((Ss78*Ss940*
     &            (3*Ss941 + 12*Ss943 + 
     &              (-6 + 2*Ss268 - 12*Ss928 + 
     &                 Ss263*(-12 + 24*Ss928) + 2*Ss930)*MStau2(1))
     &            )/Ss935 + 
     &         (-3*(Sf7581 + 4*Sf7584)*Ss78*Ss940 + 
     &            2*(Sf9301*Ss46*Ss47*
     &                (Ss484 + 12*Ss474*Ss928 + Ss930) + 
     &               3*Ss78*Ss940 + 6*Ss474*Ss78*Ss940 + 
     &               6*Ss78*Ss928*Ss940)*MStau2(1))/Sf39))/
     &     (CB2*CW**3*Ss44*SW2)

	gm2_2L = gm2_2L + (SA*Sf7990*Ss981*
     &      (-(1/360.D0*(Ss268*Ss46*Ss47)/Ss44) + 
     &        (-(1/48.D0*(Ss44*Ss78)/Sf82) + 
     &           Ss46*Ss47*(-(1/30.D0*(Ss263*Ss969)) - 1/360.D0*Ss971))
     &          /Ss44)*MStop2(1))/(CB*CW**3*SB*Ss976*SW2)

	gm2_2L = gm2_2L + (SA*Sf10333*Sf8520*Ss78*
     &      (1/120.D0*Sf10360/Sf10296 - 
     &        1/120.D0*(Sf10267*Sf10413)/(Sf10276*Sf10286))*
     &      MStop2(1)*MStop2(2))/
     &    (CB*CW**3*MZ2*SB*Sf10267*Sf9131*Sf9132*SW2*MStop2(3))

	gm2_2L = gm2_2L + 19/135.D0*(Sf7990*Ss78) + 19/135.D0*(Sf8532*Ss78)

	gm2_2L = gm2_2L + 19/270.D0*(Sf7501*Ss78) + 19/270.D0*(Sf7521*Ss78)

	gm2_2L = gm2_2L + (SA*Sf7509*Sf9108*Ss78*
     &      (-(1/240.D0*Sf9140/Sf9122) + 
     &        1/240.D0*(Sf8904*Sf9214)/(Sf8946*Sf9160))*MSbot2(1)*
     &      MSbot2(2))/
     &    (CB2*CW**3*MZ2*Sf8904*Sf9131*Sf9132*SW2*MSbot2(3))

	gm2_2L = gm2_2L + 19/90.D0*(Sf7579*Ss78) + 19/90.D0*(Sf7713*Ss78)

	end

