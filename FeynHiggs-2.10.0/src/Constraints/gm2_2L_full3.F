* gm2_2L_full3.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full3
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp5, tmp6, tmp7, tmp8

	Sf4529 = 79*MSbot2(1) + 71*MStop2(1)

	Sf6594 = 71*MSbot2(1) + 79*MStop2(2)

	Sf6468 = 71*MSbot2(1) + 79*MStop2(1)

	Sf8293 = 67*MSbot2(1) + 64*MStop2(2)

	Sf7999 = 67*MSbot2(1) + 40*MStop2(2)

	Sf8269 = 67*MSbot2(1) + 64*MStop2(1)

	Sf7791 = 67*MSbot2(1) + 40*MStop2(1)

	Sf8443 = 62*MSbot2(1) + 345*MStop2(2)

	Sf8419 = 62*MSbot2(1) + 345*MStop2(1)

	Sf8064 = 61*MSbot2(1) + 55*MStop2(2)

	Sf8241 = 61*MSbot2(1) + 39*MStop2(2)

	Sf7855 = 61*MSbot2(1) + 55*MStop2(1)

	Sf8229 = 61*MSbot2(1) + 39*MStop2(1)

	Sf7041 = 57*MSbot2(1) + 43*MStop2(2)

	Sf4534 = 57*MSbot2(1) + 43*MStop2(1)

	Sf8155 = 53*MSbot2(1) + 47*MStop2(2)

	Sf7946 = 53*MSbot2(1) + 47*MStop2(1)

	Sf6600 = 49*MSbot2(1) + 51*MStop2(2)

	Sf6474 = 49*MSbot2(1) + 51*MStop2(1)

	Sf8001 = 47*MSbot2(1) + 101*MStop2(2)

	Sf7793 = 47*MSbot2(1) + 101*MStop2(1)

	Sf6632 = 45*MSbot2(1) + 6*MStop2(2)

	Sf6514 = 45*MSbot2(1) + 6*MStop2(1)

	Sf8044 = 44*MSbot2(1) - 64*MStop2(2)

	Sf7835 = 44*MSbot2(1) - 64*MStop2(1)

	Sf6524 = 43*MSbot2(1) + 31*MStop2(2)

	Sf4275 = 43*MSbot2(1) + 31*MStop2(1)

	Sf7997 = 37*MSbot2(1) + 55*MStop2(2)

	Sf6637 = 37*MSbot2(1) + 43*MStop2(2)

	Sf8023 = 37*MSbot2(1) + 7*MStop2(2)

	Sf7789 = 37*MSbot2(1) + 55*MStop2(1)

	Sf4578 = 37*MSbot2(1) + 43*MStop2(1)

	Sf7814 = 37*MSbot2(1) + 7*MStop2(1)

	Sf7042 = 36*MSbot2(1) + 75*MStop2(2)

	Sf4535 = 36*MSbot2(1) + 75*MStop2(1)

	Sf5255 = 35*MSbot2(1) + 26*MStop2(2)

	Sf2799 = 35*MSbot2(1) + 26*MStop2(1)

	Sf8154 = 34*MSbot2(1) + 25*MStop2(2)

	Sf7945 = 34*MSbot2(1) + 25*MStop2(1)

	Sf6593 = 31*MSbot2(1) - 4*MStop2(2)

	Sf6467 = 31*MSbot2(1) - 4*MStop2(1)

	Sf5789 = 27*MSbot2(1) + 29*MStop2(2)

	Sf3794 = 27*MSbot2(1) + 29*MStop2(1)

	Sf5291 = 26*MSbot2(1) + 35*MStop2(2)

	Sf2851 = 26*MSbot2(1) + 35*MStop2(1)

	Sf8167 = 24*MSbot2(1) - 26*MStop2(2)

	Sf7958 = 24*MSbot2(1) - 26*MStop2(1)

	Sf5039 = 23*MSbot2(1) + 8*MStop2(2)

	Sf2574 = 23*MSbot2(1) + 8*MStop2(1)

	Sf8685 = 22*MSbot2(1) + 25*MStop2(2)

	Sf8645 = 22*MSbot2(1) + 25*MStop2(1)

	Sf8439 = 20*MSbot2(1) - 95*MStop2(2)

	Sf8643 = 20*MSbot2(1) - 95*MStop2(1)

	Sf5816 = 19*MSbot2(1) + 15*MStop2(2)

	Sf6599 = 19*MSbot2(1) + 8*MStop2(2)

	Sf3994 = 19*MSbot2(1) + 15*MStop2(1)

	Sf6473 = 19*MSbot2(1) + 8*MStop2(1)

	Sf5302 = 17*MSbot2(1) + 15*MStop2(2)

	Sf2862 = 17*MSbot2(1) + 15*MStop2(1)

	Sf5615 = 15*MSbot2(1) + 19*MStop2(2)

	Sf5475 = 15*MSbot2(1) + 14*MStop2(2)

	Sf6620 = 15*MSbot2(1) + 2*MStop2(2)

	Sf3416 = 15*MSbot2(1) + 19*MStop2(1)

	Sf3078 = 15*MSbot2(1) + 14*MStop2(1)

	Sf6502 = 15*MSbot2(1) + 2*MStop2(1)

	Sf8060 = 13*MSbot2(1) + 7*MStop2(2)

	Sf8240 = 13*MSbot2(1) - 6*MStop2(2)

	Sf7851 = 13*MSbot2(1) + 7*MStop2(1)

	Sf8228 = 13*MSbot2(1) - 6*MStop2(1)

	Sf5539 = 11*MSbot2(1) + 10*MStop2(2)

	Sf5448 = 11*MSbot2(1) + 9*MStop2(2)

	Sf5812 = 11*MSbot2(1) + 6*MStop2(2)

	Sf5815 = 11*MSbot2(1) + 5*MStop2(2)

	Sf3178 = 11*MSbot2(1) + 10*MStop2(1)

	Sf3059 = 11*MSbot2(1) + 9*MStop2(1)

	Sf5719 = 11*MSbot2(1) + 6*MStop2(1)

	Sf3993 = 11*MSbot2(1) + 5*MStop2(1)

	Sf6984 = 10*MSbot2(1) + 9*MStop2(2)

	Sf4271 = 10*MSbot2(1) + 9*MStop2(1)

	Sf4959 = 9*MSbot2(1) + 13*MStop2(2)

	Sf5238 = 9*MSbot2(1) + 11*MStop2(2)

	Sf6523 = 9*MSbot2(1) + 7*MStop2(2)

	Sf5811 = 9*MSbot2(1) + 2*MStop2(2)

	Sf2431 = 9*MSbot2(1) + 13*MStop2(1)

	Sf2769 = 9*MSbot2(1) + 11*MStop2(1)

	Sf4274 = 9*MSbot2(1) + 7*MStop2(1)

	Sf5718 = 9*MSbot2(1) + 2*MStop2(1)

	Sf8062 = 8*MSbot2(1) + 26*MStop2(2)

	Sf6295 = 8*MSbot2(1) + 9*MStop2(2)

	Sf8245 = 8*MSbot2(1) - 6*MStop2(2)

	Sf7853 = 8*MSbot2(1) + 26*MStop2(1)

	Sf3997 = 8*MSbot2(1) + 9*MStop2(1)

	Sf8233 = 8*MSbot2(1) - 6*MStop2(1)

	Sf5273 = 7*MSbot2(1) + 5*MStop2(2)

	Sf4900 = 7*MSbot2(1) - 4*MStop2(2)

	Sf2821 = 7*MSbot2(1) + 5*MStop2(1)

	Sf2321 = 7*MSbot2(1) - 4*MStop2(1)

	Sf7071 = 6*MSbot2(1) + 45*MStop2(2)

	Sf6950 = 6*MSbot2(1) + 4*MStop2(2)

	Sf4573 = 6*MSbot2(1) + 45*MStop2(1)

	Sf4232 = 6*MSbot2(1) + 4*MStop2(1)

	Sf6520 = 5*MSbot2(1) + 8*MStop2(2)

	Sf4901 = 5*MSbot2(1) + 4*MStop2(2)

	Sf4853 = 5*MSbot2(1) + MStop2(2)

	Sf6420 = 5*MSbot2(1) + 8*MStop2(1)

	Sf2322 = 5*MSbot2(1) + 4*MStop2(1)

	Sf2247 = 5*MSbot2(1) + MStop2(1)

	Sf4101 = 5*MSbot2(1) - 18*MStop2(1)

	Sf6983 = 4*MSbot2(1) + 3*MStop2(2)

	Sf5254 = 4*MSbot2(1) + MStop2(2)

	Sf4845 = 4*MSbot2(1) - 2*MStop2(2)

	Sf8683 = 4*MSbot2(1) - 19*MStop2(2)

	Sf7035 = 4*MSbot2(1) - 31*MStop2(2)

	Sf4270 = 4*MSbot2(1) + 3*MStop2(1)

	Sf2798 = 4*MSbot2(1) + MStop2(1)

	Sf8415 = 4*MSbot2(1) - 19*MStop2(1)

	Sf4528 = 4*MSbot2(1) - 31*MStop2(1)

	Sf5614 = 3*MSbot2(1) + 7*MStop2(2)

	Sf6621 = 3*MSbot2(1) + 4*MStop2(2)

	Sf4821 = 3*MSbot2(1) - MStop2(2)

	Sf5538 = 3*MSbot2(1) - 8*MStop2(2)

	Sf3415 = 3*MSbot2(1) + 7*MStop2(1)

	Sf6503 = 3*MSbot2(1) + 4*MStop2(1)

	Sf2328 = 3*MSbot2(1) - MStop2(1)

	Sf3177 = 3*MSbot2(1) - 8*MStop2(1)

	Sf7057 = 2*MSbot2(1) + 15*MStop2(2)

	Sf5021 = 2*MSbot2(1) + 3*MStop2(2)

	Sf4855 = 2*MSbot2(1) + MStop2(2)

	Sf8017 = 2*MSbot2(1) - MStop2(2)

	Sf4561 = 2*MSbot2(1) + 15*MStop2(1)

	Sf2554 = 2*MSbot2(1) + 3*MStop2(1)

	Sf2249 = 2*MSbot2(1) + MStop2(1)

	Sf2217 = 2*MSbot2(1) - MStop2(1)

	Sf8990 = 2*MSbot2(1) + MSbot2(2)

	Sf8992 = MSbot2(1) - MSbot2(3)

	Sf5290 = MSbot2(1) + 4*MStop2(2)

	Sf4839 = MSbot2(1) + 3*MStop2(2)

	Sf4843 = MSbot2(1) + 2*MStop2(2)

	Sf91 = MSbot2(1) + MStop2(2)

	Sf8011 = MSbot2(1) - 5*MStop2(2)

	Sf5193 = MSbot2(1) - 7*MStop2(2)

	Sf5030 = MSbot2(1) - 8*MStop2(2)

	Sf2850 = MSbot2(1) + 4*MStop2(1)

	Sf2207 = MSbot2(1) + 3*MStop2(1)

	Sf2219 = MSbot2(1) + 2*MStop2(1)

	Sf69 = MSbot2(1) + MStop2(1)

	Sf7808 = MSbot2(1) - 5*MStop2(1)

	Sf3736 = MSbot2(1) - 7*MStop2(1)

	Sf2565 = MSbot2(1) - 8*MStop2(1)

	Sf8981 = MSbot2(1) + 2*MSbot2(2)

	Sf8983 = -MSbot2(2) - MSbot2(3)

	Sf22 = MSbot(1) + MSbot(2)

	Sf20 = MSbot(1) - MSbot(2)

	Sf10415 = Log(MStop2(1)) + Log(MStop2(2))

	Sf8814 = Log(MStop2(1)) - Log(MStop2(2))

	Sf9840 = Log(MStau2(1)) + Log(MStau2(2))

	Sf8784 = Log(MStau2(1)) - Log(MStau2(2))

	Sf1519 = Log(MSneu2(1)) + Log(MStau2(2))

	Sf1518 = Log(MSneu2(1)) - Log(MStau2(2))

	Sf624 = Log(MSneu2(1)) + Log(MStau2(1))

	Sf623 = Log(MSneu2(1)) - Log(MStau2(1))

	Sf5164 = Log(MSbot2(2)) + Log(MStop2(2))

	Sf5163 = Log(MSbot2(2)) - Log(MStop2(2))

	Sf2757 = Log(MSbot2(2)) + Log(MStop2(1))

	Sf2756 = Log(MSbot2(2)) - Log(MStop2(1))

	Sf5014 = Log(MSbot2(1)) + Log(MStop2(2))

	Sf5013 = Log(MSbot2(1)) - Log(MStop2(2))

	Sf2581 = Log(MSbot2(1)) + Log(MStop2(1))

	Sf2580 = Log(MSbot2(1)) - Log(MStop2(1))

	Sf9077 = Log(MSbot2(1)) + Log(MSbot2(2))

	Sf8766 = Log(MSbot2(1)) - Log(MSbot2(2))

	Sf2225 = reMUE*TB - Re(Abot)

	Sf2196 = reMUE*TB + Re(Atop)

	Sf215 = reMUE*TB - Re(Atau)

	Sf152 = MW2*S2B*SW2 - CW2*ML2*TB

	Sf125 = MW2*S2B*SW2 - 3*CW2*MBot2*TB

	Sf10269 = MTop*Sf10268*Ss978 + 2*CA*MTop2*UStop12*UStop22

	Sf8700 = 8*CW2*Ss78 - 9*Sf8532*Ss922

	Sf8663 = 8*CW2*Ss78 - 9*Sf7990*Ss922

	Sf7746 = 4*CW2*Ss78 - 9*Sf7713*Ss922

	Sf7724 = 4*CW2*Ss78 - 9*Sf7579*Ss922

	Sf8403 = 4*CW2*Ss78 - 9*Sf7521*Ss922

	Sf8307 = 4*CW2*Ss78 - 9*Sf7501*Ss922

	Sf3008 = 52*Ss30 + 63*C2B*MSbot2(2)

	Sf2876 = 52*Ss30 + 63*C2B*MSbot2(1)

	Sf1748 = CB*Ss25**2 + 3*MW2*SB*Ss115*TB

	Sf1739 = CA*Ss25**2 - 3*MW2*SA*Ss115*TB

	Sf10252 = Ss145 - 4*MStop2(2)

	Sf10232 = Ss145 - 4*MStop2(1)

	Sf9308 = Ss145 - 4*MStau2(2)

	Sf9296 = Ss145 - 4*MStau2(1)

	Sf8891 = Ss145 - 4*MSbot2(2)

	Sf8879 = Ss145 - 4*MSbot2(1)

	Sf10247 = Ss143 - 4*MStop2(2)

	Sf10227 = Ss143 - 4*MStop2(1)

	Sf9304 = Ss143 - 4*MStau2(2)

	Sf9292 = Ss143 - 4*MStau2(1)

	Sf8887 = Ss143 - 4*MSbot2(2)

	Sf8875 = Ss143 - 4*MSbot2(1)

	Sf4201 = Ss115 + 3*MSbot2(2)

	Sf6137 = Ss115 - 3*MSbot2(2)

	Sf4121 = Ss115 + 3*MSbot2(1)

	Sf5754 = Ss115 - 3*MSbot2(1)

	Sf8808 = 21*MZ2*Sf99 + 20*MSq2Diff(2,2)**2

	Sf101 = 7*MZ2*Sf99 + 8*MSq2Diff(2,2)**2

	Sf5203 = 4*MW2*Sf99 + 7*MSq2Diff(2,2)**2

	Sf5199 = 2*MW2*Sf99 + 3*MSq2Diff(2,2)**2

	Sf5202 = MW2*Sf99 + MSq2Diff(2,2)**2

	Sf5084 = 5*Sf99*MHiggs2(4) + 2*MSq2Diff(2,2)**2

	Sf9776 = Sf9577 + Sf9578*MStau2(2)

	Sf9579 = Sf9577 + Sf9578*MStau2(1)

	Sf9781 = Sf9572 - 2*Sf9573*MStau2(2)

	Sf9574 = Sf9572 - 2*Sf9573*MStau2(1)

	Sf9782 = Sf9571 - 2*Sf9781*MStau2(2)

	Sf9575 = Sf9571 - 2*Sf9574*MStau2(1)

	Sf9784 = 4*Sf9566 + Sf9567*MStau2(2)

	Sf9568 = 4*Sf9566 + Sf9567*MStau2(1)

	Sf9785 = Sf9565 + Sf9784*MStau2(2)

	Sf9569 = Sf9565 + Sf9568*MStau2(1)

	Sf9790 = Sf9559 - Sf9560*MStau2(2)

	Sf9561 = Sf9559 - Sf9560*MStau2(1)

	Sf9791 = Sf9558 - 2*Sf9790*MStau2(2)

	Sf9562 = Sf9558 - 2*Sf9561*MStau2(1)

	Sf9787 = Sf9551 - Sf9552*MStau2(2)

	Sf9553 = Sf9551 - Sf9552*MStau2(1)

	Sf9788 = Sf9550 - 2*Sf9787*MStau2(2)

	Sf9554 = Sf9550 - 2*Sf9553*MStau2(1)

	Sf947 = SB*Sf946 + CAB*CB2*Ss443*UStau11

	Sf9352 = Sf9350 - MHiggs2(2)*MStau2(1)

	Sf9320 = Sf9319 - MHiggs2(2)*MStau2(2)

	Sf10364 = Sf9222 + 2*Ss988*MHiggs2(1)

	Sf10334 = Sf9222 + 2*Ss969*MHiggs2(1)

	Sf10162 = Sf9222 + 2*Ss949*MHiggs2(1)

	Sf10176 = Sf9222 + 2*Ss928*MHiggs2(1)

	Sf9227 = Sf9222 + 2*Ss906*MHiggs2(1)

	Sf9242 = Sf9222 + 2*Ss885*MHiggs2(1)

	Sf9150 = Sf9148 - MHiggs2(1)*MSbot2(1)

	Sf10409 = Sf9136 + 2*MStop2(2)**2

	Sf10357 = Sf9136 + 2*MStop2(1)**2

	Sf10145 = Sf9136 + 2*MStau2(2)**2

	Sf10076 = Sf9136 + 2*MStau2(1)**2

	Sf9210 = Sf9136 + 2*MSbot2(2)**2

	Sf9137 = Sf9136 + 2*MSbot2(1)**2

	Sf9111 = Sf9110 - MHiggs2(1)*MSbot2(2)

	Sf8802 = 21*MZ2*Sf91 + 20*MSq2Diff(1,2)**2

	Sf93 = 7*MZ2*Sf91 + 8*MSq2Diff(1,2)**2

	Sf6959 = 27*MW2*Sf91 + 22*MSq2Diff(1,2)**2

	Sf5191 = 4*MW2*Sf91 + 7*MSq2Diff(1,2)**2

	Sf5187 = 2*MW2*Sf91 + 3*MSq2Diff(1,2)**2

	Sf5190 = MW2*Sf91 + MSq2Diff(1,2)**2

	Sf4958 = 5*Sf91*MHiggs2(4) + 2*MSq2Diff(1,2)**2

	Sf933 = Sf906 + 18*C2B*Sf905*MSl2Diff(1,1)

	Sf932 = 115*Sf902 + 27*C2B*Sf899*MSl2Diff(1,1)

	Sf10523 = Sf9017 + 2*Ss988*MHiggs2(2)

	Sf10537 = Sf9017 + 2*Ss969*MHiggs2(2)

	Sf9612 = Sf9017 + 2*Ss949*MHiggs2(2)

	Sf9593 = Sf9017 + 2*Ss928*MHiggs2(2)

	Sf9038 = Sf9017 + 2*Ss906*MHiggs2(2)

	Sf9018 = Sf9017 + 2*Ss885*MHiggs2(2)

	Sf8960 = Sf8958 - MHiggs2(2)*MSbot2(1)

	Sf10514 = Sf8934 + 2*MStop2(2)**2

	Sf10452 = Sf8934 + 2*MStop2(1)**2

	Sf9012 = Sf8934 + 2*MSbot2(2)**2

	Sf8935 = Sf8934 + 2*MSbot2(1)**2

	Sf10466 = Sf8909 + 2*MZ2*Ss988

	Sf10456 = Sf8909 + 2*MZ2*Ss969

	Sf9345 = Sf8909 + 2*MZ2*Ss949

	Sf9322 = Sf8909 + 2*MZ2*Ss928

	Sf8953 = Sf8909 + 2*MZ2*Ss906

	Sf8910 = Sf8909 + 2*MZ2*Ss885

	Sf8907 = Sf8906 - MHiggs2(2)*MSbot2(2)

	Sf8749 = MZ2*Sf8748 + MStop2(3)*UStop12*UStop21*UStop22

	Sf8747 = MZ2*Sf8746 + 
     &   13*MStop2(3)*UStop12*UStop21*UStop22

	Sf8745 = MZ2*Sf8744 + 
     &   19*MStop2(3)*UStop12*UStop21*UStop22

	Sf8609 = MZ2*Sf8606 + 15*Sf8608

	Sf8551 = Sf7314*MSneu2(1)

	Sf8701 = 3*Sf7722*Sf8533 + 2*Sf8532*Sf8700*SW2

	Sf8507 = 3*Sf8490 - 4*MW2*Sf8491*USbot22**2

	Sf840 = 4*Sf821 - 3*C2B*Sf819*MSl2Diff(1,1)**2

	Sf839 = C2B*Sf820 + 2*Sf816*MSl2Diff(1,1)

	Sf843 = Sf817 + 9*C2B*Sf815*MSl2Diff(1,1)

	Sf8141 = MW2*Sf8140 - 59*MSq2Diff(2,2)**2

	Sf8139 = MW2*Sf8138 - 16*MSq2Diff(2,2)**2

	Sf8063 = MW2*Sf8062 - 59*MSq2Diff(1,2)**2

	Sf8061 = MW2*Sf8060 - 16*MSq2Diff(1,2)**2

	Sf8664 = 3*Sf7722*Sf7991 + 2*Sf7990*Sf8663*SW2

	Sf7932 = MW2*Sf7931 - 59*MSq2Diff(2,1)**2

	Sf7930 = MW2*Sf7929 - 16*MSq2Diff(2,1)**2

	Sf7854 = MW2*Sf7853 - 59*MSq2Diff(1,1)**2

	Sf7852 = MW2*Sf7851 - 16*MSq2Diff(1,1)**2

	Sf7778 = MZ2*Sf7777 + MStau2(3)*UStau12*UStau21*UStau22

	Sf7776 = MZ2*Sf7775 + 
     &   7*MStau2(3)*UStau12*UStau21*UStau22

	Sf7774 = MZ2*Sf7773 + 
     &   23*MStau2(3)*UStau12*UStau21*UStau22

	Sf7772 = MZ2*Sf7771 + 
     &   71*MStau2(3)*UStau12*UStau21*UStau22

	Sf7747 = 3*Sf7714*Sf7722 + 2*Sf7713*Sf7746*SW2

	Sf7725 = 3*Sf7580*Sf7722 + 2*Sf7579*Sf7724*SW2

	Sf8404 = 3*Sf7522*Sf7722 + 2*Sf7521*Sf8403*SW2

	Sf8308 = 3*Sf7503*Sf7722 + 2*Sf7501*Sf8307*SW2

	Sf8795 = 21*MZ2*Sf77 + 20*MSq2Diff(2,1)**2

	Sf79 = 7*MZ2*Sf77 + 8*MSq2Diff(2,1)**2

	Sf5183 = 4*MW2*Sf77 + 7*MSq2Diff(2,1)**2

	Sf3743 = 2*MW2*Sf77 + 3*MSq2Diff(2,1)**2

	Sf5182 = MW2*Sf77 + MSq2Diff(2,1)**2

	Sf2644 = 5*Sf77*MHiggs2(4) + 2*MSq2Diff(2,1)**2

	Sf7629 = 3*MW2*Sf7614 - 25*MSl2Diff(2,1)**2

	Sf8716 = Sf7579 + Sf7713

	Sf7571 = 3*MW2*Sf7556 - 25*MSl2Diff(1,1)**2

	Sf7397 = 4*MW2 - MZ2*Sf7396

	Sf7353 = 4*MW2 - MZ2*Sf7352

	Sf10220 = ML*SA + Ss937*UStau21*UStau22

	Sf7324 = 2*MW2 - MZ2*Sf7323

	Sf9300 = ML*SA + Ss937*UStau11*UStau12

	Sf7309 = 2*MW2 - MZ2*Sf7308

	Sf9288 = MBot*SA + Ss894*USbot21*USbot22

	Sf7306 = 2*MW2 + MZ2*Sf7304

	Sf8883 = MBot*SA + Ss894*USbot11*USbot12

	Sf7300 = 2*MW2 + MZ2*Sf7298

	Sf7287 = 28 - 6*Sf7286

	Sf7291 = 4 + 32*SW2

	Sf7489 = 45 - 32*UStop22**2

	Sf7486 = 45 - 4*Sf7485

	Sf6955 = MW2*Sf6954 - MSq2Diff(2,2)**2

	Sf6951 = MW2*Sf6950 - MSq2Diff(1,2)**2

	Sf8789 = 21*MZ2*Sf69 + 20*MSq2Diff(1,1)**2

	Sf71 = 7*MZ2*Sf69 + 8*MSq2Diff(1,1)**2

	Sf5177 = 4*MW2*Sf69 + 7*MSq2Diff(1,1)**2

	Sf3734 = 2*MW2*Sf69 + 3*MSq2Diff(1,1)**2

	Sf5176 = MW2*Sf69 + MSq2Diff(1,1)**2

	Sf2430 = 5*Sf69*MHiggs2(4) + 2*MSq2Diff(1,1)**2

	Sf6773 = MW2*Sf5381 + 3*MSq2Diff(2,2)**2

	Sf6518 = MW2*Sf5273 + 3*MSq2Diff(1,2)**2

	Sf5229 = -3*CB*MTop*SA*Sf5225 + CA*SB*Sf5228

	Sf5087 = -2*Sf5086 + Sf5085*MHiggs2(4)

	Sf4961 = -2*Sf4960 + Sf4959*MHiggs2(4)

	Sf6967 = 4*MW2*Sf4889 - MSq2Diff(2,2)**2

	Sf4890 = 2*MW2*Sf4889 + 3*MSq2Diff(2,2)**2

	Sf6966 = MW2*Sf4889 - MSq2Diff(2,2)**2

	Sf5047 = 4*MW2*Sf4879 + 7*MSq2Diff(2,2)**2

	Sf4880 = 2*MW2*Sf4879 + MSq2Diff(2,2)**2

	Sf6957 = 4*MW2*Sf4855 - MSq2Diff(1,2)**2

	Sf4856 = 2*MW2*Sf4855 + 3*MSq2Diff(1,2)**2

	Sf6956 = MW2*Sf4855 - MSq2Diff(1,2)**2

	Sf4846 = MW2*Sf4845 + 5*Sf91*MSq2Diff(1,2)

	Sf4921 = 4*MW2*Sf4843 + 7*MSq2Diff(1,2)**2

	Sf4844 = 2*MW2*Sf4843 + MSq2Diff(1,2)**2

	Sf489 = Sf477 + 3*C2B*Sf476*MSl2Diff(1,1)

	Sf488 = Sf474 + 9*C2B*Sf471*MSl2Diff(1,1)

	Sf8778 = 21*MZ2*Sf47 + 20*MSl2Diff(2,1)**2

	Sf49 = 7*MZ2*Sf47 + 8*MSl2Diff(2,1)**2

	Sf1374 = 27*MW2*Sf47 + 22*MSl2Diff(2,1)**2

	Sf1381 = 9*MW2*Sf47 + 10*MSl2Diff(2,1)**2

	Sf1984 = 4*MW2*Sf47 + 3*MSl2Diff(2,1)**2

	Sf1983 = MW2*Sf47 + 2*MSl2Diff(2,1)**2

	Sf1271 = 5*Sf47*MHiggs2(4) + 2*MSl2Diff(2,1)**2

	Sf463 = MW2*Sf462 + MSl2Diff(1,1)**2

	Sf456 = 2*MW2*Sf455 - MSl2Diff(1,1)**2

	Sf4237 = MW2*Sf4236 - MSq2Diff(2,1)**2

	Sf4233 = MW2*Sf4232 - MSq2Diff(1,1)**2

	Sf838 = -(C2B*Sf818) - 16*Sf34*MSl2Diff(1,1)

	Sf8771 = 21*MZ2*Sf34 + 20*MSl2Diff(1,1)**2

	Sf36 = 7*MZ2*Sf34 + 8*MSl2Diff(1,1)**2

	Sf457 = 27*MW2*Sf34 + 22*MSl2Diff(1,1)**2

	Sf464 = 9*MW2*Sf34 + 10*MSl2Diff(1,1)**2

	Sf1925 = 4*MW2*Sf34 + 3*MSl2Diff(1,1)**2

	Sf1924 = MW2*Sf34 + 2*MSl2Diff(1,1)**2

	Sf311 = 5*Sf34*MHiggs2(4) + 2*MSl2Diff(1,1)**2

	Sf4268 = MW2*Sf2973 + 3*MSq2Diff(2,1)**2

	Sf4266 = MW2*Sf2821 + 3*MSq2Diff(1,1)**2

	Sf2647 = -2*Sf2646 + Sf2645*MHiggs2(4)

	Sf2433 = -2*Sf2432 + Sf2431*MHiggs2(4)

	Sf4248 = 4*MW2*Sf2300 - MSq2Diff(2,1)**2

	Sf2301 = 2*MW2*Sf2300 + 3*MSq2Diff(2,1)**2

	Sf4247 = MW2*Sf2300 - MSq2Diff(2,1)**2

	Sf2593 = 4*MW2*Sf2282 + 7*MSq2Diff(2,1)**2

	Sf2283 = 2*MW2*Sf2282 + MSq2Diff(2,1)**2

	Sf4239 = 4*MW2*Sf2249 - MSq2Diff(1,1)**2

	Sf2250 = 2*MW2*Sf2249 + 3*MSq2Diff(1,1)**2

	Sf4238 = MW2*Sf2249 - MSq2Diff(1,1)**2

	Sf2354 = 4*MW2*Sf2219 + 7*MSq2Diff(1,1)**2

	Sf2220 = 2*MW2*Sf2219 + MSq2Diff(1,1)**2

	Sf1749 = CA*Sf1748 + SA*SB*Ss435

	Sf1735 = SB*Sf1734 + CAB*CB2*Ss443*UStau21

	Sf1722 = Sf1705 + 18*C2B*Sf1704*MSl2Diff(2,1)

	Sf1721 = 115*Sf1701 + 27*C2B*Sf1698*MSl2Diff(2,1)

	Sf1657 = 4*Sf1646 - 3*C2B*Sf1644*MSl2Diff(2,1)**2

	Sf1656 = C2B*Sf1645 + 2*Sf1641*MSl2Diff(2,1)

	Sf1655 = C2B*Sf1643 + 16*Sf47*MSl2Diff(2,1)

	Sf1660 = Sf1642 + 9*C2B*Sf1640*MSl2Diff(2,1)

	Sf1406 = Sf1394 + 3*C2B*Sf1393*MSl2Diff(2,1)

	Sf1405 = Sf1391 + 9*C2B*Sf1388*MSl2Diff(2,1)

	Sf1380 = MW2*Sf1379 + MSl2Diff(2,1)**2

	Sf1373 = 2*MW2*Sf1372 - MSl2Diff(2,1)**2

	Sf10473 = Sf10471 - MHiggs2(2)*MStop2(1)

	Sf10439 = Sf10438 - MHiggs2(2)*MStop2(2)

	Sf10372 = Sf10371 - MHiggs2(1)*MStop2(1)

	Sf10345 = Sf10344 - MHiggs2(1)*MStop2(2)

	Sf10221 = CB*MW*MZ*SAB*Sf10219 - 2*CW*ML*Sf10220

	Sf10088 = Sf10086 - MHiggs2(1)*MStau2(1)

	Sf10058 = Sf10057 - MHiggs2(1)*MStau2(2)

	Sf942 = ML2*SA*SB + CAB*CB2*MZ2*SW2

	Sf130 = 3*MBot2*SA*SB + CAB*CB2*MZ2*SW2

	Sf7190 = 3*MTop2*SA + 2*CAB*MZ2*SB*SW2

	Sf1744 = ML2*SA - CB*MZ2*SAB*SW2

	Sf4764 = 3*MBot2*SA - CB*MZ2*SAB*SW2

	Sf1153 = CB2*SA + CAB*SB

	Sf1740 = 2*SA - 3*CA*TB

	Sf1741 = SA + 3*CA*TB

	Sf2192 = MW2*S2B - MBot2*TB

	Sf127 = CA*S2B + C2B*SA

	Sf2316 = S2B + 2*CW2*S2B

	Sf6 = S2B - 2*CW2*S2B

	Sf9936 = MZ2**2*(-57*MZ2 + 357*MHiggs2(2))

	Sf10472 = MZ2**2 - MHiggs2(2)*MStop2(2)

	Sf10330 = MZ2**2 - MHiggs2(1)*MStop2(2)

	Sf10437 = MZ2**2 - MHiggs2(2)*MStop2(1)

	Sf10321 = MZ2**2 - MHiggs2(1)*MStop2(1)

	Sf9351 = MZ2**2 - MHiggs2(2)*MStau2(2)

	Sf10087 = MZ2**2 - MHiggs2(1)*MStau2(2)

	Sf9318 = MZ2**2 - MHiggs2(2)*MStau2(1)

	Sf10056 = MZ2**2 - MHiggs2(1)*MStau2(1)

	Sf8959 = MZ2**2 - MHiggs2(2)*MSbot2(2)

	Sf9149 = MZ2**2 - MHiggs2(1)*MSbot2(2)

	Sf8905 = MZ2**2 - MHiggs2(2)*MSbot2(1)

	Sf9109 = MZ2**2 - MHiggs2(1)*MSbot2(1)

	Sf8931 = MZ2*(-3*MZ2**2 + 3*MHiggs2(2)**2)

	Sf10005 = 831*MZ2*MHiggs2(2) - 192*MStau2(2)**2

	Sf9904 = 831*MZ2*MHiggs2(2) - 192*MStau2(1)**2

	Sf9945 = MZ2*(-222*MZ2 + 202*MHiggs2(2))

	Sf9978 = 5*MZ2*MHiggs2(2) - 8*MStau2(2)**2

	Sf9811 = 5*MZ2*MHiggs2(2) + 37*Ss143*MStau2(2)

	Sf9870 = 5*MZ2*MHiggs2(2) - 8*MStau2(1)**2

	Sf9678 = 5*MZ2*MHiggs2(2) + 37*Ss143*MStau2(1)

	Sf10248 = MZ2*MHiggs2(2) - 2*Sf10247*MStop2(2)

	Sf10249 = MZ2*MHiggs2(2) - 6*Sf10247*MStop2(2)

	Sf10228 = MZ2*MHiggs2(2) - 2*Sf10227*MStop2(1)

	Sf10229 = MZ2*MHiggs2(2) - 6*Sf10227*MStop2(1)

	Sf9693 = MZ2*MHiggs2(2) - Ss143*MStau2(2)

	Sf9832 = MZ2*MHiggs2(2) - 7*Ss143*MStau2(2)

	Sf9305 = MZ2*MHiggs2(2) - 2*Sf9304*MStau2(2)

	Sf9306 = MZ2*MHiggs2(2) - 6*Sf9304*MStau2(2)

	Sf9398 = MZ2*MHiggs2(2) - Ss143*MStau2(1)

	Sf9647 = MZ2*MHiggs2(2) - 7*Ss143*MStau2(1)

	Sf9293 = MZ2*MHiggs2(2) - 2*Sf9292*MStau2(1)

	Sf9294 = MZ2*MHiggs2(2) - 6*Sf9292*MStau2(1)

	Sf8888 = MZ2*MHiggs2(2) - 2*Sf8887*MSbot2(2)

	Sf8889 = MZ2*MHiggs2(2) - 6*Sf8887*MSbot2(2)

	Sf8876 = MZ2*MHiggs2(2) - 2*Sf8875*MSbot2(1)

	Sf8877 = MZ2*MHiggs2(2) - 6*Sf8875*MSbot2(1)

	Sf9133 = MZ2*(-3*MZ2**2 + 3*MHiggs2(1)**2)

	Sf10253 = MZ2*MHiggs2(1) - 2*Sf10252*MStop2(2)

	Sf10254 = MZ2*MHiggs2(1) - 6*Sf10252*MStop2(2)

	Sf10233 = MZ2*MHiggs2(1) - 2*Sf10232*MStop2(1)

	Sf10234 = MZ2*MHiggs2(1) - 6*Sf10232*MStop2(1)

	Sf9309 = MZ2*MHiggs2(1) - 2*Sf9308*MStau2(2)

	Sf9310 = MZ2*MHiggs2(1) - 6*Sf9308*MStau2(2)

	Sf9297 = MZ2*MHiggs2(1) - 2*Sf9296*MStau2(1)

	Sf9298 = MZ2*MHiggs2(1) - 6*Sf9296*MStau2(1)

	Sf8892 = MZ2*MHiggs2(1) - 2*Sf8891*MSbot2(2)

	Sf8893 = MZ2*MHiggs2(1) - 6*Sf8891*MSbot2(2)

	Sf8880 = MZ2*MHiggs2(1) - 2*Sf8879*MSbot2(1)

	Sf8881 = MZ2*MHiggs2(1) - 6*Sf8879*MSbot2(1)

	Sf10036 = 44*MZ2 - 8*MStau2(2)

	Sf9960 = 44*MZ2 - 8*MStau2(1)

	Sf7758 = 41*MZ2 + 4*MSl2Diff(2,1)

	Sf8728 = 23*MZ2 + 8*MSq2Diff(2,2)

	Sf8513 = 23*MZ2 - 8*MSq2Diff(2,2)

	Sf8510 = 23*MZ2 - 8*MSq2Diff(2,1)

	Sf8725 = 23*MZ2 + 8*MSq2Diff(1,2)

	Sf8497 = 23*MZ2 - 8*MSq2Diff(1,2)

	Sf8494 = 23*MZ2 - 8*MSq2Diff(1,1)

	Sf7640 = 23*MZ2 + 8*MSl2Diff(2,1)

	Sf7637 = 23*MZ2 + 8*MSl2Diff(1,1)

	Sf8807 = 23*MZ2 + 16*Sf99

	Sf8801 = 23*MZ2 + 16*Sf91

	Sf8794 = 23*MZ2 + 16*Sf77

	Sf8788 = 23*MZ2 + 16*Sf69

	Sf8777 = 23*MZ2 + 16*Sf47

	Sf8770 = 23*MZ2 + 16*Sf34

	Sf7759 = 21*MZ2 + 20*MSl2Diff(2,1)

	Sf7770 = 21*MZ2 + 20*MSl2Diff(1,1)

	Sf9825 = 15*MZ2 + 28*MStau2(2)

	Sf9644 = 15*MZ2 + 28*MStau2(1)

	Sf9980 = 7*MZ2 + 8*MStau2(2)

	Sf9887 = 7*MZ2 + 8*MStau2(1)

	Sf9833 = 6*MZ2 + 20*MStau2(2)

	Sf9804 = 6*MZ2 + 19*MStau2(2)

	Sf9649 = 6*MZ2 + 20*MStau2(1)

	Sf9674 = 6*MZ2 + 19*MStau2(1)

	Sf10464 = 4*MZ2 - MStop2(2)

	Sf8894 = 4*MZ2 + MSbot2(2)

	Sf8951 = 4*MZ2 - MSbot2(2)

	Sf10463 = 3*MZ2 + MStop2(2)

	Sf10460 = 3*MZ2 + MStop2(1)

	Sf9698 = 3*MZ2 - 4*MStau2(2)

	Sf9694 = 3*MZ2 - 8*MStau2(2)

	Sf9812 = 3*MZ2 - 11*MStau2(2)

	Sf9410 = 3*MZ2 - 4*MStau2(1)

	Sf9401 = 3*MZ2 - 8*MStau2(1)

	Sf9680 = 3*MZ2 - 11*MStau2(1)

	Sf8727 = 3*MZ2 + 5*MSq2Diff(2,2)

	Sf8514 = 3*MZ2 - 5*MSq2Diff(2,2)

	Sf8511 = 3*MZ2 - 5*MSq2Diff(2,1)

	Sf8724 = 3*MZ2 + 5*MSq2Diff(1,2)

	Sf8498 = 3*MZ2 - 5*MSq2Diff(1,2)

	Sf8495 = 3*MZ2 - 5*MSq2Diff(1,1)

	Sf7641 = 3*MZ2 + 5*MSl2Diff(2,1)

	Sf7638 = 3*MZ2 + 5*MSl2Diff(1,1)

	Sf8950 = 3*MZ2 + MSbot2(2)

	Sf8947 = 3*MZ2 + MSbot2(1)

	Sf8679 = 2*MZ2 - 13*MStop2(2)

	Sf7993 = 2*MZ2 - 13*MStop2(1)

	Sf7735 = 2*MZ2 - 13*MStau2(2)

	Sf7582 = 2*MZ2 - 13*MStau2(1)

	Sf7498 = 2*MZ2 - 13*MSneu2(1)

	Sf7524 = 2*MZ2 - 13*MSbot2(2)

	Sf7505 = 2*MZ2 - 13*MSbot2(1)

	Sf100 = 2*MZ2 + 13*Sf99

	Sf92 = 2*MZ2 + 13*Sf91

	Sf78 = 2*MZ2 + 13*Sf77

	Sf70 = 2*MZ2 + 13*Sf69

	Sf48 = 2*MZ2 + 13*Sf47

	Sf35 = 2*MZ2 + 13*Sf34

	Sf8521 = MZ2 + MStop2(4)

	Sf104 = MZ2 - MStop2(4)

	Sf7219 = MZ2 + MStop2(3)

	Sf7199 = MZ2 - MStop2(3)

	Sf8678 = MZ2 + 2*MStop2(2)

	Sf10274 = MZ2 + MStop2(2)

	Sf10273 = MZ2 - MStop2(2)

	Sf8699 = MZ2 - 2*MStop2(2)

	Sf113 = MZ2 - 4*MStop2(2)

	Sf7992 = MZ2 + 2*MStop2(1)

	Sf10265 = MZ2 + MStop2(1)

	Sf10264 = MZ2 - MStop2(1)

	Sf8660 = MZ2 - 2*MStop2(1)

	Sf82 = MZ2 - 4*MStop2(1)

	Sf7702 = MZ2 + MStau2(4)

	Sf52 = MZ2 - MStau2(4)

	Sf8536 = MZ2 - 3*MStau2(4)

	Sf1775 = MZ2 + MStau2(3)

	Sf1755 = MZ2 - MStau2(3)

	Sf9987 = MZ2 + 8*MStau2(2)

	Sf7734 = MZ2 + 2*MStau2(2)

	Sf9342 = MZ2 + MStau2(2)

	Sf9341 = MZ2 - MStau2(2)

	Sf7750 = MZ2 - 2*MStau2(2)

	Sf8538 = MZ2 - 3*MStau2(2)

	Sf61 = MZ2 - 4*MStau2(2)

	Sf10016 = MZ2 - 16*MStau2(2)

	Sf9877 = MZ2 + 8*MStau2(1)

	Sf7581 = MZ2 + 2*MStau2(1)

	Sf9315 = MZ2 + MStau2(1)

	Sf9314 = MZ2 - MStau2(1)

	Sf7728 = MZ2 - 2*MStau2(1)

	Sf39 = MZ2 - 4*MStau2(1)

	Sf8534 = MZ2 - 6*MStau2(1)

	Sf9923 = MZ2 - 16*MStau2(1)

	Sf7499 = MZ2 + 2*MSneu2(1)

	Sf1 = MZ2 - 4*MSneu2(1)

	Sf1173 = MZ2 + 4*MSl2Diff(2,1)

	Sf1171 = MZ2 + 4*MSl2Diff(1,1)

	Sf7510 = MZ2 + MSbot2(4)

	Sf17 = MZ2 - MSbot2(4)

	Sf4795 = MZ2 + MSbot2(3)

	Sf4773 = MZ2 - MSbot2(3)

	Sf7523 = MZ2 + 2*MSbot2(2)

	Sf8944 = MZ2 + MSbot2(2)

	Sf8943 = MZ2 - MSbot2(2)

	Sf8354 = MZ2 - 2*MSbot2(2)

	Sf26 = MZ2 - 4*MSbot2(2)

	Sf7504 = MZ2 + 2*MSbot2(1)

	Sf8902 = MZ2 + MSbot2(1)

	Sf8901 = MZ2 - MSbot2(1)

	Sf8257 = MZ2 - 2*MSbot2(1)

	Sf11 = MZ2 - 4*MSbot2(1)

	Sf58 = MZ2 - Sf57**2

	Sf9709 = MZ2 + Sf55*Sf57

	Sf9418 = MZ2 - Sf55*Sf57

	Sf56 = MZ2 - Sf55**2

	Sf23 = MZ2 - Sf22**2

	Sf21 = MZ2 - Sf20**2

	Sf110 = MZ2 - Sf109**2

	Sf108 = MZ2 - Sf107**2

	Sf10246 = MZ + MStop(2)

	Sf10245 = MZ - MStop(2)

	Sf10226 = MZ + MStop(1)

	Sf10225 = MZ - MStop(1)

	Sf9303 = MZ + MStau(2)

	Sf9302 = MZ - MStau(2)

	Sf9291 = MZ + MStau(1)

	Sf9290 = MZ - MStau(1)

	Sf8886 = MZ + MSbot(2)

	Sf8885 = MZ - MSbot(2)

	Sf8874 = MZ + MSbot(1)

	Sf8873 = MZ - MSbot(1)

	Sf1639 = 324*C2B*MW2**6 - 132*MW2*MSl2Diff(2,1)**5

	Sf814 = 324*C2B*MW2**6 - 132*MW2*MSl2Diff(1,1)**5

	Sf2132 = 3*MW2**2 - MStau2(2)**2

	Sf2052 = 3*MW2**2 - MStau2(1)**2

	Sf4518 = 3*MW2**2 - MSbot2(2)**2

	Sf4202 = 3*MW2**2 + Sf4201*MSbot2(2)

	Sf4414 = 3*MW2**2 - MSbot2(1)**2

	Sf4122 = 3*MW2**2 + Sf4121*MSbot2(1)

	Sf6138 = MW2*MHiggs2(4) + Sf6137*MSbot2(2)

	Sf5755 = MW2*MHiggs2(4) + Sf5754*MSbot2(1)

	Sf6089 = 26*MW2 + 27*MSbot2(2)

	Sf5699 = 26*MW2 + 27*MSbot2(1)

	Sf5847 = 19*MW2 + 30*MStop2(2)

	Sf4063 = 19*MW2 + 30*MStop2(1)

	Sf3665 = 16*MW2 + 21*MSbot2(2)

	Sf3385 = 16*MW2 + 21*MSbot2(1)

	Sf6096 = 13*MW2 + 18*MSbot2(2)

	Sf5707 = 13*MW2 + 18*MSbot2(1)

	Sf3570 = 10*MW2 + 7*MSbot2(2)

	Sf3218 = 10*MW2 + 7*MSbot2(1)

	Sf6310 = 9*MW2 + 10*MSbot2(1)

	Sf3978 = 8*MW2 + 9*MSbot2(2)

	Sf3881 = 8*MW2 + 9*MSbot2(1)

	Sf1669 = 7*MW2 + MStau2(2)

	Sf853 = 7*MW2 + MStau2(1)

	Sf3607 = 7*MW2 + MSbot2(2)

	Sf3280 = 7*MW2 + MSbot2(1)

	Sf1294 = 5*MW2 + 11*MStau2(2)

	Sf1207 = 5*MW2 + 4*MStau2(2)

	Sf2134 = 5*MW2 + 3*MStau2(2)

	Sf1652 = 5*MW2 - 21*MStau2(2)

	Sf342 = 5*MW2 + 11*MStau2(1)

	Sf211 = 5*MW2 + 4*MStau2(1)

	Sf2054 = 5*MW2 + 3*MStau2(1)

	Sf835 = 5*MW2 - 21*MStau2(1)

	Sf3048 = 5*MW2 + 11*MSbot2(2)

	Sf4522 = 5*MW2 + 3*MSbot2(2)

	Sf4734 = 5*MW2 - 21*MSbot2(2)

	Sf2921 = 5*MW2 + 11*MSbot2(1)

	Sf4418 = 5*MW2 + 3*MSbot2(1)

	Sf4639 = 5*MW2 - 21*MSbot2(1)

	Sf1805 = 4*MW2 + 11*MStau2(2)

	Sf1424 = 4*MW2 + 5*MStau2(2)

	Sf1186 = 4*MW2 - MStau2(2)

	Sf993 = 4*MW2 + 11*MStau2(1)

	Sf513 = 4*MW2 + 5*MStau2(1)

	Sf156 = 4*MW2 - MStau2(1)

	Sf4520 = 4*MW2 + MSbot2(2)

	Sf4416 = 4*MW2 + MSbot2(1)

	Sf86 = 4*MW2 - MZ2

	Sf5507 = 3*MW2 + MStop2(2)

	Sf3128 = 3*MW2 + MStop2(1)

	Sf1592 = 3*MW2 - MStau2(2)

	Sf731 = 3*MW2 - MStau2(1)

	Sf6120 = 3*MW2 + 4*MSbot2(2)

	Sf6726 = 3*MW2 - MSbot2(2)

	Sf5733 = 3*MW2 + 4*MSbot2(1)

	Sf6459 = 3*MW2 - MSbot2(1)

	Sf1718 = 2*MW2 + 43*MStau2(2)

	Sf1594 = 2*MW2 + 3*MStau2(2)

	Sf1428 = 2*MW2 + MStau2(2)

	Sf1413 = 2*MW2 - 5*MStau2(2)

	Sf928 = 2*MW2 + 43*MStau2(1)

	Sf734 = 2*MW2 + 3*MStau2(1)

	Sf517 = 2*MW2 + MStau2(1)

	Sf496 = 2*MW2 - 5*MStau2(1)

	Sf1039 = 2*MW2 + MSneu2(1)

	Sf1768 = 2*MW2 - MSl2Diff(2,1)

	Sf1738 = 2*MW2 - MSl2Diff(1,1)

	Sf3575 = 2*MW2 + 43*MSbot2(2)

	Sf6728 = 2*MW2 + 3*MSbot2(2)

	Sf3487 = 2*MW2 - 5*MSbot2(2)

	Sf3223 = 2*MW2 + 43*MSbot2(1)

	Sf6461 = 2*MW2 + 3*MSbot2(1)

	Sf3104 = 2*MW2 - 5*MSbot2(1)

	Sf15 = 2*MW2 + MZ2

	Sf5167 = MW2 + 3*MStop2(2)

	Sf5038 = MW2 + MStop2(2)

	Sf5029 = MW2 - MStop2(2)

	Sf2584 = MW2 + 3*MStop2(1)

	Sf2573 = MW2 + MStop2(1)

	Sf2563 = MW2 - MStop2(1)

	Sf1196 = MW2 - MStau2(2)

	Sf167 = MW2 - MStau2(1)

	Sf7257 = MW2 - (MSbot(2) + MStop(2))**2

	Sf7252 = MW2 - (MSbot(2) + MStop(1))**2

	Sf7255 = MW2 - (MSbot(1) + MStop(2))**2

	Sf7250 = MW2 - (MSbot(1) + MStop(1))**2

	Sf7256 = MW2 - (MSbot(2) - MStop(2))**2

	Sf7251 = MW2 - (MSbot(2) - MStop(1))**2

	Sf7254 = MW2 - (MSbot(1) - MStop(2))**2

	Sf7249 = MW2 - (MSbot(1) - MStop(1))**2

	Sf4803 = MW2 + MSq2Diff(2,2)

	Sf4802 = MW2 - MSq2Diff(2,2)

	Sf3064 = MW2 + MSq2Diff(2,1)

	Sf2741 = MW2 - MSq2Diff(2,1)

	Sf4782 = MW2 + MSq2Diff(1,2)

	Sf4781 = MW2 - MSq2Diff(1,2)

	Sf3056 = MW2 + MSq2Diff(1,1)

	Sf2564 = MW2 - MSq2Diff(1,1)

	Sf7261 = MW2 - 3*MSneu2(1)

	Sf7246 = MW2 - (MSneu(1) + MStau(2))**2

	Sf7243 = MW2 - (MSneu(1) + MStau(1))**2

	Sf7245 = MW2 - (-MSneu(1) + MStau(2))**2

	Sf7242 = MW2 - (-MSneu(1) + MStau(1))**2

	Sf1458 = MW2 + MSl2Diff(2,1)

	Sf1459 = MW2 - MSl2Diff(2,1)

	Sf555 = MW2 + MSl2Diff(1,1)

	Sf556 = MW2 - MSl2Diff(1,1)

	Sf4200 = MW2 + MSbot2(2)

	Sf2722 = MW2 - MSbot2(2)

	Sf6115 = MW2 - 5*MSbot2(2)

	Sf5464 = MW2 + 6*MSbot2(1)

	Sf4120 = MW2 + MSbot2(1)

	Sf2543 = MW2 - MSbot2(1)

	Sf5727 = MW2 - 5*MSbot2(1)

	Sf4872 = MW2 + 2*Sf99

	Sf7266 = MW2 - 3*Sf99

	Sf5201 = MW2 - 8*Sf99

	Sf4837 = MW2 + 2*Sf91

	Sf7265 = MW2 - 3*Sf91

	Sf5189 = MW2 - 8*Sf91

	Sf2275 = MW2 + 2*Sf77

	Sf7263 = MW2 - 3*Sf77

	Sf5181 = MW2 - 8*Sf77

	Sf2212 = MW2 + 2*Sf69

	Sf7262 = MW2 - 3*Sf69

	Sf5175 = MW2 - 8*Sf69

	Sf5850 = MW + MStop(2)

	Sf5849 = MW - MStop(2)

	Sf4066 = MW + MStop(1)

	Sf4065 = MW - MStop(1)

	Sf1801 = MW + MStau(2)

	Sf1800 = MW - MStau(2)

	Sf977 = MW + MStau(1)

	Sf976 = MW - MStau(1)

	Sf3561 = MW + MSbot(2)

	Sf3560 = MW - MSbot(2)

	Sf3205 = MW + MSbot(1)

	Sf3204 = MW - MSbot(1)

	Sf2194 = reMUE + TB*Re(Abot)

	Sf2245 = reMUE - TB*Re(Atop)

	Sf160 = reMUE + TB*Re(Atau)

	Sf3009 = 10*MBot2*MTop2 - 3*Sf3008*MSbot2(2)

	Sf2877 = 10*MBot2*MTop2 - 3*Sf2876*MSbot2(1)

	Sf2318 = 6*CW2*MTop2 - 4*MW2*S2B*SW2*TB

	Sf2314 = 3*CW2*MTop2 + 2*C2B*MW2*SW2

	Sf7194 = 3*CA*MTop2 - 2*MZ2*SAB*SB*SW2

	Sf4714 = 84*MTop2 + 31*MSbot2(2)

	Sf4619 = 84*MTop2 + 31*MSbot2(1)

	Sf3028 = 44*MTop2 - 23*MSbot2(2)

	Sf2901 = 44*MTop2 - 23*MSbot2(1)

	Sf3029 = 10*MTop2 - 21*MSbot2(2)

	Sf4715 = 10*MTop2 - 109*MSbot2(2)

	Sf2902 = 10*MTop2 - 21*MSbot2(1)

	Sf4620 = 10*MTop2 - 109*MSbot2(1)

	Sf3920 = 6*MTop2 - 4*MSbot2(2)

	Sf3796 = 6*MTop2 - 4*MSbot2(1)

	Sf123 = 6*MTop2 - 4*MW2*S2B*TB

	Sf4711 = 5*MTop2 + 114*MSbot2(2)

	Sf3025 = 5*MTop2 - 78*MSbot2(2)

	Sf4615 = 5*MTop2 + 114*MSbot2(1)

	Sf2898 = 5*MTop2 - 78*MSbot2(1)

	Sf2193 = MTop2 - Sf2192*TB

	Sf345 = ML2**2 + C2B**2*MW2**2

	Sf150 = CW2*ML2 - C2B*MW2*SW2

	Sf1750 = CA*ML2 - CAB*CB*MZ2*SW2

	Sf1622 = 9*ML2 + 10*MSl2Diff(2,1)

	Sf796 = 9*ML2 + 10*MSl2Diff(1,1)

	Sf1247 = 7*ML2 + 6*MStau2(2)

	Sf273 = 7*ML2 + 6*MStau2(1)

	Sf1256 = 5*ML2 - 22*MStau2(2)

	Sf289 = 5*ML2 - 22*MStau2(1)

	Sf1555 = 3*ML2 - 26*MStau2(2)

	Sf685 = 3*ML2 - 26*MStau2(1)

	Sf2109 = 2*ML2 - MStau2(2)

	Sf1974 = 2*ML2 - 39*MStau2(2)

	Sf2025 = 2*ML2 - MStau2(1)

	Sf1914 = 2*ML2 - 39*MStau2(1)

	Sf1560 = ML2 - 4*C2B*MStau2(2)

	Sf1549 = ML2 + 42*MStau2(2)

	Sf1426 = ML2 + MStau2(2)

	Sf1572 = ML2 - 3*MStau2(2)

	Sf1497 = ML2 - 4*MStau2(2)

	Sf2110 = ML2 - 7*MStau2(2)

	Sf1480 = ML2 - 13*MStau2(2)

	Sf1553 = ML2 - 23*MStau2(2)

	Sf691 = ML2 - 4*C2B*MStau2(1)

	Sf676 = ML2 + 42*MStau2(1)

	Sf515 = ML2 + MStau2(1)

	Sf706 = ML2 - 3*MStau2(1)

	Sf597 = ML2 - 4*MStau2(1)

	Sf2026 = ML2 - 7*MStau2(1)

	Sf579 = ML2 - 13*MStau2(1)

	Sf681 = ML2 - 23*MStau2(1)

	Sf1559 = ML2 - MSl2Diff(2,1)

	Sf1485 = ML2 - 3*MSl2Diff(2,1)

	Sf1530 = ML2 - 9*MSl2Diff(2,1)

	Sf689 = ML2 - MSl2Diff(1,1)

	Sf584 = ML2 - 3*MSl2Diff(1,1)

	Sf636 = ML2 - 9*MSl2Diff(1,1)

	Sf241 = ML2 - C2B*MW2

	Sf6197 = -64*MW2**3*Sf99 + 20*MHiggs2(4)**4

	Sf5819 = -64*MW2**3*Sf91 + 20*MHiggs2(4)**4

	Sf4128 = -64*MW2**3*Sf77 + 20*MHiggs2(4)**4

	Sf4008 = -64*MW2**3*Sf69 + 20*MHiggs2(4)**4

	Sf4302 = 16*MW2**3 + 9*MHiggs2(4)**3

	Sf638 = 2*MW2**3 + MHiggs2(4)**3

	Sf4286 = 8*MW2**2 + 9*MHiggs2(4)**2

	Sf4007 = 4*MW2**2 + 5*MHiggs2(4)**2

	Sf1137 = -160*MW2 + 137*MHiggs2(4)

	Sf1110 = -151*MW2 + 96*MHiggs2(4)

	Sf1133 = -92*MW2 + 73*MHiggs2(4)

	Sf1125 = -68*MW2 + 59*MHiggs2(4)

	Sf1090 = -77*MW2 + 52*MHiggs2(4)

	Sf1105 = -79*MW2 + 44*MHiggs2(4)

	Sf1103 = -29*MW2 + 8*MHiggs2(4)

	Sf203 = -17*MW2 + 7*MHiggs2(4)

	Sf4832 = -4*MW2 + 6*MHiggs2(4)

	Sf3450 = 4*MW2 + 5*MHiggs2(4)

	Sf221 = -8*MW2 + 5*MHiggs2(4)

	Sf2491 = -5*MW2 + 4*MHiggs2(4)

	Sf1114 = -13*MW2 + 4*MHiggs2(4)

	Sf4044 = 5*MW2 + 3*MHiggs2(4)

	Sf733 = -MW2 + 3*MHiggs2(4)

	Sf5171 = -2*MW2 + 3*MHiggs2(4)

	Sf540 = -5*MW2 + 3*MHiggs2(4)

	Sf191 = -7*MW2 + 3*MHiggs2(4)

	Sf2201 = MW2 + 2*MHiggs2(4)

	Sf3731 = -MW2 + 2*MHiggs2(4)

	Sf4833 = MHiggs2(4) + MStop2(2)

	Sf5172 = MHiggs2(4) + MStop2(1)

	Sf1235 = MHiggs2(4) + 18*MStau2(2)

	Sf1873 = MHiggs2(4) + 5*MStau2(2)

	Sf1185 = MHiggs2(4) - MStau2(2)

	Sf253 = MHiggs2(4) + 18*MStau2(1)

	Sf1122 = MHiggs2(4) + 5*MStau2(1)

	Sf155 = MHiggs2(4) - MStau2(1)

	Sf4894 = MHiggs2(4) + MSq2Diff(2,2)

	Sf4893 = MHiggs2(4) - MSq2Diff(2,2)

	Sf2305 = MHiggs2(4) + MSq2Diff(2,1)

	Sf2304 = MHiggs2(4) - MSq2Diff(2,1)

	Sf4860 = MHiggs2(4) + MSq2Diff(1,2)

	Sf4859 = MHiggs2(4) - MSq2Diff(1,2)

	Sf2254 = MHiggs2(4) + MSq2Diff(1,1)

	Sf2253 = MHiggs2(4) - MSq2Diff(1,1)

	Sf1040 = MHiggs2(4) + MSneu2(1)

	Sf1220 = MHiggs2(4) + MSl2Diff(2,1)

	Sf1221 = MHiggs2(4) - MSl2Diff(2,1)

	Sf229 = MHiggs2(4) + MSl2Diff(1,1)

	Sf230 = MHiggs2(4) - MSl2Diff(1,1)

	Sf3740 = MHiggs2(4) - 3*MSbot2(2)

	Sf2605 = MHiggs2(4) - 18*MSbot2(2)

	Sf3730 = MHiggs2(4) - 3*MSbot2(1)

	Sf2367 = MHiggs2(4) - 18*MSbot2(1)

	Sf4224 = MW2*TB2 + MHiggs2(4)

	Sf1234 = -2*Sf47 + MHiggs2(4)

	Sf252 = -2*Sf34 + MHiggs2(4)

	Sf1953 = 5*MW2 + MHiggs2(4)

	Sf613 = 3*MW2 + MHiggs2(4)

	Sf2586 = 2*MW2 + MHiggs2(4)

	Sf173 = -2*MW2 + MHiggs2(4)

	Sf180 = -3*MW2 + MHiggs2(4)

	Sf1948 = -4*MW2 + MHiggs2(4)

	Sf1239 = MHiggs(4) + MStau(2)

	Sf1238 = MHiggs(4) - MStau(2)

	Sf260 = MHiggs(4) + MStau(1)

	Sf259 = MHiggs(4) - MStau(1)

	Sf5055 = MHiggs(4) + MSbot(2) + MStop(2)

	Sf5054 = MHiggs(4) - MSbot(2) - MStop(2)

	Sf2601 = MHiggs(4) + MSbot(2) + MStop(1)

	Sf2600 = MHiggs(4) - MSbot(2) - MStop(1)

	Sf4929 = MHiggs(4) + MSbot(1) + MStop(2)

	Sf4928 = MHiggs(4) - MSbot(1) - MStop(2)

	Sf2363 = MHiggs(4) + MSbot(1) + MStop(1)

	Sf2362 = MHiggs(4) - MSbot(1) - MStop(1)

	Sf5052 = MHiggs(4) + MSbot(2) - MStop(2)

	Sf5053 = MHiggs(4) - MSbot(2) + MStop(2)

	Sf2598 = MHiggs(4) + MSbot(2) - MStop(1)

	Sf2599 = MHiggs(4) - MSbot(2) + MStop(1)

	Sf4926 = MHiggs(4) + MSbot(1) - MStop(2)

	Sf4927 = MHiggs(4) - MSbot(1) + MStop(2)

	Sf2360 = MHiggs(4) + MSbot(1) - MStop(1)

	Sf2361 = MHiggs(4) - MSbot(1) + MStop(1)

	Sf1270 = MHiggs(4) + MSneu(1) + MStau(2)

	Sf1269 = MHiggs(4) - MSneu(1) - MStau(2)

	Sf310 = MHiggs(4) + MSneu(1) + MStau(1)

	Sf309 = MHiggs(4) - MSneu(1) - MStau(1)

	Sf1267 = MHiggs(4) - MSneu(1) + MStau(2)

	Sf1268 = MHiggs(4) + MSneu(1) - MStau(2)

	Sf307 = MHiggs(4) - MSneu(1) + MStau(1)

	Sf308 = MHiggs(4) + MSneu(1) - MStau(1)

	Sf2265 = MHiggs(4) + MSbot(2)

	Sf2264 = MHiggs(4) - MSbot(2)

	Sf2200 = MHiggs(4) + MSbot(1)

	Sf2199 = MHiggs(4) - MSbot(1)

	Sf9054 = -8*MZ2**2 + 6*MHiggs2(2)**2

	Sf8979 = -2*MZ2**2 + MHiggs2(2)**2

	Sf9441 = 633*MZ2 + 1175*MHiggs2(2)

	Sf9943 = -936*MZ2 + 904*MHiggs2(2)

	Sf9433 = 695*MZ2 + 777*MHiggs2(2)

	Sf9946 = -604*MZ2 + 596*MHiggs2(2)

	Sf9938 = -684*MZ2 + 596*MHiggs2(2)

	Sf9467 = -578*MZ2 + 254*MHiggs2(2)

	Sf9449 = 5*MZ2 + 227*MHiggs2(2)

	Sf9518 = -276*MZ2 + 156*MHiggs2(2)

	Sf9932 = 193*MZ2 + 139*MHiggs2(2)

	Sf9942 = -MZ2 + 107*MHiggs2(2)

	Sf9453 = 5*MZ2 + 76*MHiggs2(2)

	Sf9526 = -110*MZ2 + 74*MHiggs2(2)

	Sf9937 = -129*MZ2 + 67*MHiggs2(2)

	Sf9469 = -83*MZ2 + 56*MHiggs2(2)

	Sf9511 = -97*MZ2 + 52*MHiggs2(2)

	Sf9797 = 44*MHiggs2(2) + 47*MStau2(2)

	Sf9661 = 44*MHiggs2(2) + 47*MStau2(1)

	Sf9445 = 11*MZ2 + 43*MHiggs2(2)

	Sf9503 = -78*MZ2 + 42*MHiggs2(2)

	Sf9949 = -38*MZ2 + 36*MHiggs2(2)

	Sf9941 = -13*MZ2 + 29*MHiggs2(2)

	Sf9931 = 23*MZ2 + 27*MHiggs2(2)

	Sf9459 = 5*MZ2 + 13*MHiggs2(2)

	Sf9712 = 12*MHiggs2(2) + 49*MStau2(2)

	Sf9426 = 12*MHiggs2(2) + 49*MStau2(1)

	Sf10562 = 11*MHiggs2(2) + 8*MStop2(2)

	Sf10554 = 11*MHiggs2(2) + 8*MStop2(1)

	Sf9069 = 11*MHiggs2(2) + 8*MSbot2(2)

	Sf9060 = 11*MHiggs2(2) + 8*MSbot2(1)

	Sf9741 = 9*MHiggs2(2) + 8*MStau2(2)

	Sf9486 = 9*MHiggs2(2) + 8*MStau2(1)

	Sf9544 = -5*MZ2 + 9*MHiggs2(2)

	Sf9081 = -6*MZ2 + 9*MHiggs2(2)

	Sf10043 = 8*MHiggs2(2) - 7*MStau2(2)

	Sf9964 = 8*MHiggs2(2) - 7*MStau2(1)

	Sf9437 = 19*MZ2 + 8*MHiggs2(2)

	Sf9528 = -11*MZ2 + 8*MHiggs2(2)

	Sf9523 = -13*MZ2 + 8*MHiggs2(2)

	Sf8925 = -4*MZ2 + 6*MHiggs2(2)

	Sf9055 = -8*MZ2 + 6*MHiggs2(2)

	Sf9024 = -4*MZ2 + 5*MHiggs2(2)

	Sf9933 = -81*MZ2 + 5*MHiggs2(2)

	Sf9818 = 4*MHiggs2(2) + MStau2(2)

	Sf9631 = 4*MHiggs2(2) + MStau2(1)

	Sf9464 = 23*MZ2 + 4*MHiggs2(2)

	Sf9496 = -7*MZ2 + 4*MHiggs2(2)

	Sf9768 = 3*MHiggs2(2) + 19*MStau2(2)

	Sf9539 = 3*MHiggs2(2) + 19*MStau2(1)

	Sf8991 = -5*Sf8990 + 3*MHiggs2(2)

	Sf8982 = -5*Sf8981 + 3*MHiggs2(2)

	Sf10492 = -5*Sf10298 + 3*MHiggs2(2)

	Sf10486 = -5*Sf10288 + 3*MHiggs2(2)

	Sf9430 = 5*MZ2 + 3*MHiggs2(2)

	Sf9490 = MZ2 + 3*MHiggs2(2)

	Sf8923 = -2*MZ2 + 3*MHiggs2(2)

	Sf9057 = -4*MZ2 + 3*MHiggs2(2)

	Sf10467 = 2*MHiggs2(2) + 5*MStop2(2)

	Sf10582 = 2*MHiggs2(2) + 3*MStop2(2)

	Sf10457 = 2*MHiggs2(2) + 5*MStop2(1)

	Sf10570 = 2*MHiggs2(2) + 3*MStop2(1)

	Sf9346 = 2*MHiggs2(2) + 5*MStau2(2)

	Sf9856 = 2*MHiggs2(2) + 3*MStau2(2)

	Sf9686 = 2*MHiggs2(2) - MStau2(2)

	Sf9338 = 2*MHiggs2(2) + 5*MStau2(1)

	Sf9845 = 2*MHiggs2(2) + 3*MStau2(1)

	Sf9394 = 2*MHiggs2(2) - MStau2(1)

	Sf8954 = 2*MHiggs2(2) + 5*MSbot2(2)

	Sf9094 = 2*MHiggs2(2) + 3*MSbot2(2)

	Sf8939 = 2*MHiggs2(2) + 5*MSbot2(1)

	Sf9083 = 2*MHiggs2(2) + 3*MSbot2(1)

	Sf9613 = -(Sf9302*Sf9303) + 2*MHiggs2(2)

	Sf9595 = -(Sf9290*Sf9291) + 2*MHiggs2(2)

	Sf9039 = -(Sf8885*Sf8886) + 2*MHiggs2(2)

	Sf9020 = -(Sf8873*Sf8874) + 2*MHiggs2(2)

	Sf10539 = -(Sf10245*Sf10246) + 2*MHiggs2(2)

	Sf10520 = -(Sf10225*Sf10226) + 2*MHiggs2(2)

	Sf9934 = -3*MZ2 + 2*MHiggs2(2)

	Sf10533 = MHiggs2(2) + MStop2(4)

	Sf10574 = MHiggs2(2) - MStop2(4)

	Sf10575 = MHiggs2(2) + MStop2(3)

	Sf10576 = MHiggs2(2) - MStop2(3)

	Sf10476 = MHiggs2(2) + MStop2(2)

	Sf10244 = MHiggs2(2) - MStop2(2)

	Sf10442 = MHiggs2(2) + MStop2(1)

	Sf10224 = MHiggs2(2) - MStop2(1)

	Sf9607 = MHiggs2(2) + MStau2(4)

	Sf9849 = MHiggs2(2) - MStau2(4)

	Sf9850 = MHiggs2(2) + MStau2(3)

	Sf9851 = MHiggs2(2) - MStau2(3)

	Sf10044 = MHiggs2(2) + 56*MStau2(2)

	Sf9994 = MHiggs2(2) + 8*MStau2(2)

	Sf9360 = MHiggs2(2) + MStau2(2)

	Sf9359 = MHiggs2(2) - MStau2(2)

	Sf10008 = MHiggs2(2) - 2*MStau2(2)

	Sf10007 = MHiggs2(2) - 16*MStau2(2)

	Sf9971 = MHiggs2(2) + 56*MStau2(1)

	Sf9897 = MHiggs2(2) + 8*MStau2(1)

	Sf9329 = MHiggs2(2) + MStau2(1)

	Sf9328 = MHiggs2(2) - MStau2(1)

	Sf9906 = MHiggs2(2) - 2*MStau2(1)

	Sf9913 = MHiggs2(2) - 16*MStau2(1)

	Sf9033 = MHiggs2(2) + MSbot2(4)

	Sf9087 = MHiggs2(2) - MSbot2(4)

	Sf9088 = MHiggs2(2) + MSbot2(3)

	Sf9089 = MHiggs2(2) - MSbot2(3)

	Sf8968 = MHiggs2(2) + MSbot2(2)

	Sf8967 = MHiggs2(2) - MSbot2(2)

	Sf8918 = MHiggs2(2) + MSbot2(1)

	Sf8917 = MHiggs2(2) - MSbot2(1)

	Sf9507 = 7*MZ2 + MHiggs2(2)

	Sf9075 = 2*MZ2 + MHiggs2(2)

	Sf8911 = -2*MZ2 + MHiggs2(2)

	Sf8940 = -3*MZ2 + MHiggs2(2)

	Sf10480 = MHiggs(2) + MStop(2)

	Sf10479 = MHiggs(2) - MStop(2)

	Sf10446 = MHiggs(2) + MStop(1)

	Sf10445 = MHiggs(2) - MStop(1)

	Sf9364 = MHiggs(2) + MStau(2)

	Sf9363 = MHiggs(2) - MStau(2)

	Sf9333 = MHiggs(2) + MStau(1)

	Sf9332 = MHiggs(2) - MStau(1)

	Sf8972 = MHiggs(2) + MSbot(2)

	Sf8971 = MHiggs(2) - MSbot(2)

	Sf8922 = MHiggs(2) + MSbot(1)

	Sf8921 = MHiggs(2) - MSbot(1)

	Sf8930 = MZ + MHiggs(2)

	Sf8929 = -MZ + MHiggs(2)

	Sf9184 = -8*MZ2**2 + 6*MHiggs2(1)**2

	Sf9169 = -2*MZ2**2 + MHiggs2(1)**2

	Sf10315 = 11*MHiggs2(1) + 8*MStop2(2)

	Sf10306 = 11*MHiggs2(1) + 8*MStop2(1)

	Sf10154 = 11*MHiggs2(1) + 8*MStau2(2)

	Sf10124 = 11*MHiggs2(1) + 8*MStau2(1)

	Sf9218 = 11*MHiggs2(1) + 8*MSbot2(2)

	Sf9189 = 11*MHiggs2(1) + 8*MSbot2(1)

	Sf9261 = -6*MZ2 + 9*MHiggs2(1)

	Sf9127 = -4*MZ2 + 6*MHiggs2(1)

	Sf9186 = -8*MZ2 + 6*MHiggs2(1)

	Sf9229 = -4*MZ2 + 5*MHiggs2(1)

	Sf9177 = -5*Sf8990 + 3*MHiggs2(1)

	Sf9171 = -5*Sf8981 + 3*MHiggs2(1)

	Sf10299 = -5*Sf10298 + 3*MHiggs2(1)

	Sf10289 = -5*Sf10288 + 3*MHiggs2(1)

	Sf10115 = -5*Sf10114 + 3*MHiggs2(1)

	Sf10107 = -5*Sf10106 + 3*MHiggs2(1)

	Sf9125 = -2*MZ2 + 3*MHiggs2(1)

	Sf9191 = -4*MZ2 + 3*MHiggs2(1)

	Sf10431 = 2*MHiggs2(1) + 3*MStop2(2)

	Sf10420 = 2*MHiggs2(1) + 3*MStop2(1)

	Sf10083 = 2*MHiggs2(1) + 5*MStau2(2)

	Sf10207 = 2*MHiggs2(1) + 3*MStau2(2)

	Sf10080 = 2*MHiggs2(1) + 5*MStau2(1)

	Sf10195 = 2*MHiggs2(1) + 3*MStau2(1)

	Sf9145 = 2*MHiggs2(1) + 5*MSbot2(2)

	Sf9275 = 2*MHiggs2(1) + 3*MSbot2(2)

	Sf9141 = 2*MHiggs2(1) + 5*MSbot2(1)

	Sf9263 = 2*MHiggs2(1) + 3*MSbot2(1)

	Sf10178 = -(Sf9302*Sf9303) + 2*MHiggs2(1)

	Sf10159 = -(Sf9290*Sf9291) + 2*MHiggs2(1)

	Sf9244 = -(Sf8885*Sf8886) + 2*MHiggs2(1)

	Sf9224 = -(Sf8873*Sf8874) + 2*MHiggs2(1)

	Sf10365 = -(Sf10245*Sf10246) + 2*MHiggs2(1)

	Sf10338 = -(Sf10225*Sf10226) + 2*MHiggs2(1)

	Sf10398 = MHiggs2(1) + MStop2(4)

	Sf10424 = MHiggs2(1) - MStop2(4)

	Sf10425 = MHiggs2(1) + MStop2(3)

	Sf10426 = MHiggs2(1) - MStop2(3)

	Sf10284 = MHiggs2(1) + MStop2(2)

	Sf10283 = MHiggs2(1) - MStop2(2)

	Sf10294 = MHiggs2(1) + MStop2(1)

	Sf10293 = MHiggs2(1) - MStop2(1)

	Sf10172 = MHiggs2(1) + MStau2(4)

	Sf10199 = MHiggs2(1) - MStau2(4)

	Sf10200 = MHiggs2(1) + MStau2(3)

	Sf10201 = MHiggs2(1) - MStau2(3)

	Sf10096 = MHiggs2(1) + MStau2(2)

	Sf10095 = MHiggs2(1) - MStau2(2)

	Sf10066 = MHiggs2(1) + MStau2(1)

	Sf10065 = MHiggs2(1) - MStau2(1)

	Sf9238 = MHiggs2(1) + MSbot2(4)

	Sf9267 = MHiggs2(1) - MSbot2(4)

	Sf9268 = MHiggs2(1) + MSbot2(3)

	Sf9269 = MHiggs2(1) - MSbot2(3)

	Sf9158 = MHiggs2(1) + MSbot2(2)

	Sf9157 = MHiggs2(1) - MSbot2(2)

	Sf9120 = MHiggs2(1) + MSbot2(1)

	Sf9119 = MHiggs2(1) - MSbot2(1)

	Sf9257 = 2*MZ2 + MHiggs2(1)

	Sf9113 = -2*MZ2 + MHiggs2(1)

	Sf9142 = -3*MZ2 + MHiggs2(1)

	Sf10367 = MHiggs(1) + MStop(2)

	Sf10366 = MHiggs(1) - MStop(2)

	Sf10336 = MHiggs(1) + MStop(1)

	Sf10335 = MHiggs(1) - MStop(1)

	Sf10100 = MHiggs(1) + MStau(2)

	Sf10099 = MHiggs(1) - MStau(2)

	Sf10070 = MHiggs(1) + MStau(1)

	Sf10069 = MHiggs(1) - MStau(1)

	Sf9162 = MHiggs(1) + MSbot(2)

	Sf9161 = MHiggs(1) - MSbot(2)

	Sf9124 = MHiggs(1) + MSbot(1)

	Sf9123 = MHiggs(1) - MSbot(1)

	Sf9132 = MZ + MHiggs(1)

	Sf9131 = -MZ + MHiggs(1)

	Sf121 = 6*CW2*MBot2 - 2*C2B*MW2*SW2

	Sf4768 = 3*CA*MBot2 - CAB*CB*MZ2*SW2

	Sf4725 = 84*MBot2 - 31*MSbot2(2)

	Sf4630 = 84*MBot2 - 31*MSbot2(1)

	Sf3039 = 44*MBot2 + 23*MSbot2(2)

	Sf2912 = 44*MBot2 + 23*MSbot2(1)

	Sf4726 = 10*MBot2 + 109*MSbot2(2)

	Sf3040 = 10*MBot2 + 21*MSbot2(2)

	Sf4631 = 10*MBot2 + 109*MSbot2(1)

	Sf2913 = 10*MBot2 + 21*MSbot2(1)

	Sf2615 = 7*MBot2 - 2*MSbot2(2)

	Sf2388 = 7*MBot2 - 2*MSbot2(1)

	Sf3036 = 5*MBot2 + 78*MSbot2(2)

	Sf4722 = 5*MBot2 - 114*MSbot2(2)

	Sf2909 = 5*MBot2 + 78*MSbot2(1)

	Sf4627 = 5*MBot2 - 114*MSbot2(1)

	Sf3929 = 3*MBot2 + 2*MSbot2(2)

	Sf3805 = 3*MBot2 + 2*MSbot2(1)

	Sf2864 = MBot2 - C2B*MW2

	Sf4257 = CAB*CB - SAB*SB

	Sf1157 = 2*CB - 3*SB*TB

	Sf1163 = CB + 3*SB*TB

	Sf1156 = CAB - SA*SB

	Sf149 = C2B + 4*CW2

	Sf825 = 753 - 60*C2B

	Sf756 = 747 - 132*C2B

	Sf4616 = 585 - 60*C2B

	Sf919 = 543 - 4*C2B

	Sf922 = 531 + 4*C2B

	Sf4614 = 459 + 60*C2B

	Sf832 = 457 + 10*C2B

	Sf4293 = 354 + 9*C2B

	Sf2787 = 351 - 20*C2B

	Sf4296 = 350 - 9*C2B

	Sf829 = 349 - 10*C2B

	Sf2872 = 348 - 99*C2B

	Sf425 = 327 - 10*C2B

	Sf2844 = 303 + 10*C2B

	Sf641 = 296 + 7*C2B

	Sf4344 = 289 + 21*C2B

	Sf826 = 287 + 20*C2B

	Sf643 = 276 - 7*C2B

	Sf376 = 265 + 92*C2B

	Sf2839 = 261 + 44*C2B

	Sf2793 = 257 - 10*C2B

	Sf4346 = 257 - 21*C2B

	Sf335 = 253 + 54*C2B

	Sf799 = 249 - 44*C2B

	Sf4034 = 248 + 54*C2B

	Sf2788 = 247 - 44*C2B

	Sf6613 = 246 - 32*C2B

	Sf4556 = 234 + 32*C2B

	Sf2081 = 233 + 4*C2B

	Sf774 = 233 - 4*C2B

	Sf383 = 231 + 10*C2B

	Sf3780 = 229 + 6*C2B

	Sf5248 = 222 - 92*C2B

	Sf783 = 219 + 38*C2B

	Sf2087 = 219 - 38*C2B

	Sf334 = 217 - 54*C2B

	Sf4545 = 209 + 4*C2B

	Sf6490 = 203 - 4*C2B

	Sf2884 = 203 - 54*C2B

	Sf4356 = 192 - 7*C2B

	Sf830 = 191 - 74*C2B

	Sf4360 = 188 + 7*C2B

	Sf4027 = 181 + 9*C2B

	Sf2835 = 181 - 12*C2B

	Sf4554 = 179 - 38*C2B

	Sf656 = 178 + 3*C2B

	Sf6608 = 178 - 44*C2B

	Sf6496 = 173 + 38*C2B

	Sf416 = 173 + 12*C2B

	Sf372 = 173 - 12*C2B

	Sf6605 = 171 - 90*C2B

	Sf2080 = 167 + 12*C2B

	Sf921 = 167 - 56*C2B

	Sf3777 = 163 - 6*C2B

	Sf753 = 163 - 44*C2B

	Sf775 = 161 - 12*C2B

	Sf5282 = 158 - 4*C2B

	Sf4348 = 157 + 21*C2B

	Sf831 = 155 + 74*C2B

	Sf4589 = 153 + 20*C2B

	Sf3270 = 149 - 9*C2B

	Sf664 = 141 - 7*C2B

	Sf4596 = 139 + 74*C2B

	Sf3190 = 136 + C2B

	Sf916 = 131 + 4*C2B

	Sf4055 = 130 - 27*C2B

	Sf4295 = 128 - 6*C2B

	Sf2077 = 127 + 44*C2B

	Sf666 = 125 + 7*C2B

	Sf4342 = 125 - 21*C2B

	Sf279 = 124 - 33*C2B

	Sf6494 = 123 - 16*C2B

	Sf3192 = 122 - C2B

	Sf3271 = 121 + 9*C2B

	Sf2476 = 119 + 27*C2B

	Sf3774 = 119 + 12*C2B

	Sf642 = 116 - 5*C2B

	Sf2868 = 116 - 141*C2B

	Sf2532 = 115 - 27*C2B

	Sf3319 = 115 - 84*C2B

	Sf2786 = 111 - 46*C2B

	Sf4550 = 109 + 44*C2B

	Sf426 = 109 + 22*C2B

	Sf329 = 109 + 20*C2B

	Sf4598 = 109 - 74*C2B

	Sf4072 = 107 - 24*C2B

	Sf3765 = 104 + 191*C2B

	Sf277 = 104 + 33*C2B

	Sf3768 = 104 - 229*C2B

	Sf1082 = 103 - 2*C2B

	Sf2782 = 103 - 16*C2B

	Sf6485 = 103 - 44*C2B

	Sf3235 = 101 + 63*C2B

	Sf371 = 101 + 16*C2B

	Sf4321 = 101 + 6*C2B

	Sf2527 = 101 + 3*C2B

	Sf3150 = 101 - 9*C2B

	Sf2422 = 101 - 48*C2B

	Sf336 = 97 + 10*C2B

	Sf4030 = 97 + 9*C2B

	Sf2473 = 97 - 3*C2B

	Sf2845 = 97 - 22*C2B

	Sf4549 = 97 - 28*C2B

	Sf4551 = 95 + 22*C2B

	Sf2374 = 94 + 63*C2B

	Sf5280 = 94 + 52*C2B

	Sf3787 = 91 + 108*C2B

	Sf2444 = 91 + 15*C2B

	Sf2446 = 91 + 9*C2B

	Sf773 = 91 + 8*C2B

	Sf3789 = 91 - 116*C2B

	Sf2886 = 89 - 10*C2B

	Sf2510 = 89 - 15*C2B

	Sf6484 = 89 - 22*C2B

	Sf2792 = 87 + 22*C2B

	Sf1100 = 85 - 8*C2B

	Sf382 = 85 - 22*C2B

	Sf3366 = 83 + 84*C2B

	Sf920 = 83 + 56*C2B

	Sf3090 = 82 - C2B

	Sf2836 = 81 + 16*C2B

	Sf2392 = 80 - 33*C2B

	Sf4555 = 79 + 30*C2B

	Sf2834 = 79 - 2*C2B

	Sf3779 = 79 - 6*C2B

	Sf745 = 79 - 9*C2B

	Sf1099 = 79 - 29*C2B

	Sf914 = 77 - 4*C2B

	Sf4358 = 76 + 5*C2B

	Sf3086 = 76 - 6*C2B

	Sf377 = 75 + 4*C2B

	Sf4013 = 74 + 6*C2B

	Sf4292 = 74 - 3*C2B

	Sf368 = 73 + 35*C2B

	Sf3151 = 73 + 9*C2B

	Sf660 = 73 - 7*C2B

	Sf4332 = 73 - 9*C2B

	Sf6495 = 73 - 30*C2B

	Sf2846 = 71 + 28*C2B

	Sf749 = 71 + 9*C2B

	Sf296 = 71 - 16*C2B

	Sf3116 = 71 - 39*C2B

	Sf4289 = 70 + 3*C2B

	Sf2441 = 70 - 3*C2B

	Sf4582 = 69 + 10*C2B

	Sf3795 = 69 - 12*C2B

	Sf2414 = 68 - 6*C2B

	Sf420 = 68 - 11*C2B

	Sf3261 = 68 - 99*C2B

	Sf3360 = 67 + 152*C2B

	Sf1079 = 67 + 2*C2B

	Sf4019 = 67 - 15*C2B

	Sf3759 = 66 + 163*C2B

	Sf3762 = 66 - 229*C2B

	Sf2841 = 65 + 46*C2B

	Sf4329 = 65 + 9*C2B

	Sf4050 = 65 + 6*C2B

	Sf2840 = 65 + 4*C2B

	Sf4026 = 65 - 9*C2B

	Sf3262 = 64 + 99*C2B

	Sf3089 = 64 + C2B

	Sf4035 = 64 - 27*C2B

	Sf421 = 63 - 4*C2B

	Sf378 = 62 + 11*C2B

	Sf2084 = 61 + 8*C2B

	Sf2783 = 61 + 4*C2B

	Sf746 = 61 + 2*C2B

	Sf4600 = 61 - 2*C2B

	Sf422 = 61 - 92*C2B

	Sf3264 = 59 - 21*C2B

	Sf2780 = 59 - 26*C2B

	Sf653 = 58 - C2B

	Sf2832 = 57 + 70*C2B

	Sf659 = 57 + 7*C2B

	Sf4542 = 57 - 2*C2B

	Sf1076 = 57 - 4*C2B

	Sf741 = 57 - 24*C2B

	Sf6618 = 57 - 30*C2B

	Sf3832 = 57 - 104*C2B

	Sf5651 = 57 - 148*C2B

	Sf4032 = 56 + 13*C2B

	Sf4359 = 56 - 5*C2B

	Sf373 = 55 - 2*C2B

	Sf6489 = 55 - 4*C2B

	Sf3316 = 55 - 64*C2B

	Sf823 = 54 + 5*C2B

	Sf822 = 54 - 5*C2B

	Sf3119 = 53 + 39*C2B

	Sf2521 = 53 - 3*C2B

	Sf3863 = 53 - 6*C2B

	Sf998 = 51 + 6*C2B

	Sf3788 = 51 + 4*C2B

	Sf6482 = 51 + 2*C2B

	Sf4354 = 50 - 9*C2B

	Sf257 = 50 - 21*C2B

	Sf4546 = 49 + 4*C2B

	Sf690 = 49 + 2*C2B

	Sf692 = 49 - 2*C2B

	Sf3436 = 49 - 3*C2B

	Sf5770 = 48 + 49*C2B

	Sf5769 = 48 - 11*C2B

	Sf4586 = 48 - 585*C2B

	Sf2831 = 47 + 26*C2B

	Sf2881 = 47 + 10*C2B

	Sf6431 = 47 + 3*C2B

	Sf417 = 47 - 16*C2B

	Sf2875 = 47 - 20*C2B

	Sf2387 = 47 - 27*C2B

	Sf264 = 46 + 21*C2B

	Sf4353 = 46 + 9*C2B

	Sf3838 = 46 - 191*C2B

	Sf4076 = 45 + 8*C2B

	Sf4290 = 44 + 3*C2B

	Sf3228 = 44 - 75*C2B

	Sf293 = 43 + 16*C2B

	Sf4022 = 43 + 15*C2B

	Sf4543 = 43 + 2*C2B

	Sf1057 = 41 + 16*C2B

	Sf2085 = 41 + 10*C2B

	Sf2483 = 41 + 9*C2B

	Sf872 = 41 + 3*C2B

	Sf415 = 41 + C2B

	Sf3111 = 41 - 3*C2B

	Sf412 = 41 - 26*C2B

	Sf2395 = 40 + 39*C2B

	Sf4093 = 39 + 9*C2B

	Sf4591 = 39 - 4*C2B

	Sf6446 = 39 - 33*C2B

	Sf2474 = 39 - 51*C2B

	Sf3229 = 38 + 75*C2B

	Sf1111 = 38 + 13*C2B

	Sf912 = 38 + 3*C2B

	Sf649 = 38 + C2B

	Sf5482 = 38 - 3*C2B

	Sf2442 = 38 - 30*C2B

	Sf6491 = 37 + 4*C2B

	Sf4114 = 37 + 3*C2B

	Sf873 = 37 - 3*C2B

	Sf2504 = 37 - 6*C2B

	Sf3238 = 37 - 7*C2B

	Sf2538 = 37 - 9*C2B

	Sf1054 = 37 - 17*C2B

	Sf892 = 37 - 21*C2B

	Sf3114 = 35 + 3*C2B

	Sf2828 = 35 - 2*C2B

	Sf786 = 35 - 4*C2B

	Sf3250 = 35 - 9*C2B

	Sf981 = 34 + 45*C2B

	Sf369 = 34 + 13*C2B

	Sf2478 = 34 + 3*C2B

	Sf655 = 34 - C2B

	Sf960 = 34 - 87*C2B

	Sf5648 = 34 - 89*C2B

	Sf2777 = 33 + 2*C2B

	Sf650 = 33 + C2B

	Sf675 = 33 - C2B

	Sf1904 = 33 - 2*C2B

	Sf5738 = 33 - 9*C2B

	Sf785 = 33 - 10*C2B

	Sf4328 = 33 - 24*C2B

	Sf2480 = 33 - 39*C2B

	Sf2530 = 32 - 3*C2B

	Sf276 = 32 - 13*C2B

	Sf890 = 31 + 21*C2B

	Sf364 = 31 + 20*C2B

	Sf6454 = 31 + 9*C2B

	Sf886 = 31 + 7*C2B

	Sf3771 = 31 + 6*C2B

	Sf3408 = 31 + 3*C2B

	Sf409 = 31 + 2*C2B

	Sf684 = 31 - 3*C2B

	Sf4070 = 31 - 9*C2B

	Sf870 = 31 - 33*C2B

	Sf1053 = 30 - 11*C2B

	Sf381 = 29 + 7*C2B

	Sf2404 = 29 + 6*C2B

	Sf1075 = 29 + 4*C2B

	Sf3435 = 29 + 3*C2B

	Sf4309 = 29 - C2B

	Sf2508 = 29 - 3*C2B

	Sf887 = 29 - 7*C2B

	Sf3432 = 29 - 9*C2B

	Sf3766 = 28 + 101*C2B

	Sf2396 = 28 + 11*C2B

	Sf4056 = 28 + 9*C2B

	Sf4377 = 28 + 3*C2B

	Sf5831 = 28 - 6*C2B

	Sf3767 = 28 - 139*C2B

	Sf1081 = 27 + 2*C2B

	Sf2776 = 27 - 2*C2B

	Sf3144 = 27 - 5*C2B

	Sf915 = 27 - 26*C2B

	Sf2370 = 26 + 21*C2B

	Sf255 = 26 + 7*C2B

	Sf646 = 26 + 3*C2B

	Sf651 = 26 - C2B

	Sf661 = 25 + 7*C2B

	Sf531 = 25 + 3*C2B

	Sf2829 = 25 + 2*C2B

	Sf366 = 25 - 2*C2B

	Sf662 = 25 - 7*C2B

	Sf3122 = 25 - 21*C2B

	Sf2071 = 24 + C2B

	Sf979 = 23 + 87*C2B

	Sf3330 = 23 + 84*C2B

	Sf5674 = 23 + 58*C2B

	Sf3346 = 23 + 28*C2B

	Sf2389 = 23 + 27*C2B

	Sf330 = 23 + 20*C2B

	Sf4104 = 23 + 9*C2B

	Sf682 = 23 + 8*C2B

	Sf508 = 23 + 7*C2B

	Sf2466 = 23 + 6*C2B

	Sf5664 = 23 + 2*C2B

	Sf3186 = 23 + C2B

	Sf4045 = 23 - 2*C2B

	Sf1024 = 23 - 3*C2B

	Sf3773 = 23 - 4*C2B

	Sf3288 = 23 - 12*C2B

	Sf2779 = 23 - 14*C2B

	Sf2775 = 23 - 20*C2B

	Sf3298 = 23 - 28*C2B

	Sf3311 = 23 - 84*C2B

	Sf3755 = 22 + 69*C2B

	Sf2073 = 22 + 15*C2B

	Sf5724 = 22 + 3*C2B

	Sf645 = 22 - 3*C2B

	Sf2454 = 22 - 6*C2B

	Sf4054 = 22 - 9*C2B

	Sf2373 = 22 - 21*C2B

	Sf3756 = 22 - 119*C2B

	Sf1058 = 21 + 8*C2B

	Sf607 = 21 + C2B

	Sf652 = 21 - C2B

	Sf505 = 21 - 13*C2B

	Sf3308 = 21 - 68*C2B

	Sf2394 = 20 + 39*C2B

	Sf2525 = 20 + 21*C2B

	Sf1126 = 20 + 3*C2B

	Sf1116 = 20 - 3*C2B

	Sf2536 = 20 - 18*C2B

	Sf3108 = 20 - 27*C2B

	Sf2425 = 19 + 24*C2B

	Sf2514 = 19 + 21*C2B

	Sf5689 = 19 + 12*C2B

	Sf1061 = 19 + 5*C2B

	Sf3778 = 19 + 2*C2B

	Sf2467 = 19 + C2B

	Sf1106 = 19 - 3*C2B

	Sf6423 = 19 - 6*C2B

	Sf878 = 19 - 7*C2B

	Sf6480 = 19 - 10*C2B

	Sf5635 = 19 - 12*C2B

	Sf427 = 19 - 14*C2B

	Sf950 = 19 - 15*C2B

	Sf961 = 19 - 53*C2B

	Sf3333 = 19 - 100*C2B

	Sf4298 = 18 + C2B

	Sf2372 = 18 - 7*C2B

	Sf884 = 18 - 25*C2B

	Sf5677 = 17 + 18*C2B

	Sf2452 = 17 + 12*C2B

	Sf272 = 17 + 9*C2B

	Sf2033 = 17 + 7*C2B

	Sf966 = 17 + 2*C2B

	Sf2235 = 17 - 6*C2B

	Sf4336 = 17 - 7*C2B

	Sf3231 = 17 - 9*C2B

	Sf3242 = 17 - 21*C2B

	Sf4583 = 16 + 153*C2B

	Sf3138 = 16 + 3*C2B

	Sf640 = 16 + C2B

	Sf764 = 16 - C2B

	Sf1023 = 16 - 3*C2B

	Sf2472 = 16 - 21*C2B

	Sf503 = 15 + 13*C2B

	Sf876 = 15 + 7*C2B

	Sf2076 = 15 + 4*C2B

	Sf523 = 15 - C2B

	Sf1080 = 15 - 2*C2B

	Sf2791 = 15 - 4*C2B

	Sf256 = 15 - 7*C2B

	Sf1052 = 15 - 9*C2B

	Sf3295 = 15 - 26*C2B

	Sf3109 = 14 + 27*C2B

	Sf2485 = 14 + 9*C2B

	Sf2230 = 14 + 3*C2B

	Sf1129 = 14 - 3*C2B

	Sf206 = 14 - 8*C2B

	Sf2377 = 14 - 9*C2B

	Sf2385 = 14 - 15*C2B

	Sf3363 = 13 + 68*C2B

	Sf3823 = 13 + 56*C2B

	Sf2031 = 13 + 11*C2B

	Sf2870 = 13 + 10*C2B

	Sf4401 = 13 + 9*C2B

	Sf526 = 13 + 5*C2B

	Sf6486 = 13 + 4*C2B

	Sf603 = 13 + 3*C2B

	Sf2450 = 13 + 2*C2B

	Sf301 = 13 + C2B

	Sf3191 = 13 - 7*C2B

	Sf292 = 13 - 8*C2B

	Sf1049 = 13 - 9*C2B

	Sf2439 = 13 - 12*C2B

	Sf2459 = 13 - 21*C2B

	Sf261 = 12 + 7*C2B

	Sf931 = 11 + 108*C2B

	Sf3826 = 11 + 18*C2B

	Sf3849 = 11 + 12*C2B

	Sf333 = 11 + 10*C2B

	Sf5749 = 11 + 9*C2B

	Sf4540 = 11 + 8*C2B

	Sf2462 = 11 + 6*C2B

	Sf2383 = 11 + 3*C2B

	Sf2074 = 11 + 2*C2B

	Sf365 = 11 + C2B

	Sf588 = 11 - C2B

	Sf1903 = 11 - 2*C2B

	Sf765 = 11 - 3*C2B

	Sf680 = 11 - 8*C2B

	Sf5741 = 11 - 9*C2B

	Sf2878 = 11 - 20*C2B

	Sf2379 = 10 + 21*C2B

	Sf1095 = 10 - C2B

	Sf4291 = 10 - 3*C2B

	Sf501 = 10 - 9*C2B

	Sf968 = 10 - 57*C2B

	Sf5659 = 10 - 87*C2B

	Sf3187 = 9 + 26*C2B

	Sf2827 = 9 + 20*C2B

	Sf4402 = 9 + 11*C2B

	Sf284 = 9 + 2*C2B

	Sf300 = 9 + C2B

	Sf499 = 9 - C2B

	Sf285 = 9 - 2*C2B

	Sf2086 = 9 - 4*C2B

	Sf288 = 9 - 5*C2B

	Sf2350 = 9 - 6*C2B

	Sf5640 = 9 - 13*C2B

	Sf8659 = 8 - 32*SW2

	Sf989 = 8 + 41*C2B

	Sf883 = 8 + 25*C2B

	Sf2513 = 8 + 9*C2B

	Sf3400 = 8 + 6*C2B

	Sf4112 = 8 + 3*C2B

	Sf2023 = 8 + C2B

	Sf410 = 8 - C2B

	Sf3084 = 8 - 3*C2B

	Sf5744 = 8 - 6*C2B

	Sf3423 = 8 - 9*C2B

	Sf6445 = 8 - 15*C2B

	Sf2376 = 8 - 21*C2B

	Sf3146 = 8 - 30*C2B

	Sf2407 = 7 + 15*C2B

	Sf4325 = 7 + 8*C2B

	Sf4376 = 7 + 6*C2B

	Sf2044 = 7 + 5*C2B

	Sf784 = 7 + 4*C2B

	Sf712 = 7 + 3*C2B

	Sf2882 = 7 + 2*C2B

	Sf299 = 7 + C2B

	Sf532 = 7 - C2B

	Sf2867 = 7 - 2*C2B

	Sf1060 = 7 - 3*C2B

	Sf209 = 7 - 4*C2B

	Sf186 = 7 - 5*C2B

	Sf859 = 7 - 6*C2B

	Sf4106 = 7 - 9*C2B

	Sf2424 = 7 - 24*C2B

	Sf3761 = 6 + 79*C2B

	Sf5696 = 6 + 19*C2B

	Sf270 = 6 + 5*C2B

	Sf326 = 6 + C2B

	Sf703 = 6 - C2B

	Sf327 = 6 - 5*C2B

	Sf262 = 6 - 7*C2B

	Sf2457 = 6 - 9*C2B

	Sf3233 = 6 - 84*C2B

	Sf7297 = 5 + Ss78**2

	Sf3860 = 5 + 16*C2B

	Sf3843 = 5 + 12*C2B

	Sf3429 = 5 + 9*C2B

	Sf2041 = 5 + 7*C2B

	Sf911 = 5 + 6*C2B

	Sf2499 = 5 + 4*C2B

	Sf2356 = 5 + 3*C2B

	Sf1930 = 5 + 2*C2B

	Sf759 = 5 + C2B

	Sf522 = 5 - C2B

	Sf3184 = 5 - 2*C2B

	Sf1050 = 5 - 3*C2B

	Sf190 = 5 - 4*C2B

	Sf207 = 5 - 6*C2B

	Sf713 = 5 - 7*C2B

	Sf6479 = 5 - 8*C2B

	Sf5645 = 5 - 12*C2B

	Sf3305 = 5 - 16*C2B

	Sf970 = 5 - 18*C2B

	Sf5656 = 5 - 36*C2B

	Sf3327 = 5 - 61*C2B

	Sf4399 = 4 + 15*C2B

	Sf2417 = 4 + 6*C2B

	Sf1950 = 4 + 5*C2B

	Sf2371 = 4 + 3*C2B

	Sf438 = 4 + C2B

	Sf234 = 4 - C2B

	Sf1047 = 4 - 3*C2B

	Sf956 = 4 - 5*C2B

	Sf4404 = 4 - 6*C2B

	Sf747 = 4 - 9*C2B

	Sf278 = 4 - 13*C2B

	Sf868 = 4 - 15*C2B

	Sf3266 = 4 - 63*C2B

	Sf7222 = 3 - UStop22**2

	Sf1776 = 3 - 4*SW2*UStau22**2

	Sf1779 = 3 - UStau22**2

	Sf4798 = 3 - USbot22**2

	Sf227 = 3 - Log(MSneu2(1))

	Sf3820 = 3 + 25*C2B

	Sf3749 = 3 + 17*C2B

	Sf3846 = 3 + 11*C2B

	Sf3340 = 3 + 8*C2B

	Sf171 = 3 + 7*C2B

	Sf991 = 3 + 5*C2B

	Sf168 = 3 + 4*C2B

	Sf720 = 3 + 2*C2B

	Sf118 = 3 + C2B

	Sf181 = 3 - C2B

	Sf704 = 3 - 2*C2B

	Sf3253 = 3 - 4*C2B

	Sf954 = 3 - 5*C2B

	Sf2348 = 3 - 6*C2B

	Sf3357 = 3 - 10*C2B

	Sf3750 = 3 - 31*C2B

	Sf3343 = 3 - 32*C2B

	Sf1893 = 2 - Log(MSneu2(1))

	Sf2174 = (-1 + 2*CW2)/CW2

	Sf972 = 2 + 15*C2B

	Sf858 = 2 + 6*C2B

	Sf602 = 2 + 3*C2B

	Sf243 = 2 + C2B

	Sf2 = 2 - C2B

	Sf1091 = 2 - 3*C2B

	Sf528 = 2 - 5*C2B

	Sf3422 = 2 - 6*C2B

	Sf413 = 2 - 7*C2B

	Sf6434 = 2 - 9*C2B

	Sf3760 = 2 - 19*C2B

	Sf877 = 2 - 21*C2B

	Sf7370 = 1 - UStop22**2

	Sf7313 = 1 - UStau22**2

	Sf7401 = -Sf8583

	Sf2022 = 1 - Log(MSneu2(1))

	Sf8715 = 1 - Log(MM2)

	Sf4772 = 1 + 2*CW2

	Sf1754 = 1 - 2*CW2

	Sf7198 = 1 - 4*CW2

	Sf3260 = 1 + 45*C2B

	Sf2410 = 1 + 15*C2B

	Sf199 = 1 + 13*C2B

	Sf294 = 1 + 8*C2B

	Sf721 = 1 + 7*C2B

	Sf4394 = 1 + 6*C2B

	Sf197 = 1 + 5*C2B

	Sf569 = 1 + 4*C2B

	Sf169 = 1 + 3*C2B

	Sf566 = 1 + 2*C2B

	Sf174 = 1 + C2B

	Sf175 = 1 - C2B

	Sf184 = 1 - 2*C2B

	Sf601 = 1 - 3*C2B

	Sf195 = 1 - 4*C2B

	Sf627 = 1 - 5*C2B

	Sf504 = 1 - 6*C2B

	Sf861 = 1 - 9*C2B

	Sf710 = 1 - 11*C2B

	Sf3324 = 1 - 12*C2B

	Sf891 = 1 - 14*C2B

	Sf408 = 1 - 20*C2B

	Sf3240 = 1 - 21*C2B

	Sf8633 = Sf7314*MSneu2(1) - 
     &   4*MStau2(2)*UStau22**2*(-1 + UStau22**2)

	Sf8630 = 502 + UStau22**2*(60 + 12*UStau22**2)

	Sf9774 = Sf9770*MStau2(1)**4 - 
     &   8*Sf9363**3*Sf9364**3*MHiggs2(2)*MStau2(2)**3 + 
     &   MStau2(1)*(Sf9771*MStau2(2)**2 + 
     &      MStau2(1)*(Sf9773*MStau2(2) + 
     &         MStau2(1)*(MStau2(1)**2*
     &             (-2*Sf9769 + 
     &               MStau2(1)*(-4*Sf9768 + 16*MStau2(1))) - 
     &            Sf9772*MStau2(2))))

	Sf9747 = Sf9745*MStau2(1)**4 - 
     &   2*Sf9363**3*Sf9364**3*MHiggs2(2)*MStau2(2)**3 + 
     &   MStau2(1)*(Sf9363**2*Sf9364**2*Sf9742*MStau2(2)**2 + 
     &      MStau2(1)*(-(Sf9363*Sf9364*Sf9744*MStau2(2)) + 
     &         MStau2(1)*(MStau2(1)**2*
     &             (-Sf9743 + (Sf9741 - 2*MStau2(1))*MStau2(1)) - 
     &            Sf9746*MStau2(2))))

	Sf9718 = Sf9714*MStau2(1)**4 - 
     &   32*Sf9363**3*Sf9364**3*MHiggs2(2)*MStau2(2)**3 + 
     &   MStau2(1)*(Sf9715*MStau2(2)**2 + 
     &      MStau2(1)*(Sf9717*MStau2(2) + 
     &         MStau2(1)*(MStau2(1)**2*
     &             (-Sf9713 + MStau2(1)*(-4*Sf9712 + 64*MStau2(1)))
     &              - Sf9716*MStau2(2))))

	Sf9692 = 2*MHiggs2(2)**4*
     &    (MStau2(2)**3 + MStau2(1)*
     &       (MStau2(1)*(MStau2(1) - 2*MStau2(2)) + 
     &         10*MStau2(2)**2)) - 
     &   2*MHiggs2(2)**3*MStau2(1)*
     &    (21*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(3*MStau2(1) - 9*MStau2(2)) + 
     &         17*MStau2(2)**2)) + 
     &   2*MHiggs2(2)**2*(-MStau2(2)**5 + 
     &      MStau2(1)*(16*MStau2(2)**4 + 
     &         MStau2(1)*(9*MStau2(2)**3 + 
     &            MStau2(1)*
     &             (MStau2(1)*(MStau2(1) - 12*MStau2(2)) + 
     &               19*MStau2(2)**2)))) + 
     &   MHiggs2(2)*(MStau2(2)**4 + 
     &      MStau2(1)*(-11*MStau2(2)**3 + 
     &         MStau2(1)*(-5*MStau2(2)**2 + 
     &            MStau2(1)*(4*MStau2(1) + 43*MStau2(2)))))*
     &    MStau2(3)**2 - 2*MStau2(1)*
     &    (-2*MStau2(2)**2 + MStau2(1)*(MStau2(1) + 5*MStau2(2)))*
     &    MStau2(3)**4 - MHiggs2(2)**5*MStau2(2)*MStau2(4)

	Sf9824 = Sf9822*MStau2(1)**2 + 
     &   4*Sf9360*Sf9363**3*Sf9364**3*MHiggs2(2)*MStau2(2)**2 + 
     &   MStau2(1)*(MStau2(1)**2*
     &       (-Sf9821 + MStau2(1)*
     &          (Sf9820 + MStau2(1)*
     &             (Sf9819 + MStau2(1)*(-10*Sf9818 + 2*MStau2(1))))
     &         ) + Sf9823*MStau2(2))

	Sf9803 = Sf9801*MStau2(1)**2 + 
     &   16*Sf9360*Sf9363**3*Sf9364**3*MHiggs2(2)*MStau2(2)**2 + 
     &   MStau2(1)*(MStau2(1)**2*
     &       (-Sf9800 + MStau2(1)*
     &          (Sf9799 + MStau2(1)*
     &             (-Sf9798 + (-2*Sf9797 - 10*MStau2(1))*MStau2(1))
     &            )) + Sf9802*MStau2(2))

	Sf10051 = Sf10049*MStau2(1)**2 + 
     &   4*Sf10043*Sf9363**3*Sf9364**3*MHiggs2(2)*MStau2(2)**3 + 
     &   MStau2(1)*(Sf10050*MStau2(2) + 
     &      MStau2(1)**2*(MStau2(1)*
     &          (-2*Sf10047 + 
     &            MStau2(1)*
     &             (Sf10046 + 
     &               MStau2(1)*(-Sf10045 + 4*Sf10044*MStau2(1))))+
     &           2*Sf10048*MStau2(2)))

	Sf10015 = Sf10013*MStau2(1)**2 + 
     &   16*Sf10008*Sf9363**3*Sf9364**3*MHiggs2(2)*MStau2(2)**3 + 
     &   MStau2(1)*(MStau2(1)**2*
     &       (-2*Sf10012 + 
     &         MStau2(1)*(2*Sf10011 + 
     &            MStau2(1)*
     &             (-Sf10010 + 
     &               MStau2(1)*(Sf10009 - 16*Sf10007*MStau2(1)))))-
     &        Sf10014*MStau2(2))

	Sf10001 = Sf10000*MStau2(1)**2 - 
     &   2*Sf9363**3*Sf9364**3*Sf9686*MHiggs2(2)*MStau2(2)**3 + 
     &   MStau2(1)*(MStau2(1)**2*
     &       (-2*Sf9999 + MStau2(1)*
     &          (2*Sf9997 + 
     &            MStau2(1)*
     &             (-Sf9996 + 
     &               MStau2(1)*(Sf9995 - 2*Sf9994*MStau2(1))))) + 
     &      Sf9363*Sf9364*Sf9998*MStau2(2))

	Sf5979 = MW2*(3*Sf5978 + MW2*(MW2 - 3*Sf99)) - 
     &   MStop2(2)**3 + MSbot2(2)*
     &    (42*MStop2(2)**2 + MSbot2(2)*(-MSbot2(2) + 42*MStop2(2)))

	Sf5923 = MW2*(3*Sf5922 + MW2*(MW2 - 3*Sf99)) - 
     &   MStop2(2)**3 + MSbot2(2)*
     &    (MSbot2(2)*(-MSbot2(2) - 6*MStop2(2)) - 6*MStop2(2)**2)

	Sf5553 = MW2*(3*Sf5552 + MW2*(MW2 - 3*Sf91)) - 
     &   MStop2(2)**3 + MSbot2(1)*
     &    (42*MStop2(2)**2 + MSbot2(1)*(-MSbot2(1) + 42*MStop2(2)))

	Sf5490 = MW2*(3*Sf5489 + MW2*(MW2 - 3*Sf91)) - 
     &   MStop2(2)**3 + MSbot2(1)*
     &    (MSbot2(1)*(-MSbot2(1) - 6*MStop2(2)) - 6*MStop2(2)**2)

	Sf3554 = MW2*(3*Sf3553 + MW2*(MW2 - 3*Sf77)) - 
     &   MStop2(1)**3 + MSbot2(2)*
     &    (42*MStop2(1)**2 + MSbot2(2)*(-MSbot2(2) + 42*MStop2(1)))

	Sf3479 = MW2*(3*Sf3478 + MW2*(MW2 - 3*Sf77)) - 
     &   MStop2(1)**3 + MSbot2(2)*
     &    (MSbot2(2)*(-MSbot2(2) - 6*MStop2(1)) - 6*MStop2(1)**2)

	Sf483 = MW2*(MW2*(MW2 - 3*Sf34) + 3*Sf482) - 
     &   MStau2(1)**3 + MSneu2(1)*
     &    (MSneu2(1)*(-MSneu2(1) - 6*MStau2(1)) - 6*MStau2(1)**2)

	Sf3198 = MW2*(3*Sf3197 + MW2*(MW2 - 3*Sf69)) - 
     &   MStop2(1)**3 + MSbot2(1)*
     &    (42*MStop2(1)**2 + MSbot2(1)*(-MSbot2(1) + 42*MStop2(1)))

	Sf3096 = MW2*(3*Sf3095 + MW2*(MW2 - 3*Sf69)) - 
     &   MStop2(1)**3 + MSbot2(1)*
     &    (MSbot2(1)*(-MSbot2(1) - 6*MStop2(1)) - 6*MStop2(1)**2)

	Sf1400 = MW2*(3*Sf1399 + MW2*(MW2 - 3*Sf47)) - 
     &   MStau2(2)**3 + MSneu2(1)*
     &    (MSneu2(1)*(-MSneu2(1) - 6*MStau2(2)) - 6*MStau2(2)**2)

	Sf9986 = MZ2**4*(MStau2(1)*(2*MStau2(1) - 5*MStau2(2)) + 
     &      5*MStau2(2)**2) - 
     &   MZ2**3*(9*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(17*MStau2(1) - 41*MStau2(2)) + 
     &         31*MStau2(2)**2)) + 
     &   MZ2**2*(7*MStau2(2)**4 + 
     &      MStau2(1)*(33*MStau2(2)**3 + 
     &         MStau2(1)*(MStau2(1)*
     &             (22*MStau2(1) - 45*MStau2(2)) + 15*MStau2(2)**2)
     &         )) - MZ2*(7*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(7*MStau2(1) - 41*MStau2(2)) + 
     &         11*MStau2(2)**2))*MStau2(3)**2 - 
     &   4*MStau2(2)*(MStau2(1) - MStau2(3))*MStau2(3)**4

	Sf10042 = MZ2**4*
     &    (29*MStau2(2)**2 + MStau2(1)*(8*MStau2(1) + MStau2(2)))-
     &     2*MZ2**3*(33*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(34*MStau2(1) - 43*MStau2(2)) + 
     &         128*MStau2(2)**2)) + 
     &   MZ2**2*(73*MStau2(2)**4 + 
     &      MStau2(1)*(345*MStau2(2)**3 + 
     &         MStau2(1)*(51*MStau2(2)**2 + 
     &            MStau2(1)*(16*MStau2(1) + 123*MStau2(2))))) + 
     &   2*MZ2*(-32*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(22*MStau2(1) - 17*MStau2(2)) - 
     &         13*MStau2(2)**2))*MStau2(3)**2 - 
     &   4*MStau2(2)*(-5*MStau2(2) - 2*MStau2(3))*MStau2(3)**4

	Sf9540 = Sf9535*MStau2(1) - Sf9536*MStau2(1)*MStau2(2) + 
     &   MStau2(2)**2*(Sf9537 + 
     &      MStau2(2)*(-2*Sf9538 + 
     &         MStau2(2)*(-4*Sf9539 + 16*MStau2(2))))

	Sf9427 = Sf9422*MStau2(1) - Sf9423*MStau2(1)*MStau2(2) + 
     &   MStau2(2)**2*(Sf9424 + 
     &      MStau2(2)*(-Sf9425 + 
     &         MStau2(2)*(-4*Sf9426 + 64*MStau2(2))))

	Sf9411 = Sf9406*MStau2(1) - Sf9407*MStau2(2) + 
     &   MStau2(2)**2*(Sf9408 + 
     &      MStau2(2)*(-Sf9409 + MStau2(2)*(Sf9410 + 2*MStau2(2))))

	Sf9395 = 2*Sf9390*MStau2(1) - 2*Sf9391*MStau2(2) + 
     &   MStau2(2)**2*(Sf9392 + 
     &      MStau2(2)*(-Sf9393 + 
     &         MStau2(2)*(-2*Sf9394 + 2*MStau2(2))))

	Sf9972 = Sf9966 + 2*Sf9967*MStau2(1)*MStau2(2) + 
     &   MStau2(2)**2*(-2*Sf9968 + 
     &      MStau2(2)*(Sf9969 + 
     &         MStau2(2)*(-Sf9970 + 4*Sf9971*MStau2(2))))

	Sf9924 = Sf9918 - Sf9919*MStau2(2) + 
     &   MStau2(2)**2*(2*Sf9920 + 
     &      MStau2(2)*(-Sf9921 + 
     &         MStau2(2)*(Sf9922 - 16*Sf9923*MStau2(2))))

	Sf9914 = Sf9908 - 2*Sf9909*MStau2(2) + 
     &   MStau2(2)**2*(2*Sf9910 + 
     &      MStau2(2)*(-Sf9911 + 
     &         MStau2(2)*(Sf9912 - 16*Sf9913*MStau2(2))))

	Sf9898 = Sf9892 - 2*Sf9893*MStau2(2) + 
     &   MStau2(2)**2*(2*Sf9894 + 
     &      MStau2(2)*(-Sf9895 + 
     &         MStau2(2)*(Sf9896 - 2*Sf9897*MStau2(2))))

	Sf9675 = Sf9670 - Sf9671*MStau2(2) + 
     &   MStau2(2)**2*(Sf9672 + 
     &      MStau2(2)*(-Sf9673 + 
     &         MStau2(2)*(-16*Sf9674 + 44*MStau2(2))))

	Sf6101 = 6*MTop2*MW2**2*Sf4200 - 72*MSbot2(2)**4 + 
     &   MW2*(116*MSbot2(2)**3 + 
     &      MW2*(-33*MSbot2(2)**2 + 
     &         MW2*(-(MW2*Sf504) + 2*Sf169*MSbot2(2))))

	Sf5712 = 6*MTop2*MW2**2*Sf4120 - 72*MSbot2(1)**4 + 
     &   MW2*(116*MSbot2(1)**3 + 
     &      MW2*(-33*MSbot2(1)**2 + 
     &         MW2*(-(MW2*Sf504) + 2*Sf169*MSbot2(1))))

	Sf6099 = 12*MTop2*MW2**3 + 180*MSbot2(2)**4 + 
     &   MW2*(-152*MSbot2(2)**3 + 
     &      MW2*(27*MSbot2(2)**2 + MW2*(MW2*Sf5689 + 26*MSbot2(2)))
     &      )

	Sf5710 = 12*MTop2*MW2**3 + 180*MSbot2(1)**4 + 
     &   MW2*(-152*MSbot2(1)**3 + 
     &      MW2*(27*MSbot2(1)**2 + MW2*(MW2*Sf5689 + 26*MSbot2(1)))
     &      )

	Sf1813 = 4*ML2*MW2**3 - 75*MStau2(2)**4 + 
     &   MW2*(82*MStau2(2)**3 + 
     &      MW2*(MW2*(-4*MW2*Sf174 - 4*MStau2(2)) - 
     &         11*MStau2(2)**2))

	Sf1002 = 4*ML2*MW2**3 - 75*MStau2(1)**4 + 
     &   MW2*(82*MStau2(1)**3 + 
     &      MW2*(MW2*(-4*MW2*Sf174 - 4*MStau2(1)) - 
     &         11*MStau2(1)**2))

	tmp6 = -86*MZ2**2 + MHiggs2(2)*(-108*MZ2 + 2*MHiggs2(2))

	Sf9777 = 85*MZ2**4 + 
     &   MHiggs2(2)*(147*MZ2**3 + tmp6*MHiggs2(2)) + 
     &   Sf9776*MStau2(2)

	Sf9580 = 85*MZ2**4 + 
     &   MHiggs2(2)*(147*MZ2**3 + tmp6*MHiggs2(2)) + 
     &   Sf9579*MStau2(1)

	Sf9487 = Sf9483*MStau2(1) - Sf9484*MStau2(2) + 
     &   MStau2(2)**2*(Sf9485 + MStau2(2)*(-Sf9486 + 2*MStau2(2)))

	Sf6917 = 10*MBot2*MSbot2(2)**2 + 
     &   (Sf4726 - 305*MStop2(2))*MStop2(2)**2 + 
     &   MSbot2(2)*(227*MSbot2(2)**2 + Sf4725*MStop2(2))

	Sf4727 = 10*MBot2*MSbot2(2)**2 + 
     &   (Sf4726 - 305*MStop2(1))*MStop2(1)**2 + 
     &   MSbot2(2)*(227*MSbot2(2)**2 + Sf4725*MStop2(1))

	Sf6666 = 10*MBot2*MSbot2(1)**2 + 
     &   (Sf4631 - 305*MStop2(2))*MStop2(2)**2 + 
     &   MSbot2(1)*(227*MSbot2(1)**2 + Sf4630*MStop2(2))

	Sf4632 = 10*MBot2*MSbot2(1)**2 + 
     &   (Sf4631 - 305*MStop2(1))*MStop2(1)**2 + 
     &   MSbot2(1)*(227*MSbot2(1)**2 + Sf4630*MStop2(1))

	Sf6161 = 6*MTop2*Sf99**2 - 229*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(191*MSbot2(2) - 11*MStop2(2)) + 
     &      49*MStop2(2)**2)

	Sf6177 = 6*MBot2*Sf99**2 + 229*MStop2(2)**3 + 
     &   MSbot2(2)*(-49*MStop2(2)**2 + 
     &      MSbot2(2)*(-191*MSbot2(2) + 11*MStop2(2)))

	Sf9645 = Sf9641 - Sf9642*MStau2(2) + 
     &   MStau2(2)**2*(Sf9643 + 
     &      (4*Sf9644 - 20*MStau2(2))*MStau2(2))

	Sf5782 = 6*MTop2*Sf91**2 - 229*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(191*MSbot2(1) - 11*MStop2(2)) + 
     &      49*MStop2(2)**2)

	Sf5799 = 6*MBot2*Sf91**2 + 229*MStop2(2)**3 + 
     &   MSbot2(1)*(-49*MStop2(2)**2 + 
     &      MSbot2(1)*(-191*MSbot2(1) + 11*MStop2(2)))

	Sf3913 = 6*MTop2*Sf77**2 - 229*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(191*MSbot2(2) - 11*MStop2(1)) + 
     &      49*MStop2(1)**2)

	Sf3932 = 6*MBot2*Sf77**2 + 229*MStop2(1)**3 + 
     &   MSbot2(2)*(-49*MStop2(1)**2 + 
     &      MSbot2(2)*(-191*MSbot2(2) + 11*MStop2(1)))

	Sf3785 = 6*MTop2*Sf69**2 - 229*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(191*MSbot2(1) - 11*MStop2(1)) + 
     &      49*MStop2(1)**2)

	Sf3808 = 6*MBot2*Sf69**2 + 229*MStop2(1)**3 + 
     &   MSbot2(1)*(-49*MStop2(1)**2 + 
     &      MSbot2(1)*(-191*MSbot2(1) + 11*MStop2(1)))

	Sf6899 = 6*MBot2*Sf6898 - 915*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(681*MSbot2(2) - 93*MStop2(2)) + 
     &      327*MStop2(2)**2)

	Sf6648 = 6*MBot2*Sf6647 - 915*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(681*MSbot2(1) - 93*MStop2(2)) + 
     &      327*MStop2(2)**2)

	Sf5418 = 2*MBot2*Sf5417 - 121*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(123*MSbot2(2) - 23*MStop2(2)) + 
     &      21*MStop2(2)**2)

	Sf5312 = 2*MBot2*Sf5311 - 121*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(123*MSbot2(1) - 23*MStop2(2)) + 
     &      21*MStop2(2)**2)

	Sf4703 = 6*MBot2*Sf4702 - 915*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(681*MSbot2(2) - 93*MStop2(1)) + 
     &      327*MStop2(1)**2)

	Sf4606 = 6*MBot2*Sf4605 - 915*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(681*MSbot2(1) - 93*MStop2(1)) + 
     &      327*MStop2(1)**2)

	Sf3018 = 2*MBot2*Sf3017 - 121*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(123*MSbot2(2) - 23*MStop2(1)) + 
     &      21*MStop2(1)**2)

	Sf2891 = 2*MBot2*Sf2890 - 121*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(123*MSbot2(1) - 23*MStop2(1)) + 
     &      21*MStop2(1)**2)

	Sf10031 = 2*MZ2**3*MHiggs2(2) - Sf9943*MStau2(2)**3 + 
     &   MZ2*(-4*Sf9942*MStau2(2)**2 + 
     &      MZ2*(-2*MZ2**2 + 7*Sf9941*MStau2(2)))

	Sf9944 = 2*MZ2**3*MHiggs2(2) - Sf9943*MStau2(1)**3 + 
     &   MZ2*(-4*Sf9942*MStau2(1)**2 + 
     &      MZ2*(-2*MZ2**2 + 7*Sf9941*MStau2(1)))

	Sf910 = 27*C2B*MW2**5*Sf907 + 203*MSl2Diff(1,1)**6 + 
     &   MW2**2*(MW2**2*(230*MW2**2 + Sf908) + 
     &      2*Sf909*MSl2Diff(1,1)**2)

	Sf1709 = 27*C2B*MW2**5*Sf1706 + 203*MSl2Diff(2,1)**6 + 
     &   MW2**2*(MW2**2*(230*MW2**2 + Sf1707) + 
     &      2*Sf1708*MSl2Diff(2,1)**2)

	Sf481 = 9*C2B*MW2**5*Sf478 + 89*MSl2Diff(1,1)**6 + 
     &   MW2**2*(MW2**2*(2*MW2**2 + 5*Sf479) + 
     &      2*Sf480*MSl2Diff(1,1)**2)

	Sf1398 = 9*C2B*MW2**5*Sf1395 + 89*MSl2Diff(2,1)**6 + 
     &   MW2**2*(MW2**2*(2*MW2**2 + 5*Sf1396) + 
     &      2*Sf1397*MSl2Diff(2,1)**2)

	Sf2684 = 3*MW2**3*MHiggs2(4) + Sf173*MSbot2(2)**3 + 
     &   MW2*(-(Sf2491*MSbot2(2)**2) + MW2**2*(-4*MW2 + MSbot2(2)))

	Sf2492 = 3*MW2**3*MHiggs2(4) + Sf173*MSbot2(1)**3 + 
     &   MW2*(-(Sf2491*MSbot2(1)**2) + MW2**2*(-4*MW2 + MSbot2(1)))

	Sf6097 = 2*MTop2*MW2**2 + 80*MSbot2(2)**3 + 
     &   MW2*(24*MSbot2(2)**2 + MW2*(MW2*Sf5664 + 29*MSbot2(2)))

	Sf5708 = 2*MTop2*MW2**2 + 80*MSbot2(1)**3 + 
     &   MW2*(24*MSbot2(1)**2 + MW2*(MW2*Sf5664 + 29*MSbot2(1)))

	Sf1810 = 6*ML2*MW2**2 - 300*MStau2(2)**3 + 
     &   MW2*(MW2*(-(MW2*Sf998) - 57*MStau2(2)) - 36*MStau2(2)**2)

	Sf999 = 6*ML2*MW2**2 - 300*MStau2(1)**3 + 
     &   MW2*(MW2*(-(MW2*Sf998) - 57*MStau2(1)) - 36*MStau2(1)**2)

	tmp5 = MW2**2*(-7*MW2 + 3*MHiggs2(4))

	Sf1216 = tmp5*MHiggs2(4) - Sf216*MStau2(2)**2 + 
     &   MW2**2*(4*MW2**2 + 2*Ss241*Ss242*MStau2(2))

	Sf224 = tmp5*MHiggs2(4) - Sf216*MStau2(1)**2 + 
     &   MW2**2*(4*MW2**2 + 2*Ss241*Ss242*MStau2(1))

	Sf1678 = Sf1399 + MW2*(MW2 - 2*Sf47)

	Sf865 = MW2*(MW2 - 2*Sf34) + Sf482

	tmp8 = 4*MW2**2 + MHiggs2(4)*(-4*Sf1039 + 6*MHiggs2(4))

	Sf1042 = Sf1041*Sf559 + tmp8*MHiggs2(4)**2 - 
     &   4*Sf1040*Sf237*MStau2(1)

	Sf1831 = Sf1041*Sf1462 + tmp8*MHiggs2(4)**2 - 
     &   4*Sf1040*Sf237*MStau2(2)

	Sf854 = MHiggs2(4)*(-3*Sf853 + 6*MHiggs2(4)) + 
     &   6*MStau2(1)**2 + MW2*(21*MW2 + 7*MStau2(1))

	Sf3608 = MHiggs2(4)*(-3*Sf3607 + 6*MHiggs2(4)) + 
     &   6*MSbot2(2)**2 + MW2*(21*MW2 + 7*MSbot2(2))

	Sf3281 = MHiggs2(4)*(-3*Sf3280 + 6*MHiggs2(4)) + 
     &   6*MSbot2(1)**2 + MW2*(21*MW2 + 7*MSbot2(1))

	Sf1670 = MHiggs2(4)*(-3*Sf1669 + 6*MHiggs2(4)) + 
     &   6*MStau2(2)**2 + MW2*(21*MW2 + 7*MStau2(2))

	tmp7 = -43*MZ2**2 + MHiggs2(2)*(-122*MZ2 + 4*MHiggs2(2))

	Sf9775 = 170*MZ2**3 + tmp7*MHiggs2(2) + Sf9581*MStau2(2)

	Sf9582 = 170*MZ2**3 + tmp7*MHiggs2(2) + Sf9581*MStau2(1)

	Sf7490 = 1 + UStop22**2*(7 + 8*UStop22**2)

	Sf7478 = 4 - 3*UStop22**2

	Sf7415 = 32 - 30*UStop22**2

	Sf7412 = 8 - 6*UStop22**2

	Sf8756 = -(Sf8755*MStop2(1)) + 
     &   72*(MSbot2(2) - MSbot2(3)*USbot22**2)*UStop11 - 
     &   52*MStop2(2)*UStop12*UStop21*UStop22

	Sf7448 = 2 - UStau22**2

	Sf7334 = 8 - 6*UStau22**2

	Sf7331 = 4 - 2*UStau22**2

	Sf8825 = 124 - 144*USbot22**2

	Sf8840 = -8 + 6*USbot22**2

	Sf8864 = -208 - 618*USbot22**2

	Sf8855 = -464 + 360*USbot22**2

	Sf7363 = 8 - 6*USbot22**2

	Sf9646 = Sf9639*MStau2(1)*MStau2(2) - 
     &   Sf9645*MStau2(2)**3 + 
     &   MZ2*(4*MZ2*Sf9290**3*Sf9291**3*MStau2(1)**2 + 
     &      Sf9640*MStau2(2)**2)

	Sf9697 = MStau2(1)*
     &    (MZ2*(MZ2 - 2*MStau2(2)) - 10*MStau2(2)**2) + 
     &   MStau2(2)*(MZ2*(MZ2 - 3*MStau2(2)) + 4*MStau2(2)**2) + 
     &   MStau2(1)**2*(-3*MZ2 + 6*MStau2(2) + 2*MStau2(3))

	Sf1650 = Sf829*MSneu2(1)**3 + Sf832*MStau2(2)**3 + 
     &   MSneu2(1)*(Sf830*MSneu2(1)*MStau2(2) + Sf831*MStau2(2)**2)

	Sf833 = Sf829*MSneu2(1)**3 + Sf832*MStau2(1)**3 + 
     &   MSneu2(1)*(Sf830*MSneu2(1)*MStau2(1) + Sf831*MStau2(1)**2)

	Sf1620 = Sf790*MSneu2(1)**3 + Sf793*MStau2(2)**3 + 
     &   MSneu2(1)*(3*Sf791*MSneu2(1)*MStau2(2) + 
     &      6*Sf792*MStau2(2)**2)

	Sf794 = Sf790*MSneu2(1)**3 + Sf793*MStau2(1)**3 + 
     &   MSneu2(1)*(3*Sf791*MSneu2(1)*MStau2(1) + 
     &      6*Sf792*MStau2(1)**2)

	Sf1351 = Sf425*MSneu2(1)**3 + 4*Sf427*MStau2(2)**3 + 
     &   MSneu2(1)*(294*MSneu2(1)*MStau2(2) + 3*Sf426*MStau2(2)**2)

	Sf428 = Sf425*MSneu2(1)**3 + 4*Sf427*MStau2(1)**3 + 
     &   MSneu2(1)*(294*MSneu2(1)*MStau2(1) + 3*Sf426*MStau2(1)**2)

	Sf1841 = Sf1060*MSneu2(1)**3 - Sf199*MStau2(2)**3 + 
     &   MSneu2(1)*(Sf650*MSneu2(1)*MStau2(2) + 
     &      3*Sf1061*MStau2(2)**2)

	Sf1062 = Sf1060*MSneu2(1)**3 - Sf199*MStau2(1)**3 + 
     &   MSneu2(1)*(Sf650*MSneu2(1)*MStau2(1) + 
     &      3*Sf1061*MStau2(1)**2)

	Sf6152 = 2*Sf3766*MSbot2(2)**3 + Sf3768*MStop2(2)**3 + 
     &   MSbot2(2)*(-60*C2B*MSbot2(2)*MStop2(2) - 
     &      2*Sf3767*MStop2(2)**2)

	Sf3904 = 2*Sf3766*MSbot2(2)**3 + Sf3768*MStop2(1)**3 + 
     &   MSbot2(2)*(-60*C2B*MSbot2(2)*MStop2(1) - 
     &      2*Sf3767*MStop2(1)**2)

	Sf4888 = 2*MW2*MSbot2(2) + 
     &   MSbot2(2)*(7*MSbot2(2) - 4*MStop2(2)) - 3*MStop2(2)**2

	Sf2299 = 2*MW2*MSbot2(2) + 
     &   MSbot2(2)*(7*MSbot2(2) - 4*MStop2(1)) - 3*MStop2(1)**2

	Sf3769 = 2*Sf3766*MSbot2(1)**3 + Sf3768*MStop2(1)**3 + 
     &   MSbot2(1)*(-60*C2B*MSbot2(1)*MStop2(1) - 
     &      2*Sf3767*MStop2(1)**2)

	Sf4854 = 2*MW2*MSbot2(1) + 
     &   MSbot2(1)*(7*MSbot2(1) - 4*MStop2(2)) - 3*MStop2(2)**2

	Sf2248 = 2*MW2*MSbot2(1) + 
     &   MSbot2(1)*(7*MSbot2(1) - 4*MStop2(1)) - 3*MStop2(1)**2

	Sf10235 = 2*Ss984 + 1/2.D0*Ss982 + (-1 - 2*Ss969)*MStop2(1)

	Sf10255 = 2*Ss1000 + 1/2.D0*Ss998 + (-1 - 2*Ss988)*MStop2(2)

	Sf9961 = Sf9957 - Sf9958*MStau2(2) + 
     &   MStau2(2)**2*(2*Sf9959 + Sf9960*MStau2(2))

	Sf9888 = Sf9884 - Sf9885*MStau2(2) + 
     &   MStau2(2)**2*(Sf9886 - Sf9887*MStau2(2))

	Sf1434 = 3*ML2*Sf47 + 
     &   MSneu2(1)*(11*MSneu2(1) - 4*MStau2(2)) - 7*MStau2(2)**2

	Sf1545 = 3*ML2*Sf47 + 
     &   MSneu2(1)*(-23*MSneu2(1) - 4*MStau2(2)) + 27*MStau2(2)**2

	Sf1503 = ML2*Sf47 + 
     &   MSneu2(1)*(-11*MSneu2(1) - 4*MStau2(2)) + 15*MStau2(2)**2

	Sf6287 = Sf3976 - 4*Sf3977*MStop2(2) + 
     &   (8*Sf3978 - 12*MStop2(2))*MStop2(2)**2

	Sf3979 = Sf3976 - 4*Sf3977*MStop2(1) + 
     &   (8*Sf3978 - 12*MStop2(1))*MStop2(1)**2

	Sf6273 = Sf3879 - 4*Sf3880*MStop2(2) + 
     &   (8*Sf3881 - 12*MStop2(2))*MStop2(2)**2

	Sf3882 = Sf3879 - 4*Sf3880*MStop2(1) + 
     &   (8*Sf3881 - 12*MStop2(1))*MStop2(1)**2

	Sf525 = 3*ML2*Sf34 + 
     &   MSneu2(1)*(11*MSneu2(1) - 4*MStau2(1)) - 7*MStau2(1)**2

	Sf671 = 3*ML2*Sf34 + 
     &   MSneu2(1)*(-23*MSneu2(1) - 4*MStau2(1)) + 27*MStau2(1)**2

	Sf606 = ML2*Sf34 + 
     &   MSneu2(1)*(-11*MSneu2(1) - 4*MStau2(1)) + 15*MStau2(1)**2

	Sf10363 = 2*MZ2**2*MHiggs2(1) + 3*Sf9113*MStop2(2)**2 + 
     &   MZ2*(-3*MZ2**2 + 2*Sf9131*Sf9132*MStop2(2))

	Sf10343 = 2*MZ2**2*MHiggs2(1) + 3*Sf9113*MStop2(1)**2 + 
     &   MZ2*(-3*MZ2**2 + 2*Sf9131*Sf9132*MStop2(1))

	Sf10515 = -MZ2**3 + 2*Sf10514*MStop2(2) + 
     &   MHiggs2(2)*(MZ2*MHiggs2(2) - 2*MStop2(2)**2)

	Sf10453 = -MZ2**3 + 2*Sf10452*MStop2(1) + 
     &   MHiggs2(2)*(MZ2*MHiggs2(2) - 2*MStop2(1)**2)

	Sf9013 = -MZ2**3 + 2*Sf9012*MSbot2(2) + 
     &   MHiggs2(2)*(MZ2*MHiggs2(2) - 2*MSbot2(2)**2)

	Sf8936 = -MZ2**3 + 2*Sf8935*MSbot2(1) + 
     &   MHiggs2(2)*(MZ2*MHiggs2(2) - 2*MSbot2(1)**2)

	Sf10410 = -MZ2**3 + 2*Sf10409*MStop2(2) + 
     &   MHiggs2(1)*(MZ2*MHiggs2(1) - 2*MStop2(2)**2)

	Sf10358 = -MZ2**3 + 2*Sf10357*MStop2(1) + 
     &   MHiggs2(1)*(MZ2*MHiggs2(1) - 2*MStop2(1)**2)

	Sf10146 = -MZ2**3 + 2*Sf10145*MStau2(2) + 
     &   MHiggs2(1)*(MZ2*MHiggs2(1) - 2*MStau2(2)**2)

	Sf10077 = -MZ2**3 + 2*Sf10076*MStau2(1) + 
     &   MHiggs2(1)*(MZ2*MHiggs2(1) - 2*MStau2(1)**2)

	Sf9211 = -MZ2**3 + 2*Sf9210*MSbot2(2) + 
     &   MHiggs2(1)*(MZ2*MHiggs2(1) - 2*MSbot2(2)**2)

	Sf9138 = -MZ2**3 + 2*Sf9137*MSbot2(1) + 
     &   MHiggs2(1)*(MZ2*MHiggs2(1) - 2*MSbot2(1)**2)

	Sf1214 = 7*MW2**2*MHiggs2(4) + 6*Sf173*MStau2(2)**2 + 
     &   MW2*(-10*MW2**2 + Sf221*MStau2(2))

	Sf222 = 7*MW2**2*MHiggs2(4) + 6*Sf173*MStau2(1)**2 + 
     &   MW2*(-10*MW2**2 + Sf221*MStau2(1))

	Sf2685 = 7*MW2**2*MHiggs2(4) + 3*Sf173*MSbot2(2)**2 + 
     &   MW2*(-10*MW2**2 + Sf180*MSbot2(2))

	Sf2493 = 7*MW2**2*MHiggs2(4) + 3*Sf173*MSbot2(1)**2 + 
     &   MW2*(-10*MW2**2 + Sf180*MSbot2(1))

	Sf6186 = -MSbot2(2)**2 + MW2*(MW2 - 2*MStop2(2)) + MStop2(2)**2

	Sf6104 = -MSbot2(2)**2 + MW2*(MW2 - 2*MStop2(1)) + MStop2(1)**2

	Sf5808 = -MSbot2(1)**2 + MW2*(MW2 - 2*MStop2(2)) + MStop2(2)**2

	Sf5715 = -MSbot2(1)**2 + MW2*(MW2 - 2*MStop2(1)) + MStop2(1)**2

	Sf1802 = ML2*(ML2 + 2*MW2*Sf2) - MW2**2*Sf983 + 
     &   23*Sf241*MStau2(2)

	Sf984 = ML2*(ML2 + 2*MW2*Sf2) - MW2**2*Sf983 + 
     &   23*Sf241*MStau2(1)

	Sf1531 = 6*MHiggs2(4)**2 + 
     &   MSneu2(1)*(-13*MSneu2(1) - 62*MStau2(2)) - 21*MStau2(2)**2

	Sf637 = 6*MHiggs2(4)**2 + 
     &   MSneu2(1)*(-13*MSneu2(1) - 62*MStau2(1)) - 21*MStau2(1)**2

	Sf6783 = 6*MHiggs2(4)**2 + 
     &   MSbot2(2)*(-13*MSbot2(2) - 62*MStop2(2)) - 21*MStop2(2)**2

	Sf4425 = 6*MHiggs2(4)**2 + 
     &   MSbot2(2)*(-13*MSbot2(2) - 62*MStop2(1)) - 21*MStop2(1)**2

	end

