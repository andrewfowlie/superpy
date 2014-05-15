* gm2_2L_full5.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full5
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp14, tmp15, tmp16, tmp17, tmp18, tmp19, tmp20
	ComplexType tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27
	ComplexType tmp28, tmp29

	Sf8051 = 243*MSq2Diff(1,2)**8 + 
     &   MW2*(-(Sf8045*MSq2Diff(1,2)**6) + 
     &      MW2*(2*Sf8047*MSq2Diff(1,2)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(8*MW2 + Sf8044) + 2*Sf8046) - 
     &                  Sf8048) + Sf8050) - 
     &            2*Sf8049*MSq2Diff(1,2)**2)))

	Sf7920 = 243*MSq2Diff(2,1)**8 + 
     &   MW2*(-(Sf7914*MSq2Diff(2,1)**6) + 
     &      MW2*(2*Sf7916*MSq2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(8*MW2 + Sf7913) + 2*Sf7915) - 
     &                  Sf7917) + Sf7919) - 
     &            2*Sf7918*MSq2Diff(2,1)**2)))

	Sf7842 = 243*MSq2Diff(1,1)**8 + 
     &   MW2*(-(Sf7836*MSq2Diff(1,1)**6) + 
     &      MW2*(2*Sf7838*MSq2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(8*MW2 + Sf7835) + 2*Sf7837) - 
     &                  Sf7839) + Sf7841) - 
     &            2*Sf7840*MSq2Diff(1,1)**2)))

	Sf7606 = -97*MSl2Diff(2,1)**8 + 
     &   MW2*(4*Sf7600*MSl2Diff(2,1)**6 + 
     &      MW2*(-(Sf7602*MSl2Diff(2,1)**4) + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(2*MW2 - Sf7599) + Sf7601) + 
     &                  2*Sf7603) - Sf7605) + 
     &            4*Sf7604*MSl2Diff(2,1)**2)))

	Sf7548 = -97*MSl2Diff(1,1)**8 + 
     &   MW2*(4*Sf7542*MSl2Diff(1,1)**6 + 
     &      MW2*(-(Sf7544*MSl2Diff(1,1)**4) + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(2*MW2 - Sf7541) + Sf7543) + 
     &                  2*Sf7545) - Sf7547) + 
     &            4*Sf7546*MSl2Diff(1,1)**2)))

	Sf8402 = 222*MSq2Diff(2,2)**8 + 
     &   MW2*(-(Sf8396*MSq2Diff(2,2)**6) + 
     &      MW2*(2*Sf8398*MSq2Diff(2,2)**4 + 
     &         MW2*(MW2*(Sf8401 + 
     &               MW2*(-2*Sf8399 + MW2*(2*Sf8397 - 24*MW2*Sf99))
     &               ) - Sf8400*MSq2Diff(2,2)**2)))

	Sf8306 = 222*MSq2Diff(1,2)**8 + 
     &   MW2*(-(Sf8300*MSq2Diff(1,2)**6) + 
     &      MW2*(2*Sf8302*MSq2Diff(1,2)**4 + 
     &         MW2*(MW2*(Sf8305 + 
     &               MW2*(-2*Sf8303 + MW2*(2*Sf8301 - 24*MW2*Sf91))
     &               ) - Sf8304*MSq2Diff(1,2)**2)))

	Sf8137 = -1121*MSq2Diff(2,2)**8 + 
     &   MW2*(81*Sf8130*MSq2Diff(2,2)**6 + 
     &      MW2*(-2*Sf8133*MSq2Diff(2,2)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*Sf8131 - 4*Sf8132) + Sf8134) - 
     &               Sf8136) + Sf8135*MSq2Diff(2,2)**2)))

	Sf8121 = -674*MSq2Diff(2,2)**8 + 
     &   MW2*(Sf8115*MSq2Diff(2,2)**6 + 
     &      MW2*(-2*Sf8117*MSq2Diff(2,2)**4 + 
     &         MW2*(MW2*(2*Sf8120 + 
     &               MW2*(-2*Sf8119 + 
     &                  MW2*(2*Sf8116 - 168*MW2*Sf99))) - 
     &            9*Sf8118*MSq2Diff(2,2)**2)))

	Sf8059 = -1121*MSq2Diff(1,2)**8 + 
     &   MW2*(81*Sf8052*MSq2Diff(1,2)**6 + 
     &      MW2*(-2*Sf8055*MSq2Diff(1,2)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*Sf8053 - 4*Sf8054) + Sf8056) - 
     &               Sf8058) + Sf8057*MSq2Diff(1,2)**2)))

	Sf8043 = -674*MSq2Diff(1,2)**8 + 
     &   MW2*(Sf8037*MSq2Diff(1,2)**6 + 
     &      MW2*(-2*Sf8039*MSq2Diff(1,2)**4 + 
     &         MW2*(MW2*(2*Sf8042 + 
     &               MW2*(-2*Sf8041 + 
     &                  MW2*(2*Sf8038 - 168*MW2*Sf91))) - 
     &            9*Sf8040*MSq2Diff(1,2)**2)))

	Sf7928 = -1121*MSq2Diff(2,1)**8 + 
     &   MW2*(81*Sf7921*MSq2Diff(2,1)**6 + 
     &      MW2*(-2*Sf7924*MSq2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*Sf7922 - 4*Sf7923) + Sf7925) - 
     &               Sf7927) + Sf7926*MSq2Diff(2,1)**2)))

	Sf7850 = -1121*MSq2Diff(1,1)**8 + 
     &   MW2*(81*Sf7843*MSq2Diff(1,1)**6 + 
     &      MW2*(-2*Sf7846*MSq2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*Sf7844 - 4*Sf7845) + Sf7847) - 
     &               Sf7849) + Sf7848*MSq2Diff(1,1)**2)))

	Sf7912 = 674*MSq2Diff(2,1)**8 + 
     &   MW2*(-(Sf7906*MSq2Diff(2,1)**6) + 
     &      MW2*(2*Sf7908*MSq2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(168*MW2*Sf77 - 2*Sf7907) + 2*Sf7910) - 
     &               2*Sf7911) + 9*Sf7909*MSq2Diff(2,1)**2)))

	Sf8378 = -222*MSq2Diff(2,1)**8 + 
     &   MW2*(Sf8372*MSq2Diff(2,1)**6 + 
     &      MW2*(-2*Sf8374*MSq2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(24*MW2*Sf77 - 2*Sf8373) + 2*Sf8375) - 
     &               Sf8377) + Sf8376*MSq2Diff(2,1)**2)))

	Sf7628 = -427*MSl2Diff(2,1)**8 + 
     &   MW2*(2*Sf7622*MSl2Diff(2,1)**6 + 
     &      MW2*(-3*Sf7624*MSl2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*Sf7621 - Sf7623) + 6*Sf7625) - 
     &               Sf7627) + Sf7626*MSl2Diff(2,1)**2)))

	Sf7570 = -427*MSl2Diff(1,1)**8 + 
     &   MW2*(2*Sf7564*MSl2Diff(1,1)**6 + 
     &      MW2*(-3*Sf7566*MSl2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*Sf7563 - Sf7565) + 6*Sf7567) - 
     &               Sf7569) + Sf7568*MSl2Diff(1,1)**2)))

	Sf7834 = 674*MSq2Diff(1,1)**8 + 
     &   MW2*(-(Sf7828*MSq2Diff(1,1)**6) + 
     &      MW2*(2*Sf7830*MSq2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(168*MW2*Sf69 - 2*Sf7829) + 2*Sf7832) - 
     &               2*Sf7833) + 9*Sf7831*MSq2Diff(1,1)**2)))

	Sf8282 = -222*MSq2Diff(1,1)**8 + 
     &   MW2*(Sf8276*MSq2Diff(1,1)**6 + 
     &      MW2*(-2*Sf8278*MSq2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(24*MW2*Sf69 - 2*Sf8277) + 2*Sf8279) - 
     &               Sf8281) + Sf8280*MSq2Diff(1,1)**2)))

	Sf7613 = 310*MSl2Diff(2,1)**8 + 
     &   MW2*(-3*Sf7607*MSl2Diff(2,1)**6 + 
     &      MW2*(3*Sf7609*MSl2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(84*MW2*Sf47 - Sf7608) + 6*Sf7611) - 
     &               6*Sf7612) + Sf7610*MSl2Diff(2,1)**2)))

	Sf7691 = -106*MSl2Diff(2,1)**8 + 
     &   MW2*(9*Sf7685*MSl2Diff(2,1)**6 + 
     &      MW2*(-3*Sf7687*MSl2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(12*MW2*Sf47 - Sf7686) + 6*Sf7688) - 
     &               6*Sf7690) + Sf7689*MSl2Diff(2,1)**2)))

	Sf7555 = 310*MSl2Diff(1,1)**8 + 
     &   MW2*(-3*Sf7549*MSl2Diff(1,1)**6 + 
     &      MW2*(3*Sf7551*MSl2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(84*MW2*Sf34 - Sf7550) + 6*Sf7553) - 
     &               6*Sf7554) + Sf7552*MSl2Diff(1,1)**2)))

	Sf7677 = -106*MSl2Diff(1,1)**8 + 
     &   MW2*(9*Sf7671*MSl2Diff(1,1)**6 + 
     &      MW2*(-3*Sf7673*MSl2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(12*MW2*Sf34 - Sf7672) + 6*Sf7674) - 
     &               6*Sf7676) + Sf7675*MSl2Diff(1,1)**2)))

	Sf8460 = 243*MSq2Diff(2,2)**7 + 
     &   MW2*(-(Sf8455*MSq2Diff(2,2)**5) + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(55*MW2 - 2*Sf8454) - Sf8456) + Sf8458)-
     &              Sf8459) + 2*Sf8457*MSq2Diff(2,2)**3))

	Sf8449 = 243*MSq2Diff(1,2)**7 + 
     &   MW2*(-(Sf8444*MSq2Diff(1,2)**5) + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(55*MW2 - 2*Sf8443) - Sf8445) + Sf8447)-
     &              Sf8448) + 2*Sf8446*MSq2Diff(1,2)**3))

	Sf8436 = 243*MSq2Diff(2,1)**7 + 
     &   MW2*(-(Sf8431*MSq2Diff(2,1)**5) + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(55*MW2 - 2*Sf8430) - Sf8432) + Sf8434)-
     &              Sf8435) + 2*Sf8433*MSq2Diff(2,1)**3))

	Sf8425 = 243*MSq2Diff(1,1)**7 + 
     &   MW2*(-(Sf8420*MSq2Diff(1,1)**5) + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(55*MW2 - 2*Sf8419) - Sf8421) + Sf8423)-
     &              Sf8424) + 2*Sf8422*MSq2Diff(1,1)**3))

	Sf7745 = -97*MSl2Diff(2,1)**7 + 
     &   MW2*(Sf7740*MSl2Diff(2,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(35*MW2 - Sf7739) - 6*Sf7741) + Sf7743)-
     &              Sf7744) - 3*Sf7742*MSl2Diff(2,1)**3))

	Sf7721 = -97*MSl2Diff(1,1)**7 + 
     &   MW2*(Sf7716*MSl2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(35*MW2 - Sf7715) - 6*Sf7717) + Sf7719)-
     &              Sf7720) - 3*Sf7718*MSl2Diff(1,1)**3))

	Sf5970 = 2*MSq2Diff(2,2)**7 + 
     &   MW2*(-2*Sf5965*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(5*MW2 + 6*Sf5964) - 3*Sf5967) + 
     &               4*Sf5968) - Sf5969) + 
     &         18*Sf5966*MSq2Diff(2,2)**3))

	Sf5544 = 2*MSq2Diff(1,2)**7 + 
     &   MW2*(-2*Sf5539*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(5*MW2 + 6*Sf5538) - 3*Sf5540) + 
     &               4*Sf5542) - Sf5543) + 
     &         18*Sf5541*MSq2Diff(1,2)**3))

	Sf3545 = 2*MSq2Diff(2,1)**7 + 
     &   MW2*(-2*Sf3540*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(5*MW2 + 6*Sf3539) - 3*Sf3542) + 
     &               4*Sf3543) - Sf3544) + 
     &         18*Sf3541*MSq2Diff(2,1)**3))

	Sf3183 = 2*MSq2Diff(1,1)**7 + 
     &   MW2*(-2*Sf3178*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*
     &                (MW2*(5*MW2 + 6*Sf3177) - 3*Sf3180) + 
     &               4*Sf3181) - Sf3182) + 
     &         18*Sf3179*MSq2Diff(1,1)**3))

	Sf8395 = -243*MSq2Diff(2,2)**7 + 
     &   MW2*(25*Sf8389*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 - Sf8390) + Sf8391) - 
     &               Sf8393) + Sf8394) - 2*Sf8392*MSq2Diff(2,2)**3)
     &      )

	Sf8371 = -243*MSq2Diff(2,1)**7 + 
     &   MW2*(25*Sf8365*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 - Sf8366) + Sf8367) - 
     &               Sf8369) + Sf8370) - 2*Sf8368*MSq2Diff(2,1)**3)
     &      )

	Sf8299 = -243*MSq2Diff(1,2)**7 + 
     &   MW2*(25*Sf8293*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 - Sf8294) + Sf8295) - 
     &               Sf8297) + Sf8298) - 2*Sf8296*MSq2Diff(1,2)**3)
     &      )

	Sf8275 = -243*MSq2Diff(1,1)**7 + 
     &   MW2*(25*Sf8269*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 - Sf8270) + Sf8271) - 
     &               Sf8273) + Sf8274) - 2*Sf8272*MSq2Diff(1,1)**3)
     &      )

	Sf7698 = -97*MSl2Diff(2,1)**7 + 
     &   MW2*(Sf7693*MSl2Diff(2,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 + Sf7692) - 6*Sf7694) + 
     &               Sf7696) - Sf7697) - 3*Sf7695*MSl2Diff(2,1)**3)
     &      )

	Sf7684 = -97*MSl2Diff(1,1)**7 + 
     &   MW2*(Sf7679*MSl2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 + Sf7678) - 6*Sf7680) + 
     &               Sf7682) - Sf7683) - 3*Sf7681*MSl2Diff(1,1)**3)
     &      )

	Sf5914 = 2*MSq2Diff(2,2)**7 + 
     &   MW2*(-2*Sf5909*MSq2Diff(2,2)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 + 2*Sf5149) - Sf5910) + 
     &               4*Sf5912) - Sf5913) + 
     &         2*Sf5911*MSq2Diff(2,2)**3))

	Sf5480 = 2*MSq2Diff(1,2)**7 + 
     &   MW2*(-2*Sf5475*MSq2Diff(1,2)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 + 2*Sf5030) - Sf5476) + 
     &               4*Sf5478) - Sf5479) + 
     &         2*Sf5477*MSq2Diff(1,2)**3))

	Sf3470 = 2*MSq2Diff(2,1)**7 + 
     &   MW2*(-2*Sf3465*MSq2Diff(2,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 + 2*Sf2742) - Sf3466) + 
     &               4*Sf3468) - Sf3469) + 
     &         2*Sf3467*MSq2Diff(2,1)**3))

	Sf3083 = 2*MSq2Diff(1,1)**7 + 
     &   MW2*(-2*Sf3078*MSq2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*(MW2 + 2*Sf2565) - Sf3079) + 
     &               4*Sf3081) - Sf3082) + 
     &         2*Sf3080*MSq2Diff(1,1)**3))

	Sf8332 = -25*MSl2Diff(2,1)**7 + 
     &   MW2*(MSl2Diff(2,1)**5*(103*MSneu2(1) + 140*MStau2(2)) + 
     &      MW2*(-(MSl2Diff(2,1)**3*
     &            (325*MStau2(2)**2 + 
     &              MSneu2(1)*(162*MSneu2(1) + 225*MStau2(2)))) + 
     &         MW2*(400*MStau2(2)**4 + 
     &            MW2*(MSneu2(1)*
     &                (MW2*(12*MW2 - 36*MSneu2(1)) + 
     &                  67*MSneu2(1)**2) - 
     &               3*(-11*MSneu2(1)**2 + 
     &                  MW2*(5*MW2 + 14*MSneu2(1)))*MStau2(2) + 
     &               5*(20*MW2 + 27*MSneu2(1))*MStau2(2)**2 - 
     &               275*MStau2(2)**3) + 
     &            MSneu2(1)*
     &             (-433*MStau2(2)**3 + 
     &               MSneu2(1)*
     &                (75*MStau2(2)**2 + 
     &                  MSneu2(1)*(-127*MSneu2(1) + 85*MStau2(2))))
     &            )))

	Sf8325 = -25*MSl2Diff(1,1)**7 + 
     &   MW2*(Sf8320*MSl2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(MW2*(3*MW2*Sf8319 - 2*Sf8321) + 
     &               Sf8323) - Sf8324) - Sf8322*MSl2Diff(1,1)**3))

	Sf8251 = 59*MSq2Diff(1,2)**7 + 
     &   MW2*(-(Sf8246*MSq2Diff(1,2)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf8245 - Sf8247) + Sf8249) - 
     &            Sf8250) + 2*Sf8248*MSq2Diff(1,2)**3))

	Sf8239 = 59*MSq2Diff(1,1)**7 + 
     &   MW2*(-(Sf8234*MSq2Diff(1,1)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf8233 - Sf8235) + Sf8237) - 
     &            Sf8238) + 2*Sf8236*MSq2Diff(1,1)**3))

	Sf8222 = 59*MSq2Diff(2,2)**7 + 
     &   MW2*(-(Sf8217*MSq2Diff(2,2)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf8216 - Sf8218) + Sf8220) - 
     &            Sf8221) + 2*Sf8219*MSq2Diff(2,2)**3))

	Sf8210 = 59*MSq2Diff(2,1)**7 + 
     &   MW2*(-(Sf8205*MSq2Diff(2,1)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf8204 - Sf8206) + Sf8208) - 
     &            Sf8209) + 2*Sf8207*MSq2Diff(2,1)**3))

	Sf8193 = 59*MSq2Diff(2,2)**7 + 
     &   MW2*(-(Sf8188*MSq2Diff(2,2)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf8187 - Sf8189) + Sf8191) - 
     &            Sf8192) + 2*Sf8190*MSq2Diff(2,2)**3))

	Sf8173 = 59*MSq2Diff(1,2)**7 + 
     &   MW2*(-(Sf8168*MSq2Diff(1,2)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf8167 - Sf8169) + Sf8171) - 
     &            Sf8172) + 2*Sf8170*MSq2Diff(1,2)**3))

	Sf8153 = 59*MSq2Diff(2,2)**8 + 
     &   MW2*(-3*Sf8147*MSq2Diff(2,2)**6 + 
     &      MW2*(10*Sf8149*MSq2Diff(2,2)**4 + 
     &         MW2*(MW2*(MW2*(MW2*Sf8148 - Sf8150) + Sf8152) - 
     &            10*Sf8151*MSq2Diff(2,2)**2)))

	Sf8075 = 59*MSq2Diff(1,2)**8 + 
     &   MW2*(-3*Sf8069*MSq2Diff(1,2)**6 + 
     &      MW2*(10*Sf8071*MSq2Diff(1,2)**4 + 
     &         MW2*(MW2*(MW2*(MW2*Sf8070 - Sf8072) + Sf8074) - 
     &            10*Sf8073*MSq2Diff(1,2)**2)))

	Sf7984 = 59*MSq2Diff(2,1)**7 + 
     &   MW2*(-(Sf7979*MSq2Diff(2,1)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf7978 - Sf7980) + Sf7982) - 
     &            Sf7983) + 2*Sf7981*MSq2Diff(2,1)**3))

	Sf7964 = 59*MSq2Diff(1,1)**7 + 
     &   MW2*(-(Sf7959*MSq2Diff(1,1)**5) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2*Sf7958 - Sf7960) + Sf7962) - 
     &            Sf7963) + 2*Sf7961*MSq2Diff(1,1)**3))

	Sf7944 = 59*MSq2Diff(2,1)**8 + 
     &   MW2*(-3*Sf7938*MSq2Diff(2,1)**6 + 
     &      MW2*(10*Sf7940*MSq2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*(MW2*Sf7939 - Sf7941) + Sf7943) - 
     &            10*Sf7942*MSq2Diff(2,1)**2)))

	Sf7866 = 59*MSq2Diff(1,1)**8 + 
     &   MW2*(-3*Sf7860*MSq2Diff(1,1)**6 + 
     &      MW2*(10*Sf7862*MSq2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*(MW2*Sf7861 - Sf7863) + Sf7865) - 
     &            10*Sf7864*MSq2Diff(1,1)**2)))

	Sf7667 = 25*MSl2Diff(2,1)**7 + 
     &   MW2*(-(Sf7662*MSl2Diff(2,1)**5) + 
     &      MW2*(Sf7664*MSl2Diff(2,1)**3 + 
     &         MW2*(Sf7666 + 
     &            MW2*(-Sf7665 + MW2*(Sf7663 - 3*MW2*MSneu2(1))))))

	Sf7655 = 25*MSl2Diff(1,1)**7 + 
     &   MW2*(-(Sf7650*MSl2Diff(1,1)**5) + 
     &      MW2*(Sf7652*MSl2Diff(1,1)**3 + 
     &         MW2*(Sf7654 + 
     &            MW2*(-Sf7653 + MW2*(Sf7651 - 3*MW2*MSneu2(1))))))

	Sf7636 = 25*MSl2Diff(2,1)**8 + 
     &   MW2*(-4*Sf7630*MSl2Diff(2,1)**6 + 
     &      MW2*(5*Sf7632*MSl2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*(3*MW2*Sf7631 - Sf7633) + Sf7635) - 
     &            20*Sf7634*MSl2Diff(2,1)**2)))

	Sf7620 = 25*Sf47*MSl2Diff(2,1)**6 + 
     &   MW2*(-2*Sf7616*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(3*MW2*Sf7614 - 10*Sf7615) + 
     &               Sf7617) - 2*Sf7619) + Sf7618*MSl2Diff(2,1)**2)
     &      )

	Sf7578 = 25*MSl2Diff(1,1)**8 + 
     &   MW2*(-4*Sf7572*MSl2Diff(1,1)**6 + 
     &      MW2*(5*Sf7574*MSl2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*(3*MW2*Sf7573 - Sf7575) + Sf7577) - 
     &            20*Sf7576*MSl2Diff(1,1)**2)))

	Sf7562 = 25*Sf34*MSl2Diff(1,1)**6 + 
     &   MW2*(-2*Sf7558*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(3*MW2*Sf7556 - 10*Sf7557) + 
     &               Sf7559) - 2*Sf7561) + Sf7560*MSl2Diff(1,1)**2)
     &      )

	Sf576 = 15*MSl2Diff(1,1)**7 + 
     &   MW2*(-3*Sf565*MSl2Diff(1,1)**5 + 
     &      MW2*(30*Sf568*MSl2Diff(1,1)**3 + 
     &         MW2*(6*Sf575 + 
     &            MW2*(-3*Sf574 + MW2*(Sf573 + 6*MW2*MStau2(1))))))

	Sf486 = -3*MSl2Diff(1,1)**7 + 
     &   MW2*(3*Sf471*MSl2Diff(1,1)**5 + 
     &      MW2*(MW2*(MW2*(3*Sf475 + MW2*(2*MW2*Sf470 - Sf484)) - 
     &            2*Sf485) - 6*Sf472*MSl2Diff(1,1)**3))

	Sf1477 = 15*MSl2Diff(2,1)**7 + 
     &   MW2*(-3*Sf1468*MSl2Diff(2,1)**5 + 
     &      MW2*(30*Sf1470*MSl2Diff(2,1)**3 + 
     &         MW2*(6*Sf1476 + 
     &            MW2*(-3*Sf1475 + MW2*(Sf1474 + 6*MW2*MStau2(2))))
     &         ))

	Sf1403 = -3*MSl2Diff(2,1)**7 + 
     &   MW2*(3*Sf1388*MSl2Diff(2,1)**5 + 
     &      MW2*(MW2*(MW2*(3*Sf1392 + 
     &               MW2*(2*MW2*Sf1387 - Sf1401)) - 2*Sf1402) - 
     &         6*Sf1389*MSl2Diff(2,1)**3))

	Sf8106 = 160*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf8102*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(54*MW2 - 18*Sf8101) + Sf8104) - 
     &            Sf8105) + 3*Sf8103*MSq2Diff(2,2)**2))

	Sf8028 = 160*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf8024*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(54*MW2 - 18*Sf8023) + Sf8026) - 
     &            Sf8027) + 3*Sf8025*MSq2Diff(1,2)**2))

	Sf7897 = 160*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf7893*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(54*MW2 - 18*Sf7892) + Sf7895) - 
     &            Sf7896) + 3*Sf7894*MSq2Diff(2,1)**2))

	Sf7819 = 160*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf7815*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(54*MW2 - 18*Sf7814) + Sf7817) - 
     &            Sf7818) + 3*Sf7816*MSq2Diff(1,1)**2))

	Sf8186 = 118*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf8182*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(32*MW2 - 2*Sf8181) + Sf8184) - 
     &            3*Sf8185) + 2*Sf8183*MSq2Diff(2,2)**2))

	Sf8166 = 118*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf8162*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(32*MW2 - 2*Sf8161) + Sf8164) - 
     &            3*Sf8165) + 2*Sf8163*MSq2Diff(1,2)**2))

	Sf7977 = 118*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf7973*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(32*MW2 - 2*Sf7972) + Sf7975) - 
     &            3*Sf7976) + 2*Sf7974*MSq2Diff(2,1)**2))

	Sf7957 = 118*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf7953*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(32*MW2 - 2*Sf7952) + Sf7955) - 
     &            3*Sf7956) + 2*Sf7954*MSq2Diff(1,1)**2))

	Sf363 = 3*MSl2Diff(1,1)**6 + 
     &   MW2*(-(Sf358*MSl2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(31*MW2 - 2*Sf359) + Sf361) - 
     &            Sf362) + Sf360*MSl2Diff(1,1)**2))

	Sf1314 = 3*MSl2Diff(2,1)**6 + 
     &   MW2*(-(Sf1309*MSl2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(31*MW2 - 2*Sf1310) + Sf1312) - 
     &            Sf1313) + Sf1311*MSl2Diff(2,1)**2))

	Sf5351 = 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-3*Sf5346*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(23*MW2 - Sf5347) + Sf5349) - 
     &            Sf5350) + Sf5348*MSq2Diff(2,2)**2))

	Sf5243 = 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-3*Sf5238*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(23*MW2 - Sf5239) + Sf5241) - 
     &            Sf5242) + Sf5240*MSq2Diff(1,2)**2))

	Sf2940 = 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-3*Sf2935*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(23*MW2 - Sf2936) + Sf2938) - 
     &            Sf2939) + Sf2937*MSq2Diff(2,1)**2))

	Sf2774 = 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-3*Sf2769*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(23*MW2 - Sf2770) + Sf2772) - 
     &            Sf2773) + Sf2771*MSq2Diff(1,1)**2))

	Sf8698 = 50*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf8694*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(15*MW2 - 4*Sf8693) + Sf8695) - 
     &            Sf8697) + Sf8696*MSq2Diff(2,2)**2))

	Sf8690 = 50*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf8686*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(15*MW2 - 4*Sf8685) + Sf8687) - 
     &            Sf8689) + Sf8688*MSq2Diff(1,2)**2))

	Sf8658 = 50*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf8654*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(15*MW2 - 4*Sf8653) + Sf8655) - 
     &            Sf8657) + Sf8656*MSq2Diff(2,1)**2))

	Sf8650 = 50*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf8646*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(15*MW2 - 4*Sf8645) + Sf8647) - 
     &            Sf8649) + Sf8648*MSq2Diff(1,1)**2))

	Sf7661 = 50*MSl2Diff(2,1)**6 + 
     &   MW2*(-3*Sf7657*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(15*MW2 - 4*Sf7656) + 6*Sf7658) - 
     &            Sf7660) + Sf7659*MSl2Diff(2,1)**2))

	Sf7649 = 50*MSl2Diff(1,1)**6 + 
     &   MW2*(-3*Sf7645*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(15*MW2 - 4*Sf7644) + 6*Sf7646) - 
     &            Sf7648) + Sf7647*MSl2Diff(1,1)**2))

	Sf8100 = -48*MSq2Diff(2,2)**6 + 
     &   MW2*(Sf8096*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(14*MW2 + 6*Sf8095) - Sf8097) + 
     &            Sf8099) - 9*Sf8098*MSq2Diff(2,2)**2))

	Sf8016 = -48*MSq2Diff(1,2)**6 + 
     &   MW2*(Sf8012*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(14*MW2 + 6*Sf8011) - Sf8013) + 
     &            Sf8015) - 9*Sf8014*MSq2Diff(1,2)**2))

	Sf7891 = -48*MSq2Diff(2,1)**6 + 
     &   MW2*(Sf7887*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(14*MW2 + 6*Sf7886) - Sf7888) + 
     &            Sf7890) - 9*Sf7889*MSq2Diff(2,1)**2))

	Sf7813 = -48*MSq2Diff(1,1)**6 + 
     &   MW2*(Sf7809*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(14*MW2 + 6*Sf7808) - Sf7810) + 
     &            Sf7812) - 9*Sf7811*MSq2Diff(1,1)**2))

	Sf7086 = 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf7081*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(11*MW2 - Sf7082) + Sf7084) - 
     &            Sf7085) + Sf7083*MSq2Diff(2,2)**2))

	Sf7046 = 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf7041*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(11*MW2 - Sf7042) + Sf7044) - 
     &            Sf7045) + Sf7043*MSq2Diff(1,2)**2))

	Sf4655 = 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf4650*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(11*MW2 - Sf4651) + Sf4653) - 
     &            Sf4654) + Sf4652*MSq2Diff(2,1)**2))

	Sf4539 = 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf4534*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(11*MW2 - Sf4535) + Sf4537) - 
     &            Sf4538) + Sf4536*MSq2Diff(1,1)**2))

	Sf5386 = 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-5*Sf5381*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(9*MW2 - Sf5382) + Sf5384) - 
     &            Sf5385) + Sf5383*MSq2Diff(2,2)**2))

	Sf5278 = 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-5*Sf5273*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(9*MW2 - Sf5274) + Sf5276) - 
     &            Sf5277) + Sf5275*MSq2Diff(1,2)**2))

	Sf2978 = 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-5*Sf2973*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(9*MW2 - Sf2974) + Sf2976) - 
     &            Sf2977) + Sf2975*MSq2Diff(2,1)**2))

	Sf2826 = 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-5*Sf2821*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(9*MW2 - Sf2822) + Sf2824) - 
     &            Sf2825) + Sf2823*MSq2Diff(1,1)**2))

	Sf325 = -6*MSl2Diff(1,1)**6 + 
     &   MW2*(3*Sf321*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(6*MW2 - 12*Sf320) + Sf322) - 
     &            Sf324) - Sf323*MSl2Diff(1,1)**2))

	Sf1852 = 30*MSl2Diff(2,1)**6 + 
     &   MW2*(-3*Sf1848*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(6*MW2 - 24*Sf1847) + 3*Sf1849) - 
     &            Sf1851) + 3*Sf1850*MSl2Diff(2,1)**2))

	Sf1285 = -6*MSl2Diff(2,1)**6 + 
     &   MW2*(3*Sf1281*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(6*MW2 - 12*Sf1280) + Sf1282) - 
     &            Sf1284) - Sf1283*MSl2Diff(2,1)**2))

	Sf1073 = 30*MSl2Diff(1,1)**6 + 
     &   MW2*(-3*Sf1069*MSl2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*(6*MW2 - 24*Sf1068) + 3*Sf1070) - 
     &            Sf1072) + 3*Sf1071*MSl2Diff(1,1)**2))

	Sf6859 = 4*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf6855*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(5*MW2 - 3*Sf6854) + Sf6856) - 
     &            Sf6858) + Sf6857*MSq2Diff(2,2)**2))

	Sf6745 = 4*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf6741*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(5*MW2 - 3*Sf6740) + Sf6742) - 
     &            Sf6744) + Sf6743*MSq2Diff(2,1)**2))

	Sf6604 = 4*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf6600*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(5*MW2 - 3*Sf6599) + Sf6601) - 
     &            Sf6603) + Sf6602*MSq2Diff(1,2)**2))

	Sf6478 = 4*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf6474*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(5*MW2 - 3*Sf6473) + Sf6475) - 
     &            Sf6477) + Sf6476*MSq2Diff(1,1)**2))

	Sf407 = 3*MSl2Diff(1,1)**6 + 
     &   MW2*(-(Sf403*MSl2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2 - 2*Sf402) + Sf405) - Sf406) + 
     &         Sf404*MSl2Diff(1,1)**2))

	Sf2070 = 3*MSl2Diff(1,1)**6 + 
     &   MW2*(-(Sf2066*MSl2Diff(1,1)**4) + 
     &      MW2*(MW2*(-Sf2069 + 
     &            MW2*(Sf2067 + MW2*(MW2 - 22*Sf455))) + 
     &         Sf2068*MSl2Diff(1,1)**2))

	Sf2150 = 3*MSl2Diff(2,1)**6 + 
     &   MW2*(-(Sf2146*MSl2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2 - 22*Sf1372) + Sf2147) - 
     &            Sf2149) + Sf2148*MSl2Diff(2,1)**2))

	Sf1344 = 3*MSl2Diff(2,1)**6 + 
     &   MW2*(-(Sf1340*MSl2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*(MW2*(MW2 - 2*Sf1339) + Sf1342) - 
     &            Sf1343) + Sf1341*MSl2Diff(2,1)**2))

	Sf5152 = MSq2Diff(2,2)**5 - MW2*(MW2 - 2*MStop2(2))*MStop2(2)**3

	Sf2745 = MSbot2(2)**5 - 5*MSbot2(2)**4*MStop2(1) + 
     &   10*MSbot2(2)**3*MStop2(1)**2 - 
     &   10*MSbot2(2)**2*MStop2(1)**3 - 
     &   (MW2 - MStop2(1))**2*MStop2(1)**3 + 
     &   5*MSbot2(2)*MStop2(1)**4

	Sf5033 = -(Sf5029**2*MStop2(2)**3) + 
     &   MSbot2(1)*(5*MStop2(2)**4 + 
     &      MSbot2(1)*(-10*MStop2(2)**3 + 
     &         MSbot2(1)*(MSbot2(1)*(MSbot2(1) - 5*MStop2(2)) + 
     &            10*MStop2(2)**2)))

	Sf2568 = -(Sf2563**2*MStop2(1)**3) + 
     &   MSbot2(1)*(5*MStop2(1)**4 + 
     &      MSbot2(1)*(-10*MStop2(1)**3 + 
     &         MSbot2(1)*(MSbot2(1)*(MSbot2(1) - 5*MStop2(1)) + 
     &            10*MStop2(1)**2)))

	Sf8383 = 96*MSq2Diff(2,2)**6 + 
     &   MW2*(-2*Sf8379*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(-6*Sf8382 + MW2*(Sf8380 - 12*MW2*Sf99)) + 
     &         Sf8381*MSq2Diff(2,2)**2))

	Sf8287 = 96*MSq2Diff(1,2)**6 + 
     &   MW2*(-2*Sf8283*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(-6*Sf8286 + MW2*(Sf8284 - 12*MW2*Sf91)) + 
     &         Sf8285*MSq2Diff(1,2)**2))

	Sf8094 = -256*MSq2Diff(2,2)**6 + 
     &   MW2*(2*Sf8090*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(2*Sf8093 + MW2*(Sf8092 - 84*MW2*Sf99)) - 
     &         3*Sf8091*MSq2Diff(2,2)**2))

	Sf8022 = -160*MSq2Diff(1,2)**6 + 
     &   MW2*(Sf8018*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(MW2*(54*MW2*Sf8017 - Sf8019) + Sf8021) - 
     &         3*Sf8020*MSq2Diff(1,2)**2))

	Sf8010 = -256*MSq2Diff(1,2)**6 + 
     &   MW2*(2*Sf8006*MSq2Diff(1,2)**4 + 
     &      MW2*(MW2*(Sf8009 + MW2*(Sf8008 - 84*MW2*Sf91)) - 
     &         3*Sf8007*MSq2Diff(1,2)**2))

	Sf7885 = 256*MSq2Diff(2,1)**6 + 
     &   MW2*(-2*Sf7881*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(84*MW2*Sf77 - Sf7883) - 2*Sf7884) + 
     &         3*Sf7882*MSq2Diff(2,1)**2))

	Sf8359 = -96*MSq2Diff(2,1)**6 + 
     &   MW2*(2*Sf8355*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(12*MW2*Sf77 - Sf8356) + 6*Sf8358) - 
     &         Sf8357*MSq2Diff(2,1)**2))

	Sf7807 = 256*MSq2Diff(1,1)**6 + 
     &   MW2*(-2*Sf7803*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(84*MW2*Sf69 - Sf7805) - 2*Sf7806) + 
     &         3*Sf7804*MSq2Diff(1,1)**2))

	Sf8263 = -96*MSq2Diff(1,1)**6 + 
     &   MW2*(2*Sf8259*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(12*MW2*Sf69 - Sf8260) + 6*Sf8262) - 
     &         Sf8261*MSq2Diff(1,1)**2))

	Sf8089 = -160*MSq2Diff(2,2)**6 + 
     &   MW2*(Sf8085*MSq2Diff(2,2)**4 + 
     &      MW2*(MW2*(MW2*(54*MW2*Sf4877 - Sf8086) + Sf8088) - 
     &         3*Sf8087*MSq2Diff(2,2)**2))

	Sf7880 = -160*MSq2Diff(2,1)**6 + 
     &   MW2*(Sf7876*MSq2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(54*MW2*Sf2280 - Sf7877) + Sf7879) - 
     &         3*Sf7878*MSq2Diff(2,1)**2))

	Sf7802 = -160*MSq2Diff(1,1)**6 + 
     &   MW2*(Sf7798*MSq2Diff(1,1)**4 + 
     &      MW2*(MW2*(MW2*(54*MW2*Sf2217 - Sf7799) + Sf7801) - 
     &         3*Sf7800*MSq2Diff(1,1)**2))

	Sf1910 = 5*MSl2Diff(1,1)**6 + 
     &   MW2*(-(Sf1899*MSl2Diff(1,1)**4) + 
     &      MW2*(MW2*(-Sf1909 + MW2*(Sf1908 - MW2*Sf320)) + 
     &         Sf1901*MSl2Diff(1,1)**2))

	Sf1970 = -5*MSl2Diff(2,1)**6 + 
     &   MW2*(Sf1961*MSl2Diff(2,1)**4 + 
     &      MW2*(MW2*(MW2*(MW2*Sf1280 - Sf1968) + Sf1969) - 
     &         Sf1963*MSl2Diff(2,1)**2))

	Sf8388 = -48*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(32*MW2 - Sf8384) + Sf8386) - 
     &         Sf8387) + Sf8385*MSq2Diff(2,2)**3)

	Sf8364 = -48*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(32*MW2 - Sf8360) + Sf8362) - 
     &         Sf8363) + Sf8361*MSq2Diff(2,1)**3)

	Sf8292 = -48*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(32*MW2 - Sf8288) + Sf8290) - 
     &         Sf8291) + Sf8289*MSq2Diff(1,2)**3)

	Sf8268 = -48*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(32*MW2 - Sf8264) + Sf8266) - 
     &         Sf8267) + Sf8265*MSq2Diff(1,1)**3)

	Sf806 = -6*MSl2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(24*MW2 - 6*Sf803) - Sf804) + 
     &         6*Sf805) + 12*Sf802*MSl2Diff(1,1)**3)

	Sf1631 = -6*MSl2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(24*MW2 - 6*Sf1628) - Sf1629) + 
     &         6*Sf1630) + 12*Sf1627*MSl2Diff(2,1)**3)

	Sf8692 = 48*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(Sf8452 + 
     &         MW2*(-Sf8451 + MW2*(14*MW2 + Sf8691))) - 
     &      Sf8385*MSq2Diff(2,2)**3)

	Sf8453 = 48*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(14*MW2 + 5*Sf8450) - Sf8451) + 
     &         Sf8452) - Sf8385*MSq2Diff(2,2)**3)

	Sf8684 = 48*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(Sf8441 + 
     &         MW2*(-Sf8440 + MW2*(14*MW2 + 5*Sf8683))) - 
     &      Sf8289*MSq2Diff(1,2)**3)

	Sf8442 = 48*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(14*MW2 + Sf8439) - Sf8440) + 
     &         Sf8441) - Sf8289*MSq2Diff(1,2)**3)

	Sf8652 = 48*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(Sf8428 + 
     &         MW2*(-Sf8427 + MW2*(14*MW2 + Sf8651))) - 
     &      Sf8361*MSq2Diff(2,1)**3)

	Sf8429 = 48*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(14*MW2 + 5*Sf8426) - Sf8427) + 
     &         Sf8428) - Sf8361*MSq2Diff(2,1)**3)

	Sf8644 = 48*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(Sf8417 + 
     &         MW2*(-Sf8416 + MW2*(14*MW2 + Sf8643))) - 
     &      Sf8265*MSq2Diff(1,1)**3)

	Sf8418 = 48*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(14*MW2 + 5*Sf8415) - Sf8416) + 
     &         Sf8417) - Sf8265*MSq2Diff(1,1)**3)

	Sf810 = -3*MSl2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(12*MW2 - 3*Sf803) - Sf808) + Sf809) + 
     &      6*Sf802*MSl2Diff(1,1)**3)

	Sf7111 = 3*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(-6*Sf7098 + MW2*(12*MW2 - Sf7109)) + 
     &         Sf7110) - 6*Sf6987*MSq2Diff(2,2)**3)

	Sf7073 = 3*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(-6*Sf7064 + MW2*(12*MW2 - Sf7071)) + 
     &         Sf7072) - 6*Sf6983*MSq2Diff(1,2)**3)

	Sf6885 = -3*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(6*Sf6873 + MW2*(12*MW2 - Sf6883)) - 
     &         Sf6884) + 6*Sf6872*MSq2Diff(2,2)**3)

	Sf6771 = -3*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(6*Sf6759 + MW2*(12*MW2 - Sf6769)) - 
     &         Sf6770) + 6*Sf6758*MSq2Diff(2,1)**3)

	Sf6634 = -3*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(6*Sf6622 + MW2*(12*MW2 - Sf6632)) - 
     &         Sf6633) + 6*Sf6621*MSq2Diff(1,2)**3)

	Sf6516 = -3*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(6*Sf6504 + MW2*(12*MW2 - Sf6514)) - 
     &         Sf6515) + 6*Sf6503*MSq2Diff(1,1)**3)

	Sf4680 = 3*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(-6*Sf4667 + MW2*(12*MW2 - Sf4678)) + 
     &         Sf4679) - 6*Sf4278*MSq2Diff(2,1)**3)

	Sf4575 = 3*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(-6*Sf4562 + MW2*(12*MW2 - Sf4573)) + 
     &         Sf4574) - 6*Sf4270*MSq2Diff(1,1)**3)

	Sf2170 = 3*MSl2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(6*Sf2164 + MW2*(12*MW2 - Sf2168)) - 
     &         Sf2169) - 6*Sf1495*MSl2Diff(2,1)**3)

	Sf2101 = 3*MSl2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(6*Sf2095 + MW2*(12*MW2 - Sf2099)) - 
     &         Sf2100) - 6*Sf595*MSl2Diff(1,1)**3)

	Sf1635 = -3*MSl2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(12*MW2 - 3*Sf1628) - Sf1633) + 
     &         Sf1634) + 6*Sf1627*MSl2Diff(2,1)**3)

	Sf5400 = MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf5397) + 2*Sf5398) - 
     &         Sf5399) - 2*Sf5396*MSq2Diff(2,2)**3)

	Sf5366 = -MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf5363) + 2*Sf5364) - 
     &         Sf5365) + 2*Sf5362*MSq2Diff(2,2)**3)

	Sf5294 = MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf5291) + 2*Sf5292) - 
     &         Sf5293) - 2*Sf5290*MSq2Diff(1,2)**3)

	Sf5258 = -MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf5255) + 2*Sf5256) - 
     &         Sf5257) + 2*Sf5254*MSq2Diff(1,2)**3)

	Sf434 = -MSl2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf431) + 2*Sf432) - Sf433) + 
     &      2*Sf430*MSl2Diff(1,1)**3)

	Sf392 = MSl2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf389) + 2*Sf390) - Sf391) - 
     &      2*Sf388*MSl2Diff(1,1)**3)

	Sf2992 = MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf2989) + 2*Sf2990) - 
     &         Sf2991) - 2*Sf2988*MSq2Diff(2,1)**3)

	Sf2955 = -MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf2952) + 2*Sf2953) - 
     &         Sf2954) + 2*Sf2951*MSq2Diff(2,1)**3)

	Sf2854 = MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf2851) + 2*Sf2852) - 
     &         Sf2853) - 2*Sf2850*MSq2Diff(1,1)**3)

	Sf2802 = -MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf2799) + 2*Sf2800) - 
     &         Sf2801) + 2*Sf2798*MSq2Diff(1,1)**3)

	Sf1357 = -MSl2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf1354) + 2*Sf1355) - 
     &         Sf1356) + 2*Sf1353*MSl2Diff(2,1)**3)

	Sf1329 = MSl2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(10*MW2 - Sf1326) + 2*Sf1327) - 
     &         Sf1328) - 2*Sf1325*MSl2Diff(2,1)**3)

	Sf7100 = MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf7097) - 2*Sf7098) + 
     &         Sf7099) - 2*Sf6987*MSq2Diff(2,2)**3)

	Sf7065 = MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(Sf7059 + 
     &         MW2*(MW2*(4*MW2 - Sf7057) - 2*Sf7064)) - 
     &      2*Sf6983*MSq2Diff(1,2)**3)

	Sf7060 = MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf7057) - Sf7058) + Sf7059) - 
     &      2*Sf6983*MSq2Diff(1,2)**3)

	Sf6875 = -MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf6871) + 2*Sf6873) - 
     &         Sf6874) + 2*Sf6872*MSq2Diff(2,2)**3)

	Sf6761 = -MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf6757) + 2*Sf6759) - 
     &         Sf6760) + 2*Sf6758*MSq2Diff(2,1)**3)

	Sf6624 = -MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf6620) + 2*Sf6622) - 
     &         Sf6623) + 2*Sf6621*MSq2Diff(1,2)**3)

	Sf6506 = -MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf6502) + 2*Sf6504) - 
     &         Sf6505) + 2*Sf6503*MSq2Diff(1,1)**3)

	Sf4669 = MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf4666) - 2*Sf4667) + 
     &         Sf4668) - 2*Sf4278*MSq2Diff(2,1)**3)

	Sf4564 = MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf4561) - 2*Sf4562) + 
     &         Sf4563) - 2*Sf4270*MSq2Diff(1,1)**3)

	Sf2166 = MSl2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf2163) + 2*Sf2164) - 
     &         Sf2165) - 2*Sf1495*MSl2Diff(2,1)**3)

	Sf2097 = MSl2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*(4*MW2 - Sf2094) + 2*Sf2095) - 
     &         Sf2096) - 2*Sf595*MSl2Diff(1,1)**3)

	Sf1481 = -35*MStau2(2)**4 + 
     &   MSneu2(1)*(-2*Sf1480*MStau2(2)**2 + 
     &      MSneu2(1)*(MSneu2(1)*(25*MSneu2(1) - 14*MStau2(2)) - 
     &         2*Sf1426*MStau2(2)))

	Sf580 = -35*MStau2(1)**4 + 
     &   MSneu2(1)*(-2*Sf579*MStau2(1)**2 + 
     &      MSneu2(1)*(MSneu2(1)*(25*MSneu2(1) - 14*MStau2(1)) - 
     &         2*Sf515*MStau2(1)))

	Sf1819 = -(Sf1196**3*MStau2(2)) + 
     &   MSneu2(1)*(-4*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(MSneu2(1) - 4*MStau2(2)) + 
     &         6*MStau2(2)**2))

	Sf1009 = -(Sf167**3*MStau2(1)) + 
     &   MSneu2(1)*(-4*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(MSneu2(1) - 4*MStau2(1)) + 
     &         6*MStau2(1)**2))

	Sf5153 = -(Sf5029**2*MStop2(2)**2) + 
     &   MSbot2(2)*(4*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(5*MSbot2(2) - 12*MStop2(2)) + 
     &         6*MStop2(2)**2))

	Sf2746 = -(Sf2563**2*MStop2(1)**2) + 
     &   MSbot2(2)*(4*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(5*MSbot2(2) - 12*MStop2(1)) + 
     &         6*MStop2(1)**2))

	Sf5160 = -7*Sf5029**2*MStop2(2)**2 + 
     &   MSbot2(2)*(12*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(3*MSbot2(2) - 4*MStop2(2)) - 
     &         6*MStop2(2)**2))

	Sf2753 = -7*Sf2563**2*MStop2(1)**2 + 
     &   MSbot2(2)*(12*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(3*MSbot2(2) - 4*MStop2(1)) - 
     &         6*MStop2(1)**2))

	Sf5034 = -(Sf5029**2*MStop2(2)**2) + 
     &   MSbot2(1)*(4*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(5*MSbot2(1) - 12*MStop2(2)) + 
     &         6*MStop2(2)**2))

	Sf2569 = -(Sf2563**2*MStop2(1)**2) + 
     &   MSbot2(1)*(4*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(5*MSbot2(1) - 12*MStop2(1)) + 
     &         6*MStop2(1)**2))

	Sf5042 = -7*Sf5029**2*MStop2(2)**2 + 
     &   MSbot2(1)*(12*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(3*MSbot2(1) - 4*MStop2(2)) - 
     &         6*MStop2(2)**2))

	Sf2577 = -7*Sf2563**2*MStop2(1)**2 + 
     &   MSbot2(1)*(12*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(3*MSbot2(1) - 4*MStop2(1)) - 
     &         6*MStop2(1)**2))

	Sf10032 = -72*Sf9934*MStau2(2)**4 + 
     &   MZ2*(4*Sf9933*MStau2(2)**3 + 
     &      MZ2*(Sf9932*MStau2(2)**2 + 
     &         MZ2*(25*MZ2*Ss49 - 2*Sf9931*MStau2(2))))

	Sf9935 = -72*Sf9934*MStau2(1)**4 + 
     &   MZ2*(4*Sf9933*MStau2(1)**3 + 
     &      MZ2*(Sf9932*MStau2(1)**2 + 
     &         MZ2*(25*MZ2*Ss49 - 2*Sf9931*MStau2(1))))

	Sf8244 = 16*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*Sf8240 - Sf8242) + 3*Sf8243) - 
     &      Sf8241*MSq2Diff(1,2)**3)

	Sf8334 = 16*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*Sf8228 - Sf8230) + Sf8333) - 
     &      Sf8229*MSq2Diff(1,1)**3)

	Sf8232 = 16*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*Sf8228 - Sf8230) + 3*Sf8231) - 
     &      Sf8229*MSq2Diff(1,1)**3)

	Sf8338 = 16*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*Sf8211 - Sf8213) + Sf8337) - 
     &      Sf8212*MSq2Diff(2,2)**3)

	Sf8215 = 16*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(MW2*(MW2*Sf8211 - Sf8213) + 3*Sf8214) - 
     &      Sf8212*MSq2Diff(2,2)**3)

	Sf8336 = 16*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*Sf8199 - Sf8201) + Sf8335) - 
     &      Sf8200*MSq2Diff(2,1)**3)

	Sf8203 = 16*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(MW2*(MW2*Sf8199 - Sf8201) + 3*Sf8202) - 
     &      Sf8200*MSq2Diff(2,1)**3)

	Sf8180 = -16*MSq2Diff(2,2)**5 + 
     &   MW2*(5*Sf5346*MSq2Diff(2,2)**3 + 
     &      MW2*(-3*Sf8179 + MW2*(Sf8178 + 7*MW2*MStop2(2))))

	Sf8160 = -16*MSq2Diff(1,2)**5 + 
     &   MW2*(5*Sf5238*MSq2Diff(1,2)**3 + 
     &      MW2*(-Sf8159 + MW2*(Sf8158 + 7*MW2*MStop2(2))))

	Sf8146 = 16*MSq2Diff(2,2)**6 + 
     &   MW2*(-(Sf8142*MSq2Diff(2,2)**4) + 
     &      MW2*(MW2*(MW2*Sf8143 - Sf8145) + 
     &         3*Sf8144*MSq2Diff(2,2)**2))

	Sf8068 = 16*MSq2Diff(1,2)**6 + 
     &   MW2*(-(Sf8064*MSq2Diff(1,2)**4) + 
     &      MW2*(MW2*(MW2*Sf8065 - Sf8067) + 
     &         3*Sf8066*MSq2Diff(1,2)**2))

	Sf7971 = -16*MSq2Diff(2,1)**5 + 
     &   MW2*(5*Sf2935*MSq2Diff(2,1)**3 + 
     &      MW2*(-3*Sf7970 + MW2*(Sf7969 + 7*MW2*MStop2(1))))

	Sf7951 = -16*MSq2Diff(1,1)**5 + 
     &   MW2*(5*Sf2769*MSq2Diff(1,1)**3 + 
     &      MW2*(-3*Sf7950 + MW2*(Sf7949 + 7*MW2*MStop2(1))))

	Sf7937 = 16*MSq2Diff(2,1)**6 + 
     &   MW2*(-(Sf7933*MSq2Diff(2,1)**4) + 
     &      MW2*(MW2*(MW2*Sf7934 - Sf7936) + 
     &         3*Sf7935*MSq2Diff(2,1)**2))

	Sf7859 = 16*MSq2Diff(1,1)**6 + 
     &   MW2*(-(Sf7855*MSq2Diff(1,1)**4) + 
     &      MW2*(MW2*(MW2*Sf7856 - Sf7858) + 
     &         3*Sf7857*MSq2Diff(1,1)**2))

	Sf5458 = -3*MSq2Diff(2,2)**5 + 
     &   MW2*(MW2*(-3*Sf5457 + MW2*(3*Sf5456 - 2*MW2*MSbot2(2))) + 
     &      Sf5455*MSq2Diff(2,2)**3)

	Sf5451 = -3*MSq2Diff(1,2)**5 + 
     &   MW2*(MW2*(-3*Sf5450 + MW2*(Sf5449 - 2*MW2*MSbot2(1))) + 
     &      Sf5448*MSq2Diff(1,2)**3)

	Sf4917 = -(Sf4912*MSq2Diff(2,2)**4) + 
     &   MW2*(Sf4916*MSq2Diff(2,2)**2 + 
     &      MW2*(-Sf4914 + MW2*(Sf4913 - 2*MW2*MStop2(2))))

	Sf4826 = -(Sf4821*MSq2Diff(1,2)**4) + 
     &   MW2*(Sf4825*MSq2Diff(1,2)**2 + 
     &      MW2*(-Sf4823 + MW2*(Sf4822 - 2*MW2*MStop2(2))))

	Sf3565 = 7*MSbot2(2)**5 + 
     &   MW2*(-20*MSbot2(2)**4 + 
     &      MW2*(-12*MSbot2(2)**3 + 
     &         MW2*(MW2*Sf3210 + 2*Sf3564*MSbot2(2))))

	Sf3212 = 7*MSbot2(1)**5 + 
     &   MW2*(-20*MSbot2(1)**4 + 
     &      MW2*(-12*MSbot2(1)**3 + 
     &         MW2*(MW2*Sf3210 + 2*Sf3211*MSbot2(1))))

	Sf3070 = -3*MSq2Diff(2,1)**5 + 
     &   MW2*(MW2*(-3*Sf3069 + MW2*(3*Sf3068 - 2*MW2*MSbot2(2))) + 
     &      Sf3067*MSq2Diff(2,1)**3)

	Sf3062 = -3*MSq2Diff(1,1)**5 + 
     &   MW2*(MW2*(-3*Sf3061 + MW2*(3*Sf3060 - 2*MW2*MSbot2(1))) + 
     &      Sf3059*MSq2Diff(1,1)**3)

	Sf2346 = -(Sf2341*MSq2Diff(2,1)**4) + 
     &   MW2*(Sf2345*MSq2Diff(2,1)**2 + 
     &      MW2*(-Sf2343 + MW2*(Sf2342 - 2*MW2*MStop2(1))))

	Sf2333 = -(Sf2328*MSq2Diff(1,1)**4) + 
     &   MW2*(Sf2332*MSq2Diff(1,1)**2 + 
     &      MW2*(-Sf2330 + MW2*(Sf2329 - 2*MW2*MStop2(1))))

	tmp17 = MZ2**3*(-112*MZ2 + 62*MHiggs2(2))

	Sf9789 = Sf9788*MStau2(2)**3 + 
     &   MHiggs2(2)*(-(MZ2*Sf9549*MStau2(2)**2) + 
     &      MHiggs2(2)*(tmp17*MHiggs2(2) + 
     &         3*MZ2**2*Sf9548*MStau2(2)))

	Sf9555 = Sf9554*MStau2(1)**3 + 
     &   MHiggs2(2)*(-(MZ2*Sf9549*MStau2(1)**2) + 
     &      MHiggs2(2)*(tmp17*MHiggs2(2) + 
     &         3*MZ2**2*Sf9548*MStau2(1)))

	tmp16 = MZ2**3*(-56*MZ2 + 38*MHiggs2(2))

	Sf9792 = -(Sf9791*MStau2(2)**3) + 
     &   MHiggs2(2)*(-2*MZ2*Sf9557*MStau2(2)**2 + 
     &      MHiggs2(2)*(tmp16*MHiggs2(2) - MZ2**2*Sf9556*MStau2(2))
     &      )

	Sf9563 = -(Sf9562*MStau2(1)**3) + 
     &   MHiggs2(2)*(-2*MZ2*Sf9557*MStau2(1)**2 + 
     &      MHiggs2(2)*(tmp16*MHiggs2(2) - MZ2**2*Sf9556*MStau2(1))
     &      )

	tmp15 = MZ2**2*(52*MZ2 + 52*MHiggs2(2))

	Sf9786 = -(Sf9785*MStau2(2)**2) + 
     &   MHiggs2(2)*(MHiggs2(2)*(-196*MZ2**4 + tmp15*MHiggs2(2)) + 
     &      MZ2*Sf9564*MStau2(2))

	Sf9570 = -(Sf9569*MStau2(1)**2) + 
     &   MHiggs2(2)*(MHiggs2(2)*(-196*MZ2**4 + tmp15*MHiggs2(2)) + 
     &      MZ2*Sf9564*MStau2(1))

	tmp14 = MZ2*(196*MZ2 + 66*MHiggs2(2))

	tmp28 = -234*MZ2**4 + 
     &   MHiggs2(2)*(-168*MZ2**3 + tmp14*MHiggs2(2))

	Sf9783 = tmp28*MHiggs2(2) - Sf9782*MStau2(2)

	Sf9576 = tmp28*MHiggs2(2) - Sf9575*MStau2(1)

	Sf8177 = 32*MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(13*MW2 - 2*Sf8174) + 3*Sf8176) - 
     &      2*Sf8175*MSq2Diff(2,2)**2)

	Sf8157 = 32*MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(13*MW2 - 2*Sf8154) + 3*Sf8156) - 
     &      2*Sf8155*MSq2Diff(1,2)**2)

	Sf7968 = 32*MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(13*MW2 - 2*Sf7965) + 3*Sf7967) - 
     &      2*Sf7966*MSq2Diff(2,1)**2)

	Sf7948 = 32*MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(13*MW2 - 2*Sf7945) + 3*Sf7947) - 
     &      2*Sf7946*MSq2Diff(1,1)**2)

	Sf896 = MSl2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(10*MW2 - 25*Sf34) + 7*Sf895) - 
     &      7*Sf34*MSl2Diff(1,1)**2)

	Sf1695 = MSl2Diff(2,1)**4 + 
     &   MW2*(MW2*(7*Sf1694 + MW2*(10*MW2 - 25*Sf47)) - 
     &      7*Sf47*MSl2Diff(2,1)**2)

	Sf6781 = 14*MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(6*MW2 - 3*Sf6778) + 2*Sf6780) - 
     &      Sf6779*MSq2Diff(2,2)**2)

	Sf6526 = 14*MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(6*MW2 - 3*Sf6523) + 2*Sf6525) - 
     &      Sf6524*MSq2Diff(1,2)**2)

	Sf4285 = 14*MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(6*MW2 - 3*Sf4282) + 2*Sf4284) - 
     &      Sf4283*MSq2Diff(2,1)**2)

	Sf4277 = 14*MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(6*MW2 - 3*Sf4274) + 2*Sf4276) - 
     &      Sf4275*MSq2Diff(1,1)**2)

	Sf1835 = MW2*(MW2*(MW2*(3*MW2 - 3*Sf1832) + Sf1833) - 
     &      Sf1834) + 8*MSl2Diff(2,1)**3*MStau2(2)

	Sf1046 = MW2*(MW2*(MW2*(3*MW2 - 3*Sf1043) + Sf1044) - 
     &      Sf1045) + 8*MSl2Diff(1,1)**3*MStau2(1)

	Sf6990 = -MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6988) + Sf6989) - 
     &      Sf6987*MSq2Diff(2,2)**2)

	Sf6986 = -MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6984) + Sf6985) - 
     &      Sf6983*MSq2Diff(1,2)**2)

	Sf6777 = -MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6775) + Sf6776) + 
     &      Sf5149*MSq2Diff(2,2)**2)

	Sf6695 = -MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6693) + Sf6694) + 
     &      Sf2742*MSq2Diff(2,1)**2)

	Sf6522 = -MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6520) + Sf6521) + 
     &      Sf5030*MSq2Diff(1,2)**2)

	Sf6422 = -MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6420) + Sf6421) + 
     &      Sf2565*MSq2Diff(1,1)**2)

	Sf6196 = 6*MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6193) + 2*Sf6195) - 
     &      Sf6194*MSq2Diff(2,2)**2)

	Sf6192 = 3*MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6189) + Sf6191) - 
     &      Sf6190*MSq2Diff(2,2)**2)

	Sf6110 = 3*MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf6107) + Sf6109) - 
     &      Sf6108*MSq2Diff(2,1)**2)

	Sf5818 = 6*MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf5815) + 2*Sf5817) - 
     &      Sf5816*MSq2Diff(1,2)**2)

	Sf5814 = 3*MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf5811) + Sf5813) - 
     &      Sf5812*MSq2Diff(1,2)**2)

	Sf5721 = 3*MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf5718) + Sf5720) - 
     &      Sf5719*MSq2Diff(1,1)**2)

	Sf4910 = MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf4907) + 3*Sf4909) - 
     &      Sf4908*MSq2Diff(2,2)**2)

	Sf4903 = MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf4900) + Sf4902) - 
     &      Sf4901*MSq2Diff(1,2)**2)

	Sf4281 = -MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf4279) + Sf4280) - 
     &      Sf4278*MSq2Diff(2,1)**2)

	Sf4273 = -MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf4271) + Sf4272) - 
     &      Sf4270*MSq2Diff(1,1)**2)

	Sf4003 = 6*MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf4000) + 2*Sf4002) - 
     &      Sf4001*MSq2Diff(2,1)**2)

	Sf3996 = 6*MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf3993) + 2*Sf3995) - 
     &      Sf3994*MSq2Diff(1,1)**2)

	Sf512 = MSl2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - 9*Sf34) + Sf511) - 
     &      7*Sf34*MSl2Diff(1,1)**2)

	Sf2338 = MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf2335) + 3*Sf2337) - 
     &      Sf2336*MSq2Diff(2,1)**2)

	Sf2324 = MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*(2*MW2 - Sf2321) + 3*Sf2323) - 
     &      Sf2322*MSq2Diff(1,1)**2)

	Sf1423 = MSl2Diff(2,1)**4 + 
     &   MW2*(MW2*(Sf1422 + MW2*(2*MW2 - 9*Sf47)) - 
     &      7*Sf47*MSl2Diff(2,1)**2)

	Sf5204 = MHiggs2(4)*
     &    (-(MW2*Sf5203) + 
     &      MHiggs2(4)*(6*Sf5202 + 
     &         MHiggs2(4)*(Sf5201 + 2*MHiggs2(4)))) + 
     &   4*MW2**2*MSq2Diff(2,2)**2

	Sf5192 = MHiggs2(4)*
     &    (-(MW2*Sf5191) + 
     &      MHiggs2(4)*(6*Sf5190 + 
     &         MHiggs2(4)*(Sf5189 + 2*MHiggs2(4)))) + 
     &   4*MW2**2*MSq2Diff(1,2)**2

	Sf5184 = MHiggs2(4)*
     &    (-(MW2*Sf5183) + 
     &      MHiggs2(4)*(6*Sf5182 + 
     &         MHiggs2(4)*(Sf5181 + 2*MHiggs2(4)))) + 
     &   4*MW2**2*MSq2Diff(2,1)**2

	Sf5178 = MHiggs2(4)*
     &    (-(MW2*Sf5177) + 
     &      MHiggs2(4)*(6*Sf5176 + 
     &         MHiggs2(4)*(Sf5175 + 2*MHiggs2(4)))) + 
     &   4*MW2**2*MSq2Diff(1,1)**2

	Sf4891 = 2*MW2*MSq2Diff(2,2)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf4888 - Sf4887*MHiggs2(4) + MHiggs2(4)**2) - 
     &      Sf4890*MSq2Diff(2,2))

	Sf4857 = 2*MW2*MSq2Diff(1,2)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf4854 - Sf4853*MHiggs2(4) + MHiggs2(4)**2) - 
     &      Sf4856*MSq2Diff(1,2))

	Sf2302 = 2*MW2*MSq2Diff(2,1)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf2299 - Sf2298*MHiggs2(4) + MHiggs2(4)**2) - 
     &      Sf2301*MSq2Diff(2,1))

	Sf2251 = 2*MW2*MSq2Diff(1,1)**3 + 
     &   MHiggs2(4)*(MHiggs2(4)*
     &       (Sf2248 - Sf2247*MHiggs2(4) + MHiggs2(4)**2) - 
     &      Sf2250*MSq2Diff(1,1))

	Sf1985 = MHiggs2(4)*
     &    (-2*MW2*Sf1983 + 
     &      MHiggs2(4)*(Sf1984 - 4*Sf47*MHiggs2(4) + MHiggs2(4)**2)
     &      ) + 2*MW2**2*MSl2Diff(2,1)**2

	Sf1926 = MHiggs2(4)*
     &    (-2*MW2*Sf1924 + 
     &      MHiggs2(4)*(Sf1925 - 4*Sf34*MHiggs2(4) + MHiggs2(4)**2)
     &      ) + 2*MW2**2*MSl2Diff(1,1)**2

	tmp20 = -(CA*reMUE) + SA*Re(Atop)

	Sf7214 = 6*MTop2*SA*UStop21 + CAB*SB*Sf86*UStop21 + 
     &   3*MTop*tmp20*UStop22

	Sf7192 = 6*MTop2*SA*UStop11 + CAB*SB*Sf86*UStop11 + 
     &   3*MTop*tmp20*UStop12

	Sf4790 = (6*MBot2*SA - CB*SAB*Sf15)*USbot21 + 
     &   3*MBot*(CA*reMUE + SA*Re(Abot))*USbot22

	Sf4766 = (6*MBot2*SA - CB*SAB*Sf15)*USbot11 + 
     &   3*MBot*(CA*reMUE + SA*Re(Abot))*USbot12

	Sf1965 = Sf1904*MStau2(2)**3 + 
     &   MSneu2(1)*(-8*MStau2(2)**2 + 
     &      MSneu2(1)*(50*MSneu2(1) - Sf1903*MStau2(2)))

	Sf1975 = Sf1974*MStau2(2)**2 + 
     &   MSneu2(1)*(MSneu2(1)*(50*MSneu2(1) - 11*MStau2(2)) + 
     &      2*ML2*MStau2(2))

	Sf1905 = Sf1904*MStau2(1)**3 + 
     &   MSneu2(1)*(-8*MStau2(1)**2 + 
     &      MSneu2(1)*(50*MSneu2(1) - Sf1903*MStau2(1)))

	Sf1915 = Sf1914*MStau2(1)**2 + 
     &   MSneu2(1)*(MSneu2(1)*(50*MSneu2(1) - 11*MStau2(1)) + 
     &      2*ML2*MStau2(1))

	Sf1471 = -45*MStau2(2)**3 + 
     &   MSneu2(1)*(Sf190*MStau2(2)**2 + 
     &      MSneu2(1)*(25*MSneu2(1) - Sf569*MStau2(2)))

	Sf570 = -45*MStau2(1)**3 + 
     &   MSneu2(1)*(Sf190*MStau2(1)**2 + 
     &      MSneu2(1)*(25*MSneu2(1) - Sf569*MStau2(1)))

	Sf1498 = -6*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(3*MSneu2(1) - MStau2(2)) - 
     &      Sf1497*MStau2(2))

	Sf598 = -6*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(3*MSneu2(1) - MStau2(1)) - 
     &      Sf597*MStau2(1))

	Sf1840 = Sf1058*MStau2(2)**3 + 
     &   MSneu2(1)*(-(Sf1057*MStau2(2)**2) + 
     &      MSneu2(1)*(MSneu2(1) - Sf292*MStau2(2)))

	Sf1225 = -(Sf1185**2*MStau2(2)) + 
     &   MSneu2(1)*(MSneu2(1)*(MSneu2(1) - 3*MStau2(2)) + 
     &      3*MStau2(2)**2)

	Sf1059 = Sf1058*MStau2(1)**3 + 
     &   MSneu2(1)*(-(Sf1057*MStau2(1)**2) + 
     &      MSneu2(1)*(MSneu2(1) - Sf292*MStau2(1)))

	Sf238 = -(Sf155**2*MStau2(1)) + 
     &   MSneu2(1)*(MSneu2(1)*(MSneu2(1) - 3*MStau2(1)) + 
     &      3*MStau2(1)**2)

	Sf5151 = 10*Sf5029**2*MStop2(2) + 
     &   MSbot2(2)*(MSbot2(2)*(7*MSbot2(2) - 21*MStop2(2)) + 
     &      7*MStop2(2)**2)

	Sf2744 = 10*Sf2563**2*MStop2(1) + 
     &   MSbot2(2)*(MSbot2(2)*(7*MSbot2(2) - 21*MStop2(1)) + 
     &      7*MStop2(1)**2)

	Sf5618 = -28*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(14*MSbot2(1) - 6*MStop2(2)) + 
     &      2*Sf5617*MStop2(2))

	Sf5032 = 10*Sf5029**2*MStop2(2) + 
     &   MSbot2(1)*(MSbot2(1)*(7*MSbot2(1) - 21*MStop2(2)) + 
     &      7*MStop2(2)**2)

	Sf2567 = 10*Sf2563**2*MStop2(1) + 
     &   MSbot2(1)*(MSbot2(1)*(7*MSbot2(1) - 21*MStop2(1)) + 
     &      7*MStop2(1)**2)

	Sf6300 = -3*MSq2Diff(2,2)**4 + 
     &   MW2*(MW2*(MW2*Sf5140 - Sf6299) + Sf6298*MSq2Diff(2,2)**2)

	Sf6297 = -3*MSq2Diff(1,2)**4 + 
     &   MW2*(MW2*(MW2*Sf5021 - Sf6296) + Sf6295*MSq2Diff(1,2)**2)

	Sf3559 = -2*MSbot2(2)**4 + 
     &   MW2*(16*MSbot2(2)**3 + MW2*(MW2*Sf3558 - 54*MSbot2(2)**2))

	Sf5576 = 5*MStop2(2)**4 + 
     &   MW2*(-18*MStop2(2)**3 + 
     &      MW2*(9*MW2*Sf3254 - 2*Sf5575*MStop2(2)))

	Sf3256 = 5*MStop2(1)**4 + 
     &   MW2*(-18*MStop2(1)**3 + 
     &      MW2*(9*MW2*Sf3254 - 2*Sf3255*MStop2(1)))

	Sf3203 = -2*MSbot2(1)**4 + 
     &   MW2*(16*MSbot2(1)**3 + MW2*(MW2*Sf3202 - 54*MSbot2(1)**2))

	Sf5512 = -5*MStop2(2)**4 + 
     &   MW2*(18*MStop2(2)**3 + 
     &      MW2*(3*MW2*Sf3132 + 2*Sf5511*MStop2(2)))

	Sf3134 = -5*MStop2(1)**4 + 
     &   MW2*(18*MStop2(1)**3 + 
     &      MW2*(3*MW2*Sf3132 + 2*Sf3133*MStop2(1)))

	Sf4006 = -3*MSq2Diff(2,1)**4 + 
     &   MW2*(MW2*(MW2*Sf2732 - Sf4005) + Sf4004*MSq2Diff(2,1)**2)

	Sf3999 = -3*MSq2Diff(1,1)**4 + 
     &   MW2*(MW2*(MW2*Sf2554 - Sf3998) + Sf3997*MSq2Diff(1,1)**2)

	tmp27 = 2*MZ2**2 + 3*MStop2(2)**2

	Sf10592 = tmp27*MHiggs2(2) + 
     &   MZ2*(-MZ2**2 + 2*Sf8929*Sf8930*MStop2(2))

	tmp26 = 2*MZ2**2 + 3*MStop2(1)**2

	Sf10586 = tmp26*MHiggs2(2) + 
     &   MZ2*(-MZ2**2 + 2*Sf8929*Sf8930*MStop2(1))

	tmp25 = 2*MZ2**2 + 3*MStau2(2)**2

	Sf9868 = tmp25*MHiggs2(2) + 
     &   MZ2*(-MZ2**2 + 2*Sf8929*Sf8930*MStau2(2))

	tmp24 = 2*MZ2**2 + 3*MStau2(1)**2

	Sf9860 = tmp24*MHiggs2(2) + 
     &   MZ2*(-MZ2**2 + 2*Sf8929*Sf8930*MStau2(1))

	tmp23 = 2*MZ2**2 + 3*MSbot2(2)**2

	Sf9106 = tmp23*MHiggs2(2) + 
     &   MZ2*(-MZ2**2 + 2*Sf8929*Sf8930*MSbot2(2))

	tmp22 = 2*MZ2**2 + 3*MSbot2(1)**2

	Sf9098 = tmp22*MHiggs2(2) + 
     &   MZ2*(-MZ2**2 + 2*Sf8929*Sf8930*MSbot2(1))

	Sf10331 = tmp27*MHiggs2(1) + 
     &   MZ2*(-MZ2**2 + 2*Sf9131*Sf9132*MStop2(2))

	Sf10322 = tmp26*MHiggs2(1) + 
     &   MZ2*(-MZ2**2 + 2*Sf9131*Sf9132*MStop2(1))

	Sf10217 = tmp25*MHiggs2(1) + 
     &   MZ2*(-MZ2**2 + 2*Sf9131*Sf9132*MStau2(2))

	Sf10211 = tmp24*MHiggs2(1) + 
     &   MZ2*(-MZ2**2 + 2*Sf9131*Sf9132*MStau2(1))

	Sf9285 = tmp23*MHiggs2(1) + 
     &   MZ2*(-MZ2**2 + 2*Sf9131*Sf9132*MSbot2(2))

	Sf9279 = tmp22*MHiggs2(1) + 
     &   MZ2*(-MZ2**2 + 2*Sf9131*Sf9132*MSbot2(1))

	Sf5963 = 203*MSq2Diff(2,2)**6 + 
     &   MW2**2*(MW2**2*(230*MW2**2 + Sf5961) + 
     &      2*Sf5962*MSq2Diff(2,2)**2)

	Sf5537 = 203*MSq2Diff(1,2)**6 + 
     &   MW2**2*(MW2**2*(230*MW2**2 + Sf5535) + 
     &      2*Sf5536*MSq2Diff(1,2)**2)

	Sf3538 = 203*MSq2Diff(2,1)**6 + 
     &   MW2**2*(MW2**2*(230*MW2**2 + Sf3536) + 
     &      2*Sf3537*MSq2Diff(2,1)**2)

	Sf3176 = 203*MSq2Diff(1,1)**6 + 
     &   MW2**2*(MW2**2*(230*MW2**2 + Sf3174) + 
     &      2*Sf3175*MSq2Diff(1,1)**2)

	Sf5908 = 89*MSq2Diff(2,2)**6 + 
     &   MW2**2*(MW2**2*(2*MW2**2 + 5*Sf5906) + 
     &      2*Sf5907*MSq2Diff(2,2)**2)

	Sf5474 = 89*MSq2Diff(1,2)**6 + 
     &   MW2**2*(MW2**2*(2*MW2**2 + 5*Sf5472) + 
     &      2*Sf5473*MSq2Diff(1,2)**2)

	Sf3464 = 89*MSq2Diff(2,1)**6 + 
     &   MW2**2*(MW2**2*(2*MW2**2 + 5*Sf3462) + 
     &      2*Sf3463*MSq2Diff(2,1)**2)

	Sf3077 = 89*MSq2Diff(1,1)**6 + 
     &   MW2**2*(MW2**2*(2*MW2**2 + 5*Sf3075) + 
     &      2*Sf3076*MSq2Diff(1,1)**2)

	Sf8729 = MW2*(4*MZ2*Sf8727 + MW2*(8*MW2 - Sf8728)) - 
     &   9*MZ2**2*MSq2Diff(2,2)

	Sf8726 = MW2*(4*MZ2*Sf8724 + MW2*(8*MW2 - Sf8725)) - 
     &   9*MZ2**2*MSq2Diff(1,2)

	Sf8515 = MW2*(MW2*(8*MW2 - Sf8513) + 4*MZ2*Sf8514) + 
     &   9*MZ2**2*MSq2Diff(2,2)

	Sf8512 = MW2*(MW2*(8*MW2 - Sf8510) + 4*MZ2*Sf8511) + 
     &   9*MZ2**2*MSq2Diff(2,1)

	Sf8499 = MW2*(MW2*(8*MW2 - Sf8497) + 4*MZ2*Sf8498) + 
     &   9*MZ2**2*MSq2Diff(1,2)

	Sf8496 = MW2*(MW2*(8*MW2 - Sf8494) + 4*MZ2*Sf8495) + 
     &   9*MZ2**2*MSq2Diff(1,1)

	Sf7642 = MW2*(MW2*(8*MW2 - Sf7640) + 4*MZ2*Sf7641) - 
     &   9*MZ2**2*MSl2Diff(2,1)

	Sf7639 = MW2*(MW2*(8*MW2 - Sf7637) + 4*MZ2*Sf7638) - 
     &   9*MZ2**2*MSl2Diff(1,1)

	Sf94 = MW2*(MW2*(5*MW2 - Sf92) + Sf93) - 5*MZ2*MSq2Diff(1,2)**2

	Sf80 = MW2*(MW2*(5*MW2 - Sf78) + Sf79) - 5*MZ2*MSq2Diff(2,1)**2

	Sf72 = MW2*(MW2*(5*MW2 - Sf70) + Sf71) - 5*MZ2*MSq2Diff(1,1)**2

	Sf50 = MW2*(MW2*(5*MW2 - Sf48) + Sf49) - 5*MZ2*MSl2Diff(2,1)**2

	Sf1687 = MW2*(MW2*(5*MW2 - 11*Sf47) + 7*Sf47**2) - 
     &   Sf47*MSl2Diff(2,1)**2

	Sf37 = MW2*(MW2*(5*MW2 - Sf35) + Sf36) - 5*MZ2*MSl2Diff(1,1)**2

	Sf881 = MW2*(MW2*(5*MW2 - 11*Sf34) + 7*Sf34**2) - 
     &   Sf34*MSl2Diff(1,1)**2

	Sf102 = MW2*(MW2*(5*MW2 - Sf100) + Sf101) - 
     &   5*MZ2*MSq2Diff(2,2)**2

	Sf2129 = MW2*(MW2*(3*MW2 - Sf2127) + Sf2128) + MSl2Diff(2,1)**3

	Sf2049 = MW2*(MW2*(3*MW2 - Sf2047) + Sf2048) + MSl2Diff(1,1)**3

	Sf1987 = -2*MSl2Diff(2,1)**3 + 
     &   MW2*(-Sf1986 + MW2*(MW2 + 2*MStau2(2)))

	Sf1928 = -2*MSl2Diff(1,1)**3 + 
     &   MW2*(-Sf1927 + MW2*(MW2 + 2*MStau2(1)))

	tmp29 = 8*MW2**2 + MHiggs2(4)*(-Sf2236 + 12*MHiggs2(4))

	Sf2291 = tmp29*MHiggs2(4) + 3*Sf2201*MSbot2(2)**2

	Sf2237 = tmp29*MHiggs2(4) + 3*Sf2201*MSbot2(1)**2

	Sf4881 = MHiggs2(4)*
     &    (Sf4878 + MHiggs2(4)*(-Sf4876 + 3*MHiggs2(4))) - 
     &   Sf4880*MSq2Diff(2,2)

	Sf6774 = MHiggs2(4)*
     &    (Sf6773 + MHiggs2(4)*(-3*Sf4872 + 3*MHiggs2(4))) - 
     &   4*MW2*MSq2Diff(2,2)**2

	Sf4847 = MHiggs2(4)*
     &    (Sf4846 + MHiggs2(4)*(-Sf4842 + 3*MHiggs2(4))) - 
     &   Sf4844*MSq2Diff(1,2)

	Sf6519 = MHiggs2(4)*
     &    (Sf6518 + MHiggs2(4)*(-3*Sf4837 + 3*MHiggs2(4))) - 
     &   4*MW2*MSq2Diff(1,2)**2

	Sf2284 = MHiggs2(4)*
     &    (Sf2281 + MHiggs2(4)*(-Sf2279 + 3*MHiggs2(4))) - 
     &   Sf2283*MSq2Diff(2,1)

	Sf4269 = MHiggs2(4)*
     &    (Sf4268 + MHiggs2(4)*(-3*Sf2275 + 3*MHiggs2(4))) - 
     &   4*MW2*MSq2Diff(2,1)**2

	Sf2221 = MHiggs2(4)*
     &    (Sf2218 + MHiggs2(4)*(-Sf2216 + 3*MHiggs2(4))) - 
     &   Sf2220*MSq2Diff(1,1)

	Sf4267 = MHiggs2(4)*
     &    (Sf4266 + MHiggs2(4)*(-3*Sf2212 + 3*MHiggs2(4))) - 
     &   4*MW2*MSq2Diff(1,1)**2

	Sf5200 = MHiggs2(4)*
     &    (Sf5199 - 4*Sf99*MHiggs2(4) + MHiggs2(4)**2) - 
     &   2*MW2*MSq2Diff(2,2)**2

	Sf5188 = MHiggs2(4)*
     &    (Sf5187 - 4*Sf91*MHiggs2(4) + MHiggs2(4)**2) - 
     &   2*MW2*MSq2Diff(1,2)**2

	Sf3744 = MHiggs2(4)*
     &    (Sf3743 - 4*Sf77*MHiggs2(4) + MHiggs2(4)**2) - 
     &   2*MW2*MSq2Diff(2,1)**2

	Sf3735 = MHiggs2(4)*
     &    (Sf3734 - 4*Sf69*MHiggs2(4) + MHiggs2(4)**2) - 
     &   2*MW2*MSq2Diff(1,1)**2

	Sf10581 = Sf9081*MStop2(2)**2 + 
     &   MHiggs2(2)*(MHiggs2(2)*(-2*MZ2 + 3*MHiggs2(2)) + 
     &      2*Sf8929*Sf8930*MStop2(2))

	Sf10569 = Sf9081*MStop2(1)**2 + 
     &   MHiggs2(2)*(MHiggs2(2)*(-2*MZ2 + 3*MHiggs2(2)) + 
     &      2*Sf8929*Sf8930*MStop2(1))

	Sf9855 = Sf9081*MStau2(2)**2 + 
     &   MHiggs2(2)*(MHiggs2(2)*(-2*MZ2 + 3*MHiggs2(2)) + 
     &      2*Sf8929*Sf8930*MStau2(2))

	Sf9844 = Sf9081*MStau2(1)**2 + 
     &   MHiggs2(2)*(MHiggs2(2)*(-2*MZ2 + 3*MHiggs2(2)) + 
     &      2*Sf8929*Sf8930*MStau2(1))

	Sf9093 = Sf9081*MSbot2(2)**2 + 
     &   MHiggs2(2)*(MHiggs2(2)*(-2*MZ2 + 3*MHiggs2(2)) + 
     &      2*Sf8929*Sf8930*MSbot2(2))

	Sf9082 = Sf9081*MSbot2(1)**2 + 
     &   MHiggs2(2)*(MHiggs2(2)*(-2*MZ2 + 3*MHiggs2(2)) + 
     &      2*Sf8929*Sf8930*MSbot2(1))

	Sf10430 = Sf9261*MStop2(2)**2 + 
     &   MHiggs2(1)*(MHiggs2(1)*(-2*MZ2 + 3*MHiggs2(1)) + 
     &      2*Sf9131*Sf9132*MStop2(2))

	Sf10419 = Sf9261*MStop2(1)**2 + 
     &   MHiggs2(1)*(MHiggs2(1)*(-2*MZ2 + 3*MHiggs2(1)) + 
     &      2*Sf9131*Sf9132*MStop2(1))

	Sf10206 = Sf9261*MStau2(2)**2 + 
     &   MHiggs2(1)*(MHiggs2(1)*(-2*MZ2 + 3*MHiggs2(1)) + 
     &      2*Sf9131*Sf9132*MStau2(2))

	Sf10194 = Sf9261*MStau2(1)**2 + 
     &   MHiggs2(1)*(MHiggs2(1)*(-2*MZ2 + 3*MHiggs2(1)) + 
     &      2*Sf9131*Sf9132*MStau2(1))

	Sf9274 = Sf9261*MSbot2(2)**2 + 
     &   MHiggs2(1)*(MHiggs2(1)*(-2*MZ2 + 3*MHiggs2(1)) + 
     &      2*Sf9131*Sf9132*MSbot2(2))

	Sf9262 = Sf9261*MSbot2(1)**2 + 
     &   MHiggs2(1)*(MHiggs2(1)*(-2*MZ2 + 3*MHiggs2(1)) + 
     &      2*Sf9131*Sf9132*MSbot2(1))

	Sf7203 = Sf7199*Sf86*UStop11*UStop21**2 - 
     &   MZ2*(MZ2*Sf7202 - 
     &      4*SW2*MStop2(3)*UStop12*UStop21*UStop22)

	Sf6970 = 36*MTop2*Sf6966*TB*UStop21 + 
     &   Sf6969*TB*UStop21 - 36*MTop*Sf2245*Sf6966*UStop22

	Sf4904 = MTop2*TB*UStop21 + C2B*MW2*TB*UStop21 - 
     &   MTop*Sf2245*UStop22

	Sf7216 = 3*MTop*(reMUE*SA + CA*Re(Atop))*UStop21 + 
     &   2*Sf7194*UStop22

	Sf6979 = 6*MTop2*SA*UStop21 + 
     &   SB*Sf1156*Sf86*UStop21 + 
     &   3*MTop*SA*Re(Atop)*UStop22

	tmp19 = CA*reMUE - SA*Re(Atop)

	Sf7213 = 3*MTop*tmp19*UStop21 - 2*Sf7190*UStop22

	Sf6943 = MTop2*UStop21 - MW2*S2B*TB*UStop21 + 
     &   MTop*Sf2196*UStop22

	Sf4251 = 36*MTop2*Sf4247*TB*UStop11 + 
     &   Sf4250*TB*UStop11 - 36*MTop*Sf2245*Sf4247*UStop12

	Sf4242 = 36*MTop2*Sf4238*TB*UStop11 + 
     &   Sf4241*TB*UStop11 - 36*MTop*Sf2245*Sf4238*UStop12

	Sf2325 = MTop2*TB*UStop11 + C2B*MW2*TB*UStop11 - 
     &   MTop*Sf2245*UStop12

	Sf7195 = 3*MTop*(reMUE*SA + CA*Re(Atop))*UStop11 + 
     &   2*Sf7194*UStop12

	Sf4262 = 6*MTop2*SA*UStop11 + 
     &   SB*Sf1156*Sf86*UStop11 + 
     &   3*MTop*SA*Re(Atop)*UStop12

	Sf7191 = 3*MTop*tmp19*UStop11 - 2*Sf7190*UStop12

	Sf4225 = MTop2*UStop11 - MW2*S2B*TB*UStop11 + 
     &   MTop*Sf2196*UStop12

	Sf1759 = Sf1755*Ss443*UStau11*UStau21**2 - 
     &   MZ2*(MZ2*Sf1758 - 
     &      2*SW2*MStau2(3)*UStau12*UStau21*UStau22)

	tmp18 = reMUE*SA - CA*Re(Atau)

	Sf1772 = ML*tmp18*UStau21 - 2*Sf1750*UStau22

	Sf1769 = ML*(CA*reMUE + SA*Re(Atau))*UStau21 + 
     &   2*Sf1744*UStau22

	Sf1751 = ML*tmp18*UStau11 - 2*Sf1750*UStau12

	Sf1745 = ML*(CA*reMUE + SA*Re(Atau))*UStau11 + 
     &   2*Sf1744*UStau12

	Sf4777 = Sf15*Sf4773*USbot11*USbot21**2 - 
     &   MZ2*(MZ2*Sf4776 - 
     &      2*SW2*MSbot2(3)*USbot12*USbot21*USbot22)

	Sf6061 = (Sf5636 + 28*MSbot2(2))*USbot21 - 
     &   12*MBot*Sf2225*USbot22

	Sf6102 = (6*MBot2*MW2**2*Sf4200 - Sf6101)*USbot21 - 
     &   6*MBot*MW2**2*Sf2225*Sf4200*USbot22

	Sf3615 = (Sf3289 + 32*MSbot2(2))*USbot21 - 
     &   12*MBot*Sf2225*USbot22

	Sf3939 = (Sf3289 + 28*MSbot2(2))*USbot21 - 
     &   12*MBot*Sf2225*USbot22

	tmp21 = reMUE*SA - CA*Re(Abot)

	Sf4792 = 3*MBot*tmp21*USbot21 - 2*Sf4768*USbot22

	Sf145 = 6*MBot2*SA*USbot21 - CB*SAB*Sf15*USbot21 + 
     &   3*MBot*SA*Re(Abot)*USbot22

	Sf6100 = (12*MBot2*MW2**3 - Sf6099)*USbot21 - 
     &   12*MBot*MW2**3*Sf2225*USbot22

	Sf6098 = (2*MBot2*MW2**2 - Sf6097)*USbot21 - 
     &   2*MBot*MW2**2*Sf2225*USbot22

	Sf4789 = 3*MBot*(CA*reMUE + SA*Re(Abot))*USbot21 + 
     &   2*Sf4764*USbot22

	Sf5637 = (Sf5636 + 28*MSbot2(1))*USbot11 - 
     &   12*MBot*Sf2225*USbot12

	Sf5713 = (6*MBot2*MW2**2*Sf4120 - Sf5712)*USbot11 - 
     &   6*MBot*MW2**2*Sf2225*Sf4120*USbot12

	Sf3290 = (Sf3289 + 32*MSbot2(1))*USbot11 - 
     &   12*MBot*Sf2225*USbot12

	Sf3815 = (Sf3289 + 28*MSbot2(1))*USbot11 - 
     &   12*MBot*Sf2225*USbot12

	Sf4769 = 3*MBot*tmp21*USbot11 - 2*Sf4768*USbot12

	Sf134 = 6*MBot2*SA*USbot11 - CB*SAB*Sf15*USbot11 + 
     &   3*MBot*SA*Re(Abot)*USbot12

	Sf5711 = (12*MBot2*MW2**3 - Sf5710)*USbot11 - 
     &   12*MBot*MW2**3*Sf2225*USbot12

	Sf5709 = (2*MBot2*MW2**2 - Sf5708)*USbot11 - 
     &   2*MBot*MW2**2*Sf2225*USbot12

	Sf4765 = 3*MBot*(CA*reMUE + SA*Re(Abot))*USbot11 + 
     &   2*Sf4764*USbot12

	Sf9488 = -(Sf9487*MStau2(2)**2) + 
     &   Sf9332*(Sf9332*Sf9333**2*Sf9481*MStau2(1)**2 - 
     &      Sf9333*Sf9482*MStau2(1)*MStau2(2))

	Sf1964 = 4*Sf181*MStau2(2)**2 + 
     &   MSneu2(1)*(25*MSneu2(1) + 3*MStau2(2))

	Sf1902 = 4*Sf181*MStau2(1)**2 + 
     &   MSneu2(1)*(25*MSneu2(1) + 3*MStau2(1))

	Sf1469 = -15*MStau2(2)**2 + 
     &   MSneu2(1)*(5*MSneu2(1) - 2*Sf566*MStau2(2))

	Sf1526 = MSneu2(1)*(5*MSneu2(1) - 10*MStau2(2)) - 
     &   Sf1413*MStau2(2)

	Sf567 = -15*MStau2(1)**2 + 
     &   MSneu2(1)*(5*MSneu2(1) - 2*Sf566*MStau2(1))

	Sf632 = MSneu2(1)*(5*MSneu2(1) - 10*MStau2(1)) - Sf496*MStau2(1)

	Sf1496 = -4*MStau2(2)**2 + 
     &   MSneu2(1)*(MSneu2(1) - Sf174*MStau2(2))

	Sf596 = -4*MStau2(1)**2 + 
     &   MSneu2(1)*(MSneu2(1) - Sf174*MStau2(1))

	Sf707 = 9*MSneu2(1)**2 + Sf706*(MSneu2(1) + 2*MStau2(1))

	Sf1573 = 9*MSneu2(1)**2 + Sf1572*(MSneu2(1) + 2*MStau2(2))

	Sf5393 = 2*Sf2846*MStop2(2)**2 + 
     &   MSbot2(2)*(288*MSbot2(2) + 3*Sf2845*MStop2(2))

	Sf2985 = 2*Sf2846*MStop2(1)**2 + 
     &   MSbot2(2)*(288*MSbot2(2) + 3*Sf2845*MStop2(1))

	Sf6036 = -10*MStop2(2)**2 + 
     &   MSbot2(2)*(MSbot2(2) - 3*Sf175*MStop2(2))

	Sf3695 = -10*MStop2(1)**2 + 
     &   MSbot2(2)*(MSbot2(2) - 3*Sf175*MStop2(1))

	Sf5287 = 2*Sf2846*MStop2(2)**2 + 
     &   MSbot2(1)*(288*MSbot2(1) + 3*Sf2845*MStop2(2))

	Sf2847 = 2*Sf2846*MStop2(1)**2 + 
     &   MSbot2(1)*(288*MSbot2(1) + 3*Sf2845*MStop2(1))

	Sf5881 = Sf5880*MStop2(2) + 
     &   MSbot2(1)*(18*MSbot2(1) + Sf601*MStop2(2))

	Sf5879 = MSbot2(1)*(16*MSbot2(1) - MStop2(2)) - 
     &   3*Sf5878*MStop2(2)

	Sf5616 = -10*MStop2(2)**2 + 
     &   MSbot2(1)*(MSbot2(1) - 3*Sf175*MStop2(2))

	Sf3417 = -10*MStop2(1)**2 + 
     &   MSbot2(1)*(MSbot2(1) - 3*Sf175*MStop2(1))

	Sf10333 = 3*CW*MTop*Sf10268*Ss978 - 
     &   Ss977*UStop12*UStop22 + 
     &   2*Ss979*UStop12*UStop22

	Sf10435 = 3*CW*MTop*Sf10268*Ss966 + 
     &   Ss965*UStop12*UStop22 - 
     &   2*Ss967*UStop12*UStop22

	Sf10055 = CW*ML*Sf9311*Ss937 - 
     &   Ss936*UStau12*UStau22 + 
     &   2*Ss938*UStau12*UStau22

	Sf9312 = CW*ML*Sf9311*Ss924 + 
     &   Ss923*UStau12*UStau22 - 
     &   2*Ss925*UStau12*UStau22

	Sf9108 = 3*CW*MBot*Sf8896*Ss894 - 
     &   Ss893*USbot12*USbot22 + 
     &   2*Ss895*USbot12*USbot22

	Sf8897 = 3*CW*MBot*Sf8896*Ss881 + 
     &   Ss880*USbot12*USbot22 - 
     &   2*Ss882*USbot12*USbot22

	Sf10127 = MZ2*(3*MZ2*Sf9131*Sf9132 - 
     &      Sf10123*MStau2(2)) + 6*MStau2(1)*MStau2(2)*MStau2(4)

	Sf10156 = MZ2*(3*MZ2*Sf9131*Sf9132 - 
     &      Sf10151*MStau2(1)) + 6*MStau2(1)*MStau2(2)*MStau2(4)

	Sf9192 = MZ2*(3*MZ2*Sf9131*Sf9132 - Sf9188*MSbot2(2)) + 
     &   6*MSbot2(1)*MSbot2(2)*MSbot2(4)

	Sf9220 = MZ2*(3*MZ2*Sf9131*Sf9132 - Sf9216*MSbot2(1)) + 
     &   6*MSbot2(1)*MSbot2(2)*MSbot2(4)

	Sf8529 = MStop2(3)**4 + 
     &   MZ2*(MZ2*Sf8528 - 2*MStop2(3)**2*MStop2(4))

	Sf8084 = 1121*MSq2Diff(2,2)**4 + 
     &   MW2*(24*MW2*Sf8081 + 2*Sf8083*MSq2Diff(2,2)**2)

	Sf8082 = 160*MSq2Diff(2,2)**4 + 
     &   MW2*(12*MW2*Sf8081 + Sf8080*MSq2Diff(2,2)**2)

	Sf8005 = 1121*MSq2Diff(1,2)**4 + 
     &   MW2*(24*MW2*Sf8002 + 2*Sf8004*MSq2Diff(1,2)**2)

	Sf8003 = 160*MSq2Diff(1,2)**4 + 
     &   MW2*(12*MW2*Sf8002 + Sf8001*MSq2Diff(1,2)**2)

	Sf7875 = 1121*MSq2Diff(2,1)**4 + 
     &   MW2*(24*MW2*Sf7872 + 2*Sf7874*MSq2Diff(2,1)**2)

	Sf7873 = 160*MSq2Diff(2,1)**4 + 
     &   MW2*(12*MW2*Sf7872 + Sf7871*MSq2Diff(2,1)**2)

	Sf7797 = 1121*MSq2Diff(1,1)**4 + 
     &   MW2*(24*MW2*Sf7794 + 2*Sf7796*MSq2Diff(1,1)**2)

	Sf7795 = 160*MSq2Diff(1,1)**4 + 
     &   MW2*(12*MW2*Sf7794 + Sf7793*MSq2Diff(1,1)**2)

	Sf7710 = MStau2(3)**4 + 
     &   MZ2*(MZ2*Sf7709 - 2*MStau2(3)**2*MStau2(4))

	Sf7590 = 427*MSl2Diff(2,1)**4 + 
     &   MW2*(12*MW2*Sf7589 + 3*Sf7588*MSl2Diff(2,1)**2)

	Sf7532 = 427*MSl2Diff(1,1)**4 + 
     &   MW2*(12*MW2*Sf7531 + 3*Sf7530*MSl2Diff(1,1)**2)

	Sf7518 = MSbot2(3)**4 + 
     &   MZ2*(MZ2*Sf7517 - 2*MSbot2(3)**2*MSbot2(4))

	Sf1537 = MHiggs2(4)*(Sf655*MHiggs2(4) - 180*MStau2(2)) - 
     &   Sf656*MStau2(2)**2

	Sf657 = MHiggs2(4)*(Sf655*MHiggs2(4) - 180*MStau2(1)) - 
     &   Sf656*MStau2(1)**2

	Sf2019 = MZ2*(Sf56*Sf58*Sf60**2 + MZ2*Sf61*Sf65**2) + 
     &   Sf56*Sf58*Sf60**2*MStau2(3)

	Sf6117 = -20*MSbot2(2)**3 + MW2*(MW2*Sf5728 + Sf6116*MSbot2(2))

	Sf5730 = -20*MSbot2(1)**3 + MW2*(MW2*Sf5728 + Sf5729*MSbot2(1))

	Sf7529 = -97*MSl2Diff(1,1)**4 + 
     &   MW2*(12*MW2*Sf559 - 3*Sf7528*MSl2Diff(1,1)**2)

	Sf5574 = 5*MStop2(2)**3 + MW2*(MW2*Sf5573 + 6*MStop2(2)**2)

	Sf5572 = MStop2(2)**3 + MW2*(MW2*Sf5571 - 9*MStop2(2)**2)

	Sf5510 = -5*MStop2(2)**3 + MW2*(MW2*Sf5509 - 6*MStop2(2)**2)

	Sf8079 = -243*MSq2Diff(2,2)**4 + 
     &   MW2*(24*MW2*Sf5141 - 2*Sf8078*MSq2Diff(2,2)**2)

	Sf8077 = -48*MSq2Diff(2,2)**4 + 
     &   MW2*(12*MW2*Sf5141 - Sf8076*MSq2Diff(2,2)**2)

	Sf5088 = 8*MSq2Diff(2,2)**4 + 
     &   MHiggs2(4)*(Sf5087*MHiggs2(4) + Sf4887*MSq2Diff(2,2)**2)

	Sf8000 = -243*MSq2Diff(1,2)**4 + 
     &   MW2*(24*MW2*Sf5022 - 2*Sf7999*MSq2Diff(1,2)**2)

	Sf7998 = -48*MSq2Diff(1,2)**4 + 
     &   MW2*(12*MW2*Sf5022 - Sf7997*MSq2Diff(1,2)**2)

	Sf4962 = 8*MSq2Diff(1,2)**4 + 
     &   MHiggs2(4)*(Sf4961*MHiggs2(4) + Sf4853*MSq2Diff(1,2)**2)

	Sf1788 = MZ2*(MZ2*Sf39*Sf43**2 + Sf56*Sf58*Sf60**2) - 
     &   Sf56*Sf58*Sf60**2*MStau2(3)

	Sf4426 = MHiggs2(4)*
     &    (Sf4292*MHiggs2(4) - 360*MSbot2(2)) - Sf4293*MSbot2(2)**2

	Sf4294 = MHiggs2(4)*
     &    (Sf4292*MHiggs2(4) - 360*MSbot2(1)) - Sf4293*MSbot2(1)**2

	Sf5854 = 40*MStop2(2)**3 + MW2*(MW2*Sf4071 + Sf5853*MStop2(2))

	Sf4074 = 40*MStop2(1)**3 + MW2*(MW2*Sf4071 + Sf4073*MStop2(1))

	Sf3252 = 5*MStop2(1)**3 + MW2*(MW2*Sf3251 + 6*MStop2(1)**2)

	Sf3249 = MStop2(1)**3 + MW2*(MW2*Sf3248 - 9*MStop2(1)**2)

	Sf3131 = -5*MStop2(1)**3 + MW2*(MW2*Sf3130 - 6*MStop2(1)**2)

	Sf7870 = -243*MSq2Diff(2,1)**4 + 
     &   MW2*(24*MW2*Sf2733 - 2*Sf7869*MSq2Diff(2,1)**2)

	Sf7868 = -48*MSq2Diff(2,1)**4 + 
     &   MW2*(12*MW2*Sf2733 - Sf7867*MSq2Diff(2,1)**2)

	Sf2648 = 8*MSq2Diff(2,1)**4 + 
     &   MHiggs2(4)*(Sf2647*MHiggs2(4) + Sf2298*MSq2Diff(2,1)**2)

	Sf7792 = -243*MSq2Diff(1,1)**4 + 
     &   MW2*(24*MW2*Sf2555 - 2*Sf7791*MSq2Diff(1,1)**2)

	Sf7790 = -48*MSq2Diff(1,1)**4 + 
     &   MW2*(12*MW2*Sf2555 - Sf7789*MSq2Diff(1,1)**2)

	Sf2434 = 8*MSq2Diff(1,1)**4 + 
     &   MHiggs2(4)*(Sf2433*MHiggs2(4) + Sf2247*MSq2Diff(1,1)**2)

	Sf2608 = MHiggs2(4)*
     &    (Sf2373*MHiggs2(4) - 120*MSbot2(2)) - Sf2374*MSbot2(2)**2

	Sf2375 = MHiggs2(4)*
     &    (Sf2373*MHiggs2(4) - 120*MSbot2(1)) - Sf2374*MSbot2(1)**2

	Sf1241 = MHiggs2(4)*(7*Sf2*MHiggs2(4) - 60*MStau2(2)) - 
     &   Sf264*MStau2(2)**2

	Sf265 = MHiggs2(4)*(7*Sf2*MHiggs2(4) - 60*MStau2(1)) - 
     &   Sf264*MStau2(1)**2

	Sf516 = MW2*(MW2*Sf174 - Sf515) + MStau2(1)**2

	Sf7587 = -97*MSl2Diff(2,1)**4 + 
     &   MW2*(12*MW2*Sf1462 - 3*Sf7586*MSl2Diff(2,1)**2)

	Sf1427 = MW2*(Sf1426 - MW2*Sf174) - MStau2(2)**2

	Sf10270 = 3*CW*Sf10269 - Ss977*UStop12*UStop22 - 
     &   4*MW*MZ*SAB*SB*SW2*UStop12*UStop22

	Sf1742 = MW2**2*Sf1741 + 
     &   MHiggs2(4)*(-(MW2*Sf1740) + SA*MHiggs2(4))

	Sf10278 = -3*Sf9113*MStop2(2)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf9131*Sf9132*MStop2(2))

	Sf10508 = -3*Sf8911*MStop2(2)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf8929*Sf8930*MStop2(2))

	Sf10271 = -3*Sf9113*MStop2(1)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf9131*Sf9132*MStop2(1))

	Sf10501 = -3*Sf8911*MStop2(1)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf8929*Sf8930*MStop2(1))

	Sf10139 = -3*Sf9113*MStau2(2)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf9131*Sf9132*MStau2(2))

	Sf9386 = -3*Sf8911*MStau2(2)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf8929*Sf8930*MStau2(2))

	Sf10133 = -3*Sf9113*MStau2(1)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf9131*Sf9132*MStau2(1))

	Sf9378 = -3*Sf8911*MStau2(1)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf8929*Sf8930*MStau2(1))

	Sf9204 = -3*Sf9113*MSbot2(2)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf9131*Sf9132*MSbot2(2))

	Sf9009 = -3*Sf8911*MSbot2(2)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf8929*Sf8930*MSbot2(2))

	Sf9198 = -3*Sf9113*MSbot2(1)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf9131*Sf9132*MSbot2(1))

	Sf9001 = -3*Sf8911*MSbot2(1)**2 + 
     &   MZ2*(MZ2**2 - 2*Sf8929*Sf8930*MSbot2(1))

	Sf8809 = MZ2*(12*MZ2 + 43*Sf99) + 8*MSq2Diff(2,2)**2

	Sf8803 = MZ2*(12*MZ2 + 43*Sf91) + 8*MSq2Diff(1,2)**2

	Sf8796 = MZ2*(12*MZ2 + 43*Sf77) + 8*MSq2Diff(2,1)**2

	Sf8790 = MZ2*(12*MZ2 + 43*Sf69) + 8*MSq2Diff(1,1)**2

	Sf8779 = MZ2*(12*MZ2 + 43*Sf47) + 8*MSl2Diff(2,1)**2

	Sf8772 = MZ2*(12*MZ2 + 43*Sf34) + 8*MSl2Diff(1,1)**2

	Sf5419 = MW2**4*(49*MW2**2 + 6*Sf5416) - 9*MTop2*MW2**2*Sf5418

	Sf5313 = MW2**4*(49*MW2**2 + 6*Sf5310) - 9*MTop2*MW2**2*Sf5312

	Sf3019 = MW2**4*(49*MW2**2 + 6*Sf3016) - 9*MTop2*MW2**2*Sf3018

	Sf2892 = MW2**4*(49*MW2**2 + 6*Sf2889) - 9*MTop2*MW2**2*Sf2891

	Sf6900 = MW2**4*(7*MW2**2 + 6*Sf6897) - MTop2*MW2**2*Sf6899

	Sf6649 = MW2**4*(7*MW2**2 + 6*Sf6646) - MTop2*MW2**2*Sf6648

	Sf4704 = MW2**4*(7*MW2**2 + 6*Sf4701) - MTop2*MW2**2*Sf4703

	Sf4607 = MW2**4*(7*MW2**2 + 6*Sf4604) - MTop2*MW2**2*Sf4606

	Sf4519 = -(Sf4518*MHiggs2(4)) + MW2*(2*MW2**2 - Ss25*MSbot2(2))

	Sf4415 = -(Sf4414*MHiggs2(4)) + MW2*(2*MW2**2 - Ss25*MSbot2(1))

	Sf2133 = -(Sf2132*MHiggs2(4)) + MW2*(2*MW2**2 - Ss25*MStau2(2))

	Sf2053 = -(Sf2052*MHiggs2(4)) + MW2*(2*MW2**2 - Ss25*MStau2(1))

	Sf2724 = 3*Sf173*MSbot2(2)**2 + MW2*(MW2**2 - Sf180*MSbot2(2))

	Sf2545 = 3*Sf173*MSbot2(1)**2 + MW2*(MW2**2 - Sf180*MSbot2(1))

	Sf1155 = MW2*(8*MW2 - MZ2*Sf1152) - MZ2**2*Sf1154

	Sf1174 = MW2*(5*MW2 - 2*Sf1173) + 5*MZ2*MSl2Diff(2,1)

	Sf1172 = MW2*(5*MW2 - 2*Sf1171) + 5*MZ2*MSl2Diff(1,1)

	Sf1193 = MW2*(2*MW2 + 4*Sf47) + MSl2Diff(2,1)**2

	Sf164 = MW2*(2*MW2 + 4*Sf34) + MSl2Diff(1,1)**2

	Sf4906 = MW2*(MW2 - 2*Sf99) + MSq2Diff(2,2)**2

	Sf4820 = MW2*(MW2 - 2*Sf91) + MSq2Diff(1,2)**2

	Sf2334 = MW2*(MW2 - 2*Sf77) + MSq2Diff(2,1)**2

	Sf2320 = MW2*(MW2 - 2*Sf69) + MSq2Diff(1,1)**2

	Sf1195 = MW2*(MW2 - 2*Sf47) + MSl2Diff(2,1)**2

	Sf166 = MW2*(MW2 - 2*Sf34) + MSl2Diff(1,1)**2

	Sf988 = ML2*(2*ML2 - 4*MW2*Sf174) - MW2**2*Sf987

	Sf978 = C2B*MW2**2*Sf118 + ML2*(ML2 - MW2*Sf720)

	Sf5198 = -(Sf5171*Sf5197) + 
     &   MHiggs2(4)**2*(-2*Sf5038 + 3*MHiggs2(4))

	Sf5186 = -(Sf5171*Sf5185) + 
     &   MHiggs2(4)**2*(-2*Sf5038 + 3*MHiggs2(4))

	Sf5180 = -(Sf5171*Sf5179) + 
     &   MHiggs2(4)**2*(-2*Sf2573 + 3*MHiggs2(4))

	Sf5174 = -(Sf5171*Sf5173) + 
     &   MHiggs2(4)**2*(-2*Sf2573 + 3*MHiggs2(4))

	Sf2232 = 8*MW2**2 + MHiggs2(4)*(Sf2231 + 9*MHiggs2(4))

	Sf5206 = MHiggs2(4)*(2*Sf5205 + 7*MHiggs2(4)) + 
     &   7*MSq2Diff(2,2)**2

	Sf5194 = MHiggs2(4)*(2*Sf5193 + 7*MHiggs2(4)) + 
     &   7*MSq2Diff(1,2)**2

	Sf3746 = MHiggs2(4)*(2*Sf3745 + 7*MHiggs2(4)) + 
     &   7*MSq2Diff(2,1)**2

	Sf3737 = MHiggs2(4)*(2*Sf3736 + 7*MHiggs2(4)) + 
     &   7*MSq2Diff(1,1)**2

	Sf5207 = MHiggs2(4)*(-2*Sf4868 + 3*MHiggs2(4)) + 
     &   3*MSq2Diff(2,2)**2

	Sf5195 = MHiggs2(4)*(-2*Sf4839 + 3*MHiggs2(4)) + 
     &   3*MSq2Diff(1,2)**2

	Sf3747 = MHiggs2(4)*(-2*Sf2271 + 3*MHiggs2(4)) + 
     &   3*MSq2Diff(2,1)**2

	Sf3738 = MHiggs2(4)*(-2*Sf2207 + 3*MHiggs2(4)) + 
     &   3*MSq2Diff(1,1)**2

	Sf157 = 2*MW2**2 + MHiggs2(4)*(-Sf156 + 3*MHiggs2(4))

	Sf1187 = 2*MW2**2 + MHiggs2(4)*(-Sf1186 + 3*MHiggs2(4))

	Sf1240 = Sf262*MStau2(2)**2 + 
     &   MHiggs2(4)*(2*MHiggs2(4) + Sf261*MStau2(2))

	Sf263 = Sf262*MStau2(1)**2 + 
     &   MHiggs2(4)*(2*MHiggs2(4) + Sf261*MStau2(1))

	Sf4867 = -2*Sf99*MHiggs2(4) + MHiggs2(4)**2 + MSq2Diff(2,2)**2

	Sf4831 = -2*Sf91*MHiggs2(4) + MHiggs2(4)**2 + MSq2Diff(1,2)**2

	Sf2270 = -2*Sf77*MHiggs2(4) + MHiggs2(4)**2 + MSq2Diff(2,1)**2

	Sf2206 = -2*Sf69*MHiggs2(4) + MHiggs2(4)**2 + MSq2Diff(1,1)**2

	Sf1184 = -2*Sf47*MHiggs2(4) + MHiggs2(4)**2 + MSl2Diff(2,1)**2

	Sf154 = -2*Sf34*MHiggs2(4) + MHiggs2(4)**2 + MSl2Diff(1,1)**2

	Sf9375 = -(Sf8923*Sf9374) + 
     &   MHiggs2(2)**2*(-2*Sf9342 + 3*MHiggs2(2))

	Sf8998 = -(Sf8923*Sf8997) + 
     &   MHiggs2(2)**2*(-2*Sf8944 + 3*MHiggs2(2))

	Sf10498 = -(Sf10497*Sf8923) + 
     &   MHiggs2(2)**2*(-2*Sf10274 + 3*MHiggs2(2))

	Sf10482 = MHiggs2(2)*(-2*MZ2 + 4*MHiggs2(2)) + Sf8925*MStop2(2)

	Sf10448 = MHiggs2(2)*(-2*MZ2 + 4*MHiggs2(2)) + Sf8925*MStop2(1)

	Sf9366 = MHiggs2(2)*(-2*MZ2 + 4*MHiggs2(2)) + Sf8925*MStau2(2)

	Sf9335 = MHiggs2(2)*(-2*MZ2 + 4*MHiggs2(2)) + Sf8925*MStau2(1)

	Sf8974 = MHiggs2(2)*(-2*MZ2 + 4*MHiggs2(2)) + Sf8925*MSbot2(2)

	Sf8926 = MHiggs2(2)*(-2*MZ2 + 4*MHiggs2(2)) + Sf8925*MSbot2(1)

	Sf10130 = -(Sf10129*Sf9125) + 
     &   MHiggs2(1)**2*(-2*Sf9342 + 3*MHiggs2(1))

	Sf9195 = -(Sf9125*Sf9194) + 
     &   MHiggs2(1)**2*(-2*Sf8944 + 3*MHiggs2(1))

	Sf10389 = -(Sf10388*Sf9125) + 
     &   MHiggs2(1)**2*(-2*Sf10274 + 3*MHiggs2(1))

	Sf10380 = MHiggs2(1)*(-2*MZ2 + 4*MHiggs2(1)) + Sf9127*MStop2(2)

	Sf10353 = MHiggs2(1)*(-2*MZ2 + 4*MHiggs2(1)) + Sf9127*MStop2(1)

	Sf10102 = MHiggs2(1)*(-2*MZ2 + 4*MHiggs2(1)) + Sf9127*MStau2(2)

	Sf10072 = MHiggs2(1)*(-2*MZ2 + 4*MHiggs2(1)) + Sf9127*MStau2(1)

	Sf9164 = MHiggs2(1)*(-2*MZ2 + 4*MHiggs2(1)) + Sf9127*MSbot2(2)

	Sf9128 = MHiggs2(1)*(-2*MZ2 + 4*MHiggs2(1)) + Sf9127*MSbot2(1)

	Sf10243 = Ss965*UStop21**2 - 2*Sf10242*UStop22

	Sf7271 = Sf86**2*UStop21**2 + 16*MZ2**2*SW2**2*UStop22**2

	Sf117 = -4*MZ2*SW2*UStop22**2 + Sf86*(1 - UStop22**2)

	Sf5056 = MBot*Sf2225*TB*USbot22*UStop21 - 
     &   MTop*Sf2245*USbot21*UStop22

	Sf4930 = MBot*Sf2225*TB*USbot12*UStop21 - 
     &   MTop*Sf2245*USbot11*UStop22

	Sf5143 = 3*Sf99*TB*MSq2Diff(2,2)*UStop21 + 
     &   4*MTop*Sf2245*Sf5140*UStop22

	Sf5024 = 3*Sf91*TB*MSq2Diff(1,2)*UStop21 + 
     &   4*MTop*Sf2245*Sf5021*UStop22

	Sf7218 = Sf7217*UStop21 + Sf7216*UStop22

	Sf7215 = Sf7214*UStop21 - Sf7213*UStop22

	Sf4816 = Sf102*MSq2Diff(2,2)*UStop21**2 + 
     &   Sf80*MSq2Diff(2,1)*UStop22**2

	Sf4811 = Sf94*MSq2Diff(1,2)*UStop21**2 + 
     &   Sf72*MSq2Diff(1,1)*UStop22**2

	Sf10223 = -2*Sf10222*UStop12 + Ss965*UStop22**2

	Sf7253 = 16*MZ2**2*SW2**2*UStop12**2 + Sf86**2*UStop22**2

	Sf87 = -4*MZ2*SW2*UStop12**2 + Sf86*UStop22**2

	Sf8516 = Sf8512*UStop22**2 + Sf8515*(1 - UStop22**2)

	Sf8500 = Sf8496*UStop22**2 + Sf8499*(1 - UStop22**2)

	Sf4817 = Sf80*UStop22**2 + Sf102*(1 - UStop22**2)

	Sf4812 = Sf72*UStop22**2 + Sf94*(1 - UStop22**2)

	Sf2602 = MBot*Sf2225*TB*USbot22*UStop11 - 
     &   MTop*Sf2245*USbot21*UStop12

	Sf2364 = MBot*Sf2225*TB*USbot12*UStop11 - 
     &   MTop*Sf2245*USbot11*UStop12

	Sf2735 = 3*Sf77*TB*MSq2Diff(2,1)*UStop11 + 
     &   4*MTop*Sf2245*Sf2732*UStop12

	Sf2557 = 3*Sf69*TB*MSq2Diff(1,1)*UStop11 + 
     &   4*MTop*Sf2245*Sf2554*UStop12

	Sf7197 = Sf7196*UStop11 + Sf7195*UStop12

	Sf7193 = Sf7192*UStop11 - Sf7191*UStop12

	Sf8601 = 3*Sf8600*MStau2(1) + 
     &   7*Sf8538*UStau12*UStau21*UStau22

	Sf8577 = Sf8576*MStau2(1) + 
     &   Sf8538*UStau12*UStau21*UStau22

	Sf8540 = 3*Sf8539*MStau2(1) + 
     &   Sf8538*UStau12*UStau21*UStau22

	Sf2171 = Sf345*UStau21**2 + ML2*Sf215**2*UStau22**2

	Sf1182 = MW2*Sf149*UStau21**2 + 2*Sf150*UStau22**2

	Sf1183 = MW2*S2B*UStau21**2 - 2*Sf152*UStau22**2

	Sf1636 = Sf215**2*UStau22**2 + ML2*(1 - UStau22**2)

	Sf1362 = 3*MStau2(2)*UStau21 + ML*Sf215*UStau22

	Sf1190 = ML*TB*UStau21 + Sf160*UStau22

	Sf1227 = Sf241*UStau21 - ML*Sf215*UStau22

	Sf1226 = Sf239*UStau21 + 2*ML*Sf215*UStau22

	Sf1774 = Sf1773*UStau21 - Sf1772*UStau22

	Sf1771 = Sf1770*UStau21 + Sf1769*UStau22

	Sf1882 = Sf1140*UStau21 + 6*ML*Sf215*UStau22

	Sf1881 = Sf1138*UStau21 + 3*ML*Sf215*UStau22

	Sf1879 = Sf1134*UStau21 - 3*ML*Sf215*UStau22

	Sf1877 = Sf1130*UStau21 + 3*ML*Sf215*UStau22

	Sf1876 = Sf1127*UStau21 - 3*ML*Sf215*UStau22

	Sf1866 = Sf1107*UStau21 + 3*ML*Sf215*UStau22

	Sf1297 = 2*ML**3*UStau21 - ML2*Sf215*UStau22

	Sf1293 = ML*UStau21 - Sf215*UStau22

	Sf8602 = 2*Sf8534*Sf8598*UStau12**2 + 
     &   Sf8536*Sf8597*Ss443**2*UStau21**2

	Sf8541 = 2*Sf8534*Sf8535*UStau12**2 + 
     &   Sf8536*Ss443**2*UStau21**2

	Sf1177 = Sf50*MSl2Diff(2,1)*UStau21**2 + 
     &   Sf37*MSl2Diff(1,1)*UStau22**2

	Sf8819 = 2*Sf8536*UStau21**2 + Sf8534*UStau22**2

	Sf7699 = (Sf7667*UStau21**2 + 
     &     (Sf1195**3*Sf7655*UStau22**2)/Sf166**3)/Sf1195**3

	Sf1178 = Sf37*UStau22**2 + Sf50*(1 - UStau22**2)

	Sf2102 = ML2*Sf215**2*UStau12**2 + Sf345*UStau22**2

	Sf151 = 2*Sf150*UStau12**2 + MW2*Sf149*UStau22**2

	Sf1175 = Sf1172*UStau22**2 + Sf1174*(1 - UStau22**2)

	Sf153 = -2*Sf152*UStau12**2 + MW2*S2B*UStau22**2

	Sf811 = Sf215**2*UStau12**2 + ML2*UStau22**2

	Sf444 = 3*MStau2(1)*UStau11 + ML*Sf215*UStau12

	Sf161 = ML*TB*UStau11 + Sf160*UStau12

	Sf242 = Sf241*UStau11 - ML*Sf215*UStau12

	Sf240 = Sf239*UStau11 + 2*ML*Sf215*UStau12

	Sf1753 = Sf1752*UStau11 - Sf1751*UStau12

	Sf1747 = Sf1746*UStau11 + Sf1745*UStau12

	Sf1141 = Sf1140*UStau11 + 6*ML*Sf215*UStau12

	Sf1139 = Sf1138*UStau11 + 3*ML*Sf215*UStau12

	Sf1135 = Sf1134*UStau11 - 3*ML*Sf215*UStau12

	Sf1131 = Sf1130*UStau11 + 3*ML*Sf215*UStau12

	Sf1128 = Sf1127*UStau11 - 3*ML*Sf215*UStau12

	Sf1108 = Sf1107*UStau11 + 3*ML*Sf215*UStau12

	Sf346 = 2*ML**3*UStau11 - ML2*Sf215*UStau12

	Sf341 = ML*UStau11 - Sf215*UStau12

	Sf3578 = Sf2807*USbot21**2 + MBot2*Sf2225**2*USbot22**2

	Sf7241 = Sf15**2*USbot21**2 + 4*MZ2**2*SW2**2*USbot22**2

	Sf30 = -2*MZ2*SW2*USbot22**2 + Sf15*(1 - USbot22**2)

	end

