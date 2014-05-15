* gm2_2L_full1.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 10 Aug 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full1
	implicit none

#include "gm2_2L_full.h"

	Sf7402 = USbot22**2

	Sf7314 = UStau22**2

	Sf7368 = UStop22**2

	Sf5 = S2B

	Sf8599 = UStau11*
     &   (42*MZ2 + 6*MSneu2(1) - 177*MStau2(1) - 111*MStau2(2) - 
     &     37*(MZ2 - 3*MStau2(1) - 3*MStau2(2))*UStau22**2)

	Sf8537 = UStau11*
     &   (-9*MZ2 + 70*MStau2(1) + 42*MStau2(2) + 
     &     14*(MZ2 - 3*MStau2(1) - 3*MStau2(2))*UStau22**2)

	Sf7780 = UStau11*
     &   (-8*MSl2Diff(1,1) + MStau2(3)*(26 - 26*UStau22**2))

	Sf8629 = 15*(4 - USbot22**2)

	Sf8821 = -300 + 372*USbot22**2

	Sf8838 = 7*(-7 + 9*USbot22**2)

	Sf7455 = 5*(-3 + 4*UStau22**2)

	Sf7406 = 3*UStop22**2

	Sf7325 = 3*UStau22**2

	Sf7359 = 3*USbot22**2

	Sf7379 = 3*(1 + USbot22**2)

	Sf8754 = UStop11*(478 - 52*UStop22**2)

	Sf8751 = UStop11*(235 - 97*UStop22**2)

	Sf8746 = UStop11*(43 - 13*UStop22**2)

	Sf8744 = UStop11*(49 - 19*UStop22**2)

	Sf7320 = UStau12*UStau22

	Sf7771 = UStau11*(205 - 71*UStau22**2)

	Sf7773 = UStau11*(72 - 23*UStau22**2)

	Sf7775 = UStau11*(25 - 7*UStau22**2)

	Sf8600 = UStau11*(13 - 7*UStau22**2)

	Sf8607 = UStau11*(12 - 7*UStau22**2)

	Sf8603 = UStau11*(19 - 14*UStau22**2)

	Sf8542 = UStau11*(-3 + 7*UStau22**2)

	Sf8576 = UStau11*(5 - 3*UStau22**2)

	Sf8543 = UStau11*(3 - 2*UStau22**2)

	Sf8578 = UStau11*(5 - 4*UStau22**2)

	Sf8606 = UStau11*(36 - 35*UStau22**2)

	Sf8836 = -25*(-1 + USbot22**2)

	Sf1211 = MSneu2(1)*(MSneu2(1)**3 - 4*MStau2(2)**3)

	Sf218 = MSneu2(1)*(MSneu2(1)**3 - 4*MStau2(1)**3)

	Sf9052 = MZ2*MHiggs2(2)*(-5*MZ2**2 + 5*MHiggs2(2)**2)

	Sf8977 = MZ2*MHiggs2(2)*(-MZ2**2 + MHiggs2(2)**2)

	Sf9182 = MZ2*MHiggs2(1)*(-5*MZ2**2 + 5*MHiggs2(1)**2)

	Sf9167 = MZ2*MHiggs2(1)*(-MZ2**2 + MHiggs2(1)**2)

	Sf8909 = MZ2*(1 + 2*Ss474)

	Sf9017 = (1 + 2*Ss269)*MHiggs2(2)

	Sf9222 = (1 + 2*Ss263)*MHiggs2(1)

	Sf8755 = 4*UStop11*(31 - 13*UStop22**2)

	Sf7393 = 3*UStop12*UStop22

	Sf5225 = MStop2(3)*UStop21*UStop22

	Sf8748 = -(UStop11*(-4 + UStop22**2))

	Sf7779 = 2*UStau11*(54 - 13*UStau22**2)

	Sf7777 = -(UStau11*(-4 + UStau22**2))

	Sf8539 = -(UStau11*(-2 + UStau22**2))

	Sf7301 = -3*USbot12*USbot22

	Sf5224 = MSbot2(3)*USbot21*USbot22

	Sf7394 = -(MZ2*Sf7393)

	Sf7321 = -(MZ2*Sf7320)

	Sf8520 = -3*UStop12*UStop22

	Sf7701 = -(Sf7320*(Ss922 + 2*SW2))

	Sf60 = -(Sf7320*(Ss443 + 2*MZ2*SW2))

	Sf7509 = -3*USbot12*USbot22

	Of1 = 27*MBot2 - 27*MTop2 + MW2 - 27*C2B*MW2

	Of2 = 18*MBot2 - 18*MTop2 - MW2 - 18*C2B*MW2

	Of3 = 3*MBot2 - 3*MTop2 + MW2 - 3*C2B*MW2

	Sf7286 = 4 + 48*SW2

	Sf10022 = 296*MStau2(2)**6 + 
     &   MZ2*(-652*MStau2(2)**5 + 
     &      MZ2*(595*MStau2(2)**4 + 
     &         MZ2*(-268*MStau2(2)**3 + 
     &            MZ2*(-38*MStau2(2)**2 + 
     &               MZ2*(19*MZ2 + 48*MStau2(2))))))

	Sf9917 = 296*MStau2(1)**6 + 
     &   MZ2*(-652*MStau2(1)**5 + 
     &      MZ2*(595*MStau2(1)**4 + 
     &         MZ2*(-268*MStau2(1)**3 + 
     &            MZ2*(-38*MStau2(1)**2 + 
     &               MZ2*(19*MZ2 + 48*MStau2(1))))))

	Sf10021 = 1428*MStau2(2)**6 + 
     &   MZ2*(-3198*MStau2(2)**5 + 
     &      MZ2*(1413*MStau2(2)**4 + 
     &         MZ2*(-2942*MStau2(2)**3 + 
     &            MZ2*(MZ2*(5*MZ2 - 332*MStau2(2)) + 
     &               2562*MStau2(2)**2))))

	Sf9918 = 1428*MStau2(1)**6 + 
     &   MZ2*(-3198*MStau2(1)**5 + 
     &      MZ2*(1413*MStau2(1)**4 + 
     &         MZ2*(-2942*MStau2(1)**3 + 
     &            MZ2*(MZ2*(5*MZ2 - 332*MStau2(1)) + 
     &               2562*MStau2(1)**2))))

	Sf10050 = -56*MStau2(2)**6 + 
     &   MHiggs2(2)*(232*MStau2(2)**5 + 
     &      MHiggs2(2)*(-527*MStau2(2)**4 + 
     &         MHiggs2(2)*(516*MStau2(2)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(25*MHiggs2(2) - 132*MStau2(2)) - 
     &               58*MStau2(2)**2))))

	Sf9965 = -56*MStau2(1)**6 + 
     &   MHiggs2(2)*(232*MStau2(1)**5 + 
     &      MHiggs2(2)*(-527*MStau2(1)**4 + 
     &         MHiggs2(2)*(516*MStau2(1)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(25*MHiggs2(2) - 132*MStau2(1)) - 
     &               58*MStau2(1)**2))))

	Sf10014 = -296*MStau2(2)**6 + 
     &   MHiggs2(2)*(820*MStau2(2)**5 + 
     &      MHiggs2(2)*(-1145*MStau2(2)**4 + 
     &         MHiggs2(2)*(1164*MStau2(2)**3 + 
     &            MHiggs2(2)*
     &             (-610*MStau2(2)**2 + 
     &               MHiggs2(2)*(19*MHiggs2(2) + 48*MStau2(2))))))

	Sf9907 = -296*MStau2(1)**6 + 
     &   MHiggs2(2)*(820*MStau2(1)**5 + 
     &      MHiggs2(2)*(-1145*MStau2(1)**4 + 
     &         MHiggs2(2)*(1164*MStau2(1)**3 + 
     &            MHiggs2(2)*
     &             (-610*MStau2(1)**2 + 
     &               MHiggs2(2)*(19*MHiggs2(2) + 48*MStau2(1))))))

	Sf10013 = -780*MStau2(2)**6 + 
     &   MHiggs2(2)*(-1526*MStau2(2)**5 + 
     &      MHiggs2(2)*(8305*MStau2(2)**4 + 
     &         MHiggs2(2)*(-8574*MStau2(2)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(5*MHiggs2(2) - 332*MStau2(2)) + 
     &               2830*MStau2(2)**2))))

	Sf9908 = -780*MStau2(1)**6 + 
     &   MHiggs2(2)*(-1526*MStau2(1)**5 + 
     &      MHiggs2(2)*(8305*MStau2(1)**4 + 
     &         MHiggs2(2)*(-8574*MStau2(1)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(5*MHiggs2(2) - 332*MStau2(1)) + 
     &               2830*MStau2(1)**2))))

	Sf10000 = 36*MStau2(2)**6 + 
     &   MHiggs2(2)*(8*MStau2(2)**5 + 
     &      MHiggs2(2)*(5*MStau2(2)**4 + 
     &         MHiggs2(2)*(-78*MStau2(2)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(MHiggs2(2) - 16*MStau2(2)) + 
     &               44*MStau2(2)**2))))

	Sf10049 = 420*MStau2(2)**6 + 
     &   MHiggs2(2)*(-1234*MStau2(2)**5 + 
     &      MHiggs2(2)*(1013*MStau2(2)**4 + 
     &         MHiggs2(2)*(-1110*MStau2(2)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(MHiggs2(2) - 268*MStau2(2)) + 
     &               1250*MStau2(2)**2))))

	Sf9892 = 36*MStau2(1)**6 + 
     &   MHiggs2(2)*(8*MStau2(1)**5 + 
     &      MHiggs2(2)*(5*MStau2(1)**4 + 
     &         MHiggs2(2)*(-78*MStau2(1)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(MHiggs2(2) - 16*MStau2(1)) + 
     &               44*MStau2(1)**2))))

	Sf9966 = 420*MStau2(1)**6 + 
     &   MHiggs2(2)*(-1234*MStau2(1)**5 + 
     &      MHiggs2(2)*(1013*MStau2(1)**4 + 
     &         MHiggs2(2)*(-1110*MStau2(1)**3 + 
     &            MHiggs2(2)*
     &             (MHiggs2(2)*(MHiggs2(2) - 268*MStau2(1)) + 
     &               1250*MStau2(1)**2))))

	Sf9809 = -104*MStau2(2)**5 + 
     &   MZ2*(498*MStau2(2)**4 + 
     &      MZ2*(-1494*MStau2(2)**3 + 
     &         MZ2*(MZ2*(163*MZ2 - 944*MStau2(2)) + 
     &            1881*MStau2(2)**2)))

	Sf9669 = -104*MStau2(1)**5 + 
     &   MZ2*(498*MStau2(1)**4 + 
     &      MZ2*(-1494*MStau2(1)**3 + 
     &         MZ2*(MZ2*(163*MZ2 - 944*MStau2(1)) + 
     &            1881*MStau2(1)**2)))

	Sf9808 = 288*MStau2(2)**5 + 
     &   MZ2*(-918*MStau2(2)**4 + 
     &      MZ2*(-819*MStau2(2)**3 + 
     &         MZ2*(MZ2*(63*MZ2 - 849*MStau2(2)) + 
     &            1675*MStau2(2)**2)))

	Sf9670 = 288*MStau2(1)**5 + 
     &   MZ2*(-918*MStau2(1)**4 + 
     &      MZ2*(-819*MStau2(1)**3 + 
     &         MZ2*(MZ2*(63*MZ2 - 849*MStau2(1)) + 
     &            1675*MStau2(1)**2)))

	Sf10020 = 2496*MStau2(2)**5 + 
     &   MZ2*(-2268*MStau2(2)**4 + 
     &      MZ2*(-2110*MStau2(2)**3 + 
     &         MZ2*(MZ2*(36*MZ2 - 1826*MStau2(2)) + 
     &            6908*MStau2(2)**2)))

	Sf9919 = 2496*MStau2(1)**5 + 
     &   MZ2*(-2268*MStau2(1)**4 + 
     &      MZ2*(-2110*MStau2(1)**3 + 
     &         MZ2*(MZ2*(36*MZ2 - 1826*MStau2(1)) + 
     &            6908*MStau2(1)**2)))

	Sf9830 = -8*MStau2(2)**5 + 
     &   MZ2*(66*MStau2(2)**4 + 
     &      MZ2*(-270*MStau2(2)**3 + 
     &         MZ2*(MZ2*(25*MZ2 - 164*MStau2(2)) + 
     &            351*MStau2(2)**2)))

	Sf9639 = -8*MStau2(1)**5 + 
     &   MZ2*(66*MStau2(1)**4 + 
     &      MZ2*(-270*MStau2(1)**3 + 
     &         MZ2*(MZ2*(25*MZ2 - 164*MStau2(1)) + 
     &            351*MStau2(1)**2)))

	Sf8723 = -99*MZ2**5 + 
     &   MW2*(211*MZ2**4 + 
     &      MW2*(-208*MZ2**3 + 
     &         MW2*(MW2*(256*MW2 - 640*MZ2) + 480*MZ2**2)))

	Sf9802 = -232*MStau2(2)**5 + 
     &   MHiggs2(2)*(550*MStau2(2)**4 + 
     &      MHiggs2(2)*(-1496*MStau2(2)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (163*MHiggs2(2) - 1034*MStau2(2)) + 
     &            2049*MStau2(2)**2)))

	Sf9656 = -232*MStau2(1)**5 + 
     &   MHiggs2(2)*(550*MStau2(1)**4 + 
     &      MHiggs2(2)*(-1496*MStau2(1)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (163*MHiggs2(2) - 1034*MStau2(1)) + 
     &            2049*MStau2(1)**2)))

	Sf9801 = 794*MStau2(2)**5 + 
     &   MHiggs2(2)*(-84*MStau2(2)**4 + 
     &      MHiggs2(2)*(-2832*MStau2(2)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (63*MHiggs2(2) - 1070*MStau2(2)) + 
     &            3133*MStau2(2)**2)))

	Sf9657 = 794*MStau2(1)**5 + 
     &   MHiggs2(2)*(-84*MStau2(1)**4 + 
     &      MHiggs2(2)*(-2832*MStau2(1)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (63*MHiggs2(2) - 1070*MStau2(1)) + 
     &            3133*MStau2(1)**2)))

	Sf9823 = -40*MStau2(2)**5 + 
     &   MHiggs2(2)*(70*MStau2(2)**4 + 
     &      MHiggs2(2)*(-248*MStau2(2)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (25*MHiggs2(2) - 182*MStau2(2)) + 
     &            375*MStau2(2)**2)))

	Sf9626 = -40*MStau2(1)**5 + 
     &   MHiggs2(2)*(70*MStau2(1)**4 + 
     &      MHiggs2(2)*(-248*MStau2(1)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (25*MHiggs2(2) - 182*MStau2(1)) + 
     &            375*MStau2(1)**2)))

	Sf10012 = -208*MStau2(2)**5 + 
     &   MHiggs2(2)*(-1306*MStau2(2)**4 + 
     &      MHiggs2(2)*(-3089*MStau2(2)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (18*MHiggs2(2) - 977*MStau2(2)) + 
     &            3258*MStau2(2)**2)))

	Sf9909 = -208*MStau2(1)**5 + 
     &   MHiggs2(2)*(-1306*MStau2(1)**4 + 
     &      MHiggs2(2)*(-3089*MStau2(1)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (18*MHiggs2(2) - 977*MStau2(1)) + 
     &            3258*MStau2(1)**2)))

	Sf9822 = 158*MStau2(2)**5 + 
     &   MHiggs2(2)*(60*MStau2(2)**4 + 
     &      MHiggs2(2)*(-528*MStau2(2)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (9*MHiggs2(2) - 146*MStau2(2)) + 
     &            475*MStau2(2)**2)))

	Sf9627 = 158*MStau2(1)**5 + 
     &   MHiggs2(2)*(60*MStau2(1)**4 + 
     &      MHiggs2(2)*(-528*MStau2(1)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (9*MHiggs2(2) - 146*MStau2(1)) + 
     &            475*MStau2(1)**2)))

	Sf9999 = 40*MStau2(2)**5 + 
     &   MHiggs2(2)*(13*MStau2(2)**4 + 
     &      MHiggs2(2)*(-61*MStau2(2)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (3*MHiggs2(2) - 34*MStau2(2)) + 72*MStau2(2)**2)
     &         ))

	Sf9893 = 40*MStau2(1)**5 + 
     &   MHiggs2(2)*(13*MStau2(1)**4 + 
     &      MHiggs2(2)*(-61*MStau2(1)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*
     &             (3*MHiggs2(2) - 34*MStau2(1)) + 72*MStau2(1)**2)
     &         ))

	Sf9998 = 8*MStau2(2)**5 + 
     &   MHiggs2(2)*(10*MStau2(2)**4 + 
     &      MHiggs2(2)*(-45*MStau2(2)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*(MHiggs2(2) - 5*MStau2(2)) + 
     &            39*MStau2(2)**2)))

	Sf9891 = 8*MStau2(1)**5 + 
     &   MHiggs2(2)*(10*MStau2(1)**4 + 
     &      MHiggs2(2)*(-45*MStau2(1)**3 + 
     &         MHiggs2(2)*(MHiggs2(2)*(MHiggs2(2) - 5*MStau2(1)) + 
     &            39*MStau2(1)**2)))

	Sf7479 = 1 + 3*UStop22**2

	Sf7449 = 1 + UStau22**2

	Sf7597 = 21613*MStau2(2)**4 + 
     &   MSneu2(1)*(-4382*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (13621*MSneu2(1) - 3518*MStau2(2)) - 
     &         4662*MStau2(2)**2))

	Sf7539 = 21613*MStau2(1)**4 + 
     &   MSneu2(1)*(-4382*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (13621*MSneu2(1) - 3518*MStau2(1)) - 
     &         4662*MStau2(1)**2))

	Sf7627 = 14539*MStau2(2)**4 + 
     &   MSneu2(1)*(10*MStau2(2)**3 + 
     &      MSneu2(1)*(870*MStau2(2)**2 + 
     &         MSneu2(1)*(11515*MSneu2(1) + 442*MStau2(2))))

	Sf7569 = 14539*MStau2(1)**4 + 
     &   MSneu2(1)*(10*MStau2(1)**3 + 
     &      MSneu2(1)*(870*MStau2(1)**2 + 
     &         MSneu2(1)*(11515*MSneu2(1) + 442*MStau2(1))))

	Sf7697 = -2060*MStau2(2)**4 + 
     &   MSneu2(1)*(863*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (2381*MSneu2(1) - 1019*MStau2(2)) - 
     &         165*MStau2(2)**2))

	Sf7683 = -2060*MStau2(1)**4 + 
     &   MSneu2(1)*(863*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (2381*MSneu2(1) - 1019*MStau2(1)) - 
     &         165*MStau2(1)**2))

	Sf1705 = 2350*MStau2(2)**4 + 
     &   MSneu2(1)*(718*MStau2(2)**3 + 
     &      MSneu2(1)*(792*MStau2(2)**2 + 
     &         MSneu2(1)*(2350*MSneu2(1) + 718*MStau2(2))))

	Sf906 = 2350*MStau2(1)**4 + 
     &   MSneu2(1)*(718*MStau2(1)**3 + 
     &      MSneu2(1)*(792*MStau2(1)**2 + 
     &         MSneu2(1)*(2350*MSneu2(1) + 718*MStau2(1))))

	Sf7605 = 3379*MStau2(2)**4 + 
     &   MSneu2(1)*(-1754*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (1867*MSneu2(1) - 1466*MStau2(2)) - 
     &         954*MStau2(2)**2))

	Sf7547 = 3379*MStau2(1)**4 + 
     &   MSneu2(1)*(-1754*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (1867*MSneu2(1) - 1466*MStau2(1)) - 
     &         954*MStau2(1)**2))

	Sf7744 = -3209*MStau2(2)**4 + 
     &   MSneu2(1)*(1919*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (1484*MSneu2(1) - 467*MStau2(2)) + 273*MStau2(2)**2
     &         ))

	Sf7720 = -3209*MStau2(1)**4 + 
     &   MSneu2(1)*(1919*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*
     &          (1484*MSneu2(1) - 467*MStau2(1)) + 273*MStau2(1)**2
     &         ))

	Sf7612 = 75*MStau2(2)**4 + 
     &   MSneu2(1)*(1212*MStau2(2)**3 + 
     &      MSneu2(1)*(-368*MStau2(2)**2 + 
     &         MSneu2(1)*(489*MSneu2(1) + 384*MStau2(2))))

	Sf7554 = 75*MStau2(1)**4 + 
     &   MSneu2(1)*(1212*MStau2(1)**3 + 
     &      MSneu2(1)*(-368*MStau2(1)**2 + 
     &         MSneu2(1)*(489*MSneu2(1) + 384*MStau2(1))))

	Sf7666 = -187*MStau2(2)**4 + 
     &   MSneu2(1)*(145*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(280*MSneu2(1) - 253*MStau2(2)) + 
     &         15*MStau2(2)**2))

	Sf7654 = -187*MStau2(1)**4 + 
     &   MSneu2(1)*(145*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(280*MSneu2(1) - 253*MStau2(1)) + 
     &         15*MStau2(1)**2))

	Sf1394 = 275*MStau2(2)**4 + 
     &   MSneu2(1)*(53*MStau2(2)**3 + 
     &      MSneu2(1)*(72*MStau2(2)**2 + 
     &         MSneu2(1)*(275*MSneu2(1) + 53*MStau2(2))))

	Sf477 = 275*MStau2(1)**4 + 
     &   MSneu2(1)*(53*MStau2(1)**3 + 
     &      MSneu2(1)*(72*MStau2(1)**2 + 
     &         MSneu2(1)*(275*MSneu2(1) + 53*MStau2(1))))

	Sf1716 = -199*MStau2(2)**4 + 
     &   MSneu2(1)*(-146*MStau2(2)**3 + 
     &      MSneu2(1)*(54*MStau2(2)**2 + 
     &         MSneu2(1)*(253*MSneu2(1) + 38*MStau2(2))))

	Sf926 = -199*MStau2(1)**4 + 
     &   MSneu2(1)*(-146*MStau2(1)**3 + 
     &      MSneu2(1)*(54*MStau2(1)**2 + 
     &         MSneu2(1)*(253*MSneu2(1) + 38*MStau2(1))))

	Sf1402 = -211*MStau2(2)**4 + 
     &   MSneu2(1)*(46*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(193*MSneu2(1) - 58*MStau2(2)) + 
     &         30*MStau2(2)**2))

	Sf485 = -211*MStau2(1)**4 + 
     &   MSneu2(1)*(46*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(193*MSneu2(1) - 58*MStau2(1)) + 
     &         30*MStau2(1)**2))

	Sf7635 = 275*MStau2(2)**4 + 
     &   MSneu2(1)*(-214*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(155*MSneu2(1) - 118*MStau2(2)) - 
     &         18*MStau2(2)**2))

	Sf7577 = 275*MStau2(1)**4 + 
     &   MSneu2(1)*(-214*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(155*MSneu2(1) - 118*MStau2(1)) - 
     &         18*MStau2(1)**2))

	Sf7619 = 200*MStau2(2)**4 + 
     &   MSneu2(1)*(77*MStau2(2)**3 + 
     &      MSneu2(1)*(42*MStau2(2)**2 + 
     &         MSneu2(1)*(140*MSneu2(1) + 77*MStau2(2))))

	Sf7561 = 200*MStau2(1)**4 + 
     &   MSneu2(1)*(77*MStau2(1)**3 + 
     &      MSneu2(1)*(42*MStau2(1)**2 + 
     &         MSneu2(1)*(140*MSneu2(1) + 77*MStau2(1))))

	Sf8331 = -400*MStau2(2)**4 + 
     &   MSneu2(1)*(433*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(127*MSneu2(1) - 85*MStau2(2)) - 
     &         75*MStau2(2)**2))

	Sf8324 = -400*MStau2(1)**4 + 
     &   MSneu2(1)*(433*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(127*MSneu2(1) - 85*MStau2(1)) - 
     &         75*MStau2(1)**2))

	Sf7690 = 91*MStau2(2)**4 + 
     &   MSneu2(1)*(28*MStau2(2)**3 + 
     &      MSneu2(1)*(24*MStau2(2)**2 + 
     &         MSneu2(1)*(97*MSneu2(1) + 16*MStau2(2))))

	Sf7676 = 91*MStau2(1)**4 + 
     &   MSneu2(1)*(28*MStau2(1)**3 + 
     &      MSneu2(1)*(24*MStau2(1)**2 + 
     &         MSneu2(1)*(97*MSneu2(1) + 16*MStau2(1))))

	Sf1472 = -70*MStau2(2)**4 + 
     &   MSneu2(1)*(52*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(50*MSneu2(1) - 28*MStau2(2)) - 
     &         4*MStau2(2)**2))

	Sf571 = -70*MStau2(1)**4 + 
     &   MSneu2(1)*(52*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(50*MSneu2(1) - 28*MStau2(1)) - 
     &         4*MStau2(1)**2))

	Sf1476 = -35*MStau2(2)**4 + 
     &   MSneu2(1)*(26*MStau2(2)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(25*MSneu2(1) - 14*MStau2(2)) - 
     &         2*MStau2(2)**2))

	Sf575 = -35*MStau2(1)**4 + 
     &   MSneu2(1)*(26*MStau2(1)**3 + 
     &      MSneu2(1)*(MSneu2(1)*(25*MSneu2(1) - 14*MStau2(1)) - 
     &         2*MStau2(1)**2))

	Sf8113 = 49443*MStop2(2)**4 + 
     &   MSbot2(2)*(-12984*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (37455*MSbot2(2) - 11688*MStop2(2)) - 
     &         14786*MStop2(2)**2))

	Sf7904 = 49443*MStop2(1)**4 + 
     &   MSbot2(2)*(-12984*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (37455*MSbot2(2) - 11688*MStop2(1)) - 
     &         14786*MStop2(1)**2))

	Sf8136 = 37815*MStop2(2)**4 + 
     &   MSbot2(2)*(2352*MStop2(2)**3 + 
     &      MSbot2(2)*(3730*MStop2(2)**2 + 
     &         MSbot2(2)*(33279*MSbot2(2) + 3000*MStop2(2))))

	Sf7927 = 37815*MStop2(1)**4 + 
     &   MSbot2(2)*(2352*MStop2(1)**3 + 
     &      MSbot2(2)*(3730*MStop2(1)**2 + 
     &         MSbot2(2)*(33279*MSbot2(2) + 3000*MStop2(1))))

	Sf8394 = -4592*MStop2(2)**4 + 
     &   MSbot2(2)*(1899*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (6261*MSbot2(2) - 3047*MStop2(2)) - 
     &         521*MStop2(2)**2))

	Sf8370 = -4592*MStop2(1)**4 + 
     &   MSbot2(2)*(1899*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (6261*MSbot2(2) - 3047*MStop2(1)) - 
     &         521*MStop2(1)**2))

	Sf8128 = 7181*MStop2(2)**4 + 
     &   MSbot2(2)*(-4792*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (4913*MSbot2(2) - 4360*MStop2(2)) - 
     &         2830*MStop2(2)**2))

	Sf7919 = 7181*MStop2(1)**4 + 
     &   MSbot2(2)*(-4792*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (4913*MSbot2(2) - 4360*MStop2(1)) - 
     &         2830*MStop2(1)**2))

	Sf8459 = -7503*MStop2(2)**4 + 
     &   MSbot2(2)*(4397*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (3728*MSbot2(2) - 1305*MStop2(2)) + 
     &         683*MStop2(2)**2))

	Sf8435 = -7503*MStop2(1)**4 + 
     &   MSbot2(2)*(4397*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*
     &          (3728*MSbot2(2) - 1305*MStop2(1)) + 
     &         683*MStop2(1)**2))

	Sf8120 = 1763*MStop2(2)**4 + 
     &   MSbot2(2)*(6515*MStop2(2)**3 + 
     &      MSbot2(2)*(-3941*MStop2(2)**2 + 
     &         MSbot2(2)*(3626*MSbot2(2) + 2789*MStop2(2))))

	Sf7911 = 1763*MStop2(1)**4 + 
     &   MSbot2(2)*(6515*MStop2(1)**3 + 
     &      MSbot2(2)*(-3941*MStop2(1)**2 + 
     &         MSbot2(2)*(3626*MSbot2(2) + 2789*MStop2(1))))

	Sf5960 = 2350*MStop2(2)**4 + 
     &   MSbot2(2)*(718*MStop2(2)**3 + 
     &      MSbot2(2)*(792*MStop2(2)**2 + 
     &         MSbot2(2)*(2350*MSbot2(2) + 718*MStop2(2))))

	Sf3535 = 2350*MStop2(1)**4 + 
     &   MSbot2(2)*(718*MStop2(1)**3 + 
     &      MSbot2(2)*(792*MStop2(1)**2 + 
     &         MSbot2(2)*(2350*MSbot2(2) + 718*MStop2(1))))

	Sf8401 = 886*MStop2(2)**4 + 
     &   MSbot2(2)*(526*MStop2(2)**3 + 
     &      MSbot2(2)*(302*MStop2(2)**2 + 
     &         MSbot2(2)*(940*MSbot2(2) + 418*MStop2(2))))

	Sf8377 = 886*MStop2(1)**4 + 
     &   MSbot2(2)*(526*MStop2(1)**3 + 
     &      MSbot2(2)*(302*MStop2(1)**2 + 
     &         MSbot2(2)*(940*MSbot2(2) + 418*MStop2(1))))

	Sf8221 = -389*MStop2(2)**4 + 
     &   MSbot2(2)*(337*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(670*MSbot2(2) - 701*MStop2(2)) + 
     &         83*MStop2(2)**2))

	Sf8209 = -389*MStop2(1)**4 + 
     &   MSbot2(2)*(337*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(670*MSbot2(2) - 701*MStop2(1)) + 
     &         83*MStop2(1)**2))

	Sf8152 = 555*MStop2(2)**4 + 
     &   MSbot2(2)*(-504*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(375*MSbot2(2) - 360*MStop2(2)) - 
     &         50*MStop2(2)**2))

	Sf7943 = 555*MStop2(1)**4 + 
     &   MSbot2(2)*(-504*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(375*MSbot2(2) - 360*MStop2(1)) - 
     &         50*MStop2(1)**2))

	Sf8192 = -850*MStop2(2)**4 + 
     &   MSbot2(2)*(971*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(299*MSbot2(2) - 247*MStop2(2)) - 
     &         173*MStop2(2)**2))

	Sf7983 = -850*MStop2(1)**4 + 
     &   MSbot2(2)*(971*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(299*MSbot2(2) - 247*MStop2(1)) - 
     &         173*MStop2(1)**2))

	Sf5905 = 275*MStop2(2)**4 + 
     &   MSbot2(2)*(53*MStop2(2)**3 + 
     &      MSbot2(2)*(72*MStop2(2)**2 + 
     &         MSbot2(2)*(275*MSbot2(2) + 53*MStop2(2))))

	Sf3461 = 275*MStop2(1)**4 + 
     &   MSbot2(2)*(53*MStop2(1)**3 + 
     &      MSbot2(2)*(72*MStop2(1)**2 + 
     &         MSbot2(2)*(275*MSbot2(2) + 53*MStop2(1))))

	Sf5969 = -203*MStop2(2)**4 + 
     &   MSbot2(2)*(-115*MStop2(2)**3 + 
     &      MSbot2(2)*(24*MStop2(2)**2 + 
     &         MSbot2(2)*(245*MSbot2(2) + 49*MStop2(2))))

	Sf3544 = -203*MStop2(1)**4 + 
     &   MSbot2(2)*(-115*MStop2(1)**3 + 
     &      MSbot2(2)*(24*MStop2(1)**2 + 
     &         MSbot2(2)*(245*MSbot2(2) + 49*MStop2(1))))

	Sf5913 = -183*MStop2(2)**4 + 
     &   MSbot2(2)*(45*MStop2(2)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(169*MSbot2(2) - 47*MStop2(2)) + 
     &         16*MStop2(2)**2))

	Sf3469 = -183*MStop2(1)**4 + 
     &   MSbot2(2)*(45*MStop2(1)**3 + 
     &      MSbot2(2)*(MSbot2(2)*(169*MSbot2(2) - 47*MStop2(1)) + 
     &         16*MStop2(1)**2))

	Sf8035 = 49443*MStop2(2)**4 + 
     &   MSbot2(1)*(-12984*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (37455*MSbot2(1) - 11688*MStop2(2)) - 
     &         14786*MStop2(2)**2))

	Sf7826 = 49443*MStop2(1)**4 + 
     &   MSbot2(1)*(-12984*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (37455*MSbot2(1) - 11688*MStop2(1)) - 
     &         14786*MStop2(1)**2))

	Sf8058 = 37815*MStop2(2)**4 + 
     &   MSbot2(1)*(2352*MStop2(2)**3 + 
     &      MSbot2(1)*(3730*MStop2(2)**2 + 
     &         MSbot2(1)*(33279*MSbot2(1) + 3000*MStop2(2))))

	Sf7849 = 37815*MStop2(1)**4 + 
     &   MSbot2(1)*(2352*MStop2(1)**3 + 
     &      MSbot2(1)*(3730*MStop2(1)**2 + 
     &         MSbot2(1)*(33279*MSbot2(1) + 3000*MStop2(1))))

	Sf8298 = -4592*MStop2(2)**4 + 
     &   MSbot2(1)*(1899*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (6261*MSbot2(1) - 3047*MStop2(2)) - 
     &         521*MStop2(2)**2))

	Sf8274 = -4592*MStop2(1)**4 + 
     &   MSbot2(1)*(1899*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (6261*MSbot2(1) - 3047*MStop2(1)) - 
     &         521*MStop2(1)**2))

	Sf8050 = 7181*MStop2(2)**4 + 
     &   MSbot2(1)*(-4792*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (4913*MSbot2(1) - 4360*MStop2(2)) - 
     &         2830*MStop2(2)**2))

	Sf7841 = 7181*MStop2(1)**4 + 
     &   MSbot2(1)*(-4792*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (4913*MSbot2(1) - 4360*MStop2(1)) - 
     &         2830*MStop2(1)**2))

	Sf8448 = -7503*MStop2(2)**4 + 
     &   MSbot2(1)*(4397*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (3728*MSbot2(1) - 1305*MStop2(2)) + 
     &         683*MStop2(2)**2))

	Sf8424 = -7503*MStop2(1)**4 + 
     &   MSbot2(1)*(4397*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*
     &          (3728*MSbot2(1) - 1305*MStop2(1)) + 
     &         683*MStop2(1)**2))

	Sf8042 = 1763*MStop2(2)**4 + 
     &   MSbot2(1)*(6515*MStop2(2)**3 + 
     &      MSbot2(1)*(-3941*MStop2(2)**2 + 
     &         MSbot2(1)*(3626*MSbot2(1) + 2789*MStop2(2))))

	Sf7833 = 1763*MStop2(1)**4 + 
     &   MSbot2(1)*(6515*MStop2(1)**3 + 
     &      MSbot2(1)*(-3941*MStop2(1)**2 + 
     &         MSbot2(1)*(3626*MSbot2(1) + 2789*MStop2(1))))

	Sf5534 = 2350*MStop2(2)**4 + 
     &   MSbot2(1)*(718*MStop2(2)**3 + 
     &      MSbot2(1)*(792*MStop2(2)**2 + 
     &         MSbot2(1)*(2350*MSbot2(1) + 718*MStop2(2))))

	Sf3173 = 2350*MStop2(1)**4 + 
     &   MSbot2(1)*(718*MStop2(1)**3 + 
     &      MSbot2(1)*(792*MStop2(1)**2 + 
     &         MSbot2(1)*(2350*MSbot2(1) + 718*MStop2(1))))

	Sf8305 = 886*MStop2(2)**4 + 
     &   MSbot2(1)*(526*MStop2(2)**3 + 
     &      MSbot2(1)*(302*MStop2(2)**2 + 
     &         MSbot2(1)*(940*MSbot2(1) + 418*MStop2(2))))

	Sf8281 = 886*MStop2(1)**4 + 
     &   MSbot2(1)*(526*MStop2(1)**3 + 
     &      MSbot2(1)*(302*MStop2(1)**2 + 
     &         MSbot2(1)*(940*MSbot2(1) + 418*MStop2(1))))

	Sf8250 = -389*MStop2(2)**4 + 
     &   MSbot2(1)*(337*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(670*MSbot2(1) - 701*MStop2(2)) + 
     &         83*MStop2(2)**2))

	Sf8238 = -389*MStop2(1)**4 + 
     &   MSbot2(1)*(337*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(670*MSbot2(1) - 701*MStop2(1)) + 
     &         83*MStop2(1)**2))

	Sf8074 = 555*MStop2(2)**4 + 
     &   MSbot2(1)*(-504*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(375*MSbot2(1) - 360*MStop2(2)) - 
     &         50*MStop2(2)**2))

	Sf7865 = 555*MStop2(1)**4 + 
     &   MSbot2(1)*(-504*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(375*MSbot2(1) - 360*MStop2(1)) - 
     &         50*MStop2(1)**2))

	Sf8172 = -850*MStop2(2)**4 + 
     &   MSbot2(1)*(971*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(299*MSbot2(1) - 247*MStop2(2)) - 
     &         173*MStop2(2)**2))

	Sf7963 = -850*MStop2(1)**4 + 
     &   MSbot2(1)*(971*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(299*MSbot2(1) - 247*MStop2(1)) - 
     &         173*MStop2(1)**2))

	Sf5471 = 275*MStop2(2)**4 + 
     &   MSbot2(1)*(53*MStop2(2)**3 + 
     &      MSbot2(1)*(72*MStop2(2)**2 + 
     &         MSbot2(1)*(275*MSbot2(1) + 53*MStop2(2))))

	Sf3074 = 275*MStop2(1)**4 + 
     &   MSbot2(1)*(53*MStop2(1)**3 + 
     &      MSbot2(1)*(72*MStop2(1)**2 + 
     &         MSbot2(1)*(275*MSbot2(1) + 53*MStop2(1))))

	Sf5543 = -203*MStop2(2)**4 + 
     &   MSbot2(1)*(-115*MStop2(2)**3 + 
     &      MSbot2(1)*(24*MStop2(2)**2 + 
     &         MSbot2(1)*(245*MSbot2(1) + 49*MStop2(2))))

	Sf3182 = -203*MStop2(1)**4 + 
     &   MSbot2(1)*(-115*MStop2(1)**3 + 
     &      MSbot2(1)*(24*MStop2(1)**2 + 
     &         MSbot2(1)*(245*MSbot2(1) + 49*MStop2(1))))

	Sf5479 = -183*MStop2(2)**4 + 
     &   MSbot2(1)*(45*MStop2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(169*MSbot2(1) - 47*MStop2(2)) + 
     &         16*MStop2(2)**2))

	Sf3082 = -183*MStop2(1)**4 + 
     &   MSbot2(1)*(45*MStop2(1)**3 + 
     &      MSbot2(1)*(MSbot2(1)*(169*MSbot2(1) - 47*MStop2(1)) + 
     &         16*MStop2(1)**2))

	Sf9807 = 68*MStau2(2)**4 + 
     &   MZ2*(-614*MStau2(2)**3 + 
     &      MZ2*(MZ2*(337*MZ2 - 1599*MStau2(2)) - 257*MStau2(2)**2)
     &      )

	Sf9671 = 68*MStau2(1)**4 + 
     &   MZ2*(-614*MStau2(1)**3 + 
     &      MZ2*(MZ2*(337*MZ2 - 1599*MStau2(1)) - 257*MStau2(1)**2)
     &      )

	Sf9828 = -100*MStau2(2)**4 + 
     &   MZ2*(-122*MStau2(2)**3 + 
     &      MZ2*(MZ2*(55*MZ2 - 249*MStau2(2)) - 179*MStau2(2)**2))

	Sf9641 = -100*MStau2(1)**4 + 
     &   MZ2*(-122*MStau2(1)**3 + 
     &      MZ2*(MZ2*(55*MZ2 - 249*MStau2(1)) - 179*MStau2(1)**2))

	Sf10019 = 812*MStau2(2)**4 + 
     &   MZ2*(188*MStau2(2)**3 + 
     &      MZ2*(MZ2*(47*MZ2 - 1591*MStau2(2)) + 1435*MStau2(2)**2)
     &      )

	Sf9920 = 812*MStau2(1)**4 + 
     &   MZ2*(188*MStau2(1)**3 + 
     &      MZ2*(MZ2*(47*MZ2 - 1591*MStau2(1)) + 1435*MStau2(1)**2)
     &      )

	Sf9829 = -126*MStau2(2)**4 + 
     &   MZ2*(-153*MStau2(2)**3 + 
     &      MZ2*(MZ2*(9*MZ2 - 111*MStau2(2)) + 205*MStau2(2)**2))

	Sf9640 = -126*MStau2(1)**4 + 
     &   MZ2*(-153*MStau2(1)**3 + 
     &      MZ2*(MZ2*(9*MZ2 - 111*MStau2(1)) + 205*MStau2(1)**2))

	Sf10041 = 200*MStau2(2)**4 + 
     &   MZ2*(-46*MStau2(2)**3 + 
     &      MZ2*(51*MStau2(2)**2 + MZ2*(8*MZ2 + 86*MStau2(2))))

	Sf9957 = 200*MStau2(1)**4 + 
     &   MZ2*(-46*MStau2(1)**3 + 
     &      MZ2*(51*MStau2(1)**2 + MZ2*(8*MZ2 + 86*MStau2(1))))

	Sf9702 = 12*MStau2(2)**4 + 
     &   MZ2*(-15*MStau2(2)**3 + 
     &      MZ2*(-3*MStau2(2)**2 + MZ2*(6*MZ2 + 16*MStau2(2))))

	Sf9406 = 12*MStau2(1)**4 + 
     &   MZ2*(-15*MStau2(1)**3 + 
     &      MZ2*(-3*MStau2(1)**2 + MZ2*(6*MZ2 + 16*MStau2(1))))

	Sf9984 = 24*MStau2(2)**4 + 
     &   MZ2*(-3*MStau2(2)**3 + 
     &      MZ2*(-33*MStau2(2)**2 + MZ2*(5*MZ2 + 31*MStau2(2))))

	Sf9883 = 24*MStau2(1)**4 + 
     &   MZ2*(-3*MStau2(1)**3 + 
     &      MZ2*(-33*MStau2(1)**2 + MZ2*(5*MZ2 + 31*MStau2(1))))

	Sf9703 = 26*MStau2(2)**4 + 
     &   MZ2*(28*MStau2(2)**3 + 
     &      MZ2*(7*MStau2(2)**2 + MZ2*(4*MZ2 + 24*MStau2(2))))

	Sf9407 = 26*MStau2(1)**4 + 
     &   MZ2*(28*MStau2(1)**3 + 
     &      MZ2*(7*MStau2(1)**2 + MZ2*(4*MZ2 + 24*MStau2(1))))

	Sf9985 = 56*MStau2(2)**4 + 
     &   MZ2*(56*MStau2(2)**3 + 
     &      MZ2*(15*MStau2(2)**2 + MZ2*(2*MZ2 + 41*MStau2(2))))

	Sf9884 = 56*MStau2(1)**4 + 
     &   MZ2*(56*MStau2(1)**3 + 
     &      MZ2*(15*MStau2(1)**2 + MZ2*(2*MZ2 + 41*MStau2(1))))

	Sf9991 = -12*MStau2(2)**4 + 
     &   MZ2*(-40*MStau2(2)**3 + 
     &      MZ2*(MZ2*(MZ2 - 12*MStau2(2)) + 3*MStau2(2)**2))

	Sf10040 = -120*MStau2(2)**4 + 
     &   MZ2*(102*MStau2(2)**3 + 
     &      MZ2*(MZ2*(MZ2 - 256*MStau2(2)) + 345*MStau2(2)**2))

	Sf9874 = -12*MStau2(1)**4 + 
     &   MZ2*(-40*MStau2(1)**3 + 
     &      MZ2*(MZ2*(MZ2 - 12*MStau2(1)) + 3*MStau2(1)**2))

	Sf9956 = -120*MStau2(1)**4 + 
     &   MZ2*(102*MStau2(1)**3 + 
     &      MZ2*(MZ2*(MZ2 - 256*MStau2(1)) + 345*MStau2(1)**2))

	Sf8722 = 13*MZ2**4 + 
     &   MW2*(-97*MZ2**3 + 
     &      MW2*(MW2*(64*MW2 - 208*MZ2) + 228*MZ2**2))

	Sf3662 = 210*MSbot2(2)**4 + 
     &   MW2*(20*MSbot2(2)**3 + 
     &      MW2*(60*MSbot2(2)**2 + MW2*(50*MW2 + 61*MSbot2(2))))

	Sf3382 = 210*MSbot2(1)**4 + 
     &   MW2*(20*MSbot2(1)**3 + 
     &      MW2*(60*MSbot2(1)**2 + MW2*(50*MW2 + 61*MSbot2(1))))

	Sf8489 = -13*MZ2**4 + 
     &   MW2*(19*MZ2**3 + MW2*(MW2*(32*MW2 - 80*MZ2) + 42*MZ2**2))

	Sf7757 = 13*MZ2**4 + 
     &   MW2*(-71*MZ2**3 + 
     &      MW2*(MW2*(32*MW2 - 112*MZ2) + 138*MZ2**2))

	Sf8598 = 7*MZ2**4 + 
     &   MW2*(-37*MZ2**3 + MW2*(MW2*(16*MW2 - 56*MZ2) + 70*MZ2**2))

	Sf10048 = -448*MStau2(2)**4 + 
     &   MHiggs2(2)*(128*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (415*MHiggs2(2) - 1278*MStau2(2)) + 
     &         193*MStau2(2)**2))

	Sf9967 = -448*MStau2(1)**4 + 
     &   MHiggs2(2)*(128*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (415*MHiggs2(2) - 1278*MStau2(1)) + 
     &         193*MStau2(1)**2))

	Sf9800 = 866*MStau2(2)**4 + 
     &   MHiggs2(2)*(1296*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (284*MHiggs2(2) - 2163*MStau2(2)) + 
     &         1352*MStau2(2)**2))

	Sf9658 = 866*MStau2(1)**4 + 
     &   MHiggs2(2)*(1296*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (284*MHiggs2(2) - 2163*MStau2(1)) + 
     &         1352*MStau2(1)**2))

	Sf9717 = -196*MStau2(2)**4 + 
     &   MHiggs2(2)*(-1226*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (173*MHiggs2(2) - 716*MStau2(2)) + 
     &         1725*MStau2(2)**2))

	Sf9422 = -196*MStau2(1)**4 + 
     &   MHiggs2(2)*(-1226*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (173*MHiggs2(2) - 716*MStau2(1)) + 
     &         1725*MStau2(1)**2))

	Sf9715 = 64*MStau2(2)**4 + 
     &   MHiggs2(2)*(144*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (101*MHiggs2(2) - 314*MStau2(2)) + 5*MStau2(2)**2))

	Sf9421 = 64*MStau2(1)**4 + 
     &   MHiggs2(2)*(144*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (101*MHiggs2(2) - 314*MStau2(1)) + 5*MStau2(1)**2))

	Sf10011 = 108*MStau2(2)**4 + 
     &   MHiggs2(2)*(512*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (47*MHiggs2(2) - 2127*MStau2(2)) + 
     &         2589*MStau2(2)**2))

	Sf9910 = 108*MStau2(1)**4 + 
     &   MHiggs2(2)*(512*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (47*MHiggs2(2) - 2127*MStau2(1)) + 
     &         2589*MStau2(1)**2))

	Sf9821 = 230*MStau2(2)**4 + 
     &   MHiggs2(2)*(360*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (44*MHiggs2(2) - 309*MStau2(2)) + 104*MStau2(2)**2)
     &      )

	Sf9628 = 230*MStau2(1)**4 + 
     &   MHiggs2(2)*(360*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (44*MHiggs2(2) - 309*MStau2(1)) + 104*MStau2(1)**2)
     &      )

	Sf9771 = 16*MStau2(2)**4 + 
     &   MHiggs2(2)*(36*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (23*MHiggs2(2) - 74*MStau2(2)) - MStau2(2)**2))

	Sf9773 = -76*MStau2(2)**4 + 
     &   MHiggs2(2)*(-266*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (23*MHiggs2(2) - 80*MStau2(2)) + 303*MStau2(2)**2))

	Sf9534 = 16*MStau2(1)**4 + 
     &   MHiggs2(2)*(36*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (23*MHiggs2(2) - 74*MStau2(1)) - MStau2(1)**2))

	Sf9535 = -76*MStau2(1)**4 + 
     &   MHiggs2(2)*(-266*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*
     &          (23*MHiggs2(2) - 80*MStau2(1)) + 303*MStau2(1)**2))

	Sf9997 = 60*MStau2(2)**4 + 
     &   MHiggs2(2)*(25*MStau2(2)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*(7*MHiggs2(2) - 54*MStau2(2)) + 
     &         60*MStau2(2)**2))

	Sf9894 = 60*MStau2(1)**4 + 
     &   MHiggs2(2)*(25*MStau2(1)**3 + 
     &      MHiggs2(2)*(MHiggs2(2)*(7*MHiggs2(2) - 54*MStau2(1)) + 
     &         60*MStau2(1)**2))

	Sf9558 = 142*MZ2**4 + 
     &   MHiggs2(2)*(-345*MZ2**3 + 
     &      MHiggs2(2)*(218*MZ2**2 + 
     &         MHiggs2(2)*(17*MZ2 + 6*MHiggs2(2))))

	Sf9571 = 366*MZ2**4 + 
     &   MHiggs2(2)*(-1475*MZ2**3 + 
     &      MHiggs2(2)*(928*MZ2**2 + 
     &         MHiggs2(2)*(-303*MZ2 + 6*MHiggs2(2))))

	Sf10047 = -252*MStau2(2)**4 + 
     &   MHiggs2(2)*(-682*MStau2(2)**3 + 
     &      MHiggs2(2)*(-615*MStau2(2)**2 + 
     &         MHiggs2(2)*(5*MHiggs2(2) + 573*MStau2(2))))

	Sf9968 = -252*MStau2(1)**4 + 
     &   MHiggs2(2)*(-682*MStau2(1)**3 + 
     &      MHiggs2(2)*(-615*MStau2(1)**2 + 
     &         MHiggs2(2)*(5*MHiggs2(2) + 573*MStau2(1))))

	Sf9565 = 280*MZ2**4 + 
     &   MHiggs2(2)*(-53*MZ2**3 + 
     &      MHiggs2(2)*(-518*MZ2**2 + 
     &         MHiggs2(2)*(525*MZ2 + 4*MHiggs2(2))))

	Sf9550 = 58*MZ2**4 + 
     &   MHiggs2(2)*(83*MZ2**3 + 
     &      MHiggs2(2)*(290*MZ2**2 + 
     &         MHiggs2(2)*(39*MZ2 + 4*MHiggs2(2))))

	Sf9690 = 13*MStau2(2)**4 + 
     &   MHiggs2(2)*(-9*MStau2(2)**3 + 
     &      MHiggs2(2)*(-9*MStau2(2)**2 + 
     &         MHiggs2(2)*(2*MHiggs2(2) + 17*MStau2(2))))

	Sf9390 = 13*MStau2(1)**4 + 
     &   MHiggs2(2)*(-9*MStau2(1)**3 + 
     &      MHiggs2(2)*(-9*MStau2(1)**2 + 
     &         MHiggs2(2)*(2*MHiggs2(2) + 17*MStau2(1))))

	Sf9691 = 31*MStau2(2)**4 + 
     &   MHiggs2(2)*(21*MStau2(2)**3 + 
     &      MHiggs2(2)*(19*MStau2(2)**2 + 
     &         MHiggs2(2)*(MHiggs2(2) + 9*MStau2(2))))

	Sf9391 = 31*MStau2(1)**4 + 
     &   MHiggs2(2)*(21*MStau2(1)**3 + 
     &      MHiggs2(2)*(19*MStau2(1)**2 + 
     &         MHiggs2(2)*(MHiggs2(2) + 9*MStau2(1))))

	Sf8863 = 2278 - 1449*USbot22**2

	Sf8848 = 286 - 120*USbot22**2

	Sf8638 = 28 - 21*USbot22**2

	Sf8854 = -1 + 144*USbot22**2

	Sf8556 = -4 + 3*USbot22**2

	Sf8631 = 3 - 4*USbot22**2

	Sf1770 = (2*ML2*SA - CB*SAB*Ss443)*UStau21 + 
     &   ML*(CA*reMUE + SA*Re(Atau))*UStau22

	Sf1746 = (2*ML2*SA - CB*SAB*Ss443)*UStau11 + 
     &   ML*(CA*reMUE + SA*Re(Atau))*UStau12

	Sf7377 = 2 - 6*USbot22**2

	Sf8849 = -126*USbot22**2

	Sf7357 = -4 + 6*USbot22**2

	Sf10313 = MStop2(2)**4 + 
     &   MStop2(1)*(-4*MStop2(2)**3 + 
     &      MStop2(1)*(MStop2(1)**2 + 6*MStop2(2)**2))

	Sf10304 = MStop2(2)**4 + 
     &   MStop2(1)**2*(6*MStop2(2)**2 + 
     &      MStop2(1)*(-3*MStop2(2) - MStop2(3)))

	Sf10153 = MStau2(2)**4 + 
     &   MStau2(1)*(-4*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)**2 + 6*MStau2(2)**2))

	Sf10126 = MStau2(2)**4 + 
     &   MStau2(1)**2*(6*MStau2(2)**2 + 
     &      MStau2(1)*(-3*MStau2(2) - MStau2(3)))

	Sf1279 = -353*MStau2(2)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(353*MSneu2(1) - 188*MStau2(2)) + 
     &      188*MStau2(2)**3)

	Sf319 = -353*MStau2(1)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(353*MSneu2(1) - 188*MStau2(1)) + 
     &      188*MStau2(1)**3)

	Sf1646 = -65*MStau2(2)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(65*MSneu2(1) - 32*MStau2(2)) + 
     &      32*MStau2(2)**3)

	Sf821 = -65*MStau2(1)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(65*MSneu2(1) - 32*MStau2(1)) + 
     &      32*MStau2(1)**3)

	Sf1960 = -17*MStau2(2)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(35*MSneu2(1) - 23*MStau2(2)) + 
     &      5*MStau2(2)**3)

	Sf1898 = -17*MStau2(1)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(35*MSneu2(1) - 23*MStau2(1)) + 
     &      5*MStau2(1)**3)

	Sf1846 = -13*MStau2(2)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(13*MSneu2(1) - 7*MStau2(2)) + 
     &      7*MStau2(2)**3)

	Sf1067 = -13*MStau2(1)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(13*MSneu2(1) - 7*MStau2(1)) + 
     &      7*MStau2(1)**3)

	Sf1464 = 5*MStau2(2)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(5*MSneu2(1) - MStau2(2)) - 
     &      MStau2(2)**3)

	Sf561 = 5*MStau2(1)**4 + 
     &   MSneu2(1)*(MSneu2(1)**2*(5*MSneu2(1) - MStau2(1)) - 
     &      MStau2(1)**3)

	Sf7626 = 14555*MStau2(2)**3 + 
     &   MSneu2(1)*(14109*MStau2(2)**2 + 
     &      MSneu2(1)*(12719*MSneu2(1) + 13353*MStau2(2)))

	Sf7568 = 14555*MStau2(1)**3 + 
     &   MSneu2(1)*(14109*MStau2(1)**2 + 
     &      MSneu2(1)*(12719*MSneu2(1) + 13353*MStau2(1)))

	Sf7595 = 5033*MStau2(2)**3 + 
     &   MSneu2(1)*(4716*MStau2(2)**2 + 
     &      MSneu2(1)*(4061*MSneu2(1) + 4230*MStau2(2)))

	Sf7537 = 5033*MStau2(1)**3 + 
     &   MSneu2(1)*(4716*MStau2(1)**2 + 
     &      MSneu2(1)*(4061*MSneu2(1) + 4230*MStau2(1)))

	Sf7596 = 6643*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(2269*MSneu2(1) - 1638*MStau2(2)) - 
     &      180*MStau2(2)**2)

	Sf7538 = 6643*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(2269*MSneu2(1) - 1638*MStau2(1)) - 
     &      180*MStau2(1)**2)

	Sf7696 = -1319*MStau2(2)**3 + 
     &   MSneu2(1)*(-195*MStau2(2)**2 + 
     &      MSneu2(1)*(1705*MSneu2(1) + 297*MStau2(2)))

	Sf7682 = -1319*MStau2(1)**3 + 
     &   MSneu2(1)*(-195*MStau2(1)**2 + 
     &      MSneu2(1)*(1705*MSneu2(1) + 297*MStau2(1)))

	Sf7625 = 1449*MStau2(2)**3 + 
     &   MSneu2(1)*(670*MStau2(2)**2 + 
     &      MSneu2(1)*(963*MSneu2(1) + 508*MStau2(2)))

	Sf7567 = 1449*MStau2(1)**3 + 
     &   MSneu2(1)*(670*MStau2(1)**2 + 
     &      MSneu2(1)*(963*MSneu2(1) + 508*MStau2(1)))

	Sf7610 = -2015*MStau2(2)**3 + 
     &   MSneu2(1)*(5883*MStau2(2)**2 + 
     &      MSneu2(1)*(901*MSneu2(1) + 4911*MStau2(2)))

	Sf7552 = -2015*MStau2(1)**3 + 
     &   MSneu2(1)*(5883*MStau2(1)**2 + 
     &      MSneu2(1)*(901*MSneu2(1) + 4911*MStau2(1)))

	Sf7743 = -2857*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(707*MSneu2(1) - 165*MStau2(2)) - 
     &      477*MStau2(2)**2)

	Sf7719 = -2857*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(707*MSneu2(1) - 165*MStau2(1)) - 
     &      477*MStau2(1)**2)

	Sf7604 = 935*MStau2(2)**3 + 
     &   MSneu2(1)*(648*MStau2(2)**2 + 
     &      MSneu2(1)*(683*MSneu2(1) + 558*MStau2(2)))

	Sf7546 = 935*MStau2(1)**3 + 
     &   MSneu2(1)*(648*MStau2(1)**2 + 
     &      MSneu2(1)*(683*MSneu2(1) + 558*MStau2(1)))

	Sf1338 = 100*MStau2(2)**3 + 
     &   MSneu2(1)*(90*MStau2(2)**2 + 
     &      MSneu2(1)*(601*MSneu2(1) + 237*MStau2(2)))

	Sf401 = 100*MStau2(1)**3 + 
     &   MSneu2(1)*(90*MStau2(1)**2 + 
     &      MSneu2(1)*(601*MSneu2(1) + 237*MStau2(1)))

	Sf7689 = 701*MStau2(2)**3 + 
     &   MSneu2(1)*(495*MStau2(2)**2 + 
     &      MSneu2(1)*(593*MSneu2(1) + 531*MStau2(2)))

	Sf7675 = 701*MStau2(1)**3 + 
     &   MSneu2(1)*(495*MStau2(1)**2 + 
     &      MSneu2(1)*(593*MSneu2(1) + 531*MStau2(1)))

	Sf1337 = 10*MStau2(2)**3 + 
     &   MSneu2(1)*(189*MStau2(2)**2 + 
     &      MSneu2(1)*(574*MSneu2(1) + 255*MStau2(2)))

	Sf400 = 10*MStau2(1)**3 + 
     &   MSneu2(1)*(189*MStau2(1)**2 + 
     &      MSneu2(1)*(574*MSneu2(1) + 255*MStau2(1)))

	Sf7660 = 587*MStau2(2)**3 + 
     &   MSneu2(1)*(9*MStau2(2)**2 + 
     &      MSneu2(1)*(407*MSneu2(1) + 69*MStau2(2)))

	Sf7648 = 587*MStau2(1)**3 + 
     &   MSneu2(1)*(9*MStau2(1)**2 + 
     &      MSneu2(1)*(407*MSneu2(1) + 69*MStau2(1)))

	Sf1645 = 457*MStau2(2)**3 + 
     &   MSneu2(1)*(155*MStau2(2)**2 + 
     &      MSneu2(1)*(349*MSneu2(1) + 191*MStau2(2)))

	Sf820 = 457*MStau2(1)**3 + 
     &   MSneu2(1)*(155*MStau2(1)**2 + 
     &      MSneu2(1)*(349*MSneu2(1) + 191*MStau2(1)))

	Sf7611 = 215*MStau2(2)**3 + 
     &   MSneu2(1)*(583*MStau2(2)**2 + 
     &      MSneu2(1)*(341*MSneu2(1) + 205*MStau2(2)))

	Sf7553 = 215*MStau2(1)**3 + 
     &   MSneu2(1)*(583*MStau2(1)**2 + 
     &      MSneu2(1)*(341*MSneu2(1) + 205*MStau2(1)))

	Sf1343 = 76*MStau2(2)**3 + 
     &   MSneu2(1)*(327*MStau2(2)**2 + 
     &      MSneu2(1)*(327*MSneu2(1) + 294*MStau2(2)))

	Sf406 = 76*MStau2(1)**3 + 
     &   MSneu2(1)*(327*MStau2(1)**2 + 
     &      MSneu2(1)*(327*MSneu2(1) + 294*MStau2(1)))

	Sf1411 = -387*MStau2(2)**3 + 
     &   MSneu2(1)*(-165*MStau2(2)**2 + 
     &      MSneu2(1)*(303*MSneu2(1) + 9*MStau2(2)))

	Sf494 = -387*MStau2(1)**3 + 
     &   MSneu2(1)*(-165*MStau2(1)**2 + 
     &      MSneu2(1)*(303*MSneu2(1) + 9*MStau2(1)))

	Sf2145 = -230*MStau2(2)**3 + 
     &   MSneu2(1)*(-30*MStau2(2)**2 + 
     &      MSneu2(1)*(280*MSneu2(1) + 84*MStau2(2)))

	Sf2065 = -230*MStau2(1)**3 + 
     &   MSneu2(1)*(-30*MStau2(1)**2 + 
     &      MSneu2(1)*(280*MSneu2(1) + 84*MStau2(1)))

	Sf7618 = 325*MStau2(2)**3 + 
     &   MSneu2(1)*(435*MStau2(2)**2 + 
     &      MSneu2(1)*(265*MSneu2(1) + 399*MStau2(2)))

	Sf7560 = 325*MStau2(1)**3 + 
     &   MSneu2(1)*(435*MStau2(1)**2 + 
     &      MSneu2(1)*(265*MSneu2(1) + 399*MStau2(1)))

	Sf1704 = 199*MStau2(2)**3 + 
     &   MSneu2(1)*(345*MStau2(2)**2 + 
     &      MSneu2(1)*(253*MSneu2(1) + 291*MStau2(2)))

	Sf905 = 199*MStau2(1)**3 + 
     &   MSneu2(1)*(345*MStau2(1)**2 + 
     &      MSneu2(1)*(253*MSneu2(1) + 291*MStau2(1)))

	Sf7603 = 841*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(247*MSneu2(1) - 522*MStau2(2)) - 
     &      324*MStau2(2)**2)

	Sf7545 = 841*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(247*MSneu2(1) - 522*MStau2(1)) - 
     &      324*MStau2(1)**2)

	Sf2149 = 219*MStau2(2)**3 + 
     &   MSneu2(1)*(54*MStau2(2)**2 + 
     &      MSneu2(1)*(244*MSneu2(1) + 123*MStau2(2)))

	Sf2069 = 219*MStau2(1)**3 + 
     &   MSneu2(1)*(54*MStau2(1)**2 + 
     &      MSneu2(1)*(244*MSneu2(1) + 123*MStau2(1)))

	Sf1313 = 231*MStau2(2)**3 + 
     &   MSneu2(1)*(306*MStau2(2)**2 + 
     &      MSneu2(1)*(232*MSneu2(1) + 255*MStau2(2)))

	Sf362 = 231*MStau2(1)**3 + 
     &   MSneu2(1)*(306*MStau2(1)**2 + 
     &      MSneu2(1)*(232*MSneu2(1) + 255*MStau2(1)))

	Sf1851 = 309*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(201*MSneu2(1) - 45*MStau2(2)) - 
     &      81*MStau2(2)**2)

	Sf1072 = 309*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(201*MSneu2(1) - 45*MStau2(1)) - 
     &      81*MStau2(1)**2)

	Sf1393 = 211*MStau2(2)**3 + 
     &   MSneu2(1)*(165*MStau2(2)**2 + 
     &      MSneu2(1)*(193*MSneu2(1) + 135*MStau2(2)))

	Sf476 = 211*MStau2(1)**3 + 
     &   MSneu2(1)*(165*MStau2(1)**2 + 
     &      MSneu2(1)*(193*MSneu2(1) + 135*MStau2(1)))

	Sf1703 = -177*MStau2(2)**3 + 
     &   MSneu2(1)*(-167*MStau2(2)**2 + 
     &      MSneu2(1)*(181*MSneu2(1) + 83*MStau2(2)))

	Sf904 = -177*MStau2(1)**3 + 
     &   MSneu2(1)*(-167*MStau2(1)**2 + 
     &      MSneu2(1)*(181*MSneu2(1) + 83*MStau2(1)))

	Sf2144 = 160*MStau2(2)**3 + 
     &   MSneu2(1)*(99*MStau2(2)**2 + 
     &      MSneu2(1)*(160*MSneu2(1) + 105*MStau2(2)))

	Sf2064 = 160*MStau2(1)**3 + 
     &   MSneu2(1)*(99*MStau2(1)**2 + 
     &      MSneu2(1)*(160*MSneu2(1) + 105*MStau2(1)))

	Sf7617 = 275*MStau2(2)**3 + 
     &   MSneu2(1)*(219*MStau2(2)**2 + 
     &      MSneu2(1)*(155*MSneu2(1) + 171*MStau2(2)))

	Sf7665 = -79*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(155*MSneu2(1) - 51*MStau2(2)) + 
     &      15*MStau2(2)**2)

	Sf7559 = 275*MStau2(1)**3 + 
     &   MSneu2(1)*(219*MStau2(1)**2 + 
     &      MSneu2(1)*(155*MSneu2(1) + 171*MStau2(1)))

	Sf7653 = -79*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(155*MSneu2(1) - 51*MStau2(1)) + 
     &      15*MStau2(1)**2)

	Sf1392 = -129*MStau2(2)**3 + 
     &   MSneu2(1)*(-55*MStau2(2)**2 + 
     &      MSneu2(1)*(101*MSneu2(1) + 3*MStau2(2)))

	Sf475 = -129*MStau2(1)**3 + 
     &   MSneu2(1)*(-55*MStau2(1)**2 + 
     &      MSneu2(1)*(101*MSneu2(1) + 3*MStau2(1)))

	Sf1308 = 601*MStau2(2)**3 + 
     &   MSneu2(1)*(237*MStau2(2)**2 + 
     &      MSneu2(1)*(100*MSneu2(1) + 90*MStau2(2)))

	Sf357 = 601*MStau2(1)**3 + 
     &   MSneu2(1)*(237*MStau2(1)**2 + 
     &      MSneu2(1)*(100*MSneu2(1) + 90*MStau2(1)))

	Sf8330 = -275*MStau2(2)**3 + 
     &   MSneu2(1)*(135*MStau2(2)**2 + 
     &      MSneu2(1)*(67*MSneu2(1) + 33*MStau2(2)))

	Sf8323 = -275*MStau2(1)**3 + 
     &   MSneu2(1)*(135*MStau2(1)**2 + 
     &      MSneu2(1)*(67*MSneu2(1) + 33*MStau2(1)))

	Sf7688 = 53*MStau2(2)**3 + 
     &   MSneu2(1)*(49*MStau2(2)**2 + 
     &      MSneu2(1)*(59*MSneu2(1) + 31*MStau2(2)))

	Sf7674 = 53*MStau2(1)**3 + 
     &   MSneu2(1)*(49*MStau2(1)**2 + 
     &      MSneu2(1)*(59*MSneu2(1) + 31*MStau2(1)))

	Sf1634 = -48*MStau2(2)**3 + 
     &   MSneu2(1)*(-45*MStau2(2)**2 + 
     &      MSneu2(1)*(57*MSneu2(1) + 36*MStau2(2)))

	Sf809 = -48*MStau2(1)**3 + 
     &   MSneu2(1)*(-45*MStau2(1)**2 + 
     &      MSneu2(1)*(57*MSneu2(1) + 36*MStau2(1)))

	Sf1969 = 33*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(50*MSneu2(1) - 11*MStau2(2)) - 
     &      8*MStau2(2)**2)

	Sf1909 = 33*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(50*MSneu2(1) - 11*MStau2(1)) - 
     &      8*MStau2(1)**2)

	Sf2169 = -57*MStau2(2)**3 + 
     &   MSneu2(1)*(-36*MStau2(2)**2 + 
     &      MSneu2(1)*(48*MSneu2(1) + 45*MStau2(2)))

	Sf2100 = -57*MStau2(1)**3 + 
     &   MSneu2(1)*(-36*MStau2(1)**2 + 
     &      MSneu2(1)*(48*MSneu2(1) + 45*MStau2(1)))

	Sf7633 = 100*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(40*MSneu2(1) - 42*MStau2(2)) - 
     &      54*MStau2(2)**2)

	Sf7575 = 100*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(40*MSneu2(1) - 42*MStau2(1)) - 
     &      54*MStau2(1)**2)

	Sf1475 = -45*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(25*MSneu2(1) - MStau2(2)) + 
     &      5*MStau2(2)**2)

	Sf574 = -45*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(25*MSneu2(1) - MStau2(1)) + 
     &      5*MStau2(1)**2)

	Sf1630 = -16*MStau2(2)**3 + 
     &   MSneu2(1)*(-15*MStau2(2)**2 + 
     &      MSneu2(1)*(19*MSneu2(1) + 12*MStau2(2)))

	Sf805 = -16*MStau2(1)**3 + 
     &   MSneu2(1)*(-15*MStau2(1)**2 + 
     &      MSneu2(1)*(19*MSneu2(1) + 12*MStau2(1)))

	Sf2165 = -19*MStau2(2)**3 + 
     &   MSneu2(1)*(-12*MStau2(2)**2 + 
     &      MSneu2(1)*(16*MSneu2(1) + 15*MStau2(2)))

	Sf2096 = -19*MStau2(1)**3 + 
     &   MSneu2(1)*(-12*MStau2(1)**2 + 
     &      MSneu2(1)*(16*MSneu2(1) + 15*MStau2(1)))

	Sf7634 = 20*MStau2(2)**3 + 
     &   MSneu2(1)*(9*MStau2(2)**2 + 
     &      MSneu2(1)*(14*MSneu2(1) + 9*MStau2(2)))

	Sf7576 = 20*MStau2(1)**3 + 
     &   MSneu2(1)*(9*MStau2(1)**2 + 
     &      MSneu2(1)*(14*MSneu2(1) + 9*MStau2(1)))

	Sf1284 = -97*MStau2(2)**3 + 
     &   MSneu2(1)*(253*MStau2(2)**2 + 
     &      MSneu2(1)*(11*MSneu2(1) + 217*MStau2(2)))

	Sf1845 = -11*MStau2(2)**3 + 
     &   MSneu2(1)*(-MStau2(2)**2 + 
     &      MSneu2(1)*(11*MSneu2(1) + MStau2(2)))

	Sf324 = -97*MStau2(1)**3 + 
     &   MSneu2(1)*(253*MStau2(1)**2 + 
     &      MSneu2(1)*(11*MSneu2(1) + 217*MStau2(1)))

	Sf1066 = -11*MStau2(1)**3 + 
     &   MSneu2(1)*(-MStau2(1)**2 + 
     &      MSneu2(1)*(11*MSneu2(1) + MStau2(1)))

	Sf1307 = 574*MStau2(2)**3 + 
     &   MSneu2(1)*(255*MStau2(2)**2 + 
     &      MSneu2(1)*(10*MSneu2(1) + 189*MStau2(2)))

	Sf356 = 574*MStau2(1)**3 + 
     &   MSneu2(1)*(255*MStau2(1)**2 + 
     &      MSneu2(1)*(10*MSneu2(1) + 189*MStau2(1)))

	Sf1677 = -7*MStau2(2)**3 + 
     &   MSneu2(1)*(-5*MStau2(2)**2 + 
     &      MSneu2(1)*(8*MSneu2(1) + 4*MStau2(2)))

	Sf864 = -7*MStau2(1)**3 + 
     &   MSneu2(1)*(-5*MStau2(1)**2 + 
     &      MSneu2(1)*(8*MSneu2(1) + 4*MStau2(1)))

	Sf1834 = 13*MStau2(2)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(3*MSneu2(1) - MStau2(2)) - 
     &      15*MStau2(2)**2)

	Sf1045 = 13*MStau2(1)**3 + 
     &   MSneu2(1)*(MSneu2(1)*(3*MSneu2(1) - MStau2(1)) - 
     &      15*MStau2(1)**2)

	Sf8135 = 38197*MStop2(2)**3 + 
     &   MSbot2(2)*(38839*MStop2(2)**2 + 
     &      MSbot2(2)*(35443*MSbot2(2) + 37705*MStop2(2)))

	Sf7926 = 38197*MStop2(1)**3 + 
     &   MSbot2(2)*(38839*MStop2(1)**2 + 
     &      MSbot2(2)*(35443*MSbot2(2) + 37705*MStop2(1)))

	Sf8111 = 24494*MStop2(2)**3 + 
     &   MSbot2(2)*(23135*MStop2(2)**2 + 
     &      MSbot2(2)*(21578*MSbot2(2) + 21677*MStop2(2)))

	Sf7902 = 24494*MStop2(1)**3 + 
     &   MSbot2(2)*(23135*MStop2(1)**2 + 
     &      MSbot2(2)*(21578*MSbot2(2) + 21677*MStop2(1)))

	Sf8134 = 22179*MStop2(2)**3 + 
     &   MSbot2(2)*(12127*MStop2(2)**2 + 
     &      MSbot2(2)*(17805*MSbot2(2) + 10669*MStop2(2)))

	Sf7925 = 22179*MStop2(1)**3 + 
     &   MSbot2(2)*(12127*MStop2(1)**2 + 
     &      MSbot2(2)*(17805*MSbot2(2) + 10669*MStop2(1)))

	Sf8112 = 27099*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*
     &       (13977*MSbot2(2) - 11195*MStop2(2)) - 
     &      6821*MStop2(2)**2)

	Sf7903 = 27099*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*
     &       (13977*MSbot2(2) - 11195*MStop2(1)) - 
     &      6821*MStop2(1)**2)

	Sf8393 = -2751*MStop2(2)**3 + 
     &   MSbot2(2)*(-263*MStop2(2)**2 + 
     &      MSbot2(2)*(4637*MSbot2(2) + 689*MStop2(2)))

	Sf8369 = -2751*MStop2(1)**3 + 
     &   MSbot2(2)*(-263*MStop2(1)**2 + 
     &      MSbot2(2)*(4637*MSbot2(2) + 689*MStop2(1)))

	Sf8127 = 4230*MStop2(2)**3 + 
     &   MSbot2(2)*(2833*MStop2(2)**2 + 
     &      MSbot2(2)*(3474*MSbot2(2) + 2563*MStop2(2)))

	Sf7918 = 4230*MStop2(1)**3 + 
     &   MSbot2(2)*(2833*MStop2(1)**2 + 
     &      MSbot2(2)*(3474*MSbot2(2) + 2563*MStop2(1)))

	Sf8105 = 2179*MStop2(2)**3 + 
     &   MSbot2(2)*(-363*MStop2(2)**2 + 
     &      MSbot2(2)*(3313*MSbot2(2) + 15*MStop2(2)))

	Sf7896 = 2179*MStop2(1)**3 + 
     &   MSbot2(2)*(-363*MStop2(1)**2 + 
     &      MSbot2(2)*(3313*MSbot2(2) + 15*MStop2(1)))

	Sf8119 = 1915*MStop2(2)**3 + 
     &   MSbot2(2)*(2684*MStop2(2)**2 + 
     &      MSbot2(2)*(2482*MSbot2(2) + 983*MStop2(2)))

	Sf7910 = 1915*MStop2(1)**3 + 
     &   MSbot2(2)*(2684*MStop2(1)**2 + 
     &      MSbot2(2)*(2482*MSbot2(2) + 983*MStop2(1)))

	Sf8458 = -6365*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(1833*MSbot2(2) - 277*MStop2(2)) - 
     &      959*MStop2(2)**2)

	Sf8434 = -6365*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(1833*MSbot2(2) - 277*MStop2(1)) - 
     &      959*MStop2(1)**2)

	Sf8126 = 3289*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(1507*MSbot2(2) - 3001*MStop2(2)) - 
     &      2407*MStop2(2)**2)

	Sf7917 = 3289*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(1507*MSbot2(2) - 3001*MStop2(1)) - 
     &      2407*MStop2(1)**2)

	Sf8088 = 775*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(1207*MSbot2(2) - 693*MStop2(2)) - 
     &      801*MStop2(2)**2)

	Sf7879 = 775*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(1207*MSbot2(2) - 693*MStop2(1)) - 
     &      801*MStop2(1)**2)

	Sf8400 = 949*MStop2(2)**3 + 
     &   MSbot2(2)*(681*MStop2(2)**2 + 
     &      MSbot2(2)*(787*MSbot2(2) + 735*MStop2(2)))

	Sf8376 = 949*MStop2(1)**3 + 
     &   MSbot2(2)*(681*MStop2(1)**2 + 
     &      MSbot2(2)*(787*MSbot2(2) + 735*MStop2(1)))

	Sf8099 = 409*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(643*MSbot2(2) - 147*MStop2(2)) - 
     &      129*MStop2(2)**2)

	Sf7890 = 409*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(643*MSbot2(2) - 147*MStop2(1)) - 
     &      129*MStop2(1)**2)

	Sf5380 = 100*MStop2(2)**3 + 
     &   MSbot2(2)*(90*MStop2(2)**2 + 
     &      MSbot2(2)*(601*MSbot2(2) + 237*MStop2(2)))

	Sf2972 = 100*MStop2(1)**3 + 
     &   MSbot2(2)*(90*MStop2(1)**2 + 
     &      MSbot2(2)*(601*MSbot2(2) + 237*MStop2(1)))

	Sf5379 = 10*MStop2(2)**3 + 
     &   MSbot2(2)*(189*MStop2(2)**2 + 
     &      MSbot2(2)*(574*MSbot2(2) + 255*MStop2(2)))

	Sf2971 = 10*MStop2(1)**3 + 
     &   MSbot2(2)*(189*MStop2(1)**2 + 
     &      MSbot2(2)*(574*MSbot2(2) + 255*MStop2(1)))

	Sf8387 = -316*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(497*MSbot2(2) - 250*MStop2(2)) + 
     &      69*MStop2(2)**2)

	Sf8363 = -316*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(497*MSbot2(2) - 250*MStop2(1)) + 
     &      69*MStop2(1)**2)

	Sf8697 = 527*MStop2(2)**3 + 
     &   MSbot2(2)*(29*MStop2(2)**2 + 
     &      MSbot2(2)*(467*MSbot2(2) + 49*MStop2(2)))

	Sf8657 = 527*MStop2(1)**3 + 
     &   MSbot2(2)*(29*MStop2(1)**2 + 
     &      MSbot2(2)*(467*MSbot2(2) + 49*MStop2(1)))

	Sf8452 = -389*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(406*MSbot2(2) - 141*MStop2(2)) + 
     &      124*MStop2(2)**2)

	Sf8428 = -389*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(406*MSbot2(2) - 141*MStop2(1)) + 
     &      124*MStop2(1)**2)

	Sf8220 = -163*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(375*MSbot2(2) - 215*MStop2(2)) + 
     &      11*MStop2(2)**2)

	Sf8208 = -163*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(375*MSbot2(2) - 215*MStop2(1)) + 
     &      11*MStop2(1)**2)

	Sf8399 = 313*MStop2(2)**3 + 
     &   MSbot2(2)*(290*MStop2(2)**2 + 
     &      MSbot2(2)*(340*MSbot2(2) + 209*MStop2(2)))

	Sf8375 = 313*MStop2(1)**3 + 
     &   MSbot2(2)*(290*MStop2(1)**2 + 
     &      MSbot2(2)*(340*MSbot2(2) + 209*MStop2(1)))

	Sf8185 = 413*MStop2(2)**3 + 
     &   MSbot2(2)*(-23*MStop2(2)**2 + 
     &      MSbot2(2)*(323*MSbot2(2) + 7*MStop2(2)))

	Sf7976 = 413*MStop2(1)**3 + 
     &   MSbot2(2)*(-23*MStop2(1)**2 + 
     &      MSbot2(2)*(323*MSbot2(2) + 7*MStop2(1)))

	Sf5385 = 142*MStop2(2)**3 + 
     &   MSbot2(2)*(291*MStop2(2)**2 + 
     &      MSbot2(2)*(303*MSbot2(2) + 288*MStop2(2)))

	Sf2977 = 142*MStop2(1)**3 + 
     &   MSbot2(2)*(291*MStop2(1)**2 + 
     &      MSbot2(2)*(303*MSbot2(2) + 288*MStop2(1)))

	Sf8118 = -199*MStop2(2)**3 + 
     &   MSbot2(2)*(1645*MStop2(2)**2 + 
     &      MSbot2(2)*(287*MSbot2(2) + 1483*MStop2(2)))

	Sf7909 = -199*MStop2(1)**3 + 
     &   MSbot2(2)*(1645*MStop2(1)**2 + 
     &      MSbot2(2)*(287*MSbot2(2) + 1483*MStop2(1)))

	Sf6853 = -230*MStop2(2)**3 + 
     &   MSbot2(2)*(-30*MStop2(2)**2 + 
     &      MSbot2(2)*(280*MSbot2(2) + 84*MStop2(2)))

	Sf6739 = -230*MStop2(1)**3 + 
     &   MSbot2(2)*(-30*MStop2(1)**2 + 
     &      MSbot2(2)*(280*MSbot2(2) + 84*MStop2(1)))

	Sf6858 = 173*MStop2(2)**3 + 
     &   MSbot2(2)*(20*MStop2(2)**2 + 
     &      MSbot2(2)*(246*MSbot2(2) + 73*MStop2(2)))

	Sf6744 = 173*MStop2(1)**3 + 
     &   MSbot2(2)*(20*MStop2(1)**2 + 
     &      MSbot2(2)*(246*MSbot2(2) + 73*MStop2(1)))

	Sf7079 = -280*MStop2(2)**3 + 
     &   MSbot2(2)*(-84*MStop2(2)**2 + 
     &      MSbot2(2)*(230*MSbot2(2) + 30*MStop2(2)))

	Sf4648 = -280*MStop2(1)**3 + 
     &   MSbot2(2)*(-84*MStop2(1)**2 + 
     &      MSbot2(2)*(230*MSbot2(2) + 30*MStop2(1)))

	Sf5350 = 257*MStop2(2)**3 + 
     &   MSbot2(2)*(296*MStop2(2)**2 + 
     &      MSbot2(2)*(210*MSbot2(2) + 261*MStop2(2)))

	Sf2939 = 257*MStop2(1)**3 + 
     &   MSbot2(2)*(296*MStop2(1)**2 + 
     &      MSbot2(2)*(210*MSbot2(2) + 261*MStop2(1)))

	Sf7085 = 234*MStop2(2)**3 + 
     &   MSbot2(2)*(79*MStop2(2)**2 + 
     &      MSbot2(2)*(179*MSbot2(2) + 20*MStop2(2)))

	Sf4654 = 234*MStop2(1)**3 + 
     &   MSbot2(2)*(79*MStop2(1)**2 + 
     &      MSbot2(2)*(179*MSbot2(2) + 20*MStop2(1)))

	Sf6852 = 160*MStop2(2)**3 + 
     &   MSbot2(2)*(99*MStop2(2)**2 + 
     &      MSbot2(2)*(160*MSbot2(2) + 105*MStop2(2)))

	Sf7080 = 160*MStop2(2)**3 + 
     &   MSbot2(2)*(105*MStop2(2)**2 + 
     &      MSbot2(2)*(160*MSbot2(2) + 99*MStop2(2)))

	Sf6738 = 160*MStop2(1)**3 + 
     &   MSbot2(2)*(99*MStop2(1)**2 + 
     &      MSbot2(2)*(160*MSbot2(2) + 105*MStop2(1)))

	Sf4649 = 160*MStop2(1)**3 + 
     &   MSbot2(2)*(105*MStop2(1)**2 + 
     &      MSbot2(2)*(160*MSbot2(2) + 99*MStop2(1)))

	Sf8191 = -555*MStop2(2)**3 + 
     &   MSbot2(2)*(341*MStop2(2)**2 + 
     &      MSbot2(2)*(145*MSbot2(2) + 61*MStop2(2)))

	Sf7982 = -555*MStop2(1)**3 + 
     &   MSbot2(2)*(341*MStop2(1)**2 + 
     &      MSbot2(2)*(145*MSbot2(2) + 61*MStop2(1)))

	Sf5345 = 601*MStop2(2)**3 + 
     &   MSbot2(2)*(237*MStop2(2)**2 + 
     &      MSbot2(2)*(100*MSbot2(2) + 90*MStop2(2)))

	Sf2934 = 601*MStop2(1)**3 + 
     &   MSbot2(2)*(237*MStop2(1)**2 + 
     &      MSbot2(2)*(100*MSbot2(2) + 90*MStop2(1)))

	Sf8150 = 189*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(99*MSbot2(2) - 125*MStop2(2)) - 
     &      143*MStop2(2)**2)

	Sf7941 = 189*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(99*MSbot2(2) - 125*MStop2(1)) - 
     &      143*MStop2(1)**2)

	Sf8337 = -30*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(87*MSbot2(2) - 72*MStop2(2)) + 
     &      15*MStop2(2)**2)

	Sf8335 = -30*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(87*MSbot2(2) - 72*MStop2(1)) + 
     &      15*MStop2(1)**2)

	Sf5968 = -61*MStop2(2)**3 + 
     &   MSbot2(2)*(-51*MStop2(2)**2 + 
     &      MSbot2(2)*(68*MSbot2(2) + 39*MStop2(2)))

	Sf8093 = -121*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(68*MSbot2(2) - 929*MStop2(2)) - 
     &      362*MStop2(2)**2)

	Sf3543 = -61*MStop2(1)**3 + 
     &   MSbot2(2)*(-51*MStop2(1)**2 + 
     &      MSbot2(2)*(68*MSbot2(2) + 39*MStop2(1)))

	Sf7884 = -121*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(68*MSbot2(2) - 929*MStop2(1)) - 
     &      362*MStop2(1)**2)

	Sf8151 = 85*MStop2(2)**3 + 
     &   MSbot2(2)*(34*MStop2(2)**2 + 
     &      MSbot2(2)*(67*MSbot2(2) + 34*MStop2(2)))

	Sf7942 = 85*MStop2(1)**3 + 
     &   MSbot2(2)*(34*MStop2(1)**2 + 
     &      MSbot2(2)*(67*MSbot2(2) + 34*MStop2(1)))

	Sf7110 = -48*MStop2(2)**3 + 
     &   MSbot2(2)*(-45*MStop2(2)**2 + 
     &      MSbot2(2)*(57*MSbot2(2) + 36*MStop2(2)))

	Sf4679 = -48*MStop2(1)**3 + 
     &   MSbot2(2)*(-45*MStop2(1)**2 + 
     &      MSbot2(2)*(57*MSbot2(2) + 36*MStop2(1)))

	Sf8145 = 37*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(55*MSbot2(2) - 27*MStop2(2)) - 
     &      9*MStop2(2)**2)

	Sf7936 = 37*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(55*MSbot2(2) - 27*MStop2(1)) - 
     &      9*MStop2(1)**2)

	Sf6884 = -57*MStop2(2)**3 + 
     &   MSbot2(2)*(-36*MStop2(2)**2 + 
     &      MSbot2(2)*(48*MSbot2(2) + 45*MStop2(2)))

	Sf6770 = -57*MStop2(1)**3 + 
     &   MSbot2(2)*(-36*MStop2(1)**2 + 
     &      MSbot2(2)*(48*MSbot2(2) + 45*MStop2(1)))

	Sf8382 = 35*MStop2(2)**3 + 
     &   MSbot2(2)*(-6*MStop2(2)**2 + 
     &      MSbot2(2)*(32*MSbot2(2) + 3*MStop2(2)))

	Sf5912 = -41*MStop2(2)**3 + 
     &   MSbot2(2)*(-15*MStop2(2)**2 + 
     &      MSbot2(2)*(32*MSbot2(2) + 3*MStop2(2)))

	Sf8358 = 35*MStop2(1)**3 + 
     &   MSbot2(2)*(-6*MStop2(1)**2 + 
     &      MSbot2(2)*(32*MSbot2(2) + 3*MStop2(1)))

	Sf3468 = -41*MStop2(1)**3 + 
     &   MSbot2(2)*(-15*MStop2(1)**2 + 
     &      MSbot2(2)*(32*MSbot2(2) + 3*MStop2(1)))

	Sf8214 = -10*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(29*MSbot2(2) - 24*MStop2(2)) + 
     &      5*MStop2(2)**2)

	Sf8202 = -10*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(29*MSbot2(2) - 24*MStop2(1)) + 
     &      5*MStop2(1)**2)

	Sf7099 = -16*MStop2(2)**3 + 
     &   MSbot2(2)*(-15*MStop2(2)**2 + 
     &      MSbot2(2)*(19*MSbot2(2) + 12*MStop2(2)))

	Sf4668 = -16*MStop2(1)**3 + 
     &   MSbot2(2)*(-15*MStop2(1)**2 + 
     &      MSbot2(2)*(19*MSbot2(2) + 12*MStop2(1)))

	Sf6874 = -19*MStop2(2)**3 + 
     &   MSbot2(2)*(-12*MStop2(2)**2 + 
     &      MSbot2(2)*(16*MSbot2(2) + 15*MStop2(2)))

	Sf6760 = -19*MStop2(1)**3 + 
     &   MSbot2(2)*(-12*MStop2(1)**2 + 
     &      MSbot2(2)*(16*MSbot2(2) + 15*MStop2(1)))

	Sf8179 = -23*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(14*MSbot2(2) - 7*MStop2(2)) + 
     &      16*MStop2(2)**2)

	Sf7970 = -23*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(14*MSbot2(2) - 7*MStop2(1)) + 
     &      16*MStop2(1)**2)

	Sf5344 = 574*MStop2(2)**3 + 
     &   MSbot2(2)*(255*MStop2(2)**2 + 
     &      MSbot2(2)*(10*MSbot2(2) + 189*MStop2(2)))

	Sf2933 = 574*MStop2(1)**3 + 
     &   MSbot2(2)*(255*MStop2(1)**2 + 
     &      MSbot2(2)*(10*MSbot2(2) + 189*MStop2(1)))

	Sf4914 = -5*MStop2(2)**3 + 
     &   MSbot2(2)*(21*MStop2(2)**2 + 
     &      MSbot2(2)*(9*MSbot2(2) + 7*MStop2(2)))

	Sf2343 = -5*MStop2(1)**3 + 
     &   MSbot2(2)*(21*MStop2(1)**2 + 
     &      MSbot2(2)*(9*MSbot2(2) + 7*MStop2(1)))

	Sf5457 = -3*MStop2(2)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(5*MSbot2(2) - 3*MStop2(2)) + 
     &      MStop2(2)**2)

	Sf3069 = -3*MStop2(1)**3 + 
     &   MSbot2(2)*(MSbot2(2)*(5*MSbot2(2) - 3*MStop2(1)) + 
     &      MStop2(1)**2)

	Sf9067 = MSbot2(2)**4 + 
     &   MSbot2(1)*(-4*MSbot2(2)**3 + 
     &      MSbot2(1)*(MSbot2(1)**2 + 6*MSbot2(2)**2))

	Sf9058 = MSbot2(2)**4 + 
     &   MSbot2(1)**2*(6*MSbot2(2)**2 + 
     &      MSbot2(1)*(-3*MSbot2(2) - MSbot2(3)))

	Sf8057 = 38197*MStop2(2)**3 + 
     &   MSbot2(1)*(38839*MStop2(2)**2 + 
     &      MSbot2(1)*(35443*MSbot2(1) + 37705*MStop2(2)))

	Sf7848 = 38197*MStop2(1)**3 + 
     &   MSbot2(1)*(38839*MStop2(1)**2 + 
     &      MSbot2(1)*(35443*MSbot2(1) + 37705*MStop2(1)))

	Sf8033 = 24494*MStop2(2)**3 + 
     &   MSbot2(1)*(23135*MStop2(2)**2 + 
     &      MSbot2(1)*(21578*MSbot2(1) + 21677*MStop2(2)))

	Sf7824 = 24494*MStop2(1)**3 + 
     &   MSbot2(1)*(23135*MStop2(1)**2 + 
     &      MSbot2(1)*(21578*MSbot2(1) + 21677*MStop2(1)))

	Sf8056 = 22179*MStop2(2)**3 + 
     &   MSbot2(1)*(12127*MStop2(2)**2 + 
     &      MSbot2(1)*(17805*MSbot2(1) + 10669*MStop2(2)))

	Sf7847 = 22179*MStop2(1)**3 + 
     &   MSbot2(1)*(12127*MStop2(1)**2 + 
     &      MSbot2(1)*(17805*MSbot2(1) + 10669*MStop2(1)))

	Sf8034 = 27099*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*
     &       (13977*MSbot2(1) - 11195*MStop2(2)) - 
     &      6821*MStop2(2)**2)

	Sf7825 = 27099*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*
     &       (13977*MSbot2(1) - 11195*MStop2(1)) - 
     &      6821*MStop2(1)**2)

	Sf8297 = -2751*MStop2(2)**3 + 
     &   MSbot2(1)*(-263*MStop2(2)**2 + 
     &      MSbot2(1)*(4637*MSbot2(1) + 689*MStop2(2)))

	Sf8273 = -2751*MStop2(1)**3 + 
     &   MSbot2(1)*(-263*MStop2(1)**2 + 
     &      MSbot2(1)*(4637*MSbot2(1) + 689*MStop2(1)))

	Sf8049 = 4230*MStop2(2)**3 + 
     &   MSbot2(1)*(2833*MStop2(2)**2 + 
     &      MSbot2(1)*(3474*MSbot2(1) + 2563*MStop2(2)))

	Sf7840 = 4230*MStop2(1)**3 + 
     &   MSbot2(1)*(2833*MStop2(1)**2 + 
     &      MSbot2(1)*(3474*MSbot2(1) + 2563*MStop2(1)))

	Sf8027 = 2179*MStop2(2)**3 + 
     &   MSbot2(1)*(-363*MStop2(2)**2 + 
     &      MSbot2(1)*(3313*MSbot2(1) + 15*MStop2(2)))

	Sf7818 = 2179*MStop2(1)**3 + 
     &   MSbot2(1)*(-363*MStop2(1)**2 + 
     &      MSbot2(1)*(3313*MSbot2(1) + 15*MStop2(1)))

	Sf8041 = 1915*MStop2(2)**3 + 
     &   MSbot2(1)*(2684*MStop2(2)**2 + 
     &      MSbot2(1)*(2482*MSbot2(1) + 983*MStop2(2)))

	Sf7832 = 1915*MStop2(1)**3 + 
     &   MSbot2(1)*(2684*MStop2(1)**2 + 
     &      MSbot2(1)*(2482*MSbot2(1) + 983*MStop2(1)))

	Sf8447 = -6365*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(1833*MSbot2(1) - 277*MStop2(2)) - 
     &      959*MStop2(2)**2)

	Sf8423 = -6365*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(1833*MSbot2(1) - 277*MStop2(1)) - 
     &      959*MStop2(1)**2)

	Sf8048 = 3289*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(1507*MSbot2(1) - 3001*MStop2(2)) - 
     &      2407*MStop2(2)**2)

	Sf7839 = 3289*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(1507*MSbot2(1) - 3001*MStop2(1)) - 
     &      2407*MStop2(1)**2)

	Sf8021 = 775*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(1207*MSbot2(1) - 693*MStop2(2)) - 
     &      801*MStop2(2)**2)

	Sf7801 = 775*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(1207*MSbot2(1) - 693*MStop2(1)) - 
     &      801*MStop2(1)**2)

	Sf8304 = 949*MStop2(2)**3 + 
     &   MSbot2(1)*(681*MStop2(2)**2 + 
     &      MSbot2(1)*(787*MSbot2(1) + 735*MStop2(2)))

	Sf8280 = 949*MStop2(1)**3 + 
     &   MSbot2(1)*(681*MStop2(1)**2 + 
     &      MSbot2(1)*(787*MSbot2(1) + 735*MStop2(1)))

	Sf8015 = 409*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(643*MSbot2(1) - 147*MStop2(2)) - 
     &      129*MStop2(2)**2)

	Sf7812 = 409*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(643*MSbot2(1) - 147*MStop2(1)) - 
     &      129*MStop2(1)**2)

	Sf5272 = 100*MStop2(2)**3 + 
     &   MSbot2(1)*(90*MStop2(2)**2 + 
     &      MSbot2(1)*(601*MSbot2(1) + 237*MStop2(2)))

	Sf2820 = 100*MStop2(1)**3 + 
     &   MSbot2(1)*(90*MStop2(1)**2 + 
     &      MSbot2(1)*(601*MSbot2(1) + 237*MStop2(1)))

	Sf5271 = 10*MStop2(2)**3 + 
     &   MSbot2(1)*(189*MStop2(2)**2 + 
     &      MSbot2(1)*(574*MSbot2(1) + 255*MStop2(2)))

	Sf2819 = 10*MStop2(1)**3 + 
     &   MSbot2(1)*(189*MStop2(1)**2 + 
     &      MSbot2(1)*(574*MSbot2(1) + 255*MStop2(1)))

	Sf8291 = -316*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(497*MSbot2(1) - 250*MStop2(2)) + 
     &      69*MStop2(2)**2)

	Sf8267 = -316*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(497*MSbot2(1) - 250*MStop2(1)) + 
     &      69*MStop2(1)**2)

	Sf8689 = 527*MStop2(2)**3 + 
     &   MSbot2(1)*(29*MStop2(2)**2 + 
     &      MSbot2(1)*(467*MSbot2(1) + 49*MStop2(2)))

	Sf8649 = 527*MStop2(1)**3 + 
     &   MSbot2(1)*(29*MStop2(1)**2 + 
     &      MSbot2(1)*(467*MSbot2(1) + 49*MStop2(1)))

	Sf8441 = -389*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(406*MSbot2(1) - 141*MStop2(2)) + 
     &      124*MStop2(2)**2)

	Sf8417 = -389*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(406*MSbot2(1) - 141*MStop2(1)) + 
     &      124*MStop2(1)**2)

	Sf8249 = -163*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(375*MSbot2(1) - 215*MStop2(2)) + 
     &      11*MStop2(2)**2)

	Sf8237 = -163*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(375*MSbot2(1) - 215*MStop2(1)) + 
     &      11*MStop2(1)**2)

	Sf8303 = 313*MStop2(2)**3 + 
     &   MSbot2(1)*(290*MStop2(2)**2 + 
     &      MSbot2(1)*(340*MSbot2(1) + 209*MStop2(2)))

	Sf8279 = 313*MStop2(1)**3 + 
     &   MSbot2(1)*(290*MStop2(1)**2 + 
     &      MSbot2(1)*(340*MSbot2(1) + 209*MStop2(1)))

	Sf8165 = 413*MStop2(2)**3 + 
     &   MSbot2(1)*(-23*MStop2(2)**2 + 
     &      MSbot2(1)*(323*MSbot2(1) + 7*MStop2(2)))

	Sf7956 = 413*MStop2(1)**3 + 
     &   MSbot2(1)*(-23*MStop2(1)**2 + 
     &      MSbot2(1)*(323*MSbot2(1) + 7*MStop2(1)))

	Sf5277 = 142*MStop2(2)**3 + 
     &   MSbot2(1)*(291*MStop2(2)**2 + 
     &      MSbot2(1)*(303*MSbot2(1) + 288*MStop2(2)))

	Sf2825 = 142*MStop2(1)**3 + 
     &   MSbot2(1)*(291*MStop2(1)**2 + 
     &      MSbot2(1)*(303*MSbot2(1) + 288*MStop2(1)))

	Sf8040 = -199*MStop2(2)**3 + 
     &   MSbot2(1)*(1645*MStop2(2)**2 + 
     &      MSbot2(1)*(287*MSbot2(1) + 1483*MStop2(2)))

	Sf7831 = -199*MStop2(1)**3 + 
     &   MSbot2(1)*(1645*MStop2(1)**2 + 
     &      MSbot2(1)*(287*MSbot2(1) + 1483*MStop2(1)))

	Sf6598 = -230*MStop2(2)**3 + 
     &   MSbot2(1)*(-30*MStop2(2)**2 + 
     &      MSbot2(1)*(280*MSbot2(1) + 84*MStop2(2)))

	Sf6472 = -230*MStop2(1)**3 + 
     &   MSbot2(1)*(-30*MStop2(1)**2 + 
     &      MSbot2(1)*(280*MSbot2(1) + 84*MStop2(1)))

	Sf6603 = 173*MStop2(2)**3 + 
     &   MSbot2(1)*(20*MStop2(2)**2 + 
     &      MSbot2(1)*(246*MSbot2(1) + 73*MStop2(2)))

	Sf6477 = 173*MStop2(1)**3 + 
     &   MSbot2(1)*(20*MStop2(1)**2 + 
     &      MSbot2(1)*(246*MSbot2(1) + 73*MStop2(1)))

	Sf7039 = -280*MStop2(2)**3 + 
     &   MSbot2(1)*(-84*MStop2(2)**2 + 
     &      MSbot2(1)*(230*MSbot2(1) + 30*MStop2(2)))

	Sf4532 = -280*MStop2(1)**3 + 
     &   MSbot2(1)*(-84*MStop2(1)**2 + 
     &      MSbot2(1)*(230*MSbot2(1) + 30*MStop2(1)))

	Sf5242 = 257*MStop2(2)**3 + 
     &   MSbot2(1)*(296*MStop2(2)**2 + 
     &      MSbot2(1)*(210*MSbot2(1) + 261*MStop2(2)))

	Sf2773 = 257*MStop2(1)**3 + 
     &   MSbot2(1)*(296*MStop2(1)**2 + 
     &      MSbot2(1)*(210*MSbot2(1) + 261*MStop2(1)))

	Sf7045 = 234*MStop2(2)**3 + 
     &   MSbot2(1)*(79*MStop2(2)**2 + 
     &      MSbot2(1)*(179*MSbot2(1) + 20*MStop2(2)))

	Sf4538 = 234*MStop2(1)**3 + 
     &   MSbot2(1)*(79*MStop2(1)**2 + 
     &      MSbot2(1)*(179*MSbot2(1) + 20*MStop2(1)))

	Sf6597 = 160*MStop2(2)**3 + 
     &   MSbot2(1)*(99*MStop2(2)**2 + 
     &      MSbot2(1)*(160*MSbot2(1) + 105*MStop2(2)))

	Sf7040 = 160*MStop2(2)**3 + 
     &   MSbot2(1)*(105*MStop2(2)**2 + 
     &      MSbot2(1)*(160*MSbot2(1) + 99*MStop2(2)))

	Sf6471 = 160*MStop2(1)**3 + 
     &   MSbot2(1)*(99*MStop2(1)**2 + 
     &      MSbot2(1)*(160*MSbot2(1) + 105*MStop2(1)))

	Sf4533 = 160*MStop2(1)**3 + 
     &   MSbot2(1)*(105*MStop2(1)**2 + 
     &      MSbot2(1)*(160*MSbot2(1) + 99*MStop2(1)))

	Sf8171 = -555*MStop2(2)**3 + 
     &   MSbot2(1)*(341*MStop2(2)**2 + 
     &      MSbot2(1)*(145*MSbot2(1) + 61*MStop2(2)))

	Sf7962 = -555*MStop2(1)**3 + 
     &   MSbot2(1)*(341*MStop2(1)**2 + 
     &      MSbot2(1)*(145*MSbot2(1) + 61*MStop2(1)))

	Sf8009 = -242*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(136*MSbot2(1) - 1858*MStop2(2)) - 
     &      724*MStop2(2)**2)

	Sf5237 = 601*MStop2(2)**3 + 
     &   MSbot2(1)*(237*MStop2(2)**2 + 
     &      MSbot2(1)*(100*MSbot2(1) + 90*MStop2(2)))

	Sf2768 = 601*MStop2(1)**3 + 
     &   MSbot2(1)*(237*MStop2(1)**2 + 
     &      MSbot2(1)*(100*MSbot2(1) + 90*MStop2(1)))

	Sf8072 = 189*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(99*MSbot2(1) - 125*MStop2(2)) - 
     &      143*MStop2(2)**2)

	Sf7863 = 189*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(99*MSbot2(1) - 125*MStop2(1)) - 
     &      143*MStop2(1)**2)

	Sf8333 = -30*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(87*MSbot2(1) - 72*MStop2(1)) + 
     &      15*MStop2(1)**2)

	Sf5542 = -61*MStop2(2)**3 + 
     &   MSbot2(1)*(-51*MStop2(2)**2 + 
     &      MSbot2(1)*(68*MSbot2(1) + 39*MStop2(2)))

	Sf3181 = -61*MStop2(1)**3 + 
     &   MSbot2(1)*(-51*MStop2(1)**2 + 
     &      MSbot2(1)*(68*MSbot2(1) + 39*MStop2(1)))

	Sf7806 = -121*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(68*MSbot2(1) - 929*MStop2(1)) - 
     &      362*MStop2(1)**2)

	Sf8073 = 85*MStop2(2)**3 + 
     &   MSbot2(1)*(34*MStop2(2)**2 + 
     &      MSbot2(1)*(67*MSbot2(1) + 34*MStop2(2)))

	Sf7864 = 85*MStop2(1)**3 + 
     &   MSbot2(1)*(34*MStop2(1)**2 + 
     &      MSbot2(1)*(67*MSbot2(1) + 34*MStop2(1)))

	Sf7072 = -48*MStop2(2)**3 + 
     &   MSbot2(1)*(-45*MStop2(2)**2 + 
     &      MSbot2(1)*(57*MSbot2(1) + 36*MStop2(2)))

	Sf4574 = -48*MStop2(1)**3 + 
     &   MSbot2(1)*(-45*MStop2(1)**2 + 
     &      MSbot2(1)*(57*MSbot2(1) + 36*MStop2(1)))

	Sf8067 = 37*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(55*MSbot2(1) - 27*MStop2(2)) - 
     &      9*MStop2(2)**2)

	Sf7858 = 37*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(55*MSbot2(1) - 27*MStop2(1)) - 
     &      9*MStop2(1)**2)

	Sf6633 = -57*MStop2(2)**3 + 
     &   MSbot2(1)*(-36*MStop2(2)**2 + 
     &      MSbot2(1)*(48*MSbot2(1) + 45*MStop2(2)))

	Sf6515 = -57*MStop2(1)**3 + 
     &   MSbot2(1)*(-36*MStop2(1)**2 + 
     &      MSbot2(1)*(48*MSbot2(1) + 45*MStop2(1)))

	Sf8159 = -69*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(42*MSbot2(1) - 21*MStop2(2)) + 
     &      48*MStop2(2)**2)

	Sf8286 = 35*MStop2(2)**3 + 
     &   MSbot2(1)*(-6*MStop2(2)**2 + 
     &      MSbot2(1)*(32*MSbot2(1) + 3*MStop2(2)))

	Sf5478 = -41*MStop2(2)**3 + 
     &   MSbot2(1)*(-15*MStop2(2)**2 + 
     &      MSbot2(1)*(32*MSbot2(1) + 3*MStop2(2)))

	Sf8262 = 35*MStop2(1)**3 + 
     &   MSbot2(1)*(-6*MStop2(1)**2 + 
     &      MSbot2(1)*(32*MSbot2(1) + 3*MStop2(1)))

	Sf3081 = -41*MStop2(1)**3 + 
     &   MSbot2(1)*(-15*MStop2(1)**2 + 
     &      MSbot2(1)*(32*MSbot2(1) + 3*MStop2(1)))

	Sf8243 = -10*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(29*MSbot2(1) - 24*MStop2(2)) + 
     &      5*MStop2(2)**2)

	Sf8231 = -10*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(29*MSbot2(1) - 24*MStop2(1)) + 
     &      5*MStop2(1)**2)

	Sf7059 = -16*MStop2(2)**3 + 
     &   MSbot2(1)*(-15*MStop2(2)**2 + 
     &      MSbot2(1)*(19*MSbot2(1) + 12*MStop2(2)))

	Sf4563 = -16*MStop2(1)**3 + 
     &   MSbot2(1)*(-15*MStop2(1)**2 + 
     &      MSbot2(1)*(19*MSbot2(1) + 12*MStop2(1)))

	Sf6623 = -19*MStop2(2)**3 + 
     &   MSbot2(1)*(-12*MStop2(2)**2 + 
     &      MSbot2(1)*(16*MSbot2(1) + 15*MStop2(2)))

	Sf6505 = -19*MStop2(1)**3 + 
     &   MSbot2(1)*(-12*MStop2(1)**2 + 
     &      MSbot2(1)*(16*MSbot2(1) + 15*MStop2(1)))

	Sf7950 = -23*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(14*MSbot2(1) - 7*MStop2(1)) + 
     &      16*MStop2(1)**2)

	Sf5236 = 574*MStop2(2)**3 + 
     &   MSbot2(1)*(255*MStop2(2)**2 + 
     &      MSbot2(1)*(10*MSbot2(1) + 189*MStop2(2)))

	Sf2767 = 574*MStop2(1)**3 + 
     &   MSbot2(1)*(255*MStop2(1)**2 + 
     &      MSbot2(1)*(10*MSbot2(1) + 189*MStop2(1)))

	Sf4823 = -5*MStop2(2)**3 + 
     &   MSbot2(1)*(21*MStop2(2)**2 + 
     &      MSbot2(1)*(9*MSbot2(1) + 7*MStop2(2)))

	Sf2330 = -5*MStop2(1)**3 + 
     &   MSbot2(1)*(21*MStop2(1)**2 + 
     &      MSbot2(1)*(9*MSbot2(1) + 7*MStop2(1)))

	Sf5450 = -3*MStop2(2)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(5*MSbot2(1) - 3*MStop2(2)) + 
     &      MStop2(2)**2)

	Sf3061 = -3*MStop2(1)**3 + 
     &   MSbot2(1)*(MSbot2(1)*(5*MSbot2(1) - 3*MStop2(1)) + 
     &      MStop2(1)**2)

	Sf9806 = -528*MStau2(2)**3 + 
     &   MZ2*(MZ2*(605*MZ2 - 1521*MStau2(2)) - 562*MStau2(2)**2)

	Sf9672 = -528*MStau2(1)**3 + 
     &   MZ2*(MZ2*(605*MZ2 - 1521*MStau2(1)) - 562*MStau2(1)**2)

	Sf10018 = -264*MStau2(2)**3 + 
     &   MZ2*(MZ2*(116*MZ2 - 3333*MStau2(2)) - 1768*MStau2(2)**2)

	Sf9921 = -264*MStau2(1)**3 + 
     &   MZ2*(MZ2*(116*MZ2 - 3333*MStau2(1)) - 1768*MStau2(1)**2)

	Sf9827 = -240*MStau2(2)**3 + 
     &   MZ2*(MZ2*(95*MZ2 - 351*MStau2(2)) - 190*MStau2(2)**2)

	Sf9642 = -240*MStau2(1)**3 + 
     &   MZ2*(MZ2*(95*MZ2 - 351*MStau2(1)) - 190*MStau2(1)**2)

	Sf10039 = 160*MStau2(2)**3 + 
     &   MZ2*(MZ2*(68*MZ2 - 123*MStau2(2)) - 86*MStau2(2)**2)

	Sf9958 = 160*MStau2(1)**3 + 
     &   MZ2*(MZ2*(68*MZ2 - 123*MStau2(1)) - 86*MStau2(1)**2)

	Sf9983 = 64*MStau2(2)**3 + 
     &   MZ2*(100*MStau2(2)**2 + MZ2*(17*MZ2 + 45*MStau2(2)))

	Sf9701 = 24*MStau2(2)**3 + 
     &   MZ2*(68*MStau2(2)**2 + MZ2*(17*MZ2 + 24*MStau2(2)))

	Sf9885 = 64*MStau2(1)**3 + 
     &   MZ2*(100*MStau2(1)**2 + MZ2*(17*MZ2 + 45*MStau2(1)))

	Sf9408 = 24*MStau2(1)**3 + 
     &   MZ2*(68*MStau2(1)**2 + MZ2*(17*MZ2 + 24*MStau2(1)))

	Sf9990 = -8*MStau2(2)**3 + 
     &   MZ2*(MZ2*(4*MZ2 - 41*MStau2(2)) + 2*MStau2(2)**2)

	Sf9875 = -8*MStau2(1)**3 + 
     &   MZ2*(MZ2*(4*MZ2 - 41*MStau2(1)) + 2*MStau2(1)**2)

	Sf9992 = -24*MStau2(2)**3 + 
     &   MZ2*(MZ2*(MZ2 - 4*MStau2(2)) + 37*MStau2(2)**2)

	Sf9700 = -2*MStau2(2)**3 + 
     &   MZ2*(MZ2*(MZ2 - 15*MStau2(2)) + 5*MStau2(2)**2)

	Sf9873 = -24*MStau2(1)**3 + 
     &   MZ2*(MZ2*(MZ2 - 4*MStau2(1)) + 37*MStau2(1)**2)

	Sf9405 = -2*MStau2(1)**3 + 
     &   MZ2*(MZ2*(MZ2 - 15*MStau2(1)) + 5*MStau2(1)**2)

	Sf8490 = 99*MZ2**4 + 
     &   MW2*(-112*MZ2**3 + MW2*(64*MW2**2 - 60*MZ2**2))

	Sf3568 = 70*MSbot2(2)**3 + 
     &   MW2*(90*MSbot2(2)**2 + MW2*(203*MW2 + 108*MSbot2(2)))

	Sf3216 = 70*MSbot2(1)**3 + 
     &   MW2*(90*MSbot2(1)**2 + MW2*(203*MW2 + 108*MSbot2(1)))

	Sf3976 = 240*MSbot2(2)**3 + 
     &   MW2*(192*MSbot2(2)**2 + MW2*(160*MW2 + 191*MSbot2(2)))

	Sf3879 = 240*MSbot2(1)**3 + 
     &   MW2*(192*MSbot2(1)**2 + MW2*(160*MW2 + 191*MSbot2(1)))

	Sf8735 = -13*MZ2**3 + MW2*(MW2*(64*MW2 - 144*MZ2) + 84*MZ2**2)

	Sf8491 = 13*MZ2**3 + MW2*(MW2*(32*MW2 - 48*MZ2) - 6*MZ2**2)

	Sf7764 = -13*MZ2**3 + MW2*(MW2*(32*MW2 - 80*MZ2) + 58*MZ2**2)

	Sf7783 = -11*MZ2**3 + MW2*(MW2*(16*MW2 - 40*MZ2) + 35*MZ2**2)

	Sf7784 = -13*MZ2**3 + MW2*(MW2*(16*MW2 - 48*MZ2) + 45*MZ2**2)

	Sf3663 = 21*MSbot2(2)**3 + 
     &   MW2*(17*MSbot2(2)**2 + MW2*(8*MW2 + 14*MSbot2(2)))

	Sf3383 = 21*MSbot2(1)**3 + 
     &   MW2*(17*MSbot2(1)**2 + MW2*(8*MW2 + 14*MSbot2(1)))

	Sf8611 = -7*MZ2**3 + MW2*(MW2*(8*MW2 - 24*MZ2) + 23*MZ2**2)

	Sf9716 = -208*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*
     &       (1210*MHiggs2(2) - 1159*MStau2(2)) - 1790*MStau2(2)**2
     &      )

	Sf9423 = -208*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*
     &       (1210*MHiggs2(2) - 1159*MStau2(1)) - 1790*MStau2(1)**2
     &      )

	Sf9799 = 164*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*
     &       (399*MHiggs2(2) - 2016*MStau2(2)) + 636*MStau2(2)**2)

	Sf9659 = 164*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*
     &       (399*MHiggs2(2) - 2016*MStau2(1)) + 636*MStau2(1)**2)

	Sf9451 = 32*MZ2**3 + 
     &   MHiggs2(2)*(4489*MZ2**2 + 
     &      MHiggs2(2)*(-2809*MZ2 + 292*MHiggs2(2)))

	Sf9772 = -160*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(226*MHiggs2(2) - 157*MStau2(2)) - 
     &      434*MStau2(2)**2)

	Sf9536 = -160*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(226*MHiggs2(2) - 157*MStau2(1)) - 
     &      434*MStau2(1)**2)

	Sf9443 = -10*MZ2**3 + 
     &   MHiggs2(2)*(-1661*MZ2**2 + 
     &      MHiggs2(2)*(437*MZ2 + 226*MHiggs2(2)))

	Sf9714 = -152*MStau2(2)**3 + 
     &   MHiggs2(2)*(-686*MStau2(2)**2 + 
     &      MHiggs2(2)*(208*MHiggs2(2) + 1303*MStau2(2)))

	Sf9424 = -152*MStau2(1)**3 + 
     &   MHiggs2(2)*(-686*MStau2(1)**2 + 
     &      MHiggs2(2)*(208*MHiggs2(2) + 1303*MStau2(1)))

	Sf9564 = 287*MZ2**3 + 
     &   MHiggs2(2)*(1031*MZ2**2 + 
     &      MHiggs2(2)*(-845*MZ2 + 165*MHiggs2(2)))

	Sf10010 = -216*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*
     &       (116*MHiggs2(2) - 3813*MStau2(2)) - 112*MStau2(2)**2)

	Sf9911 = -216*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*
     &       (116*MHiggs2(2) - 3813*MStau2(1)) - 112*MStau2(1)**2)

	Sf9820 = 140*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(57*MHiggs2(2) - 360*MStau2(2)) + 
     &      240*MStau2(2)**2)

	Sf9629 = 140*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(57*MHiggs2(2) - 360*MStau2(1)) + 
     &      240*MStau2(1)**2)

	Sf9770 = -200*MStau2(2)**3 + 
     &   MHiggs2(2)*(-266*MStau2(2)**2 + 
     &      MHiggs2(2)*(52*MHiggs2(2) + 229*MStau2(2)))

	Sf9537 = -200*MStau2(1)**3 + 
     &   MHiggs2(2)*(-266*MStau2(1)**2 + 
     &      MHiggs2(2)*(52*MHiggs2(2) + 229*MStau2(1)))

	Sf9509 = -80*MZ2**3 + 
     &   MHiggs2(2)*(119*MZ2**2 + 
     &      MHiggs2(2)*(-119*MZ2 + 44*MHiggs2(2)))

	Sf9516 = -66*MZ2**3 + 
     &   MHiggs2(2)*(-78*MZ2**2 + 
     &      MHiggs2(2)*(62*MZ2 + 42*MHiggs2(2)))

	Sf9557 = -47*MZ2**3 + 
     &   MHiggs2(2)*(616*MZ2**2 + 
     &      MHiggs2(2)*(-351*MZ2 + 30*MHiggs2(2)))

	Sf9457 = -246*MZ2**3 + 
     &   MHiggs2(2)*(2619*MZ2**2 + 
     &      MHiggs2(2)*(-3715*MZ2 + 30*MHiggs2(2)))

	Sf9501 = -50*MZ2**3 + 
     &   MHiggs2(2)*(-31*MZ2**2 + 
     &      MHiggs2(2)*(55*MZ2 + 26*MHiggs2(2)))

	Sf10046 = 504*MStau2(2)**3 + 
     &   MHiggs2(2)*(116*MStau2(2)**2 + 
     &      MHiggs2(2)*(20*MHiggs2(2) + 1101*MStau2(2)))

	Sf9969 = 504*MStau2(1)**3 + 
     &   MHiggs2(2)*(116*MStau2(1)**2 + 
     &      MHiggs2(2)*(20*MHiggs2(2) + 1101*MStau2(1)))

	Sf9996 = 120*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(16*MHiggs2(2) - 63*MStau2(2)) + 
     &      34*MStau2(2)**2)

	Sf9895 = 120*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(16*MHiggs2(2) - 63*MStau2(1)) + 
     &      34*MStau2(1)**2)

	Sf9435 = 38*MZ2**3 + 
     &   MHiggs2(2)*(-635*MZ2**2 + 
     &      MHiggs2(2)*(-689*MZ2 + 16*MHiggs2(2)))

	Sf9746 = 14*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(13*MHiggs2(2) - 13*MStau2(2)) - 
     &      20*MStau2(2)**2)

	Sf9483 = 14*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(13*MHiggs2(2) - 13*MStau2(1)) - 
     &      20*MStau2(1)**2)

	Sf9745 = 16*MStau2(2)**3 + 
     &   MHiggs2(2)*(-2*MStau2(2)**2 + 
     &      MHiggs2(2)*(7*MHiggs2(2) + 19*MStau2(2)))

	Sf9484 = 16*MStau2(1)**3 + 
     &   MHiggs2(2)*(-2*MStau2(1)**2 + 
     &      MHiggs2(2)*(7*MHiggs2(2) + 19*MStau2(1)))

	Sf9689 = 68*MStau2(2)**3 + 
     &   MHiggs2(2)*(87*MStau2(2)**2 + 
     &      MHiggs2(2)*(6*MHiggs2(2) + 24*MStau2(2)))

	Sf9392 = 68*MStau2(1)**3 + 
     &   MHiggs2(2)*(87*MStau2(1)**2 + 
     &      MHiggs2(2)*(6*MHiggs2(2) + 24*MStau2(1)))

	Sf9551 = 58*MZ2**3 + 
     &   MHiggs2(2)*(7*MZ2**2 + MHiggs2(2)*(43*MZ2 + 4*MHiggs2(2)))

	Sf9494 = -7*MZ2**3 + 
     &   MHiggs2(2)*(MZ2**2 + MHiggs2(2)*(MZ2 + 4*MHiggs2(2)))

	Sf9053 = -4*MZ2**3 + 
     &   MHiggs2(2)*(-7*MZ2**2 + MHiggs2(2)*(7*MZ2 + 3*MHiggs2(2)))

	Sf9559 = 84*MZ2**3 + 
     &   MHiggs2(2)*(67*MZ2**2 + 
     &      MHiggs2(2)*(42*MZ2 + 2*MHiggs2(2)))

	Sf9572 = 86*MZ2**3 + 
     &   MHiggs2(2)*(139*MZ2**2 + MHiggs2(2)*(-MZ2 + 2*MHiggs2(2)))

	Sf9577 = 196*MZ2**3 + 
     &   MHiggs2(2)*(-109*MZ2**2 + 
     &      MHiggs2(2)*(-185*MZ2 + 2*MHiggs2(2)))

	Sf9744 = 8*MStau2(2)**3 + 
     &   MHiggs2(2)*(-9*MStau2(2)**2 + 
     &      MHiggs2(2)*(MHiggs2(2) + 6*MStau2(2)))

	Sf9688 = -4*MStau2(2)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(MHiggs2(2) - 18*MStau2(2)) + 
     &      5*MStau2(2)**2)

	Sf9482 = 8*MStau2(1)**3 + 
     &   MHiggs2(2)*(-9*MStau2(1)**2 + 
     &      MHiggs2(2)*(MHiggs2(2) + 6*MStau2(1)))

	Sf9389 = -4*MStau2(1)**3 + 
     &   MHiggs2(2)*(MHiggs2(2)*(MHiggs2(2) - 18*MStau2(1)) + 
     &      5*MStau2(1)**2)

	Sf9549 = 119*MZ2**3 + 
     &   MHiggs2(2)*(MHiggs2(2)**2 + 
     &      MZ2*(203*MZ2 + 211*MHiggs2(2)))

	Sf9566 = -69*MZ2**3 + 
     &   MHiggs2(2)*(MHiggs2(2)**2 + MZ2*(-66*MZ2 + 26*MHiggs2(2)))

	Sf8978 = -2*MZ2**3 + 
     &   MHiggs2(2)*(MHiggs2(2)**2 + MZ2*(-MZ2 + MHiggs2(2)))

	Sf9183 = -4*MZ2**3 + 
     &   MHiggs2(1)*(-7*MZ2**2 + MHiggs2(1)*(7*MZ2 + 3*MHiggs2(1)))

	Sf9168 = -2*MZ2**3 + 
     &   MHiggs2(1)*(MHiggs2(1)**2 + MZ2*(-MZ2 + MHiggs2(1)))

	Sf1734 = 2*ML2*SA*UStau21 - CB*SAB*Ss443*UStau21 + 
     &   ML*SA*Re(Atau)*UStau22

	Sf8544 = Sf8543*MStau2(1) + 
     &   2*UStau11*(MSneu2(1) + MStau2(2)*(1 - UStau22**2))

	Sf946 = 2*ML2*SA*UStau11 - CB*SAB*Ss443*UStau11 + 
     &   ML*SA*Re(Atau)*UStau12

	Sf2310 = (reMUE*TB + Re(Atop))*USbot21 + 
     &   MBot*Ss199*USbot22

	Sf8866 = 78 - 241*(-1 + USbot22**2)

	Sf8636 = 52 - 7*USbot22**2

	Sf8560 = 135 + 14*(-1 + USbot22**2)

	Sf8588 = 11 - 2*USbot22**2

	Sf5228 = 3*MBot*MSbot2(3)*USbot21*USbot22 + 
     &   ML*MStau2(2)*UStau21*UStau22

	Sf2259 = (reMUE*TB + Re(Atop))*USbot11 + 
     &   MBot*Ss199*USbot12

	Sf8530 = MStop2(2)**2 + MStop2(1)*(MStop2(1) + 12*MStop2(2))

	Sf9417 = -10*MStau2(2)**2 + 
     &   MStau2(1)*(793*MStau2(1) + 879*MStau2(2))

	Sf9479 = -2*MStau2(2)**2 + 
     &   MStau2(1)*(17*MStau2(1) + 15*MStau2(2))

	Sf9708 = -793*MStau2(2)**2 + 
     &   MStau2(1)*(-869*MStau2(2) - 10*MStau2(3))

	Sf9739 = -17*MStau2(2)**2 + 
     &   MStau2(1)*(-13*MStau2(2) - 2*MStau2(3))

	Sf7711 = MStau2(2)**2 + MStau2(1)*(MStau2(1) + 12*MStau2(2))

	Sf1328 = MSneu2(1)**2*(28*MSneu2(1) - 33*MStau2(2)) + 
     &   5*MStau2(2)**3

	Sf391 = MSneu2(1)**2*(28*MSneu2(1) - 33*MStau2(1)) + 
     &   5*MStau2(1)**3

	Sf1959 = -8*MStau2(2)**3 + 
     &   MSneu2(1)**2*(25*MSneu2(1) + 3*MStau2(2))

	Sf1897 = -8*MStau2(1)**3 + 
     &   MSneu2(1)**2*(25*MSneu2(1) + 3*MStau2(1))

	Sf1356 = 28*MStau2(2)**3 + 
     &   MSneu2(1)*(5*MSneu2(1)**2 - 33*MStau2(2)**2)

	Sf433 = 28*MStau2(1)**3 + 
     &   MSneu2(1)*(5*MSneu2(1)**2 - 33*MStau2(1)**2)

	Sf7594 = 10981*MStau2(2)**2 + 
     &   MSneu2(1)*(10009*MSneu2(1) + 14368*MStau2(2))

	Sf7536 = 10981*MStau2(1)**2 + 
     &   MSneu2(1)*(10009*MSneu2(1) + 14368*MStau2(1))

	Sf7624 = 2913*MStau2(2)**2 + 
     &   MSneu2(1)*(2697*MSneu2(1) + 3764*MStau2(2))

	Sf7566 = 2913*MStau2(1)**2 + 
     &   MSneu2(1)*(2697*MSneu2(1) + 3764*MStau2(1))

	Sf7602 = 2323*MStau2(2)**2 + 
     &   MSneu2(1)*(1999*MSneu2(1) + 2776*MStau2(2))

	Sf7544 = 2323*MStau2(1)**2 + 
     &   MSneu2(1)*(1999*MSneu2(1) + 2776*MStau2(1))

	Sf1707 = MSneu2(1)*(1907*MSneu2(1) - 80*MStau2(2)) + 
     &   1907*MStau2(2)**2

	Sf908 = MSneu2(1)*(1907*MSneu2(1) - 80*MStau2(1)) + 
     &   1907*MStau2(1)**2

	Sf7623 = 2887*MStau2(2)**2 + 
     &   MSneu2(1)*(1375*MSneu2(1) + 1420*MStau2(2))

	end

