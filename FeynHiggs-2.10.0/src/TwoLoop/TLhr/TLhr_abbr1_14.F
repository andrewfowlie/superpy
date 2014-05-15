* TLhr_abbr1_14.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_14
	implicit none

#include "TLhr.h"


        Cd(7404) = -16*Cd(5439) + 
     &   3*MTy2*(3*MStop2(1)**3 + 
     &      MTy2*(MTy2*(20*MTy2 - 8*MStop2(1)) + 6*MStop2(1)**2))-
     &     Cd(7401)*MStop2(2) + 3*Cd(7402)*MStop2(2)**2 + 
     &   Cd(7403)*MStop2(2)**3

        Cd(7887) = 30*MTy2**4 + 2*Cd(7886) - 
     &   6*MTy2**3*MStop2(1) + 13*MTy2**2*MStop2(1)**2 + 
     &   5*MTy2*MStop2(1)**3 + 
     &   MStop2(2)*(-Cd(7883) + Cd(3770)*Cd(6005)*MStop2(2) + 
     &      Cd(7884)*MStop2(2)**2)

        Cd(7918) = 5*Cd(5439) + 
     &   2*(MTy2*(MTy2 - MStop2(1))*
     &       (5*MStop2(1)**2 + MTy2*(5*MTy2 + 7*MStop2(1))) + 
     &      MStop2(2)*(Cd(7916) - 
     &         MStop2(2)*(Cd(7917) + 5*Cd(45)*MStop2(2))))

        Cd(8003) = -2*MGl2*Cd(7997) + 
     &   MGl2**2*(Cd(7988) + 
     &      MGl2*(2*Cd(7975) + MGl2*(6*MGl2*Cd(7978) - Cd(8002))))+
     &     (MGl2**6*(4*MGl2 - 2*Cd(5861)) + 
     &      Cd(44)**2*Cd(48)**2*Cd(7968))*UStop2(1,1)

        Cd(8141) = 3*Cd(7886) - 
     &   2*MStop2(2)*(Cd(8139) + Cd(8140)*MStop2(2)) + 
     &   3*MTy2*(MTy2*(2*MStop2(1)**2 + 
     &         MTy2*(5*MTy2 + 10*MStop2(1))) - 6*MStop2(2)**3 + 
     &      6*MStop2(2)**2*MStop2(3) - 
     &      MStop2(3)**2*(MStop2(3) + 3*MStop2(4)))

        Cd(10160) = Cd(9948)*
     &    (Cd(10101) - 3*Cd(9951)*UStop2(1,1)**2) + 
     &   UStop2(1,1)*(Cd(10101) + 
     &      UStop2(1,1)*(2*Cd(9953) - Cd(9951)*UStop2(1,1))) + 
     &   4*Cd(9953)*UStop2(1,1)*UStop2(1,2) + 
     &   (2*Cd(9953) - 3*Cd(9951)*UStop2(1,1))*UStop2(1,2)**2 - 
     &   Cd(9951)*UStop2(1,2)**3

	Cd(10221) = Opt(5) + Opt(6)

	Cd(10263) = Opt(5) - Opt(6)

        Cd(10881) = MTy2*Cd(10751) + 2*M_3*MTy*Cd(10880) + 
     &   2*M_3C*MTy*(-3*MUE**2*UCStop(1,3)**3 + 
     &      2*MUE2*UCStop(1,3)*UStop2(1,1)**2 + 
     &      MUEC**2*UCStopC(1,3)*UStop2(1,1)**2) + 
     &   MUEC**2*MStop2(1)*UCStopC(1,3)**2*(1 + UStop2(1,2)) + 
     &   MUE**2*UCStop(1,3)**2*(MStop2(1) - MStop2(3)*UStop2(1,2))

        Cd(13705) = MUEC*TB*Cd(12700)*Cd(13697) + 
     &   2*MGl2*MGlmT2*(MTy2 - MSbot2(1))*MSbot2(1) + 
     &   At*(AtC*Cd(13697) + 2*MGl2*MGlmT2*MUE*TB*MSbot2(1)) - 
     &   MStop2(2)*(Cd(13700) + 2*Cd(13703)*MStop2(2) + 
     &      MStop2(2)**2*(-Cd(13704) + 2*MStop2(2)))

        Cd(13966) = MUEC*TB*Cd(12700)*Cd(13959) + 
     &   2*MGlmT2*MTy2*(MTy2 - MSbot2(1))*MSbot2(1) + 
     &   At*(AtC*Cd(13959) + 2*MGlmT2*MTy2*MUE*TB*MSbot2(1)) - 
     &   Cd(13961)*MStop2(2) - 2*Cd(13964)*MStop2(2)**2 - 
     &   Cd(13965)*MStop2(2)**3 + 2*MStop2(2)**4

        Cd(13978) = MStop2(2)*
     &    (2*MGl2*Cd(13976) + MGl2**2*Cd(13977) + 
     &      2*MSbot2(1)*(Cd(13974) - Cr(53)*MSbot2(1))*MStop2(2) - 
     &      2*Cd(13971)*MStop2(2)**3) - 
     &   MStop2(2)**2*(At*MUE*TB*Cd(13734)*Cr(56) + 
     &      Cd(13972)*MStop2(2))*UStop2(1,1) + 
     &   (Cd(13966) - AtC*MUEC*TB*Cd(13734)*MStop2(2)**2)*
     &    UStop2(1,1)**2

        Cd(14055) = -(MGl2**2*TB2*Cd(14054)) + 
     &   MTy*MUE*TB2*Cd(9948)*(MTy2**2 + MStop2(1)**2)*
     &    UCStop(1,3) + (-(At*MUE*TB2*Cd(14052)) + 
     &      MGl2*TB*Cd(14053) + 
     &      MTy2*TB*(MUE2*TB2 - MSq2Diff(1,1))*(MTy2 - MStop2(1)))*
     &    UStop2(1,2)**2

        Cd(1920) = -2*(-4*MTy2**3 + 
     &      MGl2*(MGl2**2 + 6*MTy2**2)) + 8*MGl2**2*Cd(49) + 
     &   (MGl2 - 4*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(-2*Cd(1918) + 
     &      MStop2(2)*(-Cd(1919) + 6*MStop2(2)))

        Cd(297) = MGl2*(Cd(289) + MGl2*Cd(291)) - 
     &   MTy2**2*Cd(292) + MTy2*Cd(294) - MTy*Cd(296) + 
     &   (MTy2**3 - Cd(44)*Cd(48)**2)*UStop2(1,1) + 
     &   MStop2(1)*MStop2(2)*(-MStop2(1) + MStop2(2)*UStop2(1,2))

        Cd(3228) = -(MGl2**3*Cd(1718)) + MGl2**2*Cd(3225) + 
     &   MTy2*(-(MGl2*Cd(3223)) + Cd(3227) + 
     &      MTy2*(8*MTy2**2 - 3*Cd(637)*MStop2(2))) + 
     &   (-5*MTy2**3 - 3*MGl2*MStop2(2)**2 + MStop2(2)**3)*
     &    MStop2(3)

        Cd(4929) = 30*MTy2**4 - Cd(5439) - 
     &   12*MTy2**3*MStop2(1) + 17*MTy2**2*MStop2(1)**2 + 
     &   9*MTy2*MStop2(1)**3 + 
     &   MStop2(2)*(-Cd(4926) + 
     &      MStop2(2)*(Cd(4927) + Cd(4928)*MStop2(2)))

	Cd(5295) = -Cd(5439) + Opt(9)

        Cd(5301) = -2*Cd(4924)*MStop2(2)**2 + 
     &   Cd(1755)*MStop2(3)**3 + Opt(9)

        Cd(5434) = 8*MTy2**4 - Cd(5439) - 16*MTy2**3*MStop2(1) + 
     &   14*MTy2**2*MStop2(1)**2 - MTy2*MStop2(1)**3 - 
     &   MStop2(2)*(Cd(5431) - Cd(5432)*MStop2(2) + 
     &      Cd(5433)*MStop2(2)**2)

        Cd(5606) = 162*MTy2**4 - Cd(5439) + 
     &   40*MTy2**3*MStop2(1) + 38*MTy2**2*MStop2(1)**2 + 
     &   50*MTy2*MStop2(1)**3 + 
     &   2*MStop2(2)*(Cd(5603) + 
     &      MStop2(2)*(Cd(5604) + Cd(5605)*MStop2(2)))

        Cd(5789) = 15*MTy2**5 + 43*MTy2**4*MStop2(1) - 
     &   12*MTy2**2*MStop2(1)**3 - 3*MTy2*MStop2(1)**4 - 
     &   Cd(5786)*MStop2(2) + 
     &   2*MStop2(2)**2*(Cd(5787) + 
     &      MStop2(2)*(Cd(5788) + 2*Cd(44)*MStop2(2))) - 
     &   Cd(5044)*MStop2(3)

        Cd(5938) = -Cd(5439) + 
     &   2*MTy2*(-4*MStop2(1)**3 + 
     &      MTy2*(MTy2*(15*MTy2 - 6*MStop2(1)) - 9*MStop2(1)**2))-
     &     2*MStop2(2)*(2*Cd(5936) + 
     &      MStop2(2)*(Cd(5937) + 2*Cd(2092)*MStop2(2)))

        Cd(6383) = MGl2*(MGl2*(MGl2 - Cd(44)) - Cd(3573)) + 
     &   Cd(47)*MGlpTmSt2(2) + 
     &   MTy2*(-MStop2(1)**2 + MTy2*(MTy2 + MStop2(1))) - 
     &   Cd(6382)*MStop2(2) - Cd(6381)*MStop2(3)

        Cd(6540) = MTy2*(-3*MStop2(1)**4 + 
     &      MTy2*(14*MStop2(1)**3 + 
     &         MTy2**2*(30*MTy2 + 27*MStop2(1)))) - 
     &   Cd(6538)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(6539) + 
     &      MStop2(2)*(-2*Cd(6130) + Cd(6131)*MStop2(2))) - 
     &   Cd(5044)*MStop2(3)

        Cd(6778) = MTy2**2*(5*MTy2 + 2*MStop2(1))*
     &    (MTy2*(9*MTy2 - 4*MStop2(1)) + 4*MStop2(1)**2) - 
     &   Cd(6775)*MStop2(2) + 2*Cd(6776)*MStop2(2)**2 + 
     &   2*Cd(6777)*MStop2(2)**3 + Cd(3856)*MStop2(2)**4 - 
     &   Cd(5044)*MStop2(3)

        Cd(7184) = 36*MTy2**5 - 3*Cd(7183) + 
     &   6*MTy2**4*MStop2(1) + 4*MTy2**3*MStop2(1)**2 + 
     &   22*MTy2**2*MStop2(1)**3 - 4*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(7179) - 2*Cd(7180)*MStop2(2) + 
     &      2*Cd(7181)*MStop2(2)**2)

        Cd(7535) = 14*MTy2**4 + Cd(5439) + 4*MTy2**3*MStop2(1) + 
     &   90*MTy2**2*MStop2(1)**2 - 28*MTy2*MStop2(1)**3 + 
     &   2*MStop2(2)*(2*Cd(7532) + 
     &      MStop2(2)*(Cd(7533) - 2*Cd(7534)*MStop2(2)))

        Cd(7703) = 83*MTy2**4 + Cd(5439) - 
     &   20*MTy2**3*MStop2(1) + 10*MTy2**2*MStop2(1)**2 + 
     &   3*MTy2*MStop2(1)**3 + 
     &   MStop2(2)*(-Cd(7701) + 
     &      MStop2(2)*(Cd(7702) + Cd(1684)*MStop2(2)))

	Cd(8205) = -Cd(5439) + Opt(10)

        Cd(8278) = -2*Cd(4924)*MStop2(2)**2 + 
     &   Cd(1755)*MStop2(3)**3 + Opt(10)

        Cd(8297) = 6*MTy2**2*(MTy2 - MStop2(1))**2*
     &    (4*MTy2 - MStop2(1)) - Cd(8296)*MStop2(2) + 
     &   2*Cd(6690)*MStop2(2)**2 - 2*Cd(6691)*MStop2(2)**3 + 
     &   Cd(1532)*MStop2(2)**4 - Cd(5044)*MStop2(3)

	Cd(8349) = -Cd(5439) + Opt(11)

        Cd(8384) = -2*Cd(4924)*MStop2(2)**2 + 
     &   Cd(1755)*MStop2(3)**3 + Opt(11)

        Cd(8409) = 10*MTy2**4 - Cd(5439) - 4*MTy2**3*MStop2(1) + 
     &   6*MTy2**2*MStop2(1)**2 + 2*Cd(8408)*MStop2(2)**2 + 
     &   4*MTy2*(-(Cd(8407)*MStop2(2)) + 
     &      (2*MStop2(1) + MStop2(3))*
     &       (MStop2(1)**2 + MStop2(2)*MStop2(3)))

        Cd(10254) = Cd(425) + 2*At*Cd(10253) - 4*MStop2(2) - 
     &   4*(-2*M_3C + AtC)*MTy*UCStop(1,3) + 
     &   8*M_3*MTy*UCStopC(1,3) - 
     &   8*AtC*(M_3C*UCStop(1,3)**2 + 
     &      M_3*UStop2(1,1)*
     &       (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Cd(10661) = -(At*Cd(10660)) + Cd(10355)*MStop2(3) + 
     &   At**2*UCStopC(1,3)**2 + 
     &   AtC*(-8*M_3C + AtC)*UCStop(1,3)**2*UStop2(1,1) + 
     &   4*(M_3*AtC*UStop2(1,1)**3 + 
     &      (2*M_3*AtC - 3*MStop2(1))*UStop2(1,1)**2*
     &       UStop2(1,2) + M_3*AtC*UStop2(1,1)*UStop2(1,2)**2 + 
     &      M_3*AtC*UStop2(1,2)**3)

        Cd(2311) = MStop2(1)*
     &    (22*MTy2**2 + MGl2*(-10*MTy2 - 7*MStop2(1)) + 
     &      9*MStop2(1)**2) + 2*Cd(2310)*MStop2(2) + 
     &   Cd(1606)*MStop2(2)**2 + Opt(12)

        Cd(2514) = -(MStop2(1)*
     &      (-16*MTy2**2 + MGl2*(4*MTy2 + MStop2(1)))) - 
     &   Cd(2512)*MStop2(2) + Cd(2513)*MStop2(2)**2 - 
     &   9*MStop2(2)**3 + Opt(12)

        Cd(3172) = 2*MTy2**2*(MTy2*Cd(3167) - Cd(3168))*
     &    UStop(1,1) + MStop2(3)*
     &    (-4*MTy*Cd(3171)*MStop2(2) - 
     &      2*(-3*MTy2*Cd(138)*Cd(141)*Cd(3169) + 
     &         Cd(607)*MStop2(2)**4 + 
     &         MGl2*(-2*Cd(1421)*MStop2(2)**3 + 
     &            MGl2*(12*Cd(843)*MStop2(2)**2 + 
     &               MGl2*(MGl2*Cd(874) - 2*Cd(1122)*MStop2(2)))))*
     &       UStop(1,1) + Cd(3166)*UStop2(1,2))

        Cd(3393) = -(MTy*Cd(3392)) - MTy2*Cd(3385)*MStop2(3) + 
     &   (-6*MTy2**4*Cd(3376) + MTy2**2*Cd(3379) + 
     &      MTy2**3*Cd(3381)*MStop2(3) + 
     &      (MGl2 - MStop2(2))**4*MStop2(2)*MStop2(3)**2)*
     &    UUStop(1,1)

        Cd(4935) = Cd(1272)*Cd(4934) - 
     &   2*MTy2*(6*MStop2(1)**4 + 
     &      MTy2*(-9*MStop2(1)**3 + 
     &         MTy2*(MTy2*(12*MTy2 - 18*MStop2(1)) + 
     &            14*MStop2(1)**2))) + 
     &   MStop2(2)*(Cd(4930) - 
     &      MStop2(2)*(2*Cd(4931) - 6*Cd(4932)*MStop2(2) + 
     &         Cd(1541)*MStop2(2)**2))

        Cd(5045) = MTy2*(-25*MStop2(1)**4 + 
     &      MTy2*(52*MStop2(1)**3 + 
     &         MTy2*(MTy2*(36*MTy2 - 29*MStop2(1)) + 
     &            2*MStop2(1)**2))) - Cd(5040)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(5041) + 
     &      MStop2(2)*(-2*Cd(5042) + Cd(1545)*MStop2(2))) + 
     &   Cd(5044)*MStop2(3)

        Cd(5097) = MTy2*(MTy2 - MStop2(1))*
     &    (23*MStop2(1)**3 + 
     &      MTy2*(MTy2*(27*MTy2 - 23*MStop2(1)) - 7*MStop2(1)**2))+
     &     MStop2(2)*(Cd(5093) + 
     &      MStop2(2)*(-2*Cd(5094) + 
     &         MStop2(2)*(-2*Cd(5095) + Cd(5096)*MStop2(2)))) + 
     &   Cd(5044)*MStop2(3)

        Cd(5244) = 36*MTy2**5 + 10*MTy2**4*MStop2(1) - 
     &   4*MTy2**3*MStop2(1)**2 - 30*MTy2**2*MStop2(1)**3 + 
     &   20*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(5240) + 
     &      MStop2(2)*(-2*Cd(5241) - 2*Cd(5242)*MStop2(2) + 
     &         Cd(5243)*MStop2(2)**2)) - Cd(5044)*MStop2(3)

        Cd(5275) = 72*MTy2**5 - Cd(1272)*Cd(4934) - 
     &   108*MTy2**4*MStop2(1) + 76*MTy2**3*MStop2(1)**2 - 
     &   34*MTy2**2*MStop2(1)**3 + 20*MTy2*MStop2(1)**4 + 
     &   MStop2(2)*(-Cd(5271) + 2*Cd(5272)*MStop2(2) + 
     &      MStop2(2)**2*(-2*Cd(5273) + Cd(5274)*MStop2(2)))

        Cd(5449) = 24*MTy2**5 + 3*Cd(1272)*Cd(4934) - 
     &   36*MTy2**4*MStop2(1) + 17*MTy2**3*MStop2(1)**2 - 
     &   2*MTy2**2*MStop2(1)**3 - 6*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(5445) - Cd(5446)*MStop2(2) + 
     &      MStop2(2)**2*(2*Cd(5447) + Cd(5448)*MStop2(2)))

        Cd(5508) = -24*MTy2**5 + 3*Cd(1272)*Cd(4934) + 
     &   36*MTy2**4*MStop2(1) - 31*MTy2**3*MStop2(1)**2 + 
     &   14*MTy2**2*MStop2(1)**3 - 14*MTy2*MStop2(1)**4 + 
     &   MStop2(2)*(Cd(5504) - 
     &      MStop2(2)*(Cd(5505) - 2*Cd(5506)*MStop2(2) + 
     &         Cd(5507)*MStop2(2)**2))

        Cd(5542) = 27*MTy2**5 + 23*MTy2**4*MStop2(1) - 
     &   10*MTy2**3*MStop2(1)**2 - 24*MTy2**2*MStop2(1)**3 + 
     &   20*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(5538) + 
     &      MStop2(2)*(-2*Cd(5539) - 2*Cd(5540)*MStop2(2) + 
     &         Cd(5541)*MStop2(2)**2)) - Cd(5044)*MStop2(3)

        Cd(5658) = 60*MTy2**5 + 378*MTy2**4*MStop2(1) + 
     &   76*MTy2**3*MStop2(1)**2 + 182*MTy2**2*MStop2(1)**3 + 
     &   80*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(5654) + 
     &      MStop2(2)*(2*Cd(5655) + 2*Cd(5656)*MStop2(2) + 
     &         Cd(5657)*MStop2(2)**2)) + Cd(5044)*MStop2(3)

        Cd(5678) = 108*MTy2**5 + 230*MTy2**4*MStop2(1) - 
     &   65*MTy2**3*MStop2(1)**2 - 218*MTy2**2*MStop2(1)**3 + 
     &   150*MTy2*MStop2(1)**4 - Cd(5674)*MStop2(2) + 
     &   MStop2(2)**2*(Cd(5675) + 2*Cd(5676)*MStop2(2) - 
     &      3*Cd(5677)*MStop2(2)**2) - 17*Cd(5044)*MStop2(3)

        Cd(5854) = 3*Cd(1272)*Cd(4934) + 
     &   2*MTy2*(MStop2(1)**4 + 
     &      MTy2*(3*MStop2(1)**3 + 
     &         MTy2*(MTy2*(18*MTy2 - 9*MStop2(1)) - 
     &            12*MStop2(1)**2))) - Cd(5850)*MStop2(2) - 
     &   2*Cd(5851)*MStop2(2)**2 + 2*Cd(5852)*MStop2(2)**3 + 
     &   Cd(5853)*MStop2(2)**4

        Cd(6035) = MTy2*(22*MStop2(1)**4 + 
     &      MTy2*(-40*MStop2(1)**3 + 
     &         MTy2*(-2*MStop2(1)**2 + 
     &            MTy2*(9*MTy2 + 15*MStop2(1))))) - 
     &   Cd(2092)*Cd(6032)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(6033) + 
     &      MStop2(2)*(2*Cd(6034) - Cd(5541)*MStop2(2))) - 
     &   Cd(5044)*MStop2(3)

        Cd(6132) = MTy2*(-MStop2(1)**4 + 
     &      MTy2*(22*MStop2(1)**3 + 
     &         MTy2*(8*MStop2(1)**2 + 
     &            MTy2*(60*MTy2 + 11*MStop2(1))))) - 
     &   Cd(6128)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(6129) + 
     &      MStop2(2)*(-2*Cd(6130) + Cd(6131)*MStop2(2))) - 
     &   Cd(5044)*MStop2(3)

        Cd(6149) = MTy2*(-8*MStop2(1)**4 + 
     &      MTy2*(10*MStop2(1)**3 + 
     &         MTy2*(MTy2*(24*MTy2 - 30*MStop2(1)) + 
     &            20*MStop2(1)**2))) - Cd(6145)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(6146) + 
     &      MStop2(2)*(-2*Cd(6147) + Cd(6148)*MStop2(2))) - 
     &   Cd(5044)*MStop2(3)

        Cd(6491) = 2*MTy2*
     &    (-6*MStop2(1)**4 + 
     &      MTy2*(9*MStop2(1)**3 + 
     &         MTy2*(MTy2*(6*MTy2 - 3*MStop2(1)) + 2*MStop2(1)**2))
     &      ) - Cd(6489)*MStop2(2) + 2*Cd(6490)*MStop2(2)**2 - 
     &   2*Cd(6147)*MStop2(2)**3 + Cd(6148)*MStop2(2)**4 - 
     &   Cd(5044)*MStop2(3)

        Cd(6549) = 12*MTy2**5 + 4*MTy2**4*MStop2(1) - 
     &   13*MTy2**3*MStop2(1)**2 + 22*MTy2**2*MStop2(1)**3 - 
     &   7*MTy2*MStop2(1)**4 - Cd(6546)*MStop2(2) + 
     &   MStop2(2)**2*(Cd(6547) - 2*Cd(6136)*MStop2(2) + 
     &      Cd(6548)*MStop2(2)**2) - Cd(5044)*MStop2(3)

        Cd(6635) = 12*MTy2**4*Cd(3378) + Cd(1272)*Cd(4934) - 
     &   6*MTy2*(2*MStop2(1)**4 + 
     &      MTy2*(-3*MStop2(1)**3 + 
     &         MTy2*(6*MTy2**2 + 6*MStop2(1)**2))) + 
     &   Cd(6633)*MStop2(1)*MStop2(2) - 2*Cd(6634)*MStop2(2)**2 + 
     &   2*Cd(6147)*MStop2(2)**3 - Cd(6148)*MStop2(2)**4

        Cd(6692) = 12*MTy2**4*Cd(286) + Cd(1272)*Cd(4934) - 
     &   2*MTy2*(6*MStop2(1)**4 + 
     &      MTy2*(-9*MStop2(1)**3 + 
     &         MTy2*(6*MTy2**2 + 10*MStop2(1)**2))) + 
     &   Cd(6689)*MStop2(1)*MStop2(2) - 2*Cd(6690)*MStop2(2)**2 + 
     &   2*Cd(6691)*MStop2(2)**3 - Cd(1532)*MStop2(2)**4

        Cd(6786) = 2*MTy2*
     &    (-2*MStop2(1)**4 + 
     &      MTy2*(5*MStop2(1)**3 + 
     &         MTy2*(MTy2*(9*MTy2 - 8*MStop2(1)) + 3*MStop2(1)**2))
     &      ) - Cd(6783)*MStop2(2) + 2*Cd(6784)*MStop2(2)**2 - 
     &   2*Cd(6785)*MStop2(2)**3 + Cd(6553)*MStop2(2)**4 - 
     &   Cd(5044)*MStop2(3)

        Cd(6880) = -(Cd(6878)*MStop2(2)) + 
     &   MStop2(2)**3*(-6*Cd(4932) + Cd(1541)*MStop2(2)) + 
     &   MTy2*(-6*MStop2(1)**4 + 
     &      MTy2*(6*MStop2(1)**3 + 
     &         MTy2*(MTy2*(18*MTy2 - 30*MStop2(1)) + 
     &            20*MStop2(1)**2)) + 2*Cd(6879)*MStop2(2)**2) - 
     &   Cd(5044)*MStop2(3)

        Cd(6937) = 60*MTy2**5 - 59*MTy2**4*MStop2(1) + 
     &   32*MTy2**3*MStop2(1)**2 - 6*MTy2**2*MStop2(1)**3 + 
     &   5*MTy2*MStop2(1)**4 + 
     &   MStop2(2)*(Cd(6934) + 
     &      MStop2(2)*(2*Cd(6935) + 2*Cd(6936)*MStop2(2) - 
     &         Cd(3770)*MStop2(2)**2)) + Cd(5044)*MStop2(3)

        Cd(6941) = 24*MTy2**5 - 52*MTy2**4*MStop2(1) + 
     &   45*MTy2**3*MStop2(1)**2 - 30*MTy2**2*MStop2(1)**3 + 
     &   11*MTy2*MStop2(1)**4 - Cd(6938)*MStop2(2) + 
     &   MStop2(2)**2*(Cd(6939) + 2*Cd(6940)*MStop2(2) - 
     &      Cd(5034)*MStop2(2)**2) + Cd(5044)*MStop2(3)

        Cd(6951) = 2*MTy2*(MTy2 - MStop2(1))**2*
     &    (8*MStop2(1)**2 + MTy2*(12*MTy2 + 3*MStop2(1))) - 
     &   MStop2(2)*(Cd(6949) + 
     &      MStop2(2)*(-2*Cd(6950) + 
     &         MStop2(2)*(-2*Cd(6691) + Cd(1532)*MStop2(2)))) + 
     &   Cd(5044)*MStop2(3)

        Cd(7072) = MTy2*(MStop2(1)**4 + 
     &      MTy2*(-22*MStop2(1)**3 + 
     &         MTy2*(MTy2*(27*MTy2 - 46*MStop2(1)) + 
     &            8*MStop2(1)**2))) + 
     &   MStop2(2)*(Cd(7068) - 
     &      MStop2(2)*(2*Cd(7069) + 
     &         MStop2(2)*(-2*Cd(7070) + Cd(7071)*MStop2(2)))) + 
     &   3*Cd(5044)*MStop2(3)

        Cd(7259) = Cd(1272)*Cd(4934) + 
     &   2*MTy2*(2*MStop2(1)**4 + 
     &      MTy2*(-7*MStop2(1)**3 + 
     &         MTy2*(MTy2*(36*MTy2 - 54*MStop2(1)) + 
     &            34*MStop2(1)**2))) - Cd(7255)*MStop2(2) + 
     &   2*Cd(7256)*MStop2(2)**2 - 2*Cd(7257)*MStop2(2)**3 + 
     &   Cd(7258)*MStop2(2)**4

        Cd(7410) = 24*MTy2**5 - 7*Cd(1272)*Cd(4934) - 
     &   36*MTy2**4*MStop2(1) + 47*MTy2**3*MStop2(1)**2 - 
     &   22*MTy2**2*MStop2(1)**3 + 14*MTy2*MStop2(1)**4 + 
     &   MStop2(2)*(-Cd(7406) + 
     &      MStop2(2)*(Cd(7407) - 2*Cd(7408)*MStop2(2) + 
     &         Cd(7409)*MStop2(2)**2))

        Cd(7457) = MTy2*(-4*MStop2(1)**4 + 
     &      MTy2*(28*MStop2(1)**3 + 
     &         MTy2*(-2*MStop2(1)**2 + 
     &            MTy2*(27*MTy2 + 19*MStop2(1))))) - 
     &   Cd(7453)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(7454) + 
     &      MStop2(2)*(-2*Cd(7455) + Cd(7456)*MStop2(2))) - 
     &   3*Cd(5044)*MStop2(3)

        Cd(7569) = 6*MTy2**5 + 12*MTy2**4*MStop2(1) - 
     &   8*MTy2**3*MStop2(1)**2 - 2*MTy2**2*MStop2(1)**3 + 
     &   4*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(7565) + 
     &      MStop2(2)*(-2*Cd(7566) + 2*Cd(7567)*MStop2(2) + 
     &         Cd(7568)*MStop2(2)**2)) + Cd(5044)*MStop2(3)

        Cd(7589) = 15*MTy2**5 + 91*MTy2**4*MStop2(1) - 
     &   8*MTy2**3*MStop2(1)**2 + 12*MTy2**2*MStop2(1)**3 - 
     &   3*MTy2*MStop2(1)**4 - Cd(7586)*MStop2(2) + 
     &   2*MStop2(2)**2*(Cd(7587) - Cd(7588)*MStop2(2) + 
     &      2*Cd(44)*MStop2(2)**2) - Cd(5044)*MStop2(3)

        Cd(7640) = Cd(1272)*Cd(4934) + 
     &   2*MTy2*(5*MStop2(1)**4 + 
     &      MTy2*(-3*MStop2(1)**3 + 
     &         MTy2*(MTy2*(18*MTy2 - 9*MStop2(1)) + 8*MStop2(1)**2)
     &         )) - Cd(7637)*MStop2(2) + 2*Cd(7638)*MStop2(2)**2 - 
     &   2*Cd(7639)*MStop2(2)**3 + Cd(6625)*MStop2(2)**4

        Cd(7787) = 9*MTy2**5 + 11*MTy2**4*MStop2(1) + 
     &   6*MTy2**3*MStop2(1)**2 + 12*MTy2**2*MStop2(1)**3 - 
     &   2*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(7784) + 
     &      MStop2(2)*(2*Cd(7785) + 2*Cd(7786)*MStop2(2) - 
     &         Cd(7456)*MStop2(2)**2)) - 3*Cd(5044)*MStop2(3)

        Cd(7874) = 60*MTy2**5 + 51*MTy2**4*MStop2(1) - 
     &   4*MTy2**3*MStop2(1)**2 + 26*MTy2**2*MStop2(1)**3 + 
     &   7*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(7871) + 
     &      MStop2(2)*(-2*Cd(7872) + 2*Cd(7873)*MStop2(2) + 
     &         Cd(2763)*MStop2(2)**2)) - Cd(5044)*MStop2(3)

        Cd(7891) = Cd(7183) + 
     &   2*MTy2*(4*MStop2(1)**4 + 
     &      MTy2*(-5*MStop2(1)**3 + 
     &         MTy2*(MTy2*(12*MTy2 - 9*MStop2(1)) + 6*MStop2(1)**2)
     &         )) - MStop2(2)*
     &    (Cd(7888) - 2*MStop2(2)*
     &       (Cd(7889) - 3*Cd(7890)*MStop2(2)))

        Cd(8163) = 12*MTy2**5 + Cd(7183) + 6*MTy2**4*MStop2(1) - 
     &   4*MTy2**3*MStop2(1)**2 - 2*MTy2**2*MStop2(1)**3 + 
     &   4*MTy2*MStop2(1)**4 - 
     &   MStop2(2)*(Cd(8161) - 2*Cd(8162)*MStop2(2) + 
     &      6*Cd(7890)*MStop2(2)**2)

        Cd(9361) = -2*M_3*Cd(9325)*MStop2(1)**3*MStop2(3) + 
     &   MTy2*(2*Cd(9360) + 
     &      MTy2*(2*M_3*Cd(9351) + 
     &         MTy2*(2*M_3*Cd(9337) + 
     &            MTy2*(-Cd(9341) + 16*M_3*MTy2*UUStop(1,3)))))+
     &     MGl2*(2*M_3*Cd(9331)*MStop2(1)**2*MStop2(3) + 
     &      MGl2*(-6*M_3*Cd(9327)*MStop2(1)*MStop2(3) + 
     &         MGl2*(2*M_3*Cd(9329)*MStop2(3) - 
     &            2*M_3*MGl2*MStop2(3)*UUStop(2,3))))

        Cd(9625) = -2*M_3*Cd(9591)*MStop2(2)**3*MStop2(3) + 
     &   MTy2*(2*Cd(9624) + 
     &      MTy2*(2*M_3*Cd(9616) + 
     &         MTy2*(2*M_3*Cd(9602) + 
     &            MTy2*(-Cd(9606) + 16*M_3*MTy2*UUStop(1,3)))))+
     &     MGl2*(2*M_3*Cd(9595)*MStop2(2)**2*MStop2(3) + 
     &      MGl2*(-6*M_3*Cd(9593)*MStop2(2)*MStop2(3) + 
     &         MGl2*(2*M_3*Cd(9597)*MStop2(3) - 
     &            2*M_3*MGl2*MStop2(3)*UUStop(2,3))))

        Cd(10416) = Cd(9948)*Cd(10415) - Opt(16) - 
     &   3*M_3C*AtC*UCStop(1,3)**3 + 
     &   UCStop(1,3)*(-(MTy2*UStop2(1,1)*UStop2(1,2)) - 
     &      MStop2(1)*UStop2(1,2)**2 + 
     &      AtC*UStop2(1,1)*(M_3*UStop2(1,1) - M_3*UStop2(1,2))
     &      )

        Cd(12335) = Opt(14) + Opt(15) - 
     &   2*MTy*Cd(10234)*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Cd(12366) = Opt(14) - Opt(15) + 
     &   2*MTy*Cd(10234)*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2)

        Cd(12452) = Cd(9948)*Cd(12451) - Opt(16) - 
     &   3*M_3C*AtC*UCStop(1,3)**3 + 
     &   UCStop(1,3)*(MStop2(1)*UStop2(1,1)*UStop2(1,2) + 
     &      MTy2*UStop2(1,2)**2 + 
     &      AtC*UStop2(1,1)*(M_3*UStop2(1,1) - M_3*UStop2(1,2))
     &      )

        Cd(13656) = -(MTy*MUE*TB*Cd(13628)*Cd(13641)*Cr(9)) + 
     &   MUEC*TB*Cd(13641)*Cd(13649)*Cr(9) + 
     &   At*Cd(13641)*Cd(13651)*Cr(9) - 
     &   AtC*Cd(13641)*Cd(13653)*Cr(9) - 
     &   Cd(13643)*UStop2(1,1)*UStop2(1,2) + 
     &   2*Cd(12700)*Cd(12701)*Cd(13641)*Cr(9)*UStop2(1,1)*
     &    UStop2(1,2) - Cd(13638)*MSq2Diff(1,2)*UStop2(1,1)*
     &    UStop2(1,2) - 2*Cd(491)*MSq2Diff(1,1)*MSq2Diff(1,2)*
     &    UStop2(1,1)*UStop2(1,2) - 
     &   Cd(13640)*UStop2(1,1)**2*UStop2(1,2) - 
     &   Cd(13646)*UCStop(3,1)*UStopC(1,1) + 
     &   Cd(13655)*MSq2Diff(1,1)*UStop2(1,2)*UUStopC(3,3)

        Cd(1891) = 2*MGl2**3 - 6*MGl2**2*(MTy2 + MStop2(1)) + 
     &   MGl2*(4*MTy2**2 + 5*MStop2(1)**2) + 
     &   (MTy2 + MStop2(1))*
     &    (MStop2(1)**2 + MTy2*(4*MTy2 + 8*MStop2(1))) + 
     &   MStop2(2)*(-2*Cd(1889) + Cd(1890)*MStop2(2))

        Cd(1926) = 6*(MGl2 - 2*MTy2)*
     &    (MGl2*(MGl2 - 2*MTy2) + 2*MTy2**2) + 
     &   3*(MTy2**2 + MGl2*(-MGl2 + MTy2))*MStop2(1) - 
     &   2*MTy2*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(1924) + 
     &      MStop2(2)*(-Cd(1925) + 15*MStop2(2)))

        Cd(5217) = -(MGl2*Cd(5212)) + 2*MGl2**2*Cd(5213) + 
     &   MGl2**3*(-(MGl2*Cd(5210)) + Cd(5216)) + 
     &   (MGl2**5 - MStop2(1)*MStop2(2)**4 + 
     &      MTy2*(Cd(594)*MStop2(2)**3 + 
     &         MTy2*(-2*Cd(595)*MStop2(2)**2 + 
     &            MTy2*(MTy2*(MTy2 - Cd(552)) + 
     &               2*Cd(2802)*MStop2(2)))))*MStop2(3)

        Cd(9214) = MTy*Cd(9213) - 
     &   (2*MTy2*(Cd(9202) + MTy2*(8*Cd(9198) + MTy2*Cd(9205))) + 
     &      MTy2**4*(-24*MTy2 + Cd(9197))*MStop2(3) + 
     &      (MGl2 - MStop2(1))**4*MStop2(3)**2)*UUStop(1,1)

        Cd(10884) = MTy2*Cd(10751) + 2*M_3*MTy*Cd(10883) + 
     &   MUE**2*MStop2(1)*UCStop(1,3)**2*UStop2(1,1) + 
     &   MUEC**2*UCStopC(1,3)**2*
     &    (MStop2(1) + MStop2(3)*UStop2(1,2)) + 
     &   2*M_3C*MTy*(-3*MUE**2*UCStop(1,3)**3 + 
     &      2*MUE2*UCStop(1,3)*
     &       (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &         UStop2(1,2)**2) + 
     &      MUEC**2*UCStopC(1,3)*
     &       (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &         UStop2(1,2)**2))

        Cd(11786) = At*Cd(11785) - 
     &   2*MUE*(AtC*Cd(47)*(MGlmT2 + MStop2(1)) - 
     &      Cd(152)*(2*M_3C*Cd(152) + AtC*MStop2(1))*MStop2(3))*
     &    UCStop(1,3)**2 + 
     &   2*Cd(152)*Cd(9948)*Cd(11784)*MStop2(3)*UStop2(1,1) + 
     &   MUEC*(4*M_3*Cd(152)**2*MStop2(3)*UCStopC(1,3)**2 + 
     &      AtC*(2*Cd(152)*MStop2(1)*MStop2(3)*UStop2(1,1)*
     &          UStop2(1,2)**2 + 
     &         Cd(47)*(MGlmT2 + MStop2(1))*
     &          (UStop2(1,1)**2 + UStop2(1,2)**2)))

        Cd(11927) = -Opt(16) + 
     &   M_3*At*UCStopC(1,3)*UStop2(1,2)**2 + 
     &   UCStop(1,3)*(AtC*Cd(10021) + 
     &      Cd(787)*MStop2(3)*UStop2(1,1) - 
     &      MTy2*UStop2(1,1)*UStop2(1,2) + 
     &      M_3C*At*UStop2(1,2)**2 + MGlpT2*UStop2(1,2)**2)

        Cd(12402) = 4*Cd(12370) + 
     &   (MGlmT2*Cd(47)*Cd(11750))/MStop2(3) - 
     &   UStop2(1,1)**2*(-5*MGlmT2*Cd(47) + 3*Cd(9948)*Cd(12391) + 
     &      (4*MGlmT2*Cd(47) + Cd(12371))*UStop2(1,1)) + 
     &   10*MGlmT2*Cd(47)*UStop2(1,1)*UStop2(1,2) + 
     &   (5*MGlmT2*Cd(47) - 
     &      3*(4*MGlmT2*Cd(47) + Cd(12371))*UStop2(1,1))*
     &    UStop2(1,2)**2 - 
     &   (4*MGlmT2*Cd(47) + Cd(12371))*UStop2(1,2)**3

	Cd(12449) = Cd(10416) - MStop2(3)*UCStop(1,3)

        Cd(14064) = MTy2*Cd(14063)*MSq2Diff(1,1)*MStop2(1) + 
     &   MGl2*Cd(14058)*MSq2Diff(1,1)*MStop2(3)*UStop2(1,2) + 
     &   MGl2**2*TB*MSq2Diff(1,1)**2*MStop2(3)*UStop2(1,2)**2 + 
     &   At*MUE*TB2*Cd(47)*UStop2(1,2)*
     &    (MStop2(1)*MStop2(3)*UStop2(1,2) + 
     &      MSq2Diff(1,1)*(-(MStop2(1)*UStop2(1,1)) + 
     &         MStop2(3)*UStop2(1,2))) + 
     &   MTy*((At*TB + MUEC*TB2)*Cd(47)*MSq2Diff(1,1)*MStop2(1)*
     &       UCStopC(1,3)*UStop2(1,2) + 
     &      MUE*TB2*Cd(47)*UCStop(1,3)*
     &       (Cd(9948)*MSq2Diff(1,1)*MStop2(3) + 
     &         MStop2(1)*(-(MSq2Diff(1,1)*UStop2(1,1)) + 
     &            MStop2(3)*UStop2(1,2))) - 
     &      TB*Cd(13469)*Cd(14007)*MStop2(3)*UStop2(1,1)*
     &       UStopC(1,1))

        Cd(2054) = MStop2(1)*
     &    (-11*MStop2(1)**2 + MTy2*(8*MTy2 + 9*MStop2(1)) + 
     &      MGl2*(-8*MGl2 + 13*MStop2(1))) - 
     &   2*Cd(2053)*MStop2(2) + Cd(1890)*MStop2(2)**2 + Opt(17)

        Cd(6153) = -4*(-20*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - 6*MTy2) + 15*MTy2**2)) + 
     &   2*Cd(6152) + 10*MGlmT2*(MGl2 - 2*MTy2)*MStop2(1) + 
     &   (-9*MGl2 + 20*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(2*Cd(6150) - Cd(6151)*MStop2(2))

        Cd(6971) = -4*(-20*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - 6*MTy2) + 15*MTy2**2)) + 
     &   2*Cd(6970) + 14*MGlmT2*(MGl2 - 2*MTy2)*MStop2(1) + 
     &   (-19*MGl2 + 28*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(2*Cd(6969) - Cd(6953)*MStop2(2))

        Cd(14179) = At*Cd(14142)*Cd(14178) - 
     &   2*M_3C*AtC*(Cd(47) - MGlpTmSt2(1)*MSq2Diff(1,1))*
     &    UCStop(1,3)**2 + 
     &   Cd(9948)*Cd(14047)*Cd(14142)*UStop2(1,1) + 
     &   2*Cd(47)*MSq2Diff(1,1)*UStop2(1,1) + 
     &   Cd(164)*MSq2Diff(1,1)**2*UStop2(1,1) - 
     &   Cd(47)*MSq2Diff(1,1)*UStop2(1,1)**2 + 
     &   Cd(164)*MSq2Diff(1,1)**2*UStop2(1,2) - 
     &   2*M_3*AtC*Cd(47)*UStop2(1,1)*UStop2(1,2) + 
     &   2*M_3*AtC*MGlpTmSt2(1)*MSq2Diff(1,1)*UStop2(1,1)*
     &    UStop2(1,2) + MUE2*TB2*Cd(47)*UStop2(1,2)**2 - 
     &   MUE2*TB2*MGlpTmSt2(1)*MSq2Diff(1,1)*UStop2(1,2)**2 + 
     &   MGlpTmSt2(1)*MSq2Diff(1,1)**2*UStop2(1,2)**2

        Cd(5838) = 4*MTy2**3*Cd(1272)*Cd(5814) + 
     &   MGl2*(MGl2*(-2*Cd(5824) + 
     &         MGl2*(MGl2*(MGl2*(-2*MGl2 + Cd(5817)) - 
     &               2*Cd(5818)) + 2*Cd(5827))) + 2*Cd(5835)) + 
     &   MTy2*Cd(5837) - 2*
     &    (Cd(1280)*MStop2(1)**2*MStop2(2)**2 + 
     &      MTy2**2*(11*Cd(831)*MStop2(1)*MStop2(2) + 
     &         MTy2**2*(Cd(5816) + MTy2*(9*MTy2 - 23*MStop2(4)))))

        Cd(9626) = MTy*Cd(9625) + 
     &   MTy2*(Cd(9589) + MTy2*
     &       (Cd(9574) + MTy2*
     &          (-Cd(9581) + 
     &            MTy2*(Cd(9568) + 
     &               MTy2*(Cd(9561) - 12*MTy2*UUStop(1,1)))))) + 
     &   MStop2(3)*(-(Cd(3937)*MStop2(2)**4) + 
     &      MGl2*(Cd(9562)*MStop2(2)**3 + 
     &         MGl2*(-2*Cd(9563)*MStop2(2)**2 + 
     &            MGl2*(2*Cd(9564)*MStop2(2) + 
     &               MGl2*(-Cd(9565) + MGl2*UUStop(2,1))))))

        Cd(11813) = 2*Cd(9948)*Cd(11810) - 
     &   MGlmT2*(MStop2(3) - 2*M_3C*At*UStop2(1,1)**2) - 
     &   2*(Cd(11812)*MStop2(1) + 
     &      3*M_3*At*MGlmT2*UCStopC(1,3)**2 + 
     &      UStop2(1,2)*(-(Cd(9955)*MStop2(1)**2) + 
     &         M_3C*At*MGlmT2*UStop2(1,1) + 
     &         MGl2*MGlmT2*UStop2(1,1) - MGlmT2*MTy2*UStop2(1,2) + 
     &         Cd(47)*UStop2(1,2)))

        Cd(10061) = At*Cd(10057) - 4*Cd(9948)*Cd(10060) - 
     &   2*AtC*MUE*UCStop(1,3)**2 + 
     &   8*M_3C*MUE*UCStop(1,3)**2*(1 + UStop2(1,2)) + 
     &   8*M_3*MUEC*UCStopC(1,3)**2*(1 + UStop2(1,2)) + 
     &   AtC*MUEC*(UStop2(1,1)**2 + UStop2(1,2)**2) - 
     &   4*(M_3*MUE + M_3C*MUEC)*UStop2(1,1)*
     &    (UStop2(1,1)**2 + UStop2(1,2)**2)

        Cd(5918) = 4*MGl2*Cd(5904) + 
     &   MGl2**2*(2*Cd(5914) + 
     &      MGl2*(4*Cd(5895) + 
     &         MGl2*(-4*MGl2*
     &             (MGl2*(2*MGl2 - Cd(5861)) + Cd(5890)) + Cd(5917)
     &            ))) + 2*Cd(1272)*
     &    (MTy2**2*(Cd(5886) + 5*MTy2**2*Cd(5888)) + 
     &      Cd(5062)*MStop2(1)**2*MStop2(2)**2) - 
     &   4*MTy2*(MTy2**2*Cd(5881) + Cd(5883)*MStop2(1)*MStop2(2) + 
     &      MTy2**4*(2*Cd(5879) + MTy2*(2*MTy2 - 7*MStop2(4)))) + 
     &   Cd(302)*Cd(5877)*UStop2(1,2)

        Cd(14139) = At*Cd(14135) - MTy2*Cd(14138) + 
     &   MGl2*(AtC*MUEC*TB2 + MUE2*TB*TB2 - TB*MSq2Diff(1,1))*
     &    UStop2(1,2)**2 + 
     &   MTy*UCStopC(1,3)*(MUEC*TB2*(MGl2 - MStop2(1)) - 
     &      2*M_3*(AtC*MUEC*TB2 + MUE2*TB*TB2 - 
     &         TB*MSq2Diff(1,1))*UStop2(1,2)) + 
     &   MTy*UCStop(1,3)*((AtC*TB + MUE*TB2)*(MGl2 - MStop2(1)) - 
     &      2*M_3C*(AtC*MUEC*TB2 + MUE2*TB*TB2 - 
     &         TB*MSq2Diff(1,1))*UStop2(1,2))

        Cd(11820) = 2*Cd(9948)*Cd(11819) - 
     &   2*Cd(11816)*MStop2(1) - 
     &   6*M_3C*AtC*MGlmT2*UCStop(1,3)**3 + 
     &   MGlmT2*UCStop(1,3)*
     &    (-MStop2(3) + 2*M_3C*At*UStop2(1,1)**2 + 
     &      2*M_3*AtC*UStop2(1,1)**2) + 
     &   2*M_3*At*MGlmT2*UCStopC(1,3)*UStop2(1,1)*
     &    (UStop2(1,1) - UStop2(1,2)) + 
     &   2*UCStop(1,3)*UStop2(1,2)*
     &    (Cd(9955)*MStop2(1)**2 - 
     &      (M_3C*At + M_3*AtC + MGl2)*MGlmT2*UStop2(1,1) + 
     &      (MGlmT2*MTy2 - Cd(47))*UStop2(1,2))

        Cd(4951) = 2*(-3*MTy2**4 + 
     &      MGl2*(2*MTy2**3 + 
     &         MGl2*(MGl2*(MGl2 - 3*MTy2) + 2*MTy2**2))) - 
     &   (8*MTy2**3 + MGl2*(-30*MTy2**2 + MGl2*(7*MGl2 + 8*MTy2)))*
     &    MStop2(1) + 2*(11*MTy2**2 + MGl2*(5*MGl2 + 11*MTy2))*
     &    MStop2(1)**2 - (7*MGl2 + 18*MTy2)*MStop2(1)**3 + 
     &   2*MStop2(1)**4 - MStop2(2)*
     &    (Cd(4946) + MStop2(2)*(-Cd(4949) + Cd(4950)*MStop2(2)))

        Cd(12383) = Cd(9948)*Cd(12382) - 
     &   2*M_3C*AtC*Cd(47)*UCStop(1,3)**2 - 
     &   At*Cd(47)*Cd(10001)*UCStopC(1,3)**2 + 
     &   Cd(12369)*UStop2(1,1) - 
     &   AtC**2*Cd(47)*UCStop(1,3)**2*UStop2(1,1) - 
     &   MGlmT2*Cd(47)*UStop2(1,1)**2 - 
     &   4*Cd(12370)*UStop2(1,1)**2 - 
     &   AtC**2*Cd(47)*UCStop(1,3)**2*UStop2(1,2) - 
     &   2*M_3C*At*Cd(47)*UStop2(1,1)*UStop2(1,2) - 
     &   2*M_3*AtC*Cd(47)*UStop2(1,1)*UStop2(1,2) - 
     &   2*MGlmT2*Cd(47)*UStop2(1,1)*UStop2(1,2) - 
     &   8*Cd(12370)*UStop2(1,1)*UStop2(1,2) - 
     &   2*At*AtC*Cd(47)*UStop2(1,1)**2*UStop2(1,2) - 
     &   MGlmT2*Cd(47)*UStop2(1,2)**2 - 
     &   4*Cd(12370)*UStop2(1,2)**2 - 
     &   2*At*AtC*Cd(47)*UStop2(1,1)*UStop2(1,2)**2 - 
     &   Cd(47)*Cd(118)*(4*UStop2(1,1)**2 + 4*UStop2(1,2)**2 + 
     &      UStop2(1,1)*(-5 + 8*UStop2(1,2)))

        Cd(12333) = -(At*Cd(164)*Cd(2597)*Cd(12332)) - 
     &   AtC*MGl2*Cd(164)*Cd(2597)*
     &    (MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)) + 
     &   MTy2*Cd(164)*Cd(2597)*
     &    (Cd(9948)*Cd(11800)*UStop2(1,1)**2 + 
     &      4*M_3*(MUEC*UCStopC(1,3)**2 + 
     &         MUE*UStop2(1,1)*UStop2(1,2)**2) - 
     &      (-4*M_3C + AtC)*
     &       (MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)**2)
     &      ) + 2*MTy*(MUEC*Cd(12329)*UCStopC(1,3) + 
     &      UCStop(1,3)*(MUE*Cd(47)**2 - MUE*Cd(12298)*MStop2(1) + 
     &         Cd(164)*Cd(2597)*
     &          (-(MGl2*MUE) - MTy2*MUE + MUE*MStop2(1) + 
     &            AtC*(M_3C*MUEC*UStop2(1,1)*UStop2(1,2) + 
     &               MUE*(M_3C*UCStop(1,3)**2 + 
     &                  M_3*UStop2(1,1)*UStop2(1,2))))))

        Cd(12353) = AtC**2*MGl2*UCStop(1,3)**2 + 
     &   At*(Cd(12313)*UCStopC(1,3)**2 + 
     &      2*Cd(12351)*UStop2(1,1)) + 
     &   MTy2*(AtC*(-8*M_3C + AtC)*UCStop(1,3)**2 + 
     &      4*(3*Cd(138)*Cd(9948)*UStop2(1,1)**2 + 
     &         (MGl2 - MStop2(2))*UStop2(1,1)**3 - 
     &         2*M_3*AtC*UStop2(1,1)*UStop2(1,2) + 
     &         3*(MGl2 - MStop2(2))*UStop2(1,1)*UStop2(1,2)**2 + 
     &         (MGl2 - MStop2(2))*UStop2(1,2)**3)) + 
     &   2*MTy*(Cd(12352)*UCStopC(1,3) + 
     &      M_3*At**2*UCStopC(1,3)**3 + 
     &      UCStop(1,3)*(4*M_3C*MTy2 + 
     &         AtC*(-2*MGlpTmSt2(2) + 
     &            AtC*(M_3C*UCStop(1,3)**2 + 
     &               M_3*UStop2(1,1)*UStop2(1,2)))))

        Cd(12364) = At*Cd(251)*Cd(2681)*Cd(12363) + 
     &   AtC*MGl2*Cd(251)*Cd(2681)*
     &    (MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)) + 
     &   MTy2*Cd(251)*Cd(2681)*
     &    (-(Cd(9948)*Cd(11800)*UStop2(1,1)**2) - 
     &      4*M_3*(MUEC*UCStopC(1,3)**2 + 
     &         MUE*UStop2(1,1)*UStop2(1,2)**2) + 
     &      (-4*M_3C + AtC)*
     &       (MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)**2)
     &      ) + 2*MTy*(MUEC*Cd(12360)*UCStopC(1,3) + 
     &      UCStop(1,3)*(MUE*Cd(50)**2 - MUE*Cd(4187)*MStop2(2) + 
     &         Cd(251)*Cd(2681)*
     &          (-(MGl2*MUE) - MTy2*MUE + MUE*MStop2(2) + 
     &            AtC*(M_3C*MUEC*UStop2(1,1)*UStop2(1,2) + 
     &               MUE*(M_3C*UCStop(1,3)**2 + 
     &                  M_3*UStop2(1,1)*UStop2(1,2))))))

        Cd(12315) = 4*Cd(12298)*Cd(12314)*MStop2(1) + 
     &   AtC**2*MGl2*Cd(164)*Cd(2597)*UCStop(1,3)**2 - 
     &   2*MTy*Cd(12312)*UCStopC(1,3) - 
     &   2*M_3*At**2*MTy*Cd(164)*Cd(2597)*UCStopC(1,3)**3 + 
     &   At*Cd(164)*Cd(2597)*
     &    (Cd(12313)*UCStopC(1,3)**2 + 2*Cd(12310)*UStop2(1,1)) + 
     &   MTy2*(AtC*(-8*M_3C + AtC)*Cd(164)*Cd(2597)*
     &       UCStop(1,3)**2 + 
     &      4*UStop2(1,1)**2*
     &       (Cd(164)*Cd(2597)*(MGl2 - MStop2(1))*
     &          (3*Cd(9948) + UStop2(1,1)) + 
     &         Cd(47)**2*(-3 + 2*UStop2(1,1))) - 
     &      8*M_3*AtC*Cd(164)*Cd(2597)*UStop2(1,1)*UStop2(1,2) - 
     &      12*(Cd(47)**2 - MGl2*Cd(164)*Cd(2597) + 
     &         Cd(164)*Cd(2597)*MStop2(1))*UStop2(1,1)*
     &       UStop2(1,2)**2 - 
     &      4*(Cd(47)**2 - MGl2*Cd(164)*Cd(2597) + 
     &         Cd(164)*Cd(2597)*MStop2(1))*UStop2(1,2)**3) - 
     &   2*MTy*UCStop(1,3)*
     &    (4*M_3C*MTy2*Cd(164)*Cd(2597) + 
     &      AtC*(2*Cd(47)**2 + 
     &         Cd(164)*Cd(2597)*
     &          (-2*MGlpTmSt2(1) + 
     &            AtC*(M_3C*UCStop(1,3)**2 + 
     &               M_3*UStop2(1,1)*UStop2(1,2)))))

        Cd(10168) = -2*Cd(47)*Cd(50)**2*Cd(9947)*Cd(10166) - 
     &   2*Cd(47)**2*Cd(50)*Cd(9949)*Cd(10167) + 
     &   8*M_3C*AtC*Cd(50)**2*Cd(9953)*UCStop(1,3)**2 + 
     &   8*M_3C*AtC*Cd(47)**2*Cd(9954)*UCStop(1,3)**2 + 
     &   At*((-8*M_3 + At)*
     &       (Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952)) + 
     &      8*M_3*Cd(50)**2*Cd(9953) + 8*M_3*Cd(47)**2*Cd(9954)
     &      )*UCStopC(1,3)**2 + 
     &   8*M_3C*At*Cd(50)**2*Cd(9948)*Cd(9953)*UStop2(1,1) + 
     &   8*M_3*AtC*Cd(50)**2*Cd(9948)*Cd(9953)*UStop2(1,1) + 
     &   8*M_3C*At*Cd(47)**2*Cd(9948)*Cd(9954)*UStop2(1,1) + 
     &   8*M_3*AtC*Cd(47)**2*Cd(9948)*Cd(9954)*UStop2(1,1) - 
     &   8*M_3C*AtC*Cd(50)**2*Cd(9951)*UCStop(1,3)**2*
     &    UStop2(1,1) + AtC**2*Cd(50)**2*Cd(9951)*UCStop(1,3)**2*
     &    UStop2(1,1) - 8*M_3C*AtC*Cd(47)**2*Cd(9952)*
     &    UCStop(1,3)**2*UStop2(1,1) + 
     &   AtC**2*Cd(47)**2*Cd(9952)*UCStop(1,3)**2*UStop2(1,1) - 
     &   8*M_3C*At*Cd(50)**2*Cd(9948)*Cd(9951)*UStop2(1,1)**2 - 
     &   8*M_3*AtC*Cd(50)**2*Cd(9948)*Cd(9951)*UStop2(1,1)**2 - 
     &   8*M_3C*At*Cd(47)**2*Cd(9948)*Cd(9952)*UStop2(1,1)**2 - 
     &   8*M_3*AtC*Cd(47)**2*Cd(9948)*Cd(9952)*UStop2(1,1)**2 - 
     &   8*M_3C*AtC*Cd(50)**2*Cd(9951)*UCStop(1,3)**2*
     &    UStop2(1,2) + AtC**2*Cd(50)**2*Cd(9951)*UCStop(1,3)**2*
     &    UStop2(1,2) - 8*M_3C*AtC*Cd(47)**2*Cd(9952)*
     &    UCStop(1,3)**2*UStop2(1,2) + 
     &   AtC**2*Cd(47)**2*Cd(9952)*UCStop(1,3)**2*UStop2(1,2) + 
     &   2*At*AtC*Cd(50)**2*Cd(9951)*UStop2(1,1)**2*UStop2(1,2) + 
     &   2*At*AtC*Cd(47)**2*Cd(9952)*UStop2(1,1)**2*UStop2(1,2) - 
     &   8*M_3C*At*Cd(50)**2*Cd(9951)*UStop2(1,1)*
     &    UStop2(1,2)**2 - 
     &   8*M_3*AtC*Cd(50)**2*Cd(9951)*UStop2(1,1)*
     &    UStop2(1,2)**2 + 
     &   2*At*AtC*Cd(50)**2*Cd(9951)*UStop2(1,1)*UStop2(1,2)**2 - 
     &   8*M_3C*At*Cd(47)**2*Cd(9952)*UStop2(1,1)*
     &    UStop2(1,2)**2 - 
     &   8*M_3*AtC*Cd(47)**2*Cd(9952)*UStop2(1,1)*
     &    UStop2(1,2)**2 + 
     &   2*At*AtC*Cd(47)**2*Cd(9952)*UStop2(1,1)*UStop2(1,2)**2

        Opt(8) = -(M_3*MTy2*Cd(4757)) + 
     &   2*MStop2(3)*(-(MTy*Cd(4770)*UCStop(1,3)*UStop2(1,2)) + 
     &      M_3*MStop2(2)*
     &       (-(MTy2**2*Cd(4758)) + 
     &         MGl2*MStop2(2)*
     &          (MGl2*(MGl2 - 3*MStop2(2)) + 3*MStop2(2)**2) + 
     &         MTy2**3*(Cd(4759) - MTy2*UStop2(1,2)**2)))

        Opt(21) = (MTy2**2*Cd(8366) + 
     &     (2*MTy2**6 + MGl2**2*
     &         (MGl2*(2*MGl2**3 - 4*Cd(1272)*Cd(5444)) + Cd(8361))-
     &          2*MGl2*Cd(1272)*Cd(5397)*MStop2(1)*MStop2(2) + 
     &        Cd(523)*MStop2(1)**2*MStop2(2)**2)*MStop2(3))*
     &   UStop(1,1)

	Cd(9940) = Cd(9932)*Cr(34)*MStop2(2) - 2*Cd(9939)*MStop2(3)

        Cd(9941) = 3/(64.D0*Pi)*
     &   (Alfa2L*(Cd(9940) + 
     &        8*Cd(402)*Cd(403)*Cr(13)*MStop2(1)*MStop2(2)))/
     &    (SB2*SW2*MStop2(2)*MStop2(3)**2)

	Cd(9944) = Cd(9939)*Cr(34) - 2*Cd(9943)*MStop2(3)

	Cd(9945) = 1/2.D0*(A0delStop(2)*Cd(9932)) - Cd(9944)*MStop2(3)

        Cd(9946) = 3/(64.D0*Pi)*
     &   (Alfa2L*(Cd(9945) + 
     &        4*A0delStop(1)*Cd(402)*Cd(403)*MStop2(2)))/
     &    (SB2*SW2*MStop2(2))

        Cd(9965) = 2*Cd(50)*Cd(9947)*Cd(9964)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(9961)*MGlpTmSt2(2)

        Cd(9971) = Cd(50)*Cd(9947)*Cd(9969)*Cr(28)*
     &    MGlpTmSt2(1) - Cd(9949)*Cd(9950)*Cd(9970)*MGlpTmSt2(2)

        Cd(9973) = Cd(50)*Cd(9947)*Cd(9972)*MGlpTmSt2(1) + 
     &   Cd(9949)*Cd(9950)*Cd(9969)*MGlpTmSt2(2)

	Cd(9974) = At*Cd(9971) - 2*AtC*Cd(9973)*UCStop(1,3)**2

	Cd(9980) = Cd(50)**2*Cd(9951) + Cd(47)**2*Cd(9952)

        Cd(9985) = Cd(50)*Cd(9947)*MGlpTmSt2(1) - 
     &   Cd(47)*Cd(9949)*MGlpTmSt2(2)

	Cd(9986) = At*Cd(9984) + 2*AtC*Cd(9985)*UCStop(1,3)**2

        Cd(9989) = Cd(50)*Cd(9947)*MGlpTmSt2(1) - 
     &   Cd(47)*Cd(9949)*Cd(9972)*MGlpTmSt2(2)

	Cd(9990) = AtC*Cd(9984) + 2*At*Cd(9989)*UCStopC(1,3)**2

        Cd(10000) = 2*Cd(50)*Cd(9947)*Cd(9996)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(9999)*MGlpTmSt2(2)

        Cd(10003) = Cd(50)*Cd(9947)*Cd(9972)*MGlpTmSt2(1) + 
     &   Cd(9949)*Cd(9950)*MGlpTmSt2(2)

	Cd(10004) = AtC*Cd(9971) - 2*At*Cd(10003)*UCStopC(1,3)**2

        Cd(10009) = Cd(47)*Cd(9949)*Cr(75)*MGlpTmSt2(2) + 
     &   Cd(50)*Cd(9947)*MGlpTmSt2(1)*UStop2(1,1)

	Cd(10012) = Cd(9951) - Cd(9953)

	Cd(10013) = Cd(9952) - Cd(9954)

        Cd(10022) = Cd(50)*Cd(9947)*Cd(10019)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(10021)*MGlpTmSt2(2)

        Cd(10023) = Cd(50)*Cd(9947)*Cd(10018)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(10020)*MGlpTmSt2(2)

	Cd(10031) = -(Cd(50)**2*Cd(9951)) + Cd(47)**2*Cd(9952)

	Cd(995) = AtC*Cd(967) - 2*MTy*Cd(989)*MStop2(3)*UCStopC(1,3)

        Cd(10044) = Cd(50)*Cd(9947)*Cd(10042)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(10043)*MGlpTmSt2(2)

        Cd(10045) = At*Cd(10044) + 
     &   Cd(50)*Cd(9947)*Cd(9955)*MGlpTmSt2(1)*MStop2(1)*
     &    UStop2(1,2)

	Cd(10046) = AtC*Cd(10023)*Cr(5) + Cd(10045)*UCStopC(1,3)

	Cd(10062) = 4*Cd(164)*Cd(10056) - Cd(10061)*MGlpTmSt2(1)

        Cd(10071) = Cd(50)*Cd(9947)*Cd(9969)*Cr(28)*
     &    MGlpTmSt2(1) + Cd(9949)*Cd(9950)*Cd(10070)*MGlpTmSt2(2)

        Cd(10072) = (MUE2*
     &     (2*Cd(47)**2*Cd(10071) - 
     &       2*Cd(50)*Cd(9950)*Cd(9951)*Cd(9981)*MStop2(3)*
     &        UStop2(1,1)))/(Cd(47)**2*Cd(50)*Cd(9950)*MStop2(3))

        Cd(10073) = Cd(50)*Cd(9947)*MGlpTmSt2(1) - 
     &   Cd(47)*Cd(9949)*Cd(9955)*MGlpTmSt2(2)

        Cd(10074) = (2*Cd(47)**2*Cd(10003) + 
     &     Cd(50)*Cd(9950)*Cd(9951)*MStop2(3))/
     &   (Cd(47)**2*Cd(50)*Cd(9950)*MStop2(3))

        Cd(10076) = Cd(50)*Cd(9947)*MGlpTmSt2(1) + 
     &   Cd(9949)*Cd(9950)*MGlpTmSt2(2)

        Cd(10078) = (-2*Cd(47)**2*Cd(10076)*Cd(10077) + 
     &     Cd(50)*Cd(9950)*Cd(9951)*MStop2(3))/
     &   (Cd(47)**2*Cd(50)*Cd(9950)*MStop2(3))

        Cd(10080) = Cd(50)*Cd(9947)*Cd(9969)*Cr(28)*
     &    MGlpTmSt2(1) + Cd(47)*Cd(9949)*Cd(9970)*MGlpTmSt2(2)

        Cd(10081) = Cd(50)*Cd(9947)*Cd(10077)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(9969)*MGlpTmSt2(2)

        Cd(10089) = Cd(50)*Cd(9947)*Cd(10088)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(10087)*MGlpTmSt2(2)

        Cd(10093) = Cd(50)*Cd(9947)*Cd(10091)*MGlpTmSt2(1) - 
     &   Cd(47)*Cd(9949)*Cd(10092)*MGlpTmSt2(2)

        Cd(10097) = Cd(50)*Cd(9947)*Cd(10096)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*Cd(10095)*MGlpTmSt2(2)

	Cd(10104) = Cd(50)**2*Cd(9953) + Cd(47)**2*Cd(9954)

	Cd(10107) = Cd(50)**2*Cd(10101) - 3*Cd(9980)*UStop2(1,1)**2

	Cd(10113) = Cd(251)*Cd(10103) - Cd(50)*Cd(10112)*MGlpTmSt2(2)

        Cd(10115) = -2*Cd(50)*Cd(9949)*Cd(10113) + 
     &   Cd(9952)*Cd(10103)*Cd(10114)

        Cd(10116) = Cd(50)*Cd(9947)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*MGlpTmSt2(2)

	Cd(10117) = 2*Cd(10086)*Cd(10116) + AtC*Cd(10073)*UCStop(1,3)**2

        Cd(10118) = Cd(47)*Cd(9949)*MGlpTmSt2(2) + 
     &   2*Cd(50)*Cd(9947)*MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2)

	Cd(10119) = 2*Cd(10039)*Cd(10116) + AtC*Cd(10118)

        Cd(10123) = Cd(50)*Cd(9947)*Cd(10121) + 
     &   Cd(47)*Cd(9949)*Cd(10122)

        Cd(10125) = Cd(50)*Cd(9947)*Cd(10077)*MGlpTmSt2(1) + 
     &   Cd(47)*Cd(9949)*MGlpTmSt2(2)

        Cd(10128) = Cd(50)*Cd(9947)*Cd(10077)*MGlpTmSt2(1) - 
     &   Cd(9949)*Cd(9950)*MGlpTmSt2(2)

        Cd(10131) = Cd(9950)*Cd(9951)*Cd(10114) - 
     &   2*Cd(47)*Cd(9947)*Cd(10130)

        Cd(10144) = Cd(50)*Cd(9947)*Cd(10143)*MGlpTmSt2(1) + 
     &   2*Cd(47)*Cd(9949)*Cr(5)*Cr(75)*MGlpTmSt2(2)

	Cd(10158) = Cd(10102) - 3*Cd(9952)*UStop2(1,1)**2

        Cd(10161) = (Cd(50)*Cd(10160)*MStop2(3) + 
     &     6*Cd(47)**2*Cd(9949)*Cd(9955)*MGlpTmSt2(2)*
     &      UStop2(1,1)**2)/(Cd(47)**2*Cd(50)*MStop2(3))

        Cd(10175) = -(Cd(50)*Cd(9947)*Cd(10171)) + 
     &   Cd(47)*Cd(9949)*Cd(10174)

        Cd(10184) = -(Cd(50)*Cd(9947)*Cd(10179)) + 
     &   Cd(47)*Cd(9949)*Cd(10183)

        Cd(10236) = (Cd(10232) - Cd(152)*Cd(2597)*Cd(10235))/
     &   (Cd(47)**2*MStop2(1)**2)

	Cd(1016) = Cd(1014) + 2*Cd(1015)*UStop2(1,1)

        Cd(10258) = Cd(50)**2*Cd(9981) - 
     &   Cd(179)*Cd(2681)*MStop2(2)*UStop2(1,2)

	Cd(10259) = Cd(50)**2 - Cd(179)*Cd(2681)*MStop2(2)

        Cd(10269) = Cd(179)*Cd(2681)*Cd(9981)*MStop2(2) - 
     &   Cd(50)**2*UStop2(1,2)

	Cd(10302) = 4*Cd(47)*Cd(10301)*MStop2(1) + Cd(10296)*MStop2(3)

        Cd(10335) = (4*(-2*Cd(47)*Cd(10334) + 
     &       Cd(10327)*MStop2(3)))/(Cd(47)**2*MStop2(3))

        Cd(10374) = (2*Cd(47)**2*Cd(10370) + 
     &     Cd(10373)*MStop2(3))/(Cd(47)**2*MStop2(3))

	Cd(1029) = -2*Cd(958) + Cd(1024)

        Cd(10380) = (2*Cd(47)**2*Cd(10377) + 
     &     Cd(10379)*MStop2(3))/(Cd(47)**2*MStop2(3))

	Cd(1030) = -4*Cd(960) + Cd(1028)

	Cd(10392) = AtC*Cd(10387) + At*Cd(10391)

	Cd(10405) = MGlmT2*Cd(10399) - Cd(10404)*MStop2(1)

	Cd(10412) = Cd(13)*Cd(10408) + 2*Cd(10280)*Cd(10411)

	Cd(10417) = 2*Cd(10416) + Cd(10413)*MStop2(3)*UCStop(1,3)

	Cd(10424) = MUE*Cd(10417) + MUEC*Cd(10423)

	Cd(10427) = -2*Cd(10426) + Cd(10025)*MStop2(3)*UCStop(1,3)

	Cd(10431) = MUE*Cd(10427) - MUEC*Cd(10430)

	Cd(10432) = MGlmT2*Cd(10431) - Cd(10424)*MStop2(1)

        Cd(10464) = (Cd(10442) + Cd(10463))/
     &   (Cd(47)**2*MStop2(1)**2*MStop2(3))

	Cd(10469) = Cd(47)*Cd(492) - MGlmT2*Cd(10468)

        Cd(10470) = (-8*Cd(47)*Cd(471) - 8*Cd(10469)*MStop2(3))/
     &   (Cd(47)*MStop2(3))

	Cd(10498) = M_3C*Cd(10497) + Cd(10495)*UCStopC(1,3)

	Cd(10504) = M_3C*Cd(10499) - Cd(10503)*UCStopC(1,3)

	Cd(10505) = MGlmT2*Cd(10498) + Cd(10504)*MStop2(1)

	Cd(10525) = AtC*Cd(10524) - M_3C*Cd(152)*MStop2(3)

	Cd(10528) = At*Cd(10524) - M_3*Cd(152)*MStop2(3)

	Cd(10544) = MGlmT2*Cd(10399) - Cd(10543)*MStop2(1)

        Cd(10559) = Cd(164)*Cd(10558)*MStop2(3) + 
     &   Cd(13)*Cd(47)*Cd(152)*Cd(9955)*UStop2(1,1)

	Cd(10560) = 2*Cd(10559) + Cd(13)*Cd(152)*Cd(164)*MStop2(3)

        Cd(10568) = 2*AtC*Cd(47)*Cd(10553)*Cr(5) + 
     &   Cd(10567)*UCStopC(1,3)

	Cd(10609) = Cd(50)*Cr(28) + MGlmT2*Cd(9970)*MStop2(1)

	Cd(10610) = Cd(50) - MGlmT2*Cd(9969)*MStop2(1)

	Cd(10611) = Cd(50) - MGlmT2*MStop2(1)

	Cd(10624) = MGlmT2*Cd(10623) + Cd(10622)*MStop2(2)

	Cd(10625) = 4*Cd(50)*Cd(10624)*MStop2(2) + Cd(10621)*MStop2(3)

        Cd(10626) = (-2*MGlmT2*Cd(50)*Cd(9955)*Cr(75) + 
     &     Cd(10602)*MStop2(3))/MStop2(3)

	Cd(10627) = Cd(10604) + Cd(10626)*MStop2(2)

	Cd(10628) = Cd(50)*Cd(10306) + Cd(10602)

	Cd(10633) = -2*MGlmT2*Cd(50)*Cd(10025) + Cd(10603)

        Cd(10636) = -2*Cd(50)**2*UCStop(1,3) + 
     &   2*At*Cd(10604)*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(10637) = -2*Cd(50)**2 + At*Cd(10039)*Cd(10604)

        Cd(10640) = (Cd(10602)*MStop2(3) + 
     &     3*MGlmT2*Cd(50)*Cd(9955)*UStop2(1,1)**2)/MStop2(3)

	Cd(10642) = 13*Cd(50) + 4*Cd(10602)

        Cd(10656) = -(MGlmT2*Cd(50)*Cd(10655)) + 
     &   4*Cd(10015)*Cd(10602)*MStop2(3)*UStop2(1,1)

	Cd(10683) = MUEC*Cd(10682) + 2*MUE*Cd(10680)*UCStop(1,3)

	Cd(10684) = MGlmT2*Cd(10683) + Cd(10676)*MStop2(2)

	Cd(10695) = MUE*Cd(10682) + 2*MUEC*Cd(10694)*UCStopC(1,3)

	Cd(10700) = MGlmT2*Cd(10695) + Cd(10699)*MStop2(2)

	Cd(10726) = MGlmT2*Cd(10725) + Cd(10719)*MStop2(2)

	Cd(10730) = Cd(50) + Cd(9948)*Cd(10729)*MStop2(1)

	Cd(10731) = 4*MGlmT2*Cd(50) + Cd(10604)

	Cd(10732) = Cd(50)*Cd(10730) - Cd(9948)*Cd(10731)*MStop2(3)

	Cd(10735) = Cd(50)*MStop2(3) - Cd(10734)*MStop2(2)*UStop2(1,1)

	Cd(10737) = 3*MGlmT2*Cd(50) + Cd(10604)

	Cd(10738) = Cd(50)*Cd(179)*MStop2(1) + Cd(10737)*MStop2(3)

	Cd(10739) = Cd(50)*Cd(10736) - Cd(10738)*UStop2(1,1)

        Cd(10750) = -4*Cd(50)*Cd(179)*Cd(10749)*MStop2(2) + 
     &   2*MUE2*Cd(10603)*MStop2(3)*UStop2(1,1)**2

	Cd(10757) = -(MGlmT2*Cd(10756)) + Cd(10751)*MStop2(2)

	Cd(10777) = Cd(179)*Cd(251) + Cd(50)*Cd(10306)

        Cd(10778) = (-2*MGlmT2*Cd(50)*Cd(9955)*Cr(75) + 
     &     Cd(179)*Cd(251)*MStop2(3))/MStop2(3)

        Cd(10789) = (Cd(251)*Cd(10788)*MStop2(3) + 
     &     12*MGlmT2*Cd(50)*Cd(9955)*UStop2(1,1)**2)/MStop2(3)

	Cd(10790) = (Cd(10789) + 4*Cd(251)*Cd(9969)*MStop2(2))/MStop2(2)

	Cd(10791) = 5*MGlmT2*Cd(50) - 2*Cd(179)*Cd(251)*Cd(10015)

	Cd(10795) = AtC*Cd(9970) + 4*Cd(10794)

        Cd(10799) = (8*(Cd(10798) + 
     &       12*Cd(9955)*MStop2(1)*UStop2(1,1)**2))/Cd(50)

	Cd(10809) = Cd(50) + MStop2(1)*MStop2(2)*UStop2(1,1)

	Cd(10810) = AtC*Cd(10809) - M_3C*MStop2(2)*MStop2(3)

	Cd(10813) = MGlmT2*Cd(10812) + Cd(10810)*UCStop(1,3)**2

	Cd(10817) = Cd(50) + Cd(10816)*MStop2(2)

	Cd(10818) = At*Cd(10817) - M_3*MStop2(2)*MStop2(3)

	Cd(10821) = MGlmT2*Cd(10820) + Cd(10818)*UCStopC(1,3)**2

	Cd(10827) = At*Cd(10825) + 2*AtC*Cd(10826)*UCStop(1,3)**2

	Cd(10829) = AtC*Cd(10825) + 2*At*Cd(10828)*UCStopC(1,3)**2

	Cd(10830) = MUE*Cd(10827) + MUEC*Cd(10829)

	Cd(105) = -(Cd(50)*Cd(95)) + Cd(47)*Cd(104)

        Cd(10848) = Cd(50)*Cd(179)*Cr(75)**2 - 
     &   2*Cd(251)*MStop2(2)*MStop2(3)*UStop2(1,2)

        Cd(10849) = Cd(50)*Cd(179)*Cd(9948)*Cd(9983) + 
     &   Cd(10848)*UStop2(1,1)

	Cd(1077) = Cd(1076)*UStop2(1,1) + 2*Cd(1073)*UStop2(1,2)

	Cd(10858) = At*Cd(10849) - 2*Cd(10857)

	Cd(10869) = AtC*Cd(10849) - 2*Cd(10868)

	Cd(10873) = Cd(50) - MGlmT2*MStop2(2)

	Cd(10874) = MGlmT2*Cd(774)*Cd(9948) + Cd(10873)*UStop2(1,1)

	Cd(10875) = Cd(50) - MGlmT2*Cd(814)

	Cd(10876) = Cd(50) - MGlmT2*Cd(652)

	Cd(10877) = MGlmT2*Cd(9948)*MStop2(3) - Cd(10876)*UStop2(1,1)

        Cd(10890) = Cd(50)*Cd(10889) + 
     &   2*MTy2*MUE2*Cd(48)*Cd(251)*MStop2(3)*UStop2(1,1)**2

	Cd(1083) = -Cd(967) + Cd(1039)

	Cd(106) = Cd(85)*MStop2(3)*UStop2(1,2) - Cd(105)*UStopC(1,1)

	Cd(1086) = -2*Cd(984) + Cd(1068)

	Cd(1087) = 2*MTy*Cd(1085) + AtC*Cd(1086)*UCStop(1,3)

	Cd(1088) = Cd(989) - Cd(1077)

	Cd(1089) = 2*Cd(993) - Cd(1082)

        Cd(10972) = 4*MUE2*Cd(10970)*MStop2(1)*UStop2(1,1)*
     &    UStopC(2,2) - MUEC**2*Cd(10971)*UUStopC(1,1)

	Cd(10973) = 2*Cd(10945) - Cd(10972)

	Cd(11004) = Cd(11003)*UStop2(1,2) + Cd(11002)*UUStopC(2,2)

        Cd(11005) = -2*Cd(11004)*MStop2(1)**2 + 
     &   6*Cd(9794)*MStop2(3)**2*UStop(1,1)

	Cd(11007) = MUEC*Cd(10999) + MUE*Cd(11005)

	Cd(11045) = Cd(11044)*MStop2(3) - 2*Cd(11042)*UStopC(1,2)

	Cd(11096) = -(Cd(11095)*MStop2(3)) + 2*Cd(11093)*UStopC(1,2)

	Cd(11100) = Cd(11096) + Cd(11099)

	Cd(11148) = Cd(11141) + Cd(11147)

        Cd(11174) = 4*Cd(11173)*MStop2(3)*UStopC(1,1) - 
     &   2*Cd(11007)*UUStopC(1,3)

	Cd(1110) = Cd(1108)*UStop2(1,1) + Cd(1109)*UStop2(1,2)

        Cd(11190) = Cd(11189)*MStop2(1)**2 + 
     &   6*Cd(9794)*MStop2(3)**2*UStop(1,2)

	Cd(11192) = Cd(11190) - 2*Cd(11191)*MStop2(3)

	Cd(11194) = AtC*Cd(11034) + 2*Cd(11048)

	Cd(11209) = Cd(11074) + AtC*Cd(11099)

        Cd(11210) = 4*Cd(11209)*UStopC(1,2) - 
     &   AtC**2*Cd(11087)*UUStopC(1,2)

	Cd(11211) = 2*Cd(11131) + AtC*Cd(11147)

        Cd(11212) = Cd(11211)*UStopC(1,1) + 
     &   4*AtC**2*Cd(11112)*UCStop(1,3)*UStop2(1,1)*UStopC(1,2)

	Cd(1116) = MUEC*Cd(1096) + 2*MUE*Cd(1105)*UCStop(1,3)**2

	Cd(1117) = MUE*Cd(1100)*UCStop(1,3) - MUEC*Cd(1110)*UCStopC(1,3)

	Cd(1118) = MUE*Cd(1096) - MUEC*Cd(1115)*UCStopC(1,3)**2

	Cd(11259) = Cd(11251) + 4*Cd(9794)*Cd(11248)*UStop2(1,1)

        Cd(11313) = (-4*Cd(11292)*MStop2(2) + 
     &     Cd(11312)*MStop2(3)**2)/MStop2(3)**2

        Cd(11314) = (Cd(11313)*MStop2(3)**2 - 
     &     4*Cd(11278)*MStop2(2)**2*UUStopC(1,2))/
     &   (MStop2(2)**2*MStop2(3)**2)

        Cd(11368) = Cd(11360)*MStop2(2)*UUStopC(1,1)*
     &    UUStopC(2,3) - Cd(11364)*UUStopC(3,1)

	Cd(11373) = Cd(11347) + Cd(9809)*Cd(11372)*MStop2(3)*UStop2(1,2)

	Cd(11419) = Cd(11418) + 2*Cd(11397)*MStop2(3)

	Cd(11420) = Cd(11380)*UUStopC(1,1) + Cd(11384)*UUStopC(1,2)

        Cd(11432) = Cd(11431)*MStop2(3) - 
     &   2*Cd(11288)*MStop2(2)*UUStopC(3,1)

	Cd(11441) = 2*Cd(11294) - Cd(11439)

        Cd(11442) = Cd(11440)*MStop2(3)*UCStopC(1,3) + 
     &   2*Cd(11441)*MStop2(2)*UStopC(1,1)

	Cd(11550) = 2*Cd(11397) + Cd(11533)

        Cd(11562) = Cd(9948)*Cd(11560)*MStop2(3)*UStop2(1,1) - 
     &   2*Cd(11561)*MStop2(2)*UUStopC(3,3)

	Cd(11604) = Cd(11590) - Cd(9809)*Cd(11603)*UStop2(1,2)

	Cd(11645) = MUEC*Cd(11639) + MUE*Cd(11644)*UCStop(1,3)**2

	Cd(11650) = MUE*Cd(11639) + MUEC*Cd(11649)*UCStopC(1,3)**2

        Cd(11676) = 4*Cd(11675)*MStop2(3)*UStopC(1,1) - 
     &   Cd(11674)*UUStopC(1,3)

	Cd(11686) = Cd(11664) + Cd(11676)

        Cd(11687) = (Cd(11686) - 
     &     Cd(11685)*MStop2(2)*MStop2(3)**2*UUStopC(1,1))/
     &   (MStop2(2)*MStop2(3)**2)

        Cd(11688) = (Cd(11687)*MStop2(3) + 
     &     4*MUEC*Cd(11652)*UCStopC(1,3)*UStop2(1,1)*UStopC(1,1))/
     &   (MStop2(1)*MStop2(3))

	Cd(11689) = SA2*Cd(11636) + CA2*Cd(11688)

        Cd(11704) = (4*MTy*Cd(13)*Cd(11647) + Cd(11703))/
     &   (MStop2(1)*MStop2(2)*MStop2(3))

	Cd(11733) = (Cd(11732) + 4*Cd(9969)*MStop2(1))/MStop2(1)

	Cd(11768) = M_3C*Cd(11767) + Cd(11766)*UCStopC(1,3)

	Cd(11799) = -2*Cd(11798) + At*Cd(11790)*Cr(28)

	Cd(11806) = -2*Cd(11805) + AtC*Cd(11790)*Cr(28)

        Cd(11814) = 2*AtC*Cd(164)*Cd(10018)*Cr(5) + 
     &   Cd(11813)*UCStopC(1,3)

	Cd(11821) = MUEC*Cd(11814) + MUE*Cd(11820)

	Cd(11830) = -(MGlmT2*Cd(474)) + Cd(11829)*MStop2(1)

        Cd(11831) = -(Cd(47)*Cd(11830)) + 
     &   MGl2*Cd(152)*Cd(1201)*MStop2(3)

	Cd(11835) = Cd(50) - Cd(138)**2*Cd(9969)

        Cd(11836) = MGl2*Cd(179)*MStop2(3) + 
     &   3*MGlmT2*Cd(50)*Cd(9955)*UStop2(1,1)**2

	Cd(11837) = Cd(50) - MGl2*Cd(179)

	Cd(11838) = Cd(50) + MGl2*Cd(179)

	Cd(11840) = AtC*Cd(11837)*Cd(11838) + Cd(11839)*MStop2(2)

	Cd(11842) = At*Cd(11837)*Cd(11838) + Cd(11841)*MStop2(2)

	Cd(11843) = Cd(11840)*UCStop(1,3) + Cd(11842)*UCStopC(1,3)

	Cd(11849) = (Cd(11848) + 4*Cd(9969)*MStop2(2))/MStop2(2)

	Cd(11853) = Cd(50)*Cr(28) + MGlmT2*Cd(10070)*MStop2(1)

	Cd(11854) = Cd(50) + MGlmT2*MStop2(1)

	Cd(11855) = Cd(50) + MGlmT2*Cd(9972)*MStop2(1)

	Cd(11862) = 8*(Cd(11859) - Cd(11861))

        Cd(11885) = (8*(Cd(11884) + 
     &       12*Cd(9955)*MStop2(1)*UStop2(1,1)**2))/Cd(50)

	Cd(11895) = Cd(251)*Cd(10803) + Cd(50)*Cr(28)

	Cd(11898) = Cd(50) + Cd(251)*MStop2(1)*UStop2(1,1)

	Cd(11899) = AtC*Cd(11898) - M_3C*Cd(251)*MStop2(3)

        Cd(11900) = 2*Cd(138)*Cd(10169)*MStop2(2) - 
     &   Cd(11899)*UCStop(1,3)**2

	Cd(11903) = Cd(50) + Cd(251)*Cd(10816)

	Cd(11904) = At*Cd(11903) - M_3*Cd(251)*MStop2(3)

        Cd(11905) = 2*Cd(138)*Cd(10176)*MStop2(2) - 
     &   Cd(11904)*UCStopC(1,3)**2

	Cd(11914) = At*Cd(9970) + 2*Cd(11913)

        Cd(11922) = -(Cd(50)*Cd(251)*Cd(10696)) + 
     &   2*Cd(179)*Cd(10388)*MStop2(2)*MStop2(3)

	Cd(11929) = AtC*Cd(10020)*Cr(5) + Cd(11928)*UCStopC(1,3)

	Cd(11930) = MUE*Cd(11927) + MUEC*Cd(11929)

	Cd(11932) = 2*Cd(11931) + Cd(10548)*MStop2(3)*UCStop(1,3)

	Cd(11936) = MUE*Cd(11932) + MUEC*Cd(11935)

	Cd(11943) = Cd(50) + Cd(251)*Cd(11942)

	Cd(11944) = Cd(50) + Cd(251)*Cd(9972)*MStop2(1)

        Cd(11951) = Cd(50)*Cd(11950) - 
     &   2*MGl2**2*MUE2*Cd(179)*MStop2(3)*UStop2(1,1)**2

	Cd(11954) = MTy2*Cd(10751) + 2*MTy*Cd(11953)

	Cd(11960) = 2*Cd(50) - Cd(138)**2

	Cd(11961) = -(Cd(50)*Cd(118)**2) + Cd(47)*Cd(11960)

	Cd(11962) = Cd(50)**2*MGlpTmSt2(1) + Cd(47)**2*MGlpTmSt2(2)

        Cd(11963) = (Cd(10193)*Cd(11962)*MStop2(3) + 
     &     6*Cd(47)**2*Cd(50)*Cd(9955)*UStop2(1,1)**2)/
     &   (Cd(47)**2*MStop2(3))

	Cd(11964) = -(Cd(50)*Cd(118)) + Cd(47)*Cd(138)

        Cd(11966) = Cd(50)**2*Cd(152)*MGlpTmSt2(1) - 
     &   Cd(47)**2*Cd(179)*MGlpTmSt2(2)

        Cd(11967) = Cd(50)**2*MGlpTmSt2(1)*MStop2(1) - 
     &   Cd(47)**2*MGlpTmSt2(2)*MStop2(2)

	Cd(1190) = Cd(1188)*UStop2(1,1) + Cd(1189)*UStop2(1,2)

	Cd(11974) = Cd(50)*Cd(11972) - Cd(11973)*MGlpTmSt2(2)

	Cd(11977) = Cd(50)*Cr(28) + Cd(47)*UStop2(1,1)**2

	Cd(11980) = Cd(47) - Cd(50)*Cd(9972)

	Cd(11984) = Cd(47) - Cd(50)

	Cd(11985) = Cd(50) - Cd(47)*Cd(9972)

	Cd(11987) = Cd(50) - 2*At*AtC*MGlpTmSt2(2)*UStop2(1,1)**2

	Cd(11992) = -(Cd(50)*Cd(11989)) + Cd(11991)*MGlpTmSt2(2)

	Cd(11997) = Cd(47)*Cd(6096) + Cd(50)*Cd(152)*MGlpTmSt2(1)

	Cd(11998) = 2*Cd(47) + Cd(50)

	Cd(12000) = 3*Cd(47) + 2*Cd(50)

	Cd(12002) = Cd(50)*Cd(11997) + Cd(47)**2*Cd(179)*MGlpTmSt2(2)

        Cd(12006) = 2*Cd(12002)*Cd(12005) + 
     &   At*Cd(47)*Cd(50)*Cd(12004)*UCStopC(1,3)**2

        Cd(12011) = Cd(50)*Cd(12010) - 
     &   2*At*AtC*MGlpTmSt2(2)*UStop2(1,1)**2

	Cd(12022) = Cd(50)*Cd(12020) - Cd(12021)*MGlpTmSt2(2)

	Cd(12033) = Cd(47)*Cd(12028) + Cd(50)*Cd(12032)

	Cd(12038) = -(Cd(50)*Cd(12034)) + Cd(47)*Cd(12037)

	Cd(12039) = AtC*Cd(11983) + 2*Cd(12038)

	Cd(12042) = Cd(47)*Cr(75) + Cd(50)*UStop2(1,1)

        Cd(1197) = MUE*Cd(1181)*MStop2(3)*UCStop(1,3) + 
     &   MUEC*Cd(1190)*UCStopC(1,3)

	Cd(12045) = Cd(50)**2*MGlpTmSt2(1) - Cd(47)**2*MGlpTmSt2(2)

	Cd(12050) = Cd(50)*Cd(10019) + Cd(47)*Cd(10021)

	Cd(12051) = Cd(50)*Cd(10018) + Cd(47)*Cd(10020)

	Cd(12052) = Cd(50)*Cd(10042) + Cd(47)*Cd(10043)

	Cd(12066) = Cd(50)*Cd(12065) + Cd(10782)*MGlpTmSt2(2)

        Cd(12089) = (Cd(13)**2*MGlpTmSt2(2)*MStop2(3)*
     &      UCStopC(1,3) + 2*Cd(50)*Cd(12087)*UCStopC(3,3))/
     &   MStop2(3)

        Cd(12091) = Cd(47)*Cd(12085)*UStop2(1,2)*UStopC(1,2) + 
     &   Cd(50)*Cd(12086)*Cr(91)*UUStopC(1,3)

        Cd(12093) = Cd(9948)*MGlpTmSt2(2)*MStop2(3) + 
     &   2*Cd(50)*UStop2(1,1)

        Cd(12096) = (2*Cd(47)*UCStopC(1,3)*UStop2(1,2) + 
     &     2*Cd(50)*Cr(91)*UUStopC(3,1))/(Cd(47)*Cd(50))

        Cd(12097) = (-(Cd(47)**2*Cd(12096)) + 
     &     MGlpTmSt2(1)*MStop2(3)*UCStopC(1,3)*UStop2(1,1))/
     &   (Cd(47)**2*MStop2(3))

        Cd(12098) = (2*Cd(50)*UCStopC(1,3)*UStop2(1,1) - 
     &     2*Cd(47)*Cr(91)*UUStopC(3,1))/(Cd(47)*Cd(50))

        Cd(12099) = (Cd(50)**2*Cd(12098) + 
     &     MGlpTmSt2(2)*MStop2(3)*UCStopC(1,3)*UStop2(1,2))/
     &   (Cd(50)**2*MStop2(3))

        Cd(12100) = (Cd(50)*UStop2(1,1) - Cd(47)*UStop2(1,2))/
     &   (Cd(47)*Cd(50))

	Cd(12101) = Cd(12099)*UStop2(1,1) - Cd(12097)*UStop2(1,2)

	Cd(12102) = (-6*Cd(50) + MGlpTmSt2(2)*MStop2(3))/MStop2(3)

        Cd(12103) = (-(Cd(47)**2*Cd(12102)) + 
     &     Cd(50)**2*MGlpTmSt2(1))/(Cd(47)**2*Cd(50)**2)

        Cd(12104) = Cd(12093)*UCStop(1,3)*UStop2(1,1) + 
     &   2*Cd(50)**2*Cd(12100)*Cr(91)*UUStopC(3,2)

        Cd(12112) = Cd(47)*Cd(9948)*UStop2(1,1)**2 - 
     &   Cd(12111)*UStop2(1,2)**2

        Cd(12113) = MStop2(2)**2*MStop2(3) + 
     &   2*Cd(50)*Cd(1272)*UStop2(1,1)

	Cd(12118) = 2*Cd(50)*Cd(12117) - MStop2(2)**2*MStop2(3)

        Cd(12119) = Cd(47)**2*Cd(12118) + 
     &   Cd(50)**2*Cd(6458)*MStop2(3)*UStop2(1,1)

	Cd(12120) = Cd(47)**2 + Cd(50)**2*UStop2(1,2)

	Cd(12127) = Cd(50)*Cd(12124) - Cd(12126)*MStop2(2)*MStop2(3)

        Cd(12133) = (2*Cd(47)*Cd(12081)*MGlpTmSt2(2)*
     &      UStopC(1,1) + 2*Cd(50)*Cd(12078)*MGlpTmSt2(1)*
     &      UUStopC(1,3))/(Cd(47)*Cd(50))

	Cd(12167) = MUE*Cd(12166) - 2*MUEC*Cd(12157)*UStop(1,2)

        Cd(12171) = Cd(12170)*MStop2(3)*UStop(2,1) + 
     &   2*MTy*Cd(50)*Cd(9809)*UStop2(1,1)

        Cd(12178) = Cd(12177)*UStop(1,2) - 
     &   AtC*MTy*Cd(47)*Cr(9)*UUStop(1,1)

	Cd(12181) = -(Cd(50)*Cd(152)) + Cd(47)*Cd(12180)

        Cd(12182) = 2*MTy2*Cd(47)*UStop2(1,1)*UUStop(3,3) + 
     &   Cd(12181)*UUStop(3,4)

	Cd(12183) = (2*MTy*Cd(47)*Cd(12179) + Cd(12182))/(Cd(47)*Cd(50))

        Cd(12184) = (2*Cd(12178)*UStopC(1,1) - 
     &     Cd(47)**2*Cd(12183)*UStopC(2,1))/Cd(47)**2

        Cd(12185) = (Cd(47)**2*Cd(12184) - 
     &     MTy*Cd(9794)*MGlpTmSt2(1)*MStop2(3)*UCStopC(1,3)*
     &      UStop2(1,1))/(Cd(47)**2*MStop2(3))

	Cd(12188) = Cd(12187)*MStop2(3) - 2*AtC*MTy*Cd(47)*UCStop(1,3)

        Cd(12189) = Cd(12188)*UStop(1,1) + 
     &   2*MTy2*Cd(47)*Cr(91)*UStopC(2,2)

        Cd(12191) = (Cd(12189)*UStop(2,2) + 
     &     2*MTy*Cd(47)*Cd(12190)*UStopC(2,1))/Cd(47)

        Cd(12192) = (Cd(9948)*Cd(12170)*MStop2(3)*UStop(2,2) + 
     &     Cd(50)*Cd(12191)*UStopC(1,1))/MStop2(3)

        Cd(12193) = (-(Cd(47)*Cd(12192)*MStop2(3)) + 
     &     2*MTy*Cd(50)**2*Cd(10987)*UStop2(1,1)*UStopC(1,1))/
     &   (Cd(47)*Cd(50)**2*MStop2(3))

	Cd(12200) = MGl2*Cd(12199) + Cd(44)*Cd(12197)*UStop(1,1)

	Cd(12207) = MGl2*Cd(12206) + Cd(44)*Cd(12204)*UStop(1,1)

	Cd(12208) = -(MUEC*Cd(10991)*Cd(12202)) + MUE*Cd(12207)

	Cd(12211) = -(Cd(50)*Cd(152)) + Cd(47)*Cd(12210)

        Cd(12213) = (-(MTy*Cd(47)*Cd(12212)) + 
     &     Cd(12211)*UStop(1,1))/(Cd(47)*Cd(50))

        Cd(12214) = (-(Cd(47)*Cd(12213)*UStop(2,1)) + 
     &     MTy*Cd(698)*UStop2(1,1))/Cd(47)

	Cd(12215) = Cd(12185)*UStopC(1,2) + Cd(12193)*UStopC(2,2)

        Cd(12216) = (Cd(12171)*UCStop(1,3)*UStopC(2,2) + 
     &     2*Cd(50)**2*Cd(12214)*UUStopC(3,2))/Cd(50)**2

        Cd(12221) = (MTy*Cd(47)*Cd(12220)*UStopC(2,2) + 
     &     MUEC*Cd(12187)*MStop2(3)*UUStop(3,2))/Cd(47)

	Cd(12222) = Cd(12221) - MTy*MUEC*Cd(12218)*UCStopC(1,3)

        Cd(12223) = (MUE*Cd(12189)*UStop2(1,1) + 
     &     2*Cd(47)*Cd(12222)*UStopC(2,1))/Cd(47)

	Cd(12224) = MTy*Cd(12202)*UStop(1,1) + At*Cd(1985)*UStop(1,2)

        Cd(12227) = MUEC*Cd(47)*Cd(12226)*UCStopC(1,3) - 
     &   MUE*Cd(12224)*UStop(1,1)*UStop2(1,1)

	Cd(12238) = Cd(44)*Cd(12174) + MGl2*Cd(12237)

	Cd(12239) = MUE*Cd(12238) - 2*MTy*Cd(47)*Cd(12234)*UCStopC(1,3)

        Cd(12240) = AtC*MTy*Cd(47)*Cd(10946)*Cr(9) - 
     &   Cd(12239)*UStop(1,2)

        Cd(12244) = 2*MTy2*Cd(47)*Cd(11161)*UStop2(1,1) + 
     &   Cd(12187)*MStop2(3)*UUStop(3,4)

        Cd(12245) = (Cd(12244) - 
     &     2*AtC*MTy*Cd(47)*Cr(91)*UCStop(1,3))/Cd(47)

	Cd(12246) = MUE*Cd(12245) + 2*MTy*Cd(12243)*UStopC(2,1)

	Cd(12250) = Cd(12223)*UStopC(1,1) + Cd(12246)*UUStopC(3,4)

	Cd(12265) = Cd(50)*Cd(118)**2 - Cd(47)*Cd(138)**2

	Cd(12269) = Cd(47)*Cd(12268) + 2*MGl2*Cd(50)*MStop2(1)

	Cd(12270) = Cd(9948)*Cd(12267) - Cd(12269)*UStop2(1,2)

        Cd(12271) = -(MUEC*Cd(12270)*UCStopC(1,3)) + 
     &   MUE*Cd(12265)*UCStop(1,3)*UStop2(1,1)

	Cd(12275) = MTy2*Cd(50) + Cd(47)*Cd(12274)

	Cd(12281) = MTy2*Cd(50) + Cd(47)*Cd(12280)

	Cd(12282) = Cd(47)*Cd(9948)*Cd(12278) - Cd(12281)*UStop2(1,2)

	Cd(12327) = Cd(12326)*UCStopC(1,3)**2 + Cd(12324)*UStop2(1,1)

	Cd(12328) = MUE*Cd(12322) + MUEC*Cd(12327)

        Cd(12336) = (Cd(1201)*Cd(12323) - 
     &     Cd(164)*Cd(2597)*Cd(12335))/(Cd(47)**2*MStop2(1)**2)

	Cd(12338) = 2*Cd(4187) - Cd(4190)*Cd(9969)*MGlpTmSt2(2)

	Cd(12339) = Cd(50)**2 - Cd(138)*Cd(4190)*Cd(9969)*MGlpTmSt2(2)

	Cd(12340) = Cd(50)**2 - MGl2*Cd(4187)

	Cd(12341) = AtC*Cd(4187) + 2*M_3C*Cd(4190)*MGlpTmSt2(2)

	Cd(12342) = AtC*Cd(12340) + Cd(12341)*MStop2(2)

	Cd(12343) = At*Cd(4187) + 2*M_3*Cd(4190)*MGlpTmSt2(2)

	Cd(12344) = At*Cd(12340) + Cd(12343)*MStop2(2)

	Cd(12345) = 4*(Cd(12342)*UCStop(1,3) + Cd(12344)*UCStopC(1,3))

	Cd(12346) = AtC*Cd(251)*Cd(2681) + 4*M_3C*Cd(4187)

	Cd(12347) = AtC*Cd(251)*Cd(2681) + 8*M_3C*Cd(4187)

	Cd(12356) = Cd(50)**2 - Cd(138)*Cd(4187)

	Cd(12357) = Cd(50)**2 + Cd(251)*Cd(2681)*MStop2(2)

        Cd(12358) = 4*Cd(4187)*Cd(10224)*MStop2(2) + 
     &   Cd(424)*Cd(12357)*UStop2(1,1)

        Cd(12367) = (Cd(1201)*Cd(12357) - 
     &     Cd(251)*Cd(2681)*Cd(12366))/(Cd(50)**2*MStop2(2)**2)

	Cd(12387) = 4*M_3*AtC*Cd(12385) + At*Cd(12386)

	Cd(12388) = 36*Cd(47)**2 + Cd(12387)*UStop2(1,2)

	Cd(12390) = 4*AtC*Cd(47)**2 + Cd(12389)*MStop2(3)

        Cd(12421) = (-2*Cd(47)**2*Cd(10370) + 
     &     Cd(12420)*MStop2(3))/(Cd(47)**2*MStop2(3))

        Cd(12430) = (-2*Cd(47)**2*Cd(10377) + 
     &     Cd(12429)*MStop2(3))/(Cd(47)**2*MStop2(3))

	Cd(12444) = MUE*Cd(12437) + MUEC*Cd(12443)

	Cd(12450) = MUEC*Cd(12448) + MUE*Cd(12449)

	Cd(12453) = 2*Cd(12452) + Cd(10025)*MStop2(3)*UCStop(1,3)

	Cd(12459) = MUE*Cd(12453) + MUEC*Cd(12458)

	Cd(12460) = MGlmT2*Cd(12459) + 2*Cd(12450)*MStop2(1)

        Cd(1240) = (Cd(1239) + 
     &     4*MTy*Cd(585)*Cd(1238)*MStop2(3)*MStop2(4))/
     &   (MStop2(1)*MStop2(3)**2)

        Cd(12479) = (Cd(12464) + Cd(12478))/
     &   (Cd(47)**2*MStop2(1)**2*MStop2(3))

	Cd(12482) = MGl2*Cd(179)**2 + Cd(4187)

	Cd(12484) = Cd(50)*Cd(138) + Cd(12481)

        Cd(12485) = (12*Cd(50)**2 - 
     &     4*Cd(138)*Cd(9969)*Cd(12484))/Cd(50)**2

        Cd(12486) = Cd(12483)*MStop2(3) + 
     &   12*MGlmT2*Cd(50)*Cd(9955)*MStop2(2)*UStop2(1,1)**2

	Cd(12487) = 2*Cd(50)**2 - MGl2*Cd(12482)

	Cd(12488) = Cd(50)*Cd(138) + Cd(12482)

	Cd(12489) = 3*Cd(50)*Cd(138) + Cd(12481)

	Cd(12490) = AtC*Cd(12488) + M_3C*Cd(12489)

	Cd(12491) = AtC*Cd(12487) + Cd(12490)*MStop2(2)

	Cd(12492) = At*Cd(12488) + M_3*Cd(12489)

	Cd(12493) = At*Cd(12487) + Cd(12492)*MStop2(2)

	Cd(12494) = 8*(Cd(12491)*UCStop(1,3) + Cd(12493)*UCStopC(1,3))

	Cd(12497) = MGlmT2*Cd(50) + 4*Cd(12482)

        Cd(12498) = (-(MGlmT2*Cd(50)*Cd(11863)) + 
     &     Cd(10112)*Cd(12497)*MStop2(3))/MStop2(3)

	Cd(12500) = MGlmT2*Cd(50) + 2*Cd(12482)

	Cd(12502) = Cd(50) + MGlmT2*Cd(11942)

	Cd(1243) = 2*CA2*Cd(1240) + S2A*Cd(1242)

	Cd(12503) = -4*Cd(50)*Cd(12502) + Cd(12483)*MStop2(3)

	Cd(12504) = AtC*Cd(12503) + 8*M_3C*Cd(12500)*MStop2(3)

	Cd(12505) = -4*Cd(50)*Cd(11855) + Cd(12483)*MStop2(3)

	Cd(12508) = AtC*Cd(12507) + 4*Cd(10039)*Cd(12500)*MStop2(3)

        Cd(12511) = 2*Cd(50)**2*UCStop(1,3) + 
     &   At*Cd(12483)*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(12512) = 4*Cd(50)**2 + At*Cd(10039)*Cd(12483)

        Cd(1244) = 1/(128.D0*Pi**3)*
     &   (Alfas2L*(4*SA2*Cd(1220) + 
     &        Cd(1243)*MStop2(1)*MStop2(3)**2))/
     &    (SB2*SW2*MStop2(1)*MStop2(3)**2)

	Cd(12520) = 8*Cd(50)**2 + Cd(12400)*Cd(12483)

	Cd(12522) = 8*MGlmT2*Cd(50) + Cd(12483)

	Cd(12529) = 2*MGlmT2*Cd(50)*Cd(9958) + Cd(12483)*UStop2(1,2)**2

	Cd(12533) = Cd(50)*Cr(28) - MGlmT2*Cd(9970)*MStop2(1)

	Cd(12534) = Cd(50) + MGlmT2*Cd(9969)*MStop2(1)

	Cd(12535) = At*Cd(12533) - 2*AtC*Cd(12534)*UCStop(1,3)**2

	Cd(12541) = AtC*Cd(12539) + At*Cd(12540)

	Cd(12545) = AtC*Cd(12533) - 2*At*Cd(11854)*UCStopC(1,3)**2

        Cd(12554) = -4*Cd(50)*Cd(9955)*Cr(75)*MStop2(2)**2 + 
     &   Cd(12483)*MStop2(3)

        Cd(12555) = -4*Cd(50)*Cd(251)*Cd(9955)*Cr(75)*
     &    MStop2(2) + Cd(12483)*MStop2(3)

	Cd(12558) = Cd(10033)*Cd(12483) + 2*Cd(50)**2*UCStop(1,3)

	Cd(12566) = 2*Cd(50)**2 + At*Cd(10039)*Cd(12483)

        Cd(12567) = Cd(12566)*UCStopC(1,3) + 
     &   AtC*Cd(12483)*Cr(5)*UStop2(1,1)*UStop2(1,2)

        Cd(12576) = 2*Cd(50)*Cd(11950)*MStop2(2) - 
     &   MGl2*MUE2*Cd(12483)*MStop2(3)*UStop2(1,1)**2

        Cd(12581) = -(Cd(50)**2*Cd(118)*Cd(164)) + 
     &   Cd(47)**2*Cd(138)*Cd(251)

	Cd(12599) = Cd(50)**2*Cd(164) + Cd(47)**2*Cd(251)

        Cd(12608) = (Cd(13)**2*Cd(47)**2*MGlpTmSt2(2) - 
     &     2*MUE2*Cd(50)**2*MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12610) = (Cd(12609)*MGlpTmSt2(1) - 
     &     Cd(47)**2*Cd(12608)*UCStopC(1,3))/Cd(47)**2

        Cd(12613) = (Cd(50)**2*Cd(10585)*MGlpTmSt2(1) + 
     &     Cd(47)**2*Cd(12612)*MGlpTmSt2(2)*UStopC(2,2))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12614) = Cd(47)**2*MStop2(2)**2 + 
     &   Cd(50)**2*MStop2(1)**2*UStop2(1,2)

        Cd(12615) = Cd(47)**2*MStop2(2)**2 + 
     &   Cd(50)**2*MGlpTmSt4(1)*UStop2(1,1)

        Cd(12616) = MUE**2*Cd(12614)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(12615)*UCStopC(1,3)**2

        Cd(12635) = 2*MTy2*Cd(167) - 
     &   MTy*Cd(10985)*MGlpTmSt2(1)*UStopC(1,2)

        Cd(12638) = 2*MTy2*Cd(167) - 
     &   MTy*Cd(9794)*MGlpTmSt2(1)*UStopC(1,1)

        Cd(12642) = MTy*Cd(12641)*UStop2(1,1) - 
     &   2*MTy2*MUE*Cd(167)*UUStop(1,1)

        Cd(12643) = (Cd(47)**2*Cd(9948)*Cd(12639) - 
     &     Cd(50)**2*Cd(12638)*UStop2(1,2))/(Cd(47)**2*Cd(50)**2)

        Cd(12644) = MTy*MUE*Cd(10985)*MGlpTmSt2(1)*
     &    UCStop(1,3)**2*UStopC(1,2) + Cd(12642)*UUStopC(1,2)

        Cd(1258) = (Cd(1257) + 
     &     8*MTy2*Cd(1252)*MStop2(3)*UStop2(1,1))/
     &   (MStop2(1)*MStop2(3))

	Cd(12662) = Cd(50)*Cd(118)**2 - Cd(47)*Cd(11960)

	Cd(12663) = 2*M_3*Cd(12002) + At*Cd(47)*Cd(50)*Cd(12004)

	Cd(12668) = Cd(50)*Cd(12667) - Cd(11973)*MGlpTmSt2(2)

	Cd(1260) = -4*SA2*Cd(1250) + S2A*Cd(1259)

        Cd(1261) = 1/(128.D0*Pi**3)*
     &   (Alfas2L*(Cd(1260) + 2*CA2*Cd(1258)*MStop2(1)*MStop2(3)))/
     &    (SB2*SW2*MStop2(1)*MStop2(3))

        Cd(12681) = MUEC*Cd(12270)*UCStopC(1,3) - 
     &   MUE*Cd(12265)*UCStop(1,3)*UStop2(1,1)

	Cd(12683) = Cd(47)*Cd(9948)*Cd(12682) - Cd(12281)*UStop2(1,2)

        Cd(12743) = 3/(64.D0*Pi)*
     &   (Alfa2L*(4*A0delStop(1)*Cd(12707)*Cd(12708)*MStop2(2) + 
     &        Cd(12742)*MStop2(3)**2))/
     &    (SW2*TB2*MStop2(2)*MStop2(3)**2)

	Cd(12752) = 8*Cd(16)*Cd(12748) - 4*Cd(12751)*MStop2(3)

        Cd(12774) = Cd(12769)*MStop2(2) - 
     &   2*Cd(12700)*Cd(12770)*UUStopC(2,3)

	Cd(12792) = Cd(17)*Cd(12791) - 2*Cd(12761)*MStop2(3)

	Cd(12799) = Cd(22)*Cd(12795) + 4*Cd(12798)*MStop2(3)

        Cd(12817) = Cd(12816)*UStop2(1,1) - 
     &   4*Cd(12809)*MStop2(1)*MStop2(3)*UStopC(1,1)

        Cd(12884) = Cd(22)*Cd(12883) - 
     &   4*Cd(12876)*Cr(13)*MStop2(1)*MStop2(3)

	Cd(126) = 2*Cd(117) + Cd(125)*UStop2(1,1)

	Cd(12922) = Cd(12701)*Cd(12906) + Cd(12897)*UUStopC(3,2)

        Cd(12953) = Cd(16)*Cd(12952) - 
     &   4*Cd(12948)*Cr(9)*Cr(34)*MStop2(3)

	Cd(12972) = M_3*Cd(12965) + M_3C*Cd(12971)

	Cd(12973) = Cd(12959) - 2*Cd(12972)

        Cd(12991) = 1/(64.D0*Pi**4)*
     &   (-Cd(12973) + 64*Pi**4*
     &       (Cd(12990)*MGlpTmStxGlT4(1)*MStop2(1)))/
     &    (MGlpTmStxGlT4(1)*MStop2(1))

        Cd(12992) = 1/(64.D0*Pi**4)*
     &   (64*Pi**4*(Cd(12991)*MGlpTmStxGlT4(1)**2*MStop2(1)) + 
     &      MTy2*Cd(164)*Cd(12836)**2*Cr(45)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(1)**2*MStop2(1)*MStop2(3))

        Cd(13069) = 1/(64.D0*Pi**4)*
     &   (2*MTy2*Cd(13068) + Cd(13065)*MGlpTmSt2(2)*MStop2(3))/MTy2

	Cd(13071) = (MGl*Cd(13069) + Cd(13070))/(MGl*TB2*MStop2(3))

        Cd(13072) = 1/Pi**4*
     &   (Cd(13046) - 64*Pi**4*(TB2*Cd(13071)*MGlpTmStxGlT4(2)))/
     &    TB2

        Cd(13091) = -(Cd(251)*Cr(28)*MGlpTmSt2(2)) + 
     &   2*Cd(96)*Cd(254)*UStop2(1,1)*UStop2(1,2)

	Cd(13092) = MTy2*Cd(13091) - MGlpTmSt4(2)*MStop2(2)

	Cd(13093) = Cd(251)*MStop2(3) - Cd(96)*Cd(254)*UStop2(1,2)

        Cd(13094) = Cd(251)*Cd(13021)*UStopC(1,1) + 
     &   Cd(96)*Cd(254)*Cd(12983)*UStopC(1,2)

        Cd(13095) = Cd(13094)*UUStopC(2,1) - 
     &   Cd(96)*Cd(254)*Cd(12984)*UCStop(3,1)*UUStopC(3,1)

        Cd(13096) = Cd(13095) + 
     &   Cd(96)*Cd(254)*Cd(12700)*Cd(13019)*UStop2(1,1)*UStopC(1,1)

        Cd(13097) = Cd(12700)*Cd(13093)*UCStop(1,3)**2 - 
     &   Cd(251)*Cd(13018)*UStop2(1,1)*UStopC(1,2)

        Cd(13149) = 1/(256.D0*Pi**4)*
     &   (Pi**4*(Cd(12700)*Cd(13148)*MStop2(1)**2*MStop2(3)) - 
     &      2*Cd(13144)*UStop2(1,1)*UStop2(1,2)*UUStopC(1,1))/
     &    (MStop2(1)**2*MStop2(3))

        Cd(13150) = Cd(12700)**2*Cd(13142)*UCStop(1,3) + 
     &   Cd(12701)**2*Cd(13138)*UStop(1,2)

	Cd(13152) = Cd(13132) + Cd(13119)*Cr(13)

        Cd(13223) = 1/(256.D0*Pi**4)*
     &   (Cd(13222) + 256*Pi**4*
     &       (Cd(12701)*Cd(13221)*MStop2(1)*MStop2(2)*MStop2(3)))/
     &    (MStop2(1)*MStop2(2)*MStop2(3))

	Cd(13231) = 2*M_3*Cd(13230) + Cd(13225)*MGlpTmSt2(1)

        Cd(13232) = 1/(64.D0*Pi**4)*
     &   (MGl*Cd(13231) + 64*Pi**4*(M_3C*Cd(13228)))/
     &    (MGl*TB2*MStop2(3))

        Cd(13233) = 1/(128.D0*Pi**4)*
     &   (Cd(1985)*Cd(12836)**2*Cr(45) - 
     &      128*Pi**4*(TB2*Cd(13232)*MGlpTmStxGlT4(1)))/
     &    (TB2*MGlpTmStxGlT4(1)**2)

	Cd(13272) = M_3C*Cd(13244) + M_3*Cd(13250)

	Cd(13288) = 2*Cd(13287) + Cd(12992)*Cr(13)

        Cd(13292) = Cd(13013)*MGlpTmSt2(2) + 
     &   Cd(179)*Cd(12701)**2*UCStopC(1,3)**2

        Cd(13301) = (MStop2(3)*
     &     (Cd(13300)*MStop2(2) - 
     &       Cd(13299)*UStopC(1,1)*UUStopC(2,1)))/MStop2(2)**2

        Cd(13302) = 1/(128.D0*Pi**4)*
     &   (Cd(13301)*MStop2(2) + 
     &      Cd(12700)*Cr(91)*MStop2(3)*UStop2(1,1)*UUStopC(3,3))/
     &    (MStop2(2)*MStop2(3))

        Cd(13318) = 1/(512.D0*Pi**4)*
     &   (-(Cd(13317)*MGlpTmStxGlT4(2)**2) + 
     &      Cd(179)*Cd(2681)*Cd(12999)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(2)**2*MStop2(3))

        Cd(13364) = 1/(512.D0*Pi**4)*
     &   (Cd(13317)*MGlpTmStxGlT4(2)**2 + 
     &      Cd(251)*Cd(2681)*Cd(12999)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(2)**2*MStop2(3))

        Cd(13410) = 1/(512.D0*Pi**4)*
     &   (Pi**4*(Cd(12700)*Cd(13403)*MStop2(2)**2) - 
     &      4*Cd(13409)*UUStopC(2,1))/MStop2(2)**2

	Cd(13411) = Cd(12701)**2*Cd(13390) + Cd(12700)**2*Cd(13396)

	Cd(13436) = Cd(13378) + Cd(13435)*Cr(34)

	Cd(13456) = M_3*Cd(13447) + M_3C*Cd(13448)

	Cd(13457) = Cd(13455) - 2*Cd(13456)

	Cd(13458) = M_3*Cd(13250) + M_3C*Cd(13438)

        Cd(13459) = 1/(64.D0*Pi**4)*
     &   (Cd(13457) + 64*Pi**4*(Cd(13445)*MGlpTmStxGlT4(2)))/
     &    MGlpTmStxGlT4(2)

	Cd(13462) = M_3C*Cd(13448) + M_3*Cd(13461)

	Cd(13463) = Cd(13455) - 2*Cd(13462)

        Cd(13464) = 1/(128.D0*Pi**4)*
     &   (-2*Cd(13463)*MGlpTmStxGlT4(2) + 
     &      Cd(2008)*Cd(12999)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(2)**2*MStop2(3))

	Cd(1341) = -(Cd(424)*Cd(1332)) + 2*Cd(1340)

	Cd(13477) = -(Cd(13476)*MSbot2(1)) + MStop2(1)**2*UStop2(1,2)

        Cd(13482) = 2*Cd(13472)*Cd(13477) + 
     &   Cd(13478)*Cd(13481)*MStop2(1)

        Cd(13511) = Cd(22)*Cd(13510) - 
     &   2*Cd(13516)*Cr(13)*MSq2Diff(1,1)*MStop2(1)

        Cd(13597) = -(3/(4.D0*Pi)*
     &     (Alfa2L*(Cd(13596) - 
     &          8*MTy2*MSq2Diff(1,1)*MSq2Diff(1,2)))/
     &      (MTy*SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2)))

        Cd(13605) = Cd(16)*Cd(13603) + 
     &   2*Cd(13604)*Cr(34)*MSq2Diff(1,2)*MStop2(2)

        Cd(13657) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(A0delStop(2)*Cd(13656) + 
     &        A0delq(MSbot2(1),Q)*Cd(13891)*Cr(9)*MStop2(3)))/
     &    (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2)*MStop2(3))

        Cd(13658) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(Cd(13656)*MStop2(2) + 
     &        Cd(13891)*Cr(9)*MSbot2(1)*MStop2(3)))/
     &    (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2)*MStop2(3))

        Cd(13659) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(Cd(13656)*Cr(34)*MStop2(2) + 
     &        Cd(13466)*Cd(13891)*Cr(9)*MSbot2(1)*MStop2(3)))/
     &    (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2)*MStop2(3))

	Cd(1364) = Cd(1289) - Cd(1356)

	Cd(1365) = Cd(1293) - Cd(1363)

        Cd(13709) = -(Cd(13708)*UCStop(1,3)) + 
     &   MTy*Cd(179)*Cd(12701)*UStop2(1,1)*UStop2(1,2)

        Cd(13757) = AtC*Cd(13734)*MStop2(2)**2*UStop2(1,1)**2 + 
     &   Cd(13713)*UStopC(1,1)

	Cd(13823) = Cd(13811)*MSq2Diff(1,1)**2 - Cd(13822)*MSq2Diff(1,2)

        Cd(13876) = 1/2.D0*(A0delq(MSbot2(1),Q)*Cd(13823)) - 
     &   Cd(13466)*Cd(13875)*MSbot2(1)*MSq2Diff(1,2)

	Cd(13880) = Cd(13466)*Cd(13823) - 2*Cd(13875)*MSq2Diff(1,2)

        Cd(13886) = A0delq(MSbot2(1),Q)*Cd(13885) - 
     &   A0delStop(2)*Cd(13891)*MSq2Diff(1,1)

        Cd(13887) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(Cd(13885)*MSbot2(1) - 
     &        Cd(13891)*MSq2Diff(1,1)*MStop2(2)))/
     &    (SW2*TB2*MSq2Diff(1,1)**2*MSq2Diff(1,2))

	Cd(1383) = 2*Cd(1304) + Cd(1381)

	end
