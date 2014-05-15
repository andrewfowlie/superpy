* TLhr_abbr1_9.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_9
	implicit none

#include "TLhr.h"


        Cd(4670) = 6*MTy2*(MGl2*(MGl2 - 3*MTy2) + 3*MTy2**2) + 
     &   MStop2(2)*(-Cd(4668) + Cd(4669)*MStop2(2)) + Opt(7)

        Cd(4733) = 4*MTy2*(MGl2*(MGl2 - 3*MTy2) + 3*MTy2**2) + 
     &   MStop2(2)*(Cd(4731) - Cd(4732)*MStop2(2)) - Opt(7)

        Cd(4925) = MTy2*(MStop2(1)**3 + 
     &      MTy2*(MTy2*(4*MTy2 - 8*MStop2(1)) + 3*MStop2(1)**2)) + 
     &   3*MTy2*Cd(2193)*MStop2(2)**2 + Cd(1592)*MStop2(2)**3 - 
     &   MStop2(2)*(Cd(4923) + 2*Cd(4924)*MStop2(2)) + 
     &   Cd(1755)*MStop2(3)**3

        Cd(4963) = 5*MTy2**4 + 22*MTy2**3*MStop2(1) - 
     &   (-2*MTy2**2 + Cd(4962))*MStop2(1)**2 - 
     &   12*MTy2*MStop2(1)**3 - 
     &   2*MStop2(2)*(Cd(4960) + 
     &      MStop2(2)*(Cd(4961) + 3*Cd(44)*MStop2(2))) + 
     &   Cd(607)*MStop2(3)**3

        Cd(5039) = 20*MTy2**4 + 4*MTy2**3*MStop2(1) - 
     &   2*MTy2**2*MStop2(1)**2 - 32*MTy2*MStop2(1)**3 + 
     &   4*Cd(5036)*MStop2(2) - 
     &   2*MStop2(2)**2*(3*Cd(5037) + Cd(5038) - 
     &      2*Cd(67)*MStop2(2)) + 9*Cd(1755)*MStop2(3)**3

        Cd(5092) = 15*MTy2**4 - 40*MTy2**3*MStop2(1) - 
     &   16*MTy2**2*MStop2(1)**2 + 4*Cd(5089)*MStop2(2) + 
     &   (2*Cd(5090) - Cd(5091))*MStop2(2)**2 + 
     &   7*Cd(1755)*MStop2(3)**3 - 
     &   6*MTy2*(2*MStop2(2)**3 - 18*MStop2(2)**2*MStop2(3) + 
     &      3*MStop2(3)**2*(MStop2(3) + 3*MStop2(4)))

        Cd(5326) = 330*MTy2**4 + 100*MTy2**3*MStop2(1) + 
     &   62*MTy2**2*MStop2(1)**2 + 44*MTy2*MStop2(1)**3 + 
     &   4*Cd(5323)*MStop2(2) + 
     &   2*MStop2(2)**2*(-Cd(4924) + Cd(5324) + 
     &      2*Cd(5325)*MStop2(2)) + Cd(1755)*MStop2(3)**3

        Cd(5496) = 8*MTy2**4 - 16*MTy2**3*MStop2(1) + 
     &   10*MTy2**2*MStop2(1)**2 - 7*MTy2*MStop2(1)**3 - 
     &   Cd(5493)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(4924) + Cd(5494) - 
     &      Cd(5495)*MStop2(2)) - Cd(1755)*MStop2(3)**3

        Cd(5537) = 15*MTy2**4 + 52*MTy2**3*MStop2(1) + 
     &   10*MTy2**2*MStop2(1)**2 - 
     &   MStop2(2)**2*(2*Cd(5535) + Cd(5536) + 
     &      12*Cd(2763)*MStop2(2)) + 
     &   4*MTy2*(6*MStop2(1)**3 - Cd(5534)*MStop2(2)) - 
     &   4*Cd(1755)*MStop2(3)**3

        Cd(5764) = 1/4.D0*(MStop2(4)**4 - 
     &      MStop2(3)**3*(MStop2(3) + 4*MStop2(4))) - 
     &   Cd(5763)*MStop2(2) + 
     &   MTy2*(-MStop2(1)**3 + 
     &      MTy2*(-3*MStop2(1)**2 + MTy2*(MTy2 + MStop2(1))) + 
     &      6*Cd(1735)*MStop2(2)**2)

        Cd(5926) = 43*MTy2**4 - 8*MTy2**3*MStop2(1) + 
     &   6*MTy2**2*MStop2(1)**2 - 5*MTy2*MStop2(1)**3 - 
     &   Cd(5924)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(4924) + Cd(5925) - 
     &      5*Cd(44)*MStop2(2)) - Cd(1755)*MStop2(3)**3

        Cd(6031) = 5*MTy2**4 + 20*MTy2**3*MStop2(1) + 
     &   2*MTy2**2*MStop2(1)**2 + 24*MTy2*MStop2(1)**3 - 
     &   MStop2(2)*(4*Cd(6027) + 
     &      MStop2(2)*(2*Cd(6028) + Cd(6030) + 
     &         4*Cd(6029)*MStop2(2))) - 6*Cd(1755)*MStop2(3)**3

        Cd(6219) = 47*MTy2**4 - 6*MTy2**3*MStop2(1) + 
     &   2*MTy2**2*MStop2(1)**2 - 6*MTy2*MStop2(1)**3 - 
     &   Cd(6216)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(4924) + Cd(6217) - 
     &      Cd(6218)*MStop2(2)) - Cd(1755)*MStop2(3)**3

        Cd(6296) = -(Cd(2417)*
     &      (MStop2(3)**3*(-4*MStop2(2) + MStop2(3)) + 
     &        2*MStop2(2)**2*
     &         (MStop2(2)*(MStop2(2) - 2*MStop2(3)) + 
     &           3*MStop2(3)**2))) + 
     &   2*MStop2(1)*MStop2(2)*
     &    (-2*MStop2(2)*(Cd(6289)*MStop2(1) + 
     &         Cd(6290)*MStop2(2)) + Cd(6290)*MStop2(3)*MStop2(4))+
     &     2*MTy2*(Cd(6295)*MStop2(4) + 
     &      MTy2*(Cd(6292) + 
     &         MTy2*(5*MTy2 + 2*MStop2(4))*UStop2(1,1)))

        Cd(6469) = 15*MTy2**4 + 26*MTy2**3*MStop2(1) + 
     &   10*MTy2**2*MStop2(1)**2 - 2*MTy2*MStop2(1)**3 - 
     &   2*Cd(6466)*MStop2(2) - 
     &   MStop2(2)**2*(2*Cd(6467) + Cd(6468) - 
     &      4*Cd(1735)*MStop2(2)) + 3*Cd(1755)*MStop2(3)**3

        Cd(6518) = 30*MTy2**4 + 32*MTy2**3*MStop2(1) - 
     &   10*MTy2**2*MStop2(1)**2 - 24*MTy2*MStop2(1)**3 - 
     &   2*MStop2(2)*(2*Cd(6515) + 
     &      MStop2(2)*(Cd(6516) - Cd(6517) + 2*Cd(5863)*MStop2(2)))
     &     + 9*Cd(1755)*MStop2(3)**3

        Cd(6563) = 15*MTy2**4 + 4*MTy2**3*MStop2(1) + 
     &   8*MTy2**2*MStop2(1)**2 - 4*MTy2*MStop2(1)**3 + 
     &   4*Cd(2092)*Cd(6560)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(6561) + Cd(6562) + 
     &      4*Cd(3718)*MStop2(2)) + 3*Cd(1755)*MStop2(3)**3

        Cd(6612) = 15*MTy2**4 - 24*MTy2**3*MStop2(1) - 
     &   2*MTy2**2*MStop2(1)**2 - 8*MTy2*MStop2(1)**3 + 
     &   4*Cd(6609)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(6610) + Cd(6611) + 
     &      12*Cd(2092)*MStop2(2)) + 2*Cd(1755)*MStop2(3)**3

        Cd(6637) = 2*MTy2**3*Cd(3318) - 
     &   MTy2*(MStop2(1)**3 + MTy2*(6*MTy2**2 + 3*MStop2(1)**2)) - 
     &   3*Cd(6636)*MStop2(1)*MStop2(2) + 
     &   (-3*MTy2*Cd(4362) + 2*Cd(4924))*MStop2(2)**2 + 
     &   Cd(2221)*MStop2(2)**3 - Cd(1755)*MStop2(3)**3

        Cd(6688) = -15*MTy2**4 + 4*MTy2**3*Cd(286) - 
     &   11*MTy2**2*MStop2(1)**2 - 9*MTy2*MStop2(1)**3 + 
     &   MStop2(2)*(Cd(6685)*MStop2(1) - 
     &      MStop2(2)*(Cd(6686) - Cd(6687) + Cd(3718)*MStop2(2)))-
     &     Cd(1755)*MStop2(3)**3

        Cd(6694) = 2*MTy2**3*Cd(607) - 
     &   MTy2*(MStop2(1)**3 + MTy2*(2*MTy2**2 + 3*MStop2(1)**2)) - 
     &   3*Cd(6693)*MStop2(1)*MStop2(2) + 
     &   (3*MTy2*Cd(45) + 2*Cd(4924))*MStop2(2)**2 - 
     &   Cd(2200)*MStop2(2)**3 - Cd(1755)*MStop2(3)**3

        Cd(6758) = 15*MTy2**4 - 78*MTy2**3*MStop2(1) - 
     &   38*MTy2**2*MStop2(1)**2 - 30*MTy2*MStop2(1)**3 + 
     &   2*Cd(6756)*MStop2(2) + 
     &   MStop2(2)**2*(Cd(6611) + 2*Cd(6757) + 
     &      4*Cd(3846)*MStop2(2)) + 2*Cd(1755)*MStop2(3)**3

        Cd(6832) = 2*(15*MTy2**4 - 7*MTy2**3*MStop2(1) + 
     &      (MTy2**2 + Cd(6831))*MStop2(1)**2 - 
     &      3*MTy2*MStop2(1)**3 + 
     &      MStop2(2)*(-Cd(6830) + 
     &         MStop2(2)*(Cd(2080) + 2*Cd(155)*MStop2(2)))) + 
     &   Cd(607)*MStop2(3)**3

        Cd(6875) = 45*MTy2**4 - 20*MTy2**3*MStop2(1) + 
     &   14*MTy2**2*MStop2(1)**2 - 6*MTy2*MStop2(1)**3 - 
     &   2*Cd(6873)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(1532)*Cd(2193) + Cd(6874) + 
     &      2*Cd(4928)*MStop2(2)) - 5*Cd(1755)*MStop2(3)**3

        Cd(7178) = 5*MTy2**4 + 6*MTy2**3*MStop2(1) - 
     &   MTy2**2*MStop2(1)**2 - MTy2*MStop2(1)**3 - 
     &   MStop2(2)*(Cd(7176) + 
     &      MStop2(2)*(Cd(6562) + Cd(7177) + 3*Cd(6529)*MStop2(2)))
     &     - 3*Cd(1755)*MStop2(3)**3

        Cd(7467) = 30*MTy2**4 - 12*MTy2**3*MStop2(1) + 
     &   17*MTy2**2*MStop2(1)**2 + 9*MTy2*MStop2(1)**3 - 
     &   Cd(4926)*MStop2(2) + 
     &   MStop2(2)**2*(-2*Cd(4924) + Cd(4927) + 
     &      Cd(4928)*MStop2(2)) + Cd(1755)*MStop2(3)**3

        Cd(7491) = 2*MTy2*
     &    (9*MStop2(1)**3 + 
     &      MTy2*(51*MStop2(1)**2 + MTy2*(81*MTy2 + 20*MStop2(1))))
     &     + 2*MStop2(2)*(Cd(7488) + 
     &      MStop2(2)*(-Cd(4924) + 3*Cd(7489) + Cd(7490)*MStop2(2))
     &      ) + Cd(1755)*MStop2(3)**3

        Cd(7564) = 15*MTy2**4 + 8*MTy2**3*MStop2(1) + 
     &   12*MTy2**2*MStop2(1)**2 + 2*MTy2*MStop2(1)**3 - 
     &   2*Cd(7561)*MStop2(2) + 
     &   MStop2(2)**2*(Cd(6030) - 2*Cd(7562) + 
     &      2*Cd(7563)*MStop2(2)) + 6*Cd(1755)*MStop2(3)**3

        Cd(7908) = MGl2**2*Cd(1952) - 2*MGl2*MTy2*Cd(6168) + 
     &   3*MTy2*(-2*MStop2(1)**2 + MTy2*(5*MTy2 + 7*MStop2(1))) + 
     &   MStop2(2)*(-(Cd(2193)*Cd(2200)) + Cd(7907)*MStop2(2)) + 
     &   MStop2(3)*(-4*Cd(832) + 3*MGl2*MStop2(4))

        Cd(7975) = Cd(782)*
     &    (MStop2(3)**3*(-4*MStop2(2) + MStop2(3)) + 
     &      2*MStop2(2)**2*
     &       (MStop2(2)*(MStop2(2) - 2*MStop2(3)) + 3*MStop2(3)**2)
     &      ) + 2*MStop2(1)*MStop2(2)*
     &    (6*Cd(172)*MStop2(1)*MStop2(2) + 
     &      2*Cd(7969)*MStop2(2)**2 - Cd(7969)*MStop2(3)*MStop2(4))
     &     + 2*MTy2*(Cd(7974)*MStop2(4) + 
     &      MTy2*(Cd(7971) - 
     &         MTy2*(5*MTy2 + 2*MStop2(4))*UStop2(1,1)))

        Cd(8041) = 3*MTy2*(MTy2 - MStop2(1))*
     &    (MTy2*(4*MTy2 - 4*MStop2(1)) + 3*MStop2(1)**2) - 
     &   Cd(8039)*MStop2(2) + 2*Cd(4924)*MStop2(2)**2 + 
     &   3*MTy2*Cd(8040)*MStop2(2)**2 - Cd(4520)*MStop2(2)**3 - 
     &   Cd(1755)*MStop2(3)**3

        Cd(8070) = 60*MTy2**5 + 3*MTy2**4*MStop2(1) + 
     &   4*MTy2**3*MStop2(1)**2 + 
     &   MTy2*(5*MStop2(1)**4 - Cd(8067)*MStop2(2)) - 
     &   MStop2(2)**2*(-2*Cd(8068) + 
     &      MStop2(2)*(-2*Cd(8069) + Cd(2065)*MStop2(2) + 
     &         MStop2(2)**2))

        Cd(8117) = 59*MTy2**4 + 1/4.D0*(Cd(3514)**2*Cd(8116)) - 
     &   26*MTy2**3*MStop2(1) + 40*MTy2**2*MStop2(1)**2 - 
     &   7*MTy2*MStop2(1)**3 + 
     &   MStop2(2)*(-Cd(8114) + 
     &      MStop2(2)*(Cd(8115) + 3*Cd(1684)*MStop2(2))) + 
     &   Cd(607)*MStop2(3)**3

        Cd(8160) = MTy2*(MStop2(1)**3 + 
     &      MTy2*(9*MStop2(1)**2 + MTy2*(15*MTy2 + 2*MStop2(1))))-
     &     Cd(8158)*MStop2(2) + 
     &   MStop2(2)**2*(Cd(6562) + Cd(8159) + Cd(7884)*MStop2(2)) + 
     &   3*Cd(1755)*MStop2(3)**3

        Cd(8611) = 2*(-2*Cd(8610)*MStop2(2) - MStop2(2)**4 + 
     &      MTy2*(-4*MStop2(1)**3 + 
     &         MTy2*(MTy2*(MTy2 - 2*MStop2(1)) + 5*MStop2(1)**2) + 
     &         Cd(155)*MStop2(2)**2)) + 
     &   6*MStop2(2)**2*MStop2(3)**2 - Cd(874)*MStop2(3)**3

	Cd(9281) = MGl2*(MGl2*(2*MGl2 - Cd(1591)) + Cd(1593)) - Cd(9272)

        Cd(10235) = MGl2*Cd(10233) + Opt(15) + 
     &   MGlpT2*MUEC**2*UCStopC(1,3)**2 - 
     &   2*MTy*Cd(10234)*UCStopC(1,3)*UStop2(1,1)*UStop2(1,2) + 
     &   MTy2*(MUE**2*UCStop(1,3)**2 + 
     &      2*MUE2*UStop2(1,1)*
     &       (Cd(9948)*UStop2(1,1) + UStop2(1,2)**2))

        Cd(10282) = MGlmT2*(Cr(28) + 2*Cd(9948)*UStop2(1,1)) + 
     &   AtC*(-2*M_3C*UCStop(1,3)**2 + 
     &      AtC*UCStop(1,3)**2*(-UStop2(1,1) - UStop2(1,2)) - 
     &      2*M_3*UStop2(1,1)*UStop2(1,2)) + 
     &   At*(-(Cd(10001)*UCStopC(1,3)**2) - 
     &      2*M_3C*UStop2(1,1)*UStop2(1,2) + 
     &      AtC*UStop2(1,1)*
     &       (-2*UStop2(1,1)*UStop2(1,2) - 2*UStop2(1,2)**2))

        Cd(11947) = MTy2*Cd(1201) + 
     &   MUEC**2*(MGl2 - MStop2(2))*UCStopC(1,3)**2 + 
     &   2*M_3*MTy*MUEC**2*UCStopC(1,3)**3 - 
     &   2*MUE2*Cd(11946)*UStop2(1,1)*UStop2(1,2) + 
     &   2*M_3C*MTy*MUEC**2*UCStopC(1,3)*UStop2(1,1)*
     &    UStop2(1,2) + MUE**2*UCStop(1,3)*
     &    (UCStop(1,3)*(MGl2 - MStop2(2) + 
     &         2*M_3C*MTy*UCStop(1,3)) + 
     &      2*M_3*MTy*UStop2(1,1)*UStop2(1,2))

        Cd(13913) = -2*M_3*MTy2*TB*UCStopC(1,3)**2 + 
     &   MTy*UCStopC(1,3)*(TB*MGlpTmSt2(1) - 
     &      2*M_3*MUE*TB2*UStop2(1,2)) + 
     &   UStop2(1,2)*(AtC*TB*Cd(13912) + 
     &      MUE*TB2*(MGlpT2 - 2*M_3C*MTy*UCStop(1,3) - 
     &         MStop2(1)*UStop2(1,2)))

        Cd(2042) = 2*(-13*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - 3*MTy2) + 9*MTy2**2)) - 
     &   (-5*MTy2**2 + MGl2*(MGl2 + 3*MTy2))*MStop2(1) - 
     &   MStop2(2)*(Cd(2040) + 
     &      MStop2(2)*(-3*Cd(2041) + 5*MStop2(2)))

        Cd(3329) = 4*MGlmT2**3 + 
     &   (MGl2*(5*MGl2 - MTy2) - MTy2**2)*MStop2(1) + 
     &   MStop2(2)*(-Cd(3327) + 
     &      MStop2(2)*(Cd(3328) + 11*MStop2(2)))

        Cd(4000) = 4*MGlmT2**3 + 
     &   (MGl2*(5*MGl2 - MTy2) - MTy2**2)*MStop2(1) - 
     &   Cd(3999)*MStop2(2) + Cd(3403)*MStop2(2)**2 - 
     &   5*MStop2(2)**3

        Cd(4326) = 10*MTy2*(MGl2*(MGl2 - 3*MTy2) + 3*MTy2**2) - 
     &   MTy2*(17*MGl2 + 7*MTy2)*MStop2(1) - Cd(4325)*MStop2(2) + 
     &   Opt(13)

        Cd(4333) = 6*MTy2*(MGl2*(MGl2 - 3*MTy2) + 3*MTy2**2) - 
     &   3*MTy2*(3*MGl2 + MTy2)*MStop2(1) - Cd(4332)*MStop2(2) + 
     &   Opt(13)

        Cd(4972) = 18*MTy2**5 + 11*MTy2**4*MStop2(1) - 
     &   4*MTy2**3*MStop2(1)**2 + 56*MTy2**2*MStop2(1)**3 - 
     &   24*MTy2*MStop2(1)**4 - MStop2(1)**5 + 
     &   MStop2(2)*(-Cd(4969) + 
     &      2*MStop2(2)*(Cd(4970) - Cd(4971)*MStop2(2) + 
     &         Cd(1592)*MStop2(2)**2))

        Cd(6365) = 3*MGl2*(MGl2**2 - 2*MTy2**2) - 
     &   8*(3*MTy2**2 + MGl2*(MGl2 + 3*MTy2))*MStop2(1) + 
     &   (7*MGl2 + 18*MTy2)*MStop2(1)**2 - 2*MStop2(1)**3 + 
     &   MStop2(2)*(-2*Cd(6363) + Cd(6364)*MStop2(2))

        Cd(6736) = -4*MTy2**3 + 
     &   MGl2*(MGl2*(MGl2 - Cd(98)) + Cd(1944)) + 
     &   MTy2**2*MStop2(1) + (-10*MTy2 + Cd(5058))*MStop2(1)**2 + 
     &   MStop2(2)*(Cd(6735) - Cd(3127)*MStop2(2)) + 
     &   Cd(286)*MStop2(3)**2

        Cd(6835) = MStop2(1)**5 - 4*Cd(6834)*MStop2(2)**3 + 
     &   MTy2*(-6*MStop2(1)**4 + 
     &      MTy2*(6*MStop2(1)**3 + 
     &         MTy2*(MTy2*(12*MTy2 - 21*MStop2(1)) + 
     &            12*MStop2(1)**2)) + 4*Cd(2049)*MStop2(2)**2) + 
     &   Cd(45)*(-(Cd(6833)*MStop2(2)) + 2*MStop2(2)**4)

        Cd(8053) = MGl2*(MGl2*(-MGl2 + Cd(44)) + Cd(3573)) - 
     &   2*Cd(47)*MGlpTmSt2(2) - 
     &   MTy2*(-MStop2(1)**2 + MTy2*(MTy2 + MStop2(1))) + 
     &   MStop2(2)*(Cd(8052) + MStop2(2)*(-Cr(86) + MStop2(2)))

        Cd(10772) = -(Cd(9948)*Cd(10606)) - MGlmT2*Cr(28) + 
     &   2*M_3*AtC*UStop2(1,1)*UStop2(1,2) + 
     &   AtC*UCStop(1,3)**2*(2*M_3C + AtC*UStop2(1,2)) + 
     &   At*(AtC + 2*M_3*UCStopC(1,3)**2 + 
     &      2*M_3C*UStop2(1,1)*UStop2(1,2))

        Cd(10879) = MGl2*Cd(1201) + MTy2*Cd(10448) - Opt(15) - 
     &   MUEC**2*MStop2(2)*UCStopC(1,3)**2 + 
     &   2*MTy*(2*M_3*MUE2 + M_3C*MUEC**2)*UCStopC(1,3)*
     &    UStop2(1,1)*UStop2(1,2) - 
     &   MStop2(2)*(MUE**2*UCStop(1,3)**2 + 
     &      2*MUE2*UStop2(1,1)*UStop2(1,2)**2)

        Cd(11801) = At*Cd(10050) - 
     &   Cd(9948)*Cd(11800)*UStop2(1,1)**2 - 
     &   4*M_3*(MUEC*UCStopC(1,3)**2 + 
     &      MUE*UStop2(1,1)*UStop2(1,2)**2) + 
     &   (-4*M_3C + AtC)*
     &    (MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)**2)

        Cd(1518) = 2*(MGl2 - 2*MTy2)*
     &    (MGl2*(MGl2 - 2*MTy2) + 2*MTy2**2) + 
     &   (MTy2**2 + MGl2*(-MGl2 + MTy2))*MStop2(1) - 
     &   2*MTy2*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(1516) + MStop2(2)*(-Cd(1517) + 5*MStop2(2)))

        Cd(1966) = 2*(-13*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - 3*MTy2) + 9*MTy2**2)) - 
     &   (MTy2**2 + MGl2*(7*MGl2 + 3*MTy2))*MStop2(1) + 
     &   3*(3*MGl2 + 8*MTy2)*MStop2(1)**2 - 5*MStop2(1)**3 - 
     &   Cd(1965)*MStop2(2)

        Cd(2021) = MGl2**3 - 4*MGl2**2*MTy2 + 6*MGl2*MTy2**2 - 
     &   4*MTy2**3 + 2*(MGl2*(MGl2 - MTy2) - MTy2**2)*MStop2(1) + 
     &   4*MTy2*MStop2(1)**2 - 
     &   2*MStop2(2)*(3*Cd(2019) + 
     &      MStop2(2)*(-3*Cd(2020) + 5*MStop2(2)))

        Cd(2059) = MStop2(1)*
     &    (-8*MGl2**2 + 8*MTy2**2 + 13*MGlpT2*MStop2(1) - 
     &      11*MStop2(1)**2) - 2*Cd(2057)*MStop2(2) + 
     &   Cd(2058)*MStop2(2)**2 + Opt(17)

        Cd(2771) = 4*MGlmT2**3 + 
     &   MStop2(1)*(-13*MStop2(1)**2 + 
     &      MTy2*(11*MTy2 + 14*MStop2(1)) + 
     &      MGl2*(-15*MGl2 + 8*MTy2 + 21*MStop2(1))) - 
     &   Cd(2770)*MStop2(2)

        Cd(2813) = 3*MGlmT2**3 + 
     &   MStop2(1)*(-9*MStop2(1)**2 + 
     &      MTy2*(9*MTy2 + 11*MStop2(1)) + 
     &      MGl2*(-11*MGl2 + 5*MTy2 + 15*MStop2(1))) - 
     &   Cd(2770)*MStop2(2)

        Cd(2867) = 10*MGlmT2**3 + 
     &   MStop2(1)*(-19*MStop2(1)**2 + 
     &      MTy2*(7*MTy2 + 23*MStop2(1)) + 
     &      MGl2*(-39*MGl2 + 31*MTy2 + 48*MStop2(1))) + 
     &   Cd(2866)*MStop2(2)

        Cd(6436) = 2*MTy2*(MGl2*(MGl2 - 3*MTy2) + 3*MTy2**2) + 
     &   (-9*MTy2**2 + MGl2*(5*MGl2 + MTy2))*MStop2(1) + 
     &   (-3*MGl2 + 4*MTy2)*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(6434) + MStop2(2)*(-Cd(6435) + MStop2(2)))

        Cd(6525) = 54*MTy2**5 - 27*MTy2**4*MStop2(1) - 
     &   6*MTy2**3*MStop2(1)**2 + 60*MTy2**2*MStop2(1)**3 - 
     &   23*MTy2*MStop2(1)**4 - Cd(6519)*MStop2(2) + 
     &   MStop2(2)**2*(2*Cd(6520) + 
     &      MStop2(2)*(-2*Cd(6521) + Cd(6524) + Cd(6522)*MStop2(2))
     &      ) + Cd(6523)*MStop2(3)**3

        Cd(6657) = -6*(-20*MTy2**3 + 
     &      MGl2*(MGl2**2 + 15*MTy2**2)) + 4*MGl2**2*Cd(6624) + 
     &   2*Cd(6656) + 16*MTy2*(-3*MGl2 + 2*MTy2)*MStop2(1) + 
     &   (-13*MGl2 + 36*MTy2)*MStop2(1)**2 - 
     &   MStop2(2)*(2*Cd(6654) + Cd(6655)*MStop2(2))

        Cd(6762) = 54*MTy2**5 - 134*MTy2**4*MStop2(1) + 
     &   56*MTy2**3*MStop2(1)**2 - 30*MTy2**2*MStop2(1)**3 - 
     &   20*MTy2*MStop2(1)**4 + 
     &   MStop2(2)*(Cd(6759) + 
     &      MStop2(2)*(-2*Cd(6760) + 
     &         MStop2(2)*(Cd(6524) + 2*Cd(6761) - 
     &            Cd(2922)*MStop2(2)))) + Cd(6523)*MStop2(3)**3

        Cd(9558) = 2*(MGl2 - 2*MTy2)*
     &    (MGl2*(MGl2 - 2*MTy2) + 2*MTy2**2) + 
     &   (MTy2**2 + MGl2*(-MGl2 + MTy2))*MStop2(1) + 
     &   2*MTy2*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(9556) + MStop2(2)*(-Cd(9557) + 5*MStop2(2)))

	Cd(10607) = Cd(10772) - Cd(10605)*MStop2(2)

        Cd(11816) = -(Cd(9948)*Cd(11815)) + Opt(16) + 
     &   UCStop(1,3)*(MStop2(3) + 3*M_3C*AtC*UCStop(1,3)**2 + 
     &      UStop2(1,1)*(-(MGlmT2*UStop2(1,2)) + 
     &         MTy2*UStop2(1,2) + 
     &         AtC*(-(M_3*UStop2(1,1)) + M_3*UStop2(1,2))))

        Cd(12445) = -(MGlpT2*Cd(13)) - Opt(19) + 
     &   MStop2(1)*(MUE*UCStop(1,3) + MUEC*UCStopC(1,3))

        Cd(14190) = MGl2*MUEC*TB2*UStop2(1,1) + 
     &   2*M_3*MTy*MUEC*TB2*UCStopC(1,3)*UStop2(1,1) + 
     &   MTy2*(-2*M_3C*TB*UCStop(1,3)**2 + 
     &      MUEC*TB2*UStop2(1,1)) + 
     &   MTy*UCStop(1,3)*(-(TB*
     &         (MGlpT2 - MSbot2(1) + MSq2Diff(1,2))) + 
     &      2*M_3C*MUEC*TB2*UStop2(1,1)) + 
     &   MUEC*TB2*UStop2(1,1)*
     &    (MSq2Diff(1,2) - MSbot2(1)*UStop2(1,1))

        Cd(14192) = -2*M_3*MTy2*TB*UCStopC(1,3)**2 + 
     &   MTy*UCStopC(1,3)*(-(TB*
     &         (MGlpT2 - MSbot2(1) + MSq2Diff(1,2))) + 
     &      2*M_3*MUE*TB2*UStop2(1,1)) + 
     &   UStop2(1,1)*(AtC*TB*Cd(14191) + 
     &      MUE*TB2*(MGlpT2 + MSq2Diff(1,2) + 
     &         2*M_3C*MTy*UCStop(1,3) - MSbot2(1)*UStop2(1,1)))

        Cd(1498) = 2*MGl2**3 - 8*MTy2**3 - 7*MGl2**2*MStop2(1) + 
     &   9*MGl2*MStop2(1)**2 - 5*MStop2(1)**3 + 
     &   MTy2**2*(12*MGl2 + 7*MStop2(1)) - Cd(1497)*MStop2(2) + 
     &   MTy2*(MGl2*(-8*MGl2 + 7*MStop2(1)) + 2*MStop2(2)**2 - 
     &      4*MStop2(3)*MStop2(4))

        Cd(2323) = MStop2(1)*
     &    (-8*MGl2**2 + 8*MTy2**2 + 11*MGlpT2*MStop2(1) - 
     &      5*MStop2(1)**2) - 
     &   MStop2(2)*(2*Cd(2322) + Cd(2058)*MStop2(2)) + Opt(17)

        Cd(2328) = MGl2**3 - 4*MTy2**3 + 3*MTy2**2*MStop2(1) - 
     &   MStop2(1)**3 - MGl2**2*(4*MTy2 + 3*MStop2(1)) + 
     &   3*MGl2*(MStop2(1)**2 + MTy2*(2*MTy2 + MStop2(1))) - 
     &   Cd(2327)*MStop2(2) - MTy2*MStop2(3)*MStop2(4)

        Cd(2535) = -(MStop2(1)*
     &      (2*MGl2**2 + 4*MTy2**2 + MGlpT2*MStop2(1))) - 
     &   Cd(2533)*MStop2(2) + Cd(2534)*MStop2(2)**2 - 
     &   5*MStop2(2)**3 + Opt(17)

        Cd(6388) = MGl2**3 - 8*MGl2**2*MTy2 + 14*MGl2*MTy2**2 - 
     &   8*MTy2**3 + (15*MTy2**2 + MGl2*(-3*MGl2 + MTy2))*
     &    MStop2(1) + 2*(MGl2 - 2*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(-Cd(6386) + MStop2(2)*(-Cd(6387) + MStop2(2)))

        Cd(7291) = 2*(95*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - MTy2) - 23*MTy2**2)) - 
     &   5*Cd(1272)*Cd(1280) + 
     &   (MGl2*(-7*MGl2 - 14*MTy2) + 111*MTy2**2)*MStop2(1) + 
     &   (9*MGl2 + 25*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(-Cd(7289) + Cd(7290)*MStop2(2))

        Cd(7574) = MGl2**3 - 6*MGl2**2*MTy2 + 15*MGl2*MTy2**2 - 
     &   20*MTy2**3 - 2*Cd(7573) + 
     &   4*MGlmT2*(MGl2 - 2*MTy2)*MStop2(1) + 
     &   2*(-7*MGl2 + 4*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(-2*Cd(7570) + Cd(7571)*MStop2(2))

        Cd(7896) = 4*(-20*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - 6*MTy2) + 15*MTy2**2)) - 
     &   2*Cd(7895) - 6*MGlmT2*(MGl2 - 2*MTy2)*MStop2(1) - 
     &   (MGl2 + 12*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(-2*Cd(7892) + Cd(7893)*MStop2(2))

        Cd(8058) = MGl2**3 - 8*MGl2**2*MTy2 + 22*MGl2*MTy2**2 - 
     &   24*MTy2**3 - (MGl2*(3*MGl2 - MTy2) + MTy2**2)*MStop2(1) + 
     &   2*(MGl2 - 2*MTy2)*MStop2(1)**2 + 
     &   MStop2(2)*(-Cd(8057) + MStop2(2)*(-Cd(6387) + MStop2(2)))

        Cd(8107) = 2*(5*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - MTy2) - 3*MTy2**2)) + 
     &   (MGl2*(MGl2 - 3*MTy2) - 5*MTy2**2)*MStop2(1) - 
     &   (MGl2 - 4*MTy2)*MStop2(1)**2 - 
     &   MStop2(2)*(Cd(8105) + MStop2(2)*(-Cd(8106) + MStop2(2)))

        Cd(11724) = -MGlmT2 + Cd(118)*(5*Cd(9948) - Cd(11723)) + 
     &   AtC*(-2*M_3C*UCStop(1,3)**2 - AtC*UCStop(1,3)**2 - 
     &      2*M_3*UStop2(1,1)*UStop2(1,2)) + 
     &   At*(-(Cd(10001)*UCStopC(1,3)**2) - 
     &      2*M_3C*UStop2(1,1)*UStop2(1,2) - 
     &      2*AtC*UStop2(1,1)*UStop2(1,2))

        Cd(14136) = MStop2(1)*UStop2(1,1) - 
     &   MGl2*UStop2(1,1)**2 - 
     &   MTy*UCStopC(1,3)*(At - 2*M_3*UStop2(1,1)) - 
     &   MTy*UCStop(1,3)*(AtC - 2*M_3C*UStop2(1,1)) + 
     &   (-(MUE2*TB2) + MSq2Diff(1,1))*UStop2(1,2)**2 + 
     &   2*At*(M_3*UCStopC(1,3)**2 + 
     &      M_3C*UStop2(1,1)*UStop2(1,2)) + 
     &   AtC*(2*M_3C*UCStop(1,3)**2 + 
     &      2*M_3*UStop2(1,1)*UStop2(1,2) - At*UStop2(1,2)**2)

        Cd(6954) = 4*(-20*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - 6*MTy2) + 15*MTy2**2)) - 
     &   14*MGlmT2*(MGl2 - 2*MTy2)*MStop2(1) + 
     &   (19*MGl2 - 28*MTy2)*MStop2(1)**2 - 2*Cd(6952)*MStop2(2) + 
     &   (-2*Cd(1357) + Cd(6953))*MStop2(2)**2 - 
     &   5*Cd(872)*MStop2(3)**2

        Cd(8300) = 4*(-20*MTy2**3 + 
     &      MGl2*(MGl2*(MGl2 - 6*MTy2) + 15*MTy2**2)) - 
     &   18*MGlmT2*(MGl2 - 2*MTy2)*MStop2(1) + 
     &   (29*MGl2 - 36*MTy2)*MStop2(1)**2 - 2*Cd(8298)*MStop2(2) - 
     &   (2*Cd(4513) + Cd(8299))*MStop2(2)**2 - 
     &   11*Cd(872)*MStop2(3)**2

        Cd(10682) = Cd(9948)*Cd(10681) + MStop2(3) - 
     &   UStop2(1,1)*(4*MTy*
     &       (M_3C*UCStop(1,3) + M_3*UCStopC(1,3))*UStop2(1,1)+
     &        MTy2*(UStop2(1,1) - UStop2(1,2))**2 - 
     &      MStop2(1)*(UStop2(1,1) - UStop2(1,2))**2) - 
     &   MGl2*(UStop2(1,2)**3 + 
     &      UStop2(1,1)*(UStop2(1,2)**2 + 
     &         UStop2(1,1)*(UStop2(1,1) + 2*UStop2(1,2))))

        Cd(10831) = Opt(20) - 
     &   MStop2(2)*(MUE*UCStop(1,3) + MUEC*UCStopC(1,3))

        Cd(11808) = Opt(20) - 
     &   MStop2(1)*(MUE*UCStop(1,3) + MUEC*UCStopC(1,3))

        Cd(4978) = 2*(-3*MTy2**4 + 
     &      MGl2*(2*MTy2**3 + 
     &         MGl2*(MGl2*(MGl2 - 3*MTy2) + 2*MTy2**2))) - 
     &   7*MGl2*(MGl2**2 - 2*MTy2**2)*MStop2(1) + 
     &   2*(MGl2*(3*MGl2 - 7*MTy2) - 7*MTy2**2)*MStop2(1)**2 - 
     &   (MGl2 - 6*MTy2)*MStop2(1)**3 - 
     &   MStop2(2)*(Cd(4975) + 
     &      MStop2(2)*(Cd(4976) + 
     &         MStop2(2)*(-Cd(4977) + 2*MStop2(2))))

        Cd(10065) = At*Cd(10064) - 
     &   2*MUE*UCStop(1,3)**2*
     &    (AtC - 4*M_3C*Cd(9948) - 4*M_3C*UStop2(1,1)) + 
     &   4*M_3*MUE*((Cd(9948) - UStop2(1,1))*UStop2(1,1)**2 - 
     &      2*UStop2(1,1)**2*UStop2(1,2) - 
     &      UStop2(1,1)*UStop2(1,2)**2 - UStop2(1,2)**3) + 
     &   MUEC*(8*M_3*UCStopC(1,3)**2 + 
     &      UStop2(1,1)**2*
     &       (4*M_3C*Cd(9948) - 4*M_3C*UStop2(1,1) + 
     &         AtC*(3*Cd(9948) + UStop2(1,1))) - 
     &      2*(4*M_3C + AtC)*UStop2(1,1)**2*UStop2(1,2) + 
     &      (-4*M_3C*UStop2(1,1) + AtC*(Cd(9948) + UStop2(1,1)))*
     &       UStop2(1,2)**2 - 4*M_3C*UStop2(1,2)**3)

        Cd(10676) = MGl2*Cd(10674) + Cd(9948)*Cd(10675) + 
     &   2*MUE*MStop2(1)*UCStop(1,3)**2*UStop2(1,1) + 
     &   MTy2*(-2*MUE*UCStop(1,3)**2 + 
     &      MUEC*UStop2(1,1)*(UStop2(1,1) - UStop2(1,2))**2) + 
     &   4*MTy*(-3*M_3C*MUE*UCStop(1,3)**3 + 
     &      (M_3*MUE + M_3C*MUEC)*UCStop(1,3)*
     &       (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &         UStop2(1,2)**2) + 
     &      M_3*MUEC*UCStopC(1,3)*
     &       (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &         UStop2(1,2)**2)) + 
     &   MUEC*(-(MStop2(2)*UStop2(1,1)**3) + 
     &      UStop2(1,2)*((-2*MStop2(1) + MStop2(3))*
     &          UStop2(1,1)**2 - 
     &         (MStop2(2) + 2*MStop2(3))*UStop2(1,1)*UStop2(1,2) - 
     &         MStop2(1)*UStop2(1,2)**2))

        Cd(14193) = AtC*Cd(14190) + At*Cd(14192) + 
     &   MGl2*TB*UStop2(1,1)*
     &    (MUE2*TB2 - MSq2Diff(1,2)*UStop2(1,1)) - 
     &   MTy*UCStopC(1,3)*(MGl2*MUEC*TB2 + 
     &      MUEC*TB2*(-MSbot2(1) + MSq2Diff(1,2)) + 
     &      2*M_3*TB*(-(MUE2*TB2) + MSq2Diff(1,2))*UStop2(1,1) + 
     &      MTy2*(MUEC*TB2 - 2*M_3*TB*UStop2(1,2))) - 
     &   MTy*UCStop(1,3)*(MGl2*MUE*TB2 + 
     &      MUE*TB2*(-MSbot2(1) + MSq2Diff(1,2)) + 
     &      2*M_3C*TB*(-(MUE2*TB2) + MSq2Diff(1,2))*
     &       UStop2(1,1) + MTy2*(MUE*TB2 - 2*M_3C*TB*UStop2(1,2))
     &      ) + MTy2*(-2*TB2*
     &       (M_3C*MUE*UCStop(1,3)**2 + 
     &         M_3*MUEC*UCStopC(1,3)**2) + 
     &      MUE2*TB*TB2*UStop2(1,1) + 
     &      TB*(MGl2*UStop2(1,2) - MSbot2(1)*UStop2(1,2)**2 + 
     &         MSq2Diff(1,2)*(-UStop2(1,1) + UStop2(1,2))))

        Cd(10699) = MGl2*Cd(10696) + Cd(9948)*Cd(10698) - 
     &   12*M_3*MTy*MUEC*UCStopC(1,3)**3 + 
     &   2*MUEC*UCStopC(1,3)**2*
     &    (MStop2(1) + MTy2*(-1 - UStop2(1,2)) + 
     &      MStop2(3)*UStop2(1,2)) + 
     &   4*MTy*(M_3*MUE + M_3C*MUEC)*UCStopC(1,3)*
     &    (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &      UStop2(1,2)**2) + 
     &   MUE*(-(MStop2(2)*UStop2(1,1)**3) + 
     &      MTy2*UStop2(1,1)*(UStop2(1,1) - UStop2(1,2))**2 + 
     &      4*M_3C*MTy*UCStop(1,3)*
     &       (UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &         UStop2(1,2)**2) + 
     &      UStop2(1,2)*((-2*MStop2(1) + MStop2(3))*
     &          UStop2(1,1)**2 - 
     &         (MStop2(2) + 2*MStop2(3))*UStop2(1,1)*UStop2(1,2) - 
     &         MStop2(1)*UStop2(1,2)**2))

        Opt(1) = UStop2(1,1)*
     &    (Cd(9948)*Cd(10496) - 4*MTy2*UCStop(1,3)*UStop2(1,2)) - 
     &   At**2*UCStopC(1,3)*(UStop2(1,1)**2 + UStop2(1,2)**2)

        Opt(5) = MGlpT2*MUE*UCStop(1,3)**2 + 
     &   MUEC*UStop2(1,1)*(MGl2*Cd(9948)*UStop2(1,1) + 
     &      UStop2(1,2)*(MTy2 + MGl2*UStop2(1,2)))

        Opt(9) = 6*MTy2**4 - 12*MTy2**3*MStop2(1) + 
     &   62*MTy2**2*MStop2(1)**2 - 24*MTy2*MStop2(1)**3 + 
     &   2*MStop2(2)*(-2*Cd(5292) + 
     &      MStop2(2)*(Cd(5293) - 2*Cd(5294)*MStop2(2)))

        Opt(10) = 2*MTy2*(MTy2 - MStop2(1))**2*
     &    (3*MTy2 + 4*MStop2(1)) - 4*Cd(8204)*MStop2(2) - 
     &   2*MTy2*Cd(1852)*MStop2(2)**2 + 4*Cd(2092)*MStop2(2)**3

        Opt(11) = 2*MTy2*(MTy2 + 3*MStop2(1))*
     &    (MTy2*(MTy2 - 5*MStop2(1)) + 2*MStop2(1)**2) - 
     &   2*MStop2(2)*(2*Cd(7151) + 
     &      MStop2(2)*(Cd(7152) - 2*Cd(2086)*MStop2(2)))

        Opt(12) = 2*MGl2**2*Cd(1952) + 
     &   MTy2*(-8*MGl2*MTy2 + 12*MTy2**2 - 3*MStop2(1)**2)

        Opt(14) = MGlpT2*
     &    (MUE**2*UCStop(1,3)**2 + MUEC**2*UCStopC(1,3)**2) + 
     &   2*MUE2*UStop2(1,1)*
     &    (MGl2*Cd(9948)*UStop2(1,1) + 
     &      UStop2(1,2)*(MTy2 + MGl2*UStop2(1,2)))

	Cd(9) = 2*Cd(7) + Cd(8)*Cr(13)

	Cd(96) = MGl2 + Cd(48)

	Cd(985) = 2*(MStop2(1) - Cd(594)*UStop2(1,2))

        Cd(9935) = 2*SA2*MStop2(2)*UStop(2,1) - 
     &   2*Cd(9808)*Cd(9810)*UStopC(2,2)

        Cd(9936) = 2*Cd(9792)*Cd(9929)*MStop2(2) + 
     &   Cd(9935)*MStop2(3)*UStopC(2,1)

	Cd(9937) = Cd(9934)*UStopC(2,2) + Cd(9807)**2*UUStopC(2,1)

	Cd(9938) = 4*Cd(9808)*UUStopC(3,2) - 4*Cd(9810)*UUStopC(3,3)

	Cd(9942) = Cd(9808)*UUStopC(2,2) - 2*Cd(9810)*UUStopC(2,3)

	Cd(9943) = Cd(9808)*Cd(9942) + Cd(9807)**2*UUStopC(2,1)

	Cd(988) = Cd(985) + 2*MStop2(2)*UStop2(1,1)

	Cd(992) = 4*(3*MStop2(1)*MStop2(3) + Cd(523)*UStop2(1,2))

	Cd(10019) = M_3*Cd(10018) - 3*M_3C*UCStop(1,3)**2

	Cd(10021) = M_3*Cd(10020) - 3*M_3C*UCStop(1,3)**2

	Cd(994) = MTy*Cd(981) - AtC*Cd(984)*UCStop(1,3)

	Cd(10026) = Cd(164) - Cd(10025)*MGlpTmSt2(1)

	Cd(10030) = Cd(251) + Cd(10029)*MGlpTmSt2(2)

	Cd(97) = Cd(96)*UStop(1,1) + 2*M_3*MTy*UStop(1,2)

	Cd(10038) = Cd(251) - Cd(10037)*MGlpTmSt2(2)

	Cd(10042) = M_3C*Cd(10018) - 3*M_3*UCStopC(1,3)**2

	Cd(10043) = M_3C*Cd(10020) - 3*M_3*UCStopC(1,3)**2

	Cd(10057) = MUE*Cd(9966) - 2*MUEC*UCStopC(1,3)**2

	Cd(10066) = 4*Cd(251)*Cd(10056) - Cd(10065)*MGlpTmSt2(2)

	Cd(1001) = 3*MStop2(3)**2 - Cd(975)*UStop2(1,2)

	Cd(10121) = Cd(164) + 2*MGlpTmSt2(1)*UStop2(1,1)

	Cd(10122) = Cd(251) + 2*MGlpTmSt2(2)*UStop2(1,1)

	Cd(10124) = Cr(28)*MGlpTmSt2(2) + Cd(251)*UStop2(1,1)

	Cd(10130) = Cd(164)*Cd(9950) - Cd(47)*Cd(10112)*MGlpTmSt2(1)

	Cd(10148) = 3*M_3 - At*Cd(10025)

	Cd(10149) = 3*M_3C*UCStop(1,3) + Cd(10148)*UCStopC(1,3)

	Cd(10151) = Cd(251) - Cd(10025)*MGlpTmSt2(2)

	Cd(1007) = Cd(947) + MStop2(3)*MStop2(4)

	Cd(1008) = MStop2(1)**2 + Cd(1007)*UStop2(1,2)

	Cd(10166) = Cd(164)*Cd(10163) + Cd(10165)*MGlpTmSt2(1)

	Cd(10167) = Cd(251)*Cd(10163) + Cd(10165)*MGlpTmSt2(2)

	Cd(10171) = 4*Cd(164)*Cd(10169) - Cd(10170)*MGlpTmSt2(1)

	Cd(10173) = At*Cd(9970) + 4*Cd(10172)

	Cd(10174) = 4*Cd(251)*Cd(10169) + Cd(10173)*MGlpTmSt2(2)

	Cd(10178) = 4*M_3C*Cd(9966) + Cd(10177)*UCStopC(1,3)**2

	Cd(10179) = 4*Cd(164)*Cd(10176) + Cd(10178)*MGlpTmSt2(1)

	Cd(10183) = 4*Cd(251)*Cd(10176) + Cd(10182)*MGlpTmSt2(2)

	Cd(10192) = Cd(164) + 2*MGlpTmSt2(1)

        Cd(10194) = (-(Cd(152)*Cd(2597)*Cd(10193)) + 
     &     4*Cd(10192)*MStop2(1)**2)/MStop2(1)**2

	Cd(10195) = Cd(47)**2 - 4*MGlpTmSt4(1)*MStop2(1)**2

        Cd(10196) = At*Cd(10195) + 
     &   16*M_3*MTy2*MGlpTmSt2(1)*MStop2(1)**2

	Cd(10202) = M_3C*Cd(10200) + Cd(10201)*UCStopC(1,3)

	Cd(10205) = Cd(164) + 4*MGlpTmSt2(1)*UStop2(1,1)

	Cd(10212) = 4*Cd(10211) - AtC*Cd(9981)*UStop2(1,1)

	Cd(1013) = MStop2(1)**2 + Cd(286)*MStop2(3)

        Cd(10217) = Cd(47)**2*Cd(9981) - 
     &   Cd(152)*Cd(2597)*MStop2(1)*UStop2(1,2)

	Cd(10218) = Cd(47)**2 - Cd(152)*Cd(2597)*MStop2(1)

	Cd(10223) = MGl2*Cd(9948)*UStop2(1,1) + Cd(10222)*UStop2(1,2)

        Cd(10226) = 2*M_3C*MTy2*Cd(795)*MGlpTmSt2(1) - 
     &   MTy*MUE*MGlpTmSt4(1)*UCStop(1,3)

        Cd(10231) = Cd(152)*Cd(2597)*Cd(9981)*MStop2(1) - 
     &   Cd(47)**2*UStop2(1,2)

	Cd(1015) = Cd(1008) - 2*MStop2(1)**2*UStop2(1,1)

	Cd(10238) = -4*MTy2 + 4*MTy*Cd(585)

        Cd(10240) = Cd(138)*Cd(251) - 
     &   2*Cd(10239)*MGlpTmSt2(2)*UStop2(1,1)

	Cd(10247) = MTy2 - MTy*Cd(585)

	Cd(10253) = -2*Cd(10252) + AtC*Cd(9981)*UStop2(1,1)

	Cd(1017) = 10*MStop2(1)*MStop2(3) + Cd(1012)*UStop2(1,1)

	Cd(10261) = Cd(179)*Cd(2681) - 4*MGlpTmSt2(2)*MStop2(2)**2

	Cd(10265) = MGl2*Cd(9948)*UStop2(1,1) + Cd(10264)*UStop2(1,2)

	Cd(1019) = 4*MStop2(1)**2 + Cd(286)*MStop2(3)

	Cd(10279) = 2*Cd(152)*Cd(2597) + 4*MTy2*Cd(164)*MStop2(1)

	Cd(10280) = Cd(152)*Cd(2597) + 2*MTy2*Cd(164)*MStop2(1)

	Cd(10283) = Cd(47) - MGlmT2*Cd(9948)*Cd(10015)

	Cd(10284) = 2*MGlmT2*Cd(47) + Cd(10280)

	Cd(1020) = Cd(951) + MStop2(3)*MStop2(4)

	Cd(10285) = Cd(47)**2 + At*AtC*Cd(9948)*Cd(10284)

	Cd(10286) = At*AtC*Cd(47)**2 + Cd(10285)*MStop2(3)

	Cd(10287) = Cd(47)**2 + Cd(10284)*MStop2(3)*UStop2(1,1)

	Cd(10288) = At*AtC*Cd(10287) - Cd(47)**2*MStop2(3)

	Cd(10289) = 4*M_3C*MGlmT2*Cd(47) + AtC*Cd(10284)

	Cd(10290) = 2*AtC*Cd(47)**2 - Cd(10289)*MStop2(3)

	Cd(10291) = 4*M_3*MGlmT2*Cd(47) + At*Cd(10284)

	Cd(10292) = 2*At*Cd(47)**2 - Cd(10291)*MStop2(3)

	Cd(1021) = MStop2(1)**2 + Cd(1020)*UStop2(1,2)

	Cd(10301) = MGlmT2*Cd(10297) + Cd(10300)*MStop2(1)

	Cd(10304) = 2*Cd(47)*Cd(9955)*Cd(10303) + Cd(10278)*MStop2(3)

	Cd(10305) = Cd(10304)*MStop2(1) + Cd(10280)*MStop2(3)

	Cd(10307) = Cd(10278) + Cd(47)*Cd(10306)

	Cd(10308) = 3*M_3 + At*Cd(10306)

	Cd(10309) = 3*M_3C*UCStop(1,3) + Cd(10308)*UCStopC(1,3)

        Cd(10314) = (Cd(10278)*Cd(10311)*MStop2(3) - 
     &     2*MGlmT2*Cd(47)*Cd(9955)*Cd(10313)*UStop2(1,2))/
     &   MStop2(3)

	Cd(1023) = Cd(1021) - Cd(1019)*UStop2(1,1)

        Cd(10317) = -2*Cd(47)**2*UCStop(1,3) + 
     &   2*At*Cd(10280)*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(10318) = -2*Cd(47)**2 + At*Cd(10039)*Cd(10280)

	Cd(10326) = -2*MGlmT2*Cd(47)*Cd(10025) + Cd(10279)

	Cd(10328) = 2*Cd(47)*UStop2(1,1) + M_3C*At*MGlmT2*UStop2(1,2)

        Cd(10331) = At*MGlmT2*Cr(5)*Cr(75) - 
     &   Cd(47)*UCStop(1,3)*UStop2(1,2)

        Cd(10332) = Cd(47)*UCStop(1,3)*UStop2(1,1) + 
     &   At*MGlmT2*Cr(5)*UStop2(1,2)

	Cd(10333) = Cd(9948)*Cd(10332) + Cd(10331)*UStop2(1,1)

	Cd(100) = Cd(48)**2 + MGl2*Cd(99)

	Cd(10337) = Cd(10278) - Cd(164)*Cd(9969)*MStop2(1)

	Cd(10338) = Cd(10280) + 2*Cd(10337)*MStop2(1)

	Cd(10339) = 13*Cd(47) + 4*Cd(10278)

	Cd(10345) = 8*M_3C + AtC*Cd(9972)

	Cd(10346) = 4*M_3*Cd(10344) + Cd(10345)*UCStop(1,3)**2

        Cd(10351) = -8*Cd(10211)*Cd(10280) + 
     &   AtC*Cd(9981)*Cd(10279)*UStop2(1,1)

	Cd(10354) = Cd(10279) + 4*Cd(10278)*Cd(10353)

	Cd(1027) = Cd(832)*UStop2(1,1)**2 + Cd(1013)*UStop2(1,2)**2

        Cd(10363) = -(MGlmT2*Cd(47)*Cd(10362)) + 
     &   4*Cd(10015)*Cd(10278)*MStop2(3)*UStop2(1,1)

	Cd(10368) = 2*(AtC*UCStop(1,3)**2 + At*Cd(10367)*UStop2(1,1))

	Cd(10370) = At*Cd(9969)*Cr(28) + 2*AtC*Cd(10077)*UCStop(1,3)**2

	Cd(10377) = AtC*Cd(9969)*Cr(28) + 2*At*Cd(10077)*UCStopC(1,3)**2

        Cd(10385) = MUEC*Cd(10383)*Cr(28) + 
     &   2*MUE*Cd(10384)*UCStop(1,3)**2

        Cd(10386) = MGlmT2*Cd(47)*Cd(10058)*MStop2(1) - 
     &   MUEC*Cd(10280)*MStop2(3)*UStop2(1,1)**2

        Cd(10389) = MUE*Cd(10383)*Cr(28) + 
     &   2*MUEC*Cd(10384)*UCStopC(1,3)**2

        Cd(10390) = MGlmT2*Cd(47)*Cd(10059)*MStop2(1) - 
     &   MUE*Cd(10280)*MStop2(3)*UStop2(1,1)**2

	Cd(10395) = MGlmT2*Cd(10055) - Cd(10057)*MStop2(1)

        Cd(10396) = 2*Cd(10280)*Cd(10388) + 
     &   MUE*Cd(10279)*UStop2(1,1)**2*UStop2(1,2)

	Cd(10397) = 2*Cd(47)*Cd(10395)*MStop2(1) - Cd(10396)*MStop2(3)

	Cd(10400) = MUEC*Cd(9966) - 2*MUE*UCStop(1,3)**2

        Cd(10402) = MUEC*Cd(10401)*Cr(28) - 
     &   2*MUE*Cd(9969)*UCStop(1,3)**2

        Cd(10403) = MUE*Cd(10401)*Cr(28) - 
     &   2*MUEC*Cd(9969)*UCStopC(1,3)**2

        Cd(10408) = Cd(44)*Cd(10279) + 
     &   2*Cd(10278)*MGlpTmSt2(1)*MStop2(1)

        Cd(10409) = M_3C*Cd(795)*UCStop(1,3) + 
     &   M_3*Cd(13)*UStop2(1,1)*UStop2(1,2)

        Cd(10410) = M_3*Cd(10219)*UCStopC(1,3) + 
     &   M_3C*Cd(13)*UStop2(1,1)*UStop2(1,2)

	Cd(10422) = 2*Cd(10421) + Cd(10413)*MStop2(3)

	Cd(10423) = 2*AtC*Cd(10018)*Cr(5) + Cd(10422)*UCStopC(1,3)

	Cd(10425) = MTy2*UStop2(1,1) - Cd(118)*UStop2(1,2)

	Cd(10428) = At*Cd(10096) - Cd(9955)*Cd(10425)

	Cd(10429) = -2*Cd(10428) + Cd(10025)*MStop2(3)

	Cd(10430) = 2*AtC*Cd(10091)*Cr(5) - Cd(10429)*UCStopC(1,3)

	Cd(101) = -2*MGl2 + Cd(98)

	Cd(10436) = -2*MUE2*Cd(485) + Cd(427)*Cd(1272)

	Cd(10439) = Cd(152)*MStop2(2) + MGlmT2*MStop2(1)*UStop2(1,2)

        Cd(10444) = -(Cd(47)*Cd(10443)*Cr(28)*MStop2(1)) + 
     &   Cd(10280)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(1036) = 2*Cd(1035)*MStop2(3) - 8*Cd(1013)*UStop2(1,2)**2

	Cd(10445) = MGlmT2 - Cd(9972)*MStop2(1)

	Cd(10446) = 4*Cd(47)*Cd(10445)*MStop2(1) + Cd(10279)*MStop2(3)

	Cd(10449) = -(Cd(427)*Cr(85)) + MUE2*Cr(28)*UStop2(1,1)

        Cd(10451) = -4*Cd(47)*Cd(152)*Cd(10450)*MStop2(1) + 
     &   2*MUE2*Cd(10279)*MStop2(3)*UStop2(1,1)**2

	Cd(1037) = Cd(1034) + 2*Cd(1015)*UStop2(1,1)

	Cd(1038) = 10*MStop2(1)*MStop2(3) - Cd(1036)*UStop2(1,1)

	Cd(10467) = M_3C*UCStop(1,3)**2 + Cd(10466)*UStop2(1,1)

	Cd(10477) = (Cd(10476) + 4*Cd(9969)*MStop2(1))/MStop2(1)

	Cd(10478) = 5*MGlmT2*Cd(47) - 2*Cd(152)*Cd(164)*Cd(10015)

	Cd(10481) = -4*Cd(10480) + AtC*Cd(9969)*Cr(28)

	Cd(10490) = M_3C*Cd(10488) + Cd(10489)*UCStopC(1,3)

	Cd(10491) = Cd(152)*Cd(10487) + MTy2*Cd(10490)

	Cd(1041) = 5*Cd(968) + 24*MTy2*UStop2(1,2)

	Cd(10494) = 2*Cd(9955)*Cd(10493) + Cd(10025)*MStop2(3)

	Cd(10503) = At*Cd(10501) - M_3*Cd(10502)

	Cd(10509) = 2*At*Cd(10091)*Cr(5) + Cd(10508)*UCStop(1,3)

	Cd(10513) = 2*Cd(10512) + Cd(10025)*MStop2(3)*UCStop(1,3)

	Cd(10514) = MGlmT2*Cd(10513) + Cd(10509)*MStop2(1)

	Cd(10519) = Cd(47) + MStop2(1)**2

	Cd(10521) = Cd(47) + Cd(10520)*MStop2(1)

	Cd(102) = Cd(48)*Cd(89) - MGl2*Cd(101)

	Cd(10534) = -(MGlmT2*Cd(9972)) + MStop2(1)

	Cd(10535) = At*Cd(10533)*Cr(28) + 2*AtC*Cd(10534)*UCStop(1,3)**2

        Cd(10536) = AtC*Cd(10533)*Cr(28) + 
     &   2*At*Cd(10534)*UCStopC(1,3)**2

	Cd(10537) = MUE*Cd(10535) + MUEC*Cd(10536)

	Cd(10540) = -2*Cd(152)*Cd(10052) + AtC*Cd(10539)

        Cd(10542) = Cd(47)*Cd(152)*Cd(10055) - 
     &   2*Cd(164)*Cd(10541)*MStop2(3)

	Cd(1046) = Cd(1041) + 2*MTy2*Cd(1045)*UStop2(1,1)

        Cd(10546) = -2*Cd(47)*Cd(9955)*Cd(10453) + 
     &   Cd(152)*Cd(164)*MStop2(3)

        Cd(10547) = Cd(10306)*MStop2(3) + 
     &   2*Cd(9955)*MStop2(1)*UStop2(1,2)

        Cd(10549) = Cd(10548)*MStop2(3) - 
     &   2*Cd(9955)*MStop2(1)*UStop2(1,2)

	Cd(10550) = MGlmT2*Cd(10549) + Cd(10547)*MStop2(1)

	Cd(1047) = Cd(1046)*UStop2(1,1) + Cd(1044)*UStop2(1,2)

        Cd(10557) = Cd(10556)*UCStopC(1,3) + 
     &   AtC*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(10561) = At*Cd(10042) + Cd(9955)*MStop2(1)*UStop2(1,2)

	Cd(10562) = 2*Cd(10561) + Cd(10306)*MStop2(3)

	Cd(1048) = 4*MStop2(2)**2 - Cd(595)*MStop2(3)

	Cd(10564) = At*Cd(10096) - Cd(9955)*MStop2(1)*UStop2(1,2)

	Cd(10565) = 2*Cd(10564) - Cd(10563)*MStop2(3)

	Cd(10566) = MGlmT2*Cd(10565) + Cd(10562)*MStop2(1)

	Cd(10572) = Cd(47) - Cd(152)*MStop2(1)

	Cd(1049) = Cd(973) - MStop2(2)*MStop2(3)

	Cd(10573) = MGlmT2*Cd(1281)*Cd(9948) - Cd(10572)*UStop2(1,1)

	Cd(10574) = Cd(47) - Cd(152)*Cd(1281)

	Cd(10575) = MGlmT2*Cd(9948)*MStop2(1) - Cd(10574)*UStop2(1,1)

	Cd(10581) = MGl2*Cd(1201) - Cd(10448)*MStop2(1)

        Cd(10582) = Cd(47)*Cd(152)*Cd(10450) + 
     &   2*MUE2*Cd(164)*MStop2(1)*MStop2(3)*UStop2(1,1)**2

	Cd(1050) = 16*MStop2(1)**2 - 4*Cd(1049)*UStop2(1,2)

	Cd(10587) = M_3C*Cd(10585) + M_3*Cd(10586)

	Cd(10591) = M_3C*Cd(10588) + M_3*Cd(10590)

	Cd(10595) = M_3C*Cd(10592) + M_3*Cd(10594)

	Cd(10596) = -(MGlmT2*Cd(10595)) + Cd(10591)*MStop2(1)

	Cd(10603) = 2*Cd(179)*Cd(2681) + 4*MTy2*Cd(251)*MStop2(2)

	Cd(10604) = Cd(179)*Cd(2681) + 2*MTy2*Cd(251)*MStop2(2)

        Cd(103) = 2*(Cd(100)*MStop2(3)*UStop(1,1) + 
     &     M_3*MTy*Cd(102)*UStop(1,2))

        Cd(1056) = -2*Cd(1055)*MStop2(1) + 
     &   2*Cd(1054)*MStop2(2)*UStop2(1,2)

        Cd(10645) = -8*Cd(10211)*Cd(10604) + 
     &   AtC*Cd(9981)*Cd(10603)*UStop2(1,1)

	Cd(10651) = 4*M_3*Cd(10649) + Cd(10650)*UCStop(1,3)**2

	Cd(10652) = M_3C*Cd(10649) + 2*M_3*Cd(9969)*UCStopC(1,3)**2

	Cd(1057) = Cd(1050) - 4*Cd(1048)*UStop2(1,1)

	Cd(10653) = AtC*Cd(9970) + 4*Cd(10652)

	Cd(10657) = 4*Cd(10353)*Cd(10602) + Cd(10603)

	Cd(10659) = M_3C*Cd(10658) - 2*M_3*UCStopC(1,3)**2

	Cd(10660) = AtC*Cd(10658) - 4*Cd(10659)

	Cd(1058) = Cd(1053) - Cd(1057)*UStop2(1,1)

	Cd(10667) = MUEC*Cd(10666) - MUE*Cd(843)*UCStop(1,3)**2

        Cd(10672) = MUE*Cd(10670)*UCStop(1,3) - 
     &   MUEC*Cd(10671)*UStop2(1,1)*UStop2(1,2)

	Cd(1059) = Cd(1056)*MStop2(3) - Cd(1058)*UStop2(1,1)

	Cd(10686) = MUE*Cd(10666) - MUEC*Cd(843)*UCStopC(1,3)**2

        Cd(10690) = -(MUEC*Cd(10689)*UCStopC(1,3)) + 
     &   MUE*Cd(10671)*UStop2(1,1)*UStop2(1,2)

	Cd(10693) = 2*M_3C*MTy*Cr(75) - Cd(44)*UCStopC(1,3)

	Cd(10706) = MGlmT2*Cd(1952) - Cd(10705)*MStop2(2)

	Cd(10713) = MUEC*Cd(10712) - 2*MUE*Cd(491)*UCStop(1,3)**2

	Cd(10714) = 2*MTy2*Cd(10713) - MUE*MStop2(3)*UCStop(1,3)**2

	Cd(10716) = Cd(10715)*MStop2(3) - 2*MGlpT2*Cr(75)*UStop2(1,2)

	Cd(10717) = MUE*Cd(10712) - 2*MUEC*UCStopC(1,3)**2

	Cd(10718) = 2*MTy2*Cd(10717) - MUEC*MStop2(3)*UCStopC(1,3)**2

	Cd(10720) = 2*MTy2*Cd(10674) + Cd(795)*MStop2(3)

	Cd(10723) = Cd(10722)*MStop2(3) - 2*MTy2*Cr(75)*UStop2(1,2)

        Cd(104) = Cd(64)*Cd(96)*Cd(97)*MStop2(3) - 
     &   Cd(103)*T134(MGl2,MTy2,MStop2(2))

	Cd(1065) = 2*Cd(1064)*MStop2(1) + 4*Cd(1020)*UStop2(1,2)

	Cd(10734) = 2*Cd(251)*MStop2(1) - MGlmT2*MStop2(3)

	Cd(1066) = Cd(1065) - 2*Cd(1061)*UStop2(1,1)

	Cd(10746) = MGlmT2*Cd(10745) - Cd(9972)*MStop2(1)*MStop2(2)

	Cd(10760) = MGlmT2*Cd(10021) - Cd(10759)*MStop2(2)

	Cd(10762) = MGlmT2*Cd(10020) - MStop2(2)*UStop2(1,1)**2

	Cd(10765) = MGlmT2*Cd(10043) - Cd(10764)*MStop2(2)

	Cd(1069) = 3*MStop2(1) - 2*Cd(286)*UStop2(1,2)

	Cd(10773) = -8*Cd(10772) + 8*Cd(10605)*MStop2(2)

        Cd(10781) = -(Cd(179)*UCStop(1,3)) + 
     &   2*At*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(1070) = 22*MStop2(3) - 4*Cd(1069)*UStop2(1,2)

	Cd(1071) = 3*Cd(985) - 2*MStop2(3)*UStop2(1,2)

	Cd(1072) = 8*MStop2(3) - Cd(874)*UStop2(1,2)**2

	Cd(10803) = Cd(10658)*MStop2(1) + Cd(9958)*MStop2(3)*UStop2(1,1)

        Cd(10812) = -(Cd(10811)*UCStop(1,3)**2) + 
     &   M_3*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(10828) = MGlmT2 - Cd(9972)*MStop2(2)

	Cd(1075) = Cd(1071) + 2*Cd(1074)*UStop2(1,1)

	Cd(10832) = -(Cd(10715)*MStop2(3)) + 2*MGlpT2*Cr(75)*UStop2(1,2)

        Cd(10836) = M_3C*Cd(10835)*UCStop(1,3) + 
     &   M_3*Cd(13)*UStop2(1,1)**2

	Cd(10837) = MUEC*Cd(10020) - 3*MUE*UCStop(1,3)**2

        Cd(10838) = M_3*Cd(13)*Cd(10020) + 
     &   M_3C*Cd(10837)*UCStop(1,3)

	Cd(1076) = Cd(1070) - Cd(1075)*UStop2(1,1)

        Cd(10842) = M_3*Cd(10841)*UCStopC(1,3) + 
     &   M_3C*Cd(13)*UStop2(1,1)**2

	Cd(10843) = MUE*Cd(10020) - 3*MUEC*UCStopC(1,3)**2

        Cd(10844) = M_3C*Cd(13)*Cd(10020) + 
     &   M_3*Cd(10843)*UCStopC(1,3)

	Cd(10850) = 2*Cd(179)*Cd(10086) + AtC*MStop2(2)*UCStop(1,3)**2

	Cd(10854) = AtC + 4*M_3C*Cd(491)

	Cd(10855) = 2*M_3*Cd(10712) - Cd(10854)*UCStop(1,3)**2

	Cd(10856) = MGlmT2*Cd(10853) + Cd(10855)*MStop2(2)

	Cd(10859) = AtC*Cd(10049) + At*Cd(10050)

	Cd(10863) = 2*Cd(179)*Cd(10039) + At*MStop2(2)*UCStopC(1,3)**2

	Cd(10866) = 2*M_3C*Cd(10712) - Cd(10865)*UCStopC(1,3)**2

	Cd(10867) = MGlmT2*Cd(10864) + Cd(10866)*MStop2(2)

	Cd(10882) = 2*MUE2*Cd(10020) - 3*MUEC**2*UCStopC(1,3)**2

        Cd(10883) = MUE**2*Cd(10020)*UCStop(1,3) + 
     &   Cd(10882)*UCStopC(1,3)

        Cd(10903) = AtC*MUE*UCStop(1,3)**2 + 
     &   Cd(424)*UStop2(1,1)*UStop2(1,2)

        Cd(10905) = At*MUEC*Cd(10904)*UCStopC(1,3)**2 + 
     &   2*Cd(10903)*UStop2(1,2)

	Cd(10914) = S2A*Cd(10905) + 2*SA2*Cd(10902)*UStop2(1,2)

	Cd(10915) = CA2*Cd(10913) + 2*Cd(10914)*UStop2(1,1)

	Cd(10916) = MStop2(1) - Cd(286)*UStop2(1,1)

	Cd(10917) = 6*MStop2(3) + 2*Cd(10916)*UStop2(1,1)

	Cd(1084) = -2*Cd(972) + Cd(1047)

	Cd(10) = A0delStop(1)*dMSfsq1133eps(-1) - Cd(9)*MStop2(1)

	Cd(1085) = Cd(981) + Cd(1059)

        Cd(10936) = Cd(10934) + 
     &   3*Cd(10921)*MStop2(1)*MStop2(3)*UCStop(3,1)*UStop2(1,1)

        Cd(10937) = 2*Cd(10935)*MStop2(3)*UStop(1,2) - 
     &   2*Cd(10936)*UStopC(1,1)

        Cd(10966) = Cd(10965)*UUStop(1,2) - 
     &   2*Cr(91)**2*MStop2(1)**2*UUStopC(2,1)

	Cd(10969) = 2*(Cd(10966) - 2*Cd(10968)*MStop2(3))

	Cd(10978) = MUEC*Cd(10977) + MUE*Cd(9794)*UUStop(1,1)

	Cd(10986) = MUE*Cd(9809)*UStop(1,1) + MUEC*Cd(10985)*UStop(2,2)

	Cd(10988) = 2*MUEC*Cd(10987)*UStopC(2,1) + Cd(10986)*UStopC(2,2)

	Cd(10990) = AtC*Cd(565)*UStop(1,1) + Cd(10989)*UStop(1,2)

	Cd(10992) = MUEC*Cd(10991) + MUE*Cd(9794)*UStop(1,1)

	Cd(10996) = MTy*Cd(10995) + 2*At*UUStop(1,2)*UUStop(2,3)

        Cd(10997) = Cd(10994)*UCStop(3,1)*UCStopC(1,3) + 
     &   Cd(10996)*UUStopC(2,1)

	Cd(10998) = 2*Cd(10997)*MStop2(1)**2 + MTy*Cd(10993)*UUStop(1,2)

	Cd(10999) = -(AtC*Cd(10938)) + 2*Cd(10998)

	Cd(11003) = Cd(698)*UCStop(3,1) - 2*Cd(11000)*UCStop(3,3)

	Cd(1093) = -Cd(948) + Cd(1008)

        Cd(11014) = 2*Cd(11013)*MStop2(1) - 
     &   MUE*Cd(1143)*MStop2(3)*UCStop(1,3)*UStop(1,1)

	Cd(11015) = MUE*Cd(9794)*UStop(2,1) + 2*MTy*MUEC*UUStop(3,2)

	Cd(11016) = Cd(11015)*UStopC(2,1) + 2*MTy*MUE*Cr(91)*UStopC(2,2)

	Cd(11021) = Cd(11016) + Cd(11020)*UStop2(1,2)**2

        Cd(11022) = Cd(11021)*MStop2(1) - 
     &   MTy*MUE*Cd(1143)*MStop2(3)*UStop(1,1)*UStop2(1,2)

	Cd(11039) = 2*MTy*Cd(10995) - At*Cr(91)*UUStop(3,2)

	Cd(11040) = MUE*Cd(11039) - 4*MTy*MUEC*UUStop(1,2)*UUStop(2,2)

	Cd(11041) = Cd(11038)*UUStopC(2,1) - Cd(11040)*UUStopC(2,3)

	Cd(11043) = 5*Cd(11016) + 3*Cd(11020)*UStop2(1,2)**2

        Cd(11044) = Cd(11043)*MStop2(1) + 
     &   8*MTy*MUE*MStop2(3)*UStop(1,1)*UStop2(1,2)**3

        Cd(11046) = Cd(11041)*MStop2(1)**2 - 
     &   MTy*MUE*Cd(11035)*UStop(1,2)

	Cd(11047) = 3*Cd(11016) + Cd(11020)*UStop2(1,2)**2

	Cd(11049) = AtC*Cd(11014) + 2*Cd(11022)

	Cd(11052) = Cd(10984) - Cd(10978)*UStop2(1,1)*UStopC(1,2)

        Cd(11061) = Cd(11059)*MStop2(1) + 
     &   MTy2*Cd(1143)*MStop2(3)*UStop2(1,2)

        Cd(11062) = Cd(11061)*UStop(1,2) + 
     &   2*MTy*Cd(11060)*MStop2(1)*UStop2(1,2)**2*UUStop(3,2)

	Cd(11064) = Cd(9794)*Cd(9948) + 2*MTy*Cr(91)*UStopC(2,2)

	Cd(11065) = Cd(11064) + Cd(11063)*UStop2(1,2)**2*UStopC(2,2)

        Cd(11066) = Cd(11065)*MStop2(1) - 
     &   MTy*Cd(1143)*MStop2(3)*UStop(1,1)*UStop2(1,2)

	Cd(11070) = 8*MTy*Cd(11069)*UCStop(3,1) - Cd(11068)*UStopC(2,1)

        Cd(11071) = Cd(11070)*MStop2(1)**2*UCStopC(1,3) - 
     &   2*MTy2*Cd(11067)*UUStop(1,2)

	Cd(11074) = -2*Cd(11073)*MStop2(3) + Cd(11071)*UStopC(1,2)

        Cd(11076) = Cd(11070)*MStop2(1)**2*UCStopC(1,3) - 
     &   2*MTy2*Cd(11075)*UUStop(1,2)

	Cd(11079) = Cd(11078)*MStop2(3) - Cd(11076)*UStopC(1,2)

	Cd(11081) = 4*MStop2(1) + Cd(1143)*MStop2(3)

        Cd(11082) = 4*Cd(9948)*Cd(11080)*MStop2(1)*UStop(1,1) - 
     &   Cd(11081)*UUStop(1,1)

        Cd(11092) = Cd(11037)*UStop2(1,2)*UStopC(2,1) - 
     &   Cd(11039)*UUStopC(2,3)

	Cd(11094) = 5*Cd(11064) + 3*Cd(11063)*UStop2(1,2)**2*UStopC(2,2)

        Cd(11095) = Cd(11094)*MStop2(1) + 
     &   8*MTy*MStop2(3)*UStop(1,1)*UStop2(1,2)**3

	Cd(11097) = Cd(11092)*MStop2(1)**2 - MTy*Cd(11035)*UStop(1,2)

	Cd(11098) = 3*Cd(11064) + Cd(11063)*UStop2(1,2)**2*UStopC(2,2)

	Cd(1102) = -Cd(947) + 2*Cd(1019)

	Cd(11101) = Cd(11062) - AtC*Cd(11066)

        Cd(11115) = 6*At*MTy*UStop(1,2)*UStop2(1,1) - 
     &   Cd(11114)*UStopC(2,1)

	Cd(11116) = Cd(11113) + Cd(11115)*UStop2(1,2)**2

	Cd(11117) = 12*Cd(9794) - 5*At*UStop(1,2)*UStop2(1,2)**2

        Cd(11118) = Cd(11116)*MStop2(1) + 
     &   MTy*Cd(11117)*MStop2(3)*UStop2(1,2)

	Cd(1104) = 2*Cd(832) - Cd(842)

        Cd(11120) = 6*Cd(9794)*MStop2(3)**2 - 
     &   4*Cd(11119)*MStop2(1)**2*UStop2(1,1)*UStopC(2,2)

	Cd(11124) = Cd(11120) + 5*Cd(9794)*MStop2(3)**2

	Cd(11125) = 2*Cd(11121)*UUStopC(2,1) + Cd(11123)*UUStopC(2,3)

	Cd(11126) = Cd(11125)*MStop2(1)**2 + MTy*Cd(11124)*UStop(1,2)

	Cd(11127) = 5*Cd(11113) + 3*Cd(11115)*UStop2(1,2)**2

	Cd(11129) = Cd(11125)*MStop2(1)**2 + MTy*Cd(11120)*UStop(1,2)

	Cd(11130) = 3*Cd(11113) + Cd(11115)*UStop2(1,2)**2

	Cd(11133) = 2*MTy*Cd(11001) + At*Cr(91)*UUStop(3,2)

	Cd(11134) = Cd(698)*UCStop(3,1) - 2*Cd(11132)*UCStop(3,3)

	Cd(11135) = Cd(11134)*UStop2(1,2) + Cd(11133)*UUStopC(2,2)

        Cd(11136) = 2*Cd(11135)*MStop2(1)**2 - 
     &   11*Cd(9794)*MStop2(3)**2*UStop(1,1)

	Cd(1106) = 22*MStop2(3) - 4*Cd(535)*UStop2(1,2)**2

        Cd(11140) = 5*Cd(9809)*UCStop(3,3) + 
     &   3*Cd(11139)*UStop2(1,2)*UUStopC(3,2)

	Cd(11142) = 12*Cd(9794) + 5*Cd(10980)*UStop2(1,2)**2

        Cd(11143) = Cd(9809)*UCStop(3,3) + 
     &   Cd(11139)*UStop2(1,2)*UUStopC(3,2)

        Cd(11144) = -2*Cd(11143)*MStop2(1) + 
     &   Cd(11142)*MStop2(3)*UCStop(1,3)

        Cd(11145) = -2*Cd(11135)*MStop2(1)**2 + 
     &   6*Cd(9794)*MStop2(3)**2*UStop(1,1)

        Cd(11146) = 3*Cd(9809)*UCStop(3,3) + 
     &   Cd(11139)*UStop2(1,2)*UUStopC(3,2)

	Cd(11149) = 2*Cd(11118) + AtC*Cd(11144)

	Cd(11154) = 2*Cd(11054) + Cd(9794)*UStop2(1,1)

	Cd(11155) = Cd(11058) - Cd(9794)*Cd(11154)*UStop2(1,1)

        Cd(11159) = -(MUE*Cd(9794)*UCStop(1,3)*UStop(1,1)) + 
     &   MUEC*Cd(11158)*UStop2(1,2)

	Cd(11164) = Cd(11163)*MStop2(1) + Cd(11159)*MStop2(3)

	Cd(1109) = 4*Cd(1003) - Cd(1106)

        Cd(11172) = Cd(11171)*MStop2(1) - 
     &   Cd(10990)*MStop2(3)*UStop2(1,2)

	Cd(11173) = 3*Cd(10988)*MStop2(1) + Cd(11172)*UUStopC(1,2)

	Cd(11183) = Cd(698)*Cd(9948) - Cd(11182)*UStopC(2,2)

        Cd(11185) = -(Cd(11184)*MStop2(3)*UCStop(1,3)) + 
     &   Cd(11183)*MStop2(1)**2*UCStop(3,3)

	Cd(11188) = MTy*Cd(11001) - 2*At*UUStop(1,3)*UUStop(2,2)

        Cd(11191) = 2*Cd(9948)*Cd(11186)*MStop2(1) - 
     &   Cd(10980)*MStop2(3)*UStop(1,2)*UStop2(1,2)**2

        Cd(11207) = -2*Cd(11206)*MStop2(1)**2 + 
     &   6*Cd(9794)**2*MStop2(3)**2

	Cd(1113) = -Cd(982) + Cd(1063)

	Cd(11223) = Cd(11221)*UUStopC(3,1) + MUE*Cd(11222)*UUStopC(3,3)

	Cd(11225) = -(Cd(11224)*MStop2(1)**2) + 3*Cr(26)*MStop2(2)

	Cd(11226) = 2*MUE2*Cd(10020) - MUEC**2*UCStopC(1,3)**2

	Cd(11232) = 2*Cd(11163)*MStop2(1) + 4*Cd(11159)*MStop2(3)

        Cd(11234) = 4*Cd(9948)*Cd(11186)*MStop2(1) - 
     &   2*Cd(11233)*MStop2(3)

	Cd(11238) = MUEC*Cd(11234) + 2*Cd(11237)*UCStop(1,3)

        Cd(11239) = Cd(11171)*MStop2(1) - 
     &   2*Cd(10990)*MStop2(3)*UStop2(1,2)

	Cd(11240) = Cd(10988)*MStop2(1) - Cd(11239)*UUStopC(1,2)

	Cd(11244) = Cd(11238) - Cd(11232)*UStop2(1,1)

	Cd(11249) = 5*Cd(9794)**2 + 2*Cd(11056)*UStop2(1,2)**2

	Cd(11262) = Cd(10917) - Cd(10920)

	Cd(11264) = Cd(10937) - Cd(10941)

	Cd(11267) = Cd(10964) - Cd(10969)

        Cd(11268) = 4*MUE2*Cd(11266)*MStop2(1)*UStop2(1,1)*
     &    UStopC(2,2) - MUEC**2*Cd(11267)*UUStopC(1,1)

	Cd(11270) = 2*Cd(11223)*MStop2(1) - Cd(11269)*MStop2(3)

        Cd(11289) = 2*Cd(11287)*MStop2(3) - 
     &   2*Cd(11288)*MStop2(2)*UUStopC(3,1)

        Cd(11291) = Cd(11289) + 
     &   2*MUE**2*Cr(91)*MStop2(3)*UCStop(1,3)*UStop2(1,1)

	Cd(11292) = Cd(11290)*UStopC(1,2) + Cd(11291)*UUStopC(3,2)

        Cd(11294) = Cd(970)*UStop(1,2) - 
     &   UCStop(3,1)*UStop(2,2)*UStop2(1,2)

        Cd(11295) = Cd(970)*MStop2(3)*UCStopC(1,3) - 
     &   2*Cd(11294)*MStop2(2)*UStopC(1,1)

        Cd(11296) = 4*Cd(11295)*MStop2(3)*UCStopC(1,3) - 
     &   4*Cr(91)**2*MStop2(2)**2*UUStopC(1,1)*UUStopC(2,1)

        Cd(11298) = 2*(MUE**2*Cd(970)*UUStop(2,1) + 
     &     Cd(11297)*UUStop(2,2))

        Cd(11302) = Cd(11298)*MStop2(3) + 
     &   4*Cd(11301)*MStop2(2)*UUStopC(3,1)

        Cd(11303) = -2*Cd(11302)*MStop2(3) + 
     &   4*Cd(11300)*MStop2(2)**2*UUStopC(1,1)

	Cd(11304) = 2*(MUE*Cd(970) + MUEC*UCStopC(1,3)**2)

        Cd(11308) = Cd(11305)*MStop2(2)*UCStop(3,1)*UStopC(2,2)*
     &    UUStopC(1,1) + Cd(11307)*MStop2(3)*UUStopC(3,3)

        Cd(11309) = -4*Cd(11308)*MStop2(2) + 
     &   2*Cd(9948)*Cd(11304)*MStop2(3)**2*UStop2(1,1)

	Cd(11310) = MUEC*Cd(11296) + 2*Cd(11309)

	Cd(11311) = MUEC*Cd(11310) - Cd(11303)*UUStopC(2,2)

	Cd(11335) = 2*Cd(698)*MStop2(2) - 2*Cd(11334)*MStop2(3)

	Cd(11336) = MTy*Cd(10995) + At*UUStop(1,2)*UUStop(2,3)

	Cd(11337) = At*Cd(11336) - 2*MTy2*Cr(58)*UUStop(3,1)

        Cd(11338) = AtC*Cd(11335)*UCStopC(3,4) - 
     &   2*At*Cd(11333)*MStop2(3)*UStopC(1,1)

        Cd(11339) = Cd(11338)*UStop2(1,1) - 
     &   4*Cd(11337)*MStop2(2)*UUStopC(1,3)

        Cd(11345) = Cd(11344)*UUStop(2,3) + 
     &   Cd(652)*MStop2(2)*UCStopC(3,1)*UStop2(1,2)*UUStop(3,3)

        Cd(11346) = 8*MTy2*UCStopC(3,1)*UStop2(1,2)*
     &    UUStop(3,1) + Cd(698)**2*UUStopC(1,1)

	Cd(1127) = 2*Cd(1123)*UStop2(1,1)**2 + Cd(1126)*UStop2(1,2)

        Cd(11352) = 2*Cd(11351)*MStop2(3) - 
     &   2*Cd(11133)*MStop2(2)*UUStopC(3,2)

        Cd(11353) = Cd(11352)*UStopC(1,1) + 
     &   2*Cd(11182)*MStop2(2)*UStop2(1,1)*UUStopC(3,3)

	Cd(11354) = 4*MTy*Cd(11350) + Cd(11353)*MStop2(2)

	Cd(11360) = 2*Cd(11205) + AtC*Cd(698)*UUStop(3,1)

	Cd(11362) = 3*Cd(9794)*Cd(9809) + Cd(11114)*UStop2(1,1)**2

	Cd(1129) = 4*Cd(652) - Cd(1122)*Cr(75)**2

	Cd(11370) = AtC*Cd(11326) + 4*Cd(11354)

	Cd(11371) = 4*Cd(11322) + AtC*Cd(11370)

	Cd(11372) = 2*Cd(11329) + Cd(9809)*MStop2(3)*UStop2(1,2)

	Cd(11378) = -2*Cd(11188) - AtC*Cr(91)*UUStop(3,1)

	Cd(11379) = MUEC*Cd(11378) - MUE*Cd(698)*UUStop(3,1)

	Cd(11396) = MUEC*Cd(11387) - Cd(11390)*UUStopC(2,2)

	Cd(11399) = MUEC*Cd(11398)*UStop(2,2) - MUE*Cd(9809)*UUStop(3,1)

	Cd(11402) = 2*Cd(10996) - AtC*Cr(91)*UUStop(3,1)

	Cd(11403) = MUEC*Cd(11402) - MUE*Cd(698)*UUStop(3,1)

        Cd(11404) = Cd(11403)*MStop2(2)*UStopC(1,1) + 
     &   2*Cd(11401)*MStop2(3)*UStopC(2,2)

	Cd(11407) = MUE*Cd(11133) + 4*MTy*MUEC*UUStop(1,2)*UUStop(2,2)

	Cd(11411) = Cd(11366)*MStop2(3) - Cd(11410)*UCStop(1,3)

	Cd(11416) = 2*MUE*Cd(11411) - MUEC*Cd(11415)

	Cd(11417) = Cd(11404) + Cd(11399)*MStop2(3)*UStop2(1,2)

        Cd(11428) = Cd(11427)*MStop2(3) - 
     &   2*Cd(11284)*MStop2(2)*UUStopC(3,1)

	Cd(11430) = Cd(11429) + MUE2*Cd(11286)*UStop2(1,2)

	Cd(11433) = Cd(11285) - Cd(11428)

	Cd(11437) = 3*Cd(970) + Cd(1136)

        Cd(11438) = Cd(2364)*UStop(1,2) - 
     &   UCStop(3,1)*UStop(2,2)*UStop2(1,2)

	Cd(11440) = -2*Cd(970) + Cd(11437)

        Cd(11446) = Cd(11445)*MStop2(3) + 
     &   2*Cd(11301)*MStop2(2)*UUStopC(3,1)

	Cd(11448) = 3*Cd(11304) + 2*Cd(11447)

	Cd(11451) = 2*Cd(11304) - Cd(11448)

	Cd(11452) = 2*Cd(11307) - Cd(11450)

        Cd(11453) = Cd(9948)*Cd(11451)*MStop2(3)*UStop2(1,1) - 
     &   2*Cd(11452)*MStop2(2)*UUStopC(3,3)

	Cd(1138) = MStop2(1) + MTy2*Cd(1137)*UStop2(1,2)

	Cd(11462) = 4*MTy2*Cd(11461)*MStop2(2) + Cd(11460)*MStop2(3)**2

        Cd(11465) = Cd(11464)*MStop2(3) - 
     &   4*Cd(11119)*MStop2(2)*UStop2(1,1)*UStopC(1,2)

        Cd(11469) = Cd(11468)*MStop2(3)**2 - 
     &   4*Cd(1281)*MStop2(2)*UStop2(1,1)**2

        Cd(11472) = Cd(11471)*MStop2(3) + 
     &   4*MStop2(2)*UStop(1,1)*UStop2(1,1)**2

	Cd(11475) = 6*Cd(11329) - 2*Cd(11474)*MStop2(3)

	Cd(1140) = Cd(1139)*UStop2(1,1) + Cd(1138)*UStop2(1,2)

        Cd(11494) = Cd(11493)*MStop2(3) - 
     &   4*Cr(91)*MStop2(2)*UStop2(1,1)*UStop2(1,2)

        Cd(11495) = Cd(11491)*MStop2(3)**2*UStop(2,1)*
     &    UStop2(1,1) + Cd(11494)*MStop2(2)*UStopC(1,2)

        Cd(11498) = Cd(11497)*MStop2(3) - 
     &   2*Cd(11133)*MStop2(2)*UUStopC(3,2)

        Cd(11499) = Cd(11498)*UStopC(1,1) + 
     &   2*Cd(11182)*MStop2(2)*UStop2(1,1)*UUStopC(3,3)

	Cd(11500) = 2*MTy*Cd(11495) + Cd(11499)*MStop2(2)

	Cd(11501) = Cd(640)*Cd(1657) + MStop2(3)*MStop2(4)*UStop2(1,1)

        Cd(11508) = Cd(11504)*MStop2(3)*UCStopC(3,1) + 
     &   4*Cd(11069)*MStop2(2)*UStop2(1,1)*UStop2(1,2)*UStopC(1,2)

	Cd(11509) = 2*MTy*Cd(11508)*MStop2(2) - Cd(11507)*MStop2(3)**2

	Cd(11511) = Cd(9794)*Cd(9809) + Cd(11114)*UStop2(1,1)**2

        Cd(11513) = 2*Cd(11510)*UStop(1,1) - 
     &   AtC*Cd(814)*UCStop(1,3)*UUStop(3,1)

	Cd(11527) = Cd(1428) + 3*Cd(11375)

	Cd(112) = Cd(110) + MGlpTmSt2(1)*UStop2(1,1)

	Cd(1145) = Cd(1144)*MStop2(1) - Cd(1143)*MStop2(2)*UStop2(1,2)

        Cd(11534) = Cd(11403)*MStop2(2)*UStopC(1,1) + 
     &   3*Cd(11401)*MStop2(3)*UStopC(2,2)

        Cd(11535) = MUE*Cd(11471)*MStop2(3)*UStop(2,1) - 
     &   Cd(11405)*MStop2(2)*UStop2(1,1)*UStopC(1,1)

	Cd(11539) = Cd(11519)*MStop2(3) - 2*Cd(11538)*UCStop(1,3)

        Cd(11542) = Cd(11541)*MStop2(3) + 
     &   12*MStop2(2)*UStop2(1,1)**2*UUStop(3,4)

	Cd(11545) = AtC*Cd(11542) + 2*Cd(11544)*UUStop(3,2)

	Cd(11546) = MUE*Cd(11539) + MUEC*Cd(11545)

	Cd(11547) = 2*Cd(11534) + 3*Cd(11399)*MStop2(3)*UStop2(1,2)

	Cd(1147) = Cd(1145) - Cd(1146)*UStop2(1,1)

	Cd(11549) = 2*Cd(11375) - Cd(11527)

        Cd(11577) = 4*MTy2*Cd(9948)*Cd(11574)*MStop2(2) + 
     &   Cd(11576)*MStop2(3)

        Cd(11578) = 2*MTy*Cd(11463)*MStop2(2)*UCStopC(1,3) + 
     &   Cd(11577)*UStop(2,2)

        Cd(11581) = Cd(11579)*MStop2(3)*UStop(2,1)*UStop2(1,1) + 
     &   2*Cd(11580)*MStop2(2)*UStopC(1,2)

	Cd(11582) = MTy*Cd(11581) + Cd(11496)*MStop2(2)*UStopC(1,1)

	Cd(11583) = 2*Cd(11329) - 3*Cd(11474)*MStop2(3)

	Cd(11591) = 2*Cd(640) - Cd(1657)*UStop2(1,1)

	Cd(11594) = MStop2(2)*UStop(2,1) + Cd(11593)*UStop2(1,1)

        Cd(11595) = Cd(11594)*MStop2(3) - 
     &   2*MTy*Cd(11503)*MStop2(2)*UCStopC(3,1)

	Cd(11597) = Cd(11511) - AtC*Cd(11596)*UCStop(1,3)

        Cd(11599) = 2*Cd(11518)*MStop2(2) - 
     &   Cd(11598)*MStop2(3)*UStop(1,1)

	Cd(11607) = 3*Cd(1428) + 2*Cd(11375)

	Cd(11608) = MUEC*Cd(11528) - Cd(11390)*UUStopC(2,2)

	Cd(11612) = Cd(11518) + 2*AtC*Cd(10925)*UCStop(1,3)**2

	Cd(11613) = -(MUEC*Cd(11611)) + MUE*Cd(11612)

        Cd(11625) = (2*Cr(75)**2*MStop2(1) + 
     &     Cd(11624)*MStop2(2))/(MStop2(1)*MStop2(2))

        Cd(113) = M_3C*Cd(112)*UCStop(1,3) + 
     &   M_3*Cd(110)*UCStopC(1,3)

        Cd(11638) = Cd(11626)*MStop2(2) + 
     &   2*Cr(75)**2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(11639) = Cd(11638)*MStop2(1) - 
     &   2*Cd(11637)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(11641) = -(Cd(11640)*MStop2(2)) + 2*Cr(75)**2*MStop2(3)

        Cd(11643) = Cd(11642)*MStop2(2) + 
     &   4*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(11647) = Cr(75)**2*MStop2(1) - Cd(11646)*MStop2(2)

        Cd(11648) = Cd(11642)*MStop2(2) + 
     &   2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(11657) = 4*(Cd(11656)*MStop2(1)*MStop2(2) + 
     &     Cd(5519)*UUStop(3,1)*UUStop(3,2))

        Cd(11658) = 4*MUE*MStop2(1)*MStop2(2)*UUStop(2,1) - 
     &   MUEC*Cd(5519)*UUStop(2,2)

	Cd(11660) = 4*Cd(11658)*UUStop(1,3) - 4*Cd(11659)*UUStop(2,3)

	Cd(11661) = MUE*Cd(11657)*UUStopC(2,2) - Cd(11660)*UUStopC(2,3)

	Cd(11662) = MUE*Cd(11661) + Cd(11655)*UUStopC(2,1)

        Cd(11663) = Cd(13)*Cr(9)*MStop2(1) - 
     &   MUE*MStop2(2)*UCStop(1,3)*UStop2(1,2)

        Cd(11664) = 4*Cd(11216)*Cd(11663)*MStop2(3)*
     &    UStopC(1,2) - Cd(11662)*UUStopC(1,2)

	Cd(11665) = Cd(11221)*UStopC(2,1) + MUE*Cd(11222)*UStopC(2,2)

	Cd(1159) = 2*Cd(1158)*MStop2(3) - Cd(874)*UStop2(1,2)**2

        Cd(11667) = MUE*Cd(5519)*UUStop(1,1) - 
     &   4*MUEC*MStop2(1)*MStop2(2)*UUStop(1,2)

	Cd(11668) = 4*(Cd(11666)*UUStop(1,3) + Cd(11667)*UUStop(2,3))

        Cd(11670) = 4*MUE*MStop2(1)*MStop2(2)*UUStop(1,1) - 
     &   MUEC*Cd(5519)*UUStop(1,2)

	Cd(11671) = 4*Cd(11669)*UUStop(1,3) - 4*Cd(11670)*UUStop(2,3)

        Cd(11673) = 4*Cd(5519)*Cd(11672) - 
     &   8*MUE2*Cd(5945)*UUStop(3,1)*UUStop(3,2)

        Cd(11675) = Cd(13)*Cd(11652)*Cr(9)*MStop2(1) - 
     &   Cd(11665)*MStop2(2)*UUStopC(1,2)

        Cd(11679) = MUE*Cd(5519)*UUStop(2,1) - 
     &   4*MUEC*MStop2(1)*MStop2(2)*UUStop(2,2)

	Cd(11681) = 4*(Cd(11679)*UUStop(1,3) + Cd(11680)*UUStop(2,3))

	Cd(11682) = MUEC*Cd(11657)*UUStopC(2,1) + Cd(11681)*UUStopC(2,3)

	Cd(11683) = MUEC*Cd(11682) + Cd(11678)*UUStopC(2,2)

	Cd(11684) = Cd(11221)*UUStopC(3,2) - MUEC*Cd(11222)*UUStopC(3,4)

        Cd(11685) = (Cd(11683) + 
     &     4*Cd(11684)*MStop2(2)*MStop2(3))/
     &   (MStop2(2)*MStop2(3)**2)

        Cd(11695) = Cd(11694)*MStop2(2) + 
     &   8*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(11701) = Cd(9948)*MStop2(2) - Cd(674)*UStop2(1,2)

	Cd(11702) = 2*Cd(11700)*MStop2(3) - Cd(11701)*UStop2(1,1)

        Cd(11705) = Cd(11216)*Cd(11663)*UStopC(1,2) + 
     &   2*Cd(10908)*MStop2(3)*UUStopC(1,2)

        Cd(11706) = Cd(11675)*UStopC(1,1) - 
     &   2*Cd(10912)*MStop2(3)*UUStopC(1,3)

        Cd(11707) = (-2*Cd(11684)*MStop2(2) + 
     &     4*Cd(10910)*MStop2(3))/(MStop2(2)*MStop2(3))

	Cd(11708) = Cd(11705) + Cd(11706)

        Cd(11709) = (2*Cd(11708) + 
     &     Cd(11707)*MStop2(2)*MStop2(3)*UUStopC(1,1))/
     &   (MStop2(2)*MStop2(3))

        Cd(11710) = (Cd(11709)*MStop2(3) + 
     &     2*MUEC*Cd(11652)*UCStopC(1,3)*UStop2(1,1)*UStopC(1,1))/
     &   (MStop2(1)*MStop2(3))

	Cd(11712) = Cd(47) - Cd(118)**2*Cd(9969)

	Cd(11713) = Cd(47) - MGl2*Cd(152)

	Cd(11714) = Cd(47) + MGl2*Cd(152)

	Cd(1164) = 6*Cd(512)*MStop2(1) - Cd(1163)*UStop2(1,2)

	Cd(11716) = AtC*Cd(11713)*Cd(11714) + Cd(11715)*MStop2(1)

	Cd(11718) = At*Cd(11713)*Cd(11714) + Cd(11717)*MStop2(1)

	Cd(11719) = 16*(Cd(11716)*UCStop(1,3) + Cd(11718)*UCStopC(1,3))

	Cd(11720) = AtC*MGlmT2 + 4*M_3C*Cd(118)

	Cd(11721) = AtC*MGlmT2 + 8*M_3C*Cd(118)

	Cd(1165) = 5*Cd(607)*UStop2(1,1)**2 + Cd(1141)*UStop2(1,2)**2

	Cd(11729) = Cd(47)*Cd(11726) - MGlmT2*Cd(11728)*MStop2(3)

	Cd(11731) = 4*Cd(152) + At*AtC*UStop2(1,1)

	Cd(11743) = -16*M_3C + AtC*Cd(9972)

	Cd(11744) = 4*M_3*Cd(11742) + Cd(11743)*UCStop(1,3)**2

	Cd(11756) = (-Cd(11755) + 4*Cd(11752)*MStop2(1))/MStop2(1)

	Cd(11765) = 2*Cd(9955)*Cd(10425) + Cd(10025)*MStop2(3)

        Cd(11770) = -(At*Cd(10091)*Cr(5)) + 
     &   Cd(9955)*Cd(10425)*UCStop(1,3)

	Cd(11771) = 2*Cd(11770) + Cd(10025)*MStop2(3)*UCStop(1,3)

	Cd(11774) = MGlmT2*Cd(11771) - 2*Cd(11773)

        Cd(11780) = (-4*Cd(47)*Cd(10370) + 
     &     8*MGlmT2*Cd(9964)*MStop2(3))/(Cd(47)*MStop2(3))

	Cd(11784) = 2*Cd(152)*Cd(10224) + Cd(424)*MStop2(1)*UStop2(1,1)

        Cd(11785) = Cd(47)*Cd(164)*Cd(10055) + 
     &   2*Cd(152)*Cd(10388)*MStop2(1)*MStop2(3)

        Cd(11788) = (-4*Cd(47)*Cd(10377) + 
     &     8*MGlmT2*Cd(9996)*MStop2(3))/(Cd(47)*MStop2(3))

	Cd(11794) = AtC*Cd(11792) - 4*M_3C*Cd(11793)

	Cd(11802) = At*Cd(11792) - 4*M_3*Cd(11793)

	Cd(11810) = MGl2*MGlmT2*UStop2(1,1) + Cd(11809)*UStop2(1,2)

	Cd(11818) = MGlmT2*Cd(11817) + Cd(47)*UCStop(1,3)

        Cd(11819) = MGl2*MGlmT2*UCStop(1,3)*UStop2(1,1) + 
     &   Cd(11818)*UStop2(1,2)

	Cd(11824) = Cd(164)*Cd(10522) + Cd(47)*Cr(28)

	Cd(11825) = Cd(47) + Cd(164)*MStop2(2)

	Cd(11845) = 8*M_3*Cd(138) + At*Cd(9972)*MStop2(1)

	Cd(11847) = 4*Cd(179) + At*AtC*UStop2(1,1)

	Cd(11860) = M_3*AtC*UStop2(1,1) + 2*Cd(138)*UStop2(1,2)

	Cd(11868) = 3*M_3 + At*Cd(10029)

	Cd(11871) = 2*Cd(10249) - Cd(179)*UCStop(1,3)

	Cd(11875) = (-Cd(11874) + 4*Cd(11752)*MStop2(2))/MStop2(2)

	Cd(11880) = 4*M_3*Cd(11878) + Cd(11879)*UCStop(1,3)**2

	Cd(11882) = M_3C*Cd(11878) - 2*M_3*Cd(11881)*UCStopC(1,3)**2

	Cd(11883) = AtC*Cd(9970) + 4*Cd(11882)

	Cd(11889) = AtC*Cd(10070) - 4*Cd(11888)

	Cd(11910) = AtC - 4*M_3C*Cd(11909)

	Cd(11911) = 2*M_3*Cd(11908) + Cd(11910)*UCStop(1,3)**2

	Cd(11916) = -8*M_3 + At*Cd(9972)

	Cd(11917) = 2*M_3C*Cd(11908) + Cd(11916)*UCStopC(1,3)**2

	Cd(11918) = 2*M_3C*Cd(10658) + Cd(10140)*UCStopC(1,3)**2

	Cd(11923) = AtC*MUEC*Cd(9983) + At*Cd(10697)

	Cd(11933) = At*Cd(10043) - Cd(118)*Cr(75)*UStop2(1,2)

	Cd(11934) = 2*Cd(11933) + Cd(10548)*MStop2(3)

	Cd(11935) = 2*AtC*Cd(10020)*Cr(5) + Cd(11934)*UCStopC(1,3)

	Cd(11948) = MTy2*Cd(10749) + 2*MTy*Cd(10584)*Cr(75)

	Cd(11949) = -(MTy2*Cd(10754)) + 2*MTy*Cd(10584)*Cr(75)

	Cd(11953) = M_3*Cd(10880) + M_3C*Cd(11952)

	Cd(11955) = MTy2*Cd(10756) + 2*MTy*Cd(11953)

	Cd(11972) = Cd(11971)*UStop2(1,1) + Cd(290)*UStop2(1,2)

	Cd(1191) = Cd(1080) + Cd(1164)

	Cd(11996) = Cd(11995)*UStop2(1,1) + Cd(2361)*UStop2(1,2)

	Cd(12027) = AtC + M_3C*Cd(12026)

	Cd(12028) = M_3*Cd(12025) + Cd(12027)*UCStop(1,3)**2

	Cd(12031) = AtC + 4*M_3C*Cd(12030)

	Cd(12032) = 2*M_3*Cd(12029)*Cr(28) - Cd(12031)*UCStop(1,3)**2

        Cd(12034) = M_3C*Cd(12029)*Cr(28) - 
     &   2*M_3*Cd(12030)*UCStopC(1,3)**2

	Cd(12037) = M_3C*Cd(12036) - 2*M_3*Cd(10025)*UCStopC(1,3)**2

        Cd(12044) = MGlpTmSt2(1)*MStop2(3) - 
     &   2*Cd(47)*Cd(9955)*UStop2(1,2)

	Cd(12057) = At*Cd(12055) - M_3*Cd(12056)

        Cd(12059) = -(Cd(152)*UCStop(1,3)) + 
     &   2*At*Cr(5)*UStop2(1,1)*UStop2(1,2)

	Cd(12064) = -(M_3*Cd(12056)) + At*Cd(12063)

        Cd(12074) = 2*Cd(12072)*Cd(12073) + 
     &   Cd(9948)*Cd(12070)*Cd(12071)*MStop2(3)

        Cd(12082) = 2*Cd(47)*UStop2(1,1) + 
     &   MStop2(1)*MStop2(3)*UStop2(1,2)

        Cd(12083) = MUE*Cd(12082)*UStop(1,1) - 
     &   3*MUEC*Cd(47)*UCStopC(1,3)*UStop(1,2)

        Cd(12084) = (-Cd(12083) + 
     &     MGlpT2*MUE*MStop2(3)*UStop(1,1)*UStop2(1,2))/MStop2(3)

        Cd(12088) = MGlpTmSt2(1)*MStop2(3)*UStop2(1,1) - 
     &   2*Cd(47)*UStop2(1,2)

        Cd(12090) = MUEC*MGlpTmSt2(1)*UCStopC(1,3)**2 + 
     &   2*Cd(12084)*UStopC(1,1)

        Cd(12095) = (Cd(47)*Cd(12094) + 
     &     MUEC*MGlpTmSt2(1)*MStop2(3)*UStop(1,2)*UStop2(1,1))/
     &   MStop2(3)

	Cd(1203) = Cd(1202) - Cd(491)*UStop2(1,1)

        Cd(12141) = MUE*Cd(9794)*UCStop(1,3)*UStop(1,1) + 
     &   MUEC*Cd(12140)*UStop2(1,2)

	Cd(1207) = 2*Cd(426)*Cd(505) + Cd(424)*Cd(1203)*UStop2(1,1)

        Cd(12144) = Cd(12072)*Cd(12143)*UCStop(3,1) + 
     &   MTy*Cd(12142)*MStop2(3)*UUStop(1,3)

        Cd(12151) = 2*Cd(10925)*Cd(12072)*UCStop(3,1) + 
     &   Cd(12142)*MStop2(3)*UUStop(1,3)

	Cd(12153) = Cd(12150) + AtC*Cd(12152)

	Cd(12154) = -2*MTy*MUE*Cd(12151) + Cd(12072)*Cd(12153)

	Cd(12155) = Cd(12142)*MStop2(3) - 4*Cd(12072)*UStop2(1,1)

	Cd(12159) = 3*MStop2(3) + 2*Cd(1592)*UStop2(1,1)

	Cd(12160) = 2*Cd(45)*MStop2(3) - Cd(44)*Cd(1592)*UStop2(1,1)

	Cd(12161) = MStop2(3) - 2*Cd(44)*UStop2(1,1)

        Cd(12163) = Cd(44)**2*Cd(12161) + 
     &   6*Cd(12162)*MGlpTmSt2(1)*UStop2(1,2)

	Cd(12169) = AtC*UCStop(1,3) - Cd(9809)*UStopC(2,1)

	Cd(12170) = 2*Cd(12168) - MTy*Cd(12169)*MGlpTmSt2(2)

	Cd(12172) = Cd(1667)*UStop(1,1) + At*MTy*UStop(1,2)

	Cd(12173) = 2*MTy2*Cd(44) + MStop2(1)*MStop2(3)

	Cd(12174) = Cd(12173)*UStop(1,1) + At*MTy*Cd(44)*UStop(1,2)

	Cd(12176) = MGl2*Cd(45) + Cd(47)*UStop2(1,1)

	Cd(12195) = Cd(1647)*UStop(1,1) - At*MTy*UStop(1,2)

	Cd(12196) = 2*MTy2*Cd(44) - MStop2(1)*MStop2(3)

	Cd(12197) = Cd(12196)*UStop(1,1) - At*MTy*Cd(44)*UStop(1,2)

	Cd(12198) = MTy2*Cd(88) + 2*MStop2(1)*MStop2(3)

	Cd(12204) = 2*MTy*Cd(1561)*UStop(1,1) - At*Cd(1701)*UStop(1,2)

	Cd(12212) = 3*AtC*UCStop(1,3)*UStop(1,1) - Cd(698)*UStopC(2,1)

        Cd(12230) = 2*Cd(12229) - 
     &   AtC*MTy*MGlpTmSt2(1)*MStop2(3)*UCStop(1,3)

        Cd(12232) = 2*MTy*Cd(12231)*UUStop(1,1) + 
     &   At*MGlpTmSt2(1)*MStop2(3)*UStop2(1,2)*UUStop(1,3)

	Cd(12236) = MTy2*Cd(12235) - 2*MStop2(1)*MStop2(3)

        Cd(12247) = MUEC*Cd(9794)*UCStopC(1,3)**2*UStopC(1,1) + 
     &   MUE*Cd(10985)*UStop2(1,1)*UStop2(1,2)*UStopC(1,2)

        Cd(12248) = (Cd(12247)*MGlpTmSt2(1)*MStop2(3) - 
     &     2*MUE*Cd(47)*Cd(12241)*UStop2(1,2)*UUStopC(3,4))/
     &   MStop2(3)

	Cd(12254) = MUEC*Cd(12253) + MUE*Cd(698)*UStop(2,1)

	Cd(12255) = -3*MGl2 + Cd(198)

	Cd(12256) = Cd(44)*Cd(48) + MGl2*Cd(12255)

        Cd(12257) = Cd(12256)*Cr(91)*MStop2(3) - 
     &   2*Cd(47)*Cd(11106)*MGlpTmSt2(2)*UStop2(1,1)

        Cd(12258) = MTy*Cd(12257) + 
     &   AtC*Cd(47)*Cd(10925)*MGlpTmSt2(2)*UCStop(1,3)

        Cd(12263) = 2*MTy*Cd(11394)*Cd(12256)*MStop2(3) - 
     &   Cd(47)*Cd(12262)*MGlpTmSt2(2)*UStopC(2,2)

	Cd(1220) = 2*MUE2*Cd(1215) - Cd(427)*Cd(1219)

	Cd(12284) = MStop2(2) + Cd(1281)*UStop2(1,1)

	Cd(12287) = Cd(12286)*UStop2(1,1) - Cd(12284)*UStop2(1,2)

        Cd(12288) = -2*MUEC*Cr(85)*MStop2(2)*UCStopC(1,3)**2 + 
     &   MUE*Cd(12287)*UStop2(1,1)

        Cd(12289) = MUEC*Cd(9794)*UCStopC(1,3)*UStop2(1,1)*
     &    UStopC(1,1) + MUE*Cd(10985)*UCStop(1,3)*UStop2(1,2)*
     &    UStopC(1,2)

	Cd(12298) = MGlmT2**3 - Cd(4092)*MStop2(1)

	Cd(12299) = 2*Cd(12298) - Cd(4095)*Cd(9969)*MGlpTmSt2(1)

	Cd(12300) = Cd(47)**2 - Cd(118)*Cd(4095)*Cd(9969)*MGlpTmSt2(1)

	Cd(12301) = Cd(47)**2 - MGl2*Cd(12298)

	Cd(12302) = AtC*Cd(12298) + 2*M_3C*Cd(4095)*MGlpTmSt2(1)

	Cd(12303) = AtC*Cd(12301) + Cd(12302)*MStop2(1)

	Cd(12304) = At*Cd(12298) + 2*M_3*Cd(4095)*MGlpTmSt2(1)

	Cd(12305) = At*Cd(12301) + Cd(12304)*MStop2(1)

	Cd(12306) = Cd(12303)*UCStop(1,3) + Cd(12305)*UCStopC(1,3)

	Cd(12307) = AtC*Cd(164)*Cd(2597) + 4*M_3C*Cd(12298)

	Cd(12308) = AtC*Cd(164)*Cd(2597) + 8*M_3C*Cd(12298)

	Cd(12310) = AtC*Cd(10223) - 4*M_3C*MTy2*UStop2(1,2)

	Cd(12311) = Cd(47)**2 - Cd(164)*Cd(2597)*MGlpTmSt2(1)

	Cd(12317) = Cd(47)**2 - Cd(118)*Cd(12298)

	Cd(12318) = 4*M_3*Cd(12298) + At*Cd(164)*Cd(2597)*UStop2(1,1)

	Cd(12319) = Cd(47)**2 + Cd(164)*Cd(2597)*MStop2(1)*UStop2(1,2)

        Cd(12320) = Cd(9948)*Cd(12318)*MStop2(1) + 
     &   At*Cd(12319)*UStop2(1,2)

	Cd(12321) = AtC*Cd(47)**2 + Cd(12307)*MStop2(1)

	Cd(12322) = Cd(12321)*UCStop(1,3)**2 + Cd(12320)*UStop2(1,1)

        Cd(120) = 2*(Cd(58)*UStop(1,1) + 
     &     M_3*MTy*MGlpTmSt2(1)*UStop(1,2))

	Cd(12323) = Cd(47)**2 + Cd(164)*Cd(2597)*MStop2(1)

	Cd(12325) = At*Cd(164)*Cd(2597) + 4*M_3*Cd(12298)

	Cd(12326) = At*Cd(47)**2 + Cd(12325)*MStop2(1)

	Cd(12331) = MTy2*Cd(9948)*UStop2(1,1) + Cd(12330)*UStop2(1,2)

	Cd(12351) = AtC*Cd(10265) - 4*M_3C*MTy2*UStop2(1,2)

	Cd(12362) = MTy2*Cd(9948)*UStop2(1,1) + Cd(12361)*UStop2(1,2)

	Cd(12369) = 4*MGl2*MTy2*Cd(152) + 2*Cd(4095)*MGlpTmSt2(1)

	Cd(12370) = MGl2*Cd(152)**2 + Cd(12298)

	Cd(12371) = 2*Cd(164)*Cd(2597) + 4*MGl2*Cd(152)*MStop2(1)

	Cd(12372) = Cd(47)*Cd(118) + Cd(12369)

        Cd(1230) = 4*(Cd(523)*UStop2(1,1) + 
     &     4*MStop2(1)*MStop2(2)*UStop2(1,2))

        Cd(12373) = (12*Cd(47)**2 - 
     &     4*Cd(118)*Cd(9969)*Cd(12372))/Cd(47)**2

	Cd(12374) = 2*Cd(47)**2 - MGl2*Cd(12370)

	Cd(12375) = Cd(47)*Cd(118) + Cd(12370)

	Cd(12376) = 3*Cd(47)*Cd(118) + Cd(12369)

	Cd(12377) = AtC*Cd(12375) + M_3C*Cd(12376)

	Cd(12378) = AtC*Cd(12374) + Cd(12377)*MStop2(1)

	Cd(12379) = At*Cd(12375) + M_3*Cd(12376)

	Cd(12380) = At*Cd(12374) + Cd(12379)*MStop2(1)

	Cd(12381) = 8*(Cd(12378)*UCStop(1,3) + Cd(12380)*UCStopC(1,3))

	Cd(12382) = 5*Cd(47)*Cd(118) + Cd(12369)

	Cd(12385) = MGlmT2*Cd(47) + 2*Cd(12370)

	Cd(12391) = 4*MGlmT2*Cd(47) + Cd(12371)

	Cd(12392) = 4*Cd(47)**2 + Cd(12391)*MStop2(3)

	Cd(12395) = MGlmT2*Cd(47) + 4*Cd(12370)

        Cd(12396) = (-(MGlmT2*Cd(47)*Cd(10294)) + 
     &     Cd(10112)*Cd(12395)*MStop2(3))/MStop2(3)

	Cd(12401) = 8*Cd(47)**2*Cd(10208) + Cd(12371)*Cd(12400)

	Cd(12407) = Cd(47)*Cd(12405) + 2*Cd(12406)*MStop2(1)**2

        Cd(12411) = Cd(10198)*Cd(12371) + 
     &   2*Cd(12370)*MStop2(1)*UCStop(1,3)

	Cd(12413) = MGlmT2*Cd(11771) - 2*Cd(12412)*MStop2(1)

        Cd(12424) = 2*Cd(47)*Cd(164)*Cd(10054)*MStop2(1) - 
     &   Cd(795)*Cd(12371)*MStop2(3)

        Cd(12425) = 2*Cd(47)*Cd(164)*Cd(10055)*MStop2(1) - 
     &   Cd(10219)*Cd(12371)*MStop2(3)

	Cd(12426) = AtC*Cd(12424) + At*Cd(12425)

	Cd(12432) = AtC - 4*M_3C*Cd(11881)

	end
