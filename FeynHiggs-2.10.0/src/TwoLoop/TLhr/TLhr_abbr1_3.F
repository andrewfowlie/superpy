* TLhr_abbr1_3.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_3
	implicit none

#include "TLhr.h"


	Cd(5401) = 13*MStop2(2) - 10*MStop2(3)

	Cd(5402) = 13*MStop2(1)**2 + Cd(5401)*MStop2(3)

	Cd(5403) = 3*MStop2(2) + 2*MStop2(3)

	Cd(5404) = 3*MStop2(1)**2 + Cd(5403)*MStop2(3)

	Cd(5405) = 36*MStop2(3) + MStop2(4)

	Cd(5406) = 2*MStop2(2)**2 - Cd(5405)*MStop2(3)

	Cd(5410) = 3*MStop2(2) - 10*MStop2(3)

	Cd(5411) = 3*MStop2(1)**2 + Cd(5410)*MStop2(3)

	Cd(5414) = 4*MStop2(3) - 3*MStop2(4)

	Cd(5415) = 6*MStop2(1)**2 - Cd(5414)*MStop2(3)

	Cd(5416) = 4*MStop2(3) + 67*MStop2(4)

	Cd(5417) = 134*MStop2(2)**2 - Cd(5416)*MStop2(3)

	Cd(5419) = 45*MStop2(2) - 34*MStop2(3)

	Cd(5420) = 45*MStop2(1)**2 + Cd(5419)*MStop2(3)

	Cd(5421) = 44*MStop2(3) + 95*MStop2(4)

	Cd(5422) = 190*MStop2(2)**2 - Cd(5421)*MStop2(3)

	Cd(5424) = 4*MStop2(3) - 9*MStop2(4)

	Cd(5425) = 18*MStop2(1)**2 - Cd(5424)*MStop2(3)

	Cd(535) = 3*MStop2(2) - MStop2(3)

	Cd(5428) = 20*MStop2(2) - 3*MStop2(3)

	Cd(5429) = 20*MStop2(2)**2 - Cd(5428)*MStop2(3)

	Cd(5433) = MTy2 + 10*MStop2(1)

	Cd(536) = 2*(MStop2(2) + MStop2(4)*UStop2(1,2))

	Cd(5438) = MTy2 + 34*MStop2(1)

	Cd(5441) = 10*MStop2(2) - MStop2(3)

	Cd(5442) = 10*MStop2(2)**2 - Cd(5441)*MStop2(3)

	Cd(537) = 2*UStop2(1,1)**2 + UStop2(1,2)**2

	Cd(5448) = 6*MTy2 + 31*MStop2(1)

	Cd(5451) = 12*MTy2*Cd(198) + 5*MStop2(3)**2

	Cd(5452) = 4*MTy2 - 6*MStop2(4)

	Cd(5455) = 7*MStop2(2) - 2*MStop2(3)

	Cd(5456) = 7*MStop2(2)**2 - Cd(5455)*MStop2(3)

	Cd(5457) = MTy2**2 + Cd(5456)

	Cd(538) = UStop2(1,1)*(UStop2(1,1)**3 - 2*UStop2(1,2)**3)

	Cd(5458) = 2*Cd(5457) - 5*MTy2*MStop2(4)

	Cd(5459) = 2*MTy2*Cd(5458) - Cd(1227)*MStop2(4)

	Cd(5460) = 3*MGl2 + Cd(5452)

	Cd(5464) = 60*MStop2(2) - 19*MStop2(3)

	Cd(5477) = 10*MTy2 - 3*MStop2(1)

        Cd(540) = -2*MStop2(1)*
     &   (Cr(28)*MStop2(1) - 2*Cd(535)*UStop2(1,1)*UStop2(1,2))

	Cd(541) = MStop2(2) + MStop2(4)*UStop2(1,2)

	Cd(5489) = 7*MStop2(3) - 12*MStop2(4)

	Cd(5490) = 24*MStop2(2)**2 + Cd(5489)*MStop2(3)

	Cd(5495) = 7*MTy2 - 18*MStop2(1)

	Cd(5501) = 17*MTy2 + 66*MStop2(1)

	Cd(5507) = 14*MTy2 + 39*MStop2(1)

	Cd(5510) = 12*MTy2*Cd(198) - 5*MStop2(3)**2

	Cd(5517) = 6*MStop2(1)**2 - MStop2(3)**2

	Cd(5518) = 2*Cd(5517)*MStop2(2)**2 + MStop2(3)**3*MStop2(4)

	Cd(5519) = 2*MStop2(2)**2 - Cd(2114)*MStop2(3)

	Cd(5525) = 11*MTy2 + 9*MStop2(1)

	Cd(5526) = 9*MStop2(2)**2 - MStop2(3)*MStop2(4)

	Cd(545) = Cr(28)*MStop2(1) - 2*Cd(535)*UStop2(1,1)*UStop2(1,2)

	Cd(5531) = 40*MTy2 + 63*MStop2(1)

	Cd(5532) = 5*MStop2(2) + 4*MStop2(3)

	Cd(546) = Cd(427)*Cd(540) + 2*MUE2*Cd(545)*Cr(28)*MStop2(1)

	Cd(5541) = 23*MTy2 + 16*MStop2(1)

	Cd(547) = 2*Cd(519)*MStop2(2) + 4*MStop2(1)*UStop2(1,2)

	Cd(5549) = 60*MStop2(2) - 17*MStop2(3)

	Cd(5550) = 60*MStop2(2)**2 - Cd(5549)*MStop2(3)

	Cd(5552) = 20*MStop2(2) - 7*MStop2(3)

	Cd(5554) = 12*MStop2(2)**2 - Cd(4662)*MStop2(3)

	Cd(548) = 2*Cr(85)*MStop2(2) + 3*MStop2(1)*UStop2(1,2)

	Cd(5559) = 9*MStop2(2)**2 - Cd(1570)*MStop2(3)

	Cd(5560) = 12*MStop2(2)**2 - Cd(2344)*MStop2(3)

	Cd(53) = 2*M_3*MTy*MGlpTmSt2(1) - MGlpTmSt4(1)*UCStop(1,3)

	Cd(549) = Cd(548) - 2*MStop2(4)*UStop2(1,1)

	Cd(5581) = 11*MTy2 + 7*MStop2(1)

	Cd(5584) = 8*MStop2(3) + 5*MStop2(4)

	Cd(5585) = 10*MStop2(2)**2 - Cd(5584)*MStop2(3)

	Cd(5586) = 21*MStop2(3) + 10*MStop2(4)

	Cd(5587) = 20*MStop2(2)**2 - Cd(5586)*MStop2(3)

	Cd(5596) = 4*MStop2(3) + 7*MStop2(4)

	Cd(552) = 5*MStop2(2) - MStop2(3)

	Cd(5599) = 340*MStop2(2) - 83*MStop2(3)

	Cd(5600) = 340*MStop2(2)**2 - Cd(5599)*MStop2(3)

	Cd(5601) = 7*MStop2(2)**2 - Cd(1060)*MStop2(3)

	Cd(5605) = 25*MTy2 - 8*MStop2(1)

	Cd(553) = 2*MStop2(2) + MStop2(1)*UStop2(1,2)

	Cd(5617) = 20*MTy2 - 7*MStop2(4)

	Cd(5620) = 2*MTy2 + 9*MStop2(4)

	Cd(5621) = 2*MGl2 + Cd(5617)

	Cd(555) = 7*MStop2(2) - 3*MStop2(3)

	Cd(556) = -(Cr(115)*MStop2(1)) + Cd(555)*UStop2(1,1)*UStop2(1,2)

	Cd(5641) = 151*MTy2 + 72*MStop2(1)

	Cd(5643) = 14*MStop2(2)**2 - Cd(5596)*MStop2(3)

	Cd(557) = MStop2(2) + MStop2(3)*UStop2(1,2)

	Cd(5647) = 47*MTy2 - 17*MStop2(1)

	Cd(5648) = 11*MTy2 - MStop2(1)

	Cd(5651) = 303*MTy2 + 397*MStop2(1)

	Cd(5652) = 191*MTy2 - 25*MStop2(1)

	Cd(5657) = 92*MTy2 - 15*MStop2(1)

	Cd(5659) = 216*MStop2(2) + 50*MStop2(3)

	Cd(5660) = 97*MStop2(2) - 15*MStop2(3)

	Cd(5663) = MStop2(2) + 46*MStop2(3)

	Cd(54) = -(M_3C*UCStop(1,3)**2) + M_3*UStop2(1,1)**2

	Cd(558) = 4*MStop2(1)*UStop(1,2)

	Cd(5669) = 4*MStop2(3) - 17*MStop2(4)

	Cd(5670) = 30*MStop2(2)**2 - Cd(5669)*MStop2(3)

	Cd(5671) = 220*MStop2(3) - 7*MStop2(4)

	Cd(5677) = 54*MTy2 + 61*MStop2(1)

	Cd(5687) = 48*MTy2**4 - 5*Cd(5442)*MStop2(3)**2

	Cd(5689) = 12*MStop2(2) + MStop2(3)

	Cd(5690) = 12*MStop2(2)**2 - Cd(5689)*MStop2(3)

	Cd(5691) = 22*MStop2(2) + 3*MStop2(3)

	Cd(5694) = 7*MStop2(3) + 24*MStop2(4)

	Cd(5695) = 48*MStop2(2)**2 - Cd(5694)*MStop2(3)

	Cd(5699) = 16*MTy2**4 - 5*MStop2(1)*MStop2(2)*MStop2(3)**2

	Cd(5704) = 9*MTy2 + 2*Cd(3378)

	Cd(5706) = 7*MStop2(2)**2 - Cd(1226)*MStop2(3)

	Cd(5707) = 6*MStop2(2)**2 - 4*Cd(1755)*MStop2(3)

	Cd(5708) = 9*MStop2(2) - 20*MStop2(3)

	Cd(5709) = 9*MStop2(3) + 5*MStop2(4)

	Cd(5710) = Cd(5708)*MStop2(2)**2 + Cd(5709)*MStop2(3)**2

	Cd(5711) = 17*MStop2(2) + MStop2(3)

	Cd(5712) = 45*MStop2(2) - 67*MStop2(3)

	Cd(5713) = Cd(5712)*MStop2(2)**2 + Cd(5711)*MStop2(3)**2

	Cd(5714) = 9*MStop2(2) + 4*MStop2(3)

	Cd(5716) = 7*MStop2(2) - 9*MStop2(3)

	Cd(5717) = 15*MStop2(2) - MStop2(3)

	Cd(5718) = 5*Cd(5716)*MStop2(2)**2 + Cd(5717)*MStop2(3)**2

	Cd(5723) = 9*MStop2(2)**2 - Cd(5689)*MStop2(3)

	Cd(5724) = 7*MStop2(2) - 3*MStop2(3)

	Cd(5726) = 15*MStop2(2)**2 - 2*Cd(1755)*MStop2(3)

	Cd(5727) = 5*MTy2 + Cd(5714)

	Cd(5728) = Cd(5726) + MTy2*Cd(5727)

	Cd(5729) = 7*MStop2(2) + 2*MStop2(3)

	Cd(5731) = 2*MTy2 - 3*MStop2(4)

	Cd(563) = AtC*UUStop(1,1) - 2*At*UUStop(1,2)

	Cd(564) = MUEC*Cd(563) + At*MUE*UUStop(1,1)

	Cd(565) = 2*MUE*UUStop(1,1) - MUEC*UUStop(1,2)

	Cd(566) = AtC*Cd(565) - At*MUE*UUStop(1,2)

	Cd(5768) = 14*MTy2 + 3*MStop2(1)

	Cd(5772) = Cd(990)*MStop2(2)**2 - Cd(535)*MStop2(3)**2

	Cd(567) = Cd(564)*UUStopC(1,1) - Cd(566)*UUStopC(1,2)

	Cd(5779) = 5*MTy2 + 16*MStop2(1)

	Cd(569) = MStop2(3)**2 - 4*MTy2*MStop2(1)*UStop2(1,1)

	Cd(5796) = 15*MStop2(2) - 4*MStop2(3)

	Cd(5804) = MStop2(4) + 2*Cd(302)*UStop2(1,2)**2

	Cd(570) = Cd(569) + 4*MTy2*MStop2(1)*UStop2(1,1)

	Cd(5806) = 7*MStop2(3) + 2*MStop2(4)

	Cd(5807) = 4*MStop2(2)**2 - Cd(5806)*MStop2(3)

	Cd(5809) = 2*MGlmT2 - 5*MStop2(1)

	Cd(5813) = 25*MStop2(2) - 23*MStop2(3)

	Cd(5814) = 25*MStop2(1)**2 + Cd(5813)*MStop2(3)

	Cd(5815) = 95*MStop2(2) - 18*MStop2(3)

	Cd(571) = Cd(570)*UStop2(1,1) + MStop2(3)**2*UStop2(1,2)

	Cd(5816) = 95*MStop2(2)**2 - Cd(5815)*MStop2(3)

	Cd(5817) = 4*MTy2 + 6*MStop2(4)

	Cd(5820) = 13*MStop2(2) - MStop2(3)

	Cd(5821) = 13*MStop2(2)**2 - Cd(5820)*MStop2(3)

	Cd(5822) = 19*MStop2(2) - 16*MStop2(3)

	Cd(5823) = 19*MStop2(1)**2 + Cd(5822)*MStop2(3)

	Cd(5825) = 2*MStop2(3) - 5*MStop2(4)

	Cd(572) = Cd(571) + Cd(471)*MStop2(1)

	Cd(5826) = 10*MStop2(1)**2 - Cd(5825)*MStop2(3)

	Cd(5828) = 3*MStop2(1)**2 + Cd(1898)*MStop2(3)

	Cd(5831) = 20*MStop2(3) + 3*MStop2(4)

	Cd(5832) = 12*MStop2(2)**2 - 2*Cd(5831)*MStop2(3)

	Cd(5833) = 36*MStop2(3) - 23*MStop2(4)

	Cd(5834) = 46*MStop2(1)**2 - Cd(5833)*MStop2(3)

	Cd(573) = 8*CA2*Cd(572) + Cd(475)*MStop2(1)

	Cd(5836) = 7*MStop2(1)**2 + Cd(5184)*MStop2(3)

	Cd(5841) = 42*MStop2(2) - 5*MStop2(3)

	Cd(5844) = 10*MStop2(2)**2 - Cd(1022)*MStop2(3)

	Cd(5853) = 2*MTy2 + 7*MStop2(1)

	Cd(5861) = 10*MTy2 + 7*MStop2(4)

	Cd(5863) = 5*MTy2 + 4*MStop2(1)

	Cd(5869) = 6*MStop2(2)**2 + Cd(5207)*MStop2(3)

	Cd(5870) = -(3/8.D0*(Cd(3514)**2*Cd(3662)**2)) + MStop2(3)**4

	Cd(5872) = 9*MTy2 + 8*MStop2(1)

	Cd(577) = Cr(28)*MStop2(2) + MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(5878) = 21*MStop2(2) - 5*MStop2(3)

	Cd(5879) = 21*MStop2(2)**2 - Cd(5878)*MStop2(3)

	Cd(5884) = 16*MStop2(1) + MStop2(3)

	Cd(578) = MStop2(2) + 2*MTy2*UStop2(1,2)

	Cd(5887) = 4*MStop2(3) - 7*MStop2(4)

	Cd(5888) = 14*MStop2(1)**2 - Cd(5887)*MStop2(3)

	Cd(5889) = 8*MStop2(3) + 21*MStop2(4)

	Cd(5891) = 5*MStop2(2) - 7*MStop2(3)

	Cd(5892) = 5*MStop2(1)**2 + Cd(5891)*MStop2(3)

	Cd(5894) = 17*MStop2(3) + 24*MStop2(4)

	Cd(579) = 4*MTy2 + MStop2(3)

	Cd(5899) = 4*MStop2(3) - 15*MStop2(4)

	Cd(5900) = 30*MStop2(2)**2 + Cd(5899)*MStop2(3)

	Cd(5902) = 13*MStop2(3) + 6*MStop2(4)

	Cd(5905) = 3*MStop2(2) + 4*MStop2(3)

	Cd(580) = Cd(579)*MStop2(2) - 4*MTy2*MStop2(1)*UStop2(1,2)

	Cd(5906) = 3*MStop2(2)**2 - Cd(5905)*MStop2(3)

	Cd(5907) = 28*MStop2(2) - 25*MStop2(3)

	Cd(5912) = 10*MStop2(3) - 9*MStop2(4)

	Cd(5913) = 18*MStop2(1)**2 - Cd(5912)*MStop2(3)

	Cd(5915) = 8*MStop2(3) - 35*MStop2(4)

	Cd(581) = Cd(580) - 2*MTy2*MStop2(4)*UStop2(1,1)

	Cd(5916) = 70*MStop2(1)**2 - Cd(5915)*MStop2(3)

	Cd(5920) = 2*MStop2(3) + 5*MStop2(4)

	Cd(5921) = 10*MStop2(2)**2 - Cd(5920)*MStop2(3)

	Cd(582) = Cd(581)*UStop2(1,1) + Cd(578)*MStop2(3)*UStop2(1,2)

	Cd(5927) = 14*MStop2(2) + MStop2(3)

	Cd(5930) = 28*MStop2(3) + 17*MStop2(4)

	Cd(5931) = 34*MStop2(2)**2 - Cd(5930)*MStop2(3)

	Cd(5933) = 4*MTy2 + 3*MStop2(4)

	Cd(583) = At*AtC*Cd(577) + Cd(582)

	Cd(5945) = 6*MStop2(2)**2 - Cd(1569)*MStop2(3)

	Cd(584) = 4*Cd(583) - 2*Cd(425)*MStop2(4)

	Cd(5954) = 6*MTy2 + MStop2(4)

	Cd(585) = AtC*UCStop(1,3) + At*UCStopC(1,3)

	Cd(5957) = 2*MStop2(3) - MStop2(4)

	Cd(5958) = 2*MStop2(2)**2 + Cd(5957)*MStop2(3)

	Cd(6) = At*CA + MUEC*SA

	Cd(57) = 2*MGlpT2 - MStop2(1)

	Cd(586) = (-Cd(584) + 8*MTy*Cd(585)*MStop2(3))/MStop2(3)

	Cd(5969) = MTy2*Cd(5370) + 3*MStop2(1)*MStop2(2)

	Cd(5975) = -MGl2 + Cd(5954)

	Cd(587) = 2*Cd(424)*Cd(577) - 2*Cd(426)*MStop2(4)

	Cd(588) = (-Cd(587) + 4*MTy*Cd(13)*MStop2(3))/MStop2(3)

	Cd(589) = MStop2(3) + MStop2(1)*UStop2(1,1)

	Cd(5996) = 3*MStop2(2)**2 - Cd(5211)*MStop2(3)

	Cd(5997) = 50*MStop2(2) - 7*MStop2(3)

	Cd(5998) = 50*MStop2(2)**2 - Cd(5997)*MStop2(3)

	Cd(6003) = 23*MTy2 + 17*MStop2(1)

	Cd(6005) = 7*MTy2 - 3*MStop2(1)

	Cd(590) = Cd(589)*UStop2(1,1) + MStop2(2)*UStop2(1,2)**2

	Cd(591) = 4*MUE2*Cd(590) - 2*Cd(427)*MStop2(4)

	Cd(6018) = MStop2(2)**2 + MStop2(3)*MStop2(4)

	Cd(6024) = 22*MTy2 + 39*MStop2(1)

	Cd(6025) = 5*MStop2(2) + 36*MStop2(3)

	Cd(592) = CA2*Cd(586) + S2A*Cd(588)

	Cd(6029) = 7*MTy2 + 10*MStop2(1)

        Cd(593) = 1/(256.D0*Pi**3)*
     &   (Alfas2L*Cr(9)**2*(-(SA2*Cd(591)) + Cd(592)*MStop2(3)))/
     &    (SB2*SW2*MStop2(3))

	Cd(594) = 5*MStop2(1) + MStop2(3)

	Cd(6050) = 7*MStop2(2) + 3*MStop2(3)

	Cd(6052) = 25*MStop2(2) - 6*MStop2(3)

	Cd(6053) = 7*MStop2(2) + 9*MStop2(3)

	Cd(6054) = -5*Cd(6053)*MStop2(2)**2 + 3*Cd(6052)*MStop2(3)**2

	Cd(595) = 5*MStop2(2) - 3*MStop2(3)

	Cd(6058) = 7*MStop2(2)**2 - Cd(5207)*MStop2(3)

	Cd(6059) = 2*MStop2(3) + 3*MStop2(4)

	Cd(6060) = Cd(1124)*MStop2(2)**2 - Cd(6059)*MStop2(3)**2

	Cd(6062) = 7*MStop2(3) + 9*MStop2(4)

	Cd(6064) = 5*MStop2(2) + MStop2(3)

	Cd(6065) = 5*MStop2(2)**2 - Cd(6064)*MStop2(3)

	Cd(596) = 3 + UStop2(1,2)

	Cd(6066) = 2*MStop2(3) - 7*MStop2(4)

	Cd(6067) = 3*MStop2(2)**2 - Cd(6066)*MStop2(3)

	Cd(6069) = 7*MStop2(2) - 10*MStop2(3)

	Cd(6070) = 9*MStop2(3) + 11*MStop2(4)

        Cd(6071) = Cd(1124)*Cd(6069)*MStop2(2)**2 + 
     &   Cd(6070)*MStop2(3)**3

	Cd(6073) = 9*MStop2(2) - 19*MStop2(3)

	Cd(6074) = 15*MStop2(2)**2 - Cd(5920)*MStop2(3)

	Cd(597) = 2*Cd(596)*MStop2(1) + 4*MStop2(2)*UStop2(1,2)

	Cd(6076) = 9*MStop2(2) + MStop2(3)

	Cd(6077) = 5*MStop2(2) + 3*MStop2(3)

	Cd(6083) = 9*MStop2(2) + 23*MStop2(3)

	Cd(6084) = 101*MStop2(3) + 33*MStop2(4)

	Cd(6085) = 2*Cd(6083)*MStop2(2)**2 - Cd(6084)*MStop2(3)**2

	Cd(598) = Cd(597) - 2*Cd(595)*UStop2(1,1)

	Cd(599) = 2*UStop2(1,1)**3 - UStop2(1,2)**3

	Cd(6104) = (-8*MStop2(2) + 4*MStop2(3))/MStop2(3)

	Cd(600) = 2*Cd(599)*MStop2(2) + Cd(598)*UStop2(1,1)*UStop2(1,2)

	Cd(602) = 5*MStop2(1) + MStop2(3)

	Cd(6131) = 5*MTy2 - 4*MStop2(1)

	Cd(603) = MStop2(1) + MStop2(4)*UStop2(1,2)

	Cd(6141) = 13*MTy2 + 15*MStop2(1)

	Cd(604) = Cd(603) - MStop2(4)*UStop2(1,1)

	Cd(6148) = 16*MTy2 + 5*MStop2(1)

	Cd(605) = Cd(599)*MStop2(2) + 2*Cd(604)*UStop2(1,1)*UStop2(1,2)

	Cd(6157) = 11*MTy2 + 32*MStop2(1)

	Cd(6169) = 10*MStop2(1) + MStop2(3)

	Cd(6171) = 5*MTy2 + Cd(286)

	Cd(6172) = Cd(1216) + MTy2*Cd(6171)

	Cd(6174) = 3*MStop2(2)**2 - Cd(1141)*MStop2(3)

	Cd(607) = 4*MStop2(1) + MStop2(3)

	Cd(6176) = 3*MStop2(2)**2 - 2*Cd(2114)*MStop2(3)

	Cd(6177) = 15*MStop2(1)**2 - 2*Cd(5825)*MStop2(3)

	Cd(6179) = 21*MStop2(2) - 4*MStop2(3)

	Cd(6181) = 14*MStop2(2) - MStop2(3)

	Cd(6182) = 14*MStop2(2)**2 - Cd(6181)*MStop2(3)

	Cd(608) = 3*MStop2(1) + Cd(535)*UStop2(1,2)

	Cd(6191) = 2*MStop2(3) + 9*MStop2(4)

	Cd(6192) = 18*MStop2(2)**2 - Cd(6191)*MStop2(3)

	Cd(6203) = -2*MGl2 + Cd(5861)

	Cd(610) = 2*MStop2(1) + MStop2(2)*UStop2(1,2)

	Cd(6213) = 15*MStop2(2) - 2*MStop2(3)

	Cd(6214) = 15*MStop2(2)**2 - Cd(6213)*MStop2(3)

	Cd(6218) = 6*MTy2 - 5*MStop2(1)

	Cd(6225) = 3*MGlpT2 - 7*MStop2(1)

	Cd(613) = 2*Cr(85)*MStop2(1) + 3*MStop2(2)*UStop2(1,2)

	Cd(6240) = 8*MStop2(3) + MStop2(4)

	Cd(6241) = 2*MStop2(2)**2 - Cd(6240)*MStop2(3)

	Cd(6242) = MStop2(1)**2 + Cd(256)*MStop2(3)

	Cd(614) = MStop2(2) + Cd(594)*UStop2(1,2)

	Cd(6259) = 10*MStop2(1)**2 + Cd(6169)*MStop2(3)

	Cd(6260) = 3 - 4*UStop2(1,2)**2

	Cd(6261) = 25 - 8*UStop2(1,2)**2

	Cd(6262) = 15 - 8*UStop2(1,2)**2

	Cd(6263) = 5 - 8*UStop2(1,2)**2

	Cd(6265) = 7 - 6*UStop2(1,2)**2

	Cd(6266) = 7 - 10*UStop2(1,2)**2

	Cd(6268) = 1 - 2*UStop2(1,2)**2

	Cd(6269) = 3 + 2*UStop2(1,2)**2

	Cd(6270) = 5 + 2*UStop2(1,2)**2

	Cd(6271) = 3 - 10*UStop2(1,2)**2

	Cd(6275) = 1 + 6*UStop2(1,2)**2

	Cd(617) = 7*MStop2(2) - 4*MStop2(3)

	Cd(6277) = 5 + 6*UStop2(1,2)**2

	Cd(6280) = 1 - 4*UStop2(1,2)**2

	Cd(6281) = 1 + 4*UStop2(1,2)**2

	Cd(6282) = 5 + 12*UStop2(1,2)**2

	Cd(6284) = 3 + 4*UStop2(1,2)**2

	Cd(6285) = 7 + 20*UStop2(1,2)**2

	Cd(6289) = 4 - 7*UStop2(1,2)

	Cd(6290) = 10 - 9*UStop2(1,2)

	Cd(6291) = 2 - 9*UStop2(1,2)

	Cd(6293) = 9 - 14*UStop2(1,2)

	Cd(6294) = 4 - 9*UStop2(1,2)

	Cd(6297) = 9 - 5*UStop2(1,2)

	Cd(6299) = 9 - 11*UStop2(1,2)

	Cd(6302) = 26 - 20*UStop2(1,2)

	Cd(6303) = 8 - 11*UStop2(1,2)

	Cd(6304) = Cd(6303)*MStop2(2) + 9*MTy2*UStop2(1,1)

	Cd(6306) = 3 - 4*UStop2(1,2)

	Cd(6308) = 19 - 22*UStop2(1,2)

	Cd(6310) = 9 - 8*UStop2(1,2)

	Cd(6312) = 11 - 14*UStop2(1,2)

	Cd(6313) = 6 - 9*UStop2(1,2)

	Cd(6315) = 20 - 19*UStop2(1,2)

	Cd(6318) = 11 - 9*UStop2(1,2)

	Cd(6320) = 1 - 4*UStop2(1,2)

	Cd(6322) = 1 - 11*UStop2(1,2)

	Cd(6324) = 6 - 22*UStop2(1,2)

	Cd(6325) = 3 - 11*UStop2(1,2)

	Cd(622) = MStop2(1) + Cd(595)*UStop2(1,2)

	Cd(6326) = 17 - 19*UStop2(1,2)

	Cd(6333) = 16*MTy2 + 3*MStop2(1)

	Cd(6334) = 13*MGlmT2 - 9*MStop2(1)

	Cd(6338) = 8*MStop2(2) - 17*MStop2(3)

	Cd(6339) = 24*MStop2(2) - 17*MStop2(3)

	Cd(6341) = 5*MStop2(2) - 16*MStop2(3)

	Cd(624) = 2 - 6*UStop2(1,2)

	Cd(6348) = 7*MTy2 + 13*MStop2(1)

	Cd(6352) = 9*MStop2(1)**2 + Cd(5184)*MStop2(3)

	Cd(6353) = 9*MStop2(2)**2 - Cd(6181)*MStop2(3)

	Cd(6354) = 6*MStop2(2) + MStop2(3)

	Cd(6355) = -4*Cd(6354)*MStop2(1)**2 + 3*MStop2(3)**3

	Cd(6360) = 8*MGl2*MTy2 - MGlpTmSt4(2)

	Cd(6362) = 5*MGl2 + 6*MTy2

	Cd(626) = MStop2(3)**2 - MTy2*MStop2(1)*UStop2(1,1)

	Cd(627) = Cd(626) + MTy2*MStop2(1)*UStop2(1,1)

	Cd(628) = Cd(627)*UStop2(1,1) + MStop2(3)**2*UStop2(1,2)

	Cd(6393) = 16*MTy2 + 21*MStop2(1)

	Cd(6405) = 10*MTy2 - MStop2(1)

	Cd(6415) = 35*MTy2 + 36*MStop2(1)

	Cd(6421) = 5*MTy2 - 8*MStop2(1)

	Cd(6425) = 5*MGlmT2 - 13*MStop2(1)

	Cd(632) = 2*MStop2(2) + MStop2(4)*UStop2(1,2)

	Cd(6433) = 7*MGl2 + 10*MTy2

	Cd(6435) = 4*MGl2 + 3*MStop2(1)

	Cd(6440) = 4*MTy2 + 19*MStop2(1)

	Cd(6445) = MTy2 + 15*MStop2(1)

	Cd(634) = UStop2(1,1)**3 - 2*UStop2(1,2)**2

	Cd(6452) = 8*MStop2(2)**2 - MStop2(3)*MStop2(4)

        Cd(635) = 4*Cd(634)*MStop2(2)*UStop2(1,1) - 
     &   2*Cd(632)*UStop2(1,1)**3

	Cd(6461) = 25*MTy2 + 37*MStop2(1)

	Cd(6464) = 13*MTy2 + 35*MStop2(1)

	Cd(637) = 8*MStop2(2) - MStop2(3)

	Cd(640) = MStop2(2) - MStop2(1)*UStop2(1,1)

	Cd(6507) = 8*MTy2 + 13*MStop2(1)

	Cd(6510) = 5*MTy2 + 11*MStop2(1)

	Cd(6511) = 3*MStop2(2) - 13*MStop2(3)

	Cd(641) = 2*Cd(640)*MStop2(2) - 2*MStop2(1)**2*UStop2(1,2)

	Cd(6522) = MTy2 + 12*MStop2(1)

	Cd(642) = Cd(641) - 2*MStop2(2)*MStop2(3)*UStop2(1,1)

	Cd(6529) = MTy2 + 7*MStop2(1)

	Cd(6531) = 3*(Cd(529)**2 - MStop2(2)**2)

	Cd(6532) = MTy2 - 16*MStop2(1)

	Cd(6534) = 15*MStop2(2)**2 - Cd(5796)*MStop2(3)

	Cd(6543) = 11*MTy2 - 15*MStop2(1)

	Cd(6548) = 11*MTy2 + 19*MStop2(1)

	Cd(6553) = 10*MTy2 + 19*MStop2(1)

	Cd(6554) = 5*MStop2(1)**2 + MStop2(3)*MStop2(4)

	Cd(645) = 3*MStop2(2) - 4*MStop2(3)

	Cd(6557) = MTy2 - 34*MStop2(1)

	Cd(6558) = 34*MStop2(3) - 5*MStop2(4)

	Cd(63) = 3 - 2*Log(MStop2(2)/Q)

	Cd(646) = 4*UStop2(1,1)**2 + Cd(596)*UStop2(1,2)

	Cd(6573) = 10*MTy2**2 - 3*MStop2(1)**2

	Cd(647) = -(Cd(646)*MStop2(2)) + Cd(645)*UStop2(1,1)*UStop2(1,2)

	Cd(648) = Cr(85)*MStop2(2) - Cd(256)*UStop2(1,1)

	Cd(649) = 2*MStop2(2)*UStop2(1,1)**2 + Cd(648)*UStop2(1,2)

	Cd(6605) = 21*MStop2(2) - 17*MStop2(3)

	Cd(650) = 3*MStop2(2) + Cr(26)*UStop2(1,2)

	Cd(6607) = 9*MTy2 + 11*MStop2(1)

	Cd(651) = Cd(650) - 3*MStop2(2)*UStop2(1,1)

	Cd(6622) = 19*MTy2 + 4*MStop2(1)

	Cd(6625) = 10*MTy2 + 9*MStop2(1)

	Cd(652) = MStop2(2) + 2*MStop2(3)

	Cd(6628) = 22*MTy2 + 31*MStop2(1)

	Cd(6630) = 31*MTy2 - 7*MStop2(1)

	Cd(653) = Cd(652) - MStop2(2)

	Cd(654) = 4*Cd(640)*MStop2(2) - 2*Cd(523)*UStop2(1,2)

	Cd(6653) = 27*MGl2 + 4*MTy2

	Cd(655) = Cd(654) - 4*MStop2(2)*MStop2(3)*UStop2(1,1)

	Cd(6656) = 3/4.D0*(Cd(3662)**2*MStop2(4)) + Cd(535)*MStop2(3)**2

	Cd(6664) = 12*MTy2 + MStop2(3)

        Cd(656) = -2*MUEC*Cd(655)*UCStopC(1,3)**2 + 
     &   2*MUE*Cd(649)*MStop2(3)*UStop2(1,2)

        Cd(657) = 2*MUE*Cd(651)*UCStop(1,3)**2 + 
     &   MUEC*Cd(647)*UStop2(1,2)

	Cd(6681) = 40*MStop2(2) - 9*MStop2(3)

        Cd(658) = 2*MUE*MStop2(3)*UCStop(1,3) + 
     &   MUEC*Cd(653)*UCStopC(1,3)

	Cd(659) = AtC*Cd(657) + 4*MTy*Cd(658)*UStop2(1,2)

	Cd(660) = At*Cd(656) - Cd(659)*MStop2(3)

	Cd(6713) = 6*MTy2 + Cd(607)

	Cd(661) = Cd(596) - 3*UStop2(1,1)

	Cd(6716) = -4*MTy2 + Cd(607)

	Cd(662) = Cd(661)*UStop2(1,1)**3 + 4*UStop2(1,1)*UStop2(1,2)**3

	Cd(6730) = MTy2 + 4*MStop2(1)

	Cd(6733) = MStop2(2) + 2*MStop2(3)

        Cd(663) = Cd(662)*MStop2(2) + 
     &   4*MStop2(3)*UStop2(1,1)**2*UStop2(1,2)**2

	Cd(6740) = 15*MTy2 + 4*MStop2(1)

	Cd(664) = MStop2(2) + Cd(286)*UStop2(1,2)

	Cd(6747) = 11*MTy2 + 39*MStop2(1)

	Cd(6750) = 22*MTy2 - 103*MStop2(1)

	Cd(6751) = 6*MStop2(2) - 13*MStop2(3)

	Cd(6754) = 4*MTy2 - 7*MStop2(1)

	Cd(665) = Cd(664) - MStop2(2)*UStop2(1,1)

	Cd(65) = 2*M_3*MTy*MGlpTmSt2(2) + MGlpTmSt4(2)*UCStop(1,3)

	Cd(666) = 1 - 7*UStop2(1,2)

	Cd(6773) = 8*MTy2 + 17*MStop2(1)

	Cd(667) = Cd(666)*MStop2(2) + 4*MStop2(1)*UStop2(1,2)

	Cd(6788) = 3*MTy2 - 2*MStop2(3)

	Cd(6802) = 54*MStop2(2) - 35*MStop2(3)

	Cd(6803) = 21*MStop2(2) - 5*MStop2(3)

	Cd(670) = Cr(85)*MStop2(2) - 4*MStop2(1)*UStop2(1,2)

	Cd(6806) = 8*MTy2 - 53*MStop2(1)

	Cd(6807) = 25*MStop2(2) - 21*MStop2(3)

	Cd(671) = UStop2(1,1)**2 - 2*Cr(85)*UStop2(1,2)

	Cd(6821) = 14*MTy2 - 11*MStop2(1)

	Cd(6824) = 15*MTy2 + 13*MStop2(1)

	Cd(6828) = 7*MStop2(2) - 15*MStop2(3)

	Cd(673) = MTy2 + MStop2(3)

	Cd(674) = MStop2(2) - 4*MStop2(3)

	Cd(6848) = 7*MTy2 + 19*MStop2(1)

	Cd(6851) = 11*MTy2 - 26*MStop2(1)

	Cd(6852) = 5*MStop2(2) - 6*MStop2(3)

	Cd(675) = Cd(673)*MStop2(2) - MTy2*Cd(674)*UStop2(1,2)

	Cd(6861) = 8*MTy2 - 7*MStop2(1)

	Cd(676) = Cd(675) - MTy2*Cr(26)*UStop2(1,1)

	Cd(6871) = 8*MTy2**2 + 11*MStop2(1)**2

	Cd(677) = Cd(676)*UStop2(1,1) + Cd(578)*MStop2(3)*UStop2(1,2)

	Cd(6876) = 5*MTy2 - 9*MStop2(1)

	Cd(6894) = 17*MTy2 + 12*MStop2(1)

	Cd(6905) = MTy2 + Cd(594)

	Cd(6909) = 10*MTy2 + MStop2(4)

	Cd(6911) = MStop2(2)**2 - Cr(2)*MStop2(3)

	Cd(6912) = 2*MStop2(2)**2 - Cd(990)*MStop2(3)

	Cd(6913) = -6*MTy2**2 + Cd(6912)

	Cd(6916) = 5*MTy2**4 + MStop2(1)**2*MStop2(2)**2

	Cd(6917) = MTy2**2*Cd(6912) - 3*Cd(6916)

	Cd(685) = -(MTy2*UStop(1,2)) + AtC*MTy*UStop(1,1)*UStop2(1,2)

	Cd(6958) = 15*MGl2 - 4*Cd(1074)

	Cd(6962) = 13*MGl2 - 21*MStop2(1)

	Cd(67) = 4*MTy2 + MStop2(1)

	Cd(689) = 2*CA*MTy*UStop(2,1) + Cd(6)*UStop(2,2)

	Cd(7001) = 9*MGlmT2 - MStop2(1)

	Cd(690) = Cd(5)*UStop(2,1) + 2*CA*MTy*UStop(2,2)

	Cd(7007) = MGlmT2 - 3*MStop2(1)

	Cd(691) = Cd(689)*UStopC(1,1) + Cd(690)*UStopC(1,2)

        Cd(692) = -2*Cd(689)*Cd(691)*MStop2(1) + 
     &   2*CA2*MStop2(2)*MStop2(3)*UStop(1,1)

	Cd(7032) = 3*MStop2(2) - 5*MStop2(3)

	Cd(693) = 2*CA*MTy*UStop(1,1) + Cd(6)*UStop(1,2)

	Cd(694) = Cd(5)*UStop(1,1) + 2*CA*MTy*UStop(1,2)

	Cd(7050) = 4*MTy2 + 9*MStop2(1)

	Cd(7051) = 3*MStop2(1)**2 + MStop2(3)**2

	Cd(695) = MUE*SA*Cr(91) + 4*CA*MTy*UUStop(3,2)

	Cd(7054) = 3*MGlmT2 - 5*MStop2(1)

	Cd(7059) = MTy2 + 18*MStop2(1)

	Cd(7062) = 5*MStop2(2) - 33*MStop2(3)

	Cd(68) = 3*MGlpT2 - MStop2(4)

	Cd(696) = MUE*SA*UUStop(3,1) + At*CA*UUStop(3,2)

	Cd(7071) = 4*MTy2 - 21*MStop2(1)

	Cd(698) = At*Cr(91) + 4*MTy*UUStop(3,1)

	Cd(7088) = 3*MStop2(2)**2 + MStop2(3)**2

	Cd(699) = CA*Cd(698) + MUEC*SA*Cr(91)

	Cd(7103) = 49*MStop2(2) - 34*MStop2(3)

	Cd(7105) = 30*MStop2(1) + MStop2(3)

	Cd(7106) = 49*MStop2(1)**2 + Cd(7105)*MStop2(3)

	Cd(7108) = 17*MStop2(2) - 2*MStop2(3)

	Cd(7114) = Cd(5184)*MStop2(1)**2 - 3*MStop2(3)**2*MStop2(4)

	Cd(7116) = 79*MStop2(3) + 98*MStop2(4)

	Cd(7118) = 49*MStop2(2) - 102*MStop2(3)

	Cd(7119) = 101*MStop2(3) + 150*MStop2(4)

	Cd(7120) = 5*Cd(7118)*MStop2(2)**2 + Cd(7119)*MStop2(3)**2

	Cd(7121) = 42*MStop2(2)**2 - 8*Cd(1570)*MStop2(3)

	Cd(7122) = 35*MStop2(2) - 78*MStop2(3)

	Cd(7123) = 7*MStop2(3) + 6*MStop2(4)

	Cd(702) = Cd(693)*UStopC(1,1) + Cd(694)*UStopC(1,2)

	Cd(7124) = Cd(7122)*MStop2(2)**2 + 3*Cd(7123)*MStop2(3)**2

	Cd(7126) = 63*MStop2(2) - 30*MStop2(3)

	Cd(7127) = 7*MStop2(3) + 33*MStop2(4)

	Cd(7128) = 105*MStop2(2)**2 - 2*Cd(7127)*MStop2(3)

	Cd(7130) = 68*MStop2(2) + 5*MStop2(3)

	Cd(7131) = 21*MStop2(2) - 46*MStop2(3)

	Cd(7132) = 3*Cd(7131)*MStop2(2)**2 + Cd(7130)*MStop2(3)**2

	Cd(7134) = 89*MStop2(3) + 86*MStop2(4)

	Cd(7141) = 3*MStop2(3) - 5*MStop2(4)

	Cd(7146) = 2*(2*MTy2 + 3*Cd(1272))*MStop2(3)

	Cd(7147) = MTy2 + 8*MStop2(1)

	Cd(705) = -(AtC*UCStop(1,3)) + 4*MTy*UStop2(1,1)

	Cd(7171) = 46*MTy2 + 99*MStop2(1)

	Cd(7172) = 20*MStop2(2) + 51*MStop2(3)

	Cd(708) = AtC*S2A + MUE*SA2

	Cd(7189) = 29*MTy2 - 9*MStop2(1)

	Cd(709) = -(Cd(708)*UCStop(1,3)) + 2*MTy*S2A*UStop2(1,1)

        Cd(710) = MUE*Cd(709)*UCStopC(3,4) + 
     &   Cd(689)*Cd(690)*UUStopC(3,3)

	Cd(7205) = MStop2(3) + 18*MStop2(4)

	Cd(7206) = 72*MStop2(2)**2 - 2*Cd(7205)*MStop2(3)

	Cd(7208) = 24*MStop2(2) - 7*MStop2(3)

	Cd(7209) = 24*MStop2(2)**2 - Cd(7208)*MStop2(3)

	Cd(7211) = 8*MStop2(2)**2 - Cd(1121)*MStop2(3)

	Cd(7213) = 6*MStop2(2) - 7*MStop2(3)

	Cd(7214) = 41*MStop2(3) + 18*MStop2(4)

	Cd(7215) = 36*MStop2(2)**2 - Cd(7214)*MStop2(3)

	Cd(7217) = 19*MStop2(2) + MStop2(3)

	Cd(7220) = 5*MStop2(3)**4 - 9*MStop2(4)**4

	Cd(7221) = 36*MStop2(2)**2 - Cd(7205)*MStop2(3)

	Cd(712) = CA2*MStop2(3)*UCStopC(3,1) - Cd(694)**2*UUStopC(3,2)

        Cd(713) = CA2*MStop2(3)*UCStopC(3,1) + 
     &   Cd(693)*Cd(694)*UUStopC(3,3)

	Cd(7234) = 50*MStop2(2) - 11*MStop2(3)

	Cd(7235) = 50*MStop2(2)**2 - Cd(7234)*MStop2(3)

	Cd(7238) = 30*MStop2(2) - MStop2(3)

	Cd(7239) = 30*MStop2(2)**2 - Cd(7238)*MStop2(3)

	Cd(714) = Cd(712) - Cd(713)

        Cd(715) = 2*(Cd(714)*MStop2(2) + 
     &     Cd(690)*Cd(694)*MStop2(1)*UUStopC(2,2))

	Cd(7258) = 4*MTy2 + 15*MStop2(1)

	Cd(7273) = 2*MStop2(3) + 15*MStop2(4)

	Cd(7274) = 30*MStop2(2)**2 - Cd(7273)*MStop2(3)

	Cd(7277) = 34*MStop2(2)**2 + Cd(5669)*MStop2(3)

	Cd(7278) = 15*MStop2(2)**2 - Cd(5717)*MStop2(3)

        Cd(719) = 2*Conjugate(dY33eps(0))*MStop2(1) + 
     &   Conjugate(dY33eps(-1))*MStop2(3)

	Cd(7301) = -MGl2 + Cd(5349)

	Cd(720) = 2*AtC*CA2 + MUE*S2A

	Cd(7307) = 19*MStop2(1)**2 + Cd(5190)*MStop2(3)

	Cd(7308) = 74*MStop2(2) - 17*MStop2(3)

	Cd(7309) = 74*MStop2(2)**2 - Cd(7308)*MStop2(3)

	Cd(7310) = 155*MStop2(2) - 126*MStop2(3)

	Cd(7311) = 155*MStop2(1)**2 + Cd(7310)*MStop2(3)

	Cd(7312) = 92*MStop2(2) - 7*MStop2(3)

	Cd(7315) = 50*MTy2 + 17*MStop2(4)

	Cd(7317) = 11*MStop2(2) - 2*MStop2(3)

	Cd(7318) = 11*MStop2(1)**2 + Cd(7317)*MStop2(3)

	Cd(7319) = 4*MStop2(3) - MStop2(4)

	Cd(7320) = 2*MStop2(2)**2 + Cd(7319)*MStop2(3)

	Cd(7322) = 36*MStop2(3) - 19*MStop2(4)

	Cd(7323) = 38*MStop2(2)**2 + Cd(7322)*MStop2(3)

	Cd(7325) = 3*MStop2(3) - 10*MStop2(4)

	Cd(7328) = 12*MStop2(3) + 91*MStop2(4)

	Cd(7329) = 182*MStop2(2)**2 - Cd(7328)*MStop2(3)

	Cd(7331) = 17*MStop2(2) - 18*MStop2(3)

	Cd(7332) = 17*MStop2(1)**2 + Cd(7331)*MStop2(3)

	Cd(7333) = 9*MStop2(2) + 2*MStop2(3)

	Cd(7334) = 9*MStop2(1)**2 + Cd(7333)*MStop2(3)

	Cd(7336) = 7*MStop2(3) - 10*MStop2(4)

	Cd(7338) = 68*MStop2(3) + 179*MStop2(4)

	Cd(7339) = 358*MStop2(2)**2 - Cd(7338)*MStop2(3)

	Cd(7341) = 26*MStop2(2) - 17*MStop2(3)

	Cd(7342) = 26*MStop2(1)**2 + Cd(7341)*MStop2(3)

	Cd(7346) = 29*MTy2 + 5*MStop2(1)

	Cd(7355) = 7*MTy2 - 74*MStop2(1)

	Cd(7358) = 48*MStop2(2) - MStop2(3)

	Cd(7359) = 48*MStop2(2)**2 - Cd(7358)*MStop2(3)

	Cd(7360) = 30*MStop2(2) + 7*MStop2(3)

	Cd(70) = 2*MGlpT2 - MStop2(2)

	Cd(7363) = 11*MStop2(2) + 3*MStop2(3)

	Cd(7373) = 23*MStop2(2) - 9*MStop2(3)

	Cd(7374) = 23*MStop2(2)**2 - Cd(7373)*MStop2(3)

	Cd(7379) = 19*MTy2 + 15*MStop2(1)

	Cd(7386) = 22*MTy2 + 15*MStop2(4)

        Cd(729) = 2*Conjugate(dY33eps(1))*MStop2(1) + 
     &   Conjugate(dY33eps(0))*MStop2(3)

	Cd(7393) = 5*MTy2 - 27*MStop2(1)

	Cd(7398) = 31*MTy2 - 10*MStop2(1)

	Cd(7403) = 9*MTy2 + 106*MStop2(1)

	Cd(7409) = 14*MTy2 + 59*MStop2(1)

	Cd(7412) = 41*MStop2(3) + 30*MStop2(4)

	Cd(7425) = 28*MTy2 - 33*MStop2(1)

	Cd(7426) = 2*MTy2 + 33*MStop2(1)

	Cd(7429) = 22*MStop2(2) + MStop2(3)

	Cd(7433) = MGl2 - 3*Cd(5225)

	Cd(734) = Cd(689)*Cd(691) - CA2*MStop2(3)*UStop(1,1)

	Cd(7442) = 44*MTy2 + 81*MStop2(1)

	Cd(7443) = 5*MStop2(2) + 18*MStop2(3)

	Cd(735) = Cd(694)*UUStopC(3,2) + Cd(693)*UUStopC(3,3)

	Cd(7456) = MTy2 - 30*MStop2(1)

        Cd(736) = Cd(694)*
     &   (-2*Cd(735)*MStop2(2) + 2*Cd(690)*MStop2(1)*UUStopC(2,2))

	Cd(7464) = 19*MStop2(3) + 30*MStop2(4)

	Cd(737) = Cd(690)*UUStopC(3,2) + Cd(689)*UUStopC(3,3)

	Cd(738) = Cd(690)*Cd(737) + CA2*MStop2(3)*UCStop(3,1)

	Cd(7482) = 5*MStop2(3) + 18*MStop2(4)

	Cd(7483) = 36*MStop2(2)**2 - Cd(7482)*MStop2(3)

	Cd(7485) = 260*MStop2(2) - 43*MStop2(3)

	Cd(7486) = 260*MStop2(2)**2 - Cd(7485)*MStop2(3)

	Cd(7490) = 9*MTy2 - 8*MStop2(1)

	Cd(7502) = 4*MTy2 - 7*MStop2(4)

	Cd(7504) = 9*MStop2(2) + 7*MStop2(3)

	Cd(7505) = 9*MStop2(2)**2 - Cd(7504)*MStop2(3)

	Cd(7506) = 12*MStop2(3) - 7*MStop2(4)

	Cd(7507) = 14*MStop2(2)**2 + Cd(7506)*MStop2(3)

	Cd(7509) = MTy2 + Cd(5349)

	Cd(7510) = 2*Cd(7509) + 9*MStop2(4)

	Cd(7511) = 2*MGl2 + Cd(7502)

	Cd(7528) = 70*MStop2(2) - 11*MStop2(3)

	Cd(743) = Cd(694)*Cd(735) + Cd(693)*Cd(702)*UStopC(2,1)

	Cd(7534) = 7*MTy2 - 6*MStop2(1)

	Cd(7536) = 13*MStop2(2) - 3*MStop2(3)

	Cd(7537) = 13*MStop2(2)**2 - Cd(7536)*MStop2(3)

	Cd(7559) = 3*MTy2**2 - MStop2(1)**2

	Cd(7563) = MTy2 + 13*MStop2(1)

	Cd(7568) = 2*MTy2 - 9*MStop2(1)

	Cd(7572) = 22*MStop2(3) + MStop2(4)

	Cd(7573) = Cd(7572)*MStop2(2)**2 - 7*Cd(535)*MStop2(3)**2

	Cd(7579) = 5*MTy2 + 32*MStop2(1)

	Cd(7603) = 4*MGlmT2 - 7*MStop2(1)

	Cd(7607) = 35*MStop2(2) - 29*MStop2(3)

	Cd(7608) = 35*MStop2(1)**2 + Cd(7607)*MStop2(3)

	Cd(7609) = 27*MStop2(2) - 2*MStop2(3)

	Cd(7613) = 21*MStop2(2) - 22*MStop2(3)

	Cd(7614) = 21*MStop2(1)**2 + Cd(7613)*MStop2(3)

	Cd(7616) = 18*MStop2(2) + MStop2(3)

	Cd(7617) = 18*MStop2(2)**2 - Cd(7616)*MStop2(3)

	Cd(7619) = 6*MStop2(2)**2 + Cd(1810)*MStop2(3)

	Cd(7621) = 11*MStop2(3) + 8*MStop2(4)

	Cd(7623) = 40*MStop2(3) - 27*MStop2(4)

	Cd(7624) = 54*MStop2(1)**2 - Cd(7623)*MStop2(3)

	Cd(7662) = 70*MStop2(2) - 57*MStop2(3)

	Cd(7663) = 70*MStop2(1)**2 + Cd(7662)*MStop2(3)

	Cd(7670) = 10*MGl2 + 7*MStop2(4)

	Cd(757) = 2*MTy - AtC*UCStop(1,3)

	Cd(7672) = 5*MStop2(1)**2 + Cd(1325)*MStop2(3)

	Cd(7675) = 15*MStop2(2) - 8*MStop2(3)

	Cd(7676) = 15*MStop2(2)**2 - Cd(7675)*MStop2(3)

	Cd(7677) = 9*MStop2(2) - 10*MStop2(3)

	Cd(758) = -(MTy*UCStopC(1,3)) + AtC*UStop2(1,1)*UStop2(1,2)

	Cd(7682) = 26*MStop2(2) + 7*MStop2(3)

	Cd(7686) = 2*MStop2(1)**2 + Cd(2114)*MStop2(3)

	Cd(7687) = 14*MStop2(1)**2 - Cd(6066)*MStop2(3)

	Cd(7689) = 6*MStop2(2)**2 - Cd(2985)*MStop2(3)

	Cd(7693) = 16*MStop2(3) - 9*MStop2(4)

	Cd(7694) = 18*MStop2(1)**2 - Cd(7693)*MStop2(3)

	Cd(7709) = 8*MTy2 + 3*MStop2(4)

	Cd(763) = Cd(169)*UStop2(1,1) - 2*UStop2(1,2)

	Cd(7757) = 118*MStop2(2) - MStop2(3)

	Cd(7758) = 118*MStop2(2)**2 - Cd(7757)*MStop2(3)

	Cd(766) = Cd(426) + AtC*MUEC*UStop2(1,1)*UStop2(1,2)

	Cd(7765) = 22*MTy2 + 27*MStop2(4)

	Cd(767) = Cd(766)*MStop2(4) + MTy*MUEC*MStop2(3)*UCStopC(1,3)

	Cd(7774) = 26*MTy2 + 57*MStop2(1)

	Cd(7775) = 5*MStop2(2) + 78*MStop2(3)

	Cd(7781) = MTy2 + 16*MStop2(1)

	Cd(769) = UCStop(3,1)**2 + UCStopC(3,1)**2

	Cd(7802) = 5*MGl2 + 7*Cr(1)

	Cd(7803) = 135*MStop2(2) - 32*MStop2(3)

	Cd(7804) = 35*Cr(2)*MStop2(2)**2 - Cd(7803)*MStop2(3)**2

	Cd(7805) = 9*MStop2(2) - 23*MStop2(3)

	Cd(770) = 4*AtC*MTy*Cr(75)*UCStop(1,3)

	Cd(7811) = MStop2(3) - 7*MStop2(4)

	Cd(7813) = 7*MStop2(2)**2 - Cd(6066)*MStop2(3)

	Cd(7816) = Cd(3504)*MStop2(2)**2 - Cd(5920)*MStop2(3)**2

	Cd(7817) = 53*MStop2(2) + 8*MStop2(3)

	Cd(7818) = 9*MStop2(2) + 43*MStop2(3)

	Cd(7819) = Cd(7818)*MStop2(2)**2 - Cd(7817)*MStop2(3)**2

	Cd(7820) = 157*MStop2(2) - 158*MStop2(3)

	Cd(7821) = 45*MStop2(2) + 302*MStop2(3)

	Cd(7822) = Cd(7821)*MStop2(1)**2 + Cd(7820)*MStop2(3)**2

	Cd(7825) = MStop2(2) + 13*MStop2(3)

	Cd(7826) = 36*MStop2(1) + MStop2(3)

	Cd(7827) = 5*MStop2(2) + 37*MStop2(3)

	Cd(7828) = Cd(7827)*MStop2(1)**2 + Cd(7826)*MStop2(3)**2

	Cd(7829) = 2*MStop2(3) - 11*MStop2(4)

	Cd(7830) = 3*MStop2(2)**2 - Cd(7829)*MStop2(3)

	Cd(7832) = 9*MStop2(2) + 5*MStop2(3)

	Cd(7849) = (-4*MStop2(2) + 2*MStop2(3))/MStop2(3)

	Cd(774) = MStop2(2) - 5*MStop2(3)

	Cd(775) = MStop2(3) - 2*Cd(774)*UStop2(1,2)

	Cd(7865) = 3*MTy2 + 25*MStop2(1)

	Cd(776) = Cd(596)*MStop2(1) - Cd(519)*MStop2(2)

	Cd(7884) = 3*MTy2 + 17*MStop2(1)

	Cd(778) = 2*MStop2(2) - 3*MStop2(3)

	Cd(7894) = 8*MStop2(2) + 9*MStop2(3)

	Cd(7895) = Cd(7894)*MStop2(2)**2 - 3*Cd(535)*MStop2(3)**2

	Cd(779) = 1 + 4*UStop2(1,2)

	Cd(7900) = 9*MTy2 - 32*MStop2(1)

	Cd(7901) = 7*MStop2(1)**2 - MStop2(3)*MStop2(4)

	Cd(7907) = 13*MTy2 - 3*MStop2(1)

	Cd(780) = Cd(779)*MStop2(1) - 3*MStop2(2)

	Cd(7914) = 14*MStop2(2)**2 - Cd(5927)*MStop2(3)

	Cd(781) = 5 - UStop2(1,2)

	Cd(7923) = 16*MStop2(3) + 9*MStop2(4)

	Cd(7924) = 18*MStop2(2)**2 - Cd(7923)*MStop2(3)

	Cd(782) = 3 - 7*UStop2(1,2)

	Cd(7930) = 3*MTy2 - 23*MStop2(1)

	Cd(7939) = 21*MStop2(2) + MStop2(3)

	Cd(7940) = 21*MStop2(2)**2 - Cd(7939)*MStop2(3)

	Cd(7947) = 3*MGlmT2 - 7*MStop2(1)

	Cd(7) = dMSfsq1133eps(0) - 2*dMSfsq1133eps(1)

	Cd(786) = 6*MStop2(1) + MStop2(3)

	Cd(7969) = 10 - 11*UStop2(1,2)

	Cd(7970) = 10 - 3*UStop2(1,2)

	Cd(7971) = 7*Cd(523) - 2*Cd(7970)*MStop2(1)*MStop2(2)

	Cd(7972) = 4 + UStop2(1,2)

	Cd(7973) = 9 - 4*UStop2(1,2)

	Cd(7976) = 8 - 6*UStop2(1,2)

	Cd(7977) = Cd(7976)*MStop2(2) + 6*MTy2*UStop2(1,1)

	Cd(787) = 1 + 2*UStop2(1,2)

	Cd(7980) = 8 - 3*UStop2(1,2)

	Cd(7981) = 7 + 3*UStop2(1,2)

	Cd(7983) = 5 - 6*UStop2(1,2)

	Cd(7985) = 7 + UStop2(1,2)

	Cd(7986) = 27 - 25*UStop2(1,2)

	Cd(788) = 3 - 6*UStop2(1,2)

	Cd(7991) = 13 - 14*UStop2(1,2)

	Cd(7993) = 7 - 6*UStop2(1,2)

	Cd(7995) = 24 - 25*UStop2(1,2)

	Cd(7999) = 11 - 5*UStop2(1,2)

	Cd(8000) = 13 - 25*UStop2(1,2)

	Cd(8007) = 11*MGlmT2 - 7*MStop2(1)

	Cd(790) = 2*UStop2(1,1)**3 + Cr(85)*UStop2(1,2)**2

	Cd(8011) = 17*MStop2(2) - 48*MStop2(3)

	Cd(8013) = 5*MStop2(2) - 12*MStop2(3)

	Cd(791) = Cd(790)*MStop2(4) - Cd(786)*UStop2(1,1)**2*UStop2(1,2)

	Cd(8023) = 4*MTy2 - 27*MStop2(1)

	Cd(8025) = 17*MStop2(1) + MStop2(3)

	Cd(8026) = 15*MStop2(1)**2 + Cd(8025)*MStop2(3)

	Cd(8027) = 20*MStop2(2) - 21*MStop2(3)

	Cd(8028) = 20*MStop2(2) + 11*MStop2(3)

	Cd(8029) = 2*Cd(8028)*MStop2(1)**2 + Cd(8027)*MStop2(3)**2

	Cd(8033) = 53*MTy2 - 18*MStop2(1)

	Cd(8034) = MTy2 - 15*MStop2(1)

	Cd(8040) = 7*MTy2 - 5*MStop2(1)

	Cd(8052) = 4*MGlpT2**2 - MStop2(1)**2

	Cd(8056) = 5*MGl2 + 26*MTy2

	Cd(795) = MUE*UCStop(1,3)**2 + MUEC*UStop2(1,1)*UStop2(1,2)

	Cd(8074) = 11*MTy2 - 5*MStop2(1)

	Cd(8084) = 19*MTy2 + 40*MStop2(1)

	Cd(8092) = 7*MTy2 - 4*MStop2(1)

	Cd(8094) = 7*MGlmT2 - 11*MStop2(1)

	Cd(799) = MStop2(2) - Cd(256)*UStop2(1,2)

	Cd(8102) = MTy2 + Cd(44)

	Cd(8103) = 2*Cd(8102) - MStop2(4)

	Cd(8104) = 7*MGl2 + 18*MTy2

	Cd(800) = 2*UStop2(1,1) - UStop2(1,2)

	Cd(8111) = 20*MTy2 + 27*MStop2(1)

	Cd(8112) = 26*MStop2(3) - 7*MStop2(4)

	Cd(8116) = 5*MStop2(3)**2 + MStop2(4)**2

	Cd(802) = MStop2(2) + 4*MStop2(3)

	Cd(8132) = 23*MTy2 + 39*MStop2(1)

	Cd(8133) = 8*MStop2(1)**2 + MStop2(3)*MStop2(4)

	Cd(8136) = 17*MTy2 + 33*MStop2(1)

	Cd(8137) = 5*MStop2(2) - 27*MStop2(3)

	Cd(803) = 2*Cr(1) + 3*MStop2(3)*UStop2(1,2)

	Cd(8140) = 7*MTy2**2 + 33*MStop2(1)**2

	Cd(804) = MStop2(3) - MStop2(1)*UStop2(1,2)

	Cd(8151) = 21*MTy2 + 8*MStop2(1)

	Cd(8153) = 14*MTy2**2 - 5*MStop2(1)**2

	Cd(8154) = 4*MStop2(2) + 21*MStop2(3)

	Cd(806) = 2 - 3*UStop2(1,2)

	Cd(807) = 4 - UStop2(1,2)

	Cd(8184) = 8*MStop2(2) - 19*MStop2(3)

	Cd(8185) = Cd(8184)*MStop2(2)**2 + Cd(7536)*MStop2(3)**2

	Cd(8186) = 6*MStop2(3) - 7*MStop2(4)

	Cd(8187) = 8*MStop2(2)**2 + 2*Cd(8186)*MStop2(3)

	Cd(8189) = (10*MTy2 + 3*Cd(1272))*MStop2(3)

	Cd(809) = 3 + 2*UStop2(1,2)

	Cd(811) = 2*Cd(512)*MStop2(1) - Cd(523)*UStop2(1,2)

        Cd(812) = Cd(607)*UStop2(1,1)**3 - 
     &   Cr(75)*MStop2(4)*UStop2(1,2)**2

	Cd(8232) = MTy2*Cd(315) - MStop2(2)*MStop2(4)

	Cd(8233) = MGl2 - 2*MStop2(4)

	Cd(813) = Cd(812)*MStop2(3) + Cd(811)*UStop2(1,1)**2

	Cd(814) = 2*MStop2(2) - 3*MStop2(3)

	Cd(815) = 2*MStop2(3) + Cd(814)*UStop2(1,2)

	Cd(8263) = (MTy2*Cd(181) + Cd(2567)*MStop2(1))*MStop2(3)

        Cd(8292) = Cd(1430)*Cr(26)*MStop2(2)**2 - 
     &   2*Cd(1755)*MStop2(3)**3

	Cd(8310) = 7*MTy2 - MStop2(1)

	Cd(8311) = 17*MStop2(2) - 3*MStop2(3)

	Cd(8314) = 19*MGlmT2 - 7*MStop2(1)

	Cd(8324) = 7*MGlmT2 + MStop2(1)

	Cd(8329) = 22*MStop2(3) - 7*MStop2(4)

	Cd(8332) = 5*MGlmT2 - 7*MStop2(1)

	Cd(824) = 8*MTy2*UStop2(1,2)*UUStop(1,1) - MStop2(3)*UUStop(2,1)

	Cd(825) = Cd(769) - UCStop(3,1)**2

	Cd(8360) = 30*MStop2(2)**2 - Cd(5104)*MStop2(3)

	Cd(8363) = 7*MStop2(2) - 4*MStop2(3)

	Cd(8364) = 7*MStop2(1)**2 + Cd(8363)*MStop2(3)

	Cd(8365) = 10*MStop2(2)**2 + Cd(5361)*MStop2(3)

	Cd(8399) = (2*MGl2 + 3*Cd(1272))*MStop2(3)

	Cd(8400) = 8*MStop2(3) + 15*MStop2(4)

	Cd(8401) = 30*MStop2(2)**2 - Cd(8400)*MStop2(3)

	Cd(831) = 5*MStop2(2)**2 - Cd(552)*MStop2(3)

	Cd(8422) = 5*MTy2**2 - MStop2(1)**2

	Cd(8425) = 17*MStop2(2) - MStop2(3)

	Cd(832) = 3*MStop2(2)**2 - Cd(535)*MStop2(3)

	Cd(8426) = 17*MStop2(2)**2 - Cd(8425)*MStop2(3)

	Cd(8428) = 7*MStop2(2) + MStop2(3)

	Cd(8429) = 7*MStop2(2)**2 - Cd(8428)*MStop2(3)

	Cd(8430) = 11*MStop2(2) - 5*MStop2(3)

	Cd(8431) = 11*MStop2(2)**2 - Cd(8430)*MStop2(3)

	Cd(8432) = 6*MStop2(2)**2 + Cd(5414)*MStop2(3)

	Cd(834) = 22*MStop2(2) - 5*MStop2(3)

	Cd(835) = 22*MStop2(2)**2 - Cd(834)*MStop2(3)

	Cd(836) = Cd(596)*MStop2(1) + Cd(519)*MStop2(2)

	Cd(837) = 1 + 8*UStop2(1,2)

	Cd(8479) = MTy2**2 - MStop2(1)*MStop2(2)

	Cd(838) = 2 + 20*UStop2(1,2)

	Cd(8502) = 2*MTy2 + 7*MStop2(4)

	Cd(8503) = 4*MTy2**2 - MStop2(4)**2

	Cd(8511) = MGl2 - 9*MTy2

	Cd(8513) = MGlmT2*MGlpT2 + 2*MTy2*MStop2(2)

	Cd(8514) = MStop2(2)**4 + Cd(8513)*MStop2(3)**2

	Cd(8515) = 2*MGl2 + 3*MTy2

	Cd(8516) = MGlpT2*Cd(8511) + 2*Cd(8515)*MStop2(2)

	Cd(8517) = 4*MGlmT2*MGlpT2**2 - Cd(8516)*MStop2(2)

	Cd(8518) = MStop2(2)**5 + Cd(8517)*MStop2(3)**2

	Cd(8519) = 4*MGlpT2 - MStop2(2)

	Cd(842) = 6*MStop2(1)**2 + Cd(607)*MStop2(3)

	Cd(843) = 2*MStop2(1) + MStop2(3)

	Cd(8536) = -4*MGl2 + Cd(8502)

	Cd(8537) = 4*Cd(8503) + MGl2*Cd(8536)

	Cd(844) = Cd(843)*MStop2(3) - MStop2(1)*MStop2(2)*UStop2(1,1)

	Cd(845) = Cd(779)*MStop2(1) + 2*MStop2(2)*UStop2(1,1)

	Cd(8560) = (10*MTy2 + 3*Cd(843))*MStop2(3)

	Cd(846) = 2*Cd(845)*MStop2(1) + 6*MStop2(2)**2

	Cd(8568) = 4*MStop2(3) + 5*MStop2(4)

	Cd(847) = 8*MStop2(2)**2 - 3*MStop2(3)*MStop2(4)

	Cd(848) = 10*MStop2(1)**2 + Cd(607)*MStop2(3)

	Cd(8584) = 2 + 3*UStop2(1,2)

	Cd(8585) = 4 + 3*UStop2(1,2)

	Cd(8588) = 1 - 5*UStop2(1,2)

	Cd(8589) = 7 + 5*UStop2(1,2)

	Cd(8590) = 2*MTy2*Cd(8588) + Cd(8589)*MStop2(4)

	Cd(849) = 2*(Cd(847) - Cd(848)*UStop2(1,2))

	Cd(8591) = 3 + 7*UStop2(1,2)

	Cd(850) = Cd(849)*UStop2(1,1) - 4*Cd(844)*UStop2(1,2)

	Cd(852) = MStop2(2)**2 - Cd(535)*MStop2(3)

	Cd(8623) = MGlpT2 + 4*MStop2(1)

	Cd(854) = 6 - 12*UStop2(1,2)

	Cd(856) = 2*MStop2(1)**2 - Cd(607)*MStop2(3)

	Cd(8662) = 2*MStop2(2) + 7*MStop2(3)

	Cd(8663) = 2*MStop2(1)**2 + Cd(4267)*MStop2(3)

	Cd(857) = 2*MStop2(1)*MStop2(2) + Cd(856)*UStop2(1,2)

	Cd(8672) = 2*MStop2(2)**2 + 11*Cd(1807)*MStop2(3)

	Cd(8674) = 3*MStop2(2)**2 + Cd(4471)*MStop2(3)

	Cd(8675) = 11*MStop2(2) - 9*MStop2(3)

	Cd(8676) = 2*MStop2(2)**2 + Cd(8675)*MStop2(3)

	Cd(8678) = 11*MStop2(2) - 15*MStop2(3)

	Cd(8679) = 2*MStop2(2)**2 + Cd(8678)*MStop2(3)

	Cd(8690) = 4*MStop2(1)**2 + Cd(786)*MStop2(3)

	Cd(860) = 4*MStop2(3) + 9*MStop2(4)

	Cd(861) = 17*MStop2(2)**2 - Cd(860)*MStop2(3)

	Cd(8713) = -(MGl2*Cd(1592)) + Cd(45)*Cd(2763)

	Cd(8717) = 2*MTy2 - MStop2(2)

	Cd(8718) = MStop2(1)**2 - Cd(8717)*MStop2(2)

	Cd(862) = 4*MStop2(1)*MStop2(2) + 2*Cd(861)*UStop2(1,2)

	Cd(8722) = Cd(45)*Cd(607) - MGl2*MStop2(1)

	Cd(8730) = 3*MStop2(2) + MStop2(3)

	Cd(8731) = -(Cd(8730)*MStop2(1)) + Cd(2036)*MStop2(2)

	Cd(8770) = MStop2(2) - MTy2*UStop2(1,1)

	Cd(8771) = Cd(169)*MStop2(2) + 3*MTy2*UStop2(1,1)

	Cd(8774) = 5 - 7*UStop2(1,2)

	Cd(869) = 7*MStop2(1)**2 + Cd(786)*MStop2(3)

	Cd(8792) = 6*MStop2(1)**2 + Cd(4267)*MStop2(3)

	Cd(870) = MStop2(2)**2 + Cd(869)*UStop2(1,2)

	Cd(8802) = 4*MStop2(2) - 13*MStop2(3)

	Cd(8803) = 4*MTy2 + 7*MStop2(3)

	Cd(8805) = 8*MStop2(2) - 21*MStop2(3)

	Cd(8809) = UUStop(1,1) - 2*UUStop(2,1)

	Cd(8810) = -(Cd(8809)*MStop2(1)) + MStop2(2)*UUStop(1,1)

	Cd(8811) = -(Cd(3705)*MStop2(1)) + MStop2(2)*UUStop(1,1)

	Cd(8812) = 122*UUStop(1,1) + UUStop(2,1)

	Cd(8813) = 148*UUStop(1,1) + UUStop(2,1)

	Cd(8815) = 3*UUStop(1,1) - 2*UUStop(2,1)

	Cd(8816) = -(Cd(8815)*MStop2(1)) + 3*MStop2(2)*UUStop(1,1)

	Cd(8817) = 37*UUStop(1,1) + UUStop(2,1)

	Cd(8818) = 38*UUStop(1,1) + UUStop(2,1)

	Cd(872) = MStop2(3) + 3*MStop2(4)

	Cd(8822) = 71*UUStop(1,1) - 2*UUStop(2,1)

	Cd(8823) = 59*UUStop(1,1) - UUStop(2,1)

	Cd(8825) = 8*UUStop(1,1) + UUStop(2,1)

	Cd(8826) = 12*UUStop(1,1) + UUStop(2,1)

	Cd(8828) = 65*UUStop(1,1) - 3*UUStop(2,1)

	Cd(873) = 13*MStop2(2)**2 - 2*Cd(872)*MStop2(3)

	Cd(8829) = 130*UUStop(1,1) - 3*UUStop(2,1)

	Cd(8832) = 5*MStop2(3) - 2*MStop2(4)

	Cd(8833) = 6*MStop2(2)**2 - Cd(8832)*MStop2(3)

	Cd(8834) = 15*UUStop(1,1) - 2*UUStop(2,1)

	Cd(8835) = 21*UUStop(1,1) - UUStop(2,1)

	Cd(8837) = 2*UUStop(1,1) + UUStop(2,1)

	Cd(874) = 4*MStop2(2) - MStop2(3)

	Cd(875) = 3*MStop2(2)**2 - MStop2(3)*MStop2(4)

	Cd(8851) = 3*MGl2 - 4*MTy2

        Cd(8854) = 2*Cd(3718)*MStop2(3)*UUStop(1,2) - 
     &   Cd(44)*MStop2(1)*UUStop(2,2)

	Cd(876) = 4*Cd(875) + 6*Cd(831)*UStop2(1,2)

	Cd(877) = Cd(876) - 2*Cd(873)*UStop2(1,1)

	Cd(8874) = 3*MGl2 + Cd(3842)

	Cd(878) = Cd(832) + Cd(874)*MStop2(1)*UStop2(1,2)

	Cd(8892) = 1/2.D0*(Cd(872)*MStop2(3)) + MStop2(2)**2

	Cd(8895) = 9*MTy2 - 4*MStop2(1)

	Cd(8897) = 4*Cd(49)*Cr(1) - MGl2*MStop2(2)

	Cd(880) = 13*MStop2(2) - 6*MStop2(3)

	Cd(8901) = MGl2*MStop2(2) - 2*Cd(49)*MStop2(4)

	Cd(881) = MStop2(2)**2 - Cd(880)*MStop2(3)

	Cd(8909) = 3*MStop2(2)**2 + Cd(4662)*MStop2(3)

	Cd(8910) = MStop2(2)**2 + Cd(872)*MStop2(3)

	Cd(8911) = MStop2(2)**2 + Cd(1757)*MStop2(3)

	Cd(8913) = MStop2(2)**2 + Cd(1569)*MStop2(3)

	Cd(8919) = 2*Cd(44)*MStop2(2) - Cd(843)*MStop2(3)

	Cd(8920) = 4*MTy2**2 - Cd(315)*MStop2(4)

	Cd(8921) = Cd(8920)*MStop2(2) + 2*MTy2*MStop2(3)**2

	Cd(8923) = 11*MTy2 + 2*MStop2(1)

	Cd(883) = 4 - 6*UStop2(1,2)

	Cd(8943) = 2*MGl2 - 5*MTy2

	Cd(8944) = 2*MTy2*Cd(2617) - Cd(8943)*MStop2(1)

	Cd(8) = dMSfsq1133eps(-1) - 2*dMSfsq1133eps(0)

	Cd(86) = MGl2 + Cd(44)

	Cd(8953) = -(Cd(49)*Cd(874)) + MGl2*MStop2(2)

	Cd(8962) = 9*MTy2 + 4*MStop2(1)

	Cd(8964) = MTy2*Cd(579) - Cd(1141)*MStop2(2)

	Cd(8965) = 12*MTy2 + 7*MStop2(3)

	Cd(8966) = -2*MTy2*Cd(8965) + 2*Cd(3318)*MStop2(2)

	Cd(887) = 17*MStop2(2) - 5*MStop2(3)

	Cd(8975) = 7*MTy2 + 8*MStop2(1)

	Cd(888) = 17*MStop2(2)**2 - Cd(887)*MStop2(3)

	Cd(890) = 6 + 4*UStop2(1,2)

	Cd(8996) = 2*MGl2*Cd(45) - MTy2*Cd(5294)

	Cd(9005) = Cd(169)*MStop2(2) + MTy2*UStop2(1,1)

	Cd(9006) = Cd(806)*MStop2(2) + 3*MTy2*UStop2(1,1)

	Cd(9027) = 3*MTy2 - 4*MStop2(1)

	Cd(9030) = 9*MTy2 + 2*MStop2(1)

	Cd(894) = 4*MStop2(2)**2 - 6*MStop2(3)*MStop2(4)

	Cd(9041) = 4*MTy2 - 7*MStop2(3)

	Cd(9042) = 4*MStop2(2) + 9*MStop2(3)

	Cd(87) = Cd(86)*UStop(1,1) - 2*M_3*MTy*UStop(1,2)

	Cd(895) = 2*Cd(172)*MStop2(1)**2 - 2*MStop2(2)**2*UStop2(1,2)

	Cd(9043) = 8*MStop2(2) + 13*MStop2(3)

	Cd(9048) = -(Cd(8809)*MStop2(2)) + MStop2(1)*UUStop(1,1)

	Cd(9049) = -(Cd(3705)*MStop2(2)) + MStop2(1)*UUStop(1,1)

	Cd(9050) = 52*UUStop(1,1) + UUStop(2,1)

	Cd(9052) = -(Cd(8815)*MStop2(2)) + 3*MStop2(1)*UUStop(1,1)

	Cd(9054) = 25*UUStop(1,1) - UUStop(2,1)

	Cd(9055) = 17*UUStop(1,1) + 2*UUStop(2,1)

	Cd(9058) = 8*UUStop(1,1) - UUStop(2,1)

	Cd(9060) = 41*UUStop(1,1) - UUStop(2,1)

	Cd(9061) = 83*UUStop(1,1) - 2*UUStop(2,1)

        Cd(897) = 2*MStop2(1)*MStop2(2)*UCStopC(3,1)**2 + 
     &   8*MTy2*MStop2(4)*UStop2(1,1)*UStop2(1,2)**2

	Cd(9063) = 106*UUStop(1,1) - 3*UUStop(2,1)

	Cd(9064) = 89*UUStop(1,1) - 3*UUStop(2,1)

	Cd(9067) = 6*MStop2(2)**2 + Cd(1121)*MStop2(3)

	Cd(9069) = 15*UUStop(1,1) - UUStop(2,1)

	Cd(9070) = 39*UUStop(1,1) - 2*UUStop(2,1)

	Cd(898) = 4*MStop2(1)**2 + MStop2(3)*MStop2(4)

	Cd(9075) = MTy2 + 5*MStop2(2)

	Cd(9076) = 3*MTy2 + 7*MStop2(2)

        Cd(9080) = 2*Cd(3969)*MStop2(3)*UUStop(1,2) + 
     &   Cd(48)*MStop2(2)*UUStop(2,2)

        Cd(899) = 16*MTy2*Cd(523)*UUStop(1,1) - 
     &   MStop2(1)*MStop2(2)*MStop2(3)*UUStop(2,1)

	Cd(9085) = 9*MTy2 + MStop2(2)

	Cd(900) = Cd(899) - 8*MTy2*Cd(898)*UStop2(1,2)*UUStop(1,1)

	Cd(9102) = 3*MGl2 + Cd(4852)

        Cd(901) = Cd(900)*UStop2(1,2) - 
     &   MStop2(1)*MStop2(2)*MStop2(3)*UUStop(2,1)

	Cd(902) = 3*MStop2(2)**2 - Cd(307)*MStop2(3)

	Cd(903) = 8*MTy2*UUStop(1,1) - 3*MStop2(3)*UUStop(2,1)

        Cd(904) = Cd(903)*MStop2(1)*MStop2(2) + 
     &   8*MTy2*Cd(902)*UStop2(1,2)*UUStop(1,1)

        Cd(905) = 4*MTy2*UStop2(1,1)**4 + 
     &   MStop2(3)*UCStopC(3,1)**2*UStop2(1,2)

	Cd(9143) = Cd(86)*UStop(1,1) - M_3*MTy*UStop(1,2)

	Cd(9145) = Cd(44) + 3*Cd(2200)

	Cd(906) = Cd(901) - Cd(904)*UStop2(1,1)

	Cd(9157) = Cd(48) + 3*Cd(2254)

	Cd(9170) = Cd(86)*UStop(1,1) + M_3*MTy*UStop(1,2)

	Cd(9181) = 3*Cd(86)*UStop(1,1) - 2*M_3*MTy*UStop(1,2)

	Cd(9187) = Cd(3718)*MStop2(1) - Cd(1592)*MStop2(2)

	Cd(9197) = 96*MStop2(2) - 65*MStop2(3)

	Cd(9199) = 48*MStop2(2) - 37*MStop2(3)

	Cd(9200) = 24*MStop2(2) - 11*MStop2(3)

	Cd(9201) = 12*MStop2(2) + 5*MStop2(3)

	Cd(911) = Cd(693)*UUStopC(3,1) + Cd(694)*UUStopC(3,3)

	Cd(9203) = 48*MStop2(2) - 67*MStop2(3)

	Cd(9204) = 72*MStop2(2) - 25*MStop2(3)

	Cd(9206) = 12*MStop2(2) - 29*MStop2(3)

	Cd(9207) = 12*MStop2(2) - 13*MStop2(3)

	Cd(9210) = 4*MStop2(2) - 17*MStop2(3)

        Cd(9211) = MStop2(3)*
     &   (Cd(9210)*MStop2(1)**2 + 4*Cd(47)*MStop2(3)*UStop2(1,2))

        Cd(912) = -2*Cd(689)*Cd(911)*MStop2(2) + 
     &   2*CA2*MStop2(3)**2*UStop2(1,1)

	Cd(913) = Cd(694)*UUStopC(3,2) + Cd(693)*UUStopC(3,4)

	Cd(9229) = 11*MTy2 + 14*MStop2(1)

        Cd(914) = -2*Cd(690)*Cd(913)*MStop2(2) + 
     &   2*CA2*MStop2(3)**2*UStop2(1,2)

	Cd(9241) = -(MTy2*Cd(786)) + MStop2(1)*MStop2(3)

	Cd(89) = 2*MTy2 - MStop2(4)

	Cd(915) = Cd(912) + Cd(914)

	Cd(9242) = 5*MTy2 + 8*MStop2(1)

	Cd(9245) = 11*MTy2 + 3*MStop2(1)

	Cd(9250) = -3*MGl2 + Cd(9245)

	Cd(9283) = 2*UUStop(1,1) - 5*UUStop(2,1)

	Cd(9284) = Cd(9283)*MStop2(1) - 2*MStop2(2)*UUStop(1,1)

	Cd(9285) = 12*UUStop(1,1) - UUStop(2,1)

	Cd(9287) = 8*UUStop(1,1) - 5*UUStop(2,1)

	Cd(9288) = -(Cd(9287)*MStop2(1)) + 8*MStop2(2)*UUStop(1,1)

	Cd(9289) = 6*UUStop(1,1) - 5*UUStop(2,1)

	Cd(9290) = -(Cd(9289)*MStop2(1)) + 6*MStop2(2)*UUStop(1,1)

	Cd(9291) = 4*UUStop(1,1) - 5*UUStop(2,1)

	Cd(9292) = -(Cd(9291)*MStop2(1)) + 4*MStop2(2)*UUStop(1,1)

	Cd(9293) = 10*UUStop(1,1) + UUStop(2,1)

	Cd(9294) = 148*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9295) = 16*UUStop(1,1) + UUStop(2,1)

	Cd(9296) = 32*UUStop(1,1) + UUStop(2,1)

	Cd(9297) = Cd(9295)*MStop2(1) - Cd(9296)*MStop2(2)

	Cd(9299) = 48*UUStop(1,1) - 2*UUStop(2,1)

	Cd(9300) = Cd(9299)*MStop2(1) - 2*Cd(3950)*MStop2(2)

	Cd(921) = Cd(693)**2*UUStopC(1,1) + Cd(694)**2*UUStopC(1,2)

	Cd(9301) = 24*UUStop(1,1) - UUStop(2,1)

	Cd(9303) = Cd(3938)*MStop2(1) - Cd(3682)*MStop2(2)

	Cd(9304) = 6*UUStop(1,1) + UUStop(2,1)

	Cd(9305) = 24*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9308) = 32*UUStop(1,1) - UUStop(2,1)

	Cd(9310) = 60*UUStop(1,1) + UUStop(2,1)

	Cd(9311) = -2*Cd(9310)*MStop2(1) + 2*Cd(3682)*MStop2(2)

	Cd(9312) = 76*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9315) = 4*UUStop(1,1) - UUStop(2,1)

	Cd(9316) = 16*UUStop(1,1) - UUStop(2,1)

	Cd(9317) = Cd(9315)*MStop2(1) - Cd(9316)*MStop2(2)

	Cd(9318) = 22*UUStop(1,1) + UUStop(2,1)

	Cd(9320) = 8*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9324) = 3*UUStop(1,3) + UUStop(2,3)

	Cd(9325) = Cd(9324)*MStop2(1) - 3*MStop2(2)*UUStop(1,3)

	Cd(9326) = 3*UUStop(1,3) + 2*UUStop(2,3)

	Cd(9327) = Cd(9326)*MStop2(1) - 3*MStop2(2)*UUStop(1,3)

	Cd(9328) = 3*UUStop(1,3) + 4*UUStop(2,3)

	Cd(9329) = Cd(9328)*MStop2(1) - 3*MStop2(2)*UUStop(1,3)

	Cd(9330) = 9*UUStop(1,3) + 4*UUStop(2,3)

	Cd(9331) = Cd(9330)*MStop2(1) - 9*MStop2(2)*UUStop(1,3)

	Cd(9332) = 33*UUStop(1,3) - 4*UUStop(2,3)

	Cd(9333) = 9*UUStop(1,3) + 2*UUStop(2,3)

	Cd(9334) = 7*UUStop(1,3) - UUStop(2,3)

	Cd(9335) = UUStop(1,3) + UUStop(2,3)

	Cd(9336) = Cd(9334)*MStop2(1) + Cd(9335)*MStop2(2)

	Cd(9338) = 28*UUStop(1,3) - UUStop(2,3)

	Cd(9339) = 4*UUStop(1,3) + UUStop(2,3)

	Cd(9340) = Cd(9338)*MStop2(1) + Cd(9339)*MStop2(2)

	Cd(90) = MGl2 - 2*MStop2(1)

	Cd(9341) = 2*M_3*(Cd(9340) + 32*MGl2*UUStop(1,3))

	Cd(9342) = 7*UUStop(1,3) + 2*UUStop(2,3)

	Cd(9343) = 35*UUStop(1,3) + 4*UUStop(2,3)

	Cd(9345) = 5*UUStop(1,3) + UUStop(2,3)

	Cd(9346) = 11*UUStop(1,3) - 6*UUStop(2,3)

	Cd(9348) = 2*UUStop(1,3) - 3*UUStop(2,3)

	Cd(9349) = 28*UUStop(1,3) + 6*UUStop(2,3)

	Cd(9350) = Cd(9349)*MStop2(1) + 2*Cd(9348)*MStop2(2)

	Cd(9352) = 9*MStop2(3) + 20*MStop2(4)

	Cd(9353) = 24*MStop2(2)**2 - Cd(9352)*MStop2(3)

	Cd(9354) = UUStop(1,3) - UUStop(2,3)

	Cd(9355) = Cd(9345)*MStop2(1) + Cd(9354)*MStop2(2)

	Cd(9356) = 11*UUStop(1,3) + 2*UUStop(2,3)

	Cd(9357) = 5*UUStop(1,3) + 4*UUStop(2,3)

	Cd(9363) = MGl2**2 - 3*MTy2**2

	Cd(928) = Cd(689)*UUStopC(3,1) + Cd(690)*UUStopC(3,4)

        Cd(929) = -4*Cd(693)*Cd(928)*MStop2(1) + 
     &   4*CA2*MStop2(3)**2*UStop2(1,2)

	Cd(9392) = 8*MTy2*UUStop(1,3) + MStop2(3)*UUStop(2,3)

	Cd(9406) = MTy2**2 + Cd(286)*MStop2(1)

	Cd(9407) = 12*MTy2*Cd(9406) - 6*Cd(4243)*MStop2(3)

	Cd(932) = 2*CA2*MStop2(2)*UStop(2,2) + Cd(690)**2*UStopC(2,2)

	Cd(9416) = 84*MStop2(2) - 31*MStop2(3)

	Cd(933) = CA2*MStop2(2)*UStop(2,1) + Cd(689)*Cd(690)*UStopC(2,2)

	Cd(9427) = 8*MTy2*UUStop(1,1) - MStop2(3)*UUStop(2,1)

        Cd(934) = Cd(693)*Cd(928)*MStop2(2) + 
     &   Cd(933)*MStop2(3)*UStopC(2,1)

	Cd(9435) = 8*MTy2*UUStop(1,2) + MStop2(3)*UUStop(2,2)

	Cd(91) = Cd(44)**2 + MGl2*Cd(90)

	Cd(935) = Cd(932)*UStopC(2,2) + Cd(689)**2*UUStopC(2,1)

	Cd(9472) = 48*MStop2(2) - 31*MStop2(3)

	Cd(9473) = 24*MStop2(2) + 19*MStop2(3)

	Cd(9474) = 36*MStop2(2) - 47*MStop2(3)

	Cd(9479) = 12*MStop2(2) + 17*MStop2(3)

	Cd(9485) = 7*MTy2 - MStop2(2)

	Cd(9499) = 4*MTy2 - 11*MStop2(1)

	Cd(9500) = 19*MTy2 - 14*MStop2(1)

	Cd(9510) = MTy2*Cd(1226) + MStop2(2)*MStop2(3)

	Cd(9515) = 11*MTy2 + 3*MStop2(2)

	Cd(9517) = -(MGl2*Cd(3233)) + 4*Cd(9510)

	Cd(9520) = -3*MGl2 + Cd(9515)

	Cd(945) = 1 - 3*Cr(15)

	Cd(9560) = 84*UUStop(1,1) + UUStop(2,1)

	Cd(947) = 8*MStop2(1)**2 + Cd(607)*MStop2(3)

	Cd(9562) = -(Cd(9287)*MStop2(2)) + 8*MStop2(1)*UUStop(1,1)

	Cd(9563) = -(Cd(9289)*MStop2(2)) + 6*MStop2(1)*UUStop(1,1)

	Cd(9564) = -(Cd(9291)*MStop2(2)) + 4*MStop2(1)*UUStop(1,1)

	Cd(9565) = -(Cd(9283)*MStop2(2)) + 2*MStop2(1)*UUStop(1,1)

	Cd(9566) = 52*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9567) = 194*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9569) = Cd(9301)*MStop2(1) - Cd(3950)*MStop2(2)

	Cd(948) = MStop2(1)**2 + Cd(947)*UStop2(1,2)

	Cd(9571) = 24*UUStop(1,1) - 5*UUStop(2,1)

	Cd(9572) = 42*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9575) = 34*UUStop(1,1) - UUStop(2,1)

	Cd(9577) = -2*Cd(3938)*MStop2(1) + 2*Cd(8826)*MStop2(2)

	Cd(9578) = 4*UUStop(1,1) + 5*UUStop(2,1)

	Cd(9579) = 98*UUStop(1,1) + 5*UUStop(2,1)

	Cd(949) = 2*MStop2(2)**2 - MStop2(3)**2

	Cd(9582) = 20*UUStop(1,1) - UUStop(2,1)

	Cd(9583) = Cd(9058)*MStop2(1) - Cd(9582)*MStop2(2)

	Cd(9584) = 34*UUStop(1,1) + UUStop(2,1)

	Cd(9586) = 32*UUStop(1,1) - 5*UUStop(2,1)

	Cd(9590) = 3*UUStop(1,3) - UUStop(2,3)

	Cd(9591) = -(Cd(9590)*MStop2(2)) + 3*MStop2(1)*UUStop(1,3)

	Cd(9592) = 3*UUStop(1,3) - 2*UUStop(2,3)

	Cd(9593) = -(Cd(9592)*MStop2(2)) + 3*MStop2(1)*UUStop(1,3)

	Cd(9594) = 9*UUStop(1,3) - 4*UUStop(2,3)

	Cd(9595) = -(Cd(9594)*MStop2(2)) + 9*MStop2(1)*UUStop(1,3)

	Cd(9596) = 3*UUStop(1,3) - 4*UUStop(2,3)

	Cd(9597) = -(Cd(9596)*MStop2(2)) + 3*MStop2(1)*UUStop(1,3)

	Cd(9598) = 9*UUStop(1,3) - 2*UUStop(2,3)

	Cd(9599) = 33*UUStop(1,3) + 4*UUStop(2,3)

	Cd(9600) = 7*UUStop(1,3) + UUStop(2,3)

	Cd(951) = MStop2(1)**2 + Cd(607)*MStop2(3)

	Cd(9601) = Cd(9354)*MStop2(1) + Cd(9600)*MStop2(2)

	Cd(9603) = 4*UUStop(1,3) - UUStop(2,3)

	Cd(9604) = 28*UUStop(1,3) + UUStop(2,3)

	Cd(9605) = Cd(9603)*MStop2(1) + Cd(9604)*MStop2(2)

	Cd(9606) = 2*M_3*(Cd(9605) + 32*MGl2*UUStop(1,3))

	Cd(9607) = 35*UUStop(1,3) - 4*UUStop(2,3)

	Cd(9608) = 7*UUStop(1,3) - 2*UUStop(2,3)

	Cd(9610) = 14*UUStop(1,3) - 3*UUStop(2,3)

	Cd(952) = MStop2(1)**2 + Cd(951)*UStop2(1,2)

	Cd(9611) = 4*UUStop(1,3) + 6*UUStop(2,3)

	Cd(9612) = Cd(9611)*MStop2(1) + 2*Cd(9610)*MStop2(2)

	Cd(9613) = 5*UUStop(1,3) - UUStop(2,3)

	Cd(9614) = 11*UUStop(1,3) + 6*UUStop(2,3)

	Cd(9617) = 9*MStop2(3) + 4*MStop2(4)

	Cd(9618) = 24*MStop2(2)**2 - Cd(9617)*MStop2(3)

	Cd(9619) = Cd(9335)*MStop2(1) + Cd(9613)*MStop2(2)

	Cd(9620) = 5*UUStop(1,3) - 4*UUStop(2,3)

	Cd(953) = 12*MStop2(2)*MStop2(3) - 4*Cd(952)*UStop2(1,2)

	Cd(9621) = 11*UUStop(1,3) - 2*UUStop(2,3)

	Cd(954) = Cd(948) - 2*MStop2(1)**2*UStop2(1,1)

	Cd(9649) = 14*MTy2 + MStop2(1)

	Cd(956) = 12*MStop2(1)*MStop2(3) - Cd(953)*UStop2(1,1)

	Cd(9655) = 13*MTy2 + 16*MStop2(1)

	Cd(9672) = 13*MTy2 - 7*MStop2(1)

	Cd(9675) = MTy2**2 + Cd(535)*MStop2(2)

	Cd(9676) = 6*(2*MTy2*Cd(9675) + Cd(4591)*MStop2(3))

	Cd(959) = 3 + 4*UStop2(1,2)

	Cd(9684) = 11*MStop2(3) + 6*MStop2(4)

	Cd(9693) = 7*MTy2 - 8*MStop2(1)

	Cd(962) = 12*MStop2(2)*MStop2(3) - 4*Cd(949)*UStop2(1,2)**2

	Cd(9716) = 2*MTy*SA + Cd(368)*UCStop(1,3)

	Cd(9717) = Cd(368)*(UStop2(1,1)**2 + UStop2(1,2)**2)

	Cd(9719) = 2*MTy*SA - Cd(368)*UCStop(1,3)

	Cd(9722) = 2*MTy*SA + Cd(369)*UCStopC(1,3)

	Cd(9723) = Cd(369)*(UStop2(1,1)**2 + UStop2(1,2)**2)

	Cd(9725) = 2*MTy*SA - Cd(369)*UCStopC(1,3)

	Cd(9730) = MUE*S2A - 2*AtC*SA2

	Cd(9731) = MUEC*S2A - 2*At*SA2

	Cd(9739) = Cd(9716) + Cd(369)*UCStopC(1,3)

	Cd(9741) = Cd(9719) - Cd(369)*UCStopC(1,3)

	Cd(9747) = Cd(9722) + Cd(368)*UCStop(1,3)

	Cd(9749) = Cd(9725) - Cd(368)*UCStop(1,3)

	Cd(966) = 12*MStop2(1)*MStop2(3) - Cd(962)*UStop2(1,1)

	Cd(9767) = -(S2A*Cd(473)) + 8*CA2*Cd(474)

	Cd(9768) = 8*SA2*Cd(472) + Cd(9767)*MStop2(2)

	Cd(9769) = -(Cd(424)*Cd(485)) + Cd(426)*MStop2(4)

	Cd(968) = MStop2(1) + 4*MTy2*UStop2(1,2)

	Cd(9770) = (2*Cd(9769) + 4*MTy*Cd(13)*MStop2(3))/MStop2(3)

	Cd(9778) = MStop2(1) - Cd(1281)*UStop2(1,2)

	Cd(9779) = Cd(9778) - MStop2(1)*UStop2(1,1)

	Cd(969) = 6 + 4*UStop2(1,2)**2

	Cd(9780) = Cd(9779)*UStop2(1,1) + 2*MStop2(3)*UStop2(1,2)**2

	Cd(9783) = MTy*Cd(9780) + AtC*MStop2(4)*UCStop(1,3)*UStop2(1,2)

	Cd(970) = 3 - UStop2(1,2)**2

	Cd(9792) = 2*MTy*SA*UStop(1,1) - Cd(369)*UStop(1,2)

	Cd(9793) = Cd(368)*UStop(1,1) - 2*MTy*SA*UStop(1,2)

	Cd(9794) = 2*MTy*UStop(1,1) + At*UStop(1,2)

	Cd(9795) = SA*Cd(9794) - CA*MUEC*UStop(1,2)

	Cd(9798) = UUStopC(3,3) + UUStopC(3,4)

	Cd(971) = 3*MStop2(1) + 2*MTy2*Cd(970)*UStop2(1,2)

	Cd(9802) = Cd(9793)*UUStopC(3,2) - Cd(9795)*UUStopC(3,4)

	Cd(9803) = Cd(9793)*Cd(9802) - SA2*MStop2(3)*UCStopC(3,1)

	Cd(9807) = 2*MTy*SA*UStop(2,1) - Cd(369)*UStop(2,2)

	Cd(9808) = Cd(368)*UStop(2,1) - 2*MTy*SA*UStop(2,2)

	Cd(9809) = 2*MTy*UStop(2,1) + At*UStop(2,2)

	Cd(9810) = SA*Cd(9809) - CA*MUEC*UStop(2,2)

	Cd(9814) = Cd(9808)*UUStopC(3,2) - Cd(9810)*UUStopC(3,3)

	Cd(9815) = Cd(9808)*Cd(9814) + SA2*MStop2(3)*UCStop(3,1)

	Cd(973) = 2*MStop2(1)**2 - MStop2(3)*MStop2(4)

	Cd(974) = 16*MStop2(1)**2 - 4*Cd(973)*UStop2(1,2)

	Cd(9839) = Cd(775) - 5*MStop2(3)*UStop2(1,1)

	Cd(975) = 2*MStop2(1)**2 - Cd(523)*UStop2(1,2)

	Cd(9840) = Cd(9839)*UStop2(1,1) + Cd(776)*UStop2(1,2)

	Cd(9841) = Cd(780) - Cd(814)*UStop2(1,1)

	Cd(9842) = Cd(9841)*UStop2(1,1)**2 - 2*MStop2(4)*UStop2(1,2)**2

	Cd(976) = 3*MStop2(2)*MStop2(3) - Cd(975)*UStop2(1,2)

	Cd(9852) = 1 + 10*UStop2(1,2)

	Cd(9856) = 2 + 8*UStop2(1,2)

	Cd(977) = 3 + UStop2(1,2)**3

	Cd(9862) = 4*MStop2(1)**2 + 2*Cd(1009)*MStop2(3)

	Cd(9863) = 4*MStop2(1)*MStop2(2) + Cd(9862)*UStop2(1,2)

	Cd(9864) = Cd(9863) - 4*Cd(852)*UStop2(1,1)

	Cd(9866) = 2*MStop2(1)*MStop2(2) + Cd(861)*UStop2(1,2)

	Cd(978) = 4*Cd(977)*MStop2(1) - 4*Cd(970)*MStop2(2)*UStop2(1,2)

	Cd(9872) = MTy2 - AtC*MTy*UCStop(1,3)

	Cd(9875) = 2*Cd(9872) + MStop2(2)

	Cd(979) = Cd(974) - 4*Cd(898)*UStop2(1,1)

	Cd(9878) = Cd(424)*Cd(577) - Cd(426)*MStop2(4)

        Cd(9879) = ((2*Cd(9878) - 4*MTy*Cd(13)*MStop2(3))*
     &     UStop2(1,1)*UStop2(1,2))/MStop2(3)

	Cd(9880) = 2*MUE2*Cd(577) - Cd(427)*MStop2(4)

	Cd(9883) = Cr(28)*MStop2(2) - 2*Cd(286)*UStop2(1,1)*UStop2(1,2)

        Cd(9884) = -2*MStop2(2)*
     &   (Cr(28)*MStop2(2) - 2*Cd(286)*UStop2(1,1)*UStop2(1,2))

	Cd(9885) = Cd(425)*Cd(9884) + 2*At*AtC*Cd(9883)*Cr(28)*MStop2(2)

        Cd(9886) = (Cd(9885) + 
     &     8*MStop2(2)*MStop2(3)**2*UStop2(1,1)*UStop2(1,2))/
     &   (MZ2*SW2*MStop2(3)**2)

	Cd(980) = 4*Cd(976) + Cd(979)*UStop2(1,1)

	Cd(9887) = Cd(9883)*(2*Cd(426) - Cd(424)*Cr(28))

	Cd(9888) = Cd(427)*Cd(9884) + 2*MUE2*Cd(9883)*Cr(28)*MStop2(2)

	Cd(9889) = CA2*Cd(9888) + MZ2*SA2*SW2*Cd(9886)*MStop2(3)**2

        Cd(9890) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(9889) + S2A*Cd(9887)*MStop2(2)))/
     &    (SB2*SW2*MStop2(2)**2*MStop2(3)**2)

        Cd(9892) = (8*Cr(75)**2*
     &     (Cd(471)*MStop2(1) + MStop2(3)**2))/MStop2(3)**2

	Cd(9895) = Cr(28) - 2*UStop2(1,1)*UStop2(1,2)

	Cd(981) = Cd(978)*MStop2(3) - Cd(980)*UStop2(1,1)

	Cd(9900) = Cr(28)*MStop2(2) - MStop2(4)*UStop2(1,1)*UStop2(1,2)

	Cd(9901) = 2*Cd(9900) - Cr(28)*MStop2(3)

        Cd(9903) = Cr(28)*MStop2(4) - 
     &   2*MStop2(1)*UStop2(1,1)*UStop2(1,2)

	Cd(9904) = Cd(1213)*MStop2(2) - Cd(9903)*UStop2(1,1)*UStop2(1,2)

	Cd(982) = 6*MStop2(1)**2 + MStop2(3)*MStop2(4)

        Cd(9911) = Cd(9792)*Cd(9807)*UUStopC(3,1) - 
     &   Cd(9793)*Cd(9810)*UUStopC(3,3)

        Cd(9912) = -2*Cd(9911)*MStop2(2) + 
     &   2*SA2*MStop2(3)**2*UStop2(1,1)

	Cd(9913) = Cd(9793)*UUStopC(3,2) - Cd(9792)*UUStopC(3,4)

        Cd(9914) = -2*Cd(9808)*Cd(9913)*MStop2(2) + 
     &   2*SA2*MStop2(3)**2*UStop2(1,2)

	Cd(983) = MStop2(1)*MStop2(2) + Cd(951)*UStop2(1,2)

	Cd(9915) = Cd(9912) + Cd(9914)

	Cd(9916) = Cd(9915) + 2*Cd(402)*Cd(403)*MStop2(2)

	Cd(9917) = Cd(9793)*UUStopC(1,2) - 2*Cd(9795)*UUStopC(1,3)

	Cd(9918) = Cd(9793)*Cd(9917) + Cd(9792)**2*UUStopC(1,1)

	Cd(9922) = Cd(9792)**2*UUStopC(1,1) + Cd(9793)**2*UUStopC(1,2)

	Cd(9929) = Cd(9807)*UUStopC(3,1) - Cd(9808)*UUStopC(3,4)

        Cd(9930) = -4*Cd(9792)*Cd(9929)*MStop2(1) + 
     &   4*SA2*MStop2(3)**2*UStop2(1,2)

	Cd(9931) = 8*Cd(9808)*UUStopC(3,2) - 8*Cd(9810)*UUStopC(3,3)

	Cd(9934) = 2*SA2*MStop2(2)*UStop(2,2) + Cd(9808)**2*UStopC(2,2)

        Cd(986) = 4*(3*MStop2(3) + 
     &     UStop2(1,2)*(-2*MStop2(1) + Cd(607)*UStop2(1,2)))

        Cd(9958) = UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &   3*UStop2(1,2)**2

        Cd(9961) = Cd(9948)*Cd(9960) + 2*M_3C*UCStop(1,3)**2 + 
     &   At*Cd(9958)*UStop2(1,1)

        Cd(9964) = Cd(9963)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(9948)*Cd(9962) + At*UStop2(1,2)**2)

	Cd(9972) = 1 - Cd(9955)

        Cd(9996) = Cd(9995)*UCStopC(1,3)**2 + 
     &   UStop2(1,1)*(Cd(9948)*Cd(9994) + AtC*UStop2(1,2)**2)

        Cd(10017) = Cd(10015)*Cd(10016)*UCStop(1,3) - 
     &   3*M_3C*AtC*UCStop(1,3)**3 + 
     &   M_3*At*Cd(10016)*UCStopC(1,3)

        Cd(10018) = UStop2(1,1)**2 + 
     &   (Cd(9948) - UStop2(1,1))*UStop2(1,2)

        Cd(10020) = UStop2(1,1)*(UStop2(1,1) - UStop2(1,2)) + 
     &   UStop2(1,2)**2

	Cd(10025) = Cd(290) - 2*UStop2(1,1)

	Cd(10028) = 2*Cd(9948) + Cd(10721)

        Cd(10029) = -(Cd(10028)*UStop2(1,1)) + 
     &   2*UStop2(1,1)**2 - Cd(10027)*UStop2(1,2)

        Cd(10036) = M_3C*AtC*Cd(10016)*UCStop(1,3) + 
     &   Cd(10015)*Cd(10016)*UCStopC(1,3) - 
     &   3*M_3*At*UCStopC(1,3)**3

	Cd(10049) = Cd(10382) + MUEC*Cd(9948)*UStop2(1,1)**2

	Cd(10050) = Cd(10388) + MUE*Cd(9948)*UStop2(1,1)**2

	Cd(10064) = Cd(9948)*Cd(10063) + Cd(10696)

	Cd(10077) = Cd(11622) - 2*UStop2(1,2)

	Cd(10087) = -Cd(10759) + M_3*Cd(9948)*Cr(75)

	Cd(10088) = -Cd(10454) + M_3*Cd(9948)*UStop2(1,1)

        Cd(10091) = (Cd(9948) - UStop2(1,1))*UStop2(1,1) - 
     &   UStop2(1,2)**2

	Cd(10095) = -Cd(10764) + M_3C*Cd(9948)*Cr(75)

	Cd(10096) = -Cd(10321) + M_3C*Cd(9948)*UStop2(1,1)

        Cd(10112) = 2*Cd(9948)*UStop2(1,1) + UStop2(1,1)**2 + 
     &   UStop2(1,2)**2

	Cd(10114) = Cd(425) + 2*At*AtC*Cd(9981)*UStop2(1,1)

	Cd(98) = Cd(2930) + 4*MStop2(2)

        Cd(10138) = AtC*Cd(10136)*UCStop(1,3) + 
     &   At*UCStopC(1,3)*(Cd(10016)*Cd(10137) - 
     &      3*M_3*At*UCStopC(1,3)**2)

        Cd(10143) = Cd(10142)*UCStopC(1,3) + 
     &   UCStop(1,3)*(2*M_3C*UStop2(1,1) + AtC*UStop2(1,2))

        Cd(10153) = AtC**2*Cd(10032) + 
     &   At*(At*Cd(10039)*UCStopC(1,3) + 
     &      2*AtC*Cr(5)*UStop2(1,1)*UStop2(1,2))

        Cd(10163) = (Cd(10162) - 4*UStop2(1,1))*UStop2(1,1) + 
     &   Cd(2361)*UStop2(1,2)

        Cd(10172) = M_3*Cd(9948)*Cd(9983) - 
     &   2*M_3C*UCStop(1,3)**2 + M_3*Cr(75)**2*UStop2(1,1)

        Cd(10182) = 4*Cd(9948)*Cd(10181) + 
     &   Cd(10180)*UCStopC(1,3)**2 + 4*M_3C*Cr(75)**2*UStop2(1,1)

	Cd(10191) = Cd(425) + 2*At*AtC*UStop2(1,1)*UStop2(1,2)

        Cd(10201) = 4*M_3*MTy2 + 
     &   At*(-2*MGlpTmSt2(1) + M_3*At*UCStopC(1,3)**2)

	end
