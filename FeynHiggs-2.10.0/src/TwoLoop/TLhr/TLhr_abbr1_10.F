* TLhr_abbr1_10.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_10
	implicit none

#include "TLhr.h"


	Cd(12436) = -2*Cd(12435) + At*Cd(9969)*Cr(28)

	Cd(12438) = At - 4*M_3*Cd(11881)

	Cd(12442) = -2*Cd(12441) + AtC*Cd(9969)*Cr(28)

	Cd(12446) = Cd(12371)*Cd(12445) - 2*Cd(13)*Cd(12370)*MStop2(1)

	Cd(12448) = AtC*Cd(10018)*Cr(5) - Cd(12447)*UCStopC(1,3)

	Cd(12457) = 2*Cd(12456) + Cd(10025)*MStop2(3)

	Cd(12458) = 2*AtC*Cd(10018)*Cr(5) + Cd(12457)*UCStopC(1,3)

	Cd(12463) = MUE2*Cd(10522) - Cd(427)*MStop2(2)

        Cd(12465) = 2*Cd(47)*Cd(164)*Cr(28)*MStop2(1) - 
     &   Cd(12371)*MStop2(3)*UStop2(1,1)**2

        Cd(12466) = 2*Cd(47)*Cd(164)*Cr(28)*MStop2(1) - 
     &   Cd(12371)*MStop2(3)*UStop2(1,2)**2

	Cd(12467) = Cd(9948)*Cd(12465) + Cd(12466)*UStop2(1,1)

        Cd(12468) = 4*Cd(47)*Cd(164)*Cd(9972)*MStop2(1) + 
     &   Cd(12371)*MStop2(3)

        Cd(12471) = 2*Cd(47)*Cd(12470)*Cr(28)*MStop2(1) - 
     &   Cd(12371)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(12472) = MGlmT2 + Cd(9972)*MStop2(1)

	Cd(12473) = 4*Cd(47)*Cd(12472)*MStop2(1) + Cd(12371)*MStop2(3)

	Cd(12481) = 4*MGl2*MTy2*Cd(179) + 2*Cd(4190)*MGlpTmSt2(2)

	Cd(1241) = Cd(424)*Cd(1229) - 4*Cd(426)*Cd(1236)

	Cd(12483) = 2*Cd(251)*Cd(2681) + 4*MGl2*Cd(179)*MStop2(2)

        Cd(1242) = (Cd(1241) + 
     &     4*MTy*Cd(13)*Cd(1238)*MStop2(3)*MStop2(4))/
     &   (MStop2(1)*MStop2(3)**2)

	Cd(12527) = At*Cd(9970) + 2*AtC*Cd(9969)*UCStop(1,3)**2

	Cd(12528) = 4*Cd(138)*Cd(10169) - Cd(12527)*MStop2(1)

	Cd(12542) = AtC*Cd(9970) + 2*At*UCStopC(1,3)**2

	Cd(12543) = 4*Cd(138)*Cd(10176) - Cd(12542)*MStop2(1)

	Cd(12549) = At*Cd(10070) - 2*Cd(11911)

	Cd(12551) = AtC*Cd(10070) - 2*Cd(11917)

	Cd(1248) = Cd(1246) - Cd(256)*UStop2(1,1)

	Cd(12556) = Cd(11907)*UStop2(1,1) + Cd(9948)*UStop2(1,2)

	Cd(12557) = At*Cd(10020)*Cr(5) + AtC*Cd(10021)*UCStop(1,3)

	Cd(12563) = AtC*Cd(10020)*Cr(5) + At*Cd(10043)*UCStopC(1,3)

	Cd(12572) = 2*MUE2*Cd(577) - Cd(427)*Cd(1272)

	Cd(12583) = Cd(152)*Cd(4095) + Cd(10114)*MGlpTmSt2(1)*MStop2(1)

	Cd(12584) = Cd(179)*Cd(4190) + Cd(10114)*MGlpTmSt2(2)*MStop2(2)

	Cd(12590) = 2*Cd(164)*Cd(10112) + MGlpTmSt2(1)

	Cd(12591) = 2*Cd(251)*Cd(10112) + MGlpTmSt2(2)

	Cd(12594) = Cd(164) + MGlpTmSt2(1)

	Cd(12596) = Cd(251) + MGlpTmSt2(2)

	Cd(12600) = 8*M_3*Cd(164) - Cd(12016)*MGlpTmSt2(1)

        Cd(12601) = Cd(12600)*UCStopC(1,3)**2 + 
     &   8*M_3C*MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2)

	Cd(12602) = 8*M_3*Cd(251) - Cd(12016)*MGlpTmSt2(2)

        Cd(12603) = Cd(12602)*UCStopC(1,3)**2 + 
     &   8*M_3C*MGlpTmSt2(2)*UStop2(1,1)*UStop2(1,2)

	Cd(1254) = Cd(787)*MStop2(1) - 3*MStop2(2)

        Cd(12620) = Cd(427)*Cd(12619) - 
     &   2*MUE2*Cd(12618)*UStop2(1,1)*UStop2(1,2)

	Cd(1255) = Cd(1254) + MStop2(4)*UStop2(1,1)

	Cd(12627) = 2*MTy*Cd(12625) - Cd(12626)*MGlpTmSt4(2)

        Cd(12628) = Cd(12141)*MGlpTmSt4(1) + 
     &   2*MTy*MUEC*Cd(164)*MGlpTmSt2(1)*UUStop(1,3)

	Cd(12630) = 2*MTy2*Cd(12629) + MTy*Cd(12626)*MGlpTmSt2(2)

	Cd(12631) = 2*MTy2*Cd(110) + AtC*MTy*MGlpTmSt2(1)*UCStop(1,3)

	Cd(1256) = Cd(1255)*UStop2(1,1) - Cr(85)*MStop2(3)*UStop2(1,2)

        Cd(12634) = MUE*Cd(12633)*MGlpTmSt2(1) + 
     &   MUEC*Cd(12631)*UUStop(1,2)

	Cd(12639) = 2*MTy2*Cd(12629) - MTy*Cd(12169)*MGlpTmSt2(2)

	Cd(12641) = MUE*Cd(12640) + MUEC*Cd(10991)*MGlpTmSt2(1)

        Cd(12646) = MUEC*Cd(9794)*UCStopC(1,3)*UStop2(1,1) + 
     &   Cd(12149)*UUStopC(1,2)

        Cd(12647) = Cd(12646)*UStopC(1,1) + 
     &   MUE*Cd(10985)*UCStop(1,3)*UStop2(1,2)*UStopC(1,2)

        Cd(12655) = Cd(12654)*UCStopC(1,3)**2 + 
     &   8*M_3C*MGlpTmSt2(1)*UStop2(1,1)*UStop2(1,2)

        Cd(12657) = Cd(12656)*UCStopC(1,3)**2 + 
     &   8*M_3C*MGlpTmSt2(2)*UStop2(1,1)*UStop2(1,2)

        Cd(12722) = 3/(64.D0*Pi)*
     &   (Alfa2L*(4*A0delStop(2)*Cd(12707)*Cd(12708)*MStop2(1) + 
     &        Cd(12721)*MStop2(3)**2))/
     &    (SW2*TB2*MStop2(1)*MStop2(3)**2)

        Cd(12725) = 3/(16.D0*Pi)*
     &   (Alfa2L*Cd(12723)*
     &      (Cd(17)*MStop2(1)*MStop2(2)*MStop2(3) + 
     &        Cd(12724)*MStop2(4)))/
     &    (SW2*TB2*MStop2(1)*MStop2(2)*MStop2(3)**2)

        Cd(12727) = 3/(8.D0*Pi)*
     &   (Alfa2L*Cd(12723)*(2*MStop2(3) - Cd(12726)*MStop2(4)))/
     &    (SW2*TB2*MStop2(3)**2)

        Cd(12729) = 3/(16.D0*Pi)*
     &   (Alfa2L*Cd(12728)*
     &      (Cd(17)*MStop2(1)*MStop2(2)*MStop2(3) + 
     &        Cd(12724)*MStop2(4)))/
     &    (SW2*TB2*MStop2(1)*MStop2(2)*MStop2(3)**2)

        Cd(12730) = 3/(8.D0*Pi)*
     &   (Alfa2L*Cd(12728)*(2*MStop2(3) - Cd(12726)*MStop2(4)))/
     &    (SW2*TB2*MStop2(3)**2)

        Cd(12737) = Cd(1272)*Cd(12700)**2*UCStop(1,3)**2 + 
     &   2*Cd(12736)*MStop2(2)*UStop2(1,1)

	Cd(12740) = Cd(12737) + 2*Cd(12739)

        Cd(12741) = (2*Cd(12740) - 
     &     2*Cd(12701)**2*MStop2(3)*UCStopC(1,3)**2)/MStop2(3)

        Cd(12745) = (MStop2(2)*
     &     (-(Cd(12735)*Cr(34)*MStop2(2)) + 
     &       2*Cd(12741)*MStop2(3)**2))/MStop2(3)**2

        Cd(12746) = 3/(64.D0*Pi)*
     &   (Alfa2L*(-8*Cd(12707)*Cd(12708)*Cr(13)*MStop2(1)*
     &         MStop2(2)**2 + Cd(12745)*MStop2(3)**2))/
     &    (SW2*TB2*MStop2(2)**2*MStop2(3)**2)

	Cd(12748) = -(Cd(12747)*MStop2(1)) + MStop2(3)**2

	Cd(12751) = 2*Cd(12750)*Cr(34) - Cd(12720)*MStop2(3)

	Cd(1269) = 6*MStop2(2)**2 - Cd(874)*MStop2(3)

	Cd(12759) = Cd(12757)*UCStop(1,3)**2 + Cd(12758)*UCStopC(1,3)**2

        Cd(12760) = Cd(12759)*UStop2(1,2) - 
     &   Cd(12700)**2*UCStop(1,3)*UUStop(3,4)*UUStopC(3,2)

	Cd(12766) = 4*(MStop2(2) - Cd(535)*UStop2(1,1)*UStop2(1,2))

	Cd(12786) = Cd(12780) + 4*Cd(12700)*Cd(12782)*UUStopC(2,3)

        Cd(12788) = Cd(12757)*MStop2(4)*UCStop(1,3)**2 + 
     &   Cd(12787)*UCStopC(1,3)**2

        Cd(12789) = Cd(12788)*UStop2(1,2) - 
     &   Cd(12700)**2*MStop2(4)*UCStop(1,3)*UUStop(3,4)*
     &    UUStopC(3,2)

	Cd(12795) = -8*Cd(12794)*MStop2(2) + 8*MStop2(3)**2

	Cd(12798) = Cd(12797)*Cr(13) + Cd(12720)*MStop2(3)

        Cd(12800) = -(Cd(12757)*UCStop(1,3)**2) + 
     &   4*Cd(12700)*Cd(12701)*UStop2(1,1)*UStop2(1,2)

        Cd(12801) = Cd(12800)*UStop(1,1) - 
     &   Cd(12758)*UCStopC(1,3)*UStop(1,2)*UStop2(1,1)

        Cd(12805) = -(Cd(12701)*Cd(12803)*UCStop(3,1)) + 
     &   Cd(12801)*UStop2(1,2)

        Cd(12811) = 2*Cd(286)*UCStopC(1,3)*UCStopC(3,1) - 
     &   MStop2(3)*UCStopC(3,3)*UStop2(1,2)

        Cd(12814) = 2*Cd(12700)*Cd(12813) + 
     &   Cd(12701)*Cd(12811)*UCStop(3,1)*UCStopC(1,3)

        Cd(12815) = Cd(12814) - 
     &   4*Cd(12700)*MStop2(3)*UStop2(1,1)**2*UStop2(1,2)**2

	Cd(125) = 4*Cd(124) - Cd(120)*UStopC(1,1)

        Cd(12826) = Cd(12825) + 
     &   2*Cd(12700)*Cr(28)*MStop2(1)*UStop2(1,1)*UStop2(1,2)

	Cd(1277) = 10*MStop2(2)**2 - 3*Cd(535)*MStop2(3)

	Cd(1278) = MStop2(2)**2 + Cd(1277)*UStop2(1,2)

	Cd(12850) = -2*Cd(12849) + 2*Cd(589)*Cd(12701)*UStop2(1,1)

	Cd(12851) = Cd(12848) + Cd(12850)

        Cd(12853) = 1/2.D0*(A0delStop(2)*Cd(12850)) - 
     &   Cd(12846)*Cr(34)*MStop2(2)*MStop2(3)

        Cd(12854) = 3/(16.D0*Pi)*
     &   (Alfa2L*(Cd(12853)*MStop2(1) + Cd(12852)*MStop2(2)))/
     &    (SW2*TB2*MStop2(1)*MStop2(2)*MStop2(3))

	Cd(12856) = Cd(12850)*Cr(34) - 2*Cd(12846)*MStop2(3)

	Cd(12857) = Cd(12855) + Cd(12856)

	Cd(1279) = 3*MStop2(1)*MStop2(3) + Cd(1278)*UStop2(1,2)

	Cd(12860) = Cd(12858) + Cd(12859)

        Cd(12861) = -(3/(8.D0*Pi)*
     &     (Alfa2L*(Cd(12860) - 16*MTy2*MStop2(3)))/
     &      (MTy*SW2*TB2*MStop2(3)))

        Cd(12863) = 1/2.D0*(A0delStop(2)*Cd(12858)) - 
     &   Cd(12720)*Cr(34)*MStop2(2)*MStop2(3)

        Cd(12864) = 1/2.D0*(A0delStop(1)*Cd(12859)) - 
     &   Cd(12720)*Cr(13)*MStop2(1)*MStop2(3)

	Cd(12865) = Cd(12864) + 4*Cd(12862)*MStop2(1)*MStop2(3)

	Cd(12866) = Cd(12863)*MStop2(1) + Cd(12865)*MStop2(2)

	Cd(12868) = Cd(12858)*Cr(34) - 2*Cd(12720)*MStop2(3)

	Cd(12869) = Cd(12859)*Cr(13) - 2*Cd(12720)*MStop2(3)

        Cd(12871) = Cd(12701)*Cr(75)*UStop2(1,1)*UStop2(1,2) + 
     &   Cd(9798)*Cd(12700)*UCStop(1,3)*UUStop(3,1)

        Cd(12874) = -(Cd(12700)*Cd(12873)*UStop2(1,2)) + 
     &   Cd(9798)*Cd(12701)*UCStopC(1,3)*UUStop(3,2)

        Cd(12875) = Cd(12874) - 
     &   Cd(12701)*Cr(91)*UCStopC(1,3)*UUStopC(3,1)

        Cd(12877) = Cd(286)*Cd(12701)*Cr(75)*UStop2(1,1)*
     &    UStop2(1,2) + Cd(9798)*Cd(12700)*MStop2(4)*UCStop(1,3)*
     &    UUStop(3,1)

        Cd(12881) = Cd(12700)*Cd(12880)*UStop2(1,2) + 
     &   Cd(9798)*Cd(12701)*MStop2(4)*UCStopC(1,3)*UUStop(3,2)

	Cd(12885) = 12*MStop2(2)**2 - 5*Cd(535)*MStop2(3)

	Cd(12891) = 2*Cd(523)*UUStop(3,3) + 4*MStop2(1)**2*UUStop(3,4)

        Cd(12893) = 2*Cd(12892)*MStop2(1)*MStop2(2) + 
     &   Cd(12891)*UStop2(1,2)*UStopC(2,1)

	Cd(12894) = 4*MStop2(2)**2*UUStop(3,3) + 2*Cd(523)*UUStop(3,4)

        Cd(12896) = 2*Cd(12895)*MStop2(1)*MStop2(2) + 
     &   Cd(12894)*UCStop(1,3)

        Cd(12897) = Cd(12701)*Cd(12893)*UStop(2,2) + 
     &   Cd(12700)*Cd(12896)*UStop2(1,2)

        Cd(12898) = 2*Cd(523)*UStop(1,1)*UStop2(1,2) - 
     &   4*Cr(91)*MStop2(1)*MStop2(2)*UStopC(2,2)

        Cd(12899) = 4*Cd(11080)*MStop2(1)*MStop2(2)*
     &    UStop2(1,2) + Cd(12894)*UStop2(1,1)*UStopC(2,2)

        Cd(12901) = Cd(12900)*UStopC(2,1) - 
     &   UStop(2,1)*UStop2(1,1)*UUStop(1,2)

        Cd(12904) = Cd(12899)*UStop2(1,1) + 
     &   2*Cd(523)*Cd(12901)*UUStopC(3,2)

	Cd(12907) = 2*(Cd(10916)*MStop2(2) + MStop2(1)**2*UStop2(1,2))

	Cd(12911) = Cd(12909) + Cd(12910)

        Cd(12912) = 2*Cr(28)*MStop2(1)*MStop2(2)*UCStopC(1,3) - 
     &   Cd(12911)*UUStopC(3,1)

	Cd(12913) = 4*MStop2(1)**2*UUStop(3,3) + 2*Cd(523)*UUStop(3,4)

	Cd(12914) = 2*Cd(523)*UUStop(3,3) + 4*MStop2(2)**2*UUStop(3,4)

	Cd(12915) = Cd(12913)*UUStopC(3,3) + Cd(12914)*UUStopC(3,4)

        Cd(1285) = 3*MStop2(2)*MStop2(3) - 
     &   Cd(535)*MStop2(1)*UStop2(1,1)**3

        Cd(12917) = Cd(12916)*UStop2(1,2) - 
     &   Cd(11106)*UStop2(1,1)**2*UUStopC(3,4)

        Cd(12918) = 2*Cd(12917)*MStop2(1)*MStop2(2) + 
     &   Cd(12915)*UStop2(1,1)*UStop2(1,2)

        Cd(12919) = Cd(12700)*Cd(12918) - 
     &   Cd(12701)*Cd(12912)*UCStopC(1,3)

	Cd(12920) = MUEC*TB2*Cd(12907) + 2*At*TB*Cd(12908)

        Cd(12921) = Cd(12919) + 
     &   2*Cd(11328)*Cd(12701)*Cr(91)*MStop2(4)*UCStopC(1,3)*
     &    UUStopC(3,1)

        Cd(12930) = UStop2(1,1)*UStop2(1,2)**3 - 
     &   Cd(12901)*UStopC(1,1)*UUStopC(3,2)

        Cd(12932) = Cd(12930)*MStop2(3) - 
     &   Cd(12931)*UStop2(1,1)*UUStopC(3,3)

        Cd(1287) = 6*MStop2(3)*MStop2(4) + 
     &   4*Cd(286)*MStop2(2)*UStop2(1,1)**2

        Cd(12942) = Cd(12701)*Cd(12932) + 
     &   Cd(12929)*UStop2(1,2)*UUStopC(3,2)

        Cd(12945) = 2*Cd(12701)*Cr(75)*UStop2(1,1)*UStop2(1,2) + 
     &   Cd(9798)*Cd(12700)*UCStop(1,3)*UUStop(3,1)

        Cd(12946) = Cd(12700)*Cd(12878)*UStop2(1,2) + 
     &   Cd(9798)*Cd(12701)*UCStopC(1,3)*UUStop(3,2)

	Cd(1288) = Cd(1264)*UStop2(1,1) - Cd(307)*MStop2(2)*UStop2(1,2)

        Cd(12947) = Cd(12946) - 
     &   Cd(12701)*Cr(91)*UCStopC(1,3)*UUStopC(3,1)

        Cd(12949) = Cd(12700)*Cd(12879)*UStop2(1,2) + 
     &   Cd(9798)*Cd(12701)*MStop2(4)*UCStopC(1,3)*UUStop(3,2)

        Cd(12950) = (2*Cd(12949)*Cr(9) - 
     &     2*Cd(12700)*Cr(9)*Cr(91)*MStop2(2)*UStop2(1,2)*
     &      UUStopC(3,3))/MStop2(2)

        Cd(12951) = 2*Cr(9)*
     &   (Cd(12700)*Cd(12877) + Cd(12882)*Cr(91)*MStop2(4))

        Cd(12952) = (-Cd(12951) + Cd(12701)*Cd(12950)*MStop2(2))/
     &   MStop2(2)

	Cd(12954) = Cd(152) + 2*MStop2(1)

        Cd(12959) = -(Cd(152)*Cd(12958)) + 
     &   Cd(12954)*MGlpTmSt2(1)*MStop2(3)

        Cd(12964) = MTy*Cd(12954)*MStop2(3) - 
     &   Cd(152)*Cd(12701)*Cd(12961)*UStopC(1,1)

        Cd(12971) = Cd(152)*Cd(12970) + 
     &   MTy*Cd(12954)*MStop2(3)*UCStop(1,3)

        Cd(12975) = 4*Cd(12700)*MStop2(4)*UUStop(1,1) - 
     &   4*Cd(778)*Cd(12701)*UUStop(1,2)

	Cd(12977) = Cd(12975)*UStop2(1,1) + 4*Cd(12974)*UStop2(1,2)

        Cd(12978) = Cd(12977)*UStopC(1,2) + 
     &   4*Cd(12976)*MStop2(1)*UStop2(1,1)*UStopC(2,2)

	Cd(12980) = AtC*Cd(778) + MUE*TB*Cd(814)

        Cd(12981) = Cd(12980)*UUStop(1,1) - 
     &   Cd(12701)*MStop2(2)*UUStop(1,2)

        Cd(12985) = Cd(12979)*MStop2(2)*UStop2(1,1) - 
     &   Cd(12981)*UStop2(1,2)

        Cd(12986) = Cd(12983)*UStopC(1,2)*UUStopC(2,1) - 
     &   Cd(12984)*UCStop(3,1)*UUStopC(3,1)

	Cd(1292) = Cd(1290) - Cd(875)*UStop2(1,1)

	Cd(12987) = -4*Cd(12986)*MStop2(1) + 4*Cd(12985)*UUStopC(1,1)

        Cd(12988) = 1/(256.D0*Pi**4)*
     &   (Cd(12987) - 4*Cd(12700)*Cr(91)*MStop2(1)*UStop2(1,1)*
     &       UUStopC(3,3))/(MStop2(1)*MStop2(3))

	Cd(12989) = 4*MStop2(3)**2 - Cd(12700)*Cd(12978)*UStopC(1,2)

	Cd(13000) = Cd(179) + 2*MStop2(2)

        Cd(13005) = -(Cd(179)*Cd(12701)*Cd(13004)) + 
     &   MTy*Cd(13000)*MStop2(3)

        Cd(13010) = MTy*Cd(13000)*MStop2(3) + 
     &   Cd(179)*Cd(12700)*Cd(13009)*UCStop(1,3)

	Cd(13012) = MStop2(2) + Cd(179)*UStop2(1,2)

	Cd(1295) = 5*MStop2(2)**2 - Cd(874)*MStop2(3)

        Cd(13022) = Cd(13021)*UStopC(1,1) - 
     &   Cd(12983)*MGlpTmSt2(2)*UStopC(1,2)

        Cd(13023) = Cd(179)*
     &   (Cd(13022)*UUStopC(2,1) + 
     &     Cd(12984)*MGlpTmSt2(2)*UCStop(3,1)*UUStopC(3,1))

        Cd(13024) = Cd(13023) - 
     &   Cd(179)*Cd(12700)*Cd(13019)*MGlpTmSt2(2)*UStop2(1,1)*
     &    UStopC(1,1)

	Cd(1296) = 6*MStop2(2)*MStop2(3) + 4*Cd(1295)*UStop2(1,2)**2

        Cd(13029) = Cd(12700)*MStop2(4)*UUStop(1,1) - 
     &   Cd(814)*Cd(12701)*UUStop(1,2)

	Cd(13031) = Cd(13029)*UStop2(1,1) + Cd(12974)*UStop2(1,2)

        Cd(13032) = 2*Cd(13031)*UStopC(1,2) + 
     &   Cd(13030)*MStop2(1)*UStop2(1,1)*UStopC(2,2)

        Cd(13033) = Cd(778)*Cd(12700)*UUStop(1,1) - 
     &   Cd(12701)*MStop2(2)*UUStop(1,2)

	Cd(13035) = Cd(13034)*UStopC(1,1) - Cd(12983)*UStopC(1,2)

        Cd(13036) = Cd(13035)*UUStopC(2,1) + 
     &   Cd(12984)*UCStop(3,1)*UUStopC(3,1)

        Cd(13037) = -2*Cd(13036)*MStop2(1) + 
     &   2*Cd(13033)*UStop2(1,2)*UUStopC(1,1)

        Cd(13038) = (Cd(13037) - 
     &     2*Cd(12979)*MStop2(2)*UStop2(1,1)*UUStopC(1,1))/
     &   MStop2(2)

        Cd(13039) = 1/(128.D0*Pi**4)*
     &   (Cd(13038)*MStop2(2) + 
     &      2*Cd(12700)*Cr(91)*MStop2(1)*UStop2(1,1)*UUStopC(3,3))/
     &    MStop2(2)

        Cd(13040) = 1/(128.D0*Pi**4)*
     &   (128*Pi**4*(Cd(12701)*Cd(13039)*MStop2(2)) - 
     &      Cd(12700)*Cd(13032)*UStopC(1,2))/MStop2(2)

        Cd(13041) = 1/(64.D0*Pi**4)*
     &   (-(64*Pi**4*(Cd(13040)*MStop2(2))) + 
     &      Cd(13028)*MStop2(3)**2)/(MStop2(2)*MStop2(3))

	Cd(13045) = -(M_3C*Cd(12993)) + M_3*Cd(12995)

        Cd(13046) = 2*MTy*Cd(13045)*MGlpTmSt2(2) + 
     &   Cd(12998)*MGlpTmSt4(2)

        Cd(13048) = 1/(64.D0*Pi**4)*
     &   (-2*MTy2*Cd(12958) + Cd(13047)*MGlpTmSt2(1)*MStop2(3))/
     &    MTy2

        Cd(13049) = 1/2.D0*
     &   (MGl*(Cd(118) - MTy2*Cd(10025))*MStop2(3)*UCStop(1,3))/MTy

	Cd(13052) = Cd(13049) + MGl*Cd(13051)

        Cd(13053) = 1/2.D0*
     &   (MGl*(Cd(118) - MTy2*Cd(10025))*MStop2(3)*UCStopC(1,3))/
     &    MTy

	Cd(13055) = Cd(13053) + MGl*Cd(13054)

	Cd(13056) = M_3C*Cd(13052) + M_3*Cd(13055)

        Cd(13057) = -(1/(16.D0*Pi**4)*
     &     (Cd(13056) - 16*Pi**4*(MGl*Cd(13048)))/
     &      (MGl*TB2*MStop2(3)))

        Cd(13058) = 1/2.D0*
     &   (MGl*(2*MTy2 + Cd(179))*MStop2(3)*UCStopC(1,3))/MTy

	Cd(13060) = Cd(13058) + MGl*Cd(13059)

	Cd(13061) = 2*MTy2 + Cd(179)

        Cd(13063) = 1/(32.D0*Pi**4)*
     &   (MGl*(-2*MTy*Cd(13062) + 
     &        Cd(13061)*MStop2(3)*UCStop(1,3)))/MTy

        Cd(13066) = Cd(13022)*UUStopC(2,1) + 
     &   Cd(12984)*MGlpTmSt2(2)*UCStop(3,1)*UUStopC(3,1)

        Cd(13067) = Cd(13066) - 
     &   Cd(12700)*Cd(13019)*MGlpTmSt2(2)*UStop2(1,1)*UStopC(1,1)

	Cd(13070) = M_3C*Cd(13063) + 1/(16.D0*Pi**4)*(M_3*Cd(13060))

	Cd(13075) = MTy2*Cd(164) + MGlpTmSt2(1)*MStop2(1)

        Cd(13076) = 1/(64.D0*Pi**4)*
     &   (MTy2*Cd(164)*Cd(12958) - 
     &      Cd(13075)*MGlpTmSt2(1)*MStop2(3))/(MTy2*MStop2(1))

        Cd(13080) = 1/(32.D0*Pi**4)*
     &   (MGl*Cd(164)*Cd(12970) + Pi**4*(Cd(13079)*MStop2(1)))/
     &    MStop2(1)

        Cd(13081) = 1/Pi**4*
     &   (MGl*(MTy2*Cd(164) + MGlpTmSt2(1)*MStop2(1))*MStop2(3)*
     &      UCStopC(1,3))/(MTy*MStop2(1))

        Cd(13082) = 1/(32.D0*Pi**4)*
     &   (MGl*Cd(164)*Cd(13054) + Pi**4*(Cd(13081)*MStop2(1)))/
     &    MStop2(1)

	Cd(13083) = M_3C*Cd(13080) + M_3*Cd(13082)

	Cd(13084) = (MGl*Cd(13076) + Cd(13083))/(MGl*TB2*MStop2(3))

        Cd(13087) = 1/Pi**4*
     &   (MGl*(MTy2*Cd(251) + MGlpTmSt2(2)*MStop2(2))*MStop2(3)*
     &      UCStopC(1,3))/(MTy*MStop2(2))

        Cd(13088) = 1/(32.D0*Pi**4)*
     &   (MGl*Cd(251)*Cd(13059) + Pi**4*(Cd(13087)*MStop2(2)))/
     &    MStop2(2)

	Cd(13089) = MTy2*Cd(251) + MGlpTmSt2(2)*MStop2(2)

        Cd(13090) = 1/(32.D0*Pi**4)*
     &   (MGl*(MTy*Cd(251)*Cd(13062) - 
     &        Cd(13089)*MStop2(3)*UCStop(1,3)))/(MTy*MStop2(2))

	Cd(13100) = M_3*Cd(13088) - M_3C*Cd(13090)

        Cd(13103) = 1/(64.D0*Pi**4)*
     &   (-(64*Pi**4*(MTy2*Cd(13040)*MStop2(2))) + 
     &      Cd(13102)*MStop2(3)**2)/(MTy2*MStop2(2)*MStop2(3))

	Cd(1304) = 2*MTy2*Cd(1303) - 3*MStop2(2)*UStop2(1,2)

        Cd(13110) = ((-(Cd(13109)*MStop2(1)) + 
     &       Cd(13108)*MStop2(2))*UCStopC(1,3)*UStopC(1,1))/
     &   MStop2(1)**2

	Cd(13116) = Cd(607) + MStop2(3)

        Cd(13117) = 1/2.D0*(Cd(12700)*Cd(13116)*UUStop(1,1)) + 
     &   Cd(12701)*MStop2(4)*UUStop(1,2)

        Cd(13124) = 1/2.D0*
     &   ((Cd(13123)*MStop2(1) + Cd(13371)*MStop2(2))*UCStopC(1,3)*
     &      UStopC(1,1))/MStop2(1)**2

        Cd(13131) = Cd(12700)**2*Cd(13122)*UCStop(1,3) - 
     &   Cd(12701)**2*Cd(13124)*UStop(1,2)

        Cd(13136) = 2*Cd(13135)*UStopC(1,1) + 
     &   UCStopC(3,1)*UStop2(1,2)*UStopC(2,1)

        Cd(13137) = 2*(Cd(13136)*UCStopC(1,3) + 
     &     Cd(10925)*UStop2(1,2)*UStopC(2,1)*UUStopC(1,1))

        Cd(13139) = Cd(1657)*Cd(12700)*MStop2(3)*UUStop(1,1) + 
     &   Cd(951)*Cd(12701)*UUStop(1,2)

	Cd(13141) = 2*Cd(800)*UCStop(1,3) + Cd(11106)*UUStopC(3,2)

	Cd(13143) = (AtC + MUE*TB)*Cd(951)

        Cd(13144) = Cd(13143)*UUStop(1,1) + 
     &   Cd(607)*Cd(12701)*MStop2(3)*UUStop(1,2)

	Cd(13145) = 3*MStop2(3)**2 + Cd(5706)*UStop2(1,1)*UStop2(1,2)

        Cd(1308) = 4*MStop2(2)*MStop2(4) - 
     &   4*Cd(874)*MStop2(1)*UStop2(1,1)

	Cd(1309) = Cd(1308) + 2*Cd(523)*UStop2(1,2)

        Cd(13163) = Cd(5519)*UStop2(1,1) - 
     &   2*Cd(286)*MStop2(2)*UStop2(1,2)

        Cd(13164) = Cd(13163)*UStop(1,1) + 
     &   2*MStop2(1)*MStop2(3)*UCStop(3,1)*UStop(2,1)

        Cd(13165) = -(Cr(1)*MStop2(2)*UStop(1,1)*
     &      UStop2(1,1)**2) + Cd(13164)*UStop2(1,2)

	Cd(13166) = MStop2(1)*MStop2(3) + Cd(814)*MStop2(2)*UStop2(1,2)

        Cd(13167) = Cd(13166)*UStop(1,2)*UStop2(1,1) + 
     &   Cr(91)*MStop2(1)*MStop2(3)*UStop2(1,2)*UStopC(2,1)

        Cd(13168) = Cd(5519)*UUStop(3,3) + 
     &   MStop2(1)*MStop2(2)*UUStop(3,4)

        Cd(13169) = 2*MStop2(2)*MStop2(3)*UStop(1,2)*
     &    UStop2(1,2) + Cd(13168)*UStopC(2,1)

	Cd(13171) = Cd(13167) - Cd(13169)*UStop2(1,2)

        Cd(13172) = Cd(13170)*MStop2(3)*UCStopC(1,3)*
     &    UStop2(1,1) - Cd(13171)*UStopC(1,1)

	Cd(13175) = MStop2(1)**2*UUStop(3,3) + Cd(5519)*UUStop(3,4)

        Cd(13176) = Cd(12701)*Cd(13175)*UStop(1,2) + 
     &   2*Cd(12700)*MStop2(1)*MStop2(2)*UStop(2,1)*UUStop(1,1)

	Cd(13181) = Cd(13180) + 2*Cd(13176)*UStopC(1,1)

        Cd(13182) = 2*Cd(13173)*MStop2(2)*UStop2(1,2)*
     &    UUStopC(1,2) + Cd(13181)*UUStopC(3,2)

        Cd(13183) = Cd(13182)*UStop2(1,1) + 
     &   Cd(13174)*MStop2(1)*MStop2(3)*UCStopC(3,1)*UStopC(2,2)*
     &    UUStopC(2,2)

        Cd(13186) = 2*Cd(12700)*Cr(1)*MStop2(2)*UUStop(1,1) + 
     &   Cd(874)*Cd(12701)*MStop2(3)*UUStop(1,2)

        Cd(13192) = Cd(13186)*UStop2(1,1)*UUStopC(1,1) + 
     &   Cd(13190)*UStopC(1,1)*UUStopC(2,1)

	Cd(13196) = Cd(13195)*MStop2(1) - Cd(13192)*UStop2(1,2)

        Cd(13197) = 1/(128.D0*Pi**4)*
     &   (Cd(13196) + 2*Cd(12700)*Cd(13184)*MStop2(2)*
     &       UStop2(1,1)**2*UStop2(1,2)**2)/(MStop2(2)*MStop2(3))

        Cd(13198) = Cd(13183) + 
     &   Cd(12700)*Cd(13165)*UCStop(1,3)*UStopC(1,2)

        Cd(13199) = Cd(12700)*Cd(13198) - 
     &   Cd(12701)**2*Cd(13172)*UCStopC(1,3)

        Cd(13200) = 1/(128.D0*Pi**4)*
     &   (-Cd(13199) + 128*Pi**4*
     &       (Cd(12701)*Cd(13197)*MStop2(2)*MStop2(3)))/
     &    (MStop2(2)*MStop2(3))

        Cd(13205) = -(1/2.D0*
     &      (Cr(91)*MStop2(1)*MStop2(3)*UStop2(1,2)*UStopC(2,1)))+
     &     Cd(13204)*UStop(1,2)*UStop2(1,1)

	Cd(1314) = Cd(1312) - Cd(874)*UStop2(1,1)**2

	Cd(13206) = Cd(13205) + Cd(13203)*MStop2(3)*UStop2(1,2)

        Cd(13207) = Cd(13170)*MStop2(3)*UCStopC(1,3)*
     &    UStop2(1,1) + 2*Cd(13206)*UStopC(1,1)

        Cd(13214) = Cd(13208)*UUStopC(1,2) - 
     &   Cd(12701)*Cd(13213)*UUStopC(3,3)

	Cd(13215) = 2*Cd(13214)*UStop2(1,2) + Cd(13210)*UUStopC(3,2)

        Cd(13216) = 1/2.D0*
     &    (Cd(13174)*MStop2(1)*MStop2(3)*UCStopC(3,1)*UStopC(2,2)*
     &      UUStopC(2,2)) + Cd(13215)*UStop2(1,1)

	Cd(13226) = Cd(164) + MGlpTmSt2(1)*UStop2(1,1)

        Cd(13227) = 1/Pi**4*
     &   (MGl*(MTy*Cd(13226) + Cd(13078)*MGlpTmSt2(1))*MStop2(3)*
     &     UCStop(1,3))

        Cd(13228) = 1/32.D0*
     &   (Cd(13227) + 1/Pi**4*(MGl*Cd(12970)*MGlpTmSt2(1)))

        Cd(13229) = MTy*Cd(12594)*MStop2(3) - 
     &   Cd(12701)*Cd(12961)*MGlpTmSt2(1)*UStopC(1,1)

	Cd(1317) = MStop2(3) - Cd(535)*UStop2(1,2)

	Cd(1318) = 3*MStop2(2) - 2*Cd(1317)*UStop2(1,2)

        Cd(13246) = 2*MTy*Cd(800)*Cd(12700)*UStop2(1,2) - 
     &   Cd(13245)*UStopC(1,1)

	Cd(13249) = Cd(13246) - 3*MTy*Cd(12701)*UCStopC(1,3)**2

        Cd(13253) = Cd(302)*Cd(12700)*UUStop(1,1) - 
     &   Cd(12701)*Cd(13252)*UUStop(1,2)

	Cd(13254) = Cd(13253)*UStop2(1,1) - Cd(13251)*UStop2(1,2)

	Cd(1319) = 6*MStop2(1)**2 - 2*Cd(1318)*MStop2(2)

	Cd(13262) = AtC*Cd(13252) + MUE*TB*Cd(13261)

        Cd(13263) = Cd(13262)*UUStop(1,1) - 
     &   Cd(12701)*MGlpTmSt2(2)*UUStop(1,2)

        Cd(13264) = Cd(12979)*MGlpTmSt2(2)*UStop2(1,1) - 
     &   Cd(13263)*UStop2(1,2)

	Cd(13265) = Cd(12986)*MGlpTmSt2(1) - Cd(13264)*UUStopC(1,1)

        Cd(13267) = Cd(13265) - 
     &   Cd(12700)*Cd(13266)*MGlpTmSt2(1)*UStop2(1,1)

        Cd(13268) = 2*Cd(12700)*MGlpTmSt2(1)*UCStop(1,3)**2*
     &    UStop2(1,1) - Cd(13254)*UUStopC(1,2)

        Cd(13269) = Cd(13267) + 
     &   Cd(11106)*Cd(12701)*MGlpTmSt2(1)*UCStopC(1,3)*UUStopC(3,1)

	Cd(13270) = Cd(12700)*Cd(13268) + Cd(12701)*Cd(13269)

        Cd(13271) = 1/(64.D0*Pi**4)*
     &   (Cd(13270) + Cd(13260)*MStop2(3))/MStop2(3)

	Cd(13275) = MTy*Cd(10279)*Cr(5) - Cd(10280)*MGlpTmSt2(1)

	Cd(13289) = 2*M_3*MTy*Cd(12995) + Cd(12998)*MGlpTmSt2(2)

	Cd(13290) = M_3C*MTy*Cd(10603)*Cd(12993) - Cd(10604)*Cd(13289)

	Cd(1324) = MStop2(1)**2 - Cd(256)*MStop2(3)

        Cd(13314) = 1/(512.D0*Pi**4)*
     &   (-(512*Pi**4*(Cd(13313)*MGlpTmStxGlT4(1)**2)) + 
     &      Cd(152)*Cd(2597)*Cd(12836)**2*Cr(45)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(1)**2*MStop2(3))

        Cd(13317) = Cd(13316)*MStop2(3) + 
     &   2*Cd(12700)*Cd(12701)*Cd(13120)*UStop2(1,1)**2*UStop2(1,2)

        Cd(13321) = 1/(128.D0*Pi**4)*
     &   (2*MTy2*Cd(12958) - Cd(13047)*MGlpTmSt2(1)*MStop2(3))/MTy2

        Cd(13322) = 1/4.D0*
     &   (MGl*(Cd(118) - MTy2*Cd(10025))*MStop2(3)*UCStop(1,3))/MTy

	Cd(13323) = MGl*Cd(13051) + 2*Cd(13322)

        Cd(13324) = 1/4.D0*
     &   (MGl*(Cd(118) - MTy2*Cd(10025))*MStop2(3)*UCStopC(1,3))/
     &    MTy

	Cd(13325) = MGl*Cd(13054) + 2*Cd(13324)

	Cd(13326) = M_3C*Cd(13323) + M_3*Cd(13325)

        Cd(13327) = 1/(32.D0*Pi**4)*
     &   (Cd(13326) + 32*Pi**4*(MGl*Cd(13321)))/(MGl*TB2*MStop2(3))

        Cd(13328) = 1/4.D0*
     &   (MGl*(2*MTy2 + Cd(179))*MStop2(3)*UCStopC(1,3))/MTy

	Cd(13329) = MGl*Cd(13059) + 2*Cd(13328)

        Cd(13332) = 1/(64.D0*Pi**4)*
     &   (2*MTy*Cd(13331) - MGl*Cd(13061)*MStop2(3)*UCStop(1,3))/
     &    MTy

	Cd(13333) = Cd(13067) - Cd(12701)*MStop2(3)*UCStopC(1,3)**2

	Cd(13334) = -(Cd(12700)*Cd(13025)) + Cd(12701)*Cd(13333)

        Cd(13335) = 1/(128.D0*Pi**4)*
     &   (2*MTy2*Cd(13334) - Cd(13065)*MGlpTmSt2(2)*MStop2(3))/MTy2

	Cd(13336) = M_3C*Cd(13332) - 1/(32.D0*Pi**4)*(M_3*Cd(13329))

	Cd(13337) = (MGl*Cd(13335) + Cd(13336))/(MGl*TB2*MStop2(3))

        Cd(13340) = 1/(256.D0*Pi**4)*
     &   (MTy2*Cd(13283) - 2*Cd(45)*MStop2(1)*MStop2(3))/
     &    (MTy2*MStop2(1)**2)

        Cd(13341) = 1/(128.D0*Pi**4)*
     &   (-(MTy2*Cd(164)*Cd(12958)) + 
     &      Cd(13075)*MGlpTmSt2(1)*MStop2(3))/(MTy2*MStop2(1))

        Cd(13343) = 1/(64.D0*Pi**4)*
     &   (MGl*Cd(164)*Cd(12970) + 2*Pi**4*(Cd(13342)*MStop2(1)))/
     &    MStop2(1)

        Cd(13344) = 1/(2.D0*Pi**4)*
     &   (MGl*(MTy2*Cd(164) + MGlpTmSt2(1)*MStop2(1))*MStop2(3)*
     &      UCStopC(1,3))/(MTy*MStop2(1))

        Cd(13345) = 1/(64.D0*Pi**4)*
     &   (MGl*Cd(164)*Cd(13054) + 2*Pi**4*(Cd(13344)*MStop2(1)))/
     &    MStop2(1)

	Cd(13346) = M_3C*Cd(13343) + M_3*Cd(13345)

	Cd(13347) = (MGl*Cd(13341) - Cd(13346))/(MGl*TB2*MStop2(3))

        Cd(13350) = 1/(2.D0*Pi**4)*
     &   (MGl*(MTy2*Cd(251) + MGlpTmSt2(2)*MStop2(2))*MStop2(3)*
     &      UCStopC(1,3))/(MTy*MStop2(2))

        Cd(13351) = 1/(64.D0*Pi**4)*
     &   (MGl*Cd(251)*Cd(13059) + 2*Pi**4*(Cd(13350)*MStop2(2)))/
     &    MStop2(2)

        Cd(13352) = 1/(64.D0*Pi**4)*
     &   (MGl*(-(MTy*Cd(251)*Cd(13062)) + 
     &        Cd(13089)*MStop2(3)*UCStop(1,3)))/(MTy*MStop2(2))

	Cd(1329) = Cd(1263) - Cd(1326)

	Cd(13354) = M_3*Cd(13351) + M_3C*Cd(13352)

        Cd(13360) = 1/(128.D0*Pi**4)*
     &   (Cd(13046)*MGlpTmStxGlT4(1)**2 + 
     &      Cd(12836)**2*Cd(13044)*MGlpTmStxGlT4(2)**2)/
     &    (TB2*MGlpTmStxGlT4(1)**2*MGlpTmStxGlT4(2)**2)

        Cd(13362) = 1/(512.D0*Pi**4)*
     &   (512*Pi**4*(Cd(13313)*MGlpTmStxGlT4(1)**2) + 
     &      Cd(164)*Cd(2597)*Cd(12836)**2*Cr(45)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(1)**2*MStop2(3))

	Cd(1330) = Cd(1262) - 2*Cd(1324)

        Cd(13384) = 8*MStop2(2)*UCStop(3,1)*UStop2(1,1)*
     &    UStopC(1,1) + 4*Cd(970)*MStop2(3)*UStopC(2,1)

        Cd(13385) = Cd(13384)*UStop(2,2) + 
     &   4*MStop2(3)*UCStopC(1,3)*UStop2(1,1)**2

        Cd(13387) = Cd(13386)*UCStopC(1,3) + 
     &   Cd(10925)*UStop2(1,2)*UUStopC(3,1)

        Cd(13391) = MStop2(3)**2*UStop2(1,1)**2 + 
     &   2*Cd(535)*MStop2(2)*UStop2(1,2)**2

	Cd(13393) = Cd(13391) + Cd(13392)*MStop2(3)

        Cd(13402) = (3*Cr(28)*MStop2(3) + 
     &     Cd(13401)*UStop2(1,1)*UStop2(1,2))/MStop2(3)

        Cd(13431) = 1/Pi**4*
     &   (2*Cd(13430)*MStop2(2) - 
     &      2*Cd(13428)*UStop2(1,1)*UStop2(1,2))/MStop2(2)**2

        Cd(13433) = 1/(512.D0*Pi**4)*
     &   (Pi**4*(Cd(12700)*Cd(13431)*MStop2(2)**2) + 
     &      2*Cd(13432)*UUStopC(2,1))/MStop2(2)**2

	Cd(1337) = -2*Cd(832) + Cd(1269)

        Cd(13442) = -Cd(13265) + 
     &   Cd(12700)*Cd(13266)*MGlpTmSt2(1)*UStop2(1,1)

        Cd(13443) = Cd(13442) - 
     &   Cd(11106)*Cd(12701)*MGlpTmSt2(1)*UCStopC(1,3)*UUStopC(3,1)

	Cd(13444) = Cd(12700)*Cd(13268) - Cd(12701)*Cd(13443)

        Cd(13445) = -(1/(64.D0*Pi**4)*
     &     (Cd(13444) - Cd(13441)*MStop2(3))/MStop2(3))

        Cd(13446) = Cd(12701)*Cd(13004)*MGlpTmSt2(2) - 
     &   MTy*Cd(12596)*MStop2(3)

        Cd(13448) = -(Cd(13062)*MGlpTmSt2(2)) + 
     &   MTy*Cd(12596)*MStop2(3)*UCStop(1,3)

	Cd(13449) = Cr(28)*MGlpTmSt2(2) + Cd(251)*UStop2(1,2)

	Cd(13450) = Cd(251)*MGlpTmSt2(2) + 2*MGlpTmSt4(2)*UStop2(1,2)

	Cd(13451) = Cd(13449)*MGlpTmSt2(2) + Cd(13450)*UStop2(1,1)

        Cd(13453) = Cd(13066)*MGlpTmSt2(2) - 
     &   Cd(12700)*Cd(13452)*MGlpTmSt4(2)*UStop2(1,1)

        Cd(13454) = Cd(13451) + 
     &   Cd(12701)**2*MGlpTmSt2(2)*UCStopC(1,3)**2

        Cd(13481) = MSq2Diff(1,1)*
     &   (Cd(13480) - MSq2Diff(1,1)*UStop2(1,2))

	Cd(13484) = 2*A0delq(MSbot2(1),Q) - Cd(13483)*MSbot2(1)

	Cd(13498) = Cd(13497)*MSq2Diff(1,1) + Cd(13488)*MSq2Diff(1,2)

        Cd(13499) = (MSbot2(1)*
     &     (Cd(13484)*Cd(13492) - 
     &       4*Cd(13466)*Cd(13498)*MSbot2(1)*MSq2Diff(1,1)*
     &        MSq2Diff(1,2)))/(MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

	Cd(13510) = Cd(13509)*MSbot2(1) - MSq2Diff(1,1)**2*UStop2(1,2)

	Cd(132) = Cd(130) + MGlpTmSt2(2)*UStop2(1,1)

	Cd(13514) = Cd(13512) + 2*Cd(13513)*MSbot2(1)

        Cd(13517) = A0delStop(2)*Cd(13514) - 
     &   2*Cd(13516)*Cr(34)*MSq2Diff(1,1)*MStop2(2)

        Cd(13518) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(A0delq(MSbot2(1),Q)*Cd(13469)*Cd(13470)*
     &         Cr(9)**2*MStop2(1) - 
     &        2*Cd(13517)*UStop2(1,1)*UStop2(1,2)))/
     &    (SW2*TB2*MSq2Diff(1,1)**2*MStop2(1))

        Cd(13519) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(Cd(13469)*Cd(13470)*Cr(9)**2*MSbot2(1)*
     &         MStop2(1) - 
     &        2*Cd(13514)*MStop2(2)*UStop2(1,1)*UStop2(1,2)))/
     &    (SW2*TB2*MSq2Diff(1,1)**2*MStop2(1))

	Cd(13520) = Cd(13514)*Cr(34) - 2*Cd(13516)*MSq2Diff(1,1)

        Cd(13521) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(Cd(13466)*Cd(13469)*Cd(13470)*Cr(9)**2*
     &         MSbot2(1)*MStop2(1) - 
     &        2*Cd(13520)*MStop2(2)*UStop2(1,1)*UStop2(1,2)))/
     &    (SW2*TB2*MSq2Diff(1,1)**2*MStop2(1))

        Cd(13533) = MSq2Diff(1,2)*
     &   (-Cd(13532) + MSq2Diff(1,2)*UStop2(1,1))

        Cd(13534) = 2*Cd(13525)*Cd(13530) + 
     &   Cd(13531)*Cd(13533)*MStop2(2)

        Cd(13545) = Cd(13544)*Conjugate(dY33eps(-1)) + 
     &   Cd(13542)*dY33eps(-1)

        Cd(13546) = Cd(13544)*Conjugate(dY33eps(0)) + 
     &   Cd(13542)*dY33eps(0)

        Cd(13547) = Cd(13544)*Conjugate(dY33eps(1)) + 
     &   Cd(13542)*dY33eps(1)

	Cd(13548) = 2*Cd(13547) + Cd(13546)*Cr(13)

	Cd(13549) = 1/2.D0*(A0delStop(1)*Cd(13545)) + Cd(13548)*MStop2(1)

        Cd(13553) = Cd(13552)*Conjugate(dY33eps(-1)) + 
     &   Cd(13551)*dY33eps(-1)

        Cd(13554) = Cd(13552)*Conjugate(dY33eps(0)) + 
     &   Cd(13551)*dY33eps(0)

        Cd(13555) = Cd(13552)*Conjugate(dY33eps(1)) + 
     &   Cd(13551)*dY33eps(1)

	Cd(13556) = 2*Cd(13555) + Cd(13554)*Cr(34)

	Cd(13557) = 1/2.D0*(A0delStop(2)*Cd(13553)) + Cd(13556)*MStop2(2)

        Cd(13568) = (Cd(13567)*MSq2Diff(1,1) + 
     &     Cd(13564)*MSq2Diff(1,2))/(MSq2Diff(1,1)*MSq2Diff(1,2))

	Cd(13569) = Cd(13559)*Cd(13565) + Cd(13566)*Cd(13568)*MSbot2(1)

	Cd(13589) = Cd(13587)*MSq2Diff(1,1) + Cd(13588)*MSq2Diff(1,2)

	Cd(13594) = A0delStop(1)*Cd(13591) - 4*Cd(12862)*MSq2Diff(1,1)

        Cd(13598) = 2*Cd(12867)*MSq2Diff(1,1) + 
     &   Cd(13591)*Cr(13)*MStop2(1)

        Cd(13603) = -(Cd(13602)*MSbot2(1)) + 
     &   MSq2Diff(1,2)**2*UStop2(1,1)

        Cd(13606) = -(A0delStop(2)*Cd(13603)) - 
     &   Cd(13604)*Cr(34)*MSq2Diff(1,2)*MStop2(2)

        Cd(13607) = Cd(13606) - 
     &   A0delq(MSbot2(1),Q)*Cd(13504)*Cd(13505)*MStop2(2)

	Cd(13608) = -Cd(13603) - Cd(13504)*Cd(13505)*MSbot2(1)

	Cd(13609) = -2*(Cd(13603)*Cr(34) + Cd(13604)*MSq2Diff(1,2))

        Cd(133) = M_3C*Cd(132)*UCStop(1,3) + 
     &   M_3*Cd(130)*UCStopC(1,3)

        Cd(13610) = Cd(13609) - 
     &   2*Cd(13466)*Cd(13504)*Cd(13505)*MSbot2(1)

        Cd(13631) = Cr(9)*
     &   (MTy*Cd(13629) - Cd(13630)*UStop2(1,1)*UStop2(1,2))

	Cd(13640) = 3*MTy2*Cd(13639) - 3*MSq2Diff(1,1)*MSq2Diff(1,2)

        Cd(13645) = MTy2*Cd(13644) - 
     &   3*MSq2Diff(1,1)*MSq2Diff(1,2)*UUStop(3,3)

        Cd(13646) = 2*MTy2*MSq2Diff(1,2)*UStop(2,1) - 
     &   Cd(13645)*UStopC(1,2)

	Cd(13677) = Cd(13481) - Cd(13660)*Cd(13668)

	Cd(13678) = Cd(13670)*Cd(13676) + 4*Cd(13677)*MStop2(1)

	Cd(1362) = Cd(1358) - 4*Cd(902)*UStop2(1,1)

	Cd(13693) = -Cd(13533) - Cd(13679)*Cd(13685)

	Cd(13694) = Cd(13687)*Cd(13692) + 4*Cd(13693)*MStop2(2)

	Cd(13700) = MGlmT2*MGlpT2*Cd(13698) - 2*Cd(13699)*MSbot2(1)

	Cd(13707) = MTy2 + Cd(13706)*UStop2(1,1)

        Cd(13710) = 2*MUE*TB*MSq2Diff(1,2)*UStop(1,1) - 
     &   MTy*Cd(48)*UStop(1,2)

	Cd(13712) = MGl2*MUE*TB*Cd(13711) - AtC*Cd(194)*MSq2Diff(1,2)

	Cd(13715) = Cd(13714)*UStop2(1,1) + MTy*Cd(13504)*UStopC(2,1)

	Cd(13717) = Cd(179)*Cd(13715) - 2*MSq2Diff(1,2)*MStop2(2)

        Cd(13718) = -(Cd(13717)*UCStopC(1,3)) + 
     &   MTy*Cd(13716)*UStop2(1,1)

        Cd(13720) = 2*Cd(194)*MSq2Diff(1,2) - 
     &   Cd(13719)*MStop2(2)**2*UStop2(1,1)

	Cd(13723) = 2*Cr(56)*MStop2(2) + Cd(13722)*UStop2(1,1)

	Cd(1367) = 5*MStop2(1)**2 + Cd(607)*MStop2(3)

	Cd(13729) = Cd(13725)*MStop2(2) - Cd(13728)*UStop2(1,1)

	Cd(13732) = Cd(13725)*MStop2(2) - Cd(13731)*UStop2(1,1)

        Cd(13735) = 2*Cr(87)*MSq2Diff(1,2)**2 + 
     &   MUE2*TB2*Cd(13734)*UStop2(1,1)**2

	Cd(1368) = 5*MStop2(2)*MStop2(3) + 4*Cd(1367)*UStop2(1,2)**2

	Cd(13740) = -2*MGl2*Cd(13739) + Cd(13738)*MStop2(2)

	Cd(13741) = -2*Cd(13737)*MSbot2(1) + Cd(2365)*MStop2(2)

	Cd(13747) = -(Cd(13746)*MStop2(2)) + Cd(13728)*UStop2(1,1)

	Cd(1369) = Cd(1296) - Cd(1368)

	Cd(13749) = -(Cd(13748)*MStop2(2)) + 2*At*MUE*TB*UStop2(1,1)

	Cd(1370) = -2*Cd(852) + Cd(1264)

        Cd(13760) = -2*Cd(13759)*MSbot2(1) + 
     &   Cd(152)*MGlpTmSt2(1)*MStop2(1)

	Cd(13761) = Cd(47)*MStop2(1) + Cd(13760)*UStop2(1,2)**2

	Cd(13762) = -2*Cd(13759)*MSbot2(1) + 2*MTy2*Cd(164)*MStop2(1)

	Cd(13763) = Cd(13760) - Cd(13762)

        Cd(13766) = 2*MUE*TB*Cd(13765)*UStop(1,1) - 
     &   MTy*Cd(44)*UStop(1,2)

	Cd(13768) = -(MGl2*Cd(46)*MStop2(1)) + Cd(13760)*UStop2(1,2)

        Cd(13769) = Cd(44)*Cd(13766)*MStop2(1) + 
     &   Cd(13077)*Cd(13760)*UStop2(1,1)

        Cd(13773) = 2*MSq2Diff(1,1)*MStop2(1) - 
     &   Cd(152)*Cd(13772)*UStopC(1,2)

	Cd(13774) = Cd(13771) - Cd(12701)*UStop2(1,1)

        Cd(13775) = -(Cd(13773)*UCStop(1,3)) + 
     &   MTy*Cd(152)*Cd(13774)*UStop2(1,1)

	Cd(13777) = 2*MSq2Diff(1,1) + Cd(13541)*UStop2(1,2)

        Cd(13778) = MUE*TB*MStop2(1)*UStop2(1,1) - 
     &   Cd(152)*Cd(13776)*UUStopC(1,1)

        Cd(13783) = -(At*Cd(12700)*Cd(13782)*MStop2(1)) + 
     &   2*Cd(13781)*UStop2(1,1)

	Cd(13785) = -2*Cd(10367)*MSbot2(1) + Cd(13784)*MStop2(1)

	Cd(13789) = -(At*AtC*Cd(13785)) + 2*Cd(13788)*MSq2Diff(1,1)**2

        Cd(13790) = At*Cd(12700)*Cd(13785) + 
     &   2*Cd(10367)*MSq2Diff(1,1)**2

	Cd(13791) = Cd(10367)*MSbot2(1) - MStop2(1)*UStop2(1,2)

	Cd(13795) = -2*Cd(13794)*MSbot2(1) + At*Cd(12700)*UStop2(1,2)

	Cd(13798) = Cd(13797)*MStop2(1) - MSbot2(1)*UStop2(1,1)

	Cd(13799) = At*AtC*Cd(13791) - Cd(13798)*MSq2Diff(1,1)

        Cd(13802) = MSq2Diff(1,1)*
     &   (M_3C*Cd(13775) - Cd(585)*MStop2(1)**2)

	Cd(13826) = MSq2Diff(1,2)*(MGl2*MGlmT2 - Cd(13825)*MStop2(2))

	Cd(13830) = (AtC + MUE*TB)*Cd(13826)

	Cd(13832) = MGl2*MGlmT2 - Cd(48)*MStop2(2)

	Cd(13833) = AtC*Cd(194) - MUE*TB*Cd(13832)

	Cd(13834) = MSq2Diff(1,2)*(Cd(13833) + 2*MGl2*MUE*TB*MStop2(2))

	Cd(13835) = MSq2Diff(1,2)*(Cd(194) - MStop2(2)**2*UStop2(1,1))

        Cd(13836) = Cd(13834) - 
     &   AtC*MSq2Diff(1,2)*MStop2(2)**2*UStop2(1,1)

	Cd(13838) = Cd(13725)*MSq2Diff(1,2) + 2*Cd(13728)*UStop2(1,1)

	Cd(13853) = MSq2Diff(1,1)*(MGl2*MGlmT2 - Cd(13852)*MStop2(1))

	Cd(1380) = -(MTy2*Cd(1379)) + MStop2(2)*UStop2(1,2)

	Cd(13860) = 2*At*Cd(10367)*Cd(12700) + Cd(13786)*MSq2Diff(1,1)

        Cd(13874) = Cd(13567)*MSq2Diff(1,1) + 
     &   Cd(12700)*Cd(13560)*MSq2Diff(1,2)

	Cd(13875) = Cd(13873) + Cd(13874)*MSq2Diff(1,1)

        Cd(13882) = Cd(13637)*UStopC(2,1) - 
     &   Cd(12700)*Cd(13633)*UStopC(2,2)

	Cd(13883) = MTy*Cd(13632)*UStopC(2,1) + Cd(13636)*UStopC(2,2)

	Cd(13889) = MTy*Cd(13632)*UUStopC(3,1) + Cd(13636)*UUStopC(3,3)

        Cd(13890) = Cd(12700)*Cd(13633)*UUStopC(3,2) - 
     &   Cd(13637)*UUStopC(3,4)

	Cd(13891) = Cd(13889) - Cd(13890)

        Cd(13892) = 2*Cd(13504)*Cd(13505)*Cr(9)*MSq2Diff(1,1) + 
     &   Cd(13891)*MSq2Diff(1,2)

	Cd(13897) = 2*Cd(9948)*Cd(13896) + MSq2Diff(1,1)

	Cd(1384) = Cd(1299) - Cd(1373)

	Cd(13902) = TB*Cd(13901) + 2*TB2*Cd(424)*UStop2(1,1)

	Cd(13904) = MSq2Diff(1,2) + 2*Cd(13903)*UStop2(1,1)

	Cd(13905) = TB*Cd(13904) + 2*TB2*Cd(424)*UStop2(1,1)

	Cd(13908) = MSq2Diff(1,1) - 2*Cd(13907)*UStop2(1,1)

	Cd(13909) = TB*Cd(13908) - 4*TB2*Cd(10224)*UStop2(1,1)

	Cd(13919) = TB*Cd(13918) + MUE*TB2*Cd(13914)*UStop2(1,1)

        Cd(13937) = -(Cd(13715)*UCStopC(1,3)) + 
     &   MTy*Cd(13936)*UStop2(1,1)

	Cd(13938) = Cd(48)*Cd(13934) - MGl2*Cd(13935)

	Cd(13939) = MTy2*Cd(251)*Cd(13937) - Cd(13938)*UStopC(1,1)

	Cd(13943) = Cd(13941)*MSq2Diff(1,2) - Cd(13942)*UStop2(1,1)

	Cd(13958) = 2*Cd(49)*MSbot2(1) - MGl2*MStop2(2)

	Cd(13959) = Cd(48)**2*Cd(13734) + MGl2*Cd(13958)

	Cd(13960) = 3*MGl2*MTy2 + At*MUE*TB*Cd(1496)

	Cd(1391) = -Cd(523) + Cd(1390)

        Cd(13968) = MTy*MUE*TB*Cd(13734)*UStop(1,1) + 
     &   MSq2Diff(1,2)*MStop2(2)*UStop(1,2)

	Cd(13974) = Cd(13973)*MStop2(2) + Cd(13972)*UStop2(1,1)

	Cd(13975) = 2*Cd(13941)*MStop2(2) + Cd(13942)*UStop2(1,1)

	Cd(13981) = 2*Cd(13788)*MStop2(2) + Cd(13942)*UStop2(1,1)

	Cd(13985) = Cd(13984)*MStop2(2) + 2*Cd(13942)*UStop2(1,1)

	Cd(13987) = Cd(13986)*MStop2(2) + 2*Cd(13942)*UStop2(1,1)

	Cd(13990) = Cd(13989)*MStop2(2) + Cd(13942)*UStop2(1,1)

	Cd(13993) = -(Cd(12273)*MSbot2(1)) + Cd(10025)*MStop2(2)

	Cd(1394) = 6*Cd(874) - 2*MStop2(3)

        Cd(13998) = MUE*TB*Cd(13734)*UCStop(1,3) - 
     &   Cd(13997)*UCStopC(1,3)

        Cd(14000) = MUE*TB*MStop2(2)*UCStop(1,3) + 
     &   Cd(13999)*UCStopC(1,3)

	Cd(14006) = -2*Cd(45)*MSbot2(1) + MGl2*MStop2(1)

	Cd(14007) = Cd(44)**2*Cd(14005) - MGl2*Cd(14006)

        Cd(14009) = (MTy*Cd(9948)*Cd(13469)*Cd(14007)*
     &      MStop2(3) - Cd(47)*Cd(14008)*MSq2Diff(1,1)*MStop2(1)*
     &      UStop(1,1))/(MTy*Cd(47)*MStop2(1)*MStop2(3))

	Cd(14011) = MTy*Cd(13473) + Cd(14010)

        Cd(14012) = MTy2*Cd(164)*Cd(14011) - 
     &   MGlpTmSt2(1)*MSq2Diff(1,1)*MStop2(1)

	Cd(14021) = 2*MStop2(2) + Cd(14020)*UStop2(1,2)

	Cd(14038) = MGlmT2*MTy2 + Cd(116)*MStop2(1)

        Cd(14039) = -2*Cd(14038)*MSbot2(1) + 
     &   Cd(164)*MGlpTmSt2(1)*MStop2(1)

	Cd(14042) = MGlmT2*MGlpT2*Cd(13698) + Cd(14041)*MStop2(1)

	Cd(14048) = -Cd(47) + Cd(118)**2

        Cd(14050) = 2*Cd(9950)*MSq2Diff(1,1) + 
     &   Cd(47)*Cd(14047)*UStop2(1,2)

        Cd(14051) = (2*TB*Cd(14012)*Cr(5) + 
     &     MUE*TB2*Cd(47)*MStop2(1)*UCStop(1,3))/
     &   (TB*TB2*Cd(47)*MStop2(1))

	Cd(14053) = Cd(44)*MSq2Diff(1,1) + MUE2*TB2*MStop2(1)

	Cd(1400) = -4*Cd(874) + Cd(1394)

	Cd(14060) = Cd(14059)*MStop2(3) + Cd(47)*UStop2(1,1)

        Cd(14063) = -(MUE*TB2*Cd(14062)*MStop2(3)) + 
     &   TB*Cd(14060)*UStop2(1,2)

        Cd(14066) = 2*MUE2*TB2*UStop2(1,2)**2 - 
     &   2*AtC*Cd(14009)*UStopC(1,2)

	Cd(14068) = MSq2Diff(1,1)*(MGlmT2*MTy2 + Cd(116)*MStop2(1))

	Cd(14072) = -(Cd(14038)*MSbot2(1)) + Cd(14071)*MStop2(1)

	Cd(14073) = (MGlmT2*MTy2*Cd(13698) + Cd(14072))*MSq2Diff(1,1)

	Cd(14078) = MUE2*TB*Cd(14077) + At*MGlmT2*MUE*MGlpTmSt2(1)

	Cd(14079) = TB2*Cd(14078) - MGlmT2*TB*Cd(44)*MSq2Diff(1,1)

	Cd(14082) = TB*Cd(44)*MSq2Diff(1,1) + TB2*Cd(14081)*MStop2(1)

	Cd(14083) = 2*TB*MSq2Diff(1,1)*MStop2(1) - Cd(14082)*UStop2(1,2)

	Cd(14085) = TB*MSq2Diff(1,1) - Cd(14084)*UStop2(1,2)

        Cd(14086) = -2*MTy*MUE*TB2*Cd(9948)*UCStop(1,3) + 
     &   Cd(14085)*UStop2(1,2)

        Cd(14090) = 2*MSq2Diff(1,1)*MStop2(1)**2 - 
     &   Cd(47)*Cd(14047)*UStop2(1,2)

	Cd(14099) = 2*MSq2Diff(1,2) - Cd(14098)*UStop2(1,2)

        Cd(14123) = Cd(47)*Cd(9948)*Cd(13896)*MGlpTmSt2(1) + 
     &   Cd(164)*Cd(9950)*MSq2Diff(1,1)

	Cd(1408) = Cd(778) - 2*MStop2(1)

	Cd(14138) = TB*Cd(14136) - TB2*Cd(14137)

	Cd(14142) = Cd(47) - MGlpTmSt2(1)*MSq2Diff(1,1)

	Cd(14145) = Cd(164)*MSq2Diff(1,1) + 2*Cd(47)*UStop2(1,1)

        Cd(14152) = Cd(9948)*Cd(14149)*MSq2Diff(1,1) + 
     &   At*Cd(14151)*UStop2(1,2)

        Cd(14154) = Cd(9948)*Cd(14153)*MSq2Diff(1,1) + 
     &   AtC*Cd(14151)*UStop2(1,2)

        Cd(14155) = MUE2*TB2*Cd(47) - 
     &   Cd(14047)*MGlpTmSt2(1)*MSq2Diff(1,1)

        Cd(14159) = TB*Cd(164)*MSq2Diff(1,1)*UStop2(1,1) - 
     &   Cd(14158)*MGlpTmSt2(1)*UStop2(1,2)**2

	Cd(14161) = TB*Cd(14160) + TB2*Cd(424)*UStop2(1,2)**2

        Cd(14164) = Cd(164)*MSq2Diff(1,1) - 
     &   Cd(13896)*MGlpTmSt2(1)*UStop2(1,2)

	Cd(14165) = Cd(47)*Cd(14163) + Cd(14164)*MSq2Diff(1,1)

        Cd(14169) = Cd(47)*Cd(14168) + 
     &   MGlpTmSt2(1)*MSq2Diff(1,1)*MStop2(3)*UStop2(1,2)

        Cd(14170) = 2*AtC*Cd(14169)*Cr(5) - 
     &   Cd(14142)*MStop2(1)*MStop2(3)*UCStopC(1,3)

        Cd(14176) = (2*Cd(47)*Cd(10169)*MSq2Diff(1,1) + 
     &     Cd(14142)*Cd(14175)*MStop2(3))/(Cd(47)*MStop2(3))

        Cd(14180) = (2*Cd(47)*Cd(10176)*MSq2Diff(1,1) + 
     &     Cd(14142)*Cd(14178)*MStop2(3))/(Cd(47)*MStop2(3))

        Cd(14188) = MGl2*TB*Cd(251)*MSq2Diff(1,2) + 
     &   Cd(14187)*MGlpTmSt2(2)*UStop2(1,1)

        Cd(14254) = -(Cd(14253)*MGlpTmSt2(2)) + 
     &   2*Cd(251)*Cr(5)*MSq2Diff(1,2)

        Cd(14276) = Cd(14275)*UCStop(1,3) + 
     &   MTy*Cd(12701)*Cd(14274)*UStop2(1,2)

	Cd(14281) = MUEC*TB*Cd(14275) - At*Cd(14280)

        Cd(14286) = -(Cd(14285)*UCStop(1,3)) + 
     &   MTy*Cd(12701)*Cd(14284)*UStop2(1,1)

        Cd(140) = 2*Cd(71)*UStop(1,1) - 
     &   2*M_3*MTy*MGlpTmSt2(2)*UStop(1,2)

	Cd(1425) = 2*Cd(1423) - Cd(1424)*UStop2(1,1)

	Cd(14306) = Cd(14302) - Cd(14304)

	Cd(14307) = Cd(14305)*MSq2Diff(1,2) - AtC*Cd(14306)*UStopC(1,2)

	Cd(14323) = Cd(14318) + Cd(14322)

	Cd(14324) = Cd(14307) + Cd(14323)*UStopC(1,1)

	Cd(14325) = Cd(14310)*MSq2Diff(1,1) + Cd(14324)*MSq2Diff(1,2)

	Cd(14330) = Cd(13526)*Cd(14329) + 2*Cd(14328)*MSq2Diff(1,2)

        Cd(14336) = Cd(14318)*MSq2Diff(1,1)*UStopC(1,1) - 
     &   2*Cd(13469)*Cd(13663)*MSbot2(1)*MSq2Diff(1,2)*UUStopC(1,3)

	Cd(14339) = Cd(14336) + Cd(14338)*MSq2Diff(1,1)

	Cd(14340) = Cd(14335)*MSq2Diff(1,1)**2 + Cd(14339)*MSq2Diff(1,2)

        Cd(14342) = (-Cd(14340) + 
     &     2*Cd(13469)*Cd(14341)*MSbot2(1)*MSq2Diff(1,2)**2)/
     &   (MSbot2(1)*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

	Cd(1431) = 2*Cd(1430) - 3*Cd(786)*UStop2(1,2)**2

	Cd(1436) = Cd(1435) + Cd(1431)*UStop2(1,1)

	Cd(1438) = Cd(1436)*UStop2(1,1) + Cd(1437)*UStop2(1,2)

	Cd(1444) = 3*MStop2(1) - Cd(1360)*MStop2(2)

	Cd(1457) = Cd(1375)*MStop2(2) - 2*Cd(814)*UStop2(1,2)

	Cd(1458) = 4*MTy2*Cd(1457) + MStop2(2)*MStop2(3)

	Cd(1461) = MTy2*Cd(1460) - MStop2(2)*UStop2(1,2)

	Cd(145) = 4*Cd(144) - Cd(140)*UStopC(1,1)

        Cd(1490) = Cd(52)*MGlpTmSt2(1)*MStop2(3) + 
     &   Cd(56)*T134(MGl2,MTy2,MStop2(1))

        Cd(1491) = Cd(64)*MGlpTmSt2(2)*MStop2(3) - 
     &   Cd(69)*T134(MGl2,MTy2,MStop2(2))

	Cd(1495) = Cd(44)**2 + MGl2*Cd(1494)

	Cd(1499) = Cd(44)**3*Cd(89) + MGl2*Cd(1498)

	Cd(1500) = Cd(47)*Cd(56) + Cd(1499)*UStop2(1,2)

	Cd(1501) = -((-2*MTy2 + Cd(843))*MStop2(3))

	Cd(1505) = Cd(44)**2*Cd(1501) + MGl2*Cd(1504)*MStop2(3)

        Cd(1506) = 2*Cd(47)*MGlpTmSt2(1)*MStop2(3) + 
     &   Cd(1505)*UStop2(1,2)

	Cd(1510) = Cd(48)**2 + MGl2*Cd(1509)

	Cd(1514) = Cd(48)**2*Cd(1511) - MGl2*Cd(1513)*MStop2(3)

	Cd(1519) = Cd(48)**3*Cd(89) + MGl2*Cd(1518)

	Cd(1527) = 2*M_3*MTy*UCStopC(1,3) + Cd(96)*UStop2(1,1)

	Cd(1528) = Cd(96)**2*Cd(1527) + Cd(1526)*UStop2(1,2)

	Cd(1537) = Cd(48)*Cd(1534) + 2*MGl2*Cd(1536)

	Cd(1538) = 2*Cd(48)**3*Cd(1529) - MGl2*Cd(1537)

	Cd(1549) = Cd(48)*Cd(1534) - MGl2*Cd(1548)

	Cd(1550) = 2*Cd(48)**3*Cd(1529) - MGl2*Cd(1549)

	Cd(1553) = Cd(1547) - 3*MGl2*Cd(1552)

	Cd(1554) = 3*Cd(48)**3*Cd(89) - MGl2*Cd(1553)

	Cd(1564) = Cd(48)*Cd(1558) + MGl2*Cd(1563)

	Cd(1565) = Cd(48)**3*Cd(1562) + MGl2*Cd(1564)

	Cd(1573) = Cd(1559)*Cd(1560) + MGl2*Cd(1572)

	Cd(1574) = Cd(48)*Cd(1558) - MGl2*Cd(1573)

	Cd(1575) = Cd(48)**3*Cd(1561) + MGl2*Cd(1574)

	Cd(1576) = -6*MGl2 + Cd(1568)

	Cd(1577) = 2*Cd(1571) + MGl2*Cd(1576)

	Cd(1578) = Cd(48)**2*Cd(1567) - MGl2*Cd(1577)

	Cd(14) = CA*Cd(11)**2*Cd(12) - 4*Pi*(Alfas2L*SA*Cd(13)*Cr(9)**2)

	Cd(1588) = Cd(44)*Cd(1586) + MGl2*Cd(1587)

	Cd(1589) = Cd(44)**3*Cd(1562) - MGl2*Cd(1588)

	Cd(1590) = 4*Cd(47)*Cd(58) - 2*Cd(1589)*UStop2(1,2)

	Cd(157) = -MGl2 + Cd(45)

	Cd(1599) = 2*Cd(44)*Cd(1586) - MGl2*Cd(1598)

	Cd(1600) = 2*Cd(44)**3*Cd(1562) - MGl2*Cd(1599)

	Cd(1601) = Cd(1590) - Cd(1600)*UStop2(1,1)

	Cd(1608) = -2*M_3*MTy*UCStopC(1,3) + Cd(86)*UStop2(1,1)

	Cd(1609) = Cd(86)**2*Cd(1608) + Cd(1607)*UStop2(1,2)

	Cd(1616) = Cd(44)*Cd(1612) - MGl2*Cd(1615)

	Cd(1617) = Cd(44)**3*Cd(1529) - MGl2*Cd(1616)

	Cd(1618) = 2*Cd(47)*MGlpTmSt4(1) - Cd(1617)*UStop2(1,2)

	Cd(1623) = Cd(1592)*Cd(1613) + MGl2*Cd(1622)

	Cd(1624) = Cd(44)*Cd(1612) + MGl2*Cd(1623)

	Cd(1625) = Cd(44)**3*Cd(1529) - MGl2*Cd(1624)

        Cd(1633) = -2*M_3C*MTy*Cd(1495)*UCStop(1,3) + 
     &   Cd(1607)*UStop2(1,2)

	Cd(1639) = Cd(44)*Cd(1635) - MGl2*Cd(1638)

	Cd(1640) = Cd(44)**3*Cd(1636) - 3*MGl2*Cd(1639)

	Cd(1641) = 3*Cd(47)*MGlpTmSt4(1) - Cd(1640)*UStop2(1,2)**2

        Cd(1642) = Cd(1641)*MStop2(3) + 
     &   2*M_3C*MTy*Cd(1506)*UCStop(1,3)

        Cd(1643) = Cd(1642)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(47)*MGlpTmSt2(1)*MStop2(3)*UStop(1,2)*
     &    UStop2(1,2)

	Cd(162) = MGl2*Cd(161) + Cd(44)*MStop2(1)

	Cd(1648) = MTy2*Cd(1647) + Cd(786)*MStop2(1)

	Cd(163) = MGlpTmSt2(1)*MStop2(1) - Cd(162)*UStop2(1,2)

        Cd(1656) = 2*M_3C*MTy*Cd(1510)*UCStop(1,3) + 
     &   Cd(1526)*UStop2(1,2)

	Cd(1661) = Cd(48)*Cd(1658) + MGl2*Cd(1660)

	Cd(1662) = Cd(48)**3*Cd(1659) + 3*MGl2*Cd(1661)

	Cd(1671) = 3*Cd(48)*Cd(1666) - 2*MGl2*Cd(1670)

	Cd(1672) = 2*Cd(48)**3*Cd(319) - MGl2*Cd(1671)

	Cd(1677) = -MGl2 + Cd(607)

	Cd(1678) = Cd(89)*Cd(1592) + MGl2*Cd(1677)

	Cd(1680) = Cd(44) + Cd(1592)

	Cd(167) = -MGl2 + Cd(44)

	Cd(1695) = -2*MGl2 + Cd(1591)

	Cd(1696) = Cd(1593) - MGl2*Cd(1695)

	Cd(1705) = Cd(48)*Cd(1703) - MGl2*Cd(1704)

	Cd(1706) = Cd(48)**3*Cd(1701) + MGl2*Cd(1705)

	Cd(1707) = Cd(48) + Cd(1559)

	Cd(1723) = -2*MGl2 + Cd(1540)

	Cd(1724) = Cd(1543) - MGl2*Cd(1723)

        Cd(1730) = MGlpTmSt2(2)*MGlpTmSt4(2)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(96)**2*UStop(1,2)*UStop2(1,2)

	Cd(1742) = Cd(48)**2*Cd(1733) - MGl2*Cd(1741)

	Cd(1743) = Cd(48)**4*Cd(187) - MGl2*Cd(1742)

	Cd(1744) = -(MGl2*Cd(1552)) + Cd(1740)

	Cd(1745) = Cd(48)**3*Cd(89) - MGl2*Cd(1744)

	Cd(1746) = Cd(1743)*UStop(1,1) + 4*M_3*MTy*Cd(1745)*UStop(1,2)

	Cd(1753) = 5*MStop2(2)**2 + 3*Cd(637)*MStop2(3)

        Cd(174) = M_3C*Cd(152)*UCStop(1,3) - 
     &   M_3*Cd(157)*UCStopC(1,3)

	Cd(1775) = MGl2**2 + Cd(48)**2

	Cd(1777) = -2*MGl2 + Cd(1773)

	Cd(1778) = 2*Cd(1774) + MGl2*Cd(1777)

	Cd(1779) = Cd(48)**2*Cd(1511) + MGl2*Cd(1778)*MStop2(3)

	Cd(16) = 2*A0delStop(2) - Cd(15)*MStop2(2)

	Cd(1790) = 2*Cd(1786) - MGl2*Cd(1789)

        Cd(1817) = 2*M_3C*Cd(47)*Cd(56)*UCStop(1,3)*
     &    UStop(1,1) - M_3*Cd(1816)*UStop(1,2)

        Cd(1819) = MGlpTmSt2(1)*MGlpTmSt4(1)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(86)**2*UStop(1,2)*UStop2(1,2)

	Cd(1821) = -2*MGl2 + Cd(1503)

	Cd(1822) = 2*Cd(1820) + MGl2*Cd(1821)

	Cd(1823) = Cd(44)**2*Cd(1501) - MGl2*Cd(1822)*MStop2(3)

        Cd(1826) = Cd(44)**3*Cd(1825) + 
     &   8*M_3C*MTy*Cd(47)*MGlpTmSt2(1)*UCStop(1,3)

	Cd(1828) = MGlmT2**2*MGlpT2 - Cd(1827)*MStop2(1)

	Cd(1829) = 2*Cd(1828)*MStop2(3) + Cd(1823)*UStop2(1,2)

        Cd(1831) = Cd(1592)*Cd(1830)*UStop(1,1) - 
     &   12*M_3*MTy*MStop2(1)*UStop(1,2)

        Cd(1834) = Cd(44)*Cd(1833)*UStop(1,1) - 
     &   8*M_3*MTy*MTy2*MStop2(1)*UStop(1,2)

	Cd(1843) = Cd(96)*UCStop(1,3) + 2*M_3*MTy*UStop2(1,2)

	Cd(183) = MGl2*Cd(182) + Cd(48)*MStop2(2)

	Cd(1856) = Cd(1846) - 2*MGl2*Cd(1855)

	Cd(1857) = Cd(1853) + 2*MGl2*Cd(1856)

	Cd(1858) = 2*Cd(48)*Cd(1850) - MGl2*Cd(1857)

	Cd(1859) = 2*Cd(48)**3*Cd(1847) - MGl2*Cd(1858)

	Cd(1862) = -(Cd(86)*UCStop(1,3)) + 2*M_3*MTy*UStop2(1,2)

	Cd(1863) = 3*Cd(45) + 2*MStop2(2)

	Cd(184) = -MGl2 + Cd(49)

	Cd(1873) = -MGl2 + Cd(1863)

	Cd(1874) = Cd(1865) - 2*MGl2*Cd(1873)

	Cd(1875) = Cd(1871) + 2*MGl2*Cd(1874)

	Cd(1876) = Cd(44)*Cd(1868) - MGl2*Cd(1875)

	Cd(1877) = 2*Cd(44)**3*Cd(1847) - MGl2*Cd(1876)

	Cd(1890) = MGl2 - 3*Cd(45)

	Cd(1896) = 6*Cd(44)**3*Cd(89) + 2*MGl2*Cd(1895)

	Cd(1900) = -MGl2 + Cd(1897)

	Cd(1901) = Cd(1899) + MGl2*Cd(1900)

	Cd(1902) = Cd(44)**2*Cd(48) - MGl2*Cd(1901)

	Cd(1909) = MGl2 - 2*Cd(49)

	Cd(1910) = Cd(48)**2 + MGl2*Cd(1909)

	Cd(1913) = Cd(48)**2*Cd(1567) + MGl2*Cd(1912)

	Cd(1914) = Cd(96)*UStop(1,1) - 2*M_3*MTy*UStop(1,2)

        Cd(1915) = 2*M_3C*MTy*Cd(1913)*UCStop(1,3) + 
     &   3*Cd(1910)*MGlpTmSt2(2)*MStop2(3)*UStop2(1,2)

	Cd(1918) = 4*MGlpT2*MTy2 + Cd(1917)*MStop2(1)

	Cd(190) = Cd(189)*MStop2(2) - Cd(48)*Cd(187)*UStop2(1,1)

	Cd(1924) = 21*Cd(1922) - Cd(1923)*MStop2(1)

        Cd(192) = M_3C*Cd(179)*UCStop(1,3) - 
     &   M_3*Cd(184)*UCStopC(1,3)

	Cd(1946) = Cd(44)**3*Cd(1942) + MGl2*Cd(1945)

	Cd(1948) = 2*MTy2*Cd(1947) + MStop2(3)*MStop2(4)

	Cd(1970) = 2*Cd(44)*Cd(86) - Cd(1495)*UStop2(1,1)

	Cd(1971) = Cd(109) - Cd(1970)*UStop2(1,1)

	Cd(1973) = MTy2*Cd(1972) + MStop2(1)*MStop2(3)

	Cd(195) = MGlpTmSt2(2)*MStop2(2) - Cd(194)*UStop2(1,2)

	Cd(1982) = Cd(44)**2*Cd(1644) + MGl2*Cd(1981)

	Cd(18) = CA*Cd(12) - SA*Cd(13)

	Cd(1985) = MGlmT2**2 - Cd(57)*MStop2(1)

	Cd(1992) = Cd(47)*Cd(1983) + Cd(1982)*UStop2(1,1)**3

	Cd(197) = M_3C*MTy*Cd(196)*UCStop(1,3) + Cd(195)*UStop2(1,2)

	Cd(1996) = Cd(48)**2*Cd(1501) + MGl2*Cd(1513)*MStop2(3)

	Cd(2001) = MTy2*Cd(2000) - 3*MStop2(2)**2*MStop2(3)

	Cd(2008) = MGlmT2**2 - Cd(70)*MStop2(2)

	Cd(2017) = 2*Cd(48)*Cd(96) - Cd(1510)*UStop2(1,1)

	Cd(2018) = -Cd(129) - Cd(2017)*UStop2(1,1)

	Cd(2022) = Cd(48)**3*Cd(1942) + MGl2*Cd(2021)

	Cd(200) = Cd(44)*Cd(48) - MGl2*Cd(199)

	Cd(2029) = -(MGl2*Cd(637)) + 3*Cd(874)*MStop2(2)

        Cd(201) = 2*Cd(200)*UStop(1,1) + 
     &   4*M_3*MTy*MStop2(3)*UStop(1,2)

	Cd(2060) = Cd(44)*Cd(2050) - MGl2*Cd(2059)

	Cd(204) = MGl2*Cd(203) - 2*Cd(44)*Cd(48)*MStop2(2)

	Cd(2071) = Cd(44)*Cd(2068) - MGl2*Cd(2070)

	Cd(2072) = Cd(44)**3*Cd(2063) - MGl2*Cd(2071)

        Cd(205) = Cd(204)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(184)*MStop2(3)*UStop(1,2)

	Cd(2075) = Cd(44)*Cd(2066) - MGl2*Cd(2074)

	Cd(2090) = Cd(48)*Cd(2083) - 2*MGl2*Cd(2089)

	Cd(2091) = 2*Cd(48)**3*Cd(2048) - MGl2*Cd(2090)

	Cd(2095) = Cd(48)*Cd(2082) - MGl2*Cd(2094)

	Cd(208) = -MGl2**2 + Cd(48)*Cd(206)

	Cd(2104) = Cd(48)*Cd(2101) + MGl2*Cd(2103)

	Cd(2105) = Cd(48)**3*Cd(1701) + MGl2*Cd(2104)

	Cd(2108) = Cd(48)*Cd(2099) + MGl2*Cd(2107)

	Cd(209) = MGl2*Cd(207) - Cd(208)*MStop2(2)

        Cd(210) = 2*Cd(209)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(184)*MStop2(3)*UStop(1,2)

	Cd(2124) = -MGl2 + Cd(2113)

	Cd(2125) = Cd(2117) - 4*MGl2*Cd(2124)

	Cd(2126) = 2*Cd(2119) + MGl2*Cd(2125)

	Cd(2127) = 2*Cd(44)*Cd(2121) - MGl2*Cd(2126)

	Cd(2128) = Cd(44)**3*Cd(2115) - MGl2*Cd(2127)

	Cd(2130) = Cd(2123) + 2*MGl2*Cd(2129)

	Cd(2131) = Cd(44)**2*Cd(2122) - 2*MGl2*Cd(2130)

        Cd(2132) = Cd(2128)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(2131)*MStop2(3)*UStop(1,2)

	Cd(211) = Cd(197) + M_3C*MTy*Cd(179)*UCStop(1,3)*UStop2(1,1)

	Cd(2141) = Cd(1592)*Cd(2137) + MGl2*Cd(2140)

	Cd(2142) = Cd(44)*Cd(2136) + MGl2*Cd(2141)

	Cd(2159) = -MGl2 + Cd(2088)

	Cd(2160) = Cd(2148) - 4*MGl2*Cd(2159)

	Cd(2161) = 2*Cd(2152) + MGl2*Cd(2160)

	Cd(2162) = 2*Cd(48)*Cd(2156) - MGl2*Cd(2161)

	Cd(2163) = Cd(48)**3*Cd(2115) - MGl2*Cd(2162)

	Cd(2165) = Cd(2158) + 2*MGl2*Cd(2164)

	Cd(2166) = Cd(48)**2*Cd(2157) - 2*MGl2*Cd(2165)

        Cd(2167) = Cd(2163)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(2166)*MStop2(3)*UStop(1,2)

	Cd(2176) = Cd(1559)*Cd(2170) + MGl2*Cd(2175)

	Cd(2177) = Cd(48)*Cd(2172) + MGl2*Cd(2176)

	Cd(2178) = Cd(48)**3*Cd(2168) - MGl2*Cd(2177)

	Cd(2179) = -4*MGl2 + Cd(2173)

	Cd(2180) = 2*Cd(2174) + MGl2*Cd(2179)

	Cd(2181) = Cd(48)**2*Cd(1529) - MGl2*Cd(2180)

	Cd(2182) = Cd(2178)*UStop(1,1) - 2*M_3*MTy*Cd(2181)*UStop(1,2)

        Cd(2183) = Cd(64)*Cd(2182) + 
     &   2*Cd(96)**2*Cd(97)*Cd(1508)*MStop2(3)

        Cd(2184) = Cd(2183)*MStop2(3) - 
     &   Cd(2167)*T134(MGl2,MTy2,MStop2(2))

        Cd(2188) = Cd(44)**3*UStop(1,1) - 
     &   M_3*MTy*Cd(2187)**2*UStop(1,2)

	Cd(2189) = Cd(2188) + MGl2*Cd(2186)*UStop(1,1)

	Cd(2190) = MTy*Cd(86)*UStop(1,1) - 2*M_3*MTy2*UStop(1,2)

        Cd(2191) = -(Cd(86)*Cd(167)*Cd(2190)) + 
     &   2*MTy*Cd(2189)*UStop2(1,2)

	Cd(2209) = Cd(2188) + Cd(2208)*UStop(1,1)

	Cd(2212) = Cd(44)**2*Cd(2206) - MGl2*Cd(2211)

	Cd(2213) = Cd(2209)*MStop2(3) + M_3*MTy*Cd(2212)*UStop(1,2)

	Cd(2214) = MGl2*Cd(2125) + 2*Cd(2194)

	Cd(2215) = 2*Cd(44)*Cd(2121) - MGl2*Cd(2214)

	Cd(2216) = Cd(44)**3*Cd(2115) - MGl2*Cd(2215)

	Cd(220) = AtC*Cd(219)*UCStop(1,3) + At*Cd(154)*UCStopC(1,3)

	Cd(2229) = 2*Cd(2188) + Cd(2228)*UStop(1,1)

	Cd(2230) = Cd(2229)*MStop2(3) + M_3*MTy*Cd(2212)*UStop(1,2)

	Cd(2232) = Cd(2223) + MGl2*Cd(2231)

	Cd(2233) = Cd(44)*Cd(2222) + MGl2*Cd(2232)

	Cd(2234) = Cd(2224) + MGl2*MStop2(3)

	Cd(2235) = Cd(2226) + MGl2*Cd(2234)

	Cd(2236) = Cd(44)**2*Cd(2225) - MGl2*Cd(2235)

        Cd(2242) = Cd(48)**3*UStop(1,1) + 
     &   M_3*MTy*Cd(2241)**2*UStop(1,2)

	Cd(2243) = Cd(2242) + MGl2*Cd(2240)*UStop(1,1)

	Cd(2244) = MTy*Cd(96)*UStop(1,1) + 2*M_3*MTy2*UStop(1,2)

	Cd(2262) = Cd(2242) + Cd(2261)*UStop(1,1)

	Cd(223) = M_3C*Cd(179)*UCStop(1,3) - M_3*Cd(49)*UCStopC(1,3)

	Cd(2265) = Cd(48)**2*Cd(2256) - MGl2*Cd(2264)

	Cd(2266) = Cd(2262)*MStop2(3) + M_3*MTy*Cd(2265)*UStop(1,2)

	Cd(2267) = MGl2*Cd(2160) + 2*Cd(2249)

	Cd(2268) = 2*Cd(48)*Cd(2156) - MGl2*Cd(2267)

	Cd(2269) = Cd(48)**3*Cd(2115) - MGl2*Cd(2268)

	Cd(2271) = Cd(48)**3*Cd(2196) - MGl2*Cd(2270)

	Cd(2279) = MTy2*Cd(1560) - MStop2(2)*MStop2(3)

	Cd(2284) = 2*Cd(2242) + Cd(2283)*UStop(1,1)

	Cd(2285) = Cd(2284)*MStop2(3) + M_3*MTy*Cd(2265)*UStop(1,2)

	Cd(2287) = Cd(2276) + MGl2*Cd(2286)

	Cd(2288) = Cd(48)*Cd(2278) + MGl2*Cd(2287)

	Cd(2289) = Cd(48)**3*Cd(2168) - MGl2*Cd(2288)

	Cd(2290) = MGl2*Cd(2280) + Cd(2281)

	Cd(226) = 2*Cd(183)*UStop(1,1) + M_3*MTy*Cd(184)*UStop(1,2)

	Cd(2299) = 4*M_3*MTy2*Cd(44) - MTy*Cd(86)**2*UCStop(1,3)

	Cd(2300) = Cd(44)**2*Cd(86) - MGl2**2*MGlpTmSt2(1)

	Cd(2301) = Cd(2188)*UStop(1,2) + MGl2*Cd(2186)*UUStop(1,3)

	Cd(2302) = Cd(2301) - M_3C*MTy*Cd(1495)*UUStop(1,1)

	Cd(227) = Cd(183)*UStop(1,1) + M_3*MTy*Cd(184)*UStop(1,2)

	Cd(2304) = Cd(86)*Cd(2299) - 2*M_3*MTy2*Cd(1495)*UStop2(1,1)

	Cd(2306) = MGl2 + Cd(1494)*UStop2(1,2)

	Cd(2307) = MGl2*Cd(2306) - Cd(44)**2*UStop2(1,1)

	Cd(2315) = -MGl2 + 4*Cd(45)

	Cd(2316) = 2*Cd(2314) - MGl2*Cd(2315)

	Cd(2320) = Cd(44)**3*Cd(1962) + MGl2*Cd(1966)

	Cd(2324) = Cd(44)*Cd(2121) - MGl2*Cd(2323)

	Cd(2325) = Cd(44)**3*Cd(2115) - 2*MGl2*Cd(2324)

	Cd(2329) = Cd(44)**3*Cd(48) + MGl2*Cd(2328)

	Cd(2332) = MTy2 + 2*Cd(256)

	Cd(2334) = 3*MTy2*Cd(2332) - Cd(2333)*MStop2(1)

	Cd(2335) = -3*MGl2 + Cd(2330)

	Cd(2336) = Cd(2334) + MGl2*Cd(2335)

	Cd(2337) = Cd(44)**2*Cd(2331) - MGl2*Cd(2336)

	Cd(2338) = Cd(2188)*UStop(1,2) + Cd(2208)*UUStop(1,3)

	Cd(2339) = Cd(2338)*MStop2(3) + M_3*MTy*Cd(2337)*UUStop(1,2)

        Cd(2340) = Cd(2339)*MStop2(3) - 
     &   4*M_3C*MTy*Cd(2329)*UUStop(1,1)

	Cd(2341) = Cd(2325)*UCStop(1,3) + 2*Cd(2340)*UUStopC(1,2)

	Cd(2342) = -(MTy*Cd(2341)) + M_3*Cd(2320)*MStop2(3)

	Cd(2357) = 2*M_3*MStop2(1) - M_3C*Cd(607)*UCStop(1,3)**2

        Cd(2363) = M_3*Cd(2359)*MStop2(3) - 
     &   2*M_3C*Cd(2362)*UCStop(1,3)**2

	Cd(2381) = Cd(508)*MStop2(1) + Cd(169)*MStop2(2)

        Cd(2382) = M_3*Cd(2380)*MStop2(3) - 
     &   4*M_3C*Cd(2381)*UCStop(1,3)**2

	Cd(22) = 2*A0delStop(1) - Cd(21)*MStop2(1)

        Cd(2390) = 6*M_3*Cd(1136)*MStop2(1)*MStop2(3) + 
     &   M_3C*Cd(2389)*UCStop(1,3)**2

        Cd(2395) = 4*M_3*Cd(1179)*MStop2(1)*MStop2(3) - 
     &   M_3C*Cd(2394)*UCStop(1,3)**2

	Cd(2399) = Cd(44)*Cd(2136) + MGl2*Cd(2398)

	Cd(2400) = Cd(44)**3*Cd(2134) - MGl2*Cd(2399)

	Cd(2405) = Cd(44)*Cd(2222) - MGl2*Cd(2404)

	Cd(2406) = Cd(44)**3*Cd(2134) - MGl2*Cd(2405)

	Cd(2408) = 2*Cd(2188)*UStop(1,2) + Cd(2228)*UUStop(1,3)

	Cd(2409) = Cd(2408)*MStop2(3) + M_3*MTy*Cd(2337)*UUStop(1,2)

	Cd(2414) = Cd(2413)*MStop2(1) + Cd(508)*MStop2(2)

	Cd(2416) = Cd(2415)*MStop2(1) + Cd(519)*MStop2(2)

	Cd(239) = Cd(44)**2 - MGl2*Cd(45)

	Cd(2430) = MGl2**3*Cd(2365) + Cd(1982)*UStop2(1,1)**3

	Cd(2433) = Cd(788)*MStop2(1) + Cd(2432)*MStop2(2)

	Cd(2438) = M_3*Cd(970) - 4*M_3C*Cd(172)*UCStop(1,3)**2

	Cd(2450) = 4*M_3*MTy2*Cd(48) + MTy*Cd(96)**2*UCStop(1,3)

	Cd(2451) = Cd(48)**2*Cd(96) - MGl2**2*MGlpTmSt2(2)

	Cd(242) = 2*MGl2*MStop2(1) - Cd(239)*UStop2(1,2)

	Cd(2452) = Cd(2242)*UStop(1,2) + MGl2*Cd(2240)*UUStop(1,3)

	Cd(2453) = Cd(2452) + M_3C*MTy*Cd(1510)*UUStop(1,1)

	Cd(2455) = Cd(96)*Cd(2450) - 2*M_3*MTy2*Cd(1510)*UStop2(1,1)

	Cd(2457) = MGl2 + Cd(1509)*UStop2(1,2)

	Cd(2458) = MGl2*Cd(2457) - Cd(48)**2*UStop2(1,1)

	Cd(2463) = Cd(48)*Cd(2172) + MGl2*Cd(2462)

	Cd(2464) = Cd(48)**3*Cd(2168) - MGl2*Cd(2463)

	Cd(2470) = Cd(48)*Cd(2278) - MGl2*Cd(2469)

	Cd(2471) = Cd(48)**3*Cd(2168) - MGl2*Cd(2470)

	Cd(2473) = 3*MTy2*Cd(306) - Cd(2472)*MStop2(2)

	Cd(2476) = -3*MGl2 + Cd(2474)

	Cd(2477) = Cd(2473) + MGl2*Cd(2476)

	Cd(2478) = Cd(48)**2*Cd(2475) - MGl2*Cd(2477)

	Cd(2479) = 2*Cd(2242)*UStop(1,2) + Cd(2283)*UUStop(1,3)

	Cd(2480) = Cd(2479)*MStop2(3) + M_3*MTy*Cd(2478)*UUStop(1,2)

	Cd(2481) = Cd(2480) + M_3C*MTy*Cd(1910)*MStop2(3)*UUStop(1,1)

	Cd(2482) = Cd(2471)*UCStop(1,3) + 2*Cd(2481)*UUStopC(1,2)

	Cd(2485) = Cd(508)*MStop2(1) + Cd(2413)*MStop2(2)

	Cd(2486) = Cd(519)*MStop2(1) + Cd(2415)*MStop2(2)

	Cd(2494) = M_3*Cd(970) - 4*M_3C*UCStop(1,3)**2

	Cd(2499) = Cd(290)*MStop2(1) - Cd(2498)*MStop2(2)

        Cd(2500) = 2*M_3*Cd(1202)*MStop2(2) + 
     &   M_3C*Cd(2499)*UCStop(1,3)**2

        Cd(247) = M_3C*Cd(164)*UCStop(1,3) - 
     &   M_3*Cd(167)*UCStopC(1,3)

	Cd(2517) = MTy2*Cd(579) - Cd(1657)*MStop2(2)

	Cd(2524) = Cd(48)*Cd(2517) - MGl2*Cd(2523)

	Cd(2525) = Cd(44)*Cd(48)**3 - MGl2*Cd(2524)

	Cd(2526) = -2*MGl2*Cd(2159) + Cd(2520)

	Cd(2527) = Cd(2522) + 2*MGl2*Cd(2526)

	Cd(2528) = 2*Cd(48)*Cd(2156) - MGl2*Cd(2527)

	Cd(2529) = 2*Cd(48)**3*Cd(2521) - MGl2*Cd(2528)

        Cd(2530) = -(Cd(2529)*UCStop(1,3)) + 
     &   2*M_3*MTy*Cd(2525)*UStop2(1,2)

	Cd(2532) = Cd(48)**3*Cd(1962) + MGl2*Cd(2042)

	Cd(2536) = Cd(48)*Cd(2156) - MGl2*Cd(2535)

	Cd(2537) = Cd(48)**3*Cd(2115) - 2*MGl2*Cd(2536)

	Cd(2540) = Cd(48)**4 - MGl2*Cd(2539)

	Cd(2541) = Cd(2242)*UStop(1,2) + Cd(2261)*UUStop(1,3)

	Cd(2542) = Cd(2541)*MStop2(3) + M_3*MTy*Cd(2478)*UUStop(1,2)

        Cd(2543) = Cd(2542)*MStop2(3) - 
     &   4*M_3C*MTy*Cd(2540)*UUStop(1,1)

	Cd(2555) = 2*M_3*MStop2(2) - M_3C*Cd(874)*UCStop(1,3)**2

	Cd(2556) = 4*M_3*MStop2(2) - M_3C*Cd(637)*UCStop(1,3)**2

        Cd(2559) = M_3*Cd(2359)*MStop2(3) + 
     &   2*M_3C*Cd(2558)*UCStop(1,3)**2

	Cd(2570) = MStop2(1) + Cd(290)*MStop2(2)

        Cd(2571) = M_3*Cd(2366)*MStop2(3) + 
     &   8*M_3C*Cd(2570)*UCStop(1,3)**2

	Cd(253) = Cd(48)**2 - MGl2*Cd(49)

        Cd(2574) = 3*M_3*Cd(2364)*MStop2(2)*MStop2(3) - 
     &   M_3C*Cd(2573)*UCStop(1,3)**2

	Cd(24) = SA*Cd(13) + CA*Cd(23)

	Cd(254) = -MGl2 + Cd(48)

        Cd(2586) = 4*M_3*Cd(1179)*MStop2(2)*MStop2(3) + 
     &   M_3C*Cd(2585)*UCStop(1,3)**2

	Cd(2602) = 2*Cd(44)**2 + MGl2*Cd(2601)

	Cd(2605) = -(MGl2*Cd(1789)) + Cd(2598)

	Cd(2607) = -(MGl2*Cd(1789)) + Cd(2604)

	Cd(2608) = Cd(44)*Cd(45)*Cd(2193) - MGl2*Cd(2607)

	Cd(257) = M_3C*Cd(251)*UCStop(1,3) - M_3*Cd(48)*UCStopC(1,3)

        Cd(2615) = Cd(2614)*UStop(1,1) + 
     &   M_3*MTy*Cd(2608)*UStop(1,2)*UStop2(1,2)

	Cd(2618) = Cd(2598) + MGl2*Cd(2617)

	Cd(2621) = Cd(44)*Cd(45)*Cd(2193) - MGl2*Cd(2620)

	Cd(2628) = 9*MGl2**2 + 2*Cd(44)**2

	Cd(2630) = -(MGl2*Cd(1789)) + Cd(2604)

	Cd(2631) = Cd(44)*Cd(45)*Cd(2193) - MGl2*Cd(2630)

	Cd(2637) = Cd(2635) - Cd(152)*Cd(2597)*UStop2(1,1)

	Cd(2638) = Cd(2636) - Cd(2637)*UStop2(1,1)

	Cd(260) = 2*MGl2*MStop2(2) - Cd(253)*UStop2(1,2)

	Cd(2645) = 4*Cd(86)*MStop2(1)**2 - Cd(152)*Cd(2597)*UStop2(1,2)

	Cd(2647) = 2*Cd(86)**2*UStop(1,1) + M_3*MTy*Cd(167)*UStop(1,2)

	Cd(2653) = 4*MGlpTmSt2(1)*MStop2(1)**2 + Cd(2608)*UStop2(1,2)

        Cd(2654) = Cd(2652)*UUStop(1,1) + 
     &   2*M_3*MTy*Cd(2653)*UUStop(1,3)

	Cd(262) = M_3C*Cd(261)*UCStop(1,3) + MTy*Cd(260)*UStop2(1,2)

	Cd(2662) = Cd(44)*Cd(45)*Cd(2193) - MGl2*Cd(2661)

	Cd(2666) = 2*M_3*MTy2*MGlpTmSt2(1) - MTy*Cd(86)**2*UCStop(1,3)

	Cd(264) = -MGl2 + 2*Cd(263)

	Cd(2684) = 2*Cd(48)**2 + MGl2*Cd(2601)

	Cd(2688) = -(MGl2*Cd(1855)) + Cd(2682)

	Cd(2690) = -(MGl2*Cd(1855)) + Cd(2687)

	Cd(265) = Cd(48)**2 - MGl2*Cd(264)

	Cd(2691) = Cd(48)*Cd(49)*Cd(2686) - MGl2*Cd(2690)

        Cd(2696) = Cd(2695)*UStop(1,1) - 
     &   M_3*MTy*Cd(2691)*UStop(1,2)*UStop2(1,2)

	Cd(2698) = MGl2*Cd(2617) + Cd(2682)

        Cd(266) = 2*MTy*Cd(265)*UStop(1,1) - 
     &   2*M_3*Cd(96)*MStop2(3)*UStop(1,2)

	Cd(2704) = 9*MGl2**2 + 2*Cd(48)**2

	Cd(2706) = -(MGl2*Cd(1855)) + Cd(2687)

	Cd(2707) = Cd(48)*Cd(49)*Cd(2686) - MGl2*Cd(2706)

	Cd(267) = 2*Cd(49)*MStop2(1) - MGl2*MStop2(2)

	Cd(2713) = Cd(2711) - Cd(179)*Cd(2681)*UStop2(1,1)

	Cd(2714) = Cd(2712) - Cd(2713)*UStop2(1,1)

	Cd(268) = -(MGl2*Cd(267)) + Cd(48)**2*Cr(26)

	Cd(2721) = 4*Cd(96)*MStop2(2)**2 - Cd(179)*Cd(2681)*UStop2(1,2)

	Cd(2723) = Cd(48)*Cd(49)*Cd(2686) - MGl2*Cd(2722)

	Cd(2725) = 2*Cd(96)**2*UStop(1,1) - M_3*MTy*Cd(254)*UStop(1,2)

	Cd(2729) = 4*MGlpTmSt2(2)*MStop2(2)**2 + Cd(2691)*UStop2(1,2)

        Cd(269) = 2*Cd(268)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(254)*MStop2(3)*UStop(1,2)

        Cd(2730) = Cd(2728)*UUStop(1,1) - 
     &   2*M_3*MTy*Cd(2729)*UUStop(1,3)

	Cd(2738) = Cd(48)*Cd(49)*Cd(2686) - MGl2*Cd(2737)

	Cd(270) = MGl2*MStop2(2) - Cd(49)*MStop2(4)

	Cd(2742) = 2*M_3*MTy2*MGlpTmSt2(2) + MTy*Cd(96)**2*UCStop(1,3)

	Cd(271) = MGl2*Cd(270) + Cd(48)**2*MStop2(1)

	Cd(2756) = MGl2*Cd(1890) + Cd(2314)

	Cd(2757) = Cd(44)**2*Cd(45) - MGl2*Cd(2756)

	Cd(2758) = MTy2*Cd(673) - MStop2(1)**2

        Cd(272) = 2*Cd(271)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(254)*MStop2(3)*UStop(1,2)

	Cd(273) = Cd(262) + M_3C*MTy2*Cd(251)*UCStop(1,3)*UStop2(1,1)

	Cd(2772) = -(MGl2*Cd(2771)) + Cd(44)**2*Cd(1561)*MStop2(1)

	Cd(2777) = Cd(2769) - MGl2*Cd(2776)

	Cd(2778) = MGl2*Cd(2777) + Cd(44)**2*Cd(1561)*MStop2(1)

	Cd(274) = Cd(272)*UStop2(1,1) + Cd(269)*UStop2(1,2)

	Cd(2779) = 3*MGl2*Cd(1890) + Cd(2773)

	Cd(2780) = Cd(44)*Cd(2775) - MGl2*Cd(2779)

	Cd(2784) = 2*Cd(45)*MStop2(1) + Cd(155)*MStop2(2)

	Cd(2787) = Cd(2785)*MStop2(1) + Cd(2786)*MStop2(2)

	Cd(275) = -(MTy*Cd(274)) + Cd(266)*MStop2(2)

	Cd(2790) = 2*Cd(45)*Cd(2788)*MStop2(1) + Cd(2789)*MStop2(2)

	Cd(2791) = -(MGl2*Cd(286)) + 3*Cd(2784)

	Cd(2792) = 3*Cd(2787) - MGl2*Cd(2791)

	Cd(2793) = Cd(2790) - MGl2*Cd(2792)

	Cd(2794) = MGl2*Cd(2793) + 3*Cd(44)**2*Cd(48)*MStop2(1)**2

        Cd(2795) = Cd(2794)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(2757)*MStop2(1)*UStop(1,2)

	Cd(2797) = 2*(6*Cd(45)*MStop2(1) + Cd(1684)*MStop2(3))

	Cd(2798) = 4*Cd(2759)*MStop2(1) + 2*Cd(2604)*MStop2(3)

	Cd(276) = 2*Cd(273)*MStop2(3) - Cd(275)*UStopC(1,1)

	Cd(2799) = Cd(2797) - 2*MGl2*MStop2(4)

	Cd(2800) = Cd(2798) - MGl2*Cd(2799)

	Cd(2801) = Cd(44)*Cd(2796) + MGl2*Cd(2800)

	Cd(2804) = (MGl2*Cd(594) + 5*MStop2(1)**2)*MStop2(3)

	Cd(2806) = MGl2*Cd(594) - Cd(2567)*MStop2(1)

	Cd(2814) = -(MGl2*Cd(2813)) + Cd(44)**2*Cd(2812)*MStop2(1)

	Cd(2818) = 2*M_3C*MTy*Cd(152)*UStop(1,1) + Cd(156)*UStop(1,2)

	Cd(2823) = Cd(2822)*MStop2(1) + Cd(2789)*MStop2(2)

	Cd(2824) = 10*MGl2**3 + Cd(44)**2*Cd(1529)

	Cd(2825) = 3*MGl2**3*Cd(2617) - Cd(2824)*MStop2(1)

	Cd(2828) = 2*MTy2**2*Cd(786) - 6*MStop2(1)**2*MStop2(4)

	Cd(2832) = Cd(2830)*MStop2(1) + Cd(2831)*MStop2(2)

	Cd(2833) = Cd(44)*Cd(2828) - 2*MGl2*Cd(2832)

	Cd(2836) = -(MGl2*Cd(607)) + Cd(2835)

	Cd(2837) = 3*Cd(2821) - MGl2*Cd(2836)

	Cd(2838) = Cd(2823) - MGl2*Cd(2837)

	Cd(2839) = MGl2*Cd(2838) + Cd(44)**2*Cd(1529)*MStop2(1)**2

	Cd(2841) = Cd(162)*UStop(1,1) - M_3*MTy*Cd(157)*UStop(1,2)

        Cd(2847) = 2*Cd(2757)*MStop2(1) + 
     &   Cd(152)*Cd(2597)*MStop2(3)*UStop2(1,1)

	Cd(2849) = 2*(Cd(2768)*MStop2(1) + Cd(2598)*MStop2(3))

	Cd(2850) = 2*(Cd(2769)*MStop2(1) + Cd(44)*Cd(2599)*MStop2(3))

	Cd(2852) = -(MGl2**2*Cd(2851)) + Cd(44)**2*Cd(1529)*MStop2(1)

	Cd(2856) = 4*Cd(2773)*MStop2(1) + 2*Cd(2604)*MStop2(3)

	Cd(2857) = -(MGl2*Cd(2849)) + Cd(2850)

	Cd(282) = Cd(44)*Cd(48) + MGl2*Cd(281)

	Cd(2859) = -(MGl2**3*Cd(2858)) + Cd(44)**2*Cd(1529)*MStop2(1)

	Cd(2860) = Cd(2859)*MStop2(1) - MGl2**3*Cd(2826)*MStop2(2)

	Cd(2868) = 2*MTy2*Cd(2635) - Cd(118)**3*MStop2(1)

	Cd(2869) = MTy2*Cd(2867)*MStop2(1) + Cd(2868)*MStop2(3)

	Cd(2873) = 4*Cd(2870)*MStop2(1) + 2*Cd(2872)*MStop2(2)

	Cd(2884) = Cd(2869) - MTy2*Cd(2873)*UStop2(1,1)

	Cd(2891) = Cd(2888) - MGl2*Cd(2890)

	Cd(2892) = MGl2*Cd(2891) + Cd(44)**2*Cd(2063)*MStop2(1)

	Cd(2897) = -(MGl2**3*Cd(2896)) + Cd(44)**2*Cd(2063)*MStop2(1)

	Cd(2901) = 2*MTy2**2*Cd(607) - 2*Cd(1634)*MStop2(1)**2

	Cd(2903) = 4*Cd(2902)*MStop2(1) + 2*Cd(2604)*MStop2(3)

	Cd(2912) = Cd(2911)*MStop2(1) - Cd(2786)*MStop2(2)

	Cd(2914) = Cd(2913)*MStop2(1) - Cd(2789)*MStop2(2)

	Cd(2915) = 24*Cd(45)*MStop2(1) + 2*Cd(1684)*MStop2(3)

	Cd(2917) = -2*MGl2*Cd(607) + Cd(2915)

	Cd(2918) = Cd(2903) - MGl2*Cd(2917)

	Cd(2919) = Cd(44)*Cd(2901) - MGl2*Cd(2918)

	Cd(2927) = MTy2*Cd(2925) - Cd(2926)*MStop2(1)*MStop2(3)

	Cd(2929) = MTy2*Cd(2928) + Cd(1532)*MStop2(1)*MStop2(3)

	Cd(2934) = -2*MGl2 + Cd(2921)

	Cd(2937) = Cd(2929) - MGl2*Cd(2936)

	Cd(2938) = Cd(44)*Cd(2927) - MGl2*Cd(2937)

	Cd(2945) = 12*Cd(45)*MStop2(1) - 2*Cd(1684)*MStop2(3)

	Cd(2947) = 3*Cd(2946) - MTy2*MStop2(2)

	Cd(2948) = 4*Cd(2947)*MStop2(1) - 2*Cd(2604)*MStop2(3)

	Cd(2949) = -(MGl2*Cd(2943)) + Cd(2945)

	Cd(2950) = Cd(2948) - MGl2*Cd(2949)

	Cd(2951) = Cd(44)*Cd(2944) - MGl2*Cd(2950)

	Cd(2952) = 4*Cd(2794)*UStop(1,1) - M_3*MTy*Cd(2951)*UStop(1,2)

	Cd(2961) = Cd(44)*Cd(2954) + MGl2*Cd(2960)

	Cd(2962) = Cd(2645) - 2*MTy2*Cd(164)*MStop2(1)*UStop2(1,2)

	Cd(2963) = Cd(2961)*MStop2(1) + 4*MTy2*Cd(2962)*MStop2(3)

	Cd(2975) = 2*(Cd(2974)*MStop2(1) + Cd(2598)*MStop2(3))

	Cd(2977) = 2*(Cd(2976)*MStop2(1) + Cd(44)*Cd(2599)*MStop2(3))

	Cd(2979) = -(MGl2**3*Cd(2896)) + Cd(44)**2*Cd(2978)*MStop2(1)

	Cd(2982) = Cd(2981)*MStop2(1) - Cd(2831)*MStop2(2)

	Cd(2983) = Cd(44)*Cd(2796) - 2*MGl2*Cd(2982)

	Cd(2989) = 4*MTy2*Cd(44) + MStop2(1)**2

	Cd(2990) = 3*MTy2*Cd(2988) - 4*Cd(2989)*MStop2(1)*MStop2(3)

	Cd(2998) = Cd(2990) + MGl2*Cd(2997)

	Cd(295) = -2*MGl2*Cd(199) + Cd(288)

	Cd(2999) = MGl2*Cd(2998) + Cd(44)**2*Cd(2986)*MStop2(1)

	Cd(3001) = Cd(44)*Cd(2994) - MGl2*Cd(3000)

        Cd(3002) = Cd(2999)*UStop(1,1) - 
     &   M_3*MTy*Cd(3001)*MStop2(1)*MStop2(3)*UStop(1,2)

	Cd(3007) = -(MGl2*Cd(1789)) + Cd(3005)

	Cd(3008) = Cd(44)*Cd(3006) - MGl2*Cd(3007)

	Cd(296) = -2*Cd(295)*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))

	Cd(3012) = 2*Cd(3011)*MStop2(1) - 2*Cd(2872)*MStop2(2)

	Cd(3013) = (MGl2*Cd(607) + 4*MStop2(1)**2)*MStop2(3)

        Cd(3029) = MStop2(3)*
     &   (Cd(3028)*MStop2(1)**3 + 4*Cd(3008)*MStop2(3)*UStop2(1,2))

	Cd(3036) = Cd(44)**2 + MGl2*Cd(3035)

	Cd(3041) = MGlpT2*MTy2 + Cd(3040)*MStop2(1)

	Cd(3046) = 2*(Cd(3045)*MStop2(1) + Cd(2872)*MStop2(2))

	Cd(3051) = -(MGl2*MGlpT2) + Cd(3049)

	Cd(3052) = Cd(3050) - MGl2*Cd(3051)

	Cd(3053) = 4*MGl2**3 + Cd(44)**2*Cd(2064)

	Cd(3054) = MGl2*Cd(3052) - Cd(3053)*MStop2(1)

	Cd(3066) = 14*MGl2**3 + Cd(44)**2*Cd(2134)

	Cd(3067) = 3*MGl2**3*Cd(2617) - Cd(3066)*MStop2(1)

	Cd(3072) = -2*MGl2 + Cd(3070)

	Cd(3078) = Cd(2924) - MGl2*Cd(3077)

	Cd(303) = MGl2 + Cd(89)

	Cd(3079) = MGl2*Cd(3078) + Cd(44)**2*Cd(2134)*MStop2(1)

	Cd(3082) = 3*MGl2*Cd(1890) + Cd(3081)

	Cd(3083) = Cd(44)*Cd(3080) - MGl2*Cd(3082)

	Cd(3085) = MGl2*Cd(3084) - Cd(44)**2*MStop2(1)

	Cd(3097) = 2*Cd(47)*Cd(2874)*MStop2(1) + Cd(2873)*UStop2(1,1)**3

	Cd(3100) = Cd(2314) - MGl2*Cd(3099)

	Cd(3101) = MGl2*Cd(3100) - Cd(44)**2*MStop2(1)

	Cd(3102) = 2*MGl2*Cd(1890) + Cd(2902)

	Cd(3103) = Cd(44)*Cd(3062) - MGl2*Cd(3102)

        Cd(3104) = Cd(3101)*MStop2(3)*UStop(1,1) - 
     &   M_3*MTy*Cd(3103)*UStop(1,2)

	Cd(3114) = -MGl2 + 3*Cd(49)

	Cd(3115) = Cd(3113) - MGl2*Cd(3114)

	Cd(3116) = Cd(48)**2*Cd(49) - MGl2*Cd(3115)

	Cd(3120) = 4*(Cd(3116) - Cd(3119)*UStop2(1,2))

	Cd(3121) = 4*Cd(1523)*Cd(3116) - Cd(3120)*UStop2(1,1)

	Cd(3134) = Cd(48)**2*Cd(1529) + MGl2*Cd(3133)

	Cd(309) = -MGl2 + Cd(285)

	Cd(3142) = Cd(3131) - MGl2*Cd(3141)

	Cd(3143) = MGl2*Cd(3142) + Cd(48)**2*Cd(1529)*MStop2(2)

	Cd(3144) = -3*MGl2*Cd(3114) + Cd(3140)

	Cd(3145) = Cd(48)*Cd(3138) + MGl2*Cd(3144)

	Cd(310) = Cd(308) + MGl2*Cd(309)

	Cd(3154) = Cd(3151) + MGl2*Cd(3153)

	Cd(3157) = 12*Cd(49)*MStop2(2) - 2*Cd(1845)*MStop2(3)

	Cd(3158) = -4*Cd(3118)*MStop2(2) + 2*Cd(2687)*MStop2(3)

	Cd(311) = Cd(44)*Cd(48)**2 - MGl2*Cd(310)

	Cd(3160) = 4*Cd(3159)*MStop2(2) + 2*Cd(49)*Cd(2686)*MStop2(3)

	Cd(3163) = Cd(3157) - 2*MGl2*MStop2(4)

	Cd(3164) = Cd(3158) + MGl2*Cd(3163)

	Cd(3165) = Cd(48)*Cd(3160) - MGl2*Cd(3164)

	end
