* TLhr_abbr1_2.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_2
	implicit none

#include "TLhr.h"


	Cd(13367) = 20*MStop2(2) - 22*MStop2(3)

	Cd(13368) = Cd(1281) + 3*MStop2(4)

	Cd(13369) = 1/2.D0*(Cd(13367) - 4*Cd(13368))

	Cd(13370) = (Cd(13369) - 1/2.D0*(Cr(28)*MStop2(3)))*UCStop(1,3)**2

	Cd(13371) = Cd(1136) - UStop2(1,1)**2

        Cd(13372) = 1/2.D0*
     &   ((8*MStop2(2) - Cd(13371)*MStop2(3))*UCStopC(1,3)**2)/
     &    MStop2(2)**2

        Cd(13373) = 1/Pi**4*
     &   (4*Cr(28)*MStop2(2) + 
     &      2*Cd(13125)*UStop2(1,1)*UStop2(1,2))/MStop2(2)**2

	Cd(13374) = Cd(12700)*UUStop(2,1) + Cd(12701)*UUStop(2,2)

        Cd(13375) = 1/(512.D0*Pi**4)*
     &   (Pi**4*(Cd(12700)*Cd(13373)*MStop2(2)**2) - 
     &      Cd(13374)*MStop2(3)*UStop2(1,1)*UStop2(1,2)*
     &       UUStopC(2,1))/MStop2(2)**2

        Cd(13376) = Cd(12701)*UStop2(1,1)**3*UStop2(1,2) + 
     &   Cd(12700)*UCStop(1,3)*UCStopC(3,1)*UUStop(3,4)*
     &    UUStopC(2,2)

	Cd(13377) = Cd(12700)*Cd(13370) + Cd(13376)*MStop2(3)

        Cd(13379) = 2*Cd(12700)*MStop2(3)*UUStop(1,1) - 
     &   Cd(12701)*MStop2(2)*UUStop(1,2)

        Cd(13380) = Cd(12700)*Cr(91)*MStop2(2)*UStop(2,1) + 
     &   Cd(12701)*MStop2(1)*UStop(1,2)*UUStop(2,2)

        Cd(13381) = Cd(12700)*Cd(12771)*UStop(2,1) + 
     &   Cd(12701)*MStop2(2)*UStop(1,2)*UUStop(2,2)

        Cd(13382) = Cd(13380)*UStop2(1,1)*UStopC(1,2) - 
     &   Cd(13381)*UCStopC(3,1)*UStopC(2,2)

        Cd(13383) = Cd(13379)*MStop2(2)*UStop2(1,1)*UStop2(1,2)*
     &    UUStopC(1,2) + Cd(13382)*MStop2(3)*UUStopC(2,2)

	Cd(13388) = 3 + UStop2(1,1)*UStop2(1,2)

        Cd(13389) = (Cd(13388)*MStop2(3)*UCStopC(1,3) - 
     &     2*Cd(12783)*UStop2(1,2)*UUStopC(3,1))/MStop2(3)

	Cd(13392) = 3*MStop2(4) + MStop2(3)*UStop2(1,2)**2

	Cd(1333) = MStop2(3)**2 + Cd(1281)*MStop2(2)*UStop2(1,2)

	Cd(13394) = 5*UStop2(1,1) - 2*UStop2(1,2)

	Cd(13395) = Cd(13394)*UCStop(1,3) + Cd(11106)*UUStopC(3,2)

	Cd(13397) = MStop2(2)**2 - 2*MStop2(3)*MStop2(4)

	Cd(13398) = Cr(26)*UUStop(3,3) + 2*MStop2(2)*UUStop(3,4)

	Cd(13399) = 3*MStop2(3)**2 + Cd(13397)*UStop2(1,1)*UStop2(1,2)

	Cd(1334) = 2*MStop2(2) + 5*MStop2(3)

        Cd(13404) = Cd(12700)*MStop2(3)*UUStop(2,1) - 
     &   Cd(12701)*MStop2(2)*UUStop(2,2)

        Cd(13405) = Cd(13404)*UStop(1,1) - 
     &   Cd(12701)*MStop2(2)*UStop(2,1)*UUStop(3,2)

        Cd(13406) = Cd(12700)*MStop2(2)*UUStop(2,1) - 
     &   Cd(12701)*MStop2(3)*UUStop(2,2)

        Cd(13407) = Cd(13406)*UStop(1,1) + 
     &   Cd(12701)*MStop2(2)*UStop(2,1)*UUStop(3,2)

	Cd(13408) = Cd(13405) - Cd(13407)

        Cd(13409) = Cd(13408)*UStop2(1,2)*UStopC(1,1) + 
     &   Cd(13191)*MStop2(2)*UUStopC(3,3)

	Cd(13413) = 1/2.D0*Cd(13367) - 6*MStop2(4)

        Cd(131) = MGlpTmSt4(2)*UStop(1,1) - 
     &   2*M_3*MTy*MGlpTmSt2(2)*UStop(1,2)

	Cd(1335) = Cd(1268) - 6*Cd(1280)

	Cd(13414) = MStop2(3)*UStop2(1,1)**2 - MStop2(2)*UStop2(1,2)**2

	Cd(13415) = Cd(13414) + MStop2(3)*UStop2(1,2)**2

	Cd(13416) = Cd(13413) + 2*Cd(13415)

        Cd(13417) = (UCStop(1,3)**2*
     &     (-Cd(13416) + 2*MStop2(2)*UStop2(1,2)**2))/MStop2(2)**2

        Cd(13418) = Cd(13380) - 
     &   Cd(12701)*MStop2(3)*UStop(1,2)*UUStop(2,2)

        Cd(13419) = Cd(13381) - 
     &   Cd(12700)*MStop2(3)*UStop(2,1)*UUStop(3,4)

        Cd(13420) = Cd(13418)*UStop2(1,1)*UStopC(1,2) - 
     &   Cd(13419)*UCStopC(3,1)*UStopC(2,2)

        Cd(13421) = 1/2.D0*
     &   (MStop2(3)*(2*Cd(12700)*MStop2(2)*UCStop(1,3)**2*
     &        UStop2(1,1)*UStop2(1,2) + Cd(13420)*UUStopC(2,2)))

	Cd(13422) = 5 - 2*UStop2(1,2)**2

        Cd(13423) = 4*MStop2(2)*UCStop(3,1)*UStop2(1,1)*
     &    UStopC(1,1) + 2*Cd(13422)*MStop2(3)*UStopC(2,1)

        Cd(13424) = Cd(13423)*UStop(2,2) + 
     &   4*MStop2(3)*UCStopC(1,3)*UStop2(1,1)**2

	Cd(13425) = 1 - UStop2(1,1)*UStop2(1,2)

        Cd(13426) = Cd(13425)*UCStopC(1,3) + 
     &   2*UStop2(1,2)*UUStop(3,3)*UUStopC(3,1)

        Cd(13427) = 1/2.D0*
     &   ((Cd(13424) - 4*Cd(13426)*MStop2(2))*UCStopC(1,3))/
     &    MStop2(2)**2

	Cd(13428) = 5*MStop2(3) - 4*MStop2(1)*UStop2(1,1)*UStop2(1,2)

        Cd(13429) = UStop2(1,2)*
     &   (UStop2(1,2)**2 + UUStop(3,3)*UUStopC(3,3))

        Cd(13432) = Cd(13409) + 
     &   Cd(13374)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(13434) = Cd(12700)**2*Cd(13417) - Cd(12701)**2*Cd(13427)

	Cd(13440) = 1 - Cr(28)

	Cd(13452) = 3*UStop2(1,1)*UStop2(1,2) + Cr(91)*UUStopC(3,3)

	Cd(1339) = 4*Cr(8)*MStop2(1) - 10*MStop2(3)

	Cd(13466) = 1 - Log(MSbot2(1)/Q)

	Cd(13467) = Cd(13466)*dMSfsq1133eps(0) + 2*dMSfsq1133eps(1)

        Cd(13468) = 1/2.D0*
     &    (A0delq(MSbot2(1),Q)*dMSfsq1133eps(-1)) + 
     &   Cd(13467)*MSbot2(1)

	Cd(13469) = MTy*UStop(1,1) + Cd(12701)*UStop(1,2)

	Cd(13470) = MTy*UStopC(1,1) + Cd(12700)*UStopC(1,2)

	Cd(13471) = Cr(13)*dMSfsq1133eps(0) + 2*dMSfsq1133eps(1)

        Cd(13472) = 1/2.D0*(A0delStop(1)*dMSfsq1133eps(-1)) + 
     &   Cd(13471)*MStop2(1)

	Cd(13473) = Cd(12701)*UCStopC(1,3) + MTy*UStop2(1,1)

	Cd(13475) = Cd(13473) + MUE*TB*UCStop(1,3)

	Cd(13478) = Cr(13)*dMSfsq1133eps(-1) + 2*dMSfsq1133eps(0)

	Cd(13479) = At*MUE*TB + MUE2*TB2

	Cd(1342) = 12*MStop2(1)**2 - 4*MStop2(3)**2

        Cd(13485) = (2*MSq2Diff(1,1)**2*MStop2(2)*UStop2(1,1) + 
     &     2*MSq2Diff(1,2)**2*MStop2(1)*UStop2(1,2))/
     &   (MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

        Cd(13486) = (MSq2Diff(1,2)**2*MStop2(1) - 
     &     MSq2Diff(1,1)**2*MStop2(2))/
     &   (MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

	Cd(13487) = Cd(12701)*Cd(13485) + 2*MTy*Cd(13486)*UCStop(1,3)

	Cd(13488) = MTy*Cd(12701)*UCStopC(1,3) + MTy2*UStop2(1,1)

	Cd(13489) = 2 - Cd(12700)*Cd(13487)

	Cd(13490) = -(MTy*Cd(12701)*UCStopC(1,3)) + MTy2*UStop2(1,2)

	Cd(13491) = Cd(13489)*MSq2Diff(1,2)**2 - 2*Cd(13490)*MStop2(2)

	Cd(1343) = MStop2(3) + 9*MStop2(4)

        Cd(13492) = Cd(13491)*MSq2Diff(1,1)**2 - 
     &   2*Cd(13488)*MSq2Diff(1,2)**2*MStop2(1)

        Cd(13493) = (MSq2Diff(1,1)*UStop2(1,1) + 
     &     MSq2Diff(1,2)*UStop2(1,2))/(MSq2Diff(1,1)*MSq2Diff(1,2))

        Cd(13494) = (-MSq2Diff(1,1) + MSq2Diff(1,2))/
     &   (MSq2Diff(1,1)*MSq2Diff(1,2))

	Cd(13495) = Cd(12701)*Cd(13493) + MTy*Cd(13494)*UCStop(1,3)

	Cd(13496) = 1 + Cd(12700)*Cd(13495)

	Cd(13500) = A0delStop(1) + Cd(13466)*Cr(13)*MStop2(1)

	Cd(13501) = Cd(13500)*MSbot2(1) + A0delq(MSbot2(1),Q)*MStop2(1)

	Cd(1344) = 40*MStop2(2)**2 - 3*Cd(1343)*MStop2(3)

	Cd(13502) = A0delq(MSbot2(1),Q) + Cd(13466)*Cr(34)*MSbot2(1)

	Cd(13503) = A0delStop(2)*MSbot2(1) + Cd(13502)*MStop2(2)

	Cd(13504) = MTy*UStop(2,1) + Cd(12701)*UStop(2,2)

	Cd(13505) = MTy*UStopC(2,1) + Cd(12700)*UStopC(2,2)

	Cd(13506) = MTy*UCStop(1,3) + MUEC*TB*UStop2(1,2)

	Cd(13507) = MTy*UCStopC(1,3) + MUE*TB*UStop2(1,2)

	Cd(13508) = MTy*UCStopC(1,3) + Cd(12700)*UStop2(1,2)

	Cd(1345) = 6*MStop2(2) + 11*MStop2(3)

        Cd(13512) = 2*MTy2*MSbot2(1)*UStop2(1,1) - 
     &   2*MSq2Diff(1,1)**2*UStop2(1,2)

	Cd(13515) = MTy2*UStop2(1,1) - MSq2Diff(1,1)*UStop2(1,2)

	Cd(13522) = Cd(13466)*dMSfsq1233eps(0) + 2*dMSfsq1233eps(1)

        Cd(13523) = 1/2.D0*
     &    (A0delq(MSbot2(1),Q)*dMSfsq1233eps(-1)) + 
     &   Cd(13522)*MSbot2(1)

	Cd(13524) = Cr(34)*dMSfsq1233eps(0) + 2*dMSfsq1233eps(1)

        Cd(13525) = 1/2.D0*(A0delStop(2)*dMSfsq1233eps(-1)) + 
     &   Cd(13524)*MStop2(2)

	Cd(13526) = -(Cd(12701)*UCStopC(1,3)) + MTy*UStop2(1,2)

	Cd(13527) = Cd(13526) - MUE*TB*UCStop(1,3)

	Cd(13528) = MUE2*TB2 + 2*MStop2(2)

	Cd(13529) = Cd(13528)*UStop2(1,1) + AtC*Cd(13504)*UStopC(2,2)

	Cd(1347) = 5*MStop2(3) - 3*MStop2(4)

	Cd(13531) = Cr(34)*dMSfsq1233eps(-1) + 2*dMSfsq1233eps(0)

        Cd(13535) = Cd(13469)*Cd(13505)*Conjugate(dY33eps(-1)) + 
     &   Cd(13470)*Cd(13504)*dY33eps(-1)

        Cd(13536) = Cd(13469)*Cd(13505)*Conjugate(dY33eps(0)) + 
     &   Cd(13470)*Cd(13504)*dY33eps(0)

        Cd(13537) = Cd(13469)*Cd(13505)*Conjugate(dY33eps(1)) + 
     &   Cd(13470)*Cd(13504)*dY33eps(1)

	Cd(13538) = Cd(13466)*Cd(13536) + 2*Cd(13537)

        Cd(13539) = 1/2.D0*(A0delq(MSbot2(1),Q)*Cd(13535)) + 
     &   Cd(13538)*MSbot2(1)

	Cd(13540) = MTy*UCStopC(3,1) + Cd(12701)*UCStopC(3,3)

	Cd(1348) = 4*MStop2(2)**2 + Cd(1347)*MStop2(3)

	Cd(13543) = MTy*UCStop(3,1) + Cd(12701)*UCStop(3,4)

	Cd(1349) = 2*MStop2(2) - 9*MStop2(3)

	Cd(13558) = Cd(13466)*dMSfsq1143eps(0) + 2*dMSfsq1143eps(1)

        Cd(13559) = 1/2.D0*
     &    (A0delq(MSbot2(1),Q)*dMSfsq1143eps(-1)) + 
     &   Cd(13558)*MSbot2(1)

	Cd(13560) = MTy*UCStop(1,3) + Cd(12701)*UStop2(1,2)

	Cd(1350) = Cd(1348) + 2*Cd(1349)*MStop2(2)*UStop2(1,2)

	Cd(13561) = -(MTy*UCStop(1,3)) + Cd(12701)*UStop2(1,1)

	Cd(13562) = MTy*Cd(13526) + Cd(12700)*Cd(13561)

	Cd(13563) = -2*MSq2Diff(1,2)**2 + 2*Cd(13562)*MStop2(2)

	Cd(13564) = MTy*Cd(13473) + Cd(12700)*Cd(13560)

        Cd(13565) = (Cd(13563)*MSq2Diff(1,1)**2 + 
     &     2*Cd(13564)*MSq2Diff(1,2)**2*MStop2(1))/
     &   (MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

	Cd(13566) = Cd(13466)*dMSfsq1143eps(-1) + 2*dMSfsq1143eps(0)

	Cd(13570) = Cr(13)*dMSfsq1143eps(0) + 2*dMSfsq1143eps(1)

        Cd(13571) = 1/2.D0*(A0delStop(1)*dMSfsq1143eps(-1)) + 
     &   Cd(13570)*MStop2(1)

	Cd(13572) = Cr(34)*dMSfsq1143eps(0) + 2*dMSfsq1143eps(1)

        Cd(13573) = 1/2.D0*(A0delStop(2)*dMSfsq1143eps(-1)) + 
     &   Cd(13572)*MStop2(2)

	Cd(13574) = Cd(13561)*MSq2Diff(1,1) + Cd(13560)*MSq2Diff(1,2)

	Cd(13578) = -(MTy*UCStopC(1,3)) + Cd(12700)*UStop2(1,1)

	Cd(13579) = Cd(13578)*MSq2Diff(1,1) + Cd(13508)*MSq2Diff(1,2)

	Cd(13583) = 3*MTy*UCStop(1,3) + 2*Cd(12701)*UStop2(1,2)

	Cd(13584) = 3*Cd(12701)*UCStopC(1,3) + 4*MTy*UStop2(1,1)

	Cd(13585) = -3*MTy*UCStop(1,3) + 2*Cd(12701)*UStop2(1,1)

	Cd(13586) = -3*Cd(12701)*UCStopC(1,3) + 4*MTy*UStop2(1,2)

	Cd(13588) = Cd(12700)*Cd(13583) + MTy*Cd(13584)

	Cd(13590) = 3*MTy*UStop(1,1) + 2*Cd(12701)*UStop(1,2)

	Cd(1353) = 10*MStop2(2)*MStop2(3) - Cd(1342)*UStop2(1,1)**3

	Cd(13592) = 3*MTy*UStop(2,1) + 2*Cd(12701)*UStop(2,2)

	Cd(13600) = -(MTy*UCStop(1,3)) + MUEC*TB*UStop2(1,1)

	Cd(13601) = -(MTy*UCStopC(1,3)) + MUE*TB*UStop2(1,1)

	Cd(13611) = 2*At + MTy*UCStop(1,3)

	Cd(13612) = 2*Cd(12700)*UStop(1,1) - MTy*UStop(1,2)

	Cd(13613) = UStop(1,2) + Cr(91)*UStopC(2,1)

	Cd(13614) = MTy*Cd(13613) - 4*Cd(12700)*UStop(1,1)

	Cd(13615) = UStop(1,1) - Cr(91)*UStopC(2,2)

	Cd(13616) = MTy2 + 3*MSq2Diff(1,1)

	Cd(13617) = MTy*MUE*TB*Cr(91) - Cd(13616)*UUStop(3,2)

	Cd(13618) = 3*MTy2 + MSq2Diff(1,1)

	Cd(13619) = -2*MTy*MUE*TB*UCStop(1,3) + Cd(13618)*UStop2(1,2)

        Cd(13620) = 4*MUE*TB*UStop2(1,1) - 
     &   2*MTy*UCStopC(1,3)*UStop2(1,2)

        Cd(13621) = MUEC*TB*Cd(13614) + 
     &   2*Cd(12701)*Cd(13612)*UStop2(1,1)

	Cd(13622) = 2*Cd(13611)*UStop2(1,1) - MTy*Cd(13615)*UStopC(1,2)

	Cd(13624) = 3*MTy2 + MSq2Diff(1,2)

	Cd(13625) = MTy2 + 3*MSq2Diff(1,2)

	Cd(13626) = Cd(13624) + Cd(13625)

	Cd(13627) = Cr(9)*UCStopC(1,3) + Cr(91)*UUStopC(3,1)

	Cd(13628) = Cr(9)*UCStop(1,3) - Cr(91)*UUStopC(3,2)

	Cd(1357) = 8*MStop2(2) - 15*MStop2(3)

	Cd(13630) = 4*Cd(12700)*Cd(12701) - Cd(13626)

	Cd(13633) = MTy*Cr(91) + 2*Cd(12701)*UUStop(3,2)

	Cd(13634) = MUE*TB*UStop(2,1) - MTy*UStop(2,2)

	Cd(13635) = MTy2 - At*MUE*TB

	Cd(13638) = MTy2 + MSq2Diff(1,1)

	Cd(1358) = 8*MStop2(3)**2 + 2*Cd(1357)*MStop2(2)*UStop2(1,2)

	Cd(13641) = 2*MSbot2(1) - MStop2(4)

	Cd(13642) = MTy2 + MSq2Diff(1,2)

        Cd(13643) = Cd(13638)*MSq2Diff(1,2) + 
     &   Cd(13642)*MSq2Diff(1,1)*UStop2(1,2)

	Cd(13647) = -(MTy*UCStopC(1,3)) + 2*MUE*TB*UStop2(1,1)

	Cd(13648) = UCStopC(1,3)*UStop2(1,1) - Cr(91)*UUStopC(3,1)

	Cd(1359) = 1 - 2*UStop2(1,1)*UStop2(1,2)

	Cd(13649) = MTy*Cd(13648) + Cd(13647)*UStop2(1,2)

	Cd(13650) = MTy*UStopC(2,1) + 2*Cd(12700)*UStopC(2,2)

	Cd(13651) = MTy*Cd(13648) + Cd(13650)*UStop(2,2)*UStop2(1,2)

	Cd(13652) = MTy*Cr(91) + 2*MUEC*TB*UUStop(3,2)

	Cd(13653) = MTy*Cr(9)*UCStop(1,3) - Cd(13652)*UUStopC(3,2)

	Cd(13654) = 6*UUStop(3,3) + UUStop(3,4)

	Cd(13655) = Cd(13654)*MSq2Diff(1,2) + MTy2*UUStop(3,4)

	Cd(13660) = 2 - Log(MStop2(1)/Q)

	Cd(13661) = 7 - 2*UStop2(1,1)

	Cd(13662) = 5*MSbot2(1) - Cd(13661)*MStop2(1)

	Cd(13663) = Cd(12700)*UStop(1,1) - MTy*UStop(1,2)

        Cd(13664) = 5*MTy*MSq2Diff(1,1) + 
     &   2*Cd(13663)*MStop2(1)*UStopC(1,2)

	Cd(13665) = 5*MSq2Diff(1,1) - 2*MStop2(1)*UStop2(1,2)

        Cd(13666) = MUE*TB*Cd(13665)*UCStop(1,3) + 
     &   2*Cd(13469)*MStop2(1)*UStop2(1,1)*UStopC(1,1)

	Cd(13667) = MUE2*TB2*UStop2(1,2) + AtC*Cd(13469)*UStopC(1,2)

	Cd(13669) = 3 - Log(MStop2(1)/Q)

	Cd(13670) = 2*A0delStop(1) + Cd(13669)**2*MStop2(1)

	Cd(13671) = 4 - UStop2(1,1)

	Cd(13672) = -6*MSbot2(1) + 2*Cd(13671)*MStop2(1)

        Cd(13673) = 6*MTy*MSq2Diff(1,1) + 
     &   2*Cd(13663)*MStop2(1)*UStopC(1,2)

	Cd(13674) = 3*MSq2Diff(1,1) - MStop2(1)*UStop2(1,2)

        Cd(13675) = 2*(MUE*TB*Cd(13674)*UCStop(1,3) + 
     &     Cd(13469)*MStop2(1)*UStop2(1,1)*UStopC(1,1))

	Cd(13679) = 2 - Log(MStop2(2)/Q)

	Cd(13680) = 5*MSq2Diff(1,2) - 2*MStop2(2)*UStop2(1,1)

	Cd(13681) = Cd(12700)*UStop(2,1) - MTy*UStop(2,2)

        Cd(13682) = 5*MTy*MSq2Diff(1,2) + 
     &   2*Cd(13681)*MStop2(2)*UStopC(2,2)

        Cd(13683) = MUE*TB*Cd(13680)*UCStop(1,3) - 
     &   2*Cd(13504)*MStop2(2)*UStop2(1,2)*UStopC(2,1)

	Cd(13684) = MUE2*TB2*UStop2(1,1) + AtC*Cd(13504)*UStopC(2,2)

	Cd(13686) = 3 - Log(MStop2(2)/Q)

	Cd(13687) = 2*A0delStop(2) + Cd(13686)**2*MStop2(2)

	Cd(13688) = 6*MSq2Diff(1,2) - 2*MStop2(2)*UStop2(1,1)

        Cd(13689) = 6*MTy*MSq2Diff(1,2) + 
     &   2*Cd(13681)*MStop2(2)*UStopC(2,2)

	Cd(13690) = 3*MSq2Diff(1,2) - MStop2(2)*UStop2(1,1)

        Cd(13691) = 2*MUE*TB*Cd(13690)*UCStop(1,3) - 
     &   2*Cd(13504)*MStop2(2)*UStop2(1,2)*UStopC(2,1)

	Cd(13696) = MGl2 + MSbot2(1)

	Cd(13698) = MTy2 + At*MUE*TB

	Cd(134) = 2*MGlpTmSt4(2) - Cd(131)*UStopC(1,1)

	Cd(13711) = MSq2Diff(1,2)*(MGlmT2 - 2*MStop2(2))

	Cd(13719) = 2*MSbot2(2) - 2*MSbot2(3) - MStop2(2)

	Cd(13721) = 2*Cr(56)*MSbot2(1) - Cr(49)*MStop2(2)

	Cd(13726) = AtC*MUEC*TB + MUE2*TB2

	Cd(13727) = MUE*TB*Cr(56) + AtC*UStop2(1,1)

	Cd(13728) = At*Cd(13727) + Cd(13726)*UStop2(1,1)

	Cd(13730) = MUE*TB*Cr(49) + AtC*UStop2(1,1)

	Cd(13731) = At*Cd(13730) + Cd(13726)*UStop2(1,1)

	Cd(13737) = UStop2(1,1) + UStop2(1,2)**2

	Cd(13739) = Cd(13737)*MSbot2(1) - MStop2(2)*UStop2(1,2)

	Cd(13752) = MUEC*TB*MSbot2(1) + At*MSq2Diff(1,2)

	Cd(13765) = MTy2 - MSbot2(1)

	Cd(13767) = MSbot2(1) + MStop2(1)

	Cd(13780) = -MSq2Diff(1,1)

	Cd(13781) = -(MGl2*MSq2Diff(1,1)) + Cd(13780)*MStop2(1)

	Cd(13792) = At*Cd(12700) - MSbot2(1)

	Cd(1374) = 3 + 2*UStop2(1,1)*UStop2(1,2)

	Cd(13805) = -MSbot2(1) + 2*MStop2(1)

        Cd(13806) = Cd(13805)*MSq2Diff(1,2) - 
     &   MSbot2(1)*MSq2Diff(1,1)*UStop2(1,1)

	Cd(1375) = 1 + Cd(1374)*UStop2(1,2)

	Cd(13808) = MSbot2(2) - MSbot2(3) - 2*MStop2(2)

	Cd(13809) = Cr(56)*MSbot2(1) - 2*MStop2(2)

	Cd(13810) = MTy*Cd(13809) - Cd(12700)*MSbot2(1)*UCStop(1,3)

	Cd(13812) = Cd(12700)*UStop(2,1) + MTy*UStop(2,2)

        Cd(13813) = MTy*Cd(13805)*MSq2Diff(1,2) - 
     &   Cd(13812)*MSbot2(1)*MSq2Diff(1,1)*UStopC(2,2)

        Cd(13814) = Cd(13813)*UStop(1,2) - 
     &   MTy*Cr(91)*MSbot2(1)*MSq2Diff(1,1)*UStopC(2,1)

	Cd(13815) = Cd(12700)*UStop(1,1) + MTy*UStop(1,2)

	Cd(13816) = Cd(9948)*MStop2(1) - MSbot2(1)*UStop2(1,2)

	Cd(13817) = MStop2(2) + MStop2(1)*UStop2(1,2)

	Cd(1376) = Cd(1375)*MStop2(2) - 4*MStop2(1)*UStop2(1,2)

        Cd(13819) = At*AtC*Cd(13816)*UStop2(1,1) + 
     &   Cd(13469)*Cd(13815)*MSq2Diff(1,2)*UUStopC(1,3)

        Cd(13821) = MTy*UStop2(1,1)*UStopC(1,1) + 
     &   Cd(12700)*UStop2(1,2)*UStopC(1,2)

	Cd(1377) = 4*MTy2*Cd(1376) + MStop2(2)*MStop2(3)

	Cd(13828) = 2*MGl2 + MSbot2(1)

	Cd(13839) = Cr(87)*MSq2Diff(1,2) - MUE2*TB2*UStop2(1,1)**2

	Cd(13840) = MSq2Diff(1,2)*(MGl2 + MStop2(2))

	Cd(13841) = At*AtC + MUE2*TB2

	Cd(13847) = MSq2Diff(1,2)

	Cd(13848) = MGl2*MSq2Diff(1,2) + Cd(13847)*MStop2(2)

	Cd(13854) = MSq2Diff(1,1)*(MTy2 + 2*MStop2(1))

	Cd(13857) = MSq2Diff(1,1)*(MGl2 + MStop2(1))

	Cd(13863) = MSbot2(1)**2 - At*Cd(12700)*MSq2Diff(1,1)

	Cd(13866) = -(MTy2*MSq2Diff(1,1)) + 2*Cd(13780)*MStop2(1)

	Cd(13871) = -(MUEC*TB*MSq2Diff(1,1)) + At*MStop2(1)

	Cd(13872) = (Cd(13871) - At*MSbot2(1))*MSq2Diff(1,2)

        Cd(13873) = -(MTy*Cd(13872)*UCStopC(1,3)) + 
     &   MTy2*MSq2Diff(1,1)*MSq2Diff(1,2)*UStop2(1,1)

	Cd(1382) = 4*Cr(75)**2*MStop2(1) + 10*MStop2(3)

	Cd(13893) = 1 - 2*Cd(9969)*UStop2(1,1)

	Cd(13901) = MSq2Diff(1,2) + 2*Cd(13841)*UStop2(1,1)

	Cd(13912) = MGlpT2 - 2*MTy*Cr(5)

	Cd(1387) = 1 + 3*UStop2(1,1)

        Cd(13934) = MTy*MTy2*MUE*TB*UStop(1,1) - 
     &   MSq2Diff(1,2)*MStop2(2)*UStop(1,2)

        Cd(13935) = MTy*MTy2*MUE*TB*UStop(1,1) + 
     &   MSq2Diff(1,2)*MStop2(2)*UStop(1,2)

	Cd(1388) = 3 - 7*UStop2(1,1)

	Cd(13936) = AtC - Cd(12700)*UStop2(1,1)

	Cd(13942) = At*Cd(12700)*Cr(56) + Cd(13726)*UStop2(1,1)

	Cd(13944) = 2*Cr(56)*MStop2(2) + MUE2*TB2*UStop2(1,1)

	Cd(1390) = 2*MStop2(1)**2 + MStop2(2)*MStop2(3)

	Cd(13955) = MGl2*MSq2Diff(1,2) + 2*Cd(13847)*MStop2(2)

	Cd(13962) = MGl2 + At*MUE*TB

        Cd(13967) = 2*MTy*MUE*TB*MSbot2(1)*UStop(1,1) - 
     &   2*MSq2Diff(1,2)*MStop2(2)*UStop(1,2)

	Cd(13969) = MTy*MUE*TB*UStop(1,1) - MSq2Diff(1,2)*UStop(1,2)

	Cd(13972) = At*AtC*Cr(56) + MUE2*TB2*UStop2(1,1)

	Cd(13996) = MUE*TB*UCStop(1,3) + Cd(12701)*UCStopC(1,3)

	Cd(1395) = 6 - 11*UStop2(1,2)

	Cd(14005) = -2*MSbot2(1) + MStop2(1)

	Cd(14008) = MStop2(3) - MTy*Cd(13469)*UStopC(1,1)

	Cd(14013) = MStop2(3) - MTy*Cd(13812)*UStopC(2,2)

	Cd(14014) = Cd(14013)*UStop(1,2) - MTy2*Cr(91)*UStopC(2,1)

	Cd(1396) = -(Cd(1395)*MStop2(2)) + 4*MStop2(1)*UStop2(1,2)

	Cd(14015) = Cd(13780)*MStop2(1) + MSq2Diff(1,1)*MStop2(2)

	Cd(14016) = 2 + Cd(9948)

	Cd(14018) = 5 - Cd(9948)

	Cd(1397) = 22*MStop2(3) - 2*Cd(1396)*UStop2(1,2)

        Cd(14031) = MUE*TB*UCStop(1,3)*UStop(1,1) - 
     &   At*UStop(1,2)*UStop2(1,1)

	Cd(14032) = MStop2(3) - MStop2(1)*UStop2(1,1)

        Cd(14033) = MUE*TB*MStop2(1)*UCStop(1,3)*UStop(1,1) + 
     &   At*Cd(14032)*UStop(1,2)

	Cd(1398) = 3*Cd(1312) + 2*MStop2(3)

	Cd(14037) = MUE*TB*UStop(1,1) + MTy*UStop(1,2)

	Cd(14040) = MUE*TB*UUStop(1,1) + At*UUStop(1,2)

	Cd(14047) = MUE2*TB2 - MSq2Diff(1,1)

        Cd(14054) = MTy*MUE*Cd(9948)*UCStop(1,3) + 
     &   MUE2*TB*UStop2(1,2)**2

	Cd(14056) = Cr(54)*MSq2Diff(1,1) - MUE2*TB2*UStop2(1,2)

	Cd(14057) = Cd(14056)*MStop2(1) - MTy2*MSq2Diff(1,1)*UStop2(1,2)

	Cd(14058) = TB*Cd(14057) - At*MUE*TB2*MStop2(1)*UStop2(1,2)

	Cd(14059) = Cr(85)*MSq2Diff(1,1) - MUE2*TB2*UStop2(1,2)

	Cd(14061) = 2*Cd(9948) - UStop2(1,2)

	Cd(14062) = 2*MTy*Cd(14061)*UCStop(1,3) + At*UStop2(1,2)**2

	Cd(14069) = 2*MTy2 + At*MUE*TB

	Cd(14070) = At*MUE*TB - MGlpTmSt2(1)

	Cd(14076) = Cd(9950) - MStop2(1)**2

	Cd(14081) = At*MUE + MUE2*TB

	Cd(1403) = 6*MStop2(2)*MStop2(3) - 4*Cd(1327)*UStop2(1,2)

	Cd(1404) = MStop2(2) + 11*MStop2(3)

	Cd(14096) = MSq2Diff(1,2) + MTy*Cd(13812)*UStopC(2,2)

	Cd(14097) = Cd(14096)*UStop(1,2) + MTy2*Cr(91)*UStopC(2,1)

	Cd(14103) = 2*MTy*UStop(2,1) + Cd(12701)*UStop(2,2)

	Cd(138) = MGl2 - MStop2(2)

	Cd(1405) = MStop2(3) - MStop2(4)*UStop2(1,2)

        Cd(14104) = Cd(14103)*MSq2Diff(1,2)*UStopC(2,1) + 
     &   MTy*Cd(13504)*Cd(13812)*UUStopC(2,3)

	Cd(14107) = MSq2Diff(1,2) + MTy*Cd(13469)*UStopC(1,1)

	Cd(14108) = MSq2Diff(1,2)*UStop(1,1) + MTy*Cd(13633)*UStopC(2,2)

	Cd(1406) = MStop2(2) + 6*Cd(1405)*UStop2(1,2)

	Cd(1407) = Cd(1404)*MStop2(1) - Cd(1406)*MStop2(2)

	Cd(14125) = AtC*TB + MUE*TB2

	Cd(14131) = MTy*TB*UCStopC(1,3) + MUE*TB2*UStop2(1,2)**2

	Cd(14132) = TB*MStop2(1) + 2*M_3*MUE*TB2*UStop2(1,2)

        Cd(14133) = Cd(14132)*UCStopC(1,3) + 
     &   2*M_3C*MUE*TB2*UCStop(1,3)*UStop2(1,2)

	Cd(14134) = -2*MTy*Cr(5) + MGl2*UStop2(1,2)

	Cd(14147) = At*TB + MUEC*TB2

	Cd(14148) = Cd(14125)*UCStop(1,3) + Cd(14147)*UCStopC(1,3)

	Cd(14150) = MStop2(1) + MSq2Diff(1,1)*UStop2(1,2)

        Cd(14151) = Cd(9950)*Cd(14150) + 
     &   MGlpTmSt2(1)*MSq2Diff(1,1)*MStop2(1)*UStop2(1,2)

        Cd(14160) = 2*MSq2Diff(1,1)*UStop2(1,1) + 
     &   Cd(13841)*UStop2(1,2)**2

	Cd(14163) = MSq2Diff(1,1) + Cd(13841)*UStop2(1,2)

	Cd(14168) = Cr(75)*MSq2Diff(1,1) - MStop2(3)*UStop2(1,2)

	Cd(14171) = At - 2*M_3*UStop2(1,1)

	Cd(14175) = At*Cd(9981) - 2*Cd(10086)

	Cd(14177) = AtC*Cd(10169) + At*Cd(10176)

	Cd(14178) = AtC*Cd(9981) - 2*Cd(10039)

	Cd(1413) = 1 - 8*UStop2(1,2)**2

        Cd(14186) = Cd(9969)*Cd(13841)*MSbot2(1)*UStop2(1,1) - 
     &   MTy2**2*UStop2(1,2)

	Cd(1414) = Cd(1413)*MStop2(2) + 5*MStop2(3)*UStop2(1,2)

	Cd(14200) = MGl2*Cd(9948) - MStop2(2)*UStop2(1,2)

	Cd(1415) = MStop2(3) - 4*MStop2(4)

	Cd(14212) = 2*MSq2Diff(1,2) + MUE2*TB2*UStop2(1,1)

	Cd(1416) = 12*MStop2(3) - Cd(1415)*UStop2(1,2)**2

        Cd(14221) = Cr(28)*MSq2Diff(1,2) + 
     &   MSq2Diff(1,1)*UStop2(1,1)*UStop2(1,2)

	Cd(14222) = MSq2Diff(1,1) - 2*MSq2Diff(1,2)

	Cd(14223) = M_3*Cd(14221) + M_3C*Cd(14222)*UCStop(1,3)**2

	Cd(14228) = M_3C*Cd(14221) + M_3*Cd(14222)*UCStopC(1,3)**2

	Cd(14235) = MSbot2(1) - MSq2Diff(1,2)

	Cd(14236) = -(Cr(75)*MSq2Diff(1,2)) + MSq2Diff(1,1)*UStop2(1,1)

	Cd(14237) = -(Cr(75)*MSq2Diff(1,2)**2) + Cd(14236)*MStop2(3)

	Cd(14240) = At - 2*M_3*UStop2(1,2)

        Cd(14241) = -(Cd(14240)*UCStopC(1,3)) + 
     &   2*M_3C*UCStop(1,3)*UStop2(1,2)

	Cd(14242) = Cr(75)*MSq2Diff(1,2) + MSq2Diff(1,1)*UStop2(1,2)

	Cd(14243) = 2*M_3*Cd(14242) - At*MSq2Diff(1,1)

        Cd(14244) = 2*M_3C*Cd(14242)*UCStop(1,3) + 
     &   Cd(14243)*UCStopC(1,3)

	Cd(14247) = MSq2Diff(1,2) + MUE2*TB2*UStop2(1,1)

        Cd(14248) = Cd(14247)*MSq2Diff(1,1) - 
     &   MUE2*TB2*Cr(75)*MSq2Diff(1,2)

	Cd(14249) = 2*M_3*Cd(14248) - At*MSq2Diff(1,1)*MSq2Diff(1,2)

        Cd(14250) = 2*M_3C*Cd(14248)*UCStop(1,3) + 
     &   Cd(14249)*UCStopC(1,3)

	Cd(14251) = MUE2*TB2 - MSq2Diff(1,2)

	Cd(14252) = At*MSq2Diff(1,2) - 2*M_3*Cd(14251)*UStop2(1,1)

        Cd(14253) = -(Cd(14252)*UCStopC(1,3)) + 
     &   2*M_3C*Cd(14251)*UCStop(1,3)*UStop2(1,1)

	Cd(14255) = MSq2Diff(1,2)*UCStop(1,3) - 2*At*Cr(5)*UStop2(1,1)

        Cd(14256) = -(Cd(14255)*MSq2Diff(1,1)) - 
     &   2*At*Cr(5)*Cr(75)*MSq2Diff(1,2)

	Cd(1421) = 8*MStop2(2) - 5*MStop2(3)

	Cd(14265) = Cd(13638)*UStop(1,1) + MTy*Cd(12701)*UStop(1,2)

	Cd(14266) = Cd(13638)*UStopC(1,1) + MTy*Cd(12700)*UStopC(1,2)

	Cd(14268) = MTy2*UCStop(1,3) + MTy*Cd(12701)*UStop2(1,2)

	Cd(14269) = -(MTy2*UCStop(1,3)) + MTy*Cd(12701)*UStop2(1,1)

	Cd(14270) = MTy*Cd(13526) + MSq2Diff(1,2)

	Cd(14271) = MTy*Cd(12700)*Cd(14269) + MTy2*Cd(14270)

	Cd(14272) = MTy2*Cd(13473) + Cd(12700)*Cd(14268)

	Cd(1422) = 5*MStop2(3) - 8*MStop2(2)*UStop2(1,2)

        Cd(14273) = Cd(14271)*MSq2Diff(1,1) + 
     &   MTy*Cd(14272)*MSq2Diff(1,2)

	Cd(14277) = MTy2 + 2*MSq2Diff(1,1)

	Cd(14278) = MTy*Cd(14277)*UStop(1,1) + At*Cd(13638)*UStop(1,2)

	Cd(14288) = MTy2 + 2*MSq2Diff(1,2)

	Cd(14289) = MTy*Cd(14288)*UStop(2,1) + At*Cd(13642)*UStop(2,2)

	Cd(14290) = MTy2 - 2*MSbot2(1)

	Cd(14291) = Cd(14290)*MSbot2(1) - 2*Cd(13765)*MStop2(2)

	Cd(14292) = MTy2*Cd(13808) - MSbot2(1)*MSq2Diff(1,2)

	Cd(14299) = Cd(13642)*UStop(2,1) + MTy*Cd(12701)*UStop(2,2)

	Cd(14300) = Cd(13642)*UStopC(2,1) + MTy*Cd(12700)*UStopC(2,2)

	Cd(14301) = MSbot2(2) - MSbot2(3) - 5*MStop2(2)

	Cd(14303) = MSbot2(2) - MSbot2(3) - 3*MStop2(2)

	Cd(14305) = MTy*MUE*TB*UCStop(1,3) + Cd(13479)*UStop2(1,2)

        Cd(14308) = 6*MTy*MSq2Diff(1,2) - 
     &   2*Cd(13681)*MSbot2(1)*UStopC(2,2)

        Cd(14309) = 5*MTy*MSq2Diff(1,2) - 
     &   2*Cd(13681)*MSbot2(1)*UStopC(2,2)

	Cd(14311) = UCStop(3,1)*UStop(2,2) + UStop(1,2)*UStop2(1,2)

	Cd(14312) = -2*MTy*Cd(14311)*MSbot2(1) + Cd(14308)*UStop(1,2)

	Cd(14313) = At*AtC - MTy2

        Cd(14314) = 3*MTy2*MSq2Diff(1,2) + 
     &   Cd(14313)*MSbot2(1)*UStop2(1,2)

	Cd(14315) = 3*MSq2Diff(1,2) + MSbot2(1)*UStop2(1,1)

        Cd(14316) = Cd(12700)*UCStop(1,3)*UStop(1,1) - 
     &   Cd(11166)*UStopC(2,1)

	Cd(14319) = 2*MTy*Cd(14311)*MSbot2(1) - Cd(14309)*UStop(1,2)

        Cd(14320) = 5*MTy2*MSq2Diff(1,2) + 
     &   2*Cd(14313)*MSbot2(1)*UStop2(1,2)

	Cd(14321) = 5*MSq2Diff(1,2) + 2*MSbot2(1)*UStop2(1,1)

	Cd(14328) = Cd(12700)*Cd(13561) + MSq2Diff(1,2)

	Cd(14331) = MTy2*UStop(1,1) + MTy*Cd(12701)*UStop(1,2)

	Cd(14337) = 3*MUE*TB*UStop(1,1) - MTy*UStop(1,2)

        Cd(14341) = MTy*UStop2(1,1)*UStopC(1,1) - 
     &   Cd(12700)*UStop2(1,2)*UStopC(1,2)

	Cd(14345) = MTy*Cd(12700)*UCStop(1,3) + MTy2*UStop2(1,1)

	Cd(14346) = -(Cd(13642)*MSq2Diff(1,1)) + Cd(14345)*MStop2(3)

	Cd(14347) = Cd(9948)*Cd(12700) + MTy*UCStopC(1,3)

	Cd(14348) = Cd(13578)*MSq2Diff(1,1) + Cd(14347)*MSq2Diff(1,2)

        Cd(14349) = Cd(14270)*MSq2Diff(1,1) + 
     &   MTy*Cd(13473)*MSq2Diff(1,2)

	Cd(14350) = MGlpT2 - MSbot2(1)

	Cd(1430) = MStop2(2) - 6*MStop2(3)

	Cd(1432) = 2*MStop2(2) - 5*MStop2(3)

	Cd(1433) = 2 - 6*UStop2(1,2)**2

	Cd(1434) = -(Cd(1433)*MStop2(2)) + Cd(1432)*UStop2(1,2)

	Cd(141) = MGl2 + MStop2(2)

	Cd(1435) = 12*MStop2(3) + Cd(1434)*UStop2(1,2)

	Cd(1437) = 2*MStop2(2) + Cd(1421)*UStop2(1,1)**3

	Cd(1439) = 3*MStop2(1) - 2*Cr(53)*MStop2(2)

	Cd(1440) = Cd(1421)*UStop2(1,1) - 4*MStop2(2)*UStop2(1,2)

	Cd(1442) = 4*MStop2(2) - 5*MStop2(3)

	Cd(142) = 2*Cd(141) + MGlpTmSt2(2)*UStop2(1,2)

	Cd(1451) = MStop2(2) - Cd(1421)*UStop2(1,2)**2

        Cd(143) = 2*(M_3*Cd(142)*UCStopC(1,3) + 
     &     M_3C*MGlpTmSt2(2)*UCStop(1,3)*UStop2(1,1))

	Cd(1462) = 2 - Cr(75)**2

	Cd(1467) = 2*MStop2(2) - 5*MStop2(4)*UStop2(1,1)

	Cd(1468) = 4*MStop2(2) - 5*MStop2(3)

	Cd(1469) = 4*Cd(1439) + MStop2(3)*UStop2(1,1)**2

	Cd(1473) = 4*MStop2(2) - 5*MStop2(4)*UStop2(1,2)

	Cd(1474) = 12*MStop2(3) + Cd(1473)*UStop2(1,2)

	Cd(1475) = 6*MStop2(3) + Cr(85)*MStop2(2)*UStop2(1,2)

	Cd(1480) = 2*MStop2(2) - Cd(1415)*UStop2(1,2)

	Cd(1481) = 2*MStop2(2) + 3*MStop2(3)

	Cd(1482) = 2*MStop2(2) + Cr(85)*MStop2(3)*UStop2(1,2)

	Cd(1483) = Cd(1481)*MStop2(1) - Cd(1482)*MStop2(2)

	Cd(13) = MUE*UCStop(1,3) + MUEC*UCStopC(1,3)

	Cd(1493) = MGlpT2 + MStop2(1)

	Cd(1496) = MGl2 + 2*MTy2

	Cd(1508) = MGlpT2 + MStop2(2)

	Cd(1511) = (-2*MTy2 + Cd(1272))*MStop2(3)

	Cd(1523) = 2*UStop2(1,1)**2 + UStop2(1,2)

	Cd(1532) = 8*MTy2 - 3*MStop2(1)

	Cd(1533) = 11*MTy2 - 3*MStop2(1)

	Cd(1535) = 8*MTy2 - 3*MStop2(3)

	Cd(1541) = 12*MTy2 + MStop2(1)

	Cd(1542) = 7*MTy2 + 3*MStop2(1)

	Cd(152) = MGlmT2 - MStop2(1)

	Cd(1545) = 21*MTy2 + 10*MStop2(1)

	Cd(1546) = 16*MTy2 - 9*MStop2(1)

	Cd(153) = 2*MGl2 + MTy2

	Cd(1556) = 16*MStop2(2) - 3*MStop2(3)

	Cd(1557) = 16*MStop2(2) - 9*MStop2(3)

	Cd(1559) = MTy2 + 3*MStop2(2)

	Cd(1569) = 6*MStop2(2) + MStop2(3)

	Cd(1570) = 9*MStop2(2) - MStop2(3)

	Cd(1572) = -4*MGl2 + Cd(1556)

	Cd(155) = MTy2 + 2*MStop2(1)

	Cd(1583) = 8*MStop2(1) + MStop2(3)

	Cd(1584) = 8*MTy2 + 3*MStop2(3)

	Cd(1585) = 2*MTy2*Cd(1584) + 6*Cd(1583)*MStop2(1)

	Cd(1592) = MTy2 + 3*MStop2(1)

	Cd(1594) = 8*MStop2(2) - 5*MStop2(3)

	Cd(1595) = 3*MStop2(2) - 4*MStop2(3)

	Cd(1596) = 12*MStop2(1)**2 + Cd(1595)*MStop2(3)

	Cd(1610) = 16*MStop2(2) - 13*MStop2(3)

	Cd(1611) = 19*MTy2 - 7*MStop2(1)

	Cd(1612) = Cd(1611)*MStop2(1) - 3*Cd(1592)*MStop2(2)

	Cd(1622) = -4*MGl2 + Cd(1610)

	Cd(160) = 2*MTy2 + 3*MStop2(1)

	Cd(1634) = 4*MStop2(2) - MStop2(3)

	Cd(161) = -MGl2 + Cd(155)

	Cd(1635) = MTy2*Cd(1442) - Cd(1634)*MStop2(1)

	Cd(1636) = 3*(MTy2 - MStop2(2))

	Cd(1637) = 5*MTy2 - 3*MStop2(4)

	Cd(1647) = 2*MTy2 + MStop2(3)

	Cd(1657) = 4*MStop2(1) + MStop2(3)

	Cd(1658) = -(MTy2*Cd(1141)) + Cd(1657)*MStop2(2)

	Cd(1659) = 3*(MTy2 - MStop2(1))

	Cd(164) = MGlmT2 + MStop2(1)

	Cd(1667) = 2*MTy2 - MStop2(3)

	Cd(1668) = MTy2*Cd(1667) + Cd(1226)*MStop2(2)

	Cd(165) = -2*MGl2 + Cd(160)

	Cd(1684) = 3*MTy2 + 5*MStop2(1)

	Cd(1702) = 8*MStop2(2) - 9*MStop2(3)

	Cd(169) = 3 - UStop2(1,2)

	Cd(1715) = 21*MTy2 + 4*MStop2(1)

	Cd(1718) = 8*MTy2 + MStop2(3)

	Cd(1719) = MTy2*Cd(1718) - Cd(1122)*MStop2(2)

	Cd(1732) = 5*MStop2(2) - 8*MStop2(3)

	Cd(1735) = MTy2 - 2*MStop2(1)

	Cd(1738) = 7*MTy2 + 6*MStop2(1)

	Cd(1739) = 4*MTy2 - 3*MStop2(1)

	Cd(172) = 2 - UStop2(1,2)

	Cd(1747) = 5*MStop2(2) + 8*MStop2(3)

	Cd(1749) = 16*MStop2(2) - MStop2(3)

	Cd(1750) = 5*MStop2(2)**2 + Cd(1749)*MStop2(3)

	Cd(1751) = MStop2(2) - 8*MStop2(3)

	Cd(1752) = MStop2(2)**2 + MStop2(3)**2

	Cd(173) = -(Cd(172)*MStop2(1)) + MStop2(2)

	Cd(1755) = MStop2(3) + 2*MStop2(4)

	Cd(1756) = MStop2(2)**2 + 2*Cd(1755)*MStop2(3)

	Cd(1757) = 3*MStop2(3) + 4*MStop2(4)

	Cd(1758) = 5*MStop2(2)**2 + 3*Cd(1757)*MStop2(3)

	Cd(1759) = 5*MStop2(3) + 8*MStop2(4)

	Cd(1760) = 5*MStop2(2)**2 + 2*Cd(1759)*MStop2(3)

	Cd(1761) = 7*MStop2(3) + 8*MStop2(4)

	Cd(1763) = MStop2(2) + 16*MStop2(3)

	Cd(1764) = 7*MStop2(3) + 4*MStop2(4)

	Cd(1765) = 7*MStop2(2)**2 + 2*Cd(1764)*MStop2(3)

	Cd(1768) = 8*MStop2(2) - 2*MStop2(3)

	Cd(1789) = -MGl2 + Cd(1684)

	Cd(1795) = 20*MStop2(1) + MStop2(3)

	Cd(1797) = 3*MGl2 + 5*Cr(1)

	Cd(1798) = 14*MStop2(3) - 15*MStop2(4)

	Cd(1799) = 5*MStop2(2)**2 - Cd(1798)*MStop2(3)

	Cd(1800) = 9*MStop2(3) - 5*MStop2(4)

	Cd(1801) = 5*MStop2(2)**2 - Cd(1800)*MStop2(3)

	Cd(1802) = 3*MStop2(3) - 4*MStop2(4)

	Cd(1803) = MStop2(2)**2 - Cd(1802)*MStop2(3)

	Cd(1804) = 6*MStop2(3) - 5*MStop2(4)

	Cd(1805) = 5*MStop2(2)**2 - 2*Cd(1804)*MStop2(3)

	Cd(1806) = MStop2(2) - 11*MStop2(3)

	Cd(1807) = MStop2(1)

	Cd(1809) = MStop2(2) + 19*MStop2(3)

	Cd(1810) = 16*MStop2(3) - 3*MStop2(4)

	Cd(1811) = 7*MStop2(2)**2 - Cd(1810)*MStop2(3)

	Cd(1814) = 8*MStop2(2) - 6*MStop2(3)

	Cd(179) = MGlmT2 - MStop2(2)

	Cd(1824) = 20*MStop2(2) - 17*MStop2(3)

	Cd(1832) = 23*MTy2 - 11*MStop2(1)

	Cd(1833) = Cd(1832)*MStop2(1) - 3*Cd(1592)*MStop2(2)

	Cd(181) = MTy2 + 2*MStop2(2)

	Cd(182) = -MGl2 + Cd(181)

	Cd(1845) = 3*MTy2 + 5*MStop2(2)

	Cd(1849) = 17*MTy2 - 6*MStop2(1)

	Cd(1852) = 5*MTy2 - 2*MStop2(1)

	Cd(1855) = -MGl2 + Cd(1845)

	Cd(1864) = 4*MTy2 + 7*MStop2(1)

	Cd(1866) = 2*MTy2**2 - MStop2(1)**2

	Cd(17) = Cr(13) + Cr(34)

	Cd(186) = 2*(MTy2*MStop2(1) + Cr(26)*MStop2(2))

	Cd(1888) = 7*MGl2 + 3*MTy2

	Cd(1898) = 3*MStop2(2) - MStop2(3)

	Cd(188) = MStop2(2) - MStop2(1)*UStop2(1,2)

	Cd(1908) = 2*M_3C*MTy*UCStop(1,3) + MGlpTmSt2(2)*UStop2(1,2)

	Cd(189) = MStop2(2)*UStop2(1,1) + 2*MStop2(3)*UStop2(1,2)

	Cd(1923) = 9*MGl2 + 10*MTy2

	Cd(191) = 2*Cd(189)*MStop2(2) - Cd(186)*UStop2(1,1)

	Cd(1952) = MTy2 - MStop2(3)

	Cd(1953) = 7*MGl2 - 9*MStop2(1)

	Cd(1954) = MGl2 - 3*MStop2(1)

	Cd(1955) = Cd(1953)*MStop2(1) + Cd(1954)*MStop2(2)

	Cd(1962) = 12*MTy2 - MStop2(4)

	Cd(1974) = 17*MStop2(2) - 12*MStop2(3)

	Cd(1975) = 3*(4*MStop2(2) - 7*MStop2(3))

	Cd(2) = dMSfsq1233eps(-1) - 2*dMSfsq1233eps(0)

	Cd(196) = 2*MStop2(2) + Cd(179)*UStop2(1,2)

	Cd(1998) = 20*MTy2 - 9*MStop2(1)

	Cd(1999) = 17*MTy2 - 10*MStop2(1)

	Cd(2002) = MGl2 - 3*MStop2(2)

	Cd(198) = 2*MTy2 + MStop2(4)

	Cd(199) = -MGl2 + Cd(198)

	Cd(2020) = 2*MGlpT2 + MStop2(1)

	Cd(2026) = 19*MTy2 + MStop2(1)

	Cd(2031) = 6*MStop2(2)**2 + MStop2(3)*MStop2(4)

	Cd(2036) = 2*MTy2 + MStop2(2)

	Cd(202) = 2*(MTy2*MStop2(1) + Cd(155)*MStop2(2))

	Cd(203) = Cd(202) - 2*MGl2*MStop2(1)

	Cd(2058) = MGlpT2 - 3*MStop2(1)

	Cd(2064) = 13*MTy2 - MStop2(1)

	Cd(2065) = 5*MTy2 + 7*MStop2(1)

	Cd(2066) = Cd(2064)*MStop2(1) - Cd(2065)*MStop2(2)

	Cd(2067) = 11*MTy2 + MStop2(1)

	Cd(2068) = Cd(2067)*MStop2(1) - 3*Cd(1592)*MStop2(2)

	Cd(2081) = 7*MTy2 - 2*MStop2(1)

	Cd(19) = A0delStop(1) + Cr(13)*Cr(34)*MStop2(1)

	Cd(2084) = 2*MTy2 - 3*MStop2(1)

	Cd(2086) = 3*MTy2 - 2*MStop2(1)

	Cd(2092) = 2*MTy2 - MStop2(1)

	Cd(207) = Cd(186) - 2*MGl2*MStop2(1)

	Cd(2093) = 13*MGlpT2 - 6*MStop2(1)

	Cd(2098) = 8*MStop2(1) + MStop2(3)

	Cd(2100) = 11*MTy2 - 9*MStop2(1)

	Cd(2106) = 7*MGlpT2 - 5*MStop2(1)

	Cd(2114) = 2*MStop2(3) + MStop2(4)

	Cd(2116) = 10*MStop2(2) - 11*MStop2(3)

	Cd(2122) = 4*(MTy2 - MStop2(1))

	Cd(2129) = -MGl2 + Cd(1592)

	Cd(2133) = 8*MStop2(2) - 11*MStop2(3)

	Cd(2135) = 5*MTy2 - 17*MStop2(1)

	Cd(2136) = Cd(2135)*MStop2(1) + 3*Cd(1592)*MStop2(2)

	Cd(2140) = -2*MGl2 + Cd(2133)

	Cd(2151) = 5*MTy2 - 6*MStop2(1)

	Cd(2155) = 13*MTy2 - 14*MStop2(1)

	Cd(2157) = 4*(MTy2 - MStop2(2))

	Cd(2164) = -MGl2 + Cd(1559)

	Cd(2169) = 8*MStop2(2) + 3*MStop2(3)

	Cd(2171) = 5*MTy2 + 9*MStop2(1)

	Cd(2175) = -2*MGl2 + Cd(2169)

	Cd(20) = A0delStop(2)*MStop2(1) + Cd(19)*MStop2(2)

	Cd(2187) = MTy2 + Cd(118)

	Cd(2193) = MTy2 - 3*MStop2(1)

	Cd(2196) = 10*MTy2 - MStop2(4)

	Cd(2200) = 3*MTy2 + MStop2(1)

	Cd(2201) = 4*MStop2(2) - 7*MStop2(3)

	Cd(2202) = 2*MTy2 + Cd(2201)

	Cd(2203) = 6*MStop2(2) - 7*MStop2(3)

	Cd(2204) = MTy2*Cd(2202) - Cd(2203)*MStop2(1)

	Cd(2205) = 4*Cd(1592) - MStop2(3)

	Cd(2210) = -4*MGl2 + Cd(2205)

	Cd(2211) = 2*Cd(2204) + MGl2*Cd(2210)

	Cd(2220) = 8*MStop2(2) - 9*MStop2(3)

	Cd(2221) = MTy2 - 5*MStop2(1)

	Cd(2222) = 3*Cd(2221)*MStop2(1) + Cd(2065)*MStop2(2)

	Cd(2231) = -2*MGl2 + Cd(2220)

	Cd(2241) = MTy2 + Cd(138)

	Cd(2248) = 11*MTy2 - 18*MStop2(1)

	Cd(2251) = 24*MTy2 - 5*MStop2(1)

	Cd(2252) = 6*MTy2 - MStop2(1)

	Cd(2254) = 3*MTy2 + MStop2(2)

	Cd(2255) = 4*Cd(1559) + MStop2(3)

	Cd(2257) = 4*MStop2(2) + 3*MStop2(3)

	Cd(2258) = 2*MTy2 + Cd(2257)

	Cd(2259) = MTy2*Cd(2258) - Cd(1569)*MStop2(2)

	Cd(2263) = -4*MGl2 + Cd(2255)

	Cd(2264) = 2*Cd(2259) + MGl2*Cd(2263)

	Cd(2273) = 8*MStop2(2) + MStop2(3)

	Cd(2274) = 4*MTy2 - MStop2(1)

	Cd(2275) = 9*MTy2 + 5*MStop2(1)

	Cd(2277) = 3*MTy2 + 7*MStop2(1)

	Cd(2286) = -2*MGl2 + Cd(2273)

	Cd(2309) = MGl2 - 7*MTy2

	Cd(2321) = 5*MGl2 + MTy2

	Cd(2326) = 3*MGl2 + 4*MTy2

	Cd(2333) = 9*MStop2(2) - 11*MStop2(3)

	Cd(2344) = MStop2(3) + 6*MStop2(4)

	Cd(2345) = 10*MStop2(2)**2 - Cd(2344)*MStop2(3)

	Cd(2346) = 5*MStop2(3) + 12*MStop2(4)

	Cd(2347) = 14*MStop2(2)**2 - Cd(2346)*MStop2(3)

	Cd(2349) = 3*MStop2(3) + 14*MStop2(4)

	Cd(2350) = 20*MStop2(2)**2 - Cd(2349)*MStop2(3)

	Cd(2351) = 13 - 9*UStop2(1,2)

	Cd(2355) = M_3 - 2*M_3C*UCStop(1,3)**2

	Cd(2356) = 4*M_3*MStop2(1) - M_3C*Cd(1421)*UCStop(1,3)**2

	Cd(2358) = 4*M_3*MStop2(1) - M_3C*Cd(1583)*UCStop(1,3)**2

	Cd(2359) = 3 - 2*UStop2(1,2)**2

	Cd(2360) = 11 - 12*UStop2(1,2)

	Cd(2361) = 5 - 4*UStop2(1,2)

	Cd(2364) = 1 - UStop2(1,2)**2

	Cd(2365) = 1 + UStop2(1,2)**2

	Cd(2366) = 3 + UStop2(1,2)**2

	Cd(2368) = 7 - 4*UStop2(1,2)

	Cd(2370) = 5 - 3*UStop2(1,2)

	Cd(2371) = -(Cd(2370)*MStop2(1)) + Cr(85)*MStop2(2)

	Cd(2373) = 17 - 8*UStop2(1,2)

	Cd(2377) = 5 - 3*UStop2(1,2)**2

	Cd(2378) = 17 - 24*UStop2(1,2)

	Cd(2380) = 5 - UStop2(1,2)**2

	Cd(2384) = 6 + UStop2(1,2)**2

	Cd(2385) = 5 + 2*UStop2(1,2)

	Cd(2386) = 13 - 6*UStop2(1,2)

	Cd(2387) = -(Cd(2386)*MStop2(1)) + Cd(2385)*MStop2(2)

        Cd(2388) = M_3*Cd(2384)*MStop2(3) - 
     &   2*M_3C*Cd(2387)*UCStop(1,3)**2

	Cd(2391) = 5 - 8*UStop2(1,2)

	Cd(2392) = 35 - 8*UStop2(1,2)

	Cd(2393) = 1 - 16*UStop2(1,2)

	Cd(2413) = 7 - 5*UStop2(1,2)

	Cd(2415) = 7 - 11*UStop2(1,2)

	Cd(2417) = 7 - 3*UStop2(1,2)

	Cd(2418) = 5 - 9*UStop2(1,2)

	Cd(2419) = Cd(2417)*MStop2(1) + Cd(2418)*MStop2(2)

	Cd(2421) = 7 - 9*UStop2(1,2)

	Cd(2422) = Cr(85)*MStop2(1) + Cd(2421)*MStop2(2)

	Cd(2423) = 13 - 11*UStop2(1,2)

	Cd(2424) = -(Cd(2423)*MStop2(1)) + Cd(2370)*MStop2(2)

	Cd(2425) = MTy2*Cd(2424) + Cd(2422)*MStop2(1)

	Cd(2429) = 1 - 3*UStop2(1,2)**2

	Cd(2432) = 5 - 2*UStop2(1,2)

	Cd(2436) = 21 - 22*UStop2(1,2)

	Cd(2437) = Cd(2436)*MStop2(1) + Cd(1052)*MStop2(2)

	Cd(2440) = M_3*Cd(2366) - 4*M_3C*Cr(85)*UCStop(1,3)**2

	Cd(2441) = 5 - 20*UStop2(1,2)

	Cd(2442) = Cd(2441)*MStop2(1) - Cd(2361)*MStop2(2)

        Cd(2443) = 4*M_3*Cd(2365)*MStop2(1) + 
     &   M_3C*Cd(2442)*UCStop(1,3)**2

	Cd(2444) = 17 - 20*UStop2(1,2)

	Cd(2445) = Cd(2444)*MStop2(1) + Cd(2368)*MStop2(2)

        Cd(2446) = 6*M_3*Cd(2364)*MStop2(1) - 
     &   M_3C*Cd(2445)*UCStop(1,3)**2

	Cd(2468) = 9*MGlpT2 + 5*MStop2(1)

	Cd(244) = -(Cd(169)*MStop2(1)) + 2*MStop2(2)

	Cd(2472) = 9*MStop2(2) + 2*MStop2(3)

	Cd(245) = 2*(MGl2*Cd(173) + Cd(244)*MStop2(1))

	Cd(2487) = Cd(2418)*MStop2(1) + Cd(2417)*MStop2(2)

	Cd(2489) = Cd(2421)*MStop2(1) + Cr(85)*MStop2(2)

	Cd(2490) = Cd(2370)*MStop2(1) - Cd(2423)*MStop2(2)

	Cd(2491) = MTy2*Cd(2490) + Cd(2489)*MStop2(2)

	Cd(2498) = 17 - 2*UStop2(1,2)

	Cd(2501) = 7 - 2*UStop2(1,2)

	Cd(2502) = Cd(1052)*MStop2(1) - Cd(2501)*MStop2(2)

        Cd(2503) = M_3*Cd(2429)*MStop2(2) - 
     &   M_3C*Cd(2502)*UCStop(1,3)**2

	Cd(2505) = M_3*Cd(2366) - 8*M_3C*UCStop(1,3)**2

        Cd(2506) = 6*M_3*Cd(2364)*MStop2(2) + 
     &   M_3C*Cd(2368)*MStop2(3)*UCStop(1,3)**2

	Cd(2507) = 11 + 4*UStop2(1,2)

	Cd(2508) = Cd(2361)*MStop2(1) + Cd(2507)*MStop2(2)

        Cd(2509) = 4*M_3*Cd(2365)*MStop2(2) - 
     &   M_3C*Cd(2508)*UCStop(1,3)**2

	Cd(2518) = 6*MTy2 + MStop2(1)

	Cd(2534) = 11*MGlpT2 - 18*MStop2(1)

	Cd(2538) = 6*MGl2 + 5*MTy2

	Cd(2546) = 5*MStop2(3) + 2*MStop2(4)

	Cd(2547) = 14*MStop2(2)**2 - Cd(2546)*MStop2(3)

	Cd(2548) = 10*MStop2(2)**2 - Cd(1121)*MStop2(3)

	Cd(2550) = 20*MStop2(2)**2 - 3*Cd(1755)*MStop2(3)

	Cd(251) = MGlmT2 + MStop2(2)

	Cd(2554) = 4*M_3*MStop2(2) - M_3C*Cd(1122)*UCStop(1,3)**2

	Cd(2557) = 17 - 16*UStop2(1,2)

	Cd(2560) = 7*MStop2(2) - MStop2(3)

	Cd(2561) = 55 - 48*UStop2(1,2)

	Cd(2562) = Cd(2560)*MStop2(1) - Cd(2561)*MStop2(2)**2

	Cd(2563) = MStop2(1) - Cd(1052)*MStop2(2)

        Cd(2564) = M_3*Cd(2380)*MStop2(3) + 
     &   8*M_3C*Cd(2563)*UCStop(1,3)**2

	Cd(2567) = 5*MStop2(2) - 3*MStop2(3)

	Cd(2568) = 13 - 8*UStop2(1,2)

	Cd(2569) = Cd(2567)*MStop2(1) - Cd(2568)*MStop2(2)**2

	Cd(2572) = 13 - 16*UStop2(1,2)

	Cd(2575) = 9*MStop2(2) + MStop2(3)

	Cd(2576) = 7 - 16*UStop2(1,2)

	Cd(2577) = Cd(2575)*MStop2(1) + Cd(2576)*MStop2(2)**2

        Cd(2578) = 6*M_3*Cd(2365)*MStop2(2)*MStop2(3) + 
     &   3*M_3C*Cd(2577)*UCStop(1,3)**2

	Cd(2580) = 1 - 8*UStop2(1,2)

	Cd(2581) = 7*MStop2(1) + Cd(2580)*MStop2(2)

        Cd(2582) = M_3*Cd(2384)*MStop2(3) + 
     &   2*M_3C*Cd(2581)*UCStop(1,3)**2

	Cd(2583) = 13 - 40*UStop2(1,2)

	Cd(2584) = 7 + 8*UStop2(1,2)

	Cd(2587) = 41 - 32*UStop2(1,2)

	Cd(2588) = 3*Cd(1898)*MStop2(1) - Cd(2587)*MStop2(2)**2

        Cd(2589) = 6*M_3*Cd(1136)*MStop2(2)*MStop2(3) + 
     &   M_3C*Cd(2588)*UCStop(1,3)**2

	Cd(2596) = 4*MGlpT2 - 3*MStop2(1)

	Cd(2597) = MGlmT2**2 - Cd(2596)*MStop2(1)

	Cd(2600) = 5*MGl2 + 3*MTy2

	Cd(2601) = 9*MGl2 + 5*MTy2

	Cd(256) = MStop2(2) - 2*MStop2(3)

	Cd(2613) = 7*MGl2 + 4*MTy2

	Cd(2617) = MGl2 - 3*MTy2

	Cd(2627) = 5*MGl2 - 2*MTy2

        Cd(2636) = 4*MGlpTmSt2(1)*MStop2(1)**2 - 
     &   Cd(152)*Cd(2597)*UStop2(1,2)**2

	Cd(261) = MGlpTmSt2(2)*MStop2(2) + MTy2*Cd(251)*UStop2(1,2)

	Cd(2651) = 29*MGl2 + 20*MTy2

	Cd(2655) = Cd(164) + 2*MGlpTmSt2(1)*UStop2(1,2)

        Cd(2656) = 4*Cd(2655)*MStop2(1)**2 - 
     &   Cd(152)*Cd(2597)*UStop2(1,2)**2

        Cd(2659) = 4*MGlpTmSt2(1)*MStop2(1)**2 - 
     &   Cd(152)*Cd(2597)*UStop2(1,2)

	Cd(2680) = 4*MGlpT2 - 3*MStop2(2)

	Cd(25) = A0delStop(1) - A0delStop(2)

	Cd(2681) = MGlmT2**2 - Cd(2680)*MStop2(2)

	Cd(2686) = MTy2 - 3*MStop2(2)

        Cd(2712) = 4*MGlpTmSt2(2)*MStop2(2)**2 - 
     &   Cd(179)*Cd(2681)*UStop2(1,2)**2

	Cd(2731) = Cd(251) + 2*MGlpTmSt2(2)*UStop2(1,2)

        Cd(2732) = -4*Cd(2731)*MStop2(2)**2 + 
     &   Cd(179)*Cd(2681)*UStop2(1,2)**2

        Cd(2735) = 4*MGlpTmSt2(2)*MStop2(2)**2 - 
     &   Cd(179)*Cd(2681)*UStop2(1,2)

	Cd(2763) = 3*MTy2 + 4*MStop2(1)

	Cd(2774) = 3*MTy2 + MStop2(3)

	Cd(2775) = MTy2*Cd(2774) - 3*MStop2(1)**2

	Cd(2796) = 2*Cd(1064)*MStop2(1)**2 - 2*MTy2**2*MStop2(4)

	Cd(2802) = 5*MStop2(2) - 2*MStop2(3)

	Cd(2803) = 10*MStop2(2) - 7*MStop2(3)

	Cd(2808) = 3*MGl2 + MTy2

	Cd(2829) = 13*MGl2 - 8*MTy2

	Cd(2834) = 9*MTy2 + 10*MStop2(1)

	Cd(2835) = Cd(2834)*MStop2(1) + 3*Cd(155)*MStop2(2)

	Cd(280) = 6*MTy2 - MStop2(4)

	Cd(281) = MGl2 + Cd(280)

	Cd(283) = MTy2 - MStop2(4)

	Cd(2874) = 2*MStop2(1) + Cd(152)*UStop2(1,2)**2

	Cd(2875) = 6*MTy2 - MStop2(3)

	Cd(2876) = 26*MStop2(2) - 27*MStop2(3)

	Cd(2877) = 18*MTy2 + Cd(2876)

	Cd(2878) = MTy2*Cd(2877) - 3*MStop2(1)*MStop2(3)

	Cd(2887) = MTy2**2 + 7*MStop2(1)**2

	Cd(2893) = 2*Cd(155)*MStop2(1)*MStop2(2) + MTy2**2*MStop2(4)

	Cd(2894) = 4*MStop2(2) + MStop2(3)

	Cd(2900) = Cd(1684) + MStop2(1)

	Cd(286) = 3*MStop2(1) + MStop2(3)

	Cd(2909) = 9*MTy2 + 14*MStop2(1)

	Cd(2910) = Cd(2909)*MStop2(1) - 3*Cd(155)*MStop2(2)

	Cd(2922) = 2*MTy2 + 5*MStop2(1)

	Cd(2926) = 8*MTy2 - MStop2(1)

	Cd(2930) = 4*MTy2 - MStop2(3)

	Cd(2931) = 20*MStop2(2) - 21*MStop2(3)

	Cd(2932) = 12*MTy2 + Cd(2931)

	Cd(2933) = MTy2*Cd(2932) - 3*MStop2(1)*MStop2(3)

	Cd(2936) = -(MGl2*Cd(2930)) + Cd(2933)

	Cd(2943) = 4*MStop2(2) - 6*MStop2(3)

	Cd(290) = 1 - 2*UStop2(1,2)

	Cd(2953) = MTy2**2 + 3*MStop2(1)**2

	Cd(2955) = 12*MTy2 + Cd(1148)

	Cd(2956) = MTy2*Cd(2955) + 9*MStop2(1)*MStop2(3)

	Cd(2957) = 4*MTy2 + 3*MStop2(3)

	Cd(2958) = MTy2 + 6*MStop2(1)

	Cd(2964) = 4*MTy2 - 3*MStop2(3)

	Cd(2965) = 12*MStop2(2) - 17*MStop2(3)

	Cd(2966) = 12*MTy2 + Cd(2965)

	Cd(2967) = MTy2*Cd(2966) - 9*MStop2(1)*MStop2(3)

	Cd(2970) = 3*MGl2 - MStop2(1)

	Cd(2971) = 16*MTy2*MGlpTmSt2(1) - 3*Cd(2970)*MStop2(1)

	Cd(2973) = 4*MTy2 + 3*MStop2(1)

	Cd(3) = 2*Cd(1) + Cd(2)*Cr(34)

	Cd(293) = 2*Cr(54)*MStop2(1) - Cr(47)*MStop2(2)

	Cd(2985) = 6*MStop2(2) + 7*MStop2(3)

	Cd(2987) = MTy2 + 5*MStop2(1)

	Cd(2991) = 24*MStop2(1) + MStop2(3)

	Cd(2992) = 18*MTy2 + Cd(2991)

	Cd(2993) = MTy2*Cd(2992) - 4*MStop2(1)*MStop2(3)

	Cd(2994) = 5*MTy2 - 3*MStop2(1)

	Cd(2995) = 5*MTy2 + 17*MStop2(1)

	Cd(3000) = -3*MGl2 + 2*Cd(2973)

	Cd(3006) = MTy2**2 - 3*MStop2(1)**2

	Cd(3009) = 11*MGl2 + 3*MTy2

	Cd(3014) = 72*MStop2(2) - 31*MStop2(3)

	Cd(3015) = 9*MStop2(2) - 7*MStop2(3)

	Cd(3016) = 8*MStop2(2) - MStop2(3)

	Cd(3018) = 48*MStop2(2) - 23*MStop2(3)

	Cd(3019) = 8*MStop2(2) + 3*MStop2(3)

	Cd(3021) = 16*MStop2(2) - 17*MStop2(3)

	Cd(3022) = 72*MStop2(2) - 41*MStop2(3)

	Cd(3024) = 4*MStop2(2) - 17*MStop2(3)

	Cd(3026) = 6*MStop2(2) - 11*MStop2(3)

	Cd(3028) = 4*MStop2(2) - 19*MStop2(3)

	Cd(3033) = 7*MGlpT2*MTy2 + 3*MGl2*MGlpTmSt2(1)

	Cd(3035) = 4*MGl2 + MTy2

	Cd(3062) = MTy2*Cd(1647) - 2*MStop2(1)**2

	Cd(302) = 2*MGlpT2 - MStop2(4)

	Cd(304) = MStop2(4) + 2*Cd(302)*UStop2(1,2)

	Cd(3099) = -MGl2 + Cd(2973)

	Cd(306) = MTy2 + 2*Cr(1)

	Cd(307) = 3*MStop2(1) + MStop2(3)

	Cd(3122) = 3*MTy2 + 4*MStop2(2)

	Cd(3127) = 8*MTy2 + 3*MStop2(1)

	Cd(308) = MTy2*Cd(306) - Cd(307)*MStop2(2)

	Cd(3129) = 11*MTy2 - 2*MStop2(1)

	Cd(3130) = 14*MTy2 - 3*MStop2(1)

	Cd(3132) = 3*MGl2 - 2*MTy2

	Cd(3137) = 3*MTy2 - MStop2(3)

	Cd(3138) = -(MTy2*Cd(3137)) + 3*MStop2(2)**2

	Cd(3139) = 9*MTy2 + MStop2(1)

	Cd(3149) = 2*Cd(1657)*MStop2(2) + 3*MTy2*MStop2(4)

	Cd(3150) = 2*Cd(181)*MStop2(1)*MStop2(2) + MTy2**2*MStop2(4)

	Cd(3152) = Cd(3149) - MGl2*MStop2(1)

	Cd(3153) = 3*Cd(3150) - MGl2*Cd(3152)

	Cd(3159) = -(MTy2*Cd(1952)) + MStop2(2)**2

	Cd(3161) = -MGl2 + Cd(3122)

	Cd(3189) = 9*MGl2 + 2*MTy2

	Cd(315) = 2*MStop2(2) + MStop2(3)

	Cd(3208) = 2*MTy2**2*Cd(1226) - 6*MStop2(2)**2*MStop2(4)

	Cd(3221) = 24*MTy2 + 5*MStop2(1)

	Cd(3222) = 27*MTy2 - 2*MStop2(1)

	Cd(3226) = 16*MStop2(2) + MStop2(3)

	Cd(3229) = 2*MStop2(2) + Cd(179)*UStop2(1,2)**2

	Cd(3230) = 26*MStop2(2) + MStop2(3)

	Cd(3231) = 18*MTy2 + Cd(3230)

	Cd(3232) = MTy2*Cd(3231) + 3*MStop2(2)*MStop2(3)

	Cd(3233) = 6*MTy2 + MStop2(3)

	Cd(3234) = 3*MGl2 - MTy2

	Cd(3235) = 6*MTy2*Cd(3234) + MGlmT2*MStop2(1)

	Cd(3236) = 27*MGl2 - MTy2

	Cd(3237) = MTy2*Cd(3236) - Cd(153)*MStop2(1)

	Cd(3240) = 3*MGl2 - MStop2(2)

	Cd(3250) = 6*MTy2 - 7*MStop2(1)

	Cd(3252) = 4*MTy2 + 5*MStop2(1)

	Cd(3256) = 3*Cd(3150) + MGl2*Cd(2617)*MStop2(1)

	Cd(3258) = MGl2 - 3*Cd(155)

	Cd(3262) = Cd(1845) + MStop2(2)

	Cd(31) = 8*A0delT*CA - 4*CA*MTy2*Cr(15)

	Cd(3286) = 12*MStop2(2) + 5*MStop2(3)

	Cd(3288) = 6*MTy2 + 7*MStop2(1)

	Cd(3290) = 23*MTy2 - 2*MStop2(1)

	Cd(3292) = 8*MTy2 - 3*MStop2(2)

	Cd(3294) = 6*MStop2(2) - 3*MStop2(3)

	Cd(3296) = 20*MStop2(2) + MStop2(3)

	Cd(3297) = 12*MTy2 + Cd(3296)

	Cd(3298) = MTy2*Cd(3297) + 3*MStop2(2)*MStop2(3)

	Cd(3308) = 2*MTy2**2*Cd(315) - 2*Cd(1160)*MStop2(2)**2

	Cd(3316) = MTy2**2 + 3*MStop2(2)**2

	Cd(3318) = 12*MStop2(2) - 5*MStop2(3)

	Cd(3319) = 12*MTy2 + Cd(3318)

	Cd(3320) = MTy2*Cd(3319) - 9*MStop2(2)*MStop2(3)

	Cd(3326) = 10*MGl2 + 13*MTy2

	Cd(3342) = MGl2 - 7*MStop2(1)

	Cd(3346) = 24*MStop2(2) - MStop2(3)

	Cd(3347) = 18*MTy2 + Cd(3346)

	Cd(3348) = MTy2*Cd(3347) + 4*MStop2(2)*MStop2(3)

	Cd(3349) = 8*MTy2 - MStop2(2)

	Cd(3351) = 12*MTy2 - 5*MStop2(1)

	Cd(3355) = 4*MTy2 + 3*MStop2(2)

	Cd(3356) = 5*MTy2 - 3*MStop2(2)

	Cd(3357) = 18*MTy2 + 5*MStop2(1)

	Cd(3358) = 13*MTy2 + 17*MStop2(1)

	Cd(3362) = -3*MGl2 + 2*Cd(3355)

	Cd(3368) = MTy2**2 - 3*MStop2(2)**2

	Cd(3372) = 3*MGl2*MGlmT2 + MTy2*MStop2(1)

	Cd(3373) = 3*MGl2 + 2*MTy2

	Cd(3377) = 36*MStop2(2) - 41*MStop2(3)

	Cd(3378) = 9*MStop2(2) - 4*MStop2(3)

	Cd(3380) = 36*MStop2(2) - 31*MStop2(3)

	Cd(3382) = 4*MStop2(2) - 11*MStop2(3)

	Cd(3383) = 24*MStop2(2) - 25*MStop2(3)

	Cd(3386) = 4*MStop2(2) + 13*MStop2(3)

	Cd(3387) = 4*MStop2(2) + 15*MStop2(3)

	Cd(3389) = 6*MStop2(2) + 5*MStop2(3)

	Cd(3394) = 7*MGlpT2*MTy2 + 3*MGl2*MGlpTmSt2(2)

	Cd(3399) = 12*MTy2 + Cd(3286)

	Cd(3400) = MTy2*Cd(3399) + 9*MStop2(2)*MStop2(3)

	Cd(3401) = 4*MTy2*Cd(3234) - MGlpT2*MStop2(1)

	Cd(3402) = MTy2*Cd(2601) - 3*Cd(153)*MStop2(1)

	Cd(3409) = 13*MTy2 - MStop2(2)

	Cd(3426) = 4*MTy2*Cd(3234) - MGlmT2*MStop2(1)

	Cd(3427) = 7*MGl2 - 11*MTy2

	Cd(3428) = MTy2*Cd(3427) - Cd(153)*MStop2(1)

	Cd(3432) = Cd(2735) - 2*MTy2*Cd(251)*MStop2(2)*UStop2(1,2)

	Cd(3434) = 11*MGl2 - 9*MTy2

	Cd(3435) = 2*MGl2 - 3*MTy2

	Cd(3436) = 3*MGl2*Cd(3434) + Cd(3435)*MStop2(1)

	Cd(3439) = 5*MTy2**3 + Cd(138)**3

	Cd(3451) = 8*MTy2*MGlpTmSt2(2) - Cd(3240)*MStop2(2)

	Cd(3470) = -(MTy2*Cd(1667)) + 2*MStop2(2)**2

	Cd(3476) = -MGl2 + Cd(3355)

	Cd(3491) = 7*MTy2 + 9*MStop2(1)

	Cd(3504) = MStop2(2) + 9*MStop2(3)

	Cd(3505) = 2*MStop2(1)**2 + Cd(1795)*MStop2(3)

	Cd(3506) = 2*MStop2(2) + 13*MStop2(3)

	Cd(3507) = 2*MStop2(1)**2 + Cd(1078)*MStop2(3)

	Cd(3509) = 15*MStop2(1) + MStop2(3)

	Cd(3510) = 3*MStop2(1)**2 + Cd(3509)*MStop2(3)

	Cd(3511) = 8*MStop2(3) - 3*MStop2(4)

	Cd(3512) = 2*MStop2(2)**2 - Cd(3511)*MStop2(3)

	Cd(3514) = MStop2(3) - MStop2(4)

	Cd(3515) = MStop2(2)**2 - 4*Cd(3514)*MStop2(3)

	Cd(3516) = 7*MStop2(3) - 8*MStop2(4)

	Cd(3517) = 2*MStop2(2)**2 - Cd(3516)*MStop2(3)

	Cd(3546) = 2*(MTy2 - MStop2(2))

	Cd(3562) = MGl2 - 4*MTy2

	Cd(3569) = MStop2(3) - 5*MStop2(4)

	Cd(3570) = 8*MStop2(2)**2 + Cd(3569)*MStop2(3)

	Cd(3591) = 8*MStop2(2)**2 - Cd(2273)*MStop2(3)

	Cd(3662) = MStop2(3) + MStop2(4)

	Cd(3665) = MTy2 - 6*MStop2(1)

	Cd(3669) = -(MTy2*Cd(1952)) + MStop2(1)**2

	Cd(35) = 1 - 2*Log(MTy2/Q)

	Cd(3680) = 2*UUStop(1,1) - UUStop(2,1)

	Cd(3681) = -(Cd(3680)*MStop2(1)) + 2*MStop2(2)*UUStop(1,1)

	Cd(3682) = 48*UUStop(1,1) + UUStop(2,1)

	Cd(3683) = Cd(3682)*MStop2(2) - MStop2(1)*UUStop(2,1)

	Cd(3684) = Cd(3683)*MStop2(1) + 48*MGl2*MStop2(2)*UUStop(1,1)

	Cd(3685) = 25*UUStop(1,1) + UUStop(2,1)

	Cd(3686) = 25*UUStop(1,1) + 2*UUStop(2,1)

	Cd(3688) = 45*UUStop(1,1) + UUStop(2,1)

	Cd(3689) = 45*UUStop(1,1) + 2*UUStop(2,1)

	Cd(3692) = 39*UUStop(1,1) - UUStop(2,1)

	Cd(3693) = 66*UUStop(1,1) - UUStop(2,1)

	Cd(3695) = 19*UUStop(1,1) - 3*UUStop(2,1)

	Cd(3696) = 14*UUStop(1,1) - 3*UUStop(2,1)

	Cd(3698) = 25*UUStop(1,1) + 3*UUStop(2,1)

	Cd(3699) = 26*UUStop(1,1) + 3*UUStop(2,1)

	Cd(3702) = 17*UUStop(1,1) - UUStop(2,1)

	Cd(3703) = 29*UUStop(1,1) - 2*UUStop(2,1)

	Cd(3705) = UUStop(1,1) - UUStop(2,1)

	Cd(3706) = 5*UUStop(1,1) + 2*UUStop(2,1)

	Cd(3718) = 5*MTy2 - MStop2(1)

	Cd(3724) = MGlpT2 - 4*MStop2(1)

	Cd(368) = CA*MUE - AtC*SA

	Cd(3741) = -3*MGl2 + 2*Cd(160)

	Cd(369) = CA*MUEC - At*SA

	Cd(3757) = 3*MTy2 + 2*MStop2(1)

	Cd(3761) = 5*MTy2 + MStop2(1)

	Cd(3770) = MTy2 - 4*MStop2(1)

	Cd(3781) = MStop2(2) + 8*MStop2(3)

	Cd(3782) = 2*MStop2(2)**2 + Cd(1749)*MStop2(3)

	Cd(3784) = MStop2(2) + 6*MStop2(3)

	Cd(3785) = MStop2(2)**2 + Cd(1226)*MStop2(3)

	Cd(3787) = 2*MStop2(2)**2 + Cd(1759)*MStop2(3)

	Cd(3788) = 2*MStop2(2)**2 + Cd(1764)*MStop2(3)

	Cd(3789) = 5*MStop2(3) + 6*MStop2(4)

	Cd(3790) = 3*MStop2(2)**2 + Cd(3789)*MStop2(3)

	Cd(3814) = 2*(MTy2 - MStop2(1))

	Cd(3830) = 8*MStop2(2)**2 + Cd(1032)*MStop2(3)

	Cd(3842) = MTy2 - 9*MStop2(1)

	Cd(3846) = 2*MTy2 - 7*MStop2(1)

	Cd(3852) = 4*MStop2(2)**2 - Cd(1141)*MStop2(3)

	Cd(3856) = 3*MTy2 - 7*MStop2(1)

	Cd(3858) = 7*MTy2 + 2*MStop2(1)

	Cd(3885) = MGl2 - 2*MTy2

	Cd(3888) = MGl2 + 18*MStop2(1)

	Cd(3937) = -(Cd(3680)*MStop2(2)) + 2*MStop2(1)*UUStop(1,1)

	Cd(3938) = 24*UUStop(1,1) + UUStop(2,1)

	Cd(3939) = 72*UUStop(1,1) + UUStop(2,1)

	Cd(3940) = Cd(3938)*MStop2(1) - Cd(3939)*MStop2(2)

	Cd(3941) = -(Cd(3940)*MStop2(2)) + 24*MGl2*Cd(315)*UUStop(1,1)

	Cd(3942) = 7*UUStop(1,1) + UUStop(2,1)

	Cd(3943) = 61*UUStop(1,1) + 2*UUStop(2,1)

	Cd(3945) = 33*UUStop(1,1) + UUStop(2,1)

	Cd(3946) = 69*UUStop(1,1) + 2*UUStop(2,1)

	Cd(3949) = 45*UUStop(1,1) - UUStop(2,1)

	Cd(3950) = 48*UUStop(1,1) - UUStop(2,1)

	Cd(3952) = 31*UUStop(1,1) - 3*UUStop(2,1)

	Cd(3953) = 22*UUStop(1,1) + 3*UUStop(2,1)

	Cd(3955) = 10*UUStop(1,1) - 3*UUStop(2,1)

	Cd(3956) = 37*UUStop(1,1) + 3*UUStop(2,1)

	Cd(3959) = 5*UUStop(1,1) - UUStop(2,1)

	Cd(3961) = 7*UUStop(1,1) - UUStop(2,1)

	Cd(3969) = 5*MTy2 - MStop2(2)

	Cd(4) = A0delStop(2)*dMSfsq1233eps(-1) - Cd(3)*MStop2(2)

	Cd(3975) = MGlpT2 - 4*MStop2(2)

	Cd(3991) = 2*MTy2 + 3*MStop2(2)

	Cd(3992) = 13*MTy2 + 9*MStop2(1)

	Cd(3996) = -3*MGl2 + 2*Cd(3991)

	Cd(4032) = 3*MGl2 + 5*MTy2

	Cd(4033) = 4*MGl2 + 7*MTy2

	Cd(4081) = 3*MTy2 - MStop2(1)

	Cd(4089) = 2*MGl2 + 9*MTy2

	Cd(402) = Cd(368)*UCStop(3,3) + Cd(369)*UCStop(3,4)

	Cd(403) = Cd(369)*UCStopC(3,3) + Cd(368)*UCStopC(3,4)

	Cd(4117) = MTy2 - 2*MStop2(2)

	Cd(404) = Cd(402)*Conjugate(dY33eps(-1)) + Cd(403)*dY33eps(-1)

	Cd(405) = Cd(402)*Conjugate(dY33eps(0)) + Cd(403)*dY33eps(0)

	Cd(406) = Cd(402)*Conjugate(dY33eps(1)) + Cd(403)*dY33eps(1)

	Cd(407) = 2*Cd(406) + Cd(405)*Cr(13)

	Cd(408) = 1/2.D0*(A0delStop(1)*Cd(404)) + Cd(407)*MStop2(1)

	Cd(409) = 2*Cd(406) + Cd(405)*Cr(34)

	Cd(4178) = 3*MTy2 - MStop2(2)

	Cd(410) = 1/2.D0*(A0delStop(2)*Cd(404)) + Cd(409)*MStop2(2)

	Cd(4193) = 5*MTy2 + 3*MStop2(2)

	Cd(4216) = -MGl2 + Cd(2200)

	Cd(4225) = 12*MStop2(1) + MStop2(3)

	Cd(4240) = MTy2*Cd(286) + 3*MStop2(1)**2

	Cd(4241) = 2*MTy2*MStop2(1) + Cd(2193)*MStop2(2)

	Cd(4243) = MTy2**2 + MStop2(1)**2

	Cd(4244) = Cd(286)*Cd(4243) + 2*MTy2*MStop2(1)**2

	Cd(4245) = -(MGl2*Cd(4240)) + 3*Cd(4244)

	Cd(4260) = 10*MStop2(2) - 7*MStop2(3)

	Cd(4263) = MGl2 + 3*MTy2

	Cd(4267) = 9*MStop2(2) - 7*MStop2(3)

	Cd(4269) = 9*MStop2(1) + MStop2(3)

	Cd(41) = 4 - 3*Log(MTy2/Q)

	Cd(4279) = Cd(4081)*MStop2(1) + Cd(2193)*MStop2(2)

	Cd(4310) = MGl2*MTy2 - 3*Cd(4243)

        Cd(423) = 2*(Cr(28)*MStop2(2) + 
     &     MStop2(3)*UStop2(1,1)*UStop2(1,2))

	Cd(424) = At*MUE + AtC*MUEC

	Cd(4324) = 11*MTy2 + 5*MStop2(1)

	Cd(4330) = MGlpTmSt2(1)*MStop2(1) + MTy2*Cd(164)*UStop2(1,2)**2

	Cd(425) = AtC**2*UCStop(1,3)**2 + At**2*UCStopC(1,3)**2

	Cd(426) = AtC*MUE*UCStop(1,3)**2 + At*MUEC*UCStopC(1,3)**2

	Cd(4347) = 4*MTy2*MStop2(1)**2 + 3*Cd(4243)*MStop2(4)

	Cd(4348) = Cd(4347) - 2*MGl2*MStop2(1)**2

	Cd(427) = MUE**2*UCStop(1,3)**2 + MUEC**2*UCStopC(1,3)**2

	Cd(4352) = 2*(MTy2 - 3*MStop2(1))

	Cd(4362) = 3*MTy2 - 5*MStop2(1)

	Cd(4363) = Cd(4362)*MStop2(1) - Cd(2193)*MStop2(2)

	Cd(42) = 3 - Cd(41)*Log(MTy2/Q)

	Cd(4377) = MTy2*Cd(2086) + Cd(2973)*MStop2(2)

	Cd(4379) = 12*MStop2(2) - MStop2(3)

	Cd(4394) = -(MTy2*Cd(2930)) + Cd(1442)*MStop2(1)

	Cd(4395) = 12*MStop2(2) - 7*MStop2(3)

	Cd(432) = 2*At*CA2 + MUEC*S2A

	Cd(4422) = Cd(1202)*MStop2(1) + 3*MGl2*UStop2(1,2)**2

	Cd(4424) = Cd(172)*MStop2(1) - MStop2(2)*UStop2(1,2)

	Cd(4428) = -(Cd(2432)*MStop2(1)) + 3*MStop2(2)

	Cd(4429) = MTy2*Cd(4424) - Cd(4428)*MStop2(1)

        Cd(435) = 2*Cr(28)*MStop2(1) - 
     &   2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(4430) = 1 + 6*UStop2(1,2)

	Cd(4431) = 7*MStop2(1) - Cd(4430)*MStop2(2)

	Cd(4446) = 13*MTy2 + 2*MStop2(1)

	Cd(4462) = Cd(155)*MStop2(1) + MTy2*MStop2(2)

	Cd(4467) = MGl2 + 4*MTy2

	Cd(43) = 6*A0delT + MTy2*Cd(42)

	Cd(4471) = 9*MStop2(2) - 11*MStop2(3)

	Cd(4495) = 48*MStop2(2) - 43*MStop2(3)

	Cd(4496) = 24*MStop2(2) - 37*MStop2(3)

	Cd(4498) = 72*MStop2(2) - 89*MStop2(3)

	Cd(4499) = 48*MStop2(2) - 83*MStop2(3)

	Cd(4500) = 16*MStop2(2) - MStop2(3)

	Cd(4502) = 32*MStop2(2) - 17*MStop2(3)

	Cd(4504) = 72*MStop2(2) - 5*MStop2(3)

	Cd(4505) = 36*MStop2(2) - 11*MStop2(3)

	Cd(4507) = 12*MStop2(2) - 19*MStop2(3)

	Cd(4508) = 4*MTy2 + 11*MStop2(3)

	Cd(4509) = 12*MStop2(2) - 25*MStop2(3)

	Cd(4510) = 4*MStop2(2) + 25*MStop2(3)

	Cd(4512) = 2*MStop2(2) + 9*MStop2(3)

	Cd(4513) = 8*MStop2(2) - 33*MStop2(3)

	Cd(4515) = 4*MStop2(2) - 9*MStop2(3)

	Cd(4520) = 9*MTy2 - 5*MStop2(1)

	Cd(4528) = 9*MTy2 - MStop2(1)

	Cd(4530) = 7*MTy2 + MStop2(1)

	Cd(4531) = 5*MTy2 + 3*MStop2(1)

	Cd(4561) = -MGl2 + Cd(2254)

	Cd(44) = MTy2 - MStop2(1)

	Cd(4570) = -2*MTy2*Cd(2774) + Cd(1226)*MStop2(2)

	Cd(4589) = 2*MStop2(2)**2 + MTy2*MStop2(4)

	Cd(4591) = MTy2**2 + MStop2(2)**2

	Cd(4592) = 4*MTy2*MStop2(2)**2 + 3*Cd(4591)*MStop2(4)

	Cd(4593) = -(MGl2*Cd(4589)) + Cd(4592)

	Cd(4612) = 9*MStop2(2) - 2*MStop2(3)

	Cd(4625) = -(MGl2*MTy2) + 3*Cd(4591)

	Cd(4660) = 24*MTy2 - 11*MStop2(1)

	Cd(4662) = 12*MStop2(2) - MStop2(3)

	Cd(4666) = MGlpTmSt2(2)*MStop2(2) + MTy2*Cd(251)*UStop2(1,2)**2

	Cd(45) = MTy2 + MStop2(1)

	Cd(4685) = Cd(4592) - 2*MGl2*MStop2(2)**2

	Cd(4689) = 2*(MTy2 - 3*MStop2(2))

	Cd(4700) = 3*Cd(315)*Cd(4591) + 4*MTy2*MStop2(2)**2

	Cd(463) = -7*MTy2 + 6*Cr(66)

	Cd(4715) = 12*MStop2(1) + MStop2(3)

	Cd(4730) = 17*MTy2 - 2*MStop2(1)

	Cd(4732) = 7*Cd(1496) - MStop2(1)

	Cd(4748) = 9*MTy2 - 2*MStop2(1)

	Cd(4759) = Cd(1202)*MStop2(2) + 3*MGl2*UStop2(1,2)**2

	Cd(468) = MStop2(3)**2 + 4*MTy2*MStop2(2)*UStop2(1,1)

	Cd(4760) = Cd(172)*MStop2(2) - MStop2(1)*UStop2(1,2)

	Cd(4763) = 2*MStop2(1) - Cd(169)*MStop2(2)

	Cd(4766) = 3*MStop2(1) - Cd(2432)*MStop2(2)

	Cd(4767) = MTy2*Cd(4760) - Cd(4766)*MStop2(2)

	Cd(4768) = -(Cd(4430)*MStop2(1)) + 7*MStop2(2)

	Cd(46) = -MGl2 + 2*Cd(45)

	Cd(469) = Cd(468) - 4*MTy2*MStop2(2)*UStop2(1,1)

	Cd(470) = Cd(469)*UStop2(1,1) + MStop2(3)**2*UStop2(1,2)

	Cd(471) = Cd(425) - At*AtC*Cr(28)

	Cd(4798) = MGl2 + 3*MStop2(1)

	Cd(472) = Cd(470) + Cd(471)*MStop2(2)

	Cd(4802) = 31*MTy2**2 - 3*MStop2(2)**2

	Cd(4803) = 9*MTy2**2 - Cd(2472)*MStop2(2)

	Cd(4804) = 2*MTy2*Cd(4803) + Cd(4802)*MStop2(3)

	Cd(4805) = 4*MTy2 - MStop2(2)

	Cd(473) = 8*Cd(426) - 4*Cd(424)*Cr(28)

	Cd(474) = Cd(427) - MUE2*Cr(28)

	Cd(4829) = 36*MStop2(2) + 17*MStop2(3)

	Cd(475) = S2A*Cd(473) + 8*SA2*Cd(474)

	Cd(4830) = 24*MStop2(2) + 35*MStop2(3)

	Cd(4832) = 16*MStop2(1) + MStop2(3)

	Cd(4833) = Cd(4832)*MStop2(2) + 6*MGl2*MStop2(4)

	Cd(4834) = 36*MStop2(2) - 67*MStop2(3)

	Cd(4835) = 18*MStop2(2) - 25*MStop2(3)

	Cd(4837) = 12*MStop2(2) + 13*MStop2(3)

	Cd(4838) = 24*MStop2(2) - 5*MStop2(3)

	Cd(476) = 8*CA2*Cd(472) + Cd(475)*MStop2(2)

	Cd(4841) = 4*MTy2 - 11*MStop2(3)

	Cd(4842) = 4*MStop2(2) + 5*MStop2(3)

	Cd(4843) = 12*MStop2(2) + 7*MStop2(3)

	Cd(4844) = 2*MStop2(2) - 11*MStop2(3)

	Cd(4845) = 8*MStop2(2) + 25*MStop2(3)

	Cd(4847) = 4*MStop2(2) - 29*MStop2(3)

	Cd(477) = Cr(115)*MStop2(1) - MStop2(2)*UStop2(1,1)*UStop2(1,2)

	Cd(4852) = MTy2 - 9*MStop2(2)

	Cd(478) = MStop2(1) + 2*MTy2*UStop2(1,2)

	Cd(4862) = MTy2 - 4*MStop2(2)

	Cd(4863) = 4*MTy2 + 5*MStop2(2)

	Cd(4867) = 9*MTy2 - MStop2(2)

	Cd(4868) = 7*MTy2 + MStop2(2)

	Cd(47) = Cd(44)**2 - MGl2*Cd(46)

	Cd(479) = 4*MTy2*Cd(188) + MStop2(1)*MStop2(3)

	Cd(480) = Cd(479) - 2*MTy2*MStop2(4)*UStop2(1,1)

	Cd(481) = 2*AtC*MTy*UCStop(1,3) + Cd(478)*UStop2(1,2)

	Cd(482) = Cd(481)*MStop2(3) + Cd(480)*UStop2(1,1)

	Cd(4905) = 5*MGl2 - 11*MTy2

	Cd(4906) = 6*MTy2*Cd(2617) - Cd(4905)*MStop2(1)

	Cd(4918) = MStop2(3) - 6*MStop2(4)

	Cd(4919) = 12*MStop2(2)**2 + Cd(4918)*MStop2(3)

	Cd(4920) = 20*MStop2(2) - MStop2(3)

	Cd(4921) = 20*MStop2(2)**2 - Cd(4920)*MStop2(3)

	Cd(4928) = 9*MTy2 + 7*MStop2(1)

	Cd(485) = Cr(28)*MStop2(1) - MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(4936) = 60*MStop2(2) - 11*MStop2(3)

	Cd(486) = 2*Cd(424)*Cd(485) - 2*Cd(426)*MStop2(4)

	Cd(4944) = 26*MTy2**2 + 9*MStop2(1)**2

	Cd(4948) = 3*MGl2 + 14*MTy2

	Cd(487) = (Cd(486) - 4*MTy*Cd(13)*MStop2(3))/MStop2(3)

	Cd(4950) = MGl2 - 2*Cd(2200)

	Cd(4955) = 7*MTy2 + 11*MStop2(1)

	Cd(4958) = 26*MTy2 + 5*MStop2(1)

	Cd(488) = 4*MUE2*Cd(477) - 2*Cd(427)*MStop2(4)

	Cd(5) = AtC*CA + MUE*SA

	Cd(48) = MTy2 - MStop2(2)

	Cd(4974) = 15*MGl2 + 22*MTy2

	Cd(4980) = 60*MStop2(2) - 13*MStop2(3)

	Cd(4983) = 38*MTy2 + 51*MStop2(1)

	Cd(4985) = 27*MStop2(2) - MStop2(3)

	Cd(4987) = 35*MStop2(3) + 18*MStop2(4)

	Cd(4988) = 36*MStop2(2)**2 - Cd(4987)*MStop2(3)

	Cd(491) = 2 - UStop2(1,1)

	Cd(4989) = 55*MStop2(3) + 18*MStop2(4)

	Cd(4990) = 36*MStop2(2)**2 - Cd(4989)*MStop2(3)

	Cd(4992) = 24*MStop2(2)**2 - Cd(4838)*MStop2(3)

	Cd(4998) = MStop2(3) - 18*MStop2(4)

	Cd(492) = Cd(491)*UStop2(1,1) + UStop2(1,2)**2

	Cd(4999) = 36*MStop2(2)**2 + Cd(4998)*MStop2(3)

	Cd(5005) = -MGl2 + 3*Cd(198)

	Cd(493) = MTy2*Cd(492) + AtC*MTy*UCStop(1,3)

	Cd(5016) = 7*MTy2 - 10*MStop2(1)

	Cd(5017) = 2*MTy2 - 5*MStop2(1)

	Cd(5020) = 12*MTy2 + 7*MStop2(1)

	Cd(5023) = 7*MGlmT2 - 17*MStop2(1)

	Cd(5029) = 11*MTy2 + 24*MStop2(1)

	Cd(5034) = 3*MTy2 + 11*MStop2(1)

	Cd(496) = Cd(427) + 2*MUE2*UStop2(1,1)*UStop2(1,2)

	Cd(5053) = 72*MStop2(2) - 39*MStop2(3)

	Cd(5057) = 12*MTy2 + 17*MStop2(1)

	Cd(5058) = 5*MStop2(2) - 2*MStop2(3)

	Cd(498) = 2*MTy + AtC*UCStop(1,3)

	Cd(5062) = 2*MStop2(1)**2 + MStop2(3)*MStop2(4)

	Cd(49) = MTy2 + MStop2(2)

	Cd(500) = MTy2*UStop(1,2) + AtC*MTy*UStop(1,1)*UStop2(1,2)

	Cd(5080) = 14*MStop2(1) + MStop2(3)

	Cd(5081) = 21*MStop2(1)**2 + Cd(5080)*MStop2(3)

	Cd(5087) = 7*MTy2 - 27*MStop2(1)

	Cd(5096) = 20*MTy2 + 7*MStop2(1)

	Cd(5104) = 30*MStop2(2) - 7*MStop2(3)

	Cd(5105) = 5*MTy2 + 3*MStop2(4)

	Cd(5108) = 8*MStop2(3) + 3*MStop2(4)

	Cd(5109) = 6*MStop2(2)**2 - Cd(5108)*MStop2(3)

	Cd(5110) = 3*MStop2(2)**2 - Cd(1074)*MStop2(3)

	Cd(5112) = 3*MTy2**2 + Cd(5110)

	Cd(505) = Cr(85) - UStop2(1,1)

	Cd(506) = Cd(505)*UStop2(1,1)**3 + 2*UStop2(1,1)*UStop2(1,2)**3

	Cd(5141) = 3*MStop2(2)**2 - Cd(1022)*MStop2(3)

	Cd(5142) = MTy2 + Cr(84)

	Cd(5143) = Cd(5141) + MTy2*Cd(5142)

	Cd(5144) = 5*MStop2(2) + MStop2(3)

	Cd(5146) = 5*MStop2(3) + 4*MStop2(4)

	Cd(5147) = 7*MStop2(2)**2 - 2*Cd(5146)*MStop2(3)

        Cd(507) = Cd(506)*MStop2(1) - 
     &   2*MStop2(3)*UStop2(1,1)**2*UStop2(1,2)**2

	Cd(5150) = 5*MStop2(2)**2 - 6*MStop2(3)*MStop2(4)

	Cd(5152) = 5*MStop2(2) - MStop2(3)

	Cd(5154) = 8*MTy2 + MStop2(1)

	Cd(5158) = Cd(1657)*MStop2(2)**2 + MStop2(3)**3

	Cd(508) = 1 - 3*UStop2(1,2)

	Cd(5168) = 49*MStop2(2) - 30*MStop2(3)

	Cd(509) = Cd(508)*MStop2(1) + 2*MStop2(2)*UStop2(1,2)

	Cd(5170) = 7*MStop2(2) - 16*MStop2(3)

	Cd(5171) = Cd(5170)*MStop2(2)**2 + Cd(1570)*MStop2(3)**2

	Cd(5172) = 41*MStop2(2) - 7*MStop2(3)

	Cd(5173) = 49*MStop2(2)**2 - Cd(5172)*MStop2(3)

	Cd(5176) = 31*MStop2(2) - 2*MStop2(3)

	Cd(5180) = 5*MStop2(3) - 6*MStop2(4)

	Cd(5182) = 7*MStop2(3) + 15*MStop2(4)

	Cd(5183) = 49*MStop2(2)**2 - 2*Cd(5182)*MStop2(3)

	Cd(5184) = 7*MStop2(1) + MStop2(3)

	Cd(5185) = 35*MStop2(2) - 66*MStop2(3)

	Cd(5186) = Cd(5185)*MStop2(2)**2 + 3*Cd(1764)*MStop2(3)**2

	Cd(511) = Cr(1) - MStop2(1)

	Cd(5189) = 63*MStop2(2) - 34*MStop2(3)

	Cd(5190) = 19*MStop2(2) - 14*MStop2(3)

	Cd(5191) = 49*MStop2(1)**2 + 3*Cd(5190)*MStop2(3)

	Cd(5192) = 7*MStop2(3) + 31*MStop2(4)

	Cd(5193) = 105*MStop2(2)**2 - 2*Cd(5192)*MStop2(3)

	Cd(5195) = 63*MStop2(2) - 118*MStop2(3)

	Cd(5196) = 55*MStop2(3) + 24*MStop2(4)

	Cd(5197) = Cd(5195)*MStop2(2)**2 + Cd(5196)*MStop2(3)**2

	Cd(512) = MStop2(1) - MStop2(2)*UStop2(1,1)

	Cd(5199) = 61*MStop2(3) + 34*MStop2(4)

	Cd(5202) = 126*MStop2(2) - 76*MStop2(3)

	Cd(5207) = 2*MStop2(3) - 3*MStop2(4)

	Cd(513) = 2*Cd(512)*MStop2(1) - 2*MStop2(2)**2*UStop2(1,2)

	Cd(5211) = 3*MStop2(2) + 2*MStop2(3)

	Cd(5214) = MStop2(2)**2 - 2*MTy2*MStop2(3)

	Cd(5215) = MStop2(2)**2 + 2*MStop2(3)*MStop2(4)

	Cd(514) = Cd(513) + 2*MStop2(1)*MStop2(3)*UStop2(1,1)

	Cd(5224) = 2*MTy2 + 3*MStop2(4)

	Cd(5225) = MTy2 + MStop2(4)

	Cd(5227) = -2*MGl2 + Cd(5224)

	Cd(5230) = 13*MTy2 + 12*MStop2(1)

	Cd(5233) = 14*MTy2 + 27*MStop2(1)

	Cd(5234) = 20*MStop2(2) + 9*MStop2(3)

	Cd(5243) = 24*MTy2 + 19*MStop2(1)

	Cd(518) = 2*MUE2*Cd(485) - Cd(427)*MStop2(4)

	Cd(519) = 1 + 3*UStop2(1,2)

	Cd(5269) = 17*MTy2 + 23*MStop2(1)

	Cd(5274) = 20*MTy2 + 9*MStop2(1)

	Cd(520) = Cd(519) - UStop2(1,1)

	Cd(5287) = 90*MStop2(2) - 19*MStop2(3)

	Cd(521) = Cr(85)*MStop2(1) - Cr(1)*UStop2(1,1)

	Cd(5290) = 5*MTy2 + 13*MStop2(1)

	Cd(5294) = 6*MTy2 + 5*MStop2(1)

        Cd(522) = Cd(520)*MStop2(1)*UStop2(1,1)**3 + 
     &   2*Cd(521)*UStop2(1,1)*UStop2(1,2)**2

	Cd(523) = 2*MStop2(2)**2 - MStop2(3)*MStop2(4)

	Cd(5308) = -MGl2 + Cd(5224)

	Cd(524) = 4*Cd(512)*MStop2(1) - 2*Cd(523)*UStop2(1,2)

	Cd(5320) = 2*MStop2(3) + 7*MStop2(4)

	Cd(5321) = 14*MStop2(2)**2 - Cd(5320)*MStop2(3)

	Cd(5325) = 11*MTy2 - 4*MStop2(1)

	Cd(525) = Cd(524) + 4*MStop2(1)*MStop2(3)*UStop2(1,1)

	Cd(526) = MStop2(1) + Cr(1)*UStop2(1,2)

	Cd(5340) = 9*MGl2 + 65*MTy2

	Cd(5343) = MGl2 + 5*MTy2

	Cd(527) = UStop2(1,1)**2 + 4*UStop2(1,2)

	Cd(5349) = 10*MTy2 + 3*MStop2(4)

	Cd(5355) = 22*MStop2(2) - 3*MStop2(3)

	Cd(5356) = 22*MStop2(2)**2 - Cd(5355)*MStop2(3)

	Cd(5357) = MStop2(2) + 4*MStop2(3)

	Cd(5358) = MStop2(1)**2 + Cd(5357)*MStop2(3)

	Cd(5359) = -7*MTy2**2 + Cd(5358)

	Cd(5361) = 8*MStop2(3) - 5*MStop2(4)

	Cd(5362) = 20*MStop2(2)**2 + 2*Cd(5361)*MStop2(3)

	Cd(51) = 3 - 2*Log(MStop2(1)/Q)

	Cd(529) = 2*MTy2 + MStop2(1)

	Cd(5370) = 3*MTy2 + MStop2(4)

	Cd(5371) = 3*Cd(523) + 2*MTy2*Cd(5370)

	Cd(5373) = -MGl2 + 2*Cd(198)

	Cd(5374) = Cd(5371) - MGl2*Cd(5373)

	Cd(5375) = MTy2 - 2*MStop2(4)

	Cd(5376) = 3*Cd(523) + MTy2*Cd(5375)

	Cd(5384) = 33*MStop2(2) - 26*MStop2(3)

	Cd(5385) = 33*MStop2(1)**2 + Cd(5384)*MStop2(3)

	Cd(5386) = 75*MStop2(2) - 62*MStop2(3)

	Cd(5387) = 75*MStop2(1)**2 + Cd(5386)*MStop2(3)

        Cd(531) = AtC*Cd(522) - 
     &   2*At*MStop2(4)*UCStopC(1,3)**2*UStop2(1,1)*UStop2(1,2)

	Cd(5388) = 44*MStop2(2) - 3*MStop2(3)

	Cd(5391) = 84*MTy2 + 34*MStop2(4)

	Cd(5392) = 78*MStop2(2) - 19*MStop2(3)

	Cd(5394) = 28*MStop2(2) - 3*MStop2(3)

	Cd(5397) = 3*MStop2(1)**2 + Cd(286)*MStop2(3)

	Cd(5399) = 28*MStop2(3) + 51*MStop2(4)

	Cd(5400) = 102*MStop2(2)**2 - Cd(5399)*MStop2(3)

	end
