* TLhr_abbr1_16.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_16
	implicit none

#include "TLhr.h"


        Cd(9095) = 3*Cd(9094)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(4565)*UStop(1,2)

        Cd(9108) = M_3*MTy*Cd(9107)*UUStop(1,2) + 
     &   2*Cd(50)**2*MStop2(2)*UUStop(2,3)

	Cd(9109) = Cd(9105) + 4*MTy*MTy2*Cd(9095)*UStop2(1,2)**2

	Cd(9128) = Cd(44)**2*Cd(1561) - MGl2*Cd(2143)

	Cd(9153) = Cd(52)*Cd(9152) - 2*Cd(1493)*Cd(9144)*MStop2(3)

	Cd(9159) = Cd(9155)*MGlpTmSt4(2) + Cd(96)*Cd(9158)*UStop(1,1)

	Cd(9160) = Cd(9159)*MStop2(3) + M_3*MTy*Cd(2265)*UStop(1,2)

	Cd(9162) = 2*Cd(9155)*MGlpTmSt4(2) + Cd(96)*Cd(9158)*UStop(1,1)

	Cd(9163) = Cd(9162)*MStop2(3) + M_3*MTy*Cd(2265)*UStop(1,2)

	Cd(908) = MTy*Cd(886)*UCStopC(1,3) + AtC*Cd(871)*UStop2(1,1)

	Cd(9194) = 3*Cd(9190)*UStop(1,1) - M_3*MTy*Cd(9193)*UStop(1,2)

	Cd(9248) = 4*Cd(44)*Cd(9243) + MGl2*Cd(9247)

	Cd(9249) = Cd(44)**3*Cd(2986) - MGl2*Cd(9248)

        Cd(9253) = Cd(9249)*UStop(1,1) - 
     &   M_3*MTy*Cd(9252)*MStop2(3)*UStop(1,2)

	Cd(9263) = Cd(9231) + MGl2*Cd(9262)

	Cd(9264) = Cd(44)**2*Cd(9227) - MGl2*Cd(9263)

	Cd(9265) = Cd(9253)*MStop2(3) - MTy2*Cd(9240)*UStop2(1,2)**2

        Cd(9266) = M_3C*Cd(9261)*UUStop(1,1) - 
     &   M_3*Cd(9264)*UUStop(1,2)

        Cd(9267) = Cd(9265)*UStopC(1,1) + 
     &   MTy*Cd(9266)*MStop2(3)*UUStopC(1,3)

	Cd(9277) = -(MGl2*Cd(1602)) + Cd(9272)

	Cd(9278) = Cd(44)**3*Cd(89) - MGl2*Cd(9277)

	Cd(9279) = Cd(9276)*UStop(1,1) - 2*M_3*MTy*Cd(9278)*UStop(1,2)

        Cd(9280) = 3*Cd(1791)*UStop(1,1) - 
     &   8*M_3*MGl2*MTy*MTy2*MStop2(3)**2*UStop(1,2)

	Cd(923) = Cd(922)*Cr(13) - Cd(917)*MStop2(3)

        Cd(9323) = -(Cd(118)**2*Cd(9322)) + 
     &   4*Cd(9280)*UStop(1,1)*UStop2(1,2)

	Cd(924) = A0delStop(1)*Cd(915) - Cd(923)*MStop2(3)

        Cd(925) = 3/(32.D0*Pi)*
     &   (Alfa2L*(Cd(924) + 
     &        2*A0delStop(2)*Cd(359)*Cd(360)*MStop2(1)))/
     &    (SB2*SW2*MStop2(1))

	Cd(926) = Cd(915)*Cr(13)*MStop2(1) - Cd(922)*MStop2(3)

        Cd(927) = 3/(16.D0*Pi)*
     &   (Alfa2L*(Cd(926) + 
     &        2*Cd(359)*Cd(360)*Cr(34)*MStop2(1)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(3)**2)

	Cd(9394) = Cd(9391) - MGl2*Cd(9393)

	Cd(9395) = Cd(44)**2*Cd(9384) - MGl2*Cd(9394)

	Cd(9397) = Cd(9396)*UCStop(1,3) + Cd(9395)*UUStopC(1,2)

        Cd(9404) = 3*Cd(9403)*MStop2(3)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(1793)*UStop(1,2)

	Cd(9411) = 2*Cd(44)*Cd(9409) - MGl2*Cd(9410)

	Cd(9412) = Cd(44)**3*Cd(3736) + MGl2*Cd(9411)

        Cd(9415) = 2*Cd(9412)*UStop(1,1) + 
     &   M_3*MTy*Cd(9414)*MStop2(3)*UStop(1,2)

	Cd(9431) = Cd(44)**2*Cd(9422) - MGl2*Cd(9430)

	Cd(9437) = Cd(9433) + MGl2*Cd(9436)

	Cd(9438) = Cd(9434) - MGl2*Cd(9437)

	Cd(9439) = Cd(44)**2*Cd(9432) - MGl2*Cd(9438)

        Cd(9440) = M_3*MTy*Cd(9439) - 
     &   Cd(47)**2*MGlpTmSt2(1)*MStop2(3)*UUStop(2,3)

	Cd(9441) = Cd(9415)*MStop2(3) + MTy2*Cd(9404)*UStop2(1,2)**2

	Cd(9442) = -(M_3C*MTy*Cd(9431)) + Cd(9440)

	Cd(9443) = Cd(9441)*UStopC(1,1) - Cd(9442)*UUStopC(1,3)

        Cd(9450) = Cd(9448)*UUStop(1,1) + 
     &   4*M_3*MTy*Cd(9449)*UUStop(1,3)

        Cd(9468) = 3*Cd(96)*Cd(1708)*MStop2(3)*UStop(1,1) + 
     &   M_3*MTy*Cd(9467)*UStop(1,2)

        Cd(9469) = 3*Cd(1776)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(50)*MStop2(3)*UStop(1,2)

        Cd(9478) = 4*Cd(9469)*MStop2(3)*UStop(1,1)*UStop2(1,2) + 
     &   2*Cd(9477)*UUStop(1,1)

        Cd(9523) = Cd(9519)*UStop(1,1) + 
     &   M_3*MTy*Cd(9522)*MStop2(3)*UStop(1,2)

	Cd(9531) = Cd(9529) + MGl2*Cd(9530)

	Cd(9534) = Cd(9502) + MGl2*Cd(9533)

	Cd(9535) = Cd(48)**2*Cd(9497) - MGl2*Cd(9534)

	Cd(9536) = Cd(9523)*MStop2(3) + MTy2*Cd(9509)*UStop2(1,2)**2

	Cd(9660) = Cd(9658) - MGl2*Cd(9659)

	Cd(9661) = Cd(48)**2*Cd(9648) - MGl2*Cd(9660)

	Cd(9663) = Cd(9662)*UCStop(1,3) + Cd(9661)*UUStopC(1,2)

        Cd(9670) = 3*Cd(9669)*MStop2(3)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(1745)*UStop(1,2)

	Cd(9679) = 2*Cd(48)*Cd(9674) - MGl2*Cd(9678)

	Cd(9680) = Cd(48)**3*Cd(3989) + MGl2*Cd(9679)

        Cd(9683) = 2*Cd(9680)*UStop(1,1) + 
     &   M_3*MTy*Cd(9682)*MStop2(3)*UStop(1,2)

	Cd(9700) = Cd(48)**2*Cd(9690) - MGl2*Cd(9699)

	Cd(9705) = Cd(9702) + MGl2*Cd(9704)

	Cd(9706) = Cd(9703) - MGl2*Cd(9705)

	Cd(9707) = Cd(48)**2*Cd(9701) - MGl2*Cd(9706)

        Cd(9708) = M_3*MTy*Cd(9707) - 
     &   Cd(50)**2*MGlpTmSt2(2)*MStop2(3)*UUStop(2,3)

	Cd(9709) = Cd(9683)*MStop2(3) + MTy2*Cd(9670)*UStop2(1,2)**2

	Cd(9710) = -(M_3C*MTy*Cd(9700)) + Cd(9708)

	Cd(9711) = Cd(9709)*UStopC(1,1) - Cd(9710)*UUStopC(1,3)

	Cd(9785) = Cd(424)*Cd(9777)*MStop2(3) + Cd(9784)*UStop2(1,1)

	Cd(9789) = CA2*Cd(546) + SA2*Cd(560)*MStop2(3)**2

        Cd(9790) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(9789) + S2A*Cd(9788)*MStop2(1)))/
     &    (SB2*SW2*MStop2(1)**2*MStop2(3)**2)

	Cd(9800) = Cd(402)*Cd(9739)*MStop2(1) - Cd(9799)*MStop2(2)

	Cd(9806) = 2*Cd(9800) - Cd(9805)

	Cd(9813) = Cd(9812)*MStop2(1) - Cd(403)*Cd(9741)*MStop2(2)

	Cd(9819) = 2*Cd(9813) + Cd(9818)

        Cd(9820) = A0delStop(1)*Cd(9800) + 
     &   Cd(9799)*Cr(13)*MStop2(1)*MStop2(3)

        Cd(9822) = 3/(16.D0*Pi)*
     &   (Alfa2L*(-(Cd(9821)*MStop2(1)) + Cd(9820)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(2)*MStop2(3)**2)

	Cd(9823) = Cd(9800)*Cr(13) + Cd(9799)*MStop2(3)

	Cd(9825) = 2*Cd(9823) - Cd(9824)

        Cd(9826) = A0delStop(2)*Cd(9813) + 
     &   Cd(9812)*Cr(34)*MStop2(2)*MStop2(3)

        Cd(9828) = 3/(16.D0*Pi)*
     &   (Alfa2L*(Cd(9826)*MStop2(1) + Cd(9827)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(2)*MStop2(3)**2)

	Cd(9829) = Cd(9813)*Cr(34) + Cd(9812)*MStop2(3)

	Cd(9831) = 2*Cd(9829) + Cd(9830)

	Cd(9848) = -2*SA2*Cd(819) + S2A*Cd(9847)

        Cd(9849) = 1/(128.D0*Pi**3)*
     &   (Alfas2L*(Cd(9848) + 
     &        4*CA2*Cd(13)**2*Cr(75)**2*MStop2(3)**2))/
     &    (SB2*SW2*MStop2(3)**2)

	Cd(9850) = S2A*Cd(823) - 4*CA2*Cd(13)**2*Cr(75)*MStop2(4)

        Cd(9851) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*Cr(9)*(-Cd(9850) + SA2*Cd(828)*MStop2(2)))/
     &    (SB2*SW2*MStop2(2))

	Cd(9868) = MUE*Cd(9860)*UCStop(1,3) + MUEC*Cd(9865)*UCStopC(1,3)

        Cd(9869) = 4*MUEC*Cd(9867)*UCStopC(1,3)**2 + 
     &   MUE*Cd(9855)*UStop2(1,1)*UStop2(1,2)

	Cd(9924) = Cd(9923)*Cr(13) - Cd(9918)*MStop2(3)

	Cd(9925) = A0delStop(1)*Cd(9915) - Cd(9924)*MStop2(3)

        Cd(9926) = 3/(32.D0*Pi)*
     &   (Alfa2L*(Cd(9925) + 
     &        2*A0delStop(2)*Cd(402)*Cd(403)*MStop2(1)))/
     &    (SB2*SW2*MStop2(1))

	Cd(9927) = Cd(9915)*Cr(13)*MStop2(1) - Cd(9923)*MStop2(3)

        Cd(9928) = 3/(16.D0*Pi)*
     &   (Alfa2L*(Cd(9927) + 
     &        2*Cd(402)*Cd(403)*Cr(34)*MStop2(1)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(3)**2)

        Cd(9977) = (-(Cd(9974)*MStop2(1)) + 
     &     Cd(50)*Cd(9950)*
     &      (Cd(9968)*T134(MGl2,MTy2,MStop2(1)) - 
     &        Cd(9976)*T134(MGl2,MTy2,MStop2(2))))/
     &   (Cd(50)*Cd(9950))

        Cd(9987) = -((Cd(9980)*Cd(9982))/
     &      (Cd(47)**2*Cd(50)**2)) + 
     &   (At*Cd(47)*Cd(50)*Cd(9955)*Cd(9979) + Cd(9986)*MStop2(3))/
     &    (Cd(47)*Cd(50)*MStop2(3)**2)

        Cd(9991) = -((Cd(9980)*Cd(9988))/
     &      (Cd(47)**2*Cd(50)**2)) + 
     &   (AtC*Cd(47)*Cd(50)*Cd(9955)*Cd(9979) + 
     &      Cd(9990)*MStop2(3))/(Cd(47)*Cd(50)*MStop2(3)**2)

        Cd(10007) = (-(Cd(10004)*MStop2(1)) + 
     &     Cd(50)*Cd(9950)*
     &      (Cd(10002)*T134(MGl2,MTy2,MStop2(1)) - 
     &        Cd(10006)*T134(MGl2,MTy2,MStop2(2))))/
     &   (Cd(50)*Cd(9950))

        Cd(10010) = 2*Cd(47)*Cd(9955)*Cd(10009) + 
     &   Cd(50)*(-Cd(9951) + Cd(9953))*MStop2(3)

        Cd(10011) = Cd(50)*Cd(10010) + 
     &   Cd(47)**2*(Cd(9952) - Cd(9954))*MStop2(3)

        Cd(10024) = At*Cd(10023)*Cr(5) + 
     &   UCStop(1,3)*(AtC*Cd(10022) + 
     &      Cd(50)*Cd(9947)*Cd(9955)*MGlpTmSt2(1)*MStop2(1)*
     &       UStop2(1,2))

        Cd(10034) = Cd(10031)*Cd(10033) + 
     &   Cd(47)*Cd(50)*(Cd(50)*Cd(9947)*Cd(10026) - 
     &      Cd(47)*Cd(9949)*Cd(10030))*UCStop(1,3)

        Cd(10040) = Cd(47)*Cd(50)*
     &    (Cd(50)*Cd(9947)*Cd(10026) - Cd(47)*Cd(9949)*Cd(10038))+
     &     At*Cd(10031)*Cd(10039)

        Cd(10051) = Cd(47)**2*(Cd(9952) - Cd(9954)) + 
     &   Cd(50)**2*Cd(10012)

        Cd(10053) = AtC*Cd(9980)*Cd(10049) + 
     &   At*Cd(9980)*Cd(10050) - 4*Cd(10051)*Cd(10052)

        Cd(10067) = Cd(47)*Cd(50)*
     &    (-(Cd(50)*Cd(9947)*Cd(10062)) + 
     &      Cd(47)*Cd(9949)*Cd(10066)) + Cd(10053)*MStop2(3)

        Cd(10068) = Cd(10067)*MStop2(3) - 
     &   Cd(47)**2*Cd(50)**2*Cd(424)*Cd(9955)*
     &    (Cr(28)*T134(MGl2,MTy2,MStop2(2)) + 
     &      T134(MGl2,MTy2,MStop2(1))*UStop2(1,2)**2)

        Cd(10082) = (2*Cd(47)*Cd(50)*Cd(10081) + 
     &     Cd(9980)*MStop2(3))/MStop2(3)

        Cd(10083) = 2*Cd(47)*Cd(50)*
     &    (Cd(50)*Cd(9947)*Cd(10077)*MGlpTmSt2(1) + 
     &      Cd(47)*Cd(9949)*MGlpTmSt2(2)) + Cd(9980)*MStop2(3)

        Cd(10090) = (Cd(10031)*Cd(10086))/
     &    (Cd(47)**2*Cd(50)**2) + 
     &   (2*(-2*Cd(47)*Cd(50)*Cd(10014)*Cd(10085) + 
     &        Cd(10089)*MStop2(3)))/(Cd(47)*Cd(50)*MStop2(3)**2)

        Cd(10098) = (Cd(10031)*Cd(10039))/
     &    (Cd(47)**2*Cd(50)**2) + 
     &   (2*(-2*Cd(47)*Cd(50)*Cd(10014)*Cd(10094) + 
     &        Cd(10097)*MStop2(3)))/(Cd(47)*Cd(50)*MStop2(3)**2)

        Cd(10105) = 2*Cd(50)**2*Cd(9953) + 
     &   2*Cd(47)**2*Cd(9954) - Cd(9980)*UStop2(1,1)

        Cd(10106) = Cd(50)**2*Cd(10101) + 4*Cd(9948)*Cd(10104) + 
     &   Cd(10105)*UStop2(1,1)

        Cd(10108) = 2*Cd(50)**2*Cd(9953) + 
     &   2*Cd(47)**2*Cd(9954) - 3*Cd(9980)*UStop2(1,1)

        Cd(10120) = AtC*Cd(10117) + 
     &   At*(Cd(10119) + At*Cd(50)*Cd(9947)*MGlpTmSt2(1)*
     &       UCStopC(1,3)**2)

        Cd(10126) = AtC**2*Cd(10081)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(10080) + At*Cd(10125)*UCStopC(1,3)**2)

        Cd(10129) = AtC**2*Cd(10128)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(10071) + 
     &      At*Cd(10076)*Cd(10077)*UCStopC(1,3)**2)

        Cd(10134) = (-2*Cd(47)*Cd(9955)*Cd(10009) + 
     &     Cd(50)*Cd(10012)*MStop2(3))/(Cd(50)*MStop2(3))

        Cd(10139) = Cd(47)**2*(Cd(9952) - Cd(9954)) - 
     &   Cd(50)**2*Cd(10012)

        Cd(10154) = 2*Cd(47)*Cd(50)*
     &    (Cd(50)*Cd(9947)*Cd(10150) - Cd(47)*Cd(9949)*Cd(10152))+
     &     Cd(10031)*Cd(10153)

        Cd(10155) = AtC**2*Cd(10089)*UCStop(1,3) + 
     &   At*(2*AtC*Cd(10093)*Cr(5) + At*Cd(10097)*UCStopC(1,3))

        Cd(10185) = At*Cd(10184) + 
     &   AtC*(Cd(10175) + AtC*Cd(10081)*UCStop(1,3)**2)

        Cd(10215) = MTy2*(4*MTy2*Cd(10194) + Cd(10214)) + 
     &   (2*MTy*Cd(10203))/MStop2(1)**2

        Cd(10229) = Cd(152)*Cd(2597)*
     &    (AtC*Cd(10221) + At*Cd(10225)) + 2*Cd(10228)

        Cd(10246) = MTy2*
     &    (Cd(9948)*Cd(10240) + Cd(251)*Cd(10241)) + 
     &   Cd(10245)*MGlpTmSt2(2)

        Cd(1018) = 2*Cd(1016)*UStop2(1,1)**2 + 
     &   Cd(1017)*UStop2(1,2) - 4*MStop2(1)**2*UStop2(1,2)**4

        Cd(10270) = -(Cd(10259)*
     &      (MUE**2*UCStop(1,3)**2 + MUEC**2*UCStopC(1,3)**2)) + 
     &   2*MUE2*Cd(10269)*UStop2(1,1)

        Cd(10406) = At*Cd(10397) + 
     &   2*Cd(47)*Cd(10405)*MStop2(1) + Cd(10394)*MStop2(3)

        Cd(10433) = 2*Cd(47)*
     &    (Cd(13)*Cd(47)*Cd(10407) + Cd(10432)*MStop2(1)) + 
     &   Cd(10412)*MStop2(3)

        Cd(10473) = (-4*MGl2 - Cd(10472) + 
     &     4*Cd(10471)*MStop2(1))/MStop2(1)

        Cd(10545) = At*Cd(10542) + Cd(47)*Cd(10544) - 
     &   2*Cd(164)*Cd(10540)*MStop2(3)

        Cd(10578) = 8*(MUE2*Cd(10576) + 
     &     MUE**2*Cd(10577)*UCStop(1,3)**2 + 
     &     MUEC**2*Cd(10577)*UCStopC(1,3)**2)

        Cd(10612) = -(AtC**2*Cd(10610)*UCStop(1,3)**2) + 
     &   At*(AtC*Cd(10609) - At*Cd(10611)*UCStopC(1,3)**2)

        Cd(10629) = Cd(50)*Cd(10309) + 
     &   AtC*Cd(10628)*UCStop(1,3) + At*Cd(10602)*UCStopC(1,3)

        Cd(10638) = AtC*(Cd(10636) + 
     &      AtC*Cd(10086)*Cd(10604)*UCStop(1,3)) + 
     &   At*Cd(10637)*UCStopC(1,3)

        Cd(10641) = -2*Cd(251)*Cd(9969) + 
     &   (Cd(10604) + 2*Cd(10640)*MStop2(2))/MStop2(2)**2

        Cd(10643) = -(Cd(50)*Cd(10342)) + 4*Cd(10602)*Cr(28) + 
     &   2*Cd(9948)*Cd(10642)*UStop2(1,1)

        Cd(10662) = Cd(9948)*Cd(10656) + 
     &   MGlmT2*Cd(50)*Cd(10661) + Cd(10657)*MStop2(3)

        Cd(10685) = Cd(50)*
     &    (2*Cd(50)*Cd(10667) - 2*Cd(10684)*MStop2(2)) + 
     &   Cd(10673)*MStop2(3)

        Cd(10701) = Cd(50)*
     &    (2*Cd(50)*Cd(10686) - 2*Cd(10700)*MStop2(2)) + 
     &   Cd(10691)*MStop2(3)

        Cd(10733) = Cd(10732)*UStop2(1,1)**2 + 
     &   Cd(50)*MStop2(2)*(MStop2(3) - Cd(152)*UStop2(1,1)**3)

        Cd(10752) = Cd(9948)*Cd(10750) + 
     &   4*Cd(50)*Cd(179)*Cd(10751)*MStop2(2) + 
     &   Cd(10448)*Cd(10603)*MStop2(3)

        Cd(10758) = 4*Cd(50)*(Cd(9948)*Cd(10755) - Cd(10757))*
     &    MStop2(2) + Cd(10753)*MStop2(3)

	Cd(10768) = MTy2*Cd(10752) + MGl2*Cd(10758) + 4*MTy*Cd(10767)

        Cd(10774) = (Cd(10612) + Cd(50)*Cd(492)*MStop2(3) - 
     &     MGlmT2*Cd(10618)*MStop2(3))/MStop2(3)

        Cd(10776) = (-4*MGl2 - Cd(10775) + 
     &     4*Cd(10471)*MStop2(2))/MStop2(2)

	Cd(10814) = At*Cd(10804) + Cd(9948)*Cd(10808) - 2*Cd(10813)

	Cd(10822) = AtC*Cd(10804) + Cd(9948)*Cd(10815) - 2*Cd(10821)

        Cd(10870) = 2*MGl2*Cd(251)*Cd(10859) + 
     &   (MUE*Cd(10858) + MUEC*Cd(10869))/MStop2(3)

        Cd(10885) = -(MGl2*Cd(10757)) + MGlmT2*Cd(10884) - 
     &   Cd(10881)*MStop2(2)

        Cd(10975) = (Cd(10973) - 
     &     MUE**2*MStop2(3)**2*UCStop(1,3)**2*UStop2(1,2)**2 + 
     &     Cd(10974)*MStop2(3)**2*UUStopC(1,2))/MStop2(3)**2

	Cd(11006) = Cd(11007) + 5*Cd(10992)*MStop2(3)**2

	Cd(11193) = 4*AtC*Cd(11179) + 4*MUE*Cd(11185) - 2*MUEC*Cd(11192)

	Cd(11213) = Cd(11210) - 2*Cd(11212) + Cd(11208)*UUStopC(1,1)

        Cd(11214) = (Cd(11213) + 
     &     4*Cd(9794)*MStop2(3)*UStop2(1,1)*
     &      (-Cd(11200) + Cd(9794)*MStop2(3)*UStop2(1,1))*
     &      UUStopC(1,1))/MStop2(3)**2

	Cd(11228) = 2*Cd(11218) - Cd(11220) - 4*Cd(11227)*UUStopC(1,2)

        Cd(11229) = (-Cd(11228) + 
     &     4*MUE*MStop2(3)*UCStop(1,3)*UStop2(1,2)*
     &      (MUE*MStop2(3)*UCStop(1,3)*UStop2(1,2) + 
     &        Cd(11217)*UStopC(1,2)))/(MStop2(1)**2*MStop2(3)**2)

	Cd(11252) = Cd(11074) + Cd(11079) + AtC*Cd(11100)

	Cd(11255) = 2*Cd(11128) - 2*Cd(11131) - AtC*Cd(11148)

	Cd(1119) = AtC*Cd(1116) + At*Cd(1118) - 2*MTy*Cd(1117)*MStop2(3)

	Cd(1130) = -Cd(957) + Cd(1018) - Cd(1127)*MStop2(3)

	Cd(1132) = Cd(1030) + Cd(1129)*MStop2(3)

        Cd(1133) = 2*MUE2*Cd(1130) - 
     &   MUE**2*Cd(1131)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(1132)*UCStopC(1,3)**2

        Cd(11421) = 2*MStop2(2)*
     &    (Cd(11419) - 2*Cd(11420)*MStop2(2)) - 
     &   4*Cd(12)*Cd(13)*Cd(11375)*MStop2(3)**2

        Cd(11435) = Cd(11289) - Cd(11432) - 
     &   MUE**2*Cr(91)*MStop2(3)*UCStop(1,3)*UStop2(1,1)

        Cd(11436) = UStopC(1,2)*
     &    (MUEC*Cd(9948)*Cd(11433) - MUE*Cd(11434)*UUStopC(2,2)) + 
     &   Cd(11435)*UUStopC(3,2)

	Cd(11523) = 4*Cd(11322) - 4*Cd(11466) + AtC*Cd(11522)

        Cd(11525) = 2*Cd(11347) - Cd(11489) + 
     &   Cd(9809)*Cd(11524)*MStop2(3)*UStop2(1,2)

	Cd(11551) = Cd(11418) - Cd(11548) + Cd(11550)*MStop2(3)

	Cd(11557) = Cd(11438) + Cd(11441)

	Cd(11567) = Cd(11289) - Cd(11432) + Cd(11429)*MStop2(3)

        Cd(11569) = MUEC*Cd(9948)*Cd(11566)*UStopC(1,2) - 
     &   MUE*Cd(11568)*UUStopC(2,2) + Cd(11567)*UUStopC(3,2)

	Cd(11605) = 4*Cd(11347) - 2*Cd(11489) - Cd(11604)*MStop2(3)

        Cd(11651) = AtC*Cd(11645) + At*Cd(11650) - 
     &   2*MTy*Cd(13)*Cd(11647)*MStop2(3)

        Cd(11697) = (4*MTy*Cd(585)*Cd(11625) + 
     &     (2*Cd(11696) - 4*Cd(10191)*Cr(75)**2*MStop2(1))/
     &      (MStop2(1)*MStop2(2)))/MStop2(3)

	Cd(11711) = -(SA2*Cd(11697)) + S2A*Cd(11704) + CA2*Cd(11710)

	Cd(1168) = Cd(1083) - Cd(1127)*MStop2(3)

        Cd(11769) = MGlmT2*Cd(11768) + 
     &   M_3C*Cd(11763)*MStop2(1) + Cd(11762)*UCStopC(1,3)

	Cd(1170) = -4*Cd(984) + 2*Cd(1068) - Cd(1152)

	Cd(1171) = Cd(1088) + Cd(1162)

	Cd(1172) = Cd(1089) + Cd(1166)

        Cd(11807) = -2*MGl2*Cd(152)*Cd(11801) + 
     &   (Cd(47)*(MUE*Cd(11799) + MUEC*Cd(11806)))/MStop2(3)

        Cd(11822) = Cd(47)*Cd(11821) + 
     &   MGl2*Cd(152)*(Cd(13)*Cd(152) + 2*Cd(11808))*MStop2(3)

        Cd(1178) = 2*Cd(1175)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(1177)*UStop2(1,1) - 
     &      4*MStop2(1)*MStop2(3)*UStop2(1,2))

        Cd(11856) = -(AtC**2*Cd(11854)*UCStop(1,3)**2) + 
     &   At*(AtC*Cd(11853) - At*Cd(11855)*UCStopC(1,3)**2)

        Cd(11901) = At*Cd(11895) + Cd(251)*Cd(9948)*Cd(11897) + 
     &   2*Cd(11900)

        Cd(11906) = AtC*Cd(11895) + Cd(251)*Cd(9948)*Cd(11902) + 
     &   2*Cd(11905)

        Cd(11915) = (-(MGlmT2*Cd(11914)) + Cd(12549)*MStop2(2))/
     &   MStop2(2)

        Cd(11921) = (2*MGl2*Cd(179)*Cd(11801))/MStop2(2) + 
     &   (Cd(50)*(MUE*Cd(11915) + MUEC*Cd(11920)))/MStop2(3)

        Cd(11937) = MGlmT2*Cd(11936) + 
     &   2*Cd(13)*Cd(9948)*Cd(11926)*Cr(75) + 2*Cd(11930)*MStop2(2)

        Cd(11938) = -(Cd(50)*Cd(11937)) + 
     &   MGl2*Cd(179)*(Cd(13)*Cd(179) + 2*Cd(10831))*MStop2(3)

        Cd(11945) = 8*MUE2*Cd(11941) - 
     &   8*(MUE**2*Cd(11943)*UCStop(1,3)**2 + 
     &      MUEC**2*Cd(11944)*UCStopC(1,3)**2)

        Cd(11956) = MGl2*Cd(251)*Cd(10751) + MGlmT2*Cd(11955) + 
     &   Cd(11954)*MStop2(2)

        Cd(11957) = -(Cd(9948)*Cd(11951)) + Cd(50)*Cd(11956) + 
     &   MGl2*Cd(179)*Cd(11947)*MStop2(3)

        Cd(11978) = -(Cd(47)*Cd(10658)) + Cd(9948)*Cd(11977) + 
     &   Cd(50)*Cr(28)*UStop2(1,1)

        Cd(11981) = AtC**2*Cd(11979)*UCStop(1,3)**2 + 
     &   At*(AtC*Cd(11978) + At*Cd(11980)*UCStopC(1,3)**2)

        Cd(11982) = (-(Cd(50)*Cd(11732)*MGlpTmSt2(1)*
     &        MStop2(3)) + 
     &     Cd(47)*(2*Cd(11981) + Cd(50)*Cd(11976)*MStop2(3)))/
     &   (Cd(50)*MStop2(3))

        Cd(11986) = -(AtC**2*Cd(11984)*UCStop(1,3)**2) + 
     &   At*(AtC*Cd(11983) + At*Cd(11985)*UCStopC(1,3)**2)

        Cd(11988) = Cd(47)**2*Cd(11987) + 
     &   Cd(50)**2*(Cd(47) - 2*At*AtC*MGlpTmSt2(1)*UStop2(1,1)**2)

        Cd(11999) = M_3*Cd(50)*Cd(11997) + 
     &   Cd(47)*(M_3*Cd(47)*Cd(179)*MGlpTmSt2(2) + 
     &      At*Cd(50)*Cd(11998)*UStop2(1,1))

        Cd(12001) = Cd(47)*UStop2(1,1)*
     &    (UStop2(1,1) - UStop2(1,2)) + Cd(12000)*UStop2(1,2)**2

        Cd(12012) = Cd(47)**2*Cd(12011) + 
     &   Cd(50)**2*(Cd(47)*Cd(12009) - 
     &      2*At*AtC*MGlpTmSt2(1)*UStop2(1,1)**2)

        Cd(12040) = AtC*Cd(12033) - 
     &   At*(Cd(12039) + At*Cd(11985)*UCStopC(1,3)**2)

        Cd(12043) = (2*Cd(47)*Cd(50)*Cd(9955)*Cd(12042) + 
     &     Cd(12045)*MStop2(3))/MStop2(3)

        Cd(12048) = Cd(50)*Cd(12047) + 
     &   AtC*Cd(12042)*UCStop(1,3) + At*Cd(47)*Cr(75)*UCStopC(1,3)

        Cd(12053) = AtC**2*Cd(12050)*UCStop(1,3) + 
     &   At*(2*AtC*Cd(12051)*Cr(5) + At*Cd(12052)*UCStopC(1,3))

        Cd(12067) = Cd(47)**2*Cd(12066) - 
     &   Cd(50)**2*(Cd(47)*Cd(12058) + Cd(12060)*MGlpTmSt2(1))

        Cd(12121) = (-(MUE**2*Cd(12120)*UCStop(1,3)**2) - 
     &     MUEC**2*Cd(47)**2*UCStopC(1,3)**2 + 
     &     2*MUE2*Cd(47)**2*Cd(12108)*UStop2(1,1))/Cd(47)**2

	Cd(12131) = Cd(12127) - MGl2*Cd(12130)*MStop2(3)

        Cd(12138) = -2*MTy*Cd(12135) + Cd(6680)*Cd(12137) - 
     &   AtC*Cd(12136)*UCStop(1,3)

        Cd(1208) = 2*CA2*Cd(1206) - S2A*Cd(1207) - 
     &   4*SA2*Cd(1201)*UStop2(1,2)

        Cd(12201) = AtC*MTy*Cd(47)*Cd(12194) + MUE*Cd(12200) - 
     &   MTy2*MUEC*Cd(47)*UUStop(1,2)

        Cd(12233) = MTy*(2*Cd(12227) + MUE*Cd(12232)) - 
     &   MUEC*Cd(12230)*UUStop(1,2)

        Cd(12276) = MGlpT2*Cd(47)*Cd(424)*Cd(12272) + 
     &   MUE*UCStop(1,3)*(-(MTy*Cd(12275)) + 
     &      2*AtC*MGlpT2*Cd(47)*UCStop(1,3))

        Cd(12283) = -(MUEC*UCStopC(1,3)*
     &      (MTy*Cd(12282) + 
     &        2*At*MGlpT2*Cd(47)*Cr(85)*UCStopC(1,3))) + 
     &   Cd(12276)*UStop2(1,1)

        Cd(12290) = Cd(50)*Cd(12070)*Cd(12289)*MStop2(3) + 
     &   Cd(47)**2*(At*Cd(12288) + AtC*Cd(12285)*UStop2(1,1))

        Cd(12334) = (MTy*(Cd(12328) + Cd(12333)) + 
     &     2*Cd(13)*Cd(12317)*MStop2(1))/
     &   (MTy*Cd(47)**2*MStop2(1)**2)

        Cd(12354) = -4*Cd(50)**2*Cd(12349) + 
     &   Cd(251)*Cd(2681)*Cd(12353) + 
     &   4*Cd(4187)*Cd(12350)*MStop2(2)

        Cd(12409) = MGlmT2*Cd(11768)*MStop2(1) + 
     &   M_3C*Cd(11763)*MStop2(1)**2 + Cd(12408)*UCStopC(1,3)

        Cd(12461) = 2*Cd(13)*Cd(47)**2*Cd(10407) - 
     &   2*Cd(47)*Cd(12460)*MStop2(1) + Cd(12446)*MStop2(3)

        Cd(12501) = 9*Cd(50)**2*Cd(491)*UStop2(1,1) + 
     &   UStop2(1,2)*(2*M_3*AtC*Cd(12500)*UStop2(1,1) + 
     &      9*Cd(50)**2*UStop2(1,2))

        Cd(12513) = AtC*(2*Cd(12511) + 
     &      AtC*Cd(10086)*Cd(12483)*UCStop(1,3)) + 
     &   At*Cd(12512)*UCStopC(1,3)

	Cd(12514) = -2*MGlmT2*Cd(50)*Cd(10029) + 2*Cd(12482) - Cd(12483)

        Cd(12531) = 2*MGlmT2*Cd(50)*Cd(12530) + 
     &   UStop2(1,1)*(8*M_3*Cd(12482) + At*Cd(12483)*UStop2(1,1))

        Cd(12532) = 4*M_3C*MGlmT2*Cd(50) + 
     &   8*M_3C*Cd(12482) + AtC*Cd(12483)

	Cd(12538) = Cd(12487) + Cd(12488)*MStop2(2)

        Cd(12544) = 4*MGlmT2*Cd(50)*Cd(9997) + 
     &   8*M_3*Cd(12482) + At*Cd(12483)

        Cd(12546) = 2*MGlmT2*Cd(50)*Cd(9998) + 
     &   UStop2(1,1)*(8*M_3C*Cd(12482) + 
     &      AtC*Cd(12483)*UStop2(1,1))

        Cd(12550) = (Cd(12434)*Cd(12483) + 
     &     (2*Cd(50)*MStop2(2)*
     &        (-(MGlmT2*Cd(11914)) + Cd(12549)*MStop2(2)))/
     &      MStop2(3))/MStop2(2)**2

	Cd(12559) = 2*MGlmT2*Cd(50)*Cd(10025) + 2*Cd(12482) - Cd(12483)

        Cd(12577) = 2*Cd(9948)*Cd(12576) - 
     &   4*Cd(50)*Cd(11956)*MStop2(2) - 
     &   Cd(12483)*Cd(12575)*MStop2(3)

        Cd(12582) = 2*MTy2*Cd(585)*Cd(12581) - 
     &   Cd(50)**2*Cd(152)*Cd(4095)*Cr(5) + 
     &   Cd(47)**2*Cd(179)*Cd(4190)*Cr(5)

        Cd(12586) = 8*MTy2**2*
     &    (Cd(50)**2*Cd(164) + Cd(47)**2*Cd(251)) + 
     &   Cd(11962)*Cd(12585)

        Cd(12592) = 2*Cd(50)**2*Cd(164) + 2*Cd(47)**2*Cd(251) + 
     &   Cd(10193)*Cd(11962)

        Cd(12593) = 4*M_3*Cd(50)**2*Cd(164) + 
     &   4*M_3*Cd(47)**2*Cd(251) - At*Cd(11962)*UStop2(1,1)

        Cd(12598) = Cd(50)**2*Cd(12595) + Cd(47)**2*Cd(12597) + 
     &   Cd(9948)*Cd(12593)*UStop2(1,1)

        Cd(12604) = Cd(50)**2*Cd(12601) + Cd(47)**2*Cd(12603) + 
     &   8*M_3C*Cd(9948)*Cd(12599)*UStop2(1,1)

        Cd(12605) = -(At*Cd(12604)) + 
     &   AtC*(-2*Cd(12598) + AtC*Cd(11962)*UCStop(1,3)**2)

        Cd(12622) = -2*MTy*Cd(47)**2*Cd(50)**2*
     &    (M_3*Cd(12610) + M_3C*Cd(12613)) + MGlpT2*Cd(12621)

        Cd(12636) = MUE*Cd(12635)*UStop2(1,1)*UStop2(1,2) - 
     &   MTy*MUEC*Cd(9794)*MGlpTmSt2(1)*UCStopC(1,3)**2*
     &    UStopC(1,1) - Cd(12634)*UUStopC(1,1)

        Cd(12645) = (Cd(50)**2*Cd(12644) + 
     &     MUE*Cd(47)**2*Cd(12639)*UCStop(1,3)**2 + 
     &     MUEC*Cd(47)**2*Cd(50)**2*Cd(12643)*UStop2(1,1))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12653) = 2*Cd(47)**2*Cd(12597) + 
     &   Cd(50)**2*Cd(12652) + 2*Cd(9948)*Cd(12593)*UStop2(1,1)

        Cd(12658) = Cd(50)**2*Cd(12655) + Cd(47)**2*Cd(12657) + 
     &   8*M_3C*Cd(9948)*Cd(12599)*UStop2(1,1)

        Cd(12659) = At*Cd(12658) + 
     &   AtC*(Cd(12653) - AtC*Cd(11962)*UCStop(1,3)**2)

        Cd(12670) = (-(Cd(50)*Cd(11732)*MGlpTmSt2(1)*
     &        MStop2(3)) + 
     &     Cd(47)*(2*Cd(11981) + Cd(50)*Cd(12669)*MStop2(3)))/
     &   (Cd(50)*MStop2(3))

        Cd(12684) = -(MUEC*UCStopC(1,3)*
     &      (MTy*Cd(12683) + 
     &        2*At*MGlpT2*Cd(47)*Cr(85)*UCStopC(1,3))) + 
     &   Cd(12276)*UStop2(1,1)

        Cd(12777) = Cd(12701)*Cd(12774) + 
     &   8*MStop2(2)*(Cd(12776) - 
     &      Cd(12701)**2*Cd(12763)*UCStopC(1,3)**2*UStop2(1,2))

        Cd(12818) = Cd(12817)*UStop2(1,2) + 
     &   Cd(12701)*(Cd(12815) + 
     &      2*Cd(12701)*MStop2(1)*UCStopC(1,3)**2*UStop2(1,1)*
     &       UStop2(1,2))

	Cd(13027) = -2*M_3*Cd(13006) + 2*M_3C*Cd(13011) + Cd(13026)

        Cd(13042) = 1/(64.D0*Pi**4)*
     &   ((Cd(13027) - 64*Pi**4*
     &          (Cd(13041)*MGlpTmStxGlT4(2)*MStop2(2)))*
     &       MGlpTmStxGlT4(2) + MTy2*Cd(251)*Cd(12999)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(2)**2*MStop2(2)*MStop2(3))

        Cd(13073) = 1/Pi**4*
     &    (Cd(12836)**2*Cd(13044) - 
     &       64*Pi**4*(TB2*Cd(13057)*MGlpTmStxGlT4(1)))/
     &     (TB2*MGlpTmStxGlT4(1)**2) + 
     &   Cd(13072)/MGlpTmStxGlT4(2)**2

        Cd(13086) = 1/(64.D0*Pi**4)*
     &    (-(64*Pi**4*(Cd(13085)*MGlpTmStxGlT4(1)**2*MStop2(1))) + 
     &       MGl2*Cd(152)*Cd(12836)**2*Cr(45)*MStop2(3))/
     &     (TB2*MGlpTmStxGlT4(1)**2*MStop2(1)*MStop2(3)) + 
     &   Cd(13084)/MGlpTmStxGlT4(1)

        Cd(13098) = -(Cd(12700)*Cd(13097)) + 
     &   Cd(12701)*(Cd(13096) - 
     &      Cd(251)*Cd(12701)*MStop2(3)*UCStopC(1,3)**2)

        Cd(13224) = (Cd(20)*Cd(13200) + 4*Cd(13162)*MStop2(1) + 
     &     2*Cd(17)*Cd(13223)*MStop2(1)**2*MStop2(2))/MStop2(1)

        Cd(13293) = Cd(12701)*Cd(13024) - 
     &   Cd(179)*Cd(12700)*Cd(13025) - Cd(13292)*MStop2(3)

	Cd(1323) = S2A*Cd(1276) + SA2*Cd(1294) + 2*CA2*Cd(1322)

	Cd(13294) = 2*M_3C*Cd(13011) - 2*M_3*Cd(13291) + Cd(13293)

        Cd(13304) = Cd(12701)*Cd(13302) - 
     &   1/(256.D0*Pi**4)*(Cd(13303) - 
     &       2*Cd(13028)*MStop2(2)*MStop2(3))/MStop2(2)**2

        Cd(13305) = 1/(256.D0*Pi**4)*
     &   (-2*(Cd(13294) - 128*Pi**4*
     &          (Cd(13304)*MGlpTmStxGlT4(2)*MStop2(2)))*
     &       MGlpTmStxGlT4(2)*MStop2(2) + Cd(13290)*MStop2(3))/
     &    (TB2*MGlpTmStxGlT4(2)**2*MStop2(2)**2*MStop2(3))

        Cd(13315) = 1/2.D0*
     &   (4*Cd(13314) + A0delStop(1)*Cd(12992)*MStop2(1) + 
     &      2*Cd(13287)*Cr(13)*MStop2(1)**2)/MStop2(1)

        Cd(13319) = 1/2.D0*
     &   (4*Cd(13318) + A0delStop(2)*Cd(13042)*MStop2(2) + 
     &      2*Cd(13305)*Cr(34)*MStop2(2)**2)/MStop2(2)

        Cd(13356) = Cd(12701)*Cd(13302) - 
     &   1/(256.D0*Pi**4)*(MTy2*Cd(13303) - 
     &       2*Cd(13102)*MStop2(2)*MStop2(3))/(MTy2*MStop2(2)**2)

        Cd(13363) = 1/2.D0*
     &   (4*Cd(13362) + A0delStop(1)*Cd(13086)*MStop2(1) + 
     &      2*Cd(13348)*Cr(13)*MStop2(1)**2)/MStop2(1)

        Cd(13460) = 1/32.D0*
     &   (1/Pi**4*(MGl2*MTy2*Cd(12999))/MGlpTmStxGlT4(2)**2 + 
     &      (-(1/Pi**4*Cd(13458)) + 32*Cd(13459)*MStop2(3))/
     &       MStop2(3)**2)/TB2

        Cd(13465) = 1/(64.D0*Pi**4)*
     &   (TB2*(64*Pi**4*(Cd(64)*Cd(13464)*MGlpTmStxGlT4(2)**2) + 
     &         64*Pi**4*(Cd(13460)*MGlpTmStxGlT4(2)**2*
     &            T134(MGl2,MTy2,MStop2(2)))) - 
     &      Cd(1508)*Cd(12999)*MGlpTmSt4(2))/
     &    (TB2*MGlpTmStxGlT4(2)**2)

        Cd(13868) = Cd(13801)*Cd(13867) + 
     &   (-(M_3C*Cd(13775)) + M_3*Cd(13779))*MSq2Diff(1,1)

        Cd(13895) = Cd(47)*Cd(138) - Cd(13894)*UStop2(1,1) + 
     &   Cd(50)*Cd(118)*UStop2(1,2)

        Cd(1386) = Cd(1383)*MStop2(3) - 
     &   AtC*MTy*Cd(1382)*UCStop(1,3) - Cd(1385)*UStop2(1,1)

        Cd(13923) = Cd(47)*Cd(138) - Cd(13894)*UStop2(1,1) + 
     &   Cd(50)*Cd(118)*UStop2(1,2)

        Cd(14112) = (Cd(13469)*Cd(13821))/MSq2Diff(1,1)**2 + 
     &   (MTy*Cd(14106)*MSq2Diff(1,1)**2 + 
     &      MTy2*Cd(14111)*MSq2Diff(1,2))/
     &    (MTy*MTy2*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2)

        Cd(14140) = Cd(9948)*Cd(14129) - 
     &   Cd(14139)*MGlpTmSt2(1) + 
     &   TB*Cd(164)*Cd(14130)*MSq2Diff(1,1)

        Cd(14162) = Cd(9948)*Cd(14157) + Cd(47)*Cd(14161) + 
     &   Cd(14159)*MSq2Diff(1,1)

        Cd(14207) = 2*M_3C*At*Cd(50) + AtC*Cd(14205) - 
     &   Cd(14206)*MSq2Diff(1,2)

	Cd(14220) = Cd(50)*Cd(9948) - Cd(14201)

        Cd(14245) = Cd(50)*Cd(14244) - 
     &   MSq2Diff(1,1)*(Cd(14241)*MGlpTmSt2(2)*MSq2Diff(1,2) + 
     &      AtC*Cd(14234)*UCStop(1,3))

        Cd(1420) = 2*Cd(1417)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(1419)*UStop2(1,1) + 
     &      4*MStop2(2)*MStop2(3)*UStop2(1,2))

        Cd(14258) = Cd(50)*Cd(14250) + AtC*Cd(14257) + 
     &   Cd(14254)*MSq2Diff(1,1)*MSq2Diff(1,2)

        Cd(14294) = -(Cd(14293)*UCStopC(1,3)) + 
     &   MTy*Cd(14291)*UStop2(1,2) + MGl2*Cd(14289)*UStopC(2,1)

        Cd(14295) = MTy*(Cd(12700)*Cd(14286) + Cd(14294)) + 
     &   Cd(14283)*MSq2Diff(1,2)**2

        Cd(14351) = MTy*(Cd(14102) + Cd(14110)*MSq2Diff(1,1)) + 
     &   MTy2*Cd(13469)*Cd(13821)*MSq2Diff(1,2)

	Cd(1449) = 2*Cd(1293) - 2*Cd(1363) - Cd(1446)

        Cd(1450) = 2*MUE2*Cd(1447) - 
     &   MUE**2*Cd(1448)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(1449)*UCStopC(1,3)**2

        Cd(1479) = Cd(1313)*MStop2(3) - Cd(1476)*UStop2(1,1) - 
     &   Cd(1478)*UStop2(1,2)

        Cd(1521) = Cd(64)*Cd(1515) - 
     &   2*Cd(1520)*T134(MGl2,MTy2,MStop2(2)) - 
     &   2*Cd(1508)*Cd(1510)*MStop2(3)**2*UStop2(1,2)

        Cd(1579) = Cd(1566) - 
     &   2*M_3*MTy*Cd(1578)*UCStopC(1,3) - Cd(1575)*UStop2(1,1)

        Cd(1628) = Cd(1618) + 
     &   2*M_3*MTy*Cd(1627)*UCStopC(1,3) - Cd(1625)*UStop2(1,1)

        Cd(1655) = Cd(52)*Cd(1643) - 
     &   Cd(1654)*T134(MGl2,MTy2,MStop2(1)) - 
     &   2*Cd(1493)*Cd(1633)*MStop2(3)**2*UStop(1,1)*UStop2(1,2)

        Cd(1676) = Cd(64)*Cd(1665) - 
     &   Cd(1675)*T134(MGl2,MTy2,MStop2(2)) + 
     &   2*Cd(1508)*Cd(1656)*MStop2(3)**2*UStop(1,1)*UStop2(1,2)

        Cd(1861) = Cd(64)*Cd(1844) - 
     &   2*Cd(96)**2*Cd(1508)*Cd(1843)*MStop2(3)**2 - 
     &   Cd(1860)*T134(MGl2,MTy2,MStop2(2))

        Cd(1881) = -(Cd(1625)*MStop2(3)*UCStop(1,3)) + 
     &   MTy*(4*M_3*Cd(47)*MGlpTmSt2(1)*MStop2(3) + 
     &      2*M_3*Cd(1880)*UStop2(1,2))

        Cd(1882) = Cd(52)*Cd(1881) + 
     &   2*Cd(86)**2*Cd(1493)*Cd(1862)*MStop2(3)**2 + 
     &   Cd(1878)*T134(MGl2,MTy2,MStop2(1))

        Cd(1960) = Cd(44)**3*Cd(1948) - MGl2*Cd(1959) - 
     &   2*MTy2*Cd(47)**2*UStop2(1,2)

        Cd(1969) = -(Cd(1968)*MStop2(3)) + 
     &   UStop2(1,1)**2*(Cd(1960) - 2*MTy2*Cd(1946)*UStop2(1,1))

        Cd(2035) = Cd(48)**3*Cd(1948) - MGl2*Cd(2034) - 
     &   2*MTy2*Cd(50)*Cd(2028)*UStop2(1,2)

        Cd(2045) = Cd(2044)*MStop2(3) + 
     &   UStop2(1,1)**2*(Cd(2035) - 2*MTy2*Cd(2022)*UStop2(1,1))

        Cd(2046) = Cd(64)*Cd(2016) - 
     &   Cd(2045)*T134(MGl2,MTy2,MStop2(2)) - 
     &   4*MTy2*Cd(1508)*Cd(2018)*MStop2(3)**2*UStop2(1,1)

        Cd(2061) = 2*Cd(44)**3*Cd(2048) - 2*MGl2*Cd(2060) - 
     &   Cd(2056)*UStop2(1,2)

        Cd(2145) = (Cd(44)**3*Cd(2134) - MGl2*Cd(2142))*
     &    UStop(1,1) + 2*M_3*MTy*Cd(2144)*UStop(1,2)

        Cd(2219) = -(M_3*Cd(2218)*MStop2(3)*UStop(1,2)) + 
     &   MTy*(Cd(2216)*UStop(1,1) + 
     &      2*Cd(2213)*MStop2(3)*UStop2(1,2))

        Cd(237) = 1/2.D0*Cd(236) - Cd(222)*Cr(13)*MStop2(1) + 
     &   Cd(232)*Cr(34)*MStop2(2)

        Cd(2411) = M_3*Cd(2402) + 
     &   MTy*(Cd(2406)*UCStop(1,3) - 2*Cd(2410)*UUStopC(1,2))

        Cd(2428) = M_3C*MStop2(1)**3*MStop2(3)*
     &    UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(2411) - Cd(2401)*UStop2(1,1))

        Cd(2493) = -(M_3C*MStop2(2)**3*MStop2(3)*
     &      UCStop(1,3)**2) + 
     &   UStop2(1,1)*(Cd(2483) - Cd(2466)*UStop2(1,1))

        Cd(249) = (4*Cd(243)*MStop2(3) - 
     &     AtC*Cd(248)*UCStop(1,3) + 
     &     2*At*Cd(238)*MStop2(3)*UCStopC(1,3))/MStop2(3)

        Cd(2609) = -(M_3*MTy*Cd(2608)*UCStopC(1,3)) + 
     &   Cd(2606)*UStop2(1,1) - Cd(2603)*UStop2(1,2)

        Cd(2616) = Cd(2606)*UStop2(1,1)**2 + 
     &   Cd(2611)*UStop2(1,2) + Cd(2615)*UStopC(1,1)

        Cd(2632) = -(M_3*MTy*Cd(2631)*UCStopC(1,3)) + 
     &   Cd(2619)*UStop2(1,1) + Cd(2629)*UStop2(1,2)

        Cd(2641) = 2*Cd(86)*Cd(2190)*MStop2(1)**2 + 
     &   MTy*(Cd(2623)*UStop2(1,1) + Cd(2640)*UStop2(1,2))

        Cd(2649) = MTy*(2*M_3C*Cd(2645)*UCStop(1,3)*
     &       UStop2(1,1) + 
     &      8*M_3*MGlpTmSt2(1)*MStop2(1)**2*UCStopC(1,3)*
     &       UStop2(1,2)) - Cd(2648)*UStopC(1,1)

        Cd(2670) = Cd(2669)*UStop2(1,1) + 
     &   UCStop(1,3)*(M_3C*MTy2*Cd(2659)*UCStop(1,3) - 
     &      MTy*Cd(2639)*UStop2(1,2))

        Cd(2674) = 2*Cd(86)*Cd(2190)*MStop2(1)**2 + 
     &   MTy*(Cd(2623)*UStop2(1,1) + Cd(2673)*UStop2(1,2))

        Cd(2692) = M_3*MTy*Cd(2691)*UCStopC(1,3) + 
     &   Cd(2689)*UStop2(1,1) - Cd(2685)*UStop2(1,2)

        Cd(2697) = Cd(2689)*UStop2(1,1)**2 + 
     &   Cd(2694)*UStop2(1,2) + Cd(2696)*UStopC(1,1)

        Cd(2708) = M_3*MTy*Cd(2707)*UCStopC(1,3) + 
     &   Cd(2699)*UStop2(1,1) + Cd(2705)*UStop2(1,2)

        Cd(2717) = 2*Cd(96)*Cd(2244)*MStop2(2)**2 + 
     &   MTy*(Cd(2701)*UStop2(1,1) + Cd(2716)*UStop2(1,2))

        Cd(2727) = MTy*(2*M_3C*Cd(2721)*UCStop(1,3)*
     &       UStop2(1,1) + 
     &      8*M_3*MGlpTmSt2(2)*MStop2(2)**2*UCStopC(1,3)*
     &       UStop2(1,2)) + Cd(2726)*UStopC(1,1)

        Cd(2746) = Cd(2745)*UStop2(1,1) + 
     &   UCStop(1,3)*(M_3C*MTy2*Cd(2735)*UCStop(1,3) + 
     &      MTy*Cd(2715)*UStop2(1,2))

        Cd(2750) = 2*Cd(96)*Cd(2244)*MStop2(2)**2 + 
     &   MTy*(Cd(2701)*UStop2(1,1) + Cd(2749)*UStop2(1,2))

        Cd(2815) = 3*Cd(2766) + 
     &   2*M_3C*MTy*Cd(2761)*UCStop(1,3) - 
     &   Cd(2814)*UStop2(1,2)**2

        Cd(2940) = Cd(2939)*MStop2(1) + 
     &   MTy*(-(Cd(2920)*UStop2(1,1)) + Cd(2908)*UStop2(1,2))

        Cd(3004) = MTy*(M_3C*Cd(2963)*UCStop(1,3)*
     &       UStop2(1,1) - 
     &      M_3*Cd(2972)*MStop2(1)*UCStopC(1,3)*UStop2(1,2)) + 
     &   Cd(3003)*UStopC(1,1)

        Cd(3092) = 2*MTy*Cd(3076) + 
     &   M_3C*Cd(3065)*UCStop(1,3)**2 - Cd(3091)*UStop2(1,1)

        Cd(3106) = Cd(2939)*MStop2(1) + 
     &   MTy*(-(Cd(2920)*UStop2(1,1)) + Cd(3105)*UStop2(1,2))

        Cd(3178) = 3*Cd(3125) - 
     &   2*M_3C*MTy*Cd(3177)*UCStop(1,3) - 
     &   Cd(3176)*UStop2(1,2)**2

        Cd(3305) = Cd(3304)*MStop2(2) + 
     &   MTy*(Cd(3284)*UStop2(1,1) - Cd(3271)*UStop2(1,2))

        Cd(3366) = MTy*(M_3C*Cd(3325)*UCStop(1,3)*
     &       UStop2(1,1) + 
     &      M_3*Cd(3330)*MStop2(2)*UCStopC(1,3)*UStop2(1,2)) - 
     &   Cd(3365)*UStopC(1,1)

        Cd(335) = Cd(330)*UCStop(1,3) + 
     &   UStop2(1,1)*(Cd(334) - 2*M_3*MTy*Cd(50)*UStop2(1,1))

        Cd(3419) = Cd(3366) + Cd(3418) + 
     &   MTy2*Cd(3315)*UStop2(1,1)**2*UStopC(1,1)

        Cd(3466) = 2*MTy*Cd(3450) + 
     &   M_3C*Cd(3433)*UCStop(1,3)**2 - Cd(3465)*UStop2(1,1)

        Cd(342) = 2*Cd(243) + AtC*Cd(341)*UCStop(1,3) + 
     &   At*Cd(238)*UCStopC(1,3)

        Cd(3483) = Cd(3304)*MStop2(2) + 
     &   MTy*(Cd(3475)*UStop2(1,1) - Cd(3482)*UStop2(1,2))

        Cd(352) = 2*Cd(351) + AtC*Cd(255)*UCStop(1,3) + 
     &   At*Cd(345)*UCStopC(1,3)

        Cd(3577) = MGl2*(Cd(3574) - MGl2*Cd(3576)) - 
     &   Cd(44)**2*Cd(2048)*MStop2(1)

        Cd(3606) = Cd(3605) + 
     &   MTy*(-2*Cd(3590)*UStop2(1,1) + 
     &      8*Cd(3104)*MStop2(3)*UStop2(1,2))

        Cd(3618) = MGl2*(Cd(3561) - MGl2*Cd(3617)) - 
     &   2*Cd(44)**2*Cd(2048)*MStop2(1)

        Cd(3657) = -(Cd(44)**2*Cd(2115)*MStop2(1)) + 
     &   MGl2*(2*Cd(3656) - 4*MGl2**3*MStop2(2))

        Cd(3729) = 2*Cd(3728)*MStop2(3) - 
     &   2*M_3C*MTy*Cd(3717)*UCStop(1,3) + Cd(3727)*UUStopC(1,2)

        Cd(3881) = MGl2*(Cd(3824) - MGl2*Cd(3880)) - 
     &   2*Cd(48)**2*Cd(2048)*MStop2(2)

        Cd(3890) = MGl2*(Cd(3857) - 2*MGl2*Cd(3889)) - 
     &   Cd(48)**2*Cd(2115)*MStop2(2)

        Cd(3981) = 2*Cd(3980)*MStop2(3) + 
     &   2*M_3C*MTy*Cd(3407)*UCStop(1,3) + Cd(3979)*UUStopC(1,2)

        Cd(390) = 1/2.D0*Cd(389) - Cd(384)*Cr(13)*MStop2(1) + 
     &   Cd(385)*Cr(34)*MStop2(2)

        Cd(4004) = 2*M_3C*MTy*Cd(3982)*UCStop(1,3)*
     &    UStop2(1,1) - 2*Cd(4003)*UStopC(1,1) + 
     &   Cd(4002)*UUStopC(1,3)

        Cd(4026) = -(M_3*MTy*Cd(4025)*UCStopC(1,3)) + 
     &   Cd(4023)*UStop2(1,1) + Cd(4020)*UStop2(1,2)

        Cd(4037) = Cd(4023)*UStop2(1,1)**2 - 
     &   Cd(4031)*UStop2(1,2) + Cd(4036)*UStopC(1,1)

        Cd(4044) = -(M_3*MTy*Cd(4043)*UCStopC(1,3)) + 
     &   Cd(4023)*UStop2(1,1) + Cd(4020)*UStop2(1,2)

        Cd(4126) = M_3*MTy*Cd(4125)*UCStopC(1,3) + 
     &   Cd(4123)*UStop2(1,1) + Cd(4120)*UStop2(1,2)

        Cd(4132) = Cd(4123)*UStop2(1,1)**2 - 
     &   Cd(4128)*UStop2(1,2) + Cd(4131)*UStopC(1,1)

	Cd(4199) = Cd(4183) - Cd(4197) - MTy*Cd(4198)

        Cd(4277) = Cd(4275) + 
     &   2*Cd(4249)*UStop2(1,1)**2*UStopC(1,1) + 
     &   2*Cd(47)*Cd(4276)*MStop2(1)*UStop2(1,2)*UStopC(1,2)

        Cd(4486) = MTy2*(M_3C*Cd(4452)*UCStop(1,3)*
     &       UStop2(1,1) - 
     &      M_3*Cd(4457)*MStop2(1)*UCStopC(1,3)*UStop2(1,2)) + 
     &   Cd(4485)*UStopC(1,1)

        Cd(4538) = Cd(4486) + Cd(4537) + 
     &   MTy*MTy2*Cd(4444)*UStop2(1,1)**2*UStopC(1,1)

        Cd(4623) = Cd(4621) + 
     &   2*Cd(4596)*UStop2(1,1)**2*UStopC(1,1) - 
     &   2*Cd(50)*Cd(4622)*MStop2(2)*UStop2(1,2)*UStopC(1,2)

        Cd(4819) = MTy2*(M_3C*Cd(4787)*UCStop(1,3)*
     &       UStop2(1,1) - 
     &      M_3*Cd(4792)*MStop2(2)*UCStopC(1,3)*UStop2(1,2)) - 
     &   Cd(4818)*UStopC(1,1)

        Cd(4876) = Cd(4819) + Cd(4875) + 
     &   MTy*MTy2*Cd(4780)*UStop2(1,1)**2*UStopC(1,1)

        Cd(4901) = Cd(4621) + 
     &   2*Cd(4899)*UStop2(1,1)**2*UStopC(1,1) - 
     &   2*Cd(50)*Cd(4622)*MStop2(2)*UStop2(1,2)*UStopC(1,2)

        Cd(5004) = Cd(44)**3*Cd(48)**2*Cd(4954) + 
     &   MGl2*(-2*Cd(44)*Cd(4968) + MGl2*Cd(5003))

        Cd(5265) = 4*Cd(47)*Cd(50)*Cd(5229) + 
     &   4*M_3*MTy*Cd(5264)*UCStopC(1,3) - 2*Cd(5259)*UStop2(1,1)

        Cd(5281) = MTy*(2*M_3C*Cd(5280)*UUStop(1,1) - 
     &      8*M_3*Cd(5130)*UUStop(1,2)) - 2*Cd(5259)*UUStop(1,3)

        Cd(5282) = Cd(5265)*UCStop(1,3) - 
     &   8*M_3*MTy*Cd(47)*Cd(50)*Cd(295)*UStop2(1,1)**2 + 
     &   Cd(5281)*UUStopC(1,2)

	Cd(5381) = Cd(5360) - Cd(5365) - MGl2*Cd(5380)

        Cd(534) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(S2A*Cd(517) + 
     &        MStop2(3)*(CA2*Cd(533)*MStop2(3) + 
     &           2*SA2*Cd(518)*UStop2(1,1)*UStop2(1,2))))/
     &    (SB2*SW2*MStop2(3)**2)

        Cd(5736) = Cd(44)**4*Cd(48)**3 + 
     &   MGl2*(Cd(44)**2*Cd(48)*Cd(5725) - MGl2*Cd(5735))

        Cd(5803) = 2*Cd(44)*Cd(48)*Cd(5785) - MGl2*Cd(5802) - 
     &   2*MTy2*Cd(47)**2*Cd(50)*UStop2(1,2)

        Cd(5840) = Cd(5839)*MStop2(3) + 
     &   2*UStop2(1,1)**2*(Cd(5803) - 2*MTy2*Cd(5776)*UStop2(1,1))

        Cd(5978) = 2*Cd(44)*Cd(48)*
     &    (Cd(44)*Cd(48) + MGl2*Cd(89)) - MGl2*Cd(5977)

        Cd(6202) = Cd(44)**4*Cd(48)**2 + 
     &   MGl2*(2*Cd(44)**2*Cd(6175) - MGl2*Cd(6201))

        Cd(6347) = UStop2(1,1)*
     &    (Cd(6346) - Cd(47)*Cd(6337)*UStop2(1,1)) + 
     &   4*Cd(47)*Cd(50)*Cd(295)*UStop2(1,2)

	Cd(631) = S2A*Cd(612) + 4*SA2*Cd(620) + 2*CA2*Cd(630)

        Cd(6432) = 2*Cd(47)*Cd(6431)*MStop2(3) + 
     &   UStop2(1,1)*(2*Cd(6430) - Cd(6390)*UStop2(1,1))

        Cd(6673) = Cd(6996) + 
     &   2*M_3C*MTy*Cd(6660)*UCStop(1,3)*UStop(1,1)

        Cd(6709) = Cd(6708)*MStop2(3)*UCStop(1,3) + 
     &   MTy*(4*M_3*Cd(47)*Cd(50)**2 - 
     &      2*M_3*Cd(6704)*UStop2(1,2))

        Cd(7014) = Cd(7013)*UStop2(1,1)*UStopC(1,1) - 
     &   2*Cd(47)*Cd(50)*Cd(6675)*UStop2(1,2)*UStopC(1,2) - 
     &   Cd(7012)*UUStopC(1,3)

        Cd(7042) = Cd(44)**4*Cd(48)**2 + 
     &   MGl2*(-2*Cd(44)*Cd(7033) - MGl2*Cd(7041))

	Cd(7165) = -Cd(5079) + 2*Cd(7164)

        Cd(7204) = 2*Cd(47)*Cd(50)*Cd(5229) + 
     &   2*M_3*MTy*Cd(7203)*UCStopC(1,3) - Cd(7198)*UStop2(1,1)

        Cd(711) = CA2*Cd(707) - Cd(710)*MStop2(1) + 
     &   Cd(360)*Cd(704)*MStop2(2)

        Cd(7229) = MTy*(2*M_3C*Cd(7228)*UUStop(1,1) - 
     &      4*M_3*Cd(5130)*UUStop(1,2)) - Cd(7198)*UUStop(1,3)

        Cd(7230) = Cd(7204)*UCStop(1,3) - 
     &   4*M_3*MTy*Cd(47)*Cd(50)*Cd(295)*UStop2(1,1)**2 + 
     &   Cd(7229)*UUStopC(1,2)

        Cd(7306) = Cd(44)**2*Cd(48)**2*Cd(288) + 2*Cd(7300) + 
     &   MGl2*Cd(7305)

	Cd(7547) = Cd(8440) - Cd(5739)*UStop2(1,2)

        Cd(745) = Cd(741)*Conjugate(dY33eps(0)) + 
     &   Cd(740)*dY33eps(0) + Cd(744)*MStop2(3)

        Cd(7553) = Cd(44)**2*Cd(48)**2*Cd(5649) - 
     &   MGl2*Cd(7552) + 4*M_3C*MTy*MTy2*Cd(7526)*UCStop(1,3)

        Cd(747) = Cd(741)*Conjugate(dY33eps(1)) + 
     &   Cd(740)*dY33eps(1) + Cd(746)*MStop2(3)

        Cd(7601) = Cd(44)*Cd(48)*Cd(7585) - MGl2*Cd(7600) - 
     &   MTy2*Cd(47)**2*Cd(50)*UStop2(1,2)

        Cd(7743) = 2*Cd(7736) - Cd(7739)*MStop2(3)*UStop(1,1) - 
     &   M_3*MTy*Cd(7742)*UStop(1,2)

        Cd(773) = S2A*Cd(768) + 2*CA2*Cd(772) + 
     &   8*SA2*Cd(13)**2*Cr(75)*MStop2(4)

        Cd(7906) = Cd(44)*Cd(48)*Cd(7870) - MGl2*Cd(7905) - 
     &   MTy2*Cd(7899)*UStop2(1,2)

        Cd(7963) = 2*Cd(7736)*UStop(1,2) - 
     &   M_3*MTy*Cd(7962)*UUStop(1,2) - 
     &   Cd(7739)*MStop2(3)*UUStop(1,3)

	Cd(8038) = Cd(6367) - 2*Cd(8037)

        Cd(8101) = -(Cd(47)*Cd(8100)*MStop2(3)) + 
     &   UStop2(1,1)*(2*Cd(8099) - Cd(8060)*UStop2(1,1))

	Cd(8223) = -Cd(8261) + 2*MTy*Cd(8215)*UStop2(1,2)

        Cd(8224) = 2*Cd(8223) + 
     &   MTy*(Cd(8178)*UStop2(1,1) + 2*Cd(6586)*UStop2(1,2))

        Cd(8236) = Cd(8321) + 
     &   2*M_3C*MTy*Cd(8230)*MStop2(3)*UCStop(1,3)

        Cd(8251) = Cd(8250)*MStop2(3) + 
     &   MTy*(4*M_3*Cd(47)*Cd(50)**2 - 
     &      2*M_3*Cd(6704)*UStop2(1,2))

        Cd(8268) = Cd(44)**2*Cd(48)**4 + 
     &   MGl2*(-2*Cd(44)*Cd(48)**2*Cd(8262) - MGl2*Cd(8267))

        Cd(8276) = Cd(44)**2*Cd(48)**4 + 
     &   MGl2*(-(MGl2*Cd(8208)) - 2*Cd(44)*Cd(48)**2*Cd(8262))

        Cd(8287) = Cd(44)*Cd(48)*Cd(89)*Cd(5556) - 
     &   MGl2*Cd(8286) + Cd(8283)*UStop2(1,1)

        Cd(8337) = Cd(6996)*UStop(1,2) - 
     &   2*(M_3C*MTy*Cd(8336)*UUStop(1,1) + 
     &      Cd(8321)*MStop2(3)*UUStop(1,3))

        Cd(8339) = Cd(8338)*UStop2(1,1)*UStopC(1,1) - 
     &   Cd(47)*Cd(50)*Cd(6675)*UStop2(1,2)*UStopC(1,2) - 
     &   Cd(8337)*UUStopC(1,3)

        Cd(8396) = 2*M_3*MTy*Cd(8238)*MStop2(3)*UCStopC(1,3) + 
     &   Cd(8395)*UStop2(1,1) + Cd(8390)*UStop2(1,2)

	Cd(8420) = Cd(8440) - Cd(5463)*UStop2(1,2)

        Cd(8572) = Cd(44)**2*Cd(48)**2*Cd(89) + 
     &   MGl2*(2*Cd(44)*Cd(48)*Cd(5995) - MGl2*Cd(8571))

        Cd(8582) = 2*M_3*MTy2*Cd(8579) + 
     &   MTy*(-(Cd(8539)*MStop2(3)*UCStop(1,3)) + 
     &      2*Cd(8581)*UUStopC(1,2))

        Cd(8603) = Cd(8037) + 
     &   2*M_3*MTy*Cd(7044)*MStop2(3)*UCStopC(1,3) + 
     &   Cd(8602)*UStop2(1,1)

        Cd(8637) = 2*M_3*MTy2*Cd(8222)*UCStopC(1,3) + 
     &   MTy*(Cd(8635)*MStop2(3) + Cd(8636)*UStopC(1,1))

        Cd(8640) = Cd(8321) + 
     &   2*M_3C*MTy*Cd(8238)*MStop2(3)*UCStop(1,3)

        Cd(8707) = Cd(8703) + 
     &   4*M_3*MTy*Cd(8706)*UCStopC(1,3) - 2*Cd(4234)*UStop2(1,1)

        Cd(8754) = Cd(8745) - 
     &   4*MTy*(Cd(8753)*UStop2(1,1) + 
     &      2*Cd(8729)*MStop2(3)*UStop2(1,2))

        Cd(8765) = MTy*(2*M_3C*Cd(4400)*UUStop(1,1) + 
     &      M_3*Cd(8764)*UUStop(1,2)) + 
     &   2*Cd(4415)*MStop2(3)*UUStop(1,3)

        Cd(868) = -2*MTy*Cd(866) + At*Cd(867) + 
     &   8*AtC*Cd(795)*Cd(833)*Cr(75)

        Cd(8788) = Cd(8745) - 
     &   4*MTy*(Cd(8787)*UStop2(1,1) + 
     &      2*Cd(8786)*MStop2(3)*UStop2(1,2))

        Cd(8880) = 2*M_3C*MTy2*Cd(8864)*UCStop(1,3)*
     &    UStop2(1,1) - Cd(8879)*UStopC(1,1) + 
     &   MTy*Cd(8878)*UUStopC(1,3)

        Cd(8938) = Cd(8935) - 
     &   4*M_3*MTy*Cd(8937)*UCStopC(1,3) - 2*Cd(4581)*UStop2(1,1)

        Cd(8948) = Cd(48)**3*Cd(8941) + 
     &   MGl2*(Cd(48)*Cd(8942) - 4*MGl2*Cd(8947))

        Cd(9110) = 2*M_3C*MTy2*Cd(4786)*UCStop(1,3)*
     &    UStop2(1,1) + Cd(9109)*UStopC(1,1) - 
     &   MTy*Cd(9108)*UUStopC(1,3)

        Cd(9129) = (Cd(44)**3*Cd(2134) - MGl2*Cd(2142))*
     &    UStop(1,1) + 2*M_3*MTy*Cd(9128)*UStop(1,2)

        Cd(9149) = -(M_3*Cd(2218)*MStop2(3)*UStop(1,2)) + 
     &   MTy*(Cd(2216)*UStop(1,1) + 
     &      2*Cd(9148)*MStop2(3)*UStop2(1,2))

        Cd(9161) = -(M_3*Cd(2271)*MStop2(3)*UStop(1,2)) + 
     &   MTy*(Cd(2269)*UStop(1,1) + 
     &      2*Cd(9160)*MStop2(3)*UStop2(1,2))

        Cd(9164) = -(M_3*Cd(2291)*UStop(1,2)) + 
     &   MTy*(Cd(2289)*UStop(1,1) + 2*Cd(9163)*UStop2(1,2))

        Cd(9532) = Cd(48)**2*Cd(9527) - MGl2*Cd(9531) + 
     &   4*MTy2*Cd(2008)*MStop2(3)*UStop2(1,2)

        Cd(9786) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(S2A*Cd(9785) + 
     &        MStop2(3)*(SA2*Cd(533)*MStop2(3) + 
     &           2*CA2*Cd(518)*UStop2(1,1)*UStop2(1,2))))/
     &    (SB2*SW2*MStop2(3)**2)

        Cd(9838) = -2*SA2*Cd(772) + S2A*Cd(9837) - 
     &   8*CA2*Cd(13)**2*Cr(75)*MStop2(4)

        Cd(9870) = -2*MTy*Cd(9868) + At*Cd(9869) + 
     &   8*AtC*Cd(795)*Cd(833)*Cr(75)

        Cd(996) = 4*Cd(972)*MStop2(3)**2 - 
     &   2*AtC*Cd(994)*UCStop(1,3) + 
     &   At*(-Cd(995) + At*Cd(993)*UCStopC(1,3)**2)

        Cd(10075) = (2*Cd(50)*Cd(10073)*MStop2(3) + 
     &     Cd(47)*(-(Cd(9952)*MStop2(2)*MStop2(3)) + 
     &        Cd(50)**2*(2*Cd(10014) - 
     &           Cd(10074)*MStop2(1)*MStop2(3))))/
     &   (Cd(47)*Cd(50)**2*MStop2(3))

        Cd(10079) = (2*Cd(50)**2*Cd(9947)*MGlpTmSt2(1)*
     &      MStop2(3) + Cd(47)*
     &      (-(Cd(9952)*MStop2(2)*MStop2(3)) + 
     &        Cd(50)**2*(2*Cd(10014) - 
     &           Cd(10078)*MStop2(1)*MStop2(3))))/
     &   (Cd(47)*Cd(50)**2*MStop2(3))

        Cd(10084) = MUE**2*Cd(10082)*UCStop(1,3)**2 + 
     &   (2*MUE2*Cd(47)*Cd(50)*Cd(10080) + 
     &      MUEC**2*Cd(10083)*UCStopC(1,3)**2)/MStop2(3) + 
     &   2*MUE2*Cd(9980)*UStop2(1,1)*UStop2(1,2)

        Cd(10133) = Cd(10012) + 
     &   (Cd(47)*(-(Cd(47)*Cd(10013)*MStop2(3)) + 
     &        2*Cd(50)**2*Cd(9947)*Cd(9955)*MGlpTmSt2(1)*
     &         UStop2(1,2)))/(Cd(50)**2*MStop2(3))

        Cd(10216) = (-(Cd(47)**2*Cd(10190)) - 
     &     MGl2*Cd(152)*Cd(2597)*Cd(10191) + 
     &     MStop2(1)*(Cd(152)*Cd(2597)*Cd(10114) + 
     &        Cd(10215)*MStop2(1)))/(Cd(47)**2*MStop2(1)**2)

        Cd(10255) = 4*MGl2*Cd(10247) + 2*MTy*Cd(10250) + 
     &   MTy2*(4*MTy2*Cd(10193) + Cd(10254))

        Cd(10260) = At*(MUEC*Cd(10259)*UCStopC(1,3)**2 + 
     &      MUE*Cd(10258)*UStop2(1,1)) + 
     &   AtC*(MUE*Cd(10259)*UCStop(1,3)**2 + 
     &      MUEC*Cd(10258)*UStop2(1,1))

        Cd(1031) = MUE2*(-2*Cd(957) + 2*Cd(1018)) - 
     &   2*MUE**2*Cd(1029)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(1030)*UCStopC(1,3)**2

        Cd(10527) = At*Cd(10523) - Cd(9948)*Cd(10526) - 
     &   2*Cd(10525)*UCStop(1,3)**2 + 
     &   2*M_3*MGlmT2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(10530) = AtC*Cd(10523) - Cd(9948)*Cd(10529) - 
     &   2*Cd(10528)*UCStopC(1,3)**2 + 
     &   2*M_3C*MGlmT2*MStop2(3)*UStop2(1,1)*UStop2(1,2)

        Cd(10619) = -(Cd(10114)*Cd(10604)) + 
     &   (2*Cd(50)*(Cd(10612) + Cd(50)*Cd(492)*MStop2(3) - 
     &        MGlmT2*Cd(10618)*MStop2(3)))/MStop2(3)

        Cd(10635) = (2*MGlmT2*Cd(50)*
     &      (Cd(10634) + Cd(10149)*MStop2(3)))/MStop2(3) + 
     &   AtC*Cd(10633)*UCStop(1,3) + At*Cd(10603)*UCStopC(1,3)

        Cd(10727) = 2*Cd(50)*(Cd(9948)*Cd(10708) + Cd(10726))*
     &    MStop2(2) + (2*MTy*Cd(13)*Cd(50)**2 + Cd(10711))*
     &    MStop2(3)

        Cd(10740) = Cd(50)*Cd(10735)*UStop2(1,1)*UStop2(1,2) + 
     &   MStop2(2)*(Cd(10733) + 
     &      UStop2(1,2)**2*
     &       (Cd(10739) - MGlmT2*Cd(50)*MStop2(1)*UStop2(1,2)))

        Cd(10796) = AtC*Cd(10793) + 
     &   At*(Cd(10795) + At*UCStopC(1,3)**2) + 
     &   12*Cd(9955)*MStop2(1)*UStop2(1,1)**2

        Cd(10797) = 4*Cd(179)*Cd(251)*Cd(10353) + 
     &   (MGlmT2*Cd(50)*(Cd(10796) + Cd(10355)*MStop2(3)))/
     &    MStop2(3) - 2*Cd(9948)*Cd(10791)*UStop2(1,1)

        Cd(10871) = MTy2*
     &    (2*MTy2*Cd(251)*Cd(10532) + Cd(10870)) - 
     &   (MGl2*Cd(50)*Cd(10830) + 2*MTy*Cd(10847))/MStop2(3)

        Cd(10878) = Cd(10874)*UStop2(1,1)**2 + 
     &   UStop2(1,2)*(Cd(10875)*UStop2(1,1)**2 + 
     &      UStop2(1,2)*(-Cd(10877) + Cd(10611)*UStop2(1,2)))

        Cd(10891) = Cd(50)*(-(Cd(50)*Cd(427)) + Cd(10885)) + 
     &   Cd(9948)*Cd(10890) + MTy2*Cd(251)*Cd(10879)*MStop2(3)

        Cd(10892) = (8*(Cd(10891) + 
     &       MUE2*Cd(50)*(Cd(10878) + 
     &          MStop2(2)*(Cd(10070)*MStop2(1) + MStop2(3)))))/
     &   (Cd(50)**2*MStop2(2)*MStop2(3))

	Cd(11102) = Cd(11252) - Cd(11101)*MStop2(3)

	Cd(11150) = Cd(11255) - Cd(11149)*MStop2(3)

	Cd(11243) = AtC*Cd(11030) + 2*Cd(11045) - Cd(11194)

        Cd(11369) = 4*Cd(11368)*MStop2(2) - 
     &   2*Cd(11359)*UStopC(2,1) + 
     &   AtC*(4*AtC*Cd(640)*MStop2(4)*UCStop(1,3)**2*UStop2(1,2) - 
     &      Cd(11367)*UUStopC(3,4))

        Cd(11454) = 2*MUEC*
     &    (Cd(11453) - MUEC*Cd(11442)*UCStopC(1,3)) + 
     &   4*MUE*Cd(11443)*UCStop(1,3)*UStop2(1,1)*UStopC(2,2) - 
     &   2*Cd(11446)*UUStopC(2,2)

        Cd(11526) = 2*Cd(11369) - 2*Cd(11521) - 
     &   Cd(11523)*UStopC(2,2) + 2*Cd(11525)*UUStopC(2,1)

	Cd(11602) = Cd(11523) + Cd(11601)*MStop2(3)

	Cd(11617) = Cd(11550) + Cd(11609) + Cd(11615)

        Cd(11832) = MGl2**2*Cd(152)*Cd(10579) - 
     &   2*MTy*Cd(11828) + 
     &   (-(MGl2*Cd(11827)) + MTy2*Cd(11831))/MStop2(3)

        Cd(11857) = (Cd(11856) + 5*Cd(50)*Cd(492)*MStop2(3) + 
     &     AtC*MGlmT2*Cd(11851)*MStop2(3) + 
     &     At*MGlmT2*Cd(11852)*MStop2(3))/MStop2(3)

        Cd(11993) = Cd(9948)*Cd(11988) - Cd(47)**2*Cd(11992) + 
     &   Cd(50)**2*(Cd(47)*Cd(11989) - Cd(11990)*MGlpTmSt2(1))

        Cd(11994) = (-2*Cd(47)*Cd(50)*Cd(11986) + 
     &     (Cd(11993) - 4*(Cd(50)**2*MGlpTmSt2(1)*MStop2(1) + 
     &           Cd(47)**2*MGlpTmSt2(2)*MStop2(2)))*MStop2(3))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(12003) = 2*M_3C*Cd(12002)*UCStop(1,3)**2 + 
     &   (2*Cd(9948)*Cd(11999) + At*Cd(47)*Cd(50)*Cd(12001))*
     &    UStop2(1,1) + At*Cd(47)**2*Cd(50)*Cd(9948)*UStop2(1,2)**2

        Cd(12023) = Cd(9948)*Cd(12012) + Cd(47)**2*Cd(12022) + 
     &   Cd(50)**2*(Cd(47)*Cd(12014) - Cd(12018)*MGlpTmSt2(1))

        Cd(12092) = (-Cd(12089) + 
     &     (Cd(50)*(MUEC*Cd(50)*Cd(12090)*MStop2(3)*UCStopC(1,3) + 
     &          MUE*(-2*Cd(47)*Cd(12091) + 
     &             MUE*Cd(50)*Cd(12088)*UCStop(1,3)*UStop2(1,2))))/
     &      (Cd(47)**2*MStop2(3)))/Cd(50)**2

        Cd(12105) = MUEC**2*Cd(12101) + 
     &   (2*MUE2*Cd(12104))/(Cd(50)**2*MStop2(3)) + 
     &   (MUE*UCStop(1,3)*(-(MUE*Cd(47)**2*Cd(12103)*
     &           UCStop(1,3)**2) - 2*Cd(12095)*UStopC(1,2)))/
     &    Cd(47)**2

        Cd(12116) = Cd(47)**2*Cd(12113) + 
     &   Cd(50)*(2*Cd(47)*Cd(50)*Cd(12115) + 
     &      MStop2(1)*(2*Cd(47)*Cd(6096)*Cd(9948) + 
     &         Cd(50)*MStop2(1)*MStop2(3)*UStop2(1,2)))

        Cd(12128) = Cd(47)*
     &    (Cd(47)*(Cd(9948)*Cd(12123) + Cd(12127)) - 
     &      Cd(50)**2*Cd(12122)*Cr(75)) + 
     &   MUE**2*Cd(50)**2*MStop2(1)*MStop2(3)*UCStop(1,3)**2*
     &    UStop2(1,2)

        Cd(12132) = Cd(47)**2*Cd(9948)*Cd(12129) + 
     &   Cd(47)**2*Cd(12131) - 
     &   Cd(50)**2*(Cd(47)*Cd(12122)*Cr(75) + 
     &      MUE**2*Cd(118)*MStop2(3)*UCStop(1,3)**2*UStop2(1,2))

        Cd(12499) = 4*Cd(50)*Cd(11734) + 
     &   (Cd(10114)*Cd(12483) + 4*Cd(12498)*MStop2(2))/
     &    MStop2(2)**2 - (4*Cd(50)*Cd(11863))/MStop2(3)

        Cd(12509) = 4*Cd(12501)*MStop2(3) + 
     &   AtC*Cd(12504)*UCStop(1,3)**2 + 
     &   At*(2*Cd(12508) + At*Cd(12505)*UCStopC(1,3)**2)

        Cd(12536) = (2*Cd(50)*Cd(12535) + 
     &     Cd(9948)*Cd(12531)*MStop2(3) + 
     &     Cd(12532)*MStop2(3)*UCStop(1,3)**2 + 
     &     At*Cd(12529)*MStop2(3)*UStop2(1,1))/MStop2(3)

        Cd(12537) = (2*Cd(9982))/MStop2(2)**2 + 
     &   (Cd(12536)/MStop2(2) + 
     &      (2*Cd(50)*(Cd(12528) + Cd(9961)*MStop2(3)))/MStop2(3))/
     &    Cd(50)**2

        Cd(12547) = (2*Cd(50)*Cd(12545) + 
     &     Cd(9948)*Cd(12546)*MStop2(3) + 
     &     Cd(12544)*MStop2(3)*UCStopC(1,3)**2 + 
     &     AtC*Cd(12529)*MStop2(3)*UStop2(1,1))/MStop2(3)

        Cd(12548) = (2*Cd(9988))/MStop2(2)**2 + 
     &   (Cd(12547)/MStop2(2) + 
     &      (2*Cd(50)*(Cd(12543) + Cd(9999)*MStop2(3)))/MStop2(3))/
     &    Cd(50)**2

	Cd(12587) = -(Cd(50)**2*Cd(164)) + Cd(47)**2*Cd(251) - Cd(12045)

        Cd(12606) = (4*MGl2*Cd(12592) + Cd(12605) - 
     &     4*Cd(50)**2*Cd(12590)*MStop2(1) - 
     &     4*Cd(47)**2*Cd(12591)*MStop2(2))/(Cd(47)**2*Cd(50)**2)

        Cd(12624) = (Cd(12616) + Cd(12622) + 
     &     Cd(50)**2*Cd(12623)*MGlpTmSt4(1) + 
     &     2*MUE2*Cd(47)**2*Cd(9981)*MStop2(2)**2*UStop2(1,1))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12660) = (-4*MGl2*Cd(12592) + Cd(12659) + 
     &     4*Cd(50)**2*Cd(12590)*MStop2(1) + 
     &     4*Cd(47)**2*Cd(12591)*MStop2(2))/(Cd(47)**2*Cd(50)**2)

        Cd(12671) = Cd(9948)*Cd(11988) - Cd(47)**2*Cd(11992) + 
     &   Cd(50)**2*(Cd(47)*Cd(11989) - Cd(11990)*MGlpTmSt2(1))

        Cd(12672) = (2*Cd(47)*Cd(50)*Cd(11986) + 
     &     (-Cd(12671) + 4*Cd(50)**2*MGlpTmSt2(1)*MStop2(1) + 
     &        4*Cd(47)**2*MGlpTmSt2(2)*MStop2(2))*MStop2(3))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(12673) = Cd(9948)*Cd(12012) + Cd(47)**2*Cd(12022) + 
     &   Cd(50)**2*(Cd(47)*Cd(12014) - Cd(12018)*MGlpTmSt2(1))

        Cd(12676) = (Cd(12616) + Cd(12622) + 
     &     Cd(50)**2*Cd(12623)*MGlpTmSt4(1) + 
     &     2*MUE2*Cd(47)**2*Cd(9981)*MStop2(2)**2*UStop2(1,1))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(13151) = 1/512.D0*
     &   (512*Cd(12701)*Cd(13149) + 
     &     1/Pi**4*(Cd(13150)*MStop2(1)**2*MStop2(3) + 
     &         12*MStop2(1)*MStop2(3)**2 - 
     &         4*Cd(12700)*Cd(13139)*UStop2(1,1)*UStop2(1,2)*
     &          UUStopC(1,2))/(MStop2(1)**2*MStop2(3)))

        Cd(13273) = 1/32.D0*
     &   (1/Pi**4*(MGl2*MTy2*Cd(12836)**2*Cr(45))/
     &       (TB2*MGlpTmStxGlT4(1)**2) + 
     &     1/Pi**4*(Cd(13272) + 32*Pi**4*(Cd(13271)*MStop2(3)))/
     &       (TB2*MStop2(3)**2) + (32*Cd(13232))/MGlpTmStxGlT4(1))

        Cd(13412) = 1/512.D0*
     &   (512*Cd(12701)*Cd(13410) + 
     &     1/Pi**4*(Cd(13411) + 
     &        (4*(Cd(12700)*Cd(13383) + 3*MStop2(2)*MStop2(3)**2))/
     &         (MStop2(2)**2*MStop2(3))))

        Cd(1366) = MUE2*(4*Cd(1286) - Cd(1354)) - 
     &   MUE**2*Cd(1364)*UCStop(1,3)**2 + 
     &   MUEC**2*Cd(1365)*UCStopC(1,3)**2

        Cd(13956) = MTy2**2*Cd(13846)*Cd(13955) + 
     &   MSq2Diff(1,2)*(Cd(13954) + 
     &      MTy2*(-(MTy2**2*Cd(13846)) + 
     &         Cd(138)*Cd(13846)*MStop2(2)))

        Cd(14088) = MTy*Cd(14087) + MTy2*Cd(14086)*MStop2(1) + 
     &   MGl2*UStop2(1,2)*(Cd(14083) + 
     &      MGl2*TB*MSq2Diff(1,1)*UStop2(1,2))

        Cd(14141) = ((TB*Cd(14124) + 
     &        TB2*Cd(47)*Cd(424)*Cd(9981)*MGlpTmSt2(1))*MStop2(1)+
     &       Cd(9950)*(-Cd(14140) + 
     &        MTy2**2*TB*Cd(9969)*MGlpTmSt2(1)*UStop2(1,1)))/
     &   (TB*TB2*Cd(47)*Cd(9950)*MSq2Diff(1,1))

        Cd(14195) = (TB*Cd(14185)*MSq2Diff(1,2)*MStop2(2) - 
     &     Cd(10103)*(Cd(14194) + 
     &        TB*MGlpTmSt2(2)*
     &         (Cd(14186) - Cd(13841)*MSq2Diff(1,2)*UStop2(1,1))))/
     &   (TB*TB2*Cd(50)*Cd(10103)*MSq2Diff(1,2))

        Cd(14246) = (2*Cd(50)*Cr(5)*Cr(75)*MSq2Diff(1,2)**2 + 
     &      Cd(14245)*MStop2(3))/
     &    (TB2*Cd(50)*MSq2Diff(1,1)*MSq2Diff(1,2)**2*MStop2(3)) - 
     &   (Cd(14234)*(MUE*UCStop(1,3) + MUEC*UCStopC(1,3)))/
     &    (TB*Cd(50)*MSq2Diff(1,2)**2)

        Cd(1429) = Cd(424)*Cd(1420) - 
     &   2*MTy*Cd(13)*Cd(1428)*MStop2(3)**2 - 
     &   Cd(1427)*(AtC*MUE*UCStop(1,3)**2 + 
     &      At*MUEC*UCStopC(1,3)**2)

	Cd(1465) = Cd(1383) + Cd(1461) + AtC*MTy*Cd(1462)*UCStop(1,3)

        Cd(2658) = Cd(2649) - Cd(2657) + 
     &   Cd(2644)*UStop2(1,1)**2*UStopC(1,1) + 
     &   Cd(2654)*UUStopC(1,1)

        Cd(2734) = Cd(2727) + Cd(2733) - 
     &   Cd(2720)*UStop2(1,1)**2*UStopC(1,1) - 
     &   Cd(2730)*UUStopC(1,1)

        Cd(3044) = MGl2**2*(-(MGl2*Cd(2964)) + Cd(2967)) - 
     &   MTy2*Cd(3042) - Cd(3043)*MStop2(1)*MStop2(3)

        Cd(3146) = 2*Cd(3125) + 
     &   2*M_3*MTy*Cd(3145)*UCStopC(1,3) - 
     &   Cd(3143)*UStop2(1,1) + Cd(3136)*UStop2(1,2)

        Cd(3198) = -(Cd(3197)*UCStop(1,3)) + 
     &   Cd(3195)*UStop2(1,1) + 
     &   4*MStop2(2)*(-2*M_3*MTy*Cd(3116)*UStop2(1,1)**2 + 
     &      Cd(50)*Cd(3196)*UStopC(1,2))

        Cd(3219) = 2*Cd(3214)*MStop2(2) - 
     &   M_3*MTy*Cd(3218)*UCStopC(1,3) - 
     &   MGl2*Cd(3211)*UStop2(1,1) - Cd(3207)*UStop2(1,2)

        Cd(3550) = 8*MTy*Cd(3549)*MStop2(3) - 
     &   4*Cd(3545)*UStop2(1,1) + 
     &   4*Cd(3534)*UCStop(1,3)*UStop2(1,2) - 
     &   4*Cd(47)*Cd(3548)*UStopC(1,2)

        Cd(3744) = MGl2**2*(-(MGl2*Cd(2964)) + Cd(2967)) - 
     &   MTy2*Cd(3042) + 3*Cd(2970)*MStop2(1)**2*MStop2(3)

        Cd(3818) = 2*MTy*Cd(3817)*MStop2(3) - 
     &   Cd(3813)*UStop2(1,1) + Cd(3803)*UCStop(1,3)*UStop2(1,2) + 
     &   Cd(50)*Cd(3816)*UStopC(1,2)

        Cd(4077) = M_3C*Cd(4064)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(4072) + Cd(4073)*UStop2(1,1)) + 
     &   MTy*Cd(4076)*UCStop(1,3)*UStop2(1,2)

        Cd(4104) = Cd(4087) + Cd(4103) + 
     &   MTy*(MTy2*Cd(4078)*UStop2(1,1)**2*UStopC(1,1) + 
     &      Cd(4094)*UUStopC(1,1))

        Cd(4175) = M_3C*Cd(4162)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(4170) + Cd(4171)*UStop2(1,1)) - 
     &   MTy*Cd(4174)*UCStop(1,3)*UStop2(1,2)

        Cd(4211) = M_3C*Cd(4162)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(4170) + Cd(4210)*UStop2(1,1)) - 
     &   MTy*Cd(4174)*UCStop(1,3)*UStop2(1,2)

        Cd(4237) = 2*Cd(4223) + 
     &   2*M_3*MTy*Cd(4236)*UCStopC(1,3) - 
     &   Cd(4234)*UStop2(1,1) - Cd(4229)*UStop2(1,2)

        Cd(4321) = 4*Cd(4223)*MStop2(1) + 
     &   M_3*MTy*Cd(4320)*UCStopC(1,3) - Cd(4318)*UStop2(1,1) - 
     &   Cd(4313)*UStop2(1,2)

        Cd(4434) = 2*M_3*MStop2(1)**5 - 
     &   M_3C*Cd(4402)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(4419) - Cd(4407)*UStop2(1,1))

        Cd(4657) = 4*Cd(4569)*MStop2(2) - 
     &   M_3*MTy*Cd(4656)*UCStopC(1,3) - Cd(4653)*UStop2(1,1) - 
     &   Cd(4648)*UStop2(1,2)

        Cd(4771) = 2*M_3*MStop2(2)**5 + 
     &   M_3C*Cd(4736)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(4756) + Cd(4741)*UStop2(1,1))

        Cd(4912) = 2*M_3*MStop2(2)**5 + 
     &   M_3C*Cd(4736)*UCStop(1,3)**2 + 
     &   UStop2(1,1)*(Cd(4756) - Cd(4911)*UStop2(1,1))

        Cd(5011) = Cd(50)*Cd(4952) - 
     &   2*M_3*MTy*Cd(5010)*UCStopC(1,3) - 
     &   Cd(5004)*UStop2(1,1) - Cd(47)*Cd(4979)*UStop2(1,2)

        Cd(5382) = 2*Cd(44)*Cd(48)*Cd(5352) + 2*Cd(5369) + 
     &   8*MTy2*Cd(5377) - MGl2*Cd(5381)

        Cd(5469) = 2*Cd(44)*Cd(48)*Cd(5435) + 
     &   MGl2*(-2*Cd(5450) + MGl2*Cd(5468)) - 
     &   8*MTy2*Cd(5463)*MStop2(3)**2*UStop2(1,2)

        Cd(6109) = 2*M_3*MTy*Cd(6107) - 
     &   (-2*MGl2*Cd(6100) + Cd(6105))*UCStop(1,3) + 
     &   (2*MGl2**2*Cd(6108)*UCStop(1,3))/(Cd(47)**2*Cd(50)**2)

        Cd(6120) = (Cd(47)*Cd(50)*
     &      (16*M_3*MTy*Cd(295)*UStop2(1,1)**2 + 
     &        Cd(47)*Cd(50)*MStop2(3)*
     &         (-2*M_3C*MTy*Cd(6097)*UCStop(1,3)**2 + 
     &           Cd(6109)*UStop2(1,1))) - Cd(6119)*UStopC(1,2))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(7045) = Cd(50)*Cd(7027) + 
     &   2*M_3*MTy*Cd(7044)*MStop2(3)*UCStopC(1,3) + 
     &   Cd(7042)*UStop2(1,1) + Cd(47)*Cd(7037)*UStop2(1,2)

        Cd(7371) = Cd(44)*Cd(48)*Cd(7352) + 
     &   MGl2*(-2*Cd(7366) + MGl2*Cd(7370)) + 
     &   4*MTy2*Cd(5463)*MStop2(3)**2*UStop2(1,2)

        Cd(7541) = MGl2*(MGl2*(2*MGl2**2 - Cd(7529)) + 
     &      4*Cd(7531)) + Cd(7535)

        Cd(7853) = 2*M_3*MTy*Cd(7852) - 
     &   (MGl2**2*Cd(6108)*UCStop(1,3))/(Cd(47)**2*Cd(50)**2) + 
     &   (-2*MGl2*Cd(7845) + Cd(7850))*UCStop(1,3)

        Cd(7863) = (Cd(47)*Cd(50)*
     &      (8*M_3*MTy*Cd(295)*UStop2(1,1)**2 + 
     &        Cd(47)*Cd(50)*MStop2(3)*
     &         (2*M_3C*MTy*Cd(7843)*UCStop(1,3)**2 + 
     &           Cd(7853)*UStop2(1,1))) - Cd(7862)*UStopC(1,2))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(8558) = (-(MGl2*(MGl2*Cd(8551) + 2*Cd(8552))) + 
     &      Cd(8555))*UStop(1,1) + 2*M_3*MTy*Cd(8557)*UStop(1,2)

        Cd(8701) = -2*MTy*Cd(8700)*MStop2(3) + 
     &   Cd(8697)*UStop2(1,1) + Cd(8693)*UCStop(1,3)*UStop2(1,2) + 
     &   Cd(47)*Cd(8699)*UStopC(1,2)

        Cd(8933) = -2*MTy*Cd(8932)*MStop2(3) + 
     &   Cd(8930)*UStop2(1,1) + Cd(8926)*UCStop(1,3)*UStop2(1,2) - 
     &   Cd(50)*Cd(8931)*UStopC(1,2)

        Cd(909) = 2*Cd(907) - AtC**2*Cd(879)*UCStop(1,3)**2 + 
     &   At*(-2*Cd(908) + At*Cd(893)*UCStopC(1,3)**2)

        Cd(9978) = -(Cd(9965)/(Cd(47)*Cd(50))) + 
     &   (Cd(9952)*Cd(9957)*MStop2(2))/Cd(50)**2 + 
     &   Cd(9977)/MStop2(3) + 
     &   (MStop2(1)*(At*Cd(47)**2*Cd(9955)*Cd(9956) + 
     &        Cd(9951)*Cd(9957)*MStop2(3)**2))/
     &    (Cd(47)**2*MStop2(3)**2)

        Cd(10008) = -(Cd(10000)/(Cd(47)*Cd(50))) + 
     &   (Cd(9952)*Cd(9993)*MStop2(2))/Cd(50)**2 + 
     &   Cd(10007)/MStop2(3) + 
     &   (MStop2(1)*(AtC*Cd(47)**2*Cd(9955)*Cd(9956) + 
     &        Cd(9951)*Cd(9993)*MStop2(3)**2))/
     &    (Cd(47)**2*MStop2(3)**2)

        Cd(10035) = (4*Cd(10014)*Cd(10017) + 
     &     (MStop2(3)*(2*Cd(47)*Cd(50)*Cd(10024) - 
     &          MStop2(3)*(Cd(10034) + 
     &             (-(Cd(50)**2*Cd(10012)*MStop2(1)) + 
     &                Cd(47)**2*Cd(10013)*MStop2(2))*UCStop(1,3))))
     &       /(Cd(47)**2*Cd(50)**2))/MStop2(3)**2

        Cd(10099) = MUE**2*Cd(10090)*UCStop(1,3) + 
     &   MUEC**2*Cd(10098)*UCStopC(1,3) + 
     &   (2*MUE2*Cr(5)*(-4*Cd(10014)*Cd(10016) + 
     &        (MStop2(3)*(2*Cd(47)*Cd(50)*Cd(10093) + 
     &             Cd(10031)*MStop2(3)*UStop2(1,1)*UStop2(1,2)))/
     &         (Cd(47)**2*Cd(50)**2)))/MStop2(3)**2

        Cd(10109) = Cd(47)**2*Cd(10102) + 
     &   Cd(10106)*UStop2(1,1) + 
     &   UStop2(1,2)*(Cd(10107) + 
     &      UStop2(1,2)*(Cd(10108) - Cd(9980)*UStop2(1,2)))

        Cd(10375) = (4*Cd(9964))/Cd(47) + 
     &   (Cd(10368) - Cd(10374)*MStop2(1))/MStop2(1)**2 + 
     &   (2*MGlmT2*Cd(10369) - 2*Cd(10370)*MStop2(1))/
     &    (Cd(47)*MStop2(3))

        Cd(10381) = (4*Cd(9996))/Cd(47) + 
     &   (2*Cd(9988) - Cd(10380)*MStop2(1))/MStop2(1)**2 + 
     &   (2*MGlmT2*Cd(10376) - 2*Cd(10377)*MStop2(1))/
     &    (Cd(47)*MStop2(3))

        Cd(10434) = MUE*Cd(10375) + MUEC*Cd(10381) + 
     &   (2*MGl2*Cd(10392) - MTy2*Cd(10406) + 2*MTy*Cd(10433))/
     &    (Cd(47)**2*MStop2(1)**2*MStop2(3))

        Cd(10516) = -(AtC*Cd(10515)) + 
     &   Cd(164)*Cd(10491)*MStop2(3) + 
     &   AtC**2*Cd(10485)*UCStop(1,3) + 
     &   Cd(47)*(Cd(10505) - 
     &      2*At*Cd(47)*Cd(9955)*UCStopC(1,3)*UStop2(1,1))

        Cd(10517) = MTy2*(MTy2*Cd(164)*Cd(10473) + Cd(10484)) + 
     &   Cd(47)*Cd(10295)*MStop2(1) + 
     &   (-(MGl2*Cd(47)*Cd(10301)) + 2*MTy*Cd(10516))/
     &    (MStop2(1)*MStop2(3))

        Cd(10518) = (8*Cd(47)*Cd(10282) + 8*Cd(10517) + 
     &     (-32*MTy2**3*Cd(164)*Cd(10193) + Cd(47)**2*Cd(10470))/
     &      MStop2(1) - (8*MGlmT2*Cd(47)*Cd(10294))/MStop2(3))/
     &   Cd(47)**2

        Cd(10569) = 2*MTy2**2*Cd(13)*Cd(164) + 
     &   (MGl2*Cd(13)*Cd(10546) - MUE*Cd(10554) + MTy2*Cd(10560) - 
     &      MUEC*Cd(10568))/MStop2(3)

        Cd(10800) = (8*Cd(10797))/(Cd(50)**2*MStop2(2)) + 
     &   (8*(Cd(10787) - MGl2*Cd(10790) + 
     &        4*Cd(251)*Cd(9969)*MStop2(2)))/Cd(50)**2 - 
     &   Cd(10799)/MStop2(3)

        Cd(1090) = -4*Cd(1084)*MStop2(3)**2 - 
     &   2*AtC*Cd(1087)*UCStop(1,3) + 
     &   At*(2*AtC*Cd(1083) + 
     &      4*MTy*Cd(1088)*MStop2(3)*UCStopC(1,3) + 
     &      At*Cd(1089)*UCStopC(1,3)**2)

	Cd(11050) = Cd(11243) - Cd(11049)*MStop2(3)

        Cd(11606) = 2*Cd(11369) - 2*Cd(11521) - 
     &   2*Cd(11600)*MStop2(3) - Cd(11602)*UStopC(2,2) + 
     &   Cd(11605)*UUStopC(2,1)

        Cd(11776) = Cd(47)*Cd(11769) - AtC*Cd(11775) + 
     &   MGl2*Cd(152)*(Cd(152)*Cd(10141) + Cd(10202))*MStop2(3) + 
     &   AtC**2*Cd(11758)*UCStop(1,3)

        Cd(11777) = MGl2*
     &    (-(MGl2*Cd(152)*Cd(11733)) + Cd(11737)) - 
     &   MTy2*Cd(11757) + (2*MTy*Cd(11776) + 
     &      Cd(47)*Cd(11730)*MStop2(1)**2)/(MStop2(1)*MStop2(3))

        Cd(11892) = (8*Cd(11741))/Cd(50) - 
     &   (8*(MGlmT2*Cd(50)*Cd(11891) + 
     &        MGl2*Cd(179)*
     &         (-4*MGl2*Cd(10193) + Cd(11875)*MStop2(2))))/
     &    (Cd(50)**2*MStop2(2)) - Cd(11885)/MStop2(3)

        Cd(11970) = (2*(MGl2*Cd(585)*Cd(11966) - 
     &       4*MTy2*Cd(11967)*Cr(5) - 
     &       Cd(50)**2*Cd(152)*Cd(585)*MGlpTmSt2(1)*MStop2(1) + 
     &       Cd(47)*(Cd(50)*Cd(11964)*Cd(11969) + 
     &          Cd(47)*Cd(179)*Cd(585)*MGlpTmSt2(2)*MStop2(2))))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12041) = (2*Cd(47)*Cd(50)*Cd(12040) + 
     &     (Cd(12023) + 4*Cd(50)**2*Cd(10471)*MGlpTmSt2(1)*
     &         MStop2(1) + 
     &        4*Cd(47)**2*(-(MGl2*Cd(11963)) + 
     &           Cd(10471)*MGlpTmSt2(2)*MStop2(2)))*MStop2(3))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(12393) = -(AtC*Cd(12390)*UCStop(1,3)**2) - 
     &   At*(At*Cd(12392) + 8*M_3*Cd(12385)*MStop2(3))*
     &    UCStopC(1,3)**2 + 4*At*AtC*Cd(47)**2*UStop2(1,2)**2 - 
     &   36*Cd(47)**2*MStop2(3)*UStop2(1,2)**2

        Cd(12403) = (4*MGl2*Cd(10193)*Cd(12371) - Cd(12401) + 
     &      4*Cd(12402)*MStop2(1))/MStop2(1)**2 + 
     &   (4*Cd(47)*(Cd(12398) + Cd(11741)*MStop2(3)))/MStop2(3)

        Cd(12422) = (4*Cd(9964))/Cd(47) + 
     &   (Cd(10368) + Cd(12421)*MStop2(1))/MStop2(1)**2 - 
     &   (2*(Cd(11779) + Cd(10370)*MStop2(1)))/(Cd(47)*MStop2(3))

        Cd(12431) = (4*Cd(9996))/Cd(47) + 
     &   (2*Cd(9988) + Cd(12430)*MStop2(1))/MStop2(1)**2 - 
     &   (2*(Cd(11787) + Cd(10377)*MStop2(1)))/(Cd(47)*MStop2(3))

        Cd(12665) = (MGl2*Cd(585)*Cd(11966) - 
     &     4*MTy2*Cd(11967)*Cr(5) - 
     &     Cd(50)**2*Cd(152)*Cd(585)*MGlpTmSt2(1)*MStop2(1) + 
     &     Cd(47)*(Cd(50)*Cd(11964)*Cd(11969) + 
     &        Cd(47)*Cd(179)*Cd(585)*MGlpTmSt2(2)*MStop2(2)))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(12674) = (2*Cd(47)*Cd(50)*Cd(12040) + 
     &     (Cd(12673) + 4*Cd(50)**2*Cd(10471)*MGlpTmSt2(1)*
     &         MStop2(1) + 
     &        4*Cd(47)**2*(-(MGl2*Cd(11963)) + 
     &           Cd(10471)*MGlpTmSt2(2)*MStop2(2)))*MStop2(3))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

        Cd(1411) = 2*(-2*Cd(1386)*MStop2(3) - 
     &     AtC**2*Cd(1393)*UCStop(1,3)**2 + 
     &     At*(AtC*Cd(1372) - 
     &        2*MTy*Cd(1402)*MStop2(3)*UCStopC(1,3) + 
     &        At*Cd(1410)*UCStopC(1,3)**2))

        Cd(2864) = 4*Cd(2766)*MStop2(1) + 
     &   M_3*MTy*Cd(2863)*UCStopC(1,3) - 
     &   (MGl2*Cd(2857) + 2*Cd(2860))*UStop2(1,1) - 
     &   Cd(2855)*UStop2(1,2)

        Cd(301) = MGlmT2*(2*MGlmT2**2 - Cd(298)*MStop2(1)) - 
     &   Cd(299)*MStop2(2) + Cd(164)*MStop2(2)**2 + 
     &   Cd(300)*UStop2(1,2)

        Cd(4200) = AtC*(4*Cd(4175) - AtC*Cd(4135)*UCStop(1,3)) + 
     &   At*(4*Cd(4161) - At*Cd(4143)*UCStopC(1,3) + 
     &      2*AtC*Cd(4132)*UStop2(1,1))

        Cd(4885) = 2*Cd(4883)*MStop2(2) - 
     &   M_3*MTy*(2*Cd(48)**2*Cd(4649) + MGl2*Cd(4655))*
     &    UCStopC(1,3) - Cd(4882)*UStop2(1,1) - 
     &   2*Cd(4884)*UStop2(1,2)

        Cd(5947) = MGl2*(MGl2*
     &       (MGl2*(2*MGl2 - 6*Cd(198)) + Cd(5934)) - 
     &      4*Cd(198)*Cd(5935)) + Cd(5938)

        Cd(6161) = 2*Cd(44)*Cd(48)*Cd(6127) + 
     &   MGl2*(-2*Cd(6138) + MGl2*(2*Cd(6132) - MGl2*Cd(6160))) - 
     &   MTy2*Cd(6156)*UStop2(1,2)

        Cd(6595) = 2*M_3*Cd(44)*Cd(48)*Cd(6545)*UCStopC(1,3) + 
     &   MTy*Cd(44)**2*Cd(48)**2*Cd(6514)*UStop2(1,1) - 
     &   MGl2*MTy*Cd(6590)*UStop2(1,1) - 
     &   2*(M_3*MGl2*Cd(6594)*UCStopC(1,3) + 
     &      MTy*Cd(6586)*UStop2(1,2)*UStopC(1,1))

        Cd(6728) = M_3C*MTy*Cd(5280)*UCStop(1,3)**2 + 
     &   Cd(6709)*UStop2(1,1) + 
     &   Cd(47)*(Cd(47)*Cd(6680)*MStop2(3)*UCStop(1,3) + 
     &      4*M_3*MTy*Cd(50)*Cd(6679)*UStop2(1,1)**2) + 
     &   Cd(50)*Cd(6727)*UStop2(1,2)*UStopC(1,2)

	Cd(7046) = Cd(5011) - 2*Cd(7045)*MStop2(3)

        Cd(728) = -(Cd(724)*Conjugate(dY33eps(0))) + 
     &   2*Cd(711)*dY33eps(0) + Cd(727)*MStop2(3) + 
     &   Cd(692)*dY33eps(0)*UStopC(2,1) - 
     &   Cd(689)*Cd(693)*Cd(719)*UUStopC(2,1)

        Cd(731) = -(Cd(724)*Conjugate(dY33eps(1))) + 
     &   2*Cd(711)*dY33eps(1) + Cd(730)*MStop2(3) + 
     &   Cd(692)*dY33eps(1)*UStopC(2,1) - 
     &   Cd(689)*Cd(693)*Cd(729)*UUStopC(2,1)

        Cd(8254) = 2*M_3C*MTy*Cd(7228)*UCStop(1,3)**2 + 
     &   Cd(8251)*UStop2(1,1) + 
     &   Cd(47)*(Cd(47)*Cd(8242)*MStop2(3)*UCStop(1,3) + 
     &      4*M_3*MTy*Cd(50)*Cd(6679)*UStop2(1,1)**2) + 
     &   Cd(50)*Cd(8253)*UStop2(1,2)*UStopC(1,2)

        Cd(8350) = 2*Cd(44)*Cd(48)*Cd(7155) + 
     &   MGl2*(MGl2*(MGl2*(-2*MGl2**2 + Cd(7148)) - 4*Cd(7149)) + 
     &      Cd(8349))

        Cd(9268) = 2*M_3C*MTy*Cd(9236) + 2*Cd(9267) + 
     &   Cd(9225)*MStop2(3)**2 + 
     &   4*MTy2*Cd(9194)*MStop2(3)*UStop2(1,1)**2*UStopC(1,1) - 
     &   Cd(9214)*UUStopC(1,1)

        Cd(9462) = Cd(9454) + 2*M_3C*MTy*Cd(9456) - Cd(9461) - 
     &   2*Cd(96)**2*Cd(9446)*UStop2(1,1)**2*UStopC(1,1) + 
     &   Cd(9450)*UUStopC(1,1)

        Cd(9712) = -2*M_3C*MTy*Cd(9663) + 2*Cd(9711) + 
     &   Cd(9641)*MStop2(3) + 
     &   4*MTy2*Cd(9551)*UStop2(1,1)**2*UStopC(1,1) + 
     &   Cd(9626)*UUStopC(1,1)

        Cd(10639) = 2*Cd(10629) - MTy2*Cd(10632) + 
     &   (-2*MGl2*Cd(585)*Cd(10627) + Cd(10638) + 
     &      Cd(10635)*MStop2(2))/MStop2(2)**2 + 
     &   (2*Cd(50)*Cd(10630))/MStop2(3)

        Cd(10663) = 8*MGl2*Cd(10641) - 4*Cd(10643) + 
     &   (Cd(10646) + 16*Cd(251)*Cd(9969)*MStop2(2)**3 + 
     &      (4*MStop2(2)*(-Cd(10662) + 
     &           Cd(50)*Cd(10654)*MStop2(2)))/MStop2(3))/
     &    MStop2(2)**2

        Cd(1173) = 4*Cd(1167)*MStop2(3)**2 + 
     &   AtC*UCStop(1,3)*(4*MTy*Cd(1169) + 
     &      AtC*Cd(1170)*UCStop(1,3)) - 
     &   At**2*Cd(1172)*UCStopC(1,3)**2 - 
     &   2*At*(AtC*Cd(1168) + 
     &      2*MTy*Cd(1171)*MStop2(3)*UCStopC(1,3))

        Cd(12348) = 8*Cd(50)**2*Cd(492) + 
     &   AtC*Cd(12347)*UCStop(1,3)**2 + 
     &   At*(At*Cd(251)*Cd(2681) + 8*M_3*Cd(4187))*
     &    UCStopC(1,3)**2 + 
     &   8*M_3*AtC*Cd(4187)*UStop2(1,1)*UStop2(1,2) + 
     &   2*At*Cd(12346)*UStop2(1,1)*UStop2(1,2)

        Cd(12462) = -(MUE*Cd(12422)) - MUEC*Cd(12431) + 
     &   (MTy2*Cd(12444))/Cd(47)**2 - 
     &   (4*Cd(13)*Cd(12423))/(MTy*Cd(47)**2*MStop2(1)) - 
     &   (MGl2*Cd(12426) + 2*MTy*Cd(12461))/
     &    (Cd(47)**2*MStop2(1)**2*MStop2(3))

        Cd(14181) = (MUE*Cd(14176) + MUEC*Cd(14180))/TB + 
     &   (2*Cd(47)*Cd(14177)*MSq2Diff(1,1) + 
     &      MStop2(3)*(Cd(14179) + 
     &         Cd(14142)*(MGl2*Cd(9969) - MStop2(1))*UStop2(1,1)))/
     &    (TB2*Cd(47)*MStop2(3))

        Cd(1783) = (MGl2*
     &       (Cd(48)*Cd(1558) + 
     &         MGl2*(MGl2*(4*MGl2 - Cd(1556)) - Cd(1559)*Cd(1560)))
     &        + Cd(1781))*UStop(1,1) + 
     &   4*M_3*MTy*Cd(1782)*MGlpTmSt2(2)*UStop(1,2)

        Cd(4546) = 4*Cd(4223)*MStop2(1) + 
     &   M_3*MTy*(Cd(44)**2*Cd(4314) - MGl2*Cd(4319))*
     &    UCStopC(1,3) + 2*MGl2**3*Cd(4291)*UStop2(1,1) - 
     &   Cd(4544)*UStop2(1,1) - 2*Cd(4545)*UStop2(1,2)

        Cd(7932) = MGl2*(MGl2*
     &       (MGl2*(MGl2*(4*MGl2 - 2*Cd(5861)) + 12*Cd(7913)) - 
     &         Cd(7915)) - 2*Cd(7918)) + Cd(7925)

        Cd(8225) = MTy*Cd(8131)*MStop2(3) + 
     &   2*M_3*(-(Cd(44)*Cd(48)*Cd(6545)) + MGl2*Cd(6594))*
     &    UCStopC(1,3) + MTy*
     &    (-(Cd(44)**2*Cd(48)**2*Cd(6514)) + MGl2*Cd(6590))*
     &    UStop2(1,1) + Cd(8224)*UStopC(1,1)

        Cd(8379) = Cd(44)**2*Cd(48)**2*Cd(2048) + 
     &   MGl2*(-2*Cd(44)*Cd(48)*Cd(7155) + 
     &      MGl2*(MGl2*(MGl2*(2*MGl2**2 - Cd(7148)) + 
     &            4*Cd(7149)) - Cd(8349)))

        Cd(8462) = 4*Cd(44)**2*Cd(48)**2*Cd(198) + 
     &   MGl2*(-Cd(8409) + 
     &      MGl2*(MGl2*(MGl2*(2*MGl2 - 4*Cd(198)) + Cd(8455)) + 
     &         4*Cd(89)*MStop2(3)**2))

        Cd(8499) = MGlmT2*Cd(8496) - Cd(8498)*MStop2(2) + 
     &   Cd(164)*MStop2(2)**2*
     &    (2*Cd(8497) - 4*MGlpT2*MStop2(2) + MStop2(2)**2) + 
     &   Cd(302)*Cd(8494)*UStop2(1,2)**2

        Cd(8522) = MStop2(1)**2*
     &    (Cd(8518) - Cd(8514)*MStop2(1)) - 
     &   MStop2(1)**4*MStop2(2)**2*MStop2(3) - 
     &   (-2*MGlmT2**3*Cd(8509) + Cd(8521) + Cd(8512)*MStop2(1))*
     &    MStop2(3)**2

        Cd(8562) = -(MGl2**4*Cd(8560)) + 
     &   (MGl2*(2*MGl2*(MGl2*(MGl2**2 + Cd(8190)) - Cd(8191)) + 
     &         Cd(8192)) + Cd(8561))*MStop2(3)

        Cd(10267) = AtC*Cd(179)*Cd(2681)*Cd(10263) + 
     &   At*Cd(179)*Cd(2681)*Cd(10266) - 
     &   2*(MTy*Cd(10262)*(MUE*UCStop(1,3) + MUEC*UCStopC(1,3)) + 
     &      2*MTy2*Cd(10261)*
     &       (M_3C*MUE*UCStop(1,3)**2 + 
     &         M_3*MUEC*UCStopC(1,3)**2 + 
     &         Cd(10224)*UStop2(1,1)*UStop2(1,2)))

        Cd(10336) = (8*MGl2*Cd(585)*Cd(10305) + 
     &     8*Cd(47)*Cd(10325)*MStop2(1)**2 - 
     &     (4*Cd(10319) + Cd(47)**2*Cd(10335)*MStop2(1) + 
     &        4*(2*Cd(10310) - MTy2*Cd(10316))*MStop2(1)**2)*
     &      MStop2(3) + 16*Cd(47)*Cd(585)*Cd(9955)*MStop2(1)**3*
     &      UStop2(1,1))/(Cd(47)**2*MStop2(1)**2*MStop2(3))

        Cd(11778) = (8*MTy*MTy2*Cd(47)*
     &      (Cd(11729) + Cd(11722)*MStop2(1)) - 
     &     (MTy2*Cd(11719) + 
     &        8*MTy*(-4*MGl2*MTy2**3*Cd(152) + 
     &           (Cd(47)*(Cd(11712) - MTy2*Cd(11724)) + 
     &              MTy2*Cd(11777))*MStop2(1)))*MStop2(3))/
     &   (MTy*MTy2*Cd(47)**2*MStop2(1)*MStop2(3))

        Cd(11823) = (4*MGl2*MTy*Cd(11786) - 8*MTy2*Cd(11822) + 
     &     (8*Cd(13)*Cd(11782) - 8*MGl2**2*MTy*Cd(152)*Cd(11783) + 
     &        4*MTy*MTy2*Cd(11807) + 
     &        MTy*Cd(47)**2*(MUE*Cd(11781) + MUEC*Cd(11789))*
     &         MStop2(1))*MStop2(3))/
     &   (MTy*Cd(47)**2*MStop2(1)*MStop2(3))

        Cd(12316) = (-4*Cd(12299) - (4*Cd(12300))/MTy2 + 
     &     (MTy*(4*MTy2**2*Cd(164)*Cd(2597) - 
     &           Cd(47)**2*Cd(10190) + Cd(12315)) - 
     &        (4*Cd(12306) + MTy*Cd(12309))*MStop2(1))/
     &      (MTy*MStop2(1)**2))/Cd(47)**2

        Cd(12355) = (-4*Cd(12338) - (4*Cd(12339))/MTy2 + 
     &     (MTy*(4*MTy2**2*Cd(251)*Cd(2681) - 
     &           Cd(50)**2*Cd(10190) + Cd(12354)) + 
     &        (Cd(12345) - MTy*Cd(12348))*MStop2(2))/
     &      (MTy*MStop2(2)**2))/Cd(50)**2

        Cd(12923) = Cd(12701)*Cd(12921) + 
     &   Cd(12757)*Cd(12887)*UCStop(1,3)**2 - 
     &   2*Cd(12888)*MStop2(1)*MStop2(2)*MStop2(3)*
     &    UCStopC(3,1)**2 + 
     &   Cd(9798)*(At**2*Cd(12907) + MUEC*Cd(12920))*UCStopC(1,3)*
     &    UUStop(3,2) + Cd(12700)*
     &    (Cd(12922) - Cd(12700)*Cd(12890)*Cr(91)*UStopC(2,2)*
     &       UUStopC(1,2))

        Cd(127) = -(MTy2*Cd(123)) + MTy*Cd(126) + 
     &   M_3C*Cd(119)*UCStop(1,3) + 
     &   M_3*(3*MTy2**2 + (MGl2 - MStop2(1))**2)*UCStopC(1,3)

        Cd(13945) = MStop2(2)**2*
     &    (-(Cd(13944)*MSbot2(1)) + 
     &      Cr(56)*(MSbot2(1)**2 - At*AtC*MSq2Diff(1,2)) + 
     &      MStop2(2)**2) + 
     &   (Cd(13930) + MUE2*TB2*MStop2(2)**3 + MStop2(2)**4)*
     &    UStop2(1,1)

	end
