* TLhr_abbr0.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr0
	implicit none

#include "TLhr.h"


	Cr(1) = MStop2(2) + MStop2(3)

	Cr(101) = (dMTfin*(-MUE + AtC*TB))/TB

	Cr(102) = MUE - AtC*TB

	Cr(108) = Log(Q) + Log(MStop2(1)/Q)

	Cr(11) = 1 - 2*Log(MStop2(2)/Q)*UStop2(1,2)

	Cr(109) = Log(Q) + 2*Re(B0q(MTy2,0.D0,MTy2,1.D0))

	Cr(112) = 6 - 7*UStop2(1,1)

	Cr(115) = UStop2(1,1)**2 + UStop2(1,2)

	Cr(119) = M_3C*MTy + MGlpTmSt2(2)*UCStopC(1,3)

	Cr(122) = UCStop(3,4)*UCStopC(3,1) + UCStop(3,1)*UCStopC(3,3)

	Cr(13) = 1 - Log(MStop2(1)/Q)

	Cr(131) = Re(A0delStop(1)) + Re(A0delStop(2))

	Cr(132) = Re(A0delStop(1)) - Re(A0delStop(2))

	Cr(138) = M_3C*UCStop(3,3) + M_3*UCStop(3,4)

	Cr(14) = 1 - Log(MGl2/Q)

	Cr(15) = 1 - Log(MTy2/Q)

	Cr(16) = Log(Q) + Re(B0q(MStop2(1),MGl2,MTy2,1.D0))

	Cr(2) = MStop2(2) + 3*MStop2(3)

	Cr(20) = Re(A0delGl) + Re(A0delT)

	Cr(21) = Cr(20) - 1/2.D0*Re(A0delStop(1))

	Cr(22) = Re(A0delStop(1)) - 2*Re(A0delStop(2))

	Cr(23) = 3/2.D0*Re(A0delStop(1)) - 2*Re(A0delStop(2))

	Cr(24) = Cr(22) - Cr(23)*UStop2(1,1)

        Cr(25) = -(1/2.D0*(Re(A0delStop(1))*UStop2(1,2)**2)) + 
     &   Cr(24)*UStop2(1,1)

	Cr(26) = MStop2(2) - 2*MStop2(3)

	Cr(27) = 2*(MStop2(2) - 2*MStop2(3))

	Cr(28) = UStop2(1,1)**2 + UStop2(1,2)**2

	Cr(30) = MGlpT2 + MGlpTmSt2(2)

	Cr(32) = 1 - 2*Log(MStop2(1)/Q)

	Cr(34) = 1 - Log(MStop2(2)/Q)

	Cr(36) = Log(Q) + Re(B0q(MStop2(2),MGl2,MTy2,1.D0))

	Cr(4) = MGlpT2 + MGlpTmSt2(1)

	Cr(40) = Cr(20) - 1/2.D0*Re(A0delStop(2))

	Cr(41) = -(1/2.D0*Re(A0delStop(2))) + Re(A0delStop(1))

        Cr(42) = 1/2.D0*(Cr(28)*Re(A0delStop(2))) + 
     &   2*Cr(41)*UStop2(1,1)*UStop2(1,2)

	Cr(47) = 2 + 4*UStop2(1,2)

	Cr(48) = Cr(47)*UStop(1,1) - UCStop(3,1)*UStop(2,1)

	Cr(49) = 2 + UStop2(1,1)

	Cr(5) = M_3C*UCStop(1,3) + M_3*UCStopC(1,3)

	Cr(50) = 2*Cr(49)*UStop(1,2) - UCStop(3,1)*UStop(2,2)

	Cr(53) = 1 + UStop2(1,1)**2

	Cr(54) = 2 + UStop2(1,2)

	Cr(55) = 2*Cr(54)*UStop(1,1) - UCStop(3,1)*UStop(2,1)

	Cr(56) = 1 + UStop2(1,1)

	Cr(57) = 2*Cr(56)*UStop(1,2) - UCStop(3,1)*UStop(2,2)

	Cr(58) = UUStop(3,3) - UUStop(3,4)

	Cr(60) = Cr(53) - Cr(58)*UUStopC(3,4)

	Cr(62) = MGlpTmSt2(1) + 2*MStop2(1)

	Cr(63) = MGlpTmSt2(2) + 2*MStop2(2)

	Cr(65) = MGl2 - MGl2*Log(MGl2/Q)

	Cr(66) = MTy2 - MTy2*Log(MTy2/Q)

	Cr(68) = Log(Q) + Re(B0q(MTy2,MGl2,MStop2(1),1.D0))

	Cr(69) = Log(Q) + Re(B0q(MTy2,MGl2,MStop2(2),1.D0))

	Cr(7) = 2 - 3*UStop2(1,1)

	Cr(70) = Cr(65) + Cr(66)

	Cr(72) = Cr(16)*MGlpTmSt2(1) + 2*Log(Q)*MStop2(1)

	Cr(75) = UStop2(1,1) - UStop2(1,2)

	Cr(8) = Cr(7)*UStop2(1,1) - UStop2(1,2)**2

	Cr(84) = MStop2(2) - 7*MStop2(3)

	Cr(85) = 1 + UStop2(1,2)

	Cr(86) = 4*MGlpT2 - MStop2(1)

	Cr(87) = UStop2(1,1)**2 - UStop2(1,2)

	Cr(88) = 1 - 3*UStop2(1,1)*UStop2(1,2)

	Cr(9) = 1 - 2*UStop2(1,1)

	Cr(91) = UUStop(3,3) + UUStop(3,4)

	Cr(92) = 6 + UStop2(1,1)

	Cr(93) = Cr(91)*UCStop(3,1) - Cr(92)*UUStop(1,3)

	Cr(94) = Cr(93)*MStop2(3) - 8*M_3*MTy*UUStop(1,2)

	Cr(95) = AtC + 2*M_3C*UStop2(1,2)**2

	Cr(96) = 2*(MUE - TB*Cr(95))

	Cr(10) = -Cr(75) + Cr(8)*Log(MStop2(1)/Q)

	Cr(103) = Cr(14) - Log(Q)

	Cr(110) = Cr(69) + Log(MStop2(2)/Q)

	Cr(113) = 1 - Cr(112)*UStop2(1,1) + UStop2(1,2)**2

        Cr(114) = MGlpTmSt2(1) - 
     &   2*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))

	Cr(117) = Cr(114) - MGlpTmSt2(1) - MStop2(3)

        Cr(12) = Cr(10)*MStop2(1) + 
     &   2*(Cr(11)*MStop2(2) + Cr(9)*MStop2(3))*UStop2(1,1)

        Cr(139) = Log(MStop2(1)/Q)*MStop2(1) - 
     &   Log(MStop2(2)/Q)*MStop2(2) + MStop2(3)

	Cr(17) = MGl2*Cr(14) + MTy2*Cr(15) + Cr(16)*MGlpTmSt2(1)

        Cr(18) = Cr(13) - 4*Log(Q) - 
     &   4*Re(B0q(MStop2(1),0.D0,MStop2(1),1.D0))

	Cr(19) = -4*MTy*Cr(5)*Cr(16) + 2*Cr(17) - Cr(18)*MStop2(1)

        Cr(29) = 2*Cr(28)*MStop2(2) + 4*Cr(26)*UStop2(1,1) - 
     &   2*Cr(27)*UStop2(1,1)**2

        Cr(3) = 4*Cr(1)*UStop2(1,1) - 2*Cr(2)*UStop2(1,1)**2 + 
     &   2*MStop2(1)*UStop2(1,2)**2

	Cr(31) = 4*MTy*Cr(5) + 2*Cr(30) + 3*MStop2(2)

        Cr(33) = Cr(32)*MStop2(1) + Log(MStop2(2)/Q)*MStop2(2) - 
     &   MStop2(3)

	Cr(37) = MGl2*Cr(14) + MTy2*Cr(15) + Cr(36)*MGlpTmSt2(2)

        Cr(38) = Cr(34) - 4*Log(Q) - 
     &   4*Re(B0q(MStop2(2),0.D0,MStop2(2),1.D0))

	Cr(39) = 4*MTy*Cr(5)*Cr(36) + 2*Cr(37) - Cr(38)*MStop2(2)

        Cr(51) = 3*UStop2(1,2)**2 + Cr(48)*UStopC(1,1) - 
     &   Cr(50)*UStopC(1,2)

	Cr(52) = -(Cr(51)*UStop2(1,1)) + 3*UStop2(1,1)**3 - UStop2(1,2)

        Cr(59) = 3*UStop2(1,2)**2 + Cr(55)*UStopC(1,1) - 
     &   Cr(57)*UStopC(1,2)

	Cr(6) = 4*Cr(4) - 8*MTy*Cr(5) + 6*MStop2(1)

        Cr(67) = 2*Log(Q) + Re(B0q(MStop2(1),MGl2,MTy2,1.D0)) + 
     &   Re(B0q(MStop2(2),MGl2,MTy2,1.D0))

        Cr(73) = Cr(36)*MGlpTmSt2(2) + 
     &   2*MStop2(2)*(Log(Q) + Re(B0q(MStop2(2),0.D0,MStop2(2),1.D0)))

	Cr(80) = Cr(34) - Log(Q)

        Cr(81) = MGlpTmSt2(2) + 
     &   2*MTy*(M_3C*UCStop(1,3) + M_3*UCStopC(1,3))

        Cr(61) = 1 - Cr(59)*UStop2(1,1) + 4*UStop2(1,1)**3 - 
     &   Cr(60)*UStop2(1,2)

	Cr(83) = -Cr(131) + 2*Re(A0delGl) + 2*Re(A0delT)

        Cr(89) = Cr(87)*MStop2(2) - Cr(88)*MStop2(3) + 
     &   8*M_3*MTy*UCStopC(1,3) + Cr(86)*UStop2(1,2)

        Cr(97) = 4*MGl2 + 4*MTy2 - Cr(85)*MStop2(1) - 
     &   MStop2(2)*UStop2(1,1)

	Cr(104) = MTy2*Cr(15) + MGl2*Cr(103)

	Cr(105) = MGlpT2 + 2*MTy*Cr(5)

	Cr(107) = 2*MTy*Cr(5) + MStop2(1)

	Cr(116) = Cr(13)*MStop2(1) - Cr(34)*MStop2(2)

	Cr(118) = Cr(113)*Cr(116) - 2*Cr(117)*Log(Q)

	Cr(120) = Cr(114) + 2*M_3*MTy*UCStopC(1,3)

	Cr(121) = M_3C*MTy*Cr(28) - Cr(120)*UCStopC(1,3)

	Cr(133) = 2*Cr(75)*Cr(122) + Cr(113)*UCStopC(1,3)

	Cr(134) = -2*Cr(20) + Cr(131)

	Cr(135) = (2*MTy*TB*Cr(132)*Cr(133) + Cr(102)*Cr(134))/(MTy*TB)

        Cr(136) = -2*Re(PhiAt*Cr(20)*Cr(102)) + 
     &   Re(PhiAt*Cr(102)*Cr(131))

        Cr(137) = 2*Pi*(Alfas2L*
     &     (Cr(136)/(MTy*TB) + 2*Re(PhiAt*Cr(132)*Cr(133))))

        Cr(35) = Cr(28)*Cr(34)*MStop2(2) + 
     &   2*Cr(33)*UStop2(1,1)*UStop2(1,2)

	Cr(43) = Cr(9) + 2*UStop2(1,1)

	Cr(45) = -2*MTy*Cr(5) + MGlpTmSt2(1)

	Cr(46) = 2*MTy*Cr(5) + MGlpTmSt2(2)

	Cr(74) = Cr(67) - 2*Log(Q)

        Cr(90) = -(Cr(89)*UCStopC(1,3)) + 
     &   4*M_3C*MTy*Cr(85)*UStop2(1,1)

	Cr(98) = Cr(84)*UStop2(1,1)**2 + Cr(97)*UStop2(1,2)

	Cr(99) = Cr(90) + Cr(98)*UCStopC(1,3) - Cr(94)*UUStopC(1,1)

        Cr(106) = 2*Cr(104) - 
     &   MGlpT2*Re(B0q(MTy2,MGl2,MStop2(1),1.D0)) - 
     &   Cr(105)*Re(B0q(MTy2,MGl2,MStop2(2),1.D0))

	Cr(44) = 1 - UStop2(1,1) - Cr(43)*UStop2(1,2)

        Cr(76) = Cr(74)*Cr(75) - 
     &   2*Re(B0q(MStop2(1),MGl2,MTy2,1.D0))*UStop2(1,1) + 
     &   2*Re(B0q(MStop2(2),MGl2,MTy2,1.D0))*UStop2(1,2)

        Cr(79) = 2*Cr(20)*Cr(44) + Cr(52)*Re(A0delStop(1)) - 
     &   Cr(61)*Re(A0delStop(2))

        Cr(123) = Cr(75)*Cr(116)*Cr(122) - 
     &   2*UCStopC(1,3)*(MStop2(1)*
     &       Re(B0q(MStop2(1),0.D0,MStop2(1),1.D0)) - 
     &      MStop2(2)*Re(B0q(MStop2(2),0.D0,MStop2(2),1.D0)) + 
     &      M_3*MTy*Log(Q)*UCStopC(1,3))

        Cr(111) = 2*MTy2*Cr(109) + Cr(108)*MStop2(1) + 
     &   Cr(110)*MStop2(2) - MStop2(4) + 
     &   Cr(107)*Re(B0q(MTy2,MGl2,MStop2(1),1.D0))

        Cr(124) = 2*(Cr(123) + 2*M_3C*MTy*Cr(115)*Log(Q) + 
     &      Cr(121)*Re(B0q(MStop2(1),MGl2,MTy2,1.D0)) + 
     &      Cr(119)*Re(B0q(MStop2(2),MGl2,MTy2,1.D0))) + 
     &   Cr(118)*UCStopC(1,3)

        Cr(71) = -(Cr(45)*Cr(68)) - Cr(46)*Cr(69) + 
     &   2*Cr(44)*Cr(70) + Cr(13)*Cr(52)*MStop2(1) - 
     &   Cr(34)*Cr(61)*MStop2(2)

        Cr(77) = -(MTy*Cr(5)*Cr(76)) + 
     &   2*MTy2*Re(B0q(MTy2,0.D0,MTy2,1.D0)) - Cr(72)*UStop2(1,1) - 
     &   2*MStop2(1)*Re(B0q(MStop2(1),0.D0,MStop2(1),1.D0))*
     &    UStop2(1,1) - Cr(73)*UStop2(1,2)

        Cr(64) = -4*MTy2 - 2*MGlpT2*Cr(44) + Cr(45) + Cr(46) - 
     &   Cr(52)*MStop2(1) + Cr(61)*MStop2(2) + 
     &   2*Cr(62)*UStop2(1,1) + 2*Cr(63)*UStop2(1,2)

        Cr(82) = -2*MGl2*(-1 + Log(MGl2/Q) + Log(Q)) + 
     &   (-1 + Log(Q) + Log(MStop2(1)/Q))*MStop2(1) - 
     &   Cr(80)*MStop2(2) + 
     &   2*MTy2*(1 - Log(MTy2/Q) + Log(Q) + 
     &      2*Re(B0q(MTy2,0.D0,MTy2,1.D0))) - 
     &   Cr(45)*Re(B0q(MTy2,MGl2,MStop2(1),1.D0)) - 
     &   Cr(81)*Re(B0q(MTy2,MGl2,MStop2(2),1.D0))

	Cr(100) = MTy*Cr(96) - TB*Cr(99)

	Cr(125) = Cr(106) + Cr(111)

	Cr(126) = (2*MTy*TB*Cr(124) - Cr(102)*Cr(125))/(MTy*TB)

	Cr(127) = Im(PhiAt*Cr(101) - 1/(6.D0*Pi)*(Alfas2L*PhiAt*Cr(126)))

	Cr(128) = Re(PhiAt*Cr(102)*Cr(106)) + Re(PhiAt*Cr(102)*Cr(111))

        Cr(129) = 4*Pi*(Alfas2L*
     &     (Cr(128)/(MTy*TB) - 2*Re(PhiAt*Cr(124))))

	Cr(130) = 1/(24.D0*Pi**2)*Cr(129) + Re(PhiAt*Cr(101))

        Cr(78) = Alfas2L*(Cr(71) + 2*Cr(77)) - 
     &   6*Pi*(dMTfin*MTy) + 4*Alfas2L*MTy2*Log(Q)

	end
