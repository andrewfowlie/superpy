* TLhr_self1.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_self1
	implicit none

#include "TLhr.h"


        tadh0 = 1/(128.D0*Pi**4)*
     &   (EL2L*(-(4*Pi*(MTy2*MStop2(1)*MStop2(2)*
     &             (3*Pi*(Cd(47)*Cd(50)*(Cd(365) - Cd(367))) + 
     &               Alfas2L*
     &                (-2*SA*Cd(108) + 
     &                  Cd(40)*Cd(47)*Cd(50)*Cr(9)*
     &                   (-2*Cd(20) + Cd(22)*MStop2(1) + 
     &                     Cd(16)*MStop2(2)))))) + 
     &        (2*Pi*(MTy2*(-4*Alfas2L*MTy2*Cd(358) + 
     &                3*Pi*(Cd(4)*Cd(47)*Cd(50)*Cd(440)*MStop2(1)))
     &              ) + (6*Pi**2*
     &               (MTy2*Cd(10)*Cd(47)*Cd(50)*Cd(442)) + 
     &              (-(MTy2*Cd(14)*Cd(16)*Cd(47)*Cd(50)) + 
     &                 4*Pi*
     &                  (Alfas2L*MTy2*
     &                     (2*CA*
     &                       (-4*MTy*Cd(43)*Cd(47)*Cd(50) + 
     &                       Cd(151)) + 
     &                       Cd(47)*Cd(50)*
     &                       (3*Cd(16)*Cd(18) - 
     &                       Cd(22)*Cd(24)*(-3 + Cr(9)**2)) + 
     &                       2*Cd(355)*Cr(66)) + 
     &                    Pi*
     &                     (MTy*Cd(47)*Cd(50)*
     &                       (3*Cd(37) + 
     &                       32*Pi**2*
     &                       (MTy2*SB*SW*
     &                       (-2*MGl2*Cd(237) + Cd(235)*Cr(65))))))
     &                 )*MStop2(1))*MStop2(2))*MStop2(3) + 
     &        2*MTy2*(-(32*Pi**4*
     &              (A0delGl*MTy*SB*SW*Cd(47)*Cd(50)*Cd(217)*
     &                MStop2(1)*MStop2(2))) - 
     &           4*Pi*(A0delT*Alfas2L*Cd(338)*MStop2(1)*
     &              MStop2(2)) + 
     &           Cd(47)*Cd(50)*
     &            (4*Pi*(Alfas2L*
     &                 (MStop2(1)*MStop2(2)*
     &                    (Cd(24)*(2 + (1 + Cr(9)**2)*Cr(13))*
     &                       MStop2(1) + 
     &                      Cd(18)*(2 + Cr(34))*MStop2(2)) + 
     &                   2*(Cd(18)*MStop2(1) + Cd(24)*MStop2(2))*
     &                    (-Cd(20) + Cd(17)*MStop2(1)*MStop2(2))*
     &                    UStop2(1,1)*UStop2(1,2))) + 
     &              MStop2(1)*MStop2(2)*
     &               (6*Pi**2*
     &                  (CA*
     &                    (Cd(26)*UCStop(1,3) + 
     &                      Cd(27)*UCStopC(1,3))) + 
     &                 Cd(14)*Cr(34)*MStop2(2))))*MStop2(3)))/
     &    (MTy*MW*SB*SW*Cd(47)*Cd(50)*MStop2(1)*MStop2(2)*
     &      MStop2(3))

	DTLHR 'tadh0 =', tadh0 ENDL

        tadHH = 1/(128.D0*Pi**4)*
     &   (EL2L*(4*Pi*(MTy2*MStop2(1)*MStop2(2)*
     &           (3*Pi*(Cd(47)*Cd(50)*(Cd(408) - Cd(410))) - 
     &             Alfas2L*
     &              (2*CA*Cd(108) + 
     &                Cd(47)*Cd(50)*Cd(380)*Cr(9)*
     &                 (-2*Cd(20) + Cd(22)*MStop2(1) + 
     &                   Cd(16)*MStop2(2))))) + 
     &        (4*Pi*(Alfas2L*MTy2*
     &              (-2*MTy2*Cd(401) + 
     &                (Cd(47)*Cd(50)*
     &                    (-8*MTy*SA*Cd(43) + 3*Cd(16)*Cd(371) + 
     &                      Cd(22)*Cd(372)*(-3 + Cr(9)**2)) + 
     &                   2*(SA*Cd(151) + Cd(398)*Cr(66)))*
     &                 MStop2(1)*MStop2(2))) + 
     &           Cd(47)*Cd(50)*
     &            (6*Pi**2*(MTy2*Cd(4)*Cd(9739)*MStop2(1)) + 
     &              (6*Pi**2*(MTy2*Cd(10)*Cd(9741)) - 
     &                 (MTy2*Cd(16)*Cd(370) + 
     &                    12*Pi**2*(MTy*Cd(379)) - 
     &                    128*Pi**4*
     &                     (MTy*MTy2*SB*SW*
     &                       (-2*MGl2*Cd(390) + Cd(388)*Cr(65))))*
     &                  MStop2(1))*MStop2(2)))*MStop2(3) - 
     &        2*MTy2*(32*Pi**4*
     &            (A0delGl*MTy*SB*SW*Cd(47)*Cd(50)*Cd(383)*
     &              MStop2(1)*MStop2(2)) + 
     &           4*Pi*(A0delT*Alfas2L*Cd(394)*MStop2(1)*
     &              MStop2(2)) - 
     &           Cd(47)*Cd(50)*
     &            (4*Pi*(Alfas2L*
     &                 (MStop2(1)*MStop2(2)*
     &                    (-(Cd(372)*(2 + (1 + Cr(9)**2)*Cr(13))*
     &                       MStop2(1)) + 
     &                      Cd(371)*(2 + Cr(34))*MStop2(2)) + 
     &                   2*(Cd(371)*MStop2(1) - Cd(372)*MStop2(2))*
     &                    (-Cd(20) + Cd(17)*MStop2(1)*MStop2(2))*
     &                    UStop2(1,1)*UStop2(1,2))) + 
     &              MStop2(1)*MStop2(2)*
     &               (6*Pi**2*
     &                  (SA*
     &                    (Cd(26)*UCStop(1,3) + 
     &                      Cd(27)*UCStopC(1,3))) + 
     &                 Cd(370)*Cr(34)*MStop2(2))))*MStop2(3)))/
     &    (MTy*MW*SB*SW*Cd(47)*Cd(50)*MStop2(1)*MStop2(2)*
     &      MStop2(3))

	DTLHR 'tadHH =', tadHH ENDL

        selfh0h0 = 1/(64.D0*Pi**2)*
     &   (MTy2*(MStop2(1)*(256*Pi**3*
     &            (Alfa2L*SB2*SW2*Cd(16)*Cd(830)*MStop2(1)*
     &              MStop2(2)**2) + 
     &           Alfa2L*(-2*Alfas2L*Cd(20)*Cd(910) + 
     &              (Alfas2L*
     &                  (2*Cd(20)*Cd(631) + Cd(16)*Cd(1323)) + 
     &                 12*Pi*Cd(733))*MStop2(1) - 
     &              4*Alfas2L*Cd(17)*Cd(682)*MStop2(1)**2) + 
     &           (Alfa2L*(Alfas2L*Cd(22)*Cd(1006) + 
     &                 12*Pi*Cd(749)) + 
     &              (64*Pi**2*
     &                  (SB2*SW2*
     &                    (Cd(925)*dMSfsq1133eps(-1) + 
     &                      Cd(942)*dMSfsq1233eps(-1))) - 
     &                 Alfa2L*
     &                  (Alfas2L*
     &                     (4*Cd(1200) - 4*Cd(1489) + 
     &                       8*(SA2*Cd(1939) + S2A*Cd(2595)) + 
     &                       8*CA2*Cd(9715) + 
     &                       Cr(9)*
     &                       (Cd(16)*Cd(573)*Cr(9) + 
     &                       Cd(22)*(-2*Cd(773) + Cd(476)*Cr(9)))-
     &                        4*Cd(1120)*Cr(13) - 
     &                       4*Cd(1412)*Cr(34)) - 
     &                    6*Pi*
     &                     (4*Cd(916)*dMSfsq1133eps(1) + 
     &                       Cd(931)*dMSfsq1233eps(1))))*MStop2(1))
     &             *MStop2(2)) + 
     &        4*(3*Pi*(Alfa2L*CA*Cd(452)*Conjugate(dAf133eps(0))*
     &              MStop2(1)**2*MStop2(2)*MStop2(3)) + 
     &           6*Pi*(Alfa2L*CA*Cd(422)*Conjugate(dAf133eps(1))*
     &              MStop2(1)**2*MStop2(2)*MStop2(3)) + 
     &           MStop2(3)*
     &            (Alfa2L*MStop2(1)**2*MStop2(2)*
     &               (CA*(3*Pi*(Cd(444)*dAf133eps(0)) + 
     &                    6*Pi*(Cd(416)*dAf133eps(1))) + 
     &                 Alfas2L*
     &                  (2*Cd(756)*Cr(9)*Cr(13)*MStop2(1) + 
     &                    4*Cd(761)*Cr(9)*Cr(34)*MStop2(2))) - 
     &              Alfa2L*
     &               (64*Pi**3*
     &                  (SB2*SW2*Cd(20)*Cd(1244)*MStop2(1)**2) + 
     &                 (64*Pi**3*(SB2*SW2*Cd(20)*Cd(562)) + 
     &                    (-8*Alfas2L*CA2*Cd(946) + 
     &                       Alfas2L*(Cd(497) + Cd(576))*
     &                      Cr(9)**2 + 
     &                       64*Pi**3*
     &                       (SB2*SW2*
     &                       (2*Cd(17)*Cd(821) - 2*MGl2*Cd(4012) + 
     &                       2*MTy2*Cd(9123) + 4*Cd(490)*Cr(13) + 
     &                       4*Cd(593)*Cr(34) - Cd(4015)*Cr(65) + 
     &                       Cd(9125)*Cr(66))))*MStop2(1)**2)*
     &                  MStop2(2))*MStop2(3)) + 
     &           MStop2(1)*
     &            (32*Pi**3*
     &               (A0delGl*Alfa2L*SB2*SW2*Cd(4013)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) - 
     &              64*Pi**3*
     &               (A0delT*Alfa2L*SB2*SW2*Cd(9124)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) - 
     &              16*Pi**2*
     &               (SB2*SW2*Cd(455)*Conjugate(dAf133eps(-1))*
     &                 MStop2(1)*MStop2(2)*MStop2(3)**2) - 
     &              16*Pi**2*
     &               (SB2*SW2*Cd(447)*dAf133eps(-1)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) + 
     &              16*Pi**2*
     &               (SB2*SW2*Cd(462)*dMf133eps(-1)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) - 
     &              16*Pi**2*
     &               (SB2*SW2*Cd(467)*dMf133eps(0)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) + 
     &              6*Pi*(Alfa2L*Cd(439)*dMf133eps(1)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) + 
     &              16*Pi**2*
     &               (SB2*SW2*Cd(927)*dMSfsq1133eps(0)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) + 
     &              16*Pi**2*
     &               (SB2*SW2*Cd(938)*dMSfsq1233eps(0)*MStop2(1)*
     &                 MStop2(2)*MStop2(3)**2) + 
     &              128*Pi**3*
     &               (Alfa2L*SB2*SW2*Cd(17)*Cd(1261)*MStop2(1)**2*
     &                 MStop2(2)*MStop2(3)**2) - 
     &              128*Pi**3*
     &               (Alfa2L*SB2*SW2*Cd(17)*Cd(534)*MStop2(2)**2*
     &                 MStop2(3)**2) - 
     &              Alfa2L*Alfas2L*Cd(688)*MStop2(1)**2*
     &               MStop2(3)**2*UStop2(1,1) - 
     &              2*Alfa2L*Alfas2L*Cd(1208)*MStop2(1)*MStop2(2)*
     &               MStop2(3)**2*UStop2(1,1) - 
     &              4*Alfa2L*Alfas2L*Cd(504)*MStop2(2)**2*
     &               MStop2(3)**2*UStop2(1,1) - 
     &              4*Alfa2L*Alfas2L*Cd(568)*MStop2(1)*MStop2(2)*
     &               (Cd(22)*MStop2(1) + Cd(16)*MStop2(2))*
     &               UStop2(1,1)*UStop2(1,2)))))/
     &    (MW2*SB2*SW2*MStop2(1)**2*MStop2(2)*MStop2(3)**2)

	DTLHR 'selfh0h0 =', selfh0h0 ENDL

        selfHHHH = 1/(64.D0*Pi**2)*
     &   (MTy2*(-2*(Alfa2L*Alfas2L*
     &            (-16*SA2*Cd(946) + 2*Cd(9899) + 
     &              4*S2A*Cd(10069) + 4*CA2*Cd(10100) + 
     &              4*Cd(10915) - 2*Cd(11157) - MTy2*Cd(12696) + 
     &              2*Cd(9773)*Cr(9)**2 + Cd(12699)*Cr(66)) + 
     &           4*Alfa2L*SW2*
     &            (-(Alfas2L*MZ2*SA2*Cd(10188)) + 
     &              32*Pi**3*
     &               (SB2*(2*Cd(17)*Cd(9849) + Cd(16)*Cd(9894) - 
     &                   2*MGl2*Cd(10897) + 
     &                   Cd(20)*
     &                    (Cd(9790) + Cd(9890) + Cd(11690)) + 
     &                   4*Cd(9772)*Cr(13) - 4*Cd(9910)*Cr(34) + 
     &                   Cd(10901)*Cr(65)))) - 
     &           32*Pi**2*(SB2*SW2*
     &              (Cd(9751)*Conjugate(dAf133eps(-1)) - 
     &                Cd(9822)*Conjugate(dY33eps(-1)) + 
     &                Cd(9743)*dAf133eps(-1) + 
     &                Cd(9762)*dMf133eps(-1) - 
     &                Cd(9766)*dMf133eps(0) + 
     &                Cd(9738)*dMf133eps(1) + 
     &                Cd(9928)*dMSfsq1133eps(0) + 
     &                Cd(9941)*dMSfsq1233eps(0) + 
     &                Cd(9828)*dY33eps(-1)))) + 
     &        (Alfa2L*Alfas2L*
     &            (-2*Cd(20)*Cd(9871) - 
     &              Cd(22)*Cr(9)*(2*Cd(9838) + Cd(9768)*Cr(9))*
     &               MStop2(1)*MStop2(2)) + 
     &           2*(-(6*Pi*
     &                 (Alfa2L*Cd(9825)*Conjugate(dY33eps(0))*
     &                   MStop2(1)*MStop2(2))) - 
     &              12*Pi*(Alfa2L*Cd(9806)*Conjugate(dY33eps(1))*
     &                 MStop2(1)*MStop2(2)) + 
     &              Pi*((3*Alfa2L*
     &                    (4*Cd(9916)*dMSfsq1133eps(1) + 
     &                      Cd(9933)*dMSfsq1233eps(1) + 
     &                      2*Cd(9831)*dY33eps(0) + 
     &                      4*Cd(9819)*dY33eps(1)) + 
     &                   32*Pi*
     &                    (SB2*SW2*
     &                      (Cd(9926)*dMSfsq1133eps(-1) + 
     &                       Cd(9946)*dMSfsq1233eps(-1) + 
     &                       4*Pi*
     &                       (Alfa2L*Cd(16)*Cd(9851)*MStop2(2)))))*
     &                 MStop2(1)*MStop2(2)) + 
     &              Alfa2L*
     &               (-(6*Pi*
     &                    (SA*Cd(9754)*Conjugate(dAf133eps(0))*
     &                      MStop2(1)*MStop2(2)*MStop2(3))) - 
     &                 12*Pi*
     &                  (SA*Cd(9727)*Conjugate(dAf133eps(1))*
     &                    MStop2(1)*MStop2(2)*MStop2(3)) + 
     &                 2*MStop2(1)*
     &                  (-(128*Pi**3*
     &                       (SB2*SW2*Cd(17)*Cd(9882)*MStop2(1)*
     &                       MStop2(3))) + 
     &                    MStop2(2)*
     &                     (SA*
     &                       (-(3*Pi*(Cd(9746)*dAf133eps(0))) - 
     &                       6*Pi*(Cd(9721)*dAf133eps(1))) + 
     &                       Alfas2L*
     &                       (-2*Cd(9833)*Cr(9)*Cr(13)*MStop2(1) + 
     &                       4*Cd(9832)*Cr(9)*Cr(34)*MStop2(2))))*
     &                  MStop2(3) - 
     &                 MStop2(2)*MStop2(3)**2*
     &                  (SB2*SW2*
     &                     (128*Pi**3*
     &                       (Cd(22)*Cd(11231)*MStop2(1)**2) - 
     &                       512*Pi**3*
     &                       (Cd(11274)*Cr(13)*MStop2(1)**3) + 
     &                       256*Pi**3*
     &                      (Cd(17)*Cd(9786)*MStop2(2)) + 
     &                       128*Pi**3*
     &                       (Cd(16)*Cd(11423)*MStop2(1)*MStop2(2))
     &                        + 512*Pi**3*
     &                       (Cd(11621)*MStop2(1)*MStop2(2)**2) - 
     &                       512*Pi**3*
     &                       (Cd(11554)*Cr(34)*MStop2(1)*
     &                       MStop2(2)**2)) - 
     &                    A0delGl*Alfas2L*Cd(10898)*MStop2(1) + 
     &                    A0delT*Alfas2L*Cd(12297)*MStop2(1) + 
     &                    Alfas2L*
     &                     (-2*Cd(17)*Cd(11711)*MStop2(1)**2*
     &                       MStop2(2) + 
     &                       8*Cd(9775)*MStop2(2)*UStop2(1,1))) - 
     &                 8*Alfas2L*MStop2(1)*
     &                  (Cd(9791)*MStop2(2)*
     &                     (Cd(22)*MStop2(1) + Cd(16)*MStop2(2)) + 
     &                    Cd(9874)*MStop2(1)*MStop2(3)**2)*
     &                  UStop2(1,1)*UStop2(1,2))))/
     &         (MStop2(1)*MStop2(2)*MStop2(3)**2)))/(MW2*SB2*SW2)

	DTLHR 'selfHHHH =', selfHHHH ENDL

        selfA0A0 = 1/(64.D0*Pi**2)*
     &   (MTy2*(12*Pi*(Alfa2L*MTy*Cd(12857)*
     &           Conjugate(dAf133eps(0))*MStop2(1)*MStop2(2)**2*
     &           MStop2(3)) + 
     &        24*Pi*(Alfa2L*MTy*Cd(12851)*Conjugate(dAf133eps(1))*
     &           MStop2(1)*MStop2(2)**2*MStop2(3)) + 
     &        Alfa2L*(6*Pi*
     &            (MStop2(2)*
     &              (2*Cd(12866)*dMf133eps(-1)*MStop2(3) + 
     &                MStop2(1)*MStop2(2)*
     &                 (4*MTy*Cd(12709)*dMSfsq1133eps(1) + 
     &                   MTy*Cd(12744)*dMSfsq1233eps(1) - 
     &                   2*
     &                    (MTy*
     &                       (Cd(12844)*dAf133eps(0) + 
     &                       2*Cd(12838)*dAf133eps(1)) - 
     &                      Cd(12870)*dMf133eps(0))*MStop2(3)))) + 
     &           Alfas2L*MTy*
     &            (-(Cd(12793)*MStop2(1)) + 
     &              MStop2(2)*
     &               (4*Cd(12944) + 
     &                 MStop2(2)*
     &                  (Cd(12831) + 
     &                    MStop2(1)*
     &                     (Cr(9)*
     &                       (4*Cd(12884) - Cd(12752)*Cr(9) - 
     &                       Cd(12799)*Cr(9)) - 
     &                       4*Cd(12953)*MStop2(2) + 
     &                       32*MStop2(3)*
     &                       (32*Pi**4*
     &                       (-Cd(13224) + Cd(13153)*MStop2(1) + 
     &                       Cd(13437)*MStop2(2) + 
     &                       TB2*
     &                       (Cd(13274) - 2*MGl2*Cd(13320) + 
     &                       2*MTy2*Cd(13366) + Cd(13465) - 
     &                       Cd(13307)*Cr(65) - Cd(13359)*Cr(66))*
     &                       MStop2(3)) - 3*MTy2*MStop2(3))))))) - 
     &        32*MTy*MStop2(1)*MStop2(2)**2*
     &         (32*Pi**4*(A0delGl*Alfa2L*Alfas2L*TB2*Cd(13043)*
     &              MStop2(3)**2) + 
     &           2*A0delT*Alfa2L*Alfas2L*
     &            (3 + 16*Pi**4*(TB2*Cd(13105)))*MStop2(3)**2 + 
     &           (2*Pi**2*(SW2*TB2*
     &                 (-(Cd(12854)*Conjugate(dAf133eps(-1))) + 
     &                   Cd(12725)*Conjugate(dY33eps(-1)) + 
     &                   Cd(12727)*Conjugate(dY33eps(0)) + 
     &                   Cd(12841)*dAf133eps(-1) + 
     &                   Cd(12861)*dMf133eps(1) - 
     &                   Cd(12722)*dMSfsq1133eps(-1) + 
     &                   Cd(12719)*dMSfsq1133eps(0) - 
     &                   Cd(12743)*dMSfsq1233eps(-1) + 
     &                   Cd(12746)*dMSfsq1233eps(0) + 
     &                   Cd(12729)*dY33eps(-1) + 
     &                   Cd(12730)*dY33eps(0))) - 
     &              Alfa2L*Alfas2L*MTy2*Cd(41)*Log(MTy2/Q))*
     &            MStop2(3)**2 + 
     &           Alfa2L*Alfas2L*Cd(12707)*Cd(12708)*
     &            (Cd(22)*MStop2(1) + Cd(16)*MStop2(2))*
     &            UStop2(1,1)*UStop2(1,2))))/
     &    (MTy*MW2*SW2*TB2*MStop2(1)*MStop2(2)**2*MStop2(3)**2)

	DTLHR 'selfA0A0 =', selfA0A0 ENDL

        selfHmHp = 1/(16.D0*Pi**2)*
     &   (MTy2*(Pi*(-(16*Pi*
     &              (Cd(13577)*Conjugate(dAf133eps(0)) + 
     &                Cd(13576)*Conjugate(dAf133eps(1)) + 
     &                Cd(13582)*dAf133eps(0) + 
     &                Cd(13581)*dAf133eps(1) + 
     &                Cd(13597)*dMf133eps(1))) + 
     &           (3*Alfa2L*
     &              (MTy*(Cd(13575)*Conjugate(dAf133eps(-1)) + 
     &                   Cd(13580)*dAf133eps(-1)) + 
     &                Cd(13595)*dMf133eps(-1) + 
     &                2*Cd(13599)*dMf133eps(0)))/
     &            (MTy*SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2))) + 
     &        3*Pi*(Alfa2L*
     &            (2*Cd(13504)*Cd(13505)*Cd(13573)*MSbot2(1)*
     &               MSq2Diff(1,1)**2 + 
     &              2*Cd(13469)*Cd(13470)*Cd(13571)*MSbot2(1)*
     &               MSq2Diff(1,2)**2 - 
     &              Cd(13569)*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2))/
     &          (SW2*TB2*MSbot2(1)*MSq2Diff(1,1)**2*
     &            MSq2Diff(1,2)**2) + 
     &        3*Pi*(Alfa2L*
     &            (Cd(13482) + 
     &              2*Cd(13468)*Cd(13469)*Cd(13470)*MStop2(1)))/
     &          (SW2*TB2*MSq2Diff(1,1)**2*MStop2(1)) - 
     &        3*Pi*(Alfa2L*
     &            (Cd(13534) - 
     &              2*Cd(13504)*Cd(13505)*Cd(13523)*MStop2(2)))/
     &          (SW2*TB2*MSq2Diff(1,2)**2*MStop2(2)) + 
     &        6*Pi*(Alfa2L*
     &            (-(Cd(13557)*MSq2Diff(1,1)) + 
     &              Cd(13549)*MSq2Diff(1,2) + Cd(13539)*MStop2(3)))
     &           /(SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2)*MStop2(3))+
     &          (Alfa2L*Alfas2L*
     &           (-2*Cd(13503)*Cd(13504)*Cd(13505)*MSbot2(1)*
     &              MSq2Diff(1,1)**2 - 
     &             2*Cd(13469)*Cd(13470)*Cd(13501)*MSbot2(1)*
     &              MSq2Diff(1,2)**2 - 
     &             Cd(13499)*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2))/
     &         (SW2*TB2*MSbot2(1)*MSq2Diff(1,1)**2*
     &           MSq2Diff(1,2)**2) + 
     &        (Alfa2L*(Alfas2L*MSq2Diff(1,1)*
     &              (2*Cd(13503)*Cd(13892) - 
     &                Cd(13694)*MSq2Diff(1,1)) - 
     &             (Alfas2L*Cd(13678) + 
     &                (Alfas2L*(8*Cd(43) + TB2*Cd(14356)) + 
     &                   64*Pi**3*
     &                    (SW2*TB2*
     &                      (-2*MGl2*Cd(13877) + 
     &                       2*MTy2*Cd(14115) + Cd(13881)*Cr(65) + 
     &                       Cd(14122)*Cr(66))) - 
     &                   32*Pi**3*
     &                    (SW2*TB2*
     &                      (A0delGl*Cd(13824) + 
     &                       A0delStop(1)*Cd(13887) - 
     &                       A0delT*Cd(14117) + 
     &                       2*Cd(13888)*Cr(13)*MStop2(1))))*
     &                 MSq2Diff(1,1)**2)*MSq2Diff(1,2)**2 + 
     &             2*Alfas2L*Cd(13886)*MSq2Diff(1,2)*MStop2(1)))/
     &         (SW2*TB2*MSq2Diff(1,1)**2*MSq2Diff(1,2)**2) + 
     &        (2*Alfa2L*(MSq2Diff(1,2)*
     &              (16*Pi**3*
     &                 (SW2*TB2*MSq2Diff(1,1)*
     &                   (A0delStop(1)*Cd(13658) + 
     &                     2*(2*Cd(13657) + Cd(13659)*Cr(13))*
     &                      MStop2(1))*MStop2(3)) + 
     &                Alfas2L*Cd(22)*Cd(13623)*Cr(9)*MStop2(1)) + 
     &             Alfas2L*
     &              (Cd(16)*Cd(13631)*MSq2Diff(1,1)*MStop2(2) - 
     &                Cd(13503)*Cd(13891)*Cr(9)*MStop2(3))))/
     &         (SW2*TB2*MSq2Diff(1,1)*MSq2Diff(1,2)*MStop2(3)) + 
     &        (-2*Alfa2L*Alfas2L*
     &            (Cd(13503)*Cd(13504)*Cd(13505) + Cd(13605))*
     &            Cr(9)**2*MStop2(2) + 
     &           4*Alfa2L*Alfas2L*
     &            (2*A0delStop(1)*Cd(13608)*MStop2(2) + 
     &              MStop2(1)*
     &               (2*Cd(13607) + Cd(13610)*Cr(13)*MStop2(2)))*
     &            UStop2(1,1)*UStop2(1,2))/
     &         (SW2*TB2*MSq2Diff(1,2)**2*MStop2(2)) - 
     &        (2*Alfa2L*(16*Pi**3*
     &              (SW2*TB2*MSq2Diff(1,1)**2*
     &                (A0delStop(1)*Cd(13519) + 
     &                  2*(2*Cd(13518) + Cd(13521)*Cr(13))*
     &                   MStop2(1))) + 
     &             Alfas2L*
     &              (-(Cd(13511)*Cr(9)**2) + 
     &                4*Cd(13469)*Cd(13470)*Cd(13503)*UStop2(1,1)*
     &                 UStop2(1,2))))/(SW2*TB2*MSq2Diff(1,1)**2)))/
     &    MW2

	DTLHR 'selfHmHp =', selfHmHp ENDL

	end
