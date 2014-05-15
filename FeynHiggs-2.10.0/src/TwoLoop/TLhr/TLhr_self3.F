* TLhr_self3.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:06

#include "externals.h"
#include "types.h"


	subroutine TLhr_self3
	implicit none

#include "TLhr.h"


        tadA0 = 1/(128.D0*Pi**4)*
     &   (EL2L*MTy2*(16*cI*Pi*
     &         (A0delGl*Alfas2L*Co(2507))/(Cd(47)*Cd(50)) + 
     &        (8*Pi*(MTy2*TB*Cd(47)*Cd(50)*
     &              (32*Pi**3*(MGl2*SW*Co(2511)) + 
     &                cI*(Alfas2L*Co(2515)))) - 
     &           cI*(8*Pi*(Alfas2L*MTy2*
     &                 (2*Co(2509)*Cr(65) - Co(2526)*Cr(66))) + 
     &              Cd(47)*Cd(50)*Cd(12836)*
     &               (MTy2*Cd(11)**2*(Co(2497) - Co(2504)) + 
     &                 4*Pi*
     &                  (Alfas2L*MTy2*(Ce(893) - Ce(896)) + 
     &                    3*Pi*(MTy*Co(2502))))))/
     &         (MTy2*Cd(47)*Cd(50)) + 
     &        2*Pi*(32*Pi**3*(A0delT*SW*TB*Co(2521)) - 
     &           12*cI*Pi*(Cd(12707)*Conjugate(dY33eps(1))) + 
     &           12*cI*Pi*(Cd(12708)*dY33eps(1)) + 
     &           cI*(-8*Alfas2L*MTy2*Co(2524) + 
     &               3*Pi*(Cd(47)*Cd(50)*Cd(12836)*
     &                  (Cd(4)*MStop2(1) - Cd(10)*MStop2(2))))/
     &             (Cd(47)*Cd(50)*MStop2(1)*MStop2(2)) + 
     &           3*cI*Pi*(Cd(25)*Cd(12707)*Conjugate(dY33eps(-1)))/
     &             MStop2(3) + 
     &           6*cI*Pi*(Cd(12707)*Conjugate(dY33eps(0))*Cr(116))/
     &             MStop2(3) - 
     &           cI*(3*Pi*(Cd(12708)*
     &                  (Cd(25)*dY33eps(-1) + 2*Cr(116)*dY33eps(0))
     &                  ) + 
     &               4*Alfas2L*Cd(12836)*Cr(9)*Cr(75)*
     &                (-2*Cd(20) + Cd(22)*MStop2(1) + 
     &                  Cd(16)*MStop2(2)))/MStop2(3) - 
     &           6*cI*Pi*(Cd(26)*UCStop(1,3)) + 
     &           6*cI*Pi*(Cd(27)*UCStopC(1,3)) + 
     &           4*cI*(Alfas2L*Cd(12836)*Co(2501)*MStop2(3)*
     &               UStop2(1,1)*UStop2(1,2))/(MStop2(1)*MStop2(2))
     &           )))/(MTy*MW*SW*TB)

	DTLHR 'tadA0 =', tadA0 ENDL

        selfh0A0 = 1/(64.D0*Pi**2)*
     &   (MTy2*(cI*(Alfa2L*Alfas2L*
     &           (TB*Cd(20)*Co(94)*MStop2(1)**2 + 
     &             2*MStop2(1)*
     &              (Cd(20)*Co(137) + 
     &                TB*(2*Co(66) + Cd(17)*Co(83))*MStop2(1)**2)*
     &              MStop2(2) - 4*TB*Cd(20)*Co(30)*MStop2(2)**2))+
     &          2*MStop2(1)*MStop2(2)*
     &         (-(3*cI*Pi*(Alfa2L*Co(236)*Conjugate(dY33eps(0))*
     &                MStop2(1)*MStop2(2))) + 
     &           3*cI*Pi*(Alfa2L*Co(155)*Conjugate(dAf133eps(0))*
     &              MStop2(1)*MStop2(2)*MStop2(3)) + 
     &           (Pi*((32*Pi*
     &                    (SB*SW2*TB*
     &                      (Co(223)*dMSfsq1133eps(-1) + 
     &                       Co(210)*dMSfsq1233eps(-1))) + 
     &                   3*cI*(Alfa2L*Co(255)*dY33eps(0)))*
     &                 MStop2(1)) + 
     &              4*cI*(Alfa2L*Alfas2L*TB*
     &                 (2*Co(23) - Cd(17)*Co(36))*MStop2(2)))*
     &            MStop2(2) + 
     &           2*(3*cI*Pi*
     &               (Alfa2L*Co(144)*Conjugate(dAf133eps(1))*
     &                 MStop2(1)*MStop2(2)*MStop2(3)) + 
     &              16*Pi**2*
     &               (SB*SW2*TB*Co(152)*Conjugate(dAf133eps(-1))*
     &                 MStop2(1)*MStop2(2)*MStop2(3)**2) + 
     &              16*Pi**2*
     &               (SB*SW2*TB*Co(239)*Conjugate(dY33eps(-1))*
     &                 MStop2(1)*MStop2(2)*MStop2(3)**2) + 
     &              16*Pi**2*
     &               (SB*SW2*TB*
     &                 (Co(179)*dAf133eps(-1) + 
     &                   Co(182)*dAf133eps(0) + 
     &                   Co(162)*dAf133eps(1) + 
     &                   Co(190)*dMf133eps(-1) + 
     &                   Co(225)*dMSfsq1133eps(0) + 
     &                   Co(207)*dMSfsq1233eps(0) + 
     &                   Co(258)*dY33eps(-1))*MStop2(1)*MStop2(2)*
     &                 MStop2(3)**2) + 
     &              Alfa2L*
     &               (4*Cd(16)*
     &                  (16*Pi**3*
     &                     (SB*SW2*TB*
     &                       (Co(50) + Co(118)*MStop2(2)**2)*
     &                       MStop2(3)**2) + 
     &                    cI*
     &                     (Alfas2L*Co(48)*MStop2(2)**2*
     &                       UStop2(1,1)*UStop2(1,2)))*MStop2(1) + 
     &                 (Pi*
     &                     (SW2*TB*
     &                       (-(3*cI*
     &                       (MZ2*
     &                       (Co(193)*dMf133eps(0) - 
     &                       2*Co(172)*dMf133eps(1)))) + 
     &                       64*Pi**2*
     &                       (SB*
     &                       (4*Co(10) + 4*Co(52) + 
     &                       2*Cd(17)*Co(128) + SA*Co(1030) + 
     &                       CA*Co(1902) + 
     &                       4*Cr(13)*
     &                       (Co(18) + Co(103)*MStop2(1)**2) + 
     &                       4*Cr(34)*
     &                       (Co(59) + Co(106)*MStop2(2)**2))))*
     &                       MStop2(1)*MStop2(3)**2) + 
     &                    64*Pi**3*
     &                     (SB*SW2*TB*Cd(22)*
     &                       (Co(6) + Co(114)*MStop2(1)**2)*
     &                       MStop2(3)**2) + 
     &                    4*cI*
     &                     (Alfas2L*Cd(22)*Co(48)*MStop2(1)**2*
     &                       UStop2(1,1)*UStop2(1,2)))*MStop2(2))))
     &        ))/
     &    (MW2*SB*SW2*TB*MStop2(1)**2*MStop2(2)**2*MStop2(3)**2)

	DTLHR 'selfh0A0 =', selfh0A0 ENDL

        selfHHA0 = 1/(64.D0*Pi**2)*
     &   (MTy2*(-(cI*(Alfa2L*Alfas2L*
     &             (TB*Cd(20)*Co(1922)*MStop2(1)**2 + 
     &               2*MStop2(1)*
     &                (Cd(20)*Co(1928) + 
     &                  TB*(2*Co(1918) + Cd(17)*Co(1921))*
     &                   MStop2(1)**2)*MStop2(2) + 
     &               4*TB*Cd(20)*Co(1911)*MStop2(2)**2))) + 
     &        2*MStop2(1)*MStop2(2)*
     &         (-(3*cI*Pi*(Alfa2L*Co(2010)*Conjugate(dY33eps(0))*
     &                MStop2(1)*MStop2(2))) - 
     &           3*cI*Pi*(Alfa2L*Co(1946)*Conjugate(dAf133eps(0))*
     &              MStop2(1)*MStop2(2)*MStop2(3)) + 
     &           (Pi*((32*Pi*
     &                    (SB*SW2*TB*
     &                      (Co(2042)*dMSfsq1133eps(-1) + 
     &                       Co(1999)*dMSfsq1233eps(-1))) - 
     &                   3*cI*(Alfa2L*Co(2028)*dY33eps(0)))*
     &                 MStop2(1)) - 
     &              4*cI*(Alfa2L*Alfas2L*TB*
     &                 (2*Co(1907) + Cd(17)*Co(1913))*MStop2(2)))*
     &            MStop2(2) + 
     &           2*(-(3*cI*Pi*
     &                 (Alfa2L*Co(1935)*Conjugate(dAf133eps(1))*
     &                   MStop2(1)*MStop2(2)*MStop2(3))) + 
     &              16*Pi**2*
     &               (SB*SW2*TB*Co(1943)*Conjugate(dAf133eps(-1))*
     &                 MStop2(1)*MStop2(2)*MStop2(3)**2) + 
     &              16*Pi**2*
     &               (SB*SW2*TB*Co(2013)*Conjugate(dY33eps(-1))*
     &                 MStop2(1)*MStop2(2)*MStop2(3)**2) + 
     &              16*Pi**2*
     &               (SB*SW2*TB*
     &                 (Co(1967)*dAf133eps(-1) + 
     &                   Co(1970)*dAf133eps(0) + 
     &                   Co(1952)*dAf133eps(1) + 
     &                   Co(1978)*dMf133eps(-1) + 
     &                   Co(2044)*dMSfsq1133eps(0) + 
     &                   Co(1996)*dMSfsq1233eps(0) + 
     &                   Co(2031)*dY33eps(-1))*MStop2(1)*MStop2(2)*
     &                 MStop2(3)**2) + 
     &              Alfa2L*
     &               (4*Cd(16)*
     &                  (16*Pi**3*
     &                     (SB*SW2*TB*
     &                       (Co(1915) + Co(1926)*MStop2(2)**2)*
     &                       MStop2(3)**2) + 
     &                    cI*
     &                     (Alfas2L*Co(1914)*MStop2(2)**2*
     &                       UStop2(1,1)*UStop2(1,2)))*MStop2(1) + 
     &                 (Pi*
     &                     (SW2*TB*
     &                       (3*cI*
     &                       (MZ2*
     &                       (Co(1981)*dMf133eps(0) - 
     &                       2*Co(1960)*dMf133eps(1))) + 
     &                       64*Pi**2*
     &                       (SB*
     &                       (4*Co(1904) + 4*Co(1916) + 
     &                       2*Cd(17)*Co(1927) + CA*Co(2142) + 
     &                       SA*Co(2495) + 
     &                       4*Cr(13)*
     &                       (Co(1905) + Co(1923)*MStop2(1)**2) + 
     &                       4*Cr(34)*
     &                       (Co(1917) + Co(1924)*MStop2(2)**2))))*
     &                       MStop2(1)*MStop2(3)**2) + 
     &                    64*Pi**3*
     &                     (SB*SW2*TB*Cd(22)*
     &                       (Co(1903) + Co(1925)*MStop2(1)**2)*
     &                       MStop2(3)**2) + 
     &                    4*cI*
     &                     (Alfas2L*Cd(22)*Co(1914)*MStop2(1)**2*
     &                       UStop2(1,1)*UStop2(1,2)))*MStop2(2))))
     &        ))/
     &    (MW2*SB*SW2*TB*MStop2(1)**2*MStop2(2)**2*MStop2(3)**2)

	DTLHR 'selfHHA0 =', selfHHA0 ENDL

	end
