* gm2_2L_2.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:39
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_2
	implicit none

#include "gm2_2L.h"

	ComplexType tmp1, tmp2, tmp3, tmp4

	tmp1 = ML2**2*(9/40.D0 - 33/40.D0*ML2/MZ2)

	tmp2 = SA2/Ss44 + CA2/Ss49

	tmp3 = -3 + Ss324 - 6*Ss474 + Ss306*(-6 + 12*Ss474) + Ss484

	tmp4 = 8*Ss45*Ss46*Ss47*Ss48*Ss628 + 
     &   Ss12*(Ss12*Ss63 + 4*Ss482*Ss63)

	gm2_2L = (1/2880.D0*(ML2*(-1657*ML2 + 1420*MW2 - 931*Ss31) + 
     &         MW2*Ss31) + 
     &      1/1440.D0*(-49*ML2**4*Ss31 - 3*Ss430*Ss431)/MW2**3)/
     &    (Ss28**2*SW2) + (ML2**2*Ss325*
     &      (MW2*(1/240.D0*Ss324 + 1/20.D0*(Ss306*Ss398) + 
     &           1/240.D0*Ss455) - 1/20.D0*Ss534))/(MW2**3*Ss28*SW2) + 
     &   (1/1440.D0*(49*ML2**3 + 3*ML2*Ss320*Ss321 + 
     &         18*ML2**2*Ss325*Ss327 + 
     &         3*(MBot2 - MTop2)**2*(49*Ss29 + 228*Ss371) + 
     &         216*Ss508*Ss515) - 1/960.D0*Ss615/Ss37**3)/
     &    (MW2**3*SW2) + 1/90.D0*
     &    (MBot2*(23 - 36*Ss141**2 + 60*Ss474 + 6*Ss484 + 
     &         (43*Ss10 + MBot2*(46 + 120*Ss141 + 12*Ss214) + 
     &            72*Ss141*Ss475 + 12*Ss479 - 72*Ss551 + 18*Ss719)/
     &          MZ2))/Ss7 + 
     &   1/30.D0*(ML2*(23 - 36*Ss306**2 + 60*Ss474 + 6*Ss484 + 
     &         (43*Ss11 + ML2*(46 + 120*Ss306 + 12*Ss324) + 
     &            72*Ss306*Ss476 + 12*Ss481 - 72*Ss625 - 36*Ss800)/
     &          MZ2))/Ss8 + 
     &   (MTop2*(2/45.D0*(-23 + 36*Ss335**2 - 60*Ss474 - 6*Ss484) - 
     &        2/45.D0*(43*Ss12 + MTop2*(46 + 120*Ss335 + 12*Ss375) + 
     &            72*Ss335*Ss477 + 12*Ss482 - 72*Ss628 + 36*Ss804)/
     &          MZ2))/Ss9 + 
     &   1/60.D0*((5 - 4*Ss885)*
     &       ((CA*Ss884)/Ss878 + (SA*Ss897)/Ss892))/(CB2*CW) + 
     &   1/60.D0*((5 - 4*Ss906)*
     &       ((CA*Ss905)/Ss903 + (SA*Ss915)/Ss913))/(CB2*CW) + 
     &   3/20.D0*((5 - 4*Ss928)*
     &       ((CA*Ss927)/Ss921 + (SA*Ss940)/Ss935))/(CB2*CW) + 
     &   3/20.D0*((5 - 4*Ss949)*
     &       ((CA*Ss948)/Ss946 + (SA*Ss958)/Ss956))/(CB2*CW) + 
     &   1/15.D0*((CA*Ss968*(5 - 4*Ss969))/Ss964 + 
     &       SA*((-5*Ss981)/Ss976 + (4*Ss969*Ss981)/Ss976))/
     &     (CB*CW*SB) + 1/15.D0*
     &    ((CA*Ss987*(5 - 4*Ss988))/Ss986 + 
     &       SA*((-5*Ss997)/Ss995 + (4*Ss988*Ss997)/Ss995))/
     &     (CB*CW*SB) + 1/5760.D0*
     &    ((ML2**2*(1954*ML2 + 537*Ss31))/MW2 + 
     &       (ML2*(12*Ss306*Ss308 + ML2**2*(248*ML2 + 79*Ss31) + 
     &            12*Ss308*Ss398))/MW2**2)/(Ss28**2*SW2) + 
     &   1/720.D0*(MBot2*(-6*MBot2*Ss211*Ss46*Ss47*
     &          (Ss141*(6*Ss141 - 12*Ss474) - Ss484) + 
     &         23*MBot2*Ss56 + 
     &         6*(5*MZ2 + Ss10 + 2*Ss474*Ss478)*Ss56 + 
     &         (12*MBot2*Ss211*Ss46*Ss47*
     &             (MBot2*(6*Ss141**2 + Ss214) - 6*Ss551) - 
     &            (MBot2*(74*MBot2 - 31*Ss10) + Ss141*Ss142 - 
     &               12*Ss479*Ss480)*Ss56)/MZ2))/(CW2*Ss7**2*SW2)+
     &     1/5760.D0*Ss780/(MW2**3*Ss28*Ss37**3*SW2)

	gm2_2L = gm2_2L + (CA*MTop2**2*SA*Ss63*
     &      (11/20.D0*MTop2/(MZ2*Ss17) + 3/20.D0/Ss18))/
     &    (CB*CW2*SB*Ss9*SW2) + 
     &   (CA*MTop2*SA*Ss45*Ss46*Ss47*Ss48*
     &      (-3 + Ss375 - 6*Ss474 + Ss335*(-6 + 12*Ss474) + Ss484)*
     &      (1/120.D0 - 1/120.D0*Ss44/Ss49))/(CB*CW2*SB*Ss44*Ss9*SW2)+
     &     (Ss78**2*TB2*(-(1/480.D0*
     &           (ML2*(43*Ss11 + 
     &                6*MZ2*
     &                 (Ss324 + 10*Ss474 + Ss306*(10 + 12*Ss474))+
     &                  2*Ss627))/Ss62) + 
     &        1/240.D0*(ML2**2*
     &            (23 + Ss306*(60 + 36*Ss306) + 6*Ss324))/Ss68))/
     &    (CW2*Ss8*SW2) + (1/480.D0*Ss877/(MM2*MZ2) + 
     &      1/360.D0*(MTop2*(36*MTop2*Ss335**2*Ss46*Ss47*Ss48 - 
     &            (MTop2*(23 + 60*Ss474) + 
     &               6*(Ss12 + MZ2*(5 + 6*Ss335 + 6*Ss474) - 
     &                  12*Ss482))*Ss63))/Ss9**2)/(CW2*SW2) + 
     &   1/240.D0*(Ss53**2*(ML2*
     &          (6*(Ss11 + MZ2*(5 + 6*Ss306 + 6*Ss474) - 
     &               12*Ss481) + 
     &            ML2*(23 - 36*Ss306**2 + 60*Ss474 + 6*Ss484)) + 
     &         (ML2**2*(31*Ss11 + 
     &              ML2*(-74 - 168*Ss306 + 12*Ss324) + 
     &              72*Ss306*Ss476 + 156*Ss481 - 72*Ss625))/MZ2))/
     &     (CW2*MM2*Ss8**2*SW2) + 
     &   3/160.D0*(ML2*(Ss51**2 + (Ss49*Ss55**2)/Ss44)*
     &       (Ss11*(Ss11 + 4*Ss481) + 8*Ss52*Ss625))/
     &     (CB2*CW2*MZ2*Ss49*Ss8*SW2) + 
     &   1/360.D0*(MTop2**2*((-6*Ss46*Ss47*Ss48*
     &             (12*Ss335*Ss477 + Ss629) + 
     &            (-31*Ss12 + MTop2*(74 + 168*Ss335) - 156*Ss482)*
     &             Ss63)/(MZ2*Ss9**2) + 
     &         ((23 + Ss335*(60 + 36*Ss335) + 6*Ss375)*Ss63)/
     &          (Ss5*Ss9)))/(CW2*SW2) + 
     &   (1/240.D0*(MTop2*(600*MBot2*Ss432 + 
     &            Ss32*Ss33*(121*Ss34 + 228*Ss432) + MBot2*Ss34*Ss80)
     &           + MBot2*(MW2*(-137*Ss34 - 132*Ss432) + 
     &            372*MTop2*Ss432 + MBot2*(188*Ss34 + 240*Ss432))*TB2)+
     &        1/80.D0*(-2*MTop2**2*Ss32*Ss33*(13*Ss34 + 36*Ss432) + 
     &          MBot2*MTop2**2*(Ss34*Ss79 + Ss432*(36 - 144*TB2)) - 
     &          MBot2**3*(17*Ss34 + 36*Ss432)*TB2 + 
     &          MBot2**2*MTop2*
     &           (Ss34*(1 + 61*TB2) + Ss432*(36 + 180*TB2)))/MW2)/
     &    (Ss25*Ss37*SW2) + 
     &   1/120.D0*((MBot2**2*Ss108 - 6*Ss398*Ss423 + 
     &          MTop2*(3*Ss375*Ss382 + 
     &             Ss335*(-6*Ss360 + 36*Ss382*Ss398) - 
     &             Ss32*Ss33*Ss43 + 3*Ss382*Ss455) + 
     &          MBot2*(3*Ss214*Ss262 + 
     &             Ss141*(-6*Ss175 + 36*Ss262*Ss398) + 
     &             3*Ss262*Ss455 + Ss88) + 27*MBot2**4*Ss72*TB2 - 
     &          6*MBot2**3*Ss77*TB2)/(MW2*Ss37) + 
     &       (-6*Ss275*Ss293 + 
     &          MTop2*(-6*Ss335*Ss349 - Ss26*Ss27*Ss40) + 
     &          MBot2*(6*Ss141*Ss159 + Ss94 + 
     &             MBot2*(-Ss102 - 3*MBot2*Ss74*TB2)))/Ss24)/
     &     (Ss25**3*SW2)

	gm2_2L = gm2_2L + (Ss51**2*tmp1)/(CB2*CW2*Ss3*Ss8*SW2) + 
     &   (Ss55**2*tmp1)/(CB2*CW2*Ss15*Ss8*SW2) - 
     &   1/240.D0*(MBot2*(Ss211*
     &          (Ss214*Ss46*Ss47 + 
     &            Ss46*Ss47*(12*Ss141*Ss474 + Ss484)) - 3*Ss56 - 
     &         6*Ss141*Ss56 - 6*Ss474*Ss56)*Ss58*tmp2)/
     &     (CB2*CW2*Ss7*SW2) + 
     &   1/160.D0*(MBot2*(Ss10*(Ss10*Ss56 + 4*Ss479*Ss56) + 
     &         8*Ss211*Ss46*Ss47*Ss551*Ss58)*tmp2)/
     &     (CB2*CW2*MZ2*Ss7*SW2) + 
     &   (1/240.D0*(3*Ss449*Log(MBot2) - 2*Ss453*Log(ML2))/CB + 
     &      1/80.D0*(Ss441*Log(MTop2)))/(MW2**3*MZ2*Ss25**3*SW2**2) + 
     &   1/40.D0*Log(MZ2)/SW2**2 + 
     &   1/720.D0*(12*Ss128*Ss632 + 18*Ss503*Ss648 + Ss708 + 
     &       6*Ss498*Sqrt(Ss496)*Log(Abs(Ss500/Ss499)) + 
     &       18*Ss503*Sqrt(Ss501)*Log(Abs(Ss505/Ss504)))/
     &     (MZ2**3*SW2**2) + 
     &   (1/240.D0*(-2*Ss636*Ss639 - 6*Ss642*Ss647 - 
     &         Ss485*Log(Ss28**2)) + 
     &      1/80.D0*(MBot*MTop*Ss489*Ss495*
     &          Log(1/4.D0*Sq(Ss488)/(MBot2*MTop2)))/
     &        (Ss490*Ss491))/(MW2**3*MZ2*SW2**2) + 
     &   1/30.D0*(MBot2*SA2*(33 + 12*(-2 + Ss141)*Ss141 + 12*Ss263 - 
     &         2*Ss268 + (-3*Ss14 - 24*Ss141*Ss265 - 12*Ss266 + 
     &            MBot2*(-156 + 72*Ss141 - 4*Ss214 - 24*Ss263 + 
     &               4*Ss268) + 24*Ss547)/MHiggs2(1)))/(CB2*Ss13)+
     &     (MBot2**2*SA2*(-(3/40.D0*Ss56/(CW2*Ss7*SW2)) + 
     &        1/15.D0*(-3*Ss14 + MBot2*(66 + 72*Ss141 + 4*Ss214) + 
     &            24*Ss141*Ss264 - 12*Ss266 - 24*Ss547)/
     &          MHiggs2(1)**2))/(CB2*Ss13) + 
     &   (CA*MTop2**2*SA*(-(3/20.D0*Ss63/(CW2*Ss9*SW2)) + 
     &        4/15.D0*(-3*Ss19 + 24*Ss335*Ss362 - 12*Ss369 + 
     &            MTop2*(66 + 72*Ss335 + 4*Ss375) - 24*Ss619)/
     &          MHiggs2(1)**2))/(CB*SB*Ss17) + 
     &   (MBot2*(1/240.D0*(SA2*
     &            (Ss211*(Ss214*Ss46*Ss47 + 
     &                 (12*Ss141*Ss263 + Ss268)*Ss46*Ss47) - 
     &              3*Ss56 - 6*Ss141*Ss56 - 6*Ss263*Ss56)*Ss60)/
     &          (CB2*CW2*Ss44*SW2) - 2/5.D0*Ss715/MHiggs2(1)**2))/
     &    Ss13 + (ML2*(1/80.D0*
     &         ((-3 - 6*Ss306 + Ss263*(-6 + 12*Ss306) + Ss324)*
     &            Ss54*Ss55**2)/(CB2*CW2*Ss44*SW2) + 
     &        6/5.D0*Ss789/MHiggs2(1)**2))/Ss15 + 
     &   (MTop2*(-(1/120.D0*(CA*SA*
     &              (Ss268*Ss46*Ss47*Ss48 + 
     &                (12*Ss263*Ss335 + Ss375)*Ss46*Ss47*Ss48 - 
     &                3*Ss63 - 6*Ss263*Ss63 - 6*Ss335*Ss63)*Ss64)/
     &            (CB*CW2*SB*Ss44*SW2)) + 8/5.D0*Ss798/MHiggs2(1)**2)
     &      )/Ss17 + (Ss55**2*
     &      (-(1/80.D0*(ML2*Ss52*tmp3)/Ss44) - 
     &        1/20.D0*(ML2**3*(6*Ss306**2 + Ss324))/(Ss15*MHiggs2(1))
     &        ))/(CB2*CW2*Ss8*SW2)

	gm2_2L = gm2_2L + (1/240.D0*
     &       (CA2*MBot2*(Ss211*
     &             (Ss214*Ss46*Ss47 + 
     &               (12*Ss141*Ss269 + Ss274)*Ss46*Ss47) - 
     &            3*Ss56 - 6*Ss141*Ss56 - 6*Ss269*Ss56)*Ss57)/
     &        (Ss1*Ss49) + 
     &      33/40.D0*(ML2**3*Ss55**2)/(MZ2*Ss15*MHiggs2(1)))/
     &    (CB2*CW2*SW2) + (ML2*SA2*
     &      (1/10.D0*(-12*Ss263 + 2*Ss268 - 
     &           3*(11 + Ss306*(-8 + 4*Ss306))) + 
     &        1/10.D0*(3*Ss16 + 24*Ss306*Ss315 + 12*Ss318 + 
     &            ML2*(156 - 72*Ss306 + 4*Ss324) - 24*Ss617)/
     &          MHiggs2(1)))/(CB2*Ss15) + 
     &   (CA*MTop2*SA*(2/15.D0*
     &         (12*Ss263 - 2*Ss268 + 3*(11 + Ss335*(-8 + 4*Ss335)))
     &          - 2/15.D0*(3*Ss19 + 24*Ss335*Ss364 + 12*Ss369 + 
     &            MTop2*(156 + 24*Ss263 - 4*Ss268 - 72*Ss335 + 
     &               4*Ss375) - 24*Ss619)/MHiggs2(1)))/(CB*SB*Ss17)
     &     + (CA*SA*(-(11/20.D0*(MTop2**3*Ss63)/Ss18) + 
     &        1/30.D0*(MTop2**2*
     &            (MTop2*Ss366*(6*Ss335**2 + Ss375)*Ss46*Ss47*Ss48 - 
     &              6*Ss335*Ss338*Ss63))/(Ss17*MHiggs2(1))))/
     &    (CB*CW2*MZ2*SB*Ss9*SW2) - 
     &   3/5.D0*(ML2**2*((SA2*
     &            (-(Ss16/MHiggs2(1)**2) + 
     &              (8*Ss306*Ss313)/MHiggs2(1)**2 - 
     &              (4*Ss618)/MHiggs2(1)**2))/Ss15 + 
     &         (22*ML2*SA2)/(Ss15*MHiggs2(1)**2) + 
     &         CA2*(((8*Ss306*Ss314)/MHiggs2(2)**2 - 
     &               Ss4/MHiggs2(2)**2 - (4*Ss622)/MHiggs2(2)**2)/
     &             Ss3 + (22*ML2)/(Ss3*MHiggs2(2)**2))))/CB2 - 
     &   1/60.D0*(MBot2**3*Ss211*(6*Ss141**2 + Ss214)*Ss46*Ss47*
     &       (SA2/MHiggs2(1) + CA2/MHiggs2(2)))/
     &     (CB2*CW2*MZ2*Ss7*SW2) + 
     &   1/60.D0*(MBot2**3*Ss211*(6*Ss141**2 + Ss214)*Ss46*Ss47*
     &       (SA2/(Ss13*MHiggs2(1)) + CA2/(Ss1*MHiggs2(2))))/
     &     (CB2*CW2*Ss7*SW2) + 
     &   (CA2*MBot2**2*(-(3/40.D0*Ss56/(CW2*Ss7*SW2)) + 
     &        1/15.D0*(-3*Ss2 + MBot2*(66 + 72*Ss141 + 4*Ss214) + 
     &            24*Ss141*Ss270 - 12*Ss272 - 24*Ss549)/
     &          MHiggs2(2)**2))/(CB2*Ss1) + 
     &   (Ss51**2*(1/80.D0*(ML2*
     &            (-3 - 6*Ss306 + Ss269*(-6 + 12*Ss306) + Ss324)*
     &            Ss50)/Ss49 + 33/40.D0*ML2**3/(MZ2*MHiggs2(2))))/
     &    (CB2*CW2*Ss3*SW2) + 
     &   (Ss51**2*(-(1/80.D0*(ML2*Ss52*tmp3)/Ss49) - 
     &        1/20.D0*(ML2**3*(6*Ss306**2 + Ss324))/(Ss3*MHiggs2(2)))
     &      )/(CB2*CW2*Ss8*SW2) + 
     &   (Ss55**2*(3/10.D0*(ML2**2*Ss306*Ss311)/Ss15 - 
     &        1/20.D0*(ML2**3*(6*Ss306**2 + Ss324))))/
     &    (CB2*CW2*MZ2*Ss8*SW2*MHiggs2(1)) + 
     &   (CA*SA*(-(11/20.D0*(MTop2**3*Ss63)/MZ2) - 
     &        1/80.D0*(MTop2*(Ss19*(Ss19*Ss63 + 4*Ss369*Ss63) - 
     &              8*Ss46*Ss47*Ss48*Ss619*Ss64))/Ss44))/
     &    (CB*CW2*SB*Ss17*SW2*MHiggs2(1))

	gm2_2L = gm2_2L - 1/40.D0*
     &    (MBot2**2*Ss56*((4*SA2*Ss141*Ss146)/(Ss13*MHiggs2(1)) + 
     &         (11*MBot2*SA2*Ss61)/(Ss13*MHiggs2(1)) + 
     &         CA2*((4*Ss141*Ss144)/(Ss1*MHiggs2(2)) + 
     &            (11*MBot2*Ss59)/(Ss1*MHiggs2(2)))))/
     &     (CB2*CW2*MZ2*Ss7*SW2) + 
     &   1/30.D0*(CA2*MBot2*(33 + 12*(-2 + Ss141)*Ss141 + 12*Ss269 - 
     &         2*Ss274 + (-3*Ss2 - 24*Ss141*Ss271 - 12*Ss272 + 
     &            MBot2*(-156 + 72*Ss141 - 4*Ss214 - 24*Ss269 + 
     &               4*Ss274) + 24*Ss549)/MHiggs2(2)))/(CB2*Ss1) - 
     &   1/160.D0*(MBot2*((SA2*
     &            (Ss14*(Ss14*Ss56 + 4*Ss266*Ss56) + 
     &              8*Ss211*Ss46*Ss47*Ss547*Ss60))/
     &          (Ss13*Ss44*MHiggs2(1)) + 
     &         (CA2*(Ss2*(Ss2*Ss56 + 4*Ss272*Ss56) + 
     &              8*Ss211*Ss46*Ss47*Ss549*Ss57))/
     &          (Ss1*Ss49*MHiggs2(2))))/(CB2*CW2*SW2) + 
     &   1/160.D0*(ML2*((Ss55**2*
     &            (Ss16*(3*Ss16 + 12*Ss318) + 2*Ss54*Ss616))/
     &          (Ss15*Ss44*MHiggs2(1)) + 
     &         (Ss51**2*(12*Ss319*Ss4 + 3*Ss4**2 + 2*Ss50*Ss620))/
     &          (Ss3*Ss49*MHiggs2(2))))/(CB2*CW2*SW2) - 
     &   2/15.D0*(CA*MTop2*SA*(12*Ss269 - 2*Ss274 + 
     &         3*(11 + Ss335*(-8 + 4*Ss335)) + 
     &         (-3*Ss20 - 24*Ss335*Ss365 - 12*Ss370 + 
     &            MTop2*(-156 + 72*Ss335 - 4*Ss375) + 24*Ss623)/
     &          MHiggs2(2)))/(CB*SB*Ss18) + 
     &   (CA2*ML2*(1/10.D0*(-12*Ss269 + 2*Ss274 - 
     &           3*(11 + Ss306*(-8 + 4*Ss306))) + 
     &        1/10.D0*(24*Ss306*Ss316 + 12*Ss319 + 
     &            ML2*(156 - 72*Ss306 + 4*Ss324) + 3*Ss4 - 24*Ss621
     &            )/MHiggs2(2)))/(CB2*Ss3) + 
     &   (CA*SA*(-(1/80.D0*(MTop2*tmp4)/Ss44) - 
     &        1/30.D0*(MTop2**2*
     &            (MTop2*Ss367*(6*Ss335**2 + Ss375)*Ss46*Ss47*Ss48 - 
     &              6*Ss335*Ss339*Ss63))/(Ss18*MHiggs2(2))))/
     &    (CB*CW2*MZ2*SB*Ss9*SW2) + 
     &   (CA*MTop2*SA*(1/120.D0*
     &         ((Ss274*Ss46*Ss47*Ss48 + 
     &             (12*Ss269*Ss335 + Ss375)*Ss46*Ss47*Ss48 - 
     &             3*Ss63 - 6*Ss269*Ss63 - 6*Ss335*Ss63)*Ss65) + 
     &        1/80.D0*(Ss20*(Ss20*Ss63 + 4*Ss370*Ss63) - 
     &            8*Ss46*Ss47*Ss48*Ss623*Ss65)/MHiggs2(2)))/
     &    (CB*CW2*SB*Ss18*Ss49*SW2) + 
     &   (CA*SA*(CW2*MTop2**2*
     &         (-(88/5.D0*MTop2/MHiggs2(2)**2) + 
     &           4/5.D0*Ss20/MHiggs2(2)**2 - 
     &           32/5.D0*(Ss335*Ss363)/MHiggs2(2)**2 + 
     &           16/5.D0*Ss624/MHiggs2(2)**2) + 
     &        11/20.D0*(MTop2**3*Ss63)/(MZ2*SW2*MHiggs2(2))))/
     &    (CB*CW2*SB*Ss18) - 
     &   2/45.D0*(MTop2*(23 + 60*Ss138 + 6*Ss140 - 36*Ss335**2 + 
     &         (72*Ss335*Ss361 + 12*Ss368 + 
     &            MTop2*(46 + 120*Ss335 + 12*Ss375) - 72*Ss543 + 
     &            43*Ss6)/MHiggs2(3)))/Ss5 + 
     &   (Ss51**2*(3/10.D0*(ML2**2*Ss306*Ss310)/Ss3 - 
     &        1/20.D0*(ML2**3*(6*Ss306**2 + Ss324))))/
     &    (CB2*CW2*MZ2*Ss8*SW2*MHiggs2(2))

	gm2_2L = gm2_2L - 1/90.D0*
     &    (MBot2*TB2*(23 + 60*Ss138 + 6*Ss140 - 36*Ss141**2 + 
     &         (72*Ss141*Ss210 + 12*Ss212 + 
     &            MBot2*(46 + 120*Ss141 + 12*Ss214) - 72*Ss536 + 
     &            43*Ss67)/MHiggs2(3)))/Ss66 - 
     &   1/30.D0*(ML2*TB2*(23 + 60*Ss138 + 6*Ss140 - 36*Ss306**2 + 
     &         (72*Ss306*Ss312 + 12*Ss317 + 
     &            ML2*(46 + 120*Ss306 + 12*Ss324) - 72*Ss540 + 
     &            43*Ss69)/MHiggs2(3)))/Ss68 + 
     &   (TB2*(1/720.D0*(MBot2**2*
     &            (23 + Ss141*(60 + 36*Ss141) + 6*Ss214)*Ss56)/
     &          (Ss66*Ss7) - 
     &        1/480.D0*(ML2*Ss78**2*
     &            (2*Ss542 + 43*Ss69 + 
     &              6*(10*Ss306 + Ss138*(10 + 12*Ss306) + Ss324)*
     &               MHiggs2(3)))/(Ss62*Ss68)))/(CW2*SW2) - 
     &   1/1440.D0*(MBot2*Ss56*TB2*
     &       ((43*Ss10 + 6*MZ2*
     &             (Ss214 + 10*Ss474 + Ss141*(10 + 12*Ss474)) + 
     &            2*Ss554)/Ss7 + 
     &         (2*Ss539 + 43*Ss67 + 
     &            6*(10*Ss141 + Ss138*(10 + 12*Ss141) + Ss214)*
     &             MHiggs2(3))/Ss66))/(CW2*Ss62*SW2) - 
     &   1/720.D0*(MTop2*Ss63*((43*Ss12 + 
     &            6*MZ2*(Ss375 + 10*Ss474 + 
     &               Ss335*(10 + 12*Ss474)) + 2*Ss631)/Ss9 + 
     &         (2*Ss545 + 43*Ss6 + 
     &            6*(10*Ss335 + Ss138*(10 + 12*Ss335) + Ss375)*
     &             MHiggs2(3))/Ss5))/(CW2*Ss62*SW2) - 
     &   3/10.D0*(ML2**2*TB2*(Ss534 + (MW2**2*Ss529)/MHiggs2(4)**2))/
     &     (MW2**2*Ss25*SW2) + 
     &   3/10.D0*((Ss574*Ss578)/Ss37 + 
     &       (MW2**2*Ss563*Ss568)/(Ss24*MHiggs2(4)**2))/
     &     (MW2**2*Ss25*SW2) + 
     &   (ML2*TB2*(-(1/4.D0*Ss431/Ss25) + 
     &        3/40.D0*(ML2**2*Ss327)/MHiggs2(4)**2))/(MW2*SW2) + 
     &   (1/160.D0*(MBot2*(MTop2*(-(Ss334*Ss336) + 2*Ss388*Ss389) + 
     &             Ss397*Ss399 + 2*Ss462*Ss463) + 
     &          MTop2*(-(Ss407*Ss408) - 2*Ss472*Ss473))/Ss37**3 + 
     &      1/40.D0*(Ss515*Ss523)/(Ss30*MHiggs2(4)**2))/(MW2**2*SW2)+
     &     (MBot2*MTop2*(-(1/120.D0*Ss120) + 
     &        1/40.D0*(Ss214*Ss238 - 2*Ss209*Ss335 + 
     &            Ss141*(-2*Ss209 + 12*Ss238*Ss335) + Ss238*Ss375)/
     &          (MW2*MHiggs2(4))))/(Ss24*Ss30*Ss37*SW2) + 
     &   (1/1920.D0*(343*ML2**2 - 1101*MTop2*Ss29 + 1188*MTop2*Ss371 + 
     &         MBot2*(885*Ss29 + 4860*Ss371) - 16*Ss398*Ss634) + 
     &      3/40.D0*(ML2**3*Ss327*TB2)/MHiggs2(4))/(MW2**2*SW2) + 
     &   (ML2*TB2*(1/4.D0*Ss322/Ss25 + 
     &        1/240.D0*(204*Ss321 + ML2*(23 + 6*Ss324 - 18*Ss327))/
     &          MW2))/(SW2*MHiggs2(4)) + 
     &   (1/40.D0*(MBot2*(-(Ss214*Ss252) + 
     &             Ss252*(-12*Ss141*Ss275 - Ss295)) - 
     &          MTop2*Ss295*Ss305 - 12*MTop2*Ss275*Ss305*Ss335 - 
     &          MTop2*Ss305*Ss375)/Ss24 - 
     &      1/40.D0*(ML2*(ML2*Ss245*Ss295 + 4*Ss275*Ss307)*TB2))/
     &    (Ss25**3*SW2*MHiggs2(4))

	gm2_2L = gm2_2L + (-(1/720.D0*Ss137/(MZ2**3*Ss25**2*SW2**2)) + 
     &      1/120.D0*(-(Ss429*Log(MW2)) + 
     &          18*MW2**2*(Ss124*TB*Log(MW2)**2 + 
     &             Log(MHiggs2(4))*
     &              (Ss373 - Ss124*TB*Log(MHiggs2(4))))*MHiggs2(4))
     &         /(Ss25**3*SW2))/(CB*MW2**2*SB) + 
     &   (TB2*(1/40.D0*(ML2**2*
     &            (-6*Ss398*Ss409 + Ss255*Ss455 + 
     &              3*(1 + Ss306*(-2 + 4*Ss398))*MHiggs2(4)**2))/
     &          MW2 - 1/40.D0*
     &         (ML2*(-4*Ss398*(-MW2 + MHiggs2(4))*
     &              (-4*MW2 - 3*ML2*Ss306 + 4*MHiggs2(4)) + 
     &             3*ML2*(MW2*(-3 - 2*Ss306) + 4*MHiggs2(4) + 
     &                2*Ss275*MHiggs2(4))))))/(Ss25**3*SW2) + 
     &   (1/80.D0*(-(MTop2*(26*MTop2*Ss21*Ss26*Ss27 + 
     &               72*MTop2*Ss26*Ss27*Ss372 + 
     &               MBot2*(36*MTop2*Ss372 + MBot2*(Ss21 + 36*Ss372) + 
     &                  MTop2*Ss21*Ss79))) + 
     &          MBot2*(144*MTop2**2*Ss372 + 
     &             MBot2*(-61*MTop2*Ss21 - 180*MTop2*Ss372 + 
     &                MBot2*(17*Ss21 + 36*Ss372)))*TB2)/MHiggs2(4) + 
     &      1/240.D0*(MTop2*(228*Ss26*Ss27*Ss372 + 
     &            Ss21*(121*Ss26*Ss27 - MBot2*Ss80) + 
     &            MBot2*Ss372*(-600 - 372*TB2)) + 
     &         MBot2*TB2*(MBot2*(-188*Ss21 - 240*Ss372) + 
     &            137*Ss21*MHiggs2(4) + 132*Ss372*MHiggs2(4))))/
     &    (Ss24*Ss25*SW2) + 
     &   1/60.D0*(SA*Ss897*((MSbot2(1)*
     &            (-Ss898 + 8*Ss885*Ss899 - 4*Ss902 + 22*MSbot2(1))
     &            )/MHiggs2(1)**2 + 
     &         (Ss898 + 4*Ss900 + 
     &            (-27 + (8 - 4*Ss885)*Ss885)*MSbot2(1))/MHiggs2(1)
     &         ))/(CB2*CW*Ss892) + 
     &   1/180.D0*(CA*Ss884*((3*Ss888 + 12*Ss889 + 
     &            (-81 + 2*Ss274 + (24 - 12*Ss885)*Ss885 + 
     &               Ss269*(-12 + 24*Ss885))*MSbot2(1))/MHiggs2(2)+
     &           (MSbot2(1)*
     &            (-3*(Ss888 + 4*Ss889 + 8*Ss891) + 
     &              2*(33 + Ss885*(36 + 12*Ss885) + 2*Ss887)*
     &               MSbot2(1)))/MHiggs2(2)**2))/(CB2*CW*Ss878) + 
     &   1/60.D0*(SA*Ss915*((MSbot2(2)*
     &            (-Ss916 + 8*Ss906*Ss917 - 4*Ss920 + 22*MSbot2(2))
     &            )/MHiggs2(1)**2 + 
     &         (Ss916 + 4*Ss918 + 
     &            (-27 + (8 - 4*Ss906)*Ss906)*MSbot2(2))/MHiggs2(1)
     &         ))/(CB2*CW*Ss913) + 
     &   1/180.D0*(CA*Ss905*((3*Ss909 + 12*Ss910 + 
     &            (-81 + 2*Ss274 + (24 - 12*Ss906)*Ss906 + 
     &               Ss269*(-12 + 24*Ss906))*MSbot2(2))/MHiggs2(2)+
     &           (MSbot2(2)*
     &            (-3*(Ss909 + 4*Ss910 + 8*Ss912) + 
     &              2*(33 + Ss906*(36 + 12*Ss906) + 2*Ss908)*
     &               MSbot2(2)))/MHiggs2(2)**2))/(CB2*CW*Ss903) + 
     &   (-(23/240.D0*(ML2*Ss71*TB2)/Ss25) + 
     &      1/80.D0*(MTop2**2*(-26*Ss29 + 72*Ss371) + 
     &          MBot2*(MTop2*(108*Ss371 + Ss29*(-1 - 44*TB2)) + 
     &             MBot2*(17*Ss29*TB2 + 180*Ss371*TB2)))/Ss30)/
     &    (MW2*SW2*MHiggs2(4))

	gm2_2L = gm2_2L + (SA*Ss940*
     &      (-(3/20.D0*((Ss941 - 8*Ss928*Ss942 + 4*Ss945 - 
     &                22*MStau2(1))*MStau2(1))/MHiggs2(1)**2) + 
     &        3/20.D0*(Ss941 + 4*Ss943 + 
     &            (-27 + (8 - 4*Ss928)*Ss928)*MStau2(1))/MHiggs2(1)
     &        ))/(CB2*CW*Ss935) + 
     &   1/20.D0*(CA*Ss927*((3*Ss931 + 12*Ss932 + 
     &            (-81 + 2*Ss274 + (24 - 12*Ss928)*Ss928 + 
     &               Ss269*(-12 + 24*Ss928))*MStau2(1))/MHiggs2(2)+
     &           (MStau2(1)*
     &            (-3*(Ss931 + 4*Ss932 + 8*Ss934) + 
     &              2*(33 + Ss928*(36 + 12*Ss928) + 2*Ss930)*
     &               MStau2(1)))/MHiggs2(2)**2))/(CB2*CW*Ss921) + 
     &   (SA*Ss958*(-(3/20.D0*
     &           ((Ss959 - 8*Ss949*Ss960 + 4*Ss963 - 22*MStau2(2))*
     &              MStau2(2))/MHiggs2(1)**2) + 
     &        3/20.D0*(Ss959 + 4*Ss961 + 
     &            (-27 + (8 - 4*Ss949)*Ss949)*MStau2(2))/MHiggs2(1)
     &        ))/(CB2*CW*Ss956) + 
     &   1/20.D0*(CA*Ss948*((3*Ss952 + 12*Ss953 + 
     &            (-81 + 2*Ss274 + (24 - 12*Ss949)*Ss949 + 
     &               Ss269*(-12 + 24*Ss949))*MStau2(2))/MHiggs2(2)+
     &           (MStau2(2)*
     &            (-3*(Ss952 + 4*Ss953 + 8*Ss955) + 
     &              2*(33 + Ss949*(36 + 12*Ss949) + 2*Ss951)*
     &               MStau2(2)))/MHiggs2(2)**2))/(CB2*CW*Ss946) + 
     &   (CA*(1/80.D0*(MTop2*SA*tmp4)/(CW2*MZ2*Ss49*Ss9*SW2) + 
     &        1/45.D0*(Ss968*
     &            (3*Ss972 + 12*Ss973 + 
     &              (-81 + 2*Ss274 + (24 - 12*Ss969)*Ss969 + 
     &                 Ss269*(-12 + 24*Ss969))*MStop2(1)))/
     &          (CW*Ss964*MHiggs2(2))))/(CB*SB) + 
     &   1/45.D0*(SA*Ss981*((-3*(Ss982 + 4*Ss984) + 
     &            (12*Ss263 - 2*Ss268 + 
     &               3*(27 + Ss969*(-8 + 4*Ss969)))*MStop2(1))/
     &          MHiggs2(1) + 
     &         (MStop2(1)*(3*Ss982 + 12*Ss984 + 24*Ss985 + 
     &              Ss969*(-24*Ss983 - 72*MStop2(1)) + 
     &              (-66 - 4*Ss971)*MStop2(1)))/MHiggs2(1)**2))/
     &     (CB*CW*SB*Ss976) + 
     &   1/45.D0*(SA*Ss997*MStop2(2)*
     &       (12*Ss1000 + 24*Ss1001 + 3*Ss998 + 
     &         Ss988*(-24*Ss999 - 72*MStop2(2)) + 
     &         (-66 - 4*Ss990)*MStop2(2)))/
     &     (CB*CW*SB*Ss995*MHiggs2(1)**2) + 
     &   1/45.D0*((CA*Ss987*(3*Ss991 + 12*Ss992 + 
     &            (-81 + 2*Ss274 + (24 - 12*Ss988)*Ss988 + 
     &               Ss269*(-12 + 24*Ss988))*MStop2(2)))/
     &        (Ss986*MHiggs2(2)) + 
     &       (SA*Ss997*(-3*(4*Ss1000 + Ss998) + 
     &            (12*Ss263 - 2*Ss268 + 
     &               3*(27 + Ss988*(-8 + 4*Ss988)))*MStop2(2)))/
     &        (Ss995*MHiggs2(1)))/(CB*CW*SB) + 
     &   1/45.D0*(CA*((Ss968*MStop2(1)*
     &            (-3*(Ss972 + 4*Ss973 + 8*Ss975) + 
     &              2*(33 + Ss969*(36 + 12*Ss969) + 2*Ss971)*
     &               MStop2(1)))/Ss964 + 
     &         (Ss987*MStop2(2)*
     &            (-3*(Ss991 + 4*Ss992 + 8*Ss994) + 
     &              2*(33 + Ss988*(36 + 12*Ss988) + 2*Ss990)*
     &               MStop2(2)))/Ss986))/(CB*CW*SB*MHiggs2(2)**2)

	end

