* gm2_2L_full4.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 30 Nov 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full4
	implicit none

#include "gm2_2L_full.h"

	ComplexType tmp10, tmp11, tmp12, tmp13, tmp9

	Sf6529 = 6*MHiggs2(4)**2 + 
     &   MSbot2(1)*(-13*MSbot2(1) - 62*MStop2(2)) - 21*MStop2(2)**2

	Sf4288 = 6*MHiggs2(4)**2 + 
     &   MSbot2(1)*(-13*MSbot2(1) - 62*MStop2(1)) - 21*MStop2(1)**2

	Sf2266 = Sf2209 - Sf2201*MSbot2(2)

	Sf2202 = Sf2209 - Sf2201*MSbot2(1)

	Sf2267 = 2*MW2**2 + MHiggs2(4)*(-3*MW2 + 4*MHiggs2(4)) + 
     &   Sf2201*MSbot2(2)

	Sf2203 = 2*MW2**2 + MHiggs2(4)*(-3*MW2 + 4*MHiggs2(4)) + 
     &   Sf2201*MSbot2(1)

	Sf4873 = 2*MW2*Sf99 + 
     &   MHiggs2(4)*(-2*Sf4872 + 3*MHiggs2(4)) + MSq2Diff(2,2)**2

	Sf4838 = 2*MW2*Sf91 + 
     &   MHiggs2(4)*(-2*Sf4837 + 3*MHiggs2(4)) + MSq2Diff(1,2)**2

	Sf2276 = 2*MW2*Sf77 + 
     &   MHiggs2(4)*(-2*Sf2275 + 3*MHiggs2(4)) + MSq2Diff(2,1)**2

	Sf2213 = 2*MW2*Sf69 + 
     &   MHiggs2(4)*(-2*Sf2212 + 3*MHiggs2(4)) + MSq2Diff(1,1)**2

	Sf4874 = MW2*Sf4868 - Sf4872*MHiggs2(4) + 
     &   MHiggs2(4)**2 + MSq2Diff(2,2)**2

	Sf4869 = -2*Sf4868*MHiggs2(4) + MHiggs2(4)**2 + 
     &   MSq2Diff(2,2)**2 + 4*MW2*MStop2(2)

	Sf4840 = MW2*Sf4839 - Sf4837*MHiggs2(4) + 
     &   MHiggs2(4)**2 + MSq2Diff(1,2)**2

	Sf2277 = MW2*Sf2271 - Sf2275*MHiggs2(4) + 
     &   MHiggs2(4)**2 + MSq2Diff(2,1)**2

	Sf2272 = -2*Sf2271*MHiggs2(4) + MHiggs2(4)**2 + 
     &   MSq2Diff(2,1)**2 + 4*MW2*MStop2(1)

	Sf2214 = MW2*Sf2207 - Sf2212*MHiggs2(4) + 
     &   MHiggs2(4)**2 + MSq2Diff(1,1)**2

	Sf2208 = -2*Sf2207*MHiggs2(4) + MHiggs2(4)**2 + 
     &   MSq2Diff(1,1)**2 + 4*MW2*MStop2(1)

	Sf9719 = -332*MZ2**2 + 
     &   MHiggs2(2)*(62*MZ2 + 224*MHiggs2(2)) + Sf9467*MStau2(2)

	Sf9468 = -332*MZ2**2 + 
     &   MHiggs2(2)*(62*MZ2 + 224*MHiggs2(2)) + Sf9467*MStau2(1)

	Sf9748 = -44*MZ2**2 + 
     &   MHiggs2(2)*(-4*MZ2 + 32*MHiggs2(2)) + Sf9526*MStau2(2)

	Sf9527 = -44*MZ2**2 + 
     &   MHiggs2(2)*(-4*MZ2 + 32*MHiggs2(2)) + Sf9526*MStau2(1)

	Sf10517 = -MZ2**2 + Sf8678*MHiggs2(2) + MHiggs2(2)**2 + 
     &   6*MStop2(2)**2

	Sf10454 = -MZ2**2 + Sf7992*MHiggs2(2) + MHiggs2(2)**2 + 
     &   6*MStop2(1)**2

	Sf9015 = -MZ2**2 + Sf7523*MHiggs2(2) + MHiggs2(2)**2 + 
     &   6*MSbot2(2)**2

	Sf8937 = -MZ2**2 + Sf7504*MHiggs2(2) + MHiggs2(2)**2 + 
     &   6*MSbot2(1)**2

	Sf10412 = -MZ2**2 + Sf8678*MHiggs2(1) + MHiggs2(1)**2 + 
     &   6*MStop2(2)**2

	Sf10359 = -MZ2**2 + Sf7992*MHiggs2(1) + MHiggs2(1)**2 + 
     &   6*MStop2(1)**2

	Sf10148 = -MZ2**2 + Sf7734*MHiggs2(1) + MHiggs2(1)**2 + 
     &   6*MStau2(2)**2

	Sf10078 = -MZ2**2 + Sf7581*MHiggs2(1) + MHiggs2(1)**2 + 
     &   6*MStau2(1)**2

	Sf9213 = -MZ2**2 + Sf7523*MHiggs2(1) + MHiggs2(1)**2 + 
     &   6*MSbot2(2)**2

	Sf9139 = -MZ2**2 + Sf7504*MHiggs2(1) + MHiggs2(1)**2 + 
     &   6*MSbot2(1)**2

	Sf5409 = 4*MBot2 - 4*MTop2 + C2B*(35*MSbot2(2) + 13*MStop2(2))

	Sf3001 = 4*MBot2 - 4*MTop2 + C2B*(35*MSbot2(2) + 13*MStop2(1))

	Sf6889 = 4*MBot2 - 4*MTop2 + C2B*(23*MSbot2(2) + 25*MStop2(2))

	Sf4684 = 4*MBot2 - 4*MTop2 + C2B*(23*MSbot2(2) + 25*MStop2(1))

	Sf5303 = 4*MBot2 - 4*MTop2 + C2B*(35*MSbot2(1) + 13*MStop2(2))

	Sf2863 = 4*MBot2 - 4*MTop2 + C2B*(35*MSbot2(1) + 13*MStop2(1))

	Sf6638 = 4*MBot2 - 4*MTop2 + C2B*(23*MSbot2(1) + 25*MStop2(2))

	Sf4579 = 4*MBot2 - 4*MTop2 + C2B*(23*MSbot2(1) + 25*MStop2(1))

	Sf7487 = 16 + UStop22**2*(-24 + 9*UStop22**2)

	Sf7459 = 1 + UStau22**2*(-2 + 5*UStau22**2)

	Sf7456 = 4 + UStau22**2*(-8 + 5*UStau22**2)

	Sf8604 = -7*Sf8543*MStau2(1) + 10*MSneu2(1)*UStau11 + 
     &   14*MStau2(2)*UStau12*UStau21*UStau22

	Sf7483 = -8 - USbot22**2 + 9*UStop22**2

	Sf7475 = -Sf7483

	Sf8741 = -MStop2(1) + MSbot2(2)*USbot21**2 + 
     &   MSbot2(1)*USbot22**2

	Sf7267 = 3*MStop2(2) + Sf7265*USbot22**2 + 
     &   Sf7266*(1 - USbot22**2)

	Sf10583 = Sf10582*MStop2(1)**2 - MStop2(1)**3 - 
     &   Sf10244*MStop2(2)**2

	Sf10432 = Sf10431*MStop2(1)**2 - MStop2(1)**3 - 
     &   Sf10283*MStop2(2)**2

	Sf8527 = -(3/2.D0*MStop2(4)) + Log(MStop2(1))*MStop2(1) + 
     &   Log(MStop2(2))*MStop2(2)

	Sf10563 = Sf10562*MStop2(1) + 30*MHiggs2(2)*MStop2(2) + 
     &   8*MStop2(2)**2

	Sf10497 = -Sf8528 + 2*Sf10476*MStop2(1)

	Sf10465 = 2*Sf10463*MStop2(1) - MStop2(1)**2 + Sf10464*MStop2(2)

	Sf10316 = Sf10315*MStop2(1) + 30*MHiggs2(1)*MStop2(2) + 
     &   8*MStop2(2)**2

	Sf10388 = -Sf8528 + 2*Sf10284*MStop2(1)

	Sf10555 = 30*MHiggs2(2)*MStop2(1) + 8*MStop2(1)**2 + 
     &   Sf10554*MStop2(2)

	Sf10307 = 30*MHiggs2(1)*MStop2(1) + 8*MStop2(1)**2 + 
     &   Sf10306*MStop2(2)

	Sf7264 = 3*MStop2(1) + Sf7262*USbot22**2 + 
     &   Sf7263*(1 - USbot22**2)

	Sf9857 = Sf9856*MStau2(1)**2 - MStau2(1)**3 - 
     &   Sf9359*MStau2(2)**2

	Sf10208 = Sf10207*MStau2(1)**2 - MStau2(1)**3 - 
     &   Sf10095*MStau2(2)**2

	Sf7708 = -(3/2.D0*MStau2(4)) + Log(MStau2(1))*MStau2(1) + 
     &   Log(MStau2(2))*MStau2(2)

	Sf9374 = -Sf7709 + 2*Sf9360*MStau2(1)

	Sf10155 = Sf10154*MStau2(1) + 30*MHiggs2(1)*MStau2(2) + 
     &   8*MStau2(2)**2

	Sf10129 = -Sf7709 + 2*Sf10096*MStau2(1)

	Sf10125 = 30*MHiggs2(1)*MStau2(1) + 8*MStau2(1)**2 + 
     &   Sf10124*MStau2(2)

	Sf1615 = Sf778*MSneu2(1)**2 + 
     &   Sf779*MSneu2(1)*MStau2(2) + Sf780*MStau2(2)**2

	Sf781 = Sf778*MSneu2(1)**2 + Sf779*MSneu2(1)*MStau2(1) + 
     &   Sf780*MStau2(1)**2

	Sf1504 = Sf499*MSneu2(1)**2 + 6*MSneu2(1)*MStau2(2) + 
     &   Sf607*MStau2(2)**2

	Sf608 = Sf499*MSneu2(1)**2 + 6*MSneu2(1)*MStau2(1) + 
     &   Sf607*MStau2(1)**2

	Sf7436 = MSneu2(1) + MStau2(2)*UStau21**2 + 
     &   MStau2(1)*UStau22**2

	Sf2111 = Sf2109*MSneu2(1) + 8*MSneu2(1)**2 + Sf2110*MStau2(2)

	Sf2027 = Sf2025*MSneu2(1) + 8*MSneu2(1)**2 + Sf2026*MStau2(1)

	Sf6911 = Sf4596*MSbot2(2)**2 + 
     &   Sf4598*MSbot2(2)*MStop2(2) + 5*Sf4600*MStop2(2)**2

	Sf4718 = Sf4596*MSbot2(2)**2 + 
     &   Sf4598*MSbot2(2)*MStop2(1) + 5*Sf4600*MStop2(1)**2

	Sf6163 = Sf3787*MSbot2(2)**2 + 
     &   2*Sf3788*MSbot2(2)*MStop2(2) + Sf3789*MStop2(2)**2

	Sf3915 = Sf3787*MSbot2(2)**2 + 
     &   2*Sf3788*MSbot2(2)*MStop2(1) + Sf3789*MStop2(1)**2

	Sf6178 = 3*Sf3778*MSbot2(2)**2 + 
     &   Sf3779*MSbot2(2)*MStop2(2) - Sf3780*MStop2(2)**2

	Sf3933 = 3*Sf3778*MSbot2(2)**2 + 
     &   Sf3779*MSbot2(2)*MStop2(1) - Sf3780*MStop2(1)**2

	Sf6150 = 3*Sf3760*MSbot2(2)**2 - 
     &   Sf3761*MSbot2(2)*MStop2(2) - Sf3762*MStop2(2)**2

	Sf3902 = 3*Sf3760*MSbot2(2)**2 - 
     &   Sf3761*MSbot2(2)*MStop2(1) - Sf3762*MStop2(1)**2

	Sf6045 = Sf3435*MSbot2(2)**2 + 14*MSbot2(2)*MStop2(2) + 
     &   Sf3436*MStop2(2)**2

	Sf3709 = Sf3435*MSbot2(2)**2 + 14*MSbot2(2)*MStop2(1) + 
     &   Sf3436*MStop2(1)**2

	Sf5179 = -Sf2733 + 2*Sf5172*MSbot2(2)

	Sf5197 = -Sf5141 + 2*Sf4833*MSbot2(2)

	Sf9276 = Sf9275*MSbot2(1)**2 - MSbot2(1)**3 - 
     &   Sf9157*MSbot2(2)**2

	Sf9095 = Sf9094*MSbot2(1)**2 - MSbot2(1)**3 - 
     &   Sf8967*MSbot2(2)**2

	Sf5249 = Sf5248*MSbot2(1)**2 + 
     &   Sf2787*MSbot2(1)*MStop2(2) + Sf2788*MStop2(2)**2

	Sf6660 = Sf4596*MSbot2(1)**2 + 
     &   Sf4598*MSbot2(1)*MStop2(2) + 5*Sf4600*MStop2(2)**2

	Sf4623 = Sf4596*MSbot2(1)**2 + 
     &   Sf4598*MSbot2(1)*MStop2(1) + 5*Sf4600*MStop2(1)**2

	Sf3790 = Sf3787*MSbot2(1)**2 + 
     &   2*Sf3788*MSbot2(1)*MStop2(1) + Sf3789*MStop2(1)**2

	Sf5800 = 3*Sf3778*MSbot2(1)**2 + 
     &   Sf3779*MSbot2(1)*MStop2(2) - Sf3780*MStop2(2)**2

	Sf3809 = 3*Sf3778*MSbot2(1)**2 + 
     &   Sf3779*MSbot2(1)*MStop2(1) - Sf3780*MStop2(1)**2

	Sf5767 = 3*Sf3760*MSbot2(1)**2 - 
     &   Sf3761*MSbot2(1)*MStop2(2) - Sf3762*MStop2(2)**2

	Sf3763 = 3*Sf3760*MSbot2(1)**2 - 
     &   Sf3761*MSbot2(1)*MStop2(1) - Sf3762*MStop2(1)**2

	Sf3437 = Sf3435*MSbot2(1)**2 + 14*MSbot2(1)*MStop2(1) + 
     &   Sf3436*MStop2(1)**2

	Sf5583 = Sf3270*MSbot2(1)**2 + 174*MSbot2(1)*MStop2(2) + 
     &   Sf3271*MStop2(2)**2

	Sf5521 = Sf3150*MSbot2(1)**2 + 78*MSbot2(1)*MStop2(2) + 
     &   Sf3151*MStop2(2)**2

	Sf7516 = -(3/2.D0*MSbot2(4)) + Log(MSbot2(1))*MSbot2(1) + 
     &   Log(MSbot2(2))*MSbot2(2)

	Sf9219 = Sf9218*MSbot2(1) + 30*MHiggs2(1)*MSbot2(2) + 
     &   8*MSbot2(2)**2

	Sf9194 = -Sf7517 + 2*Sf9158*MSbot2(1)

	Sf9070 = Sf9069*MSbot2(1) + 30*MHiggs2(2)*MSbot2(2) + 
     &   8*MSbot2(2)**2

	Sf8997 = -Sf7517 + 2*Sf8968*MSbot2(1)

	Sf8952 = 2*Sf8950*MSbot2(1) - MSbot2(1)**2 + Sf8951*MSbot2(2)

	Sf8895 = 2*Sf8885*Sf8886*MSbot2(1) + MSbot2(1)**2 + 
     &   Sf8894*MSbot2(2)

	Sf5173 = -Sf2555 + 2*Sf5172*MSbot2(1)

	Sf9061 = 30*MHiggs2(2)*MSbot2(1) + 8*MSbot2(1)**2 + 
     &   Sf9060*MSbot2(2)

	Sf9190 = 30*MHiggs2(1)*MSbot2(1) + 8*MSbot2(1)**2 + 
     &   Sf9189*MSbot2(2)

	Sf8815 = (MZ2*Sf10415 - Sf8814*MStop2(3))/MZ2

	Sf8785 = (MZ2*Sf9840 - Sf8784*MStau2(3))/MZ2

	Sf8781 = (MW2*Sf1519 - Sf1518*MSl2Diff(2,1))/MW2

	Sf8774 = (MW2*Sf624 - Sf623*MSl2Diff(1,1))/MW2

	Sf8811 = (MW2*Sf5164 + Sf5163*MSq2Diff(2,2))/MW2

	Sf8798 = (MW2*Sf2757 + Sf2756*MSq2Diff(2,1))/MW2

	Sf8805 = (MW2*Sf5014 + Sf5013*MSq2Diff(1,2))/MW2

	Sf8792 = (MW2*Sf2581 + Sf2580*MSq2Diff(1,1))/MW2

	Sf8767 = (MZ2*Sf9077 - Sf8766*MSbot2(3))/MZ2

	Sf6946 = Sf6945*TB2 + reMUE*TB*UStop21 + MTop*UStop22

	Sf4228 = Sf4227*TB2 + reMUE*TB*UStop11 + MTop*UStop12

	Sf10470 = Ss484 + 12*Ss474*Ss988 + Ss990

	Sf10436 = Ss484 + 12*Ss474*Ss969 + Ss971

	Sf9349 = Ss484 + 12*Ss474*Ss949 + Ss951

	Sf9313 = Ss484 + 12*Ss474*Ss928 + Ss930

	Sf8957 = Ss484 + 12*Ss474*Ss906 + Ss908

	Sf8900 = Ss484 + 12*Ss474*Ss885 + Ss887

	Sf5999 = 45*Ss30 - Sf3261*MSbot2(2) + Sf3262*MStop2(2)

	Sf3594 = 45*Ss30 - Sf3261*MSbot2(2) + Sf3262*MStop2(1)

	Sf5580 = 45*Ss30 - Sf3261*MSbot2(1) + Sf3262*MStop2(2)

	Sf3263 = 45*Ss30 - Sf3261*MSbot2(1) + Sf3262*MStop2(1)

	Sf5064 = 21*Ss30 - Sf2383*MSbot2(2) + Sf765*MStop2(2)

	Sf2613 = 21*Ss30 - Sf2383*MSbot2(2) + Sf765*MStop2(1)

	Sf4938 = 21*Ss30 - Sf2383*MSbot2(1) + Sf765*MStop2(2)

	Sf2384 = 21*Ss30 - Sf2383*MSbot2(1) + Sf765*MStop2(1)

	Sf6208 = 18*Ss30 - Sf4035*MSbot2(2) - 39*Sf234*MStop2(2)

	Sf4147 = 18*Ss30 - Sf4035*MSbot2(2) - 39*Sf234*MStop2(1)

	Sf5830 = 18*Ss30 - Sf4035*MSbot2(1) - 39*Sf234*MStop2(2)

	Sf4036 = 18*Ss30 - Sf4035*MSbot2(1) - 39*Sf234*MStop2(1)

	Sf6837 = 15*Ss30 - Sf6446*MSbot2(2) - Sf3231*MStop2(2)

	Sf6715 = 15*Ss30 - Sf6446*MSbot2(2) - Sf3231*MStop2(1)

	Sf7024 = 15*Ss30 + Sf4401*MSbot2(2) + 3*Sf4402*MStop2(2)

	Sf4508 = 15*Ss30 + Sf4401*MSbot2(2) + 3*Sf4402*MStop2(1)

	Sf6582 = 15*Ss30 - Sf6446*MSbot2(1) - Sf3231*MStop2(2)

	Sf6447 = 15*Ss30 - Sf6446*MSbot2(1) - Sf3231*MStop2(1)

	Sf7002 = 15*Ss30 + Sf4401*MSbot2(1) + 3*Sf4402*MStop2(2)

	Sf4403 = 15*Ss30 + Sf4401*MSbot2(1) + 3*Sf4402*MStop2(1)

	Sf5094 = 12*Ss30 - Sf2441*MSbot2(2) - Sf2442*MStop2(2)

	Sf2654 = 12*Ss30 - Sf2441*MSbot2(2) - Sf2442*MStop2(1)

	Sf4968 = 12*Ss30 - Sf2441*MSbot2(1) - Sf2442*MStop2(2)

	Sf2443 = 12*Ss30 - Sf2441*MSbot2(1) - Sf2442*MStop2(1)

	Sf10475 = Ss274 + 12*Ss269*Ss988 + Ss990

	Sf10441 = Ss274 + 12*Ss269*Ss969 + Ss971

	Sf9358 = Ss274 + 12*Ss269*Ss949 + Ss951

	Sf9327 = Ss274 + 12*Ss269*Ss928 + Ss930

	Sf8966 = Ss274 + 12*Ss269*Ss906 + Ss908

	Sf8916 = Ss274 + 12*Ss269*Ss885 + Ss887

	Sf10378 = Ss268 + 12*Ss263*Ss988 + Ss990

	Sf10351 = Ss268 + 12*Ss263*Ss969 + Ss971

	Sf10094 = Ss268 + 12*Ss263*Ss949 + Ss951

	Sf10064 = Ss268 + 12*Ss263*Ss928 + Ss930

	Sf9156 = Ss268 + 12*Ss263*Ss906 + Ss908

	Sf9118 = Ss268 + 12*Ss263*Ss885 + Ss887

	Sf10236 = Ss982 + 4*Ss984 + 4*Ss263*MHiggs2(1)

	Sf5068 = 42*Ss215*Ss218 + Sf2395*MSbot2(2) + 9*Sf2396*MStop2(2)

	Sf2619 = 42*Ss215*Ss218 + Sf2395*MSbot2(2) + 9*Sf2396*MStop2(1)

	Sf4942 = 42*Ss215*Ss218 + Sf2395*MSbot2(1) + 9*Sf2396*MStop2(2)

	Sf2397 = 42*Ss215*Ss218 + Sf2395*MSbot2(1) + 9*Sf2396*MStop2(1)

	Sf5990 = 30*Ss215*Ss218 - Sf3228*MSbot2(2) + Sf3229*MStop2(2)

	Sf3580 = 30*Ss215*Ss218 - Sf3228*MSbot2(2) + Sf3229*MStop2(1)

	Sf5564 = 30*Ss215*Ss218 - Sf3228*MSbot2(1) + Sf3229*MStop2(2)

	Sf3230 = 30*Ss215*Ss218 - Sf3228*MSbot2(1) + Sf3229*MStop2(1)

	Sf5106 = 21*Ss215*Ss218 - Sf2473*MSbot2(2) - Sf2474*MStop2(2)

	Sf2673 = 21*Ss215*Ss218 - Sf2473*MSbot2(2) - Sf2474*MStop2(1)

	Sf4980 = 21*Ss215*Ss218 - Sf2473*MSbot2(1) - Sf2474*MStop2(2)

	Sf2475 = 21*Ss215*Ss218 - Sf2473*MSbot2(1) - Sf2474*MStop2(1)

	Sf6788 = 12*Ss215*Ss218 + 31*MSbot2(2) + 76*MStop2(2)

	Sf4433 = 12*Ss215*Ss218 + 31*MSbot2(2) + 76*MStop2(1)

	Sf6533 = 12*Ss215*Ss218 + 31*MSbot2(1) + 76*MStop2(2)

	Sf4305 = 12*Ss215*Ss218 + 31*MSbot2(1) + 76*MStop2(1)

	Sf5122 = 9*Ss215*Ss218 + Sf2514*MSbot2(2) + 41*MStop2(2)

	Sf2701 = 9*Ss215*Ss218 + Sf2514*MSbot2(2) + 41*MStop2(1)

	Sf6830 = 9*Ss215*Ss218 - Sf2459*MSbot2(2) + Sf4394*MStop2(2)

	Sf6705 = 9*Ss215*Ss218 - Sf2459*MSbot2(2) + Sf4394*MStop2(1)

	Sf7018 = 9*Ss215*Ss218 - Sf2348*MSbot2(2) + 7*Sf169*MStop2(2)

	Sf4497 = 9*Ss215*Ss218 - Sf2348*MSbot2(2) + 7*Sf169*MStop2(1)

	Sf4996 = 9*Ss215*Ss218 + Sf2514*MSbot2(1) + 41*MStop2(2)

	Sf2515 = 9*Ss215*Ss218 + Sf2514*MSbot2(1) + 41*MStop2(1)

	Sf6575 = 9*Ss215*Ss218 - Sf2459*MSbot2(1) + Sf4394*MStop2(2)

	Sf6435 = 9*Ss215*Ss218 - Sf2459*MSbot2(1) + Sf4394*MStop2(1)

	Sf6996 = 9*Ss215*Ss218 - Sf2348*MSbot2(1) + 7*Sf169*MStop2(2)

	Sf4388 = 9*Ss215*Ss218 - Sf2348*MSbot2(1) + 7*Sf169*MStop2(1)

	Sf6824 = 6*Ss215*Ss218 - 2*Sf6423*MSbot2(2) - 3*Sf1095*MStop2(2)

	Sf6696 = 6*Ss215*Ss218 - 2*Sf6423*MSbot2(2) - 3*Sf1095*MStop2(1)

	Sf7013 = 6*Ss215*Ss218 + Sf4377*MSbot2(2) + 4*Sf4112*MStop2(2)

	Sf4487 = 6*Ss215*Ss218 + Sf4377*MSbot2(2) + 4*Sf4112*MStop2(1)

	Sf5934 = 6*Ss215*Ss218 - Sf3108*MSbot2(2) + Sf3109*MStop2(2)

	Sf3491 = 6*Ss215*Ss218 - Sf3108*MSbot2(2) + Sf3109*MStop2(1)

	Sf6569 = 6*Ss215*Ss218 - 2*Sf6423*MSbot2(1) - 3*Sf1095*MStop2(2)

	Sf6424 = 6*Ss215*Ss218 - 2*Sf6423*MSbot2(1) - 3*Sf1095*MStop2(1)

	Sf6991 = 6*Ss215*Ss218 + Sf4377*MSbot2(1) + 4*Sf4112*MStop2(2)

	Sf4378 = 6*Ss215*Ss218 + Sf4377*MSbot2(1) + 4*Sf4112*MStop2(1)

	Sf5501 = 6*Ss215*Ss218 - Sf3108*MSbot2(1) + Sf3109*MStop2(2)

	Sf3110 = 6*Ss215*Ss218 - Sf3108*MSbot2(1) + Sf3109*MStop2(1)

	Sf6031 = 3*Ss215*Ss218 - Sf765*MSbot2(2) + Sf1050*MStop2(2)

	Sf3690 = 3*Ss215*Ss218 - Sf765*MSbot2(2) + Sf1050*MStop2(1)

	Sf6040 = 3*Ss215*Ss218 - Sf747*MSbot2(2) - Sf3423*MStop2(2)

	Sf3700 = 3*Ss215*Ss218 - Sf747*MSbot2(2) - Sf3423*MStop2(1)

	Sf5941 = 3*Ss215*Ss218 + Sf1116*MSbot2(2) - Sf3138*MStop2(2)

	Sf3505 = 3*Ss215*Ss218 + Sf1116*MSbot2(2) - Sf3138*MStop2(1)

	Sf5611 = 3*Ss215*Ss218 - Sf765*MSbot2(1) + Sf1050*MStop2(2)

	Sf3412 = 3*Ss215*Ss218 - Sf765*MSbot2(1) + Sf1050*MStop2(1)

	Sf5620 = 3*Ss215*Ss218 - Sf747*MSbot2(1) - Sf3423*MStop2(2)

	Sf3424 = 3*Ss215*Ss218 - Sf747*MSbot2(1) - Sf3423*MStop2(1)

	Sf5516 = 3*Ss215*Ss218 + Sf1116*MSbot2(1) - Sf3138*MStop2(2)

	Sf3139 = 3*Ss215*Ss218 + Sf1116*MSbot2(1) - Sf3138*MStop2(1)

	Sf6811 = 2*Ss215*Ss218 + Sf4359*MSbot2(2) + Sf4360*MStop2(2)

	Sf4472 = 2*Ss215*Ss218 + Sf4359*MSbot2(2) + Sf4360*MStop2(1)

	Sf6556 = 2*Ss215*Ss218 + Sf4359*MSbot2(1) + Sf4360*MStop2(2)

	Sf4361 = 2*Ss215*Ss218 + Sf4359*MSbot2(1) + Sf4360*MStop2(1)

	Sf10256 = 4*Ss1000 + Ss998 + 4*Ss263*MHiggs2(1)

	Sf6046 = 9*MBot2*Sf99**2 - 9*MTop2*Sf99**2 + Sf6045*MSq2Diff(2,2)

	Sf6153 = Sf3765*MSbot2(2)**4 + 
     &   12*MBot2*MTop2*Sf99*MSq2Diff(2,2) - Sf6152*MStop2(2)

	Sf6912 = -5*Sf4543*MSbot2(2)**3 + 
     &   20*MTop2*Sf99*MSq2Diff(2,2) - Sf6911*MStop2(2)

	Sf6179 = -(Sf3777*MSbot2(2)**3) + 
     &   12*MBot2*Sf99*MSq2Diff(2,2) + Sf6178*MStop2(2)

	Sf5139 = 3/2.D0*Sf99 - Log(MSbot2(2))*MSbot2(2) - 
     &   Log(MStop2(2))*MStop2(2)

	Sf9879 = Sf9875 - Sf9876*MStau2(2) + 2*Sf9877*MStau2(2)**2

	Sf9681 = Sf9679 - 4*Sf9680*MStau2(2) + 10*MStau2(2)**2

	Sf9584 = Sf9580 - Sf9582*MStau2(2) + Sf9583*MStau2(2)**2

	Sf1752 = CA*ML*Sf945 - CAB*CB*Ss443*UStau11 - 
     &   ML*reMUE*SA*UStau12

	Sf5584 = 63*MBot2*Sf91**2 - 63*MTop2*Sf91**2 - Sf5583*MSq2Diff(1,2)

	Sf5522 = 15*MBot2*Sf91**2 - 15*MTop2*Sf91**2 - Sf5521*MSq2Diff(1,2)

	Sf6661 = -5*Sf4543*MSbot2(1)**3 + 
     &   20*MTop2*Sf91*MSq2Diff(1,2) - Sf6660*MStop2(2)

	Sf5801 = -(Sf3777*MSbot2(1)**3) + 
     &   12*MBot2*Sf91*MSq2Diff(1,2) + Sf5800*MStop2(2)

	Sf5020 = 3/2.D0*Sf91 - Log(MSbot2(1))*MSbot2(1) - 
     &   Log(MStop2(2))*MStop2(2)

	Sf8681 = 3/2.D0*Sf8678 - MZ2*Log(MZ2) - 2*Log(MStop2(2))*MStop2(2)

	Sf1667 = 2*Sf850 - Sf191*MStau2(2) - 4*MStau2(2)**2

	Sf851 = 2*Sf850 - Sf191*MStau2(1) - 4*MStau2(1)**2

	Sf3605 = 2*Sf850 - Sf191*MSbot2(2) - 4*MSbot2(2)**2

	Sf3278 = 2*Sf850 - Sf191*MSbot2(1) - 4*MSbot2(1)**2

	Sf7995 = 3/2.D0*Sf7992 - MZ2*Log(MZ2) - 2*Log(MStop2(1))*MStop2(1)

	Sf7737 = 3/2.D0*Sf7734 - MZ2*Log(MZ2) - 2*Log(MStau2(2))*MStau2(2)

	Sf3710 = 9*MBot2*Sf77**2 - 9*MTop2*Sf77**2 + Sf3709*MSq2Diff(2,1)

	Sf3905 = Sf3765*MSbot2(2)**4 + 
     &   12*MBot2*MTop2*Sf77*MSq2Diff(2,1) - Sf3904*MStop2(1)

	Sf4719 = -5*Sf4543*MSbot2(2)**3 + 
     &   20*MTop2*Sf77*MSq2Diff(2,1) - Sf4718*MStop2(1)

	Sf3934 = -(Sf3777*MSbot2(2)**3) + 
     &   12*MBot2*Sf77*MSq2Diff(2,1) + Sf3933*MStop2(1)

	Sf2730 = 3/2.D0*Sf77 - Log(MSbot2(2))*MSbot2(2) - 
     &   Log(MStop2(1))*MStop2(1)

	Sf7584 = 3/2.D0*Sf7581 - MZ2*Log(MZ2) - 2*Log(MStau2(1))*MStau2(1)

	Sf7526 = 3/2.D0*Sf7523 - MZ2*Log(MZ2) - 2*Log(MSbot2(2))*MSbot2(2)

	Sf7507 = 3/2.D0*Sf7504 - MZ2*Log(MZ2) - 2*Log(MSbot2(1))*MSbot2(1)

	Sf8717 = Sf7501 + Sf7521 + 3*Sf8716

	Sf7500 = 3/2.D0*Sf7499 - MZ2*Log(MZ2) - 2*Log(MSneu2(1))*MSneu2(1)

	Sf1600 = ML2*Sf741 - 2*Sf742*MSneu2(1) + Sf743*MStau2(2)

	Sf744 = ML2*Sf741 - 2*Sf742*MSneu2(1) + Sf743*MStau2(1)

	Sf7217 = 3*CA*MTop*Sf6978 - SAB*SB*Sf86*UStop21 + 
     &   3*MTop*reMUE*SA*UStop22

	Sf6969 = 36*C2B*MW2*Sf6966 + 27*MW2*Sf99*MSq2Diff(2,2) + 
     &   22*MSq2Diff(2,2)**3

	Sf3438 = 9*MBot2*Sf69**2 - 9*MTop2*Sf69**2 + Sf3437*MSq2Diff(1,1)

	Sf3770 = Sf3765*MSbot2(1)**4 + 
     &   12*MBot2*MTop2*Sf69*MSq2Diff(1,1) - Sf3769*MStop2(1)

	Sf4624 = -5*Sf4543*MSbot2(1)**3 + 
     &   20*MTop2*Sf69*MSq2Diff(1,1) - Sf4623*MStop2(1)

	Sf3810 = -(Sf3777*MSbot2(1)**3) + 
     &   12*MBot2*Sf69*MSq2Diff(1,1) + Sf3809*MStop2(1)

	Sf2551 = 3/2.D0*Sf69 - Log(MSbot2(1))*MSbot2(1) - 
     &   Log(MStop2(1))*MStop2(1)

	Sf6311 = Sf6309 - 2*Sf6310*MStop2(2) + 5*MStop2(2)**2

	Sf6224 = Sf5848 - 2*Sf5847*MSbot2(2) + 20*MSbot2(2)**2

	Sf5858 = Sf5848 - 2*Sf5847*MSbot2(1) + 20*MSbot2(1)**2

	Sf1449 = 2*Sf545 - Sf191*MStau2(2) - 4*MStau2(2)**2

	Sf546 = 2*Sf545 - Sf191*MStau2(1) - 4*MStau2(1)**2

	Sf3521 = 2*Sf545 - Sf191*MSbot2(2) - 4*MSbot2(2)**2

	Sf3159 = 2*Sf545 - Sf191*MSbot2(1) - 4*MSbot2(1)**2

	Sf5226 = 3*MBot*SA*SB*Sf5224 + 3*CA*CB*MTop*Sf5225 + 
     &   ML*SA*SB*MStau2(2)*UStau21*UStau22

	Sf4878 = 2*MW2*Sf4877 + 5*MSbot2(2)**2 - 5*MStop2(2)**2

	Sf1378 = 3/2.D0*Sf47 - Log(MSneu2(1))*MSneu2(1) - 
     &   Log(MStau2(2))*MStau2(2)

	Sf7196 = 3*CA*MTop*Sf4261 - SAB*SB*Sf86*UStop11 + 
     &   3*MTop*reMUE*SA*UStop12

	Sf4250 = 36*C2B*MW2*Sf4247 + 27*MW2*Sf77*MSq2Diff(2,1) + 
     &   22*MSq2Diff(2,1)**3

	Sf4241 = 36*C2B*MW2*Sf4238 + 27*MW2*Sf69*MSq2Diff(1,1) + 
     &   22*MSq2Diff(1,1)**3

	Sf4165 = Sf4064 - 2*Sf4063*MSbot2(2) + 20*MSbot2(2)**2

	Sf4080 = Sf4064 - 2*Sf4063*MSbot2(1) + 20*MSbot2(1)**2

	Sf6018 = Sf3664 - Sf3665*MStop2(2) + 3*MStop2(2)**2

	Sf3672 = Sf3664 - Sf3665*MStop2(1) + 3*MStop2(1)**2

	Sf461 = 3/2.D0*Sf34 - Log(MSneu2(1))*MSneu2(1) - 
     &   Log(MStau2(1))*MStau2(1)

	Sf5598 = Sf3384 - Sf3385*MStop2(2) + 3*MStop2(2)**2

	Sf3392 = Sf3384 - Sf3385*MStop2(1) + 3*MStop2(1)**2

	Sf2281 = 2*MW2*Sf2280 + 5*MSbot2(2)**2 - 5*MStop2(1)**2

	Sf2218 = 2*MW2*Sf2217 + 5*MSbot2(1)**2 - 5*MStop2(1)**2

	Sf1765 = 3*Sf1738 + Sf556 - MSl2Diff(1,1)

	Sf1773 = CA*ML*Sf1733 - CAB*CB*Ss443*UStau21 - 
     &   ML*reMUE*SA*UStau22

	Sf1785 = Sf1459 + 3*Sf1768 - MSl2Diff(2,1)

	Sf4793 = 3*CA*MBot*Sf144 - CAB*CB*Sf15*USbot21 - 
     &   3*MBot*reMUE*SA*USbot22

	Sf4770 = 3*CA*MBot*Sf133 - CAB*CB*Sf15*USbot11 - 
     &   3*MBot*reMUE*SA*USbot12

	Sf1160 = 2*CB2*SA + 2*CAB*SB + 3*SB*Sf1156

	Sf10516 = MZ2**2*(-MZ2 + MHiggs2(2)) - 2*MStop2(2)**3

	Sf10451 = MZ2**2*(-MZ2 + MHiggs2(2)) - 2*MStop2(1)**3

	Sf9014 = MZ2**2*(-MZ2 + MHiggs2(2)) - 2*MSbot2(2)**3

	Sf8933 = MZ2**2*(-MZ2 + MHiggs2(2)) - 2*MSbot2(1)**3

	Sf10411 = MZ2**2*(-MZ2 + MHiggs2(1)) - 2*MStop2(2)**3

	Sf10356 = MZ2**2*(-MZ2 + MHiggs2(1)) - 2*MStop2(1)**3

	Sf10147 = MZ2**2*(-MZ2 + MHiggs2(1)) - 2*MStau2(2)**3

	Sf10075 = MZ2**2*(-MZ2 + MHiggs2(1)) - 2*MStau2(1)**3

	Sf9212 = MZ2**2*(-MZ2 + MHiggs2(1)) - 2*MSbot2(2)**3

	Sf9135 = MZ2**2*(-MZ2 + MHiggs2(1)) - 2*MSbot2(1)**3

	Sf10383 = 2*MZ2**2 + Sf8528

	Sf9369 = 2*MZ2**2 + Sf7709

	Sf8986 = 2*MZ2**2 + Sf7517

	Sf10507 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 3*Sf8911*MStop2(2)

	Sf10598 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 2*Sf8911*MStop2(2)

	Sf10502 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 3*Sf8911*MStop2(1)

	Sf10595 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 2*Sf8911*MStop2(1)

	Sf9383 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 3*Sf8911*MStau2(2)

	Sf9354 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 2*Sf8911*MStau2(2)

	Sf9379 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 3*Sf8911*MStau2(1)

	Sf9324 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 2*Sf8911*MStau2(1)

	Sf9006 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 3*Sf8911*MSbot2(2)

	Sf8962 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 2*Sf8911*MSbot2(2)

	Sf9002 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 3*Sf8911*MSbot2(1)

	Sf8913 = MZ2*(-5*MZ2 + 3*MHiggs2(2)) + 2*Sf8911*MSbot2(1)

	Sf10026 = MZ2*(-MZ2 + MHiggs2(2)) - Sf9949*MStau2(2)

	Sf9950 = MZ2*(-MZ2 + MHiggs2(2)) - Sf9949*MStau2(1)

	Sf10277 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 3*Sf9113*MStop2(2)

	Sf10374 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 2*Sf9113*MStop2(2)

	Sf10272 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 3*Sf9113*MStop2(1)

	Sf10347 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 2*Sf9113*MStop2(1)

	Sf10138 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 3*Sf9113*MStau2(2)

	Sf10090 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 2*Sf9113*MStau2(2)

	Sf10134 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 3*Sf9113*MStau2(1)

	Sf10061 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 2*Sf9113*MStau2(1)

	Sf9203 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 3*Sf9113*MSbot2(2)

	Sf9152 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 2*Sf9113*MSbot2(2)

	Sf9199 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 3*Sf9113*MSbot2(1)

	Sf9115 = MZ2*(-5*MZ2 + 3*MHiggs2(1)) + 2*Sf9113*MSbot2(1)

	Sf83 = Sf8660 - MZ*Sqrt(Sf82)

	Sf62 = Sf7750 - MZ*Sqrt(Sf61)

	Sf40 = Sf7728 - MZ*Sqrt(Sf39)

	Sf27 = Sf8354 - MZ*Sqrt(Sf26)

	Sf114 = Sf8699 - MZ*Sqrt(Sf113)

	Sf12 = Sf8257 - MZ*Sqrt(Sf11)

	Sf8 = MZ2 - MZ*Sqrt(Sf1) - 2*MSneu2(1)

	Sf1425 = ML2*MW2 + MW2**2*Sf175 + 3*MStau2(2)**2

	Sf514 = ML2*MW2 + MW2**2*Sf175 + 3*MStau2(1)**2

	Sf1212 = MW2*(-8*MW2 + 5*MHiggs2(4)) + 4*Sf173*MStau2(2)

	Sf219 = MW2*(-8*MW2 + 5*MHiggs2(4)) + 4*Sf173*MStau2(1)

	Sf2686 = MW2*(-8*MW2 + 5*MHiggs2(4)) + 3*Sf173*MSbot2(2)

	Sf2494 = MW2*(-8*MW2 + 5*MHiggs2(4)) + 3*Sf173*MSbot2(1)

	Sf2725 = MW2*(-5*MW2 + 4*MHiggs2(4)) + 3*Sf173*MSbot2(2)

	Sf2546 = MW2*(-5*MW2 + 4*MHiggs2(4)) + 3*Sf173*MSbot2(1)

	Sf1213 = MW2*(-4*MW2 + 3*MHiggs2(4)) - Sf173*MStau2(2)

	Sf220 = MW2*(-4*MW2 + 3*MHiggs2(4)) - Sf173*MStau2(1)

	Sf2723 = MW2*(-4*MW2 + 3*MHiggs2(4)) - Sf173*MSbot2(2)

	Sf2544 = MW2*(-4*MW2 + 3*MHiggs2(4)) - Sf173*MSbot2(1)

	Sf4804 = 3*MW2 - Sf3064*UStop22**2 + 
     &   Sf4803*(-1 + UStop22**2)

	Sf4783 = 3*MW2 - Sf3056*UStop22**2 + 
     &   Sf4782*(-1 + UStop22**2)

	Sf7247 = 3*MW2 - Sf556*UStau22**2 + 
     &   Sf1459*(-1 + UStau22**2)

	Sf7212 = 3*MW2 - Sf4781*USbot22**2 + 
     &   Sf4802*(-1 + USbot22**2)

	Sf7189 = 3*MW2 - Sf2564*USbot22**2 + 
     &   Sf2741*(-1 + USbot22**2)

	Sf4876 = 2*MW2 + 7*MSbot2(2) - MStop2(2)

	Sf2279 = 2*MW2 + 7*MSbot2(2) - MStop2(1)

	Sf4842 = 2*MW2 + 7*MSbot2(1) - MStop2(2)

	Sf2216 = 2*MW2 + 7*MSbot2(1) - MStop2(1)

	Sf1192 = MW2 + Sf47

	Sf163 = MW2 + Sf34

	Sf44 = MW2 - Sf47

	Sf31 = MW2 - Sf34

	Sf5147 = MW2 + Sf99

	Sf2739 = MW2 + Sf77

	Sf96 = MW2 - Sf99

	Sf74 = MW2 - Sf77

	Sf5028 = MW2 + Sf91

	Sf2561 = MW2 + Sf69

	Sf88 = MW2 - Sf91

	Sf66 = MW2 - Sf69

	Sf4691 = 30*MBot2*MTop2 + 684*Ss30*MSbot2(2) + 417*C2B*MSbot2(2)**2

	Sf4590 = 30*MBot2*MTop2 + 684*Ss30*MSbot2(1) + 417*C2B*MSbot2(1)**2

	Sf1791 = 15*ML2 + MW2*Sf950 + 28*MStau2(2)

	Sf951 = 15*ML2 + MW2*Sf950 + 28*MStau2(1)

	Sf1251 = 14*ML2 + Sf277*MSneu2(1) - Sf278*MStau2(2)

	Sf281 = 14*ML2 + Sf277*MSneu2(1) - Sf278*MStau2(1)

	Sf1245 = 7*ML2 + Sf175*MSneu2(1) - Sf174*MStau2(2)

	Sf269 = 7*ML2 + Sf175*MSneu2(1) - Sf174*MStau2(1)

	Sf1554 = 7*ML2 + 189*MSneu2(1) + Sf682*MStau2(2)

	Sf683 = 7*ML2 + 189*MSneu2(1) + Sf682*MStau2(1)

	Sf2115 = 5*ML2 + Sf2031*MSneu2(1) + 3*Sf174*MStau2(2)

	Sf2032 = 5*ML2 + Sf2031*MSneu2(1) + 3*Sf174*MStau2(1)

	Sf1576 = 5*ML2 + Sf169*MSneu2(1) - Sf710*MStau2(2)

	Sf711 = 5*ML2 + Sf169*MSneu2(1) - Sf710*MStau2(1)

	Sf1522 = 5*ML2 + MW2*Sf627 + 4*MStau2(2)

	Sf628 = 5*ML2 + MW2*Sf627 + 4*MStau2(1)

	Sf1523 = 5*ML2 - 5*C2B*MW2 + 2*MStau2(2)

	Sf629 = 5*ML2 - 5*C2B*MW2 + 2*MStau2(1)

	Sf1205 = 4*ML2 - Sf207*MHiggs2(4) + 5*MStau2(2)

	Sf208 = 4*ML2 - Sf207*MHiggs2(4) + 5*MStau2(1)

	Sf1583 = 3*ML2 + Sf720*MSneu2(1) + Sf721*MStau2(2)

	Sf722 = 3*ML2 + Sf720*MSneu2(1) + Sf721*MStau2(1)

	Sf1571 = 2*ML2 - Sf703*MSneu2(1) - 2*Sf704*MStau2(2)

	Sf705 = 2*ML2 - Sf703*MSneu2(1) - 2*Sf704*MStau2(1)

	Sf2108 = 2*ML2 + 4*Sf118*MSneu2(1) + Sf2023*MStau2(2)

	Sf2024 = 2*ML2 + 4*Sf118*MSneu2(1) + Sf2023*MStau2(1)

	Sf1199 = 2*ML2 + 2*MW2*Sf181 + 5*MStau2(2)

	Sf182 = 2*ML2 + 2*MW2*Sf181 + 5*MStau2(1)

	Sf1539 = ML2 - Sf661*MSneu2(1) + Sf662*MStau2(2)

	Sf663 = ML2 - Sf661*MSneu2(1) + Sf662*MStau2(1)

	Sf1433 = ML2 + Sf522*MSneu2(1) + Sf523*MStau2(2)

	Sf524 = ML2 + Sf522*MSneu2(1) + Sf523*MStau2(1)

	Sf1259 = ML2 + Sf293*MSneu2(1) + Sf294*MStau2(2)

	Sf295 = ML2 + Sf293*MSneu2(1) + Sf294*MStau2(1)

	Sf1492 = ML2 + Sf181*MSneu2(1) - Sf522*MStau2(2)

	Sf592 = ML2 + Sf181*MSneu2(1) - Sf522*MStau2(1)

	Sf1500 = ML2 + Sf169*MSneu2(1) + 3*Sf174*MStau2(2)

	Sf600 = ML2 + Sf169*MSneu2(1) + 3*Sf174*MStau2(1)

	Sf1127 = 3*ML2 - MW2*Sf1126 + 19*MHiggs2(4)

	Sf1138 = -3*ML2 - MW2*Sf684 + 19*MHiggs2(4)

	Sf1140 = -6*ML2 - 2*MW2*Sf1129 + 7*MHiggs2(4)

	Sf1130 = -3*ML2 - MW2*Sf1129 + 6*MHiggs2(4)

	Sf1134 = 3*ML2 + MW2*Sf1129 + 4*MHiggs2(4)

	Sf1529 = 3*MHiggs2(4) - 5*MSneu2(1) - 7*MStau2(2)

	Sf635 = 3*MHiggs2(4) - 5*MSneu2(1) - 7*MStau2(1)

	Sf6782 = 3*MHiggs2(4) - 5*MSbot2(2) - 7*MStop2(2)

	Sf4424 = 3*MHiggs2(4) - 5*MSbot2(2) - 7*MStop2(1)

	Sf6528 = 3*MHiggs2(4) - 5*MSbot2(1) - 7*MStop2(2)

	Sf4287 = 3*MHiggs2(4) - 5*MSbot2(1) - 7*MStop2(1)

	Sf4483 = MW2 + 3*MHiggs2(4) + 4*MSbot2(2)

	Sf4372 = MW2 + 3*MHiggs2(4) + 4*MSbot2(1)

	Sf4482 = Sf733 - 2*MSbot2(2)

	Sf4371 = Sf733 - 2*MSbot2(1)

	Sf1444 = Sf191 - 4*MStau2(2)

	Sf539 = Sf191 - 4*MStau2(1)

	Sf3525 = Sf191 - 4*MSbot2(2)

	Sf3163 = Sf191 - 4*MSbot2(1)

	Sf1222 = Sf47 + MHiggs2(4)

	Sf231 = Sf34 + MHiggs2(4)

	Sf1219 = -Sf47 + MHiggs2(4)

	Sf228 = -Sf34 + MHiggs2(4)

	Sf1232 = Sf1185 - 3*MSneu2(1)

	Sf250 = Sf155 - 3*MSneu2(1)

	Sf4875 = Sf99 + MHiggs2(4)

	Sf2278 = Sf77 + MHiggs2(4)

	Sf4892 = -Sf99 + MHiggs2(4)

	Sf6213 = -Sf4868 + MHiggs2(4)

	Sf2303 = -Sf77 + MHiggs2(4)

	Sf4152 = -Sf2271 + MHiggs2(4)

	Sf5058 = Sf3740 - MStop2(2)

	Sf2604 = Sf3740 - MStop2(1)

	Sf4841 = Sf91 + MHiggs2(4)

	Sf2215 = Sf69 + MHiggs2(4)

	Sf4858 = -Sf91 + MHiggs2(4)

	Sf5836 = -Sf4839 + MHiggs2(4)

	Sf2252 = -Sf69 + MHiggs2(4)

	Sf4041 = -Sf2207 + MHiggs2(4)

	Sf4932 = Sf3730 - MStop2(2)

	Sf2366 = Sf3730 - MStop2(1)

	Sf1508 = Sf613 + 4*MStau2(2)

	Sf612 = Sf613 + 4*MStau2(1)

	Sf3716 = Sf6120 + MHiggs2(4)

	Sf3444 = Sf5733 + MHiggs2(4)

	Sf6139 = -Sf6728 + MHiggs2(4)

	Sf5756 = -Sf6461 + MHiggs2(4)

	Sf239 = -2*ML2 + 2*C2B*MW2 + MHiggs2(4)

	Sf1107 = -3*ML2 - MW2*Sf1106 + MHiggs2(4)

	Sf10560 = -4*MZ2**2 + 3*MHiggs2(2)**2 + Sf9055*MStop2(2)

	Sf10552 = -4*MZ2**2 + 3*MHiggs2(2)**2 + Sf9055*MStop2(1)

	Sf9066 = -4*MZ2**2 + 3*MHiggs2(2)**2 + Sf9055*MSbot2(2)

	Sf9056 = -4*MZ2**2 + 3*MHiggs2(2)**2 + Sf9055*MSbot2(1)

	Sf10553 = 7*MHiggs2(2) + 4*MStop2(1) - 4*MStop2(3)

	Sf10561 = 7*MHiggs2(2) + 4*MStop2(2) + 4*MStop2(3)

	Sf10564 = 7*MHiggs2(2) - 11*MStop2(1) - 30*MStop2(2)

	Sf10556 = 7*MHiggs2(2) - 30*MStop2(1) - 11*MStop2(2)

	Sf9059 = 7*MHiggs2(2) + 4*MSbot2(1) - 4*MSbot2(3)

	Sf9068 = 7*MHiggs2(2) + 4*MSbot2(2) + 4*MSbot2(3)

	Sf9071 = 7*MHiggs2(2) - 11*MSbot2(1) - 30*MSbot2(2)

	Sf9062 = 7*MHiggs2(2) - 30*MSbot2(1) - 11*MSbot2(2)

	Sf10543 = -2*MZ2 + 4*MHiggs2(2) + 3*MStop2(2)

	Sf10529 = -2*MZ2 + 4*MHiggs2(2) + 3*MStop2(1)

	Sf9617 = -2*MZ2 + 4*MHiggs2(2) + 3*MStau2(2)

	Sf9603 = -2*MZ2 + 4*MHiggs2(2) + 3*MStau2(1)

	Sf9043 = -2*MZ2 + 4*MHiggs2(2) + 3*MSbot2(2)

	Sf9029 = -2*MZ2 + 4*MHiggs2(2) + 3*MSbot2(1)

	Sf10491 = 3*MHiggs2(2) + 2*MStop2(1) - 2*MStop2(3)

	Sf10485 = 3*MHiggs2(2) + 2*MStop2(2) + 2*MStop2(3)

	Sf8989 = 3*MHiggs2(2) + 2*MSbot2(1) - 2*MSbot2(3)

	Sf8980 = 3*MHiggs2(2) + 2*MSbot2(2) + 2*MSbot2(3)

	Sf10540 = Sf8923 - MStop2(2)

	Sf10519 = Sf8923 - MStop2(1)

	Sf9614 = Sf8923 - MStau2(2)

	Sf9594 = Sf8923 - MStau2(1)

	Sf9040 = Sf8923 - MSbot2(2)

	Sf9019 = Sf8923 - MSbot2(1)

	Sf10312 = -4*MZ2**2 + 3*MHiggs2(1)**2 + Sf9186*MStop2(2)

	Sf10303 = -4*MZ2**2 + 3*MHiggs2(1)**2 + Sf9186*MStop2(1)

	Sf10152 = -4*MZ2**2 + 3*MHiggs2(1)**2 + Sf9186*MStau2(2)

	Sf10122 = -4*MZ2**2 + 3*MHiggs2(1)**2 + Sf9186*MStau2(1)

	Sf9217 = -4*MZ2**2 + 3*MHiggs2(1)**2 + Sf9186*MSbot2(2)

	Sf9187 = -4*MZ2**2 + 3*MHiggs2(1)**2 + Sf9186*MSbot2(1)

	Sf10305 = 7*MHiggs2(1) + 4*MStop2(1) - 4*MStop2(3)

	Sf10314 = 7*MHiggs2(1) + 4*MStop2(2) + 4*MStop2(3)

	Sf10317 = 7*MHiggs2(1) - 11*MStop2(1) - 30*MStop2(2)

	Sf10308 = 7*MHiggs2(1) - 30*MStop2(1) - 11*MStop2(2)

	Sf10121 = 7*MHiggs2(1) + 4*MStau2(1) - 4*MStau2(3)

	Sf10150 = 7*MHiggs2(1) + 4*MStau2(2) + 4*MStau2(3)

	Sf10151 = 7*MHiggs2(1) - 11*MStau2(1) - 30*MStau2(2)

	Sf10123 = 7*MHiggs2(1) - 30*MStau2(1) - 11*MStau2(2)

	Sf9185 = 7*MHiggs2(1) + 4*MSbot2(1) - 4*MSbot2(3)

	Sf9215 = 7*MHiggs2(1) + 4*MSbot2(2) + 4*MSbot2(3)

	Sf9216 = 7*MHiggs2(1) - 11*MSbot2(1) - 30*MSbot2(2)

	Sf9188 = 7*MHiggs2(1) - 30*MSbot2(1) - 11*MSbot2(2)

	Sf10401 = -2*MZ2 + 4*MHiggs2(1) + 3*MStop2(2)

	Sf10394 = -2*MZ2 + 4*MHiggs2(1) + 3*MStop2(1)

	Sf10182 = -2*MZ2 + 4*MHiggs2(1) + 3*MStau2(2)

	Sf10168 = -2*MZ2 + 4*MHiggs2(1) + 3*MStau2(1)

	Sf9248 = -2*MZ2 + 4*MHiggs2(1) + 3*MSbot2(2)

	Sf9234 = -2*MZ2 + 4*MHiggs2(1) + 3*MSbot2(1)

	Sf10297 = 3*MHiggs2(1) + 2*MStop2(1) - 2*MStop2(3)

	Sf10287 = 3*MHiggs2(1) + 2*MStop2(2) + 2*MStop2(3)

	Sf10113 = 3*MHiggs2(1) + 2*MStau2(1) - 2*MStau2(3)

	Sf10105 = 3*MHiggs2(1) + 2*MStau2(2) + 2*MStau2(3)

	Sf9176 = 3*MHiggs2(1) + 2*MSbot2(1) - 2*MSbot2(3)

	Sf9170 = 3*MHiggs2(1) + 2*MSbot2(2) + 2*MSbot2(3)

	Sf10368 = Sf9125 - MStop2(2)

	Sf10337 = Sf9125 - MStop2(1)

	Sf10179 = Sf9125 - MStau2(2)

	Sf10158 = Sf9125 - MStau2(1)

	Sf9245 = Sf9125 - MSbot2(2)

	Sf9223 = Sf9125 - MSbot2(1)

	Sf10251 = Sf113 + MHiggs2(1)

	Sf10231 = Sf82 + MHiggs2(1)

	Sf2807 = MBot2**2 + MTop2**2 + C2B**2*MW2**2

	Sf3839 = 191*MBot2 - 191*MTop2 + MW2*Sf3838

	Sf3361 = 152*MBot2 - 152*MTop2 - MW2*Sf3360

	Sf5652 = 148*MBot2 - 148*MTop2 + MW2*Sf5651

	Sf3833 = 104*MBot2 - 104*MTop2 + MW2*Sf3832

	Sf3334 = 100*MBot2 - 100*MTop2 + MW2*Sf3333

	Sf5649 = 89*MBot2 - 89*MTop2 + MW2*Sf5648

	Sf5660 = 87*MBot2 - 87*MTop2 + MW2*Sf5659

	Sf3367 = 84*MBot2 - 84*MTop2 - MW2*Sf3366

	Sf3331 = 84*MBot2 - 84*MTop2 - MW2*Sf3330

	Sf3320 = 84*MBot2 - 84*MTop2 + MW2*Sf3319

	Sf3312 = 84*MBot2 - 84*MTop2 + MW2*Sf3311

	Sf3581 = 75*MBot2 - 75*MTop2 - 7*Sf3231*MSbot2(2)

	Sf3232 = 75*MBot2 - 75*MTop2 - 7*Sf3231*MSbot2(1)

	Sf3364 = 68*MBot2 - 68*MTop2 - MW2*Sf3363

	Sf3309 = 68*MBot2 - 68*MTop2 + MW2*Sf3308

	Sf4143 = 66*MBot2 - 66*MTop2 - Sf4027*MSbot2(2)

	Sf4028 = 66*MBot2 - 66*MTop2 - Sf4027*MSbot2(1)

	Sf3317 = 64*MBot2 - 64*MTop2 + MW2*Sf3316

	Sf3328 = 61*MBot2 - 61*MTop2 + MW2*Sf3327

	Sf3897 = 58*MBot2 - 58*MTop2 - 24*C2B*MSbot2(2)

	Sf3752 = 58*MBot2 - 58*MTop2 - 24*C2B*MSbot2(1)

	Sf5675 = 58*MBot2 - 58*MTop2 - MW2*Sf5674

	Sf3824 = 56*MBot2 - 56*MTop2 - MW2*Sf3823

	Sf2675 = 51*MBot2 - 51*MTop2 - 4*Sf2478*MSbot2(2)

	Sf2479 = 51*MBot2 - 51*MTop2 - 4*Sf2478*MSbot2(1)

	Sf2618 = 42*MBot2 - 42*MTop2 - 3*Sf2392*MSbot2(2)

	Sf2393 = 42*MBot2 - 42*MTop2 - 3*Sf2392*MSbot2(1)

	Sf2679 = 39*MBot2 - 39*MTop2 - 2*Sf2485*MSbot2(2)

	Sf2632 = 39*MBot2 - 39*MTop2 - Sf2356*MSbot2(2)

	Sf3500 = 39*MBot2 - 39*MTop2 - 7*Sf169*MSbot2(2)

	Sf2486 = 39*MBot2 - 39*MTop2 - 2*Sf2485*MSbot2(1)

	Sf2413 = 39*MBot2 - 39*MTop2 - Sf2356*MSbot2(1)

	Sf3124 = 39*MBot2 - 39*MTop2 - 7*Sf169*MSbot2(1)

	Sf5657 = 36*MBot2 - 36*MTop2 + MW2*Sf5656

	Sf4513 = 33*MBot2 - 33*MTop2 + 4*Sf1047*MSbot2(2)

	Sf2634 = 33*MBot2 - 33*MTop2 - 38*MSbot2(2)

	Sf4409 = 33*MBot2 - 33*MTop2 + 4*Sf1047*MSbot2(1)

	Sf2416 = 33*MBot2 - 33*MTop2 - 38*MSbot2(1)

	Sf3344 = 32*MBot2 - 32*MTop2 + MW2*Sf3343

	Sf2656 = 30*MBot2 - 30*MTop2 - Sf2446*MSbot2(2)

	Sf2447 = 30*MBot2 - 30*MTop2 - Sf2446*MSbot2(1)

	Sf3347 = 28*MBot2 - 28*MTop2 - MW2*Sf3346

	Sf3299 = 28*MBot2 - 28*MTop2 + MW2*Sf3298

	Sf2625 = 27*MBot2 - 27*MTop2 + 4*Sf601*MSbot2(2)

	Sf4144 = 27*MBot2 - 27*MTop2 + Sf3270*MSbot2(2)

	Sf3495 = 27*MBot2 - 27*MTop2 - Sf3116*MSbot2(2)

	Sf2403 = 27*MBot2 - 27*MTop2 + 4*Sf601*MSbot2(1)

	Sf4029 = 27*MBot2 - 27*MTop2 + Sf3270*MSbot2(1)

	Sf3117 = 27*MBot2 - 27*MTop2 - Sf3116*MSbot2(1)

	Sf3296 = 26*MBot2 - 26*MTop2 + MW2*Sf3295

	Sf3821 = 25*MBot2 - 25*MTop2 - MW2*Sf3820

	Sf2660 = 24*MBot2 - 24*MTop2 - Sf2452*MSbot2(2)

	Sf2453 = 24*MBot2 - 24*MTop2 - Sf2452*MSbot2(1)

	Sf3492 = 21*MBot2 - 21*MTop2 - Sf3111*MSbot2(2)

	Sf4443 = 21*MBot2 - 21*MTop2 - 4*MSbot2(2)

	Sf4442 = 21*MBot2 - 21*MTop2 - 19*MSbot2(2)

	Sf3112 = 21*MBot2 - 21*MTop2 - Sf3111*MSbot2(1)

	Sf4316 = 21*MBot2 - 21*MTop2 - 4*MSbot2(1)

	Sf4315 = 21*MBot2 - 21*MTop2 - 19*MSbot2(1)

	Sf5436 = 20*MBot2 + Sf2881*MSbot2(2) - 27*Sf2882*MStop2(2)

	Sf3044 = 20*MBot2 + Sf2881*MSbot2(2) - 27*Sf2882*MStop2(1)

	Sf5330 = 20*MBot2 + Sf2881*MSbot2(1) - 27*Sf2882*MStop2(2)

	Sf2917 = 20*MBot2 + Sf2881*MSbot2(1) - 27*Sf2882*MStop2(1)

	Sf5697 = 19*MBot2 - 19*MTop2 - MW2*Sf5696

	Sf5678 = 18*MBot2 - 18*MTop2 - MW2*Sf5677

	Sf3827 = 18*MBot2 - 18*MTop2 - MW2*Sf3826

	Sf3861 = 16*MBot2 - 16*MTop2 - MW2*Sf3860

	Sf3306 = 16*MBot2 - 16*MTop2 + MW2*Sf3305

	Sf4503 = 15*MBot2 - 15*MTop2 - Sf4394*MSbot2(2)

	Sf2666 = 15*MBot2 - 15*MTop2 - Sf2462*MSbot2(2)

	Sf2629 = 15*MBot2 - 15*MTop2 - 58*MSbot2(2)

	Sf4395 = 15*MBot2 - 15*MTop2 - Sf4394*MSbot2(1)

	Sf2463 = 15*MBot2 - 15*MTop2 - Sf2462*MSbot2(1)

	Sf2409 = 15*MBot2 - 15*MTop2 - 58*MSbot2(1)

	Sf5641 = 13*MBot2 - 13*MTop2 + MW2*Sf5640

	Sf4493 = 12*MBot2 - 12*MTop2 + Sf684*MSbot2(2)

	Sf4384 = 12*MBot2 - 12*MTop2 + Sf684*MSbot2(1)

	Sf5646 = 12*MBot2 - 12*MTop2 + MW2*Sf5645

	Sf5636 = 12*MBot2 - 12*MTop2 + MW2*Sf5635

	Sf3850 = 12*MBot2 - 12*MTop2 - MW2*Sf3849

	Sf3844 = 12*MBot2 - 12*MTop2 - MW2*Sf3843

	Sf3325 = 12*MBot2 - 12*MTop2 + MW2*Sf3324

	Sf3289 = 12*MBot2 - 12*MTop2 + MW2*Sf3288

	Sf3847 = 11*MBot2 - 11*MTop2 - MW2*Sf3846

	Sf3358 = 10*MBot2 - 10*MTop2 + MW2*Sf3357

	Sf3686 = 9*MBot2 - 9*MTop2 - Sf858*MSbot2(2)

	Sf3706 = 9*MBot2 - 9*MTop2 - 2*Sf747*MSbot2(2)

	Sf6721 = 9*MBot2 - 9*MTop2 - 4*Sf602*MSbot2(2)

	Sf3496 = 9*MBot2 - 9*MTop2 - 2*Sf504*MSbot2(2)

	Sf4453 = 9*MBot2 - 9*MTop2 + Sf4329*MSbot2(2)

	Sf4194 = 9*MBot2 - 9*MTop2 - 2*Sf4112*MSbot2(2)

	Sf4138 = 9*MBot2 - 9*MTop2 - Sf4019*MSbot2(2)

	Sf3705 = 9*MBot2 - 9*MTop2 + Sf3429*MSbot2(2)

	Sf4512 = 9*MBot2 - 9*MTop2 + Sf3250*MSbot2(2)

	Sf3685 = 9*MBot2 - 9*MTop2 - Sf2508*MSbot2(2)

	Sf4174 = 9*MBot2 - 9*MTop2 - Sf2462*MSbot2(2)

	Sf4177 = 9*MBot2 - 9*MTop2 - Sf2356*MSbot2(2)

	Sf2624 = 9*MBot2 - 9*MTop2 - Sf1903*MSbot2(2)

	Sf4459 = 9*MBot2 - 9*MTop2 + 49*MSbot2(2)

	Sf3507 = 9*MBot2 - 9*MTop2 + 10*MSbot2(2)

	Sf4460 = 9*MBot2 - 9*MTop2 + 4*MSbot2(2)

	Sf4190 = 9*MBot2 - 9*MTop2 - 20*MSbot2(2)

	Sf3506 = 9*MBot2 - 9*MTop2 - 44*MSbot2(2)

	Sf4139 = 9*MBot2 - 9*MTop2 - 46*MSbot2(2)

	Sf3407 = 9*MBot2 - 9*MTop2 - Sf858*MSbot2(1)

	Sf3431 = 9*MBot2 - 9*MTop2 - 2*Sf747*MSbot2(1)

	Sf6453 = 9*MBot2 - 9*MTop2 - 4*Sf602*MSbot2(1)

	Sf3118 = 9*MBot2 - 9*MTop2 - 2*Sf504*MSbot2(1)

	Sf4330 = 9*MBot2 - 9*MTop2 + Sf4329*MSbot2(1)

	Sf4113 = 9*MBot2 - 9*MTop2 - 2*Sf4112*MSbot2(1)

	Sf4020 = 9*MBot2 - 9*MTop2 - Sf4019*MSbot2(1)

	Sf3430 = 9*MBot2 - 9*MTop2 + Sf3429*MSbot2(1)

	Sf4408 = 9*MBot2 - 9*MTop2 + Sf3250*MSbot2(1)

	Sf3406 = 9*MBot2 - 9*MTop2 - Sf2508*MSbot2(1)

	Sf4089 = 9*MBot2 - 9*MTop2 - Sf2462*MSbot2(1)

	Sf4092 = 9*MBot2 - 9*MTop2 - Sf2356*MSbot2(1)

	Sf2402 = 9*MBot2 - 9*MTop2 - Sf1903*MSbot2(1)

	Sf4338 = 9*MBot2 - 9*MTop2 + 49*MSbot2(1)

	Sf3141 = 9*MBot2 - 9*MTop2 + 10*MSbot2(1)

	Sf4339 = 9*MBot2 - 9*MTop2 + 4*MSbot2(1)

	Sf4108 = 9*MBot2 - 9*MTop2 - 20*MSbot2(1)

	Sf3140 = 9*MBot2 - 9*MTop2 - 44*MSbot2(1)

	Sf4021 = 9*MBot2 - 9*MTop2 - 46*MSbot2(1)

	Sf4071 = 9*MBot2 - 9*MTop2 + MW2*Sf4070

	Sf4077 = 8*MBot2 - 8*MTop2 - MW2*Sf4076

	Sf3341 = 8*MBot2 - 8*MTop2 - MW2*Sf3340

	Sf3338 = 8*MBot2 - 8*MTop2 - MW2*Sf294

	Sf4451 = 7*MBot2 - 7*MTop2 + Sf4325*MSbot2(2)

	Sf2669 = 7*MBot2 - 7*MTop2 - Sf2467*MSbot2(2)

	Sf6710 = 7*MBot2 - 7*MTop2 - Sf186*MSbot2(2)

	Sf3493 = 7*MBot2 - 7*MTop2 + 2*Sf175*MSbot2(2)

	Sf2710 = 7*MBot2 - 7*MTop2 + 17*Sf174*MSbot2(2)

	Sf2706 = 7*MBot2 - 7*MTop2 + 5*Sf174*MSbot2(2)

	Sf4499 = 7*MBot2 - 7*MTop2 - Sf169*MSbot2(2)

	Sf4326 = 7*MBot2 - 7*MTop2 + Sf4325*MSbot2(1)

	Sf2468 = 7*MBot2 - 7*MTop2 - Sf2467*MSbot2(1)

	Sf6440 = 7*MBot2 - 7*MTop2 - Sf186*MSbot2(1)

	Sf3113 = 7*MBot2 - 7*MTop2 + 2*Sf175*MSbot2(1)

	Sf2526 = 7*MBot2 - 7*MTop2 + 17*Sf174*MSbot2(1)

	Sf2520 = 7*MBot2 - 7*MTop2 + 5*Sf174*MSbot2(1)

	Sf4390 = 7*MBot2 - 7*MTop2 - Sf169*MSbot2(1)

	Sf3303 = 7*MBot2 - 7*MTop2 + MW2*Sf413

	Sf5093 = 6*MBot2 - 6*MTop2 + 5*MSq2Diff(2,2)

	Sf2653 = 6*MBot2 - 6*MTop2 + 5*MSq2Diff(2,1)

	Sf4967 = 6*MBot2 - 6*MTop2 + 5*MSq2Diff(1,2)

	Sf2440 = 6*MBot2 - 6*MTop2 + 5*MSq2Diff(1,1)

	Sf6112 = 6*MBot2 - 6*MTop2 + Sf859*MSbot2(2)

	Sf4464 = 6*MBot2 - 6*MTop2 + Sf4344*MSbot2(2)

	Sf6711 = 6*MBot2 - 6*MTop2 + Sf4106*MSbot2(2)

	Sf4146 = 6*MBot2 - 6*MTop2 + Sf4032*MSbot2(2)

	Sf4502 = 6*MBot2 - 6*MTop2 - Sf3084*MSbot2(2)

	Sf2691 = 6*MBot2 - 6*MTop2 + Sf2073*MSbot2(2)

	Sf6707 = 6*MBot2 - 6*MTop2 + Sf207*MSbot2(2)

	Sf4173 = 6*MBot2 - 6*MTop2 - Sf1116*MSbot2(2)

	Sf4489 = 6*MBot2 - 6*MTop2 + 5*MSbot2(2)

	Sf4132 = 6*MBot2 - 6*MTop2 - 4*MSbot2(2)

	Sf4130 = 6*MBot2 - 6*MTop2 - 74*MSbot2(2)

	Sf5723 = 6*MBot2 - 6*MTop2 + Sf859*MSbot2(1)

	Sf4345 = 6*MBot2 - 6*MTop2 + Sf4344*MSbot2(1)

	Sf6441 = 6*MBot2 - 6*MTop2 + Sf4106*MSbot2(1)

	Sf4033 = 6*MBot2 - 6*MTop2 + Sf4032*MSbot2(1)

	Sf4393 = 6*MBot2 - 6*MTop2 - Sf3084*MSbot2(1)

	Sf2500 = 6*MBot2 - 6*MTop2 + Sf2073*MSbot2(1)

	Sf6437 = 6*MBot2 - 6*MTop2 + Sf207*MSbot2(1)

	Sf4088 = 6*MBot2 - 6*MTop2 - Sf1116*MSbot2(1)

	Sf4380 = 6*MBot2 - 6*MTop2 + 5*MSbot2(1)

	Sf4012 = 6*MBot2 - 6*MTop2 - 4*MSbot2(1)

	Sf4010 = 6*MBot2 - 6*MTop2 - 74*MSbot2(1)

	Sf5662 = 6*MBot2 - 6*MTop2 - MW2*Sf4394

	Sf2236 = 6*MBot2 - 6*MTop2 + MW2*Sf2235

	Sf2659 = 5*MBot2 - 5*MTop2 + Sf2450*MSbot2(2)

	Sf2697 = 5*MBot2 - 5*MTop2 + 4*Sf243*MSbot2(2)

	Sf2451 = 5*MBot2 - 5*MTop2 + Sf2450*MSbot2(1)

	Sf2507 = 5*MBot2 - 5*MTop2 + 4*Sf243*MSbot2(1)

	Sf2693 = 4*MBot2 - 4*MTop2 + Sf438*MSbot2(2)

	Sf2502 = 4*MBot2 - 4*MTop2 + Sf438*MSbot2(1)

	Sf3132 = 4*MBot2 - 4*MTop2 - MW2*Sf569

	Sf3254 = 4*MBot2 - 4*MTop2 + MW2*Sf3253

	Sf3207 = 4*MBot2 - 4*MTop2 + MW2*Sf195

	Sf5700 = 4*MBot2 - 4*MTop2 - MW2*Sf168

	Sf5880 = 3*MBot2 - 3*MTop2 + Sf603*MStop2(2)

	Sf5617 = 3*MBot2 - 3*MTop2 + 10*MStop2(2)

	Sf5099 = 3*MBot2 - 3*MTop2 + MSq2Diff(2,2)

	Sf6829 = 3*MBot2 - 3*MTop2 - MSq2Diff(2,2)

	Sf5933 = 3*MBot2 - 3*MTop2 - 5*MSq2Diff(2,2)

	Sf6245 = 3*MBot2 - 3*MTop2 - 11*MSq2Diff(2,2)

	Sf6786 = 3*MBot2 - 3*MTop2 - 20*MSq2Diff(2,2)

	Sf2663 = 3*MBot2 - 3*MTop2 + MSq2Diff(2,1)

	Sf4496 = 3*MBot2 - 3*MTop2 - MSq2Diff(2,1)

	Sf3490 = 3*MBot2 - 3*MTop2 - 5*MSq2Diff(2,1)

	Sf4197 = 3*MBot2 - 3*MTop2 - 11*MSq2Diff(2,1)

	Sf4430 = 3*MBot2 - 3*MTop2 - 20*MSq2Diff(2,1)

	Sf4973 = 3*MBot2 - 3*MTop2 + MSq2Diff(1,2)

	Sf6574 = 3*MBot2 - 3*MTop2 - MSq2Diff(1,2)

	Sf5500 = 3*MBot2 - 3*MTop2 - 5*MSq2Diff(1,2)

	Sf5870 = 3*MBot2 - 3*MTop2 - 11*MSq2Diff(1,2)

	Sf6527 = 3*MBot2 - 3*MTop2 - 20*MSq2Diff(1,2)

	Sf2458 = 3*MBot2 - 3*MTop2 + MSq2Diff(1,1)

	Sf4387 = 3*MBot2 - 3*MTop2 - MSq2Diff(1,1)

	Sf3107 = 3*MBot2 - 3*MTop2 - 5*MSq2Diff(1,1)

	Sf4117 = 3*MBot2 - 3*MTop2 - 11*MSq2Diff(1,1)

	Sf4301 = 3*MBot2 - 3*MTop2 - 20*MSq2Diff(1,1)

	Sf4181 = 3*MBot2 - 3*MTop2 - Sf765*MSbot2(2)

	Sf3689 = 3*MBot2 - 3*MTop2 - Sf712*MSbot2(2)

	Sf4447 = 3*MBot2 - 3*MTop2 + 4*Sf601*MSbot2(2)

	Sf6128 = 3*MBot2 - 3*MTop2 - Sf588*MSbot2(2)

	Sf4193 = 3*MBot2 - 3*MTop2 - Sf532*MSbot2(2)

	Sf4189 = 3*MBot2 - 3*MTop2 - Sf4106*MSbot2(2)

	Sf4492 = 3*MBot2 - 3*MTop2 + Sf3436*MSbot2(2)

	Sf2713 = 3*MBot2 - 3*MTop2 + 4*Sf2530*MSbot2(2)

	Sf4131 = 3*MBot2 - 3*MTop2 - Sf2508*MSbot2(2)

	Sf2674 = 3*MBot2 - 3*MTop2 - Sf2476*MSbot2(2)

	Sf2655 = 3*MBot2 - 3*MTop2 - Sf2444*MSbot2(2)

	Sf2639 = 3*MBot2 - 3*MTop2 - Sf2422*MSbot2(2)

	Sf4498 = 3*MBot2 - 3*MTop2 + Sf1091*MSbot2(2)

	Sf4438 = 3*MBot2 - 3*MTop2 + 86*MSbot2(2)

	Sf3701 = 3*MBot2 - 3*MTop2 + 32*MSbot2(2)

	Sf3682 = -3*MTop2 + Sf3929

	Sf6698 = 3*MBot2 - 3*MTop2 - MSbot2(2)

	Sf4184 = 3*MBot2 - 3*MTop2 - 13*MSbot2(2)

	Sf4432 = 3*MBot2 - 3*MTop2 - 17*MSbot2(2)

	Sf4452 = 3*MBot2 - 3*MTop2 - 53*MSbot2(2)

	Sf4097 = 3*MBot2 - 3*MTop2 - Sf765*MSbot2(1)

	Sf3411 = 3*MBot2 - 3*MTop2 - Sf712*MSbot2(1)

	Sf4320 = 3*MBot2 - 3*MTop2 + 4*Sf601*MSbot2(1)

	Sf5743 = 3*MBot2 - 3*MTop2 - Sf588*MSbot2(1)

	Sf4111 = 3*MBot2 - 3*MTop2 - Sf532*MSbot2(1)

	Sf4107 = 3*MBot2 - 3*MTop2 - Sf4106*MSbot2(1)

	Sf4383 = 3*MBot2 - 3*MTop2 + Sf3436*MSbot2(1)

	Sf2531 = 3*MBot2 - 3*MTop2 + 4*Sf2530*MSbot2(1)

	Sf4011 = 3*MBot2 - 3*MTop2 - Sf2508*MSbot2(1)

	Sf2477 = 3*MBot2 - 3*MTop2 - Sf2476*MSbot2(1)

	Sf2445 = 3*MBot2 - 3*MTop2 - Sf2444*MSbot2(1)

	Sf2423 = 3*MBot2 - 3*MTop2 - Sf2422*MSbot2(1)

	Sf4389 = 3*MBot2 - 3*MTop2 + Sf1091*MSbot2(1)

	Sf4311 = 3*MBot2 - 3*MTop2 + 86*MSbot2(1)

	Sf3425 = 3*MBot2 - 3*MTop2 + 32*MSbot2(1)

	Sf3403 = -3*MTop2 + Sf3805

	Sf6426 = 3*MBot2 - 3*MTop2 - MSbot2(1)

	Sf4100 = 3*MBot2 - 3*MTop2 - 13*MSbot2(1)

	Sf4304 = 3*MBot2 - 3*MTop2 - 17*MSbot2(1)

	Sf4327 = 3*MBot2 - 3*MTop2 - 53*MSbot2(1)

	Sf5728 = 3*MBot2 - 3*MTop2 - MW2*Sf2371

	Sf2231 = 3*MBot2 - 3*MTop2 - MW2*Sf2230

	Sf2223 = 3*MBot2 - 3*MTop2 + MW2*Sf1060

	Sf3816 = 3*MBot2 - 3*MTop2 + MW2*Sf1047

	Sf6701 = 2*MBot2 - 2*MTop2 - Sf499*MSbot2(2)

	Sf4471 = 2*MBot2 - 2*MTop2 - Sf4356*MSbot2(2)

	Sf6429 = 2*MBot2 - 2*MTop2 - Sf499*MSbot2(1)

	Sf4357 = 2*MBot2 - 2*MTop2 - Sf4356*MSbot2(1)

	Sf3293 = 2*MBot2 - 2*MTop2 + MW2*Sf184

	Sf5878 = MBot2 - MTop2 + 5*MStop2(2)

	Sf6836 = MBot2 - MTop2 - MSq2Diff(2,2)

	Sf5048 = MBot2 - MTop2 - 2*MSq2Diff(2,2)

	Sf6025 = MBot2 - MTop2 - 3*MSq2Diff(2,2)

	Sf6214 = MBot2 - MTop2 - 4*MSq2Diff(2,2)

	Sf4507 = MBot2 - MTop2 - MSq2Diff(2,1)

	Sf2594 = MBot2 - MTop2 - 2*MSq2Diff(2,1)

	Sf3679 = MBot2 - MTop2 - 3*MSq2Diff(2,1)

	Sf4153 = MBot2 - MTop2 - 4*MSq2Diff(2,1)

	Sf6581 = MBot2 - MTop2 - MSq2Diff(1,2)

	Sf4922 = MBot2 - MTop2 - 2*MSq2Diff(1,2)

	Sf5605 = MBot2 - MTop2 - 3*MSq2Diff(1,2)

	Sf5837 = MBot2 - MTop2 - 4*MSq2Diff(1,2)

	Sf4400 = MBot2 - MTop2 - MSq2Diff(1,1)

	Sf2355 = MBot2 - MTop2 - 2*MSq2Diff(1,1)

	Sf3399 = MBot2 - MTop2 - 3*MSq2Diff(1,1)

	Sf4042 = MBot2 - MTop2 - 4*MSq2Diff(1,1)

	Sf4437 = MBot2 - MTop2 + Sf4309*MSbot2(2)

	Sf3510 = MBot2 - MTop2 - Sf3144*MSbot2(2)

	Sf6702 = MBot2 - MTop2 - Sf300*MSbot2(2)

	Sf4446 = MBot2 - MTop2 - Sf2828*MSbot2(2)

	Sf2698 = MBot2 - MTop2 + Sf2508*MSbot2(2)

	Sf3681 = MBot2 - MTop2 - Sf1050*MSbot2(2)

	Sf3723 = MBot2 - MTop2 - 6*C2B*MSbot2(2)

	Sf3702 = MBot2 - MTop2 + 2*MSbot2(2)

	Sf3696 = MBot2 - MTop2 - MSbot2(2)

	Sf6132 = MBot2 - MTop2 - 4*MSbot2(2)

	Sf6697 = MBot2 - MTop2 - 5*MSbot2(2)

	Sf4454 = MBot2 - MTop2 - 6*MSbot2(2)

	Sf4488 = MBot2 - MTop2 - 11*MSbot2(2)

	Sf4310 = MBot2 - MTop2 + Sf4309*MSbot2(1)

	Sf3145 = MBot2 - MTop2 - Sf3144*MSbot2(1)

	Sf6430 = MBot2 - MTop2 - Sf300*MSbot2(1)

	Sf4319 = MBot2 - MTop2 - Sf2828*MSbot2(1)

	Sf2509 = MBot2 - MTop2 + Sf2508*MSbot2(1)

	Sf3402 = MBot2 - MTop2 - Sf1050*MSbot2(1)

	Sf3452 = MBot2 - MTop2 - 6*C2B*MSbot2(1)

	Sf3426 = MBot2 - MTop2 + 2*MSbot2(1)

	Sf3418 = MBot2 - MTop2 - MSbot2(1)

	Sf5748 = MBot2 - MTop2 - 4*MSbot2(1)

	Sf6425 = MBot2 - MTop2 - 5*MSbot2(1)

	Sf4331 = MBot2 - MTop2 - 6*MSbot2(1)

	Sf4379 = MBot2 - MTop2 - 11*MSbot2(1)

	Sf3210 = MBot2 - MTop2 + MW2*Sf651

	Sf3206 = MBot2 - MTop2 + MW2*Sf410

	Sf3213 = MBot2 - MTop2 + MW2*Sf3192

	Sf4067 = MBot2 - MTop2 - MW2*Sf243

	Sf5638 = MBot2 - MTop2 + MW2*Sf175

	Sf2552 = -MTop2 + Sf2864

	Sf1154 = CA*CAB*CB - CA2*SB2 + SA*Sf1153

	Sf1894 = 12 + Pi**2 - 6*Sf1893*Log(MSneu2(1))

	Sf8561 = 72 + 3*UStau22**2

	Sf8637 = 396 + UStau22**2*(72 + 21*UStau22**2)

	Sf8828 = 245 + UStau22**2*(-54 + 54*UStau22**2)

	Sf8846 = 28 + UStau22**2*(-6 + 15*UStau22**2)

	Sf8622 = 96 + 3*UStau22**2

	Sf8584 = 38 + 3*UStau22**2

	Sf7472 = Sf5396 + MSbot2(1) + 4*MStop2(1)

	Sf7464 = Sf99 - MSbot2(3)*USbot22**2 - 
     &   MStop2(3)*UStop22**2

	Sf10348 = 2*MZ2**2*MHiggs2(1) + 
     &   MZ2*(-3*MZ2**2 - Sf9142*MStop2(1)) - Sf10347*MStop2(2) + 
     &   Sf9113*(MStop2(1)**2 + MStop2(2)**2)

	Sf10377 = 2*MZ2**2*MHiggs2(1) - Sf10374*MStop2(1) + 
     &   MZ2*(-3*MZ2**2 - Sf9142*MStop2(2)) + 
     &   Sf9113*(MStop2(1)**2 + MStop2(2)**2)

	Sf1875 = -(Sf1874*MSneu2(1)) + 2*Sf1873*MSneu2(1)**2 - 
     &   2*MStau2(2)**3 + MHiggs2(4)*
     &    (MHiggs2(4)*(2*MHiggs2(4) - 11*MStau2(2)) + 
     &      11*MStau2(2)**2)

	Sf1124 = -(Sf1123*MSneu2(1)) + 2*Sf1122*MSneu2(1)**2 - 
     &   2*MStau2(1)**3 + MHiggs2(4)*
     &    (MHiggs2(4)*(2*MHiggs2(4) - 11*MStau2(1)) + 
     &      11*MStau2(1)**2)

	Sf9979 = 8*(-3*MStau2(1)*MStau2(2)**2 + MStau2(2)**3 + 
     &      MStau2(1)**2*(2*MStau2(2) + MStau2(3))) + 
     &   15*MZ2*MHiggs2(2)*MStau2(4)

	Sf10006 = Sf10005*MStau2(1) + 
     &   831*MZ2*MHiggs2(2)*MStau2(2) + 64*MStau2(2)**3 + 
     &   MStau2(1)**2*(128*MStau2(2) + 64*MStau2(3))

	Sf9905 = -64*MStau2(3)**3 + 831*MZ2*MHiggs2(2)*MStau2(4)

	Sf9871 = -8*MStau2(3)**3 + 15*MZ2*MHiggs2(2)*MStau2(4)

	Sf5165 = Sf99*(Ss509 - 3*Log(MStop2(2)/MSbot2(2))**2) - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2) - 
     &   3*Sf5163*Sf5164*MSq2Diff(2,2) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf5015 = Sf91*(Ss509 - 3*Log(MSbot2(1)/MStop2(2))**2) - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   3*Sf5013*Sf5014*MSq2Diff(1,2) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf2758 = Sf77*(Ss509 - 3*Log(MStop2(1)/MSbot2(2))**2) - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2) - 
     &   3*Sf2756*Sf2757*MSq2Diff(2,1) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1)

	Sf2582 = Sf69*(Ss509 - 3*Log(MStop2(1)/MSbot2(1))**2) - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   3*Sf2580*Sf2581*MSq2Diff(1,1) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1)

	Sf1520 = Sf47*(Ss509 - 3*Log(MSneu2(1)/MStau2(2))**2) + 
     &   3*Sf1518*Sf1519*MSl2Diff(2,1) - 
     &   12*Sf227*Log(MSneu2(1))*MSneu2(1) - 
     &   12*Ss954*Log(MStau2(2))*MStau2(2)

	Sf625 = Sf34*(Ss509 - 3*Log(MSneu2(1)/MStau2(1))**2) + 
     &   3*Sf623*Sf624*MSl2Diff(1,1) - 
     &   12*Sf227*Log(MSneu2(1))*MSneu2(1) - 
     &   12*Ss933*Log(MStau2(1))*MStau2(1)

	Sf7420 = 2*MW2 + MSq2Diff(2,2) - MSbot2(3)*USbot22**2

	Sf1511 = MW2**3 - Sf613*MStau2(2)**2 + 4*MStau2(2)**3 + 
     &   MHiggs2(4)*(3*MW2**2 + 2*Sf180*MStau2(2))

	Sf616 = MW2**3 - Sf613*MStau2(1)**2 + 4*MStau2(1)**3 + 
     &   MHiggs2(4)*(3*MW2**2 + 2*Sf180*MStau2(1))

	Sf3713 = MW2**3 - Sf613*MSbot2(2)**2 + 4*MSbot2(2)**3 + 
     &   MHiggs2(4)*(3*MW2**2 + 2*Sf180*MSbot2(2))

	Sf3441 = MW2**3 - Sf613*MSbot2(1)**2 + 4*MSbot2(1)**3 + 
     &   MHiggs2(4)*(3*MW2**2 + 2*Sf180*MSbot2(1))

	Sf254 = -2*Sf253*MSneu2(1) - 22*MSneu2(1)**2 - 
     &   6*MStau2(1)**2 + MHiggs2(4)*(4*MHiggs2(4) + 2*MStau2(1))

	Sf1236 = -2*Sf1235*MSneu2(1) - 22*MSneu2(1)**2 - 
     &   6*MStau2(2)**2 + MHiggs2(4)*(4*MHiggs2(4) + 2*MStau2(2))

	Sf6337 = 3*MBot2 - 3*MTop2 + (1 - 3*C2B)*MSbot2(2) + 
     &   Sf4104*MStop2(2)

	Sf4188 = 3*MBot2 - 3*MTop2 + (1 - 3*C2B)*MSbot2(2) + 
     &   Sf4104*MStop2(1)

	Sf6314 = 3*MBot2 - 3*MTop2 + (1 - 3*C2B)*MSbot2(1) + 
     &   Sf4104*MStop2(2)

	Sf4105 = 3*MBot2 - 3*MTop2 + (1 - 3*C2B)*MSbot2(1) + 
     &   Sf4104*MStop2(1)

	Sf7431 = 2 + 30*UStop22**2

	Sf7428 = 2 + 6*UStop22**2

	Sf8757 = MZ2*(MZ2**3*(-(MZ2*Sf8754) + Sf8756) + 
     &      8*Sf8722*UStop11*UStop12**2) + 
     &   Sf7199*Sf7756*Sf86**2*UStop11*UStop21**2

	Sf7346 = 2 + 6*UStau22**2

	Sf7343 = 2 + 2*UStau22**2

	Sf7781 = MZ2*(MZ2**3*(-(MZ2*Sf7779) + Sf7780) + 
     &      4*Sf7757*UStau11*UStau12**2) + 
     &   Sf1755*Sf7756*Ss443**2*UStau11*UStau21**2

	Sf8842 = -2 - 6*USbot22**2

	Sf7380 = 2 + 6*USbot22**2

	Sf9608 = 3/2.D0*Sf9607 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   Log(MStau2(1))*MStau2(1) - Log(MStau2(2))*MStau2(2)

	Sf948 = CA2*ML*SB*Sf945 + SA*Sf947 - 
     &   CA*CB*Sf127*Ss443*UStau11 + CB*ML*reMUE*UStau12

	Sf9239 = 3/2.D0*Sf9238 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   Log(MSbot2(1))*MSbot2(1) - Log(MSbot2(2))*MSbot2(2)

	Sf9034 = 3/2.D0*Sf9033 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   Log(MSbot2(1))*MSbot2(1) - Log(MSbot2(2))*MSbot2(2)

	Sf8526 = 3/2.D0*Sf8521 - MZ2*Log(MZ2) - 
     &   Log(MStop2(1))*MStop2(1) - Log(MStop2(2))*MStop2(2)

	Sf7707 = 3/2.D0*Sf7702 - MZ2*Log(MZ2) - 
     &   Log(MStau2(1))*MStau2(1) - Log(MStau2(2))*MStau2(2)

	Sf7515 = 3/2.D0*Sf7510 - MZ2*Log(MZ2) - 
     &   Log(MSbot2(1))*MSbot2(1) - Log(MSbot2(2))*MSbot2(2)

	Sf7399 = 2*MW2 + MSq2Diff(2,1) - MSbot2(3)*USbot22**2

	Sf7386 = 2*MW2 - MSq2Diff(2,2) - MStop2(3)*UStop22**2

	Sf7367 = 2*MW2 - MSq2Diff(1,2) - MStop2(3)*UStop22**2

	Sf7311 = Sf1428 - MSneu2(1) - MStau2(3)*UStau22**2

	Sf5156 = 3/2.D0*Sf5147 - MW2*Log(MW2) - 
     &   Log(MSbot2(2))*MSbot2(2) - Log(MStop2(2))*MStop2(2)

	Sf5037 = 3/2.D0*Sf5028 - MW2*Log(MW2) - 
     &   Log(MSbot2(1))*MSbot2(1) - Log(MStop2(2))*MStop2(2)

	Sf4882 = 3/2.D0*Sf4875 - Log(MHiggs2(4))*MHiggs2(4) - 
     &   Log(MSbot2(2))*MSbot2(2) - Log(MStop2(2))*MStop2(2)

	Sf4848 = 3/2.D0*Sf4841 - Log(MHiggs2(4))*MHiggs2(4) - 
     &   Log(MSbot2(1))*MSbot2(1) - Log(MStop2(2))*MStop2(2)

	Sf2749 = 3/2.D0*Sf2739 - MW2*Log(MW2) - 
     &   Log(MSbot2(2))*MSbot2(2) - Log(MStop2(1))*MStop2(1)

	Sf2572 = 3/2.D0*Sf2561 - MW2*Log(MW2) - 
     &   Log(MSbot2(1))*MSbot2(1) - Log(MStop2(1))*MStop2(1)

	Sf454 = 3/2.D0*Sf231 - Log(MHiggs2(4))*MHiggs2(4) - 
     &   Log(MSneu2(1))*MSneu2(1) - Log(MStau2(1))*MStau2(1)

	Sf2285 = 3/2.D0*Sf2278 - Log(MHiggs2(4))*MHiggs2(4) - 
     &   Log(MSbot2(2))*MSbot2(2) - Log(MStop2(1))*MStop2(1)

	Sf2222 = 3/2.D0*Sf2215 - Log(MHiggs2(4))*MHiggs2(4) - 
     &   Log(MSbot2(1))*MSbot2(1) - Log(MStop2(1))*MStop2(1)

	Sf1736 = CA2*ML*SB*Sf1733 + SA*Sf1735 - 
     &   CA*CB*Sf127*Ss443*UStau21 + CB*ML*reMUE*UStau22

	Sf469 = 3/2.D0*Sf163 - MW2*Log(MW2) - 
     &   Log(MSneu2(1))*MSneu2(1) - Log(MStau2(1))*MStau2(1)

	Sf1371 = 3/2.D0*Sf1222 - Log(MHiggs2(4))*MHiggs2(4) - 
     &   Log(MSneu2(1))*MSneu2(1) - Log(MStau2(2))*MStau2(2)

	Sf1386 = 3/2.D0*Sf1192 - MW2*Log(MW2) - 
     &   Log(MSneu2(1))*MSneu2(1) - Log(MStau2(2))*MStau2(2)

	Sf10534 = 3/2.D0*Sf10533 - Log(MHiggs2(2))*MHiggs2(2) - 
     &   Log(MStop2(1))*MStop2(1) - Log(MStop2(2))*MStop2(2)

	Sf10399 = 3/2.D0*Sf10398 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   Log(MStop2(1))*MStop2(1) - Log(MStop2(2))*MStop2(2)

	Sf10173 = 3/2.D0*Sf10172 - Log(MHiggs2(1))*MHiggs2(1) - 
     &   Log(MStau2(1))*MStau2(1) - Log(MStau2(2))*MStau2(2)

	Sf5853 = 24*MBot2 - 24*MTop2 + MW2*Sf4072 + 244*MStop2(2)

	Sf4073 = 24*MBot2 - 24*MTop2 + MW2*Sf4072 + 244*MStop2(1)

	Sf6238 = 9*MBot2 - 9*MTop2 - Sf5738*MSbot2(2) + Sf601*MStop2(2)

	Sf6125 = 9*MBot2 - 9*MTop2 - Sf5738*MSbot2(2) + Sf601*MStop2(1)

	Sf5100 = 9*MBot2 - 9*MTop2 - 39*MSbot2(2) - Sf2459*MStop2(2)

	Sf2664 = 9*MBot2 - 9*MTop2 - 39*MSbot2(2) - Sf2459*MStop2(1)

	Sf5868 = 9*MBot2 - 9*MTop2 - Sf5738*MSbot2(1) + Sf601*MStop2(2)

	Sf5739 = 9*MBot2 - 9*MTop2 - Sf5738*MSbot2(1) + Sf601*MStop2(1)

	Sf4974 = 9*MBot2 - 9*MTop2 - 39*MSbot2(1) - Sf2459*MStop2(2)

	Sf2460 = 9*MBot2 - 9*MTop2 - 39*MSbot2(1) - Sf2459*MStop2(1)

	Sf5855 = 9*MBot2 - 9*MTop2 + MW2*Sf4054 + 62*MStop2(2)

	Sf4075 = 9*MBot2 - 9*MTop2 + MW2*Sf4054 + 62*MStop2(1)

	Sf5573 = 9*MBot2 - 9*MTop2 + MW2*Sf3250 - 10*MStop2(2)

	Sf3251 = 9*MBot2 - 9*MTop2 + MW2*Sf3250 - 10*MStop2(1)

	Sf2590 = 6*MBot2 - 6*MTop2 + MW2*Sf2350 - 6*MSbot2(2)

	Sf2351 = 6*MBot2 - 6*MTop2 + MW2*Sf2350 - 6*MSbot2(1)

	Sf2589 = 6*MBot2 - 6*MTop2 + MW2*Sf2348 - 7*MSbot2(2)

	Sf2349 = 6*MBot2 - 6*MTop2 + MW2*Sf2348 - 7*MSbot2(1)

	Sf6240 = 3*MBot2 - 3*MTop2 - Sf5741*MSbot2(2) + 3*Sf175*MStop2(2)

	Sf6127 = 3*MBot2 - 3*MTop2 - Sf5741*MSbot2(2) + 3*Sf175*MStop2(1)

	Sf5080 = 3*MBot2 - 3*MTop2 + Sf2425*MSbot2(2) + 3*Sf293*MStop2(2)

	Sf2640 = 3*MBot2 - 3*MTop2 + Sf2425*MSbot2(2) + 3*Sf293*MStop2(1)

	Sf6326 = 3*MBot2 - 3*MTop2 + 7*MSbot2(2) + 3*Sf712*MStop2(2)

	Sf4172 = 3*MBot2 - 3*MTop2 + 7*MSbot2(2) + 3*Sf712*MStop2(1)

	Sf6198 = 3*MBot2 - 3*MTop2 - 29*MSbot2(2) + 25*MStop2(2)

	Sf4129 = 3*MBot2 - 3*MTop2 - 29*MSbot2(2) + 25*MStop2(1)

	Sf5871 = 3*MBot2 - 3*MTop2 - Sf5741*MSbot2(1) + 3*Sf175*MStop2(2)

	Sf5742 = 3*MBot2 - 3*MTop2 - Sf5741*MSbot2(1) + 3*Sf175*MStop2(1)

	Sf4954 = 3*MBot2 - 3*MTop2 + Sf2425*MSbot2(1) + 3*Sf293*MStop2(2)

	Sf2426 = 3*MBot2 - 3*MTop2 + Sf2425*MSbot2(1) + 3*Sf293*MStop2(1)

	Sf6301 = 3*MBot2 - 3*MTop2 + 7*MSbot2(1) + 3*Sf712*MStop2(2)

	Sf4087 = 3*MBot2 - 3*MTop2 + 7*MSbot2(1) + 3*Sf712*MStop2(1)

	Sf5820 = 3*MBot2 - 3*MTop2 - 29*MSbot2(1) + 25*MStop2(2)

	Sf4009 = 3*MBot2 - 3*MTop2 - 29*MSbot2(1) + 25*MStop2(1)

	Sf6118 = 3*MBot2 - 3*MTop2 - MW2*Sf603 - 14*MSbot2(2)

	Sf5731 = 3*MBot2 - 3*MTop2 - MW2*Sf603 - 14*MSbot2(1)

	Sf5509 = 3*MBot2 - 3*MTop2 - MW2*Sf2383 - 2*MStop2(2)

	Sf3130 = 3*MBot2 - 3*MTop2 - MW2*Sf2383 - 2*MStop2(1)

	Sf3558 = 2*MBot2 - 2*MTop2 + MW2*Sf3184 + 33*MSbot2(2)

	Sf3202 = 2*MBot2 - 2*MTop2 + MW2*Sf3184 + 33*MSbot2(1)

	Sf6803 = MBot2 - MTop2 + Sf4336*MSbot2(2) - Sf2033*MStop2(2)

	Sf4458 = MBot2 - MTop2 + Sf4336*MSbot2(2) - Sf2033*MStop2(1)

	Sf6548 = MBot2 - MTop2 + Sf4336*MSbot2(1) - Sf2033*MStop2(2)

	Sf4337 = MBot2 - MTop2 + Sf4336*MSbot2(1) - Sf2033*MStop2(1)

	Sf8856 = Sf8855*MSbot2(1) - 
     &   96*SW2**2*(7*Sf8534 + 50*MStau2(1)) + 
     &   3*(-4087*MStau2(1) + 
     &      64*SW2**2*(7*Sf8534 - 7*Sf8536 + 25*MStau2(1)) - 
     &      2977*MStau2(2))*UStau22**2 + 
     &   (-672*(Sf8534 - 2*Sf8536)*SW2**2 + 
     &      9*(787*MStau2(1) + 899*MStau2(2)))*UStau22**4 - 
     &   8*(30*MSneu2(1) + 69*MStau2(2) + 232*MStop2(1) + 
     &      52*MStop2(2) + MSbot2(2)*(13 + 45*USbot22**2) - 
     &      180*(MStop2(1) - MStop2(2))*UStop22**2)

	Sf8558 = 6*Sf7314*MSneu2(1) - 
     &   3*MStau2(2)*UStau22**2 + 3*MStau2(2)*UStau22**4 - 
     &   2*Sf8534*SW2**2*(-1 + UStau22**2)**2

	Sf8869 = 48656 + 
     &   USbot22**2*(1533 + 3*Sf8868 - 678*USbot22**2) + 
     &   UStau22**2*(-9252 + 15174*UStau22**2)

	Sf8822 = 12221 + UStau22**2*(-2454 + 3342*UStau22**2)

	Sf8614 = 939 + UStau22**2*(66 + 6*UStau22**2)

	Sf8554 = 60 + UStau22**2*(6 - 21*UStau22**2)

	Sf8589 = 48 + UStau22**2*(9 + 3*UStau22**2)

	Sf8837 = 847 + UStau22**2*(-180 + 207*UStau22**2)

	Sf8721 = 3*Sf20**2*Sf22**2*Sf7497*Sf7509**2 + 
     &   9*Sf55**2*Sf57**2*Sf7497*Sf7701**2 + 
     &   3*Sf107**2*Sf109**2*Sf7497*Sf8520**2 + 
     &   MZ2*Sf7497*(-3*Sf7503*MSbot2(1) - 3*Sf7522*MSbot2(2) - 
     &      9*MSneu2(1) - 9*Sf7580*MStau2(1) - 
     &      9*Sf7714*MStau2(2) - 3*Sf7991*MStop2(1) - 
     &      3*Sf8533*MStop2(2)) + 
     &   CW2*MZ2*(16*Sf7501*Ss78*SW2*MSbot2(1) + 
     &      16*Sf7521*Ss78*SW2*MSbot2(2) + 
     &      48*Sf7579*Ss78*SW2*MStau2(1) + 
     &      48*Sf7713*Ss78*SW2*MStau2(2) + 
     &      32*Sf7990*Ss78*SW2*MStop2(1) + 
     &      32*Sf8532*Ss78*SW2*MStop2(2))

	Sf7480 = MSbot2(1) + 2*Sf7478*MStop2(1) + 
     &   2*Sf7479*MStop2(2) + MSbot2(3)*(-1 + 3*USbot22**2)

	Sf9214 = Sf9212*MHiggs2(1) - Sf9211*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf9213 - 2*Ss916*MSbot2(1) + 
     &      MSbot2(1)**2) - Sf9133*MSbot2(2) + 
     &   MSbot2(2)**2*(Sf9134 + MSbot2(2)**2)

	Sf9140 = Sf9135*MHiggs2(1) - Sf9133*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf9134 + MSbot2(1)**2) - Sf9138*MSbot2(2) + 
     &   MSbot2(2)**2*(Sf9139 - 2*Ss898*MSbot2(2) + MSbot2(2)**2)

	Sf9016 = Sf9014*MHiggs2(2) - Sf9013*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf9015 - 2*Ss909*MSbot2(1) + 
     &      MSbot2(1)**2) - Sf8931*MSbot2(2) + 
     &   MSbot2(2)**2*(Sf8932 + MSbot2(2)**2)

	Sf8938 = Sf8933*MHiggs2(2) - Sf8931*MSbot2(1) + 
     &   MSbot2(1)**2*(Sf8932 + MSbot2(1)**2) - Sf8936*MSbot2(2) + 
     &   MSbot2(2)**2*(Sf8937 - 2*Ss888*MSbot2(2) + MSbot2(2)**2)

	Sf10518 = Sf10516*MHiggs2(2) - Sf10515*MStop2(1) + 
     &   MStop2(1)**2*(Sf10517 - 2*Ss991*MStop2(1) + 
     &      MStop2(1)**2) - Sf8931*MStop2(2) + 
     &   MStop2(2)**2*(Sf8932 + MStop2(2)**2)

	tmp11 = MStop2(2)**2 + MStop2(1)*(5*MStop2(1) + 10*MStop2(2))

	tmp13 = MZ2**2 - MZ2*(3*MStop2(1) + MStop2(2)) + MStop2(3)**2

	Sf10455 = -(MZ2**2*tmp11) + tmp13*MHiggs2(2)**2 + 
     &   MZ2**3*(3*MStop2(1) + MStop2(2)) + MStop2(3)**4 + 
     &   MHiggs2(2)*(MZ2*(-MZ2**2 + tmp11) - 
     &      2*MStop2(3)**2*MStop2(4))

	tmp12 = 5*MStop2(2)**2 + MStop2(1)*(MStop2(1) + 10*MStop2(2))

	Sf10413 = -(MZ2**2*tmp12) + 
     &   MZ2**3*(MStop2(1) + 3*MStop2(2)) + MStop2(3)**4 + 
     &   MHiggs2(1)**2*(MZ2**2 - MZ2*(MStop2(1) + 3*MStop2(2)) + 
     &      MStop2(3)**2) + 
     &   MHiggs2(1)*(MZ2*(-MZ2**2 + tmp12) - 
     &      2*MStop2(3)**2*MStop2(4))

	Sf10360 = -(MZ2**2*tmp11) + tmp13*MHiggs2(1)**2 + 
     &   MZ2**3*(3*MStop2(1) + MStop2(2)) + MStop2(3)**4 + 
     &   MHiggs2(1)*(MZ2*(-MZ2**2 + tmp11) - 
     &      2*MStop2(3)**2*MStop2(4))

	tmp10 = 5*MStau2(2)**2 + MStau2(1)*(MStau2(1) + 10*MStau2(2))

	Sf10149 = -(MZ2**2*tmp10) + 
     &   MZ2**3*(MStau2(1) + 3*MStau2(2)) + MStau2(3)**4 + 
     &   MHiggs2(1)**2*(MZ2**2 - MZ2*(MStau2(1) + 3*MStau2(2)) + 
     &      MStau2(3)**2) + 
     &   MHiggs2(1)*(MZ2*(-MZ2**2 + tmp10) - 
     &      2*MStau2(3)**2*MStau2(4))

	tmp9 = MStau2(2)**2 + MStau2(1)*(5*MStau2(1) + 10*MStau2(2))

	Sf10079 = -(MZ2**2*tmp9) + 
     &   MZ2**3*(3*MStau2(1) + MStau2(2)) + MStau2(3)**4 + 
     &   MHiggs2(1)**2*(MZ2**2 - MZ2*(3*MStau2(1) + MStau2(2)) + 
     &      MStau2(3)**2) + 
     &   MHiggs2(1)*(MZ2*(-MZ2**2 + tmp9) - 
     &      2*MStau2(3)**2*MStau2(4))

	Sf7227 = Sf4803 + 3*Sf7212 + MSq2Diff(2,2) + 
     &   (Sf4782 - Sf4803 - MSbot2(3))*USbot22**2

	Sf7208 = Sf3064 + 3*Sf7189 + MSq2Diff(2,1) + 
     &   (Sf3056 - Sf3064 - MSbot2(3))*USbot22**2

	Sf4807 = Sf4802 + 3*Sf4804 - MSq2Diff(2,2) + 
     &   (Sf2741 - Sf4802 - MStop2(3))*UStop22**2

	Sf4786 = Sf4781 + 3*Sf4783 - MSq2Diff(1,2) + 
     &   (Sf2564 - Sf4781 - MStop2(3))*UStop22**2

	Sf9993 = Sf9991*MStau2(1)**2 - Sf9990*MStau2(1)**3 + 
     &   Sf9988*MStau2(1)**4 - 2*Sf9987*MStau2(1)**5 + 
     &   MZ2*(Sf9992*MStau2(1)*MStau2(2) - 2*Sf9989*MStau2(2)**3)

	Sf9835 = MZ2*(MStau2(1)*(6*MStau2(1) - 68*MStau2(2)) + 
     &      6*MStau2(2)**2) - 
     &   2*(-7*MStau2(2) - MStau2(3))*MStau2(3)**2 + 
     &   7*MZ2**2*MStau2(4)

	Sf9814 = -4*MZ2*(3*MStau2(2)**2 + 
     &      MStau2(1)*(3*MStau2(1) + 68*MStau2(2))) + 
     &   2*(5*MStau2(1) + 32*MStau2(2))*MStau2(3)**2 + 
     &   37*MZ2**2*MStau2(4)

	Sf10416 = 12*Ss974*Log(MStop2(1))*MStop2(1) + 
     &   12*Ss993*Log(MStop2(2))*MStop2(2) - 
     &   3*Sf10415*Sf8814*MStop2(3) - Ss509*MStop2(4) + 
     &   3*Log(MStop2(1)/MStop2(2))**2*MStop2(4)

	Sf9841 = 12*Ss933*Log(MStau2(1))*MStau2(1) + 
     &   12*Ss954*Log(MStau2(2))*MStau2(2) - 
     &   3*Sf8784*Sf9840*MStau2(3) - Ss509*MStau2(4) + 
     &   3*Log(MStau2(1)/MStau2(2))**2*MStau2(4)

	Sf7268 = MW2 + MSneu2(1) + 
     &   MSbot2(2)*(-3 + 3*USbot21**2 + 3*USbot22**2)

	Sf9078 = 12*Ss890*Log(MSbot2(1))*MSbot2(1) + 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2) - 
     &   3*Sf8766*Sf9077*MSbot2(3) - Ss509*MSbot2(4) + 
     &   3*Log(MSbot2(1)/MSbot2(2))**2*MSbot2(4)

	Sf10250 = 1/12.D0*(Sf8678*Ss509) - 
     &   1/2.D0*(MZ2*Log(MZ2/MStop2(2))**2) - MZ2*Ss550*Log(MZ2) - 
     &   2*Ss993*Log(MStop2(2))*MStop2(2) + 
     &   TF(MZ2,MStop2(2),MStop2(2))

	Sf10230 = 1/12.D0*(Sf7992*Ss509) - 
     &   1/2.D0*(MZ2*Log(MZ2/MStop2(1))**2) - MZ2*Ss550*Log(MZ2) - 
     &   2*Ss974*Log(MStop2(1))*MStop2(1) + 
     &   TF(MZ2,MStop2(1),MStop2(1))

	Sf9307 = 1/12.D0*(Sf7734*Ss509) - 
     &   1/2.D0*(MZ2*Log(MZ2/MStau2(2))**2) - MZ2*Ss550*Log(MZ2) - 
     &   2*Ss954*Log(MStau2(2))*MStau2(2) + 
     &   TF(MZ2,MStau2(2),MStau2(2))

	Sf9295 = 1/12.D0*(Sf7581*Ss509) - 
     &   1/2.D0*(MZ2*Log(MZ2/MStau2(1))**2) - MZ2*Ss550*Log(MZ2) - 
     &   2*Ss933*Log(MStau2(1))*MStau2(1) + 
     &   TF(MZ2,MStau2(1),MStau2(1))

	Sf8890 = 1/12.D0*(Sf7523*Ss509) - 
     &   1/2.D0*(MZ2*Log(MZ2/MSbot2(2))**2) - MZ2*Ss550*Log(MZ2) - 
     &   2*Ss911*Log(MSbot2(2))*MSbot2(2) + 
     &   TF(MZ2,MSbot2(2),MSbot2(2))

	Sf8878 = 1/12.D0*(Sf7504*Ss509) - 
     &   1/2.D0*(MZ2*Log(MZ2/MSbot2(1))**2) - MZ2*Ss550*Log(MZ2) - 
     &   2*Ss890*Log(MSbot2(1))*MSbot2(1) + 
     &   TF(MZ2,MSbot2(1),MSbot2(1))

	Sf4834 = Sf5022 + MHiggs2(4)**2 - 2*Sf4833*MSbot2(1) - 
     &   Sf4832*MStop2(2)

	Sf10524 = Sf8528 + MHiggs2(2)**2 - Sf8925*MStop2(1) - 
     &   2*Sf10442*MStop2(2)

	Sf10538 = -Sf10497 + MHiggs2(2)**2 - Sf8925*MStop2(2)

	Sf9611 = -Sf9374 + MHiggs2(2)**2 - Sf8925*MStau2(2)

	Sf9598 = Sf7709 + MHiggs2(2)**2 - Sf8925*MStau2(1) - 
     &   2*Sf9329*MStau2(2)

	Sf9037 = -Sf8997 + MHiggs2(2)**2 - Sf8925*MSbot2(2)

	Sf9023 = Sf7517 + MHiggs2(2)**2 - Sf8925*MSbot2(1) - 
     &   2*Sf8918*MSbot2(2)

	Sf10341 = Sf8528 + MHiggs2(1)**2 - Sf9127*MStop2(1) - 
     &   2*Sf10294*MStop2(2)

	Sf10361 = -Sf10388 + MHiggs2(1)**2 - Sf9127*MStop2(2)

	Sf10163 = Sf7709 + MHiggs2(1)**2 - Sf9127*MStau2(1) - 
     &   2*Sf10066*MStau2(2)

	Sf10177 = -Sf10129 + MHiggs2(1)**2 - Sf9127*MStau2(2)

	Sf9243 = -Sf9194 + MHiggs2(1)**2 - Sf9127*MSbot2(2)

	Sf9228 = Sf7517 + MHiggs2(1)**2 - Sf9127*MSbot2(1) - 
     &   2*Sf9120*MSbot2(2)

	Sf5575 = Of1 + 15*MStop2(2)

	Sf3255 = Of1 + 15*MStop2(1)

	Sf5571 = Of2 + 27*MStop2(2)

	Sf3248 = Of2 + 27*MStop2(1)

	Sf3564 = -7*MTop2 - 17*MW2 - 7*C2B*MW2 + Sf2615

	Sf3211 = -7*MTop2 - 17*MW2 - 7*C2B*MW2 + Sf2388

	Sf4918 = Of3 + 8*MStop2(2)

	Sf5511 = 3*MBot2 - 3*MTop2 - 3*C2B*MW2 + Sf5167

	Sf2347 = Of3 + 8*MStop2(1)

	Sf3133 = 3*MBot2 - 3*MTop2 - 3*C2B*MW2 + Sf2584

	Sf6116 = Of3 + 17*MSbot2(2)

	Sf5729 = Of3 + 17*MSbot2(1)

	Sf1152 = 6 - CA*CAB*CB - CB2*SA2 - CAB*SA*SB + CA2*SB2

	Sf8826 = 48*Sf8534*SW2**2 + Sf8825*MSbot2(1) + 
     &   352*SW2**2*MStau2(1) + 
     &   (3058*MStau2(1) - 
     &      32*SW2**2*(3*Sf8534 - 3*Sf8536 + 11*MStau2(1)) + 
     &      1883*MStau2(2))*UStau22**2 + 
     &   3*(16*Sf8534*SW2**2 - 32*Sf8536*SW2**2 - 545*MStau2(1) - 
     &      569*MStau2(2))*UStau22**4 + 
     &   4*(4*MSneu2(1) + 49*MStau2(2) + 124*MStop2(1) + 
     &      40*MStop2(2) + MSbot2(2)*(-5 + 36*USbot22**2) - 
     &      84*(MStop2(1) - MStop2(2))*UStop22**2)

	Sf7782 = 2*MW2**3*MZ2**3*Sf7772 - 
     &   12*MW2**4*MZ2**2*Sf7774 + 32*MW2**5*MZ2*Sf7776 - 
     &   64*MW2**6*Sf7778 + MW2**2*Sf7781 + 
     &   MZ2**5*(MW2*Sf7770 - 9*MZ2*MSl2Diff(1,1))*UStau11

	Sf7284 = 3*MSbot2(2)*USbot21**2 + 
     &   3*MSbot2(1)*USbot22**2 + MStau2(2)*UStau21**2 + 
     &   MStau2(1)*UStau22**2 + 3*MStop2(2)*UStop21**2 + 
     &   3*MStop2(1)*UStop22**2

	Sf6207 = 4*MW2**3 - 27*Ss215*Ss218*MSbot2(2)**2 + 
     &   Sf4026*MSbot2(2)**3 - Sf4143*MSbot2(2)*MStop2(2) - 
     &   Sf4144*MStop2(2)**2 - Sf4030*MStop2(2)**3

	Sf4145 = 4*MW2**3 - 27*Ss215*Ss218*MSbot2(2)**2 + 
     &   Sf4026*MSbot2(2)**3 - Sf4143*MSbot2(2)*MStop2(1) - 
     &   Sf4144*MStop2(1)**2 - Sf4030*MStop2(1)**3

	Sf5829 = 4*MW2**3 - 27*Ss215*Ss218*MSbot2(1)**2 + 
     &   Sf4026*MSbot2(1)**3 - Sf4028*MSbot2(1)*MStop2(2) - 
     &   Sf4029*MStop2(2)**2 - Sf4030*MStop2(2)**3

	Sf4031 = 4*MW2**3 - 27*Ss215*Ss218*MSbot2(1)**2 + 
     &   Sf4026*MSbot2(1)**3 - Sf4028*MSbot2(1)*MStop2(1) - 
     &   Sf4029*MStop2(1)**2 - Sf4030*MStop2(1)**3

	Sf8850 = Sf8849*MSbot2(1) - 42*MSneu2(1) - 
     &   80*SW2**2*(Sf8534 + 6*MStau2(1)) + 
     &   126*MSbot2(2)*(-1 + USbot22**2) + 
     &   2*(80*SW2**2*(Sf8534 - Sf8536 + 3*MStau2(1)) - 
     &      3*(67*MStau2(1) + 53*MStau2(2)))*UStau22**2 - 
     &   40*(2*Sf8534*SW2**2 - 4*Sf8536*SW2**2 - 6*MStau2(1) - 
     &      9*MStau2(2))*UStau22**4 - 
     &   42*(MStau2(2) + 3*(MStop2(2) - MStop2(3)*UStop22**2))

	Sf8834 = -2*MSbot2(4) - 6*MStau2(2) - 8*MStop2(4) - 
     &   3*(20*MStau2(1) + 9*MStau2(2))*UStau22**2 + 
     &   27*MStau2(4)*UStau22**4

	Sf7281 = 21*MSbot2(2) + 7*MSneu2(1) + 7*MStau2(2) + 
     &   21*MStop2(2) - 21*MSbot2(3)*USbot22**2 - 
     &   7*MStau2(3)*UStau22**2 - 21*MStop2(3)*UStop22**2

	Sf9831 = MZ2**5*(4*MStau2(2)**2 + 
     &      MStau2(1)*(9*MStau2(1) + 25*MStau2(2))) - 
     &   MZ2**4*(12*MStau2(2)**3 + 
     &      MStau2(1)*(164*MStau2(2)**2 + 
     &         MStau2(1)*(55*MStau2(1) + 111*MStau2(2)))) + 
     &   MZ2**3*(12*MStau2(2)**4 + 
     &      MStau2(1)*(351*MStau2(2)**3 + 
     &         MStau2(1)*(205*MStau2(2)**2 + 
     &            MStau2(1)*(95*MStau2(1) + 249*MStau2(2))))) - 
     &   MZ2**2*(4*MStau2(2)**5 + 
     &      MStau2(1)*(270*MStau2(2)**4 + 
     &         MStau2(1)*(153*MStau2(2)**3 + 
     &            MStau2(1)*
     &             (-179*MStau2(2)**2 + 
     &               MStau2(1)*(9*MStau2(1) + 351*MStau2(2)))))) - 
     &   2*MZ2*MStau2(1)*(-33*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(30*MStau2(1) - 34*MStau2(2)) - 
     &         3*MStau2(2)**2))*MStau2(3)**2 + 
     &   4*MStau2(1)*(MStau2(1)*(5*MStau2(1) - 8*MStau2(2)) - 
     &      2*MStau2(2)**2)*MStau2(3)**4

	Sf8469 = Sf7436 + 3*MSbot2(2)*USbot21**2 + 
     &   3*MSbot2(1)*USbot22**2 + 3*MStop2(2)*UStop21**2 + 
     &   3*MStop2(1)*UStop22**2

	Sf9852 = Sf9607*Ss509 - 
     &   3*Sf9851*Log(MHiggs2(2)/MStau2(1))**2 - 
     &   3*Sf9850*Log(MHiggs2(2)/MStau2(2))**2 + 
     &   3*Sf9849*Log(MStau2(1)/MStau2(2))**2 - 
     &   12*Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   12*Ss933*Log(MStau2(1))*MStau2(1) - 
     &   12*Ss954*Log(MStau2(2))*MStau2(2)

	Sf9270 = Sf9238*Ss509 - 
     &   3*Sf9269*Log(MHiggs2(1)/MSbot2(1))**2 - 
     &   3*Sf9268*Log(MHiggs2(1)/MSbot2(2))**2 + 
     &   3*Sf9267*Log(MSbot2(1)/MSbot2(2))**2 - 
     &   12*Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2)

	Sf9090 = Sf9033*Ss509 - 
     &   3*Sf9089*Log(MHiggs2(2)/MSbot2(1))**2 - 
     &   3*Sf9088*Log(MHiggs2(2)/MSbot2(2))**2 + 
     &   3*Sf9087*Log(MSbot2(1)/MSbot2(2))**2 - 
     &   12*Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2)

	Sf10327 = Sf8521*Ss509 - 12*MZ2*Ss550*Log(MZ2) - 
     &   3*Sf7199*Log(MZ2/MStop2(1))**2 - 
     &   3*Sf7219*Log(MZ2/MStop2(2))**2 + 
     &   3*Sf104*Log(MStop2(1)/MStop2(2))**2 - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf9865 = Sf7702*Ss509 - 12*MZ2*Ss550*Log(MZ2) - 
     &   3*Sf1755*Log(MZ2/MStau2(1))**2 - 
     &   3*Sf1775*Log(MZ2/MStau2(2))**2 + 
     &   3*Sf52*Log(MStau2(1)/MStau2(2))**2 - 
     &   12*Ss933*Log(MStau2(1))*MStau2(1) - 
     &   12*Ss954*Log(MStau2(2))*MStau2(2)

	Sf9103 = Sf7510*Ss509 - 12*MZ2*Ss550*Log(MZ2) - 
     &   3*Sf4773*Log(MZ2/MSbot2(1))**2 - 
     &   3*Sf4795*Log(MZ2/MSbot2(2))**2 + 
     &   3*Sf17*Log(MSbot2(1)/MSbot2(2))**2 - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2)

	Sf5453 = Sf5147*Ss509 - 12*MW2*Ss530*Log(MW2) - 
     &   3*Sf4803*Log(MW2/MSbot2(2))**2 - 
     &   3*Sf4802*Log(MW2/MStop2(2))**2 + 
     &   3*Sf96*Log(MStop2(2)/MSbot2(2))**2 - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf5446 = Sf5028*Ss509 - 12*MW2*Ss530*Log(MW2) - 
     &   3*Sf4782*Log(MW2/MSbot2(1))**2 - 
     &   3*Sf4781*Log(MW2/MStop2(2))**2 + 
     &   3*Sf88*Log(MSbot2(1)/MStop2(2))**2 - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf4895 = Sf4875*Ss509 - 
     &   3*Sf4894*Log(MHiggs2(4)/MSbot2(2))**2 - 
     &   3*Sf4893*Log(MHiggs2(4)/MStop2(2))**2 + 
     &   3*Sf4892*Log(MStop2(2)/MSbot2(2))**2 - 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4) - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf4861 = Sf4841*Ss509 - 
     &   3*Sf4860*Log(MHiggs2(4)/MSbot2(1))**2 - 
     &   3*Sf4859*Log(MHiggs2(4)/MStop2(2))**2 + 
     &   3*Sf4858*Log(MSbot2(1)/MStop2(2))**2 - 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4) - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf3065 = Sf2739*Ss509 - 12*MW2*Ss530*Log(MW2) - 
     &   3*Sf3064*Log(MW2/MSbot2(2))**2 - 
     &   3*Sf2741*Log(MW2/MStop2(1))**2 + 
     &   3*Sf74*Log(MStop2(1)/MSbot2(2))**2 - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1)

	Sf3057 = Sf2561*Ss509 - 12*MW2*Ss530*Log(MW2) - 
     &   3*Sf3056*Log(MW2/MSbot2(1))**2 - 
     &   3*Sf2564*Log(MW2/MStop2(1))**2 + 
     &   3*Sf66*Log(MStop2(1)/MSbot2(1))**2 - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1)

	Sf232 = Sf231*Ss509 - 
     &   3*Sf230*Log(MHiggs2(4)/MSneu2(1))**2 - 
     &   3*Sf229*Log(MHiggs2(4)/MStau2(1))**2 + 
     &   3*Sf228*Log(MSneu2(1)/MStau2(1))**2 - 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4) - 
     &   12*Sf227*Log(MSneu2(1))*MSneu2(1) - 
     &   12*Ss933*Log(MStau2(1))*MStau2(1)

	Sf2306 = Sf2278*Ss509 - 
     &   3*Sf2305*Log(MHiggs2(4)/MSbot2(2))**2 - 
     &   3*Sf2304*Log(MHiggs2(4)/MStop2(1))**2 + 
     &   3*Sf2303*Log(MStop2(1)/MSbot2(2))**2 - 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4) - 
     &   12*Ss911*Log(MSbot2(2))*MSbot2(2) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1)

	Sf2255 = Sf2215*Ss509 - 
     &   3*Sf2254*Log(MHiggs2(4)/MSbot2(1))**2 - 
     &   3*Sf2253*Log(MHiggs2(4)/MStop2(1))**2 + 
     &   3*Sf2252*Log(MStop2(1)/MSbot2(1))**2 - 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4) - 
     &   12*Ss890*Log(MSbot2(1))*MSbot2(1) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1)

	Sf557 = Sf163*Ss509 - 12*MW2*Ss530*Log(MW2) - 
     &   3*Sf556*Log(MW2/MSneu2(1))**2 - 
     &   3*Sf555*Log(MW2/MStau2(1))**2 + 
     &   3*Sf31*Log(MSneu2(1)/MStau2(1))**2 - 
     &   12*Sf227*Log(MSneu2(1))*MSneu2(1) - 
     &   12*Ss933*Log(MStau2(1))*MStau2(1)

	Sf1223 = Sf1222*Ss509 - 
     &   3*Sf1221*Log(MHiggs2(4)/MSneu2(1))**2 - 
     &   3*Sf1220*Log(MHiggs2(4)/MStau2(2))**2 + 
     &   3*Sf1219*Log(MSneu2(1)/MStau2(2))**2 - 
     &   12*Ss524*Log(MHiggs2(4))*MHiggs2(4) - 
     &   12*Sf227*Log(MSneu2(1))*MSneu2(1) - 
     &   12*Ss954*Log(MStau2(2))*MStau2(2)

	Sf1460 = Sf1192*Ss509 - 12*MW2*Ss530*Log(MW2) - 
     &   3*Sf1459*Log(MW2/MSneu2(1))**2 - 
     &   3*Sf1458*Log(MW2/MStau2(2))**2 + 
     &   3*Sf44*Log(MSneu2(1)/MStau2(2))**2 - 
     &   12*Sf227*Log(MSneu2(1))*MSneu2(1) - 
     &   12*Ss954*Log(MStau2(2))*MStau2(2)

	Sf10577 = Sf10533*Ss509 - 
     &   3*Sf10576*Log(MHiggs2(2)/MStop2(1))**2 - 
     &   3*Sf10575*Log(MHiggs2(2)/MStop2(2))**2 + 
     &   3*Sf10574*Log(MStop2(1)/MStop2(2))**2 - 
     &   12*Ss548*Log(MHiggs2(2))*MHiggs2(2) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf10427 = Sf10398*Ss509 - 
     &   3*Sf10426*Log(MHiggs2(1)/MStop2(1))**2 - 
     &   3*Sf10425*Log(MHiggs2(1)/MStop2(2))**2 + 
     &   3*Sf10424*Log(MStop2(1)/MStop2(2))**2 - 
     &   12*Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   12*Ss974*Log(MStop2(1))*MStop2(1) - 
     &   12*Ss993*Log(MStop2(2))*MStop2(2)

	Sf10202 = Sf10172*Ss509 - 
     &   3*Sf10201*Log(MHiggs2(1)/MStau2(1))**2 - 
     &   3*Sf10200*Log(MHiggs2(1)/MStau2(2))**2 + 
     &   3*Sf10199*Log(MStau2(1)/MStau2(2))**2 - 
     &   12*Ss546*Log(MHiggs2(1))*MHiggs2(1) - 
     &   12*Ss933*Log(MStau2(1))*MStau2(1) - 
     &   12*Ss954*Log(MStau2(2))*MStau2(2)

	Sf8843 = 2*Sf8840*MSbot2(1) + 2*Sf8842*MSbot2(2) - 
     &   36*MStau2(2) - 64*MStop2(1) - 40*MStop2(2) - 
     &   3*(141*MStau2(1) + 73*MStau2(2))*UStau22**2 + 
     &   207*(MStau2(1) + MStau2(2))*UStau22**4 + 
     &   24*(MStop2(1) - MStop2(2))*UStop22**2

	Sf9704 = MZ2**4*(3*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(4*MStau2(1) - 6*MStau2(2)) + 
     &         17*MStau2(2)**2)) - 
     &   MZ2**3*(3*MStau2(2)**4 + 
     &      MStau2(1)*(36*MStau2(2)**3 + 
     &         MStau2(1)*(MStau2(1)*
     &             (17*MStau2(1) - 24*MStau2(2)) + 16*MStau2(2)**2)
     &         )) + MZ2**2*
     &    (MStau2(2)**5 + MStau2(1)*
     &       (27*MStau2(2)**4 + 
     &         MStau2(1)*(3*MStau2(2)**3 + 
     &            MStau2(1)*
     &             (MStau2(1)*(18*MStau2(1) - 24*MStau2(2)) + 
     &               7*MStau2(2)**2)))) - 
     &   MZ2*MStau2(1)*(9*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(3*MStau2(1) - 31*MStau2(2)) + 
     &         3*MStau2(2)**2))*MStau2(3)**2 - 
     &   2*MStau2(1)*(-MStau2(2)**2 + 
     &      MStau2(1)*(MStau2(1) + 2*MStau2(2)))*MStau2(3)**4 - 
     &   MZ2**5*MStau2(2)*MStau2(4)

	Sf9810 = MZ2**5*(16*MStau2(2)**2 + 
     &      MStau2(1)*(63*MStau2(1) + 163*MStau2(2))) - 
     &   MZ2**4*(48*MStau2(2)**3 + 
     &      MStau2(1)*(944*MStau2(2)**2 + 
     &         MStau2(1)*(337*MStau2(1) + 849*MStau2(2)))) + 
     &   MZ2**3*(48*MStau2(2)**4 + 
     &      MStau2(1)*(1881*MStau2(2)**3 + 
     &         MStau2(1)*(1675*MStau2(2)**2 + 
     &            MStau2(1)*(605*MStau2(1) + 1599*MStau2(2))))) - 
     &   MZ2**2*(16*MStau2(2)**5 + 
     &      MStau2(1)*(1494*MStau2(2)**4 + 
     &         MStau2(1)*(819*MStau2(2)**3 + 
     &            MStau2(1)*
     &             (-257*MStau2(2)**2 + 
     &               MStau2(1)*(279*MStau2(1) + 1521*MStau2(2))))))
     &     - 2*MZ2*MStau2(1)*
     &    (-249*MStau2(2)**3 + 
     &      MStau2(1)*(MStau2(1)*(48*MStau2(1) - 136*MStau2(2)) - 
     &         39*MStau2(2)**2))*MStau2(3)**2 + 
     &   4*MStau2(1)*(MStau2(1)*(11*MStau2(1) - 32*MStau2(2)) - 
     &      26*MStau2(2)**2)*MStau2(3)**4

	Sf10023 = -4*MZ2**5*
     &    (-4*MStau2(2)**3 + 
     &      MStau2(1)*(12*MStau2(2)**2 + 
     &         MStau2(1)*(9*MStau2(1) + 83*MStau2(2)))) + 
     &   2*MZ2**4*(-40*MStau2(2)**4 + 
     &      MStau2(1)*(19*MStau2(2)**3 + 
     &         MStau2(1)*(1281*MStau2(2)**2 + 
     &            MStau2(1)*(47*MStau2(1) + 913*MStau2(2))))) - 
     &   2*MZ2**3*(-72*MStau2(2)**5 + 
     &      MStau2(1)*(-134*MStau2(2)**4 + 
     &         MStau2(1)*(1471*MStau2(2)**3 + 
     &            MStau2(1)*
     &             (3454*MStau2(2)**2 + 
     &               MStau2(1)*(58*MStau2(1) + 1591*MStau2(2))))))+
     &     MZ2**2*(-112*MStau2(2)**6 + 
     &      MStau2(1)*(-595*MStau2(2)**5 + 
     &         MStau2(1)*(1413*MStau2(2)**4 + 
     &            MStau2(1)*
     &             (2110*MStau2(2)**3 + 
     &               MStau2(1)*
     &                (2870*MStau2(2)**2 + 
     &                  MStau2(1)*(69*MStau2(1) + 3333*MStau2(2))))
     &            ))) + MZ2**6*MStau2(1)*
     &    (-14*MStau2(2) - 5*MStau2(3)) - 
     &   2*MZ2*(-16*MStau2(2)**5 + 
     &      MStau2(1)*(-358*MStau2(2)**4 + 
     &         MStau2(1)*(899*MStau2(2)**3 + 
     &            MStau2(1)*
     &             (1022*MStau2(2)**2 + 
     &               MStau2(1)*(8*MStau2(1) + 957*MStau2(2))))))*
     &    MStau2(3)**2 + 4*MStau2(1)*MStau2(2)*
     &    (-74*MStau2(2)**2 + 
     &      MStau2(1)*(64*MStau2(1) + 61*MStau2(2)))*MStau2(3)**4

	Sf8865 = Sf8864*MSbot2(1) - 414*MSneu2(1) - 
     &   414*MStau2(2) - 832*MStop2(1) - 826*MStop2(2) + 
     &   MSbot2(2)*(-826 + 618*USbot22**2) + 
     &   12*SW2**2*(-1 + UStau22**2)*
     &    (535*MStau2(1) + 161*Sf8536*UStau22**2) + 
     &   UStau22**2*(-8361*MStau2(1) - 6717*MStau2(2) + 
     &      9*(562*MStau2(1) + 723*MStau2(2))*UStau22**2) - 
     &   966*Sf8534*SW2**2*
     &    (1 + UStau22**2*(-2 + UStau22**2)) + 
     &   6*(MStop2(1) - MStop2(2))*UStop22**2

	Sf7278 = Sf7245*Sf7246 + 3*Sf7256*Sf7257 + 
     &   3*Sf4804*MSbot2(2) + Sf7247*MSneu2(1) + 
     &   Sf1768*MStau2(2) + 3*Sf7212*MStop2(2) + 
     &   (Sf7242*Sf7243 - Sf7245*Sf7246 + Sf1738*MStau2(1) - 
     &      Sf1768*MStau2(2))*UStau22**2 + 
     &   3*(Sf7251*Sf7252 - Sf7256*Sf7257 + Sf7189*MStop2(1) - 
     &      Sf7212*MStop2(2))*UStop22**2 + 
     &   3*USbot22**2*(Sf4783*MSbot2(1) - Sf4804*MSbot2(2) + 
     &      Sf7249*Sf7250*UStop22**2 - 
     &      Sf7251*Sf7252*UStop22**2 + 
     &      Sf7254*Sf7255*(1 - UStop22**2) + 
     &      Sf7256*Sf7257*(-1 + UStop22**2))

	Sf8719 = Sf7501**2*MSbot2(1) - Sf7503*MSbot2(1) + 
     &   Sf7521**2*MSbot2(2) - Sf7522*MSbot2(2) + 
     &   Sf7509**2*MSbot2(4) + 3*Sf7579**2*MStau2(1) - 
     &   3*Sf7580*MStau2(1) + 3*Sf7713**2*MStau2(2) - 
     &   3*Sf7714*MStau2(2) + 3*Sf7701**2*MStau2(4) + 
     &   Sf7990**2*MStop2(1) - Sf7991*MStop2(1) + 
     &   Sf8532**2*MStop2(2) - Sf8533*MStop2(2) + 
     &   Sf8520**2*MStop2(4)

	Sf7202 = UStop11*(Sf7201 + 4*SW2*(-1 + UStop22**2))

	Sf1758 = UStau11*(Sf1757 + 2*SW2*(-1 + UStau22**2))

	Sf8508 = USbot21*
     &   (MZ2*Sf8507 + 2*MW2*Sf4795*Sf8491*USbot22**2)

	Sf4776 = USbot11*(Sf4775 + 2*SW2*(-1 + USbot22**2))

	Sf7422 = Sf7401*MSbot2(2)

	Sf7429 = MZ2*Sf7428

	Sf7426 = 2*MZ2

	Sf7413 = MZ2*Sf7412

	Sf7409 = 2*MZ2

	Sf7384 = 2*MZ2

	Sf7378 = MZ2*Sf7377

	Sf7361 = 2*MZ2

	Sf7358 = MZ2*Sf7357

	Sf7344 = MZ2*Sf7343

	Sf7341 = 2*MZ2

	Sf7332 = MZ2*Sf7331

	Sf7328 = 2*MZ2

	Sf7303 = MZ2*Sf7301

	Sf8733 = (-(Sf7756*Sf86**2) + 4*Sf8722)*UStop12*UStop22

	Sf7762 = Sf7320*(2*Sf7757 - Sf7756*Ss443**2)

	Sf8820 = (-Sf8534 + Sf8536)*UStau21*UStau22

	Sf7388 = -(Sf7370*MStop2(2))

	Sf7372 = Sf7388

	Sf112 = -((Sf86 + 4*MZ2*SW2)*UStop12*UStop22)

	Sf25 = -((Sf15 + 2*MZ2*SW2)*USbot12*USbot22)

	Sf8640 = -3*(10*Sf7314*MSneu2(1) - 
     &     7*MStau2(2)*UStau22**2*(-1 + UStau22**2))

	Sf8114 = -1121*MSq2Diff(2,2)**8 + 
     &   MW2*(81*Sf8108*MSq2Diff(2,2)**6 + 
     &      MW2*(-2*Sf8110*MSq2Diff(2,2)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(348*MW2 - 12*Sf8107) + Sf8109) + 
     &                  Sf8112) - Sf8113) + 
     &            2*Sf8111*MSq2Diff(2,2)**2)))

	Sf8036 = -1121*MSq2Diff(1,2)**8 + 
     &   MW2*(81*Sf8030*MSq2Diff(1,2)**6 + 
     &      MW2*(-2*Sf8032*MSq2Diff(1,2)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*
     &                   (MW2*(348*MW2 - 12*Sf8029) + 2*Sf8031) + 
     &                  Sf8034) - Sf8035) + 
     &            2*Sf8033*MSq2Diff(1,2)**2)))

	Sf7905 = -1121*MSq2Diff(2,1)**8 + 
     &   MW2*(81*Sf7899*MSq2Diff(2,1)**6 + 
     &      MW2*(-2*Sf7901*MSq2Diff(2,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(348*MW2 - 12*Sf7898) + Sf7900) + 
     &                  Sf7903) - Sf7904) + 
     &            2*Sf7902*MSq2Diff(2,1)**2)))

	Sf7827 = -1121*MSq2Diff(1,1)**8 + 
     &   MW2*(81*Sf7821*MSq2Diff(1,1)**6 + 
     &      MW2*(-2*Sf7823*MSq2Diff(1,1)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(348*MW2 - 12*Sf7820) + Sf7822) + 
     &                  Sf7825) - Sf7826) + 
     &            2*Sf7824*MSq2Diff(1,1)**2)))

	Sf7598 = -427*MSl2Diff(2,1)**8 + 
     &   MW2*(4*Sf7592*MSl2Diff(2,1)**6 + 
     &      MW2*(-(Sf7594*MSl2Diff(2,1)**4) + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(98*MW2 - Sf7591) + 5*Sf7593) + 
     &                  2*Sf7596) - Sf7597) + 
     &            4*Sf7595*MSl2Diff(2,1)**2)))

	Sf7540 = -427*MSl2Diff(1,1)**8 + 
     &   MW2*(4*Sf7534*MSl2Diff(1,1)**6 + 
     &      MW2*(-(Sf7536*MSl2Diff(1,1)**4) + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(98*MW2 - Sf7533) + 5*Sf7535) + 
     &                  2*Sf7538) - Sf7539) + 
     &            4*Sf7537*MSl2Diff(1,1)**2)))

	Sf8129 = 243*MSq2Diff(2,2)**8 + 
     &   MW2*(-(Sf8123*MSq2Diff(2,2)**6) + 
     &      MW2*(2*Sf8125*MSq2Diff(2,2)**4 + 
     &         MW2*(MW2*(MW2*
     &                (MW2*(MW2*(8*MW2 + Sf8122) + 2*Sf8124) - 
     &                  Sf8126) + Sf8128) - 
     &            2*Sf8127*MSq2Diff(2,2)**2)))

	end

