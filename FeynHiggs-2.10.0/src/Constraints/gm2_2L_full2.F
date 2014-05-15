* gm2_2L_full2.F
* this file is part of FeynHiggs
* generated 22 Dec 2005 10:40
* last modified 10 Aug 11 th

#include "externals.h"
#include "types.h"


	subroutine gm2_2L_full2
	implicit none

#include "gm2_2L_full.h"

	Sf7565 = 2887*MStau2(1)**2 + 
     &   MSneu2(1)*(1375*MSneu2(1) + 1420*MStau2(1))

	Sf1702 = 913*MStau2(2)**2 + 
     &   MSneu2(1)*(913*MSneu2(1) + 1096*MStau2(2))

	Sf903 = 913*MStau2(1)**2 + 
     &   MSneu2(1)*(913*MSneu2(1) + 1096*MStau2(1))

	Sf1708 = 886*MStau2(2)**2 + 
     &   MSneu2(1)*(886*MSneu2(1) + 149*MStau2(2))

	Sf909 = 886*MStau2(1)**2 + 
     &   MSneu2(1)*(886*MSneu2(1) + 149*MStau2(1))

	Sf7608 = 599*MStau2(2)**2 + 
     &   MSneu2(1)*(599*MSneu2(1) + 818*MStau2(2))

	Sf7550 = 599*MStau2(1)**2 + 
     &   MSneu2(1)*(599*MSneu2(1) + 818*MStau2(1))

	Sf7695 = 554*MStau2(2)**2 + 
     &   MSneu2(1)*(589*MSneu2(1) + 705*MStau2(2))

	Sf7681 = 554*MStau2(1)**2 + 
     &   MSneu2(1)*(589*MSneu2(1) + 705*MStau2(1))

	Sf7609 = 1007*MStau2(2)**2 + 
     &   MSneu2(1)*(539*MSneu2(1) + 342*MStau2(2))

	Sf7551 = 1007*MStau2(1)**2 + 
     &   MSneu2(1)*(539*MSneu2(1) + 342*MStau2(1))

	Sf7742 = 685*MStau2(2)**2 + 
     &   MSneu2(1)*(470*MSneu2(1) + 693*MStau2(2))

	Sf7718 = 685*MStau2(1)**2 + 
     &   MSneu2(1)*(470*MSneu2(1) + 693*MStau2(1))

	Sf7659 = 535*MStau2(2)**2 + 
     &   MSneu2(1)*(427*MSneu2(1) + 462*MStau2(2))

	Sf7647 = 535*MStau2(1)**2 + 
     &   MSneu2(1)*(427*MSneu2(1) + 462*MStau2(1))

	Sf1278 = 371*MStau2(2)**2 + 
     &   MSneu2(1)*(371*MSneu2(1) + 414*MStau2(2))

	Sf318 = 371*MStau2(1)**2 + 
     &   MSneu2(1)*(371*MSneu2(1) + 414*MStau2(1))

	Sf1644 = 449*MStau2(2)**2 + 
     &   MSneu2(1)*(341*MSneu2(1) + 330*MStau2(2))

	Sf819 = 449*MStau2(1)**2 + 
     &   MSneu2(1)*(341*MSneu2(1) + 330*MStau2(1))

	Sf1391 = 307*MStau2(2)**2 + 
     &   MSneu2(1)*(307*MSneu2(1) + 418*MStau2(2))

	Sf474 = 307*MStau2(1)**2 + 
     &   MSneu2(1)*(307*MSneu2(1) + 418*MStau2(1))

	Sf1397 = 274*MStau2(2)**2 + 
     &   MSneu2(1)*(274*MSneu2(1) + 23*MStau2(2))

	Sf480 = 274*MStau2(1)**2 + 
     &   MSneu2(1)*(274*MSneu2(1) + 23*MStau2(1))

	Sf1342 = 61*MStau2(2)**2 + 
     &   MSneu2(1)*(272*MSneu2(1) + 315*MStau2(2))

	Sf405 = 61*MStau2(1)**2 + 
     &   MSneu2(1)*(272*MSneu2(1) + 315*MStau2(1))

	Sf1312 = 248*MStau2(2)**2 + 
     &   MSneu2(1)*(265*MSneu2(1) + 375*MStau2(2))

	Sf7664 = 210*MStau2(2)**2 + 
     &   MSneu2(1)*(265*MSneu2(1) + 237*MStau2(2))

	Sf361 = 248*MStau2(1)**2 + 
     &   MSneu2(1)*(265*MSneu2(1) + 375*MStau2(1))

	Sf7652 = 210*MStau2(1)**2 + 
     &   MSneu2(1)*(265*MSneu2(1) + 237*MStau2(1))

	Sf1643 = 287*MStau2(2)**2 + 
     &   MSneu2(1)*(251*MSneu2(1) + 326*MStau2(2))

	Sf818 = 287*MStau2(1)**2 + 
     &   MSneu2(1)*(251*MSneu2(1) + 326*MStau2(1))

	Sf1336 = 81*MStau2(2)**2 + 
     &   MSneu2(1)*(235*MSneu2(1) + 174*MStau2(2))

	Sf399 = 81*MStau2(1)**2 + 
     &   MSneu2(1)*(235*MSneu2(1) + 174*MStau2(1))

	Sf1335 = 51*MStau2(2)**2 + 
     &   MSneu2(1)*(208*MSneu2(1) + 171*MStau2(2))

	Sf398 = 51*MStau2(1)**2 + 
     &   MSneu2(1)*(208*MSneu2(1) + 171*MStau2(1))

	Sf7687 = 233*MStau2(2)**2 + 
     &   MSneu2(1)*(173*MSneu2(1) + 170*MStau2(2))

	Sf7673 = 233*MStau2(1)**2 + 
     &   MSneu2(1)*(173*MSneu2(1) + 170*MStau2(1))

	Sf2148 = 176*MStau2(2)**2 + 
     &   MSneu2(1)*(167*MSneu2(1) + 233*MStau2(2))

	Sf2068 = 176*MStau2(1)**2 + 
     &   MSneu2(1)*(167*MSneu2(1) + 233*MStau2(1))

	Sf2147 = 116*MStau2(2)**2 + 
     &   MSneu2(1)*(165*MSneu2(1) + 127*MStau2(2))

	Sf2067 = 116*MStau2(1)**2 + 
     &   MSneu2(1)*(165*MSneu2(1) + 127*MStau2(1))

	Sf1341 = 47*MStau2(2)**2 + 
     &   MSneu2(1)*(164*MSneu2(1) + 173*MStau2(2))

	Sf404 = 47*MStau2(1)**2 + 
     &   MSneu2(1)*(164*MSneu2(1) + 173*MStau2(1))

	Sf7593 = -809*MStau2(2)**2 + 
     &   MSneu2(1)*(163*MSneu2(1) + 364*MStau2(2))

	Sf1277 = 163*MStau2(2)**2 + 
     &   MSneu2(1)*(163*MSneu2(1) + 180*MStau2(2))

	Sf7535 = -809*MStau2(1)**2 + 
     &   MSneu2(1)*(163*MSneu2(1) + 364*MStau2(1))

	Sf317 = 163*MStau2(1)**2 + 
     &   MSneu2(1)*(163*MSneu2(1) + 180*MStau2(1))

	Sf8329 = 325*MStau2(2)**2 + 
     &   MSneu2(1)*(162*MSneu2(1) + 225*MStau2(2))

	Sf8322 = 325*MStau2(1)**2 + 
     &   MSneu2(1)*(162*MSneu2(1) + 225*MStau2(1))

	Sf1401 = MSneu2(1)*(111*MSneu2(1) - 66*MStau2(2)) - 
     &   177*MStau2(2)**2

	Sf484 = MSneu2(1)*(111*MSneu2(1) - 66*MStau2(1)) - 
     &   177*MStau2(1)**2

	Sf1282 = -23*MStau2(2)**2 + 
     &   MSneu2(1)*(109*MSneu2(1) + 202*MStau2(2))

	Sf322 = -23*MStau2(1)**2 + 
     &   MSneu2(1)*(109*MSneu2(1) + 202*MStau2(1))

	Sf7686 = 103*MStau2(2)**2 + 
     &   MSneu2(1)*(103*MSneu2(1) + 82*MStau2(2))

	Sf7672 = 103*MStau2(1)**2 + 
     &   MSneu2(1)*(103*MSneu2(1) + 82*MStau2(1))

	Sf1311 = 110*MStau2(2)**2 + 
     &   MSneu2(1)*(101*MSneu2(1) + 173*MStau2(2))

	Sf1283 = 113*MStau2(2)**2 + 
     &   MSneu2(1)*(101*MSneu2(1) + 74*MStau2(2))

	Sf7694 = -69*MStau2(2)**2 + 
     &   MSneu2(1)*(101*MSneu2(1) + 9*MStau2(2))

	Sf360 = 110*MStau2(1)**2 + 
     &   MSneu2(1)*(101*MSneu2(1) + 173*MStau2(1))

	Sf323 = 113*MStau2(1)**2 + 
     &   MSneu2(1)*(101*MSneu2(1) + 74*MStau2(1))

	Sf7680 = -69*MStau2(1)**2 + 
     &   MSneu2(1)*(101*MSneu2(1) + 9*MStau2(1))

	Sf1844 = 91*MStau2(2)**2 + 
     &   MSneu2(1)*(91*MSneu2(1) + 102*MStau2(2))

	Sf1065 = 91*MStau2(1)**2 + 
     &   MSneu2(1)*(91*MSneu2(1) + 102*MStau2(1))

	Sf1642 = 85*MStau2(2)**2 + 
     &   MSneu2(1)*(85*MSneu2(1) + 98*MStau2(2))

	Sf817 = 85*MStau2(1)**2 + 
     &   MSneu2(1)*(85*MSneu2(1) + 98*MStau2(1))

	Sf1629 = MSneu2(1)*(84*MSneu2(1) - 132*MStau2(2)) - 
     &   132*MStau2(2)**2

	Sf804 = MSneu2(1)*(84*MSneu2(1) - 132*MStau2(1)) - 
     &   132*MStau2(1)**2

	Sf1306 = 235*MStau2(2)**2 + 
     &   MSneu2(1)*(81*MSneu2(1) + 174*MStau2(2))

	Sf355 = 235*MStau2(1)**2 + 
     &   MSneu2(1)*(81*MSneu2(1) + 174*MStau2(1))

	Sf1850 = 99*MStau2(2)**2 + 
     &   MSneu2(1)*(79*MSneu2(1) + 78*MStau2(2))

	Sf1390 = 79*MStau2(2)**2 + 
     &   MSneu2(1)*(79*MSneu2(1) + 64*MStau2(2))

	Sf1071 = 99*MStau2(1)**2 + 
     &   MSneu2(1)*(79*MSneu2(1) + 78*MStau2(1))

	Sf473 = 79*MStau2(1)**2 + 
     &   MSneu2(1)*(79*MSneu2(1) + 64*MStau2(1))

	Sf2143 = 76*MStau2(2)**2 + 
     &   MSneu2(1)*(69*MSneu2(1) + 89*MStau2(2))

	Sf2063 = 76*MStau2(1)**2 + 
     &   MSneu2(1)*(69*MSneu2(1) + 89*MStau2(1))

	Sf7601 = -403*MStau2(2)**2 + 
     &   MSneu2(1)*(65*MSneu2(1) + 332*MStau2(2))

	Sf7543 = -403*MStau2(1)**2 + 
     &   MSneu2(1)*(65*MSneu2(1) + 332*MStau2(1))

	Sf7616 = 70*MStau2(2)**2 + 
     &   MSneu2(1)*(64*MSneu2(1) + 109*MStau2(2))

	Sf7558 = 70*MStau2(1)**2 + 
     &   MSneu2(1)*(64*MSneu2(1) + 109*MStau2(1))

	Sf1396 = MSneu2(1)*(61*MSneu2(1) - 16*MStau2(2)) + 
     &   61*MStau2(2)**2

	Sf479 = MSneu2(1)*(61*MSneu2(1) - 16*MStau2(1)) + 
     &   61*MStau2(1)**2

	Sf1958 = 36*MStau2(2)**2 + 
     &   MSneu2(1)*(55*MSneu2(1) + 51*MStau2(2))

	Sf1896 = 36*MStau2(1)**2 + 
     &   MSneu2(1)*(55*MSneu2(1) + 51*MStau2(1))

	Sf7632 = 65*MStau2(2)**2 + 
     &   MSneu2(1)*(53*MSneu2(1) + 68*MStau2(2))

	Sf7574 = 65*MStau2(1)**2 + 
     &   MSneu2(1)*(53*MSneu2(1) + 68*MStau2(1))

	Sf1305 = 208*MStau2(2)**2 + 
     &   MSneu2(1)*(51*MSneu2(1) + 171*MStau2(2))

	Sf354 = 208*MStau2(1)**2 + 
     &   MSneu2(1)*(51*MSneu2(1) + 171*MStau2(1))

	Sf1963 = 39*MStau2(2)**2 + 
     &   MSneu2(1)*(50*MSneu2(1) + 39*MStau2(2))

	Sf1901 = 39*MStau2(1)**2 + 
     &   MSneu2(1)*(50*MSneu2(1) + 39*MStau2(1))

	Sf2142 = -25*MStau2(2)**2 + 
     &   MSneu2(1)*(45*MSneu2(1) + 14*MStau2(2))

	Sf2062 = -25*MStau2(1)**2 + 
     &   MSneu2(1)*(45*MSneu2(1) + 14*MStau2(1))

	Sf1633 = MSneu2(1)*(42*MSneu2(1) - 66*MStau2(2)) - 
     &   66*MStau2(2)**2

	Sf808 = MSneu2(1)*(42*MSneu2(1) - 66*MStau2(1)) - 
     &   66*MStau2(1)**2

	Sf1389 = 39*MStau2(2)**2 + 
     &   MSneu2(1)*(40*MSneu2(1) + 54*MStau2(2))

	Sf7663 = MSneu2(1)*(40*MSneu2(1) - 6*MStau2(2)) - 
     &   12*MStau2(2)**2

	Sf472 = 39*MStau2(1)**2 + 
     &   MSneu2(1)*(40*MSneu2(1) + 54*MStau2(1))

	Sf7651 = MSneu2(1)*(40*MSneu2(1) - 6*MStau2(1)) - 
     &   12*MStau2(1)**2

	Sf7658 = 59*MStau2(2)**2 + 
     &   MSneu2(1)*(37*MSneu2(1) + 34*MStau2(2))

	Sf1410 = MSneu2(1)*(37*MSneu2(1) - 22*MStau2(2)) - 
     &   59*MStau2(2)**2

	Sf7646 = 59*MStau2(1)**2 + 
     &   MSneu2(1)*(37*MSneu2(1) + 34*MStau2(1))

	Sf493 = MSneu2(1)*(37*MSneu2(1) - 22*MStau2(1)) - 
     &   59*MStau2(1)**2

	Sf1849 = 57*MStau2(2)**2 + 
     &   MSneu2(1)*(29*MSneu2(1) + 10*MStau2(2))

	Sf1070 = 57*MStau2(1)**2 + 
     &   MSneu2(1)*(29*MSneu2(1) + 10*MStau2(1))

	Sf1641 = 25*MStau2(2)**2 + 
     &   MSneu2(1)*(25*MSneu2(1) + 28*MStau2(2))

	Sf1968 = 12*MStau2(2)**2 + 
     &   MSneu2(1)*(25*MSneu2(1) + 3*MStau2(2))

	Sf816 = 25*MStau2(1)**2 + 
     &   MSneu2(1)*(25*MSneu2(1) + 28*MStau2(1))

	Sf1908 = 12*MStau2(1)**2 + 
     &   MSneu2(1)*(25*MSneu2(1) + 3*MStau2(1))

	Sf1327 = 11*MStau2(2)**2 + 
     &   MSneu2(1)*(23*MSneu2(1) + 5*MStau2(2))

	Sf390 = 11*MStau2(1)**2 + MSneu2(1)*(23*MSneu2(1) + 5*MStau2(1))

	Sf1957 = -3*MStau2(2)**2 + 
     &   MSneu2(1)*(19*MSneu2(1) + 6*MStau2(2))

	Sf1466 = MSneu2(1)*(19*MSneu2(1) - 16*MStau2(2)) + 
     &   19*MStau2(2)**2

	Sf1895 = -3*MStau2(1)**2 + 
     &   MSneu2(1)*(19*MSneu2(1) + 6*MStau2(1))

	Sf563 = MSneu2(1)*(19*MSneu2(1) - 16*MStau2(1)) + 
     &   19*MStau2(1)**2

	Sf8328 = -50*MStau2(2)**2 + 
     &   MSneu2(1)*(18*MSneu2(1) + 21*MStau2(2))

	Sf8321 = -50*MStau2(1)**2 + 
     &   MSneu2(1)*(18*MSneu2(1) + 21*MStau2(1))

	Sf1543 = 19*MStau2(2)**2 + 
     &   MSneu2(1)*(15*MSneu2(1) + 66*MStau2(2))

	Sf1474 = MSneu2(1)*(15*MSneu2(1) - 6*MStau2(2)) - 
     &   45*MStau2(2)**2

	Sf669 = 19*MStau2(1)**2 + 
     &   MSneu2(1)*(15*MSneu2(1) + 66*MStau2(1))

	Sf573 = MSneu2(1)*(15*MSneu2(1) - 6*MStau2(1)) - 45*MStau2(1)**2

	Sf1422 = 13*MStau2(2)**2 + 
     &   MSneu2(1)*(13*MSneu2(1) + 12*MStau2(2))

	Sf511 = 13*MStau2(1)**2 + 
     &   MSneu2(1)*(13*MSneu2(1) + 12*MStau2(1))

	Sf1463 = 11*MStau2(2)**2 + 
     &   MSneu2(1)*(11*MSneu2(1) + 14*MStau2(2))

	Sf2164 = -7*MStau2(2)**2 + 
     &   MSneu2(1)*(11*MSneu2(1) + 11*MStau2(2))

	Sf1355 = 23*MStau2(2)**2 + 
     &   MSneu2(1)*(11*MSneu2(1) + 5*MStau2(2))

	Sf2001 = MSneu2(1)*(11*MSneu2(1) - 24*MStau2(2)) + 
     &   13*MStau2(2)**2

	Sf560 = 11*MStau2(1)**2 + 
     &   MSneu2(1)*(11*MSneu2(1) + 14*MStau2(1))

	Sf2095 = -7*MStau2(1)**2 + 
     &   MSneu2(1)*(11*MSneu2(1) + 11*MStau2(1))

	Sf432 = 23*MStau2(1)**2 + MSneu2(1)*(11*MSneu2(1) + 5*MStau2(1))

	Sf1945 = MSneu2(1)*(11*MSneu2(1) - 24*MStau2(1)) + 
     &   13*MStau2(1)**2

	Sf1833 = 17*MStau2(2)**2 + 
     &   MSneu2(1)*(9*MSneu2(1) + 22*MStau2(2))

	Sf7741 = MSneu2(1)*(9*MSneu2(1) - 99*MStau2(2)) - 
     &   239*MStau2(2)**2

	Sf1044 = 17*MStau2(1)**2 + 
     &   MSneu2(1)*(9*MSneu2(1) + 22*MStau2(1))

	Sf7717 = MSneu2(1)*(9*MSneu2(1) - 99*MStau2(1)) - 
     &   239*MStau2(1)**2

	Sf7589 = 7*MStau2(2)**2 + MSneu2(1)*(7*MSneu2(1) + 12*MStau2(2))

	Sf1701 = 7*MStau2(2)**2 + MSneu2(1)*(7*MSneu2(1) + 10*MStau2(2))

	Sf1465 = MSneu2(1)*(7*MSneu2(1) - MStau2(2)) + 7*MStau2(2)**2

	Sf7531 = 7*MStau2(1)**2 + MSneu2(1)*(7*MSneu2(1) + 12*MStau2(1))

	Sf902 = 7*MStau2(1)**2 + MSneu2(1)*(7*MSneu2(1) + 10*MStau2(1))

	Sf562 = MSneu2(1)*(7*MSneu2(1) - MStau2(1)) + 7*MStau2(1)**2

	Sf1470 = 6*MStau2(2)**2 + MSneu2(1)*(5*MSneu2(1) + 6*MStau2(2))

	Sf2128 = MSneu2(1)*(5*MSneu2(1) - 3*MStau2(2)) - 2*MStau2(2)**2

	Sf568 = 6*MStau2(1)**2 + MSneu2(1)*(5*MSneu2(1) + 6*MStau2(1))

	Sf2048 = MSneu2(1)*(5*MSneu2(1) - 3*MStau2(1)) - 2*MStau2(1)**2

	Sf7615 = 10*MStau2(2)**2 + MSneu2(1)*(4*MSneu2(1) + 7*MStau2(2))

	Sf1700 = 3*MStau2(2)**2 + MSneu2(1)*(4*MSneu2(1) + 6*MStau2(2))

	Sf7557 = 10*MStau2(1)**2 + MSneu2(1)*(4*MSneu2(1) + 7*MStau2(1))

	Sf901 = 3*MStau2(1)**2 + MSneu2(1)*(4*MSneu2(1) + 6*MStau2(1))

	Sf1694 = 3*MStau2(2)**2 + MSneu2(1)*(3*MSneu2(1) + 4*MStau2(2))

	Sf1548 = 3*MStau2(2)**2 + MSneu2(1)*(3*MSneu2(1) + 2*MStau2(2))

	Sf1986 = MSneu2(1)*(3*MSneu2(1) - 4*MStau2(2)) + MStau2(2)**2

	Sf895 = 3*MStau2(1)**2 + MSneu2(1)*(3*MSneu2(1) + 4*MStau2(1))

	Sf674 = 3*MStau2(1)**2 + MSneu2(1)*(3*MSneu2(1) + 2*MStau2(1))

	Sf1927 = MSneu2(1)*(3*MSneu2(1) - 4*MStau2(1)) + MStau2(1)**2

	Sf1399 = MStau2(2)**2 + MSneu2(1)*(MSneu2(1) + 3*MStau2(2))

	Sf7631 = MSneu2(1)*(MSneu2(1) - 4*MStau2(2)) + 5*MStau2(2)**2

	Sf482 = MStau2(1)**2 + MSneu2(1)*(MSneu2(1) + 3*MStau2(1))

	Sf7573 = MSneu2(1)*(MSneu2(1) - 4*MStau2(1)) + 5*MStau2(1)**2

	Sf5365 = MSbot2(2)**2*(28*MSbot2(2) - 33*MStop2(2)) + 
     &   5*MStop2(2)**3

	Sf2954 = MSbot2(2)**2*(28*MSbot2(2) - 33*MStop2(1)) + 
     &   5*MStop2(1)**3

	Sf5399 = 28*MStop2(2)**3 + 
     &   MSbot2(2)*(5*MSbot2(2)**2 - 33*MStop2(2)**2)

	Sf2991 = 28*MStop2(1)**3 + 
     &   MSbot2(2)*(5*MSbot2(2)**2 - 33*MStop2(1)**2)

	Sf8110 = 13854*MStop2(2)**2 + 
     &   MSbot2(2)*(13125*MSbot2(2) + 18434*MStop2(2))

	Sf7901 = 13854*MStop2(1)**2 + 
     &   MSbot2(2)*(13125*MSbot2(2) + 18434*MStop2(1))

	Sf8133 = 11490*MStop2(2)**2 + 
     &   MSbot2(2)*(11004*MSbot2(2) + 15221*MStop2(2))

	Sf7924 = 11490*MStop2(1)**2 + 
     &   MSbot2(2)*(11004*MSbot2(2) + 15221*MStop2(1))

	Sf8125 = 2752*MStop2(2)**2 + 
     &   MSbot2(2)*(2509*MSbot2(2) + 3362*MStop2(2))

	Sf7916 = 2752*MStop2(1)**2 + 
     &   MSbot2(2)*(2509*MSbot2(2) + 3362*MStop2(1))

	Sf8392 = 1952*MStop2(2)**2 + 
     &   MSbot2(2)*(2271*MSbot2(2) + 2595*MStop2(2))

	Sf8368 = 1952*MStop2(1)**2 + 
     &   MSbot2(2)*(2271*MSbot2(2) + 2595*MStop2(1))

	Sf8104 = 853*MStop2(2)**2 + 
     &   MSbot2(2)*(2203*MSbot2(2) + 534*MStop2(2))

	Sf7895 = 853*MStop2(1)**2 + 
     &   MSbot2(2)*(2203*MSbot2(2) + 534*MStop2(1))

	Sf5961 = MSbot2(2)*(1907*MSbot2(2) - 80*MStop2(2)) + 
     &   1907*MStop2(2)**2

	Sf3536 = MSbot2(2)*(1907*MSbot2(2) - 80*MStop2(1)) + 
     &   1907*MStop2(1)**2

	Sf8457 = 2487*MStop2(2)**2 + 
     &   MSbot2(2)*(1763*MSbot2(2) + 2568*MStop2(2))

	Sf8433 = 2487*MStop2(1)**2 + 
     &   MSbot2(2)*(1763*MSbot2(2) + 2568*MStop2(1))

	Sf8391 = -799*MStop2(2)**2 + 
     &   MSbot2(2)*(1743*MSbot2(2) + 282*MStop2(2))

	Sf8367 = -799*MStop2(1)**2 + 
     &   MSbot2(2)*(1743*MSbot2(2) + 282*MStop2(1))

	Sf8117 = 2738*MStop2(2)**2 + 
     &   MSbot2(2)*(1685*MSbot2(2) + 501*MStop2(2))

	Sf7908 = 2738*MStop2(1)**2 + 
     &   MSbot2(2)*(1685*MSbot2(2) + 501*MStop2(1))

	Sf8132 = 1781*MStop2(2)**2 + 
     &   MSbot2(2)*(1214*MSbot2(2) + 1128*MStop2(2))

	Sf7923 = 1781*MStop2(1)**2 + 
     &   MSbot2(2)*(1214*MSbot2(2) + 1128*MStop2(1))

	Sf8109 = -6116*MStop2(2)**2 + 
     &   MSbot2(2)*(1174*MSbot2(2) + 7470*MStop2(2))

	Sf7900 = -6116*MStop2(1)**2 + 
     &   MSbot2(2)*(1174*MSbot2(2) + 7470*MStop2(1))

	Sf5959 = 913*MStop2(2)**2 + 
     &   MSbot2(2)*(913*MSbot2(2) + 1096*MStop2(2))

	Sf3534 = 913*MStop2(1)**2 + 
     &   MSbot2(2)*(913*MSbot2(2) + 1096*MStop2(1))

	Sf5962 = 886*MStop2(2)**2 + 
     &   MSbot2(2)*(886*MSbot2(2) + 149*MStop2(2))

	Sf3537 = 886*MStop2(1)**2 + 
     &   MSbot2(2)*(886*MSbot2(2) + 149*MStop2(1))

	Sf8103 = 731*MStop2(2)**2 + 
     &   MSbot2(2)*(857*MSbot2(2) + 876*MStop2(2))

	Sf7894 = 731*MStop2(1)**2 + 
     &   MSbot2(2)*(857*MSbot2(2) + 876*MStop2(1))

	Sf8116 = 697*MStop2(2)**2 + 
     &   MSbot2(2)*(697*MSbot2(2) + 622*MStop2(2))

	Sf7907 = 697*MStop2(1)**2 + 
     &   MSbot2(2)*(697*MSbot2(2) + 622*MStop2(1))

	Sf8086 = MSbot2(2)*(565*MSbot2(2) - 360*MStop2(2)) + 
     &   133*MStop2(2)**2

	Sf7877 = MSbot2(2)*(565*MSbot2(2) - 360*MStop2(1)) + 
     &   133*MStop2(1)**2

	Sf8184 = 718*MStop2(2)**2 + 
     &   MSbot2(2)*(520*MSbot2(2) + 366*MStop2(2))

	Sf7975 = 718*MStop2(1)**2 + 
     &   MSbot2(2)*(520*MSbot2(2) + 366*MStop2(1))

	Sf8183 = 589*MStop2(2)**2 + 
     &   MSbot2(2)*(508*MSbot2(2) + 507*MStop2(2))

	Sf7974 = 589*MStop2(1)**2 + 
     &   MSbot2(2)*(508*MSbot2(2) + 507*MStop2(1))

	Sf8386 = -147*MStop2(2)**2 + 
     &   MSbot2(2)*(465*MSbot2(2) + 46*MStop2(2))

	Sf8362 = -147*MStop2(1)**2 + 
     &   MSbot2(2)*(465*MSbot2(2) + 46*MStop2(1))

	Sf8696 = 499*MStop2(2)**2 + 
     &   MSbot2(2)*(463*MSbot2(2) + 462*MStop2(2))

	Sf8656 = 499*MStop2(1)**2 + 
     &   MSbot2(2)*(463*MSbot2(2) + 462*MStop2(1))

	Sf8087 = 365*MStop2(2)**2 + 
     &   MSbot2(2)*(437*MSbot2(2) + 326*MStop2(2))

	Sf7878 = 365*MStop2(1)**2 + 
     &   MSbot2(2)*(437*MSbot2(2) + 326*MStop2(1))

	Sf8398 = 520*MStop2(2)**2 + 
     &   MSbot2(2)*(385*MSbot2(2) + 219*MStop2(2))

	Sf8374 = 520*MStop2(1)**2 + 
     &   MSbot2(2)*(385*MSbot2(2) + 219*MStop2(1))

	Sf5411 = 371*MStop2(2)**2 + 
     &   MSbot2(2)*(371*MSbot2(2) + 414*MStop2(2))

	Sf3003 = 371*MStop2(1)**2 + 
     &   MSbot2(2)*(371*MSbot2(2) + 414*MStop2(1))

	Sf8219 = 229*MStop2(2)**2 + 
     &   MSbot2(2)*(315*MSbot2(2) + 258*MStop2(2))

	Sf8381 = 297*MStop2(2)**2 + 
     &   MSbot2(2)*(315*MSbot2(2) + 244*MStop2(2))

	Sf8207 = 229*MStop2(1)**2 + 
     &   MSbot2(2)*(315*MSbot2(2) + 258*MStop2(1))

	Sf8357 = 297*MStop2(1)**2 + 
     &   MSbot2(2)*(315*MSbot2(2) + 244*MStop2(1))

	Sf5904 = 307*MStop2(2)**2 + 
     &   MSbot2(2)*(307*MSbot2(2) + 418*MStop2(2))

	Sf3460 = 307*MStop2(1)**2 + 
     &   MSbot2(2)*(307*MSbot2(2) + 418*MStop2(1))

	Sf8091 = MSbot2(2)*(299*MSbot2(2) - 220*MStop2(2)) + 
     &   137*MStop2(2)**2

	Sf7882 = MSbot2(2)*(299*MSbot2(2) - 220*MStop2(1)) + 
     &   137*MStop2(1)**2

	Sf5907 = 274*MStop2(2)**2 + 
     &   MSbot2(2)*(274*MSbot2(2) + 23*MStop2(2))

	Sf3463 = 274*MStop2(1)**2 + 
     &   MSbot2(2)*(274*MSbot2(2) + 23*MStop2(1))

	Sf8695 = 310*MStop2(2)**2 + 
     &   MSbot2(2)*(266*MSbot2(2) + 204*MStop2(2))

	Sf8655 = 310*MStop2(1)**2 + 
     &   MSbot2(2)*(266*MSbot2(2) + 204*MStop2(1))

	Sf5384 = 130*MStop2(2)**2 + 
     &   MSbot2(2)*(261*MSbot2(2) + 325*MStop2(2))

	Sf2976 = 130*MStop2(1)**2 + 
     &   MSbot2(2)*(261*MSbot2(2) + 325*MStop2(1))

	Sf8456 = MSbot2(2)*(259*MSbot2(2) - 1092*MStop2(2)) - 
     &   2985*MStop2(2)**2

	Sf8432 = MSbot2(2)*(259*MSbot2(2) - 1092*MStop2(1)) - 
     &   2985*MStop2(1)**2

	Sf8097 = MSbot2(2)*(257*MSbot2(2) - 30*MStop2(2)) + 
     &   95*MStop2(2)**2

	Sf7888 = MSbot2(2)*(257*MSbot2(2) - 30*MStop2(1)) + 
     &   95*MStop2(1)**2

	Sf8451 = MSbot2(2)*(237*MSbot2(2) - 28*MStop2(2)) - 
     &   285*MStop2(2)**2

	Sf8427 = MSbot2(2)*(237*MSbot2(2) - 28*MStop2(1)) - 
     &   285*MStop2(1)**2

	Sf8092 = 235*MStop2(2)**2 + 
     &   MSbot2(2)*(235*MSbot2(2) + 874*MStop2(2))

	Sf5378 = 81*MStop2(2)**2 + 
     &   MSbot2(2)*(235*MSbot2(2) + 174*MStop2(2))

	Sf7883 = 235*MStop2(1)**2 + 
     &   MSbot2(2)*(235*MSbot2(2) + 874*MStop2(1))

	Sf2970 = 81*MStop2(1)**2 + 
     &   MSbot2(2)*(235*MSbot2(2) + 174*MStop2(1))

	Sf6891 = 305*MStop2(2)**2 + 
     &   MSbot2(2)*(227*MSbot2(2) + 196*MStop2(2))

	Sf4686 = 305*MStop2(1)**2 + 
     &   MSbot2(2)*(227*MSbot2(2) + 196*MStop2(1))

	Sf5349 = 247*MStop2(2)**2 + 
     &   MSbot2(2)*(222*MSbot2(2) + 351*MStop2(2))

	Sf2938 = 247*MStop2(1)**2 + 
     &   MSbot2(2)*(222*MSbot2(2) + 351*MStop2(1))

	Sf5377 = 51*MStop2(2)**2 + 
     &   MSbot2(2)*(208*MSbot2(2) + 171*MStop2(2))

	Sf2969 = 51*MStop2(1)**2 + 
     &   MSbot2(2)*(208*MSbot2(2) + 171*MStop2(1))

	Sf8190 = 360*MStop2(2)**2 + 
     &   MSbot2(2)*(193*MSbot2(2) + 249*MStop2(2))

	Sf7981 = 360*MStop2(1)**2 + 
     &   MSbot2(2)*(193*MSbot2(2) + 249*MStop2(1))

	Sf6162 = 229*MStop2(2)**2 + 
     &   MSbot2(2)*(191*MSbot2(2) + 180*MStop2(2))

	Sf3914 = 229*MStop2(1)**2 + 
     &   MSbot2(2)*(191*MSbot2(2) + 180*MStop2(1))

	Sf6856 = 91*MStop2(2)**2 + 
     &   MSbot2(2)*(178*MSbot2(2) + 103*MStop2(2))

	Sf6742 = 91*MStop2(1)**2 + 
     &   MSbot2(2)*(178*MSbot2(2) + 103*MStop2(1))

	Sf6857 = 148*MStop2(2)**2 + 
     &   MSbot2(2)*(165*MSbot2(2) + 203*MStop2(2))

	Sf6743 = 148*MStop2(1)**2 + 
     &   MSbot2(2)*(165*MSbot2(2) + 203*MStop2(1))

	Sf5416 = 163*MStop2(2)**2 + 
     &   MSbot2(2)*(163*MSbot2(2) + 180*MStop2(2))

	Sf3016 = 163*MStop2(1)**2 + 
     &   MSbot2(2)*(163*MSbot2(2) + 180*MStop2(1))

	Sf7083 = 147*MStop2(2)**2 + 
     &   MSbot2(2)*(160*MSbot2(2) + 209*MStop2(2))

	Sf4652 = 147*MStop2(1)**2 + 
     &   MSbot2(2)*(160*MSbot2(2) + 209*MStop2(1))

	Sf5383 = 81*MStop2(2)**2 + 
     &   MSbot2(2)*(158*MSbot2(2) + 181*MStop2(2))

	Sf2975 = 81*MStop2(1)**2 + 
     &   MSbot2(2)*(158*MSbot2(2) + 181*MStop2(1))

	Sf6850 = -75*MStop2(2)**2 + 
     &   MSbot2(2)*(135*MSbot2(2) + 42*MStop2(2))

	Sf6736 = -75*MStop2(1)**2 + 
     &   MSbot2(2)*(135*MSbot2(2) + 42*MStop2(1))

	Sf5410 = 121*MStop2(2)**2 + 
     &   MSbot2(2)*(123*MSbot2(2) + 100*MStop2(2))

	Sf3002 = 121*MStop2(1)**2 + 
     &   MSbot2(2)*(123*MSbot2(2) + 100*MStop2(1))

	Sf8397 = 107*MStop2(2)**2 + 
     &   MSbot2(2)*(107*MSbot2(2) + 74*MStop2(2))

	Sf8373 = 107*MStop2(1)**2 + 
     &   MSbot2(2)*(107*MSbot2(2) + 74*MStop2(1))

	Sf5348 = 134*MStop2(2)**2 + 
     &   MSbot2(2)*(103*MSbot2(2) + 183*MStop2(2))

	Sf2937 = 134*MStop2(1)**2 + 
     &   MSbot2(2)*(103*MSbot2(2) + 183*MStop2(1))

	Sf8218 = MSbot2(2)*(99*MSbot2(2) - 54*MStop2(2)) - 
     &   35*MStop2(2)**2

	Sf8206 = MSbot2(2)*(99*MSbot2(2) - 54*MStop2(1)) - 
     &   35*MStop2(1)**2

	Sf7084 = 190*MStop2(2)**2 + 
     &   MSbot2(2)*(97*MSbot2(2) + 109*MStop2(2))

	Sf4653 = 190*MStop2(1)**2 + 
     &   MSbot2(2)*(97*MSbot2(2) + 109*MStop2(1))

	Sf6890 = 85*MStop2(2)**2 + 
     &   MSbot2(2)*(85*MSbot2(2) + 98*MStop2(2))

	Sf4685 = 85*MStop2(1)**2 + 
     &   MSbot2(2)*(85*MSbot2(2) + 98*MStop2(1))

	Sf5343 = 235*MStop2(2)**2 + 
     &   MSbot2(2)*(81*MSbot2(2) + 174*MStop2(2))

	Sf2932 = 235*MStop2(1)**2 + 
     &   MSbot2(2)*(81*MSbot2(2) + 174*MStop2(1))

	Sf5903 = 79*MStop2(2)**2 + 
     &   MSbot2(2)*(79*MSbot2(2) + 64*MStop2(2))

	Sf3459 = 79*MStop2(1)**2 + 
     &   MSbot2(2)*(79*MSbot2(2) + 64*MStop2(1))

	Sf7078 = 69*MStop2(2)**2 + 
     &   MSbot2(2)*(76*MSbot2(2) + 89*MStop2(2))

	Sf4647 = 69*MStop2(1)**2 + 
     &   MSbot2(2)*(76*MSbot2(2) + 89*MStop2(1))

	Sf8098 = 59*MStop2(2)**2 + 
     &   MSbot2(2)*(73*MSbot2(2) + 60*MStop2(2))

	Sf7889 = 59*MStop2(1)**2 + 
     &   MSbot2(2)*(73*MSbot2(2) + 60*MStop2(1))

	Sf8189 = -189*MStop2(2)**2 + 
     &   MSbot2(2)*(71*MSbot2(2) + 108*MStop2(2))

	Sf8380 = 71*MStop2(2)**2 + 
     &   MSbot2(2)*(71*MSbot2(2) + 50*MStop2(2))

	Sf7980 = -189*MStop2(1)**2 + 
     &   MSbot2(2)*(71*MSbot2(2) + 108*MStop2(1))

	Sf8356 = 71*MStop2(1)**2 + 
     &   MSbot2(2)*(71*MSbot2(2) + 50*MStop2(1))

	Sf6851 = 76*MStop2(2)**2 + 
     &   MSbot2(2)*(69*MSbot2(2) + 89*MStop2(2))

	Sf6737 = 76*MStop2(1)**2 + 
     &   MSbot2(2)*(69*MSbot2(2) + 89*MStop2(1))

	Sf8149 = 72*MStop2(2)**2 + 
     &   MSbot2(2)*(63*MSbot2(2) + 76*MStop2(2))

	Sf7940 = 72*MStop2(1)**2 + 
     &   MSbot2(2)*(63*MSbot2(2) + 76*MStop2(1))

	Sf5906 = MSbot2(2)*(61*MSbot2(2) - 16*MStop2(2)) + 
     &   61*MStop2(2)**2

	Sf3462 = MSbot2(2)*(61*MSbot2(2) - 16*MStop2(1)) + 
     &   61*MStop2(1)**2

	Sf8213 = MSbot2(2)*(55*MSbot2(2) - 14*MStop2(2)) - 
     &   13*MStop2(2)**2

	Sf8201 = MSbot2(2)*(55*MSbot2(2) - 14*MStop2(1)) - 
     &   13*MStop2(1)**2

	Sf5911 = 53*MStop2(2)**2 + 
     &   MSbot2(2)*(54*MSbot2(2) + 73*MStop2(2))

	Sf3467 = 53*MStop2(1)**2 + 
     &   MSbot2(2)*(54*MSbot2(2) + 73*MStop2(1))

	Sf5342 = 208*MStop2(2)**2 + 
     &   MSbot2(2)*(51*MSbot2(2) + 171*MStop2(2))

	Sf2931 = 208*MStop2(1)**2 + 
     &   MSbot2(2)*(51*MSbot2(2) + 171*MStop2(1))

	Sf5967 = MSbot2(2)*(46*MSbot2(2) - 9*MStop2(2)) - 
     &   52*MStop2(2)**2

	Sf3542 = MSbot2(2)*(46*MSbot2(2) - 9*MStop2(1)) - 
     &   52*MStop2(1)**2

	Sf8176 = 33*MStop2(2)**2 + 
     &   MSbot2(2)*(43*MSbot2(2) + 12*MStop2(2))

	Sf7967 = 33*MStop2(1)**2 + 
     &   MSbot2(2)*(43*MSbot2(2) + 12*MStop2(1))

	Sf5910 = MSbot2(2)*(42*MSbot2(2) - 27*MStop2(2)) - 
     &   76*MStop2(2)**2

	Sf3466 = MSbot2(2)*(42*MSbot2(2) - 27*MStop2(1)) - 
     &   76*MStop2(1)**2

	Sf8144 = 23*MStop2(2)**2 + 
     &   MSbot2(2)*(29*MSbot2(2) + 20*MStop2(2))

	Sf7935 = 23*MStop2(1)**2 + 
     &   MSbot2(2)*(29*MSbot2(2) + 20*MStop2(1))

	Sf6897 = 25*MStop2(2)**2 + 
     &   MSbot2(2)*(25*MSbot2(2) + 28*MStop2(2))

	Sf6780 = 16*MStop2(2)**2 + 
     &   MSbot2(2)*(25*MSbot2(2) + 7*MStop2(2))

	Sf7077 = MSbot2(2)*(25*MSbot2(2) - 14*MStop2(2)) - 
     &   45*MStop2(2)**2

	Sf4701 = 25*MStop2(1)**2 + 
     &   MSbot2(2)*(25*MSbot2(2) + 28*MStop2(1))

	Sf4284 = 16*MStop2(1)**2 + 
     &   MSbot2(2)*(25*MSbot2(2) + 7*MStop2(1))

	Sf4646 = MSbot2(2)*(25*MSbot2(2) - 14*MStop2(1)) - 
     &   45*MStop2(1)**2

	Sf5364 = 11*MStop2(2)**2 + 
     &   MSbot2(2)*(23*MSbot2(2) + 5*MStop2(2))

	Sf2953 = 11*MStop2(1)**2 + 
     &   MSbot2(2)*(23*MSbot2(2) + 5*MStop2(1))

	Sf5158 = MSbot2(2)*(17*MSbot2(2) - 3*MStop2(2)) + MStop2(2)**2

	Sf2751 = MSbot2(2)*(17*MSbot2(2) - 3*MStop2(1)) + MStop2(1)**2

	Sf6191 = MSbot2(2)*(15*MSbot2(2) - 2*MStop2(2)) + 3*MStop2(2)**2

	Sf6109 = MSbot2(2)*(15*MSbot2(2) - 2*MStop2(1)) + 3*MStop2(1)**2

	Sf6989 = 11*MStop2(2)**2 + 
     &   MSbot2(2)*(13*MSbot2(2) + 8*MStop2(2))

	Sf8178 = MSbot2(2)*(13*MSbot2(2) - 4*MStop2(2)) - 
     &   37*MStop2(2)**2

	Sf8143 = MSbot2(2)*(13*MSbot2(2) - 6*MStop2(2)) + 7*MStop2(2)**2

	Sf4280 = 11*MStop2(1)**2 + 
     &   MSbot2(2)*(13*MSbot2(2) + 8*MStop2(1))

	Sf7969 = MSbot2(2)*(13*MSbot2(2) - 4*MStop2(1)) - 
     &   37*MStop2(1)**2

	Sf7934 = MSbot2(2)*(13*MSbot2(2) - 6*MStop2(1)) + 7*MStop2(1)**2

	Sf6873 = -7*MStop2(2)**2 + 
     &   MSbot2(2)*(11*MSbot2(2) + 11*MStop2(2))

	Sf5398 = 23*MStop2(2)**2 + 
     &   MSbot2(2)*(11*MSbot2(2) + 5*MStop2(2))

	Sf6195 = MSbot2(2)*(11*MSbot2(2) - MStop2(2)) + 6*MStop2(2)**2

	Sf6759 = -7*MStop2(1)**2 + 
     &   MSbot2(2)*(11*MSbot2(2) + 11*MStop2(1))

	Sf2990 = 23*MStop2(1)**2 + 
     &   MSbot2(2)*(11*MSbot2(2) + 5*MStop2(1))

	Sf4002 = MSbot2(2)*(11*MSbot2(2) - MStop2(1)) + 6*MStop2(1)**2

	Sf4916 = -5*MStop2(2)**2 + MSbot2(2)*(9*MSbot2(2) + 2*MStop2(2))

	Sf8124 = MSbot2(2)*(9*MSbot2(2) - 473*MStop2(2)) + 
     &   360*MStop2(2)**2

	Sf2345 = -5*MStop2(1)**2 + MSbot2(2)*(9*MSbot2(2) + 2*MStop2(1))

	Sf7915 = MSbot2(2)*(9*MSbot2(2) - 473*MStop2(1)) + 
     &   360*MStop2(1)**2

	Sf8148 = MSbot2(2)*(8*MSbot2(2) - 30*MStop2(2)) + 
     &   26*MStop2(2)**2

	Sf7939 = MSbot2(2)*(8*MSbot2(2) - 30*MStop2(1)) + 
     &   26*MStop2(1)**2

	Sf5086 = 11*MStop2(2)**2 + 
     &   MSbot2(2)*(7*MSbot2(2) + 26*MStop2(2))

	Sf8081 = 7*MStop2(2)**2 + MSbot2(2)*(7*MSbot2(2) + 12*MStop2(2))

	Sf5958 = 7*MStop2(2)**2 + MSbot2(2)*(7*MSbot2(2) + 10*MStop2(2))

	Sf7098 = MSbot2(2)*(7*MSbot2(2) - 11*MStop2(2)) - 
     &   11*MStop2(2)**2

	Sf2646 = 11*MStop2(1)**2 + 
     &   MSbot2(2)*(7*MSbot2(2) + 26*MStop2(1))

	Sf7872 = 7*MStop2(1)**2 + MSbot2(2)*(7*MSbot2(2) + 12*MStop2(1))

	Sf3533 = 7*MStop2(1)**2 + MSbot2(2)*(7*MSbot2(2) + 10*MStop2(1))

	Sf4667 = MSbot2(2)*(7*MSbot2(2) - 11*MStop2(1)) - 
     &   11*MStop2(1)**2

	Sf5966 = 5*MStop2(2)**2 + MSbot2(2)*(6*MSbot2(2) + 9*MStop2(2))

	Sf3541 = 5*MStop2(1)**2 + MSbot2(2)*(6*MSbot2(2) + 9*MStop2(1))

	Sf6898 = 5*MStop2(2)**2 + MSbot2(2)*(5*MSbot2(2) + 42*MStop2(2))

	Sf5417 = MSbot2(2)*(5*MSbot2(2) - 22*MStop2(2)) + 5*MStop2(2)**2

	Sf4702 = 5*MStop2(1)**2 + MSbot2(2)*(5*MSbot2(2) + 42*MStop2(1))

	Sf3017 = MSbot2(2)*(5*MSbot2(2) - 22*MStop2(1)) + 5*MStop2(1)**2

	Sf6776 = 15*MStop2(2)**2 + 
     &   MSbot2(2)*(3*MSbot2(2) + 14*MStop2(2))

	Sf4913 = MStop2(2)**2 + MSbot2(2)*(3*MSbot2(2) + 12*MStop2(2))

	Sf6203 = 3*MStop2(2)**2 + MSbot2(2)*(3*MSbot2(2) + 2*MStop2(2))

	Sf4909 = MSbot2(2)*(3*MSbot2(2) - 2*MStop2(2)) - MStop2(2)**2

	Sf6694 = 15*MStop2(1)**2 + 
     &   MSbot2(2)*(3*MSbot2(2) + 14*MStop2(1))

	Sf2342 = MStop2(1)**2 + MSbot2(2)*(3*MSbot2(2) + 12*MStop2(1))

	Sf4137 = 3*MStop2(1)**2 + MSbot2(2)*(3*MSbot2(2) + 2*MStop2(1))

	Sf2337 = MSbot2(2)*(3*MSbot2(2) - 2*MStop2(1)) - MStop2(1)**2

	Sf6215 = 9*MStop2(2)**2 + MSbot2(2)*(MSbot2(2) + 22*MStop2(2))

	Sf5922 = MStop2(2)**2 + MSbot2(2)*(MSbot2(2) + 3*MStop2(2))

	Sf5978 = MSbot2(2)*(MSbot2(2) - 13*MStop2(2)) + MStop2(2)**2

	Sf4154 = 9*MStop2(1)**2 + MSbot2(2)*(MSbot2(2) + 22*MStop2(1))

	Sf3478 = MStop2(1)**2 + MSbot2(2)*(MSbot2(2) + 3*MStop2(1))

	Sf3553 = MSbot2(2)*(MSbot2(2) - 13*MStop2(1)) + MStop2(1)**2

	Sf5257 = MSbot2(1)**2*(28*MSbot2(1) - 33*MStop2(2)) + 
     &   5*MStop2(2)**3

	Sf2801 = MSbot2(1)**2*(28*MSbot2(1) - 33*MStop2(1)) + 
     &   5*MStop2(1)**3

	Sf5293 = 28*MStop2(2)**3 + 
     &   MSbot2(1)*(5*MSbot2(1)**2 - 33*MStop2(2)**2)

	Sf2853 = 28*MStop2(1)**3 + 
     &   MSbot2(1)*(5*MSbot2(1)**2 - 33*MStop2(1)**2)

	Sf8032 = 13854*MStop2(2)**2 + 
     &   MSbot2(1)*(13125*MSbot2(1) + 18434*MStop2(2))

	Sf7823 = 13854*MStop2(1)**2 + 
     &   MSbot2(1)*(13125*MSbot2(1) + 18434*MStop2(1))

	Sf8055 = 11490*MStop2(2)**2 + 
     &   MSbot2(1)*(11004*MSbot2(1) + 15221*MStop2(2))

	Sf7846 = 11490*MStop2(1)**2 + 
     &   MSbot2(1)*(11004*MSbot2(1) + 15221*MStop2(1))

	Sf8047 = 2752*MStop2(2)**2 + 
     &   MSbot2(1)*(2509*MSbot2(1) + 3362*MStop2(2))

	Sf7838 = 2752*MStop2(1)**2 + 
     &   MSbot2(1)*(2509*MSbot2(1) + 3362*MStop2(1))

	Sf8296 = 1952*MStop2(2)**2 + 
     &   MSbot2(1)*(2271*MSbot2(1) + 2595*MStop2(2))

	Sf8272 = 1952*MStop2(1)**2 + 
     &   MSbot2(1)*(2271*MSbot2(1) + 2595*MStop2(1))

	Sf8026 = 853*MStop2(2)**2 + 
     &   MSbot2(1)*(2203*MSbot2(1) + 534*MStop2(2))

	Sf7817 = 853*MStop2(1)**2 + 
     &   MSbot2(1)*(2203*MSbot2(1) + 534*MStop2(1))

	Sf5535 = MSbot2(1)*(1907*MSbot2(1) - 80*MStop2(2)) + 
     &   1907*MStop2(2)**2

	Sf3174 = MSbot2(1)*(1907*MSbot2(1) - 80*MStop2(1)) + 
     &   1907*MStop2(1)**2

	Sf8446 = 2487*MStop2(2)**2 + 
     &   MSbot2(1)*(1763*MSbot2(1) + 2568*MStop2(2))

	Sf8422 = 2487*MStop2(1)**2 + 
     &   MSbot2(1)*(1763*MSbot2(1) + 2568*MStop2(1))

	Sf8295 = -799*MStop2(2)**2 + 
     &   MSbot2(1)*(1743*MSbot2(1) + 282*MStop2(2))

	Sf8271 = -799*MStop2(1)**2 + 
     &   MSbot2(1)*(1743*MSbot2(1) + 282*MStop2(1))

	Sf8039 = 2738*MStop2(2)**2 + 
     &   MSbot2(1)*(1685*MSbot2(1) + 501*MStop2(2))

	Sf7830 = 2738*MStop2(1)**2 + 
     &   MSbot2(1)*(1685*MSbot2(1) + 501*MStop2(1))

	Sf8054 = 1781*MStop2(2)**2 + 
     &   MSbot2(1)*(1214*MSbot2(1) + 1128*MStop2(2))

	Sf7845 = 1781*MStop2(1)**2 + 
     &   MSbot2(1)*(1214*MSbot2(1) + 1128*MStop2(1))

	Sf7822 = -6116*MStop2(1)**2 + 
     &   MSbot2(1)*(1174*MSbot2(1) + 7470*MStop2(1))

	Sf5533 = 913*MStop2(2)**2 + 
     &   MSbot2(1)*(913*MSbot2(1) + 1096*MStop2(2))

	Sf3172 = 913*MStop2(1)**2 + 
     &   MSbot2(1)*(913*MSbot2(1) + 1096*MStop2(1))

	Sf5536 = 886*MStop2(2)**2 + 
     &   MSbot2(1)*(886*MSbot2(1) + 149*MStop2(2))

	Sf3175 = 886*MStop2(1)**2 + 
     &   MSbot2(1)*(886*MSbot2(1) + 149*MStop2(1))

	Sf8025 = 731*MStop2(2)**2 + 
     &   MSbot2(1)*(857*MSbot2(1) + 876*MStop2(2))

	Sf7816 = 731*MStop2(1)**2 + 
     &   MSbot2(1)*(857*MSbot2(1) + 876*MStop2(1))

	Sf8038 = 697*MStop2(2)**2 + 
     &   MSbot2(1)*(697*MSbot2(1) + 622*MStop2(2))

	Sf7829 = 697*MStop2(1)**2 + 
     &   MSbot2(1)*(697*MSbot2(1) + 622*MStop2(1))

	Sf8031 = -3058*MStop2(2)**2 + 
     &   MSbot2(1)*(587*MSbot2(1) + 3735*MStop2(2))

	Sf8019 = MSbot2(1)*(565*MSbot2(1) - 360*MStop2(2)) + 
     &   133*MStop2(2)**2

	Sf7799 = MSbot2(1)*(565*MSbot2(1) - 360*MStop2(1)) + 
     &   133*MStop2(1)**2

	Sf8164 = 718*MStop2(2)**2 + 
     &   MSbot2(1)*(520*MSbot2(1) + 366*MStop2(2))

	Sf7955 = 718*MStop2(1)**2 + 
     &   MSbot2(1)*(520*MSbot2(1) + 366*MStop2(1))

	Sf8163 = 589*MStop2(2)**2 + 
     &   MSbot2(1)*(508*MSbot2(1) + 507*MStop2(2))

	Sf7954 = 589*MStop2(1)**2 + 
     &   MSbot2(1)*(508*MSbot2(1) + 507*MStop2(1))

	Sf8290 = -147*MStop2(2)**2 + 
     &   MSbot2(1)*(465*MSbot2(1) + 46*MStop2(2))

	Sf8266 = -147*MStop2(1)**2 + 
     &   MSbot2(1)*(465*MSbot2(1) + 46*MStop2(1))

	Sf8688 = 499*MStop2(2)**2 + 
     &   MSbot2(1)*(463*MSbot2(1) + 462*MStop2(2))

	Sf8648 = 499*MStop2(1)**2 + 
     &   MSbot2(1)*(463*MSbot2(1) + 462*MStop2(1))

	Sf8020 = 365*MStop2(2)**2 + 
     &   MSbot2(1)*(437*MSbot2(1) + 326*MStop2(2))

	Sf7800 = 365*MStop2(1)**2 + 
     &   MSbot2(1)*(437*MSbot2(1) + 326*MStop2(1))

	Sf8302 = 520*MStop2(2)**2 + 
     &   MSbot2(1)*(385*MSbot2(1) + 219*MStop2(2))

	Sf8278 = 520*MStop2(1)**2 + 
     &   MSbot2(1)*(385*MSbot2(1) + 219*MStop2(1))

	Sf5305 = 371*MStop2(2)**2 + 
     &   MSbot2(1)*(371*MSbot2(1) + 414*MStop2(2))

	Sf2866 = 371*MStop2(1)**2 + 
     &   MSbot2(1)*(371*MSbot2(1) + 414*MStop2(1))

	Sf8248 = 229*MStop2(2)**2 + 
     &   MSbot2(1)*(315*MSbot2(1) + 258*MStop2(2))

	Sf8285 = 297*MStop2(2)**2 + 
     &   MSbot2(1)*(315*MSbot2(1) + 244*MStop2(2))

	Sf8236 = 229*MStop2(1)**2 + 
     &   MSbot2(1)*(315*MSbot2(1) + 258*MStop2(1))

	Sf8261 = 297*MStop2(1)**2 + 
     &   MSbot2(1)*(315*MSbot2(1) + 244*MStop2(1))

	Sf5470 = 307*MStop2(2)**2 + 
     &   MSbot2(1)*(307*MSbot2(1) + 418*MStop2(2))

	Sf3073 = 307*MStop2(1)**2 + 
     &   MSbot2(1)*(307*MSbot2(1) + 418*MStop2(1))

	Sf8007 = MSbot2(1)*(299*MSbot2(1) - 220*MStop2(2)) + 
     &   137*MStop2(2)**2

	Sf7804 = MSbot2(1)*(299*MSbot2(1) - 220*MStop2(1)) + 
     &   137*MStop2(1)**2

	Sf5473 = 274*MStop2(2)**2 + 
     &   MSbot2(1)*(274*MSbot2(1) + 23*MStop2(2))

	Sf3076 = 274*MStop2(1)**2 + 
     &   MSbot2(1)*(274*MSbot2(1) + 23*MStop2(1))

	Sf8687 = 310*MStop2(2)**2 + 
     &   MSbot2(1)*(266*MSbot2(1) + 204*MStop2(2))

	Sf8647 = 310*MStop2(1)**2 + 
     &   MSbot2(1)*(266*MSbot2(1) + 204*MStop2(1))

	Sf5276 = 130*MStop2(2)**2 + 
     &   MSbot2(1)*(261*MSbot2(1) + 325*MStop2(2))

	Sf2824 = 130*MStop2(1)**2 + 
     &   MSbot2(1)*(261*MSbot2(1) + 325*MStop2(1))

	Sf8445 = MSbot2(1)*(259*MSbot2(1) - 1092*MStop2(2)) - 
     &   2985*MStop2(2)**2

	Sf8421 = MSbot2(1)*(259*MSbot2(1) - 1092*MStop2(1)) - 
     &   2985*MStop2(1)**2

	Sf8013 = MSbot2(1)*(257*MSbot2(1) - 30*MStop2(2)) + 
     &   95*MStop2(2)**2

	Sf7810 = MSbot2(1)*(257*MSbot2(1) - 30*MStop2(1)) + 
     &   95*MStop2(1)**2

	Sf8440 = MSbot2(1)*(237*MSbot2(1) - 28*MStop2(2)) - 
     &   285*MStop2(2)**2

	Sf8416 = MSbot2(1)*(237*MSbot2(1) - 28*MStop2(1)) - 
     &   285*MStop2(1)**2

	Sf8008 = 235*MStop2(2)**2 + 
     &   MSbot2(1)*(235*MSbot2(1) + 874*MStop2(2))

	Sf5270 = 81*MStop2(2)**2 + 
     &   MSbot2(1)*(235*MSbot2(1) + 174*MStop2(2))

	Sf7805 = 235*MStop2(1)**2 + 
     &   MSbot2(1)*(235*MSbot2(1) + 874*MStop2(1))

	Sf2818 = 81*MStop2(1)**2 + 
     &   MSbot2(1)*(235*MSbot2(1) + 174*MStop2(1))

	Sf6640 = 305*MStop2(2)**2 + 
     &   MSbot2(1)*(227*MSbot2(1) + 196*MStop2(2))

	Sf4581 = 305*MStop2(1)**2 + 
     &   MSbot2(1)*(227*MSbot2(1) + 196*MStop2(1))

	Sf5241 = 247*MStop2(2)**2 + 
     &   MSbot2(1)*(222*MSbot2(1) + 351*MStop2(2))

	Sf2772 = 247*MStop2(1)**2 + 
     &   MSbot2(1)*(222*MSbot2(1) + 351*MStop2(1))

	Sf5269 = 51*MStop2(2)**2 + 
     &   MSbot2(1)*(208*MSbot2(1) + 171*MStop2(2))

	Sf2817 = 51*MStop2(1)**2 + 
     &   MSbot2(1)*(208*MSbot2(1) + 171*MStop2(1))

	Sf8170 = 360*MStop2(2)**2 + 
     &   MSbot2(1)*(193*MSbot2(1) + 249*MStop2(2))

	Sf7961 = 360*MStop2(1)**2 + 
     &   MSbot2(1)*(193*MSbot2(1) + 249*MStop2(1))

	Sf5783 = 229*MStop2(2)**2 + 
     &   MSbot2(1)*(191*MSbot2(1) + 180*MStop2(2))

	Sf3786 = 229*MStop2(1)**2 + 
     &   MSbot2(1)*(191*MSbot2(1) + 180*MStop2(1))

	Sf6601 = 91*MStop2(2)**2 + 
     &   MSbot2(1)*(178*MSbot2(1) + 103*MStop2(2))

	Sf6475 = 91*MStop2(1)**2 + 
     &   MSbot2(1)*(178*MSbot2(1) + 103*MStop2(1))

	Sf6602 = 148*MStop2(2)**2 + 
     &   MSbot2(1)*(165*MSbot2(1) + 203*MStop2(2))

	Sf6476 = 148*MStop2(1)**2 + 
     &   MSbot2(1)*(165*MSbot2(1) + 203*MStop2(1))

	Sf5310 = 163*MStop2(2)**2 + 
     &   MSbot2(1)*(163*MSbot2(1) + 180*MStop2(2))

	Sf2889 = 163*MStop2(1)**2 + 
     &   MSbot2(1)*(163*MSbot2(1) + 180*MStop2(1))

	Sf7043 = 147*MStop2(2)**2 + 
     &   MSbot2(1)*(160*MSbot2(1) + 209*MStop2(2))

	Sf4536 = 147*MStop2(1)**2 + 
     &   MSbot2(1)*(160*MSbot2(1) + 209*MStop2(1))

	Sf5275 = 81*MStop2(2)**2 + 
     &   MSbot2(1)*(158*MSbot2(1) + 181*MStop2(2))

	Sf2823 = 81*MStop2(1)**2 + 
     &   MSbot2(1)*(158*MSbot2(1) + 181*MStop2(1))

	Sf6469 = -75*MStop2(1)**2 + 
     &   MSbot2(1)*(135*MSbot2(1) + 42*MStop2(1))

	Sf5304 = 121*MStop2(2)**2 + 
     &   MSbot2(1)*(123*MSbot2(1) + 100*MStop2(2))

	Sf2865 = 121*MStop2(1)**2 + 
     &   MSbot2(1)*(123*MSbot2(1) + 100*MStop2(1))

	Sf8301 = 107*MStop2(2)**2 + 
     &   MSbot2(1)*(107*MSbot2(1) + 74*MStop2(2))

	Sf8277 = 107*MStop2(1)**2 + 
     &   MSbot2(1)*(107*MSbot2(1) + 74*MStop2(1))

	Sf5240 = 134*MStop2(2)**2 + 
     &   MSbot2(1)*(103*MSbot2(1) + 183*MStop2(2))

	Sf2771 = 134*MStop2(1)**2 + 
     &   MSbot2(1)*(103*MSbot2(1) + 183*MStop2(1))

	Sf8247 = MSbot2(1)*(99*MSbot2(1) - 54*MStop2(2)) - 
     &   35*MStop2(2)**2

	Sf8235 = MSbot2(1)*(99*MSbot2(1) - 54*MStop2(1)) - 
     &   35*MStop2(1)**2

	Sf7044 = 190*MStop2(2)**2 + 
     &   MSbot2(1)*(97*MSbot2(1) + 109*MStop2(2))

	Sf4537 = 190*MStop2(1)**2 + 
     &   MSbot2(1)*(97*MSbot2(1) + 109*MStop2(1))

	Sf6639 = 85*MStop2(2)**2 + 
     &   MSbot2(1)*(85*MSbot2(1) + 98*MStop2(2))

	Sf4580 = 85*MStop2(1)**2 + 
     &   MSbot2(1)*(85*MSbot2(1) + 98*MStop2(1))

	Sf5235 = 235*MStop2(2)**2 + 
     &   MSbot2(1)*(81*MSbot2(1) + 174*MStop2(2))

	Sf2766 = 235*MStop2(1)**2 + 
     &   MSbot2(1)*(81*MSbot2(1) + 174*MStop2(1))

	Sf5469 = 79*MStop2(2)**2 + 
     &   MSbot2(1)*(79*MSbot2(1) + 64*MStop2(2))

	Sf3072 = 79*MStop2(1)**2 + 
     &   MSbot2(1)*(79*MSbot2(1) + 64*MStop2(1))

	Sf7038 = 69*MStop2(2)**2 + 
     &   MSbot2(1)*(76*MSbot2(1) + 89*MStop2(2))

	Sf4531 = 69*MStop2(1)**2 + 
     &   MSbot2(1)*(76*MSbot2(1) + 89*MStop2(1))

	Sf8014 = 59*MStop2(2)**2 + 
     &   MSbot2(1)*(73*MSbot2(1) + 60*MStop2(2))

	Sf7811 = 59*MStop2(1)**2 + 
     &   MSbot2(1)*(73*MSbot2(1) + 60*MStop2(1))

	Sf8169 = -189*MStop2(2)**2 + 
     &   MSbot2(1)*(71*MSbot2(1) + 108*MStop2(2))

	Sf8284 = 71*MStop2(2)**2 + 
     &   MSbot2(1)*(71*MSbot2(1) + 50*MStop2(2))

	Sf7960 = -189*MStop2(1)**2 + 
     &   MSbot2(1)*(71*MSbot2(1) + 108*MStop2(1))

	Sf8260 = 71*MStop2(1)**2 + 
     &   MSbot2(1)*(71*MSbot2(1) + 50*MStop2(1))

	Sf6596 = 76*MStop2(2)**2 + 
     &   MSbot2(1)*(69*MSbot2(1) + 89*MStop2(2))

	Sf6470 = 76*MStop2(1)**2 + 
     &   MSbot2(1)*(69*MSbot2(1) + 89*MStop2(1))

	Sf8071 = 72*MStop2(2)**2 + 
     &   MSbot2(1)*(63*MSbot2(1) + 76*MStop2(2))

	Sf7862 = 72*MStop2(1)**2 + 
     &   MSbot2(1)*(63*MSbot2(1) + 76*MStop2(1))

	Sf5472 = MSbot2(1)*(61*MSbot2(1) - 16*MStop2(2)) + 
     &   61*MStop2(2)**2

	Sf3075 = MSbot2(1)*(61*MSbot2(1) - 16*MStop2(1)) + 
     &   61*MStop2(1)**2

	Sf8242 = MSbot2(1)*(55*MSbot2(1) - 14*MStop2(2)) - 
     &   13*MStop2(2)**2

	Sf8230 = MSbot2(1)*(55*MSbot2(1) - 14*MStop2(1)) - 
     &   13*MStop2(1)**2

	Sf5477 = 53*MStop2(2)**2 + 
     &   MSbot2(1)*(54*MSbot2(1) + 73*MStop2(2))

	Sf3080 = 53*MStop2(1)**2 + 
     &   MSbot2(1)*(54*MSbot2(1) + 73*MStop2(1))

	Sf5234 = 208*MStop2(2)**2 + 
     &   MSbot2(1)*(51*MSbot2(1) + 171*MStop2(2))

	Sf2765 = 208*MStop2(1)**2 + 
     &   MSbot2(1)*(51*MSbot2(1) + 171*MStop2(1))

	Sf5540 = MSbot2(1)*(46*MSbot2(1) - 9*MStop2(2)) - 
     &   52*MStop2(2)**2

	Sf3180 = MSbot2(1)*(46*MSbot2(1) - 9*MStop2(1)) - 
     &   52*MStop2(1)**2

	Sf6595 = -25*MStop2(2)**2 + 
     &   MSbot2(1)*(45*MSbot2(1) + 14*MStop2(2))

	Sf8156 = 33*MStop2(2)**2 + 
     &   MSbot2(1)*(43*MSbot2(1) + 12*MStop2(2))

	Sf7947 = 33*MStop2(1)**2 + 
     &   MSbot2(1)*(43*MSbot2(1) + 12*MStop2(1))

	Sf5476 = MSbot2(1)*(42*MSbot2(1) - 27*MStop2(2)) - 
     &   76*MStop2(2)**2

	Sf3079 = MSbot2(1)*(42*MSbot2(1) - 27*MStop2(1)) - 
     &   76*MStop2(1)**2

	Sf8066 = 23*MStop2(2)**2 + 
     &   MSbot2(1)*(29*MSbot2(1) + 20*MStop2(2))

	Sf7857 = 23*MStop2(1)**2 + 
     &   MSbot2(1)*(29*MSbot2(1) + 20*MStop2(1))

	Sf6646 = 25*MStop2(2)**2 + 
     &   MSbot2(1)*(25*MSbot2(1) + 28*MStop2(2))

	Sf6525 = 16*MStop2(2)**2 + 
     &   MSbot2(1)*(25*MSbot2(1) + 7*MStop2(2))

	Sf7037 = MSbot2(1)*(25*MSbot2(1) - 14*MStop2(2)) - 
     &   45*MStop2(2)**2

	Sf4604 = 25*MStop2(1)**2 + 
     &   MSbot2(1)*(25*MSbot2(1) + 28*MStop2(1))

	Sf4276 = 16*MStop2(1)**2 + 
     &   MSbot2(1)*(25*MSbot2(1) + 7*MStop2(1))

	Sf4530 = MSbot2(1)*(25*MSbot2(1) - 14*MStop2(1)) - 
     &   45*MStop2(1)**2

	Sf5256 = 11*MStop2(2)**2 + 
     &   MSbot2(1)*(23*MSbot2(1) + 5*MStop2(2))

	Sf2800 = 11*MStop2(1)**2 + 
     &   MSbot2(1)*(23*MSbot2(1) + 5*MStop2(1))

	Sf5040 = MSbot2(1)*(17*MSbot2(1) - 3*MStop2(2)) + MStop2(2)**2

	Sf2575 = MSbot2(1)*(17*MSbot2(1) - 3*MStop2(1)) + MStop2(1)**2

	Sf5813 = MSbot2(1)*(15*MSbot2(1) - 2*MStop2(2)) + 3*MStop2(2)**2

	Sf5720 = MSbot2(1)*(15*MSbot2(1) - 2*MStop2(1)) + 3*MStop2(1)**2

	Sf7058 = MSbot2(1)*(14*MSbot2(1) - 22*MStop2(2)) - 
     &   22*MStop2(2)**2

	Sf6985 = 11*MStop2(2)**2 + 
     &   MSbot2(1)*(13*MSbot2(1) + 8*MStop2(2))

	Sf8158 = MSbot2(1)*(13*MSbot2(1) - 4*MStop2(2)) - 
     &   37*MStop2(2)**2

	Sf8065 = MSbot2(1)*(13*MSbot2(1) - 6*MStop2(2)) + 7*MStop2(2)**2

	Sf4272 = 11*MStop2(1)**2 + 
     &   MSbot2(1)*(13*MSbot2(1) + 8*MStop2(1))

	Sf7949 = MSbot2(1)*(13*MSbot2(1) - 4*MStop2(1)) - 
     &   37*MStop2(1)**2

	Sf7856 = MSbot2(1)*(13*MSbot2(1) - 6*MStop2(1)) + 7*MStop2(1)**2

	Sf6622 = -7*MStop2(2)**2 + 
     &   MSbot2(1)*(11*MSbot2(1) + 11*MStop2(2))

	Sf5292 = 23*MStop2(2)**2 + 
     &   MSbot2(1)*(11*MSbot2(1) + 5*MStop2(2))

	Sf5817 = MSbot2(1)*(11*MSbot2(1) - MStop2(2)) + 6*MStop2(2)**2

	Sf6504 = -7*MStop2(1)**2 + 
     &   MSbot2(1)*(11*MSbot2(1) + 11*MStop2(1))

	Sf2852 = 23*MStop2(1)**2 + 
     &   MSbot2(1)*(11*MSbot2(1) + 5*MStop2(1))

	Sf3995 = MSbot2(1)*(11*MSbot2(1) - MStop2(1)) + 6*MStop2(1)**2

	Sf4825 = -5*MStop2(2)**2 + MSbot2(1)*(9*MSbot2(1) + 2*MStop2(2))

	Sf4902 = MSbot2(1)*(9*MSbot2(1) - 6*MStop2(2)) - 3*MStop2(2)**2

	Sf8046 = MSbot2(1)*(9*MSbot2(1) - 473*MStop2(2)) + 
     &   360*MStop2(2)**2

	Sf2332 = -5*MStop2(1)**2 + MSbot2(1)*(9*MSbot2(1) + 2*MStop2(1))

	Sf7837 = MSbot2(1)*(9*MSbot2(1) - 473*MStop2(1)) + 
     &   360*MStop2(1)**2

	Sf8070 = MSbot2(1)*(8*MSbot2(1) - 30*MStop2(2)) + 
     &   26*MStop2(2)**2

	Sf7861 = MSbot2(1)*(8*MSbot2(1) - 30*MStop2(1)) + 
     &   26*MStop2(1)**2

	Sf4960 = 11*MStop2(2)**2 + 
     &   MSbot2(1)*(7*MSbot2(1) + 26*MStop2(2))

	Sf8002 = 7*MStop2(2)**2 + MSbot2(1)*(7*MSbot2(1) + 12*MStop2(2))

	Sf5532 = 7*MStop2(2)**2 + MSbot2(1)*(7*MSbot2(1) + 10*MStop2(2))

	Sf7064 = MSbot2(1)*(7*MSbot2(1) - 11*MStop2(2)) - 
     &   11*MStop2(2)**2

	Sf2432 = 11*MStop2(1)**2 + 
     &   MSbot2(1)*(7*MSbot2(1) + 26*MStop2(1))

	Sf7794 = 7*MStop2(1)**2 + MSbot2(1)*(7*MSbot2(1) + 12*MStop2(1))

	Sf3171 = 7*MStop2(1)**2 + MSbot2(1)*(7*MSbot2(1) + 10*MStop2(1))

	Sf4562 = MSbot2(1)*(7*MSbot2(1) - 11*MStop2(1)) - 
     &   11*MStop2(1)**2

	Sf5541 = 5*MStop2(2)**2 + MSbot2(1)*(6*MSbot2(1) + 9*MStop2(2))

	Sf3179 = 5*MStop2(1)**2 + MSbot2(1)*(6*MSbot2(1) + 9*MStop2(1))

	Sf6647 = 5*MStop2(2)**2 + MSbot2(1)*(5*MSbot2(1) + 42*MStop2(2))

	Sf5311 = MSbot2(1)*(5*MSbot2(1) - 22*MStop2(2)) + 5*MStop2(2)**2

	Sf4605 = 5*MStop2(1)**2 + MSbot2(1)*(5*MSbot2(1) + 42*MStop2(1))

	Sf2890 = MSbot2(1)*(5*MSbot2(1) - 22*MStop2(1)) + 5*MStop2(1)**2

	Sf6521 = 15*MStop2(2)**2 + 
     &   MSbot2(1)*(3*MSbot2(1) + 14*MStop2(2))

	Sf4822 = MStop2(2)**2 + MSbot2(1)*(3*MSbot2(1) + 12*MStop2(2))

	Sf5825 = 3*MStop2(2)**2 + MSbot2(1)*(3*MSbot2(1) + 2*MStop2(2))

	Sf6421 = 15*MStop2(1)**2 + 
     &   MSbot2(1)*(3*MSbot2(1) + 14*MStop2(1))

	Sf2329 = MStop2(1)**2 + MSbot2(1)*(3*MSbot2(1) + 12*MStop2(1))

	Sf4018 = 3*MStop2(1)**2 + MSbot2(1)*(3*MSbot2(1) + 2*MStop2(1))

	Sf2323 = MSbot2(1)*(3*MSbot2(1) - 2*MStop2(1)) - MStop2(1)**2

	Sf5838 = 9*MStop2(2)**2 + MSbot2(1)*(MSbot2(1) + 22*MStop2(2))

	Sf5489 = MStop2(2)**2 + MSbot2(1)*(MSbot2(1) + 3*MStop2(2))

	Sf5552 = MSbot2(1)*(MSbot2(1) - 13*MStop2(2)) + MStop2(2)**2

	Sf4043 = 9*MStop2(1)**2 + MSbot2(1)*(MSbot2(1) + 22*MStop2(1))

	Sf3095 = MStop2(1)**2 + MSbot2(1)*(MSbot2(1) + 3*MStop2(1))

	Sf3197 = MSbot2(1)*(MSbot2(1) - 13*MStop2(1)) + MStop2(1)**2

	Sf7519 = MSbot2(2)**2 + MSbot2(1)*(MSbot2(1) + 12*MSbot2(2))

	Sf8531 = -9*MStop2(3)**4 + 
     &   MZ2*(MZ2*Sf8530 + 10*MStop2(3)**2*MStop2(4))

	Sf7712 = -9*MStau2(3)**4 + 
     &   MZ2*(MZ2*Sf7711 + 10*MStau2(3)**2*MStau2(4))

	Sf7520 = -9*MSbot2(3)**4 + 
     &   MZ2*(MZ2*Sf7519 + 10*MSbot2(3)**2*MSbot2(4))

	Sf8553 = 65 + 7*USbot22**2

	Sf8613 = 5 + 4*USbot22**2

	Sf7474 = 23 - 14*UStop22**2

	Sf7482 = 14 - 5*UStop22**2

	Sf564 = 23*MSl2Diff(1,1)**6 + 
     &   MW2**2*(MW2**2*Sf563 + 10*Sf562*MSl2Diff(1,1)**2)

	Sf1467 = 23*MSl2Diff(2,1)**6 + 
     &   MW2**2*(MW2**2*Sf1466 + 10*Sf1465*MSl2Diff(2,1)**2)

	Sf10282 = -3*MHiggs2(1)*MStop2(2)**2 + 
     &   MZ2*(MZ2**2 - 2*Ss44*MStop2(2))

	Sf10281 = -3*MHiggs2(1)*MStop2(1)**2 + 
     &   MZ2*(MZ2**2 - 2*Ss44*MStop2(1))

	Sf9805 = MZ2*(279*MZ2 - 464*MStau2(2)) - 672*MStau2(2)**2

	Sf9673 = MZ2*(279*MZ2 - 464*MStau2(1)) - 672*MStau2(1)**2

	Sf10017 = MZ2*(69*MZ2 - 1882*MStau2(2)) - 780*MStau2(2)**2

	Sf9922 = MZ2*(69*MZ2 - 1882*MStau2(1)) - 780*MStau2(1)**2

	Sf9813 = MZ2*(37*MZ2 - 272*MStau2(2)) - 118*MStau2(2)**2

	Sf9679 = MZ2*(37*MZ2 - 272*MStau2(1)) - 118*MStau2(1)**2

	Sf10037 = MZ2*(29*MZ2 - 8*MStau2(2)) + 28*MStau2(2)**2

	Sf9955 = MZ2*(29*MZ2 - 8*MStau2(1)) + 28*MStau2(1)**2

	Sf9982 = 36*MStau2(2)**2 + MZ2*(22*MZ2 + 55*MStau2(2))

	Sf9886 = 36*MStau2(1)**2 + MZ2*(22*MZ2 + 55*MStau2(1))

	Sf9699 = 6*MStau2(2)**2 + MZ2*(18*MZ2 + 37*MStau2(2))

	Sf9409 = 6*MStau2(1)**2 + MZ2*(18*MZ2 + 37*MStau2(1))

	Sf9826 = MZ2*(9*MZ2 - 188*MStau2(2)) - 240*MStau2(2)**2

	Sf9643 = MZ2*(9*MZ2 - 188*MStau2(1)) - 240*MStau2(1)**2

	Sf10038 = MZ2*(8*MZ2 - 61*MStau2(2)) + 30*MStau2(2)**2

	Sf9959 = MZ2*(8*MZ2 - 61*MStau2(1)) + 30*MStau2(1)**2

	Sf9834 = MZ2*(7*MZ2 - 68*MStau2(2)) - 34*MStau2(2)**2

	Sf9648 = MZ2*(7*MZ2 - 68*MStau2(1)) - 34*MStau2(1)**2

	Sf9981 = 4*MStau2(2)**2 + MZ2*(5*MZ2 + MStau2(2))

	Sf9988 = MZ2*(5*MZ2 - 42*MStau2(2)) + 20*MStau2(2)**2

	Sf9882 = 4*MStau2(1)**2 + MZ2*(5*MZ2 + MStau2(1))

	Sf9876 = MZ2*(5*MZ2 - 42*MStau2(1)) + 20*MStau2(1)**2

	Sf8522 = 9*MStop2(3)**2 + MZ2*(4*MZ2 - 13*MStop2(4))

	Sf7703 = 9*MStau2(3)**2 + MZ2*(4*MZ2 - 13*MStau2(4))

	Sf7511 = 9*MSbot2(3)**2 + MZ2*(4*MZ2 - 13*MSbot2(4))

	Sf9989 = MZ2*(2*MZ2 - 3*MStau2(2)) + MStau2(2)**2

	Sf9872 = MZ2*(2*MZ2 - 3*MStau2(1)) + MStau2(1)**2

	Sf8813 = MStop2(3)**2 + MZ2*(MZ2 - 2*MStop2(4))

	Sf8783 = MStau2(3)**2 + MZ2*(MZ2 - 2*MStau2(4))

	Sf9695 = MZ2*(MZ2 - 2*MStau2(2)) - 10*MStau2(2)**2

	Sf9696 = MZ2*(MZ2 - 3*MStau2(2)) + 4*MStau2(2)**2

	Sf9400 = MZ2*(MZ2 - 2*MStau2(1)) - 10*MStau2(1)**2

	Sf9399 = MZ2*(MZ2 - 3*MStau2(1)) + 4*MStau2(1)**2

	Sf8765 = MSbot2(3)**2 + MZ2*(MZ2 - 2*MSbot2(4))

	Sf10471 = -MZ2**2 + MHiggs2(2)*(3*MZ2 + 2*MStop2(2))

	Sf10438 = -MZ2**2 + MHiggs2(2)*(3*MZ2 + 2*MStop2(1))

	Sf9350 = -MZ2**2 + MHiggs2(2)*(3*MZ2 + 2*MStau2(2))

	Sf9319 = -MZ2**2 + MHiggs2(2)*(3*MZ2 + 2*MStau2(1))

	Sf8958 = -MZ2**2 + MHiggs2(2)*(3*MZ2 + 2*MSbot2(2))

	Sf8906 = -MZ2**2 + MHiggs2(2)*(3*MZ2 + 2*MSbot2(1))

	Sf10371 = -MZ2**2 + MHiggs2(1)*(3*MZ2 + 2*MStop2(2))

	Sf10344 = -MZ2**2 + MHiggs2(1)*(3*MZ2 + 2*MStop2(1))

	Sf10086 = -MZ2**2 + MHiggs2(1)*(3*MZ2 + 2*MStau2(2))

	Sf10057 = -MZ2**2 + MHiggs2(1)*(3*MZ2 + 2*MStau2(1))

	Sf9148 = -MZ2**2 + MHiggs2(1)*(3*MZ2 + 2*MSbot2(2))

	Sf9110 = -MZ2**2 + MHiggs2(1)*(3*MZ2 + 2*MSbot2(1))

	Sf6095 = 180*MSbot2(2)**2 + MW2*(89*MW2 + 148*MSbot2(2))

	Sf5706 = 180*MSbot2(1)**2 + MW2*(89*MW2 + 148*MSbot2(1))

	Sf3664 = 63*MSbot2(2)**2 + MW2*(35*MW2 + 65*MSbot2(2))

	Sf3977 = 45*MSbot2(2)**2 + MW2*(35*MW2 + 52*MSbot2(2))

	Sf3384 = 63*MSbot2(1)**2 + MW2*(35*MW2 + 65*MSbot2(1))

	Sf3880 = 45*MSbot2(1)**2 + MW2*(35*MW2 + 52*MSbot2(1))

	Sf1809 = 75*MStau2(2)**2 + MW2*(29*MW2 + 53*MStau2(2))

	Sf997 = 75*MStau2(1)**2 + MW2*(29*MW2 + 53*MStau2(1))

	Sf5570 = 5*MStop2(2)**2 + MW2*(28*MW2 + 17*MStop2(2))

	Sf3247 = 5*MStop2(1)**2 + MW2*(28*MW2 + 17*MStop2(1))

	Sf6309 = 30*MSbot2(1)**2 + MW2*(24*MW2 + 23*MSbot2(1))

	Sf5508 = 5*MStop2(2)**2 + MW2*(16*MW2 + 17*MStop2(2))

	Sf3129 = 5*MStop2(1)**2 + MW2*(16*MW2 + 17*MStop2(1))

	Sf7756 = MW2*(16*MW2 - 32*MZ2) + 13*MZ2**2

	Sf3569 = 7*MSbot2(2)**2 + MW2*(15*MW2 + 13*MSbot2(2))

	Sf3217 = 7*MSbot2(1)**2 + MW2*(15*MW2 + 13*MSbot2(1))

	Sf8597 = MW2*(8*MW2 - 16*MZ2) + 7*MZ2**2

	Sf4521 = 3*MSbot2(2)**2 + MW2*(7*MW2 + 3*MSbot2(2))

	Sf4417 = 3*MSbot2(1)**2 + MW2*(7*MW2 + 3*MSbot2(1))

	Sf8535 = MW2*(2*MW2 - 3*MZ2) + MZ2**2

	Sf5848 = MW2*(MW2 - 88*MStop2(2)) + 40*MStop2(2)**2

	Sf4064 = MW2*(MW2 - 88*MStop2(1)) + 40*MStop2(1)**2

	Sf850 = 14*MW2**2 + MHiggs2(4)*(-21*MW2 + 9*MHiggs2(4))

	Sf2209 = 4*MW2**2 + MHiggs2(4)*(-7*MW2 + 6*MHiggs2(4))

	Sf1041 = 4*MW2**2 + MHiggs2(4)*(-8*MW2 + 6*MHiggs2(4))

	Sf216 = 8*MW2**2 + MHiggs2(4)*(-13*MW2 + 5*MHiggs2(4))

	Sf1874 = MHiggs2(4)*(4*MHiggs2(4) - 15*MStau2(2)) + 
     &   8*MStau2(2)**2

	Sf1123 = MHiggs2(4)*(4*MHiggs2(4) - 15*MStau2(1)) + 
     &   8*MStau2(1)**2

	Sf237 = 2*MW2**2 + MHiggs2(4)*(-4*MW2 + 3*MHiggs2(4))

	Sf1008 = 3*MW2**2 + MHiggs2(4)*(-4*MW2 + 2*MHiggs2(4))

	Sf4870 = MSq2Diff(2,2)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 2*MStop2(2))

	Sf2273 = MSq2Diff(2,1)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 2*MStop2(1))

	Sf4835 = MSq2Diff(1,2)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 2*MStop2(2))

	Sf2210 = MSq2Diff(1,1)**2 + 
     &   MHiggs2(4)*(MHiggs2(4) - 2*MStop2(1))

	Sf1923 = MW2*(-2*MW2 - 2*MHiggs2(4)) + MHiggs2(4)**2

	Sf217 = MW2*(2*MW2 - 3*MHiggs2(4)) + MHiggs2(4)**2

	Sf545 = MW2*(6*MW2 - 5*MHiggs2(4)) + MHiggs2(4)**2

	Sf9450 = -935*MZ2**2 + MHiggs2(2)*(2272*MZ2 + 1305*MHiggs2(2))

	Sf9444 = -96*MZ2**2 + MHiggs2(2)*(-1306*MZ2 + 420*MHiggs2(2))

	Sf9434 = 249*MZ2**2 + MHiggs2(2)*(1600*MZ2 + 349*MHiggs2(2))

	Sf9556 = -401*MZ2**2 + MHiggs2(2)*(-130*MZ2 + 245*MHiggs2(2))

	Sf9713 = -208*MStau2(2)**2 + 
     &   MHiggs2(2)*(224*MHiggs2(2) + 6*MStau2(2))

	Sf9425 = -208*MStau2(1)**2 + 
     &   MHiggs2(2)*(224*MHiggs2(2) + 6*MStau2(1))

	Sf9458 = -107*MZ2**2 + MHiggs2(2)*(289*MZ2 + 161*MHiggs2(2))

	Sf9463 = 80*MZ2**2 + MHiggs2(2)*(71*MZ2 + 82*MHiggs2(2))

	Sf9798 = MHiggs2(2)*(80*MHiggs2(2) - 298*MStau2(2)) - 
     &   244*MStau2(2)**2

	Sf9660 = MHiggs2(2)*(80*MHiggs2(2) - 298*MStau2(1)) - 
     &   244*MStau2(1)**2

	Sf10009 = MHiggs2(2)*(69*MHiggs2(2) - 1418*MStau2(2)) - 
     &   620*MStau2(2)**2

	Sf9912 = MHiggs2(2)*(69*MHiggs2(2) - 1418*MStau2(1)) - 
     &   620*MStau2(1)**2

	Sf9522 = -88*MZ2**2 + MHiggs2(2)*(-22*MZ2 + 52*MHiggs2(2))

	Sf9508 = -43*MZ2**2 + MHiggs2(2)*(-16*MZ2 + 33*MHiggs2(2))

	Sf9452 = 21*MZ2**2 + MHiggs2(2)*(-491*MZ2 + 33*MHiggs2(2))

	Sf9436 = 76*MZ2**2 + MHiggs2(2)*(-262*MZ2 + 32*MHiggs2(2))

	Sf9769 = MHiggs2(2)*(28*MHiggs2(2) - 33*MStau2(2)) - 
     &   80*MStau2(2)**2

	Sf9538 = MHiggs2(2)*(28*MHiggs2(2) - 33*MStau2(1)) - 
     &   80*MStau2(1)**2

	Sf9567 = 622*MZ2**2 + MHiggs2(2)*(-502*MZ2 + 20*MHiggs2(2))

	Sf9819 = -20*MStau2(2)**2 + 
     &   MHiggs2(2)*(16*MHiggs2(2) + 34*MStau2(2))

	Sf9630 = -20*MStau2(1)**2 + 
     &   MHiggs2(2)*(16*MHiggs2(2) + 34*MStau2(1))

	Sf9495 = -28*MZ2**2 + MHiggs2(2)*(-2*MZ2 + 16*MHiggs2(2))

	Sf10045 = 700*MStau2(2)**2 + 
     &   MHiggs2(2)*(15*MHiggs2(2) + 766*MStau2(2))

	Sf9970 = 700*MStau2(1)**2 + 
     &   MHiggs2(2)*(15*MHiggs2(2) + 766*MStau2(1))

	Sf9743 = 14*MStau2(2)**2 + 
     &   MHiggs2(2)*(14*MHiggs2(2) + 15*MStau2(2))

	Sf9485 = 14*MStau2(1)**2 + 
     &   MHiggs2(2)*(14*MHiggs2(2) + 15*MStau2(1))

	Sf9500 = -6*MZ2**2 + MHiggs2(2)*(-8*MZ2 + 11*MHiggs2(2))

	Sf9578 = 591*MZ2**2 + MHiggs2(2)*(-587*MZ2 + 10*MHiggs2(2))

	Sf9995 = 68*MStau2(2)**2 + 
     &   MHiggs2(2)*(9*MHiggs2(2) + 8*MStau2(2))

	Sf9896 = 68*MStau2(1)**2 + 
     &   MHiggs2(2)*(9*MHiggs2(2) + 8*MStau2(1))

	Sf9510 = -15*MZ2**2 + MHiggs2(2)*(11*MZ2 + 9*MHiggs2(2))

	Sf9502 = -18*MZ2**2 + MHiggs2(2)*(4*MZ2 + 9*MHiggs2(2))

	Sf9560 = 129*MZ2**2 + MHiggs2(2)*(-74*MZ2 + 7*MHiggs2(2))

	Sf9581 = 353*MZ2**2 + MHiggs2(2)*(-354*MZ2 + 7*MHiggs2(2))

	Sf9742 = MHiggs2(2)*(5*MHiggs2(2) - MStau2(2)) - 2*MStau2(2)**2

	Sf9481 = MHiggs2(2)*(5*MHiggs2(2) - MStau2(1)) - 2*MStau2(1)**2

	Sf9493 = -3*MZ2**2 + MHiggs2(2)*(-16*MZ2 + 5*MHiggs2(2))

	Sf9573 = 267*MZ2**2 + MHiggs2(2)*(-252*MZ2 + 5*MHiggs2(2))

	Sf9543 = 5*MZ2**2 + MHiggs2(2)*(-13*MZ2 + 4*MHiggs2(2))

	Sf9517 = -28*MZ2**2 + MHiggs2(2)*(-46*MZ2 + 4*MHiggs2(2))

	Sf9442 = 117*MZ2**2 + MHiggs2(2)*(-328*MZ2 + 4*MHiggs2(2))

	Sf9687 = 32*MStau2(2)**2 + 
     &   MHiggs2(2)*(2*MHiggs2(2) + 35*MStau2(2))

	Sf9393 = 32*MStau2(1)**2 + 
     &   MHiggs2(2)*(2*MHiggs2(2) + 35*MStau2(1))

	Sf9552 = 29*MZ2**2 + MHiggs2(2)*(-5*MZ2 + 2*MHiggs2(2))

	Sf9583 = 85*MZ2**2 + MHiggs2(2)*(-85*MZ2 + 2*MHiggs2(2))

	Sf9548 = MHiggs2(2)**2 + MZ2*(61*MZ2 + 32*MHiggs2(2))

	Sf9515 = MHiggs2(2)**2 + MZ2*(-3*MZ2 + 8*MHiggs2(2))

	Sf8932 = MHiggs2(2)**2 + MZ2*(-5*MZ2 + 5*MHiggs2(2))

	Sf8934 = MZ2*(5*MZ2 - 5*MHiggs2(2)) + MHiggs2(2)**2

	Sf9134 = MHiggs2(1)**2 + MZ2*(-5*MZ2 + 5*MHiggs2(1))

	Sf9136 = MZ2*(5*MZ2 - 5*MHiggs2(1)) + MHiggs2(1)**2

	Sf793 = 320 + (483 - 18*C2B)*C2B

	Sf790 = 320 + (-546 - 24*C2B)*C2B

	Sf752 = 230 + (-657 - 57*C2B)*C2B

	Sf779 = 89 + C2B*(19 + 2*C2B)

	Sf778 = 76 + (-61 - C2B)*C2B

	Sf780 = 69 + (42 - C2B)*C2B

	Sf791 = 66 + C2B*(-51 + 10*C2B)

	Sf8823 = 470 - 411*USbot22**2

	Sf769 = 45 + (249 - 22*C2B)*C2B

	Sf757 = 35 + (105 - 6*C2B)*C2B

	Sf792 = 35 + C2B*(36 + 2*C2B)

	Sf743 = 31 + (330 - 45*C2B)*C2B

	Sf754 = 28 + (99 - 15*C2B)*C2B

	Sf767 = 25 + (-140 - 14*C2B)*C2B

	Sf967 = 16 + (-51 - 3*C2B)*C2B

	Sf5665 = 16 + (-69 - 3*C2B)*C2B

	Sf8841 = 14 - 63*USbot22**2

	Sf768 = 14 + (163 - 22*C2B)*C2B

	Sf8829 = -16 + 27*USbot22**2

	Sf8832 = 11 - 27*USbot22**2

	Sf7502 = 11 + SW2*(-80 + 160*SW2)

	Sf4775 = 9 + SW2*(-6 + 4*USbot12**2)

	Sf7201 = 9 - 4*SW2*(1 + 2*UStop22**2)

	Sf8615 = 22 - 15*USbot22**2

	Sf8617 = 7 + 15*USbot22**2

	Sf3864 = 6 + (53 - 3*C2B)*C2B

	Sf5680 = 6 + (-2 - 3*C2B)*C2B

	Sf7458 = 5 - 20*UStau22**2

	Sf760 = 5 + (-21 - 6*C2B)*C2B

	Sf7411 = -3 + 6*UStop22**2

	Sf7414 = -6 + 9*UStop22**2

	Sf7424 = 3 - 3*UStop22**2

	Sf7427 = -Sf7411

	Sf7430 = 3 - 9*UStop22**2

	Sf7330 = -1 + 4*UStau22**2

	Sf7333 = -2 + 5*UStau22**2

	Sf7339 = 3 - 3*UStau22**2

	Sf7342 = 3 - 4*UStau22**2

	Sf7345 = 3 - 5*UStau22**2

	Sf7362 = 6 - 3*USbot22**2

	Sf7382 = 3 - 3*USbot22**2

	Sf1757 = 3 + SW2*(-6 + 4*UStau12**2)

	Sf5690 = 3 + (-19 - 6*C2B)*C2B

	Sf740 = 2 + (57 - 12*C2B)*C2B

	Sf983 = 2 + (4 - C2B)*C2B

	Sf742 = 2 + C2B*(-96 + 3*C2B)

	Sf7445 = 6 - 5*UStau22**2

	Sf7452 = 4 - 3*UStau22**2

	Sf8590 = -Sf8556

	Sf8569 = -2 + 3*USbot22**2

	Sf8623 = -5 + 6*USbot22**2

	Sf8592 = 1 + 3*USbot22**2

	Sf8571 = 1 - 3*USbot22**2

	Sf8625 = 1 - 6*USbot22**2

	Sf7722 = 1 + (2 - 4*SW2)*SW2

	Sf7497 = 1 + SW2*(-4 + 8*SW2)

	Sf3349 = 1 + C2B*(4 + 2*C2B)

	Sf987 = 1 + (-4 - 2*C2B)*C2B

	Sf8752 = -(MZ2*Sf8751) - 
     &   97*MStop2(3)*UStop12*UStop21*UStop22

	Sf8533 = Ss211**2*UStop21**2 + 16*SW2**2*UStop22**2

	Sf8532 = -4*SW2*UStop22**2 + Ss211*(1 - UStop22**2)

	Sf7396 = Sf7479

	Sf10268 = UStop12*UStop21 + UStop11*UStop22

	Sf6945 = Re(Abot)*UStop21 + MTop*UStop22

	Sf10242 = 3*CW*MTop*Ss966*UStop21 - Ss967*UStop22

	Sf6975 = reMUE*SB*UStop21 + 2*CB*MTop*UStop22

	Sf6978 = 2*MTop*UStop21 + Re(Atop)*UStop22

	Sf7485 = 9 - 9*UStop22**2

	Sf7477 = 6 - 6*UStop22**2

	Sf7469 = Sf7424

	Sf7465 = 2 - 2*UStop22**2

	Sf8549 = MStop2(2)*UStop21**2 + MStop2(1)*UStop22**2

	Sf7991 = 16*SW2**2*UStop12**2 + Ss211**2*UStop22**2

	Sf7990 = -4*SW2 + 3*UStop22**2

	Sf7352 = 4 - 3*UStop22**2

	Sf4227 = Re(Abot)*UStop11 + MTop*UStop12

	Sf10222 = 3*CW*MTop*Ss966*UStop11 - Ss967*UStop12

	Sf4256 = reMUE*SB*UStop11 + 2*CB*MTop*UStop12

	Sf4261 = 2*MTop*UStop11 + Re(Atop)*UStop12

	Sf8608 = -(Sf8607*MStau2(1)) + 
     &   7*MStau2(2)*UStau12*UStau21*UStau22

	Sf7714 = Ss922**2*UStau21**2 + 4*SW2**2*UStau22**2

	Sf10219 = Ss922*UStau21**2 + 2*SW2*UStau22**2

	Sf7713 = -2*SW2*UStau22**2 + Ss922*(1 - UStau22**2)

	Sf7248 = Ss443**2*UStau21**2 + 4*MZ2**2*SW2**2*UStau22**2

	Sf65 = -2*MZ2*SW2*UStau22**2 + Ss443*(1 - UStau22**2)

	Sf7323 = Sf7449

	Sf9311 = UStau12*UStau21 + UStau11*UStau22

	Sf1730 = ML*reMUE*SA*UStau21 - 2*MZ2*SAB*SB*SW2*UStau22

	Sf1729 = CB*reMUE*UStau21 + 2*ML*SB*UStau22

	Sf1733 = 2*ML*UStau21 + Re(Atau)*UStau22

	Sf7447 = 4 - 4*UStau22**2

	Sf7440 = Sf7339

	Sf7437 = 2 - 2*UStau22**2

	Sf7580 = 4*SW2**2*UStau12**2 + Ss922**2*UStau22**2

	Sf9299 = 2*SW2*UStau12**2 + Ss922*UStau22**2

	Sf7579 = -2*SW2 + UStau22**2

	Sf7244 = 4*MZ2**2*SW2**2*UStau12**2 + Ss443**2*UStau22**2

	Sf43 = -2*MZ2*SW2*UStau12**2 + Ss443*UStau22**2

	Sf7308 = 2 - UStau22**2

	Sf941 = ML*reMUE*SA*UStau11 - 2*MZ2*SAB*SB*SW2*UStau12

	Sf940 = CB*reMUE*UStau11 + 2*ML*SB*UStau12

	Sf945 = 2*ML*UStau11 + Re(Atau)*UStau12

	Sf7522 = Ss879**2*USbot21**2 + 4*SW2**2*USbot22**2

	Sf9287 = Ss879*USbot21**2 + 2*SW2*USbot22**2

	Sf7521 = -2*SW2*USbot22**2 + Ss879*(1 - USbot22**2)

	Sf8867 = 57 + 52*USbot21**2

	Sf7304 = Sf8571

	Sf8896 = USbot12*USbot21 + USbot11*USbot22

	Sf141 = 3*MBot*reMUE*SA*USbot21 - 2*MZ2*SAB*SB*SW2*USbot22

	Sf140 = CB*reMUE*USbot21 + 2*MBot*SB*USbot22

	Sf144 = 2*MBot*USbot21 + Re(Abot)*USbot22

	Sf8868 = -5*Sf8867 + 34*USbot12**2

	Sf8621 = 7 - 7*USbot22**2

	Sf8845 = 2 - 2*USbot22**2

	Sf7466 = Sf8845

	Sf8583 = -1 + USbot22**2

	Sf7503 = 4*SW2**2*USbot12**2 + Ss879**2*USbot22**2

	Sf8882 = 2*SW2*USbot12**2 + Ss879*USbot22**2

	Sf7501 = -2*SW2 + 3*USbot22**2

	Sf7298 = Sf8569

	Sf129 = 3*MBot*reMUE*SA*USbot11 - 2*MZ2*SAB*SB*SW2*USbot12

	Sf128 = CB*reMUE*USbot11 + 2*MBot*SB*USbot12

	Sf133 = 2*MBot*USbot11 + Re(Abot)*USbot12

	Sf8524 = 9*MStop2(3)**2 + 4*MZ2*MStop2(4)

	Sf8523 = 6*MStop2(3)**2 + MZ2*MStop2(4)

	Sf8661 = 2*MStop2(3)**2 - MZ2*MStop2(4)

	Sf8525 = MStop2(3)**2 - MZ2*MStop2(4)

	Sf8528 = MStop2(1)**2 + MStop2(2)**2

	Sf10298 = 2*MStop2(1) + MStop2(2)

	Sf10300 = MStop2(1) - MStop2(3)

	Sf10288 = MStop2(1) + 2*MStop2(2)

	Sf10290 = -MStop2(2) - MStop2(3)

	Sf109 = MStop(1) + MStop(2)

	Sf107 = MStop(1) - MStop(2)

	Sf7705 = 9*MStau2(3)**2 + 4*MZ2*MStau2(4)

	Sf7704 = 6*MStau2(3)**2 + MZ2*MStau2(4)

	Sf7723 = 2*MStau2(3)**2 - MZ2*MStau2(4)

	Sf7706 = MStau2(3)**2 - MZ2*MStau2(4)

	Sf9878 = 2*Sf9872*MStau2(1)**3 - Sf9873*MStau2(1)*MStau2(2)

	Sf7709 = MStau2(1)**2 + MStau2(2)**2

	Sf8579 = Sf8543*MStau2(1) - 
     &   2*MStau2(2)*UStau12*UStau21*UStau22

	Sf10114 = 2*MStau2(1) + MStau2(2)

	Sf10116 = MStau2(1) - MStau2(3)

	Sf10106 = MStau2(1) + 2*MStau2(2)

	Sf10108 = -MStau2(2) - MStau2(3)

	Sf57 = MStau(1) + MStau(2)

	Sf55 = MStau(1) - MStau(2)

	Sf6058 = 3*MSq2Diff(2,2)**2 - MW2*MStop2(2)

	Sf3727 = 3*MSq2Diff(2,1)**2 - MW2*MStop2(1)

	Sf3456 = 3*MSq2Diff(1,1)**2 - MW2*MStop2(1)

	Sf1544 = MSneu2(1)**3 - MStau2(2)**3

	Sf670 = MSneu2(1)**3 - MStau2(1)**3

	Sf1462 = MSneu2(1)**2 + MStau2(2)**2

	Sf559 = MSneu2(1)**2 + MStau2(1)**2

	Sf1610 = Sf768*MSneu2(1) + Sf769*MStau2(2)

	Sf770 = Sf768*MSneu2(1) + Sf769*MStau2(1)

	Sf7622 = 1411*MSneu2(1) + 1465*MStau2(2)

	Sf7564 = 1411*MSneu2(1) + 1465*MStau2(1)

	Sf7591 = 826*MSneu2(1) - 254*MStau2(2)

	Sf7533 = 826*MSneu2(1) - 254*MStau2(1)

	Sf7592 = 797*MSneu2(1) + 824*MStau2(2)

	Sf7534 = 797*MSneu2(1) + 824*MStau2(1)

	Sf7693 = 662*MSneu2(1) + 655*MStau2(2)

	Sf7679 = 662*MSneu2(1) + 655*MStau2(1)

	Sf7740 = 619*MSneu2(1) + 698*MStau2(2)

	Sf7716 = 619*MSneu2(1) + 698*MStau2(1)

	Sf7607 = 475*MSneu2(1) + 547*MStau2(2)

	Sf7549 = 475*MSneu2(1) + 547*MStau2(1)

	Sf1334 = 208*MSneu2(1) + 119*MStau2(2)

	Sf397 = 208*MSneu2(1) + 119*MStau2(1)

	Sf1333 = 181*MSneu2(1) + 101*MStau2(2)

	Sf396 = 181*MSneu2(1) + 101*MStau2(1)

	Sf7600 = 179*MSneu2(1) + 188*MStau2(2)

	Sf7542 = 179*MSneu2(1) + 188*MStau2(1)

	Sf7662 = 128*MSneu2(1) + 115*MStau2(2)

	Sf7650 = 128*MSneu2(1) + 115*MStau2(1)

	Sf1304 = 119*MSneu2(1) + 208*MStau2(2)

	Sf353 = 119*MSneu2(1) + 208*MStau2(1)

	Sf8327 = 103*MSneu2(1) + 140*MStau2(2)

	Sf8320 = 103*MSneu2(1) + 140*MStau2(1)

	Sf1303 = 101*MSneu2(1) + 181*MStau2(2)

	Sf352 = 101*MSneu2(1) + 181*MStau2(1)

	Sf7739 = 92*MSneu2(1) + 367*MStau2(2)

	Sf7715 = 92*MSneu2(1) + 367*MStau2(1)

	Sf7621 = 89*MSneu2(1) + 413*MStau2(2)

	Sf7563 = 89*MSneu2(1) + 413*MStau2(1)

	Sf7692 = 79*MSneu2(1) - 52*MStau2(2)

	Sf7678 = 79*MSneu2(1) - 52*MStau2(1)

	Sf7657 = 77*MSneu2(1) + 85*MStau2(2)

	Sf1706 = 77*MSneu2(1) - 131*MStau2(2)

	Sf7645 = 77*MSneu2(1) + 85*MStau2(1)

	Sf907 = 77*MSneu2(1) - 131*MStau2(1)

	Sf1310 = 73*MSneu2(1) + 68*MStau2(2)

	Sf359 = 73*MSneu2(1) + 68*MStau2(1)

	Sf2141 = 71*MSneu2(1) + 79*MStau2(2)

	Sf2061 = 71*MSneu2(1) + 79*MStau2(1)

	Sf2146 = 48*MSneu2(1) + 55*MStau2(2)

	Sf2066 = 48*MSneu2(1) + 55*MStau2(1)

	Sf7599 = 46*MSneu2(1) - 26*MStau2(2)

	Sf7541 = 46*MSneu2(1) - 26*MStau2(1)

	Sf1848 = 45*MSneu2(1) + 49*MStau2(2)

	Sf2168 = 45*MSneu2(1) + 6*MStau2(2)

	Sf1069 = 45*MSneu2(1) + 49*MStau2(1)

	Sf2099 = 45*MSneu2(1) + 6*MStau2(1)

	Sf1339 = 41*MSneu2(1) + 10*MStau2(2)

	Sf402 = 41*MSneu2(1) + 10*MStau2(1)

	Sf7685 = 39*MSneu2(1) + 47*MStau2(2)

	Sf7671 = 39*MSneu2(1) + 47*MStau2(1)

	Sf7588 = 37*MSneu2(1) + MStau2(2)

	Sf1395 = 37*MSneu2(1) - 59*MStau2(2)

	Sf7530 = 37*MSneu2(1) + MStau2(1)

	Sf478 = 37*MSneu2(1) - 59*MStau2(1)

	Sf1326 = 35*MSneu2(1) + 26*MStau2(2)

	Sf389 = 35*MSneu2(1) + 26*MStau2(1)

	Sf7630 = 32*MSneu2(1) + 35*MStau2(2)

	Sf7572 = 32*MSneu2(1) + 35*MStau2(1)

	Sf1340 = 31*MSneu2(1) + 16*MStau2(2)

	Sf2140 = 31*MSneu2(1) - 4*MStau2(2)

	Sf403 = 31*MSneu2(1) + 16*MStau2(1)

	Sf2060 = 31*MSneu2(1) - 4*MStau2(1)

	Sf1354 = 26*MSneu2(1) + 35*MStau2(2)

	Sf431 = 26*MSneu2(1) + 35*MStau2(1)

	Sf1468 = 25*MSneu2(1) + 27*MStau2(2)

	Sf1961 = 25*MSneu2(1) + 22*MStau2(2)

	Sf565 = 25*MSneu2(1) + 27*MStau2(1)

	Sf1899 = 25*MSneu2(1) + 22*MStau2(1)

	Sf1640 = 23*MSneu2(1) + 27*MStau2(2)

	Sf1967 = 23*MSneu2(1) + 19*MStau2(2)

	Sf7586 = 23*MSneu2(1) + 11*MStau2(2)

	Sf815 = 23*MSneu2(1) + 27*MStau2(1)

	Sf1907 = 23*MSneu2(1) + 19*MStau2(1)

	Sf7528 = 23*MSneu2(1) + 11*MStau2(1)

	Sf1309 = 22*MSneu2(1) + 25*MStau2(2)

	Sf358 = 22*MSneu2(1) + 25*MStau2(1)

	Sf1388 = 21*MSneu2(1) + 19*MStau2(2)

	Sf1281 = 21*MSneu2(1) + 17*MStau2(2)

	Sf471 = 21*MSneu2(1) + 19*MStau2(1)

	Sf321 = 21*MSneu2(1) + 17*MStau2(1)

	Sf7656 = 19*MSneu2(1) + 28*MStau2(2)

	Sf7644 = 19*MSneu2(1) + 28*MStau2(1)

	Sf2163 = 15*MSneu2(1) + 2*MStau2(2)

	Sf2094 = 15*MSneu2(1) + 2*MStau2(1)

	Sf1409 = 14*MSneu2(1) - 32*MStau2(2)

	Sf492 = 14*MSneu2(1) - 32*MStau2(1)

	Sf1698 = 13*MSneu2(1) + 11*MStau2(2)

	Sf899 = 13*MSneu2(1) + 11*MStau2(1)

	Sf1681 = 11*MSneu2(1) + 7*MStau2(2)

	Sf869 = 11*MSneu2(1) + 7*MStau2(1)

	Sf2127 = 7*MSneu2(1) + 2*MStau2(2)

	Sf1387 = 7*MSneu2(1) - 16*MStau2(2)

	Sf2047 = 7*MSneu2(1) + 2*MStau2(1)

	Sf470 = 7*MSneu2(1) - 16*MStau2(1)

	Sf1699 = 5*MSneu2(1) + 76*MStau2(2)

	Sf1280 = 5*MSneu2(1) + MStau2(2)

	Sf900 = 5*MSneu2(1) + 76*MStau2(1)

	Sf320 = 5*MSneu2(1) + MStau2(1)

	Sf1379 = 4*MSneu2(1) + 6*MStau2(2)

	Sf1627 = 4*MSneu2(1) + 3*MStau2(2)

	Sf1325 = 4*MSneu2(1) + MStau2(2)

	Sf8326 = 4*MSneu2(1) - 5*MStau2(2)

	Sf462 = 4*MSneu2(1) + 6*MStau2(1)

	Sf802 = 4*MSneu2(1) + 3*MStau2(1)

	Sf388 = 4*MSneu2(1) + MStau2(1)

	Sf8319 = 4*MSneu2(1) - 5*MStau2(1)

	Sf1832 = 3*MSneu2(1) + 5*MStau2(2)

	Sf1495 = 3*MSneu2(1) + 4*MStau2(2)

	Sf1836 = 3*MSneu2(1) - MStau2(2)

	Sf1043 = 3*MSneu2(1) + 5*MStau2(1)

	Sf595 = 3*MSneu2(1) + 4*MStau2(1)

	Sf1048 = 3*MSneu2(1) - MStau2(1)

	Sf1628 = 2*MSneu2(1) + 15*MStau2(2)

	Sf1384 = 2*MSneu2(1) + 3*MStau2(2)

	Sf1372 = 2*MSneu2(1) + MStau2(2)

	Sf803 = 2*MSneu2(1) + 15*MStau2(1)

	Sf467 = 2*MSneu2(1) + 3*MStau2(1)

	Sf455 = 2*MSneu2(1) + MStau2(1)

	Sf7614 = MSneu2(1) + 5*MStau2(2)

	Sf1353 = MSneu2(1) + 4*MStau2(2)

	Sf1847 = MSneu2(1) + 2*MStau2(2)

	Sf47 = MSneu2(1) + MStau2(2)

	Sf1233 = MSneu2(1) - 3*MStau2(2)

	Sf7556 = MSneu2(1) + 5*MStau2(1)

	Sf430 = MSneu2(1) + 4*MStau2(1)

	Sf1068 = MSneu2(1) + 2*MStau2(1)

	Sf34 = MSneu2(1) + MStau2(1)

	Sf251 = MSneu2(1) - 3*MStau2(1)

	Sf1527 = Ss115*MSl2Diff(2,1)**2 - MW2*MHiggs2(4)*MStau2(2)

	Sf633 = Ss115*MSl2Diff(1,1)**2 - MW2*MHiggs2(4)*MStau2(1)

	Sf7513 = 9*MSbot2(3)**2 + 4*MZ2*MSbot2(4)

	Sf7512 = 6*MSbot2(3)**2 + MZ2*MSbot2(4)

	Sf8258 = 2*MSbot2(3)**2 - MZ2*MSbot2(4)

	Sf7514 = MSbot2(3)**2 - MZ2*MSbot2(4)

	Sf6221 = MSbot2(2)**3 - MStop2(2)**3

	Sf4162 = MSbot2(2)**3 - MStop2(1)**3

	Sf6299 = 7*MSbot2(2)**2 + 9*MStop2(2)**2

	Sf4005 = 7*MSbot2(2)**2 + 9*MStop2(1)**2

	Sf5456 = 3*MSbot2(2)**2 - MStop2(2)**2

	Sf3068 = 3*MSbot2(2)**2 - MStop2(1)**2

	Sf5141 = MSbot2(2)**2 + MStop2(2)**2

	Sf2733 = MSbot2(2)**2 + MStop2(1)**2

	Sf7387 = Sf7368*MSbot2(2) - MStop2(2)

	Sf8115 = 3113*MSbot2(2) + 3437*MStop2(2)

	Sf7906 = 3113*MSbot2(2) + 3437*MStop2(1)

	Sf8123 = 1789*MSbot2(2) + 1843*MStop2(2)

	Sf7914 = 1789*MSbot2(2) + 1843*MStop2(1)

	Sf8455 = 1546*MSbot2(2) + 1729*MStop2(2)

	Sf8431 = 1546*MSbot2(2) + 1729*MStop2(1)

	Sf8102 = 1009*MSbot2(2) + 955*MStop2(2)

	Sf7893 = 1009*MSbot2(2) + 955*MStop2(1)

	Sf8085 = 721*MSbot2(2) + 667*MStop2(2)

	Sf7876 = 721*MSbot2(2) + 667*MStop2(1)

	Sf8396 = 655*MSbot2(2) + 763*MStop2(2)

	Sf8372 = 655*MSbot2(2) + 763*MStop2(1)

	Sf8182 = 547*MSbot2(2) + 583*MStop2(2)

	Sf7973 = 547*MSbot2(2) + 583*MStop2(1)

	Sf8131 = 483*MSbot2(2) + 969*MStop2(2)

	Sf7922 = 483*MSbot2(2) + 969*MStop2(1)

	Sf8090 = 433*MSbot2(2) + 379*MStop2(2)

	Sf7881 = 433*MSbot2(2) + 379*MStop2(1)

	Sf8217 = 303*MSbot2(2) + 262*MStop2(2)

	Sf8205 = 303*MSbot2(2) + 262*MStop2(1)

	Sf8096 = 299*MSbot2(2) + 281*MStop2(2)

	Sf7887 = 299*MSbot2(2) + 281*MStop2(1)

	Sf8390 = 258*MSbot2(2) - 92*MStop2(2)

	Sf8366 = 258*MSbot2(2) - 92*MStop2(1)

	Sf8385 = 251*MSbot2(2) + 233*MStop2(2)

	Sf8361 = 251*MSbot2(2) + 233*MStop2(1)

	Sf8188 = 244*MSbot2(2) + 321*MStop2(2)

	Sf7979 = 244*MSbot2(2) + 321*MStop2(1)

	Sf8694 = 239*MSbot2(2) + 247*MStop2(2)

	Sf8654 = 239*MSbot2(2) + 247*MStop2(1)

	Sf5376 = 208*MSbot2(2) + 119*MStop2(2)

	Sf2968 = 208*MSbot2(2) + 119*MStop2(1)

	Sf8384 = 203*MSbot2(2) + 16*MStop2(2)

	Sf8360 = 203*MSbot2(2) + 16*MStop2(1)

	Sf5375 = 181*MSbot2(2) + 101*MStop2(2)

	Sf8107 = 181*MSbot2(2) + 46*MStop2(2)

	Sf2967 = 181*MSbot2(2) + 101*MStop2(1)

	Sf7898 = 181*MSbot2(2) + 46*MStop2(1)

	Sf8379 = 139*MSbot2(2) + 121*MStop2(2)

	Sf8355 = 139*MSbot2(2) + 121*MStop2(1)

	Sf5341 = 119*MSbot2(2) + 208*MStop2(2)

	Sf2930 = 119*MSbot2(2) + 208*MStop2(1)

	Sf5347 = 115*MSbot2(2) + 118*MStop2(2)

	Sf2936 = 115*MSbot2(2) + 118*MStop2(1)

	Sf8108 = 103*MSbot2(2) + 105*MStop2(2)

	Sf7899 = 103*MSbot2(2) + 105*MStop2(1)

	Sf5340 = 101*MSbot2(2) + 181*MStop2(2)

	Sf8147 = 101*MSbot2(2) + 107*MStop2(2)

	Sf2929 = 101*MSbot2(2) + 181*MStop2(1)

	Sf7938 = 101*MSbot2(2) + 107*MStop2(1)

	Sf5382 = 94*MSbot2(2) + 57*MStop2(2)

	Sf2974 = 94*MSbot2(2) + 57*MStop2(1)

	Sf8130 = 93*MSbot2(2) + 95*MStop2(2)

	Sf7921 = 93*MSbot2(2) + 95*MStop2(1)

	Sf8083 = 89*MSbot2(2) + 8*MStop2(2)

	Sf7874 = 89*MSbot2(2) + 8*MStop2(1)

	Sf8181 = 85*MSbot2(2) + 112*MStop2(2)

	Sf7972 = 85*MSbot2(2) + 112*MStop2(1)

	Sf7076 = 79*MSbot2(2) + 71*MStop2(2)

	Sf4645 = 79*MSbot2(2) + 71*MStop2(1)

	Sf6849 = 71*MSbot2(2) + 79*MStop2(2)

	Sf6735 = 71*MSbot2(2) + 79*MStop2(1)

	Sf8389 = 67*MSbot2(2) + 64*MStop2(2)

	Sf8078 = 67*MSbot2(2) + 40*MStop2(2)

	Sf8365 = 67*MSbot2(2) + 64*MStop2(1)

	Sf7869 = 67*MSbot2(2) + 40*MStop2(1)

	Sf8454 = 62*MSbot2(2) + 345*MStop2(2)

	Sf8430 = 62*MSbot2(2) + 345*MStop2(1)

	Sf8142 = 61*MSbot2(2) + 55*MStop2(2)

	Sf8212 = 61*MSbot2(2) + 39*MStop2(2)

	Sf7933 = 61*MSbot2(2) + 55*MStop2(1)

	Sf8200 = 61*MSbot2(2) + 39*MStop2(1)

	Sf7081 = 57*MSbot2(2) + 43*MStop2(2)

	Sf4650 = 57*MSbot2(2) + 43*MStop2(1)

	Sf8175 = 53*MSbot2(2) + 47*MStop2(2)

	Sf7966 = 53*MSbot2(2) + 47*MStop2(1)

	Sf6855 = 49*MSbot2(2) + 51*MStop2(2)

	Sf6741 = 49*MSbot2(2) + 51*MStop2(1)

	Sf8080 = 47*MSbot2(2) + 101*MStop2(2)

	Sf7871 = 47*MSbot2(2) + 101*MStop2(1)

	Sf6883 = 45*MSbot2(2) + 6*MStop2(2)

	Sf6769 = 45*MSbot2(2) + 6*MStop2(1)

	Sf8122 = 44*MSbot2(2) - 64*MStop2(2)

	Sf7913 = 44*MSbot2(2) - 64*MStop2(1)

	Sf6779 = 43*MSbot2(2) + 31*MStop2(2)

	Sf4283 = 43*MSbot2(2) + 31*MStop2(1)

	Sf8076 = 37*MSbot2(2) + 55*MStop2(2)

	Sf6888 = 37*MSbot2(2) + 43*MStop2(2)

	Sf8101 = 37*MSbot2(2) + 7*MStop2(2)

	Sf7867 = 37*MSbot2(2) + 55*MStop2(1)

	Sf4683 = 37*MSbot2(2) + 43*MStop2(1)

	Sf7892 = 37*MSbot2(2) + 7*MStop2(1)

	Sf7082 = 36*MSbot2(2) + 75*MStop2(2)

	Sf4651 = 36*MSbot2(2) + 75*MStop2(1)

	Sf5363 = 35*MSbot2(2) + 26*MStop2(2)

	Sf2952 = 35*MSbot2(2) + 26*MStop2(1)

	Sf8174 = 34*MSbot2(2) + 25*MStop2(2)

	Sf7965 = 34*MSbot2(2) + 25*MStop2(1)

	Sf6848 = 31*MSbot2(2) - 4*MStop2(2)

	Sf6734 = 31*MSbot2(2) - 4*MStop2(1)

	Sf6167 = 27*MSbot2(2) + 29*MStop2(2)

	Sf3919 = 27*MSbot2(2) + 29*MStop2(1)

	Sf5397 = 26*MSbot2(2) + 35*MStop2(2)

	Sf2989 = 26*MSbot2(2) + 35*MStop2(1)

	Sf8187 = 24*MSbot2(2) - 26*MStop2(2)

	Sf7978 = 24*MSbot2(2) - 26*MStop2(1)

	Sf5157 = 23*MSbot2(2) + 8*MStop2(2)

	Sf2750 = 23*MSbot2(2) + 8*MStop2(1)

	Sf8693 = 22*MSbot2(2) + 25*MStop2(2)

	Sf8653 = 22*MSbot2(2) + 25*MStop2(1)

	Sf8691 = 20*MSbot2(2) - 95*MStop2(2)

	Sf8651 = 20*MSbot2(2) - 95*MStop2(1)

	Sf6194 = 19*MSbot2(2) + 15*MStop2(2)

	Sf6854 = 19*MSbot2(2) + 8*MStop2(2)

	Sf4001 = 19*MSbot2(2) + 15*MStop2(1)

	Sf6740 = 19*MSbot2(2) + 8*MStop2(1)

	Sf5408 = 17*MSbot2(2) + 15*MStop2(2)

	Sf3000 = 17*MSbot2(2) + 15*MStop2(1)

	Sf6035 = 15*MSbot2(2) + 19*MStop2(2)

	Sf5909 = 15*MSbot2(2) + 14*MStop2(2)

	Sf6871 = 15*MSbot2(2) + 2*MStop2(2)

	Sf3694 = 15*MSbot2(2) + 19*MStop2(1)

	Sf3465 = 15*MSbot2(2) + 14*MStop2(1)

	Sf6757 = 15*MSbot2(2) + 2*MStop2(1)

	Sf8138 = 13*MSbot2(2) + 7*MStop2(2)

	Sf8211 = 13*MSbot2(2) - 6*MStop2(2)

	Sf7929 = 13*MSbot2(2) + 7*MStop2(1)

	Sf8199 = 13*MSbot2(2) - 6*MStop2(1)

	Sf5965 = 11*MSbot2(2) + 10*MStop2(2)

	Sf5455 = 11*MSbot2(2) + 9*MStop2(2)

	Sf6190 = 11*MSbot2(2) + 6*MStop2(2)

	Sf6193 = 11*MSbot2(2) + 5*MStop2(2)

	Sf3540 = 11*MSbot2(2) + 10*MStop2(1)

	Sf3067 = 11*MSbot2(2) + 9*MStop2(1)

	Sf6108 = 11*MSbot2(2) + 6*MStop2(1)

	Sf4000 = 11*MSbot2(2) + 5*MStop2(1)

	Sf6988 = 10*MSbot2(2) + 9*MStop2(2)

	Sf4279 = 10*MSbot2(2) + 9*MStop2(1)

	Sf5085 = 9*MSbot2(2) + 13*MStop2(2)

	Sf5346 = 9*MSbot2(2) + 11*MStop2(2)

	Sf6778 = 9*MSbot2(2) + 7*MStop2(2)

	Sf6189 = 9*MSbot2(2) + 2*MStop2(2)

	Sf2645 = 9*MSbot2(2) + 13*MStop2(1)

	Sf2935 = 9*MSbot2(2) + 11*MStop2(1)

	Sf4282 = 9*MSbot2(2) + 7*MStop2(1)

	Sf6107 = 9*MSbot2(2) + 2*MStop2(1)

	Sf8140 = 8*MSbot2(2) + 26*MStop2(2)

	Sf6298 = 8*MSbot2(2) + 9*MStop2(2)

	Sf8216 = 8*MSbot2(2) - 6*MStop2(2)

	Sf7931 = 8*MSbot2(2) + 26*MStop2(1)

	Sf4004 = 8*MSbot2(2) + 9*MStop2(1)

	Sf8204 = 8*MSbot2(2) - 6*MStop2(1)

	Sf5381 = 7*MSbot2(2) + 5*MStop2(2)

	Sf4907 = 7*MSbot2(2) - 4*MStop2(2)

	Sf2973 = 7*MSbot2(2) + 5*MStop2(1)

	Sf2335 = 7*MSbot2(2) - 4*MStop2(1)

	Sf7109 = 6*MSbot2(2) + 45*MStop2(2)

	Sf6954 = 6*MSbot2(2) + 4*MStop2(2)

	Sf4678 = 6*MSbot2(2) + 45*MStop2(1)

	Sf4236 = 6*MSbot2(2) + 4*MStop2(1)

	Sf6775 = 5*MSbot2(2) + 8*MStop2(2)

	Sf4908 = 5*MSbot2(2) + 4*MStop2(2)

	Sf4887 = 5*MSbot2(2) + MStop2(2)

	Sf6334 = 5*MSbot2(2) - 18*MStop2(2)

	Sf6693 = 5*MSbot2(2) + 8*MStop2(1)

	Sf2336 = 5*MSbot2(2) + 4*MStop2(1)

	Sf2298 = 5*MSbot2(2) + MStop2(1)

	Sf4185 = 5*MSbot2(2) - 18*MStop2(1)

	Sf6987 = 4*MSbot2(2) + 3*MStop2(2)

	Sf5362 = 4*MSbot2(2) + MStop2(2)

	Sf8450 = 4*MSbot2(2) - 19*MStop2(2)

	Sf7075 = 4*MSbot2(2) - 31*MStop2(2)

	Sf4278 = 4*MSbot2(2) + 3*MStop2(1)

	Sf2951 = 4*MSbot2(2) + MStop2(1)

	Sf8426 = 4*MSbot2(2) - 19*MStop2(1)

	Sf4644 = 4*MSbot2(2) - 31*MStop2(1)

	Sf6034 = 3*MSbot2(2) + 7*MStop2(2)

	Sf6872 = 3*MSbot2(2) + 4*MStop2(2)

	Sf4912 = 3*MSbot2(2) - MStop2(2)

	Sf5964 = 3*MSbot2(2) - 8*MStop2(2)

	Sf3693 = 3*MSbot2(2) + 7*MStop2(1)

	Sf6758 = 3*MSbot2(2) + 4*MStop2(1)

	Sf2341 = 3*MSbot2(2) - MStop2(1)

	Sf3539 = 3*MSbot2(2) - 8*MStop2(1)

	Sf7097 = 2*MSbot2(2) + 15*MStop2(2)

	Sf5140 = 2*MSbot2(2) + 3*MStop2(2)

	Sf4889 = 2*MSbot2(2) + MStop2(2)

	Sf4877 = 2*MSbot2(2) - MStop2(2)

	Sf4666 = 2*MSbot2(2) + 15*MStop2(1)

	Sf2732 = 2*MSbot2(2) + 3*MStop2(1)

	Sf2300 = 2*MSbot2(2) + MStop2(1)

	Sf2280 = 2*MSbot2(2) - MStop2(1)

	Sf7421 = MSbot2(2) - Sf7402*MStop2(2)

	Sf5396 = MSbot2(2) + 4*MStop2(2)

	Sf4868 = MSbot2(2) + 3*MStop2(2)

	Sf4879 = MSbot2(2) + 2*MStop2(2)

	Sf99 = MSbot2(2) + MStop2(2)

	Sf8095 = MSbot2(2) - 5*MStop2(2)

	Sf5205 = MSbot2(2) - 7*MStop2(2)

	Sf5149 = MSbot2(2) - 8*MStop2(2)

	Sf2988 = MSbot2(2) + 4*MStop2(1)

	Sf2271 = MSbot2(2) + 3*MStop2(1)

	Sf2282 = MSbot2(2) + 2*MStop2(1)

	Sf77 = MSbot2(2) + MStop2(1)

	Sf7886 = MSbot2(2) - 5*MStop2(1)

	Sf3745 = MSbot2(2) - 7*MStop2(1)

	Sf2742 = MSbot2(2) - 8*MStop2(1)

	Sf5844 = MSbot2(1)**3 - MStop2(2)**3

	Sf4060 = MSbot2(1)**3 - MStop2(1)**3

	Sf5449 = 9*MSbot2(1)**2 - 3*MStop2(2)**2

	Sf6296 = 7*MSbot2(1)**2 + 9*MStop2(2)**2

	Sf3998 = 7*MSbot2(1)**2 + 9*MStop2(1)**2

	Sf3060 = 3*MSbot2(1)**2 - MStop2(1)**2

	Sf5022 = MSbot2(1)**2 + MStop2(2)**2

	Sf2555 = MSbot2(1)**2 + MStop2(1)**2

	Sf7517 = MSbot2(1)**2 + MSbot2(2)**2

	Sf7369 = Sf7368*MSbot2(1) - MStop2(2)

	Sf5017 = 4*MW2*MSbot2(1) - 6*MSq2Diff(1,2)**2

	Sf5185 = 2*MHiggs2(4)*MSbot2(1) - MSq2Diff(1,2)**2

	Sf8037 = 3113*MSbot2(1) + 3437*MStop2(2)

	Sf7828 = 3113*MSbot2(1) + 3437*MStop2(1)

	Sf8045 = 1789*MSbot2(1) + 1843*MStop2(2)

	Sf7836 = 1789*MSbot2(1) + 1843*MStop2(1)

	Sf8444 = 1546*MSbot2(1) + 1729*MStop2(2)

	Sf8420 = 1546*MSbot2(1) + 1729*MStop2(1)

	Sf8024 = 1009*MSbot2(1) + 955*MStop2(2)

	Sf7815 = 1009*MSbot2(1) + 955*MStop2(1)

	Sf8018 = 721*MSbot2(1) + 667*MStop2(2)

	Sf7798 = 721*MSbot2(1) + 667*MStop2(1)

	Sf8300 = 655*MSbot2(1) + 763*MStop2(2)

	Sf8276 = 655*MSbot2(1) + 763*MStop2(1)

	Sf8162 = 547*MSbot2(1) + 583*MStop2(2)

	Sf7953 = 547*MSbot2(1) + 583*MStop2(1)

	Sf8053 = 483*MSbot2(1) + 969*MStop2(2)

	Sf7844 = 483*MSbot2(1) + 969*MStop2(1)

	Sf8006 = 433*MSbot2(1) + 379*MStop2(2)

	Sf7803 = 433*MSbot2(1) + 379*MStop2(1)

	Sf8246 = 303*MSbot2(1) + 262*MStop2(2)

	Sf8234 = 303*MSbot2(1) + 262*MStop2(1)

	Sf8012 = 299*MSbot2(1) + 281*MStop2(2)

	Sf7809 = 299*MSbot2(1) + 281*MStop2(1)

	Sf8294 = 258*MSbot2(1) - 92*MStop2(2)

	Sf8270 = 258*MSbot2(1) - 92*MStop2(1)

	Sf8289 = 251*MSbot2(1) + 233*MStop2(2)

	Sf8265 = 251*MSbot2(1) + 233*MStop2(1)

	Sf8168 = 244*MSbot2(1) + 321*MStop2(2)

	Sf7959 = 244*MSbot2(1) + 321*MStop2(1)

	Sf8686 = 239*MSbot2(1) + 247*MStop2(2)

	Sf8646 = 239*MSbot2(1) + 247*MStop2(1)

	Sf5268 = 208*MSbot2(1) + 119*MStop2(2)

	Sf2816 = 208*MSbot2(1) + 119*MStop2(1)

	Sf8288 = 203*MSbot2(1) + 16*MStop2(2)

	Sf8264 = 203*MSbot2(1) + 16*MStop2(1)

	Sf5267 = 181*MSbot2(1) + 101*MStop2(2)

	Sf8029 = 181*MSbot2(1) + 46*MStop2(2)

	Sf2815 = 181*MSbot2(1) + 101*MStop2(1)

	Sf7820 = 181*MSbot2(1) + 46*MStop2(1)

	Sf8283 = 139*MSbot2(1) + 121*MStop2(2)

	Sf8259 = 139*MSbot2(1) + 121*MStop2(1)

	Sf5233 = 119*MSbot2(1) + 208*MStop2(2)

	Sf2764 = 119*MSbot2(1) + 208*MStop2(1)

	Sf5239 = 115*MSbot2(1) + 118*MStop2(2)

	Sf2770 = 115*MSbot2(1) + 118*MStop2(1)

	Sf8030 = 103*MSbot2(1) + 105*MStop2(2)

	Sf7821 = 103*MSbot2(1) + 105*MStop2(1)

	Sf5232 = 101*MSbot2(1) + 181*MStop2(2)

	Sf8069 = 101*MSbot2(1) + 107*MStop2(2)

	Sf2763 = 101*MSbot2(1) + 181*MStop2(1)

	Sf7860 = 101*MSbot2(1) + 107*MStop2(1)

	Sf5274 = 94*MSbot2(1) + 57*MStop2(2)

	Sf2822 = 94*MSbot2(1) + 57*MStop2(1)

	Sf8052 = 93*MSbot2(1) + 95*MStop2(2)

	Sf7843 = 93*MSbot2(1) + 95*MStop2(1)

	Sf8004 = 89*MSbot2(1) + 8*MStop2(2)

	Sf7796 = 89*MSbot2(1) + 8*MStop2(1)

	Sf8161 = 85*MSbot2(1) + 112*MStop2(2)

	Sf7952 = 85*MSbot2(1) + 112*MStop2(1)

	Sf7036 = 79*MSbot2(1) + 71*MStop2(2)

	end

