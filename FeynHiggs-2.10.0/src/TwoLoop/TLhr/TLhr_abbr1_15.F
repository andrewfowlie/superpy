* TLhr_abbr1_15.F
* this file is part of FeynHiggs
* generated 30 Nov 2011 16:05

#include "externals.h"
#include "types.h"


	subroutine TLhr_abbr1_15
	implicit none

#include "TLhr.h"


        Cd(13888) = 1/(32.D0*Pi**3)*
     &   (Alfas2L*(Cd(13466)*Cd(13885)*MSbot2(1) - 
     &        Cd(13891)*Cr(34)*MSq2Diff(1,1)*MStop2(2)))/
     &    (SW2*TB2*MSq2Diff(1,1)**2*MSq2Diff(1,2))

	Cd(13894) = -(Cd(50)*Cd(118)) + 2*Cd(47)*Cd(9969)*Cd(13841)

	Cd(13899) = Cd(47)*Cd(9948)*MSq2Diff(1,1) + Cd(50)*MSq2Diff(1,2)

        Cd(13900) = Cd(13899)*UStop2(1,1) + 
     &   Cd(47)*MSq2Diff(1,1)*UStop2(1,2)**2

        Cd(13911) = Cd(47)*Cd(13906)*MSq2Diff(1,1) + 
     &   Cd(50)*Cd(13910)*MSq2Diff(1,2)

        Cd(13920) = Cd(47)*Cd(13919)*MSq2Diff(1,1) - 
     &   Cd(50)*Cd(13913)*MSq2Diff(1,2)

	Cd(13930) = (Cd(13926) + Cd(13929))*MSq2Diff(1,2)

        Cd(13932) = -(Cd(13931)*UCStop(1,3)) + 
     &   MTy*MTy2*Cd(251)*Cd(12701)*UStop2(1,1)*UStop2(1,2)

	Cd(13933) = (AtC + MUE*TB)*Cd(253)

        Cd(13940) = MSq2Diff(1,2)*
     &   (Cd(13933) + AtC*MStop2(2)**2*UStop2(1,1))

	Cd(13954) = M_3C*Cd(13932) + M_3*Cd(13939)

        Cd(14046) = Cd(14036)*MSq2Diff(1,1)*MStop2(1)*
     &    UStopC(1,1) + MTy2*Cd(14045)*MStop2(3)*UUStopC(1,3)

        Cd(14049) = Cd(14046) - 
     &   Cd(14048)*MSq2Diff(1,1)**2*MStop2(1)*MStop2(3)*UStop2(1,2)

        Cd(14065) = Cd(9950)*Cd(14064) - 
     &   Cd(47)*Cd(14055)*MSq2Diff(1,1)*MStop2(3)

        Cd(14075) = 2*MGl2*Cd(14015)*MSq2Diff(1,1)*MStop2(1)*
     &    UStop2(1,1) + Cd(14074)*MStop2(3)*UUStopC(1,3)

        Cd(14087) = MUE*TB2*Cd(47)*Cd(9948)*UCStop(1,3) + 
     &   2*TB*Cd(239)*Cd(13469)*UStop2(1,1)*UStopC(1,1)

        Cd(14091) = Cd(14090)*UStop2(1,2) - 
     &   2*AtC*Cd(239)*Cd(9948)*Cd(13469)*UStopC(1,2)

	Cd(14092) = 2*Cd(14075) - Cd(14091)*MSq2Diff(1,1)*MStop2(3)

        Cd(14093) = -(Cd(9950)*Cd(14092)) + 
     &   Cd(47)*Cd(14047)*MSq2Diff(1,1)*MStop2(1)**2*MStop2(3)*
     &    UStop2(1,2)**2

        Cd(14102) = Cd(14101)*MSq2Diff(1,1)*UStopC(1,1) + 
     &   MTy*Cd(13469)*Cd(13815)*MSq2Diff(1,2)*UUStopC(1,3)

	Cd(14111) = Cd(14102) + Cd(14110)*MSq2Diff(1,1)

        Cd(14129) = MGl2*TB*Cd(164)*MSq2Diff(1,1) - 
     &   Cd(14128)*MGlpTmSt2(1)*UStop2(1,1)

	Cd(14157) = TB*Cd(164)*MSq2Diff(1,1)**2 + Cd(14156)*UStop2(1,1)

        Cd(14185) = Cd(251)*Cd(10103) - 
     &   Cd(50)*Cd(9969)*MGlpTmSt2(2)*UStop2(1,1)

        Cd(14197) = Cd(50)*MSbot2(1) - 
     &   Cd(14196)*MGlpTmSt2(2)*MSq2Diff(1,2)

	Cd(14198) = At*Cd(14197) + 2*M_3*Cd(50)*MSq2Diff(1,2)

	Cd(14199) = AtC*Cd(14197) + 2*M_3C*Cd(50)*MSq2Diff(1,2)

	Cd(14201) = MGlpTmSt2(2)*MSq2Diff(1,2) - Cd(50)*UStop2(1,1)

	Cd(14203) = TB*Cd(251)*MSq2Diff(1,2)**2 + Cd(14202)*UStop2(1,1)

	Cd(14205) = 2*M_3*Cd(50) + At*MGlpTmSt2(2)*MSq2Diff(1,2)

	Cd(14209) = Cd(50) - Cd(9969)*MGlpTmSt2(2)*MSq2Diff(1,2)

        Cd(14210) = Cd(9969)*MGlpTmSt2(2)*MSq2Diff(1,2) - 
     &   Cd(50)*UStop2(1,1)

	Cd(14211) = At*AtC*Cd(14209) - MUE2*TB2*Cd(14210)

	Cd(14213) = Cd(14211)*UStop2(1,1) + Cd(50)*Cd(14212)*UStop2(1,2)

	Cd(14214) = MGlpTmSt2(2)*MSq2Diff(1,2) - Cd(50)*UStop2(1,2)

	Cd(14215) = Cd(50)*Cd(9948)*UStop2(1,1) - Cd(14214)*UStop2(1,2)

	Cd(14216) = MTy2**2*Cd(14215) - Cd(14213)*MSbot2(1)

	Cd(14218) = Cd(50) - MGlpTmSt2(2)*MSq2Diff(1,2)*UStop2(1,2)

        Cd(14219) = Cd(9948)*MGlpTmSt2(2)*MSq2Diff(1,2)*
     &    UStop2(1,1) - Cd(14218)*UStop2(1,2)

        Cd(14227) = MUE2*TB2*Cd(50) + 
     &   2*M_3*AtC*MGlpTmSt2(2)*MSq2Diff(1,2)

	Cd(14234) = Cd(50) - MGlpTmSt2(2)*MSq2Diff(1,2)

        Cd(14238) = Cd(50)*Cd(14237) - 
     &   MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)*MStop2(3)*
     &    UStop2(1,1)

        Cd(14239) = 2*AtC*Cd(14238)*Cr(5) + 
     &   Cd(14234)*Cd(14235)*MSq2Diff(1,1)*MStop2(3)*UCStopC(1,3)

        Cd(14257) = Cd(50)*Cd(14256) + 
     &   Cd(14255)*MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)

        Cd(14261) = Cd(50)*Cd(14260) - 
     &   Cd(14259)*MGlpTmSt2(2)*MSq2Diff(1,1)*MSq2Diff(1,2)*
     &    MStop2(3)

	Cd(14293) = MUEC*TB*Cd(14287) + At*Cd(14292)

	Cd(14296) = Cd(12700)*Cd(14276) - Cd(14282)

        Cd(14327) = 2*Cd(14305)*MSq2Diff(1,2) - 
     &   AtC*Cd(14326)*UStopC(1,2)

	Cd(14333) = Cd(14327) + Cd(14332)*UStopC(1,1)

	Cd(14334) = Cd(14330)*MSq2Diff(1,1) + Cd(14333)*MSq2Diff(1,2)

	Cd(14343) = Cd(13466)*Cd(14325) + Cd(14334)

        Cd(14344) = 4*Cd(14343) - 
     &   Cd(13484)*Cd(14342)*MSq2Diff(1,1)*MSq2Diff(1,2)**2

	Cd(1446) = Cd(1445)*MStop2(3) - 2*Cd(1443)*UStop2(1,1)

	Cd(146) = 2*Cd(137) + Cd(145)*UStop2(1,1)

	Cd(1492) = (Cd(50)*Cd(1490) + Cd(47)*Cd(1491))/(Cd(47)*Cd(50))

        Cd(1515) = 2*Cd(50)*MGlpTmSt2(2)*MStop2(3) - 
     &   Cd(1514)*UStop2(1,2)

	Cd(1520) = Cd(50)*Cd(69) + Cd(1519)*UStop2(1,2)

	Cd(1539) = 4*Cd(50)*Cd(71) - Cd(1538)*UStop2(1,2)

	Cd(1551) = Cd(1539) - Cd(1550)*UStop2(1,1)

        Cd(1555) = Cd(1551)*MStop2(3) - 
     &   4*M_3*MTy*Cd(1554)*UCStopC(1,3)

	Cd(1566) = 2*Cd(50)*MGlpTmSt4(2) - Cd(1565)*UStop2(1,2)

	Cd(1603) = Cd(1597) - 3*MGl2*Cd(1602)

	Cd(1604) = 3*Cd(44)**3*Cd(89) - MGl2*Cd(1603)

        Cd(1605) = Cd(1601)*MStop2(3) - 
     &   4*M_3*MTy*Cd(1604)*UCStopC(1,3)

	Cd(159) = MUE*Cd(154)*UCStop(1,3) - MUEC*Cd(158)*UCStopC(1,3)

	Cd(1627) = Cd(44)**2*Cd(1620) - MGl2*Cd(1626)

	Cd(1650) = 3*Cd(44)*Cd(1645) - MGl2*Cd(1649)

	Cd(1651) = 2*Cd(44)**3*Cd(1644) - 2*MGl2*Cd(1650)

	Cd(1652) = 6*Cd(47)*Cd(58) - Cd(1651)*UStop2(1,2)**2

        Cd(1653) = Cd(1652)*MStop2(3) - 
     &   4*M_3C*MTy*Cd(1500)*UCStop(1,3)

        Cd(1654) = Cd(1653)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(47)*Cd(56)*UStop(1,2)*UStop2(1,2)

	Cd(1663) = 3*Cd(50)*MGlpTmSt4(2) - Cd(1662)*UStop2(1,2)**2

        Cd(1664) = Cd(1663)*MStop2(3) - 
     &   2*M_3C*MTy*Cd(1515)*UCStop(1,3)

        Cd(1665) = Cd(1664)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(50)*MGlpTmSt2(2)*MStop2(3)*UStop(1,2)*
     &    UStop2(1,2)

	Cd(1673) = 6*Cd(50)*Cd(71) - Cd(1672)*UStop2(1,2)**2

        Cd(1674) = Cd(1673)*MStop2(3) - 
     &   4*M_3C*MTy*Cd(1520)*UCStop(1,3)

        Cd(1675) = Cd(1674)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(50)*Cd(69)*UStop(1,2)*UStop2(1,2)

        Cd(1682) = 3*Cd(1679)*UStop(1,1) - 
     &   8*M_3*MTy*Cd(86)*Cd(1681)*UStop(1,2)

	Cd(1693) = Cd(44)*Cd(1689) - MGl2*Cd(1692)

	Cd(1694) = Cd(44)**3*Cd(1683) - MGl2*Cd(1693)

	Cd(1698) = Cd(1694)*UStop(1,1) - 8*M_3*MTy*Cd(1697)*UStop(1,2)

        Cd(1699) = Cd(52)*Cd(1682) + 
     &   2*Cd(86)**3*Cd(1493)*MStop2(3)*UStop(1,1)

        Cd(1700) = Cd(1699)*MStop2(3) + 
     &   Cd(1698)*T134(MGl2,MTy2,MStop2(1))

        Cd(1709) = Cd(1706)*UStop(1,1) + 
     &   8*M_3*MTy*Cd(96)*Cd(1708)*UStop(1,2)

	Cd(1721) = Cd(48)*Cd(1716) - MGl2*Cd(1720)

	Cd(1722) = Cd(48)**3*Cd(1710) - MGl2*Cd(1721)

	Cd(1726) = Cd(1722)*UStop(1,1) - 8*M_3*MTy*Cd(1725)*UStop(1,2)

        Cd(1727) = Cd(64)*Cd(1709) - 
     &   2*Cd(96)**3*Cd(1508)*MStop2(3)*UStop(1,1)

        Cd(1728) = Cd(1727)*MStop2(3) + 
     &   Cd(1726)*T134(MGl2,MTy2,MStop2(2))

        Cd(1729) = (Cd(50)**2*Cd(1700) - Cd(47)**2*Cd(1728))/
     &   (Cd(47)**2*Cd(50)**2)

        Cd(1771) = 2*M_3C*Cd(50)*Cd(69)*UCStop(1,3)*
     &    UStop(1,1) - M_3*Cd(1770)*UStop(1,2)

        Cd(1780) = Cd(1776)*MStop2(3)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(1779)*UStop(1,2)

        Cd(1781) = Cd(48)**3*Cd(1562) - 
     &   8*M_3C*MTy*Cd(50)*MGlpTmSt2(2)*UCStop(1,3)

        Cd(176) = (AtC*Cd(175)*UCStop(1,3) - 
     &     2*At*Cd(154)*MStop2(3)*UCStopC(1,3))/MStop2(3)

	Cd(1792) = -(MGl2*Cd(1602)) + Cd(1788)

	Cd(1793) = Cd(44)**3*Cd(89) - MGl2*Cd(1792)

	Cd(1794) = Cd(1791)*UStop(1,1) + 4*M_3*MTy*Cd(1793)*UStop(1,2)

	Cd(177) = (4*MTy*Cd(168) + Cd(176))/MTy

        Cd(178) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(-2*SA*Cd(159) + CA*MTy*Cd(177)))/
     &    (MTy*SB*SW*Cd(47)*MStop2(1))

	Cd(1838) = Cd(50)*Cd(1490) + Cd(47)*Cd(1491)

        Cd(1839) = 4*M_3C*MTy*Cd(1838)*UStop(1,1) - 
     &   Cd(85)*MStop2(3)*UStop(1,2)

	Cd(1854) = Cd(50)*Cd(69) + Cd(1745)*UStop2(1,2)

	Cd(1860) = 4*M_3*MTy*Cd(1854) - Cd(1859)*UCStop(1,3)

	Cd(1872) = Cd(47)*Cd(56) + Cd(1793)*UStop2(1,2)

	Cd(1878) = 4*M_3*MTy*Cd(1872) - Cd(1877)*UCStop(1,3)

	Cd(1880) = Cd(44)**2*Cd(1501) - MGl2*Cd(1879)*MStop2(3)

	Cd(1884) = Cd(44)**2*Cd(1620) - MGl2*Cd(1626)

        Cd(1885) = -2*M_3C*MTy*Cd(1884)*UCStop(1,3) + 
     &   Cd(1617)*UStop2(1,2)

	Cd(1886) = 2*Cd(47)*Cd(86)*Cd(87) - Cd(1885)*UStop(1,1)

	Cd(1892) = Cd(44)*Cd(1868) - 2*MGl2*Cd(1891)

	Cd(1893) = 2*Cd(44)**3*Cd(1847) - MGl2*Cd(1892)

	Cd(1903) = 2*Cd(1902)*UStop(1,1) - 2*M_3*MTy*Cd(93)*UStop(1,2)

        Cd(1904) = -2*M_3C*MTy*Cd(1896)*UCStop(1,3) + 
     &   Cd(1893)*UStop2(1,2)

	Cd(1905) = 2*Cd(47)*Cd(1903) - Cd(1904)*UStop(1,1)

        Cd(1906) = Cd(52)*Cd(1886) - 
     &   2*Cd(1493)*Cd(1633)*MStop2(3)*UStop(1,1)

        Cd(1907) = Cd(1906)*MStop2(3) + 
     &   Cd(1905)*T134(MGl2,MTy2,MStop2(1))

	Cd(1916) = 2*Cd(50)*Cd(96)*Cd(1914) + Cd(1915)*UStop(1,1)

	Cd(1921) = Cd(48)**2*Cd(1847) - MGl2*Cd(1920)

	Cd(1927) = 6*Cd(48)**3*Cd(89) + 2*MGl2*Cd(1926)

	Cd(1929) = MGl2*Cd(309) + Cd(1928)

	Cd(1930) = Cd(44)*Cd(48)**2 - MGl2*Cd(1929)

        Cd(1931) = 2*Cd(1930)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(102)*UStop(1,2)

        Cd(1932) = -(M_3C*MTy*Cd(1927)*UCStop(1,3)) + 
     &   Cd(1921)*MGlpTmSt2(2)*UStop2(1,2)

	Cd(1933) = Cd(50)*Cd(1931) - Cd(1932)*UStop(1,1)

        Cd(1934) = Cd(64)*Cd(1916) - 
     &   2*Cd(1508)*Cd(1510)*Cd(1908)*MStop2(3)*UStop(1,1)

        Cd(1935) = Cd(1934)*MStop2(3) - 
     &   2*Cd(1933)*T134(MGl2,MTy2,MStop2(2))

        Cd(1937) = (Cd(50)**2*Cd(1907) + Cd(47)**2*Cd(1935))/
     &   (Cd(47)**2*Cd(50)**2)

	Cd(1959) = Cd(1951) - MGl2*Cd(1958)

	Cd(1968) = Cd(47)*Cd(1961) - Cd(1967)*UStop2(1,1)

	Cd(1979) = Cd(44)**2*Cd(1973) - MGl2*Cd(1978)

	Cd(1991) = Cd(1990) - Cd(1979)*UStop2(1,1)

	Cd(1993) = 2*MTy2*Cd(1992) + Cd(1991)*UStop2(1,1)

        Cd(1994) = Cd(52)*Cd(1993) + 
     &   4*MTy2*Cd(1493)*Cd(1971)*MStop2(3)*UStop2(1,1)

	Cd(2014) = Cd(2013) - Cd(2005)*UStop2(1,1)

	Cd(2015) = Cd(50)*Cd(2006)*MStop2(3) + Cd(1996)*UStop2(1,1)**3

	Cd(2016) = 2*MTy2*Cd(2015) + Cd(2014)*MStop2(3)*UStop2(1,1)

	Cd(2034) = Cd(2027) - MGl2*Cd(2033)

	Cd(2044) = Cd(50)*Cd(2037) + Cd(2043)*UStop2(1,1)

	Cd(2055) = Cd(44)*Cd(2052) - 2*MGl2*Cd(2054)

	Cd(2056) = 2*Cd(44)**3*Cd(2048) - MGl2*Cd(2055)

	Cd(2077) = Cd(52)*Cd(2076) + 2*Cd(1493)*Cd(2062)*MStop2(3)

	Cd(2110) = Cd(64)*Cd(2109) - 2*Cd(1508)*Cd(2097)*MStop2(3)

        Cd(2111) = Cd(2110)*MStop2(3) - 
     &   Cd(2096)*T134(MGl2,MTy2,MStop2(2))

	Cd(2144) = Cd(44)**2*Cd(1562) - MGl2*Cd(2143)

	Cd(213) = 2*Cd(211)*MStop2(3) + Cd(212)*UStopC(1,1)

        Cd(214) = (2*AtC*Cd(185)*MStop2(3)*UCStop(1,3) + 
     &     At*Cd(193)*UCStopC(1,3))/MStop2(3)

	Cd(215) = (2*MTy*Cd(213) + Cd(214)*MStop2(3))/(MTy*MStop2(3))

        Cd(216) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(-2*SA*Cd(13)*Cd(180) + CA*MTy*Cd(215)))/
     &    (MTy*SB*SW*Cd(50)*MStop2(2))

	Cd(217) = 2*Cd(178)*MStop2(1) - 2*Cd(216)*MStop2(2)

	Cd(2218) = Cd(44)**3*Cd(2196) - MGl2*Cd(2217)

	Cd(221) = (2*MTy*Cd(168) - Cd(220))/MTy

	Cd(2238) = Cd(52)*Cd(2237) - 2*Cd(1493)*Cd(2191)*MStop2(3)

        Cd(222) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(-(SA*Cd(159)) + CA*MTy*Cd(221)))/
     &    (MTy*SB*SW*Cd(47)*MStop2(1))

        Cd(2245) = -(Cd(96)*Cd(254)*Cd(2244)) + 
     &   2*MTy*Cd(2243)*UStop2(1,2)

	Cd(2293) = Cd(64)*Cd(2292) - 2*Cd(1508)*Cd(2245)*MStop2(3)

        Cd(2294) = Cd(2293)*MStop2(3) - 
     &   Cd(2272)*T134(MGl2,MTy2,MStop2(2))

	Cd(2312) = Cd(1951) - MGl2*Cd(2311)

	Cd(2313) = Cd(44)**3*Cd(1948) - MGl2*Cd(2312)

        Cd(2318) = -(Cd(2128)*UCStop(1,3)) + 
     &   2*M_3*MTy*Cd(2317)*UStop2(1,2)

	Cd(2319) = M_3*Cd(2313) - MTy*Cd(2318)

	Cd(230) = AtC*Cd(185)*UCStop(1,3) - At*Cd(225)*UCStopC(1,3)

	Cd(231) = (2*MTy*Cd(229) + Cd(230))/MTy

        Cd(232) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(-(SA*Cd(13)*Cd(180)) + CA*MTy*Cd(231)))/
     &    (MTy*SB*SW*Cd(50)*MStop2(2))

	Cd(233) = 2*Cd(222) - Cd(178)*Cr(13)

	Cd(234) = 2*Cd(232) - Cd(216)*Cr(34)

	Cd(2376) = 2*M_3*Cd(2367) + M_3C*Cd(2375)*UCStop(1,3)**2

	Cd(235) = Cd(233)*MStop2(1) - Cd(234)*MStop2(2)

	Cd(236) = A0delStop(1)*Cd(178) - A0delStop(2)*Cd(216)

	Cd(2401) = M_3*Cd(1979) + MTy*Cd(2400)*UCStop(1,3)

	Cd(2407) = Cd(44)**2*Cd(1562) - MGl2*Cd(2143)

	Cd(2410) = Cd(2409) + M_3C*MTy*Cd(2407)*UUStop(1,1)

        Cd(2435) = 2*M_3*Cd(2431) + 
     &   M_3C*Cd(118)*Cd(2434)*MStop2(3)*UCStop(1,3)**2

	Cd(241) = MUE*Cd(238)*UCStop(1,3) + MUEC*Cd(240)*UCStopC(1,3)

	Cd(2461) = Cd(48)**2*Cd(1997) - MGl2*Cd(2460)

        Cd(2465) = -(Cd(2464)*UCStop(1,3)) + 
     &   2*M_3*MTy*Cd(96)*Cd(1708)*UStop2(1,2)

	Cd(2466) = M_3*Cd(2461) + MTy*Cd(2465)

	Cd(2483) = M_3*Cd(2467) - MTy*Cd(2482)

	Cd(2515) = Cd(2027) - MGl2*Cd(2514)

	Cd(2516) = Cd(48)**3*Cd(1948) - MGl2*Cd(2515)

	Cd(2531) = M_3*Cd(2516) - MTy*Cd(2530)

        Cd(2611) = -(M_3C*MTy*Cd(152)*Cd(2597)*UCStop(1,3)) + 
     &   Cd(2603)*UStop2(1,2)

        Cd(2622) = -(M_3C*MTy*Cd(2621)*UCStop(1,3)) + 
     &   Cd(2619)*UStop2(1,2)

	Cd(2623) = Cd(2606)*UStop(1,1) - M_3*MTy*Cd(2608)*UStop(1,2)

        Cd(2624) = Cd(2622)*UCStop(1,3) + 
     &   Cd(2623)*UStop2(1,1)*UStopC(1,2)

	Cd(2639) = 2*Cd(86)**2*MStop2(1)**2 + Cd(2619)*UStop2(1,2)

	Cd(2640) = 2*Cd(2606)*UStop(1,1) - M_3*MTy*Cd(2631)*UStop(1,2)

	Cd(2644) = 2*Cd(2606)*UStop(1,1) + M_3*MTy*Cd(2608)*UStop(1,2)

	Cd(2646) = 3*Cd(2619)*UStop(1,1) - M_3*MTy*Cd(2621)*UStop(1,2)

	Cd(2648) = 4*Cd(2647)*MStop2(1)**2 + Cd(2646)*UStop2(1,2)**2

	Cd(2660) = M_3*MTy2*Cd(2635) + MTy*Cd(2619)*UCStop(1,3)

	Cd(2667) = 2*Cd(2666)*MStop2(1)**2 - MTy*Cd(2665)*UUStopC(1,2)

	Cd(2668) = Cd(2660) - M_3*MTy2*Cd(152)*Cd(2597)*UStop2(1,1)

	Cd(2669) = Cd(2667) - Cd(2668)*UStop2(1,1)

	Cd(2673) = 2*Cd(2606)*UStop(1,1) - M_3*MTy*Cd(2608)*UStop(1,2)

        Cd(2694) = M_3C*MTy*Cd(179)*Cd(2681)*UCStop(1,3) + 
     &   Cd(2685)*UStop2(1,2)

        Cd(2700) = -(M_3C*MTy*Cd(179)*Cd(2681)*UCStop(1,3)) + 
     &   Cd(2699)*UStop2(1,2)

	Cd(2701) = Cd(2689)*UStop(1,1) + M_3*MTy*Cd(2691)*UStop(1,2)

        Cd(2702) = Cd(2700)*UCStop(1,3) + 
     &   Cd(2701)*UStop2(1,1)*UStopC(1,2)

	Cd(2715) = 2*Cd(96)**2*MStop2(2)**2 + Cd(2699)*UStop2(1,2)

	Cd(2716) = 2*Cd(2689)*UStop(1,1) + M_3*MTy*Cd(2707)*UStop(1,2)

	Cd(2720) = 2*Cd(2689)*UStop(1,1) - M_3*MTy*Cd(2691)*UStop(1,2)

	Cd(2724) = 3*Cd(2699)*UStop(1,1) + M_3*MTy*Cd(2723)*UStop(1,2)

	Cd(2726) = 4*Cd(2725)*MStop2(2)**2 + Cd(2724)*UStop2(1,2)**2

	Cd(2736) = M_3*MTy2*Cd(2711) - MTy*Cd(2689)*UCStop(1,3)

	Cd(2743) = 2*Cd(2742)*MStop2(2)**2 - MTy*Cd(2741)*UUStopC(1,2)

	Cd(2744) = Cd(2736) - M_3*MTy2*Cd(179)*Cd(2681)*UStop2(1,1)

	Cd(2745) = Cd(2743) - Cd(2744)*UStop2(1,1)

	Cd(2749) = 2*Cd(2689)*UStop(1,1) + M_3*MTy*Cd(2691)*UStop(1,2)

	Cd(2761) = Cd(2757) - Cd(2760)*UStop2(1,2)

	Cd(2762) = -4*Cd(1523)*Cd(2757) + 4*Cd(2761)*UStop2(1,1)

	Cd(2782) = M_3C*MTy*Cd(2762) + 2*Cd(2781)*UCStopC(1,3)

	Cd(277) = 2*Cd(276) - At*Cd(259)*UCStopC(1,3)

        Cd(2810) = 8*M_3C*Cd(2757)*MStop2(1)*MStop2(3)*
     &    UCStop(1,3)*UStop(1,1) + M_3*Cd(2809)*UStop(1,2)

        Cd(278) = (Cd(277) + 
     &     2*AtC*Cd(255)*MStop2(3)*UCStop(1,3))/MStop2(3)

	Cd(279) = -2*SA*Cd(13)*Cd(252) + CA*Cd(278)

	Cd(2845) = 4*(Cd(2757) - Cd(2760)*UStop2(1,2))

	Cd(2846) = Cd(2845) - 8*Cd(2757)*UStop2(1,1)

        Cd(2848) = Cd(2846)*MStop2(1)*UStop2(1,1) - 
     &   2*Cd(2847)*UStop2(1,2)

	Cd(2863) = Cd(44)*Cd(2828) - MGl2*Cd(2862)

	Cd(2883) = Cd(2882)*MStop2(1) + MTy2*Cd(2873)*UStop2(1,2)**2

	Cd(2885) = Cd(2883) + Cd(2884)*UStop2(1,1)

        Cd(2886) = 2*MTy2*Cd(47)*Cd(2874)*MStop2(1) - 
     &   Cd(2885)*UStop2(1,1)

	Cd(284) = 2*MTy*Cd(295)*Cr(5) + Cd(282)*MStop2(3)

	Cd(2899) = Cd(2892)*MStop2(1) + Cd(2898)*UStop2(1,2)

	Cd(2907) = Cd(44)*Cd(2901) - MGl2*Cd(2906)

        Cd(2908) = 4*Cd(2904)*MStop2(3)*UStop(1,1) - 
     &   M_3*MTy*Cd(2907)*UStop(1,2)

	Cd(2920) = 2*Cd(2916)*UStop(1,1) + M_3*MTy*Cd(2919)*UStop(1,2)

	Cd(2939) = 2*MTy*Cd(2935)*UStop(1,1) + M_3*Cd(2938)*UStop(1,2)

	Cd(2984) = 3*Cd(2980)*UStop(1,1) + M_3*MTy*Cd(2983)*UStop(1,2)

	Cd(3003) = 2*Cd(3002)*MStop2(1) - MTy2*Cd(2984)*UStop2(1,2)**2

        Cd(3031) = M_3C*MTy2*Cd(3012)*MStop2(3)*UCStop(1,3)*
     &    UUStop(1,1) + M_3*Cd(3030)*UUStop(1,3)

	Cd(300) = Cd(47)*MGlpTmSt2(2) + Cd(50)*MGlpTmSt2(1)*UStop2(1,2)

	Cd(3055) = Cd(3054)*MStop2(1) + MTy2*Cd(3037)*UStop2(1,2)

	Cd(3056) = Cd(3034)*MStop2(1)**2 + Cd(3055)*UStop2(1,2)**2

	Cd(3064) = Cd(3061) + 4*MTy2*Cd(3063)*UStop2(1,2)

	Cd(3065) = -(Cd(3064)*MStop2(1)) + 2*MTy2*Cd(2659)*MStop2(3)

	Cd(3069) = M_3*Cd(2869) - 2*MTy*Cd(3068)*UCStop(1,3)

	Cd(3074) = Cd(3073)*UStop(1,1) - M_3*MTy*Cd(2757)*UStop(1,2)

        Cd(3075) = Cd(2898)*UCStop(1,3)*UStop2(1,2) + 
     &   Cd(3074)*MStop2(1)*UStopC(1,2)

        Cd(3076) = 2*M_3*MTy*Cd(47)*MStop2(1)**2 + 
     &   Cd(3075)*UStop2(1,2)

	Cd(3087) = 2*(Cd(3086)*MStop2(1) + Cd(2665)*MStop2(3))

	Cd(3088) = M_3*Cd(2882) + 2*MTy*Cd(3079)*UCStop(1,3)

	Cd(3089) = Cd(3088)*MStop2(1) + MTy*Cd(3087)*UUStopC(1,2)

	Cd(3090) = Cd(3069) - M_3*MTy2*Cd(2873)*UStop2(1,1)

	Cd(3091) = Cd(3089) + Cd(3090)*UStop2(1,1)

	Cd(3096) = Cd(2883) + Cd(2869)*UStop2(1,1)

	Cd(3098) = MTy2*Cd(3097) - Cd(3096)*UStop2(1,1)

	Cd(3105) = 4*Cd(3104)*MStop2(1) + 2*Cd(2673)*MStop2(3)

        Cd(3156) = Cd(3155)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(3116)*MStop2(2)*UStop(1,2)

	Cd(3176) = -(MGl2*Cd(3175)) + Cd(48)**2*Cd(319)*MStop2(2)

	Cd(3182) = Cd(3172) + 2*Cd(3156)*MStop2(3)*UStop2(1,1)

        Cd(314) = M_3C*Cd(305)*UCStop(1,3) + 
     &   M_3*MStop2(1)*MStop2(2)*UCStopC(1,3)

	Cd(3192) = 4*Cd(3191)*MStop2(2) - 2*Cd(179)*Cd(2681)*MStop2(3)

        Cd(3197) = M_3C*MTy*Cd(3192)*UCStop(1,3) + 
     &   2*Cd(3188)*MStop2(3)*UStop2(1,2)

	Cd(3246) = Cd(3228) + 4*MTy2*Cd(3119)*UStop2(1,1)

	Cd(3247) = Cd(3245) - Cd(3246)*UStop2(1,1)

	Cd(3248) = Cd(3247)*MStop2(2) + 2*MTy2*Cd(2714)*MStop2(3)

        Cd(3249) = 2*MTy2*Cd(50)*Cd(3229)*MStop2(2) - 
     &   Cd(3248)*UStop2(1,1)

	Cd(3261) = Cd(3255)*MStop2(2) + Cd(3260)*UStop2(1,2)

	Cd(3284) = 2*Cd(3280)*UStop(1,1) - M_3*MTy*Cd(3283)*UStop(1,2)

	Cd(3304) = 2*MTy*Cd(3300)*UStop(1,1) - M_3*Cd(3303)*UStop(1,2)

        Cd(329) = Cd(328)*UCStopC(1,3) + 
     &   2*M_3C*MTy*Cd(326)*UStop2(1,1)

	Cd(3344) = Cd(48)*Cd(3336) + 2*MGl2*Cd(3343)

	Cd(3345) = 3*Cd(3340)*UStop(1,1) - M_3*MTy*Cd(3344)*UStop(1,2)

        Cd(3364) = Cd(3361)*UStop(1,1) - 
     &   M_3*MTy*Cd(3363)*MStop2(2)*MStop2(3)*UStop(1,2)

	Cd(3365) = 2*Cd(3364)*MStop2(2) - MTy2*Cd(3345)*UStop2(1,2)**2

	Cd(334) = 2*M_3*MTy*Cd(50) - Cd(333)*UStopC(1,2)

	Cd(3416) = Cd(3415)*MStop2(2) + MTy2*Cd(3397)*UStop2(1,2)

	Cd(3417) = Cd(3395)*MStop2(2)**2 + Cd(3416)*UStop2(1,2)**2

	Cd(3418) = Cd(3417)*MStop2(3) + M_3C*MTy*Cd(3408)*UCStop(1,3)

	Cd(3423) = 4*Cd(3422)*MStop2(2) - 2*Cd(179)*Cd(2681)*MStop2(3)

	Cd(3433) = Cd(3431)*MStop2(2) - 2*MTy2*Cd(3432)*MStop2(3)

	Cd(3448) = Cd(3447)*UStop(1,1) + M_3*MTy*Cd(3116)*UStop(1,2)

        Cd(3449) = Cd(3260)*UCStop(1,3)*UStop2(1,2) + 
     &   Cd(3448)*MStop2(2)*UStopC(1,2)

        Cd(3450) = 2*M_3*MTy*Cd(50)*MStop2(2)**2 - 
     &   Cd(3449)*UStop2(1,2)

	Cd(3461) = 2*(Cd(3460)*MStop2(2) + Cd(2741)*MStop2(3))

	Cd(3462) = M_3*Cd(3452) - 2*MTy*Cd(3455)*UCStop(1,3)

	Cd(3463) = Cd(3462)*MStop2(2) - MTy*Cd(3461)*UUStopC(1,2)

	Cd(3464) = Cd(3444) - M_3*MTy2*Cd(3423)*UStop2(1,1)

	Cd(3465) = Cd(3463) + Cd(3464)*UStop2(1,1)

	Cd(3475) = 2*Cd(3280)*UStop(1,1) - M_3*MTy*Cd(3474)*UStop(1,2)

	Cd(3482) = 4*Cd(3481)*MStop2(2) + 2*Cd(2749)*MStop2(3)

        Cd(3500) = Cd(3499)*UStop(1,1) + 
     &   8*M_3*MTy*Cd(2757)*MStop2(3)*UStop(1,2)

        Cd(3503) = Cd(2317)*MStop2(1)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(3502)*MStop2(3)*UStop(1,2)

	Cd(3533) = Cd(3528) + MGl2*Cd(3532)

	Cd(3534) = MGl2*Cd(3533) - Cd(44)**2*Cd(1847)*MStop2(1)

	Cd(3543) = Cd(3528) + MGl2*Cd(3542)

	Cd(3544) = MGl2*Cd(3543) - Cd(44)**2*Cd(1847)*MStop2(1)

	Cd(3545) = 2*M_3*MTy*Cd(3540)*MStop2(3) - Cd(3544)*UCStop(1,3)

	Cd(3555) = Cd(3554) + 2*MTy2*Cd(3553)*UStop2(1,2)**2

	Cd(3557) = Cd(3555) + Cd(3556)*UStop2(1,1)

	Cd(3558) = 2*MTy2*Cd(47)*Cd(2874) - Cd(3557)*UStop2(1,1)

	Cd(3567) = Cd(3561) + 2*MGl2*Cd(3566)

	Cd(3568) = -(MGl2*Cd(3567)) + 2*Cd(44)**2*Cd(2048)*MStop2(1)

	Cd(3576) = Cd(3571) - MGl2*Cd(3575)

	Cd(3588) = Cd(3585) - MGl2*Cd(3587)

	Cd(3589) = 2*MGl2*Cd(3588) - Cd(44)**2*Cd(2115)*MStop2(1)

        Cd(3590) = Cd(3589)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(3103)*MStop2(3)*UStop(1,2)

        Cd(3605) = 2*(MTy*Cd(3602)*UStop(1,1) + 
     &     M_3*Cd(3604)*MStop2(3)*UStop(1,2))

	Cd(3611) = Cd(3585) - MGl2*Cd(3610)

	Cd(3612) = 2*MGl2*Cd(3611) - Cd(44)**2*Cd(2115)*MStop2(1)

	Cd(3613) = M_3*Cd(3552)*MStop2(3) - MTy*Cd(3612)*UCStop(1,3)

	Cd(3617) = Cd(3615) - MGl2*Cd(3616)

	Cd(3627) = M_3*Cd(3554) + 2*MTy*Cd(3626)*UUStopC(1,2)

	Cd(3628) = Cd(3627)*MStop2(3) + MTy*Cd(3625)*UCStop(1,3)

        Cd(3629) = Cd(3613) - 
     &   2*M_3*MTy2*Cd(3553)*MStop2(3)*UStop2(1,1)

	Cd(3630) = Cd(3628) + Cd(3629)*UStop2(1,1)

        Cd(3631) = M_3C*Cd(3064)*MStop2(3)*UCStop(1,3)**2 + 
     &   Cd(3630)*UStop2(1,1)

	Cd(3636) = 4*Cd(1523)*Cd(2757) - Cd(2845)*UStop2(1,1)

	Cd(3648) = M_3C*MTy*Cd(3636) - 2*Cd(3647)*UCStopC(1,3)

	Cd(3649) = Cd(3555) + Cd(3552)*UStop2(1,1)

	Cd(3651) = 2*MTy2*Cd(3650) - Cd(3649)*UStop2(1,1)

	Cd(3656) = Cd(3585) - MGl2*Cd(3655)

        Cd(3673) = Cd(3670)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(3672)*MStop2(3)*UStop(1,2)

        Cd(3674) = 2*MStop2(1)*
     &   (3*Cd(2317)*MStop2(1)*UStop(1,1) - 
     &     4*M_3*MTy*MTy2*Cd(167)*MStop2(3)**2*UStop(1,2))

        Cd(3709) = Cd(118)**2*Cd(3708)*MStop2(1) + 
     &   Cd(3674)*UStop(1,1)*UStop2(1,2)

	Cd(3716) = Cd(3715) + 4*MTy2*Cd(2760)*UStop2(1,2)

	Cd(3717) = 2*Cd(3714)*MStop2(1)*MStop2(3) + Cd(3716)*UStop2(1,2)

	Cd(3726) = Cd(3723) - MGl2*Cd(3725)

	Cd(3727) = MGl2*Cd(3726) - Cd(44)**2*Cd(3719)*MStop2(1)

	Cd(3728) = Cd(3711)*MStop2(1) - 2*MTy2*Cd(3713)*UStop2(1,2)**3

        Cd(3743) = 2*Cd(3740)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(3742)*MStop2(1)*MStop2(3)*UStop(1,2)

	Cd(3746) = Cd(3743) - 2*MTy2*Cd(3735)*UStop2(1,2)**2

        Cd(3768) = Cd(3765)*UStop(1,1) - 
     &   8*M_3*MTy*Cd(3767)*MStop2(3)*UStop(1,2)

	Cd(3769) = -2*Cd(3768)*MStop2(2) + MGl2**2*Cd(3759)*UStop(1,1)

        Cd(3780) = Cd(3777)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(3779)*MStop2(3)*UStop(1,2)

	Cd(3813) = 2*M_3*MTy*Cd(3809)*MStop2(3) - Cd(3812)*UCStop(1,3)

	Cd(3821) = 2*MTy2*Cd(50)*Cd(3229) - Cd(3247)*UStop2(1,1)

	Cd(3839) = Cd(3829) - 2*Cd(3838)*UStop2(1,2)

	Cd(3849) = MGl2**3*Cd(3840) + 2*Cd(3848)*MStop2(2)

        Cd(3850) = Cd(3849)*UStop(1,1) - 
     &   8*M_3*MTy*Cd(3480)*MStop2(2)*MStop2(3)*UStop(1,2)

        Cd(3867) = 2*(MTy*Cd(3864)*UStop(1,1) + 
     &     M_3*Cd(3866)*MStop2(3)*UStop(1,2))

	Cd(3868) = Cd(3867) + 8*MTy*Cd(3481)*MStop2(3)*UStop2(1,2)

	Cd(3869) = -(Cd(3868)*MStop2(2)) + MTy*Cd(3850)*UStop2(1,1)

        Cd(3870) = 2*MTy*Cd(3839)*MStop2(2)*UStop2(1,2) + 
     &   Cd(3869)*UStopC(1,1)

        Cd(3871) = 2*M_3C*Cd(3821)*MStop2(2)*MStop2(3) + 
     &   Cd(3870)*UCStopC(1,3)

	Cd(3872) = Cd(3431) + 4*MTy2**2*Cd(251)*MStop2(3)*UStop2(1,2)

        Cd(3874) = Cd(3849)*UCStop(1,3) + 
     &   4*M_3*MTy*Cd(3116)*MStop2(2)*MStop2(3)*UStop2(1,2)

	Cd(3875) = MTy*Cd(3874) - 2*M_3*Cd(3873)*MStop2(2)*MStop2(3)

	Cd(3880) = Cd(3878) - 2*MGl2*Cd(3879)

	Cd(380) = SA*Cd(39) + 2*CA*Cd(13)*Cr(75)

	Cd(3892) = M_3*Cd(3241) - 2*MTy*Cd(3891)*UUStopC(1,2)

        Cd(381) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(2*CA*Cd(159) + MTy*SA*Cd(177)))/
     &    (MTy*SB*SW*Cd(47)*MStop2(1))

        Cd(382) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(2*CA*Cd(13)*Cd(180) + MTy*SA*Cd(215)))/
     &    (MTy*SB*SW*Cd(50)*MStop2(2))

	Cd(3907) = M_3C*MTy*Cd(3121) + 2*Cd(3906)*UCStopC(1,3)

	Cd(3911) = Cd(3867) + 8*MTy*Cd(3910)*MStop2(3)*UStop2(1,2)

	Cd(3912) = -(Cd(3911)*MStop2(2)) + MTy*Cd(3850)*UStop2(1,1)

        Cd(3913) = 2*MTy*Cd(3839)*MStop2(2)*UStop2(1,2) + 
     &   Cd(3912)*UStopC(1,1)

	Cd(383) = 2*Cd(381)*MStop2(1) - 2*Cd(382)*MStop2(2)

        Cd(3914) = 2*M_3C*Cd(3821)*MStop2(2)*MStop2(3) + 
     &   Cd(3913)*UCStopC(1,3)

        Cd(384) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(CA*Cd(159) + MTy*SA*Cd(221)))/
     &    (MTy*SB*SW*Cd(47)*MStop2(1))

        Cd(3928) = Cd(3925)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(3927)*MStop2(3)*UStop(1,2)

        Cd(3931) = 3*Cd(3930)*UStop(1,1) - 
     &   4*M_3*MTy*MTy2*Cd(254)*MStop2(3)**2*UStop(1,2)

        Cd(385) = 1/(16.D0*Pi**3)*
     &   (Alfas2L*(CA*Cd(13)*Cd(180) + MTy*SA*Cd(231)))/
     &    (MTy*SB*SW*Cd(50)*MStop2(2))

        Cd(3936) = 4*M_3C*MTy2*Cd(3933)*MStop2(3)*UCStop(1,3)*
     &    UUStop(1,1) - M_3*Cd(3935)*UUStop(1,3)

	Cd(386) = 2*Cd(384) - Cd(381)*Cr(13)

	Cd(387) = 2*Cd(385) - Cd(382)*Cr(34)

	Cd(388) = Cd(386)*MStop2(1) - Cd(387)*MStop2(2)

        Cd(3964) = -(Cd(138)**2*Cd(3963)) + 
     &   2*Cd(3931)*UStop(1,1)*UStop2(1,2)

	Cd(389) = A0delStop(1)*Cd(381) - A0delStop(2)*Cd(382)

	Cd(3978) = Cd(3974) - MGl2*Cd(3977)

	Cd(3979) = MGl2*Cd(3978) - Cd(48)**2*Cd(3970)*MStop2(2)

	Cd(3980) = Cd(3966)*MStop2(2) - 2*MTy2*Cd(3968)*UStop2(1,2)**3

        Cd(3998) = 2*Cd(3995)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(3997)*MStop2(2)*MStop2(3)*UStop(1,2)

        Cd(4002) = M_3*MTy*Cd(4001)*UUStop(1,2) + 
     &   2*Cd(50)**2*MStop2(2)*UUStop(2,3)

	Cd(392) = 2*CA*Cd(13)*Cd(252) + SA*Cd(278)

	Cd(4003) = Cd(3998) - 2*MTy2*Cd(3988)*UStop2(1,2)**2

	Cd(4019) = Cd(44)**2*Cd(2763) - MGl2*Cd(4021)

	Cd(4020) = Cd(44)**3*Cd(1735) - MGl2*Cd(4019)

        Cd(4036) = Cd(4035)*UStop(1,1) + 
     &   M_3*MTy*Cd(4025)*UStop(1,2)*UStop2(1,2)

	Cd(4062) = MTy*Cd(4055)*UStop2(1,2) + Cd(4061)*UStopC(1,1)

        Cd(4063) = -(Cd(4062)*UCStopC(1,3)) + 
     &   M_3C*Cd(4053)*UStop2(1,1)

	Cd(4072) = Cd(4071)*MStop2(1) + MTy*Cd(4070)*UUStopC(1,2)

	Cd(4076) = Cd(4075)*MStop2(1) + Cd(4035)*UStop2(1,2)

        Cd(4094) = MTy2*Cd(4091)*UUStop(1,1) - 
     &   2*M_3*MTy*Cd(4093)*UUStop(1,3)

	Cd(4109) = MTy*Cd(4055)*UStop2(1,2) + Cd(4108)*UStopC(1,1)

        Cd(4110) = -(Cd(4109)*UCStopC(1,3)) + 
     &   M_3C*Cd(4053)*UStop2(1,1)

	Cd(4119) = Cd(48)**2*Cd(3122) - MGl2*Cd(4121)

	Cd(4120) = Cd(48)**3*Cd(4117) - MGl2*Cd(4119)

        Cd(4128) = M_3C*MTy*Cd(4116)*UCStop(1,3) + 
     &   Cd(4120)*UStop2(1,2)

        Cd(4131) = Cd(4130)*UStop(1,1) - 
     &   M_3*MTy*Cd(4125)*UStop(1,2)*UStop2(1,2)

        Cd(4143) = Cd(4142)*UCStopC(1,3) + 
     &   M_3C*MTy*Cd(4116)*UStop2(1,1)*UStop2(1,2)

	Cd(4160) = MTy*Cd(4153)*UStop2(1,2) + Cd(4159)*UStopC(1,1)

	Cd(4161) = Cd(4160)*UCStopC(1,3) + M_3C*Cd(4151)*UStop2(1,1)

	Cd(4170) = Cd(4169)*MStop2(2) - MTy*Cd(4168)*UUStopC(1,2)

	Cd(4174) = Cd(4173)*MStop2(2) + Cd(4130)*UStop2(1,2)

        Cd(4189) = MTy2*Cd(4185)*UUStop(1,1) + 
     &   2*M_3*MTy*Cd(4188)*UUStop(1,3)

        Cd(4198) = MTy2*Cd(4176)*UStop2(1,1)**2*UStopC(1,1) + 
     &   Cd(4189)*UUStopC(1,1)

	Cd(4204) = MTy*Cd(4153)*UStop2(1,2) + Cd(4203)*UStopC(1,1)

	Cd(4205) = Cd(4204)*UCStopC(1,3) + M_3C*Cd(4151)*UStop2(1,1)

	Cd(4236) = 3*Cd(44)**3 - MGl2*Cd(4235)

	Cd(416) = Cd(413) - Cd(415)

        Cd(4265) = 8*M_3C*Cd(4217)*MStop2(1)*MStop2(3)*
     &    UCStop(1,3)*UStop(1,1) + M_3*Cd(4264)*UStop(1,2)

	Cd(4274) = Cd(4273)*MStop2(1) - Cd(4031)*MStop2(3)*UStop2(1,2)

        Cd(4275) = 2*Cd(4274)*UStop2(1,1) - 
     &   4*M_3*MTy*Cd(4217)*MStop2(1)*UCStopC(1,3)*UStop2(1,2)

	Cd(422) = Cd(419) - Cd(421)

	Cd(4317) = Cd(44)*Cd(4306) - MGl2*Cd(4316)

	Cd(4318) = Cd(44)**3*Cd(4307) - MGl2*Cd(4317)

	Cd(4351) = Cd(4344)*MStop2(1) - Cd(4350)*UStop2(1,2)

	Cd(4391) = 2*MTy*Cd(4351)*UStop2(1,2) - Cd(4390)*UStopC(1,1)

	Cd(4392) = M_3C*Cd(4339) + Cd(4391)*UCStopC(1,3)

	Cd(4405) = 3*Cd(44)*Cd(4365) - MGl2*Cd(4404)

	Cd(4406) = Cd(44)**3*Cd(4363) - MGl2*Cd(4405)

	Cd(4407) = M_3*Cd(4403) + 2*MTy*Cd(4406)*UCStop(1,3)

	Cd(4417) = 2*(Cd(4416)*MStop2(1) + Cd(4070)*MStop2(3))

	Cd(4418) = M_3*Cd(4408) + 2*MTy*Cd(4411)*UCStop(1,3)

	Cd(4419) = Cd(4418)*MStop2(1) - MTy*Cd(4417)*UUStopC(1,2)

	Cd(4452) = Cd(4450)*MStop2(1) + 4*Cd(4451)*MStop2(3)

	Cd(4470) = 3*Cd(4465)*UStop(1,1) - M_3*MTy*Cd(4469)*UStop(1,2)

	Cd(4480) = Cd(44)**3*Cd(4473) - MGl2*Cd(4479)

        Cd(4484) = 2*MTy*Cd(4480)*UStop(1,1) - 
     &   M_3*Cd(4483)*MStop2(1)*MStop2(3)*UStop(1,2)

	Cd(4485) = Cd(4484)*MStop2(1) - MTy*MTy2*Cd(4470)*UStop2(1,2)**2

	Cd(441) = Cd(413)*Cr(34) + 2*Cd(440)*MStop2(3)*UCStopC(1,3)

	Cd(443) = Cd(415)*Cr(13) + 2*Cd(442)*MStop2(3)*UCStopC(1,3)

        Cd(4518) = M_3C*MTy2*Cd(4494)*MStop2(3)*UCStop(1,3)*
     &    UUStop(1,1) + M_3*Cd(4517)*UUStop(1,3)

	Cd(444) = Cd(441) - Cd(443)

	Cd(4523) = Cd(4457)*MStop2(1) - MTy2*Cd(4522)*UStop2(1,2)

        Cd(4524) = Cd(4521)*MStop2(1)**2*MStop2(3) + 
     &   MTy2*Cd(4523)*UStop2(1,2)

        Cd(445) = 1/2.D0*(A0delStop(2)*Cd(413)) + 
     &   Cd(440)*Cr(34)*MStop2(2)*MStop2(3)*UCStopC(1,3)

	Cd(4537) = MTy*Cd(4536)*MStop2(3) - M_3C*Cd(4524)*UCStop(1,3)

        Cd(446) = 1/2.D0*(A0delStop(1)*Cd(415)) + 
     &   Cd(442)*Cr(13)*MStop2(1)*MStop2(3)*UCStopC(1,3)

	Cd(4545) = Cd(4229)*MStop2(1) + Cd(4020)*MStop2(3)

        Cd(447) = 3/(16.D0*Pi)*
     &   (Alfa2L*CA*(-(Cd(445)*MStop2(1)) + Cd(446)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(2)*MStop2(3))

	Cd(4550) = 2*Cd(44)**3 - MGl2*Cd(4549)

        Cd(4551) = Cd(4415)*MStop2(3)*UStop(1,1) + 
     &   M_3*MTy*Cd(4550)*UStop(1,2)

	Cd(4552) = 4*Cd(4551)*MStop2(1) - 2*Cd(4107)*MStop2(3)

	Cd(4553) = Cd(4374)*UStop2(1,1) + Cd(4552)*UStop2(1,2)

	Cd(4554) = -(MTy*Cd(4553)) + Cd(4389)*MStop2(1)

	Cd(4555) = 2*MTy*Cd(4351)*UStop2(1,2) - Cd(4554)*UStopC(1,1)

	Cd(4556) = M_3C*Cd(4339) + Cd(4555)*UCStopC(1,3)

	Cd(449) = Cd(419)*Cr(34) + 2*Cd(448)*MStop2(3)*UCStop(1,3)

	Cd(4586) = M_3C*MTy*Cd(4567) - 2*Cd(4585)*UCStopC(1,3)

	Cd(451) = Cd(421)*Cr(13) + 2*Cd(450)*MStop2(3)*UCStop(1,3)

	Cd(452) = Cd(449) - Cd(451)

        Cd(4604) = 2*Cd(4569)*MStop2(2)*UStop(1,1) + 
     &   M_3*MTy*Cd(4603)*UStop(1,2)

        Cd(453) = 1/2.D0*(A0delStop(2)*Cd(419)) + 
     &   Cd(448)*Cr(34)*MStop2(2)*MStop2(3)*UCStop(1,3)

        Cd(454) = 1/2.D0*(A0delStop(1)*Cd(421)) + 
     &   Cd(450)*Cr(13)*MStop2(1)*MStop2(3)*UCStop(1,3)

	Cd(4620) = Cd(4619)*MStop2(2) + Cd(4128)*MStop2(3)*UStop2(1,2)

        Cd(4621) = 2*Cd(4620)*UStop2(1,1) + 
     &   4*M_3*MTy*Cd(4562)*MStop2(2)*UCStopC(1,3)*UStop2(1,2)

        Cd(455) = 3/(16.D0*Pi)*
     &   (Alfa2L*CA*(-(Cd(453)*MStop2(1)) + Cd(454)*MStop2(2)))/
     &    (SB2*SW2*MStop2(1)*MStop2(2)*MStop2(3))

	Cd(4652) = 6*Cd(48)*Cd(4642) - MGl2*Cd(4651)

	Cd(4653) = Cd(48)**3*Cd(4643) - MGl2*Cd(4652)

	Cd(4676) = Cd(4665) - 4*MTy2*Cd(4565)*UStop2(1,1)

	Cd(4677) = Cd(4675) - Cd(4676)*UStop2(1,1)

	Cd(4678) = Cd(4677)*MStop2(2) + 2*Cd(4151)*MStop2(3)

	Cd(4679) = 2*Cd(50)*Cd(4666)*MStop2(2) - Cd(4678)*UStop2(1,1)

	Cd(4688) = Cd(4684)*MStop2(2) - Cd(4687)*UStop2(1,2)

	Cd(4727) = 2*MTy*Cd(4688)*UStop2(1,2) - Cd(4726)*UStopC(1,1)

	Cd(4728) = M_3C*Cd(4679) - Cd(4727)*UCStopC(1,3)

	Cd(4754) = 2*(Cd(4753)*MStop2(2) + Cd(4168)*MStop2(3))

	Cd(4755) = M_3*Cd(4742) - 2*MTy*Cd(4746)*UCStop(1,3)

	Cd(4756) = Cd(4755)*MStop2(2) - MTy*Cd(4754)*UUStopC(1,2)

	Cd(4786) = Cd(4785) + 8*MGl2*MTy2*Cd(179)*MStop2(3)*UStop2(1,2)

	Cd(4787) = Cd(4786)*MStop2(2) - 4*Cd(4162)*MStop2(3)

	Cd(4813) = Cd(48)**3*Cd(4806) - MGl2*Cd(4812)

        Cd(4817) = 2*MTy*Cd(4813)*UStop(1,1) - 
     &   M_3*Cd(4816)*MStop2(2)*MStop2(3)*UStop(1,2)

	Cd(4818) = Cd(4817)*MStop2(2) + MTy*MTy2*Cd(4801)*UStop2(1,2)**2

        Cd(4825) = 6*Cd(4569)*MStop2(2)*UStop(1,1) - 
     &   M_3*MTy*Cd(4824)*MStop2(3)*UStop(1,2)

        Cd(4840) = 2*Cd(4825)*UStop(1,1)*UStop2(1,2) - 
     &   Cd(138)*Cd(4839)*MStop2(2)*UUStop(1,1)

	Cd(4859) = Cd(4858) - 4*MTy2*Cd(4565)*UStop2(1,2)

        Cd(4860) = Cd(4854)*MStop2(2)*MStop2(3) - 
     &   MTy2*Cd(4859)*UStop2(1,2)

	Cd(4861) = Cd(4860)*MStop2(2) + 2*Cd(4192)*MStop2(3)

	Cd(4875) = MTy*Cd(4874)*MStop2(3) + M_3C*Cd(4861)*UCStop(1,3)

	Cd(4883) = 2*Cd(4569) + MGl2**3*Cd(3122)*UStop2(1,1)

	Cd(4897) = 2*MTy*Cd(4688)*UStop2(1,2) - Cd(4896)*UStopC(1,1)

	Cd(4898) = M_3C*Cd(4679) - Cd(4897)*UCStopC(1,3)

	Cd(4939) = 2*Cd(4929) - MGl2*Cd(4938)

	Cd(4940) = Cd(4935) + MGl2*Cd(4939)

	Cd(4941) = 2*(Cd(44)*Cd(48)*Cd(4925) + MGl2*Cd(4940))

        Cd(4942) = -4*Cd(47)*Cd(50)*Cd(295)*Cr(28) + 
     &   Cd(4941)*UStop2(1,1)*UStop2(1,2)

	Cd(4952) = Cd(44)**2*Cd(48)*Cd(4943) + MGl2*Cd(4951)

	Cd(5009) = Cd(5002) - MGl2*Cd(5008)

	Cd(5010) = Cd(44)*Cd(48)*Cd(4995) - MGl2*Cd(5009)

	Cd(5025) = Cd(5022) - MGl2*Cd(5024)

	Cd(5026) = Cd(44)*Cd(5018) - MGl2*Cd(5025)

	Cd(5027) = Cd(47)**2*Cd(5015) + Cd(50)*Cd(5026)*UStop2(1,2)

	Cd(5069) = Cd(5035) - MGl2*Cd(5068)

	Cd(5070) = Cd(5039) + MGl2*Cd(5069)

	Cd(5071) = Cd(5045) - MGl2*Cd(5070)

	Cd(5072) = Cd(5051) - MGl2*Cd(5071)

	Cd(5073) = Cd(44)**2*Cd(48)**2*Cd(5032) - MGl2*Cd(5072)

	Cd(5079) = Cd(5073)*UStop(1,1) + 2*M_3*MTy*Cd(5078)*UStop(1,2)

	Cd(5123) = 2*Cd(5097) - MGl2*Cd(5122)

	Cd(5124) = Cd(5103) - MGl2*Cd(5123)

	Cd(5125) = 2*Cd(44)**2*Cd(48)**2*Cd(5084) - MGl2*Cd(5124)

	Cd(5128) = Cd(5113) + MGl2*Cd(5127)

	Cd(5129) = 2*Cd(5118) + MGl2*Cd(5128)

	Cd(5130) = Cd(44)**2*Cd(48)**2*Cd(288) - MGl2*Cd(5129)

	Cd(5131) = Cd(5125)*UStop(1,1) - 8*M_3*MTy*Cd(5130)*UStop(1,2)

	Cd(5139) = Cd(5137) - MGl2*Cd(5138)

	Cd(5140) = Cd(48)*Cd(5134) - 2*MGl2*Cd(5139)

	Cd(5161) = Cd(48)*Cd(5151) + MGl2*Cd(5160)

	Cd(5162) = Cd(44)*Cd(48)**3*Cd(187) + MGl2*Cd(5161)

	Cd(5165) = 2*Cd(5157) - MGl2*Cd(5164)

	Cd(5166) = Cd(48)*Cd(5159) + MGl2*Cd(5165)

        Cd(5167) = 2*(Cd(5162)*UStop(1,1) + 
     &     M_3*MTy*Cd(5166)*UStop(1,2))

        Cd(5218) = 2*M_3C*Cd(47)*Cd(5140)*UCStop(1,3)*
     &    UStop(1,1) - 4*M_3*Cd(5217)*UStop(1,2)

	Cd(517) = Cd(424)*Cd(507)*MStop2(3) - Cd(516)*UStop2(1,1)

	Cd(5257) = Cd(5244) - MGl2*Cd(5256)

	Cd(5258) = Cd(5250) - MGl2*Cd(5257)

	Cd(5259) = Cd(44)**2*Cd(48)**2*Cd(89)*Cd(1847) - MGl2*Cd(5258)

	Cd(5262) = 2*Cd(4929) - MGl2*Cd(5261)

	Cd(5263) = Cd(4935) + MGl2*Cd(5262)

	Cd(5264) = Cd(44)*Cd(48)*Cd(5252) + MGl2*Cd(5263)

	Cd(5279) = Cd(5275) - MGl2*Cd(5278)

	Cd(5280) = 2*Cd(44)*Cd(48)*Cd(4995) - 2*MGl2*Cd(5279)

	Cd(5284) = -4*Cd(47)*Cd(50)*Cd(295) - Cd(4941)*UStop2(1,2)

        Cd(5285) = -4*Cd(47)*Cd(50)*Cd(295)*Cd(1523) - 
     &   Cd(5284)*UStop2(1,1)

	Cd(5304) = Cd(5301) + MGl2*Cd(5303)

	Cd(5305) = Cd(5300) + MGl2*Cd(5304)

	Cd(5306) = 2*Cd(44)**2*Cd(48)**2*Cd(5286) - 2*MGl2*Cd(5305)

	Cd(5307) = 4*Cd(47)*Cd(50)*Cd(282) - Cd(5306)*UStop2(1,2)

	Cd(5311) = Cd(5295) + MGl2*Cd(5310)

	Cd(5312) = Cd(5300) + MGl2*Cd(5311)

	Cd(5313) = Cd(44)**2*Cd(48)**2*Cd(5286) - MGl2*Cd(5312)

	Cd(5314) = Cd(5307) - 2*Cd(5313)*UStop2(1,1)

	Cd(5316) = -Cd(5275) + MGl2*Cd(5315)

	Cd(5317) = Cd(44)*Cd(48)*Cd(4995) + MGl2*Cd(5316)

        Cd(5318) = Cd(5314)*MStop2(3) - 
     &   4*M_3*MTy*Cd(5317)*UCStopC(1,3)

	Cd(5319) = 2*M_3C*MTy*Cd(5285) - Cd(5318)*UCStopC(1,3)

	Cd(5347) = 2*Cd(5337) - MGl2*Cd(5346)

	Cd(5348) = Cd(44)**2*Cd(48)**2*Cd(5322) - MGl2*Cd(5347)

        Cd(5435) = MTy2*Cd(5434) + 
     &   5*MStop2(1)*MStop2(2)*MStop2(3)**2*MStop2(4)

	Cd(5443) = 4*MTy2*Cd(5440) + 10*Cd(5442)*MStop2(3)**2*MStop2(4)

        Cd(5450) = 2*(MTy2*Cd(5449) + 
     &     5*Cd(5444)*MStop2(1)*MStop2(2)*MStop2(3)**2)

	Cd(5467) = 2*Cd(5430) + MGl2*Cd(5466)

	Cd(5468) = Cd(5443) - 2*MGl2*Cd(5467)

	Cd(5485) = Cd(5472) - 2*MGl2*Cd(5484)

	Cd(5486) = 2*Cd(5476) + MGl2*Cd(5485)

	Cd(5487) = Cd(5482) - MGl2*Cd(5486)

	Cd(5488) = Cd(44)*Cd(48)*Cd(5470) - MGl2*Cd(5487)

        Cd(5509) = 2*(MTy2*Cd(5508) + 
     &     5*Cd(5444)*MStop2(1)*MStop2(2)*MStop2(3)**2)

	Cd(5515) = 2*Cd(5509) + MGl2*Cd(5514)

	Cd(5516) = 2*Cd(44)*Cd(48)*Cd(5497) - MGl2*Cd(5515)

	Cd(5521) = -(MGl2*Cd(4939)) + Cd(5520)

	Cd(5522) = Cd(44)*Cd(48)*Cd(5252) - MGl2*Cd(5521)

	Cd(5523) = Cd(5516) + 4*MTy2*Cd(5522)*UStop2(1,2)

	Cd(5524) = Cd(5488)*MStop2(3)**2 + 2*MTy2*Cd(5523)*UStop2(1,2)

	Cd(5566) = Cd(5542) - MGl2*Cd(5565)

	Cd(5567) = Cd(5548) - 2*MGl2*Cd(5566)

	Cd(5568) = 2*Cd(44)**2*Cd(48)**2*Cd(5529) - MGl2*Cd(5567)

	Cd(5573) = Cd(5561) + MGl2*Cd(5572)

	Cd(5574) = 2*Cd(44)*Cd(48)*Cd(89)**2*Cd(5556) + MGl2*Cd(5573)

	Cd(5575) = Cd(5568)*UStop(1,1) - M_3*MTy*Cd(5574)*UStop(1,2)

	Cd(5593) = 4*Cd(5588) - MGl2*Cd(5592)

	Cd(5594) = 2*Cd(44)**2*Cd(48)**2*Cd(89)**2 - MGl2*Cd(5593)

        Cd(5595) = 3*Cd(5594)*MStop2(3)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(5264)*UStop(1,2)

	Cd(5624) = 2*Cd(5606) - MGl2*Cd(5623)

	Cd(5625) = Cd(5611) - MGl2*Cd(5624)

	Cd(5626) = Cd(5616) - MGl2*Cd(5625)

	Cd(5627) = Cd(44)**2*Cd(48)**2*Cd(5598) - MGl2*Cd(5626)

	Cd(5629) = -Cd(5472) + MGl2*Cd(5628)

	Cd(5630) = Cd(5619) - MGl2*Cd(5629)

	Cd(5631) = Cd(5482) - MGl2*Cd(5630)

	Cd(5632) = Cd(44)*Cd(48)*Cd(5470) - MGl2*Cd(5631)

        Cd(5633) = 2*MTy*Cd(5627)*UStop(1,1) - 
     &   M_3*Cd(5632)*MStop2(3)*UStop(1,2)

	Cd(5635) = 2*Cd(198)*Cd(5553) - MGl2*Cd(5634)

	Cd(5636) = Cd(5555) - MGl2*Cd(5635)

	Cd(5637) = Cd(5561) + 2*MGl2*Cd(5636)

	Cd(5638) = 2*Cd(44)*Cd(48)*Cd(89)**2*Cd(5556) + MGl2*Cd(5637)

        Cd(5679) = MTy2*Cd(5678) + 
     &   3*Cd(5601)*MStop2(1)*MStop2(2)*MStop2(3)*MStop2(4)

        Cd(560) = (Cd(559)*MStop2(1) + 
     &     2*Cd(558)*MStop2(3)**2*UStop2(1,1)*UStopC(1,2))/
     &   MStop2(3)**2

	Cd(5688) = 2*MTy2*Cd(5686) - Cd(5687)*MStop2(4)

        Cd(562) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*(Cd(561) + CA2*Cd(560)*MStop2(3)**2))/
     &    (SB2*SW2*MStop2(3)**2)

	Cd(5733) = Cd(5730) - MGl2*Cd(5732)

	Cd(5734) = Cd(5722) + MGl2*Cd(5733)

	Cd(5735) = Cd(5715) + MGl2*Cd(5734)

	Cd(5742) = Cd(5653) - MGl2*Cd(5741)

	Cd(5743) = Cd(5658) + MGl2*Cd(5742)

	Cd(5744) = Cd(5679) - MGl2*Cd(5743)

	Cd(5745) = Cd(5673) - MGl2*Cd(5744)

	Cd(5746) = Cd(44)**2*Cd(48)**2*Cd(5649) - MGl2*Cd(5745)

	Cd(5748) = 4*Cd(5684) - MGl2*Cd(5747)

	Cd(5749) = Cd(5688) - MGl2*Cd(5748)

	Cd(5750) = 2*Cd(5698) - MGl2*Cd(5749)

	Cd(5751) = Cd(44)*Cd(48)*Cd(5702) - MGl2*Cd(5750)

	Cd(5775) = Cd(5769) - MGl2*Cd(5774)

	Cd(5776) = Cd(44)*Cd(48)*Cd(5764) - MGl2*Cd(5775)

	Cd(5801) = 2*Cd(5789) - MGl2*Cd(5800)

	Cd(5802) = Cd(5795) - MGl2*Cd(5801)

	Cd(5839) = 2*Cd(47)*Cd(50)*Cd(5805) - 2*Cd(5838)*UStop2(1,1)

	Cd(5977) = 2*Cd(5959) - MGl2*Cd(5976)

	Cd(5983) = Cd(44)*Cd(48)*Cd(5980) - MGl2*Cd(5982)

	Cd(5986) = Cd(5845) - 2*MGl2*Cd(5985)

	Cd(5987) = 2*Cd(5849) + MGl2*Cd(5986)

	Cd(5988) = Cd(5854) - MGl2*Cd(5987)

	Cd(5989) = Cd(5860) - MGl2*Cd(5988)

	Cd(5990) = Cd(44)**2*Cd(48)**2*Cd(89)*Cd(2521) - MGl2*Cd(5989)

	Cd(6017) = -2*(2*Cd(47)*Cd(50)*Cd(295) + Cd(5522)*UStop2(1,2))

	Cd(6045) = Cd(6035) - MGl2*Cd(6044)

	Cd(6046) = Cd(6040) - 2*MGl2*Cd(6045)

	Cd(6047) = Cd(44)**2*Cd(48)**2*Cd(6022) - MGl2*Cd(6046)

        Cd(6048) = Cd(6047)*UStop(1,1) - 
     &   16*M_3*MTy*Cd(5130)*UStop(1,2)

        Cd(6090) = -2*M_3C*Cd(47)*Cd(50)*Cd(295)*UCStop(1,3)*
     &    UStop(1,1) + M_3*Cd(6087)*Cd(6088)*Cd(6089)*UStop(1,2)

        Cd(6093) = M_3C*Cd(6017)*UCStop(1,3)*UStop2(1,1) + 
     &   4*M_3*Cd(47)*Cd(50)*Cd(295)*UCStopC(1,3)*UStop2(1,2)

	Cd(6106) = 2*MGl2*Cd(6095) + Cd(6099)

        Cd(6107) = (-4*Cd(295) - 
     &     Cd(47)*Cd(50)*Cd(6106)*MStop2(3)*UStop2(1,2))/
     &   (Cd(47)*Cd(50)*MStop2(3))

	Cd(6115) = Cd(5244) - MGl2*Cd(6114)

	Cd(6116) = Cd(5250) - MGl2*Cd(6115)

	Cd(6117) = Cd(44)**2*Cd(48)**2*Cd(89)*Cd(1847) - MGl2*Cd(6116)

        Cd(6119) = 4*Cd(47)*Cd(50)*Cd(6118) - 
     &   2*Cd(6117)*UStop(1,1)*UStop2(1,2)

	Cd(6154) = 2*Cd(6142) - MGl2*Cd(6153)

	Cd(6155) = Cd(6149) - MGl2*Cd(6154)

	Cd(6156) = 2*Cd(44)*Cd(48)*Cd(6144) - 2*MGl2*Cd(6155)

	Cd(6165) = Cd(5789) - MGl2*Cd(6164)

	Cd(6166) = Cd(6163) - MGl2*Cd(6165)

	Cd(6167) = Cd(44)*Cd(48)*Cd(5785) - MGl2*Cd(6166)

	Cd(6200) = 4*Cd(6173) - MGl2*Cd(6199)

	Cd(6201) = Cd(6178) + MGl2*Cd(6200)

	Cd(6206) = 2*Cd(6186) - MGl2*Cd(6205)

	Cd(6207) = 2*Cd(6193) + MGl2*Cd(6206)

	Cd(6208) = 2*Cd(6198) - MGl2*Cd(6207)

	Cd(6209) = 2*Cd(44)**2*Cd(48)**2*Cd(89)*Cd(2521) - MGl2*Cd(6208)

	Cd(6235) = Cd(5854) - MGl2*Cd(6234)

	Cd(6236) = Cd(5860) - MGl2*Cd(6235)

        Cd(6237) = 2*Cd(44)**2*Cd(48)**2*Cd(89)*Cd(2521) - 
     &   2*MGl2*Cd(6236)

	Cd(6255) = Cd(44)*Cd(48)*Cd(6252) - MGl2*Cd(6254)

	Cd(6336) = Cd(5157) - MGl2*Cd(6335)

	Cd(6337) = Cd(48)*Cd(5159) + 2*MGl2*Cd(6336)

	Cd(6371) = 2*Cd(6351) - MGl2*Cd(6370)

	Cd(6372) = Cd(6359) - 2*MGl2*Cd(6371)

	Cd(6373) = 2*Cd(48)*Cd(6356) + MGl2*Cd(6372)

	Cd(6374) = 2*Cd(44)**2*Cd(48)**3*Cd(1644) - MGl2*Cd(6373)

	Cd(6375) = 2*Cd(6367) - Cd(6374)*UStop2(1,1)

	Cd(6376) = 2*MGl2*Cd(5316) + Cd(44)*Cd(48)*Cd(6369)

        Cd(6377) = Cd(6375)*MStop2(3) - 
     &   2*M_3*MTy*Cd(6376)*UCStopC(1,3)

        Cd(6390) = Cd(47)*Cd(6389)*MStop2(3) - 
     &   4*MTy2*Cd(5522)*UStop2(1,2)

	Cd(6418) = Cd(6400) + MGl2*Cd(6417)

	Cd(6419) = Cd(6406) - 2*MGl2*Cd(6418)

	Cd(6420) = 2*Cd(44)*Cd(6412) - MGl2*Cd(6419)

	Cd(6429) = Cd(6420) - MTy2*Cd(50)*Cd(6428)*UStop2(1,2)

	Cd(6430) = 2*Cd(47)**2*Cd(6391)*MStop2(3) - Cd(6429)*UStop2(1,2)

        Cd(6431) = MTy2*Cd(6383)*UStop2(1,1)**3 + 
     &   Cd(50)*Cd(282)*UStop2(1,2)**2

        Cd(62) = Cd(52)*Cd(55)*MStop2(3) + 
     &   2*Cd(61)*T134(MGl2,MTy2,MStop2(1))

	Cd(6502) = Cd(6486) - MGl2*Cd(6501)

	Cd(6503) = Cd(6491) - MGl2*Cd(6502)

	Cd(6504) = Cd(44)*Cd(48)*Cd(6488) - MGl2*Cd(6503)

	Cd(6505) = Cd(6498)*UStop(1,1) + 4*M_3*MTy*Cd(6504)*UStop(1,2)

        Cd(6550) = MTy2*Cd(6549) - 
     &   Cd(1272)*Cd(5828)*MStop2(1)*MStop2(2)*MStop2(3)

	Cd(6577) = Cd(6559) - MGl2*Cd(6576)

	Cd(6578) = Cd(6563) + MGl2*Cd(6577)

	Cd(6579) = Cd(6572) - MGl2*Cd(6578)

	Cd(6580) = Cd(48)*Cd(6567) - MGl2*Cd(6579)

	Cd(6581) = Cd(44)**2*Cd(48)**3*Cd(6552) - MGl2*Cd(6580)

	Cd(6584) = Cd(6491) - MGl2*Cd(6583)

	Cd(6585) = Cd(44)*Cd(48)*Cd(6488) - MGl2*Cd(6584)

	Cd(6586) = Cd(6581)*UStop(1,1) + M_3*MTy*Cd(6585)*UStop(1,2)

	Cd(6588) = Cd(6518) + MGl2*Cd(6587)

	Cd(6589) = Cd(6525) - MGl2*Cd(6588)

	Cd(6590) = Cd(1644)*Cd(6530) - MGl2*Cd(6589)

	Cd(6592) = Cd(6537) - MGl2*Cd(6591)

	Cd(6593) = Cd(6540) - MGl2*Cd(6592)

	Cd(6594) = Cd(6550) - MGl2*Cd(6593)

	Cd(6596) = Cd(6460)*MStop2(3) + Cd(6505)*UStop2(1,1)*UStopC(1,1)

	Cd(6646) = 2*Cd(6629) - MGl2*Cd(6645)

	Cd(6647) = 2*Cd(6632) + MGl2*Cd(6646)

	Cd(6648) = Cd(6635) - MGl2*Cd(6647)

	Cd(6649) = Cd(44)*Cd(48)*Cd(6637) - MGl2*Cd(6648)

	Cd(6650) = Cd(6644)*UStop(1,1) - 2*M_3*MTy*Cd(6649)*UStop(1,2)

	Cd(6652) = -(Cd(47)*Cd(6651)*UStop(1,1)) + Cd(6650)*UStop2(1,2)

	Cd(6659) = Cd(6635) - MGl2*Cd(6658)

	Cd(6660) = Cd(44)*Cd(48)*Cd(6637) - MGl2*Cd(6659)

        Cd(6668) = Cd(6667)*MStop2(3)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(2317)*UStop(1,2)

        Cd(6671) = (2*MGl2*MTy*Cd(47)*MStop2(3)*UStop(1,1) + 
     &     M_3*Cd(6670)*UStop(1,2))/MGl

        Cd(6672) = 4*MGl*MTy*Cd(6671) - 
     &   Cd(47)*MGlpTmSt4(2)*MStop2(3)*UStop(1,1)

        Cd(6675) = 4*M_3C*MTy*Cd(6674)*UStop(1,1) + 
     &   Cd(282)*MStop2(3)*UStop(1,2)*UStop2(1,2)

	Cd(6676) = Cd(6652) + Cd(47)*Cd(50)*Cd(6603)*UStop2(1,1)

	Cd(6703) = Cd(6692) - MGl2*Cd(6702)

	Cd(6704) = Cd(44)*Cd(48)*Cd(6694) - MGl2*Cd(6703)

	Cd(6705) = -(MGl2*Cd(6370)) + 2*Cd(6695)

	Cd(6706) = Cd(6698) - MGl2*Cd(6705)

	Cd(6707) = 2*Cd(48)*Cd(6696) - MGl2*Cd(6706)

	Cd(6708) = Cd(44)**2*Cd(48)**3*Cd(1644) + MGl2*Cd(6707)

        Cd(6787) = MTy2*Cd(6786) - 
     &   2*Cd(832)*MStop2(1)*MStop2(2)*MStop2(3)*MStop2(4)

	Cd(6790) = Cd(6752) - 2*MGl2*Cd(6789)

	Cd(6791) = 2*Cd(6758) + MGl2*Cd(6790)

	Cd(6792) = Cd(6762) - MGl2*Cd(6791)

	Cd(6793) = 2*Cd(6768) - MGl2*Cd(6792)

	Cd(6794) = Cd(44)**2*Cd(48)**2*Cd(6755) - MGl2*Cd(6793)

	Cd(6796) = 2*Cd(6774) - MGl2*Cd(6795)

	Cd(6797) = Cd(6778) - MGl2*Cd(6796)

	Cd(6798) = Cd(6787) - MGl2*Cd(6797)

	Cd(6799) = Cd(44)*Cd(48)*Cd(6782) - MGl2*Cd(6798)

        Cd(6800) = 2*M_3*Cd(6799)*UCStopC(1,3) + 
     &   MTy*Cd(6794)*UStop2(1,1)

	Cd(6837) = 2*Cd(48)**2*Cd(6812) - MGl2*Cd(6836)

	Cd(6838) = Cd(44)**2*Cd(48)**4*Cd(6823) - MGl2*Cd(6837)

	Cd(6840) = 2*Cd(6829) - MGl2*Cd(6839)

	Cd(6841) = Cd(6832) - MGl2*Cd(6840)

	Cd(6842) = Cd(6835) - MGl2*Cd(6841)

	Cd(6843) = Cd(44)**2*Cd(48)*Cd(6826) - MGl2*Cd(6842)

	Cd(6844) = Cd(6838)*UStop(1,1) + 2*M_3*MTy*Cd(6843)*UStop(1,2)

        Cd(6845) = Cd(47)**2*Cd(6745)*MStop2(3) + 
     &   Cd(6844)*UStop2(1,1)*UStopC(1,1)

	Cd(6846) = Cd(6800) - MTy*Cd(6845)

	Cd(6882) = Cd(6854) - MGl2*Cd(6881)

	Cd(6883) = Cd(6862) + MGl2*Cd(6882)

	Cd(6884) = Cd(6867) - MGl2*Cd(6883)

	Cd(6885) = Cd(48)**2*Cd(6857) - MGl2*Cd(6884)

	Cd(6886) = Cd(44)**2*Cd(48)**4*Cd(2331) - MGl2*Cd(6885)

	Cd(6888) = 2*Cd(6872) - MGl2*Cd(6887)

	Cd(6889) = Cd(6875) - MGl2*Cd(6888)

	Cd(6890) = Cd(6880) - MGl2*Cd(6889)

	Cd(6891) = Cd(44)*Cd(48)*Cd(6877) - MGl2*Cd(6890)

	Cd(6892) = Cd(6886)*UStop(1,1) + M_3*MTy*Cd(6891)*UStop(1,2)

	Cd(6920) = 4*Cd(6897) - MGl2*Cd(6919)

	Cd(6921) = Cd(6906) - MGl2*Cd(6920)

	Cd(6922) = 2*Cd(6902) - MGl2*Cd(6921)

	Cd(6923) = Cd(44)**2*Cd(48)**2*Cd(288) - MGl2*Cd(6922)

	Cd(6925) = Cd(6910) - MGl2*Cd(6924)

	Cd(6926) = Cd(6914) - MGl2*Cd(6925)

	Cd(6927) = Cd(6918) + MGl2*Cd(6926)

	Cd(6928) = Cd(44)**2*Cd(48)**2*Cd(6907) - MGl2*Cd(6927)

	Cd(6929) = MTy*Cd(6923)*UStop(1,1) + M_3*Cd(6928)*UStop(1,2)

        Cd(6942) = MTy2*Cd(6941) + 
     &   Cd(902)*MStop2(1)*MStop2(2)*MStop2(3)*MStop2(4)

	Cd(6956) = Cd(6951) - MGl2*Cd(6955)

	Cd(6957) = Cd(44)*Cd(48)*Cd(6944) - MGl2*Cd(6956)

	Cd(6964) = 2*MGl2**3*Cd(6963) + Cd(6957)*UStop2(1,1)

	Cd(6966) = Cd(6929)*MStop2(3) - 2*MTy*Cd(6892)*UStop2(1,1)

	Cd(6972) = Cd(6968) - MGl2*Cd(6971)

	Cd(6973) = Cd(6951) - MGl2*Cd(6972)

	Cd(6974) = 2*Cd(44)*Cd(48)*Cd(6944) - 2*MGl2*Cd(6973)

        Cd(6981) = Cd(50)*Cd(6980)*MStop2(3) + 
     &   M_3C*MTy*Cd(6974)*UUStop(1,1)

        Cd(6994) = -3*Cd(47)**2*Cd(5013)*UStop(1,1) + 
     &   Cd(6650)*UStop2(1,2)

        Cd(6995) = 2*(Cd(6994) - 
     &     M_3C*MTy*Cd(47)*Cd(6993)*UCStop(1,3)*UStop(1,1))

	Cd(6996) = Cd(47)*Cd(6672) + Cd(50)*Cd(6668)*UStop2(1,2)

	Cd(7003) = Cd(7000) - MGl2*Cd(7002)

	Cd(7004) = Cd(48)*Cd(6999) - 2*MGl2*Cd(7003)

	Cd(7009) = Cd(7006) - MGl2*Cd(7008)

	Cd(7010) = Cd(44)*Cd(7005) + 2*MGl2*Cd(7009)

	Cd(7011) = Cd(47)*Cd(7004) - Cd(50)*Cd(7010)*UStop2(1,2)

        Cd(7012) = 2*Cd(6996)*UStop(1,2) - 
     &   2*M_3C*MTy*Cd(7011)*UUStop(1,1)

	Cd(7013) = Cd(6995) + 2*Cd(47)*Cd(50)*Cd(6603)*UStop2(1,1)

	Cd(7020) = -(MGl2*Cd(5635)) + 2*Cd(7018)

	Cd(7021) = Cd(7019) - MGl2*Cd(7020)

	Cd(7022) = Cd(44)*Cd(48)*Cd(89)**2*Cd(5556) - MGl2*Cd(7021)

        Cd(7023) = -2*Cd(47)*Cd(50)*Cd(295)*Cr(28) + 
     &   Cd(7022)*UStop2(1,1)*UStop2(1,2)

	Cd(7027) = Cd(44)**2*Cd(48)**2 + MGl2*Cd(7026)

	Cd(7037) = Cd(44)**2*Cd(48)**2 - MGl2*Cd(7036)

	Cd(7040) = 4*Cd(7029) - MGl2*Cd(7039)

	Cd(7041) = Cd(7031) - MGl2*Cd(7040)

	Cd(7082) = Cd(7067) + MGl2*Cd(7081)

	Cd(7083) = Cd(7072) - MGl2*Cd(7082)

	Cd(7084) = Cd(7078) - MGl2*Cd(7083)

	Cd(7085) = 3*Cd(44)**2*Cd(48)**2*Cd(7064) - MGl2*Cd(7084)

	Cd(7086) = Cd(7085)*UStop(1,1) - 4*M_3*MTy*Cd(5130)*UStop(1,2)

	Cd(7102) = Cd(5162)*UStop(1,1) - 2*M_3*MTy*Cd(7101)*UStop(1,2)

        Cd(7144) = M_3C*Cd(47)*Cd(7093)*UCStop(1,3)*
     &    UStop(1,1) + M_3*Cd(5217)*UStop(1,2)

        Cd(7164) = Cd(7160)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(7163)*MStop2(3)*UStop(1,2)

	Cd(7196) = Cd(7184) - MGl2*Cd(7195)

	Cd(7197) = 2*Cd(7190) - MGl2*Cd(7196)

	Cd(7198) = Cd(44)**2*Cd(48)**2*Cd(7175) - MGl2*Cd(7197)

	Cd(7202) = Cd(7019) - MGl2*Cd(7201)

	Cd(7203) = Cd(44)*Cd(48)*Cd(89)**2*Cd(5556) - MGl2*Cd(7202)

	Cd(7227) = Cd(7222) - MGl2*Cd(7226)

	Cd(7228) = Cd(44)*Cd(48)*Cd(7212) - MGl2*Cd(7227)

	Cd(7232) = -2*Cd(47)*Cd(50)*Cd(295) - Cd(7022)*UStop2(1,2)

        Cd(7233) = 2*Cd(47)*Cd(50)*Cd(295)*Cd(1523) + 
     &   Cd(7232)*UStop2(1,1)

	Cd(7250) = -Cd(7243) + MGl2*Cd(7249)

	Cd(7251) = 2*Cd(7247) - MGl2*Cd(7250)

	Cd(7252) = Cd(44)**2*Cd(48)**2*Cd(7237) - MGl2*Cd(7251)

	Cd(7253) = 2*Cd(47)*Cd(50)*Cd(282) - Cd(7252)*UStop2(1,2)

	Cd(7262) = Cd(7243) - MGl2*Cd(7261)

	Cd(7263) = 2*Cd(7247) + MGl2*Cd(7262)

	Cd(7264) = Cd(44)**2*Cd(48)**2*Cd(7237) - MGl2*Cd(7263)

	Cd(7265) = Cd(7253) - Cd(7264)*UStop2(1,1)

	Cd(7269) = Cd(7259) - MGl2*Cd(7268)

	Cd(7270) = Cd(44)*Cd(48)*Cd(7212) - MGl2*Cd(7269)

        Cd(7271) = Cd(7265)*MStop2(3) - 
     &   2*M_3*MTy*Cd(7270)*UCStopC(1,3)

	Cd(7272) = 2*M_3C*MTy*Cd(7233) + Cd(7271)*UCStopC(1,3)

	Cd(7298) = Cd(5365) + MGl2*Cd(7297)

	Cd(7299) = 2*Cd(44)*Cd(48)*Cd(5352) + MGl2*Cd(7298)

	Cd(7300) = Cd(44)**2*Cd(48)**2*Cd(288) + MGl2*Cd(7299)

	Cd(7304) = Cd(5360) + MGl2*Cd(7303)

	Cd(7305) = 2*Cd(5369) - MGl2*Cd(7304)

        Cd(7352) = MTy2*Cd(7351) + 
     &   5*MStop2(1)*MStop2(2)*MStop2(3)**2*MStop2(4)

	Cd(725) = Cd(723) - Cd(689)*Cd(691)*UStopC(2,1)

	Cd(7357) = 2*MTy2*Cd(7356) + 5*Cd(5442)*MStop2(3)**2*MStop2(4)

	Cd(7369) = 4*Cd(7348) + MGl2*Cd(7368)

	Cd(7370) = Cd(7357) - MGl2*Cd(7369)

        Cd(727) = -(Cd(726)*Conjugate(dY33eps(-1))) + 
     &   Cd(725)*dY33eps(-1)

        Cd(7400) = MTy2*Cd(7399) - 
     &   5*MStop2(1)*MStop2(2)*MStop2(3)**2*MStop2(4)

	Cd(730) = -(Cd(726)*Conjugate(dY33eps(0))) + Cd(725)*dY33eps(0)

	Cd(7405) = 2*MTy2*Cd(7404) - 5*Cd(5442)*MStop2(3)**2*MStop2(4)

        Cd(7411) = MTy2*Cd(7410) - 
     &   5*Cd(5444)*MStop2(1)*MStop2(2)*MStop2(3)**2

	Cd(7415) = 4*Cd(7395) + MGl2*Cd(7414)

	Cd(7416) = Cd(7405) - MGl2*Cd(7415)

	Cd(7417) = 2*Cd(7411) - MGl2*Cd(7416)

	Cd(7418) = Cd(44)*Cd(48)*Cd(7400) - MGl2*Cd(7417)

	Cd(7419) = -Cd(7019) + MGl2*Cd(7020)

	Cd(7420) = Cd(44)*Cd(48)*Cd(89)**2*Cd(5556) + MGl2*Cd(7419)

	Cd(7421) = Cd(7418) - 2*MTy2*Cd(7420)*UStop2(1,2)

	Cd(7422) = Cd(7391)*MStop2(3)**2 - 2*MTy2*Cd(7421)*UStop2(1,2)

	Cd(7437) = 2*Cd(7432) - MGl2*Cd(7436)

	Cd(7438) = Cd(44)**2*Cd(48)**2*Cd(7423) - MGl2*Cd(7437)

        Cd(7439) = 3*Cd(7438)*MStop2(3)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(7203)*UStop(1,2)

	Cd(7470) = Cd(7445) - MGl2*Cd(7469)

	Cd(7471) = Cd(7452) - MGl2*Cd(7470)

	Cd(7472) = Cd(7457) + MGl2*Cd(7471)

	Cd(7473) = Cd(7463) - MGl2*Cd(7472)

	Cd(7474) = 3*Cd(44)**2*Cd(48)**2*Cd(7447) - MGl2*Cd(7473)

	Cd(7479) = Cd(44)*Cd(48)*Cd(5252) - MGl2*Cd(7478)

	Cd(7480) = Cd(7474)*UStop(1,1) - M_3*MTy*Cd(7479)*UStop(1,2)

	Cd(739) = 2*Cd(738)*MStop2(1) - 2*Cd(360)*Cd(704)*MStop2(2)

	Cd(740) = Cd(739) + 2*Cd(734)*MStop2(1)*UStopC(2,1)

	Cd(7519) = Cd(7508) + MGl2*Cd(7518)

	Cd(7520) = Cd(7385) - MGl2*Cd(7519)

	Cd(742) = Cd(741)*Conjugate(dY33eps(-1)) + Cd(740)*dY33eps(-1)

	Cd(7521) = Cd(44)*Cd(48)*Cd(7372) - MGl2*Cd(7520)

        Cd(7522) = 2*MTy*Cd(7517)*UStop(1,1) - 
     &   M_3*Cd(7521)*MStop2(3)*UStop(1,2)

	Cd(7524) = 2*Cd(4929) - MGl2*Cd(7523)

	Cd(7525) = Cd(5520) - MGl2*Cd(7524)

	Cd(7526) = Cd(44)*Cd(48)*Cd(5252) - MGl2*Cd(7525)

        Cd(744) = Cd(743)*Conjugate(dY33eps(-1)) + 
     &   Cd(360)*Cd(704)*dY33eps(-1)

	Cd(7546) = Cd(7540) + MGl2*Cd(7545)

	Cd(7550) = Cd(5658) + MGl2*Cd(7549)

	Cd(7551) = Cd(5679) - MGl2*Cd(7550)

	Cd(7552) = Cd(5673) - MGl2*Cd(7551)

        Cd(746) = Cd(743)*Conjugate(dY33eps(0)) + 
     &   Cd(360)*Cd(704)*dY33eps(0)

	Cd(7576) = Cd(7569) - MGl2*Cd(7575)

	Cd(7577) = Cd(44)*Cd(48)*Cd(7560) - MGl2*Cd(7576)

	Cd(7599) = Cd(7589) - MGl2*Cd(7598)

	Cd(7600) = Cd(7595) - MGl2*Cd(7599)

	Cd(7605) = Cd(7602) + MGl2*Cd(7604)

	Cd(7606) = Cd(44)*Cd(48)*Cd(89)**2 + MGl2*Cd(7605)

        Cd(7626) = 3*Cd(1272)*Cd(5397)*MStop2(1)*MStop2(2) + 
     &   Cd(302)*Cd(7606)*UStop2(1,2)**2

	Cd(7658) = 4*Cd(198)*Cd(7649) - MGl2*Cd(7657)

	Cd(7659) = 3*Cd(7652) - MGl2*Cd(7658)

	Cd(7660) = 2*Cd(7655) - MGl2*Cd(7659)

	Cd(7661) = Cd(44)**2*Cd(48)**2*Cd(7648) - MGl2*Cd(7660)

	Cd(7724) = Cd(7714) - MGl2*Cd(7723)

	Cd(7725) = 2*Cd(7719) - MGl2*Cd(7724)

	Cd(7726) = Cd(44)**2*Cd(48)**2*Cd(7647) - MGl2*Cd(7725)

        Cd(7730) = Cd(7726)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(7729)*MStop2(3)*UStop(1,2)

	Cd(7738) = 2*Cd(5959) - MGl2*Cd(7737)

	Cd(7739) = Cd(44)**2*Cd(48)**2 - MGl2*Cd(7738)

	Cd(7747) = Cd(7640) - MGl2*Cd(7746)

	Cd(7748) = 2*Cd(7646) - MGl2*Cd(7747)

	Cd(7749) = Cd(44)**2*Cd(48)**2*Cd(89)*Cd(7647) - MGl2*Cd(7748)

	Cd(7768) = Cd(7759) + MGl2*Cd(7767)

	Cd(7769) = Cd(7764) - MGl2*Cd(7768)

	Cd(7770) = 3*Cd(44)**2*Cd(48)**2*Cd(6438) - MGl2*Cd(7769)

	Cd(7771) = -2*Cd(47)*Cd(50)*Cd(295) - Cd(7420)*UStop2(1,2)

	Cd(7797) = Cd(7787) - MGl2*Cd(7796)

	Cd(7798) = Cd(7792) - MGl2*Cd(7797)

	Cd(7799) = Cd(44)**2*Cd(48)**2*Cd(7778) - MGl2*Cd(7798)

	Cd(7800) = Cd(7799)*UStop(1,1) - 8*M_3*MTy*Cd(5130)*UStop(1,2)

	Cd(7840) = 3*Cd(47)*Cd(50)*Cd(282) - Cd(7770)*UStop2(1,2)**2

        Cd(7841) = 2*M_3C*Cd(7771)*UCStop(1,3)*UStop2(1,1) + 
     &   4*M_3*Cd(47)*Cd(50)*Cd(295)*UCStopC(1,3)*UStop2(1,2)

	Cd(7851) = MGl2*Cd(6095) + Cd(7844)

        Cd(7852) = (-2*Cd(295) + 
     &     Cd(47)*Cd(50)*Cd(7851)*MStop2(3)*UStop2(1,2))/
     &   (Cd(47)*Cd(50)*MStop2(3))

	Cd(7859) = -Cd(7184) + MGl2*Cd(7858)

        Cd(75) = Cd(64)*Cd(66)*MStop2(3) + 
     &   Cd(74)*T134(MGl2,MTy2,MStop2(2))

	Cd(7860) = 2*Cd(7190) + MGl2*Cd(7859)

	Cd(7861) = Cd(44)**2*Cd(48)**2*Cd(7175) - MGl2*Cd(7860)

        Cd(7862) = 2*Cd(47)*Cd(50)*Cd(6118) - 
     &   Cd(7861)*UStop(1,1)*UStop2(1,2)

	Cd(7897) = 2*Cd(7887) + MGl2*Cd(7896)

	Cd(7898) = Cd(7891) - MGl2*Cd(7897)

	Cd(7899) = Cd(44)*Cd(48)*Cd(7882) - MGl2*Cd(7898)

	Cd(7904) = Cd(7874) - MGl2*Cd(7903)

	Cd(7905) = Cd(7880) - MGl2*Cd(7904)

	Cd(7910) = Cd(7589) - MGl2*Cd(7909)

	Cd(7911) = Cd(7595) - MGl2*Cd(7910)

	Cd(7912) = Cd(44)*Cd(48)*Cd(7585) - MGl2*Cd(7911)

	Cd(7957) = Cd(7640) - MGl2*Cd(7956)

	Cd(7958) = 2*Cd(7646) - MGl2*Cd(7957)

	Cd(7959) = Cd(44)**2*Cd(48)**2*Cd(89)*Cd(7647) - MGl2*Cd(7958)

	Cd(76) = Cd(50)*Cd(62) - Cd(47)*Cd(75)

	Cd(8037) = Cd(47)**2*MGlpTmSt4(2) + Cd(50)*Cd(7027)*UStop2(1,2)

	Cd(8043) = Cd(44)**2*Cd(48)**3*Cd(8036) - MGl2*Cd(8042)

	Cd(8048) = Cd(7259) - MGl2*Cd(8047)

	Cd(8049) = Cd(44)*Cd(48)*Cd(8041) - MGl2*Cd(8048)

        Cd(8060) = Cd(47)*Cd(8059)*MStop2(3) - 
     &   2*MTy2*Cd(7420)*UStop2(1,2)

	Cd(8087) = Cd(8070) - MGl2*Cd(8086)

	Cd(8088) = Cd(8082) - MGl2*Cd(8087)

	Cd(798) = 2*Cd(796) - At*Cd(797)*MStop2(3)

	Cd(8089) = Cd(44)*Cd(8076) - MGl2*Cd(8088)

	Cd(8098) = Cd(8089) - MTy2*Cd(50)*Cd(8097)*UStop2(1,2)

	Cd(8099) = Cd(47)**2*Cd(8062)*MStop2(3) - Cd(8098)*UStop2(1,2)

	Cd(8128) = Cd(44)**2*Cd(48)**2*Cd(8109) - MGl2*Cd(8127)

	Cd(8130) = Cd(8128) - Cd(50)**2*Cd(8129)*UStop2(1,2)

	Cd(8131) = Cd(47)*Cd(8108) - Cd(8130)*UStop2(1,2)

	Cd(8168) = 2*Cd(8141) + MGl2*Cd(8167)

	Cd(8169) = Cd(8150) - MGl2*Cd(8168)

	Cd(8170) = 2*Cd(48)*Cd(8145) - MGl2*Cd(8169)

	Cd(8171) = 3*Cd(44)**2*Cd(48)**3*Cd(89)*Cd(187) - MGl2*Cd(8170)

	Cd(8176) = Cd(8163) - MGl2*Cd(8175)

	Cd(8177) = Cd(44)*Cd(48)*Cd(8157) + MGl2*Cd(8176)

	Cd(8178) = Cd(8171)*UStop(1,1) - 2*M_3*MTy*Cd(8177)*UStop(1,2)

	Cd(8207) = 4*Cd(8194) - MGl2*Cd(8206)

	Cd(8208) = Cd(8199) - MGl2*Cd(8207)

	Cd(8209) = 2*Cd(48)*Cd(8195) + MGl2*Cd(8208)

	Cd(8210) = Cd(44)**2*Cd(48)**4 - MGl2*Cd(8209)

	Cd(8212) = Cd(8205) - 2*MGl2*Cd(8211)

	Cd(8213) = Cd(44)*Cd(48)*Cd(89)*Cd(5556) - MGl2*Cd(8212)

	Cd(8214) = MGl2**4*Cd(8200) - Cd(8213)*MStop2(3)

        Cd(8215) = Cd(8210)*MStop2(3)*UStop(1,1) + 
     &   M_3*MTy*Cd(8214)*UStop(1,2)

	Cd(8217) = 2*Cd(44)**2*Cd(48)*Cd(8182) - MGl2*Cd(8216)

	Cd(8218) = Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(8217)

	Cd(8220) = Cd(8192) - 2*MGl2*Cd(8219)

	Cd(8221) = Cd(44)*Cd(48)*Cd(89)*Cd(5556) - MGl2*Cd(8220)

	Cd(8222) = MGl2**4*Cd(8189) + Cd(8221)*MStop2(3)

        Cd(8239) = Cd(8210)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(8238)*MStop2(3)*UStop(1,2)

        Cd(8253) = Cd(8252)*MStop2(3)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(6726)*UStop(1,2)

	Cd(8258) = 2*Cd(8257)*MStop2(3) + Cd(6737)*UStop2(1,1)

        Cd(8261) = MTy*Cd(8218)*MStop2(3)*UStop(1,1) - 
     &   2*M_3*MTy2*Cd(8222)*UStop(1,2)

	Cd(8266) = -(MGl2*Cd(8206)) + 4*Cd(8244)

	Cd(8267) = Cd(8199) - MGl2*Cd(8266)

	Cd(8275) = Cd(44)**2*Cd(48)**2*Cd(288) - MGl2*Cd(8217)

	Cd(8282) = Cd(8278) - MGl2*Cd(8281)

	Cd(8283) = -(Cd(44)*Cd(48)*Cd(6826)) + MGl2*Cd(8282)

	Cd(8302) = Cd(8297) - MGl2*Cd(8301)

	Cd(8303) = Cd(44)*Cd(48)*Cd(8295) - MGl2*Cd(8302)

        Cd(8305) = Cd(50)*Cd(8304)*MStop2(3) + 
     &   2*M_3C*MTy*Cd(8303)*UUStop(1,1)

	Cd(820) = -(S2A*Cd(798)) + 2*CA2*Cd(819)

        Cd(821) = -(1/(128.D0*Pi**3)*
     &     (Alfas2L*(Cd(820) - 
     &          4*SA2*Cd(13)**2*Cr(75)**2*MStop2(3)**2))/
     &      (SB2*SW2*MStop2(3)**2))

	Cd(8318) = Cd(6650) + 2*Cd(8239)*MStop2(3)

	Cd(8320) = -(Cd(47)*Cd(8319)*UStop(1,1)) + Cd(8318)*UStop2(1,2)

	Cd(8321) = Cd(47)**2*MGlpTmSt4(2) + Cd(50)*Cd(8235)*UStop2(1,2)

	Cd(8338) = Cd(8320) + Cd(47)*Cd(50)*Cd(6603)*UStop2(1,1)

        Cd(8372) = 2*(Cd(47)*Cd(8371)*MStop2(3) + 
     &     M_3*Cd(5217)*UStop(1,2))

	Cd(828) = (Cd(827) + 2*Cd(826)*MStop2(2))/MStop2(2)

	Cd(8388) = Cd(8384) - MGl2*Cd(8387)

	Cd(8389) = 2*Cd(44)*Cd(48)*Cd(7155) + MGl2*Cd(8388)

	Cd(8390) = Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(8389)

	Cd(8393) = Cd(8349) - MGl2*Cd(8392)

	Cd(8394) = 2*Cd(44)*Cd(48)*Cd(7155) + MGl2*Cd(8393)

	Cd(8395) = Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(8394)

	Cd(829) = S2A*Cd(823) + 4*SA2*Cd(13)**2*Cr(75)*MStop2(4)

        Cd(8398) = 2*M_3C*MTy*Cd(8238)*MStop2(3)*UCStop(1,3) + 
     &   Cd(8390)*UStop2(1,2)

        Cd(830) = 1/(64.D0*Pi**3)*
     &   (Alfas2L*Cr(9)*(Cd(829) + CA2*Cd(828)*MStop2(2)))/
     &    (SB2*SW2*MStop2(2))

	Cd(8412) = Cd(8409) - MGl2*Cd(8411)

	Cd(8413) = 4*Cd(44)**2*Cd(48)**2*Cd(198) - MGl2*Cd(8412)

	Cd(8414) = Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(8413)

	Cd(8437) = Cd(8427) + MGl2*Cd(8436)

	Cd(8438) = Cd(8433) - MGl2*Cd(8437)

	Cd(8439) = Cd(44)*Cd(48)*Cd(8423) - MGl2*Cd(8438)

	Cd(8440) = Cd(44)*Cd(48)*Cd(6438) - MGl2*Cd(7546)

	Cd(8441) = 2*Cd(8440) - Cd(8238)*UStop2(1,2)

	Cd(8442) = Cd(8439) - 4*MTy2**2*Cd(8441)*UStop2(1,2)

	Cd(8447) = Cd(8384) - MGl2*Cd(8446)

	Cd(8448) = 2*Cd(44)*Cd(48)*Cd(7155) + MGl2*Cd(8447)

	Cd(8449) = Cd(44)**2*Cd(48)**2*Cd(2048) - MGl2*Cd(8448)

        Cd(8450) = 3*Cd(8449)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(8238)*MStop2(3)*UStop(1,2)

        Cd(81) = Cd(80)*T134(MGl2,MTy2,MStop2(2)) + 
     &   Cd(64)*MGlpTmSt2(2)*MStop2(3)*UStop2(1,1)*UStop2(1,2)

	Cd(8466) = -Cd(8427) + MGl2*Cd(8465)

	Cd(8467) = Cd(8461) + MGl2*Cd(8466)

	Cd(8468) = Cd(44)*Cd(48)*Cd(8457) - MGl2*Cd(8467)

	Cd(8488) = Cd(8205) - MGl2*Cd(8281)

	Cd(8489) = Cd(44)*Cd(48)*Cd(89)*Cd(5556) - MGl2*Cd(8488)

	Cd(8493) = Cd(8490) - MGl2*Cd(8492)

	Cd(8494) = Cd(44)*Cd(48)*Cd(5556) - MGl2*Cd(8493)

        Cd(8521) = MGlmT2*Cd(8520)*MStop2(2) - 
     &   Cd(302)*Cd(8508)*UStop2(1,2)

	Cd(8529) = Cd(8490) - MGl2*Cd(8528)

	Cd(8530) = Cd(44)*Cd(48)*Cd(5556) - MGl2*Cd(8529)

        Cd(8531) = Cd(8527)*MStop2(3)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(8530)*UStop(1,2)

	Cd(8540) = Cd(8531)*UStop2(1,1) - 2*Cd(7736)*UStop2(1,2)

        Cd(8541) = Cd(8539)*MStop2(3)*UStop2(1,1) - 
     &   Cd(302)*Cd(8540)*UStopC(1,1)

        Cd(8550) = (UStop(1,1)*
     &     (2*M_3C*MTy*Cd(47)**2*Cd(50)**2*Cd(8549)*MStop2(3)*
     &        UCStop(1,3) + Cd(8548)*UStop2(1,2)))/
     &   (Cd(47)**2*Cd(50)**2*MStop2(3))

	Cd(82) = Cd(50)*Cd(78) - Cd(47)*Cd(81)

	Cd(8557) = MGl2*Cd(6095) + Cd(8556)

	Cd(8561) = -(Cd(44)*Cd(48)*Cd(6826)) + Cd(8489)*UStop2(1,2)

	Cd(8565) = Cd(44)*Cd(48)*Cd(89)*Cd(5556) + MGl2**2*Cd(8564)

	Cd(8566) = -4*MGl2*Cd(8483) + Cd(8565)*MStop2(3)

	Cd(8571) = Cd(89)*Cd(8569) + MGl2*Cd(8570)

	Cd(8581) = Cd(302)*Cd(8580) + M_3C*MTy*Cd(8489)*UUStop(1,1)

	Cd(8601) = 2*Cd(48)*Cd(8195) + MGl2*Cd(8267)

	Cd(8602) = Cd(44)**2*Cd(48)**4 - MGl2*Cd(8601)

	Cd(8615) = Cd(8614) + Cd(50)*Cd(8346)*UStop2(1,2)

	Cd(8616) = Cd(47)**2*Cd(251) + Cd(8615)*UStop2(1,2)

	Cd(8618) = Cd(8616) + Cd(8617)*UStop2(1,1)

	Cd(8620) = Cd(44)**2*Cd(48)**2 - MGl2*Cd(8619)

	Cd(8630) = Cd(47)*Cd(8620) - Cd(8629)*UStop2(1,2)

	Cd(8631) = -2*MGl2*Cd(8211) + Cd(8278)

	Cd(8632) = -(Cd(44)*Cd(48)*Cd(6826)) + MGl2*Cd(8631)

	Cd(8633) = MGl2**4*Cd(8200) + Cd(8632)*MStop2(3)

        Cd(8634) = Cd(8602)*MStop2(3)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(8633)*UStop(1,2)

	Cd(8635) = Cd(8630) - Cd(8218)*UStop2(1,1)

	Cd(8636) = Cd(8634)*UStop2(1,1) + 2*Cd(8215)*UStop2(1,2)

        Cd(8645) = 2*M_3C*MTy*Cd(8283)*UUStop(1,1) + 
     &   Cd(50)**2*MGlpTmSt4(1)*UUStop(1,3)

        Cd(8649) = -2*M_3C*MTy*Cd(47)*Cd(8356)*UCStop(1,3)*
     &    UStop(1,1) + Cd(8239)*UStop2(1,2)

        Cd(8651) = 2*M_3C*MTy*Cd(8650)*UUStop(1,1) + 
     &   Cd(8321)*UUStop(1,3)

        Cd(8652) = Cd(8649)*UStop2(1,1)*UStopC(1,1) + 
     &   Cd(8651)*UUStopC(1,3)

        Cd(8657) = Cd(4273)*UStop2(1,1) - 
     &   2*M_3*MTy*Cd(4217)*UCStopC(1,3)*UStop2(1,2)

        Cd(8667) = Cd(8666)*UStop(1,1) + 
     &   8*M_3*MTy*Cd(4217)*MStop2(3)*UStop(1,2)

        Cd(8671) = Cd(2317)*MStop2(1)*UStop(1,1) + 
     &   4*M_3*MTy*Cd(8670)*MStop2(3)*UStop(1,2)

	Cd(8686) = Cd(8657) + Cd(47)*Cd(4276)*UStop2(1,2)*UStopC(1,2)

	Cd(8695) = Cd(44)**3 - MGl2*Cd(8694)

	Cd(8696) = Cd(4217) + Cd(8695)*UStop2(1,2)

	Cd(8697) = 2*M_3*MTy*Cd(8696)*MStop2(3) + Cd(8693)*UCStop(1,3)

	Cd(8703) = 4*Cd(4223) - 2*Cd(4229)*UStop2(1,2)

	Cd(8725) = Cd(8716) + 4*Cd(8724)*UStop2(1,2)

        Cd(8745) = MTy*Cd(8741)*UStop(1,1) + 
     &   M_3*Cd(8744)*MStop2(3)*UStop(1,2)

        Cd(8753) = Cd(8749)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(8752)*MStop2(3)*UStop(1,2)

	Cd(8758) = M_3*Cd(4327)*MStop2(3) + 2*MTy*Cd(8757)*UCStop(1,3)

	Cd(866) = MUE*Cd(851)*UCStop(1,3) + MUEC*Cd(859)*UCStopC(1,3)

	Cd(8764) = 3*Cd(44)**3 - MGl2*Cd(8763)

	Cd(8783) = M_3C*MTy*Cd(4221) + 2*Cd(8782)*UCStopC(1,3)

        Cd(8787) = Cd(8749)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(8728)*MStop2(3)*UStop(1,2)

        Cd(8800) = Cd(8797)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(8799)*MStop2(3)*UStop(1,2)

        Cd(8801) = 2*MStop2(1)*
     &   (3*Cd(2317)*MStop2(1)*UStop(1,1) + 
     &     4*M_3*MGl2*MTy*Cd(157)*MStop2(3)**2*UStop(1,2))

        Cd(8840) = -(Cd(118)*Cd(8839)*MStop2(1)) + 
     &   Cd(8801)*UStop(1,1)*UStop2(1,2)

	Cd(8848) = Cd(8847) + 8*MTy2*Cd(4219)*UStop2(1,2)

        Cd(8849) = 2*Cd(8843)*MStop2(1)*MStop2(3) - 
     &   MTy2*Cd(8848)*UStop2(1,2)

	Cd(8862) = 4*Cd(8861)*MStop2(3) - Cd(8860)*UUStopC(1,2)

	Cd(8863) = MTy*Cd(8862) - M_3C*Cd(8849)*UCStop(1,3)

	Cd(8864) = Cd(4450) - 8*MGl2*MTy2*Cd(152)*MStop2(3)*UStop2(1,2)

        Cd(8868) = 3*Cd(8867)*UStop(1,1) - 
     &   2*M_3*MTy*Cd(4219)*UStop(1,2)

	Cd(8873) = 2*Cd(44)**3*Cd(3736) + MGl2*Cd(8872)

        Cd(8877) = 2*(MTy*Cd(8873)*UStop(1,1) + 
     &     M_3*Cd(8876)*MStop2(1)*MStop2(3)*UStop(1,2))

        Cd(8878) = M_3*MTy*Cd(8847)*UUStop(1,2) + 
     &   2*Cd(47)**2*MStop2(1)*UUStop(2,3)

	Cd(8879) = Cd(8877) + 4*MTy*MTy2*Cd(8868)*UStop2(1,2)**2

        Cd(8890) = 2*Cd(4619)*UStop2(1,1) + 
     &   4*M_3*MTy*Cd(4562)*UCStopC(1,3)*UStop2(1,2)

        Cd(8900) = Cd(8899)*UStop(1,1) - 
     &   8*M_3*MTy*Cd(254)*Cd(4739)*MStop2(3)*UStop(1,2)

        Cd(8907) = Cd(8903)*UStop(1,1) - 
     &   4*M_3*MTy*Cd(8906)*MStop2(3)*UStop(1,2)

	Cd(8917) = Cd(8890) - 2*Cd(50)*Cd(4622)*UStop2(1,2)*UStopC(1,2)

	Cd(8928) = Cd(48)**3 - MGl2*Cd(8927)

	Cd(8929) = Cd(4562) + Cd(8928)*UStop2(1,2)

	Cd(8930) = 2*M_3*MTy*Cd(8929)*MStop2(3) + Cd(8926)*UCStop(1,3)

	Cd(8935) = 4*Cd(4569) - 2*Cd(4576)*UStop2(1,2)

	Cd(8940) = 2*Cd(50)*Cd(4666) - Cd(4677)*UStop2(1,1)

        Cd(8974) = MTy*Cd(8970)*UStop(1,1) + 
     &   M_3*Cd(8973)*MStop2(3)*UStop(1,2)

        Cd(8979) = Cd(8978)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(4894)*MStop2(3)*UStop(1,2)

        Cd(8980) = Cd(8979)*UStop2(1,1) + 
     &   2*Cd(8957)*MStop2(3)*UStop2(1,2)

	Cd(8981) = Cd(8974) - 4*MTy*Cd(8980)

        Cd(8984) = Cd(8899)*UStop(1,1) - 
     &   8*M_3*MTy*Cd(4562)*MStop2(3)*UStop(1,2)

        Cd(8985) = -(Cd(1756)*MStop2(2)**4*UStop(1,1)) + 
     &   Cd(8907)*MStop2(2)*UStop2(1,2)

        Cd(8987) = -2*Cd(8917)*MStop2(3) + 
     &   2*Cd(8984)*UStop2(1,1)**2*UStopC(1,1)

        Cd(8994) = Cd(8978)*UCStop(1,3) + 
     &   M_3*MTy*Cd(4562)*MStop2(3)*UStop2(1,2)

	Cd(8995) = 2*MTy*Cd(8994) + M_3*Cd(8993)*MStop2(3)

	Cd(9001) = M_3*Cd(4671) - 2*MTy*Cd(4753)*UUStopC(1,2)

	Cd(9002) = Cd(9001)*MStop2(3) + 2*MTy*Cd(9000)*UCStop(1,3)

	Cd(9017) = M_3C*MTy*Cd(4567) - 2*Cd(9016)*UCStopC(1,3)

        Cd(9021) = Cd(8979)*UStop2(1,1) + 
     &   2*Cd(9020)*MStop2(3)*UStop2(1,2)

	Cd(9022) = Cd(8974) - 4*MTy*Cd(9021)

	Cd(9025) = 2*MTy*Cd(8994) + M_3*Cd(4737)*MStop2(3)

        Cd(9039) = Cd(9036)*UStop(1,1) + 
     &   2*M_3*MTy*Cd(9038)*MStop2(3)*UStop(1,2)

        Cd(9040) = 3*Cd(8903)*UStop(1,1) + 
     &   4*M_3*MGl2*MTy*Cd(184)*MStop2(3)**2*UStop(1,2)

        Cd(9073) = -(Cd(138)*Cd(9072)) + 
     &   2*Cd(9040)*UStop(1,1)*UStop2(1,2)

	Cd(9089) = 4*Cd(9088)*MStop2(3) - Cd(9087)*UUStopC(1,2)

	Cd(9090) = MTy*Cd(9089) + 2*M_3C*Cd(4860)*UCStop(1,3)

	end
