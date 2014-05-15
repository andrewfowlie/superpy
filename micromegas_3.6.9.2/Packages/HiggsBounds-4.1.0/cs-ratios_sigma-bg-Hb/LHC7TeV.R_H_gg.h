!******************************************************
        real*8 function lhc7_rH_gg(x)
!******************************************************
!* ggH process : consensus numbers used (PDF: MSTW 2008 NNLO)
!*		 see: arXiv:1101.0593 [hep-ph] 
!* 		 numbers taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt7TeV
!* bbH process : numbers in range [70:1020] generated with bbh@nnlo v. 1.3 (PDF: MSTW 2008 NNLO) 
!* 		 with input parameters according to arXiv:1101.0593 [hep-ph]
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:1000], deviations from data below 0.012%
!*      slight extrapolation to range [70:1020] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,x
        real*8 a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2

        i1 = 593099.89872571
        j1 = -39709.0334757866
        k1 = 1143.7934255263
        a1 = -17.7227682955945
        b1 = 0.191890374015844
        c1 = -0.00128789621013637
        d1 = 5.73969299012758e-06
        e1 = -1.68044503342273e-08
        f1 = 3.10098872124758e-11
        g1 = -3.265795492463e-14
        h1 = 1.49463610979039e-17

        i2 = 62461446.2253099
        j2 = -1030866.64181057
        k2 = 7410.55699727643
        a2 = -29.6033746719771
        b2 = 0.0804362310250513
        c2 = -0.000140562830181593
        d2 = 1.65223884829233e-07
        e2 = -1.28738488029277e-10
        f2 = 6.34186901246237e-14
        g2 = -1.77318371831275e-17
        h2 = 2.1150324470805e-21

        lhc7_rH_gg=0d0

        if(x .lt. 70d0) then
        write(*,*)'function lhc7_rH_gg might not be a good fit (m_H < 70 GeV)'
        endif
        if((x .ge. 70d0) .and. (x .le. 350d0)) then
        lhc7_rH_gg=i1/x**3+j1/x**2+k1/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        endif
        if((x .ge. 350d0) .and. (x .le. 1020d0)) then
        lhc7_rH_gg=i2/x**3+j2/x**2+k2/x+a2+b2*x+c2*x**2+d2*x**3+e2*x**4+f2*x**5+g2*x**6+h2*x**7
        endif
        if(x .gt. 1020d0) then
        write(*,*)'function lhc7_rH_gg might not be a good fit (m_H > 1020 GeV)'
        endif
                                                                                                                                        
        end function
        
